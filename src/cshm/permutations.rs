//! MAIN FUNCTIONS OF CSHM, OPTIMAL PERMUTATION FINDING FOR A REFERENCE SHAPE GIVEN A NON-ALIGNED PROBLEM SHAPE.

use crate::cshm::bounds::*;
use crate::cshm::linalg::*;
use crate::geometry::center_and_normalise;
use nalgebra::{Matrix3, Vector3};

/// The two centred, normalised point sets and the tables
/// precomputed from them once per search.
pub(super) struct SearchTables<'a> {
    pub(super) reference: &'a [Vector3<f64>],
    pub(super) problem: &'a [Vector3<f64>],
    /// `hi[r][p]` is the correlation block `reference[r] * problem[p]^T`.
    pub(super) hi: Vec<Vec<Matrix3<f64>>>,
    /// Norm of every reference point.
    pub(super) ref_norms: Vec<f64>,
    /// `prob_suffix[p]` is the sum of the problem-point norms from `p` to the end.
    pub(super) prob_suffix: Vec<f64>,
}

impl<'a> SearchTables<'a> {
    /// Precomputes the hi correlation blocks, the reference norms and the
    /// problem norms sum.
    pub(super) fn new(reference: &'a [Vector3<f64>], problem: &'a [Vector3<f64>]) -> Self {
        Self {
            reference,
            problem,
            hi: precompute_correlation_blocks(reference, problem),
            ref_norms: precompute_norms(reference),
            prob_suffix: precompute_suffix_sums(&precompute_norms(problem)),
        }
    }
}

/// Encodes the recursion state.
struct SearchState {
    // Currently assigned points.
    assigned: Vec<bool>,
    // Current permutation being explored.
    current_perm: Vec<usize>,
    /// Correlation matrix of the partial permutation, kept up to date incrementally.
    h_partial: Matrix3<f64>,
    // Best scores found yet.
    best_s: f64,
    best_perm: Vec<usize>,
    best_rot_matrix: Matrix3<f64>,
}

/// Recursively finds the best permutation of a given reference shape so that its points align
/// to the problem shape. Will prune non-optimal permutations using the partial sum of the singular
/// values.
/// Fixes the permutation of the central atom for centered structures.
pub(crate) fn find_best_permutation(
    reference: &mut [Vector3<f64>],
    problem: &mut [Vector3<f64>],
    has_centre: bool,
) -> (f64, Vec<usize>, Vec<Vector3<f64>>, Matrix3<f64>) {
    let n = problem.len();
    debug_assert_eq!(n, reference.len());

    center_and_normalise(reference);
    let (problem_centroid, normalisation_constant) = center_and_normalise(problem);

    // Precompute the correlation matrices and norms for all points once.
    let tables = SearchTables::new(reference, problem);

    // Initialise the recursion state and the return values.
    let mut state = SearchState {
        assigned: vec![false; n],
        current_perm: Vec::with_capacity(n),
        h_partial: Matrix3::zeros(),
        best_s: f64::INFINITY,
        best_perm: Vec::new(),
        best_rot_matrix: Matrix3::zeros(),
    };

    // Fixes permutation of the central atom if found.
    if has_centre {
        state.assigned[0] = true;
        state.current_perm.push(0);
        state.h_partial = tables.hi[0][0];
    }

    branch(&tables, &mut state);

    let reconstructed: Vec<Vector3<f64>> = reference
        .iter()
        .map(|p| {
            (state.best_rot_matrix.transpose() * p) / normalisation_constant + problem_centroid
        })
        .collect();

    (
        state.best_s,
        state.best_perm,
        reconstructed,
        state.best_rot_matrix,
    )
}

fn branch(tables: &SearchTables, state: &mut SearchState) {
    let n = tables.reference.len();
    debug_assert_eq!(n, tables.problem.len());

    if state.current_perm.len() == n {
        // If a permutation is complete then:
        let reordered: Vec<Vector3<f64>> = state
            .current_perm
            .iter()
            .map(|&p| tables.reference[p])
            .collect();
        let h = correlation_matrix(tables.problem, &reordered);
        let (rot_matrix, a_i) = optimal_rotation(h);
        let s = shape_measure(&a_i, n).max(0.0); // max 0.0 makes sure the s value doesn't go below 0 because floating point errors.

        if s < state.best_s {
            state.best_s = s;
            state.best_perm = state.current_perm.clone();
            state.best_rot_matrix = rot_matrix;
        }
        return;
    }

    let pos = state.current_perm.len(); // Next problem-point index to assign

    for ref_idx in 0..n {
        if state.assigned[ref_idx] {
            continue;
        }

        state.h_partial += tables.hi[ref_idx][pos]; // Sum the corresponding point to the partial correlation matrix
        state.assigned[ref_idx] = true; // Mark the point as assigned.

        let a_partial: f64 = nuclear_norm(state.h_partial); // Calculate the partial SV sum,
        // Calculate the estimated remaining contributions
        // to the correlation matrix measure of the rest of points.
        let remaining_bound =
            max_unassigned_norm(&tables.ref_norms, &state.assigned) * tables.prob_suffix[pos + 1];
        let a_bound = a_partial + remaining_bound;
        let s_bound = (1.0 - a_bound.powi(2) / ((n as f64).powi(2))) * 100.0;

        if s_bound < state.best_s {
            // If we have found a better s
            state.current_perm.push(ref_idx); // Add the matrix when pushing new point to list.
            // Recursively call the branch function again.
            branch(tables, state);
            state.current_perm.pop();
        }
        state.assigned[ref_idx] = false;
        state.h_partial -= tables.hi[ref_idx][pos]; // Subtract the matrix when backtracking the current
    }
}
