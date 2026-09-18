//! AUTOMORPHISM FINDING FOR REFERENCE SHAPES
use nalgebra::{Matrix3, Vector3};

use crate::cshm::bounds::*;
use crate::cshm::linalg::*;
use crate::cshm::permutations::SearchTables;

/// Finds a Reference Shape's Automorphisms.
/// Expects centered and normalised Shapes.
/// `has_centre` marks whether `reference[0]` is the shape's centre point (as opposed
/// to a vertex): the centre can only ever map to itself, so it is fixed up front
/// instead of being explored as a branch.
pub fn find_automorphisms(reference: &[Vector3<f64>], has_centre: bool) -> Vec<Vec<usize>> {
    let n = reference.len();
    let mut automorphisms: Vec<Vec<usize>> = Vec::new();
    // A shape's automorphisms are its zero-deviation permutations against itself, so the search
    // tables are those of a reference-vs-reference search.
    let tables = SearchTables::new(reference, reference);

    let mut assigned = vec![false; n];
    let mut current_perm = Vec::with_capacity(n);
    let mut h_partial: Matrix3<f64> = Matrix3::zeros();

    if has_centre {
        assigned[0] = true;
        current_perm.push(0);
        h_partial = tables.hi[0][0];
    }

    const EPS: f64 = 1e-6;

    automorphism_branch(
        &tables,
        EPS,
        &mut assigned,
        &mut current_perm,
        &mut h_partial,
        &mut automorphisms,
    );
    automorphisms
}

fn automorphism_branch(
    tables: &SearchTables,
    epsilon: f64,
    assigned: &mut [bool],
    current_perm: &mut Vec<usize>,
    h_partial: &mut Matrix3<f64>,
    automorphisms: &mut Vec<Vec<usize>>,
) {
    let reference = tables.reference;
    let n = reference.len();
    // If the permutation is complete -> Score it using full SVD
    if current_perm.len() == n {
        // Build the reordered shape
        let reordered: Vec<Vector3<f64>> = current_perm.iter().map(|&i| reference[i]).collect();
        let h = correlation_matrix(reference, &reordered);
        let (_, a_i) = optimal_rotation(h);
        let s = shape_measure(&a_i, n);

        // If permutation scores a CShM ~ 0 it is an automorphism.
        if s.abs() < epsilon {
            automorphisms.push(current_perm.clone());
        }
        return;
    }

    let pos = current_perm.len();

    for ref_idx in 0..n {
        // Main loop to look for perms.
        if assigned[ref_idx] {
            // Skip the currently assigned points in the permutation
            continue;
        }
        // Choosing the next point starts here:

        //
        *h_partial += tables.hi[ref_idx][pos];
        assigned[ref_idx] = true;

        let a_partial: f64 = nuclear_norm(*h_partial);
        let remaining_bound =
            max_unassigned_norm(&tables.ref_norms, assigned) * tables.prob_suffix[pos + 1];
        let a_bound = a_partial + remaining_bound;
        let s_bound = (1.0 - a_bound.powi(2) / (n as f64).powi(2)) * 100.0;

        if s_bound < epsilon {
            current_perm.push(ref_idx);
            automorphism_branch(
                tables,
                epsilon,
                assigned,
                current_perm,
                h_partial,
                automorphisms,
            );

            current_perm.pop();
        }
        assigned[ref_idx] = false;
        *h_partial -= tables.hi[ref_idx][pos];
    }
}
