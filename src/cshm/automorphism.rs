use nalgebra::{Matrix3, Vector3};

use crate::cshm::geometry::*;
use crate::cshm::bounds::*;


/// Finds a Reference Shape's Automorphisms.
/// Expects centered and normalised Shapes.
pub fn find_automorphisms(reference: &[Vector3<f64>]) -> Vec<Vec<usize>> {
    let n = reference.len();
    let mut automorphisms: Vec<Vec<usize>> = Vec::new();
    let hi = precompute_correlation_blocks(reference, reference);

    let mut assigned = vec![false; n];
    let mut current_perm = Vec::with_capacity(n);
    let mut h_partial: Matrix3<f64> = Matrix3::zeros();

    const EPS:f64 = 1e-6;

    automorphism_branch(
        reference,
        &hi,
        &mut assigned,
        &mut current_perm,
        &mut h_partial,
        EPS,
        &mut automorphisms,
    );
    automorphisms
}

fn automorphism_branch(
    reference: &[Vector3<f64>],
    hi: &Vec<Vec<Matrix3<f64>>>,
    assigned: &mut [bool],
    current_perm: &mut Vec<usize>,
    h_partial: &mut Matrix3<f64>,
    epsilon: f64,
    automorphisms: &mut Vec<Vec<usize>>,
) {
    let n = reference.len();
    // If the permutation is complete -> Score it using full SVD
    if current_perm.len() == n {
        // Build the reordered shape
        let reordered: Vec<Vector3<f64>> = current_perm.iter().map(|&i| reference[i]).collect();
        let h = crate::cshm::correlation_matrix(reference, &reordered);
        let (_, a_i) = crate::cshm::optimal_rotation(h);
        let s = crate::cshm::shape_measure(&a_i, n);

        // If permutation scores a CShM ~ 0 it is an automorphism.
        if s.abs() < epsilon {
            automorphisms.push(current_perm.clone());
        }
        return;
    }

    let pos = current_perm.len();

    for ref_idx in 0..n {
        // Main loop to look for perms.
        if assigned[ref_idx] { // Skip the currently asigned points in the permutation
            continue;
        }
        // Choosing the next point starts here:

        //
        *h_partial += hi[ref_idx][pos];
        assigned[ref_idx] = true;

        let a_partial: f64 = singular_values(*h_partial).iter().sum();
        let remaining_bound =
            max_unassigned_norm(reference, assigned) *
                unassigned_norms_sum(&reference[pos + 1..]);
        let a_bound = a_partial + remaining_bound;
        let s_bound = (1.0 - a_bound.powi(2) / (n as f64).powi(2)) * 100.0;

        if s_bound < epsilon {
            current_perm.push(ref_idx);
            automorphism_branch(
                reference, &hi, assigned, current_perm,
                h_partial, epsilon, automorphisms,
            );

            current_perm.pop();
        }
        assigned[ref_idx] = false;
        *h_partial -= hi[ref_idx][pos];
    }
}

