/// MAIN FUNCTIONS OF CSHM, OPTIMAL PERMUTATION FINDING FOR A REFERENCE SHAPE GIVEN A NON-ALIGNED PROBLEM SHAPE.

use std::collections::HashSet;
use std::time::Instant;
use nalgebra::{Matrix3, Vector3};

use crate::cshm::geometry::*;
use crate::cshm::bounds::*;
use crate::cshm::automorphism::*;

pub fn best_permutation_branch_and_bound(
    reference: &mut [Vector3<f64>],
    problem: &mut [Vector3<f64>],
) -> (f64, Vec<usize>) {
    let n = problem.len();
    debug_assert_eq!(n, reference.len());

    center_and_normalise(reference);
    center_and_normalise(problem);

    let start = Instant::now();
    let ref_automorphisms: Vec<Vec<usize>> = find_automorphisms(reference);
    let time = start.elapsed();
    println!("Automorphisms Branch Time: {:?}", time);

    let mut visited: HashSet<Vec<usize>> = std::collections::HashSet::new();

    let mut best_s = f64::INFINITY;
    let mut best_perm: Vec<usize> = Vec::new();

    let mut assigned = vec![false; n];
    let mut current_perm: Vec<usize> = Vec::with_capacity(n);

    let mut h_partial = Matrix3::zeros();

    let hi = precompute_correlation_blocks(&reference, &problem);

    branch(
        &reference,
        &problem,
        &hi,
        &ref_automorphisms,
        &mut visited,
        &mut assigned,
        &mut current_perm,
        &mut h_partial,
        &mut best_s,
        &mut best_perm
    );

    (best_s, best_perm)
}

fn branch(
    reference: &[Vector3<f64>],
    problem: &[Vector3<f64>],
    hi: &Vec<Vec<Matrix3<f64>>>,
    ref_automorphisms: &Vec<Vec<usize>>,
    visited: &mut std::collections::HashSet<Vec<usize>>,
    assigned: &mut [bool],
    current_perm: &mut Vec<usize>,
    h_partial: &mut Matrix3<f64>,
    best_s: &mut f64,
    best_perm: &mut Vec<usize>,

) {

    let n = reference.len();
    debug_assert_eq!(n, problem.len());


    if current_perm.len() == n { // If a permutation is complete then:
        if visited.contains(current_perm) {
            return;
        }

        let reordered: Vec<Vector3<f64>> = current_perm.iter().map(|&p| reference[p]).collect();
        let h = correlation_matrix(problem, &reordered);
        let (_, a_i) = optimal_rotation(h);
        let s = shape_measure(&a_i, n);

        for a in ref_automorphisms {
            let equiv: Vec<usize> = (0..n).map(|i| a[current_perm[i]]).collect();
            visited.insert(equiv);
        }

        if s < *best_s {
            *best_s = s;
            *best_perm = current_perm.clone();
        }
        return;
    }

    let pos = current_perm.len(); // Next problem-point index to assign

    for ref_idx in 0..n {
        if assigned[ref_idx] {
            continue;
        }

        *h_partial += hi[ref_idx][pos]; // Sum the corresponding point to the partial correlation matrix
        assigned[ref_idx] = true; // Mark the point as assigned.

        let a_partial: f64 = singular_values(*h_partial).iter().sum(); // Calculate the partial SV sum,
        // Calculate the estimated remaining contributions
        // to the correlation matrix measure of the rest of points.
        let remaining_bound = max_unassigned_norm(reference, assigned) * unassigned_norms_sum(&problem[pos + 1..]);
        let a_bound = a_partial + remaining_bound;
        let s_bound = (1.0 - a_bound.powi(2)/((n as f64).powi(2))) * 100.0;


        if s_bound < *best_s { // If we have found a better s
            current_perm.push(ref_idx); // Add the matrix when pushing new point to list.
            // Recursively call the branch function again.
            branch(
                reference, problem, hi, ref_automorphisms, visited,
                assigned, current_perm, h_partial, best_s, best_perm,
            );
            current_perm.pop();

        }
        assigned[ref_idx] = false;
        *h_partial -= hi[ref_idx][pos]; // Subtract the matrix when backtracking the current
    }
}