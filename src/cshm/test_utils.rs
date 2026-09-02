#![cfg(test)]
use crate::cshm::linalg::*;
use crate::geometry::center_and_normalise;
use itertools::Itertools;
use nalgebra::Vector3;
/// BRUTE FORCE FUNCTIONS USED FOR SANITY-CHECKS AND TESTING. NOT USING THEM IN FINAL BUILDS.
use std::time::Instant;

/// Brute-force search all n! automorphisms of a Reference Shape. Expects centered and normalised coordinates to work.
pub(crate) fn naive_find_automorphism(points: &[Vector3<f64>]) -> Vec<Vec<usize>> {
    // Expects centered and normalised coordinates.
    let n = points.len();
    let mut automorphisms: Vec<Vec<usize>> = Vec::new();

    for perm in (0..n).permutations(n) {
        // Build the reordered point set using this permutation:
        let reordered: Vec<Vector3<f64>> = perm.iter().map(|&p| points[p]).collect();

        let h = correlation_matrix(points, &reordered);
        let (_, a_i) = optimal_rotation(h);
        let s = shape_measure(&a_i, n);

        if s.abs() < 1e-6 {
            automorphisms.push(perm);
        }
    }
    automorphisms
}


/// Finds the best permutation by iterating all over the n! permutation list
/// but prunes by automorphisms using naive_find_automorphism.
pub(crate) fn best_permutation_brute_force(
    reference: &mut [Vector3<f64>],
    problem: &mut [Vector3<f64>],
) -> (f64, Vec<usize>) {
    let n = problem.len();
    debug_assert_eq!(n, reference.len());

    // Normalise inputs
    center_and_normalise(reference);
    center_and_normalise(problem);

    let mut best_s = f64::INFINITY;
    let mut best_perm: Vec<usize> = Vec::new();

    let naive_start = Instant::now();
    let ref_automorphisms: Vec<Vec<usize>> = naive_find_automorphism(reference);
    let naive_time = naive_start.elapsed();
    println!("naive time: {:?}", naive_time);
    let mut visited: std::collections::HashSet<Vec<usize>> = std::collections::HashSet::new();

    for perm in (0..n).permutations(n) {
        // Build the reordered point set using this permutation:
        let reordered: Vec<Vector3<f64>> = perm.iter().map(|&p| reference[p]).collect();

        // If this permutation of the reference shape is an automorphism of a visited order, skip it.
        if visited.contains(&perm) {
            continue;
        }

        // Mark every permutation equivalent to this one (perm composed with
        // each automorphism) as visited, so we don't redo them later.
        for a in &ref_automorphisms {
            let equiv: Vec<usize> = (0..n).map(|i| a[perm[i]]).collect();
            visited.insert(equiv);
        }


        // Shape calculation
        let h = correlation_matrix(problem, &reordered);
        let (_, a_i) = optimal_rotation(h);
        let s = shape_measure(&a_i, n);

        if s < best_s {
            best_s = s;
            best_perm = perm;
        }
    }
    (best_s, best_perm)
}