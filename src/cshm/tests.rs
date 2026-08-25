/// TESTS FOR CSHM MODULE
use itertools::Itertools;
use nalgebra::{Vector3, Matrix3};
use std::time::Instant;

use crate::cshm::{geometry::*, automorphism::*, permutations::*, test_utils::*};

#[test]
fn centre_and_normalises_correctly() {
    let mut points = Vec::from(
        [Vector3::new(1.0, 0.0, 0.0), // Regular octahedron centered at 1 0 0
            Vector3::new(2.0, 0.0, 1.0),
            Vector3::new(1.0, 1.0, 0.0),
            Vector3::new(1.0, 0.0, 1.0),
            Vector3::new(-0.0, 0.0, 1.0),
            Vector3::new(1.0, -1.0, 0.0),
            Vector3::new(1.0, 0.0, -1.0),
        ]
    );

    center_and_normalise(&mut points);

    let mut centroid = [0.0; 3];

    for point in points.iter() {
        centroid[0] += point[0];
        centroid[1] += point[1];
        centroid[2] += point[2];
    }

    centroid[0] /= points.len() as f64;
    centroid[1] /= points.len() as f64;
    centroid[2] /= points.len() as f64;

    for point in centroid.iter() {assert!(point.abs() < 1e-10)}

    let mut s2: f64 = 0.0;
    for point in points.iter() {s2 += point[0] * point[0] + point[1] * point[1] + point[2] * point[2]};
    assert!((s2 - points.len() as f64).abs() < 1e-10);

}


#[test]
fn recovers_from_rotation() {
    let reference = [
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(2.0, 0.0, 1.0),
        Vector3::new(1.0, 1.0, 0.0),
        Vector3::new(1.0, 0.0, 1.0),
        Vector3::new(-0.0, 0.0, 1.0),
    ];

    let know_rotation = Matrix3::new(
        0.0, -1.0, 0.0,
        1.0,  0.0, 0.0,
        0.0,  0.0, 1.0,
    );

    let rotated = reference.iter().map(|&v| know_rotation * v).collect::<Vec<_>>();

    let h = correlation_matrix(&reference, &rotated);
    let (recovered_r, _) = optimal_rotation(h);


    recovered_r.iter().zip(know_rotation.iter()).for_each(|(r1, r2)| {assert!((r1-r2).abs() < 1e-10)});

}

#[test]
fn perfect_gives_zero() {
    let mut points = [
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(2.0, 0.0, 0.0),
        Vector3::new(1.0, 1.0, 0.0),
        Vector3::new(1.0, 1.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
    ];

    center_and_normalise(&mut points);

    let h = correlation_matrix(&points, &points);
    let (_, s) = optimal_rotation(h);
    let s = shape_measure(&s, points.len());
    println!("s is: {}", s);
    assert!(s.abs() < 1e-10);
}



fn octahedron() -> [Vector3<f64>; 7] {
    [
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(0.0, 0.0, -1.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(-1.0, 0.0, 0.0),
        Vector3::new(-0.0, -1.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
    ]
}
#[test]
fn bnb_automorphisms_matches_naive() {
    let mut reference = octahedron();

    center_and_normalise(&mut reference);
    let mut bnb_automorphisms = find_automorphisms(&reference);
    let mut all_automorphisms = naive_find_automorphism(&reference);
    let perms = bnb_automorphisms.len();
    bnb_automorphisms.sort();
    all_automorphisms.sort();
    assert_eq!(perms, 48);    // Finds all the 48 symmetry elements of the octahedron.
    assert_eq!(bnb_automorphisms, all_automorphisms, "Expected {:?}, found: {:?}", all_automorphisms.len(), bnb_automorphisms.len());
}

#[test]
fn water_matches_result_from_shape21 () {
    // Coords from problem taken from water dataset from Olex2
    let mut problem = [
        Vector3::new(0.0,       8.0648,     0.0),       // Central Mn
        Vector3::new(-1.332417, 6.56007,   -1.099508),  // N2'
        Vector3::new(1.3324,    9.5695,     1.0995),    // N2
        Vector3::new(0.521462,  6.437162,   1.333904),  // O4
        Vector3::new(-0.5215,   9.6924,    -1.3339),    // O4'
        Vector3::new(1.619207,  7.429778,  -1.386425),  // O5
        Vector3::new(-1.6192,   8.6998,     1.3864),    // O5'
    ];


    // This list of atoms is manually set in order to correspond
    // to the correct assignation of point pairs.
    let mut reference = [
        Vector3::new(0.0, 0.0, 0.0), // Central atom first.
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(-1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(0.0, -1.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
        Vector3::new(0.0, 0.0, -1.0),
    ];

    center_and_normalise(&mut reference);
    center_and_normalise(&mut problem);

    let h = correlation_matrix(&reference, &problem);
    let (_, s) = optimal_rotation(h);
    let s = shape_measure(&s, problem.len());
    println!("s is: {}", s);
    assert!((s - 0.18).abs() < 1e-2);
}

fn square() -> Vec<Vector3<f64>> {
    vec![
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(1.0, 1.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
    ]
}

#[test]
fn identity_problem_gives_zero_score() {
    let mut reference = square();
    let mut problem = square();


    let (best_s, best_perm) = best_permutation_brute_force(&mut reference, &mut problem);
    assert!(best_s.abs() < 1e-10, "expected near-zero shape measure, got {best_s}");
    assert_eq!(best_perm.len(), 4)
}

#[test]
fn permuted_problem_gives_zero_score() {
    let mut reference = square();
    // problem is reference shuffled — brute force should still find perm
    // that undoes the shuffle and scores ~0.
    let mut problem = vec![reference[2], reference[0], reference[3], reference[1]];

    let (best_s, _best_perm) = best_permutation_brute_force(&mut reference, &mut problem);
    assert!(best_s.abs() < 1e-10, "expected near-zero shape measure, got {best_s}");
}

#[test]
fn automorphism_does_not_change_best_score() {
    let mut reference = square();
    let mut problem = [
        Vector3::new(0.1, 0.0, 0.0),
        Vector3::new(1.0, 0.1, 0.0),
        Vector3::new(0.0, 0.9, 0.0),
        Vector3::new(0.95, 1.0, 0.0),
    ]; // Noisy square.
    let (best_s, _) = best_permutation_brute_force(&mut reference, &mut problem);

    // Brute force all n! permutations.
    let n = reference.len();
    let mut true_best_s = f64::INFINITY;

    for perm in (0..n).permutations(n) {
        let reordered: Vec<Vector3<f64>> = perm.iter().map(|&p| reference[p]).collect();
        let h = correlation_matrix(&problem, &reordered);
        let (_, a_i) = optimal_rotation(h);
        let s = shape_measure(&a_i, n);
        if s < true_best_s {
            true_best_s = s;
        }
    }
    assert!(
        (best_s - true_best_s).abs() < 1e-9,
        "dedup changed the result: got {best_s}, expected {true_best_s}"
    );
}

#[test]
/// Tests main permutation finding algorithm against a benzene ring found in the Olex2 'water' dataset.
fn unordered_coordinates_matches_shape21() {
    let mut problem = [
        Vector3::new(-4.09549, 5.2296, -3.52751),
        Vector3::new(-3.62959, 2.88708, -3.07315),
        Vector3::new(-2.30454, 4.56031, -1.99152),
        Vector3::new(-4.37536, 3.8733, -3.71495),
        Vector3::new(-2.57894, 3.214, -2.20586),
        Vector3::new(-3.05614, 5.527, -2.6526),
    ];

    let mut reference_hexagon = [ // Regular hexagon reference from standard_shapes.rs
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.5, 0.8660254, 0.0),
        Vector3::new(-0.5, 0.8660254, 0.0),
        Vector3::new(-1.0, 0.0, 0.0),
        Vector3::new(-0.5, -0.8660254, 0.0),
        Vector3::new(0.5, -0.8660254, 0.0),
    ];

    let mut reference_octahedron = [
        Vector3::new(0.0, 0.0, -1.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(-1.0, 0.0, 0.0),
        Vector3::new(-0.0, -1.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0)
    ];


    center_and_normalise(&mut reference_hexagon);
    center_and_normalise(&mut reference_octahedron);
    center_and_normalise(&mut problem);

    let shape21_hex = 0.035;
    let shape21_oc = 33.342;
    let (best_s_hex, _) = best_permutation_brute_force(&mut reference_hexagon, &mut problem);
    let (best_s_oc, _) = best_permutation_brute_force(&mut reference_octahedron, &mut problem);

    assert!((best_s_hex - shape21_hex).abs() < 1e-3,
            "Result doesn't match SHAPE 2.1. Expected {shape21_hex}, got {best_s_hex}");
    assert!((best_s_oc - shape21_oc).abs() < 1e-3,
            "Result doesn't match SHAPE 2.1. Expected {shape21_oc}, got {best_s_oc}");
}

#[test]
fn branch_and_bound_matches_brute_force_results() {
    let mut reference = square();
    let mut problem = [
        Vector3::new(0.1, 0.0, 0.0),
        Vector3::new(1.0, 0.1, 0.0),
        Vector3::new(0.0, 0.9, 0.0),
        Vector3::new(0.95, 1.0, 0.0),
    ]; // Noisy square.

    let (s_bf,_) = best_permutation_brute_force(&mut reference, &mut problem);
    let (s_bnb, _, _) = find_best_permutation(&mut reference, &mut problem);

    assert!((s_bf - s_bnb).abs() < 1e-10, "true optimal value was pruned. Expected {s_bf}, found {s_bnb}.")
}

#[test]
fn bnb_matches_bf_matches_shape21_hard() {
    let mut problem = [ // Coordinates from refined FeHS dataset.
        Vector3::new(4.92991, 10.3899, 12.9237),
        Vector3::new(6.20468, 10.6922, 14.7747),
        Vector3::new(5.00034, 8.51382, 13.9503),
        Vector3::new(6.80263, 10.0749, 11.7725),
        Vector3::new(5.47238, 12.2809, 12.1201),
        Vector3::new(3.30536, 11.3249, 13.862),
        Vector3::new(3.74731, 9.47729, 11.4687),
    ];

    let mut reference = octahedron();

    let shape21_result = 2.109;

    let (s_bf,_) = best_permutation_brute_force(&mut reference, &mut problem);
    let (s_bnb, _, _) = find_best_permutation(&mut reference, &mut problem);

    assert!((s_bf - s_bnb).abs() < 1e-10, "true optimal value was pruned. Expected {s_bf}, found {s_bnb}.");
    assert!((s_bnb - shape21_result).abs() < 1e-3, "Calculation doesn't match SHAPE 2.1 output: Expected {shape21_result}, found {s_bnb}.")
}

#[test]
fn bnb_is_faster_than_bf() {
    let problem = [ // From Eu7 dataset.
        Vector3::new(5.44844, 5.38278, 7.85016),
        Vector3::new(6.16348, 8.57381, 7.75539),
        Vector3::new(4.83935, 2.1473, 8.09115),
        Vector3::new(2.93143, 5.57194, 8.26729),
        Vector3::new(7.99552, 5.18822, 8.38373),
        Vector3::new(5.22877, 5.61787, 10.3971),
        Vector3::new(4.39975, 6.62039, 5.71296),
        Vector3::new(6.48067, 4.00061, 5.84709),
    ];

    let reference = [ // Capped trigonal prism
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
        Vector3::new(0.68689018, 0.68689018, 0.23741035),
        Vector3::new(-0.68689018, 0.68689018, 0.23741035),
        Vector3::new(0.68689018, -0.68689018, 0.23741035),
        Vector3::new(-0.68689018, -0.68689018, 0.23741035),
        Vector3::new(0.61748904, 0.0, -0.78657949),
        Vector3::new(-0.61748904, 0.0, -0.78657949),
    ];

    let mut ref_bf = reference.clone();
    let mut prob_bf = problem.clone();

    let shape21_result = 2.630;

    let start_bf = Instant::now();
    let (s_bf,_) = best_permutation_brute_force(&mut ref_bf, &mut prob_bf);
    let time_bf = start_bf.elapsed();

    let mut ref_bnb = reference.clone();
    let mut prob_bnb = problem.clone();

    let start_bnb = Instant::now();
    let (s_bnb, _, _) = find_best_permutation(&mut ref_bnb, &mut prob_bnb);
    let time_bnb = start_bnb.elapsed();

    assert!((s_bf - s_bnb).abs() < 1e-10, "true optimal value was pruned. Expected {s_bf}, found {s_bnb}.");
    assert!((s_bnb - shape21_result).abs() < 1e-3, "Calculation doesn't match SHAPE 2.1 output: Expected {shape21_result}, found {s_bnb}.");

    println!("Brute force: {:?}, Branch & bound: {:?}", time_bf, time_bnb);
    assert!(time_bnb < time_bf, "expected B&B to be faster: bf={:?}, bnb={:?}", time_bf, time_bnb);
}

#[test]
fn strain_12_point_test() {
    let mut problem = [ // 11 coordinate Lanthanum complex
        Vector3::new(4.95508, 11.3487, 7.16088),
        Vector3::new(5.71619, 10.8511, 9.51126),
        Vector3::new(2.62287, 11.6697, 8.08526),
        Vector3::new(5.03435, 13.813, 6.40506),
        Vector3::new(5.95314, 11.5481, 4.72175),
        Vector3::new(4.64156, 13.2917, 8.93972),
        Vector3::new(5.48014, 8.89859, 6.26662),
        Vector3::new(3.88241, 9.16413, 8.32509),
        Vector3::new(3.27692, 12.3607, 5.2059),
        Vector3::new(3.09533, 9.74408, 5.83297),
        Vector3::new(7.15277, 12.9425, 7.81747),
        Vector3::new(7.65921, 10.4947, 6.88486),
    ];

    let mut reference = [
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.5),
        Vector3::new(0.80901699, 0.58778525, -0.5),
        Vector3::new(0.30901699, 0.95105652, 0.5),
        Vector3::new(-0.30901699, 0.95105652, -0.5),
        Vector3::new(-0.80901699, 0.58778525, 0.5),
        Vector3::new(-1.0, 0.0, -0.5),
        Vector3::new(-0.80901699, -0.58778525, 0.5),
        Vector3::new(-0.30901699, -0.95105652, -0.5),
        Vector3::new(0.30901699, -0.95105652, 0.5),
        Vector3::new(0.80901699, -0.58778525, -0.5),
        Vector3::new(0.0, 0.0, -1.11803399),
    ];

    let start = Instant::now();
    let (s, _, _) = find_best_permutation(&mut reference, &mut problem);
    let time = start.elapsed();
    println!("Find Optimal Permutation Time: {:?}", time);

    let shape21_result = 7.288;

    assert!((s - shape21_result).abs() < 1e-3, "Calculation doesn't match SHAPE 2.1. Expected: {shape21_result}, found {s}.");
}
