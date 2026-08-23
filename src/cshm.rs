use itertools::Itertools;
use nalgebra::{Matrix3, Vector3};

pub fn center_and_normalise(points: &mut [[f64; 3]]) {
    let n = points.len() as f64;
    let mut centroid = [0.0; 3];

    for point in points.iter() {
        centroid[0] += point[0];
        centroid[1] += point[1];
        centroid[2] += point[2];
    }

    centroid[0] /= n;
    centroid[1] /= n;
    centroid[2] /= n;

    for point in points.iter_mut() {
        point[0] -= centroid[0];
        point[1] -= centroid[1];
        point[2] -= centroid[2];
    }

    let mut s2 = 0.0;
    for point in points.iter_mut() {
        s2 += point[0] * point[0] + point[1] * point[1] + point[2] * point[2]
    }

    // Normalization
    // Σ A²·|P_i|² = N
    // A = sqrt(N / Σ|P_i|²)

    let scale_factor = (n / s2).sqrt();

    for point in points.iter_mut() {
        point[0] *= scale_factor;
        point[1] *= scale_factor;
        point[2] *= scale_factor;
    }
}

pub fn correlation_matrix(reference: &[[f64; 3]], problem: &[[f64; 3]]) -> Matrix3<f64> {
    /// Expects centered and normalised points.

    let mut h = Matrix3::<f64>::zeros();

    for (p, q) in reference.iter().zip(problem.iter()) {
        let p_vec = Vector3::new(p[0], p[1], p[2]);
        let q_vec = Vector3::new(q[0], q[1], q[2]);

        h += p_vec * q_vec.transpose();
    }
    h
}

pub fn optimal_rotation(h: Matrix3<f64>) -> (Matrix3<f64>, Vector3<f64>)  {
    let svd = h.svd(true, true); // computes U and V^T
    let u = svd.u.unwrap();
    let mut v_t = svd.v_t.unwrap();
    let s = svd.singular_values;

    // The SVD algorithm can produce the reflected shape instead of the actual one.
    // To convert back to the desired rotation, flip one of the rows of the v_t matrix.
    if v_t.transpose().determinant() * u.transpose().determinant() < 0.0 {
        // TODO Flip one the rows so Det becomes positive.
        v_t.set_row(2, &(-v_t.row(2)));
    };
    (v_t.transpose() * u.transpose(), s)
}

pub fn shape_measure(singular_values: &Vector3<f64>, n: usize) -> f64 {
    let a: f64 = singular_values.iter().sum();
    (1.0 - a*a/(n as f64 * n as f64)) * 100.0
}

/// Brute-force search all n! automorphisms of a Reference Shape. Expects centered and normalised coordinates to work.
pub fn naive_find_automorphism(points: &[[f64; 3]]) -> Vec<Vec<usize>> {
    // Expects centered and normalised coordinates.
    let n = points.len();
    let mut automorphisms: Vec<Vec<usize>> = Vec::new();

    for perm in (0..n).permutations(n) {
        // Build the reordered point set using this permutation:
        let reordered: Vec<[f64;3]> = perm.iter().map(|&p| points[p]).collect();

        let h = correlation_matrix(points, &reordered);
        let (_, a_i) = optimal_rotation(h);
        let s = shape_measure(&a_i, n);

        if s.abs() < 1e-6 {
            automorphisms.push(perm);
        }
    }
    automorphisms
}

/// Finds the best permutation but prunes by automorphisms using naive_find_automorphism.
pub fn best_permutation_brute_force(
    reference: &[[f64; 3]],
    problem: &[[f64; 3]],
) -> (f64, Vec<usize>) {
    let n = problem.len();
    debug_assert_eq!(n, reference.len());

    let mut best_s = f64::INFINITY;
    let mut best_perm: Vec<usize> = Vec::new();

    let ref_automorphisms: Vec<Vec<usize>> = naive_find_automorphism(reference);
    let mut visited: std::collections::HashSet<Vec<usize>> = std::collections::HashSet::new();

    for perm in (0..n).permutations(n) {
        // Build the reordered point set using this permutation:
        let reordered: Vec<[f64;3]> = perm.iter().map(|&p| reference[p]).collect();

        // If this permutation of the reference shape is an automorphism of a visited order, skip it.
        if visited.contains(&perm) {
            continue;
        }

        // Mark every permutation equivalent to this one (perm composed with
        // each automorphism) as visited, so we don't redo them later.
        for a in &ref_automorphisms {
            let equiv: Vec<usize> = a.iter().map(|&i| perm[i]).collect();
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


mod tests {
    use crate::data::standard_shapes;
    use crate::shapes::ReferenceShape;
    use super::*;

    #[test]
    fn centre_and_normalises_correctly() {
        let mut points = Vec::from(
            [[1.0, 0.0, 0.0], // Regular octahedron centered at 1 0 0
                [2.0, 0.0, 1.0],
                [1.0, 1.0, 0.0],
                [1.0, 0.0, 1.0],
                [-0.0, 0.0, 1.0],
                [1.0, -1.0, 0.0],
                [1.0, 0.0, -1.0],
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
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 1.0],
            [1.0, 1.0, 0.0],
            [1.0, 0.0, 1.0],
            [-0.0, 0.0, 1.0],
        ];

        let know_rotataion = Matrix3::new(
            0.0, -1.0, 0.0,
            1.0,  0.0, 0.0,
            0.0,  0.0, 1.0,
        );

        let rotated: Vec<[f64; 3]> = reference.iter().map(|p| {
            let v = Vector3::new(p[0], p[1], p[2]);
            let r = know_rotataion * v;
            [r[0], r[1], r[2]]
        }).collect();

        let h = correlation_matrix(&reference, &rotated);
        let (recovered_r, _) = optimal_rotation(h);


        recovered_r.iter().zip(know_rotataion.iter()).for_each(|(r1, r2)| {assert!((r1-r2).abs() < 1e-10)});

    }

    #[test]
    fn perfect_gives_zero() {
        let mut points = [
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];

        center_and_normalise(&mut points);

        let h = correlation_matrix(&points, &points);
        let (_, s) = optimal_rotation(h);
        let s = shape_measure(&s, points.len());
        println!("s is: {}", s);
        assert!(s.abs() < 1e-10);
    }

    #[test]
    fn water_matches_result_from_shape21 () {
        // Coords from problem taken from water dataset from Olex2
        let mut problem = [
            [0.0,       8.0648,     0.0],       // Central Mn
            [-1.332417, 6.56007,   -1.099508],  // N2'
            [1.3324,    9.5695,     1.0995],    // N2
            [0.521462,  6.437162,   1.333904],  // O4
            [-0.5215,   9.6924,    -1.3339],    // O4'
            [1.619207,  7.429778,  -1.386425],  // O5
            [-1.6192,   8.6998,     1.3864],    // O5'
        ];


        // This list of atoms is manually set in order to correspond
        // to the correct assignation of point pairs.
        let mut reference = [
            [0.0, 0.0, 0.0], // Central atom first.
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ];

        center_and_normalise(&mut reference);
        center_and_normalise(&mut problem);

        let h = correlation_matrix(&reference, &problem);
        let (_, s) = optimal_rotation(h);
        let s = shape_measure(&s, problem.len());
        println!("s is: {}", s);
        assert!((s - 0.18).abs() < 1e-2);
    }

    fn square() -> Vec<[f64; 3]> {
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ]
    }

    #[test]
    fn identity_problem_gives_zero_score() {
        let mut reference = square();
        let mut problem = square();

        center_and_normalise(&mut reference);
        center_and_normalise(&mut problem);

        let (best_s, best_perm) = best_permutation_brute_force(&reference, &problem);
        assert!(best_s.abs() < 1e-10, "expected near-zero shape measure, got {best_s}");
        assert_eq!(best_perm.len(), 4)
    }

    #[test]
    fn permuted_problem_gives_zero_score() {
        let mut reference = square();
        // problem is reference shuffled — brute force should still find perm
        // that undoes the shuffle and scores ~0.
        let mut problem = vec![reference[2], reference[0], reference[3], reference[1]];

        center_and_normalise(&mut reference);
        center_and_normalise(&mut problem);

        let (best_s, _best_perm) = best_permutation_brute_force(&reference, &problem);
        assert!(best_s.abs() < 1e-10, "expected near-zero shape measure, got {best_s}");
    }

    #[test]
    fn automorphism_does_not_change_best_score() {
        let reference = square();
        let problem = [
            [0.1, 0.0, 0.0],
            [1.0, 0.1, 0.0],
            [0.0, 0.9, 0.0],
            [0.95, 1.0, 0.0],
        ]; // Noisy square.
        let (best_s, _) = best_permutation_brute_force(&reference, &problem);

        // Brute force all n! permutations.
        let n = reference.len();
        let mut true_best_s = f64::INFINITY;

        for perm in (0..n).permutations(n) {
            let reordered: Vec<[f64; 3]> = perm.iter().map(|&p| reference[p]).collect();
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
            [-4.09549, 5.2296, -3.52751],
            [-3.62959, 2.88708, -3.07315],
            [-2.30454, 4.56031, -1.99152],
            [-4.37536, 3.8733, -3.71495],
            [-2.57894, 3.214, -2.20586],
            [-3.05614, 5.527, -2.6526],
        ];

        let mut reference_hexagon = [ // Regular hexagon reference from standard_shapes.rs
            [1.0, 0.0, 0.0],
            [0.5, 0.8660254, 0.0],
            [-0.5, 0.8660254, 0.0],
            [-1.0, 0.0, 0.0],
            [-0.5, -0.8660254, 0.0],
            [0.5, -0.8660254, 0.0]
        ];

        let mut reference_octahedron = [
            [0.0, 0.0, -1.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
            [-0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0]
        ];


        center_and_normalise(&mut reference_hexagon);
        center_and_normalise(&mut reference_octahedron);
        center_and_normalise(&mut problem);

        let shape21_hex = 0.035;
        let shape21_oc = 33.342;
        let (best_s_hex, _) = best_permutation_brute_force(&reference_hexagon, &problem);
        let (best_s_oc, _) = best_permutation_brute_force(&reference_octahedron, &problem);

        assert!((best_s_hex - shape21_hex).abs() < 1e-3,
                "Result doesn't match SHAPE 2.1. Expected {shape21_hex}, got {best_s_hex}");
        assert!((best_s_oc - shape21_oc).abs() < 1e-3,
                "Result doesn't match SHAPE 2.1. Expected {shape21_oc}, got {best_s_oc}");
    }
}