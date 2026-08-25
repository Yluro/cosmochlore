use std::collections::HashSet;
use itertools::Itertools;
use nalgebra::{Matrix3, Vector3};
use std::time::Instant;

pub fn center_and_normalise(points: &mut [Vector3<f64>]) {
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

/// Calculates the correlation matrix H of two given sets of points.
/// H = Sum (P^T x Q)
pub fn correlation_matrix(reference: &[Vector3<f64>], problem: &[Vector3<f64>]) -> Matrix3<f64> {
    // Expects centered and normalised points.
    let mut h = Matrix3::<f64>::zeros();
    for (p, q) in reference.iter().zip(problem.iter()) {
        h += p * q.transpose();
    }
    h
}

pub fn optimal_rotation(h: Matrix3<f64>) -> (Matrix3<f64>, Vector3<f64>)  {
    let svd = h.svd(true, true); // computes U and V^T
    let u = svd.u.unwrap();
    let mut v_t = svd.v_t.unwrap();
    let a_i = svd.singular_values;

    // The SVD algorithm can produce the reflected shape instead of the actual one.
    // To convert back to the desired rotation, flip one of the rows of the v_t matrix.
    if v_t.transpose().determinant() * u.transpose().determinant() < 0.0 {
        // TODO Flip one the rows so Det becomes positive.
        v_t.set_row(2, &(-v_t.row(2)));
    };
    (v_t.transpose() * u.transpose(), a_i)
}

pub fn shape_measure(singular_values: &Vector3<f64>, n: usize) -> f64 {
    let a: f64 = singular_values.iter().sum();
    (1.0 - a*a/(n as f64 * n as f64)) * 100.0
}


/// Brute-force search all n! automorphisms of a Reference Shape. Expects centered and normalised coordinates to work.
fn naive_find_automorphism(points: &[Vector3<f64>]) -> Vec<Vec<usize>> {
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

/// Finds a Reference Shape's Automorphisms.
/// Expects centered and normalised Shapes.
fn find_automorphisms(reference: &[Vector3<f64>]) -> Vec<Vec<usize>> {
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
        if assigned[ref_idx] { // Skip the currently asigned points in the permutation
            continue;
        }
        // Choosing the next point starts here:

        //
        *h_partial += hi[ref_idx][pos];
        assigned[ref_idx] = true;

        let a_partial = singular_value_sum(*h_partial);
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





/// Finds the best permutation by iterating all over the n! permutation list
/// but prunes by automorphisms using naive_find_automorphism.
pub fn best_permutation_brute_force(
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

        let a_partial = singular_value_sum(*h_partial); // Calculate the partial SV sum,
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

fn max_unassigned_norm(reference: &[Vector3<f64>], assigned: &[bool]) -> f64 {
    reference
        .iter()
        .zip(assigned.iter())
        .filter(|&(_, is_assigned)| !is_assigned )
        .map(|(point, _)| {
            point.norm()
        })
        .fold(0.0, f64::max)
}

fn unassigned_norms_sum(problem_remaining: &[Vector3<f64>]) -> f64 {
    problem_remaining.iter().map(|p| p.norm()).sum()
}


/// Precomputes the correlation block for the single point analysis.
/// Returns a vector of vectors: hi[ref_id][pro_id] that stores the matrix.
fn precompute_correlation_blocks(
    reference: &[Vector3<f64>],
    problem: &[Vector3<f64>],
) -> Vec<Vec<Matrix3<f64>>> {
    let n = reference.len();
    debug_assert_eq!(n, problem.len());

    let mut hi: Vec<Vec<Matrix3<f64>>> = vec![vec![Matrix3::zeros(); n]; n];

    for ref_idx in 0..n {
        for problem_idx in 0..n {
            hi[ref_idx][problem_idx] = reference[ref_idx] * problem[problem_idx].transpose();
        }
    }
    hi
}


/// Reimplementation of cshm.f90's suml routine.
/// Calculates the singular value sum of a given correlation matrix using
/// nalgebra's .symmetric_eigen() method. Masks eigenvalues > 0 to 0.0 to avoid floating point errors.
fn singular_value_sum(h: Matrix3<f64>) -> f64 {
    let m = h.transpose() * h;
    let eig = m.symmetric_eigen();

    eig.eigenvalues.iter().map(|&v| v.max(0.0).sqrt()).sum()
}


mod tests {
    use super::*;

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
        let mut bnb_autom = find_automorphisms(&reference);
        let mut all_autom = naive_find_automorphism(&reference);
        let perms = bnb_autom.len();
        bnb_autom.sort();
        all_autom.sort();
        assert_eq!(perms, 48);    // Finds all the 48 symmetry elements of the octahedron.
        assert_eq!(bnb_autom, all_autom, "Expected {:?}, found: {:?}", all_autom.len(), bnb_autom.len());
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
        let (s_bnb, _) = best_permutation_branch_and_bound(&mut reference, &mut problem);

        assert!((s_bf - s_bnb).abs() < 1e-10, "true optimal value was pruned. Expected {s_bf}, found {s_bnb}.")
    }

    #[test]
    fn bnb_mathes_bf_matches_shape21_hard() {
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
        let (s_bnb, _) = best_permutation_branch_and_bound(&mut reference, &mut problem);

        assert!((s_bf - s_bnb).abs() < 1e-10, "true optimal value was pruned. Expected {s_bf}, found {s_bnb}.");
        assert!((s_bnb - shape21_result).abs() < 1e-3, "Calculation doesn't match SHAPE 2.1 output: Expected {shape21_result}, found {s_bnb}.")
    }

    #[test]
    fn bnb_is_faster_than_bf() {
        let mut problem = [ // From Eu7 dataset.
            Vector3::new(5.44844, 5.38278, 7.85016),
            Vector3::new(6.16348, 8.57381, 7.75539),
            Vector3::new(4.83935, 2.1473, 8.09115),
            Vector3::new(2.93143, 5.57194, 8.26729),
            Vector3::new(7.99552, 5.18822, 8.38373),
            Vector3::new(5.22877, 5.61787, 10.3971),
            Vector3::new(4.39975, 6.62039, 5.71296),
            Vector3::new(6.48067, 4.00061, 5.84709),
        ];

        let mut reference = [ // Capped trigonal prism
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
        let (s_bnb, _) = best_permutation_branch_and_bound(&mut ref_bnb, &mut prob_bnb);
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
        let (s, _) = best_permutation_branch_and_bound(&mut reference, &mut problem);
        let time = start.elapsed();
        println!("Find Optimal Permutation Time: {:?}", time);

        let shape21_result = 7.288;

        assert!((s - shape21_result).abs() < 1e-3, "Calculation doesn't match SHAPE 2.1. Expected: {shape21_result}, found {s}.");
    }
}