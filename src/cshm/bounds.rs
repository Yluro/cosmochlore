use nalgebra::{Matrix3, Vector3};
/// HELPER FUNCTIONS TO DEFINE THE BOUNDS FOR THE BRANCHING ALGORITHMS

/// Precomputes the correlation block for the single point analysis.
/// Returns a vector of vectors: hi[ref_id][pro_id] that stores the matrix.
pub(crate) fn precompute_correlation_blocks(
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


pub(crate) fn max_unassigned_norm(reference: &[Vector3<f64>], assigned: &[bool]) -> f64 {
    reference
        .iter()
        .zip(assigned.iter())
        .filter(|&(_, is_assigned)| !is_assigned )
        .map(|(point, _)| {
            point.norm()
        })
        .fold(0.0, f64::max)
}

pub(crate) fn unassigned_norms_sum(problem_remaining: &[Vector3<f64>]) -> f64 {
    problem_remaining.iter().map(|p| p.norm()).sum()
}



