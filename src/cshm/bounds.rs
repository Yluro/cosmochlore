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


/// Precomputes the Euclidean norm of every point once, so the branch-and-bound
/// loop can look these up instead of recomputing `.norm()` at every node.
pub(crate) fn precompute_norms(points: &[Vector3<f64>]) -> Vec<f64> {
    points.iter().map(|p| p.norm()).collect()
}

/// Precomputes suffix sums of `norms`: `suffix[k] = norms[k..].sum()`, with
/// `suffix[norms.len()] = 0.0`.
pub(crate) fn precompute_suffix_sums(norms: &[f64]) -> Vec<f64> {
    let mut suffix = vec![0.0; norms.len() + 1];
    for k in (0..norms.len()).rev() {
        suffix[k] = suffix[k + 1] + norms[k];
    }
    suffix
}

pub(crate) fn max_unassigned_norm(ref_norms: &[f64], assigned: &[bool]) -> f64 {
    ref_norms
        .iter()
        .zip(assigned.iter())
        .filter(|&(_, is_assigned)| !is_assigned )
        .map(|(&norm, _)| norm)
        .fold(0.0, f64::max)
}



