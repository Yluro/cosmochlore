use nalgebra::{Normed, Vector3};



/// Shape deviation for two given shapes. It returns the
///
/// Assumes centered and normalized points.
/// Assumes correct permutation of points.
pub(crate) fn shape_deviation(reference: &[Vector3<f64>], problem: &[Vector3<f64>]) -> f64 {

    let n = reference.len() as f64;
    let problem_centroid = problem.iter().sum::<Vector3<f64>>() / n;

    let denominator = problem.iter()
        .map(|p| (p - problem_centroid).norm_squared())
        .sum::<f64>();

    let numerator = problem.iter().zip(reference.iter())
        .map(
            |(q , p)| (q - p).norm_squared()
        ).sum::<f64>();


    100.0 * numerator / denominator

}