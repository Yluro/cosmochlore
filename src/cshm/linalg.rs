/// GEOMETRY RELATED FUNCTIONS
use nalgebra::{Matrix3, Vector3};


/// Calculates the correlation matrix H of two given sets of points.
/// Returns H = Sum (P^T x Q).
pub fn correlation_matrix(reference: &[Vector3<f64>], problem: &[Vector3<f64>]) -> Matrix3<f64> {
    // Expects centered and normalized points.
    let mut h = Matrix3::<f64>::zeros();
    for (p, q) in reference.iter().zip(problem.iter()) {
        h += p * q.transpose();
    }
    h
}

/// Finds the optimal rotation given a correlation matrix H using
/// Kabsch's SVD algorithm.
/// Returns the Orientation Matrix and a vector of eigenvalues.
pub fn optimal_rotation(h: Matrix3<f64>) -> (Matrix3<f64>, Vector3<f64>)  {
    let svd = h.svd(true, true); // computes U and V^T
    let u = svd.u.unwrap();
    let v_t = svd.v_t.unwrap();
    let a_i = svd.singular_values;

    (v_t.transpose() * u.transpose(), a_i)
}

/// Computes the shape measure given a list of eigenvalues from the SVD.
pub fn shape_measure(singular_values: &Vector3<f64>, n: usize) -> f64 {
    let a: f64 = singular_values.iter().sum();
    (1.0 - a*a/(n as f64 * n as f64)) * 100.0
}

/// Computes the singular values of the SVD using M = H^T * H.
///Reimplementation of cshm.f90's suml routine.
/// Calculates the singular value sum of a given correlation matrix using
/// nalgebra's .symmetric_eigen() method. Masks eigenvalues > 0 to 0.0 to avoid floating point errors.
pub(crate) fn singular_values(h: Matrix3<f64>) -> Vec<f64> {
    let m = h.transpose() * h;
    let eig = m.symmetric_eigen();

    eig.eigenvalues.iter().map(|&v| v.max(0.0).sqrt()).collect()
}