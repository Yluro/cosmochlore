/// GEOMETRY RELATED FUNCTIONS
use nalgebra::{Matrix3, Vector3};


/// Puts the shape centroid in the origin of coordinates and calculate the A
/// normalisation factor for a cloud of points A²Σ|P_i|² = N
pub fn center_and_normalise(points: &mut [Vector3<f64>]) -> (Vector3<f64>, f64)  {
    let n = points.len() as f64;
    let mut centroid: Vector3<f64> = Vector3::new(0.0, 0.0, 0.0);

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
    (centroid, scale_factor)

}


/// Calculates the correlation matrix H of two given sets of points.
/// Returns H = Sum (P^T x Q).
pub fn correlation_matrix(reference: &[Vector3<f64>], problem: &[Vector3<f64>]) -> Matrix3<f64> {
    // Expects centered and normalised points.
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
    let mut v_t = svd.v_t.unwrap();
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