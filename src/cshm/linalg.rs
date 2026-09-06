//! LINEAR ALGEBRA RELATED FUNCTIONS
use nalgebra::{Matrix3, Vector3};
use std::f64::consts::PI;

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

/// Eigenvalues of M = H^T*H of a symmetric 3x3 matrix via the trigonometric solution
/// of its characteristic cubic (Smith, 1961).
#[inline]
fn symmetric_eigenvalues_3x3(m: &Matrix3<f64>) -> [f64; 3] {
    let (a, b, c) = (m[(0, 0)], m[(1, 1)], m[(2, 2)]); // Diagonal values
    let (d, e, f) = (m[(0, 1)], m[(0, 2)], m[(1, 2)]); // Off diagonal

    let p1 = d * d + e * e + f * f;
    if p1 <= 0.0 {
        return [a, b, c]; // m is already diagonal.
    }

    let q = (a + b + c) / 3.0; // mean of the diagonal == trace / 3
    let (da, db, dc) = (a - q, b - q, c - q);
    let p2 = da * da + db * db + dc * dc + 2.0 * p1;
    let p = (p2 / 6.0).sqrt();

    // r = det((m - q*I) / p) / 2, clamped for floating-point safety before acos.
    let inv = 1.0 / p;
    let (ba, bb, bc) = (da * inv, db * inv, dc * inv);
    let (bd, be, bf) = (d * inv, e * inv, f * inv);
    let det_b = ba * (bb * bc - bf * bf) - bd * (bd * bc - bf * be) + be * (bd * bf - bb * be);
    let r = (det_b / 2.0).clamp(-1.0, 1.0);

    let phi = r.acos() / 3.0;
    let eig1 = q + 2.0 * p * phi.cos();
    let eig3 = q + 2.0 * p * (phi + 2.0 * PI / 3.0).cos();
    let eig2 = 3.0 * q - eig1 - eig3; // trace = eig1 + eig2 + eig3 = 3q

    [eig1, eig2, eig3]
}

/// Computes the nuclear norm (sum of singular values) of a correlation matrix H,
/// via the closed-form eigenvalues of M = H^T * H.
pub(crate) fn nuclear_norm(h: Matrix3<f64>) -> f64 {
    let m = h.transpose() * h;
    symmetric_eigenvalues_3x3(&m)
        .iter()
        .map(|&v| v.max(0.0).sqrt())
        .sum()
}