/// General geometric tools.

use nalgebra::{Vector3, Matrix3};

/// Puts the shape centroid in the origin of coordinates and calculate the A
/// normalization factor for a cloud of points A²Σ|P_i|² = N
///
/// Mutates the input so its normalized and its centroid sits at the origin.
/// Returns the
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

/// Places the centroid of the points at the origin.
pub fn center_by_centroid(points: &mut [Vector3<f64>]) -> Vector3<f64> {

    let n = points.len() as f64;
    let centroid = points.iter().sum::<Vector3<f64>>() / n;

    for p in points.iter_mut() {
        *p -= centroid;
    }

    centroid
}

/// Places the first point in the list at the origin.
pub fn center_by_first_point(points: &mut [Vector3<f64>]) -> Vector3<f64> {
    let p0 = points[0];
    for p in points.iter_mut() {
        *p -= p0;
    }
    p0
}

/// Places the centre of a structure in a given point.
pub fn center_by_coordinate(points: &mut [Vector3<f64>], centre: Vector3<f64>) -> Vector3<f64> {
    for p in points.iter_mut() {
        *p -= centre;
    }
    centre
}

/// Assumes centered points.
pub(crate) fn normalise(points: &mut [Vector3<f64>]) -> f64 {

    let n = points.len() as f64;
    let sq = points.iter().map(|v| v.norm_squared()).sum::<f64>();

    let scale_factor = (n / sq).sqrt();

    for p in points.iter_mut() {
        *p *= scale_factor;
    }

    scale_factor
}

/// Returns the rotation matrix given the axis of rotation and the angle in degrees.
pub fn rotation_matrix(axis: Vector3<f64>, angle: f64 ) -> Matrix3<f64> {

    let axis = axis.normalize();

    let k: Matrix3<f64> = Matrix3::new(
        0.0,    -axis[2], axis[1],
        axis[2],0.0, -axis[0],
        -axis[1],axis[0], 0.0,
    );
    let rads = angle.to_radians();

    let rot: Matrix3<f64> = Matrix3::identity() + rads.sin() * k + (1.0 - rads.cos()) * k * k;
    rot.transpose()
}