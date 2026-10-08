//! General geometric tools.

use itertools::Itertools;
use nalgebra::{Matrix3, Vector3};

/// The tetrahedral angle `θ = arccos(-1/3)` in degrees.
pub const TETRAHEDRAL_ANGLE: f64 = 109.471_220_634_490_69;

/// Puts the shape centroid in the origin of coordinates and calculate the A
/// normalization factor for a cloud of points A²Σ|P_i|² = N
///
/// Mutates the input so its normalized and its centroid sits at the origin.
/// Returns the
pub fn center_and_normalise(points: &mut [Vector3<f64>]) -> (Vector3<f64>, f64) {
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

/// Returns the rotation matrix given the axis of rotation and the angle in degrees
/// using the Rodriges formula.
pub fn rotation_matrix(axis: Vector3<f64>, angle: f64) -> Matrix3<f64> {
    let axis = axis.normalize();

    let k: Matrix3<f64> = Matrix3::new(
        0.0, -axis.z, axis.y, axis.z, 0.0, -axis.x, -axis.y, axis.x, 0.0,
    );
    let rads = angle.to_radians();

    let rot: Matrix3<f64> = Matrix3::identity() + rads.sin() * k + (1.0 - rads.cos()) * k * k;
    rot.transpose()
}

/// Rodrigues' rotation formula from a rotation vector `v`: direction is
/// the axis, `|v|` (radians) is the angle. `v = 0` maps to identity.
pub fn rotation_matrix_from_vector(v: Vector3<f64>) -> Matrix3<f64> {
    let angle = v.norm();
    if angle < 1e-10 {
        return Matrix3::identity();
    }
    let axis = v / angle;
    let k = Matrix3::new(
        0.0, -axis.z, axis.y, axis.z, 0.0, -axis.x, -axis.y, axis.x, 0.0,
    );
    Matrix3::identity() + angle.sin() * k + (1.0 - angle.cos()) * (k * k)
}

/// Compute the angle in degrees between two vectors, in `[0, 180]`.
///
/// NaN if either vector has zero length.
pub fn angle_btw_vectors(v1: Vector3<f64>, v2: Vector3<f64>) -> f64 {
    v1.normalize()
        .dot(&v2.normalize())
        .clamp(-1.0, 1.0)
        .acos()
        .to_degrees()
}

/// Every angle in degrees between two of the vectors, in the order of the index pairs
/// `(0, 1), (0, 2), ..., (0, n-1), (1, 2), ..., (n-2, n-1)`.
pub fn pairwise_angles(vectors: &[Vector3<f64>]) -> Vec<f64> {
    vectors
        .iter()
        .array_combinations()
        .map(|[vi, vj]| angle_btw_vectors(*vi, *vj))
        .collect()
}

/// Compute the angle in degrees between two vectors, signed by whether `v1`, `v2` and `direct`
/// form a right- or left-handed set (i.e. whether going from `v1` to `v2` is CW or CCW when
/// viewed from the `direct` side).
pub fn angle_sign(v1: Vector3<f64>, v2: Vector3<f64>, direct: Vector3<f64>) -> f64 {
    let v1 = v1.normalize();
    let v2 = v2.normalize();
    let angle = v1.dot(&v2).clamp(-1.0, 1.0).acos().to_degrees();

    let matrix = Matrix3::new(
        v1.x, v1.y, v1.z, v2.x, v2.y, v2.z, direct.x, direct.y, direct.z,
    );

    if matrix.determinant() < 0.0 {
        -angle
    } else {
        angle
    }
}

/// The index of the vector in `vectors` that has the largest angle with `from`, that is, the one
/// most nearly anti-parallel to it (the first one, when several tie).
pub fn most_trans_index(from: Vector3<f64>, vectors: &[Vector3<f64>]) -> usize {
    let mut best_angle = 0.0;
    let mut best_index = 0;
    for (n, v) in vectors.iter().enumerate() {
        let angle = angle_btw_vectors(from, *v);
        if angle > best_angle {
            best_angle = angle;
            best_index = n;
        }
    }
    best_index
}

/// Coefficients `(a, b, c, d)` of the plane `a*x + b*y + c*z = d` through three points.
pub fn find_eq_of_plane(x: Vector3<f64>, y: Vector3<f64>, z: Vector3<f64>) -> (f64, f64, f64, f64) {
    let cross = (z - x).cross(&(y - x));
    let d = cross.dot(&z);
    (cross.x, cross.y, cross.z, d)
}

/// Orthogonal projection of point `p` onto the plane `a*x + b*y + c*z = d`.
pub fn project_atom_onto_plane(p: Vector3<f64>, a: f64, b: f64, c: f64, d: f64) -> Vector3<f64> {
    let plane = Vector3::new(a, b, c);
    let lambda = (d - plane.dot(&p)) / plane.dot(&plane);
    p + lambda * plane
}

/// Six times the signed volume of the tetrahedron with vertices `a`, `b`, `c`, `d`.
pub fn tetrahedron_volume_x6(
    a: Vector3<f64>,
    b: Vector3<f64>,
    c: Vector3<f64>,
    d: Vector3<f64>,
) -> f64 {
    (a - d).dot(&(b - d).cross(&(c - d)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn tetrahedral_angle_constant_is_arccos_minus_one_third() {
        let exact = (-1.0_f64 / 3.0).acos().to_degrees();

        assert!((TETRAHEDRAL_ANGLE - exact).abs() < 1e-12);
    }

    #[test]
    fn angle_btw_vectors_gives_the_angle_in_degrees() {
        let x = Vector3::new(2.0, 0.0, 0.0);

        assert!((angle_btw_vectors(x, Vector3::new(0.0, 3.0, 0.0)) - 90.0).abs() < 1e-12);
        assert!((angle_btw_vectors(x, Vector3::new(1.0, 1.0, 0.0)) - 45.0).abs() < 1e-12);
        assert!((angle_btw_vectors(x, Vector3::new(-5.0, 0.0, 0.0)) - 180.0).abs() < 1e-6);
        // Parallel vectors: acos is ill-conditioned at 1, so only demand a tiny angle.
        assert!(angle_btw_vectors(x, Vector3::new(7.0, 0.0, 0.0)) < 1e-5);
    }

    #[test]
    fn angle_btw_vectors_is_nan_for_a_zero_length_vector() {
        assert!(angle_btw_vectors(Vector3::zeros(), Vector3::x()).is_nan());
    }

    #[test]
    fn pairwise_angles_lists_every_pair_in_index_order() {
        let vectors = [Vector3::x(), Vector3::y(), -Vector3::x(), Vector3::z()];
        let angles = pairwise_angles(&vectors);

        // (0,1) (0,2) (0,3) (1,2) (1,3) (2,3)
        let expected = [90.0, 180.0, 90.0, 90.0, 90.0, 90.0];
        assert_eq!(angles.len(), expected.len());
        for (angle, expected) in angles.iter().zip(expected) {
            assert!((angle - expected).abs() < 1e-6);
        }
    }

    #[test]
    fn angle_sign_follows_the_handedness_of_the_viewing_direction() {
        let (x, y, z) = (Vector3::x(), Vector3::y(), Vector3::z());

        assert!((angle_sign(x, y, z) - 90.0).abs() < 1e-12);
        assert!((angle_sign(x, y, -z) + 90.0).abs() < 1e-12);
        assert!((angle_sign(y, x, z) + 90.0).abs() < 1e-12);
    }

    #[test]
    fn most_trans_index_picks_the_most_opposite_vector() {
        let vectors = [
            Vector3::x(),
            Vector3::y(),
            Vector3::new(-1.0, 0.1, 0.0),
            Vector3::z(),
        ];

        assert_eq!(most_trans_index(Vector3::x(), &vectors), 2);
        assert_eq!(most_trans_index(Vector3::new(0.0, -1.0, 0.0), &vectors), 1);
    }

    #[test]
    fn plane_through_three_points_and_projection_onto_it() {
        let (a, b, c, d) = find_eq_of_plane(Vector3::x(), Vector3::y(), Vector3::z());

        // The plane x + y + z = 1, up to the scale of its normal.
        assert!((a - b).abs() < 1e-12 && (b - c).abs() < 1e-12);
        assert!((d / a - 1.0).abs() < 1e-12);

        let projected = project_atom_onto_plane(Vector3::zeros(), a, b, c, d);
        let third = 1.0 / 3.0;
        assert!((projected - Vector3::new(third, third, third)).norm() < 1e-12);
        // A point already on the plane does not move.
        assert!((project_atom_onto_plane(Vector3::x(), a, b, c, d) - Vector3::x()).norm() < 1e-12);
    }

    #[test]
    fn tetrahedron_volume_x6_is_signed() {
        let (x, y, z, o) = (Vector3::x(), Vector3::y(), Vector3::z(), Vector3::zeros());

        // The corner of the unit cube: volume 1/6.
        assert!((tetrahedron_volume_x6(x, y, z, o) - 1.0).abs() < 1e-12);
        assert!((tetrahedron_volume_x6(y, x, z, o) + 1.0).abs() < 1e-12);
    }

    #[test]
    fn zero_returns_identity() {
        let v = Vector3::zeros();
        let rot_mat = rotation_matrix_from_vector(v);

        assert_eq!(
            rot_mat,
            Matrix3::identity(),
            "expected Identity, found {}",
            rot_mat
        );
    }

    #[test]
    fn pi_half_gives_quarter_rotation() {
        let v = Vector3::new(0.0, 0.0, PI / 2.0);
        let rot_mat = rotation_matrix_from_vector(v);

        let expected_rot = Matrix3::new(0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0);

        assert!(
            (expected_rot - rot_mat).abs().max() < 1e-10,
            "expected near-zero difference, found {}",
            (rot_mat - expected_rot).abs().max()
        );
    }
}
