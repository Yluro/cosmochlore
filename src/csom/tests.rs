use nalgebra::Vector3;
use crate::csom::dev::best_permutation;
use crate::geometry::{normalise, rotation_matrix};

#[test]
fn can_rotate_and_assign() {

    let mut shape: Vec<Vector3<f64>> = Vec::from( // Trigonal pyramid centered at origo.
        [Vector3::zeros(),           // Point 0
        Vector3::new(0.0, 0.0, -1.0), // 1
        Vector3::new(1.0, 0.0, 0.0), // 2
        Vector3::new(-0.5, 0.8660254, 0.0), // 3
        Vector3::new(-0.5, -0.8660254, 0.0), // 4
        Vector3::new(0.0, 0.0, 1.), // 5
    ]);

    normalise(&mut shape);

    let axis = Vector3::new(0.0, 0.0, 1.0);
    let angle = 115.0; // Approximate C3 rotation.

    let rot_mat = rotation_matrix(axis, angle);

    let rotated: Vec<Vector3<f64>> = shape.iter().map(|v| rot_mat * *v).collect::<Vec<Vector3<f64>>>();

    let (_, perm) = best_permutation(&shape, rotated.as_slice());
    let expected_perm = vec![0, 1, 3, 4, 2, 5];

    assert_eq!(perm, expected_perm)


}