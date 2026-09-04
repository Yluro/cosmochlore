use nalgebra::{Matrix3, Vector3};
use crate::csom::dev::*;
use crate::csom::io::strip_label;
use crate::geometry::{center_by_centroid, center_by_first_point, normalise, rotation_matrix};


fn octahedron() -> Vec<Vector3<f64>> {
    Vec::from([ // Regular octahedron centered at origo
        Vector3::zeros(),           // Point 0
        Vector3::new(0.0, 0.0, -1.0), // 1
        Vector3::new(1.0, 0.0, 0.0), // 2
        Vector3::new(0.0, 1.0, 0.0), // 3
        Vector3::new(-1.0, 0.0, 0.0), // 4
        Vector3::new(0.0, -1.0, 0.0), // 5
        Vector3::new(0.0, 0.0, 1.0), //6
    ])
}

#[test]
fn can_rotate_and_assign() {

    let mut shape: Vec<Vector3<f64>> = octahedron();
    normalise(&mut shape);

    let axis = Vector3::new(0.0, 0.0, 1.0);
    let angle = 85.0; // Approximate C4 rotation.

    let rot_mat = rotation_matrix(axis, angle);

    let rotated: Vec<Vector3<f64>> = shape.iter().map(|v| rot_mat * *v).collect::<Vec<Vector3<f64>>>();

    let (_, perm) = best_permutation(&shape, rotated.as_slice());
    let expected_perm = vec![0, 1, 3, 4, 5, 2, 6];

    assert_eq!(perm, expected_perm)
}


#[test]
fn mismatched_gives_nonzero_dev() {

    let mut shape: Vec<Vector3<f64>> = octahedron();
    normalise(&mut shape);

    let labels = [
        "Fe".to_string(),
        "N12".to_string(),
        "Cl00A".to_string(),
        "Cl3".to_string(),
        "ClA".to_string(),
        "OA2".to_string(),
        "N2'".to_string(),
    ];

    let stripped = labels.iter().map(|l| strip_label(l).to_string() ).collect::<Vec<String>>();
    assert_eq!(vec!["Fe", "N", "Cl", "Cl", "Cl", "O", "N"], stripped);


    let axis = Vector3::new(0.0, 0.0, 1.0);
    let angle = 90.0; // Approximate C4 rotation.

    let rot_mat = rotation_matrix(axis, angle);

    let rotated: Vec<Vector3<f64>> = shape.iter().map(|v| rot_mat * *v).collect::<Vec<Vector3<f64>>>();

    let (a, b) = best_permutation_multiple_atoms(&shape, rotated.as_slice(), &stripped);

    // Sanity: same length, all atoms retained.
    assert_eq!(a.len(), shape.len());
    assert_eq!(b.len(), shape.len());

    let sds = sds_dev(&a, &b);
    assert!(sds > 1e-6, "expected nonzero deviation due to mismatched Cl group, got {sds}");
}

#[test]
fn matches_expected_sds_after_permutation() {
    let mut square = vec![
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(1.0, 1.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
    ];

    let labels = [
        "K1".to_string(),
        "K2".to_string(),
        "K3".to_string(),
        "O7".to_string()
    ];

    let stripped = labels.iter().map(|l| strip_label(l).to_string() ).collect::<Vec<String>>();

    let plane_of_sym_x = Matrix3::new(
        -1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
    );

    let _ = center_by_centroid(&mut square);
    let operated = square.iter().map(|v| plane_of_sym_x * *v).collect::<Vec<Vector3<f64>>>();

    let (a, b) = best_permutation_multiple_atoms(&square, operated.as_slice(), &stripped);
    let sds = sds_dev(&a, &b);
    assert!((sds - 100.0).abs() < 1e-6, "expected sds = 100, got {sds}");
}

#[test]
fn oc_for_oh_point_group() {

    let mut shape: Vec<Vector3<f64>> = octahedron();
    let labels = [
        "Fe".to_string(),
        "N".to_string(),
        "N".to_string(),
        "N".to_string(),
        "N".to_string(),
        "N".to_string(),
        "N".to_string()
    ];
    let _ = center_by_first_point(&mut shape);
    let _ = normalise(&mut shape);

    let pg_name = "Oh".to_string();

    let average_dev = point_group_dev(&shape, &labels, pg_name);
    assert!(average_dev.is_ok());

    let average_dev = average_dev.unwrap();
    assert!(average_dev.abs() < 1e-3, "expected average_dev = 0.0, found {}",  average_dev);
}

