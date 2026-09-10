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

    let (a, b, _) = best_permutation_multiple_atoms(&shape, rotated.as_slice(), &stripped, false, false);

    // Sanity: same length, all atoms retained.
    assert_eq!(a.len(), shape.len());
    assert_eq!(b.len(), shape.len());

    let sds = sds_dev(&a, &b);
    assert!(sds > 1e-6, "expected nonzero deviation due to mismatched Cl group, got {sds}");
}

#[test]
fn ignoring_labels_finds_the_perfect_match_across_mismatched_labels() {

    let mut shape: Vec<Vector3<f64>> = octahedron();
    normalise(&mut shape);

    // Deliberately mislabeled compared to `mismatched_gives_nonzero_dev`: with labels
    // honoured, the lone "Cl" ligand can't be matched into the "N" group, so a perfect
    // rotation still shows nonzero deviation. Ignoring labels lifts that restriction.
    let labels = [
        "Fe".to_string(),
        "N".to_string(),
        "Cl".to_string(),
        "N".to_string(),
        "N".to_string(),
        "N".to_string(),
        "N".to_string(),
    ];

    let axis = Vector3::new(0.0, 0.0, 1.0);
    let angle = 90.0; // Exact C4 rotation: a true symmetry operation of the octahedron.

    let rot_mat = rotation_matrix(axis, angle);
    let rotated: Vec<Vector3<f64>> = shape.iter().map(|v| rot_mat * *v).collect::<Vec<Vector3<f64>>>();

    let (_, honouring_labels, _) = best_permutation_multiple_atoms(&shape, rotated.as_slice(), &labels, false, false);
    assert!(sds_dev(&shape, &honouring_labels) > 1e-6, "mismatched label should prevent a perfect match");

    let (a, b, _) = best_permutation_multiple_atoms(&shape, rotated.as_slice(), &labels, true, false);
    let sds = sds_dev(&a, &b);
    assert!(sds < 1e-6, "expected near-zero deviation once labels are ignored, got {sds}");
}

#[test]
fn pinned_centre_atom_is_never_reassigned() {

    let mut shape: Vec<Vector3<f64>> = octahedron(); // index 0 sits at the origin
    normalise(&mut shape);

    // All ligands share one label, so without pinning the Hungarian search would be free
    // (and, for a symmetry operation, equally correct) to route the origin point anywhere.
    let labels = vec!["N".to_string(); shape.len()];

    let axis = Vector3::new(0.0, 0.0, 1.0);
    let angle = 90.0;
    let rot_mat = rotation_matrix(axis, angle);
    let rotated: Vec<Vector3<f64>> = shape.iter().map(|v| rot_mat * *v).collect::<Vec<Vector3<f64>>>();

    let (a, b, _) = best_permutation_multiple_atoms(&shape, rotated.as_slice(), &labels, false, true);

    // The pinned centre (origin) must appear paired with itself among the returned pairs.
    let centre_pair = a.iter().zip(b.iter()).find(|(pa, _)| pa.norm() < 1e-9);
    let (_, centre_b) = centre_pair.expect("origin point should be present in the output");
    assert!(centre_b.norm() < 1e-9, "pinned centre should be matched to itself, got {centre_b}");
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

    let (a, b, _) = best_permutation_multiple_atoms(&square, operated.as_slice(), &stripped, false, false);
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

    let pg_name = "Oh";

    let average_dev = point_group_dev(&shape, &labels, pg_name, false, true);
    assert!(average_dev.is_ok());

    let average_dev = average_dev.unwrap();
    assert!(average_dev.abs() < 1e-3, "expected average_dev = 0.0, found {}",  average_dev);
}

#[test]
fn operated_image_keeps_every_atom_with_its_own_label() {

    let mut shape: Vec<Vector3<f64>> = octahedron();
    let _ = center_by_first_point(&mut shape);
    let _ = normalise(&mut shape);

    // One ligand deliberately mislabelled, so ignoring labels lets it pair with the image of
    // a nitrogen -- the case that used to make the reconstructed .xyz misleading, because the
    // image was filed under the label of the atom it was *matched to* rather than its own.
    let labels: Vec<String> = ["Fe", "N", "Cl", "N", "N", "N", "N"]
        .iter()
        .map(|l| l.to_string())
        .collect();

    let ops = point_group_operation_deviations(&shape, &labels, "Oh", true, true)
        .expect("Oh is a valid point group");

    for op in &ops {
        for (i, point) in op.image.iter().enumerate() {
            let expected = op.matrix * shape[i];
            assert!(
                (point - expected).norm() < 1e-12,
                "{}: image[{i}] should be the operation applied to atom {i}, got {point}",
                op.name
            );
        }

        // The pairing has to stay a permutation of the atoms, one image per atom.
        let mut seen = vec![false; labels.len()];
        for &j in &op.pairing {
            assert!(!seen[j], "{}: atom {j} was paired with twice", op.name);
            seen[j] = true;
        }
    }

    // The mislabelled ligand must really be paired across elements for some operation,
    // otherwise the check above never sees the case that broke.
    let crossed = ops.iter().any(|op| {
        op.pairing.iter().enumerate().any(|(i, &j)| labels[i] != labels[j])
    });
    assert!(crossed, "expected ignoring labels to pair the Cl with a nitrogen image somewhere");
}
