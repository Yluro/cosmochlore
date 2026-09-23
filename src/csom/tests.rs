use crate::csom::assignment::*;
use crate::csom::deviation::*;
use crate::csom::optimize::{refine_axis_from_seed, search_best_axis};
use crate::csom::prepare::{CsomStructure, strip_label};
use crate::csom::seeding::{fibonacci_hemisphere, in_plane_period, seed_rotation_vectors};
use crate::csom::types::{CsomError, OptimiserSettings};
use crate::data::pgs::get_pointgroup;
use crate::geometry::{
    center_by_centroid, center_by_first_point, normalise, rotation_matrix,
    rotation_matrix_from_vector,
};
use nalgebra::{Matrix3, Rotation3, Vector3};
use std::f64::consts::{PI, TAU};

fn octahedron() -> Vec<Vector3<f64>> {
    Vec::from([
        // Regular octahedron centered at origo
        Vector3::zeros(),             // Point 0
        Vector3::new(0.0, 0.0, -1.0), // 1
        Vector3::new(1.0, 0.0, 0.0),  // 2
        Vector3::new(0.0, 1.0, 0.0),  // 3
        Vector3::new(-1.0, 0.0, 0.0), // 4
        Vector3::new(0.0, -1.0, 0.0), // 5
        Vector3::new(0.0, 0.0, 1.0),  //6
    ])
}

#[test]
fn can_rotate_and_assign() {
    let mut shape: Vec<Vector3<f64>> = octahedron();
    normalise(&mut shape);

    let axis = Vector3::new(0.0, 0.0, 1.0);
    let angle = 85.0; // Approximate C4 rotation.

    let rot_mat = rotation_matrix(axis, angle);

    let rotated: Vec<Vector3<f64>> = shape
        .iter()
        .map(|v| rot_mat * *v)
        .collect::<Vec<Vector3<f64>>>();

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

    let stripped = labels
        .iter()
        .map(|l| strip_label(l).to_string())
        .collect::<Vec<String>>();
    assert_eq!(vec!["Fe", "N", "Cl", "Cl", "Cl", "O", "N"], stripped);

    let axis = Vector3::new(0.0, 0.0, 1.0);
    let angle = 90.0; // Approximate C4 rotation.

    let rot_mat = rotation_matrix(axis, angle);

    let rotated: Vec<Vector3<f64>> = shape
        .iter()
        .map(|v| rot_mat * *v)
        .collect::<Vec<Vector3<f64>>>();

    let groups = group_by_label(&stripped);
    let (a, b, _) =
        best_permutation_multiple_atoms(&shape, rotated.as_slice(), &groups, false, false);

    // Sanity: same length, all atoms retained.
    assert_eq!(a.len(), shape.len());
    assert_eq!(b.len(), shape.len());

    let sds = sds_dev(&a, &b);
    assert!(
        sds > 1e-6,
        "expected nonzero deviation due to mismatched Cl group, got {sds}"
    );
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
    let rotated: Vec<Vector3<f64>> = shape
        .iter()
        .map(|v| rot_mat * *v)
        .collect::<Vec<Vector3<f64>>>();

    let groups = group_by_label(&labels);
    let (_, honouring_labels, _) =
        best_permutation_multiple_atoms(&shape, rotated.as_slice(), &groups, false, false);
    assert!(
        sds_dev(&shape, &honouring_labels) > 1e-6,
        "mismatched label should prevent a perfect match"
    );

    let (a, b, _) =
        best_permutation_multiple_atoms(&shape, rotated.as_slice(), &groups, true, false);
    let sds = sds_dev(&a, &b);
    assert!(
        sds < 1e-6,
        "expected near-zero deviation once labels are ignored, got {sds}"
    );
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
    let rotated: Vec<Vector3<f64>> = shape
        .iter()
        .map(|v| rot_mat * *v)
        .collect::<Vec<Vector3<f64>>>();

    let groups = group_by_label(&labels);
    let (a, b, _) =
        best_permutation_multiple_atoms(&shape, rotated.as_slice(), &groups, false, true);

    // The pinned centre (origin) must appear paired with itself among the returned pairs.
    let centre_pair = a.iter().zip(b.iter()).find(|(pa, _)| pa.norm() < 1e-9);
    let (_, centre_b) = centre_pair.expect("origin point should be present in the output");
    assert!(
        centre_b.norm() < 1e-9,
        "pinned centre should be matched to itself, got {centre_b}"
    );
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
        "O7".to_string(),
    ];

    let stripped = labels
        .iter()
        .map(|l| strip_label(l).to_string())
        .collect::<Vec<String>>();

    let plane_of_sym_x = Matrix3::new(-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

    let _ = center_by_centroid(&mut square);
    let operated = square
        .iter()
        .map(|v| plane_of_sym_x * *v)
        .collect::<Vec<Vector3<f64>>>();

    let groups = group_by_label(&stripped);
    let (a, b, _) =
        best_permutation_multiple_atoms(&square, operated.as_slice(), &groups, false, false);
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
        "N".to_string(),
    ];
    let _ = center_by_first_point(&mut shape);
    let _ = normalise(&mut shape);

    let pg_name = "Oh";

    let groups = group_by_label(&labels);
    let average_dev = point_group_dev(&shape, &groups, pg_name, false, true);
    assert!(average_dev.is_ok());

    let average_dev = average_dev.unwrap();
    assert!(
        average_dev.abs() < 1e-3,
        "expected average_dev = 0.0, found {}",
        average_dev
    );
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

    let groups = group_by_label(&labels);
    let ops = point_group_operation_deviations(&shape, &groups, "Oh", true, true)
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
        op.pairing
            .iter()
            .enumerate()
            .any(|(i, &j)| labels[i] != labels[j])
    });
    assert!(
        crossed,
        "expected ignoring labels to pair the Cl with a nitrogen image somewhere"
    );
}

/// A tighter search than the CLI default, so the optimiser tests converge well inside their
/// tolerances: `seeds` seeds, 1000 iterations, tolerance 1e-8.
fn thorough_search(seeds: usize) -> OptimiserSettings {
    OptimiserSettings {
        seeds,
        iterations: 1000,
        tolerance: 1e-8,
    }
}

fn octahedron_structure() -> CsomStructure {
    let points = vec![
        Vector3::new(0.0, 0.0, -1.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(-1.0, 0.0, 0.0),
        Vector3::new(0.0, -1.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
    ];
    let labels = vec!["N".to_string(); 6];
    let groups = group_by_label(&labels);
    CsomStructure {
        points,
        has_centre: false,
        groups,
    }
}

/// A water molecule (C2v) with O–H = 0.9584 A, H–O–H = 104.45°
fn water_structure() -> CsomStructure {
    let bond = 0.9584;
    let half_angle = 104.45f64.to_radians() / 2.0;

    let mut points = vec![
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(bond * half_angle.sin(), 0.0, -bond * half_angle.cos()),
        Vector3::new(-bond * half_angle.sin(), 0.0, -bond * half_angle.cos()),
    ];
    center_by_centroid(&mut points);

    let labels = vec!["O".to_string(), "H".to_string(), "H".to_string()];
    let groups = group_by_label(&labels);
    CsomStructure {
        points,
        has_centre: false,
        groups,
    }
}

#[test]
fn refine_axis_from_seed_converges_for_perfect_octahedron() {
    let structure = octahedron_structure();
    // Start from a slightly off-identity guess so the simplex has real work to do.
    let axis0 = Vector3::new(0.01, 0.02, 0.03);

    let (_, cost) = refine_axis_from_seed(axis0, &structure, "Oh", 1000, 1e-8, false)
        .expect("Oh is a valid point group");

    assert!(
        cost.abs() < 1e-3,
        "expected near-zero deviation, got {cost}"
    );
}

#[test]
fn refine_axis_from_seed_rejects_unknown_point_group() {
    let structure = octahedron_structure();
    let axis0 = Vector3::zeros();

    let result = refine_axis_from_seed(axis0, &structure, "NotAGroup", 1000, 1e-8, false);

    assert!(matches!(result, Err(CsomError::WrongSpaceGroup { .. })));
}

#[test]
fn search_best_axis_converges_for_rotated_octahedron() {
    let mut structure = octahedron_structure();

    // Rotate the octahedron off its canonical axes so the search has to find the
    // symmetry axis rather than starting right on top of it.
    let rot_mat = rotation_matrix_from_vector(Vector3::new(1.0, 1.0, 1.0));
    structure.points = structure.points.iter().map(|p| rot_mat * p).collect();

    let (_, cost) = search_best_axis(&structure, "Oh", thorough_search(8), false)
        .expect("Oh is a valid point group");

    assert!(
        cost.abs() < 1e-3,
        "expected near-zero deviation, got {cost}"
    );
}

#[test]
fn search_best_axis_recovers_the_c2_axis_of_a_rotated_water_molecule() {
    let mut structure = water_structure();

    // An arbitrary rotation with no special relation to water's own C2 axis.
    let applied_rotation = rotation_matrix_from_vector(Vector3::new(0.4, -0.3, 0.9));
    structure.points = structure
        .points
        .iter()
        .map(|p| applied_rotation * p)
        .collect();

    let (axis, cost) = search_best_axis(&structure, "C2v", thorough_search(20), false)
        .expect("C2v is a valid point group");
    assert!(
        cost.abs() < 1e-3,
        "expected near-zero deviation, got {cost}"
    );

    // The rotation the optimiser found should undo `applied_rotation` well enough
    // that composing the two carries water's own C2 axis (z, in the untouched
    // molecule) back onto itself.
    let recovered_rotation = rotation_matrix_from_vector(axis);
    let recovered_z_axis = recovered_rotation * applied_rotation * Vector3::z();

    assert!(
        (recovered_z_axis - Vector3::z()).norm() < 1e-2,
        "expected the recovered rotation to realign water's C2 axis with z, got {recovered_z_axis}"
    );
}

#[test]
fn search_best_axis_rejects_unknown_point_group() {
    let structure = octahedron_structure();

    let result = search_best_axis(&structure, "NotAGroup", thorough_search(8), false);

    assert!(matches!(result, Err(CsomError::WrongSpaceGroup { .. })));
}

/// The pyridine ring of `tests/Co110_pyridine.xyz` (six atoms, no centre), as `csom --nc`
/// prepares it. Its normal is tilted 68 degrees from `z`.
fn pyridine_ring_structure() -> CsomStructure {
    let mut points = vec![
        Vector3::new(2.178776, 6.987902, 12.271005),
        Vector3::new(0.680230, 6.204588, 14.464450),
        Vector3::new(1.483930, 5.816450, 12.233931),
        Vector3::new(0.735126, 5.411294, 13.325930),
        Vector3::new(1.384024, 7.394160, 14.494005),
        Vector3::new(2.124076, 7.759082, 13.373102),
    ];
    center_by_centroid(&mut points);
    normalise(&mut points);

    let labels: Vec<String> = ["N", "C", "C", "C", "C", "C"]
        .iter()
        .map(|l| l.to_string())
        .collect();
    let groups = group_by_label(&labels);
    CsomStructure {
        points,
        has_centre: false,
        groups,
    }
}

#[test]
fn fibonacci_hemisphere_samples_unit_vectors_on_the_upper_hemisphere() {
    let directions = fibonacci_hemisphere(20);

    assert_eq!(directions.len(), 20);
    assert!(
        (directions[0] - Vector3::z()).norm() < 1e-12,
        "the pole should be the first sample, got {}",
        directions[0]
    );
    for d in &directions {
        assert!((d.norm() - 1.0).abs() < 1e-12, "not a unit vector: {d}");
        assert!(d.z > 0.0, "sample below the equator: {d}");
    }
}

#[test]
fn in_plane_period_is_the_spacing_of_the_vertical_elements() {
    let period = |pg: &str| in_plane_period(get_pointgroup(pg).expect("valid point group"));
    let expect = |pg: &str, q: f64| {
        let found = period(pg).unwrap_or_else(|| panic!("{pg} should have an in-plane period"));
        assert!(
            (found - TAU / q).abs() < 1e-9,
            "{pg}: expected a period of 2pi/{q}, got 2pi/{}",
            TAU / found
        );
    };

    // Dihedral-type groups: vertical planes / perpendicular C2 axes every pi/n.
    expect("C2v", 4.0);
    expect("C3v", 6.0);
    expect("D2h", 4.0);
    expect("D3h", 6.0);
    expect("D6h", 12.0);
    expect("D7h", 14.0);
    expect("D8h", 16.0);
    expect("D2d", 4.0);
    expect("D4d", 8.0);
    // Cubic groups: z is a C2 axis of T/Td and a C4 axis of O/Oh, but both only repeat every
    // 90 degrees (the in-plane axes of O alternate between C4 and C2).
    expect("Td", 4.0);
    expect("Th", 4.0);
    expect("O", 4.0);
    expect("Oh", 4.0);
    // Icosahedral: z is a C5 axis. Its five equatorial C2 axes repeat every 36 degrees, but
    // the other five C5 axes sit at azimuths 0, 72, 144, ... and only repeat every 72.
    expect("I", 5.0);
    expect("Ih", 5.0);

    // Every operation is about z: the in-plane angle is irrelevant.
    for pg in ["E", "Cs", "Ci", "C2", "C6", "C7", "S4", "C3h"] {
        assert!(
            period(pg).is_none(),
            "{pg} should be invariant under any spin about z"
        );
    }
}

#[test]
fn every_seed_aligns_its_direction_with_z() {
    let n = 20;
    let seeds = seed_rotation_vectors(n, get_pointgroup("C2v").expect("valid point group"));

    // C2v's period of 90 degrees needs three in-plane spins at 20 directions (lattice spacing
    // ~32 degrees), so every direction contributes three seeds.
    let per_direction = seeds.len() / n;
    assert_eq!(seeds.len(), n * per_direction);
    assert_eq!(per_direction, 3);

    for (d, chunk) in fibonacci_hemisphere(n)
        .iter()
        .zip(seeds.chunks(per_direction))
    {
        let first = rotation_matrix_from_vector(chunk[0]);
        for (k, seed) in chunk.iter().enumerate() {
            let rotation = rotation_matrix_from_vector(*seed);
            let aligned = rotation * d;
            assert!(
                (aligned - Vector3::z()).norm() < 1e-9,
                "seed for {d} should carry it onto z, got {aligned}"
            );

            // The in-plane spins of one direction differ from the first only by a rotation
            // about z of k * period / 3.
            let spin = Rotation3::from_matrix_unchecked(rotation * first.transpose());
            let expected = Rotation3::from_axis_angle(&Vector3::z_axis(), k as f64 * PI / 6.0);
            assert!(
                (spin.matrix() - expected.matrix()).norm() < 1e-9,
                "in-plane spin {k} of {d} should be {} degrees about z",
                30 * k
            );
        }
    }

    // With no in-plane dependence (C6 is invariant under any spin about z) there is exactly
    // one seed per direction, and the first seed is the identity: the input orientation is
    // always tried.
    let seeds = seed_rotation_vectors(n, get_pointgroup("C6").expect("valid point group"));
    assert_eq!(seeds.len(), n);
    assert!(seeds[0].norm() < 1e-12, "first seed should be the identity");
}

#[test]
fn default_search_finds_the_c6_axis_of_a_tilted_ring() {
    // Regression: seeding with rotations by a fixed 1 rad about lattice axes could not bring
    // this ring's normal (68 degrees from z) anywhere near z, and every refinement ended in
    // the in-plane local minimum at a deviation of 34 instead of 0.03.
    let structure = pyridine_ring_structure();

    let (_, cost) = search_best_axis(&structure, "D6h", OptimiserSettings::default(), true)
        .expect("D6h is a valid point group");

    assert!(
        cost < 0.1,
        "expected the near-perfect C6 axis of the ring to be found, got a deviation of {cost}"
    );
}
