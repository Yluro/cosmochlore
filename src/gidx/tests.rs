use crate::cshm::test_utils::structure_from_shape;
use crate::geometry::{TETRAHEDRAL_ANGLE, angle_btw_vectors};
use crate::gidx::calc::{calc_tau4, calc_tau4_prime, calc_tau6, calculate_gidx, trans_angles};
use crate::gidx::{GidxError, GidxResult};
use crate::xyz::{Atom, Structure, parse_xyz};
use itertools::Itertools;
use nalgebra::{Rotation3, Unit, Vector3};

/// Positions of the built-in shapes in their vertex-count list (`id` in `data/shapes` minus 1).
const SP_4: usize = 0;
const T_4: usize = 1;
const SS_4: usize = 2;
const V_TBPY_4: usize = 3;
const TBPY_5: usize = 2;
const SPY_5: usize = 3;
const HP_6: usize = 0;
const OC_6: usize = 2;
const TPR_6: usize = 3;

/// Values for `tests/ML4.xyz`, `tests/ML5.xyz` and `tests/FeHS.xyz`, worked out from their
/// coordinates by an independent calculation.
const ML4_ALPHA: f64 = 113.23;
const ML4_BETA: f64 = 144.34;
const ML4_TAU4: f64 = 0.7262;
const ML4_TAU4_PRIME: f64 = 0.6298;
const ML5_ALPHA: f64 = 150.18;
const ML5_BETA: f64 = 177.91;
const ML5_TAU5: f64 = 0.4621;
const FEHS_ALPHAS: [f64; 3] = [162.81, 161.60, 161.17];
const FEHS_TAU6: f64 = 0.3023;

/// `(alpha, beta, tau4, tau4')` of a four-coordinate result.
fn four(result: &GidxResult) -> (f64, f64, f64, f64) {
    match *result {
        GidxResult::Four {
            alpha,
            beta,
            tau4,
            tau4_prime,
        } => (alpha, beta, tau4, tau4_prime),
        _ => panic!("expected a four-coordinate result"),
    }
}

/// `(alpha, beta, tau5)` of a five-coordinate result.
fn five(result: &GidxResult) -> (f64, f64, f64) {
    match *result {
        GidxResult::Five { alpha, beta, tau5 } => (alpha, beta, tau5),
        _ => panic!("expected a five-coordinate result"),
    }
}

/// `([alpha1, alpha2, alpha3], tau6)` of a six-coordinate result.
fn six(result: &GidxResult) -> ([f64; 3], f64) {
    match *result {
        GidxResult::Six {
            alpha1,
            alpha2,
            alpha3,
            tau6,
        } => ([alpha1, alpha2, alpha3], tau6),
        _ => panic!("expected a six-coordinate result"),
    }
}

/// A structure whose ligands lie in one plane through the centre, at the given azimuths
/// (degrees) and distances, so every angle between two of them is known by hand.
fn planar_structure(centre: Vector3<f64>, ligands: &[(f64, f64)]) -> Structure {
    Structure {
        centre: Some(Atom {
            label: "M".to_string(),
            coords: centre,
        }),
        ligands: ligands
            .iter()
            .map(|&(azimuth, distance)| Atom {
                label: "L".to_string(),
                coords: centre
                    + distance
                        * Vector3::new(azimuth.to_radians().cos(), azimuth.to_radians().sin(), 0.0),
            })
            .collect(),
    }
}

#[test]
fn square_plane_gives_zero_tau4_and_tau4_prime() {
    let result = calculate_gidx(&structure_from_shape(4, SP_4)).unwrap();
    let (alpha, beta, tau4, tau4_prime) = four(&result);

    assert!((alpha - 180.0).abs() < 1e-6);
    assert!((beta - 180.0).abs() < 1e-6);
    assert!(tau4.abs() < 1e-6);
    assert!(tau4_prime.abs() < 1e-6);
}

#[test]
fn tetrahedron_gives_unit_tau4_and_tau4_prime() {
    // The built-in T-4 coordinates are tabulated to 8 digits, so its angles are
    // arccos(-1/3) = 109.4712...° to about 10⁻⁷ degrees.
    let result = calculate_gidx(&structure_from_shape(4, T_4)).unwrap();
    let (alpha, beta, tau4, tau4_prime) = four(&result);

    assert!((alpha - TETRAHEDRAL_ANGLE).abs() < 1e-6);
    assert!((beta - TETRAHEDRAL_ANGLE).abs() < 1e-6);
    assert!((tau4 - 1.0).abs() < 1e-6);
    assert!((tau4_prime - 1.0).abs() < 1e-6);
}

#[test]
fn exact_tetrahedron_gives_exactly_one_not_the_papers_rounded_value() {
    // Alternate corners of a cube: every angle is arccos(-1/3) to machine precision.
    let corners = [
        Vector3::new(1.0, 1.0, 1.0),
        Vector3::new(1.0, -1.0, -1.0),
        Vector3::new(-1.0, 1.0, -1.0),
        Vector3::new(-1.0, -1.0, 1.0),
    ];
    let structure = Structure {
        centre: Some(Atom {
            label: "M".to_string(),
            coords: Vector3::zeros(),
        }),
        ligands: corners
            .iter()
            .map(|&coords| Atom {
                label: "L".to_string(),
                coords,
            })
            .collect(),
    };
    let result = calculate_gidx(&structure).unwrap();
    let (_, _, tau4, tau4_prime) = four(&result);

    // With the rounded 141° and 109.5° this would be 1.0004.
    assert!((tau4 - 1.0).abs() < 1e-12);
    assert!((tau4_prime - 1.0).abs() < 1e-12);
}

#[test]
fn seesaw_matches_hand_computed_tau4_and_tau4_prime() {
    // Axial-axial 180°, and every other angle is 90°, so alpha = 90° and beta = 180°.
    let result = calculate_gidx(&structure_from_shape(4, SS_4)).unwrap();
    let (alpha, beta, tau4, tau4_prime) = four(&result);

    assert!((alpha - 90.0).abs() < 1e-6);
    assert!((beta - 180.0).abs() < 1e-6);
    assert!((tau4 - 90.0 / (360.0 - 2.0 * TETRAHEDRAL_ANGLE)).abs() < 1e-6);
    assert!((tau4_prime - 90.0 / (360.0 - TETRAHEDRAL_ANGLE)).abs() < 1e-6);
}

#[test]
fn trigonal_pyramid_matches_hand_computed_tau4_and_tau4_prime() {
    // Vacant trigonal bipyramid: axial-equatorial 90°, equatorial-equatorial 120°.
    let result = calculate_gidx(&structure_from_shape(4, V_TBPY_4)).unwrap();
    let (alpha, beta, tau4, tau4_prime) = four(&result);

    assert!((alpha - 120.0).abs() < 1e-6);
    assert!((beta - 120.0).abs() < 1e-6);
    assert!((tau4 - 120.0 / (360.0 - 2.0 * TETRAHEDRAL_ANGLE)).abs() < 1e-6);
    assert!((tau4_prime - 60.0 / (180.0 - TETRAHEDRAL_ANGLE)).abs() < 1e-6);
}

#[test]
fn tau4_prime_separates_angle_pairs_that_share_a_tau4() {
    // (180°, 90°) and (135°, 135°) both add up to 270°, so tau4 cannot tell them apart.
    assert!((calc_tau4(90.0, 180.0) - calc_tau4(135.0, 135.0)).abs() < 1e-12);
    assert!(
        (calc_tau4_prime(90.0, 180.0) - calc_tau4_prime(135.0, 135.0)).abs() > 0.2,
        "tau4' should differ for a seesaw and a flattened tetrahedron"
    );
}

#[test]
fn trigonal_bipyramid_gives_unit_tau5() {
    // Axial-axial 180° and equatorial-equatorial 120°: (180 - 120) / 60.
    let result = calculate_gidx(&structure_from_shape(5, TBPY_5)).unwrap();
    let (alpha, beta, tau5) = five(&result);

    assert!((alpha - 120.0).abs() < 1e-6);
    assert!((beta - 180.0).abs() < 1e-6);
    assert!((tau5 - 1.0).abs() < 1e-6);
}

#[test]
fn square_pyramid_gives_zero_tau5() {
    // The two trans-basal angles are equal (arccos(-0.875) = 151.04°), so beta = alpha.
    let result = calculate_gidx(&structure_from_shape(5, SPY_5)).unwrap();
    let (alpha, beta, tau5) = five(&result);

    assert!((alpha - 151.0442).abs() < 1e-3);
    assert!((beta - 151.0442).abs() < 1e-3);
    assert!(tau5.abs() < 1e-6);
}

#[test]
fn octahedron_gives_zero_tau6() {
    let result = calculate_gidx(&structure_from_shape(6, OC_6)).unwrap();
    let (alphas, tau6) = six(&result);

    for alpha in alphas {
        assert!((alpha - 180.0).abs() < 1e-6);
    }
    assert!(tau6.abs() < 1e-6);
}

#[test]
fn trigonal_prism_matches_the_hand_computed_tau6() {
    // In the uniform prism (equal edges) the vector from the centre to a top vertex and to a
    // diagonally opposite bottom one have cos = -5/7, i.e. 135.58°, not the rounded 135° the
    // paper uses for 0.75: tau6 = (540 - 3 * 135.58) / 180 = 0.7403.
    let result = calculate_gidx(&structure_from_shape(6, TPR_6)).unwrap();
    let (alphas, tau6) = six(&result);
    let trans = (-5.0_f64 / 7.0).acos().to_degrees();

    for alpha in alphas {
        assert!((alpha - trans).abs() < 1e-3);
    }
    assert!((tau6 - (540.0 - 3.0 * trans) / 180.0).abs() < 1e-3);
    assert!((calc_tau6(135.0, 135.0, 135.0) - 0.75).abs() < 1e-12);
}

#[test]
fn planar_hexagon_is_indistinguishable_from_the_octahedron() {
    // Three 180° angles again, and tau6 sees nothing else.
    let result = calculate_gidx(&structure_from_shape(6, HP_6)).unwrap();

    assert!(six(&result).1.abs() < 1e-6);
}

#[test]
fn six_ligands_pair_up_instead_of_taking_the_three_largest_angles() {
    // Azimuths 0, 175, 185, 90, 270 and 45°. The three largest angles overall are 180°
    // (90-270) and two 175° (0-175 and 0-185), but those two share ligand 0. The pairing that
    // uses every ligand once is 90-270 (180°), 0-175 (175°) and 185-45 (140°), whose sum,
    // 495°, beats every other: tau6 = (540 - 495) / 180 = 0.25.
    let ligands = [
        (0.0, 2.0),
        (175.0, 2.1),
        (185.0, 1.9),
        (90.0, 2.0),
        (270.0, 2.2),
        (45.0, 2.0),
    ];

    for order in [[0, 1, 2, 3, 4, 5], [5, 3, 1, 0, 4, 2], [2, 4, 0, 5, 1, 3]] {
        let shuffled: Vec<(f64, f64)> = order.iter().map(|&i| ligands[i]).collect();
        let structure = planar_structure(Vector3::new(0.5, 1.0, -2.0), &shuffled);
        let (alphas, tau6) = six(&calculate_gidx(&structure).unwrap());

        assert!((alphas[0] - 180.0).abs() < 1e-9);
        assert!((alphas[1] - 175.0).abs() < 1e-9);
        assert!((alphas[2] - 140.0).abs() < 1e-9);
        assert!((tau6 - 0.25).abs() < 1e-9);
    }
}

#[test]
fn trans_angles_match_a_brute_force_search_over_every_ordering() {
    // Deterministic pseudo-random vectors; the pairing of (p0,p1), (p2,p3), (p4,p5) is tried
    // for all 720 orderings p of the six ligands.
    let mut state = 0x2545_f491_4f6c_dd1d_u64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
    };

    for _ in 0..50 {
        let vectors: Vec<Vector3<f64>> = (0..6)
            .map(|_| Vector3::new(next(), next(), next()) + Vector3::new(0.1, 0.0, 0.0))
            .collect();

        let brute_force = (0..6)
            .permutations(6)
            .map(|p| {
                angle_btw_vectors(vectors[p[0]], vectors[p[1]])
                    + angle_btw_vectors(vectors[p[2]], vectors[p[3]])
                    + angle_btw_vectors(vectors[p[4]], vectors[p[5]])
            })
            .fold(f64::NEG_INFINITY, f64::max);

        let found: f64 = trans_angles(&vectors).unwrap().iter().sum();
        assert!(
            (found - brute_force).abs() < 1e-9,
            "{found} vs {brute_force}"
        );
    }
}

#[test]
fn four_ligands_use_the_two_largest_angles_wherever_they_come_in_the_file() {
    // Azimuths 0, 100, 180 and 260°: the six angles are 180, 160, 100, 100, 80 and 80, so
    // beta = 180° and alpha = 160°, whichever pair of ligands they come from.
    let ligands = [(0.0, 2.0), (100.0, 2.1), (180.0, 1.9), (260.0, 2.0)];
    let expected_tau4 = 20.0 / (360.0 - 2.0 * TETRAHEDRAL_ANGLE);
    let expected_tau4_prime = 20.0 / (360.0 - TETRAHEDRAL_ANGLE);

    for order in [[0, 1, 2, 3], [3, 1, 0, 2], [2, 3, 1, 0]] {
        let shuffled: Vec<(f64, f64)> = order.iter().map(|&i| ligands[i]).collect();
        let structure = planar_structure(Vector3::new(1.0, -2.0, 3.0), &shuffled);
        let result = calculate_gidx(&structure).unwrap();
        let (alpha, beta, tau4, tau4_prime) = four(&result);

        assert!((alpha - 160.0).abs() < 1e-9);
        assert!((beta - 180.0).abs() < 1e-9);
        assert!((tau4 - expected_tau4).abs() < 1e-9);
        assert!((tau4_prime - expected_tau4_prime).abs() < 1e-9);
    }
}

#[test]
fn five_ligands_use_the_two_largest_angles_wherever_they_come_in_the_file() {
    // Azimuths 0, 50, 130, 190 and 250°: the ten angles are 170, 160, 140, 130, 120, 110,
    // 80, 60, 60 and 50, so beta = 170° and alpha = 160°.
    let ligands = [
        (0.0, 2.0),
        (50.0, 2.2),
        (130.0, 2.0),
        (190.0, 1.8),
        (250.0, 2.1),
    ];

    for order in [[0, 1, 2, 3, 4], [4, 2, 0, 3, 1], [1, 3, 4, 0, 2]] {
        let shuffled: Vec<(f64, f64)> = order.iter().map(|&i| ligands[i]).collect();
        let structure = planar_structure(Vector3::new(-1.5, 0.5, 2.0), &shuffled);
        let result = calculate_gidx(&structure).unwrap();
        let (alpha, beta, tau5) = five(&result);

        assert!((alpha - 160.0).abs() < 1e-9);
        assert!((beta - 170.0).abs() < 1e-9);
        assert!((tau5 - 10.0 / 60.0).abs() < 1e-9);
    }
}

#[test]
fn indices_do_not_change_under_rotation_scaling_or_translation() {
    let rotation =
        Rotation3::from_axis_angle(&Unit::new_normalize(Vector3::new(1.0, 2.0, -0.5)), 0.7);
    let shift = Vector3::new(4.0, -1.0, 9.0);

    for (vertices, index) in [(4, SS_4), (5, TBPY_5), (6, TPR_6)] {
        let reference = calculate_gidx(&structure_from_shape(vertices, index)).unwrap();

        let mut moved = structure_from_shape(vertices, index);
        for atom in moved.ligands.iter_mut().chain(moved.centre.iter_mut()) {
            atom.coords = 2.5 * (rotation * atom.coords) + shift;
        }
        let moved = calculate_gidx(&moved).unwrap();

        match (&reference, &moved) {
            (GidxResult::Four { .. }, GidxResult::Four { .. }) => {
                let (a, b) = (four(&reference), four(&moved));
                assert!((a.0 - b.0).abs() < 1e-9 && (a.1 - b.1).abs() < 1e-9);
                assert!((a.2 - b.2).abs() < 1e-9 && (a.3 - b.3).abs() < 1e-9);
            }
            (GidxResult::Five { .. }, GidxResult::Five { .. }) => {
                let (a, b) = (five(&reference), five(&moved));
                assert!((a.0 - b.0).abs() < 1e-9 && (a.1 - b.1).abs() < 1e-9);
                assert!((a.2 - b.2).abs() < 1e-9);
            }
            _ => {
                let (a, b) = (six(&reference), six(&moved));
                for k in 0..3 {
                    assert!((a.0[k] - b.0[k]).abs() < 1e-9);
                }
                assert!((a.1 - b.1).abs() < 1e-9);
            }
        }
    }
}

#[test]
fn matches_the_hand_checked_four_coordinate_file() {
    // tests/ML4.xyz: a tetrahedron distorted toward a seesaw, with unequal bond lengths.
    let structure = parse_xyz(
        concat!(env!("CARGO_MANIFEST_DIR"), "/tests/ML4.xyz"),
        Some(0),
    )
    .unwrap();
    let result = calculate_gidx(&structure).unwrap();
    let (alpha, beta, tau4, tau4_prime) = four(&result);

    println!("{:?}", result);
    assert!((alpha - ML4_ALPHA).abs() < 1e-2);
    assert!((beta - ML4_BETA).abs() < 1e-2);
    assert!((tau4 - ML4_TAU4).abs() < 1e-4);
    assert!((tau4_prime - ML4_TAU4_PRIME).abs() < 1e-4);
}

#[test]
fn matches_the_hand_checked_five_coordinate_file() {
    // tests/ML5.xyz: a trigonal bipyramid bent about halfway toward a square pyramid.
    let structure = parse_xyz(
        concat!(env!("CARGO_MANIFEST_DIR"), "/tests/ML5.xyz"),
        Some(0),
    )
    .unwrap();
    let result = calculate_gidx(&structure).unwrap();
    let (alpha, beta, tau5) = five(&result);

    println!("{:?}", result);
    assert!((alpha - ML5_ALPHA).abs() < 1e-2);
    assert!((beta - ML5_BETA).abs() < 1e-2);
    assert!((tau5 - ML5_TAU5).abs() < 1e-3);
}

#[test]
fn matches_the_hand_checked_six_coordinate_file() {
    // tests/FeHS.xyz: the high-spin iron(II) complex of the odis example.
    let structure = parse_xyz(
        concat!(env!("CARGO_MANIFEST_DIR"), "/tests/FeHS.xyz"),
        Some(0),
    )
    .unwrap();
    let result = calculate_gidx(&structure).unwrap();
    let (alphas, tau6) = six(&result);

    println!("{:?}", result);
    for (alpha, expected) in alphas.iter().zip(FEHS_ALPHAS) {
        assert!((alpha - expected).abs() < 1e-2);
    }
    assert!((tau6 - FEHS_TAU6).abs() < 1e-4);
}

#[test]
fn reproduces_the_papers_tau6_for_its_six_coordinate_cadmium() {
    // tests/CdONCl4.xyz: Cd2 of Stoeckli-Evans et al., rebuilt from the fractional coordinates
    // of their supporting information. The paper gives trans angles of 170.34, 161.49 and
    // 159.61 degrees and tau6 = [540 - (159.61 + 161.49 + 170.34)] / 180 = 0.27.
    let structure = parse_xyz(
        concat!(env!("CARGO_MANIFEST_DIR"), "/tests/CdONCl4.xyz"),
        Some(0),
    )
    .unwrap();
    let (alphas, tau6) = six(&calculate_gidx(&structure).unwrap());

    for (alpha, published) in alphas.iter().zip([170.34, 161.49, 159.61]) {
        assert!((alpha - published).abs() < 1e-2, "{alpha} vs {published}");
    }
    assert!((tau6 - 0.27).abs() < 5e-3);
    assert!((tau6 - 0.2698).abs() < 1e-4);
}

#[test]
fn reproduces_the_papers_tau5_for_its_five_coordinate_cadmium() {
    // tests/CdONCl3.xyz: Cd1 of the same paper. It gives 161.00 and 142.63 degrees and
    // tau5 = (161.0 - 142.63) / 60 = 0.31.
    let structure = parse_xyz(
        concat!(env!("CARGO_MANIFEST_DIR"), "/tests/CdONCl3.xyz"),
        Some(0),
    )
    .unwrap();
    let (alpha, beta, tau5) = five(&calculate_gidx(&structure).unwrap());

    assert!((alpha - 142.63).abs() < 1e-2);
    assert!((beta - 161.00).abs() < 1e-2);
    assert!((tau5 - 0.31).abs() < 5e-3);
    assert!((tau5 - 0.3063).abs() < 1e-4);
}

#[test]
fn rejects_a_structure_without_a_centre() {
    let mut structure = structure_from_shape(4, T_4);
    structure.centre = None;

    assert!(matches!(
        calculate_gidx(&structure),
        Err(GidxError::NoCentre)
    ));
}

#[test]
fn rejects_coordination_numbers_other_than_four_five_and_six() {
    // Three ligands (the smallest .xyz the parser accepts, with a centre, is 3 atoms).
    let mut three = structure_from_shape(4, T_4);
    three.ligands.truncate(3);
    assert!(matches!(
        calculate_gidx(&three),
        Err(GidxError::UnsupportedCoordination { n: 3 })
    ));

    for vertices in [7, 8] {
        assert!(matches!(
            calculate_gidx(&structure_from_shape(vertices, 1)),
            Err(GidxError::UnsupportedCoordination { n }) if n == vertices as usize
        ));
    }
}

#[test]
fn rejects_a_ligand_on_top_of_the_centre_for_every_coordination_number() {
    // With four or five ligands the NaN angles would otherwise sort to the bottom and the
    // indices would come out of the remaining angles without any sign of a problem.
    for (vertices, index) in [(4, T_4), (5, TBPY_5), (6, OC_6)] {
        let mut structure = structure_from_shape(vertices, index);
        structure.ligands[2].coords = structure.centre.as_ref().unwrap().coords;

        assert!(
            matches!(calculate_gidx(&structure), Err(GidxError::UndefinedAngle)),
            "{vertices} ligands"
        );
    }
}
