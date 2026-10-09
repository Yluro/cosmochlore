use crate::cshm::test_utils::structure_from_shape;
use crate::geometry::{TETRAHEDRAL_ANGLE, angle_btw_vectors};
use crate::gidx::calc::calculate_gidx;
use crate::gidx::calc::tau4::{calc_tau4, calc_tau4_prime};
use crate::gidx::calc::tau6::{calc_tau6, calc_tau6_prime, trans_angles};
use crate::gidx::calc::tau8::{
    SA_DIHEDRAL, TD_PHI, TD_THETA, calc_cube_delta, calc_tau8, calc_tau8_prime,
};
use crate::gidx::{GidxError, GidxResult, SquareFace};
use crate::xyz::{Atom, Structure, parse_xyz};
use itertools::Itertools;
use nalgebra::{Rotation3, Unit, Vector3};

/// Index of each built-in shape in its vertex-count list.
const SP_4: usize = 0;
const T_4: usize = 1;
const SS_4: usize = 2;
const V_TBPY_4: usize = 3;
const TBPY_5: usize = 2;
const SPY_5: usize = 3;
const HP_6: usize = 0;
const PPY_6: usize = 1;
const OC_6: usize = 2;
const TPR_6: usize = 3;
const JPPY_6: usize = 4;
const HBPY_8: usize = 2;
const CU_8: usize = 3;
const SAPR_8: usize = 4;
const TDD_8: usize = 5;
const BTPR_8: usize = 9;

/// Independently computed values for `tests/ML4.xyz`, `ML5.xyz` and `FeHS.xyz`.
const ML4_ALPHA: f64 = 113.23;
const ML4_BETA: f64 = 144.34;
const ML4_TAU4: f64 = 0.7262;
const ML4_TAU4_PRIME: f64 = 0.6298;
const ML5_ALPHA: f64 = 150.18;
const ML5_BETA: f64 = 177.91;
const ML5_TAU5: f64 = 0.4621;
const FEHS_ALPHAS: [f64; 3] = [162.81, 161.60, 161.17];
const FEHS_TAU6: f64 = 0.3023;
const FEHS_TAU6_PRIME: f64 = 0.4330;

/// `tests/tau8/` structures (SI of Turnbull et al.): file and the tau8, tau8' and delta its
/// Table S10 gives.
const TAU8_SI: [(&str, f64, f64, f64); 8] = [
    ("WF4_PH3_4_D2d", 0.92, 1.00, 0.08),
    ("WF4_PH3_4_D2", 0.36, 0.41, 0.05),
    ("WBr4_PH3_4_D2", 0.45, 0.52, 0.07),
    ("MoF4_PH3_4_D2", 0.01, 0.14, 0.13),
    ("MoCl4_PH3_4_D2", 0.55, 0.62, 0.06),
    ("NbCl4_PH3_4_D2", 0.22, 0.32, 0.09),
    ("TaCl4_PH3_4_D2", 0.01, 0.15, 0.14),
    ("TaCl4_PH3_4_D2d", 0.91, 1.00, 0.09),
];

/// `tests/cd6/` sites (CSD CIFs): file, tau6 and tau6', computed independently.
const CD6_SITES: [(&str, f64, f64); 7] = [
    ("IQATAY_Cd1", 0.0000, 0.1201),
    ("IQATAY_Cd2", 0.0586, 0.1879),
    ("IQATAY_Cd3", 0.1589, 0.3053),
    ("IQATAY_Cd4", 0.1733, 0.2934),
    ("ULESAH_Cd1", 0.7336, 0.5110),
    ("VAPSAG_Cd1", 0.8199, 1.0445),
    ("CIJBII_Cd1", 0.6820, 0.9950),
];

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
            ..
        } => ([alpha1, alpha2, alpha3], tau6),
        _ => panic!("expected a six-coordinate result"),
    }
}

/// `([gamma1, ..., gamma5], tau6')` of a six-coordinate result.
fn six_prime(result: &GidxResult) -> ([f64; 5], f64) {
    match *result {
        GidxResult::Six {
            gammas, tau6_prime, ..
        } => (gammas, tau6_prime),
        _ => panic!("expected a six-coordinate result"),
    }
}

/// `(faces, dihedral_a, dihedral_b, tau8')` of an eight-coordinate result.
fn eight(result: &GidxResult) -> ([SquareFace; 2], f64, f64, f64) {
    match *result {
        GidxResult::Eight {
            faces,
            dihedral_a,
            dihedral_b,
            tau8_prime,
            ..
        } => (faces, dihedral_a, dihedral_b, tau8_prime),
        _ => panic!("expected an eight-coordinate result"),
    }
}

/// `cube_delta` of an eight-coordinate result.
fn cube_delta(result: &GidxResult) -> f64 {
    match *result {
        GidxResult::Eight { cube_delta, .. } => cube_delta,
        _ => panic!("expected an eight-coordinate result"),
    }
}

/// A structure with the centre at the origin and the ligands at the given vectors.
fn structure_from_vectors(vectors: &[Vector3<f64>]) -> Structure {
    Structure {
        centre: Some(Atom {
            label: "M".to_string(),
            coords: Vector3::zeros(),
        }),
        ligands: vectors
            .iter()
            .map(|&coords| Atom {
                label: "L".to_string(),
                coords,
            })
            .collect(),
    }
}

/// Kepert's hard-sphere dodecahedron from its defining conditions (see `TD_THETA`).
fn hard_sphere_dodecahedron() -> Vec<Vector3<f64>> {
    let mut w = 0.64_f64;
    for _ in 0..50 {
        w -= (12.0 * w.powi(3) - 7.0 * w * w - 2.0 * w + 1.0) / (36.0 * w * w - 14.0 * w - 2.0);
    }
    let (c_a, c_b) = (w.sqrt(), (2.0 * w - 1.0) / w.sqrt());
    let (s_a, s_b) = ((1.0 - c_a * c_a).sqrt(), (1.0 - c_b * c_b).sqrt());

    vec![
        Vector3::new(-s_a, 0.0, c_a),
        Vector3::new(s_a, 0.0, c_a),
        Vector3::new(0.0, -s_b, c_b),
        Vector3::new(0.0, s_b, c_b),
        Vector3::new(-s_b, 0.0, -c_b),
        Vector3::new(s_b, 0.0, -c_b),
        Vector3::new(0.0, -s_a, -c_a),
        Vector3::new(0.0, s_a, -c_a),
    ]
}

/// The equal-edge square antiprism: tan²α = 2√2, the lower ring turned by 45°.
fn hard_sphere_antiprism() -> Vec<Vector3<f64>> {
    let alpha = (2.0 * 2.0_f64.sqrt()).sqrt().atan();
    let (s, c) = (alpha.sin(), alpha.cos());

    (0..4)
        .map(|k| {
            let azimuth = k as f64 * std::f64::consts::FRAC_PI_2;
            Vector3::new(s * azimuth.cos(), s * azimuth.sin(), c)
        })
        .chain((0..4).map(|k| {
            let azimuth = std::f64::consts::FRAC_PI_4 + k as f64 * std::f64::consts::FRAC_PI_2;
            Vector3::new(s * azimuth.cos(), s * azimuth.sin(), -c)
        }))
        .collect()
}

/// A six-coordinate site of `tests/cd6/`, with the metal as the centre.
fn cd6_structure(name: &str) -> Structure {
    let path = format!("{}/tests/cd6/{name}.xyz", env!("CARGO_MANIFEST_DIR"));
    parse_xyz(&path, Some(0)).unwrap()
}

/// Planar ligands at the given azimuths (degrees) and distances, so every angle is known.
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
fn planar_hexagon_scores_like_the_octahedron_in_tau6_but_not_in_tau6_prime() {
    // Three 180° angles again, and tau6 sees nothing else. Its five smallest angles are 60°
    // against the octahedron's 90°: (450 - 300) / 90 = 5/3.
    let result = calculate_gidx(&structure_from_shape(6, HP_6)).unwrap();

    assert!(six(&result).1.abs() < 1e-6);
    assert!((six_prime(&result).1 - 5.0 / 3.0).abs() < 1e-6);
}

#[test]
fn ideal_shapes_give_the_tau6_prime_anchors() {
    // Five smallest angles: 90° for the octahedron, 72° for the pentagonal pyramids, 60° for
    // the hexagon and arccos(1/7) = 81.79° for the prism with equal edges.
    let prism_angle = (1.0_f64 / 7.0).acos().to_degrees();
    let anchors = [
        (HP_6, 60.0, 5.0 / 3.0),
        (PPY_6, 72.0, 1.0),
        (OC_6, 90.0, 0.0),
        (TPR_6, prism_angle, (450.0 - 5.0 * prism_angle) / 90.0),
        (JPPY_6, 72.0, 1.0),
    ];

    for (index, gamma, expected) in anchors {
        let result = calculate_gidx(&structure_from_shape(6, index)).unwrap();
        let (gammas, tau6_prime) = six_prime(&result);

        for g in gammas {
            assert!((g - gamma).abs() < 1e-3, "shape {index}: {g} vs {gamma}");
        }
        assert!((tau6_prime - expected).abs() < 1e-5, "shape {index}");
    }
    assert!((calc_tau6_prime(&[90.0; 5])).abs() < 1e-12);
    assert!((calc_tau6_prime(&[72.0; 5]) - 1.0).abs() < 1e-12);
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

    for (vertices, index) in [(4, SS_4), (5, TBPY_5), (6, TPR_6), (8, TDD_8)] {
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
            (GidxResult::Eight { .. }, GidxResult::Eight { .. }) => {
                let (a, b) = (eight(&reference), eight(&moved));
                for k in 0..2 {
                    assert!((a.0[k].theta - b.0[k].theta).abs() < 1e-9);
                    assert!((a.0[k].phi - b.0[k].phi).abs() < 1e-9);
                    assert!((a.0[k].tau8 - b.0[k].tau8).abs() < 1e-9);
                }
                assert!((a.1 - b.1).abs() < 1e-9 && (a.2 - b.2).abs() < 1e-9);
                assert!((a.3 - b.3).abs() < 1e-9);
                assert!((cube_delta(&reference) - cube_delta(&moved)).abs() < 1e-9);
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
                let (a, b) = (six_prime(&reference), six_prime(&moved));
                for k in 0..5 {
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
    assert!((six_prime(&result).1 - FEHS_TAU6_PRIME).abs() < 1e-3);
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
    // The 76° bite of the chelate pulls the smallest angles down: 76.15°, 77.99°, 82.93°...
    assert!((six_prime(&calculate_gidx(&structure).unwrap()).1 - 0.4710).abs() < 1e-3);
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
fn literature_octahedra_and_prism_reproduce_the_published_tau6() {
    // Stoeckli-Evans et al. give tau6 = 0, 0.06, 0.16 and 0.17 for the four Cd atoms of IQATAY and
    // 0.73 for the trigonal prism ULESAH.
    let published = [
        ("IQATAY_Cd1", 0.0),
        ("IQATAY_Cd2", 0.06),
        ("IQATAY_Cd3", 0.16),
        ("IQATAY_Cd4", 0.17),
        ("ULESAH_Cd1", 0.73),
    ];

    for (name, value) in published {
        let (_, tau6) = six(&calculate_gidx(&cd6_structure(name)).unwrap());

        assert!((tau6 - value).abs() < 5e-3, "{name}: {tau6} vs {value}");
    }
}

#[test]
fn literature_pentagonal_pyramids_give_tau6_prime_close_to_one() {
    // The paper quotes tau6 = 1.02 for VAPSAG and 1.0 for CIJBII, which it gets from the five
    // adjacent base angles. Those add up to about 360° whatever the distortion, so that value is
    // (1 + tau6') / 2.
    for (name, published) in [("VAPSAG_Cd1", 1.02), ("CIJBII_Cd1", 1.0)] {
        let (gammas, tau6_prime) = six_prime(&calculate_gidx(&cd6_structure(name)).unwrap());

        assert!((tau6_prime - 1.0).abs() < 0.06, "{name}: {tau6_prime}");
        assert!(
            ((1.0 + tau6_prime) / 2.0 - published).abs() < 5e-3,
            "{name}"
        );
        assert!((gammas.iter().sum::<f64>() - 360.0).abs() < 5.0, "{name}");
    }
}

#[test]
fn tau6_prime_orders_the_literature_octahedra_prism_and_pyramids() {
    let tau6_prime = |name: &str| six_prime(&calculate_gidx(&cd6_structure(name)).unwrap()).1;

    let octahedra = ["IQATAY_Cd1", "IQATAY_Cd2", "IQATAY_Cd3", "IQATAY_Cd4"];
    let most_distorted_octahedron = octahedra.iter().map(|n| tau6_prime(n)).fold(0.0, f64::max);

    assert!(most_distorted_octahedron < tau6_prime("ULESAH_Cd1"));
    assert!(tau6_prime("ULESAH_Cd1") < tau6_prime("CIJBII_Cd1"));
    assert!(tau6_prime("ULESAH_Cd1") < tau6_prime("VAPSAG_Cd1"));
}

#[test]
fn literature_sites_match_the_independently_computed_values() {
    for (name, tau6, tau6_prime) in CD6_SITES {
        let result = calculate_gidx(&cd6_structure(name)).unwrap();

        assert!((six(&result).1 - tau6).abs() < 1e-3, "{name} tau6");
        assert!(
            (six_prime(&result).1 - tau6_prime).abs() < 1e-3,
            "{name} tau6'"
        );
    }
}

#[test]
fn hard_sphere_polyhedra_have_the_edges_the_constants_assume() {
    // The dodecahedron: 14 of its 18 edges (the A-A edges and both kinds of A-B edge) have one
    // length; the four B-B edges, between the upper and the lower B sites, are 1.25 times longer.
    let td = hard_sphere_dodecahedron();
    let edge = |i: usize, j: usize| (td[i] - td[j]).norm();
    let a_a = edge(0, 1);
    for (i, j) in [
        (6, 7),
        (0, 2),
        (0, 3),
        (0, 4),
        (1, 2),
        (1, 5),
        (6, 4),
        (7, 5),
    ] {
        assert!((edge(i, j) - a_a).abs() < 1e-12, "edge {i}-{j}");
    }
    assert!((edge(2, 4) / a_a - 1.25).abs() < 0.01);

    // The antiprism: the edges of a ring equal the edges between the rings.
    let sa = hard_sphere_antiprism();
    let ring = (sa[0] - sa[1]).norm();
    assert!(((sa[0] - sa[4]).norm() - ring).abs() < 1e-12);
    assert!(((sa[0] - sa[7]).norm() - ring).abs() < 1e-12);
}

#[test]
fn hard_sphere_dodecahedron_gives_exactly_one_for_both_indices() {
    let result = calculate_gidx(&structure_from_vectors(&hard_sphere_dodecahedron())).unwrap();
    let (faces, dihedral_a, dihedral_b, tau8_prime) = eight(&result);

    for face in faces {
        assert!((face.theta - TD_THETA).abs() < 1e-9, "theta {}", face.theta);
        assert!((face.phi - TD_PHI).abs() < 1e-9, "phi {}", face.phi);
        assert!((face.tau8 - 1.0).abs() < 1e-12);
    }
    assert!((dihedral_a - 90.0).abs() < 1e-9 && (dihedral_b - 90.0).abs() < 1e-9);
    assert!((tau8_prime - 1.0).abs() < 1e-12);
    assert!(cube_delta(&result).abs() < 1e-12);
}

#[test]
fn hard_sphere_antiprism_gives_exactly_zero_for_both_indices() {
    let result = calculate_gidx(&structure_from_vectors(&hard_sphere_antiprism())).unwrap();
    let (faces, dihedral_a, dihedral_b, tau8_prime) = eight(&result);

    for face in faces {
        assert!((face.theta - face.phi).abs() < 1e-9);
        assert!(face.tau8.abs() < 1e-9);
    }
    assert!((dihedral_a - SA_DIHEDRAL).abs() < 1e-9 && (dihedral_b - SA_DIHEDRAL).abs() < 1e-9);
    assert!(tau8_prime.abs() < 1e-9);
    assert!(cube_delta(&result).abs() < 1e-9);
}

#[test]
fn the_constants_are_the_roundings_the_paper_prints() {
    // 129.8, 97.1 and 65.5 (and 2 * 65.5 = 131.0) are what the exact values round to.
    assert!((TD_THETA - 129.8).abs() < 0.05);
    assert!((TD_PHI - 97.1).abs() < 0.05);
    assert!((SA_DIHEDRAL - 65.5).abs() < 0.05);
    assert!((SA_DIHEDRAL - (2.0_f64.sqrt() - 1.0).acos().to_degrees()).abs() < 1e-12);

    // The anchors of both formulas.
    assert!((calc_tau8(TD_THETA, TD_PHI) - 1.0).abs() < 1e-12);
    assert!(calc_tau8(118.5, 118.5).abs() < 1e-12);
    assert!((calc_tau8_prime(90.0, 90.0) - 1.0).abs() < 1e-12);
    assert!(calc_tau8_prime(SA_DIHEDRAL, SA_DIHEDRAL).abs() < 1e-12);
}

#[test]
fn built_in_antiprism_and_dodecahedron_give_the_ideal_indices() {
    let (faces, dihedral_a, dihedral_b, tau8_prime) =
        eight(&calculate_gidx(&structure_from_shape(8, SAPR_8)).unwrap());
    for face in faces {
        assert!(face.tau8.abs() < 1e-6);
    }
    assert!((dihedral_a - SA_DIHEDRAL).abs() < 1e-5 && (dihedral_b - SA_DIHEDRAL).abs() < 1e-5);
    assert!(tau8_prime.abs() < 1e-6);

    // SHAPE tabulates TDD-8 to 8 digits, so it hits 1 to about 10⁻⁴.
    let (faces, dihedral_a, dihedral_b, tau8_prime) =
        eight(&calculate_gidx(&structure_from_shape(8, TDD_8)).unwrap());
    for face in faces {
        assert!((face.tau8 - 1.0).abs() < 1e-3, "tau8 {}", face.tau8);
    }
    assert!((dihedral_a - 90.0).abs() < 1e-4 && (dihedral_b - 90.0).abs() < 1e-4);
    assert!((tau8_prime - 1.0).abs() < 1e-4);
}

#[test]
fn cube_has_equal_diagonals_but_orthogonal_planes() {
    // Both faces: zero tau8 like the antiprism, tau8' = 1 like the dodecahedron.
    let result = calculate_gidx(&structure_from_shape(8, CU_8)).unwrap();
    let (faces, dihedral_a, dihedral_b, tau8_prime) = eight(&result);

    for face in faces {
        assert!(face.tau8.abs() < 1e-6);
        // Opposite corners of a cube face subtend arccos(-1/3) at its centre.
        assert!((face.theta - TETRAHEDRAL_ANGLE).abs() < 1e-4);
    }
    assert!((dihedral_a - 90.0).abs() < 1e-4 && (dihedral_b - 90.0).abs() < 1e-4);
    assert!((tau8_prime - 1.0).abs() < 1e-4);
    // The cube is the one polyhedron with delta = 1.
    assert!((cube_delta(&result) - 1.0).abs() < 1e-4);
}

#[test]
fn bicapped_trigonal_prism_has_contradictory_faces() {
    // The signature of a BTP: one square face looks like a dodecahedron, the other like an
    // antiprism or cube. The SI quotes tau8 = 1.10 and 0.00 for WF6(py)2.
    let (faces, ..) = eight(&calculate_gidx(&structure_from_shape(8, BTPR_8)).unwrap());
    let (low, high) = if faces[0].tau8 < faces[1].tau8 {
        (faces[0].tau8, faces[1].tau8)
    } else {
        (faces[1].tau8, faces[0].tau8)
    };

    assert!(low.abs() < 1e-3, "{low}");
    assert!((high - 1.17).abs() < 0.01, "{high}");
}

#[test]
fn reproduces_the_supporting_information_table_for_its_computed_structures() {
    // Table S10 rounds to two decimals, and NbCl4(PH3)4 (D2) differs by 0.01 in tau8'.
    for (name, tau8, tau8_prime, delta) in TAU8_SI {
        let path = format!("{}/tests/tau8/{name}.xyz", env!("CARGO_MANIFEST_DIR"));
        let structure = parse_xyz(&path, Some(0)).unwrap();
        let result = calculate_gidx(&structure).unwrap();
        let (faces, _, _, computed_prime) = eight(&result);

        for face in faces {
            assert!(
                (face.tau8 - tau8).abs() < 0.015,
                "{name}: tau8 {}",
                face.tau8
            );
        }
        assert!(
            (computed_prime - tau8_prime).abs() < 0.015,
            "{name}: tau8' {computed_prime}"
        );
        assert!((cube_delta(&result) - delta).abs() < 0.015, "{name}: delta");
    }
}

#[test]
fn cube_delta_takes_the_larger_tau8_as_the_supporting_information_does() {
    // Its WF6(py)2: tau8 = 1.10 and 0.00, tau8' = 0.61, delta = -0.49 (in either face order).
    assert!((calc_cube_delta(1.10, 0.0, 0.61) + 0.49).abs() < 1e-12);
    assert!((calc_cube_delta(0.0, 1.10, 0.61) + 0.49).abs() < 1e-12);
    // And [UF7(NH3)]3-: tau8 = 1.11 and 0.47, tau8' = 0.76, delta = -0.36 (unrounded inputs).
    assert!((calc_cube_delta(1.11, 0.47, 0.76) + 0.36).abs() < 0.015);

    // The ideal bicapped trigonal prism: 0.408 - 1.167.
    let result = calculate_gidx(&structure_from_shape(8, BTPR_8)).unwrap();
    assert!((cube_delta(&result) + 0.759).abs() < 0.005);
}

#[test]
fn hexagonal_bipyramid_has_no_plane_through_its_axial_ligands() {
    // Its two axial ligands are antipodal, so the plane they share with the centre is undefined.
    assert!(matches!(
        calculate_gidx(&structure_from_shape(8, HBPY_8)),
        Err(GidxError::UndefinedPlane)
    ));
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
fn rejects_coordination_numbers_other_than_four_five_six_and_eight() {
    // Three ligands (the smallest .xyz the parser accepts, with a centre, is 3 atoms).
    let mut three = structure_from_shape(4, T_4);
    three.ligands.truncate(3);
    assert!(matches!(
        calculate_gidx(&three),
        Err(GidxError::UnsupportedCoordination { n: 3 })
    ));

    for vertices in [7, 9] {
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
    for (vertices, index) in [(4, T_4), (5, TBPY_5), (6, OC_6), (8, TDD_8)] {
        let mut structure = structure_from_shape(vertices, index);
        structure.ligands[2].coords = structure.centre.as_ref().unwrap().coords;

        assert!(
            matches!(calculate_gidx(&structure), Err(GidxError::UndefinedAngle)),
            "{vertices} ligands"
        );
    }
}
