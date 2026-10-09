//! tau8, tau8' and delta (eight ligands).

use super::angle_matrix;
use crate::gidx::{GidxError, GidxResult, SquareFace};
use itertools::Itertools;
use nalgebra::Vector3;

/// Face-diagonal angles in degrees of Kepert's hard-sphere trigonal dodecahedron: A-A
/// (`TD_THETA`) and B-B (`TD_PHI`). The paper rounds them to 129.8° and 97.1°.
///
/// All edges but the four B-B ones are equal, so cos²θA = w, the root of
/// 12w³ − 7w² − 2w + 1 = 0 near 0.64, and tan θA · tan θB = 2. The diagonals subtend
/// arccos(−cos²θA) and arccos(−cos²θB).
pub(in crate::gidx) const TD_THETA: f64 = 129.820_772_928_068_46;

pub(in crate::gidx) const TD_PHI: f64 = 97.071_371_178_684_94;

/// Angle in degrees between the A (or B) planes of the equal-edge square antiprism,
/// arccos(√2 − 1). The paper rounds it to 65.5°; the dodecahedron has 90°.
pub(in crate::gidx) const SA_DIHEDRAL: f64 = 65.530_199_479_297_8;

/// Relative cross-product length below which two vectors count as collinear.
const COLLINEAR: f64 = 1e-9;

/// tau8 (per square face), tau8' and delta, following the SI of Turnbull et al. (2021).
///
/// Face 1 holds the first ligand and face 2 the other four. The larger diagonal of a face
/// joins its A ligands, the smaller its B ligands.
pub(super) fn eight_coordinate(vectors: &[Vector3<f64>]) -> Result<GidxResult, GidxError> {
    let between = angle_matrix(vectors)?;

    let first = square_face_diagonals(vectors, &between, 0);
    let remaining: Vec<usize> = (0..8)
        .filter(|i| !first.iter().flatten().contains(i))
        .collect();
    let second = complement_diagonals(
        &between,
        [remaining[0], remaining[1], remaining[2], remaining[3]],
    );

    let (face_1, a_1, b_1) = square_face(&between, first);
    let (face_2, a_2, b_2) = square_face(&between, second);

    let delta_a = paired_dihedral(vectors, &between, [a_1[0], a_1[1], a_2[0], a_2[1]], false)?;
    let delta_b = paired_dihedral(vectors, &between, [b_1[0], b_1[1], b_2[0], b_2[1]], true)?;

    let tau8_prime = calc_tau8_prime(delta_a, delta_b);
    Ok(GidxResult::Eight {
        faces: [face_1, face_2],
        delta_a,
        delta_b,
        tau8_prime,
        cube_delta: calc_cube_delta(face_1.tau8, face_2.tau8, tau8_prime),
    })
}

/// Diagonals of the square face holding `first`: its partner is at the third largest of its
/// seven angles, and the other corners are the two ligands closest to both.
fn square_face_diagonals(
    vectors: &[Vector3<f64>],
    between: &[Vec<f64>],
    first: usize,
) -> [[usize; 2]; 2] {
    let mut partners: Vec<usize> = (0..vectors.len()).filter(|&i| i != first).collect();
    partners.sort_by(|&a, &b| between[first][b].total_cmp(&between[first][a]));
    let second = partners[2];

    let distance_to_both =
        |i: usize| (vectors[i] - vectors[first]).norm() + (vectors[i] - vectors[second]).norm();
    let mut corners: Vec<usize> = (0..vectors.len())
        .filter(|&i| i != first && i != second)
        .collect();
    corners.sort_by(|&a, &b| distance_to_both(a).total_cmp(&distance_to_both(b)));

    [[first, second], [corners[0], corners[1]]]
}

/// Splits four ligands into two diagonals: the pairing with the largest angle sum.
fn complement_diagonals(between: &[Vec<f64>], four: [usize; 4]) -> [[usize; 2]; 2] {
    let [a, b, c, d] = four;
    let sum = |pairs: &[[usize; 2]; 2]| {
        between[pairs[0][0]][pairs[0][1]] + between[pairs[1][0]][pairs[1][1]]
    };

    [[[a, b], [c, d]], [[a, c], [b, d]], [[a, d], [b, c]]]
        .into_iter()
        .max_by(|x, y| sum(x).total_cmp(&sum(y)))
        .expect("there are always three pairings")
}

/// A face from its two diagonals, with its A pair (larger angle) and B pair.
fn square_face(
    between: &[Vec<f64>],
    diagonals: [[usize; 2]; 2],
) -> (SquareFace, [usize; 2], [usize; 2]) {
    let [d1, d2] = diagonals;
    let (angle_1, angle_2) = (between[d1[0]][d1[1]], between[d2[0]][d2[1]]);
    let (a, b, theta, phi) = if angle_1 >= angle_2 {
        (d1, d2, angle_1, angle_2)
    } else {
        (d2, d1, angle_2, angle_1)
    };

    let face = SquareFace {
        theta,
        phi,
        tau8: calc_tau8(theta, phi),
    };
    (face, a, b)
}

/// Acute angle between the planes through the centre and two pairs of four sites: the
/// first is paired with its smallest-angle partner (A) or largest (B, `largest`).
fn paired_dihedral(
    vectors: &[Vector3<f64>],
    between: &[Vec<f64>],
    sites: [usize; 4],
    largest: bool,
) -> Result<f64, GidxError> {
    let first = sites[0];
    let by_angle = |a: &usize, b: &usize| between[first][*a].total_cmp(&between[first][*b]);
    let candidates = sites[1..].iter().copied();
    let partner = if largest {
        candidates.max_by(by_angle)
    } else {
        candidates.min_by(by_angle)
    }
    .expect("three candidates");
    let others: Vec<usize> = sites[1..]
        .iter()
        .copied()
        .filter(|&i| i != partner)
        .collect();

    plane_angle(
        vectors[first],
        vectors[partner],
        vectors[others[0]],
        vectors[others[1]],
    )
}

/// Acute angle between the planes through the centre and (`a`, `b`) and (`c`, `d`).
fn plane_angle(
    a: Vector3<f64>,
    b: Vector3<f64>,
    c: Vector3<f64>,
    d: Vector3<f64>,
) -> Result<f64, GidxError> {
    let (n1, n2) = (a.cross(&b), c.cross(&d));
    if n1.norm() <= COLLINEAR * a.norm() * b.norm() || n2.norm() <= COLLINEAR * c.norm() * d.norm()
    {
        return Err(GidxError::UndefinedPlane);
    }

    Ok((n1.dot(&n2).abs() / (n1.norm() * n2.norm()))
        .clamp(0.0, 1.0)
        .acos()
        .to_degrees())
}

/// Turnbull et al. (2021): 0 for equal diagonals (antiprism, cube), 1 for the dodecahedron.
pub(in crate::gidx) fn calc_tau8(theta: f64, phi: f64) -> f64 {
    (theta - phi).abs() / (TD_THETA - TD_PHI)
}

/// Turnbull et al. (2021): 1 for the dodecahedron and the cube, 0 for the antiprism.
pub(in crate::gidx) fn calc_tau8_prime(delta_a: f64, delta_b: f64) -> f64 {
    1.0 - (180.0 - (delta_a + delta_b)).abs() / (180.0 - 2.0 * SA_DIHEDRAL)
}

/// Δ = tau8' − tau8 with the larger tau8, as in the SI: 1 for a cube, 0 for the ideal
/// dodecahedron and antiprism.
pub(in crate::gidx) fn calc_cube_delta(tau8_1: f64, tau8_2: f64, tau8_prime: f64) -> f64 {
    tau8_prime - tau8_1.max(tau8_2)
}
