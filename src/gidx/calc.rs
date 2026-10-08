use crate::geometry::{TETRAHEDRAL_ANGLE, pairwise_angles};
use crate::gidx::{GidxError, GidxResult};
use crate::xyz::{Atom, Structure};
use itertools::Itertools;
use nalgebra::Vector3;

/// Angle between the equatorial and axial planes of the trigonal bipyramid in Addison's tau5.
const TAU5_DENOMINATOR: f64 = 60.0;

/// Calculate the geometry indices (tau4 and tau4' for four ligands, tau5 for five, tau6 for six)
/// of a structure.
///
/// The structure must contain a centre atom and exactly four, five or six ligand points. Only
/// the ligand-centre-ligand angles enter, so bond lengths play no part.
pub fn calculate_gidx(structure: &Structure) -> Result<GidxResult, GidxError> {
    let centre = structure.centre.as_ref().ok_or(GidxError::NoCentre)?;
    let vectors = centre_to_ligand_vectors(centre, &structure.ligands);

    match vectors.len() {
        4 => four_coordinate(&vectors),
        5 => five_coordinate(&vectors),
        6 => six_coordinate(&vectors),
        n => Err(GidxError::UnsupportedCoordination { n }),
    }
}

/// The vector from the centre to each ligand, in the order the ligands come in the file.
fn centre_to_ligand_vectors(centre: &Atom, ligands: &[Atom]) -> Vec<Vector3<f64>> {
    ligands
        .iter()
        .map(|ligand| ligand.coords - centre.coords)
        .collect()
}

/// tau4 and tau4' from the two largest angles between the four ligand vectors.
fn four_coordinate(vectors: &[Vector3<f64>]) -> Result<GidxResult, GidxError> {
    let (alpha, beta) = two_largest_angles(vectors)?;
    Ok(GidxResult::Four {
        alpha,
        beta,
        tau4: calc_tau4(alpha, beta),
        tau4_prime: calc_tau4_prime(alpha, beta),
    })
}

/// tau5 from the two largest angles between the five ligand vectors.
fn five_coordinate(vectors: &[Vector3<f64>]) -> Result<GidxResult, GidxError> {
    let (alpha, beta) = two_largest_angles(vectors)?;
    Ok(GidxResult::Five {
        alpha,
        beta,
        tau5: calc_tau5(alpha, beta),
    })
}

/// tau6 from the three trans angles of the six ligand vectors.
fn six_coordinate(vectors: &[Vector3<f64>]) -> Result<GidxResult, GidxError> {
    let [alpha1, alpha2, alpha3] = trans_angles(vectors)?;
    Ok(GidxResult::Six {
        alpha1,
        alpha2,
        alpha3,
        tau6: calc_tau6(alpha1, alpha2, alpha3),
    })
}

/// Fails if any of the angles is NaN. A zero-length vector (a ligand on top of the centre) gives
/// one, and a later sort could silently push it aside and leave a plausible-looking index. The
/// parser already rejects non-finite coordinates.
fn ensure_defined(angles: impl IntoIterator<Item = f64>) -> Result<(), GidxError> {
    if angles.into_iter().any(f64::is_nan) {
        return Err(GidxError::UndefinedAngle);
    }
    Ok(())
}

/// The second-largest (`alpha`) and largest (`beta`) of all the angles between pairs of vectors.
fn two_largest_angles(vectors: &[Vector3<f64>]) -> Result<(f64, f64), GidxError> {
    let mut angles = pairwise_angles(vectors);
    ensure_defined(angles.iter().copied())?;
    angles.sort_by(|a, b| b.total_cmp(a));

    Ok((angles[1], angles[0]))
}

/// The three principal (trans) angles of six vectors, largest first: the angles of the three
/// disjoint pairs of vectors whose sum is greatest, so each ligand is used exactly once.
///
/// The three largest angles overall are not always a valid choice: in a distorted structure two
/// of them can share a ligand.
pub(super) fn trans_angles(vectors: &[Vector3<f64>]) -> Result<[f64; 3], GidxError> {
    let n = vectors.len();
    let angles = pairwise_angles(vectors);
    ensure_defined(angles.iter().copied())?;

    // `pairwise_angles` lists the angles in the order of these index pairs.
    let mut between = vec![vec![0.0; n]; n];
    for ([i, j], angle) in (0..n).array_combinations().zip(angles) {
        between[i][j] = angle;
        between[j][i] = angle;
    }

    let free: Vec<usize> = (0..n).collect();
    let (_, mut pairs) = best_pairing(&between, &free);
    pairs.sort_by(|a, b| b.total_cmp(a));

    Ok([pairs[0], pairs[1], pairs[2]])
}

/// Pairs up every index in `free` so that the sum of `between[i][j]` over the pairs is greatest.
/// Returns that sum and the angle of each pair. `free` must hold an even number of indices.
fn best_pairing(between: &[Vec<f64>], free: &[usize]) -> (f64, Vec<f64>) {
    let Some((&first, others)) = free.split_first() else {
        return (0.0, Vec::new());
    };

    let mut best: Option<(f64, Vec<f64>)> = None;
    for (k, &partner) in others.iter().enumerate() {
        let rest: Vec<usize> = others
            .iter()
            .enumerate()
            .filter(|&(m, _)| m != k)
            .map(|(_, &i)| i)
            .collect();

        let angle = between[first][partner];
        let (rest_sum, mut pairs) = best_pairing(between, &rest);
        let sum = angle + rest_sum;

        if best.as_ref().is_none_or(|(best_sum, _)| sum > *best_sum) {
            pairs.push(angle);
            best = Some((sum, pairs));
        }
    }

    best.expect("an even, non-empty set always has a partner for its first index")
}

/// Yang, Powell and Houser's four-coordinate index (Dalton Trans. 2007, 955): 1 for a regular
/// tetrahedron and 0 for a square plane. It falls to about 0.64 for the seesaw.
pub(super) fn calc_tau4(alpha: f64, beta: f64) -> f64 {
    (360.0 - (alpha + beta)) / (360.0 - 2.0 * TETRAHEDRAL_ANGLE)
}

/// Okuniewski et al.'s refinement of tau4 (Polyhedron 2015, 90, 47): same limits, but it also
/// weighs the difference between the two angles. Geometries whose `alpha + beta` add up the same
/// but whose angles differ, like a seesaw (180°, 90°) and a flattened tetrahedron (135°, 135°),
/// share a tau4 yet no longer share a tau4'.
pub(super) fn calc_tau4_prime(alpha: f64, beta: f64) -> f64 {
    let theta = TETRAHEDRAL_ANGLE;
    (beta - alpha) / (360.0 - theta) + (180.0 - beta) / (180.0 - theta)
}

/// Addison et al.'s five-coordinate index (J. Chem. Soc., Dalton Trans. 1984, 1349): 0 for a
/// square pyramid and 1 for a trigonal bipyramid.
fn calc_tau5(alpha: f64, beta: f64) -> f64 {
    (beta - alpha) / TAU5_DENOMINATOR
}

/// Stoeckli-Evans et al.'s six-coordinate index (Acta Cryst. E 2025, 81, 393): 0 for an
/// octahedron, where the three trans angles are 180°, and growing as they close up. It is built
/// from the trans angles alone, so a planar hexagon (also three 180° angles) gives 0 as well.
pub(super) fn calc_tau6(alpha1: f64, alpha2: f64, alpha3: f64) -> f64 {
    (3.0 * 180.0 - (alpha1 + alpha2 + alpha3)) / 180.0
}
