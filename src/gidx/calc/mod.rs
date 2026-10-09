//! Geometry indices by coordination number, with the angle helpers they share.

pub(super) mod tau4;
pub(super) mod tau5;
pub(super) mod tau6;
pub(super) mod tau8;

use crate::geometry::pairwise_angles;
use crate::gidx::{GidxError, GidxResult};
use crate::xyz::{Atom, Structure};
use itertools::Itertools;
use nalgebra::Vector3;

/// Geometry indices of a centre with four, five, six or eight ligands, from the
/// ligand-centre-ligand angles alone.
pub fn calculate_gidx(structure: &Structure) -> Result<GidxResult, GidxError> {
    let centre = structure.centre.as_ref().ok_or(GidxError::NoCentre)?;
    let vectors = centre_to_ligand_vectors(centre, &structure.ligands);

    match vectors.len() {
        4 => tau4::four_coordinate(&vectors),
        5 => tau5::five_coordinate(&vectors),
        6 => tau6::six_coordinate(&vectors),
        8 => tau8::eight_coordinate(&vectors),
        n => Err(GidxError::UnsupportedCoordination { n }),
    }
}

/// Vectors from the centre to each ligand, in file order.
fn centre_to_ligand_vectors(centre: &Atom, ligands: &[Atom]) -> Vec<Vector3<f64>> {
    ligands
        .iter()
        .map(|ligand| ligand.coords - centre.coords)
        .collect()
}

/// Fails on a NaN angle (e.g. a ligand on the centre), which a sort could hide.
fn ensure_defined(angles: impl IntoIterator<Item = f64>) -> Result<(), GidxError> {
    if angles.into_iter().any(f64::is_nan) {
        return Err(GidxError::UndefinedAngle);
    }
    Ok(())
}

/// Symmetric matrix of the angles between every two vectors.
fn angle_matrix(vectors: &[Vector3<f64>]) -> Result<Vec<Vec<f64>>, GidxError> {
    let n = vectors.len();
    let angles = pairwise_angles(vectors);
    ensure_defined(angles.iter().copied())?;

    // `pairwise_angles` lists the angles in the order of these index pairs.
    let mut between = vec![vec![0.0; n]; n];
    for ([i, j], angle) in (0..n).array_combinations().zip(angles) {
        between[i][j] = angle;
        between[j][i] = angle;
    }

    Ok(between)
}

/// The second-largest (`alpha`) and largest (`beta`) angle between two vectors.
fn two_largest_angles(vectors: &[Vector3<f64>]) -> Result<(f64, f64), GidxError> {
    let mut angles = pairwise_angles(vectors);
    ensure_defined(angles.iter().copied())?;
    angles.sort_by(|a, b| b.total_cmp(a));

    Ok((angles[1], angles[0]))
}
