//! tau4 and tau4' (four ligands).

use super::two_largest_angles;
use crate::geometry::TETRAHEDRAL_ANGLE;
use crate::gidx::{GidxError, GidxResult};
use nalgebra::Vector3;

/// tau4 and tau4' from the two largest angles.
pub(super) fn four_coordinate(vectors: &[Vector3<f64>]) -> Result<GidxResult, GidxError> {
    let (alpha, beta) = two_largest_angles(vectors)?;
    Ok(GidxResult::Four {
        alpha,
        beta,
        tau4: calc_tau4(alpha, beta),
        tau4_prime: calc_tau4_prime(alpha, beta),
    })
}

/// Yang et al. (2007): 1 for a tetrahedron, 0 for a square plane.
pub(in crate::gidx) fn calc_tau4(alpha: f64, beta: f64) -> f64 {
    (360.0 - (alpha + beta)) / (360.0 - 2.0 * TETRAHEDRAL_ANGLE)
}

/// Okuniewski et al. (2015): tau4 that also weighs `beta - alpha`.
pub(in crate::gidx) fn calc_tau4_prime(alpha: f64, beta: f64) -> f64 {
    let theta = TETRAHEDRAL_ANGLE;
    (beta - alpha) / (360.0 - theta) + (180.0 - beta) / (180.0 - theta)
}
