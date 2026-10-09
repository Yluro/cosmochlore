//! tau5 (five ligands).

use super::two_largest_angles;
use crate::gidx::{GidxError, GidxResult};
use nalgebra::Vector3;

/// Denominator of Addison's tau5, in degrees.
const TAU5_DENOMINATOR: f64 = 60.0;

/// tau5 from the two largest angles.
pub(super) fn five_coordinate(vectors: &[Vector3<f64>]) -> Result<GidxResult, GidxError> {
    let (alpha, beta) = two_largest_angles(vectors)?;
    Ok(GidxResult::Five {
        alpha,
        beta,
        tau5: calc_tau5(alpha, beta),
    })
}

/// Addison et al. (1984): 0 for a square pyramid, 1 for a trigonal bipyramid.
fn calc_tau5(alpha: f64, beta: f64) -> f64 {
    (beta - alpha) / TAU5_DENOMINATOR
}
