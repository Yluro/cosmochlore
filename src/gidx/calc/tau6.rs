//! tau6 and tau6' (six ligands).

use super::{angle_matrix, ensure_defined};
use crate::geometry::pairwise_angles;
use crate::gidx::{GidxError, GidxResult};
use nalgebra::Vector3;

/// tau6 from the three trans angles, tau6' from the five smallest.
pub(super) fn six_coordinate(vectors: &[Vector3<f64>]) -> Result<GidxResult, GidxError> {
    let [alpha1, alpha2, alpha3] = trans_angles(vectors)?;
    let gammas = smallest_angles(vectors)?;
    Ok(GidxResult::Six {
        alpha1,
        alpha2,
        alpha3,
        tau6: calc_tau6(alpha1, alpha2, alpha3),
        gammas,
        tau6_prime: calc_tau6_prime(&gammas),
    })
}

/// The five smallest angles between two vectors, smallest first; they may share ligands.
fn smallest_angles(vectors: &[Vector3<f64>]) -> Result<[f64; 5], GidxError> {
    let mut angles = pairwise_angles(vectors);
    ensure_defined(angles.iter().copied())?;
    angles.sort_by(|a, b| a.total_cmp(b));

    Ok([angles[0], angles[1], angles[2], angles[3], angles[4]])
}

/// The three trans angles of six vectors, largest first: the disjoint pairing with the
/// greatest angle sum (the three largest angles can share a ligand).
pub(in crate::gidx) fn trans_angles(vectors: &[Vector3<f64>]) -> Result<[f64; 3], GidxError> {
    let between = angle_matrix(vectors)?;

    let free: Vec<usize> = (0..vectors.len()).collect();
    let (_, mut pairs) = best_pairing(&between, &free);
    pairs.sort_by(|a, b| b.total_cmp(a));

    Ok([pairs[0], pairs[1], pairs[2]])
}

/// Pairs up `free` with the greatest angle sum; returns the sum and the pair angles.
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

/// Stoeckli-Evans et al. (2025): 0 for an octahedron, from the trans angles alone.
pub(in crate::gidx) fn calc_tau6(alpha1: f64, alpha2: f64, alpha3: f64) -> f64 {
    (3.0 * 180.0 - (alpha1 + alpha2 + alpha3)) / 180.0
}

/// cosmochlore's own complement to tau6: 0 for an octahedron (cis angles 90°), 1 for a
/// pentagonal pyramid (base angles 72°). The paper's pyramid value is (1 + tau6') / 2.
pub(in crate::gidx) fn calc_tau6_prime(gammas: &[f64; 5]) -> f64 {
    (5.0 * 90.0 - gammas.iter().sum::<f64>()) / 90.0
}
