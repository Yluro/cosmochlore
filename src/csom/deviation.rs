use crate::csom::assignment::best_permutation_multiple_atoms;
use crate::csom::types::CsomError;
use crate::csom::types::CsomOperationResult;
use crate::data::pgs::{get_pointgroup, to_matrix3};
use nalgebra::{Matrix3, Vector3};

/// Shape deviation for two given shapes. It returns squared-distance-sum deviation (0 to 100)
///
/// Assumes centered and normalized points. Scores `reference[i]` against
/// `problem[pairing[i]]`, summing in atom order so the result is bit-stable.
pub(crate) fn sds_dev(
    reference: &[Vector3<f64>],
    problem: &[Vector3<f64>],
    pairing: &[usize],
) -> f64 {
    let n = reference.len() as f64;
    let problem_centroid = problem.iter().sum::<Vector3<f64>>() / n;

    // Permutation-invariant: needs no pairing.
    let denominator: f64 = problem
        .iter()
        .map(|p| (p - problem_centroid).norm_squared())
        .sum();

    let numerator: f64 = reference
        .iter()
        .zip(pairing)
        .map(|(p, &j)| (problem[j] - p).norm_squared())
        .sum();

    100.0 * numerator / denominator
}

/// A symmetry operation by name, its matrix converted once for scoring.
pub(crate) type Operation = (&'static str, Matrix3<f64>);

/// The operations of point group `pg` (the identity is implicit and left out).
pub(crate) fn pointgroup_operations(pg: &str) -> Result<Vec<Operation>, CsomError> {
    let ops =
        get_pointgroup(pg).ok_or_else(|| CsomError::WrongSpaceGroup { pg: pg.to_string() })?;
    Ok(ops.iter().map(|&(name, m)| (name, to_matrix3(m))).collect())
}

/// Scores `points` against one symmetry operation: its image, the pairing and the deviation.
/// Shared by the search, which keeps only `deviation`, and the per-operation report.
fn score_operation(
    points: &[Vector3<f64>],
    groups: &[Vec<usize>],
    ignore_labels: bool,
    has_centre: bool,
    &(name, matrix): &Operation,
) -> CsomOperationResult {
    let image: Vec<Vector3<f64>> = points.iter().map(|p| matrix * p).collect();

    // The assignment only decides which image point each atom is *scored against*; it
    // never changes whose image a point is, so `image` stays in the input's atom order
    // and the pairing is reported separately.
    let pairing =
        best_permutation_multiple_atoms(points, &image, groups, ignore_labels, has_centre);
    CsomOperationResult {
        name,
        matrix,
        deviation: sds_dev(points, &image, &pairing),
        image,
        pairing,
    }
}

/// Deviation of `points` from every individual symmetry operation of `pg`.
pub(crate) fn point_group_operation_deviations(
    points: &[Vector3<f64>],
    groups: &[Vec<usize>],
    pg: &str,
    ignore_labels: bool,
    has_centre: bool,
) -> Result<Vec<CsomOperationResult>, CsomError> {
    Ok(pointgroup_operations(pg)?
        .iter()
        .map(|op| score_operation(points, groups, ignore_labels, has_centre, op))
        .collect())
}

/// Average deviation of `points` over `operations`: the cost of the refinement search.
pub fn point_group_dev(
    points: &[Vector3<f64>],
    groups: &[Vec<usize>],
    operations: &[Operation],
    ignore_labels: bool,
    has_centre: bool,
) -> f64 {
    let sum: f64 = operations
        .iter()
        .map(|op| score_operation(points, groups, ignore_labels, has_centre, op).deviation)
        .sum();

    sum / operations.len() as f64
}
