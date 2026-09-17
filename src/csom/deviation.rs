use crate::csom::assignment::best_permutation_multiple_atoms;
use crate::csom::types::CsomError;
use crate::csom::types::CsomOperation;
use crate::data::pgs::{get_pointgroup, to_matrix3, SymmetryOperation};
use nalgebra::Vector3;

/// Shape deviation for two given shapes. It returns squared-distance-sum deviation (0 to 100)
///
/// Assumes centered and normalized points.
/// Assumes correct point-to-point correspondence.
pub(crate) fn sds_dev(reference: &[Vector3<f64>], problem: &[Vector3<f64>]) -> f64 {
    let n = reference.len() as f64;
    let problem_centroid = problem.iter().sum::<Vector3<f64>>() / n;

    let denominator: f64 = problem.iter()
        .map(|p| (p - problem_centroid).norm_squared())
        .sum();

    let numerator: f64 = reference.iter().zip(problem)
        .map(
            |(p , q)| (q - p).norm_squared()
        ).sum();

    100.0 * numerator / denominator
}

/// Deviation of `points` from every individual symmetry operation of `pg`.
pub(crate) fn point_group_operation_deviations(
    points: &[Vector3<f64>],
    groups: &[Vec<usize>],
    pg: &str,
    ignore_labels: bool,
    has_centre: bool,
) -> Result<Vec<CsomOperation>, CsomError> {
    let ops = get_pointgroup(pg).ok_or_else(|| CsomError::WrongSpaceGroup { pg: pg.to_string() })?;

    Ok(ops
        .iter()
        .map(|op| operation_deviation(points, groups, ignore_labels, has_centre, op))
        .collect())
}

/// Deviation of `points` from a single symmetry operation.
fn operation_deviation(
    points: &[Vector3<f64>],
    groups: &[Vec<usize>],
    ignore_labels: bool,
    has_centre: bool,
    (name, matrix): &SymmetryOperation,
) -> CsomOperation {
    let sym_op = to_matrix3(*matrix);
    let image: Vec<Vector3<f64>> = points.iter().map(|p| sym_op * p).collect();

    // The assignment only decides which image point each atom is *scored against*; it
    // never changes whose image a point is, so `image` stays in the input's atom order
    // and the pairing is reported separately.
    let (a, b, pairing) = best_permutation_multiple_atoms(points, &image, groups, ignore_labels, has_centre);

    CsomOperation {
        name: name.to_string(),
        matrix: sym_op,
        deviation: sds_dev(&a, &b),
        image,
        pairing,
    }
}

/// Average CSM deviation of `points` from every symmetry operation of point group `pg`.
pub fn point_group_dev(
    points: &[Vector3<f64>],
    groups: &[Vec<usize>],
    pg: &str,
    ignore_labels: bool,
    has_centre: bool,
) -> Result<f64, CsomError> {
    let devs = point_group_operation_deviations(points, groups, pg, ignore_labels, has_centre)?;

    Ok(devs.iter().map(|op| op.deviation).sum::<f64>() / devs.len() as f64)
}
