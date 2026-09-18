//! Starting orientations for the symmetry-axis search.
//!
//! The point-group tables put the principal axis along `z`, so the search has to find the
//! rotation that carries the structure's own principal axis there. Its seeds are built in two
//! levels:
//!
//! 1. Candidate axis directions are spread over a Fibonacci hemisphere, and each one is
//!    carried onto `z` by the smallest rotation that does so. This guarantees that, whatever
//!    the input orientation, some seed starts with the structure's true axis within the
//!    lattice spacing of `z`. (Using the lattice directions directly as rotation vectors -- a
//!    rotation by exactly 1 rad about each -- does not: it can only bring an axis to within
//!    |1 rad - tilt| of `z`, so an axis tilted far from `z` starts every refinement far from
//!    the answer, which is how the search used to miss the C6 axis of a tilted ring.)
//! 2. Aligning the axis fixes only two of the three rotational degrees of freedom. When the
//!    group also cares about the spin about `z` -- when it has vertical mirror planes or
//!    perpendicular C2 axes -- each direction is additionally spun about `z` through one
//!    [`in_plane_period`] in steps no coarser than the lattice spacing, so the seeds cover the
//!    in-plane angle at the same resolution as the axis tilt.
//!
//! Every seed is then refined by `optimize::refine_axis_from_seed`.

use crate::data::pgs::{SymmetryOperation, to_matrix3};
use nalgebra::{Matrix3, Rotation3, Vector3};
use std::f64::consts::{PI, TAU};

/// Spreads `n` unit vectors evenly over the upper hemisphere (z >= 0) with a Fibonacci
/// lattice. Each one is a candidate direction for the structure's principal symmetry axis.
///
/// Returns `n` candidate z-axes as Vec<Vector3<f64>>
pub(crate) fn fibonacci_hemisphere(n: usize) -> Vec<Vector3<f64>> {
    if n == 0 {
        return vec![Vector3::z()];
    }

    let golden_angle: f64 = PI * (3.0 - 5.0f64.sqrt());

    (0..n)
        .map(|i| {
            // z runs from 1 (the pole) down to 1/n, stopping short of the equator, where a
            // sample and its antipode would seed the same axis twice.
            let z = 1.0 - i as f64 / n as f64;
            let r = (1.0 - z * z).sqrt();
            let theta = golden_angle * i as f64;

            Vector3::new(r * theta.cos(), r * theta.sin(), z)
        })
        .collect()
}

/// In-plane period of a point group: the smallest spin about `z` that carries the group's
/// operation set onto itself (by conjugation), so that two orientations of a structure
/// differing by that spin score identically.
///
/// `None` when the group is invariant under *every* spin about `z`, which is the case when all
/// of its operations are rotations or rotoreflections about `z` (C_n, S_n, C_nh, Cs, Ci, E):
/// the in-plane angle then does not affect the deviation at all.
pub(crate) fn in_plane_period(ops: &[SymmetryOperation]) -> Option<f64> {
    const EPS: f64 = 1e-9;

    let matrices: Vec<Matrix3<f64>> = ops.iter().map(|(_, m)| to_matrix3(*m)).collect();

    // An operation commutes with every spin about z iff it is a rotation or rotoreflection
    // about z itself: [[c, -s, 0], [s, c, 0], [0, 0, +-1]].
    let is_axial = |m: &Matrix3<f64>| {
        m[(0, 2)].abs() < EPS
            && m[(1, 2)].abs() < EPS
            && m[(2, 0)].abs() < EPS
            && m[(2, 1)].abs() < EPS
            && (m[(0, 0)] - m[(1, 1)]).abs() < EPS
            && (m[(0, 1)] + m[(1, 0)]).abs() < EPS
    };
    if matrices.iter().all(is_axial) {
        return None;
    }

    // Order of the principal axis: the finest proper rotation about z in the table.
    let n = matrices
        .iter()
        .filter(|m| is_axial(m) && (m[(2, 2)] - 1.0).abs() < EPS)
        .map(|m| m[(1, 0)].atan2(m[(0, 0)]).abs())
        .filter(|&angle| angle > EPS)
        .fold(1usize, |n, angle| n.max((TAU / angle).round() as usize));

    // Some tables (I, Ih, D7h, D8h) list their off-axis matrices to three decimals only, so
    // membership is tested well above that precision.
    let maps_group_onto_itself = |angle: f64| {
        let spin = *Rotation3::from_axis_angle(&Vector3::z_axis(), angle).matrix();
        let spin_inv = spin.transpose();
        matrices.iter().all(|m| {
            let conjugated = spin * m * spin_inv;
            matrices
                .iter()
                .any(|other| (conjugated - other).norm() < 1e-2)
        })
    };

    // The period is 2*pi/q for the largest q that works: at least n, since spinning by C_n
    // itself changes nothing, and at most 2n, the spacing of the vertical planes and C2 axes
    // of the dihedral groups (it is n, not 2n, for O and Oh, whose in-plane axes alternate
    // between C4 and C2 and so only repeat every 90 degrees).
    let q = (1..=2 * n)
        .rev()
        .find(|&q| maps_group_onto_itself(TAU / q as f64))
        .expect("q = 1 is a full turn, which always maps the group onto itself");

    Some(TAU / q as f64)
}

/// The starting rotation vectors (axis-angle, Rodrigues form) of a search for the point group
/// with operations `ops`: every [`fibonacci_hemisphere`] direction of `n_directions`, carried
/// onto `z`, at every in-plane angle the group tells apart (see the module docs).
pub(crate) fn seed_rotation_vectors(
    n_directions: usize,
    ops: &[SymmetryOperation],
) -> Vec<Vector3<f64>> {
    // Angular spacing of the direction lattice: each of the n samples owns a 2*pi/n patch
    // of the hemisphere, taken as a square of that area.
    let spacing = (TAU / n_directions.max(1) as f64).sqrt();

    let (n_in_plane, step) = match in_plane_period(ops) {
        Some(period) => {
            let m = (period / spacing).ceil().max(1.0);
            (m as usize, period / m)
        }
        None => (1, 0.0),
    };

    let z = Vector3::z();
    fibonacci_hemisphere(n_directions)
        .into_iter()
        .flat_map(|d| {
            // `d` is never antiparallel to z on the upper hemisphere, so this never fails;
            // for `d = z` it is the identity.
            let align = Rotation3::rotation_between(&d, &z).unwrap_or_else(Rotation3::identity);
            (0..n_in_plane).map(move |k| {
                let spin = Rotation3::from_axis_angle(&Vector3::z_axis(), k as f64 * step);
                (spin * align).scaled_axis()
            })
        })
        .collect()
}
