use crate::csom::deviation::point_group_dev;
use crate::csom::prepare::CsomStructure;
use crate::csom::types::{CsomError, OptimiserSettings};
use crate::data::pgs::get_pointgroup_map;
use crate::geometry::rotation_matrix_from_vector;
use argmin::core::{CostFunction, Error, Executor, State};
use argmin::solver::neldermead::NelderMead;
use nalgebra::Vector3;

/// Generates N evenly-spaced points around a sphere of Radius = 1.
fn fibonacci_sphere_sampling(n: usize) -> Vec<Vector3<f64>> {
    if n == 0 {
        return vec![Vector3::new(0.0, 0.0, 1.0)];
    }

    let phi: f64 = std::f64::consts::PI * (5.0f64.sqrt() - 1f64);
    let mut points = Vec::new();

    for i in 0..n {
        let y = 1.0 - (i as f64 / (n as f64 - 1.0));
        let r = (1.0 - y * y).sqrt();
        let theta = phi * i as f64;

        points.push(Vector3::new(r * theta.cos(), y, r * theta.sin()));
    }
    points
}

fn orientation_cost(
    v: &[f64],
    structure: &[Vector3<f64>],
    deviation: &mut impl FnMut(&[Vector3<f64>]) -> f64,
) -> f64 {
    let r = rotation_matrix_from_vector(Vector3::new(v[0], v[1], v[2]));
    let rotated: Vec<Vector3<f64>> = structure.iter().map(|p| r * p).collect();
    deviation(&rotated)
}

/// Binds a structure and a target point group together so argmin can minimize
/// [`orientation_cost`] over rotation vectors: `Param` is `[x, y, z]`, an axis-angle
/// (Rodrigues) rotation vector, and the cost is the average CSM deviation of `structure`
/// from `pg_name` once rotated by that vector.
struct OrientationProblem<'a> {
    structure: &'a CsomStructure,
    pg_name: &'a str,
    ignore_labels: bool,
}

impl CostFunction for OrientationProblem<'_> {
    type Param = Vec<f64>;
    type Output = f64;

    fn cost(&self, v: &Self::Param) -> Result<Self::Output, Error> {
        let mut deviation = |rotated: &[Vector3<f64>]| {
            point_group_dev(
                rotated,
                &self.structure.groups,
                self.pg_name,
                self.ignore_labels,
                self.structure.has_centre,
            )
            .expect(
                "point group name was validated in refine_axis_from_seed before this closure runs",
            )
        };
        Ok(orientation_cost(v, &self.structure.points, &mut deviation))
    }
}

/// Refines a candidate symmetry axis `axis0` with Nelder-Mead, minimizing the average csom
/// deviation of `structure` from the `pg_name` point group.
///
/// Returns the optimised rotation vector (axis-angle, Rodrigues form) and its deviation.
pub(crate) fn refine_axis_from_seed(
    axis0: Vector3<f64>,
    structure: &CsomStructure,
    pg_name: &str,
    max_iters: usize,
    tolerance: f64,
    ignore_labels: bool,
) -> Result<(Vector3<f64>, f64), CsomError> {
    // Fail fast on an unknown point group instead of on every cost evaluation.
    if get_pointgroup_map(pg_name).is_none() {
        return Err(CsomError::WrongSpaceGroup {
            pg: pg_name.to_string(),
        });
    }

    // Initialize a small simplex around the sampled axis (n + 1 = 4 points for 3 parameters).
    let s0 = vec![axis0.x, axis0.y, axis0.z];
    let s1 = vec![axis0.x + 0.1, axis0.y, axis0.z];
    let s2 = vec![axis0.x, axis0.y + 0.1, axis0.z];
    let s3 = vec![axis0.x, axis0.y, axis0.z + 0.1];

    let solver = NelderMead::new(vec![s0, s1, s2, s3])
        .with_sd_tolerance(tolerance)
        .map_err(|e| CsomError::OptimizationFailed(e.to_string()))?;

    let problem = OrientationProblem {
        structure,
        pg_name,
        ignore_labels,
    };

    let result = Executor::new(problem, solver)
        .configure(|state| state.max_iters(max_iters as u64))
        .run()
        .map_err(|e| CsomError::OptimizationFailed(e.to_string()))?;

    let best_v = result
        .state()
        .get_best_param()
        .expect("Nelder-Mead always has a best parameter after running");
    let best_cost = result.state().get_best_cost();

    Ok((Vector3::new(best_v[0], best_v[1], best_v[2]), best_cost))
}

/// Samples `settings.seeds` candidate axes on a Fibonacci sphere and refines each one with
/// [`refine_axis_from_seed`].
///
/// Returns the best rotation vector found and its deviation score.
pub(crate) fn search_best_axis(
    structure: &CsomStructure,
    pg_name: &str,
    settings: OptimiserSettings,
    ignore_labels: bool,
) -> Result<(Vector3<f64>, f64), CsomError> {
    let mut best: Option<(Vector3<f64>, f64)> = None;

    for axis0 in fibonacci_sphere_sampling(settings.seeds) {
        let candidate = refine_axis_from_seed(
            axis0,
            structure,
            pg_name,
            settings.iterations,
            settings.tolerance,
            ignore_labels,
        )?;

        // First sample always wins since best is None.
        // Then only keep the ones that score lower S-value.
        if best.is_none() || candidate.1 < best.as_ref().unwrap().1 {
            best = Some(candidate);
        }
    }

    // fibonacci_sphere_sampling will always yield one point. This error never fires.
    best.ok_or_else(|| CsomError::OptimizationFailed("no candidate axes were sampled".to_string()))
}
