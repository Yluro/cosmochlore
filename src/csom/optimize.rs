use nalgebra::Vector3;
use crate::geometry::rotation_matrix_from_vector;
use argmin::core::{CostFunction, Error, Executor, State};
use argmin::solver::neldermead::NelderMead;
use crate::csom::io::CsomStructure;
use crate::csom::dev::point_group_dev;
use crate::csom::CsomError;
use crate::data::pgs::get_pointgroup_map;

/// Generates N evenly-spaced points around a sphere of Radius = 1.
fn fibonacci_sphere_sampling(n: usize) -> Vec<Vector3<f64>> {
    if n == 0 { return vec![Vector3::new(0.0, 0.0, 1.0)]; }

    let phi: f64 = std::f64::consts::PI * (5.0f64.sqrt() - 1f64);
    let mut points = Vec::new();

    for i in 0..n {

        let y = 1.0 - (i as f64 / (n as f64 - 1.0));
        let r = (1.0 - y*y).sqrt();
        let theta = phi * i as f64;

        points.push(Vector3::new(r * theta.cos(), y, r * theta.sin()));

    };
    points
}
pub fn orientation_cost(
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
}

impl CostFunction for OrientationProblem<'_> {
    type Param = Vec<f64>;
    type Output = f64;

    fn cost(&self, v: &Self::Param) -> Result<Self::Output, Error> {
        let mut deviation = |rotated: &[Vector3<f64>]| {
            point_group_dev(rotated, &self.structure.labels, self.pg_name.to_string())
                .expect("point group name was validated before optimisation started")
        };
        Ok(orientation_cost(v, &self.structure.points, &mut deviation))
    }
}

/// Refines a candidate symmetry axis `axis0` with Nelder-Mead, minimizing the average CSM
/// deviation of `structure` from the `pg_name` point group.
///
/// Returns the optimised rotation vector (axis-angle, Rodrigues form) and its deviation.
fn optimise_axis(
    axis0: Vector3<f64>,
    structure: &CsomStructure,
    pg_name: &str,
) -> Result<(Vector3<f64>, f64), CsomError> {
    // Fail fast on an unknown point group instead of on every cost evaluation.
    if get_pointgroup_map(pg_name).is_none() {
        return Err(CsomError::WrongSpaceGroup { pg: pg_name.to_string() });
    }

    // Initialize a small simplex around the sampled axis (n + 1 = 4 points for 3 parameters).
    let s0 = vec![axis0.x, axis0.y, axis0.z];
    let s1 = vec![axis0.x + 0.001, axis0.y, axis0.z];
    let s2 = vec![axis0.x + 0.001, axis0.y + 0.001, axis0.z];
    let s3 = vec![axis0.x + 0.001, axis0.y, axis0.z + 0.001];

    let solver = NelderMead::new(vec![s0, s1, s2, s3]);

    let problem = OrientationProblem { structure, pg_name };

    let result = Executor::new(problem, solver)
        .configure(|state| state.max_iters(1000))
        .run()
        .map_err(|e| CsomError::OptimizationFailed(e.to_string()))?;

    let best_v = result
        .state()
        .get_best_param()
        .expect("Nelder-Mead always has a best parameter after running");
    let best_cost = result.state().get_best_cost();

    Ok((Vector3::new(best_v[0], best_v[1], best_v[2]), best_cost))
}


#[cfg(test)]
mod tests {
    use super::*;

    fn octahedron() -> CsomStructure {
        let points = vec![
            Vector3::new(0.0, 0.0, -1.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(-1.0, 0.0, 0.0),
            Vector3::new(0.0, -1.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ];
        let labels = vec!["N".to_string(); 6];
        CsomStructure { labels, points }
    }

    #[test]
    fn optimise_axis_converges_for_perfect_octahedron() {
        let structure = octahedron();
        // Start from a slightly off-identity guess so the simplex has real work to do.
        let axis0 = Vector3::new(0.01, 0.02, 0.03);

        let (_, cost) = optimise_axis(axis0, &structure, "Oh").expect("Oh is a valid point group");

        assert!(cost.abs() < 1e-3, "expected near-zero deviation, got {cost}");
    }

    #[test]
    fn optimise_axis_rejects_unknown_point_group() {
        let structure = octahedron();
        let axis0 = Vector3::zeros();

        let result = optimise_axis(axis0, &structure, "NotAGroup");

        assert!(matches!(result, Err(CsomError::WrongSpaceGroup { .. })));
    }
}
