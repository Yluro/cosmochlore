use crate::cli::CsomArgs;
use crate::xyz::{parse_xyz, Structure};
use std::error::Error;
use nalgebra::{Matrix3, Vector3};
use crate::csom::dev::point_group_operation_deviations;
use crate::csom::io::{prepare_csom_structure, CenteringMode};
use crate::csom::optimize::find_best_axis;
use crate::data::pgs::POINTGROUP_NAMES;
use crate::geometry::rotation_matrix_from_vector;
use crate::out::{print_csom_table, write_csom_csv, write_csom_details_csv};

pub(crate) mod io;
mod dev;
#[cfg(test)]
mod tests;
mod optimize;

pub struct CsomResult {
    /// Point group analysed.
    pub point_group: String,

    /// Deviation from the ideal symmetry
    pub deviation: f64,

    /// Rotation matrix that defines the refined axis.
    pub rotation: Matrix3<f64>,

    /// Deviation from each individual symmetry operation.
    pub operations: Vec<(String, Matrix3<f64>, f64)>,
}


pub fn csom_main(args: CsomArgs) -> Result<(), Box<dyn Error>> {

    // 1. Parse input .xyz file and form structure.
    let structure = parse_xyz(&args.name, args.not_centered)?;

    // 2. Fetch the desired point groups.
    let point_groups = args.point_groups;

    if point_groups.is_none() {
        todo!("Auto point-group analysis is not complete yet. Please specify the --pg option")
    }

    let point_groups = point_groups.unwrap();
    let samples = args.samples.unwrap_or(20);
    let iterations = args.iterations.unwrap_or(1000);

    // 3. Prepare the structure and measure it against each point group.
    let results = calc_csom(structure, args.centering_mode, args.vector, &point_groups, samples, iterations, args.full)?;

    print_csom_table(&results, &args.name);

    // 4. If requested, write the summary table (point group, deviation, rotation) to a .csv file.
    if args.table {
        write_csom_csv(&results, &args.name)?;
    }

    // 5. If requested, write the per-operation breakdown (name, matrix, deviation) to a .csv
    //    file per point group.
    if args.full {
        write_csom_details_csv(&results, &args.name)?;
    }

    Ok(())
}

/// Prepares `structure` for CSOM analysis and, for each of `point_groups`, searches the
/// Fibonacci sphere for the best-fitting symmetry axis, returning one [`CsomResult`] per
/// point group in the same order they were requested.
///
/// When `with_operations` is true, each result's `operations` is also filled in with the
/// deviation of every individual symmetry operation at the refined axis; this costs one
/// extra measurement pass per point group, so leave it false when that detail isn't needed.
pub fn calc_csom(
    structure: Structure,
    centering_mode: CenteringMode,
    centering_vector: Option<Vec<f64>>,
    point_groups: &[String],
    samples: usize,
    iterations: usize,
    with_operations: bool,
) -> Result<Vec<CsomResult>, CsomError> {

    for pg in point_groups { if !POINTGROUP_NAMES.contains(&pg.as_str()) {return Err(CsomError::WrongSpaceGroup { pg: pg.clone() })}}

    // Prepare the structure depending on centering mode.
    // (prepare_csom_structure converts --vector's raw f64s into a Vector3 itself.)
    let (csom_structure, _scale, _original_centroid) =
        prepare_csom_structure(structure, centering_mode, centering_vector);

    let mut results: Vec<CsomResult> = Vec::new();
    for point_group in point_groups {
        let (rotation_vector, deviation) = find_best_axis(samples, &csom_structure, point_group, iterations)?;
        let rotation = rotation_matrix_from_vector(rotation_vector);

        let operations = if with_operations {
            // Re-measure at the refined axis to break the overall deviation down by operation.
            let rotated_points: Vec<Vector3<f64>> = csom_structure.points.iter().map(|p| rotation * p).collect();
            point_group_operation_deviations(&rotated_points, &csom_structure.labels, point_group)?
                .into_iter()
                .map(|(name, matrix, dev)| (name.to_string(), matrix, dev))
                .collect()
        } else {
            Vec::new()
        };

        results.push(CsomResult {
            point_group: point_group.clone(),
            deviation,
            rotation,
            operations,
        });
    }

    Ok(results)
}


#[derive(Debug)]
pub enum CsomError {
    WrongSpaceGroup { pg: String},
    OptimizationFailed(String),

}

impl std::fmt::Display for CsomError  {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            CsomError::WrongSpaceGroup { pg } => {write!(f, "wrong point group name: {}", {pg})}
            CsomError::OptimizationFailed(msg) => {write!(f, "axis optimisation failed: {}", msg)}
        }
    }
}

impl std::error::Error for CsomError {}