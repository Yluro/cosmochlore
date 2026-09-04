use crate::cli::CsomArgs;
use crate::xyz;
use std::error::Error;
use nalgebra::{Matrix3};
use crate::csom::io::prepare_csom_structure;
use crate::csom::optimize::find_best_axis;
use crate::data::pgs::POINTGROUP_NAMES;
use crate::geometry::rotation_matrix_from_vector;
use crate::out::print_csom_table;

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
}


pub fn csom_main(args: CsomArgs) -> Result<(), Box<dyn Error>> {

    // 1. Parse input .xyz file and form structure.
    let structure = xyz::parse_xyz(&args.name, !args.not_centered)?;

    // 2. Prepare the structure depending on centering mode.
    //    (prepare_csom_structure converts --vector's raw f64s into a Vector3 itself.)
    let (csom_structure, _scale, _original_centroid) =
        prepare_csom_structure(structure, args.centering_mode, args.vector);

    // 3. Fetch the desired point groups.
    let point_groups = args.point_groups;

    if point_groups.is_none() {
        todo!("Auto point-group analysis is not complete yet. Please specify the --pg option")
    }

    let point_groups = point_groups.unwrap();
    for pg in &point_groups { if !POINTGROUP_NAMES.contains(&pg.as_str()) {return Err(CsomError::WrongSpaceGroup { pg: pg.clone() }.into())}}

    let samples = args.samples.unwrap_or(20);

    // 4. For each requested point group, search the Fibonacci sphere for the
    //    best-fitting symmetry axis and record its deviation + rotation.
    let mut results: Vec<CsomResult> = Vec::new();
    for point_group in point_groups {
        let (rotation_vector, deviation) = find_best_axis(samples, &csom_structure, &point_group)?;

        let result = CsomResult {
            point_group,
            deviation,
            rotation: rotation_matrix_from_vector(rotation_vector),
        };

        // 5. Add the result of the point group to the report.
        results.push(result);

    }
    print_csom_table(&results, &args.name);
    Ok(())
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