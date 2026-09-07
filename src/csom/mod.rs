use crate::cli::CsomArgs;
use crate::xyz::{parse_xyz, resolve_center, Structure};
use std::error::Error;
use nalgebra::{Matrix3, Vector3};
use crate::csom::dev::point_group_operation_deviations;
use crate::csom::io::{prepare_csom_structure, CenteringMode};
use crate::csom::optimize::find_best_axis;
use crate::data::pgs::POINTGROUP_NAMES;
use crate::geometry::rotation_matrix_from_vector;
use crate::out::{print_csom_table, write_csom_csv, write_csom_details_csv, write_csom_operated_xyz, write_csom_merged_mol2};

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

    /// Per-operation breakdown:  name,  matrix, deviation, reconstructed struc
    pub operations: Vec<(String, Matrix3<f64>, f64, Vec<Vector3<f64>>)>,

    /// Normalization scale factor for the structure
    pub scale: f64,

    /// Centering vector used.
    pub centroid: Vector3<f64>,
}


pub fn csom_main(args: CsomArgs) -> Result<(), Box<dyn Error>> {

    // 1. Parse input .xyz file and form structure.
    // args.center is 1-based for the user get mapped to 0 based for parser.
    let center = resolve_center(args.not_centered, args.center.map(|c| c - 1));
    let structure = parse_xyz(&args.name, center)?;

    // 2. Fetch the desired point groups.
    let point_groups = args.point_groups;

    if point_groups.is_none() {
        todo!("Auto point-group analysis is not complete yet. Please specify the --pg option")
    }

    let point_groups = point_groups.unwrap();
    let samples = args.samples.unwrap_or(20);
    let iterations = args.iterations.unwrap_or(1000);

    // Capture the atom labels and original coordinates, in file order, before `structure` is
    // consumed by `calc_csom
    let has_centre_atom = structure.centre.is_some();
    let mut labels: Vec<String> = Vec::new();
    let mut original_coords: Vec<Vector3<f64>> = Vec::new();
    if let Some(ref centre) = structure.centre {
        labels.push(centre.label.clone());
        original_coords.push(centre.coords);
    }
    labels.extend(structure.ligands.iter().map(|l| l.label.clone()));
    original_coords.extend(structure.ligands.iter().map(|l| l.coords));

    // The per-operation breakdown (and the operated coordinates it now carries) is only
    // computed when either --full or --operated actually needs it.
    let with_operations = args.full || args.operated;

    // 3. Prepare the structure and measure it against each point group.
    let results = calc_csom(structure, args.centering_mode, args.vector, &point_groups, samples, iterations, with_operations, args.ignore_labels)?;

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

    // 6. If requested, write the operated coordinates (recovered rotation + each symmetry
    //    operation applied to the original structure) to a .xyz file per point group, plus a
    //    merged .mol2 overlaying every image in one 3D structure for viewers such as Mercury.
    if args.operated {
        write_csom_operated_xyz(&results, &args.name, &labels, &original_coords)?;
        write_csom_merged_mol2(&results, &args.name, &labels, &original_coords, has_centre_atom)?;
    }

    Ok(())
}

/// Prepares `structure` for CSOM analysis and, for each of `point_groups`.
///
/// When `with_operations` is true, each result's `operations` is also filled in with the
/// deviation of every individual symmetry operation at the refined axis
///
/// When `ignore_labels` is true, atom labels are not used to restrict the permutation search.
/// Returns one `CsomResult` per point group analysed.
pub fn calc_csom(
    structure: Structure,
    centering_mode: CenteringMode,
    centering_vector: Option<Vec<f64>>,
    point_groups: &[String],
    samples: usize,
    iterations: usize,
    with_operations: bool,
    ignore_labels: bool,
) -> Result<Vec<CsomResult>, CsomError> {

    for pg in point_groups { if !POINTGROUP_NAMES.contains(&pg.as_str()) {return Err(CsomError::WrongSpaceGroup { pg: pg.clone() })}}

    // Prepare the structure depending on centering mode.
    // (prepare_csom_structure converts --vector's raw f64s into a Vector3 itself.)
    let (csom_structure, scale, centroid) =
        prepare_csom_structure(structure, centering_mode, centering_vector);

    let mut results: Vec<CsomResult> = Vec::new();
    for point_group in point_groups {
        let (rotation_vector, deviation) = find_best_axis(samples, &csom_structure, point_group, iterations, ignore_labels)?;
        let rotation = rotation_matrix_from_vector(rotation_vector);

        let operations = if with_operations {
            // Re-measure at the refined axis to break the overall deviation down by operation.
            let rotated_points: Vec<Vector3<f64>> = csom_structure.points.iter().map(|p| rotation * p).collect();
            point_group_operation_deviations(&rotated_points, &csom_structure.labels, point_group, ignore_labels, csom_structure.has_centre)?
                .into_iter()
                .map(|(name, matrix, dev, xyz)| (name.to_string(), matrix, dev, xyz))
                .collect()
        } else {
            Vec::new()
        };

        results.push(CsomResult {
            point_group: point_group.clone(),
            deviation,
            rotation,
            operations,
            scale,
            centroid,
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