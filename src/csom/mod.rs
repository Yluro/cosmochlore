use crate::cli::CsomArgs;
use crate::csom::deviation::point_group_operation_deviations;
use crate::csom::prepare::{CenteringMode, prepare_csom_structure};
use crate::csom::optimize::search_best_axis;
use crate::data::pgs::POINTGROUP_NAMES;
use crate::error::Error;
use crate::geometry::rotation_matrix_from_vector;
use crate::out::{print_csom_table, write_csom_csv, write_csom_details_csv, write_csom_merged_mol2, write_csom_operated_xyz};
use crate::xyz::{Structure, parse_xyz, resolve_center};
use nalgebra::Vector3;
use types::{CsomError, CsomOperation, CsomResult};

mod assignment;
mod deviation;
pub(crate) mod prepare;
#[cfg(test)]
mod tests;
mod optimize;
pub(crate) mod types;


pub fn csom_main(args: CsomArgs) -> Result<(), Error> {

    // 1. Parse input .xyz file and form structure.
    // args.center is 1-based for the user get mapped to 0 based for parser.
    let center = resolve_center(args.not_centered, args.center.map(|c| c - 1));
    let structure = parse_xyz(&args.name, center)?;

    // 2. Fetch the desired point groups. When none are given, analyse against every supported
    //    point group.
    let point_groups = args.point_groups.unwrap_or_else(|| POINTGROUP_NAMES.iter().map(|pg| pg.to_string()).collect());

    // Capture the atom labels and original coordinates, in file order, before `structure` is
    // consumed by `calc_csom`.
    let has_centre_atom = structure.centre.is_some();
    let atoms = structure.atoms();
    let labels: Vec<String> = atoms.iter().map(|a| a.label.clone()).collect();
    let original_coords: Vec<Vector3<f64>> = atoms.iter().map(|a| a.coords).collect();

    // The per-operation breakdown (and the operated coordinates it now carries) is only
    // computed when either --full or --operated actually needs it.
    let with_operations = args.full || args.operated;

    // 3. Prepare the structure and measure it against each point group.
    let results = calc_csom(structure, args.centering_mode, args.vector, &point_groups, args.seeds, args.iterations, args.tolerance, with_operations, args.ignore_labels)?;

    print_csom_table(&results, &args.name);

    // 4. If requested, write the summary table (point group, deviation, rotation) to a .csv file.
    if args.table {
        write_csom_csv(&results, &args.name)?;
    }

    // 5. If requested, write the per-operation breakdown (name, matrix, deviation) to a .csv
    //    file per point group.
    if args.full {
        write_csom_details_csv(&results, &args.name, &labels)?;
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
    tolerance: f64,
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
        let (rotation_vector, deviation) = search_best_axis(samples, &csom_structure, point_group, iterations, tolerance, ignore_labels)?;
        let rotation = rotation_matrix_from_vector(rotation_vector);

        let operations: Vec<CsomOperation> = if with_operations {
            // Re-measure at the refined axis to break the overall deviation down by operation.
            let rotated_points: Vec<Vector3<f64>> = csom_structure.points.iter().map(|p| rotation * p).collect();
            point_group_operation_deviations(&rotated_points, &csom_structure.groups, point_group, ignore_labels, csom_structure.has_centre)?
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