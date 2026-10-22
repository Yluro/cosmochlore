use crate::cli::CsomArgs;
use crate::csom::deviation::point_group_operation_deviations;
use crate::csom::optimize::search_best_axis;
use crate::csom::prepare::prepare_csom_structure;
use crate::data::pgs::POINTGROUP_NAMES;
use crate::error::Error;
use crate::geometry::rotation_matrix_from_vector;
use crate::xyz::{Structure, parse_xyz, resolve_center};
use nalgebra::Vector3;
use types::{CsomError, CsomOperationResult, CsomResult, CsomSettings};

mod assignment;
mod deviation;
mod optimize;
pub(crate) mod prepare;
mod seeding;
#[cfg(test)]
mod tests;
pub(crate) mod types;

pub fn run(args: &CsomArgs) -> Result<Vec<CsomResult>, Error> {
    // 1. Parse input .xyz file and form structure.
    // args.center is 1-based for the user get mapped to 0 based for parser.
    let center = resolve_center(args.not_centered, args.center.map(|c| c - 1));
    let structure = parse_xyz(&args.name, center)?;

    // 2. Fetch the desired point groups. When none are given, analyse against every supported
    //    point group.
    let point_groups = args
        .point_groups
        .clone()
        .unwrap_or_else(|| POINTGROUP_NAMES.iter().map(|pg| pg.to_string()).collect());

    // 3. Measure it against each point group.
    Ok(calc_csom(&structure, &point_groups, &args.settings()?)?)
}

/// Centres `structure` and, for each of `point_groups`, searches for the orientation of minimum
/// deviation. With `settings.operations`, each result also gets the deviation of every operation
/// at the refined axis (`operations`) and the input with its image under each one (`operated`).
pub fn calc_csom(
    structure: &Structure,
    point_groups: &[String],
    settings: &CsomSettings,
) -> Result<Vec<CsomResult>, CsomError> {
    for pg in point_groups {
        if !POINTGROUP_NAMES.contains(&pg.as_str()) {
            return Err(CsomError::WrongSpaceGroup { pg: pg.clone() });
        }
    }

    // 4. Move the structure to the measuring frame (centring and normalizing)
    let (csom_structure, scale, centroid) = prepare_csom_structure(structure, settings.centering);

    let mut results: Vec<CsomResult> = Vec::new();
    for point_group in point_groups {
        // 5. Find the best orientation to apply the symmetry operations.
        let (rotation_vector, deviation) = search_best_axis(
            &csom_structure,
            point_group,
            settings.search,
            settings.ignore_labels,
        )?;
        let rotation = rotation_matrix_from_vector(rotation_vector);

        // If per-operation deviations are needed, use the rotation vector from the optimiser and
        // rerun the symmetry elements.
        let operations: Vec<CsomOperationResult> = if settings.operations {
            // Re-measure at the refined axis to break the overall deviation down by operation.
            let rotated_points: Vec<Vector3<f64>> =
                csom_structure.points.iter().map(|p| rotation * p).collect();
            point_group_operation_deviations(
                &rotated_points,
                &csom_structure.groups,
                point_group,
                settings.ignore_labels,
                csom_structure.has_centre,
            )?
        } else {
            Vec::new()
        };

        // Each image back in the input frame: `rotation` is orthogonal, so its transpose undoes
        // the alignment.
        let operated = settings.operations.then(|| {
            let to_input = |p: &Vector3<f64>| (rotation.transpose() * p) / scale + centroid;
            let mut blocks = vec![structure.clone()];
            blocks.extend(operations.iter().map(|op| {
                let coords: Vec<Vector3<f64>> = op.image.iter().map(to_input).collect();
                structure.with_coords(&coords)
            }));
            blocks
        });

        results.push(CsomResult {
            point_group: point_group.clone(),
            deviation,
            rotation,
            operations,
            operated,
        });
    }

    Ok(results)
}
