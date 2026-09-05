pub mod linalg;
pub mod automorphism;
pub mod permutations;
pub mod bounds;
#[cfg(test)]
pub mod test_utils;
#[cfg(test)]
pub mod tests;

use crate::cli::CshmArgs;
use crate::coordinates::{points_from_reference_shape, points_from_structure};
use crate::out::{print_crab, print_cshm_table, write_cshm_csv, write_cshm_reconstructed_xyz};
use crate::shapes::{ReferenceShape, check_vertex_count};
use crate::xyz::Structure;
use crate::{shapes, xyz, yaml};
use nalgebra::Vector3;
pub(crate) use permutations::find_best_permutation;

pub struct CShMResult {
    pub name: String,
    pub symbol: String,
    pub symm: String,
    pub cshm: f64,
    pub perm: Vec<usize>,
    pub xyz: Vec<Vector3<f64>>
}


pub fn calc_cshm(reference_shapes: Vec<ReferenceShape>, problem_structure: &Structure, has_centre: bool) -> Vec<CShMResult> {

    let problem = points_from_structure(&problem_structure);
    let mut results: Vec<CShMResult> = Vec::new();


    for shape in reference_shapes {
        let mut reference = points_from_reference_shape(&shape, has_centre);
        let mut problem_copy = problem.clone();

        let (s, best_perm, reconstructed, _) = find_best_permutation(&mut reference, &mut problem_copy);

        results.push(
            CShMResult {
                name: shape.name,
                symbol: shape.symbol,
                symm: shape.symm,
                cshm: s,
                perm: best_perm,
                xyz: reconstructed,
            })
    }
    results
}

pub fn cshm_main(args: CshmArgs) -> Result<(), Box<dyn std::error::Error>> {
    // 1. Parse input .xyz file and form structure.
    let structure = xyz::parse_xyz(&args.name, args.not_centered)?;

    // 2. Get vertex count.
    let n = structure.ligands.len() as u8;

    // 3. Fetch builtin-shapes for that vertex count and index selection
    let indices = &args.shapes;
    let mut ref_shapes =  shapes::resolve_shapes(n, indices.as_deref())?;

    // 4. If user has input any shapes, add them to compare list.
    if let Some(files) = &args.user_shapes {
        let mut user_shapes: Vec<ReferenceShape> = Vec::new();

        for file in files {
            let mut s = yaml::parse_custom_shapes(file)?;
            user_shapes.append(&mut s);
        }


        for shape in user_shapes { // For each parsed shape, check its vertex count and append to list if right.
            check_vertex_count(&shape, n)?;
            ref_shapes.push(shape);
        }
    }

    // 5. Compute the CShM for each reference shape against the selected problem structures.
    let has_centre = !args.not_centered;
    let results = calc_cshm(ref_shapes, &structure, has_centre);


    // 6. Output results.
    print_cshm_table(&results, &args.name);

    // 7. If --table is passed
    if args.table { write_cshm_csv(&results, &args.name)?; }

    // 8. If --ideal is passed
    if args.ideal {
        let mut labels: Vec<String> = Vec::new();
        if !args.not_centered {
            labels.push(structure.centre.unwrap().label);
        }
        for ligand in &structure.ligands {
            labels.push(ligand.label.clone());
        }
        write_cshm_reconstructed_xyz(&args.name, &results, &labels)?;
    }

    // If --crab is passed.
    if args.crab { print_crab(); }

    Ok(())
}