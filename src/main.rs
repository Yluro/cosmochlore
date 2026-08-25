mod xyz;
mod cli;
mod yaml;
mod shapes;
mod data;
mod cshm;
mod coordinates;
mod out;

use std::time::Instant;
use clap::Parser;
use crate::cli::{Cli};
use crate::out::{welcome_msg, MeasureResult, output_table};
use crate::coordinates::{points_from_reference_shape, points_from_structure};
use crate::cshm::find_best_permutation;
use crate::shapes::{check_vertex_count, ReferenceShape};

fn main() {
    let main_start = Instant::now();
    welcome_msg();
    let args = Cli::parse();

    // 1. Parse input .xyz file and form structure.
    let structure = match xyz::parse_xyz(&args.name, args.not_centered){
        Ok(s) => s,
        Err(err) => panic!("{}", err)
    };

    // 2. Get vertex count.
    let n = structure.ligands.len() as u8;

    // 3. Fetch builtin-shapes for that vertex count and index selection
    let indices = &args.shapes;
    let mut ref_shapes =  match shapes::resolve_shapes(n, indices.as_deref()) {
        Ok(s) => s,
        Err(err) => panic!("{}", err)
    };

    // 4. If user has input any shapes, add them to compare list.
    if let Some(files) = args.user_shapes {
        let mut user_shapes: Vec<ReferenceShape> = Vec::new();
        for file in files {
            match yaml::parse_custom_shapes(&file) {
                Ok(mut s) => user_shapes.append(&mut s),
                Err(err) => panic!("{}", err)
            }
        }

        for shape in user_shapes { // For each parsed shape, check its vertex count and append to list if right.
            match check_vertex_count(&shape, n) {
                Ok(()) => ref_shapes.push(shape),
                Err(err) => panic!("{}", err),
            }
        }
    }

    // 5. Compute the CShM for each reference shape against the selected problem structures.
    let has_centre = !args.not_centered;
    let problem = points_from_structure(&structure);

    let mut results: Vec<MeasureResult> = Vec::new();

    for shape in ref_shapes {
        let mut reference = points_from_reference_shape(&shape, has_centre);
        let mut problem_copy = problem.clone();

        let (s, _, rot_mat) = find_best_permutation(&mut reference, &mut problem_copy);

        results.push(
            MeasureResult {
                name: shape.name,
                symbol: shape.symbol,
                symm: shape.symm,
                cshm: s,
                rot_mat: rot_mat,
            }
        );
    }

    output_table(&results);
    let main_time = main_start.elapsed();
    println!("Calculations done in {:?}", main_time);
}
