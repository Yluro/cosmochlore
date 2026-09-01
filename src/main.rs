mod xyz;
mod cli;
mod yaml;
mod shapes;
mod data;
mod cshm;
mod csom;
mod odis;
mod coordinates;
mod out;


use crate::cli::{Cli, Command, CshmArgs, OdisArgs};
use crate::coordinates::{points_from_reference_shape, points_from_structure};
use crate::out::{print_odis_table, print_cshm_table, welcome_msg, write_cshm_csv, write_cshm_reconstructed_xyz, print_crab};
use crate::shapes::{ReferenceShape, check_vertex_count};
use clap::Parser;
use std::time::Instant;

use crate::cshm::{find_best_permutation, CShMResult};
use crate::odis::{calculate_od, ODResult};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let main_start = Instant::now();
    welcome_msg();
    let args = Cli::parse();


    let run = match args.command {
        Command::Cshm(cshm_args) => { main_cshm(&cshm_args) },
        Command::Csom(_csom_args) => { unimplemented!() },
        Command::Odis(odis_args) => { main_odis(&odis_args) },
    };

    if let Err(err) = run {
        eprintln!("Error: {err}");
        std::process::exit(1);
    }

    let main_time = main_start.elapsed();
    println!("Program finished in {:?}", main_time);
    //if args.crab {print_crab()}
    run
}

fn main_cshm(args: &CshmArgs) -> Result<(), Box<dyn std::error::Error>> {
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
    let problem = points_from_structure(&structure);

    let mut results: Vec<CShMResult> = Vec::new();

    for shape in ref_shapes {
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
            }
        );
    }


    // 6. Output results.
    print_cshm_table(&results, &args.name);

    // 7. If --table is passed

    print_cshm_table(&results, &args.name);

    if args.table {
        write_cshm_csv(&args.name, &results)?;
    }

    // If --ideal is passed
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


fn main_odis(args: &OdisArgs) -> Result<(), Box<dyn std::error::Error>> {

    // 1. Extract structure from .xyz
    let structure = xyz::parse_xyz(&args.name, false)?; // Structure should be centered by default.

    // 2. Calculate the common parameters
    let result = calculate_od(&structure)?;


    // 3. Output results
    print_odis_table(&result, &args.name);

    // 4. Calculate chsm against OC-6 and TRP-6 if --full is passed.
    let n = 6; // Number of vertices
    let indices = Some(Vec::from([2, 3])); // Indices of OC and TRP
    let ref_shapes =  shapes::resolve_shapes(n, indices.as_deref())?;
    let has_centre = true;
    let problem = points_from_structure(&structure);
    let mut cshm_results: Vec<CShMResult> = Vec::new();
    for shape in ref_shapes {
        let mut reference = points_from_reference_shape(&shape, has_centre);
        let mut problem_copy = problem.clone();

        let (s, best_perm, reconstructed, _) = find_best_permutation(&mut reference, &mut problem_copy);

        cshm_results.push(
            CShMResult {
                name: shape.name,
                symbol: shape.symbol,
                symm: shape.symm,
                cshm: s,
                perm: best_perm,
                xyz: reconstructed,
            });
    }

    print_cshm_table(&cshm_results, &args.name);
    write_cshm_csv(&args.name, &cshm_results)?;



    Ok(())
}


