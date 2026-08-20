mod xyz;
mod cli;
mod yaml;
mod shapes;
mod data;
mod cshm;
mod points;

use clap::Parser;
use crate::cli::{Cli, welcome_msg};
use crate::shapes::{check_vertex_count, ReferenceShape};

fn main() {

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

    // 5. Print the resulting shape names
    for shape in ref_shapes{
        println!("{}, {}", shape.name, shape.symbol);
    }
}
