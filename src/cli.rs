use clap::Parser;

#[derive(Parser, Debug)]
#[clap(name = "kosmochlor", about = "Continuous Shape and Symmetry Measurements in Rust.", author = "JSG (jose.serranog@ub.edu)", version = env!("CARGO_PKG_VERSION"))]
pub struct Cli {
    /// Path or name of the .xyz file containing the atom labels and coordinates of the problem shape.
    name: String,

    /// Treat the problem shape as a non-centered polyhedron.
    #[arg(short = 'n', long = "nc")]
    not_centered: bool,

    /// IDs of the standard reference shapes to compare the problem shape to as found in the SHAPE 2.1 user manual.
    #[arg(short = 's', long = "sh", num_args = 1..)]
    shapes: Option<Vec<u8>>,

    /// Name or path of YAML files containing user-defined reference shapes to compare the problem shape to.
    #[arg(short = 'r', long ="ref", num_args = 1..)]
    user_shapes: Option<Vec<String>>
}