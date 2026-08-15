use clap::Parser;

#[derive(Parser, Debug)]
#[clap(name = "kosmochlor", about = "Continuous Shape and Symmetry Measurements in Rust.")]
pub struct Cli {
    /// Path or name of the .xyz file with the problem shape.
    name: String,

    /// Treat the problem shape as a non-centered polyhedron.
    #[arg(short = 'n', long = "nc")]
    not_centered: bool,

    /// Standard reference shapes to compare to.
    #[arg(short = 's', long = "sh", num_args = 1..)]
    shapes: Option<Vec<u8>>,

    /// Name or path of YAML files containing user-defined reference shapes to compare to.
    #[arg(short = 'r', long ="ref", num_args = 1..)]
    user_shapes: Option<Vec<String>>
}