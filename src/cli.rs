use clap::{Args, Parser, Subcommand};


#[derive(Parser, Debug)]
#[clap(name = "cosmochlore", about = env!("CARGO_PKG_DESCRIPTION"), author = env!("CARGO_PKG_AUTHORS"), version = env!("CARGO_PKG_VERSION"))]
pub struct Cli {
    #[clap(subcommand)]
    pub command: Command,
}

#[derive(Subcommand, Debug)]
pub enum Command {
    /// Continuous Shape Measures
    Cshm(CshmArgs),
    /// Continuous Symmetry Operation Measures
    Csom(CsomArgs),
    /// Octahedral Distortion Analysis
    Odis(OdisArgs),
}

#[derive(Args, Debug)]
pub struct CshmArgs {
    /// Path or name of the .xyz file containing the atom labels and coordinates of the problem shape.
    pub name: String,

    /// Treat the problem shape as a non-centered.
    #[arg(short = 'n', long = "nc")]
    pub not_centered: bool,

    /// Indices of the standard reference shapes to compare the problem shape to as found in the user manual.
    #[arg(short = 's', long = "sh", num_args = 1..)]
    pub shapes: Option<Vec<usize>>,

    /// Name or path of YAML files containing user-defined reference shapes to compare the problem shape to.
    #[arg(short = 'r', long ="ref", num_args = 1..)]
    pub user_shapes: Option<Vec<String>>,

    /// Write the output table to a .csv file.
    #[arg(short = 't', long = "table")]
    pub table: bool,

    /// Write a reconstructed version of the idealised polyhedra to a .xyz file.
    #[arg(short = 'i', long = "ideal")]
    pub ideal: bool,

    /// Easter-egg
    #[arg(short = 'c', long = "crab", hide = true)]
    pub crab: bool,
}

#[derive(Args, Debug)]
pub struct  CsomArgs {
    /// Path or name of the .xyz file containing the atom labels and coordinates of the problem shape.
    pub name: String,
    
    /// Treat the structure as non-centered
    #[arg(short = 'n', long = "nc")]
    pub not_centered: bool,

    /// Space groups to measure in Schoenflies notation.
    #[arg(short = 'p', long = "pg", num_args = 1..)]
    pub point_groups: Option<Vec<String>>,
}

#[derive(Args, Debug)]
pub struct OdisArgs {
    /// Path or name of the .xyz file containing the atom labels and coordinates of the problem shape.
    pub name: String,

    /// Full analysis of the octahedron. Including CShM and CSoM values.
    #[arg(short = 'f', long = "full")]
    pub full: bool,

    /// Write the output table to a .csv file.
    #[arg(short = 't', long = "table")]
    pub table: bool,
}