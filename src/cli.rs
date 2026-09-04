use clap::{Args, Parser, Subcommand};
use crate::csom::io::CenteringMode;

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

    /// Centering mode. Defaults to Auto. If manual is passed. A centering vector is required.
    ///
    /// Auto mode centers by the first atom if the structure is centered and
    /// centers by the centroid if the structure is not centered.
    #[arg(short = 'c', long = "center", value_enum, default_value = "auto")]
    pub centering_mode: CenteringMode,

    /// Centering vector for manual centering.
    #[arg(short = 'u', long = "vector", num_args = 3, requires_if("manual", "centering_mode"))]
    pub vector: Option<Vec<f64>>,

    /// Write per-operation deviation details (name, matrix, deviation) to a .csv file for
    /// each point group analysed.
    #[arg(short = 'f', long = "full")]
    pub full: bool,

    /// Write the summary output table (point group, deviation, rotation) to a .csv file.
    #[arg(short = 't', long = "table")]
    pub table: bool,

    /// Number of samples taken f the Fibonacci sphere.
    #[arg(short = 's', long = "samples", default_value = "20")]
    pub samples: Option<usize>,

    /// Maximum number of iterations for Nelder-Mead optimization of the z-axis.
    #[arg(short = 'i', long = "iterations", default_value = "200")]
    pub iterations: Option<usize>,
    
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