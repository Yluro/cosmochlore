use crate::cli::CsomArgs;
use crate::xyz;
use std::error::Error;

mod io;
mod dev;
mod tests;
mod optimize;


#[derive(Debug)]
pub enum CsomError {
    WrongSpaceGroup { pg: String},
    OptimizationFailed(String),

}

impl std::fmt::Display for CsomError  {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            CsomError::WrongSpaceGroup { pg } => {write!(f, "wrong point group name: {}", {pg})}
            CsomError::OptimizationFailed(msg) => {write!(f, "axis optimisation failed: {}", msg)}
        }
    }
}

impl std::error::Error for CsomError {}

pub fn csom_main(args: CsomArgs) -> Result<(), Box<dyn Error>> {

    // 1. Parse input .xyz file and form structure.
    let structure = xyz::parse_xyz(&args.name, !args.not_centered)?;

    // TEST COMMIT.



    Ok(())
}