pub mod calc;

use crate::cli::OdisArgs;
use crate::cshm::calc_cshm;
use crate::out::{print_cshm_table, print_odis_table, write_cshm_csv};
use crate::{shapes, xyz};
pub use calc::calculate_od;

#[derive(Debug, Clone)]
pub struct OdisResult {
    pub d_mean: f64,
    pub zeta: f64,
    pub delta: f64,
    pub sigma: f64,
    //pub theta: f64,
    //pub vol: f64,
    pub tau: f64,
    pub mu: f64,
}

#[derive(Debug)]
pub enum OdisError {
    NoCentre,
    IncorrectNumberOfPoints{ n: usize },
}

impl std::fmt::Display for OdisError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            OdisError::NoCentre => {write!(f, "structure has no central atom")}
            OdisError::IncorrectNumberOfPoints{ n} => { write!(f, "wrong number of points, expected: 7, found: {}", n) },
        }
    }
}

impl std::error::Error for OdisError {}

pub fn main_odis(args: &OdisArgs) -> Result<(), Box<dyn std::error::Error>> {

    // 1. Extract structure from .xyz
    let structure = xyz::parse_xyz(&args.name, false)?; // Structure should be centered by default.

    // 2. Calculate the common parameters
    let odis_result = calculate_od(&structure)?;


    // 3. Output results
    print_odis_table(&odis_result, &args.name);

    // 4. Calculate cshm against OC-6 and TRP-6 if --full is passed.
    let n = 6; // Number of vertices
    let indices = Some(Vec::from([2, 3])); // Indices of OC and TRP
    let ref_shapes =  shapes::resolve_shapes(n, indices.as_deref())?;
    let has_centre = true;

    let cshm_results = calc_cshm(ref_shapes, &structure, has_centre);

    print_cshm_table(&cshm_results, &args.name);
    write_cshm_csv(&cshm_results, &args.name)?;

    // 5. Calculate csom against Oh, D4h, D3h



    Ok(())
}