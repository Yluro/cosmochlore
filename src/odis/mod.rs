pub mod calc;

use crate::cli::OdisArgs;
use crate::cshm::calc_cshm;
use crate::cshm::shape_lookup::resolve_shapes;
use crate::cshm::types::CShMResult;
use crate::csom::calc_csom;
use crate::csom::prepare::CenteringMode;
use crate::csom::types::{CsomResult, OptimiserSettings};
use crate::error::Error;
use crate::xyz;
pub use calc::calculate_od;

#[derive(Debug, Clone)]
pub struct OdisResult {
    pub d_mean: f64,
    pub zeta: f64,
    pub delta: f64,
    pub sigma: f64,
    pub theta: f64,
    pub vol: f64,
    pub tau: f64,
    pub mu: f64,
}

#[derive(Debug, thiserror::Error)]
pub enum OdisError {
    #[error("structure has no central atom")]
    NoCentre,
    #[error("wrong number of points, expected: 7, found: {n}")]
    IncorrectNumberOfPoints { n: usize },
}

/// What `odis` hands to `out`: the distortion parameters and, with `--full`, the cshm and csom
/// results against the octahedral references.
pub struct OdisRun {
    pub result: OdisResult,
    pub full: Option<OdisFull>,
}

pub struct OdisFull {
    pub cshm: Vec<CShMResult>,
    pub csom: Vec<CsomResult>,
}

/// Computes the octahedral distortion parameters, plus cshm and csom with `--full`.
pub fn run(args: &OdisArgs) -> Result<OdisRun, Error> {
    // 1. Extract structure from .xyz
    let center = xyz::resolve_center(false, args.center.map(|c| c - 1));
    let structure = xyz::parse_xyz(&args.name, center)?;

    // 2. Calculate the common parameters
    let result = calculate_od(&structure)?;
    if !args.full {
        return Ok(OdisRun { result, full: None });
    }

    // 3. Calculate cshm against OC-6 and TPR-6.
    let ref_shapes = resolve_shapes(6, Some(&[2, 3]))?;
    let cshm = calc_cshm(ref_shapes, &structure, true, false);

    // 4. Calculate csom against Oh and the octahedral distortions.
    let csom_point_groups: Vec<String> = [
        "Oh",  // Ideal octahedron.
        "D4h", // Tetragonal distortion (axial elongation/compression, e.g. Jahn-Teller).
        "D3h", // Trigonal prism face twisting.
        "D3d", // Trigonal distortion (trigonal antiprismatic twist).
        "D2h", // Rhombic distortion.
        "C4v", // One ligand distinct from the other five (square-pyramidal-like).
        "C3v", // Facial substitution pattern (fac-MA3B3).
        "C2v", // cis-disubstitution pattern (cis-MA4B2).
    ]
    .iter()
    .map(|pg| pg.to_string())
    .collect();

    let csom = calc_csom(
        structure,
        CenteringMode::First,
        None,
        &csom_point_groups,
        OptimiserSettings::default(),
        false,
        false,
    )?;

    Ok(OdisRun {
        result,
        full: Some(OdisFull { cshm, csom }),
    })
}
