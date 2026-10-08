pub mod calc;
#[cfg(test)]
mod tests;

use crate::cli::GidxArgs;
use crate::error::Error;
use crate::out::*;
use crate::xyz;
pub use calc::calculate_gidx;

/// The angles and geometry indices that apply to the coordination number of the centre. All
/// angles are in degrees.
#[derive(Debug, Clone, PartialEq)]
pub enum GidxResult {
    /// Four-coordinate centre: `tau4` (Yang 2007) and `tau4_prime` (Okuniewski 2015), built from
    /// the two largest ligand-centre-ligand angles, `beta >= alpha`. Both indices are 0 for
    /// square planar and 1 for tetrahedral.
    Four {
        alpha: f64,
        beta: f64,
        tau4: f64,
        tau4_prime: f64,
    },
    /// Five-coordinate centre: `tau5` (Addison 1984), built from the two largest angles,
    /// `beta >= alpha`. 0 for square pyramidal and 1 for trigonal bipyramidal.
    Five { alpha: f64, beta: f64, tau5: f64 },
    /// Six-coordinate centre: `tau6` (Stoeckli-Evans 2025), built from the three principal
    /// (trans) angles, `alpha1 >= alpha2 >= alpha3`. 0 for octahedral and 0.75 for the
    /// trigonal prism (when its trans angles are 135°).
    Six {
        alpha1: f64,
        alpha2: f64,
        alpha3: f64,
        tau6: f64,
    },
}

#[derive(Debug, thiserror::Error)]
pub enum GidxError {
    #[error("structure has no central atom")]
    NoCentre,
    #[error("wrong number of ligands, expected: 4, 5 or 6, found: {n}")]
    UnsupportedCoordination { n: usize },
    #[error(
        "a ligand-centre-ligand angle is undefined (NaN): check that no ligand sits on the central atom"
    )]
    UndefinedAngle,
}

pub fn main_gidx(args: GidxArgs) -> Result<(), Error> {
    // 1. Extract structure from .xyz
    let center = xyz::resolve_center(false, args.center.map(|c| c - 1));
    let structure = xyz::parse_xyz(&args.name, center)?;

    // 2. Calculate the geometry indices for the structure's coordination number
    let gidx_result = calculate_gidx(&structure)?;

    // 3. Output results
    print_gidx_table(&gidx_result, &args.name);
    if args.table {
        write_gidx_csv(&gidx_result, &args.name)?
    };

    Ok(())
}
