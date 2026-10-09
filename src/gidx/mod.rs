pub mod calc;
#[cfg(test)]
mod tests;

use crate::cli::GidxArgs;
use crate::error::Error;
use crate::out::*;
use crate::xyz;
pub use calc::calculate_gidx;

/// A square face of an eight-coordinate centre: its diagonals subtend `theta >= phi`
/// degrees and give `tau8`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SquareFace {
    pub theta: f64,
    pub phi: f64,
    pub tau8: f64,
}

/// Angles (degrees) and indices for the coordination number of the centre.
#[derive(Debug, Clone, PartialEq)]
pub enum GidxResult {
    /// `tau4` and `tau4_prime` from the two largest angles, `beta >= alpha`.
    Four {
        alpha: f64,
        beta: f64,
        tau4: f64,
        tau4_prime: f64,
    },
    /// `tau5` from the two largest angles, `beta >= alpha`.
    Five { alpha: f64, beta: f64, tau5: f64 },
    /// `tau6` from the trans angles `alpha1 >= alpha2 >= alpha3`, and cosmochlore's own
    /// `tau6_prime` from the five smallest angles `gammas`.
    Six {
        alpha1: f64,
        alpha2: f64,
        alpha3: f64,
        tau6: f64,
        gammas: [f64; 5],
        tau6_prime: f64,
    },
    /// `tau8` per square face, the A and B plane angles `dihedral_a` and `dihedral_b`, `tau8_prime`,
    /// and `cube_delta` = tau8' − tau8 (larger tau8). Faces near 1 and 0 point to a bicapped prism.
    Eight {
        faces: [SquareFace; 2],
        dihedral_a: f64,
        dihedral_b: f64,
        tau8_prime: f64,
        cube_delta: f64,
    },
}

#[derive(Debug, thiserror::Error)]
pub enum GidxError {
    #[error("structure has no central atom")]
    NoCentre,
    #[error("wrong number of ligands, expected: 4, 5, 6 or 8, found: {n}")]
    UnsupportedCoordination { n: usize },
    #[error(
        "a ligand-centre-ligand angle is undefined (NaN): check that no ligand sits on the central atom"
    )]
    UndefinedAngle,
    #[error(
        "a plane through the central atom and two ligands is undefined: the two ligands lie on one line with it"
    )]
    UndefinedPlane,
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
