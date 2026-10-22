pub mod linalg;
// Only consumer is tests::bnb_automorphisms_matches_naive; the CShM search itself
// prunes by bounds, not by automorphisms.
#[cfg(test)]
pub mod automorphism;
pub mod bounds;
pub mod permutations;
pub mod shape_lookup;
#[cfg(test)]
pub mod test_utils;
#[cfg(test)]
pub mod tests;
pub(crate) mod types;

use crate::cli::CshmArgs;
use crate::data::standard_shapes::ReferenceShape;
use crate::error::Error;
use crate::xyz::Structure;
use crate::{xyz, yaml};
use nalgebra::Vector3;
pub(crate) use permutations::find_best_permutation;
use shape_lookup::{check_vertex_count, resolve_shapes};
use types::CShMResult;

pub fn calc_cshm(
    reference_shapes: Vec<ReferenceShape>,
    problem_structure: &Structure,
    has_centre: bool,
    with_ideal: bool, // Also recover the rotation matrix to rebuild the ideal shape
) -> Vec<CShMResult> {
    let problem: Vec<Vector3<f64>> = problem_structure
        .atoms()
        .iter()
        .map(|atom| atom.coords)
        .collect();
    let mut results: Vec<CShMResult> = Vec::new();

    for shape in reference_shapes {
        let mut reference = shape.points(has_centre);
        let mut problem_copy = problem.clone();

        let (s, best_perm, reconstructed, _) =
            find_best_permutation(&mut reference, &mut problem_copy, has_centre);

        let ideal =
            with_ideal.then(|| ideal_structure(problem_structure, &best_perm, &reconstructed));
        results.push(CShMResult {
            name: shape.name,
            symbol: shape.symbol,
            symm: shape.symm,
            cshm: s,
            ideal,
        })
    }
    results
}

/// `problem` with every atom moved to its matched vertex of the fitted shape. `perm[p]` is the
/// vertex of `fitted` matched to atom `p` of `problem.atoms()`.
fn ideal_structure(problem: &Structure, perm: &[usize], fitted: &[Vector3<f64>]) -> Structure {
    let coords: Vec<Vector3<f64>> = perm.iter().map(|&vertex| fitted[vertex]).collect();
    problem.with_coords(&coords)
}

/// Measures the input structure against the selected reference shapes.
pub fn run(args: &CshmArgs) -> Result<Vec<CShMResult>, Error> {
    // 1. Parse input .xyz file and form structure.
    // args.center is 1-based for the user; the parser wants a 0-based table index.
    let center = xyz::resolve_center(args.not_centered, args.center.map(|c| c - 1));
    let structure = xyz::parse_xyz(&args.name, center)?;

    // 2. Get vertex count.
    let n = structure.ligands.len() as u8;

    // 3. Fetch builtin-shapes for that vertex count and index selection
    let mut ref_shapes = resolve_shapes(n, args.shapes.as_deref())?;

    // 4. If user has input any shapes, add them to compare list.
    for file in args.user_shapes.iter().flatten() {
        for shape in yaml::parse_custom_shapes(file)? {
            check_vertex_count(&shape, n)?;
            ref_shapes.push(shape);
        }
    }

    // 5. Compute the CShM for each reference shape against the selected problem structures.
    Ok(calc_cshm(
        ref_shapes,
        &structure,
        center.is_some(),
        args.ideal,
    ))
}
