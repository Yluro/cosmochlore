pub mod geometry;
pub mod automorphism;
pub mod permutations;
pub mod bounds;
#[cfg(test)]
pub mod test_utils;
#[cfg(test)]
pub mod tests;

use nalgebra::Vector3;
pub(crate) use permutations::find_best_permutation;
use crate::coordinates::{points_from_reference_shape, points_from_structure};
use crate::shapes::ReferenceShape;
use crate::xyz::Structure;

pub struct CShMResult {
    pub name: String,
    pub symbol: String,
    pub symm: String,
    pub cshm: f64,
    pub perm: Vec<usize>,
    pub xyz: Vec<Vector3<f64>>
}


pub fn calc_cshm(reference_shapes: Vec<ReferenceShape>, problem_structure: &Structure, has_centre: bool) -> Vec<CShMResult> {

    let problem = points_from_structure(&problem_structure);
    let mut results: Vec<CShMResult> = Vec::new();


    for shape in reference_shapes {
        let mut reference = points_from_reference_shape(&shape, has_centre);
        let mut problem_copy = problem.clone();

        let (s, best_perm, reconstructed, _) = find_best_permutation(&mut reference, &mut problem_copy);

        results.push(
            CShMResult {
                name: shape.name,
                symbol: shape.symbol,
                symm: shape.symm,
                cshm: s,
                perm: best_perm,
                xyz: reconstructed,
            })
    }
    results
}