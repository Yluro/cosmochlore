pub mod geometry;
pub mod automorphism;
pub mod permutations;
pub mod bounds;
#[cfg(test)]
pub mod test_utils;
#[cfg(test)]
pub mod tests;

use nalgebra::Vector3;
pub use permutations::find_best_permutation;

pub struct CShMResult {
    pub name: String,
    pub symbol: String,
    pub symm: String,
    pub cshm: f64,
    pub perm: Vec<usize>,
    pub xyz: Vec<Vector3<f64>>
}