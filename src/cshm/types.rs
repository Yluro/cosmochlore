use nalgebra::Vector3;

/// The continuous shape measure of a structure against one reference shape, plus what the
/// writers need to put the idealised polyhedron back onto the structure.
pub struct CShMResult {
    /// Name of the reference shape. E.g. Tetrahedron
    pub name: String,

    /// Symbol of the reference shape. E.g. TD-4
    pub symbol: String,

    /// Point group symmetry of the reference shape. E.g. Td
    pub symm: String,

    /// Continuous shape measure: 0 for an exact match, 100 at most.
    pub cshm: f64,

    /// `perm[p]` is the reference vertex matched to problem atom `p`. With a centre atom,
    /// index 0 is the centre on both sides.
    pub perm: Vec<usize>,

    /// The reference shape rotated, scaled and translated onto the problem structure, in
    /// reference-vertex order: `xyz[perm[p]]` is the ideal position of problem atom `p`.
    pub xyz: Vec<Vector3<f64>>,
}
