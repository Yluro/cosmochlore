use crate::xyz::Structure;

/// The continuous shape measure of a structure against one reference shape.
pub struct CShMResult {
    /// Name of the reference shape. E.g. Tetrahedron
    pub name: String,

    /// Symbol of the reference shape. E.g. TD-4
    pub symbol: String,

    /// Point group symmetry of the reference shape. E.g. Td
    pub symm: String,

    /// Continuous shape measure: 0 for an exact match, 100 at most.
    pub cshm: f64,

    /// The reference shape fitted onto the structure.
    pub ideal: Option<Structure>,
}
