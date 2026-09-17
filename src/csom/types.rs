use nalgebra::{Matrix3, Vector3};

#[derive(Debug, thiserror::Error)]
pub enum CsomError {
    #[error("wrong point group name: {pg}")]
    WrongSpaceGroup { pg: String},
    #[error("axis optimisation failed: {0}")]
    OptimizationFailed(String),
}

pub struct CsomResult {
    /// Point group analysed.
    pub point_group: String,

    /// Deviation from the ideal symmetry
    pub deviation: f64,

    /// Rotation matrix that defines the refined axis.
    pub rotation: Matrix3<f64>,

    /// Per-operation breakdown.
    pub operations: Vec<CsomOperation>,

    /// Normalization scale factor for the structure
    pub scale: f64,

    /// Centering vector used.
    pub centroid: Vector3<f64>,
}

/// One symmetry operation of a point group, measured against the structure at the refined axis.
pub struct CsomOperation {
    /// Schoenflies name of the operation.
    pub name: String,

    /// The operation's matrix, in the refined (rotated) frame.
    pub matrix: Matrix3<f64>,

    /// Deviation of the structure from this single operation.
    pub deviation: f64,

    /// `image[i]` is where the operation sends atom `i`, so it is still atom `i`'s own
    /// position and carries atom `i`'s label.
    pub image: Vec<Vector3<f64>>,

    /// `pairing[i]` is the atom whose image atom `i` was scored against. Honouring labels
    /// that is always an atom of the same element; with `--ignore` an atom can be paired
    /// with the image of a different element, which is what lowers the deviation.
    pub pairing: Vec<usize>,
}