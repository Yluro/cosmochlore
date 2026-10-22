use crate::xyz::Structure;
use nalgebra::{Matrix3, Vector3};

/// Knobs of the symmetry-axis search. `Default` is the one place the values live: the `csom`
/// CLI defaults are read from it, and `odis --full` runs the same search with the same values.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OptimiserSettings {
    /// Number of Fibonacci-sphere samples used to seed the search.
    pub seeds: usize,
    /// Maximum number of Nelder-Mead iterations spent refining each seed.
    pub iterations: usize,
    /// Nelder-Mead convergence tolerance: the standard deviation of the simplex's deviations.
    pub tolerance: f64,
}

impl Default for OptimiserSettings {
    fn default() -> Self {
        Self {
            seeds: 20,
            iterations: 200,
            tolerance: 1e-6,
        }
    }
}

/// Centering mode enumerator.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Centering {
    /// The centre atom when the structure has one, else the centroid.
    Auto,
    First,
    Centroid,
    /// A given point.
    Manual(Vector3<f64>),
}

/// Everything `calc_csom` needs besides the structure and the point groups.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CsomSettings {
    pub centering: Centering,
    pub search: OptimiserSettings,
    /// Match atoms regardless of their labels.
    pub ignore_labels: bool,
    /// Also fill every result's `operations` and `operated`.
    pub operations: bool,
}

impl Default for CsomSettings {
    fn default() -> Self {
        Self {
            centering: Centering::Auto,
            search: OptimiserSettings::default(),
            ignore_labels: false,
            operations: false,
        }
    }
}

#[derive(Debug, thiserror::Error)]
pub enum CsomError {
    #[error("wrong point group name: {pg}")]
    WrongSpaceGroup { pg: String },
    #[error("manual centering needs --vector X Y Z")]
    MissingVector,
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

    /// Per-operation breakdown, with `--details` or `--operated`.
    pub operations: Vec<CsomOperationResult>,

    /// Images of operated structures with `--details` or `--operated`.
    pub operated: Option<Vec<Structure>>,
}

/// One symmetry operation of a point group, measured against the structure at the refined axis.
pub struct CsomOperationResult {
    /// Schoenflies name of the operation.
    pub name: &'static str,

    /// The operation's matrix, in the refined (rotated) frame.
    pub matrix: Matrix3<f64>,

    /// Deviation of the structure from this single operation.
    pub deviation: f64,

    /// `image[i]` is where the operation sends atom `i`, so it is still atom `i`'s own
    /// position and carries atom `i`'s label.
    pub image: Vec<Vector3<f64>>,

    /// `pairing[i]` is the atom whose image atom `i` was scored against.
    pub pairing: Vec<usize>,
}
