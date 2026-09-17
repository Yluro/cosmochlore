use crate::csom::assignment::group_by_label;
use crate::geometry::{center_by_centroid, center_by_coordinate, center_by_first_point, normalise};
use crate::xyz::Structure;
use nalgebra::Vector3;

pub struct CsomStructure {
    pub points: Vec<Vector3<f64>>,
    pub has_centre: bool,

    /// Stripped element labels partitioned into same-element index groups, precomputed once
    /// here since it's invariant across the many calls a csom search makes (see
    /// `group_by_label`).
    pub groups: Vec<Vec<usize>>,
}

#[derive(Debug, Clone, clap::ValueEnum)]
pub enum CenteringMode {
    Auto,
    First,
    Centroid,
    Manual,
}

pub(crate) fn prepare_csom_structure(
    structure: Structure,
    centering_mode: CenteringMode,
    centering_vector: Option<Vec<f64>>,
) -> (CsomStructure, f64, Vector3<f64>) {
    let atoms = structure.atoms();
    let striped_labels: Vec<String> = atoms.iter().map(|a| strip_label(&a.label)).collect();
    let mut points: Vec<Vector3<f64>> = atoms.iter().map(|a| a.coords).collect();

    let (original_centroid, has_centre): (Vector3<f64>, bool) = match centering_mode {
        CenteringMode::Auto => {
            if structure.centre.is_some() {
                (center_by_first_point(&mut points), true)
            } else {
                (center_by_centroid(&mut points), false)
            }
        }
        CenteringMode::First => (center_by_first_point(&mut points), true),
        CenteringMode::Centroid => (center_by_centroid(&mut points), false),
        // No error handling here because data at this point should be trusted.
        CenteringMode::Manual => {
            //let centering_vector = centering_vector.unwrap();
            let centering_vector = Vector3::from_column_slice(&centering_vector.unwrap());
            (center_by_coordinate(&mut points, centering_vector), false)
        }
    };

    let scaling_factor = normalise(&mut points);
    let groups = group_by_label(&striped_labels);

    (
        CsomStructure {
            points,
            has_centre,
            groups,
        },
        scaling_factor,
        original_centroid,
    )
}

pub(crate) fn strip_label(label: &str) -> String {
    let end = label
        .char_indices()
        .nth(
            if label.chars().nth(1).is_some_and(|c| c.is_ascii_lowercase()) {
                2
            } else {
                1
            },
        )
        .map_or(label.len(), |(i, _)| i);

    label[..end].to_string()
}

pub(crate) fn strip_all_labels(labels: &[String]) -> Vec<String> {
    labels.iter().map(|l| strip_label(l)).collect()
}
