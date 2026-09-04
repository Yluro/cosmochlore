use crate::geometry::{center_by_centroid, center_by_coordinate, center_by_first_point, normalise};
use crate::xyz::Structure;
use nalgebra::Vector3;

pub struct CsomStructure {
    pub labels: Vec<String>,
    pub points: Vec<Vector3<f64>>,
}

#[derive(Debug, Clone, clap::ValueEnum)]
pub enum CenteringMode { Auto, First, Centroid, Manual }

pub(crate) fn prepare_csom_structure(
    structure: Structure,
    centering_mode: CenteringMode,
    centering_vector: Option<Vec<f64>>) -> (CsomStructure, f64, Vector3<f64>)
{
    let mut labels: Vec<String> = Vec::new();
    let mut points: Vec<Vector3<f64>> = Vec::new();

    match structure.centre {
        Some(ref centre) => {
            labels.push(centre.label.clone());
            points.push(centre.coords)
        } ,
        None => () ,
    };

    structure.ligands.iter().for_each(|ligand| {
        labels.push(ligand.label.clone());
        points.push(ligand.coords.clone());
    });

    let striped_labels: Vec<String> = labels.iter().map(|l| strip_label(l).to_string()).collect();


    let original_centroid: Vector3<f64> = match centering_mode {
        CenteringMode::Auto => {
            if structure.centre.is_some() {
                center_by_first_point(&mut points)
            } else {
                center_by_centroid(&mut points)
            }
        },
        CenteringMode::First => center_by_first_point(&mut points),
        CenteringMode::Centroid => center_by_centroid(&mut points),
        // No error handling here because data at this point should be trusted.
        CenteringMode::Manual => {
            //let centering_vector = centering_vector.unwrap();
            let centering_vector = Vector3::from_column_slice(&centering_vector.unwrap());
            center_by_coordinate(&mut points, centering_vector)},
    };



    let scaling_factor = normalise(&mut points);

    (CsomStructure {
        labels: striped_labels,
        points
    }, scaling_factor, original_centroid)
}

pub(crate) fn strip_label(label: &String) -> String {
    let end = label
        .char_indices()
        .nth(
            if label
                .chars()
                .nth(1)
                .is_some_and(|c| c.is_ascii_lowercase())
            {
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

