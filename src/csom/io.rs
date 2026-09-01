use nalgebra::Vector3;
use crate::csom::center::*;
use crate::xyz::Structure;

pub struct CsomStructure {
    pub labels: Vec<String>,
    pub points: Vec<Vector3<f64>>,
}

pub enum CenteringMode { Auto, First, Centroid, Manual }

pub(crate) fn prepare_csom_structure(
    structure: Structure,
    centering_mode: CenteringMode,
    centering_vector: Option<Vector3<f64>>) -> (CsomStructure, f64, Vector3<f64>)
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
        CenteringMode::Manual => center_by_coordinate(&mut points, centering_vector.unwrap()),
    };



    let scaling_factor = normalise(&mut points);

    (CsomStructure {
        labels: striped_labels,
        points
    }, scaling_factor, original_centroid)
}

fn strip_label(label: &str ) -> &str{
    let end = label
        .char_indices() // Iterate over elements of str
        .nth(
            if label.chars().nth(1).is_some_and(|c| c.is_ascii_lowercase()) {2}
            else {1}
            )
        .map_or(label.len(), |(i, _)| i);

    &label[..end]
}




