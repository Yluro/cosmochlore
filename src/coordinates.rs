use nalgebra::Vector3;
use crate::xyz::Structure;
use crate::shapes::ReferenceShape;

pub fn points_from_structure(structure: &Structure) -> Vec<Vector3<f64>> {
    let mut points: Vec<Vector3<f64>> = Vec::new();
    if let Some(centre) = &structure.centre {
        points.push(<Vector3<f64>>::from(centre.coords));
    }
    for ligand in &structure.ligands {
        points.push(<Vector3<f64>>::from(ligand.coords));
    }
    points
}

pub fn points_from_reference_shape(shape: &ReferenceShape, use_centre: bool) -> Vec<Vector3<f64>> {
    let mut points: Vec<Vector3<f64>> = Vec::new();
    if use_centre {
        points.push(<Vector3<f64>>::from(shape.centre));
    }
    for ligand in &shape.vertices {
        points.push(<Vector3<f64>>::from(*ligand));
    }
    points
}