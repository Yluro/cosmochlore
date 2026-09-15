/// Functions to convert the input data structures (Structure and ReferenceShape)
/// into Vector3 coordinates for the math modules.
use nalgebra::Vector3;
use crate::xyz::Structure;
use crate::shapes::ReferenceShape;

pub fn points_from_structure(structure: &Structure) -> Vec<Vector3<f64>> {
    structure.atoms().iter().map(|atom| atom.coords).collect()
}

pub fn points_from_reference_shape(shape: &ReferenceShape, use_centre: bool) -> Vec<Vector3<f64>> {
    let mut points: Vec<Vector3<f64>> = Vec::new();
    if use_centre {
        points.push(shape.centre);
    }
    for ligand in &shape.vertices {
        points.push(*ligand);
    }
    points
}