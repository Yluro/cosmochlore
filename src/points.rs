use crate::xyz::Structure;
use crate::shapes::ReferenceShape;

pub fn points_from_structure(structure: &Structure) -> Vec<[f64; 3]> {
    let mut points: Vec<[f64; 3]> = Vec::new();
    if let Some(centre) = &structure.centre {
        points.push(centre.coords);
    }
    for ligand in &structure.ligands {
        points.push(ligand.coords);
    }
    points
}

pub fn points_from_reference_shape(shape: &ReferenceShape, use_centre: bool) -> Vec<[f64; 3]> {
    let mut points: Vec<[f64; 3]> = Vec::new();
    if use_centre {
        points.push(shape.centre);
    }
    for ligand in &shape.vertices {
        points.push(*ligand);
    }
    points
}