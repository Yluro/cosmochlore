use crate::data;
use crate::xyz::{Atom, Structure};
use data::standard_shapes::builtin_shapes;
use nalgebra::Vector3;

#[derive(Debug, Clone, PartialEq)]
pub struct ReferenceShape {
    /// Symbol of the reference shape. E.g. TD-4
    pub symbol: String,
    /// Name of the reference shape. E.g. Tetrahedron
    pub name: String,
    /// ID of the reference shape. E.g. 1
    pub id: u8,
    /// Point group symmetry of the shape. E.g. C2v
    pub symm: String,
    /// Coordinates of the centre of the reference shape. Usually [0.0, 0.0, 0.0]
    pub centre: Vector3<f64>,
    /// Coordinates of the vertex positions the reference shape.
    pub vertices: Vec<Vector3<f64>>,
}

impl ReferenceShape {
    /// Centre and vertices as one ordered list, centre first when `use_centre` is true.
    pub fn points(&self, use_centre: bool) -> Vec<Vector3<f64>> {
        let mut points = Vec::with_capacity(self.vertices.len() + 1);
        if use_centre {
            points.push(self.centre);
        }
        points.extend(self.vertices.iter().copied());
        points
    }
}


/// Build a structure from a given reference shape vertices and index.
pub fn structure_from_shape(vertices: u8, index: usize) -> Structure {
    let shape = resolve_shapes(vertices, Some(&[index])).unwrap().remove(0);
    Structure {
        centre: Some(Atom { label: "M".to_string(), coords: shape.centre }),
        ligands: shape.vertices.iter()
            .map(|v| Atom { label: "L".to_string(), coords: *v })
            .collect(),
    }
}


fn shape_by_vertex(no_vertices: u8) -> Result<Vec<ReferenceShape>, ShapeLookupError> {
    let shapes_map = builtin_shapes();
    match shapes_map.get(&no_vertices) {
        Some(shapes) => Ok(shapes.clone()),
        None => Err(ShapeLookupError::NoShapesForVertexCount(no_vertices)),
    }
}

fn shapes_by_index(shapes: &[ReferenceShape], indices: &[usize]) -> Result<Vec<ReferenceShape>, ShapeLookupError> {
    let mut result: Vec<ReferenceShape> = Vec::new();

    for idx in indices {
        let found: Option<&ReferenceShape> = shapes.get(*idx);

        match found {
            Some(shape) => result.push(shape.clone()),
            None => return Err(ShapeLookupError::IndexOutOfBounds(*idx))
        }

    }
    Ok(result)
}

pub fn resolve_shapes(no_vertices: u8, indices: Option<&[usize]>) -> Result<Vec<ReferenceShape>, ShapeLookupError> {
    let ref_shapes = shape_by_vertex(no_vertices);

    match indices {
        Some(indices) => shapes_by_index(ref_shapes?.as_slice(), indices),
        None => Ok(ref_shapes?.to_vec()),
    }
}

#[derive(Debug, thiserror::Error)]
pub enum ShapeLookupError {
    #[error("no reference shapes with {0} vertices.")]
    NoShapesForVertexCount(u8),
    #[error("there is no reference shape with index {0}")]
    IndexOutOfBounds(usize),
    #[error("wrong number of vertices for shape {symbol}: expected {expected}, found {found}")]
    VertexCountMismatch { symbol: String, expected: u8, found: usize}, // Used for user input.
}

pub fn check_vertex_count(shape: &ReferenceShape, count: u8) -> Result<(), ShapeLookupError> {
    if shape.vertices.len() == count as usize {
        return Ok(());
    }
    Err(ShapeLookupError::VertexCountMismatch {
        symbol: shape.symbol.clone(),
        expected: count as usize as u8,
        found: shape.vertices.len()
    })
}


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn finds_correct_number_of_shapes() {
        let map = builtin_shapes();

        let mut counts: Vec<usize> = map.values().map(|shapes| shapes.len()).collect();
        counts.sort();

        let mut expected = vec![3, 4, 4, 5, 5, 7, 13, 13, 13, 7, 13, 1, 2, 1, 1];
        expected.sort();

        assert_eq!(counts, expected)
    }


    #[test]
    fn matches_shape21_library() {
        let shapes = shape_by_vertex(12).unwrap();
        let result = shapes_by_index(&shapes, &[12]).unwrap();

        assert_eq!(result[0].name, "Sphenomegacorona J88");
        assert_eq!(result[0].symbol, "JSPMC-12");

        let shapes = shape_by_vertex(6).unwrap();
        let result = shapes_by_index(&shapes, &[0, 2]).unwrap();

        assert_eq!(result[0].name, "Hexagon");
        assert_eq!(result[1].symbol, "OC-6");
    }

    #[test]
    fn no_indices_pass() {
        let expect = 13;
        let result = resolve_shapes(12, None).unwrap().len();
        assert_eq!(result, expect);
    }

    #[test]
    fn error_on_bad_vertex_count() {
        let result = resolve_shapes(99, Some([0, 1].as_slice()));
        assert!(matches!(result, Err(ShapeLookupError::NoShapesForVertexCount(99))));
    }

    #[test]
    fn error_on_bad_index() {
        let result = resolve_shapes(6, Some([0, 1, 7].as_slice()));
        assert!(matches!(result, Err(ShapeLookupError::IndexOutOfBounds(7))));
    }
}