use crate::data;
use data::standard_shapes::builtin_shapes;

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
    pub centre: [f64; 3],
    /// Coordinates of the vertex positions the reference shape.
    pub vertices: Vec<[f64; 3]>,
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

#[derive(Debug)]
pub enum ShapeLookupError {
    NoShapesForVertexCount(u8),
    IndexOutOfBounds(usize)
}

impl std::fmt::Display for ShapeLookupError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            ShapeLookupError::NoShapesForVertexCount(n) => write!(f, "no reference shapes with {} vertices.", n),
            ShapeLookupError::IndexOutOfBounds(idx) => write!(f, "there is no reference shape with index {}", idx),
        }
    }
}

impl std::error::Error for ShapeLookupError {}

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