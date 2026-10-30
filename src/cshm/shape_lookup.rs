//! Resolves the reference shapes a `cshm` run compares against: the built-in ones by vertex
//! count and index, plus the check every user-defined shape must pass.

use crate::data::standard_shapes::{BuiltinShape, ReferenceShape, BUILTIN_SHAPES};

pub(crate) fn shape_by_vertex(no_vertices: u8) -> Result<&'static [BuiltinShape], ShapeLookupError> {

    let (_, shapes) = BUILTIN_SHAPES
        .iter()
        .find(|(n, _)| *n == no_vertices)
        .ok_or(ShapeLookupError::NoShapesForVertexCount(no_vertices))?;

    Ok(shapes)

}

pub(crate) fn shapes_by_index(
    shapes: &[BuiltinShape],
    indices: &[usize],
) -> Result<Vec<ReferenceShape>, ShapeLookupError> {
    let mut result: Vec<ReferenceShape> = Vec::new();

    for idx in indices {
        let found: Option<&BuiltinShape> = shapes.get(*idx);

        match found {
            Some(shape) => result.push(ReferenceShape::from(shape)),
            None => return Err(ShapeLookupError::IndexOutOfBounds(*idx)),
        }
    }
    Ok(result)
}

pub fn resolve_shapes(
    no_vertices: u8,
    indices: Option<&[usize]>,
) -> Result<Vec<ReferenceShape>, ShapeLookupError> {
    let ref_shapes = shape_by_vertex(no_vertices)?;

    match indices {
        Some(indices) => shapes_by_index(ref_shapes, indices),
        None => Ok(ref_shapes.iter().map(ReferenceShape::from).collect()), // No indices were supplied, so all shapes are passed.
    }
}

#[derive(Debug, thiserror::Error)]
pub enum ShapeLookupError {
    #[error("no reference shapes with {0} vertices.")]
    NoShapesForVertexCount(u8),
    #[error("there is no reference shape with index {0}")]
    IndexOutOfBounds(usize),
    #[error("wrong number of vertices for shape {symbol}: expected {expected}, found {found}")]
    VertexCountMismatch {
        symbol: String,
        expected: u8,
        found: usize,
    }, // Used for user input.
}

pub fn check_vertex_count(shape: &ReferenceShape, count: u8) -> Result<(), ShapeLookupError> {
    if shape.vertices.len() == count as usize {
        return Ok(());
    }
    Err(ShapeLookupError::VertexCountMismatch {
        symbol: shape.symbol.clone(),
        expected: count as usize as u8,
        found: shape.vertices.len(),
    })
}
