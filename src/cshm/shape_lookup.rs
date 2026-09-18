//! Resolves the reference shapes a `cshm` run compares against: the built-in ones by vertex
//! count and index, plus the check every user-defined shape must pass.

use crate::data::standard_shapes::{ReferenceShape, builtin_shapes};

pub(crate) fn shape_by_vertex(no_vertices: u8) -> Result<Vec<ReferenceShape>, ShapeLookupError> {
    let shapes_map = builtin_shapes();
    match shapes_map.get(&no_vertices) {
        Some(shapes) => Ok(shapes.clone()),
        None => Err(ShapeLookupError::NoShapesForVertexCount(no_vertices)),
    }
}

pub(crate) fn shapes_by_index(
    shapes: &[ReferenceShape],
    indices: &[usize],
) -> Result<Vec<ReferenceShape>, ShapeLookupError> {
    let mut result: Vec<ReferenceShape> = Vec::new();

    for idx in indices {
        let found: Option<&ReferenceShape> = shapes.get(*idx);

        match found {
            Some(shape) => result.push(shape.clone()),
            None => return Err(ShapeLookupError::IndexOutOfBounds(*idx)),
        }
    }
    Ok(result)
}

pub fn resolve_shapes(
    no_vertices: u8,
    indices: Option<&[usize]>,
) -> Result<Vec<ReferenceShape>, ShapeLookupError> {
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
