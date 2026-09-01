use nalgebra::Vector3;

/// Places the centroid of the points at the origin.
pub(crate) fn center_by_centroid(points: &mut [Vector3<f64>]) -> Vector3<f64> {

    let n = points.len() as f64;
    let centroid = points.iter().sum::<Vector3<f64>>() / n;

    for p in points.iter_mut() {
        *p -= centroid;
    }

    centroid
}

/// Places the first point in the list at the origin.
pub(crate) fn center_by_first_point(points: &mut [Vector3<f64>]) -> Vector3<f64> {

    let n = points.len() as f64;
    let p0 = points[0];

    for p in points.iter_mut() {
        *p -= p0;
    }

    p0
}

pub(crate) fn center_by_coordinate(points: &mut [Vector3<f64>], centre: Vector3<f64>) -> Vector3<f64> {
    for p in points.iter_mut() {
        *p -= centre;
    }
    centre
}


/// Assumes centered points.
pub(crate) fn normalise(points: &mut [Vector3<f64>]) -> f64 {

    let n = points.len() as f64;
    let sq = points.iter().map(|v| v.norm_squared()).sum::<f64>();

    let scale_factor = (n / sq).sqrt();

    for p in points.iter_mut() {
        *p *= scale_factor;
    }

    scale_factor
}