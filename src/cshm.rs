use nalgebra::{Matrix3, OVector, Vector3};

pub fn center_and_normalise(points: &mut [[f64; 3]]) {
    let n = points.len() as f64;
    let mut centroid = [0.0; 3];

    for point in points.iter() {
        centroid[0] += point[0];
        centroid[1] += point[1];
        centroid[2] += point[2];
    }

    centroid[0] /= n;
    centroid[1] /= n;
    centroid[2] /= n;

    for point in points.iter_mut() {
        point[0] -= centroid[0];
        point[1] -= centroid[1];
        point[2] -= centroid[2];
    }

    let mut s2 = 0.0;
    for point in points.iter_mut() {
        s2 += point[0] * point[0] + point[1] * point[1] + point[2] * point[2]
    }

    // Normalization
    // Σ A²·|P_i|² = N
    // A = sqrt(N / Σ|P_i|²)

    let scale_factor = (n / s2).sqrt();

    for point in points.iter_mut() {
        point[0] *= scale_factor;
        point[1] *= scale_factor;
        point[2] *= scale_factor;
    }
}

pub fn correlation_matrix(reference: &[[f64; 3]], problem: &[[f64; 3]]) -> Matrix3<f64> {
    /// Expects centered and normalised points.

    let mut h = Matrix3::<f64>::zeros();

    for (p, q) in reference.iter().zip(problem.iter()) {
        let p_vec = Vector3::new(p[0], p[1], p[2]);
        let q_vec = Vector3::new(q[0], q[1], q[2]);

        h += p_vec * q_vec.transpose();
    }
    h
}

pub fn optimal_rotation(h: Matrix3<f64>) -> (Matrix3<f64>, Vector3<f64>)  {
    let svd = h.svd(true, true); // computes U and V^T
    let u = svd.u.unwrap();
    let mut v_t = svd.v_t.unwrap();
    let s = svd.singular_values;

    // The SVD algorithm can produce the reflected shape instead of the actual one.
    // To convert back to the desired rotation, flip one of the rows of the v_t matrix.
    if v_t.transpose().determinant() * u.transpose().determinant() < 0.0 {
        // TODO Flip one the rows so Det becomes positive.
        v_t.set_row(2, &(-v_t.row(2)));
    };
    (v_t.transpose() * u.transpose(), s)
}

pub fn shape_measure(singular_values: &Vector3<f64>, n: usize) -> f64 {
    let a: f64 = singular_values.iter().sum();
    (1.0 - a*a/(n as f64 * n as f64)) * 100.0
}


mod tests {
    use crate::shapes::ReferenceShape;
    use super::*;


    #[test]
    fn centre_and_normalises_correctly() {
        let mut points = Vec::from(
            [[1.0, 0.0, 0.0], // Regular octahedron centered at 1 0 0
                [2.0, 0.0, 1.0],
                [1.0, 1.0, 0.0],
                [1.0, 0.0, 1.0],
                [-0.0, 0.0, 1.0],
                [1.0, -1.0, 0.0],
                [1.0, 0.0, -1.0],
            ]
        );

        center_and_normalise(&mut points);

        let mut centroid = [0.0; 3];

        for point in points.iter() {
            centroid[0] += point[0];
            centroid[1] += point[1];
            centroid[2] += point[2];
        }

        centroid[0] /= points.len() as f64;
        centroid[1] /= points.len() as f64;
        centroid[2] /= points.len() as f64;

        for point in centroid.iter() {assert!(point.abs() < 1e-10)}

        let mut s2: f64 = 0.0;
        for point in points.iter() {s2 += point[0] * point[0] + point[1] * point[1] + point[2] * point[2]};
        assert!((s2 - points.len() as f64).abs() < 1e-10);

    }


    #[test]
    fn recovers_from_rotation() {
        let reference = [
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 1.0],
            [1.0, 1.0, 0.0],
            [1.0, 0.0, 1.0],
            [-0.0, 0.0, 1.0],
        ];

        let know_rotataion = Matrix3::new(
            0.0, -1.0, 0.0,
            1.0,  0.0, 0.0,
            0.0,  0.0, 1.0,
        );

        let rotated: Vec<[f64; 3]> = reference.iter().map(|p| {
            let v = Vector3::new(p[0], p[1], p[2]);
            let r = know_rotataion * v;
            [r[0], r[1], r[2]]
        }).collect();

        let h = correlation_matrix(&reference, &rotated);
        let (recovered_r, _) = optimal_rotation(h);


        recovered_r.iter().zip(know_rotataion.iter()).for_each(|(r1, r2)| {assert!((r1-r2).abs() < 1e-10)});

    }

    #[test]
    fn perfect_gives_zero() {
        let mut points = [
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];

        center_and_normalise(&mut points);

        let h = correlation_matrix(&points, &points);
        let (_, s) = optimal_rotation(h);
        let s = shape_measure(&s, points.len());
        println!("s is: {}", s);
        assert!(s.abs() < 1e-10);
    }

    #[test]
    fn water_matches_result_from_shape21 () {
        // Coords from problem taken from water dataset from Olex2
        let mut problem = [
            [0.0,       8.0648,     0.0],       // Central Mn
            [-1.332417, 6.56007,   -1.099508],  // N2'
            [1.3324,    9.5695,     1.0995],    // N2
            [0.521462,  6.437162,   1.333904],  // O4
            [-0.5215,   9.6924,    -1.3339],    // O4'
            [1.619207,  7.429778,  -1.386425],  // O5
            [-1.6192,   8.6998,     1.3864],    // O5'
        ];


        // This list of atoms is manually set in order to correspond
        // to the correct assignation of point pairs.
        let mut reference = [
            [0.0, 0.0, 0.0], // Central atom first.
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ];

        center_and_normalise(&mut reference);
        center_and_normalise(&mut problem);

        let h = correlation_matrix(&reference, &problem);
        let (_, s) = optimal_rotation(h);
        let s = shape_measure(&s, problem.len());
        println!("s is: {}", s);
        assert!((s - 0.18).abs() < 1e-2);
    }
}