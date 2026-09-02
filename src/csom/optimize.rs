use nalgebra::Vector3;

/// Generates N evenly-spaced points around a sphere of Radius = 1.
fn fibonacci_sphere_sampling(n: usize) -> Vec<Vector3<f64>> {

    if n == 0 { return vec![Vector3::new(0.0, 0.0, 1.0)]; }

    let phi: f64 = std::f64::consts::PI * (5.0f64.sqrt() - 1f64);
    let mut points = Vec::with_capacity(n);

    for i in 0..n {

        let y = 1.0 - (i as f64 / (n as f64 - 1.0));
        let r = (1.0 - y*y).sqrt();
        let theta = phi * i as f64;

        points[i] = Vector3::new(r * theta.cos(), y, r * theta.sin());

    };
    points
}