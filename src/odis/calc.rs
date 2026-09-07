use crate::odis::{OdisError, OdisResult};
use crate::xyz::Structure;
use itertools::Itertools;
use nalgebra::Vector3;

/// Calculate the octahedral distortion (OD) of a structure.
///
/// The structure must contain a centre atom and exactly six ligand points.
pub fn calculate_od(structure: &Structure) -> Result<OdisResult, OdisError> {
    // const N: i32 = 7;

    // 0. Extract the centre and ligand coordinates from the structure parsed by xyz.rs
    let centre= structure.centre.as_ref().ok_or(OdisError::NoCentre)?;

    let centre = Vector3::new(
        centre.coords[0],
        centre.coords[1],
        centre.coords[2],
    );

    let ligands: Vec<Vector3<f64>> = structure.ligands
        .iter()
        .map(|ligand| Vector3::new(ligand.coords[0], ligand.coords[1], ligand.coords[2]))
        .collect();

    if ligands.len() != 6 { return Err(OdisError::IncorrectNumberOfPoints {n: ligands.len()}) };

    // 1. Build an array of points and vectors of the octahedron for easy calculation later.
    let points = [
        centre,
        ligands[0],
        ligands[1],
        ligands[2],
        ligands[3],
        ligands[4],
        ligands[5],
    ];

    let vectors = [
        points[1] - points[0],
        points[2] - points[0],
        points[3] - points[0],
        points[4] - points[0],
        points[5] - points[0],
        points[6] - points[0],
    ];

    let distances = [
        vectors[0].norm(),
        vectors[1].norm(),
        vectors[2].norm(),
        vectors[3].norm(),
        vectors[4].norm(),
        vectors[5].norm(),
    ];

    // 3. D-mean, zeta (bond length distortion) and Delta (octahedral tilting) and mu (centroid-deviation) calculations.

    let d_mean = distances
        .iter()
        .sum::<f64>() / ligands.len() as f64;

    let zeta = distances
        .iter()
        .map(|di| (di - d_mean).abs())
        .sum::<f64>();

    let delta = distances
        .iter()
        .map(|di| ((di - d_mean)/(d_mean)).powi(2) )
        .sum::<f64>() / distances.len() as f64;

    let mu = vectors.iter().sum::<Vector3<f64>>().norm() / vectors.len() as f64;

    // 4. Cis and trans angle distortion calculation

    let mut angles = Vec::new();

    for comb in vectors.iter().combinations(2) {
        let vi = comb[0];
        let vj = comb[1];
        let angle = (vi.dot(&vj) / (vi.norm() * vj.norm()))
            .clamp(-1.0, 1.0)
            .acos()
            .abs()
            .to_degrees();
        angles.push(angle);
    };
    angles.sort_by(|a, b| a.total_cmp(b));
    let phis = &angles[..12]; // All 12 small angles in the list have to be cis angles.
    let psis = &angles[12..];    // All 3 others have to be trans angles.

    let sigma: f64 = phis.iter().map(|phi| {(90.0 - phi).abs()}).sum();
    let tau: f64 = psis.iter().map(|phi| {(180.0 - phi).abs()}).sum();

    // 5. Theta (face twisting) calculation, following OctaDist's algorithm: walk the eight
    // triangular faces of the (possibly distorted) octahedron and, for each, project the centre
    // and the three non-face ligands onto the face's plane, then sum the deviation of the six
    // resulting angles from the ideal 60°. The final theta is half the total over all 8 faces,
    // since each pair of opposite faces is counted twice (once from either side).
    let ligand_arr: [Vector3<f64>; 6] = [
        ligands[0], ligands[1], ligands[2], ligands[3], ligands[4], ligands[5],
    ];
    let theta = calc_theta(centre, &ligand_arr);

    // 6. Octahedron volume, via decomposition into tetrahedra (see `calc_vol`).
    let vol = calc_vol(centre, &ligand_arr);

    Ok(
        OdisResult {
            d_mean,
            zeta,
            delta,
            sigma,
            theta,
            vol,
            tau,
            mu,
        }
    )
}

/// Compute the angle in degrees between two vectors, in `[0, 180]`.
fn angle_btw_vectors(v1: Vector3<f64>, v2: Vector3<f64>) -> f64 {
    v1.normalize()
        .dot(&v2.normalize())
        .clamp(-1.0, 1.0)
        .acos()
        .to_degrees()
}

/// Compute the angle in degrees between two vectors, signed by whether `v1`, `v2` and `direct`
/// form a right- or left-handed set (i.e. whether going from `v1` to `v2` is CW or CCW when
/// viewed from the `direct` side).
fn angle_sign(v1: Vector3<f64>, v2: Vector3<f64>, direct: Vector3<f64>) -> f64 {
    let v1 = v1.normalize();
    let v2 = v2.normalize();
    let angle = v1.dot(&v2).clamp(-1.0, 1.0).acos().to_degrees();

    let matrix = nalgebra::Matrix3::new(
        v1.x, v1.y, v1.z,
        v2.x, v2.y, v2.z,
        direct.x, direct.y, direct.z,
    );

    if matrix.determinant() < 0.0 { -angle } else { angle }
}

/// Coefficients `(a, b, c, d)` of the plane `a*x + b*y + c*z = d` through three points.
fn find_eq_of_plane(x: Vector3<f64>, y: Vector3<f64>, z: Vector3<f64>) -> (f64, f64, f64, f64) {
    let cross = (z - x).cross(&(y - x));
    let d = cross.dot(&z);
    (cross.x, cross.y, cross.z, d)
}

/// Orthogonal projection of point `p` onto the plane `a*x + b*y + c*z = d`.
fn project_atom_onto_plane(p: Vector3<f64>, a: f64, b: f64, c: f64, d: f64) -> Vector3<f64> {
    let plane = Vector3::new(a, b, c);
    let lambda = (d - plane.dot(&p)) / plane.dot(&plane);
    p + lambda * plane
}

/// The index (0..6) of the ligand whose metal-to-ligand vector is most nearly anti-parallel to
/// (i.e. has the largest angle with) `from`.
fn most_trans_index(from: Vector3<f64>, metal_to_lig: &[Vector3<f64>; 6]) -> usize {
    let mut best_angle = 0.0;
    let mut best_index = 0;
    for (n, v) in metal_to_lig.iter().enumerate() {
        let angle = angle_btw_vectors(from, *v);
        if angle > best_angle {
            best_angle = angle;
            best_index = n;
        }
    }
    best_index
}

/// Reorder the six ligand coordinates into the arrangement `calc_theta` needs in order to walk
/// the octahedron's eight triangular faces, following `OctaDist`'s `determine_faces`. Each of the
/// three rounds below puts the ligand most nearly trans to `ligands[0]`, `[1]` and `[2]` into
/// position `4`, `5` and `3` respectively (`OctaDist` also uses this pass to flag non-octahedral
/// geometries; that diagnostic isn't surfaced here, so we skip straight to the reordering).
fn determine_faces(centre: Vector3<f64>, ligands: &[Vector3<f64>; 6]) -> [Vector3<f64>; 6] {
    let mut lig = *ligands;

    for &(from, to) in &[(0usize, 4usize), (1, 5), (2, 3)] {
        let metal_to_lig: [Vector3<f64>; 6] = std::array::from_fn(|i| lig[i] - centre);
        lig.swap(to, most_trans_index(metal_to_lig[from], &metal_to_lig));
    }

    lig
}

/// Cycle `coord_lig` to the next of the octahedron's eight triangular faces (kept in its first
/// three elements): N1N2N3 -> N1N4N2 -> N1N6N4 -> N1N3N6, then flip to the opposite four faces
/// and repeat. Shared by `calc_theta` and `calc_vol`, which both need to walk all 8 faces.
fn cycle_octahedron_faces(coord_lig: &mut [Vector3<f64>; 6], round: usize) {
    let tmp = coord_lig[1];
    coord_lig[1] = coord_lig[3];
    coord_lig[3] = coord_lig[5];
    coord_lig[5] = coord_lig[2];
    coord_lig[2] = tmp;

    if round == 3 {
        coord_lig.swap(0, 4);
        coord_lig.swap(1, 5);
        coord_lig.swap(2, 3);
    }
}

/// Calculate the Theta (face-twisting) octahedral distortion parameter, in degrees.
fn calc_theta(centre: Vector3<f64>, ligands: &[Vector3<f64>; 6]) -> f64 {
    let mut coord_lig = determine_faces(centre, ligands);
    let mut eight_theta = Vec::with_capacity(8);

    for r in 0..8 {
        let (a, b, c, d) = find_eq_of_plane(coord_lig[0], coord_lig[1], coord_lig[2]);

        let projected_m = project_atom_onto_plane(centre, a, b, c, d);
        let projected_lig4 = project_atom_onto_plane(coord_lig[3], a, b, c, d);
        let projected_lig5 = project_atom_onto_plane(coord_lig[4], a, b, c, d);
        let projected_lig6 = project_atom_onto_plane(coord_lig[5], a, b, c, d);

        let vector_theta = [
            coord_lig[0] - projected_m,
            coord_lig[1] - projected_m,
            coord_lig[2] - projected_m,
            projected_lig4 - projected_m,
            projected_lig5 - projected_m,
            projected_lig6 - projected_m,
        ];

        // The reference face's vertices are traversed either CW or CCW; take the direction from
        // whichever of the first two edges spans the smaller angle.
        let a12 = angle_btw_vectors(vector_theta[0], vector_theta[1]);
        let a13 = angle_btw_vectors(vector_theta[0], vector_theta[2]);
        let direction = if a12 < a13 {
            vector_theta[0].cross(&vector_theta[1])
        } else {
            vector_theta[2].cross(&vector_theta[0])
        };

        let indi_theta = [
            angle_sign(vector_theta[0], vector_theta[3], direction),
            angle_sign(vector_theta[3], vector_theta[1], direction),
            angle_sign(vector_theta[1], vector_theta[4], direction),
            angle_sign(vector_theta[4], vector_theta[2], direction),
            angle_sign(vector_theta[2], vector_theta[5], direction),
            angle_sign(vector_theta[5], vector_theta[0], direction),
        ];

        eight_theta.push(indi_theta.iter().map(|t| (t - 60.0).abs()).sum::<f64>());

        cycle_octahedron_faces(&mut coord_lig, r);
    }

    // Each pair of opposite faces contributes twice (once viewed from either side).
    eight_theta.iter().sum::<f64>() / 2.0
}

/// Six times the signed volume of the tetrahedron with vertices `a`, `b`, `c`, `d`.
fn tetrahedron_volume_x6(a: Vector3<f64>, b: Vector3<f64>, c: Vector3<f64>, d: Vector3<f64>) -> f64 {
    (a - d).dot(&(b - d).cross(&(c - d)))
}

/// Calculate the volume of the (possibly distorted) octahedron spanned by the six ligand points,
/// in cubic Angstrom. Unlike `OctaDist`, which calls out to `scipy`'s `ConvexHull`, this walks
/// the same eight triangular faces used by `calc_theta` and decomposes the octahedron into eight
/// tetrahedra, each with its base on one face and its apex at the ligands' centroid. Since the
/// centroid sits inside the (convex) coordination polyhedron, those eight tetrahedra exactly
/// tile it, so their volumes simply sum to the polyhedron's volume — no averaging or halving
/// needed here, unlike `calc_theta` (this is a straight decomposition, not a doubled count).
fn calc_vol(centre: Vector3<f64>, ligands: &[Vector3<f64>; 6]) -> f64 {
    let centroid = ligands.iter().sum::<Vector3<f64>>() / 6.0;
    let mut coord_lig = determine_faces(centre, ligands);

    let mut volume = 0.0;
    for r in 0..8 {
        volume += tetrahedron_volume_x6(centroid, coord_lig[0], coord_lig[1], coord_lig[2]).abs() / 6.0;
        cycle_octahedron_faces(&mut coord_lig, r);
    }

    volume
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::shapes;
    use crate::xyz::parse_xyz;

    #[test]
    fn matches_octadist_results() {
        let structure = parse_xyz(r".\tests\FeHS.xyz", Some(0)).unwrap();

        let calc = calculate_od(&structure);
        assert!(calc.is_ok());

        let calc = calc.unwrap();
        println!("{:?}", calc);
        assert!((calc.d_mean - 2.1623).abs() < 1e-3 );
        assert!((calc.zeta - 0.3621) < 1e-3 );
        assert!((calc.delta - 0.001006).abs() < 1e-3 );
        assert!((calc.sigma - 82.29).abs() < 1e-2 );
        assert!((calc.theta - 306.81).abs() < 1e-2 );
        // Cross-checked against `scipy.spatial.ConvexHull` on the same six ligand points.
        assert!((calc.vol - 12.9644).abs() < 1e-3 );
        assert!((calc.tau - 54.41).abs() < 1e-2 );
        assert!((calc.mu - 0.17).abs() < 1e-2 );
    }

    #[test]
    fn calc_vol_matches_known_polyhedron_volumes() {
        // Regular octahedron with vertices at unit distance along the axes: two square
        // pyramids of base area 2 and height 1, so 2 * (1/3 * 2 * 1) = 4/3.
        let octahedron = shapes::structure_from_shape(6, 2);
        assert!((calculate_od(&octahedron).unwrap().vol - 4.0 / 3.0).abs() < 1e-6);

        // Ideal (uniform) trigonal prism, cross-checked against `scipy.spatial.ConvexHull`.
        let prism = shapes::structure_from_shape(6, 3);
        assert!((calculate_od(&prism).unwrap().vol - 2.25).abs() < 1e-3);
    }

    /// Rotate the "top" three ligands (the ones on the +z side of the centre) of a trigonal
    /// prism around the z axis, tracing out the Bailar twist pathway towards an octahedron.
    fn twist_top_face(structure: &Structure, degrees: f64) -> Structure {
        let rotation = nalgebra::Rotation3::from_axis_angle(&Vector3::z_axis(), degrees.to_radians());
        let mut twisted = structure.clone();
        for atom in twisted.ligands.iter_mut() {
            if atom.coords.z > 0.0 {
                atom.coords = rotation * atom.coords;
            }
        }
        twisted
    }

    #[test]
    fn ideal_trigonal_prism_gives_the_maximum_theta() {
        let octahedron = shapes::structure_from_shape(6,2);
        let prism = shapes::structure_from_shape(6,3);

        // Theta is 0 for the perfectly staggered (Oh) octahedron...
        let octahedron_theta = calculate_od(&octahedron).unwrap().theta;
        assert!(octahedron_theta.abs() < 1e-6);

        // ...and reaches its maximum, 1170°, at the opposite (eclipsed, D3h) extreme: the ideal
        // trigonal prism.
        let prism_theta = calculate_od(&prism).unwrap().theta;
        assert!((prism_theta - 1170.0).abs() < 1e-3);
    }
}
