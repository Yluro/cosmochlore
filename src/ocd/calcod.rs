use crate::xyz::Structure;
use nalgebra::Vector3;
use itertools::Itertools;

#[derive(Debug, Clone)]
pub struct ODResult {
    pub d_mean: f64,
    pub zeta: f64,
    pub delta: f64,
    pub sigma: f64,
    //pub theta: f64,
    //pub vol: f64,
    pub tau: f64,
    pub mu: f64,
}

/// Calculate the octahedral distortion (OD) of a structure.
///
/// The structure must contain a centre atom and exactly six ligand points.
pub fn calculate_od(structure: &Structure) -> Result<ODResult, ODError> {
    const N: i32 = 7;

    // 0. Extract the centre and ligand coordinates from the structure parsed by xyz.rs
    let centre= structure.centre.as_ref().ok_or(ODError::NoCentre)?;

    let centre = Vector3::new(
        centre.coords[0],
        centre.coords[1],
        centre.coords[2],
    );

    let ligands: Vec<Vector3<f64>> = structure.ligands
        .iter()
        .map(|ligand| Vector3::new(ligand.coords[0], ligand.coords[1], ligand.coords[2]))
        .collect();

    if ligands.len() != 6 { return Err(ODError::IncorrectNumberOfPoints {n: ligands.len()}) };

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

    let a = for comb in vectors.iter().combinations(2) {
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


    // 5. Theta (face twisting) calculation ...





    Ok(
        ODResult {
            d_mean,
            zeta,
            delta,
            sigma,
            //theta,
            //vol,
            tau,
            mu,
        }
    )
}

#[derive(Debug)]
pub enum ODError {
    NoCentre,
    IncorrectNumberOfPoints{ n: usize },
}

impl std::fmt::Display for ODError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            ODError::NoCentre => {write!(f, "structure has no central atom")}
            ODError::IncorrectNumberOfPoints{ n} => { write!(f, "wrong number of points, expected: 7, found: {}", n) },
        }
    }
}

impl std::error::Error for ODError {}



mod tests {
    use super::*;
    use crate::xyz::parse_xyz;

    #[test]
    fn matches_octadist_results() {
        let structure = parse_xyz(r".\tests\FeHS.xyz", false).unwrap();

        let calc = calculate_od(&structure);
        assert!(calc.is_ok());

        let calc = calc.unwrap();
        println!("{:?}", calc);
        assert!((calc.d_mean - 2.1623).abs() < 1e-3 );
        assert!((calc.zeta - 0.3621) < 1e-3 );
        assert!((calc.delta - 0.001006).abs() < 1e-3 );
        assert!((calc.sigma - 82.29).abs() < 1e-2 );
        //assert!((calc.theta - 306.79).abs() < 1e-2
        assert!((calc.tau - 54.41).abs() < 1e-2 );
        assert!((calc.mu - 0.17).abs() < 1e-2 );
    }
}