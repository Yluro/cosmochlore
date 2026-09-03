pub mod calc;

pub use calc::calculate_od;

#[derive(Debug, Clone)]
pub struct OdisResult {
    pub d_mean: f64,
    pub zeta: f64,
    pub delta: f64,
    pub sigma: f64,
    //pub theta: f64,
    //pub vol: f64,
    pub tau: f64,
    pub mu: f64,
}

#[derive(Debug)]
pub enum OdisError {
    NoCentre,
    IncorrectNumberOfPoints{ n: usize },
}

impl std::fmt::Display for OdisError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            OdisError::NoCentre => {write!(f, "structure has no central atom")}
            OdisError::IncorrectNumberOfPoints{ n} => { write!(f, "wrong number of points, expected: 7, found: {}", n) },
        }
    }
}

impl std::error::Error for OdisError {}