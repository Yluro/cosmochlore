pub mod odis;

pub use odis::calculate_od;

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