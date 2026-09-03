pub mod csom;
mod io;
mod dev;
mod tests;
mod optimize;


#[derive(Debug)]
pub enum CsomError {
    WrongSpaceGroup { pg: String},

}

impl std::fmt::Display for CsomError  {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            CsomError::WrongSpaceGroup { pg } => {write!(f, "wrong point group name: {}", {pg})}
        }
    }
}

impl std::error::Error for CsomError {}
