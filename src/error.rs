use crate::csom::types::CsomError;
use crate::odis::OdisError;
use crate::shapes::ShapeLookupError;
use crate::xyz::XyzParseError;
use crate::yaml::YamlParseError;

/// Crate-wide error, covering every module's own error type plus file I/O.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error(transparent)]
    Xyz(#[from] XyzParseError),
    #[error(transparent)]
    Yaml(#[from] YamlParseError),
    #[error(transparent)]
    Shape(#[from] ShapeLookupError),
    #[error(transparent)]
    Csom(#[from] CsomError),
    #[error(transparent)]
    Odis(#[from] OdisError),
    #[error(transparent)]
    Io(#[from] std::io::Error),
}
