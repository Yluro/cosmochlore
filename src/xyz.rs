use std::string::ParseError;

#[derive(Debug, Clone)]
pub struct Atom {
    pub label: String,
    pub coords:  [f64; 3],
}

#[derive(Debug, Clone)]
pub struct Structure {
    pub centre: Option<Atom>,
    pub ligands: Vec<Atom>,
}

pub fn parse_xyz(path: &str) -> Result<Structure, XyzParseError> {

    unimplemented!()
}

#[derive(Debug)]
pub enum XyzParseError {
    Io(std::io::Error), // Can't read file
    Empty,              // File is empty
    BadHeader(String),  // File has wrong header
    BadLine { line_no: usize, line: String },  // File has a bad line (not enough whitespaces)
    BadCoordinate { line_no: usize, source: std::num::ParseFloatError }, // Coordinates are not floats.
    TooFewAtoms(usize), // File contains two atoms or fewer.
}

impl std::fmt::Display for XyzParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            XyzParseError::Io(err) => write!(f, "could not read file: {}", err),
            XyzParseError::Empty => write!(f, "empty file"),
            XyzParseError::BadHeader(err) => write!(f, "bad header: {}", err),
            XyzParseError::BadLine { line_no, line } => write!(f, "bad line in {}: {}", line_no, line),
            XyzParseError::BadCoordinate { line_no, source } => write!(f, "bad coordinates in line {}: {}", line_no, source),
            XyzParseError::TooFewAtoms(count) => write!(f, "too few atoms: expected at least 3, found {}", count),
        }
    }
}

impl std::error::Error for XyzParseError {}

impl From<std::io::Error> for XyzParseError {
    fn from(err: std::io::Error) -> Self {
        XyzParseError::Io(err)
    }
}