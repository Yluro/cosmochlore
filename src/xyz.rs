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

fn read_file(path: &str) -> Result<String,XyzParseError> {
    std::fs::read_to_string(path).map_err(XyzParseError::Io)
}

pub fn parse_xyz(path: &str, not_centered: bool) -> Result<Structure, XyzParseError> {
    let content = read_file(path)?;
    parse_xyz_contents(&content, not_centered)
}

fn parse_xyz_contents(content: &str, not_centered: bool) -> Result<Structure, XyzParseError> {
    // Parses through a .xyz file and gets the header, the comment and the atom table.
    let mut lines = content.lines();

    // Get the header of the line
    let header = match lines.next() {
        Some(line) => line,
        None => return Err(XyzParseError::Empty),
    };

    // get the number of atoms (first line of .xyz)
    let no_atoms = match header.parse::<usize>() {
        Ok(n) => n,
        Err(_) => return Err(XyzParseError::BadHeader(header.to_string()))
    };

    // Comment line in xyz, useless
    let _comment = lines.next();

    //
    let mut atoms: Vec<Atom> = Vec::new();
    for (i, line) in lines.enumerate() {
        if line.trim().is_empty() { continue; } // Skip blank lines

        // Split the line into a Vector of 4 strings.
        let parts = line.split_whitespace().collect::<Vec<&str>>();

        // If the split has more or less than four parts (label, x, y, z) return error.
        if parts.len() != 4 { return Err(XyzParseError::BadLine {line_no: i + 3, line: line.to_string()})}

        // Main logic, asign labels and xyz coords.
        let label = parts[0];
        let mut coords = [0.0; 3]; // Create array [0.0, 0.0, 0.0]
        for (j, part) in parts[1..4].iter().enumerate() {
            coords[j] = match part.parse::<f64>() {
                Ok(n) => n,
                Err(source) => return Err(XyzParseError::BadCoordinate {line_no: i + 3, source}),
            };
        }

        let atom = Atom {label:label.to_string(), coords};
        atoms.push(atom);
    };

    // Sanity checks.
    // If the number of atoms doesn't coincide with the header, return Err
    if atoms.len() != no_atoms {
        return Err(XyzParseError::IncorrectAtomCount{header: no_atoms, count: atoms.len()})
    }
    // Or if there are less than three atoms, return Err
    if atoms.len() < 3 {
        return Err(XyzParseError::TooFewAtoms(atoms.len()))
    }


    let structure = if !not_centered {
        let center_atom = atoms.remove(0);
        Structure {
            centre: Some(center_atom),
            ligands: atoms,
            }
    } else {
        Structure {
            centre: None,
            ligands: atoms,
        }
    };
    Ok(structure)
}

#[derive(Debug)]
pub enum XyzParseError {
    Io(std::io::Error), // Can't read file
    Empty,              // File is empty
    BadHeader(String),  // File has wrong header
    BadLine { line_no: usize, line: String },  // File has a bad line (not enough whitespaces)
    BadCoordinate { line_no: usize, source: std::num::ParseFloatError }, // Coordinates are not floats.
    TooFewAtoms(usize), // File contains two atoms or fewer.
    IncorrectAtomCount{ header: usize, count: usize}, // The xyz table has a different number of atoms than expected in the header.
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
            XyzParseError::IncorrectAtomCount { header, count} => write!(f, "incorrect number of atoms, expected {}, found {}", header, count),
        }
    }
}

impl std::error::Error for XyzParseError {}

impl From<std::io::Error> for XyzParseError {
    fn from(err: std::io::Error) -> Self {
        XyzParseError::Io(err)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn parses_valid_file() {
        let input = "5\ncomment\nC 0.0 0.0 0.0\nH1 0.5 0.1 0.9\nH2 0.2 0.8 -0.6\nH3 0.3 -0.9 -0.4\nH4 0.3 -0.5 1.2 ";
        let result = parse_xyz_contents(input, false);

        println!("{:?}", result);
        assert!(result.is_ok());
        let result = result.unwrap();
        assert!(result.centre.is_some());
        assert_eq!(result.ligands.len(), 4);
        assert_eq!(result.ligands[0].label, "H1");
        assert_eq!(result.ligands[0].coords[0], 0.5);
    }

    #[test]
    fn empty_file_error() {
        let contents = "";
        let result = parse_xyz_contents(contents, false);
        assert!(matches!(result, Err(XyzParseError::Empty)));
    }

    #[test]
    fn bad_header_error() {
        let input = "banana\ncomment\nC 0.0 0.0 0.0\nH1 0.5 0.1 0.9\nH2 0.2 0.8 -0.6\nH3 0.3 -0.9 -0.4\nH4 0.3 -0.5 1.2 ";
        let result = parse_xyz_contents(input, true);
        assert!(matches!(result, Err(XyzParseError::BadHeader(_))));
    }

    #[test]
    fn bad_line_error() {
        let input = "5\ncomment\nC 0.0 0.0\nH1 0.5 0.1 0.9\nH2 0.2 0.8 -0.6\nH3 0.3 -0.9 -0.4\nH4 0.3 -0.5 1.2 ";
        let result = parse_xyz_contents(input, true);
        assert!(matches!(result, Err(XyzParseError::BadLine {line_no: 3, line: _})));
    }

    #[test]
    fn bad_coordinate_error() {
        let input = "5\ncomment\nC banana 0.0 0.0\nH1 0.5 0.1 0.9\nH2 0.2 0.8 -0.6\nH3 0.3 -0.9 -0.4\nH4 0.3 -0.5 1.2 ";
        let result = parse_xyz_contents(input, true);
        assert!(matches!(result, Err(XyzParseError::BadCoordinate {line_no: 3, source: _})));
    }

    #[test]
    fn too_few_error() {
        let input = "2\ncomment\nC 0.0 0.0 0.0\nH1 0.5 0.1 0.9\n";
        let result = parse_xyz_contents(input, true);
        assert!(matches!(result, Err(XyzParseError::TooFewAtoms(2))));
    }

    #[test]
    fn invalid_count_error() {
        let input = "2\ncomment\nC 0.0 0.0 0.0\nH1 0.5 0.1 0.9\nH2 0.2 0.8 -0.6\nH3 0.3 -0.9 -0.4\nH4 0.3 -0.5 1.2 ";
        let result = parse_xyz_contents(input, true);
        assert!(matches!(result, Err(XyzParseError::IncorrectAtomCount {header: 2, count: 5})))
    }
}
