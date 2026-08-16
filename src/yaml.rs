use crate::shapes::ReferenceShape;

#[derive(Debug, PartialEq)]
enum Section {
    Unset,
    Vertices,
    Centre,
}


fn read_yaml(path: &str) -> Result<String, YamlParseError> {
    Ok(std::fs::read_to_string(path)?)
}

fn parse_yaml_str(content: &str, file: &str) -> Result<Vec<ReferenceShape>, YamlParseError> {

    if content.trim().is_empty() {
        return Err(YamlParseError::FileEmpty { file: file.to_string() })
    }

    let mut shapes: Vec<ReferenceShape> = Vec::new();

    let mut symbol: Option<String> = None;
    let mut id: Option<u8> = None;
    let mut symm: Option<String> = None;
    let mut name: Option<String> = None;
    let mut vertices: Vec<[f64;3]> = Vec::new();
    let mut centre: Vec<[f64;3]> = Vec::new();
    let mut current_section = Section::Unset;

    for (line_no, raw_line) in content.lines().enumerate() {

        if raw_line.trim().is_empty() { continue; }

        // New entry means not start with whitespace and end with ":"
        let is_new_entry = !raw_line.starts_with(' ') && raw_line.trim_end().ends_with(':');
        let line = raw_line.trim();

        // If it's a new entry, finalise the previous entry and get the symbol of the next one.
        if is_new_entry {
            if symbol.is_some() {
                let shape = finalise_entry(&file, &symbol, &id, &symm, &name, &vertices, &centre);
                shapes.push(shape?);
            }

            let new_symbol = line.trim_end_matches(':').to_string();
            if new_symbol.is_empty() {
                return Err(YamlParseError::BadSymbol { file: file.to_string(), line_no });
            }

            // reset everything for the new entry
            symbol = Some(line.trim_end_matches(':').to_string());
            id = None;
            symm = None;
            name = None;
            vertices.clear();
            centre.clear();
            current_section = Section::Unset;
            continue;
        }

        // This block checks for field keywords like: id, symm, name.
        // If it finds the keyword but doesn't finde value, returns YamlParseError::EmptyFieldValue
        if line.starts_with("id") {
            let value = line.split(':').nth(1).unwrap().trim();

            if value.is_empty() {
                return Err(YamlParseError::MissingField {
                    file: file.to_string(), symbol: symbol.clone().unwrap_or_default(), field: "id"
                })
            }

            id = Some(value.parse::<u8>().unwrap());
            continue;
        }

        if line.starts_with("symmetry") ||  line.starts_with("symm"){
            let value = line.split(':').nth(1).unwrap_or("").trim();
            if value.is_empty() {
                return Err(YamlParseError::NoValue {
                    file: file.to_string(), symbol: symbol.clone().unwrap_or_default(), field: "symm", line_no,
                });
            }
            symm = Some(value.to_string());
            continue;
        }

        if line.starts_with("name") {
            let value = line.split(':').nth(1).unwrap_or("").trim();
            if value.is_empty() {
                return Err(YamlParseError::NoValue {
                    file: file.to_string(), symbol: symbol.clone().unwrap_or_default(), field: "name", line_no,
                });
            }
            name = Some(value.to_string());
            continue;
        }

        if line.starts_with("centre") ||  line.starts_with("center") ||  line.starts_with("metal"){
            current_section = Section::Centre;
            continue;
        }

        if line.starts_with("vertices") || line.starts_with("ligands") {
            current_section = Section::Vertices;
            continue;
        }

        // This block parses coordinate lines.
        if line.starts_with("-") {
            let inner = line
                .trim_start_matches('-')
                .trim()
                .trim_end_matches(']')
                .trim_start_matches('[');

            let parts: Vec<&str> = inner.split(',').collect();
            if parts.len() != 3 {
                return Err(YamlParseError::BadCoordinate {file: file.to_string(), line_no})
            }


            let mut numbers: Vec<f64> = Vec::new();
            for part in parts {
                let x = part.trim().parse::<f64>()
                    .map_err(|_| YamlParseError::BadCoordinate {file: file.to_string(), line_no})?;
                numbers.push(x);
            }
            let xyz = [numbers[0], numbers[1], numbers[2]];

            match current_section {
                Section::Vertices => vertices.push(xyz),
                Section::Centre => centre.push(xyz),
                Section::Unset => return Err(YamlParseError::UnexpectedCoordinate {file: file.to_string(), line_no}),

            }
            continue;
        }
        // If a line is not caught by any if block, means it contains something unknown to the parser
        return Err(YamlParseError::UnparsableLine {file: file.to_string(), line_no })
    }

    // finalise the LAST entry (no more "new entry" line to trigger it)
    if symbol.is_some() {
        let shape = finalise_entry(file, &symbol, &id, &symm, &name, &vertices, &centre)?;
        shapes.push(shape);
    }

    Ok(shapes)
}

fn finalise_entry(
    file: &str,
    symbol: &Option<String>,
    id: &Option<u8>,
    symm: &Option<String>,
    name: &Option<String>,
    vertices: &[[f64; 3]],
    centre: &[[f64; 3]]
) -> Result<ReferenceShape, YamlParseError> {

    let symbol = symbol.clone().ok_or_else(|| YamlParseError::MissingField {
        file: file.to_string(), symbol: "?".to_string(), field: "symbol"
    })?;

    let symm = symm.clone().ok_or_else(|| YamlParseError::MissingField {
        file: file.to_string(), symbol: symbol.to_string(), field: "symm"
    })?;

    let name = name.clone().ok_or_else(|| YamlParseError::MissingField {
        file: file.to_string(), symbol: symbol.to_string(), field: "name"
    })?;

    if vertices.is_empty() {
        return Err(YamlParseError::MissingField {
            file: file.to_string(), symbol: symbol.to_string(), field: "vertices"
        })
    }

    if centre.is_empty() {
        return Err(YamlParseError::MissingField {
            file: file.to_string(), symbol: symbol.to_string(), field: "centre"
        })
    }


    Ok(ReferenceShape {
        symbol,
        name,
        id: id.unwrap_or(0),
        symm,
        centre: centre[0],
        vertices: vertices.to_vec(),
    })

}


pub fn parse_custom_shapes(path: &str) -> Result<Vec<ReferenceShape>, YamlParseError> {
    let content = read_yaml(path)?;
    parse_yaml_str(&content, &path)
}





// ERROR HANDLING
#[derive(Debug)]
pub enum YamlParseError {
    Io(std::io::Error),                                     // Can't read file
    FileEmpty { file: String },                              // File is empty
    BadSymbol { file: String, line_no: usize },                  // Entry appears but is in wrong format for example ":" meets criteria of starting with no ' ' and ending with ':'
    MissingField {file: String, symbol: String, field: &'static str }, // Some entry has missing fields
    NoValue {file: String, symbol: String, field: &'static str, line_no: usize },   // field appeared, but has no value
    UnparsableLine {file: String, line_no: usize },         // Some entry has unexpected fields. It reaches the end of the "if chain" without getting parsed.
    UnexpectedCoordinate {file: String, line_no: usize },   // Coordinate appears when Section::Unset
    BadCoordinate {file: String, line_no: usize },          // wrong coordinate format
}

impl std::fmt::Display for YamlParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            YamlParseError::Io(err) => write!(f, "could not read file: {err}"),
            YamlParseError::FileEmpty {file} => write!(f, "could not read {file} file: empty file"),
            YamlParseError::BadSymbol {file, line_no} => write!(f, "bad entry in file {file} at line {}", line_no + 1),
            YamlParseError::MissingField {file,  symbol , field  } => write!(f, "entry for {symbol} is missing {field} in file {file}"),
            YamlParseError::NoValue {file, symbol , field, line_no } => write!(f, "entry for {symbol} has no value for {field} in file {file} at line {}", line_no + 1),
            YamlParseError::UnparsableLine {file, line_no  } => write!(f, "unexpected contents in file {file} at line {}", line_no + 1),
            YamlParseError::UnexpectedCoordinate {file, line_no  } => write!(f, "found unexpected coordinate block in file {file} at line {} ", line_no + 1),
            YamlParseError::BadCoordinate {file,  line_no  } => write!(f, "bad coordinates in file {file} at line {} ", line_no + 1),
         }
    }
}

// Rust technicalities to handle errors with ? syntax.
impl std::error::Error for YamlParseError {}

impl From<std::io::Error> for YamlParseError {
    fn from(err: std::io::Error) -> Self {
        YamlParseError::Io(err)
    }
}


// TESTS
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_valid_custom_shape() {
        let input = "\
TEST-3:
  id: 1
  symmetry: D3h
  name: Test Trigonal
  vertices:
    - [1.0, 0.0, 0.0]
    - [-0.5, 0.87, 0.0]
    - [-0.5, -0.87, 0.0]
  center:
    - [0.0, 0.0, 0.0]
TEST-4:
  symmetry: D3h2
  name: Test Trigonal2
  vertices:
    - [1.0, 0.0, 0.0]
    - [-0.5, 0.87, 0.0]
    - [-0.5, -0.87, 0.0]
  center:
    - [0.0, 0.0, 0.0]
";
        let result = parse_yaml_str(input, "test.yaml");
        assert!(result.is_ok());

        let shapes = result.unwrap();
        assert_eq!(shapes.len(), 2);
        let shape1 = &shapes[0];
        let shape2 = &shapes[1];
        assert_eq!(shape1.id, 1);
        assert_eq!(shape2.id, 0);
        assert_eq!(shape1.symbol, "TEST-3");
        assert_eq!(shape2.vertices[0], [1.0, 0.0, 0.0]);
    }

    #[test]
    fn missing_field_error() { // Input has missing name field.
        let input = "\
TEST-2:
  symmetry: D3h
  vertices:
    - [1.0, 0.0, 0.0]
    - [-0.5, 0.87, 0.0]
    - [-0.5, -0.87, 0.0]
  center:
    - [0.0, 0.0, 0.0]
";
        let result = parse_yaml_str(input, "test.yaml");
        assert!(result.is_err());
        match result {
            Err(YamlParseError::MissingField { file, symbol, field }) => {
                assert_eq!(file, "test.yaml");
                assert_eq!(symbol, "TEST-2");
                assert_eq!(field, "name");
            }
            other => panic!("unexpected result: {:?}", other),
        }
    }

    #[test]
    fn unexpected_coordinate_error() { // Input has missing name field.
        let input = "\
TEST-2:
  symmetry: D3h
  name: Test Trigonal
    - [1.0, 0.0, 0.0]
  vertices:
    - [-0.5, 0.87, 0.0]
    - [-0.5, -0.87, 0.0]
  center:
    - [0.0, 0.0, 0.0]
";
        let result = parse_yaml_str(input, "test.yaml");
        assert!(result.is_err());
        match result {
            Err(YamlParseError::UnexpectedCoordinate { file, line_no }) => {
                assert_eq!(file, "test.yaml");
                assert_eq!(line_no, 3);
            }
            other => panic!("unexpected result: {:?}", other),
        }
    }

    #[test]
    fn bad_coordinate_error() { // Input has missing name field.
        let input = "\
TEST-2:
  symmetry: D3h
  name: Test Trigonal
  vertices:
    - [1.0, 0,0, 0.0]
    - [-0.5, 0.87, 0.0]
    - [-0.5, -0.87, 0.0]
  center:
    - [0.0, 0.0, 0.0]
";
        let result = parse_yaml_str(input, "test.yaml");
        assert!(result.is_err());
        match result {
            Err(YamlParseError::BadCoordinate { file, line_no }) => {
                assert_eq!(file, "test.yaml");
                assert_eq!(line_no, 4);
            }
            other => panic!("unexpected result: {:?}", other),
        }
    }

    #[test]
    fn unexpected_line_error() { // Input has missing name field.
        let input = "\
TEST-2:
  symmetry: D3h
  banana: banana
  name: Test Trigonal
  vertices:
    - [1.0, 0,0, 0.0]
    - [-0.5, 0.87, 0.0]
    - [-0.5, -0.87, 0.0]
  center:
    - [0.0, 0.0, 0.0]
";
        let result = parse_yaml_str(input, "test.yaml");
        println!("{:?}", result);
        assert!(result.is_err());
        match result {
            Err(YamlParseError::UnparsableLine { file, line_no }) => {
                assert_eq!(file, "test.yaml");
                assert_eq!(line_no, 2);
            }
            other => panic!("unexpected result: {:?}", other),
        }
    }

    #[test]
    fn no_value_error() { // Input has missing name field.
        let input = "\
TEST-2:
  symmetry:
  name: Test Trigonal
  vertices:
    - [1.0, 0,0, 0.0]
    - [-0.5, 0.87, 0.0]
    - [-0.5, -0.87, 0.0]
  center:
    - [0.0, 0.0, 0.0]
";
        let result = parse_yaml_str(input, "test.yaml");
        println!("{:?}", result);
        assert!(result.is_err());
        match result {
            Err(YamlParseError::NoValue { file, symbol, field, line_no }) => {
                assert_eq!(file, "test.yaml");
                assert_eq!(symbol, "TEST-2");
                assert_eq!(field, "symm");
                assert_eq!(line_no, 1);
            }
            other => panic!("unexpected result: {:?}", other),
        }
    }

    #[test]
    fn file_empty_error() { // Input has missing name field.
        let input = "\
\n\n   \n
";
        let result = parse_yaml_str(input, "test.yaml");
        println!("{:?}", result);
        assert!(result.is_err());
        match result {
            Err(YamlParseError::FileEmpty { file}) => {
                assert_eq!(file, "test.yaml");
            }
            other => panic!("unexpected result: {:?}", other),
        }
    }

    #[test]
    fn bad_entry_error() { // Input has missing name field.
        let input = "\
:
  symmetry: D3h
  name: Test Trigonal
  vertices:
    - [1.0, 0.0, 0.0]
    - [-0.5, 0.87, 0.0]
    - [-0.5, -0.87, 0.0]
  center:
    - [0.0, 0.0, 0.0]
";
        let result = parse_yaml_str(input, "test.yaml");
        println!("{:?}", result);
        assert!(result.is_err());
        match result {
            Err(YamlParseError::BadSymbol { file, line_no }) => {
                assert_eq!(file, "test.yaml");
                assert_eq!(line_no, 0);
            }
            other => panic!("unexpected result: {:?}", other),
        }
    }
}
