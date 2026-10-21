//! CSV files. A stopgap: each command lists its header and rows by hand until the CSV writer is
//! designed. Fields are written unquoted: the parsers reject `xyz::FORBIDDEN_CHARS`.

use super::create_output;
use std::io::Write;

/// Writes `header` and `rows`, comma-separated, to `path`.
pub fn write_csv(
    path: &str,
    what: &str,
    header: &[&str],
    rows: &[Vec<String>],
    silent: bool,
) -> std::io::Result<()> {
    let mut file = create_output(path, what, silent)?;
    writeln!(file, "{}", header.join(","))?;
    for row in rows {
        debug_assert_eq!(row.len(), header.len(), "one field per column");
        writeln!(file, "{}", row.join(","))?;
    }
    Ok(())
}
