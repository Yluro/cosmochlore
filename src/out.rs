use crate::cshm::types::CShMResult;
use crate::csom::prepare::strip_all_labels;
use crate::csom::types::CsomResult;
use crate::data::elements::covalent_radius;
use crate::gidx::GidxResult;
use crate::odis::OdisResult;
use nalgebra::{Matrix3, Vector3};
use std::fs::File;
use std::io::Write;
use std::sync::atomic::{AtomicBool, Ordering};

/// Set by `--silent`; suppresses everything the program prints to stdout.
static SILENT: AtomicBool = AtomicBool::new(false);

pub fn set_silent(silent: bool) {
    SILENT.store(silent, Ordering::Relaxed);
}

/// `println!` that does nothing when `--silent` is active.
macro_rules! say {
    ($($arg:tt)*) => {
        if !SILENT.load(Ordering::Relaxed) {
            println!($($arg)*);
        }
    };
}

/// Escapes a CSV field per RFC 4180: wraps it in double quotes, doubling any embedded quotes,
/// whenever it contains a comma, a quote or a newline.
fn csv_field(s: &str) -> String {
    if s.contains(',') || s.contains('"') || s.contains('\n') || s.contains('\r') {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

// ---------------------------------------------------------------------------------------------
// D10: `Table` and `CsvWriter`. Not wired in yet; the free `print_*_table` / `write_*_csv`
// functions below still do the actual work. Both structs own their `silent` flag, so they do
// not touch the `SILENT` static.
// ---------------------------------------------------------------------------------------------

#[allow(dead_code)]
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Align {
    Left,
    Right,
}

/// How wide a [`Table`] column is.
#[allow(dead_code)]
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Width {
    /// Exactly this many characters.
    Fixed(usize),
    /// The longest cell of the column (header excluded) plus this much padding.
    Fit(usize),
}

#[allow(dead_code)]
#[derive(Clone, Debug)]
pub struct Column {
    header: String,
    width: Width,
    align: Align,
    /// Spaces printed before the column.
    gap: usize,
}

#[allow(dead_code)]
impl Column {
    /// A left-aligned column preceded by one space.
    pub fn new(header: &str, width: Width) -> Self {
        Column {
            header: header.to_string(),
            width,
            align: Align::Left,
            gap: 1,
        }
    }

    pub fn right(mut self) -> Self {
        self.align = Align::Right;
        self
    }

    pub fn gap(mut self, gap: usize) -> Self {
        self.gap = gap;
        self
    }
}

#[allow(dead_code)]
#[derive(Clone, Debug)]
enum Line {
    /// A full-width rule made of this character (`=` or `-`).
    Rule(char),
    /// Free text, printed as is.
    Text(String),
    /// The column headers.
    Header,
    /// One cell per column; cells are formatted by the caller.
    Row(Vec<String>),
}

/// A text table built line by line and rendered once, so column widths, rules and number
/// formatting live in one place. Two shapes use it:
/// - column tables (cshm, csom): `Column`s with headers, one `row` per result;
/// - property tables (odis, gidx): three columns (label, value, unit), one `row` per measure.
#[allow(dead_code)]
#[derive(Clone, Debug)]
pub struct Table {
    silent: bool,
    columns: Vec<Column>,
    rule_width: Option<usize>,
    lines: Vec<Line>,
}

#[allow(dead_code)]
impl Table {
    pub fn new(silent: bool, columns: Vec<Column>) -> Self {
        Table {
            silent,
            columns,
            rule_width: None,
            lines: Vec::new(),
        }
    }

    /// The label / value / unit layout of the odis and gidx tables.
    pub fn properties(silent: bool) -> Self {
        Table::new(
            silent,
            vec![
                Column::new("", Width::Fixed(16)).gap(0),
                Column::new("", Width::Fixed(12)).right().gap(0),
                Column::new("", Width::Fixed(12)).gap(2),
            ],
        )
    }

    /// Forces the width of the rules instead of deriving it from the columns.
    pub fn rule_width(mut self, width: usize) -> Self {
        self.rule_width = Some(width);
        self
    }

    pub fn rule(mut self, ch: char) -> Self {
        self.lines.push(Line::Rule(ch));
        self
    }

    pub fn text(mut self, text: impl Into<String>) -> Self {
        self.lines.push(Line::Text(text.into()));
        self
    }

    pub fn header(mut self) -> Self {
        self.lines.push(Line::Header);
        self
    }

    pub fn row<S: AsRef<str>>(mut self, cells: &[S]) -> Self {
        assert_eq!(cells.len(), self.columns.len(), "one cell per column");
        self.lines.push(Line::Row(
            cells.iter().map(|c| c.as_ref().to_string()).collect(),
        ));
        self
    }

    /// A property row: `label`, `value` printed with `decimals` places, `unit` (may be empty).
    /// The label is indented by one space like the section titles.
    pub fn value(self, label: &str, value: f64, decimals: usize, unit: &str) -> Self {
        self.row(&[
            format!(" {label}"),
            format!("{value:.decimals$}"),
            unit.to_string(),
        ])
    }

    fn widths(&self) -> Vec<usize> {
        self.columns
            .iter()
            .enumerate()
            .map(|(i, col)| match col.width {
                Width::Fixed(w) => w,
                Width::Fit(pad) => {
                    let longest = self
                        .lines
                        .iter()
                        .filter_map(|l| match l {
                            Line::Row(cells) => Some(cells[i].len()),
                            _ => None,
                        })
                        .max()
                        .unwrap_or(0);
                    longest + pad
                }
            })
            .collect()
    }

    fn format_cells<'a>(&self, cells: impl Iterator<Item = &'a str>, widths: &[usize]) -> String {
        let mut out = String::new();
        for ((col, &w), cell) in self.columns.iter().zip(widths).zip(cells) {
            out.push_str(&" ".repeat(col.gap));
            match col.align {
                Align::Left => out.push_str(&format!("{cell:<w$}")),
                Align::Right => out.push_str(&format!("{cell:>w$}")),
            }
        }
        out
    }

    /// The table as text, one `\n`-terminated line per entry.
    pub fn render(&self) -> String {
        let widths = self.widths();
        let rule_width = self.rule_width.unwrap_or_else(|| {
            widths.iter().sum::<usize>() + self.columns.iter().map(|c| c.gap).sum::<usize>()
        });

        let mut out = String::new();
        for line in &self.lines {
            match line {
                Line::Rule(ch) => out.push_str(&ch.to_string().repeat(rule_width)),
                Line::Text(text) => out.push_str(text),
                Line::Header => out.push_str(
                    &self.format_cells(self.columns.iter().map(|c| c.header.as_str()), &widths),
                ),
                Line::Row(cells) => {
                    out.push_str(&self.format_cells(cells.iter().map(String::as_str), &widths))
                }
            }
            out.push('\n');
        }
        out
    }

    /// Prints the table to stdout unless `--silent`.
    pub fn print(&self) {
        if !self.silent {
            print!("{}", self.render());
        }
    }
}

/// `<file>` without a trailing `.xyz`, the stem every output file name is built from.
#[allow(dead_code)]
pub fn output_stem(file_name: &str) -> &str {
    file_name.strip_suffix(".xyz").unwrap_or(file_name)
}

/// Creates an output file and reports it ("Writing {what} to {path}...") unless `silent`.
#[allow(dead_code)]
pub fn create_output(path: &str, what: &str, silent: bool) -> Result<File, std::io::Error> {
    let file = File::create(path)?;
    if !silent {
        println!("Writing {what} to {path}...");
    }
    Ok(file)
}

/// Writes a CSV file one record at a time, escaping every field with [`csv_field`].
#[allow(dead_code)]
pub struct CsvWriter {
    file: File,
}

#[allow(dead_code)]
impl CsvWriter {
    /// Creates `path` and reports it as "Writing {what} to {path}..." unless `silent`.
    pub fn create(path: &str, what: &str, silent: bool) -> Result<Self, std::io::Error> {
        Ok(CsvWriter {
            file: create_output(path, what, silent)?,
        })
    }

    /// Creates `<file stem>_<suffix>.csv`, e.g. `("a.xyz", "cshm_table")` -> `a_cshm_table.csv`.
    pub fn create_for(
        file_name: &str,
        suffix: &str,
        what: &str,
        silent: bool,
    ) -> Result<Self, std::io::Error> {
        let path = format!("{}_{}.csv", output_stem(file_name), suffix);
        Self::create(&path, what, silent)
    }

    /// Writes one record; fields are pre-formatted by the caller (number precision lives there).
    pub fn record<S: AsRef<str>>(&mut self, fields: &[S]) -> Result<(), std::io::Error> {
        let line: Vec<String> = fields.iter().map(|f| csv_field(f.as_ref())).collect();
        writeln!(self.file, "{}", line.join(","))
    }
}

fn format_matrix3(m: &Matrix3<f64>) -> String {
    format!(
        "[{:.4} {:.4} {:.4}; {:.4} {:.4} {:.4}; {:.4} {:.4} {:.4}]",
        m[(0, 0)],
        m[(0, 1)],
        m[(0, 2)],
        m[(1, 0)],
        m[(1, 1)],
        m[(1, 2)],
        m[(2, 0)],
        m[(2, 1)],
        m[(2, 2)],
    )
}

pub fn say_finished(elapsed: std::time::Duration) {
    say!("Program finished in {:?}", elapsed);
}

pub fn welcome_msg() {
    let msg: &str = {
        r"
  _  ______   _____ __  __  ____   _____ _    _ _      ____  _____
 | |/ / __ \ / ____|  \/  |/ __ \ / ____| |  | | |    / __ \|  __ \
 | ' / |  | | (___ | \  / | |  | | |    | |__| | |   | |  | | |__) |
 |  /  |  | |\___ \| |\/| | |  | | |    |  __  | |   | |  | |  _  /
 | . \ |__| |____) | |  | | |__| | |____| |  | | |___| |__| | | \ \
 |_|\_\____/|_____/|_|  |_|\____/ \_____|_|  |_|______\____/|_|  \_\
"
    };
    say!("{}", msg);
    say!("{}", env!("CARGO_PKG_DESCRIPTION"));
    say!("Version: {}", env!("CARGO_PKG_VERSION"));
    say!("Authors: {}", env!("CARGO_PKG_AUTHORS"));
    say!("Repository: {}", env!("CARGO_PKG_REPOSITORY"));
}

pub fn print_cshm_table(results: &[CShMResult], file: &str) {
    say!("\nInput file: {}", file);

    let name_width = results.iter().map(|r| r.name.len()).max().unwrap() + 2;
    let symbol_width = results.iter().map(|r| r.symbol.len()).max().unwrap() + 2;
    let symm_width = "Symmetry".len() + 2;
    let total_width = symbol_width + name_width + symm_width + 7 + 4;

    say!("{}", "=".repeat(total_width));
    say!(
        " {:<sw$} {:<nw$} {:<syw$} {:<7}",
        "Symbol",
        "Shape",
        "Symmetry",
        "CShM",
        sw = symbol_width,
        nw = name_width,
        syw = symm_width
    );
    say!("{}", "-".repeat(total_width));

    for result in results {
        say!(
            " {:<sw$} {:<nw$} {:<syw$} {:<7.3}",
            result.symbol,
            result.name,
            result.symm,
            result.cshm,
            sw = symbol_width,
            nw = name_width,
            syw = symm_width
        );
    }

    let min_s = results
        .iter()
        .map(|r| r.cshm)
        .min_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    say!("{}", "-".repeat(total_width));
    if min_s > 10.0 {
        say!(
            "Only extremely distorted geometries were found for this shape. Make sure the .xyz file is correct."
        )
    }
}

/// Writes the cshm results table (symbol, name, symmetry, s-value)
/// to a `<file>_cshm_table.csv` file.
pub fn write_cshm_csv(results: &[CShMResult], file_name: &str) -> Result<(), std::io::Error> {
    let out_name = file_name
        .strip_suffix(".xyz")
        .unwrap_or(file_name)
        .to_owned()
        + "_cshm_table.csv";
    let mut file = File::create(&out_name)?;

    say!("Writing output table to {}...", out_name);

    writeln!(file, "Symbol,Name,Symmetry,CShM").expect("Unable to write to file.");
    for r in results {
        writeln!(
            file,
            "{},{},{},{:.3}",
            csv_field(&r.symbol),
            csv_field(&r.name),
            csv_field(&r.symm),
            r.cshm
        )?;
    }

    Ok(())
}

/// Writes all ideal reference shapes scaled and rotated to align to the problem structure
/// to a `<file>_ideal.xyz` file.
pub fn write_cshm_reconstructed_xyz(
    file_name: &str,
    results: &[CShMResult],
    labels: &[String],
) -> Result<(), std::io::Error> {
    let out_name = file_name
        .strip_suffix(".xyz")
        .unwrap_or(file_name)
        .to_owned()
        + "_ideal.xyz";
    let mut file = File::create(&out_name)?;

    say!(
        "Writing idealised polyhedra coordinates to table to {}...",
        out_name
    );

    for result in results {
        // Write the preamble to each xyz block.
        writeln!(file, "{}", labels.len())?;
        writeln!(
            file,
            "{} {} CShM = {:.3}",
            result.symbol, result.symm, result.cshm
        )?;

        let mut inverse_perm = vec![0usize; result.perm.len()];
        for (problem_idx, &ref_idx) in result.perm.iter().enumerate() {
            inverse_perm[ref_idx] = problem_idx;
        }

        for (ref_idx, point) in result.xyz.iter().enumerate() {
            let problem_idx = inverse_perm[ref_idx];
            writeln!(
                file,
                "{}  {:.6}  {:.6}  {:.6}",
                labels[problem_idx], point.x, point.y, point.z
            )?;
        }
        writeln!(file)?;
    }

    Ok(())
}

pub fn print_crab() {
    say!(
        r"
     /\
    ( /   @ @    ()
     \\ __| |__  /
      \/   v   \/
     /-|       |-\
    / /-\     /-\ \
     / /-`---'-\ \
      /         \ "
    )
}

pub fn print_odis_table(result: &OdisResult, file: &str) {
    say!("\nInput file: {}", file);
    say!("{}", "=".repeat(34));
    say!(" Octahedral distortion parameters");
    say!("{}", "-".repeat(34));
    say!(
        "{:<16}{:>12.4}  {:<12}",
        " Mean d(M-X)",
        result.d_mean,
        "Ang"
    );
    say!("{:<16}{:>12.4}  {:<12}", " Zeta", result.zeta, "Ang");
    say!("{:<16}{:>12.6}  {:<12}", " Delta", result.delta, "");
    say!("{:<16}{:>12.2}  {:<12}", " Sigma", result.sigma, "deg");
    say!("{:<16}{:>12.2}  {:<12}", " Theta", result.theta, "deg");
    say!("{:<16}{:>12.4}  {:<12}", " Volume", result.vol, "Ang^3");
    say!("{}", "-".repeat(34));
    say!("{:<16}{:>12.2}  {:<12}", " Tau", result.tau, "deg");
    say!("{:<16}{:>12.2}  {:<12}", " Mu", result.mu, "Ang");
    say!("{}", "=".repeat(34));
}

fn print_gidx_angle(label: &str, angle: f64) {
    say!("{:<16}{:>12.2}  {:<12}", label, angle, "deg");
}

fn print_gidx_index(label: &str, index: f64) {
    say!("{:<16}{:>12.3}  {:<12}", label, index, "");
}

pub fn print_gidx_table(result: &GidxResult, file: &str) {
    say!("Input file: {}", file);
    say!("{}", "=".repeat(34));
    match *result {
        GidxResult::Four {
            alpha,
            beta,
            tau4,
            tau4_prime,
        } => {
            say!(" Geometry indices (CN = 4)");
            say!("{}", "-".repeat(34));
            print_gidx_angle(" Alpha", alpha);
            print_gidx_angle(" Beta", beta);
            say!("{}", "-".repeat(34));
            print_gidx_index(" Tau4", tau4);
            print_gidx_index(" Tau4'", tau4_prime);
        }
        GidxResult::Five { alpha, beta, tau5 } => {
            say!(" Geometry indices (CN = 5)");
            say!("{}", "-".repeat(34));
            print_gidx_angle(" Alpha", alpha);
            print_gidx_angle(" Beta", beta);
            say!("{}", "-".repeat(34));
            print_gidx_index(" Tau5", tau5);
        }
        GidxResult::Six {
            alpha1,
            alpha2,
            alpha3,
            tau6,
            gammas,
            tau6_prime,
        } => {
            say!(" Geometry indices (CN = 6)");
            say!("{}", "-".repeat(34));
            print_gidx_angle(" Alpha 1", alpha1);
            print_gidx_angle(" Alpha 2", alpha2);
            print_gidx_angle(" Alpha 3", alpha3);
            for (k, gamma) in gammas.iter().enumerate() {
                print_gidx_angle(&format!(" Gamma {}", k + 1), *gamma);
            }
            say!("{}", "-".repeat(34));
            print_gidx_index(" Tau6", tau6);
            print_gidx_index(" Tau6'", tau6_prime);
        }
        GidxResult::Eight {
            faces,
            dihedral_a,
            dihedral_b,
            tau8_prime,
            cube_delta,
        } => {
            say!(" Geometry indices (CN = 8)");
            say!("{}", "-".repeat(34));
            for (k, face) in faces.iter().enumerate() {
                print_gidx_angle(&format!(" Theta {}", k + 1), face.theta);
                print_gidx_angle(&format!(" Phi {}", k + 1), face.phi);
            }
            print_gidx_angle(" Dihedral A", dihedral_a);
            print_gidx_angle(" Dihedral B", dihedral_b);
            say!("{}", "-".repeat(34));
            for (k, face) in faces.iter().enumerate() {
                print_gidx_index(&format!(" Tau8 (face {})", k + 1), face.tau8);
            }
            print_gidx_index(" Tau8'", tau8_prime);
            print_gidx_index(" Tau8'-Tau8", cube_delta);
        }
    }
    say!("{}", "=".repeat(34));
}

pub fn print_csom_table(results: &[CsomResult], file: &str) {
    say!("\nInput file: {}", file);
    say!("{}", "=".repeat(20));
    say!(" {:<12} {:<10}", "Point group", "CSoM");
    say!("{}", "-".repeat(20));
    for result in results {
        say!(" {:<12} {:<10.3}", result.point_group, result.deviation);
    }
    say!("{}", "-".repeat(20));
}

pub fn write_odis_csv(result: OdisResult, file_name: &str) -> Result<(), std::io::Error> {
    let out_name = file_name
        .strip_suffix(".xyz")
        .unwrap_or(file_name)
        .to_owned()
        + "_odis_table.csv";
    let mut file = File::create(&out_name)?;

    say!("Writing output table to {}...", out_name);

    writeln!(file, "d_mean,zeta,delta,sigma,theta,vol,tau,mu")?;
    writeln!(
        file,
        "{:.4},{:.4},{:.6},{:.2},{:.2},{:.4},{:.4},{:.4}",
        result.d_mean,
        result.zeta,
        result.delta,
        result.sigma,
        result.theta,
        result.vol,
        result.tau,
        result.mu
    )?;

    Ok(())
}

/// Writes the geometry indices (the angles they are built from and the indices for the
/// coordination number) to a `<file>_gidx_table.csv` file.
pub fn write_gidx_csv(result: &GidxResult, file_name: &str) -> Result<(), std::io::Error> {
    let out_name = file_name
        .strip_suffix(".xyz")
        .unwrap_or(file_name)
        .to_owned()
        + "_gidx_table.csv";
    let mut file = File::create(&out_name)?;

    say!("Writing output table to {}...", out_name);

    match *result {
        GidxResult::Four {
            alpha,
            beta,
            tau4,
            tau4_prime,
        } => {
            writeln!(file, "alpha,beta,tau4,tau4_prime")?;
            writeln!(
                file,
                "{:.2},{:.2},{:.4},{:.4}",
                alpha, beta, tau4, tau4_prime
            )?;
        }
        GidxResult::Five { alpha, beta, tau5 } => {
            writeln!(file, "alpha,beta,tau5")?;
            writeln!(file, "{:.2},{:.2},{:.4}", alpha, beta, tau5)?;
        }
        GidxResult::Six {
            alpha1,
            alpha2,
            alpha3,
            tau6,
            gammas,
            tau6_prime,
        } => {
            writeln!(
                file,
                "alpha1,alpha2,alpha3,tau6,gamma1,gamma2,gamma3,gamma4,gamma5,tau6_prime"
            )?;
            writeln!(
                file,
                "{:.2},{:.2},{:.2},{:.4},{:.2},{:.2},{:.2},{:.2},{:.2},{:.4}",
                alpha1,
                alpha2,
                alpha3,
                tau6,
                gammas[0],
                gammas[1],
                gammas[2],
                gammas[3],
                gammas[4],
                tau6_prime
            )?;
        }
        GidxResult::Eight {
            faces,
            dihedral_a,
            dihedral_b,
            tau8_prime,
            cube_delta,
        } => {
            writeln!(
                file,
                "theta1,phi1,tau8_1,theta2,phi2,tau8_2,dihedral_a,dihedral_b,tau8_prime,cube_delta"
            )?;
            writeln!(
                file,
                "{:.2},{:.2},{:.4},{:.2},{:.2},{:.4},{:.2},{:.2},{:.4},{:.4}",
                faces[0].theta,
                faces[0].phi,
                faces[0].tau8,
                faces[1].theta,
                faces[1].phi,
                faces[1].tau8,
                dihedral_a,
                dihedral_b,
                tau8_prime,
                cube_delta
            )?;
        }
    }

    Ok(())
}

/// Writes the csom summary table (point group, deviation, refined rotation matrix) to a
///  `<file>_csom_table.csv` file.
pub fn write_csom_csv(results: &[CsomResult], file_name: &str) -> Result<(), std::io::Error> {
    let out_name = file_name
        .strip_suffix(".xyz")
        .unwrap_or(file_name)
        .to_owned()
        + "_csom_table.csv";
    let mut file = File::create(&out_name)?;

    say!("Writing output table to {}...", out_name);

    writeln!(file, "PointGroup,Dev,Rotation Matrix")?;
    for r in results {
        writeln!(
            file,
            "{},{:.3},{}",
            csv_field(&r.point_group),
            r.deviation,
            format_matrix3(&r.rotation)
        )?;
    }

    Ok(())
}

/// Formats an operation's atom pairing as `Source>Target` entries in the input's atom order:
/// `A>B` means the operation carries atom `A` onto the site of atom `B`.
fn format_pairing(pairing: &[usize], labels: &[String]) -> String {
    // `pairing[target] = source`, so invert it to walk the atoms in input order.
    let mut target_of = vec![0usize; pairing.len()];
    for (target, &source) in pairing.iter().enumerate() {
        target_of[source] = target;
    }

    target_of
        .iter()
        .enumerate()
        .map(|(source, &target)| format!("{}>{}", labels[source], labels[target]))
        .collect::<Vec<String>>()
        .join(" ")
}

/// Writes the individual symmetry operation deviations (name, matrix, s-value, atom pairing)
/// per point group analysed to a  `<file>_<point group>_details.csv` file.
pub fn write_csom_details_csv(
    results: &[CsomResult],
    file_name: &str,
    labels: &[String],
) -> Result<(), std::io::Error> {
    let stem = file_name.strip_suffix(".xyz").unwrap_or(file_name);

    for result in results {
        let out_name = format!("{}_{}_details.csv", stem, result.point_group);
        let mut file = File::create(&out_name)?;

        say!("Writing operation details to {}...", out_name);

        writeln!(file, "name,op_matrix,dev,pairing")?;
        for op in &result.operations {
            writeln!(
                file,
                "{},{},{:.3},{}",
                csv_field(&op.name),
                format_matrix3(&op.matrix),
                op.deviation,
                csv_field(&format_pairing(&op.pairing, labels)),
            )?;
        }
    }

    Ok(())
}

/// Writes, for each point group, a single `<file>_<point group>_operated.xyz` with all
/// reconstructed operated structures by each symmetry element.
pub fn write_csom_operated_xyz(
    results: &[CsomResult],
    file_name: &str,
    labels: &[String],
    original_coords: &[Vector3<f64>],
) -> Result<(), std::io::Error> {
    let stem = file_name.strip_suffix(".xyz").unwrap_or(file_name);

    for result in results {
        let out_name = format!("{}_{}_operated.xyz", stem, result.point_group);
        let mut file = File::create(&out_name)?;

        say!("Writing operated coordinates to {}...", out_name);

        // The identity isn't stored among `result.operations` (the point-group tables in
        // data/pgs.rs omit E -- it trivially gives zero deviation for any structure), so write
        // the untouched original structure as its own "E" block first.
        writeln!(file, "{}", labels.len())?;
        writeln!(file, "{} E dev = 0.000", result.point_group)?;
        for (label, point) in labels.iter().zip(original_coords) {
            writeln!(
                file,
                "{}  {:.6}  {:.6}  {:.6}",
                label, point.x, point.y, point.z
            )?;
        }
        writeln!(file)?;

        // `rotation` is orthogonal, so its transpose is its inverse -- undoes the
        // CSOM-alignment rotation without an explicit matrix inversion.
        let rotation_inv = result.rotation.transpose();

        // `op.image[i]` is where the operation sends atom `i`, so every atom keeps its own
        // label: an atom whose image lands on a *different* atom's site shows up as such.
        for op in &result.operations {
            writeln!(file, "{}", labels.len())?;
            writeln!(
                file,
                "{} {} dev = {:.3}",
                result.point_group, op.name, op.deviation
            )?;
            for (label, point) in labels.iter().zip(&op.image) {
                let real_point = (rotation_inv * point) / result.scale + result.centroid;
                writeln!(
                    file,
                    "{}  {:.6}  {:.6}  {:.6}",
                    label, real_point.x, real_point.y, real_point.z
                )?;
            }
            writeln!(file)?;
        }
    }

    Ok(())
}

/// Two atoms are bonded in the merged .mol2 when they are closer than this factor times the
/// sum of their covalent radii (the same criterion as RDKit's connectivity perception).
const BOND_TOLERANCE: f64 = 1.3;

/// Bonds of one block of the merged .mol2, as 0-based atom index pairs `(i, j)` with `i < j`:
/// every pair within `BOND_TOLERANCE` times the sum of their covalent radii, plus, when the
/// structure has a centre atom (index 0), the centre to every ligand whatever its distance,
/// so the coordination polyhedron is always drawn.
///
/// Each block is a rigid image of the original structure, so one list serves all of them.
fn perceive_bonds(
    elements: &[String],
    coords: &[Vector3<f64>],
    has_centre_atom: bool,
) -> Vec<(usize, usize)> {
    let radii: Vec<f64> = elements.iter().map(|e| covalent_radius(e)).collect();

    let mut bonds = Vec::new();
    for i in 0..coords.len() {
        for j in (i + 1)..coords.len() {
            let cutoff = BOND_TOLERANCE * (radii[i] + radii[j]);
            if (has_centre_atom && i == 0) || (coords[i] - coords[j]).norm() <= cutoff {
                bonds.push((i, j));
            }
        }
    }
    bonds
}

/// Writes, for each point group, a single `<file>_<point group>_merged.mol2` with all
/// the reconstructed operated structures. Bonds come from [`perceive_bonds`].
pub fn write_csom_merged_mol2(
    results: &[CsomResult],
    file_name: &str,
    labels: &[String],
    original_coords: &[Vector3<f64>],
    has_centre_atom: bool,
) -> Result<(), std::io::Error> {
    let stem = file_name.strip_suffix(".xyz").unwrap_or(file_name);
    let elements = strip_all_labels(labels);
    let n_per_block = labels.len();
    let bonds = perceive_bonds(&elements, original_coords, has_centre_atom);

    for result in results {
        let out_name = format!("{}_{}_merged.mol2", stem, result.point_group);
        let mut file = File::create(&out_name)?;

        say!("Writing merged mol2 to {}...", out_name);

        // Each block (the original "E" structure, then one per symmetry operation) as
        // (substructure name, that block's atom coordinates recovered to the original frame).
        let rotation_inv = result.rotation.transpose();
        let mut blocks: Vec<(&str, Vec<Vector3<f64>>)> = vec![("E", original_coords.to_vec())];
        blocks.extend(result.operations.iter().map(|op| {
            let real_points = op
                .image
                .iter()
                .map(|p| (rotation_inv * p) / result.scale + result.centroid)
                .collect();
            (op.name.as_str(), real_points)
        }));

        writeln!(file, "@<TRIPOS>MOLECULE")?;
        writeln!(file, "{}_{}_merged", stem, result.point_group)?;
        writeln!(
            file,
            "{} {} {} 0 0",
            n_per_block * blocks.len(),
            bonds.len() * blocks.len(),
            blocks.len()
        )?;
        writeln!(file, "SMALL")?;
        writeln!(file, "NO_CHARGES")?;
        writeln!(file)?;

        writeln!(file, "@<TRIPOS>ATOM")?;
        let mut atom_id = 0usize;
        for (block_idx, (name, points)) in blocks.iter().enumerate() {
            let subst_id = block_idx + 1;
            for ((label, element), point) in labels.iter().zip(&elements).zip(points) {
                atom_id += 1;
                writeln!(
                    file,
                    "{} {}.{} {:.6} {:.6} {:.6} {} {} {} 0.0000",
                    atom_id, label, name, point.x, point.y, point.z, element, subst_id, name
                )?;
            }
        }

        if !bonds.is_empty() {
            writeln!(file, "@<TRIPOS>BOND")?;
            let mut bond_id = 0usize;
            for block in 0..blocks.len() {
                let base = block * n_per_block;
                for &(i, j) in &bonds {
                    bond_id += 1;
                    writeln!(file, "{} {} {} 1", bond_id, base + i + 1, base + j + 1)?;
                }
            }
        }

        writeln!(file, "@<TRIPOS>SUBSTRUCTURE")?;
        for (block, (name, _)) in blocks.iter().enumerate() {
            writeln!(file, "{} {} {}", block + 1, name, block * n_per_block + 1)?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn csv_field_leaves_plain_text_untouched() {
        assert_eq!(csv_field("vOC-5"), "vOC-5");
    }

    #[test]
    fn csv_field_quotes_a_comma() {
        assert_eq!(csv_field("My Shape, Custom"), "\"My Shape, Custom\"");
    }

    #[test]
    fn csv_field_doubles_embedded_quotes() {
        assert_eq!(csv_field("6\" wide"), "\"6\"\" wide\"");
    }

    #[test]
    fn property_table_matches_the_odis_layout() {
        let table = Table::properties(false)
            .rule_width(34)
            .rule('=')
            .text(" Octahedral distortion parameters")
            .rule('-')
            .value("Mean d(M-X)", 2.0, 4, "Ang")
            .value("Delta", 0.5, 6, "")
            .rule('=')
            .render();
        let expected = format!(
            "{}\n Octahedral distortion parameters\n{}\n{:<16}{:>12.4}  {:<12}\n{:<16}{:>12.6}  {:<12}\n{}\n",
            "=".repeat(34),
            "-".repeat(34),
            " Mean d(M-X)",
            2.0,
            "Ang",
            " Delta",
            0.5,
            "",
            "=".repeat(34),
        );
        assert_eq!(table, expected);
    }

    #[test]
    fn column_table_matches_the_cshm_layout() {
        let table = Table::new(
            false,
            vec![
                Column::new("Symbol", Width::Fit(2)),
                Column::new("Shape", Width::Fit(2)),
                Column::new("Symmetry", Width::Fixed(10)),
                Column::new("CShM", Width::Fixed(7)),
            ],
        )
        .rule('=')
        .header()
        .rule('-')
        .row(&["OC-6", "Octahedron", "Oh", &format!("{:.3}", 0.0)])
        .render();
        let (sw, nw) = (6, 12);
        let expected = format!(
            "{}\n {:<sw$} {:<nw$} {:<10} {:<7}\n{}\n {:<sw$} {:<nw$} {:<10} {:<7.3}\n",
            "=".repeat(sw + nw + 10 + 7 + 4),
            "Symbol",
            "Shape",
            "Symmetry",
            "CShM",
            "-".repeat(sw + nw + 10 + 7 + 4),
            "OC-6",
            "Octahedron",
            "Oh",
            0.0,
        );
        assert_eq!(table, expected);
    }

    #[test]
    fn csv_writer_escapes_fields_and_reports_nothing_when_silent() {
        let path = std::env::temp_dir().join("cosmochlore_csv_writer_test.csv");
        let path = path.to_str().unwrap();
        let mut csv = CsvWriter::create(path, "test table", true).unwrap();
        csv.record(&["a", "b"]).unwrap();
        csv.record(&["x, y", "6\" wide"]).unwrap();
        drop(csv);
        assert_eq!(
            std::fs::read_to_string(path).unwrap(),
            "a,b\n\"x, y\",\"6\"\" wide\"\n"
        );
        std::fs::remove_file(path).unwrap();
    }

    fn labels(symbols: &[&str]) -> Vec<String> {
        symbols.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn perceive_bonds_joins_a_ring_without_a_centre() {
        // A regular hexagon of carbons with 1.39 A edges, like a benzene ring.
        let coords: Vec<Vector3<f64>> = (0..6)
            .map(|k| {
                let angle = k as f64 * std::f64::consts::FRAC_PI_3;
                Vector3::new(1.39 * angle.cos(), 1.39 * angle.sin(), 0.0)
            })
            .collect();

        let bonds = perceive_bonds(&labels(&["C"; 6]), &coords, false);

        // Neighbours only: the 1,3 (2.41 A) and 1,4 (2.78 A) pairs are well past the cutoff.
        assert_eq!(bonds, vec![(0, 1), (0, 5), (1, 2), (2, 3), (3, 4), (4, 5)]);
    }

    #[test]
    fn perceive_bonds_always_joins_the_centre_to_every_ligand() {
        // An octahedron with 3.0 A bonds: past the covalent cutoff for Fe-N (2.9 A), so
        // only the centre rule can draw the polyhedron.
        let mut coords = vec![Vector3::zeros()];
        for axis in [Vector3::x(), Vector3::y(), Vector3::z()] {
            coords.push(3.0 * axis);
            coords.push(-3.0 * axis);
        }
        let elements = labels(&["Fe", "N", "N", "N", "N", "N", "N"]);

        let bonds = perceive_bonds(&elements, &coords, true);
        assert_eq!(bonds, vec![(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6)]);

        // Without a declared centre the same coordinates are simply seven loose atoms.
        assert!(perceive_bonds(&elements, &coords, false).is_empty());
    }

    #[test]
    fn perceive_bonds_uses_the_fallback_radius_for_unknown_labels() {
        // 2.5 A is a bond for a metal-sized dummy atom next to a nitrogen ((1.5 + 0.71) * 1.3
        // = 2.87 A) but not between two nitrogens (1.85 A).
        let coords = vec![Vector3::zeros(), Vector3::new(2.5, 0.0, 0.0)];

        assert_eq!(
            perceive_bonds(&labels(&["Xx", "N"]), &coords, false),
            vec![(0, 1)]
        );
        assert!(perceive_bonds(&labels(&["N", "N"]), &coords, false).is_empty());
    }
}
