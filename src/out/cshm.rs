//! `cshm` output: the shape table, `--table` CSV, `--ideal` xyz and the crab.

use super::csv::write_csv;
use super::table::{Row, Table, Width};
use super::{create_output, output_stem, table_path};
use crate::cli::CshmArgs;
use crate::cshm::types::CShMResult;
use std::io::{Result, Write};

pub fn cshm(results: &[CShMResult], args: &CshmArgs, silent: bool) -> Result<()> {
    print(results, &args.name, silent);
    if args.table {
        write_table_csv(results, &args.name, silent)?;
    }
    if args.ideal {
        write_ideal(results, &args.name, silent)?;
    }
    if args.crab {
        crab(silent);
    }
    Ok(())
}

/// The shape table, with a warning when every shape is far off.
pub(super) fn print(results: &[CShMResult], file_name: &str, silent: bool) {
    let widths = [
        Width::Dynamic,
        Width::Dynamic,
        Width::Fixed(10),
        Width::Fixed(7),
    ];
    let mut table = Table::new(&widths);
    table.row(Row::Text(format!("\nInput file: {file_name}")));
    table.row(Row::DoubleRule);
    table.row(Row::from(["Symbol", "Shape", "Symmetry", "CShM"]));
    table.row(Row::Rule);
    for r in results {
        let cshm = format!("{:.3}", r.cshm);
        table.row(Row::from([
            r.symbol.clone(),
            r.name.clone(),
            r.symm.clone(),
            cshm,
        ]));
    }
    table.row(Row::Rule);
    if results.iter().all(|r| r.cshm > 10.0) {
        table.row(Row::Text(
            "Only extremely distorted geometries were found for this shape. Make sure the .xyz file is correct."
                .to_string(),
        ));
    }
    table.print(silent);
}

/// `<file>_cshm_table.csv`: symbol, name, symmetry and CShM per shape.
pub(super) fn write_table_csv(results: &[CShMResult], file_name: &str, silent: bool) -> Result<()> {
    let rows: Vec<Vec<String>> = results
        .iter()
        .map(|r| {
            let cshm = format!("{:.3}", r.cshm);
            vec![r.symbol.clone(), r.name.clone(), r.symm.clone(), cshm]
        })
        .collect();
    let header = ["Symbol", "Name", "Symmetry", "CShM"];
    let path = table_path(file_name, "cshm");
    write_csv(&path, "output table", &header, &rows, silent)
}

/// `<file>_ideal.xyz`: one block per shape with its fitted ideal polyhedron (`CShMResult::ideal`).
fn write_ideal(results: &[CShMResult], file_name: &str, silent: bool) -> Result<()> {
    let path = format!("{}_ideal.xyz", output_stem(file_name));
    let mut file = create_output(&path, "idealised polyhedra coordinates to table", silent)?;
    for result in results {
        let Some(ideal) = &result.ideal else {
            continue;
        };
        let atoms = ideal.atoms();
        writeln!(file, "{}", atoms.len())?;
        writeln!(
            file,
            "{} {} CShM = {:.3}",
            result.symbol, result.symm, result.cshm
        )?;
        for atom in &atoms {
            let p = atom.coords;
            writeln!(file, "{}  {:.6}  {:.6}  {:.6}", atom.label, p.x, p.y, p.z)?;
        }
        writeln!(file)?;
    }
    Ok(())
}

fn crab(silent: bool) {
    say!(
        silent,
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
