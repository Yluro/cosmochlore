//! `odis` output: the distortion parameters, plus the cshm and csom tables with `--full`.

use super::csv::write_csv;
use super::table::{Row, Table, Width};
use super::{cshm, csom, table_path};
use crate::cli::OdisArgs;
use crate::odis::{OdisResult, OdisRun};
use std::io::Result;

pub fn odis(run: &OdisRun, args: &OdisArgs, silent: bool) -> Result<()> {
    print(&run.result, &args.name, silent);
    if args.table || args.full {
        write_table_csv(&run.result, &args.name, silent)?;
    }
    if let Some(full) = &run.full {
        cshm::print(&full.cshm, &args.name, silent);
        cshm::write_table_csv(&full.cshm, &args.name, silent)?;
        csom::print(&full.csom, &args.name, silent);
        csom::write_table_csv(&full.csom, &args.name, silent)?;
    }
    Ok(())
}

fn print(r: &OdisResult, file_name: &str, silent: bool) {
    let mut table = Table::new(&[Width::Fixed(15), Width::Fixed(13), Width::Fixed(6)]);
    table.row(Row::Text(format!("\nInput file: {file_name}")));
    table.row(Row::DoubleRule);
    table.row(Row::Text(" Octahedral distortion parameters".to_string()));
    let row = |label: &str, value: String, unit: &str| {
        Row::from([label.to_string(), value, unit.to_string()])
    };
    table.row(Row::Rule);
    table.row(row("Mean d(M-X)", format!("{:>12.4}", r.d_mean), "Ang"));
    table.row(row("Zeta", format!("{:>12.4}", r.zeta), "Ang"));
    table.row(row("Delta", format!("{:>12.6}", r.delta), ""));
    table.row(row("Sigma", format!("{:>12.2}", r.sigma), "deg"));
    table.row(row("Theta", format!("{:>12.2}", r.theta), "deg"));
    table.row(row("Volume", format!("{:>12.4}", r.vol), "Ang^3"));
    table.row(Row::Rule);
    table.row(row("Tau", format!("{:>12.2}", r.tau), "deg"));
    table.row(row("Mu", format!("{:>12.2}", r.mu), "Ang"));
    table.row(Row::DoubleRule);
    table.print(silent);
}

/// `<file>_odis_table.csv`: one row with every parameter.
fn write_table_csv(r: &OdisResult, file_name: &str, silent: bool) -> Result<()> {
    let header = [
        "d_mean", "zeta", "delta", "sigma", "theta", "vol", "tau", "mu",
    ];
    let row = vec![
        format!("{:.4}", r.d_mean),
        format!("{:.4}", r.zeta),
        format!("{:.6}", r.delta),
        format!("{:.2}", r.sigma),
        format!("{:.2}", r.theta),
        format!("{:.4}", r.vol),
        format!("{:.4}", r.tau),
        format!("{:.4}", r.mu),
    ];
    let path = table_path(file_name, "odis");
    write_csv(&path, "output table", &header, &[row], silent)
}
