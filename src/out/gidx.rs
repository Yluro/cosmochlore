//! `gidx` output: the geometry indices table and its `--table` CSV.

use super::csv::write_csv;
use super::table::{Row, Table, Width};
use super::table_path;
use crate::cli::GidxArgs;
use crate::gidx::GidxResult;
use std::io::Result;

pub fn gidx(result: &GidxResult, args: &GidxArgs, silent: bool) -> Result<()> {
    print(result, &args.name, silent);
    if args.table {
        write_table_csv(result, &args.name, silent)?;
    }
    Ok(())
}

fn print(result: &GidxResult, file_name: &str, silent: bool) {
    let angle = |label: &str, v: f64| {
        Row::from([label.to_string(), format!("{v:>12.2}"), "deg".to_string()])
    };
    let index =
        |label: &str, v: f64| Row::from([label.to_string(), format!("{v:>12.3}"), String::new()]);
    let cn = match result {
        GidxResult::Four { .. } => 4,
        GidxResult::Five { .. } => 5,
        GidxResult::Six { .. } => 6,
        GidxResult::Eight { .. } => 8,
    };

    let mut table = Table::new(&[Width::Fixed(15), Width::Fixed(13), Width::Fixed(6)]);
    table.row(Row::Text(format!("\nInput file: {file_name}")));
    table.row(Row::DoubleRule);
    table.row(Row::Text(format!(" Geometry indices (CN = {cn})")));
    table.row(Row::Rule);
    match *result {
        GidxResult::Four {
            alpha,
            beta,
            tau4,
            tau4_prime,
        } => {
            table.row(angle("Alpha", alpha));
            table.row(angle("Beta", beta));
            table.row(Row::Rule);
            table.row(index("Tau4", tau4));
            table.row(index("Tau4'", tau4_prime));
        }
        GidxResult::Five { alpha, beta, tau5 } => {
            table.row(angle("Alpha", alpha));
            table.row(angle("Beta", beta));
            table.row(Row::Rule);
            table.row(index("Tau5", tau5));
        }
        GidxResult::Six {
            alpha1,
            alpha2,
            alpha3,
            tau6,
            gammas,
            tau6_prime,
        } => {
            table.row(angle("Alpha 1", alpha1));
            table.row(angle("Alpha 2", alpha2));
            table.row(angle("Alpha 3", alpha3));
            for (n, gamma) in (1..).zip(gammas) {
                table.row(angle(&format!("Gamma {n}"), gamma));
            }
            table.row(Row::Rule);
            table.row(index("Tau6", tau6));
            table.row(index("Tau6'", tau6_prime));
        }
        GidxResult::Eight {
            faces,
            dihedral_a,
            dihedral_b,
            tau8_prime,
            cube_delta,
        } => {
            for (n, face) in (1..).zip(faces) {
                table.row(angle(&format!("Theta {n}"), face.theta));
                table.row(angle(&format!("Phi {n}"), face.phi));
            }
            table.row(angle("Dihedral A", dihedral_a));
            table.row(angle("Dihedral B", dihedral_b));
            table.row(Row::Rule);
            for (n, face) in (1..).zip(faces) {
                table.row(index(&format!("Tau8 (face {n})"), face.tau8));
            }
            table.row(index("Tau8'", tau8_prime));
            table.row(index("Tau8'-Tau8", cube_delta));
        }
    }
    table.row(Row::DoubleRule);
    table.print(silent);
}

/// `<file>_gidx_table.csv`: one row with the angles (2 decimals) and indices (4 decimals).
fn write_table_csv(result: &GidxResult, file_name: &str, silent: bool) -> Result<()> {
    let a = |v: f64| format!("{v:.2}");
    let i = |v: f64| format!("{v:.4}");
    let (header, row): (Vec<&str>, Vec<String>) = match *result {
        GidxResult::Four {
            alpha,
            beta,
            tau4,
            tau4_prime,
        } => (
            vec!["alpha", "beta", "tau4", "tau4_prime"],
            vec![a(alpha), a(beta), i(tau4), i(tau4_prime)],
        ),
        GidxResult::Five { alpha, beta, tau5 } => (
            vec!["alpha", "beta", "tau5"],
            vec![a(alpha), a(beta), i(tau5)],
        ),
        GidxResult::Six {
            alpha1,
            alpha2,
            alpha3,
            tau6,
            gammas,
            tau6_prime,
        } => {
            let mut header = vec!["alpha1", "alpha2", "alpha3", "tau6"];
            header.extend([
                "gamma1",
                "gamma2",
                "gamma3",
                "gamma4",
                "gamma5",
                "tau6_prime",
            ]);
            let mut row = vec![a(alpha1), a(alpha2), a(alpha3), i(tau6)];
            row.extend(gammas.map(a));
            row.push(i(tau6_prime));
            (header, row)
        }
        GidxResult::Eight {
            faces,
            dihedral_a,
            dihedral_b,
            tau8_prime,
            cube_delta,
        } => {
            let mut header = vec!["theta1", "phi1", "tau8_1", "theta2", "phi2", "tau8_2"];
            header.extend(["dihedral_a", "dihedral_b", "tau8_prime", "cube_delta"]);
            let mut row = Vec::new();
            for face in faces {
                row.extend([a(face.theta), a(face.phi), i(face.tau8)]);
            }
            row.extend([a(dihedral_a), a(dihedral_b), i(tau8_prime), i(cube_delta)]);
            (header, row)
        }
    };
    let path = table_path(file_name, "gidx");
    write_csv(&path, "output table", &header, &[row], silent)
}
