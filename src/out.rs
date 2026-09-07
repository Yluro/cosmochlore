use crate::cshm::CShMResult;
use crate::odis::OdisResult;
use crate::csom::CsomResult;
use nalgebra::Matrix3;
use std::fs::File;
use std::io::Write;

fn format_matrix3(m: &Matrix3<f64>) -> String {
    format!(
        "[{:.4} {:.4} {:.4}; {:.4} {:.4} {:.4}; {:.4} {:.4} {:.4}]",
        m[(0, 0)], m[(0, 1)], m[(0, 2)],
        m[(1, 0)], m[(1, 1)], m[(1, 2)],
        m[(2, 0)], m[(2, 1)], m[(2, 2)],
    )
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
    println!("{}", msg);
    println!("{}", env!("CARGO_PKG_DESCRIPTION"));
    println!("Version: {}", env!("CARGO_PKG_VERSION"));
    println!("Authors: {}", env!("CARGO_PKG_AUTHORS"));
    println!("Repository: {}", env!("CARGO_PKG_REPOSITORY"));
}



pub fn print_cshm_table(results: &[CShMResult], file: &str) {

    println!("\nInput file: {}", file);

    let name_width = results.iter().map(|r| {r.name.len()}).max().unwrap() + 2;
    let symbol_width = results.iter().map(|r| {r.symbol.len()}).max().unwrap() + 2;
    let symm_width = "Symmetry".len() + 2;
    let total_width = symbol_width + name_width + symm_width + 7 + 4;

    println!("{}", "=".repeat(total_width));
    println!(" {:<sw$} {:<nw$} {:<syw$} {:<7}", "Symbol", "Shape", "Symmetry", "CShM", sw=symbol_width, nw = name_width, syw=symm_width);
    println!("{}", "-".repeat(total_width));

    for result in results {
        println!(" {:<sw$} {:<nw$} {:<syw$} {:<7.3}", result.symbol, result.name, result.symm, result.cshm, sw = symbol_width, nw = name_width, syw=symm_width);
    }

    let min_s = results.iter().map(|r| { r.cshm }).min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    println!("{}", "-".repeat(total_width));
    if min_s > 10.0 {println!("Only extremely distorted geometries were found for this shape. Make sure the .xyz file is correct.")}
}

/// Writes the cshm results table (symbol, name, symmetry, s-value)
/// to a `<file>_cshm_table.csv` file.
pub fn write_cshm_csv (results: &[CShMResult], file_name: &String) -> Result<(), std::io::Error> {
    let out_name = file_name.strip_suffix(".xyz").unwrap_or(file_name).to_owned() + "_cshm_table.csv";
    let mut file = File::create(&out_name)?;

    println!("Writing output table to {}...", out_name);

    writeln!(file, "Symbol,Name,Symmetry,CShM")
        .expect("Unable to write to file.");
    for r in results {
        writeln!(file, "{},{},{},{:.3}", r.symbol, r.name, r.symm, r.cshm)?;
    }

    Ok(())
}

/// Writes all ideal reference shapes scaled and rotated to align to the problem structure
/// to a `<file>_ideal.xyz` file.
pub fn write_cshm_reconstructed_xyz (file_name: &String, results: &[CShMResult], labels: &[String]) -> Result<(), std::io::Error> {
    let out_name = file_name.strip_suffix(".xyz").unwrap_or(file_name).to_owned() + "_ideal.xyz";
    let mut file = File::create(&out_name)
        .expect("Unable to create file.");

    println!("Writing idealised polyhedra coordinates to table to {}...", out_name);

    for result in results {
        // Write the preamble to each xyz block.
        writeln!(file, "{}", labels.len())?;
        writeln!(file, "{} {} CShM = {:.3}", result.symbol, result.symm, result.cshm)?;

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
        writeln!(file, "")?;
    };

    Ok(())
}

pub fn print_crab() {
    println!(r"
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
    println!("\nInput file: {}", file);
    println!("{}", "=".repeat(34));
    println!(" Octahedral distortion parameters");
    println!("{}", "-".repeat(34));
    println!("{:<16}{:>12.4}  {:<12}", " Mean d(M-X)", result.d_mean, "Ang");
    println!("{:<16}{:>12.4}  {:<12}", " Zeta", result.zeta, "Ang");
    println!("{:<16}{:>12.6}  {:<12}", " Delta", result.delta, "");
    println!("{:<16}{:>12.2}  {:<12}", " Sigma", result.sigma, "deg");
    //println!("{:<16}{:>12.2}  {:<12}", " Theta", result.theta, "deg");
    //println!("{:<16}{:>12.4}  {:<12}", " Volume", result.volume, "Ang^3");
    println!("{}", "-".repeat(34));
    println!("{:<16}{:>12.2}  {:<12}", " Tau", result.tau, "deg");
    println!("{:<16}{:>12.2}  {:<12}", " Mu", result.mu, "Ang");
    println!("{}", "=".repeat(34));
}

pub fn print_csom_table(results: &[CsomResult], file: &str) {
    println!("\nInput file: {}", file);
    println!("{}", "=".repeat(24));
    println!(" {:<12} {:<10}", "Point group", "Deviation");
    println!("{}", "-".repeat(24));
    for result in results {
        println!(" {:<12} {:<10.3}", result.point_group, result.deviation);
    }
    println!("{}", "-".repeat(24));
}

pub fn write_odis_csv(result: OdisResult, file_name: &str) -> Result<(), std::io::Error> {
    let out_name = file_name.strip_suffix(".xyz").unwrap_or(file_name).to_owned() + "_odis_table.csv";
    let mut file = File::create(&out_name)?;

    println!("Writing output table to {}...", out_name);

    writeln!(file, "d_mean,zeta,delta,sigma,tau,mu")?;
    writeln!(file, "{:.4},{:.4},{:.6},{:.2},{:.4},{:.4}",
        result.d_mean,
        result.zeta,
        result.delta,
        result.sigma,
        result.tau,
        result.mu
    )?;

    Ok(())
}

/// Writes the csom summary table (point group, deviation, refined rotation matrix) to a
///  `<file>_csom_table.csv` file.
pub fn write_csom_csv(results: &[CsomResult], file_name: &str) -> Result<(), std::io::Error> {
    let out_name = file_name.strip_suffix(".xyz").unwrap_or(file_name).to_owned() + "_csom_table.csv";
    let mut file = File::create(&out_name)?;

    println!("Writing output table to {}...", out_name);

    writeln!(file, "PointGroup,Dev,Rotation Matrix")?;
    for r in results {
        writeln!(file, "{},{:.3},{}", r.point_group, r.deviation, format_matrix3(&r.rotation))?;
    }

    Ok(())
}

/// Writes the individual symmetry operation deviations (name, matrix, s-value) per point group analysed
/// to a  `<file>_<point group>_details.csv` file.
pub fn write_csom_details_csv(results: &[CsomResult], file_name: &str) -> Result<(), std::io::Error> {
    let stem = file_name.strip_suffix(".xyz").unwrap_or(file_name);

    for result in results {
        let out_name = format!("{}_{}_details.csv", stem, result.point_group);
        let mut file = File::create(&out_name)?;

        println!("Writing operation details to {}...", out_name);

        writeln!(file, "name,op_matrix,dev")?;
        for (name, matrix, dev) in &result.operations {
            writeln!(file, "{},{},{:.3}", name, format_matrix3(matrix), dev)?;
        }
    }

    Ok(())
}