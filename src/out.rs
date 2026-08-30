use std::fs::File;
use std::io::Write;
use nalgebra::Vector3;
pub struct CShMResult {
    pub name: String,
    pub symbol: String,
    pub symm: String,
    pub cshm: f64,
    pub perm: Vec<usize>,
    pub xyz: Vec<Vector3<f64>>
}

pub fn welcome_msg() {
    let msg: &str = {r"
  _  ______   _____ __  __  ____   _____ _    _ _      ____  _____
 | |/ / __ \ / ____|  \/  |/ __ \ / ____| |  | | |    / __ \|  __ \
 | ' / |  | | (___ | \  / | |  | | |    | |__| | |   | |  | | |__) |
 |  <| |  | |\___ \| |\/| | |  | | |    |  __  | |   | |  | |  _  /
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



pub fn output_table(results: &[CShMResult]) {

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

pub fn write_cshm_csv (file_name: &String, results: &[CShMResult]) {
    let out_name = file_name.strip_suffix(".xyz").unwrap_or(file_name).to_owned() + "_table.csv";
    let mut file = File::create(&out_name)
        .expect("Unable to create file.");

    println!("Writing output table to {}...", out_name);

    writeln!(file, "Symbol,Name,Symmetry,CShM")
        .expect("Unable to write to file.");
    for r in results {
        writeln!(file, "{},{},{},{:.3}", r.symbol, r.name, r.symm, r.cshm)
            .expect("Unable to write to file.");
    }
}

pub fn write_cshm_reconstructed_xyz (file_name: &String, results: &[CShMResult], labels: &[String]) {
    let out_name = file_name.strip_suffix(".xyz").unwrap_or(file_name).to_owned() + "_ideal.xyz";
    let mut file = File::create(&out_name)
        .expect("Unable to create file.");

    println!("Writing idealised polyhedra coordinates to table to {}...", out_name);

    for result in results {
        // Write the preamble to each xyz block.
        writeln!(file, "{}", labels.len()).expect("Unable to write to file.");
        writeln!(file, "{} {} CShM = {:.3}", result.symbol, result.symm, result.cshm).expect("Unable to write to file.");

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
            ).expect("Unable to write to file.");
        }
        writeln!(file, "").expect("Unable to write to file.");
    };
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