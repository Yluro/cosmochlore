use nalgebra::Matrix3;
pub struct MeasureResult {
    pub name: String,
    pub symbol: String,
    pub symm: String,
    pub cshm: f64,
    pub rot_mat: Matrix3<f64>
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



pub fn output_table(results: &[MeasureResult]) {

    let name_width = results.iter().map(|r| {r.name.len()}).max().unwrap() + 2;
    let symbol_width = results.iter().map(|r| {r.symbol.len()}).max().unwrap() + 2;
    let total_width = symbol_width + name_width + 7 + 3;

    println!("{}", "=".repeat(total_width));
    println!(" {:<sw$} {:<nw$} {:<7}", "Symbol", "Shape", "CShM", sw=symbol_width, nw = name_width, );
    println!("{}", "-".repeat(total_width));

    for result in results {
        println!(" {:<sw$} {:<nw$} {:<7.3}", result.symbol, result.name, result.cshm, sw = symbol_width, nw = name_width);
    }

    let min_s = results.iter().map(|r| { r.cshm }).min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    println!("{}", "-".repeat(total_width));
    if min_s > 10.0 {println!("Only extremely distorted geometries were found for this shape. Make sure the .xyz file is correct.")}
}