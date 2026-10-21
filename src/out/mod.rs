//! Everything cosmochlore prints or writes. `main` hands each subcommand's result here; each
//! command's output lives in its own submodule, built on `table` (terminal) and `csv` (files).

/// `println!` unless `$silent` (`--silent`).
macro_rules! say {
    ($silent:expr, $($arg:tt)*) => {
        if !$silent {
            println!($($arg)*);
        }
    };
}

mod cshm;
mod csom;
mod csv;
mod gidx;
mod odis;
mod table;

pub use self::cshm::cshm;
pub use self::csom::csom;
pub use self::gidx::gidx;
pub use self::odis::odis;
use std::fs::File;

pub fn banner(silent: bool) {
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
    say!(silent, "{}", msg);
    say!(silent, "{}", env!("CARGO_PKG_DESCRIPTION"));
    say!(silent, "Version: {}", env!("CARGO_PKG_VERSION"));
    say!(silent, "Authors: {}", env!("CARGO_PKG_AUTHORS"));
    say!(silent, "Repository: {}", env!("CARGO_PKG_REPOSITORY"));
}

pub fn finished(elapsed: std::time::Duration, silent: bool) {
    say!(silent, "Program finished in {:?}", elapsed);
}

/// `<file>` without a trailing `.xyz`, the stem every output file name is built from.
fn output_stem(file_name: &str) -> &str {
    file_name.strip_suffix(".xyz").unwrap_or(file_name)
}

/// Creates an output file and reports it as "Writing {what} to {path}..." unless `silent`.
fn create_output(path: &str, what: &str, silent: bool) -> std::io::Result<File> {
    let file = File::create(path)?;
    say!(silent, "Writing {what} to {path}...");
    Ok(file)
}

/// `<file stem>_<kind>_table.csv`, the name of a subcommand's `--table` CSV.
fn table_path(file_name: &str, kind: &str) -> String {
    format!("{}_{kind}_table.csv", output_stem(file_name))
}
