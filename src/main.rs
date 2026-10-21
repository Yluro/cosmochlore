mod cli;
mod cshm;
mod csom;
mod data;
mod error;
pub mod geometry;
mod gidx;
mod odis;
mod out;
mod xyz;
mod yaml;

use crate::cli::{Cli, Command};
use crate::error::Error;
use clap::Parser;
use std::time::Instant;

fn main() {
    let start = Instant::now();
    let cli = Cli::parse();
    out::banner(cli.silent);
    if let Err(err) = dispatch(cli.command, cli.silent) {
        eprintln!("Error: {err}");
        std::process::exit(1);
    }
    out::finished(start.elapsed(), cli.silent);
}

/// Runs one subcommand and hands its result to `out`.
fn dispatch(command: Command, silent: bool) -> Result<(), Error> {
    match command {
        Command::Cshm(args) => Ok(out::cshm(&cshm::run(&args)?, &args, silent)?),
        Command::Csom(args) => Ok(out::csom(&csom::run(&args)?, &args, silent)?),
        Command::Odis(args) => Ok(out::odis(&odis::run(&args)?, &args, silent)?),
        Command::Gidx(args) => Ok(out::gidx(&gidx::run(&args)?, &args, silent)?),
    }
}
