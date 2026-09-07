mod xyz;
mod cli;
mod yaml;
mod shapes;
mod data;
mod cshm;
mod csom;
mod odis;
mod coordinates;
mod out;
pub mod geometry;

use crate::cli::{Cli, Command};
use crate::out::welcome_msg;
use clap::Parser;
use std::time::Instant;


fn main() -> Result<(), Box<dyn std::error::Error>> {
    let main_start = Instant::now();
    let args = Cli::parse();
    welcome_msg();
    
    let run = match args.command {
        Command::Cshm(cshm_args) => { cshm::cshm_main(cshm_args) },
        Command::Csom(csom_args) => { csom::csom_main(csom_args) },
        Command::Odis(odis_args) => { odis::main_odis(odis_args) },
    };

    if let Err(err) = run {
        eprintln!("Error: {err}");
        std::process::exit(1);
    }

    let main_time = main_start.elapsed();
    println!("Program finished in {:?}", main_time);
    //if args.crab {print_crab()}
    run
}


