mod cli;
mod cmd;
mod data;
mod error;
mod logger;
mod utils;

use clap::Parser;
use crate::logger::init_logger;
use crate::utils::helpers::format_time;
use log::info;
use std::time::Instant;
use std::process;



fn main() {
    let start_time = Instant::now();

    // initialize logger
    init_logger();

    // print env version and contact
    let pkg_name = env!("CARGO_PKG_NAME");
    let pkg_version = env!("CARGO_PKG_VERSION");
    let pkg_authors = env!("CARGO_PKG_AUTHORS");

    info!("{} version v{} ({})", pkg_name, pkg_version, pkg_authors);

    // run the application
    if let Err(_err) = run_app() {
        // since subcommands use error!() before returning an Err,
        // the error message has already been printed via init_logger().
        // simply exiting with a non-zero code.
        process::exit(1);
    }

    let elapsed = start_time.elapsed();
    info!("total execution time: {}", format_time(elapsed));
}

/// Helper function to dispatch CLI commands and propagate errors with `?`
fn run_app() -> Result<(), Box<dyn std::error::Error>> {
    let cli = cli::Cli::parse();

    match cli.command {
        cli::Commands::FaLength { fname } => {
            cmd::falength::run(&fname)?;
        }

        cli::Commands::FaStats { fname } => {
            cmd::fastats::run(&fname)?;
        }

        cli::Commands::FaFilter {
            fname,
            min_len,
            max_len,
            min_gc,
        } => {
            cmd::fafilter::run(&fname, min_len, max_len, min_gc)?;
        }

        cli::Commands::FaFai { fname } => {
            cmd::fafai::run(&fname)?;
        }
    }

    Ok(())
}