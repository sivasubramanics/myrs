mod cli;
mod commands;
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

    // Initialize your custom logger format
    init_logger();

    // Run the application logic and catch any returned errors
    if let Err(_err) = run_app() {
        // Since subcommands use error!() before returning an Err,
        // the error message has already been printed via init_logger().
        // Simply exit with a non-zero code.
        process::exit(1);
    }

    let elapsed = start_time.elapsed();
    info!("Total execution time: {}", format_time(elapsed));
}

/// Helper function to dispatch CLI commands and propagate errors with `?`
fn run_app() -> Result<(), Box<dyn std::error::Error>> {
    let cli = cli::Cli::parse();

    match cli.command {
        cli::Commands::FaLength { fname } => {
            commands::falength::run(&fname)?;
        }

        cli::Commands::FaStats { fname } => {
            commands::fastats::run(&fname)?;
        }

        cli::Commands::FaFilter {
            fname,
            min_len,
            max_len,
            min_gc,
        } => {
            commands::fafilter::run(&fname, min_len, max_len, min_gc)?;
        }

        cli::Commands::FaFai { fname } => {
            commands::fafai::run(&fname)?;
        }
    }

    Ok(())
}