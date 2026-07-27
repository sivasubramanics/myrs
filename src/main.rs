mod cli;
mod commands;
mod data;
mod error;
mod logger;
mod utils;

use clap::Parser;
use crate::logger::init_logger;
use log::info;
use std::time::Instant;

/// Formats a Duration into HH:MM:SS.mmm string format
fn format_duration(duration: std::time::Duration) -> String {
    let total_millis = duration.as_millis();
    let millis = total_millis % 1000;
    let total_secs = total_millis / 1000;
    let seconds = total_secs % 60;
    let minutes = (total_secs / 60) % 60;
    let hours = total_secs / 3600;

    format!("{:02}:{:02}:{:02}.{:03}", hours, minutes, seconds, millis)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let start_time = Instant::now();

    init_logger();
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
    }

    let elapsed = start_time.elapsed();
    info!("Total execution time: {}", format_duration(elapsed));

    Ok(())
}