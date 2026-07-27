mod cli;
mod commands;
mod data;
mod utils;
mod error;

use clap::Parser;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = cli::Cli::parse();

    match cli.command {
        cli::Commands::Length { fname } => {
            commands::length::run(&fname)?;
        }

        cli::Commands::Summary { fname } => {
            commands::summary::run(&fname)?;
        }

        cli::Commands::Filter { fname , min_len, max_len, min_gc } => {
            commands::filter::run(&fname, min_len, max_len, min_gc)?;
        }
    }

    Ok(())
}