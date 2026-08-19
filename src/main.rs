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
use std::process;
use std::time::Instant;

fn main() {
    let start_time = Instant::now();

    // initialize logger
    init_logger();

    // run the application
    if let Err(_err) = run_app() {
        log::error!("{}", _err);
        process::exit(1);
    }

    let elapsed = start_time.elapsed();
    info!("total execution time: {}", format_time(elapsed));
}

/// Helper function to dispatch CLI commands and propagate errors with `?`
fn run_app() -> Result<(), Box<dyn std::error::Error>> {
    let cli = cli::Cli::parse();

    // print env version and contact
    let pkg_name = env!("CARGO_PKG_NAME");
    let pkg_version = env!("CARGO_PKG_VERSION");
    let pkg_authors = env!("CARGO_PKG_AUTHORS");

    info!("{} version v{} ({})", pkg_name, pkg_version, pkg_authors);

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

        cli::Commands::FaOneRecord {
            fname,
            name,
            output,
            fold_width,
        } => {
            cmd::faonerecord::run(&fname, &name, output, fold_width)?;
        }

        cli::Commands::FaSomeRecords {
            fname,
            names_file,
            names_list,
            output,
            fold_width,
        } => {
            cmd::fasomerecords::run(
                &fname,
                names_file.as_deref(),
                names_list,
                output,
                fold_width,
            )?;
        }

        cli::Commands::FqStats { fname, threads } => {
            cmd::fqstats::run(&fname, threads)?;
        }

        cli::Commands::DumpKMC {
            fname,
            output,
            in_memory,
        } => {
            cmd::dumpkmc::run(&fname, output, in_memory)?;
        }

        cli::Commands::FaKmerCOV {
            fname,
            kmc_db,
            output,
            in_memory,
            threads,
        } => {
            cmd::fakmercov::run(&fname, &kmc_db, &output, in_memory, threads)?;
        }

        cli::Commands::HicContactMatrix(args) => {
            let filters = cmd::hicmatrix::FilterOptions {
                min_mapq: args.min_mapq,
                max_nm: Some(args.max_nm),
                min_as: Some(args.min_as),
            };
            cmd::hicmatrix::run(&args.input, args.binsize, args.threads, filters)?;
        }
    }

    Ok(())
}