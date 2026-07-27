use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(name = "rustlearn")]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    FaLength {
        #[arg(short, long)]
        fname: String,
    },

    FaStats {
        #[arg(short, long)]
        fname: String,
    },

    FaFilter {
        /// input fasta
        #[arg(short, long)]
        fname: String,
        /// min length
        #[arg(short = 'm', long)]
        min_len: Option<usize>,
        /// max length
        #[arg(short = 'M', long)]
        max_len: Option<usize>,
        /// minimum GC percentage
        #[arg(short = 'g', long)]
        min_gc: Option<f64>,
    },

    FaFai {
        #[arg(short, long)]
        fname: String,
    }
}