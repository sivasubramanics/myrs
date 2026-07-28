use clap::{Parser, Subcommand};

#[derive(Parser, Debug)]
#[command(
name = "myrs",
author,
version,
about = "(learning) small rust project with utility tools.",
arg_required_else_help = true,
// Custom formatted help screen for top-level `./myrs --help`
// override_help = "(learning) small rust project with utility tools.
//
// Usage: myrs <SUBCOMMAND>
//
// \x1b[4mFASTA Utilities:\x1b[0m
//   fa-length        Calculate and report sequence lengths for records in a FASTA file
//   fa-stats         Generate summary statistics (GC content, N50, sequence counts)
//   fa-filter        Filter FASTA records by sequence length and GC percentage
//   fa-fai           Generate a .faidx index file for fast random access
//   fa-one-record    Extract a single sequence record from an indexed FASTA file
//   fa-some-records  Extract multiple sequence records from a FASTA file
//
// \x1b[4mFASTQ Utilities:\x1b[0m
//   fq-stats         Generate summary statistics for FASTQ files
//
// Options:
// -h, --help       Print help
// -V, --version    Print version"
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Calculate and report sequence lengths for records in a FASTA file
    #[command(arg_required_else_help = true)]
    FaLength {
        /// Input FASTA file path (supports .gz)
        #[arg(short, long)]
        fname: String,
    },

    /// Generate summary statistics (GC content, N50, sequence counts)
    #[command(arg_required_else_help = true)]
    FaStats {
        /// Input FASTA file path (supports .gz)
        #[arg(short, long)]
        fname: String,
    },

    /// Filter FASTA records by sequence length and GC percentage
    #[command(arg_required_else_help = true)]
    FaFilter {
        /// Input FASTA file path (supports .gz)
        #[arg(short, long)]
        fname: String,

        /// Minimum sequence length threshold (bp)
        #[arg(short = 'm', long)]
        min_len: Option<usize>,

        /// Maximum sequence length threshold (bp)
        #[arg(short = 'M', long)]
        max_len: Option<usize>,

        /// Minimum GC percentage threshold (e.g. 45.0 for 45%)
        #[arg(short = 'g', long)]
        min_gc: Option<f64>,
    },

    /// Generate a .faidx index file for fast random access
    #[command(arg_required_else_help = true)]
    FaFai {
        /// Input FASTA file path (uncompressed)
        #[arg(short, long)]
        fname: String,
    },

    /// Extract a single sequence record from an indexed FASTA file
    #[command(arg_required_else_help = true)]
    FaOneRecord {
        /// Input FASTA file path
        #[arg(short, long)]
        fname: String,

        /// Sequence ID or header name to extract
        #[arg(short = 'n', long = "name")]
        name: String,

        /// Output file path (default: <input>.<name>.fa)
        #[arg(short = 'o', long = "out")]
        output: Option<String>,

        /// Line folding/wrapping width for sequence output
        #[arg(short = 'w', long = "width")]
        fold_width: Option<usize>,
    },

    /// Extract multiple sequence records from a FASTA file
    #[command(arg_required_else_help = true)]
    FaSomeRecords {
        /// Input FASTA file path
        #[arg(short, long)]
        fname: String,

        /// Text file containing sequence IDs to extract (one per line)
        #[arg(short = 'q', long = "names-file")]
        names_file: Option<String>,

        /// Sequence IDs to extract (supports space-separated and/or comma-separated values)
        #[arg(
        short = 'n',
        long = "names",
        value_delimiter = ',',
        num_args = 1..,
        value_name = "NAMES"
        )]
        names_list: Option<Vec<String>>,

        /// Output file path (default: <input>.some.fa)
        #[arg(short = 'o', long = "out")]
        output: Option<String>,

        /// Line folding/wrapping width for sequence output
        #[arg(short = 'w', long = "width")]
        fold_width: Option<usize>,
    },

    /// Summary stats for FASTQ files
    #[command(arg_required_else_help = true)]
    FqStats {
        /// Input FASTQ file path (supports .gz)
        #[arg(short, long)]
        fname: String,
        /// Number of threads
        #[arg(short = 't', long = "threads", default_value_t = 2)]
        threads: usize,
    },

    /// Dump KMC database k-mers and their counts to text format
    #[command(arg_required_else_help = true)]
    DumpKMC {
        /// Base prefix name of KMC database (e.g. 'db' for db.kmc_pre / db.kmc_suf)
        #[arg(short = 'p', long = "kmc")]
        fname: String,

        /// Output file path (defaults to stdout if omitted)
        #[arg(short, long)]
        output: Option<String>,

        /// Pre-load suffix buffers completely into memory instead of memory mapping
        #[arg(short, long, default_value_t = false)]
        in_memory: bool,
    },

    #[command(arg_required_else_help = true)]
    FaKmerCOV {
        #[arg(short = 'r', long = "reference")]
        fname: String,
        #[arg(short = 'k', long = "kmc")]
        kmc_db: String,
        #[arg(short = 'o', long = "output")]
        output: String,
        #[arg(short = 'm', long = "memory", default_value_t = false)]
        in_memory: bool,
        #[arg(short = 't', long = "threads", default_value_t = 2)]
        threads: usize,
    }


}