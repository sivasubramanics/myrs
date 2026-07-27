/// Default buffer size for buffered I/O readers/writers (512 KB).
pub const DEFAULT_BUF_SIZE: usize = 512 * 1024;

/// Default line wrap width for FASTA output sequence lines.
pub const DEFAULT_FOLD_WIDTH: usize = 80;

/// Standard FASTA file extensions recognized across bioinformatics tools.
pub const VALID_FASTA_EXTENSIONS: &[&str] = &["fasta", "fa", "fna", "faa", "ffn", "frn"];

/// Standard FASTQ file extensions.
pub const VALID_FASTQ_EXTENSIONS: &[&str] = &["fastq", "fq"];

/// column width
pub const DEFAULT_COLUMN_WIDTH: usize = 30;