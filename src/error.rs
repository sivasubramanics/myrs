use std::io;
use std::path::PathBuf;
use thiserror::Error;

/// Project-wide result alias.
pub type Result<T> = std::result::Result<T, AppError>;

#[derive(Error, Debug)]
pub enum AppError {
    #[error("Unrecognized file extension '{ext}' for file '{path}'.\nSupported extensions: .fasta, .fa, .fna, .faa, .ffn, .frn (optionally with .gz)")]
    InvalidExtension {
        ext: String,
        path: PathBuf,
    },

    #[error("Failed to access file '{path}': {source}")]
    FileIo {
        path: PathBuf,
        #[source]
        source: io::Error,
    },

    #[error("I/O error during FASTA processing: {0}")]
    Io(#[from] io::Error),

    #[error("Gzip decoding failed: {0}")]
    GzipDecode(String),
}