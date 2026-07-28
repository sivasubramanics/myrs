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

    #[error("I/O error during processing: {0}")]
    Io(#[from] io::Error),

    #[error("Gzip decoding failed: {0}")]
    GzipDecode(String),

    #[error("Memory mapping failed: {0}")]
    Mmap(io::Error),

    #[error("Input file is compressed. Please decompress before indexing.")]
    CompressedFileNotSupported,

    #[error("Invalid input: {0}")]
    InvalidInput(String),

    #[error("KMC database error: {0}")]
    KmcError(String),

    #[error("Threading error: {0}")]
    ThreadingError(String),

    #[error("Generic error: {0}")]
    Generic(String),
}