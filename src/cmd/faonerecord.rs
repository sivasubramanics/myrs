use crate::data::fasta::FastaWriter;
use crate::data::faidx::FaidxIndex;
use crate::error::{AppError, Result};
use log::{error, info};
use std::path::{Path, PathBuf};
use crate::cmd::fafai;
use crate::utils::helpers::has_valid_extension;
use crate::utils::VALID_FASTA_EXTENSIONS;

pub fn run(
    fname: &str,
    name: &str,
    output: Option<String>,
    fold_width: Option<usize>,
) -> Result<()> {
    let fasta_path = Path::new(fname);

    if !has_valid_extension(fname, VALID_FASTA_EXTENSIONS, false) {
       error!("input file '{}' does not have a valid FASTA extension. Supported extensions: .fasta, .fa, .fna, .faa, .ffn, .frn (not compressed)", fname);
       return Err(AppError::InvalidInput("invalid file extension".into()));
    }

    let fai_path = PathBuf::from(format!("{}.faidx", fname));

    if !fai_path.exists() {
        info!("index file {:?} not found. Generating .faidx index...", fai_path);
        fafai::run(fname)?;
    }

    let out_path = match output {
        Some(path_str) => PathBuf::from(path_str),
        None => derive_output_path(fasta_path, name)?,
    };

    info!(
        "extracting record '{}' from {:?} to {:?}",
        name, fasta_path, out_path
    );

    let index = FaidxIndex::from_path(&fai_path, &fasta_path)?;
    let record = index.fetch(name)?;

    let mut writer = FastaWriter::create_path(&out_path)?;
    writer.write_record(&record, fold_width)?;
    writer.flush()?;

    info!("successfully wrote record '{}' to {:?}", name, out_path);
    Ok(())
}

/// Helper function to construct default output filename (e.g. input.fa -> input.chr01.fa)
fn derive_output_path(input_path: &Path, record_name: &str) -> Result<PathBuf> {
    let parent = input_path.parent().unwrap_or_else(|| Path::new(""));
    let file_name = input_path
        .file_name()
        .and_then(|f| f.to_str())
        .ok_or_else(|| AppError::InvalidInput(format!("invalid input path '{:?}'", input_path)))?;

    let new_filename = if let Some((stem, ext)) = file_name.rsplit_once('.') {
        format!("{}.{}.{}", stem, record_name, ext)
    } else {
        format!("{}.{}", file_name, record_name)
    };

    Ok(parent.join(new_filename))
}