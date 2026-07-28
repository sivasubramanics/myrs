use crate::data::kmc::Kmc;
use crate::error::Result;
use crate::utils::helpers::create_file;
use log::info;
use std::io::BufWriter;
use std::path::{Path, PathBuf};

pub fn run(
    kmc_prefix: &str,
    output: Option<String>,
    in_memory: bool,
) -> Result<()> {
    // Strip `.kmc_pre` or `.kmc_suf` suffix if passed by the user
    let db_base_path = kmc_prefix
        .strip_suffix(".kmc_pre")
        .or_else(|| kmc_prefix.strip_suffix(".kmc_suf"))
        .unwrap_or(kmc_prefix);

    // Determine target output path: use `output` if supplied, else default to `<db_base_path>.dump.tsv`
    let out_path = output
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(format!("{}.dump.tsv", db_base_path)));

    info!("loading KMC database: '{}'", db_base_path);
    let kmc = Kmc::new(db_base_path, in_memory)?;

    info!("dumping k-mer table to output file: '{}'", out_path.display());
    let file = create_file(&out_path)?;
    let writer = BufWriter::new(file);
    kmc.dump(writer)?;

    info!("successfully dumped KMC database to '{}'.", out_path.display());
    Ok(())
}