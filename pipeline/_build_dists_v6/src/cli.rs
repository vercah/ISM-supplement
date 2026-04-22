use clap::Parser;
use std::path::PathBuf;

#[derive(Debug, Clone, Parser)]
#[command(author, version, about)]
pub struct Cli {
    /// Directory containing Fulgor dump files.
    pub dumpdir: PathBuf,

    /// Dataset basename.
    pub dataset: String,

    /// k-mer length.
    pub k: usize,

    /// Output directory.
    pub out_prefix: PathBuf,

    /// Number of worker threads.
    pub n_cores: usize,

    /// Optional block size. If absent, it is derived from --memory-gb and thread count.
    #[arg(long)]
    pub block_size: Option<usize>,

    /// Approximate RAM budget used to size block-local accumulators.
    #[arg(long, default_value_t = 10)]
    pub memory_gb: usize,

    /// Number of rows processed per Rayon task.
    #[arg(long, default_value_t = 4096)]
    pub rows_per_task: usize,
}
