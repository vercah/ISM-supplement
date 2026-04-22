mod cli;
mod compute;
mod model;
mod parse;
mod write;

use anyhow::Result;
use clap::Parser;
use cli::Cli;
use compute::{compute_block, compute_cooc_marginals, pick_block_size};
use parse::{classify_rows, load_dataset};
use rayon::ThreadPoolBuilder;
use write::Writers;

fn main() -> Result<()> {
    let cli = Cli::parse();

    ThreadPoolBuilder::new()
        .num_threads(cli.n_cores.max(1))
        .build_global()
        .expect("failed to initialize Rayon thread pool");

    let mut ds = load_dataset(&cli.dumpdir, &cli.dataset, cli.k, &cli.out_prefix)?;
    classify_rows(&mut ds);
    let marg = compute_cooc_marginals(&ds);

    let block_size = cli
        .block_size
        .unwrap_or_else(|| pick_block_size(ds.n_colors, cli.n_cores.max(1), cli.memory_gb));

    eprintln!(
        "loaded {} colors, {} color sets, total_kmers={}, block_size={}",
        ds.n_colors,
        ds.rows.len(),
        ds.total_kmers,
        block_size
    );

    let mut writers = Writers::create(&ds)?;

    let mut block_i0 = 0usize;
    while block_i0 < ds.n_colors {
        let block_i1 = (block_i0 + block_size).min(ds.n_colors);
        eprintln!("computing block [{}:{})", block_i0, block_i1);
        let block = compute_block(&ds, block_i0, block_i1, cli.rows_per_task);
        writers.write_block(&ds, &marg, &block)?;
        block_i0 = block_i1;
    }

    writers.flush()?;
    Ok(())
}
