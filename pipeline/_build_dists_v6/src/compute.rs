use crate::model::{BlockAcc, CoocMarginals, Dataset, RowMode};
use crate::parse::{colors_slice, lower_bound};
use rayon::prelude::*;

pub fn compute_cooc_marginals(ds: &Dataset) -> CoocMarginals {
    let mut m = CoocMarginals::new(ds.n_colors);
    for row in &ds.rows {
        if row.mode != RowMode::Cooc {
            continue;
        }
        let colors = colors_slice(&ds.colors_flat, &row.colors);
        for &sample in colors {
            let s = sample as usize;
            m.unitig[s] += row.unitig_weight;
            m.kmer[s] += row.kmer_weight;
            m.uniq[s] += 1;
        }
    }
    m
}

pub fn pick_block_size(n_colors: usize, threads: usize, memory_gb: usize) -> usize {
    let budget_bytes = (memory_gb as u128) * 1024u128 * 1024u128 * 1024u128;
    let reserve_fraction = 60u128;
    let usable = budget_bytes * reserve_fraction / 100u128;

    // We allocate 6 block-local u64 arrays per task accumulator.
    // Assume up to `threads` concurrent tasks plus one merge accumulator.
    let consumers = (threads.max(1) + 1) as u128;
    let bytes_per_row = 6u128 * (n_colors as u128) * (std::mem::size_of::<u64>() as u128);
    let derived = if bytes_per_row == 0 || consumers == 0 {
        256usize
    } else {
        (usable / (bytes_per_row * consumers)) as usize
    };

    derived.clamp(64, 1024)
}

pub fn compute_block(ds: &Dataset, block_i0: usize, block_i1: usize, rows_per_task: usize) -> BlockAcc {
    let n = ds.n_colors;
    let init = || BlockAcc::new(block_i0, block_i1, n);

    ds.rows
        .par_chunks(rows_per_task.max(1))
        .map(|chunk| {
            let mut acc = init();
            for row in chunk {
                let colors = colors_slice(&ds.colors_flat, &row.colors);
                match row.mode {
                    RowMode::Cooc => add_cooc_row(&mut acc, colors, row.unitig_weight, row.kmer_weight),
                    RowMode::Diff => add_diff_row(&mut acc, colors, row.unitig_weight, row.kmer_weight),
                }
            }
            acc
        })
        .reduce(init, |mut a, b| {
            a.merge_from(b);
            a
        })
}

fn add_cooc_row(acc: &mut BlockAcc, colors: &[u32], unitig_weight: u64, kmer_weight: u64) {
    if colors.len() < 2 {
        return;
    }
    let lo = lower_bound(colors, acc.block_i0 as u32);
    if lo == colors.len() {
        return;
    }
    let hi = lower_bound(colors, acc.block_i1 as u32);
    if lo >= hi {
        return;
    }

    for p in lo..hi {
        let i = colors[p] as usize;
        let base = acc.row_base(i);
        for &j_u32 in &colors[p + 1..] {
            let j = j_u32 as usize;
            if unitig_weight != 0 {
                acc.same_unitig[base + j] += unitig_weight;
            }
            if kmer_weight != 0 {
                acc.same_kmer[base + j] += kmer_weight;
            }
            acc.same_uniq[base + j] += 1;
        }
    }
}

fn add_diff_row(acc: &mut BlockAcc, colors: &[u32], unitig_weight: u64, kmer_weight: u64) {
    let n = acc.n;
    if colors.is_empty() || colors.len() == n {
        return;
    }

    let mut missing = Vec::with_capacity(n - colors.len());
    let mut expected = 0usize;
    for &c_u32 in colors {
        let c = c_u32 as usize;
        while expected < c {
            missing.push(expected as u32);
            expected += 1;
        }
        expected = c + 1;
    }
    while expected < n {
        missing.push(expected as u32);
        expected += 1;
    }

    let colors_lo = lower_bound(colors, acc.block_i0 as u32);
    let colors_hi = lower_bound(colors, acc.block_i1 as u32);
    for &i_u32 in &colors[colors_lo..colors_hi] {
        let i = i_u32 as usize;
        let base = acc.row_base(i);
        let miss_start = lower_bound(&missing, (i + 1) as u32);
        for &j_u32 in &missing[miss_start..] {
            let j = j_u32 as usize;
            if unitig_weight != 0 {
                acc.diff_unitig[base + j] += unitig_weight;
            }
            if kmer_weight != 0 {
                acc.diff_kmer[base + j] += kmer_weight;
            }
            acc.diff_uniq[base + j] += 1;
        }
    }

    let miss_lo = lower_bound(&missing, acc.block_i0 as u32);
    let miss_hi = lower_bound(&missing, acc.block_i1 as u32);
    for &i_u32 in &missing[miss_lo..miss_hi] {
        let i = i_u32 as usize;
        let base = acc.row_base(i);
        let color_start = lower_bound(colors, (i + 1) as u32);
        for &j_u32 in &colors[color_start..] {
            let j = j_u32 as usize;
            if unitig_weight != 0 {
                acc.diff_unitig[base + j] += unitig_weight;
            }
            if kmer_weight != 0 {
                acc.diff_kmer[base + j] += kmer_weight;
            }
            acc.diff_uniq[base + j] += 1;
        }
    }
}
