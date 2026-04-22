use crate::model::{Dataset, Row, RowMode};
use anyhow::{anyhow, bail, Context, Result};
use std::cmp::Ordering;
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, Write};
use std::ops::Range;
use std::path::{Path, PathBuf};

pub fn load_dataset(dumpdir: &Path, dataset: &str, k: usize, outdir: &Path) -> Result<Dataset> {
    let dumpdir = dumpdir.to_path_buf();
    let outdir = outdir.to_path_buf();
    fs::create_dir_all(&outdir)?;

    let color_sets_file = dumpdir.join(format!("{dataset}_k{k}.color_sets.txt"));
    let metadata_file = dumpdir.join(format!("{dataset}_k{k}.metadata.txt"));
    let unitigs_file = dumpdir.join(format!("{dataset}_k{k}.unitigs.fa"));
    let names_file = resolve_names_file(&dumpdir, dataset, k)?;

    for p in [&color_sets_file, &metadata_file, &unitigs_file, &names_file] {
        if !p.exists() {
            bail!("missing required input: {}", p.display());
        }
    }

    let names = load_names(&names_file)?;
    let n_colors = names.len();
    let (mut rows_colors, mut colors_flat) = load_color_sets(&color_sets_file, n_colors)?;
    let total_kmers = accumulate_unitig_weights(&unitigs_file, k, &mut rows_colors)?;
    append_num_kmers_if_missing(&metadata_file, total_kmers)?;

    let rows = rows_colors
        .into_iter()
        .map(|(range, unitig_weight, kmer_weight)| Row {
            colors: range,
            unitig_weight,
            kmer_weight,
            mode: RowMode::Cooc,
        })
        .collect();

    if colors_flat.capacity() == 0 {
        colors_flat.shrink_to_fit();
    }

    Ok(Dataset {
        names,
        colors_flat,
        rows,
        n_colors,
        total_kmers,
        dumpdir,
        dataset: dataset.to_string(),
        k,
        outdir,
    })
}

pub fn classify_rows(ds: &mut Dataset) {
    let n = ds.n_colors;
    for row in &mut ds.rows {
        let s = row.colors.len();
        let c = n.saturating_sub(s);
        let cooc_cost = s.saturating_mul(s.saturating_sub(1)) / 2;
        let diff_cost = s.saturating_mul(c);
        row.mode = if diff_cost < cooc_cost {
            RowMode::Diff
        } else {
            RowMode::Cooc
        };
    }
}

fn resolve_names_file(dumpdir: &Path, dataset: &str, k: usize) -> Result<PathBuf> {
    let candidates = [
        dumpdir.join(format!("{dataset}_k{k}.filenames.txt")),
        PathBuf::from(format!("01_datasets/{dataset}.txt")),
        PathBuf::from(format!("../01_datasets/{dataset}.txt")),
    ];
    for p in candidates {
        if p.exists() {
            return Ok(p);
        }
    }
    bail!(
        "could not locate a names file; expected one of: {dataset}_k{k}.filenames.txt, 01_datasets/{dataset}.txt, ../01_datasets/{dataset}.txt"
    )
}

fn load_names(path: &Path) -> Result<Vec<String>> {
    let file = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut names = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        names.push(canonical_sample_name(trimmed));
    }
    Ok(names)
}

fn canonical_sample_name(value: &str) -> String {
    let stripped = value.trim().trim_end_matches('/');
    let name = stripped.rsplit('/').next().unwrap_or(stripped);
    const FASTA_SUFFIXES: [&str; 6] = [".fa.gz", ".fasta.gz", ".fna.gz", ".fa", ".fasta", ".fna"];
    for suffix in FASTA_SUFFIXES {
        if let Some(prefix) = name.strip_suffix(suffix) {
            return prefix.to_string();
        }
    }
    match name.rfind('.') {
        Some(0) | None => name.to_string(),
        Some(pos) if pos + 1 == name.len() => name[..pos].to_string(),
        Some(pos) => name[..pos].to_string(),
    }
}

type RowWeights = (Range<usize>, u64, u64);

fn load_color_sets(path: &Path, n_colors: usize) -> Result<(Vec<RowWeights>, Vec<u32>)> {
    let file = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut ranges: Vec<Option<RowWeights>> = Vec::new();
    let mut flat: Vec<u32> = Vec::new();
    let mut implicit_id = 0usize;

    for raw_line in reader.lines() {
        let line = raw_line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let mut parts = line.split_whitespace();
        let first = parts.next().ok_or_else(|| anyhow!("empty color-set line"))?;

        let (color_set_id, remainder_tokens) = if let Some(id) = first.strip_prefix("color_set_id=") {
            let parsed_id: usize = id.parse().with_context(|| format!("bad color_set_id in line: {line}"))?;
            let rest: Vec<&str> = parts.collect();
            (parsed_id, rest)
        } else {
            let mut rest = Vec::with_capacity(1 + parts.clone().count());
            rest.push(first);
            rest.extend(parts);
            let parsed_id = implicit_id;
            implicit_id += 1;
            (parsed_id, rest)
        };

        if ranges.len() <= color_set_id {
            ranges.resize_with(color_set_id + 1, || None);
        }

        let start = flat.len();
        let mut colors: Vec<u32> = Vec::new();
        for token in remainder_tokens {
            if token.starts_with("size=") {
                continue;
            }
            let value: usize = token.parse().with_context(|| format!("bad color value in line: {line}"))?;
            if value >= n_colors {
                bail!("color id {value} out of range [0,{n_colors}) in {}", path.display());
            }
            colors.push(value as u32);
        }

        colors.sort_unstable();
        colors.dedup();
        flat.extend(colors);
        let end = flat.len();
        ranges[color_set_id] = Some((start..end, 0, 0));
    }

    let mut out = Vec::with_capacity(ranges.len());
    for (idx, maybe_row) in ranges.into_iter().enumerate() {
        match maybe_row {
            Some(row) => out.push(row),
            None => bail!("missing color_set_id={idx} in {}", path.display()),
        }
    }

    Ok((out, flat))
}

fn accumulate_unitig_weights(path: &Path, k: usize, rows: &mut [RowWeights]) -> Result<u64> {
    let file = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut total_kmers = 0u64;
    let mut current_color_set_id: Option<usize> = None;
    let mut current_len = 0usize;

    let mut flush_current = |color_set_id: Option<usize>, seq_len: usize, total_kmers: &mut u64, rows: &mut [RowWeights]| -> Result<()> {
        if let Some(csid) = color_set_id {
            if csid >= rows.len() {
                bail!("unitigs.fa refers to out-of-range color_set_id={csid}");
            }
            let n_kmers = seq_len.saturating_sub(k).saturating_add(1);
            if n_kmers > 0 {
                *total_kmers += n_kmers as u64;
                rows[csid].1 += 1;
                rows[csid].2 += n_kmers as u64;
            }
        }
        Ok(())
    };

    for raw_line in reader.lines() {
        let line = raw_line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(rest) = line.strip_prefix('>') {
            flush_current(current_color_set_id, current_len, &mut total_kmers, rows)?;
            current_color_set_id = Some(parse_color_set_id_from_header(rest)?);
            current_len = 0;
        } else {
            current_len += line.len();
        }
    }

    flush_current(current_color_set_id, current_len, &mut total_kmers, rows)?;
    Ok(total_kmers)
}

fn parse_color_set_id_from_header(rest: &str) -> Result<usize> {
    for token in rest.split_whitespace() {
        if let Some(value) = token.strip_prefix("color_set_id=") {
            return value
                .parse::<usize>()
                .with_context(|| format!("bad color_set_id in FASTA header: {rest}"));
        }
    }
    bail!("could not find color_set_id=... in FASTA header: {rest}")
}

fn append_num_kmers_if_missing(metadata_file: &Path, total_kmers: u64) -> Result<()> {
    let file = File::open(metadata_file).with_context(|| format!("open {}", metadata_file.display()))?;
    let reader = BufReader::new(file);
    for line in reader.lines() {
        if line?.starts_with("num_kmers=") {
            return Ok(());
        }
    }

    let mut out = OpenOptions::new()
        .append(true)
        .open(metadata_file)
        .with_context(|| format!("append {}", metadata_file.display()))?;
    writeln!(out, "num_kmers={total_kmers}")?;
    Ok(())
}

pub fn colors_slice<'a>(flat: &'a [u32], range: &Range<usize>) -> &'a [u32] {
    &flat[range.clone()]
}

pub fn lower_bound(slice: &[u32], target: u32) -> usize {
    let mut lo = 0usize;
    let mut hi = slice.len();
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        match slice[mid].cmp(&target) {
            Ordering::Less => lo = mid + 1,
            _ => hi = mid,
        }
    }
    lo
}
