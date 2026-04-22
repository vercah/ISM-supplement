use std::ops::Range;
use std::path::PathBuf;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RowMode {
    Cooc,
    Diff,
}

#[derive(Debug, Clone)]
pub struct Row {
    pub colors: Range<usize>,
    pub unitig_weight: u64,
    pub kmer_weight: u64,
    pub mode: RowMode,
}

#[derive(Debug, Clone)]
pub struct Dataset {
    pub names: Vec<String>,
    pub colors_flat: Vec<u32>,
    pub rows: Vec<Row>,
    pub n_colors: usize,
    pub total_kmers: u64,
    pub dumpdir: PathBuf,
    pub dataset: String,
    pub k: usize,
    pub outdir: PathBuf,
}

#[derive(Debug, Clone)]
pub struct CoocMarginals {
    pub unitig: Vec<u64>,
    pub kmer: Vec<u64>,
    pub uniq: Vec<u64>,
}

impl CoocMarginals {
    pub fn new(n: usize) -> Self {
        Self {
            unitig: vec![0; n],
            kmer: vec![0; n],
            uniq: vec![0; n],
        }
    }
}

#[derive(Debug, Clone)]
pub struct BlockAcc {
    pub block_i0: usize,
    pub block_i1: usize,
    pub n: usize,
    pub same_unitig: Vec<u64>,
    pub same_kmer: Vec<u64>,
    pub same_uniq: Vec<u64>,
    pub diff_unitig: Vec<u64>,
    pub diff_kmer: Vec<u64>,
    pub diff_uniq: Vec<u64>,
}

impl BlockAcc {
    pub fn new(block_i0: usize, block_i1: usize, n: usize) -> Self {
        let area = (block_i1 - block_i0) * n;
        Self {
            block_i0,
            block_i1,
            n,
            same_unitig: vec![0; area],
            same_kmer: vec![0; area],
            same_uniq: vec![0; area],
            diff_unitig: vec![0; area],
            diff_kmer: vec![0; area],
            diff_uniq: vec![0; area],
        }
    }

    #[inline]
    pub fn row_base(&self, i: usize) -> usize {
        (i - self.block_i0) * self.n
    }

    pub fn merge_from(&mut self, other: Self) {
        debug_assert_eq!(self.same_unitig.len(), other.same_unitig.len());
        add_into(&mut self.same_unitig, &other.same_unitig);
        add_into(&mut self.same_kmer, &other.same_kmer);
        add_into(&mut self.same_uniq, &other.same_uniq);
        add_into(&mut self.diff_unitig, &other.diff_unitig);
        add_into(&mut self.diff_kmer, &other.diff_kmer);
        add_into(&mut self.diff_uniq, &other.diff_uniq);
    }
}

fn add_into(dst: &mut [u64], src: &[u64]) {
    for (d, s) in dst.iter_mut().zip(src) {
        *d += *s;
    }
}
