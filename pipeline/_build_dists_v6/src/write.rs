use crate::model::{BlockAcc, CoocMarginals, Dataset};
use anyhow::{Context, Result};
use itoa::Buffer;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

pub struct Writers {
    unitig: BufWriter<File>,
    kmer: BufWriter<File>,
    uniq: BufWriter<File>,
}

impl Writers {
    pub fn create(ds: &Dataset) -> Result<Self> {
        let unitig = BufWriter::with_capacity(16 * 1024 * 1024, create_file(&ds.outdir.join(format!("{}_k{}_unitig.dists.txt", ds.dataset, ds.k)))?);
        let kmer = BufWriter::with_capacity(16 * 1024 * 1024, create_file(&ds.outdir.join(format!("{}_k{}_kmer.dists.txt", ds.dataset, ds.k)))?);
        let uniq = BufWriter::with_capacity(16 * 1024 * 1024, create_file(&ds.outdir.join(format!("{}_k{}_uniqrow.dists.txt", ds.dataset, ds.k)))?);
        Ok(Self { unitig, kmer, uniq })
    }

    pub fn flush(&mut self) -> Result<()> {
        self.unitig.flush()?;
        self.kmer.flush()?;
        self.uniq.flush()?;
        Ok(())
    }

    pub fn write_block(&mut self, ds: &Dataset, marg: &CoocMarginals, block: &BlockAcc) -> Result<()> {
        let n = ds.n_colors;
        let mut buf_unitig = Vec::<u8>::with_capacity(8 * 1024 * 1024);
        let mut buf_kmer = Vec::<u8>::with_capacity(8 * 1024 * 1024);
        let mut buf_uniq = Vec::<u8>::with_capacity(8 * 1024 * 1024);
        let name_prefixes: Vec<String> = ds.names.iter().map(|s| format!("{s}\t")).collect();

        for i in block.block_i0..block.block_i1 {
            let base = block.row_base(i);
            for j in (i + 1)..n {
                let unitig_dist = block.diff_unitig[base + j]
                    + marg.unitig[i]
                    + marg.unitig[j]
                    - 2 * block.same_unitig[base + j];
                let kmer_dist = block.diff_kmer[base + j]
                    + marg.kmer[i]
                    + marg.kmer[j]
                    - 2 * block.same_kmer[base + j];
                let uniq_dist = block.diff_uniq[base + j]
                    + marg.uniq[i]
                    + marg.uniq[j]
                    - 2 * block.same_uniq[base + j];

                append_line(&mut buf_unitig, &name_prefixes[i], &ds.names[j], unitig_dist);
                append_line(&mut buf_kmer, &name_prefixes[i], &ds.names[j], kmer_dist);
                append_line(&mut buf_uniq, &name_prefixes[i], &ds.names[j], uniq_dist);
            }
        }

        self.unitig.write_all(&buf_unitig)?;
        self.kmer.write_all(&buf_kmer)?;
        self.uniq.write_all(&buf_uniq)?;
        Ok(())
    }
}

fn create_file(path: &Path) -> Result<File> {
    File::create(path).with_context(|| format!("create {}", path.display()))
}

fn append_line(buf: &mut Vec<u8>, left_prefix: &str, right: &str, value: u64) {
    buf.extend_from_slice(left_prefix.as_bytes());
    buf.extend_from_slice(right.as_bytes());
    buf.push(b'\t');
    let mut itoa = Buffer::new();
    buf.extend_from_slice(itoa.format(value).as_bytes());
    buf.push(b'\n');
}
