# _build_dists_v6

Rust implementation of the Fulgor dump distance builder.

## What it computes

Given a Fulgor dump, it writes three pairwise distance files:

- `{dataset}_k{k}_unitig.dists.txt`
- `{dataset}_k{k}_kmer.dists.txt`
- `{dataset}_k{k}_uniqrow.dists.txt`

For a color set row `S` and weight `w`, the contribution to the Hamming distance between colors `i` and `j` is `w` iff exactly one of `i,j` belongs to `S`.

The three weightings are:

- `unitig` – one per dumped unitig
- `kmer` – number of kmers implied by unitig length
- `uniqrow` – one per distinct dumped color set

## Why this version exists

The Python v1/v2 implementation uses dense `num_colors × num_colors` matrices and updates rows against their complements, which is fundamentally memory-heavy and slow at 10k colors.

The C++ v4/v5 code switches to blockwise output and the identity

`distance(i,j) = marginal(i) + marginal(j) - 2 * cooccurrence(i,j)`

which removes the need for a full global dense matrix, but it still pays too much on dense rows and can multiply memory by the number of parallel workers.

This Rust v6 uses an adaptive hybrid:

- **cooccurrence mode** for sparse rows – enumerate only `S choose 2`
- **direct-difference mode** for dense rows – enumerate only `S × complement(S)`
- **streamed block writing** – never materialize a full `n × n` matrix

That is the right shape for datasets where rows are either very sparse or nearly universal.

## Input format

Supported color-set line formats:

- legacy / repo style: `color_set_id=17 size=3 1 4 9`
- current Fulgor dump style: `size=3 1 4 9`

Supported unitig header styles:

- `> color_set_id=17`
- `> unitig_id=123 color_set_id=17`

The code prefers `{dumpdir}/{dataset}_k{k}.filenames.txt` when present and otherwise falls back to `01_datasets/{dataset}.txt` and `../01_datasets/{dataset}.txt`.

## Build

```bash
cd pipeline/_build_dists_v6
make
```

## Run

```bash
./_build_dists_v6 05_dumps ngono 31 06_ham_dists 10
```

Optional knobs:

```bash
./_build_dists_v6 05_dumps ngono 31 06_ham_dists 10 \
  --block-size 512 \
  --rows-per-task 4096 \
  --memory-gb 10
```

## Notes

- For 10k colors, the output itself is very large: `10,000 * 9,999 / 2 = 49,995,000` pairs per metric.
- Text formatting and compression are a substantial part of end-to-end time.
- The implementation is designed to stay under the requested RAM budget rather than to maximize temporary in-memory buffering.
