# ISM-supplement


<!-- vim-markdown-toc GFM -->

* [Prerequisites](#prerequisites)
    * [Setting up the environment](#setting-up-the-environment)
    * [Concorde TSP solver](#concorde-tsp-solver)
* [Running the pipeline](#running-the-pipeline)
    * [Step 1: Prepare dataset file(s)](#step-1-prepare-dataset-files)
        * [Minigono example](#minigono-example)
        * [Using your own data](#using-your-own-data)
    * [Step 2: Configure pipeline parameters](#step-2-configure-pipeline-parameters)
        * [Output order names](#output-order-names)
    * [Step 3: Run the pipeline](#step-3-run-the-pipeline)
* [Precomputed results](#precomputed-results)
* [Pipeline workarounds and implementation notes](#pipeline-workarounds-and-implementation-notes)
    * [Concorde TSP solver](#concorde-tsp-solver)
    * [Phylogenetic tree inference and ordering extraction](#phylogenetic-tree-inference-and-ordering-extraction)
    * [Fulgor indexing](#fulgor-indexing)
    * [Randomized ordering](#randomized-ordering)

<!-- vim-markdown-toc -->


## Intro

Supplementary material and experimental pipelines for the paper:

> Veronika Hendrychová and Karel Břinda (2026). **Why phylogenies compress so well: combinatorial guarantees
> under the Infinite Sites Model**. *bioRxiv* 2026.03.18.712055, 2026.
> [https://doi.org/10.64898/2026.03.18.712055](https://doi.org/10.64898/2026.03.18.712055)

BibTeX:

```bibtex
@article{hendrychova2026ism,
  author = {Veronika Hendrychová and Karel Břinda},
  title = {Why phylogenies compress so well: combinatorial guarantees under the Infinite Sites Model},
  journal = {bioRxiv},
  volume = {2026.03.18.712055},
  year = {2026},
  doi = {10.64898/2026.03.18.712055},
  url = {https://doi.org/10.64898/2026.03.18.712055}
}
```

## Prerequisites

### Setting up the environment

Create a dedicated Conda environment with all dependencies using [Bioconda](https://bioconda.github.io/):

```bash
git clone https://github.com/vercah/ism-supplement
cd ism-supplement

# create the ism environment with conda
conda env create -f environment.yml
# alternative faster route: via mamba
mamba env create -f environment.yml

conda activate ism
```

Now test if Fulgor functions correctly by `fulgor help`. If you see an error like `GLIBCXX_3.4.21' not found`
(happens on some older Linux systems), run additionally

```bash
conda install libstdcxx-ng
```

The repository already vendors `galitime` `0.3.0` from the upstream GitHub release as the standalone
executable `bin/galitime`, so no separate `galitime` installation is required.


### Concorde TSP solver

[Concorde](https://www.math.uwaterloo.ca/tsp/concorde.html) is installed via the helper script in `bin/`. The pipeline calls `bin/concorde` directly, so no `PATH` changes are required.

```bash
make -C bin
```

This downloads, builds, and places the `concorde` binary in `bin/concorde`.

## Running the pipeline

The example pipeline is located in `pipeline/`.


### Step 1: Prepare dataset file(s)

The pipeline auto-detects every `*.txt` file in `pipeline/01_datasets/`. Each such file is treated as one dataset and must contain absolute paths to input genomes, one per line.

Supported genome inputs are `.fa`, `.fasta`, `.fna`, and their `.gz` variants.

#### Minigono example

The bundled example dataset is `minigono`, stored in `data/minigono/`. `prepare_minigono.sh` scans that directory for supported FASTA files and writes a dataset manifest:

```bash
cd pipeline
./prepare_minigono.sh
```

This creates `pipeline/01_datasets/minigono.txt`.

**Using your own data**

`link_genomes.sh` supports three input modes:

- a dataset folder name under `data/`
- a directory path containing FASTA files
- one or more explicit FASTA paths, including shell-expanded globs

If you place a folder with genome files in `data/`, you can generate the corresponding file-of-filenames in `01_datasets/` like this:

```bash
cd pipeline/01_datasets
./link_genomes.sh my_dataset
```

This creates `pipeline/01_datasets/my_dataset.txt`.

You can also point `link_genomes.sh` directly at an external directory or an explicit shell-expanded file list:

```bash
cd pipeline/01_datasets
./link_genomes.sh ~/tmp/neisseria_gonorrhoeae__01
./link_genomes.sh ~/tmp/neisseria_gonorrhoeae__01/SAM*.fa
```

In that mode, the dataset name is derived from the common parent directory, so the examples above create `pipeline/01_datasets/neisseria_gonorrhoeae__01.txt`.

Any dataset `.txt` file present in `01_datasets/` will be picked up automatically when the pipeline runs.

### Step 2: Configure pipeline parameters

Edit the single workflow config file `pipeline/config.yaml`. This is the only intended user-facing pipeline config file.

```yaml
# ====================================================================
# Pipeline configuration
# ====================================================================
# The Snakefile loads this file via `configfile:` and validates every
# setting at startup. Unknown keys under any mapping are rejected.
# All bool-maps require at least one entry to be true.
# ====================================================================

# k-mer sizes to build Fulgor indices for.
kmer_sizes:
  - 31

# Seed for every deterministic shuffle step:
# initial dataset randomization and the `randomized` baseline order.
random_seed: 1

# Geometric subsampling density for automatically derived N values.
#
# 0 - evaluate only the full dataset size.
# 1 - evaluate approximately one value per decade: 1, 10, 100, ...
# 2 - evaluate approximately two values per decade: 1, 3, 10, 32, ...
#
# The full dataset size is always included.
sampling_conf: 0

# Matrix types to build from the Fulgor dump.
matrices:
  kmer: true
  unitig: true
  uniqrow: true

# Genome orders to generate and evaluate.
orders:
  best_tsp: true
  worst_tsp: true
  nj_attotree: true
  upgma_attotree: true
  nj_hamming: true
  upgma_hamming: true
  randomized: true

# Passed to every `attotree` invocation as `-k` and `-s`.
attotree_k: 31
attotree_s: 10000
```

Configuration keys:

- **`kmer_sizes`** — list of k-mer sizes to evaluate
- **`random_seed`** — base seed controlling reproducible dataset shuffling and randomized output ordering
- **`sampling_conf`** — geometric subsampling density for automatically derived `N` values
- **`matrices`** — enable or disable the `kmer`, `unitig`, and `uniqrow` matrix families
- **`orders`** — enable or disable the final ordering families:
  `best_tsp`, `worst_tsp`, `nj_attotree`, `upgma_attotree`, `nj_hamming`, `upgma_hamming`, `randomized`
- **`attotree_k`** — attotree k-mer size passed as `attotree -k`
- **`attotree_s`** — attotree sketch size passed as `attotree -s`

`sampling_conf` controls which subset sizes are evaluated for each dataset:

- `0` preserves the previous behavior and evaluates only the full dataset size.
- `1` evaluates one logarithmic sample per decade: `1, 10, 100, ...`, plus the full dataset size.
- `2` evaluates two logarithmic samples per decade: `1, round(sqrt(10)), 10, round(sqrt(10) * 10), 100, ...`, plus the full dataset size.

More generally, a positive value `c` means `c` logarithmic steps per order of magnitude. Candidate values are computed with standard Python `round(...)`, clamped to the dataset size, deduplicated after rounding, and returned in ascending order. The final full dataset size is always included. For example, dataset size `7` gives `[7]` when `sampling_conf: 0` and `[1, 7]` when `sampling_conf: 1`; dataset size `100` gives `[1, 10, 100]` when `sampling_conf: 1` and `[1, 3, 10, 32, 100]` when `sampling_conf: 2`.

Selections remain nested: `02_order_randomization/{dataset}.txt` is the reproducibly shuffled master list, and `03_selection/{dataset}.N{N}.txt` is the first `N` entries from that list, sorted for easier debugging. Therefore larger `N` selections still contain the smaller selected set. When `N = 1`, the pipeline writes the single selected genome directly as the order file in `09_orders/` for each enabled ordering family, skips Concorde and attotree, and still runs the normal final evaluation so `10_runs/..._N1_...runs` is produced.

Two common edits:

Disable worst-case TSP orders:

```yaml
orders:
  worst_tsp: false
```

Disable the unitig matrix family:

```yaml
matrices:
  unitig: false
```

The pipeline generates results for the Cartesian product of the enabled settings across all auto-detected dataset files in `01_datasets/`. `N_values` are derived automatically from the number of lines in each dataset file and the `sampling_conf` setting.

### Output order names

The `method` column in `final_data.tsv` and the method component in filenames use the following order tokens:

| Method token | Meaning | Source distances |
|---|---|---|
| `optimal` | Minimum-runs ordering from Concorde TSP | Pipeline Hamming distances |
| `worst` | Maximum-runs ordering from Concorde TSP on inverted distances | Pipeline Hamming distances |
| `randomized` | Deterministic seed-based random ordering | None |
| `nj-attotree` | Neighbor-Joining tree order from `attotree` | AttoTree sketch distances from FASTA files |
| `upgma-attotree` | UPGMA tree order from `attotree` | AttoTree sketch distances from FASTA files |
| `nj-hamming` | Neighbor-Joining tree order from `quicktree` | Pipeline Hamming distance matrix |
| `upgma-hamming` | UPGMA tree order from `quicktree` | Pipeline Hamming distance matrix |

The final per-run filename pattern is:

```text
10_runs/{dataset}_{method}_k{k}_N{N}_{type}.runs
```

where `{type}` is one of `kmer`, `unitig`, or `uniqrow`.

Use hyphens, not underscores, inside method tokens. Underscores delimit fields in generated filenames.

### Step 3: Run the pipeline

After editing `pipeline/config.yaml`, run the pipeline through `make` from `pipeline/`:

```bash
cd pipeline
make
```

Per-order run results will be generated in the `10_runs/` directory. The aggregate run-count TSV is generated as `final_data.tsv`, and the aggregated timing TSV is generated as `final_time.tsv`.

## Precomputed results

The file `final_data_1k.tsv` contains the aggregated results used to produce the figures and tables in the paper.


## Pipeline workarounds and implementation notes

The pipeline includes several workarounds to handle limitations of the tools it relies on, as well as implementation choices that are worth documenting.

### Concorde TSP solver

#### Cycle-to-path conversion via a dummy separator node

Concorde solves the *cycle* version of TSP, but the pipeline needs a *path*. To convert between the two, a dummy node called `_DUMMY_CITY_SEPARATOR_` with zero-weight edges to all other nodes is prepended to every instance (`_export_tsp_instance.py`). This guarantees that the optimal cycle passes through the dummy for free, so removing it yields an optimal open path. The extraction script (`_extract_path_from_tsp_solution.py`) asserts the dummy is the first node in the solution and strips it.

#### Padding small instances by duplicating the first genome

For instances with fewer than 35 nodes, Concorde switches internally to a Held-Karp dynamic-programming routine (`CCheldkarp_small`) that imposes strict limits on edge-weight magnitude, which can cause solver failures. To force Concorde onto its general branch-and-cut code path, for small selected sets the pipeline repeats the first real genome until the final TSP instance, including the dummy separator node, has at least 35 nodes (`_export_tsp_instance.py`). Each copy has zero distance to itself and to other copies, and inherits the original genome's distances to all other real nodes. In any optimal tour, the copies therefore cluster together as consecutive neighbors and can be collapsed back into a single node after solving, recovering the original ordering (`_extract_path_from_tsp_solution.py`). Importantly, the duplicated node must be a real genome, not the dummy separator node -- duplicating the separator would create additional zero-cost "teleportation" shortcuts between arbitrary points in the tour, distorting the optimal path.

#### Distance scaling to prevent arithmetic overflow

When maximum pairwise Hamming distances exceed 10,000, the resulting optimal tour lengths can overflow Concorde's internal integer arithmetic (error `OVERFLOW in CCbigguy_addmult`). To prevent this, when any of the computed distances is `max_dist > 10000`, the pipeline uniformly scales positive distances as `d' = max(1, round(d / (max_dist / 10000)))` and leaves zero distances unchanged (`_export_tsp_instance.py`). The `max(1, ...)` floor ensures no positive distance collapses to zero. This approximation may merge close distances into the same integer, but does not significantly affect results given the large original range.

#### Worst-case instance by distance inversion

To find the *worst* (maximum-runs) ordering using Concorde's minimization, all distances are inverted: `d' = D_max - d` (`_export_tsp_instance.py`). This converts the maximization problem into a minimization problem. The theoretical justification is given in the supplement of the preprint (Note S2).


### Phylogenetic tree inference and ordering extraction

The pipeline evaluates two independent tree-based ordering families.

The `nj-attotree` and `upgma-attotree` orders run `attotree` directly on the selected genome FASTA files. These orders do not use the Hamming distance matrices generated by the pipeline; `attotree` computes its own sketch-based distances internally. The method passed to `attotree -m` is `nj` or `upgma`.

The `nj-hamming` and `upgma-hamming` orders first convert the selected pairwise Hamming distance matrix to PHYLIP format with `_dists_to_phylip.py`, then run `quicktree` on that matrix. The `upgma-hamming` order passes `-upgma` to `quicktree`; the `nj-hamming` order uses quicktree's Neighbor-Joining mode.

For both tree families, the final genome order is the left-to-right leaf order extracted from the Newick tree by `_tree_leaf_order.py`. The pipeline does not apply rooting, leaf-order normalization, or multifurcation resolution.


### Fulgor indexing

#### File descriptor limit

Fulgor may open many files simultaneously during index construction. The pipeline raises the soft file-descriptor limit to 4,096 before calling `fulgor build` (`Snakefile`).

#### Distance-builder wrapper

The workflow calls `pipeline/_build_dists`. This file is a symlink created by `make` and points to the compiled Rust v6 implementation at `pipeline/_build_dists_v6/_build_dists_v6`.


### Randomized ordering

The pipeline now uses deterministic seed-based shuffling in `_shuffle_lines_with_seed.py` rather than `sort -R`. A single `random_seed` config value controls both dataset shuffling and randomized output ordering.
