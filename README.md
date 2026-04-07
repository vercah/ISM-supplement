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
    * [Step 3: Run the pipeline](#step-3-run-the-pipeline)
* [Precomputed results](#precomputed-results)

<!-- vim-markdown-toc -->


## Intro

Supplementary material and experimental pipelines for the paper:

> **Why phylogenies compress so well: combinatorial guarantees under the Infinite Sites Model**
>
> bioRxiv preprint: [https://www.biorxiv.org/content/10.64898/2026.03.18.712055v1](https://www.biorxiv.org/content/10.64898/2026.03.18.712055v1)

## Prerequisites

### Setting up the environment

Create a dedicated Conda environment with all dependencies using [Bioconda](https://bioconda.github.io/):

```bash
conda create -n ism -c bioconda -c conda-forge \
    python=3.12 \
    snakemake fulgor attotree \
    numpy pandas joblib xopen progressbar2 ete3
conda activate ism
```

### Concorde TSP solver

[Concorde](https://www.math.uwaterloo.ca/tsp/concorde.html) is installed via the helper script in `bin/`. The pipeline calls `bin/concorde` directly, so no `PATH` changes are required.

```bash
make -C bin
```

This downloads, builds, and places the `concorde` binary in `bin/concorde`.

## Running the pipeline

The example pipeline is located in `experiments/pipeline/`.


### Step 1: Prepare dataset file(s)

The pipeline auto-detects every `*.txt` file in `experiments/pipeline/01_datasets/`. Each such file is treated as one dataset and must contain absolute paths to input genomes, one per line.

#### Minigono example

The bundled example dataset is `minigono`, stored in `experiments/data/minigono/`. To prepare its file-of-filenames, run:

```bash
cd experiments/pipeline
./prepare_minigono.sh
```

This creates `experiments/pipeline/01_datasets/minigono.txt`.

**Using your own data**

You can add your own datasets by placing a folder with genome files (`.fa`, `.fasta`, or `.fna`) into `experiments/data/`. Then generate the corresponding file-of-filenames in `01_datasets/`:

```bash
cd experiments/pipeline/01_datasets
./generate_data.sh my_dataset
```

This creates `experiments/pipeline/01_datasets/my_dataset.txt`. Any dataset `.txt` file present in `01_datasets/` will be picked up automatically when the pipeline runs.

### Step 2: Configure pipeline parameters

The pipeline parameters are defined at the top of `experiments/pipeline/Snakefile`:

- **`k_values`** — list of k-mer sizes (default: `[31]`)
- **`matrix_types`** — types of presence/absence matrices (default: `["kmer", "unitig"]`)
- **`N_values`** — number of genomes to evaluate; by default derived from the number of lines in each dataset file

The pipeline generates results for the Cartesian product of all parameter combinations across all auto-detected dataset files in `01_datasets/`. For example, if `01_datasets/` contains `minigono.txt` and `my_dataset.txt`, both datasets will be processed automatically.

### Step 3: Run the pipeline

Run the pipeline through `make`:

```bash
cd experiments/pipeline
make
```

All results will be generated in the `11_runs/` directory.

## Precomputed results

The file `final_data_1k.tsv` contains the aggregated results used to produce the figures and tables in the paper.
