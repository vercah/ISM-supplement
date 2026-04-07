# ISM-supplement


<!-- vim-markdown-toc GFM -->

* [Prerequisites](#prerequisites)
    * [Setting up the environment](#setting-up-the-environment)
    * [Concorde TSP solver](#concorde-tsp-solver)
* [Running the pipeline](#running-the-pipeline)
    * [Step 1: Generate the dataset file](#step-1-generate-the-dataset-file)
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
    snakemake fulgor attotree \
    numpy pandas joblib tqdm xopen progressbar2 ete3
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


### Step 1: Generate the dataset file

The pipeline needs a text file listing absolute paths to input genomes (one per line).

**Minigono example:** The included example dataset is `minigono` (located in `experiments/data/minigono/`).

To generate the dataset file, run from the `01_datasets/` directory:

```bash
cd experiments/pipeline/01_datasets
./generate_data.sh minigono
```

This creates `minigono.txt` with absolute paths to the genome files.

**Using your own data**

You can add your own datasets by placing a folder with genome files (`.fa`, `.fasta`, or `.fna`) into `experiments/data/`. Then generate the dataset file by passing the folder name:

```bash
cd experiments/pipeline/01_datasets
./generate_data.sh my_dataset
```

### Step 2: Configure pipeline parameters

The pipeline parameters are defined at the top of `experiments/pipeline/Snakefile`:

- **`k_values`** — list of k-mer sizes (default: `[31]`)
- **`g_values`** — list of dataset names (default: `["minigono"]`)
- **`matrix_types`** — types of presence/absence matrices (default: `["kmer", "unitig"]`)
- **`N_values`** — number of genomes to evaluate; by default derived from the dataset size

The pipeline generates results for the Cartesian product of all parameter combinations, so specifying e.g. `k_values = [15, 23, 31]` and `g_values = ["minigono", "my_dataset"]` will produce results for all six combinations. If you added a custom dataset in Step 1, make sure to add its name to the `g_values` list.

### Step 3: Run the pipeline

Run the pipeline through `make`:

```bash
cd experiments/pipeline
make
```

All results will be generated in the `11_runs/` directory.

## Precomputed results

The file `final_data_1k.tsv` contains the aggregated results used to produce the figures and tables in the paper.
