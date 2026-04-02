# ISM-supplement

Supplementary material and experimental pipelines for the paper:

> **Why phylogenies compress so well: combinatorial guarantees under the Infinite Sites Model**
>
> bioRxiv preprint: [https://www.biorxiv.org/content/10.64898/2026.03.18.712055v1](https://www.biorxiv.org/content/10.64898/2026.03.18.712055v1)

## Prerequisites

### External tools

The following tools must be installed and available in your `PATH`:

- **[Snakemake](https://snakemake.readthedocs.io/)** — workflow management system

  ```bash
  pip install snakemake-minimal
  ```

- **[fulgor](https://github.com/jermp/fulgor)** — colored compacted de Bruijn graph index

  Build from source following the instructions in the [fulgor repository](https://github.com/jermp/fulgor).

- **[attotree](https://github.com/karel-brinda/attotree)** — rapid phylogenetic tree construction

  ```bash
  pip install attotree
  ```

- **[Concorde](https://www.math.uwaterloo.ca/tsp/concorde.html)** — TSP solver

  Must be installed manually. Download from the [Concorde website](https://www.math.uwaterloo.ca/tsp/concorde.html) and ensure the `concorde` binary is in your `PATH`.

### Python packages

```bash
pip install numpy pandas joblib tqdm xopen progressbar2 ete3
```

### Using a Conda environment (recommended)

It may be useful to create a dedicated Conda environment to avoid dependency conflicts:

```bash
conda create -n ism python=3.10
conda activate ism
pip install snakemake-minimal attotree numpy pandas joblib tqdm xopen progressbar2 ete3
```

## Running the experimental pipeline

The example pipeline is located in `experiments/example_pipeline/`.

### Step 0: Clone the repository

```bash
git clone https://github.com/vercah/ISM-supplement.git
cd ISM-supplement
```

If using a Conda environment, activate it before proceeding:

```bash
conda activate ism
```

### Step 1: Generate the dataset file

The pipeline needs a text file listing absolute paths to input genomes (one per line). The included example dataset is `minigono` (located in `experiments/datasets/minigono/`).

To generate the dataset file, run from the `01_datasets/` directory:

```bash
cd experiments/example_pipeline/01_datasets
./generate_data.sh minigono
```

This creates `minigono.txt` with absolute paths to the genome files.

#### Using your own data

You can add your own datasets by placing a folder with genome files (`.fa`, `.fasta`, or `.fna`) into `experiments/datasets/`. Then generate the dataset file by passing the folder name:

```bash
cd experiments/example_pipeline/01_datasets
./generate_data.sh my_dataset
```

### Step 2: Configure pipeline parameters

The pipeline parameters are defined at the top of `experiments/example_pipeline/Snakefile`:

- **`k_values`** — list of k-mer sizes (default: `[31]`)
- **`g_values`** — list of dataset names (default: `["minigono"]`)
- **`matrix_types`** — types of presence/absence matrices (default: `["kmer", "unitig"]`)
- **`N_values`** — number of genomes to evaluate; by default derived from the dataset size

The pipeline generates results for the Cartesian product of all parameter combinations, so specifying e.g. `k_values = [15, 23, 31]` and `g_values = ["minigono", "my_dataset"]` will produce results for all six combinations. If you added a custom dataset in Step 1, make sure to add its name to the `g_values` list.

### Step 3: Run the pipeline

From the `experiments/example_pipeline/` directory, either run the provided script:

```bash
cd experiments/example_pipeline
./RunScript
```

Or run Snakemake directly with the same flags:

```bash
cd experiments/example_pipeline
snakemake -j 32 --resources concurrency=1 --latency-wait 30 --rerun-incomplete --keep-going --show-failed-logs --reason --use-conda
```

All results will be generated in the `11_runs/` directory.
