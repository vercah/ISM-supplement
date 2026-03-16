# ISM-supplement

Supplementary material and experimental pipelines for the paper "Why phylogenies compress so well: combinatorial guarantees under the Infinite Sites Model".

## Running the experimental pipeline

The example pipeline is located in `experiments/example_pipeline/`.

### Prerequisites

The pipeline requires the following tools to be installed and available in your `PATH`:

- [Snakemake](https://snakemake.readthedocs.io/) (with Conda support)
- [fulgor](https://github.com/jermp/fulgor)
- [attotree](https://github.com/karel-brinda/attotree)
- [concorde](https://www.math.uwaterloo.ca/tsp/concorde.html) (TSP solver)

### Steps

1. **Generate dataset file with absolute paths.** Navigate to the `01_datasets/` directory and run the provided script:

   ```bash
   cd experiments/example_pipeline/01_datasets
   bash generate_data.sh
   ```

   This creates a text file listing absolute paths to the input genomes (one per line), which is required by the pipeline.

2. **Run the pipeline.** From the `example_pipeline/` directory, execute:

   ```bash
   cd experiments/example_pipeline
   bash RunScript
   ```
