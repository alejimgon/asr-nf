# Ancestral Sequence Reconstruction Pipeline (asr-nf)

This repository provides a reproducible pipeline for ancestral sequence reconstruction using [Nextflow](https://www.nextflow.io/) and [IQ-TREE](http://www.iqtree.org/). The workflow supports iTOL visualization, taxonomy annotation, and extraction of ancestral FASTA sequences for selected nodes.

## Features

- **Automated environment setup** with Conda
- **Ancestral state reconstruction** using IQ-TREE
- **iTOL-compatible outputs** for visualization
- **Flexible FASTA extraction** for all or selected ancestral nodes
- **Taxonomy annotation**  for visualization (optional)

## Requirements

- [Miniconda/Anaconda](https://docs.conda.io/en/latest/)
- [Nextflow](https://www.nextflow.io/) (installed via Conda environment)

## Installation

1. **Clone the repository:**
    ```bash
    git clone https://github.com/alejimgon/asr-nf.git
    cd asr-nf
    ```

2. **Set up the Conda environment:**
    ```bash
    source scripts/setup.sh
    ```
    This will create and activate the `asr-nf` environment with all dependencies.

## Usage

1. **Prepare your input files:**
    - Place your alignment files in `data/alignments/`
    - Place your tree files in `data/phylogenetic_trees/`
    - Edit `data/input_table.txt` with your analysis parameters (see below).

2. **Run the pipeline:**
    ```bash
    nextflow run main.nf
    ```

3. **Optional parameters:**
    - Fetch taxonomy and generate iTOL labels:
        ```bash
        nextflow run main.nf --tax true --email your_email@example.com
        ```
    - Extract FASTA for all or selected nodes:
        ```bash
        nextflow run main.nf --extract_fasta true
        nextflow run main.nf --extract_fasta true --fasta_nodes Node13,Node25
        nextflow run main.nf --extract_fasta true --fasta_nodes "Nodes13|Nodes25"
        nextflow run main.nf --extract_fasta true --fasta_nodes my_nodes.txt
        nextflow run main.nf --extract_fasta true --fasta_missing X
        ```

## Input Table Format

`data/input_table.txt` should be a tab-delimited file with the following columns (no header):

```
alignment_file   model   position   aa   asr_min
```

Example:
```
example.aln   JTT+G   42   K   0.7
```

## Output

- Results are saved in the `outputs/` directory.
- iTOL files and extracted FASTA sequences are organized by analysis.

## Troubleshooting

- **Conda environment not activated:**  
  Always use `source scripts/setup.sh` to ensure the environment is activated in your current shell.
- **Missing dependencies:**  
  Re-run the setup script or check `env/env.yaml`.

## License

This project is for non-commercial use.  
See [LICENSE](LICENSE) for details.

## Citation

If you use this pipeline, please cite [IQ-TREE](http://www.iqtree.org/) and [Nextflow](https://www.nextflow.io/).

---

## Developed
Developed by Alejandro Jiménez-González