# Documentation

## Overview
The 'solution_run.py' script is a Python-based tool designed to conduct statistical analyses on gene knockout and mutation data. It is best run as part of a pipeline also including 'solution_prep.py' and 'solution_split.py'
It generates a clustered heatmap to provide a visual representation of the results. 
This script reads input data from TSV files, performs one-way ANOVA tests on all gene-mutation pairs, and creates a heatmap to illustrate the associations between genes and mutations.

The script also includes a set of unit tests to ensure the correctness of its functions. These tests are defined in the TestSolution class within the 'test_soluttion.py' script.

## Features
Input Data: Accepts transformed mutation data and gene knockout data in TSV file formats.
Statistical Analysis: Conducts one-way ANOVA tests on gene expression changes between mutants and wild-types for all gene-mutation pairs.
Heatmap Generation: Creates a clustered heatmap to visualize the statistical results.


## Usage
To run the script, follow the steps below:
    Ensure you have Docker installed on your system.
    Create the Docker image using the provided 'Dockerfile'.
    Run the Docker container with the required input files and output directory.
    **See the 'solution_installation.txt' file for Docker installation instructions.


Alternatively (if you are in a haste), execute it from the command line with the command below (fingers crossed):
    python solution_run.py -m mutations.tsv -g gene_kos.tsv -o output_directory


## Required Arguments
-m, --mutations: Path to the transposed mutations TSV file.
-g, --gene_kos: Path to the gene knockout TSV file.
-o, --outdir: Output directory where results and heatmap will be saved.

## Optional Arguments
-v, --value_column: Name of the column containing values for plotting ('P-value', 'stats') (default: 'P-value').
-i, --index_column: Name of the column to be used as the heatmap's index (default: 'mutation_id').
-c, --columns_column: Name of the column to be used as the heatmap's columns (default: 'gene_id').
-l, --log-level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL) (default: INFO).

## Output
The script produces the following output:
Statistical Results: A TSV file (gene_vs_mutation_results.tsv) containing statistical results for gene-mutation comparisons.
Heatmap: A PNG image file (gene_vs_mutation_heatmap.png) showing the clustered heatmap.
There could be multiple outputs for these depending on if the input data is split

## Dependencies
Python 3.8 or higher
Required Python packages (specified in requirements.txt):
    argparse
    logging
    os
    scipy
    pandas
    seaborn
    matplotlib


## Docker Support
This script can be containerized using Docker for easy deployment.
A Dockerfile and usage instructions are provided for packaging and running the scripts as a cextflow pipeline within a Docker container.

## Author
Joachim Nwezeobi
