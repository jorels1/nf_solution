#!/usr/bin/python

import argparse
import logging
import os
import sys
import scipy
from scipy.stats import f_oneway
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from solution_prep import setup_logging



#set default constants
DEFAULT_COLORMAP = 'coolwarm_r'
DEFAULT_CLUSTERING_METHOD = 'complete'



def get_input_dataframes(mutations, gene_kos):

    """
    This function reads the mutation and gene knockout data from TSV files and return them as dataframes.

    Parameters:
    - mutations (str): The file path to the transposed mutations data in TSV format.
    - gene_kos (str): The file path to the gene knockout data in TSV format.

    Returns:
    - tuple: A tuple containing two pandas dataframes: the mutation and gene knockout data.
    """

    logging.info(f'Loading the input data: {mutations} and {gene_kos}')

    # Read the mutation data
    mutations_T_df = pd.read_csv(mutations, sep='\t')

    # Read the gene knockout data
    genes_df = pd.read_csv(gene_kos, sep='\t')

    return mutations_T_df, genes_df


def generate_gene_and_mutation_lists(genes_df, mutations_T_df):
    
    """
    This function extracts gene and mutation IDs from the dataframes and returns them as lists

    Parameters:
    - genes_df (pd.DataFrame): DataFrame containing gene knockout data.
    - mutations_T_df (pd.DataFrame): Transposed DataFrame containing mutation status data.

    Returns:
    - tuple: A tuple containing two lists: (gene_list, mutation_list).
    """

    logging.info(f'Generating the gene and mutation lists from the input data')

    # Generate the gene list (excluding the first column which is assumed to be an identifier)
    gene_list = list(genes_df.columns[1:])

    # Generate the mutation list (excluding the first column which is assumed to be an identifier)
    mutation_list = list(mutations_T_df.columns[1:])

    return gene_list, mutation_list


def perform_one_way_anova(genes_df, mutations_T_df, gene_list, mutation_list, outdir):
    
    """
    This function performs a statistical analysis comparing the gene expression changes between mutants and wild-types
    for multiple pairs of genes and mutations.

    Parameters:
    - genes_df (pd.DataFrame): DataFrame containing gene expression data.
    - mutations_T_df (pd.DataFrame): DataFrame containing transposed mutation data.
    - gene_list (list): List of gene IDs to analyze.
    - mutation_list (list): List of mutation IDs to analyze.
    - outdir (str): Directory to write output file

    Returns:
    - pd.DataFrame: DataFrame containing statistical results for gene-mutation comparisons.
    """

    logging.info(f'Performing a one way anova')

    # Create an empty list to store the results
    output = []

    # Iterate through each gene and mutation combination
    for gene in gene_list:
        for mutation in mutation_list:
            mutants_fold_changes = genes_df.loc[mutations_T_df[mutation] == 1][gene]
            wt_fold_changes = genes_df.loc[mutations_T_df[mutation] == 0][gene]
            statistic, p_value = f_oneway(mutants_fold_changes, wt_fold_changes)

            # Append results as dictionaries to the output list
            output.append({
                'gene_id': gene,
                'mutation_id': mutation,
                'stats': statistic,
                'P-value': p_value
            })

    # Create a DataFrame from the output list
    output_df = pd.DataFrame(output)

    return output_df


def create_clustered_heatmap(output_df, value_column, index_column, columns_column, cmap=DEFAULT_COLORMAP, 
                             method=DEFAULT_CLUSTERING_METHOD, title=None):
    
    """
    This function takes in a dataframe, converts it into a pivot table, and creates a clustered heatmap from the pivot table.

    Parameters:
    - output_df (pd.DataFrame): The input DataFrame containing the statistical results computed above.
    - value_column (str): The name of the column with the values to be plotted.
    - index_column (str): The name of the column to be used as the index for the heatmap.
    - columns_column (str): The name of the column to be used as the columns for the heatmap.
    - cmap (str, optional): The colormap to be used for coloring the heatmap. Default is DEFAULT_COLORMAP.
    - method (str, optional): The clustering method to be used. Default is DEFAULT_CLUSTERING_METHOD.
    - title (str, optional): The title of the heatmap. Default is None.

    Returns:
    - plt.figure: The generated heatmap figure.
    """

    logging.info(f'Creating a clustered heatmap using the analysis result')

    # Pivot the data to create a matrix
    heatmap_data = output_df.pivot_table(values=value_column, index=index_column, columns=columns_column)

    # Set font scale (adjust font size as needed)
    sns.set(font_scale=1)

    # Create a clustered heatmap
    sns.clustermap(heatmap_data, cmap=cmap, method=method, figsize=(10, 6), annot=True, fmt='.3f')

    # Customize the title
    if title:
        title = f'Clustered Heatmap {value_column} for the Gene-Mutation Analysis'
        plt.title(title, loc='left')

    return plt


def solution(args):

    """
    This function wraps all the above function to load the gene knockout and mutation data, perform statistical analysis, and creates a heatmap.

    Parameters:
    - args (argparse.Namespace): An argparse namespace containing command-line arguments.

    Returns:
    - None
    """

    # Set up logging and create output directories
    setup_logging(log_level=args.log_level)
    os.makedirs(args.outdir, exist_ok=True)

    # Load the input data
    mutations_T_df, genes_df = get_input_dataframes(mutations=args.mutations, gene_kos=args.gene_kos)

    # Extract the gene and mutation ID lists from the genes_df and mutations_T_df
    gene_list, mutation_list = generate_gene_and_mutation_lists(genes_df=genes_df, mutations_T_df=mutations_T_df)

    # Perform the one-way ANOVA on the dataset
    output_df = perform_one_way_anova(genes_df=genes_df, mutations_T_df=mutations_T_df, gene_list=gene_list, 
                mutation_list=mutation_list, outdir=args.outdir)

    # Plot the heatmap
    heatmap_figure = create_clustered_heatmap(output_df=output_df, value_column=args.value_column, 
                    index_column=args.index_column, columns_column=args.columns_column)

    # Extract the chunk number for mutations and gene_kos from filename
    mutation_chunk_number = os.path.basename(args.mutations).split('_')[-1].replace('.tsv', '')
    gene_chunk_number = os.path.basename(args.gene_kos).split('_')[-1].replace('.tsv', '')

    # Create a combined unique identifier to add as a suffix to the output filename
    combined_chunk_number = f"mutations_{mutation_chunk_number}_vs_genes_{gene_chunk_number}"


    # Save the statistcial results and heatmap
    try:
        # Save the statistical results CSV
        results_csv_path = os.path.join(args.outdir, f'gene_vs_mutation_results_{combined_chunk_number}.tsv')
        output_df.to_csv(results_csv_path, index=False, sep='\t')
        logging.info(f"Statistical results saved to {results_csv_path}")

        # Save the heatmap
        heatmap_path = os.path.join(args.outdir, f'gene_vs_mutation_heatmap_{combined_chunk_number}.png')
        heatmap_figure.savefig(heatmap_path, dpi=300)
        logging.info(f"Heatmap saved to {heatmap_path}")

    except Exception as e:
        error_message = f"Error: An exception occurred while saving output files: {e}"
        logging.error(error_message)
        raise RuntimeError(error_message)

    logging.info('solution.py is complete')



def main():

    """
    Entry point for the gene-mutation analysis script.

    This function parses command-line arguments to specify input files, output settings, and logging options.
    It then calls the 'solution_run()' function to perform the analysis and generate results.

    Command-Line Arguments:
    - -m, --mutations: Mutations TSV file inherited from solution_prep.py (required).
    - -g, --gene_kos: Gene KOs TSV file path (required).
    - -v, --value_column: Name of the column for plotting values ('P-value', 'stats') (default: 'P-value').
    - -i, --index_column: Name of the index column (default: 'mutation_id').
    - -c, --columns_column: Name of the columns for the heatmap (default: 'gene_id').
    - -o, --outdir: Output directory path (default: 'output').
    - -l, --log-level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL) (default: INFO).

    Returns:
    - None
    """
    
    parser = argparse.ArgumentParser(description="A python script to perform one way ANOVA tests across all (gene, mutation) pairs")
    parser.add_argument('-m', '--mutations', help='Mutations tsv file', required=True)
    parser.add_argument('-g', '--gene_kos', help='Gene KOs tsv file', required=True)
    parser.add_argument('-v', '--value_column', help="The name of th column with the values to be plotted ('P-value', 'stats'). \
                        Default: 'P-value'", default='P-value')
    parser.add_argument('-i', '--index_column', help='The name of the column to be used as the index. \
                        Default: mutation_id', default='mutation_id')
    parser.add_argument('-c', '--columns_column', help='The name of the column to be used as the columns. \
                        Default: gene_id', default='gene_id')
    parser.add_argument('-o', '--outdir', help='Output directory. Default: output', default='output')
    parser.add_argument('-l', '--log-level', help='Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL). \
                        Default: INFO', default='INFO')


    args = parser.parse_args()

    # Validate input file paths
    if not os.path.exists(args.mutations):
        error_message = f"Error: Mutations file '{args.mutations}' does not exist."
        logging.error(error_message)
        raise ValueError(error_message)
    
    if not os.path.exists(args.gene_kos):
        error_message = f"Error: Gene KOs file '{args.gene_kos}' does not exist."
        logging.error(error_message)
        raise ValueError(error_message)

    try:
        solution(args)

    except Exception as e:
        error_message = f"Error: An exception occurred during analysis: {e}"
        logging.error(error_message)
        raise RuntimeError(error_message)


if __name__ == '__main__':
    main()
