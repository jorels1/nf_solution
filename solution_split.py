#!/usr/bin/python

import pandas as pd
import logging
import argparse
import os
import sys
import tempfile
from solution_prep import setup_logging



def split_tsv_columns(path, splitsize, outdir, prefix=""):

    """
    This function splits an input TSV file into chunks based on a specified number of columns.

    Parameters:
    - path (str): File path of the TSV to be split.
    - splitsize (int): Number of columns for each chunk.
    - outdir (str): Directory to write output file
    - prefix (str): Prefix for the output chunk filenames.


    Returns:
    - list: A list containing file paths of all chunks.
    """

    logging.info(f'Splitting TSV file into chunks: {path}')
    
    # Read the TSV file into a pandas DataFrame
    df = pd.read_csv(path, sep='\t')
    
    # Identify how many chunks are needed
    n_cols = df.shape[1] - 1

    # Check if the number of columns is less than the splitsize
    if n_cols < splitsize:
        error_message = (f"Error: The input TSV has {n_cols} columns (excluding the first column), "
                         f"which is fewer than the specified splitsize of {splitsize}.")
        logging.error(error_message)
        raise ValueError(error_message)

    chunks = n_cols // splitsize + (n_cols % splitsize > 0)
    
    # Create an empty list for the fles
    output_files = []
    
    # Split the DataFrame into chunks and save each chunk
    for chunk_index in range(chunks):
        start_col = chunk_index * splitsize + 1  # +1 to skip the first column
        end_col = min((chunk_index + 1) * splitsize + 1, df.shape[1])
        
        # Select the required columns for this chunk
        chunk_df = df.iloc[:, [0] + list(range(start_col, end_col))]
        
        # Save the chunk to a file
        try:
            result_tsv_path = os.path.join(outdir, f"{prefix}_chunk_{chunk_index + 1}.tsv")
            chunk_df.to_csv(result_tsv_path, sep='\t', index=False)
            output_files.append(result_tsv_path)

        except Exception as e:
            error_message = f"Error: An exception occurred while saving output files: {e}"
            logging.error(error_message)
            raise RuntimeError(error_message)

    
    logging.info(f"TSV file has been split into {len(output_files)} chunks.")

    return output_files



def solution_split(args):

    """
    This function wraps the split_tsv_columns function to load the mutation data and log the file chunks created.

    Parameters:
    - args (argparse.Namespace): An argparse namespace containing command-line arguments.

    Returns:
    - None
    """

    # Set up logging and create the output directory
    setup_logging(log_level=args.log_level)
    os.makedirs(args.outdir, exist_ok=True)

    # Split the transposed mutations file into chunks
    mutation_chunks = split_tsv_columns(path=args.mutations, splitsize=args.splitsize, outdir=args.outdir, prefix="mutations")
    gene_chunks = split_tsv_columns(path=args.gene_kos, splitsize=args.splitsize, outdir=args.outdir, prefix="genes")

    # Log the chunks created
    for chunk in mutation_chunks:
        logging.info(f"Chunk created at: {chunk}")
    for chunk in gene_chunks:  # Fixed error in this line
        logging.info(f"Chunk created at: {chunk}")



def main():

    """
    Split the transposed mutations and the gene files into specified splitsizes.

    This function parses command-line arguments to specify the input file, splitsize, output settings, and logging options.
    It then calls the 'solution_split()' function to perform the data splitting.

    Command-Line Arguments:
    - -m, --mutations: Mutations TSV file path (required).
    - -g, genes_df (pd.DataFrame): DataFrame containing gene expression data.
    - -c, --splitsize Number of columns for each chunk
    - -l, --log-level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL) (default: INFO).
    - -o, --outdir: Output directory path (default: 'output').  

    Returns:
    - None
    """

    parser = argparse.ArgumentParser(description="Split mutation and gene files")
    parser.add_argument('-m', '--mutations', help='Mutations tsv file', required=True)
    parser.add_argument('-g', '--gene_kos', help='Gene KOs tsv file', required=True)
    parser.add_argument('-c', '--splitsize', help='Number of columns for each chunk. Default: 10', default=10, type=int)
    parser.add_argument('-o', '--outdir', help='Output directory. Default: output', default='output')
    parser.add_argument('-l', '--log-level', help='Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL). \
                        Default: INFO', default='INFO')


    # Parse the command-line arguments passed to the script.
    args = parser.parse_args()


    # Validate the input file path
    if not os.path.exists(args.mutations):
        error_message = f"Error: Mutations file '{args.mutations}' does not exist."
        logging.error(error_message)
        raise ValueError(error_message)
    
    if not os.path.exists(args.gene_kos):
        error_message = f"Error: Gene KOs file '{args.gene_kos}' does not exist."
        logging.error(error_message)
        raise ValueError(error_message)


    # Execute the solution_prep function: in case of error, print and exit with a failure
    try:
        solution_split(args)

    except Exception as e:
        error_message = f"Error: An exception occurred during analysis: {e}"
        logging.error(error_message)
        raise RuntimeError(error_message)


if __name__ == '__main__':
    main()
