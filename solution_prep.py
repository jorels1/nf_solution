#!/usr/bin/python

import pandas as pd
import logging
import argparse
import os
import sys



def setup_logging(log_level):

    """
    This function configures the logging system.

    Parameters:
    - log_level (str): The log level to set (DEBUG, INFO, WARNING, ERROR, CRITICAL).

    Returns:
    - None
    """

    # Define a dictionary to map log level strings to their corresponding constants
    log_level_map = {
        'DEBUG': logging.DEBUG,
        'INFO': logging.INFO,
        'WARNING': logging.WARNING,
        'ERROR': logging.ERROR,
        'CRITICAL': logging.CRITICAL
    }

    # Validate the log level argument
    if log_level not in log_level_map:
        error_message = (f"Invalid log level '{log_level}'. Use one of: DEBUG, INFO, WARNING, ERROR, CRITICAL.")
        logging.error(error_message)
        raise ValueError(error_message)

    # Clear existing log handlers to prevent duplicate messages when configuring new ones in setup_logging."
    try:
        del logging.root.handlers[:]
    except:
        pass

    # Configure logging based on the specified log level
    logging.basicConfig(level=log_level_map[log_level], format='[%(levelname)s] %(message)s', stream=sys.stdout)


def get_dataframe_and_transpose(mutations):

    """
    This function reads the mutation data and transpose.

    Parameters:
    - mutations (str): The file path to the mutations data in TSV format.

    Returns:
    - tuple: A tuple containing two pandas dataframes: the mutation and gene knockout data.
    """

    logging.info(f'Loading the input data: {mutations}')

    # Read the mutation data
    if isinstance(mutations, pd.DataFrame):
        mutations_df=mutations
    else:
        mutations_df = pd.read_csv(mutations, sep='\t')

    # Transpose the mutation data
    mutations_T = mutations_df.transpose()

    #make the first row the header
    mutations_T.columns = mutations_T.iloc[0]

    #Drop the first row after assigning as header
    mutations_T_df = mutations_T.drop(mutations_T.index[0])

    # Rename the 'index' column to 'Model'
    mutations_T_df.index.name = 'Model'

    return mutations_T_df


def solution_prep(args):

    """
    This function wraps the above function to tranpose mutation data and generate an output dataframe.

    Parameters:
    - args (argparse.Namespace): An argparse namespace containing command-line arguments.

    Returns:
    - None
    """

    # Set up logging and create the output directory
    setup_logging(log_level=args.log_level)
    os.makedirs(args.outdir, exist_ok=True)

    # Load the input data
    mutations_df = get_dataframe_and_transpose(mutations=args.mutations)

    # Save the dataframe to a TSV file
    try:
        # Save the dataframe to a TSV file
        result_tsv_path = os.path.join(args.outdir, 'mutations_processed.tsv')
        mutations_df.to_csv(result_tsv_path, sep='\t', header=True, index=True)
        logging.info(f"Transposed mutations file saved to {result_tsv_path}")

    except Exception as e:
        error_message = f"Error: An exception occurred while saving output files: {e}"
        logging.error(error_message)
        raise RuntimeError(error_message)


def main():

    """
    Prep the mutation data by transposing it to assume the same orientation as the gene_kos file.

    This function parses command-line arguments to specify the input file, output settings, and logging options.
    It then calls the 'solution_prep()' function to transpose the data and generate results.

    Command-Line Arguments:
    - -m, --mutations: Mutations TSV file path (required).
    - -l, --log-level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL) (default: INFO).
    - -o, --outdir: Output directory path (default: 'output').  

    Returns:
    - None
    """

    parser = argparse.ArgumentParser("Transpose the mutations TSV file")
    parser.add_argument('-m', '--mutations', help='Mutations TSV file', required=True)
    parser.add_argument('-o', '--outdir', help='Output directory. Default: output', default='output')
    parser.add_argument('-l', '--log-level', help='Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL). \
                        Default: INFO', default='INFO')


    # Parse the command-line arguments passed to the script.
    args = parser.parse_args()


    # Validate the input file path
    if not os.path.exists(args.mutations):
        error_message = "Error: Mutations file '{args.mutations}' does not exist."
        logging.error(error_message)
        raise ValueError(error_message)


    # Execute the solution_prep function: in case of error, print and exit with a failure
    try:
        solution_prep(args)

    except Exception as e:
        error_message = f"Error: An exception occurred during analysis: {e}"
        logging.error(error_message)
        raise RuntimeError(error_message)


if __name__ == '__main__':
    main()
