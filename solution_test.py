#!/usr/bin/python

import unittest
import os
import shutil
import pandas as pd
import numpy as np
from solution_prep import (setup_logging, get_dataframe_and_transpose)
from solution_split import split_tsv_columns
from solution_run import (get_input_dataframes, generate_gene_and_mutation_lists, 
                            perform_one_way_anova, create_clustered_heatmap)



class TestSolution(unittest.TestCase):



    def test_setup_logging(self):

        """
        This test checks if the setup_logging function returns None.
        """

        # Test setup_logging function
        result = setup_logging('INFO')
        self.assertIsNone(result, "setup_logging did not return None as expected")


    def test_get_dataframe_and_transpose(self):

        """
        This test checks if the transpose_dataframe function correctly transposes a DataFrame and sets column names.
        """

        # Test transpose_dataframe function
        mutations_df = pd.DataFrame({
            'Mutation': [1, 2, 3], 
            'Gene1': [0, 1, 0], 
            'Gene2': [1, 0, 1]})

        transposed_df = get_dataframe_and_transpose(mutations_df)

        self.assertIsInstance(transposed_df, pd.DataFrame, "transposed_df is not a Dataframe.")
        self.assertEqual(transposed_df.index.name, 'Model', "The first column name is not 'Model' and so the index rename did not work")


    def test_split_tsv_columns(self):

        """
        This test checks if the split_tsv_columns function correctly splits a dataframe into chunks and also checks for edge cases
        """

        # Create a temp file for the test
        temp_file_path = 'temp_test_file.tsv'
        df = pd.DataFrame({
            'Model': ['ModelM', 'ModelO', 'ModelS', 'ModelA'],
            'GeneA_KO': [0.0756, -0.1628, -0.1107, -0.4571],
            'GeneB_KO': [-0.8701, -0.3590, -0.9904, -0.1207],
            'GeneC_KO': [-2.0091, -0.8265, -1.1181, -0.2318]
        })

        df.to_csv(temp_file_path, sep='\t', index=False)

        result_files = split_tsv_columns(temp_file_path, 2, ".", prefix="test")

        # Test the functionality of split_tsv_columns()
        # Check the contents of the first output file
        df_chunk1 = pd.read_csv(result_files[0], sep='\t')
        self.assertListEqual(df_chunk1.columns.tolist(), ['Model', 'GeneA_KO', 'GeneB_KO'])
        
        # Check the contents of the second output file
        df_chunk2 = pd.read_csv(result_files[1], sep='\t')
        self.assertListEqual(df_chunk2.columns.tolist(), ['Model', 'GeneC_KO'])

        # Test the edge case: where splitsize is larger than the number of columns
        with self.assertRaises(ValueError) as context:
            split_tsv_columns(temp_file_path, 5, ".", prefix="test")
        self.assertTrue("which is fewer than the specified splitsize" in str(context.exception))

        # Clean up the temp test files
        os.remove(temp_file_path)
        for file in result_files:
            os.remove(file)


    def test_get_input_dataframes(self):

        """
        This test checks if the get_input_dataframes function correctly reads input files and returns two DataFrames.
        """

        # Test get_input_dataframes function
        mutations_file = 'test_mutations.tsv'
        gene_kos_file = 'test_gene_kos.tsv'

        self.assertTrue(os.path.isfile(mutations_file), f"{mutations_file} does not exist.")
        self.assertTrue(os.path.isfile(gene_kos_file), f"{gene_kos_file} does not exist.")

        mutations_df, genes_df = get_input_dataframes(mutations_file, gene_kos_file)

        self.assertIsInstance(mutations_df, pd.DataFrame, "mutations_df is not a DataFrame")
        self.assertIsInstance(genes_df, pd.DataFrame, "genes_df is not a DataFrame")
    
        

    def test_generate_gene_and_mutation_lists(self):

        """
        This test case checks if the function correctly extracts gene and mutation IDs from given dataframes.
        """

        # Create sample dataframes for testing
        genes_df = pd.DataFrame({'Model': ['Model1', 'Model2', 'Model3'], 'Gene1': [0.0756, -0.1628, -0.1107], 'Gene2': [-0.8701, -0.3590, -0.9904], 'Gene3': [-2.0091, -0.8265, -1.1181]})
        mutations_T_df = pd.DataFrame({'Model': ['Model1', 'Model2', 'Model3'], 'Mutation1': [0, 1, 0], 'Mutation2': [1, 0, 1], 'Mutation3': [0, 0, 1]})
        
        # Call the function to generate gene and mutation lists
        gene_list, mutation_list = generate_gene_and_mutation_lists(genes_df, mutations_T_df)

        # Check if the returned values are lists
        self.assertIsInstance(gene_list, list, "gene_list is not a list.")
        self.assertIsInstance(mutation_list, list, "mutation_list is not a list")

        # Check if the extracted lists match the expected values
        self.assertEqual(gene_list, ['Gene1', 'Gene2', 'Gene3'], "Incorrect gene_list")
        self.assertEqual(mutation_list, ['Mutation1', 'Mutation2', 'Mutation3'], "Incorrect mutation_list.")


    def test_perform_one_way_anova(self):

        """
        This function tests the perform_one_way_anova function by providing sample data and comparing the output.
        """

        # Create sample data
        genes_df = pd.DataFrame({'Model': ['Model1', 'Model2', 'Model3'], 'Gene1': [0.0756, -0.1628, -0.1107], 'Gene2': [-0.8701, -0.3590, -0.9904], 'Gene3': [-2.0091, -0.8265, -1.1181]})
        mutations_T_df = pd.DataFrame({'Model': ['Model1', 'Model2', 'Model3'], 'Mutation1': [0, 1, 0], 'Mutation2': [1, 0, 1], 'Mutation3': [0, 0, 1]})
        gene_list = ['Gene1', 'Gene2', 'Gene3']
        mutation_list = ['Mutation1', 'Mutation2', 'Mutation3']
        outdir = 'test_output'

        # Call the function being tested
        output_df = perform_one_way_anova(genes_df, mutations_T_df, gene_list, mutation_list, outdir)

        # Check the output
        self.assertIsInstance(output_df, pd.DataFrame, "output_df is not a DataFrame")
        self.assertEqual(len(output_df), len(gene_list) * len(mutation_list), "Incorrect output_df length")

    def test_create_clustered_heatmap(self):

        """
        This function tests the create_clustered_heatmap function by providing sample data and checking the generated heatmap.
        """

        # Create sample data
        genes = ['Gene1', 'Gene2', 'Gene3']
        mutations = ['Mutation1', 'Mutation2', 'Mutation3']

        # Create all combinations of genes and mutations
        combinations = [(gene, mutation) for gene in genes for mutation in mutations]

        # Create a DataFrame
        data = {
            'gene_id': [combo[0] for combo in combinations],
            'mutation_id': [combo[1] for combo in combinations],
            'stats': np.random.rand(len(combinations)),  # Filling with random numbers
            'P-value': np.random.rand(len(combinations))  # Filling with random numbers
        }

        output_df = pd.DataFrame(data)

        #add the extra settings
        value_column = 'stats'
        index_column = 'gene_id'
        columns_column = 'mutation_id'

        # Call the function being tested
        heatmap_figure = create_clustered_heatmap(output_df, value_column, index_column, columns_column)

        # Check the output
        self.assertIsNotNone(heatmap_figure, "heatmap_figure is None")


if __name__ == '__main__':
    unittest.main()