#!/usr/bin/env nextflow

/*
 * A pipeline for automating the analysis of gene-mutations comparison with Nextflow
 *
 * Author:
 * - Joachim Nwezeobi <elochukwujoachim@yahoo.com>
 */

/* 
 * enables modules 
 */
nextflow.enable.dsl=2

/*
 * Set up pipeline parameters. They can be overriden on the command line.
 */

params.mutations_file = params.mutations_file ?: "Mutations.tsv"
params.gene_kos_file = params.gene_kos_file ?: "Gene_KOs.tsv"
params.outdir = params.outdir ?: "output"
params.splitsize = params.splitsize ?: 5
params.split = params.split ?: false
params.value_column = params.value_column ?: "P-value"
params.index_column = params.index_column ?: "mutation_id"
params.columns_column = params.columns_column ?: "gene_id"



// this prevents a warning of undefined parameter
params.help             = false


log.info """\
        mutations_file  : ${params.mutations_file}
        gene_kos_file   : ${params.gene_kos_file}
        outdir          : ${params.outdir}
        splitsize       : ${params.splitsize}
        split           : ${params.split}
        value_column    : ${params.value_column}
        index_column    : ${params.index_column}
        columns_column  : ${params.columns_column}
        """

        .stripIndent()

// this prints the help in case you use --help parameter in the command line and it stops the pipeline
if (params.help) {
    log.info 'This is a Nextflow a pipeline for analysing gene-mutation data'
    log.info 'Please define mutations_file, gene_kos_file and output.\n'
    log.info '\n'
    exit 1
}


// Process to transpose dataset1.csv
process PrepSolution {

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path "mutations_file"

    output:
    path("mutations_processed.tsv"), emit: transposed_mutations_file

    script:
    """
    python /bin/solution_prep.py -m ${mutations_file} -o .

    """
}

process SplitSolution {

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path transposed_mutations_file
    path gene_kos_file

    output:
    path("mutations_chunk_*.tsv"), emit: mutation_chunk
    path("genes_chunk_*.tsv"), emit: gene_chunk

    script:
    """
    python /bin/solution_split.py -m ${transposed_mutations_file} -g ${gene_kos_file} \
                                    -c ${params.splitsize} -o .

    """
}



process RunSolution {

    publishDir "${params.outdir}", mode: 'copy', overwrite: true, pattern: 'gene_vs_mutation_*.{tsv,png}'

    input:
    tuple path(mutation_chunks, stageAs: "dir1/*"), path(gene_chunks, stageAs: "dir2/*")


    output:
    path('gene_vs_mutation_results_*.tsv')
    path('gene_vs_mutation_heatmap_*.png')


    script:
    """
    python /bin/solution_run.py -m ${mutation_chunks} -g ${gene_chunks} -v ${params.value_column} \
                                -i ${params.index_column} -c ${params.columns_column} -o .

    """
}


workflow {

    // Define the channels for the input mutations and genes data sets
    mutations_ch = Channel.fromPath(params.mutations_file, checkIfExists: true)
    genes_ch = Channel.fromPath(params.gene_kos_file, checkIfExists: true)


    // Use PrepSolution to generate the transposed mutations file
    mutation_prep_ch = PrepSolution(mutations_ch)
    .transposed_mutations_file


    // Use the transposed mutations file to split it into chunks if params.split is true
    if (params.split) {
        // Run the SplitSolution process to generate split dataframes depending on the number of columns specified
        splitMutations = SplitSolution(mutation_prep_ch, genes_ch)
        // Create a cartesian product of the different chunks to ensure all possible combination is run
        combinedChunks = splitMutations.mutation_chunk.flatten()combine(splitMutations.gene_chunk.flatten())
    }
    
    else {
        // Simply generate a list of input data as needed by the RunSolution process
        combinedChunks = mutation_prep_ch.combine(genes_ch)
    }


    // Run the RunSolution process
    RunSolution(combinedChunks)

}