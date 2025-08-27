#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
USAGE:

1. Create the conda environment (only needs to be done once):
    conda env create -f env.yaml
2. Activate the environment:
    conda activate ancestral_reconstruction
3. Run the pipeline (Stage 1: iTOL visualization):
    nextflow run main.nf
    nextflow run main.nf --tax true --email your_email@example.com   # to fetch taxonomy for leaf IDs and generate LABELS file
4. Visualize the iTOL files and identify interesting nodes.
5. (Optional) Extract FASTA sequences for all or selected nodes (Stage 2):
    nextflow run main.nf --extract_fasta true
    # To extract only specific nodes:
    nextflow run main.nf --extract_fasta true --fasta_nodes Node13,Node25
    nextflow run main.nf --extract_fasta true --fasta_nodes Node13|Node25
    nextflow run main.nf --extract_fasta true --fasta_nodes my_nodes.txt
    # To change the missing character in FASTA output:
    nextflow run main.nf --extract_fasta true --fasta_missing X

Parameters:
    --extract_fasta   Set to true to extract FASTA files for nodes (default: false)
    --fasta_nodes     Node filter: comma-separated list, range (NodeX|NodeY), or file with node names
    --fasta_missing   Character for missing states in FASTA (default: '-')
    --tax             Fetch taxonomy for leaf IDs and generate LABELS file (default: false)
    --email           Email address for NCBI Entrez (required if --tax is used)
*/

params.input_table     = 'data/input_table.txt'
params.alignments_dir  = 'data/alignments/'
params.trees_dir       = 'data/phylogenetic_trees/'
params.outputs_dir     = 'outputs/'
params.email           = null
params.tax             = false
params.extract_fasta   = false
params.fasta_nodes     = null
params.fasta_missing   = '-'

// Sanitize output file names
def sanitizePrefix(aln_name, model, asr_min) {
    return "${aln_name}_${model}_${asr_min}".replaceAll(/[^A-Za-z0-9_]/, "_")
}

// --- ENVIRONMENT CHECK ---
process CHECK_CONDA_ENV {
    tag "check_conda_env"
    output:
        path 'check_env.ok'
    script:
    """
    if [ -z "\$CONDA_DEFAULT_ENV" ]; then
        echo "ERROR: You must run this pipeline inside a conda environment!" >&2
        exit 1
    fi
    python -c "import Bio" || (echo "ERROR: biopython is not installed in this environment!" >&2; exit 1)
    iqtree --version || (echo "ERROR: iqtree is not installed in this environment!" >&2; exit 1)
    touch check_env.ok
    """
}

// --- IQTREE ANCESTRAL RECONSTRUCTION ---
process IQTREE_ANCESTRAL {
    tag { sanitizePrefix(aln_name, model, asr_min) }
    publishDir { "${params.outputs_dir}/${sanitizePrefix(aln_name, model, asr_min)}" }, mode: 'copy'

    input:
        tuple val(aln_name), path(alignment), path(tree_file), val(model), val(asr_min)

    output:
        tuple val(aln_name), val(model), val(asr_min), path("*.state"), path(alignment), path(tree_file)

    script:
    def out_prefix = sanitizePrefix(aln_name, model, asr_min)
    def asr_min_opt = asr_min ? "--asr-min ${asr_min}" : ""
    def model_opt = model ? "-m ${model}" : ""
    """
    iqtree -s ${alignment} -te ${tree_file} --ancestral -asr ${asr_min_opt} ${model_opt} -nt AUTO --prefix ${out_prefix}
    """
}
// --- EXTRACT AA TO ITOL ---
process EXTRACT_AA_TO_ITOL {
    tag { sanitizePrefix(alignment.simpleName, model, asr_min) }

    publishDir "${params.outputs_dir}/itol/${alignment.simpleName}_${model}_${asr_min}/itol_${aa}${position}", mode: 'copy'

    input:
        path extract_aa_to_itol_script
        tuple path(alignment), path(state_file), val(model), val(asr_min), val(position), val(aa)

    output:
        path "${alignment.simpleName}_itol_*", emit: itol_outputs
        path "versions.yml", emit: versions

    script:
    def tax_opt = params.tax ? "--taxonomy" : ""
    def email_opt = (params.tax && params.email) ? "--email ${params.email}" : ""
    """
    python ${extract_aa_to_itol_script} ${state_file} ${params.outputs_dir}/itol/${alignment.simpleName} -al ${alignment} -p ${position} -a ${aa} ${tax_opt} ${email_opt}
    """
}

// --- EXTRACT FASTA FROM STATE ---
process EXTRACT_FASTA_FROM_STATE {
    tag { state_file.getBaseName() + "_" + aa + position }
    publishDir { "${params.outputs_dir}/${state_file.getBaseName()}_${aa}${position}" }, mode: 'copy'

    input:
        tuple path(state_file), path(alignment_path), path(tree_file), val(model), val(position), val(aa), val(asr_min), val(out_prefix), path(state_to_fasta_script)

    output:
        path "*.fasta"

    script:
    def node_opt = params.fasta_nodes ? "-n ${params.fasta_nodes}" : ""
    def missing_opt = params.fasta_missing ? "-m ${params.fasta_missing}" : ""
    """
    python ${state_to_fasta_script} ${state_file} ${out_prefix}_nodes.fasta ${node_opt} ${missing_opt}
    """
}


// --- MAIN WORKFLOW ---
workflow {
    // 1. Check environment
    CHECK_CONDA_ENV()

    // 2. Parse input table and prepare unique IQTREE jobs
    ch_full_input = Channel
    .fromPath(params.input_table)
    .splitCsv(header: false, sep: '\t')
    .map { row ->
        def (alignment, model, position, aa, asr_min) = row*.trim()
        def aln_name = alignment.replace('.aln','')
        def treefile = alignment.replace('.aln', '.treefile')
        tuple(aln_name, file("${params.alignments_dir}/${alignment}"), file("${params.trees_dir}/${treefile}"), model, asr_min, position, aa)
    }

    ch_iqtree_input = ch_full_input
        .map { a, b, c, d, e, f, g -> tuple(a, b, c, d, e) }
        .distinct() // Ensures unique (aln_name, alignment, tree_file, model, asr_min) tuples

    // 3. Run IQTREE_ANCESTRAL for each unique combination
    IQTREE_ANCESTRAL(ch_iqtree_input)

}