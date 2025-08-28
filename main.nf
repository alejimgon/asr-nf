#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
USAGE:

1. Set up and activate the Conda environment
    source scripts/setup.sh

2. Run the pipeline (Stage 1: iTOL visualization):
    nextflow run main.nf
    nextflow run main.nf --tax true --email your_email@example.com   # to fetch taxonomy for leaf IDs and generate LABELS file

3. Visualize the iTOL files and identify interesting nodes.

4. (Optional) Extract FASTA sequences for all or selected nodes (Stage 2):
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
        tuple val(aln_name), path(alignment), path(tree_file), val(model), val(asr_min), val(position), val(aa)

    output:
        tuple val(aln_name), val(model), val(asr_min), path("*.state"), path(alignment), path(tree_file), val(position), val(aa)

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
    publishDir "outputs/${sanitizePrefix(alignment.simpleName, model, asr_min)}/itol_${aa}${position}/", mode: 'copy'

    input:
        tuple path(alignment), path(state_file), val(model), val(asr_min), val(position), val(aa), val(aln_name), path(tree_file)

    output:
        path "${alignment.simpleName}_treecolors.txt"
        path "${alignment.simpleName}_node_prob.txt"
        path "${alignment.simpleName}_leaf_labels.txt", optional: true
        path "${alignment.simpleName}_taxonomy.txt", optional: true

    script:
    def tax_opt = params.tax ? "--taxonomy" : ""
    def email_opt = (params.tax && params.email) ? "--email ${params.email}" : ""
    """
    python ${projectDir}/scripts/extract_aa_to_itol.py ${state_file} ${alignment.simpleName} -al ${alignment} -p ${position} -a ${aa} ${tax_opt} ${email_opt}
    """
}

// --- EXTRACT FASTA FROM STATE ---
process EXTRACT_FASTA_FROM_STATE {
    tag { sanitizePrefix(aln_name, model, asr_min) }
    publishDir { "${params.outputs_dir}/${sanitizePrefix(aln_name, model, asr_min)}" }, mode: 'copy'

    input:
        tuple val(aln_name), val(model), val(asr_min), path(state_file), path(alignment), path(tree_file), val(position), val(aa)

    output:
        path "${sanitizePrefix(aln_name, model, asr_min)}_nodes.fasta"

    script:
    def out_prefix = sanitizePrefix(aln_name, model, asr_min)
    def node_opt = params.fasta_nodes ? "-n ${params.fasta_nodes}" : ""
    def missing_opt = params.fasta_missing ? "-m ${params.fasta_missing}" : ""
    """
    python ${projectDir}/scripts/state_to_fasta.py ${state_file} ${out_prefix}_nodes.fasta ${node_opt} ${missing_opt}
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

    // 3. Run IQTREE_ANCESTRAL for each unique combination
    IQTREE_ANCESTRAL_OUT = IQTREE_ANCESTRAL(ch_full_input)

    // 4. Pass values and run EXTRACT_AA_TO_ITOL
    EXTRACT_AA_TO_ITOL(
        IQTREE_ANCESTRAL_OUT.map { a, b, c, state_file, alignment, tree_file, position, aa ->
            tuple(alignment, state_file, b, c, position, aa, a, tree_file)
        }
    )

    // 5. If needed, run EXTRACT_FASTA_FROM_STATE
    if (params.extract_fasta) {
        EXTRACT_FASTA_FROM_STATE(
            IQTREE_ANCESTRAL_OUT.map { a, b, c, state_file, alignment, tree_file, position, aa ->
                tuple(a, b, c, state_file, alignment, tree_file, position, aa)
            }
        )
    }
}