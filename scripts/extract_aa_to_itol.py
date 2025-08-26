#!/usr/bin/env python3
"""
This script extract the predicted amino acid at a given position for each node from an IQ-TREE sequence ancestral reconstruction file (.state file).
It also extract the amino acid at the same position for each leaf from a multiple sequence alignment (MSA) file (FASTA format).
Using the extracted information, it generates several iTOL-compatible files for visualization.

Usage:

    python extract_aa_to_itol.py input.state output [prefix] -al alignment.fasta -p [NUM] -a [AA]

Optional arguments:

    -tax, --taxonomy   Fetch taxonomy for leaf IDs and generate LABELS file (requires Biopython)
    -e, --email        Email address for NCBI Entrez (required if -tax is used)

Requirements:

    - Biopython

"""

import argparse
import csv
import os
from Bio import Entrez, SeqIO
import re

parser = argparse.ArgumentParser(description="Extracts the amino acid at a given position for each node from an IQ-TREE .state file and generates an iTOL binary dataset file.")

parser.add_argument("input_file", help="Input .state file from IQ-TREE")
parser.add_argument("output", help="Output prefix for iTOL files")
parser.add_argument("-p", "--position", type=int, help="1-based position to check")
parser.add_argument("-a", "--aa", type=str, help="Amino acid of interest")
parser.add_argument("-al", "--alignment", type=str, help="Alignment file (FASTA format) for leaf label coloring")
parser.add_argument("-tax", action="store_true", help="Fetch taxonomy for leaf IDs and generate LABELS file (requires Biopython)")
parser.add_argument("-e", "--email", type=str, help="Email address for NCBI Entrez (required if -tax is used)")

args = parser.parse_args()

input_file = args.input_file
output_file = args.output
position = args.position
aa_of_interest = args.aa.upper()

# Color for the amino acid of interest and for others
color_aa = "#1f77b4"   # Blue for aa of interest
color_other = "#d62728" # Red for others

nodes = []

def branch_colors(nodes, output_file, aa_of_interest, color_aa, color_other):
    """Write a file with branch colors based on amino acid states."""

    with open(output_file, "w") as out:
        out.write("TREE_COLORS\n")
        out.write("SEPARATOR TAB\n")
        out.write("DATA\n")
        for node, state, pAA in nodes:
            if state == aa_of_interest:
                color = color_aa
            elif state:
                color = color_other
            else:
                color = "#cccccc"
            out.write(f"{node}\tbranch\t{color}\tnormal\t1\n")


def node_prob_aa(nodes, output_file):
    """Write a file with metadata: node probability and amino acid."""

    with open(output_file, "w") as out:
        out.write("METADATA\n")
        out.write("SEPARATOR COMMA\n")
        out.write("FIELD_LABELS,probability,AA\n")
        out.write("DATA\n")
        for node, state, pAA in nodes:
            out.write(f"{node},{pAA},{state}\n")


def leaf_label_colors_from_alignment(alignment_file, output_file, position, aa_of_interest, color_aa, color_other):
    """Read a FASTA alignment, extract AA at position for each leaf, and write TREE_COLORS label lines."""

    with open(output_file, "w") as out:
        out.write("TREE_COLORS\n")
        out.write("SEPARATOR TAB\n")
        out.write("DATA\n")
        with open(alignment_file) as f:
            for record in SeqIO.parse(f, "fasta"):
                seq = str(record.seq)
                if len(seq) >= position:
                    aa = seq[position-1].upper()
                    color = color_aa if aa == aa_of_interest else color_other
                    out.write(f"{record.id}\tbranch\t{color}\tnormal\t1\n")
                    out.write(f"{record.id}\tlabel\t{color}\tnormal\t1\n")


def taxonomy_label(alignment_file, output_file, email):
    """For each accession in the alignment FASTA, fetch the organism name from NCBI and write an iTOL LABELS file."""

    Entrez.email = email
    accessions = []
    with open(alignment_file) as f:
        for record in SeqIO.parse(f, "fasta"):
            accessions.append(record.id)
    with open(output_file, "w") as out:
        out.write("LABELS\n")
        out.write("SEPARATOR TAB\n")
        out.write("DATA\n")
        for acc in accessions:
            try:
                handle = Entrez.efetch(db="protein", id=acc, rettype="gb", retmode="text")
                rec = SeqIO.read(handle, "genbank")
                # Try to get organism from source feature
                species = None
                for feature in rec.features:
                    if feature.type == "source":
                        species = feature.qualifiers.get("organism", [None])[0]
                        break
                # Fallback to description [Organism Name]
                if not species and rec.description:
                    m = re.search(r'\[([^\[\]]+)\]$', rec.description)
                    if m:
                        species = m.group(1)
                # Fallback to annotations
                if not species:
                    species = rec.annotations.get("organism", None)
                out.write(f"{acc}\t{acc}_{species}\n")
            except Exception as e:
                out.write(f"{acc}\tError: {e}\n")


with open(input_file) as f:
    # Skip comment lines
    lines = [line for line in f if not line.startswith("#")]
    reader = csv.DictReader(lines, delimiter="\t")
    for row in reader:
        if int(row["Site"]) == position:
            node = row["Node"]
            state = row["State"]
            # Get probability of amino acid in each node
            pAA = 0.0
            for key in row:
                if key.startswith("p_") and key.endswith(state):
                    pAA = float(row[key])
            nodes.append((node, state, pAA))

# Write TREE_COLORS branch coloring file
treecolors_file = os.path.splitext(output_file)[0] + "_treecolors.txt"
branch_colors(nodes, treecolors_file, aa_of_interest, color_aa, color_other)

# Write TREE_COLORS branch coloring file with probabilities
node_prob_file = os.path.splitext(output_file)[0] + "_node_prob.txt"
node_prob_aa(nodes, node_prob_file)

# Write TREE_COLORS leaf label coloring file from alignment (if alignment exists)
if args.alignment:
    leaf_label_file = os.path.splitext(output_file)[0] + "_leaf_labels.txt"
    leaf_label_colors_from_alignment(args.alignment, leaf_label_file, position, aa_of_interest, color_aa, color_other)
    if args.tax:
        if not args.email:
            raise ValueError("Email address is required for taxonomy fetching with -tax option.")
        taxonomy_label_file = os.path.splitext(output_file)[0] + "_taxonomy.txt"
        taxonomy_label(args.alignment, taxonomy_label_file, args.email)
