#!/usr/bin/env python3
"""
Extract reconstructed sequences for all nodes from an IQ-TREE .state file and write them in FASTA format.

Usage:
    python state_to_fasta.py input.state output.fasta [--missing MISSING_CHAR] [--nodes NODE_FILTER]

Arguments:
    input.state      Input .state file from IQ-TREE
    output.fasta     Output FASTA file

Options:
    -m, --missing MISSING_CHAR  Character to use for missing states (default: -)
    -n, --nodes NODE_FILTER      Node filter: comma-separated list, range 'NodeX|NodeY', or file with node names
"""

import argparse
import csv
from collections import defaultdict

def parse_node_filter(node_filter):
    """Parse the node filter argument and return a set of node names to include."""
    if node_filter is None:
        return None
    # If it's a file
    import os
    if os.path.isfile(node_filter):
        with open(node_filter) as f:
            return set(line.strip() for line in f if line.strip())
    # If it's a range: Node13|Node25
    if '|' in node_filter:
        start, end = node_filter.split('|')
        # Extract prefix and numbers
        import re
        m1 = re.match(r'(\D+)(\d+)', start)
        m2 = re.match(r'(\D+)(\d+)', end)
        if m1 and m2 and m1.group(1) == m2.group(1):
            prefix = m1.group(1)
            n1 = int(m1.group(2))
            n2 = int(m2.group(2))
            return set(f"{prefix}{i}" for i in range(min(n1, n2), max(n1, n2)+1))
    # If it's a comma-separated list
    nodes = [n.strip() for n in node_filter.split(',') if n.strip()]
    return set(nodes)

def parse_ancestral_state_file(state_file, missing_char='-', node_filter=None):
    """Parse the .state file and return a dict of node:sequence, optionally filtering nodes."""
    node_sequences = defaultdict(list)
    with open(state_file, 'r') as f:
        lines = [line for line in f if not line.startswith('#') and line.strip()]
        reader = csv.DictReader(lines, delimiter='\t')
        for row in reader:
            node = row['Node']
            if node_filter is not None and node not in node_filter:
                continue
            state = row['State']
            node_sequences[node].append(state if state != '-' else missing_char)
    # Join the states to form sequences
    return {node: ''.join(seq) for node, seq in node_sequences.items()}

def write_fasta(sequences, output_file):
    """Write node sequences to a FASTA file."""
    with open(output_file, 'w') as f:
        for node, seq in sequences.items():
            f.write(f">{node}\n{seq}\n")

def main():
    parser = argparse.ArgumentParser(description="Extract reconstructed sequences from IQ-TREE .state file and write FASTA.")
    parser.add_argument("input_file", help="Input .state file from IQ-TREE")
    parser.add_argument("output_file", help="Output FASTA file")
    parser.add_argument("-m", "--missing", default="-", help="Character for missing states (default: -)")
    parser.add_argument("-n", "--nodes", default=None, help="Node filter: comma-separated list, range 'NodeX|NodeY', or file with node names")
    args = parser.parse_args()

    node_filter = parse_node_filter(args.nodes) if args.nodes else None
    sequences = parse_ancestral_state_file(args.input_file, missing_char=args.missing, node_filter=node_filter)
    write_fasta(sequences, args.output_file)

if __name__ == "__main__":
    main()