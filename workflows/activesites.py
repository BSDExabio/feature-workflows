#!/usr/bin/env python3
""" Given a set of Alphafold2 proteins, find any active sites as identified
by similar proteins found via tmalign and uniprot.
"""
import argparse
import rich

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Utility for finding active sites for proteins processed by alphafold2')
    parser.add_argument('model_quality', help='Model quality pickle file')
    parser.add_argument('alignment_results', help='Alignment results pickle file')
    parser.add_argument('pdb_to_uniprot', help='PDB to UniProt mapping pickle file')
    parser.add_argument('UniProt_metadata', help='UniProt metadata pickle file')

    args = parser.parse_args()
