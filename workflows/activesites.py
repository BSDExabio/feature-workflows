#!/usr/bin/env python3
""" Given a set of Alphafold2 proteins, find any active sites as identified
by similar proteins found via tmalign and uniprot.

Exits with a 0 if ok, 1 if the model quality file is wrong, 2 if the alignment
results files is wrong, and 3 if the UniProp metadata file is wrong.
"""
import argparse
import pickle
import logging
import sys

import toolz
import pandas as pd

from rich.logging import RichHandler

# Create unique logger for this namespace
rich_handler = RichHandler(rich_tracebacks=True,
                           markup=True)
logging.basicConfig(level='INFO', format='%(message)s',
                    datefmt="[%Y/%m/%d %H:%M:%S]",
                    handlers=[rich_handler])
logger = logging.getLogger(__name__)


def load_pickle_file(filename):
    logger.info(f'Loading {filename}')
    with open(filename, 'rb') as pickle_file:
        return pickle.load(pickle_file)


def validate_input_files(proteinID_model_qualities,
                         proteinID_aln_results,
                         PDBID_to_UniProt_map,
                         UniProt_metadata_dict):
    """ Ensure that the *right* pickle files have been loaded.

    We do that by checking for certain data structure features by file.

    We will immediately exit if we don't find certain features

    :param proteinID_model_qualities:
    :param proteinID_aln_results:
    :param PDBID_to_UniProt_map:
    :param UniProt_metadata_dict:
    :return: None
    """
    # first check the alphafold model qualities dict
    key = toolz.first(proteinID_model_qualities.keys())

    if not set(['model_name', 'ptms', 'plddts', 'iterations']) == \
           set(proteinID_model_qualities[key].keys()):
        logger.critical(f'Invalid model quality file')
        sys.exit(1)

    # second, check the alignment results
    key = toolz.first(proteinID_aln_results.keys())
    if not set(['Len1', 'Len2', 'RMSD', 'TMscore1', 'TMscore2']) <= \
           set(proteinID_aln_results[key].keys()):
        logger.critical(f'Invalid alignment results file')
        sys.exit(2)

    # third, check the PDBID to UniProt mapping dict
    if len(PDBID_to_UniProt_map.keys()) < 84000:
        logger.critical(f'Invalid PDB to UniProt mapping file')
        sys.exit(3)

    # finally, check the UniProt metadata
    key = toolz.first(UniProt_metadata_dict.keys())
    if not set(['entry_name', 'status', 'features', 'sequence']) <= \
           set(UniProt_metadata_dict[key].keys()):
        logger.critical(f'Invalid UniProp metadata file')
        sys.exit(4)

    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Utility for finding active '
                                                 'sites for proteins processed '
                                                 'by alphafold2')

    parser.add_argument('model_quality', help='Model quality pickle file')
    parser.add_argument('alignment_results',
                        help='Alignment results pickle file')
    parser.add_argument('pdb_to_uniprot',
                        help='PDB to UniProt mapping pickle file')
    parser.add_argument('uniprot_metadata', help='UniProt metadata pickle file')

    args = parser.parse_args()

    proteinID_model_qualities = load_pickle_file(args.model_quality)
    proteinID_aln_results = load_pickle_file(args.alignment_results)
    PDBID_to_UniProt_map = load_pickle_file(args.pdb_to_uniprot)
    UniProt_metadata_dict = load_pickle_file(args.uniprot_metadata)

    # Ensure we got the right data loaded from those pickle files by doing
    # things like checking for dictionary keys and the number of rows. I.e., we
    # don't want to get too far into this having one or more of these be wrong.
    validate_input_files(proteinID_model_qualities,
                         proteinID_aln_results,
                         PDBID_to_UniProt_map,
                         UniProt_metadata_dict)

    # Find all of the alphafold proteins that have high quality predicted
    # confidence as denoted by ptms scores > .7
    # TODO Make .7 an optional command line argument
    af_candidates = toolz.valfilter(lambda x: x['ptms'] > 0.7,
                                    proteinID_model_qualities)


    logger.info('Done')