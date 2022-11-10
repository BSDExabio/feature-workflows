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
from pathlib import Path
import re

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


def get_ids(af_candidates, proteinID_aln_results, PDBID_to_UniProt_map):
    """ Look up the PDB IDs for all the AF candidate proteins

        :param af_candidates: dict of alphafold proteins
        :param proteinID_aln_results: dict of proteins containing PDB IDs
        :param PDBID_to_UniProt_map: to look up UniProt ID for PDB IDs
        :returns: updated af_candidates dict with ['pdbids'] a list of PDB Is
    """
    total_skipped = 0 # total skipped because no UniProt ID could be found
    for protein in af_candidates.keys():
        af_candidates[protein]['pdbids'] = {}
        af_candidates[protein]['maxTMscore'] = proteinID_aln_results[protein]['maxTMscore']
        for pdbid_path in proteinID_aln_results[protein]['Target Path']:
            # Target Path is just a capture of the fully qualified path to a
            # PDB file for a protein.  So we use path string manipulation to
            # extract the PDB ID, which is just the path stem.
            pdbid = Path(pdbid_path).stem
            uniprotid = PDBID_to_UniProt_map.get(pdbid, None)
            if uniprotid is None:
                logger.warning(f'{pdbid} has no corresponding UniProt ID ...skipping')
                total_skipped += 1
                continue
            af_candidates[protein]['pdbids'][pdbid] = {'ecIDs'   : set([]),
                                                       'UniProt' : uniprotid,
                                                       'ACT_SITE': None,
                                                       'BINDING' : None}
    if total_skipped > 0:
        logger.warning(f'Skipped {total_skipped} PDB IDs because no corresponding UniProt ID could be found')

    return toolz.valfilter(lambda x: x['pdbids'] != {}, af_candidates)


def extract_uniprot_info(af_pdbids, UniProt_metadata_dict):
    """ Look up interesting UniProt information by PDB ID

    :param af_pdbids: from which to cross-reference
    :param UniProt_metadata_dict: of UniProt in which to look up IDs
    :return: interesting UniProt info
    """
    for af_value in af_pdbids.values():
        for pdb_value in af_value['pdbids'].values():
            uniprot = pdb_value['UniProt']
            if UniProt_metadata_dict[uniprot]['ecIDs'] == []:
                # Skip any UniProt proteins that are not enzymes
                continue
            # Add all the enzymes to 'ecIDs' list
            for enzyme in UniProt_metadata_dict[uniprot]['ecIDs']:
                # EC code three dot separated numbers and then a space, dash,
                # number, or some other character
                found_enzyme = re.search('EC=(\d+.\d+.\d+.[ \-\w+])', enzyme)
                if found_enzyme: # extract EC code
                    pdb_value['ecIDs'].add(found_enzyme.group(1))

    return af_pdbids



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

    logger.info(f'Have {len(af_candidates)} AF candidates out of '
                f'{len(proteinID_model_qualities)} proteins')

    # Next we grab all the corresponding PDB IDs for those candidates to later
    # use for looking up info in UniProt.
    af_with_pdbids = get_ids(af_candidates, proteinID_aln_results, PDBID_to_UniProt_map)

    # Now we extract all the goodness from UniProt for each PDB ID
    results = extract_uniprot_info(af_with_pdbids, UniProt_metadata_dict)

    logger.info('Done')
