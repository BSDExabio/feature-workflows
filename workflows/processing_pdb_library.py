
import argparse
import time
import platform
from pathlib import Path
import logging
import csv
import pandas
import pickle
import MDAnalysis

import dask
import dask.config
from distributed import Client, Worker, as_completed, get_worker

# codes housed in structural_DLFA repo
import rcsb_query
import uniprot_query

#######################################
### LOGGING FUNCTIONS
#######################################

def append_timings(csv_writer, file_object, hostname, worker_id, start_time, stop_time,
                   query, task_type, return_code):
    """ append the task timings to the CSV timings file
    :param csv_writer: CSV to which to append timings
    :param file_object: file object associated with the CSV
    :param hostname: on which the processing took place
    :param worker_id: of the dask worker that did the processing
    :param start_time: start time in *NIX epoch seconds
    :param stop_time: stop time in same units
    :param query: query used as input to a task
    :param task_type: integer used to denote which step of the workflow has been performed
    :param return_code: result of the task; 1=successful, 0=failure
    """
    csv_writer.writerow({'hostname'   : hostname,
                         'worker_id'  : worker_id,
                         'start_time' : start_time,
                         'stop_time'  : stop_time,
                         'query'      : query,
                         'task_type'  : task_type,
                         'return_code': return_code})
    file_object.flush()


def setup_logger(name, log_file, level=logging.INFO):
    """To setup as many loggers as you want"""
    formatter = logging.Formatter('%(asctime)s    %(levelname)s       %(message)s')
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger


def clean_logger(logger):
    """To cleanup the logger instances once we are done with them"""
    for handle in logger.handlers:
        handle.flush()
        handle.close()
        logger.removeHandler(handle)


#######################################
### DASK RELATED FUNCTIONS
#######################################

def get_num_workers(client):
    """ Get the number of active workers
    :param client: active dask client
    :return: the number of workers registered to the scheduler
    """
    scheduler_info = client.scheduler_info()

    return len(scheduler_info['workers'].keys())


#######################################
### SUB-FUNCTIONS
#######################################

def get_pdbid_chainid(pdb_file):
    """Gather the PDBID_CHAINID string from a pdb_file
    NOTE: this function might need to change based on where the PDBID_CHAINID
    information is stored within the pdb_file.
    :param pdb_file: path object that points to a pdb file, which has an 
                     associated PDBID_CHAINID.
    :return: string of format AAAA_B where AAAA is the four character RCSB PDB
             ID and B is the chain ID of interest
    
    For the PDB70 structure set, pdb files have the PBDID_CHAINID as the stem.
    For the SCOPe structure set, PDBID_CHAINID are stored in the REMARK lines.
    """
    # FOR PARSING PDB70 structures
    return pdb_file.stem
    # FOR PARSING SCOPe structures
    #return ....

#######################################
### TASK FUNCTIONS
#######################################
# Since tasks 0 to 2 are gonna be strung together and parsed in the same 
# `as_completed' for-loop, we need their function's output to be 
# similarly formatted.
# Output format: 
#       task_num: a number to denote which task in the workflow has been 
#                 completed.
#       results: a dictionary of a dictionary. The inner dictionary can be 
#                filled with any key:value pairs; the sole key for the outer
#                dictionary should be the relevant ID that was parsed.
#       hostname: a string that describes the platform/compute node the task
#                 was performed on.
#       start,stop: epoch times for when the task began and finished.
#       return_code: integer value; 0 if the job successfully completed all
#                    necessary tasks; !=0 if the task errored out.

# step0
def _parse_pdb_file(pdb_file):
    """Load the structure file to gather sequence information of the structure
    :param pdb_file: string or path object, points to PDB file that has some form of PDBID_CHAINID identifier. The get_pdbid_chainid() is used to gather that information. 
    :return: step_number, subdictionary, host_information, worker_ID, start_time, stop_time, return_code

    NOTE: this function assumes only a single chain is present in the structure file; will need to rework this is applied to a multi-chained structure
    """
    start_time = time.time()
    worker = get_worker()
    pdb_dict = {}
    pdb_dict['structure_path'] = str(pdb_file)
    pdbid_chainid = get_pdbid_chainid(pdb_file)
    return_code = -9999
    try:
        u = MDAnalysis.Universe(str(pdb_file))
        sel = u.select_atoms('protein')
        start_resid = sel.residues[-1].resid
        end_resid   = sel.residues[-1].resid
        # need to check for non-standard aa resnames that MDAnalysis cannot parse into a sequence.
        # replace these nonstandard resnames with their closest analog 
        # MSE -> MET, PYL -> LYS, SEC -> CYS; list may be incomplete and so this may fail for rare cases
        # standard residues that MDAnalysis comprehends: https://userguide.mdanalysis.org/1.1.1/standard_selections.html
        for resid in range(sel.n_residues):
            if sel.residues[resid].resname == 'MSE':
                sel.residues[resid].resname = 'MET'
            elif sel.residues[resid].resname == 'PYL':
                sel.residues[resid].resname = 'LYS'
            elif sel.residues[resid].resname == 'SEC':
                sel.residues[resid].resname = 'CYS'
        seq = sel.residues.sequence().seq
        pdb_dict['struct_first_resid'] = sel.residues[0].resid
        pdb_dict['struct_last_resid']  = sel.residues[-1].resid
        pdb_dict['struct_seq']  = sel.residues.sequence(format='string')
        return_code = 0
    except Exception as e:
        print(f'failed to load {pdb_file} into MDAnalysis and gather residue indices and sequence. Exception: {e}', file=sys.stdout, flush=True)
    
    stop_time = time.time()
    return 0, {pdbid_chainid: pdb_dict}, platform.node(), worker.id, start_time, stop_time, return_code

# step 1
def _query_rcsb(step0_results):
    """Gather UniProt ID of a pdbID-chainID string using RCSB graphql request
    :param step0_results: tuple, element 1 of this tuple is the results dictionary output from step0, with the sole key should have the expected format AAAA_B (e.g. 2JLR_A) where AAAA is a PDB ID and B is the chain. 
    :return: step_number, subdictionary, host_information, worker_ID, start_time, stop_time, return_code
    """
    start_time = time.time()
    worker = get_worker()
    uniprotid = None
    return_code = -9999
    try:
        uniprotid = rcsb_query.query_uniprot_str(pdbid_chainid)
        return_code = 0
    except Exception as e:
        print(f'failed to pull the uniprot accession id associated with {pdbid_chainid}. Exception: {e}', file=sys.stdout, flush=True)

    stop_time = time.time()
    return 1, {pdbid_chainid: uniprotid}, platform.node(), worker.id, start_time, stop_time, return_code


# step 2
def _query_uniprot_flat_file(uniprotid):
    """Gather UniProt flat file metadata associated with a uniprotid
    :param uniprotid: string, numerous formats; expected to be associated with an entry in the UniProtKB database
    :return: step_number, subdictionary, host_information, worker_ID, start_time, stop_time, return_code
    """
    start_time = time.time()
    worker = get_worker()
    meta_dict = {}
    return_code = -9999
    try:
        meta_dict = uniprot_query.request_uniprot_metadata(uniprotid)
        return_code = 0
    except Exception as e:
        print(f'failed to pull {uniprotid} flat file. Exception: {e}', file=sys.stdout, flush=True)

    stop_time = time.time()
    return 2, {uniprotid: meta_dict}, platform.node(), worker.id, start_time, stop_time, return_code


# step 3
def _blast_aln(structure_dict,uniprotid_dict,blastp_path,ids):
    """Run and parse a blastp sequence alignment between the structure and flat file sequences
    """
    start_time = time.time()
    worker = get_worker()
    seq_aln_dict = {}
    return_code = -9999
    try:
        # unfortunately need to write fasta files for blastp to run
        # query will always be the PDBID_CHAINID structure's sequence
        # subject will always be the uniprot flat file's sequence
        with open(f'{worker.id}_query.fasta','w') as fsta:
            fsta.write(f">query\n{structure_dict['struct_seq']}")
        with open(f'{worker.id}_sbjct.fasta','w') as fsta:
            fsta.write(f"{uniprotid_dict['sequence']}")

        xml_output = NcbiblastpCommandline(query=f'{worker.id}_query.fasta',subject='{worker.id}_sbjct.fasta',cmd=blastp_path,outfmt=5)()[0]
        blast_record = NCBIXML.read(StringIO(xml_output))
        seq_aln_dict['query_start'] = blast_record.alignments[0].hsps[0].query_start
        seq_aln_dict['query_end']   = blast_record.alignments[0].hsps[0].query_end
        seq_aln_dict['sbjct_start'] = blast_record.alignments[0].hsps[0].sbjct_start
        seq_aln_dict['sbjct_end']   = blast_record.alignments[0].hsps[0].sbjct_end
        seq_aln_dict['alignment']   = list(zip(blast_record.alignments[0].hsps[0].query,blast_record.alignments[0].hsps[0].match,blast_record.alignments[0].hsps[0].sbjct))
        seq_aln_dict['e-value']     = blast_record.alignments[0].hsps[0].expect
        seq_aln_dict['score']       = blast_record.alignments[0].hsps[0].score
        seq_aln_dict['bits']        = blast_record.alignments[0].hsps[0].bits
        seq_aln_dict['n_gaps']      = blast_record.alignments[0].hsps[0].gaps
        return_code = 0        
        
    except Exception as e:
        print(f'failed to run the blastp alignment. Exception: {e}')

    stop_time = time.time()
    return 3, seq_aln_dict, ids, platform.node(), worker.id, start_time, stop_time, return_code


#######################################
### MAIN
#######################################

if __name__ == '__main__':
    # read command line arguments.
    parser = argparse.ArgumentParser(description='Parsing and requesting metadata associated with structural alignment hits.')
    parser.add_argument('--scheduler-file', '-s', required=True, help='dask scheduler file')
    parser.add_argument('--input-list-file', '-inp', required=True, help='list file that contains the paths to pdb files. Stem of the file paths is assumed to be a PDBID_CHAINID format.')
    parser.add_argument('--timings-file', '-ts', required=True, help='CSV file for protein processing timings')
    parser.add_argument('--tskmgr-log-file', '-log', required=True, help='string that will be used to store logging info for this run')
    args = parser.parse_args()

    # start dask client.
    client = Client(scheduler_file=args.scheduler_file,timeout=5000,name='AlignmentTaskMgr')

    # set up the main logger file and list all relevant parameters.
    main_logger = setup_logger('tskmgr_logger', args.tskmgr_log_file)
    main_logger.info(f'Starting dask pipeline and setting up logging. Time: {time.time()}')
    main_logger.info(f'Scheduler file: {args.scheduler_file}')
    main_logger.info(f'Timing file: {args.timings_file}')
    dask_parameter_string = ''
    for key, value in dask.config.config.items():
        dask_parameter_string += f"'{key}': '{value}'\n"
    dask_parameter_string += f'Client information: {client}\n'
    dask_parameter_string += '################################################################################'
    main_logger.info(f'Dask parameters:\n{dask_parameter_string}')
   
    with open(args.input_list_file,'r') as list_file:
        pdbid_chainid_files = [Path(line.strip()) for line in list_file.readlines()]
    main_logger.info(f'{len(pdbid_chainid_files)} PDBID_CHAINID files will be parsed.')

    # set up timing log file.
    timings_file = open( args.timings_file, 'w')
    timings_csv = csv.DictWriter(timings_file,['hostname','worker_id','start_time','stop_time','query','task_type','return_code'])
    timings_csv.writeheader()
   
    # dictionary within which the structural information associated with each PDBID_CHAINID will be stored
    # if a pdbid_chainid structure file fails to be read then an incomplete dictionary will be left instead; this incomplete dictionary will only contain "structure_path" key; missing "struct_seq" and numerous other important keys
    pdbid_chainid_metadata_dict = {}
    
    # list within which pdbid_chainid strings will be stored that were successfully parsed.
    pdbid_chainid_seqs = []

    # dictionary within which the pdbid_chainid key will map to the UniProt Accession ID value
    # if a pdbid_chainid does not have an associated UniprotID then a None object will be left instead
    pdbid_to_uniprot_dict= {}
    
    # list within which uniprotid strings will be stored that were successfully queried.
    uniprotid_list = []
    
    # dictionary within which the metadata associated with UniprotIDs (keys) will be stored. 
    # will be a dictionary of dictionaries, because there are many fields/types of information
    # in the Uniprot flat files. 
    uniprot_metadata_dict= {}

    ### do the thing; step0 and step1 are chained together
    step0_futures = client.map(_parse_pdb_file, pdbid_chainid_files)
    step1_futures = client.map(_query_rcsb, step0_futures)
    futures_bucket = step0_futures + step1_futures
    ac = as_completed(futures_bucket)
    for finished_task in ac:
        task_num, results, hostname, workerid, start, stop, return_code = finished_task.result()
        # step0 query_string -> pdbid_chainid
        # step1 query_string -> pdbid_chainid
        # step2 query_string -> uniprotid
        query_string = list(results.keys())[0]
        append_timings(timings_csv,timings_file,hostname,workerid,start,stop,query_string,task_num,return_code)
        # handling step 0 results:
        if task_num == 0:
            pdbid_chainid_metadata_dict.update(results)
            if return_code == 0:
                main_logger.info(f'The pdb structure file associated with {query_string} has been parsed. Return code: {return_code}. Took {stop-start} seconds.')
                pdbid_chainid_seqs.append(query_string)
            else:
                main_logger.info(f'The structure associated with {query_string} failed to be parsed. Still passing PDBID_CHAINID onto step 1 but no step 3 will be run. Return code: {return_code}. Took {stop-start} seconds.')

        # handling step 1 results:
        if task_num == 1:
            pdbid_to_uniprot_dict.update(results)
            if return_code == 0 and results[query_string] and results[query_string] not in uniprotid_list:
                main_logger.info(f'The UniProt accession ID associated with {query_string} has been queried. Return code: {return_code}. Took {stop-start} seconds.')
                pdbid_chainid_metadata_dict[query_string]['uniprotid_aln'] = {results[query_string] : None}
                uniprotid_list.append(results[query_string])
                ac.add(client.submit(_query_uniprot_flat_file,results[query_string]))
            elif return_code == 0 and results[query_string]:
                main_logger.info(f'The UniProt accession ID associated with {query_string} has been queried. The UniprotID has already been seen. Return code: {return_code}. Took {stop-start} seconds.')
                pdbid_chainid_metadata_dict[query_string]['uniprotid_aln'] = {results[query_string] : None}
            else:
                main_logger.info(f'Failed to query {query_string} to gather its UniProtID. Return code: {return_code}. Took {stop-start} seconds.')
        
        # handling step 2 results:
        elif task_num == 2:
            uniprot_metadata_dict.update(results)
            if return_code == 0:
                main_logger.info(f'The flat file associated with {query_string} has been parsed. Return code: {return_code}. Took {stop-start} seconds.')
            else:
                main_logger.info(f'The flat file associated with {query_string} failed to be parsed. Return code: {return_code}. Took {stop-start} seconds.')

    main_logger.info(f'Done parsing PDBID_CHAINID files, querying RCSB, and parsing UniProt flat files. Writing map and uniprot metadata dictionaries to file. {time.time()}')

    # save dictionary of the pdbid_chainid to uniprotid mapping
    with open('pdbid_to_uniprotid_map.pkl', 'wb') as out:
        pickle.dump(pdbid_to_uniprot_dict,out)

    # save dictionary of uniprot accession id meta data
    with open('uniprot_metadata.pkl', 'wb') as out:
        pickle.dump(uniprot_metadata_dict,out)

    main_logger.info(f'Beginning to run blastp alignments between structure and uniprot sequences to get residue index mapping. {time.time()}')
    task3_futures = []
    for pdbid_chainid in pdbid_to_uniprot_dict.keys():
        uniprotid = pdbid_to_uniprot_dict[pdbid_chainid]
        # check to see that the pdbid_chainid maps to a real uniprotid
        # and
        # check to see if the 'struct_seq' key is present in the pdbid_chainid_metadata_dict subdictionary
        if not uniprotid and 'struct_seq' in pdbid_chainid_metadata_dict[pdbid_chainid].keys():
            task3_futures.append(client.submit(_blast_aln,pdbid_chainid_metadata_dict[pdbid_chainid],uniprot_metadata_dict[uniprotid], args.blastp_path, [pdbid_chainid,uniprotid]))

    ac = as_completed(task3_futures)
    for finished_task in ac:
        task_num, results, ids, hostname, workerid, start, stop, return_code = finished_task.result()
        pdbid_chainid = ids[0]
        uniprotid     = ids[1]
        pdbid_chainid_metadata_dict[pdbid_chainid]['uniprotID_aln'][uniprotid] = results

    # close log files and shut down the cluster.
    timings_file.close()
    main_logger.info(f'Done. Shutting down the cluster. Time: {time.time()}')
    clean_logger(main_logger)

