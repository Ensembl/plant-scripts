#!/usr/bin/env python3 

# Python3 script to run RepeatDetector (Red v2) and optionally 
# feed results into an Ensembl core database
#
# Tested with: 
# pyenv local 3.7.6
# pip install --user sqlalchemy_utils pymysql
# 
# Bruno Contreras Moreira, Carlos García Girón EBI-EMBL 2020
#
# See https://www.ebi.ac.uk/seqdb/confluence/display/EnsGen/Repeat+masking+with+Red

import argparse
import os
import re
import errno
import subprocess

import sqlalchemy as db
import sqlalchemy_utils as db_utils

# sqlalchemy requires MySQLdb but MySQLdb doesn't support Python 3.x
# pymysql can be imported and used instead
import pymysql
pymysql.install_as_MySQLdb()


def parse_FASTA_sequences( genome_file , dirname ):
    '''Takes FASTA genome file name, parses individual sequences and 
       saves them in multiple files in named directory.
       Returns: list of successfully parsed sequence names'''

    seq_names = []
    try:
        file = open(genome_file)
    except OSError as error:
        print("# ERROR: cannot open/read file:", genome_file, error)
        return num_seqs

    seq_filepath = ''
    prev_filepath = ''
    for line in file:
       header = re.search(r'^>', line) 
       if header:
           # check previous file was open
           prev_filepath = seq_filepath
           if prev_filepath:
               seqfile.close()

           # open temp FASTA file for this sequence only
           seq_name_match = re.search(r'^>(\S+)', line)
           if seq_name_match:
               seq_name = seq_name_match.group(1)
               seq_filename = seq_name + '.fa'
               seq_filepath = os.path.join(dirname, seq_filename)

               try:
                   seqfile = open(seq_filepath,"w")
               except OSError as error:
                   print("# ERROR: cannot create file ", seq_filepath, error)

               seqfile.write(">%s\n" % seq_name)
               seq_names.append(seq_name)
           else:
               print("# ERROR: cannot parse FASTA header:", header)
       else:
           if seqfile:
               seqfile.write(line)
    
    if seqfile: seqfile.close()
    file.close()

    return seq_names


def parse_params_from_log( log_filename ):
    '''Parses Red stdout log and returns two strings: 
       i) Red version ii) Parameters of this job'''

    version = 'NA'
    params = ''
    try:
        logfile = open(log_filename)
    except OSError as error:
        print("# ERROR: cannot open/read file:", log_filename, error)
        return version

    for line in logfile:
        versionre = re.search(r'^Version: (\S+)', line)
        if versionre:
            version = versionre.group(1)
        else:
            paramre = re.search(r'^(-\w+: \S+)', line)
            if paramre:
                params = params + ' ' + paramre.group(1)

    logfile.close

    return version, params


def _parse_rptfiles_from_log( log_filename ):
    '''Parses Red stdout log and returns a list with
       the names of output files with repeat coordinates'''

    rpt_files = []

    try:
        logfile = open(log_filename)
    except OSError as error:
        print("# ERROR: cannot open/read file:", log_filename, error)
        return rpt_files

    job_done = False
    repeats_ok = True
    for line in logfile:
        repeats = re.search(r'locations to: (\S+)', line)
        if repeats:
            rpt_filename = repeats.group(1)
            if not(os.path.isfile(rpt_filename)):
                repeats_ok = False
                break
            else:
                rpt_files.append(rpt_filename)
        else: # last line in log
            summary = re.search(r'Genome length: \d+', line)
            if summary:
                job_done = True

    logfile.close

    if repeats_ok and job_done:
        return rpt_files
    else:
        return []


def run_red( red_exe, cores, gnmdirname, rptdirname, log_filepath):
    '''Calls Red, waits for job completion and logs stdout.
       Returns list of TSV repeat filenames.
       If repeat outfiles are in place then Red is skipped.
       Note: repeats are requested in format 3, 1-based inclusive (Red2)'''

    rpt_files = []

    # check whether previous results exist
    if(os.path.isfile(log_filepath)):
        rpt_files = _parse_rptfiles_from_log(log_filepath)
        if rpt_files:
            print("# re-using previous Red results")				
            return rpt_files

    # open new log file 
    try:
        logfile = open(log_filepath,"w")
    except OSError as error:
        print("# ERROR: cannot create file ", log_filepath, error)

    cmd = red_exe + \
            ' -cor '+ cores + \
            ' -frm 3'+ \
            ' -gnm ' + gnmdirname + \
            ' -rpt ' + rptdirname 

    # check Red binary
    if not(os.path.isfile(red_exe)):
        raise FileNotFoundError(errno.ENOENT,os.strerror(errno.ENOENT),red_exe)

    # run Red and capture stdout
    try:
        print("# Red command: ", cmd)
        osresponse = subprocess.check_call(cmd.split(),stdout=logfile)
    except subprocess.CalledProcessError as err:
        print("# ERROR: cannot run Red ", err.returncode)
    finally:
        logfile.close()

    # parse log and capture repeat filenames
    rpt_files = _parse_rptfiles_from_log(log_filepath)
    return rpt_files


def store_repeats_database( rptdir, seq_name_list, rpt_file_list,\
    red_path,red_version, red_params, logic_name, db_url):
    '''Store parsed Red repeats in Ensembl core database
       accessible from passed URL. Note that the analysis logic name
       and software details are also passed in order to
       fill the analysis table.
       Returns number of inserted repeats.'''

    name_to_seqregion = {}

    # core database handles
    engine = db.create_engine(db_url)
    connection = engine.connect()
    metadata = db.MetaData()
    
    # relevant db tables 
    analysis_table = db.Table('analysis',metadata,autoload=True,autoload_with=engine)
    meta_table = db.Table('meta',metadata,autoload=True,autoload_with=engine)
    repeat_consensus_table = \
        db.Table('repeat_consensus',metadata,autoload=True,autoload_with=engine)
    repeat_feature_table = \
        db.Table('repeat_feature',metadata,autoload=True,autoload_with=engine)
    seq_region_table = db.Table('seq_region',metadata,autoload=True,autoload_with=engine)
    seq_syn_table = \
        db.Table('seq_region_synonym',metadata,autoload=True,autoload_with=engine)

    # fetch seq_region_ids of sequences
    for seq_name in seq_name_list:
        seq_query = db.select([seq_region_table.columns.seq_region_id])
        seq_query = seq_query.where(seq_region_table.columns.name == seq_name)
        seq_results = connection.execute(seq_query).fetchall()
        if seq_results:
            seq_region_id = seq_results[0][0]
        else:
            # try synonyms if that failed
            syn_query = db.select([seq_syn_table.columns.seq_region_id])
            syn_query = syn_query.where(seq_syn_table.columns.synonym == seq_name)
            syn_results = connection.execute(syn_query).fetchall()
            if syn_results:
                seq_region_id = syn_results[0][0]
            else:
                print("# ERROR: cannot find seq_region_id for sequence %s\n" % seq_name)
                return 0              

        print("# sequence %s corresponds to seq_region_id %d" % (seq_name, seq_region_id))	
        name_to_seqregion[seq_name] = seq_region_id

    # insert Red analysis, fails if logic_name exists
    analysis_insert = analysis_table.insert().values({ \
        'created':db.sql.func.now(), \
        'logic_name':logic_name, \
        'program':'Red', \
        'program_version':red_version, \
        'program_file':red_path,
        'parameters': red_params,
        'gff_source':logic_name,
        'gff_feature':'repeat' })
    connection.execute(analysis_insert)

    # fetch the assigned analysis_id for the new Red analysis
    analysis_query = db.select([analysis_table.columns.analysis_id])
    analysis_query = \
        analysis_query.where(analysis_table.columns.logic_name == logic_name)
    analysis_results = connection.execute(analysis_query).fetchall()
    analysis_id = analysis_results[0][0]

    # insert repeat analysis meta keys, will fails if exists
    meta_insert = meta_table.insert().values({ \
        'species_id':1, \
        'meta_key':'repeat.analysis', \
        'meta_value':logic_name })
    connection.execute(meta_insert)

    # insert dummy repeat consensus, will fail if it exists
    # Note: Red repeats are not annotated by default, 
    # thus they are linked to a dummy repeat consensus
    repeat_consensus_insert = repeat_consensus_table.insert().values({ \
        'repeat_name':logic_name, \
        'repeat_class':logic_name, \
        'repeat_type':logic_name, \
        'repeat_consensus':'N' })
    connection.execute(repeat_consensus_insert)

    # fetch the repeat_consensus_id of the new dummy consensus
    repeat_consensus_query = \
        db.select([repeat_consensus_table.columns.repeat_consensus_id])
    repeat_consensus_query = \
        repeat_consensus_query.where( \
            repeat_consensus_table.columns.repeat_name == logic_name)
    repeat_consensus_results = connection.execute(repeat_consensus_query).fetchall()
    dummy_consensus_id = repeat_consensus_results[0][0]

    # parse repeats and produce a TSV file to be loaded in repeat table
    TSVfilename =_parse_repeats(rptdir, rpt_file_list, name_to_seqregion,\
        analysis_id, dummy_consensus_id)

    # actually insert repeat features
    repeat_query = "LOAD DATA LOCAL INFILE '" + TSVfilename +\
        "' INTO TABLE repeat_feature FIELDS TERMINATED BY '\\t' " +\
        "LINES TERMINATED BY '\\n' (seq_region_id,seq_region_start," + \
        "seq_region_end,repeat_start,repeat_end,repeat_consensus_id,analysis_id)"
    repeat_result = connection.execute(repeat_query).rowcount

    return repeat_result, name_to_seqregion


def _parse_repeats(rptdir, rpt_file_list, name2region, analysis_id, repeat_consensus_id):
    '''Parses 1-based inclusive coords produced by Red in rpt dir and 
       creates TSV file ready to be loaded in Ensembl core database.
       Returns TSV filename.''' 
	   
    if not rpt_file_list:
        print("# ERROR: got no repeat files")

    # open new TSV file
    outfilename = os.path.join(rptdir, 'ensembl.tsv')
    try:
        tsvfile = open(outfilename,"w")
    except OSError as error:
        print("# ERROR: cannot create file ", outfilename, error)

    # parse repeat coord files, one per sequence
    for filename in rpt_file_list:
        try:
            rptfile = open(filename)
        except OSError as error:
            print("# ERROR: cannot open/read file:", filename, error)

        for line in rptfile:
            column = line.split()
             
            if column[0] in name2region: 
                seq_region_id = name2region[column[0]]
            else:
                seq_region_id = column[0]
                print("# ERROR: cannot fetch seq_region_id for sequence", column[0])

            seq_region_start = column[1]
            seq_region_end = column[2]
            repeat_start = 1
            repeat_end = int(seq_region_end) - int(seq_region_start) + 1

            print("%s\t%d\t%d\t%d\t%d\t%d\t%d" % (\
                seq_region_id,\
                int(seq_region_start),\
                int(seq_region_end),\
                int(repeat_start),\
                repeat_end,\
                int(repeat_consensus_id),\
                int(analysis_id)), \
				file=tsvfile)

        rptfile.close()

    tsvfile.close()
    
    return outfilename


def main():

    parser=argparse.ArgumentParser()
    
    parser.add_argument("fasta_file",
        help="path to FASTA file with top-level genomic sequences")
    parser.add_argument("outdir",
        help="path to directory to store Red results")
    parser.add_argument("--exe", default="Red",
        help="path to Red executable, default: Red")
    parser.add_argument("--cor", default=1,
        help="number of cores for Red, default: 1")
    parser.add_argument("--host",
        help="name of the database host, required to store repeats in Ensembl core")
    parser.add_argument("--user",
        help="host user, required to store repeats in Ensembl core")
    parser.add_argument("--pw",
        help="host password, required to store repeats in Ensembl core")
    parser.add_argument("--port", type=int,
        help="host port, required to store repeats in Ensembl core")
    parser.add_argument("--db",
        help="name of the core database, required to store repeats in Ensembl core")
    parser.add_argument("--logic_name", default="repeatdetector",
        help="logic name of Ensembl analysis, default: repeatdetector")

    args = parser.parse_args()

    # create output directory and subdirs if required,
    # these follow Red nomenclature
    gnmdir = args.outdir
    rptdir = os.path.join(args.outdir, 'rpt')
    outdirs = [ gnmdir, rptdir ]
    for dir in outdirs:
        try:
            os.makedirs( dir, mode=0o777, exist_ok=True)
        except OSError as error:
            print("# ERROR: cannot create ", dir)
            print(error) 

    log_filepath = os.path.join(gnmdir, 'log.txt')

    # save individual sequences to output directory, 
    # this allows for multi-threaded Red jobs
    print("# parsing FASTA file")
    sequence_names = parse_FASTA_sequences( args.fasta_file, gnmdir)
    if len(sequence_names) == 0:
        print("# ERROR: cannot parse ", args.fasta_file)
    else:
        print("# number of input sequences = %d\n\n" % len(sequence_names))

    # run Red, or else re-use previous results
    print("# running Red")
    repeat_filenames = run_red( args.exe, args.cor, \
	    gnmdir, rptdir, log_filepath) 
    print("# TSV files with repeat coords: %s\n\n" % rptdir)

    # (optionally) store repeat features in core database
    if args.user and args.pw and args.host and args.port and args.db:
        db_url = 'mysql://' + \
                args.user + ':' + \
                args.pw + '@' + \
                args.host + ':' + \
                str(args.port) + '/' + \
                args.db + '?' + \
               'local_infile=1'

        red_version, red_params = parse_params_from_log(log_filepath)

        num_repeats, name2seqregion  = \
            store_repeats_database( rptdir, sequence_names, repeat_filenames, \
                args.exe, red_version, red_params, args.logic_name,\
                db_url)
        print("\n# stored %d repeats\n" % num_repeats);

        # print sequence name synonyms to file
        syn_filepath = os.path.join(gnmdir, 'synonyms.tsv')
        try:
            synfile = open(syn_filepath, "w")
            for name, seqregion in name2seqregion.items():
                synfile.write("%s\t%s\n" % (name, seqregion))
            synfile.close()
        except OSError as error:
            print("# ERROR: cannot create file:", syn_filepath, error)

    elif args.user or args.pw or args.host or args.port or args.db:
        print("# ERROR: make sure you set all Ensembl core params")


if __name__ == "__main__":
    main()


