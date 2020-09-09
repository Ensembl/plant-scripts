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
       Returns: integer with number of parsed sequences'''
    num_seqs = 0
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
                   seqfile = open(seq_filepath,"w+")
               except OSError as error:
                   print("# ERROR: cannot create file ", seq_filepath, error)

               seqfile.write(">%s\n" % seq_name)
               num_seqs = num_seqs + 1
           else:
               print("# ERROR: cannot parse FASTA header:", header)
       else:
           if seqfile:
               seqfile.write(line)
    
    if seqfile: 
        seqfile.close()
    file.close()

    return num_seqs

def parse_version_from_log( log_filename ):
    '''Parses Red stdout log and returns a string with version'''
    version = 'NA'
    try:
        logfile = open(log_filename)
    except OSError as error:
        print("# ERROR: cannot open/read file:", log_filename, error)
        return version

    for line in logfile:
        versionre = re.search(r'^Version: (\S+)', line)
        if versionre:
            version = versionre.group(1)
            break

    logfile.close

    return version


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
        logfile = open(log_filepath,"w+")
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
        print("# ERROR: cannot run Red " + err.returncode)
    finally:
        logfile.close()

    # parse log and capture repeat filenames
    rpt_files = _parse_rptfiles_from_log(log_filepath)
    return rpt_files


def store_repeats_database( rpt_file_list, red_path, red_version, logic_name, db_url):
    '''Store parsed Red repeats in Ensembl core database
       accessible from passed URL. Note that the analysis logic name
       and the software version are also passed in order to
       fill the analysis table'''

    engine = db.create_engine(db_url)
    connection = engine.connect()
    metadata = db.MetaData()
    
    # handles for relevant db tables 
    analysis_table = db.Table('analysis',metadata,autoload=True,autoload_with=engine)
    meta_table = db.Table('meta',metadata,autoload=True,autoload_with=engine)
    repeat_consensus_table = \
        db.Table('repeat_consensus',metadata,autoload=True,autoload_with=engine)
    repeat_feature_table = \
        db.Table('repeat_feature',metadata,autoload=True,autoload_with=engine)

    # insert Red analysis
    analysis_insert = analysis_table.insert().values({ \
        'created':db.sql.func.now(), \
        'logic_name':logic_name, \
        'program':'Red', \
        'program_version':red_version, \
        'program_file':red_path })
    connection.execute(analysis_insert)

    # fetch the assigned analysis_id for the new Red analysis
    analysis_query = db.select([analysis_table.columns.analysis_id])
    analysis_query = \
        analysis_query.where(analysis_table.columns.logic_name == logic_name)
    analysis_results = connection.execute(analysis_query).fetchall()
    analysis_id = analysis_results[0][0]

    # insert repeat analysis meta keys
    meta_insert = meta_table.insert().values({ \
        'species_id':1, \
        'meta_key':'repeat.analysis', \
        'meta_value':logic_name })
    connection.execute(meta_insert)

    # insert dummy repeat consensus
    # Note: Red repeats are not annotated by default 
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
    repeat_consensus_id = repeat_consensus_results[0][0]

    # _parse_repeats( rpt_file_list, 1, 1)


def _parse_repeats(rpt_file_list, repeat_consensus_id, analysis_id):
    '''Parse the 1-based inclusive coords produced by Red in rpt dir and 
       create TSV file to be loaded in Ensembl core database''' 
    if not rpt_file_list:
        print("# ERROR: got no repeat files")

    for filename in rpt_file_list:
        try:
            rptfile = open(filename)
        except OSError as error:
            print("# ERROR: cannot open/read file:", filename, error)
        
        print(filename)


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
        help="path to Red executable, default: repeatdetector")

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
    n_of_sequences = parse_FASTA_sequences( args.fasta_file, gnmdir)
    if n_of_sequences == 0:
        print("# ERROR: cannot parse ", args.fasta_file)
    else:
        print("# number of input sequences = %d\n\n" % n_of_sequences)

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
                args.port + '/' + \
                args.db + '?' + \
               'local_infile=1'

        store_repeats_database(repeat_filenames, \
            arg.exe, parse_version_from_log(log_filepath), \
            args.logic_name, db_url)


if __name__ == "__main__":
    main()


