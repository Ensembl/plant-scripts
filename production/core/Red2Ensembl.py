#!/usr/bin/env python3 

# Python3 script to run RepeatDetector (Red v2) and optionally 
# feed results into an Ensembl core database
# Bruno Contreras Moreira, Carlos García Girón EBI-EMBL 2020
#
# See https://www.ebi.ac.uk/seqdb/confluence/display/EnsGen/Repeat+masking+with+Red

import argparse
import os
import re
import errno
import subprocess

#import sqlalchemy as db
#import sqlalchemy_utils as db_utils

# sqlalchemy requires MySQLdb but MySQLdb doesn't support Python 3.x
# pymysql can be imported and used instead
#import pymysql
#pymysql.install_as_MySQLdb()

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


def run_red( red_exe, cores, gnmdirname, rptdirname ):
    '''Calls Red, waits for jobs and returns log file name.
       If repeat outfiles are in place then Red is skipped.
       Note: repeats are requested in format 3, 1-based inclusive (Red2)'''

    # set log file name
    log_filepath = os.path.join(gnmdirname, 'log.txt')
    logfile = open(log_filepath)

    # check whether previous results exist
    if logfile:
        job_done = False
        repeats_ok = True
        for line in logfile:
            #Printing locations to: outdir/rpt/1.tsv
            repeats = re.search(r'locations to (\S+)', line)
            summary = re.search(r'Genome length: \d+', line)
            if repeats:
                rptfile = repeats.group(1)
                if not(os.path.isfile(rptfile)):
                    repeats_ok = False
            elif summary:
                job_done = True

            if repeats_ok == False:
                break
        
        logfile.close()
        
        if job_done:
            print("# re-using previous Red results")				
            return log_filepath

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

    return log_filepath

def save_repeats_core( target_db ):
    '''Store parsed Red repeats in Ensembl core database'''    

def parse_repeats(repeat_consensus_id,analysis_id):
    '''Parse the BED-like format produced by Red in rpt dir and 
	create new file with 1-based inclusive coords for Ensembl''' 


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
    parser.add_argument("--pass",
        help="host password, required to store repeats in Ensembl core")
    parser.add_argument("--port", type=int,
        help="host port, required to store repeats in Ensembl core")
    parser.add_argument("--db",
        help="name of the core database, required to store repeats in Ensembl core")

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

    # save individual sequences to output directory, 
    # this allows for multi-threaded Red jobs
    print("# parsing FASTA file")
    n_of_sequences = parse_FASTA_sequences( args.fasta_file, gnmdir)
    if n_of_sequences == 0:
        print("# ERROR: cannot parse ", args.fasta_file)
        exit(-1)
    else:
        print("# number of input sequences = %d\n" % n_of_sequences)

    # run Red
    print("# running Red")
    log_filename = run_red( args.exe, args.cor, gnmdir, rptdir) 

    # optionally parse output and feed into core   



if __name__ == "__main__":
    main()


