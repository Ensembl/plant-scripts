#!/usr/bin/env python3 

# Python3 script to run RepeatDetector (Red) and feed results into an Ensembl core database
# Bruno Contreras Moreira, Carlos Garcia Giron EBI-EMBL 2020
#
# See https://www.ebi.ac.uk/seqdb/confluence/display/EnsGen/Repeat+masking+with+Red

import argparse
import os
import re
#import filecmp
#import errno
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
               seqfile = open(seq_filepath,"w+")
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


def run_red( red_exe, gnmdirname, rptdirname ):
    '''Calls Red and waits for it to terminate.
       Note repeats are requested in format 2: chrName start end'''
    cmd = red_exe + \
            ' -frm 2'+ \
            ' -gnm ' + gnmdirname + \
            ' -rpt ' + rptdirname
    print(cmd)

    try:
        response = subprocess.check_call(cmd.split())
    except subprocess.CalledProcessError as err:
        print("# ERROR: cannot run Red " + err.returncode)


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
            print("# ERROR: cannot create " + dir)
            print(error) 

    # save individual sequences to output directory, 
    # this allows for multi-threaded Red jobs
    print("# parsing FASTA file")
    n_of_sequences = parse_FASTA_sequences( args.fasta_file, gnmdir )
    print("# number of input sequences = %d\n" % n_of_sequences)

    # run Red
    print("# running Red")
    run_red( args.exe, gnmdir, rptdir ) 

    # optionally parse output and feed into core   



if __name__ == "__main__":
    main()


