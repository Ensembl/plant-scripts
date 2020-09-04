#!/usr/bin/env python3 

# Python3 script to run RepeatDetector (Red) and feed results into an Ensembl core database
# Bruno Contreras Moreira, Carlos Garcia Giron EBI-EMBL 2020
#
# See https://www.ebi.ac.uk/seqdb/confluence/display/EnsGen/Repeat+masking+with+Red

import argparse
import os
#import filecmp
#import errno
#import subprocess

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


def run_red( red_exe, red_args ):
    ''' calls Red execute and waits for it to terminate'''
        # output format: 1 (chrName:start-end) or 2 (chrName start end)
        # Note that chrName includes the '>' character
        #cmd = self.param('red_path')+ \
        #      ' -frm 2'+ \
        #      ' -gnm '+self.param('gnm')+ \
        #      ' -rpt '+self.param('rpt')
        #try:
        #    response = subprocess.check_call(cmd.split())
        #except subprocess.CalledProcessError as err:
        #    print("Could not run Red. Return code "+err.returncode)


def save_repeats_core( target_db ):
    '''Store parsed Red repeats in Ensembl core database'''    

def parse_repeats(repeat_consensus_id,analysis_id):
    '''Parse the BED-like format produced by Red in rpt dir and 
	create new file with 1-based inclusive coords for Ensembl''' 


def main():

    parser=argparse.ArgumentParser()
    
    parser.add_argument("fasta_file",
        help="path to an input FASTA file with genomic sequences")
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

    # check script params

    # create output directory and subdirs if required,
    # these follow Red nomenclature
	gnmdir = args.outdir
	rptdir = args.outdir+'/rpt'
    outdirs = [ gnmdir, rptdir ]
    for dir in outdirs:
        try:
            os.makedirs( dir, mode=0o777, exist_ok=True)
        except OSError as error:
            print("# ERROR: cannot create " + dir)
            print(error) 

    # save individual sequences to output directory, 
    # this allows for multi-threaded Red jobs
    n_of_sequences = parse_FASTA_sequences( args.fasta_file, gnmdir )


    # run Red
    # optionally parse output and feed into core   



if __name__ == "__main__":
    main()


