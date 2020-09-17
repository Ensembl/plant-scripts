#!/usr/bin/env python3 

# Python3 script to annotate Red repeats and optionally 
# feed the new consensus_repeats into an Ensembl core database.
#
# Tested with: 
# pyenv local 3.7.6
# pip install --user sqlalchemy_utils pymysql
# 
# Bruno Contreras Moreira EBI-EMBL 2020
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

def fetch_repeats_FASTA( logpath, synpath, annotdir, minlen ):
    '''Parses previous log and synonyms to retrieve repeat coords and puts 
       the corresponding sequence segments in a new FASTA file.
       Only repeats with length >= minlen are fetched.
       Returns name of FASTA file.'''

    #repeat_seq_query = "SELECT r.seq_region_id, r.seq_region_start, r.seq_region_end, " +\
    #    "SUBSTR(d.sequence,r.seq_region_start,r.seq_region_end-r.seq_region_start+1) " +\
    #    "FROM repeat_feature r INNER JOIN dna d USING (seq_region_id)"
    #repeat_seq_result = connection.execute(repeat_seq_query).fetchall()
    #print(repeat_seq_result[0][0])

    repeat_FASTA_file = ''
    name_to_seqregion = {}

    # check whether previous files and synonyms exist
    if(os.path.isfile(logpath)):
        rpt_files = _parse_rptfiles_from_log(logpath)
        if not rpt_files:
            print("# ERROR: cannot find previous repeat files from", logpath)
            return repeat_FASTA_file
        
        fa_files = _parse_fafiles_from_log(logpath)
        if not rpt_files:
            print("# ERROR: cannot find previous sequence files from", logpath)
            return repeat_FASTA_file
	
        try:
            synfile = open(synpath)
        except OSError as error:
            print("# ERROR: cannot open/read file:", synpath, error)
            return repeat_FASTA_file
        for line in synfile:
            synre = re.search(r'^(\S+)\t(\d+)', line)
            if synre:
                name_to_seqregion[ synre.group(1) ] = synre.group(2)
 
    # create output FASTA file
    repeat_FASTA_file = os.path.join(annotdir, 'Red_repeats.fna')
    try:
        outfile = open(repeat_FASTA_file,"w")
    except OSError as error:
        print("# ERROR: cannot create file ", repeat_FASTA_file, error)

    # cut repeat sequences from coordinates, chr by chr
    for file in fa_files:

        chrsequence = ''

        # build rpt_filename
        rpt_filename = os.path.basename(file)
        rpt_filename = rpt_filename.replace(".fa",".tsv")

        # open sequence file
        try:
            fafile = open(file)
        except OSError as error:
            print("# ERROR: cannot open/read file:", file, error)
            outfile.close()
            os.remove(repeat_FASTA_file)
            return ''

        for line in fafile:
            namere = re.search(r'^>(\S+)', line)
            if namere:
                seq_name = namere.group(1)
                if seq_name in name_to_seqregion:
                    seq_region_id = name_to_seqregion[ seq_name ]
                else:
                    print("# ERROR: cannot find seg_region_id for sequence ", 
                        seq_name, error)
                    outfile.close()
                    os.remove(repeat_FASTA_file)
                    return ''
            else:
                chrsequence = chrsequence + line.rstrip()
        fafile.close() 
    
	    # open corresponding repeat file
        for rfilename in rpt_files:
            rpt_basename = os.path.basename(rfilename)
            if rpt_basename != rpt_filename: 
                continue
            try:
                rfile = open(rfilename)
            except OSError as error:
                 print("# ERROR: cannot open/read file:", rfilename, error)
                 outfile.close()
                 os.remove(repeat_FASTA_file)
                 return ''
            for line in rfile:
                ( name, start, end ) = line.split()
                # convert coord strings to substring indexes (0-based excl)
                idx_start = int(start)-1
                idx_end = int(end)
                if(idx_end-idx_start+1 < minlen): 
                    continue
                if chrsequence[idx_start:idx_end] == '':
                    print("# ERROR: cannot cut seq_region %s:%s-%s" 
                        % (seq_region_id, start, end)) 
                    outfile.close()
                    os.remove(repeat_FASTA_file)
                    return ''
                else:
                    # note that coords in FASTA header are still 1-based incl
                    outfile.write(">%s:%s:%s\n%s\n" % 
                        (seq_region_id, start, end,
                        chrsequence[idx_start:idx_end]))
    
    outfile.close()
    return repeat_FASTA_file


def format_reference_minimap( miniexe, cores, fasta_file, outdir):
    '''Takes FASTA library and formats it with minimap2.
       Returns name of formatted library file.'''

    # output for nrTEplantsJune2020.fna (708MB) -> requires 3.7GB RAM
	#[M::main::44.595*0.98] loaded/built the index for 171104 target sequence(s)
	#[M::mm_idx_stat] kmer size: 15; skip: 10; is_hpc: 0; #seq: 171104
	#[M::main] Version: 2.17-r974-dirty
	#[M::main] CMD: github/minimap2/minimap2 --idx-no-seq -d repeat_lib.mmi
	#[M::main] Real time: 47.044 sec; CPU: 44.924 sec; Peak RSS: 3.713 GB

    log_filepath = os.path.join(outdir, 'repeat_lib.mmi.log')
    formatted_filename = os.path.join(outdir, 'repeat_lib.mmi')

    # check whether previous results exist
    if(os.path.isfile(formatted_filename) and
        os.path.getsize(formatted_filename) > 0 and
        os.path.isfile(log_filepath)):
        try:
            logfile = open(log_filepath)
        except OSError as error:
            print("# ERROR: cannot open/read file ", log_filepath, error)

        for line in logfile:
            # jobs might die due insufficient RAM; make sure
            # they completed by checking the final memory report
            ramre = re.search(r'Peak RSS+', line)
            if ramre:
                print("# re-using previously formatted repeat library ", 
                    formatted_filename)
                return formatted_filename

    # open new log file
    try:
        logfile = open(log_filepath,"w")
    except OSError as error:
        print("# ERROR: cannot create file ", log_filepath, error)

    cmd = miniexe + ' -t ' + cores + ' --idx-no-seq -d '\
        + formatted_filename + ' ' + fasta_file
																	
    # check minimap2 binary
    if not(os.path.isfile(miniexe)):
        raise FileNotFoundError(errno.ENOENT,os.strerror(errno.ENOENT),miniexe)
	
    # run minimap2 and capture stderr
    try:
        osresponse = subprocess.check_call(cmd.split(),stderr=logfile)
    except subprocess.CalledProcessError as err:
        print("# ERROR: cannot run minimap2 ", cmd, err.returncode)
        # does not seem to work on bshell jobs
        #except OSError as err:
        #    print("# ERROR: need more RAM to run minimap2 ", err.returncode)
    finally:
        logfile.close()

    return formatted_filename


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

def _parse_fafiles_from_log( log_filename ):
    '''Parses names of FASTA (.fa) sequence files from log. 
       Returns a list with the filenames'''

    fa_files = []

    try:
        logfile = open(log_filename)
    except OSError as error:
        print("# ERROR: cannot open/read file:", log_filename, error)
        return fa_files

    fafiles_ok = True
    for line in logfile:
        fare = re.search(r'^Counting k-mers in (\S+) ...', line)
        if fare:
            fa_filename = fare.group(1)
            if not(os.path.isfile(fa_filename)):
                fafiles_ok = False
                break
            else:
                fa_files.append(fa_filename)

    logfile.close

    if fafiles_ok:
        return fa_files
    else:
        return []


def run_minimap( miniexe, cores, lib_filename, fasta_filename, outdir):
    '''Calls minimap2,  waits for job completion and logs stderr.
       Returns name of file with mappings.
       If previous output exist minimap2 is skipped.'''

    output_filename = os.path.join(outdir, 'repeat_mappings.paf') 
    log_filepath = os.path.join(outdir, 'repeat_mappings.log')

    # check whether previous results exist
    if(os.path.isfile(output_filename) and
        os.path.getsize(output_filename) > 0 and
        os.path.isfile(log_filepath)):
        try:
            logfile = open(log_filepath)
        except OSError as error:
            print("# ERROR: cannot open/read file ", log_filepath, error)
            return ''

        for line in logfile:
            # jobs might die due to insufficient RAM; make sure
            # they completed by checking the final memory report
            ramre = re.search(r'Peak RSS+', line)
            if ramre:
                print("# re-using previous mappings file ",output_filename)
																
        return output_filename

    # open new log & output files 
    try:
        logfile = open(log_filepath,"w")
    except OSError as error:
        print("# ERROR: cannot create file ", log_filepath, error)
        return ''
    try:
        outfile = open(output_filename,"w")
    except OSError as error:
        print("# ERROR: cannot create file ", output_filename, error)
        return ''

    cmd = miniexe + \
            ' -K100M --score-N 0 -x map-ont ' +\
            ' -t '+ cores + ' ' +\
            lib_filename + ' ' +\
            fasta_filename

    # check binary
    if not(os.path.isfile(miniexe)):
        raise FileNotFoundError(errno.ENOENT,os.strerror(errno.ENOENT),miniexe)

    # run minimap2 and capture stdout & stderr
    try:
        osresponse = subprocess.check_call(cmd.split(),stdout=outfile,stderr=logfile)
    except subprocess.CalledProcessError as err:
        print("# ERROR: cannot run minimap2 ", cmd,  err.returncode)
    finally:
        logfile.close()
        outfile.close()

    return output_filename


def store_repeats_database( rptdir, seq_name_list, rpt_file_list,\
    red_path,red_version, red_params, logic_name, db_url):
    '''Store parsed Red repeats in Ensembl core database
       accessible from passed URL. Note that the analysis logic name
       and software details are also passed in order to
       fill the analysis table.
       Returns number of inserted repeats.'''

    num_repeats = 0
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
                return num_repeats              

        print("# sequence %s corresponds to seq_region_id %d" % (seq_name, seq_region_id))	
        name_to_seqregion[seq_name] = seq_region_id

    # insert Red analysis, will fails if logic_name exists
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

    return repeat_result


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
   
    parser.add_argument("repeat_fasta_file",
        help="path to FASTA file with repeat sequences in RepBase format")
    parser.add_argument("outdir",
        help="path to directory with stored Red results")
    parser.add_argument("--exe", default="minimap2",
        help="path to minimap2 executable, default: minimap2")
    parser.add_argument("--cor", default=1,
        help="number of cores for minimap2, default: 1")
    parser.add_argument("--minlen", default=90,
        help="min length of Red repeats to be annotated, default: 90")
    parser.add_argument("--host",
        help="name of the database host")
    parser.add_argument("--user",
        help="host user")
    parser.add_argument("--pw",
        help="host password")
    parser.add_argument("--port", type=int,
        help="host port")
    parser.add_argument("--db",
        help="name of the core database")
    parser.add_argument("--logic_name", default="repeatdetector",
        help="logic name of Ensembl analysis, default: repeatdetector")

    args = parser.parse_args()

    # create output subdir if required,
    # these follow Red nomenclature
    gnmdir = args.outdir
    annotdir = os.path.join(args.outdir, 'annot') 
    try:
        os.makedirs( annotdir, mode=0o777, exist_ok=True)
    except OSError as error:
        print("# ERROR: cannot create ", annotdir)
        print(error) 

    log_filepath = os.path.join(gnmdir, 'log.txt')
    syn_filepath = os.path.join(gnmdir, 'synonyms.tsv')

    # fetch sequences of Red repeats and save in FASTA file
    repeats_filename = fetch_repeats_FASTA( log_filepath, syn_filepath,\
        annotdir, int(args.minlen) )
    print("# FASTA file with repeat sequences (length>%s): %s\n\n"\
        % (args.minlen, repeats_filename))

    # format repeat library for minimap2
    formatted_lib_filename = format_reference_minimap( args.exe, args.cor,\
        args.repeat_fasta_file, annotdir)

    # run minimap2, or else re-use previous results
    print("# running minimap2")

    map_filename = run_minimap( args.exe, args.cor, \
	    formatted_lib_filename, repeats_filename, annotdir)


    # make URL to connect to core database
    if args.user and args.pw and args.host and args.port and args.db:
        db_url = 'mysql://' + \
        args.user + ':' + \
        args.pw + '@' + \
        args.host + ':' + \
        str(args.port) + '/' + \
        args.db + '?' + \
        'local_infile=1'

#        num_annot = store_repeat_annot_database( map_filename, \
#            args.repeat_fasta_file, args.logic_name, db_url)
#        print("\n# stored %d repeat annotations\n" % num_annot);



if __name__ == "__main__":
    main()


