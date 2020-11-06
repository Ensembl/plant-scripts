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
            print("# WARNING: cannot find synonyms file, will use real names")
        else:    
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
                    seq_region_id = seq_name
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
       Returns name of file with sorted mappings (system sort).
       If previous output exist minimap2 & sort are skipped.'''

    output_filename = os.path.join(outdir, 'repeat_mappings.paf')
    sorted_filename = os.path.join(outdir, 'repeat_mappings.sort.paf')
    log_filepath = os.path.join(outdir, 'repeat_mappings.log')

    # check whether previous results exist
    if(os.path.isfile(sorted_filename) and
        os.path.getsize(sorted_filename) > 0 and
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
                print("# re-using previous mappings file ",sorted_filename)
																
        return sorted_filename

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

    try:
        sortfile = open(sorted_filename,"w")
    except OSError as error:
        print("# ERROR: cannot create file ", sorted_filename, error)
        return ''

    # put together minimap2 command
    # map-ont is actually the default as of Sep2020
    # Note: also tested map-pab, asm20, little difference
    cmd = miniexe + \
            ' -K100M --score-N 0 ' +\
            ' -x map-ont ' +\
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
        print("# ERROR: cannot run minimap2 ", cmd, err.returncode)
    finally:
        logfile.close()
        outfile.close()
    print(cmd)
    
    # PAF format equivalent in BLAST+:
    # blastn -task blastn -query test.fna -db nrTEplants.fna -outfmt '6 qseqid qlen qstart qend strand sseqid slen  sstart send positive length bitscore evalue'

    # sort results on repeat, start coord, end coord
    # See: https://github.com/lh3/miniasm/blob/master/PAF.md
    # Note: this allows updating the overlap coords while parsing
    # in function make_annotation_report 
    cmd = 'sort -k1,1 -k3,3n -k4,4n ' + output_filename

    try:
        osresponse = subprocess.check_call(cmd.split(),stdout=sortfile)
    except subprocess.CalledProcessError as err:
        print("# ERROR: cannot run sort ", cmd, err.returncode)
    finally:
        sortfile.close()

    os.remove(output_filename)

    return sorted_filename


def make_annotation_report( map_filename, log_filename, 
    minlen, verbose=False):
    '''Parse file with sorted mappings and print repeat annotation stats.
       Only alignments > minlen are considered.
       Returns dictionary with mapped repeats.'''

    rep_match = {}
    matched_repeats = {}
    annotated_length = 0
    stats = {}

    try:
        mapfile = open(map_filename)
    except OSError as error:
        print("# ERROR: cannot open/read file ", map_filename, error)
        return {}

    for line in mapfile:
        paf = line.split("\t") 
        # See https://github.com/lh3/miniasm/blob/master/PAF.md
        #  0 : masked seq_region
        #  2 : 0-based start coord of seq_region
        #  3 : 0-based, exclusive end coord of seq_region
        #  5 : repeat name (ie RLG_159077:mipsREdat_9.3p_ALL#LTR/Gypsy)

        # work out length of annotated repeat
        start = int(paf[2])
        end = int(paf[3])
        annotlen = end-start

        if annotlen >= minlen:

            # update last aligned coord of this repeat,
            # compute overlap with previous repeat,
            # skip short or redundant alignments

            # -----------
            #               -------  
            if paf[0] not in rep_match or start >= rep_match[paf[0]]:
                rep_match[paf[0]] = end
            elif start < rep_match[paf[0]] and end <= rep_match[paf[0]]:
                # -------------- redundant
                #     ------
                #print("> %s %s %d %d %d" % (paf[0], paf[5], start, end, annotlen))
                continue
            else:
                # -----------
                #     -------***
                overlap = rep_match[paf[0]] - start
                annotlen = annotlen - overlap
                if annotlen >= minlen:
                    rep_match[paf[0]] = end
                    start = start + overlap
                else:
                    #print(">> %s %s %d %d %d" % (paf[0], paf[5], start, end, annotlen))
                    continue

            # record this annotated segment
            if paf[0] not in matched_repeats: 
                matched_repeats[paf[0]] = [ paf[5] ]
            else:
                matched_repeats[paf[0]].append(paf[5])
				
            # collect stats
            annotated_length = annotated_length + annotlen
            classre = re.search(r'#(\S+)', paf[5])
            if classre:
                repclass = classre.group(1)
                if repclass in stats:
                    stats[repclass] = stats[repclass] + annotlen
                else: 
                    stats[repclass] = annotlen

            if verbose:
                print("# %s %s %d %d %d" % 
                    (paf[0], paf[5], start, end, annotlen))
 
        #else: 
        #     print(">>> %s %s %d %d %d" % (paf[0], paf[5], start, end, annotlen))

    mapfile.close()

    # fetch summary from log and print it with annotation stats
    try:
        logfile = open(log_filename)
    except OSError as error:
        print("# ERROR: cannot open/read file:", log_filename, error)
        return {}

    for line in logfile:
        statsre = re.search(r'^Genome length: (\d+) - Repeats length: (\d+)', line)
        if statsre:
            print("# Genome length: %s Repeated content: %s %2.1f%% Annotated: %d %2.1f%%\n" %
                (statsre.group(1), statsre.group(2),\
                 100*int(statsre.group(2))/int(statsre.group(1)),
                 annotated_length, 100*annotated_length/int(statsre.group(1))))
   
    logfile.close()

    # print repeat class stats
    for repclass in sorted(stats.keys()): 
        print("%s\t%d" %
            (repclass, stats[repclass]))


    # parse repeat FASTA file to fetch sequences of matched repeats
    # TO BE DONE if those sequences/consensi are valuable

    return matched_repeats 


def store_repeat_annot_database( workdir, matched_repeats, 
    repeat_fasta_file, logic_name, db_url):
    '''Store parsed repeat annotations in Ensembl core database
       accessible from passed URL. First param is dictionary
       with seq_region_id:start:end keys and annotations as values.
	   Note analysis table entry with logic name is updated.
       Returns number of inserted annotations.'''

#repeat_seq_query = "SELECT r.seq_region_id, r.seq_region_start, r.seq_region_end, " +\
    #    "SUBSTR(d.sequence,r.seq_region_start,r.seq_region_end-r.seq_region_start+1) " +\
	    #    "FROM repeat_feature r INNER JOIN dna d USING (seq_region_id)"
		    #repeat_seq_result = connection.execute(repeat_seq_query).fetchall()
			    #print(repeat_seq_result[0][0])


    num_annot = 0

    # core database handles
    engine = db.create_engine(db_url)
    connection = engine.connect()
    metadata = db.MetaData()
    
    # relevant db tables 
    analysis_table = db.Table('analysis',metadata,autoload=True,autoload_with=engine)
    repeat_feature_table = \
        db.Table('repeat_feature',metadata,autoload=True,autoload_with=engine)

    # create temporary table with annotated repeat features
    #for repeat in matched_repeats:
        #(seq_region_id,seq_region_start,seq_region_end) = repeat.split(":")
        #print("%s %s %s" % (seq_region_id,seq_region_start,seq_region_end))
        
    # TO BE DONE
    return 0
    

    #update_stmt = analysis_table.update()\
    #    .where(analysis_table.logic_name=logic_name)\
    #    .values(db_file=repeat_fasta_file)
    #connection.execute(update_stmt)

    # fetch the assigned analysis_id for logic_name
    analysis_query = db.select([analysis_table.columns.analysis_id])
    analysis_query = \
        analysis_query.where(analysis_table.columns.logic_name == logic_name)
    analysis_results = connection.execute(analysis_query).fetchall()
    analysis_id = analysis_results[0][0]

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
        help="min length of repeats to be annotated, default: 90bp")
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

    # format reference sequences for minimap2
    # Note: a precomputed repeat library is used as reference
    formatted_lib_filename = format_reference_minimap( args.exe, args.cor,\
        args.repeat_fasta_file, annotdir)

    # run minimap2, or else re-use previous results
    # Note: Red-called repeats are handled as long reads
    print("# running minimap2")
    map_filename = run_minimap( args.exe, args.cor, \
        formatted_lib_filename, repeats_filename, annotdir)
    print("# mapped repeats: ", map_filename)

    matched_repeats = make_annotation_report( map_filename, log_filepath,\
        int(args.minlen))

    # make URL to connect to core database
    if args.user and args.pw and args.host and args.port and args.db:
        db_url = 'mysql://' + \
        args.user + ':' + \
        args.pw + '@' + \
        args.host + ':' + \
        str(args.port) + '/' + \
        args.db + '?' + \
        'local_infile=1'

        num_annot = store_repeat_annot_database( annotdir, matched_repeats, \
            args.repeat_fasta_file, args.logic_name, db_url)

        print("\n# stored %d repeat annotations\n" % num_annot)



if __name__ == "__main__":
    main()


