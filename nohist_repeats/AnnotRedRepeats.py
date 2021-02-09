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

def _parse_minimap_log( log_filename ):
    '''Parses minimap log file and returns
       i string) version, float ii) peak RAM'''

    version = ''
    RAM = 0.0

    try:
        logfile = open(log_filename)
    except OSError as error:
        print("# ERROR: cannot open/read file:", log_filename, error)
        return version

    for line in logfile:
        vre = re.search(r'Version: (\S+)', line)
        if vre: 
            version = vre.group(1)
        
        ramre = re.search(r'Peak RSS+: (\S+)', line)
        if ramre: 
            RAM = float(ramre.group(1))
            break
						       
    logfile.close
    return version, RAM

def run_minimap( miniexe, cores, lib_filename, fasta_filename, outdir):
    '''Calls minimap, waits for job completion and logs stderr.
       If previous output exist minimap & sort are skipped.
       Returns: i string) name of file with sorted mappings (system sort).
                ii string) minimap version'''

    output_filename = os.path.join(outdir, 'repeat_mappings.paf')
    sorted_filename = os.path.join(outdir, 'repeat_mappings.sort.paf')
    log_filepath = os.path.join(outdir, 'repeat_mappings.log')

    # check whether previous results exist
    if(os.path.isfile(sorted_filename) and
        os.path.getsize(sorted_filename) > 0 and
        os.path.isfile(log_filepath)):

        # jobs might die due to insufficient RAM; make sure
        # they completed by checking the final memory report
        (minimap_version, RAM) = _parse_minimap_log(log_filepath)
        if RAM > 0:
            print("# re-using previous mappings file ",sorted_filename)
																
        return sorted_filename, minimap_version

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
    
    (minimap_version, RAM) = _parse_minimap_log(log_filepath)

    # PAF format equivalent in BLAST+:
    # blastn -task blastn -query test.fna -db nrTEplants.fna -outfmt 
    # '6 qseqid qlen qstart qend strand sseqid slen sstart send positive length bitscore'

    # perl -lane 'print join("\t",@F[0,3,6,7,2,1,3,8,9,11,3,10])' megablast>paf

    # sort results on repeat, start coord, align length, end coord
    # See: https://github.com/lh3/miniasm/blob/master/PAF.md
    # Note: this allows updating the overlap coords while parsing
    # in function make_annotation_report 
    cmd = 'sort -k1,1 -k3,3n -k10,10nr -k4,4n ' + output_filename

    try:
        osresponse = subprocess.check_call(cmd.split(),stdout=sortfile)
    except subprocess.CalledProcessError as err:
        print("# ERROR: cannot run sort ", cmd, err.returncode)
    finally:
        sortfile.close()

    os.remove(output_filename)

    return sorted_filename, minimap_version


def make_annotation_report( map_filename, log_filename, 
    minlen, verbose=False):
    '''Parses file with sorted mappings and print repeat annotation stats.
       Only alignments > minlen are considered.
       Returns dictionary with actual mapped repeat segment coords (tuple)
       as keys and a tuples of matched repeats as values'''

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
        #  7 : 0-based start coord of matched repeat
        #  8 : 0-based, exclusive end coord of matched repeat

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

            # record actual annotated repeated segment
            # Note: using Ensembl 1-based exclusive format
            repeat_coords = line.split(":")
            repeat_coords[1] = int(repeat_coords[1]) + start + 1
            repeat_coords[2] = repeat_coords[1] + (end-start) + 1

            matched_repeats[ (repeat_coords[0], repeat_coords[1], repeat_coords[2]) ] = \
                (paf[5], int(paf[7])+1, int(paf[8])+1)
				
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

    return matched_repeats 


def _annotated_repeats_to_files( workdir, matched_repeats, sequences,
    analysis_id, last_consensus_id ):
    '''Creates two TSV files in workdir to be loaded in Ensembl core db:
       i) 1-based matched repeats for repeat_feature table
       ii) annotation for repeat_consensus table
       Note1: keys in matched_repeats are tuples (seq_region, start, end)
       Note2: values in matched_repeats are tuples (repeat name, start, end)
       Returns i and ii filenames.'''

    repname_to_id = {}

    # open TSV file I (repeats)
    repfilename = os.path.join(workdir, 'repeat_feature.tsv')
    try:
        tsvfileI = open(repfilename,"w")
    except OSError as error:
        print("# ERROR: cannot create file ", repfilename, error)

    # open TSV file II (annotations/consensi)
    consfilename = os.path.join(workdir, 'repeat_consensus.tsv')
    try:
        tsvfileII = open(consfilename,"w")
    except OSError as error:
        print("# ERROR: cannot create file ", consfilename, error)

    # loop along matched repeats
    new_consensus_id = last_consensus_id + 1

    for mrep in matched_repeats.keys():
        seq_region_id = mrep[0]
        seq_region_start = mrep[1]
        seq_region_end = mrep[2]
        repeat_start = matched_repeats[mrep][1]
        repeat_end = matched_repeats[mrep][2]

        # parse repeat classification and assign repeat_consensus_id
        # example: RLG_43695:mipsREdat_9.3p_ALL#LTR/Gypsy
        repclass = matched_repeats[mrep][0].split("#")
        repeat_name = repclass[0]

        if repeat_name in repname_to_id:
            repeat_consensus_id = repname_to_id[repeat_name]
        else:
            repeat_consensus_id = new_consensus_id
            if matched_repeats[mrep][0] in sequences:
                repeat_consensus = sequences[matched_repeats[mrep][0]]
            else:
                repeat_consensus = 'N' # default
            classtype = repclass[1].split("/")
            if len(classtype) == 2:
                repeat_class = repclass[1]
                repeat_type = classtype[1]
            else:
                repeat_type = repclass[1]
                repeat_class = repclass[1]

            # only new consensi are printed to file
            print("%d\t%s\t%s\t%s\t%s" % (\
                int(repeat_consensus_id),\
                repeat_name,\
                repeat_class,\
                repeat_type,\
                repeat_consensus),\
                file=tsvfileII)

            # and stored in dictionary
            repname_to_id[repeat_name] = repeat_consensus_id		
            new_consensus_id = new_consensus_id + 1

        print("%s\t%d\t%d\t%d\t%d\t%d\t%d" % (\
            seq_region_id,\
            seq_region_start,\
            seq_region_end,\
            repeat_start,\
            repeat_end,\
            repeat_consensus_id,\
            analysis_id), \
            file=tsvfileI)

    tsvfileI.close()
    tsvfileII.close()

    return repfilename, consfilename

def _get_FASTA_sequences( filename ):
    '''Takes a FASTA filename and reads sequences.
       Return a dictionary with first non-blank word as keys and 
       sequence as values'''

    sequences = {}

    try:
        file = open(filename)
    except OSError as error:
        print("# ERROR: cannot open/read file:", filename, error)
        return sequences

    for line in file:
        header = re.search(r'^>(\S+)', line) 
        if header:
            name = header.group(1)
        else:
            if name in sequences:
                sequences[name] = sequences[name] + line.rstrip()
            else:
                sequences[name] = line.rstrip()
 
    file.close()
    return sequences


def store_annotated_repeat_database( workdir, matched_repeats, 
    exe, minimap_version, repeat_fasta_file, logic_name, db_url):
    '''Stores parsed repeat annotations in Ensembl core database
       accessible from passed URL. Second param is dictionary
       with tuples (seq_region_id,start,end) and annotations as values.
       Note that the analysis logic name and software details are 
       also passed in order to fill the analysis table.
       Returns number of inserted annotations.'''

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
    #don't read synonyms, repeat file used as input for minimap already used them
    #seq_syn_table = \
    #    db.Table('seq_region_synonym',metadata,autoload=True,autoload_with=engine)

    # insert new analysis, fails if logic_name exists
    analysis_insert = analysis_table.insert().values({ \
        'created': db.sql.func.now(), \
        'logic_name': logic_name, \
        'program':'minimap', \
        'program_version': minimap_version, \
        'program_file': exe,
        'parameters': repeat_fasta_file,
        'gff_source': logic_name,
        'gff_feature':'repeat' })
    connection.execute(analysis_insert)

    # fetch the assigned analysis_id for logic_name
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

    # insert dummy, default repeat consensus, will fail if it exists
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

    # parse repeat_fasta_file to fetch full sequences of matched repeats
    sequences = _get_FASTA_sequences(repeat_fasta_file)

    # parse repeats and produce TSV files to be loaded in db
    #analysis_id=1
    #dummy_consensus_id=1
    #repeat_result=1
    (repeat_filename, cons_filename) = _annotated_repeats_to_files(workdir,\
        matched_repeats, sequences, analysis_id, dummy_consensus_id)

    # insert repeat annotations in repeat_consensus table
    cons_query = "LOAD DATA LOCAL INFILE '" + cons_filename +\
        "' INTO TABLE repeat_consensus FIELDS TERMINATED BY '\\t' " +\
        "LINES TERMINATED BY '\\n' (repeat_consensus_id,repeat_name," + \
        "repeat_class,repeat_type,repeat_consensus)"
    consensus_result = connection.execute(cons_query).rowcount

    # insert repeat features
    repeat_query = "LOAD DATA LOCAL INFILE '" + repeat_filename +\
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
    parser.add_argument("--logic_name", default="repeatdetector_annotated",
        help="logic name of Ensembl analysis, default: repeatdetector_annotated")

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

    # run minimap, or else re-use previous results
    # Note: Red-called repeats are handled as long reads
    print("# running minimap")
    (map_filename, version) = run_minimap( args.exe, args.cor, \
        formatted_lib_filename, repeats_filename, annotdir)
    print("# mapped repeats: ", map_filename)

    matched_repeats = make_annotation_report( map_filename,\
        log_filepath,int(args.minlen))

    # make URL to connect to core database
    if args.user and args.pw and args.host and args.port and args.db:
        db_url = 'mysql://' + \
        args.user + ':' + \
        args.pw + '@' + \
        args.host + ':' + \
        str(args.port) + '/' + \
        args.db + '?' + \
        'local_infile=1'

        num_annot = store_annotated_repeat_database( annotdir, matched_repeats,\
            args.exe, version, args.repeat_fasta_file, args.logic_name, db_url)

        print("\n# stored %d repeat annotations\n" % num_annot)

    elif args.user or args.pw or args.host or args.port or args.db:
        print("# ERROR: make sure you set all Ensembl core params")


if __name__ == "__main__":
    main()


