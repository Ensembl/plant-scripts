#!/usr/bin/env python3

# Python3 script to run RepeatDetector (a fork of Red v2) to mask repeats,
# and optionally feed results into an Ensembl core database
#
# Tested with:
# pyenv local 3.7.6
# pip install --user sqlalchemy==1.3.23 sqlalchemy_utils pymysql
#
# Copyright [2020-24] EMBL-European Bioinformatics Institute

import argparse
import sys
import os
import re
import gzip
import bz2
import math
import errno
import subprocess

from typing import TYPE_CHECKING

if TYPE_CHECKING: # False at runtime
    import sqlalchemy as db
    import sqlalchemy_utils as db_utils
    # sqlalchemy requires MySQLdb but MySQLdb doesn't support Python 3.x
    # pymysql can be imported and used instead
    import pymysql
    pymysql.install_as_MySQLdb()


def _is_gz_file(filepath):
    """Checks a file is GZIP compressed by checking magic number"""
    with open(filepath, "rb") as test_f:
        return test_f.read(2) == b"\x1f\x8b"


def _is_bz2_file(filepath):
    """Checks a file is BZIP2 compressed by checking magic number"""
    with open(filepath, "rb") as test_f:
        return test_f.read(2) == b"BZ"


def parse_FASTA_sequences(genome_file, dirname):
    """Takes a FASTA genome file name, which might be GZIP/BZIP2-compressed,
    parses individual sequences and saves them in multiple files in named directory.
    Also outputs the estimated RAM needs for a genome.
    Note: forbidden filesystem chars are removed from sequence names.
    Returns: list of successfully parsed sequence names"""

    seq_names = []

    if _is_gz_file(genome_file) == True:
        try:
            file = gzip.open(genome_file, "rt")
        except OSError as error:
            print("# ERROR: cannot open/read file:", genome_file, error)
            return num_seqs
    elif _is_bz2_file(genome_file) == True:
        try:
            file = bz2.open(genome_file, "rt")
        except OSError as error:
            print("# ERROR: cannot open/read file:", genome_file, error)
            return num_seqs
    else:
        try:
            file = open(genome_file)
        except OSError as error:
            print("# ERROR: cannot open/read file:", genome_file, error)
            return num_seqs

    seq_filepath = ""
    prev_filepath = ""
    genome_length = 0
    for line in file:
        header = re.search(r"^>", line)
        if header:
            # check previous file was open
            prev_filepath = seq_filepath
            if prev_filepath:
                seqfile.close()

            # open temp FASTA file for this sequence only
            seq_name_match = re.search(r"^>(\S+)", line)
            if seq_name_match:
                raw_seq_name = seq_name_match.group(1)
                seq_name = re.sub(r"[\/:*<>|]", "_", raw_seq_name)

                seq_filename = seq_name + ".fa"
                seq_filepath = os.path.join(dirname, seq_filename)

                try:
                    seqfile = open(seq_filepath, "w")
                except OSError as error:
                    print("# ERROR: cannot create file ", seq_filepath, error)

                seqfile.write(">%s\n" % seq_name)
                seq_names.append(seq_name)
            else:
                print("# ERROR: cannot parse FASTA header:", header)
        else:
            genome_length = genome_length + len(line) - 1
            if seqfile:
                seqfile.write(line)

    if seqfile:
        seqfile.close()
    file.close()

    # estimate RAM needed for this genome using fitted linear function
    RAM = (13.9 * math.log10(genome_length)) - 115
    print("# genome length = %d bp" % genome_length)
    if RAM > 0:
        print("# estimated RAM needed to process this genome = %1.2f GB" % RAM)

    return seq_names


def parse_params_from_log(log_filename):
    """Parses Red stdout log and returns 3 strings:
    i) Red version ii) Parameters of this job
    iii) summary of masked repeats"""

    version = "NA"
    params = ""
    summary = ""
    try:
        logfile = open(log_filename)
    except OSError as error:
        print("# ERROR: cannot open/read file:", log_filename, error)
        return version

    for line in logfile:
        versionre = re.search(r"^Version: (\S+)", line)
        if versionre:
            version = versionre.group(1)
        else:
            paramre = re.search(r"^(-\w+: \S+)", line)
            if paramre:
                params = params + " " + paramre.group(1)
            else:
                sumre = re.search(r"^Genome length: \d+ - (Repeats length: .*)", line)
                if sumre:
                    summary = sumre.group(1)

    logfile.close

    return (version, params, summary)


def _parse_files_from_log(regex, log_filename, order=[]):
    """Parses Red stdout log and returns an (optionally orderred) list with
    the names of output files matched with regex"""

    parsed_files = []
    order_files = []

    try:
        logfile = open(log_filename)
    except OSError as error:
        print("# ERROR: cannot open/read file:", log_filename, error)
        return parsed_files

    job_done = False
    parsing_ok = True
    for line in logfile:
        match = re.search(regex, line)
        if match:
            filename = match.group(1)
            if not (os.path.isfile(filename)):
                parsing_ok = False
                break
            else:
                parsed_files.append(filename)
        else:  # last line in log
            summary = re.search(r"Genome length: \d+", line)
            if summary:
                job_done = True

    logfile.close

    if parsing_ok and job_done:
    
        # sort parsed files if required
        if order:
            for elem in order:
                elemOK = False
                for file in parsed_files:
                    thisregex = "rpt/" + elem + "\." 
                    match = re.search(thisregex, file)
                    if match:
                        order_files.append(file)
                        elemOK = True
                        break

                if elemOK == False:
                    print("# ERROR: cannot match ordered contig: ", elem) 
                                
            return order_files 

        else:
            return parsed_files
    else:
        return []


def run_red(red_exe, cores, outmskfilename, gnmdirname, rptdirname, 
    log_filepath, seqnames):
    """Calls Red, waits for job completion and logs stdout.
    Returns list of TSV repeat filenames in same order as input sequences.
    If repeat outfiles are in place then Red is skipped.
    Note: repeats are requested in format 3, 1-based inclusive (Red2)"""

    rpt_files = []

    # check whether previous results exist
    if os.path.isfile(log_filepath):
        rpt_files = _parse_files_from_log(
          "locations to: (\S+)", log_filepath, seqnames)
        if rpt_files:
            print("# re-using previous Red results")
            return rpt_files

    # open new log file
    try:
        logfile = open(log_filepath, "w")
    except OSError as error:
        print("# ERROR: cannot create file ", log_filepath, error)

    cmd = (
        red_exe
        + " -cor "
        + str(cores)
        + " -msk "
        + rptdirname
        + " -frm 3"
        + " -gnm "
        + gnmdirname
        + " -rpt "
        + rptdirname
    )

    # check Red binary
    if not (os.path.isfile(red_exe)):
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), red_exe)

    # run Red and capture stdout
    try:
        print("# Red command: ", cmd)
        osresponse = subprocess.check_call(cmd.split(), stdout=logfile)
    except subprocess.CalledProcessError as err:
        print("# ERROR: cannot run Red ", err.returncode)
        logfile.close()
        return rpt_files
    finally:
        logfile.close()

    # merge masked chromosomes if requested
    if outmskfilename:
        msk_files = _parse_files_from_log(
          "masked sequence to: (\S+)", log_filepath, seqnames)
        with open(outmskfilename, "w") as outmskfile:
            for fname in msk_files:
                with open(fname) as infile:
                    for line in infile:
                        outmskfile.write(line)
        
        print("# FASTA file with soft-masked sequences: %s" % (outmskfilename))

    # parse log and capture repeat filenames
    rpt_files = _parse_files_from_log(
      "locations to: (\S+)", log_filepath, seqnames)
    return rpt_files


def produce_BED(rpt_file_list, bed_filename):
    """Parses rpt files and produces BED file with sorted repeated ranges.
    Note: requires system sort.
    Returns number of lines in BED file."""

    if not rpt_file_list:
        print("# ERROR: got no repeat files")

    # open new BED file
    try:
        bed_filename_raw = bed_filename + ".raw"
        bedfile = open(bed_filename_raw, "w")
        num_lines = 0
    except OSError as error:
        print("# ERROR: cannot create file ", bed_filename_raw, error)
        return 0

    # parse repeat coord files, one per sequence
    for filename in rpt_file_list:
        try:
            rptfile = open(filename)
        except OSError as error:
            print("# ERROR: cannot open/read file:", filename, error)
            return 0

        for line in rptfile:
            column = line.split()
            seq_region_id = column[0]
            bed_start = int(column[1]) - 1
            bed_end = column[2]
            num_lines = num_lines + 1

            print("%s\t%s\t%s" % (seq_region_id, bed_start, bed_end), file=bedfile)

        rptfile.close()

    bedfile.close()

    # sort BED file
    cmd = "sort -k1,1 -k2,2g " + bed_filename_raw

    try:
        sortfile = open(bed_filename, "w")
    except OSError as error:
        print("# ERROR: cannot create file ", bed_filename, error)
        return 0

    try:
        osresponse = subprocess.check_call(cmd.split(), stdout=sortfile)
    except subprocess.CalledProcessError as err:
        print("# ERROR: cannot run sort ", cmd, err.returncode)
    finally:
        sortfile.close()

    os.remove(bed_filename_raw)

    return num_lines


def store_repeats_database(
    rptdir,
    seq_name_list,
    rpt_file_list,
    red_path,
    red_version,
    red_params,
    logic_name,
    description,
    label,
    db_url,
):
    """Store parsed Red repeats in Ensembl core database
    accessible from passed URL. Note that the analysis logic name
    and software details are also passed in order to
    feed the analysis & analysis_description tables.
    Returns number of inserted repeats."""

    name_to_seqregion = {}

    # core database handles
    engine = db.create_engine(db_url)
    connection = engine.connect()
    metadata = db.MetaData()

    # relevant db tables
    analysis_table = db.Table("analysis", metadata, autoload=True, autoload_with=engine)
    analysis_desc_table = db.Table(
        "analysis_description", metadata, autoload=True, autoload_with=engine
    )
    meta_table = db.Table("meta", metadata, autoload=True, autoload_with=engine)
    repeat_consensus_table = db.Table(
        "repeat_consensus", metadata, autoload=True, autoload_with=engine
    )
    repeat_feature_table = db.Table(
        "repeat_feature", metadata, autoload=True, autoload_with=engine
    )
    seq_region_table = db.Table(
        "seq_region", metadata, autoload=True, autoload_with=engine
    )
    seq_syn_table = db.Table(
        "seq_region_synonym", metadata, autoload=True, autoload_with=engine
    )

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

        print(
            "# sequence %s corresponds to seq_region_id %d" % (seq_name, seq_region_id)
        )
        name_to_seqregion[seq_name] = seq_region_id

    # insert Red analysis, fails if logic_name exists
    #   (make sure binary_path is within the limits of column 'program_file')
    binary_path = red_path
    if len(binary_path) > 80:
        binary_path = '...' + binary_path[75:]

    analysis_insert = analysis_table.insert().values(
        {
            "created": db.sql.func.now(),
            "logic_name": logic_name,
            "program": "Red",
            "program_version": red_version,
            "program_file": binary_path,
            "parameters": red_params,
            "gff_source": logic_name,
            "gff_feature": "repeat",
        }
    )
    connection.execute(analysis_insert)

    # fetch the assigned analysis_id for the new Red analysis
    analysis_query = db.select([analysis_table.columns.analysis_id])
    analysis_query = analysis_query.where(
        analysis_table.columns.logic_name == logic_name
    )
    analysis_results = connection.execute(analysis_query).fetchall()
    analysis_id = analysis_results[0][0]

    # insert repeat analysis meta keys, will fails if exists
    meta_insert = meta_table.insert().values(
        {"species_id": 1, "meta_key": "repeat.analysis", "meta_value": logic_name}
    )
    connection.execute(meta_insert)

    # insert Red analysis_description, fails if logic_name exists
    analysis_desc_insert = analysis_desc_table.insert().values(
        {"analysis_id": analysis_id, "description": description, "display_label": label}
    )
    connection.execute(analysis_desc_insert)

    # insert dummy repeat consensus, will fail if it exists
    # Note: Red repeats are not annotated by default,
    # thus they are linked to a dummy repeat consensus
    repeat_consensus_insert = repeat_consensus_table.insert().values(
        {
            "repeat_name": logic_name,
            "repeat_class": logic_name,
            "repeat_type": logic_name,
            "repeat_consensus": "N",
        }
    )
    connection.execute(repeat_consensus_insert)

    # fetch the repeat_consensus_id of the new dummy consensus
    repeat_consensus_query = db.select(
        [repeat_consensus_table.columns.repeat_consensus_id]
    )
    repeat_consensus_query = repeat_consensus_query.where(
        repeat_consensus_table.columns.repeat_name == logic_name
    )
    repeat_consensus_results = connection.execute(repeat_consensus_query).fetchall()
    dummy_consensus_id = repeat_consensus_results[0][0]

    # parse repeats and produce a TSV file to be loaded in repeat table
    TSVfilename = _parse_repeats(
        rptdir, rpt_file_list, name_to_seqregion, analysis_id, dummy_consensus_id
    )

    # actually insert repeat features
    repeat_query = (
        "LOAD DATA LOCAL INFILE '"
        + TSVfilename
        + "' INTO TABLE repeat_feature FIELDS TERMINATED BY '\\t' "
        + "LINES TERMINATED BY '\\n' (seq_region_id,seq_region_start,"
        + "seq_region_end,repeat_start,repeat_end,repeat_consensus_id,analysis_id)"
    )
    repeat_result = connection.execute(repeat_query).rowcount

    return repeat_result, name_to_seqregion


def _parse_repeats(
    rptdir, rpt_file_list, name2region, analysis_id, repeat_consensus_id
):
    """Parses 1-based inclusive coords produced by Red in rpt dir and
    creates TSV file ready to be loaded in Ensembl core database.
    Returns TSV filename."""

    if not rpt_file_list:
        print("# ERROR: got no repeat files")

    # open new TSV file
    outfilename = os.path.join(rptdir, "ensembl.tsv")
    try:
        tsvfile = open(outfilename, "w")
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

            print(
                "%s\t%d\t%d\t%d\t%d\t%d\t%d"
                % (
                    seq_region_id,
                    int(seq_region_start),
                    int(seq_region_end),
                    int(repeat_start),
                    repeat_end,
                    int(repeat_consensus_id),
                    int(analysis_id),
                ),
                file=tsvfile,
            )

        rptfile.close()

    tsvfile.close()

    return outfilename


def citation_string():
    """Return string with citation information."""

    citation = (
        "Contreras-Moreira et al (2021) https://doi.org/10.1002/tpg2.20143\n"
        + "Girgis HZ (2015) BMC Bioinformatics 16:227. doi: 10.1186/s12859-015-0654-5\n"
    )
    return citation


def main():

    default_exe = os.path.join(os.path.dirname(__file__), "../lib/Red/bin/Red")

    default_description = (
        'Repeats detected using <a href="https://bmcbioinformatics.biomedcentral.com'
        + '/articles/10.1186/s12859-015-0654-5">Red (REPeatDetector)</a>'
    )

    parser = argparse.ArgumentParser(
        description="Script to run RepeatDetector (a fork of Red v2) to mask repeats,\n"
        + "and optionally feed results into an Ensembl core database.",
        epilog="Citation:\n" + citation_string(),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "fasta_file", help="path to FASTA file with top-level genomic sequences"
    )
    parser.add_argument("outdir", help="path to directory to store Red temp results")
    parser.add_argument(
        "--exe",
        default=default_exe,
        help="path to Red executable, default: " + default_exe,
    )
    parser.add_argument("--cor", default=1, help="number of cores for Red, default: 1")
    parser.add_argument(
        "--msk_file",
        default="",
        help="name of output FASTA file with soft-masked sequences",
    )
    parser.add_argument(
        "--bed_file",
        default="",
        help="name of output BED file with repeated ranges, uses original sequence names",
    )
    parser.add_argument(
        "--host",
        help="name of the database host, required to store repeats in Ensembl core",
    )
    parser.add_argument(
        "--user", help="host user, required to store repeats in Ensembl core"
    )
    parser.add_argument(
        "--pw", help="host password, required to store repeats in Ensembl core"
    )
    parser.add_argument(
        "--port", type=int, help="host port, required to store repeats in Ensembl core"
    )
    parser.add_argument(
        "--db",
        help="name of the core database, required to store repeats in Ensembl core",
    )
    parser.add_argument(
        "--logic_name",
        default="repeatdetector",
        help="logic name of Ensembl analysis, default: repeatdetector",
    )
    parser.add_argument(
        "--description",
        default=default_description,
        help="quoted string with Ensembl analysis description, default: "
        + default_description,
    )
    parser.add_argument(
        "--displaylabel",
        default="Repeats:Red",
        help="string with Ensembl analysis display label, default: Repeats:Red",
    )

    args = parser.parse_args()

    # create output directory & subdirs if required,
    # these follow Red nomenclature
    gnmdir = args.outdir
    rptdir = os.path.join(args.outdir, "rpt")
    outdirs = [gnmdir, rptdir]
    for dir in outdirs:
        try:
            os.makedirs(dir, mode=0o777, exist_ok=True)
        except OSError as error:
            print("# ERROR: cannot create ", dir)
            print(error)

    log_filepath = os.path.join(gnmdir, "log.txt")

    # save individual sequences to output directory,
    # this allows for multi-threaded Red jobs
    print("# parsing FASTA file")
    sequence_names = parse_FASTA_sequences(args.fasta_file, gnmdir)
    if len(sequence_names) == 0:
        print("# ERROR: cannot parse ", args.fasta_file)
    else:
        print("# number of input sequences = %d\n\n" % len(sequence_names))

    # run Red, or else re-use previous results
    print("# running Red")
    repeat_filenames = run_red(
        args.exe, args.cor, args.msk_file, 
        gnmdir, rptdir, 
        log_filepath, sequence_names
    )
    
    if repeat_filenames:
        print("# TSV files with repeat coords: %s\n\n" % rptdir)
    else:
        print("# Red process interrupted\n")
        sys.exit(-1)

    # output BED if requested
    if args.bed_file:
        num_bed_lines = produce_BED(repeat_filenames, args.bed_file)
        print("# BED file with repeat coords: %s\n\n" % args.bed_file)

    # (optionally) store repeat features in core database
    if args.user and args.pw and args.host and args.port and args.db:
        db_url = (
            "mysql://"
            + args.user
            + ":"
            + args.pw
            + "@"
            + args.host
            + ":"
            + str(args.port)
            + "/"
            + args.db
            + "?"
            + "local_infile=1"
        )

        (red_version, red_params, red_summary) = parse_params_from_log(log_filepath)

        num_repeats, name2seqregion = store_repeats_database(
            rptdir,
            sequence_names,
            repeat_filenames,
            args.exe,
            red_version,
            red_params,
            args.logic_name,
            args.description,
            args.displaylabel,
            db_url,
        )
        print("\n# stored %d repeats\n" % num_repeats)

        # text report
        print(
            "# summary report:\nRepeated sequences were called with the Repeat Detector,"
            + " which is part of the [Ensembl Genomes repeat feature pipelines]"
            + "(http://plants.ensembl.org/info/genome/annotation/repeat_features.html). %s\n"
            % red_summary
        )

        # print sequence name synonyms to file
        syn_filepath = os.path.join(gnmdir, "synonyms.tsv")
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
