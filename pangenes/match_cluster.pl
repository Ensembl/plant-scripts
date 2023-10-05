#!/usr/bin/env perl

# This script matches input nucleotide sequences to clusters produced by get_pangenes.pl .
# It creates a sequence index with nucleotide sequences from clusters and 
# uses GMAP to scan them. Can use cdna [default ] or CDS sequences.

# Copyright [2023]
# EMBL-European Bioinformatics Institute & Estacion Experimental Aula Dei-CSIC

$|=1;

use strict;
use warnings;
use Getopt::Std;
use File::Temp qw/ tempfile /;
use File::Basename;
use FindBin '$Bin';
use lib "$Bin/lib";
use pangeneTools qw( check_installed_features feature_is_installed 
                     parse_sequence_FASTA_file );

my @FEATURES2CHECK = (
  'EXE_GMAP'
);

my $GMAPBIN = $ENV{'EXE_GMAP'};
my $GMAPBUILDBIN = $ENV{'EXE_GMAP_BUILD'};

my ($INP_dir, $INP_seqfile, $INP_isCDS, $INP_idfile) = ('','',0,'');
my ($INP_outfile,$INP_threads, $INP_ow, $INP_verbose) = ('',4, 0, 0);
my ($cluster_regex, $isCDS, $gmapdb_path, $gmapdb, $seq) = ( '.cdna.fna', 0 );
my ($max_number_isoforms, $qlen, $tlen, $identity, $cover) = ( 0 );
my ($cluster_file, $gene_id, $isof_id, $cmd, $index_file);
my ($cluster_id, $seq_id, $cluster_seq_id, $coords, $taxon);
my (%opts, %matches, @order_geneid);

getopts('hvcCwI:d:s:t:o:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0))
{
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-c print credits and checks installation\n";
  print "-d directory with pangene clusters.pl             (required, should contain *.cdna.fna or *.cds.fna files)\n";
  print "-s nucleotide sequence file in FASTA format       (required, example: -s transcripts.fna,\n";
  print "                                                   useful to have genomic coords in header ie chr1:12-1200)\n";
  print "-o output file in TSV format                      (required)\n";
  print "-C use CDS sequences                              (optional, cDNA sequences are scanned by default)\n";
  print "-t threads                                        (optional, default: -t $INP_threads)\n";
  print "-w overwrite GMAP index                           (optional, by default index is re-used if possible)\n";
  print "-v verbose                                        (optional)\n";
  #print "-I TSV file matching cluster names to pangene ids (optional, example: -I cluster2id.tsv)\n\n";
  exit(0);
}

if(defined($opts{'c'})) {
  print "\nPrimary citation:\n https://doi.org/10.1186/s13059-023-03071-z\n";
  print "\nThis software uses external algorithms, please cite them accordingly:\n";
  print " gmap https://doi.org/10.1093/bioinformatics/bti310\n";

  # check all binaries needed by this program and print diagnostic info
  print check_installed_features(@FEATURES2CHECK);
  exit(0);
}

if(defined($opts{'d'})) { 
  $INP_dir = $opts{'d'};
  $gmapdb = basename($INP_dir) . '.gmap';
  $gmapdb_path = dirname($INP_dir);

  if(defined($opts{'C'})){
    $INP_isCDS = 1;
    $cluster_regex = '.cds.fna';
    $gmapdb = basename($INP_dir) . '.cds.gmap';
  }
}
else{ die "# EXIT : need a -d directory\n" }

if(defined($opts{'s'})){  
  $INP_seqfile = $opts{'s'};

  if(!-e $INP_seqfile) {
    die "# EXIT : need a valid -s file\n"
  }
} 
else{ die "# EXIT : need parameter -s\n" }

if(defined($opts{'o'})){
  $INP_outfile = $opts{'o'};
} 
else{ die "# EXIT : need parameter -o\n" }

if(defined($opts{'I'})){ 
  $INP_idfile = $opts{'I'};  
}

if(defined($opts{'w'})){
  $INP_ow = 1;
}

if(defined($opts{'v'})){
  $INP_verbose = 1;
}

if(defined($opts{'t'}) && $opts{'t'} >= 1){
  $INP_threads = int($opts{'t'});
}

print "# $0 -d $INP_dir -s $INP_seqfile -o $INP_outfile -C $INP_isCDS " .
  "-w $INP_ow -t $INP_threads -v $INP_verbose\n\n";

# 1) parse pre-computed cluster files, write sequences to temp file & make GMAP index
opendir(INPDIR,$INP_dir) || 
  die "# ERROR: cannot list $INP_dir , please check -d argument is a valid folder\n";
my @files = grep {/$cluster_regex$/} readdir(INPDIR);
closedir(INPDIR);

if(scalar(@files) == 0) {
  die "# ERROR: cannot find any valid clusters in $INP_dir\n";
}

# check whether previous index should be re-used
$index_file = "$gmapdb.ref153positions";

if($INP_ow ||
  !-d "$gmapdb_path/$gmapdb" ||
  !-e "$gmapdb_path/$gmapdb/$index_file") {

  my ($fh, $filename) = tempfile( '/tmp/tempfnaXXXXX', UNLINK => 1);

  foreach $cluster_file (@files) {

    my ( $ref_geneid, $ref_fasta, $ref_isof_coords, $ref_taxon ) = 
      parse_sequence_FASTA_file( "$INP_dir/$cluster_file" , 0);

    foreach $gene_id (@$ref_geneid) {  
      foreach $seq (split(/\n/,$ref_fasta->{$gene_id})) {

        if($seq =~ /^>(\S+) (\S+) (\S+) (\S+)/) {
          # FASTA headers in clusters produced by get_pangenes.pl look like this:
          #>transcript:Os01t0100100-01 gene:Os01g0100100 1:2983-10815(+) [Oryza_sativa.IRGSP-1.0.chr1]

          # prepend cluster filename to header and print to temp file
          ($isof_id, $gene_id, $coords, $taxon) = ($1, $2, $3, $4);
          print $fh ">$cluster_file|$isof_id|$gene_id|$coords|$taxon\n";

        } else {
          print $fh "$seq\n";
        }
      }
    }
  } 


  # actually create GMAP index from temp file
  $cmd = "$GMAPBUILDBIN -e 0 -d $gmapdb -D $gmapdb_path $filename 2>&1";

  print "# building GMAP index $gmapdb_path/$gmapdb\n\n";
  system($cmd);
  if ( $? != 0 ) {
    die "# ERROR: failed running $GMAPBUILDBIN ($cmd)\n";
  } else {
    print "\n# created GMAP index $gmapdb_path/$gmapdb\n\n";
  }

} else {
  print "\n# re-using GMAP index $gmapdb_path/$gmapdb\n\n";
}

## 2) map input sequence(s) to indexed cluster sequences with GMAP
$cmd = "$GMAPBIN --no-chimeras -t $INP_threads -n 100 -S -d $gmapdb -D $gmapdb_path $INP_seqfile 2>&1";
open(GMAP, "$cmd |") ||
  die "# ERROR: cannot run $cmd\n"; 
while(<GMAP>) { 

  print if($INP_verbose);
  chomp;	

  #>transcript:Os01t0147300-01 gene:Os01g0147300 1:2568051-2570029(+) [Oryza_sativa.IRGSP-1.0.chr1]
  if(/^>(.*)/){ 
    # input FASTA headers could be in any format 
    $seq_id = $1;
    push(@order_geneid, $seq_id);

  } elsif(/^  Path \d+: .*?(\d+) bp\)$/){
    #Path 1: query 1..1768 (1768 bp) => genome gene:...:1..1,820 (1820 bp)
    ($tlen) = ($1); 

  } elsif(/^    Coverage: (\S+) \(query length: (\d+) bp/){	  
    ($cover,$qlen) = ($1, $2);

  } elsif(/^    Percent identity: (\S+)/) {
    #Percent identity: 87.4 (546 matches, 62 mismatches, 17 indels, 0 unknowns)
    ($identity) = ($1);

  } elsif(/^    \+(\S+)\s+\(\d+-\d+\)/){
    # expected to be on the plus/+ strand as these are sequences cut from GFF in the right strand:
    # +cluster.cdna.fna|Os01t0829900-01|gene:Os01g0829900|1:35501996-35505052(-)|[taxon]:1-1820  (1-1768)   97%

    ($cluster_seq_id) = ($1);
    if($cluster_seq_id =~ m/^([^\|]+)\|([^\|]+)\|([^\|]+)\|([^\|]+)\|\[([^\]]+)\]:/){
      ($cluster_id, $isof_id, $gene_id, $coords, $taxon) = ($1, $2, $3, $4, $5);
    }

    # compile match stats
    $matches{$seq_id}{$cluster_id}{'total'}++;

    # take stats from hit with best cover
    if(!defined($matches{$seq_id}{$cluster_id}{'cover'}) || 
      $cover > $matches{$seq_id}{$cluster_id}{'cover'}) {
      $matches{$seq_id}{$cluster_id}{'cover'} = $cover;
      $matches{$seq_id}{$cluster_id}{'identity'} = $identity;
      $matches{$seq_id}{$cluster_id}{'coords'} = $coords;
      $matches{$seq_id}{$cluster_id}{'qlength'} = $qlen;
      $matches{$seq_id}{$cluster_id}{'tlength'} = $tlen;
    }
  }
}
close(GMAP);


## 3) print results table

open(TSV,">",$INP_outfile) ||
  die "# ERROR: cannot create $INP_outfile\n";

print TSV "#query\tqlength\tpangene\tlength\t".
  "matches\tperc_qcover\tperc_identity\tcoords\n";

foreach $seq_id (@order_geneid) {

  # no matches
  if(!defined($matches{$seq_id})) {
      print TSV "$seq_id\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
      next;
  }
	
  # note there might be 1+ rows per input sequence, 
  # a sequence can potentially match several clusters
  foreach $cluster_id (keys(%{ $matches{$seq_id} })) {

    printf(TSV "%s\t%d\t%s\t%d\t%d\t%d\t%1.1f\t%s\n",
      $seq_id,
      $matches{$seq_id}{$cluster_id}{'qlength'},
      $cluster_id,
      $matches{$seq_id}{$cluster_id}{'tlength'},
      $matches{$seq_id}{$cluster_id}{'total'},      
      $matches{$seq_id}{$cluster_id}{'cover'},
      $matches{$seq_id}{$cluster_id}{'identity'},
      $matches{$seq_id}{$cluster_id}{'coords'}
    );
  }
}

close(TSV);

print "# results in TSV format: $INP_outfile\n\n";
