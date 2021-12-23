#!/usr/bin/env perl 

# This program computes whole genome alignments (WGA) to define clusters 
# of collinear, orthologous genes/features annotated in GFF files. Such 
# clusters might help define pan-genes across a pangenome.
# Several WGA algorithms are available and some parameters are customizable. 
# It is designed to process (in a multicore computer or HPC cluster) files 
# contained in a directory (-d), so that new .fna & .gff files can be added 
# while conserving previous results. 

# Copyright [2021-22] 
# EMBL-European Bioinformatics Institute & Estacion Experimental Aula Dei-CSIC

$|=1;

use strict;
use warnings;
use Getopt::Std;
use Fcntl qw(:flock);
use File::Basename;
use Cwd;
use FindBin '$Bin';
use lib "$Bin/lib";
use HPCluster;
use pangeneTools;

my $VERSION = '1.x';

## global variables that control some algorithmic choices
my $NOFSAMPLESREPORT = 20;        # number of samples used for the generation of pan/core genomes
my $MINQUAL    = 50;              # minimap2 mapping quality
my $MINOVERLAP = 0.50;

## list of features/binaries required by this program (do not edit)
my @FEATURES2CHECK = ('EXE_MINIMAP','EXE_SAMTOOLS','EXE_BEDTOOLS','EXE_WFMASH','EXE_GFFREAD',
  'EXE_COLLINEAR');

my ($newDIR,$output_mask,$pancore_mask,$include_file,%included_input_files,%opts) = ('','','',0);
my ($dowfmash,$reference_string) = (0,0);
my ($onlywga,$inputDIR,$cluster_list_file) = (0);
my ($min_overlap,$min_map_qual,$split_chr_regex) = ($MINOVERLAP,$MINQUAL,'');
my ($n_of_cpus,$do_soft) = (4,0);
my ($do_ANIb_matrix,$do_POCP_matrix) = (0,0);
my ($samtools_path,$bedtools_path) = ('','');

my ($min_cluster_size,$runmode,$do_genome_composition,$saveRAM,$ANIb_matrix_file,$POCP_matrix_file);
#my ($evalue_cutoff,$pi_cutoff,$pmatch_cutoff) = 
#  ($BLAST_PVALUE_CUTOFF_DEFAULT,$PERCENT_IDENTITY_CUTOFF_EST_DEFAULT,$PERCENT_MATCH_CUTOFF_DEFAULT);

my $random_number_generator_seed = 0;
my $pwd = getcwd(); $pwd .= '/';

getopts('hvcsPoAzWn:m:d:r:t:I:C:R:B:S:O:Q:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0))
{
  print   "\nusage: $0 [options]\n\n";
  print   "-h this message\n";
  print   "-v print version, credits and checks installation\n";
  print   "-d directory with input files (.fna & .gff pairs)        ";
  print   "   (new files can be added later, creates \n".
    "                                                           ".
    "  output folder named 'directory_pangenes')\n";
  print   "\nOptional parameters:\n";
  print   "-o only run WGA jobs and exit\n";
  print   "-c report pangene set growth analysis                       ".
    "(follows order in -I file if enforced,\n".
    "                                                               ".
    "with -t N skips clusters occup<N)\n";
  print   "-R set random seed for genome composition analysis          ".
    "(optional, requires -c, example -R 1234)\n";
  print   "-m runmode [local|cluster|dryrun]                           ".
    "(default: -m local)\n";
  print   "-n nb of threads in 'local' runmode                         ".
    "(default=$n_of_cpus)\n";
  print   "-I file with .fna files in -d to be included                ".
    "(takes all by default, requires -d)\n";
  print   "-B path to bedtools binary                                  ".
    "(optional, default: -B bedtools)\n";
 
  print   "\nAlgorithms instead of default minimap2 (Mmap):\n";
  print   "-W use wfmash aligner (Wmsh), requires samtools\n";
  print   "-S path to samtools binary                                  ".
    "(optional, default: -S samtools)\n";
  print   "\nOptions that control alignments:\n";
  print   "-O min overlap of genes                                     ".
    "(optional, range [0-1], default: -O $MINOVERLAP)\n";
  print   "-Q min mapping quality, minimap2 only                       ".
    "(optional, default: -Q $MINQUAL)\n";
  print   "-s split genome in chrs                                     ".
    '(optional, requires regex to match chr names ie: -S \'^\d+$\')'."\n";  

  print   "\nOptions that control clustering:\n";
  print   "-t report sequence clusters including at least t taxa       ".
    "(default: t=numberOfTaxa,\n".
    "                                                             ".
    "t=0 reports all clusters)\n";
  print   "-r reference transcriptome .fna file                        ".
    "(by default takes file with\n".
    "                                                            ".   
    " least annotated genes/features)\n";
  # this would require computing BLASTN all vs all
  #print   "-A calculate average identity of clustered sequences,          ".
  #  "(optional, creates tab-separated matrix,\n";
  #print   " uses blastn results                                           ".
  #  " [OMCL])\n";
  print   "-P calculate percentage of conserved sequences (POCS),      ".
    "(optional, creates tab-separated matrix)\n";
  print   "-z add soft-core to genome composition analysis\n";
    
  print "\n".
    "This program computes whole genome alignments (WGA) to define clusters\n".
    "of collinear, orthologous genes/features annotated in GFF files. Such\n".
    "clusters might help define pan-genes across a pangenome.\n".
    "Several WGA algorithms are available and some parameters are customizable.\n".
    "It is designed to process (in a multicore computer or HPC cluster) files\n".
    "contained in a directory (-d), so that new .fna & .gff files can be added\n".
    "while conserving previous results.\n";


  exit;
}

# read version number from CHANGES.txt
open(CHANGES,"$Bin/CHANGES.txt");
while(<CHANGES>) {
  if(eof && /^(\d+):/){ $VERSION = $1 } 
}
close(CHANGES);

if(defined($opts{'v'})) {

  print "\n$0 version $VERSION\n";
  print "\nPrimary citation:\n\n";
  # ...
  print "\nThis software uses genome mappers, please cite them accordingly:\n";
  print " minimap2 https://doi.org/10.1093/bioinformatics/bty191\n";
  print " wfmash   https://github.com/ekg/wfmash\n";

  # check all binaries and data needed by this program and print diagnostic info
  print check_installed_features(@FEATURES2CHECK);
  exit(0);
}

if(defined($opts{'d'})) {

  $inputDIR = $opts{'d'};
  die "# EXIT : need a valid input directory\n" if(!-e $inputDIR);
  if(basename($inputDIR) =~ /([\+])/)
  {
    die "# EXIT : need a valid input directory name, offending char: '$1'\n";
  }
  if(substr($inputDIR,length($inputDIR)-1,1) eq '/'){ chop $inputDIR }
  $newDIR = $pwd.basename($inputDIR)."_pangenes";

} else{ die "# EXIT : need a -d directory with .fasta/.fna files as input\n"; }

if(defined($opts{'m'})) {
  $runmode = $opts{'m'};
  if($runmode ne 'local' && $runmode ne 'cluster' && $runmode ne 'dryrun'){ 
    $runmode = 'cluster'
  }
} else{ $runmode = 'local'; }

if($runmode eq 'local' && defined($opts{'n'}) && $opts{'n'} > 0) {
  $n_of_cpus = $opts{'n'};
}

if($runmode eq 'cluster') {
  # check whether file 'cluster.conf' exists & parse it
  read_cluster_config( "$Bin/HPC.conf" );

  if(!cluster_is_available())
  {       
    print "# EXIT : cluster is not available, please create/edit file cluster.conf\n";
        print_cluster_config();
    die "# EXIT : or choose runmode -m local\n";
  }
}

if(defined($opts{'o'})){ $onlywga = 1; }

if($opts{'r'}){ $reference_string = $opts{'r'}; }
else{ $reference_string = 0; }

if(defined($opts{'t'}) && $opts{'t'} >= 0)
{
  $min_cluster_size = $opts{'t'};
  $output_mask .= $min_cluster_size."taxa_";
  $pancore_mask .= '_'.$min_cluster_size."taxa";
}
else{ $min_cluster_size = 'all'; $output_mask .= "alltaxa_"; }

if($opts{'I'} && $inputDIR)
{
  $include_file = $opts{'I'};
  $output_mask .= basename($include_file)."_";
  $pancore_mask = "_".basename($include_file);
}

check_installed_features(@FEATURES2CHECK);

#if(defined($opts{'A'})){ $do_ANIb_matrix = 1 }

if(defined($opts{'P'})){ $do_POCP_matrix = 1 }

if(defined($opts{'w'})) {
  if(feature_is_installed('WFMASH'))
  {
    $dowfmash = 1;
    $output_mask .= "algWmsh_";
    $pancore_mask .= "_algWmsh";
  }
  else{ warn_missing_soft('WFMASH') }

} else {
  $output_mask .= "algMmap_";
  $pancore_mask .= "_algMmap";

  if(defined($opts{'Q'})){ $min_map_qual = $opts{'Q'} }
}

if(defined($opts{'c'})) {
  $do_genome_composition = 1;
  
  if(defined($opts{'z'})){ $do_soft = 1 }
  
  if($opts{'R'}) {
    $random_number_generator_seed = $opts{'R'};
  }  
  
  srand($random_number_generator_seed);

} else{ $do_genome_composition = 0 }

#if(defined($opts{'E'}))
#{
#  $evalue_cutoff = $opts{'E'};
#  #if($evalue_cutoff > $MAXEVALUEBLASTSEARCH){ $evalue_cutoff = $MAXEVALUEBLASTSEARCH }
#  $output_mask .= "E$evalue_cutoff\_"; $pancore_mask .= "_E$evalue_cutoff";
#}
#if(defined($opts{'C'}))
#{
#  $pmatch_cutoff = $opts{'C'}; # BDBH|OMCL
#  if($pmatch_cutoff < 1){ $pmatch_cutoff = 1 }
#  elsif($pmatch_cutoff > 100){ $pmatch_cutoff = 100 }
#  $output_mask .= "C$pmatch_cutoff\_"; $pancore_mask .= "_C$pmatch_cutoff";
#}

#if(defined($opts{'S'}))
#{
#  $pi_cutoff = $opts{'S'};
#  if($pi_cutoff < 1){ $pi_cutoff = 1 }
#  elsif($pi_cutoff > 100){ $pi_cutoff = 100 }
#  $output_mask .= "S$pi_cutoff\_"; $pancore_mask .= "_S$pi_cutoff";
#}

print "# $0 -d $inputDIR -o $onlywga -r $reference_string ".
  "-t $min_cluster_size -c $do_genome_composition -z $do_soft -I $include_file -m $runmode -n $n_of_cpus -w $dowfmash ".
  "-O $min_overlap -Q $min_map_qual ".
  "-R $random_number_generator_seed -P $do_POCP_matrix\n\n";

if($runmode eq 'cluster')
{
  print "# computer cluster settings\n";
  print_cluster_config();
}

###############################################################

