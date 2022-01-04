# adapted from https://github.com/eead-csic-compbio/get_homologues 

package pangeneTools;
require Exporter;

@ISA = qw(Exporter);
@EXPORT = qw(
  set_phyTools_env 
  feature_is_installed 
  check_installed_features 
  constructDirectory
  count_GFF_genes
  parse_GFF_regex
  get_string_with_previous_genomes

  $merged_tsv_file
  $lockfile
  $selected_genomes_file
  $TMP_DIR
);

use strict;

our $WORKING_DIR   = '';
our $TMP_DIR       = '';
our $merged_tsv_file       = 'mergedpairs.tsv';
our $lockfile              = '.lock';
our $selected_genomes_file = 'selected.genomes'; # names of genomes previously used 

#################################################################

use FindBin '$Bin';
my $PANGENEPATH = "$Bin";

my %feature_output = (

  # default output of binaries when corrected installed, add as required for
  # check_installed_features
  'EXE_BEDTOOLS'=>'sage',
  'EXE_SAMTOOLS'=>'Usage',
  'EXE_MINIMAP'=>'Usage',
  'EXE_WFMASH'=>'OPTIONS',
  'EXE_GFFREAD'=>'Usage',
  'EXE_COLLINEAR'=>'usage',
  'EXE_CUTSEQUENCES'=>'usage',
  'EXE_CLUSTANALYSIS'=>'ERROR',
  'EXE_GZIP'=>'help'
);

my %ubuntu_packages = (
  'EXE_SAMTOOLS' => 'samtools', 
  'EXE_BEDTOOLS' => 'bedtools'
);

################################################################

set_pangeneTools_env();

##########################################################################################

sub set_pangeneTools_env {

  if( ! defined($ENV{'MARFIL_MISSING_BINARIES'}) ) { $ENV{'MARFIL_MISSING_BINARIES'} = '' }
  if( ! defined($ENV{'PANGENES'}) ) { $ENV{'PANGENE'} =  $PANGENEPATH .'/' }

  # installed in this repo
  if( ! defined($ENV{"EXE_MINIMAP"}) ){ 
    $ENV{"EXE_MINIMAP"} = $ENV{'PANGENE'}.'bin/minimap2-2.17/minimap2' 
  }
  if( ! defined($ENV{"EXE_GFFREAD"}) ){ 
    $ENV{"EXE_GFFREAD"} = $ENV{'PANGENE'}.'bin/gffread-0.12.7.Linux_x86_64/gffread' 
  }
  if( ! defined($ENV{"EXE_WFMASH"}) ){ 
    $ENV{"EXE_WFMASH"} = $ENV{'PANGENE'}.'bin/wfmash/build/bin/wfmash' 
  }

  # should be pre-installed in most settings
  if( ! defined($ENV{"EXE_SAMTOOLS"}) ){ $ENV{"EXE_SAMTOOLS"} = 'samtools' }
  if( ! defined($ENV{'EXE_BEDTOOLS'}) ){ $ENV{'EXE_BEDTOOLS'} = 'bedtools' }
  if( ! defined($ENV{'EXE_ZCAT'}) ){ $ENV{'EXE_ZCAT'} = 'zcat' }

  # scripts from this repo
  if( ! defined($ENV{"EXE_COLLINEAR"}) ){ 
    $ENV{"EXE_COLLINEAR"} = $ENV{'PANGENE'}."_collinear_genes.pl" 
  }
  if( ! defined($ENV{"EXE_CUTSEQUENCES"}) ){ 
    $ENV{"EXE_CUTSEQUENCES"} = $ENV{'PANGENE'}."_cut_sequences.pl" 
  }
  if( ! defined($ENV{"EXE_CLUSTANALYSIS"}) ){
    $ENV{"EXE_CLUSTANALYSIS"} = $ENV{'PANGENE'}."_cluster_analysis.pl"
  }

}

########################################################################################

# Check needed binaries and data sources required by functions here,
# fills $ENV{"MISSING_BINARIES"} with missing binaries
sub check_installed_features {

  my (@to_be_checked) = @_;
  my ($check_summary,$output) = 
    ("\nChecking required binaries and data sources, set in pangeneTools.pm or in command line:\n");

  foreach my $bin (@to_be_checked) {
    $check_summary .= sprintf("%18s : ",$bin);
    if($ENV{$bin}) {
      $output = `$ENV{$bin} 2>&1`;
      if(!$output){ $output = '' }
    }

    if($output =~ /$feature_output{$bin}/) {
      $output = "OK (path:$ENV{$bin})"
    }
    else {

      $ENV{"PANGENE_MISSING_BINARIES"} .= "$bin,";
      $output = " wrong path:$ENV{$bin} or needs to be installed";
      if($ubuntu_packages{$bin}){ 
        $output .= ", ie 'sudo apt install $ubuntu_packages{$bin}'" 
      } 
    }
    $check_summary .= $output."\n";
  }

  return $check_summary;
}

# checks whether a passed software can be used before calling it
sub feature_is_installed {

  my ($feature) = @_;

  my $env_missing = $ENV{"PANGENE_MISSING_BINARIES"};

  if($feature eq 'MINIMAP') {
    if($env_missing =~ /MINIMAP/){ return 0 }
  } elsif($feature eq 'SAMTOOLS') {
    if($env_missing =~ /SAMTOOLS /){ return 0 }
  } elsif($feature eq 'BEDTOOLS') {
    if($env_missing =~ /BEDTOOLS/){ return 0 }
  } elsif($feature eq 'WFMASH') {
    if($env_missing =~ /WFMASH/){ return 0 }
  } elsif($feature eq 'GFFREAD') {
    if($env_missing =~ /GFFREAD/){ return 0 }
  } elsif($feature eq 'COLLINEAR') {
    if($env_missing =~ /COLLINEAR/){ return 0 }
  } elsif($feature eq 'CUTSEQUENCES') {
    if($env_missing =~ /CUTSEQUENCES/){ return 0 }
  } elsif($feature eq 'ZCAT') {
    if($env_missing =~ /ZCAT/){ return 0 }
  }

  return 1;
}

# This subroutine is used to construct the directories to store
# the intermediate files and the final files.
# Arguments: 1 (string) name of desired directory
# Returns:  boolean, 1 if successful, else 0
sub constructDirectory {
  my ($dirname) = @_;

  $WORKING_DIR = $dirname . '/';
  if(!-e $WORKING_DIR) { 
    mkdir($WORKING_DIR) || return 0
  }

  $TMP_DIR = $WORKING_DIR."tmp/";
  if(!-e $TMP_DIR){ mkdir($TMP_DIR); }

  $lockfile              = $TMP_DIR.$lockfile;
  $selected_genomes_file = $TMP_DIR.$selected_genomes_file;
  $merged_tsv_file       = $TMP_DIR.$merged_tsv_file;       
  
  return 1
}

# Takes string with name of GFF file and returns number of genes
sub count_GFF_genes {

  my ($gffile) = @_;

  my $num_genes = 0;
  
  open(GFF, "<", $gffile) ||
    die "# ERROR(count_GFF_genes): cannot read $gffile\n";
  while(<GFF>) {
    my @data = split(/\t/,$_);
    if($data[0] !~ /^#/ && $data[2] eq 'gene'){ $num_genes++ }
  }
  close(GFF);

  return $num_genes
}

# Takes two params:
# 1) string with name of GFF file 
# 2) regular expression
# 3) 0-based column of GFF to be parsed
# Returns hash refeference with number of occurrences of unique matching strings 
sub parse_GFF_regex {

  my ($gffile, $regex, $column) = @_;

  my %count;

  open(GFF, "<", $gffile) ||
    die "# ERROR(parse_GFF_regex): cannot read $gffile\n";
  while(<GFF>) {
    my @data = split(/\t/,$_);
    if($data[$column] =~ m/^($regex)$/){ $count{$1}++ }
  }
  close(GFF);

  return \%count
}


# Check genomes used in previous run stored in $selected_genomes_file
sub get_string_with_previous_genomes {

  my ($selected_genomes_file) = @_;

  my ($previous_genomes,@genomes) = ('');

  open(SEL,$selected_genomes_file) || 
    die "# ERROR(get_string_with_previous_genomes): cannot read $selected_genomes_file\n";
  while(<SEL>) {
    chomp $_;
    $_ =~ s/\s+//g;
    push(@genomes,$_);
  }
  close(SEL);

  $previous_genomes = join('',sort(@genomes));

  return $previous_genomes;
}



1;
