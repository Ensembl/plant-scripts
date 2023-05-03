package pangeneTools;
require Exporter;

# Subroutines use by get_pangenes.pl and related scripts.
# Some adapted from https://github.com/eead-csic-compbio/get_homologues
#
# Copyright [2021-23]
# EMBL-European Bioinformatics Institute & Estacion Experimental de Aula Dei-CSIC

@ISA = qw(Exporter);
@EXPORT = qw(
  set_phyTools_env 
  feature_is_installed 
  check_installed_features 
  constructDirectory
  count_GFF_genes
  parse_GFF_regex
  get_string_with_previous_genomes
  parse_sequence_FASTA_file
  read_FAI_regex2hash
  extract_isoforms_FASTA
  calc_median
  calc_mode
  calc_stdev
  N50
  get_outlier_cutoffs

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

  # default output keywords of binaries when corrected installed, 
  # add as required for sub check_installed_features
  'EXE_BEDTOOLS'=>'sage',
  'EXE_SAMTOOLS'=>'Usage',
  'EXE_MINIMAP'=>'Usage',
  'EXE_WFMASH'=>'OPTIONS',
  'EXE_GSALIGN'=> 'Usage',
  'EXE_GFFREAD'=>'Usage',
  'EXE_GMAP' => 'Usage',
  'EXE_CLUSTALO' => 'FATAL',
  'EXE_ALISTAT' => 'Syntax',
  'EXE_COLLINEAR'=>'usage',
  'EXE_CUTSEQUENCES'=>'usage',
  'EXE_CLUSTANALYSIS'=>'ERROR',
  'EXE_GZIP'=>'help',
  'EXE_BZIP2'=>'help',
  'EXE_SORT'=>'Usage'
);

my %ubuntu_packages = (
  'EXE_SAMTOOLS' => 'samtools', 
  'EXE_BEDTOOLS' => 'bedtools'
);

################################################################

set_pangeneTools_env();

##########################################################################################

sub set_pangeneTools_env {

  if( ! defined($ENV{'PANGENE_MISSING_BINARIES'}) ) { $ENV{'PANGENE_MISSING_BINARIES'} = '' }
  if( ! defined($ENV{'PANGENE'}) ) { $ENV{'PANGENE'} =  $PANGENEPATH .'/' }

  # installed in this repo
  if( ! defined($ENV{"EXE_MINIMAP"}) ){ 
    $ENV{"EXE_MINIMAP"} = $ENV{'PANGENE'}.'../lib/minimap2/minimap2' 
  }
  if( ! defined($ENV{"EXE_GFFREAD"}) ){ 
    $ENV{"EXE_GFFREAD"} = $ENV{'PANGENE'}.'bin/gffread/gffread' 
  }
  if( ! defined($ENV{"EXE_WFMASH"}) ){ 
    $ENV{"EXE_WFMASH"} = $ENV{'PANGENE'}.'bin/wfmash/build/bin/wfmash' 
  }
  if( ! defined($ENV{"EXE_GSAPATH"}) ){
    $ENV{"EXE_GSAPATH"} = $ENV{'PANGENE'}.'bin/GSAlign/bin/';
    $ENV{"EXE_GSALIGN"} = $ENV{'EXE_GSAPATH'}.'GSAlign';
  }
  if( ! defined($ENV{"EXE_GMAP"}) ){
    $ENV{"EXE_GMAP"} = $ENV{'PANGENE'}.'bin/gmap/exe/bin/gmap';
  }
  if( ! defined($ENV{"EXE_CLUSTALO"}) ){
    $ENV{"EXE_CLUSTALO"} = $ENV{'PANGENE'}.'bin/clustalo';
  }
  if( ! defined($ENV{"EXE_ALISTAT"}) ){
    $ENV{"EXE_ALISTAT"} = $ENV{'PANGENE'}.'bin/AliStat/alistat';
  }


  # should be pre-installed in most settings
  if( ! defined($ENV{"EXE_SAMTOOLS"}) ){ $ENV{"EXE_SAMTOOLS"} = 'samtools' }
  if( ! defined($ENV{'EXE_BEDTOOLS'}) ){ $ENV{'EXE_BEDTOOLS'} = 'bedtools' }
  if( ! defined($ENV{'EXE_GZIP'}) ){ $ENV{'EXE_GZIP'} = 'gzip' }
  if( ! defined($ENV{'EXE_BZIP2'}) ){ $ENV{'EXE_BZIP2'} = 'bzip2' }
  if( ! defined($ENV{'EXE_SORT'}) ){ $ENV{'EXE_SORT'} = 'sort' }

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
# fills $ENV{"PANGENE_MISSING_BINARIES"} with missing binaries
sub check_installed_features {

  my (@to_be_checked) = @_;
  my ($check_summary,$output) = 
    ("\nChecking required binaries and data sources, ".
      "set in pangeneTools.pm or in command line:\n");

  foreach my $bin (@to_be_checked) {
    $check_summary .= sprintf("%18s : ",$bin);
    if($ENV{$bin}) {
      if($bin eq 'EXE_GZIP' || $bin eq 'EXE_BZIP2'){ 
        $output = `$ENV{$bin} -h 2>&1 ` 
      } elsif($bin eq 'EXE_SORT'){
        $output = `$ENV{$bin} --help`
      }else {
        $output = `$ENV{$bin} 2>&1`;
      }
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
  } elsif($feature eq 'GSALIGN') {
    if($env_missing =~ /GSALIGN/){ return 0 }
  } elsif($feature eq 'GFFREAD') {
    if($env_missing =~ /GFFREAD/){ return 0 }
  } elsif($feature eq 'CLUSTALO') {
    if($env_missing =~ /CLUSTALO/){ return 0 }
  } elsif($feature eq 'ALISTAT') {
    if($env_missing =~ /ALISTAT/){ return 0 }
  } elsif($feature eq 'COLLINEAR') {
    if($env_missing =~ /COLLINEAR/){ return 0 }
  } elsif($feature eq 'CUTSEQUENCES') {
    if($env_missing =~ /CUTSEQUENCES/){ return 0 }
  } elsif($feature eq 'GZIP') {
    if($env_missing =~ /GZIP/){ return 0 }
  } elsif($feature eq 'BZIP2') {
    if($env_missing =~ /BZIP2/){ return 0 }
  } elsif($feature eq 'SORT') {
    if($env_missing =~ /SORT/){ return 0 }
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
# Returns hash reference with number of occurrences of unique matching strings 
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



# Takes:
# i) name string of a FASTA file created by _cut_sequences.pl
# ii optional) boolean to return also hash ref mapping id to production_name 
# and parses the sequences in there. Assumes the following header formats:
# >mrnaid geneid coords [production_name]
# >chr01:11217-12435(+) [oryza_sativa_RAPDB]
# and thus supports the same gene having several associated sequences.
# Returns:
# i)   ref to list with gene ids (usually coord-sorted if parsed from GFF)
# ii)  ref to hash with FASTA strings with genes as keys,
#       might contain 2+ seqs for the same gene id
# iii) ref to hash mapping gene ids to chr coordinates
# iv optional)  ref to hash mapping gene ids to [production_name]
sub parse_sequence_FASTA_file {

  my ( $fname, $add_production_name ) = @_;
  my ( $geneid, $coords, $chr, $start, $end, $strand, $prod_name );
  my ( @geneids, %chr_coords, %fasta, %prod_names );

  open(FASTA,"<",$fname) ||
    die "# ERROR(parse_sequence_FASTA_file): cannot read $fname\n";
  while(<FASTA>) {
    #>transcript:Os01t0100100-01 gene:Os01g0100100 1:2983-10815(+) [Oryza_sativa.IRGSP-1.0.chr1]
    if(/^>\S+\s+(\S+)\s+(\S+)\s+\[([^\]]+)\]/) {
      ($geneid, $coords, $prod_name) = ($1, $2, $3);

      if($coords =~ m/^(\S+)?:(\d+)-(\d+)\(([+-])\)/) {

        ($chr, $start, $end, $strand) = ($1, $2, $3, $4);

        # conserved gene id order and chr coords, take 1st
        if(!$fasta{$geneid}) {
          push(@geneids,$geneid);
          $chr_coords{$geneid} = [$chr, $start, $end, $strand];

          if($add_production_name) {
            $prod_names{$geneid} = $prod_name
          }
        }

        # concat in case this is the second isoform
        $fasta{$geneid} .= $_;
      }

    } elsif(/^>(\S+)\s+\[([^\]]+)\]/) { #>chr01:11217-12435(+) [oryza_sativa]

      # geneid is actually segment interval
      ($geneid, $prod_name) = ($1, $2); 

      if($geneid =~ m/^(\S+)?:(\d+)-(\d+)\(([+-])\)/) {

        ($chr, $start, $end, $strand) = ($1, $2, $3, $4);

        push(@geneids,$geneid);
        $chr_coords{$geneid} = [$chr, $start, $end, $strand];
        if($add_production_name) {
          $prod_names{$geneid} = $prod_name
        }

        $fasta{$geneid} = $_;
      }

    # concat sequence
    } else { 
      $fasta{$geneid} .= $_;
    }
  }
  close(FASTA);
  
  if($add_production_name) {
    return ( \@geneids, \%fasta, \%chr_coords, \%prod_names )
  } else {
    return ( \@geneids, \%fasta, \%chr_coords )
  }
}

# Takes a string with gene sequence in FASTA format, as those in \%fasta produced  
# by sub parse_sequence_FASTA_file, and returns an array with individual
# isoform sequences, if any
sub extract_isoforms_FASTA {
  my ($all_isof_FASTA) = @_;

  my @isoform_seqs;
  my $isoform = '';

  my @seqs = split(/\n/,$all_isof_FASTA);
  foreach my $s (0 .. $#seqs) {
    if($seqs[$s] =~ m/^>/) {
      if($isoform) {
        push(@isoform_seqs, $isoform)
      }
      $isoform = "$seqs[$s]\n"; # header
    } else {
      $isoform .= "$seqs[$s]";
    }
  }

  # last isof
  push(@isoform_seqs, $isoform);

  return @isoform_seqs;
}


# Takes ref to list of numbers and returns the median
sub calc_median {

  my ($dataref) = @_;

  my $mid = int(scalar(@$dataref)/2);
  my @sorted = sort {$a<=>$b} (@$dataref);

  if(scalar(@sorted) % 2) {
    return $sorted[ $mid ]
  }
  else {
    return sprintf("%1.1f",($sorted[$mid-1] + $sorted[$mid])/2)
  }
}

# Takes ref to list of numbers and returns a list with mode(s)
sub calc_mode {

  my ($dataref) = @_;

  my ($max, $elem, @modes, %obs) = (0);

  foreach $elem (@$dataref) {
    $obs{$elem}++;
    if($obs{$elem} > $max) { 
      $max = $obs{$elem}
    }
  }

  foreach $elem (sort {$b<=>$a} keys(%obs)) {
    if($obs{$elem} == $max) {
      push(@modes, $elem);
    }
  }    

  return @modes;
}

# Takes ref to list of numbers and returns the mean
sub calc_mean {

  my ($dataref) = @_;

  my $mean = 0;
  foreach (@$dataref) { $mean += $_ }
  return $mean / scalar(@$dataref);
}

# Takes ref to list of numbers and returns the standard deviation
sub calc_stdev {

  my ($dataref) = @_;
  my $mean = calc_mean($dataref);
  my $sum_squares = 0;

  foreach (@$dataref) { 
    $sum_squares += (($_ - $mean) **2);  
  }

  return sqrt( $sum_squares / (scalar(@$dataref)-1) );
}


# Takes ref to list of numbers and returns N50
sub N50 {

  my ($dataref) = @_;

  my ($total_len,$cumul_len,$N50,$seqlen) = (0,0,-1);
  my @sorted = sort {$b<=>$a} (@$dataref);

  foreach $seqlen (@sorted) {
    $total_len += $seqlen
  }

  foreach $seqlen (@sorted) {
    $cumul_len += $seqlen;
    if($cumul_len>$total_len/2){ 
      $N50 = $seqlen;
      last 
    }
  }
    
  return $N50;
}


# Takes ref to list of numbers and returns median, 
# lower and upper cutoff (scalars) to call outliers:
# i)   median
# ii)  Q1 - 1.5 IQR 
# iii) Q3 + 1.5 IQR
sub get_outlier_cutoffs {

  my ($dataref, $verbose) = @_;

  my @values = sort {$a<=>$b} (@$dataref);

  # 25% percentile (Q1)
  my $Q1 = $values[sprintf("%.0f",(0.25*($#values)))];

  # 50% median
  my $median = $values[sprintf("%.0f",(0.5*($#values)))];

  # 75% percentile (Q3)
  my $Q3 = $values[sprintf("%.0f",(0.75*($#values)))];

  my $IQR = $Q3-$Q1; 

  print "# Q1 $Q1 Q3 $Q3 IQR $IQR\n" if($verbose);

  return (
    $median,
    sprintf("%1.1f", $Q1 - (1.5 * $IQR)),
    sprintf("%1.1f", $Q3 + (1.5 * $IQR)),
  )
}

# Takes 2 strings:
# 1) name of FASTA .fai index file
# 2) (optional) regex to match chromosome names, applied to 1st column
# Returns ref to hash with chr and/or 'unplaced' as keys and BED strings as value
sub read_FAI_regex2hash {

  my ($faifile,$regex) = @_;

  my ($seqname,$size);
  my %bed;

  open(FAI,"<$faifile") ||
    die "# ERROR(read_FAI_regex2hash): cannot read $faifile $!\n";

  while (<FAI>) {
    #1A      602900890       60      60      61
    #1B      697493198       612949359       60      61
    if(/^(\S+)\t(\d+)/) {
      ($seqname, $size) = ($1, $2);
      if($regex && $seqname !~ m/^$regex$/) {
        $bed{'unplaced'} .= "$seqname\t0\t$size\n";
      } else {
        $bed{$seqname} = "$seqname\t0\t$size\n";
      }
    }
  }

  close(FAI);

  return \%bed;
}

1;
