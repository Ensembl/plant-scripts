# adapted from https://github.com/eead-csic-compbio/get_homologues 

package pangeneTools;
require Exporter;

@ISA = qw(Exporter);
@EXPORT = qw(
  set_phyTools_env feature_is_installed check_installed_features warn_missing_soft
);

use strict;

#################################################################

use FindBin '$Bin';
my $PANGENEPATH = "$Bin";

my %feature_output = (

  # default output of binaries when corrected installed, add as required for
  # check_installed_features
  'EXE_BEDTOOLS'=>'Usage',
  'EXE_SAMTOOLS'=>'Usage',
  'EXE_MINIMAP'=>'Usage',
  'EXE_WFMASH'=>'OPTIONS',
  'EXE_GFFREAD'=>'Usage',
  'EXE_COLLINEAR'=>'usage'
);


################################################################

set_pangeneTools_env();

##########################################################################################

sub set_pangeneTools_env {

  if( ! defined($ENV{'MARFIL_MISSING_BINARIES'}) ) { $ENV{'MARFIL_MISSING_BINARIES'} = '' }
  if( ! defined($ENV{'PANGENES'}) ) { $ENV{'PANGENE'} =  $PANGENEPATH .'/' }

  # installed in this repo
  if( ! defined($ENV{"EXE_MINIMAP"}) ){ $ENV{"EXE_MINIMAP"} = $ENV{'PANGENE'}.'bin/minimap2-2.17/minimap2' }
  if( ! defined($ENV{"EXE_GFFREAD"}) ){ $ENV{"EXE_GFFREAD"} = $ENV{'PANGENE'}.'bin/gffread-0.12.7.Linux_x86_64/gffread' }
  if( ! defined($ENV{"EXE_WFMASH"}) ){ $ENV{"EXE_WFMASH"} = $ENV{'PANGENE'}.'bin/wfmash/build/bin/wfmash' }

  # should be pre-installed in most settings
  if( ! defined($ENV{"EXE_SAMTOOLS"}) ){ $ENV{"EXE_SAMTOOLS"} = 'samtools' }
  if( ! defined($ENV{'EXE_BEDTOOLS'}) ){ $ENV{'EXE_BEDTOOLS'} = 'bedtools' }

  # scripts from this repo
  if( ! defined($ENV{"EXE_COLLINEAR"}) ){ $ENV{"EXE_COLLINEAR"} = $ENV{'PANGENE'}."_collinear_genes.pl"; }
}

########################################################################################

# Check needed binaries and data sources required by functions here,
# fills $ENV{"MISSING_BINARIES"} with missing binaries
sub check_installed_features {

  my (@to_be_checked) = @_;
  my ($check_summary,$output) = 
    ("\nChecking required binaries and data sources, all set in pangeneTools.pm :\n");

  foreach my $bin (@to_be_checked) {
    $check_summary .= sprintf("%18s : ",$bin);
    if($ENV{$bin}) {
      $output = `$ENV{$bin} 2>&1`;
      if(!$output){ $output = '' }
    }

    if($output =~ /$feature_output{$bin}/) {
      $output = "OK (path:$ENV{$bin})"
    }
    else
    {
      $ENV{"PANGENE_MISSING_BINARIES"} .= "$bin,";
      $output = " needs to be installed (wrong path:$ENV{$bin}) or correct path set";
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
  }

  return 1;
}



1;
