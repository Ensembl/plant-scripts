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
use File::Temp qw(tempfile);
use File::Basename;
use File::Copy "cp";
use Cwd;
use FindBin '$Bin';
use lib "$Bin/lib";
use HPCluster;
use pangeneTools;

my $VERSION = '1.x';

## global variables that control some algorithmic choices
my $MINQUAL    = 50;       # minimap2 mapping quality
my $MINOVERLAP = 0.50;     # used by bedtools to call overlapping features    
my $MINGFFLEN  = 100;      # used by gffread to extract GFF features
my $NOFSAMPLESREPORT = 20; # number of samples used for the generation of pan/core genomes

## list of features/binaries required by this program (do not edit)
my @FEATURES2CHECK = (
  'EXE_MINIMAP','EXE_SAMTOOLS','EXE_BEDTOOLS','EXE_WFMASH','EXE_GFFREAD',
  'EXE_COLLINEAR','EXE_CUTSEQUENCES',
  'EXE_ZCAT'
);

my ($newDIR,$output_mask,$pancore_mask,$include_file,%included_input_files,%opts) = ('','','',0);
my ($dowfmash,$reference_string) = (0,0);
my ($onlywga,$inputDIR,$cluster_list_file) = (0);
my ($min_overlap,$min_map_qual,$split_chr_regex) = ($MINOVERLAP,$MINQUAL,'');
my ($n_of_cpus,$do_soft) = (4,0);
my ($do_ANIb_matrix,$do_POCP_matrix) = (0,0);
my ($samtools_path,$bedtools_path) = ('','');
my ($min_cluster_size,$runmode,$do_genome_composition,$ANIb_matrix_file,$POCP_matrix_file);

my $random_number_generator_seed = 0;
my $pwd = getcwd(); $pwd .= '/';

getopts('hvcPoAzWn:m:d:r:t:I:C:R:B:S:O:Q:s:', \%opts);

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
  print   "-r reference genome .fna file                               ".
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

if(defined($opts{'W'})) {
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

  if(defined($opts{'Q'})) { 
    $min_map_qual = $opts{'Q'}; 
    if($min_map_qual < 0){ $min_map_qual = 0 }
    $output_mask .= "Q$min_map_qual\_"; 
    $pancore_mask .= "_Q$min_map_qual";
  }
}

if(defined($opts{'s'}) && $opts{'s'} ne '') {
  $split_chr_regex = $opts{'s'};
  $output_mask .= "split_";
  $pancore_mask .= "_split"; 
}

if(defined($opts{'c'})) {
  $do_genome_composition = 1;
  
  if(defined($opts{'z'})){ $do_soft = 1 }
  
  if($opts{'R'}) {
    $random_number_generator_seed = $opts{'R'};
  }  
  
  srand($random_number_generator_seed);

} else{ $do_genome_composition = 0 }

if(defined($opts{'O'})) {
  $min_overlap = $opts{'O'};
  if($min_overlap > 1 || $min_overlap < 0.01){ $min_overlap = $MINOVERLAP }
  $output_mask .= "O$min_overlap\_"; 
  $pancore_mask .= "_O$min_overlap";
}

if(defined($opts{'B'})) {
  $bedtools_path = $opts{'B'};
  $ENV{"EXE_BEDTOOLS"} = $bedtools_path;
}

if(defined($opts{'S'})) {
  $samtools_path = $opts{'S'};
  $ENV{"EXE_SAMTOOLS"} = $samtools_path;
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
  print "\nThis software uses external algorithms, please cite them accordingly:\n";
  print " minimap2 https://doi.org/10.1093/bioinformatics/bty191\n";
  print " wfmash   https://github.com/ekg/wfmash\n";
  #print " gffread  https://f1000research.com/articles/9-304/v2\n";

  # check all binaries and data needed by this program and print diagnostic info
  print check_installed_features(@FEATURES2CHECK);
  exit(0);
}

print "# $0 -d $inputDIR -o $onlywga -r $reference_string ".
  "-t $min_cluster_size -c $do_genome_composition -z $do_soft -I $include_file -m $runmode ".
  "-n $n_of_cpus -W $dowfmash -O $min_overlap -Q $min_map_qual -s '$split_chr_regex' ".
  "-B '$bedtools_path' -S '$samtools_path' ".
  "-R $random_number_generator_seed -P $do_POCP_matrix\n\n";

if($runmode eq 'cluster')
{
  print "# computer cluster settings\n";
  print_cluster_config();
}

###############################################################

## 0) declare most important variables 

my ($total_dry, $refOK, $total_genes, $n_of_taxa) = (0,0,0,0);
my ($min_proteome_size, $reference_proteome, $infile,$command);
my ($order, $taxon, $previous_files, $current_files);
my (@taxa);

#my ($infile,$new_infile,$prot_new_infile,$p2oinfile,$seq,$seqL,$comma_input_files,@newfiles,%ressize,%orth_taxa);
#my ($label,%orthologues,$gene,$orth,$orth2,$para,%inparalogues,%paralogues,$FASTAresultsDIR,$order,$minlog);
#my ($n_of_similar_length_orthologues,$clusterfile,$prot_clusterfile,$previous_files,$current_files,$inpara);
#my ($min_proteome_size,$reference_proteome,$smallest_proteome,$proteome_size,%seq_length,$cluster,$annot);
#my ($smallest_proteome_name,$reference_proteome_name,%psize,$taxon,$taxon2,$n_of_clusters,$n_of_taxa,$n_of_residues);
#my ($pname,$n_of_sequences,$refOK,$genbankOK,$cluster_size,$prot_cluster_size,$prot_cluster);
#my (%orth_registry,%LSE_registry,$LSE_reference,$LSE_t,$redo_inp,$redo_orth,%idclusters,%names,$generef);
#my ($reparse_all,$n_of_similar_length_paralogues,$pfam_annot,$protOK,$n_of_taxa_cluster) = (0);
#my ($BDBHdone,$PARANOIDdone,$orthoMCLdone,$n_of_parsed_lines,$n_of_pfam_parsed_lines) = (0,0,0,0,0);
#my ($diff_BDBH_params,$diff_INP_params,$diff_HOM_params,$diff_OMCL_params,$lockcapableFS,$total_dry) = (0,0,0,0,0,0);
#my ($diff_ISO_params,$redo_iso,$partial_sequences,$isof,%full_length_file,%redundant_isoforms,%total_redundant) = (0);
#my ($total_clustersOK,$clgene,$clorth,%ANIb_matrix,%POCP_matrix,%GIclusters,$clustersOK,%cluster_ids) = (0);

constructDirectory($newDIR) || die "# EXIT : cannot create directory $newDIR , check permissions\n";

# 0.1) make sure there is only 1 instance writing to $newDIR
my $lockcapableFS = 0;
my ($fhtest,$testlockfilename) = tempfile( DIR => $newDIR );
if(flock($fhtest,LOCK_EX|LOCK_NB)) {
  $lockcapableFS = 1;
}
else {
  print "# WARNING : cannot lock files in $newDIR ,\n".
    "# please ensure that no other instance of the program is running at this location\n\n";
}
unlink($testlockfilename);

open(my $fhlock,">$pangeneTools::lockfile") ||
  die "# EXIT : cannot create lockfile $pangeneTools::lockfile\n";
if($lockcapableFS) {
  flock($fhlock, LOCK_EX|LOCK_NB) ||
    die "# EXIT : cannot run another instance of the program with same input data while previous is running\n";
}

# 0.2) open important files
my $input_order_file = $newDIR."/input_order.txt";
my $dryrun_file = $newDIR."/dryrun.txt";

print "# version $VERSION\n";
print "# results_directory=$newDIR\n";

# 0.3) prepare dryrun file if required
if($runmode eq 'dryrun')
{
    open(DRYRUNLOG,">",$dryrun_file);
}

## 1) read all input files, identify formats and generate required files in temporary directory

print "\n# checking input files...\n";
$min_proteome_size = -1;
$reference_proteome = $refOK = $total_genes = $n_of_taxa = 0;
$previous_files = $current_files = '';

# 1.1) open and read directory, a pair of 
# FASTA (.fa .faa .fasta) + GFF (.gff .gff3) files per genome is expected

opendir(DIR,$inputDIR) || die "# EXIT : cannot list $inputDIR\n";
my @inputfiles = sort grep {
  /\.fna$/i || /\.fna\.gz$/i || 
  /\.fa$/i || /\.fa\.gz$/i || 
  /\.fasta$/i || /\.fasta\.gz$/i  
} readdir(DIR);
closedir(DIR); 

# 1.2) sort input files and put new files towards the end of @inputfiles: LILO
if(-s $input_order_file) {
  my (@new_order_input_files,%previous_input_file,$n_of_new_infiles);

  open(ORDER,$input_order_file) || die "# EXIT : cannot read $input_order_file\n";
  while(<ORDER>) {
    chomp;
    ($order,$infile) = split(/\t/);
    if(!-s $inputDIR."/".$infile) { 
      die "# EXIT : cannot find previous input file $infile, please re-run everything\n"; 
    }

    $previous_input_file{$infile} = 1;
    $new_order_input_files[$order] = $infile;
  }
  close(ORDER);

  $n_of_new_infiles=0;
  foreach $infile (@inputfiles) {
    next if($previous_input_file{$infile});

    $new_order_input_files[++$order] = $infile;
    print "# order of new input file $infile = $order\n";
    $n_of_new_infiles++;
  }

  if($n_of_new_infiles){ 
    print "# updating $input_order_file with $n_of_new_infiles new input files\n"; 
  }
  open(ORDER,">$input_order_file") || die "# EXIT : cannot write $input_order_file\n";
  
  $order=0;
  foreach $infile (@new_order_input_files) {
    print ORDER "$order\t$infile\n";
    $order++;
  }
  close(ORDER);

  @inputfiles = @new_order_input_files;
} else {
  $order=0;
  open(ORDER,">$input_order_file") || die "# EXIT : cannot write $input_order_file\n";

  foreach $infile (@inputfiles) {
    print ORDER "$order\t$infile\n";
    $order++;
  }
  close(ORDER);
}

# 1.3) iteratively parse input files
my ($dnafile,$gffile,$plain_dnafile,$plain_gffile,$num_genes);
my ($outcDNA,$outCDS,$outpep,$clusteroutfile);
my (%cluster_PIDs,%ngenes,@gff_outfiles,@todelete);

foreach $infile (@inputfiles) {

  ++$n_of_taxa;

  if($infile =~ m/(\S+?)\.f/){
    $taxon = $1;
    push(@taxa,$taxon);
  }

  # check whether matching GFF file exists (expects same taxon preffix)
  $dnafile = $inputDIR."/$infile";
  $gffile  = $inputDIR."/$taxon";
  $outcDNA = $newDIR."/$taxon.cdna.fna";
  $outCDS  = $newDIR."/$taxon.cds.fna";
  $outpep  = $newDIR."/$taxon.cds.faa";
  push(@gff_outfiles, $outcDNA, $outCDS, $outpep);

  if(-s $gffile.'.gff'){
    $gffile .= '.gff'
  } elsif(-s $gffile.'.gff3'){
    $gffile .= '.gff3'
  } elsif(-s $gffile.'.gff.gz'){
    $gffile .= '.gff.gz'
  } elsif(-s $gffile.'.gff3.gz'){
    $gffile .= '.gff3.gz'
  } else {
    die "ERROR: cannot find matching GFF file for $dnafile\n".
      "A valid filename would be $inputDIR/$taxon.gff\n"
  }

  # make temporary copies of uncompressed FASTA & GFF files
  $plain_dnafile  = $newDIR ."/_$taxon.fna";
  $plain_gffile   = $newDIR ."/_$taxon.gff";
  $clusteroutfile = $newDIR ."/_$infile.queue";

  if(!-s $plain_dnafile) {
    if($dnafile =~ m/\.gz/) {
      print "# uncompressing $dnafile\n";
      system("$ENV{'EXE_ZCAT'} $dnafile > $plain_dnafile")     
    } else {
      cp($dnafile,$plain_dnafile)
    }
  } else {
    print "# re-using $plain_dnafile\n"
  }

  if(!-s $plain_gffile) {
    if($gffile =~ m/\.gz/) {
      print "# uncompressing $gffile\n";
      system("$ENV{'EXE_ZCAT'} $gffile > $plain_gffile")     
    } else {
      cp($gffile,$plain_gffile)
    }
  } else {
    print "# re-using $plain_gffile\n"
  }

  # work out sequence stats
  $num_genes = count_GFF_genes( $plain_gffile );
  $ngenes{$taxon} = $num_genes;
  $total_genes += $num_genes;
  print "# $dnafile ngenes=$num_genes\n";

  # extract cDNA and CDS sequences
  if(-s $outpep && -s $outCDS && -s $outcDNA) {
    #skip if already run
    next;
  } else {
    $command = "$ENV{'EXE_CUTSEQUENCES'} -sp $taxon -fa $plain_dnafile ".
      "-gf $plain_gffile -p $ENV{'EXE_GFFREAD'} -l $MINGFFLEN -o $newDIR";
    #die $command; 
    if($runmode eq 'cluster') {
      submit_cluster_job($infile,$command,$clusteroutfile,$newDIR,\%cluster_PIDs);
    } elsif($runmode eq 'dryrun') {
          $command =~ s/\\//g;
          print DRYRUNLOG "$command\n";
          $total_dry++;
    } else { # 'local' runmode
      $command = "$command > /dev/null"; 
      system("$command");
      if($? != 0) {
        die "# EXIT: failed while extracting GFF features ($command)\n";
      }
    }
  }

  #my $refOK, $n_of_genes

  # update included taxa, 


}

# wait until GFF jobs are done
if($runmode eq 'cluster') {
  check_cluster_jobs($newDIR,\%cluster_PIDs);
} elsif($runmode eq 'dryrun' && $total_dry > 0) {
  close(DRYRUNLOG);
  print "# EXIT: check the list of pending commands at $dryrun_file\n";
  exit;
} 

# confirm outfiles
foreach $gffile (@gff_outfiles) {
  if(!-e $gffile){
    die "# EXIT, $gffile does not exist, GFF job search might failed ".
      "or hard drive is still writing it (please re-run)\n";
  }
}

print "# done\n\n";

## 2) compute whole-genome alignments (WGA) and call collinear genes


