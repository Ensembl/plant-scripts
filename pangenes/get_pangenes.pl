#!/usr/bin/env perl 

# This program computes whole genome alignments (WGA) to define clusters 
# of collinear, orthologous genes/features annotated in GFF files. Such 
# clusters define pan-genes across a pangenome.
# Several WGA algorithms are available and some parameters are customizable. 
# It is designed to process (in a multicore computer or HPC cluster) files 
# contained in a directory (-d), so that new .fna & .gff files can be added 
# while conserving previous results. 
#
# This script calls _cut_sequences.pl, _collinear_genes.pl & _cluster_analysis.pl
# and produces different types of output:
# 1) clusters of CDS (nucl & pep), cDNA and gDNA sequences of collinear genes 
# 2) pangenome matrices that summarize the genome occupancy of clusters
# 3) matrix of % conserved sequences summarizing cDNA shared clusters across genomes
# 4) optionally (-c) matrices with core- and pangene set growth simulations

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
my $NUMCPUS    = 4;
my $MINQUAL    = 50;       # minimap2 mapping quality
my $MINOVERLAP = 0.50;     # used by bedtools to call overlapping features    
my $MINGFFLEN  = 100;      # used by gffread to extract GFF features
my $NOFSAMPLESREPORT = 20; # number of samples while simulating pangene growth

## list of features/binaries required by this program (do not edit)
my @FEATURES2CHECK = (
  'EXE_MINIMAP','EXE_BEDTOOLS','EXE_GFFREAD',
  'EXE_COLLINEAR','EXE_CUTSEQUENCES','EXE_CLUSTANALYSIS',
  'EXE_GZIP', 'EXE_BZIP2'
);

my (%opts,%included_input_files);
my ($newDIR,$output_mask,$pancore_mask,$include_file) = ('','','',0);
my ($dowfmash,$dogsalign,$reference_string) = (0,0,0);
my ($onlywga,$inputDIR,$alg) = (0);
my ($min_overlap,$min_map_qual,$split_chr_regex) = ($MINOVERLAP,$MINQUAL,'');
my ($n_of_cpus,$do_soft, $highly_repetitive) = ($NUMCPUS,0,0);
my ($bedtools_path,$samtools_path,$wfmash_path,$gsalign_path) = ('','','','');
my ($min_cluster_size,$runmode,$do_genome_composition);

my $random_number_generator_seed = 0;
my $pwd = getcwd(); 
$pwd .= '/';

getopts('hvcoAHzgwG:W:n:m:d:r:t:I:C:R:B:S:O:Q:s:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)) {

  print   "\nusage: $0 [options]\n\n";
  print   "-h this message\n";
  print   "-v print version, credits and checks installation\n";
  print   "-d directory with input files (.fna & .gff pairs)        ";
  print   "   (new files can be added later, creates \n".
    "                                                           ".
    "  output dir named 'directory_pangenes')\n";
  print   "\nOptional parameters:\n";
  print   "-o only run WGA jobs and exit\n";
  print   "-c report pangene set growth analysis                       ".
  #  "(follows order in -I file if enforced,\n".
  #  "                                                               ".
    "(with -t N skips clusters occup<N)\n";
  print   "-R set random seed for pangene set growth analysis          ".
    "(optional, requires -c, example -R 1234)\n";
  print   "-m runmode [local|cluster|dryrun]                           ".
    "(default: -m local)\n";
  print   "-n nb of threads                                            ".
    "(default=$n_of_cpus)\n";
  print   "-I file with filenames in -d to be included                 ".
    "(takes all by default)\n";
  print   "-B path to bedtools binary                                  ".
    "(optional, default: -B bedtools)\n";
 
  print   "\nAlgorithms instead of default minimap2 (Mmap):\n";
  print   "-w use wfmash aligner (Wmsh)\n";
  print   "-W path to wfmash binary                                    ".
    "(optional, default: -wf wfmash)\n";
  print   "-g use GSAlign aligner (GSal)                               ".
    "(optional, produces also ANI matrix)\n";
  print   "-G path to GSalign bin/ folder                              ".
    "(optional)\n";
  print   "-S path to samtools binary, required with -w -s             ".
    "(optional, default: -S samtools)\n";
  print   "\nOptions that control alignments:\n";
  print   "-O min overlap of genes                                     ".
    "(optional, range [0-1], default: -O $MINOVERLAP)\n";
  print   "-Q min mapping quality, minimap2 only                       ".
    "(optional, default: -Q $MINQUAL)\n";
  print   "-s split genome in chrs, align only homologous chrs         ".
    "(optional, requires regex to match chr names\n".
    "                                                            ".
    ' ie: -s \'^\d+$\' , remove tmp/ if new regex used)'."\n";
 
  print   "-H genome is highly repetitive                              ".
    "(optional, <minimap RAM, masks intergenes)\n"; 
  print   "\nOptions that control clustering:\n";
  print   "-t report sequence clusters including at least t taxa       ".
    "(default: t=numberOfTaxa,\n".
    "                                                             ".
    "t=0 reports all clusters)\n";
  print   "-r reference genome .fna file                               ".
    "(by default takes file with\n".
    "                                                            ".   
    " least annotated genes/features)\n";
  #print   "-z add soft-core to genome composition analysis\n";
    
  print "\n".
    "This program computes whole genome alignments (WGA) to define clusters\n".
    "of collinear, orthologous genes/features annotated in GFF files. Such\n".
    "clusters might help define pan-genes across a pangenome.\n".
    "Several WGA algorithms are available and some parameters are customizable.\n".
    "It is designed to process (in a multicore computer or HPC cluster) files\n".
    "contained in a directory (-d), so that new .fna & .gff files can be added\n".
    "while conserving previous results. Produces different types of output:\n\n".
    " 1) clusters of CDS (nucl & pep), cDNA and gDNA sequences of collinear genes\n".
    " 2) pangenome matrices that summarize the genome occupancy of clusters\n".
    " 3) matrix of % conserved sequences summarizing shared cDNA clusters across genomes\n".
    " 4) optionally (-c) matrices with core- and pan-gene set growth simulations\n";
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
  print "\nPrimary citation:\n https://github.com/Ensembl/plant-scripts/pangenes\n";
  print "\nThis software uses external algorithms, please cite them accordingly:\n";
  print " minimap2 https://doi.org/10.1093/bioinformatics/bty191\n";
  print " wfmash   https://github.com/ekg/wfmash\n";
  print " GSAlign  https://doi.org/10.1186/s12864-020-6569-1\n";

  # check all binaries needed by this program and print diagnostic info
  print check_installed_features(@FEATURES2CHECK);
  exit(0);
}

if(defined($opts{'d'})) {

  $inputDIR = $opts{'d'};
  die "# EXIT : need a valid input directory\n" if(!-e $inputDIR);
  if(basename($inputDIR) =~ /([\+])/) {
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
} else{ $runmode = 'local' }

if($runmode eq 'local' && defined($opts{'n'}) && $opts{'n'} > 0) {
  $n_of_cpus = $opts{'n'};
}

if($runmode eq 'cluster') {
  # check whether file 'cluster.conf' exists & parse it
  read_cluster_config( "$Bin/HPC.conf" );

  if(!cluster_is_available()) {       
    print "# EXIT : cluster is not available, please create/edit file cluster.conf\n";
    print_cluster_config();
    die "# EXIT : or choose runmode -m local\n";
  }
}

if(defined($opts{'H'})) { 
  $highly_repetitive = 1;
  $output_mask .= "highrep_";
  $pancore_mask .= '_highrep';
}

if(defined($opts{'o'})){ $onlywga = 1 }

if($opts{'r'}){ $reference_string = $opts{'r'} }
else{ $reference_string = 0 }

if(defined($opts{'t'}) && $opts{'t'} >= 0) {
  $min_cluster_size = $opts{'t'};
  $output_mask .= $min_cluster_size."taxa_";
  $pancore_mask .= '_'.$min_cluster_size."taxa";
} else { 
  $min_cluster_size = 'all'; 
  $output_mask .= "alltaxa_"; 
}

if($opts{'I'} && $inputDIR) {
  $include_file = $opts{'I'};
  $output_mask .= basename($include_file)."_";
  $pancore_mask = "_".basename($include_file);
}

check_installed_features(@FEATURES2CHECK);

if(defined($opts{'w'})) {

  if(defined($opts{'W'})) {
    $wfmash_path = $opts{'W'};
    $ENV{"EXE_WFMASH"} = $wfmash_path;
  }

  check_installed_features('EXE_WFMASH');
  if(feature_is_installed('WFMASH')) {
    $alg = 'Wmsh';
    $dowfmash = 1;
    $output_mask .= "alg$alg\_";
    $pancore_mask .= "_alg$alg";
  }
  else{ 
    die "# EXIT : cannot find wfmash binary, ".
      "see dependency instructions or set path with -W\n";
  }
} elsif(defined($opts{'g'})) { 

  if(defined($opts{'G'})) {
    $gsalign_path = $opts{'W'};
    $ENV{"EXE_GSAPATH"} = $gsalign_path;
  }
  
  check_installed_features('EXE_GSALIGN');
  if(feature_is_installed('GSALIGN')) {
    $alg = 'GSal';
    $dogsalign = 1;
    $output_mask .= "alg$alg\_";
    $pancore_mask .= "_alg$alg";
  }
  else{
    die "# EXIT : cannot find GSAlign binary, ".
      "see dependency instructions or set path to GSAlign/bin/ with -G\n";
  }

} else {
  $alg = 'Mmap';
  $output_mask .= "alg$alg\_";
  $pancore_mask .= "_alg$alg";

  if(defined($opts{'Q'})) { 
    $min_map_qual = $opts{'Q'}; 
    if($min_map_qual < 0){ $min_map_qual = 0 }
    $output_mask .= "Q$min_map_qual\_"; 
    $pancore_mask .= "_Q$min_map_qual";
  }
}

if(defined($opts{'S'})) {
  $samtools_path = $opts{'S'};
  $ENV{"EXE_SAMTOOLS"} = $samtools_path;
}

if(defined($opts{'s'}) && $opts{'s'} ne '') {

  if($alg eq 'Wmsh') {
    check_installed_features('EXE_SAMTOOLS');
    if(!feature_is_installed('SAMTOOLS')) {
      print "# EXIT : cannot find samtools binary, ".
        "set path with -S\n";
    } 
  }

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

# to test other minimap versions
#$ENV{"EXE_MINIMAP"} = '~/soft/minimap2-2.24_x64-linux/minimap2';

print "# $0 -d $inputDIR -o $onlywga -r $reference_string ".
  "-t $min_cluster_size -c $do_genome_composition -z $do_soft -I $include_file ".
  "-m $runmode -w $dowfmash -g $dogsalign -O $min_overlap -Q $min_map_qual ".
  "-s '$split_chr_regex' -H $highly_repetitive ".
  "-W '$wfmash_path' -G '$gsalign_path' -B '$bedtools_path' -S '$samtools_path' ".
  "-n $n_of_cpus -R $random_number_generator_seed\n\n";

if($runmode eq 'cluster') {
  print "# computer cluster settings\n";
  print_cluster_config();
}

###############################################################

## 0) declare most important variables 

my ($total_dry, $refOK, $total_genes, $n_of_taxa) = (0,0,0,0);
my ($min_geneome_size, $reference_proteome, $infile, $command);
my ($order, $taxon, $taxon2, $previous_files, $current_files);
my (@taxa);

constructDirectory($newDIR) || 
  die "# EXIT : cannot create directory $newDIR , check permissions\n";

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
    die "# EXIT : cannot run another instance of the program with " .
      "same input data while previous is running\n";
}

# 0.2) open important files
my $input_order_file = $newDIR."/input_order.txt";
my $dryrun_file = $newDIR."/dryrun.txt";

print "# version $VERSION\n";
print "# results_directory=$newDIR\n";
print "# parameters: MINGFFLEN=$MINGFFLEN\n";

# 0.3) prepare dryrun file if required
if($runmode eq 'dryrun') {
    open(DRYRUNLOG,">",$dryrun_file);
}

## 1) read all input files, identify formats and generate required files 
##    in temporary directory

print "\n# checking input files...\n";
$min_geneome_size = -1;
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
my ($dnafile,$gffile,$plain_dnafile,$plain_gffile,$num_genes,$Mb);
my ($outcDNA,$outCDS,$outpep,$clusteroutfile,$tx1,$tx2);
my (%cluster_PIDs,%ngenes,@gff_outfiles,@gff_logfiles,@to_be_deleted);

foreach $infile (@inputfiles) {

  # save this taxon
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
      system("$ENV{'EXE_GZIP'} -dc $dnafile > $plain_dnafile")     
    } else {
      cp($dnafile,$plain_dnafile)
    }
  } else {
    print "# re-using $plain_dnafile\n"
  }

  if(!-s $plain_gffile) {
    if($gffile =~ m/\.gz/) {
      print "# uncompressing $gffile\n";
      system("$ENV{'EXE_GZIP'} -dc $gffile > $plain_gffile")     
    } else {
      cp($gffile,$plain_gffile)
    }
  } else {
    print "# re-using $plain_gffile\n"
  }

  # work out sequence stats and make sure split regex works
  $Mb = (-s $plain_dnafile) / (1024 * 1024);
  $num_genes = count_GFF_genes( $plain_gffile );
  $ngenes{$taxon} = $num_genes;
  $total_genes += $num_genes;

  if($split_chr_regex) {
    my $ref_parsed_chrs = parse_GFF_regex($plain_gffile, $split_chr_regex, 0);
    printf("# %s %1.2fMB genes=%d chrs/contigs=%d\n",
      $dnafile,
      $Mb,
      $num_genes,
      scalar(keys(%$ref_parsed_chrs)));
    if(scalar(keys(%$ref_parsed_chrs)) < 1) {
      die "# ERROR: regex '$split_chr_regex' does not match chr names in $gffile, please edit\n"
    }
  } else {
    printf("# %s %1.2fMB genes=%d\n",
      $dnafile,
      $Mb,
      $num_genes);
  }

  # extract cDNA and CDS sequences
  # Note: this also creates a FASTA index file that might be used by wfmash later on
  if(-s $outpep && -s $outCDS && -s $outcDNA) {
    #print "# re-using $outCDS\n";
    
  } else {
    $command = "$ENV{'EXE_CUTSEQUENCES'} -sp $taxon -fa $plain_dnafile ".
      "-gf $plain_gffile -p $ENV{'EXE_GFFREAD'} -l $MINGFFLEN -o $newDIR"; #die $command; 

    if($runmode eq 'cluster') {
      unlink($clusteroutfile);
      submit_cluster_job("cut$infile",$command,$clusteroutfile,$newDIR,\%cluster_PIDs);
    } elsif($runmode eq 'dryrun') {
          $command =~ s/\\//g;
          print DRYRUNLOG "$command > $clusteroutfile\n";
          $total_dry++;
    } else { # 'local' runmode
      $command = "$command > $clusteroutfile"; 
      system("$command");
      if($? != 0) {
        die "# EXIT: failed while extracting GFF features ($command)\n";
      }
    }

    push(@gff_logfiles, $clusteroutfile);
  } 
}

# size of gene clusters
if($min_cluster_size eq 'all'){ 
  $min_cluster_size = $n_of_taxa 
} 

print "\n# $n_of_taxa genomes, $total_genes genes\n\n";

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

# check errors in logfiles 
foreach $gffile (@gff_logfiles) {
  open(LOG,"<",$gffile) || die "# ERROR: cannot open $gffile\n";
  while(<LOG>) {
    if(/^# ERROR/) {
      print;
      exit;
    }
  }
  close(LOG);
}

print "# done\n\n";

# 1.4) correct list of included files
if($include_file) {

  my ($included,$includedfull,$n_of_matched_included,@Itaxa);

  open(INCL,$include_file) || die "# EXIT : cannot read $include_file\n";
  while(<INCL>) {
    next if(/^#/ || /^$/);
    $included_input_files{(split)[0]} = $.;
  }
  close(INCL);

  print "# included input files (".scalar(keys(%included_input_files))."):\n";

  $refOK = $total_genes = $n_of_matched_included = 0;
  TAXON: foreach $included 
    (sort {$included_input_files{$a}<=>$included_input_files{$b}}
    keys(%included_input_files)) {
    foreach $taxon (@taxa) {
      if($taxon =~ /^$included$/) {
        push(@Itaxa,$taxon);
        $total_genes += $ngenes{$taxon};
        print ": $taxon $ngenes{$taxon}\n";
        $n_of_matched_included++;
        next TAXON;
      }
    }
  }
  print "\n";

  if($n_of_matched_included < scalar(keys(%included_input_files)))
  {
    die "# EXIT : failed to match taxa included in $include_file ($n_of_matched_included), ".
      "please make sure their names match those of input files\n";
  }

  # update @taxa, $min_cluster_size
  @taxa = @Itaxa;
  $n_of_taxa = scalar(@Itaxa);
  if($min_cluster_size eq 'all') { 
    $min_cluster_size = $n_of_taxa 
  }
  elsif($min_cluster_size > $n_of_taxa) { 
    $min_cluster_size = $n_of_taxa
  }
}

printf("# taxa considered = %d genes = %d\n\n",$n_of_taxa,$total_genes);

if($n_of_taxa<2){ die "# EXIT: need at least two taxa\n" }

# 1.5) set reference proteome index and mask, 
# by default takes genomes with least genes)

my ($geneome_size,$smallest_geneome,$smallest_geneome_name);
my ($reference_genome, $reference_name, $reference_prefix);

for($taxon=0;$taxon<scalar(@taxa);$taxon++) {

  $geneome_size = $ngenes{$taxa[$taxon]};

  # update minimal proteome size
  if((defined($min_geneome_size) && $min_geneome_size== -1) ||
    (defined($min_geneome_size) && defined($geneome_size) 
    && $geneome_size < $min_geneome_size)) {

    $min_geneome_size = $geneome_size;
    $smallest_geneome_name = $taxa[$taxon];
    $smallest_geneome = $taxon;
  }

  # check user-defined reference proteome
  if($reference_string ne '0' && 
    !$refOK && $taxa[$taxon] =~ /$reference_string/) {
    $reference_genome = $taxon;
    $reference_name = $taxa[$taxon];
    $refOK=1;
  }
}

if(!$refOK && $reference_string ne '0') {
  print "# WARNING: cannot find reference genome ($reference_string),".
    " taking default\n";
}

if($reference_string eq '' || !$refOK) {
  $reference_name = $smallest_geneome_name;
  $reference_genome = $smallest_geneome;
}

$reference_prefix = $reference_name;
$reference_prefix =~ s/[\W+]//g;
$output_mask = $reference_prefix."\_" . $output_mask;

print "# mask=$output_mask ($pancore_mask)\n\n" if(!$onlywga);

# 1.6) check previously processed input files to decide whether new WGAs
# are needed and update $pangeneTools::selected_genomes_file
$current_files = join('',sort(@taxa));
if(-s $selected_genomes_file) {
  $previous_files = get_string_with_previous_genomes($selected_genomes_file);
}

open(SEL,">$selected_genomes_file") || die "# cannot create $selected_genomes_file\n";
foreach $taxon (@taxa){ print SEL "$taxon\n" }
close(SEL);


## 2) compute pairwise whole-genome alignments (WGA) and call collinear genes

my ($outTSVfile,$outANIfile,@tmp_wga_output_files,%ANIfiles);

# remove previous merged results to make sure they are updated
unlink($merged_tsv_file);

print "\n# indexing genomes ...\n";
foreach $tx1 (0 .. $#taxa) {

  $taxon = $taxa[$tx1];
  next if($include_file && !$included_input_files{$taxon});

  $clusteroutfile = $newDIR . "/$taxon.index.queue"; 

  $command = "$ENV{'EXE_COLLINEAR'} -t $n_of_cpus ".
    "-sp1 $taxon ".
    "-fa1 $newDIR/_$taxon.fna -gf1 $newDIR/_$taxon.gff ".
    "-sp2 $taxon ".
    "-fa2 $newDIR/_$taxon.fna -gf2 $newDIR/_$taxon.gff ".
    "-T $TMP_DIR -i -r ";
    
  if($split_chr_regex) {
    $command .= "-s '$split_chr_regex' ";
  }
  if($dowfmash) {
    $command .= "-wf -W $ENV{'EXE_WFMASH'} ";
  } elsif($dogsalign) {
    $command .= "-gs -G $ENV{'EXE_GSAPATH'} ";
  } else {
    $command .= "-M $ENV{'EXE_MINIMAP'} "
  }
  
  if($highly_repetitive) {
      $command .= '-H '
  } #print "$taxon $command\n";

  if($runmode eq 'cluster') {
    unlink($clusteroutfile);
    submit_cluster_job("idx$taxon",$command,$clusteroutfile,$newDIR,\%cluster_PIDs);
  } elsif($runmode eq 'dryrun') {
    $command =~ s/\\//g;
    print DRYRUNLOG "$command > $clusteroutfile\n";
    $total_dry++;
  } else { # 'local' runmode
    $command = "$command > $clusteroutfile";
    system("$command");
    if($? != 0) {
      die "# EXIT: failed while indexing genomes ($command)\n";
    }
  }
}

# wait until cluster jobs are done
if($runmode eq 'cluster') {
  check_cluster_jobs($newDIR,\%cluster_PIDs);

} elsif($runmode eq 'dryrun' && $total_dry > 0) {
  close(DRYRUNLOG);
  print "# EXIT: check the list of pending commands at $dryrun_file\n";
  exit;
}

print "# done\n\n"; 

print "\n# running pairwise genome alignments ...\n";
foreach $tx1 (0 .. $#taxa-1) {

  $taxon = $taxa[$tx1];
  next if($include_file && !$included_input_files{$taxon});

  foreach $tx2 ($tx1+1 .. $#taxa) {

    $taxon2 = $taxa[$tx2];
    next if($include_file && !$included_input_files{$taxon2});
    #print $taxon.$taxon2."\n";

    $outTSVfile = $newDIR ."/_$taxon.$taxon2.alg$alg.overlap$min_overlap";
    $outANIfile = $newDIR ."/_$taxon.$taxon2.alg$alg";
    if($split_chr_regex) {
      $outTSVfile .= '.split';
      $outANIfile .= '.split'
    }
    if($highly_repetitive) {
      $outTSVfile .= '.highrep';
      $outANIfile .= '.highrep'
    }
    $outTSVfile .= '.tsv';
    $outANIfile .= '.ANI.tsv';

    # skip job if already computed 
    if(-s $outTSVfile) {
      push(@tmp_wga_output_files,$outTSVfile);
      $ANIfiles{$taxon}{$taxon2} = $outANIfile;
      next;
    }

    $clusteroutfile = $outTSVfile.'.queue';

    $command = "$ENV{'EXE_COLLINEAR'} -ovl $min_overlap -t $n_of_cpus ".
        "-sp1 $taxon ".
        "-fa1 $newDIR/_$taxon.fna -gf1 $newDIR/_$taxon.gff ".
        "-sp2 $taxon2 ".
        "-fa2 $newDIR/_$taxon2.fna -gf2 $newDIR/_$taxon2.gff ".
        "-B $ENV{'EXE_BEDTOOLS'} -T $TMP_DIR ".
        "-add -out $outTSVfile -r "; # reuse index & tmp files

    if($split_chr_regex) {
      $command .= "-s '$split_chr_regex' ";
    } 
    if($dowfmash) {
      $command .= "-wf -W $ENV{'EXE_WFMASH'} ";
    } elsif($dogsalign) {
      $command .= "-gs -G $ENV{'EXE_GSAPATH'} -ANI $outANIfile "
    } else {
      $command .= "-M $ENV{'EXE_MINIMAP'} "
    }
    if($highly_repetitive) {
        $command .= '-H '
    } #die $command;

    if($runmode eq 'cluster') {
      unlink($clusteroutfile);
      submit_cluster_job($taxon.$taxon2,$command,$clusteroutfile,$newDIR,\%cluster_PIDs);
    } elsif($runmode eq 'dryrun') {
      $command =~ s/\\//g;
      print DRYRUNLOG "$command > $clusteroutfile\n";
      $total_dry++;
    } else { # 'local' runmode
      $command = "$command > $clusteroutfile";
      system("$command");
      if($? != 0) {
        die "# EXIT: failed while doing genome alignment ($command)\n";
      }
    }

    push(@tmp_wga_output_files,$outTSVfile);
    $ANIfiles{$taxon}{$taxon2} = $outANIfile;
  }  
}  

# wait until cluster jobs are done
if($runmode eq 'cluster') {
  check_cluster_jobs($newDIR,\%cluster_PIDs);

} elsif($runmode eq 'dryrun' && $total_dry > 0) {
  close(DRYRUNLOG);
  print "# EXIT: check the list of pending commands at $dryrun_file\n";
  exit;
}

print "# done\n\n";

# concat alignment results to $merged_tsv_file (global var)
if(@tmp_wga_output_files) {
  print "# concatenating WGA results...\n";
  $command = 'cat ';
  foreach $outTSVfile (@tmp_wga_output_files) {
    if(!-e $outTSVfile) {
      die "# EXIT, $outTSVfile does not exist, WGA might have failed ".
        "or hard drive is still writing it (please re-run)\n";
    } else {
      $command .= "$outTSVfile ";
    }
  }

  $command .= "> $merged_tsv_file";
  system("$command"); #print $command;
  if($? != 0) {
    die "# EXIT: failed while concatenating WGA results\n";
  }
}

if($onlywga) {
  print "\n# terminating after WGA (-o)\n";
  exit(0);
}

# parse paired ANI files and produce ANI matrix (-g only)
my %ANI;
if($dogsalign && %ANIfiles) {
  my ($totalsize,$value,$size,@sizes,@ANIs);
  foreach $taxon (keys(%ANIfiles)) {
    foreach $taxon2 (keys(%{ $ANIfiles{$taxon} })) {
  
      # parse ANI values, note with -s there will be 1 value per chr 
      $outANIfile = $ANIfiles{$taxon}{$taxon2}; 
      $totalsize = 0;
      @ANIs = @sizes = ();
      open(ANIPAIR,"<",$outANIfile) ||
        warn "# WARN, $outANIfile does not exist, skip it\n";
      while(<ANIPAIR>) { 
        if(/^\S+\t(\S+)\t(\d+)/) {
          ($value,$size) = ($1, $2);
          push(@ANIs,$value);
          push(@sizes,$size);
          $totalsize += $size;
        }
      } 
      close(ANIPAIR);

      if(scalar(@sizes) > 1) {
        $value = 0;
        foreach my $v (0 .. $#ANIs) {
          $value += ($ANIs[$v] * $sizes[$v]/$totalsize);
        } 
        $value = sprintf("%1.2f",$value); 
      } else {
        $value = $ANIs[0];
      }
 
      $ANI{$taxon}{$taxon2} = $value;
      $ANI{$taxon2}{$taxon} = $value; print "$value\n";
    }
  }
} exit;

## 3) extract clusters of collinear sequences and produce pangene set matrices

print "\n# clustering sequences ";
if($do_genome_composition) {
  print "and simulating pangene set growth "
}
print "...\n";

  $newDIR = $pwd.basename($inputDIR)."_pangenes";

my $outfolder = "$newDIR/$output_mask"; 
$clusteroutfile = $outfolder.'.queue';

$command = "$ENV{'EXE_CLUSTANALYSIS'} ".
  "-T $merged_tsv_file -r $reference_name ".
  "-f $outfolder -t $min_cluster_size -s $newDIR ".
  "-B $ENV{'EXE_BEDTOOLS'} ";

if($do_genome_composition) {
  $command .= "-g $NOFSAMPLESREPORT -R $random_number_generator_seed ";
} 

if( -e $outfolder ) {
  printf("\n# WARNING : folder '%s' exists, clusters might be overwritten\n\n",
    basename($outfolder));
} else {
  if ( !mkdir($outfolder) ) {
    die "# ERROR: cannot create $outfolder\n";
  }
}  

if($runmode eq 'cluster') {
  unlink($clusteroutfile);
  submit_cluster_job("clust$reference_name",$command,$clusteroutfile,$newDIR,\%cluster_PIDs);
} elsif($runmode eq 'dryrun') {
  $command =~ s/\\//g;
  print DRYRUNLOG "$command > $clusteroutfile\n";
  $total_dry++;
} else { # 'local' runmode
  $command = "$command > $clusteroutfile";
  system("$command");
  if($? != 0) {
    die "# EXIT: failed while clustering sequences ($command)\n";
  }
}

# wait until cluster jobs are done
if($runmode eq 'cluster') {
  check_cluster_jobs($newDIR,\%cluster_PIDs);

} elsif($runmode eq 'dryrun' && $total_dry > 0) {
  close(DRYRUNLOG);
  print "# EXIT: check the list of pending commands at $dryrun_file\n";
  exit;
}

print "# done\n\n";

my $printOK = 0;
open(CLUSTERSLOG,'<', $clusteroutfile) ||
  die "# ERROR: cannot read $clusteroutfile\n";
while(<CLUSTERSLOG>) {

  if(/^# number of clusters/){ 
    $printOK = 1 
  }
  print if($printOK);
}
close(CLUSTERSLOG);

# print ANI matrix in output folder, 
# Note: order might be different in POCS matrix
my $ANI_matrix_file = $outfolder . '/ANI.tsv';
open(ANIMATRIX, ">", $ANI_matrix_file)
  || die "# EXIT: cannot create $ANI_matrix_file\n";

print ANIMATRIX "genomes";
foreach $tx1 ( 0 .. $#taxa ) {
    print ANIMATRIX "\t$taxa[$tx1]";
}
print ANIMATRIX "\n";

foreach $tx1 ( 0 .. $#taxa ) {
    print ANIMATRIX "$taxa[$tx1]";
    foreach $tx2 ( 0 .. $#taxa ) {
        if ( $tx1 == $tx2 ) {
            print ANIMATRIX "\t100.00"
        } else {
            if(defined($ANI{$taxa[$tx1]}{$taxa[$tx2]})){
                print ANIMATRIX "\t$ANI{ $taxa[$tx1] }{ $taxa[$tx2] }";
            } else {
                print ANIMATRIX "\tNA";
            }
        }
    }
    print ANIMATRIX "\n";
}
close(ANIMATRIX);

print "\n# Average Nucleotide Identity file = $ANI_matrix_file\n\n";


