#!/usr/bin/env perl

# This script assigns stable IDs to clusters and transposed pangene matrices made by get_pangenes.pl
# TODO: If a reference set of pangenes is passed, these will guide the nomenclature, need rules

# Copyright [2023]
# EMBL-European Bioinformatics Institute & Estacion Experimental Aula Dei-CSIC

$|=1;

use strict;
use warnings;
use Getopt::Std;
use File::Basename;
use File::Copy;
use FindBin '$Bin';
use lib "$Bin/lib";
use pangeneTools qw( parse_sequence_FASTA_file );

my $GENEIDIGITS = 6;

my ($INP_dir, $INP_refdir, $INP_outdir, $INP_verbose) = ('','','',0);
my ($cluster_folder, $cluster_folder_name, $cluster_list_file ) = ( '','','' );
my (%opts, @matrix_files);

getopts('hvcd:r:o:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0))
{
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-c print credits and checks installation\n";
  print "-d input directory with get_pangenes.pl results (required, example: -d /path/data_pangenes/..._algMmap_)\n";
  print "-o output directory and prefix ID               (required, example: -o Os4530.POR.1)\n";
  print "-r directory with reference pangenes            (optional, these help renaming the input ones)\n";
  print "-v verbose                                      (optional)\n";
  exit(0);
}

if(defined($opts{'c'})) {
  print "\nPrimary citation:\n https://doi.org/10.1186/s13059-023-03071-z\n";
  exit(0);
}

if(defined($opts{'d'})) { 
  $INP_dir = $opts{'d'};

  if(defined($opts{'r'})){
    $INP_refdir = $opts{'r'};
  }
}
else{ die "# EXIT : need -d directory\n" }

if(defined($opts{'o'})){
  $INP_outdir = $opts{'o'};

  my @IDparts = split(/\./,basename($INP_outdir));
  if(scalar(@IDparts) != 3) {
    die "# EXIT : need a valid -o directory/prefixID ie [clade].[group].[version]\n"
  }  
} 
else{ die "# EXIT : need -o directory/prefixID\n" }

if(defined($opts{'v'})){
  $INP_verbose = 1;
}

print "# $0 -d $INP_dir -r $INP_refdir -o $INP_outdir -v $INP_verbose\n\n";


# 0) locate .cluster_list file and guess actual folder with clusters
opendir(INPDIR,$INP_dir) ||
  die "# ERROR: cannot list $INP_dir , please check -d argument is a valid folder\n";
my @files = grep {/\.cluster_list/} readdir(INPDIR);
closedir(INPDIR);

if(@files) {
  $cluster_list_file = $files[0];
  $cluster_folder_name = (split(/\.cluster_list/,$cluster_list_file))[0];

  $cluster_list_file = $INP_dir . '/' . $cluster_list_file;
  $cluster_folder = $INP_dir . '/' . $cluster_folder_name;

} else {
  die "# ERROR: cannot locate folder with clusters within $INP_dir\n";
}

# 1) locate matrix files (transposed only .tr) 
opendir(INPDIR,$INP_dir) ||
  die "# ERROR: cannot list $INP_dir , please check -d argument is a valid folder\n";
@matrix_files = grep {/pangene_matrix.*?\.tr/} readdir(INPDIR);
closedir(INPDIR);

if(!@matrix_files) {
  die "# ERROR: cannot locate pangene matrices within $INP_dir\n";
}



# 2) create output directory
if ( -e $INP_outdir ) {
  print "\n# WARNING : folder '$INP_outdir' exists, ".
   "files will be overwritten\n\n";

} else {
  if ( !mkdir($INP_outdir) ) {
    die "# ERROR: cannot create $INP_outdir\n";
  }
}

# 3) parse cluster list file, rename clusters and create copies in output folder

my ($clustname,$size,$taxa,$cdnafile,$cdsfile,$pepfile,$gdnafile);
my ($file, $newfile);
my ($new_clustname, $clust_num, $clust_num1) = ('', 0, 0);	
my (%file2ID);

if(!$INP_refdir) {

  mkdir("$INP_outdir/$cluster_folder_name");

  open(CLUSTLIST,"<",$cluster_list_file) || 
    die "# ERROR: cannot read $cluster_list_file\n";
  while(<CLUSTLIST>) {

    if(/^cluster (\S+) size=(\d+) taxa=(\d+) taxa\(gdna\)=\S+ cdnafile: (\S+) cdsfile: (\S+) pepfile: (\S+) gdnafile: (\S+)/) {
      ($clustname,$size,$taxa,$cdnafile,$cdsfile,$pepfile,$gdnafile) =
        ($1,$2,$3,$4,$5,$6,$7);

      if($taxa > 1) {
        $clust_num++;
        $new_clustname = sprintf("%s.pan0%0${GENEIDIGITS}d", 
          $INP_outdir,$clust_num);

      } else {
        $clust_num1++;	      
        $new_clustname = sprintf("%s.pan1%0${GENEIDIGITS}d",             
          $INP_outdir,$clust_num1);
      }

      $file2ID{ $clustname } = $new_clustname;
      #print "$new_clustname $clustname $taxa\n";

      # copy and rename FASTA files of this cluster 
      foreach $file ($cdnafile,$cdsfile,$pepfile,$gdnafile) {
        next if($file eq 'void');
	
        if(!-s "$cluster_folder/$file") {
          die "# ERROR: cannot find $file , skip copy\n";
        }

        $newfile = $file;
        $newfile =~ s/\Q$clustname\E/$file2ID{$clustname}/;	
        print "$file -> $newfile\n" if($INP_verbose);

        copy( "$cluster_folder/$file", "$INP_outdir/$cluster_folder_name/$newfile");          
      }
    } 
  }
  close(CLUSTLIST);

} else {

  # TODO, steps:
  # i) read ref cluster list
  # ii) read input cluster list
  # iii) match clusters ii -> i and decide appropriate names case-by-case
}


# 1) parse pre-computed cluster files, write sequences to temp file & make GMAP index
#@files = ();
#opendir(CLUSTDIR,$cluster_folder) || 
#  die "# ERROR: cannot list $cluster_folder , please check -d argument is a valid folder\n";
#@files = grep {/$cluster_regex$/} readdir(CLUSTDIR);
#closedir(CLUSTDIR);

#if(scalar(@files) == 0) {
#  die "# ERROR: cannot find any valid clusters in $cluster_folder\n";
#}

#==> test_rice_pangenes/Oryza_nivara_v1chr1_alltaxa_5neigh_algMmap_split_/pangene_matrix_genes.tr.tab <==
#source:/home/contrera/github/plant-scripts/pangenes/test_rice_pangenes/Oryza_nivara_v1chr1_alltaxa_5neigh_algMmap_split_/Oryzanivarav1.chr1	Oryza_nivara_v1.chr1	Oryza_sativa.IRGSP-1.0.chr1	Oryza_indica.ASM465v1.chr1
#chr:1	NA	NA	NA
#gene:ONIVA01G00100	gene:ONIVA01G00100	gene:Os01g0100100	gene:BGIOSGA002569

#==> test_rice_pangenes/Oryza_nivara_v1chr1_alltaxa_5neigh_algMmap_split_/pangene_matrix.tr.bed <==
##1	NA	NA	gene:BGIOSGA002568	1	0	NA	NA	gene:BGIOSGA002568
#1	4848	11824	gene:ONIVA01G00010	1	+	gene:ONIVA01G00010	NA	NA

#==> test_rice_pangenes/Oryza_nivara_v1chr1_alltaxa_5neigh_algMmap_split_/pangene_matrix.tr.tab <==
#source:/home/contrera/github/plant-scripts/pangenes/test_rice_pangenes/Oryza_nivara_v1chr1_alltaxa_5neigh_algMmap_split_/Oryzanivarav1.chr1	Oryza_nivara_v1.chr1	Oryza_sativa.IRGSP-1.0.chr1	Oryza_indica.ASM465v1.chr1
#chr1	NA	NA	NA
#gene:ONIVA01G00100	1	1	1

