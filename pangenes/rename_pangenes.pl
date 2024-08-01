#!/usr/bin/env perl

# This script assigns stable IDs to clusters and transposed pangene matrices made by get_pangenes.pl
# 
# Inspired by proposal previously discussed in AgBioData group:
# [clade].[group].[version].panddddddd
# where:
# clade -> 2-letter for species, 1-letter for genus, followed by NCBI Taxon ID ie
# Os4530 for Oryza sativa
# O4527  for Oryza genus
# group -> unique 3-letter code for group or consortium that made the pan genes
# version -> Unique integer version (starting from 1) referring to the build of the clade by that group
# panddddddd -> 'pan' followed by numerical digits as local pan gene identifier; 
# 0000001 onwards for pan-gene clusters with 2 or more members
# 1000001 for singletons (genes found in only one genome)
#
# TODO: If a reference set of pangenes is passed, these will guide the nomenclature, need rules

# Copyright [2023-24]
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

my ($INP_dir, $INP_refdir, $INP_outdir, $INP_verbose) = ('','','',1);
my ($cluster_folder, $cluster_folder_name, $cluster_list_file ) = ( '','','' );
my (%opts, @matrix_files);

getopts('hScd:r:o:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0))
{
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-c print credits and checks installation\n";
  print "-d input directory with get_pangenes.pl results (required, example: -d /path/data_pangenes/..._algMmap_)\n";
  print "-o output directory and prefix ID               (required, example: -o Os4530.POR.1)\n";
#  print "-r directory with reference pangenes            (optional, these help renaming the input ones)\n";
  print "-S silent                                       (optional, by default prints mappings old -> renamed)\n";
  print 
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

if(defined($opts{'S'})){
  $INP_verbose = 0;
}

printf("# %s -d %s -r %s -o %s -S %d\n\n",
	$0, $INP_dir, $INP_refdir, $INP_outdir, !$INP_verbose);



## 0) locate .cluster_list file and guess actual folder with clusters
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

## 1) locate matrix files (transposed only .tr) 
opendir(INPDIR,$INP_dir) ||
  die "# ERROR: cannot list $INP_dir , please check -d argument is a valid folder\n";
@matrix_files = grep {/pangene_matrix.*?\.tr/} readdir(INPDIR);
closedir(INPDIR);

if(!@matrix_files) {
  die "# ERROR: cannot locate pangene matrices within $INP_dir\n";
}



## 2) create output directory
if ( -e $INP_outdir ) {
  print "\n# WARNING : folder '$INP_outdir' exists, ".
   "files will be overwritten\n\n";

} else {
  if ( !mkdir($INP_outdir) ) {
    die "# ERROR: cannot create $INP_outdir\n";
  }
}

## 3) parse cluster list file, rename clusters and create copies in output folder

my ($clustname,$size,$taxa,$gtaxa,$cdnafile,$cdsfile,$pepfile,$gdnafile);
my ($file, $newfile, $namecol);
my ($new_clustname, $clust_num, $clust_num1) = ('', 0, 0);	
my (%file2ID, @data);

# 3.1 no reference clusters
if(!$INP_refdir) {

  mkdir("$INP_outdir/$cluster_folder_name");

  # 3.1.1 take care of cluster_list file and FASTA cluster files

  open(NEWCLUSTLIST,">","$INP_outdir/$cluster_folder_name.cluster_list") ||
    die "# ERROR: cannot create $INP_outdir/$cluster_folder_name.cluster_list\n";

  open(CLUSTLIST,"<",$cluster_list_file) || 
    die "# ERROR: cannot read $cluster_list_file\n";
  while(<CLUSTLIST>) {

    if(/^cluster (\S+) size=(\d+) taxa=(\d+) taxa\(gdna\)=(\S+) cdnafile: (\S+) cdsfile: (\S+) pepfile: (\S+) gdnafile: (\S+)/) {
      ($clustname,$size,$taxa,$gtaxa,$cdnafile,$cdsfile,$pepfile,$gdnafile) =
        ($1,$2,$3,$4,$5,$6,$7,$8);

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

      # mappings
      print "cluster: $clustname -> $new_clustname\n" if($INP_verbose);

      # add this renamed cluster to new cluster list file
      print NEWCLUSTLIST "cluster $new_clustname size=$size taxa=$taxa taxa(gdna)=$gtaxa";

      $newfile = $cdnafile;
      if($newfile ne 'void') { $newfile =~ s/\Q$clustname\E/$new_clustname/ }
      print NEWCLUSTLIST " cdnafile: $newfile";

      $newfile = $cdsfile;
      if($newfile ne 'void') { $newfile =~ s/\Q$clustname\E/$new_clustname/ }
      print NEWCLUSTLIST " cdsfile: $newfile";

      $newfile = $pepfile;
      if($newfile ne 'void') { $newfile =~ s/\Q$clustname\E/$new_clustname/ }
      print NEWCLUSTLIST " pepfile: $newfile";

      $newfile = $gdnafile;
      if($newfile ne 'void') { $newfile =~ s/\Q$clustname\E/$new_clustname/ }
      print NEWCLUSTLIST " gdnafile: $newfile\n";

      # actually copy and rename FASTA files of this cluster 
      foreach $file ($cdnafile,$cdsfile,$pepfile,$gdnafile) {
        next if($file eq 'void');
	
        if(!-s "$cluster_folder/$file") {
          die "# ERROR: cannot find $file , stop copying\n";
        }

        $newfile = $file;
        $newfile =~ s/\Q$clustname\E/$file2ID{$clustname}/;	

		# mappings
        print "file: $file -> $newfile\n" if($INP_verbose);

        if(!copy( "$cluster_folder/$file", "$INP_outdir/$cluster_folder_name/$newfile")) {
          die "# ERROR: cannot cp $cluster_folder/$file to $INP_outdir/$cluster_folder_name/$newfile, stop copying\n";
        }		
      }
    } else {
      print NEWCLUSTLIST;
    } 
  }
  close(CLUSTLIST);

  close(NEWCLUSTLIST);

  # 3.1.2 take care of (transposed) pangene matrix files
  foreach $file (@matrix_files) {
    
    if($file =~ m/\.bed/) {
      $namecol = 3
    } else {
      $namecol = 0
    }

    open(NEWMATFILE,">","$INP_outdir/$file") ||
      die "# ERROR: cannot create $INP_outdir/$file\n";

    open(MATFILE,"<","$INP_dir/$file") ||
      die "# ERROR: cannot read $INP_dir/$file\n";

    while(<MATFILE>) {
      if(/^source/ || /^chr.+?\tNA/){ 
        print NEWMATFILE

      } else {
        @data = split(/\t/, $_);

        if(!defined($file2ID{ $data[$namecol] })) {
          die "# ERROR: cannot find ID for cluster $data[$namecol] , stop copying $file\n"
        }

	$data[$namecol] = $file2ID{ $data[$namecol] };
	
        print NEWMATFILE join("\t",@data);
      }   
    }
	  
    close(MATFILE);

    close(NEWMATFILE);   
  } 

} else {

  # 3.2

  # TODO, steps:
  # i) read ref cluster list
  # ii) read input cluster list
  # iii) match clusters ii -> i and decide appropriate names case-by-case
}
