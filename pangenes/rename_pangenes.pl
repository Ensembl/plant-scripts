#!/usr/bin/env perl

# This script assigns stable IDs to clusters and transposed pangene matrices made by get_pangenes.pl
# Dependencies: uses system sort & join
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
# If a reference set of pangenes is passed, these will guide the nomenclature, current rules 
# to go from version n to version n+1:
# * if reference gene membership is unchanged and >=50% of other existing members stay in the same cluster, the same identifier as before is used
# * if cluster did not contain a reference gene in version n, and does not in version n+1, then >=50% members rule is maintained
# * singletons get assigned the codes starting at 1000001 so they don’t interfere with numbering system of others
# * If a reference gene starts as a singleton, but later clusters with other genes, then it will be assigned a non-singleton (0000001 type) identifier

# Copyright [2023-24]
# EMBL-European Bioinformatics Institute & Estacion Experimental Aula Dei-CSIC

$|=1;

use strict;
use warnings;
use Getopt::Std;
use File::Basename;
use File::Copy;
use File::Temp qw/ tempfile /;
use FindBin '$Bin';
use lib "$Bin/lib";
use pangeneTools;

my $JOINBIN      = $ENV{'EXE_JOIN'} || 'join';
my $SORTBIN      = $ENV{'EXE_SORT'} || 'sort';
my $SORTPARS     = "--buffer-size=1G ";
$ENV{'LC_ALL'}   = 'POSIX';


my $GENEIDIGITS = 6;

my ($INP_dir, $INP_refdir, $INP_outdir, $INP_verbose) = ('','','',1);
my ($cluster_folder, $cluster_folder_name, $cluster_list_file, $matfile) = ('','','','');
my ($ref_cluster_folder, $ref_cluster_folder_name, $ref_cluster_list_file, $ref_matfile) = ('','','','');
my ($f,$colr,$col,$row,$rowr,$cmd);
my (%opts, @matrix_files, @ref_matrix_files);

getopts('hScd:r:o:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0))
{
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-c print credits and checks installation\n";
  print "-d input directory with get_pangenes.pl results    (required, example: -d /path/data_pangenes/..._algMmap_)\n";
  print "-o output directory for pangene set and prefix ID  (required, example: -o Os4530.POR.1)\n";
  print "-r folder with previous reference pangene set      (optional, expects common taxa)\n";
  print "-S silent                                          (optional, by default prints mappings old -> renamed)\n";
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
} else {
  foreach $f (@matrix_files) {
    if($f eq 'pangene_matrix_genes.tr.tab') {
      $matfile = $INP_dir . '/' . $f;
	  last; 
    }
  }
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

} else { # 3.2 reference clusters passed to guide the nomenclature of new ones

  my (@ref_taxa, @ref_clusters, @ref_cluster_names);
  my (@taxa, @clusters, @cluster_names, %input2ref, %input2ref_cluster);
  my ($n_common_taxa, $n_ref_clusters, $n_identical_clusters) = (0,0,0);

  ## 3.2.0) locate .cluster_list file and guess actual folder with reference clusters
  opendir(REFDIR,$INP_refdir) ||
    die "# ERROR: cannot list $INP_refdir , please check -r argument is a valid folder made with this script\n";
  my @files = grep {/\.cluster_list/} readdir(REFDIR);
  closedir(REFDIR);

  if(@files) {
    $ref_cluster_list_file = $files[0];
    $ref_cluster_folder_name = (split(/\.cluster_list/,$cluster_list_file))[0];

    $ref_cluster_list_file = $INP_refdir . '/' . $cluster_list_file;
    $ref_cluster_folder = $INP_refdir . '/' . $cluster_folder_name;

  } else {
    die "# ERROR: cannot locate folder with clusters within $INP_refdir\n";
  }

  ## 3.2.1) locate matrix files (transposed only .tr) 
  opendir(REFDIR,$INP_refdir) ||
    die "# ERROR: cannot list $INP_refdir , please check -r argument is a valid folder made with this script\n";
  @ref_matrix_files = grep {/pangene_matrix_genes.tr.tab/} readdir(REFDIR);
  closedir(REFDIR);

  if(!@ref_matrix_files) {
    die "# ERROR: cannot locate pangene matrices within $INP_refdir\n";
  } else {
    $ref_matfile = $INP_refdir . '/' . $ref_matrix_files[0];
  }

  # 3.2.2) read reference pangene_matrix_genes file 
  open(GENEMATR,"<",$ref_matfile) ||
    die "# ERROR: cannot read $ref_matfile\n";
  while(<GENEMATR>) {  
    #source:.. taxon1 taxon2 .. taxonN
    #chr:unsorted    NA      NA      ..      
    #Os4530.POR.1.pan0000001 Os01g0115300,LOC_Os01g02530     -       ..
    chomp;
    my @data = split(/\t/, $_);
    if($data[0] =~ /^source:/) {		
      @ref_taxa = @data[1 .. $#data];

	} elsif($data[0] =~ /^chr:/) { 

    } else {
      push(@ref_cluster_names, shift(@data));
      push(@ref_clusters, \@data);	   
      $n_ref_clusters++;
    }
  }
  close(GENEMATR);

  # 3.2.3) read input pangene_matrix_genes file
  open(GENEMAT,"<",$matfile) ||
    die "# ERROR: cannot read $matfile\n";
  while(<GENEMAT>) {
    #source:.. taxon1 taxon2 .. taxonN
    #chr:unsorted    NA      NA      ..      
    #Os4530.POR.1.pan0000001 Os01g0115300,LOC_Os01g02530     -       ..
    chomp;
    my @data = split(/\t/, $_);
    if($data[0] =~ /^source:/) {
      @taxa = @data[1 .. $#data];

	  # check common taxa
	  foreach $col (0 .. $#taxa) {
        foreach $colr (0 .. $#ref_taxa) {
          next if(defined($input2ref{$colr})); 			
          if($ref_taxa[$colr] eq $taxa[$col]) {
            $input2ref{$colr} = $col; #print "$colr $col $ref_taxa[$colr] $taxa[$col]\n";
            $n_common_taxa++;			
            last;		 	
		  }
        }		  
      }
      if($n_common_taxa == 0) {
        die "# ERROR: no common taxa between $matfile and $ref_matfile, exit\n";
      }

    } elsif($data[0] =~ /^chr:/) {

    } else {
      push(@cluster_names, shift(@data));

      # save matrix row in same order as reference matrix
      my @row;	   
	  foreach $colr (0 .. $#ref_taxa) {
        push(@row,$data[$input2ref{$colr}]);
      }		  
      push(@clusters, \@row);
    }
  }
  close(GENEMAT);

  # 3.2.4) match input to reference clusters (sort & join)
  my ($infh, $tempmatfile) = tempfile(UNLINK => 1);
  my ($infhr, $tempmatfiler) = tempfile(UNLINK => 1);
  my ($infhs, $tempmatfiles) = tempfile(UNLINK => 1);
  my ($infhrs, $tempmatfilers) = tempfile(UNLINK => 1);

  foreach $row (0 .. $#clusters) { 
    print $infh join(':',@{$clusters[$row]})."\t$cluster_names[$row]\n" 
  }

  foreach $rowr (0 .. $#ref_clusters) { 
    print $infhr join(':',@{$ref_clusters[$rowr]})."\t$ref_cluster_names[$rowr]\n" 
  }

  system("$SORTBIN $tempmatfile > $tempmatfiles");
  system("$SORTBIN $tempmatfiler > $tempmatfilers");

  $cmd = "$JOINBIN -o 1.2,2.2 -a 1 $tempmatfiles $tempmatfilers"; # > $tempjoinfile";

  open(SYSJOIN,"$cmd |") ||
    die "# ERROR: cannot run $cmd\n";
  while(<SYSJOIN>) {
    #  gene:OsGoSa_12g0016870 Os4530.POR.1.pan0040985
    #  OsMH63_12G017040 
    if(/^(\S+)\s+(\S+)/){ 
      $input2ref_cluster{ $1 } = $2;
	  $n_identical_clusters++;
    } elsif(/^(\S+)/){
      $input2ref_cluster{ $1 } = 'NA';
    }		
  }
  close(SYSJOIN);  


  # 3.2.5) match input to reference clusters and decide appropriate names, case-by-case
  #my %seen;
  #foreach $row (0 .. $#clusters) {
  #  foreach $rowr (0 .. $#ref_clusters) {
  #    next if(defined($seen{$rowr}));		
  #    #print join("\t",@{$clusters[$row]})."\n";		
  #    if(join("\t",@{$clusters[$row]}) eq join("\t",@{$ref_clusters[$rowr]})) {
  #      #print "$cluster_names[$row] $ref_cluster_names[$rowr]\n";
#		$n_identical_clusters++;
   #     $seen{$rowr}++;		
  #      last;		
  #    }
  #  }
  #}
  
  print "# identical clusters: $n_identical_clusters / $n_ref_clusters\n";
}
