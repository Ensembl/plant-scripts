#!/usr/bin/env perl

# This script takes a cluster produced by get_pangenes.pl and retrieves
# the collinearity data evidence supporting it, inferred from whole genome
# alignments, contained in mergedpairs.tsv.gz (sorted with -k1,1 -k4,4nr)

# Copyright [2022]
# EMBL-European Bioinformatics Institute & Estacion Experimental Aula Dei-CSIC

$|=1;

use strict;
use warnings;
use Getopt::Std;
use DB_File;
use Compress::Zlib qw(compress uncompress);
use FindBin '$Bin';
use lib "$Bin/lib";
use pangeneTools qw( parse_sequence_FASTA_file );

$ENV{'LC_ALL'}   = 'POSIX';

my $GZIPBIN = $ENV{'EXE_GZIP'} || 'gzip';

my (%opts,$INP_dir,$INP_clusterfile);
my ($cluster_list_file,$cluster_folder);
my ($gene_id, $hom_gene_id, $homology_type, $species, $hom_species);
my ($overlap, $coords, $hom_coords, $full_id, $hom_full_id);
my ($line, $segment, $hom_segment, $dummy, $TSVdata);
my (%TSVdb, @sorted_ids, @pairs);
my (%seen, %overlap, %cluster_gene_id, %fullid2id, %gene_length);

getopts('hd:i:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0))
{
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-d output directory produced by get_pangenes.pl    (example: -d /path/to/data_pangenes/..._algMmap_)\n";
  print "-i cluster filename as shown in .cluster_list file (example: -i gene:ONIVA01G52180.cdna.fna)\n\n";
  print "Note: reads the compressed merged TSV file in -d\n"; 
  exit(0);
}

if(defined($opts{'d'})){  $INP_dir = $opts{'d'} }
else{ die "# EXIT : need a -d directory\n" }

if(defined($opts{'i'})){  $INP_clusterfile = $opts{'i'} }
else{ die "# EXIT : need parameter -i\n" }

# 1) locate .cluster_list file to check clusterfile is there
opendir(INPDIR,$INP_dir) || 
  die "# ERROR: cannot list $INP_dir , please check -d argument is a valid folder\n";
my @files = grep {/\.cluster_list/} readdir(INPDIR);
closedir(INPDIR);

if(@files) {
  $cluster_list_file = $files[0];
  $cluster_folder = (split(/\.cluster_list/,$cluster_list_file))[0]
} else {
  die "# ERROR: cannot find .cluster_list file in $INP_dir\n";
}

my $clusternameOK = 0;
open(LIST,"<","$INP_dir/$cluster_list_file") ||
  die "# ERROR: cannot read $INP_dir/$cluster_list_file, ".
    "please check -d argument is a valid folder\n";

while(<LIST>) {
  if(/$INP_clusterfile/) {
    $clusternameOK = 1;
  }
}
close(LIST);

if($clusternameOK == 0) {
  die "# ERROR: cannot find $INP_clusterfile in $INP_dir/$cluster_list_file, please correct\n";
}

# 2) parse FASTA headers of input cluster to extract gene names
my ( $ref_geneid, $ref_fasta, $ref_coords, $ref_taxon ) = 
  parse_sequence_FASTA_file( "$INP_dir/$cluster_folder/$INP_clusterfile" , 1);

foreach $gene_id (sort @$ref_geneid) {
  $cluster_gene_id{$gene_id} = 1;
  push(@sorted_ids, $gene_id);
} 

printf("\n# cluster %s genes = %d\n",
  $INP_clusterfile, 
  scalar(keys(%cluster_gene_id)));


# 3) parse compressed merged TSV file and feed Berkeley DB
my $mergedTSV = "$INP_dir/mergedpairs.tsv.gz";
my $TSVdb_file = "$INP_dir/mergedpairs.tsv.bdb";

# only first time
if(!-s $TSVdb_file) {

  print "\n# creating database (only first time)\n";

  tie(%TSVdb, 'DB_File', $TSVdb_file, 
    O_RDWR|O_CREAT, 0666, $DB_BTREE) ||
    die "# ERROR: cannot create file $TSVdb_file: $!\n";

  open(TSV, "$GZIPBIN -dc $mergedTSV |") ||
    die "# ERROR: cannot find file $mergedTSV, please check -d argument\n";

  my ($block, $prev_gene_id) = ('', '');
  while(<TSV>) {
    if(/^([^\t]+)/) {
      $gene_id = $1;

      if($gene_id ne $prev_gene_id) {

        if($block) {
          $TSVdb{$prev_gene_id} = compress($block);
        }
        
        $block = $_;
        $prev_gene_id = $gene_id;

      } else {
        $block .= $_;
      }
    }
  }   

  $TSVdb{$gene_id} = compress($block);

  close(TSV);

  print "# done\n";  

} else {

  print "# re-using database\n";

  tie(%TSVdb, 'DB_File', $TSVdb_file,
    O_RDWR, 0666, $DB_BTREE) ||
    die "# ERROR: cannot read file $TSVdb_file: $!\n";
}

# 4) parse TSV blocks for selected genes  
foreach $gene_id (@sorted_ids) {

  next if(!$TSVdb{$gene_id});

  $TSVdata = uncompress($TSVdb{$gene_id});
 
foreach $line (split(/\n/,$TSVdata)) { 

  #LOC_Os01g01019 LOC_Os01g01019 oryza_sativa_MSU 1065 ortholog_collinear
  #gene:Osir64_01g0000020 gene:Osir64_01g0000020 oryza_sativa_ir64 1065	
  #NULL NULL NULL 100.00 1 Chr1:11217-12435(+);1:19441-20506(+)

  ( $gene_id, $segment, $species, $overlap, $homology_type, 
    $hom_gene_id, $hom_segment, $hom_species, $dummy, 
    $dummy, $dummy, $dummy, $dummy, $dummy, $coords   
  ) = split( /\t/, $line );
 
  if($homology_type eq 'ortholog_collinear') {

    next if(!$cluster_gene_id{$gene_id} || !$cluster_gene_id{$hom_gene_id});
    next if($ref_taxon->{$gene_id} ne $species || $ref_taxon->{$hom_gene_id} ne $hom_species);

    # concat species & gene_id in case there are repeated gene ids
    $full_id = $species.$gene_id;
    $hom_full_id = $hom_species.$hom_gene_id;
    $fullid2id{$full_id} = $gene_id;
    $fullid2id{$hom_full_id} = $hom_gene_id;

    $seen{$full_id}++;
    $seen{$hom_full_id}++;

    # overlap only considered for gene pairs (not segments)
    $overlap{$full_id} += $overlap;
    $overlap{$hom_full_id} += $overlap; 

    # compute gene length
    ($coords, $hom_coords) = split(/;/,$coords); 

    if(!$gene_length{$full_id}) {
      if($coords =~ m/^\S+?:(\d+)-(\d+)\([+-]\)/) {
        $gene_length{$full_id} = 1+$2-$1;
      } else {
        die "# ERROR: cannot parse $coords $_\n";
      }
    }

    if(!$gene_length{$hom_full_id}) {
      if($hom_coords =~ m/^\S+?:(\d+)-(\d+)\([+-]\)/) {
        $gene_length{$hom_full_id} = 1+$2-$1;
      } else { 
        die "# ERROR: cannot parse $hom_coords\n";
      }
    }

    push(@pairs, "$line\n");

  } elsif($homology_type eq 'segment_collinear') {

    next if(!$cluster_gene_id{$gene_id} && !$cluster_gene_id{$hom_gene_id});
    next if($cluster_gene_id{$gene_id} && $ref_taxon->{$gene_id} ne $species);
    next if($cluster_gene_id{$hom_gene_id} && $ref_taxon->{$hom_gene_id} ne $hom_species);

    # segment-gene pairs not considered for stats 

    push(@pairs, "$line\n");
  }

}

}

untie(%TSVdb);

# 5) print raw collinear TSV evidence
if(scalar(keys(%seen)) != scalar(keys(%cluster_gene_id))) {
  die "# ERROR: cannot find collinear evidence for cluster $INP_clusterfile ,".
    "please re-run get_pangenes.pl with same arguments\n";
}
else { 
  print "\n#gene_stable_id\tprotein_stable_id\tspecies\toverlap\thomology_type\t" .
    "homology_gene_stable_id\thomology_protein_stable_id\thomology_species\t".
    "overlap\tdn\tds\tgoc_score\twga_coverage\tis_high_confidence\tcoordinates\n";
  print @pairs;
}

# 6) print summary stats
print "\n#gene\tlength\tpairs\tgene_overlap\tspecies\n";
foreach $full_id (sort {$seen{$b} <=> $seen{$a}} (keys(%seen))){

  $gene_id = $fullid2id{$full_id};

  printf("%s\t%d\t%d\t%1.0f\t%s\n",
    $gene_id,
    $gene_length{$full_id},
    $seen{$full_id},
    $overlap{$full_id},
    $ref_taxon->{$gene_id}
  );
}

