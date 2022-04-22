#!/usr/bin/env perl

# This script takes a cluster produced by get_pangenes.pl and retrieves
# the collinearity data evidence supporting it, inferred from whole genome
# alignments

# Copyright [2022]
# EMBL-European Bioinformatics Institute & Estacion Experimental Aula Dei-CSIC

$|=1;

use strict;
use warnings;
use Getopt::Std;
use FindBin '$Bin';
use lib "$Bin/lib";
use pangeneTools qw(parse_sequence_FASTA_file);

my $GZIPBIN = $ENV{'EXE_GZIP'} || 'gzip';

my (%opts,$INP_dir,$INP_clusterfile);
my ($cluster_list_file,$cluster_folder);
my ($gene_id, $hom_gene_id, $homology_type, $species, $hom_species);
my ($overlap, $coords, $hom_coords, $full_id, $hom_full_id);
my (%seen, %overlap, @pairs, %cluster_gene_id, %fullid2id, %gene_length);

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
  die "# ERROR: cannot read $INP_dir/$cluster_list_file, please check -d argument is a valid folder\n";

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

foreach $gene_id (@$ref_geneid) {
  $cluster_gene_id{$gene_id} = 1;
}

printf("# cluster genes = %d\n", 
  scalar(keys(%cluster_gene_id)));

# 3) parse compressed merged TSV file
my $mergedTSV = "$INP_dir/mergedpairs.tsv.gz";

open(TSV, "$GZIPBIN -dc $mergedTSV |") ||
  die "# ERROR: cannot find file $mergedTSV, please check -d argument\n";

while(<TSV>) {
  #gene_stable_id\tprotein_stable_id\tspecies\toverlap\thomology_type\thomology_gene_stable_id\t
  # homology_protein_stable_id\thomology_species\toverlap\tdn\tds\tgoc_score\twga_coverage\t
  # is_high_confidence\tcoordinates

  my @data = split(/\t/,$_);
  next if(scalar(@data) < 15);

  ($gene_id, $species, $overlap, $homology_type, $hom_gene_id, $hom_species, $coords) = 
    ($data[0],$data[2],$data[3],$data[4],$data[5],$data[7],$data[14]);

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
        die "# ERROR: cannot parse $coords\n";
      }
    }

    if(!$gene_length{$hom_full_id}) {
      if($hom_coords =~ m/^\S+?:(\d+)-(\d+)\([+-]\)/) {
        $gene_length{$hom_full_id} = 1+$2-$1;
      } else { 
        die "# ERROR: cannot parse $hom_coords\n";
      }
    }

    push(@pairs, $_);

  } elsif($homology_type eq 'segment_collinear') {

    next if(!$cluster_gene_id{$gene_id} && !$cluster_gene_id{$hom_gene_id});
    next if($cluster_gene_id{$gene_id} && $ref_taxon->{$gene_id} ne $species);
    next if($cluster_gene_id{$hom_gene_id} && $ref_taxon->{$hom_gene_id} ne $hom_species);

    # segment-gene pairs not considered for stats 
    #if($data[1] ne 'segment'){ $seen{$gene_id}++ }
    #elsif($data[6] ne 'segment'){ $seen{$hom_gene_id}++ }

    push(@pairs, $_)
  }

}
close(TSV);

if(scalar(keys(%seen)) != scalar(keys(%cluster_gene_id))) {
  die "# ERROR: cannot find collinear evidence for cluster $INP_clusterfile ,".
    "please re-run get_pangenes.pl with same arguments\n";
}
else { 
  print "\n#gene_stable_id\tprotein_stable_id\tspecies\toverlap\thomology_type\thomology_gene_stable_id\t" .
    "homology_protein_stable_id\thomology_species\toverlap\tdn\tds\tgoc_score\twga_coverage\t" .
    "is_high_confidence\tcoordinates\n";
  print @pairs;
}

print "\n#species\tgene\tlength\tpairs\tgene_overlap\tspecies\n";
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

