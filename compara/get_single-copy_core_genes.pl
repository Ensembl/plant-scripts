#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);
use Net::FTP;
use JSON;
use Data::Dumper;
use Benchmark;

use HTTP::Tiny;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::MetaData::DBSQL::GenomeInfoAdaptor;

# Retrieves all single-copy orthologous genes/proteins shared by (plant) species in clade 
# by querying pre-computed data from Ensembl Genomes Compara. 
# Multiple copies are allowed only in polyploid species.
#
# Based on scripts at
# https://github.com/Ensembl/ensembl-compara/blob/release/97/scripts/examples/
#
# NOTE: use previous release of ensembl-compara, ie relase/97 to query Ensembl Plants 98
#
# Bruno Contreras Moreira 2019

# Ensembl Genomes 
my @divisions  = qw( Plants Bacteria Fungi Vertebrates Protists Metazoa );
my $FTPURL     = 'ftp.ensemblgenomes.org'; 
my $COMPARADIR = '/pub/plants/current/tsv/ensembl-compara/homologies';
my $RESTURL    = 'http://rest.ensembl.org';
my $INFOPOINT  = $RESTURL.'/info/genomes/division/';
my $TREEPOINT  = $RESTURL.'/genetree/member/id/';

# http://ensemblgenomes.org/info/access/mysql
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
   -host    => 'mysql-eg-publicsql.ebi.ac.uk',
   -user    => 'anonymous',
   -port    => 4157,
);

my $verbose    = 0;
my $division   = 'Plants';
my $taxonid    = 3700; # NCBI Taxonomy id, Viridiplantae=33090, Arabidopsis=3701
my $ref_genome = 'arabidopsis_thaliana'; # should be diploid and contained in $taxonid;
my $seqtype    = 'protein'; 
my $outfolder  = '';
my $out_genome = '';

my ($help,$sp,$show_supported,$request,$response,$request_time,$last_request_time);
my (@poly_species, %polyploid, %division_supported);

GetOptions(	
	"help|?"       => \$help,
	"verbose|v"    => \$verbose,
	"supported|l"  => \$show_supported,
	"division|d=s" => \$division, 
	"clade|c=i"    => \$taxonid,
	"reference|r=s"=> \$ref_genome,
	"outgroup|o=s" => \$out_genome,
	"polyploid|p=s"=> \@poly_species,
	"type|t=s"     => \$seqtype,
	"folder|f=s"   => \$outfolder
) || help_message(); 

sub help_message {
	print "\nusage: $0 [options]\n\n".
		"-l list supported species_names       (optional, example: -l)\n".
		"-d Ensembl division                   (optional, default: -d $division)\n".
		"-c NCBI Taxonomy clade of interest    (optional, default: -c $taxonid)\n".
		"-r reference species_name             (optional, default: -r $ref_genome)\n".
		"-o outgroup species_name              (optional, example: -o brachypodium_distachyon)\n".
		"-p polyploid species_name(s)          (optional, example: -p triticum_aestivum -s triticum_turgidum)\n".
		"-f folder to output FASTA files       (optional, example: -f myfolder)\n".
		"-t sequence type [protein|cdna]       (optional, requires -f, default: -t protein)\n".
		"-v verbose                            (optional, example: -v\n\n";
		exit(0);
}

if($help){ help_message() }

if($division && !grep(/$division/,@divisions)){
	die "# ERROR: accepted values for division are: ".join(',',@divisions)."\n"
}

if(@poly_species){
	foreach my $sp (@poly_species){
		$polyploid{ $sp } = 1;
	}
	printf("# polyploid species : %d\n\n",scalar(keys(%polyploid)));
}

if($outfolder){
	if(-e $outfolder){ print "# WARNING : folder $outfolder exists, files might be overwritten\n" }
	else { 
		if(!mkdir($outfolder)){ die "# ERROR: cannot create $outfolder\n" }
	}

	if($seqtype ne 'protein' && $seqtype ne 'cdna'){
		die "# ERROR: accepted values for seqtype are: protein|cdna\n"
	}
}	

if($show_supported){ print "# $0 -l \n\n" }
else {
	print "# $0 -d $division -c $taxonid -r $ref_genome -o $out_genome -f $outfolder -t $seqtype\n\n";
}

## 0) check supported species in division ##################################################

my $start_time = new Benchmark();

my $http = HTTP::Tiny->new();

$request = $INFOPOINT."Ensembl$division?";
$response = $http->get( $request, { headers => { 'Content-type' => 'application/json' } });

unless ($response->{success}) {
	die "# ERROR: failed REST request $request\n";
}

my $infodump = decode_json($response->{content});

foreach $sp (@{ $infodump }) {
	if($sp->{'has_peptide_compara'}){
		$division_supported{ $sp->{'name'} } = 1;
	}	
}

if($show_supported){

	foreach $sp (sort(keys(%division_supported))){
		print "$sp\n";
	}
	exit;
}

# check outgroup is supported
if($out_genome && !$division_supported{ $out_genome }){
	die "# ERROR: genome $out_genome is not supported\n";
}

## 1) check species in clade ##################################################################

my $n_of_species = 0;
my (@supported_species, %supported, %core );
my ($ftp, $compara_file, $stored_compara_file);

# columns of TSV file 
my ($gene_stable_id,$prot_stable_id,$species,$identity,$homology_type,$hom_gene_stable_id,
   $hom_prot_stable_id,$hom_species,$hom_identity,$dn,$ds,$goc_score,$wga_coverage,
	   $high_confidence,$homology_id);

# https://rest.ensembl.org/info/genomes/taxonomy/Homo%20sapiens?content-type=application/json

# get a metadata adaptor
my $e_gdba = Bio::EnsEMBL::MetaData::DBSQL::GenomeInfoAdaptor->build_ensembl_genomes_adaptor();

# find and iterate over all species supported Ensembl Plants within $taxonid
for $sp (@{$e_gdba->fetch_all_by_taxonomy_branch($taxonid)}) {

	if($division_supported{ $sp->name() }){
		push(@supported_species, $sp->name());
		$supported{ $sp->name() } = 1;
		print "# ".$sp->name()."\n" if($verbose);
	}
}

printf("# supported species in NCBI taxon %d : %d\n\n", $taxonid, scalar(@supported_species));

# check reference genome is supported and is not polyploid
if(!grep(/$ref_genome/,@supported_species)){
	die "# ERROR: cannot find $ref_genome in \$taxonid=$taxonid\n";
}
elsif($polyploid{ $ref_genome }){
	   die "# ERROR: $ref_genome is polyploid; reference genome must be diploid\n";
}

# add outgroup if required
if($out_genome){
	push(@supported_species,$out_genome);
	$supported{ $out_genome } = 1;
	print "# outgenome: $out_genome\n";
}

$n_of_species = scalar( @supported_species );
print "# total selected species : $n_of_species\n\n";

## 2) get orthologous (plant) genes shared by $ref_genome and other species ####################

print "# connecting to $FTPURL ...\n";

if($ftp = Net::FTP->new($FTPURL,Passive=>1,Debug =>0,Timeout=>60)){
	$ftp->login("anonymous",'-anonymous@') || 
		die "# cannot login ". $ftp->message();
	$ftp->cwd($COMPARADIR) || 
		die "# ERROR: cannot change working directory to $COMPARADIR ". $ftp->message();
	$ftp->cwd($ref_genome) || 
		die "# ERROR: cannot find $ref_genome in $COMPARADIR ". $ftp->message();

	# find out which file is to be downloaded and 
	# work out its final name with $ref_genome in it
	$compara_file = '';
	foreach my $file ( $ftp->ls() ){
		if($file =~ m/protein_default.homologies.tsv.gz/){
			$compara_file = $file;
			$stored_compara_file = $compara_file;
			$stored_compara_file =~ s/tsv.gz/$ref_genome.tsv.gz/;
			last;
		}
	}

	# download that TSV file
	unless(-s $stored_compara_file){
		$ftp->binary();
		my $downsize = $ftp->size($compara_file);
		$ftp->hash(\*STDOUT,$downsize/20) if($downsize);
		printf("# downloading %s (%1.1fMb) ...\n",$compara_file,$downsize/(1024*1024));
		print "# [        50%       ]\n# ";
		if(!$ftp->get($compara_file)){
			die "# ERROR: failed downloading $compara_file\n";
		}
	
		# rename file to final name
		rename($compara_file, $stored_compara_file);
		print "# using $stored_compara_file\n\n";
	} else {
		print "# re-using $stored_compara_file\n\n";
	}
} else { die "# ERROR: cannot connect to $FTPURL , please try later\n" }

# parse TSV file 
my (@sorted_ids);
open(TSV,"gzip -dc $stored_compara_file |") || 
	die "# ERROR: cannot open $stored_compara_file\n";
while(<TSV>){
	
	($gene_stable_id,$prot_stable_id,$species,$identity,$homology_type,$hom_gene_stable_id,
	$hom_prot_stable_id,$hom_species,$hom_identity,$dn,$ds,$goc_score,$wga_coverage,    
	$high_confidence,$homology_id) = split(/\t/);

	next if(!$supported{ $hom_species } || $hom_species eq $ref_genome);

	#https://www.ensembl.org/info/genome/compara/Ortholog_qc_manual.html#wga
	#next if($wga_coverage eq 'NULL' || $wga_coverage < 100);

	#https://www.ensembl.org/info/genome/compara/Ortholog_qc_manual.html#gic
	#next if($goc_score eq 'NULL' || $goc_score < 100);

	if($homology_type eq 'ortholog_one2one' || 
		$polyploid{ $hom_species } && $homology_type eq 'ortholog_one2many'){

		# add $ref_genome protein 
		if(!$core{ $gene_stable_id }){ 

			push(@{ $core{ $gene_stable_id }{ $ref_genome } }, $prot_stable_id );

			push(@sorted_ids, $gene_stable_id); # save cluster order
		}

		push(@{ $core{ $gene_stable_id }{ $hom_species } }, $hom_prot_stable_id );
	}
}
close(TSV);

## 3) print summary matrix of single-copy / core genes and compile sequence clusters #################

my $total_core_clusters = 0;
my ($pruned_species,$acc,$seq,$line,$filename);

# prepare param to prune species in REST requests
foreach $hom_species (@supported_species){
	$pruned_species .= "prune_species=$hom_species;";
} 

# print header
print "cluster";
foreach $hom_species (@supported_species){ 
	print "\t$hom_species";
} print "\n";

foreach $gene_stable_id (@sorted_ids){

	next if(scalar(keys(%{ $core{ $gene_stable_id } })) < $n_of_species); 

	my (%valid_prots, %align);

	$filename = $gene_stable_id;
	if($outfolder){
		if($seqtype eq 'protein'){ $filename .= '.faa' }
		else{ $filename .= '.fna' }
	}

	# matrix
	print $filename;
	foreach $hom_species (@supported_species){ 

		printf("\t%s", join(',',@{ $core{ $gene_stable_id }{ $hom_species } }) );

		# store which prots come from each species
		foreach $hom_prot_stable_id (@{ $core{ $gene_stable_id }{ $hom_species } }){
			$valid_prots{ $hom_prot_stable_id } = $hom_species;
		}
	} print "\n";


	if($outfolder){

		# retrieve cluster sequences
		$request = "$TREEPOINT$gene_stable_id?compara=$division;aligned=1;sequence=$seqtype;$pruned_species";
		$response = $http->get( $request, { headers => { 'Content-type' => 'application/json' } });

		unless ($response->{success}) {
			die "# ERROR: failed REST request $request\n"; 
		}

		if(length($response->{content})){
			my $treedump = Dumper( decode_json($response->{content}) ); 

			foreach $line (split(/\n/, $treedump ) ){
				if($line =~ m/'sequence' =>/){
					($seq,$acc) = ('','');
				}
				elsif($line =~ m/'seq' => '([^']+)'/){
					$seq = $1;
				}
				elsif($line =~ m/'accession' => '([^']+)'/ && $acc eq ''){
					$acc = $1;

					if($valid_prots{ $acc } ){
						$align{ $valid_prots{$acc} }{ $acc } = $seq;
					}
				}
			}
		}

		# save cluster to file
		if(scalar(keys(%align)) == $n_of_species){

			open(FASTA,">","$outfolder/$filename") || 
				die "# ERROR: cannot create $outfolder/$filename\n";

			foreach $hom_species (@supported_species){
				foreach $hom_prot_stable_id (@{ $core{ $gene_stable_id }{ $hom_species } }){
					print FASTA ">$hom_species $hom_prot_stable_id\n$align{ $hom_species }{ $hom_prot_stable_id }\n";
				}
			}	

			close(FASTA);
		} else { die "# ERROR: cannot retrieve sequences for $gene_stable_id\n" }
	}
	
	$total_core_clusters++; #last if($total_core_clusters == 1); # debug
}

print "\n# total single-copy core clusters : $total_core_clusters\n\n";

my $end_time = new Benchmark();
print "\n# runtime: ".timestr(timediff($end_time,$start_time),'all')."\n";
