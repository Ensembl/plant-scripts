#!/usr/bin/env perl
use strict;
use warnings;
use Net::FTP;
use HTTP::Tiny;
use JSON;
use Data::Dumper;

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

my $DIVISION  = 'Plants';
my $NCBICLADE = 3700; # NCBI Taxonomy id, Viridiplantae=33090, Arabidopsis=3701
my $REFGENOME = 'arabidopsis_thaliana'; # should be diploid and contained in $NCBICLADE
my $OUTGENOME = ''; # in case you need an outgroup
my $SEQTYPE   = 'cdna'; # cdna, protein

my $VERBOSE   = 1;

# Ensembl Genomes 
my $FTPURL     = 'ftp.ensemblgenomes.org'; 
my $COMPARADIR = '/pub/plants/current/tsv/ensembl-compara/homologies';
my $RESTURL    = 'http://rest.ensembl.org';
my $TREEPOINT  = $RESTURL.'/genetree/member/id/';


# http://ensemblgenomes.org/info/access/mysql
my $DBSERVER  = 'mysql-eg-publicsql.ebi.ac.uk';
my $DBUSER    = 'anonymous';
my $DBPORT    = 4157;

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
   -host    => $DBSERVER,
   -user    => $DBUSER,
   -port    => $DBPORT,
);

my $http = HTTP::Tiny->new();

## 1) check species in clade ##################################################################

my $n_of_species = 0;
my (@supported_species, %polyploid, %supported, %core );
my ($ftp, $compara_file, $stored_compara_file);

# columns of TSV file 
my ($gene_stable_id,$prot_stable_id,$species,$identity,$homology_type,$hom_gene_stable_id,
   $hom_prot_stable_id,$hom_species,$hom_identity,$dn,$ds,$goc_score,$wga_coverage,
	   $high_confidence,$homology_id);

# https://rest.ensembl.org/info/genomes/taxonomy/Homo%20sapiens?content-type=application/json

# get a metadata adaptor
my $e_gdba = Bio::EnsEMBL::MetaData::DBSQL::GenomeInfoAdaptor->build_ensembl_genomes_adaptor();

# find and iterate over all species suupported Ensembl Plants within $NCBICLADE
for my $genome (@{$e_gdba->fetch_all_by_taxonomy_branch($NCBICLADE)}) {
	push(@supported_species, $genome->name());
	$supported{ $genome->name() } = 1;
	print "# ".$genome->name()."\n" if($VERBOSE);
}

printf("# supported species in NCBICLADE %d : %d\n\n", $NCBICLADE, scalar(@supported_species));

# check reference genome is supported
if(!grep(/$REFGENOME/,@supported_species)){
	die "# ERROR: cannot find 'arabidopsis_thaliana' in \$NCBICLADE=$NCBICLADE\n";
}

# add outgroup if required
if($OUTGENOME){
	push(@supported_species,$OUTGENOME);
	$supported{ $OUTGENOME } = 1;
	print "# outgenome: $OUTGENOME\n";
}

$n_of_species = scalar( @supported_species );
print "# total selected species : $n_of_species\n\n";

# check for polyploid species and find reference genome index
my $genome_db_adaptor = $registry->get_adaptor('Plants', 'compara', 'GenomeDB');

#list of supported species, assembly_level 
#https://rest.ensembl.org/info/genomes/division/EnsemblPlants?content-type=application/json

my @gdbs = @{ $genome_db_adaptor->fetch_all_by_mixed_ref_lists(
	-SPECIES_LIST => \@supported_species ) };

foreach my $sp (@gdbs){
	if($sp->is_polyploid()){	
		$polyploid{$sp->name()} = 1;
		printf("%s is polyploid\n",$sp->name()) if($VERBOSE);
	}
}

if($polyploid{$REFGENOME}){
	die "# ERROR: $REFGENOME is polyploid; \$REFGENOME must be diploid\n";
}

printf("# polyploid species : %d\n\n",scalar(keys(%polyploid)));

## 2) get orthologous (plant) genes shared by $REFGENOME and other clade species ################

print "# connecting to $FTPURL ...\n";

if($ftp = Net::FTP->new($FTPURL,Passive=>1,Debug =>0,Timeout=>60)){
	$ftp->login("anonymous",'-anonymous@') || 
		die "# cannot login ". $ftp->message();
	$ftp->cwd($COMPARADIR) || 
		die "# ERROR: cannot change working directory to $COMPARADIR ". $ftp->message();
	$ftp->cwd($REFGENOME) || 
		die "# ERROR: cannot find $REFGENOME in $COMPARADIR ". $ftp->message();

	# find out which file is to be downloaded and 
	# work out its final name with $REFGENOME in it
	$compara_file = '';
	foreach my $file ( $ftp->ls() ){
		if($file =~ m/protein_default.homologies.tsv.gz/){
			$compara_file = $file;
			$stored_compara_file = $compara_file;
			$stored_compara_file =~ s/tsv.gz/$REFGENOME.tsv.gz/;
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

	next if(!$supported{ $hom_species } || $hom_species eq $REFGENOME);

	#https://www.ensembl.org/info/genome/compara/Ortholog_qc_manual.html#wga
	#next if($wga_coverage eq 'NULL' || $wga_coverage < 100);

	#https://www.ensembl.org/info/genome/compara/Ortholog_qc_manual.html#gic
	#next if($goc_score eq 'NULL' || $goc_score < 100);

	if($homology_type eq 'ortholog_one2one' || 
		$polyploid{ $hom_species } && $homology_type eq 'ortholog_one2many'){

		# add $REFGENOME protein 
		if(!$core{ $gene_stable_id }){ 

			push(@{ $core{ $gene_stable_id }{ $REFGENOME } }, $prot_stable_id );

			push(@sorted_ids, $gene_stable_id); # save cluster order
		}

		push(@{ $core{ $gene_stable_id }{ $hom_species } }, $hom_prot_stable_id );
	}
}
close(TSV);

## 3) print summary matrix of single-copy / core genes and compile sequence clusters #################

my $total_core_clusters = 0;
my ($pruned_species,$acc,$seq,$line);

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

	# matrix
	print "$gene_stable_id";
	foreach $hom_species (@supported_species){ 

		printf("\t%s", join(',',@{ $core{ $gene_stable_id }{ $hom_species } }) );

		# store which prots come from each species
		foreach $hom_prot_stable_id (@{ $core{ $gene_stable_id }{ $hom_species } }){
			$valid_prots{ $hom_prot_stable_id } = $hom_species;
		}
	} print "\n";

	# retrieve cluster sequences
	my $response = $http->get( "$TREEPOINT$gene_stable_id?compara=$DIVISION;aligned=1;sequence=$SEQTYPE;$pruned_species", 
		{ headers => { 'Content-type' => 'application/json' } });
			 
	die "Failed!\n" unless $response->{success};

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

		foreach $hom_species (@supported_species){
			foreach $hom_prot_stable_id (@{ $core{ $gene_stable_id }{ $hom_species } }){
				print ">$hom_species $hom_prot_stable_id\n$align{ $hom_species }{ $hom_prot_stable_id }\n";
			}
		}	

	} else {
		die "# ERROR: cannot retrieve sequences for $gene_stable_id\n";
	}

	$total_core_clusters++;

	last if($total_core_clusters == 1);
}

print "# total single-copy core clusters : $total_core_clusters\n";
