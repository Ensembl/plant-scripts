#!/usr/bin/env perl
use strict;
use warnings;
use Net::FTP;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::MetaData::DBSQL::GenomeInfoAdaptor;

# Retrieves all single-copy orthologous genes shared by (plant) species in clade 
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
my $NCBICLADE = 3701; # NCBI Taxonomy id, Viridiplantae=33090
my $REFGENOME = 'arabidopsis_thaliana'; # should be diploid and contained in $NCBICLADE
my $OUTGENOME = ''; # in case you need an outgroup

my $VERBOSE   = 1;

# Ensembl Genomes URLs
my $FTPURL     = 'ftp.ensemblgenomes.org'; 
my $COMPARADIR = '/pub/plants/current/tsv/ensembl-compara/homologies';
my $RESTURL    = 'rest.ensembl.org';

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

## 1) check species in clade ##################################################################

my (@supported_species, %polyploid, %supported, %core_genes );
my ($ftp, $compara_file, $stored_compara_file );

my ($sp, $tree, $leaf, $treeOK, $tree_stable_id, $align, $align_string );

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

# check for polyploid species and find reference genome index
my $genome_db_adaptor = $registry->get_adaptor('Plants', 'compara', 'GenomeDB');

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

printf("# polyploid species in clade : %d\n\n",scalar(keys(%polyploid)));

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
my ($gene_stable_id,$prot_stable_id,$species,$identity,$homology_type,$hom_gene_stable_id,
	$hom_protein_stable_id,$hom_species,$hom_identity,$dn,$ds,$goc_score,$wga_coverage,
	$high_confidence,$homology_id);

open(TSV,"gzip -dc $stored_compara_file |") || 
	die "# ERROR: cannot open $stored_compara_file\n";
while(<TSV>){
	
	($gene_stable_id,$prot_stable_id,$species,$identity,$homology_type,$hom_gene_stable_id,
	$hom_protein_stable_id,$hom_species,$hom_identity,$dn,$ds,$goc_score,$wga_coverage,    
	$high_confidence,$homology_id) = split(/\t/);

	next if(!$supported{ $hom_species } || $hom_species eq $REFGENOME);

	if($homology_type eq 'ortholog_one2one' || 
		$polyploid{ $hom_species } && $homology_type eq 'ortholog_one2many'){
	
		print;
	}
}
close(TSV);

