#!/usr/bin/env perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::MetaData::DBSQL::GenomeInfoAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

# Script to retrieve all single-copy genes shared by all (plant) species in clade 
# Queries the Compara db at the public Ensembl Genomes server by default
#
# Based on scripts at
# https://github.com/Ensembl/ensembl-compara/blob/release/97/scripts/examples/
#
# NOTE: requires the previous release of ensembl-compara; use relase/97 to query Ensembl Plants 98
#
# Bruno Contreras Moreira 2019

my $DIVISION  = 'Plants';
my $NCBICLADE = 3701; # NCBI Taxonomy id, Viridiplantae=33090
my $REFGENOME = 'arabidopsis_thaliana'; # should be included in $NCBICLADE
my $OUTGENOME = ''; # in case you need an outgroup
my $ONE2ONE   = 1; # set to 1 to take only single-copy genes in diploid species

my $VERBOSE   = 1;

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

## 1) check species in clade ####################################################################

my (@supported_species, %polyploid);
my ($sp, $tree, $tree_newick, $leaf, $treeOK);

# get a metadata adaptor
my $e_gdba = Bio::EnsEMBL::MetaData::DBSQL::GenomeInfoAdaptor->build_ensembl_genomes_adaptor();

# find and iterate over all genomes from Ensembl Plants
for my $genome (@{$e_gdba->fetch_all_by_taxonomy_branch($NCBICLADE)}) {
	push(@supported_species, $genome->name());
	print $genome->name()."\n" if($VERBOSE);
}

printf("# supported species in NCBICLADE %d : %d\n\n", $NCBICLADE, scalar(@supported_species));

# check reference genome is supported
if(!grep(/$REFGENOME/,@supported_species)){
	die "# ERROR: cannot find 'arabidopsis_thaliana' in \$NCBICLADE=$NCBICLADE\n";
}

# add outgroup if required
if($OUTGENOME){
	push(@supported_species,$OUTGENOME);
	print "# outgenome: $OUTGENOME\n";
}

# check for polyploid species
my $genome_db_adaptor = $registry->get_adaptor('Plants', 'compara', 'GenomeDB');

my @gdbs = @{ $genome_db_adaptor->fetch_all_by_mixed_ref_lists(-SPECIES_LIST => \@supported_species) };
foreach my $sp (@gdbs){
	printf("%s polyploid=%d\n",$sp->name(),$sp->is_polyploid()) if($VERBOSE);
	if($sp->is_polyploid()){	
		$polyploid{$sp->name()} = 1;
	}
}

printf("# polyploid species in clade : %d\n\n",scalar(keys(%polyploid)));


## 2) select (plant) trees including all clade species ############################################

my $gene_tree_adaptor = $registry->get_adaptor('Plants', 'compara', 'GeneTree');

my $all_protein_trees = $gene_tree_adaptor->fetch_all(
    -CLUSTERSET_ID => 'default',
    -MEMBER_TYPE => 'protein',
    -TREE_TYPE => 'tree',
);

foreach $tree (@$all_protein_trees){

	my (%included);
	$treeOK = 1;

	$tree_newick = $tree->newick_format("species");
  	
	# find clade species in this tree
	foreach $sp (@supported_species){
		if($tree_newick =~ /$sp/ig){
			$included{$sp}++;
		}
	}

   # check this tree
	foreach $sp (@supported_species){

		if(!$included{$sp} || 
			($ONE2ONE==1 && !$polyploid{$sp} && $included{$sp} > 1)){
			$treeOK = 0;
		}

		printf("%s %s seqs=%d\n",
			$tree->root_id(),$sp,$included{$sp} || 0) if($VERBOSE)
	}
   
   next if($treeOK == 0);

		

	exit;
}




