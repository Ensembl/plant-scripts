#!/usr/bin/env perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::MetaData::DBSQL::GenomeInfoAdaptor;

# Script to retrieve all single-copy genes shared by all plants in Compara 
# Queries the public Ensembl Genomes server by default
# Bruno Contreras Moreira 2019

my $NCBICLADE = 33090; # NCBI Taxonomy id, Viridiplantae
my $REFGENOME = 'arabidopsis_thaliana'; 

# http://ensemblgenomes.org/info/access/mysql
my $DBSERVER  = 'mysql-eg-publicsql.ebi.ac.uk';
my $DBUSER    = 'anonymous';
my $DBPORT    = 4157;

my $VERBOSE = 0;

my (@supported_species);

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_url("mysql://$DBUSER\@$DBSERVER:$DBPORT/");

## 1) check species in clade ####################################################################

# get a metadata adaptor
my $e_gdba = Bio::EnsEMBL::MetaData::DBSQL::GenomeInfoAdaptor->build_ensembl_genomes_adaptor();

# find and iterate over all genomes from Ensembl Plants
for my $genome (@{$e_gdba->fetch_all_by_taxonomy_branch($NCBICLADE)}) {
	push(@supported_species, $genome->name());
	print $genome->name()."\n" if($VERBOSE);
}

printf("# number species in NCBICLADE %d : %d\n\n", $NCBICLADE, scalar(@supported_species));

# check reference genome is supported
if(!grep(/$REFGENOME/,@supported_species)){
	die "# ERROR: cannot find 'arabidopsis_thaliana' in \$NCBICLADE=$NCBICLADE\n";
}

## 

#my $genetree_adaptor = $registry->get_adaptor('Multi','compara','GeneTree');

#my @trees = $genetree_adaptor->fetch_all();


## 2) get gene names of reference genome ########################################################

my @db_adaptors = @{ $registry->get_all_DBAdaptors() };

foreach my $db_adaptor (@db_adaptors) {
    my $db_connection = $db_adaptor->dbc();

    printf(
        "species/group\t%s/%s\ndatabase\t%s\nhost:port\t%s:%s\n\n",
        $db_adaptor->species(),   $db_adaptor->group(),
        $db_connection->dbname(), $db_connection->host(),
        $db_connection->port()
    );
}

#my $gene_adaptor = $registry->get_adaptor( $REFGENOME, "core", "Gene");

#$all_genes = $gene_adaptor->fetch_all_by_external_name


#my $genemember_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi','compara','GeneMember');


 
    ##fetch a gene by its stable identifier
    #my $gene_adaptor = $registry->get_adaptor("oryza_sativa", "core", "gene");
    #my $gene_adaptor = $registry->get_adaptor("triticum_aestivum_cadenza", "core", "gene");
    #my $gene = $gene_adaptor->fetch_by_stable_id('TRIAE_CS42_5AL_TGACv1_375099_AA1216080.2.path1');

    #my $transcript = $gene->canonical_transcript;
    #my $translation = $transcript->translation;
    #say $translation->seq;


