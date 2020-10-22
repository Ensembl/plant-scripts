#!/usr/bin/env perl

# Examples of queries to Ensembl Plants using the native Perl API
# Check the tutorials for more examples:
# https://m.ensembl.org/info/docs/api/core/core_tutorial.html
# https://www.ensembl.org/info/docs/api/compara/compara_tutorial.html
#
# Copyright [2017-2020] EMBL-European Bioinformatics Institute


# Install the Ensembl Perl API and updated env as explained in
# http://www.ensembl.org/info/docs/api/api_installation.html
# http://www.ensembl.org/info/docs/api/api_git.html
# https://m.ensembl.org/info/docs/api/debug_installation_guide.html
# Note you might need some dependencies, such as libmysqlclient-dev
# or Perl modules DBI & DBD::mysql


## A1) Load the Registry object with details of genomes available

use warnings;
use strict;
use List::Util qw(all);
use Bio::EnsEMBL::Registry;

Bio::EnsEMBL::Registry->load_registry_from_db(
	-USER => 'anonymous',
	-HOST => 'mysql-eg-publicsql.ebi.ac.uk',
	-PORT => '4157',
	#-VERBOSE => 1 # uncomment to see dbs loaded
);

## A2) Check which analyses are available for a species

# Note: logic_names are printed for each analysis

my $division = 'plants';
my $species = 'arabidopsis_thaliana';

my $analysis_adaptor = Bio::EnsEMBL::Registry->
	get_adaptor( $species, "core", "analysis" );

foreach my $analysis (sort 
	{$a->logic_name() cmp $b->logic_name()} 
		@{ $analysis_adaptor->fetch_all() }){
		print $analysis->logic_name(), "\n";
} 

## A3) Get soft masked sequences from Arabidopsis thaliana

my $slice_adaptor = Bio::EnsEMBL::Registry->
    get_adaptor($species, 'core', 'Slice');

my ($total,$masked, $softseq) = (0,0);
foreach my $slice (@{ $slice_adaptor->fetch_all('toplevel') }){

	# for brevity consider only the plastome
	next if($slice->seq_region_name() ne 'Pt'); 

    # note Ensembl 1-based inclusive coordinates
	printf(">%s %s %d-%d\n",
		$slice->seq_region_name(),
		$slice->coord_system_name(),
		$slice->start(),
		$slice->end());
	
	# By default repeatmask* analyses, see recipe A2 to list others
	# Repeat analyses include 'repeatmask_redat', 
	# 'repeatmask_nrplants' or 'repeatdetector_curated'
	# $slice->get_repeatmasked_seq( ['repeatmask_redat'], 1 )
	# only print a 50b segment for brevity
	print substr($slice->get_repeatmasked_seq( undef, 1 )->seq(),80,50), "\n";
}

## A4) Get BED file with repeats in chr4

my $chrname = 'chr4';

my $slice = $slice_adaptor->
	fetch_by_region( 'toplevel', $chrname );

	my @repeats = @{ $slice->get_all_RepeatFeatures() };
	my $total_repeats = 0;

foreach my $repeat (@repeats) {

	# for brevity
	last if($total_repeats++ > 10);

	printf("%s\t%d\t%d\t%s\t%s\t%s\n",
		$chrname,
		$repeat->start()-1,
		$repeat->end(),
		$repeat->analysis()->logic_name(),
		$repeat->repeat_consensus()->repeat_class(),
		$repeat->repeat_consensus()->repeat_type() );
}

## A5) Find the DEAR3 gene

# gene of interest and species
my $gene_name = 'DEAR3';

# get a gene adaptor to work with genes from
# the species
my $gene_adaptor = Bio::EnsEMBL::Registry->
	get_adaptor($species, 'core', 'gene');

# find the gene with the specified name using
# the adaptor
my ($gene_obj) = @{$gene_adaptor->
   fetch_all_by_external_name($gene_name)};

## A6) Get the transcript used in Compara analyses

# The canonical transcript is used in the gene tree analysis,
# which usually is the longest translation with no stop codons

printf(">DEAR3 %s\n%s\n",
	$gene_obj->canonical_transcript()->stable_id(),
	$gene_obj->canonical_transcript()->spliced_seq() );

printf(">DEAR3 %s CDS\n%s\n",
    $gene_obj->canonical_transcript()->stable_id(),
	$gene_obj->canonical_transcript()->translateable_seq() );

printf(">DEAR3 %s\n%s\n\n",
    $gene_obj->canonical_transcript()->translation->stable_id(),
    $gene_obj->canonical_transcript()->translate->seq() );

## A7) Find all orthologues of a gene

# get an adaptor to work with genes from compara
my $gene_member_adaptor = Bio::EnsEMBL::Registry->
	get_adaptor($division, 'compara', 'GeneMember');

# find the corresponding gene in compara
my $gene_member = $gene_member_adaptor->
	fetch_by_stable_id($gene_obj->stable_id());

# get an adaptor to work with homologues in compara
my $homology_adaptor = Bio::EnsEMBL::Registry->
	get_adaptor($division, 'compara', 'Homology');

# find all homologues of the gene
my @homologies = @{$homology_adaptor->
	fetch_all_by_Member($gene_member)};

# filter out homologues based on type
@homologies = grep {
	$_->description =~ m/ortholog/
} @homologies;

foreach my $homology (@homologies) {

	# get the protein from the target
	my $target = $homology->get_all_Members->[1];	
	
	printf("%s\t%s\t%s\t%s\n",
		$gene_obj->stable_id(), 
		$species, 
		$target->stable_id(),
		$target->genome_db->name() );
}

## A8) Get markers mapped on chr1D of bread wheat

# Note: only a few plants have markers
# As of release EG47/100:
# triticum_aestivum, oryza_indica, brassica_rapa
#
# Coordinates are returned in BED format

$species = 'triticum_aestivum';
$chrname = '1D';

$slice_adaptor = Bio::EnsEMBL::Registry->
	get_adaptor( $species, 'Core', 'Slice' );

$slice = $slice_adaptor->
	fetch_by_region( 'chromosome', $chrname );

my $total_markers = 0;
foreach my $mf (@{ $slice->get_all_MarkerFeatures() }) {

	last if($total_markers++ > 10); #for brevity

	my $marker = $mf->marker(); 

	printf("%s\t%d\t%d\t%s\t%s\t%s\t%d\n",
		$mf->seq_region_name(),
		$mf->start()-1,      
		$mf->end(), 
		$mf->display_id(),
		$marker->left_primer(),
		$marker->right_primer(),
		$marker->max_primer_dist() );
}


## A9) Find all syntelogues among rices

# Note: GOC=Gene Order Conservation score
# Read more at 
# https://www.ensembl.org/info/genome/compara/Ortholog_qc_manual.html

# get an adaptor to work with comparative sets from compara
my $mlss_adaptor = Bio::EnsEMBL::Registry->
	get_adaptor($division, 'compara', 'MethodLinkSpeciesSet');

# find the mlss that describes orthologies between these
# two rice species
my $mlss = $mlss_adaptor->fetch_by_method_link_type_registry_aliases(
    'ENSEMBL_ORTHOLOGUES', ['oryza_sativa', 'oryza_indica']);

# find all homologues between these two rice species
@homologies = @{$homology_adaptor->
    fetch_all_by_MethodLinkSpeciesSet($mlss)};

# filter out homologues based on local gene order conservation
@homologies = grep {
	$_->goc_score && $_->goc_score >= 75
} @homologies;

my $count = 1;
foreach my $homology (@homologies) {

	# get one orthologue
	my $prot = $homology->get_all_Members->[1];

    # find all orthologues in rice
    my @rice_homologies = @{$homology_adaptor->
        fetch_all_by_Member($prot,
            -METHOD_LINK_TYPE => 'ENSEMBL_ORTHOLOGUES',
            -TARGET_TAXON => 'Oryza')};

    if (all {$_->goc_score && $_->goc_score >= 75}
                @rice_homologies) {
        foreach my $rh (@rice_homologies) {
            printf("%s\t%s\t%s\t%s\t%d\n",
                $rh->get_all_Members->[0]->genome_db->name,
                $rh->get_all_Members->[0]->stable_id,
                $rh->get_all_Members->[1]->genome_db->name,
                $rh->get_all_Members->[1]->stable_id,
                $rh->goc_score);
        }
        print "\n";

        # Only print the first 10 groups
        last if $count++ == 10;
    }
}

## A10) Print all translations for otherfeatures genes

# Note: otherfeatures dbs are Ensembl databases that
# usually contain additional annotation tracks

$count = 0;
$species = 'triticum_aestivum';

$gene_adaptor = Bio::EnsEMBL::Registry->
	get_adaptor($species, "otherfeatures", "gene");
my $genes = $gene_adaptor->fetch_all_by_biotype('protein_coding');

for my $gene (@$genes){
    my $transcripts = $gene->get_all_Transcripts;
    for my $t (@$transcripts){
        if ($t->biotype ne 'protein_coding'){
            next;
        }
        $count++;
        print ">",$gene->stable_id,"\n";

        my $translation = $t->translation;
        my $sequence = $translation->seq;
        print $translation->seq, "\n";
    }

	# Print only first 10, comment for real use
    last if ($count == 10);
}

