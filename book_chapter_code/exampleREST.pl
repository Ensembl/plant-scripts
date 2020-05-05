# Example in Perl of queries to the Ensembl REST endpoints 
# 
# The full set of endpoints is documented 
# at http://rest.ensembl.org

use strict;
use warnings;
use JSON;
use HTTP::Tiny;
use Data::Dumper;

# 1. Create an HTTP client and a helper function for invoking a
# REST endpoint:

# create an HTTP client
my $http = HTTP::Tiny->new;
my $server = 'http://rest.ensembl.org';

# function for invoking endpoint, see other options at
# https://github.com/Ensembl/ensembl-rest/wiki
sub call_endpoint {
	my ($url, $verbose) = @_;
	
	print "Invoking $url\n" if($verbose);

	my $response = $http-> get($url, {headers =>	
		{'Content-type' => 'application/json'} });
	
	return decode_json($response->{content});
}

# 2. Find homologues of A. thaliana DEAR3 gene:
my $gene = 'DEAR3';
my $species = 'arabidopsis_thaliana';
my $division = 'plants';
my $homoltype = 'ortholog';

# Check https://www.ncbi.nlm.nih.gov/taxonomy for taxonomy ids
my $target_taxon = 4479; #Poaceae

my $url =
	join('/', $server, 'homology/symbol', $species, $gene).
	"?content-type=application/json&compara=$division";

if(defined($target_taxon)){
	$url .= "&target_taxon=$target_taxon";
}

# call url endpoint and get a hash back
my $homology_data = call_endpoint($url);

# parse the homologue list from the response
my @homologies = @{$homology_data->{data}[0]{homologies}};

# filter out homologues based on type
if(defined($homoltype)){
	@homologies = grep {
		$_->{type} =~ m/$homoltype/
	} @homologies;
}

# 3. Print some information about the orthologous protein:
for my $homolog (@homologies) {

	my $target_species = $homolog->{target}{species};
	my $target_id = $homolog->{target}{protein_id};
	print "$target_species $homolog->{type} $target_id\n";

	# 4. For each translation, print information about GO annotation
	# using the xrefs/id endpoint:
	$url = join('/', $server, "xrefs/id/$target_id").
		"?content-type=application/json;external_db=GO;all_levels=1";

	my $go_data = call_endpoint($url);
	for my $go (@{$go_data}) {
		print $go->{display_id}, ' ', $go->{description} || 'NA',
		' Evidence: ', join(', ', @{$go->{linkage_types}}),
		"\n";
	}
}
