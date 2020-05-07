#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);
use Benchmark;
use HTTP::Tiny;
use JSON qw(decode_json);
use FindBin '$Bin';
use lib $Bin;
use ComparaUtils qw(
	perform_rest_action $REQUEST_COUNT @DIVISIONS $BASEURL $ANNOTURL
);

# Gets genome annotations of species in a (plant) clade from Ensembl Genomes
#
# Bruno Contreras Moreira 2020

# Ensembl Genomes / NV
my $RESTURL    = 'http://rest.ensembl.org';
my $INFOPOINT  = $RESTURL.'/info/genomes/division/';
my $TAXOPOINT  = $RESTURL.'/info/genomes/taxonomy/';

my $division   = 'Plants';
my $taxonid    = ''; # NCBI Taxonomy id, Brassicaceae=3700, Asterids=71274, Poaceae=4479
my ($outfolder,$out_genome) = ('','');

my ($help,$sp,$line,$id,$show_supported,$request,$response,$baseurl);
my (@ignore_species, %ignore, %division_supported);

GetOptions(	
	"help|?"       => \$help,
	"supported|l"  => \$show_supported,
	"division|d=s" => \$division, 
	"clade|c=s"    => \$taxonid,
	"outgroup|o=s" => \$out_genome,
	"ignore|i=s"   => \@ignore_species,
	"outfile|f=s"  => \$outfolder
) || help_message(); 

sub help_message {
	print "\nusage: $0 [options]\n\n".
      "-c NCBI Taxonomy clade of interest         (required, example: -c Brassicaceae or -c 3700)\n".
      "-f output folder to save HTML files        (required, example: -f Brassicaceae_annot)\n".
      "-l list supported species_names            (optional, example: -l)\n".
      "-d Ensembl division                        (optional, default: -d $division)\n".
      "-o outgroup species_name                   (optional, example: -o brachypodium_distachyon)\n".
      "-i ignore species_name(s)                  (optional, example: -i selaginella_moellendorffii -i ...)\n\n".

    print "Example calls:\n\n".
        " perl $0 -c Liliopsida -o arabidopsis_thaliana -f Liliopsoda.Atha.EG44.annot\n";
    exit(0);
}

if($help){ help_message() }

if($division){
   if(!grep(/^$division$/,@ComparaUtils::DIVISIONS)){
      die "# ERROR: accepted values for division are: ".join(',',@ComparaUtils::DIVISIONS)."\n"
   } else {
      $baseurl = $ComparaUtils::BASEURL;
      my $lcdiv = lc($division);
      $baseurl =~ s/xxx/$lcdiv/;	
   }
}

if($show_supported){ 
	print "# $0 -d $division -l \n\n" 
}
else { 

	if($taxonid eq ''){
		print "# ERROR: need a valid NCBI Taxonomy clade, such as -c Brassicaceae or -c 3700\n\n";
		print "# Check https://www.ncbi.nlm.nih.gov/taxonomy\n";
		exit;
	} else {
		$taxonid =~ s/\s+/%20/g;
	}
	
	if(@ignore_species){
		foreach my $sp (@ignore_species){
			$ignore{ $sp } = 1;
		}
		printf("\n# ignored species : %d\n\n",scalar(keys(%ignore)));
	}

	if(!$outfolder){
		print "# ERROR: need a valid output file, such as -f Brassicaceae.tsv\n\n";
	    exit;
	} 

	print "# $0 -d $division -c $taxonid -o $out_genome -f $outfolder\n\n";
}

my $start_time = new Benchmark();

# new object and params for REST requests
my $http = HTTP::Tiny->new();
my $global_headers = { 'Content-Type' => 'application/json' };
$ComparaUtils::REQUEST_COUNT = 0; 

## 0) check supported species in division ##################################################

$request = $INFOPOINT."Ensembl$division?";

$response = perform_rest_action( $http, $request, $global_headers );
my $infodump = decode_json($response);

foreach $sp (@{ $infodump }) {
	if($sp->{'assembly_accession'}){ 
		$division_supported{ $sp->{'name'} } = 1;
	}	
}

# list supported species and exit
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

my ($n_of_species, $n_of_sequences) = ( 0, 0 );
my (@supported_species, %supported);

$request = $TAXOPOINT."$taxonid?";

$response = perform_rest_action( $http, $request, $global_headers );
$infodump = decode_json($response);

foreach $sp (@{ $infodump }) {
   if($sp->{'name'} && $division_supported{ $sp->{'name'} }){

		next if($ignore{ $sp->{'name'} });

		# add sorted clade species
		$supported{ $sp->{'name'} } = 1;
		push(@supported_species, $sp->{'name'});
   }
}

printf("# supported species in NCBI taxon %s : %d\n\n", $taxonid, scalar(@supported_species));

# add outgroup if required
if($out_genome){
	push(@supported_species,$out_genome);
	$supported{ $out_genome } = 1;
	print "# outgenome: $out_genome\n";
}

$n_of_species = scalar( @supported_species );
print "# total selected species : $n_of_species\n\n";

## 2) get annotation (static content) for selected (plant) species ####################

mkdir($outfolder) if not(-e $outfolder) ;

# iteratively get ans save HTML files
my ($url,$outfile);
$global_headers = { 'Content-Type' => 'text/html' };
foreach $sp ( @supported_species ){

	# build URL
	$url = $baseurl;
    $url .= "/".ucfirst($sp);
    $url .= $ComparaUtils::ANNOTURL;

	# build filename
	$outfile = $outfolder.'/'.$sp.'.html';

	open(OUTFILE,">",$outfile) || die "# ERROR: cannot create $outfile\n";
	$response = perform_rest_action( $http, $url, $global_headers );
	print OUTFILE $response;
	close(OUTFILE);

	print "# written $outfile\n";
}
