#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);
use JSON qw(decode_json);
use Time::HiRes;
use HTTP::Tiny;
use FindBin '$Bin';
use lib $Bin;
use ComparaUtils qw(
	perform_rest_action  
   $REQUEST_COUNT 
);

# Retrieves annotations for an Ensembl gene/protein transfered from its orthologues 
# Bruno Contreras Moreira 2019

# Ensembl Genomes 
my $RESTURL    = 'http://rest.ensembl.org';
my $HOMOLPOINT = $RESTURL.'/homology/id/';
my $XREFPOINT  = $RESTURL.'/xrefs/';

my ($help,$xref,$orth,$orth_id,$orth_sp,$request,$response,$dump);
my ($id,$nondirect,$request_time,$last_request_time) = ('',0,0,0);
my (@orths,%annot);

GetOptions(	
	"help|?"       => \$help,
	"nondirect|n"  => \$nondirect,
	"id|i=s"       => \$id,
) || help_message(); 

sub help_message {
	print "\nusage: $0 [options]\n\n".
      "-i gene/protein id                 (required, example: -i MARPO_0110s0028 or -i HORVU5Hr1G095630 -n)\n".
		"-n allow non-direct annotations    (optional, example: -n\n";
		exit(0);
}

if($help){ help_message() }

if($id eq ''){
	die "# ERROR: need a valid gene/protein id, such as -i MARPO_0110s0028 or -i HORVU5Hr1G095630 -n\n";
}

print "# $0 -n $nondirect -i $id\n\n";

# new object and params for REST requests
my $http = HTTP::Tiny->new();
my $global_headers = { 'Content-Type' => 'application/json' };
$ComparaUtils::REQUEST_COUNT = 0;


# get orthologues & their annotations
$request = $HOMOLPOINT."$id?type=orthologues";
$response = perform_rest_action( $http, $request, $global_headers );
$dump = decode_json($response);
@orths = @{ $dump->{'data'}->[0]->{'homologies'} };
if(scalar(@orths == 0)){
	die "# no orthologues found for $id\n\n";
}

foreach $orth (@orths) {

	$orth_sp = $orth->{'target'}->{'species'};
	$orth_id = $orth->{'target'}->{'id'};
	
	# check annotations of this ortholog
	$request = $XREFPOINT."id/$orth_id?";
	$response = perform_rest_action( $http, $request, $global_headers );
	my $xref_dump = decode_json($response);
	
	foreach $xref (@{ $xref_dump }) {

		next if($nondirect == 0 && $xref->{'info_type'} ne "DIRECT");

		# pick reactome annotations only once
		if($xref->{'dbname'} eq "Plant_Reactome_Pathway" && 
				$xref->{'primary_id'} =~ m/R-\w+-(\d+)/){
			$orth_sp = $xref->{'dbname'};
			$orth_id = $1; 
		} elsif($xref->{'dbname'} eq "Plant_Reactome_Reaction"){
			next;
		}


		if($xref->{'description'} && $xref->{'description'} ne '~'){
			$annot{$orth_sp}{$orth_id}{$xref->{'description'}}++;
		}
		elsif($xref->{'display_id'} && $xref->{'display_id'} ne $orth_id){
         $annot{$orth_sp}{$orth_id}{$xref->{'display_id'}}++;
      }
	}
}

if(scalar(keys(%annot)) == 0){
	printf("$id\tNA\tNA\tNA\n");
}
else{
	foreach $orth_sp (keys(%annot)){
		foreach $orth_id (keys(%{ $annot{$orth_sp} })){
			foreach $xref (keys(%{ $annot{$orth_sp}{$orth_id} })){
				printf("%s\t%s\t%s\t%s\n",$id,$orth_sp,$orth_id,$xref);
			}
		}
	}
}

