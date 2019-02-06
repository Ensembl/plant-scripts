#!/usr/bin/env perl
# gets available EVA studies for a user-provided GCA_accession
# 2018-9 Bruno Contreras Moreira EMBL-EBI

use strict;
use warnings;
use Getopt::Std;
use REST::Client;
use JSON qw(from_json);

my $EVASTUDYLIST = 'https://www.ebi.ac.uk/eva/webservices/rest/v1/meta/studies/all';

my ($GCA_acc,$TSV_file,$species,%wantedGCA,%opts);

getopts('ha:f:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)){
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-a GCA accession                              (example: -a GCA_000001735.2)\n";
  print "-f TSV file with GCA accessions in 1st column (example: -f Viridiplantae.2019-01-22.report.tsv)\n\n";
  exit(0);
}

if($opts{'a'}){ 
	$GCA_acc = $opts{'a'};
	$wantedGCA{$GCA_acc} = 1;
}
elsif($opts{'f'}){
	$TSV_file = $opts{'f'};
	open(TSV,"<",$TSV_file) || die "# ERROR: cannot read $TSV_file\n";
	while(<TSV>) {
		next if(/^#/);
		my @data = split(/\t/,$_);
		$GCA_acc = $data[0];
		$species = $data[1];
		$wantedGCA{$GCA_acc} = $species;
	}
	close(TSV);
}
else {
	die "# ERROR: please use -a OR -f\n\n";
}

my $cl=REST::Client->new(); 
$cl->GET($EVASTUDYLIST,{Accept => 'application/json'}); 
my $json = from_json($cl->responseContent());

foreach my $study (@{ $json->{'response'}->[0]->{'result'} }){

	$GCA_acc = $study->{'assemblyAccession'};

	next if(!$wantedGCA{ $GCA_acc });

	$species = $wantedGCA{$GCA_acc};
	if($species eq '1'={ $species = 'NA' }
	
	print "$study->{'id'}\t$study->{'assemblyAccession'}\t$species\t$study->{'name'}\n";
}

