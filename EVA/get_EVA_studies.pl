#!/usr/bin/env perl
# gets available EVA studies for a user-provided GCA_accession
# 2018 Bruno Contreras Moreira EMBL-EBI

use strict;
use warnings;
use REST::Client;
use JSON qw(from_json);

my $EVASTUDYLIST = 'https://www.ebi.ac.uk/eva/webservices/rest/v1/meta/studies/all';

my $cl=REST::Client->new(); 
$cl->GET($EVASTUDYLIST,{Accept => 'application/json'}); 
my $json = from_json($cl->responseContent());

foreach my $study (@{ $json->{'response'}->[0]->{'result'} }){
	print "$study->{'id'}\t$study->{'assemblyAccession'}\t$study->{'name'}\n";
}

