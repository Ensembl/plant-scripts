#!/usr/bin/env perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;

# compute percentage of masked toplevel sequence

my ($dbname, $species, $servercmd, $analysis, $server_details);
my ($host,$user,$port) = ('','',''); 
my $usage = "# usage: $0 <dbname> <species> <server cmd> [analysis name]";

if(!$ARGV[2]){ die "$usage\n" }
else {
	($dbname, $species, $servercmd, $analysis) = @ARGV;

	# get db server connection details
	chomp( $server_details = `$servercmd -details script` );
	if($server_details =~ m/--host (\S+) --port (\d+) --user (\S+)/){
		($host,$port,$user) = ($1,$2,$3);
	} 
}

# open db connection, assume it's a core schema
my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
  -host    => $host,
  -user    => $user,
  -port    => $port,
  -species => $species,
  -dbname  => $dbname
);


my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'Core', 'Slice');
my @slices = @{ $slice_adaptor->fetch_all('toplevel') };

my ($total,$masked, $softseq) = (0,0);
foreach my $slice (@slices){

	if(defined($analysis)){
		$softseq = $slice->get_repeatmasked_seq( ["$analysis"] , 1 )->seq();
	} else {
		$softseq = $slice->get_repeatmasked_seq( undef, 1 )->seq();
		$analysis = 'all';
	}

	$total += length($softseq);
	$masked += $softseq =~ tr/a-z//;
}

printf("# %s %d/%d %2.2f\n", 
	$analysis,
	$masked,
	$total,
	(100*$masked)/$total);
