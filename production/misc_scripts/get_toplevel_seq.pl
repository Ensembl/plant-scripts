#!/usr/bin/env perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;

# prints FASTA sequence of toplevel coord_system

my ($dbname, $species, $servercmd, $server_details);
my ($host,$user,$port) = ('','',''); 

if(!$ARGV[2]){ die "# usage: $0 <dbname> <species> <server cmd>\n" }
else {
	($dbname, $species, $servercmd) = @ARGV;

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

foreach my $slice (@slices){
	printf(">%s %s %d-%d\n",	
		$slice->seq_region_name(),
		$slice->coord_system_name(),
		$slice->start(), $slice->end());
	print $slice->seq();
	print "\n";
}
