#!/usr/bin/env perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $FASTAWIDTH = 60;

# prints FASTA sequence of toplevel coord_system

my ($dbname, $species, $servercmd, $mask, $server_details);
my ($host,$user,$port,$logic_name) = ('','','',''); 
my $usage = "# usage: $0 <dbname> <species> <server cmd> <mask: none, hard, soft> [repeats logic_name, default=RepeatMasker]";

if(!$ARGV[3]){ die "$usage\n" }
else {
	($dbname, $species, $servercmd, $mask, $logic_name) = @ARGV;

	# check mask
	if($mask ne 'none' && 
			$mask ne 'hard' &&
			$mask ne 'soft'){
		die "$usage\n"
	}

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
	printf(">%s dna:%s %s %s\n",	
		$slice->seq_region_name(),
		$slice->coord_system_name(),
		$slice->name(),
		$mask );

	my $seq;
	if($mask eq 'hard'){
		if($logic_name ne ''){
			$seq = $slice->get_repeatmasked_seq( [$logic_name] )->seq();
		} else {
			$seq = $slice->get_repeatmasked_seq()->seq();
		}
	} elsif($mask eq 'soft'){
		if($logic_name ne ''){
			$seq = $slice->get_repeatmasked_seq( [$logic_name], 1 )->seq();
		} else {
			$seq = $slice->get_repeatmasked_seq( undef, 1 )->seq();	
		}
	}
	else {
		$seq = $slice->seq();
	}
	
	while($seq =~ m/(.{1,$FASTAWIDTH})/g){
		print "$1\n"
	}
}
