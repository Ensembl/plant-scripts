#!/usr/bin/env perl
use strict;
use warnings;

# Removes genes (or any user-selected logic names/biotypes) from a core db. 
# If -file is passed only genes with stable_id in the list are removed. Adapted from 
# https://www.ebi.ac.uk/seqdb/confluence/display/EnsGen/Core+database+production

use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
 
my ($help, $host, $port, $user, $pass, $dbname);
my (@logic_name, @biotype, $listfile, %stable_ids); 

GetOptions(
  "help|?" => \$help,
  "host=s", \$host,
  "P|port=i", \$port,
  "user=s", \$user,
  "p|pass=s", \$pass,
  "dbname=s", \$dbname,
  "logic_name=s", \@logic_name,
  "biotype=s", \@biotype,
  "file=s", \$listfile
) || help_message(); 

if($help){ help_message() }

sub help_message {
	print "\nusage: $0 [options]\n\n".
	"-host          (required)\n".
	"-port          (required)\n".
	"-user          (required)\n".
	"-pass          (required)\n".
	"-dbname        (required)\n".
	"-logic_name    (required, example: -logic_name maker)\n".
	"-biotype       (optional, example: -biotype protein_coding)\n".
	"-file          (optional, example: -file stable_ids.txt)\n";
	exit(0);
}
 
my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new (
    -host   => $host,
    -port   => $port,
    -user   => $user,
    -pass   => $pass,
    -dbname => $dbname,
);
 
if($listfile) {
	open(LIST,"<",$listfile) || die "# ERROR: cannot open $listfile\n";
	while(<LIST>){
		chomp;
		my $stableid = (split)[0];
		$stable_ids{ $stableid }++;
	}
	close(LIST);

	printf("\n# stable ids read from %s : %d\n\n", $listfile, scalar(keys(%stable_ids)));
}

my $ga = $dba->get_adaptor("Gene");
 
my %biotype = map { $_ => 1 } @biotype;

my $deleted = 0;

foreach my $logic_name (@logic_name) {
	my $genes = $ga->fetch_all_by_logic_name($logic_name);
	foreach my $gene (@$genes) {
    	
		if (scalar(@biotype) == 0 || exists $biotype{$gene->biotype}) {
			if(scalar(keys(%stable_ids)) > 0) {
				next if(!defined($stable_ids{ $gene->stable_id }));
			}
			$ga->remove($gene);
			$deleted++;
		}
	}
}

print "# deleted $deleted genes\n\n";

