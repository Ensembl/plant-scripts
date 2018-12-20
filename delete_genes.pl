#!/usr/bin/env perl
use strict;
use warnings;

# Removes genes (or any user-selected logic names/biotypes) from a core db; adapted from 
# https://www.ebi.ac.uk/seqdb/confluence/display/EnsGen/Core+database+production

use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
 
my ($help, $host, $port, $user, $pass, $dbname, @logic_name, @biotype);
 
GetOptions(
  "help|?" => \$help,
  "host=s", \$host,
  "P|port=i", \$port,
  "user=s", \$user,
  "p|pass=s", \$pass,
  "dbname=s", \$dbname,
  "logic_name=s", \@logic_name,
  "biotype=s", \@biotype,
) || help_message(); 

if($help){ help_message() }

sub help_message {
	print "\nusage: $0 [options]\n\n".
	"-host                 (required)\n".
	"-port                 (required)\n".
	"-user                 (required)\n".
	"-pass                 (required)\n".
	"-dbname               (required)\n".
	"-biotype              (optional, example: -biotype protein_coding)\n".
	"-logic_names          (optional, example: -logic_names maker)\n";
	exit(0);
}
 
my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new
(
    -host   => $host,
    -port   => $port,
    -user   => $user,
    -pass   => $pass,
    -dbname => $dbname,
);
 
my $ga = $dba->get_adaptor("Gene");
 
my %biotype = map { $_ => 1 } @biotype;

my $deleted = 0;

foreach my $logic_name (@logic_name) {
  my $genes = $ga->fetch_all_by_logic_name($logic_name);
  foreach my $gene (@$genes) {
    if (scalar(@biotype) == 0 || exists $biotype{$gene->biotype}) {
      $ga->remove($gene);
      $deleted++;
    }
  }
}

print "# deleted $deleted genes\n\n";

