#!/usr/bin/env perl
use strict;
use warnings;

# from https://www.ebi.ac.uk/seqdb/confluence/display/EnsGen/Core+database+production#Coredatabaseproduction-Reusinganassembly

use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
 
my ($host, $port, $user, $pass, $dbname, @logic_name, @biotype);
 
GetOptions(
  "host=s", \$host,
  "P|port=i", \$port,
  "user=s", \$user,
  "p|pass=s", \$pass,
  "dbname=s", \$dbname,
  "logic_name=s", \@logic_name,
  "biotype=s", \@biotype,
);
 
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
 
foreach my $logic_name (@logic_name) {
  my $genes = $ga->fetch_all_by_logic_name($logic_name);
  foreach my $gene (@$genes) {
    if (scalar(@biotype) == 0 || exists $biotype{$gene->biotype}) {
      $ga->remove($gene);
    }
  }
}
