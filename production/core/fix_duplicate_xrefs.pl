#!/usr/bin/env perl
use warnings;
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Data::Dumper;
use Getopt::Std;

# Takes a db server details plus the name of core db and eliminates duplicated xrefs.
#
# Example call: fix_duplicate_xrefs.pl -p eg-p3-w -d solanum_lycopersicum_core_42_95_3
#
# Adapted from Dan Staines's original by Bruno Contreras Moreira EMBL-EBI 2019

my $VERBOSE = 0;

my ($db_name,$prod_server,$server_args,$dba,%opts);
getopts('hp:d:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)){
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-d species database name     (required, example: -d solanum_lycopersicum_core_42_95_3)\n";
  print "-p production db server rw   (required, example: -p eg-p3-w)\n\n";
  exit(0);
}

if($opts{'d'}){ $db_name = $opts{'d'} }
else{ die "# EXIT : need a valid -d species db name, such as -f solanum_lycopersicum_core_42_95_3\n" }

if($opts{'p'}){
        $prod_server = $opts{'p'};
        chomp( $server_args = `$prod_server details` );
	if($server_args =~ m/--host=(\S+) --port=(\S+) --user=(\S+) --pass=(\S+)/){
		$dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
			-HOST=>$1,-PORT=>$2,-USER=>$3,-PASS=>$4,-DBNAME=>$db_name
		);
		print "# host=$1\n# port=$2\n# user=$3\n# pass=$4\n" if($VERBOSE);	
	}
	else{ die "# EXIT : cannot parse $prod_server details\n" }
}
else{ die "# EXIT : need a valid -p production db server, such as -p eg-p3-w\n" }


my $sql = q/select xx.xref_id,xx.external_db_id,xx.dbprimary_acc,xx.version from (select x.dbprimary_acc acc, x.external_db_id db, count(*) from xref x group by x.dbprimary_acc, x.external_db_id having count(*)>1) as dups join xref xx on (xx.dbprimary_acc=dups.acc and xx.external_db_id=dups.db) order by xx.external_db_id,xx.dbprimary_acc,xx.version/;

my $tables = {
    object_xref=>["xref_id"],
    associated_xref=>["xref_id"],
    associated_xref=>["source_xref_id"],
    ontology_xref=>["source_xref_id"],
    dependent_xref=>["dependent_xref_id","master_xref_id"],
    external_synonym=>["xref_id"],
};

print "# Checking xref duplicates in $db_name\n";
my $dups = {};
$dba->dbc()->sql_helper()->execute_no_return(
-SQL=>$sql,
    -CALLBACK=>sub {
        my ($row) = @_;
        my ($xref_id,$db_id,$acc,$version) = @$row;
        my $key = join ('-',$db_id,$acc);
        push @{$dups->{$key}}, $xref_id; 
	print "$key $xref_id\n" if($VERBOSE);
        return;
    }
);

print "# Found ".scalar(keys (%$dups))." duplicates\n";

my $n= 0;
while(my ($key,$xrefs) = each %$dups) {
    my $keep;
    for my $xref (@$xrefs) {
        if(!defined $keep) {
            $keep = $xref;
        } else {
            $n++;
            # replace
            while(my ($table,$cols) = each %$tables) {
                for my $col (@$cols) {
                    eval {
                        $dba->dbc()->sql_helper()->execute_update(
                            -SQL=>"update $table set $col=$keep where $col=$xref"
                            );
                    };
                    if($@) {
                        $dba->dbc()->sql_helper()->execute_update(
                            -SQL=>"delete from $table where $col=$xref"
                            );                        
                    }
                }
            }                  
            # delete
            $dba->dbc()->sql_helper()->execute_update(
                -SQL=>"delete from xref where xref_id=$xref"
                );
        }
    } 
}

print "# Removed $n duplicates\n";
