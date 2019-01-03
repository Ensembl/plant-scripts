#!/usr/bin/env perl
use warnings;
use strict;

# Adapted from Dan Staines's original by B Contreras Moreira 
# EMBL-EBI 2019

use Bio::EnsEMBL::Utils::CliHelper;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Data::Dumper;

my $cli_helper = Bio::EnsEMBL::Utils::CliHelper->new();

# get the basic options for connecting to a database server
my $optsd = [@{$cli_helper->get_dba_opts()}];
my $opts = $cli_helper->process_args($optsd, \&pod2usage);

my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -USER=>$opts->{user},
    -PASS=>$opts->{pass},
    -PORT=>$opts->{port},
    -HOST=>$opts->{host},
    -DBNAME=>$opts->{dbname}
);


my $sql = q/select xx.xref_id,xx.external_db_id,xx.dbprimary_acc,xx.version from (select  x.dbprimary_acc acc, x.external_db_id db, count(*) from xref x group by x.dbprimary_acc, x.external_db_id having count(*)>1) as dups join xref xx on (xx.dbprimary_acc=dups.acc and xx.external_db_id=dups.db) order by xx.external_db_id,xx.dbprimary_acc,xx.version/;

my $tables = {
    object_xref=>["xref_id"],
    associated_xref=>["xref_id"],
    associated_xref=>["source_xref_id"],
    ontology_xref=>["source_xref_id"],
    dependent_xref=>["dependent_xref_id","master_xref_id"],
    external_synonym=>["xref_id"],
};
print "Processing $opts->{dbname}\n";
my $dups = {};
$dba->dbc()->sql_helper()->execute_no_return(
-SQL=>$sql,
    -CALLBACK=>sub {
        my ($row) = @_;
        my ($xref_id,$db_id,$acc,$version) = @$row;
        my $key = join ('-',$db_id,$acc);
        push @{$dups->{$key}}, $xref_id;
        return;
    }
);

print "Found ".scalar(keys (%$dups))." duplicates\n";

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

print "Removed $n duplicates\n";
