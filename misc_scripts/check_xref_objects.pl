#!/usr/bin/env perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use Bio::SeqIO;
use feature qw /say/;
use DBI;

#my $db = 'arabidopsis_thaliana_core_38_91_11';
my $db = 'physcomitrella_patens_core_38_91_11';

{
    my $dsn = "DBI:mysql:database=$db;host=mysql-eg-prod-1.ebi.ac.uk;port=4238";
    my $dbh = DBI->connect($dsn, 'ensrw', 'writ3rp1');

    my @obj_xref_ids = get_ids($dbh);
    
    my $dsn = "DBI:mysql:database=$db;host=mysql-eg-staging-2.ebi.ac.uk;port=4275";
    my $dbh = DBI->connect($dsn, 'ensro', '');

    my $count = 0;
    my $sql;
    for my $id (@obj_xref_ids){
        $sql = "select logic_name from object_xref join analysis using (analysis_id) where object_xref_id = $id";
        my $sth = $dbh->prepare($sql);
        $sth->execute();
        
        my $ref = $sth->fetchrow_hashref();
        my $logic_name = $ref->{logic_name};
        #say $logic_name;
        if ($logic_name ne 'goa_import'){
            say $logic_name;
        }
        $count++;
        if ($count % 1000 == 0){
            warn $count;
        }
    }


}

#======================================== 
sub get_ids {
#======================================== 
    my ($dbh) = @_; 
    my @ids;

    my $sql = qq{
    SELECT distinct(ontology_xref.object_xref_id) FROM ontology_xref LEFT JOIN object_xref ON 
    ontology_xref.object_xref_id = object_xref.object_xref_id WHERE object_xref.object_xref_id IS NULL;
    };

    my $count = 0;

    my $sth = $dbh->prepare($sql);
    $sth->execute();
    while (my $ref = $sth->fetchrow_hashref()) {
        my $obj_xref_id   = $ref->{object_xref_id};
        push (@ids, $obj_xref_id);
        $count++;
        #say $count;
        #if ($count == 100){
        #    last;
        #}
    }
    $sth->finish();
    return (@ids);
}


