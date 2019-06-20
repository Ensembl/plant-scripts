#!/usr/bin/env perl

# adapted from eg-pipelines/scripts/xrefs_pipeline/load_sgd_gene_names.pl
# Bruno Contreras 2019

use strict;
use warnings;
use Readonly;


use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Utils::CliHelper;

use Carp;
use Data::Dumper;
use Log::Log4perl qw(:easy);

my $logger = get_logger();

my $cli_helper = Bio::EnsEMBL::Utils::CliHelper->new();

my $optsd = [ @{ $cli_helper->get_dba_opts() } ];
push( @{$optsd}, "file:s" );
push( @{$optsd}, "verbose" );

my $opts = $cli_helper->process_args( $optsd, \&pod2usage );
if ( $opts->{verbose} ) {
  Log::Log4perl->easy_init($DEBUG);
}
else {
  Log::Log4perl->easy_init($INFO);
}

$logger->info( "Loading " . $opts->{dbname} );
my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                     -USER   => $opts->{user},
									   -PASS   => $opts->{pass},
									   -HOST   => $opts->{host},
									   -PORT   => $opts->{port},
									   -DBNAME => $opts->{dbname} );

my $file = $opts->{file};
open my ($fh), '<', $file or die "Could not open $file: $!";

my $gene_adaptor = $dba->get_GeneAdaptor;
my $dea = $dba->get_DBEntryAdaptor;
while ( my $line = <$fh> ) {
    chomp($line);
    my @row = split /\t/, $line;
        
    if ( length $row[3] == 0 ) {
    	next;
    }
    
    my $gene = $gene_adaptor->fetch_by_stable_id($row[3]);
    if ( !$gene ) {
        next;
    }
    
    my $new_name = 0;
    my $new_description = 0;
    
    my $old_display_xref = $gene->display_xref();
    my $new_display_xref = '';
    if ( $old_display_xref ) {
        if ( $old_display_xref->display_id() ne $row[4] and length $row[4] > 0 ) {
            $logger->info( "Modifying " . $row[3] . " - " . $old_display_xref->display_id() . " to " . $row[4]);
            my @synonyms = $old_display_xref->get_all_synonyms();
            
            $new_display_xref = Bio::EnsEMBL::DBEntry -> new (
                -PRIMARY_ID  => $row[0],
                -DBNAME      => 'SGD_GENE',
                -DISPLAY_ID  => $row[4],
                -DESCRIPTION => $row[15],
                -INFO_TYPE   => 'NONE',
            );
            
            foreach my $synonym (@{$synonyms[0]}) {
                #$logger->info( "  Old Synonyms");
                #$logger->info(Dumper $synonym);
                $new_display_xref->add_synonym($synonym);
            }
            
            my $dbRef = $dea->_check_external_db($new_display_xref,0);
            my $xref_id = $dea->_store_or_fetch_xref($new_display_xref,$dbRef);
            $new_display_xref->dbID($xref_id);
            $gene->display_xref($new_display_xref);
            
            $new_name = 1;
        }
    } else {
        if ( length $row[4] > 0 ) {
            $logger->info( "New name " . $row[3] . " - " . $row[4]);
            $new_display_xref = Bio::EnsEMBL::DBEntry -> new (
                -PRIMARY_ID  => $row[0],
                -DBNAME      => 'SGD_GENE',
                -DISPLAY_ID  => $row[4],
                -DESCRIPTION => $row[15],
                -INFO_TYPE   => 'NONE',
            );
            
            my $dbRef = $dea->_check_external_db($new_display_xref,0);
            my $xref_id = $dea->_store_or_fetch_xref($new_display_xref,$dbRef);
            $new_display_xref->dbID($xref_id);
            $gene->display_xref($new_display_xref);
            
            $gene->display_xref($new_display_xref);
            
            $new_name = 1;
        }
    }
    
    my $new_description = $row[15] . ' [Source:SGD;Acc:' . $row[0] . ']';
    if ( $gene->description() ne $new_description ) {
        $logger->info( "New Description " . $row[3] . " - " . $new_description);
        $gene->description($new_description);
        
        $new_description = 1;
    }
    
    if ( $new_name eq 1 or $new_description eq 1 ) {
        $gene_adaptor->update($gene);
    }
}
