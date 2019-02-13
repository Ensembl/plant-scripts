#!/usr/bin/env perl
use strict;
use warnings;
use feature qw/say/;

# For debugging
use Data::Dumper;
use Time::Piece;

use Getopt::Long;

# Connect to Ensembl and friends
use Bio::EnsEMBL::Registry;

use Bio::EnsEMBL::Utils::Exception;

# Mapping between object type and table name

my $reg_conf;
my $species;
my $db_type; #e.g. "core"
my @feature_types; #e.g. "gene, transcript"
my $sql_out_dir;
my $logfile;

GetOptions("reg_conf=s" => \$reg_conf,
	   "species=s" => \$species,
	   "db_type=s" => \$db_type,
);


Bio::EnsEMBL::Registry->load_all($reg_conf);

print "Getting a database adaptor ($db_type:$species)\n";
my $db_adaptor = Bio::EnsEMBL::Registry->
  get_DBAdaptor($species, $db_type)
    or die "could not get db adaptor for sp: $species db_type: $db_type\n$!\n";

#print $db_adaptor->dbname;
#
@feature_types = ('exon');


foreach my $feature_type (@feature_types) {
    my $f_adaptor = $db_adaptor->get_adaptor($feature_type)
        or die "couldn't get feature adaptor for a $feature_type\n$!\n";#
        #
    #  print LOGFILE "Counting features ($feature_type)\n" if ($logfile);
    #  print LOGFILE "There are ", $f_adaptor->generic_count, " ", $feature_type, "s\n" if ($logfile);
    
    my $items = $f_adaptor->fetch_all();
    #my $exon = $f_adaptor->fetch_by_stable_id('Traes_7BL_4955BA7E7.E4');
    #my $exon = $f_adaptor->fetch_by_stable_id('Traes_4DS_4BD85B5C7.E17');

    #print ">", $exon->stable_id, "\n";
    #print $exon->seq->seq, "\n";
    #die('');
    
    #$items = [$exon];
    my $count = 0;
    ##Get all exons
    while (my $exon = shift @{$items}) {
        my $sid = $exon->stable_id;
        if ($sid !~ /^TRAES3B/){
            next;   
        }
        print ">", $sid, "\n";
        print $exon->seq->seq, "\n";
        #if ($count++ > 100){
        #    die('');
        #}
    }
    die('');
}
