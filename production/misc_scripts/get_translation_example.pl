#!/usr/bin/env perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use Bio::SeqIO;
use feature qw /say/;

##Use the relevant reg conf
my $reg_conf = '/homes/gnaamati/registries/prod2.reg';
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_all($reg_conf);

##fetch a gene by its stable identifier
#my $transcript_adaptor = $registry->get_adaptor("triticum_aestivum", "otherfeatures", "transcript");
#my $transcript_adaptor = $registry->get_adaptor("triticum_dicoccoides", "core", "transcript");

##Get adaptor for core
my $transcript_adaptor = $registry->get_adaptor("triticum_aestivum", "core", "transcript");

#my $gene = $gene_adaptor->fetch_by_stable_id('');
#my $trans = $transcript_adaptor->fetch_by_dbID(2);
#my $trans = $transcript_adaptor->fetch_by_stable_id('TRIDC1AG012410.5');
my $trans = $transcript_adaptor->fetch_by_stable_id('TraesCS4A02G403700.1');
#print Dumper $trans;

my $temp = $trans->translate();
say $temp->seq;
#say "\n\n";

#say $trans->seq->translate->seq;


#my $exons = $trans->get_all_Exons;

#for my $e (@$exons){
#    say "e\n",$e->seq->seq;
#}
#say $trans->seq->translateable_seq;
#print $transcript->translateable_seq(), "\n";

#my $translation = $trans->translation;
#print $translation->seq, "\n";


