#!/usr/bin/env perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use feature qw /say/;
use Data::Dumper;

{
    ##Load the registry for prod3_test
    my $registry = 'Bio::EnsEMBL::Registry';
    #my $reg_conf = 'prod3_test.reg';
    my $reg_conf = '/homes/gnaamati/registries/prod2.reg';
    $registry->load_all($reg_conf);
    
    ##fetch a gene by its stable identifier
    #my $gene_adaptor = $registry->get_adaptor("oryza_sativa", "core", "gene");
    my $gene_adaptor = $registry->get_adaptor("triticum_aestivum_cadenza", "core", "gene");
    my $gene = $gene_adaptor->fetch_by_stable_id('TRIAE_CS42_5AL_TGACv1_375099_AA1216080.2.path1');

    my $transcript = $gene->canonical_transcript;
    my $translation = $transcript->translation;
    say $translation->seq;

}
    #say $transcript->seq->seq;
    #
__END__

    my $slice_adaptor = $registry->get_adaptor("oryza_sativa", "core", "slice");
    #my $slice         = $slice_adaptor->fetch_by_region('toplevel','BA000029.3');
    my $slice         = $slice_adaptor->fetch_by_region('toplevel','X15901.1');

    #my $slice = $transcript->slice;
    #print Dumper $slice;

    my $sub_slice = $slice->sub_Slice(1,100);
    #say $sub_slice->seq;
    say $sub_slice->seq;

    #print Dumper $translation;
    say $translation->seq;

    ##transform the gene
    #my $proj_gene = $gene->transform('toplevel');
    #if (!$proj_gene){
    #    say 'could not project ', $gene->stable_id;
    #}
}

