use strict;
use warnings;
use feature qw/say/;
use Data::Dumper;
use lib '/nfs/production/panda/ensemblgenomes/development/gnaamati/lib';
use FileReader qw(slurp read_file file2hash file2hash_tab);
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception;

{
    my ($names) = @ARGV;
    if (@ARGV < 1){
        usage();
    }
    #my $reg_file = '/homes/gnaamati/registries/prod3_gn.reg';
    my $reg_file = '/homes/gnaamati/registries/prod3_3b.reg';

    my $slice_adaptor = get_slice_adaptor($reg_file);
    my @names = slurp($names);
    
    my $count = 0;
    for my $name (@names){
        say ">$name";
        my $slice = $slice_adaptor->fetch_by_region('toplevel',$name);
        #say $slice;
        #my $seq_region_id = $slice_adaptor->get_seq_region_id($slice);
        say $slice->seq;

        #say "seq_region id = $seq_region_id\n";

        if ($count++ == 100){
            die;
        }
    }
}

#========================================
sub get_slice_adaptor {
#========================================
    my ($reg_file) = @_;
    Bio::EnsEMBL::Registry->load_all($reg_file);
    my $species = 'triticum_aestivum';
    my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'Core', 'Slice');
    return $slice_adaptor;
}

sub usage {
    say "Usage perl parse_exons.pl [a] [b]";
    exit 0;
}

