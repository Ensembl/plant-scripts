use strict;
use warnings;
use Test::More tests => 1;

ok( eval{ `perl ens_single-copy_core_genes.pl -c Oryza -r oryza_sativa` } =~ /# total single-copy core clusters/ , 'ens_single-copy_core_genes.pl' );
