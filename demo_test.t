use strict;
use warnings;
use Test::More;

my $number_of_tests = 7;

ok( eval{ `recipes/exampleFTP.sh --spider test 2>&1` } =~ /Brachypodium_distachyon/ , 
	'exampleFTP.sh' );

ok( eval{ `recipes/exampleMySQL.sh test` } =~ /_core_/ , 'exampleMySQL.sh' );

ok( eval{ `python recipes/exampleREST.py test` } =~ /hordeum_vulgare/ , 'exampleREST.py' );

ok( eval{ `perl recipes/exampleREST.pl test` } =~ /hordeum_vulgare/ , 'exampleREST.pl' );

ok( eval{ `Rscript recipes/exampleREST.R test` } =~ /hordeum_vulgare/ , 'exampleREST.R' );

ok( eval{ `perl recipes/exampleAPI.pl test` } =~ /xref/ , 'exampleAPI.pl' );

ok( eval{ `perl recipes/exampleCRAM.pl test` } =~ /subgroup/ , 'exampleCRAM.pl' );

if($ARGV[0] && $ARGV[0] eq 'biomart'){
	ok( eval{ `Rscript recipes/exampleBiomart.R test` } =~ /IWGSC/ , 'exampleBiomaRt.R' );
	$number_of_tests++;
}

done_testing( $number_of_tests );
