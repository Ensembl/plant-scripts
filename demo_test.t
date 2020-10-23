use strict;
use warnings;
use Test::More;

my $number_of_tests = 6;

ok( eval{ `./exampleFTP.sh --spider test 2>&1` } =~ /Brachypodium_distachyon/ , 
	'exampleFTP.sh' );

ok( eval{ `./exampleMySQL.sh test` } =~ /_core_/ , 'exampleMySQL.sh' );

ok( eval{ `python exampleREST.py test` } =~ /hordeum_vulgare/ , 'exampleREST.py' );

ok( eval{ `perl exampleREST.pl test` } =~ /hordeum_vulgare/ , 'exampleREST.pl' );

ok( eval{ `Rscript exampleREST.R test` } =~ /hordeum_vulgare/ , 'exampleREST.R' );

ok( eval{ `perl exampleAPI.pl test` } =~ /xref/ , 'exampleAPI.pl' );

if($ARGV[0] && $ARGV[0] eq 'biomart'){
	ok( eval{ `Rscript exampleBiomart.R test` } =~ /IWGSC/ , 'exampleBiomaRt.R' );
	$number_of_tests++;
}

done_testing( $number_of_tests );
