use strict;
use warnings;
use Test::More;

my $number_of_tests = 1;

ok( eval{ `recipes/exampleMySQL.sh test` } =~ /_core_/ , 'exampleMySQL.sh' );

if(!$ARGV[0] || $ARGV[0] ne 'travis'){

	# FTP/REST/API tests might timeout from Travis
	
	$number_of_tests += 5;
	    
	ok( eval{ `recipes/exampleFTP.sh --spider test 2>&1` } =~ /Brachypodium_distachyon/ ,
		'exampleFTP.sh' );

	ok( eval{ `python recipes/exampleREST.py test` } =~ /hordeum_vulgare/ , 'exampleREST.py' );

	ok( eval{ `perl recipes/exampleREST.pl test` } =~ /hordeum_vulgare/ , 'exampleREST.pl' );

	ok( eval{ `Rscript recipes/exampleREST.R test` } =~ /hordeum_vulgare/ , 'exampleREST.R' );

	ok( eval{ `perl recipes/exampleCRAM.pl test` } =~ /subgroup/ , 'exampleCRAM.pl' );
}

# requires perl API to be installed ie 'make install_ensembl'
if($ARGV[0] && $ARGV[0] eq 'API'){
	ok( eval{ `perl recipes/exampleAPI.pl test` } =~ /xref/ , 'exampleAPI.pl' );
    $number_of_tests++;
}

# requires BiomaRt R library ie 'make install_biomart_r'
if($ARGV[0] && $ARGV[0] eq 'biomart'){
	ok( eval{ `Rscript recipes/exampleBiomart.R test` } =~ /IWGSC/ , 'exampleBiomaRt.R' );
	$number_of_tests++;
}

done_testing( $number_of_tests );
