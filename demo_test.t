use strict;
use warnings;
use Test::More tests => 5;

ok( eval{ `./exampleFTP.sh --spider test 2>&1` } =~ /Brachypodium_distachyon/ , 
	'exampleFTP.sh' );

ok( eval{ `./exampleMySQL.sh test` } =~ /_core_/ , 'exampleMySQL.sh' );

ok( eval{ `python exampleREST.py test` } =~ /hordeum_vulgare/ , 'exampleREST.py' );

ok( eval{ `perl exampleREST.pl test` } =~ /hordeum_vulgare/ , 'exampleREST.pl' );

ok( eval{ `perl exampleAPI.pl test` } =~ /xref/ , 'exampleAPI.pl' );

#ok( eval{ `Rscript exampleBiomart.R test` } =~ /IWGSC/ , 'exampleBiomaRt.R' );
