use strict;
use warnings;
use Test::More tests => 5;

ok( eval{ `bash exampleFTP.sh --spider test 2>&1` } =~ /Brachypodium_distachyon/ , 
	'exampleFTP.sh' );

ok( eval{ `bash exampleMySQL.sh test` } =~ /COUNT/ , 'exampleMySQL.sh' );

ok( eval{ `python exampleREST.py test` } =~ /hordeum_vulgare/ , 'exampleREST.py' );

ok( eval{ `perl exampleREST.pl test` } =~ /hordeum_vulgare/ , 'exampleREST.pl' );

ok( eval{ `perl exampleAPI.pl test` } =~ /xref/ , 'exampleAPI.pl' );

#takes too long and BiomaRt package not maintained by Ensembl
#ok( eval{ `Rscript exampleBiomart.R test` } =~ /IWGSC/ , 'exampleBiomaRt.R' );
