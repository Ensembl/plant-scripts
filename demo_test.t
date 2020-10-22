use strict;
use warnings;
use Test::More tests => 5;

ok( eval{ `bash exampleFTP.sh --spider 2>&1` } !~ /No such file/ , 'exampleFTP.sh --spider ' );

ok( eval{ `bash exampleMySQL.sh test` } =~ /stable_id/ , 'exampleMySQL.sh' );

ok( eval{ `Rscript exampleBiomart.R test` } =~ /IWGSC/ , 'exampleBiomaRt.R' );

ok( eval{ `python exampleREST.py test` } =~ /hordeum_vulgare/ , 'exampleREST.py' );

ok( eval{ `perl exampleREST.pl test` } =~ /hordeum_vulgare/ , 'exampleREST.pl' );
