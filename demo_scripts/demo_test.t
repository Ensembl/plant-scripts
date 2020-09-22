use strict;
use warnings;
use Test::More tests => 5;

ok( eval{ `bash exampleFTP.sh --spider 2>&1` } !~ /No such file/ , 'exampleFTP.sh --spider ' );

ok( eval{ `bash exampleMySQL.sh` } =~ /COUNT/ , 'exampleMySQL.sh' );

ok( eval{ `perl exampleAPI.pl` } =~ /Cadenza1734/ , 'exampleAPI.pl' );

ok( eval{ `Rscript exampleBiomart.R` } =~ /IWGSC/ , 'exampleBiomaRt.R' );

ok( eval{ `python exampleREST.py` } =~ /TraesCS4B02G042700.1\t812/ , 'exampleREST.py' );
