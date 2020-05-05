use strict;
use warnings;
use Test::More tests => 3;

ok( eval{ `perl exampleREST.pl` } =~ /zea_mays/ , 'exampleREST.pl' );

ok( eval{ `bash exampleFTP.sh` } =~ /# got Compara/ , 'exampleFTP.sh' );

ok( eval{ `bash exampleFTP.sh` } =~ /COUNT/ , 'exampleMySQL.sh' );
