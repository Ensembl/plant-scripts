#!/usr/bin/perl

# adapted from 
# https://ben.lobaugh.net/blog/202527/mysql-who-see-who-is-connected-to-your-mysql-server-with-this-script
# by Bruno Contreras EMBL-EBI 2019

use strict;
use warnings;
use Getopt::Long;
use DBI;

# Set defaults
my $db = 'information_schema';
my ($host,$user,$port,$pass,$help,$host_params) = ('','','',''); 

# Get commandline options
GetOptions (
	"help|?" => \$help,
	"host=s" => \$host,
) || help_message(); 

if($help){ help_message() }

sub help_message {
	print "\nusage: $0 [options]\n\n".
	"--host       (required, example: --host mysql-eg-hive\n";
	exit(0);
}

if($host eq '') {
	die "# ERROR: please indicate a valid mysql host\n"
}

chomp( $host_params = `$host details script` );
if($host_params =~ m/--host (\S+) --port (\d+) --user (\S+) --pass (\S+)/){
	($host,$port,$user,$pass) = ($1,$2,$3,$4);
}
elsif($host_params =~ m/--host (\S+) --port (\d+) --user (\S+)/){
        ($host,$port,$user) = ($1,$2,$3);
}

# Create database connection
my $dbh = DBI->connect( "DBI:mysql:$db:$host:$port", $user, $pass, {'PrintError'=>0} ) ||
	die "# ERROR: cannot connect to $host : $DBI::errstr\n";

# Get the current processes to check user
my $sql = "show processlist";
my $sth = $dbh->prepare($sql);
$sth->execute() || die "# SQL ERROR: $DBI::errstr\n";

# Create an array of all the users
my %user;
while (my @row = $sth->fetchrow_array()) {
    $user{$row[1]}++;
}

# Display User details
print "# Users logged in $host:\n";
for my $x (keys %user) {
    print "$x - $user{$x} connections\n";
}
