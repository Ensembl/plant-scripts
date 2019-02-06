#!/usr/bin/perl

# adapted from 
# https://ben.lobaugh.net/blog/202527/mysql-who-see-who-is-connected-to-your-mysql-server-with-this-script
# by Bruno Contreras EMBL-EBI 2019

use strict;
use warnings;
use DBI;
use 5.010;
use Getopt::Long;

# Set defaults
my $host = 'localhost';
my $db = 'information_schema';
my $user = ''; 
my $pass = ''; 

# Get commandline options
GetOptions (
    "host=s"    => \$host,
    "db=s"      => \$db,
    "user=s"    => \$user,
    "pass=s"    => \$pass,
) or die ( "Valid options\n--host\n--user\n--password\n--db database" );

# Create database connection
my $dbh = DBI->connect( "DBI:mysql:$db:$host", $user, $pass, {'PrintError'=>0} )
        or die "** Connection error!\nTry different options:\n--host\n--user\n--password\n--db database\n\n\n";

# Get the current processes to check user
my $sql = "show processlist";
my $sth = $dbh->prepare($sql);
$sth->execute() or die "SQL error: $DBI::errstr\n";

# Create an array of all the users
my %user;
while (my @row = $sth->fetchrow_array()) {
    $user{$row[1]}++;
}

# Display User details
print "MySQL - Logged in users\n";
print "=" x 15, "\n";

for my $x (keys %user) {
    print "$x - $user{$x} connections\n";
}


print "\n";
