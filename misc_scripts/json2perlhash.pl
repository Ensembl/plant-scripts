#!/usr/bin/perl
use strict;
use warnings;
use JSON;
use Data::Dumper;

# I used to convert webdata from productiondb 

if(!$ARGV[0]){ die "# usage: $0 'json string'\n" }

my $json = $ARGV[0];
print  Dumper(decode_json($json));
