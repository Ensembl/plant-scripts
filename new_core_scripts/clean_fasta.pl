#!/usr/bin/env perl
## Cleans headers in fasta file leaving only chr number or name
## NOTE: edit the regex below if needed

use strict;
use warnings;
use feature qw/say/;

{
    my ($fasta_file) = @ARGV;
    if (@ARGV != 1){
        usage();
    }
    open IN, "<", $fasta_file or die "can't open $fasta_file\n";
    while (my $line = <IN>){
        chomp($line);
        if ($line =~ /^>chr(.*?)__/ || # wheat
			  $line =~ /^>Chr(\S+)/ || $line =~ /^>chromosome\s*(\S+)/ #generic
			  ){
            $line = ">$1";
        }
        say $line;
    }
}

sub usage {
    say "Usage perl parse_fasta [fasta_file]";
    exit 0;
}
 

