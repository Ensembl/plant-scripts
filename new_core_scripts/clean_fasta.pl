##Parses a fasta file for various purposes
use strict;
use warnings;
use feature qw/say/;
use Data::Dumper;
use lib '/nfs/production/panda/ensemblgenomes/development/gnaamati/lib';
use FileReader qw(slurp read_file file2hash file2hash_tab);
use constant {
    TRUE  => 1,
    FALSE => 2,
};

#my $fasta_dir = '/nfs/nobackup/ensemblgenomes/grabmuel/ensgen-sequences/ftp';

{
    my ($fasta_file) = @ARGV;
    if (@ARGV != 1){
        usage();
    }
    open IN, "<", $fasta_file or die "can't open $fasta_file\n";
    while (my $line = <IN>){
        chomp($line);
        if ($line =~ /^>chr(.*?)__/){
            $line = ">$1";
        }
        say $line;
    }
}

sub usage {
    say "Usage perl parse_fasta [fasta_file]";
    exit 0;
}
 

