use 5.14.0;
use warnings;
use Data::Dumper;
use lib '/nfs/production/panda/ensemblgenomes/development/gnaamati/lib';
use FileReader qw(slurp slurp_hash_list read_file file2hash file2hash_tab line2hash);

{
    my ($file) = @ARGV;
    if (@ARGV < 1){
        usage();
    }
    my $header;
    my $header_count;
    my $count;
    my $print = 0;
    open IN, "<", $file or die "can't open $file\n";
    while (my $line = <IN>){
        chomp($line);
        
        ##Create a new chunk each time there is a header
        if ($line =~ /(TuUngrouped_contig_\d+)/){
            $header = ">$1";
            ##Print just the ungrouped contigs
            say $header;
            $print = 1;
            next;
        }
        if ($print){ 
            say $line;
        }
    }
}

sub usage {
    say "Usage perl break_fasta.pl [a] [b]";
    exit 0;
}
 
