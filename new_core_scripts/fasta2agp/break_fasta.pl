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
    open IN, "<", $file or die "can't open $file\n";
    while (my $line = <IN>){
        chomp($line);
        
        ##Create a new chunk each time there is a header
        if ($line =~ /chromosome (\d)\,/){
            
            ##Header is the chrom number with the chunk count
            $header = ">$1";
            $header_count = 1;
            say $header."_$header_count";
            $count=0;
            next;
        }
        say $line;
        $count++;

        ##Make each chunk 600,000 bases long
        if ($count % 10_000 == 0){
            $header_count++;
            say $header."_$header_count";
        }
    }
}

sub usage {
    say "Usage perl break_fasta.pl [a] [b]";
    exit 0;
}
 
