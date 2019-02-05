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
    open IN, "<", $file or die "can't open $file\n";

    ##This is the starting chrom header, everytime it changes total_length becomes 0
    my $current_chr_header = 1;
    
    my ($chr_header,$chunk_header,$count,$chunk_length,$chunk_count,$total_length);
    while (my $line = <IN>){
        chomp($line);
        if ($line =~ />(\d+)_(\d+)/){
            if ($chunk_length){
                
                ##Get the coordinates for the start and end of the ASM part
                my $asm_start = $total_length+1;
                my $asm_end   = $total_length+$chunk_length;
                
                ##Print the output for the AGP file
                print "$chr_header\t$asm_start\t$asm_end\t$chunk_count\tW\t";
                say "$chunk_header\t1\t$chunk_length\t+";
                $total_length = $total_length + $chunk_length;
            }
            $chunk_length = 0;
            $chr_header   = $1;
            $chunk_count  = $2;
            
            ##If the chrom header changes, reset the total length
            if ($chr_header != $current_chr_header){
                $total_length = 0;
                $current_chr_header = $chr_header;
            }

            ##The new chunk header for the new chunk
            $chunk_header = $chr_header."_".$chunk_count;
            next;
        }
        $chunk_length += length($line);
    } 
    
    ##Printout for the last chunk
    ##Get the coordinates for the start and end of the ASM part
    my $asm_start = $total_length+1;
    my $asm_end   = $total_length+$chunk_length;
    
    ##Print the output for the AGP file
    print "$chr_header\t$asm_start\t$asm_end\t$chunk_count\tW\t";
    say "$chunk_header\t1\t$chunk_length\t+";
    $total_length = $total_length + $chunk_length;

}


sub usage {
    say "Usage perl fasta2agp.pl [a] [b]";
    exit 0;
}
 
