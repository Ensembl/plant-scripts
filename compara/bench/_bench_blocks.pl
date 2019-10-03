use strict;
use warnings;

if(!$ARGV[1]){
	die "# usage: $0 <block_file> <gtf_file.gz>\n"; 
}

my ($block_file, $gtf_file) = @ARGV;

my ($bkid, $first_gene, $last_gene, $geneid);
my ($found, $total, $ids, @data);

# read pre-defined blocks of synthenic gnes
open(BLOCK, "<", $block_file) || die "# ERROR: cannot read $block_file\n";
while(<BLOCK>){
	#A       AT1G01010       AT1G19840
	if(/^(\S+)\s+(\S+)\s+(\S+)/){
		($bkid, $first_gene, $last_gene) = ($1, $2, $3);

		print "$bkid [$first_gene,$last_gene]";

		($total, $found, $ids) = (0, 0, '');
		open(GTF, "gzip -dc $gtf_file |") || die "# ERROR: cannot read $gtf_file\n";
		while(<GTF>){
			next if(/^#/);
			#1       araport11       gene    3631    5899    .       +       .       gene_id "AT1G01010";
			@data = split(/\t/,$_);
			
			next if($data[2] ne 'gene');
	
			$geneid = '';	
			if($data[8] =~ m/gene_id "([^;"]+)/){ $geneid = $1 }

			if($geneid eq $first_gene){ $found = 1 }

			if($found){ 
				$ids .= "\t$geneid";
				$total++;	
			}

			if($geneid eq $last_gene){
				print "\t$total$ids\n";
				last;
			}
		}
		close(GTF);

		#if($total == 0){
		#	print "$total\n";
		#}	
	}
}
close(BLOCK);
