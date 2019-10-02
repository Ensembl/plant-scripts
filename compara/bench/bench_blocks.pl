use strict;
use warnings;

my ($block_file, $gff_file) = @ARGV;

my ($bkid, $first_gene, $last_gene, $geneid);
my ($found, $total, $ids, @data);

open(BLOCK, "<", $block_file) || die "# ERROR: cannot read $block_file\n";
while(<BLOCK>){
	#A       AT1G01010       AT1G19840
	if(/^(\S+)\s+(\S+)\s+(\S+)/){
		($bkid, $first_gene, $last_gene) = ($1, $2, $3);

		print "$bkid [$first_gene,$last_gene]";

		($total, $found, $ids) = (0, 0, '');
		open(GFF, "<", $gff_file) || die "# ERROR: cannot read $gff_file\n";
		while(<GFF>){
			next if(/^#/);
			#1	araport11	gene	2846894	2848100	.	-	.	ID=gene:AT1G08880;
			@data = split(/\t/,$_);
			
			next if($data[2] ne 'gene');
	
			$geneid = '';	
			if($data[8] =~ m/ID=gene:([^;]+)/){ $geneid = $1 }

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
		close(GFF);

		#if($total == 0){
		#	print "$total\n";
		#}	
	}
}
close(BLOCK);
