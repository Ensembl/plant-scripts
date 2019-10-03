use strict;
use warnings;

if(!$ARGV[1]){
	die "# usage: $0 <block_file> <gtf_file.gz>\n"; 
}

my ($block_file, $gtf_file_gz) = @ARGV;

my $ref_blocks = populate_gene_blocks_GTF( $block_file, $gtf_file_gz );

foreach my $blockid (keys(%$ref_blocks)){
	print "$blockid $ref_blocks->{$blockid}{'interval'} $ref_blocks->{$blockid}{'total'}\n";
}


sub populate_gene_blocks_GTF {

	my ($block_file, $gtf_file) = @_;

	my ($bkid, $first_gene, $last_gene, $geneid);
	my ($found, $total, @data, %blocks);

	# reads pre-defined blocks of synthenic genes and computes list of genes for each block
	open(BLOCK, "<", $block_file) || die "# ERROR: cannot read $block_file\n";
	while(<BLOCK>){
		#A       AT1G01010       AT1G19840
		if(/^(\S+)\s+(\S+)\s+(\S+)/){
			($bkid, $first_gene, $last_gene) = ($1, $2, $3);

			# start a new block
			$blocks{$bkid}{'interval'} = "[$first_gene,$last_gene]";

			($total, $found) = (0, 0);
			open(GTF, "gzip -dc $gtf_file |") || die "# ERROR: cannot read $gtf_file\n";
			while(<GTF>){
				next if(/^#/);
				#1       araport11       gene    3631    5899    .       +       .       gene_id "AT1G01010";...
				@data = split(/\t/,$_);
				next if($data[2] ne 'gene');
	
				$geneid = '';	
				if($data[8] =~ m/gene_id "([^;"]+)/){ $geneid = $1 }

				if($geneid eq $first_gene){ $found = 1 }

				if($found){ # add sorted gene names to current interval
					push(@{ $blocks{$bkid}{'genes'} }, $geneid );
					$total++;	
				}

				if($geneid eq $last_gene){ 
					$blocks{$bkid}{'total'} = $total;
					last;
				}
			}
			close(GTF);
		}
	}
	close(BLOCK);

	return \%blocks;
}	
