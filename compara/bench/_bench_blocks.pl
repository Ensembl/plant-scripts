use strict;
use warnings;

# Bruno Contreras Moreira 2019

if(!$ARGV[2]){
	die "# usage: $0 <block_file> <gtf_file.gz> <synths TSV file> \n"; 
}

my ($block_file, $gtf_file_gz, $synths_file) = @ARGV;

# parse genes spanning blocks and return ref to hash
my ($ref_blocks, $ref_gene2block) = populate_gene_blocks_GTF( $block_file, $gtf_file_gz );

# loop through list of synthelogs and compute stats
my ($ref_gene, $c, $bkid);
my (@species, %stats);
open(TSV, "<", $synths_file) || die "# cannot read $synths_file\n";
while(<TSV>){
	next if(/^[#\s+]/);	
	chomp;
	my @col = split(/\t/,$_);

	#gene    location        arabidopsis_lyrata      brassica_rapa   ....
	if($col[0] eq 'gene'){ 
		foreach $c (2 .. $#col){ $species[$c] = $col[$c] }
	} 
	else { #AT1G01070       1:38444-41017:- scaffold_100145.1       NA      NA      AT1G01070.1     NA

		# find out to which block(s) and species this (reference genome) synthelog belongs
		$ref_gene = $col[0];
		if($ref_gene2block->{$ref_gene}){
			foreach $bkid (@{ $ref_gene2block->{$ref_gene} }){

				foreach $c (2 .. $#col){
					if($col[$c] ne 'NA'){
						$stats{$species[$c]}{$bkid}++;
					}
				}
			}
		} else {
			$bkid = 'noblock'; # not in predefined block
         
			foreach $c (2 .. $#col){
         	if($col[$c] ne 'NA'){
            	$stats{$species[$c]}{$bkid}++;
            }
         }
		}
	}

}
close(TSV);

# print stats summary
my @allblocks = ( keys(%$ref_blocks) , 'noblock'  );

# header
print "block\tinterval\tgenes\t%synthelogs";
foreach $c (2 .. $#species){ print "\t$species[$c]" }
print "\n";

foreach $bkid ( @allblocks  ){
	printf("%s\t%s\t%d",
		$bkid,$ref_blocks->{$bkid}{'interval'} || '[NA,NA]',
		$ref_blocks->{$bkid}{'total'} || 0 );

	# all synthelogs including arabidopsis_thaliana
	if( $ref_blocks->{$bkid}{'total'} ){
	   printf("\t%1.1f",100*$stats{'arabidopsis_thaliana'}{$bkid} / $ref_blocks->{$bkid}{'total'}); 
	} else {
		print "\tNA";
	}

	foreach $c (2 .. $#species){
		printf("\t%d",$stats{$species[$c]}{$bkid} || 0); 
	}

	print "\n";
}

##########################################


sub populate_gene_blocks_GTF {

	my ($block_file, $gtf_file) = @_;

	my ($bkid, $first_gene, $last_gene, $geneid);
	my ($found, $total, @data);
	my (%blocks, %member);

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

				if($found){ 
					
					# add sorted gene names to current interval
					push(@{ $blocks{$bkid}{'genes'} }, $geneid ); 

					# store to which block a random gene belongs
					push(@{ $member{$geneid} }, $bkid);
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

	return (\%blocks, \%member);
}	
