#!/usr/bin/perl
use strict;
use warnings;

# 2) produce a list of Pfam domains of clusters that contain both TEs and cDNAs,
# considered ambiguous sequence clusters

my ($inputlog,$pos_pfam_list,$neg_pfam_list);
if(!$ARGV[2]){ die "usage: $0 <log.annot> <control+> <control->\n" }
else { 
	($inputlog,$pos_pfam_list,$neg_pfam_list) = @ARGV;
}

## parse list of positive Pfam domains, know to be part of TEs
my %pos_pfam;
open(POS,"<",$pos_pfam_list) || die "# ERROR: cannot read $pos_pfam_list\n";
while(<POS>){
	if(/^(PF\d+)/){ $pos_pfam{ $1 } = 1 }
}
close(POS);
printf("# positive Pfam domains = %d\n",scalar(keys(%pos_pfam)));

## parse list of negative Pfam domains, know to be part of non-TE protein coding genes
my %neg_pfam;
open(NEG,"<",$neg_pfam_list) || die "# ERROR: cannot read $neg_pfam_list\n";
while(<NEG>){
   if(/^(PF\d+)/){ $neg_pfam{ $1 } = 1 }
}
close(NEG);
printf("# negative Pfam domains = %d\n",scalar(keys(%neg_pfam)));

## parse input log and compute stats
my $n_of_TE_clusters = 0;
my $n_of_mixed_clusters = 0;
my (%TEcol,%stats,%pfam,$col,$noTE,$TE,$tseqs,$dom,$file,$fullfile,$num_potential_host_gene);

open(LOG,"<",$inputlog) || die "# ERRO: cannot read $inputlog\n";
while(<LOG>){

	# get columns that correspond to TEs
	if(/^# TE column (\d+)/){
		$TEcol{$1} = 1;
	}

	# init stats and parse cluster file 
	if(/^# \S+ -f (\S+)/){
		$fullfile = $1;	
		%pfam = ();
      ($noTE, $TE, $tseqs, $num_potential_host_gene) = (0,0,0,0);

		open(FASTA,"<",$fullfile) || die "# ERROR: cannot read $fullfile\n";
		while(my $line = <FASTA>){
			if($line =~ m/PotentialHostGene/ || $line =~ m/^>PHG/){ $num_potential_host_gene++ }
		}
		close(FASTA);
	}	

	# check number of sequences in cluster
	if(/^# total\s+sequences: (\d+) taxa/){
		$tseqs = $1;
	}

	# parse domains (only if this cluster was 'annotated', 
	# that is, contained > 1 sequence and both TE & cDNA seqs,
	# as opposed to 'not annotated')
	if(/^# Pfam domains: $/){ $pfam{'unknown'}++ } 
	while(/(PF[\d\.]+)/g){ $pfam{$1}++; }

	#if(!$pfam{'PF13041'}) { %pfam = (); $noTE=0; $TE=0; next; } # debugging

	# check cluster composition
	if(/^(\d+_\S+?.fna)/){

		$n_of_TE_clusters++;

		# nothing else to do for clusters 'not annotated'
		next if(scalar(keys(%pfam)) < 1);

		$file = $1;
		my @data = split(/\t/,$_);
		foreach $col (1 .. $#data-2){
			if(defined($TEcol{$col}) && $data[$col] ne '-'){
            $TE++;
         }
			elsif(!defined($TEcol{$col}) && $data[$col] ne '-'){
				$noTE++;
			}
		}

		# this is a mixed cluster, contains both TEs and cDNAs (noTE)
		if($noTE > 0){ 
			foreach $dom (keys(%pfam)){
				$stats{$dom}{'total'} += $pfam{$dom}; # total domain occurrences, redundant
				$stats{$dom}{'seqs'} += $tseqs;
				$stats{$dom}{'noTE'} += $noTE;
				$stats{$dom}{'TE'} += $TE;
				$stats{$dom}{'totalclusters'}++;
				$stats{$dom}{'potential_genes'} += $num_potential_host_gene;
				$stats{$dom}{'clusters'} .= "$file,";
			}

			$n_of_mixed_clusters++;
		}

		%pfam = ();
		($noTE, $TE, $tseqs) = (0,0,0);
	}		
}
close(LOG);

# print overall stats
my ($note, $frac_potgenes);

print "# TEclusters=$n_of_TE_clusters\n";
print "# mixedclusters=$n_of_mixed_clusters\n";

print "domain\ttotclusters\toccurrences\ttotseqs\tTElibs\tcDNAlibs\tpotgenes\tfrac_potgenes\tnotes\tclusters\n";
foreach $dom (sort {$stats{$b}{'totalclusters'}<=>$stats{$a}{'totalclusters'}} keys(%stats)){

	$frac_potgenes = 0;
	if($dom ne 'unknown' && $stats{$dom}{'seqs'} > 0){ 
		$frac_potgenes = $stats{$dom}{'potential_genes'}/$stats{$dom}{'seqs'} 
	} 

	$note = '';
	if(defined($pos_pfam{ $dom })){ $note = 'positive_control' }
	elsif(defined($neg_pfam{ $dom })){ $note = 'negative_control' }

	printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%1.2f\t%s\t%s\n",
		$dom,
		$stats{$dom}{'totalclusters'},$stats{$dom}{'total'},$stats{$dom}{'seqs'},
		$stats{$dom}{'TE'},$stats{$dom}{'noTE'},
		$stats{$dom}{'potential_genes'},
		$frac_potgenes,
		$note,
		$stats{$dom}{'clusters'});
}


