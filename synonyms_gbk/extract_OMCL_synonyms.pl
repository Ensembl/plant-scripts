#!/usr/bin/perl -w

use strict;
use warnings;
use File::Basename;
use Bio::SeqIO;

my $VERBOSE = 0;
my $MINANNOTCHARS = 3;

# reads i) a FASTA file with query sequences with full headers starting with a unique integer and 
# ii) a folder with get_homologues clusters including sequences from i)
#
# Produces a TSV file with gene synonyms based on those clusters
#
# Bruno Contreras Moreira 2019

die "# usage: $0 <query .fna.gz> <cluster dir>\n" unless($ARGV[1]);
 
my ($fastafile, $clusterdir) = @ARGV;

my $fastabasename = basename($fastafile);
my ($queryid, $subjectid, $geneid, $acc, $full);
my (%annot);

# 1) read query sequence full FASTA headers
open(FASTA,"zcat $fastafile |") || die "# cannot read $fastafile\n";
while(<FASTA>){
	if(/^>(\d+)\|([^\|]+)\|([^\|]*)/){
		($queryid, $acc, $full) = ($1,$2,$3);
	
		next if(length($full) < $MINANNOTCHARS || 
			$full =~ m/predicted/); # what's the point of adding predicted stuff?

		$annot{$queryid}{'acc'} = $acc;
		$annot{$queryid}{'full_name'} = $full || 'NA';
	}	
}
close(FASTA);

# 2) read clusters of sequences (OMCL)
opendir(DIR,$clusterdir) || die "# ERROR: cannot list $clusterdir\n";
my @clusters = grep{/\.fna/} readdir(DIR);
closedir(DIR);

foreach my $clustfile (@clusters){

	my (%gene, %CDS, %clust_annot);

	open(FNA,"<","$clusterdir/$clustfile") ||	die "# cannot read $clusterdir/$clustfile\n";
	while(<FNA>){

		#>HORVU4Hr1G080410.7 [Hordeum_vulgare.IBSC_v2.cdna.all.fa.gz] | aligned:1-638 (638)
		#>18906|BAJ93979.1|predicted [barley_cds.fna.gz] | aligned:1-648 (648)

		next if not(/^>/);

		# is this a query sequence (a CDS)
		if(/\[$fastabasename\]/){
		
			if(/^>(\d+)/){ # only 1st integer is wanted 
				$queryid = $1;
				$CDS{$queryid}++;
			}
			else{ next }

		} else { # this a ref transcript 
		
			if(/^>(\S+)/){
            $subjectid = $1;
				$geneid = $subjectid;
				$geneid =~ s/\.\d+$//;
            $gene{$geneid}++;
         }		
		}
	}
	close(FNA);

	# skip clusters with a mix of genes
	if(scalar(keys(%gene)) > 1){
		print "# skip mixed cluster: $clusterdir/$clustfile ".join(',',keys(%gene))."\n" if($VERBOSE);
		next;
	}

	# skip clusters with divergent annotations
	foreach $queryid (keys(%CDS)){

		next if(!defined($annot{$queryid}{'acc'}) || !defined($annot{$queryid}{'full_name'}));

		$clust_annot{$annot{$queryid}{'full_name'}}++;
	}

	if(scalar(keys(%clust_annot)) > 1){
		print "# skip mixed annotations: $clusterdir/$clustfile ".join(',',keys(%clust_annot))."\n" if($VERBOSE);
      next;
	}

	# print TSV synonyms
	foreach $queryid (keys(%CDS)){

		next if(!defined($annot{$queryid}{'acc'}) || !defined($annot{$queryid}{'full_name'}));

		printf("%s\t%s\t%s",
         $geneid, $annot{$queryid}{'acc'}, $annot{$queryid}{'full_name'});
		
		if($VERBOSE){ print "\t$queryid" }

		print "\n";
	}		
}

