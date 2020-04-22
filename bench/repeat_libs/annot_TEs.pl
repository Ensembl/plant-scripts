#!/usr/bin/perl 
use strict;
use warnings;
use File::Copy;

# 1) align nucleotides sequences in clusters and annotate any Pfam domains;
# note that original FASTA headers are shortened and stored in $short_dir
# with lower-case sequences

my $align_dir = 'TE_alignments';
my $short_dir = 'TE_clusters'; # shortened FASTA headers

my %TElibs = (
	'SINEs.plants.fna.gz.nucl'         => 1,
	'mipsREdat_9.3p_ALL.fasta.gz.nucl' => 1,
	'repetDB.Mar2020.fna.gz.nucl'      => 1,
	'trep-db_nr_Rel-19.fasta.gz.nucl'  => 1,
	'rice6.9.5.liban.fna.gz.nucl'      => 1,
	'maizeTE11122019.fna.gz.nucl'      => 1,
	'SoyBaseTE.fna.gz.nucl'            => 1,
   'TAIR10_TE.fna.gz.nucl'            => 1
);

my $annotEXE = '~/soft/github/get_homologues/annotate_cluster.pl';

my $panmatfile;
if(!defined($ARGV[0])){ die "# usage: $0 <pangenome_matrix_genes_t0.tr.tab>\n" }
else{ 
	$panmatfile = $ARGV[0];
}

my (%TEcolumns);
my ($cluster_folder,$cluster_name,$col,$last_col,$TEok,$num_seqs);
my ($seqid,$taxon);

open(MAT,"<",$panmatfile) || die "# error: cannot read $panmatfile\n";
while(<MAT>){
	chomp;

	if(/^source:(\S+)/){
		$cluster_folder = "$1/";
		my @names = split(/\t/,$_);
		$last_col = $#names;
		foreach $col ( 1 .. $last_col ){
			if($TElibs{ $names[$col] }){
				$TEcolumns{ $col } = 1; print "# TE column $col\n";
			}
		}
	} else {

		# check this clusters composition, does it include TEs?
		$TEok = $num_seqs = 0;
		my @data = split(/\t/,$_);
		$cluster_name = $data[0];

		#next if($cluster_name !~ '539866_'); #debug

		foreach $col ( 1 .. $last_col ){

			if($data[$col] eq '-'){ next }

			my @seqids = split(/,/,$data[$col]);
			$num_seqs += scalar(@seqids);		

			if($TEcolumns{ $col }){ 
				$TEok++;
			} 
      }

		next if($TEok == 0); # skip clusters without TE sequences

		# shorthen FASTA headers
		my $short_cluster_name = "$short_dir/$cluster_name";
		open(SHORT,">",$short_cluster_name) || die "# error: cannot write to $short_cluster_name\n";
		
		open(CLUSTER,"<",$cluster_folder.$cluster_name) ||
			die "# error: cannot read $cluster_folder$cluster_name\n";
		while(my $line = <CLUSTER>){
			if($line =~ m/^>(\S+).*?\[(\S+?)\]$/ || $line =~ m/^>(\S+).*?\[(\S+?)\]\s\|/ ||
				$line =~ m/^>(\S+).*?\[(.+?)\]\s\|/ || $line =~ m/^>(\S+).*?\[(.+?)\]/)
			{
				($seqid, $taxon) = ($1,$2);
				# repetDB
				# PotentialHostGene-DHX-RG_VvTEdenovoV2-B-R17515-Map5
				if($seqid =~ /PotentialHostGene.*?([^\-_]+-[^\-]+-[^\-]+-Map\d+)/){ 
					$seqid = "PHG_$1"; 
				}
				elsif($seqid =~ m/\S+?_([^\-_]+-[^\-]+-[^\-]+-Map\d+)/){ $seqid = $1 }
				else{ $seqid = (split(/\|/,$seqid))[0] } # mips
				print SHORT ">$seqid [$taxon]\n";
			}
			else{ print SHORT lc($line) }
		}
		close(CLUSTER);

		close(SHORT);

		# align and annotate TE sequences
		if($num_seqs > 1) {
			system("$annotEXE -f $short_dir/$cluster_name -D -o $align_dir/$cluster_name");
		} else {
			copy("$short_dir/$cluster_name" , "$align_dir/$cluster_name");				
		}

		print "$_\t$num_seqs\n"; #exit;
	}	
}
close(MAT);

