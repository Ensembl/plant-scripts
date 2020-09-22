#!/bin/env perl
use strict;
use warnings;

# add genes to a GFF file containing only mRNA features with no explicit parent
# by Bruno Contreras Moreira EMBL-EBI 2019-20

if(!$ARGV[0]){ die "# usage: $0 <GFF>\n" }

my ($tparent, $tname, $gene, $transcript, $other, $type, $id, %seen);

open(GFF,"<",$ARGV[0]) || die "# ERROR: cannot read $ARGV[0]\n";
while(<GFF>){

	#CcLG11  GLEAN   mRNA    9854    10924   0.575745        +       .       ID=C.cajan_00001;evid_id=C.cajan_GLEAN_10046622;
	#CcLG11  GLEAN   CDS     9854    10924   .       +       0       Parent=C.cajan_00001;
	#CcLG11  GLEAN   mRNA    20087   31342   0.0164527       -       .       ID=C.cajan_00004;evid_id=C.cajan_GLEAN_10046625;
	#CcLG11  GLEAN   CDS     31146   31342   .       -       0       Parent=C.cajan_00004;
	#CcLG11  GLEAN   CDS     30427   30541   .       -       1       Parent=C.cajan_00004;
	#CcLG11  GLEAN   CDS     29974   30238   .       -       0       Parent=C.cajan_00004;
	#CcLG11  GLEAN   CDS     29375   29848   .       -       2       Parent=C.cajan_00004;

	next if(/^#/);

	my @col = split(/\t/,$_);

	if(defined($col[2]) && $col[2] eq 'mRNA'){
		
		$gene = $_;
		$transcript = $_;

		if($col[8] =~ m/ID=([^;]+);/){

			# create a new gene
			$tparent = $1;
			$tname = $tparent.'.t';
			$gene =~ s/\tmRNA\t/\tgene\t/;	
			

			# first print the new gene feature
			if(!$seen{$tparent}){
				print $gene;
				$seen{$tparent}++;
			}

			# now update and print transcript
			$transcript =~ s/ID=[^;]+;/ID=$tname;Parent=$tparent;/;
			print $transcript;
		}
	
	} else {
		$other = $_;
		chomp($other);
		$type = lc($col[2]);
	
		# update parent
		$other =~ s/Parent=[^;]+(;*)/Parent=$tparent.t$1\n/;

		# add ID if missing
		if($other !~ /ID=/){
			#count feature
			$seen{"$tparent.t.$type"}++;
			$id = "$tparent.t.$type".$seen{"$tparent.t.$type"};
			$other =~ s/Parent=/ID=$id;Parent=/;
		}

		print $other;
	}

}
close(GFF);
