use strict;
use warnings;
use File::Temp qw/ tempfile /;
use Set::IntSpan::Fast::XS;

# Takes four files and computes stats on their masked, lower case segments.
# Files 1 & 2 are obtained from ftp://ftp.ensemblgenomes.org/pub/plants/release-XX/fasta
# File 3 was masked from the original in 1
#
# 1) FASTA file with control raw genome
# 2) FASTA file with control masked genome 
# 3) FASTA file with experiment masked genome
# 4) GFF   file with mapped genes/features 

# Bruno Contreras Moreira 2019

my $BUFSIZE = 500_000; # in bytes
my $VERBOSE = 0;
my $TEMPRFIX= 'tmp1lineXXXX';
my $FEATTYPE= 'gene';

my $GUNZIPEXE = "gzip -cd ";
# Note hard-masked Ns are converted to lower-case ns
my $FASTA21LINEXE = "perl -ne 'if(!/^>/){ chomp; s/N/n/g; print }' ";
my $FASTA2CHREXE = "perl -ne 'if(/^>/){ print \"\n\" if(\$.>1)} else { chomp; s/N/n/g; print }' ";

###############################################################################

my ($raw_file, $reference_file, $test_file, $feature_file) = @ARGV;

if(!$ARGV[3]){
	die "# usage: $0 <FASTA raw> <FASTA control> <FASTA experiment> <GFF/BED>\n\n";
}

my ($raw_buffer, $raw_read, $ref_buffer, $test_buffer, $ref_read, $test_read);
my ($ord_ref, $ord_test, $ref_lc, $test_lc, $raw_bs);
my ($chr, $pos, $corrpos, $endpos, $off, $feat_overlap);
my ($nb_chr,$N,$mism,$total,$tp,$tn,$fp,$fn) = (0,0,0,0,0,0,0,0);
my ($ref_feat, $test_feat, $ref_mask, $test_mask) = (0,0,0,0);
my (%int2chr,%chr2int,%offset);
my $feat_set = Set::IntSpan::Fast::XS->new();


## 1) read raw FASTA file to compute chr sizes to be used as offsets for BED/GFF

my ($fh_len, $len_filename) = tempfile( $TEMPRFIX , UNLINK => !$VERBOSE );
system("$GUNZIPEXE $raw_file | $FASTA2CHREXE > $len_filename");

open(RAW,"$GUNZIPEXE $raw_file |") || die "# ERROR: cannot read $raw_file\n";
while(<RAW>){
	if(/^>(\S+)/){
		$nb_chr++;
		$chr = $1;
		$int2chr{$nb_chr} = $chr;
		$chr2int{$chr} = $nb_chr;
		$offset{$chr} = 0;
	}
}
close(RAW);

my $accum_offset = 0;
open($fh_len,"<",$len_filename) || die "# ERROR: cannot read $len_filename\n";
while(<$fh_len>){

	chomp;

	# length of current chr is offset for next 
	next if not defined($int2chr{ $. + 1 });
	$chr = $int2chr{ $. + 1 };

	if($VERBOSE){
		printf("# %s %d offset %d\n", $chr, length($_), length($_) + $accum_offset);
	}	

   $offset{ $chr } = length($_) + $accum_offset;
	$accum_offset += length($_);
}
close($fh_len);

## 2) read BED/GFF features and correct coordinates with offsets ###############

open(BEDGFF,"<",$feature_file) || die "# ERROR: cannot read $feature_file\n";
while(<BEDGFF>){
	next if(/^#/);
	#1	araport11	gene	3549920	3550848	.	-	.	ID=gene:AT1G10690;Name=SMR8;biotype=pr...
	my @feat = split(/\t/,$_);

	next if($feat[2] ne $FEATTYPE); # comment if needed

	($chr, $pos, $endpos) = @feat[0,3,4]; 
	$off = $offset{ $chr } - 1; # zero-based 

	$pos += $off;
	$endpos += $off;

	$feat_set->add_range( $pos, $endpos );
}
close(BEDGFF);

## 3) convert FASTA files to single-line raw sequences and compute size in bytes

my ($fh_raw, $raw_filename) = tempfile( $TEMPRFIX , UNLINK => !$VERBOSE );
system("$GUNZIPEXE $raw_file | $FASTA21LINEXE > $raw_filename");

open($fh_raw,"<",$raw_filename) || die "# ERROR: cannot read $raw_filename\n";
my $raw_size = -s $raw_filename;

print "# flat size: $raw_size ($raw_file)\n";

my ($fh_ref, $ref_filename) = tempfile( $TEMPRFIX , UNLINK => !$VERBOSE );
system("$GUNZIPEXE $reference_file | $FASTA21LINEXE > $ref_filename");

open($fh_ref,"<",$ref_filename) || die "# ERROR: cannot read $ref_filename\n";
my $ref_size = -s $ref_filename;

print "# flat size: $ref_size ($reference_file)\n";

my ($fh_test, $test_filename) = tempfile( $TEMPRFIX , UNLINK => !$VERBOSE );
system("$GUNZIPEXE $test_file | $FASTA21LINEXE > $test_filename");

open($fh_test,"<",$test_filename) || die "# ERROR: cannot read $test_filename\n";
my $test_size = -s $test_filename;

print "# flat size: $test_size ($test_file)\n\n";


## 4) read files using a buffer #################################################
# https://www.perlmonks.org/?node_id=457032

my $accum_pos = -1;
while( $raw_read = sysread($fh_raw , $raw_buffer , $BUFSIZE) ) {

	$ref_read = sysread($fh_ref , $ref_buffer , $BUFSIZE);

	$test_read = sysread($fh_test , $test_buffer , $BUFSIZE);

	if($raw_read != $test_read){
		die "# ERROR: sequences have different length\n";
	} 

	# take reference as gold standard and keep ading tp,tn,fp,fn
	my @raw = split(//,$raw_buffer);
	my @ref = split(//,$ref_buffer);
	my @test = split(//,$test_buffer);

	foreach $pos (0 .. $ref_read-1){

		# overall position after adding chromosomes
		$accum_pos++;

		$raw_bs = lc($raw[$pos]);

		# skip N-segments in raw sequence
      if($raw_bs eq 'n'){
         $N++;
         next;
      }

		# debugging
		#if($VERBOSE && $feat_set->contains( $accum_pos )) { print "$raw[$pos] $accum_pos\n" }

		$ord_ref  = ord($ref[$pos]);
		$ord_test = ord($test[$pos]);
		$ref_lc   = lc($ref[$pos]);
		$test_lc  = lc($test[$pos]);

		# skip mismatches
		if($test_lc ne 'n' && $ref_lc ne 'n' && $test_lc ne $ref_lc){ 
			$mism++;
			next; 
		} 
		elsif($VERBOSE){ print "$ref[$pos] $test[$pos]\n" } 

		# case of bases in buffer
   	# A 65 C 67 G 71 N 78 T 84
   	# a 97 c 99 g 103 n 110 t 116

		if( $feat_set->contains( $accum_pos )) { $feat_overlap = 1 }
		else{ $feat_overlap = 0 }

		if($ord_ref > 84){ # ref masked 

			$ref_mask++;

			if( $feat_overlap ){	$ref_feat++ }

			if($ord_test > 84){
				
				$tp++;						
				$test_mask++;
				
				if( $feat_overlap ){ $test_feat++ }

			} else { 
				$fn++;
			}
		} else { # ref non masked

			if($ord_test > 84){
            $fp++;
				$test_mask++;

				if( $feat_overlap ){ $test_feat++ }

         } else {
				$tn++;
			}
		}

		$total++;
	}
}

## 5) print stats summary ###########################################

print "#N\tmismatch\tmatch\tfeat\tref_overlap\tref_overlap_perc\ttest_overlap\ttest_overlap_perc\t".
	"ref_mask\tref_mask_perc\ttest_mask\ttest_mask_perc\t".
	"TP\tTN\tFP\tFN\tspecificity\tsensitivity\tFDR\tFOR\n";

printf("%d\t%d\t%d\t%s\t", $N, $mism, $total, $FEATTYPE );

printf("%d\t%1.2f\t", $ref_feat, 100*$ref_feat/$total);

printf("%d\t%1.2f\t", $test_feat, 100*$test_feat/$total);

printf("%d\t%1.2f\t", $ref_mask, 100*$ref_mask/$total);

printf("%d\t%1.2f\t", $test_mask, 100*$test_mask/$total);

printf("%d\t%d\t%d\t%d\t", $tp, $tn, $fp, $fn);

if( ($tn+$fp) > 0){
	printf("%1.3f\t",$tn/($tn+$fp));
} else { 
	printf("NA\t");
}

if( ($tp+$fn) > 0){
	printf("%1.3f\t",$tp/($tp+$fn));
} else { 
	print "NA\t" 
}

if( ($tp+$fp) > 0){
	printf("%1.3f\t", $fp / ($fp+$tp) );
} else {
	print "NA\t";
}

if( ($tn+$fn) > 0){
   printf("%1.3f\n", $fn / ($fn+$tn) );
} else {
   print "NA\n";
}



