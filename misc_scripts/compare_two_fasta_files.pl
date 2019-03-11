#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Std;
use File::Temp qw/ tempfile /;
use FindBin '$Bin';
use lib "$Bin/../Tools";
use FASTATools;

# Takes two FASTA nucleotide files and checks whether they contain the same headers & sequences
# I use to compare assemblies from ENA and other sources, such as JGI, to see whether
# GFF files map to the same sequences
# Bruno Contreras Moreira EMBL-EBI 2019 

# Example:
# perl ~/plant_tools/misc_scripts/compare_two_fasta_files.pl -i GCA_000263155.2.fasta.gz -j Sitalica_312_v2.fa.gz -a
#
# reading GCA_000263155.2.fasta.gz ...
# reading Sitalica_312_v2.fa.gz ...
# Aligning sequence fragments (size=10000) ...
# ENA|CM003528|CM003528.1	scaffold_1	identical	42145699	42145699	
# ENA|CM003529|CM003529.1	scaffold_2	different	49199997	49200776	
# ENA|CM003530|CM003530.1	scaffold_3	different	50651951	50652576	
# ENA|CM003531|CM003531.1	scaffold_4	different	40407879	40408058	
# ENA|CM003532|CM003532.1	scaffold_5	different	47252588	47253416	
# ENA|CM003533|CM003533.1	scaffold_6	different	36014550	36015257	
# ENA|CM003534|CM003534.1	scaffold_7	different	35964315	35964515	
# ENA|CM003535|CM003535.1	scaffold_8	different	40689132	40690061	
# ENA|CM003536|CM003536.1	scaffold_9	different	58970307	58970518

my $SKIPPEPTIDES = 1;
my $MINIDENTITY  = 98;
my $FRAGMENTSIZE = 10_000;
my $BLASTNEXE    = 'blastn'; # change as needed, assumes BLAST+

my (%opts,$INP_file1,$INP_file2,$INP_do_align);

getopts('hai:j:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0))
{
  print   "\nusage: $0 [options]\n\n";
  print   "-h this message\n";
  print   "-i FASTA filename         (required)\n";
  print   "-j FASTA filename         (required)\n";
  print   "-a align sequences        (optional, takes fragments of size=$FRAGMENTSIZE)\n\n";
  exit(0);
}

if(defined($opts{'i'})){  $INP_file1 = $opts{'i'}; }
else{ die "# EXIT : need -i FASTA filename\n"; }

if(defined($opts{'j'})){  $INP_file2 = $opts{'j'}; }
else{ die "# EXIT : need -j FASTA filename\n"; }

if(defined($opts{'a'})){  $INP_do_align = 1 }
else{ $INP_do_align = 0 }

##########################################################

my ($seq1,$seq2,$matched_seq,$header,$sequence,$report);
my (%id1toid2,%bits);

print "# reading $INP_file1 ...\n";
my $ref_FASTA1 = read_FASTA_file_array($INP_file1, $SKIPPEPTIDES);
print "# reading $INP_file2 ...\n";
my $ref_FASTA2 = read_FASTA_file_array($INP_file2, $SKIPPEPTIDES);

if($INP_do_align) {

	print "# Aligning sequence fragments (size=$FRAGMENTSIZE) ...\n";

	# create temporary file to store uncompressed file1
	my ($fh1,$tmpfile1) = tempfile();
	foreach $seq1 ( 0 .. $#{$ref_FASTA1} ) {
		print $fh1 ">$ref_FASTA1->[$seq1][NAME]\n".
			substr($ref_FASTA1->[$seq1][SEQ],0,$FRAGMENTSIZE)."\n";
	}

	# create temporary file to store uncompressed file2
	my ($fh2,$tmpfile2) = tempfile();
        foreach $seq2 ( 0 .. $#{$ref_FASTA2} ) {
                print $fh2 ">$ref_FASTA2->[$seq2][NAME]\n".
			substr($ref_FASTA2->[$seq2][SEQ],0,$FRAGMENTSIZE)."\n";
        }

	# run BLASTN
	my $cmd = "$BLASTNEXE -query $tmpfile1 -subject $tmpfile2 -outfmt 6 ".
                "-max_target_seqs 1 -perc_identity $MINIDENTITY -soft_masking true ";

	#ENA|CM003528|CM003528.1	scaffold_1	100.00	10000	0	0	1	10000	1	10000	0.0	18467
	#ENA|CM003529|CM003529.1	scaffold_2	100.00	5289	0	0	4712	10000	4712	10000	0.0	 9768
	#ENA|CM003529|CM003529.1	scaffold_2	100.00	3788	0	0	1	3788	1	3788	0.0	 6996
	open(BLASTN,"$cmd |") || die "# ERROR: cannot run $cmd\n";
	while(<BLASTN>) {
		chomp;
		my @data = split(/\t/,$_);
		next if($id1toid2{$data[0]});
		$id1toid2{$data[0]} = $data[1];
		$bits{$data[0]} = $data[11];
	}
	close(BLASTN);
}


## compare sequences in both files
print "#sequence1\tsequence2\tcomparison\tlength1\tlength2\n";
foreach $seq1 ( 0 .. $#{$ref_FASTA1} ) {

	$header = $ref_FASTA1->[$seq1][NAME];
	$matched_seq = undef;

	# check if corresponding sequence was previously determined by BLASTN alignment
	if(defined($id1toid2{$header})){

		foreach $seq2 ( 0 .. $#{$ref_FASTA2} ) {
                	if($ref_FASTA2->[$seq2][NAME] eq $id1toid2{$header}) {
                        	$matched_seq = $seq2;
                        	last;
                	}
		}
        } 
	
	# check for identical sequences if you feel lucky
	if(!defined($matched_seq)){

		foreach $seq2 ( 0 .. $#{$ref_FASTA2} ) {
			if($ref_FASTA1->[$seq1][SEQ] eq $ref_FASTA2->[$seq2][SEQ]){
				$matched_seq = $seq2;
                                last;		
			}
                }
	} else {

		# otherwise check whether there's a sequence with identical header in 2nd FASTA
		foreach $seq2 ( 0 .. $#{$ref_FASTA2} ) { 
			if($header eq $ref_FASTA2->[$seq2][NAME]) {
				$matched_seq = $seq2;
				last;
			}
		}
	}

	# print report
	$report = "$header\t";
	if(defined($matched_seq)){
	
		$report .= "$ref_FASTA2->[$matched_seq][NAME]\t";
	
		if($ref_FASTA1->[$seq1][SEQ] eq $ref_FASTA2->[$matched_seq][SEQ]){
			$report .= "identical\t";
                } else {
                	$report .= "different\t";        
                }

		$report .= length($ref_FASTA1->[$seq1][SEQ])."\t".length($ref_FASTA2->[$matched_seq][SEQ])."\t";
		
		# add BLASTN bitscore if available
		#if($bits{$header}){ $report .= "$bits{$header}" }
		#else{ $report .= "NA" }
	}
	else {
		$report .= "NA\t".length($ref_FASTA1->[$seq1][SEQ])."\tNA";
	}	

	print "$report\n";
}
