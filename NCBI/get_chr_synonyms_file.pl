#!/usr/bin/env perl
#
# takes ENA accession and visits NCBI FTP to retrieve contig/chr synonyms from full assembly report
# 2018 Bruno Contreras Moreira EMBL-EBI
#
# Example calls: 
# $ ./get_chr_synonyms.pl GCA_000188115.3
# $ ./get_chr_synonyms.pl Viridiplantae.2018-10-09.report.tsv

use strict;
use warnings;
use Net::FTP;
use Getopt::Std;

# hard-coded paths, based on
# ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/188/115/GCA_000188115.3_SL3.0/GCA_000188115.3_SL3.0_assembly_report.txt
my $NCBIFTPURL = 'ftp.ncbi.nlm.nih.gov';
my $GGAGENPATH = '/genomes/all/GCA/'; # absolute

# EXTERNAL_ID_DB values in db schema
my %ENSEXTIDDB = (
	'RefSeq' => 50840,
	'INSDC'  => 50710
);

my (@ENA_accessions,@species_names,%opts);
my $input_is_file = 0;
my ($prod_server,$ensembl_version,$species_db_name) = ('','','');
my ($acc2path,$fullpath,$species,$GCA_version,$eg_version,$prod_db_args);
my ($dbhost,$dbport,$dbuser,$dbpass,$report_file,$ENA_accession);

getopts('hs:p:a:v:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)){
  print "\nusage: $0 [options]\n\n";
  print " -h this message\n";
  print " -a ENA assembly accession OR file with 1 accession/line   (required, example: -a GCA_000188115.3)\n";
  print " -p add synonyms to this production db server              (optional, example: -p eg-p2-w)\n\n";
  print "When -a is a file produced by check_new_ENA_assemblies.pl:\n";
  print " -v next Ensembl version                                   (optional, example: -v 95, required with -p)\n\n";
  print "When -a is a single ENA accession:\n";
  print " -s species_name, required with -p                         (optional, example: -s solanum_lycopersicum_core_42_95_3)\n\n";
  exit(0);
}

## 0) check input params
if($opts{'a'}){ 

	if(-s $opts{'a'}){ # input file instead of accession code

		# ENA_assembly_id	scientific_name	...
		# GCA_000001735.2	Arabidopsis thaliana	...
		# GCA_000004515.4	Glycine max	...

		open(INFILE,'<',$opts{'a'}) || 
			die "# ERROR: cannot read $opts{'a'}\n";
		while(<INFILE>){
			next if(/^#/); # skip commented lines
			my @tsvdata = split(/\t/,$_);
			push(@ENA_accessions,$tsvdata[0]);
			$species = lc($tsvdata[1]);
			$species =~ s/\s+/_/;
			push(@species_names,$species);
		}
		close(INFILE);
		
		$input_is_file = 1;		

	} else { 
		push(@ENA_accessions,$opts{'a'}); # single accession
	}
}
else{ die "# ERROR: nead valid ENA accession\n\n# example: $0 -a GCA_000188115.3\n"; }

if($opts{'p'}){ 
	$prod_server = $opts{'p'};
	chomp( $prod_db_args = `$prod_server -details script` ); 
	if($prod_db_args =~ m/--host (\S+) --port (\d+) --user (\S+) --pass (\S+)/){
		($dbhost,$dbport,$dbuser,$dbpass) = ($1,$2,$3,$4);
	}
	else{ die "# ERROR: cannot parse db server details\n" }

	if($input_is_file == 1){
		if($opts{'v'}){
        		$ensembl_version = $opts{'v'};
			$eg_version = $ensembl_version-53;
		}
		else{ die "# EXIT : need a valid -v version, such as -v 95\n" }
	}
	else{
		if($opts{'s'}){ 
			$species_db_name = $opts{'s'};
			if($species_db_name =~ m/_core_\d+_\d+_\d+$/){
				push(@species_names,$opts{'s'});
			}
			else{ die "# EXIT : need a valid -s species_db_name, such as -s solanum_lycopersicum_core_42_95_3\n" }
		}
		else{ die "# EXIT : need a valid -s species_db_name, such as -s solanum_lycopersicum_core_42_95_3\n" }
	}
}

printf("# %s -a %s -p %s -v %s -s %s\n\n", 
	$0, $opts{'a'}, $prod_server, $ensembl_version, $species_db_name );

## 1) connect to FTP site 
my $ftp = Net::FTP->new($NCBIFTPURL, Debug => 0) ||
	die "# ERROR: cannot connect to $NCBIFTPURL: $@";

$ftp->login("anonymous",'-anonymous@') || 
	die "# ERROR: cannot login ", $ftp->message();

## 2) check all ENA accessions
foreach my $input_sp (0 .. $#ENA_accessions){

	$ENA_accession = $ENA_accessions[$input_sp];

	# 2.0) check accession format
	# https://stackoverflow.com/questions/25523899/in-perl-best-way-to-insert-a-char-every-n-chars
	if($ENA_accession =~ /GCA_(\d+)\.(\d+)$/){

		($acc2path,$GCA_version) = ($1,$2);
                $acc2path =~ s/...\K(?=.)/\//sg;

		$species_db_name = $species_names[$input_sp];

		if($input_is_file == 1){
                	$species_db_name .= '_core_'.$eg_version.'_'.$ensembl_version.'_'.$GCA_version;
        	}
        }
        else{ die "# EXIT : bad ENA accession: $ENA_accession lacks .version\n"; }

	# 2.1) compose folder path based on ENA accession
	$fullpath = $GGAGENPATH . $acc2path;

	# 2.2) get to appropriate remote folder
	$ftp->cwd($fullpath) || 
		die "# ERROR: cannot change working directory ", $ftp->message();

	# 2.3) check right assembly subfolder and try to copy assembly report file
	$report_file = 'NA';
	foreach my $subfolder ( $ftp->ls() ){
		if($subfolder =~ $ENA_accession) {
		
			$ftp->cwd($subfolder) || 
				die "# ERROR: cannot change working directory ", $ftp->message();

			foreach my $file ( $ftp->ls() ){

				if($file =~ $ENA_accession && $file =~ /assembly_report.txt/){

					$ftp->get($file) || 
						die "# ERROR: cannot download $file ", $ftp->message;	

					if(-s $file){ $report_file = $file }
					last;
				}
			}
			last;
		}
	} 

	print "# assembly report: $report_file\n\n";

	# 2.4) parse report file and add synonyms to db if requested 
	if(defined($dbpass)){

		my ($community_name,$seqtype,$insdc_acc,$refseq_acc);

		open(REPORT,'<',$report_file) || 
			die "# ERROR: cannot read $report_file\n";
		while(<REPORT>){
			# Sequence-Name Role Assigned-Molecule Assigned-Molecule-Location/Type GenBank-Accn Relationship RefSeq-Accn
			# SL3.0ch01	assembled-molecule	1	Chromosome	CM001064.3 <>	na
			# SL3.00SC0000001	unplaced-scaffold	na	na	AEKE03000001.1 <>	na
			next if(/^#/); 
			my @tsvdata = split(/\t/,$_);
			($community_name,$seqtype,$insdc_acc,$refseq_acc) = @tsvdata[0,3,4,6];

			print "$community_name,$seqtype,$insdc_acc,$refseq_acc  $report_file\n";
		}
		close(REPORT);

	}
}

# 3) disconnect and exit
$ftp->quit();
