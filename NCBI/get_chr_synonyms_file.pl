#!/usr/bin/env perl
#
# takes ENA accessions and visits NCBI FTP to retrieve contig/chr synonyms from full assembly report
# 2018 Bruno Contreras Moreira EMBL-EBI
#
# Example calls: 
# $ ./get_chr_synonyms.pl -a GCA_000188115.3
# $ ./get_chr_synonyms.pl -a GCA_000188115.3 -R p2panreg -s solanum_lycopersicum
# $ ./get_chr_synonyms.pl -a Viridiplantae.2018-10-09.report.tsv -R p2panreg

use strict;
use warnings;
use Net::FTP;
use Getopt::Std;
use Bio::EnsEMBL::Registry;

# hard-coded paths, based on
# ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/188/115/GCA_000188115.3_SL3.0/GCA_000188115.3_SL3.0_assembly_report.txt
my $NCBIFTPURL = 'ftp.ncbi.nlm.nih.gov';
my $GGAGENPATH = '/genomes/all/GCA/'; # absolute

my (@ENA_accessions,@species_names,%opts);
my $input_is_file = 0;
my ($reg_file,$species_name) = ('','');
my ($acc2path,$fullpath,$species);
my ($report_file,$ENA_accession);

getopts('hs:R:a:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)){
  print "\nusage: $0 [options]\n\n";
  print " -h this message\n";
  print " -a ENA assembly accession OR file with 1 accession/line   (required, example: -a GCA_000188115.3)\n";
  print " -R add synonyms to db server pointed by to registry file  (optional, example: -R \$p2panreg)\n";
  print " -s species_name, required with -R & single ENA accession  (optional, example: -s solanum_lycopersicum)\n\n";
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
			$species =~ s/_+$//;	
			push(@species_names,$species);
		}
		close(INFILE);
		
		$input_is_file = 1;		

	} else { 
		push(@ENA_accessions,$opts{'a'}); # single accession
	}
}
else{ die "# ERROR: nead valid ENA accession\n\n# example: $0 -a GCA_000188115.3\n"; }

if($opts{'R'}){

	if(-e $opts{'R'}){ $reg_file = $opts{'R'} }
	else{ die "# EXIT : need a valid -R file, such as -R \$p2panreg\n" }

	if($input_is_file == 0){
		
		if($opts{'s'}){ 
			$species = $opts{'s'};
			push(@species_names,$species);
		}
		else{ die "# EXIT : need a valid -s species_db_name, such as -s solanum_lycopersicum\n" }
	}
}

printf("# %s -a %s -R %s -s %s\n\n", 
	$0, $opts{'a'}, $reg_file, $species_name );

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
	if($ENA_accession =~ /GCA_(\d+)\.\d+$/){

		$acc2path = $1;
                $acc2path =~ s/...\K(?=.)/\//sg;

		$species_name = $species_names[$input_sp];
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
	if(defined($reg_file)){

		my ($community_name,$seqrole,$ENA_acc,$insdc_acc,$refseq_acc);
		my ($insdc_db_id,$refseq_db_id,$seq_region_id);

		# connect to production db
		my $registry = 'Bio::EnsEMBL::Registry';
		$registry->load_all($reg_file);
		my $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species_names[$input_sp], "core");

		# find out ids of external dbs for synonyms
		my $external_db_sql = 
			"SELECT external_db_id FROM external_db WHERE db_name = ?;";

		my $sth = $dba->dbc->prepare($external_db_sql);
		$sth->execute('INSDC');
		my @results = $sth->fetchrow_array;		
		if(scalar(@results) > 0){ $insdc_db_id = $results[0] }
		else{ die "# ERROR: cannot find internal_db_id for INSDC\n" }

		$sth->execute('RefSeq');
                @results = $sth->fetchrow_array;             
                if(scalar(@results) > 0){ $refseq_db_id = $results[0] }  
                else{ die "# ERROR: cannot find internal_db_id for RefSeq\n" }

		# prepare seq_region_id query
		my $seq_region_id_sql =
                        "SELECT seq_region_id FROM seq_region WHERE name = ?;";
                my $sth2 = $dba->dbc->prepare($seq_region_id_sql);

		# prepare synonym sql
		my $synonym_sql =
			"INSERT IGNORE INTO seq_region_synonym ".
			"(seq_region_id, synonym, external_db_id) VALUES (?, ?, ?);";
		my $sth3 = $dba->dbc->prepare($synonym_sql);

		# actually parse and insert synonyms
		open(REPORT,'<',$report_file) || 
			die "# ERROR: cannot read $report_file\n";
		while(<REPORT>){
			# Sequence-Name Role Assigned-Molecule Assigned-Molecule-Location/Type GenBank-Accn Relationship RefSeq-Accn
			# SL3.0ch01	assembled-molecule	1	Chromosome	CM001064.3 <>	na
			# SL3.00SC0000001	unplaced-scaffold	na	na	AEKE03000001.1 <>	na
			next if(/^#/); 
			my @tsvdata = split(/\t/,$_);
			($community_name,$seqrole,$ENA_acc,$insdc_acc,$refseq_acc) = @tsvdata[0,1,2,4,6];
			$seq_region_id = '';

			next if($ENA_acc eq 'na');

			$sth2->execute($ENA_acc);

			my @results = $sth2->fetchrow_array;
			if(scalar(@results) > 0){ $seq_region_id = $results[0] }
                	else{ 
				print "# WARNING: cannot find seq_region_id for $ENA_acc, skip it\n";
				next;
			}

			# first INSDC synonyms
    			$sth3->execute($ENA_acc, $insdc_acc, $insdc_db_id);

			# now RefSeq synonyms		
			next if($refseq_acc eq 'na');

			$sth3->execute($ENA_acc, $refseq_acc, $refseq_db_id);	

		}
		close(REPORT);
	}
}

# 3) disconnect and exit
$ftp->quit();
