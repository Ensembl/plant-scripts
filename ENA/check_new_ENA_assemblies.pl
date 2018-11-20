#!/usr/bin/env perl
# checks ENA annotated assemblies not currently included in a give Ensembl division, 
# EnsemblPlants by default 
# 2018 Bruno Contreras Moreira EMBL-EBI
#
# Example run, takes a few hours: 
# $ ./check_new_ENA_assemblies.pl -t Viridiplantae 
#
# Or if you want to re-run without re-downloading ENA files: 
# $ ./check_new_ENA_assemblies.pl -t Viridiplantae -S

# TO BE DONE: check whether GFF file is available from NCBI/RefSeq

use strict;
use warnings;
use POSIX qw(strftime);
use Getopt::Std;
use HTTP::Tiny;
use JSON;
use Net::FTP;

# hard-coded default parameters
my $required_assembly_state = 'chromosome';
my $required_genome_representation = 'full';
my $ensembl_division = 'EnsemblPlants';
my $ENAfolder = './ENAfiles/'; 
my $used_saved_files = 0;
my $KEEPCHR = 1; # bool, set to zero to avoid keeping local copies of ENA files

# URLs and query paths, can change with time
my $RESTURL    = 'http://rest.ensemblgenomes.org';
my $ENAURL     = 'https://www.ebi.ac.uk/ena/portal/api/search?result=assembly&query=tax_tree';
my $ENAquery   = '&fields=assembly_level%2Cgenome_representation%2Cassembly_name'.
			'%2Cassembly_title%2Cscientific_name%2Cstrain%2Cstudy_accession';
my $ENAviewURL = 'https://www.ebi.ac.uk/ena/data/view';
# also NCBI
my $NCBIFTPURL = 'ftp.ncbi.nlm.nih.gov';
my $GGAGENPATH = '/genomes/all/GCA/'; # absolute

# binaries, edit if required
my $WGETEXE = 'wget';


## 0) handle user arguments
my (%opts, $taxon_name);
getopts('hSd:t:r:s:f:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)){
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-t taxon scientific name as in NCBI Taxonomy   (required, example: -t Arabidopsis%20thaliana)\n";
  print "-d Ensembl division                            (optional, default: -d $ensembl_division)\n";
  print "-f path to folder to save downloaded ENA files (optional, default: -f $ENAfolder)\n";
  print "-S used available local ENA files              (optional, default: download from scratch)\n";
  print "-s required ENA assembly state                 (optional, default: -s $required_assembly_state)\n";
  print "-r required ENA assembly representation        (optional, default: -r $required_genome_representation)\n\n";
  print "Please check https://www.ebi.ac.uk/ena/browse/assembly-format for valid values for -r/-s\n\n";
  exit(0);
}

if(defined($opts{'t'})){  
  $taxon_name = $opts{'t'}; 
}
else{ die "# EXIT : need a valid -t taxon_name, such as Viridiplantae or Arabidopsis%20thaliana\n"; }

if($opts{'d'}){ $ensembl_division = $opts{'d'} }

if($opts{'s'}){ $required_assembly_state = $opts{'s'} }

if($opts{'r'}){ $required_genome_representation = $opts{'r'} }

if($opts{'f'}){ $ENAfolder = $opts{'f'} }
if(!-d $ENAfolder){ mkdir($ENAfolder) }

if($opts{'S'}){ $used_saved_files = 1 }

# complete REST URLs based on arguments
my $RESTtaxonomy_query = "$RESTURL/taxonomy/name/$taxon_name?";
my $RESTdivision_query = "$RESTURL/info/genomes/division/$ensembl_division?";

printf("# %s -t %s -d %s -s %s -r %s -f %s -S %d (\$KEEPCHR=%d)\n\n",
  $0, $taxon_name, $ensembl_division, $required_assembly_state, 
  $required_genome_representation, $ENAfolder, $used_saved_files, $KEEPCHR);


## 1) get taxon ID from REST ensemblgenomes
my $taxon_id = '';
open(REST,"$WGETEXE -q --header='Content-type:application/json' $RESTtaxonomy_query -O- |") || 
	die "ERROR: cannot connect to $RESTtaxonomy_query\n";
while(<REST>) {
	# ids are listed from ancestor (Eukarya) to actual taxon of interest
	while(/"id":"(\d+)"/g) {
		$taxon_id = $1
	}
}
close(REST);
print "# taxon '$taxon_name' corresponds to taxonomy_id=$taxon_id\n\n"; 

# compose output filenames
my $datestring = strftime("%F", localtime());
my $ENAcurrent_file = 'ENA.'.$taxon_name.'.'.$datestring.'.tsv';
my $ENAlogfile = 'ENA.'.$taxon_name.'.'.$datestring.'.err';
my $report_file = $taxon_name.'.'.$datestring.'.report.tsv';

## 2) save local copy of ENAs current $taxon_name contents
system("$WGETEXE -q '$ENAURL($taxon_id)$ENAquery' -O $ENAcurrent_file -o $ENAlogfile");

## 3) check currently supported assemblies in Ensembl division  
my %current_ensembl_assembly_ids;
my $http = HTTP::Tiny->new();
my $response = $http->get($RESTdivision_query, {
	headers => { 'Content-type' => 'application/json' }
}); die "Failed!\n" unless $response->{success};
          
if(length($response->{content})) {
	my $array_ref = decode_json($response->{content});
	foreach my $sp_hash (@$array_ref) {
		if(defined($sp_hash->{assembly_id})){
			$current_ensembl_assembly_ids{ $sp_hash->{assembly_id} } = 1;
			#print "$sp_hash->{assembly_id}\n"; # debug
		}
	}
}
     
## 4) check ENA annotated assemblies which might not be currently supported in Ensembl division
my ($assembly_id,$assembly_state,$genome_representation,$strain,$fullpath,$acc2path);
my ($full_assembly_id,$chr_acc,$n_of_genes,$scientific_name,$ENAfile);

# connect to NCBI FTP site 
my $ftp = Net::FTP->new($NCBIFTPURL, Debug => 0) ||
	die "# ERROR: cannot connect to $NCBIFTPURL: $@";

$ftp->login("anonymous",'-anonymous@') || 
	die "# ERROR: cannot login ", $ftp->message();


open(REPORT,">",$report_file) || die "# ERROR: cannot create $report_file\n";
print REPORT "# ENA_assembly_id\tscientific_name\tstrain\tannotated_genes\tchromosomes\tchromosome_accessions\n";

open(LOCALENA,$ENAcurrent_file) || die "ERROR: cannot read $ENAcurrent_file\n";
while(<LOCALENA>){
	
	#accession	assembly_level	genome_representation	assembly_name	assembly_title	scientific_name	strain	study_accession
	#GCA_000001735	chromosome	full	TAIR10;TAIR10.1	TAIR10.1 assembly for Arabidopsis thaliana	Arabidopsis thaliana
	
	my @data = split(/\t/,$_);
	($assembly_id,$assembly_state,$genome_representation,$scientific_name,$strain) = @data[0,1,2,5,6];
	next if($assembly_id eq 'accession');

	if($strain eq ''){ $strain = 'NA' }

	# get this assembly most recent version attached to GCA accession + linked chr accessions
	# NOTE: when you request a view you get the most recent version (Dan Bolser)
	my @chr_accessions;

	open(XMLENA,"$WGETEXE -q '$ENAviewURL/$assembly_id&display=xml' -O- -o $ENAlogfile |") ||
		die "ERROR: cannot query $ENAviewURL/$assembly_id&display=xml\n";
	while(<XMLENA>){

		# <ASSEMBLY accession="GCA_000001735.2" alias="TAIR10.1" ... ">
		if(/<ASSEMBLY accession=\"(\S+)\"/){
			$full_assembly_id = $1;
		}
		elsif(/<CHROMOSOME accession=\"(\S+)\"/){ #<CHROMOSOME accession="CP002684.1">
			$chr_acc = $1;
			push(@chr_accessions,$chr_acc);
		}
	}
	close(XMLENA);	

	# skip assemblies non reaching the required state/genome representation or 
	# those with matching version in Ensembl division
	next if($assembly_state ne $required_assembly_state || 
		$genome_representation ne $required_genome_representation ||
		defined($current_ensembl_assembly_ids{$full_assembly_id}) );


	# check whether NCBI/RefSeq provide GFF for this assembly
	if($full_assembly_id =~ /GCA_(\d+)\.\d+$/){

		$acc2path = $1;
                $acc2path =~ s/...\K(?=.)/\//sg;
        }
        else{ die "# EXIT : bad ENA accession: $full_assembly_id lacks .version\n"; }

	# compose folder path based on ENA full assembly id
	$fullpath = $GGAGENPATH . $acc2path;

	# 2.2) get to appropriate remote folder
	$ftp->cwd($fullpath) || 
		die "# ERROR: cannot change working directory ", $ftp->message();



	# compute the number of annotated genes on top of this assembly 
	$n_of_genes = 0;
	foreach $chr_acc (@chr_accessions){

		# download a local copy of all chr files of this assembly
		# Note: XML is smaller than text in my Athaliana tests, 
		# but EMBL/text format is more readable
		$ENAfile = "$ENAfolder/$chr_acc.embl.gz";

		if($used_saved_files && -s $ENAfile){
			print "# re-using $ENAfile\n";
		}
		else {
			print "# downloading $ENAfile\n";
        		system("$WGETEXE -q -c -N '$ENAviewURL/$chr_acc&display=text&download=gzip' -O $ENAfile -o $ENAlogfile");
			if($? != 0){
                		die "ERROR: cannot download $ENAviewURL/$chr_acc&display=text&download=gzip\n";
			}
		}	

		open(EMBLENA,"zcat $ENAfile |") || die "ERROR: cannot read $ENAfile\n";
        	while(<EMBLENA>){
			if(/^FT   gene/){
				$n_of_genes++;
			}
		}
		close(EMBLENA);

		unlink($ENAfile) if($KEEPCHR == 0);
	}

	print REPORT "$full_assembly_id\t$scientific_name\t$strain\t$n_of_genes\t".
		scalar(@chr_accessions)."\t".join(',',@chr_accessions)."\n";
}
close(LOCALENA);

close(REPORT);

$ftp->quit();

print "# created report $report_file\n";
