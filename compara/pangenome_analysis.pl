#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);
use Net::FTP;
use JSON qw(decode_json);
use Data::Dumper;
use Benchmark;
use Time::HiRes;
use HTTP::Tiny;

# Simulates pangenome growth by counting clusters of orthologous genes shared by (plant) species in clade 
# by querying pre-computed Compara data from Ensembl Genomes
#
# Based on scripts at
# https://github.com/Ensembl/ensembl-compara/blob/release/97/scripts/examples/
# https://github.com/Ensembl/ensembl-rest/wiki/Example-Perl-Client
#
# Bruno Contreras Moreira 2019

# Ensembl Genomes 
my @divisions  = qw( Plants Bacteria Fungi Vertebrates Protists Metazoa );
my $FTPURL     = 'ftp.ensemblgenomes.org'; 
my $COMPARADIR = '/pub/xxx/current/tsv/ensembl-compara/homologies';
my $RESTURL    = 'http://rest.ensembl.org';
my $INFOPOINT  = $RESTURL.'/info/genomes/division/';
my $TAXOPOINT  = $RESTURL.'/info/genomes/taxonomy/';

my $verbose    = 0;
my $division   = 'Plants';
my $taxonid    = ''; # NCBI Taxonomy id, Brassicaceae=3700, Asterids=71274, Poaceae=4479
my $ref_genome = ''; # should be contained in $taxonid;
my $seqtype    = 'protein'; 
my $comparadir = '';
my $outfolder  = '';
my $out_genome = '';


my ($help,$sp,$show_supported,$request,$response);
my ($GOC,$WGA,$LOWCONF) = (0,0,0);
my ($request_time,$last_request_time) = (0,0);
my (@ignore_species, %ignore, %polyploid, %division_supported);

GetOptions(	
	"help|?"       => \$help,
	"verbose|v"    => \$verbose,
	"supported|l"  => \$show_supported,
	"division|d=s" => \$division, 
	"clade|c=s"    => \$taxonid,
	"reference|r=s"=> \$ref_genome,
	"outgroup|o=s" => \$out_genome,
	"ignore|i=s"   => \@ignore_species,
	"type|t=s"     => \$seqtype,
	"GOC|G=i"      => \$GOC,
	"WGA|W=i"      => \$WGA,
	"LC|L"         => \$LOWCONF,
	"folder|f=s"   => \$outfolder
) || help_message(); 

sub help_message {
	print "\nusage: $0 [options]\n\n".
      "-c NCBI Taxonomy clade of interest         (required, example: -c Brassicaceae or -c 3700)\n".
		"-f output folder                           (required, example: -f myfolder)\n".
		"-r reference species_name to name clusters (required, example: -r arabidopsis_thaliana)\n".
		"-l list supported species_names            (optional, example: -l)\n".
		"-d Ensembl division                        (optional, default: -d $division)\n".
		"-r reference species_name                  (optional, default: -r $ref_genome)\n".
		"-o outgroup species_name                   (optional, example: -o brachypodium_distachyon)\n".
		"-i ignore species_name(s)                  (optional, example: -i selaginella_moellendorffii -i ...)\n".
		"-t sequence type [protein|cdna]            (optional, requires -f, default: -t protein)\n".
		"-L allow low-confidence orthologues        (optional, by default these are skipped)\n".
		"-v verbose                                 (optional, example: -v\n";

	print "\nThe following options are only available for some clades:\n\n".
		"-G min Gene Order Conservation [0:100]  (optional, example: -G 75)\n".
		"   see modules/Bio/EnsEMBL/Compara/PipeConfig/EBI/Plants/ProteinTrees_conf.pm\n".
		"   at https://github.com/Ensembl/ensembl-compara\n\n".
		"-W min Whole Genome Align score [0:100] (optional, example: -W 75)\n".
		"   see ensembl-compara/scripts/pipeline/compara_plants.xml\n".
		"   at https://github.com/Ensembl/ensembl-compara\n\n";
	print "Read about GOC and WGA at:\n".
		"https://www.ensembl.org/info/genome/compara/Ortholog_qc_manual.html\n\n";

	print "Example calls:\n\n".
		" perl $0 -c Brassicaceae -f Brassicaceae\n".
		" perl $0 -c Brassicaceae -f Brassicaceae -t cdna -o theobroma_cacao\n".
		" perl $0 -f poaceae -c 4479 -r oryza_sativa -WGA 75\n".
		exit(0);
}

if($help){ help_message() }

if($taxonid eq ''){
	print "# ERROR: need a valid NCBI Taxonomy clade, such as -c Brassicaceae or -c 3700\n\n";
	print "# Check https://www.ncbi.nlm.nih.gov/taxonomy\n";
	exit;
}

if($ref_genome eq ''){
   print "# ERROR: need a valid reference species_name, such as -r arabidopsis_thaliana)\n\n";
   exit;
}

if($division){
	if(!grep(/$division/,@divisions)){
		die "# ERROR: accepted values for division are: ".join(',',@divisions)."\n"
	} else {
		$comparadir = $COMPARADIR;
		my $lcdiv = lc($division);
		$comparadir =~ s/xxx/$lcdiv/;
	}
}

if(@ignore_species){
	foreach my $sp (@ignore_species){
		$ignore{ $sp } = 1;
	}
	printf("\n# ignored species : %d\n\n",scalar(keys(%ignore)));
}

if($outfolder){
	if(-e $outfolder){ print "\n# WARNING : folder '$outfolder' exists, files might be overwritten\n\n" }
	else { 
		if(!mkdir($outfolder)){ die "# ERROR: cannot create $outfolder\n" }
	}

	if($seqtype ne 'protein' && $seqtype ne 'cdna'){
		die "# ERROR: accepted values for seqtype are: protein|cdna\n"
	}
} else {
	print "# ERROR: need a valid output folder, such as -f Brassicaceae\n\n";
   exit;
}

if($show_supported){ print "# $0 -l \n\n" }
else {
	print "# $0 -d $division -c $taxonid -r $ref_genome -o $out_genome ".
		"-f $outfolder -t $seqtype -G $GOC -W $WGA -L $LOWCONF\n\n";
}

my $start_time = new Benchmark();

# new object for REST requests
my $http = HTTP::Tiny->new();
my $global_headers = { 'Content-Type' => 'application/json' };
my $request_count = 0; # global counter to avoid overload

## 0) check supported species in division ##################################################

$request = $INFOPOINT."Ensembl$division?";

$response = perform_rest_action( $request, $global_headers );
my $infodump = decode_json($response);

foreach $sp (@{ $infodump }) {
	if($sp->{'has_peptide_compara'}){
		$division_supported{ $sp->{'name'} } = 1;
	}	
}

# list supported species and exit
if($show_supported){

	foreach $sp (sort(keys(%division_supported))){
		print "$sp\n";
	}
	exit;
}

# check outgroup is supported
if($out_genome && !$division_supported{ $out_genome }){
	die "# ERROR: genome $out_genome is not supported\n";
}

## 1) check species in clade ##################################################################

my $n_of_species = 0;
my (@supported_species, %supported, %core, %present);

$request = $TAXOPOINT."$taxonid?";

$response = perform_rest_action( $request, $global_headers );
$infodump = decode_json($response);

foreach $sp (@{ $infodump }) {
   if($sp->{'name'} && $division_supported{ $sp->{'name'} }){

		next if($ignore{ $sp->{'name'} });

		# add sorted clade species except reference
		$supported{ $sp->{'name'} } = 1;
		if( $sp->{'name'} ne $ref_genome ){
			push(@supported_species, $sp->{'name'});
		}
   }
}

# check reference genome is supported 
if(!$supported{ $ref_genome }){
	die "# ERROR: cannot find $ref_genome within NCBI taxon $taxonid\n";
} else {
	unshift(@supported_species, $ref_genome);

	if($verbose){
		foreach $sp (@supported_species){
			print "# $sp\n";
		}
	}
}

printf("# supported species in NCBI taxon %s : %d\n\n", $taxonid, scalar(@supported_species));

# add outgroup if required
if($out_genome){
	push(@supported_species,$out_genome);
	$supported{ $out_genome } = 1;
	print "# outgenome: $out_genome\n";
}

$n_of_species = scalar( @supported_species );
print "# total selected species : $n_of_species\n\n";

## 2) get orthologous (plant) genes shared by selected species ####################

# columns of TSV file 
my ($gene_stable_id,$prot_stable_id,$species,$identity,$homology_type,$hom_gene_stable_id,
   $hom_prot_stable_id,$hom_species,$hom_identity,$dn,$ds,$goc_score,$wga_coverage,
	$high_confidence,$homology_id);

# iteratively get and parse TSV files, starting with reference
# NOTE: these files are bulky and might take some time to download
my (@sorted_ids);

foreach $sp ( @supported_species ){
	
	my $stored_compara_file = download_TSV_file( $comparadir, $sp );
	open(TSV,"gzip -dc $stored_compara_file |") || die "# ERROR: cannot open $stored_compara_file\n";
	while(<TSV>){
	
		($gene_stable_id,$prot_stable_id,$species,$identity,$homology_type,$hom_gene_stable_id,
		$hom_prot_stable_id,$hom_species,$hom_identity,$dn,$ds,$goc_score,$wga_coverage,    
		$high_confidence,$homology_id) = split(/\t/);

		next if(!$supported{ $hom_species });

		next if($LOWCONF == 0 && ($high_confidence eq 'NULL' || $high_confidence == 0));

		next if($WGA && ($wga_coverage eq 'NULL' || $wga_coverage < $WGA));

		next if($GOC && ($goc_score eq 'NULL' || $goc_score < $GOC));

		if($homology_type =~ m/ortholog/) {

			# add $ref_genome protein 
			if(!$core{ $gene_stable_id }){ 

				push(@{ $core{ $gene_stable_id }{ $ref_genome } }, $prot_stable_id );

				$present{ $ref_genome }++;

				push(@sorted_ids, $gene_stable_id); # save cluster order
			}

			push(@{ $core{ $gene_stable_id }{ $hom_species } }, $hom_prot_stable_id );

			$present{ $hom_species }++;

		} elsif($homology_type =~ m/within_species_paralog/) {
		
		}
	}
	close(TSV);
}

exit;

# check GOC / WGA availability
foreach $hom_species (@supported_species){
	if(!defined( $present{ $hom_species } )&& $WGA){
		print "# WGA not available: $hom_species\n";
	} elsif(!defined( $present{ $hom_species } ) && $GOC){
		print "# GOC not available: $hom_species\n";
	}
}


## 3) print summary matrix of single-copy / core genes and compile sequence clusters #################

my $total_core_clusters = 0;
my ($pruned_species,$treedump,$acc,$seq,$line,$filename);

# prepare param to prune species in REST requests
if($outfolder){
	foreach $hom_species (@supported_species){
		$pruned_species .= "prune_species=$hom_species;";
	}
}	

foreach $gene_stable_id (@sorted_ids){

	next if(scalar(keys(%{ $core{ $gene_stable_id } })) < $n_of_species); 

	my (%valid_prots, %align);

	$filename = $gene_stable_id;
	if($outfolder){
		if($seqtype eq 'protein'){ $filename .= '.faa' }
		else{ $filename .= '.fna' }
	}

	# print matrix header
	if($total_core_clusters == 0){
		print "cluster";
		foreach $hom_species (@supported_species){
			print "\t$hom_species";
		} 
		print "\n";
	}
	
	# print a matrix row in TSV format
	print $filename;
	foreach $hom_species (@supported_species){ 

		printf("\t%s", join(',',@{ $core{ $gene_stable_id }{ $hom_species } }) );

		# store which prots come from each species
		foreach $hom_prot_stable_id (@{ $core{ $gene_stable_id }{ $hom_species } }){
			$valid_prots{ $hom_prot_stable_id } = $hom_species;
		}
	} print "\n";

	$total_core_clusters++; #last if($total_core_clusters == 1); # debug
}

print "\n# total single-copy core clusters : $total_core_clusters\n\n";

# print diagnostics
if($total_core_clusters == 0){

	print "# diagnostic stats (species\tclusters)\n\n";

	foreach $hom_species (@supported_species){ $present{ $hom_species } = 0 unless($present{ $hom_species }) }

	foreach $hom_species (sort {$present{$a}<=>$present{$b}} (@supported_species)){
		printf("%s %d\n",	$hom_species, $present{ $hom_species } );
	}
}

my $end_time = new Benchmark();
print "\n# runtime: ".timestr(timediff($end_time,$start_time),'all')."\n";

###################################################################################################

# download compressed TSV file from FTP site, renames it 
# and saves it in current folder; uses FTP globals defined above
sub download_TSV_file {

	my ($dir,$ref_genome) = @_;
	my ($compara_file,$stored_compara_file) = ('','');

	#print "# connecting to $FTPURL ...\n";

	if(my $ftp = Net::FTP->new($FTPURL,Passive=>1,Debug =>0,Timeout=>60)){
		$ftp->login("anonymous",'-anonymous@') ||
			die "# cannot login ". $ftp->message();
		$ftp->cwd($dir) ||
		   die "# ERROR: cannot change working directory to $dir ". $ftp->message();
		$ftp->cwd($ref_genome) ||
			die "# ERROR: cannot find $ref_genome in $dir ". $ftp->message();

		# find out which file is to be downloaded and 
		# work out its final name with $ref_genome in it
		foreach my $file ( $ftp->ls() ){
			if($file =~ m/protein_default.homologies.tsv.gz/){
				$compara_file = $file;
				$stored_compara_file = $compara_file;
				$stored_compara_file =~ s/tsv.gz/$ref_genome.tsv.gz/;
				last;
			}
		}
		
		# download that TSV file
		unless(-s $stored_compara_file){
			$ftp->binary();
			my $downsize = $ftp->size($compara_file);
			$ftp->hash(\*STDOUT,$downsize/20) if($downsize);
			printf("# downloading %s (%1.1fMb) ...\n",$stored_compara_file,$downsize/(1024*1024));
			print "# [        50%       ]\n# ";
			if(!$ftp->get($compara_file)){
				die "# ERROR: failed downloading $compara_file\n";
			}

			# rename file to final name
			rename($compara_file, $stored_compara_file);
			print "# using $stored_compara_file\n\n";
		} else {
			print "# re-using $stored_compara_file\n\n";
		}
	} else { die "# ERROR: cannot connect to $FTPURL , please try later\n" }

	return $stored_compara_file;
}


# uses global $request_count
sub perform_rest_action {
	my ($url, $headers) = @_;
	$headers ||= {};
	$headers->{'Content-Type'} = 'application/json' unless exists $headers->{'Content-Type'};

	if($request_count == 15) { # check every 15
		my $current_time = Time::HiRes::time();
		my $diff = $current_time - $last_request_time;

		# if less than a second then sleep for the remainder of the second
		if($diff < 1) {
			Time::HiRes::sleep(1-$diff);
		}
		# reset
		$last_request_time = Time::HiRes::time();
		$request_count = 0;
	}

	my $response = $http->get($url, {headers => $headers});
	my $status = $response->{status};
	
	if(!$response->{success}) {
		# check for rate limit exceeded & Retry-After (lowercase due to our client)
		if(($status == 429 || $status == 599) && exists $response->{headers}->{'retry-after'}) {
			my $retry = $response->{headers}->{'retry-after'};
			Time::HiRes::sleep($retry);
			# afterr sleeping see that we re-request
			return perform_rest_action($url, $headers);
		}
		else {
			my ($status, $reason) = ($response->{status}, $response->{reason});
			die "# ERROR: failed REST request $url\n# Status code: ${status}\n# Reason: ${reason}\n# Please re-run";
		}
	}

	$request_count++;

	if(length($response->{content})) { return $response->{content} } 
	else { return '' }	
}

