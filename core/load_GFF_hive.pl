#!/bin/env perl
use strict;
use warnings;
use Getopt::Std;
use POSIX qw(strftime);
use Cwd;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Slice qw(split_Slices);
use Bio::Seq;

# TO BE DONE: automatically extract a stable_gene_id\tdescription file to load descriptions


# This script takes a GFF3 & a peptide FASTA file and attempts to load the 
# features on top of a previously loaded ENA genome assembly in hive.
# This should be run after loading a genome from ENA 
# NOTE: it will update meta.genebuild dates and version
# NOTE: hive pipelines must be run in eb-cli nodes
#
# It uses env $USER to create hive job names and assumes Ensembl-version API
# is loaded in @INC / $PERL5LIB
#
# Adapted from Dan Bolser's run_the_gff_loader2.sh by B Contreras Moreira 2018-9
#
# https://www.ebi.ac.uk/seqdb/confluence/display/EnsGen/Load+GFF3+Pipeline
#
## check user arguments ######################################################
##############################################################################

my (%opts,$species,$protein_fasta_file,$gff3_file,$gene_source,$ensembl_version);
my ($prodbname,$pipeline_dir,$reg_file,$hive_args,$hive_db,$hive_url,$argsline);
my ($skip_bad_types,$rerun,$sub_chr_names,$nonzero,$synonyms,$overwrite,$max_feats) = (0,0,'',0,0,0,0);
my ($check_gff_CDS,$check_chr_ends,$add_to_previous,$names_stable,$otherfeats) = (0,0,0,0,0);
my ($new_gff3file,$short_gff3file,$synonym_file,$gbkfile,%synonym) = ('','','','');
my $logicname = '';
my @seqregion_types = ('chromosome','supercontig','contig');
my $hive_db_cmd = 'mysql-ens-hive-prod-2-ensrw';

getopts('hTNoawzyrceY:n:s:f:g:L:S:v:R:H:P:m:D:G:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)){
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-s species_name                                (required, example: -s arabidopsis_thaliana)\n";
  print "-f protein FASTA file                          (required, example: -f atha.pep.fasta)\n";
  print "-v next Ensembl version                        (required, example: -v 95)\n";
  print "-g GFF3 file                                   (required, example: -g atha.gff)\n";
  print "-R registry file, can be env variable          (required, example: -R \$p1panreg)\n";
  print "-D ensembl_production db name                  (required, example: -D ensembl_production_97)\n";
  print "-P folder to put pipeline files, can be env    (required, example: -P \$gfftmp)\n";
  print "-H hive database command                       (optional, default: $hive_db_cmd)\n";
  print "-G GenBank file to load sequence edits         (optional, example: -G atha.gbk)\n";
  print "-m max genes to load                           (optional, default: all loaded)\n";
  print "-a incrementally add GFF to prev annotation    (optional, default: delete previous)\n";
  print "-N take name instead of ID as stable_id        (optional, default: ID is stable_id)\n";
  print "-n replace chr names with Perl-regex           (optional, example: -n 'SL3.0ch0*(\\d+)' )\n";
  print "-z skip chr zero in GFF3 file                  (optional, requires -n)\n";
  print "-y saves original chr names as synonyms in db  (optional, requires -n)\n";
  print "-e print sequence of chromosome ends in db     (optional, requires -n)\n";
  print "-Y add contig synonyms from file before load   (optional, takes a TSV file: gffname\\ENAtsynonym)\n";
  print "-c check CDS coords in GFF3                    (optional, requires -m)\n";
  print "-T skip non-recognized feature types           (optional)\n";
  print "-S source of annotation for dumps, one word    (optional, example: -S SOL, default: 3rd col of GFF3)\n";
  print "-L logicname to import analysis desc & webdata (optional, example: -L SOL_genomics, should be in production db)\n";
  print "-o add features to other_features db           (optional, default: add to core db)\n";
  print "-w over-write db (hive_force_init)             (optional, useful when a previous run failed)\n";                             
  print "-r re-run jump to beekeper.pl                  (optional, default: run init script from scratch)\n\n";
  exit(0);
}

if($opts{'v'}){
        $ensembl_version = $opts{'v'};

        # check Ensembl API is in env
	if(!grep(/ensembl-$ensembl_version\/ensembl-hive\/modules/,@INC)){
                die "# EXIT : cannot find ensembl-$ensembl_version/ensembl-hive/modules in \$PERL5LIB / \@INC\n"
        } 
}
else{ die "# EXIT : need a valid -v version, such as -v 95\n" }

if($opts{'s'}){ 
	$species = $opts{'s'}; 
	$hive_db = $ENV{'USER'}."_load_gff3_$species";  
} 
else{ die "# EXIT : need a valid -s species_name, such as -s arabidopsis_thaliana\n" }

if($opts{'f'} && -e $opts{'f'}){ $protein_fasta_file = $opts{'f'} }
else{ die "# EXIT : need a valid -f file, such as -f atha.pep.fasta\n" }

if($opts{'g'} && -e $opts{'g'}){ $gff3_file = $opts{'g'} }
else{ die "# EXIT : need a valid -g file, such as -g atha.gff\n" }

if($opts{'S'}){ $gene_source = $opts{'S'} }
else{
	chomp( $gene_source = `grep -v "^#" $gff3_file | grep gene | cut -f 2 | head -1` );

	if(!$gene_source){ 
		die "# EXIT : cannot parse annotation source from $gff3_file\n";
	} elsif($gene_source =~ m/\s+/) {
		die "# EXIT : annotation source '$gene_source' parsed from $gff3_file contains blanks\n";
	}
}

if($opts{'L'}){ $logicname = $opts{'L'} }

if($opts{'o'}){ $otherfeats = 1 }

if($opts{'R'} && -e $opts{'R'}){ $reg_file = $opts{'R'} }
else{ die "# EXIT : need a valid -R file, such as -R \$p2panreg\n" }

if($opts{'D'}){ $prodbname = $opts{'D'} }
else{ die "# EXIT : need a valid -D value, such as -D ensembl_production_97\n" }

if($opts{'G'}){ 
	$gbkfile = $opts{'G'};
	if(!-e $gbkfile){ 
		die "# EXIT : need a valid -G file\n" 
	}
}

if($opts{'H'}){ $hive_db_cmd = $opts{'H'} }
chomp( $hive_args = `$hive_db_cmd details script` );
chomp( $hive_url  = `$hive_db_cmd --details url` );
$hive_url .= $hive_db;

if($opts{'P'} && -d $opts{'P'}){ $pipeline_dir = "$opts{'P'}/$species\_$ensembl_version" }
else{ die "# EXIT : need a valid -P folder to put pipeline files, such as -P \$gfftmp\n" }

if($opts{'m'} && $opts{'m'} > 0){ 
	$max_feats = int($opts{'m'}); 
	if($opts{'c'}){ $check_gff_CDS = 1 }
}

if($opts{'T'}){ $skip_bad_types = 1 }

if($opts{'r'}){ $rerun = 1 }

if($opts{'w'}){ $overwrite = 1 }

if($opts{'a'}){ $add_to_previous = 1 }

if($opts{'N'}){ $names_stable = 1 }

if($opts{'n'}){ 
	$sub_chr_names = $opts{'n'};

	if($opts{'z'}){ $nonzero = 1 }
	if($opts{'y'}){ $synonyms = 1 }
	if($opts{'e'}){ $check_chr_ends = 1 }
}
elsif($opts{'Y'}){ $synonym_file = $opts{'Y'} }

$argsline = sprintf("%s -s %s -f %s -g %s -S %s -L %s -o %d -v %s -R %s -D %s -H %s -P %s -G %s -m %d ".
		"-Y %s -n '%s' -z %d -y %d -e %d -c %d -N %d -a %d -T %d -w %d -r %d",
  $0, $species, $protein_fasta_file, $gff3_file, $gene_source, $logicname, $otherfeats,
  $ensembl_version, $reg_file, $prodbname, $hive_db_cmd, $pipeline_dir, $gbkfile, $max_feats,
  $synonym_file, $sub_chr_names, $nonzero, $synonyms, 
  $check_chr_ends, $check_gff_CDS, $names_stable, 
  $add_to_previous, $skip_bad_types, $overwrite, $rerun );

print "# $argsline\n\n";


## check ID=names in GFF3 file to warn about gene:, mRNA:,... tags, which otherwise are added as stable_ids to db
#SL3.0ch00	maker_ITAG	gene	16480	17940	.	+	.	ID=gene:Solyc00g005000.3...
#SL3.0ch00	maker_ITAG	mRNA	16480	17940	.	+	.	ID=mRNA:Solyc00g005000.3.1...
#SL3.0ch00	maker_ITAG	exon	16480	16794	.	+	.	ID=exon:Solyc00g005000.3.1.1...
#SL3.0ch00	maker_ITAG	CDS	16480	16794	.	+	0	ID=CDS:Solyc00g005000.3.1.1...
open(GFF,'<',$gff3_file) || die "# ERROR: cannot read $gff3_file\n";
while(<GFF>){
	next if(/^#/);
	my @gffdata = split(/\t/,$_);

	if($gffdata[8] && $gffdata[8] =~ m/ID=gene:/){ 

		print "# ERROR: please edit the GFF file to remove redundant ID names:\n$_\n\n";
		#print "# You can try: \$ perl -lne 's/ID=\\w+:/ID=/; print' <gff3file> \n\n";
		exit(0);
	}
}
close(GFF);

## replace chr names with natural & add original names to synonyms in db;
## also check whether coords in GFF3 are within chromosomes in db
########################################################################
if($sub_chr_names ne ''){

	# connect to production db
	my $registry = 'Bio::EnsEMBL::Registry';
	$registry->load_all($reg_file);
	my $slice_adaptor = $registry->get_adaptor($species, "core", "slice");

	my $prev_sep_line = 0; # to count ### separators
	my (%synonyms,$chr_int,$chr_orig,$chr_start,$chr_end,$chr_length);
	$new_gff3file = $gff3_file . '.edited';

	open(NEWGFF,'>',$new_gff3file) || die "# ERROR: cannot create $new_gff3file\n";

	# read original GFF3 file, chr name is 1st column
	open(GFF,'<',$gff3_file) || die "# ERROR: cannot read $gff3_file\n";
	while(<GFF>){
		if(/^$sub_chr_names/){ # /^SL3.0ch0*(\d+)/

			$chr_int = $1; # natural chr number
			next if($nonzero == 1 && $chr_int eq '0');
			
			my @gffdata = split(/\t/,$_); 
			$chr_orig = shift(@gffdata); # original chr name, 1st column

			if(!defined($synonyms{ $chr_int })){
				print "# chr $chr_orig replaced by $chr_int\n";
			}
			$synonyms{ $chr_int } = $chr_orig;
				
			print NEWGFF "$chr_int\t" . join("\t",@gffdata);
		}
		elsif(/^##sequence-region\s+$sub_chr_names\s+(\d+)\s+(\d+)/){ ##sequence-region   SL3.0ch00 16480 20797619
			$chr_int = $1; # natural chr number
                        
			next if($nonzero == 1 && $chr_int eq '0');

			# check chr end sequences
			if($check_chr_ends){
				($chr_start,$chr_end) = ($2,$3);
				my $chr_slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr_int );
				$chr_length = length($chr_slice->seq);

				print "# chr $chr_int length=$chr_length ";
				if($chr_end > $chr_length){ print "WARNING: $chr_end > $chr_length" }
				print "\n";
				print "# ".substr($chr_slice->seq,0,30) ." .. ".substr($chr_slice->seq,-30)."\n"; 
			}
			
			print NEWGFF $_;
		} 
		else{ 
			next if($prev_sep_line == 1 && /^###/);
			print NEWGFF $_; 
			
			if(/^###/){ $prev_sep_line = 1 }
			else{ $prev_sep_line = 0 }
		} 
	}
	close(GFF);

	close(NEWGFF);

	print "# created edited GFF3 file: $new_gff3file\n\n";

	# add chr name synonyms to target db
	if($synonyms == 1){
			
		for $chr_int (keys(%synonyms)){

			$chr_orig = $synonyms{ $chr_int };
			my $chr_slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr_int );

			print "# adding synonym $chr_orig ($chr_int)\n";
			$chr_slice->add_synonym( $chr_orig );
		}

		# DanBolser's SQL queries for doing just that
		# find out coord_system_id for complete chromosomes
		#SELECT coord_system_id FROM solanum_lycopersicum_core_42_95_3.coord_system WHERE name = "chromosome";

		#INSERT INTO seq_region_synonym (seq_region_id, synonym, external_db_id)
		#SELECT seq_region_id, CONCAT("SL3.0ch", name), 50691 FROM seq_region
		#WHERE coord_system_id = 2 AND name > 9;

		#find out external_db_id of original chr ids
		#SELECT * FROM external_db WHERE db_name rlike 'sgn';
	}
}
else{ $new_gff3file = $gff3_file }

if($synonym_file ne '' && -s $synonym_file){	

	print "# adding synonyms from $synonym_file\n";

	# connect to production db
	my ($orig_gff,$synonym,$contig_slice,$regtype);
	my $registry = 'Bio::EnsEMBL::Registry';
	$registry->load_all($reg_file); 
	my $slice_adaptor = $registry->get_adaptor($species, "core", "slice");

	# parse synonyms
	my $n_of_synonyms = 0;
	open(TSV,"<",$synonym_file) || 
		die "# ERROR: cannot read $synonym_file\n";
	while(<TSV>){
		chomp;
		($orig_gff,$synonym) = split(/\s+/,$_); 
		foreach $regtype (@seqregion_types){
			$contig_slice = $slice_adaptor->fetch_by_region( $regtype, $synonym );
			last if(defined($contig_slice));
		}

		if(defined($contig_slice)){
	                #print "# adding synonym $orig_gff ($synonym)\n";
			$synonym{$orig_gff} = $synonym;
                	$contig_slice->add_synonym( $orig_gff );
			$n_of_synonyms++;
		} else {
			print "# cannot add unmatched synonym $orig_gff ($synonym)\n";
		}
	}
	close(TSV);	

	print "# total=$n_of_synonyms synonyms\n";

	# filter out sequences in GFF not matched in synonym list
	my $n_filtered_genes = 0;
	my %filter_scaffolds;
	my $filter_gff3file = $gff3_file . ".syn.filt";
	my $filter_log = $synonym_file.'.log';

	open(FILTERLOG,">",$filter_log) || die "# ERROR: cannot create $filter_log\n";

	open(FILTERGFF,">",$filter_gff3file) || die "# ERROR: cannot create $filter_gff3file\n";

	open(GFF,'<',$new_gff3file) || die "# ERROR: cannot read $new_gff3file\n";
        while(<GFF>){
		if(/^#/){ print FILTERGFF $_ }
                else{
                        my @gffdata = split(/\t/,$_);

			if(!defined($synonym{$gffdata[0]})){

				if($gffdata[2] eq 'gene'){ 
					print FILTERLOG "$gffdata[0] filtered out: $gffdata[8]";

					if(!$filter_scaffolds{$gffdata[0]}){
						$n_filtered_genes++;
						$filter_scaffolds{$gffdata[0]} = 1;
					}	
				}
				next;
			}
			
			print FILTERGFF;
		}
	}
	close(GFF);

	close(FILTERGFF);

	close(FILTERLOG);

	$new_gff3file = $filter_gff3file;

        print "\n# synonym-filtered GFF3 file: $filter_gff3file\n\n";
	print "# filtered sequences=$n_filtered_genes synonym-filtered log file: $filter_log\n\n";
}

## make custom-GFF subset with user-provided number of features (genes)
## and cut CDS sequences as internal control
if($max_feats > 0){

	my $num_of_features = 0;
	my ($chr,$CDS_start,$CDS_end,$CDS_strand,$CDS_seq,$CDS_name);
	$short_gff3file = $gff3_file . ".m$max_feats";

	# connect to production db
	my $registry = 'Bio::EnsEMBL::Registry';
	$registry->load_all($reg_file);
	my $slice_adaptor = $registry->get_adaptor($species, "core", "slice");

	if($check_gff_CDS){ print "\n# CDS sequences, internal check:\n\n" } 

        open(SHORTGFF,'>',$short_gff3file) || die "# ERROR: cannot create $short_gff3file\n";

	open(GFF,'<',$new_gff3file) || die "# ERROR: cannot read $new_gff3file\n";
        while(<GFF>){
		if(/^#/){ print SHORTGFF $_ }
		else{
			my @gffdata = split(/\t/,$_);
                        if($gffdata[2] eq 'gene'){ $num_of_features++ }
			elsif($check_gff_CDS && $gffdata[2] =~ m/mRNA/i){ 
				if(defined($CDS_seq) && $CDS_seq ne ''){ # print previous CDS
					print ">$CDS_strand|$CDS_name";
						
					if($CDS_strand eq '-'){ # get reverse complement
						$CDS_seq =~ tr/ACGTacgtyrkmYRKM/TGCAtgcarymkRYMK/; 
						$CDS_seq = reverse($CDS_seq);
					} #print "$CDS_seq\n"; 

					my $seq_obj = Bio::Seq->new(-seq => $CDS_seq, -alphabet => 'dna' );
					print $seq_obj->translate()->seq()."\n";
				}

				# CDS_name is set to parent mRNA, sequence is initialized
				($CDS_name,$CDS_seq) = ($gffdata[8],'');				
			}

			print SHORTGFF $_;

			if($check_gff_CDS == 1 && $gffdata[2] eq 'CDS'){ 
				#1 maker_ITAG CDS 30927	31259 .	- 0 ID=CDS:Soly...
				($chr,$CDS_start,$CDS_end,$CDS_strand) = @gffdata[0,3,4,6];
			
				# cut always leading strand, take rev comp just before translating	
				my $chr_slice = $slice_adaptor->fetch_by_region('chromosome',$chr,$CDS_start,$CDS_end);
				$CDS_seq .= $chr_slice->seq();
			}

			if($num_of_features > $max_feats){
				
				if($check_gff_CDS && $CDS_seq ne ''){ # print last CDS
                                        print ">$CDS_strand|$CDS_name";
                                        
					if($CDS_strand eq '-'){ # get reverse complement
                                                $CDS_seq =~ tr/ACGTacgtyrkmYRKM/TGCAtgcarymkRYMK/; 
                                                $CDS_seq = reverse($CDS_seq);
                                        } #print "|$CDS_seq|\n";

                                        my $seq_obj = Bio::Seq->new(-seq => $CDS_seq, -alphabet => 'dna' );
                                        print $seq_obj->translate()->seq()."\n";
                                }

				last;
			}
		}
	}
	close(GFF);

	close(SHORTGFF);

	$new_gff3file = $short_gff3file;

	print "\n# shortened GFF3 file: $short_gff3file\n";
}

## Update genebuild metadata
#########################################################################

# connect to production db
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_all($reg_file);
my $meta_adaptor = $registry->get_adaptor($species, "core", "MetaContainer");

my $year_month = strftime("%Y-%m", localtime());
my $year_month_day = strftime("%Y-%m-%d", localtime());

my $version_start_date = "$year_month-$gene_source";

my $genebuild_method = "Generated from $gene_source annotation";

print "\n# Setting meta genebuild.method to $genebuild_method \n";

if($meta_adaptor->single_value_by_key( 'genebuild.method' )){
        if(!$meta_adaptor->key_value_exists( 'genebuild.method', $genebuild_method )) {
                $meta_adaptor->update_key_value( 'genebuild.method', $genebuild_method );
        }
} else {
                $meta_adaptor->store_key_value( 'genebuild.method', $genebuild_method );
}

print "\n# Setting meta genebuild.version to $version_start_date\n";

if($meta_adaptor->single_value_by_key( 'genebuild.version' )){
	if(!$meta_adaptor->key_value_exists( 'genebuild.version', $version_start_date )) {
		$meta_adaptor->update_key_value( 'genebuild.version', $version_start_date );
	} 
} else {
		$meta_adaptor->store_key_value( 'genebuild.version', $version_start_date );
}

print "\n# Setting meta genebuild.start_date to $version_start_date\n";
if($meta_adaptor->single_value_by_key( 'genebuild.start_date' )){
	if(!$meta_adaptor->key_value_exists( 'genebuild.start_date', $version_start_date )) {
		$meta_adaptor->update_key_value( 'genebuild.start_date', $version_start_date );
	} 
} else {
	$meta_adaptor->store_key_value( 'genebuild.start_date', $version_start_date );
}

print "\n# Setting meta genebuild.last_geneset_update to $year_month_day\n\n";
if($meta_adaptor->single_value_by_key( 'genebuild.last_geneset_update' )){
	if(!$meta_adaptor->key_value_exists( 'genebuild.last_geneset_update', $year_month_day )) {
		$meta_adaptor->update_key_value( 'genebuild.last_geneset_update', $year_month_day );
	}
} else {
	$meta_adaptor->store_key_value( 'genebuild.last_geneset_update', $year_month_day )
}

## Run init script and produce a hive_db with all tasks to be carried out
#########################################################################

my $initcmd = "init_pipeline.pl Bio::EnsEMBL::EGPipeline::PipeConfig::LoadGFF3_conf ".
    	"$hive_args ".
    	"--registry $reg_file ".
        "--production_db $prodbname ".
    	"--pipeline_dir $pipeline_dir ".
    	"--species $species ".
    	"--gff3_file $new_gff3file ".
    	"--protein_fasta_file $protein_fasta_file ".
    	"--gene_source '$gene_source' ".
	#"-fix_models 0 ".
	"--hive_force_init $overwrite ";

if($add_to_previous){ $initcmd .= "-delete_existing 0 " }

if($skip_bad_types){ $initcmd .= "-types_complete 0 " }

if($names_stable){ $initcmd .= "-use_name_field stable_id " }

if($logicname){ $initcmd .= "--logic_name $logicname " }

if($otherfeats){ $initcmd .= "-db_type otherfeatures " }

if($gbkfile){ $initcmd .= "-genbank_file $gbkfile " }

print "# $initcmd\n\n";

if($rerun == 0){

	open(INITRUN,"$initcmd |") || die "# ERROR: cannot run $initcmd\n";
	while(<INITRUN>){
		print;
	}
	close(INITRUN);
}

## Send jobs to hive 
######################################################################### 

print "# hive job URL: $hive_url";

system("beekeeper.pl -url '$hive_url;reconnect_when_lost=1' -sync");
system("runWorker.pl -url '$hive_url;reconnect_when_lost=1' -reg_conf $reg_file");
system("beekeeper.pl -url '$hive_url;reconnect_when_lost=1' -reg_conf $reg_file -loop");

print "# hive job URL: $hive_url\n\n";
print "beekeeper.pl -url '$hive_url;reconnect_when_lost=1' -sync\n\n";
print "beekeeper.pl -url '$hive_url;reconnect_when_lost=1' -reg_conf $reg_file -loop\n\n";


