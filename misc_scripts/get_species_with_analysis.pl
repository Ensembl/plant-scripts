
#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);

my $GETDBSEXE = $ENV{'ENSAPIPATH'}."/ensembl-metadata/misc_scripts/get_list_databases_for_division.pl"; # see .bashrc
my $SQL       = "SELECT logic_name FROM analysis WHERE logic_name RLIKE 'XXX'";

my $prod_db_cmd  = 'mysql-ens-meta-prod-1';
my $stag_db_cmd  = 'mysql-ens-sta-3';
my $division     = 'Plants';

my ($logic_name, $ensembl_version, $help, $schema_name, $value, $query);

GetOptions(	
	"help|?" => \$help,
	"division|v=s" => \$division,
	"version|v=s"  => \$ensembl_version,
	"logic|l=s"    => \$logic_name,
	"prodcmd|p=s"  => \$prod_db_cmd,    
	"stagcmd|s=s"  => \$stag_db_cmd
) || help_message(); 

if($help){ help_message() }

sub help_message {
	print "\nusage: $0 [options]\n\n".
		"-v next Ensembl version      (required, example: -v 95)\n".
		"-l logic_name of analysis    (required, example: -l gramene_plant_reactome)\n";
		"-d Ensembl division          (optional, default: -d Plants)\n".
		"-p production db command     (optional, default: -p $prod_db_cmd)\n".
		"-s staging db command        (optional, default: -s $stag_db_cmd)\n\n";
	exit(0);
}

if(!$ensembl_version){
	die "# ERROR: need argument -v, example -v 95\n\n";
}

if(!$logic_name){
	   die "# ERROR: need argument -l, example -l gramene_plant_reactome\n\n";
} else {
	$query = $SQL;
	$query =~ s/XXX/$logic_name/;
}


open(GETDBS,"perl $GETDBSEXE \$($prod_db_cmd details script) -release $ensembl_version -division $division |") ||
	die "# ERROR: cannot run $GETDBSEXE \$($prod_db_cmd details script) -release $ensembl_version -division $division\n";
while(<GETDBS>){
	$schema_name = (split)[0]; 
	next if($schema_name !~ /core/ && $schema_name !~ /otherfeat/);
	
	# query this schema to find out assembly provider, accession and genebuild.version
	my (@values);
	open(SQL,"$stag_db_cmd $schema_name -e \"$query\" |") ||
		die "# ERROR: cannot query $stag_db_cmd $schema_name\n\n";
	while(<SQL>){
		chomp;
		$value = $_;  
	   next if($value eq 'logic_name');
		push(@values,$value);
	}
	close(SQL);

	# print all
	printf("%s\t%s\n", 
		$schema_name, 
		join("\t",@values)) if(@values);
}
close(GETDBS);
