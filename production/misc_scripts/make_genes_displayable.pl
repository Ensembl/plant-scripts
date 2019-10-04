#!/bin/env perl
use strict;
use warnings;
use Getopt::Std;

# Adds analysis description and makes genes displayable to a core db
# Ideally this data should first be added to production db
#
# Example: perl _make_genes_displayable.pl -s glycine_max -p mysql-ens-plants-prod-1-ensrw
#
# Bruno Contreras Moreira EMBL-EBI 2019

my (%opts,$coredb,$prod_server,$description,$gene_source);

my $web_data = q( {\'colour_key\' => \'[biotype]\',\'caption\' => \'Genes\',\'name\' => \'Genes\',\'label_key\' => \'[biotype]\',\'default\' => {\'MultiTop\' => \'gene_label\',\'contigviewbottom\' => \'transcript_label\',\'MultiBottom\' => \'collapsed_label\',\'contigviewtop\' => \'gene_label\',\'cytoview\' => \'gene_label\',\'alignsliceviewbottom\' => \'as_collapsed_label\'},\'key\' => \'ensembl\'} );

my $default_description_suffix = '. See the species \"About\" page for more information.';

getopts('hs:p:d:g:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)){
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-s species_core_db                             (required, example: -s theobroma_cacao_core_43_96_2)\n";
  print "-g gene_source                                 (required, example: -g INRA)\n";
  print "-d description                                 (required, example: -d 'Genes imported from INRA')\n";
  print "-p production db server, where data is loaded  (required, example: -p mysql-ens-plants-prod-1-ensrw)\n";
  exit(0);
}

if($opts{'s'}){ $coredb = $opts{'s'} } 
else{ die "# EXIT : need a valid -s species_name, such as -s theobroma_cacao_core_43_96_2\n" }

if($opts{'g'}){ $gene_source = $opts{'g'} } 
else{ die "# EXIT : need a valid -g gene_source, such as -g INRA\n" }

if($opts{'d'}){ 
	$description = $opts{'d'}; 
	$description .= $default_description_suffix;
}
else{ die "# EXIT : need a valid -d description, such as -d 'Genes imported from INRA'\n" }

if($opts{'p'}){ $prod_server = $opts{'p'} }
else{ die "# EXIT : need a valid -p production db server, such as -p mysql-ens-plants-prod-1-ensrw\n" }


## now check current production server

my @analysis_ids = `$prod_server $coredb -e "select distinct(analysis_id) from gene where source = '$gene_source'"`;
if(scalar(@analysis_ids) > 2){
	die "# EXIT: there are several analysis_id values, skip them\n";
} if(scalar(@analysis_ids) == 0){
        die "# EXIT: cannot find gene_source: $gene_source\n";
} else {
	print "\n# matched analysis_id: $analysis_ids[1]\n";
}


foreach my $analysis_id (@analysis_ids){
	next if($analysis_id eq 'analysis_id');

	print `$prod_server $coredb -e "UPDATE analysis_description SET displayable=1 WHERE analysis_id=$analysis_id;"`;
	print `$prod_server $coredb -e "UPDATE analysis_description SET description='$description' WHERE analysis_id=$analysis_id;"`;
	print `$prod_server $coredb -e "UPDATE analysis_description SET display_label='Genes' WHERE analysis_id=$analysis_id;"`;
	print `$prod_server $coredb -e "UPDATE analysis_description SET web_data='$web_data' WHERE analysis_id=$analysis_id;"`;
}
