#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# Release Coordinator, please update this file before starting every release
# and check the changes back into GIT for everyone's benefit.

# Things that normally need updating are:
#
# 1. Release number
# 2. Check the name prefix of all databases
# 3. Possibly add entries for core databases that are still on genebuilders' servers

use strict;
use warnings;
use Bio::EnsEMBL::Utils::ConfigRegistry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyDBAdaptor;

my $curr_release = 97;
my $prev_release = $curr_release - 1;

# ---------------------- CURRENT CORE DATABASES----------------------------------

# most cores are on EG servers, but some are on ensembl's vertannot-staging
# explained at https://www.ensembl.org/info/docs/Doxygen/core-api/classBio_1_1EnsEMBL_1_1Registry.html
Bio::EnsEMBL::Registry->load_registry_from_url("mysql://ensro\@mysql-ens-plants-prod-2:4208/$curr_release");
#Bio::EnsEMBL::Registry->remove_DBAdaptor('saccharomyces_cerevisiae', 'core'); # never use EG's version of yeast
#Bio::EnsEMBL::Registry->load_registry_from_url("mysql://ensro\@mysql-ens-vertannot-staging:4573/$curr_release");

# ---------------------- PREVIOUS CORE DATABASES---------------------------------

# previous release core databases will be required by LoadMembers only
# !!! COMMENT THIS SECTION OUT FOR ALL OTHER PIPELINES (for speed) !!!

#my $suffix_separator = '__cut_here__';
#Bio::EnsEMBL::Registry->load_registry_from_db(
#    -host   => 'mysql-ens-mirror-3',
#    -port   => 4275,
#    -user   => 'ensro',
#    -pass   => '',
#    -db_version     => $prev_release,
#    -species_suffix => $suffix_separator.$prev_release,
#);
#Bio::EnsEMBL::Registry->remove_DBAdaptor('saccharomyces_cerevisiae'.$suffix_separator.$prev_release, 'core'); # never use EG's version of yeast
#Bio::EnsEMBL::Registry->load_registry_from_db(
#    -host   => 'mysql-ens-mirror-1',
#    -port   => 4240,
#    -user   => 'ensro',
#    -pass   => '',
#    -db_version     => $prev_release,
#    -species_suffix => $suffix_separator.$prev_release,
#);

#------------------------COMPARA DATABASE LOCATIONS----------------------------------


my $compara_dbs = {
    # general compara dbs
    'compara_master' => [ 'mysql-ens-plants-prod-2', 'pan_compara_master' ],
    'compara_curr'   => [ 'mysql-ens-plants-prod-2', 'ensembl_compara_pan_homology_44_97' ],
    'compara_prev'   => [ 'mysql-ens-sta-3', 'ensembl_compara_pan_homology_44_97' ],

    # homology dbs
    'compara_members' => [ 'mysql-ens-hive-prod-2', 'bcontreras_pan_load_members_97' ],
    'compara_ptrees'  => [ 'mysql-ens-hive-prod-2', 'ensembl_compara_pan_hom_44_97' ],

    # LASTZ dbs
    #'lastz_a' => [ 'mysql-ens-compara-prod-8', 'carlac_plants_lastz_batch1_95' ],
    #'lastz_b' => [ 'mysql-ens-compara-prod-6', 'muffato_plants_lastz_b_95' ],
    #'lastz_c' => [ 'mysql-ens-compara-prod-6', 'muffato_plants_lastz_c_95' ],
    #'lastz_lang_rerun' => [ 'mysql-ens-compara-prod-8', 'carlac_plants_lastz_lang_rerun_95' ],

    # synteny
    #'compara_syntenies' => [ 'mysql-ens-compara-prod-5', 'mateus_synteny_plants_42_95' ],
}; 

add_compara_dbs( $compara_dbs ); # NOTE: by default, '%_prev' dbs will have a read-only connection

# ----------------------NON-COMPARA DATABASES------------------------

# NCBI taxonomy database (also maintained by production team):
Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyDBAdaptor->new(
    -host => 'mysql-ens-sta-1',
    -user => 'ensro',
    -port => 4519,
    -group => 'taxonomy',
    -species => 'ncbi_taxonomy',
    -dbname => 'ncbi_taxonomy',
);

# -------------------------------------------------------------------

sub add_compara_dbs {
    my $compara_dbs = shift;

    foreach my $alias_name ( keys %$compara_dbs ) {
        my ( $host, $db_name ) = @{ $compara_dbs->{$alias_name} };

        my ( $user, $pass ) = ( 'ensro', '' );
        if($alias_name =~ m/compara_curr/){
		my $details = `$host-ensrw details`;
		if($details =~ m/--pass=(\S+)/){
		        ( $user, $pass ) = ( 'ensrw', $1 );
		}
	}

        Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
            -host => $host,
            -user => $user,
            -pass => $pass,
            -port => get_port($host),
            -species => $alias_name,
            -dbname  => $db_name,
        );
    }
}

sub get_port {
    my $host = shift;
    my $port = `$host port`;
    chomp $port;
    return $port;
}

1;
