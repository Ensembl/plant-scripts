#!/usr/bin/env perl

## To run: perl add_new_core --param_file=[PARAM_FILE]
## You can find param files in param_file_examples dir
#
## Guy Gnaamati, Bruno Contreras Moreira 2019-20

use 5.14.0;
use warnings;
use FindBin '$Bin';
use lib "$Bin/..";
use Tools::FileReader qw( file2hash_tab );
use File::Temp qw( tempdir );
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use DBI;
use JSON qw( decode_json );
use HTTP::Tiny;

# should contain ensembl, ensembl-production, ensembl-pipelines
my $ENSEMBLPATH = $ENV{ENSEMBL_ROOT_DIR}; 

# alias for pan production db server with rw permissions
my $PANPRODSERVER = 'mysql-eg-pan-prod-ensrw';

# REST configuration, used to get taxonomy
my $RESTURL     = 'http://rest.ensembl.org';
my $TAXOPOINT   = $RESTURL.'/taxonomy/classification/';

my $VERBOSE = 1;

# only required if meta table is to be created from scratch, not copied over
my @METAKEYS = qw( 
    provider.name provider.url
    assembly.accession assembly.name 
    species.division species.production_name species.display_name species.taxonomy_id species.strain
);

my %DERIVED_METAKEYS = (
	'assembly.name' => 
		['assembly.default'],
	'species.production_name' => 
		['species.db_name'],
	'species.display_name' => 
		['species.scientific_name','species.species_name','species.wikipedia_name','species.url','species.wikipedia_url'],
	'species.taxonomy_id' => 
		['species.species_taxonomy_id']
);

my %OPTIONAL_METAKEYS = ( 
	'species.strain' => 1,
	'provider.url' => 1 
);

############################################################################

my $core;
my $file2;
my $param_file;

{
    GetOptions (
	 	"param_file=s" => \$param_file,
      ) or die("Incorrect Usage");

    if (!$param_file){
        usage();
    }

    ## read param file
    my $h   = file2hash_tab($param_file);

    ## check db server connection details
    if($h->{'host'}){
    	if(!$h->{'user'}){ 
    		my $server_args;
    		chomp( $server_args = `$h->{'host'} details` );
    		if($server_args =~ m/--host=(\S+) --port=(\S+) --user=(\S+) --pass=(\S+)/){
				$h->{'host'} = $1;
    			$h->{'port'} = $2;
				$h->{'user'} = $3;
				$h->{'pass'} = $4;
      	} else {
    			die "# ERROR: please set port, user & password in param_file\n";
    		}
    	}
    } else {
    	die "# ERROR: please set host in param_file\n";
	 }

   ## check assembly.name is set
   if(!$h->{'assembly.name'}){
   	die "# ERROR: please set assembly.name in param_file\n";
   }

   ##connect to db server 
   my $dbh = get_dbh($h);

   ##creating db and adding tables
   create_db($h);

   ##Adding controlled vocab
   add_cv($h);

   ##Loading Fasta data
	# will issue harmless warnings:
	# -------------------- WARNING ----------------------
	# MSG: Name 2C_43 does not look like a valid accession - are you sure this is what you want?
	# FILE: ensembl-pipeline/scripts/load_seq_region.pl LINE: 258
   load_fasta($h);

   ##Loading AGP data
   if($h->{'agp_file'}){
      load_agp($h);
   }

   ##updating meta table
   if($h->{'meta_source'}){
   	copy_meta($h, $dbh); # requires manual tweaking
   } else {
   	workout_meta($h, $dbh, $TAXOPOINT);	
   }

   ##Add seq region attribs
   add_seq_region_attribs($h, $dbh);
}

#======================================== 
sub create_db {

#======================================== 
    my ($h) = @_;

    ##Creating the DB and adding tables
    warn "# create_db: creating and populating new core for $h->{core}\n";
    my $cmd = "mysqladmin -h $h->{host} -P $h->{port} -u$h->{user} -p$h->{pass} CREATE $h->{core}";
    open(SQL,"$cmd 2>&1 |") || die "# ERROR(create_db): cannot run $cmd\n";
    while(<SQL>){
    	if(/mysqladmin: CREATE DATABASE failed/){
    		die "# ERROR (create_db): $h->{core} already exists, remove it and re-run\n";
    	}
    }
    
    ##Adding tables
    $cmd = "mysql -h$h->{host} -P$h->{port} -u$h->{user} -p$h->{pass} $h->{core} < $ENSEMBLPATH/ensembl/sql/table.sql";
    system($cmd);
}


#======================================== 
sub add_cv {
#======================================== 
    my ($h) = @_;
    warn "# add_cv : adding controlled vocabulary for $h->{'core'}\n";
    
    my $path = "$ENSEMBLPATH/ensembl-production/scripts/production_database";    
	 my $tmpdir = tempdir( CLEANUP => 1 );
    my $cmd = "perl $path/populate_production_db_tables.pl ".
              "--host $h->{host} --port $h->{port} --user $h->{user} --pass $h->{pass} ".
              "\$($PANPRODSERVER details prefix_m) ".
              "--database $h->{core} ".
              "--dumppath $tmpdir --dropbaks";
    system($cmd);

	 warn "# add_cv: done\n\n";
}

#======================================== 
sub set_top_level {
#======================================== 
    my ($h) = @_;
    my $path = "$ENSEMBLPATH/ensembl-pipeline/scripts";
    my $cmd = "perl $path/set_toplevel.pl ".
              "--host $h->{host} --port $h->{port} --user $h->{user} --pass $h->{pass} ".
              "--dbname $h->{core} ";
    say $cmd;
    system($cmd);
}

#======================================== 
sub load_fasta {

# loads chunks as virtual scaffolds	
#======================================== 
    my ($h) = @_;
    my $path = "$ENSEMBLPATH/ensembl-pipeline/scripts";
    my $cmd = "perl $path/load_seq_region.pl ".
              "--host $h->{host} --port $h->{port} --user $h->{user} --pass $h->{pass} ".
              "--dbname $h->{core} ".
              "--coord_system_name scaffold --coord_system_version $h->{'assembly.name'} ".
              "--rank 2 -default_version -sequence_level ".
              "--fasta_file $h->{fasta_file}";
    say $cmd,"\n";
    system($cmd);
}

#======================================== 
sub load_agp {
#======================================== 
    my ($h) = @_;
    my $path = "$ENSEMBLPATH/ensembl-pipeline/scripts";
    
    ##Load AGP part 1
    my $cmd = "perl $path/load_seq_region.pl ".
              "--host $h->{host} --port $h->{port} --user $h->{user} --pass $h->{pass} ".
              "--dbname $h->{core} ".
              "--coord_system_name chromosome --coord_system_version $h->{'assembly.name'} ".
              "--rank 1 --default_version ".
              "-agp_file $h->{agp_file}";
    say $cmd,"\n";
    system($cmd);

    ##Load AGP part 2
    $cmd = "perl $path/load_agp.pl ".
            "--host $h->{host} --port $h->{port} --user $h->{user} --pass $h->{pass} ".
            "--dbname $h->{core} ".
            "--assembled_name chromosome -assembled_version $h->{'assembly.name'} ".
            "--component_name scaffold ".
            "-agp_file $h->{agp_file}";
    say $cmd,"\n";
    system($cmd);

}

#======================================== 
sub copy_meta {
#======================================== 
    my ($h, $dbh) = @_;
    my ($sql, $sth);
    
    ##Deleting current meta
    $sql = "delete from $h->{'core'}.meta";
    $sth = $dbh->prepare($sql);
    $sth->execute();

    ##Copying meta from other source
    $sql = qq{
    insert into $h->{'core'}.meta
    (select * from $h->{'meta_source'}.meta)
    };
    $sth = $dbh->prepare($sql);
    $sth->execute();

    ##Add test regarding the version

    ##Cleaning up meta table (repeats and patches)
    $sql = "delete from $h->{'core'}.meta where meta_key rlike ";
    my @keys = qw/patch repeat/;
    for my $key (@keys){
        my $sql_to_run = $sql."'$key'";
        $sth = $dbh->prepare($sql_to_run);
        $sth->execute();
    }

}

#========================================
sub workout_meta {

# Work out the key meta data for this core db.
# Mandatory params: accession, version, production_name, display_name, taxonomy_id, division, provider_name
# Optional params: strain, provider_url

#======================================== 

	my ($h, $dbh, $rest_entry_point) = @_;
	my ($sql, $sth, $metakey, $derived_key, $derived_value);

	foreach $metakey (@METAKEYS) {

		if(!$h->{$metakey}){
			if($OPTIONAL_METAKEYS{$metakey}){ 
				next 
			} else{
				die "# ERROR (workout_meta) : please set param '$metakey' in param_file\n";
			}
		}

	   warn "# workout_meta: $metakey=$h->{$metakey}\n" if($VERBOSE);
	   $sql = qq{INSERT INTO $h->{'core'}.meta (species_id,meta_key,meta_value) VALUES (1, '$metakey', '$h->{$metakey}');};
	   $sth = $dbh->prepare($sql);
	   $sth->execute();

		# check derived meta keys
		if($DERIVED_METAKEYS{$metakey}){
			foreach $derived_key (@{ $DERIVED_METAKEYS{$metakey}  }){

				$derived_value = $h->{$metakey};
				if($derived_key =~ m/url/){ $derived_value =~ s/\s/_/g }

				$sql = qq{INSERT INTO $h->{'core'}.meta (species_id,meta_key,meta_value) VALUES (1, '$derived_key', '$derived_value');};
				$sth = $dbh->prepare($sql);
				$sth->execute();
			}
		}	
	
		# obtain full taxonomy for passed taxonomy_id from Ensembl REST interface
		if($metakey eq 'species.taxonomy_id'){
			my $http = HTTP::Tiny->new();
			my $request = $rest_entry_point.$h->{$metakey};
			my $response = $http->get($request, {headers => {'Content-Type' => 'application/json'}});
			if($response->{success} && $response->{content}){
				my $taxondump = decode_json($response->{content});
				foreach my $taxon (@{ $taxondump }) {
					next if(!$taxon->{'name'});
					warn "# workout_meta: species.classification, $taxon->{'name'}\n" if($VERBOSE);
					$sql = qq{INSERT INTO $h->{'core'}.meta (species_id,meta_key,meta_value) VALUES (1, 'species.classification', '$taxon->{'name'}');};
					$sth = $dbh->prepare($sql);
					$sth->execute();
				}
			} else {
				warn "# ERROR (workout_meta) : $request request failed, please re-run\n";
			}
		} 
	}
}

#======================================== 
sub add_seq_region_attribs { 
#======================================== 
    my ($h, $dbh) = @_;
    my ($sql, $sth);

    print Dumper $dbh;

    my $seq_region_file = $h->{seq_region_file};
    my $core = $h->{core};
    
    ##Get the seq_regions
    $sql = "select seq_region_id, name from $core.seq_region where coord_system_id=2 order by name asc;";
    $sth = $dbh->prepare($sql);
    $sth->execute();

    my $rank = 1;
    while (my $ref = $sth->fetchrow_hashref()) {
        my ($seq_region_id, $name) = ($ref->{seq_region_id},$ref->{name});

        my $comp;
        if ($name =~ /\d(\w)/){
            $comp = $1;
        }
        elsif ($name = 'Un'){
            $comp = 'U';
        }
        else{
            say "no $comp for $name";
        }

        ##Insert polyploid value
        my $sql = qq{
            insert into $core.seq_region_attrib
                (seq_region_id, attrib_type_id, value)
            values
                ($seq_region_id, 425, '$comp');
        };
        run_sql($dbh,$sql,$core);

        ##Insert top level
        $sql = qq{
            insert into $core.seq_region_attrib
                (seq_region_id, attrib_type_id, value)
            values
                ($seq_region_id, 6, '1');
        };
        run_sql($dbh,$sql,$core);


        ##Insert karyotype rank
        $sql = qq{
            insert into $core.seq_region_attrib
                (seq_region_id, attrib_type_id, value)
            values
                ($seq_region_id, 367, $rank);
        };
        run_sql($dbh,$sql,$core);
        $rank++;
    }

    $sth->finish();
}

#======================================== 
sub run_sql {
#======================================== 
    my ($dbh, $sql, $core) = @_;
    
    say "running:\n $sql";
    my $sth = $dbh->prepare($sql);
    $sth->execute();
    $sth->finish();
}


#======================================== 
sub get_params {
#======================================== 
    my ($file) = @_;
    my $h = file2hash_tab($param_file);

    my $user = $h->{user}; 
    my $pass = $h->{pass}; 
    my $host = $h->{host}; 
    my $port = $h->{port}; 
    my $core = $h->{core};
    return ($user,$pass,$host,$port,$core);

}

#======================================== 
sub get_dbh {
#======================================== 
    my ($h) = @_;
    my $dsn = "DBI:mysql:host=$h->{host};port=$h->{port}";
    my $dbh = DBI->connect($dsn, $h->{user}, $h->{pass});
    return $dbh;
}

#======================================== 
sub usage {
#======================================== 
    say "Usage perl add_new_core --param_file=[PARAM_FILE]";
    exit 0;
}

