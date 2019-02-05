##To run: perl add_new_core --param_file=[PARAM_FILE]
##you can find param files in param_file_examples dir
use 5.14.0;
use warnings;
use Data::Dumper;
use lib '/nfs/production/panda/ensemblgenomes/development/gnaamati/lib';
use FileReader qw(slurp slurp_hash_list read_file file2hash file2hash_tab line2hash);
use Getopt::Long;
use Pod::Usage;
my $core;
my $file2;
my $param_file;

{
    GetOptions ("param_file=s" => \$param_file,
                "file2=s" => \$file2) 
    or die("Incorrect Usage");

    if (!$param_file){
        usage();
    }

    my $h   = file2hash_tab($param_file);
    my $dbh = get_dbh($h);

    ##creating db and adding tables
    create_db($h);

    ##Adding controlled vocab
    add_cv($h);

    ##Loading Fasta data
    load_fasta($h);

    ##Loading AGP data
    load_agp($h);

    ##updating top level attrib
    set_top_level($h);

    ##updating meta table
    copy_meta($h, $dbh);
    
}

#======================================== 
sub create_db {
#======================================== 
    my ($h) = @_;

    ##Creating the DB and adding tables
    warn "Creating and populating new core for $h->{core}\n";
    my $cmd = "mysqladmin -h $h->{host} -P $h->{port} -u$h->{user} -p$h->{pass} CREATE $h->{core}";
    system($cmd);
    
    ##Adding tables
    $cmd = "mysql -h$h->{host} -P$h->{port} -u$h->{user} -p$h->{pass} $h->{core} < $ENV{ENSEMBL_ROOT_DIR}/ensembl/sql/table.sql";
    system($cmd);
}


#======================================== 
sub add_cv {
#======================================== 
    my ($h) = @_;
    warn "Adding controlled vocabulary for $core\n";
    
    my $path = '/nfs/panda/ensemblgenomes/apis/ensembl/master/ensembl-production/scripts/production_database';    
    system("mkdir /tmp/prod_db_tables");
    my $cmd = "perl $path/populate_production_db_tables.pl ".
              "--host $h->{host} --port $h->{port} --user $h->{user} --pass $h->{pass} ".
              '$(mysql-pan-prod-ensrw details prefix_m) '.
              "--database $h->{core} ".
              "--dumppath /tmp/prod_db_tables --dropbaks";
    system($cmd);
}

#======================================== 
sub set_top_level {
#======================================== 
    my ($h) = @_;
    my $path = "$ENV{ENSEMBL_ROOT_DIR}/ensembl-pipeline/scripts";
    my $cmd = "perl $path/set_toplevel.pl ".
              "--host $h->{host} --port $h->{port} --user $h->{user} --pass $h->{pass} ".
              "--dbname $h->{core} ";
    say $cmd;
    #system($cmd);
}

#======================================== 
sub load_fasta {
#======================================== 
    my ($h) = @_;
    my $path = "$ENV{ENSEMBL_ROOT_DIR}/ensembl-pipeline/scripts";
    my $cmd = "perl $path/load_seq_region.pl ".
              "--host $h->{host} --port $h->{port} --user $h->{user} --pass $h->{pass} ".
              "--dbname $h->{core} ".
              "--coord_system_name scaffold --coord_system_version $h->{version} ".
              "--rank 2 -default_version -sequence_level ".
              "--fasta_file $h->{fasta_file}";
    say $cmd,"\n";
    system($cmd);
}

#======================================== 
sub load_agp {
#======================================== 
    my ($h) = @_;
    my $path = "$ENV{ENSEMBL_ROOT_DIR}/ensembl-pipeline/scripts";
    
    ##Load AGP part 1
    my $cmd = "perl $path/load_seq_region.pl ".
              "--host $h->{host} --port $h->{port} --user $h->{user} --pass $h->{pass} ".
              "--dbname $h->{core} ".
              "--coord_system_name chromosome --coord_system_version $h->{version} ".
              "--rank 1 --default_version ".
              "-agp_file $h->{agp_file}";
    say $cmd,"\n";
    system($cmd);

    ##Load AGP part 2
    $cmd = "perl $path/load_agp.pl ".
            "--host $h->{host} --port $h->{port} --user $h->{user} --pass $h->{pass} ".
            "--dbname $h->{core} ".
            "--assembled_name chromosome -assembled_version $h->{version} ".
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
    
    $sql = "delete from $h->{'core'}.meta";
    $sth = $dbh->prepare($sql);
    $sth->execute();

    $sql = qq{
    insert into $h->{'core'}.meta
    (select * from $h->{'meta_source'}.meta)
    };
    $sth = $dbh->prepare($sql);
    $sth->execute();

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

