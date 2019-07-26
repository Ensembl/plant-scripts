##A cheat sheet to add new 10+ genomes, very specific for 10+
##Please Run this from the relevant path
#Usage:perl add_10_plus.pl [species]
use 5.14.0;
use warnings;
use Data::Dumper;
use Cwd;
use DBI;
use Tools::FileReader qw(slurp slurp_hash_list read_file file2hash file2hash_tab line2hash);

{
    my ($species) = @ARGV;
    if (@ARGV < 1){
        usage();
    }
    ##see what I did here with the Uppercase first?
    my $Species = ucfirst($species);

    ##make sure we are in the right cwd
    my $path = "/nfs/production/panda/ensemblgenomes/data/Plants/Wheat/$Species";
    say $path;
    my $dir = getcwd;
    if ($dir ne $path){
        die("Please run this from $path");
    }
    
    ##Get the fasta file
    my $fasta_file = "$species.genome.fa";

    ##Clean up the file:
    my $cmd = "perl $ENV{plant_tools}/new_core_scripts/clean_fasta.pl $fasta_file > clean_$species.fa";
    say "1. Cleaning up the fasta headers:\n$cmd";
    system($cmd);

    ##Break into chunks
    my $chunk_file = $species."_chunks.fa";
    $cmd = "perl $ENV{plant_tools}/new_core_scripts/fasta2agp/break_fasta.pl  clean_$species.fa  > $chunk_file";
    say "2. Breaking into chunks:\n$cmd";
    system($cmd);
    
    ##Create AGP file
    $cmd = "perl $ENV{plant_tools}/new_core_scripts/fasta2agp/fasta2agp.pl  $chunk_file  > $species.agp";
    say "3. Creating the agp file:\n$cmd";
    system($cmd);

    ##Create conf file (using lancer as a template
    my $file = "/nfs/production/panda/ensemblgenomes/data/Plants/Wheat/Lancer/lancer.conf";
    my $conf_file = "$path/$species.conf";
    open OUT, ">", $conf_file;
    my @lines = slurp($file);
    for my $line (@lines){
        ##Use lancer as a template for the meta table
        if ($line =~ /meta_source/){
            print OUT "meta_source\ttriticum_aestivum_lancer_core_45_98_1\n";
            next;
        }
        ##All other lines: swap lancer for our species (lower and upper case)
        $line =~ s/Lancer/$Species/;
        $line =~ s/lancer/$species/;
        print OUT $line,"\n";
    }
    close OUT;

    ##Nearly there: load the new assembly!
    $cmd = "perl $ENV{plant_tools}/new_core_scripts/add_new_core.pl --param_file=$conf_file";
    say "4. Loading the assembly into the DB:\n$cmd";
    system($cmd);

    ##Final SQL update for the meta table
    update_meta($conf_file,$Species,$species);
    
}

##======================================== 
sub update_meta {
##======================================== 
    my ($conf_file,$Species,$species) = @_;
    my ($h) = file2hash_tab($conf_file);
    my $dsn = "DBI:mysql:host=$h->{host};port=$h->{port}";
    my $dbh = DBI->connect($dsn, $h->{user}, $h->{pass});

    ##Update the meta table with our species name
    my $core = $h->{core};
    my $Sql = "UPDATE $h->{core}.meta SET meta_value = REPLACE(meta_value,'Lancer','$Species')";
    my $sql = "UPDATE $h->{core}.meta SET meta_value = REPLACE(meta_value,'lancer','$species')";
    
    ##Prepare and run the SQLs
    my $sth = $dbh->prepare($Sql);
    $sth->execute();
    
    my $sth = $dbh->prepare($sql);
    $sth->execute();
    $sth->finish();
}


sub usage {
    say "Usage perl add_10_plus.pl [a] [b]";
    exit 0;
}
 
