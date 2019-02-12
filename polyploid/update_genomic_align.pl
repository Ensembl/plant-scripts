##perl update_genomic_align.pl iwgsc_db_data.txt iwgsc_alignment.dat
##perl update_genomic_align.pl zavitan_db_data.txt zavitan_alignment.dat
use 5.14.0;
use warnings;
use Data::Dumper;
use lib '/nfs/production/panda/ensemblgenomes/development/gnaamati/lib';
use FileReader qw(slurp slurp_hash_list read_file file2hash file2hash_tab line2hash);
use DBI;
my $mlss = 9822;

##The genomic align id start (prefix is mlss_id)
my $gab_id_count  = 98220000000000;
my $ga_id_count   = 98220000000000;
               
##GAB = genomic align block
##GA  = genomic align
{
    my ($db_data, $ihv_data) = @ARGV;
    if (@ARGV < 2){
        usage();
    }
    my $ihv_hash = {};
    my $gab_new_id_ref = {};
    my $new_gab_id;
    my $new_ga_id;

    
    ##Get all the group ids that contain homeoulugues (from iwgsc_alignment.dat)
    ## Call this IHV hash (interhomologous variants)
    my @lines = slurp($ihv_data);
    for my $line (@lines){
        my ($a,$b,$group_id) = split(/\t/,$line);
        $ihv_hash->{$group_id}++;
    }

    ##Go over all the genomic align and genomic align block ids
    @lines = slurp($db_data);
    my $header = shift(@lines);
    my $count = 0;
    for my $line (@lines){
        my $h = line2hash($header, $line);

        ##skip if not in IHV hash
        if (!$ihv_hash->{$h->{group_id}}){
            next;
        }

        ##Check if we already created a new GAB id, if not create a new one
        my $gab_id = $h->{genomic_align_block_id};
        if ($gab_new_id_ref->{$gab_id}){
            $new_gab_id = $gab_new_id_ref->{$gab_id};
        }
        else{   
            $new_gab_id = ++$gab_id_count;
            $gab_new_id_ref->{$gab_id} = $new_gab_id;
            
            ##Create SQL for new GAB
            my $sql = qq{
            insert into genomic_align_block 
                (genomic_align_block_id,method_link_species_set_id,score,perc_id,length,group_id,level_id)
            values
                ($new_gab_id,$mlss,$h->{score},$h->{perc_id},$h->{length},$h->{group_id},$h->{level_id});
            };
            say $sql;
        }



        ##SQL for new genmic align entry
        $new_ga_id = ++$ga_id_count;
        my $sql = qq{
        insert into genomic_align
            (genomic_align_id, genomic_align_block_id, method_link_species_set_id, dnafrag_id, dnafrag_start, dnafrag_end, dnafrag_strand, cigar_line, visible, node_id)
        values
            ($new_ga_id, $new_gab_id,$mlss,$h->{dnafrag_id}, $h->{dnafrag_start}, $h->{dnafrag_end},$h->{dnafrag_strand},'$h->{cigar_line}', $h->{visible},$h->{node_id});
        };
        say $sql;

        $count++;
        #if ($count == 6){
        #    die;
        #}
    }

}

sub get_gab_id {
    



}

sub get_max_gab_id {

    return "42";

    my $dsn = "DBI:mysql:database=ensembl_compara_plants_40_93;host=mysql-eg-prod-1.ebi.ac.uk;port=4238";

    my $dbh = DBI->connect($dsn, 'ensrw', 'writ3rp1');


    # now retrieve data from the table.
    my $sth = $dbh->prepare("select max(genomic_align_block_id) from genomic_align_block;");
    $sth->execute();
    my @row = $sth->fetchrow_array();
    my $max = $row[0];
    $sth->finish();
    $dbh->disconnect();

    return $max;

}

sub usage {
    say "Usage perl update_genomic_align.pl [a] [b]";
    exit 0;
}



 
