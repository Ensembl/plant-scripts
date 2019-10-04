##perl update_genomic_align.pl iwgsc_db_data.txt iwgsc_alignment.dat
##perl update_genomic_align.pl zavitan_db_data.txt zavitan_alignment.dat
use 5.14.0;
use warnings;
use Data::Dumper;
my $mlss = 9822;
use Carp qw(croak);

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

    }

}


##Slurps a file into a list
sub slurp {
    my ($file) = @_;
    my (@data, @data_chomped);
    open IN, "<", $file or croak "can't open $file\n";
    @data = <IN>;
    for my $line (@data){
        chomp($line);
        push (@data_chomped, $line);
    }
    close IN;
    return (@data_chomped);
}

##Splits a line and parses it into the hash according to given headers
sub line2hash {
    my ($header, $line) = @_;

    ##Split header fields and line fields
    my @header_fields = split (/\t/, $header);
    my @line_fields   = split (/\t/, $line);

    ##Fill the hash
    my $hash = {};
    my $i = 0;
    for my $fld (@header_fields){
        $fld =~ s/\s/_/g;
        $hash->{lc($fld)} = $line_fields[$i];
        $i++;
    }

    return $hash;
}





sub usage {
    say "Usage perl update_genomic_align.pl [a] [b]";
    exit 0;
}



 
