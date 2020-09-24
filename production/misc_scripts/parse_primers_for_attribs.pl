# perl parse_primers.pl primer_files/Kronos_primers.csv ems_name2id_kronos.txt > sqls_to_update.sql
use 5.14.0;
use warnings;
use Data::Dumper;
use lib '/nfs/production/panda/ensemblgenomes/development/gnaamati/lib';
use FileReader qw(slurp slurp_hash_list read_file file2hash file2hash_tab line2hash_comma);
use Tools::Chrom2Id qw(chrom2id_wheat);

##Need to update this according to the attrib table in the DB
my $attrib_ref = {
    link_to_primer          => 610,
    snp_type                => 611,
    total_contigs           => 612,
    primer_type             => 613,
    ems_genotype            => 614,
    mutant_quality          => 659,
    residual_heterogeneity  => 660,

};


{
    my ($primers, $name2id_file) = @ARGV;
    if (@ARGV != 2){
        usage();
    }
    my $id_count = 1;

    ##Get all variation ids according to variation name
    my $name2id_ref = file2hash_tab($name2id_file);
    
    my @lines = slurp($primers);
    my $header = shift(@lines);
    my $count = 0;
    for my $line (@lines){
        my $h = line2hash_comma($header, $line);

        ##Update names in input file to the attrib names
        $h->{'link_to_primer'} = $h->{marker};
        $h->{'ems_genotype'}   = $h->{'ems-gt'};
        $h->{'ems_genotype'}   = $h->{'ems.gt'};
        $h->{'mutant_quality'} = $h->{'mut_qual'};
        $h->{'residual_heterogeneity'} = $h->{'rh'};

        #Get variation dbid from variation name by primer id
        my $var_id = $name2id_ref->{$h->{marker}};
        
        ##Don't update if no var ID (probably a deletion/insertion mutant which we don't display
        if (!$var_id){
            next;
        }
        
        ##Don't link out to primer when we have no common field (with the KASP marker data)
        if (!$h->{common}){
            $h->{link_to_primer} = '';
        }
        
        ##Print SQLs
        my @types = qw/link_to_primer snp_type total_contigs primer_type ems_genotype mutant_quality residual_heterogeneity/;
        for my $type (@types){
            get_sql($type, $h, $var_id);
        }

        ##Simple counter to track progress of updates
        if ($id_count % 100_000 == 0){
            warn $id_count;
        }
        $id_count++;
    }
    
}

##Print the SQL according to the primer type
sub get_sql {
    my ($type, $h, $var_id) = @_;
    
    my $primer_type = $h->{$type};
    
    ##Special case for residual_heterogeneity
    if ($type eq 'residual_heterogeneity'){
        if ($primer_type){
            $primer_type = 'Yes';
        }
        else{
            $primer_type = 'No';
        }
    }

    ##Don't print SQL if primers are empty
    if (!$primer_type){
        return;
    }
    
    ##Longer names for Hom and Het
    if ($primer_type eq 'Hom'){
        $primer_type = 'homozygous';
    }
    if ($primer_type eq 'Het'){
        $primer_type = 'heterozygous';
    }

    my $sql = qq{
    insert into variation_attrib 
        (variation_id,attrib_id,value)
    values
        ($var_id, $attrib_ref->{$type},"$primer_type");
    };
    say $sql;
}



sub usage {
    say "Usage perl parse_primers.pl [a] [b]";
    exit 0;
}
 
