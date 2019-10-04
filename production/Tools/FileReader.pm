###======================================== 
#A set of subroutines for file/api handlind and parsing 
###======================================== 
package Tools::FileReader;
require Exporter;
# symbols to export on request
@ISA = qw(Exporter);
@EXPORT_OK = qw(slurp slurp_and_chomp read_file write_file file2hash file2hash_tab hash2file_tab hash2file shift_first get_files_from_dir line2hash line2hash_comma slurp_excel slurp_hash_list complement);  
use 5.14.0;
use warnings;
use DirHandle;
use Carp qw(croak);
use Data::Dumper;
#use Spreadsheet::ParseExcel;

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

##Slurps but doesn't chomp
sub slurp_no_chomp {
    my ($file) = @_;
    my @data; 
    open IN, "<", $file or croak "can't open $file\n";
    @data = <IN>;
    close IN;
    return (@data);
}

##Slurp to lines of hashes
sub slurp_hash_list {
    my ($file) = @_;
    my @data   = slurp($file); 
    my $header = shift(@data);
    my @hash_list;
    
    for my $line (@data){
        my $h = line2hash($header,$line);
        push (@hash_list, $h);
    }
    
    return (@hash_list);
}


##Shifts the first line of a file
sub shift_first {
    my ($file) = @_;
    my @data = slurp($file);
    my $first = shift(@data);
    return ($first);
}


##Parses the lines of a file into a hash
##Each key gets a value of 1
sub file2hash {
    my ($file) = @_;
    my $data = {};
    my @lines = slurp($file);
    
    for my $line (@lines){
        chomp ($line);
        if ($line){
            $data->{$line}++;
        }
    }
    return ($data);
}


##Parses the data into a hash with the keys and values
## tab delimited. Skips lines starting with #
sub file2hash_tab {
    my ($file) = @_;
    my $data = {};
    my @lines = slurp($file);
    for my $line (@lines){
        chomp ($line);
        if ($line =~ /^\#/){
            next;
        }
        my ($key, $value) = split(/\t/, $line);
        
        if ($key){
            $data->{$key} = $value;
        }
    }
    return ($data);
}

##parses a hash into a file (tab delimited)
sub hash2file_tab {
    my ($hash, $file) = @_;
    open OUT, ">", $file;
    for my $key (keys %$hash){
        print OUT "$key\t$hash->{$key}\n";
    }
    close OUT;
}

##parses a hash into a file (only the keys)
sub hash2file {
    my ($hash, $file) = @_;
    open OUT, ">", $file;
    for my $key (keys %$hash){
        print OUT "$key\n";
    }
    close OUT;
}


##reads data from a file into a string (doesn't chomp)
sub read_file {
    my ($file) = @_;
    my $data; 
    open IN, "<", $file or croak "can't open $file\n";
    while (my $line = <IN>){
                $data .= $line;
    }
    return ($data);
    close IN;
}

##write data from a string into a file
sub write_file {
        my ($file, $data) = @_;
        open OUT, ">", $file or croak "can't open $file\n";
        print OUT "$_\n" foreach (@$data);
        close OUT;
}

##Gets files from a directory and returns them in list format
sub get_files_from_dir {
    my ($dir) = @_;
    my $dir_handle = new DirHandle $dir or die "unable to open $dir $!";
    my @dir_files;
    while (defined(my $file = $dir_handle->read)) {
        push (@dir_files, $file);
    }
    return @dir_files;
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

##Splits a line and parses it into the hash according to given headers (comma delimited)
sub line2hash_comma {
    my ($header, $line) = @_;

    ##Split header fields and line fields
    my @header_fields = split (/,/, $header);
    my @line_fields   = split (/,/, $line);

    ##Fill the hash
    my $hash = {};
    my $i = 0;
    for my $fld (@header_fields){
        $hash->{lc($fld)} = $line_fields[$i];
        $i++;
    }

    return $hash;
}

#========================================
##Always wanted to give this name to a sub
#========================================
sub complement {
    my $nuc = shift;
    given($nuc) {
        when('A') { return 'T'; }
        when('C') { return 'G'; }
        when('G') { return 'C'; }
        when('T') { return 'A'; }
        default { return; }
    }
}


##Returns an Excel file as a list of tab delimited lines
#======================================== 
#sub slurp_excel {
#======================================== 
    #my ($file) = @_;
    
    #my $parser   = Spreadsheet::ParseExcel->new();
    #my $workbook = $parser->parse($file);

    #if ( !defined $workbook ) {
        #die $parser->error(), ".\n";
    #}

    #my @lines;
    #for my $worksheet ( $workbook->worksheets() ) {
        #my ( $row_min, $row_max ) = $worksheet->row_range();
        #my ( $col_min, $col_max ) = $worksheet->col_range();

        #for my $row ( $row_min .. $row_max ) {
            #my @line;
            #for my $col ( $col_min .. $col_max ) {
                #my $cell = $worksheet->get_cell( $row, $col );
                #if (not defined $cell){
                    #push(@line, '');
                #}
                #next unless $cell;
                
                #my $val = $cell->value();
                #push (@line, $val);
            #}
            #my $line = join("\t", @line);
            #push (@lines, $line);
        #}
    #}
    #return (@lines);
#}


=pod

=head1 NAME

Utils::FileReader;
 Module for reading/writing to a file (and other cool stuff)

=head1 SYNOPSIS
    
    use lib $ENV{UNITY_ROOT} . '/lib/Miner';
    use Utils::FileReader qw(slurp line2hash read_file);
    
    my $input_file = 'input.txt';
    my @lines = slurp($input_file);
    for my $line (@lines){
        print "$line\n";
    }


=head1 SUBROUTINES

=over 8


=item read_file ( $file )

    Reads data from a file into a string (doesn't chomp)

=item slurp( $file )
    
    Parses the file and returns the lines as an array.
    Each line is chomped to remove "\n";

=item slurp_no_chomp( $file )
    
    Same as slurp but doesn't chomp lines

=item shift_first( $file )
    
    Shifts the first line of a file
    The file remains the same without the first line

=item file2hash ( $file )
    
    Parses the lines of a file into a hash
    Each key gets a value of 1

=item file2hash_tab ( $file)
    
    Parses the data into a hash with the keys and values
    tab delimited. Skips lines starting with #

=item hash2file_tab ($hash, $file)
    
    Parses a hash into a file (tab delimited)

=item sub get_files_from_dir ( $dir )
    
    Gets files from a directory and returns them in list format

=item line2hash ($header, $line)

    Splits a line and parses it into the hash according to given headers

=item slurp_excel ($file)

    Parses an Excel file and returns the data as a list of tab delimited lines

=item get_slice_adaptor ($reg_file, $species)

    Create a slice adaptor according to the species and registry

=item get_gene_adaptor ($reg_file, $species)

    Create a slice adaptor according to the species and registry



=back

=head1 AUTHOR

Guy Naamati





