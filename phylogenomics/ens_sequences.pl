#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);
use Benchmark;
use HTTP::Tiny;
use JSON qw(decode_json);
use FindBin '$Bin';
use lib $Bin;
use PlantCompUtils qw(
  list_ensembl_mysql_dbs get_canonical_transcript_ids
  download_FASTA_file parse_isoform_FASTA_file perform_rest_action
  $REQUEST_COUNT $FASTADIR @DIVISIONS
);

# Downloads cds/pep sequences of species in a clade from Ensembl Plants.
# Uses canonical transcripts, used in the gene tree analysis,
# which usually are the longest translation with no stop codons
#
# Copyright [2019-2021] EMBL-European Bioinformatics Institute

# Ensembl Genomes
my $RESTURL   = 'http://rest.ensembl.org';
my $INFOPOINT = $RESTURL . '/info/genomes/division/';
my $TAXOPOINT = $RESTURL . '/info/genomes/taxonomy/';

my $downloadir = $Bin . '/downloads';
my $division   = 'Plants';
my $seqtype    = 'protein';
my $taxonid    = '';    # NCBI Taxonomy id Brassicaceae=3700
                        # Asterids=71274, Poaceae=4479

my ( $fastadir, $outfile, $out_genome ) = ( '', '', '' );

my ( $help, $sp, $line, $id, $show_supported, $request, $response );
my ( $filename, $dnafile, $pepfile, $seqfolder, $ext );
my ( @ignore_species, %ignore, %division_supported );

GetOptions(
    "help|?"       => \$help,
    "supported|l"  => \$show_supported,
    "division|d=s" => \$division,
    "clade|c=s"    => \$taxonid,
    "outgroup|o=s" => \$out_genome,
    "ignore|i=s"   => \@ignore_species,
    "type|t=s"     => \$seqtype,
    "outfile|f=s"  => \$outfile
) || help_message();

sub help_message {
    print "\nusage: $0 [options]\n\n"
      . "-c NCBI Taxonomy clade of interest (required, example: -c Brassicaceae or -c 3700)\n"
      . "-f output FASTA file               (required, example: -f myfile.fasta)\n"
      . "-l list supported species_names    (optional, example: -l)\n"
      . "-d Ensembl division                (optional, default: -d $division)\n"
      . "-o outgroup species_name           (optional, example: -o brachypodium_distachyon)\n"
      . "-i ignore species_name(s)          (optional, example: -i selaginella_moellendorffii -i ...)\n"
      . "-t sequence type [protein|cdna]    (optional, default: -t protein)\n\n";

    print "Example calls:\n\n"
      . " perl $0 -c Liliopsida -o arabidopsis_thaliana -f Liliopsoda.Atha.EG44.faa\n";
    exit(0);
}

if ($help) { help_message() }

if ($division) {
    if ( !grep( /^$division$/, @PlantCompUtils::DIVISIONS ) ) {
        die "# ERROR: accepted values for division are: "
          . join( ',', @PlantCompUtils::DIVISIONS ) . "\n";
    }
    else {
        my $lcdiv = lc($division);
        $fastadir = $PlantCompUtils::FASTADIR;
        $fastadir =~ s/xxx/$lcdiv/;
    }
}

if ($show_supported) {
    print "# $0 -d $division -l \n\n";
}
else {

    if ( $taxonid eq '' ) {
        print "# ERROR: need a valid NCBI Taxonomy clade, ".
               "such as -c Brassicaceae or -c 3700\n\n";
        print "# Check https://www.ncbi.nlm.nih.gov/taxonomy\n";
        exit;
    }
    else {
        $taxonid =~ s/\s+/%20/g;
    }

    if (@ignore_species) {
        foreach my $sp (@ignore_species) {
            $ignore{$sp} = 1;
        }
        printf( "\n# ignored species : %d\n\n", scalar( keys(%ignore) ) );
    }

    if ( $seqtype ne 'protein' && $seqtype ne 'cdna' ) {
        die "# ERROR: accepted values for seqtype are: protein|cdna\n";
    }
    else {
        if ( $seqtype eq 'protein' ) {
            $ext       = '.faa';
            $seqfolder = 'pep';
        }
        else {
            $ext       = '.fna';
            $seqfolder = 'cdna';
        }
    }

    if ( !$outfile ) {
        print "# ERROR: need a valid output file, such as -f Brassicaceae.fasta\n\n";
        exit;
    }

    print "# $0 -d $division -c $taxonid -o $out_genome -f $outfile -t $seqtype\n\n";
}

my $start_time = new Benchmark();

# new object and params for REST requests
my $http = HTTP::Tiny->new();
my $global_headers = { 'Content-Type' => 'application/json' };
$PlantCompUtils::REQUEST_COUNT = 0;

## 0) check supported species in division 

$request = $INFOPOINT . "Ensembl$division?";

$response = perform_rest_action( $http, $request, $global_headers );
my $infodump = decode_json($response);

foreach $sp ( @{$infodump} ) {
    if ( $sp->{'has_peptide_compara'} ) {
        $division_supported{ $sp->{'name'} } = 1;
    }
}

# list supported species and exit
if ($show_supported) {

    foreach $sp ( sort( keys(%division_supported) ) ) {
        print "$sp\n";
    }
    exit;
}

# check outgroup is supported
if ( $out_genome && !$division_supported{$out_genome} ) {
    die "# ERROR: genome $out_genome is not supported\n";
}

## 1) check species in clade

my ( $n_of_species, $n_of_sequences ) = ( 0, 0 );
my ( @supported_species, %supported );

$request = $TAXOPOINT . "$taxonid?";

$response = perform_rest_action( $http, $request, $global_headers );
$infodump = decode_json($response);

foreach $sp ( @{$infodump} ) {
    if ( $sp->{'name'} && $division_supported{ $sp->{'name'} } ) {

        next if ( $ignore{ $sp->{'name'} } );

        # add sorted clade species
        $supported{ $sp->{'name'} } = 1;
        push( @supported_species, $sp->{'name'} );
    }
}

printf( "# supported species in NCBI taxon %s : %d\n\n",
    $taxonid, scalar(@supported_species) );

# add outgroup if required
if ($out_genome) {
    push( @supported_species, $out_genome );
    $supported{$out_genome} = 1;
    print "# outgenome: $out_genome\n";
}

$n_of_species = scalar(@supported_species);
print "# total selected species : $n_of_species\n\n";

## 2) connect to public Ensembl server and 
##    find latest database schema for each species

my $ref_dbs = list_ensembl_mysql_dbs();

my %species2db;
foreach $sp (@supported_species) {
  foreach my $db (@$ref_dbs) {
    if($db =~ /$sp\_core_\d+/) {
      $species2db{$sp} = $db;
    } 
  }
}

## 3) get sequences for selected (plant) species

open( OUTFILE, ">", $outfile ) || die "# ERROR: cannot create $outfile\n";

# iteratively get and parse FASTA files
foreach $sp (@supported_species) {

    # get list of canonical transcripts for this species
    my $ref_canon_isofs = get_canonical_transcript_ids($species2db{$sp});

    printf("# %s canonical isoforms=%d\n", $sp,
      scalar(keys(%$ref_canon_isofs)));

    # now get FASTA file and parse it, selected/longest isoforms are read
    my $stored_sequence_file =
      download_FASTA_file( $fastadir, "$sp/$seqfolder", $downloadir );

    my ( $ref_sequence, $ref_header ) =
      parse_isoform_FASTA_file($stored_sequence_file, $ref_canon_isofs);

    foreach $id ( keys(%$ref_sequence) ) {
        print OUTFILE ">$id $ref_header->{$id} [$sp]\n$ref_sequence->{$id}\n";
        $n_of_sequences++;
    }
}

close(OUTFILE);

print "# created $outfile with $n_of_sequences sequences\n";
