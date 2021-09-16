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
  parse_isoform_FASTA_file download_FASTA_file 
  sort_isoforms_chr sort_clusters_by_position
  download_compara_TSV_file download_MAF_files parse_MAF_file
  perform_rest_action write_boxplot_file factorial fisher_yates_shuffle
  $REQUEST_COUNT $COMPARADIR $FASTADIR $MAFDIR @DIVISIONS
);

# Produces pan-gene analyses based on clusters of orthologous genes shared by
# species in a plant clade by querying pre-computed Compara protein trees 
# (orthology) and whole genome alignments (WGA, synteny, optional -W).
# Optionally, pan-genes might be named consistently across aligned genomes.
#
# Note: most species have only LASTZ genome alignments to a reference genome;
# some clades ie Oryza also have EPO multiple alignments, explained at
# http://plants.ensembl.org/info/genome/compara/multiple_genome_alignments.html
# Both types are used to compute WGAcoverage. 
#
# Copyright [2019-2021] EMBL-European Bioinformatics Institute

# Ensembl Genomes
my $RESTURL   = 'http://rest.ensembl.org';
my $INFOPOINT = $RESTURL . '/info/genomes/division/';
my $TAXOPOINT = $RESTURL . '/info/genomes/taxonomy/';

my $GZIPEXE   = 'gzip';
my $BEDTOOLSEXE = 'bedtools';
my $TRANSPOSEXE =
'perl -F\'\t\' -ane \'$F[$#F]=~s/\n//g;$r++;for(1 .. @F){$m[$r][$_]=$F[$_-1]};'
  . '$mx=@F;END{for(1 .. $mx){for $t(1 .. $r){print"$m[$t][$_]\t"}print"\n"}}\'';

# gene overlap in MAF alignment block
my $OVERLAP = 0.5;

# genome composition report
my $RNDSEED          = 12345;
my $NOFSAMPLESREPORT = 10;

my $downloadir = $Bin . '/downloads';
my $division   = 'Plants';
my $seqtype    = 'protein';
my $taxonid    = '';
# NCBI Taxonomy id, Brassicaceae=3700, Asterids=71274, Poaceae=4479

my $ref_genome = '';    # should be contained in $taxonid;
my ( $clusterdir, $comparadir, $fastadir) = ( '', '', '' );
my ( $beddir, $mafdir ) = ( '', '' );
my ($outfolder, $out_genome, $params) = ('', '', '');

my ( $help, $sp, $sp2, $show_supported, $request, $response );
my ( $filename, $dnafile, $pepfile, $seqfolder, $ext );
my ( $n_core_clusters, $n_cluster_sp, $n_cluster_seqs ) = ( 0, 0, 0 );
my ( $GOC, $WGA, $LOWCONF, $NOSINGLES , $GROWTH ) = ( 0, 0, 0, 0, 0 );
my ( $CHREGEX, $MAF ) = ( '', '' );
my ( $verbose, $bedtoolsexe ) = ( 0, $BEDTOOLSEXE );
my ( @ignore_species, @userTSVfiles, %ignore, %division_supported );

GetOptions(
    "help|?"        => \$help,
    "verbose|v"     => \$verbose,
    "supported|l"   => \$show_supported,
    "division|d=s"  => \$division,
    "clade|c=s"     => \$taxonid,
    "reference|r=s" => \$ref_genome,
    "outgroup|o=s"  => \$out_genome,
    "ignore|i=s"    => \@ignore_species,
    "user|u=s"      => \@userTSVfiles, 
    "type|t=s"      => \$seqtype,
    "GOC|G=i"       => \$GOC,
    "WGA|W=i"       => \$WGA,
    "LC|L"          => \$LOWCONF,
    "S|S"           => \$NOSINGLES,
	"growth|g"      => \$GROWTH,
	"MAF|M=s"       => \$MAF,     
    "folder|f=s"    => \$outfolder,
    "bedexe|exe=s"  => \$bedtoolsexe,
    "position|p=s"  => \$CHREGEX
) || help_message();

sub help_message {
    print "\nusage: $0 [options]\n\n"
      . "-c NCBI Taxonomy clade of interest         (required, example: -c Brassicaceae or -c 3700)\n"
      . "-f output folder                           (required, example: -f myfolder)\n"
      . "-r reference species_name to name clusters (required, example: -r arabidopsis_thaliana)\n"
      . "-l list supported species_names            (optional, example: -l)\n"
      . "-d Ensembl division                        (optional, default: -d $division)\n"
      . "-u parse species from user TSV files       (optional, as opposed E! species & Compara orthologues)\n"
      . "-o outgroup species_name                   (optional, example: -o brachypodium_distachyon)\n"
      . "-i ignore species_name(s)                  (optional, example: -i selaginella_moellendorffii -i ...)\n"

      # commented until I find a way of matching protein ids to transcript ids
      #. "-t sequence type [protein|cdna]            (optional, default: -t protein)\n"

      . "-g do pangene set growth simulation        (optional, produces [core|pan_gene]*.tab files)\n" 
      . "-L allow low-confidence orthologues        (optional, by default these are skipped)\n"
      . "-S skip singletons                         (optional, by default unclustered sequences are taken)\n"
      . '-p sort pangene clusters by chr position   (optional, requires regex to match chr names, example: -p \'^\d+$\';'
      . "\n                                            the example regular expression matches natural numbers,\n"
      . "                                            like the default chr names used in Ensembl Plants)\n"
      . "-b path to bedtools, useful with -M        (optional, if not in \$PATH, example: -b /path/to/bedtools)\n"
      . "-v verbose                                 (optional, example: -v\n";

    print "\nThe following options are only available for some clades:\n\n"
     #untested      
     #. "-M parse multiple genome alignments        (optional, requires bedtools, example: -M 8_rice.epo)\n\n"
      . "-G min Gene Order Conservation [0:100]  (optional, example: -G 75)\n"
      . "   see modules/Bio/EnsEMBL/Compara/PipeConfig/EBI/Plants/ProteinTrees_conf.pm\n"
      . "   at https://github.com/Ensembl/ensembl-compara\n\n"
      . "-W min Whole Genome Align score [0:100] (optional, example: -W 75)\n"
      . "   see ensembl-compara/scripts/pipeline/compara_plants.xml\n"
      . "   at https://github.com/Ensembl/ensembl-compara\n\n";
    print "Read about GOC and WGA at:\n"
      . "https://www.ensembl.org/info/genome/compara/Ortholog_qc_manual.html\n\n";

    print "Example calls:\n\n"
      . " perl $0 -c Brassicaceae -f Brassicaceae -o beta_vulgaris\n"
      . " perl $0 -f poaceae -c 4479 -r oryza_sativa -G 75 -g\n"
      . " perl $0 -f poaceae -c 4479 -r oryza_sativa -WGA 75 -M 8_rice.epo\n"
      . exit(0);
}

if ($help) { help_message() }

if ($division) {
    if ( !grep( /^$division$/, @PlantCompUtils::DIVISIONS ) ) {
        die "# ERROR: accepted values for division are: "
          . join( ',', @PlantCompUtils::DIVISIONS ) . "\n";
    }
    else {
        my $lcdiv = lc($division);

        $comparadir = $PlantCompUtils::COMPARADIR;
        $comparadir =~ s/xxx/$lcdiv/;

        $fastadir = $PlantCompUtils::FASTADIR;
        $fastadir =~ s/xxx/$lcdiv/;

        if(defined($MAF)){
            $mafdir = $PlantCompUtils::MAFDIR;
            $mafdir =~ s/xxx/$lcdiv/;
            $mafdir .= $MAF;

			# test bedtools is available
            my $bedversion = `$bedtoolsexe -version`;
            if($bedversion !~ "bedtools") {
                print "# ERROR: please set a valid -b /path/to/bedtools with option -M\n";
                exit;
            }
        }
    }
}

if ($show_supported) {
    print "# $0 -d $division -l \n\n";
}
else {

    if ( $ref_genome eq '' ) {
        print "# ERROR: need a valid reference species_name, ".
               "such as -r arabidopsis_thaliana)\n\n";
        exit;
    }
    else {
        $clusterdir = $ref_genome;
        $clusterdir =~ s/_//g;
        if ($out_genome) {
            $clusterdir .= "_plus_" . $out_genome;
        }
    }

    if ( $taxonid eq '' ) {
        print "# ERROR: need a valid NCBI Taxonomy clade, ".
               "such as -c Brassicaceae or -c 3700\n\n";
        print "# Check https://www.ncbi.nlm.nih.gov/taxonomy\n";
        exit;
    }
    else {
        $taxonid =~ s/\s+/%20/g;
        $clusterdir .= "_$taxonid\_algEnsemblCompara";
    }

    if ($CHREGEX) {
        $params .= "_chrpos";
    }

    if ($GOC) {
        $params .= "_GOC$GOC";
    }

    if ($WGA) {
        $params .= "_WGA$WGA";
    }

    if ($LOWCONF) {
        $params .= "_LC";
    }

    if ($NOSINGLES) {
        $params .= "_nosingles";
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
            die "# ERROR: currently cannot accept seqtype = cdna\n";
            # TODO: the trick would be to download both pep & cdna and pair protein stable_ids 
            # (used in Compara) to transcript ids from the FASTA header in pep
        }
    }

    if ($outfolder) {
        if ( -e $outfolder ) {
            print "\n# WARNING : folder '$outfolder' exists, files might be overwritten\n\n";
        }
        else {
            if ( !mkdir($outfolder) ) {
                die "# ERROR: cannot create $outfolder\n";
            }
        }

        # create $clusterdir with $params
        $clusterdir .= $params;
        if ( !-e "$outfolder/$clusterdir" ) {
            if ( !mkdir("$outfolder/$clusterdir") ) {
                die "# ERROR: cannot create $outfolder/$clusterdir\n";
            }
        }

        # create $beddir 
        $beddir = "$outfolder/bed";
        if ( defined($MAF) && !-e $beddir ) {
            if ( !mkdir("$beddir") ) {
                die "# ERROR: cannot create $beddir\n";
            }
        }
    }
    else {
        print
          "# ERROR: need a valid output folder, such as -f Brassicaceae\n\n";
        exit;
    }

    print "# $0 -d $division -c $taxonid -r $ref_genome -o $out_genome "
      . "-f $outfolder -t $seqtype -b $bedtoolsexe -M $MAF "
      . "-p '$CHREGEX' -G $GOC -W $WGA "
      . "-g $GROWTH -L $LOWCONF -S $NOSINGLES -v $verbose\n\n";
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

my ( $n_of_species, $cluster_id, $chr ) = ( 0, '' );
my ( @supported_species, @cluster_ids, %sorted_cluster_ids );
my ( %supported, %incluster, %cluster, %compara_isoform );
my ( %sequence, %header, %bedfiles );
my ( %totalgenes, %totalclusters, %POCP_matrix );
my ( %MAFblocks, %isoform2block, %sorted_ids, %id2chr );

$request = $TAXOPOINT . "$taxonid?";

$response = perform_rest_action( $http, $request, $global_headers );
$infodump = decode_json($response);

foreach $sp ( @{$infodump} ) {
    if ( $sp->{'name'} && $division_supported{ $sp->{'name'} } ) {

        next if ( $ignore{ $sp->{'name'} } );

        # add sorted clade species except reference
        $supported{ $sp->{'name'} } = 1;
        if ( $sp->{'name'} ne $ref_genome ) {
            push( @supported_species, $sp->{'name'} );
        }
    }
}

# check reference genome is supported
if ( !$supported{$ref_genome} ) {
    die "# ERROR: cannot find $ref_genome within NCBI taxon $taxonid\n";
}
else {
    # ref genome is first in array
    unshift( @supported_species, $ref_genome );

    if ($verbose) {
        foreach $sp (@supported_species) {
            print "# $sp\n";
        }
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

## 2) get orthologous (plant) genes shared by selected species

# columns of TSV file
# NOTE: high_conf are not be available for some divisions
my (
    $gene_stable_id,     $prot_stable_id, $species,
    $identity,           $homology_type,  $hom_gene_stable_id,
    $hom_prot_stable_id, $hom_species,    $hom_identity,
    $dn,                 $ds,             $goc_score,
    $wga_coverage,       $high_confidence
);

# Iteratively get and parse TSV files that define pairs of orthologues
# sequences computed by Ensembl Compara. After parsing all pairwise
# TSV files clusters emerge.
# Note: each sequence is identified with a protein stable_id that
# corresponds to the canonical isoform of that gene.
# Read more at
# http://plants.ensembl.org/info/genome/compara/peptide_compara.html
# http://plants.ensembl.org/info/website/glossary.html
foreach $sp (@supported_species) {

    # get TSV file; these files are bulky and might take some time to download
    my $stored_compara_file =
      download_compara_TSV_file( $comparadir, $sp, $downloadir );

    # uncompress on the fly and parse
    open( TSV, "$GZIPEXE -dc $stored_compara_file |" )
      || die "# ERROR: cannot open $stored_compara_file\n";
    while ( my $line = <TSV> ) {

        #ATMG00030 ATMG00030.1 arabidopsis_thaliana 52.3364 ortholog_one2many \
        #Tp57577 Tp57577 trifolium_pratense 16.8675 NULL NULL NULL NULL 0
        (
            $gene_stable_id,     $prot_stable_id, $species,
            $identity,           $homology_type,  $hom_gene_stable_id,
            $hom_prot_stable_id, $hom_species,    $hom_identity,
            $dn,                 $ds,             $goc_score,
            $wga_coverage,       $high_confidence
        ) = split( /\t/, $line );

        next if ( !$supported{$species} || !$supported{$hom_species} );

        if ( defined($high_confidence) ) {
            if ( $LOWCONF == 0
                && ( $high_confidence eq 'NULL' || $high_confidence == 0 ) )
            {
                #print 
                #  "# skip $prot_stable_id,$hom_prot_stable_id due to low-confidence\n"
                #  if ($verbose);
                next;
            }
        }

        next if ( $WGA && ( $wga_coverage eq 'NULL' || $wga_coverage < $WGA ) );

        next if ( $GOC && ( $goc_score eq 'NULL' || $goc_score < $GOC ) );

        if ( $homology_type =~ m/ortholog/ ) {

            # add $species protein to cluster only if not clustered yet
            if ( !$incluster{$prot_stable_id} ) {

                if ( $incluster{$hom_prot_stable_id} ) {

                    # use existing cluster_id from other species ortholog
                    $cluster_id = $incluster{$hom_prot_stable_id};
                }
                else {

                    # otherwise create a new one
                    $cluster_id = $prot_stable_id;
                    push( @cluster_ids, $cluster_id );
                }

                # record to which cluster this protein belongs
                $incluster{$prot_stable_id} = $cluster_id;

                push( @{ $cluster{$cluster_id}{$species} }, $prot_stable_id );

            }
            else {
                # set cluster for $hom_species anyway
                $cluster_id = $incluster{$prot_stable_id};
            } 

            # now add $hom_species protein to previously defined cluster
            if ( !$incluster{$hom_prot_stable_id} ) {

                # record to which cluster this protein belongs
                $incluster{$hom_prot_stable_id} = $cluster_id;

                push(
                    @{ $cluster{$cluster_id}{$hom_species} },
                    $hom_prot_stable_id
                );
            }

            # save isoforms used in compara
            $compara_isoform{$prot_stable_id} = 1;
            $compara_isoform{$hom_prot_stable_id} = 1;
        }
    }
    close(TSV); 
} 

# count how many clusters include each species
foreach $cluster_id (@cluster_ids) {
    foreach $species (@supported_species) {
        if ( $cluster{$cluster_id}{$species} ) {
            $totalclusters{$species}++;
        }
    }
}

# Get and parse FASTA files to get sequences & headers 
# of isoforms in the Compara clusters
# Note: uses %compara_isoform, created previously
foreach $sp (@supported_species) {

    my $stored_sequence_file =
        download_FASTA_file( $fastadir, "$sp/$seqfolder", $downloadir );

    my ( $ref_sequence, $ref_header ) =
       parse_isoform_FASTA_file( $stored_sequence_file, \%compara_isoform );

    my ($ref_bedfiles, $ref_sorted_ids, $ref_id2chr) =
       sort_isoforms_chr($ref_header, $beddir, $sp, $CHREGEX);

    printf("# wrote sorted isoforms of $sp in %d BED files (1-based)\n\n",
       scalar(keys(%$ref_bedfiles)));

    $bedfiles{$sp} = $ref_bedfiles;
    $sorted_ids{$sp} = $ref_sorted_ids;
    $id2chr{$sp} = $ref_id2chr;

    # count number of genes/selected isoforms in this species
    $totalgenes{$sp} = scalar( keys(%$ref_sequence) );

    # save these sequences
    foreach $prot_stable_id ( keys(%$ref_sequence) ) {
        $sequence{$sp}{$prot_stable_id} = $ref_sequence->{$prot_stable_id};
        $header{$sp}{$prot_stable_id} = $ref_header->{$prot_stable_id};
    } 
}

# Note: even using -W clusters might contain sequences from different chr
# We might want to split or duplicate these clusters
# Example: Os07t0248800-02,Os10t0102900-00 (oryza_sativa)	BGIOSGA024529-PA (oryza_indica)
# http://plants.ensembl.org/Oryza_sativa/Location/Multi?db=core;g=Os10g0102900;g1=BGIOSGA024529;r=10:227687-234770;s1=Oryza_indica;t=Os10t0102900-00;r1=7:8306227-8312670:1;time=1620306359
# http://plants.ensembl.org/Oryza_sativa/Location/Multi?db=core;g=Os07g0248800;g1=BGIOSGA024529;r=7:8264635-8271779;s1=Oryza_indica;r1=7:8306227-8312670:1;time=1620306355

# add unclustered sequences as singletons
my $total_seqs = 0;
foreach $sp (@supported_species) {

    my $singletons = 0;

    foreach $prot_stable_id ( sort keys( %{ $sequence{$sp} } ) ) {

        next if ( $NOSINGLES || $incluster{$prot_stable_id} );    # skip

        # create new cluster
        $cluster_id = $prot_stable_id;
        $incluster{$prot_stable_id} = $cluster_id;

        push( @{ $cluster{$cluster_id}{$sp} }, $prot_stable_id );
        push( @cluster_ids,                    $cluster_id );

        # add this singleton to total clusters
        $totalclusters{$sp}++;

        $singletons++;
    }

    $total_seqs += $totalgenes{$sp};

    if($CHREGEX) {
        printf( "# %s : sequences = %d clusters = %d (singletons = %d, placed = %d)\n",
            $sp, $totalgenes{$sp}, $totalclusters{$sp}, $singletons,  
            scalar(keys(%{ $id2chr{$sp} })));
    } else {
        printf( "# %s : sequences = %d clusters = %d (singletons = %d)\n",
            $sp, $totalgenes{$sp}, $totalclusters{$sp}, $singletons );
    }
}

printf( "\n# total sequences = %d\n\n", $total_seqs );

# NOTE: not used so far
# Ideas: % genome in blocks, <genes> / block, detection of overlapping non-orths
# download multiple alignment files in MAF format if required,
# parse them and get stats
if($MAF ne '') {

    my ($n_of_blocks, $maf_file, @stored_maf_files) = (0);

    print "# parsing multiple alignment MAF files\n";

	@stored_maf_files = download_MAF_files( $mafdir, $MAF, $downloadir, $verbose );
    printf( "\n# total MAF files = %d\n", scalar(@stored_maf_files) );

    #my $tmpi = 0;#debug
    foreach $maf_file (@stored_maf_files) {

        parse_MAF_file( $maf_file, \@supported_species,
            \%bedfiles, $bedtoolsexe, $OVERLAP,
            \%MAFblocks, \%isoform2block, $verbose);

		#last if($tmpi++ > 2);
    }

	printf("# total aligned blocks = %d\n", scalar(keys(%MAFblocks)));
}

## 3) write sequence clusters, summary text file and POCP matrix

# POCP=Percent Conserved Sequences (POCP) matrix
my $POCP_matrix_file = "$outfolder/POCP.matrix$params\.tab";

my $cluster_summary_file = "$outfolder/$clusterdir.cluster_list";

open( CLUSTER_LIST, ">", $cluster_summary_file )
  || die "# ERROR: cannot create $cluster_summary_file\n";

$n_core_clusters = 0;

foreach $cluster_id (@cluster_ids) {

    if ( scalar( keys( %{ $cluster{$cluster_id} } ) ) == $n_of_species ) {
        $n_core_clusters++;
    }

    # sequence cluster
    $n_cluster_sp = $n_cluster_seqs = 0;
    $filename = $cluster_id;

    # for summary, in case this was run twice (cdna & prot)
    $dnafile = $filename . '.fna';
    $pepfile = $filename . '.faa';

    # write sequences and count sequences
    my ( %cluster_stats, @cluster_species );
    open( CLUSTER, ">", "$outfolder/$clusterdir/$filename$ext" )
      || die "# ERROR: cannot create $outfolder/$clusterdir/$filename$ext\n";

	foreach $species (@supported_species) {
        next if ( !$cluster{$cluster_id}{$species} );
        $n_cluster_sp++;
        foreach $prot_stable_id ( @{ $cluster{$cluster_id}{$species} } ) {
            printf(CLUSTER ">%s [%s] %s\n", $prot_stable_id, $species, 
                $header{$species}{$prot_stable_id} || "no_header");
            if ( $sequence{$species}{$prot_stable_id} ) {
                print CLUSTER "$sequence{$species}{$prot_stable_id}\n";
            }
            else {
                print "# cannot find peptide $prot_stable_id ($species)\n"
                  if ($verbose);
            }

            $n_cluster_seqs++;
            $cluster_stats{$species}++;
        }
    }
    close(CLUSTER);

    # cluster summary
    @cluster_species = keys(%cluster_stats);
    if ( !-s "$outfolder/$clusterdir/$dnafile" ) { $dnafile = 'void' }
    if ( !-s "$outfolder/$clusterdir/$pepfile" ) { $pepfile = 'void' }

    print CLUSTER_LIST
         "cluster $cluster_id size=$n_cluster_seqs taxa=$n_cluster_sp ".
          "file: $dnafile aminofile: $pepfile\n";

    foreach $species (@cluster_species) {
        foreach $prot_stable_id ( @{ $cluster{$cluster_id}{$species} } ) {
            print CLUSTER_LIST ": $species\n";
        }
    }

    # update PCOP data
    foreach $sp ( 0 .. $#cluster_species - 1 ) {
        foreach $sp2 ( $sp + 1 .. $#cluster_species ) {

            # add the number of sequences in this cluster from a pair of species/taxa
            $POCP_matrix{ $cluster_species[$sp] }{ $cluster_species[$sp2] } +=
              $cluster_stats{ $cluster_species[$sp] };
            $POCP_matrix{ $cluster_species[$sp] }{ $cluster_species[$sp2] } +=
              $cluster_stats{ $cluster_species[$sp2] };

            # now in reverse order to make sure it all adds up
            $POCP_matrix{ $cluster_species[$sp2] }{ $cluster_species[$sp] } +=
              $cluster_stats{ $cluster_species[$sp] };
            $POCP_matrix{ $cluster_species[$sp2] }{ $cluster_species[$sp] } +=
              $cluster_stats{ $cluster_species[$sp2] };
        }
    }
}

close(CLUSTER_LIST);

printf( "\n# number_of_clusters = %d (core = %d)\n\n",
    scalar(@cluster_ids), $n_core_clusters );
print "# cluster_list = $outfolder/$clusterdir.cluster_list\n";
print "# cluster_directory = $outfolder/$clusterdir\n";

# print POCP matrix
open( POCPMATRIX, ">$POCP_matrix_file" )
  || die "# EXIT: cannot create $POCP_matrix_file\n";

print POCPMATRIX "genomes";
foreach $sp ( 0 .. $#supported_species ) {
    print POCPMATRIX "\t$supported_species[$sp]";
}
print POCPMATRIX "\n";

my (%POCP2ref,$perc);

foreach $sp ( 0 .. $#supported_species ) {
    print POCPMATRIX "$supported_species[$sp]";
    foreach $sp2 ( 0 .. $#supported_species ) {

        if ( $sp == $sp2 ) { 
            print POCPMATRIX "\t100.00"
        }
        else {
            if ( $POCP_matrix{ $supported_species[$sp] }
                { $supported_species[$sp2] } )
            {
                $perc = sprintf("\t%1.2f",
                    (
                        100 * $POCP_matrix{ $supported_species[$sp] }
                          { $supported_species[$sp2] }
                    ) / (
                        $totalgenes{ $supported_species[$sp] } +
                          $totalgenes{ $supported_species[$sp2] }
                    )
                );
                print POCPMATRIX "\t$perc";

                # save %POCP for all species vs reference
                if($sp == 0){ $POCP2ref{$supported_species[$sp2]} = $perc }
            }
            else {
                print POCPMATRIX "\tNA";
            }
        }
    }
    print POCPMATRIX "\n";
}
close(POCPMATRIX);

print "\n# percent_conserved_proteins_file = $POCP_matrix_file\n\n";

# sort species from ref down by decreasing POCP
my @supported_species_POCP;

# reference goes in 1st place
push(@supported_species_POCP, $ref_genome);

foreach $sp2 (sort {$POCP2ref{$b}<=>$POCP2ref{$a}} keys(%POCP2ref)) {
    push(@supported_species_POCP, $sp2);
}


## 4)  write pangenome matrices in output folder

## if required sort clusters following gene order of i) ref species 
## and ii) other supported species sorted by shared from close to distant
if(!$CHREGEX){
	push(@{ $sorted_cluster_ids{'unsorted'} }, @cluster_ids );
}
else {
    %sorted_cluster_ids = sort_clusters_by_position( 
        \@supported_species_POCP, \%sorted_ids, \%incluster, 
        \%cluster, \%id2chr, $CHREGEX );

    foreach $chr (sort keys(%sorted_cluster_ids)) {
	    printf("# clusters sorted by position in chr %s = %d\n", 
        $chr, scalar(@{ $sorted_cluster_ids{$chr} }));
    }
}

# set matrix filenames and write headers
my $pangenome_matrix_file = "$outfolder/pangenome_matrix$params\.tab";
my $pangenome_gene_file   = "$outfolder/pangenome_matrix_genes$params\.tab";
my $pangenome_matrix_tr   = "$outfolder/pangenome_matrix$params\.tr.tab";
my $pangenome_gene_tr     = "$outfolder/pangenome_matrix_genes$params\.tr.tab";
my $pangenome_fasta_file  = "$outfolder/pangenome_matrix$params\.fasta";

open( PANGEMATRIX, ">$pangenome_matrix_file" )
  || die "# EXIT: cannot create $pangenome_matrix_file\n";

open( PANGENEMATRIX, ">$pangenome_gene_file" )
  || die "# EXIT: cannot create $pangenome_gene_file\n";

print PANGEMATRIX "source:$outfolder/$clusterdir";
foreach $chr (keys(%sorted_cluster_ids)) {
    print PANGEMATRIX "\tchr$chr";
    foreach $cluster_id (@{ $sorted_cluster_ids{$chr} }) {
        print PANGEMATRIX "\t$cluster_id$ext"; 
    }
}	
print PANGEMATRIX "\n";

print PANGENEMATRIX "source:$outfolder/$clusterdir";
foreach $chr (keys(%sorted_cluster_ids)) {
    print PANGENEMATRIX "\tchr$chr";
    foreach $cluster_id (@{ $sorted_cluster_ids{$chr} }) {
        print PANGENEMATRIX "\t$cluster_id$ext"; 
    }
}	
print PANGENEMATRIX "\n";

open( PANGEMATRIF, ">$pangenome_fasta_file" )
  || die "# EXIT: cannot create $pangenome_fasta_file\n";

foreach $species (@supported_species_POCP) {

    print PANGEMATRIX "$species";
    print PANGENEMATRIX "$species";
    print PANGEMATRIF ">$species\n";

    foreach $chr (keys(%sorted_cluster_ids)) {

        # chr lines have no genes
        print PANGEMATRIX "\tNA";
        print PANGENEMATRIX "\tNA";

        foreach $cluster_id (@{ $sorted_cluster_ids{$chr} }) {

            if ( $cluster{$cluster_id}{$species} ) {
                printf( PANGEMATRIX "\t%d",
                    scalar( @{ $cluster{$cluster_id}{$species} } )
                );
                printf( PANGENEMATRIX "\t%s",
                    join( ',', @{ $cluster{$cluster_id}{$species} } )
                );
                print PANGEMATRIF "1";
            }
            else {    # absent genes
                print PANGEMATRIX "\t0";
                print PANGENEMATRIX "\t-";
                print PANGEMATRIF "0";
            }
        }
    }

    print PANGEMATRIX "\n";
    print PANGENEMATRIX "\n";
    print PANGEMATRIF "\n";
}

close(PANGEMATRIX);
close(PANGENEMATRIX);
close(PANGEMATRIF);

system("$TRANSPOSEXE $pangenome_matrix_file > $pangenome_matrix_tr");
system("$TRANSPOSEXE $pangenome_gene_file > $pangenome_gene_tr");

print
"# pangenome_file = $pangenome_matrix_file tranposed = $pangenome_matrix_tr\n";
print
  "# pangenome_genes = $pangenome_gene_file transposed = $pangenome_gene_tr\n";
print "# pangenome_FASTA_file = $pangenome_fasta_file\n";


## TODO: haplotypes per cluster


exit if(!$GROWTH);

## 5) optionally make genome composition analysis to simulate pangene growth
## NOTE: this is measured in clusters added/missed per genome

my ( $core_occup, $mean, $sd, $data_file, $sort, $s ); #$s = sample
my ( %previous_sorts, @sample, @clusters, @pangenome, @coregenome );
my @taxa    = @supported_species;
my @tmptaxa = @taxa;

my $n_of_permutations = sprintf( "%g", factorial($n_of_species) );
if ( $n_of_permutations < $NOFSAMPLESREPORT ) {
    $NOFSAMPLESREPORT = $n_of_permutations;
}
printf( "\n# genome composition report (samples=%d,seed=%d)\n",
    $NOFSAMPLESREPORT, $RNDSEED );

# random-sort the list of taxa $NOFSAMPLESREPORT times
for ( $s = 0 ; $s < $NOFSAMPLESREPORT ; $s++ ) {
    if ( $s > 0 ) {    # reshuffle until a new permutation is obtained
        $sort = fisher_yates_shuffle( \@tmptaxa );
        while ( $previous_sorts{$sort} ) {
            $sort = fisher_yates_shuffle( \@tmptaxa );
        }
        $previous_sorts{$sort} = 1;
    }
    push( @{ $sample[$s] }, @tmptaxa );
}

# sample taxa in random order
for ( $s = 0 ; $s < $NOFSAMPLESREPORT ; $s++ ) {
    my ( %n_of_taxa_in_cluster, $sample );
    @tmptaxa = @{ $sample[$s] };

    $sample = "## sample $s ($tmptaxa[0] | ";
    for ( $sp = 0 ; $sp < $n_of_species ; $sp++ ) {
        $sp2 = 0;
        while ( $tmptaxa[$sp] ne $taxa[$sp2] ) { $sp2++ }
        $sample .= "$sp2,";
        if ( length($sample) > 70 ) {
            $sample .= '...';    # trim it
            last;
        }
    }
    $sample .= ')';
    print "$sample\n";

    # calculate pan/core-gene size adding genomes one-by-one
    $coregenome[$s][0] = $totalclusters{ $tmptaxa[0] };
    $pangenome[$s][0]  = $coregenome[$s][0];
    print
      "# adding $tmptaxa[0]: core=$coregenome[$s][0] pan=$pangenome[$s][0]\n"
      if ($verbose);

    for ( $sp = 1 ; $sp < $n_of_species ; $sp++ ) {
        $coregenome[$s][$sp] = 0;
        $pangenome[$s][$sp]  = $pangenome[$s][ $sp - 1 ];
        $core_occup          = $sp + 1;

        foreach $chr (keys(%sorted_cluster_ids)) {
	        foreach $cluster_id (@{ $sorted_cluster_ids{$chr} }) {

                # check reference species is in this cluster (1st iteration only)
                if ( $sp == 1 && $cluster{$cluster_id}{ $tmptaxa[0] } ) {
                     $n_of_taxa_in_cluster{$cluster_id}++;
                }

                # check $sp is in this cluster
                if ( $cluster{$cluster_id}{ $tmptaxa[$sp] } ) {
                    $n_of_taxa_in_cluster{$cluster_id}++;
                }

                # check cluster occupancy
                if (   $n_of_taxa_in_cluster{$cluster_id}
                    && $cluster{$cluster_id}{ $tmptaxa[$sp] } )
                {

                    # core genes must contain all previously seen species
                    if ( $n_of_taxa_in_cluster{$cluster_id} == $core_occup ) {
                        $coregenome[$s][$sp]++;

                    }    # pan genes must be novel to this species
                    elsif ( $n_of_taxa_in_cluster{$cluster_id} == 1 ) {
                        $pangenome[$s][$sp]++;
                    }
                }
            }
        }

        print "# adding $tmptaxa[$sp]: core=$coregenome[$s][$sp] pan=$pangenome[$s][$sp]\n"
          if ($verbose);
    }
}

# write genome composition stats to boxplot files
my $pan_file  = "$outfolder/pan_gene$params\.tab";
my $core_file = "$outfolder/core_gene$params\.tab";

write_boxplot_file( $pan_file, $n_of_species, $NOFSAMPLESREPORT, \@pangenome );
print "\n# pan-gene (number of clusters) = $pan_file\n";

write_boxplot_file( $core_file, $n_of_species, $NOFSAMPLESREPORT,
    \@coregenome );
print "# core-gene (number of clusters) = $core_file\n";

my $end_time = new Benchmark();
print "\n# runtime: "
  . timestr( timediff( $end_time, $start_time ), 'all' ) . "\n";
