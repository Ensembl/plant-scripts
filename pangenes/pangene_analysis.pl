#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);

# Produces pan-gene analysis based on clusters of collinear genes shared by
# species in a pre-computed Minimap2/Wfmash synteny TSV file, produced by 
# get_collinear_genes.pl
#
# Copyright [2021] EMBL-European Bioinformatics Institute

my $TRANSPOSEXE =
'perl -F\'\t\' -ane \'$F[$#F]=~s/\n//g;$r++;for(1 .. @F){$m[$r][$_]=$F[$_-1]};'
  . '$mx=@F;END{for(1 .. $mx){for $t(1 .. $r){print"$m[$t][$_]\t"}print"\n"}}\'';

# these globals match sequence type to extension
my @SEQTYPE = qw( cds pep cdna ); # Note cdna goes last
my %SEQEXT = (
    'cdna' => '.cdna.fna',
    'cds' => '.cds.fna',
    'pep' => '.cds.faa'
);

# genome composition report
my $RNDSEED          = 12345;
my $NOFSAMPLESREPORT = 10;

my ( $ref_genome, $seqfolder, $clusterdir ) = ( '', '', '' );
my ( $outfolder, $params) = ('', '');
my ( $help, $sp, $sp2, $show_supported );
my ( $infile, $filename, $cdsfile, $pepfile );
my ( $n_core_clusters, $n_cluster_sp, $n_cluster_seqs ) = ( 0, 0, 0 );
my ( $NOSINGLES , $GROWTH ) = ( 0, 0 );
my ( $n_of_species, $verbose ) = ( 0, 0 );
my ( @infiles, @supported_species, @ignore_species);
my ( %species, %ignore, %supported );

GetOptions(
    "help|?"        => \$help,
    "verbose|v"     => \$verbose,
    "supported|l"   => \$show_supported,
    "TSV|T=s"       => \@infiles,
    "reference|r=s" => \$ref_genome,
    "ignore|i=s"    => \@ignore_species,
    "S|S"           => \$NOSINGLES,
	"growth|g"      => \$GROWTH,
    "folder|f=s"    => \$outfolder,
	"seq|s=s"       => \$seqfolder,
) || help_message();

sub help_message {
    print "\nusage: $0 [options]\n\n"
      . "-T input collinear TSV file(s)             (required, example: -T Minimap2.homologies.rice.overlap0.5.tsv -T ...)\n"
      . "-f output folder                           (required, example: -f myfolder)\n"
      . "-r reference species_name to name clusters (required, example: -r arabidopsis_thaliana)\n"
      . "-l list supported species in -T file       (optional, example: -l)\n"
      . "-i ignore species_name(s)                  (optional, example: -i selaginella_moellendorffii -i ...)\n"
      . "-g do pangene set growth simulation        (optional, produces [core|pan_gene]*.tab files)\n" 
      . "-S skip singletons                         (optional, by default unclustered sequences are taken)\n"
      . "-s folder with gene seqs of species in TSV (optional, default: \$PWD)\n"
      . "-v verbose                                 (optional, example: -v\n";

    exit(0);
}

if ($help) { help_message() }

if ($show_supported) {
    print "# $0 -l \n\n";
}
else {

    if ( scalar(@infiles) < 1 ) {
        print "# ERROR: need at least one valid input TSV file\n\n";
        exit;
    }
    else {
        $clusterdir = $ref_genome;
        $clusterdir =~ s/_//g;
    }

    if($seqfolder && $seqfolder !~ /\/$/){
        $seqfolder .= '/';
    }

    if ($NOSINGLES) {
        $params .= "_nosingles";
    }      	

    if(!$ref_genome) {
        die "# ERROR: please set -r reference_genome\n";
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

		# save list of input file names for future reference
		open(INLOG,">","$outfolder/input.log") ||
			die "# ERROR: cannot create $outfolder/input.log\n";
		foreach $infile (@infiles) {
			print INLOG "$infile\n";
		}
		close(INLOG);

        # create $clusterdir with $params
        $clusterdir .= $params;
        if ( !-e "$outfolder/$clusterdir" ) {
            if ( !mkdir("$outfolder/$clusterdir") ) {
                die "# ERROR: cannot create $outfolder/$clusterdir\n";
            }
        }
    }
    else {
        print
          "# ERROR: need a valid output folder, such as -f Brassicaceae\n\n";
        exit;
    }

    print "# $0 -r $ref_genome -f $outfolder -g $GROWTH -S $NOSINGLES -v $verbose\n";
    print "# ";
    foreach $infile (@infiles) {
        print "-T $infile ";
    } print "\n";

    if (@ignore_species) {
        print "# ";
        foreach $sp (@ignore_species) {
            $ignore{$sp} = 1;
            print "-i $sp ";
        } print "\n";
    }
}

## 1) check (supported) species in TSV file(s)

foreach $infile (@infiles){
    open(TSV,"<",$infile) || die "# ERROR: cannot read $infile\n";
    while(<TSV>){
        my @F = split(/\t/,$_);
        ($sp, $sp2) = ($F[2], $F[7]);
        $species{$sp} = 1;
        $species{$sp2} = 1;
    }
    close(TSV);
}

# add species but skip ignored ones and ref
for $sp (keys(%species)) {

    next if( $ignore{$sp} || $sp eq 'species' || $sp eq 'homology_species');

    $supported{ $sp } = 1;
    if ( !$ref_genome || $sp ne $ref_genome ) {
        push( @supported_species, $sp );
    }
}

# check reference genome is supported
if ( $ref_genome && !$supported{$ref_genome} ) {
    die "# ERROR: cannot find $ref_genome in input TSV file(s)\n";
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

if($show_supported) {
    foreach $sp (@supported_species) {
        print "$sp\n";
    }
	exit(0);
}

$n_of_species = scalar(@supported_species);
print "# total selected species : $n_of_species\n\n";

## 2) infer pairs of collinear gens and make up clusters 

my ( $cluster_id, $chr, $seqtype ) = ( 0, '' );
my ( @cluster_ids, %sorted_cluster_ids );
my ( $ref_geneid, $ref_fasta );
my ( %incluster, %cluster, %sequence );
my ( %totalgenes, %totalclusters, %POCP_matrix );
my ( %sorted_ids, %id2chr );

# columns of TSV file as produced by get_collinear_genes.pl
# Os01g0100200 Os01g0100200 oryza_sativa 1050 ortholog_collinear ONIVA01G00100 ONIVA01G00100 oryza_nivara 1050 NULL NULL NULL 100.00 1 1:11217-12435;1:2929-12267
# NOTE: identity is really alignment length
# NOTE: dn,ds,goc are not computed and have dummy values
my (
    $gene_stable_id,     $prot_stable_id, $species,
    $identity,           $homology_type,  $hom_gene_stable_id,
    $hom_prot_stable_id, $hom_species,    $hom_identity,
    $dn,                 $ds,             $goc_score,
    $wga_coverage,       $high_confidence,
	$coordinates
);

# Iteratively get and parse TSV files that define pairs of collinear 
# genes computed with get_collinear_genes. After parsing all pairwise
# TSV files clusters emerge.
foreach $infile (@infiles) {

    # uncompress on the fly and parse
    open(TSV, "<", $infile)
      || die "# ERROR: cannot open $infile\n";
    while ( my $line = <TSV> ) {

        (
            $gene_stable_id,     $prot_stable_id, $species,
            $identity,           $homology_type,  $hom_gene_stable_id,
            $hom_prot_stable_id, $hom_species,    $hom_identity,
            $dn,                 $ds,             $goc_score,
            $wga_coverage,       $high_confidence,
            $coordinates		
        ) = split( /\t/, $line );

        next if ( !$supported{$species} || !$supported{$hom_species} );

        if ( $homology_type =~ m/ortholog/ ) {

            # add $species protein to cluster only if not clustered yet
            if ( !$incluster{$gene_stable_id} ) {

                if ( $incluster{$hom_gene_stable_id} ) {

                    # use existing cluster_id from other species ortholog
                    $cluster_id = $incluster{$hom_gene_stable_id};
                }
                else {

                    # otherwise create a new one
                    $cluster_id = $gene_stable_id;
                    push( @cluster_ids, $cluster_id );
                }

                # record to which cluster this gene belongs
                $incluster{$gene_stable_id} = $cluster_id;

                push( @{ $cluster{$cluster_id}{$species} }, $gene_stable_id );

            }
            else {
                # set cluster for $hom_species anyway
                $cluster_id = $incluster{$gene_stable_id};
            } 

            # now add $hom_species protein to previously defined cluster
            if ( !$incluster{$hom_gene_stable_id} ) {

                # record to which cluster this protein belongs
                $incluster{$hom_gene_stable_id} = $cluster_id;

                push(
                    @{ $cluster{$cluster_id}{$hom_species} },
                    $hom_gene_stable_id
                );
            }
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

    foreach $seqtype (@SEQTYPE) {

        $filename = "$seqfolder$sp$SEQEXT{$seqtype}";
        if(!-s $filename){
            die "# ERROR: cannot find sequence file $filename, set path with -seq\n"; 
        }

        ( $ref_geneid, $ref_fasta ) = parse_sequence_FASTA_file( $filename );

        $sorted_ids{$sp} = $ref_geneid;

        # count number of genes in this species
        $totalgenes{$sp} = scalar(@$ref_geneid);

        # save these sequences
        foreach $gene_stable_id ( @$ref_geneid ) {
            $sequence{$sp}{$gene_stable_id}{$seqtype} = $ref_fasta->{$gene_stable_id};
        }

        # log
        open(INLOG,">>","$outfolder/input.log") ||
            die "# ERROR: cannot create $outfolder/input.log\n";
        print INLOG "$filename\n";
        close(INLOG);
    }
}

# add unaligned sequences as singletons 
my $total_seqs = 0;
foreach $sp (@supported_species) {

    my $singletons = 0;

    foreach $gene_stable_id ( @{ $sorted_ids{$sp} } ) {

        next if ( $NOSINGLES || $incluster{$gene_stable_id} );

        # create new cluster
        $cluster_id = $gene_stable_id;
        $incluster{$gene_stable_id} = $cluster_id;

        push( @{ $cluster{$cluster_id}{$sp} }, $gene_stable_id );
        push( @cluster_ids,                    $cluster_id );

        # add this singleton to total clusters
        $totalclusters{$sp}++;
        $singletons++;
    }

    $total_seqs += $totalgenes{$sp};

    printf( "# %s : sequences = %d clusters = %d (singletons = %d)\n",
        $sp, $totalgenes{$sp}, $totalclusters{$sp}, $singletons );
}

printf( "\n# total sequences = %d\n\n", $total_seqs );

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

    # create all types of sequence clusters
    foreach $seqtype (@SEQTYPE) {

        my ( %cluster_stats, @cluster_species );
        $n_cluster_sp = $n_cluster_seqs = 0;

        $filename = $cluster_id;

        # write sequences and count sequences
        open( CLUSTER, ">", "$outfolder/$clusterdir/$filename$SEQEXT{$seqtype}" )
          || die "# ERROR: cannot create $outfolder/$clusterdir/$filename$SEQEXT{$seqtype}\n";

	    foreach $species (@supported_species) {
            next if ( !$cluster{$cluster_id}{$species} );

			$n_cluster_sp++;
            foreach $gene_stable_id ( @{ $cluster{$cluster_id}{$species} } ) {
                next if(!$sequence{$species}{$gene_stable_id}{$seqtype}); 

                print CLUSTER $sequence{$species}{$gene_stable_id}{$seqtype};
		
                # use cDNA seq type to compute stats
                if($seqtype eq 'cdna'){
                    $n_cluster_seqs++;
                    $cluster_stats{$species}++;
                }
            }
        }
        close(CLUSTER);

        # cluster summary and PCOP update
		# done with cDNAs, after cds & pep have been processed
        if($seqtype eq 'cdna'){
            @cluster_species = keys(%cluster_stats);

            $cdsfile = $filename.$SEQEXT{'cds'};
            $pepfile = $filename.$SEQEXT{'pep'};

            if ( !-s "$outfolder/$clusterdir/$cdsfile" ) { $cdsfile = 'void' }
            if ( !-s "$outfolder/$clusterdir/$pepfile" ) { $pepfile = 'void' }

            print CLUSTER_LIST
              "cluster $cluster_id size=$n_cluster_seqs taxa=$n_cluster_sp ".
              "cdnafile: $filename$SEQEXT{$seqtype} cdsfile: $cdsfile pepfile: $pepfile\n";

            foreach $species (@cluster_species) {
                foreach $gene_stable_id ( @{ $cluster{$cluster_id}{$species} } ) {
                    print CLUSTER_LIST ": $species\n";
                }
            }

            # add sequences in this cluster from a pair of species/taxa
            foreach $sp ( 0 .. $#cluster_species - 1 ) {
                foreach $sp2 ( $sp + 1 .. $#cluster_species ) {

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

# Note: clusters are not necessarily ordered
push(@{ $sorted_cluster_ids{'unsorted'} }, @cluster_ids );

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
        print PANGEMATRIX "\t$cluster_id"; 
    }
}	
print PANGEMATRIX "\n";

print PANGENEMATRIX "source:$outfolder/$clusterdir";
foreach $chr (keys(%sorted_cluster_ids)) {
    print PANGENEMATRIX "\tchr$chr";
    foreach $cluster_id (@{ $sorted_cluster_ids{$chr} }) {
        print PANGENEMATRIX "\t$cluster_id"; 
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


#######################################################

# Takes the name string of a FASTA file created by cut_sequences.pl
# and parses the sequences in there. Assumes the following header
# format: ">mrnaid geneid coords [production_name]" and thus supports
# the same gene having several associated sequences.
# Returns:
# i) ref to list with coord-sorted gene ids
# ii) ref to hash with FASTA strings with genes as keys
sub parse_sequence_FASTA_file {

    my ( $fname ) = @_;
    my ( $geneid, @geneids, %fasta );

    open(FASTA,"<",$fname) || 
        die "# ERROR(parse_sequence_FASTA_file): cannot read $fname\n";
    while(<FASTA>){
        if(/^>\S+\s+(\S+)\s+\S+\s+\[[^\]]+\]/){
            $geneid = $1;

            # conserved gene id order			
            if(!$fasta{$geneid}){
                push(@geneids,$geneid);
            }

            $fasta{$geneid} .= $_;
        } else {
            $fasta{$geneid} .= $_;
        }
    }
    close(FASTA);

    return ( \@geneids, \%fasta );
}
