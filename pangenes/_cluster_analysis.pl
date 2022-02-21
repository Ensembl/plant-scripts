#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);

# Makes pan-gene analysis based on clusters of collinear genes shared by
# species in a pre-computed Minimap2/Wfmash synteny TSV file, produced by 
# _collinear_genes.pl
# Adapted from https://github.com/eead-csic-compbio/get_homologues

# Copyright [2021-22]
# EMBL-European Bioinformatics Institute & Estacion Experimental Aula Dei-CSIC

# ./_cluster_analysis.pl -T rices4_pangenes/tmp/mergedpairs.tsv -f folder \
#   -s rices4_pangenes/ -r Oryza_sativa.IRGSP-1.0 -g 8 -t 3

my $BEDTOOLSEXE = 'bedtools';
my $TRANSPOSEXE =
'perl -F\'\t\' -ane \'$F[$#F]=~s/\n//g;$r++;for(1 .. @F){$m[$r][$_]=$F[$_-1]};'
  . '$mx=@F;END{for(1 .. $mx){for $t(1 .. $r){print"$m[$t][$_]\t"}print"\n"}}\'';

# these globals match sequence type to extension
my @SEQTYPE = qw( cds pep cdna );
my %SEQEXT = (
    'gdna' => '.gdna.fna',
    'cdna' => '.cdna.fna',
    'cds'  => '.cds.fna',
    'pep'  => '.cds.faa'
);

# genome composition report
my $RNDSEED          = 12345;
my $NOFSAMPLESREPORT = 10;

my ( $ref_genome, $seqfolder, $clusterdir ) = ( '', '', '' );
my ( $outfolder, $params, $bedtools_path) = ('', '', '');
my ( $help, $sp, $sp2, $show_supported, $seed );
my ( $infile, $filename, $cdsfile, $pepfile, $gdnafile );
my ( $n_core_clusters, $n_cluster_sp, $n_cluster_seqs ) = ( 0, 0, 0 );
my ( $NOSINGLES , $dogrowth ) = ( 0, 0 );
my ( $n_of_species, $verbose, $min_taxa ) = ( 0, 0, 0 );
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
    "growth|g=i"    => \$dogrowth,
    "folder|f=s"    => \$outfolder,
    "seq|s=s"       => \$seqfolder,
    "mintaxa|t=i"   => \$min_taxa,
    #"soft|z"        => \$dosoft,
    "seed|R=i"      => \$seed,
    "bedtools|B=s"  => \$bedtools_path
) || help_message();

sub help_message {
    print "\nusage: $0 [options]\n\n"
      . "-T input collinear TSV file(s)             (required, example: -T Minimap2.homologies.rice.overlap0.5.tsv -T ...)\n"
      . "-f output folder                           (required, example: -f myfolder)\n"
      . "-r reference species_name to name clusters (required, example: -r arabidopsis_thaliana)\n"
      . "-l list supported species in -T file       (optional, example: -l)\n"
      . "-i ignore species_name(s)                  (optional, example: -i selaginella_moellendorffii -i ...)\n"
      . "-g do pangene set growth simulations       (optional, example: -g 10; makes [core|pan_gene]*.tab files, ignores -t)\n" 
      . "-S skip singletons                         (optional, by default unclustered sequences are taken)\n"
      . "-s folder with gene seqs of species in TSV (optional, default current folder, files created by _cut_sequences.pl)\n"
      . "-t consider only clusters with -t taxa     (optional, by default all clusters are taken)\n"
      . "-R random seed for genome growth analysis  (optional, requires -g, example -R 1234)\n"
      . "-B path to bedtools binary                 (optional, default: -B bedtools)\n"
      #. "-z add soft-core to genome growth analysis (optional)\n"
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

    if ( $dogrowth && $dogrowth > 0) {
        $NOFSAMPLESREPORT = $dogrowth;

        if($seed) { $RNDSEED = $seed }
    }

    if ($outfolder) {
        if ( -e $outfolder ) {
            print "\n# WARNING : folder '$outfolder' exists, files might be overwritten\n\n";
        } else {
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

    print "# $0 -r $ref_genome -f $outfolder -g $dogrowth -S $NOSINGLES ".
      "-v $verbose -t $min_taxa -R $RNDSEED -B $bedtools_path\n";

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

# check binaries
if(!$bedtools_path) { 
    $bedtools_path = $BEDTOOLSEXE
} 
if(`$bedtools_path` !~ 'sage') {
        print "# ERROR: cannot find binary file $bedtools_path , exit\n";
        exit(-1)
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

## 2) parse pairs of collinear genes and make up clusters 

my ( $cluster_id, $chr, $seqtype, $coords_id ) = ( 0, '' );
my ( $line, $stable_id, $segment_species, $coords, $segment_coords );
my ( @cluster_ids, %sorted_cluster_ids, %segment_species );
my ( $ref_geneid, $ref_fasta, $segment_cluster, $num_segments );
my ( %incluster, %cluster, %coords, %sequence, %segment, %segment_sequence );
my ( %totalgenes, %totalclusters, %POCS_matrix );
my ( %sorted_ids, %id2chr );

# columns of TSV file as produced by get_collinear_genes.pl

#gene:BGIOSGA002571 gene:BGIOSGA002571 Oryza_indica.ASM465v1.chr1 2418 ortholog_collinear
# gene:ONIVA01G00120 gene:ONIVA01G00120 Oryza_nivara_v1.chr1 2418 NULL NULL NULL 100.00 1 1:39544-42130(+);1:116435-120177(+)
#Oryza_indica.ASM465v1.chr1:1:222338-228913 segment Oryza_indica.ASM465v1.chr1 6575 segment_collinear
# gene:ONIVA01G00200 gene:ONIVA01G00200 Oryza_nivara_v1.chr1 6575 NULL NULL NULL 100.00 1 1:222338-228913(+);1:160018-166571(+)

# NOTE: wga_cov,dn,ds,goc are not computed and have dummy values

my (
    $gene_stable_id,     $prot_stable_id, $species,
    $overlap,            $homology_type,  $hom_gene_stable_id,
    $hom_prot_stable_id, $hom_species,    $hom_identity,
    $dn,                 $ds,             $goc_score,
    $wga_coverage,       $high_confidence,
    $coordinates
);

# Iteratively get & parse TSV files that define pairs of collinear genes (made with _collinear_genes.pl) 
# Note: Clusters emerge after parsing all pairwise TSV files, see heuristic below
print "# parsing TSV files\n";
foreach $infile (@infiles) {

    open(TSV, "<", $infile)
      || die "# ERROR: cannot open $infile\n";
    while ( $line = <TSV> ) {

        (
            $gene_stable_id,     $prot_stable_id, $species,
            $overlap,            $homology_type,  $hom_gene_stable_id,
            $hom_prot_stable_id, $hom_species,    $hom_identity,
            $dn,                 $ds,             $goc_score,
            $wga_coverage,       $high_confidence,
            $coordinates		
        ) = split( /\t/, $line );

        next if ( !$supported{$species} || !$supported{$hom_species} );

        if ( $homology_type =~ m/ortholog/ ) {

            # add $species gene to cluster only if not clustered yet
            if ( !$incluster{$gene_stable_id} ) {

                if ( $incluster{$hom_gene_stable_id} ) {

                    # use existing cluster_id from other species ortholog
                    $cluster_id = $incluster{$hom_gene_stable_id};
                } else {

                    # otherwise create a new one
                    $cluster_id = $gene_stable_id;
                    push( @cluster_ids, $cluster_id );
                }

                # record to which cluster this gene belongs
                $incluster{$gene_stable_id} = $cluster_id;

                push( @{ $cluster{$cluster_id}{$species} }, $gene_stable_id );

            } else {
                # set cluster for $hom_species anyway
                $cluster_id = $incluster{$gene_stable_id};
            } 

            # now add $hom_species gene 
            if ( !$incluster{$hom_gene_stable_id} ) {

                # record to which cluster this gene belongs, 
                # currently 1st pair where $hom_gene_stable_id appears (heuristic)
                $incluster{$hom_gene_stable_id} = $cluster_id;

                push( @{ $cluster{$cluster_id}{$hom_species} },
                    $hom_gene_stable_id
                );
            } else {
                # $hom_species gene already clustered (same or different cluster)
                if($cluster_id ne $incluster{$hom_gene_stable_id} && $verbose) {
                    print "# WARN: possibly conflicting clusters for $cluster_id & $incluster{$hom_gene_stable_id}\n";
                    # TODO: merge clusters? Would require a hash of merged cluster_ids
                    # downside: will inflate core subset
                }
            }
                
        } elsif($homology_type =~ m/segment_collinear/) {

            chomp($coordinates);

            # record coords of this gene-segment pair
            if($prot_stable_id eq 'segment') {
                $stable_id = $hom_gene_stable_id;
                $segment_species = $species;
                $species = $hom_species;
                ($segment_coords, $coords) = split(/;/,$coordinates);
            } else {
                $stable_id = $gene_stable_id;
                $segment_species = $hom_species;
                ($coords, $segment_coords) = split(/;/,$coordinates); 
            }

            # assign coords to stable_id
            $segment{$stable_id}{$species} = $coords;
            $segment{$stable_id}{$segment_species} = $segment_coords;
        }
    }
    close(TSV);
} print "\n";  

# count how many clusters include each species
foreach $cluster_id (@cluster_ids) {
    foreach $species (@supported_species) {
        if ( $cluster{$cluster_id}{$species} ) {
            $totalclusters{$species}++;
        }
    }
}  

# 2.1) Write BED & FASTA files with genomic segments (gdna), one per species,
# and store them in a hash with (species,coords) keys
foreach $species (@supported_species) {

    $num_segments = 0;

    # write BED
    $filename = "$seqfolder$species.gdna.bed";
    open(GDNA,">",$filename) ||
        die "# ERROR: cannot create $filename\n";

    foreach $stable_id (keys(%segment)) {

        next if(!defined($segment{$stable_id}{$species}));

        #1:125929-131075(+)
        #UN-Ctg123:12-1234(-)
        $coords_id = $segment{$stable_id}{$species};
        if($coords_id =~ m/^(\S+)?:(\d+)-(\d+)\(([+-])\)/) {
            print GDNA "$1\t$2\t$3\tNA\t0\t$4\n";
            $num_segments++;
        }
    }
    close(GDNA); 

    # extract segments with help from bedtools
    my $FASTAgenome_file = "$seqfolder\_$species.fna";
    my $outFASTAfile = "$seqfolder$species.gdna.fna";

    if($num_segments) {
        my $cmd = "$bedtools_path getfasta -fi $FASTAgenome_file -bed $filename -s -fo $outFASTAfile"; 
        system("$cmd");
        if ( $? != 0 ) {
            die "# ERROR: failed running bedtools ($cmd)\n";
        }
        elsif ( !-s $outFASTAfile ) {
            die "# ERROR: failed generating $outFASTAfile file ($cmd)\n";
        }
    }

    print "# $outFASTAfile : $num_segments genomic segments\n";

    $filename = "$seqfolder$species$SEQEXT{'gdna'}";
    if(!-s $filename){
        print "# WARN: cannot find sequence file $filename, skip it\n" if($verbose);
        next;
    }

    my $ref_fasta = parse_GETFASTA_file( $filename, $species );

    foreach $coords_id ( keys(%$ref_fasta) ) {
        $segment_sequence{$species}{$coords_id} = $ref_fasta->{$coords_id};
    }

    # log
    open(INLOG,">>","$outfolder/input.log") ||
        die "# ERROR: cannot create $outfolder/input.log\n";
    print INLOG "$filename\n";
    close(INLOG);
} print "\n";


# 2.2) Parse FASTA files to get main cluster sequences (@SEQTYPE)
# Note: there might be 1+ sequences for the same gene, same species in GFF
foreach $species (@supported_species) {

    foreach $seqtype (@SEQTYPE) {

        $filename = "$seqfolder$species$SEQEXT{$seqtype}";
        if(!-s $filename){
            die "# ERROR: cannot find sequence file $filename, set path with -seq\n"; 
        }

        ( $ref_geneid, $ref_fasta ) = parse_sequence_FASTA_file( $filename );

        $sorted_ids{$species} = $ref_geneid;

        # count number of genes in this species
        $totalgenes{$species} = scalar(@$ref_geneid);

        # save these sequences
        foreach $gene_stable_id ( @$ref_geneid ) {
            $sequence{$species}{$gene_stable_id}{$seqtype} = $ref_fasta->{$gene_stable_id};
        }

        # log
        open(INLOG,">>","$outfolder/input.log") ||
            die "# ERROR: cannot create $outfolder/input.log\n";
        print INLOG "$filename\n";
        close(INLOG);
    }
}

# 2.3) add unaligned sequences as singletons 
my $total_seqs = 0;
foreach $species (@supported_species) {

    my $singletons = 0;

    foreach $gene_stable_id ( @{ $sorted_ids{$species} } ) {
        
        next if ( $NOSINGLES || $incluster{$gene_stable_id} );

        # create new cluster
        $cluster_id = $gene_stable_id;
        $incluster{$gene_stable_id} = $cluster_id;

        push( @{ $cluster{$cluster_id}{$species} }, $gene_stable_id );
        push( @cluster_ids, $cluster_id );

        # add this singleton to total clusters
        $totalclusters{$species}++;
        $singletons++;
    }

    $total_seqs += $totalgenes{$species};

    printf( "# %s : sequences = %d clusters = %d (singletons = %d)\n",
        $species, $totalgenes{$species}, $totalclusters{$species}, $singletons );
}

printf( "\n# total sequences = %d\n\n", $total_seqs );


# 2.4) create and print shadow clusters of genomic sequences
foreach $cluster_id (@cluster_ids) {

    next if( scalar( keys( %{ $cluster{$cluster_id} } ) ) < $min_taxa);

    if($cluster{$cluster_id}{$ref_genome}) {
        $filename = $cluster{$cluster_id}{$ref_genome}[0]
    } else { $filename = $cluster_id }

    my ($new,%seen) = (0);
    $segment_cluster = '';

    foreach $species (@supported_species) {
        next if ( !$cluster{$cluster_id}{$species} );

        my $segment_seqs_OK = 0;
        foreach $stable_id (@{ $cluster{$cluster_id}{$species} }) {
            next if(!$segment{$stable_id}{$species});

            # sp1, found both in cluster and in segment pair
            if(!$seen{$species}) {
                $coords_id = $segment{$stable_id}{$species};
                $segment_cluster .= $segment_sequence{$species}{$coords_id};
                $segment_seqs_OK++;
                $seen{$species}=1;
            }

            # sp2, which might be missing from cluster
            foreach $sp (@supported_species) { 
                next if(!$segment{$stable_id}{$sp} || $seen{$sp});
                $coords_id = $segment{$stable_id}{$sp};
                $segment_cluster .= $segment_sequence{$sp}{$coords_id};
                $segment_seqs_OK++;
                $seen{$sp}=1;
                if(!$cluster{$cluster_id}{$sp}) {
                    $new++  
                }
                last;
            }
        }    

        if($segment_seqs_OK) {
            $segment_species{$cluster_id} += $segment_seqs_OK;
        }
    } 

    # these are regular clusters; genuine segments have 2+ species 
    if(!$segment_species{$cluster_id} || $segment_species{$cluster_id} < 2 || $new < 1) {
        $segment_species{$cluster_id} = 0;
        next;
    }  

    open( CLUSTER, ">", "$outfolder/$clusterdir/$filename$SEQEXT{'gdna'}" )
        || die "# ERROR: cannot create $outfolder/$clusterdir/$filename$SEQEXT{'gdna'}\n";
    print CLUSTER $segment_cluster;
    close(CLUSTER);
} 

%segment_sequence = (); # not needed anymore


# 3) write main sequence clusters, summary text file and POCS matrix

# POCS=Percent Conserved Sequences (POCS) matrix
my $POCS_matrix_file = "$outfolder/POCS.matrix$params\.tab";

my $cluster_summary_file = "$outfolder/$clusterdir.cluster_list";

open( CLUSTER_LIST, ">", $cluster_summary_file )
  || die "# ERROR: cannot create $cluster_summary_file\n";

$n_core_clusters = 0;

foreach $cluster_id (@cluster_ids) {

    next if(scalar( keys( %{ $cluster{$cluster_id} } ) ) < $min_taxa);

    if ( scalar( keys( %{ $cluster{$cluster_id} } ) ) == $n_of_species ) {
        $n_core_clusters++;
    }

    # create all types of sequence clusters
    foreach $seqtype (@SEQTYPE) {

        my ( %cluster_stats, @cluster_species );
        ($n_cluster_sp, $n_cluster_seqs) = (0, 0);

        if($cluster{$cluster_id}{$ref_genome}) {
            $filename = $cluster{$cluster_id}{$ref_genome}[0]
        } else { $filename = $cluster_id }

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
	# done with cDNAs, after gdna, cds & pep files have been created
        if($seqtype eq 'cdna'){
            @cluster_species = keys(%cluster_stats);

            $cdsfile = $filename.$SEQEXT{'cds'};
            $pepfile = $filename.$SEQEXT{'pep'};
            $gdnafile = $filename.$SEQEXT{'gdna'};

            if ( !-s "$outfolder/$clusterdir/$cdsfile" ) { $cdsfile = 'void' }
            if ( !-s "$outfolder/$clusterdir/$pepfile" ) { $pepfile = 'void' }
            if ( !-s "$outfolder/$clusterdir/$gdnafile" ) { $gdnafile = 'void' }

            $num_segments = $segment_species{$cluster_id} || 'NA'; 

            print CLUSTER_LIST
              "cluster $filename size=$n_cluster_seqs taxa=$n_cluster_sp taxa(gdna)=$num_segments ".
              "cdnafile: $filename$SEQEXT{$seqtype} cdsfile: $cdsfile pepfile: $pepfile gdnafile: $gdnafile\n"; 

            foreach $species (@cluster_species) {
                foreach $gene_stable_id ( @{ $cluster{$cluster_id}{$species} } ) {
                    print CLUSTER_LIST ": $species\n";
                }
            }

            # add sequences in this cluster from a pair of species/taxa
            foreach $sp ( 0 .. $#cluster_species - 1 ) {
                foreach $sp2 ( $sp + 1 .. $#cluster_species ) {

                    $POCS_matrix{ $cluster_species[$sp] }{ $cluster_species[$sp2] } +=
                      $cluster_stats{ $cluster_species[$sp] };
                    $POCS_matrix{ $cluster_species[$sp] }{ $cluster_species[$sp2] } +=
                      $cluster_stats{ $cluster_species[$sp2] };

                    # now in reverse order to make sure it all adds up
                    $POCS_matrix{ $cluster_species[$sp2] }{ $cluster_species[$sp] } +=
                      $cluster_stats{ $cluster_species[$sp] };
                    $POCS_matrix{ $cluster_species[$sp2] }{ $cluster_species[$sp] } +=
                      $cluster_stats{ $cluster_species[$sp2] };
                }				  
            }
        }
    }
}

close(CLUSTER_LIST);

printf( "\n# number of clusters = %d (core = %d)\n\n",
    scalar(@cluster_ids), $n_core_clusters );
print "# cluster_list = $outfolder/$clusterdir.cluster_list\n";
print "# cluster_directory = $outfolder/$clusterdir\n";


# 3.1) print POCS matrix #########################################

open( POCSMATRIX, ">$POCS_matrix_file" )
  || die "# EXIT: cannot create $POCS_matrix_file\n";

print POCSMATRIX "genomes";
foreach $sp ( 0 .. $#supported_species ) {
    print POCSMATRIX "\t$supported_species[$sp]";
}
print POCSMATRIX "\n";

my (%POCS2ref,$perc);

foreach $sp ( 0 .. $#supported_species ) {
    print POCSMATRIX "$supported_species[$sp]";
    foreach $sp2 ( 0 .. $#supported_species ) {

        if ( $sp == $sp2 ) { 
            print POCSMATRIX "\t100.00"
        } else {
            if ( $POCS_matrix{ $supported_species[$sp] }
                { $supported_species[$sp2] } )
            {
                $perc = sprintf("\t%1.2f",
                    (
                        100 * $POCS_matrix{ $supported_species[$sp] }
                          { $supported_species[$sp2] }
                    ) / (
                        $totalgenes{ $supported_species[$sp] } +
                          $totalgenes{ $supported_species[$sp2] }
                    )
                );
                print POCSMATRIX "$perc";

                # save %POCS for all species vs reference
                if($sp == 0){ $POCS2ref{$supported_species[$sp2]} = $perc }
            }
            else {
                print POCSMATRIX "\tNA";
            }
        }
    }
    print POCSMATRIX "\n";
}
close(POCSMATRIX);

print "\n# percent_conserved_proteins_file = $POCS_matrix_file\n\n";

# sort species from ref down by decreasing POCS
my @supported_species_POCS;

# reference goes in 1st place
push(@supported_species_POCS, $ref_genome);

foreach $sp2 (sort {$POCS2ref{$b}<=>$POCS2ref{$a}} keys(%POCS2ref)) {
    push(@supported_species_POCS, $sp2);
}


## 4)  write pangenome matrices in output folder

# Note: clusters are not necessarily ordered
push(@{ $sorted_cluster_ids{'unsorted'} }, @cluster_ids );

# set matrix filenames and write headers
my $pangene_matrix_file = "$outfolder/pangene_matrix$params\.tab";
my $pangene_gene_file   = "$outfolder/pangene_matrix_genes$params\.tab";
my $pangene_matrix_tr   = "$outfolder/pangene_matrix$params\.tr.tab";
my $pangene_gene_tr     = "$outfolder/pangene_matrix_genes$params\.tr.tab";
my $pangene_fasta_file  = "$outfolder/pangene_matrix$params\.fasta";

open( PANGEMATRIX, ">$pangene_matrix_file" )
  || die "# EXIT: cannot create $pangene_matrix_file\n";

open( PANGENEMATRIX, ">$pangene_gene_file" )
  || die "# EXIT: cannot create $pangene_gene_file\n";

print PANGEMATRIX "source:$outfolder/$clusterdir";
foreach $chr (keys(%sorted_cluster_ids)) {
    print PANGEMATRIX "\tchr$chr";
    foreach $cluster_id (@{ $sorted_cluster_ids{$chr} }) {
        $filename = $cluster{$cluster_id}{$ref_genome}[0] || $cluster_id;
        print PANGEMATRIX "\t$filename"; 
    }
}	
print PANGEMATRIX "\n";

print PANGENEMATRIX "source:$outfolder/$clusterdir";
foreach $chr (keys(%sorted_cluster_ids)) {
    print PANGENEMATRIX "\tchr$chr";
    foreach $cluster_id (@{ $sorted_cluster_ids{$chr} }) {
        $filename = $cluster{$cluster_id}{$ref_genome}[0] || $cluster_id;
        print PANGENEMATRIX "\t$filename"; 
    }
}	
print PANGENEMATRIX "\n";

open( PANGEMATRIF, ">$pangene_fasta_file" )
  || die "# EXIT: cannot create $pangene_fasta_file\n";

foreach $species (@supported_species_POCS) {

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

system("$TRANSPOSEXE $pangene_matrix_file > $pangene_matrix_tr");
system("$TRANSPOSEXE $pangene_gene_file > $pangene_gene_tr");

print "# pangene_file (occup) = $pangene_matrix_file tranposed = $pangene_matrix_tr\n";
print "# pangene_file (names) = $pangene_gene_file transposed = $pangene_gene_tr\n";
#print "# pangene_FASTA_file = $pangene_fasta_file\n";

exit if(!$dogrowth);


## 5) optionally make genome composition analysis to simulate pangene set growth
##    as new annotated genomes are added to the pool
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

    print "# adding $tmptaxa[0]: core=$coregenome[$s][0] pan=$pangenome[$s][0]\n"
        if ($verbose);

    for ( $sp = 1 ; $sp < $n_of_species ; $sp++ ) {

        $coregenome[$s][$sp] = 0;
        $pangenome[$s][$sp]  = $pangenome[$s][ $sp - 1 ];
        $core_occup          = $sp + 1;

        foreach $chr (keys(%sorted_cluster_ids)) {
            foreach $cluster_id (@{ $sorted_cluster_ids{$chr} }) {

                # check reference species is in this cluster (1st iteration only)
                if ( $sp == 1 && $cluster{$cluster_id}{ $tmptaxa[0] } && 
                    scalar(@{ $cluster{$cluster_id}{ $tmptaxa[0] } }) > 0 ) {
                    $n_of_taxa_in_cluster{$cluster_id}++
                }

                # check $sp is represented in this cluster
                if ( $cluster{$cluster_id}{ $tmptaxa[$sp] } &&
                    scalar(@{ $cluster{$cluster_id}{ $tmptaxa[$sp] } }) > 0 ) {
                    $n_of_taxa_in_cluster{$cluster_id}++
                }

                # work out cluster occupancy
                if ( $n_of_taxa_in_cluster{$cluster_id} &&
                    $cluster{$cluster_id}{ $tmptaxa[$sp] } &&
                    scalar(@{ $cluster{$cluster_id}{ $tmptaxa[$sp] } }) > 0 ) {

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

# Takes the name string of a FASTA file created by _cut_sequences.pl
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

# Takes the name string of a FASTA file created by bedtools getfasta
# and parses the sequences in there. Assumes the following header
# format: ">chr:start-end(strand)"
# Returns:
# i) ref to hash with FASTA strings with coord as keys
sub parse_GETFASTA_file {

    my ( $fname, $species_name ) = @_;
    my ( $geneid, %fasta );

    open(FASTA,"<",$fname) ||
        die "# ERROR(parse_GETFASTA_file): cannot read $fname\n";
    while(<FASTA>){
        if(/^>(\S+)/){
            chomp;
            $geneid = $1;
            $fasta{$geneid} = "$_ [$species_name]\n";
        } else {
            $fasta{$geneid} .= $_;
        }
    }
    close(FASTA);

    return \%fasta;
}



sub factorial {
  my $max = int($_[0]);
  my $f = 1;
  for (2..$max) { $f *= $_ }
  return $f;
}

# based on http://www.unix.org.ua/orelly/perl/cookbook/ch04_18.htm
# generates a random permutation of @array in place and returns 
# string with concatenated elements
sub fisher_yates_shuffle {
  my ($array) = (@_);
  my ($i,$j);

  for ($i = @$array; --$i; )
  {
    $j = int(rand($i+1));
    next if $i == $j;
    @$array[$i,$j] = @$array[$j,$i];
  }

  return join('',@$array);
}

# write TSV file appropriate for R boxplot function
sub write_boxplot_file {

    my ( $outfile, $n_genomes, $n_samples, $ref_data ) = @_;

    my ( $s, $sp );

    open( BOXDATA, ">", $outfile )
      || die "# ERROR(write_boxplot_file): cannot create $outfile\n";
    for ( $sp = 0 ; $sp < $n_genomes ; $sp++ ) {
        printf( BOXDATA "g%d\t", $sp + 1 );    #g = genome
    }
    print BOXDATA "\n";

    for ( $s = 0 ; $s < $n_samples ; $s++ ) {
        for ( $sp = 0 ; $sp < $n_genomes ; $sp++ ) {
            print BOXDATA "$ref_data->[$s][$sp]\t";
        }
        print BOXDATA "\n";
    }
    close(BOXDATA);

    return $outfile;
}

