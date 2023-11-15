#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);
use FindBin '$Bin';
use lib "$Bin/lib";
use pangeneTools qw(parse_sequence_FASTA_file);

# Makes pan-gene analysis based on clusters of collinear genes shared by
# species in a pre-computed Minimap2/GSAlign synteny TSV file, produced by 
# _collinear_genes.pl
# Adapted from https://github.com/eead-csic-compbio/get_homologues

# Copyright [2021-23]
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

# cluster quality control
my $MAXDISTNEIGHBORS = 5;   # neighbor genes in a cluster cannot be more than N genes away
my $MINEDGESTOMERGE  = 0.75;# ratio of edges connecting two clusters so they can be merged

# genome composition report
my $RNDSEED          = 12345;
my $NOFSAMPLESREPORT = 5;

my ( $ref_genome, $seqfolder, $clusterdir ) = ( '', '', '' );
my ( $outfolder, $params, $bedtools_path) = ('', '', '');
my ( $help, $sp, $sp2, $show_supported, $seed );
my ( $infile, $filename, $cdsfile, $pepfile, $gdnafile );
my ( $n_core_clusters, $n_cluster_sp, $n_cluster_seqs ) = ( 0, 0, 0 );
my ( $maxdistneigh ) = ( $MAXDISTNEIGHBORS );
my ( $NOSINGLES, $dogrowth, $patch, $chregex ) = ( 0, 0, 0, '' );
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
    "maxdist|m=i"   => \$maxdistneigh,
    "regex|x=s"     => \$chregex,
    "patch|p"       => \$patch,
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
      . "-m max distance among neighbor genes       (optional, example: -m 10, default: $maxdistneigh)\n"
      . "-x regex to match chromosomes in genome    (optional, ie: -x '^\\d+\$')\n"
      . "-p use patched sequences                   (optional, expects patch.SEQEXT extension)\n"
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

    if($maxdistneigh < 0) {
        print "# ERROR: -m must be a natural number\n\n";
        exit;
    } else {
        $maxdistneigh = int($maxdistneigh);
    } 

    if ($NOSINGLES) {
        $params .= "_nosingles";
    }      	

    #if ($patch) { $params .= "_patch" }

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

    print "# $0 -r $ref_genome -f $outfolder -g $dogrowth -S $NOSINGLES -m $maxdistneigh ".
      "-v $verbose -t $min_taxa -x $chregex -p $patch -R $RNDSEED -B $bedtools_path\n";

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

    print "\n";
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

    $n_of_species++;
    $supported{ $sp } = $n_of_species;
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
            print "# $sp $supported{ $sp }\n";
        }
    }
}

if($show_supported) {
    foreach $sp (@supported_species) {
        print "$sp\n";
    }
    exit(0);
}

print "\n# clustering parameters:\n";
print "# \$MAXDISTNEIGHBORS: $MAXDISTNEIGHBORS\n";
print "# \$MINEDGESTOMERGE: $MINEDGESTOMERGE\n\n";

print "# total selected species : $n_of_species\n\n";

## 2) parse pairs of collinear genes and make up clusters 

my ( $chr, $cluster_id, $cluster_id2, $seqtype, $coords_id ) = ( '', 0 );
my ( $line, $stable_id, $segment_species, $coords, $segment_coords );
my ( $ref_geneid, $ref_fasta, $ref_coords, $segment_cluster, $num_segments );
my ( @sorted_chrs, %sorted_cluster_ids, %segment_species );
my ( %incluster, %cluster, %sequence, %segment, %segment_sequence );
my ( %totalgenes, %totalclusters, %totaloverlap, %POCS_matrix );
my ( %unclustered, %sorted_ids, %id2coords );
my ( %cluster_links, %toremove, @tmp_cluster_ids, @cluster_ids );
my (
    $gene_stable_id,     $prot_stable_id, $species,
    $overlap,            $homology_type,  $hom_gene_stable_id,
    $hom_prot_stable_id, $hom_species,    $hom_identity,
    $dn,                 $ds,             $goc_score,
    $wga_coverage,       $high_confidence,
    $coordinates
);

# 2.1) Parse FASTA files to get main cluster sequences (@SEQTYPE),
#      cDNA file is is where genes get their genomic coordinates
# Note1: if a gene lacks a cDNA will be skipped
# Note2: there might be 1+ sequences for the same gene, same species in GFF
foreach $species (@supported_species) {

    # init, see below
    $unclustered{$species} = 0;

    foreach $seqtype (@SEQTYPE) {

        $filename = "$seqfolder$species$SEQEXT{$seqtype}";
        if($patch) {
            $filename = "$seqfolder$species.patch$SEQEXT{$seqtype}";
        }

        if(!-s $filename){
            die "# ERROR: cannot find sequence file $filename, set path with -seq\n";
        }

        ( $ref_geneid, $ref_fasta, $ref_coords ) = parse_sequence_FASTA_file( $filename );

        # save these sequences
        foreach $gene_stable_id ( @$ref_geneid ) {
            $sequence{$species}{$gene_stable_id}{$seqtype} = $ref_fasta->{$gene_stable_id};
        }

        # log
        open(INLOG,">>","$outfolder/input.log") ||
            die "# ERROR: cannot create $outfolder/input.log\n";
        print INLOG "$filename\n";
        close(INLOG);

        # take gene coordinates only once,
        # only genes with cDNA sequences are considered.
        # WARNING: rRNA, tRNA genes don't have cDNA features
        if($seqtype eq 'cdna') {
            $sorted_ids{$species} = $ref_geneid;
            $id2coords{$species} = $ref_coords;

            # count number of genes in this species
            $totalgenes{$species} = scalar(@$ref_geneid);
        }
    }
}

# 2.2) actually parse pairs of collinear genes to grow clusters

# columns of TSV file as produced by get_collinear_genes.pl:
#gene:BGIOSGA002571 gene:BGIOSGA002571 Oryza_indica.ASM465v1.chr1 2418 ortholog_collinear
# gene:ONIVA01G00120 gene:ONIVA01G00120 Oryza_nivara_v1.chr1 2418 NULL NULL NULL 100.00 1 1:39544-42130(+);1:116435-120177(+)
#Oryza_indica.ASM465v1.chr1:1:222338-228913 segment Oryza_indica.ASM465v1.chr1 6575 segment_collinear
# gene:ONIVA01G00200 gene:ONIVA01G00200 Oryza_nivara_v1.chr1 6575 NULL NULL NULL 100.00 1 1:222338-228913(+);1:160018-166571(+)

# NOTE: wga_cov,dn,ds,goc are not computed and have dummy values

# Iteratively get & parse TSV files that define pairs of collinear genes (made with _collinear_genes.pl) 
# Clusters are first created as pairs and grow as new collinear genes from other species are added.
# In some cases selected clusters will be merged on a second step.
# Note: Assumes TSV files have been sorted by $gene_stable_id and $overlap (see get_pangenes.pl) and
# warns (-v) about conflicting (with sequences from same species) or 
# partially overlapping clusters (disjoint species, not enough links/edges to merge) 

print "# parsing TSV files\n";
foreach $infile (@infiles) {

    open(TSV, "<", $infile)
      || die "# ERROR: cannot open $infile\n";
    while ( $line = <TSV> ) {

        (   $gene_stable_id,     $prot_stable_id, $species,
            $overlap,            $homology_type,  $hom_gene_stable_id,
            $hom_prot_stable_id, $hom_species,    $hom_identity,
            $dn,                 $ds,             $goc_score,
            $wga_coverage,       $high_confidence,
            $coordinates		
        ) = split( /\t/, $line );

        next if ( !$supported{$species} || !$supported{$hom_species} );

        if ( $homology_type =~ m/ortholog/ ) {

            if($verbose) {
                if(!defined($id2coords{$species}{$gene_stable_id}[1])) {
                    print "# WARN: skip gene model $species $gene_stable_id as it lacks cDNA (no coordinates)\n";
                    next;
                } elsif(!defined($id2coords{$hom_species}{$hom_gene_stable_id}[1])) {
                    print "# WARN: skip gene model $hom_species $hom_gene_stable_id as it lacks cDNA (no coordinates)\n";
                    next;
                }
            }

            # add $species gene to cluster only if not clustered yet
            if ( !$incluster{ $supported{$species}.':::'.$gene_stable_id } ) {

                if ( $incluster{ $supported{$hom_species}.':::'.$hom_gene_stable_id } ) {

                    # use existing cluster_id from other species ortholog
                    $cluster_id = $incluster{ $supported{$hom_species}.':::'.$hom_gene_stable_id };
                } else {

                    # otherwise create a new one
                    $cluster_id = $supported{$species}.':::'.$gene_stable_id;
                    push( @tmp_cluster_ids, $cluster_id );
                }

                # record to which cluster this gene belongs
                $incluster{ $supported{$species}.':::'.$gene_stable_id } = $cluster_id;

                push( @{ $cluster{$cluster_id}{$species} }, $gene_stable_id );

                # record overlap
                $totaloverlap{$gene_stable_id} += $overlap;

            } else {
                # set cluster for $hom_species anyway
                $cluster_id = $incluster{ $supported{$species}.':::'.$gene_stable_id };
            } 

            # now add $hom_species gene 
            if ( !$incluster{ $supported{$hom_species}.':::'.$hom_gene_stable_id } ) {

                # record to which cluster this gene belongs, 
                # currently 1st pair where $hom_gene_stable_id appears (heuristic)
                $incluster{ $supported{$hom_species}.':::'.$hom_gene_stable_id } = $cluster_id;

                push( @{ $cluster{$cluster_id}{$hom_species} }, $hom_gene_stable_id);

                # record overlap
                $totaloverlap{$hom_gene_stable_id} += $overlap;

            } else { # $hom_species gene already clustered 

                # typically when a long gene model overlaps 2+ shorther ones
                if($cluster_id ne $incluster{ $supported{$hom_species}.':::'.$hom_gene_stable_id }) {

                    $cluster_links{ $cluster_id }{ $incluster{ $supported{$hom_species}.':::'.$hom_gene_stable_id } }++;
                    $cluster_links{ $incluster{ $supported{$hom_species}.':::'.$hom_gene_stable_id } }{ $cluster_id }++;
                } 
                #else { } # do nothing, previously added to same cluster
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


# 2.3) see if disjoint clusters can be merged, can happen if long genes overlap shorther ones
my ($c1, $c2, $size1, $size2, $speciesOK);

foreach $c1 (0 .. $#tmp_cluster_ids-1) {

    $cluster_id = $tmp_cluster_ids[$c1];

    next if(defined($toremove{ $cluster_id }));

    foreach $c2 ($c1+1 .. $#tmp_cluster_ids) {	
        $cluster_id2 = $tmp_cluster_ids[$c2];
         
        if(!defined($cluster_links{ $cluster_id }{ $cluster_id2 }) ||
            defined($toremove{ $cluster_id2 })) {
            next;
        }

        # check species in clusters are different (disjoint clusters)
        ($speciesOK, $size1) = (1, 0);
        foreach $species (keys(%{ $cluster{$cluster_id} })) {
            if(defined($cluster{$cluster_id2}{$species})){
                $speciesOK = 0;					
                last;
            }
            # get size of cluster1 as you go  
            $size1 += scalar(@{ $cluster{$cluster_id}{$species} });
        }

        if($speciesOK == 0){
            print "# WARN: conflicting clusters $cluster_id & $cluster_id2 ($species)\n" if($verbose);
            next;
        }

        # get size of cluster2
        $size2 = 0;
        foreach $species (keys(%{ $cluster{$cluster_id2} })) {
            $size2 += scalar(@{ $cluster{$cluster_id2}{$species} });
        }			
  
        # is there evidence to merge these clusters?
        if($cluster_links{ $cluster_id }{ $cluster_id2 } >= $MINEDGESTOMERGE * $size1 * $size2) {

            # actually merge cluster2 to cluster1, requires updating data structures 
            foreach $species (keys(%{ $cluster{$cluster_id2} })) {
                foreach $hom_gene_stable_id (@{ $cluster{$cluster_id2}{$species} }) {

                    # copy gene ids from cluster2 to cluster, for all species
                    push(@{ $cluster{$cluster_id}{$species} }, $hom_gene_stable_id);

                    # change to which cluster_id those genes now belong
                    $incluster{ $supported{$species}.':::'.$hom_gene_stable_id } = $cluster_id;
                }
            }

            delete($cluster{$cluster_id2});
            $toremove{ $cluster_id2 } = 1;

            print "# WARN: merged clusters $cluster_id & $cluster_id2 " .
                "($cluster_links{ $cluster_id }{ $cluster_id2 },$size1,$size2)\n" if($verbose);
        } else {
            print "# WARN: partially overlapping clusters $cluster_id & $cluster_id2 " .
                "($cluster_links{ $cluster_id }{ $cluster_id2 },$size1,$size2)\n" if($verbose);
        }
    }		
} print "\n";

# remove individual clusters (merged) from cluster list
foreach $cluster_id (@tmp_cluster_ids) {
    next if($toremove{ $cluster_id });
    push(@cluster_ids,$cluster_id);
}

@tmp_cluster_ids=();
undef(%toremove);

#foreach $cluster_id (@cluster_ids) {
#	next if($cluster_id !~ '1:::BaRT2v18chr3HG163450'); && 
#	foreach $species (@supported_species) {
#		foreach $gene_stable_id ( @{ $cluster{$cluster_id}{$species} } ) {
#	print "$cluster_id $gene_stable_id\n"; }}}



# 2.4) Write BED & FASTA files with genomic segments (gdna), one per species,
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
        sleep(2);
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


# 2.5) quality control of clusters, criteria: 
# i)  genes from same species should be neighbors, else should be removed
# ii) genes from same species should be ordered along chr

my ($best_index, $best_gene_stable_id, $index, $index_dist);
my ($main_cluster_id, $new_cluster_id, $gene_id);

foreach $cluster_id (@cluster_ids) {

    #next if($cluster_id !~ 'Os01g0116250'); #1:::HORVU.MOREX.r3.3HG0311160'); #gene:BGIOSGA000009'); #debug

    # set main cluster id, note it can be updated within the loop
    $main_cluster_id = $cluster_id;

    foreach $species (@supported_species) {

        if( defined($cluster{$main_cluster_id}{$species}) &&
            scalar(@{ $cluster{$main_cluster_id}{$species} }) > 1 ) {

            my (@checked_ids);

            # rank genes in terms of cumulative overlap
            my @ranked_ids = sort {$totaloverlap{$b}<=>$totaloverlap{$a}}
                @{ $cluster{$main_cluster_id}{$species} };
 
            # get index of gene with highest overlap
            $best_gene_stable_id = shift(@ranked_ids); #$best_gene_stable_id = pop(@ranked_ids); #debug
            $best_index = _get_element_index($sorted_ids{$species},$best_gene_stable_id);
            push(@checked_ids, $best_gene_stable_id);

            # iteratively compute index distance of remaining genes 
            foreach $gene_stable_id (@ranked_ids) {

                    $index = _get_element_index($sorted_ids{$species},$gene_stable_id);
                    $index_dist = abs($index - $best_index); 

                    # remove/uncluster genes that are no neighbors
                    if($index_dist > $maxdistneigh) { 
                    #if($index_dist > $maxdistneigh || $gene_stable_id eq 'gene:BGIOSGA000009') { # debug

                        # create separate singleton cluster for this gene

                        # original cluster was named after this gene,
                        # we'll need to create a new one and set it as main cluster
                        # NOTE: this can only happen once per cluster
                        if(defined($cluster{ $supported{$species}.':::'.$gene_stable_id })) {

                            # copy original cluster to new main cluster id
                            $new_cluster_id = $supported{$species}.':::'.$best_gene_stable_id;                            

                            foreach $sp2 (@supported_species) { 
                                next if($sp2 eq $species);
                                foreach $gene_id (@{ $cluster{$cluster_id}{$sp2} }) {
                                    push(@{ $cluster{$new_cluster_id}{$sp2}}, $gene_id);
                                    $incluster{ $supported{$sp2}.':::'.$gene_id } = $new_cluster_id;
                                } 
                            } 

                            # update main cluster id and add it to list
                            $main_cluster_id = $new_cluster_id;
                            push( @cluster_ids, $main_cluster_id );

                            # replace original cluster with singleton 
                            $cluster{$cluster_id} = ();
                            push(@{ $cluster{$cluster_id}{$species}}, $gene_stable_id); 

                            if($verbose) {
                                print "# WARN: remove $gene_stable_id from cluster $cluster_id => ".
                                    "$main_cluster_id ($index_dist)\n";
                            }
 
                        } else { # new cluster with new cluster_id

                            $new_cluster_id = $supported{$species}.':::'.$gene_stable_id;
                            $incluster{ $supported{$species}.':::'.$gene_stable_id } = $new_cluster_id;
                            push( @{ $cluster{$new_cluster_id}{$species} }, $gene_stable_id );
                            push( @cluster_ids, $new_cluster_id );

                            if($verbose) {
                                print "# WARN: remove $gene_stable_id from cluster $cluster_id ($index_dist)\n";
                            }

                            $unclustered{$species}++; 
                        }

                    } else { # gene passes QC
                        push(@checked_ids, $gene_stable_id);
                    }
            }

            # rank genes in terms of chr position
            @checked_ids = sort {
                $id2coords{$species}{$a}[1] <=>
                $id2coords{$species}{$b}[1]
                }
                @checked_ids;

            # updated species gene_ids in main cluster
            $cluster{$main_cluster_id}{$species} = \@checked_ids; 
        }
    }
} 

#foreach $sp2 (@supported_species){ foreach $gene_id (@{ $cluster{'gene:BGIOSGA000010'}{$sp2} }) { print "$gene_id\n" } }

%totaloverlap = (); # not needed anymore

# count how many clusters include each species
foreach $cluster_id (@cluster_ids) {
    foreach $species (@supported_species) {
        if ( $cluster{$cluster_id}{$species} ) {
            $totalclusters{$species}++;
        }
    }
}


# 2.6) add unpaired sequences as singletons 
my $total_seqs = 0;
foreach $species (@supported_species) {

    my $singletons = 0;

    foreach $gene_stable_id ( @{ $sorted_ids{$species} } ) {
        
        next if ( $NOSINGLES || $incluster{ $supported{$species}.':::'.$gene_stable_id } );

        # create new cluster
        $cluster_id = $supported{$species}.':::'.$gene_stable_id;
        $incluster{ $supported{$species}.':::'.$gene_stable_id } = $cluster_id;

        push( @{ $cluster{$cluster_id}{$species} }, $gene_stable_id );
        push( @cluster_ids, $cluster_id );

        # add this singleton to total clusters
        $totalclusters{$species}++;
        $singletons++;
    }

    $total_seqs += $totalgenes{$species};

    printf( "# %s : sequences = %d clusters = %d (unclustered = %d , singletons = %d)\n",
        $species, $totalgenes{$species}, $totalclusters{$species}, 
        $unclustered{$species}, $singletons );
}

printf( "\n# total sequences = %d\n\n", $total_seqs );


# 2.7) create and print shadow clusters of genomic sequences
foreach $cluster_id (@cluster_ids) {

    next if( scalar( keys( %{ $cluster{$cluster_id} } ) ) < $min_taxa);

    if($cluster{$cluster_id}{$ref_genome}) {
        $filename = $cluster{$cluster_id}{$ref_genome}[0]

    } else { 
        $filename = (split(/:::/,$cluster_id))[1];
    }

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

        } else { 
            $filename = (split(/:::/,$cluster_id))[1];
        }

        # write sequences and count sequences
        open( CLUSTER, ">", "$outfolder/$clusterdir/$filename$SEQEXT{$seqtype}" )
          || die "# ERROR: cannot create $outfolder/$clusterdir/$filename$SEQEXT{$seqtype}\n";

        foreach $species (@supported_species) {
            next if ( !$cluster{$cluster_id}{$species} );

            $n_cluster_sp++;
            foreach $gene_stable_id ( @{ $cluster{$cluster_id}{$species} } ) {

                # no sequence printed if not available for this seqtype,
                # frequently you might see cDNA clusters with empty twin CDS clusters
                # for non-protein coding genes 
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

print "\n# percent_conserved_sequences_file = $POCS_matrix_file\n\n";

# sort species from ref down by decreasing POCS
my @supported_species_POCS;

# reference goes in 1st place
push(@supported_species_POCS, $ref_genome);

foreach $sp2 (sort {$POCS2ref{$b}<=>$POCS2ref{$a}} keys(%POCS2ref)) {
    push(@supported_species_POCS, $sp2);
}


## 4)  write pangenome matrices in output folder
if(!$chregex){ # unsorted clusters
    push(@{ $sorted_cluster_ids{'unsorted'} }, @cluster_ids );
}
else { # ordered along homologous chromosomes matching regex

    %sorted_cluster_ids = sort_clusters_by_position( 
        \@supported_species_POCS, \%supported, \%sorted_ids, \%id2coords, 
        $chregex, \%incluster, \%cluster );

    foreach $chr (sort keys(%sorted_cluster_ids)) {
	    printf("# clusters sorted by position in chr %s = %d\n", 
        $chr, scalar(@{ $sorted_cluster_ids{$chr} }));
    }
} 

# sort chromosome names, will be used later on
@sorted_chrs = sort {$a cmp $b} keys(%sorted_cluster_ids);

# set matrix filenames and write headers
my $pangene_matrix_file = "$outfolder/pangene_matrix$params\.tab";
my $pangene_gene_file   = "$outfolder/pangene_matrix_genes$params\.tab";
my $pangene_matrix_tr   = "$outfolder/pangene_matrix$params\.tr.tab";
my $pangene_gene_tr     = "$outfolder/pangene_matrix_genes$params\.tr.tab";
my $pangene_fasta_file  = "$outfolder/pangene_matrix$params\.fasta";
my $pangene_bed_file    = "$outfolder/pangene_matrix$params\.tr.bed";

open( PANGEMATRIX, ">$pangene_matrix_file" )
  || die "# EXIT: cannot create $pangene_matrix_file\n";

open( PANGENEMATRIX, ">$pangene_gene_file" )
  || die "# EXIT: cannot create $pangene_gene_file\n";

print PANGEMATRIX "source:$outfolder/$clusterdir";
foreach $chr (@sorted_chrs) {
    print PANGEMATRIX "\tchr$chr";
    foreach $cluster_id (@{ $sorted_cluster_ids{$chr} }) {

        next if(scalar( keys( %{ $cluster{$cluster_id} } ) ) < $min_taxa);

	$filename = $cluster{$cluster_id}{$ref_genome}[0] || 
          (split(/:::/,$cluster_id))[1];;

        print PANGEMATRIX "\t$filename"; 
    }
}	
print PANGEMATRIX "\n";

print PANGENEMATRIX "source:$outfolder/$clusterdir";
foreach $chr (@sorted_chrs) {
    print PANGENEMATRIX "\tchr:$chr";
    foreach $cluster_id (@{ $sorted_cluster_ids{$chr} }) {

        next if(scalar( keys( %{ $cluster{$cluster_id} } ) ) < $min_taxa);

        $filename = $cluster{$cluster_id}{$ref_genome}[0] || 
          (split(/:::/,$cluster_id))[1];	

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

    foreach $chr (@sorted_chrs) {

        # chr lines have no genes
        print PANGEMATRIX "\tNA";
        print PANGENEMATRIX "\tNA";

        foreach $cluster_id (@{ $sorted_cluster_ids{$chr} }) {

            # skip genes with less occupancy than required
            next if(scalar( keys( %{ $cluster{$cluster_id} } ) ) < $min_taxa);


            if( defined($cluster{$cluster_id}{$species}) && 
                scalar(@{ $cluster{$cluster_id}{$species} }) > 0 ) {

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

# transposed matrices can be more convenient, with genes as rows
system("$TRANSPOSEXE $pangene_matrix_file > $pangene_matrix_tr");
system("$TRANSPOSEXE $pangene_gene_file > $pangene_gene_tr");

print "# pangene_file (occup) = $pangene_matrix_file\n";
print "# pangene_file (occup, transposed) = $pangene_matrix_tr\n";
print "# pangene_file (names) = $pangene_gene_file\n";
print "# pangene_file (names, transposed) = $pangene_gene_tr\n";
#print "# pangene_FASTA_file = $pangene_fasta_file\n";

# 4.1) BED format matrix if clusters are chr-sorted
if($chregex) {

    my ( $start, $end, $strand, $occup, @sorted_ref_gene_ids );

    open( PANGEMATRIBED, ">$pangene_bed_file" )
        || die "# EXIT: cannot create $pangene_bed_file\n";

    foreach $chr (sort keys(%sorted_cluster_ids)) {

        foreach $cluster_id (@{ $sorted_cluster_ids{$chr} }) {

            $filename = $cluster{$cluster_id}{$ref_genome}[0] || 
              (split(/:::/,$cluster_id))[1];

            # compute occupancy (number of genomes where gene is present)
            $occup = 0;
            foreach $species (@supported_species_POCS) {
                if(defined($cluster{$cluster_id}{$species}) &&
                    scalar(@{ $cluster{$cluster_id}{$species} }) > 0 ) {
                    $occup++
                } 
            }            

            # check reference genes, if any
            if( defined($cluster{$cluster_id}{$ref_genome}) &&
                scalar(@{ $cluster{$cluster_id}{$ref_genome} }) > 0 ) {

                @sorted_ref_gene_ids = @{ $cluster{$cluster_id}{$ref_genome} };

            } else { # if none add dummy
                @sorted_ref_gene_ids = ('dummy');
            }

            foreach $gene_stable_id (@sorted_ref_gene_ids) {

                # print each reference gene in a new BED line
                if($gene_stable_id ne 'dummy' ) {

                    if($id2coords{$ref_genome}{$gene_stable_id}) {

                        ( $start, $end, $strand ) =
                            @{$id2coords{$ref_genome}{$gene_stable_id}}[1,2,3];

                        printf(PANGEMATRIBED "%s\t%d\t%d\t%s\t%d\t%s\t%s",
                            $chr,$start,$end,$filename,$occup,$strand,$gene_stable_id);

                    } else { # short reference genes lack coords, cannot be parsed

                        $strand = 0; 
                        printf(PANGEMATRIBED "#%s\tnocoords\tnocoords\t%s\t%d\t%s\t%s",
                            $chr,$filename,$occup,$strand,$gene_stable_id);    
                    }
 
                } else { # gene absent in reference genome annotation

                    $strand = 0; # missing for strand

                    printf(PANGEMATRIBED "#%s\tNA\tNA\t%s\t%d\t%s\tNA",
                        $chr,$filename,$occup,$strand);
                }
         
                # print genes from other species
                foreach $species (@supported_species_POCS) {

                    next if($species eq $ref_genome);

                    if( defined($cluster{$cluster_id}{$species}) &&
                        scalar(@{ $cluster{$cluster_id}{$species} }) > 0 ) {

                        printf(PANGEMATRIBED "\t%s",
                            join( ',', @{ $cluster{$cluster_id}{$species} } ));

                    } else {    # absent genes
                        print PANGEMATRIBED "\tNA";
                    }
                }
 
                print PANGEMATRIBED "\n";
            }
        }
    }

    close(PANGEMATRIBED);

    print "# pangene_file (BED-like) = $pangene_bed_file\n";
}

exit if(!$dogrowth);


## 5) optionally make genome composition analysis to simulate pangene set growth
## as new annotated genomes are added to the pool
## NOTE: this is measured in clusters added/missed per genome, 
## a cluster might contain 1+ gene of same species

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

        # init core- & pan- values
        $coregenome[$s][$sp] = 0;
        $pangenome[$s][$sp]  = $pangenome[$s][ $sp - 1 ];

        # update required occupancy
        $core_occup          = $sp + 1;

        foreach $chr (@sorted_chrs) {
            foreach my $cl (0 .. scalar(@{ $sorted_cluster_ids{$chr} })-1) {

                $cluster_id = $sorted_cluster_ids{$chr}->[$cl];

                # check reference species is in this cluster (1st iteration only)
                if ( $sp == 1 && defined($cluster{$cluster_id}{ $tmptaxa[0] }) && 
                    scalar(@{ $cluster{$cluster_id}{ $tmptaxa[0] } }) > 0 ) {
                    $n_of_taxa_in_cluster{$cluster_id}++
                }

                # check $sp is represented in this cluster
                if ( defined($cluster{$cluster_id}{ $tmptaxa[$sp] }) &&
                    scalar(@{ $cluster{$cluster_id}{ $tmptaxa[$sp] } }) > 0 ) {
                    $n_of_taxa_in_cluster{$cluster_id}++
                }

                # work out cluster occupancy
                if ( $n_of_taxa_in_cluster{$cluster_id} &&
                    defined($cluster{$cluster_id}{ $tmptaxa[$sp] }) &&
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

# Returns a hash with lists of clusters sorted by chr position.
# Only genes in chromosomes, named according to regex, are sorted; 
# remaining genes are added to virtual chr 'unplaced'
# 
# Params:
# i)    ref to list with species, reference in 1st position
# ii)   ref to hash mapping species name to number
# iii)  ref to 2-way hash with a list of chr-sorted gene_ids per species
# iv)   ref to 2-way hash mapping gene_ids to coordinates [ chr, start, end, strand ]
# v)    regex to match chr names; unmatched are added to 'unplaced' virtual chr
# vi)   ref to hash mapping sp:::gene_ids to clusters
# vii)  ref to 2-way hash containing clustered gene_ids
# viii) optional boolean flag to enable verbose output
sub sort_clusters_by_position {

    my ($ref_species, $ref_species_num, $ref_sorted_ids, $ref_id2coords, 
        $regex, $ref_incluster, $ref_cluster, 
        $verbose) = @_;

    my $species_seen = 0;
    my ($species, $sp, $gene, $last_gene, $gene_id, $chr, $cluster_id);
    my ($cluster_idx, $cluster_id2, $chr_name, $sortable_chr, $prev_chr);
    my ($chr2, $gene_id2, $gene2, $last_gene2, $next_cluster_id);
    my $prev_cluster_idx; # index of cluster where previous clustered gene sits
    my $next_cluster_idx; # index of cluster where next clustered gene sits
    my (%cluster_seen, %sorted_cluster_ids);

    foreach $species (@$ref_species) {

        $species_seen++;
        $prev_cluster_idx = -1;
        $prev_chr = '';
        $last_gene = scalar(@{ $ref_sorted_ids->{$species} }) - 1;

        foreach $gene (0 .. $last_gene) {

            ## initialize this gene
            $sortable_chr = 0;
            ($chr_name, $cluster_id) = ('unplaced', '');

            $gene_id = $ref_sorted_ids->{$species}[$gene];
            $chr = $ref_id2coords->{$species}{$gene_id}[0];

            # check chr matches regex, only these chromosomes
            # are guaranteed to be homologous across genomes
            if($chr =~ m/$regex/) {
                $sortable_chr = 1; 
                $chr_name = $chr; 
            }
            
            # new chr
            if($prev_chr ne '' && $chr_name ne $prev_chr) {
                $prev_cluster_idx = -1;
            }


            ## find out which cluster contains this gene
            if(!defined($ref_incluster->{ $ref_species_num->{$species}.':::'.$gene_id })){
                print "# WARNING(sort_clusters_by_position): skip unclustered $gene_id\n" 
                    if($verbose);
                next;
            } else {
                $cluster_id = $ref_incluster->{ $ref_species_num->{$species}.':::'.$gene_id };
            }  

            ## 1st time this cluster was seen (a cluster can only be added to sorted list once)
            if(!$cluster_seen{$cluster_id}) {
				
                # non-reference species, find out where to insert this cluster 
                if($species_seen > 1 && $sortable_chr == 1) {

                    ## 1) index of cluster with next clustered gene on same chr,
                    ##    might be the same for consecutive PAV clusters
                    ($next_cluster_id,$next_cluster_idx) = ('',-1);
                    foreach $gene2 ($gene .. $last_gene) {
                            
                        $gene_id2 = $ref_sorted_ids->{$species}[$gene2];
                        $chr2 = $ref_id2coords->{$species}{$gene_id2}[0];
                        last if($chr2 ne $chr);
 
                        if($ref_incluster->{$gene_id2} && 
                            $cluster_seen{$ref_incluster->{ $ref_species_num->{$species}.':::'.$gene_id2 }}){

                            # cluster containg this sequence
                            $next_cluster_id = $ref_incluster->{ $ref_species_num->{$species}.':::'.$gene_id2 };

                            # index of this cluster in sorted list
                            $next_cluster_idx = _get_element_index( 
                                $sorted_cluster_ids{$chr}, $next_cluster_id);
                            last;
                        }
                    }

                    print "$species $gene_id $cluster_id $chr $prev_cluster_idx ".
                        "$next_cluster_idx $next_cluster_id $species_seen\n" if($verbose); 

                    # 2) actually insert cluster in list 

                    # before previous clusters
                    # prev: NA
                    # next: cluster-------->[0] 
                    if($prev_cluster_idx == -1) {

                        unshift( @{ $sorted_cluster_ids{$chr_name} }, $cluster_id );
                        $cluster_idx = 0
 
                    } elsif($next_cluster_idx == -1) {
                        # after previous clusters -> 
                        # prev: [$#sorted_cluster_ids]-------->cluster
                        # next: NA
							
                        push( @{ $sorted_cluster_ids{$chr_name} }, $cluster_id );
                        $cluster_idx = scalar(@{ $sorted_cluster_ids{$chr_name} })-1;

                    }else {
                        # inserted amid previous clusters
                        # prev: [0..N]----->cluster
                        # next: cluster-----> [N+1..$#sorted_cluster_ids]

                        splice( @{ $sorted_cluster_ids{$chr_name} }, 
                            $prev_cluster_idx+1, 0, $cluster_id );
                        $cluster_idx = $prev_cluster_idx+1
                    }
   
                    # take note this cluster is already processed
                    $cluster_seen{$cluster_id}++;    

                } else {

                    push(@{ $sorted_cluster_ids{$chr_name} }, $cluster_id);
                    $cluster_seen{$cluster_id}++;

                    # debugging info, cluster_id comes from 1st cluster where this
                    # sequence was grouped, which might not be the reference
                    printf("> %s %d %d\n",
                        $cluster_id,scalar(@{ $sorted_cluster_ids{$chr_name} }),$species_seen) 
                        if($verbose);
                }
            } else { # cluster already sorted, can only happen with non-ref species

                if($sortable_chr == 1) {
                    $cluster_idx =
                        _get_element_index( $sorted_cluster_ids{$chr_name}, $cluster_id);
                }
            }

            # save this cluster as previous for next iteration
            $prev_cluster_idx = $cluster_idx;
            $prev_chr = $chr_name;
        }
    }

    return %sorted_cluster_ids;
}

# Takes 3 params:
# i) list ref
# ii) string of element to search
# iii) integer with 0-based starting index, optional
# Return index of elem in list, else -1
sub _get_element_index {

    my ($ref_list, $elem) = @_;

    foreach my $idx (0 .. $#$ref_list) {
        if($ref_list->[$idx] eq $elem) {
            return $idx
        }
    }

    return -1;
}

