package PlantCompUtils;
require Exporter;

# Copyright [2019-2021] EMBL-European Bioinformatics Institute

@ISA       = qw(Exporter);
@EXPORT_OK = qw(
  list_ensembl_mysql_dbs get_canonical_transcript_ids 
  download_FASTA_file parse_isoform_FASTA_file 
  sort_isoforms_chr sort_clusters_by_position
  download_compara_TSV_file download_GTF_file get_gene_coords_GTF_file
  download_MAF_files parse_MAF_file
  perform_rest_action transverse_tree_json
  minimize_MSA write_boxplot_file factorial fisher_yates_shuffle
  @DIVISIONS $REQUEST_COUNT $COMPARADIR $FASTADIR $GTFDIR $MAFDIR
  $FTPURL
);

use strict;
use warnings;
use Net::FTP;
use Time::HiRes;
use HTTP::Tiny;
use DBI;

# Fungi Protists Metazoa have collections and one all-vs-all TSV file
# This code won't work there
our @DIVISIONS  = qw( Plants );
our $FTPURL     = 'ftp.ensemblgenomes.org';
our $COMPARADIR = '/pub/xxx/current/tsv/ensembl-compara/homologies';
our $FASTADIR   = '/pub/current/xxx/fasta';
our $GTFDIR     = '/pub/current/xxx/gtf';
our $MAFDIR     = '/pub/current/xxx/maf/ensembl-compara/multiple_alignments/';

my  $MYSQLURL   = 'mysql-eg-publicsql.ebi.ac.uk';
my  $MYSQLPORT  = 4157; 
my  $MYSQLUSER  = 'anonymous';

my $VOIDVALUE   = -999;

our $REQUEST_COUNT = 0;

my $GZIPEXE = 'gzip';  #default system gzip

# Connects to Ensembl public mysql server,
# returns reference to list of available databases
sub list_ensembl_mysql_dbs {

  my $dbh = DBI->connect("dbi:mysql:host=$MYSQLURL;port=$MYSQLPORT",$MYSQLUSER)
    || die "# ERROR: cannot connect to MySQL server $MYSQLURL port $MYSQLPORT\n";

  my $res = $dbh->selectall_arrayref('SHOW DATABASES');
  my @dbnames = map { $_->[0] } @$res;

  $dbh->disconnect();

  return \@dbnames;
}

## Connects to Ensembl public mysql server, 
## queries species core schema,
## returns reference to hash of canonical transcript stable_ids
sub get_canonical_transcript_ids {

  my ($species_core_db, $verbose) = @_;
  
  my (@tmpdata,%canon_isofs);

  my $dbh = DBI->connect("DBI:mysql:$species_core_db:host=$MYSQLURL;port=$MYSQLPORT",$MYSQLUSER)
    || die "# ERROR: cannot connect to MySQL server $MYSQLURL port $MYSQLPORT\n";

  my $sth = $dbh->prepare( 
    "SELECT t.stable_id FROM gene g, transcript t ".
    "WHERE g.canonical_transcript_id = t.transcript_id"
  );

  $sth->execute() ||
  	die "# ERROR: cannot query $MYSQLURL port $MYSQLPORT\n";

  while(@tmpdata = $sth->fetchrow_array()) {
    $canon_isofs{$tmpdata[0]} = 1;
    if($verbose){ print "$tmpdata[0]\n" }
  }

  $sth->finish();
  $dbh->disconnect();

  return \%canon_isofs;
}

# Parse a FASTA file, either pep or cdna, downloaded with download_FASTA_file
# returns:
# i) a isoform=>sequence hash with the (optionally) selected or
#    (default) longest peptide/transcript per gene
# ii) a isoform=>description hash with chr location & description
#    extracted from the FASTA header
sub parse_isoform_FASTA_file {

    my ( $FASTA_filename, $ref_isoforms2keep ) = @_;

    my ( $stable_id, $gene_stable_id, $max, $len );
    my ( $iso_selected, $len_selected );
    my ( %sequence, %sequence4gene, %header, %header4gene );

    open( FASTA, "$GZIPEXE -dc $FASTA_filename |" )
      || die "# ERROR(parse_isoform_FASTA_file): cannot open $FASTA_filename\n";
    while ( my $line = <FASTA> ) {

        #>g00297.t1 pep supercontig:Ahal2.2:FJVB01000001.1:1390275:1393444:1 gene:g00297 ...
        if ( $line =~ m/^>(\S+).*?gene:(\S+)/ ) {

            $stable_id      = $1;    # might pep or cdna id
            $gene_stable_id = $2;

			# take genomic coordinates
            $header{$stable_id} = ( split( /\s+/, $line ) )[2];

            # add description if possible
            if ( $line =~ m/description:(.*)?\[/ ) {
                $header{$stable_id} .= " $1";
            }
            else {
                $header{$stable_id} .= " no_description ";
            }
        }
        elsif ( $line =~ m/^(\S+)/ ) {
            $sequence{$gene_stable_id}{$stable_id} .= $1;
        }
    }
    close(FASTA);

    foreach $gene_stable_id ( keys(%sequence) ) {

        # work out which isoform should be kept for this gene
        ( $max, $iso_selected, $len_selected ) = ( 0, '', '' );
        foreach $stable_id ( keys( %{ $sequence{$gene_stable_id} } ) ) {

            # find longest isoform (default), note that key order is random
            $len = length( $sequence{$gene_stable_id}{$stable_id} );
            if ( $len > $max ) {
                $max          = $len;
                $len_selected = $stable_id;
            }

            if ( $ref_isoforms2keep->{$stable_id} ) {
                $iso_selected = $stable_id;
            }
        }

        if ($iso_selected) {
            $sequence4gene{$iso_selected} =
              $sequence{$gene_stable_id}{$iso_selected};
            $header4gene{$iso_selected} = $header{$iso_selected};
        }
        elsif ($len_selected) {
            $sequence4gene{$len_selected} =
              $sequence{$gene_stable_id}{$len_selected};
            $header4gene{$len_selected} = $header{$len_selected};
        }
        else {
            print "# ERROR(parse_isoform_FASTA_file): cannot select ".
                "an isoform for gene $gene_stable_id\n";
        }
    }

    return ( \%sequence4gene, \%header4gene );
}

# reads a compressed GTF file and returns a reference to an array with
# chr-order geneids and their coordinates
sub get_gene_coords_GTF_file {

    my ($GTF_filename) = @_;

    my ( $chr, $start, $end, $strand, $geneid );
    my (@chr_sorted_gene_ids);

    open( GTF, "$GZIPEXE -dc $GTF_filename |" )
      || die "# ERROR(get_gene_coords_GTF_file): cannot open $GTF_filename\n";
    while ( my $line = <GTF> ) {

        #1 araport11 gene 3631 5899 . + . gene_id "AT1G01010";...
        if ( $line =~
               m/^([^#])\t[^\t]+\tgene\t(\d+)\t(\d+)\t[^\t]\t(\S+)\t[^\t]\tgene_id "([^";]+)/
          )
        {
            ( $chr, $start, $end, $strand, $geneid ) = ( $1, $2, $3, $4, $5 );
            push( @chr_sorted_gene_ids,
                [ $geneid, $chr, $start, $end, $strand ] );
        }
    }
    close(GTF);

    return \@chr_sorted_gene_ids;
}

# download compressed GTF file from FTP site, renames it
# and saves it in $targetdir; uses FTP globals defined above
sub download_GTF_file {

    my ( $dir, $ref_genome, $targetdir ) = @_;
    my ( $gtf_file, $stored_gtf_file ) = ( '', '' );

    if ( my $ftp =
        Net::FTP->new( $FTPURL, Passive => 1, Debug => 0, Timeout => 60 ) )
    {
        $ftp->login( "anonymous", '-anonymous@' )
          || die "# ERROR(download_GTF_file): cannot login " . $ftp->message();
        $ftp->cwd($dir)
          || die
          "# ERROR(download_GTF_file): cannot change working directory to $dir "
          . $ftp->message();
        $ftp->cwd($ref_genome)
          || die "# ERROR(download_GTF_file): cannot find $ref_genome in $dir "
          . $ftp->message();

        # find out which file is to be downloaded and
        # work out its final name
        foreach my $file ( $ftp->ls() ) {
            if ( $file =~ m/gtf.gz/ && $file !~ m/.chr.gtf/ ) {
                $gtf_file        = $file;
                $stored_gtf_file = "$targetdir/$gtf_file";
                last;
            }
        }

        unless ( -s $stored_gtf_file ) {
            $ftp->binary();
            my $downsize = $ftp->size($gtf_file);
            $ftp->hash( \*STDOUT, $downsize / 20 ) if ($downsize);
            printf( "# downloading %s (%1.1fMb) ...\n",
                $gtf_file, $downsize / ( 1024 * 1024 ) );
            print "# [        50%       ]\n# ";
            if ( !$ftp->get($gtf_file) ) {
                die
                  "# ERROR(download_GTF_file): failed downloading $gtf_file\n";
            }

            # rename file to final name
            rename( $gtf_file, $stored_gtf_file );

            print "# using $gtf_file\n\n";
        }
        else {
            print "# re-using $stored_gtf_file\n\n";
        }
    }
    else {
        die "# ERROR(download_GTF_file): cannot connect to $FTPURL , please try later\n";
    }

    return $stored_gtf_file;
}

# download compressed TSV file from FTP site, renames it
# and saves it in $targetdir; uses FTP globals defined above
# NOTE: if species file is not found it tries the bulky all-vs-all file
sub download_compara_TSV_file {

    my ( $dir, $ref_genome, $targetdir ) = @_;
    my ( $compara_file, $stored_compara_file ) = ( '', '' );

    if ( my $ftp =
        Net::FTP->new( $FTPURL, Passive => 1, Debug => 0, Timeout => 60 ) )
    {
        $ftp->login( "anonymous", '-anonymous@' )
          || die "# ERROR(download_compara_TSV_file): cannot login "
          . $ftp->message();
        $ftp->cwd($dir)
          || die "# ERROR(download_compara_TSV_file): cannot change working directory to $dir "
          . $ftp->message();

        # find out which file is to be downloaded
        if ( $ftp->cwd($ref_genome) ) {
            foreach my $file ( $ftp->ls() ) {
                if ( $file =~ m/protein_default.homologies.tsv.gz/ ) {
                    $compara_file        = $file;
                    $stored_compara_file = "$targetdir/$compara_file";
                    $stored_compara_file =~ s/tsv.gz/$ref_genome.tsv.gz/;
                    last;
                }
            }
        }
        else {    # try all-vs-all file instead (Fungi, Protists, Metazoa)

            print "# WARNING(download_compara_TSV_file): cannot find ".
			    "$ref_genome in $dir, try all-vs-all\n";

            foreach my $file ( $ftp->ls() ) {
                if ( $file =~ m/protein_default.homologies.tsv.gz/ ) {
                    $compara_file        = $file;
                    $stored_compara_file = "$targetdir/$compara_file";
                    foreach my $div (@DIVISIONS) {
                        if ( $dir =~ m/($div)/i ) {
                            $div = $1;
                            $stored_compara_file =~ s/tsv.gz/$div.tsv.gz/;
                            last;
                        }
                    }
                }
            }
        }

        if ( !$compara_file ) {
            die "# ERROR(download_compara_TSV_file): cannot find file to download\n";
        }

        # download that TSV file
        unless ( -s $stored_compara_file ) {
            $ftp->binary();
            my $downsize = $ftp->size($compara_file);
            $ftp->hash( \*STDOUT, $downsize / 20 ) if ($downsize);
            printf( "# downloading %s (%1.1fMb) ...\n",
                $stored_compara_file, $downsize / ( 1024 * 1024 ) );
            print "# [        50%       ]\n# ";
            if ( !$ftp->get($compara_file) ) {
                die "# ERROR(download_compara_TSV_file): failed downloading $compara_file\n";
            }

            # rename file to final name
            rename( $compara_file, $stored_compara_file );

            print "# using $stored_compara_file\n\n";
        }
        else {
            print "# re-using $stored_compara_file\n\n";
        }
    }
    else {
        die "# ERROR(download_compara_TSV_file): cannot connect to ".
            "$FTPURL , please try later\n";
    }

    return $stored_compara_file;
}

# download compressed FASTA file from FTP site, and saves it in $targetdir;
# uses FTP globals defined above
sub download_FASTA_file {

    my ( $dir, $genome_folder, $targetdir ) = @_;
    my ( $fasta_file, $stored_fasta_file ) = ( '', '' );

    if ( my $ftp =
        Net::FTP->new( $FTPURL, Passive => 1, Debug => 0, Timeout => 60 ) )
    {
        $ftp->login( "anonymous", '-anonymous@' )
          || die "# ERROR(download_FASTA_file): cannot login "
          . $ftp->message();
        $ftp->cwd($dir)
          || die "# ERROR(download_FASTA_file): cannot change working directory to $dir "
          . $ftp->message();
        $ftp->cwd($genome_folder)
          || die "# ERROR(download_FASTA_file): cannot find $genome_folder in $dir "
          . $ftp->message();

        # find out which file is to be downloaded and
        # work out its final name
        foreach my $file ( $ftp->ls() ) {
            if ( $file =~ m/all.fa.gz/ ) {
                $fasta_file        = $file;
                $stored_fasta_file = "$targetdir/$fasta_file";
                last;
            }
        }

        # download that FASTA file
        unless ( -s $stored_fasta_file ) {
            $ftp->binary();
            my $downsize = $ftp->size($fasta_file);
            $ftp->hash( \*STDOUT, $downsize / 20 ) if ($downsize);
            printf( "# downloading %s (%1.1fMb) ...\n",
                $fasta_file, $downsize / ( 1024 * 1024 ) );
            print "# [        50%       ]\n# ";
            if ( !$ftp->get($fasta_file) ) {
                die "# ERROR(download_FASTA_file): failed downloading $fasta_file\n";
            }

            # rename file to final name
            rename( $fasta_file, $stored_fasta_file );

            print "# using $fasta_file\n\n";
        }
        else {
            print "# re-using $stored_fasta_file\n\n";
        }
    }
    else {
        die "# ERROR(download_FASTA_file): cannot connect to ".
             "$FTPURL , please try later\n";
    }

    return $stored_fasta_file;
}

# Download compressed MAF files from FTP site, and saves them in $targetdir.
# Alignments are grouped by chromosome, and then by coordinate system.
# The files named *.other*.maf contain alignments that do not include any reference
# region. Each file contains up to 200 alignments.
# uses FTP globals defined above
sub download_MAF_files {

    my ( $dir, $suffix, $targetdir, $verbose ) = @_;
    my ( $maf_file, $stored_maf_file, @stored_files ) = ( '', '' );

    if ( my $ftp =
           Net::FTP->new( $FTPURL, Passive => 1, Debug => 0, Timeout => 60 ) )
    {
        $ftp->login( "anonymous", '-anonymous@' )
          || die "# ERROR(download_FASTA_file): cannot login "
          . $ftp->message();
        $ftp->cwd($dir)
          || die "# ERROR(download_FASTA_file): cannot change working directory to $dir "
          . $ftp->message();

        # find out which files are to be downloaded
        foreach my $file ( $ftp->ls() ) {
            if ( $file =~ m/^$suffix.*.maf.gz/ ) {
                $maf_file = $file;
                $stored_maf_file = "$targetdir/$maf_file";

                unless ( -s $stored_maf_file ) {
                    $ftp->binary();
                    # don't print progress as there are too many files
                    #my $downsize = $ftp->size($maf_file);
                    #$ftp->hash( \*STDOUT, $downsize / 20 ) if ($downsize);
                    #printf( "# downloading %s (%1.1fMb) ...\n",
                    #$maf_file, $downsize / ( 1024 * 1024 ) );
                    #print "# [        50%       ]\n# ";
                    if ( !$ftp->get($maf_file) ) {
                        die "# ERROR(download_MAF_files): failed downloading $maf_file\n";
                    }
                    # rename file to final name
                    rename( $maf_file, $stored_maf_file );
                    print "# using $maf_file\n" if($verbose);
                }
                else {
				    print "# re-using $stored_maf_file\n" if($verbose);
               }

               push(@stored_files, $stored_maf_file);
			}
		}
	}
	else {
        die "# ERROR(download_MAF_files): cannot connect to ".
            "$FTPURL , please try later\n";
    }

	return @stored_files;
}	

# Takes a ref to a list of Ensembl production names and 
# returns a hash with abbreviated names as keys (as those in MAF
# Newick trees) and production names as values
sub _shorten_production_name {
    my ($ref_long_names) = @_;

    my (%abbrev,$short);
    foreach my $name (@$ref_long_names){
        if($name =~ m/([^_]+)_([^_]+)_([^_]+)/) { # trinomial
            $short = uc(substr($1,0,1)) . substr($2,0,1) . substr($3,0,2);
        } elsif($name =~ m/([^_]+)_([^_]+)/) { # binomial
            $short = uc(substr($1,0,1)) . substr($2,0,3);
	    }		
        $abbrev{ $short } = $name;
	}

    return \%abbrev;
}

# Takes 7 args:
# i)   path to GZIP-compressed MAF file (readonly)
# ii)  ref to list of target species, only those will be parsed (readonly)
# iii) ref to hash with BED filenames  (readonly)
# iv)  path to bedtools
# v)   fraction of gene overlap [0-1] 
# vi)  ref to hash with blocks to add parsed data
# vii) ref to hash mapping isoforms to blocks
# viii)boolean to enable verbose mode
sub parse_MAF_file {

    my ($maf_filename, $ref_species, $ref_bedfiles, 
        $bedexe, $gene_overlap, $ref_blocks, 
        $ref_iso2block, $verbose) = @_;

    my ($sp, $chr, $start, $end, $strand);
    my ($species, $command, $newick);

    # shorten species names, to match those in Newick string
    my $ref_shortnames = _shorten_production_name($ref_species);

    open( MAF, "$GZIPEXE -dc $maf_filename |" )
      || die "# ERROR(parse_MAF_file): cannot open $maf_filename\n";
    while ( my $line = <MAF> ) {

        # tree: ((Oind_4_11238666_11248253[+]:0.0061515,Oniv_4_7962011_7971991[+]:...
        if($line =~ /^# tree: (\S+)/){

            # we will use Newick tree to identify this block
            $newick = $1;

            while($newick =~ m/(\w{4})_([\w_]+?)_(\d+)_(\d+)\[([+-])\]:/g) {

                #Osat_5_7890_14004[+]:1e-05
                #Oruf_5_7624_13133[+]:0.001504
                #Osat_Syng_TIGR_005_8011_14124[-]:0.000799
                #Aseq_Ancestor_9910_23859_1_5502[+]:4.8e-05
	            			   
                ($sp, $chr, $start, $end, $strand) = ($1,$2,$3,$4,$5);

                # skip unselected species
				next if not defined($ref_shortnames->{$sp});

                $species = $ref_shortnames->{$sp};

                # compute block length per species
                $ref_blocks->{$newick}{$species}{'length'} = ($end-$start)+1;

                # block occupancy
                $ref_blocks->{$newick}{'ocuppancy'}++;

				# find genes within block, 
				# should work only for scaffolds containing genes
                if($ref_bedfiles->{$species}{$chr}) {

                    $command = "echo -e \"$chr\t$start\t$end\" | " .
                        "$bedexe intersect -sorted -f $gene_overlap " .
                        "-a $ref_bedfiles->{$species}{$chr} -b stdin";
                    open(BEDTOOLS,"$command |") 
                      || die "# ERROR(parse_MAF_file): cannot run $command\n";
                    while(<BEDTOOLS>) {
                        #10	2645232	2645660	ONIVA10G02410.1	0	-
						my @bedata = split;
                        $ref_iso2block->{$species}{$bedata[3]} = $newick;
                        push( @{ $ref_blocks->{$newick}{$species}{'genes'} }, $bedata[3] );
                    }
                    close(BEDTOOLS);
                } elsif(1 || $verbose) {
                    print "# WARNING(parse_MAF_file): cannot find BED file for $species chr $chr\n";
                }
            }
        }
    }
    close(MAF); 
}

# Takes 4 args:
# i)   hash ref of headers produced by parse_isoform_FASTA_file
# ii)  folder to place BED files
# iii) species production_name
# iv)  string with regex to match chr names
# Produces 1-based inclusive BED file with canonical isoforms, one per chr/scaffold, 
# sorted by start coordinate. These can be directly compared to intervals in MAF files
# Returns:
# i)   hash ref with chr as key and BED filename as value
# ii)  hash ref with chr-sorted isoform names, one per chr
# iii) hash ref matching isoform ids to chrs (those in scaffols are excluded)
sub sort_isoforms_chr {
    my ($ref_header, $bedfolder, $species, $regex) = @_;
    my ($stable_id,$start,$end);
    my ($chr,$strand,$bedfile);
    my (%raw,%bedfiles,%sorted_ids,%id2chr);

    for $stable_id (keys(%$ref_header)){
        # chromosome:IRGSP-1.0:12:8823315:8825166:-1 , Ensembl 1-based inclusive
        # scaffold:v.1.0:scaffold_9:1012538:1013882:-1
		if($ref_header->{$stable_id} =~ m/[^:]+:[^:]+:([^:]+):([^:]+):([^:]+):([\d-]+)/) {
            ($chr,$start,$end,$strand) = ($1,$2,$3,$4); 

            push(@{ $raw{$chr} }, [$start,$end,$stable_id,$strand] );
            
            if($chr =~ m/$regex/) {
                $id2chr{ $stable_id } = $chr;
            }
        }
    } 

    # sort isoforms along chr/scaffolds
	# Note: sometimes genes on the same strand can overlap, see
    # http://plants.ensembl.org/Oryza_sativa/Gene/Summary?db=core;g=Os06g0168150;r=6:3426914-3434445
	foreach $chr (keys(%raw)) {
       
        # create new BED file for this chr
        $bedfile = "$bedfolder/$species.$chr.bed";
		open(BED,">",$bedfile) || die "# ERROR(sort_isoforms_chr): cannot create $bedfile\n";

        foreach my $isof (sort {$a->[0] <=> $b->[0]} @{ $raw{$chr} }) {
            printf(BED "%s\t%d\t%d\t%s\t0\t%s\n",
                $chr,$isof->[0],$isof->[1],$isof->[2],
                $isof->[3] == -1 ? '-' : '+');
            push(@{ $sorted_ids{$chr} }, $isof->[2]);
        }

		close(BED);

        $bedfiles{$chr} = $bedfile;
    }
    
    return (\%bedfiles, \%sorted_ids, \%id2chr);
}

# Returns a hash with lists of clusters sorted by chr position.
# Only genes in chromosomes, named according to regex, are sorted; all
# other genes are added to chr 'unplaced', and clusters containing genes from
# different chromosomes added to 'mixed' virtual chromosome
# Takes 5(6) parms:
# i)   ref to list with species production names
# ii)  ref to 2-way hash with chr-sorted lists of genes
# iii) ref to hash mapping gene isoform ids to clusters
# iv)  ref to 2-way hash containing clustered isoforms
# v)   ref to hash mapping gene isoform ids to chr
# vi)  regex to match chr names
# vii) optional boolean flag to enable verbose output
sub sort_clusters_by_position {
    my ($ref_supported_species, $ref_sorted_ids, 
        $ref_incluster, $ref_cluster, $ref_id2chr, $regex, $verbose) = @_;

    my ($species_seen, $isof, $isof2, $chr, $chr2, $sp, $cluster_id) = ( 0 );
    my ($ref_sorted_chr_ids,$last_isof,$next_cluster_id, $is_mixed);
    my ($isof_stable_id, $cluster_idx, $isof_idx, $sortable_chr, $chr_name);
    my $prev_cluster_idx; # index of cluster where previous clustered $isof sits
    my $next_cluster_idx; # index of cluster where next clustered $isof sits
    my (%cluster_seen, %sorted_cluster_ids);
    my $n_mixed = 0; # debug

    foreach my $species (@$ref_supported_species) {

        $species_seen++;
        
        my @chrs = keys(%{ $ref_sorted_ids->{$species} });
		#@chrs = qw( 10 ); # debug

        foreach $chr (@chrs) {

            # init 
            $prev_cluster_idx = $VOIDVALUE;
            $sortable_chr = 0;
            $chr_name = 'unplaced';

            if($chr =~ m/$regex/) {
                $sortable_chr = 1; 
                $chr_name = $chr;
            }

            $ref_sorted_chr_ids = $ref_sorted_ids->{$species}{$chr};
            $last_isof = scalar(@$ref_sorted_chr_ids)-1;

            for $isof (0 .. $last_isof) {

                # init 
                $cluster_id = '';
                $is_mixed = 0;
                $isof_stable_id = $ref_sorted_chr_ids->[$isof];

                # find out which cluster contains this isoform  
                if(!defined($ref_incluster->{$isof_stable_id})){
                    print "# WARNING(sort_clusters_by_position): skip unclustered ".
                        "$isof_stable_id\n" if($verbose);
                    next;
                } else {
                    $cluster_id = $ref_incluster->{$isof_stable_id};
                }

                # first time this cluster was seen (a cluster can only be added once)
                if(!$cluster_seen{$cluster_id}) {
				
                    # does cluster contain seqs from different chrs?
                    my %cluster_chrs;
                    foreach $sp (keys(%{ $ref_cluster->{$cluster_id} })) {
                        foreach $isof2 (@{ $ref_cluster->{$cluster_id}{$sp} }) { 
                            if($ref_id2chr->{$sp}{$isof2}) { 
                                # only genes in chr matching regex count
                                $chr2 = $ref_id2chr->{$sp}{$isof2};
                                print "$cluster_id $sp $isof2 $chr2\n" if($verbose);
                                $cluster_chrs{$chr2}++;
                            }
                        } 
                    }
                        
                    # add it virtual chr 'mixed' if so
                    if(scalar(keys(%cluster_chrs)) > 1) {
                        $n_mixed++;
                        $chr_name = 'mixed';
                        $is_mixed = 1;
                        foreach $chr2 (keys(%cluster_chrs)){ print "$cluster_id $chr2 $n_mixed\n"; }
				    }	

                    # non-reference species, find out where to insert this cluster 
                    if($species_seen > 1 && $sortable_chr == 1 && $is_mixed == 0) {

                        # 1) get index of cluster harborig next clustered isoform ($isof_stable_id)
                        $next_cluster_idx = $VOIDVALUE;
                        foreach $isof2 ($isof .. $last_isof) {
                            
                            if($ref_incluster->{$ref_sorted_chr_ids->[$isof2]} && 
                                   $cluster_seen{$ref_incluster->{$ref_sorted_chr_ids->[$isof2]}}){

                                # cluster containg this sequence
                                $next_cluster_id = $ref_incluster->{$ref_sorted_chr_ids->[$isof2]};

                                # index of this cluster in sorted list
                                $next_cluster_idx = _get_element_index( $sorted_cluster_ids{$chr}, $next_cluster_id);
                                last;
                            }
                        }

                        print "$species $isof_stable_id $cluster_id $chr $prev_cluster_idx ".
                            "$next_cluster_idx $next_cluster_id $ref_sorted_chr_ids->[$last_isof]\n" if($verbose); 

                        # 2) actually insert cluster in list 

                        # before previous clusters
                        # prev: NA
                        # next: cluster-------->[0] 
                        if($prev_cluster_idx == $VOIDVALUE) {

                            unshift( @{ $sorted_cluster_ids{$chr_name} }, $cluster_id );
                            $cluster_idx = 0 
                        } 
                        elsif($next_cluster_idx == $VOIDVALUE) {
                            # after previous clusters -> 
                            # prev: [$#sorted_cluster_ids]-------->cluster
                            # next: NA
							
							push( @{ $sorted_cluster_ids{$chr_name} }, $cluster_id );
                            $cluster_idx = scalar(@{ $sorted_cluster_ids{$chr_name} })-1;
                        }
                        else {
                            # inserted amid previous clusters
                            # prev: [0..N]----->cluster
                            # next: cluster-----> [N+1..$#sorted_cluster_ids]
					
                            splice( @{ $sorted_cluster_ids{$chr_name} }, 
                                $prev_cluster_idx+1, 0, $cluster_id );
                            $cluster_idx = $prev_cluster_idx+1
                        }

                    } else { 
                        push(@{ $sorted_cluster_ids{$chr_name} }, $cluster_id);
                        $cluster_seen{$cluster_id}++;

						printf("> %s %d\n",
                            $cluster_id,scalar(@{ $sorted_cluster_ids{$chr_name} })) if($verbose);
                    }

                } else {
                    # cluster already sorted, can only happen with non-ref species
                    if($sortable_chr == 1) {
                        $cluster_idx =
                            _get_element_index( $sorted_cluster_ids{$chr_name}, $cluster_id);
                    }
                }

				# save this cluster as previous for next iteration
                $prev_cluster_idx = $cluster_idx;
            }
        }
    }

    return %sorted_cluster_ids;
}

# Takes 3 params:
# i) list ref
# ii) string of element to search
# iii) integer with 0-based starting index, optional
# Return index of elem in list, else $VOIDVALUE
sub _get_element_index {

    my ($ref_list, $elem) = @_;
	
    foreach my $idx (0 .. $#$ref_list) {
        if($ref_list->[$idx] eq $elem) {
            return $idx
        }
    }

    return $VOIDVALUE;
}

# uses global $REQUEST_COUNT
# takes a HTTP::Tiny object as first param
# based on examples at 
# https://github.com/Ensembl/ensembl-rest/wiki/Example-Perl-Client
sub perform_rest_action {
    my ( $http, $url, $headers ) = @_;

    $headers ||= {};
    $headers->{'Content-Type'} = 'application/json'
      unless exists $headers->{'Content-Type'};

    my ( $diff, $last_request_time, $current_time ) = ( 0, 0 );

    if ( $REQUEST_COUNT == 15 ) {    # check every 15
        $current_time = Time::HiRes::time();
        $diff         = $current_time - $last_request_time;

        # if less than a second then sleep for the remainder of the second
        if ( $diff < 1 ) {
            Time::HiRes::sleep( 1 - $diff );
        }

        # reset
        $last_request_time = Time::HiRes::time();
        $REQUEST_COUNT     = 0;
    }

    my $response = $http->get( $url, { headers => $headers } );
    my $status = $response->{status};

    if ( !$response->{success} ) {

     # check for rate limit exceeded & Retry-After (lowercase due to our client)
        if ( ( $status == 429 || $status == 599 )
            && exists $response->{headers}->{'retry-after'} )
        {
            my $retry = $response->{headers}->{'retry-after'};
            Time::HiRes::sleep($retry);

            # after sleeping see that we re-request
            return perform_rest_action( $url, $headers );
        }
        else {
            my ( $status, $reason ) =
              ( $response->{status}, $response->{reason} );
            die "# ERROR(perform_rest_action): failed REST request $url\n".
                 "# Status code: ${status}\n# Reason: ${reason}\n# Please re-run";
        }
    }

    $REQUEST_COUNT++;

    if   ( length( $response->{content} ) ) { return $response->{content} }
    else                                    { return '' }
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

sub factorial {
    my $max = int( $_[0] );
    my $f   = 1;
    for ( 2 .. $max ) { $f *= $_ }
    return $f;
}

# generates a random permutation of @array in place
sub fisher_yates_shuffle {
    my $array = shift;
    my ( $i, $j );

    for ( $i = @$array ; --$i ; ) {
        $j = int( rand( $i + 1 ) );
        next if $i == $j;
        @$array[ $i, $j ] = @$array[ $j, $i ];
    }

    return join( '', @$array );
}

# takes 
# i) decoded JSON->{'tree'}, as returned from 
#          http://rest.ensembl.org/genetree/member/id/ and
# ii) reference to hash which gets populated with { accession, sequence } pairs
sub transverse_tree_json {

    my ( $tree, $ref_hash ) = @_;

    my ( $child, $id );

    if ( $tree->{'sequence'} ) {
        foreach $id ( @{ $tree->{'sequence'}{'id'} } ) {
            $ref_hash->{ $id->{'accession'} } =
              $tree->{'sequence'}{'mol_seq'}{'seq'};
        }
    }

    foreach $child ( @{ $tree->{'children'} } ) {
        transverse_tree_json( $child, $ref_hash );
    }
}

# takes a two-key hash with aligned sequences
# and returns a hash ref with minimized aligned seqs,
# after removing gap-only columns
sub minimize_MSA {

    my ( $ref_hash ) = @_;

    my $n_of_seqs = scalar( keys( %$ref_hash ));
    my ( %gap, %minMSA, $key1, $key2, $pos );

    # save position of observed gaps
    foreach $key1 (keys( %$ref_hash )) {
        foreach $key2 (keys( %{$ref_hash->{$key1}} )) {
            while($ref_hash->{$key1}{$key2} =~ m/-/g) {
                $gap{ pos($ref_hash->{$key1}{$key2})-1 }++;
            } 
        }
    }

    # copy sequences to output array but skip gap-only columns
    foreach $key1 (keys( %$ref_hash )) {
        foreach $key2 (keys( %{$ref_hash->{$key1}} )) {
            foreach $pos (0 .. length($ref_hash->{$key1}{$key2})-1) {
                next if($gap{ $pos } && $gap{ $pos } == $n_of_seqs);
                $minMSA{$key1}{$key2} .= substr($ref_hash->{$key1}{$key2},$pos,1); 
            }
        }
    }

    return \%minMSA;
}

1;
