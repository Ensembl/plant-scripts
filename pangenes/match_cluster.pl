#!/usr/bin/env perl

# This script matches input nucleotide sequences to clusters produced by get_pangenes.pl .
# It creates a temporary sequence index with nucleotide sequences from clusters and 
# uses GMAP to scan them. Supports cdna [default] or CDS sequences.
# Optionally, sequence index can be exported as gene-based pangenome for mapping,
# with <global pangenome positions> estimated from reference genome.
#
# Note: -F looks for ../INP_dir/_reference_genome.fna.fai (ie _pangenes folder)

# Copyright [2023-25]
# EMBL-European Bioinformatics Institute & Estacion Experimental Aula Dei-CSIC

$|=1;

use strict;
use warnings;
use Getopt::Std;
use Cwd qw/ abs_path /;
use File::Temp qw/ tempfile /;
use File::Basename;
use File::Copy;
use FindBin '$Bin';
use lib "$Bin/lib";
use pangeneTools qw( check_installed_features 
                     feature_is_installed 
                     parse_sequence_FASTA_file );

my @FEATURES2CHECK = (
  'EXE_GMAP'
);

my $GMAPBIN = $ENV{'EXE_GMAP'};
my $GMAPBUILDBIN = $ENV{'EXE_GMAP_BUILD'};

my ($INP_dir, $INP_seqfile, $INP_isCDS, $INP_idfile) = ('','',0,'');
my ($INP_outfile,$INP_threads, $INP_ow, $INP_verbose) = ('',4, 0, 0);
my ($INP_make_FASTA_reference,$INP_minident) = (0, 95.0);
my ($pangene_bed_file,$pangene_ref_file,$neigh) = ('','',0);
my ($cluster_regex, $isCDS, $cluster_folder ) = ('.cdna.fna', 0, '');
my ($gmapdb_path, $gmapdb, $seq, $regex, $gcoords);
my ($max_number_isoforms, $qlen, $tlen, $identity, $cover) = (0);
my ($cluster_file, $gene_id, $isof_id, $cmd, $index_file);
my ($cluster_id, $seq_id, $cluster_seq_id, $coords, $taxon);
my (%opts, %matches, @order_geneid, %global_coords, @pangene_files);

getopts('hFvcCwI:d:s:t:o:i:', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0)) {
  print "\nusage: $0 [options]\n\n";
  print "-h this message\n";
  print "-c print credits and checks installation\n"; 
  print "-d directory produced by get_pangenes.pl          (example: -d /path/data_pangenes/..._algMmap_ ,\n";
  print "                                                   should include folder with [cdna|cds].fna files)\n";
  print "-s nucleotide sequence file in FASTA format       (example: -s transcripts.fna,\n";
  print "                                                   helps if header has genomic coords ie chr1:12-1200)\n";
  print "-o output file in TSV format                      (required with -s)\n";
  print "-F make pangene reference FASTA format            (optional, pangenome coordinates estimated from reference,\n";
  print "                                                   requires -d 0taxa.._split dir, obtained with get_pangenes.pl -s -t 0)\n";
  print "-C use CDS sequences                              (optional, cDNA sequences are scanned by default)\n";
  print "-i min % sequence identity                        (optional, default: -i $INP_minident)\n";
  print "-t threads                                        (optional, default: -t $INP_threads)\n";
  print "-w overwrite GMAP index                           (optional, by default index is re-used if possible)\n";
  print "-v verbose                                        (optional)\n";
  #print "-I TSV file matching cluster names to pangene ids (optional, example: -I cluster2id.tsv)\n\n";
  exit(0);
}

if(defined($opts{'c'})) {
  print "\nPrimary citation:\n https://doi.org/10.1186/s13059-023-03071-z\n";
  print "\nThis software uses external algorithms, please cite them accordingly:\n";
  print " gmap https://doi.org/10.1093/bioinformatics/bti310\n";

  # check all binaries needed by this program and print diagnostic info
  print check_installed_features(@FEATURES2CHECK);
  exit(0);
}

if(defined($opts{'d'})) { 
  $INP_dir = $opts{'d'};

  if(defined($opts{'C'})){
    $INP_isCDS = 1;
    $cluster_regex = '.cds.fna';
  }

  if(defined($opts{'F'})){

    $pangene_bed_file = $INP_dir.'/pangene_matrix.tr.bed';
    $pangene_ref_file = $INP_dir.'/pangene_matrix.reference.fna';

    if($INP_dir =~ m/_0taxa_/ && $INP_dir =~ m/_split_/ && -e $pangene_bed_file) {
      $INP_make_FASTA_reference = 1;
      
      if($INP_dir =~ m/_(\d+)neigh_/) { # store neighbor distance used to compute pangenes
        $neigh = $1;
      }    

    } else {
      die "# EXIT : -F requires a -d directory with _split_ in name, obtained with get_pangenes.pl -s\n";
    }
  }    
}
else{ die "# EXIT : need a -d directory\n" }

if(defined($opts{'s'})){  
  $INP_seqfile = $opts{'s'};

  if(!-e $INP_seqfile) {
    die "# EXIT : need a valid -s file\n"
  }
} 
elsif($INP_make_FASTA_reference == 0) {  
  die "# EXIT : need parameter -s or -F\n" 
}

if(defined($opts{'o'})){
  $INP_outfile = $opts{'o'};
} 
elsif($INP_make_FASTA_reference == 0) { 
  die "# EXIT : need parameter -o\n" 
}

if(defined($opts{'I'})){ 
  $INP_idfile = $opts{'I'};  
}

if(defined($opts{'w'})){
  $INP_ow = 1;
}

if(defined($opts{'v'})){
  $INP_verbose = 1;
}

if(defined($opts{'t'}) && $opts{'t'} >= 1){
  $INP_threads = int($opts{'t'});
}

if(defined($opts{'i'}) && $opts{'i'} >= 0 && $opts{'i'} <= 100){
  $INP_minident = int($opts{'i'});
}


print "# $0 -d $INP_dir -s $INP_seqfile -o $INP_outfile -C $INP_isCDS " .
  "-w $INP_ow -t $INP_threads -i $INP_minident -F $INP_make_FASTA_reference -v $INP_verbose\n\n";


# 0) locate .cluster_list file and guess actual folder with clusters
opendir(INPDIR,$INP_dir) ||
  die "# ERROR: cannot list $INP_dir , please check -d argument is a valid folder\n";
my @dirfiles = grep {/\.cluster_list/} readdir(INPDIR);
closedir(INPDIR);

if(@dirfiles) {
  $cluster_folder = $INP_dir. '/' . (split(/\.cluster_list/,$dirfiles[0]))[0];

  $gmapdb_path = dirname($cluster_folder);  
  $gmapdb = basename($cluster_folder); 
  if($INP_make_FASTA_reference) {
    $gmapdb .= '.reference';  
  }
  if($INP_isCDS == 1) {
    $gmapdb .= '.cds';
  }	  
  $gmapdb .= '.gmap';

} else {
  die "# ERROR: cannot locate folder with clusters within $INP_dir\n";
}

print "# cluster folder: $cluster_folder\n";

# 1) parse pre-computed cluster files, write sequences to temp file & make GMAP index

# 1.1) by default parse FASTA files from pangene folder, no pangenome/graph global coords.
# Cluster files are not sorted, we noticed this can sometimes change gmap results 
if($INP_make_FASTA_reference == 0) {

  opendir(CLUSTDIR,$cluster_folder) || 
    die "# ERROR: cannot list $cluster_folder , please check -d argument is a valid folder\n";
  @pangene_files = grep {/$cluster_regex$/} readdir(CLUSTDIR);
  closedir(CLUSTDIR);

  if(scalar(@pangene_files) == 0) {
    die "# ERROR: cannot find any valid clusters in $cluster_folder\n";
  }

} else { 

  # 1.2) alternatively parse 0-based BED-file to assign pangenome/graph global coords to individual clusters

  # first it assigns start coords to non-ref pangenes, then end coords are added (two passes), rules: 
  # out-of-order ref genes (non-consecutive in same pangene) are excluded,
  # leading non-ref clusters get assigned global coords before 1st gene model,
  # middle non-ref clusters get assigned global interval defined between prev,next ref gene model on same strand,
  # trailing non-ref clusters get assigned global coords after last gene model,
  # 
  # contents of BED file look like this:
  # #1     NA    NA   gene:BGIOSGA002568      1       0       NA      ...
  # 1   4848  11824   gene:ONIVA01G00010      1       +       gene:ONIVA01G00010
  # #1     NA    NA ..
  # 1 104921  115645  gene:ONIVA01G00100      3       +       gene:ONIVA01G00100
  # #1     NA     NA  gene:BGIOSGA002570      2       0       NA      ...
  # ...
  # Note: a pangene can appear in consecutive lines if it includes neighbor ref genes:
  #1 2531186 2535515 gene:ONIVA01G03460	3 - gene:ONIVA01G03460 ..
  #1 2536276 2536905 gene:ONIVA01G03460	3 - gene:ONIVA01G03470 ..
  #1 2536922 2540925 gene:ONIVA01G03460	3 - gene:ONIVA01G03480 ..

  my ($gchr,$gstart,$gend,$gstrand,$gene);
  my ($gchr2,$gstart2,$gend2,$gstrand2,$gene2,$cluster_id2);
  my ($e,$e2,$elem,$elem2,$dist,$index,$p,$refp);
  my ($ref_genome, $fai_file) = ('','');
  my (%prev_end, %prev_start, %prev_strand);
  my (%nonconsec, %isref, %ref_chr_size);
  my (@BED,@nonref);

  ## find out last position of reference chromosomes:
  # i) find out reference genome
  open(TAB,"<", $INP_dir.'/pangene_matrix.tr.tab') ||
    die "# ERROR: cannot read $INP_dir/pangene_matrix.tr.tab\n";
  while(<TAB>) {
    $ref_genome = (split(/\t/))[1];
    last; 
  }
  close (TAB);
  
  # ii) get ref chr sizes
  $fai_file = dirname(abs_path($INP_dir))."/_$ref_genome.fna.fai";
  open(FAI,"<", $fai_file) ||
    die "# ERROR: cannot read $fai_file\n";
  while(<FAI>) {
    if(/^(\S+)\t(\d+)\t/) {	  
      $ref_chr_size{ $1 } = $2;
    }
  }
  close (FAI); 

  ## actually parse BED

  # i) slurp file into array, columns split; will be read more than once
  open(BED,"<",$pangene_bed_file) ||
    die "# ERROR: cannot read $pangene_bed_file\n";
  while(<BED>) { # 1 4848 11824 gene:ONIVA01G00010 1 + gene:ONIVA01G00010 ..	  
    my @data = split(/\t/,$_);   	  
    push(@BED, \@data);	  
  } 
  close(BED);

  # ii) identify non consecutive reference genes in same pangene,
  # should not be considered for coordinate calculations,
  for $e (0 .. $#BED) {
    next if($BED[$e]->[0] =~ m/^#/);

    $elem = $BED[$e];
    ($gchr,$gstart,$gend,$cluster_id,$gene) = 
      ($elem->[0],$elem->[1],$elem->[2],$elem->[3],$elem->[6]);

    $dist = 0;
    for $e2 ($e+1 .. $#BED) {
      next if($BED[$e2]->[0] =~ m/^#/);
      $dist++;
      last if($dist > $neigh); 

      $elem2 = $BED[$e2];
      ($gchr2,$gstart2,$gend2,$cluster_id2,$gene2) = 
        ($elem2->[0],$elem2->[1],$elem2->[2],$elem2->[3],$elem2->[6]);

      next if($gchr ne $gchr2 || $cluster_id eq $cluster_id2);

      if($gstart2 < $gstart || $gstart2 < $gend) {
        #print "$gchr,$gstart,$gend,$cluster_id,$gene -> $gchr2,$gstart2,$gend2,$cluster_id2,$gene2\n";
        $nonconsec{ $gene } = 1;
      }
    }  
  }  
  
  # iii) 1st pass, assumes pangenes are BED-sorted, 
  # skips non-consecutive, out-of-order ref genes
  for $e (0 .. $#BED) {

    $elem = $BED[$e];
    ($gchr,$gstart,$gend,$cluster_id,$gstrand,$gene) =
      ($elem->[0],$elem->[1],$elem->[2],$elem->[3],$elem->[5],$elem->[6]); 
    #print "$gchr,$gstart,$gend,$cluster_id,$gstrand,$gene\n";

    $index = scalar(@pangene_files);

    if($gchr !~ /^#/) { # ref pangene

      # skip ref genes contained in other ref genes
      if(defined($prev_strand{$gchr}) && $gstrand eq $prev_strand{$gchr} && 
        defined($prev_end{$gchr}) && $gend < $prev_end{$gchr}) {
        next;
      }

      # update pangene: contains several ref gene models
      if(defined($global_coords{ $cluster_id . $cluster_regex })) {

        if(defined($nonconsec{$gene})) { # out-of-order ref gene, skip
          next; 

        } else {	
          $global_coords{ $cluster_id . $cluster_regex }->[2] = $gend;
        }

      } else { # new ref pangene, correct start as BED is 0-based  
        $global_coords{ $cluster_id . $cluster_regex } = 
          [$gchr,$gstart+1,$gend,$gstrand,$index];
        $isref{ $cluster_id . $cluster_regex} = 1;
        push(@pangene_files, $cluster_id . $cluster_regex);  	  
      }	

      # keep track of previous pangene
      $prev_end{$gchr} = $gend;
      $prev_strand{$gchr} = $gstrand;

    } else { #non-ref pangene

      $gchr =~ s/^#//; # remove # from chr name
    
      next if($gchr eq 'unplaced');  

      push(@nonref,$cluster_id . $cluster_regex);

      if(defined($prev_end{$gchr})) { 				

        # middle/trailing non-ref pangene, end coord to be added on 2nd pass
        $global_coords{ $cluster_id . $cluster_regex } = 
          [$gchr,$prev_end{$gchr}+1,0,$prev_strand{$gchr},$index];

      } else {
        # leading non-ref, save only chr
        $global_coords{ $cluster_id . $cluster_regex } = [$gchr,0,0,'',$index]; 
      }

      push(@pangene_files, $cluster_id . $cluster_regex);
    }
  }
  undef(@BED);

  # iv) 2nd pass, complete coords of non-ref pangenes
  foreach $cluster_id (@nonref) {

    ($gchr,$gstart,$gend,$gstrand,$index) = 
      ( $global_coords{$cluster_id}->[0], $global_coords{$cluster_id}->[1],
        $global_coords{$cluster_id}->[2], $global_coords{$cluster_id}->[3],
        $global_coords{$cluster_id}->[4]+1 ); # +1 to get to next pangene
    #print "$gchr,$gstart,$gend,$cluster_id,$gstrand,$index\n";

    # find next ref pangene ($files[$p]) in BED file on same strand & chr
    $refp = -1; 
    foreach $p ($index .. $#pangene_files) {
      next if( !$isref{$pangene_files[$p]} || # nonref
        ($gchr ne $global_coords{$pangene_files[$p]}->[0]) || # wrong chr 	      
        ($gstrand ne '' && $gstrand ne $global_coords{$pangene_files[$p]}->[3]) ); # wrong strand

      $refp = $p; #print ">$cluster_id $p $pangene_files[$p]\n"; 
      last;  
    }	    

    if($gstart == 0) { # leading pangenes

      $global_coords{ $cluster_id }->[1] = 1;
      $global_coords{ $cluster_id }->[2] = $global_coords{$pangene_files[$refp]}->[1] - 1; 
      #print "$gchr,$global_coords{ $cluster_id }->[1],$global_coords{ $cluster_id }->[2],$cluster_id\n";

    } else {
      if($refp == -1) { # trailing non-ref
        $global_coords{ $cluster_id }->[2] = $ref_chr_size{ $gchr };
	#print "$gchr,$gstart,$global_coords{ $cluster_id }->[2],$cluster_id\n";

      } else {	    
        $global_coords{ $cluster_id }->[2] = $global_coords{$pangene_files[$refp]}->[1] - 1;	    
	#print "$gchr,$gstart,$global_coords{ $cluster_id }->[2],$cluster_id\n";
      }
    }   

    # QC  
    if(defined($global_coords{ $cluster_id }) && 
      $global_coords{ $cluster_id }->[0] ne 'unplaced' && 
      $gstart > $global_coords{ $cluster_id }->[2]+1) { # +1 to resist to neighbor overlapping genes ie HORVU.MOREX.r3.3HG0311710
      print "# WARN: not valid non-ref interval for pangene $cluster_id ($gchr,$gstart,$global_coords{ $cluster_id }->[2])\n";
    }
  }

}#1.2

# check whether previous index should be re-used
$index_file = "$gmapdb.ref153positions";

if($INP_ow ||
  !-d "$gmapdb_path/$gmapdb" ||
  !-e "$gmapdb_path/$gmapdb/$index_file") {

  my ($fh, $filename) = tempfile( '/tmp/tempfnaXXXXX', UNLINK => 1);

  foreach $cluster_file (@pangene_files) {

    my ( $ref_geneid, $ref_fasta, $ref_isof_coords, $ref_taxon ) = 
      parse_sequence_FASTA_file( "$cluster_folder/$cluster_file" , 0);

    foreach $gene_id (@$ref_geneid) {  
      foreach $seq (split(/\n/,$ref_fasta->{$gene_id})) {

        if($seq =~ /^>(\S+) (\S+) (\S+) (\S+)/) {

          # FASTA headers in clusters produced by get_pangenes.pl look like this:
          #>transcript:Os01t0100100-01 gene:Os01g0100100 1:2983-10815(+) [Oryza_sativa.IRGSP-1.0.chr1]

          # prepend cluster filename to header and print to temp file
          ($isof_id, $gene_id, $coords, $taxon) = ($1, $2, $3, $4);
          print $fh ">$cluster_file|$isof_id|$gene_id|$coords|$taxon|";

          # add <global coordinates> to the end of header
          if($INP_make_FASTA_reference && $global_coords{$cluster_file}->[0] ne 'unplaced') {
            printf($fh "<%s:%d-%d>\n",
              $global_coords{$cluster_file}->[0],
              $global_coords{$cluster_file}->[1],
              $global_coords{$cluster_file}->[2]);

	  } else {	    
            print $fh "<>\n";
	  }	  

        } else {
          print $fh "$seq\n";
        }
      }
    }
  } 

  # actually create GMAP index from temp file
  $cmd = "$GMAPBUILDBIN -e 0 -d $gmapdb -D $gmapdb_path $filename 2>&1";

  print "# building GMAP index $gmapdb_path/$gmapdb\n\n";
  system($cmd);
  if ( $? != 0 ) {
    die "# ERROR: failed running $GMAPBUILDBIN ($cmd)\n";
  } else {
    print "\n# created GMAP index $gmapdb_path/$gmapdb\n\n";

    if($INP_make_FASTA_reference) {
      move($filename, $pangene_ref_file);
      print "# pangene reference file: $pangene_ref_file\n\n";
    } 
  }

} else {
  print "\n# re-using GMAP index $gmapdb_path/$gmapdb\n\n";

  if(-s $pangene_ref_file) {
    print "# pangene reference file: $pangene_ref_file\n\n";
  }
}

if($INP_seqfile eq '') {
  exit(0);
}


## 2) map input sequence(s) to indexed cluster sequences with GMAP
$cmd = "$GMAPBIN --no-chimeras -t $INP_threads -n 100 -S -d $gmapdb -D $gmapdb_path $INP_seqfile 2>&1";
open(GMAP, "$cmd |") ||
  die "# ERROR: cannot run $cmd\n"; 
while(<GMAP>) { 

  print if($INP_verbose);
  chomp;	

  #>transcript:Os01t0147300-01 gene:Os01g0147300 1:2568051-2570029(+) [Oryza_sativa.IRGSP-1.0.chr1]
  if(/^>(.*)/){ 
    # input FASTA headers could be in any format 
    $seq_id = $1;
    push(@order_geneid, $seq_id);

  } elsif(/^  Path \d+: .*?(\d+) bp\)$/){
    #Path 1: query 1..1768 (1768 bp) => genome gene:...:1..1,820 (1820 bp)
    ($tlen) = ($1); 

  } elsif(/^    Coverage: (\S+) \(query length: (\d+) bp/){	  
    ($cover,$qlen) = ($1, $2);

  } elsif(/^    Percent identity: (\S+)/) {
    #Percent identity: 87.4 (546 matches, 62 mismatches, 17 indels, 0 unknowns)
    ($identity) = ($1);

  } elsif(/^    \+(\S+)\s+\(\d+-\d+\)/){
    # expected to be on the plus/+ strand as these are sequences cut from GFF in the right strand:
    # +cluster.cdna.fna|Os01t0829900-01|gene:Os01g0829900|1:35501996-35505052(-)|[taxon]|<1:35501996-35505052>:1-1820  (1-1768)   97%
    # +cluster.cdna.fna|Os01t0829900-01|gene:Os01g0829900|1:35501996-35505052(-)|[taxon]|:1-1820  (1-1768)   97%

    $cluster_seq_id = $1;

    $regex = '^([^\|]+)\|([^\|]+)\|([^\|]+)\|([^\|]+)\|\[([^\]]+)\]\|<([^>]*)';

    if($identity >= $INP_minident && $cluster_seq_id =~ m/$regex/){

      ($cluster_id, $isof_id, $gene_id, $coords, $taxon, $gcoords) = ($1, $2, $3, $4, $5, $6);
    
      # compile match stats
      $matches{$seq_id}{$cluster_id}{'total'}++;

      # take stats from hit with best cover
      if(!defined($matches{$seq_id}{$cluster_id}{'cover'}) || 
        $cover > $matches{$seq_id}{$cluster_id}{'cover'}) {

        $matches{$seq_id}{$cluster_id}{'cover'} = $cover;
        $matches{$seq_id}{$cluster_id}{'identity'} = $identity;
        $matches{$seq_id}{$cluster_id}{'taxon'} = $taxon;	
        $matches{$seq_id}{$cluster_id}{'coords'} = $coords;
        $matches{$seq_id}{$cluster_id}{'qlength'} = $qlen;
        $matches{$seq_id}{$cluster_id}{'tlength'} = $tlen;
        $matches{$seq_id}{$cluster_id}{'gcoords'} = $gcoords || 'NA';
      }	
    }
  }
}
close(GMAP);


## 3) print results table

open(TSV,">",$INP_outfile) ||
  die "# ERROR: cannot create $INP_outfile\n";

print TSV "#query\tqlength\tpangene\tlength\t".
  "matches\tperc_qcover\tperc_identity\tcoords\ttaxon\tpangenome_coords\n";

foreach $seq_id (@order_geneid) {

  # no matches
  if(!defined($matches{$seq_id})) {
      print TSV "$seq_id\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
      next;
  }
	
  # note there might be 1+ rows per input sequence, 
  # a sequence can potentially match several clusters
  foreach $cluster_id (keys(%{ $matches{$seq_id} })) {

    printf(TSV "%s\t%d\t%s\t%d\t%d\t%1.1f\t%1.1f\t%s\t%s\t%s\n",
      $seq_id,
      $matches{$seq_id}{$cluster_id}{'qlength'},
      $cluster_id,
      $matches{$seq_id}{$cluster_id}{'tlength'},
      $matches{$seq_id}{$cluster_id}{'total'},      
      $matches{$seq_id}{$cluster_id}{'cover'},
      $matches{$seq_id}{$cluster_id}{'identity'},
      $matches{$seq_id}{$cluster_id}{'coords'},
      $matches{$seq_id}{$cluster_id}{'taxon'},
      $matches{$seq_id}{$cluster_id}{'gcoords'}      
    );
  }
}

close(TSV);

print "# results in TSV format: $INP_outfile\n\n";
