# Modified from phyTools at https://github.com/eead-csic-compbio/get_homologues

package FASTATools;
require Exporter;

@ISA = qw(Exporter);
@EXPORT = qw(
  read_FASTA_file_array
  NAME SEQ 
  );

use strict;

#################################################################

use constant NAME => 0;
use constant SEQ  => 1;
use constant IDENTICALS => 2;
use constant NUMBER_IDENTICALS => 3;

sub read_FASTA_file_array
{
  # in FASTA format
  # returns a reference to a 2D array for 4 secondary keys: NAME,SEQ,IDENTICALS,NUMBER_IDENTICALS
  # first valid index (first sequence) is '0'
   
  my ( $infile, $skipamino, $skipidentical ) = @_;
  my (@FASTA,$name,$seq,$n_of_sequences,$magic);
  $n_of_sequences = -1;

  # check input file format and open it accordingly
  open(INFILE,$infile) || die "# read_FASTA_sequence_array: cannot read $infile, exit\n";
  sysread(INFILE,$magic,2);
  close(INFILE);

  if($infile =~ /\.gz$/ || $magic eq "\x1f\x8b") # GZIP compressed input
  {
    if(!open(FASTA,"gzip -dc $infile |"))
    {
      die "# read_FASTA_sequence_array: cannot read GZIP compressed $infile $!\n"
        ."# please check gzip is installed\n";
    }
  }
  elsif($infile =~ /\.bz2$/ || $magic eq "BZ") # BZIP2 compressed input
  {
    if(!open(FASTA,"bzip2 -dc $infile |"))
    {
      die "# read_FASTA_sequence_array: cannot read BZIP2 compressed $infile $!\n"
        ."# please check bzip2 is installed\n";
    }
  }
  else{ open(FASTA,"<$infile") || die "# read_FASTA_sequence_array: cannot read $infile $!\n"; }

  while(<FASTA>)
  { 
    next if(/^$/ || /^#/);
    if(/^\>(\S+)/) # first non blank
    {
      $n_of_sequences++; # first sequence ID is 0
      $name = $1;  
      $FASTA[$n_of_sequences][NAME] = $name;
    }
    elsif($n_of_sequences>-1)
    {
      $_ =~ s/[\s|\n|\-|\.]//g;
      $FASTA[$n_of_sequences][SEQ] .= $_; # conserve case
    }
  }
  close(FASTA);

  if($n_of_sequences == -1){ return \@FASTA }
  if($skipamino || $skipidentical)
  {
    my ($n_of_nr_sequences,@nrFASTA,%identical) = 0;
    if($skipidentical)
    {
      foreach $seq ( 0 .. $n_of_sequences-1 )
      {
        $FASTA[$seq][IDENTICALS] = '';
        next if($identical{$seq});
        foreach my $seq2 ( $seq + 1 .. $n_of_sequences )
        {
          if($FASTA[$seq][SEQ] eq $FASTA[$seq2][SEQ] ||  # same length
            $FASTA[$seq][SEQ] =~ /$FASTA[$seq2][SEQ]/ || # different length
            $FASTA[$seq2][SEQ] =~ /$FASTA[$seq][SEQ]/)
          {
            $identical{$seq2} = $seq;
            $FASTA[$seq][IDENTICALS] .= "$seq2,"; #print "mira $seq $FASTA[$seq][IDENTICALS]\n";
            $FASTA[$seq][NUMBER_IDENTICALS]++;
          }
        }
      }
    }
    foreach $seq ( 0 .. $n_of_sequences )
    {
      # skip protein sequences
      if($FASTA[$seq][SEQ] =~ /[QEYIDFHKLVM]/i)
      {
        print "# read_FASTA_sequence_array : skipped amino acid sequence ($FASTA[$seq][NAME])\n";
        next;
      }
      if($identical{$seq})
      {
        print "# read_FASTA_sequence_array : skipped identical sequence: ".
          "$FASTA[$seq][NAME]\n";
        next;
      }
      $nrFASTA[$n_of_nr_sequences] = $FASTA[$seq];

      # keep track of identical sequences
      if($FASTA[$seq][NUMBER_IDENTICALS])
      {
        chomp($nrFASTA[$n_of_nr_sequences][NAME]);
        $nrFASTA[$n_of_nr_sequences][NAME] .=
          " | identical sequences=$nrFASTA[$seq][NUMBER_IDENTICALS]";
      }
      $n_of_nr_sequences++;
    }
    return \@nrFASTA;
  }
  else{ return \@FASTA }
}


1;
