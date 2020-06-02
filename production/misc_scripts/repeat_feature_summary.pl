#!/usr/bin/env/perl
# Contributed by James Allen

use strict;
use warnings;
use feature 'say';

use Getopt::Long qw(:config no_ignore_case);
use Path::Tiny;

use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->no_version_check(1);

my ($reg_file, $species, @logic_name, $all, $out_dir, $bedtools_dir);

GetOptions(
  "reg_file=s",     \$reg_file,
  "species=s",      \$species,
  "logic_name=s",   \@logic_name,
  "all",            \$all,
  "out_dir=s",      \$out_dir,
  "bedtools_dir=s", \$bedtools_dir,
);

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_all($reg_file);
my $dba = $registry->get_DBAdaptor($species, 'core');
my $dbh = $dba->dbc->db_handle;
my $genome = $dba->get_GenomeContainerAdaptor;

@logic_name = () unless @logic_name;

$out_dir = '.' unless $out_dir;
path($out_dir)->mkpath;

if ($bedtools_dir) {
  $ENV{PATH} = "$ENV{PATH}:$bedtools_dir";
}

my $logic_name_sql = '
  SELECT logic_name
  FROM analysis INNER JOIN repeat_feature USING (analysis_id)
  GROUP BY logic_name
  ORDER BY logic_name
';

my $coord_sql = '
  SELECT sr.name, rf.seq_region_start - 1, rf.seq_region_end
  FROM seq_region sr
  INNER JOIN repeat_feature rf USING (seq_region_id)
  INNER JOIN analysis USING (analysis_id)
  WHERE logic_name = ?
  ORDER BY sr.name, rf.seq_region_start;
  ';

my $coord_all_sql = '
  SELECT sr.name, rf.seq_region_start - 1, rf.seq_region_end
  FROM seq_region sr
  INNER JOIN repeat_feature rf USING (seq_region_id)
  ORDER BY sr.name, rf.seq_region_start;
  ';

main();

sub main {
  say join("\t", ($species, 'genome_length', $genome->get_ref_length));

  if (!scalar(@logic_name)) {
    @logic_name = @{$dbh->selectcol_arrayref($logic_name_sql)};
  }

  foreach my $logic_name (@logic_name) {
    my $sth = $dbh->prepare($coord_sql);
    $sth->execute($logic_name);
    
    create_bed_file($sth, $logic_name);
  }

  if ($all) {
    my $sth = $dbh->prepare($coord_all_sql);
    $sth->execute();
    
    create_bed_file($sth, 'all');
  }
}

sub create_bed_file {
  my ($sth, $category) = @_;
  
  my $filename = "$out_dir/$species.$category.bed";
  my $nr_file  = path("$filename");
  my $tmp_file = path("$filename.tmp");
  
  while (my $batch = $sth->fetchall_arrayref(undef, 100000)) {
    $tmp_file->append(map { join("\t", @$_) . "\n" } @$batch);
  }
  
  my $nr_cmd = "bedtools merge -i $filename.tmp > $filename";
  if (system($nr_cmd) == 0) {
    $tmp_file->remove;
  } else {
    die "Failed to run '$nr_cmd': $?";
  }
  
  my $total = 0;
  my @lines = $nr_file->lines({ chomp => 1 });
  foreach my $line (@lines) {
    my @line = split("\t", $line);
    $total += $line[2] - $line[1];
  }
  
  say join("\t", ($species, $category, $total));
}
