use strict;
use warnings;

# This script deletes duplicated probes (same probe_seq_id) in arrays marked as probeset arrays
# and failed datacheck UniqueProbe for this reason
#
# See https://m.ensembl.org/info/docs/api/funcgen
# See https://m.ensembl.org/info/docs/api/funcgen/probe.png

# get probe_seq_id of duplicated probes
my $SQL1 =  "SELECT DISTINCT probe.probe_set_id,probe.probe_seq_id FROM array INNER JOIN array_chip USING (array_id) INNER JOIN probe USING (array_chip_id) WHERE array.is_probeset_array = true GROUP BY array.name, probe.name, probe.probe_set_id HAVING COUNT(distinct probe.probe_id) > 1";

# get probe_id of duplicated probes
my $SQL1_1 =  "SELECT probe_id FROM probe WHERE probe_seq_id=QQ AND probe_set_id=SS";

# delete probe_ids from several tables
my $SQL2_1 =  "DELETE FROM probe_transcript WHERE probe_id=";
my $SQL2_2 =  "DELETE pf, pft FROM probe_feature pf INNER JOIN probe_feature_transcript pft ON pf.probe_feature_id=pft.probe_feature_id WHERE pf.probe_id=";
my $SQL2_3 =  "DELETE FROM probe WHERE probe_id=";
my $SQL2_4 =  "DELETE t1 FROM probe_feature t1 LEFT JOIN probe t2 ON t1.probe_id = t2.probe_id WHERE t1.probe_id IS NOT NULL AND t2.probe_id IS NULL";

if(!$ARGV[1]){ die "# $0 <dbserver-w> <funcgen schema>\n" }
my ($dbserver, $fschema) = @ARGV;

my ($sql,$row,$probe_set_id,$probe_seq_id,$probe_id);

my @dup_probe_seqs = get_n_cols_sql_query($dbserver,$fschema,$SQL1,2);

foreach $row (@dup_probe_seqs){

	($probe_set_id,$probe_seq_id) = split(/,/,$row);
	
	$sql = $SQL1_1;
	$sql =~ s/=QQ/=$probe_seq_id/;
	$sql =~ s/=SS/=$probe_set_id/;

	my @probe_ids = get_n_cols_sql_query($dbserver,$fschema,$sql);
	
	my $reps = 0;
	foreach $probe_id (@probe_ids){

		# skip first probe_id, delete the rest
		if($reps > 0){
			print "# $probe_set_id $probe_seq_id $probe_id\n";

			# delete repeated probe across tables
			get_n_cols_sql_query($dbserver,$fschema,$SQL2_1.$probe_id);
			get_n_cols_sql_query($dbserver,$fschema,$SQL2_2.$probe_id);
			get_n_cols_sql_query($dbserver,$fschema,$SQL2_3.$probe_id);
		}
		$reps++;
	}
}

# finally clean orphan probe_ids
get_n_cols_sql_query($dbserver,$fschema,$SQL2_4);

# returns N column of results from SQL query
# If N > 1 values are CSV
sub get_n_cols_sql_query {
	my ($server,$schema,$SQL,$ncols) = @_;
	my (@results);

	if(!defined($ncols)){ $ncols = 1 }

	open(SQL,"$server $schema -NB -e \"$SQL\" |") ||
		die "# ERROR: cannot run $server $schema -NB -e \"$SQL\"\n";
	while(my $line = <SQL>){
		chomp($line);
		my @data = split(/\t/,$line);
		push(@results,join(',',@data[0 .. $ncols-1]));
	}
	close(SQL);

	return @results;
}
