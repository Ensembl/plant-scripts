use strict;
use warnings;

# This script deletes duplicated probes in arrays which are not marked as probeset arrays
# and failed datacheck UniqueProbe for this reason
#
# See https://m.ensembl.org/info/docs/api/funcgen
# See https://m.ensembl.org/info/docs/api/funcgen/probe.png

# get probe_seq_id of duplicated probes
my $SQL1 =  "SELECT probe.probe_seq_id FROM array INNER JOIN array_chip USING (array_id) INNER JOIN probe USING (array_chip_id) WHERE array.is_probeset_array = false GROUP BY array.name, probe.name   having COUNT(distinct probe.probe_id) > 1";

# get probe_id of duplicated probes
my $SQL1_1 =  "SELECT probe_id FROM probe WHERE probe_seq_id=";

# delete probe_ids from several tables
my $SQL2_1 =  "DELETE FROM probe_transcript WHERE probe_id=";
my $SQL2_2 =  "DELETE pf, pft FROM probe_feature pf INNER JOIN probe_feature_transcript pft ON pf.probe_feature_id=pft.probe_feature_id WHERE pf.probe_id=";
my $SQL2_3 =  "DELETE FROM probe WHERE probe_id=";
my $SQL2_4 =  "DELETE t1 FROM probe_feature t1 LEFT JOIN probe t2 ON t1.probe_id = t2.probe_id WHERE t1.probe_id IS NOT NULL AND t2.probe_id IS NULL";

if(!$ARGV[1]){ die "# $0 <dbserver-w> <funcgen schema>\n" }
my ($dbserver, $fschema) = @ARGV;

my @dup_probe_seqs = get_first_col_sql_query($dbserver,$fschema,$SQL1);

foreach my $probe_seq_id (@dup_probe_seqs){

	my @probe_ids = get_first_col_sql_query($dbserver,$fschema,$SQL1_1.$probe_seq_id);
	
	my $reps = 0;
	foreach my $probe_id (@probe_ids){

		# skip first probe_id, delete the rest
		if($reps > 0){
			print "# $probe_seq_id $probe_id\n";

			# delete repeated probe across tables
			get_first_col_sql_query($dbserver,$fschema,$SQL2_1.$probe_id);
			get_first_col_sql_query($dbserver,$fschema,$SQL2_2.$probe_id);
			get_first_col_sql_query($dbserver,$fschema,$SQL2_3.$probe_id);
		}
		$reps++;
	}
}

# finally clean orphan probe_ids
get_first_col_sql_query($dbserver,$fschema,$SQL2_4);

# returns 1st column of results from SQL query
sub get_first_col_sql_query {
	my ($server,$schema,$SQL) = @_;
	my (@results);

	open(SQL,"$server $schema -NB -e \"$SQL\" |") ||
	    die "# ERROR: cannot run $server $schema -NB -e \"$SQL\"\n";
	while(my $line = <SQL>){
	    chomp($line);
	    my @data = split(/\t/,$line);
		push(@results,$data[0]); 
	}
	close(SQL);

	return @results;
}
