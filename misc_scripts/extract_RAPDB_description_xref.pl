use strict;
use warnings;
use Encode qw(from_to);

# reads a RAP-DB GFF file and writes a 4-column TSV file with stable_id , gene name, external db and description
# for hierarchical external dbs i) RAP-DB and ii) Oryzabase
# Bruno Contreras Moreira May2019

my $inGFFfile = $ARGV[0] || die "# usage: $0 <GFF file from RAP-DB>\n";

# edit if needed
my %external_db = ('Oryzabase'=>50851, 'RAP-DB'=>50852);

my ($id,$syn,$desc,$extdb);
my (%seen);

my $gene_file = $inGFFfile.".gene.tsv";
my $mRNA_file = $inGFFfile.".transcript.tsv";

open(GENE,">",$gene_file) || die "# ERROR: cannot write to $gene_file\n";
open(MRNA,">",$mRNA_file) || die "# ERROR: cannot write to $mRNA_file\n";

open(GFF,$inGFFfile) || die "# ERROR: cannot read $inGFFfile\n";
while(<GFF>){
	chomp;
	my @col = split(/\t/);

	#chr01	irgsp1	gene	2983	10815	.	+	.	ID=Os01g0100100;Name=Os01g0100100;Note=RabGAP/TBC domain containing protein. (Os01t0100100-01);
	#chr01	irgsp1	gene	35623	41136	.	+	.	ID=Os01g0100900;Name=Os01g0100900;Note=Sphingosine-1-phosphate lyase%2C Disease resistance response (Os01t0100900-01);...;RAP-DB Gene Symbol Synonym(s)=SPL1%2C OsSPL1;RAP-DB Gene Name Synonym(s)=SPHINGOSINE-1-PHOSPHATE LYASE 1%2C Sphingosine-1-Phoshpate Lyase 1;
	#chr01	irgsp1	mRNA	27143	28644	.	+	.	ID=Os01t0100700-01;Name=Os01t0100700-01;Parent=Os01g0100700;Oryzabase Gene Symbol Synonym(s)=RPS5;Oryzabase Gene Name Synonym(s)=ribosomal protein S5%2C ribosomal protein small subunit 5;Note=Similar to 40S ribosomal protein S5-1.
	#chr01	irgsp1	mRNA	35623	41136	.	+	.	ID=Os01t0100900-01;Name=Os01t0100900-01;Parent=Os01g0100900;RAP-DB Gene Symbol Synonym(s)=SPL1%2C OsSPL1;RAP-DB Gene Name Synonym(s)=SPHINGOSINE-1-PHOSPHATE LYASE 1%2C Sphingosine-1-Phoshpate Lyase 1;Oryzabase Gene Symbol Synonym(s)=OsSPL%2C OsSPL1;Oryzabase Gene Name Synonym(s)=Sphingosine-1-phosphate lyase%2C sphingosine-1-phosphate lyase 1;Note=Sphingosine-1-phosphate lyase%2C Disease resistance response;

	if($col[2] eq 'gene' || $col[2] eq 'mRNA'){

		if($col[8] =~ /^ID=([^;]+)/){ $id=$1 }

		next if(defined($seen{$id}{$col[2]}));

		if($col[8] =~ /Note=([^;]+)/){ 
			$desc=URI2string($1); 
			$desc = (split(/. \(O/,$desc))[0];
	
		}

		if($col[8] =~ /RAP-DB Gene Symbol Synonym\(s\)=([^;\s]+)/){ 
			$syn=URI2string($1);
			$syn=clean_name($syn); 
			$extdb = $external_db{'RAP-DB'};
		}
		elsif($col[8] =~ /Oryzabase Gene Symbol Synonym\(s\)=([^;\s]+)/){ 
			$syn=URI2string($1);
			$syn=clean_name($syn);
			$extdb = $external_db{'Oryzabase'};
		}
		else{ $syn=$extdb="NULL" }

		$syn =~ s/,.*$//g;

		if($col[2] eq 'gene'){
			print GENE "$id\t$syn\t$extdb\t$desc\n";        
		}
		else{
			print MRNA "$id\t$syn\t$extdb\t$desc\n";
		}
        }

	# take only first 
	$seen{$id}{$col[2]} = 1;

	
}
close(GFF);
 
# taken from https://stackoverflow.com/questions/4510550/using-perl-how-do-i-decode-or-create-those-encodings-on-the-web
sub URI2string {

	my ($string) = @_;

	$string =~ s/%([A-Fa-f\d]{2})/chr hex $1/eg;

	# any HTML special chars as well
	$string =~ s/&apos/'/g; # this actually XML
	$string =~ s/&(\S+)?;/$1/eg;

	return $string;
}

# after reading https://github.com/Ensembl/ensj-healthcheck/blob/bf400cd6f87f542856bd73a5e38745b441bbbd07/src/org/ensembl/healthcheck/testcase/EnsTestCase.java
sub clean_name {
	my ($string) = @_;

	$string =~ s/^[\[\\:\\;\n\r\t\~\]]//g;
        $string =~ s/;.*$//; 
	$string =~ s/[\[\\:\\;\n\r\t\~\]]$//g;
        
        return $string;
}
