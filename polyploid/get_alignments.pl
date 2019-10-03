#to run: perl get_alignments.pl -release 42 > wheat9.txt
#!/usr/bin/perl

# program to calculate IHVs for wheat
# $opt_test is used flexibly when testing, it doesn't have a role in normal production
use 5.14.0;
use vars qw($opt_test $opt_release);

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Data::Dumper;
use Getopt::Long;
use DBI;
use strict;
&GetOptions("test:s", "release=i");

my $registry = 'Bio::EnsEMBL::Registry';
die "Release not set\n" if ! $opt_release;
my $egRelease = $opt_release;
my $eRelease = $egRelease + 53;
my $assembly;
my $dbId;

if ($egRelease < 25) {
	
	$assembly = 1;
	$dbId = 2000;
	
} elsif ($egRelease < 32) {
	
	$assembly = 2;
	$dbId = 2054;

} else {
	
	$assembly = 3;
	$dbId = 2054;
}

##Genome DB Id for Zavitan principal component
$dbId     = 2116;
$assembly = 1;

##mlss for LASTZ between A and B components for Zavitan
my $mlss_id  = 9741;

#$registry->load_registry_from_multiple_dbs( {
#    -host => 'mysql-eg-mirror.ebi.ac.uk',
#    -port => xxx,
#    -user => 'xxx',
#    -pass => 'xxxx',
#    -db_version => $eRelease
#});

#my $reg_conf = '/homes/gnaamati/registries/prod1-compara-merged.reg';
my $reg_conf = '/homes/gnaamati/registries/prod2-compara-staging2.reg';
$registry->load_all($reg_conf);


##
## Get all homeologues from Compara DB on Staging
##

my $dsn = 'DBI:mysql:host=mysql-eg-staging-1.ebi.ac.uk;port=xxx;' . 
          'database=ensembl_compara_plants_41_94';


my $dsn = 'DBI:mysql:host=mysql-eg-staging-2.ebi.ac.uk;port=xxx;' . 
          'database=ensembl_compara_plants_42_95';

#mysql://ensrw:writ3rp1@mysql-eg-prod-1.ebi.ac.uk:4238/
                    
print STDERR $dsn, "\n";
          
my $user = 'xxx';
my $user = 'xxx';
my $password = 'xxxx';
my $password = 'xxxx';
my $password = '';

my $dbh = DBI->connect($dsn, $user, $password)
  || die "Cannot connect to server: $DBI::errstr\n";

# find 1:1 homoeologues to detect IHVs

my $query1 = "
 select m1.stable_id AS mm1,
        m2.stable_id AS mm2, 
        m1.genome_db_id AS gdb1,
        m2.genome_db_id AS gdb2
 FROM homology h, 
     gene_member m1, 
     homology_member hm1,  
     gene_member m2, 
     homology_member hm2
 WHERE  hm1.homology_id = h.homology_id
 AND hm2.homology_id = hm1.homology_id  
 AND m1.gene_member_id = hm1.gene_member_id   
 AND m2.gene_member_id = hm2.gene_member_id
 AND m1.stable_id != m2.stable_id
 AND h.description = 'homoeolog_one2one'";


my $psmt = $dbh -> prepare($query1);
$psmt -> execute();
my %orthologues;
	
while ((my @results) = $psmt -> fetchrow()) {
	
	my $gene1Id = $results[0];
	my $gene2Id = $results[1];
	my $gene1genomeId = $results[2];
	my $gene2genomeId = $results[3];
	##Just get Zavitan genes for now
        if ($gene1Id =~ /^Traes/) {
            next;
        }

	${$orthologues{$gene1Id}}{$gene2Id} = 1;
        next;
	
        #I don't think we need this
	# only collect wheat genes (it is faster to do this as a post-SQL filter than in
	# SQL selection)
	
	#if ((($egRelease > 32) &&
	#     (($gene1genomeId == 2081) || 
	#      ($gene1genomeId == 2082) || 
	#      ($gene1genomeId == 2083)) && 
	#     (($gene2genomeId == 2081) || 
	#      ($gene2genomeId == 2082) || 
	#      ($gene2genomeId == 2083))) ||
	#    (($gene1genomeId == 2054) && ($gene2genomeId == 2054))) {
	     
		# ignore plastid genes
	
	#	if (($gene1Id !~ /^EPlTAEG/) && ($gene2Id !~ /^EPlTAEG/)) {
        #	
	#		${$orthologues{$gene1Id}}{$gene2Id} = 1;
	#	}
	#}

}


#print Dumper \%orthologues;
#use Storable; 
#store \%orthologues, 'pk_orth.str';
#die('');

#my $hashref = retrieve('pk_orth.str');
#my %orthologues = %$hashref;

#print Dumper $hashref;
#die('');



# test to see if we have a group of 3 genes that are each orthologues of each other
# and if not what the number of AB, AD, and BD pairs are
#


my ($full, $ab, $ad, $bd);

for my $orthologue (keys %orthologues) {

	my $ok = 0;
	
	if (scalar keys %{$orthologues{$orthologue}} == 2) {
	
		$ok = 1;
		my %test;
		my @test = ($orthologue);
		push @test, keys %{$orthologues{$orthologue}};
		
		for my $test (@test) {
		
			$test{$test} = 1;
		}
		
		for my $test (@test) {
		
			$ok = 0 if scalar keys %{$orthologues{$test}} != 2;
			
			for my $key (keys %{$orthologues{$test}}) {
	
				if (! $test{$key}) {
		
					$ok = 0;
				}
			}
		}
		
		if ($ok == 1) {
	
			$full ++;
		}
	} 
	
	if ($ok == 0) {
	
		for my $key (keys %{$orthologues{$orthologue}}) {
		
			    if ($egRelease > 38) {
			
				if (($key =~ /TraesCS\dA/) && ($orthologue =~ /TraesCS\dB/)) {
			    
			    	$ab++;
			
				} elsif (($key =~ /TraesCS\dA/) && ($orthologue =~ /TraesCS\dD/)) {
			    
			    	$ad++;
			
				} elsif (($key =~ /TraesCS\dB/) && ($orthologue =~ /TraesCS\dD/)) {
			    
			    	$bd++;
			    }

                    }

                        
                        
                        elsif ($egRelease > 32) {
			
				if (($key =~ /TRIAE_CS42_\dA/) && ($orthologue =~ /TRIAE_CS42_\dB/)) {
			    
			    	$ab++;
			
				} elsif (($key =~ /TRIAE_CS42_\dA/) && ($orthologue =~ /TRIAE_CS42_\dD/)) {
			    
			    	$ad++;
			
				} elsif (($key =~ /TRIAE_CS42_\dB/) && ($orthologue =~ /TRIAE_CS42_\dD/)) {
			    
			    	$bd++;
			    }
            
            } else {		
		
				if (($orthologue =~ /Traes_\dA/) &&
                    (($key =~ /Traes_\dB/) || ($key =~ /TRAES3B/))) {
                            
                    $ab++;
                        
                } elsif (($orthologue =~ /Traes_\dA/) &&
                         ($key =~ /Traes_\dD/)) {
                            
                    $ad++;
                        
                } elsif ((($orthologue =~ /Traes_\dB/) || ($key =~ /Traes_\dD/)) &&
                          ($key =~ /Traes_\d+D/)) {
                            
                    $bd++;
                }
			}
		}
	} 
}

print STDERR join "\t", "Keys:", scalar keys %orthologues,
                        "3 way: ", $full, 
                        "AB: ", $ab, 
                        "AD: ", $ad, 
                        "BD: ", $bd;
                        
print STDERR "\n";

# retrieve alignment data

my %segment2segment;
my $query2;



#AND mlss.method_link_species_set_id in (9452, 9465, 9466)
my $query2 = get_segment2segment_sql($mlss_id, $dbId, $egRelease);

#if ($opt_test =~ /^\d+$/) {
if (0) {

	$query2 .= 
		#" AND d.name = 'TGACv1_scaffold_593044_7BS' 
		#  AND d2.name='TGACv1_scaffold_621721_7DS'";
		" AND d.name = '1A' 
		  AND d2.name='1B'";
}

$query2 .= ") f";

warn "$query2\n\n";


my $psmt = $dbh -> prepare($query2);
$psmt -> execute();
	
while ((my @results) = $psmt -> fetchrow()) {
    ##record which dnafrags (hereafter called segments) are aligned to each other
    push @{${$segment2segment{$results[1]}}{$results[2]}}, \@results;
}


##===============Storable options for debugging ========================= 
#print STDERR scalar "Number of matching segments: ", scalar keys %segment2segment, "!!\n";
#die if scalar keys %segment2segment == 0;

#use Storable; 
#store \%segment2segment, 'paul_kersey_file.str'; 
#print Dumper \%segment2segment;

#my $hashref  = retrieve('paul_kersey_file.str');
#%segment2segment = %{$hashref};
#print Dumper $hashref;
##========================================== 


# get coordinates for all genes

my $dsn = 'DBI:mysql:host=mysql-eg-staging-1.ebi.ac.uk;port=4160;' . 
          'database=triticum_aestivum_core_33_86_3';

my $user='ensro';
my $password = '';		
		         
my $dsn = 'DBI:mysql:host=mysql-eg-mirror.ebi.ac.uk;port=4157;' . 
          'database=triticum_aestivum_core_' . $egRelease . '_' . $eRelease . '_' . 
          $assembly;
          
my $user = 'ensrw';
my $password = 'writ3r';
 
##Take from prod-2
#my $dsn = 'DBI:mysql:host=mysql-eg-mirror.ebi.ac.uk;port=4157;' . 
#          'database=triticum_aestivum_core_' . $egRelease . '_' . $eRelease . '_' . 
#          $assembly;

my $dsn = 'DBI:mysql:host=mysql-eg-prod-2.ebi.ac.uk;port=4239;' .
          'database=triticum_dicoccoides_core_42_95_1';

          
my $user = 'ensrw';
my $password = 'writ3rp2';
 

my $dbh = DBI->connect($dsn, $user, $password)
  || die "Cannot connect to server: $DBI::errstr\n";
  
my $query3 = "SELECT g.stable_id, s.name, g.seq_region_start, g.seq_region_end, c.name
			  FROM gene g, seq_region s, coord_system c
              WHERE s.seq_region_id = g.seq_region_id
              AND c.coord_system_id = s.coord_system_id";

my $psmt = $dbh -> prepare($query3);
$psmt -> execute();
my %gene2seqRegion;
my %gene2start;
my %gene2end;
my %seqRegion2coordSystem;
	
while ((my @results) = $psmt -> fetchrow()) {
	
	my $geneId = $results[0];
	
	# ignore plastid genes
	
	if ($geneId !~ /^EPlTAEG/) {
	
		my $seqRegionName = $results[1];
		$gene2seqRegion{$geneId} = $seqRegionName;
		$gene2start{$geneId} = $results[2];
		$gene2end{$geneId} = $results[3];
		$seqRegion2coordSystem{$geneId} = $seqRegionName;;
	}
}

#print Dumper \%seqRegion2coordSystem;
#die('');

# now the real work starts, finding alignments supporting orthology calls and looking for
# variants

#my $slice_adaptor = $registry->get_adaptor( 'wheat', 'Core', 'Slice' );
my $slice_adaptor = $registry->get_adaptor( 'triticum_dicoccoides', 'Core', 'Slice' );
my ($seqCount, $seqNoCount);
my $indel = 0;
my $id = 1;
my %idHash;
my %id2key;
my $supportingAlignments; # the number of alignments supporting pairs
my $unsupportedPairs; # the number of pairs without supporting alignments
my %hasOrthologues; 
my %geneKeys;
my %groups;
my %globalMatchLength; # for all alignments
my %matchLength; # for one alignment
my %mismatch;
my %methodCount;

open (SUP, ">support.dat") || die "Could not open support.dat\n";
open (ALIGNMENT, ">alignment.dat") || die "Could not open alignment.dat\n";

# for each gene

print STDERR "Genes:", scalar keys %gene2seqRegion, "\n";

GENE: for my $gene (keys %gene2seqRegion) {

	my $seqRegion = $gene2seqRegion{$gene};
	my $supportedOrthologues = 0;
	
	# for each orthologue of the gene
	
	PAIR: for my $orthologue (keys %{$orthologues{$gene}}) {			
		
		$hasOrthologues{$gene} = 1;
		
		# each pair will be encountered twice, but we only want to process it once
		
		my $geneKey = join '', sort ($gene, $orthologue);
                #say "$geneKey\n";
			
		if ($geneKeys{$geneKey}) {
		
			next PAIR;
		
		} else {
		
			$geneKeys{$geneKey} = 1;
		}
		
		# use form of seqRegion name to identify genome
		
		my $orthologueSeqRegion = $gene2seqRegion{$orthologue};	
		my ($genome1, $genome2);
		


		#  AND d2.name='TGACv1_scaffold_621721_7DS'";
                #  chr1A


		if ($egRelease > 38) {
                        #say "seq region = $seqRegion";
		
			($genome1) = ($seqRegion =~ /\d(\w)/);
			($genome2) = ($orthologueSeqRegion =~ /\d(\w)/);
                        #warn "genome1 = $genome1";
                        #warn "genome2 = $genome2";
                        #die;
			
		}
                elsif ($egRelease > 32) {
		
			($genome1) = ($seqRegion =~ /scaffold_\d+_\d(\w)/);
			($genome2) = ($orthologueSeqRegion =~ /scaffold_\d+_\d(\w)/) 
			
		} else {
		
			($genome1) = ($seqRegion =~ /\d(\w)/);
        	($genome2) = ($orthologueSeqRegion =~ /\d(\w)/);	
		}
		
		#print STDERR join ',', 'GR', 
	        #						$seqRegion, 
		#						$genome1, 
		#						$orthologueSeqRegion, 
		#						$genome2;
		
		#print STDERR "\n";
                #
		
		my $geneStart1 = $gene2start{$gene};
		my $geneEnd1 = $gene2end{$gene};
		my $geneStart2 = $gene2start{$orthologue};
		my $geneEnd2 = $gene2end{$orthologue};

                #say $geneStart1;
                #say $geneEnd2;
                #say $geneEnd1;
                
                
                ###Should I remove this?
                #$seqRegion = 'chr'.$seqRegion;
                #$orthologueSeqRegion = 'chr'.$orthologueSeqRegion;
                
                
                #say "seqRegion is $orthologueSeqRegion";
                #print Dumper \%segment2segment;
                #die();
				
		#  if there are DNA alignments potentially supporting the orthology
	
		if (${$segment2segment{$seqRegion}}{$orthologueSeqRegion}) {
                        

			
			my %alignmentHash;
			my %group2alignments;
			undef %matchLength;
			my ($mismatch1, $mismatch2);
			
			# retrieve the CIGAR strings defining the alignment of each orthologue to the
			# overall alignment
						
			RESULTS: for my $result 
				(@{${$segment2segment{$seqRegion}}{$orthologueSeqRegion}}) {
			
				my @results = @$result;
				my $groupId = $results[0];
				my $alignmentStart1 = $results[3];
				my $alignmentStart2 = $results[4];
				my $alignmentEnd1 = $results[5];
				my $alignmentEnd2 = $results[6];
			
				# only work with alignments overlapping genes
			
				if (($alignmentStart1 <= $geneEnd1) && 
				    ($alignmentEnd1 >= $geneStart1) &&
			    	($alignmentStart2 <= $geneEnd2) && 
			    	($alignmentEnd2 >= $geneStart2)) {
			      	 
			      	# calculate the extent to which the alignment overlaps the genes in
			      	# the orthologue pairs
			      	
			      	$alignmentHash{$groupId} +=
			      		_calculateDistance($geneStart1, 
			      						   $geneEnd1, 
			      						   $geneStart2, 
			      						   $geneEnd2, 
			      						   $alignmentStart1, 
			      						   $alignmentEnd1, 
			      						   $alignmentStart2, 
			      						   $alignmentEnd2);
			      	
			      	push @{$group2alignments{$groupId}}, join ',', @results; 
			      				      
				} 
			}
			
			if (scalar keys %alignmentHash == 0) {
			
				$unsupportedPairs ++;
				next PAIR;
			}
					
			$supportedOrthologues ++;

			# if multiple alignments pick the one with the longest cumulative
			# overlap of both genes
							
			my @sortedResults = sort {$alignmentHash{$b} <=> $alignmentHash{$a}} 
										keys %alignmentHash;
										
			my $groupId = $sortedResults[0];
			
			print ALIGNMENT join "\t", $gene, $orthologue, $groupId;
			print ALIGNMENT "\n";
											
			# only parse each group once
			
			if ($groups{$groupId}) {
			
				next PAIR;
			}
			
			$groups{$groupId} = $gene . $orthologue;
			
			# method to filter out redundant data (may no longer be needed)
			# output is a subset of input
								
			my @alignments = @{_removeRedundantGroups($group2alignments{$groupId})};
			
			# for each block in the best group
			
			my $first = 1;
			my %alignments;
			
			for my $alignment (@alignments) {
				
				# potential problem - same alignment coming twice in the same group
				
				if ($alignments{$alignment}) {
				
					print STDERR "Duplicate $alignment\n";
					print STDERR $groupId, " ", $groups{$groupId}, "\n";
					print STDERR $alignments{$alignment}, "\n";
					print STDERR $groupId, " ", $groups{$alignments{$alignment}}, "\n";
					die;
				
				} else {
				
				 	$alignments{$alignment} = $groupId;
				}
				
				print STDERR "Current alignment: $alignment\n";
				
				# parse alignment data
				
				my @alignmentData = split ',', $alignment;
				my $alignmentStart1 = $alignmentData[3];
				my $alignmentStart2 = $alignmentData[4];
				my $alignmentEnd1 = $alignmentData[5];
				my $alignmentEnd2 = $alignmentData[6];
				my $orientation1 = $alignmentData[7];
				my $orientation2 = $alignmentData[8];
				my $cigar1 = $alignmentData[9];
				my $cigar2 = $alignmentData[10];
				
				# keep track of the number of alignment groups supported by each method
				
				if ($first == 1) {
					
					my $methodId = $alignmentData[11];
					$methodCount{$methodId} ++;
					$first = 0;
				}
				
				$supportingAlignments ++;
				
				# reformat alignments for subsequent analysis
				my @cigars = ($cigar1, $cigar2);
				
				my ($matchLength, $cigarArrays) = 
					_reverseAlignments($genome1, $genome2, \@cigars);
				
				print STDERR "Match length: $matchLength\n";
					
				my @cigarArrays = @$cigarArrays;
			
				# in Ensembl, we define an alignment on opposite directions in each 
				# sequence thusly (1:2)-(1:2)-R
				# but for this program, we want to use an alternative representation i.e.
				# (1:2)-(2:1)
				
				# we also want to define a pointer to the current location, which starts
				# one position "before" the alignment (in the direction of travel) 
				# because we are going to increment at the start of each operation not 
				# after 
			
				my @current;
				my $location1;
				my $location2;
			
				if ($orientation1 == 1) {
			
					$location1 = $alignmentStart1 - 1;
				
				} elsif ($orientation1 == -1) {
			
					$location1 = $alignmentEnd1 + 1;
					$alignmentStart1 = $alignmentEnd1;
				}	
			
				if ($orientation2 == 1) {
			
					$location2 = $alignmentStart2 - 1;
				
				} elsif ($orientation2 == -1) {
			
					$location2 = $alignmentEnd2 + 1;
					$alignmentStart2 = $alignmentEnd2;
				}
				
				print STDERR "Initial locations: $location1 $location2\n";	
				
				# now iterate through the alignment
				# 2 sorts of variants - gaps between "match" runs
				# but also, mismatches exist within match sequences
				# so we have to capture all gaps, but also, examine all matches
			
				my @features1;
				my @features2;
				my @subFeatures1;
				my @subFeatures2;
				my $matchCount  = 0;
				
				for (my $n = 0; $n < $matchLength; $n ++) {
				
					# data integrity check
					# in a pairwise alignment, at least one of the sequences must be 
					# present in every position of the alignment 
					# (or else that position would not exist!)
				
					if ((${$cigarArrays[0]}[$n] ne 'M') && 
					    (${$cigarArrays[1]}[$n] ne 'M')) {
					
						die "Cigar strings misaligned\n";;
					}
					
					# both positions match
					
					if ((${$cigarArrays[0]}[$n] eq 'M') && 
					    (${$cigarArrays[1]}[$n] eq 'M')) {
					
						# increment current location
					
						$location1 = _increment($location1, $orientation1);
						$location2 = _increment($location2, $orientation2);
						
						# was there a mismatch at the previous position?
						
						if ($current[0] ne '') {
						
							# previous (non-matching) features are now completed - now
							# retrieve sequence for them
							# note we have separate features on both seq regions (i.e. we 
							# are coding the alignment twice, once wrt each reference)
						
							_closeFeatures($current[0], 
							               $current[1], 
							               \@features1, 
							               \@features2,
							               $orientation1, 
							               $orientation2);
							               
							$matchCount ++;
							$current[0] = '';
							$current[1] = '';
						}
					
					} else {
					
						if (($current[0] eq '') && ($current[1] eq '')) {
					
							# (new) mismatch, no live mismatch
							# either, end of a match region (which we want to analyse)
							# or, the alignment starts with a mistmatch (not interesting)
							# note, incrementation of location variables hasn't yet 
							# occurred
						
							print STDERR join ",", "Getting indels:",
													$n, 
							                        ${$cigarArrays[0]}[$n], 
							                        ${$cigarArrays[1]}[$n], 
							                        $location1,
							                        $location2,
							                        "\n";
						
							if ($matchCount > 0) {
					
								my ($sf1, $sf2) = _getIndels(\@features1, 
								   	       					 \@features2, 
								   	       					 $genome1, 
								   	       					 $genome2,
								   	                         $location1, 
								                             $location2, 
								                             \@subFeatures1, 
								                             \@subFeatures2,
								                             $seqRegion,
								                             $orthologueSeqRegion,
								                             $alignmentStart1,
								                             $alignmentStart2,
								                             $alignmentEnd1,
								                             $alignmentEnd2,
								                             $orientation1,
								                             $orientation2);
								           
								@subFeatures1 = @$sf1;
								@subFeatures2 = @$sf2;
							}
						}
						
						# now record mismatch - either as a new feature or an extension
						# of a current one
						
						my $genome1feature = join ',', $genome1,
						                               $seqRegion, 
						                               $location1,
						                               $orientation1;
						                               
						my $genome2feature = join ',', $genome2,
						                               $orthologueSeqRegion,  
						                               $location2,
						                               $orientation2;
						                               
						my @genomeFeatures = ($genome1feature, $genome2feature); 
						my $del = 0;
						my $ref = 1;                 
	    				
	    				if (${$cigarArrays[1]}[$n] eq 'D') {
	    				
	    					# deletion is on genome 2 not genome 1
	    				
	    					$ref = 0;
	    					$del = 1;
	    				}
	    				
	    				my ($genomeDel, $seqRegionDel, $locationDel, $orientationDel) 
	    					= split ',', $genomeFeatures[$del];
	    					
	    				my ($genomeRef, $seqRegionRef, $locationRef, $orientationRef) 
	    					= split ',', $genomeFeatures[$ref];
	    					
	    				my $oldLocationRef = $locationRef;
	    				$locationRef =_increment($locationRef, $orientationRef);
	    					
	    				print STDERR join ',', 'Feature start:',
	    										${$cigarArrays[0]}[$n],
	    										${$cigarArrays[1]}[$n],
	    										$ref,
	    										$del,
	    										$oldLocationRef,
	    										$locationRef,
	    										"\n";
	    										
	    				if ($ref == 0) {
	    				
	    					$location1 = $locationRef;
	    				
	    				} else {
	    				
	    					$location2 = $locationRef;
	    				}
	    					
	    				if (($current[0] eq '') && ($current[1] eq '')) {
	    					
	    					# new feature
	    						
	    					$current[$del] = join ',', $genomeDel,
	    										       $seqRegionDel,
	    										       $locationDel,
	    										       $locationDel,
	    										       'I';
	    						
	    					$current[$ref] = join ',', $genomeRef, 
	    										       $seqRegionRef,
	    										       $oldLocationRef, 
	    										       $locationRef, 
	    										       'D';
	    					
	    				} else {
	    				
	    					# update existing reference feature with updated location 
	    					# value 
	    					# deletion feature is unchanged
	    					
	    					my ($genome, 
	    						$seqRegion, 
	    						$start, 
	    						$end, 
	    						$type) = split ',', $current[$ref];
	    						
	    					$current[$ref] = join ',', $genome, 
	    										       $seqRegion, 
	    											   $start, 
	    											   $locationRef, 
	    											   $type; 
	    				} 
	    			}
				} # end of the alignment
				
				# close out currently open features (variants)			
							
				if ($current[0] ne '') {
						
					_closeFeatures($current[0], 
					               $current[1], 
					               \@features1, 
					               \@features2, 
					               $orientation1, 
					               $orientation2);
				
				} elsif (($current[0] eq '') && ($current[1] eq '')) {
								
					# alignment ends with a match
					# but this match region may contain indels
				
					print STDERR "Final match $location1 $location2 \n";
					
					my ($sf1, $sf2) = _getIndels(\@features1, 
							                     \@features2, 
							                     $genome1,
							                     $genome2,
							                     $location1, 
							                     $location2, 
							                     \@subFeatures1, 
							                     \@subFeatures2,
							                     $seqRegion,
							                     $orthologueSeqRegion,
							                     $alignmentStart1,
							                     $alignmentStart2,
							                     $alignmentEnd1,
							                     $alignmentEnd2,
							                     $orientation1,
							                     $orientation2);
							   
					@subFeatures1 = @$sf1;
					@subFeatures2 = @$sf2;
				}
				
				# print deletion and indel features for this alignment
				
				for my $feature (@features1) {
				
					$mismatch1 += _swapAndPrint($feature, $orientation1, $geneKey);
				}
				
				for my $feature (@features2) {
				
					$mismatch2 += _swapAndPrint($feature, $orientation2, $geneKey);
				}
				
				for my $feature (@subFeatures1) {
				
					$mismatch1 += _swapAndPrint($feature, $orientation1, $geneKey);
				}
				
				for my $feature (@subFeatures2) {
				
					$mismatch2 += _swapAndPrint($feature, $orientation2, $geneKey);
				}
				
				if ( scalar @features1 != scalar @features2) {
				
					die;
				}
				
				if ( scalar @subFeatures1 != scalar @subFeatures2) {
				
					print STDERR join '*', scalar @subFeatures1, 
					                       scalar @subFeatures2, 
					                       "\n";
					                       
					die;
				}
			}
			
			my $keyText1 = $genome1 . $genome2;
			my $keyText2 = $genome2 . $genome1;
			
			# check that mistmatch length isn't longer than total length of each alignment
			
			if (($mismatch1 > $matchLength{$keyText1}) || 
			    ($mismatch2 > $matchLength{$keyText2})) {
			
				print STDERR join ",", 'PROBLEM', 
									   $matchLength{$keyText1}, 
									   $mismatch1, 
									   $matchLength{$keyText2}, 
									   $mismatch2;
										
				print STDERR "\n";
				die;
			} 
			
			# overall length of mismatches for use in final statistics
			
			$mismatch{$keyText1} += $mismatch1;
			$mismatch{$keyText2} += $mismatch2; 
			   		
		} else {
		
			# no alignments on the right DNA segments
			
			$unsupportedPairs ++;
		}
	}
	
	# looking for fully-supported 1:1:1 triples
	
	if (scalar keys %{$orthologues{$gene}} == 2) {
	
		if ($supportedOrthologues == 2) {
	
			print SUP 'S', "\t", $gene;
		
			for my $orthologue (keys %{$orthologues{$gene}}) {
		
				print SUP "\t", $orthologue;
			}
		
			print SUP "\n";
		
		} else {
		
			print SUP 'U', "\t", $gene, "\n";
		}
	}
}

# statistical summary

print STDERR "Features with sequence: $seqCount Without sequence: $seqNoCount\n";

print STDERR join "*\t", scalar keys %hasOrthologues, 
						 scalar keys %geneKeys, 
						 $supportingAlignments, 
						 $unsupportedPairs, 
						 "\n";
                        
print STDERR "\nMismatch/Match/Percentage\n";

for my $genome (sort keys %globalMatchLength) {

	print STDERR join "\t", $genome,
							$mismatch{$genome}, 
							$globalMatchLength{$genome}, 
							($mismatch{$genome}/$globalMatchLength{$genome});
	print STDERR "\n";
}

print STDERR "Method count\n";

for my $key (%methodCount) {

	print STDERR $key, "\t", $methodCount{$key}, "\n";
}


sub _calculateDistance($$$$$$$$) {
	
	# calculate combined length of both genes covered by alignments
	
	my ($geneStart1,  
	    $geneEnd1, 
		$geneStart2, 
		$geneEnd2, 
	    $alignmentStart1, 
		$alignmentEnd1, 
		$alignmentStart2, 
		$alignmentEnd2) = @_;
			      	
	my ($start1, $start2, $end1, $end2);
			      	
	if ($alignmentStart1 <= $geneStart1) {
			      	
		$start1 = $geneStart1;
			      	
	} else {
			      	
		$start1 = $alignmentStart1;
	}
			      	
	if ($alignmentStart2 <= $geneStart2) {
			      	
		$start2 = $geneStart2;
			      	
	} else {
			      	
		$start2 = $alignmentStart2;
	}
			      	
	if ($alignmentEnd1 <= $geneEnd1) {
			      	
		$end1 = $alignmentEnd1;
			      	
	} else {
			      	
		$end1 = $geneEnd1;
	}
			      	
	if ($alignmentEnd2 <= $geneEnd2) {
			      	
		$end2 = $alignmentEnd2;
			      	
	} else {
			      	
		$end2 = $geneEnd2;
	}
			      	
	return $end1 - $start1 + $end2 - $start2;
}

sub _closeFeatures($$$$$$) {

	# find sequence of our working features once these have been extended to their 
	# maximal extent
	
	my ($current1, $current2, $features1, $features2, $orientation1, $orientation2) = @_;
	
	my ($lGenome1, 
		$lseqRegion1, 
		$lStart1, 
		$lEnd1, 
		$lType1) = split ',', $current1;

	my ($lGenome2, 
		$lseqRegion2, 
		$lStart2, 
		$lEnd2, 
		$lType2) = split ',', $current2;
		    				
	my ($s1, $s2);
		    				
	# reference bases
		    				
	$s1 = _getSequence($lseqRegion1, $lStart1, $lEnd1, $orientation1);
	$s2 = _getSequence($lseqRegion2, $lStart2, $lEnd2, $orientation2);	
		    										  									   
	if (($s1 eq '') ||  ($s2 eq '')) {
		    				
	     $seqNoCount ++;
		    				
	} else {
		    					    						    
		$current1 = join ',', $current1, 
		    				  $lGenome2,
		    				  $lseqRegion2, 
		    				  $lStart2, 
		    				  $lEnd2, 
		    				  $s1, 
		    				  $s2;
		    										  
		push @$features1, $current1;
								
		$current2 = join ',', $current2,
		    				  $lGenome1,
		    				  $lseqRegion1, 
		    				  $lStart1, 
		    				  $lEnd1,  
							  $s2, 
						      $s1;
								                      
		push @$features2, $current2;
		$seqCount ++;
	}
}

sub _getId($$$$$$) {

	# build hashes of all genes linked by orthology relationships (maximal sets, not
	# requiring universal linkage
	
	my ($lseqRegion1, $lStart1, $lEnd1, $lseqRegion2, $lStart2, $lEnd2) = @_;
	my $thisId;
	my $key1 = join ',', $lseqRegion1, $lStart1, $lEnd1;
    my $key2 = join ',', $lseqRegion2, $lStart2, $lEnd2;
    
    if ($idHash{$key1} && $idHash{$key2} && ($idHash{$key1} ne $idHash{$key2})) {
    		
    	# possible (not reciprocity: A1-B1, B2-D2, A1-D2)	
    		
    	$thisId = $idHash{$key1};
    	my $oldId = $idHash{$key2};
    	
    	for my $key (@{$id2key{$oldId}}) {
    	
    		push @{$id2key{$thisId}}, $key;
    		$idHash{$key} = $thisId;
    	}
    	
    	delete $id2key{$oldId};
    	
    } elsif ($idHash{$key1}) {
    				
    	$thisId = $idHash{$key1};
    	$idHash{$key2} = $thisId;
    	push @{$id2key{$thisId}}, $key2;
    				
    } elsif ($idHash{$key2}) {
    				
    	$thisId = $idHash{$key2};
    	$idHash{$key1} = $thisId;
    	push @{$id2key{$thisId}}, $key1;
    				
    } else {
    				
    	$thisId = $id;
    	$idHash{$key1} = $thisId;
    	$idHash{$key2} = $thisId;
    	push @{$id2key{$thisId}}, $key1;
    	push @{$id2key{$thisId}}, $key2;
    	$id ++;
    }
    
    print STDERR join '*', $key1, $key2, $idHash{$key1}, $idHash{$key2};
    print STDERR "\n";
    return $thisId;
}

sub _getIndels($$$$$$$$$$$$$$$$) {
		
	# find variations of same length within both sequences	
	# such differences (not producing gaps in the alignment) are stored within match
	# regions in CIGAR strings.
																									
	my ($f1, 
		$f2, 
		$genome1,
		$genome2,
		$location1, # location of the start of the NEXT mistmatch region on sequence 1
		$location2,  # location of the start of the NEXT mistmatch region on sequence 2
		$subFeatures1, 
		$subFeatures2, 
		$seqRegion, 
		$orthologueSeqRegion,
		$alignmentStart1,
		$alignmentStart2,
		$alignmentEnd1,
		$alignmentEnd2,
		$orientation1,
		$orientation2) = @_;
		
	my @features1 = @$f1;
	my @features2 = @$f2;
	my @subFeatures1 = @$subFeatures1;
	my @subFeatures2 = @$subFeatures2;
	my $features1size = scalar @features1;
	my $features2size = scalar @features2;
	die if $features1size != $features2size; # just a check, these two should match
	
	# we want to analyse the last match region 
	# (i.e. the region since the last mismatch "feature"
	# or alternatively, the region since the start of the alignment)
	# a match region is flanked on either side by the gene start or by  a mismatch in one 
	# sequence only
	# the span of the region is determined by the orientation of the sequences in the 
	# alignment
	# cue some tedious coordinate rearrangement...
	
	print STDERR join '^', 'Alignment region', 
							$genome1,
							$location1, 
							$alignmentStart1,
							$alignmentEnd1,
							$genome2, 
							$location2, 
							$alignmentStart2,
							$alignmentEnd2, "\n";
	
	my ($start1, $start2);
	
	if ($features1size == 0) {
	
		# begin from start of alignment
		
		print STDERR "First feature\n";
	
		if ($orientation1 == 1) {
		
			$start1 = $alignmentStart1;
		
		} elsif ($orientation1 == -1) {
		
			$start1 = $alignmentEnd1;
		}
		
		if ($orientation2 == 1) {
		
			$start2 = $alignmentStart2; # previously this was set to 1, was this a bug?
		
		} elsif ($orientation2 == -1) {
		
			$start2 = $alignmentEnd2;
		}
	
	} else {
	
		my $lastFeature1 = $features1[($features1size - 1)];
		
		my ($lGenome1, 
		    $lseqRegion1, 
		    $lStart1, 
		    $lEnd1, 
		    $lType, 
		    $lGenome2, 
		    $lseqRegion2) = split ',', $lastFeature1;
		
		# start the new feature from where the last feature ended
				
		if ($orientation1 == 1) {
		
			$start1 = $lEnd1 + 1;
				
		} elsif ($orientation1 == -1) {
		
			$start1 = $lEnd1 - 1;
		}
	
		my $lastFeature2 = $features2[($features2size - 1)];	
		
		my ($lGenome2, 
		    $lseqRegion2, 
		    $lStart2, 
		    $lEnd2, 
		    $lType, 
		    $lGenome1, 
		    $lseqRegion1) = split ',', $lastFeature2;
		
		if ($orientation2 == 1) {
		
			$start2 = $lEnd2 + 1;
		
		} elsif ($orientation2 == -1) {
		
			$start2 = $lEnd2 - 1;
		}
		
		print STDERR join ',', 'Lends:', $lEnd1, $lEnd2, "\n";
	}
	                       
	# where are we in each genome wrt the end of the last mistmatch feature
	# (or the start of both sequences, as we will be in the case of the first match 
	# feature)?    
	# in the second case, we may need to trim one of the sequences appropriately in order
	# to confine ourselves to the true match region                   
	
	my ($length1, $length2);
	
	if ($orientation1 == 1) {
	
		$length1 = $location1 - $start1 + 1;
		
	} elsif ($orientation1 == -1) {	
		
		$length1 = $start1 - $location1 + 1;
	}
	
	if ($orientation2 == 1) {
	
		$length2 = $location2 - $start2 + 1;
		
	} elsif ($orientation2 == -1) {	
		
		$length2 = $start2 - $location2 + 1;
	}
	
	print STDERR join ' ', "Initial coordinates:", $start1, $length1, $location1, $start2, $length2, $location2, "\n";
	
	# what if 2 alignment regions are of different lengths to one another?
	# e.g. 20-50 on genome 1 matches 30-61 of genome 2
	# what follows is an adjustment of the start of the alignment block in the 
	# longest sequence, which I don't understand...
				
	if ($length1 < $length2) {
		
		if ($orientation2 == 1) {
		
			$start2 = $location2 - $length1 + 1;
		
		} elsif ($orientation2 == -1) {	
		
			$start2 = $location2 + $length1 - 1;
		}
					
	} elsif ($length2 < $length1) {
		
		if ($orientation1 == 1) {
		
			$start1 = $location1 - $length2 + 1;
		
		} elsif ($orientation1 == -1) {	
		
			$start1 = $location1 + $length2 - 1;
		}
	}	
	
	print STDERR join ' ', 'Final coordinates:', $start1, $length1, $location1, $start2, $length2, $location2, "\n";
	
	# end tedious coordination manipulation, now fetch the sequences		
					
	my ($sequence1) = _getSequence($seqRegion, $start1, $location1, $orientation1);	
	
	print STDERR $sequence1, "\n";
	
	my ($sequence2) = 
		_getSequence($orthologueSeqRegion, $start2, $location2, $orientation2);	
		
	print STDERR $sequence2, "\n";
	
	die if length $sequence1 ne length $sequence2;
	
	if (($sequence1 eq '') || ($sequence2 eq '')) {
	
		# alignment exists within a sequence gap
		# this is currently not configured to crash the program, but it seems very odd
	
		$seqNoCount ++;
		return(\@subFeatures1, \@subFeatures2);;
	
	} else {
	
		$seqCount ++;
	}
	
	my @seq1 = split '', $sequence1;
	my @seq2 = split '', $sequence2;
		
	my $length = scalar @seq1;
	my $currentSub1 = '';
	my $currentSub2 = '';
	my $currentStart1 = '';
	my $currentStart2 = '';
					
	for (my $o = 0; $o < $length; $o ++) {
					
		if ($seq1[$o] eq $seq2[$o]) {
						
			if ($currentSub1 ne '') {
				
				# which way we count is determined by the orientation of the scaffold in
				# the alignments, not by the orientation of the sequence-level contigs
				
				print STDERR "SNP length $currentStart1 $currentSub1 " .
				             "$currentStart2 $currentSub2\n";	
				
				# arbiratry limit, don't know the actual limit imposed by the alignment
				# algorithms we're using
				             	
				die if length $currentSub1 > 2000;
				
				my $end1 =_increment($currentStart1, 
				                     $orientation1, 
				                     ((length $currentSub1) - 1));
				
				my $end2 = _increment($currentStart2, 
				                      $orientation2, (
				                      (length $currentSub2) - 1));
				
				if ($end1 < $currentStart1) {
				
					my $tmp = $end1;
					$end1 = $currentStart1;
					$currentStart1 = $tmp;
				}
				
				if ($end2 < $currentStart2) {
				
					my $tmp = $end2;
					$end2 = $currentStart2;
					$currentStart2 = $tmp;
				}
								
    			if (($currentStart1 < 1) ||
    			    ($end1 < 1) ||
    			    ($currentStart2 < 1) ||
    			    ($end2 < 1)) {
    			    
    			    die "Mismtatching lengths: $currentStart1 $start1 $end1 " .
    			        "$currentStart2 $start2 $end2\n";
    			}
    			
    			print STDERR "Revised SNP length $currentStart1 $currentSub1 " .
				             "$currentStart2 $currentSub2 $end1 $end2\n";	
    										
				my $currentFeature1 = join ',',
										$genome1, 
									    $seqRegion, 
										$currentStart1, 
										$end1,
										'S', 
									    $genome2, 
									    $orthologueSeqRegion,
										$currentStart2,
										$end2,
										$currentSub1,
										$currentSub2;
								
				my $currentFeature2 = join ',',
										$genome2, 
									    $orthologueSeqRegion, 
									    $currentStart2,
									    $end2,
										'S', 
										$genome1, 
									    $seqRegion,
										$currentStart1,
										$end1,
										$currentSub2,
									    $currentSub1;			           
														           
				push @subFeatures1, $currentFeature1;
				push @subFeatures2, $currentFeature2;
				$currentSub1 = '';
				$currentSub2 = '';
				$currentStart1 = '';
				$currentStart2 = '';
			}
			
		} elsif ($currentSub1 ne '') {
				
			# extension to previous indel
				
			$currentSub1 .= $seq1[$o];
			$currentSub2 .= $seq2[$o];
						
		} else {
		
			# new indel
				
			$currentSub1 = $seq1[$o];
			$currentSub2 = $seq2[$o];
			
			# current start is the start of the current feature
			
			if ($orientation1 == 1) {
			
				$currentStart1 = $start1 + $o;
			
			} else {
			
				$currentStart1 = $start1 - $o;
			}
			
			if ($orientation2 == 1) {
				
				$currentStart2 = $start2 + $o;
			
			} else {
			
				$currentStart2 = $start2 - $o;
			}
		}
	}
	
	return (\@subFeatures1, \@subFeatures2);
}

sub _getSequence($$$$) {

	# retrieve sequence from database
	
	my ($seqRegion, $start, $end, $orientation) = @_;
	
	if ($start > $end) {
	
		if ($orientation == -1) {
	
			my $tmp = $end;
			$end = $start;
			$start = $tmp;
		
		} else {
		
			print STDERR join "*:", $seqRegion, $start, $end, $orientation, "\n";
			die;
		}
	}
	
	print STDERR join '+', $seqRegion, $start, $end, $orientation, "\n";
	
	my $slice = $slice_adaptor->fetch_by_region($seqRegion2coordSystem{$seqRegion}, 
									            $seqRegion, 
									            $start, 
									            $end, 
									            $orientation);
	
	return $slice->seq();	
}

sub _increment($$$) {

	# increment locations, either up or down (dependent on orientation of alignments)

	my ($location, $increment, $quantity) = @_;
	$quantity = 1 if $quantity eq '';
	$quantity = $quantity * $increment;
	$location += $quantity;
	return $location;
}

sub _removeRedundantGroups($) {

	# there appears to be a certain amount of crap in the DB 
	# (see EG-2486, EG-2505)
	# so this bit of code is a hack to weed out redundant alignments within groups 
	# (identical or subsumed alignments)
	# note that both these JIRA issues have allegedly been solved so this section
	# of code can probably be removed

	my ($alignments) = @_;
	my @alignments = @$alignments;
	my $alignmentCount = scalar @alignments;
	my %drop;
	
	# consider all alignments pairwise
			
	for (my $n = 0; $n < $alignmentCount; $n ++) {
				
		for (my $o = ($n + 1); $o <  scalar @alignments; $o ++) {
			
			# if one group is a subset of another then mark it for deletion
			
			my $a = $alignments[$n];
			my $b = $alignments[$o];
			my @aResults = split ',', $a;
			my @bResults = split ',', $b;
					
			if (($aResults[1] eq $bResults[1]) && 
				($bResults[3] >= $aResults[3]) && 
				($bResults[5] <= $aResults[5])) {
					
				$drop{$o} = 1;
					
			} elsif (($bResults[1] eq $aResults[1]) && 
					 ($aResults[3] >= $bResults[3]) && 
					 ($aResults[5] <= $bResults[5])) {
					
				$drop{$n} = 1;
			}
		}
	}
			
	my @use;
			
	for (my $n = 0; $n < $alignmentCount; $n ++) {
			
		if (! $drop{$n}) {
				
			# alignments we want to keep
				
			push @use, $alignments[$n];
		}
	}
	
	return \@use;
}

sub _reverseSequence ($) {

	my ($sequence) = @_;
	$sequence =~ tr/ACGT/TGCA/;
	$sequence = reverse $sequence;
	return $sequence;
}

sub _reverseAlignments ($$$) {

    # we want a data structure which, for each position we want a value per sequence; 
    # whereas, in the cigar lines, for each sequence we have a value per position 
    # i.e. we need to reverse the key hierarchy on a double hash
	# also, expand representation to decompress runs i.e. M2 -> MM

	my ($genome1, $genome2, $cigars) = @_;
	print STDERR Dumper $cigars;
	my @cigars = @$cigars;
	my @cigarArrays;
	my @lengths;
	
	for (my $n = 0; $n < 2; $n ++) {
			
		my $cigar = $cigars[$n];
		print STDERR "CIGAR $cigar\n";
		my @cigarArray;
					
		while (length $cigar > 0) {
				
			my ($number, $letter, $rest) = ($cigar =~/(\d*?)(\D)(.*)/);
			$number = 1 if $number eq '';				
				
			for (my $m = 0; $m < $number; $m ++) {
				
				push @cigarArray, $letter;
							
				if ($letter eq 'M') {
							
					my $keyText;
								
					if ($n == 0) {
									
						$keyText = $genome1 . $genome2; 
								
					} elsif ($n == 1) {
								
						$keyText = $genome2 . $genome1; 
					}
								
					$globalMatchLength{$keyText} ++;
					$matchLength{$keyText} ++;
				}
			}
				
			$cigar = $rest;
		}
				
		push @lengths, scalar @cigarArray;
		push @cigarArrays, \@cigarArray;
	}
			
	if ($lengths[0] ne $lengths[1]) {
			
		die "Lengths do not map\n";
	}
	
	return ($lengths[0], \@cigarArrays);

}

sub _swapAndPrint($$$) {

	# print out final report
	
	my ($feature, $orientation, $geneKey) = @_;
	my $mismatchLength;
	
	my @components = split /,/, $feature;
	
	my $id = _getId($components[1], 
		            $components[2], 
		            $components[3],
		            $components[6],
    	            $components[7],
    	            $components[8]);
	
	unshift(@components, $id);
	
	# if sequence was reversed in alignment, need to translate back to positive strand
	
	my $nonRefSequence = pop @components;
	my $refSequence = pop @components;
	
	if ($orientation == -1) {
	
		$refSequence = _reverseSequence($refSequence);
		$nonRefSequence = _reverseSequence($nonRefSequence);
	} 
	
	push @components, $refSequence;
	push @components, $nonRefSequence;
	my $revisedFeature = join ',', @components;
	print $geneKey, ',' if $opt_test;
	print $revisedFeature, "\n";
	
	# record bases on reference not matched by non-reference
	
	if (($components[5] eq 'D') || ($components[5] eq 'S')) {
	
		my $matchLength = $components[4] - $components[3];
		
		if ($matchLength < 0) {
		
			$matchLength = $matchLength * -1;
		}
		
		$matchLength ++ if $components[5] eq 'S';
		$mismatchLength += $matchLength;
		
		print STDERR join ',', 'MISMATCH', 
								$components[1], 
								$matchLength, 
								$mismatchLength, "\n";
	}
	
	print STDERR "RETURNING $mismatchLength\n";
	return $mismatchLength;
}

sub get_segment2segment_sql {
    my ($mlss_id, $dbId, $egRelease) = @_;

    if ($egRelease > 32) {

            $query2 ="
    select distinct * from
    (select b.group_id,
            d.name, 
        d2.name AS a, 
        g1.dnafrag_start, 
        g2.dnafrag_start AS b, 
        g1.dnafrag_end, 
        g2.dnafrag_end AS c, 
        g1.dnafrag_strand,
        g2.dnafrag_strand AS d,
        g1.cigar_line, 
        g2.cigar_line AS e,
        mlss.method_link_id
    FROM dnafrag d, genomic_align g1, genomic_align g2, dnafrag d2, genomic_align_block b,
        method_link_species_set mlss
    where d.dnafrag_id = g1.dnafrag_id 
    and g2.genomic_align_block_id = g1.genomic_align_block_id
    and d2.dnafrag_id = g2.dnafrag_id
    and b.genomic_align_block_id = g1.genomic_align_block_id
    AND mlss.method_link_species_set_id = b.method_link_species_set_id
    AND mlss.method_link_species_set_id = $mlss_id
    AND d.name != d2.name";

    } else {

            $query2 = "
    select distinct * from
    (select b.group_id,
            d.name, 
        d2.name AS a, 
        g1.dnafrag_start, 
        g2.dnafrag_start AS b, 
        g1.dnafrag_end, 
        g2.dnafrag_end AS c, 
        g1.dnafrag_strand,
        g2.dnafrag_strand AS d,
        g1.cigar_line, 
        g2.cigar_line AS e,
        mlss.method_link_id
    FROM dnafrag d, genomic_align g1, genomic_align g2, dnafrag d2, genomic_align_block b,
        method_link_species_set mlss
    where d.dnafrag_id = g1.dnafrag_id 
    and g2.genomic_align_block_id = g1.genomic_align_block_id
    and d2.dnafrag_id = g2.dnafrag_id
    and b.genomic_align_block_id = g1.genomic_align_block_id
    AND d.genome_db_id = $dbId
    AND d2.genome_db_id = $dbId
    AND d.name != d2.name
    AND mlss.method_link_species_set_id = b.method_link_species_set_id";
    }

    return $query2;
}



__END__

# running the program and getting NR output

bsub -e stderr.txt -q production-rh7 -M 64000 -R "rusage[mem=64000]" 'perl /homes/pkersey/eg-analysis/src/perl/compara/wheat2.pl -release 34 > wheat9.txt'
sort -u wheat9.txt > wheat10.txt

# number of IDs; events of each type; annotated scaffolds.

m wheat10.txt | perl -e 'while (<>) {@a = split /,/, $_; $h{$a[0]}=1} print scalar keys %h'
m wheat10.txt | perl -e 'while (<>) {@a = split /,/, $_; $h{$a[5]}++} for my $k (keys %h) { print $k, "\t", $h{$k}, "\n"}'
m wheat10.txt | perl -e 'while (<>) {@a = split /,/, $_; $h{$a[2]} = 1; $h{$a[7]} = 1} print scalar keys %h, "\n"'

m wheat10.txt | perl -e 'while (<>) {@a = split /,/, $_; if (@a[1] eq "A") { $h{$a[5]}++}} for my $k (keys %h) { print $k, "\t", $h{$k}, "\n"; $s += $h{$k}} print $s, "\n"'
m wheat10.txt | perl -e 'while (<>) {@a = split /,/, $_; if (@a[1] eq "B") { $h{$a[5]}++}} for my $k (keys %h) { print $k, "\t", $h{$k}, "\n"; $s += $h{$k}} print $s, "\n"'
m wheat10.txt | perl -e 'while (<>) {@a = split /,/, $_; if (@a[1] eq "D") { $h{$a[5]}++}} for my $k (keys %h) { print $k, "\t", $h{$k}, "\n"; $s += $h{$k}} print $s, "\n"'

m wheat10.txt | perl -e 'while (<>) {@a = split /,/, $_; if (@a[1] eq "A") { $h{$a[5] .  $a[6]}++; $s{$a[6]} ++}} for my $k (keys %h) { print $k, "\t", $h{$k}, "\n"} for my $k (keys %s) { print $k, "\t", $s{$k}, "\n"}'
m wheat10.txt | perl -e 'while (<>) {@a = split /,/, $_; if (@a[1] eq "B") { $h{$a[5] .  $a[6]}++; $s{$a[6]} ++}} for my $k (keys %h) { print $k, "\t", $h{$k}, "\n"} for my $k (keys %s) { print $k, "\t", $s{$k}, "\n"}'
m wheat10.txt | perl -e 'while (<>) {@a = split /,/, $_; if (@a[1] eq "D") { $h{$a[5] .  $a[6]}++; $s{$a[6]} ++}} for my $k (keys %h) { print $k, "\t", $h{$k}, "\n"} for my $k (keys %s) { print $k, "\t", $s{$k}, "\n"}'

# transition : transversion ratio

m wheat10.txt | perl -e 'while (<>) { chomp; (@a) = split /,/, $_; if ((($a[10] eq 'A') && ($a[11] eq 'G')) || (($a[10] eq 'G') && ($a[11] eq 'A')) || (($a[10] eq 'C') && ($a[11] eq 'T')) || (($a[10] eq 'T') && ($a[11] eq 'C'))) { $ti++} elsif ((length $a[10] == 1) && (length $a[11] == 1) && ($a[10] ne 'N') && ($a[11] ne 'N')) { $tv ++}} print $ti, "\t", $tv, "\n"'
m wheat10.txt | perl -e 'while (<>) { chomp; (@a) = split /,/, $_; if (($a[10] ne 'N') && ($a[11] ne 'N') && (length $a[10] == 1) && (length $a[11] == 1)) { $k = $a[10]. $a[11]; $h{$k} ++}} for my $k (keys %h) { print $k, "\t", $h{$k}, "\n"}'

# alignments with Ns...

m wheat10.txt | perl -e 'while (<>) { chomp; if ((/,N+$/) || (/,N+,\w+$/)) { print $_, "\n"}}' > gaps

NOTE: expected 3-4% sequence difference in coding regions, transition: transversion ration >= 2:1
