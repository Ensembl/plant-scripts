#!/usr/bin/perl -w

use strict;
use warnings;
use Bio::SeqIO;

# reads GenBank (.gbk) file and extracts contained CDS nucleotide sequences
# to output FASTA file

# based on http://www.eead.csic.es/compbio/material/bioperl/node25.html
# Bruno Contreras Moreira 2019

die "# usage: $0 <file .gbk>\n" unless(@ARGV);
 
my ($gbfile) = @ARGV;

my ($dna_fasta_file,$n_of_CDSs) = extract_CDSs_from_genbank($gbfile);
print "# $0 : store $n_of_CDSs CDSs in $dna_fasta_file\n";


sub extract_CDSs_from_genbank {
	my ($infile) = @_;
	
	my ($out_dna_file,$n_of_CDS,$taxon,$gene,$gi,$crossrefs,$gbaccession,$header,$CDScoords) = ($infile,0);
	$out_dna_file =~ s/\.gbk/_cds.fna/;

	my $in = new Bio::SeqIO(-file => $infile, -format => 'genbank' );
	
	open(DNA,">",$out_dna_file) || 
		die "# extract_CDSs_from_genbank : cannot create $out_dna_file\n";
	
	while( my $seq = $in->next_seq){    
		$taxon='';
		$gbaccession = $seq->accession();   
	   foreach my $f ($seq->get_SeqFeatures){
			if($f->primary_tag() =~ /source/){
				#source          1..72630
            #         /organism="Hordeum vulgare subsp. vulgare"
            #         /mol_type="genomic DNA"
            #         /cultivar="Thibaut"
            #         /sub_species="vulgare"
            #         /db_xref="taxon:112509"
            #         /chromosome="1"
            #         /map="1S(7H)"
				if($f->has_tag('organism')){
					foreach my $element ($f->each_tag_value('organism')){ $taxon = "[$element]"; last; }
				} 
				if($f->has_tag('cultivar')){
					$taxon .= " cultivar:";
					foreach my $element ($f->each_tag_value('cultivar')){ $taxon .= "\"$element\","; }
					chop $taxon;
				}
				if($f->has_tag('map')){
					$taxon .= " map:";
					foreach my $element ($f->each_tag_value('map')){ $taxon .= "$element"; last; }
				}
				if($f->has_tag('chromosome')){
					$taxon .= " chromosome:";
					foreach my $element ($f->each_tag_value('chromosome')){ $taxon .= "$element"; last; }
				}

			} elsif($f->primary_tag() =~ /CDS/){

         	#/gene="NAC47"
            #/note="stress-responsive NAC transcription factor"
            #/codon_start=1
            #/product="NAC protein"
            #/protein_id="AMQ48929.1"

				$gene=$gi=$crossrefs='';
				$CDScoords = $f->spliced_seq(); 
								
				next if(!$CDScoords->{'seq'});

				if($f->has_tag('protein_id')){
					$gi = join(',',sort $f->each_tag_value('protein_id'));
				}

				if($gi eq '' && $f->has_tag('db_xref')){
					$crossrefs = join(',',sort $f->each_tag_value('db_xref'));
					if($crossrefs =~ /(GI\:\d+)/){ $gi = $1 } 
				}
				
				if($f->has_tag('gene')){
					$gene = join(',',sort $f->each_tag_value('gene'));
				} 
				
				if($gene eq '' && $f->has_tag('product')){
					$gene = join(',',sort $f->each_tag_value('product'));
				}
				
				$header = $n_of_CDS.'|'.$gi.'|'.$gene.'|'.$taxon; 
    
				print DNA ">$header\n$CDScoords->{'seq'}\n";
				$n_of_CDS++;	
			}
		}
   }    
   close(DNA);
	
	return ($out_dna_file,$n_of_CDS);
}
