
## Benchmark analyses

The results and plots resulting from the operations below can be browsed at:

https://docs.google.com/spreadsheets/d/1-Gaw8fP-BOESGO-kF4Jrce1lQT8VR8OhlKF7wCbIi_I/edit?usp=sharing


```

# nrplants v0.3 Jun2020
~/plant_tools/production/core/repeats_pipeline_hive.pl -v 101 -R $p1panreg -P $reptmp -D ensembl_production -s aegilops_tauschii -s arabidopsis_halleri -s arabidopsis_thaliana -s arabis_alpina -s brassica_rapa -s malus_domestica_golden -s oryza_sativa -s prunus_dulcis -s setaria_viridis -s triticum_turgidum -s vitis_vinifera -s trifolium_pratense -s olea_europaea_sylvestris -s brachypodium_distachyon -s helianthus_annuus -s zea_mays -s camelina_sativa -w -n

## compute masked fraction Ensembl production way 
## produce the corresponding BED files, which are moved to folder bed/
perl -lne '$sp=(split(/_core/,$_))[0]; $cmd="perl ~/plant_tools/production/misc_scripts/repeat_feature_summary.pl -reg_file \$p1panreg -species $sp -logic_name repeatmask_nrplants"; print $cmd; system("$cmd")' list.cores > log.nrplants.bed

perl -lne '$sp=(split(/_core/,$_))[0]; $cmd="perl ~/plant_tools/production/misc_scripts/repeat_feature_summary.pl -reg_file \$p1panreg -species $sp -logic_name repeatmask_redat"; print $cmd; system("$cmd")' list.cores > log.redat.bed

perl -lne '$sp=(split(/_core/,$_))[0]; $cmd="perl ~/plant_tools/production/misc_scripts/repeat_feature_summary.pl -reg_file \$p1panreg -species $sp -logic_name dust"; print $cmd; system("$cmd")' list.cores > log.dust.bed

perl -lne '$sp=(split(/_core/,$_))[0]; $cmd="perl ~/plant_tools/production/misc_scripts/repeat_feature_summary.pl -reg_file \$p1panreg -species $sp -logic_name trf"; print $cmd; system("$cmd")' list.cores > log.trf.bed

# wheat pangenome suppl table
perl -lne '$sp=(split(/_core/,$_))[0]; $cmd="perl ~/plant_tools/production/misc_scripts/repeat_feature_summary.pl -reg_file \$s3panreg -species $sp -logic_name repeatdetector"; print $cmd; system("$cmd")' list.cores.wheat > log.wheat.Red.bed

perl -lne '$sp=(split(/_core/,$_))[0]; $cmd="perl ~/plant_tools/production/misc_scripts/repeat_feature_summary.pl -reg_file \$s3panreg -species $sp -logic_name repeatmask_redat"; print $cmd; system("$cmd")' list.cores.wheat > log.wheat.redat.bed


# make pretty for spreadsheet
repeated_fraction]$ perl -ane 'if($F[1] =~ 'repeatmask'){ print "$F[2]\n" }elsif($F[1] =~ 'genome'){ print "$F[0]\t$F[2]\t" }' log.redat.bed
repeated_fraction]$ perl -ane 'if($F[1] =~ 'repeatmask'){ print "$F[2]\n" }elsif($F[1] =~ 'genome'){ print "$F[0]\t$F[2]\t" }' log.nrplants.bed
tail -n 1 ../Red_minimap2/*/log.txt 

# get Red repeats in BED format
perl -lne '$sp=ucfirst((split(/_core/,$_))[0]); print $sp' list.cores > list.Red

REDIR=../Red_minimap2/
while read sp; do
	perl -lane 'printf("%s\t%d\t%d\n",$F[0],$F[1]-1,$F[2])' ${REDIR}/${sp}/rpt/*.tsv | bedtools sort > bed/${sp,,}.Red.bed
done < list.Red


## gene/exon overlap stats
## http://quinlanlab.org/tutorials/bedtools/bedtools.html

# get gene & exons coords in BED format
# protein_coding genes

while read core; do
	pl1 $core -Ne "SELECT name,seq_region_start-1,seq_region_end,stable_id FROM gene INNER JOIN seq_region USING (seq_region_id) WHERE biotype='protein_coding' ORDER BY name, seq_region_start" > bed/$core.bed
done < list.cores

while read core; do
	pl1 $core -Ne "SELECT sr.name, ex.seq_region_start-1,ex.seq_region_end,g.stable_id FROM exon ex INNER JOIN seq_region sr USING (seq_region_id) INNER JOIN gene g USING (seq_region_id) INNER JOIN exon_transcript et WHERE g.biotype='protein_coding' AND g.canonical_transcript_id=et.transcript_id AND et.exon_id=ex.exon_id ORDER BY sr.name,ex.seq_region_start" > bed/$core.exon.bed
done < list.cores

# now up & downstream 500 stretches
# Note1: strand does not matter as we don't distinguish up from down
# Note2: substract genes coords to avoid nearby loci
while read core; do
	pl1 $core -Ne "SELECT name,IF(seq_region_start>499, seq_region_start-500, 0) as seq_region_start,seq_region_start FROM gene INNER JOIN seq_region USING (seq_region_id) WHERE biotype='protein_coding' ORDER BY name, seq_region_start" > bed/tmp
	pl1 $core -Ne "SELECT name,seq_region_end-1,seq_region_end+500 FROM gene INNER JOIN seq_region USING (seq_region_id) WHERE biotype='protein_coding' ORDER BY name, seq_region_start" >> bed/tmp
	bedtools sort -i bed/tmp > bed/$core.updown500.raw.bed
	bedtools subtract -sorted -a bed/$core.updown500.raw.bed -b bed/$core.bed > bed/tmp
	bedtools sort -i bed/tmp > bed/$core.updown500.bed
	
	# FASTA sequence could be produced with something like
	#-Ne "SELECT CONCAT('>',name),'\n\r',sequence FROM seq_region INNER JOIN dna USING (seq_region_id)"
done < list.cores


## compute overlaps

# ERROR: Sort order was unspecified, and file bed/aegilops_tauschii_core_48_101_3.bed is not sorted lexicographically.
bedtools sort -i bed/aegilops_tauschii_core_48_101_3.bed 
bedtools sort -i bed/aegilops_tauschii.repeatmask_nrplants.bed
bedtools sort -i bed/aegilops_tauschii.repeatmask_redat.bed
bedtools sort -i bed/aegilops_tauschii.dust.bed
bedtools sort -i bed/aegilops_tauschii.trf.bed

# 1st genes

MINOVER=50

while read -r col1 col2; do	
	red=$(bedtools intersect -a bed/${col1}.bed -b bed/${col2}.Red.bed -sorted -wo | perl -lane '$over+=$F[7]; END{print $over}')
	nrplants=$(bedtools intersect -a bed/${col1}.bed -b bed/${col2}.repeatmask_nrplants.bed -sorted -wo | perl -lane '$over+=$F[7]; END{print $over}')
	redat=$(bedtools intersect -a bed/${col1}.bed -b bed/${col2}.repeatmask_redat.bed -sorted -wo | perl -lane '$over+=$F[7]; END{print $over}')
	genes=$(bedtools intersect -a bed/${col1}.bed -b bed/${col1}.bed -sorted -wo | perl -lane '$over+=$F[8]; END{print $over}')
	read tgenes filename <<< $(wc -l bed/${col1}.bed)
	printf "%s\t%d\t%d\t%d\t%s\t%d\n" $col2 $tgenes $genes $redat $nrplants $red

	# save overlapped genes
	bedtools intersect -a bed/${col1}.bed -b bed/${col2}.Red.bed -sorted -wo > bed/${col2}.Red.genes.bed
	bedtools intersect -a bed/${col1}.bed -b bed/${col2}.repeatmask_nrplants.bed -sorted -wo > bed/${col2}.repeatmask_nrplants.genes.bed
	bedtools intersect -a bed/${col1}.bed -b bed/${col2}.repeatmask_redat.bed -sorted -wo > bed/${col2}.repeatmask_redat.genes.bed

	# MINOVER is required
	perl -lane 'print $F[3] if($F[7] > ${MINOVER})' bed/${col2}.Red.genes.bed | uniq > bed/${col2}.Red.genes.ovl${MINOVER}.list
    perl -lane 'print $F[3] if($F[7] > ${MINOVER})' bed/${col2}.repeatmask_nrplants.genes.bed | uniq > bed/${col2}.repeatmask_nrplants.genes.ovl${MINOVER}.list
    perl -lane 'print $F[3] if($F[7] > ${MINOVER})' bed/${col2}.repeatmask_redat.genes.bed | uniq > bed/${col2}.repeatmask_redat.genes.ovl${MINOVER}.list

	# tally Pfam domain occurrences
	if [ ! -f "pfam/${col2}.repeatmask_redat.genes.ovl${MINOVER}.enrich.tsv" ];
	then
		while read ln; do fgrep "$ln" pfam/${col2}.tsv; done < bed/${col2}.Red.genes.ovl${MINOVER}.list | cut -f 2 | sort | uniq -c | perl -lane 'if(/PF/){ print "$F[1]\t$F[0]" }' > pfam/${col2}.Red.genes.ovl${MINOVER}.tsv
		while read ln; do fgrep "$ln" pfam/${col2}.tsv; done < bed/${col2}.repeatmask_nrplants.genes.ovl${MINOVER}.list | cut -f 2 | sort | uniq -c | perl -lane 'if(/PF/){ print "$F[1]\t$F[0]" }' > pfam/${col2}.repeatmask_nrplants.genes.ovl${MINOVER}.tsv
		while read ln; do fgrep "$ln" pfam/${col2}.tsv; done < bed/${col2}.repeatmask_redat.genes.ovl${MINOVER}.list | cut -f 2 | sort | uniq -c | perl -lane 'if(/PF/){ print "$F[1]\t$F[0]" }' > pfam/${col2}.repeatmask_redat.genes.ovl${MINOVER}.tsv

		# compute enrichment
		Rscript pfam/enrich.R pfam/${col2}.Red.genes.ovl${MINOVER}.tsv pfam/${col2}.tsv.count.tsv pfam/${col2}.Red.genes.ovl${MINOVER}.all.tsv
		perl -lane 'print if($F[2] < 0.05)' pfam/${col2}.Red.genes.ovl${MINOVER}.all.tsv > pfam/${col2}.Red.genes.ovl${MINOVER}.enrich.tsv
		Rscript pfam/enrich.R pfam/${col2}.repeatmask_nrplants.genes.ovl${MINOVER}.tsv pfam/${col2}.tsv.count.tsv pfam/${col2}.repeatmask_nrplants.genes.ovl${MINOVER}.all.tsv
		perl -lane 'print if($F[2] < 0.05)' pfam/${col2}.repeatmask_nrplants.genes.ovl${MINOVER}.all.tsv > pfam/${col2}.repeatmask_nrplants.genes.ovl${MINOVER}.enrich.tsv
		Rscript pfam/enrich.R pfam/${col2}.repeatmask_redat.genes.ovl${MINOVER}.tsv pfam/${col2}.tsv.count.tsv pfam/${col2}.repeatmask_redat.genes.ovl${MINOVER}.all.tsv
		perl -lane 'print if($F[2] < 0.05)' pfam/${col2}.repeatmask_redat.genes.ovl${MINOVER}.all.tsv > pfam/${col2}.repeatmask_redat.genes.ovl${MINOVER}.enrich.tsv
	else
		continue
	fi

done < list.cores.sp > log.genes


# now exons
bedtools sort -i bed/aegilops_tauschii_core_48_101_3.exon.bed

while read -r col1 col2; do
    red=$(bedtools intersect -a bed/${col1}.exon.bed -b bed/${col2}.Red.bed -sorted -wo | perl -lane '$over+=$F[7]; END{print $over}')
    nrplants=$(bedtools intersect -a bed/${col1}.exon.bed -b bed/${col2}.repeatmask_nrplants.bed -sorted -wo | perl -lane '$over+=$F[7]; END{print $over}')
	redat=$(bedtools intersect -a bed/${col1}.exon.bed -b bed/${col2}.repeatmask_redat.bed -sorted -wo | perl -lane '$over+=$F[7]; END{print $over}')
	exons=$(bedtools intersect -a bed/${col1}.exon.bed -b bed/${col1}.exon.bed -sorted -wo | perl -lane '$over+=$F[8]; END{print $over}')
	read texons filename <<< $(wc -l bed/${col1}.exon.bed)
	printf "%s\t%d\t%d\t%d\t%s\t%d\n" $col2 $texons $exons $redat $nrplants $red

	# save overlapped exons
    bedtools intersect -a bed/${col1}.exon.bed -b bed/${col2}.Red.bed -sorted -wo > bed/${col2}.Red.exons.bed
	bedtools intersect -a bed/${col1}.exon.bed -b bed/${col2}.repeatmask_nrplants.bed -sorted -wo > bed/${col2}.repeatmask_nrplants.exons.bed
    bedtools intersect -a bed/${col1}.exon.bed -b bed/${col2}.repeatmask_redat.bed -sorted -wo > bed/${col2}.repeatmask_redat.exons.bed

	perl -lane 'print $F[3] if($F[7] > ${MINOVER})' bed/${col2}.Red.exons.bed | uniq > bed/${col2}.Red.exons.${MINOVER}list
	perl -lane 'print $F[3] if($F[7] > ${MINOVER})' bed/${col2}.repeatmask_nrplants.exons.bed | uniq > bed/${col2}.repeatmask_nrplants.exons.${MINOVER}list
	perl -lane 'print $F[3] if($F[7] > ${MINOVER})' bed/${col2}.repeatmask_redat.exons.bed | uniq > bed/${col2}.repeatmask_redat.exons.${MINOVER}list

done < list.cores.sp > log.exons

# now up & downstream
while read -r col1 col2; do
	red=$(bedtools intersect -a bed/${col1}.updown500.bed -b bed/${col2}.Red.bed -sorted -wo | perl -lane '$over+=$F[6]; END{print $over}')
	nrplants=$(bedtools intersect -a bed/${col1}.updown500.bed -b bed/${col2}.repeatmask_nrplants.bed -sorted -wo | perl -lane '$over+=$F[6]; END{print $over}')
	redat=$(bedtools intersect -a bed/${col1}.updown500.bed -b bed/${col2}.repeatmask_redat.bed -sorted -wo | perl -lane '$over+=$F[6]; END{print $over}')
	updown=$(bedtools intersect -a bed/${col1}.updown500.bed -b bed/${col1}.updown500.bed -sorted -wo | perl -lane '$over+=$F[6]; END{print $over}')
	printf "%s\t%d\t%d\t%s\t%d\n" $col2 $updown $redat $nrplants $red
done < list.cores.sp > log.updown500

# -> gene_exon_updown_plot

 
## check #copies, length of repeats, and overlap among methods

# repeat length stats (median) 
while read -r col1 col2; do
	printf "$col2"
	perl -lane 'print $F[2]-$F[1]' bed/${col2}.Red.bed | Rscript -e 'median(scan(file="stdin"))' 2>&1 | perl -ne 'if(/Read (\d+) items/){ print "\t$1" } elsif(/\[1\] (\d+)/){ print "\t$1" }'
	perl -lane 'print $F[2]-$F[1]' bed/${col2}.repeatmask_nrplants.bed | Rscript -e 'median(scan(file="stdin"))' 2>&1 | perl -ne 'if(/Read (\d+) items/){ print "\t$1" } elsif(/\[1\] (\d+)/){ print "\t$1" }'
	perl -lane 'print $F[2]-$F[1]' bed/${col2}.repeatmask_redat.bed | Rscript -e 'median(scan(file="stdin"))' 2>&1 | perl -ne 'if(/Read (\d+) items/){ print "\t$1" } elsif(/\[1\] (\d+)/){ print "\t$1\n" }'

done < list.cores.sp > log.repeat.length

# repeat length stats (N50)
while read -r col1 col2; do
	printf "$col2"
	perl -ane '$l=$F[2]-$F[1]; $TL+=$l; $R{$.}=$l; END{ foreach $s (sort {$R{$b}<=>$R{$a}} keys(%R)){ $t+=$R{$s}; if($t>$TL/2){ print "\t$R{$s}"; exit }}}' bed/${col2}.Red.bed
	perl -ane '$l=$F[2]-$F[1]; $TL+=$l; $R{$.}=$l; END{ foreach $s (sort {$R{$b}<=>$R{$a}} keys(%R)){ $t+=$R{$s}; if($t>$TL/2){ print "\t$R{$s}"; exit }}}' bed/${col2}.repeatmask_nrplants.bed
	perl -ane '$l=$F[2]-$F[1]; $TL+=$l; $R{$.}=$l; END{ foreach $s (sort {$R{$b}<=>$R{$a}} keys(%R)){ $t+=$R{$s}; if($t>$TL/2){ print "\t$R{$s}\n"; exit }}}' bed/${col2}.repeatmask_redat.bed

done < list.cores.sp > log.repeat.N50

# check Pfam enrichment of gene overlap

cut -f 1 pfam/*redat*enrich.tsv | sort | uniq -c | sort -nr | perl -lane 'if(/PF/){ print "$F[1]\t$F[0]" }' > pfam/enrich.redat.tsv
cut -f 1 pfam/*nrplants*enrich.tsv | sort | uniq -c | sort -nr | perl -lane 'if(/PF/){ print "$F[1]\t$F[0]" }' > pfam/enrich.nrplants.tsv
cut -f 1 pfam/*Red*enrich.tsv | sort | uniq -c | sort -nr | perl -lane 'if(/PF/){ print "$F[1]\t$F[0]" }' > pfam/enrich.Red.tsv

# Red repeat overlap vs others
while read -r col1 col2; do
	printf "$col2"
	red=$(bedtools intersect -a bed/${col2}.Red.bed -b bed/${col2}.Red.bed -sorted -wo | perl -lane '$over+=$F[6]; END{print $over}')
	nrplants=$(bedtools intersect -a bed/${col2}.Red.bed -b bed/${col2}.repeatmask_nrplants.bed -sorted -wo | perl -lane '$over+=$F[6]; END{print $over}')
	redat=$(bedtools intersect -a bed/${col2}.Red.bed -b bed/${col2}.repeatmask_redat.bed -sorted -wo | perl -lane '$over+=$F[6]; END{print $over}')
	dust=$(bedtools intersect -a bed/${col2}.Red.bed -b bed/${col2}.dust.bed -sorted -wo | perl -lane '$over+=$F[6]; END{print $over}')
	trf=$(bedtools intersect -a bed/${col2}.Red.bed -b bed/${col2}.trf.bed -sorted -wo | perl -lane '$over+=$F[6]; END{print $over}')
	totdust=$(bedtools intersect -a bed/${col2}.dust.bed -b bed/${col2}.dust.bed -sorted -wo | perl -lane '$over+=$F[6]; END{print $over}')
	tottrf=$(bedtools intersect -a bed/${col2}.trf.bed -b bed/${col2}.trf.bed -sorted -wo | perl -lane '$over+=$F[6]; END{print $over}')
	
	printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" $col2 $red $redat $nrplants $dust $trf $totdust $tottrf

done < list.cores.sp > log.repeat.overlap


## check GC content

# full genome
while read sptop; do
	echo $sptop
	perl -lne 'if(!/^>/){ $SQ=uc($_); while($SQ =~ /([ACTG])/g){ $stat{$1}++; $tot++ } } END{ printf("%s\t%1.1f\n",$ARGV,100*($stat{"G"}+$stat{"C"})/$tot) }' ${REDIR}/$sptop
done < list.toplevel > log.gc

## kmer analysis of up & downstream regions

K=31    # 21 odd, close to tallymer and DUK,2^21 (2G), 16 used by Red
MINCP=20 # min copies of K-mer to be dumped 
REDIR=../Red_minimap2/
while read -r col1 col2 col3; do
 
 	# get overlapping BED intervals   
	bedtools intersect -a bed/${col1}.updown500.bed -b bed/${col2}.Red.bed -sorted > bed/${col1}.updown500.Red.bed
	bedtools intersect -a bed/${col1}.updown500.bed -b bed/${col2}.repeatmask_nrplants.bed -sorted > bed/${col1}.updown500.nrplants.bed
	bedtools intersect -a bed/${col1}.updown500.bed -b bed/${col2}.repeatmask_redat.bed -sorted > bed/${col1}.updown500.redat.bed 

	# extract overlapping sequences
	bedtools getfasta -fi ${REDIR}/$col3 -bed bed/${col1}.updown500.Red.bed > kmer/${col1}.updown500.Red.fna
	bedtools getfasta -fi ${REDIR}/$col3 -bed bed/${col1}.updown500.nrplants.bed > kmer/${col1}.updown500.nrplants.fna
	bedtools getfasta -fi ${REDIR}/$col3 -bed bed/${col1}.updown500.redat.bed > kmer/${col1}.updown500.redat.fna

	# extract canonical kmers
	# bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction 
	~/soft/jellyfish-linux count -C -m $K -s 2G -t 4 kmer/${col1}.updown500.Red.fna -o kmer/${col1}.updown500.Red.${K}mer.jf
	red=$(~/soft/jellyfish-linux stats kmer/${col1}.updown500.Red.${K}mer.jf | perl -lane 'push(@N,$F[1]); END{ print join("\t",@N)}')
	~/soft/jellyfish-linux dump -L $MINCP -c -t kmer/${col1}.updown500.Red.${K}mer.jf | \
		sort -nrk2,2 > kmer/${col1}.updown500.Red.${K}mer.${MINCP}.tsv
	red20=$(perl -lne '$tot++; END{ print $tot}' kmer/${col1}.updown500.Red.${K}mer.${MINCP}.tsv)

	~/soft/jellyfish-linux count -C -m $K -s 2G -t 4 kmer/${col1}.updown500.nrplants.fna -o kmer/${col1}.updown500.nrplants.${K}mer.jf
	nrplants=$(~/soft/jellyfish-linux stats kmer/${col1}.updown500.nrplants.${K}mer.jf | perl -lane 'push(@N,$F[1]); END{ print join("\t",@N)}')
	~/soft/jellyfish-linux dump -L $MINCP -c -t kmer/${col1}.updown500.nrplants.${K}mer.jf | \
		sort -nrk2,2 > kmer/${col1}.updown500.nrplants.${K}mer.${MINCP}.tsv
	nrplants20=$(perl -lne '$tot++; END{ print $tot}' kmer/${col1}.updown500.nrplants.${K}mer.${MINCP}.tsv)

	~/soft/jellyfish-linux count -C -m $K -s 2G -t 4 kmer/${col1}.updown500.redat.fna -o kmer/${col1}.updown500.redat.${K}mer.jf
	redat=$(~/soft/jellyfish-linux stats kmer/${col1}.updown500.redat.${K}mer.jf | perl -lane 'push(@N,$F[1]); END{ print join("\t",@N)}')
	~/soft/jellyfish-linux dump -L $MINCP -c -t kmer/${col1}.updown500.redat.${K}mer.jf | \
		sort -nrk2,2 > kmer/${col1}.updown500.redat.${K}mer.${MINCP}.tsv
	redat20=$(perl -lne '$tot++; END{ print $tot}' kmer/${col1}.updown500.redat.${K}mer.${MINCP}.tsv)

	printf "%s\t%s\t%s\t%s\t%d\t%d\t%d\n" $col2 $redat $nrplants $red $redat20 $nrplants20 $red20

done < list.cores.sp.toplevel > log.updown500.${K}mer


## check overlap with denovo called Rgenes (NLR-annotator)

while read -r col1 col2; do
	red=$(bedtools intersect -a denovoRgenes/${col2}.nlr.bed.sorted -b bed/${col2}.Red.bed -sorted -wo | perl -lane 'next if($F[9] <= $ENV{'MINOVER'}); $over+=$F[9]; END{print $over}')
	nrplants=$(bedtools intersect -a denovoRgenes/${col2}.nlr.bed.sorted -b bed/${col2}.repeatmask_nrplants.bed -sorted -wo | perl -lane 'next if($F[9] <= $ENV{'MINOVER'}); $over+=$F[9]; END{print $over}')
	redat=$(bedtools intersect -a denovoRgenes/${col2}.nlr.bed.sorted -b bed/${col2}.repeatmask_redat.bed -sorted -wo | perl -lane 'next if($F[9] <= $ENV{'MINOVER'}); $over+=$F[9]; END{print $over}')
	genes=$(bedtools intersect -a denovoRgenes/${col2}.nlr.bed.sorted -b denovoRgenes/${col2}.nlr.bed.sorted -sorted -wo | perl -lane '$over+=$F[12]; END{print $over}')
	read tgenes filename <<< $(wc -l denovoRgenes/${col2}.nlr.bed.sorted)
	printf "%s\t%d\t%d\t%d\t%d\t%d\n" $col2 $tgenes $genes $redat $nrplants $red

done < list.cores.sp > log.Rgenes.$MINOVER

```
