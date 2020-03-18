# get cDNA sequences for for best annotated species in EG46/99 (Mar2020)
# https://docs.google.com/spreadsheets/d/1McW5NHzU6zPvYLsvQvD7ysxTxjl7iyHojLUP_WTT6Kw
cd ../../compara

# dicots
./ens_sequences.pl -t cdna -c arabidopsis_thaliana -f arabidopsis_thaliana.cdna.fna
./ens_sequences.pl -t cdna -c glycine_max -f glycine_max.cdna.fna
./ens_sequences.pl -t cdna -c solanum_lycopersicum -f solanum_lycopersicum.cdna.fna
./ens_sequences.pl -t cdna -c populus_trichocarpa -f populus_trichocarpa.cdna.fna
./ens_sequences.pl -t cdna -c vitis_vinifera -f vitis_vinifera.cdna.fna
./ens_sequences.pl -t cdna -c brassica_napus -f brassica_napus.cdna.fna
./ens_sequences.pl -t cdna -c helianthus_annuus -f helianthus_annuus.cdna.fna
./ens_sequences.pl -t cdna -c phaseolus_vulgaris -f phaseolus_vulgaris.cdna.fna
./ens_sequences.pl -t cdna -c medicago_truncatula -f medicago_truncatula.cdna.fna
# monocots
./ens_sequences.pl -t cdna -c 39947 -f oryza_sativa_japonica.cdna.fna
./ens_sequences.pl -t cdna -c zea_mays -f zea_mays.cdna.fna
./ens_sequences.pl -t cdna -c sorghum_bicolor -f sorghum_bicolor.cdna.fna
./ens_sequences.pl -t cdna -c hordeum_vulgare -f hordeum_vulgare.cdna.fna
./ens_sequences.pl -t cdna -c brachypodium_distachyon -f brachypodium_distachyon.cdna.fna

# plus TE libraries repetDB, trep, mips and SINE

contents:

arabidopsis_thaliana.cdna.fna.gz
brachypodium_distachyon.cdna.fna.gz
brassica_napus.cdna.fna.gz
cdna.list
glycine_max.cdna.fna.gz
helianthus_annuus.cdna.fna.gz
hordeum_vulgare.cdna.fna.gz
medicago_truncatula.cdna.fna.gz
mipsREdat_9.3p_ALL.fasta.gz
oryza_sativa_japonica.cdna.fna.gz
phaseolus_vulgaris.cdna.fna.gz
populus_trichocarpa.cdna.fna.gz
repetDB.Mar2020.fna.gz
repetDB.Mar2020.tsv
SINEs.plants.fna.gz
solanum_lycopersicum.cdna.fna.gz
sorghum_bicolor.cdna.fna.gz
TE.list
trep-db_nr_Rel-19.fasta.gz
vitis_vinifera.cdna.fna.gz
zea_mays.cdna.fna.gz
