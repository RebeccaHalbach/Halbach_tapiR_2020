#!/bin/bash -l
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -N  genomedistr
#$ -e ./genomedistr.err



#Author: Rebecca Halbach
#Contact: rebecca.halbach@radboudumc.nl



#requires bedtools v2.27.1 (Quinlan & Hall, Bioinformatics, 2010)
#requires GNU Awk 4.1.4
#requires samtools 1.9 (Li et al, Bioinformatics 2009, PMID: 19505943) 


#Libraries used:
#SRR5961503 (Ae. aegypti soma, oxidized)  
#SRR5961504 (Ae. aegypti germline, oxidized)  
#SRR5961505 (Ae. aegypti soma)  
#SRR5961506 (Ae. aegypti germline)  
#RDVJ106 (Aag2 total sRNA input, uninfected)  







date=$(date +%Y%m%d)
WORKDIR=$(pwd)

printf "Running piRNA analysis for Aedes aegypti on %s in the directory %s\n" "${date}" "${WORKDIR}"

mkdir "${date}".piRNAsAedesAegypti
cd "${date}".piRNAsAedesAegypti || { echo "Cannot change to specified directory"; exit 1; }

mkdir tmp
mkdir tmp/GFFs





#######Files needed#########
#Annotation file from vectorBase (BaseFeatures AaegL5 https://www.vectorbase.org/downloads)
#BASEFEATURE annotation modified with R script PrepareAnnotationFiles.Rmd
#change path!

BASEFEATURES=$HOME/AaegL5/Annotation_files/Aedes-aegypti-LVP_AGWG_BASEFEATURES_AaegL5.1_modified.gff3
REPEATFEATURES=$HOME/AaegL5/Annotation_files/aedes-aegypti-lvpagwgrepeatfeaturesaaegl5.gff3


#Handle to of bam files:
#Change path!
INFILE=$(ls $HOME/LewisMosquitosRNA/*[0-9].bam)







###################################################
#######Distribution features on genome#############
##################################################

cat "${BASEFEATURES}" "${REPEATFEATURES}"| 
	awk ' /chromosome|supercontig/' |
	sort -k1,1 -k4,5n  > tmp/chrom.info.gff3



#Split protein coding and ncRNA
sort  -k1,1 -k4,5n "${BASEFEATURES}" |
	awk   '{if ($3 ~/CDS/ ) {print $0 > "tmp/GFFs/01_protein_coding"} 
			else if (/prot/ || $9~/nontranslating_CDS/) {print $0 > "tmp/GFFs/12_protUTRIntron"} 
			else if (/intron/) {print $0 > "tmp/GFFs/13_ncRNA_intron"}
			else  if( $3 !~ /lncRNA|pseudogene/) {print $0 > "tmp/GFFs/02_ncRNAs"}}'

#Split Repeatfeatures

sort -k1,1 -k4,5n "${REPEATFEATURES}" |
	awk  '{if (/LTR/) {print $0 > "tmp/GFFs/03_LTR_Retrotransposon"} 
			else if (/SINE|LINE|Penelope/) {print $0 > "tmp/GFFs/04_non_LTR_Retrotransposon"}
			else if (/DNA/) {print $0 > "tmp/GFFs/05_DNA_transposon"}
			else if (/Helitron|MITEs/) {print $0 > "tmp/GFFs/06_otherDNA_Transposons"}
			else if (/Satellite|trf/)  {print $0 > "tmp/GFFs/07_Satellite"}
			else if (/dust/) {print $0 > "tmp/GFFs/08_dust"}
			else if (/repeat|UD|Unknown/) {print $0 > "tmp/GFFs/09_other_unknown_repeats"}}'


#calculate nucleotides covered by a feature
#subsequently these regions are then subtracted from other features to avoid double-counting

	
bedtools merge -i tmp/GFFs/01_protein_coding |
	sort -k 1,1 -k2,3n > tmp/NewFile	 #output is bed format -> zero-based
	
nucl=$(\
	awk '{$4= $3- $2 }
	{sum += $4} END{print sum}' tmp/NewFile) 

echo -e "01_protein-coding\t${nucl}" > summary_genome.tsv
	
for feature in $(find tmp/GFFs/* | tail -n +2) ; do
	
	NAME=$(basename "${feature}")
	
	bedtools subtract -a <(bedtools merge -i "${feature}") \
		-b tmp/NewFile > tmp/"${NAME}".diff
	
	nucl=$(\
		awk '{$4= $3- $2}
		{sum += $4} END{print sum}' tmp/"${NAME}".diff) #bed files -> zero-based
	
	echo -e "${NAME}\t${nucl}" >> summary_genome.tsv
	
	cat tmp/NewFile > tmp/NewFile.tmp
	
	cat tmp/NewFile.tmp tmp/"${NAME}".diff  |
		sort -k1,1 -k2,3n >> tmp/NewFile

done
	
#Regions without any feature
nucl=$(\
	cat "${BASEFEATURES}" "${REPEATFEATURES}" |
	awk '! /chromosome|supercontig|#/' |
	bedtools subtract -a tmp/chrom.info.gff3 -b stdin |
	awk '{$6= $5- $4 +1}							
	{sum += $6} END{print sum}') #output is gff3 -> 1-based 
	
echo -e "14_other_sequences\t${nucl}" >> summary_genome.tsv


#Check if sum of all features and others is equal to genome length	
SumF=$(\
	awk '{sum += $2} END{print sum}' summary_genome.tsv)

LengthGenome=$(\
	awk '{$6= $5- $4 +1}{sum += $6} END{print sum}' tmp/chrom.info.gff3)

if (( "${SumF}" == "${LengthGenome}" )); 
then
   echo "Sum of all features is equal to total genome length."
else
   echo "Something went wrong. Sum of all features is NOT equal to total genome length. Check!"
fi	
	
printf "Finished distribution of genomic features.\n"
	
	
#################################################
#########Distribution of small RNAs###############
#################################################

##Prepare Files

#grep "miRNA"  "${FEATURES} "> tmp/miRNAs.gff3

sort -k1,1 -k4,5n "${REPEATFEATURES}" |
	awk  -v OFS="\t" '{ if ($3 != "chromosome" && $3 != "supercontig") { print } }'  |
	grep -v '#' |
	awk -v OFS="\t" '{split($9, b, ";");  $10= b[1]; $11=b[2]};
		{gsub("class=", "", $11)};
		{print $1, $2, $11, $4, $5, $6, $7, $8, $9} ' |
	cat  - "${BASEFEATURES}" |
	sort -k1,1 -k4,5n > tmp/Annotation


##Running for each bam file:

for file in ${INFILE} ; do

	FILENAME=$(echo $(basename "${file}") | sed 's/\.bam$//')
	
	#Bam to bed and sort by genomic position
	#Count piRNAs (same sequence and length)
	
	bedtools bamtobed -i  "${file}"|
		sort -k1,1 -k2,2n -k3,3n |
		awk '{$4 = "x"; print}' |
		uniq -c |
		awk -v OFS="\t" '{print $2, $3, $4, $5, $1, $7}' > "${FILENAME}".bed
	
	
	#overlap with annotated features
	bedtools intersect -a "${FILENAME}".bed -b tmp/Annotation -wa -wb -sorted |	
		bedtools groupby -g 1,2,3,5,6 -c 9,15 -o distinct,distinct |
		awk '!/tRNA|rRNA/' > tmp/"${FILENAME}".overlapping
	
	#not overlapping with any annotated features
	bedtools intersect -a "${FILENAME}".bed -b tmp/Annotation -v -sorted |
		awk -v OFS="\t" '{print $1, $2, $3, $5, $6, "none", "none"}' > tmp/"${FILENAME}".nonoverlapping
	
	
	#put together 
	cat tmp/"${FILENAME}".overlapping tmp/"${FILENAME}".nonoverlapping |
		sort -k1,1 -k2,2n -k3,3n > tmp/"${FILENAME}".merged
	
	#Extract siRNAs and piRNAs
	awk '{if (($3-$2)<=22 && ($3-$2) >=19) {print}}'  tmp/"${FILENAME}".merged |
		awk -v OFS='\t' '{print $4, $6}'| 
		sort -k2,2 | 
		bedtools groupby -i stdin -g 2 -c 1 -o sum> "${FILENAME}".merged.siRNAs.genomedistribution
	
	awk '{if (($3-$2)>=23 && ($3-$2) <=33) {print}}' tmp/"${FILENAME}".merged|
		awk -v OFS='\t' '{print $4, $6}'| 
		sort -k2,2 | 
		bedtools groupby -i stdin -g 2 -c 1 -o sum > "${FILENAME}".merged.piRNAs.genomedistribution
	
	
	##################################################
	#####Make size profile for tapiR1/2################
	##################################################
	
	#tapir1 (reverse complement of first 23 nt):
	#TAAAACGACCTAGTTTTGAAGAC
	#CTAAAACGACCTATTTTGAAGAC -> tapiR1 repeat with a 1nt deletion
	
	samtools view "${file}" |
		cut -f 10 | 
		awk '/TAAAACGACCTAGTTTTGAAGAC$|CTAAAACGACCTATTTTGAAGAC$/' |
		awk '{ print length }' | 
		sort | 
		uniq -c  |
		awk -v OFS="\t" '{print $2, $1}' > "${FILENAME}".tapiR1.sizes
	
	#tapiR2 (reverse complement of first 23 nt):
	#GAATTTCTAAAACATATCCGAAA
	
	samtools view "${file}" |
		cut -f 10 | 
		awk '/GAATTTCTAAAACATATCCGAAA$/' |
		awk '{ print length }' | 
		sort | 
		uniq -c  |
		awk -v OFS="\t" '{print $2, $1}' > "${FILENAME}".tapiR2.sizes

	
	printf "Finished analysis genome distribution sRNAs for file %s.\n" "${FILENAME}"
	
done	
	
rm -r tmp/

cd "${WORKDIR}"
