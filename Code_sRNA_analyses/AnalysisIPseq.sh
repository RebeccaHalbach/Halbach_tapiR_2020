#!/bin/bash -l
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -N  analysisIPseq
#$ -e ./analysisIPseq.err

#Author: Rebecca Halbach
#Contact: rebecca.halbach@radboudumc.nl

#requires bedtools v2.27.1 (Quinlan & Hall, Bioinformatics, 2010)
#requires GNU Awk 4.1.4

#Input: bed files from the following libraries
#RDVJ106: input, uninfected
#RDVJ101: empty IP, uninfected
#RDVJ102: Ago3 IP, uninfected
#RDVJ103: Piwi4 IP, uninfected
#RDVJ104: Piwi5 IP, uninfected
#RDVJ105: Piwi6 IP, uninfected





date=$(date +%Y%m%d)
WORKDIR=$(pwd)

printf "Running analysis IP samples Aag2 cells on %s in the directory %s\n" "${date}" "${WORKDIR}"

mkdir "${date}".IPseq
cd "${date}".IPseq || { echo "Cannot change to specified directory"; exit 1; }

mkdir tmp




##Files
#change path!
InputSample=$HOME/IPseq/MappedFiles/RDVJ106.bed
IPSamples=$(ls $HOME/IPseq/MappedFiles/*10[1-5].bed)


########################################
######prepare input sample##############
########################################


#filter piRNA-sized reads, count reads in Input sample

TotalReadsInput=$(\
	wc -l "${InputSample}"  |
	awk '{print $1}')
	
	
awk -v OFS="\t" \
		'{$4= $3-$2}
		{if ($4 >=25 && $4 <= 32)
		{print $0}}' "${InputSample}"  | 
	sort  |
	uniq -c | 
	awk -v n="${TotalReadsInput}" \
		'{$1 = $1/n * 1000000}
		{if ($1 >=1) {print $2"_"$3"_"$4"_"$7, $1}}' |
	sort -k1,1 > tmp/Input.counts
	
	
	
#########################################
#######Analyse IP samples ################
##########################################

header=$(echo "1i\Chrom\tStart\tEnd\tStrand\tCount_input\tCount_IP\tEnr\tCorr_count")
	

for file in ${IPSamples} ; do

	FILENAME=$(echo $(basename "${file}") | sed 's/\.bed//')
	
	
	TotalReads=$(\
		wc -l "${file}"  |
		awk '{print $1}')
	
	#filter piRNA-sized reads, count
	
	awk -v OFS="\t" \
			'{$4= $3-$2}
			{if ($4 >=23 && $4 <= 32)
			{print $0}}' "${file}"  | 
		sort  |
		uniq -c | 
		awk -v n="${TotalReads}" \
			'{$1 = $1/n * 1000000}
			{if ($1 >=1) {print $2"_"$3"_"$4"_"$7, $1}}' |
		sort -k1,1 > tmp/"${FILENAME}".counts


	#Join with input sample, calculate enrichment and counts above background
	
	join -a1 tmp/Input.counts -a2 tmp/"${FILENAME}".counts -o auto -e 0 |
		awk -v OFS="\t" '{$4= ($3+1)/($2+1); $5= ($3-$2) }
			{split($1, b, "_");	$6= b[1]; $7=b[2]; $8=b[3]; $9=b[4]}
			{print $6, $7, $8, $9, $2, $3, $4, $5}' > tmp/"${FILENAME}".enrichment
	
	#filter enriched reads with at least 10 reads in IP sample:
	
	awk '{if ($7 >=2 && $6 >= 10) {print}}' tmp/"${FILENAME}".enrichment |
		sed -e "${header}"> "${FILENAME}".enriched
	
	#extract tapiR1/2 reads
	
	awk -v OFS="\t"  '{if($1==3 &&
			($3==41157596 ||
			$3==41158652 ||
			$3==41158499 ||
			$3==41158348 ||
			$3==41158197 ||
			$3==41158047 ||
			$3==41157894 ||
			$3==41157744 ||
			$3==41157449 ||
			$3==41157301 ||
			$3==41157150 ||
			$3==41157002 ||
			$3==41156851 ||
			$3==41156703 ||
			$3==41156553 ||
			$3==41156405 ||
			$3==41156255 ||
			$3==41156105 ||
			$3==41155955 ||
			$3==41155805)) {print}}' tmp/"${FILENAME}".enrichment | 
		sed -e "${header}" > "${FILENAME}".enrichment.tapiR1
	
	awk -v OFS="\t"  '{if($1==3 &&
			($3==41158739 ||
			$3==41158593 ||
			$3==41158440 ||
			$3==41158289 ||
			$3==41158138 ||
			$3==41157988 ||
			$3==41157835 ||
			$3==41157685 ||
			$3==41157538 ||
			$3==41157390 ||
			$3==41157242 ||
			$3==41157091 ||
			$3==41156943 ||
			$3==41156792 ||
			$3==41156644 ||
			$3==41156494 ||
			$3==41156346 ||
			$3==41156196 ||
			$3==41156046 ||
			$3==41155896 ||
			$3==41155746)) {print}}' tmp/"${FILENAME}".enrichment | 
		sed -e "${header}" > "${FILENAME}".enrichment.tapiR2
	
	printf "Finished analysis IPseq for file %s.\n" "${FILENAME}"
	
	
done

rm -r tmp/
	
	
cd "${WORKDIR}"	
	
	
	
	
	
	
