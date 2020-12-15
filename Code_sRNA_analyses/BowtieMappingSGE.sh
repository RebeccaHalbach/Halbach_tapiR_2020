#!/bin/bash -l
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -N mapping
#$ -e ./map_err.txt


#Author: Rebecca Halbach
#Contact: rebecca.halbach@radboudumc.nl


#requires: cutadapt v1.14 (Martin, EMBnet.journal 2011, doi: https://doi.org/10.14806/ej.17.1.200)
#requires: Bowtie v0.12.7 (Langmead et al., Genome Biol, 2009, PMID:19261174)

#Bowtie Index must be prepared in advance
#with bowtie-build <path/to/Aedes aegypti AaegL5 chrom fasta> Bowtieindex_AaegL5
#Needs to be done only once
#Genome fasta can be downloaded from vectorBase: https://www.vectorbase.org/downloads


date=$(date +%Y%m%d)

WORKDIR=`pwd`

mkdir "${date}".BowtieMapping
cd "${date}".BowtieMapping || { echo "Cannot change to specified directory"; exit 1; }

mkdir tmp/




#############################################
#########General settings####################
#############################################

#path to bowtie index (change!)
index=$HOME/AaegL5/Bowtieindex_AaegL5/Bowtieindex_AaegL5

#adapter to be clipped from reads
#adapter Illumina
adapter='TGGAATTCTCGGGTGCCAAGG'
#adapter used by Lewis et al, Nat. Ecol & Evol 2017:
#adapter='NNNTGGAATTCTCGGGTGCCAAGG'
#adatper NEBNext
#adapter='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'

#number of max mismatches allowed
mismatch=0




echo "This script was run on $(date)
	with
	Adapter: "${adapter}"
	Index: "${index}"
	Max. number of mismatches: "${mismatch}"
        Results will be stored in the folder: "${date}".BowtieMapping"

echo -e "\n\n"



###############################################
###Adapter clipping & mapping #################
###############################################


#Run for all fastq files in working directory#

for files in "${WORKDIR}"/*fastq.gz; do

	INFILE="${files}"
	FILENAME=$(echo $(basename "${files}") | sed 's/\.fastq.gz$//')
	CLIPPED="${FILENAME}".clipped.fq

	echo "Processing of library %s started at $(date +%H:%M:%S)" "${FILENAME}"


	#clip adapter
	#filter for minimum length of 15, max length 35
	cat "${INFILE}" | \
		gunzip -c | \
		cutadapt -a "${adapter}" -m 15 -M 35  --discard-untrimmed - > tmp/"${CLIPPED}"

	echo -e "\nDone with pre-processing of file %s at $(date +%H:%M:%S)\n\n" "${INFILE}"

	#Map with bowtie
	#options:
	#-k 1 --best: for multimappers, only give one (only the best) alignment
	#-v: number of mismatches
	#-S output is SAM file
	#samtools -F 4: don't print unmapped reads

	OUTBAM="${FILENAME}".bam
	bowtie "${index}" tmp/"${CLIPPED}" --best -k 1 --threads 24 -t -v "${mismatch}" -S |\
		samtools view -Sb -F 4 - | samtools sort - -o "${OUTBAM}"

	OUTBED="${FILENAME}".bed
	bedtools bamtobed -i "${OUTBAM}" > "${OUTBED}"

	echo -e "Done with mapping of file %s (multi- and unique mappers) at $(date +%H:%M:%S)\n\n" "${INFILE}"



	#Map again with only unique reads 
	#-m: only allow reads that can be placed m=1 (== unique) times
	#-v: number of mismatches
	#-S output is SAM file
	#samtools -F 4: don't print unmapped reads

	UNIQBAM="${FILENAME}".unique.bam
	bowtie "${index}" tmp/"${CLIPPED}" -m 1 --threads 8 -t -v "${mismatch}" -S |\
		samtools view -Sb -F 4 - | samtools sort - -o "${UNIQBAM}"

	UNIQBED="${FILENAME}".unique.bed
	bedtools bamtobed -i "${UNIQBAM}" > "${UNIQBED}"
	
	echo -e "Done with mapping of file %s (unique mappers only) at $(date +%H:%M:%S)\n\n" "${INFILE}"
done


cd "${WORKDIR}"



