#!/bin/bash -l
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -N SalmonMapping 
#$ -e ./SalmonMapping.err



#Author: Rebecca Halbach
#Contact: rebecca.halbach@radboudumc.nl


#requires: Salmon v.0.11.3 (Patro et al, Nat Methods 2017, PMID: 28263959)


#A salmon Index must be prepared in advance
#with salmon index -t <path/to/transposon fasta> Salmon_TEfam_index
#Needs to be done only once
#TEfam fasta file was downloaded from https://tefam.biochem.vt.edu/tefam


#This script quantifies reads on transposon sequences
#Output of estimated counts per transposon can be used in the downstream analysis with DEseq2.



date=$(date +%Y%m%d)

WORKDIR=`pwd`

mkdir "${date}".Salmon
cd "${date}".STARMappingSalmon || { echo "Cannot change to specified directory"; exit 1; }



#############################################
#########General settings####################
#############################################



#path to index (change!)
index=$HOME/AaegL5/Salmon_index/Salmon_TEfam_index
#path to fastq files
fastqDir="${WORKDIR}"/fastq
#library type
libType="ISR"




#############################################
#########Quantification######################
#############################################


for f in ${fastqDir}/*.fastq.gz; do 


	file=${f##*/}
	samp=${file%%.*}
	
	echo "Processing sample %s" "${samp}"
	
	salmon quant -i "${index}" \
		-l "${libType}" \
		-1 "${fastqDir}"/"${samp}"_1.fastq.gz \
		-2 "${fastqDir}"/"${samp}"_2.fastq.gz \
		-p 24 \
		-o quants/"${samp}"_quant
	
	echo "Done with quantification of %s"${samp}



done

cd ..

