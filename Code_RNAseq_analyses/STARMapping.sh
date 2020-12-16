#!/bin/bash -l
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -N STARmapping
#$ -e ./STARmapping.err




#Author: Rebecca Halbach
#Contact: rebecca.halbach@radboudumc.nl


#requires: STAR 2.7.0 (Dobin et al, Bioinformatics 2013, PMID:23104886 )


#STAR Index must be prepared in advance
#with STAR --runMode genomeGenerate \
#	   --genomeDir <path/to/directory> \
#	   --genomeFastaFiles <path/to/Aedes aegypti AaegL5 chrom fasta> \
#	   --sjdbGTFfile <path/to/Aedes BASEFEATURES_AaegL5.1.gtf> \
#	   --sjdbOverhang 
#Needs to be done only once
#Genome fasta and annotation files can be downloaded from VectorBase: https://www.vectorbase.org/downloads
#Alternatively, you can use the gff3 annotation file, then just add --sjdbGTFtagExonParentTranscript Parent


#This script will align reads to the Ae. aegypti and quantify reads per gene in 2-pass mode, meaning that
#mapping will be done once for all libraries, thenjunctions will be collected and used as annotated
#in a second round of mapping.
#Output of the second round is a bam file with aligned reads, and counts per gene to use in the downstream analysis with DEseq2.



date=$(date +%Y%m%d)

WORKDIR=`pwd`

mkdir "${date}".STARMapping
cd "${date}".STARMapping || { echo "Cannot change to specified directory"; exit 1; }

mkdir tmp/




#############################################
#########General settings####################
#############################################

#path to index (change!)
index=$HOME/AaegL5/STARAaegL5Index
fastqDir="${WORKDIR}"/fastq

echo "This script was run on $(date)
	in 2-pass mode with
	Index: ${index}
	Results will be stored in the folder: $date.STARMapping"

echo -e "\n\n"



#############################################
#########1st-pass mapping####################
#############################################


mkdir 1pass_mapping

cd 1pass_mapping || { echo "Cannot change to directory 1pass_mapping"; exit 1; }

for files in $(ls "${fastqDir}" | sed 's/\_[12].fastq.gz$//' |uniq); do
	
	INFILE="${files}"

	FILENAME=$(echo $(basename "${files}") | sed 's/\.fastq.gz$//')
		
	mkdir -p "${FILENAME}"
	cd "${FILENAME}" || { echo "Cannot change directory"; exit 1; }

	echo "Processing of library %S started at $(date +%H:%M:%S)" "${FILENAME}"

	STAR --runThreadN 24 \
		--genomeDir "${index}" \
		--readFilesIn "${INFILE}"_1.fastq.gz "${INFILE}"_2.fastq.gz \
	       	--readFilesCommand zcat \
		--outSAMtype None \
		--outMultimapperOrder Random \
		--runRNGseed 123 \
		--outSAMmultNmax 1 \
	       	--outSAMstrandField intronMotif \
		--outFileNamePrefix "${FILENAME}"

	cd ..

	
done

cd ..

echo "Done with 1st-pass mapping of all libaries in the directory %s" "${fastqDir}"



#############################################
#########process junctions###################
#############################################

#catenate all Sj (splice junction files) and then remove all junctions on mitochondria (false positives) 
#To use in 2nd-pass mapping

cat */*SJ.out.tab > tmp/SJfiles_cat.tab
sed '/Mt/d' tmp/SJfiles_cat.tab > tmp/SJfiles_cat_wo_Mt.tab





#############################################
#########2-pass mapping######################
#############################################


mkdir 2pass_mapping

cd 2pass_mapping || { echo "Cannot change to directory 2pass_mapping"; exit 1; }

for files in $(ls $WORKDIR/*.fastq.gz | sed 's/\_[12].fastq.gz$//' |uniq); do


	        INFILE="${files}"

		FILENAME=$(echo $(basename "${files}") | sed 's/\.fastq.gz$//')
		
		mkdir -p "${FILENAME}"
		cd "${FILENAME}" || { echo "Cannot change directory"; exit 1; }


	        echo "Processing of library %S started at $(date +%H:%M:%S)" "${FILENAME}"

	        STAR --runThreadN 20 \
	             --genomeDir "${index}" \
	             --readFilesIn "${INFILE}"_1.fastq.gz "${INFILE}"_2.fastq.gz \
	             --readFilesCommand zcat \
	             --outSAMtype BAM SortedByCoordinate \
	             --outMultimapperOrder Random \
	             --runRNGseed 123 \
	             --outSAMmultNmax 1 \
	             --outSAMstrandField intronMotif \
		     --quantMode GeneCounts \
		     --sjdbFileChrStartEnd "${WORKDIR}"/"${date}".STARMapping/tmp/SJfiles_cat_wo_Mt.tab \
	             --outFileNamePrefix "${FILENAME}"
	        

		cd ..
done


echo "Done with 2st-pass mapping of all libaries in the directory %s" "${fastqDir}"

cd "${WORKDIR}"

