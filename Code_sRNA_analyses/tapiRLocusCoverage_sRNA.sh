#!/bin/bash -l
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -N  sRNACoverage
#$ -e ./sRNACoverage.err



#Author: Rebecca Halbach
#Contact: rebecca.halbach@radboudumc.nl



#requires bedtools v2.27.1 (Quinlan & Hall, Bioinformatics, 2010)
#requires GNU Awk 4.1.4
#requires samtools 1.9 (Li et al, Bioinformatics 2009, PMID: 19505943) 

 
#Library: RDVJ106 (Aag2 total RNA, input uninfected)


date=$(date +%Y%m%d)
WORKDIR=$(pwd)

printf "Running analysis IP samples Aag2 cells on %s in the directory %s\n" "${date}" "${WORKDIR}"

mkdir "${date}".coverage
cd "${date}".coverage || { echo "Cannot change to specified directory"; exit 1; }

mkdir tmp


#Path to file (change!)
INFILE=$(ls $HOME/IPseq/20190926.BowtieMapping/RDVJ106.bam)

#Annotation
#From VectorBase: gene set AaegL5.1
#Or modified version from PrepareAnnotationFiles.Rmd
#change path!
BASEFEATURES=$HOME/AaegL5/Annotation_files/Aedes-aegypti-LVP_AGWG_BASEFEATURES_AaegL5.1.gff3

awk '/miRNA/' "${BASEFEATURES}" > tmp/miRNAs.gff3


#####Running for each bam file:

for file in ${INFILE} ; do

	FILENAME=$(echo $(basename "${file}") | sed 's/\.bam$//')
	
	
	###Extract reads on tapiR1 locus (as bam)
	samtools view -h "${file}" 3:41155000-41161000 > "${FILENAME}".tapiRLocus.bam
	
	NomiRNAs=$(\
	bedtools intersect -a ${file} -b tmp/miRNAs.gff3 -bed |
		uniq |
		wc -l)

	printf "The number of miRNAs in file %s is %d.\n" "${FILENAME}" "${NomiRNAs}" > "${FILENAME}".miRNAs
	
	###Get BedGraph of piRNA-sized reads (normalized to miRNAs)
	
	NormFac=$(\
	echo "scale=6;${NomiRNAs}/1000000" | bc)
	
	awk 'length($10) >= 23 &&  length($10) <= 32' || $1 ~ /^@/  "${FILENAME}".tapiRLocus | 
		samtools view -bS -| 
		bedtools genomecov -ibam stdin -bga -scale "${NormFac}" > "${FILENAME}".tapiR.coverage.norm.bed
	
done
	
rm -r tmp/

cd "${WORKDIR}"
