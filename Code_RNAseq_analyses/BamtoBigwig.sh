#!/bin/bash -l
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -N BigWigs.err
#$ -e ./BigWigs.err



#Author: Rebecca Halbach
#Contact: rebecca.halbach@radboudumc.nl


#requires bedtools v2.27.1 (Quinlan & Hall, Bioinformatics, 2010)
#requires kentUtils (https://github.com/ucscGenomeBrowser/kent)




date=$(date +%Y%m%d)
WORKDIR=$(pwd)

INFOLDER="${WORKDIR}"/bamFiles
CHROMINFO=$HOME/Annotation/Aedes-aegypti-LVP_AGWG_CHROMOSOMES_AaegL5_chromInfo.txt


printf "Convert bam file to bigwig.\n" "${date}" 


mkdir "${date}".bigwig
cd "${date}".bigwig || { echo "Cannot change to specified directory"; exit 1; }



for file in ${INFOLDER}/*.bam; do

	INFILE="${file}"

	FILENAME=$(echo $(basename ${INFILE}) | sed 's/\Aligned.sortedByCoord.out.bam$//')
	
	
	#separate files for plus and minus strand reads
	genomeCoverageBed -ibam "${INFILE}" \
			-g "${CHROMINFO}"  \
			-bg \
			-split \
			-strand + |
			sort -k1,1 -k2,2n> "${FILENAME}".plus.bedGraph
			
	bedGraphToBigWig "${FILENAME}".plus.bedGraph "${CHROMINFO}"  "${FILENAME}".plus.bw	

	
	genomeCoverageBed -ibam "${INFILE}" \
			-g "${CHROMINFO}"  \
			-bg \
			-split \
			-strand - |
			sort -k1,1 -k2,2n > "${FILENAME}".min.bedGraph

	bedGraphToBigWig "${FILENAME}".min.bedGraph "${CHROMINFO}"  "${FILENAME}".min.bw	
	
	
    echo "Converted file %s to bigwig.\n" "${FILENAME}"
	
	
done

cd "${WORKDIR}"
