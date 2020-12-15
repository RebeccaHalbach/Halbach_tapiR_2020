---
title: "Readme"
author: "Rebecca Halbach"
date: "17/04/2020"
output: md_document
---

## Transcriptome analysis in Aedes aegypti

Aim of this project is to analyse mRNAseq and small RNAseq from Aedes aegypti mosquitoes to characterize the function of a highly conserved, piRNA-producing satellite repeat during maternal-to-zygotic transition in mosquitoes. Further information can be found in:    
Halbach et al., A satellite repeat-derived piRNA controls embryonic development of Aedes. Nature 580, 274–277 (2020). 
https://doi.org/10.1038/s41586-020-2159-2  

## Prerequisites  

Bash scripts require the installation of the following tools  

* cutadapt v1.14 (Martin, EMBnet.journal 2011, doi: https://doi.org/10.14806/ej.17.1.200)
* Bowtie v0.12.7 (Langmead et al., Genome Biol, 2009, PMID:19261174)
* samtools 1.9 (Li et al, Bioinformatics 2009, PMID: 19505943)
* bedtools v2.27.1 (Quinlan, Curr. Protoc. Bioinformatics 2014, doi: https://doi.org/10.1002/0471250953.bi1112s47)
* STAR v.2.5.2b (Dobin et al., Bioinformatics 2013, PMID: 23104886)
* Salmon v.0.8.2 (Patro et al., Nat Methods, 2017, PMID: 28263959)
* requires GNU Awk 4.1.4  

R scripts require the installation of the following libraries and packages  

* R 3.6.3 (ran under: Ubuntu 18.04.4 LTS)
* tidyverse 1.2.1
* rtracklayer 1.44.2
* GenomeRanges 1.36.0
* GViz 1.28.1
* GenomicFeatures 1.36.4
* Biostrings 2.52.0
* biomaRt 2.40.4
* AnnotationForge 1.26.0
* AnnotationDbi 1.46.0
* org.Aaegypti.eg.db 0.2
* tximport 1.12.3 
* DESeq2 1.24.0
* RColorBrewer 1.1-2
* latticeExtra 0.6-28
* ggforce 0.3.1 
* venneuler 1.1-0
* clusterProfiler 3.12.0  
  

Annotation files:  
(can be downloaded from VectorBase: https://www.vectorbase.org/downloads)

* Aedes aegypti BASEFEATURES AaegL5.1 (gff3)
* Aedes aegypti REPEATFEATURES AaegL5.1 (gff3)
* Aedes aegypti genome AaegL5 assembly (fasta)
* TEfam transposon sequences (fasta)

## Sequencing libraries

Files used are deposited under BioProject number PRJNA482553 (sequencing of polyA-selected mRNAs from Aag2 cells and Ae. aegypti embryos) and PRJNA594491 (sRNAseq after PIWI immunoprecipitation).  
Above that, publicly available datasets were used from: 

* Lewis et al., Nature Ecol. & Evol. 2018 (PMID: 29203920) (: PRJNA386859 (SRR5961503-SRR5961506)  
* Akbari et al., G3 (Bethesda) 2013 (PMID:23833213): PRJNA209388 (SRR923702, SRR923826, SRR923837, SRR923853, SRR923704)  

## Author

Rebecca Halbach, van Rij lab, Radboudumc Nijmegen, Nijmegen, The Netherlands  
Website: https://vanrijlab.org/  
Contact: rebecca.halbach@radboudumc.nl


