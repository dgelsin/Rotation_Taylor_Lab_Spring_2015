#!/bin/sh

#PBS -l nodes=1:ppn=16

#PBS -M dgelsin1@jhu.edu

#tophat and cufflinks build with linuxbrew 

#PBS -l walltime=15:00:00

export PATH=/data1/taylor-fs1/linuxbrew/bin:/data1/taylor-fs1/bootstrap/usr/bin:$PATH

cd #PBS_O_WORKDIR

GENES = "/home/dgelsin1/RNA-seq_analysis/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf"
GENOME = "/home/dgelsin1/RNA-seq_analysis/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome"

#ACCESSION_DIRECTORY="cat /home/dgelsin1/RNA-seq_analysis/accession_list.txt | tail -n +${PBS_ARRAY_INDEX} | head -1"

OUTPUT_PATH1="/home/dgelsin1/RNA-seq_analysis/hum_tophat_output/"


TOPHAT_OUTPUT = $OUTPUT_PATH1

/data1/taylor-fs1/linuxbrew/bin/tophat2 -p 12 -G "/home/dgelsin1/RNA-seq_analysis/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf" "/home/dgelsin1/RNA-seq_analysis/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome" "/home/dgelsin1/RNA-seq_analysis/hum_RNA-seq/"$ACCESSION_DIRECTORY"/"$ACCESSION_DIRECTORY"_1" "/home/dgelsin1/RNA-seq_analysis/hum_RNA-seq/"$1"/"$1"_2"


#----------------------------------------------------------------------------
##################################CUFFLINKS##################################
#----------------------------------------------------------------------------
CUFFLINKS_INPUT=$TOPHAT_OUTPUT"*"
CUFFLINKS_OUTPUT=$OUTPUT_PATH2$CUFFLINKS_INPUT"/"

cufflinks -p 8 -o $CUFFLINKS_OUTPUT $TOPHAT_OUTPUT"accepted_hits.bam"
#Assemble transcripts for each sample

1.fastq.bz2

#----------------------------------------------------------------------------
##################################MAKE ASSEMBLIES.TXT##################################
#----------------------------------------------------------------------------
ASSEMBLIES_TXT = $CUFFLINKS_OUTPUT"/transcripts.gtf"

echo $ASSEMBLIES_TXT > assemblies.txt 
#creat assembly file that lists the assembly file for each sample

#----------------------------------------------------------------------------
##################################CUFFMERGE##################################
#----------------------------------------------------------------------------
CUFFMERGE_INPUT = "assemblies.txt"
GENOME_FA = "/home/dgelsin1/RNA-seq_analysis/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome.fa"
cuffmerge -g $GENOME_FILE -s GENOME_FA -p 8 $CUFFMERGE_INPUT
#cuffmerge assemblies to create a single merged transcriptome annotation





###########################################################
###########################################################
######################BETA TESTING STUFF###################
###########################################################
###########################################################

#----------------------------------------------------------------------------
#Alternative way of writing tophat?

#human muscle data
#tophat -p 8 -G hum_genes.gtf -o hM1_tophat_out SRR945190.fastq.bz2 SRR945191.fastq.bz2
#tophat -p 8 -G hum_genes.gtf -o hM2_tophat_out SRR945192.fastq.bz2 SRR94513.fastq.bz2
#tophat -p 8 -G hum_genes.gtf -o hM3_tophat_out SRR945194.fastq.bz2 SRR945195.fastq.bz2
#tophat -p 8 -G hum_genes.gtf -o hM4_tophat_out SRR945196.fastq.bz2 SRR945197.fastq.bz2
#tophat -p 8 -G hum_genes.gtf-o hM5_tophat_out SRR945198.fastq.bz2 SRR945199.fastq.bz2
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945200.fastq.bz2 SRR945201.fastq.bz2

#human kidney data
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945202.fastq.bz2 SRR945203.fastq.bz2
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945204.fastq.bz2 SRR945205.fastq.bz2
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945206.fastq.bz2 SRR945207.fastq.bz2
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945208.fastq.bz2 SRR945209.fastq.bz2
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945210.fastq.bz2 SRR945211.fastq.bz2
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945212.fastq.bz2 SRR945213.fastq.bz2

#human prefrontal cortex data
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945200.fastq.bz2 SRR945201.fastq.bz2
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945200.fastq.bz2 SRR945201.fastq.bz2
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945200.fastq.bz2 SRR945201.fastq.bz2
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945200.fastq.bz2 SRR945201.fastq.bz2
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945200.fastq.bz2 SRR945201.fastq.bz2
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945200.fastq.bz2 SRR945201.fastq.bz2

#human visual cortex data
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945200.fastq.bz2 SRR945201.fastq.bz2
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945200.fastq.bz2 SRR945201.fastq.bz2
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945200.fastq.bz2 SRR945201.fastq.bz2
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945200.fastq.bz2 SRR945201.fastq.bz2
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945200.fastq.bz2 SRR945201.fastq.bz2
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945200.fastq.bz2 SRR945201.fastq.bz2

#human cerebellar cortex data
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945200.fastq.bz2 SRR945201.fastq.bz2
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945200.fastq.bz2 SRR945201.fastq.bz2
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945200.fastq.bz2 SRR945201.fastq.bz2
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945200.fastq.bz2 SRR945201.fastq.bz2
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945200.fastq.bz2 SRR945201.fastq.bz2
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945200.fastq.bz2 SRR945201.fastq.bz2


