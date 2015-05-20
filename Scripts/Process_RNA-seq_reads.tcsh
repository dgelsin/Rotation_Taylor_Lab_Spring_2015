#!/bin/tcsh
#---------------------------------------------------------------------------
# tells the scheduler what queue to submit the job to
PBS -q batch
# requests the max wall time for this job
PBS -l walltime=10:00
#---------------------------------------------------------------------------
#
#    -l nodes=80  yields 80 processors and lets the scheduler decide where to run 
#                 the processes.  This is your best bet for getting your job out of idle
#                 in the fastest possible way as long as you don't need a large number
#                 of processors.  NOTE: you can only specify up to the maximum
#                 current count of nodes on the cluster.
#
PBS -l nodes=80
# this specifies no stderr/stdout merging should be done
PBS -j n
# tells the scheduler that we want email when the job is 
# aborted and when the job terminates normally.
PBS -m ae
# send an email to nobody@somewhere.com with job status changes
# MAKE SURE YOU CHANGE THIS LINE OR REMOVE IT!!
PBS -M dgelsin1@jhu.edu

#---------------------------------------------------------------------------

#USER SPECIFIED VARIABLES
#"User defined output PATH and also name of output folder"
#"User specified PATH to fastq files to process"

#SAMPLE_NAME="SRR945"
INPUT_DATA_FOLDER="/home/dgelsin1/RNA-seq_analysis/hum_RNA-seq/SRX3293"
SUB_INPUT_DATA = "/SRR945"
REF_GENOME = "/home/dgelsin1/RNA-seq_analysis/genes.gtf"
# want to do this: /home/dgelsin1/hum_RNA-seq/SRX3293*/SRR945* where SRX3293 represent a directory (ascending in numerical order ie: SRX29370, SRX29371, ... SRX29399) for each biological replicate RNA-seq sample, 
#followed by the paired end reads for each biological replicate (ie: /SRR945190 + /SRR945191, ... /SRR945248 + /SRR945249)


OUTPUT_PATH1="/home/dgelsin1/RNA-seq_analysis/hum_tophat_output/"
OUTPUT_PATH2="/home/dgelsin1/RNA-seq_analysis/hum_cufflinks_output/"
THREADS=16

#----------------------------------------------------------------------------
##################################TOPHAT##################################
#----------------------------------------------------------------------------
INPUT_FILES=$INPUT_DATA_FOLDER"*"SUB_INPUT_DATA"*"
GENOME_FILE=$REF_GENOME
TOPHAT_OUTPUT = $OUTPUT_PATH1$INPUT_FILES"/"
# input data, two reads per tophat run
cd $INPUT_FILES, $REF_GENOME

tophat -p 12 -G $REF_GENOME -o $TOPHAT_OUTPUT genome $INPUT_FILES
# run tophat on each pair of paired-end reads, it is installed in my local cluster bin (/home/dgelsin1/bin/)


#----------------------------------------------------------------------------
##################################CUFFLINKS##################################
#----------------------------------------------------------------------------
CUFFLINKS_INPUT=$TOPHAT_OUTPUT"*"
CUFFLINKS_OUTPUT=$OUTPUT_PATH2$CUFFLINKS_INPUT"/"

cufflinks -p 8 -o $CUFFLINKS_OUTPUT $TOPHAT_OUTPUT"accepted_hits.bam"
#Assemble transcripts for each sample



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
cuffmerge -g $GENOME_FILE -s genome.fa -p 8 $CUFFMERGE_INPUT
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
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945200.fastq.bz2 SRR945201.fastq.bz2
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945200.fastq.bz2 SRR945201.fastq.bz2
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945200.fastq.bz2 SRR945201.fastq.bz2
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945200.fastq.bz2 SRR945201.fastq.bz2
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945200.fastq.bz2 SRR945201.fastq.bz2
#tophat -p 8 -G hum_genes.gtf -o hM6_tophat_out SRR945200.fastq.bz2 SRR945201.fastq.bz2

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


