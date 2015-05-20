This supplementary dataset contains expression level estimates for lncRNAs and for Ensembl-annotated protein-coding genes. We provide the numbers of uniquely mapped reads, the raw RPKM values as well as the normalized RPKM values. The normalization procedure is the one described in Brawand et al., 2011, which uses the 1000 least-varying genes as a standard to normalize the RPKM values. 

##########################################################################

Several types of data are provided for each species and for each dataset (main dataset or second dataset with only strand specific samples). This folder contains files for the main dataset (non strand-specific samples).

1) *_NbReads_NonSS_LncRNA_ProteinCoding_*  : number of uniquely mapped reads. all samples are treated as unstranded.

2) *_NbReads_SS_LncRNA_ProteinCoding_*  : number of uniquely mapped reads, only for strand-specific samples, with the strand information taken into account. 

3) *_NormalizedRPKM_NonSS_LncRNA_ProteinCoding_*  : normalized RPKM,  all samples are treated as unstranded.

4) *_NormalizedRPKM_SS_LncRNA_ProteinCoding_*  : normalized RPKM, only for strand-specific samples, with the strand information taken into account. 

##########################################################################

In addition, we provide the normalization coefficients for the between-species comparisons, obtained with the same method. To obtain normalized RPKM, the raw RPKM values were divided by the normalization coefficients for each sample. 

Note that the chimpanzee and bonobo data are joined together in a single file named *Pan*.

##########################################################################

Finally, this dataset also contains miRNA expression values for 5 species. After adapter trimming, we mapped reads of length 15-25 on a database of miRNA precursors (Ensembl62). We then counted the number of reads whic map uniquely on each Ensembl-annotated miRNA. For between-tissues comparisons, we resampled the same number of uniquely mapped reads. 

There are 2 files per species:

miRNAs*_AllSamples_Ensembl62_NbReads.txt : number of unique reads per miRNA and per sample.

miRNAs*_AllTissues_Ensembl62_NbReads_Resampled.txt : number of unique reads per miRNA and per tissue, after resampling.

##########################################################################

