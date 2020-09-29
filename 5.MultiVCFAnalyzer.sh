#!/bin/bash

##This is an example of the config script for MultiVCFAnalzyer with parameters
##used in this study.
##By setting min allele freq for homozygous call and min allele freq for hetero call ==, 
##uncertainty in IUPAC calls are avoided and the potential heterozygous site is called as N


java -Xmx32G -jar /path/to/MultiVCFAnalyzer_0-85-1.jar \
NA \
/path/to/reference.fasta \
NA \
/path/to/output/directory \
F \ #whether to include the percentage of reads a given allele is present in in the SNP table e.g. A (70%)
30 \ #Minimal genotyping quality (GATK) - a threshold of which a SNP call falling under is 'discarded'
5 \ #Minimal coverage for base call - the minimum number of a reads that a position must be covered by to be reported
0.8 \ #Minimal allele frequency for homozygous call - the fraction of reads a base must have to be called 'homozygous'
0.8 \ #Minimal allele frequency for heterozygous call - a fraction of which whereby if a call falls above this value, and lower than the homozygous threshold, a base will be called 'heterozygous' and reported with a IUPAC uncertainity code
/path/to/sites/to/exclude.gff \
/data/stonelab/HLP/VCF_files/1126/1126_trimmed_mapped_q37_sorted_rmdup_uniq_RG.vcf \
/data/stonelab/HLP/VCF_files/1262-16/1262-16_trimmed_mapped_q37_sorted_rmdup_uniq_RG.vcf \
/data/stonelab/HLP/VCF_files/2188/2188_trimmed_mapped_q37_sorted_rmdup_uniq_RG.vcf \
/data/stonelab/HLP/VCF_files/2936/2936_trimmed_mapped_q37_sorted_rmdup_uniq_RG.vcf \
/data/stonelab/HLP/VCF_files/2DDS/2DDS_trimmed_mapped_q37_sorted_rmdup_uniq_RG.vcf \
/data/stonelab/HLP/VCF_files/3077/3077_trimmed_merged_q37_sort_rmdup_uniq_RG.vcf \
/data/stonelab/HLP/VCF_files/3208/3208_trimmed_mapped_q37_sorted_rmdup_uniq_RG.vcf \
/data/stonelab/HLP/VCF_files/511/511_trimmed_merged_q37_sort_rmdup_uniq_RG.vcf \
/data/stonelab/HLP/VCF_files/515/515_trimmed_merged_q37_sort_rmdup_uniq_RG.vcf \
/data/stonelab/HLP/VCF_files/516/516_trimmed_merged_q37_sort_rmdup_uniq_RG.vcf \
/data/stonelab/HLP/VCF_files/517/517_trimmed_merged_q37_sort_rmdup_uniq_RG.vcf \
/data/stonelab/HLP/VCF_files/518/518_trimmed_merged_q37_sort_rmdup_uniq_RG.vcf \
/data/stonelab/HLP/VCF_files/519/519_trimmed_merged_q37_sort_rmdup_uniq_RG.vcf \
/data/stonelab/HLP/VCF_files/520/520_trimmed_merged_q37_sort_rmdup_uniq_RG.vcf \
/data/stonelab/HLP/VCF_files/523A1/523A1_trimmed_merged_q37_sort_rmdup_uniq_RG.vcf \
/data/stonelab/HLP/VCF_files/536/536_trimmed_merged_q37_sort_rmdup_uniq_RG.vcf \
/data/stonelab/HLP/VCF_files/97016/97016_trimmed_mapped_q37_sorted_rmdup_uniq_RG.vcf \
/data/stonelab/HLP/Multi_VCF/outgroup_lepromatosis/outgroup.vcf 
#to be treated as an outgroup (only call variants if they are also present in the sample genomes), 
#the outgroup .vcf must be in a directory named "outgroup_[sample_name]"