#!/bin/bash

#This script adds read groups to the analysis-ready bam, and performs variant calling 
#using GATK UnifiedGenotyper

module load picard/2.18.3
module load samtools/1.9

set -ex 

#prevents the next commands from being executed after a command has failed
#also keeps track of which command is being executed and prints that to the command lin

path=/where/input/analysis/ready/bams/live
ref=/path/to/reference.fasta

#if you need to prepare .dict of reference
#picard CreateSequenceDictionary R=/path/to/reference.fasta O=/path/to/output/.fasta.dict

for a in $path/Sample*
 do
  cd $a
  
  ID=$(basename $a)
  NAME=$(echo $ID |cut -d "_" -f2)
  
  mkdir $path/$ID/5.gatk_unified_genotyper
  
  #bam files must have read groups for ANY gatk program
  picard AddOrReplaceReadGroups I=$path/"${ID}"/4.bwa/"${NAME}"_trimmed_mapped_q37_sorted_rmdup_uniq_paired.bam O=$path/"${ID}"/5.gatk_unified_genotyper/"$NAME"_trimmed_mapped_q37_sorted_rmdup_uniq_RG.bam RGID="${NAME}" RGLB=Lib RGPL=ILLUMINA RGPU=unit1 RGSM="${NAME}"

  #must re-index bam file with read groups
  samtools index $path/"${ID}"/5.gatk_unified_genotyper/"$NAME"_trimmed_mapped_q37_sorted_rmdup_uniq_RG.bam

  java -jar java -jar /path/to/GATK.jar \
    -T UnifiedGenotyper \
    -R $ref \
    -I $path/"${ID}"/5.gatk_unified_genotyper/"$NAME"_trimmed_mapped_q37_sorted_rmdup_uniq_RG.bam \
    -o $path/"${ID}"/5.gatk_unified_genotyper/"$NAME"_trimmed_mapped_q37_sorted_rmdup_uniq_RG.vcf \
    --output_mode EMIT_ALL_SITES

  done