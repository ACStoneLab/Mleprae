#!/bin/bash

#This script rescales base qualities based on expected damage patterns for ancient DNA,
#adds read groups to the rescaled bam, and performs variant calling using GATK
#UnifiedGenotyper

module load picard/2.18.3
module load samtools/1.9
module load mapdamage/2.0.9

set -ex 
#prevents the next commands from being executed after a command has failed
#also keeps track of which command is being executed and prints that to the command line

path=/where/input/data/are #input path for where raw data are
outputdir=/where/ouput/data/go #input path for output directory; make sure input does not have a final forward slash
ref=/path/to/reference.fasta

#if you need to prepare .dict of reference
#picard CreateSequenceDictionary R=/path/to/reference.fasta O=/path/to/output/.fasta.dict

for a in $path/Sample*
 do
  cd $a
  
  ID=$(basename $a)
  NAME=$(echo $ID |cut -d "_" -f2)
  
  mkdir $outputdir/"${ID}"/mapdamage
  
  mapDamage -i $path/"${ID}"/4.bwa/"${NAME}"_trimmed_merged_bwa_mapped_q37_sorted_rmdup_uniq.bam -r $ref -d $outputdir/"${ID}"/mapdamage --rescale
  samtools index $outputdir/"${ID}"/mapdamage/"{NAME}"_trimmed_merged_bwa_mapped_q37_sorted_rmdup_uniq.rescaled.bam

  mkdir 5.gatk_unified_genotyper
  
  #bam files must have read groups for ANY gatk program
  picard AddOrReplaceReadGroups I=$path/"${ID}"/mapdamage/"${NAME}"_trimmed_merged_bwa_mapped_q37_sorted_rmdup_uniq.rescaled.bam O=$path/"${ID}"/5.gatk_unified_genotyper/"$NAME"_trimmed_merged_q37_sort_rmdup_uniq_rescaled_RG.bam RGID="${NAME}" RGLB=Lib RGPL=ILLUMINA RGPU=unit1 RGSM="${NAME}"

  #must re-index bam file with read groups
  samtools index $path/"${ID}"/5.gatk_unified_genotyper/"$NAME"_trimmed_merged_q37_sort_rmdup_uniq_rescaled_RG.bam

  java -jar /path/to/GATK.jar \
    -T UnifiedGenotyper \
    -R $ref \
    -I $path/"${ID}"/5.gatk_unified_genotyper/"$NAME"_trimmed_merged_q37_sort_rmdup_uniq_rescaled_RG.bam \
    -o $path/"${ID}"/5.gatk_unified_genotyper/"$NAME"_trimmed_merged_q37_sort_rmdup_uniq_rescaled_RG.vcf \
    --output_mode EMIT_ALL_SITES

  done