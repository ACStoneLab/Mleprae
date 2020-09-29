#!/bin/bash

#Variants were filtered and aligned using MultiVCFAnalyzer, which names samples based on 
#the name of the directory they're in. This script creates a directory with the sample 
#name and copies the .vcf and vcf.idx files to this new directory. 
#Note the directory structure follows that of the previous scripts.

path=/where/input/data/are

for a in $path/Sample*
 do
  
  ID=$(basename $a)
  NAME=$(echo $ID |cut -d "_" -f2)
  
  mkdir $path/"${ID}"/Multi_VCF
  mkdir $path/"${ID}"/Multi_VCF/"${NAME}"
  
  cp $path/modern_working/"${ID}"/5.gatk_unified_genotyper/"$NAME"_*.vcf $path/Multi_VCF/"${NAME}"
  cp $path/modern_working/"${ID}"/5.gatk_unified_genotyper/"$NAME"_*.vcf.idx $path/Multi_VCF/"${NAME}"
  
  done