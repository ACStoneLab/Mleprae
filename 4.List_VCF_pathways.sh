#!/bin/bash

##This script generates a list of pathways for the VCF files to put in the
##MultiVCFAnalyzer config file
#This script may need to be modified if all of the VCF files do not have the same suffix

set -ex 
#prevents the next commands from being executed after a command has failed
#also keeps track of which command is being executed and prints that to the command line

path=/path/to/intput/data 
outputdir=/path/to/output/directory

echo -e "List of VCF pathways" > $outputdir/VCF_pathway_list.txt

for a in $path/Sample*
 do
  
  ID=$(basename $a)  #sets variable ID equal to the variables in the for loop, which are directory names
  NAME=$(echo $ID |cut -d "_" -f2)
  PWD="${path}"/"${ID}"/Multi_VCF/"${NAME}"/"${NAME}"_trimmed_mapped_q37_sorted_rmdup_uniq_RG.vcf" \\"
  
  echo $PWD >> $path2/VCF_pathway_list.txt
  
  done
  