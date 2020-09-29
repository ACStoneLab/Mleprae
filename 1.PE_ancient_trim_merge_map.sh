#!/bin/bash

#Modified from scripts written by Dr. Tanvi Honap, University of Oklahoma and Dr. Maria Nieves-Colon, Arizona State University

#Raw data should be stored in uniquely named folders with the prefix "Sample" such as Sample_1, Sample_2, Sample_3

#This script performs the following on ancient PE sequencing data 
#1.fastqc on raw data
#2.adapter removal: trimming and merging
#3.fastqc on trimmed and merged data
#4.bwa mapping and samtools filtering
#4.5 qualimap (mapping statistics)

module load fastqc/0.11.7
module load bwa/0.7.17
module load samtools/1.9
eval "$(conda shell.bash hook)"  #necessary for calling conda activate env on server
conda activate working #activates conda environment for use
qualimap=/path/to/qualimap

set -ex 
#prevents the next commands from being executed after a command has failed
#also keeps track of which command is being executed and prints that to the command line

path=/where/raw/data/are #input path for where raw data are
outputdir=/where/output/directory/is #output directory; make sure input does not have a final forward slash
ref=/where/ref/fasta/is.fasta

#indexing reference only needs to be done once
bwa index $ref
samtools faidx $ref

echo -e "Sample \t Total Reads \t Analysis-ready reads (trimmed and merged) \t Proportion of reads kept after trimming and merging \t Mapped reads \t Proportion of reads mapped (Mapped reads / Analysis-ready reads) \t Q37 mapped reads \t Unique Q37 mapped reads \t Average length of mapped reads \t Mean coverage \t Std dev of mean coverage \t % reference covered >= 1X \t % reference covered >= 5X \t Cluster factor (Mapped reads/ Unique Q37 mapped reads) \t % endogenous DNA (Unique Q37 mapped reads / Total reads)" > $outputdir/ancient_PE_leprae_pipeline_summary.txt

for a in $path/Sample_*
 do
  cd $a
  
  ID=$(basename $a)  #sets variable ID equal to the variables in the for loop, which are directory names
  NAME=$(echo $ID |cut -d "_" -f2)
  
  mkdir $outputdir/$ID
  
  echo "Running fastqc on sample "${NAME}"!" 
  
  mkdir $outputdir/"${ID}"/1.raw_fastqc 
  fastqc $path/"${ID}"/*R1.fastq.gz -o $outputdir/"${ID}"/1.raw_fastqc 
  fastqc $path/"${ID}"/*R2.fastq.gz -o $outputdir/"${ID}"/1.raw_fastqc 
  
  echo "Running adapter removal on sample "${NAME}"!"  #prints to screen which sample is being processed
  
  mkdir $outputdir/$ID/2.adapter_removal
  
  AdapterRemoval --file1 $path/"${ID}"/*R1.fastq.gz --file2 $path/"${ID}"/*R2.fastq.gz --trimns --trimqualities --minquality 20 --collapse --minlength 30 --settings $outputdir/"${ID}"/2.adapter_removal/"${NAME}"_settingz.txt --output1 $outputdir/"${ID}"/2.adapter_removal/"${NAME}"_trimmed_R1.fastq --output2 $outputdir/"${ID}"/2.adapter_removal/"${NAME}"_trimmed_R2.fastq --outputcollapsed $outputdir/"${ID}"/2.adapter_removal/"${NAME}"_trimmed_merged.fastq  --singleton $outputdir/"${ID}"/2.adapter_removal/"${NAME}"_singleton_truncated.fastq --outputcollapsedtruncated $outputdir/"${ID}"/2.adapter_removal/"${NAME}"_outputcollapsedtruncated.fastq --discarded $outputdir/"${ID}"/2.adapter_removal/"${NAME}"_discarded.fastq
  
#runs AdapterRemoval on the specified input files; parameters are defined below. 
  
#--trimqualities #trims bases with quality scores below specified quality score
#--trimns #trims ambiguous calls (Ns)
#--minquality 20 #specified quality score
#--collapse #asks AdapterRemoval to merge reads
#--minlength 30 #sets a min length for kept fragments

  echo "Running fastqc on trimmed and merged reads "${ID}"!" 
  
  mkdir $outputdir/$ID/3.post_AR_fastqc
  
  fastqc $outputdir/"${ID}"/2.adapter_removal/"${NAME}"_trimmed_merged.fastq -o $outputdir/"${ID}"/3.post_AR_fastqc

#adapter removal stats
  total_reads=`grep "Total number of read pairs" $outputdir/"${ID}"/2.adapter_removal/"${NAME}"_settingz.txt | awk '{print $NF}'` #this will not work for SE samples, but there are only 2 so I will calculate manually
  analysis_ready=`grep "Number of full-length collapsed pairs" $outputdir/Sample_"${NAME}"/2.adapter_removal/"${NAME}"_settingz.txt | awk '{print $NF}'`
  pc_merged=`echo "scale=4;$analysis_ready/$total_reads" | bc`
  
## bwa aln finds the SA coordinates of the input reads
##  -l 1000 disables seed for usage with aDNA reads; do not use if you are working with modern DNA data
##  -n denotes edit distance and allows more or less mismatches; 0.04 (default) or reduce to 0.01 (less stringent) or increase to 0.1 (to make mapping more strict)
##  -t denotes number of threads used - can be changed as desired
## bwa samse generates alignments in the SAM format given single-end reads

  mkdir $outputdir/$ID/4.bwa
  
## Here we are using only the trimmed and merged reads as input. Modify the following if you need to use unmerged R1 and R2 as well.
  bwa aln -l 1000 -n 0.1 -t 8 $ref $outputdir/Sample_"${NAME}"/2.adapter_removal/"${NAME}"_trimmed_merged.fastq > $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged.sai

  bwa samse $ref $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged.sai $outputdir/Sample_"${NAME}"/2.adapter_removal/"${NAME}"_trimmed_merged.fastq > $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged_bwa_all.sam
  
  samtools view -bSh $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged_bwa_all.sam > $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged_bwa_all.bam
  
## Filtering, sorting, removing duplicates using SAMtools
## samtools view -bSh displays SAM file as a BAM file (b); input is SAM (S),and includes the header (h)

## Create a report summary file for each sample. Counts from the report will be combined into overall summary report.

  echo " " >> $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_MappingReport.txt
  echo "${NAME}" >> $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_MappingReport.txt
  echo "---------------------------------------" >> $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_MappingReport.txt
  echo "Analysis-ready reads" >> $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_MappingReport.txt
  
## Use samtools view to count alignments and print the total number to the report summary file
  samtools view -c $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged_bwa_all.bam >> $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_MappingReport.txt
  
## Filter out unmappped and low quality reads
## samtools view -bh -F4 -q displays the previous output as a BAM file (b), and
## includes the header (h), but skips alignments with MAPQ smaller than than 37
## (-q 37), and alignments containing the 4 flag (0x4 segment unmapped)
  
  samtools view -bh -F4 $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged_bwa_all.bam > $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged_bwa_mapped.bam

  samtools view -bh -q 37 $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged_bwa_mapped.bam > $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged_bwa_mapped_q37.bam
  
## Sort alignments by leftmost coordinates

  samtools sort -o $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged_bwa_mapped_q37_sorted.bam $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged_bwa_mapped_q37.bam
  
## Print alignment statistics after removing low quality reads

  echo "Mapped reads" >> $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_MappingReport.txt
  samtools view -c $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged_bwa_mapped.bam >> $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_MappingReport.txt
  echo "Q37 mapped reads" >> $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_MappingReport.txt
  samtools view -c $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged_bwa_mapped_q37_sorted.bam >> $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_MappingReport.txt

## Remove PCR duplicates
## samtools markdup marks potential PCR duplicates, -r removes them
## if multiple read pairs have identical external coordinates, only retain the pair with highest mapping quality

  samtools markdup -r $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged_bwa_mapped_q37_sorted.bam $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged_bwa_mapped_q37_sorted_rmdup.bam

## Print results after markdup

  echo "After removing duplicates" >> $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_MappingReport.txt
  samtools view -c $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged_bwa_mapped_q37_sorted_rmdup.bam >> $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_MappingReport.txt

## Remove reads with multiple mappings
## grep -v means invert the match, it returns all non matching lines (all the reads that are not multiple mappings)
## [XT:A:U flag in the SAM file denotes unique read, and XT:A:R and XA:Z denote multiple mappings for that read]
## When working with paired end reads, it may also be useful to consider filtering out reads with the flag XT:A:M
## (one-mate recovered) which means that one of the pairs is uniquely mapped and the other isn't]
## Use awk to scan input file for lines that match the pattern: {if($0~/X1:i:0/||$0~/^@/)print $0}
## X1= Number of suboptimal hits found by BWA, if X1=0 then keep that read in file

  samtools view -h $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged_bwa_mapped_q37_sorted_rmdup.bam |grep -v 'XT:A:R'|grep -v 'XA:Z' |grep -v 'XT:A:M' |awk '{if($0~/X1:i:0/||$0~/^@/)print $0}' | samtools view -bS - > $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged_bwa_mapped_q37_sorted_rmdup_uniq.bam

## The Sample_trimmed_merged_bwa_mapped_q37_sorted_rmdup_uniq.bam file is the final BAM file that should be used for further analyses
## If working with modern DNA, use this file for further analyses

## Print results after removing reads with multiple mappings

  echo "Unique reads" >> $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_MappingReport.txt
  samtools view -c $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged_bwa_mapped_q37_sorted_rmdup_uniq.bam >> $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_MappingReport.txt

## Print average length in report using awk

  echo "Average length of mapped reads" >> $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_MappingReport.txt
  samtools view $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged_bwa_mapped_q37_sorted_rmdup_uniq.bam |awk '{SUM+=length($10);DIV++}END{print SUM/DIV}' >> $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_MappingReport.txt

## Index alignment and print additional statistics

  samtools index $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged_bwa_mapped_q37_sorted_rmdup_uniq.bam
  samtools flagstat $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged_bwa_mapped_q37_sorted_rmdup_uniq.bam >> $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_flagstat.txt

## Calculate statistics
  
## optional removal of temporary files
  rm $path/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged.sai
  rm $path/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged_bwa_all.sam
  rm $path/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged_bwa_all.bam
  rm $path/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged_bwa_mapped.bam
  rm $path/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged_bwa_mapped_q37.bam
  rm $path/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged_bwa_mapped_q37_sorted.bam
  rm $path/Sample_"${NAME}"/4.bwa/"${NAME}"_trimmed_merged_bwa_mapped_q37_sorted_rmdup.bam

## Calculate statistics
  mapped=`grep -A 1 "Mapped reads" $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_MappingReport.txt | grep -v "Mapped reads"`
  Q37=`grep -A 1 "Q37 mapped reads" $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_MappingReport.txt | grep -v "Q37 mapped reads"`
  unique=`grep -A 1 "Unique reads" $outputdir/Sample_"${NAME}"/4.bwa/"${NAME}"_MappingReport.txt | grep -v "Unique reads"`
  pc_mapped=`echo "scale=4;$mapped/$analysis_ready" | bc`
  length=`grep -A 1 "Average length of mapped reads" $outputdir/"${ID}"/4.bwa/"${NAME}"_MappingReport.txt | grep -v "Average length of mapped reads"`
  cf=`echo "scale=4;$mapped/$unique" | bc`
  endo=`echo "scale=4;$mapped/$analysis_ready" | bc`
  
  mkdir $outputdir/"${ID}"/4.5.Qualimap
  
  $qualimap bamqc -bam $outputdir/"${ID}"/4.bwa/"${NAME}"_trimmed_merged_bwa_mapped_q37_sorted_rmdup_uniq.bam -outdir $outputdir/"${ID}"/4.5.Qualimap -outformat pdf \

  mean_cov=`grep "mean coverageData" $outputdir/"${ID}"/4.5.Qualimap/genome_results.txt | cut -d "=" -f2`
  sd_cov=`grep "std coverageData" $outputdir/"${ID}"/4.5.Qualimap/genome_results.txt | cut -d "=" -f2`
  cov1=`grep "reference with a coverageData >= 1X" $outputdir/"${ID}"/4.5.Qualimap/genome_results.txt | awk '{print $4}'`
  cov5=`grep "reference with a coverageData >= 5X" $outputdir/"${ID}"/4.5.Qualimap/genome_results.txt | awk '{print $4}'`
  
  echo -e "${NAME} \t $total_reads \t $analysis_ready \t $pc_merged \t $mapped \t $pc_mapped \t $Q37 \t $unique \t $length \t $mean_cov \t $sd_cov \t $cov1 \t $cov5 \t $cf \t $endo " >> $outputdir/ancient_PE_leprae_pipeline_summary.txt

 
done



