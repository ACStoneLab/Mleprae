#!/bin/bash

##This script reads a txt file with two columns (no header) 
##Column one is Sample name and column two is the NCBI SRR# or ENA ERR#
##It fetches raw single-end data using the fastq-dump tool, gzips, and renames the
##fastq files using the name specified in the input file

module load sratoolkit/2.8.2-1

set -ex

filename=/Path/to/file.txt/with/names/and/SRRs
outputdir=/Path/to/output/directory

while read NAME SRR; do
 mkdir $outputdir/Sample_$NAME
 fastq-dump --outdir $outputdir/Sample_$NAME --gzip "${SRR}"
 mv $outputdir/Sample_$NAME/*.gz $outputdir/Sample_$NAME/$NAME.fastq.gz
 
 done < $filename
