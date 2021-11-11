#!/bin/bash
cd $1

# - .vcf.gz files have no stochastic content except a header lines with time information
#   Therefore:
# - Check md5sums for all types of files, sort

echo ".vcf files:"
for v in *.vcf.gz;do zcat $v | grep -v Date | grep -v ^##fileDate | md5sum;done | sort

