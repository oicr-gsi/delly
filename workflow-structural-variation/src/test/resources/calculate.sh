#/bin/bash
cd $1
 
find . -name "*.vcf" -exec wc -l {} +
 

