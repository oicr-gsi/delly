#/bin/bash
cd $1
 
find . -name "*.zip" -exec unzip {} \; >/dev/null
find . -type f -not -path "./*.zip" -exec md5sum {} + 
 

