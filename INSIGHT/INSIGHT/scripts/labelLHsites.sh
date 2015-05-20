# bash script for applying frequency threshold to a site summary file
# lebeling all polymorphic sites ('P') as 'L' or 'H' based on the minor allele frequency
#
# usage: bash scipts/labelLHsites.sh <summaryFile> <outFile> <freqThresold>
#
# * summaryFile   - input summary file
# * outFile       - name of file for output
# * freqThreshold - frequency threshold (in %)

set -e
set -o pipefail

summaryFile=$1
outFile=$2
freqThreshold=$3

thres=`echo $freqThreshold/100 | bc -l`

cat $summaryFile | \
   awk -v thres=$thres \
   '{ \
      if($1=="site" && $3=="P") { \
         freq=$7/($6+$7); \
         if(freq<thres) {$3="L"} \
         else {$3="H"}; \
      } \
      print;
   }' \
   > $outFile
