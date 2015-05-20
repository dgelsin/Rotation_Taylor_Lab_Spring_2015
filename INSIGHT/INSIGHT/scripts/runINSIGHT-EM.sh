# main bash script for running INSIGHT-EM on processed input files
#
# usage: bash scripts/runINSIGHT-EM.sh <dataID>  <emInputElements>  <emInputFlanks>  <outDir>  <freqThreshold>  [--pseudo]
#
# * dataID            - unique ID of data set
# * emInputElements   - EM input file for elements
# * emInputFlanks     - EM input file for flanking polymorphic sites.
# *                     can alternatively be a file with a single beta line (beta b1 b2 b3).
#                       if the first token in the file is 'beta', the file is assumed to be a beta file
# * outDir       - directory where analysis is summarized
# * freqThreshold     - frequency threshold (in %) used in distinguishing 'L' and 'H' sites
# * --pseudo          - optinal flag. If provided, then generates pseudo sites for elements and flanks
#
# ==> progress is logged (with timestamps) onto stdout.
#     successful completion is indicated by last log line being 'Done. time <time>'

set -e
set -o pipefail

dataID=$1
emInputElements=$2
emInputFlanks=$3
outDir=$4
freqThreshold=$5
usePseudocounts=$6

emInputTemp=$outDir/$dataID.tmp
logFile=$outDir/$dataID.ins.log
resultsFile=$outDir/$dataID.ins.results.txt
postCountFile_sites=$outDir/$dataID.post.sites
givenBetas="FALSE"
betaFile="NOFILE"
token=`head -n1 $emInputFlanks | awk '{print $1}'`
if [ $token == "beta" ]
then
   givenBetas="TRUE"
   betaFile=$emInputFlanks
   emInputFlanks="--noFlanks" 
fi

scriptDir=./scripts

echo '=================================================================='
echo '  Running INSIGHT-EM on data set '$dataID
echo '  Frequency thresholds: '$freqThreshold
echo '=================================================================='

##############  -A- CREATE TEMPORARY INPUT FILES ###############################

   cp   $emInputElements  $emInputTemp.ins
   if [ $givenBetas != "TRUE" ]
   then
      cp   $emInputFlanks    $emInputTemp.flankPoly.forBetas.ins
   fi

##############  -B- APPLY PSEUDO COUNTS ###############################

   if [ ""$usePseudocounts == "--pseudo" ]
   then
      echo '    adding pseudo-count sites'
      if [ $givenBetas == "TRUE" ]
      then
         bash $scriptDir/addPseudocounts.sh   $emInputTemp.ins  --noFlanks
      else
         bash $scriptDir/addPseudocounts.sh   $emInputTemp.ins  $emInputTemp.flankPoly.forBetas.ins
      fi
   fi


##############  -C- GET BETA ESTIMATES  ###############################

echo ''
   if [ $givenBetas == "TRUE" ]
   then
      betas=` \
         cat $betaFile | \
         head -n1 | \
         awk '{print $(NF-2),$(NF-1),$(NF);}'   `
   else
      bash $scriptDir/getBetaEstimates.sh   $emInputTemp.flankPoly.forBetas.ins   $freqThreshold   $logFile
      betas=` \
         cat $logFile | \
         grep "thres-"$freqThreshold"-beta" | \
         grep "Beta Line:" | \
         awk '{print $(NF-2),$(NF-1),$(NF);}'   `
   fi

##############  -B- GET SELECTION ESTIMATES  ###############################

echo ''
   bash $scriptDir/getSelectionEstimates.sh \
      $emInputTemp.ins   $freqThreshold   $betas   $logFile   $postCountFile_sites

##############  -D- PROCESS RESULTS  ###############################

echo ''
echo '=================================================================='
echo '  Processing results '
echo '=================================================================='

   bash $scriptDir/processEMresults.sh   $logFile   $dataID   $resultsFile


##############  -E- CLEAN UP  ###############################

   rm -f $emInputTemp.*

##############  -F- DONE  ###############################

echo 'Done.'

