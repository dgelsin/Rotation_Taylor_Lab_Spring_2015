# bash script for running INSIGHT-EM on input file representing a given collection of elements
# infers parameters of interest, and LRT statistics for the three tests: rho>0, eta>0, and gamma>0.
#
# usage: bash scripts/getSelectionEstimates.sh <elementSummaryFile> <freqThreshold> <beta1> <beta2> <beta3> <logFile>
#
# * elementSummaryFile   - summary file for elements to analyze (generated by processElementBed/processElements.sh script)
# * freqThreshold        - frequency threshold (in %) used in distinguishing 'L' and 'H' sites
# * beta1,2,3            - estimates for beta1, beta2, and beta3 obtained by analyzing flanking sites
# * logFile              - name of file for logging INSIGHT-EM output (appends existing file)
#
# ==> progress is logged onto stdout.
#     successful completion is indicated by last log line being 'Done.'

set -e
set -o pipefail

# path to executable
INSIGHT_EM=./bin/INSIGHT-EM-v1.1
# global options
INSIGHT_EM_OPTIONS="-v -i 100000"

elementSummaryFile=$1
freqThreshold=$2
beta1=$3
beta2=$4
beta3=$5
logFile=$6
postFile=$7

scriptDir=./scripts

echo '=================================================================='
echo '  Estimating selection parameters for freq threshold '$freqThreshold' on file '$elementSummaryFile
echo '=================================================================='

elementInputFile=${elementSummaryFile}.f${freqThreshold}

##############  -1- APPLY FREQ THRESHOLD TO INPUT FILE AND ADD BETA LINE ###############################

echo '1.  applying frequency threshold to input file and adding beta line'
   bash $scriptDir/labelLHsites.sh   $elementSummaryFile   $elementInputFile   $freqThreshold
   echo -e 'beta\t'$beta1'\t'$beta2'\t'$beta3   >>   $elementInputFile

##############  -2- ESTIMATE SELECTION PARAMETERS ###############################

echo '2.  estimating selection parameters'

   tmpLogFile=$logFile.estimate
   cmd=$INSIGHT_EM" "$INSIGHT_EM_OPTIONS
   if [ ""$postFile != "--noPost" -a ""$postFile != "" ]
   then
      cmd=$cmd" -p "$postFile
   fi
   cmd=$cmd" "$elementInputFile
   echo '    Invoking  '$cmd
   set +e
   $cmd   &>   $tmpLogFile
   if [ $? -ne 0 ]
   then
      echo '    ==> Error when executing INSIGHT-EM to estimate selection parameters.'
      echo '    ==> Command used: '$cmd
      echo '    ==> See output log file '$tmpLogFile
      exit 1
   fi
   set -e
   bash $scriptDir/getEMstatus.sh   $tmpLogFile   $cmd
   echo -e "thres-"$freqThreshold"-main\tCommandline: "$cmd >> $logFile   
   cat $tmpLogFile  |  awk -v OFS="\t" -v thres=$freqThreshold '{print "thres-"thres"-main",$0}'   >>   $logFile
   rm $tmpLogFile


##############  -3- LRT FOR RHO>0 ###############################

echo '3.  restricting rho=0 for LRT for selection'

   tmpLogFile=$logFile.rho0
   cmd=$INSIGHT_EM" "$INSIGHT_EM_OPTIONS" -fr -r 0.0 -c "$elementInputFile
   echo '    Invoking  '$cmd
   set +e
   $cmd   &>   $tmpLogFile
   if [ $? -ne 0 ]
   then
      echo '    ==> Error when executing INSIGHT-EM for LRT of rho>0.'
      echo '    ==> Command used: '$cmd
      echo '    ==> See output log file '$tmpLogFile
      exit 1
   fi
   set -e
   bash $scriptDir/getEMstatus.sh   $tmpLogFile   $cmd
   echo -e "thres-"$freqThreshold"-rho0\tCommandline: "$cmd >> $logFile   
   cat $tmpLogFile  |  awk -v OFS="\t" -v thres=$freqThreshold '{print "thres-"thres"-rho0",$0}'   >>   $logFile
   rm $tmpLogFile

##############  -4- LRT FOR ETA>0 ###############################

echo '4.  restricting eta=0 for LRT for positive selection'

   tmpLogFile=$logFile.eta0
   cmd=$INSIGHT_EM" "$INSIGHT_EM_OPTIONS" -fe -e 0.0 -c "$elementInputFile
   echo '    Invoking  '$cmd
   set +e
   $cmd   &>   $tmpLogFile
   if [ $? -ne 0 ]
   then
      echo '    ==> Error when executing INSIGHT-EM for LRT of eta>0.'
      echo '    ==> Command used: '$cmd
      echo '    ==> See output log file '$tmpLogFile
      exit 1
   fi
   set -e
   bash $scriptDir/getEMstatus.sh   $tmpLogFile   $cmd
   echo -e "thres-"$freqThreshold"-eta0\tCommandline: "$cmd >> $logFile   
   cat $tmpLogFile  |  awk -v OFS="\t" -v thres=$freqThreshold '{print "thres-"thres"-eta0",$0}'   >>   $logFile
   rm $tmpLogFile

##############  -5- LRT FOR GAMMA>0 ###############################

echo '5.  restricting gamma=0 for LRT for selection'

   tmpLogFile=$logFile.gam0
   cmd=$INSIGHT_EM" "$INSIGHT_EM_OPTIONS" -fg -g 0.0 -c "$elementInputFile
   echo '    Invoking  '$cmd
   set +e
   $cmd   &>   $tmpLogFile
   if [ $? -ne 0 ]
   then
      echo '    ==> Error when executing INSIGHT-EM for LRT of gamma>0.'
      echo '    ==> Command used: '$cmd
      echo '    ==> See output log file '$tmpLogFile
      exit 1
   fi
   set -e
   bash $scriptDir/getEMstatus.sh   $tmpLogFile   $cmd
   echo -e "thres-"$freqThreshold"-gam0\tCommandline: "$cmd >> $logFile   
   cat $tmpLogFile  |  awk -v OFS="\t" -v thres=$freqThreshold '{print "thres-"thres"-gam0",$0}'   >>   $logFile
   rm $tmpLogFile
   
##############  -5- CLEAN UP ###############################

   rm $elementInputFile

echo 'Done.'
