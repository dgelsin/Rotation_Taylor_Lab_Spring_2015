# bash script for getting exit status of INSIGHT-EM and producing the appropriate error message
#
# usage: bash scripts/getEMstatus.sh <logFile> <cmd>
#
# * logFile - log file to process
# * cmd     - commandline used

set -e
set -o pipefail

logFile=$1
shift
cmd=$@

emStatus="status_"`cat $logFile | grep "EM status:" | awk '{print $NF}'`
if [ $emStatus == "status_overshoot" ]
then
   echo '    ==> Warning: INSIGHT-EM exitted with status '$emStatus
   echo '    ==> Command used: '$cmd
   echo '    ==> See output log file '$logFile'.sav'
   echo '    ==> Consider relaxing convergence criterion by increasing the -d input option'
   cp $logFile $logFile.sav
elif [ $emStatus != "status_converged" ]
then
   echo '    ==> Error: INSIGHT-EM exitted with status '$emStatus
   echo '    ==> Command used: '$cmd
   echo '    ==> See output log file '$logFile
   exit 1
else
   echo '    INSIGHT-EM converged.'
fi

