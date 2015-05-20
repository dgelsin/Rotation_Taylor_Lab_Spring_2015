# bash script for processing results of INSIGHT-EM
#
# usage: bash scripts/processEMresults.sh <logFile> <dataID> <resultsFile>
#
# * logFile           - log file where INSIGHT-EM output is logged
#                       if logFile == --header, then only header is produced as output and not processed version of output
# * dataID            - id of data set to use in output line
# * resultsFile       - name of file where result summary line is to be output
# * --header           - when optional flag is given, then only header is produced in output file

set -e
set -o pipefail

logFile=$1
dataID=$2
resultsFile=$3


if [ $1"" == "--header" ]
then
   echo -e "dataID\tthres\trho\trho_stderr\tE[A]\tE[A]_stderr\tE[W]\tE[W]_stderr\talpha\talpha_stderr\ttau\ttau_stderr\teta\teta_stderr\tgamma\tgamma_stderr\tlnLd\tLRT[rho>0]\tLRT[eta>0]\tLRT[gamma>0]\tem_status" \
      > $resultsFile
else 
   cat $logFile | \
      sed -e 's/thres-/f=/' | \
      sed -e 's/-beta/%\tbeta/' | \
      sed -e 's/-main/%\tmain/' | \
      sed -e 's/-rho0/%\trho0/' | \
      sed -e 's/-eta0/%\teta0/' | \
      sed -e 's/-gam0/%\tgam0/' | \
      awk -v OFS="\t" -v dataID=$dataID\
      'BEGIN{  thres = "------"; } \
      { \
         if($1!=thres) { \
            if(thres != "------") { \
               print dataID, thres, rho, rho_stderr, Ea, Ea_stderr, Ew, Ew_stderr, alpha, alpha_stderr, tau, tau_stderr, eta, eta_stderr, gamma, gamma_stderr, lnLd, LRTrho, LRTeta, LRTgamma, em_status; \
            } \
            thres=$1; rho="na"; rho_stderr="na"; eta="na"; eta_stderr="na"; gamma="na"; gamma_stderr="na"; \
		      Ea="na"; Ea_stderr="na"; Ew="na"; Ew_stderr="na"; alpha="na"; alpha_stderr="na"; tau="na"; tau_stderr="na"; \
            lnLd="na"; em_status="na"; \
		      LRTrho="na";LRTeta="na";LRTgamma="na";\
         } \
         if($2=="main") { \
            if($3=="Estimates:") {  rho=$4; eta=$5; gamma=$6; Ea=$7; Ew=$8; alpha=$9; tau=$10;  } \
            else if($3=="StndrdErr:") {  rho_stderr=$4; eta_stderr=$5; gamma_stderr=$6; \
               Ea_stderr=$7; Ew_stderr=$8; alpha_stderr=$9; tau_stderr=$10;   } \
            else if($3=="EM" && $4=="status:") {  lnLd=$6; em_status="main-"$8; } \
         } else if($2=="rho0" && $3=="EM" && $4=="status:") { \
            LRTrho = 2*(lnLd - $6); em_status=em_status"-rho0-"$8; \
         } else if($2=="eta0" && $3=="EM" && $4=="status:") { \
            LRTeta = 2*(lnLd - $6); em_status=em_status"-eta0-"$8; \
         } else if($2=="gam0" && $3=="EM" && $4=="status:") { \
            LRTgamma = 2*(lnLd - $6); em_status=em_status"-gam0-"$8; \
         } \
      } \
      END{ print dataID, thres, rho, rho_stderr, Ea, Ea_stderr, Ew, Ew_stderr, alpha, alpha_stderr, tau, tau_stderr, eta, eta_stderr, gamma, gamma_stderr, lnLd, LRTrho, LRTeta, LRTgamma, em_status; }' \
      >> $resultsFile
fi


