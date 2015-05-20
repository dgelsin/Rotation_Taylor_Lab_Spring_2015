# script for adding pseuodo sites to EM input files
# add a single pseudo block with theta=0.001 and lambda=0.001 and a single H polymorphic site
#
# usage: bash scripts/addPseudocounts.sh  <emInputFile>  [<flankInputFile>]
#
# adds pseudo sites to both EM input files (one for elements and one for flanks, if given) 

set -e
set -o pipefail

elementEMInputFile=$1
flankEMInputFile=$2

theta=0.001
lambda=0.001

blockLine="block\tchrPseudo:0-500\ttheta\t"$theta"\tlambda\t"$lambda
hSiteLine="site\tchrP:0\tP\t1.0\t0.0\t50\t50"
echo -e $blockLine >> $elementEMInputFile
echo -e $hSiteLine >> $elementEMInputFile

if [ ""$flankEMInputFile != "" ]
then
   echo -e $blockLine >> $flankEMInputFile
   echo -e $hSiteLine >> $flankEMInputFile
fi

