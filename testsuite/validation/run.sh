#!/bin/bash 
#
# T-Coffee validation script. It runs a T-Coffee alignment for each entry in the Balibase RV11 dataset. 
# The result is compared with the reference alignment using the bali_score result as metric. 
# 
# If the average result is above or equals to the threshold value (default: 0.490) the validation PASSED, 
# NOT PASSED otherwise. 
# If the average result is below the threshold value by a 10% the validation is marked as FAILED and 
# the script will exit with the error code 1.
#  
# Usage: ./run.sh [dataset-against-to-test]
#
#
# Environment variables that can be used to override the default values 
# - THRESHOLD: define the validation threshold score (default: 0.490)
# - SCRATCH: the path where validation will be executed (default: ./scratch)
# - FAILED: the path where log and result files for alignment with a score value below the threshold are stored for further analysis (default: ./failed)
# - OSNAME: the system type (linux|macosx)
# - OSARCH: the system architecture (i386|x86_64) 
#
#
set -e
set -o errexit
#set -o nounset
#set -u
#set -x

##
## Defines the OSNAME variable
## 
if [ -z $OSNAME ]; then
  if [[ $OSTYPE == darwin* ]] ; then 
    OSNAME='macosx'
  fi
  if [[ $OSTYPE == linux* ]] ; then
    OSNAME='linux'
  fi 
fi
if [ -z $OSNAME ]; then
	echo 'Please define a variable $OSNAME depending your platform (linux|maxosx)'
	exit 1
fi
export OSNAME

#
# Some definition variables 
# 
ROOT=`pwd`
COUNT=0
TOTAL=0
BALISCORE=$ROOT/bali_score-$OSNAME-$HOSTTYPE
THRESHOLD=${THRESHOLD:-0.490}
FAILED=${FAILED:-$ROOT/failed}
SCRATCH=${SCRATCH:-$ROOT/scratch}


##
## The source dataset against to validate can be specified as the first argument on the command line 
## 
DATA=${1:-$ROOT/bb3/RV11}
if [[ ${DATA:0:1} != "/" ]]; then
DATA=$ROOT/$DATA 
fi 

#
# Clean before start
# 
rm -rf $SCRATCH; mkdir -p $SCRATCH
rm -rf $FAILED;  mkdir -p $FAILED
cd $SCRATCH

VERSION=`t_coffee -version`
echo "************************************************************************************************************"
echo "** Validation report for $VERSION"
echo "************************************************************************************************************"

##
## Go! 
## 

for TEST in  `ls  $DATA/*.tfa | xargs -I F basename F | cut -d '.' -f 1 | sort | uniq`
do 
  	echo -n "Validating $TEST - "
	COUNT=`echo "$COUNT +1" | bc`
 	
	## cleaning 
	set +e
	rm -rf * &> /dev/null
	rm -rf .* &> /dev/null
	set -e
	
	## creating soft linsk to target files 
	ln -s $DATA/$TEST* . 
	
	## run T-Coffee 
	t_coffee $TEST.tfa -output msf &> tcoffee.log
	
	
	## Validate
	$BALISCORE $TEST.xml $TEST.msf &> baliscore.log
	VALUE=`cat baliscore.log | grep auto | awk '{ print $4 }'`
	TOTAL=`echo "$VALUE + $TOTAL" | bc`
	RESULT=`echo "$VALUE >= $THRESHOLD" | bc`
	echo -n "bscore: $VALUE"
	
	if [ $RESULT -eq 1 ]; then 
	  echo "" 
	else 
	  echo "; warn!"
	  mkdir $FAILED/$TEST
	  mv * $FAILED/$TEST
	fi   
	
done

if [ $COUNT -gt 0 ]; then 
  AVG=`echo "scale=3; $TOTAL/$COUNT" |  bc -l` 
else 
  AVG=0
fi 

PASSED=`echo "($AVG >= $THRESHOLD)" | bc -l`
FAILED=`echo "($AVG < ($THRESHOLD * 0.9))" | bc -l`
echo ---------------------------------
echo -n "Validation average score: $AVG; "
if [ $FAILED -eq 1 ]; then
 	echo "FAILED"
	exit 1
fi

if [ $PASSED -eq 1 ]; then
 	echo "PASSED"
else
	echo "NOT PASSED"
fi

