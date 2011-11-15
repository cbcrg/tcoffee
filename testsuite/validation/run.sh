#!/bin/bash 
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
ERRORS=0
BALISCORE=$ROOT/bali_score-$OSNAME-$HOSTTYPE
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

##
## Go! 
## 

for TEST in  `ls  $DATA/*.tfa | xargs -I F basename F | cut -d '.' -f 1 | sort | uniq`
do 
  	echo -n "Validating $TEST:  "
	
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
	RESULT=`echo "$VALUE > 0.4" | bc`
	
	if [ $RESULT -eq 1 ]; then 
	  echo OK 
	else 
	  ERRORS=`echo "$ERRORS +1" | bc`
	  echo "FAILED (bali_score: $VALUE)"
	  mkdir $FAILED/$TEST
	  mv * $FAILED/$TEST
	fi   
	
done

if [ $ERRORS -ne 0 ]; then
   exit 1
fi
