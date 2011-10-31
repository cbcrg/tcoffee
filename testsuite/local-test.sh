#!/bin/bash
# set path to java using JAVA_HOME if available, otherwise assume it's on the PATH
JAVA_PATH=${JAVA_HOME:+$JAVA_HOME/bin/}java

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

##
## Define the DIR_4_TCOFFEE variable and validate the T-Coffee home path
## 
if [ -z $DIR_4_TCOFFEE ]; then
	BIN=`which t_coffee`
	DIR_4_TCOFFEE=`dirname $BIN`
	MCOFFEE_4_TCOFFEE=$DIR_4_TCOFFEE/mcoffee/
	PLUGINS_4_TCOFFEE=$DIR_4_TCOFFEE/plugins/$OSNAME/
	PERL5LIB=$DIR_4_TCOFFEE/perl/
	MAFFT_BINARIES=$BIN
fi 

if [ ! -e $DIR_4_TCOFFEE/bin ]; then 
	echo Missing path $MAFFT_BINARIES
	exit 2
fi 
if [ ! -e $DIR_4_TCOFFEE/mcoffee ]; then 
	echo Missing path $DIR_4_TCOFFEE/mcoffee
	exit 2
fi 
if [ ! -e $DIR_4_TCOFFEE/plugins/$OSNAME/ ]; then 
	echo Missing path $DIR_4_TCOFFEE/plugins/$OSNAME/
	exit 2
fi 
if [ ! -e $DIR_4_TCOFFEE/perl/ ]; then 
	echo Missing path $DIR_4_TCOFFEE/perl/
	exit 2
fi 

##
## Launch the tests
##
$JAVA_PATH -jar ./black-coffee.jar -var tcoffee.home=$DIR_4_TCOFFEE "$@"
