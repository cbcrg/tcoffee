#!/bin/bash
# set path to java using JAVA_HOME if available, otherwise assume it's on the PATH
JAVA_PATH=${JAVA_HOME:+$JAVA_HOME/bin/}java

if [ -z $DIR_4_TCOFFEE ]; then
	echo You need to define in your enviroment the variable '$DIR_4_TCOFFEE' referring the root directory of your T-Coffee installation
	exit 1
fi 

$JAVA_PATH -jar ./black-coffee.jar -var tcoffee.home=$DIR_4_TCOFFEE "$@"
