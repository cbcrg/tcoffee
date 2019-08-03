#!/bin/bash 
set -e 

#apt-get update --fix-missing && apt-get install -y openjdk-7-jre-headless 
cd /root/tcoffee/
lib/perl/lib/doc2test.pl -replay ./tests/.dumps/core.tests/
