#!/bin/bash 
set -e 

#apt-get update --fix-missing && apt-get install -y openjdk-7-jre-headless 
cd /root/tcoffee/
/root/tcoffee/lib/perl/lib/perl4makefile/doc2test.pl -replay /root/tcoffee/tests/.dumps/core.tests/
