#!/bin/bash 
set -e 

apt-get update --fix-missing && apt-get install -y openjdk-7-jre-headless 
cd /root/tcoffee/testsuite/
java -jar black-coffee.jar \
  --stop error \
  --delete never \
  --sandbox-dir /test-results/all \
  --html-path-prefix all \
  -o /test-results/index.html \
  -R ./all