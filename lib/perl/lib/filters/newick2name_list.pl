#!/usr/bin/env perl
#This script reads a newik tree and outputs a list of names:
#n seq
#name1
#name2



use Env qw(HOST);
use Env qw(HOME);
use Env qw(USER);
#Specific Files and Default Parameters


open (F,$ARGV[0]);
while ( <>)
  {
    @x=/([^:,;\)\(\s]+):[^:,;\)\(]*/g;
    @list=(@list,@x);
  }
$n=$#list+1;
foreach $n(@list){print ">$n\nsequence\n";}


close (F);
