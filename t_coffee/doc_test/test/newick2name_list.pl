#!/usr/bin/env perl


use Env qw(HOST);
use Env qw(HOME);
use Env qw(USER);


open (F,$ARGV[0]);
while ( <>)
  {
    @x=/([^:,;\)\(\s]+):[^:,;\)\(]*/g;
    @list=(@list,@x);
  }
$n=$#list+1;
foreach $n(@list){print ">$n\nsequence\n";}


close (F);

