#!/usr/bin/env perl

$dir="$ARGV[0]/";
opendir (DIR,$ARGV[0]);
@files=readdir (DIR);
closedir (DIR);
$add=$ARGV[1];

print $dir;
foreach $f (@files)
  {
    
    rename "$dir/$f", "$dir/0$f";
  }
