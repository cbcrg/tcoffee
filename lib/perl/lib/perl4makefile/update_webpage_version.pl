#!/usr/bin/env perl
#This program increments a tag VersionTag by 0.1
#It respects any symbols AROUND VersionTag = xx.xx



$file=$ARGV[0];
$version=$ARGV[1];

open (F, $file);
while ( <F>)
  {
    $l=$_;
    if ( $l=~/Version_/)
      {
	$l=~s/(Version_[^<]*)/$version/;
	print "$l";
      }
    else
      {print $l;}
  }
