#!/usr/bin/env perl
use strict;

#This program will set the version tag to whatever is in $ARGV[0] or the file ./version_number-version
#if on MASTER system, it will set this value to the Github branch value
#it will return stdout

my $file;
my $version;



if( $ARGV[0]){$file=$ARGV[0];}
else{$file="version_number.version";}


if ($ENV{'TC_MASTER_NODE'})
  {
    system ("git branch -v  > tmp.txt");
    open (F, "tmp.txt");
    while(<F>)
      {
	my $l=$_;
	if ($l=~/\* master (\w+) .*/)
	  {
	    $version=$1;
	  }
      }
    close (F);
    open (F, ">$file");
    print F "$version";
    close F;
    unlink ("tmp.txt");
  }
open (F,$file);

while (<F>)
          {
	  my $line=$_;
	  chop ($line);
	  print "$line";
	  }
close (F);

