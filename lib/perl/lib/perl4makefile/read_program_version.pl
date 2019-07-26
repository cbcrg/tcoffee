#!/usr/bin/env perl
use strict;

#This program will set the version tag to whatever is in $ARGV[0] or the file ./version_number-version
#if on MASTER system, it will set this value to the Github branch value
#it will return stdout

my $file;

my $versionF;
my $branchF;

my $version;
my $branch;



if( $ARGV[0]){$file=$ARGV[0];}
else{$file="version_number";}

$versionF=$file.".version";
$branchF=$file.".branch";


if ($ENV{'TC_MASTER_NODE'})
  {
    system ("git branch -v  > tmp.txt");
    open (F, "tmp.txt");
    while(<F>)
      {
	my $l=$_;
	if ($l=~/\* master (\w+) .*/)
	  {
	    $branch=$1;
	  }
      }
    close (F);
    open (F, ">$branchF");
    print F "$branch";
    close F;
    unlink ("tmp.txt");
  }


#get version
open (F,$versionF);
while (<F>)
          {
	  $version=$_;
	  chop ($version);
	}
close (F);



#get branch
open (F,$branchF);
while (<F>)
          {
	  $branch=$_;
	  chop ($branch);
	}
close (F);
print "$branch";
