#!/usr/bin/env perl
use Env;
use FileHandle;
use Cwd;
use File::Path;
use Sys::Hostname;
my $f = new FileHandle;

open ($f, $ARGV[1]);
$atom=$ARGV[0];

$atom=~s/PRIME/\'/;
while (<$f>)
  {
    my $l=$_;

    $l=~s/$atom/CA /;
    
    
    $l=~s/  G /GLY /g;
    $l=~s/  C /CYS /g;
    $l=~s/  T /THR /g;
    $l=~s/  A /ALA /g;
    $l=~s/  U /THR /g;
    
    $l=~s/ DG /GLY /g;
    $l=~s/ DC /CYS /g;
    $l=~s/ DT /THR /g;
    $l=~s/ DA /ALA /g;
    $l=~s/ DU /THR /g;
    
    print $l;
  }




