#!/usr/bin/env perl
use File::Copy;
use Env qw(HOST);
use Env qw(HOME);
use Env qw(USER);

foreach $v (@ARGV){$cl.=$v;}


if ( $cl=~/-k(\d+)/){$k=$1;}
else {$k=1;}
if ( $cl=~/-w(\d+)/){$w=$1;}
else {$w=-1;}
if ( $cl=~/-p(\d+)/){$p=$1;}
else {$p=-1;}

while (<STDIN>)
  {
    @l=($_=~/(\S+)/g);
    $v=$l[$k-1];
    if ( !$h{$v}){@ll=($v, @ll);}
    
    if ( $w==-1)
      {$h{$v}++;}
    else
      {$h{$v}+=$l[$w-1];}

    if ($p!=-1){$print{$v}=$l[$p-1];}

  }
foreach $v (@ll)
  {
    print "$v $print{$v} $h{$v}\n";
  }

