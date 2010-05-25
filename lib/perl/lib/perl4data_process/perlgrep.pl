#!/usr/bin/env perl
use File::Copy;
use Env qw(HOST);
use Env qw(HOME);
use Env qw(USER);




while (<STDIN>)
  {
    $line=$_;
    @list=($line=~/(\S+)/g);
    for ( $a=0; $a<=$#ARGV; $a++)
      {
	if ( $ARGV[$a]=~/\-(\d+)/)
	  {
	    $f=$1;
	    if ($ARGV[$a]=~/\>(\d+)/)
	      {if ($list[$f-1]>$1){$keep=1;}}
	    elsif ($ARGV[$a]=~/\<(\d+)/)
	      {if ($list[$f-1]<$1){$keep=1;}}
	    elsif ($ARGV[$a]=~/\=(\d+)/)
	      {if ($list[$f-1]==$1){$keep=1;}}
	    elsif ($ARGV[$a]=~/\==(\S+)/)
	      {if ($list[$f-1]=~/$1/){$keep=1;}}
	    elsif ($ARGV[$a]=~/\=(\S+)/)
	      {if ($list[$f-1] eq $1){$keep=1;}}
	    else {$keep=0;}
	    
	  }
      }
    if ( $keep==1){print $line;}
  }
