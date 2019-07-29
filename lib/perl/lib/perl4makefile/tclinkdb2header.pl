#!/usr/bin/env perl
#Version 1.01 (25/02/03)
use Cwd;
use File::Path;
$SILENT=">/dev/null 2>/dev/null";

$cl=join( " ", @ARGV);

if ( ($cl=~/-db=\s*(\S+)/)){$db=$1;}
if ( ($cl=~/-infile=\s*(\S+)/)){$infile=$1;}
if ( ($cl=~/-mode=\s*(\S+)/)){$mode=$1;}

open (F, "$db");
while (<F>)
  {
    $l=$_;
    if (($l =~/^\/\//) || ($db=~/^#/)){;}
    elsif ( !($l =~/\w/)){;}
    else
      {
	my @v=split (/\s+/, $l);
	if ( $mode eq "Perl")
	  {
	    $header.="\$$v[0]\{\"$v[1]\"\}\{\"$v[2]\"\}=\"$v[3]\";\n";
	  }
	elsif ( $mode eq "C")
	  {
	    if (($l=~/^MODE/)){;}
	    elsif ( $v[2]eq"4_TCOFFEE")
	      {
	       
		$base{$v[1]}=$v[3];
		$header.="#define $base{$v[1]}_$v[2] \"$v[1]\"\n";
	      }
	    else
	      {
		$header.="#define $base{$v[1]}_$v[2] \"$v[3]\"\n";
	      }
	  }
      }
  }
close (F);

open (F, "$infile");
$in=0;
while (<F>)
  {
    $l=$_;
    
    if ($l=~/TclinkdbStart/){$in=1;}
    if ($l=~/TclinkdbEnd/){$in=0;}
    if ($in==1)
      {
	
	print "$l\n$header";
	$in=2;
      }
    if ( $in==0){print $l;}
  }
close (F);

