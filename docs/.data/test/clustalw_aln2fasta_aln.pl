#!/usr/bin/env perl

use Env qw(HOST);
use Env qw(HOME);
use Env qw(USER);


$tmp=clean_cr ($ARGV[0]);
open (F, $tmp);

while ( <F>)
  {
    my $l=$_;
    if ( $l=~/^# STOCKHOLM/){$stockholm=1;}
    elsif ( $stockholm && $l=~/^#/)
      {
	$l=~/^#(\S+)\s+(\S+)\s+(\S*)/g;
	$l="_stockholmhasch_$1\_stockholmspace_$2 $3\n";
      }
    $file.=$l;
  }
close (F);
unlink($tmp);
$file1=$file;

$file=~s/\#/_hash_symbol_/g;
$file=~s/\@/_arobase_symbol_/g;


$file=~s/\n[\.:*\s]+\n/\n\n/g;

$file=~s/\n[ \t\r\f]+(\b)/\n\1/g;


$file=~s/(\n\S+)(\s+)(\S)/\1_blank_\3/g;

$file=~s/[ ]//g;
$file=~s/_blank_/ /g;



$file =~s/\n\s*\n/#/g;

$file.="#";
$file =~s/\n/@/g;




@blocks=split /\#/, $file;
shift (@blocks);
@s=split /\@/, $blocks[0];
$nseq=$#s+1;



$file=join '@', @blocks;
@lines=split /\@/,$file;

$c=0;

foreach $l (@lines)
  {
    if (!($l=~/\S/)){next;}
    elsif ($stockholm && ($l=~/^\/\// || $l=~/STOCKHOLM/)){next;}#get read of STOCHOLM Terminator
   
    $l=~/(\S+)\s+(\S*)/g;
    $n=$1; $s=$2;
    
    $seq[$c].=$s;
    $name[$c]=$n;
    $c++;
    
    if ( $c==$nseq){$c=0;}
    
  } 

if ( $c!=0)
      {
	print STDERR "ERROR: $ARGV[0] is NOT an MSA in Clustalw format: make sure there is no blank line within a block [ERROR]\n";
	exit (EXIT_FAILURE);
      }

for ($a=0; $a< $nseq; $a++)
  {
    $name[$a]=cleanstring ($name[$a]);
    $seq[$a]=cleanstring ($seq[$a]);
    $seq[$a]=breakstring($seq[$a], 60);
    
    $line=">$name[$a]\n$seq[$a]\n";
    
    print "$line";
  }
exit (EXIT_SUCCESS);

sub cleanstring
  {
    my $s=@_[0];
    $s=~s/_hash_symbol_/\#/g;
    $s=~s/_arobase_symbol_/\@/g;
    $s=~s/[ \t]//g;
    return $s;
  }
sub breakstring
  {
    my $s=@_[0];
    my $size=@_[1];
    my @list;
    my $n,$ns, $symbol;
    
    @list=split //,$s;
    $n=0;$ns="";
    foreach $symbol (@list)
      {
	if ( $n==$size)
	  {
	    $ns.="\n";
	    $n=0;
	  }
	$ns.=$symbol;
	$n++;
      }
    return $ns;
    }

sub clean_cr
  {
    my $f=@_[0];
    my $file;
    
    $tmp="f$.$$";
    
    
    open (IN, $f);
    open (OUT, ">$tmp");
    
    while ( <IN>)
      {
	$file=$_;
	$file=~s/\r\n/\n/g;
	$file=~s/\n\r/\n/g;
	$file=~s/\r\r/\n/g;
	$file=~s/\r/\n/g;
	print OUT "$file";
      }
    
    close (IN);
    close (OUT);
    return $tmp;
  }

