#!/usr/bin/env perl
use Cwd;
use File::Copy;
use File::Find;
use strict;
use DirHandle;
my ($source,$target,$delay,@filelist,$nfiles, $list, $sleep, $max, $nskip, $tskip);

$source="/Users/cnotredame/MusicPhotosBooksPiano/Audio/music/";
$target="/Users/cnotredame/Dropbox/personnal/Dropbox.Photos.Librairie.MusiqueLivres.Soft/Audio/music/";
$sleep="/Users/cnotredame/.sleep.txt";
$list="DefaultFileList.txt";
$delay=15;

for ($a=0; $a<@ARGV; $a++)
  {
    my $v=$ARGV[$a];
        
    if ($v eq "-s")
      {
	$source="$ARGV[++$a]/";
      }
    elsif ($v eq "-t")
      {
	$target=$ARGV[++$a];
      }
    elsif ($v eq "-delay")##ibooks to synchronyze itune lib
      {
	$delay=$ARGV[++$a];
      }
    elsif ($v eq "-list")
      {
	$list=$ARGV[++$a];
      }
    elsif ($v eq "-max")
      {
	$max=$ARGV[++$a];
      }
    elsif (int($v) eq $v)
      {
	if ($v<0){unlink($sleep);}
	else
	  {
	    open (F, ">$sleep");
	    print F "$v\n";
	    close (F);
	  }
	exit (0);
      }
  }


if    (!-d $source){die "Source: $source does not exist";}
elsif (!-d $target){die "Target: $target does not exist";}



my $n;
my (@dir) = ($source);


if (-e $list)
  {
    open (F, "$list");
    while (<F>)
      {
	my $l=$_;
	#print "SCAN: $l";
	chomp $l;
	push(@filelist,$l);
      }
    close(F);
  }
else
  {
    find(\&dir2list,@dir);
    open (F, ">$list");
    foreach my $f (@filelist){print F "$f\n";}
    close F;
  }

my $nfiles=@filelist;
my $nq;
foreach my $f (@filelist)
  {
    
    $nq+=queue_file($f, $source, $target);
    
  }

sub dir2list {
    #print $File::Find::name."\n";
    my $filename = $File::Find::name;

    $n++;
    print "DIR2LIST: $n $filename\n";
    $filename=~s/$source//;
    @filelist=(@filelist, $filename);
  }

sub queue_file
  {
    my ($f, $s, $t)=@_;
    my $to="$t/$f";
    my $from="$s/$f";
    
    if    (-d $from){if (!-d $to){mkdir "$t/$f";}}
    elsif (-f $from)
      {if (!-f $to)
	 {
	   my $size=-s $from;
	   $size/=1000000;
	   if ($max && $size>$max)
	     {
	       $nskip++;
	       $tskip+=$size;
	       print "SKIP $f [$size is Too Big][N=$nskip][T=$tskip Mb]\n";
	     }
	   else
	     {
	       my $p=int (($nq*100)/$nfiles);
	       print "QUEUE: $nq out of $nfiles [$f] : $p %\n";
	       copy ("$from", "$to");
	       my $wait=$size*$delay;
	       print "WAIT: $wait seconds (Size=$size). Edit [$sleep] to change this\n";
	       mysleep ($wait);
	     }
	 }
     }
    return 1;
  }


sub mysleep
  {
    my $d=shift;
    my $i=5;
    while ($d)
      {
	if (-e $sleep)
	  {
	    my $t=time();
	    open (F, "$sleep");
	    while (<F>){$d=$_;}
	    close (F);
	    chomp $d;
	    unlink $sleep;
	    print "FORCESLEEP: $d seconds ($t)\n";
	  }
	if   ($d<=$i){sleep $d; $d=0;}
	else {sleep $i; $d-=$i;}
      }
    return 1;
  }
	    
