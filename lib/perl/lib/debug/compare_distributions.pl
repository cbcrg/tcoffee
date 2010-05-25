#!/usr/bin/env perl

$HOME=$ENV{'HOME'};

#defaults
$archive="$HOME/untared_distributions/";
$package="T-COFFEE";
$source="t_coffee_source";
$action="diff";

$v1=$ARGV[0];
$v2=$ARGV[1];
$file=$ARGV[2];


if ( $#ARGV==-1)
  {
    print "DESCRIPTION:\n";
    print "  Compares|Reconciliates two distributions\n";
    print "  It is assumed all the distributions are\n";
    print "  in the same dir and use the same file structure\n";
    print "  The distribution source is expected to be in the form:\n";
    print "      <archive>/<package>_distribution_Version_<v>/<source>/\n";
    print "USAGE:\n";
    print "   compare_distributions <version 1> <version 2> <file>\n";
    print "PARAMETERS:\n";
    print "  <file>.................name of a source file OR file\n" ;
    print "  .......................containing source file names (NO PATHS)\n";
    print "  -diff..................makes a diff on the considered files\n";
    print "  -merge.................makes a merge on the considered files\n";
    print "  -soft..................name of the package (NO PATH!!)\n";
    print "  -source................name of the source dir (NO PATH!!)\n";
    print "  -archive...............name of the dir containing the distributions\n";
    
    print "DEFAULTS:\n";
    print "  -action=$action\n";
    print "  -soft=$package\n";
    print "  -source=$source\n";
    print "  -archive=$archive\n";
    die;
  }



foreach $arg (@ARGV)
  {
    $cl.=" $arg";
  }


 
if ( $cl=~/\-soft (\S+)/){$package =$1;}
if ( $cl=~/\-source (\S+)/){$source =$1;}
if ( $cl=~/\-archive (\S+)/){$archive =$1;}
     

if ( $cl=~/\-merge/){$action ="cp";}
elsif ($cl=~/\-diff/){$action ="diff";}

if ( -e $file)
  {
    open (F, $file);
    while (<F>)
      {
	$f= $_;
	chop $f;
	&act ($v1, $v2, $f, $action);
      }
    close (F);
  }
else
  {
    &act ($v1, $v2, $file, $action);
  }

sub act 
  {
    my $v1=@_[0];
    my $v2=@_[1];
    my $f=@_[2];
    my $action=@_[3];
    my $command;
    my $file1, $file2;
    
    print "\n###########################################\n $f\n";
   
    $file1="$archive/$package\_distribution_Version_$v1/$source/$f";
    $file2="$archive/$package\_distribution_Version_$v2/$source/$f";
    $command="$action $file1 $file2";
    print "\n$command\n";
    
    system ("$command >x");
    open (Z, "x");
    while (<Z>){print "$_";}
    close (Z);
    print "\n###########################################";
    
    
    
  }
