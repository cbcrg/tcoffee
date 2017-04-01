#!/usr/bin/env perl
#This script reads an aln in fasta format and modifies names so thathere are no duplicates
#if two sequences have the same name, the second one is renamed name_1 and so on
use strict;
use FileHandle;
use Env qw(HOST);
use Env qw(HOME);
use Env qw(USER);
my %name;
my $nseq;
my $F= new FileHandle;
open ($F, $ARGV[1]);
while(<$F>)
  {
    my $l=$_;
    if ($l=~/^#/){;}
    elsif (($l=~/\d+\s+\d+\s+(\S+)\s+(\S+)/))
      {
	my $n=$1;
	$name{$1}++;
      }
  }
close ($F);

open ($F, $ARGV[0]);
while(<$F>)
  {
    my $l=$_;
    if ($l=~/^#/){;}
    elsif ($l=~/\d+\s+\d+\s+(\S+)\s+(\S+)/)
      {
	my $n=$1;
	$name{$n}++;
	if ($name{$n}==2){$nseq++;}
      }
  }
close ($F);

print "#NAMESEQ_01\n";
print "# $nseq\n";
open ($F, $ARGV[0]);
while(<$F>)
  {
    my $l=$_;
    if ($l=~/^#/){;}
    elsif ($l=~/\d+\s+\d+\s+(\S+)\s+(\S+)/)
      {
	my $n=$1;
	if ($name{$n}==2){print "$l";}
      }
  }
close ($F);
exit (0);


