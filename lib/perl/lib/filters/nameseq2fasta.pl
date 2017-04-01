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
open ($F, $ARGV[0]);
while(<$F>)
  {
    
    my $l=$_;
    if ($l=~/^#/){;}
    elsif (($l=~/\d+\s+\d+\s+(\S+)\s+(\S+)/))
      {
	my $name=$1;
	my $seq=$2;
	print ">$name\n$seq\n";
      }
  }
close ($F);
exit (0);


