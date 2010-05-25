#!/usr/bin/env perl 
# ----------------------------------------------------------------------
#
# Copyright (C) 2007 Andreas Wilm <andreas.wilm@ucd.ie>
#
# This file is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# ----------------------------------------------------------------------





my $PROBTRESH = 0.3;# base pairs below this prob threshold will be ignored
my $WEIGHT = 100.0; # float!!
my $NUCALPH = "ACGTUNRYMKSWHBVD";
use vars qw($NUCALPH $WEIGHT);

my $myname = basename($0);

use strict;
use warnings;
#use diagnostics;

use File::Basename;
use Getopt::Long;
use File::Glob ':glob';
use File::Spec;
use File::Temp qw/ tempfile tempdir /;




my $fmsq = "";
my $flib = "";
my %OPTS;
my %seq;
my ($id, $nseq, $i);
my @nl;

# ---   option handling
#
GetOptions("in=s" => \$fmsq, "out=s" => \$flib);

if (! -s $fmsq) {
    usage("empty or non-existant file \"$fmsq\"")
}
if (length($flib)==0) {
    usage("empty out-filename")
}



# ---   read seqences
#
%seq=read_fasta_seq($fmsq);

`Gibbs 10,10,10,10 -o x`;

&gibbs2aln (x);



sub read_fasta_seq 
  {
    my $f=$_[0];
    my %hseq;
    my (@seq, @com, @name);
    my ($a, $s,$nseq);

    open (F, $f);
    while (<F>)
      {
	$s.=$_;
      }
    close (F);

    
    @name=($s=~/>(\S*).*\n[^>]*/g);
    
    @seq =($s=~/>.*.*\n([^>]*)/g);
    @com =($s=~/>(\S*)(.*)\n([^>]*)/g);


    $nseq=$#name+1;
  
    for ($a=0; $a<$nseq; $a++)
      {
	my $n=$name[$a];
	my $s;
	$hseq{$n}{name}=$n;
	$s=$seq[$a];$s=~s/\s//g;
	
	$hseq{$n}{seq}=$s;
	$hseq{$n}{com}=$com[$a];
      }
    return %hseq;
  }

sub gibbs2aln
{
  my $f=@_[0];
  
  while (<F>)
    {
      $l=$_;
      chomp $l;
      $t.="$l;";
    }
  @list=($l=~/Num Motifs(.*)Column/g);
  foreach $e (@list){print "$e\n\n";}
}
