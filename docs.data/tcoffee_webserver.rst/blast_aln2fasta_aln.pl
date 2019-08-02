#!/usr/bin/env perl

use Env qw(HOST);
use Env qw(HOME);
use Env qw(USER);

                                                        
use strict;                                             
use warnings;
use diagnostics;

my $in_hit_list, my $in_aln=0, my(%name_list)=(),my (%list)=(),my $n_seq=0; my $test=0;
my($j)=0, my $n=0, my $nom, my $lg_query, my %vu=();

open (F, ">tmp");

$/="\n";
while (<>)
{
    print F $_;
    if($_ =~ /Query=\s*(.+?)\s/i) { $nom=$1;}

    if ( /Sequences producing significant alignments/){$in_hit_list=1;}
    
    if ($_=~ /^pdb\|/i) { $_=~ s/pdb\|//g; }
    if ($_=~ /^(1_\d+)\s+\d+/) { $_=~ s/$1/QUERY/;}
      
    if ( /^(\S+).+?\s+[\d.]+\s+([\de.-]+)\s+$/ && $in_hit_list)	
    {
	my($id)=$1; # 
	$id=~ s/\|/_/g; #
	if ($id =~ /.+_$/) { chop($id) }; #
	$name_list{$n_seq++}=$id;
	$name_list{$n_seq-1}=~ s/.*\|//g;     
    }
  
    if (/query/i) {$in_aln=1;}
    if ( /^(\S+)\s+(\d+)\s+([a-zA-Z-]+)\s+(\d+)/ || /^(\S+)(\s+)(\-+)(\s+)/ && ($in_aln == 1))
    {
	my $name=$1;
	my $start=$2;
	my $seq=$3;
	my $end=$4;
		
	if ($name =~ /QUERY/i) { $lg_query=length($seq); }

	unless ($test > $n) #m
	{
	    my(@seqq)= split('',$seq);
	    my($gap_missing)= scalar(@seqq);
	    
	    while ($gap_missing != $lg_query)  { unshift (@seqq,"-"); $gap_missing= scalar(@seqq); }
	    $seq=join('',@seqq);  #m
	}
	
	if ($name =~ /QUERY/i)
	{
	    $n=0; %vu=(); $j=0;
	    $list{$n}{'real_name'}="$nom";
	}	
	else
	{
	    unless (exists $vu{$name}) { ++$j;}	
	    $list{$n}{'real_name'}=$name_list{$j-1};
	}
		
	$list{$n}{'name'}=$name;

	$seq=~tr/a-z/A-Z/;
	$list{$n}{'seq'}=$list{$n}{'seq'};
	$list{$n}{'seq'}.=$seq;

	$n++;
	$vu{$name}++;
	$test++;
   } 
    
}

my @numero=();

for (my $a=0; $a<$n; $a++) #m
{
    my $long=length($list{0}{'seq'});  
    my $long1= length($list{$a}{'seq'});
  
    while ($long1 ne $long)
    {
	$list{$a}{'seq'}.="-";
	$long1= length ($list{$a}{'seq'});
    } 
 
    push (@numero,"$list{$a}{'name'} $list{$a}{'real_name'}\n");
}

my %dejavu=();


for (my $i=0; $i<=$#numero; $i++)
{
    my $s=">$list{$i}{'real_name'}\n$list{$i}{'seq'}\n";
    my $k=0;
    
    if (exists $dejavu{$numero[$i]}) {next;}
    else
    {	
	for ($j=0; $j<$n ; $j++)
	{
	    if ("$numero[$i]" eq "$numero[$j]" && $j != $i )
	    {
		++$k;
		$s .=">$list{$j}{'real_name'}\n$list{$j}{'seq'}\n";
	    }
	}	
    }
    
    if ($k>0) 
    {
	my $cons;
	open (SOR,">tempo_aln2cons"); print SOR $s;  close SOR ;
	open (COM,"t_coffee -other_pg seq_reformat -in tempo_aln2cons -action +aln2cons +upper |") ; 
     	while (<COM>)
	{	
	    if (/^>/) { $cons =">$list{$i}{'real_name'}\n"; next;}
	    $_=~ s/\n//g;
	    $cons .=$_;
	}
	close COM; unlink ("tempo_aln2cons");
	print $cons,"\n"; print F $cons,"\n";
    }	
    else  { print $s;  print F $s; }
    
    $dejavu{$numero[$i]}++;
} #m

exit;












