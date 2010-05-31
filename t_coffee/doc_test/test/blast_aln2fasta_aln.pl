#!/usr/bin/env perl
use Env qw(HOST);
use Env qw(HOME);
use Env qw(USER);

$query_start=-1;
$query_end=-1;
open (F, ">tmp");
while (<>)
  {
    print F "$_";
    if ( /Sequences producing significant alignments/){$in_hit_list=1;}

    if ( /(\S+).*\s+(\d+)\s+([0-9e.-]+)/ && $in_hit_list)
      {	
	$name_list{$n_seq++}=$1;
	$name_list{$n_seq-1}=~s/.*\|//g;

      }
    if ( /QUERY/) {$in_aln=1;}
    if ( /(\S+)\s+(\d+)\s+([a-zA-Z-]+)\s+(\d+)/ || /(\S+)(\s+)(\-+)(\s+)/ && $in_aln)
      {
	
	$name=$1;
	$start=$2;
	$seq=$3;
	$end=$4;

	if ($1 eq "QUERY")
	  {
	    $n=0;
	    if ( $query_start==-1){$query_start=$start;}	    
	    $query_end=$end;
	  }
	else
	  {
	    $list{$n}{'real_name'}=$name_list{$n-1};
	  }

	$list{$n}{'name'}=$name;
	if ( !$list{$n}{'real_name'})
	  {
	    $list{$n}{'real_name'}=$name_list{$n+1};
	  }
	$seq=~tr/a-z/A-Z/;
	$list{$n}{'seq'}=$list{$n}{'seq'}.$seq;
	
	$n++;
	
      }
}

$list{0}{'real_name'}="QUERY_"."$query_start"."_"."$query_end";

for ($a=0; $a<$n; $a++)
  {
    print F ">$list{$a}{'real_name'}\n$list{$a}{'seq'}\n";

    print ">$list{$a}{'real_name'}\n$list{$a}{'seq'}\n";
  }
      

