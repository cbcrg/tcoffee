#!/usr/bin/env perl
use Env qw(HOST);
use Env qw(HOME);
use Env qw(USER);


$query_start=-1;
$query_end=-1;

while (<>)
  {
    if ( /\/\//){$in_aln=1;}
    elsif ( $in_aln && /(\S+)\s+(.*)/)
      {


	$name=$1;
	

	$seq=$2;
	$seq=~s/\s//g;
        $seq=~s/\~/\-/g;
	$seq=~s/\./\-/g;
	if ( $list{$n}{'name'} && $list{$n}{'name'} ne $name)
	  {
	    print "$list{$n}{'name'} Vs $name";
	    
	    exit (EXIT_FAILURE);
	  }
	else
	  {
	    $list{$n}{'name'}= $name;
	  }

	$list{$n}{'seq'}=$list{$n}{'seq'}.$seq;
	
	$nseq=++$n;
	
      }
    else
      {$n=0;}
  }


for ($a=0; $a<$nseq; $a++)
  {
    print ">$list{$a}{'name'}\n$list{$a}{'seq'}\n";
  }
      

