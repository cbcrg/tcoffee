#!/usr/bin/perl
#
#
#date : 12/10/

use Env qw(HOST);
use Env qw(HOME);
use Env qw(USER);



open (F, $ARGV[0]);
while ( <F>)
  {
    if ( /Database: All/){$record=0;}
    if ( $record){$line[$nl++]=$_;}
    if ( /ALIGNMENTS/){$record=1;}
  }

foreach $l ( @line)
  {
    @list=($l=~/(\S+)/g);
    
    if ($#list>2)
      {
	$name=$list[0];
	if (!$hasch{$name}{'index'})
	  {
	    $s=++$hasch{nseq};
	    $hasch{$name}{'index'}=$s;
	    $hasch{$s}{'name'}=$name;    
	  }
	
	$seq=substr($l,16,60);
	$seq=~s/\s/\-/g;
	$seq=~s/\d//g;
	
	if ( $name eq "Query")
	  {
	    $block=++$hasch{'nblock'};
	    @x=($seq=~/(.)/g);
	    $hasch{$block}{'block_len'}=$#x+1;
	  }
	$hasch{$name}{'block'}{$block}{'seq'}=$seq;
      }
  }


for ($a=1; $a<=$hasch{nseq}; $a++)
  {
    $name=$hasch{$a}{'name'};
    print ">$name\n";
    for ($b=1; $b<=$hasch{'nblock'}; $b++)
      {
	if ($hasch{$name}{'block'}{$b}{'seq'})
	  {

	    print "$hasch{$name}{'block'}{$b}{'seq'}\n";
	  }
	else
	  {
	    $blank=generate_gap($hasch{$b}{'block_len'});
	    print "$blank\n";
	    }
	}
    }

sub generate_gap
  {
    my $l=@_[0];
    my $blank;
    my $a;
    
    for ( $a=0; $a<$l; $a++)
      {
	$blank.='-';
      }
    return $blank;
    }
