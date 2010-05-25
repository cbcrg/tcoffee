#!/usr/bin/env perl

#Read the file
#records are in %record
#score are in %score

open F, $ARGV[0];
while (<F>)
  {
    
    if (/^\s+\n/){$nr++;}
    else
      {

	$record[$nr].=$_;
	if ( /Total prob of N-in:\s*([\d\.]+)/)
	  {
	    $score{$nr}=$1;
	  }
	$record{$nr}.=$_;	
      }
  }
close (F);

sub compare_score {$score{$b} <=> $score{$a};}
foreach $record_index (sort compare_score (keys(%score)))
  {
    print "\n$record{$record_index}";
  }

