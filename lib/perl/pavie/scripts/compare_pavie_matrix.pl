#!/usr/bin/env perl
#Takes a list of file where each file contains one sequence in fasta format
#Makes all the pairwise alignments
#outputs a matrix

%mat1=read_matrix ($ARGV[0]);
%mat2=read_matrix ($ARGV[1]);
@symbol_list=mat2symbol_list($ARGV[1]);


foreach $s1 (@symbol_list)
  {
    foreach $s2 (@symbol_list)
      {
	$n++;
	$diff=$mat1{$s1}{$s2}-$mat2{$s1}{$s2};
	$tot+=($diff*$diff);
      }
  }
$diff=sqrt ($tot/$n);
printf "%.3f",$diff;


sub read_matrix
  {
    my $file=@_[0];
    my %mat;
    
    open (F,$file);
    while ( <F>)
      {
	
	if ( /Symbol_list/){;}
	
	elsif (/GAP/)
	  {
	    
	    /GAP\s+(.+)/;
	    $mat{'gap'}{'gap'}=$1;
	  }
	else
	  {
	    
	    /(.) (.)\s+(.+)/;
	    $mat{$1}{$2}=$3;
	  }
      }
    close (F);
    return %mat;
  }
    
    
sub mat2symbol_list
  {
    my $file=@_[0];
    my @symbol_list;
    
    open (F,$file);
    while ( <F>)
      {
	
	if ( /Symbol_list/)
	  {
	    $l=$_;
	    $l=~s/Symbol_list: //;
	    
	    @symbol_list=($l=~/(\S)/g);
	  }
      }
    close (F);
    return @symbol_list;
  }
    
    
