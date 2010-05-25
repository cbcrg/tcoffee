#! /usr/bin/env perl

@file_list=@ARGV;
@code=@ARGV;

$n_code=$#file_list+1;
$date{'massachussets.txt'}=1788;
$date{'northcarolina.txt'}=1789;
$date{'kentucky.txt'}=1792;
$date{'alabama.txt'}=1819;
$date{'maine.txt'}=1820;
$date{'missouri.txt'}=1821;
$date{'iowa.txt'}=1846;
$date{'california.txt'}=1850;
$date{'oregon.txt'}=1859;
$date{'idaho.txt'}=1890;
$date{'alaska.txt'}=1959;


while ($#file_list>=0)
  {
    
    $i=$#file_list;
    $f=pop (@file_list);
    print "DIS_PROC Process $f\n";
    $text="";
    open (F, $f);
    while (<F>)
      {
	$l=$_;
	chomp $l;
	$text.=" $l";
      }
    close (F);
    $text=clean_text ($text);
    
    %list=file2di_word_list ($i,$text, %list);
  }

foreach $w (keys %{$list{TOT}})
  {
    
    for ($b="",$a=0; $a< $n_code; $a++)
      {
	if ($list{$a}{$w} && $list{$a}{$w}>0){$b.="1";}
	else {$b.="0";}
      }
    $list{PROFILE}{$w}=$b;
    #if ( $b=~/0/){print "$b $w\n";}
  }



#for ( $a=0; $a< $n_code; $a++)
#  {
#    print "\n>$code[$a]\n";
#    foreach $w (keys %{$list{TOT}})
#      {
#	if ( !($list{PROFILE}{$w}=~/0/)){next;}
#	
#	for ($b=0; $b<$n_code; $b++)
#	  {
#	    
#	    if ($b==$a)
#	      {
#		if ($list{$a}{$w}){print "W";}
#		else 
#		  {
#		    print "-";
#		  }
#	      }
#	  }
#      }
#  }
#print "\n";


print "# TC_DISTANCE_MATRIX_FORMAT_01\n";
for ( $a=0; $a<$n_code; $a++)
  {
    print "# SEQ_INDEX $code[$a] $a\n";
  }
for ($a=0; $a< $n_code-1; $a++)
  {
    for ( $b=$a+1; $b< $n_code; $b++)
      {
	$score=0;
	$best_w="";
	$best_count=0;
	
	for ($d="",$c=0; $c<$n_code; $c++)
	  {
	    if ($c==$a || $c==$b){$d.="1";}
	    else {$d.="0";}
	  }
	foreach $w (keys %{$list{TOT}})
	  {
	    $count=$list{TOT}{$w};
	   
	    
	    if ($d eq $list{PROFILE}{$w})
	      {
		$score++;
		if ( $count>$best_count && $w ne "nwords")
		   {
		     $best_count=$count;
		     $best_w=$w;
		   }
	      }
	  }
	$st1=$code[$a];
	$st2=$code[$b];
	
	$d1=$date{$st1};
	$d2=$date{$st2};
	
	$diff=($d1<$d2)?2000-$d2:2000-$d1;
	$diff = abs ((2000-$date{$st1})+(2000-$date{$st2}));
	
	
	$score=(2*$score*100)/($list{$a}{nwords}+$list{$b}{nwords});
	
	if ( $st1 ne $st2)
	  {
	    printf "DIS_DIS %-20s %-20s %5.2f %2d\nDIS_BWO $best_w $best_count\n", $code[$a],$code[$b],$score,$diff;
	    printf "BOT $a $b %.2f $code[$a] $code[$b] %.2f\n", $score*10, $score*10;
	    printf "TOP $b $a %.2f $code[$b] $code[$a] %.2f\n", $score*10, $score*10;
	  }
	
      }
  }
    

die;

for ( $c1=0; $c1<$n_code; $c1++)
  {
    for ( $c2=$c1+1; $c2<$n_code; $c2++)
      {
	$score=CompareList ($c1, $c2,%list);
	print "\n\n";
	printf "\nSCORE %d Vs %d =%.3f\n", $c1, $c2 ,$score;
      }
  }


sub clean_text
  {
    my $t=shift @_;
    
    $t=~s/\<[^\>]*\>/ /g;
    $t=~s/\d//g;
    $t=~s/["'.;,:-]/ abcdefg /g;
    $t=~s/[^\w]/ /g;
    $t=~s/\s\S\s/ /g;
    
    $t=~s/nbsp/ /g;
    $t=lc $t;

    
    return $t;
    
    }
sub file2word_list
  {
    my $i=shift (@_);
    my $l=shift (@_);
    my %wl=@_;
    
    

    $l=lc $l;
    @wordlist=($l=~/(\w+)/g);
    foreach $word1 (@wordlist)
      {
	$word=lc $word1;
	
	if ($wl{TOT}{$word}){$wl{$i}{$word}++;$wl{TOT}{$word}++; }
	else
	  {
	    $wl{$i}{$word}++;
	    $wl{TOT}{$word}++;
	    $wl{$i}{nwords}++;
	    $wl{TOT}{nwords}++;
	    
	    $nwords=$wl{TOT}{nwords};
	    $wl{wordlist}{$nwords}=$word;
	  }
      }
    return %wl;
}

sub file2di_word_list
  {
    my $i=shift (@_);
    my $l=shift (@_);
    my %wl=@_;
    my %sw;
    my $word, $word1, $word2;
    
    %sw=file2word_list ($i, $l, %sw);
    $l=lc $l;
    @wordlist=($l=~/(\w{4,})/g);

    $max_single_word=10;
    for ($a=0;$a<$#wordlist; $a++)
      {
	
	$word1=$wordlist[$a];
	$word2=$wordlist[$a+1];

	if ( $word1 eq "abcdefg" || $word2 eq "abcdefg"){next;}
	elsif ($sw{$i}{$word1}>$max_single_word || $sw{$i}{$word2}>$max_single_word ){next;}
	
	$word="$word1 $word2";	
	$word=lc $word;

	
	if ($wl{TOT}{$word}){$wl{$i}{$word}++;$wl{TOT}{$word}++; }
	else
	  {
	    $wl{$i}{$word}++;
	    $wl{TOT}{$word}++;
	    $wl{TOT}{nwords}++;
	    $wl{$i}{nwords}++;
	    $nwords=$wl{nwords};
	    $wl{wordlist}{$nwords}=$word;
	  }
      }
    return %wl;
}
sub file2period_list
  {
    my $i=shift (@_);
    my $l=shift (@_);
    my %wl=@_;
    
    
    @sentences=($l=~/#P.#(.+)/g);
    foreach $s (@sentences)
      {
	print "XXXXXX $s XXXXXXX\n";
      }
    }
