#!/usr/bin/env perl
#

# define scoring system:
# usage: pg seq1 seq2 matrix gap penalty

@arg_list=@ARGV;
if ( $#arg_list==-1)
  {
    print "Usage: pg seq1 seq2 <matrix>\n";
    die;
  }


if ( !$arg_list[2])
{
  @w=([  0,-10,-10,-10,-10,-10,-10, -3 ], 
    [-10,  0,-10,-10,-10,-10,-10, -3 ], 
    [-10,-10,  0,-10,-10,-10,-10, -3 ],
    [-10,-10,-10,  0,-10,-10,-10, -3 ],
    [-10,-10,-10,-10,  0,-10,-10, -3 ],
    [-10,-10,-10,-10,-10,  0,-10, -3 ],
    [-10,-10,-10,-10,-10,-10,  0, -3 ],
    [ -3, -3, -3, -3, -3, -3, -3, -3 ]);
}
elsif ( $arg_list[2])
  {
    open (F, $arg_list[2]);
    while (<F>)
      {
	if ( /Symbol/){;}
	elsif (/GAP/)
	  {
	  
	    /GAP\s+(.+)/;
	    for ( $a=0; $a<5; $a++){$g[$a]=$1;}
	  }
	else
	  {
	    
	    /(.) (.)\s+(.+)/;
	    
	    $mv1=uc $1;$mv2= uc $2;
	    $mv1=~ tr/ABCDEFGX/0-7/;
	    $mv2=~ tr/ABCDEFGX/0-7/;
	    
	    $w[$mv1][$mv2]=$3;
	  }
      }
    close (F);
  }

if ( !@g)
  {
    @g=(  -5, -5, -5, -5, -5, -5, -5, -5  );
  }


# input parsing

$/ = "\n>"; 

open(SEQ1, "$arg_list[0]");

while(<SEQ1>) {
  print "$_";
  /^>?([^\n]*)\n(.*)/s;
  $name = $1;
  $seq = $2;
  $seq =~ tr/[a-z]/[A-Z]/;
  $seq =~ s/[^ABCDEFGX]//gs;
  $seq =~ tr/ABCDEFGX/0-7/;
  @x = split //, " " . $seq;
  $XL = scalar(@x) - 1;
  print "\nSequence 1: $name\n\n";
  print "  ", @x, "\n";
  }
   
open(SEQ2, "$arg_list[1]");
while(<SEQ2>) {
  /^>?([^\n]*)\n(.*)/s;
  $name = $1;
  $seq = $2;
  $seq =~ tr/[a-z]/[A-Z]/;
  $seq =~ s/[^ABCDEFGX]//gs;
  $seq =~ tr/ABCDEFGX/0-7/;
  @y = split //, " " . $seq;
  $YL = scalar(@y) - 1;
  print "\nSequence 2: $name\n\n";
  print "  ", @y, "\n";
  }   

# socio- alignment algorithm fill-stage


$F[0][0] = 0;
for($i = 1; $i <= $XL; $i++) {
  $F[$i][0] = $F[$i-1][0]+$g[$x[$i]]; 
  }
for($j = 1; $j <= $YL; $j++) {
  $F[0][$j] = $F[0][$j-1]+$g[$y[$j]]; 
  }

for($i = 1; $i <= $XL; $i++) {
  for($j = 1; $j <= $YL; $j++) {
    $s[$i][$j] = $w[$x[$i]][$y[$j]];

    $F[$i][$j] = $F[$i-1][$j-1]+$s[$i][$j] ;
    if          ($F[$i-1][$j] + $g[$x[$i]] > $F[$i][$j]) {
    $F[$i][$j] = $F[$i-1][$j] + $g[$x[$i]] }
    elsif       ($F[$i][$j-1] + $g[$y[$j]] > $F[$i][$j]) {
    $F[$i][$j] = $F[$i][$j-1] + $g[$y[$j]] }
    }
  }

# print filled SW matrix

  print "\nSmith-Waterman score matrix:\n";
  for($i = 0; $i <= $XL; $i++) {
    print "\n";
    for($j = 0; $j <= $YL; $j++) {
      printf "%4i", $F[$i][$j]
      }
    }
  print "\n";

# Smith-Waterman algorithm: traceback stage 

$i = $XL; 
$j = $YL; 

while($i > 0 && $j > 0) {
  printf "\n%4i%4i", $i, $j;   
  if   ($F[$i-1][$j-1] + $s[$i][$j] == $F[$i][$j]) {
    push @upper, $x[$i];
    push @lower, $y[$j];
    $i=$i-1; $j=$j-1;
    }
  elsif($F[$i-1][$j]   + $g[$x[$i]] == $F[$i][$j]) {
    push @upper, $x[$i];
    push @lower, '-';
    $i= $i-1;
    }
  elsif($F[$i][$j-1]   + $g[$y[$j]] == $F[$i][$j]) {
    push @upper, '-';
    push @lower, $y[$j];
    $j =$j-1;
    }
  else {
    die "Error: traceback failed!\n"
    }
  } 

while($i > 0) {
  printf "\n%4i%4i", $i, $j;   
  push @upper, $x[$i];
  push @lower, '-';
  $i=$i-1;
  }
while($j > 0) {
  printf "\n%4i%4i", $i, $j;   
  push @upper, '-';
  push @lower, $y[$j];
  $j=$j-1;
  }

# print optimal alignment:

  print "\n\nOptimal alignment score: $F[$XL][$YL]\n";
  print "\nsequence 1 "; 
  while(scalar(@upper) > 0) {$seq1.=(pop @upper)}
  print "$seq1";
  print "\nsequence 2 ";
  while(scalar(@lower) > 0) {$seq2.=(pop @lower)}  
  print "$seq2";
  print "\n\n";

# print minimum output

while(scalar(@upper) > 0) {$seq1.=(pop @upper)}  
while(scalar(@lower) > 0) {$seq2.=(pop @lower)} 


$seq1=~tr/01234567/ABCDEFGX/;
$seq2=~tr/01234567/ABCDEFGX/;

print "al_seqStart\nal_seq1 $seq1\nal_seq2 $seq2\nal_seqEnd\n";


