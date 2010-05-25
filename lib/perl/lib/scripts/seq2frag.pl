#! /usr/bin/perl


$len=100;
$frag=100;
$cov=100;


$n=($len/$frag)*10+1;
$n=2;

for ( $c=0; $c<$cov; $c++)
  {
    for ($a=1; $a<=$n; $a++)
      {
	$i[$a]=rand($len);
      }
    $j[$n]=$len;
    @j=sort {$a<=>$b} @i;
    for ($a=1; $a<=$n; $a++)
      {
	@l[$j[$a]-$j[$a-1]]++;
      }
  }


for ($a=1; $a<200; $a++)
  {
    print "\n$a $l[$a]";
  }
