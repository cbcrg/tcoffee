#! /usr/bin/env perl

open (F, $ARGV[0]);
while (<F>)
  {
    $l=$_;
    @list=($l=~/([^;]+)/g);
    

    $name=$list[0];
    $size=$list[1]+$list[2]+$list[3];
    @gname=(@gname, $name);
    @gsize=(@gsize, $size);
    $total_size+=$size;
  }
close (F);
$ngenes=$#gname+1;

print "\nInitialisation Start";

for ( $a=0; $a<$ngenes; $a++)
  {
    for ($b=0; $b<$gsize[$a]; $b++, $c++)
      {
	$genome[$c]=$gname[$a];
      }
  }
print "\nInitialisation Finished";



$mut_level=$ARGV[1];
$threshold=$ARGV[2];
$n_mut=$mut_level*($total_size/1000);
$n=1000;
for ($a=0; $a<$ngenes; $a++)
      {
	$sum{$gname[$a]}=0;
	$sum2{$gname[$a]}=0;
	$extreme{$name}=0;
      }

for ( $sim=0; $sim<$n; $sim++)
  {
    for ($a=0; $a<$n_mut; $a++)
      {
	$mut_table{$genome[int (rand $total_size)]}++;
      }

    for ($a=0; $a<$ngenes; $a++)
      {
	$name=$gname[$a];
	$mutf=1000*$mut_table{$name}/$gsize[$a];
	$freq{$name}=$mutf;
	$max{$name}=($mutf>$max{$name})?$mutf:$max{$name};
	$sum{$name}+=$mutf;
	$sum2{$name}+=$mutf*$mutf;
	$mut_table{$name}=0;
	
	$extreme{$name}+=($mutf>2*$mut_level || $mutf<$mut_level/2)?1:0;
      }
    $ne=0;
    for ( $a=0; $a<$ngenes; $a++)
      {
	$name=$gname[$a];
	$mutf=$freq{$name};
	if (($mutf>2*$mut_level || $mutf<$mut_level/2) && $gsize[$a]<$threshold ){$ne++;}
      }
    $special{$ne}++;
  }

for ( $ne=0,$a=0; $a< $ngenes; $a++)
  {

    if ($gsize[$a]<$threshold ){$ne++}
}
print "\n$ne genes smaller than $threshold\n";
for ( $a=0; $a<10; $a++)
  {
    $special{$a}=($special{$a}*100)/$n;
    print "\n $a: $special{$a}";
  }

for ( $a=0; $a< $ngenes; $a++)
  {
    $sum2=$sum2{$gname[$a]};
    $sum=$sum{$gname[$a]};
    
    $sd=sqrt($sum2*$n-$sum*$sum)/$n;
        
    printf "\n%10s %10d %5.2f MAX=%.3f EXT %d", $gname[$a],$gsize[$a], $sd, $max{$gname[$a]}, $extreme{$gname[$a]};
}
