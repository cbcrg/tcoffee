#!/usr/bin/env perl
$max=5;

open (F, $ARGV[0]);

while (<F>)
  {$tmp.=$_;}


@seq=($tmp=~/(>[^>]+)/g);


if ( $#seq <=$max)
  {
    print $tmp;
    exit (EXIT_SUCCESS);
  }

for ($a=0; $a<=$#seq; $a++)
  {
    if ( ($seq[$a]=~/\d{4}/)!=1)
	 {
	   $list[$a]=1; $n++;
	 }
  }



while ($n<$max)
  {
    $x=int (rand ($#seq));
    if ($list[$x]==1){;}
    else 
      {$list[$x]=1; $n++;}
  }

for ($a=0; $a<=$#seq; $a++)
{
  if ( $list[$a]==1){print "$seq[$a]";}
}
