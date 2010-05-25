#!/usr/bin/env perl
use FileHandle;


$N=$ARGV[0];
$N11=$ARGV[1];
$N1msa=$ARGV[2];
$N1sar=$ARGV[3];



print "N=$N N11=$N11 msa=$N1msa sar=$N1sar\n";
$P=&gate_keeper_evalue ($N11, $N1msa, $N1sar, $N);


print " P=$P\n";

sub gate_keeper_evalue
  {
    my $n11=@_[0];
    my $n1msa=@_[1];
    my $n1sar=@_[2];
    my $N=@_[3];
    my $p;

    #1 Number of cases that yield EXACTLY $n11 columns:
    $f=(($n1msa-$n11)>($n1sar-$n11))?$n1msa-$n11:$n1sar-$n11;
    
    $p1=&KinN($n11, $N);
    $p2=&KinN($n1msa-$n11, $N-$n11);
    $p3=&KinN($n1sar-$n11, $N-$n11);
        
    $p_top=&KinN($n1msa, $N)*&KinN($n1sar-$n11, $N-$n1msa)*&KinN($n11, $n1msa);
    printf STDERR "2 in 8: %d %d", &KinN($n1msa-$n11, $N-$n11),&KinN(2, 8);
    print STDERR "\nP1=$p1 P2=$p2 P3=$p3 P=$p_top\n";

    #2 Divide with the number of possibilities
    $p4=&KinN($n1msa, $N);
    $p5=&KinN($n1sar, $N);
    $p_bot=$p4*$p5;

    print STDERR "\nP4=$p4 P5=$p5 P=$p_bot\n";
    $p= -log($p_top/$p_bot);
    
    return $p;
  }

sub KinN
  {
    
    my $K=@_[0];
    my $N=@_[1];
    my $bottom=$N-$K;
    my $R;
    
    
    
    if ( !$fac_lu{1} {$bottom}){$fac_lu{1}{$bottom}=&factorial(1,$bottom);}
    if ( !$fac_lu{$K+1}{$N}){$fac_lu{$K+1}{$N}=&factorial($K+1,$N);}

    $R=$fac_lu{$K+1}{$N}/$fac_lu{1}{$bottom};
    return $R;
  }
    
sub factorial
  {
    my $from=@_[0];
    my $to=@_[1];
    
    my $result=1;
    my $a;
    
    if ( $from>$to){return 1;}
    
    for ($a=$from; $a<=$to; $a++)
      {
        $result*=$a;
      }
    return $result;
  }
sub fac
  {
    my $N=@_[0];
    my $result=1;
    my $a;
    
    for ($a=1; $a<=$N; $a++)
      {
        $result*=$a;
      }
    return $result;
  }

