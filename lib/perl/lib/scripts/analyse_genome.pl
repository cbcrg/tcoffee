#! /usr/bin/perl


if ( -e $ARGV[0])
  {
    %h=&read_stat ($ARGV[0]);
  }
else
  {
    %h=&make_stat ($ARGV[0]);
  }

&display_count (%h);
&count2analyze (%h);

sub display_count 
  {
     my %h=(@_);
     
     @list1=keys (%{$h{1}});
     @list2=keys (%{$h{2}});

     foreach $m1 (@list2)
       {
	 printf "COUNT $m1 %7d %7d\n", $h{2}{$m1}{$m1}{'count_1'}, $h{2}{$m1}{$m1}{'count_2'};
       }
   }

     
sub count2analyze 
  {
    my %h=(@_);
    
    @list1=keys (%{$h{1}});
    @list2=keys (%{$h{2}});
    


    for ($a=1; $a<=2; $a++)
      {
	@list=keys (%{$h{$a}});
	foreach $m1 (@list)
	  {
	    $s{$a}{sum_1}+=$h{$a}{$m1}{$m1}{'count_1'};
	    $s{$a}{sum_2}+=$h{$a}{$m1}{$m1}{'count_2'};
	    foreach $m2 (@list)
	      {
		$s{$a}{sum}+=$h{$a}{$m1}{$m2}{'sub'};
	      }
	  }
	print "Total Count K=$a Seq1: $s{$a}{sum_1} Seq2: $s{$a}{sum_2} Tot: $s{$a}{sum}\n";
      }

    foreach $m1 (@list2)
      {
	foreach $m2 ( @list2)
	  {
	    $observed1=$h{2}{$m1}{$m2}{'sub'}/$s{2}{sum};
	    $expected1=($h{2}{$m1}{$m1}{'count_1'}/$s{2}{sum_1})*($h{2}{$m2}{$m2}{'count_2'}/$s{2}{sum_2});
	    $ratio1=$result{$m1}{$m2}{'ratio'}=$observed1/$expected1;
	    $logodd1=$result{$m1}{$m2}{'logodd'}=&mlog($ratio1);

	    @l1=split (//, $m1);
	    @l2=split (//, $m2);


	    $observed2=$h{1}{$l1[0]}{$l2[0]}{'sub'}/$s{1}{sum};
	    $observed3=$h{1}{$l1[1]}{$l2[1]}{'sub'}/$s{1}{sum};
	    $expected2=($h{1}{$l1[0]}{$l1[0]}{'count_1'}/$s{1}{sum_1})*($h{1}{$l2[0]}{$l2[0]}{'count_2'}/$s{1}{sum_2});
	    $expected3=($h{1}{$l1[1]}{$l1[1]}{'count_1'}/$s{1}{sum_1})*($h{1}{$l2[1]}{$l2[1]}{'count_2'}/$s{1}{sum_2});

	    $ratio2=$observed2/$expected2;
	    $ratio3=$observed3/$expected3;
	    $logodd2=&mlog($ratio2);
	    $logodd3=&mlog($ratio3);
	    
	    $diff1=$result{$m1}{$m2}{'diff'}=$logodd1-($logodd2+$logodd3);
	    
	    printf "#$m1 $m2 %5.2f %5.2f\n", $logodd1, $diff1;
	    $result{$m1}{$m1}{P}+=($diff1>0)?1:0;
	    $result{$m1}{$m1}{N}+=($diff1<0)?1:0;
	    $result{$m1}{$m1}{T}+=$diff1;
	  }
      }
 
    foreach $m1 (@list2)
      {
	printf "$m1 %2d %2d %5.2f\n",$result{$m1}{$m1}{P}, $result{$m1}{$m1}{N},$result{$m1}{$m1}{T};
      }
  }

sub read_stat 
  {
    my $f=@_[0];
    my ($l, @l,%h,$m1, $m2);
    open (F, "$f");
    

    while (<F>)
      {
	$l=$_;

	if ($l=~/\#/ || !($l=~/\S/)){;}
	else
	  {
	    @l=split /\s+/, $l;
	    if ($l=~/\//)
	      {
		($m1,$m2)=split /\//, $l[0];
		$h{length($m1)}{$m1}{$m2}{'sub'}=$l[1];
		#print "$m1 $m2 $h{length($m1)}{$m1}{$m2}{sub}\n"
	      }
	    else
	      {
		$m1=$l[0];
		$h{length($m1)}{$m1}{$m1}{'count_1'}=$l[1];
		$h{length($m1)}{$m1}{$m1}{'count_2'}=$l[2];
		#print "$l[0] $l[1] $l[2]\n";
	      }
	  }
      }
    return %h;
  }



# Create a genome and return the counts

sub make_stat 
  {
    my $l=@_[0];
    return &genome2stat (&make_genome2($l));
  }

sub genome2stat
  {
    my @G=@_;
    my ($a, $n1, $n2, %h);
    
    @g1=split (//,$G[0]);
    @g2=split (//,$G[1]);
    
  
    for ($a=1; $a<=$#g1; $a++)
      {
	$n1=$g1[$a];
	$n2=$g2[$a];
	
	$h{1}{$n1}{$n1}{'count_1'}++;
	$h{1}{$n2}{$n2}{'count_2'}++;
	$h{1}{$n1}{$n2}{'sub'}++;
	
	$n1=$g1[$a].$g1[$a-1];
	$n2=$g2[$a].$g2[$a-1];
	
	
	
	$h{2}{$n1}{$n1}{'count_1'}++;
	$h{2}{$n2}{$n2}{'count_2'}++;
	$h{2}{$n1}{$n2}{'sub'}++;
      }
    return %h;
  }
	
	
sub make_genome1
  {
    my $l=@_[0];
    my ($n, $a, @g, %comp, @list);
    my $id=50;
    #Make two random sequences with composition bias having %id %identity
    
    $comp{A}=10;
    $comp{G}=5;
    $comp{C}=5;
    $comp{T}=10;

    
    foreach $n (keys (%comp))
      {
	for ( $a=0; $a<$comp{$n}; $a++)
	  {
	    push @list, $n;
	  }
      }
    for ($a=0; $a<$l; $a++)
      {
	$n1=$list[rand($#list)];
	$n2=$list[rand($#list)];
	
	$g[0].=$n1;
	$g[1].=(int(rand(100))<$id)?$n1:$n2;
      }
    print "\n$g[0]\n$g[1]";
    return @g;
  }

sub make_genome2
  {
    my $l=@_[0];
    my ($n, $a, @g, %comp, %list, @G, $n1, $n2, $m1, $m2, %list2, @l2);
    my $id=50;
    #Make two random sequences with composition bias
    $A=100;
    $G=50;
    $C=50;
    $T=100;

    $comp{A}{A}=$A;
    $comp{A}{G}=$G;
    $comp{A}{C}=$C;
    $comp{A}{T}=$T;

    $comp{G}{A}=$A;
    $comp{G}{G}=$G;
    $comp{G}{C}=$C;
    $comp{G}{T}=$T;

    $comp{C}{A}=$A;
    $comp{C}{G}=$G/10;
    $comp{C}{C}=$C;
    $comp{C}{T}=$T;

    $comp{T}{A}=$A;
    $comp{T}{G}=$G;
    $comp{T}{C}=$C;
    $comp{T}{T}=$T;

    foreach $n1 (keys (%comp))
      {
	foreach $n2 (keys (%{$comp{$n1}}))
	  {
	    for ( $a=0; $a<$comp{$n1}{$n2}; $a++)
	      {
		$list{$n1}{$list{$n1}{"len"}++}=$n2;
	      }
	    $m1=$n1.$n2;
	    @l2=(@l2,$m1);
	  }
      }
    

    
    foreach $m1 (@l2)
      {
	foreach $m2 (@l2)
	  {
	    if ($m1 eq "CG"){$v=1000;}
	    else
	      {
		$v=($list2{$m1}{$m2})?$list2{$m1}{$m2}:100;
	      }
	    if ($m2 eq $m1){$v*=10;}
	    for ($a=0; $a<$v; $a++)
	      {
		$list2{$m1}{$list2{$m1}{"len"}++}=$m2;
	      }
	  }
      }

    $np1=$np2="A";
#Build first sequence    
    for ($a=0; $a<$l; $a++)
      {
	$np1=$n1=$list{$np1}{int(rand($list{$np1}{len}))};	
	
	$g[0].=$n1;
      }

#build second sequence
    @G=split (//, $g[0]);
    for ($a=0; $a<$l-1; $a++)
      {
	$m1=$G[$a].$G[$a+1];
	$m2=$list2{$m1}{int(rand($list2{$m1}{len}))};
	
	($n1,$n2)=split (//,$m2);
	$g[1].=$n1;
      }
    $g[1].=$n2;
    
    #print "\n$g[0]\n$g[1]\n";die;
    return @g;
  }
    
    
sub mlog
  {
    my $l=@_[0];

    return ($l<=0)?-9999:log($l);
  }
