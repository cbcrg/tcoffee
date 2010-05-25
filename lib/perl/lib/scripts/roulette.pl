#!/usr/bin/env perl


$attribute {'0'}="0",

$attribute {'01'}="01_o_r_C1_Q1";
$attribute {'02'}="02_e_b_C2_Q1";
$attribute {'03'}="03_o_r_C3_Q1";

$attribute {'04'}="04_e_b_C1_Q1";
$attribute {'05'}="05_o_r_C2_Q1";
$attribute {'06'}="06_e_b_C3_Q1";

$attribute {'07'}="07_o_r_C1_Q1";
$attribute {'08'}="08_e_b_C2_Q1";
$attribute {'09'}="09_o_r_C3_Q1";

$attribute {'10'}="10_e_b_C1_Q1";
$attribute {'11'}="11_o_b_C2_Q1";
$attribute {'12'}="12_e_r_C3_Q1";

$attribute {'13'}="13_o_b_C1_Q2";
$attribute {'14'}="14_e_r_C2_Q2";
$attribute {'15'}="15_o_b_C3_Q2";

$attribute {'16'}="16_e_r_C1_Q2";
$attribute {'17'}="17_o_b_C2_Q2";
$attribute {'18'}="18_e_r_C3_Q2";

$attribute {'19'}="19_o_r_C1_Q2";
$attribute {'20'}="20_e_b_C2_Q2";
$attribute {'21'}="21_o_r_C3_Q2";

$attribute {'22'}="22_e_b_C1_Q2";
$attribute {'23'}="23_o_r_C2_Q2";
$attribute {'24'}="24_e_b_C3_Q2";

$attribute {'25'}="25_o_r_C1_Q3";
$attribute {'26'}="26_e_b_C2_Q3";
$attribute {'27'}="27_o_r_C3_Q3";

$attribute {'28'}="28_e_r_C1_Q3";
$attribute {'29'}="29_o_b_C2_Q3";
$attribute {'30'}="30_e_r_C3_Q3";

$attribute {'31'}="31_o_b_C1_Q3";
$attribute {'32'}="32_e_r_C2_Q3";
$attribute {'33'}="33_o_b_C3_Q3";

$attribute {'34'}="34_e_r_C1_Q3";
$attribute {'35'}="35_o_b_C2_Q3";
$attribute {'36'}="36_e_r_C3_Q3";


$list{"00"}{"name"}="00";
$list{"00"}{"gain"}=35;

$list{"01"}{"name"}="01";
$list{"01"}{"gain"}=35;

$list{"02"}{"name"}="02";
$list{"02"}{"gain"}=35;

$list{"03"}{"name"}="03";
$list{"03"}{"gain"}=35;

$list{"04"}{"name"}="04";
$list{"04"}{"gain"}=35;


$list{"05"}{"name"}="05";
$list{"05"}{"gain"}=35;

$list{"06"}{"name"}="06";
$list{"06"}{"gain"}=35;

$list{"07"}{"name"}="07";
$list{"07"}{"gain"}=35;

$list{"08"}{"name"}="08";
$list{"08"}{"gain"}=35;

$list{"09"}{"name"}="09";
$list{"09"}{"gain"}=35;

$list{"10"}{"name"}="10";
$list{"10"}{"gain"}=35;

$list{"11"}{"name"}="11";
$list{"11"}{"gain"}=35;

$list{"12"}{"name"}="12";
$list{"12"}{"gain"}=35;

$list{"13"}{"name"}="13";
$list{"13"}{"gain"}=35;

$list{"14"}{"name"}="14";
$list{"14"}{"gain"}=35;

$list{"15"}{"name"}="15";
$list{"15"}{"gain"}=35;

$list{"16"}{"name"}="16";
$list{"16"}{"gain"}=35;

$list{"17"}{"name"}="17";
$list{"17"}{"gain"}=35;

$list{"18"}{"name"}="18";
$list{"18"}{"gain"}=35;


$list{"19"}{"name"}="19";
$list{"19"}{"gain"}=35;

$list{"20"}{"name"}="20";
$list{"20"}{"gain"}=35;


$list{"21"}{"name"}="21";
$list{"21"}{"gain"}=35;

$list{"22"}{"name"}="22";
$list{"22"}{"gain"}=35;

$list{"23"}{"name"}="23";
$list{"23"}{"gain"}=35;

$list{"24"}{"name"}="24";
$list{"24"}{"gain"}=35;

$list{"25"}{"name"}="25";
$list{"25"}{"gain"}=35;

$list{"26"}{"name"}="26";
$list{"26"}{"gain"}=35;

$list{"27"}{"name"}="27";
$list{"27"}{"gain"}=35;

$list{"28"}{"name"}="28";
$list{"28"}{"gain"}=35;

$list{"29"}{"name"}="29";
$list{"29"}{"gain"}=35;

$list{"30"}{"name"}="30";
$list{"30"}{"gain"}=35;

$list{"31"}{"name"}="31";
$list{"31"}{"gain"}=35;

$list{"32"}{"name"}="32";
$list{"32"}{"gain"}=35;


$list{"33"}{"name"}="33";
$list{"33"}{"gain"}=35;

$list{"34"}{"name"}="34";
$list{"34"}{"gain"}=35;

$list{"35"}{"name"}="35";
$list{"35"}{"gain"}=35;

$list{"36"}{"name"}="36";
$list{"36"}{"gain"}=35;

$list{"e"}{"name"}="e";
$list{"e"}{"gain"}=2;
$list{"o"}{"name"}="o";
$list{"o"}{"gain"}=2;

$list{"b"}{"name"}="b";
$list{"b"}{"gain"}=2;
$list{"r"}{"name"}="r";
$list{"r"}{"gain"}=2;

$list{"C1"}{"name"}="C1";
$list{"C1"}{"gain"}=3;

$list{"C2"}{"name"}="C2";
$list{"C2"}{"gain"}=3;

$list{"C3"}{"name"}="C3";
$list{"C3"}{"gain"}=3;

$list{"Q1"}{"name"}="Q1";
$list{"Q1"}{"gain"}=3;

$list{"Q2"}{"name"}="Q2";
$list{"Q2"}{"gain"}=3;

$list{"Q3"}{"name"}="C1";
$list{"Q3"}{"gain"}=3;



if ( $ARGV[0] eq "-bet")
  {
    $x=int rand(36);

    @bet=@ARGV;
    $v="SSS";
    $v="0".$x;
    if ( $x<=9){$v="0".$x;}
    else {sprintf $v=$x;}
    
    print "\t\t$v ";

    for ( $a=1; $a<=$#bet; $a+=2)
      {
	$b=$bet[$a];
	$vb=$bet[$a+1];

	if ($attribute{"$v"}=~/$b/)
	  {
	    $vb=$list{"$b"}{"gain"}*$vb;
	  }
	else
	  {
	    $vb=0;
	  }
	$tot_bet+=$bet[$a+1];
	$tot_gain+=$vb;
	$tot=$tot."$b $vb "; 
      }
    $delta=$tot_gain- $tot_bet;
    print "$tot_bet $tot_gain $tot\n";
  }
elsif ( $ARGV[0]="-test")
  {
    $start=$ARGV[1];
    $pr=1;
    for ( $a=0; $a<1000; $a++)
      {
	$money=$start;
	$n=0;
	$cbet=1;
	$pr=1;
	$inc=$ARGV[2];
	while ($money>0 && $money<$start+$inc)
	  {
	    
#	    if ( $gain>0){$cbet++;}
#	    else { $cbet=1;}
	    	  

	    #if ( $pr==1){$l=`roulette.pl -bet Q2 $cbet Q3 $cbet`;}
	    #elsif ( $pr==2){$l=`roulette.pl -bet Q1 $cbet Q3 $cbet`;}
	    #elsif ( $pr==3){$l=`roulette.pl -bet Q1 $cbet Q2 $cbet`;}
	    $l=`roulette.pl -bet Q1 $cbet Q2 $cbet C1 $cbet C2 $cbet o $cbet r $cbet`;
	    @r=($l=~/(\S+)/g);
	    
	    if ( $attribute{'$r[0]'}=~/Q1/){$pr=1;}
	    elsif ( $attribute{'$r[0]'}=~/Q2/){$pr=2;}
	    elsif ( $attribute{'$r[0]'}=~/Q3/){$pr=3;}
	    	    

	    $bet=$r[1];
	    $gain=$r[2];
	    
	    if ( $gain){$loose_strike=0;$win_strike++;}
	    elsif ( $gain==0){$loose_strike++;$win_strike=0;}
	    
	    $money=$money-$bet+$gain;
	    $n++;
	    
	  }
	$earnt+=$money;
	$spent+=$start;
	
	if ( $money>$start){$win++; $ngen_win+=$n}
	else
	  {
	    $loose++;
	  }
	$delta=$earnt-$spent;
	$p=($win/($win+$loose));
	print "$a $n $money Delta $delta P $p\n";
	if ( $money <=0){$a=100000;}
	$start++;
	$inc++;
      }
    $p=($win/($win+$loose));
    $n=$ngen_win/$win;
    print "$p $n $delta";
  }
