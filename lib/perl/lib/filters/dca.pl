#!/usr/bin/env perl

my $msaf="msa.in.tmp.$$";
my $msaoutf="msa.out.tmp.$$";
my $cost="blosum62.tmp.$$";

open  (F, $ARGV[0]);
open  (OUT, ">$msaf");
$nseq=0;
while (<F>)
  {
    $l=$_;
    if ( $l=~/^>(\S+)/)
      {
	my $simple="Seq$nseq";
	$s=$seqname{$nseq++}=$1;
	$translate{$simple}=$s;
	
	print OUT ">$simple\n";
	
      }
    else
      {
	$l=uc($l);
	print OUT "$l";
      }
  }
close (F);
close(OUT);

dump_blosum ($cost);
system ("dca -c $cost -q $msaf> $msaoutf 2>/dev/null");
open (F, "$msaoutf");
open (OUT, ">$ARGV[1]");

$read=0;
while (<F>)
  {
    $l=$_;
    if ($l=~/^>(\S+)/)
      {
	$read=1;
	$name=$translate{$1};
	print OUT ">$name\n";
      }
    elsif ($read && ($l=~/\S/))
      {
	print OUT "$l";
      }
    else
      {
	$read=0;
      }
  }
close (F);

unlink ($cost);
unlink ($msaf);
unlink ($msaoutf);

sub dump_blosum
  {
    my $f=shift;
    open (F, ">$f");

    print F "6\n";
    print F "- -   0\n";
    print F "W W   0\n";
    print F "Y Y   4\n";
    print F "F F   5\n";
    print F "V V   7\n";
    print F "L L   7\n";
    print F "I I   7\n";
    print F "M M   6\n";
    print F "K K   6\n";
print F "R R   6\n";
    print F "H H   3\n";
    print F "Q Q   6\n";
    print F "E E   6\n";
    print F "D D   5\n";
    print F "N N   5\n";
    print F "G G   5\n";
    print F "A A   7\n";
    print F "P P   4\n";
    print F "T T   6\n";
    print F "S S   7\n";
    print F "C C   2\n";
    print F "- C  10 \n";
    print F "- S  10\n";
    print F "- T  10 \n";
    print F "- P  10\n";
    print F "- A  10 \n";
    print F "- G  10\n";
    print F "- N  10 \n";
    print F "- D  10\n";
    print F "- E  10 \n";
    print F "- Q  10\n";
print F "- H  10 \n";
    print F "- R  10\n";
    print F "- K  10 \n";
    print F "- M  10\n";
    print F "- I  10 \n";
    print F "- L  10\n";
    print F "- V  10 \n";
    print F "- F  10\n";
    print F "- Y  10 \n";
    print F "- W  10\n";
    print F "W C  13 \n";
    print F "W S  14\n";
    print F "W T  13 \n";
    print F "W P  15\n";
    print F "W A  14 \n";
    print F "W G  13\n";
    print F "W N  15 \n";
    print F "W D  15\n";
    print F "W E  14 \n";
    print F "W Q  13\n";
    print F "W H  13 \n";
    print F "W R  14\n";
    print F "W K  14 \n";
    print F "W M  12\n";
    print F "W I  14 \n";
    print F "W L  13\n";
    print F "W V  14 \n";
    print F "W F  10\n";
    print F "W Y   9 \n";
    print F "Y C  13\n";
    print F "Y S  13 \n";
    print F "Y T  13\n";
    print F "Y P  14 \n";
    print F "Y A  13\n";
    print F "Y G  14 \n";
    print F "Y N  13\n";
    print F "Y D  14 \n";
    print F "Y E  13\n";
    print F "Y Q  12 \n";
    print F "Y H   9\n";
    print F "Y R  13 \n";
    print F "Y K  13\n";
    print F "Y M  12 \n";
    print F "Y I  12\n";
    print F "Y L  12 \n";
    print F "Y V  12\n";
    print F "Y F   8 \n";
    print F "F C  13\n";
print F "F S  13 \n";
    print F "F T  13\n";
    print F "F P  15 \n";
    print F "F A  13\n";
    print F "F G  14 \n";
    print F "F N  14\n";
    print F "F D  14 \n";
    print F "F E  14\n";
    print F "F Q  14 \n";
    print F "F H  12\n";
    print F "F R  14 \n";
    print F "F K  14\n";
    print F "F M  11 \n";
    print F "F I  11\n";
    print F "F L  11 \n";
    print F "F V  12\n";
    print F "V C  12 \n";
    print F "V S  13\n";
    print F "V T  11 \n";
    print F "V P  13\n";
    print F "V A  11 \n";
    print F "V G  14\n";
    print F "V N  14 \n";
    print F "V D  14\n";
print F "V E  13 \n";
print F "V Q  13\n";
print F "V H  14 \n";
print F "V R  14\n";
print F "V K  13 \n";
print F "V M  10\n";
print F "V I   8 \n";
print F "V L  10\n";
print F "L C  12 \n";
print F "L S  13\n";
print F "L T  12 \n";
print F "L P  14\n";
print F "L A  12 \n";
print F "L G  15\n";
print F "L N  14 \n";
print F "L D  15\n";
print F "L E  14 \n";
print F "L Q  13\n";
print F "L H  14 \n";
print F "L R  13\n";
print F "L K  13 \n";
print F "L M   9\n";
print F "L I   9 \n";
print F "I C  12\n";
print F "I S  13 \n";
print F "I T  12\n";
print F "I P  14 \n";
print F "I A  12\n";
print F "I G  15 \n";
print F "I N  14\n";
print F "I D  14 \n";
print F "I E  14\n";
print F "I Q  14 \n";
print F "I H  14\n";
print F "I R  14 \n";
print F "I K  14\n";
print F "I M  10 \n";
print F "M C  12\n";
print F "M S  12 \n";
print F "M T  12\n";
print F "M P  13 \n";
print F "M A  12\n";
print F "M G  14 \n";
print F "M N  13\n";
print F "M D  14 \n";
print F "M E  13\n";
print F "M Q  11 \n";
print F "M H  13\n";
print F "M R  12 \n";
print F "M K  12\n";
print F "K C  14 \n";
print F "K S  11\n";
print F "K T  12 \n";
print F "K P  12\n";
print F "K A  12 \n";
print F "K G  13\n";
print F "K N  11 \n";
print F "K D  12\n";
print F "K E  10 \n";
print F "K Q  10\n";
print F "K H  12 \n";
print F "K R   9\n";
print F "R C  14 \n";
print F "R S  12\n";
print F "R T  12 \n";
print F "R P  13\n";
print F "R A  12 \n";
print F "R G  13\n";
print F "R N  11 \n";
print F "R D  13\n";
print F "R E  11 \n";
print F "R Q  10\n";
print F "R H  11 \n";
print F "H C  14\n";
print F "H S  12 \n";
print F "H T  13\n";
print F "H P  13 \n";
print F "H A  13\n";
print F "H G  13 \n";
print F "H N  10\n";
print F "H D  12 \n";
print F "H E  11\n";
print F "H Q  11 \n";
print F "Q C  14\n";
print F "Q S  11 \n";
print F "Q T  12\n";
print F "Q P  12 \n";
print F "Q A  12\n";
print F "Q G  13 \n";
print F "Q N  11\n";
print F "Q D  11 \n";
print F "Q E   9\n";
print F "E C  15 \n";
print F "E S  11\n";
print F "E T  12 \n";
print F "E P  12\n";
print F "E A  12 \n";
print F "E G  13\n";
print F "E N  11 \n";
print F "E D   9\n";
print F "D C  14 \n";
print F "D S  11\n";
print F "D T  12 \n";
print F "D P  12\n";
print F "D A  13 \n";
print F "D G  12\n";
print F "D N  10 \n";
print F "N C  14\n";
print F "N S  10 \n";
print F "N T  11\n";
print F "N P  13 \n";
print F "N A  13\n";
print F "N G  11 \n";
print F "G C  14\n";
print F "G S  11 \n";
print F "G T  13\n";
print F "G P  13 \n";
print F "G A  11\n";
print F "A C  11 \n";
print F "A S  10\n";
print F "A T  11 \n";
print F "A P  12\n";
print F "P C  14 \n";
print F "P S  12\n";
print F "P T  12 \n";
print F "T C  12\n";
print F "T S  10 \n";
print F "S C  12\n";
close (F);
    return;
  }
    
