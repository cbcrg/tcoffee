#!/usr/bin/env perl

my $msaf="msa.in.tmp.$$";
my $msaoutf="msa.out.tmp.$$";
my $err="msa.out.err.$$";
open  (F, $ARGV[0]);
open  (OUT, ">$msaf");
while (<F>)
  {
    $l=$_;
    if ( $l=~/^>(\S+)/)
      {
	$s=$seqname{$nseq++}=$1;
	print OUT "$l";
	
      }
    else 
      {
	$l=uc($l);
	print OUT "$l";
      }
  }

close (F);
close(OUT);

system ("msa $msaf > $msaoutf 2>$err");
open (F, "$msaoutf");
$read=0;
$cn=0;
while (<F>)
  {
    $l=$_;
    if ($read)
      {
	if ($l=~/End gaps not penalized/){$read=0;}
	elsif (!($l=~/\S/))
	  {
	    $cn=0;
	  }
	else
	  {
	    
	    chomp ($l);
	    $seqal{$cn++}.=$l;
	    $tot++;
	  }
      }
    elsif ($l=~/Optimal Multiple Alignment/)
      {
	$read=1;
      }
  }
close (F);

if ($tot<1)
  {
    print STDERR "\nWarning: MSA returned a NULL file -- Use T-Coffee instead\n";
    open (F,$err);
    while (<F>){print "$_";}
      
    system ("t_coffee -seq $msaf -outfile $ARGV[1]  -quiet");
  }
else
  {
    open (OUT, ">$ARGV[1]");
    for ($a=0; $a<$nseq;$a++)
      {
	print OUT ">$seqname{$a}\n$seqal{$a}\n";
      }
    close (OUT);
  }
unlink ($msaf);
unlink ($msaoutf);
unlink ($err);
