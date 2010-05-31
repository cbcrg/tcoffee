#!/usr/bin/env perl
use Getopt::Long;
GetOptions("-in=s" => \$fmsq1, "-out=s" => \$outfile, "-arch=s" => \$arch,"-psv=s" => \$psv, "-hmmtop_home=s", \$hmmtop_home );
open (O, ">$outfile");

if (!$hmmtop_home){$hmmtop_home="/home/notredame/packages/hmmtop/hmmtop_2.1";}
if ($arch){$ENV{'HMMTOP_ARCH'}=$arch;}
else {$ENV{'HMMTOP_ARCH'}="$hmmtop_home/hmmtop.arch";}

if ($psv){$ENV{'HMMTOP_PSV'}=$psv;}
else{$ENV{'HMMTOP_PSV'}="$hmmtop_home/hmmtop.psv";}

$fmsq="seq2convert.$$.tmp";
system ("t_coffee -other_pg seq_reformat -in $fmsq1 -output fasta_seq > $fmsq");
%seq=read_fasta_seq($fmsq);

$tmpfile="fasta_seq2hmmtop_fasta.$$.tmp";
foreach $s (keys (%seq))
  {
    
    open F, ">$tmpfile";
    print F ">seq\n$seq{$s}{seq}\n";
    close F;

    $result=`hmmtop -if=$tmpfile -sf=FAS -pl 2>/dev/null`;
    @r=($result=~/(.+)/g);
    foreach $l (@r)
      {
	
	if ($l=~/pred(.*)/)
	  {$p.=$1;}
      }
    
    $p=~s/\s//g;
    print O ">$seq{$s}{name}\n$p\n";
    $p="";
  }
unlink "$tmpfile";
unlink "$fmsq";
close (O);

sub read_fasta_seq 
  {
    my $f=$_[0];
    my %hseq;
    my (@seq, @com, @name);
    my ($a, $s,$nseq);

    open (F, $f);
    while (<F>)
      {
	$s.=$_;
      }
    close (F);

    
    @name=($s=~/>(.*).*\n[^>]*/g);
    
    @seq =($s=~/>.*.*\n([^>]*)/g);
    @com =($s=~/>.*(.*)\n([^>]*)/g);


    $nseq=$#name+1;
    
  
    for ($a=0; $a<$nseq; $a++)
      {
	my $n=$name[$a];
	my $s;
	$hseq{$n}{name}=$n;
	$s=$seq[$a];$s=~s/\s//g;
	
	$hseq{$n}{seq}=$s;
	$hseq{$n}{com}=$com[$a];
      }
    return %hseq;
  }

