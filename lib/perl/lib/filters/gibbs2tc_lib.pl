#!/usr/bin/perl


my $PROBTRESH = 0.3;# base pairs below this prob threshold will be ignored
my $WEIGHT = 100.0; # float!!
my $NUCALPH = "ACGTUNRYMKSWHBVD";
use vars qw($NUCALPH $WEIGHT);

my $myname = basename($0);


#use diagnostics;

use File::Basename;
use Getopt::Long;
use File::Glob ':glob';
use File::Spec;
use File::Temp qw/ tempfile tempdir /;




my $fmsq = "";
my $flib = "";
my %OPTS;
my %seq;
my ($id, $nseq, $i);
my @nl;

# ---   option handling
#




# ---   read seqences
#
$fmseq=$ARGV[0];
%seq=read_fasta_seq($fmseq);

`Gibbs $fmseq 20,20,20,20 -o x`;

%tc=&gibbs2aln (x);

@sl=keys (%seq);
$nseq=$#sl+1;
print "! TC_LIB_FORMAT_01\n";
print "$nseq\n";
foreach $s (@sl)
  {
    $l=length ( $seq{$s}{seq});
    print "$seq{$s}{name} $l $seq{$s}{seq}\n";
  }
for ($a=1; $a<=$nseq-1; $a++)
  {
    $s1=$sl[$a-1];
    for ($b=$a+1; $b<=$nseq; $b++)
      {
	$s2=$sl[$b-1];
	$key="$s1\_$s2";
	
	print "# $a $b\n";
	for ($c=0; $c<$tc{$key}{n}; $c++)
	  {
	    print "\t$tc{$key}{$c}\n";
	  }
      }
  }
print "! CPU 0\n! SEQ_1_TO_N\n";

  
sub read_fasta_seq 
  {
    my $f=@_[0];
    my %hseq;
    my (@seq, @com, @name);
    my ($a, $s,$nseq);

    open (F, $f);
    while (<F>)
      {
	$s.=$_;
      }
    close (F);

    @name=($s=~/>(\S*).*\n[^>]*/g);
    
    @seq =($s=~/>.*.*\n([^>]*)/g);
    @com =($s=~/>(\S*)(.*)\n([^>]*)/g);


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

sub gibbs2aln
{
  my $f=@_[0];
  open (F, $f);
  
  
  while (<F>)
    {
      if (/Num Motifs/){$record=1;$n++;}
      elsif ($record && /\d\,\s+\d\s+\w/ )
        {$aln{$n}{$aln{$n}{n}++}=$_;}
      elsif ($record==1)
	{$record=0;}
    }
  
  for ($a=0; $a<$n; $a++)
    {
      $nseq=$aln{$a}{n};
      
      for ($b=0; $b<$nseq-1; $b++)
	{
	  @seq1=($aln{$a}{$b}=~/(\S+)/g);
	  for ($c=0; $c<$nseq; $c++)
	    {
	      @seq2=($aln{$a}{$c}=~/(\S+)/g);
	      $s1=$seq1[9];
	      $s2=$seq2[9];
	      if ( $s1 eq $s2){next;}
	      
	      $key1="$s1\_$s2";
	      $key2="$s2\_$s1";
	     
	      for ($d=$seq1[2], $e=$seq2[2]; $d<=$seq1[6]; $d++, $e++)
		{
		  $r1="$d $e 50";
		  if (!$used{$key1}{$r1})
		    {
		      $tc{$key1}{$tc{$key1}{n}++}=$r1;
		      $used{$key1}{$r1}=1;
		    }
		  $r2="$e $d 50";
		  if (!$used{$key2}{$r2})
		    {
		      $tc{$key2}{$tc{$key2}{n}++}=$r2;
		      $used{$key2}{$r2}=1;
		    }
		  
		}
	    }
	}
    }
  return %tc;
}

