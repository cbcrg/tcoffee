#!/usr/bin/env perl





$cmd=join ' ', @ARGV;
if ($cmd=~/-infile=(\S+)/){ $seqfile=$1;}
if ($cmd=~/-outfile=(\S+)/){ $libfile=$1;}



%s=read_fasta_seq ($seqfile);

open (F, ">$libfile");
foreach $name (keys (%s))
  {
    my $tclib="$name.RNAplfold_tclib";
    print (F ">$name _F_ $tclib\n");
    seq2RNAplfold2tclib ($name, $s{$name}{seq}, $tclib);
  }
close (F);
exit (EXIT_SUCCESS);

sub seq2RNAplfold2tclib
  {
    my ($name, $seq, $tclib)=@_;
    my ($tmp);
    $n++;
    $tmp="tmp4seq2RNAplfold_tclib.$$.$n.pep";
    open (RF, ">$tmp");
    print (RF ">$name\n$seq\n");
    close (RF);
    
    system "t_coffee -other_pg RNAplfold2tclib.pl -in=$tmp -out=$tclib";
    
    unlink ($tmp);
    return $tclib;
  }
    
    
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
    @com =($s=~/>\S*(.*)\n([^>]*)/g);

    
    $nseq=$#name+1;
    
    for ($a=0; $a<$nseq; $a++)
      {
	my $n=$name[$a];
	$hseq{$n}{name}=$n;
	$hseq{$n}{seq}=$seq[$a];
	$hseq{$n}{com}=$com[$a];
      }
    return %hseq;
  }

