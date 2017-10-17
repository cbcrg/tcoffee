#!/usr/bin/env perl
use strict;
use Cwd;
use File::Basename;
my $test=0;

my $tmpdir="/tmp/tco/aligners/upp/";
mymkdir ($tmpdir);



if ($ARGV[0] eq "one")
  {
    seq2msa ($ARGV[1], $ARGV[2]);
  }
elsif ($ARGV[0] eq "all")
  {
    listseq2listmsa ($ARGV[1]);
  }

sub listseq2listmsa
  {
    my $list=shift;
    my $cdir = getcwd;
    my $dir=random_string(10);
    $dir="$tmpdir/$dir/";
    my %h;
    my $n;
    mkdir  ($dir);

    open (F, $list);
    while (<F>)
      {
        my $l=$_;

        chomp($l);
        my @f=split (/\s+/, $l);
	if ( -e $f[0])
          {
            $h{$n}{in}=$f[0];
            ($h{$n}{name},$h{$n}{path})=fileparse ($f[0]);
            $h{$n}{NFin}= "$dir/$h{$n}{name}.seq";
	    
            $h{$n}{NFout}="$dir/$h{$n}{name}.aln";

            $h{$n}{out}=$f[1];

            fasta2fastaupp ($h{$n}{in}, $h{$n}{NFin});
            $n++;
          }
      }
    close (F);
    chdir ($dir);
    
    if (!$test)
      {
	system ("fbname=\$(basename `ls *.seq` .seq); \
             run_upp.py -s \${fbname}.seq -m amino --cpu 1 -d outdir -o \${fbname}.aln; \
             mv outdir/\${fbname}.aln_alignment.fasta \${fbname}.aln;");
      }
    
    foreach my $n (keys (%h))
      {
	if ($test)
	  {
	    system ("cp $h{$n}{NFin} $h{$n}{NFout}");
	    print "$h{$n}{NFin} $h{$n}{NFout} $h{$n}{out}\n";
	  }
        fastaupp2fasta ($h{$n}{NFout},$h{$n}{out});
      }
    chdir ($cdir);
  }

sub seq2msa
    {
      my ($in, $out)=@_;
      my $cdir=getcwd;
      
      
      if (!($in=~/\//)){$in=$cdir."/".$in;}
      if (!($out=~/\//)){$out=$cdir."/".$out;}
      
      my $file=random_string(10);
      $file="$tmpdir/$file";
      open (F, ">$file");
      print F "$in $out\n";
      close (F);
      
      return listseq2listmsa ($file);
    }
	
sub fasta2fastaupp
  {
    my ($in, $out)=@_;
    my ($name, $seq, $n);
    
    if (!-e $in){return;}
    
    open (IN, "$in");
    open (OUT, ">$out");
    local $/ = "\n>";  # read by FASTA record
    
    while (<IN>)
      {
	my $l=$_;
	$l=~s/>//g;
	$l=">".$l;
	
	$l=~/^>(.*)/;
	$name=$1;
	
	$l=~s/^>*.+\n//;
	$l=~s/\n//g;
	$seq=$l;
	
	$seq=~s/u/x/g;
	$seq=~s/U/X/g;
	print OUT ">$name\n$seq\n";
	$n++;
      }
    if ($n==2)
      {
	print OUT ">fake_seq4upp\n$seq\n";
      }
    close (IN);
    close (OUT);
    local $/="\n";
  }

sub fastaupp2fasta
  {
    my ($in, $out)=@_;
    my ($name, $seq, $n);
    
    if (!-e $in){return;}
    
    open (IN, "$in");
    open (OUT, ">$out");
    local $/ = "\n>";  # read by FASTA record
    
    while (<IN>)
      {
	my $l=$_;
	$l=~s/>//g;
	$l=">".$l;
	
	$l=~/^>(.*)/;
	$name=$1;
	
	$l=~s/^>*.+\n//;
	$l=~s/\n//g;
	$seq=$l;
	
	$seq=~s/x/u/g;
	$seq=~s/X/U/g;
	
	if (!($name=~/fake_seq4upp/))
	  {
	    print OUT ">$name\n$seq\n";
	  }
      }
    close (IN);
    close (OUT);
    local $/="\n";
  }

sub random_string
    {
      my $len=shift;
      my @chars = ("A".."Z", "a".."z");
      my $string;
      $string .= $chars[rand @chars] for 1..$len;
      return $string;
    }

sub mymkdir
      {
	my $d=shift;
	my $cd='/';
	
	foreach my $e (split (/\//, $d))
	  {
	    $cd.="$e/";
	    if ( !-d $cd){mkdir ($cd);}
	  }
	return;
      }
      
			  
      
