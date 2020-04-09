#!/usr/bin/env perl
#This script reads an aln in fasta format and modifies names so thathere are no duplicates
#if two sequences have the same name, the second one is renamed name_1 and so on
use strict;
use FileHandle;
use Env qw(HOST);
use Env qw(HOME);
use Env qw(USER);

my $format=file2format ($ARGV[0]);

if    ($format eq "clustalw"){clustalw2name_seq($ARGV[0]);}
elsif ($format eq "fasta")   {fasta2name_seq($ARGV[0]);}
elsif ($format eq "msf")   {msf2name_seq($ARGV[0]);}
elsif ($format eq "phylip")   {phylip2name_seq($ARGV[0]);}
elsif ($format eq "nameseq") {display_file ($ARGV[0]);}
 
exit (0);

sub file2format
  {
    my $f=shift;
    
    my $l=file2n_lines($f,2);
    
    if ( $l=~/^CLUSTAL/){return "clustalw";}
    elsif ($l=~/^SAGA/){return "clustalw";}
    elsif ($l=~/^>/){return "fasta";}
    elsif ($l=~/^PileUp/){return "msf";}
    elsif ($l=~/\s+\d+\s+\d+\s/){return "phylip";}
    elsif ($l=~/\#NAMESEQ_01/){return "nameseq";}
    else 
      {
	print STDERR "ERROR: $f FILE is NOT a supported format [ERROR]\n";
	exit (1);
      }
  }
sub display_file
    {
       my $file=shift;
       my $F= new FileHandle;
       open ($F, $file);
       while (<$F>){print "$_";}
       close ($F);
     }
sub phylip2name_seq
    {
      my $file=shift;
      my $F= new FileHandle;
      my ($seq, $name,$seq);
      my $query_start=-1;
      my $query_end=-1;
      my $in_aln=0;
      my %list;
      my ($first,$seq,$name, $cn, $nseq, $l,%len);
      
      open ($F, $file);
      <$F>;
      $l=$_;
      $l=~/\s*(\d+)\s*(\d+)/;
      $first=1;
      $cn=0;
      while (<$F>)
	{
	  my $l=$_;
	  if (!($l=~/\S/))
	    {
	      $cn=0;
	      $first=0;
	    }
	  elsif ($first==1)
	    {
	      $l=~/\s*(\S+)(.*)/;
	      my $name=$1;
	      my $seq=$2;
	      chomp ($seq);
	      $seq=~s/\s//g;
	      $list{$cn}{'name'}=$name;
	      $list{$cn}{'seq'}.=$seq;
	      $cn++;
	      $nseq++;
	    }
	  else
	    {
	      chomp ($l);
	      $l=~s/\s//g;
	      $list{$cn}{'seq'}.=$l;
	      $cn++;
	    }
	}
      close ($F);
      
      print "#NAMESEQ_01\n";
      print "# $nseq\n";
      for (my $a=0; $a<$nseq; $a++)
	{
	  my $nl=length ($list{$a}{'name'});
	  my $sl=length ($list{$a}{'seq'});
	  print ">$nl $sl $list{$a}{'name'} $list{$a}{'seq'}\n";
	}
    }
      
sub msf2name_seq
    {
      my $file=shift;
      my $F= new FileHandle;
      my ($seq, $name,$seq);
      my $query_start=-1;
      my $query_end=-1;
      my $in_aln=0;
      my %list;
      my ($seq,$name, $n, $nseq, $l,%len);
      
      open ($F, $file);
      while (<$F>)
	{
	  if ( /\/\//){$in_aln=1;}
	  elsif ( $in_aln && /(\S+)\s+(.*)/)
	    {
	      $name=$1;
	      $seq=$2;
	      $seq=~s/\s//g;
	      $seq=~s/\~/\-/g;
	      $seq=~s/\./\-/g;
	      if ( $list{$n}{'name'} && $list{$n}{'name'} ne $name)
		{
		  print "$list{$n}{'name'} Vs $name";
		  
		  exit (1);
		}
	      else
		{
		  $list{$n}{'name'}= $name;
		}
	      
	      $list{$n}{'seq'}=$list{$n}{'seq'}.$seq;
	      
	      $nseq=++$n;
	      
	    }
	  else
	    {$n=0;}
	}
      close ($F);
      print "#NAMESEQ_01\n";
      print "# $nseq\n";
      for (my $a=0; $a<$nseq; $a++)
	{
	  my $nl=length ($list{$a}{'name'});
	  my $sl=length ($list{$a}{'seq'});
	  print ">$nl $sl $list{$a}{'name'} $list{$a}{'seq'}\n";
	}
    }
    
sub fasta2name_seq
    {
      my $file=shift;
      my $F= new FileHandle;
      my ($seq, $name,$n,$l,%len);
      
      open ($F, $file);
      while (<$F>)
	{
	  if ( /^>(\S+)/){$n++;$seq="";$name=$1;}
	  else
	    {
	      $l=$_;
	      chomp ($l);
	      $seq.=$l;
	      $len{$name}=length($seq);
	    }
	}
      close ($F);
      print "#NAMESEQ_01\n";
      print "# $n";
      
      open ($F, $file);
      while (<$F>)
	{
	  my $l=$_;
	  $l=~s/\r[\n]*/\n/gm;
	  if ( ($l=~/^>(\S+)(.*)\n/))
	    {
	      my $name=$1;
	      my $comment=$2;
	      my $nl=length ($name);
	      my $sl=$len{$name};
	      if ($comment)
		{
		  $comment=~s/^\s+//g;
		  my $cl=length ($comment);
		  print "\n#$cl $comment\n";
		}
	      print "\n>$nl $sl $name ";
	    }
	  else
	    {
	      $l=$_;
	      chomp ($l);
	      print "$l";
	    }
	}
      print "\n";
      close ($F);
    }
sub clustalw2name_seq
  {
    my $fname=shift;
    my ($file1, $file);
    my (@blocks, @lines,@s, $n,$nseq, $c);
    my (@name, @seq);
    my $F= new FileHandle;
    my $stockholm;
    
    open ($F, $fname);
    
    while ( <$F>)
      {
	my $l=$_;
	$l=clean_cr($l);
	if ( $l=~/^# STOCKHOLM/){$stockholm=1;}
	elsif ( $stockholm && $l=~/^#/)
	  {
	    $l=~/^#(\S+)\s+(\S+)\s+(\S*)/g;
	    $l="_stockholmhasch_$1\_stockholmspace_$2 $3\n";
	  }
	$file.=$l;
      }
    close ($F);
        
    #Protect # and @
    $file=~s/\#/_hash_symbol_/g;
    $file=~s/\@/_arobase_symbol_/g;
    
    
    #Remove annotation
    $file=~s/\n[\.:*\s]+\n/\n\n/g;
    
    #Remove White spaces before the sequence name
    $file=~s/\n[ \t\r\f]+(\b)/\n\1/g;
    
    
    #Remove Internal Blanks
    $file=~s/(\n\S+)(\s+)(\S)/\1_blank_\3/g;
    
    $file=~s/[ ]//g;
    $file=~s/_blank_/ /g;
    
    
    #Identify Double Blank lines
    
    $file =~s/\n\s*\n/#/g;
    
    $file.="#";
    $file =~s/\n/@/g;
    
    
    
    
    #count nseq
    @blocks=split /\#/, $file;
    shift (@blocks);
    @s=split /\@/, $blocks[0];
    $nseq=$#s+1;
    
    #Merge all the sequences and split every Nseq
    
    
    $file=join '@', @blocks;
    @lines=split /\@/,$file;
    
    $c=0;
    
    foreach my $l (@lines)
      {
	my ($n, $s);
	
	if (!($l=~/\S/)){next;}
	elsif ($stockholm && ($l=~/^\/\// || $l=~/STOCKHOLM/)){next;}#get read of STOCHOLM Terminator
	
	$l=~/(\S+)\s+(\S*)/g;
	$n=$1; $s=$2;
	
	$seq[$c].=$s;
	$name[$c]=$n;
	$c++;
	
	if ( $c==$nseq){$c=0;}
	
      } 
    
    if ( $c!=0)
      {
	print STDERR "ERROR: $fname is NOT an MSA in Clustalw format: make sure there is no blank line within a block [ERROR]\n";
	exit (1);
      }
    print "#NAMESEQ_01\n";
    print "# $nseq\n";
    for (my $a=0; $a< $nseq; $a++)
      {
	$name[$a]=cleanstring ($name[$a]);
	$seq[$a]=cleanstring ($seq[$a]);
	my $ln=length ($name[$a]);
	my $ls=length ($seq[$a]);
	print ">$ln $ls $name[$a] $seq[$a]\n";
      }
  }
sub cleanstring
    {
      my $s=@_[0];
      $s=~s/_hash_symbol_/\#/g;
      $s=~s/_arobase_symbol_/\@/g;
      $s=~s/[ \t]//g;
      return $s;
    }

sub clean_cr
  {
    my $f=shift;
    $f=~s/\r\n/\n/g;
    $f=~s/\n\r/\n/g;
    $f=~s/\r\r/\n/g;
    $f=~s/\r/\n/g;
    return $f;
  }

sub file2n_lines
    {
      my $file=shift;
      my $nl=shift;
      my $ret;
      my $F=new FileHandle;
      my $n=0;
      open ($F, $file);

      while (<$F>)
	{
	  $ret.=$_;
	  $n++;
	  
	  if ($n>=$n){close ($F); return $ret;}
	}
      close ($F);
      return $ret;
    }
