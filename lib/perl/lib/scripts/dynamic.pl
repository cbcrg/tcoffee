#!/usr/bin/env perl
use Env;
use strict;
use FileHandle;
use Cwd;
use File::Path;
use Sys::Hostname;

my %method;
my $method2use;
my $infile=$ARGV[0];
my $outfile=$ARGV[1];

my $n=file2nseq($infile);
my $tree=$ENV{"child_tree_4_TCOFFEE"};
my $method_file=$ENV{dynamic_config_4_TCOFFEE};
my $blast=$ENV{blast_server_4_TCOFFEE};
my $protein_db=$ENV{protein_db_4_TCOFFEE};
my $pdb_db=$ENV{pdb_db_4_TCOFFEE};
my $cache=$ENV{cache_4_TCOFFEE};
  
if (-e $method_file)
    {
      open (F, $method_file);
      while (<F>)
	{
	  my $f=$_;
	  $f=~/(\W)+ (\d)+/;
	  $method{$1}=$2;
	}
      close(F);
    }
else
  {
    $method{"psicoffee_msa"}=50;
    $method{"famsa_msa"}=1000000000;
  }

my $set=1;
foreach my $name (sort { $method{$a} <=> $method{$b} } keys %method) 
  {
    if ($n<=$method{$name})
      {
	$method2use=$name;
	last;
      }
  }
print "\n! METHOD: $method2use $blast $protein_db\n";
if ($tree)
  {
    ;
    #set the tree if needed
  }

if ($method2use eq "tcoffee_msa")
  {
    system ("t_coffee -in $infile -outfile $outfile -output fasta_aln -quiet");    
  }
elsif ($method2use eq "psicoffee_msa")
  {
    my $cacheFlag;
    if ($cache){$cacheFlag="-cache=$cache";}
      
    if ($blast eq "LOCAL")
      {
	system ("t_coffee -mode psicoffee -in $infile -outfile $outfile -output fasta_aln -blast_server=LOCAL -protein_db=$protein_db $cacheFlag -quiet");
      }
    else
      {
      system ("t_coffee -mode psicoffee -in $infile -outfile $outfile -output fasta_aln -quiet $cacheFlag");
    }
  }
elsif ($method2use eq "famsa_msa")
  {
    system ("famsa -t 1 $infile $outfile >/dev/null 2>/dev/null");
  }

else
  {
    system ("t_coffee -in $infile -method $method2use -outfile $outfile -output fasta_aln -quiet");
  }

sub file2nseq
  {
    my ($f)=@_;
    
    open (F, $f) || return 0;
    while (<F>)
      {
	if ($_=~/^\>/){$n++;}
      }
    close (F);
    return $n;
  }

