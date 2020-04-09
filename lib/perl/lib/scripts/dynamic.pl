#!/usr/bin/env perl
use Env;
use strict;
use FileHandle;
use Cwd;
use File::Path;
use Sys::Hostname;
use File::Temp qw/ tempfile tempdir /;
my $QUIET="2>/dev/null";
my $VERBOSE=0;
our $EXIT_FAILURE=1;
our $EXIT_SUCCESS=0;

my %method;
my $method2use;
my $treeF;
my $tree=$ENV{"child_tree_4_TCOFFEE"};
my $dynamic=$ENV{dynamic_config_4_TCOFFEE};
my $blast=$ENV{blast_server_4_TCOFFEE};
my $protein_db=$ENV{protein_db_4_TCOFFEE};
my $pdb_db=$ENV{pdb_db_4_TCOFFEE};
my $cache=$ENV{cache_4_TCOFFEE};
my $template_file=$ENV{template_file_4_TCOFFEE};
my $clean;
my $treeFlag;
my $blastFlag;
my $infile;
my $outfile;
my $flush;
my ($h1, $h2);
my @tmpL;
my $tmpdir = File::Temp->newdir();
my $cdir=getcwd();




for ($a=0; $a<=$#ARGV; $a++)
  {
    if    ($ARGV[$a] eq "-seq"){$infile=$ARGV[++$a];}
    elsif ($ARGV[$a] eq "-outfile"){$outfile=$ARGV[++$a];}
    elsif ($ARGV[$a] eq "-dynamic_config"){$dynamic=$ARGV[++$a];}
    elsif ($ARGV[$a] eq "-blast_server"){$blast=$ARGV[++$a];}
    elsif ($ARGV[$a] eq "-tree") {$tree=$ARGV[++$a];}
    elsif ($ARGV[$a] eq "-protein_db") {$protein_db=$ARGV[++$a];}
    elsif ($ARGV[$a] eq "-pdb_db") {$pdb_db=$ARGV[++$a];}
    elsif ($ARGV[$a] eq "-cache") {$cache=$ARGV[++$a];}
    elsif ($ARGV[$a] eq "-method") {$method2use=$ARGV[++$a];}
    elsif ($ARGV[$a] eq "-verbose"){$VERBOSE=1; $QUIET="";} 
    elsif ($ARGV[$a] eq "-clean"){$clean=1;}
    elsif ($ARGV[$a] eq "-template_file"){$template_file=$ARGV[++$a]}
    
  }
if ($ENV{dynamic_4_TCOFFEE})
  {
    print "ERROR - $method2use resulst in infinite recursion [FATAL:dynamic.pl]\n";
    exit ($EXIT_FAILURE);
  }
else
  {
    $ENV{dynamic_4_TCOFFEE}=1;
  }

my $n=file2nseq($infile);
if ($n==0)
  {
    print "ERROR - No sequences provided [FATAL:dynamic.pl]\n";
    exit ($EXIT_FAILURE);
  }
if (!$outfile)
  {
    ($h1,$outfile)=tempfile();
    push (@tmpL,$outfile);
    $flush=1;
  }

if ($method2use){;}
else 
  {
    if (-e $dynamic)
      {
	open (F, $dynamic);
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
    
    foreach my $name (sort { $method{$a} <=> $method{$b} } keys %method) 
      {
	if ($n<=$method{$name})
	  {
	    $method2use=$name;
	    last;
	  }
      }
  }
if ($VERBOSE){print "\n! METHOD: NSEQ: $n $method2use $blast $protein_db\n";}
if ($tree)
  {
    ($h2,$treeF)=tempfile();
    push (@tmpL,$treeF);
    if ( -e $tree)
      {
	system ("cp $tree $treeF");
      }
    elsif ($method2use=~/mafft/)
      {
	system ("t_coffee -other_pg seq_reformat -in $infile -action +seq2dnd $tree -output mafttdndmatrix> $treeF");
      }
    else
      {
	system ("t_coffee -other_pg seq_reformat -in $infile -action +seq2dnd $tree -output newick> $treeF");
      }
  }

if ($clean)
  {
    
    $infile=file2abs($infile);
    $outfile=file2abs($outfile);
    $treeF=file2abs($treeF);
    $protein_db=file2abs($protein_db);
    $pdb_db=file2abs($pdb_db);
    chdir ($tmpdir);
  }
	

if ($method2use eq "tcoffee_msa")
  {
    if ($treeF){$treeFlag="-usetree $treeF "}
    my_system ("t_coffee -seq $infile -outfile $outfile -output fasta_aln $treeFlag >/dev/null  $QUIET");    
  }
elsif ($method2use eq "psicoffee_msa")
  {
    my $cacheFlag;
    if ($cache){$cacheFlag="-cache=$cache";}
    if ($treeF){$treeFlag="-usetree $treeF "}
    if ($blast eq "LOCAL"){$blastFlag="-blast_server=LOCAL -protein_db=$protein_db";}
    my_system ("t_coffee -mode psicoffee -in $infile -outfile $outfile -output fasta_aln  $cacheFlag $blastFlag>/dev/null $QUIET");
  }
elsif  ($method2use eq "accurate_msa")
  {
    my $cacheFlag;
    if ($cache){$cacheFlag="-cache=$cache";}
    if ($treeF){$treeFlag="-usetree $treeF "}
    if ($blast eq "LOCAL"){$blastFlag="-blast_server=LOCAL -protein_db=$protein_db -pdb_db=$pdb_db";}
    my_system ("t_coffee -mode accurate -in $infile -outfile $outfile -output fasta_aln  $cacheFlag $blastFlag>/dev/null  $QUIET");
  }
elsif  ($method2use eq "3dcoffee_msa")
  {
    my $cacheFlag;
    if ($cache){$cacheFlag="-cache=$cache";}
    if ($treeF){$treeFlag="-usetree $treeF "}
    if ($blast eq "LOCAL"){$blastFlag="-blast_server=LOCAL  -pdb_db=$pdb_db";}
    my_system ("t_coffee -method sap_pair TMalign_pair -template_file $template_file  -in $infile -outfile $outfile -output fasta_aln  $cacheFlag $blastFlag>/dev/null $QUIET");
  }

elsif  ($method2use eq "expresso_msa")
  {
    my $cacheFlag;
    if ($cache){$cacheFlag="-cache=$cache";}
    if ($treeF){$treeFlag="-usetree $treeF "}
    if ($blast eq "LOCAL"){$blastFlag="-blast_server=LOCAL  -pdb_db=$pdb_db";}
    my_system ("t_coffee -mode expresso -in $infile -outfile $outfile -output fasta_aln  $cacheFlag $blastFlag>/dev/null $QUIET");
  }

elsif ($method2use eq "clustalo_msa")
  {
    if ($treeF){$treeFlag="--guidetree-in=$treeF "}
    my_system ("clustalo -i $infile $treeFlag -o $outfile  --force --threads=1 $QUIET");
    }
elsif ($method2use eq "mafft_msa")
  {
    if ($treeF){$treeFlag="--treein $treeF ";}
    my_system ("mafft --anysymbol $treeFlag $infile > $outfile $QUIET");
  }
elsif ($method2use eq "mafftginsi_msa")
  {
    if ($treeF){$treeFlag="--treein $treeF ";}
    my_system ("mafft-ginsi --anysymbol $treeFlag $infile > $outfile $QUIET");
  }
elsif ($method2use eq "mafftfftns1_msa")
  {
    if ($treeF){$treeFlag="--treein $treeF ";}
    my_system ("mafft-fftns1 --anysymbol --retree1 $treeFlag $infile > $outfile $QUIET");
  }
elsif ($method2use eq "famsa_msa")
  {
    if ($treeF){$treeFlag="-gt import $treeF "}
    my_system ("famsa -t 1 $infile $outfile >/dev/null $QUIET");
  }

else
  {
    if ($treeF)
      {
	printf (STDERR "WARNING: Method $method2use CANNOT use pre-sepecified guide tree [dynamic.pl]\n");
      }
    my_system ("t_coffee -in $infile -method $method2use -outfile $outfile -output fasta_aln -quiet $QUIET");
  }
if ($clean)
  {
    chdir ($cdir);
  }

#Flush output if none provided
if ( ! -e $outfile)
  {
    print "ERROR - No MSA computed [FATAL:dynamic.pl]\n";
    exit ($EXIT_FAILURE);
  }
elsif ( $flush)
 {
   open (F, $outfile);
   while (<F>){print $_;}
   close (F);
 }
#Clean empty files
foreach my $f (@tmpL){unlink($f);}

exit ($EXIT_SUCCESS);


sub my_system 
  {
    my ($com)=@_;
    if ($VERBOSE)
      {
	print "$com\n";
      }
    system ($com);
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
sub file2abs
    {
      my $f=shift @_;
      
      
      if (!$f || $f=~/^\//){return $f;}
      return "$cdir/$f";
    }
