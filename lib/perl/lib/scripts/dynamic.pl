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
my $thread=$ENV{"child_thread_4_TCOFFEE"};
my $dynamic=$ENV{dynamic_config_4_TCOFFEE};
my $blast=$ENV{blast_server_4_TCOFFEE};
my $protein_db=$ENV{protein_db_4_TCOFFEE};
my $pdb_db=$ENV{pdb_db_4_TCOFFEE};
my $cache=$ENV{cache_4_TCOFFEE};
my $tcarg=$ENV{tcarg_4_TCOFFEE};
my $template_file=$ENV{template_file_4_TCOFFEE};
my $clean;
my $treeFlag;
my $blastFlag;
my $infile;
my $outfile;
my $flush;
my $do_exit=0;
my ($h1, $h2);
my @tmpL;
my $tmpdir = File::Temp->newdir();
my $stderrF="$tmpdir/stderr";
$QUIET="2>$stderrF";
my $cdir=getcwd();
my $threadFlag4tc;
my $threadFlag4famsa;
my $threadFlag;



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
    elsif ($ARGV[$a] eq "-thread"){$thread=$ARGV[++$a]}
    elsif ($ARGV[$a] eq "-tcarg") {$tcarg=file2string($ARGV[++$a]);}
  }

$threadFlag=($thread)?"--thread $thread ":"--thread 1 ";
$threadFlag4tc=($thread)?"-thread $thread ":"-thread 1 ";
$threadFlag4famsa=($thread)?"-t $thread ":"-t 1 ";

if ($tree eq "list")
  {
    my $f="$tmpdir/f";
    open (F, ">$f");
    print F ">a\nxxx\n>b\nyyyyy\n";
    close (F);
    print STDOUT ("**** Supported Guide tree modes:\n");
    system ("t_coffee -other_pg seq_reformat -in $f -action +seq2dnd list ");
    $do_exit=1;
  }
if ($method2use eq "list")
  {
    my %ml;
    my $listfile="$tmpdir/list";
    
    $ml{tcoffee}=1;
    $ml{psicoffee}=1;
    $ml{accurate}=1;
    $ml{'3dcoffee'}=1;
    $ml{expresso}=1;
    $ml{clustalo}=1;
    $ml{mafft}=1;
    $ml{famsa}=1;
      print STDOUT ("**** Supported MSA mode:\n");
    system ("t_coffee 2>/dev/null | grep _msa > $listfile");
    open (F, $listfile);
    while (<F>)
      {
	my $l=$_;
	$l=~/(.*_msa)\s+(.*)/;
	my $m=$1;
	my $i="$2\n";
	if ($m=~/mafftsparsescore/)
	  {
	   printf STDOUT "%-20s DOES NOT Support [-tree] -- $i", $m;
	  }
	elsif ($m=~/mafft/)
	  {
	   printf STDOUT "%-20s DOES     Support [-tree] -- $i", $m;
	  }
	elsif (!$ml{$m})
	  {
	   printf STDOUT "%-20s DOES NOT Support [-tree] -- $i", $m;
	  }
	else 
	  {
	   printf STDOUT "%-20s DOES     Support [-tree] -- $i", $m;
	  }
      }
    $do_exit=1;
  }
if ($do_exit){exit ($EXIT_SUCCESS);}

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
	system ("t_coffee -other_pg seq_reformat -in $infile -action +seq2dnd $tree -output mafftdndmatrix> $treeF");
      }
    else
      {
	system ("t_coffee -other_pg seq_reformat -in $infile -action +seq2dnd $tree -output newick> $treeF");
      }
  }
my $psiCL=get_psicl();


if ($clean)
  {
    
    $infile=file2abs($infile);
    $outfile=file2abs($outfile);
    $treeF=file2abs($treeF);
    $protein_db=file2abs($protein_db);
    $pdb_db=file2abs($pdb_db);
    chdir ($tmpdir);
  }

if ($method2use eq "tcoffee_msa" || $method2use eq "tcoffee"|| $method2use eq "t_coffee" )
  {
    
    if ($treeF){$treeFlag="-usetree $treeF "}
    my_system ("t_coffee -seq $infile -outfile $outfile -output fasta_aln $treeFlag $threadFlag4tc $tcarg>/dev/null  $QUIET");    
  }

elsif ($method2use eq "mcoffee_msa" || $method2use eq "mcoffee")
  {
    
    if ($treeF){$treeFlag="-usetree $treeF "}
    my_system ("t_coffee  -mode fcoffee -seq $infile -outfile $outfile -output fasta_aln $treeFlag $threadFlag4tc $tcarg>/dev/null  $QUIET");    
  }
elsif ($method2use eq "fmcoffee_msa" || $method2use eq "fmcoffee")
  {
    
    if ($treeF){$treeFlag="-usetree $treeF "}
    my_system ("t_coffee  -mode fmcoffee -seq $infile -outfile $outfile -output fasta_aln $treeFlag $threadFlag4tc $tcarg>/dev/null  $QUIET");    
  }
elsif ($method2use eq "rcoffee_msa" || $method2use eq "rcoffee")
  {
    
    if ($treeF){$treeFlag="-usetree $treeF "}
    my_system ("t_coffee  -mode rcoffee -seq $infile -outfile $outfile -output fasta_aln $treeFlag $threadFlag4tc $tcarg>/dev/null  $QUIET");    
  }
elsif ($method2use eq "rmcoffee_msa" || $method2use eq "rmcoffee")
  {
    
    if ($treeF){$treeFlag="-usetree $treeF "}
    my_system ("t_coffee  -mode rmcoffee -seq $infile -outfile $outfile -output fasta_aln $treeFlag $threadFlag4tc $tcarg>/dev/null  $QUIET");    
  }
elsif ($method2use eq "rsapcoffee_msa" || $method2use eq "rsapcoffee")
  {
    
    if ($treeF){$treeFlag="-usetree $treeF "}
    my_system ("t_coffee  -mode rsapcoffee -seq $infile -outfile $outfile -output fasta_aln $treeFlag $threadFlag4tc $tcarg>/dev/null  $QUIET");    
  }
elsif ($method2use eq "psicoffee_msa" || $method2use eq "psicoffee")
  {
    my $cacheFlag;
    if ($cache){$cacheFlag="-cache=$cache";}
    if ($treeF){$treeFlag="-usetree $treeF "}
    
    if ($blast eq "LOCAL"){$blastFlag="-blast_server=LOCAL -protein_db=$protein_db";}
    print ("t_coffee -mode psicoffee -in $infile -outfile $outfile -output fasta_aln  $cacheFlag $blastFlag $threadFlag4tc $psiCL>/dev/null $QUIET");
    die;
    
    my_system ("t_coffee -mode psicoffee -in $infile -outfile $outfile -output fasta_aln  $cacheFlag $blastFlag $threadFlag4tc $psiCL>/dev/null $QUIET");
  }
elsif  ($method2use eq "accurate_msa" || $method2use eq "accurate")
  {
    my $cacheFlag;
    if ($cache){$cacheFlag="-cache=$cache";}
    if ($treeF){$treeFlag="-usetree $treeF "}
    if ($blast eq "LOCAL"){$blastFlag="-blast_server=LOCAL -protein_db=$protein_db -pdb_db=$pdb_db $psiCL";}
    my_system ("t_coffee -mode accurate -in $infile -outfile $outfile -output fasta_aln  $cacheFlag $blastFlag $threadFlag4tc>/dev/null  $QUIET");
  }
elsif  ($method2use eq "3dcoffee_msa"|| $method2use eq "3dcoffee")
  {
    my $cacheFlag;
    if ($cache){$cacheFlag="-cache=$cache";}
    if ($treeF){$treeFlag="-usetree $treeF "}
    if ($blast eq "LOCAL"){$blastFlag="-blast_server=LOCAL  -pdb_db=$pdb_db";}
    my_system ("t_coffee -method sap_pair TMalign_pair -template_file $template_file  -in $infile -outfile $outfile -output fasta_aln  $cacheFlag $blastFlag $threadFlag4tc>/dev/null $QUIET");
  }
elsif  ($method2use eq "3dmcoffee_msa"|| $method2use eq "3dmcoffee")
  {
    my $cacheFlag;
    if ($cache){$cacheFlag="-cache=$cache";}
    if ($treeF){$treeFlag="-usetree $treeF "}
    if ($blast eq "LOCAL"){$blastFlag="-blast_server=LOCAL  -pdb_db=$pdb_db";}
    my_system ("t_coffee -method mustang_pair sap_pair TMalign_pair -template_file $template_file  -in $infile -outfile $outfile -output fasta_aln  $cacheFlag $blastFlag $threadFlag4tc>/dev/null $QUIET");
  }

elsif  ($method2use eq "expresso_msa" || $method2use eq "expresso")
  {
    my $cacheFlag;
    if ($cache){$cacheFlag="-cache=$cache";}
    if ($treeF){$treeFlag="-usetree $treeF "}
    if ($blast eq "LOCAL"){$blastFlag="-blast_server=LOCAL  -pdb_db=$pdb_db";}
    my_system ("t_coffee -mode expresso -in $infile -outfile $outfile -output fasta_aln  $cacheFlag $blastFlag $threadFlag4tc>/dev/null $QUIET");
  }

elsif ($method2use eq "clustalo_msa" || $method2use eq "clustalo")
  {
    if ($treeF) {$treeFlag="--guidetree-in=$treeF "}
    
    my_system ("clustalo -i $infile $treeFlag -o $outfile  --force $threadFlag $QUIET");
    }
elsif ($method2use eq "mafft_msa" || $method2use eq "mafft")
  {
    if ($treeF){$treeFlag="--treein $treeF ";}
    my_system ("mafft --anysymbol $threadFlag $treeFlag $infile > $outfile $QUIET");
  }
elsif ($method2use eq "mafftginsi_msa" || $method2use eq "mafft-ginsi")
  {
    if ($treeF){$treeFlag="--treein $treeF ";}
    my_system ("mafft-ginsi $threadFlag --anysymbol $treeFlag $infile > $outfile $QUIET");
  }
elsif ($method2use eq "mafftfftns1_msa" || $method2use eq "mafft-fftns1")
  {
   
    if ($treeF){$treeFlag="--treein $treeF ";}
   
    my_system ("mafft-fftnsi $threadFlag --retree 1 --anysymbol $treeFlag $infile > $outfile $QUIET");
  }

elsif ($method2use =~/mafft/)
   {
     if ($treeF){$treeFlag="--treein $treeF ";}
     my_system ("$method2use $threadFlag --anysymbol $treeFlag $infile > $outfile $QUIET");
    }
elsif ($method2use eq "famsa_msa")
  {
    if ($treeF){$treeFlag="-gt import $treeF "}
    my_system ("famsa $treeFlag $threadFlag4famsa $infile $outfile >/dev/null $QUIET");
  }

else
  {
    if ($treeF)
      {
	printf (STDERR "WARNING: Method $method2use CANNOT use pre-sepecified guide tree [dynamic.pl]\n");
      }
    my_system ("t_coffee -in $infile -method $method2use -outfile $outfile -output fasta_aln $tcarg -quiet $QUIET");
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


if ($VERBOSE!=-1)
  {
    open (F, "$stderrF");
    while (<F>)
      {
	my $l=$_;
	if ( $VERBOSE || $l=~/WARNING/ || $l=~/ERROR/ || $l=~/INFORNATION/){print stderr "$l";}
      }
    close (F);
  }



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
sub file2string 
    {
      my ($f)=@_;
      my $s;

      open (F, $f) || return 0;
      while (<F>)
	{
	  $s.=$_;
	}
      close (F);
      chomp($s);
      return $s;
    }   

sub get_psicl
      {
	my ($psitrim, $psitrim_mode, $pisN);
	my $cl;
	
	if ($ENV{psitrim_tree_4_TCOFFEE}){$cl.=" -psitrim_tree=".$ENV{psitrim_tree_4_TCOFFEE}." ";}
	if ($ENV{psitrim_mode_4_TCOFFEE}){$cl.=" -psitrim_mode=".$ENV{psitrim_mode_4_TCOFFEE}." ";}
	if ($ENV{psitrim_4_TCOFFEE}){$cl.=" -psitrim=".$ENV{psitrim_4_TCOFFEE}." ";}
	if ($ENV{psiJ_4_TCOFFEE}){$cl.=" -psiJ=".$ENV{psiJ_4_TCOFFEE}." ";}
	

	return $cl;
      }
      
