#!/usr/bin/env perl
use Env;

if ($argv[0] eq "-support"){$support=1;}

foreach $a (@ARGV)
  {
    $in.="$a ";
  }
print ( STDERR "IN: $in");
if ($in=~/-infile=([^-]*)/){$infile=$1;$infile=~s/\s//g;}
if ($in=~/-outfile=([^-]*)/){$outfile=$1;$outfile=~s/\s//g;}
if ($in=~/-method=([^-]*)/){$method=$1;$method=~s/\s//g;}
if ($in=~/-support/){$support=1;}
if ($in=~/-param=([^-]*)/){$param=$1;}

$param.=" >/dev/null 2>&1 ";
@method_list=("clustalw", "t_coffee","muscle");

if (!$infile)
  {
    print STDERR "generic wrapper for alignment methods\n";
    print STDERR "\t-infile=<infile> -outfile=<outfile> -method=<method> -parameter=<parmeter list>\n";
    print STDERR "Supported Methods: @method_list\n";
    exit (EXIT_SUCCESS);
  }
else
  {
    check_pg_is_installed ($method, @method_list);
  }
if ($method eq "clustalw")
  {
    $command="clustalw -infile=$infile -outfile=$outfile $param";
    `$command`;
  }
elsif ($method eq "t_coffee")
  {
    $command="t_coffee2 -infile=$infile -outfile=$outfile $param";
    `$command`;
  }
elsif ($method eq "muscle")
  {
    `muscle $infile > $outfile`;
  }
else
  {
  `$method -infile=$infile -outfile=$outfile $param`;  
  }

if ( -e $outfile){exit (EXIT_SUCCESS);}
else
  {
    print STDERR "\nCommand $command Did Not Produce File $outfile [FATAL]\n";
    exit (EXIT_FAILURE);
  }

sub check_pg_is_installed
  {
    my @ml=@_;
    my $r, $p, $m;
    my $supported=0;
    
    my $p=shift (@ml);
    $r=`which $p`;
    if ($r eq "")
      {
	print STDERR "\nProgram $p Supported but Not Installed on your system [FATAL:tc_generic_method]\n";
	exit (EXIT_FAILURE);
      }
    else
      {
	return 1;
      }
  }
$program="T-COFFEE (Version_2.36)";\n

