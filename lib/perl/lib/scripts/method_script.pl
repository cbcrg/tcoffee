#!/usr/bin/env perl
use Env;
#This is a sample script_file, it reads.
#It expets: argv[0]: -in
#           argv[1]: file_name
#           argv[2]: -out
#           argv[3]: file_name
#           argv[4]: program_name
#           argv[5]: extra parameters


$infile=$ARGV[1];
$outfile=$ARGV[3];
$method=$ARGV[5];
for ( $a=6; $a<=$#ARGV; $a++){$param.=" $ARGV[$a] ";}
$param.=" >/dev/null 2>&1 ";

if (!$infile)
  {
    print STDERR "generic wrapper for alignment methods\n";
    print STDERR "-in <infile> -out <outfile> -meth method <parmeter list>\n";
    print STDERR "Supports: clustalw\n";
  }
else if ($method eq "clustalw")
  {
    $command="clustalw -infile=$infile -outfile=$outfile $param";
    sytem ( $command);
  }
#YOUR METHOD HERE
else
  {
    print STDERR "$method Is Not a Supported Method";
  }


if ( -e $outfile){exit (EXIT_SUCCESS);}
else
  {
    print STDERR, "\nCommand $command\nDid Not Produce File $ARGV[3]\n[FATAL]\n";
    exit (EXIT_FAILURE);
  }
