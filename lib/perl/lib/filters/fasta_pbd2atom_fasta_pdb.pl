!#/bin/usr/env perl

$fasta=$ARGV[0];

open (F, $fasta);
$/='>';
while ( <F>)
  {
    print $_;
    die;
  }

