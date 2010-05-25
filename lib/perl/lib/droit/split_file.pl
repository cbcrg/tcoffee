#! /usr/bin/env perl


open ( F, $ARGV[0]);
open ( G, "0.txt");
while (<F>)
  {
    $line=$_;
    if ( $line =~/#Code/)
      {
	print $_;
	$line=~/(\S+)#$/;
	$name=$1;

	open ( G, ">$name.txt");
      }
    print G "$line";
  }
close (G);


