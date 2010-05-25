#!/usr/bin/env perl

#This program reads a file (version_number.version by default)
#and outputs its content to STDOUT

if( $ARGV[0] && -e $ARGV[0]){$file=$ARGV[0];}
else{$file="version_number.version";}

open (INFILE,$file);

while (<INFILE>)
          {
	  $line=$_;
	  chop ($line);
	  print "$line";
	  }
close (INFILE);

