#!/usr/bin/env perl
#This program increments a tag VersionTag by 0.1
#It respects any symbols AROUND VersionTag = xx.xx



$file=$ARGV[0];
$tmp="$$.tmpfile";

if ( !-e $file){die;}

open (F, $file);
open (O, ">$tmp");

while (<F>)
  {
    if ( /(.*)VersionTag\s*(.*)\s*([0-9]+)\.([0-9]+)(.*)/)
	 {
	   $main=$3;
	   $secondary=$4; 
	   if ( $secondary>=99)
	     {
	       $main++;
	       $secondary=1;
	     }
	   else
	     {
	       $secondary++;
	     }
	   
	   if ($secondary<10){$secondary="0$secondary";}
	   print O "$1\VersionTag $2 $main\.$secondary$5\n";
	   $success=1;
	 }
    else
      {
	print O "$_";
      }
  }
close (F);
close (OUT);
`mv $tmp $file`;
`chmod u+x $file`;
if ( $success){exit EXIT_SUCCESS;}
else
  {
    print STDERR "Could not change program version [ERROR:change_perl_program_version]\n";
    exit EXIT_FAILURE;
}
