#!/usr/bin/env perl

$x="x.html";
`wget http://igs-server.cnrs-mrs.fr/~cnotred/Projects_home_page/pavie.html -O$x`;
open (F, "$x");
while ( <F>)
  {

    if ( /EXAMPLE\:\<\/b\>([^<]+)/)
      {
	print "\n$1";
      }
  }
print "\n";
close (F);
