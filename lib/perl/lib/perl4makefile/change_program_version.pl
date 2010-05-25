#!/usr/bin/env perl
#This program increments the current program version
#by 0.1. The current version must be by default in version_number.version
#and it must be of the format "Version_X.X\n"


if( $ARGV[0] && -e $ARGV[0]){$file=$ARGV[0];}
else{$file="version_number.version";}


if (!(-e $file))
	{
	$line=$ARGV[1];
	}
else
	{
	open (INFILE,$file);
	while (<INFILE>)
          	{
	  	$line=$_;
		}
	close (INFILE);
	}


$line=~/Version_(.*)\.([0-9]+)/;

$main=$1;
$secondary=$2;

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


open(INFILE,">$file");
print INFILE "Version_$main\.$secondary\n";
close (INFILE);

#$new_version=$2+1;

#print INFILE "Version_$1$new_version\n";
