#!/usr/bin/env perl
use File::Copy;
use Env qw(HOST);
use Env qw(HOME);
use Env qw(USER);
#Specific Files and Default Parameters
$x_field=0;
$y_field=1;
$interval=0;
$file="stdin";

if ( !-e $ARGV[0]){open F, ">$ARGV[0]";}
else {open F, ">>$ARGV[0]";}

if ( $#ARGV==0)
  {
    while (<STDIN>)
      {
	$l=&clean_string ($_);
	print F $l;
      }
  }
elsif ( -e $ARGV[1])
  {
    open O, "$ARGV[1]";
    while (<O>)
      {
	$l=&clean_string ($_);
	print F $l;
      }
  close O;
  }
else
  {
    $l=&clean_string ($ARGV[1]);
    print F $l;
  }

close (F);

sub clean_string
  {
    my $l=@_[0];
    $l=~s/\\n/\n/g;
    $l=~s/\\t/\t/g;
    $l=~s/\#DOLLAR/\$/g;
    return $l;
  }

exit (0);
