#!/usr/bin/env perl
use Env;
use FileHandle;
use Cwd;
use File::Path;
use Sys::Hostname;
use Cwd;
use File::Copy;
use DirHandle;
use strict;
use DateTime;
use File::Basename;
use File::Find;



my $target;
my $tag;
my $versionF;
my $version;
my $mode;

my $cl=join( " ", @ARGV);
if ( ($cl=~/-target=\s*(\S+)/)){$target=$1;}
if ( ($cl=~/-tag=\s*([^-]+[^-\s])/)){$tag=$1;}
if ( ($cl=~/-version=\s*(\S+)/)){$version=$1;}
if ( ($cl=~/-versionF=\s*(\S+)/)){$versionF=$1;}
if ( ($cl=~/-mode=\s*(\S+)/)){$mode=$1;}

if ($versionF){$version=file2value($versionF);}


if (!$mode)
  {
    file2edit_file ($target, $tag, $version);
  }
elsif ($mode eq "conf.py")
  {
   
    my $in=new FileHandle();
    my $out=new FileHandle();
    my $new_file;

    open ($in, $target);
    while (<$in>)
     {
       my $l=$_;
       if ( $l=~/^version = u\'/)
	 {
	   $l="version = u\'$version\'\n";
	 }
       elsif ($l=~/^release = u\'/)
	  {
	   $l="release = u\'$version\'\n";
	 }
      $new_file.=$l;
     }
    close ($out);
    open ($in, ">$target");
    print $in "$new_file";
    close ($out);
  }

  


sub file2edit_file
      {
       my ($file, $tag, $version)=@_;
       my $new_file;
       open (F, "$file") or die "Could not open file [$file][FATAL:edit_version.pl\n";
       while (<F>)
	 {
	   my $line=$_;
	  
	   if ( $line=~/$tag/)
	     {
	       
	       $line=~s/$tag/our \$VERSION=\"$version\";#POPULATED by edit_version.pl from $tag/;
	     }
	   $new_file.=$line;
	 }
       close (F);
       open (F, ">$file") or die "Could not open file [$file][FATAL:read_program_version.pl\n";
       print F "$new_file";
       close (F);
       }

