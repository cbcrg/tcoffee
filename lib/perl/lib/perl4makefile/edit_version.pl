#!/usr/bin/env perl
use Cwd;
use File::Path;

my $target;
my $tag;
my $versionF;
my $version;


if ( ($cl=~/-target=\s*(\S+)/)){$target=$1;}
if ( ($cl=~/-tag=\s*(\S+)/)){$tag=$1;}
if ( ($cl=~/-version=\s*(\S+)/)){$version=$1;}
if ( ($cl=~/-versionF=\s*(\S+)/)){$versionF=$1;}

if ($versionF){$version=file2value($versionF);}
file2edit_file ($targe, $tag, $version);

sub file2edit_file
      {
       my ($file, $tag, $version)=@_;
       my $new_file;
       open (F, "$file") or die "Could not open file [$file][FATAL:read_program_version.pl\n";
       while (<F>)
	 {
	   my $line=$_;
	   print $line;
	   if ( $line=~/$tag/)
	     {
	       $line=~s/$tag/VERSION=$version/;
	     }
	   $new_file.=$line;
	 }
       close (F);
       open (F, ">$file") or die "Could not open file [$file][FATAL:read_program_version.pl\n";
       print F "$new_file";
       close (F);
       }

