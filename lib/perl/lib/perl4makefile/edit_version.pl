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
use File::Basename;
use File::Find;

my ($inF,$tag,$string,$outF);

for (my $i=0; $i<=$#ARGV; $i++)
  {
    my $cl=$ARGV[$i];
    if   ($cl eq "-in") {$inF=$ARGV[++$i];}
    elsif($cl eq "-tag"){$tag =$ARGV[++$i];}
    elsif($cl eq "-string") {$string   =clean($ARGV[++$i]);}
    elsif($cl eq "-out"   ) {$outF  =$ARGV[++$i];}
  }

my $in=new FileHandle();
my $out=new FileHandle();

if (!$tag)   { die "No -tag provided [FATAL:read_program_version.pl\n";}
if (!$string){ die "No -string Provided [FATAL:read_program_version.pl\n";}



open ($in, $inF) or die "Could not open file [$in][FATAL:read_program_version.pl\n";

my $new_file;
while (<$in>)
  {
    my $l=$_;
    if ( $l=~/$tag/)
      {
	if ($tag=~/#_#/){$l="$string $tag -- Populated by edit_version.pl\n";}
	else {$l="$string\n";}
      }
    $new_file.=$l;
  }
close ($in);



if (!$outF){$outF=$inF;}
open ($out, ">$outF");
print $out "$new_file";
close ($out);



sub clean
  {
    my $s= shift @_;

    if ($s=~/#_#DQ/){$s=~s/#_#DQ/\"/g;}
    if ($s=~/#_#BQ/){$s=~s/#_#BQ/\`/g;}
    if ($s=~/#_#SQ/){$s=~s/#_#SQ/\'/g;}
    if ($s=~/#_#DS/)
      {
	
	$s=~s/#_#DS/\$/g;
	
	
      }#Dollar_SIGN
    if ($s=~/#_#AS/){$s=~s/#_#AS/\@/g;}
    return $s;
  }
	

