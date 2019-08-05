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
    elsif($cl eq "-string") {$string   =$ARGV[++$i];}
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
	$l="$string $tag -- Populated by edit_version.pl\n";
      }
    $new_file.=$l;
  }
close ($in);



if (!$outF){$outF=$inF;}
open ($out, ">$outF");
print $out "$new_file";
close ($out);





