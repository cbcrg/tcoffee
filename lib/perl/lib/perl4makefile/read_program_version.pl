#!/usr/bin/env perl
use strict;

#This program will set versio_number.version
#if on MASTER system, the github branch numb er will be updated
#it will return stdout

my $file;
my $build_versionF;
my $minor_versionF;
my $major_versionF;
my $githubF;
my $version_numberF;


my $build_version;
my $minor_version;
my $major_version;
my $github;
my $version_number;

my $path;
my $circleciF;
my $pr;

my $cl=join( " ", @ARGV);

if ( ($cl=~/-path=\s*(\S+)/)){$path=$1;}
else {$path="./";}

if ( ($cl=~/-circleci=\s*(\S+)/)){$circleciF=$1;}
else {$circleciF="./.circleci/config.yml";}

if ( ($cl=~/-print/)){$pr=1;}
else {$pr=0;}

##


$build_versionF=$path."/build_version.version";
$minor_versionF=$path."/minor_version.version";
$major_versionF=$path."/major_version.version";
$githubF=$path."/github.version";
$version_numberF=$path."/version_number.version";


if (!-e $build_versionF){print STDERR "build_version File [$build_versionF] could not be opened [FATAL:read_program_version.pl]\n";die;}
if (!-e $minor_versionF){print STDERR "minor_version File [$minor_versionF] could not be opened [FATAL:read_program_version.pl]\n";die;}
if (!-e $major_versionF){print STDERR "major_version File [$major_versionF] could not be opened [FATAL:read_program_version.pl]\n";die;}



if (($cl=~/-beta/))
  {
    increase ($build_versionF, 1);
    if ($ENV{'TC_MASTER_NODE'}){system ("git branch -v  >$githubF");}
    circleci2releaseV($circleciF, "0");
  }
elsif ( ($cl=~/-stable/))
  {
    value2file(0, $build_versionF);
    increase   ($minor_versionF, 1);
    if ($ENV{'TC_MASTER_NODE'}){system ("git branch -v  >$githubF");}
    circleci2releaseV($circleciF, "1");
  }
elsif ( ($cl=~/-major/))
  {
    value2file(0, $build_versionF);
    value2file(0, $minor_versionF);
    increase   ($major_versionF, 1);
    if ($ENV{'TC_MASTER_NODE'}){system ("git branch -v  >$githubF");}
    circleci2releaseV($circleciF, "1");
  }

$build_version=file2value($build_versionF);
$minor_version=file2value($minor_versionF);
$major_version=file2value($major_versionF);
$github=github_file2value($githubF);

$version_number="Version_$major_version\.$minor_version\.$build_version\.$github\n";
value2file($version_number, $version_numberF);

if ($pr){print "$version_number"; }



sub file2value
  {
    my ($file)=@_;
    open (F, "$file") or die "Could not open file [$file][FATAL:read_program_version.pl\n";
    while (<F>)
      {
	my $l=$_;
	chomp $l;
	close (F);
	return $l;
      }
  }
sub github_file2value
    {
      my ($file)=@_;
      my $value;
      
      open (F, "$file") or die "Could not open file [$file][FATAL:read_program_version.pl\n";
      while (<F>)
	{
	  my $l=$_;
	  chomp $l;
	  if ($l=~/\* master (\w+) .*/)
	    {
	      $value=$1;
	    }
	  close (F);
	  return $value;
	}
    }

sub value2file
    {
      my ($value,$file)=@_;
      open (F, ">$file");
      print F "$value";
      close F;
    }
    

sub increase
    {
      my ($file, $inc)=@_;
      my $value=file2value($file);
      $value+=$inc;
      value2file($value, $file);
      return $value;
    }

sub circleci2releaseV
      {
       my ($file, $value)=@_;
       my $new_file;
       
       open (F, "$file") or die "Could not open file [$file][FATAL:read_program_version.pl\n";
       while (<F>)
	 {
	   my $line=$_;
	   print $line;
	   if ( $line=~/\#UPDATE_RELEASE_STATUS/)
	     {
	       $line=~s/RELEASE=\".\"/RELEASE=$value/;
	     }
	   $new_file.=$line;
	 }
       close (F);
       open (F, ">$file") or die "Could not open file [$file][FATAL:read_program_version.pl\n";
       print F "$new_file";
       close (F);
       }
      






