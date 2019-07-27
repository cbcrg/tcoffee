#!/usr/bin/env perl
use strict;

#This program will set the version tag to whatever is in $ARGV[0] or the file ./version_number-version
#if on MASTER system, it will set this value to the Github branch value
#it will return stdout

my $file;

my $major_versionF;
my $minor_versionF;
my $githubF;
my $version_numberF;

my $major_version;
my $minor_version;
my $github;
my $version_number;

my $beta
my $stable;
my $major;

my $path;

my $cl=join( " ", @ARGV);
if ( ($cl=~/-path=\s*(\S+)/))
  {
    $path=$1;
    $major_versionF=$path."/major_version.version";
    $minor_versionF=$path."/minor_version.version";
    $githubF=$path."/github.version";
    $version_number=$path."/version_number.version";
  }
if (!-e $major_versionF)return 1;
if (!-e $minor_versionF)return 1;

if (($cl=~/-beta/)){increase ($minor_versionF, 1);}
elsif ( ($cl=~/-stable/)){increase ($minor_versionF, 1);}
elsif ( ($cl=~/-major/)) {increase ($major_versionF, 1);}


if ($ENV{'TC_MASTER_NODE'}){system ("git branch -v  >$githubF");}



$minor_version=file2value($minor_versionF);
$major_version=file2value($major_versionF);
$github=github_file2value($githubF);

$version_number="$major_version\.$minor_version\.$github";
value2file($version_number, $version_numberF);
print "$version_number";



sub file2value
  {
    my ($file)=@_;
    open (F, "$file");
    while (<F>)
      {
	my $l=$_;
	chomp $l;
	close (F);
	return $l;
      }
  }
github_file2value
    {
      my ($file)=@_;
      my $value;
      
      open (F, "$file");
      while (<F>)
	{
	  my $l=$_;
	  chomp $l;
	  if ($l=~/\* master (\w+) .*/)
	    {
	      $value=$1;
	    }
	  
	  
	  close (F);
	return $l;
      }
  }
sub value2file
    {
      my ($value,$file)=@_;
      open (F, ">$file");
      print F "$vlaue\n";
      close F;
    }
    

sub increase
    {
      my ($file, $inc)=@_;
      my $value=$file2value($file);
      $value+=$inc;
      $value2file($value, $file);
      return $value;
    }








branchF;

my $version;
my $branch;



if( $ARGV[0]){$file=$ARGV[0];}
else{$file="version_number";}

$versionF=$file.".version";
$branchF=$file.".branch";


if ($ENV{'TC_MASTER_NODE'})
  {
    system ("git branch -v  > tmp.txt");
    open (F, "tmp.txt");
    while(<F>)
      {
	my $l=$_;
	if ($l=~/\* master (\w+) .*/)
	  {
	    $branch=$1;
	  }
      }
    close (F);
    open (F, ">$branchF");
    print F "$branch";
    close F;
    unlink ("tmp.txt");
  }


#get version
open (F,$versionF);
while (<F>)
          {
	  $version=$_;
	  chop ($version);
	}
close (F);



#get branch
open (F,$branchF);
while (<F>)
          {
	  $branch=$_;
	  chop ($branch);
	}
close (F);
print "$version.$branch";
