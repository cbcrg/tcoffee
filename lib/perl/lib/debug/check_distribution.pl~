#!/usr/bin/env perl

$HOME=$ENV{'HOME'};
$cache="$HOME/.t_coffee/cache4debug/";
$ENV{'CACHE_4_TCOFFEE'}=$cache;

#defaults

$data="$HOME/Dropbox/projects/c.svn/trunk/t_coffee/doc_test/data/";
$archive="$HOME/untared_distributions/";
$package="T-COFFEE";
$source="t_coffee_source";
$action="diff";

foreach $a (@ARGV){$cl.=" $a ";}
if ( $cl =~/\-h/)
  {
    print "DESCRIPTION:\n";
    print "  Checks a list of command lines on a given distribution\n";
    print "USAGE:\n";
    print "   check_distribution -file <file> ¦ -doc <doc> -tag <tag> -para <extra CL parameters>\n";
    print "PARAMETERS:\n";
    print "  -doc..................name of the msword doc file\n" ;
    print "  -tag..................name of the tag used to extract the CL from the msword file\n";
    print "  -file.................name of the file containing the command lines\n" ;
    print "  -para.................extra parameters\n";
    print "  -data.................dir containing the data to test";
    die;
  }

$doc="$HOME/projects/c/t_coffee/doc/t_coffee_doc.doc";

$tag="PROMPT:";



if ($cl=~/\-file([^-]*)/)
  {
    $file=$1;$file=~s/\s//g;
  }
if ($cl=~/\-data([^-]*)/)
  {
    $data=$1;$data=~s/\s//g;
  }

if ( $cl=~/\-doc([^-]*)/)
  {
    $doc=$1;$doc=~s/\s//g;
  }
if ($cl=~/\-tag([^-]*)/)
  {
    $tag=$1;$tag=~s/\s//g;
    if ($tag eq "none"){$tag="";}
  }
if ($cl=~/\-fo([^-]*)/)
  {
    $file_only=1;
  }

if ($cl=~/\-extra_param(.*)/)
  {$extra_param=$1;}


#1 prepare the data
if ( $data ne "current")
  {
    `mkdir check_distribution_dir`;
    `rm check_distribution_dir/*`;
    chdir "check_distribution_dir";
    
    `cp $data/* .`;
  }

#2: create the command line file
if (!$file)
  {
    $file="tmp_check_distribution.cl_file";
    print "antiword $doc | grep '$tag' | sed s/$tag//g > $file \n";
    `antiword -w0 $doc | grep '$tag' | sed s/$tag//g > $file`;
  }

if ( $file_only)
  {
    print STDERR "List of Commands in File: $file\n";
    exit (EXIT_SUCCESS);
  }



#run the check
$failed=$passed=0;
$version=`t_coffee -version`;

if ( -e $file)
  {
    open (F, $file);
    open (FJ, ">failed_list.txt");
    while (<F>)
      {
	$action= $_;
	chomp ($action);
	$n++;

	if ( $action =~/tcp/){$action=~s/tcp/t_coffee -other_pg seq_reformat/;}

	$action=~s/\[.*\]//;
	if ($action =~/\S/){$action.=" $extra_param ";}
	&act ($action);
      }
    close (F);
    close (FJ);
    $success_rate=(($passed+$failed)==0)?0:($passed*100)/($passed+$failed);
    print STDERR "\nAnalysis Finished: VERSION: $version FAILED=$failed PASSED=$passed SUCCESS=$success_rate %%\n";
  }
#clean the cache dir:
`rm -r $cache/*`;
sub act 
  {
    my $command=@_[0];
    
    if (!$command =~/\S/ || $command=~/mocca/){return;}
    $inc=$command;
    if ($command =~/>/)
      {$command .=" 2>/dev/null";}
    else
      {$command .=" >$n.log 2>&1";}

    $r=system ($command);

    if ( $r!=0)
      {
	$failed++;
	print STDERR "\n*** [FAILED] $command";
	print FJ "$inc\n";
      }
    else
      {
	$passed++;
	
	print STDERR "\n[$n:OK]: $command";
      }
    return;
  }
