#!/usr/bin/env perl
$HOME=$ENV{'HOME'};
$archive="$HOME/untared_distributions/";
if ( $#ARGV==-1)
  {
    print "\nDESCRIPTION:\n";
    print "  This program checks bug reports from T-Coffee\n";
    print "USAGE:\n";
    print "  test_error_report <file> -version <vn>\n";
    print "PARAMETERS:\n";
    print "  <file>..........T-Coffee Error Report\n";
    print "  -version........Version Number or Higher\n";
    print "  -testlog........Version Number or Higher\n";
    print "\n";
    die;
  }

$error_file=$ARGV[0];
$dir="$error_file\_dir";
`mkdir $dir`;
`cp $error_file $dir`;
chdir $dir;

if ( $#ARGV==-1)
  {
    print "This program checks bug reports from T-Coffee\n";
    print "-version version_number | current........Set version to use\n";
  }

foreach $arg (@ARGV)
  {
    $command_line.=" $arg";
  }


open (E, $error_file);
while (<E>)
  {   
   
    if (/######### FILE: (\S*)/)
	{
	 open ( F, ">$1");
	 while (<E>)
		{
		  if ( /(.*)END    \#\#\#\#\#\#/){last;}
		  print F  $_;
		}
	 close (F);
	} 
    elsif  (/######### OUTFILE: (\S*)/)
	    {
	   
	     open (F, ">ref_$1");
	     push @outfile, $1;
             while (<E>) 
		{ 
		   if ( /(.*)END    \#\#\#\#\#\#/){last;}
		   print F  $_;
		}  
	     close (F);
	    }
    elsif (/######### PROGRAM_VERSION/)
      {
       
       while (<E>)
	  {
	    if ( /\S/ && /(.*)END    \#\#\#\#\#\#/){last;}
	    /(\S+), Version_(\S+)/;
	    $program=$1;
	    $version=$2;
	  }
      }
    elsif ( /######### COMMAND_LINE/)
      {
	while (<E>)
	  {
	    if (/\S/ && /(.*)END    \#\#\#\#\#\#/){last;}
	    $command=$_;
	  }	
      }

  }
	   
@pg=($command=~/(\S+)/g);
if (($command_line=~/\-purify/))
    {
      $command=~s/$pg[0]/$pg[0]\.purify/;
    }

if (($command_line=~/\-version/))
    {
      $command_line=~/\-version (\S+)/;
      $version=$1;      
    }
	   
if ( $version eq "current")
	   {
	     $command=$command;
	   }
else
	 {
	   $c="$archive$program\_distribution_Version_$version/bin/$command";
	   if ( -e $c){$command=$c;}
	 }

if (($command_line=~/\-check_log/))
	   {
	     $command=~s/\-full_log (\S+)//;
	     chop $command;
	     $command .=" -quiet";

	     `$command`;
	     print "\n#### COMMAND: $command";
	     foreach $f (@outfile)
	       {
		 $value=`diff $f ref_$f`;
		 print "\n#### Differences on file $f:\n$value";
	       }
	   }
else
	   {
	    `$command`;
	    print "\n$command";
	  }
	    
