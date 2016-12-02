#!/usr/bin/env perl

my $address="http://www.tcoffee.org/Data/Datasets/NatureProtocolsDataset.tar.gz";
my $out="NatureProtocolsDataset.tar.gz";
&url2file ($address,$out);

if ( -e $out)
  {
    
    system ("gunzip NatureProtocolsDataset.tar.gz");
    system ("tar -xvf NatureProtocolsDataset.tar");
  	system ("rm -rf NatureProtocolsDataset.tar");  
    print "Your Data Set is in the Folder 'NatureProtocolsDataset'\n";
  }
else 
  {
    print "Could not Download Dataset --- Web site may be down -- Try again later\n";
  }




sub url2file
{
    my ($address, $out, $wget_arg, $curl_arg)=(@_);
    my ($pg, $flag, $r, $arg, $count);
    
    if (!$CONFIGURATION){&check_configuration ("wget", "INTERNET", "gzip");$CONFIGURATION=1;}
    
    if (&pg_is_installed ("wget"))   {$pg="wget"; $flag="-O";$arg=$wget_arg;}
    elsif (&pg_is_installed ("curl")){$pg="curl"; $flag="-o";$arg=$curl_arg;}
    return system ("$pg $address $flag $out>/dev/null 2>/dev/null");

}

sub pg_is_installed
  {
    my @ml=@_;
    my $r, $p, $m;
    my $supported=0;
    
    my $p=shift (@ml);
    if ($p=~/::/)
      {
	if (system ("perl -M$p -e 1")==$EXIT_SUCCESS){return 1;}
	else {return 0;}
      }
    else
      {
	$r=`which $p 2>/dev/null`;
	if ($r eq ""){return 0;}
	else {return 1;}
      }
  }
sub check_configuration 
    {
      my @l=@_;
      my $v;
      foreach my $p (@l)
	{
	  
	  if   ( $p eq "EMAIL")
	    { 
	      if ( !($EMAIL=~/@/))
		{
		  exit (EXIT_FAILURE);
		}
	    }
	  elsif( $p eq "INTERNET")
	    {
	      if ( !&check_internet_connection())
		{
		  exit (EXIT_FAILURE);
		}
	    }
	  elsif( $p eq "wget")
	    {
	      if (!&pg_is_installed ("wget") && !&pg_is_installed ("curl"))
		{
		  exit (EXIT_FAILURE);
		}
	    }
	  elsif( !(&pg_is_installed ($p)))
	    {
	      exit (EXIT_FAILURE);
	    }
	}
      return 1;
    }
sub check_internet_connection
  {
    my $internet;
    my $tmp;
    &check_configuration ( "wget"); 
    
    $tmp=&vtmpnam ();
    
    if     (&pg_is_installed    ("wget")){`wget www.google.com -O$tmp >/dev/null 2>/dev/null`;}
    elsif  (&pg_is_installed    ("curl")){`curl www.google.com -o$tmp >/dev/null 2>/dev/null`;}
    
    if ( !-e $tmp || -s $tmp < 10){$internet=0;}
    else {$internet=1;}
    if (-e $tmp){unlink $tmp;}

    return $internet;
  }

sub vtmpnam
      {
	my $r=rand(100000);
	my $f="file.$r.$$";
	while (-e $f)
	  {
	    $f=vtmpnam();
	  }
	push (@TMPFILE_LIST, $f);
	return $f;
      }


