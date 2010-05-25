#!/usr/bin/env perl
#Version 1.01 (25/02/03)
use Cwd;
$base=cwd();

$CC="kcc";
$FC="h77";

$CC=get_C_compiler($CC);
$FC=get_F_compiler($FC);

open (F, "license.txt");
while (<F>)
  {
    print "$_";
  }
close (F);
print "Scroll up to read the entire license --- Type yes if you agree, or no if you do not\n";
$input=<STDIN>;
$input=lc ($input);
if ( !$input=~/yes/)
  {
    print "--------- You have NOT accepted the terms of the license -----\nInterruption of the installation\n";
    exit (EXIT_FAILURE)
  }
else
  {
    print "--------- You have accepted the terms of the license -----\nProceeding with the installation\n";
  }
  

$cdir=cwd();
chdir "t_coffee_source";
&flush_command ("make clean");
print "\n----------Compiling T-Coffee (Should take a few minutes)\n";
&flush_command ("make -i CC=$CC CFLAGS=\"-O2\" USER_BIN=../bin/ t_coffee");
  
print "\n----------Compiling TMalign (Should take a few seconds)\n";
&flush_command ("make -i USER_BIN=../bin/ FCC=$FC TMalign");
chdir $cdir;
$i=system (" perl -MCPAN -e 'install SOAP::Lite'");

if ( $i==1)
  {
    print "------- You Need To be Root to install SOAP:Lite with CPAN\n";
    print "------- log as root: with 'su root'\n";
    print "------- Type: 'perl -MCPAN -e install SOAP::Lite' (from any location)\n";
  }
else
  {
    print "SOAP is properly installed on your system\n"
  }

if ( !-e "./bin/t_coffee" && !-e "./bin/t_coffee.exe")
  {
    print "******* ERROR ***************************\n";
    print "******* T-Coffee installation failed ****\n";
    
    if ( $CC="")
      {
	print "******* You do not have an appropriate C Compiler ****\n";
      }
  }
else
  {
    print "******* T-Coffee was successfuly installed and is now in your bin directory\n";
  }
if ( !-e "./bin/TMalign" && !-e "./bin/TMalign.exe")
  {
    print "******* ERROR ***************************\n";
    print "******* TMalign installation failed ****\n";
  
    if ( $FC="")
      {
	print "******* You do not have an appropriate C Compiler ****\n";
      }
  }
else
  {
    print "******* TMalign was successfuly installed and is now in your bin directory\n";
  }


print "!!!!!!! Remember to add the bin directory to your \$PATH when the installation is finished!\n";

exit (EXIT_SUCCESS);  

sub get_C_compiler
  {
    my $c=@_[0];
    my (@clist)=("gcc", "cc", "ic64");
    
    return get_compil ($c, @clist);
 }

sub get_F_compiler
  {
    my ($c)=@_[0];
    my @clist=("f77", "g77", "gfortran", "ifort");
    return get_compil ($c, @clist);
  } 
       
sub get_compil
  {
    my ($c,@clist)=(@_);
    
    
    if (&check_pg_is_installed ($c)){return $c;}
    else
      {
	foreach $c (@clist)
	  {
	    if  (&check_pg_is_installed ($c)){return $c;}
	  }
      }
    return "";
  }

sub check_pg_is_installed
  {
    my @ml=@_;
    my $r, $p, $m;
    my $supported=0;
    
    my $p=shift (@ml);

    
    $r=`which $p 2>/dev/null`;
    if ($r eq "")
      {
	return 0;
      }
    else
      {
	return 1;
      }
  }
sub flush_command
  {
    my $command=@_[0];
    open (COMMAND, "$command|");
    while (<COMMAND>){print "--$_";}
    close (COMMAND);
  }
