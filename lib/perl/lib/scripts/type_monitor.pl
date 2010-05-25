#! /usr/bin/env perl
system "stty cbreak </dev/tty >/dev/tty 2>&1";

if ($#ARGV==-1)
  {
    print "type_monitor -txt=<file containing text to type> -char=<max number char> -time=<max time> -log=<keystrokes> -summary=<summary> -user=<username>\n";
  }

$cl=join (" ", @ARGV);

($TXT)=(&my_get_opt ( $cl, "-txt=",1,1));
($LOG)=(&my_get_opt ( $cl, "-log=",0,0));
($SUMMARY)=(&my_get_opt ( $cl, "-summary=",0,0));
($TIME)=(&my_get_opt ( $cl, "-time=",0,0));
($CHAR)=(&my_get_opt ( $cl, "-char=",0,0));
($USER)=(&my_get_opt ( $cl, "-user=",0,0));

if (!$USER)
  {
    if (-e "run"){open (F, "run"); while (<F>){$user=$_; chomp ($user);}close (F);}
    $user++;
    open (F, ">run");
    print F "$user\n";
    close (F)
  }



if (!$TIME){$TIME=-1;}
if (!$CHAR){$CHAR=-1;}
if (!$SUMMARY){$SUMMARY="summary\_$user.txt";}
if (!$LOG){$LOG="log\_$user.txt";}



open (F, "$TXT");
open (LOG, ">$LOG");
open (SUMMARY, ">$SUMMARY");


$char=$char_error=$char_no_error=$char_correct=0;
while (<F>) 
  {
    $line=$_;
    $current="";
    &check_continue();
    
    print "$line";
    @list=split(//,$line);
    foreach $c(@list)
     {
       &check_continue();
       ($c, $value1)=&clean_char($c);
       
       $d=getc;$char++;
       $no_error=1;
       
       ($d, $value2)=&clean_char ($d);
       
       while($d ne $c)
	 {
	   print LOG "FAIL -- $value1 $value2\n";
	   print "\a";
	   print "\r$current";
	   $fail{"$value1 $value2"}++;
	   $fail{"$value1"}++;
	   
	   &check_continue();
	   $char_error++;
	   $d=getc;$char++;
	   ($d, $value2)=&clean_char ($d);
	   $no_error=0;
	  
	   
	   
	 }
       if ($no_error){$pass{$value1}++;$char_no_error++;}
       $char_correct++;
       print LOG  "PASS -- $value1\n";
       $current.="$c";
       $bg{$value1}++;
     }
  }
my_exit ();

sub check_continue
  {
    if (!$ref_time){$ref_time=time;}
    $time=time-$ref_time;
    if ($TIME != -1 && $time >$TIME) {&my_exit();}
    if ($CHAR !=-1 && $char>$CHAR){&my_exit();}
    
  }

sub my_exit
  {
    
    print SUMMARY "#type monitor\n";
    print SUMMARY "#1 - Time: $time sec.\n";
    print SUMMARY "#2 - Number of typed character (correct or incorrect) : $char\n";
    print SUMMARY "#3 - Number of mistyped target character              : $char_error\n";
    print SUMMARY "#4 - Number of target character typed without error   : $char_no_error\n";
    print SUMMARY "#5 - Number of target character typed w or wo error   : $char_correct\n";
    
    foreach $e (keys (%fail))
      {
	print SUMMARY "ERR -- $e $fail{$e}\n";
      }
    foreach $e (keys (%pass))
      {
	print SUMMARY "PAS-- $e $pass{$e}\n";
      }
    close (SUMMARY);
    close (LOG);
    print STDERR "\nSUMMARY in: $SUMMARY\n";
    print STDERR "\nLOG     in: $LOG\n";
    exit (EXIT_SUCCESS);
  }


sub clean_char
  {
    my $c=@_[0];

    if ($c eq "\n")
      {return ($c, "newline");}
    elsif ($c=~/\s/)
      {return (" ","space");}
    else 
      {return ($c, $c);}
  }


sub my_get_opt
  {
    my @list=@_;
    my $cl, $a, $argv, @argl;
    
    @argl=();
    $cl=shift @list;
    
    for ( $a=0; $a<=$#list; $a+=3)
      {
	$option=$list[$a];
	$optional=$list[$a+1];
	$status=$list[$a+2];
	$argv="";
	if ($cl=~/$option(\S+)/){$argv=$1;}
	@argl=(@argl,$argv);
	
	
	#$optional:0=>optional
	#$optional:1=>must be set
	#$status: 0=>no requirement
	#$status: 1=>must be an existing file
	#$status: 2=>must be an installed package
	

	if ($optional==0){;}
	elsif ( $optional==1 && $argv eq "")
	  {
	    print STDERR "ERROR: Option $option must be set [FATAL:$program/$mode/$method]\n";
	    exit (EXIT_FAILURE);
	  }
	if ($status==0){;}
	elsif ($status ==1 && $argv ne "" && !-e $argv)
	  {
	    print STDERR "ERROR: File $argv must exist [FATAL:$program/$mode/$method]\n";
	    exit (EXIT_FAILURE);
	  }
	elsif ( $status==2 && $argv ne "" && &check_pg_is_installed ($argv)==0)
	  {
	    print STDERR "ERROR: $argv is not installed [FATAL:$program/$mode/$method]\n";
	    exit (EXIT_FAILURE);
	  }
      }

    return @argl;
    }
