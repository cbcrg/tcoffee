#!/usr/bin/env perl

use Env qw(HOST);
use Env qw(HOME);
use Env qw(USER);
$random_tag=int (rand 10000)+1;
$unique_prefix="$$.$HOST.$random_tag";
$queue="distillery.and.mid";
$monitor=0;
$stderr_file="/dev/null";
$stdio_file="/dev/null";
$log_file="/dev/null";
$pause_time=0;
$max_sub_jobs=60;
$min_sub_jobs=30;
$output_all=0;
$var='\$';

foreach $value ( @ARGV)
    {
	if ($value ne $ARGV[$np]) 
	    {
	    ;
	    }
	elsif ($value eq "-max_sub_jobs")
	    {
	    $max_sub_jobs= $ARGV[++$np];
	    $np++;
    	    }	
	elsif ($value eq "-min_sub_jobs" )
	    {
	    $min_sub_jobs= $ARGV[++$np];
	    $np++;
    	    }
	elsif ($value eq "-para")
	    {
	    $para=1;
	    $monitor=1;
	    $np++;
    	    }
	elsif ($value eq "-monitor") 
	    {
	    $monitor=1;
	    $np++;
	    }
	elsif ($value eq "-no_monitor") 
	    {
	    $monitor=0;
	    $np++;
	    }
	elsif ($value eq "-queue")
	    {
	    $queue=$ARGV[++$np];
	    $np++;
	    }	
	elsif ($value eq "-stderr_file")
	    {
	    $stderr_file=$ARGV[++$np];
	    $np++;
	    }
	elsif ($value eq "-stdio_file")
	    {
	    $stdio_file=$ARGV[++$np];
	    $np++;
	    }
	elsif ($value eq "-output_all")
	    {
	    $output_all=1;
	    $np++;
	    }
	elsif ($value eq "-pause") 
	    {
	    $pause_time=$ARGV[++$np];
	    $np++;
	    }
	elsif ($value eq "-log")
	      {
	       $log=1;
	       
	       if ($ARGV[$np+1]=~/\-\S+/) 
	          {
		  $log_file="stderr";
	          }
	       else 
	          {
		  $log_file=$ARGV[++$np]; 
		  ++$np;
		 
	          }
	      }
	elsif ( $value eq "-com")
	    {
		
		if (!$ARGV[$np+1]=~/^\'/) { $com=$ARGV[++$np];}
		else {$com=$ARGV[++$np];}

	     $np++;
	    }
	elsif ( $value eq "-check")
	  {
	    
	    if (!$ARGV[$np+1]=~/^\'/) { $check=$ARGV[++$np];}
	    else {$check=$ARGV[++$np];}
	    $np++;
	  }
	elsif ($com eq "") 
	    {
	    $com_set=1;
	    $com=$ARGV[$np];
	    
	    $np++;
	    }
	elsif ($list eq "") 
	    {
	    $list_set=1;
	    $list=$ARGV[$np];
	    $np++;
	    }
	elsif ( $var_set eq "")
	    {
	    $var_set=1;
	    $var=$ARGV[$np];
	    $np++;
	    }
	}




if ( $com eq ""){print "You Need to Provide a Command [FATAL]\n";
	      die;
	     }



if ($list_set==0) 
    {
    $x= int (rand 100000)+1;
    $tmp_file_name="tmp_file_$x";
    open ( TMP, ">$tmp_file_name");
    while (<STDIN>)
      {
	print TMP $_;
      }
    close (TMP);
    open (F, $tmp_file_name);
    }
else 
    {
    open (F, $list);
    }

if ($para==0) 
    {

     @tc_list= <F>;
     close (F); 
     
     foreach $val(@tc_list) 
	    {
	      
	      
	      
	      $loc_com=$com;
	      if ($check){$loc_check=$check;}
	      
	      @i_val=($val=~/([^\s]+)/g);
	      
	      if ( $#i_val==0)
		{
		  if ($check){$loc_check=~s/$var/$i_val[0]/g;}
		  $loc_com=~s/$var/$i_val[0]/g;
		}
	      else
		{
		  for ($n=1; $n<=$#i_val+1;$n++ )
		    {
		      
		      $sub="$var$n";
		      
		      $loc_com=~s/$sub/$i_val[$n-1]/g;
		      if ($check){$loc_check=~s/$var/$i_val[0]/g;}
		    }
		}
	      if ( $check && -e $loc_check)
		{
		  print STDERR "skipping $loc_com...\n";
		  }
	      else
		{
		  system "$loc_com";
		}
	    }
    exit;
    }

elsif ($para==1) 
    {
    print STDERR "do parallel execution of: \"$com $list\"\n";
    
    if ($log==1) 
	{
	if ($log_file eq "stdout" || $log_file eq "stderr" ) 
		{
		$log_file="";
	        }

        else 
		{
		system "echo LOG FILE> $log_file";
		
	        }
	}
    else	
	{
	open ( OUT, ">/dev/null");
	}
	
    
    $id=0;
    $n_sub=0;
    while ($val=<F>) 
	    {	    	    
	    $job_log[$id]="$HOME/tmp/$unique_prefix.$id.log_file";
	    
	    $job=$unique_prefix."_$id";
	    open (JOB, ">$job");
	    
	    $loc_com=$com;
	    chop $val;

	    $loc_com=~s/\$/$val/g;
	 
	    print JOB "#!/bin/csh\n";
	    print JOB "#\$ -cwd\n";
	    print JOB "#\$ -N $unique_prefix\n";
	    if ($queue && !($queue eq " ")) {print JOB "#\$ -l $queue\n";}
	    print JOB "#\n";	    
            print JOB "$loc_com\n";
	    print JOB "echo FINISHED  >> $job_log[$id]\n";
	    print JOB "pwd\n";
	    
	    close (JOB);
	    if ( $output_all==1)
		{
		system "qsub $job >  $unique_prefix";		
	        }
	    else
		{system "qsub $job -e $stderr_file -o $stdio_file >$unique_prefix";	        
	        } 



	    print STDERR "$id: $output_all\n";
	    $n_sub++;
	    if ( $max_sub_jobs && $n_sub==$max_sub_jobs) 
		{
		$n_sub=monitor_process($min_sub_jobs,@job_log); 		 
		
	        }	
	   
            unlink $unique_prefix;
	    sleep $pause_time;
	    $id++;
	    }

    close (OUT);
    close (F);

    print STDERR "Your $id Jobs Have Been Submited (NAME=$unique_prefix)\n";
    monitor_process (0, @job_log);
    foreach $file(@job_log) {if (-e $file) {unlink($file);}}
    
    }

sub monitor_process ( @job_list)
    {
    my (@job_list)=@_;
    my $min_sub_jobs=shift (@job_list);
    my $n_sub_jobs;
    my $finished;
    my $n=0;

    $n_sub_jobs=-1;
    $finished=0;
    print STDERR "\nMonitor Batch: [$min_sub_jobs]";
       
    while (!$finished && (($n_sub_jobs>$min_sub_jobs)|| $n_sub_jobs==-1) ) 
	{
	$finished=1;
	$n_sub_jobs=0;
	$n=0;
	foreach $file (@job_list)
	        {
	
		if (-e $file){;}
		else 
		    {
		    $finished=0; $n_sub_jobs++;
	            }
	        }
	system "sleep 1";
        }
    
    return $n_sub_jobs;
    }
    
    
if ($tmp_file_name){unlink($tmp_file_name);}

