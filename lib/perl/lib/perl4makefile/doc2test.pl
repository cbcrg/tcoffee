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
my $RETURN=" ####RETURN#### ";
my %cl;
my %dir;
#$GIT=0: no git interaction
#$GIT=1: Commit all new files, uncommit all deleted files
#$GIT=2: blank run


my $GIT=1;

my $max=0;
my $pattern='$$';
my $log ;
my $doclog;
my $regexp;
my $reset;
my $stop_on_failed;
my $purge;
my $clean;

my $redo_failed;
my $failed;
my $clf;
my $rep;
my $cw=cwd();
my $mode="new";#will only run the new ones

my $pg="t_coffee";

$dir{examples}="$cw/examples/";
$dir{docs}    ="$cw/docs/";
$dir{tmp}     ="$cw/testsuite/validation/doc/tmp/";
$dir{ref}     ="$cw/testsuite/validation/doc/ref/";
$dir{latest}  ="$cw/testsuite/validation/doc/latest/";
$dir{log}     ="$cw/testsuite/validation/doc/log/";
$dir{failed}  ="$cw/testsuite/validation/doc/failed/";

if ($ARGV[0] eq "-help")
  {
    print "doc2check\n";
    print "Automaticly checks t_coffee command lines\n";
    print "The github dir structure is expected by default\n";
    print "tcoffee/\n";
    print "       /docs     -> contains rst doc\n";
    print "       /examples -> contains the reference files\n";
    print "       /testsuite/validation/doc/\n"; 
    print "       /testsuite/validation/doc/tmp    -> computation\n";
    print "       /testsuite/validation/doc/ref    -> succesful dumps\n";
    print "       /testsuite/validation/doc/failed -> unsuccesful dumps\n";
    print "\n";
    print "Commands are extracted from the .rst files contained in <-doc>\n";
    print "Commands are recognised as any line starting with <-pattern>\n";
    print "Duplicated commands are checked only once\n";
    print "By default the program only checks the new commands (-mode=new) \n";
    print "To check All the commands against the references, use -mode validate\n";
    print "flags:\n";
    
    print "     -command     <file containing a list of commands, for test purpose>\n";
    print "     -pattern       pattern used to recognize the command lines [def=**]\n";
    print "                    pattern will be treated as a regexp if -regexp is set\n";
    print "     -regexp        flag that causes pattern to be treated as a perl regexp\n";
    print "     -pg            specify the path of the version of T-Coffee (optional)\n";
    print "     -log           default: validation.log\n";
    print "     -doclog        default: doc.log\n";
    print "     -mode=<action> new|update|failed\n";
    print "                    new    : check ONLY CL w/o ref/dump and create ref/dump\n";
    print "                    update : check ALL  CL or create new ref/dum\n";
    print "                    failed : run   ONLY FAILURE as found /failed\n";
    print "     -reset         delete all /ref and /failed dumps before running mode\n";
    print "     -compile       stop at every FAILURE\n";
    print "     -clean         removes unused /ref/dumps, /failed/dumps and /examples files\n";
    
    
    print "     -rep           specifies a root repository\n";
    print "     -example_dir   directory containing the sample files\n"; 
    print "     -docs_dir      directory containing the .rst files\n"; 
    print "     -ref_dir       directory containing the reference dump\n"; 
    print "     -run_dir       tmp directory\n";
    
    print "     -max           max number of CL to check [DEBUG]\n";
    print "     -help          display this help message\n";
    
    
    
    print "\n";
    print "\n";
    die;
    }

@ARGV=clean_cl(@ARGV);

for (my $a=0; $a<=$#ARGV; $a++)
  {
    
    
    if ($ARGV[$a]=~/-command/)
      {
	$clf=$ARGV[++$a];
      }
    elsif ($ARGV[$a]=~/-pattern/)
      {
	$pattern=$ARGV[++$a];
      }
    elsif ($ARGV[$a]=~/-regexp/)
      {
	$regexp=1;
      }
    elsif ($ARGV[$a]=~/-reset/)
      {
	$reset=1;
      }
    elsif ($ARGV[$a]=~/-stop_on_failed/)
      {
	$stop_on_failed=1;
      }
    elsif ($ARGV[$a]=~/-clean/)
      {
	$clean=1;
      }
    elsif ($ARGV[$a]=~/-max/)
      {
	$max=$ARGV[++$a];
      }
    elsif ($ARGV[$a]=~/-pg/)
      {
	$pg=$ARGV[++$a];
      }
    elsif ($ARGV[$a]=~/-log/)
      {
	$log=$ARGV[++$a];
      }
    elsif ($ARGV[$a]=~/-doclog/)
      {
	$doclog=$ARGV[++$a];
      }
    elsif ($ARGV[$a]=~/-mode/)
      {
	#default: new      => will only run the new CL
	#         validate => will check everything
	$mode=$ARGV[++$a];
      }
    
    elsif ($ARGV[$a]=~/-example_dir/)
      {
	$dir{examples}=$ARGV[++$a];
      }
    elsif ($ARGV[$a]=~/-failed_dir/)
      {
	$dir{failed}=$ARGV[++$a];
      }
    elsif ($ARGV[$a]=~/-docs_dir/)
      {
	$dir{docs}=$ARGV[++$a];
      }
    elsif ($ARGV[$a]=~/-tmp_dir/)
      {
	$dir{tmp}=$ARGV[++$a];
      }
    elsif ($ARGV[$a]=~/-ref_dir/)
      {
	$dir{ref}=$ARGV[++$a];
      }
    elsif ($ARGV[$a]=~/-rep/)
      {
	$rep=$ARGV[++$a];
      }
    elsif ($ARGV[$a]=~/-mode/)
      {
	$mode=$ARGV[++$a];
      }
  }

#protect special char in pattern if non regexp mode is used
if (!$regexp)
  {my @char=split (//, $pattern);
   my $newP;
   foreach my $c (@char)
     {
       if ($c eq "*" | $c eq "+" || $c eq '$' || $c eq "." || $c eq "?" || $c eq "|" )
	 {
	   $newP.="\\$c";
	 }
     }
   $pattern=$newP;
 }

if ($rep)
  {
    foreach my $d (keys (%dir))
      {
	my $s=$dir{$d};
	$s=~s/$cw/$rep/;
	$dir{$d}=$s;
      }
  }
if (!$log){$log="$dir{log}/validation.log";}
open (LOG, ">$log");close (LOG);
if (!$doclog){$doclog="$dir{log}/doc.log";}


#1-minimum checks

if (!-d $dir{docs})
  {
    print STDERR "ERROR: The $dir{docs} file must be set to the dir containing the doc in .rst format [FATAL]\n";
    die;
  }
if (!-d $dir{examples})
  {
    print STDERR "ERROR: The $dir{examples} file must must be set to the dir containing the examples [FATAL]\n";
    die;
  }

# 2 - create the directory structure if needed
foreach my $d (keys (%dir))
  {
    if (!-d $dir{$d})
      {
	print "create $dir{$d}\n"; 
	system ("mkdir -p $dir{$d}");
      }
  }

# reset
if ($reset==1)
  {
    
    my @dir_list=("ref", "failed", "tmp");
    foreach my $d (@dir_list)
      {
	print STDERR "Reset $dir{$d}\n";
	my @list=dir2file_list($dir{$d});
	
	foreach my $f (@list)
	  {
	    if ($f=~/\.dump/)
	      {
		myunlink ("$dir{$d}/$f");
		print STDERR "\tRemove $dir{$d}/$f\n";
	      }
	  }
      }
    exit;
  }


# 2 - get the command line file

my @rstlist;
if (!$clf)
  {
   opendir (DIR, $dir{docs});
   my @list=readdir (DIR);
   closedir (DIR);
   foreach my $rst (@list)
     {
       my $f="$dir{docs}/$rst";
       if ($rst=~/.*\.rst$/)
	 {
	   rst2cl($f, \%cl);
	   push @rstlist, $f;
	 }
     }
 }
else
  {
    open (F, "$clf");
    while(<F>)
      {
	my $l=$_;
	$l=~/\s*(\S.*\S)\s*/;
	my $line=$1;
	
	if (($line=~/\S/))
	  {
	    $cl{$line}{0}{"doc"}=$line;
	  }
      }
    close (F);
  }


#dump the documentation

docdump (\%cl, $doclog);


# 3 - get the examples in the tmp dir
if ($clean)
  {
    clean_examples ($dir{examples}, $doclog);
    mygit_add ($dir{examples});
    exit;
  }

system ("cp $dir{examples}/* $dir{tmp}");
clean_examples ($dir{tmp}, $doclog);

# 4 - get the command lines in the reference and in the failed
my $num;
$num=dir2dump (\%cl,\%dir,"ref", $num); 
$num=dir2dump (\%cl,\%dir,"failed", $num); 

#display_command (\%cl);
#die;

# 5 - initialize the report
my $report;
my $datestring = localtime();
$report="##REPORT -- TIME -- $datestring\n";
my $version=`$pg -version`;chomp $version;
$report.="##REPORT -- SOFTWARE -- $pg\n";
$report.="##REPORT -- VERSION  -- $version\n";
print "$report";
open (LOG, ">>$log");print (LOG "$report");close (LOG);
$report="";


# 5 - run the analysis
chdir ($dir{tmp});
my (%R, %D);
my  @com_list=keys (%cl);
my $n_com=$#com_list+1;
my $count=1;
my $shell;
my ($tot_failure, $tot_success, $tot_warning, $tot_tot);
$tot_failure=$tot_warning=$tot_success=$tot_tot=0;
foreach my $com (sort (keys (%cl)))
  {
    my $dump;
    my $dump_num;
    my $ref=$cl{$com}{0}{ref};
    my $doc=$cl{$com}{0}{doc};
    my $failed=$cl{$com}{0}{failed};
    my $run=0;

    if (!($com=~/\S/)){next;}


    if (!$ref)
      {
	$num++;
	$dump="t_coffee.$num.dump";
	$dump_num=$num;
	
      }
    elsif ($doc && $ref)
      {
	$dump="$ref";
	$dump=~/t_coffee\.(\d+)\.dump/;
	$dump_num=$1;
      }


    #if ($cl{$com}{0}{doc}}{$cl{$com}{0}{doc}=$dump;}


    my $shell_com=$com;

    if ($pg ne "t_coffee")
      {
	$shell_com=~s/t_coffee/$pg/g;
      }
    
    #decide which jobs will run
    $run=0;
    if    ($mode eq "new")   {if (!$ref)   {$run=1;}}
    elsif ($mode eq "failed"){if ( $failed){$run=1;}}
    else  {$run=1;}
    
    $count++;
    if ($run)
      {
	my $success;
	if ($max && $tot_tot>=$max){last;}
	
	$tot_tot++;
	$cl{$com}{run}=1;
	%R=dump2report($com,"ref", \%cl);
	
	
	$shell=system ("export DUMP_4_TCOFFEE=$dump;$shell_com >/dev/null 2>/dev/null");
	
	%D=dump2report($com,"doc", \%cl);
	compare_reports (\%R, \%D);
	if ($shell || $D{error} || $D{MissingOutput})
	  {
	    #refresh the example files
	    system ("rm $dir{tmp}/*");
	    system ("cp $dir{examples}/* $dir{tmp}");
	    clean_examples ($dir{tmp}, $doclog);
	    $shell=system ("export DUMP_4_TCOFFEE=$dump;$shell_com >/dev/null 2>/dev/null");
	    %D=dump2report($com,"doc", \%cl);
	    compare_reports (\%R, \%D);
	  }

	
	if (! -e $dump)
	  {
	    open (F, ">$dump");
	    print F "<dumpIO><cl>$com</cl><stack>#FAILURE -- SHELL ERROR -- Check the Command Line Syntax</stack><file><stream>output</stream><name>stderr</name><content>#FAILURE -- SHELL ERROR -- Check the Command Line Syntax</content></file></dumpIO>\n";
	    close(F);
	  }
	else
	  {
	    substitute_cl4dump($dump,$com);
	  }

#Prepare the report	
	$report="";	
	if ($D{error} || $D{MissingOutput})
	  {
	    $tot_failure++;
	    $cl{$com}{status}="failure";
	    for (my $a=0; $a< $cl{$com}{n}; $a++)
	      {
		$report.="##FAILURE [$dump_num] -- COM -- $com -- Manual: $cl{$com}{$a}{file} -- Line: $cl{$com}{$a}{line} -- $dump\n";
	      }
	    if ($D{error})
	      {
		my $stack=$D{stack};
		$stack=~s/$RETURN/\n/g;
		
		if ($stack=~/\S/)
		  {
		    $report.=wraptext ($D{stack}, "##FAILURE [$dump_num] -- MSG -- ");
		  }
		else 
		  {
		    my $nstderr;
		    my $stderr=$D{stderr};
		    
		    $stderr=~s/$RETURN/\n/g;
		    my @ll=split (/\n/, $stderr);
		    
		    foreach my $l (@ll){if ($l=~/ERROR/){$nstderr.="$l\n";}}
		    $stderr=$nstderr;
		    
		    $report.=wraptext ($stderr, "##FAILURE [$dump_num] -- MSG -- ");
		  }
	      }
	  }
	elsif ($shell)
	  {
	    $tot_failure++;
	    $cl{$com}{status}="failure";
	    for (my $a=0; $a< $cl{$com}{n}; $a++)
	      {
		$report.="##FAILURE [$dump_num] -- COM -- $com -- Manual: $cl{$com}{$a}{file} -- Line: $cl{$com}{$a}{line} -- dump: $dump\n";
	      }
	    $report .="##FAILURE [$dump_num] -- SHELL ERROR -- Check the Command Line Syntax\n";
	    
	  }
	elsif ($D{warning})
	  {
	    $tot_warning++;
	    $cl{$com}{status}="warning";
	    $success=1;
	    for (my $a=0; $a< $cl{$com}{n}; $a++)
	      {
		$report.="##WARNING [$dump_num] -- COM -- $com -- Manual: $cl{$com}{$a}{file} -- Line: $cl{$com}{$a}{line} -- dump: $dump\n";
	      }
	  }
	else
	  {
	    $tot_success++;
	    $cl{$com}{status}="success";
	    $success=1;
	    for (my $a=0; $a< $cl{$com}{n}; $a++)
	      {
		$report.="##SUCCESS [$dump_num] -- COM -- $com -- Manual: $cl{$com}{$a}{file} -- Line: $cl{$com}{$a}{line} -- dump: $dump\n";
	      }
	  }

	print "##$count Out of $n_com\n$report";
	open (LOG, ">>$log");print (LOG "$report");close (LOG);


	

	if ($success)
	  {
	    if (!$ref)   {mymove ($dump,$dir{ref});}
	    else {mymove ($dump,$dir{latest});}
	    if ($failed) {purgedir4dump($dir{"failed"}, $com);}
	  }
	elsif (!$success)
	  {
	    mymove ($dump,$dir{failed});
	  }
	
	if (!$success && $stop_on_failed){exit;}
      }
  }


#Finalize the report

$report="## SUMMARY -- COMPLETE -- TOT: $tot_tot FAILURE: $tot_failure WARNING: $tot_warning SUCCESS: $tot_success\n";
foreach my $rst (@rstlist)
  {
    my ($failure, $success, $warning, $tot);
    $tot=$failure=$warning=$success=0;
    foreach my $com (keys (%cl))
      {
	if ($cl{$com}{run})
	  {
	    for (my $a=0; $a<$cl{$com}{n}; $a++)
	      {
		if ($cl{$com}{$a}{file} eq $rst)
		  {
		    if    ($cl{$com}{status} eq "failure"){$failure++;}
		    elsif ($cl{$com}{status} eq "warning"){$warning++;}
		    elsif ($cl{$com}{status} eq "success"){$success++;}
		    else {print "---Unknown status!!!$cl{$com}{status}---\n";}
		    $tot++;
		  }
	      }
	  }
      }
    $report.="## SUMMARY -- TOT: $tot FAILURE: $failure WARNING: $warning SUCCESS: $success -- $rst\n";
  }
print "$report";
open (LOG, ">>$log");print (LOG "$report");close (LOG);
purge(\%cl, \%dir);

die;

#end main block
sub docdump
  {
    my ($cl, $docdump)=@_;
    my %new;
    
    
    foreach my $c (keys (%{$cl}))
      {
	for (my $a=0; $a<$cl->{$c}{n}; $a++)
	  {
	    my $rst=$cl->{$c}{$a}{file};
	    my $line=$cl->{$c}{$a}{line};
	    $new{$rst}{$line}=$c;
	  }
      }
    open (F, ">$docdump");
    foreach my $file (keys (%new))
      {
	print F "## File: $file\n";
	foreach my $l (sort {$a <=> $b} (keys(%{$new{$file}})))
	  {
	    printf F "Line: %-4d --- $new{$file}{$l}\n", $l;
	  }
      }
    close (F);
    return;
  }
sub file2word_hash
    {
      my ($file,$h)=@_;
      open (F, $file);
      while (<F>)
	{
	  my $line=$_;
	  $line =~s/[=,;:\+\-]/ /g;
	  my @wl=split (/\s+/, $line);
	  foreach my $w (@wl)
		{
		  $w=~/.(.*)/;
		  my $sw=$1;
		  $h->{$sw}{n}++;
		  $h->{$w} {n}++;
		}
	}
      close (F);
    }
sub check_file_list
    {
      my ($dir,$hash)=@_;
      my $check=0;
      foreach my $f (keys (%{$hash}))
	{
	  if (!$hash->{$f}{checked})
	    {
	      if (!-e "$dir/$f")
		{
		  $hash->{$f}{checked}=1;
		  $hash->{$f}{exists}=-1;
		  
		}
	      else
		{
		  $check++;
		  $hash->{$f}{checked}=1;
		  $hash->{$f}{exists}=-1;
		  file2word_hash ("$dir/$f", $hash);
		}
	    }
	}
      return $check;
    }
	
sub clean_examples
    {
      my ($dir,$doc)=@_;
      my %names;
      
      if (!-d $dir){return;}
      file2word_hash ($doc, \%names);
      while (check_file_list ($dir, \%names)){;}

      opendir (DIR, $dir);
      my @fl=readdir (DIR);
      closedir (DIR);
      foreach my $f (@fl)
	{
	  if (!$names{$f})
	    {
	      #print " ------ rm $dir/$f\n";
	      myunlink("$dir/$f");
	    }
	}
    }

sub purge
  {
    my ($cl, $dir)=@_;
    my $purge;
    
    print "Purge deprecated Command Lines ...\n";
    foreach my $d ( keys (%{$dir}))
      {
	my $cd= $dir->{$d};
	print "Purge $cd\n";
	my @l=dir2file_list($cd);

	foreach my $dump (@l)
	  {
	    $dump="$cd/$dump";
	    if ($dump=~/\.dump/)
	      {
		my $com=dump2cl($dump);
		if (!$cl->{$com}{0}{doc})
		  {
		    print "-----  Purge $dump\n";
		    myunlink ("$dump");
		  }
	      }
	  }
      }
    return $purge;
  }


sub purgedir4dump
  {
    my ($d, $com)=@_;
    my (@list, $purged);
    my @list=dir2file_list ($d);
    
    foreach my $dump (@list)
      {
	
	if ($dump=~/t_coffee\.(\d+).dump/ && !($dump=~/.*dump.new.*/))
	  {
	    $dump="$d/$dump";
	    my $cd=dump2cl($dump);
	    if ($cd eq $com)
	      {
		myunlink ($dump);
		$purged++;
	      }
	  }
      }
    return $purged;
  }
sub dir2dump
  {
    my ($cl, $dir, $mode, $num)=@_;
    my ($d,@list);
    
    $d=$dir->{$mode};
    my @list=dir2file_list ($d);
        
    foreach my $dump (@list)
      {
	$dump="$d/$dump";
	if ($dump=~/t_coffee\.(\d+).dump/ && !($dump=~/.*dump.new.*/))
	  {	
	    my $cn=$1;
	    $num=($cn>$num)?$cn:$num;
	    my $com=dump2cl($dump);
	    $cl->{$com}{0}{$mode}=$dump;
	  }
      }
    return $num;
  }
sub compare_reports
  {
    my ($ref, $doc)=@_;


    if (!$ref || !$doc){return;}
    
    foreach my $f (keys (%{$ref->{file}}))
      {
	if (!$doc->{file}{$f})
	  {
	    $doc->{MissingOutputF}{$f}=1;$doc->{MissingOutput}++;
	  }
	elsif ($doc->{file}{$f} ne $ref->{file}{$f})
	  {
	    $doc->{DifferentOutput}{$f}=1;$doc->{N_DifferentOutput}++;
	  }
      }
    return;
  }
sub dump2report
  {
      my ($com, $type, $cl)=@_;
      my %R;  
      my $dump=$cl->{$com}{0}{$type};
      
      if (!$dump){return %R;}
      
      my %ref=xml2tag_list ($dump, "file");
      
      for (my $i=0; $i<$ref{n};$i++)
	{
	  
	  my $stream=xmltag2value($ref{$i}{body},"stream");
	  my $name=xmltag2value($ref{$i}{body},"name");
	  my $content=xmltag2value($ref{$i}{body},"content");
	
	if ($name  eq "stdout" || $name  eq "stderr")
	  {
	    $R{$name}=$content;
	  }
	else
	  { 
	   $R{file}{$name}{stream}=$stream;
	   $R{file}{$name}{content}=$content;
	   $R{nfiles}++;
	 }
	}
      $R{stack}=dump2stack ($dump);
      my $stderr=$R{stderr};
      my $stdout=$R{stdout};
    
      if ($stderr=~/ERROR/ || $stderr =~/FATAL/)
	{
	  $R{error}=1;
	}
      $R{warning}=($stderr=~/WARNING/g);
    
      return %R
  }
  
sub dump2stack
    {
      my ($file)=@_;
      
      my %cl=xml2tag_list ($file, "stack");
     
      my $stack= $cl{0}{body};

      if (!$stack){return "";}
      return $stack;
      
    }
sub wraptext 
    {
    my ($text, $wrap)=@_;
    
    $text=~s/$RETURN/\n/g;

    $text=~s/^\s*//;
    $text=~s/\s*$//;
      

    my @list=split ( /\n/, $text);
    my $ret;
    foreach my $l (@list)
      {
	$ret.="$wrap$l\n";
      }
    return $ret;
  }

sub dump2cl
   {
     my ($file)=@_;
     my %cl=xml2tag_list ($file, "cl");
     return clean_command ($cl{0}{body});
   }
	
   
sub rst2cl
    {
      my ($rst, $cl)=@_;
      my $lineN;
      my $clineN;
      my $wrap=0;
      my $command;
      
      open (F, $rst);
      
      while (<F>)
	{
	  my $line=$_;
	  $lineN++;
	  chomp($line);
	  if (!($line=~/\S/)){;}
	  elsif ($wrap)
	    {
	      if ($line=~/\s+(.*)\\/)
		{
		  $command.=$1;
		  $wrap=1;
		}
	      else 
		{
		  $line=~/\s+(\S.*)/;
		  $command.=$1;
		  $wrap=0;
		}
	      
	    }
	  elsif ($line=~/\s+$pattern: (.*)\\/)
	    {
	      $command=$1;
	      $clineN=$lineN;
	      $wrap=1;
	    }
	  elsif ($line=~/\s+$pattern: (.*)/)
	    {
	      $command=$1;
	      $clineN=$lineN;
	      $wrap=0;
	    }

	  if ($command && !$wrap)
	    {
	      $command=clean_command($command);
	      
	      
	      my $n=$cl->{$command}{n}++;
	      $cl->{$command}{$n}{doc}=$command;
	      $cl->{$command}{$n}{file}=$rst;
	      $cl->{$command}{$n}{line}=$clineN;
	      $command="";
	    }
	}
      close (F);
    }
sub substitute_cl4dump
      {
	my ($dump, $com)=@_;
	
	if (!-e $dump){return;}
	open (I, "$dump");
	open (O, ">$dump.tmp");
	while (<I>)
	  {
	    my $line=$_;
	    if ($line=~/\<cl\>/){print O "<cl>$com</cl>\n";}
	    else {print O "$line";}
	  }
	system ("mv $dump.tmp $dump");
      }
	
sub clean_command
      {
	my $c=shift;
	$c=~s/\s+/ /g;
	$c=~s/^\s+//g;
	$c=~s/\s+$//g;
	return $c;
      }
sub display_command
      {
	my $cl=shift;

	
	foreach my $com (keys %{$cl})
	  {
	    for (my $a=0; $a< $cl{$com}{n}; $a++)
	     {
	       print "#D [$com] -- $cl->{$com}{$a}{file} -- $cl->{$com}{$a}{line}\n";
	     }
	    if ($cl{$com}{0}{ref})
	      {
		print "#R [$com] $cl->{$com}{0}{ref}\n";
	      }
	  }
      }
#xml parsing

sub xmltag2value
  {
    my ($string_in, $tag)=@_;
    my %TAG;
    %TAG=xml2tag_list ($string_in, $tag);
    return $TAG{0}{body};
  }

sub xml2tag_list
  {
    my ($string_in,$tag)=@_;
    my ($tag_in, $tag_out, $string, $tag_in1, $tag_in2);
    my (@l, $in, $n, $t);
    my %tag;
    
    if (-e $string_in)
      {
	$string=&file2string ($string_in);
      }
    else
      {
	$string=$string_in;
      }
   
    
    $tag_in1="<$tag ";
    $tag_in2="<$tag>";
    $tag_out="/$tag>";
    $string=~s/>/>##1/g;
    $string=~s/</##2</g;
    $string=~s/##1/<#/g;
    $string=~s/##2/#>/g;
    @l=($string=~/(\<[^>]+\>)/g);
    $tag{n}=0;
    $in=0;$n=-1;



    foreach $t (@l)
      {
	
	$t=~s/<#//;
	$t=~s/#>//;

	if ( $t=~/$tag_in1/ || $t=~/$tag_in2/)
	  {

	    $in=1;
	    $tag{$tag{n}}{open}=$t;
	    $n++;

	  }
	elsif ($t=~/$tag_out/)
	  {


	    $tag{$tag{n}}{close}=$t;
	    $tag{n}++;
	    $in=0;
	  }
	elsif ($in)
	  {

	    $tag{$tag{n}}{body}.=$t;
	  }
      }

    return %tag;
  }


sub file_contains
  {
    my ($file, $tag, $max)=(@_);
    my ($n);
    $n=0;

    if ( !-e $file && ($file =~/$tag/)) {return 1;}
    elsif ( !-e $file){return 0;}
    else
      {
	open (FC, "$file");
	while ( <FC>)
	  {
	    if ( ($_=~/$tag/))
	      {
		close (FC);
		return 1;
	      }
	    elsif ($max && $n>$max)
	      {
		close (FC);
		return 0;
	      }
	    $n++;
	  }
      }
    close (FC);
    return 0;
  }


    
sub file2string
  {
    my $f=@_[0];
    my ($string, $l);
    open (F,"$f");
    while (<F>)
      {

	$l=$_;
	#chomp ($l);
	$string.=$l;
      }
    close (F);
    $string=~s/\r\n//g;
    $string=~s/\n/$RETURN/g;
    return $string;
  }


sub tag2value
  {

    my $tag=(@_[0]);
    my $word=(@_[1]);
    my $return;

    $tag=~/$word="([^"]+)"/;
    $return=$1;
    return $return;
  }
sub clean_cl
    {
      my @argl=@_;
      my $argv;
      foreach my $a (@argl)
	{
	  $a=~s/ /###SPACE###/g;
	  $argv.="$a ";
	}
      
      $argv=~s/[=,;]/ /g;
      @argl=split (/\s+/, $argv);
      for (my $a=0; $a<=$#argl; $a++){$argl[$a]=~s/###SPACE###/ /g;}
      
      return @argl;
    }

# Git functions

sub myunlink
      {
	my $f=shift;
	
	if (!-e $f){return;}
	unlink ($f);
	mygit_rm ($f);
      }
sub mymove
	{
	  my ($from, $to)=@_;
	  if (!-e $from)
	    {return;}
	  move ($from, $to);
	  mygit_rm ($from);
	  if (-d $to){$to="$to/$from";}
	  
	  if (-e $to)
	    {
	      mygit_add($to);
	    }
	}
	
sub file_is_tracked 
      {
	my $file=shift;
	if (!$GIT){return 0;}
	my $shell=system ("git ls-files $file --error-unmatch>/dev/null 2>/dev/null");
	if ($shell){return 0;}
	return 1;
      }
sub mygit_add
	{
	my @list1=@_;

	foreach my $f1 (@list1)
	  {
	    if (-d $f1)
	      {
		my @list2=dir2file_list($f1);
		foreach my $f2 (@list2)
		  {
		    mygit_add ("$f1/$f2");
		  }
	      }
	    else
	      {
		if ($GIT && !file_is_tracked ($f1)){mygit("add $f1");}
	      }
	  }
      }
sub mygit_rm
      {
	my @list1=@_;
	foreach my $f1(@list1)
	  {
	    if (-d $f1)
	      {
		my @list2=dir2file_list($f1);
		foreach my $f2 (@list2)
		  {
		    mygit_rm ("$f1/$f2");
		  }
	      }
	    else
	      {
		if ($GIT && (file_is_tracked ($f1))){mygit("rm -q -f $f1");}
	      }
	  }
      }
      
 sub mygit
      {
	my $arg=shift;
	my $com="git $arg";
	

	if (!$GIT)
	  {
	    return;
	  }
	
	elsif ($GIT==2)
	  {
	    print "###### DEBUG GIT ----------- $com\n";
	    return;
	  }
	else
	  {
	    return system ($com);
	  }
      }

sub dir2file_list
     {
       my ($cd, $pattern)=@_;
       my (@l, @nl);

       if (!-d $cd){return;}
       opendir (DIR, $cd);
       @l=readdir (DIR);
       closedir (DIR);
       
       foreach my $f (@l)
	 {
	   if ($f ne "." && $f ne "..")
	     {
	       if ($pattern && $f=~/$pattern/){push (@nl,$f);}
	       else {{push (@nl,$f);}}
	     }
	 }
       return @nl;
     }
     
