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
use DateTime;
use File::Basename;


my $RETURN=" ####RETURN#### ";
my %cl;
my %dir;
#$GIT=0: no git interaction
#$GIT=1: Commit all new files, uncommit all deleted files
#$GIT=2: blank run

my $TIMEOUT=1000;#Default Timeout in seconds
my $KEEPREPLAYED=0;
my $GIT=0;
my %FILE2IGNORE;

my $max=0;
my $pattern='$$';
my $log ;
my $docslog;
my $regexp;
my $reset;
my $stop_on_failed;
my $clean;
my $replay_dump_file;
my $replay_dump_list;
my $replay_dump_dir;
my $STRICT=0;
my $VERY_STRICT=0;
my $failed;
my $rep;
my $cw=cwd();
my $mode="new";#will only run the new ones
my $pg="t_coffee";

$dir{examples}="$cw/examples/";
$dir{docs}    ="$cw/docs/";
$dir{tmp}     ="$cw/testsuite/validation/docs/tmp/";
$dir{ref}     ="$cw/testsuite/validation/docs/ref/";
$dir{latest}  ="$cw/testsuite/validation/docs/latest/";
$dir{log}     ="$cw/testsuite/validation/docs/log/";
$dir{failed}  ="$cw/testsuite/validation/docs/failed/";

$FILE2IGNORE{'stdout'}=1;
$FILE2IGNORE{'stderr'}=1;


if ($ARGV[0] eq "-help")
  {
    print "docs2test.pl\n";
    print "Automaticly checks t_coffee command lines\n";
    print "The github dir structure is expected by default\n";
    print "tcoffee/\n";
    print "       /docs     -> contains rst docs\n";
    print "       /examples -> contains the reference files\n";
    print "       /testsuite/validation/docs/\n"; 
    print "       /testsuite/validation/docs/tmp    -> computation\n";
    print "       /testsuite/validation/docs/ref    -> succesful dumps\n";
    print "       /testsuite/validation/docs/failed -> unsuccesful dumps\n";
    print "\n";
    print "Commands are extracted from the .rst files contained in <-docs>\n";
    print "Commands are recognised as any line starting with <-pattern>\n";
    print "Duplicated commands are checked only once\n";
    print "By default the program only checks the new commands (-mode=new) \n";
    print "To check All the commands against the references, use -mode validate\n";
    print "Dumps are containers containing the CL and the input files\n";
    print "flags:\n";
    print "     -pattern       pattern used to recognize the command lines [def=**]\n";
    print "                    pattern will be treated as a regexp if -regexp is set\n";
    print "     -regexp        flag that causes pattern to be treated as a perl regexp\n";
    print "     -pg            specify the path of the version of T-Coffee (optional)\n";
    print "     -log           default: validation.log\n";
    print "     -docslog        default: docs.log\n";
    print "     -mode=<action> new|update|failed|check\n";
    print "                    new    : check ONLY CL w/o ref/dump and create ref/dump\n";
    print "                    update : check ALL  CL or create new ref/dum\n";
    print "                    failed : run   ONLY FAILURE as found /failed\n";
    print "                    check  : run   from dumps in ref";
    
    print "     -reset         delete all the dumps [CAUTION]\n";
    print "     -stop          stop at every FAILURE\n";
    print "     -clean         examples|dumps\n";
    print "                    examples: removes files in /examples/ not used by /docs [CAUTION]\n";
    print "                    dump:     removes files dumps not used by /docs [CAUTION]\n";
    
    print "     -rep           specifies a root repository\n";
    print "     -example       directory containing the sample files\n"; 
    print "     -docs          directory containing the .rst files\n";
    print "                    OR .rst file\n";
    print "                    OR file containing CLs (one per line)\n";
    print "     -ref           directory containing the reference dumps\n";
    print "                    OR dump file\n";
    print "     -failed        directory containing all the failed dumps\n";
    print "     -tmp           tmp directory\n";
    print "     -latest        latest directory\n";
    print "     -replay <file> replay an existing dump -- will produce file.replay along with diagnostic\n";
    print "     -replay_list<file> replay dumps in a file (one per line, comments start with #, files are checked)\n";
    print "     -replay_dir <dir>  replay dumps in a dir -- files are checked --\n";
    
    
    print "     -keepreplayed  Keep the dump of the replayed dump. Will be named file.replayed\n";
    print "     -strict        Will report failure if one or more replay output files are missing\n";
    print "     -very_strict   Will report failure if there is any difference between replay output\n";
    print "     -timeout       Will report failure is time is over this value [Def=$TIMEOUT sec.]\n";
    print "     -ignore        List of filesto be ignore: File1 File2 Def: -ignore stdout stderr\n";
    
    
    print "     -max           max number of CL to check [DEBUG]\n";
    print "     -help          display this help message\n";
    
    
    
    print "\n";
    print "\n";
    die;
    }

@ARGV=clean_cl(@ARGV);

for (my $a=0; $a<=$#ARGV; $a++)
  {
    
    if ($ARGV[$a]=~/-pattern/)
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
    elsif ($ARGV[$a]=~/-stop/)
      {
	$stop_on_failed=1;
      }
    elsif ($ARGV[$a]=~/-clean/)
      {
	$clean=$ARGV[++$a];;
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
    elsif ($ARGV[$a]=~/-docslog/)
      {
	$docslog=$ARGV[++$a];
      }
    elsif ($ARGV[$a]=~/-example/)
      {
	$dir{examples}=$ARGV[++$a];
      }
    elsif ($ARGV[$a]=~/-failed/)
      {
	$dir{failed}=$ARGV[++$a];
      }
    elsif ($ARGV[$a]=~/-docs/)
      {
	$dir{docs}=$ARGV[++$a];
      }
    elsif ($ARGV[$a]=~/-tmp/)
      {
	$dir{tmp}=$ARGV[++$a];
      }
    elsif ($ARGV[$a]=~/-ref/)
      {
	$dir{ref}=$ARGV[++$a];
      }
     elsif ($ARGV[$a]=~/-replay_list/)
      {
	$replay_dump_list=$ARGV[++$a];
      }
    elsif ($ARGV[$a]=~/-replay_dir/)
      {
	$replay_dump_dir=$ARGV[++$a];
      }
    
    elsif ($ARGV[$a]=~/-replay/)
      {
	$replay_dump_file=$ARGV[++$a];
      }
    elsif ($ARGV[$a]=~/-keepreplayed/)
      {
	$KEEPREPLAYED=1;
      }
    elsif ($ARGV[$a]=~/-rep/)
      {
	$rep=$ARGV[++$a];
      }
    elsif ($ARGV[$a]=~/-mode/)
      {
	$mode=$ARGV[++$a];
      }
    elsif ($ARGV[$a]=~/-timeout/)
      {
	$TIMEOUT=$ARGV[++$a];
      }
    elsif ($ARGV[$a]=~/-very_strict/)
      {
	$VERY_STRICT=1;
       
      }
    elsif ($ARGV[$a]=~/-strict/)
      {
	$STRICT=1;
      }
    elsif ($ARGV[$a]=~/-ignore/)
      {
	$FILE2IGNORE{$ARGV[++$a]}=1;
	while (!($ARGV[$a+1]=~/^-/))
	  {
	   $FILE2IGNORE{$ARGV[++$a]}=1;
	 }
      }
  }

# SIngle check mode: replay
if ($replay_dump_file)
  {
    exit(replay_dump_file ($replay_dump_file, 1));
  }
elsif ($replay_dump_list)
  {
    exit(replay_dump_list ($replay_dump_list));
  }
elsif ($replay_dump_dir)
  {
    exit(replay_dump_dir ($replay_dump_dir));
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
if (!$docslog){$docslog="$dir{log}/docs.log";}


#1-minimum checks



# 2 - create the directory structure if needed
foreach my $d (keys (%dir))
  {
    if (!-d $dir{$d} && !-e $dir{$d} && $d ne "docs")
      {
	print "create $dir{$d}\n"; 
	system ("mkdir -p $dir{$d}");	
      }
  }
# turn relative path into absolute paths
foreach my $d (keys (%dir))
  {
    $dir{$d}=rpath2apath($dir{$d});
  }

   
my $examplesD=$dir{examples};
my $docsD=$dir{docs};
my $tmpD=$dir{tmp};
my $refD=$dir{ref};
my $latestD=$dir{latest};
my $logD=$dir{log};
my $failedD=$dir{failed};

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

if ($mode eq "failed" || $mode eq "check")
  {
    $docsD=$dir{docs}="";
    $examplesD=$dir{examples}="";
  }


if ($docsD)
  {
    check_dir ($examplesD);
    
    if (-d $docsD && isrst ($docsD))
       {
	 my @list;
	 if (-d $docsD)
	   {
	     @list=dir2file_list($docsD, "\.rst\$");
	     @list=add_before_list ("$docsD/", @list);
	   }
	 else
	   {
	     @list=($docsD);
	   }
	 foreach my $f (@list)
	   {
	     rst2cl($f, \%cl);
	   }
       }
    elsif (-e $docsD)
      {
	my $lineN;
	open (F, $dir{docs});
	while(<F>)
	  {
	    my $l=$_;
	    if (!$l=~/\#\#/)
	      {
		$l=~/\s*(\S.*\S)\s*/;
		my $line=$1;
		
		if (($line=~/\S/))
		  {
		    my $n=$cl{$line}{docs}{n}++;
		    $cl{$line}{docs} {$n}{command}=$line;
		    $cl{$line}{docs} {$n}{source}=$docsD;
		    $cl{$line}{docs} {$n}{line}=$lineN;
		  }
		$lineN++;
	      }
	  }
      }
    else
      {
	#param is a CL
	$cl{$docsD}{docs}{n}++;
	$cl{$docsD}{docs}{0}{command}=$docsD;
      }
  }
if ($refD)
  {
    my @list;
    if (-d $refD)
      {
	my @l1=dir2file_list($refD);
	foreach my $dump (@l1)
	  {
	    my $dump="$refD/$dump";
	    if (isdump ($dump))
	      {
		push @list, $dump;
	      }
	  }
      }
    elsif (isdump ($refD))
      {
	push @list, $refD;
      }
    elsif (-e $refD && !isdump($refD))
      {
	print "ERROR --- $refD is not a valid dumpfile\n";
	die;
      }
    foreach my $d (@list)
      {
	refdump2cl($d, \%cl, "ref");
      }
  }

if ($failedD)
  {
    my @list;
    if (-d $failedD)
      {
	my @l1=dir2file_list($failedD);
	foreach my $dump (@l1)
	  {
	    my $dump="$failedD/$dump";
	    if (isdump ($dump))
	      {
		push @list, $dump;
	      }
	  }
      }
    elsif (isdump ($failedD))
      {
	push @list, $failedD;
      }
    elsif (-e $failedD  && !isdump($failedD))
      {
	print "ERROR --- $failedD is not a valid dumpfile\n";
	die;
      }
    foreach my $d (@list)
      {
	refdump2cl($d, \%cl, "failed");
      }
  }


#dump the documentation
docdump (\%cl, $docslog);


# 3 - get the examples in the tmp dir
if ($clean eq "examples")
  {
    clean_examples ($dir{examples}, $docslog);
    mygit_add ($dir{examples});
    exit;
  }

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
my  @com_list=keys (%cl);
my $n_com=$#com_list+1;
my $data;
my ($tot_failure, $tot_success, $tot_warning, $tot_tot);
$tot_failure=$tot_warning=$tot_success=$tot_tot=0;
my $dump_num;


my $com_index;
my $synced;
foreach my $com (sort (keys (%cl)))
  {
    my $shell_com;
    my $run;
    my $shell;
    my $etime;
    my $success=0;
    
    my $ref   =$cl{$com}{ref}{n};
    my $rdump =$cl{$com}{ref}{0}{file};
    my $doc   =$cl{$com}{doc}{n};
    my $failed=$cl{$com}{failed}{n};
    my $fdump =$cl{$com}{failed}{0}{file};
    my $ndump;
    my %R=dump2report($rdump);
    my %D;
    
    
    
    my %cl2;
    my $n2=0;
    foreach my $l ("doc", "failed", "ref")
      {
	for (my $n=0; $n<$cl{$com}{$l}{n}; $n++)
	  {
	    $n2=$cl2{n}++;
	    $cl2{$n2}{source}=$cl{$com}{$l}{$n}{source};
	    $cl2{$n2}{line}=$cl{$com}{$l}{$n}{line};
	  }
      }
    

    $com_index++;
    if ($rdump)
      {
	$ndump=path2name($rdump);
      }
    else
      {
	my $r=random_string (10);
	$ndump="t_coffee.$com_index.$r.dump";
      }
    
    
    
    $shell_com=$com;
    if ($pg ne "t_coffee"){$shell_com=~s/t_coffee/$pg/g;}
    
    
    $run=1;
    if (($mode eq "new" && !$ref) || ($mode eq "update" && !$doc))
      {
	if (!$synced){$synced=syncfiles($examplesD, $tmpD);}
	$shell=system4tc ("export DUMP_4_TCOFFEE=$ndump;$shell_com >/dev/null 2>/dev/null");
	%D=dump2report($ndump);
	if ($shell || $D{error} || $D{MissingOutput})
	  {
	    print "### Refresh TMP\n";
	    syncfiles ($examplesD, $tmpD);
	    $shell=system4tc ("export DUMP_4_TCOFFEE=$ndump;$shell_com >/dev/null 2>/dev/null");
	  }
      }
    elsif ($mode eq "check")
      {
	($shell,$etime)=dump2run ($rdump, $ndump, $shell_com);
      }
    elsif ($mode eq "failed" && $failed)
      {
	($shell, $etime)=dump2run ($fdump, $ndump, $shell_com);
      }
    else
      {
	$run=0;
      }
	
    if ($run)
      {
	$tot_tot++;
	$dump_num++;
	$cl{$com}{run}=1;
	if ($max && $tot_tot>=$max){last;}
	
	%D=dump2report($rdump);
	compare_reports (\%R, \%D);
	
	if (!-e $ndump)
	  {
	    open (F, ">$ndump");
	    print F "<DumpIO><cl>$com</cl><stack>#FAILURE -- SHELL ERROR -- Check the Command Line Syntax</stack><file><stream>output</stream><name>stderr</name><content>#FAILURE -- SHELL ERROR -- Check the Command Line Syntax</content></file></DumpIO>\n";
	    close (F);
	  }
	else
	  {
	    substitute_cl4dump($ndump,$com);
	  }
	
	$report="";	
	if ($D{error} || $D{MissingOutput})
	  {
	    $tot_failure++;
	    $cl{$com}{status}="failure";
	    for (my $a=0; $a< $cl2{n}; $a++)
	      {
		$report.="##FAILURE [$dump_num] -- COM -- $com -- Manual: $cl2{$a}{source} -- Line: $cl2{$a}{line} -- $ndump\n";
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
	    for (my $a=0; $a< $cl2{n}; $a++)
	      {
		$report.="##FAILURE [$dump_num] -- COM -- $com -- Manual: $cl2{$a}{source} -- Line: $cl{$a}{line} -- dump: $ndump\n";
	      }
	    $report .="##FAILURE [$dump_num] -- SHELL ERROR -- Check the Command Line Syntax\n";
	    
	  }
	elsif ($D{warning})
	  {
	    $tot_warning++;
	    $cl{$com}{status}="warning";
	    $success=1;
	    for (my $a=0; $a< $cl2{n}; $a++)
	      {
		$report.="##WARNING [$dump_num] -- COM -- $com -- Manual: $cl2{$a}{source} -- Line: $cl2{$a}{line} -- dump: $ndump\n";
	      }
	  }
	else
	  {
	    $tot_success++;
	    $cl{$com}{status}="success";
	    $success=1;
	    for (my $a=0; $a< $cl2{n}; $a++)
	      {
		$report.="##SUCCESS [$dump_num] -- COM -- $com -- Manual: $cl2{$a}{source} -- Line: $cl2{$a}{line} -- dump: $ndump\n";
		open (F, ">>$ndump");
		print F "<ref>\n<source>$cl2{$a}{source}</source>\n<line>$cl2{$a}{line}</line>\n</ref>\n";
		close (F);
	      }
	  }
	
	print "##$com_index Out of $n_com\n$report";
	open (LOG, ">>$log");print (LOG "$report");close (LOG);
	
	
	

	if ($success)
	  {
	    purgedir4dump($dir{"failed"}, $com);
	    if    ($mode eq "new")   {mymove ($ndump,$dir{ref});}
	    elsif ($mode eq "update")
	      {
		if (!$ref){mymove ($ndump,$dir{ref});}
		else {mymove ($ndump,$dir{latest});}
	      }
	    elsif ($mode eq "failed")
	      {
		mymove ($ndump,$dir{ref});
		
	      }
	    elsif ($mode eq "check")
	      {
		unlink ($ndump);
	      }
	  }
	elsif (!$success)
	  {
	    if (!$failed){mymove ($ndump,$dir{failed});}
	    if ( $stop_on_failed){exit;}
	  }
      }
  }


#Finalize the report

$report="## SUMMARY -- COMPLETE -- TOT: $tot_tot FAILURE: $tot_failure WARNING: $tot_warning SUCCESS: $tot_success\n";

foreach my $rst (cl2source (\%cl))
  {
    my ($failure, $success, $warning, $tot);
    $tot=$failure=$warning=$success=0;
    foreach my $com (keys (%cl))
      {
	if ($cl{$com}{run})
	  {
	    for (my $a=0; $a<$cl{$com}{n}; $a++)
	      {
		if ($cl{$com}{$a}{source} eq $rst)
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
  }
print "$report";
open (LOG, ">>$log");print (LOG "$report");close (LOG);
if ($clean=~/dump/)
  {
    purge(\%cl, \%dir);
  }
elsif ($clean=~/ref/)
  {

    foreach my $com (sort (keys (%cl)))
      {
	if (($cl{$com}{ref}{n} && !$cl{$com}{docs}{n}) )
	  {

	    for (my $i=0; $i<$cl{$com}{ref}{n}; $i++)
	      {
		my $dump=$cl{$com}{ref}{$i}{file};
		myunlink ($dump);
	      }
	  }
      }
  }
elsif ($clean=~/failed/)
  {

    foreach my $com (sort (keys (%cl)))
      {
	if (($cl{$com}{ref}{n} && $cl{$com}{failed}{n}) || (!$cl{$com}{docs}{n}) )
	  {

	    for (my $i=0; $i<$cl{$com}{failed}{n}; $i++)
	      {
		my $dump=$cl{$com}{failed}{$i}{file};
		myunlink ($dump);
	      }
	  }
      }
  }
exit;

sub replay_dump_dir
  {
    my ($replay_dir)=@_;
    my $shell=0;
    my $n=0;
    
    opendir (DIR, $replay_dir);
    my @dump_list=readdir (DIR);
    foreach my $dump (@dump_list)
      {
	if (file_is_dump($dump))
	  {
	    if (replay_dump_file ($dump, ++$n)){$shell=1;}
	  }
      }
    return $shell;
  }

sub replay_dump_list
  {
    my ($replay_list)=@_;
    my $n=0;
    my $shell=0;
    
    open (F, "$replay_list");
    while (<F>)
      {
	my $l=$_;
	chomp ($l);
	if (!$l=~/^#/ && file_is_dump($l))
	  {
	    if (replay_dump_file ($l, ++$n)){$shell=1;}
	  }
      }
    close (F);
    return $shell;
  }

      
sub replay_dump_file
  {
    my ($replay, $n)=@_;
    my $replayed=$replay.".replay";
    my ($shell,$etime)=dump2run ($replay, $replayed, "quiet");
    my $com=dump2cl ($replay);
    my ($missing, $different, $error, $warning);
    
    $missing=$warning=$error=$different=0;
    print "~ ($n) $com $etime ms ";
    my %in =dump2report($replay);
    my %out=dump2report($replayed);
    
    compare_reports (\%in, \%out, "quiet");
    $missing=$out{MissingOutput};
    $different=$out{N_DifferentOutput};
    $error=$out{error};
    
    if (!$error){$error=0;}
    if (!$warning){$warning=0;}
    if (!$different){$different=0;}
    if (!$missing){$missing=0;}
    
    
    print "MISSING_IO $missing ";
    print "DIFFERENT_IO $different ";
    print "WARNINGS $warning ";
    print "ERRORS $error ";
    
    if    ($STRICT && ($error || $missing)){$shell=1;}
    elsif ($VERY_STRICT && ($error || $missing || $warning || $different)){$shell=1;}
    
    if (!$shell)
      {
	print "PASSED\n";
      }
    else
      {
	print "FAILED\n";
      }
    if ($KEEPREPLAYED)
      {
	print "Replay File: $replayed\n";
      }
    else
      {
	unlink ($replayed);
      }
    return ($SHELL);
  }
  
#End of single replay      

sub docdump
  {
    my ($cl, $docdump)=@_;
    my %new;
    
    
    foreach my $c (keys (%{$cl}))
      {
	for (my $a=0; $a<$cl->{$c}{doc}{n}; $a++)
	  {
	    my $rst= $cl->{$c}{doc}{$a}{source};
	    my $line=$cl->{$c}{doc}{$a}{line};
	    $new{$rst}{$line}=$c;
	  }
      }
    open (F, ">$docdump");
    foreach my $file (keys (%new))
      {
	print F "## File: $file\n";
	foreach my $l (sort {$a <=> $b} (keys(%{$new{$file}})))
	  {
	    printf F "## Line: %-4d --- \n$new{$file}{$l}\n", $l;
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
    my $npurged;
    
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
		if (!$cl->{$com}{doc}{n})
		  {
		    print "-----  Purge $dump\n";
		    myunlink ($dump);
		    $npurged++;
		  }
	      }
	  }
      }
    return $npurged;
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
    my ($ref, $doc, $quiet)=@_;

    
    if (!$ref || !$doc){return;}
    
    foreach my $f (keys (%{$ref->{file}}))
      {
	if ($quiet ne "quiet")
	  {print "FILE: $f";
	   print "********\n\n";
	   print "($doc->{file}{$f}{content}\n\n\n\n($ref->{file}{$f}{content}\n\n\n";
	   print "********\n\n";
	 }


	if (!$FILE2IGNORE{$f})
	  {
	    #Get rid of the Version and CPU effects
	    my $content1=$doc->{file}{$f}{content};
	    my $content2=$doc->{file}{$f}{content};
	    
	    $content1=~s/Version_\S+ /Version_XXXX /g;
	    $content2=~s/Version_\S+ /Version_XXXX /g;
	    
	    $content1=~s/CPU\S+ /CPU=XXXX /g;
	    $content2=~s/CPU\S+ /CPU=XXXX /g;
	    
	    if (!$doc->{file}{$f})
	      {
		$doc->{MissingOutputF}{$f}=1;$doc->{MissingOutput}++;
	      }
	    elsif ($content1 ne $content2)
	      {
		$doc->{DifferentOutput}{$f}=1;$doc->{N_DifferentOutput}++;
	      }
	  }
      }
    return;
  }
sub system4tc
    {
      my $com=shift;
      system ("t_coffee -clean >/dev/null 2>/dev/null");
      return timeout_system ($com);
    }

sub dump2run
  {
    my ($idump, $odump, $com)=@_;
    my $cdir=cwd;
    my $use_stdout;
    my $shell;
    
    $idump=path2abs($idump);
    $odump=path2abs($odump);
    
    
    my $dir=random_string();

    system ("mkdir $dir");
    chdir  ($dir);
    
    if (!$com){$com=dump2cl($idump);}
    $com=dump2cl($idump);
    
    if ($com =~/.*\|(.*)/)
      {
	$com=$1;
      }
    my %ref=xml2tag_list ($idump, "file");
    

    for (my $i=0; $i<$ref{n};$i++)
	{
	  
	  my $stream=xmltag2value($ref{$i}{body},"stream");
	  my $name=xmltag2value($ref{$i}{body},"name");
	  my $content=xmltag2value($ref{$i}{body},"content");
	  $content=~s/$RETURN/\n/g;
	  
	  
	  if ($stream eq "input")
	    {
	      my $dir=dirname ($name);
	      $dir="./$dir";
	      system ("mkdir -p $dir");
	      
	      open (F, ">$name");
	      print F "$content";
	      close (F);
	      
	      if ($name eq "stdin"){$com="cat stdin | $com";}
	    }
	  if ($stream eq "output" && $name eq "stdout"){$use_stdout=1;}
	}
    
   
    my $before=time;
    unlink ($odump);
    if (!$use_stdout)
      {
	$shell=system4tc ("export DUMP_4_TCOFFEE=$odump;$com >/dev/null 2>/dev/null");
      }
    else
      {
	$shell=system4tc ("export DUMP_4_TCOFFEE=$odump;$com >stdout 2>/dev/null");
      }
    my $etime=time - $before;
    $etime*=1000;
        
    chdir($cdir);
    system ("rm -rf $dir");
    
    return ($shell, $etime);
  }
  
sub dump2report
  {
      my ($dump)=shift;
      my $cdir=cwd;
      my %ref;
      if (!$dump || !-e $dump){return %ref;}
     
      %ref=xml2tag_list ($dump, "file");
      
      for (my $i=0; $i<$ref{n};$i++)
	{
	  my $stream=xmltag2value($ref{$i}{body},"stream");
	  my $name=xmltag2value($ref{$i}{body},"name");
	  my $content=xmltag2value($ref{$i}{body},"content");
	
	if ($name  eq "stdout" || $name  eq "stderr")
	  {
	    $ref{$name}=$content;
	  }
	else
	  { 
	   $ref{file}{$name}{stream}=$stream;
	   $ref{file}{$name}{content}=$content;
	   $ref{nfiles}++;
	 }
	}
      $ref{stack}=dump2stack ($dump);
      my $stderr=$ref{stderr};
      my $stdout=$ref{stdout};
    
      if ($stderr=~/ERROR/ || $stderr =~/FATAL/)
	{
	  $ref{error}=1;
	}
      $ref{warning}=($stderr=~/WARNING/g);
    
      return %ref
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

sub refdump2cl
  {
    my ($dump, $cl, $type)=@_;
    
    if (!-e $dump || !isdump($dump)){return;}
    my $com=dump2cl($dump);
    

    my %h=xml2tag_list ($dump, "ref");
    if (!$h{n}){$h{n}=1;}
    for (my $i=0; $i<$h{n};$i++)
	{
	  my $n=$cl->{$com}{$type}{n}++;
	  my $source=xmltag2value($h{$i}{body},"source");
	  my $line=xmltag2value($h{$i}{body},"line");
	  $cl->{$com}{$type}{$n}{command}=$com;
	  $cl->{$com}{$type}{$n}{file}=$dump;
	  $cl->{$com}{$type}{$n}{line}=$line;
	  $cl->{$com}{$type}{$n}{source}=$source;
	  
	}
    return $com;
  }
sub cl2source
    {
      my $cl=shift;
      my %h;

      foreach my $com (keys (%{$cl}))
	{
	  for (my $i=0; $i< $cl->{$com}{n}; $i++)
	    {
	      my $file=$cl->{$com}{$i}{source};
	      $h{$file}=1;
	    }
	}
      return keys (%h);
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
	      
	      
	      my $n=$cl->{$command}{doc}{n}++;
	      $cl->{$command}{doc}{$n}{command}=$command;
	      $cl->{$command}{doc}{$n}{source}   =$rst;
	      $cl->{$command}{doc}{$n}{line}   =$clineN;
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
	       print "#D [$com] -- $cl->{$com}{$a}{source} -- $cl->{$com}{$a}{line}\n";
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
sub xml2tag_list_test
 { 
   my ($string_in,$tag)=@_;
   my ($tag_in, $tag_out, $string, $tag_in1, $tag_in2);
   my (@l, $in, $n, $t);
   my %tag;
   my $r1=random_string (20);
   my $r2=random_string (20);
   my $r3=random_string (20);
   my $r4=random_string (20);
   

   my $tag_in="<$tag>";
   my $tag_out="<\/$tag>";
   my $tag_inR=" $r1 ";
   my $tag_outR=" $r2 ";
   my $openR=" $r3 ";
   my $closeR=" $r4 ";
   
   
   if (-e $string_in)
     {
       $string=&file2string ($string_in);
     }
   else
     {
       $string=$string_in;
     }
   
   $string=~s/$tag_in/$tag_inR/g;
   $string=~s/$tag_out/$tag_outR/g;
   $string=~s/\</$openR/g;
   $string=~s/\>/$closeR/g;
   $string=~s/$tag_inR/\</g;
   $string=~s/$tag_outR/\>/g;
   
   my @l=($string=~/\<([^>]+)\>/g);
   foreach my $e (@l)
     {
       $e=~s/$openR/\</g;
       $e=~s/$closeR/\>/g;
       $tag{$tag{n}++}{body}=$e;
     }
   return %tag;
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
    
    my $cwd=$cw;
    
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

sub file_isdump
    {
      my $f=shift;
      if (!-e $f) {return 0;}
      open (F, "$f");
      while (<F>)
	{
	  my $l=$_;
	  close (F);
	  return $l=~/DumpIO/;
	}
    }
sub isdump
    {
      my $f1=shift;
      if (-d $f1)
	{
	  my @l=dir2file_list ($f1);
	  foreach my $f2 (@l)
	    {
	      if (file_isdump ("$f1/$f2")){return 1;}
	    }
	  return 0;
	}
      else
	{
	  return file_isdump ($f1);
	}
    }
sub isrst
    {
      my $f1=shift;
      if (-d $f1)
	{
	  my @l=dir2file_list ($f1);
	  foreach my $f2 (@l)
	    {
	      if ($f2=~/\.rst$/) {return 1;}
	    }
	  return 0;
	}
      else
	{
	  return ($f1=~/\.rst$/);
	}
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
	       if ($pattern)
		 {
		   if ($f=~/$pattern/){push (@nl,$f);}
		 }
	       else {{push (@nl,$f);}}
	     }
	 }
       return @nl;
     }
     
sub random_string
       {
	 my $l=shift;

	 if (!$l){$l=20;}
	 my $ret;
	 my $s="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890";
	 my @l=split (//,$s);
	 for (my $a=0; $a<$l; $a++)
	   {
	     my $c=int(rand($#l+1));
	     $ret.=$l[$c];
	   }
	 return $ret;
       }
sub check_dir
       {
	 my $dir=shift;

	 if (!-d $dir && !-e $dir)
	   
	   {
	     print STDERR "ERROR: $dir could not be found. In default mode run the script from the repository root [FATAL]\n";
	     die;
	   }
	 return;
       }
     
sub add_before_list
       {
	 my ($add,@list)=@_;
	 for (my $a=0; $a<=$#list; $a++)
	   {
	     $list[$a]=$add.$list[$a];
	   }
	 return @list;
       }
sub add_after_list
      {
	my (@list, $add)=@_;
	for (my $a=0; $a<=$#list; $a++)
	  {
	    $list[$a]=$list[$a].$add;
	  }
	return @list;
      }	 
sub syncfiles
    {
      my ($from, $to)=@_;
      if (!-d $from || !-d $to){return;}
      my $n=0;
      
      system ("rm $to/*");
      my @list=dir2file_list($from);
      foreach my $f (@list)
	{
	  $f="$from/$f";
	  if ( -e $f && !-d $f)
	    {$n++;
	     copy ($f, $to);
	   }
	}
      return $n;
    }
sub path2name
      {
	my $f=shift;
	if (!($f=~/\//)){return $f;}
	else
	  {
	    $f =~/.*\/([^\/]+)/;
	    return $1;
	  }
      }
sub rpath2apath
      {
	my $f=shift;
	my $cw=cwd();
	
	if ($f=~/^\//)
	  {
	    return $f;
	  }
	elsif ($f=~/^\.\/(.*)/)
	  {
	    
	    return "$cw/$1";
	  }
	else
	  {
	    return "$cw/$f";
	  }
      }
	  

sub path2abs
	{
	  my ($file)=@_;
	  if ($file=~/^^\//){return $file;}
	  my $dir=cwd;
	  $file=$dir."/".$file;
	  return $file;
	}



sub timeout_system
   {
     my ($command, $timeout)=(@_);
     my $shell;
     if (!$timeout)
       {
	 $timeout=$TIMEOUT;
       }
     
     eval 
       {
	 local $SIG{ALRM} = sub { die "alarm\n" }; # NB: \n required
	 
	 alarm($timeout);
	 $shell=system ($command);
	 alarm(0);
       };
     
     if ($@ eq "alarm\n"){return "-1";}
     else {return $shell}
   }
