#!/usr/bin/env perl
use Getopt::Long;
use File::Path;
use Env;
use FileHandle;
use Cwd;
use Sys::Hostname;
our $PIDCHILD;
our $ERROR_DONE;
our @TMPFILE_LIST;
our $EXIT_FAILURE=1;
our $EXIT_SUCCESS=0;

our $REFDIR=getcwd;
our $EXIT_SUCCESS=0;
our $EXIT_FAILURE=1;

our $PROGRAM="tc_generic_method.pl";
our $CL=$PROGRAM;

our $CLEAN_EXIT_STARTED;
our $debug_lock=$ENV{"DEBUG_LOCK"};
our $LOCKDIR=$ENV{"LOCKDIR_4_TCOFFEE"};
if (!$LOCKDIR){$LOCKDIR=getcwd();}
our $ERRORDIR=$ENV{"ERRORDIR_4_TCOFFEE"};
our $ERRORFILE=$ENV{"ERRORFILE_4_TCOFFEE"};
&set_lock ($$);
if (isshellpid(getppid())){lock4tc(getppid(), "LLOCK", "LSET", "$$\n");}
      
our $print;
my ($fmsq1, $fmsq2, $output, $outfile, $arch, $psv, $hmmtop_home, $trim, $cov, $sample, $mode, $gor_home, $gor_seq, $gor_obs);

GetOptions("-in=s" => \$fmsq1,"-output=s" =>\$output ,"-out=s" => \$outfile, "-arch=s" => \$arch,"-psv=s" => \$psv, "-hmmtop_home=s", \$hmmtop_home,"-trim=s" =>\$trim ,"-print=s" =>\$print,"-cov=s" =>\$cov , "-sample=s" =>\$sample, "-mode=s" =>\$mode, "-gor_home=s"=>\$gor_home, "-gor_seq=s"=>\$gor_seq,"-gor_obs=s"=>\$gor_obs);


if (!$mode){$mode = "hmmtop"}
elsif ($mode eq "hmmtop"){;}
elsif ($mode eq "gor"){;}
else {myexit(flush_error ("-mode=$mode is unknown"));}


our $HOME=$ENV{"HOME"};
our $MCOFFEE=($ENV{"MCOFFEE_4_TCOFFEE"})?$ENV{"MCOFFEE_4_TCOFFEE"}:"$HOME/.t_coffee/mcoffee";

if ($mode eq "hmmtop")
  {
    
    check_configuration ("hmmtop");
    if (-e $arch){$ENV{'HMMTOP_ARCH'}=$arch;}
    elsif (-e $ENV{HMMTOP_ARCH}){$arch=$ENV{HMMTOP_ARCH};}
    elsif (-e "$MCOFFEE/hmmtop.arch"){$arch=$ENV{'HMMTOP_ARCH'}="$MCOFFEE/hmmtop.arch";}
    elsif (-e "$hmmtop_home/hmmtop.arch"){$arch=$ENV{'HMMTOP_ARCH'}="$hmmtop_home/hmmtop.arch";}
    else {myexit(flush_error ( "Could not find ARCH file for hmmtop"));}
    
    
    if (-e $psv){$ENV{'HMMTOP_PSV'}=$psv;}
    elsif (-e $ENV{HMMTOP_PSV}){$psv=$ENV{HMMTOP_PSV};}
    elsif (-e "$MCOFFEE/hmmtop.psv"){$psv=$ENV{'HMMTOP_PSV'}="$MCOFFEE/hmmtop.psv";}
    elsif (-e "$hmmtop_home/hmmtop.psv"){$psv=$ENV{'HMMTOP_PSV'}="$hmmtop_home/hmmtop.psv";}
    else {myexit(flush_error ( "Could not find PSV file for hmmtop"));}

  }
elsif ($mode eq "gor")
  {
    our $GOR_SEQ;
    our $GOR_OBS;
    
    check_configuration ("gorIV");
    if (-e $gor_seq){$GOR_SEQ=$gor_seq;}
    elsif (-e $ENV{GOR_SEQ}){$GOR_SEQ=$ENV{GOR_SEQ};}
    elsif (-e "$MCOFFEE/New_KS.267.seq"){$GOR_SEQ="$MCOFFEE/New_KS.267.seq";}
    elsif (-e "$gor_home/New_KS.267.seq"){$GOR_SEQ="$gor_home/New_KS.267.seq";}
    else {myexit(flush_error ( "Could not find SEQ file for gor"));}

    if (-e $gor_obs){$GOR_OBS=$gor_obs;}
    elsif (-e $ENV{GOR_OBS}){$GOR_OBS=$ENV{GOR_OBS};}
    elsif (-e "$MCOFFEE/New_KS.267.obs"){$GOR_OBS="$MCOFFEE/New_KS.267.obs";}
    elsif (-e "$gor_home/New_KS.267.obs"){$GOR_OBS="$gor_home/New_KS.267.obs";}
    else {myexit(flush_error ( "Could not find OBS file for gor"));}
  }


if ( ! -e $fmsq1){myexit(flush_error ("Could Not Read Input file $fmsq1"));}


my $fmsq2=vtmpnam();
my $fmsq3=vtmpnam();
my $tmpfile=vtmpnam();
my $predfile=vtmpnam();

if ($trim){$trim_action=" +trim _aln_%%$trim\_K1 ";}
if ($cov) {$cov_action= " +sim_filter _aln_c$cov ";}
&safe_system("t_coffee -other_pg seq_reformat -in $fmsq1 -action +convert 'BOUJXZ-' $cov_action $trim_action -output fasta_aln -out $fmsq2");
my (%pred, %seq, %predA);


%seq=read_fasta_seq($fmsq2);
%seq=fasta2sample(\%seq, $sample);

if (1==2 &&$mode eq "hmmtop" && $output eq "cons")
  {
    fasta2hmmtop_cons($outfile,\%seq);
  }
else
  {
   
    %pred=fasta2pred(\%seq, $mode);
    %predA=pred2aln (\%pred, \%seq);
    
    
    if (!$output || $output eq "prediction"){output_fasta_seq (\%predA, $outfile);}
    elsif ($output eq "color_html"){pred2color (\%pred,\%seq, $outfile);}
    elsif ($output eq "cons"){pred2cons($outfile,\%predA);}
    else {flush_error ("$output is an unknown output mode");}
  }

sub fasta2sample
  {
    my $SR=shift;
    my $it=shift;
    my %S=%$SR;
    
    my $seq=index2seq_name (\%S, 1);
    my $l=length($S{$seq}{seq});
    my @sl=keys(%S);
    my $nseq=$#sl+1;
    my $index=$nseq;
  
    if (!$sample) {return %S;}
    for (my $a=0; $a<$it; $a++)
      {
	my $newseq="";
	my $nname="$seq\_sampled_$index";
	for (my $p=0; $p<$l; $p++)
	  {
	    my $i=int(rand($nseq));
	    
	    my $name = $sl[$i];
	    my $seq=$S{$name}{seq};
	    my $r=substr ($seq, $p, 1);
	    $newseq.=$r;
	  }
	$S{$nname}{name}=$nname;
	$S{$nname}{seq}=$newseq;
	$S{$nname}{com}="sampled";
	$S{$nname}{index}=++$index;
      }
    return %S;
  }
	      
sub fasta2pred
  {
    my $s=shift;
    my $mode=shift;

    if ( $mode eq "hmmtop"){return fasta2hmmtop_pred($s);}
    elsif ($mode eq "gor"){return fasta2gor_pred ($s);}
  }
sub fasta2hmmtop_cons
  {
    my $outfile=shift;
    my $SR=shift;
    
    my $o = new FileHandle;
    my $i = new FileHandle;
    my $tmp_in =vtmpnam();
    my $tmp_out=vtmpnam();
    my %seq=%$SR;
    my %pred;
    my $N=keys(%seq);
    
    output_fasta_seq (\%seq,$tmp_in, "seq");
    `hmmtop -pi=mpred -if=$tmp_in -sf=FAS -pl 2>/dev/null >$tmp_out`;
    open ($o, ">$outfile");
    open ($i, "$tmp_out");
    while (<$i>)
      {
	my $l=$_;
	if (($l=~/>HP\:\s+(\d+)\s+(.*)/)){my $line=">$2 NSEQ: $N\n";print $o "$line";}
	elsif ( ($l=~/.*pred(.*)/))  {my $line="$1\n";print $o "$line";}
      }
    close ($o);
    close ($i);
    return read_fasta_seq($tmp);
  }
sub fasta2hmmtop_pred
  {
    my $SR=shift;
    my $o = new FileHandle;
    my $i = new FileHandle;
    my $tmp    =vtmpnam();
    my $tmp_in =vtmpnam();
    my $tmp_out=vtmpnam();
    my %seq=%$SR;
    my %pred;
    

    output_fasta_seq (\%seq,$tmp_in, "seq");

    
    `hmmtop -if=$tmp_in -sf=FAS -pl 2>/dev/null >$tmp_out`;
    

    
    
    open ($o, ">$tmp");
    open ($i, "$tmp_out");
    while (<$i>)
      {
	my $l=$_;
	if (($l=~/>HP\:\s+(\d+)\s+(.*)/)){my $line=">$2\n";print $o "$line";}
	elsif ( ($l=~/.*pred(.*)/))  {my $line="$1\n";print $o "$line";}
      }
    close ($o);
    close ($i);
    return read_fasta_seq($tmp);
  }
    
	
	
	    
	
	

	
sub fasta2gor_pred
  {
    my $SR=shift;
    my $o = new FileHandle;
    my $i = new FileHandle;
    my $tmp    =vtmpnam();
    my $tmp_in =vtmpnam();
    my $tmp_out=vtmpnam();
    my %seq=%$SR;
    my %pred;
    

    output_fasta_seq (\%seq,$tmp_in, "seq");
    `gorIV -prd $tmp_in -seq $GOR_SEQ -obs $GOR_OBS >$tmp_out`;
    open ($o, ">$tmp");
    open ($i, "$tmp_out");
    while (<$i>)
      {
	my $l=$_;

	
	if ( $l=~/>/){print $o "$l";}
	elsif ( $l=~/Predicted Sec. Struct./){$l=~s/Predicted Sec. Struct\.//;print $o "$l";}
      }
    close ($o);
    close ($i);
    return read_fasta_seq($tmp);
  }
			
			     
sub index2seq_name
  {
    
    my $SR=shift;
    my $index=shift;
    
    
    my %S=%$SR;
    
    foreach my $s (%S)
      {
	if ( $S{$s}{index}==$index){return $s;}
      }
    return "";
  }

sub pred2cons
  {
    my $outfile=shift;
    my $predR=shift;
    my $seq=shift;
    my %P=%$predR;
    my %C;
    my ($s,@r,$nseq);
    my $f= new FileHandle;

    open ($f, ">$outfile");

    if (!$seq){$seq=index2seq_name(\%P,1);}
    foreach my $s (keys(%P))
      {
	$nseq++;
	$string= $P{$s}{seq};
	$string = uc $string;
	my @r=split (//,$string);
	for (my $a=0; $a<=$#r; $a++)
	  {
	    if (($r[$a]=~/[OHICE]/)){$C{$a}{$r[$a]}++;}
	  }
      }
    @l=keys(%C);
    
    
    $s=$P{$seq}{seq};
    print $f ">$seq pred based on $nseq\n";
    @r=split (//,$s);
    
    for (my $x=0; $x<=$#r; $x++)
      {
	if ($r[$x] ne "-")
	  {
	    my $h=$C{$x}{H};
	    my $i=$C{$x}{I};
	    my $o=$C{$x}{O};
	    my $c=$C{$x}{C};
	    my $e=$C{$x}{E};
	    my $l=$i+$o;
	    
	    if ($h>=$i && $h>=$o && $h>=$c && $h>=$e){$r[$x]='H';}
	    elsif ($i>=$o && $i>=$c && $i>=$e){$r[$x]='I';}
	    elsif ($o>=$c && $o>=$e){$r[$x]='O';}
	    elsif ($c>=$e){$r[$x]='C';}
	    else {$r[$x]='E';}
	  }
      }
    $j=join ('', @r);
    print $f "$j\n";
    close ($f);
    return $j;
  }

sub pred2aln
  {
    my $PR=shift;
    my $AR=shift;
    
    my $f=new FileHandle;
    my %P=%$PR;
    my %A=%$AR;
    my %PA;
    my $tmp=vtmpnam();
    my $f= new FileHandle;
    
    open ($f, ">$tmp");
    foreach my $s (sort{$A{$a}{index}<=>$A{$b}{index}}(keys (%A)))
      {
	my (@list, $seq, @plist, @pseq, $L, $PL, $c, $w);
	my $seq;
	my $seq=$A{$s}{seq};
	my $pred=$P{$s}{seq};
	$seq=pred2alnS($P{$s}{seq},$A{$s}{seq});
	print $f ">$s\n$seq\n";
      }
    close ($f);
    return read_fasta_seq ($tmp);
  }
sub pred2alnS
  {
    my $pred=shift;
    my $aln= shift;
    my ($j,$a,$b);
    my @P=split (//, $pred);
    my @A=split (//, $aln);
    for ($a=$b=0;$a<=$#A; $a++)
      {
	if ($A[$a] ne "-"){$A[$a]=$P[$b++];}
      }
    if ($b!= ($#P+1)){add_warning ("Could not thread sequence: $b $#P");}
    
    $j= join ('', @A);
    return $j;
  }
sub pred2color
  {
    my $predP=shift;
    my $alnP=shift;
    my $out=shift;
    my $F=new FileHandle;
    my $struc=vtmpnam();
    my $aln=vtmpnam();
    

    output_fasta_seq ($alnP, $aln);
    my %p=%$predP;
    
    open ($F, ">$struc");
    
    
    foreach my $s (keys(%p))
      {
	
	print $F ">$s\n";
	my $s=uc($p{$s}{seq});
	
	$s=~s/[Oo]/0/g;
	$s=~s/[Ee]/0/g;
	
	$s=~s/[Ii]/5/g;
	$s=~s/[Cc]/5/g;
	
	$s=~s/[Hh]/9/g;
	
	print $F "$s\n";
      }
    close ($F);
    
    
    
    safe_system ( "t_coffee -other_pg seq_reformat -in $aln -struc_in $struc -struc_in_f number_fasta -output color_html -out $out");
    return;
  }
	  
    
sub display_fasta_seq
  {
    my $SR=shift;
    my %S=%$SR;
    
    foreach my $s (sort{$S{$a}{index}<=>$S{$b}{index}}(keys (%S)))
      {
	print STDERR ">$s\n$S{$s}{seq}\n";
      }
    close ($f);
  }
sub output_fasta_seq
  {
    my $SR=shift;
    my $outfile=shift;
    my $mode =shift;
    my $f= new FileHandle;
    my %S=%$SR;
    
    
    open ($f, ">$outfile");
    foreach my $s (sort{$S{$a}{index}<=>$S{$b}{index}}(keys (%S)))
      {
	my $seq=$S{$s}{seq};
	if ( $mode eq "seq"){$seq=~s/\-//g;}
	print $f ">$s\n$seq\n";
      }
    close ($f);
  }
      
sub read_fasta_seq 
  {
    my $f=$_[0];
    my %hseq;
    my (@seq, @com, @name);
    my ($a, $s,$nseq);
    my $index;
    open (F, $f);
    while (<F>)
      {
	$s.=$_;
      }
    close (F);

    
    @name=($s=~/>(\S*).*\n[^>]*/g);
    
    @seq =($s=~/>.*.*\n([^>]*)/g);
    @com =($s=~/>.*(.*)\n([^>]*)/g);


    $nseq=$#name+1;
    
  
    for ($a=0; $a<$nseq; $a++)
      {
	my $n=$name[$a];
	my $s;
	$hseq{$n}{name}=$n;
	$s=$seq[$a];$s=~s/\s//g;
	$hseq{$n}{index}=++$index;
	$hseq{$n}{seq}=$s;
	$hseq{$n}{com}=$com[$a];
      }
    return %hseq;
  }


sub file2head
      {
	my $file = shift;
	my $size = shift;
	my $f= new FileHandle;
	my $line;
	open ($f,$file);
	read ($f,$line, $size);
	close ($f);
	return $line;
      }
sub file2tail
      {
	my $file = shift;
	my $size = shift;
	my $f= new FileHandle;
	my $line;
	
	open ($f,$file);
	seek ($f,$size*-1, 2);
	read ($f,$line, $size);
	close ($f);
	return $line;
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

sub myexit
  {
    my $code=@_[0];
    if ($CLEAN_EXIT_STARTED==1){return;}
    else {$CLEAN_EXIT_STARTED=1;}
    ### ONLY BARE EXIT
    exit ($code);
  }
sub set_error_lock
    {
      my $name = shift;
      my $pid=$$;

      
      &lock4tc ($$,"LERROR", "LSET", "$$ -- ERROR: $name $PROGRAM\n");
      return;
    }
sub set_lock
  {
    my $pid=shift;
    my $msg= shift;
    my $p=getppid();
    &lock4tc ($pid,"LLOCK","LRESET","$p$msg\n");
  }
sub unset_lock
   {
     
    my $pid=shift;
    &lock4tc ($pid,"LLOCK","LRELEASE","");
  }
sub shift_lock
  {
    my $from=shift;
    my $to=shift;
    my $from_type=shift;
    my $to_type=shift;
    my $action=shift;
    my $msg;
    
    if (!&lock4tc($from, $from_type, "LCHECK", "")){return 0;}
    $msg=&lock4tc ($from, $from_type, "LREAD", "");
    &lock4tc ($from, $from_type,"LRELEASE", $msg);
    &lock4tc ($to, $to_type, $action, $msg);
    return;
  }
sub isshellpid
  {
    my $p=shift;
    if (!lock4tc ($p, "LLOCK", "LCHECK")){return 0;}
    else
      {
	my $c=lock4tc($p, "LLOCK", "LREAD");
	if ( $c=~/-SHELL-/){return 1;}
      }
    return 0;
  }
sub isrootpid
  {
    if(lock4tc (getppid(), "LLOCK", "LCHECK")){return 0;}
    else {return 1;}
  }
sub lock4tc
	{
	  my ($pid,$type,$action,$value)=@_;
	  my $fname;
	  my $host=hostname;
	  
	  if ($type eq "LLOCK"){$fname="$LOCKDIR/.$pid.$host.lock4tcoffee";}
	  elsif ( $type eq "LERROR"){ $fname="$LOCKDIR/.$pid.$host.error4tcoffee";}
	  elsif ( $type eq "LWARNING"){ $fname="$LOCKDIR/.$pid.$host.warning4tcoffee";}
	  
	  if ($debug_lock)
	    {
	      print STDERR "\n\t---lock4tc(tcg): $action => $fname =>$value (RD: $LOCKDIR)\n";
	    }

	  if    ($action eq "LCHECK") {return -e $fname;}
	  elsif ($action eq "LREAD"){return file2string($fname);}
	  elsif ($action eq "LSET") {return string2file ($value, $fname, ">>");}
	  elsif ($action eq "LRESET") {return string2file ($value, $fname, ">");}
	  elsif ($action eq "LRELEASE") 
	    {
	      if ( $debug_lock)
		{
		  my $g=new FileHandle;
		  open ($g, ">>$fname");
		  print $g "\nDestroyed by $$\n";
		  close ($g);
		  safe_system ("mv $fname $fname.old");
		}
	      else
		{
		  unlink ($fname);
		}
	    }
	  return "";
	}
	
sub file2string
	{
	  my $file=@_[0];
	  my $f=new FileHandle;
	  my $r;
	  open ($f, "$file");
	  while (<$f>){$r.=$_;}
	  close ($f);
	  return $r;
	}
sub string2file 
    {
    my ($s,$file,$mode)=@_;
    my $f=new FileHandle;
    
    open ($f, "$mode$file");
    print $f  "$s";
    close ($f);
  }

BEGIN
    {
      srand;
    
      $SIG{'SIGUP'}='signal_cleanup';
      $SIG{'SIGINT'}='signal_cleanup';
      $SIG{'SIGQUIT'}='signal_cleanup';
      $SIG{'SIGILL'}='signal_cleanup';
      $SIG{'SIGTRAP'}='signal_cleanup';
      $SIG{'SIGABRT'}='signal_cleanup';
      $SIG{'SIGEMT'}='signal_cleanup';
      $SIG{'SIGFPE'}='signal_cleanup';
      
      $SIG{'SIGKILL'}='signal_cleanup';
      $SIG{'SIGPIPE'}='signal_cleanup';
      $SIG{'SIGSTOP'}='signal_cleanup';
      $SIG{'SIGTTIN'}='signal_cleanup';
      $SIG{'SIGXFSZ'}='signal_cleanup';
      $SIG{'SIGINFO'}='signal_cleanup';
      
      $SIG{'SIGBUS'}='signal_cleanup';
      $SIG{'SIGALRM'}='signal_cleanup';
      $SIG{'SIGTSTP'}='signal_cleanup';
      $SIG{'SIGTTOU'}='signal_cleanup';
      $SIG{'SIGVTALRM'}='signal_cleanup';
      $SIG{'SIGUSR1'}='signal_cleanup';


      $SIG{'SIGSEGV'}='signal_cleanup';
      $SIG{'SIGTERM'}='signal_cleanup';
      $SIG{'SIGCONT'}='signal_cleanup';
      $SIG{'SIGIO'}='signal_cleanup';
      $SIG{'SIGPROF'}='signal_cleanup';
      $SIG{'SIGUSR2'}='signal_cleanup';

      $SIG{'SIGSYS'}='signal_cleanup';
      $SIG{'SIGURG'}='signal_cleanup';
      $SIG{'SIGCHLD'}='signal_cleanup';
      $SIG{'SIGXCPU'}='signal_cleanup';
      $SIG{'SIGWINCH'}='signal_cleanup';
      
      $SIG{'INT'}='signal_cleanup';
      $SIG{'TERM'}='signal_cleanup';
      $SIG{'KILL'}='signal_cleanup';
      $SIG{'QUIT'}='signal_cleanup';
      
      our $debug_lock=$ENV{"DEBUG_LOCK"};
      
      
      
      
      foreach my $a (@ARGV){$CL.=" $a";}
      if ( $debug_lock ){print STDERR "\n\n\n********** START PG: $PROGRAM *************\n";}
      if ( $debug_lock ){print STDERR "\n\n\n**********(tcg) LOCKDIR: $LOCKDIR $$ *************\n";}
      if ( $debug_lock ){print STDERR "\n --- $$ -- $CL\n";}
      
	     
      
      
    }
sub flush_error
  {
    my $msg=shift;
    return add_error ($EXIT_FAILURE,$$, $$,getppid(), $msg, $CL);
  }
sub add_error 
  {
    my $code=shift;
    my $rpid=shift;
    my $pid=shift;
    my $ppid=shift;
    my $type=shift;
    my $com=shift;
    
    $ERROR_DONE=1;
    lock4tc ($rpid, "LERROR","LSET","$pid -- ERROR: $type\n");
    lock4tc ($$, "LERROR","LSET", "$pid -- COM: $com\n");
    lock4tc ($$, "LERROR","LSET", "$pid -- STACK: $ppid -> $pid\n");
   
    return $code;
  }
sub add_warning 
  {
    my $rpid=shift;
    my $pid =shift;
    my $command=shift;
    my $msg="$$ -- WARNING: $command\n";
    print STDERR "$msg";
    lock4tc ($$, "LWARNING", "LSET", $msg);
  }

sub signal_cleanup
  {
    print dtderr "\n**** $$ (tcg) was killed\n";
    &cleanup;
    exit ($EXIT_FAILURE);
  }
sub clean_dir
  {
    my $dir=@_[0];
    if ( !-d $dir){return ;}
    elsif (!($dir=~/tmp/)){return ;}#safety check 1
    elsif (($dir=~/\*/)){return ;}#safety check 2
    else
      {
	`rm -rf $dir`;
      }
    return;
  }
sub cleanup
  {
    #print stderr "\n----tc: $$ Kills $PIDCHILD\n";
    #kill (SIGTERM,$PIDCHILD);
    my $p=getppid();
    $CLEAN_EXIT_STARTED=1;
    
    
    
    if (&lock4tc($$,"LERROR", "LCHECK", ""))
      {
	my $ppid=getppid();
	if (!$ERROR_DONE) 
	  {
	    &lock4tc($$,"LERROR", "LSET", "$$ -- STACK: $p -> $$\n");
	    &lock4tc($$,"LERROR", "LSET", "$$ -- COM: $CL\n");
	  }
      }
    my $warning=&lock4tc($$, "LWARNING", "LREAD", "");
    my $error=&lock4tc($$,  "LERROR", "LREAD", "");
    #release error and warning lock if root
    
    if (isrootpid() && ($warning || $error) )
      {
	
	print STDERR "**************** Summary *************\n$error\n$warning\n";

	&lock4tc($$,"LERROR","RELEASE","");
	&lock4tc($$,"LWARNING","RELEASE","");
      } 
    
    
    foreach my $f (@TMPFILE_LIST)
      {
	if (-e $f){unlink ($f);} 
      }
    foreach my $d (@TMPDIR_LIST)
      {
	clean_dir ($d);
      }
    #No More Lock Release
    #&lock4tc($$,"LLOCK","LRELEASE",""); #release lock 

    if ( $debug_lock ){print STDERR "\n\n\n********** END PG: $PROGRAM ($$) *************\n";}
    if ( $debug_lock ){print STDERR "\n\n\n**********(tcg) LOCKDIR: $LOCKDIR $$ *************\n";}
  }
END 
  {
    
    &cleanup();
  }
   

sub safe_system 
{
  my $com=shift;
  my $ntry=shift;
  my $ctry=shift;
  my $pid;
  my $status;
  my $ppid=getppid();
  if ($com eq ""){return 1;}
  
  

  if (($pid = fork ()) < 0){return (-1);}
  if ($pid == 0)
    {
      set_lock($$, " -SHELL- $com (tcg)");
      exec ($com);
    }
  else
    {
      lock4tc ($$, "LLOCK", "LSET", "$pid\n");#update parent
      $PIDCHILD=$pid;
    }
  if ($debug_lock){printf STDERR "\n\t .... safe_system (fasta_seq2hmm)  p: $$ c: $pid COM: $com\n";}

  waitpid ($pid,WTERMSIG);

  shift_lock ($pid,$$, "LWARNING","LWARNING", "LSET");

  if ($? == $EXIT_FAILURE || lock4tc($pid, "LERROR", "LCHECK", ""))
    {
      if ($ntry && $ctry <$ntry)
	{
	  add_warning ($$,$$,"$com failed [retry: $ctry]");
	  lock4tc ($pid, "LRELEASE", "LERROR", "");
	  return safe_system ($com, $ntry, ++$ctry);
	}
      elsif ($ntry == -1)
	{
	  if (!shift_lock ($pid, $$, "LERROR", "LWARNING", "LSET"))
	    {
	      add_warning ($$,$$,"$com failed");
	    }
	  else
	    {
	      lock4tc ($pid, "LRELEASE", "LERROR", "");
	    }
	  return $?;}
      else
	{
	  if (!shift_lock ($pid,$$, "LERROR","LERROR", "LSET"))
	    {
	      myexit(add_error ($EXIT_FAILURE,$$,$pid,getppid(), "UNSPECIFIED system", $com));
	    }
	}
    }
  return $?;
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
		add_warning($$,$$,"Could Not Use EMAIL");
		myexit(add_error ($EXIT_FAILURE,$$,$$,getppid(),"EMAIL","$CL"));
	      }
	    }
	  elsif( $p eq "INTERNET")
	    {
	      if ( !&check_internet_connection())
		{
		  myexit(add_error ($EXIT_FAILURE,$$,$$,getppid(),"INTERNET","$CL"));
		}
	    }
	  elsif( $p eq "wget")
	    {
	      if (!&pg_is_installed ("wget") && !&pg_is_installed ("curl"))
		{
		  myexit(add_error ($EXIT_FAILURE,$$,$$,getppid(),"PG_NOT_INSTALLED:wget","$CL"));
		}
	    }
	  elsif( !(&pg_is_installed ($p)))
	    {
	      myexit(add_error ($EXIT_FAILURE,$$,$$,getppid(),"PG_NOT_INSTALLED:$p","$CL"));
	    }
	}
      return 1;
    }
sub pg_is_installed
  {
    my @ml=@_;
    my $r, $p, $m;
    my $supported=0;
    
    my $p=shift (@ml);
    if ($p=~/::/)
      {
	if (safe_system ("perl -M$p -e 1")==$EXIT_SUCCESS){return 1;}
	else {return 0;}
      }
    else
      {
	$r=`which $p 2>/dev/null`;
	if ($r eq ""){return 0;}
	else {return 1;}
      }
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
sub check_pg_is_installed
  {
    my @ml=@_;
    my $r=&pg_is_installed (@ml);
    if (!$r && $p=~/::/)
      {
	print STDERR "\nYou Must Install the perl package $p on your system.\nRUN:\n\tsudo perl -MCPAN -e 'install $pg'\n";
      }
    elsif (!$r)
      {
	myexit(flush_error("\nProgram $p Supported but Not Installed on your system"));
      }
    else
      {
	return 1;
      }
  }




