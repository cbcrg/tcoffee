#!/usr/bin/perl -w
#FILE:    pprun.pl
#AUTHOR:  Karsten Suhre
#DATE:    Thu Feb 9 13:56:40 CET 2006
#PURPOSE: run a list of jobs in parallel threads
#         (see perldoc perlthrtut for details on threading under perl)
#
#BUGS:    STDERR of all jobs is written to screen, not to the log file (although it should)
#MODIF:   

# "Beware of bugs in the above code; I have only proved it correct, not
# tried it."
#               -- Donald Knuth

use Config;
use threads;
use threads::shared;
use Thread::Queue;

# check if PERL is able to deal with threads
$Config{usethreads} or die "Your Perl implementation does not support threads.\n";

# global variables for threads
my $q       = new Thread::Queue;     # job queue
my $fin     = new Thread::Queue;     # queue with finished jobs
my $account = new Thread::Queue;     # queue for accountable time
my $pending : shared;                # number of pending jobs
my $njobs   : shared;                # total number of jobs
my $outdir  : shared = "outdir.pprun.$$";  # subdirectory for all pprun specific I/O
my @prior     : shared;              # prior script
my @posterior : shared;              # posterior script

########################################################################
############## SUB parse_pp ############################################
########################################################################
# read a section from the pp file
sub parse_pp ( $ ) {
  my ($section)=@_;
  my @values;
  open(IN, "$infile") or die "ERROR opening file $infile\n";
  my $nlines = 0;
  my $lread = 0;
  while (<IN>) {
    next if (/^ *$/); 
    if ( /^ *\[$section\] *$/ ) {
      $lread = 1;
      next;
    }
    if ( /^ *\[.*\] *$/ ) {
      $lread = 0;
      next;
    }
    if ( $lread ) {
      chomp; s///g;
      $values[$nlines] = $_;
      $nlines++;
    }
  }
  close IN;
  return @values;
}


########################################################################
############## SUB thread ##############################################
########################################################################
# run a thread that will process several jobs from the queue $q
sub thread ( ) {

  # count the number of tasks that this thread has executed
  my $count = 0;
  my $tid = threads->self->tid;
  my $ctid = "thread-" . $tid; # string to identify thread output
  my $clog; # the log file name, one per job
  my $tstart;
  my $tend;
  my $tdiff;

  print "$ctid: thread $tid is starting\n";
  threads->yield(); # yield to start several threads quickly

  while ( 1 ) {
    { # aquire lock on queue (to make sure $clog is unique)
      lock ($q);
      $pending = $q->pending;
      $jobid = $njobs - $pending + 1;
      $foo = $q->dequeue_nb();
    }
    if (not defined $foo) { # if there are no more jobs in the queue, return
      print "$ctid: no more jobs in queue\n";
      print "$ctid: thread is exiting\n";
      $fin->enqueue ( $tid ); # enqueue tid to signal thread's end
      return $count 
    }
    $clog = "$outdir/job.$jobid.log";
    $cjob = "$outdir/job.$jobid.sh";
    print "$ctid: starting job $jobid\n";
    open (CJOB, ">$cjob") or die "ERROR: $ctid could not create file $cjob\n";
    print CJOB <<EOF1;
# this file has been create by pprun.pl
exec 2>&1 1> $clog

# define job specific variables
export PPRUNDIR="$outdir";
export TID="$tid";
export JOBID="$jobid";

# execute prior script
@prior

echo '================================================================================'
# excute job command line
echo 'executing command: $foo'
echo '================================================================================'
echo '======================== MAIN JOB =============================================='
echo '================================================================================'
$foo
echo '================================================================================'
echo '================================================================================'
echo '================================================================================'

# execute posterior script
@posterior

EOF1
    close CJOB;

    # run the job
    $tstart = `date +"%s%N"`; $tstart = 1*$tstart;
    system "sh $cjob";
    $tend = `date +"%s%N"`; $tend = 1*$tend;
    $tdiff = ($tend - $tstart)*1E-9;

    print  "$ctid: finished job $jobid\n";
    printf  "$ctid: spent %0.2f seconds on job $jobid\n", $tdiff;
    $account->enqueue ( $tdiff ); # enqueue accountable time
    $count ++;
  }
  print "$ctid: getting next job\n";
}
########################################################################
############## MAIN ####################################################
########################################################################

# parse the input line
if ($#ARGV < 0) {
  print STDOUT <<EOF2;
usage: $0 job.pp

# pprun definition file job.pp may be as follows:

# supported sections are :
#  [threads]     number of threads to be run in parallel
#  [prior]       shell script to be run prior to every individual job
#  [posterior]   shell script to be run posterior to every individual job
#  [jobs]        list of jobs to be run in parallel; must be a shell command line
#
#  available environment variables are
#    - PPRUNDIR  : pprun specific directory for job output, can be used for other purposes
#    - TID       : id of the thread, staring with 1
#    - JOBID     : id of the job, starting with 1
#
# Here comes an example file:

---cut---
[threads]
4

[prior]
echo "jobid is \$JOBID, thread id is \$TID"
echo "PPRUNDIR is \$PPRUNDIR"

[jobs]
stresscpu -n 1000
stresscpu -n 1000
stresscpu -n 1000
stresscpu -n 1000

[posterior]
grep GFLOP \$PPRUNDIR/job.\$JOBID.log

---cut---


EOF2
  exit;
}

# get the input file
$infile = $ARGV[0];
die "ERROR, no such file: $infile\n" if (not -e $infile);

# read [threads] section
($nthreads)=parse_pp ( "threads" );
if (not defined $nthreads) {$nthreads = 1};
print "master-0: $nthreads threads to be used\n";

# read [jobs] section
@jobs=parse_pp ( "jobs" );
$njobs = $#jobs+1;
die "ERROR: no jobs found in file $infile\n" if ($njobs <1);
print "master-0: $njobs jobs to be run\n";

# read [prior] section
@prior=parse_pp ( "prior" );
if (not defined $prior[0]) {$prior[0] = "# no prior script defined"};
for (my $l=0; $l<=$#prior; $l++) {
  $prior[$l] .= "\n";
}

# read [posterior] section
@posterior=parse_pp ( "posterior" );
if (not defined $posterior[0]) {$posterior[0] = "# no posterior script defined"};
for (my $l=0; $l<=$#posterior; $l++) {
  $posterior[$l] .= "\n";
}

#############################################################

# create I/O directory
print "master-0: creating I/O directory $outdir\n";
die "ERROR: $outdir already exists\n" if ( -e $outdir );
system "mkdir $outdir";
die "ERROR: $outdir could not be created\n" if ( ! -e $outdir );

# enqueue the jobs
$q->enqueue (@jobs);
print "master-0: $njobs jobs enqueued\n";

# start global timer
$tstart = `date +"%s%N"`; $tstart = 1*$tstart;

# start the threats
print "master-0: starting threads\n";
for (my $i = 1; $i<=$nthreads; $i++) {
  $t[$i] = new threads \&thread;
}

# wait until no more jobs are pending (they may still be under execution)
$nfin = 0;
while ( $nfin < $nthreads ) {
  print "master-0: ", $nthreads - $nfin, " thread(s) are presently running\n";
  $foo = $fin->dequeue();
  $r[$foo] = $t[$foo]->join;
  $nfin++;
  print "master-0: thread $foo terminated\n";
}

# stop timer
$tend = `date +"%s%N"`; $tend = 1*$tend;
$tdiff = ($tend - $tstart)*1E-9;

#############################################################

# print some stats
for ($i = 1; $i<=$nthreads; $i++) {
  print "master-0: thread $i executed $r[$i] commands\n";
}

# get the timing
for ($j=1; $j<=$njobs; $j++) {
  $taccount += $account->dequeue_nb();
}
printf "master-0: accountable time %12.2f seconds\n", $taccount;
printf "master-0: wall clock time  %12.2f seconds\n", $tdiff;
printf "master-0: accountable/wall clock: %0.2f\n", $taccount / $tdiff;
printf "master-0: output files can be found in $outdir\n";

