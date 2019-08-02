#!/usr/bin/env perl



my $PROBTRESH = 0.3;# base pairs below this prob threshold will be ignored
my $WEIGHT = 100.0; # float!!
my $NUCALPH = "ACGTUNRYMKSWHBVD";
use vars qw($NUCALPH $WEIGHT);

my $myname = basename($0);

use strict;
use warnings;

use File::Basename;
use Getopt::Long;
use File::Glob ':glob';
use File::Spec;
use File::Temp qw/ tempfile tempdir /;




sub tcoffeelib_header($;$)
{
    my ($nseq, $fd) = @_;
    if (! defined($fd)) {
        $fd = *STDOUT;
    }
    printf $fd "! TC_LIB_FORMAT_01\n";
    printf $fd "%d\n", $nseq;
}


sub tcoffeelib_header_addseq($$;$)
{
    my ($id, $seq, $fd) = @_;
    if (! defined($fd)) {
        $fd = *STDOUT;
    }
    printf $fd "%s %d %s\n", $id, length($seq), $seq;
}


sub tcoffeelib_comment($;$)
{
    my ($comment, $fd) = @_;
    if (! defined($fd)) {
        $fd = *STDOUT;
    }
    printf $fd "!" . $comment . "\n";
}


sub tcoffeelib_struct($$$;$)
{
    my ($nseq, $len, $bpm, $fd) = @_;

    if (! defined($fd)) {
        $fd = *STDOUT;
    }

    # output basepair indices with fixed weight
    printf $fd "#%d %d\n", $nseq, $nseq;
    # output basepairs (only once) and with unit-offset
    for (my $i=0; $i<$len; $i++) {
        for (my $j=$i+1; $j<$len; $j++) {
            if (! defined($bpm->[$i][$j])) {
                print STDERR "ERROR: \$bpm->[$i][$j] undefined\n";
            }
            if ($bpm->[$i][$j]>0) {
                print $fd $i+1;
                print $fd " ";
                print $fd $j+1;
                print $fd " " . $bpm->[$i][$j] . "\n";
            }
        }
    }
}


sub tcoffeelib_footer(;$)
{
    my ($fd) = @_;
    if (! defined($fd)) {
        $fd = *STDOUT;
    }
    print $fd "! SEQ_1_TO_N\n";
}


    
sub plfold($$$)
{    
    my ($id, $seq, $probtresh) = @_;
    my (@struct);# return
    my ($templ, $fhtmp, $fnametmp, $cmd, $ctr, $window_size);
    our $ntemp++;
    
    $templ = $myname . ".pid-" . $$ .$ntemp .".XXXXXX";
    ($fhtmp, $fnametmp) = tempfile($templ, UNLINK => 1); 
    print $fhtmp ">$id\n$seq\n";

    # --- init basepair array
    #
    for (my $i=0; $i<length($seq); $i++) {
        for (my $j=$i+1; $j<length($seq); $j++) {
            $struct[$i][$j]=0;
        }
    }


    # --- call rnaplfold and drop a readme
    #
    $window_size=(length($seq)<70)?length($seq):70;
    $cmd = "RNAplfold -W $window_size < $fnametmp >/dev/null";
    system($cmd);
    
    if ($? != 0) {
        printf STDERR "ERROR: RNAplfold ($cmd) exited with error status %d\n", $? >> 8;
        return;
    }
    #unlink($fnametmp);
    my $fps = sprintf("%s_dp.ps", $id); # check long name
    
    if (! -s $fps) {
      {

	$fps = sprintf("%s_dp.ps", substr($id,0,12)); # check short name
 	if (! -s $fps)
	  {
	    die("couldn't find expected file $fps\n");
	    return;
	  }
      }
    }

    
    # --- read base pairs from created postscript
    #
    open(FH, $fps);
    while (my $line = <FH>) {
        my ($nti, $ntj, $prob);
        chomp($line);        
        # line: bp bp sqrt-prob ubox
        my @match = ($line =~ m/^([0-9]+) +([0-9]+) +([0-9\.]+) +ubox$/);
        if (scalar(@match)) {
            $nti=$1;
            $ntj=$2;
            $prob=$3*$3;# prob stored as square root

            if ($prob>$probtresh) {
                #printf STDERR "\$struct[$nti][$ntj] sqrtprob=$3 prob=$prob > $probtresh\n";
                $struct[$nti-1][$ntj-1] = $WEIGHT
            }
            # store with zero-offset
        }
    }
    close(FH);

    # remove or gzi postscript
    #
    unlink($fps);
    #
    # or gzip
    #$cmd = "gzip -qf $fps";
    #system($cmd);
    #if ($? != 0) {
    #    printf STDERR "ERROR: gzip ($cmd) exited with error status %d\n", $? >> 8;
    #}

    return \@struct;
}





sub rnaseqfmt($)
{
    my ($seq) = @_;
    # remove gaps
    $seq =~ s/-//g;
    # uppercase RNA
    $seq = uc($seq);
    # T -> U
    $seq =~ s/T/U/g;
    # check for invalid charaters
    $_ = $seq;
    s/[^$NUCALPH]//g;
    return $_;
}




sub usage(;$)
{    
    my ($errmsg) = @_;
    if ($errmsg) {
        print STDERR "ERROR: $errmsg\n";
    }
    print STDERR << "EOF";
$myname:
 Creates a T-Coffee RNA structure library from RNAplfold prediction.
 See FIXME:citation
Usage:
 $myname -in seq_file -out tcoffee_lib
EOF
    exit(1);
}

sub read_fasta_seq 
  {
    my $f=$_[0];
    my %hseq;
    my (@seq, @com, @name);
    my ($a, $s,$nseq);

    open (F, $f);
    while (<F>)
      {
	$s.=$_;
      }
    close (F);

    
    @name=($s=~/>(\S*).*\n[^>]*/g);
    
    @seq =($s=~/>.*.*\n([^>]*)/g);
    @com =($s=~/>(\S*)(.*)\n([^>]*)/g);


    $nseq=$#name+1;
  
    for ($a=0; $a<$nseq; $a++)
      {
	my $n=$name[$a];
	my $s;
	$hseq{$n}{name}=$n;
	$s=$seq[$a];$s=~s/\s//g;
	
	$hseq{$n}{seq}=$s;
	$hseq{$n}{com}=$com[$a];
      }
    return %hseq;
  }







my $fmsq = "";
my $flib = "";
my %OPTS;
my %seq;
my ($id, $nseq, $i);
my @nl;

GetOptions("in=s" => \$fmsq, "out=s" => \$flib);

if (! -s $fmsq) {
    usage("empty or non-existant file \"$fmsq\"")
}
if (length($flib)==0) {
    usage("empty out-filename")
}






%seq=read_fasta_seq($fmsq);


@nl=keys(%seq);

$nseq=$#nl+1;
open FD_LIB, ">$flib" or die "can't open $flib!";
tcoffeelib_header($nseq, *FD_LIB);
foreach $id (keys (%seq))
  {
    my ($seq, $fmtseq);
    
    $seq = $seq{$id}{seq};
    
    $fmtseq = rnaseqfmt($seq);# check here, formatting for folding important later
    if (length($seq)!=length($fmtseq)) {
        print STDERR "ERROR: invalid sequence $id is not an RNA sequence. read seq is: $seq\n";
        exit
      }
   
    tcoffeelib_header_addseq($id, uc($seq), *FD_LIB);
  }
tcoffeelib_comment("generated by $myname on " . localtime(), *FD_LIB);



$i=0;
foreach $id (keys (%seq))
  {
    my ($cleanid, $seq, $bpm);
    $seq=$seq{$id}{seq};
    $cleanid = $id;
    $cleanid =~ s,[/ ],_,g;# needed for rnaplfold
    $seq = rnaseqfmt($seq);
    
    $bpm = plfold($cleanid, rnaseqfmt($seq), $PROBTRESH);       
    
    tcoffeelib_struct($i+1, length($seq), $bpm, *FD_LIB);
    $i++;
}


tcoffeelib_footer(*FD_LIB);
close FD_LIB;
exit (0);


