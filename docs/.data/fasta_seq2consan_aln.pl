#!/usr/bin/env perl





my $FMODEL =""; 
my $TMPDIR = "/tmp";




my $NUCALPH = "ACGTUNRYMKSWHBVD";
my $PRIMNUCALPH = "ACGTUN";
use vars qw($NUCALPH $PRIMNUCALPH $TMPDIR);


my $errmsg;
use vars qw($errmsg);



use Getopt::Long;
use Cwd;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use File::Copy;
use File::Path;



sub usage(;$)
{
    my ($errmsg) = @_;
    my $myname = basename($0);

    if ($errmsg) {
        print STDERR "ERROR: $errmsg\n";
    }

    print STDERR << "EOF";
    
$myname: align two sequences by means of consan\'s sfold
Usage:
 $myname -i file -o file -d path
Options:
 -i|--in : pairwise input sequence file
 -o|--out: output alignment
 -d|--directory containing data

EOF
}

sub read_stk_aln 
  {
    my $f=$_[0];
    my ($seq, $id);
    
    my %hseq;

    open (STK, "$f");
    while (<STK>)
      {
	if ( /^#/ || /^\/\// || /^\s*$/){;}
	else
	  {
	    ($id,$seq)=/(\S+)\s+(\S+)/;
	    $hseq{$id}{'seq'}.=$seq;
	  }
      }
    close (STK);
    return %hseq;
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

    
    @name=($s=~/>(.*).*\n[^>]*/g);
    
    @seq =($s=~/>.*.*\n([^>]*)/g);
    @com =($s=~/>.*(.*)\n([^>]*)/g);

    
    $nseq=$#name+1;
    
    for ($a=0; $a<$nseq; $a++)
      {
	my $n=$name[$a];
	$hseq{$n}{name}=$n;
	$hseq{$n}{seq}=$seq[$a];
	$hseq{$n}{com}=$com[$a];
      }
    return %hseq;
  }



sub sfold_parseoutput($$)
{
    my ($frawout, $foutfa) = @_;
    my %haln;
    my ($fstk, $cmd, $id);
    open FOUTFA, ">$foutfa";
    
    $fstk = $frawout . ".stk";
    
    # first line of raw out contains info
    # remaining stuff is stockholm formatted
    $cmd = "sed -e '1d' $frawout";
    system("$cmd > $fstk");
    if ($? != 0) {
        $errmsg = "command failed with exit status $?.";
        $errmsg .=  "Command was \"$cmd\"";
        return -1;
    }

    # this gives an error message. just ignore it...
    %haln=read_stk_aln ( $fstk);
    foreach $i (keys (%haln))
      {
	my $s;
	$s=$haln{$i}{'seq'};
	$s =~ s/\./-/g;
	print FOUTFA ">$i\n$s\n";
      }
    close FOUTFA;
    return 0;
}




sub sfold_wrapper($$$$)
{
    
    my ($fs1, $fs2, $fmodel, $foutfa) = @_;
    

    my ($cmd, $frawout, $ferrlog, $freadme, $ftimelog, $fstk);

    # add  basename($fmsqin) (unknown here!)
    $frawout = "sfold.log";
    $ferrlog = "sfold.err";
    $ftimelog = "sfold.time";
    $freadme =  "sfold.README";
    $fstk = "sfold.stk";
    
    # prepare execution...
    #
    # ./tmp is essential for dswpalign
    # otherwise you'll get a segfault
    mkdir "./tmp";
    
    $cmd = "sfold -m $fmodel $fs1 $fs2";
    open(FREADME,">$freadme");
    print FREADME "$cmd\n"; 
    close(FREADME);

    # and go
    #
    system("/usr/bin/time -p -o $ftimelog $cmd >$frawout 2>$ferrlog");
    if ($? != 0) {
        $errmsg = "command failed with exit status $?";
        $errmsg .= "command was \"$cmd\". See " . getcwd . "\n";
        return -1;
    }

    return sfold_parseoutput($frawout, $foutfa);
}







my ($help, $fmsqin, $fmsaout);
GetOptions("help"  => \$help,
           "in=s" => \$fmsqin,
           "out=s" => \$fmsaout,
	   "data=s" => \$ref_dir);



if ($help) {
    usage();
    exit(0);
}
if (! defined($fmsqin)) {
    usage('missing input filename');
    exit(1);
}
if (! defined($fmsaout)) {
    usage('missing output filename');
    exit(1);

}
if (scalar(@ARGV)) {
    usage('Unknown remaining args');
    exit(1);
}

$FMODEL = "$ref_dir/mix80.mod";
if (! -e "$FMODEL") {
    die("couldn't find sfold grammar model file. Expected $FMODEL\n");
}


my %hseq=read_fasta_seq ($fmsqin);
my $id;

foreach $id (keys(%hseq))
  {
    push(@seq_array, $hseq{$id});
  }

if ( scalar(@seq_array) != 2 ) {
    die("Need *exactly* two sequences as input (pairwise alignment!).")
}



my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
my $datei = sprintf("%4d-%02d-%02d", $year+1900, $mon+1, $mday);
my $templ = basename($0) . "." . $datei . ".pid-" . $$ . ".XXXXXX";
my $wd = tempdir ( $templ, DIR => $TMPDIR);

copy($fmsqin, "$wd/" . basename($fmsqin) . ".org"); # for reproduction
copy($FMODEL, "$wd");
my $fmodel = basename($FMODEL);
my $orgwd = getcwd;
chdir $wd;



my @sepseqfiles;
foreach $id (keys(%hseq)) {
    my ($seq, $orgseq, $fname, $sout);
    $seq=$hseq{$id}{'seq'};
    
    $fname = basename($fmsqin) . "_$id.fa";
    # replace funnies in file/id name (e.g. "/" " " etc)
    $fname =~ s,[/ ],_,g;
    open (PF, ">$fname");
    print (PF ">$id\n$seq\n");
    close (PF);

    push(@sepseqfiles, $fname);
}

my ($f1, $f2, $fout);
$f1 = $sepseqfiles[0];
$f2 = $sepseqfiles[1];
$fout = $wd . basename($fmsqin) . ".out.fa";
if (sfold_wrapper($f1, $f2, $fmodel, "$fout") != 0) {
    printf STDERR "ERROR: See logs in $wd\n";
    exit(1);
} else {
    chdir $orgwd;
    copy($fout, $fmsaout);
    rmtree($wd);
   exit(0);
}

