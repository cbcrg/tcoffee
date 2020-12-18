#!/usr/bin/env perl
use Env;
use strict;
use FileHandle;
use DirHandle;
use Cwd;
use File::Path;
use Sys::Hostname;
use File::Temp qw/ tempfile tempdir /;

my $QUIET=">/dev/null 2>/dev/null";
my $VERBOSE=$ENV{VERBOSE_4_DYNAMIC};
my $FATAL="[FATAL:mmseqs2prf]";
our $EXIT_FAILURE=1;
our $EXIT_SUCCESS=0;
our $LAST_COM="";

my $tmpdir = File::Temp->newdir();

my $doquiet=0;
my ($outdir,$cachedb,$cacheq);

my $ff=new FileHandle;
my (%db, %q, %P, %H, %T);
my ($dbf, $qf, $out);
my $prot_min_cov=50;
my $prot_min_sim=0;
my $prot_max_sim=100;
#my $psitrim_mode="regtrim";
my $psitrim_mode="sorttrim";

my $psitrim_tree="codnd";
my $psitrim=100;
my $psiJ=1;
my $TF;
my $S;

my $updatedb=0;
my $updateq=0;
my $update=0;
my $qff;
my %CIRCULAR;
my $mmseqsR;
my %R;
my $split=1000000;

for ($a=0; $a<=$#ARGV; $a++)
  {
    if    ($ARGV[$a] eq "-protein_db" || $ARGV[$a] eq "-db"){$dbf=file2abs($ARGV[++$a]);}
    elsif ($ARGV[$a] eq "-q" || $ARGV[$a] eq "-i") {$qff =file2abs($ARGV[++$a]);}
    elsif ($ARGV[$a] eq "-update"){$update=1;}
    elsif ($ARGV[$a] eq "-quiet") {$doquiet=1;}
    
    elsif ($ARGV[$a] eq "-odir") {$outdir=file2abs($ARGV[++$a]);}
    elsif ($ARGV[$a] eq "-o")    {$mmseqsR=file2abs($ARGV[++$a]);}
    
    elsif ($ARGV[$a] eq "-template_file" || $ARGV[$a] eq "-tf") {$TF=($ARGV[++$a]);}
    elsif ($ARGV[$a] eq "-cachedb") {$cachedb=file2abs($ARGV[++$a]);}
    elsif ($ARGV[$a] eq "-updatedb") {$updatedb=1;}
    
    
    elsif ($ARGV[$a] eq "-cacheq")  {$cacheq=file2abs($ARGV[++$a]);}
    elsif ($ARGV[$a] eq "-updateq") {$updateq=1;}
    
    elsif ($ARGV[$a] eq "-prot_min_sim") {$prot_min_sim=($ARGV[++$a]);}
    elsif ($ARGV[$a] eq "-prot_max_sim") {$prot_max_sim=($ARGV[++$a]);}
    elsif ($ARGV[$a] eq "-prot_min_cov") {$prot_min_cov=($ARGV[++$a]);}
    
    elsif ($ARGV[$a] eq "-psitrim_mode") {$psitrim_mode=$ARGV[++$a];}
    elsif ($ARGV[$a] eq "-psiJ")         {$psiJ=$ARGV[++$a];}
    elsif ($ARGV[$a] eq "-psitrim_tree") {$psitrim_tree=$ARGV[++$a];}
    elsif ($ARGV[$a] eq "-psitrim")      {$psitrim=$ARGV[++$a];}
    elsif ($ARGV[$a] eq "-s")            {$S="-s ".$ARGV[++$a];}
    elsif ($ARGV[$a] eq "-split")        {$split=$ARGV[++$a];}
    else {die "$ARGV[$a] is an unknown argument $FATAL";}
  }

if (!$dbf){print STDERR "ERROR: mmseqs2prf required a database via -protein_db $FATAL\n";exit ($EXIT_FAILURE);}
if (!$qff){print STDERR "ERROR: mmseqs2prf required a query    -q $FATAL";exit ($EXIT_FAILURE);}
if (!is_installed("mmseqs")){print STDERR "ERROR: mmseqs2prf required mmseqs to be installed $FATAL\n";exit ($EXIT_FAILURE);}


if   ($cachedb eq "TMP"){$cachedb=$tmpdir;}
elsif(!$cachedb){$cachedb=file2path($dbf);}

if   (!$mmseqsR){$mmseqsR="./out.mmseqs";}
if   (!$cacheq || $cacheq eq "TMP"){$cacheq="$tmpdir/query";}
if   (!$outdir  ){$outdir="./R_dir";}
if   (!$doquiet) {$QUIET="";}

mymkdir ($outdir,$cachedb,$cacheq,$tmpdir);
my ($qf,%H)=filelist2h($outdir,string2fasta_list($qff));
make_output_structure($outdir,%H);


if ( ! -e $mmseqsR || $update){split_mmseqs($qf, $cacheq, $updateq, $dbf, $cachedb, $updatedb,$S,$mmseqsR,$split);}

mmseqs2prf ($mmseqsR,$outdir,$prot_min_sim,$prot_max_sim, $prot_min_cov,%H);
prf2trimprf($outdir,$psitrim_mode, $psitrim_tree, $psitrim, %H);

if ($TF){h2template_file ($TF,%H);}


sub h2template_file
  {
    my ($tf,%h)=@_;
    my $fh =new FileHandle;
    my %lu=my %lu=h2lu(%h);

    open ($fh, ">$tf");
    foreach my $s (keys{%lu})
      {
	my $f=$lu{$s}{0};
	print $fh ">$s _R_ $h{$f}{tprf}{$s}{absolute}\n";
      }
    close ($fh);
    return $tf;
  }
    
sub prf2trimprf
  {
    my ($outdir,$psitrim_mode,$psitrim_tree,$psitrim, %h)=@_;
    my $template =new FileHandle;
    my $template_file;
    my $qf=abs2file($qf);
    
    my %lu=h2lu(%h);
    foreach my $s (keys(%lu))
      {
	my $f=$lu{$s}{0};
	system ("t_coffee -other_pg seq_reformat -in  $h{$f}{prf}{$s}{absolute} -treemode=$psitrim_tree -keep 1 -action +$psitrim_mode $psitrim -output fasta_aln -out $h{$f}{tprf}{$s}{absolute}");
	foreach my $i (keys%{$lu{$s}})
	  {
	   
	    my $f2=$lu{$s}{$i};
	    if ( !-e $h{$f2}{tprf}{$s}{absolute})
	      {
		system ("cp $h{$f}{tprf}{$s}{absolute} $h{$f2}{tprf}{$s}{absolute}");
	      }
	  }
      }
    return %P;
  }




sub h2lu
  {
    my %h=@_;
    my %lu;
    my %count;
    foreach my $f (keys (%h))
      {
	foreach my $s (keys(%{$h{$f}{name}}))
	  {
	    
	    $lu{$s}{$count{$s}++}=$f;
	  }
      }
    return %lu;
  }
sub mmseqs2prf
  {
    #"query[0],target[1],qaln[2],taln[3],qstart[4],qend[5],pident[6],qcov[7],qlen[8]\" $QUIET");
    my ($out,$outdir,$min_id, $max_id,$min_cov, %h)=@_;
    my $ff  =new FileHandle;
    my $prf =new FileHandle;
    my $nn;
    my $tot;
    my $psn;
    my %lu=h2lu(%h);
    my %luf;
    
    open ($ff,$out);
    while (<$ff>)
      {
	my $l=$_;
	my @ll=split (/\s/, $l);
	my $sn=$ll[0];
	my $f =$lu{$sn}{0};
	my $cf=$h {$f}{prf}{$sn}{absolute};
	
	if ($sn ne $psn)
	  {
	    close $prf;
	    if ($luf{$cf}){open ($prf, ">>$cf");}
	    else
	      {		
		open ($prf, ">$cf");
		print $prf ">$sn\n$h{$f}{seq}{$sn}\n";
	      }
	  }
	$nn=++$luf{$cf};
	$psn=$sn;
	
	my $id=$ll[6]*100;
	my $cov=$ll[7]*100;
	my $len=$ll[8];
	print "$id $max_id $min_id\n";
	if ($id<=$max_id && $id>=$min_id && $cov>$min_cov)
	  {
	    print "Keep";
	    print $prf ">$sn\_$nn\n";
	    for (my $a=1; $a<$ll[4]; $a++){print $prf "-"}
	    
	    my @ql=split (//,$ll[2]);
	    my @tl=split (//,$ll[3]);
	    my $qlen=length($ll[2]);
	    for (my $a=0; $a<$qlen; $a++)
	      {
		if ($ql[$a] ne "-"){print $prf "$tl[$a]";}
	      }
	    for (my $a=$ll[5]; $a<$len; $a++){print $prf "-"}
	    print $prf "\n";
	  }
      }
    close($prf);
    close($ff);

    # checkout the un-used ones
    foreach my $sn (keys(%lu))
      {
	my $f=$lu{$sn}{0};

	if (!-e $h{$f}{prf}{$sn}{absolute})
	  {
	    open ($prf,">$h{$f}{prf}{$sn}{absolute}");
	    print $prf ">$sn\n$h{$f}{seq}{$sn}\n";
	    close (prf);
	  }
      }
    #duplicate prf files that are shared by different input datasets
    foreach my $sn (keys (%lu))
      {
	my $f0=$lu{$sn}{0};
	
	foreach my $i (keys(%{$lu{$sn}}))
	  {	
	    my $f=$lu{$sn}{$i};
	    if (! -e $h{$f}{prf}{$sn}{absolute}){system ("cp $h{$f0}{prf}{$sn}{absolute} $h{$f}{prf}{$sn}{absolute}");}
	  }
      }
  }


sub file2db
  {
    my ($in,$dir, $update)=(@_);
    my %f;
    my $out;
    
    if ( !-e $in)
      {
	print "$in does not exists $FATAL \n";
      }
    if ($dir)
      {
	if (!-d $dir){mkdir ($dir) or die "Could not create $dir $FATAL"; }
	$out=$dir."/".abs2file($in);}
    else {$out=$in;}
    
    
    $f{name}=$in;
    $f{db}="$out\.MMSEQSDB";
    $f{index}="$out\.MMSEQSINDEX";
    
    
    #Trigger automated update when source db younger than mmseqs file
    if (-e $f{db} && ((-M $f{db})>(-M $f{name}))){$update=1;}
    
    if (!-e $f{db} || $update)
      {
	system ("mmseqs createdb  $f{name}  $f{db} $QUIET");
      }
    
    if (!-d $f{index} || $update)
      {
	system ("mmseqs createindex $f{name} $f{index} $QUIET");
      }

    return %f;
    
  }

sub file2path
    {
      my ($f)=@_;
      $f=file2abs($f);
      $f=~/(.*\/)[^\/]*$/;
      my $cdir=$1;
      return $cdir;
    }
sub file2abs
     {
       my ($f, $mode)=@_;
       my $cdir=getcwd();
       if ($f=~/^\//){return $f;}
       return "$cdir/$f";
   }
sub abs2file
    {
      my $in=shift @_;
      my $out;
      
      if ( $in=~/\//)
	  {
	    $in=~/.*\/([^\/]*)$/;
	    $out=$1;
	  }
	else
	  {
	    $out=$in;
	  }
      
    return $out;
    }

sub read_fasta_seq
  {
    my $f=@_[0];
    my %hseq;
    my (@seq, @com, @name);
    my ($a, $s,$nseq);
    my $fh=new FileHandle;
    
    open ($fh, $f);
    while (<$fh>)
      {
	$s.=$_;
      }
    close ($fh);


    @name=($s=~/>(\S*).*\n[^>]*/g);

    @seq =($s=~/>.*.*\n([^>]*)/g);
    @com =($s=~/>\S*(.*)\n([^>]*)/g);


    $nseq=$#name+1;

    for ($a=0; $a<$nseq; $a++)
      {
	my $s;
	my $n=$name[$a];
	$hseq{$n}{name}=$n;
	$hseq{$n}{cname}=clean_file_name($n);
	
	$seq[$a]=~s/[^A-Za-z]//g;
	$hseq{$n}{order}=$a;
	$hseq{$n}{seq}=$seq[$a];
	$hseq{$n}{com}=$com[$a];

      }
    return %hseq;
  }
sub mymkdir
    {
      my @l=@_;
      foreach my $a (@l)
	{
	  if ( $a && !-d $a)
	    {
	      system ("mkdir -p $a");
	      if ( !-d $a)
		{
		  die "Could not Create $a $FATAL\n";
		}
	    }
	}

      return 1;
    }
sub clean_file_name
  {
    my $name=@_[0];

    $name=~s/[^A-Za-z1-9.-]/_/g;
    return $name;
  }

sub string2fasta_list
    {
      my $string=@_[0];
      if (!-f $string && !-d $string && !($string=~/\*/)){return ();}
      if ($CIRCULAR{$string}){print STDERR "ERROR: CIRCULAR REFERENCE $string   $FATAL\n";exit ($EXIT_FAILURE);}
      $CIRCULAR{$string}=1;
      

      my @l1=string2list($string);
      my @l2;
      
      foreach my $f (@l1)
	{

	  if (isfasta($f)){push (@l2, $f);}
	  else 
	    {
	      foreach my $string2 (file2list($f))
		{
		  my @l3=string2fasta_list($string2);
		  foreach my $string3 (@l3)
		    {
		      push (@l2, $string3);
		    }
		}
	    }
	}
      return shrinklist(@l2);
    }
sub shrinklist
      {
	my @l=@_;
	my @l2;

	foreach my $e (@l)
	  {
	    if ($e)
	      {
		print "PUSH [$e]\n";
		push (@l2,$e);
	      }
	  }
	return @l2;
      }
sub string2list
    {
      my $string=@_[0];
      my @list;

      if    (-d $string      ){@list= dir2list($string);}
      elsif (   $string=~/\*/){@list= glob    ($string);}
      elsif (-f $string      ){@list=         ($string);}
      return @list;
    }
sub dir2list
  {
    my $dir=shift;
    my @list;
    my $cdir=getcwd;
    my $dh  =new DirHandle;
    
    opendir ($dh, $dir);
    my @dlist=readdir($dh);
    closedir($dh);
    
    foreach my $f (@dlist)
      {
	if ($f eq "." || $f eq ".."){;}
	else {push (@list, string2list("$dir/$f"));}
      }
    return @list;
  }
sub file2list
{
  my $file=shift;
  my @list;
  my $fh=new FileHandle;
  
  open ($fh, "$file");
  while (<$fh>)
    {
      my $l=$_;
      chomp ($l);
      if ($l){push(@list, $l);}
    }
  close ($fh);
  return @list;
}
sub file2string
  {
    my $f=@_[0];
    my ($string, $l);
    my $fh= new FileHandle;
    open ($fh,"$f");
    while (<$fh>)
      {

	$l=$_;
	$string.=$l;
      }
    close ($fh);
    $string=~s/\r\n/\n/g;
    return $string;
  }
sub isfasta
  {
    my $file=shift;
    my $fh=new FileHandle;
    
    open ($fh, "$file");
    while (<$fh>)
      {
	my $l=$_;
	close ($fh);
	if ($l=~/^>/){return 1;}
	return 0;
      }
  }

  



sub filelist2h
  {
    my ($outdir,@list)=@_;
    my $infile="$outdir/fullseq.fa";
    my %h;
    my %lu;
    my $fh=new FileHandle;
    
    open ($fh, ">$infile");
    foreach my $f (@list)
      {
	if (!$QUIET){print "! Process $f\n";}
	my %s=read_fasta_seq($f);
	$f=abs2file ($f);

	$h{$f}{template_dir }="$outdir/$f.template_dir";
	$h{$f}{template_file}{relative}="$f.R.template_file";
	$h{$f}{template_file}{absolute}="$h{$f}{template_dir}/$h{$f}{template_file}{relative}";
		
	foreach my $s (keys(%s))
	  {
	    my $cname=$s{$s}{cname};
	    my $name =$s{$s}{ name};
	    
	    my $seq  =$s{$s}{seq};
	    my $com  =$s{$s}{com};
	    
	    #clean 
	    $seq=~s/\-//g;
	    $seq=~s/\.//g;
	    
	    if ($lu{$s}{seq} && $lu{$s}{seq} ne $seq)
	      {
		#print "$lu{$s}{seq} ---> \n$lu{$s}{source}/$s \n";
		#print "$seq ---> \n$f/$s \n";
		
		#die;
		print "ERROR: two different sequences where provided for $name: [$seq/$f] and [$lu{$s}{seq}/$lu{$s}{source}]$FATAL]\n";
		close (F);
		die;
	      }
	    else 
	      {
		$lu{$s}{seq   }=$seq;
		$lu{$s}{source}=$f;
		
		print $fh ">$name\n$seq\n";
	      }
	    $h{$f}{seq  }{$name}=$seq;
	    $h{$f}{name }{$name}=$name;
	    $h{$f}{cname}{$name}=$cname;
	    $h{$f}{prf}{$name}{relative}="$cname.R.prf";
	    $h{$f}{prf}{$name}{absolute}="$tmpdir/$h{$f}{prf}{$name}{relative}";
	    
	    $h{$f}{tprf}{$name}{relative}="$cname.R.prf";
	    $h{$f}{tprf}{$name}{absolute}="$h{$f}{template_dir}/$h{$f}{tprf}{$name}{relative}";
	    
	    $h{$f}{templates}.=">$name _R_ $h{$f}{tprf}{$name}{relative}\n";
	  }
      }
    close ($fh);
    return $infile, %h;
  }
sub make_output_structure
  {
    my ($outdir,%h)=@_;
    my $fh=new FileHandle;
    mymkdir ($outdir) || die "Could not create $outdir\n";
    
    foreach my $f (keys (%h))
	{
	  mymkdir ($h{$f}{template_dir}) or die "1 Could not create $h{$f}{template_dir}\n";
	  open  ($fh, ">$h{$f}{template_file}{absolute}") or die "Could not open $h{$f}{template_file}{absolute}";
	  print $fh "$h{$f}{templates}";
	  close $fh;
	}
    return 1;
  }
sub is_installed
  {
    my $p=@_[0];
    my $r=0;
    my $cwhich="$tmpdir/which";
    if (-e $cwhich){unlink ("$cwhich");}
    
    system ("which $p >$cwhich 2>/dev/null");
    my $w=file2string ($cwhich);
    if (($w=~/mmseqs/)){$r=1;}
    return $r;
  }

  
sub split_mmseqs
    {
      my ($qf,$cacheq,$updateq, $db,$cachedb, $updatedb, $s,$out,$split)=@_;
      
      my @dbl=splitfasta($split,(string2fasta_list($db)));
      

      if ( -e $out){unlink ($out)}


     
      foreach my $d (@dbl)
	{
	  my $uid=getuid();
	  my $lcacheo="$tmpdir/$uid/search/";
	  my $lcachedb=(($d =~/$tmpdir/))?"$tmpdir/$uid/db/":$cachedb;
	  my $lcacheq=$cacheq;
	  mymkdir ($lcacheo, $lcachedb);
	  
	  print "! Process Database $d\n";
	  
	  my %db=file2db($d ,$lcachedb,$updatedb);
	  my %q =file2db($qf,$lcacheq ,$updateq );
	  my $ld=abs2file ($d);

	  if (! -d $lcacheo){die "NO CACHE";}

	  $q{search }="$lcacheo/$ld\.MMSEQSSEARCH";
	  $q{convert}="$lcacheo/$ld\.MMSEQSCONVERT";
	 
	  system ("mmseqs search $q{db} $db{db} $q{search} $db{index} $s -a $QUIET");
	  system ("mmseqs convertalis $q{db} $db{db} $q{search} $q{convert} --format-output \"query,target,qaln,taln,qstart,qend,pident,qcov,qlen\" $QUIET");
	  system ("cat $q{convert} >> $out");
	  
	}
      return $out;
    }
    
sub splitfasta 
      {
	my ($split,@list)=@_;
	my @fl;

		
	if (!$split){return @list;}
	
	foreach my $e (@list)
	  {
	    
	    my $n=`grep -c ">" $e`;
	    
	    if ($n>$split)
	      {
		my $uid=getuid();
		my $odir="$tmpdir/$uid/";
		mymkdir ($odir);
		system ("t_coffee -other_pg seq_reformat -action +odir $odir +split $e $split");
		push (@fl,string2list ("$odir/*.split"));
	      }
	    else
	      {
		push (@fl, $e);
	      }
	  }
	return @fl;
      }
sub getuid
	{
	  my $n;
	  my $l=3;
	  my $string=randomstring ($l);
	  while ($R{$string})
	    {
	      $n++;
	      
	      if ($n==10){$l++;}
	      $string=randomstring($l);
	    }
	  $R{$string}=1;
	  return $string;
	}
		
sub randomstring
	  {
	    my $l=shift;
	    my @s;
	    my @alp=split (//, 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_');
	    
	    my $lalp=@alp;
	  
	    for ( my $i=0; $i<$l; $i++)
	      {
		$s[$i]=$alp[rand($lalp)];
	      }
	    return join ('',@s);
	  }
	  
		  
