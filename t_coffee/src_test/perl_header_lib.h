char *PerlScriptName[]={"rec_sum.pl","count.pl","p\
rocess_list.pl","make_license.pl","CCsed.script","\
msa2bootstrap.pl","t_coffee_dpa","t_coffee_dpa2","\
tc_generic_method.pl","generic_method.tc_method","\
clustalw_method.tc_method","extract_from_pdb","ins\
tall.pl","clean_cache.pl","nature_protocol.pl","mo\
cca","dalilite.pl","wublast.pl","blastpgp.pl","ncb\
iblast_lwp.pl","wublast_lwp.pl","RNAplfold2tclib.p\
l","fasta_seq2RNAplfold_templatefile.pl","fasta_se\
q2hmmtop_fasta.pl","fasta_seq2consan_aln.pl","clus\
talw_aln2fasta_aln.pl","msf_aln2fasta_aln.pl","bla\
st_aln2fasta_aln.pl","blast_xml2fasta_aln.pl","fas\
ta_aln2fasta_aln_unique_name.pl","newick2name_list\
.pl","excel2fasta.pl","any_file2unix_file.pl","End\
List"};char *PerlScriptFile[]={"use File::Copy;\nu\
se Env qw(HOST);\nuse Env qw(HOME);\nuse Env qw(US\
ER);\n$x_field=0;\n$y_field=1;\n$interval=0;\n$fil\
e=\"stdin\";\n$print_avg=1;\n$print_sd=0;\n$print_\
sum=0;\n$print_n=0;\nforeach $value ( @ARGV)\n    \
{\n	if ($value ne $ARGV[$np]) \n	    {\n	    ;\n	 \
   }\n	elsif($value eq \"-print_all\")\n	    {\n	 \
   $print_sd=$print_avg=$print_n=$print_sum=1;\n	 \
   $np++;\n	    }\n	elsif($value eq \"-print_sum\"\
)\n	    {\n	    $print_sum=1;\n	    $print_avg=0;\\
n	    $np++;\n	    }\n	elsif($value eq \"-print_n\\
")\n	    {\n	    $print_n=1;\n	    $print_avg=0;\n\
	    $np++;\n	    }\n	elsif($value eq \"-print_avg\
\")\n	    {\n	    $print_avg=1;\n	    $print_avg=0\
;\n	    $np++;\n	    }\n	elsif($value eq \"-sd\")\\
n	    {\n	    $print_sd=1;\n	    $print_avg=0;\n	 \
   $np++;\n	    }\n	elsif($value eq \"-h\")\n	    \
{\n	    $header=1;\n	    $np++;\n	    }\n	elsif ($\
value eq \"-i\")\n	    {\n	    $interval= $ARGV[++\
$np];\n	    $np++;\n    	    }\n	elsif ($value eq \
\"-r\")\n	    {\n	    $min= $ARGV[++$np];\n	    $m\
ax= $ARGV[++$np];\n	    $np++;\n    	    }\n	\n	el\
sif ($value eq \"-x\")\n	    {\n	    $x_field= $AR\
GV[++$np]-1;\n	    $np++;\n    	    }\n	elsif ($va\
lue eq \"-y\")\n	    {\n	      \n	    while ($ARGV\
[$np+1] && !($ARGV[$np+1]=~/\\-/))\n	      {\n		$y\
_field[$nyf++]=$ARGV[++$np]-1;\n		$y_field_set=1;\\
n	      }\n\n	    $np++;\n    	    }\n	elsif ($val\
ue eq \"-file\")\n	    {\n	    $file= $ARGV[++$np]\
;\n	    $file_set=1;\n	    $np++;\n    	    }     \
  \n	elsif ( $value eq \"h\" ||  $value eq \"-h\" \
|| $value eq \"-H\" || $value eq \"-help\" || $val\
ue eq \"help\")\n	  {\n	    print STDOUT \"data_an\
alyse: Analyse and discretization of data\\n\";\n	\
    print STDOUT \"       -file:    <file containi\
ng the data to analyze>,.<def=STDIN>\\n\";\n	    p\
rint STDOUT \"       -x: <field containing the X>,\
...............<Def=0>\\n\";\n	    print STDOUT \"\
       -y: <field containing the Y>,..............\
.<Def=1>\\n\";\n	    print STDOUT \"       -i:<Int\
erval size on the X>,...............<Def=0>\\n\";\\
n	    print STDOUT \"       -i:<0:only one interva\
l>\\n\";\n	    print STDOUT \"       -r:<Range of \
the X>\\n\";\n	    print STDOUT \"       -sd: prin\
t standard deviation on the Y\";\n	    print STDOU\
T \"       -h  : print column header \\n\";\n	    \
exit (0);\n	  }\n	elsif ($value=~/-/)\n	  {\n	    \
print \"$value is not a valid FLAG[FATAL]\\n\";\n	\
    exit (0);\n	   } \n	elsif ($list eq \"\") \n	 \
   {\n	    $file=$ARGV[$np];\n	    $np++;\n	    }\\
n	\n	\n      }\n\n\n\n\n\nif ($file eq \"stdin\")\\
n	{\n	$remove_file=1;\n	$file=\"tmp$$\";\n	open (F\
, \">$file\");\n	while (<STDIN>)\n		{\n		print F $\
_;\n		}\n	close (F);\n	 \n	;}\n\n\nopen(F,$file);\\
n\nif ($interval)\n  {\n    $interval_size=($max-$\
min)/$interval;\n  }\nwhile (<F>)\n  {\n    $line=\
$_;\n    if (!/\\S/){next;}\n    @list=($line=~/(\\
\S+)/g);\n    \n    if ($interval==0){$bin=0;}\n  \
  else{$bin=int (($list[$x_field]-$min)/($interval\
_size));}\n\n    \n    if ($bin && $bin==$interval\
){$bin--;}\n    for ( $a=0; $a<$nyf; $a++)\n      \
{\n	$sum{$a}{$bin}+=$list[$y_field[$a]];\n	$sum2{$\
a}{$bin}+=$list[$y_field[$a]]*$list[$y_field[$a]];\
\n	$n{$a}{$bin}++;\n      }\n  }\n\nif (!$interval\
){$interval=1;}\nfor ( $a=0; $a<$interval; $a++)\n\
  {\n    printf ( \"%3d %3d \", $interval_size*$a,\
 $interval_size*($a+1));\n    for ( $b=0; $b<$nyf;\
 $b++)	\n      {\n	$i=$interval*$a;\n	if ( $n{$b}{\
$a}==0)\n	  {\n	    $avg=0;\n	    $sd=0;\n	  }\n	e\
lse\n	  {\n	    $avg=$sum{$b}{$a}/$n{$b}{$a};\n	  \
  $sd=sqrt($sum2{$b}{$a}*$n{$b}{$a}-$sum{$b}{$a}*$\
sum{$b}{$a})/($n{$b}{$a}*$n{$b}{$a});\n	  }\n	if (\
$print_n) {printf \"%10.4f \", $n{$b}{$a};}\n	if (\
$print_sum){printf \"%10.4f \", $sum{$b}{$a};}\n	i\
f ($print_avg){printf \"%10.4f \", $avg}\n	if ($pr\
int_sd) {printf \"%10.4f \", $sd;}\n      }\n    p\
rintf (\"\\n\");\n  }\n\n\nif ( $remove_file){unli\
nk $file;}\n","use File::Copy;\nuse Env qw(HOST);\\
nuse Env qw(HOME);\nuse Env qw(USER);\n\nforeach $\
v (@ARGV){$cl.=$v;}\n\n\nif ( $cl=~/-k(\\d+)/){$k=\
$1;}\nelse {$k=1;}\nif ( $cl=~/-w(\\d+)/){$w=$1;}\\
nelse {$w=-1;}\nif ( $cl=~/-p(\\d+)/){$p=$1;}\nels\
e {$p=-1;}\n\nwhile (<STDIN>)\n  {\n    @l=($_=~/(\
\\S+)/g);\n    $v=$l[$k-1];\n    if ( !$h{$v}){@ll\
=($v, @ll);}\n    \n    if ( $w==-1)\n      {$h{$v\
}++;}\n    else\n      {$h{$v}+=$l[$w-1];}\n\n    \
if ($p!=-1){$print{$v}=$l[$p-1];}\n\n  }\nforeach \
$v (@ll)\n  {\n    print \"$v $print{$v} $h{$v}\\n\
\";\n  }\n","\nuse Env qw(HOST);\nuse Env qw(HOME)\
;\nuse Env qw(USER);\n$random_tag=int (rand 10000)\
+1;\n$unique_prefix=\"$$.$HOST.$random_tag\";\n$qu\
eue=\"distillery.and.mid\";\n$monitor=0;\n$stderr_\
file=\"/dev/null\";\n$stdio_file=\"/dev/null\";\n$\
log_file=\"/dev/null\";\n$pause_time=0;\n$max_sub_\
jobs=60;\n$min_sub_jobs=30;\n$output_all=0;\n$var=\
'\\$';\n\nforeach $value ( @ARGV)\n    {\n	if ($va\
lue ne $ARGV[$np]) \n	    {\n	    ;\n	    }\n	elsi\
f ($value eq \"-max_sub_jobs\")\n	    {\n	    $max\
_sub_jobs= $ARGV[++$np];\n	    $np++;\n    	    }	\
\n	elsif ($value eq \"-min_sub_jobs\" )\n	    {\n	\
    $min_sub_jobs= $ARGV[++$np];\n	    $np++;\n   \
 	    }\n	elsif ($value eq \"-para\")\n	    {\n	  \
  $para=1;\n	    $monitor=1;\n	    $np++;\n    	  \
  }\n	elsif ($value eq \"-monitor\") \n	    {\n	  \
  $monitor=1;\n	    $np++;\n	    }\n	elsif ($value\
 eq \"-no_monitor\") \n	    {\n	    $monitor=0;\n	\
    $np++;\n	    }\n	elsif ($value eq \"-queue\")\\
n	    {\n	    $queue=$ARGV[++$np];\n	    $np++;\n	\
    }	\n	elsif ($value eq \"-stderr_file\")\n	    \
{\n	    $stderr_file=$ARGV[++$np];\n	    $np++;\n	\
    }\n	elsif ($value eq \"-stdio_file\")\n	    {\\
n	    $stdio_file=$ARGV[++$np];\n	    $np++;\n	   \
 }\n	elsif ($value eq \"-output_all\")\n	    {\n	 \
   $output_all=1;\n	    $np++;\n	    }\n	elsif ($v\
alue eq \"-pause\") \n	    {\n	    $pause_time=$AR\
GV[++$np];\n	    $np++;\n	    }\n	elsif ($value eq\
 \"-log\")\n	      {\n	       $log=1;\n	       \n	\
       if ($ARGV[$np+1]=~/\\-\\S+/) \n	          {\
\n		  $log_file=\"stderr\";\n	          }\n	      \
 else \n	          {\n		  $log_file=$ARGV[++$np]; \
\n		  ++$np;\n		 \n	          }\n	      }\n	elsif \
( $value eq \"-com\")\n	    {\n		\n		if (!$ARGV[$n\
p+1]=~/^\\'/) { $com=$ARGV[++$np];}\n		else {$com=\
$ARGV[++$np];}\n\n	     $np++;\n	    }\n	elsif ( $\
value eq \"-check\")\n	  {\n	    \n	    if (!$ARGV\
[$np+1]=~/^\\'/) { $check=$ARGV[++$np];}\n	    els\
e {$check=$ARGV[++$np];}\n	    $np++;\n	  }\n	elsi\
f ($com eq \"\") \n	    {\n	    $com_set=1;\n	    \
$com=$ARGV[$np];\n	    \n	    $np++;\n	    }\n	els\
if ($list eq \"\") \n	    {\n	    $list_set=1;\n	 \
   $list=$ARGV[$np];\n	    $np++;\n	    }\n	elsif \
( $var_set eq \"\")\n	    {\n	    $var_set=1;\n	  \
  $var=$ARGV[$np];\n	    $np++;\n	    }\n	}\n\n\n\\
n\nif ( $com eq \"\"){print \"You Need to Provide \
a Command [FATAL]\\n\";\n	      die;\n	     }\n\n\\
n\nif ($list_set==0) \n    {\n    $x= int (rand 10\
0000)+1;\n    $tmp_file_name=\"tmp_file_$x\";\n   \
 open ( TMP, \">$tmp_file_name\");\n    while (<ST\
DIN>)\n      {\n	print TMP $_;\n      }\n    close\
 (TMP);\n    open (F, $tmp_file_name);\n    }\nels\
e \n    {\n    open (F, $list);\n    }\n\nif ($par\
a==0) \n    {\n\n     @tc_list= <F>;\n     close (\
F); \n     \n     foreach $val(@tc_list) \n	    {\\
n	      \n	      \n	      \n	      $loc_com=$com;\\
n	      if ($check){$loc_check=$check;}\n	      \n\
	      @i_val=($val=~/([^\\s]+)/g);\n	      \n	   \
   if ( $#i_val==0)\n		{\n		  if ($check){$loc_che\
ck=~s/$var/$i_val[0]/g;}\n		  $loc_com=~s/$var/$i_\
val[0]/g;\n		}\n	      else\n		{\n		  for ($n=1; $\
n<=$#i_val+1;$n++ )\n		    {\n		      \n		      $s\
ub=\"$var$n\";\n		      \n		      $loc_com=~s/$sub\
/$i_val[$n-1]/g;\n		      if ($check){$loc_check=~\
s/$var/$i_val[0]/g;}\n		    }\n		}\n	      if ( $c\
heck && -e $loc_check)\n		{\n		  print STDERR \"sk\
ipping $loc_com...\\n\";\n		  }\n	      else\n		{\\
n		  system \"$loc_com\";\n		}\n	    }\n    exit;\\
n    }\n\nelsif ($para==1) \n    {\n    print STDE\
RR \"do parallel execution of: \\\"$com $list\\\"\\
\n\";\n    \n    if ($log==1) \n	{\n	if ($log_file\
 eq \"stdout\" || $log_file eq \"stderr\" ) \n		{\\
n		$log_file=\"\";\n	        }\n\n        else \n	\
	{\n		system \"echo LOG FILE> $log_file\";\n		\n	 \
       }\n	}\n    else	\n	{\n	open ( OUT, \">/dev/\
null\");\n	}\n	\n    \n    $id=0;\n    $n_sub=0;\n\
    while ($val=<F>) \n	    {	    	    \n	    $job\
_log[$id]=\"$HOME/tmp/$unique_prefix.$id.log_file\\
";\n	    \n	    $job=$unique_prefix.\"_$id\";\n	  \
  open (JOB, \">$job\");\n	    \n	    $loc_com=$co\
m;\n	    chop $val;\n\n	    $loc_com=~s/\\$/$val/g\
;\n	 \n	    print JOB \"#!/bin/csh\\n\";\n	    pri\
nt JOB \"#\\$ -cwd\\n\";\n	    print JOB \"#\\$ -N\
 $unique_prefix\\n\";\n	    if ($queue && !($queue\
 eq \" \")) {print JOB \"#\\$ -l $queue\\n\";}\n	 \
   print JOB \"#\\n\";	    \n            print JOB\
 \"$loc_com\\n\";\n	    print JOB \"echo FINISHED \
 >> $job_log[$id]\\n\";\n	    print JOB \"pwd\\n\"\
;\n	    \n	    close (JOB);\n	    if ( $output_all\
==1)\n		{\n		system \"qsub $job >  $unique_prefix\\
";		\n	        }\n	    else\n		{system \"qsub $job\
 -e $stderr_file -o $stdio_file >$unique_prefix\";\
	        \n	        } \n\n\n\n	    print STDERR \"\
$id: $output_all\\n\";\n	    $n_sub++;\n	    if ( \
$max_sub_jobs && $n_sub==$max_sub_jobs) \n		{\n		$\
n_sub=monitor_process($min_sub_jobs,@job_log); 		 \
\n		\n	        }	\n	   \n            unlink $uniqu\
e_prefix;\n	    sleep $pause_time;\n	    $id++;\n	\
    }\n\n    close (OUT);\n    close (F);\n\n    p\
rint STDERR \"Your $id Jobs Have Been Submited (NA\
ME=$unique_prefix)\\n\";\n    monitor_process (0, \
@job_log);\n    foreach $file(@job_log) {if (-e $f\
ile) {unlink($file);}}\n    \n    }\n\nsub monitor\
_process ( @job_list)\n    {\n    my (@job_list)=@\
_;\n    my $min_sub_jobs=shift (@job_list);\n    m\
y $n_sub_jobs;\n    my $finished;\n    my $n=0;\n\\
n    $n_sub_jobs=-1;\n    $finished=0;\n    print \
STDERR \"\\nMonitor Batch: [$min_sub_jobs]\";\n   \
    \n    while (!$finished && (($n_sub_jobs>$min_\
sub_jobs)|| $n_sub_jobs==-1) ) \n	{\n	$finished=1;\
\n	$n_sub_jobs=0;\n	$n=0;\n	foreach $file (@job_li\
st)\n	        {\n	\n		if (-e $file){;}\n		else \n	\
	    {\n		    $finished=0; $n_sub_jobs++;\n	      \
      }\n	        }\n	system \"sleep 1\";\n       \
 }\n    \n    return $n_sub_jobs;\n    }\n    \n  \
  \nif ($tmp_file_name){unlink($tmp_file_name);}\n\
","\n\nforeach ($np=0; $np<=$#ARGV; $np++)\n    {\\
n    $value=$ARGV[$np];\n\n    if ($value eq \"-fi\
le\")\n      {\n      $file= $ARGV[++$np];\n      \
}\n    elsif ($value eq \"-type\")\n      {\n     \
   $type= $ARGV[++$np];\n      }\n    elsif ($valu\
e eq \"-institute\")\n      {\n        $institute=\
 $ARGV[++$np];\n      }\n    elsif ($value eq \"-a\
uthor\")\n      {\n        $author= $ARGV[++$np];\\
n      }\n    elsif ($value eq \"-date\")\n      {\
\n        $date= $ARGV[++$np];\n      }\n     elsi\
f ($value eq \"-program\")\n      {\n        $prog\
ram= $ARGV[++$np];\n      }\n    elsif ($value eq \
\"-email\")\n      {\n        $email= $ARGV[++$np]\
;\n      }\n    else\n      {\n	print \"$value is \
an unkown argument[FATAL]\\n\";\n	exit (1);\n     \
 }\n  }\n\n\n\nopen F, $file || die;\nprint $INSTI\
TUTE;\nif ( $type eq \"c\"){print \"/*************\
********************COPYRIGHT NOTICE**************\
********************/\\n\";}\nif ( $type eq \"perl\
\"){print \"#################################COPYR\
IGHT NOTICE#################################/\\n\"\
;}\nif ( $type eq \"txt\"){print \"---------------\
-------------------COPYRIGHT NOTICE---------------\
------------------/\\n\";}\n\n\nwhile (<F>)\n  {\n\
  s/\\$INSTITUTE/$institute/g;\n  s/\\$AUTHOR/$aut\
hor/g;\n  s/\\$DATE/$date/g;\n  s/\\$PROGRAM/$prog\
ram/g;  \n  s/\\$EMAIL/$email/g;  \n  if ( $type e\
q \"txt\"){print $_;}\n  elsif ($type eq \"c\"){ch\
op $_; print \"\\/*$_*\\/\\n\";}\n  elsif ($type e\
q \"perl\"){print \"\\#$_\";}\n}\nclose (F);\nif (\
 $type eq \"c\"){print \"/************************\
*********COPYRIGHT NOTICE*************************\
*********/\\n\";}\nif ( $type eq \"perl\"){print \\
"#################################COPYRIGHT NOTICE\
#################################/\\n\";}\nif ( $t\
ype eq \"txt\"){print \"--------------------------\
--------COPYRIGHT NOTICE--------------------------\
-------/\\n\";}\n\n","\nwhile (<>)	\n	{\n	s/\\=cc/\
123456789/g;\n	s/\\bcc/\\$\\(CC\\)/g;\n	s/12345678\
9/\\=cc/g;\n	print $_;\n	}\n\n","$version=\"1.00\"\
;\n$rseed= int(rand(100000))+1;\n\n\nif ( $#ARGV==\
-1)\n  {\n    print \"msa2bootstrap -i <input_file\
> -input <seq|msa|matrix|tree> -n <N-Boostrap> -o \
<outtree> -tmode <nj|upgma|parsimony|ml> -dmode <k\
imura> -alignpg <t_coffee | muscle | clustalw> -rt\
ree <file> -stype <prot|cdna|dna> -recompute -syst\
em <cygwin|unix>\";\n    print \"\\n\\t-i: input f\
ile, can be sequneces, msa, matrix, trees, type is\
 specified via -input\";\n    print \"\\n\\t-input\
: Type of input data\";\n    print \"\\n\\t\\tmsa:\
 msa in fasta format\";\n    print \"\\n\\t\\tseq:\
 compute an msa with -alignpg\";\n    print \"\\n\\
\t\\tmatrix: phylipp distance matrix fed directly \
to method -tmode [caveat: tmode=nj or upgma]\";\n \
   print \"\\n\\t\\ttree: list of newick trees dir\
ectly fed to consence in order to generate a boots\
traped tree\";\n    \n    print \"\\n\\t-n: number\
 of bootstrap replicates\";\n    print \"\\n\\t-o:\
 name of the output tree. Files are not overwritte\
n. Use -recompute to overwrite existing file\";\n \
   print \"\\n\\t-tmode: tree mode: nj|upgma|parsi\
mony|ml\";\n    print \"\\n\\t-dmode: distance mod\
e\";\n    print \"\\n\\t-alignpg: program for alig\
ning sequences (t_coffee=default)\";\n    print \"\
\\n\\t-rtree: replicate tree file (default: no fil\
e)\";\n    print \"\\n\\t-rmsa: replicate msa file\
 (default: no file)\";\n    print \"\\n\\t-rmat: r\
eplicate matrix file (default: no file)\";\n    pr\
int \"\\n\\t-stype: sequence type: protein, dna or\
 cdna\";\n    print \"\\n\\t-recompute: force file\
s to be overwritten\";\n    print \"\\n\\t-system:\
 cygwin|unix\";\n      \n\n    \n    &my_exit (EXI\
T_FAILURE);\n  }\nforeach $arg (@ARGV){$command.=\\
"$arg \";}\n\nprint \"CLINE: $command\\n\";\n$thre\
shold=100;\n$trim_msa=0;\n$stype=\"prot\";\nprint \
\"msa2bootstrap \";\n\n$system=\"cygwin\";\nif(($c\
ommand=~/\\-system (\\S+)/))\n  {\n    $system=$1;\
\n    if ( $system eq \"cygwin\")\n      {\n	$exec\
_extension=\".exe\";\n      }\n    elsif ( $system\
 eq \"unix\")\n      {\n	$exec_extension=\"\";\n	p\
rint \"system=Unix\";die;\n      }\n    else\n    \
  {\n	print \"msa2boostrap: -system=$system is an \
unknown mode [FATAL]\\n\"; die;\n      }\n    \n  \
  print \"-system $system \";\n  }\nif(($command=~\
/\\-stype (\\S+)/))\n  {\n    $stype=$1;\n  }\npri\
nt \"-stype=$stype \";\n\n\n\nif(($command=~/\\-i \
(\\S+)/))\n  {\n    $msa=$1;\n    print \"-i $msa \
\";\n  }\n\nif(($command=~/\\-rtree (\\S+)/))\n  {\
\n    $rtree=$1;\n    print \"-rtree=$rtree \";\n \
 }\n\nif(($command=~/\\-rmsa (\\S+)/))\n  {\n    $\
rmsa=$1;\n  }\nif(($command=~/\\-rmat (\\S+)/))\n \
 {\n    $rmat=$1;\n  }\n$input=\"seq\";\nif(($comm\
and=~/\\-input (\\S+)/))\n  {\n    $input=$1;\n  }\
\nprint \"-input=$input \";\n\n$dmode=\"kimura\";\\
nif(($command=~/\\-dmode (\\S+)/))\n  {\n    $dmod\
e=$1;\n  }\nprint \"-dmode=$dmode \";\n$alignpg=\"\
muscle\";\nif(($command=~/\\-alignpg (\\S+)/))\n  \
{\n    $alignpg=$1;\n  }\nprint \"-alignpg=$dmode \
\";\n\n$tmode=\"nj\";\nif(($command=~/\\-tmode (\\\
S+)/))\n  {\n    $tmode=$1;\n  }\nprint \"-tmode=$\
tmode \";\n$recompute=0;\nif(($command=~/\\-recomp\
ute/))\n  {\n    $recompute=1;\n    print \"-recom\
pute \";\n  }\n\n$out=$msa;\n$out=~s/\\..*//;\n$ou\
t.=\".bph\";\nif(($command=~/\\-o (\\S+)/))\n  {\n\
    $out=$1;\n    \n  }\nprint \"-out=$out \";\nif\
 (-e $out && !$recompute)\n  {\n    print \"\\nNo \
Computation Required $out already exists\\n\";\n  \
  &my_exit (EXIT_SUCCESS);\n    \n  }\n\n$n=100;\n\
if(($command=~/\\-n (\\d+)/))\n  {\n    $n=$1;\n  \
}\nprint \"-n=$n \";\n$seed=3;\nif(($command=~/\\-\
s (\\d+)/))\n  {\n    $seed=$1;\n  }\nprint \"-s=$\
seed \";\n\nif(($command=~/\\-run_name (\\d+)/))\n\
  {\n    $suffix=$1;\n  }\nelse\n  {\n    $msa=~/(\
[^.]+)/;\n    $suffix=$1;\n  }\nprint \"-run_name=\
$suffix\\n\";\n\n\nif ( $input eq \"seq\")\n  {\n \
   $seq=$msa;\n    $msa=\"$suffix.prot_msa\";\n   \
 \n    if ($stype eq \"cdna\")\n      {\n	$cdna_se\
q=$seq;\n	$clean_cdna_seq=&vtmpnam();\n	$seq=&vtmp\
nam();\n	`t_coffee -other_pg seq_reformat -in $cdn\
a_seq -action +clean_cdna >$clean_cdna_seq`;\n	`t_\
coffee -other_pg seq_reformat -in $clean_cdna_seq \
-action +translate >$seq`;\n	\n      }\n\n    if (\
!-e $msa || $recompute)\n      {\n	print \"\\n####\
#   Compute an MSA With $alignpg\\n\";\n	\n	if ( $\
alignpg eq \"t_coffee\")\n	  {`$alignpg $seq -outf\
ile=$msa >/dev/null 2>/dev/null`;}\n	elsif ( $alig\
npg eq \"muscle\")\n	  {\n	    `$alignpg -in $seq \
> $msa 2>/dev/null`;\n	  }\n	elsif ( $alignpg eq \\
"clustalw\")\n	  {\n	    `$alignpg -infile=$seq -o\
utfile=$msa -quicktree >/dev/null 2>/dev/null`;\n	\
  }\n	elsif ( $align eq \"mafft\")\n	  {\n	    `$a\
lignpg $seq > $msa >/dev/null 2>/dev/null`;\n	  }\\
n	else\n	  {\n	    `$alignpg -in=$seq -outfile=$ms\
a`;\n	  }\n      }\n    if (!-e $msa)\n      {\n	p\
rint \"\\nError: $alignpg Could Not produce the MS\
A $msa [FATAL]\\n\";\n      }\n\n    if ($stype eq\
 \"cdna\")\n      {\n	$msa2=\"$suffix.cdna_msa\";\\
n	`t_coffee -other_pg seq_reformat -in $clean_cdna\
_seq -in2 $msa -action +thread_dna_on_prot_aln -ou\
tput fasta_aln  >$msa2`;\n	$msa=$msa2;\n      }\n \
   \n    $input=\"msa\";\n  }\n\n\n\n$seqboot_o=&v\
tmpnam();\n$seqboot_c=&vtmpnam();\n\n$protdist_o=&\
vtmpnam();\n$protdist_c=&vtmpnam();\nif ( $input e\
q \"msa\")\n  {\n    if ($tmode eq \"nj\" || $tmod\
e eq \"upgma\"){$input=\"matrix\";}\n    \n    $lm\
sa= &vtmpnam ();\n    `t_coffee -other_pg seq_refo\
rmat -in $msa -output phylip_aln > $lmsa`;\n    \n\
    if ( -e \"outfile\"){unlink (\"outfile\");}\n \
   # run seqboot\n  \n    if ( $n>1)\n      {\n	pr\
int \"Run SeqBoot .....\";\n	open (F, \">$seqboot_\
c\");\n	print F \"$lmsa\\nR\\n$n\\nY\\n$seed\\n\";\
\n	close (F);\n	`seqboot$exec_extension  < $seqboo\
t_c`;\n	if ( -e \"outfile\"){ print \"[OK]\\n\";}\\
n	else { print \"[FAILED]\\n\";&my_exit (EXIT_FAIL\
URE);}\n	`mv outfile $seqboot_o`;\n      }\n    el\
se\n      {\n	`cp $lmsa $seqboot_o`;\n      }\n\n \
   if ($rmsa){`cp $seqboot_o $rmsa`;}\n    \n    i\
f ($tmode eq \"nj\" || $tmode eq \"upgma\")\n     \
 {\n	if ( $stype eq \"prot\")\n	  {\n	    # run pr\
otdist\n	    print \"Run Protdist [dmode=$dmode]\"\
;\n	    if ($dmode eq \"kimura\")\n	      {\n		$dm\
ode=\"P\\nP\\nP\";\n	      }\n	    else\n	      {\\
n		print \"\\n$dmode is an unknown mode for Protdi\
st [FATAL:msa2bootstrap.pl]\\n\";\n		&my_exit (EXI\
T_FAILURE);\n	      }\n	    open (F, \">$protdist_\
c\");\n	    if ($n>1){print F \"$seqboot_o\\n$dmod\
e\\nM\\nD\\n$n\\nY\\n\";}\n	    else {printf F \"$\
seqboot_o\\n$dmode\\nY\\n\";}\n	    close (F);\n	 \
   `protdist$exec_extension  < $protdist_c`;\n	   \
 if ( -e \"outfile\"){ print \"[OK]\\n\";}\n	    e\
lse { print \"[FAILED]\\n\";&my_exit (EXIT_FAILURE\
);}\n	    `mv outfile $protdist_o`;\n	 \n	  }\n	el\
sif ( $stype eq \"cdna\" || $stype eq \"dna\")\n	 \
 {\n	    print \"Run dnadist [dmode=default\";\n	 \
   open (F, \">$protdist_c\");\n	    if ($n>1){pri\
nt F \"$seqboot_o\\nM\\nD\\n$n\\nY\\n\";}\n	    el\
se {printf F \"$seqboot_o\\nY\\n\";}\n	    close (\
F);\n	    `protdist$exec_extension  < $protdist_c`\
;\n	    if ( -e \"outfile\"){ print \"[OK]\\n\";}\\
n	    else { print \"[FAILED]\\n\";&my_exit (EXIT_\
FAILURE);}\n	    `mv outfile $protdist_o`;\n	  }\n\
      }\n  }\nelsif ( $input eq \"matrix\")\n  {\n\
    $protdist_o=&vtmpnam();\n    print \"MSA: $msa\
\\n\";\n    `cp $msa $protdist_o`;\n    $n=1;\n  }\
\n\n\n\n\n\n$nb_o=&vtmpnam();\n$nb_c=&vtmpnam();\n\
if ($input eq \"matrix\" && $tmode ne \"parsimony\\
" && $tmode ne \"ml\")\n  {\n    print \"Run neigh\
bor [tmode=$tmode]\";\n\n    if ($tmode eq \"nj\")\
\n      {\n	$tmode=\"\\nN\\nN\";\n      }\n    els\
if ( $tmode eq \"upgma\")\n      {\n	$tmode = \"\\\
nN\";\n      }\n    else\n      {\n	print \"\\n ER\
ROR: $tmode is an unknown tree computation mode\\n\
\";\n	&my_exit (EXIT_FAILURE);\n      }\n\n    ope\
n (F, \">$nb_c\");\n    if ($n>1){print F \"$protd\
ist_o$tmode\\nM\\n$n\\n$seed\\nY\\n\";}\n    else \
{print F \"$protdist_o$tmode\\nY\\n\";}\n    close\
 (F);\n\n    `neighbor$exec_extension  < $nb_c`;\n\
    if ( -e \"outtree\"){ print \"[Neighbor OK]\\n\
\";}\n    else { print \"[FAILED]\\n\";&my_exit (E\
XIT_FAILURE);}\n    `mv outtree $nb_o`;\n    unlin\
k (\"outfile\");\n  }\nelsif ($input eq \"msa\" &&\
 $tmode eq \"parsimony\")\n  {\n    if ( -e \"outf\
ile\"){unlink (\"outfile\");}\n    if ( -e \"outtr\
ee\"){unlink (\"outtree\");}\n    \n    if ($stype\
 eq \"prot\")\n      {\n	print \"Run protpars [tmo\
de=$tmode]\";\n	open (F, \">$nb_c\");\n	if ($n>1){\
print F \"$seqboot_o\\nM\\nD\\n$n\\n$seed\\n10\\nY\
\\n\";}\n	else {print F \"$seqboot_o\\nY\\n\";}\n	\
close (F);\n	`protpars$exec_extension  < $nb_c`;\n\
      }\n    elsif ( $stype eq \"dna\" || $stype e\
q \"cdna\")\n      {\n	print \"Run dnapars [tmode=\
$tmode]\";\n	open (F, \">$nb_c\");\n	if ($n>1){pri\
nt F \"$seqboot_o\\nM\\nD\\n$n\\n$seed\\n10\\nY\\n\
\";}\n	else {print F \"$seqboot_o\\nY\\n\";}\n	clo\
se (F);\n	`dnapars$exec_extension  < $nb_c`;\n    \
  }\n    if ( -e \"outtree\"){ print \"[OK]\\n\";}\
\n    else { print \"[FAILED]\\n\";&my_exit (EXIT_\
FAILURE);}\n    `mv outtree $nb_o`;\n   unlink (\"\
outfile\");\n  }\nelsif ($input eq \"msa\" && $tmo\
de eq \"ml\")\n  {\n    if ( -e \"outfile\"){unlin\
k (\"outfile\");}\n    if ( -e \"outtree\"){unlink\
 (\"outtree\");}\n    \n    if ($stype eq \"prot\"\
)\n      {\n	print \"Error: ML impossible with Pro\
tein Sequences [ERROR]\";\n	&my_exit (EXIT_FAILURE\
);\n      }\n    elsif ( $stype eq \"dna\" || $sty\
pe eq \"cdna\")\n      {\n	print \"Run dnaml [tmod\
e=$tmode]\";\n	open (F, \">$nb_c\");\n	if ($n>1){p\
rint F \"$seqboot_o\\nM\\nD\\n$n\\n$seed\\n10\\nY\\
\n\";}\n	else {print F \"$seqboot_o\\nY\\n\";}\n	c\
lose (F);\n	`dnaml$exec_extension  < $nb_c`;\n    \
  }\n    if ( -e \"outtree\"){ print \"[OK]\\n\";}\
\n    else { print \"[FAILED]\\n\";&my_exit (EXIT_\
FAILURE);}\n    `mv outtree $nb_o`;\n   unlink (\"\
outfile\");\n  }\n\n\nelse\n  {\n    `cp $msa $nb_\
o`;\n    $n=2;\n  }\n\nif ($rmsa && -e $seqboot_o)\
{print \"\\nOutput List of $n Replicate MSA: $rmsa\
\\n\";`cp $seqboot_o $rmsa`;}\nif ($rmat && -e $pr\
otdist_o){print \"\\nOutput List of $n Replicate M\
ATRICES: $rmat\\n\";`cp $protdist_o $rmat`;}\nif (\
$rtree && -e $nb_o){print \"\\nOutput List of $n R\
eplicate TREES: $rtree\\n\";`cp $nb_o $rtree`;}\n\\
n\n\n$con_o=&vtmpnam();\n$con_c=&vtmpnam();\nif ($\
n >1)\n  {\n    print \"Run Consense.....\";\n    \
open (F, \">$con_c\");\n    print F \"$nb_o\\nY\\n\
\";\n    close (F);\n    `consense$exec_extension \
 < $con_c`;\n    if ( -s \"outtree\"  > 0) { print\
 \"[OK]\\n\";}\n    else { print \"[FAILED]\\n\";&\
my_exit (EXIT_FAILURE);}\n    `mv outtree $con_o`;\
\n    unlink (\"outfile\");\n  }\nelse\n  {\n    `\
cp $nb_o $con_o`;\n  }\n\n\n`cp $con_o $out`;\nif \
( !-e $out)\n  {\n    print \"Tree Computation fai\
led [FAILED]\\n\";\n    &my_exit (EXIT_FAILURE);\n\
  }\nelsif ($n>1)\n  {\n    print \"\\nOutput Boot\
strapped Tree: $out\\n\";\n    $avg=`t_coffee -oth\
er_pg seq_reformat -in $out -action +avg_bootstrap\
`;\n    $avg=~s/\\n//g;\n    print \"$avg\\n\";\n \
 }\nelse\n  {\n    print \"\\nOutput Tree: $out\\n\
\";\n  }\n\nopen (F, \"$out\");\nwhile (<F>)\n  {\\
n    \n    $tree.=$_;\n  }\nclose (F);\n$tree=~s/\\
\n//g;\nprint \"BPH: $tree\\n\";\n\n\n&my_exit (EX\
IT_SUCCESS);\n\nsub my_exit \n  {\n    my $m=@_[0]\
;\n    &clean_vtmpnam();\n    exit ($m);\n  }\nsub\
 vtmpnam \n  {\n    my $file;\n\n\n    $ntmp++;\n \
   $file=\"tmp4msa2bootstrap.$rseed.$$.$ntmp\";\n \
   \n    push (@tmpfile, $file);\n    return $file\
;\n  }\nsub clean_vtmpnam \n  {\n    my $t;\n    f\
oreach $t (@tmpfile)\n      {\n	if ( -e $t){unlink\
 ($t)};\n      }\n  }\n","use Env;\n$seq_reformat=\
\"t_coffee -other_pg seq_reformat \";\n$VersionTag\
=\"1.00\";\n$step=1;\n$unset=\"\";\n$scoreT1=$scor\
eT2=$nseqT=$dp_limit=$unset;\n@tl=();\nchomp($tc_v\
ersion=`t_coffee -version`);$tc_version=~s/PROGRAM\
: //;\n\n\nprint STDERR \"\\n*********************\
********************************************\";\np\
rint STDERR \"\\n*           HIGH LEVEL PROGRAM: T\
-COFFEE_DPA Version $VersionTag\";\nprint STDERR \\
"\\n*           LOW  LEVEL PROGRAM: $tc_version \"\
;\nprint STDERR \"\\n*****************************\
************************************\";\n\nif (!@A\
RGV)\n  {\n    print \"t_coffee_dpa accepts every \
t_coffee_flag.\\nType t_coffee to obtain a list\\n\
\";\n    print \"Requires $TC_VERSION\\n\";\n    p\
rint \"Requires \";\n    print \"t_coffee_dpa spec\
ific flags:\\n\";\n    print \"\\t-dpa_master_aln.\
...................Master alignment: provided OR c\
omputed\\n\";\n    print \"\\t-dpa_master_aln.....\
...............By default, Computed with t_coffee \
-very_fast\\n\";\n    print \"\\t-dpa_master_aln=<\
file>.............Use file, (must be an aln in Fas\
ta or ClustalW\\n\";\n    print \"\\t-dpa_master_a\
ln=<program>..........Compute aln with pg -in seq \
-out aln`\\n\";\n    print \"\\t-dpa_maxnseq......\
.................Maximum number of sequences in su\
bgroups\\n\";\n    print \"\\t-dpa_min_score1.....\
...............Minimum Id for two sequences to be \
grouped in ref_aln\\n\";\n    print \"\\t-dpa_min_\
score2....................Minimum Id within a subg\
roup\\n\";\n    print \"\\t-dpa_debug.............\
............Keep Tmp File (for debug purpose)\\n\\\
n\";\n    \n    exit (0);\n  }\nforeach $arg (@ARG\
V)\n  {\n    $arg_list.=\" $arg\";\n  }\n$arg_list\
=~s/[=,;]/ /g;\n\n\n($seq0, $arg_list)=&extract_va\
l_from_arg_list(\"^\",$arg_list, \"SPLICE\",\"unse\
t\");\n($seq1, $arg_list)=&extract_val_from_arg_li\
st(\"-seq\",$arg_list, \"SPLICE\",\"unset\");\n($s\
eq2, $arg_list)=&extract_val_from_arg_list(\"-in\"\
,$arg_list, \"KEEP\",\"unset\");\n($seq3, $arg_lis\
t)=&extract_val_from_arg_list(\"-infile\",$arg_lis\
t, \"SPLICE\",\"unset\");\n($prf,  $arg_list)=&ext\
ract_val_from_arg_list(\"-profile\",$arg_list, \"S\
PLICE\",\"unset\");\n\n$gl{'Seq'}=$seq=&vtmpnam();\
#file containing all the sequences\n\n   #1-remove\
 sequences from -in\nif ( $arg_list =~/\\-in\\b/)\\
n  {\n    my $save, $name;\n    while($arg_list=~/\
\\-in\\b[^-]+(\\bS[\\w.]+)/)\n      {\n	$name=$1;$\
name=~s/^.//;\n	if ( !-e $name){$save.=\" S$name \\
";}\n\n	$arg_list=~s/S$name/ /;\n      }\n    $arg\
_list=~s/\\-in\\b/\\-in $save /;\n  }\n   #2-prepa\
re \n\nif (!($arg_list=~/\\-outorder/))\n  {\n    \
\n    $output_cl .=\" -outorder=$seq\";\n  }\n@out\
put_flag=(\"-output\",\"-outfile\", \"-run_name\",\
 \"-outorder\"); \nforeach $v1 (@output_flag)\n  {\
\n    ($v2, $arg_list)=&extract_val_from_arg_list(\
$v1,$arg_list, \"SPLICE\",\"unset\");\n    if ($v2\
 ne \"\")\n      {\n\n	if ($v1 eq \"-run_name\"){$\
run_name=$v2;$output_cl .=\" $v1 $v2 \";}\n	elsif \
( $v1 eq \"-outorder\")\n	  {\n	    if ( $v2 eq \"\
input\"){$v2=$seq;}\n	    $outorder=$v2;$output_cl\
 .=\" $v1 $v2 \";\n	  }\n	else\n	  {\n	    $output\
_cl .=\" $v1 $v2 \";\n	  }\n      }\n }\n\n\n($dpa\
_master_aln, $arg_list)  =&extract_val_from_arg_li\
st(\"-dpa_master_aln\",$arg_list, \"SPLICE\", \"t_\
coffee\");\n$dpa_master_aln=~s/\\s//g;\n($nseqT, $\
arg_list)           =&extract_val_from_arg_list(\"\
-dpa_maxnseq\",$arg_list, \"SPLICE\", 30);\n($scor\
eT1, $arg_list)         =&extract_val_from_arg_lis\
t(\"-dpa_min_score1\",$arg_list, \"SPLICE\", 80);\\
n($scoreT2, $arg_list)         =&extract_val_from_\
arg_list(\"-dpa_min_score2\"    ,$arg_list, \"SPLI\
CE\", 30);\n($dpa_limit, $arg_list)       =&extrac\
t_val_from_arg_list(\"-dpa_limit\"        ,$arg_li\
st, \"SPLICE\", 0);\n($dpa_delta_id, $arg_list)   \
 =&extract_val_from_arg_list(\"-dpa_delta_id\"    \
    ,$arg_list, \"SPLICE\", 1);\n($dpa_debug, $arg\
_list)       =&extract_val_from_arg_list(\"-dpa_de\
bug\"           ,$arg_list, \"SPLICE\", 0);\n\n\n$\
in_seq=$seq0.\" \".$seq1.\" \".$seq2.\" \".$seq3;\\
n$in_prf=(($prf ne $unset)?\"$prf \":\"\");\n&exit\
_dpa (($in_seq eq \"\" && $in_prf eq \"\")?1:0, \"\
ERROR: You did not Provide any sequences. Use the \
-seq flag [FATAL: t_coffee_dpa]\\n\", EXIT_FAILURE\
);\n\n\nprint STDERR \"\\nSTART DPA COMPUTATION\";\
\n\n\n\nif ($in_seq=~/\\S+/)\n  {\n    \n    print\
 STDERR \"\\n Step $step: Gather all the sequences\
 into the tmp file: [$seq]\";$step++;	\n    &my_sy\
stem (\"t_coffee $in_seq -convert -quiet -output f\
asta_seq -outfile=$seq -maxnseq 0\");\n  }\n\nif (\
 !-e $seq){$seq=\"\";}\n\nif ($in_prf=~/\\S+/)\n  \
{\n    $seq_in_type=\"profile\"; \n    $seq.= $in_\
prf; \n  }\nif ($seq eq \"\"){ &exit_dpa (1, \"\\n\
ERROR: No Sequence FOund. Provide Sequences with t\
he -seq flag [FATAL: t_coffee_dpa]\", EXIT_FAILURE\
);}\n\n \n\nif ( $run_name)\n  {\n    $suffix=$run\
_name;\n  }\nelsif ($in_seq=~/\\b(S[\\w.]+\\b)/)\n\
  {\n    my $suffix1, $sufffix2;\n    $suffix1=$su\
ffix2=$1;\n    $suffix2=~s/^S//;\n    if ( -e $suf\
fix1){$suffix=$suffix1;}\n    elsif ( -e $suffix2)\
{$suffix=$suffix2;}\n    else\n      {\n	$suffix=&\
vtmpnam();	\n      }\n    $suffix=~s/\\.\\w+//;\n \
 }\n\nelse\n  {\n    $suffix=&vtmpnam();\n  }\n\n\\
nif (!$run_name){$output_cl.=\" -run_name $suffix \
\";}\n\n\n$gl{'Tree'}=&seq2dpa_tree ($seq, \"$suff\
ix.dpadnd\");\n\nprint STDERR \"\\n Step $step: Pr\
epare guide tree: $gl{'Tree'}\";$step++;\n\nprint \
STDERR \"\\n Step $step: Identify and Align Closel\
y Related Groups\";$step++;\n%gl=&make_one_pass (0\
, $scoreT1,\"Align\",%gl);\n\nprint STDERR \"\\n S\
tep $step: Make Multiple Group Alignment\";$step++\
;\nwhile (!%gl ||$gl{'Ng'}>$nseqT)\n  {\n    %gl=&\
make_one_pass ($nseqT, $scoreT2,\"t_coffee\",%gl);\
\n    if ( $gl{'Newgroups'}==0){$scoreT2--;}    \n\
  }\nprint STDERR \"\\n Step $step: Make The Final\
 Alignment\";$step++;\n\n\n$arg_list .=$output_cl;\
\n\n\n%gl=&tree2group (0,0, %gl);\n$gl{$gl{'0'}{'F\
ile'}}{'Output'}=\"\";\n$a=0;\n&align_groups (\"t_\
coffee\",'0', $arg_list, \" \", %gl);\n\n\n\nif ( \
!$dpa_keep_tmpfile){&clean_tmp_file (@tl);}\n\n\n\\
nsub seq2dpa_tree \n  {\n    my $seq=@_[0];\n    m\
y $newtree=@_[1];\n    my $aln=&vtmpnam ();\n\n   \
 &my_system (\"t_coffee -special_mode quickaln -in\
 $seq -outfile $aln -quiet\");\n    &my_system (\"\
$seq_reformat -in $aln -action +aln2tree +tree2dpa\
tree -output newick >$newtree\");\n    return $new\
tree;\n  }	\nsub seq2dpa_tree_old \n  {\n    my $a\
ln=@_[0];\n    my $newtree=@_[1];\n    \n    \n   \
 &my_system(\"$seq_reformat -in $aln -action +seq2\
dpatree -output newick > $newtree\");\n    return \
$newtree;\n  }\nsub aln2dpa_tree \n  {\n    my $al\
n=@_[0];\n    my $newtree=&vtmpnam();\n    \n    &\
my_system(\"$seq_reformat -in $aln -action +aln2tr\
ee +tree2dpatree -output newick > $newtree\");\n  \
  return $newtree;\n  }\nsub group_file2ngroups\n \
 {\n    my $file=@_[0];\n    my $n;\n    \n    ope\
n ( F, $file);\n    while (<F>)\n      {\n	$n+=/\\\
>/;\n      }\n    close (F);\n    return $n;\n  }\\
n\nsub make_one_pass\n  {\n    my ($N, $ID,$pg, %g\
l)=@_;\n    my $a;\n\n    %gl=&tree2group ($N,$ID,\
%gl);\n    if (!$gl{'Newgroups'}){return %gl;}\n  \
  else\n      {\n	for ( $a=0; $a< $ng; $a++)\n	  {\
\n	    if ($gl{$gl{$a}{'File'}}{'Ng'}>1){&display_\
group($a, %gl);}\n	    &align_groups ($pg, $a, $ar\
g_list, \" -quiet=quiet \", %gl);\n	  }\n	return %\
gl;\n      }\n  }\n\nsub tree2group \n  {\n    my \
($N, $ID, %gl)=@_;\n    my $prefix=&vtmpnam();\n  \
  my $group_file=&vtmpnam();\n    my $file;\n    m\
y $oldtree=&vtmpnam();\n    my $n;\n    my $tree;\\
n\n\n    if ( $gl{'Ng'}==1){return %gl;}\n    $tre\
e=$gl{'Tree'}; \n    \n    #1 extract the groups\n\
    &my_system (\"$seq_reformat -in $tree -action \
+tree2group $N $ID $prefix > $group_file\");\n    \
$n=group_file2ngroups($group_file);\n    \n    \n \
   $gl{'Newgroups'}=1;\n    if ( $n==$gl{'Ng'})\n \
     {\n	$gl{'Newgroups'}=0;\n	return %gl;\n      \
}\n    $gl{'Iteration'}++;\n    $gl{'MaxNseq'}=$N;\
$gl{'MinID'}=$ID;\n    $gl{'GroupFile'}=$group_fil\
e;$gl{'Ng'}=$ng=0;\n    #2 Process the group list \
into the hash\n    open (F, $group_file);\n    whi\
le (<F>)\n      {\n	$gl{'File'}.=$_;\n	if (/\\>/)\\
n	  {\n	    $line=$_;\n	    $line=~s/\\>//;\n	    \
@list=($line=~/(\\S+)/g);\n	    $file=$gl{$ng}{'Fi\
le'}=shift @list;\n	    $gl{$file}{'Output'}=$file\
;\n	    \n	    $gl{$file}{'Ng'}=$#list+1;\n	    if\
 ($gl{$file}{'Ng'}>1){ $gl{$file}{'Tlist'}=$gl{$fi\
le}{'Alist'}=\"(\";}\n	    foreach $l (@list)\n	  \
    {\n	\n		$gl{$file}{'List'}.=\" $l \";\n		\n		i\
f (!$gl{$l}{'Tlist'})\n		  {\n		    $gl{$l}{'Tlist\
'}=\"$l\";\n		    $gl{$l}{'Alist'}=\"$l\";\n		    \
$gl{$l}{'Nseq'}=1;\n		    $gl{$l}{'Ng'}=1;\n		  }\\
n		$gl{$file}{'Tlist'}.=\"$gl{$l}{'Tlist'},\";\n		\
$gl{$file}{'Alist'}.=\"$gl{$l}{'Tlist'}|\";\n		$gl\
{$file}{'Nseq'}+=$gl{$l}{'Nseq'};\n	      }\n	    \
\n\n	    chop($gl{$file}{'Tlist'});chop($gl{$file}\
{'Alist'});\n	    if ($gl{$file}{'Ng'}>1){$gl{$fil\
e}{'Tlist'}.=\")\"; $gl{$file}{'Alist'}.=\");\";}\\
n	    $ng++;\n	  }	\n      }\n    $gl{'Ng'}=$ng;\n\
    close (F);\n    \n    #3 Update the old tree w\
ith the new groups\n    $gl{'Tree'}=&vtmpnam();\n \
   &my_system (\"$seq_reformat -in $tree -action +\
collapse_tree $group_file -output newick > $gl{'Tr\
ee'}\");\n    \n    return %gl;\n  }\n\nsub displa\
y_group \n  {\n    my ($g,%gl)=@_;\n    my $f;\n  \
  \n    if ( $g==-1)\n      {\n	print STDERR \"\\n\
Iteration $gl{'Iteration'} [MaxN=$gl{'MaxNseq'}][M\
inID=$gl{'MinID'}]\";\n      }\n    else\n      {\\
n\n	$f=$gl{$g}{'File'};\n	$action=($gl{$f}{'Ng'}==\
1 || $gl{'Iteration'}==1)?\"KEEP  \":\"ALIGN \";\n\
        print STDERR \"\\n\\t[$action][MaxN=$gl{'M\
axNseq'}][MinID=$gl{'MinID'}][File $f][Nseq=$gl{$f\
}{'Nseq'}][Ngroups=$gl{$f}{'Ng'}][$gl{$f}{'Alist'}\
]\";\n      }\n  }\n      \n\n\nsub align_groups\n\
  {\n    my ($pg, $g, $arg, $extra_arg,%gl)=@_;\n \
   my $f;\n    my $Output,$Outflag;\n    \n    \n \
   $f=$gl{$g}{'File'};\n    $Output=($gl{$f}{'Outp\
ut'});\n    \n    if ( $pg eq \"Align\")\n      {\\
n	if ( !-e $f)\n	  {\n	    $command=\"$seq_reforma\
t -in $gl{'Seq'}  -action +extract_aln $gl{'GroupF\
ile'}\";\n	    if ($gl{$f}{'Ng'}>1)\n	      {\n		&\
my_system ($command);\n		$command=\"t_coffee -spec\
ial_mode quick_aln  S$f -outfile=$Output -quiet\";\
\n	      }\n	  }\n	else \n	  {$command=\"\";}\n   \
   }\n    elsif ( -e $f)\n      {	\n	$Outflag=($Ou\
tput)?\"-outfile=$Output\":\"\";\n	$command=\"$pg \
-infile $f $Outflag -quiet stdout $arg $extra_arg \
-maxnseq 0 -convert -quiet stdout\";\n      }\n   \
 elsif ( $gl{$f}{'Ng'}==1)\n      {\n	$action=($dp\
a_debug)?\"cp\":\"mv\";\n	$command=\"$action $gl{$\
f}{'List'} $Output\";\n      }\n    else\n      {\\
n	$Outflag=($Output)?\"-outfile=$Output\":\"\";\n	\
$command=\"$pg -profile $gl{$f}{'List'} $Outflag $\
arg $extra_arg -maxnseq 0\";\n      }\n    \n    &\
my_system ($command);\n    return $outfile;\n  }\n\
    \nsub my_system \n  {\n    my $command=@_[0];\\
n    my $force=@_[1];\n    my $status;\n\n    if (\
 $dpa_debug) {print STDERR \"\\nCOMMAND: $command\\
";}\n    $status=system ($command);\n\n    if (!$f\
orce)\n       {\n	 &exit_dpa (($status==1), \"Fail\
ed in Command:\\n$command\\n[FATAL: t_coffee_dpa]\\
\n\", EXIT_FAILURE);\n       }\n    \n    return $\
status;\n  }\n\nsub vtmpnam\n  {\n    my $prefix=@\
_[0];\n    my $tmp_file_name;\n\n    $tmp_prefix=(\
$prefix)?$prefix:\"dpa_tmp_file_$$\";\n   \n    $t\
mp_count++;\n    $tmp_file_name=\"$tmp_prefix\".\"\
$tmp_count\";\n    $tl[$#tl+1]=$tmp_file_name;\n  \
  return $tmp_file_name;\n  }\n\nsub clean_tmp_fil\
e\n  {\n\n    my $list;\n    my $file;\n    \n    \
if ($dpa_debug){return;}\n    $list=vtmpnam();\n  \
  `ls -1 | grep $tmp_prefix>$list`;\n    \n    ope\
n (F,$list);\n    while ( <F>)\n      {\n	$file=$_\
;\n	chop $file;\n	if ( -e $file){unlink $file;}\n \
     }\n    close (F);\n    unlink $list;\n  }\n\n\
\nsub exit_dpa\n  {\n  my $condition=@_[0];\n  my \
$error_msg=@_[1];\n  my $exit_value=@_[2];\n  if (\
 $condition)\n    {\n      print \"$error_msg\\n\"\
;\n      exit ($exit_value);\n    }\n  else\n    {\
\n      return;\n    }\n  \n}\nsub extract_val_fro\
m_arg_list\n  {\n    my $arg=@_[0];\n    my $arg_l\
ist=@_[1];\n    my $keep_flag=@_[2];\n    my $defa\
ult_value=@_[3];\n    my $val=\"\";\n    \n    #pr\
otect\n    $arg_list=~s/\\s-/ \\@/g;\n    $arg=~s/\
-/\\@/g;\n    \n    #search\n    if ($arg eq \"^\"\
)\n      {\n	$arg_list=~/^([^@]*)/;\n	$val=$1;\n  \
    }\n    else\n      {$arg_list=~/$arg ([^@]*)/;\
$val=$1;}\n    \n    #remove trailing spaces\n    \
$val=~s/\\s*$//;\n    \n    #remove the parsed seq\
uence if needed\n    if (($val ne \"\") && $keep_f\
lag ne \"KEEP\")\n      {\n	if ( $arg eq \"^\"){$a\
rg_list=~s/$val/ /;}\n	else {$arg_list=~s/($arg [^\
@]*)/ /;}\n      }\n	\n    #unprotect\n    $arg_li\
st=~s/\\@/-/g;\n    $arg=~s/\\@/-/g;\n    \n    if\
 (($val eq \"\") && $default_value ne \"unset\"){$\
val=$default_value;}\n    \n    return $val, $arg_\
list;\n  }\n$program=\"T-COFFEE (Version_8.79)\";\\
n\n","\n$DEBUG=1;\n$dpa_nseq=10;\n$dpa_sim=0;\nif \
(!@ARGV)\n  {\n    `t_coffee`;\n    exit (0);\n  }\
\nforeach $arg (@ARGV)\n  {\n    $arg_list.=\" $ar\
g\";\n  }\n$max_nseq=10;\n($seq0, $arg_list)=&extr\
act_val_from_arg_list(\"^\",$arg_list);\n($seq1, $\
arg_list)=&extract_val_from_arg_list(\"-seq\",$arg\
_list);\n($seq2, $arg_list)=&extract_val_from_arg_\
list(\"-in\",$arg_list, \"KEEP\");\n($seq3, $arg_l\
ist)=&extract_val_from_arg_list(\"-infile\",$arg_l\
ist);\n$in_seq=$seq0.\" \".$seq1.\" \".$seq2.\" \"\
.$seq3;\n\n$seq=vtmpnam();\n`t_coffee $in_seq -con\
vert -output fasta_seq -outfile=$seq`;\n\n\n($dpa_\
nseq, $arg_list)=&extract_val_from_arg_list(\"-dpa\
_nseq\",$arg_list);\n($master_aln, $arg_list)=&ext\
ract_val_from_arg_list(\"-master_aln\",$arg_list);\
\n($sim_matrix, $arg_list)=&extract_val_from_arg_l\
ist(\"-sim_matrix\",$arg_list);\n($core_seq, $arg_\
list)=&extract_val_from_arg_list(\"-core_seq\",$ar\
g_list);\n($dpa_sim, $arg_list)=&extract_val_from_\
arg_list(\"-dpa_sim\",$arg_list);\n($run_name, $ar\
g_list)=&extract_val_from_arg_list(\"-run_name\",$\
arg_list);\n($output, $arg_list)=&extract_val_from\
_arg_list(\"-output\",$arg_list);\n\n\n\nif (!$sim\
_mat && !$master_aln)#Compute the fast alignment\n\
  {\n    $ref_aln=vtmpnam();\n    `t_coffee -seq=$\
seq -very_fast -outfile=$ref_aln -quiet`;\n    \n \
 }\n\nif (!$sim_mat)\n  {\n    $sim_mat=vtmpnam();\
\n    `seq_reformat -in $ref_aln -output sim > $si\
m_mat`;\n  }\n\nif ( !$core_seq)\n  {\n    $core_s\
eq=vtmpnam();\n    `seq_reformat -in $ref_aln -act\
ion +trimTC N$max_nseq -output fasta_seq > $core_s\
eq`;\n  }\n@core_name=`seq_reformat -in $core_seq \
-output name `; \n\n@tot_name=`seq_reformat -in $s\
eq -output name `;\n\nforeach $s (@core_name){$s=~\
s/\\s//g;$hcore{$s}=1;}\nforeach $s (@tot_name){$s\
=~s/\\s//g;}\nprint STDERR \"T-Coffee_dpa:\\n\";\n\
print STDERR \"\\tTOTAL  SEQ: @tot_name\\n\";\npri\
nt STDERR \"\\tCHOSEN SEQ: @core_name\\n\";\n\n\n\\
nopen (F, $sim_mat);\nwhile ( <F>)\n  {\n    @l=($\
_=~/(\\b[\\S]+\\b)/g);\n    if (($l[0] eq \"TOP\" \
|| $l[0] eq \"BOT\"))\n      {\n	$s1=$l[1];$s2=$l[\
2];$v=$l[3];\n	if ($hcore{$s1} && !$hcore{$s2})\n	\
  {\n	    if (!$hseq{$s2}{\"sim\"} || $v>$hseq{$s2\
}{\"sim\"})\n	      {\n		$hseq{$s2}{\"sim\"}=$v;$h\
seq{$s2}{\"seq\"}=$s1;\n	      }\n	  }\n      }\n \
 }\nclose (F);\nforeach $s (@tot_name)\n  {\n\n   \
 if ( !$hseq{$s}{\"seq\"}){;}\n    else\n      {\n\
	$s2=$hseq{$s}{\"seq\"};\n	$v=$hseq{$s}{\"sim\"};\\
n		\n	if ($v>$dpa_sim)\n	  {\n	    $hseq{$s}{'used\
'}=1;\n	    $seq_list{$s2}{$seq_list{$s2}{'nseq'}+\
+}=$s;\n	  }\n      }\n  }\nforeach $s (@core_name\
){$seq_list{$s}{$seq_list{$s}{'nseq'}++}=$s;$hseq{\
$s}{'used'}=1;}\nforeach $s (@tot_name){if (!$hseq\
{$s}{'used'}){$seq_list{'unused'}{$seq_list{'unuse\
d'}{'nseq'}++}=$s;}}\n\n\n$n=0;\nforeach $s (@core\
_name)\n  {\n    $ng++;\n    $n=$seq_list{$s}{'nse\
q'};\n    for (@g_list=(), $a=0; $a<$n; $a++){@g_l\
ist=(@g_list,$seq_list{$s}{$a});}\n\n    $g_seq=vt\
mpnam();\n    $g_aln=vtmpnam();\n    \n    print S\
TDERR \"Group $ng: $#g_list Seq: @g_list: \";\n   \
 \n    \n    `seq_reformat -in $seq -action +lower\
 +keep_name +extract_seq  @g_list -output fasta_se\
q > $g_seq`;\n    \n    \n    if ( $#g_list==0)\n \
     {\n	print STDERR \"[No aln]\\n\";\n	$g_aln=$g\
_seq;\n      }\n    elsif ($#g_list<$max_nseq) \n \
     {\n	print STDERR \"[t_coffee]\\n\";\n	`t_coff\
ee $g_seq -outfile=$g_aln -quiet $arg_list`;\n    \
  }\n    else\n      {\n	print STDERR \"[t_coffee_\
dpa]\\n\";\n	`t_coffee_dpa2 $g_seq -outfile=$g_aln\
 $arg_list -sim_matrix $sim_matrix -dpa_nseq $dpa_\
nseq`;\n      }\n    @profile_list=(@profile_list,\
 $g_aln);\n  }\n\n\nprint \"UNUSED $seq_list{'unus\
ed'}{'nseq'}\";\n\nif ($seq_list{'unused'}{'nseq'}\
)\n    {\n      $prf=vtmpnam();\n      \n      `t_\
coffee -profile @profile_list $arg_list -outfile=$\
prf -quiet`;\n      $n=$seq_list{\"unused\"}{'nseq\
'};\n      $new_seq=vtmpnam();\n      $new_prf=vtm\
pnam();\n      for ($a=0; $a<$n-1; $a++)\n	{\n	  $\
s=$seq_list{\"unused\"}{$a};\n	  print STDERR \"\\\
nADD Sequence $s\";\n	  \n	  `seq_reformat -in $se\
q -action +lower +keep_name +extract_seq $s  -outp\
ut fasta_seq > $new_seq`;\n	  `t_coffee -profile $\
prf $new_seq $arg_list -outfile=$new_prf`;\n	  `cp\
 $new_prf $prf`;\n	}\n      $s=$seq_list{\"unused\\
"}{$a};\n      `seq_reformat -in $seq -action +low\
er +keep_name +extract_seq $s  -output fasta_seq >\
 $new_seq`;\n      @profile_list=($prf, $new_seq);\
\n    }\n    \n      \nif ($run_name){$arg_list.=\\
" -run_name $run_name\";}\nelse \n  {\n    $in_seq\
=~/([\\w-]+)/;\n    $arg_list.=\" -run_name $1\";\\
n  }\nif ( $output){$arg_list.=\" -output $output \
\";}\n\n`t_coffee -profile @profile_list $arg_list\
`;\n\n\n&clean (@tmp_file_list);\n\n\nsub vtmpnam\\
n  {\n    my $tmp_file_name;\n    $tmp_name_counte\
r++;\n    $tmp_file_name=\"tmp_file_$tmp_name_coun\
ter\\_Pid$$\";\n    $tmp_file_list[$ntmp_file++]=$\
tmp_file_name;\n    return $tmp_file_name;\n  }\ns\
ub clean\n  {\n  my @fl=@_;\n  my $file;\n  return\
;\n\n  foreach $file ( @fl)\n    {\n      if ( -e \
$file){unlink($file);}\n    }\n}\nsub extract_val_\
from_arg_list\n  {\n    my $arg=@_[0];\n    my $ar\
g_list=@_[1];\n    my $keep_flag=@_[2];\n    #prot\
ect\n    $arg_list=~s/\\s-/ \\@/g;\n    $arg=~s/-/\
\\@/g;\n    \n    #search\n    if ($arg eq \"^\")\\
n      {\n	$arg_list=~/^([^@]*)/;\n	$val=$1;\n    \
  }\n    else\n      {$arg_list=~/$arg ([^@]*)/;$v\
al=$1;}\n    \n    #remove the parsed sequence if \
needed\n    if ($val && $keep_flag ne \"KEEP\")\n \
     {\n	if ( $arg eq \"^\"){$arg_list=~s/$val/ /;\
}\n	else {$arg_list=~s/($arg [^@]*)/ /;}\n      }\\
n	\n    #unprotect\n    $arg_list=~s/\\@/-/g;\n   \
 $arg=~s/\\@/-/g;\n    \n    return $val, $arg_lis\
t;\n  }\n\n","use Env;\nuse FileHandle;\nuse Cwd;\\
nuse File::Path;\nuse Sys::Hostname;\n\nour $PIDCH\
ILD;\nour $ERROR_DONE;\nour @TMPFILE_LIST;\nour $E\
XIT_FAILURE=1;\nour $EXIT_SUCCESS=0;\n\nour $REFDI\
R=getcwd;\nour $EXIT_SUCCESS=0;\nour $EXIT_FAILURE\
=1;\n\nour $PROGRAM=\"tc_generic_method.pl\";\nour\
 $CL=$PROGRAM;\n\nour $CLEAN_EXIT_STARTED;\nour $d\
ebug_lock=$ENV{\"DEBUG_LOCK\"};\nour $LOCKDIR=$ENV\
{\"LOCKDIR_4_TCOFFEE\"};\nif (!$LOCKDIR){$LOCKDIR=\
getcwd();}\nour $ERRORDIR=$ENV{\"ERRORDIR_4_TCOFFE\
E\"};\nour $ERRORFILE=$ENV{\"ERRORFILE_4_TCOFFEE\"\
};\n&set_lock ($$);\nif (isshellpid(getppid())){lo\
ck4tc(getppid(), \"LLOCK\", \"LSET\", \"$$\\n\");}\
\n      \n\n\n\n\nour $BLAST_MAX_NRUNS=2;\nour $CO\
MMAND;\nour $PIDCHILD;\n\n$REF_EMAIL=\"\";\n$tmp_d\
ir=\"\";\n$init_dir=\"\";\n\n\n$test=0;\nif ($test\
==1)\n  {\n    $SERVER=\"NCBI\";\n    $query=$ARGV\
[0];\n    $hitf=$ARGV[1];\n    %s=read_fasta_seq($\
query);\n    @sl=keys(%s);\n    &blast_xml2profile\
 (\"xx\", $s{$sl[0]}{seq},$maxid,$minid,$mincov, $\
hitf);\n    myexit ($EXIT_FAILURE);\n  }\n\nforeac\
h $v(@ARGV){$cl.=\"$v \";}\n$COMMAND=$cl;\n($mode)\
=&my_get_opt ( $cl, \"-mode=\",1,0);\n\n($A)=(&my_\
get_opt ( $cl, \"-name1=\",0,0));\n($B)=(&my_get_o\
pt ( $cl, \"-name2=\",0,0));\n($TMPDIR)=(&my_get_o\
pt ( $cl, \"-tmpdir=\",0,0));\n($CACHE)=(&my_get_o\
pt ( $cl, \"-cache=\",0,0));\n($SERVER)=((&my_get_\
opt ( $cl, \"-server=\",0,0)));\n($EMAIL)=((&my_ge\
t_opt ( $cl, \"-email=\",0,0)));\n\nif (!$A){$A=\"\
A\";}\nif (!$B){$B=\"B\";}\n\n\nif (!$TMPDIR)\n  {\
\n    $HOME=$ENV{HOME};\n    if ($ENV{TMP_4_TCOFFE\
E}){$TMPDIR=$ENV{TMP_4_TCOFFEE};}\n    else{$TMPDI\
R=\"$HOME/.t_coffee/tmp/\";}\n  }\nif ( ! -d $TMPD\
IR)\n  {\n    mkdir $TMPDIR;\n  }\nif ( ! -d $TMPD\
IR)\n  {\n    print \"ERROR: Could not create temp\
orary dir: $TMPDIR\\n\";\n    myexit ($EXIT_FAILUR\
E);\n  }\n\n$EMAIL=~s/XEMAILX/\\@/g;\nif (!$EMAIL)\
\n  {\n    if ($ENV{EMAIL_4_TCOFFEE}){$EMAIL=$ENV{\
EMAIL_4_TCOFFEE};}\n    elsif ($ENV{EMAIL}){$EMAIL\
=$ENV{EMAIL};}\n    else {$EMAIL=$REF_EMAIL;}\n  }\
\n\n($maxid,$minid,$mincov)=(&my_get_opt ( $cl, \"\
-maxid=\",0,0, \"-minid=\",0,0,\"-mincov=\",0,0));\
\nif (!$cl=~/\\-maxid\\=/){$maxid=95;}\nif (!$cl=~\
/\\-minid\\=/){$minid=35;}\nif (!$cl=~/\\-mincov\\\
=/){$mincov=80;}\n\n\n\n\nif ($mode eq \"seq_msa\"\
)\n  {\n    &seq2msa($mode,&my_get_opt ( $cl, \"-i\
nfile=\",1,1, \"-method=\",1,2, \"-param=\",0,0,\"\
-outfile=\",1,0, \"-database=\",0,0));\n  }\nelsif\
 ( $mode eq \"tblastx_msa\")\n  {\n    &seq2tblast\
x_lib ($mode,&my_get_opt ( $cl, \"-infile=\",1,1, \
\"-outfile=\",1,0));\n  }\nelsif ( $mode eq \"tbla\
stpx_msa\")\n  {\n    &seq2tblastpx_lib ($mode,&my\
_get_opt ( $cl, \"-infile=\",1,1, \"-outfile=\",1,\
0));\n  }\nelsif ( $mode eq \"thread_pair\")\n  {\\
n    &seq2thread_pair($mode,&my_get_opt ( $cl, \"-\
infile=\",1,1, \"-pdbfile1=\",1,1, \"-method=\",1,\
2,\"-param=\",0,0, \"-outfile=\",1,0, ));\n  }\nel\
sif ( $mode eq \"pdbid_pair\")\n  {\n    &seq2pdbi\
d_pair($mode,&my_get_opt ( $cl, \"-pdbfile1=\",1,0\
, \"-pdbfile2=\",1,0, \"-method=\",1,2,\"-param=\"\
,0,0, \"-outfile=\",1,0, ));\n  }\nelsif ( $mode e\
q \"pdb_pair\")\n  {\n    &seq2pdb_pair($mode,&my_\
get_opt ( $cl, \"-pdbfile1=\",1,1, \"-pdbfile2=\",\
1,1, \"-method=\",1,2,\"-param=\",0,0, \"-outfile=\
\",1,0, ));\n  }\nelsif ( $mode eq \"profile_pair\\
")\n  {\n     &seq2profile_pair($mode,&my_get_opt \
( $cl, \"-profile1=\",1,1, \"-profile2=\",1,1, \"-\
method=\",1,2,\"-param=\",0,0, \"-outfile=\",1,0 )\
);\n  }\nelsif ( $mode eq \"pdb_template\")\n  {\n\
    &blast2pdb_template ($mode,&my_get_opt ( $cl, \
\"-infile=\",1,1, \"-database=\",1,0, \"-method=\"\
,1,0, \"-outfile=\",1,0,\"-pdb_type=\",1,0));\n  }\
\nelsif ( $mode eq \"profile_template\")\n  {\n   \
 \n    &psiblast2profile_template ($mode,&my_get_o\
pt ( $cl, \"-infile=\",1,1, \"-database=\",1,0, \"\
-method=\",1,0, \"-outfile=\",1,0));\n  }\nelsif (\
 $mode eq \"psiprofile_template\")\n  {\n    &psib\
last2profile_template ($mode,&my_get_opt ( $cl, \"\
-infile=\",1,1, \"-database=\",1,0, \"-method=\",1\
,0, \"-outfile=\",1,0));\n  }\nelsif ( $mode eq \"\
RNA_template\")\n  {\n    &seq2RNA_template ($mode\
,&my_get_opt ( $cl, \"-infile=\",1,1, \"-outfile=\\
",1,0));\n  }\nelsif ( $mode eq \"tm_template\")\n\
  {\n    &seq2tm_template ($mode, \"\", &my_get_op\
t ( $cl, \"-infile=\",1,1,\"-arch=\",1,1,\"-psv=\"\
,1,1, \"-outfile=\",1,0,));\n  }\nelsif ( $mode eq\
 \"psitm_template\")\n  {\n    &seq2tm_template ($\
mode,&my_get_opt ( $cl, \"-database=\",1,0, \"-inf\
ile=\",1,1, \"-arch=\",1,1,\"-psv=\",1,1, \"-outfi\
le=\",1,0,));\n  }\nelsif ( $mode eq \"ssp_templat\
e\")\n  {\n    &seq2ssp_template ($mode,&my_get_op\
t ( $cl, \"-infile=\",1,1,\"-seq=\",1,1,\"-obs=\",\
1,1, \"-outfile=\",1,0));\n  }\nelsif ( $mode eq \\
"psissp_template\")\n  {\n    &seq2ssp_template ($\
mode,&my_get_opt ( $cl, \"-infile=\",1,1,\"-seq=\"\
,1,1,\"-obs=\",1,1, \"-outfile=\",1,0));\n  }\n\ne\
lsif ( $mode eq \"rna_pair\")\n{\n    &seq2rna_pai\
r($mode,&my_get_opt ( $cl, \"-pdbfile1=\",1,1, \"-\
pdbfile2=\",1,1, \"-method=\",1,2,\"-param=\",0,0,\
 \"-outfile=\",1,0, ));\n}\nelsif ( $mode eq \"cal\
c_rna_template\")\n{\n    &calc_rna_template($mode\
,&my_get_opt ( $cl, \"-infile=\",1,1,\"-pdbfile=\"\
,1,1, \"-outfile=\",1,0));\n}\nelse\n  {\n    myex\
it(flush_error( \"$mode is an unknown mode of tc_g\
eneric_method.pl\"));\n  }\nmyexit ($EXIT_SUCCESS)\
;\n\n\nsub seq2ssp_template\n  {\n  my ($mode, $in\
file,$gor_seq,$gor_obs,$outfile)=@_;\n  my %s, %h;\
\n  my $result;\n  my (@profiles);\n  &set_tempora\
ry_dir (\"set\",$infile,\"seq.pep\");\n  %s=read_f\
asta_seq (\"seq.pep\");\n\n  \n  open (R, \">resul\
t.aln\");\n  \n  #print stdout \"\\n\";\n  foreach\
 $seq (keys(%s))\n    {\n      \n      open (F, \"\
>seqfile\");\n      $s{$seq}{seq}=uc$s{$seq}{seq};\
\n      print (F \">$s{$seq}{name}\\n$s{$seq}{seq}\
\\n\");\n      close (F);\n      $lib_name=\"$s{$s\
eq}{name}.ssp\";\n      $lib_name=&clean_file_name\
 ($lib_name);\n      \n      if ($mode eq \"ssp_te\
mplate\"){&seq2gor_prediction ($s{$seq}{name},$s{$\
seq}{seq}, \"seqfile\", $lib_name,$gor_seq, $gor_o\
bs);}\n      elsif ($mode eq \"psissp_template\")\\
n	{\n	  &seq2msa_gor_prediction ($s{$seq}{name},$s\
{$seq}{seq},\"seqfile\", $lib_name,$gor_seq, $gor_\
obs);\n	}\n    \n      if ( !-e $lib_name)\n	{\n	 \
 myexit(flush_error(\"GORIV failed to compute the \
secondary structure of $s{$seq}{name}\"));\n	  mye\
xit ($EXIT_FAILURE);\n	}\n      else\n	{\n	  print\
 stdout \"\\tProcess: >$s{$seq}{name} _E_ $lib_nam\
e \\n\";\n	  print R \">$s{$seq}{name} _E_ $lib_na\
me\\n\";\n	}\n      unshift (@profiles, $lib_name)\
;\n    }\n  close (R);\n  &set_temporary_dir (\"un\
set\",$mode, $method,\"result.aln\",$outfile, @pro\
files);\n}\n\nsub seq2tm_template\n  {\n  my ($mod\
e, $db, $infile,$arch,$psv,$outfile)=@_;\n  my %s,\
 %h;\n  my $result;\n  my (@profiles);\n  &set_tem\
porary_dir (\"set\",$infile,\"seq.pep\");\n  %s=re\
ad_fasta_seq (\"seq.pep\");\n\n  \n  open (R, \">r\
esult.aln\");\n  \n  #print stdout \"\\n\";\n  for\
each $seq (keys(%s))\n    {\n      open (F, \">seq\
file\");\n      print (F \">$s{$seq}{name}\\n$s{$s\
eq}{seq}\\n\");\n      close (F);\n      $lib_name\
=\"$s{$seq}{name}.tmp\";\n      $lib_name=&clean_f\
ile_name ($lib_name);\n\n      if ($mode eq \"tm_t\
emplate\")\n	{\n	  &safe_system (\"t_coffee -other\
_pg fasta_seq2hmmtop_fasta.pl -in=seqfile -out=$li\
b_name -arch=$arch -psv=$psv\");\n	}\n      elsif \
( $mode eq \"psitm_template\")\n	{\n	  &seq2msa_tm\
_prediction ($s{$seq}{name},$s{$seq}{seq}, $db, \"\
seqfile\", $lib_name,$arch, $psv);\n	}\n      if (\
 !-e $lib_name)\n	{\n	  myexit(flush_error(\"RNApl\
fold failed to compute the secondary structure of \
$s{$seq}{name}\"));\n	  myexit ($EXIT_FAILURE);\n	\
}\n      else\n	{\n	  print stdout \"\\tProcess: >\
$s{$seq}{name} _T_ $lib_name\\n\";\n	  print R \">\
$s{$seq}{name} _T_ $lib_name\\n\";\n	}\n      unsh\
ift (@profiles, $lib_name);\n    }\n  close (R);\n\
  &set_temporary_dir (\"unset\",$mode, $method,\"r\
esult.aln\",$outfile, @profiles);\n}\n\nsub seq2RN\
A_template\n  {\n  my ($mode, $infile,$outfile)=@_\
;\n  my %s, %h, ;\n  my $result;\n  my (@profiles)\
;\n  &set_temporary_dir (\"set\",$infile,\"seq.pep\
\");\n  %s=read_fasta_seq (\"seq.pep\");\n\n  \n  \
open (R, \">result.aln\");\n  \n  #print stdout \"\
\\n\";\n  foreach $seq (keys(%s))\n    {\n      op\
en (F, \">seqfile\");\n      print (F \">$s{$seq}{\
name}\\n$s{$seq}{seq}\\n\");\n      close (F);\n  \
    $lib_name=\"$s{$seq}{name}.rfold\";\n      $li\
b_name=&clean_file_name ($lib_name);\n      &safe_\
system (\"t_coffee -other_pg RNAplfold2tclib.pl -i\
n=seqfile -out=$lib_name\");\n      \n      if ( !\
-e $lib_name)\n	{\n	 myexit(flush_error(\"RNAplfol\
d failed to compute the secondary structure of $s{\
$seq}{name}\"));\n	  myexit ($EXIT_FAILURE);\n	}\n\
      else\n	{\n	  print stdout \"\\tProcess: >$s{\
$seq}{name} _F_ $lib_name\\n\";\n	  print R \">$s{\
$seq}{name} _F_ $lib_name\\n\";\n	}\n      unshift\
 (@profiles, $lib_name);\n    }\n  close (R);\n  &\
set_temporary_dir (\"unset\",$mode, $method,\"resu\
lt.aln\",$outfile, @profiles);\n}\n\nsub psiblast2\
profile_template \n  {\n  my ($mode, $infile, $db,\
 $method, $outfile)=@_;\n  my %s, %h, ;\n  my ($re\
sult,$psiblast_output,$profile_name,@profiles);\n \
 my $trim=0;\n  &set_temporary_dir (\"set\",$infil\
e,\"seq.pep\");\n  %s=read_fasta_seq (\"seq.pep\")\
;\n  open (R, \">result.aln\");\n  \n  #print stdo\
ut \"\\n\";\n  foreach $seq (keys(%s))\n    {\n   \
   open (F, \">seqfile\");\n      print (F \">$A\\\
n$s{$seq}{seq}\\n\");\n      close (F);\n      $ps\
iblast_output=&run_blast ($s{$seq}{name},$method, \
$db, \"seqfile\",\"outfile\");\n      \nif ( -e $p\
siblast_output)\n	{\n	  %profile=blast_xml2profile\
($s{$seq}{name}, $s{$seq}{seq},$maxid, $minid,$min\
cov,$psiblast_output);\n	  unlink ($psiblast_outpu\
t);\n	  \n	  $profile_name=\"$s{$seq}{name}.prf\";\
\n	  $profile_name=&clean_file_name ($profile_name\
);\n	  unshift (@profiles, $profile_name);\n	  out\
put_profile ($profile_name, \\%profile, $trim);\n	\
  print stdout \"\\tProcess: >$s{$seq}{name} _R_ $\
profile_name [$profile{n} Seq.] [$SERVER/blast/$db\
][$CACHE_STATUS]\\n\";\n	  print R \">$s{$seq}{nam\
e} _R_ $profile_name\\n\";\n	}\n    }\n  close (R)\
;\n  &set_temporary_dir (\"unset\",$mode, $method,\
\"result.aln\",$outfile, @profiles);\n}\n\nsub bla\
st2pdb_template \n  {\n  my ($mode, $infile, $db, \
$method, $outfile,$type)=@_;\n  my %s, %h, ;\n  my\
 ($result,$blast_output);\n  &set_temporary_dir (\\
"set\",$infile,\"seq.pep\");\n  %s=read_fasta_seq \
(\"seq.pep\");\n  open (R, \">result.aln\");\n  \n\
 \n  #print stdout \"\\n\";\n  foreach $seq (keys(\
%s))\n    {\n      my $c;\n      my $found;\n     \
 \n      open (F, \">seqfile\");\n      print (F \\
">$A\\n$s{$seq}{seq}\\n\");\n      close (F);\n   \
  \n      $blast_output=&run_blast ($s{$seq}{name}\
,$method, $db, \"seqfile\",\"outfile\");\n     \n \
     %p=blast_xml2profile($s{$seq}{name}, $s{$seq}\
{seq},$maxid, $minid,$mincov,$blast_output);\n    \
  unlink ($blast_output);\n      \n      $c=1;\n  \
    print stdout \"\\tProcess: >$s{$seq}{name} [$S\
ERVER/blast/$db][$CACHE_STATUS]\\n\";\n      while\
 (!$found && $c<$p{n})\n	{\n	  $pdbid=&id2pdbid($p\
{$c}{identifyer});\n	  if ( length ($pdbid)>5){$pd\
bid=id2pdbid($p{$c}{definition});}\n\n	  if ( leng\
th ($pdbid)>5)\n	    {\n	      myexit(add_error (E\
XIT_FAILURE,$$,$$,getppid(), \"BLAST_FAILURE::Coul\
d Not Parse PDBID ($p{$c}{identifyer},$p{$c}{defin\
ition})\"));\n	    }\n	  \n\n	  if (!&pdb_is_relea\
sed($pdbid))\n	    {\n	      print stdout \"\\t\\t\
**$pdbid [PDB NOT RELEASED or WITHDRAWN]\\n\";\n	 \
     $c++;\n	    }\n	  elsif (!&pdb_has_right_type\
 ($pdbid,$type))\n	    {\n	      my $ptype=&pdb2ty\
pe ($pdbid);\n	      my $etype=&type2etype($type);\
\n	      \n	      print stdout \"\\t\\t**$pdbid [$\
ptype cannot be used (expected: $etype)]\\n\";\n	 \
     $c++;\n	    }\n	  else\n	    {\n	      $found\
=1;\n	    }\n	}\n\n      if ($found)\n	{\n	  print\
 R \">$s{$seq}{name} _P_ $pdbid\\n\";\n	  print st\
dout \"\\t\\t >$s{$seq}{name} _P_ $pdbid\\n\";\n	}\
\n      else\n	{\n	  print R \">$s{$seq}{name}\\n\\
";\n	  print stdout \"\\t\\t >$s{$seq}{name} No Te\
mplate Selected\\n\";\n	}\n    }\n  close (R);\n  \
&set_temporary_dir (\"unset\",$mode, $method,\"res\
ult.aln\",$outfile);\n}\nsub type2etype\n  {\n    \
my $type=shift;\n    my $etype;\n    \n    if ( $t\
ype=~/n/){$etype.=\"NMR \";}\n    if ( $type=~/d/)\
{$etype.=\"diffraction \";}\n    if ( $type=~/m/){\
$etype.=\"model \";}\n    return $etype;\n  }\nsub\
 pdb2type\n  {\n     my $pdb=shift;\n     my $f=vt\
mpnam();\n     \n     my $value= &safe_system (\"t\
_coffee -other_pg extract_from_pdb -model_type $pd\
b > $f\");\n     my $r=&file2string ($f);\n     ch\
omp($r);\n     return $r;\n   }\nsub pdb_has_right\
_type\n  {\n    my $pdb=shift;\n    my $type=shift\
;\n    \n    my $f=vtmpnam();\n    \n    my $value\
= &safe_system (\"t_coffee -other_pg extract_from_\
pdb -model_type $pdb > $f\");\n    my $r=&file2str\
ing ($f);\n    chomp($r);\n\n        \n    if ( $r\
 eq \"NMR\" && $type=~/n/){return 1;}\n    elsif (\
 $r eq \"diffraction\" && $type=~/d/){return 1;}\n\
    elsif ( $r eq \"model\" && $type=~/m/){return \
1;}\n    else {return 0;}\n  }\nsub pdb_is_release\
d\n  {\n    my $pdb=shift;\n    my $f=vtmpnam();\n\
    \n    $value= &safe_system (\"t_coffee -other_\
pg extract_from_pdb -is_released_pdb_name $pdb > $\
f\");\n    my $r=&file2string ($f);\n    chomp($r)\
;\n    return $r;\n  }\nsub blast_msa\n  {\n    my\
 ($infile,$db,$outfile)=@_;\n    my ($a, %seq);\n \
   my $seqfile;\n    my $SEQ=new FileHandle;\n    \
my $seqfile=\"seqfile\";\n    my @txt;\n    \n    \
\n    %s1=&read_fasta_seq ($db);\n    \n    foreac\
h $s (keys (%s1))\n      {\n	$i=$s1{$s}{order};\n	\
$s{$i}{name}=$s;\n	$s{$i}{seq}=$s1{$s}{seq};\n	$s{\
$i}{len}=length( $s{$i}{seq});\n	$s{n}++;\n      }\
\n    \n    &safe_system (\"formatdb -i $db\");\n \
   &safe_system  (\"blastall -i $infile -d $db -m7\
 -p blastp -o io\");\n    &set_blast_type (\"io\")\
;\n    \n    %FB=&xml2tag_list (\"io\", \"Iteratio\
n\");\n    open (F, \">$outfile\");\n    print F \\
"! TC_LIB_FORMAT_01\\n\";\n    print F \"$s{n}\\n\\
";\n    for ( $a=0; $a<$s{n}; $a++)\n      {\n	pri\
nt F \"$s{$a}{name} $s{$a}{len} $s{$a}{seq}\\n\";\\
n      }\n\n\n    for ( $a=0; $a<$FB{n}; $a++)\n  \
    {\n	%p=blast_xml2profile ($s{$a}{name}, $s{$a}\
{seq},100, 0, 0, $FB{$a}{body});\n	my $query=$p{0}\
{name};\n	my $i= $s1{$query}{order}+1;\n	for ($b=1\
; $b<$p{n}; $b++)\n	  {\n	    my $l=length ($p{$b}\
{Qseq});\n	    my $hit=$p{$b}{definition};\n	    m\
y $Qstart=$p{$b}{Qstart};\n	    my $Hstart=$p{$b}{\
Hstart};\n	    my $identity=$p{$b}{identity};\n	  \
  my @lrQ=split (//,$p{$b}{Qseq});\n	    my @lrH=s\
plit (//,$p{$b}{Hseq});\n	    \n	    my $j= $s1{$h\
it}{order}+1;\n	    #if ( $j==$i){next;}\n	    pri\
ntf F \"# %d %d\\n\", $i, $j;\n	    #  print  F \"\
\\n$p{$b}{Qseq} ($Qstart)\\n$p{$b}{Hseq} ($Hstart)\
\";\n	    for ($c=0; $c<$l; $c++)\n	      {\n		my \
$rQ=$lrQ[$c];\n		my $rH=$lrH[$c];\n		my $n=0;\n		\\
n		if ($rQ ne \"-\"){$n++, $Qstart++;}\n		if ($rH \
ne \"-\"){$n++; $Hstart++;}\n		\n		if ( $n==2)\n		\
  {\n		    printf F \"\\t%d %d %d\\n\", $Qstart-1,\
 $Hstart-1,$identity;\n		  }\n	      }\n	  }\n    \
  }\n    print F \"! SEQ_1_TO_N\\n\";\n    close (\
F);\n    return $output;\n  \n  }\n\nsub blast_msa\
_old\n  {\n    my ($infile,$outfile)=@_;\n    my (\
$a, %seq);\n    %s1=&read_fasta_seq ($infile);\n  \
  foreach $s (keys (%s1))\n      {\n	$i=$s1{$s}{or\
der};\n	$s{$i}{name}=$s;\n	$s{$i}{seq}=$s1{$s}{seq\
};\n	$s{$i}{len}=length( $s{$i}{seq});\n	$s{n}++;\\
n      }\n    &safe_system (\"formatdb -i $infile\\
");\n    &safe_system (\"blastall -i $infile -d $i\
nfile -m7 -o io\");\n    &set_blast_type (\"io\");\
\n    \n    %FB=&xml2tag_list (\"io\", \"Iteration\
\");\n    \n    open (F, \">$outfile\");\n    prin\
t F \"! TC_LIB_FORMAT_01\\n\";\n    print F \"$s{n\
}\\n\";\n    for ( $a=0; $a<$s{n}; $a++)\n      {\\
n	print F \"$s{$a}{name} $s{$a}{len} $s{$a}{seq}\\\
n\";\n      }\n    for ( $a=0; $a<$FB{n}; $a++)\n \
     {\n	%p=blast_xml2profile ($s{$a}{name}, $s{$a\
}{seq},100, 0, 0, $FB{$a}{body});\n	for ($b=1; $b<\
$p{n}; $b++)\n	  {\n	    my $l=length ($p{$b}{Qseq\
});\n	    my $hit=$p{$b}{definition};\n	    my $Qs\
tart=$p{$b}{Qstart};\n	    my $Hstart=$p{$b}{Hstar\
t};\n	    my $identity=$p{$b}{identity};\n	    my \
@lrQ=split (//,$p{$b}{Qseq});\n	    my @lrH=split \
(//,$p{$b}{Hseq});\n	    my $i= $s1{$s{$a}{name}}{\
order}+1;\n	    my $j= $s1{$hit}{order}+1;\n	    #\
if ( $j==$i){next;}\n	    printf F \"# %d %d\\n\",\
 $i, $j;\n	    #  print  F \"\\n$p{$b}{Qseq} ($Qst\
art)\\n$p{$b}{Hseq} ($Hstart)\";\n	    for ($c=0; \
$c<$l; $c++)\n	      {\n		my $rQ=$lrQ[$c];\n		my $\
rH=$lrH[$c];\n		my $n=0;\n		\n		if ($rQ ne \"-\"){\
$n++, $Qstart++;}\n		if ($rH ne \"-\"){$n++; $Hsta\
rt++;}\n		\n		if ( $n==2)\n		  {\n		    printf F \\
"\\t%d %d %d\\n\", $Qstart-1, $Hstart-1,$identity;\
\n		  }\n	      }\n	  }\n      }\n    print F \"! \
SEQ_1_TO_N\\n\";\n    close (F);\n    return $outp\
ut;\n  \n  }\n\nsub seq2msa\n  {\n    my ($mode, $\
infile, $method, $param, $outfile,$database)=@_;\n\
    &set_temporary_dir (\"set\",$infile,\"seq.pep\\
", $database, \"db.pep\");\n    $param.=\" >/dev/n\
ull 2>&1 \";\n    \n    \n    #make sure test.pep \
is in FASTA\n    &safe_system (\"t_coffee -other_p\
g seq_reformat -in seq.pep -output fasta_seq > x\"\
);\n    `mv x seq.pep`;\n    \n    if ( $method eq\
 \"blastp\")\n      {\n	&blast_msa (\"seq.pep\", \\
"db.pep\",\"result.aln\");\n      }\n    elsif ( $\
method eq \"muscle\")\n      {\n	`muscle -in seq.p\
ep -out result.aln $param`;\n      }\n    elsif ( \
$method eq \"probcons\")\n      {\n	`probcons seq.\
pep >result.aln 2>/dev/null`;\n      }\n    elsif \
( $method eq \"mafft\")\n      {\n	`mafft --quiet \
--localpair --maxiterate 1000 seq.pep> result.aln \
 2>/dev/null`\n      }\n    elsif ( $method=~/pran\
k/)\n      {\n	`$method -d=seq.pep -o=result.aln -\
quiet 2>/dev/null`;\n	`mv result.aln.1.fas result.\
aln`;\n      }\n    else\n      {\n	`$method -infi\
le=seq.pep -outfile=result.aln`;\n      }\n    \n \
   &set_temporary_dir (\"unset\",$mode, $method,\"\
result.aln\",$outfile);\n    myexit ($EXIT_SUCCESS\
);\n  }\n\nsub seq2thread_pair\n  {\n    my ($mode\
, $infile, $pdbfile1, $method, $param, $outfile)=@\
_;\n    &set_temporary_dir (\"set\",$infile,\"seq.\
pep\",$pdbfile1,\"struc.pdb\");\n    if ($method e\
q \"fugueali\")\n      {\n	#Env Variable that need\
 to be defined for Fugue\n	if (!$ENV{FUGUE_LIB_LIS\
T}){$ENV{FUGUE_LIB_LIST}=\"DUMMY\";}\n	if (!$ENV{H\
OMSTRAD_PATH})  {$ENV{HOMSTRAD_PATH}=\"DUMMY\";}\n\
	if (!$ENV{HOMS_PATH}){$ENV{HOMS_PATH}=\"DUMMY\";}\
\n	\n	`joy struc.pdb >x 2>x`;\n	&check_file(\"stru\
c.tem\", \"Joy failed [FATAL:$PROGRAM/$method]\");\
\n	`melody -t struc.tem >x 2>x`;\n	&check_file(\"s\
truc.tem\", \"Melody failed [FATAL:$PROGRAM/$metho\
d]\");\n	`fugueali -seq seq.pep -prf struc.fug -pr\
int > tmp_result.aln`;\n	\n	&check_file(\"tmp_resu\
lt.aln\", \"Fugue failed [FATAL:$PROGRAM/$method]\\
");\n	&safe_system (\"t_coffee -other_pg seq_refor\
mat -in tmp_result.aln -output fasta_aln >result.a\
ln\");\n      }\n    elsif ( $method eq \"t_coffee\
\")\n      {\n	&safe_system (\"t_coffee -in Pstruc\
.pdb Sseq.pep Mslow_pair -outfile result.aln -quie\
t\");\n      }\n    else\n      {\n	&safe_system (\
\"$method -infile=seq.pep -pdbfile1=struc.pdb -out\
file=result.aln $param>x 2>x\");\n      }\n    &se\
t_temporary_dir (\"unset\",$mode,$method,\"result.\
aln\",$outfile);\n    myexit ($EXIT_SUCCESS);\n  }\
\nsub seq2pdbid_pair\n  {\n    my ($mode, $pdbfile\
1, $pdbfile2, $method, $param, $outfile)=@_;\n    \
my ($name);\n\n    \n    &set_temporary_dir (\"set\
\");\n    $name=$pdbfile1.\" \".$pdbfile2;\n\n    \
if (    &cache_file(\"GET\",\"\",\"$name\",\"$meth\
od\",\"dali\",$outfile,\"EBI\"))\n      {return $o\
utfile;}\n    else\n      {\n	if ($method eq \"dal\
iweb\")\n	  {\n	    $pdbfile1=~/(....)(.)/;\n	    \
$id1=$1; $c1=$2;\n	    \n	    $pdbfile2=~/(....)(.\
)/;\n	    $id2=$1; $c2=$2;\n	    \n	    $command=\\
"t_coffee -other_pg dalilite.pl --pdb1 $id1 --chai\
nid1 $c1 --pdb2 $id2 --chainid2 $c2 --email=$EMAIL\
  >dali_stderr 2>dali_stderr\";\n	    $dali=`$comm\
and`;\n	    \n	    open (F, \"dali_stderr\");\n	  \
  while (<F>)\n	      {\n		if ( /JobId: dalilite-(\
\\S+)/)\n		{\n		  $jobid=$1;\n		}\n	      }\n	    \
close (F);\n	    unlink (\"dali_stderr\");\n	    \\
n	    $output1=\"dalilite-$jobid.txt\";\n	    if (\
 -e $output1)\n	      {\n		unlink ($output1);\n		&\
url2file (\"http://www.ebi.ac.uk/Tools/es/cgi-bin/\
jobresults.cgi/dalilite/dalilite-$jobid/aln.html\"\
, \"output2\");\n		\n		if ( -e \"output2\")\n		  {\
\n		    my ($seq1, $seq2);\n		    $seq1=$seq2=\"\"\
;\n		    \n		    open (F, \"output2\");\n		    whi\
le (<F>)\n		      {\n			$l=$_;\n			if ( $l=~/Query\
\\s+(\\S+)/)\n			  {\n			    $seq1.=$1;\n			  }\n	\
		elsif ( $l=~/Sbjct\\s+(\\S+)/)\n			  {\n			    $\
seq2.=$1;\n			  }\n		      }\n		    close (F);\n		\
    unlink (\"output2\");\n		    if ($seq1 ne \"\"\
 && $seq2 ne \"\")\n		      {\n			$output3=\">$A\\\
n$seq1\\n>$B\\n$seq2\\n\";\n			$output3=~s/\\./-/g\
;\n			open (F, \">result.aln\");\n			print F \"$ou\
tput3\";\n			close (F);\n		      }\n		  }\n	      \
}\n	  }\n      }\n    &cache_file(\"SET\",\"\",\"$\
name\",\"$method\",\"dali\",\"result.aln\",\"EBI\"\
);\n    &set_temporary_dir (\"unset\",$mode, $meth\
od, \"result.aln\",$outfile);\n    myexit ($EXIT_S\
UCCESS);\n  }\nsub seq2pdb_pair\n  {\n    my ($mod\
e, $pdbfile1, $pdbfile2, $method, $param, $outfile\
)=@_;\n    \n    &set_temporary_dir (\"set\",$pdbf\
ile1,\"pdb1.pdb\",$pdbfile2,\"pdb2.pdb\");\n    if\
 ($method eq \"t_coffee\")\n      {\n	&safe_system\
 (\"t_coffee -in Ppdb1.pdb Ppdb2.pdb -quiet -outfi\
le=result.aln\");\n      }\n    elsif ( $method eq\
 \"DaliLite\")\n      {\n	if ( &safe_system (\"Dal\
iLite -pairwise pdb1.pdb pdb2.pdb >tmp1\")==$EXIT_\
SUCCESS)\n	  {\n	     my ($seq1, $seq2);\n	     $s\
eq1=$seq2=\"\";\n		    \n	     open (F, \"tmp1\");\
\n	     while (<F>)\n	       {\n		 $l=$_;\n		 if (\
 $l=~/Query\\s+(\\S+)/)\n		   {\n		     $seq1.=$1;\
\n		   }\n		 elsif ( $l=~/Sbjct\\s+(\\S+)/)\n		   \
{\n		     $seq2.=$1;\n		   }\n	       }\n	     clo\
se (F);\n	     unlink (\"tmp1\");\n	     if ($seq1\
 ne \"\" && $seq2 ne \"\")\n	       {\n		 my $outp\
ut3=\">$A\\n$seq1\\n>$B\\n$seq2\\n\";\n		 $output3\
=~s/\\./-/g;\n		 open (F, \">result.aln\");\n		 pr\
int F \"$output3\";\n		 close (F);\n	       }\n	  \
 }\n	else\n	  {\n	    print \"ERROR: DalLite faile\
d to align the considered structures[tc_generic_me\
thod.pl]\\n\";\n	  }    \n      }\n    elsif ( $me\
thod eq \"TMalign\")\n      {\n	if ( &safe_system \
(\"TMalign pdb1.pdb pdb2.pdb >tmp1\")==$EXIT_SUCCE\
SS)\n	  {\n	    `tail -4 tmp1 > tmp2`;\n	    \n	  \
  open (F, \"tmp2\");\n	    while (<F>)\n	      {\\
n		unshift(@l, $_);\n	      }\n	    close (F);\n	 \
   open (F, \">result.aln\");\n	    $l[3]=~s/[^a-z\
A-Z0-9-]/\\-/g;\n	    $l[1]=~s/[^a-zA-Z0-9-]/\\-/g\
;\n	    print F \">$A\\n$l[3]\\n>$B\\n$l[1]\\n\";\\
n	    close (F);\n	  }\n	else\n	  {\n	    print \"\
ERROR: TMalign failed to align the considered stru\
ctures[tc_generic_method.pl]\\n\";\n	    `rm resul\
t.aln >/dev/null 2>/dev/null`;\n	  }\n      }\n   \
 elsif ( $method eq \"mustang\")\n      {\n	if ( &\
safe_system (\"mustang -i pdb1.pdb pdb2.pdb -F fas\
ta >/dev/null 2>/dev/null\")==$EXIT_SUCCESS)\n	  {\
\n	    `mv results.afasta result.aln`;\n	  }\n	els\
e\n	  {\n	    print \"ERROR: mustang failed to ali\
gn the considered structures[tc_generic_method.pl]\
\\n\";\n	    `rm result.aln >/dev/null 2>/dev/null\
`;\n	  }\n      }\n    else\n      {\n	if ( &safe_\
system (\"$method -pdbfile1=pdb1.pep -pdbfile2=pdb\
2.pdb -outfile=result.aln $param>x 2>x\")==$EXIT_S\
UCCESS)\n	  {\n	    `mv results.afasta result.aln`\
;\n	  }\n	else\n	  {\n	    print \"ERROR: $method \
failed to align the considered structures[tc_gener\
ic_method.pl]\\n\";\n	    `rm result.aln >/dev/nul\
l 2>/dev/null`;\n	  }\n      }\n    &set_temporary\
_dir (\"unset\",$mode, $method, \"result.aln\",$ou\
tfile);\n    myexit ($EXIT_SUCCESS);\n  }\n\nsub s\
eq2profile_pair\n  {\n    my ($mode, $profile1, $p\
rofile2, $method, $param, $outfile)=@_;\n    \n   \
 \n    if ($method eq \"clustalw\")\n      {\n	&se\
t_temporary_dir (\"set\",$profile1,\"prf1.aln\",$p\
rofile2,\"prf2.aln\");\n	`clustalw -profile1=prf1.\
aln -profile2=prf2.aln -outfile=result.aln`;\n	&se\
t_temporary_dir (\"unset\",$mode, $method, \"resul\
t.aln\",$outfile);\n      }\n    elsif ( $method e\
q \"hhalign\")\n      {\n	hhalign ( $profile1,$pro\
file2,$outfile,$param);\n      }\n    else\n      \
{\n	\n	`$method -profile1=prf1.aln -profile2=prf2.\
aln -outfile=result.aln $param>x 2>x`;\n      }\n \
  \n    myexit ($EXIT_SUCCESS);\n  }\n\nsub pg_is_\
installed\n  {\n    my @ml=@_;\n    my ($r, $p, $m\
);\n    my $supported=0;\n    \n    my $p=shift (@\
ml);\n    if ($p=~/::/)\n      {\n	if (safe_system\
 (\"perl -M$p -e 1\")==$EXIT_SUCCESS){return 1;}\n\
	else {return 0;}\n      }\n    else\n      {\n	$r\
=`which $p 2>/dev/null`;\n	if ($r eq \"\"){$r=0;}\\
n	else {$r=1;}\n	\n	if ($r==0 && is_blast_package \
($p)){return pg_is_installed (\"legacy_blast.pl\")\
;}\n	else {return $r;}\n      }\n  }\n\nsub is_bla\
st_package\n  {\n    my $p=shift;\n    if ( $p=~/b\
lastp/){return 1;}\n    elsif ($p=~/blastall/){ret\
urn 1;}\n    elsif ($p=~/blastn/){return 1;}\n    \
elsif ($p=~/blastx/){return 1;}\n    elsif ($p=~/f\
ormatdb/){return 1;}\n    else {return 0;}\n  }\n \
   \nsub check_internet_connection\n  {\n    my $i\
nternet;\n    my $tmp;\n    &check_configuration (\
 \"wget\"); \n    \n    $tmp=&vtmpnam ();\n    \n \
   if     (&pg_is_installed    (\"wget\")){`wget w\
ww.google.com -O$tmp >/dev/null 2>/dev/null`;}\n  \
  elsif  (&pg_is_installed    (\"curl\")){`curl ww\
w.google.com -o$tmp >/dev/null 2>/dev/null`;}\n   \
 \n    if ( !-e $tmp || -s $tmp < 10){$internet=0;\
}\n    else {$internet=1;}\n    if (-e $tmp){unlin\
k $tmp;}\n\n    return $internet;\n  }\nsub check_\
pg_is_installed\n  {\n    my @ml=@_;\n    my $r=&p\
g_is_installed (@ml);\n    if (!$r && $p=~/::/)\n \
     {\n	print STDERR \"\\nYou Must Install the pe\
rl package $p on your system.\\nRUN:\\n\\tsudo per\
l -MCPAN -e 'install $pg'\\n\";\n      }\n    elsi\
f (!$r)\n      {\n	myexit(flush_error(\"\\nProgram\
 $p Supported but Not Installed on your system\"))\
;\n      }\n    else\n      {\n	return 1;\n      }\
\n  }\nsub set_temporary_dir\n  {\n    my @list=@_\
;\n    my $dir_mode, $a, $mode, $method;\n  \n    \
$dir_mode=shift (@list);\n\n    \n    if ( $dir_mo\
de eq \"set\")\n      {\n	$initial_dir=cwd();\n	if\
 ( !$tmp_dir)\n	  {\n	    $rand=rand (100000);\n	 \
   $tmp_dir=\"$TMPDIR/tmp4tcoffee_profile_pair_dir\
_$$\\_P_$rand\";\n	  }\n	if ( !-d $tmp_dir)\n	  {\\
n	    push (@TMPDIR_LIST, $tmp_dir);\n	    `mkdir \
$tmp_dir`;\n	  }\n	\n	for ( $a=0; $a<=$#list; $a+=\
2)\n	      {\n		if (-e $list[$a]){ `cp $list[$a] $\
tmp_dir/$list[$a+1]`;}\n	      }\n	chdir $tmp_dir;\
\n      }\n    elsif ( $dir_mode eq \"unset\")\n  \
    {\n	$mode=shift (@list);\n	$method=shift (@lis\
t);\n	\n	if (!-e $list[0])\n	  {\n	   myexit(flush\
_error(\"Program $method failed to produce $list[1\
]\" ));\n	    myexit ($EXIT_FAILURE);\n	  }\n	else\
\n	  {\n	    chdir $initial_dir;\n	    # `t_coffee\
 -other_pg seq_reformat -in $tmp_dir/$list[0] -out\
put fasta_aln -out $tmp_dir/result2.aln`;\n	    `c\
p $tmp_dir/$list[0] $tmp_dir/result2.aln`;\n	    i\
f ( $list[1] eq \"stdout\")\n	      {\n		open (F, \
\"$tmp_dir/result2.aln\");\n		while (<F>){print $_\
;}close(F);\n	      }\n	    else\n	      {\n		`mv \
$tmp_dir/result2.aln $list[1]`;\n	      }\n	    sh\
ift (@list); shift (@list);\n	    foreach $f (@lis\
t)\n	      {\n		if (-e (\"$tmp_dir/$f\")){`mv $tmp\
_dir/$f .`;}\n	      }\n	  }\n      }\n  }\n\n\n\n\
\nsub my_get_opt\n  {\n    my @list=@_;\n    my $c\
l, $a, $argv, @argl;\n    \n    @argl=();\n    $cl\
=shift @list;\n    for ( $a=0; $a<=$#list; $a+=3)\\
n      {\n	$option=$list[$a];\n	$optional=$list[$a\
+1];\n	$status=$list[$a+2];\n	$argv=\"\";\n	if ($c\
l=~/$option(\\S+)/){$argv=$1;}\n	@argl=(@argl,$arg\
v);\n	\n	\n	#$optional:0=>optional\n	#$optional:1=\
>must be set\n	#$status: 0=>no requirement\n	#$sta\
tus: 1=>must be an existing file\n	#$status: 2=>mu\
st be an installed package\n	\n\n	if ($optional==0\
){;}\n	elsif ( $optional==1 && $argv eq \"\")\n	  \
{\n	    myexit(flush_error( \"ERROR: Option $optio\
n must be set\"));\n	    myexit ($EXIT_FAILURE);\n\
	  }\n	if ($status==0){;}\n	elsif ($status ==1 && \
$argv ne \"\" && !-e $argv)\n	  {\n	    myexit(flu\
sh_error( \"File $argv must exist\"));\n	    myexi\
t ($EXIT_FAILURE);\n	  }\n	elsif ( $status==2 && $\
argv ne \"\" && &check_pg_is_installed ($argv)==0)\
\n	  {\n	    myexit(flush_error( \" $argv is not i\
nstalled\"));\n	    myexit ($EXIT_FAILURE);\n	  }\\
n      }\n\n    return @argl;\n    }\n\nsub check_\
file \n  {\n    my ($file, $msg)=@_;\n\n    if ( !\
-e $file)\n      {\n	myexit(flush_error(\"$msg\"))\
;\n      }\n    }\nsub hhalign\n  {\n    my ($aln1\
, $aln2, $outfile, $param)=@_;\n    my $h1, $h2;\n\
    \n    $h{0}{index}=0;\n    $h{1}{index}=1;\n  \
  \n    $h{0}{aln}=$aln1;\n    $h{1}{aln}=$aln2;\n\
\n   \n\n    %{$h{0}}=aln2psi_profile (%{$h{0}});\\
n    %{$h{1}}=aln2psi_profile (%{$h{1}});\n\n    $\
param=~s/#S/ /g;\n    $param=~s/#M/\\-/g;\n    $pa\
ram=~s/#E/\\=/g;\n    \n\n    \n    $command=\"hha\
lign -i $h{0}{a3m} -t $h{1}{a3m} -tc $outfile.tmp \
-rank 1 -mapt 0 $param\";\n    `$command`;\n    \n\
  #  `hhalign -i $h{0}{a3m} -t $h{1}{a3m} -tc $out\
file.tmp -rank 1 -mapt 0 -gapf 0.8 -gapg 0.8`;\n  \
  \n\n    # To run global use the following\n    \\
n    open (I, \"$outfile.tmp\");\n    open (O, \">\
$outfile\");\n    $h{0}{cons}=s/\\./x/g;\n    $h{1\
}{cons}=s/\\./x/g;\n\n    print O \"! TC_LIB_FORMA\
T_01\\n2\\n$h{0}{name} $h{0}{len} $h{0}{seq}\\n$h{\
1}{name} $h{1}{len} $h{1}{seq}\\n#1 2\\n\";\n    \\
n    while (<I>)\n      {\n	if (/(\\d+)\\s+(\\d+)\\
\s+(\\d+)/)\n	  {\n	    print O \"\\t$h{0}{$1}\\t$\
h{1}{$2}\\t$3\\n\";\n	  }\n      }\n    print O \"\
! SEQ_1_TO_N\\n\";\n\n    close (O);\n    close (I\
);\n  }\n\nsub aln2psi_profile\n  {\n    my (%h)=@\
_;\n    my ($aln,$i,$hv, $a, @c, $n);\n   \n    $i\
=$h{index};\n    $aln=$h{aln};\n\n    `cp $aln $$.\
hhh_aln`;\n    $command=\"t_coffee -other_pg seq_r\
eformat -in $aln -output hasch\";\n    $hv=`$comma\
nd`;chomp ($hv);\n    \n    $h{a2m}=\"$tmp/$hv.tmp\
4hhpred.a2m\";\n    $h{a3m}=\"$tmp/$hv.tmp4hhpred.\
a3m\";\n    if ( -e $h{a3m}){;}\n    else\n      {\
\n	`hhconsensus  -M 50 -i $h{aln} -oa2m $h{a2m}`;\\
n	if (!-e $h{a2m})\n	  {\n	    print STDERR \"Prog\
ram tc_generic_method.pl FAILED to run:\\n\\thhcon\
sensus  -M 50 -i $h{aln} -oa2m $h{a2m}\";\n	    my\
exit ($EXIT_FAILURE);\n	  }\n	\n	`hhconsensus  -M \
50 -i $h{aln} -oa3m $h{a3m}`;\n	if (!-e $h{a3m})\n\
	  {\n	    print STDERR \"Program tc_generic_metho\
d.pl FAILED to run:\\n\\thhconsensus  -M 50 -i $h{\
aln} -oa3m $h{a3m}\";\n	    myexit ($EXIT_FAILURE)\
;\n	  }\n       `buildali.pl $h{a3m} -n 1`;\n     \
 }\n    \n    \n    $h{a2m_seq}=`head -n 2 $h{a2m}\
 | grep -v \">\"`;chomp ($h{a2m_seq});\n    $h{a3m\
_seq}=`head -n 2 $h{a3m} | grep -v \">\"`;chomp ($\
h{a3m_seq});\n    $h{cons}=$h{a2m_seq};\n    $h{se\
q}=`head -n 2 $h{aln} | grep -v \">\"`;chomp ($h{s\
eq});\n    \n    \n\n    @c=split (//, $h{cons});\\
n    $h{len}=$#c+1;\n    for ($n=0,$a=0, $b=0; $a<\
$h{len};$a++)\n      {\n	if ( $c[$a]=~/[A-Z]/)\n	 \
 {\n	    $h{++$n}=++$b;\n\n	  }\n	elsif ( $c[$a]=~\
/[a-z\\.]/)\n	  {\n	    ++$b;\n	  }\n      }\n    \
\n    $name=`head -n 2 $h{aln} | grep \">\"`;\n   \
 $name=~/\\>(\\S+)/;\n    $h{name}=$1;\n    \n    \
`cp $h{a2m} $i.a2m`;\n    `cp $h{a3m} $i.a3m`;\n  \
  `cp $h{aln} $i.hh_aln`;\n    \n    return %h;\n \
 }\n\nsub read_fasta_seq \n  {\n    my $f=@_[0];\n\
    my %hseq;\n    my (@seq, @com, @name);\n    my\
 ($a, $s,$nseq);\n\n    open (F, $f);\n    while (\
<F>)\n      {\n	$s.=$_;\n      }\n    close (F);\n\
\n    \n    @name=($s=~/>(\\S*).*\\n[^>]*/g);\n   \
 \n    @seq =($s=~/>.*.*\\n([^>]*)/g);\n    @com =\
($s=~/>\\S*(.*)\\n([^>]*)/g);\n\n    \n    $nseq=$\
#name+1;\n    \n    for ($a=0; $a<$nseq; $a++)\n  \
    {\n	my $s;\n	my $n=$name[$a];\n	$hseq{$n}{name\
}=$n;\n	$seq[$a]=~s/[^A-Za-z]//g;\n	$hseq{$n}{orde\
r}=$a;\n	$hseq{$n}{seq}=$seq[$a];\n	$hseq{$n}{com}\
=$com[$a];\n	\n      }\n    return %hseq;\n  }\n\n\
sub file_contains \n  {\n    my ($file, $tag, $max\
)=(@_);\n    my ($n);\n    $n=0;\n    \n    if ( !\
-e $file && ($file =~/$tag/)) {return 1;}\n    els\
if ( !-e $file){return 0;}\n    else \n      {\n	o\
pen (FC, \"$file\");\n	while ( <FC>)\n	  {\n	    i\
f ( ($_=~/$tag/))\n	      {\n		close (FC);\n		retu\
rn 1;\n	      }\n	    elsif ($max && $n>$max)\n	  \
    {\n		close (FC);\n		return 0;\n	      }\n	    \
$n++;\n	  }\n      }\n    close (FC);\n    return \
0;\n  }\n	    \n	  \nsub file2string\n  {\n    my \
$f=@_[0];\n    my $string, $l;\n    open (F,\"$f\"\
);\n    while (<F>)\n      {\n\n	$l=$_;\n	#chomp (\
$l);\n	$string.=$l;\n      }\n    close (F);\n    \
$string=~s/\\r\\n//g;\n    $string=~s/\\n//g;\n   \
 return $string;\n  }\n\n\nsub my_get_opt\n  {\n  \
  my @list=@_;\n    my $cl, $a, $argv, @argl;\n   \
 \n    @argl=();\n    $cl=shift @list;\n    for ( \
$a=0; $a<=$#list; $a+=3)\n      {\n	$option=$list[\
$a];\n	$optional=$list[$a+1];\n	$status=$list[$a+2\
];\n	$argv=\"\";\n	if ($cl=~/$option(\\S+)/){$argv\
=$1;}\n	@argl=(@argl,$argv);\n	\n	\n	#$optional:0=\
>optional\n	#$optional:1=>must be set\n	#$status: \
0=>no requirement\n	#$status: 1=>must be an existi\
ng file\n	#$status: 2=>must be an installed packag\
e\n	\n\n	if ($optional==0){;}\n	elsif ( $optional=\
=1 && $argv eq \"\")\n	  {\n\n	    myexit(flush_er\
ror(\"Option $option must be set\"));\n	   \n	  }\\
n	if ($status==0){;}\n	elsif ($status ==1 && $argv\
 ne \"\" && !-e $argv)\n	  {\n	     myexit(flush_e\
rror(\"File $argv must exist\"));\n	   \n	  }\n	el\
sif ( $status==2 && $argv ne \"\" && &check_pg_is_\
installed ($argv)==0)\n	  {\n	    myexit(flush_err\
or(\"$argv is not installed\"));\n	   \n	  }\n    \
  }\n\n    return @argl;\n    }\n\nsub tag2value \\
n  {\n    \n    my $tag=(@_[0]);\n    my $word=(@_\
[1]);\n    my $return;\n    \n    $tag=~/$word=\"(\
[^\"]+)\"/;\n    $return=$1;\n    return $return;\\
n  }\n      \nsub hit_tag2pdbid\n  {\n    my $tag=\
(@_[0]);\n    my $pdbid;\n       \n    $tag=~/id=\\
"(\\S+)\"/;\n    $pdbid=$1;\n    $pdbid=~s/_//;\n \
   return $pdbid;\n  }\nsub id2pdbid\n  {\n    my \
$in=@_[0];\n    my $id;\n    \n    $in=~/(\\S+)/;\\
n    $id=$in;\n    $id=~s/PDB/pdb/g;\n    \n    if\
 ( $id=~/.*gnl\\|.*\\,(\\S*).*/){$id=$1;}\n    els\
if ($id =~/pdb(.*)/){$id=$1;}\n    \n    $id=~s/[:\
|��_]//g;\n    return $id;\n  }\nsub set_blast\
_type \n  {\n    my $file =@_[0];\n    if (&file_c\
ontains ($file,\"EBIApplicationResult\",100)){$BLA\
ST_TYPE=\"EBI\";}\n    elsif (&file_contains ($fil\
e,\"NCBI_BlastOutput\",100)) {$BLAST_TYPE=\"NCBI\"\
;}\n    else\n      {\n	$BLAST_TYPE=\"\";\n      }\
\n    return $BLAST_TYPE;\n  }\nsub is_valid_blast\
_xml\n    {\n      my $file=shift;\n      my $line\
;\n      \n      \n      if ( !-e $file) {return 0\
;}\n      $line=&file2tail ($file,100);\n      if \
( $line=~/<\\/EBIApplicationResult/ || $line=~/<\\\
/NCBI_BlastOutput/ || $line=~/<\\/BlastOutput/ ){r\
eturn 1;}\n      return 0;\n    }\nsub file2blast_\
flavor\n      {\n	my $file=shift;\n	if (&file_cont\
ains ($file,\"EBIApplicationResult\",100)){return \
\"EBI\";}\n	elsif (&file_contains ($file,\"NCBI_Bl\
astOutput\",100)){return \"NCBI\";}\n	else {return\
 \"UNKNOWN\";}\n      }\nsub blast_xml2profile \n \
 {\n    my ($name,$seq,$maxid, $minid, $mincov, $f\
ile)=(@_);\n    my (%p, $a, $string, $n);\n    \n\\
n\n    if ($BLAST_TYPE eq \"EBI\" || &file_contain\
s ($file,\"EBIApplicationResult\",100)){%p=ebi_bla\
st_xml2profile(@_);}\n    elsif ($BLAST_TYPE eq \"\
NCBI\" || &file_contains ($file,\"NCBI_BlastOutput\
\",100)){%p=ncbi_blast_xml2profile(@_);}\n    else\
 \n      {\n	myexit(add_error ( $$,$$,getppid(), \\
"BLAST_FAILURE::unkown XML\",$CL));\n      }\n    \
for ($a=0; $a<$p{n}; $a++)\n      {\n	my $name=$p{\
$a}{name};\n	$p{$name}{seq}=$p{$a}{seq};\n	$p{$nam\
e}{index}=$a;\n      }\n    return %p;\n  }\nsub n\
cbi_tblastx_xml2lib_file \n  {\n    my  ($outlib,$\
string)=(@_);\n    my ($L,$l, $a,$b,$c,$d,$i,$nhit\
s,@identifyerL);\n    my (%ITERATION);\n      \n  \
  open (F, \">>$outlib\");\n    \n    $seq=~s/[^a-\
zA-Z]//g;\n    $L=length ($seq);\n    \n    %ITERA\
TION=xml2tag_list ($string, \"Iteration\");\n    f\
or ($i=0; $i<$ITERATION{n};$i++)\n      {\n	my ($q\
index, $qlen, %hit, $string);\n	$string=$ITERATION\
{$i}{body};\n\n	$qindex=xmltag2value($string,\"Ite\
ration_iter-num\");\n	$qlen  =xmltag2value($string\
,\"Iteration_query-len\");\n	%hit=&xml2tag_list  (\
$string, \"Hit\");\n\n	for ($a=0; $a<$hit{n}; $a++\
)\n	  {\n	    my ($string);\n	    $string=$hit{$a}\
{body};\n	 \n	    $hindex=xmltag2value($string,\"H\
it_accession\")+1;\n	    if ($hindex<=$qindex){nex\
t;}\n	    else  {print F  \"# $qindex $hindex\\n\"\
;}\n		   \n	   \n	    $hlen=xmltag2value  ($string\
,\"Hit_len\");\n	    %HSP=&xml2tag_list  ($string,\
 \"Hsp\");\n	   \n	    for ($b=0; $b<$HSP{n}; $b++\
)\n	      {\n		my ($string, $qs,$qe,$qf,$hs,$he,$h\
f,$s, $d, $e);\n		$string=$HSP{$b}{body};\n	\n		$q\
s=xmltag2value  ($string,\"Hsp_query-from\");\n		$\
qe=xmltag2value  ($string,\"Hsp_query-to\");\n		$q\
f=xmltag2value  ($string,\"Hsp_query-frame\");\n\n\
		$hs=xmltag2value  ($string,\"Hsp_hit-from\");\n	\
	$he=xmltag2value  ($string,\"Hsp_hit-to\");\n		$h\
f=xmltag2value  ($string,\"Hsp_hit-frame\");\n		\n\
		$s=xmltag2value  ($string,\"Hsp_identity\");\n		\
$l=xmltag2value  ($string,\"Hsp_align-len\");\n		$\
s=int(($s*100)/$l);\n		\n		if ($qf>0)\n		  {$rqs=$\
qs; $rqe=$qe;}\n		else\n		  {\n		    $rqe=($qlen-$\
qs)+1;\n		    $rqs=($qlen-$qe)+1;\n		  }\n		\n		if\
 ($hf>0)\n		  {$rhs=$hs; $rhe=$he;}\n		else\n		  {\
\n		    $rhe=($hlen-$hs)+1;\n		    $rhs=($hlen-$he\
)+1;\n		  }\n		for ($d=0,$e=$rqs; $e<$rqe; $e++,$d\
++)\n		  {\n		    my ($r1,$r2);\n		    $r1=$e;\n		\
    $r2=$rhs+$d;\n		    print F \" $r1 $r2 $s 0\\n\
\";\n		  }\n	      }\n	  }\n      }\n    print F \\
"! SEQ_1_TO_N\\n\";\n    \n    close (F);\n    ret\
urn %lib;\n  }\n\nsub ncbi_tblastpx_xml2lib_file \\
n  {\n    my  ($outlib,$string,%s)=(@_);\n    my (\
$L,$l, $a,$b,$c,$d,$i,$nhits,@identifyerL);\n    m\
y (%ITERATION,%hdes, %qdes);\n      \n    open (F,\
 \">>$outlib\");\n    \n    $seq=~s/[^a-zA-Z]//g;\\
n    $L=length ($seq);\n    \n    %ITERATION=xml2t\
ag_list ($string, \"Iteration\");\n    for ($i=0; \
$i<$ITERATION{n};$i++)\n      {\n	my ($qindex, $ql\
en, %hit, $string);\n	$string=$ITERATION{$i}{body}\
;\n\n	$qdef=xmltag2value($string,\"Iteration_query\
-def\");\n	%qdes=&tblastpx_name2description($qdef,\
%s);\n	$qlen  =xmltag2value($string,\"Iteration_qu\
ery-len\");\n	%hit=&xml2tag_list  ($string, \"Hit\\
");\n\n	for ($a=0; $a<$hit{n}; $a++)\n	  {\n	    m\
y ($string);\n	    $string=$hit{$a}{body};\n	    $\
hdef=xmltag2value($string,\"Hit_def\");\n	    %hde\
s=&tblastpx_name2description($hdef,%s);\n	    if (\
$hdes{index}<=$qdes{index}){next;}\n	    else  {pr\
int F  \"# $qdes{index} $hdes{index}\\n\";}\n		   \
\n	   \n	    $hlen=xmltag2value  ($string,\"Hit_le\
n\");\n	    %HSP=&xml2tag_list  ($string, \"Hsp\")\
;\n	   \n	    for ($b=0; $b<$HSP{n}; $b++)\n	     \
 {\n		my ($string, $l,$qs,$qe,$qf,$hs,$he,$hf,$s, \
$d, $e, @s1, @s2);\n		$string=$HSP{$b}{body};\n	\n\
		$qs=xmltag2value  ($string,\"Hsp_query-from\");\\
n		$qe=xmltag2value  ($string,\"Hsp_query-to\");\n\
		$qf=$qdes{frame};\n		$qseq=xmltag2value  ($strin\
g,\"Hsp_qseq\");\n		\n		$hs=xmltag2value  ($string\
,\"Hsp_hit-from\");\n		$he=xmltag2value  ($string,\
\"Hsp_hit-to\");\n		$hf=$hdes{frame};\n		$hseq=xml\
tag2value  ($string,\"Hsp_hseq\");\n		\n		$s=xmlta\
g2value  ($string,\"Hsp_identity\");\n		$l=xmltag2\
value  ($string,\"Hsp_align-len\");\n		$s=int(($s*\
100)/$l);\n		@s1=tblastpx_hsp2coordinates($qseq,$q\
s,$qe,%qdes);\n		@s2=tblastpx_hsp2coordinates($hse\
q,$hs,$he,%hdes);\n		\n		\n		for ($f=0; $f<=$#s1; \
$f++)\n		  {\n		    if ($s1[$f]==-1 || $s2[$f]==-1\
){next;}\n		    else \n		      {\n			print F \" $s\
1[$f] $s2[$f] $s 0\\n\";\n		      }\n		  }\n	     \
 }\n	  }\n      }\n    print F \"! SEQ_1_TO_N\\n\"\
;\n    \n    close (F);\n    return %lib;\n  }\nsu\
b tblastpx_hsp2coordinates\n  {\n    my ($seq, $s,\
 $e, %des)=@_;\n    my @list;\n    my @sa;\n    my\
 @gap=(-1,-1,-1);\n    \n    $s=$des{start}+3*($s-\
1);\n  \n    if ($des{strand} eq \"d\"){;}\n    el\
se {$s=($des{length}-$s)+1;}\n    \n    foreach $c\
 (split (//,$seq))\n      {\n	if ( $c eq '-'){push\
 (@list,@gap);}\n	elsif ($des{strand} eq \"d\")\n	\
  {\n	    push(@list,$s++,$s++,$s++);\n	  }\n	else\
\n	  {\n	    push(@list, $s--,$s--,$s--);\n	  }\n \
     }\n    return @list;\n  }\n\nsub tblastpx_nam\
e2description\n  {\n    my ($name, %s)=@_;\n    my\
 @at=split(\"__\", $name);\n    my %des;\n\n    $d\
es{name}=$at[0];\n    $des{strand}=$at[1];\n    \n\
    $des{start}=$at[2];\n    $des{end}=$at[3];\n  \
  $des{length}=$at[4];\n    $des{index}=$s{$at[0]}\
{order}+1;\n    return %des;\n  }  \nsub ncbi_blas\
t_xml2profile \n  {\n    my ($name,$seq,$maxid, $m\
inid, $mincov, $string)=(@_);\n    my ($L,$l, $a,$\
b,$c,$d,$nhits,@identifyerL);\n    \n    \n    $se\
q=~s/[^a-zA-Z]//g;\n    $L=length ($seq);\n   \n  \
  %query=&xml2tag_list ($string, \"Iteration_query\
-def\");\n    $name=$query{0}{body};\n    \n    %h\
it=&xml2tag_list ($string, \"Hit\");\n    \n    \n\
    for ($nhits=0,$a=0; $a<$hit{n}; $a++)\n      {\
\n	my ($ldb,$id, $identity, $expectation, $start, \
$end, $coverage, $r);\n	my (%ID,%DE,%HSP);\n	\n	$l\
db=\"\";\n\n	%ID=&xml2tag_list ($hit{$a}{body}, \"\
Hit_id\");\n	$identifyer=$ID{0}{body};\n	\n	%DE=&x\
ml2tag_list ($hit{$a}{body}, \"Hit_def\");\n	$defi\
nition=$DE{0}{body};\n	\n	%HSP=&xml2tag_list ($hit\
{$a}{body}, \"Hsp\");\n	for ($b=0; $b<$HSP{n}; $b+\
+)\n	  {\n	    my (%START,%END,%E,%I,%Q,%M);\n\n	 \
\n	    %START=&xml2tag_list ($HSP{$b}{body}, \"Hsp\
_query-from\");\n	    %HSTART=&xml2tag_list ($HSP{\
$b}{body}, \"Hsp_hit-from\");\n	    \n	    %LEN=  \
&xml2tag_list ($HSP{$b}{body}, \"Hsp_align-len\");\
\n	    %END=  &xml2tag_list ($HSP{$b}{body}, \"Hsp\
_query-to\");\n	    %HEND=  &xml2tag_list ($HSP{$b\
}{body}, \"Hsp_hit-to\");\n	    %E=&xml2tag_list  \
   ($HSP{$b}{body}, \"Hsp_evalue\");\n	    %I=&xml\
2tag_list     ($HSP{$b}{body}, \"Hsp_identity\");\\
n	    %Q=&xml2tag_list     ($HSP{$b}{body}, \"Hsp_\
qseq\");\n	    %M=&xml2tag_list     ($HSP{$b}{body\
}, \"Hsp_hseq\");\n	    \n	    for ($e=0; $e<$Q{n}\
; $e++)\n\n	      {\n		$qs=$Q{$e}{body};\n		$ms=$M\
{$e}{body};\n		\n		$expectation=$E{$e}{body};\n		$\
identity=($LEN{$e}{body}==0)?0:$I{$e}{body}/$LEN{$\
e}{body}*100;\n		$start=$START{$e}{body};\n		$end=\
$END{$e}{body};\n		$Hstart=$HSTART{$e}{body};\n		$\
Hend=$HEND{$e}{body};\n	\n		$coverage=(($end-$star\
t)*100)/$L;\n	\n		if ($identity>$maxid || $identit\
y<$minid || $coverage<$mincov){next;}\n		@lr1=(spl\
it (//,$qs));\n		@lr2=(split (//,$ms));\n		$l=$#lr\
1+1;\n		for ($c=0;$c<$L;$c++){$p[$nhits][$c]=\"-\"\
;}\n		for ($d=0,$c=0; $c<$l; $c++)\n		  {\n		    $\
r=$lr1[$c];\n		    if ( $r=~/[A-Za-z]/)\n		      {\
\n			\n			$p[$nhits][$d + $start-1]=$lr2[$c];\n			\
$d++;\n		      }\n		  }\n		$Qseq[$nhits]=$qs;\n		$\
Hseq[$nhits]=$ms;\n		$QstartL[$nhits]=$start;\n		$\
HstartL[$nhits]=$Hstart;\n		$identityL[$nhits]=$id\
entity;\n		$endL[$nhits]=$end;\n		$definitionL[$nh\
its]=$definition;\n		$identifyerL[$nhits]=$identif\
yer;\n		$comment[$nhits]=\"$ldb|$identifyer [Eval=\
$expectation][id=$identity%][start=$Hstart end=$He\
nd]\";\n		$nhits++;\n	      }\n	  }\n      }\n    \
\n    \n    $profile{n}=0;\n    $profile{$profile{\
n}}{name}=$name;\n    $profile{$profile{n}}{seq}=$\
seq;\n    $profile {n}++;\n    \n    for ($a=0; $a\
<$nhits; $a++)\n      {\n	$n=$a+1;\n	\n	$profile{$\
n}{name}=\"$name\\_$a\";\n	$profile{$n}{seq}=\"\";\
\n	$profile{$n}{Qseq}=$Qseq[$a];\n	$profile{$n}{Hs\
eq}=$Hseq[$a];\n	$profile{$n}{Qstart}=$QstartL[$a]\
;\n	$profile{$n}{Hstart}=$HstartL[$a];\n	$profile{\
$n}{identity}=$identityL[$a];\n	$profile{$n}{defin\
ition}=$definitionL[$a];\n	$profile{$n}{identifyer\
}=$identifyerL[$a];\n	$profile{$n}{comment}=$comme\
nt[$a];\n\n	for ($b=0; $b<$L; $b++)\n	  {\n	    if\
 ($p[$a][$b])\n	      {\n		$profile{$n}{seq}.=$p[$\
a][$b];\n	      }\n	    else\n	      {\n		$profile\
{$n}{seq}.=\"-\";\n	      }\n	  }\n      }\n    \n\
    $profile{n}=$nhits+1;\n    return %profile;\n \
 }\nsub ebi_blast_xml2profile \n  {\n    my ($name\
,$seq,$maxid, $minid, $mincov, $string)=(@_);\n   \
 my ($L,$l, $a,$b,$c,$d,$nhits,@identifyerL,$ident\
ifyer);\n    \n\n    \n    $seq=~s/[^a-zA-Z]//g;\n\
    $L=length ($seq);\n    %hit=&xml2tag_list ($st\
ring, \"hit\");\n    \n    for ($nhits=0,$a=0; $a<\
$hit{n}; $a++)\n      {\n	my ($ldb,$id, $identity,\
 $expectation, $start, $end, $coverage, $r);\n	my \
(%Q,%M,%E,%I);\n	\n	$ldb=&tag2value ($hit{$a}{open\
}, \"database\");\n	$identifyer=&tag2value ($hit{$\
a}{open}, \"id\");\n\n	$description=&tag2value ($h\
it{$a}{open}, \"description\");\n	\n	%Q=&xml2tag_l\
ist ($hit{$a}{body}, \"querySeq\");\n	%M=&xml2tag_\
list ($hit{$a}{body}, \"matchSeq\");\n	%E=&xml2tag\
_list ($hit{$a}{body}, \"expectation\");\n	%I=&xml\
2tag_list ($hit{$a}{body}, \"identity\");\n	\n\n	f\
or ($b=0; $b<$Q{n}; $b++)\n	  {\n\n	    $qs=$Q{$b}\
{body};\n	    $ms=$M{$b}{body};\n	    \n	    $expe\
ctation=$E{$b}{body};\n	    $identity=$I{$b}{body}\
;\n	    \n	    	    \n	    $start=&tag2value ($Q{$\
b}{open}, \"start\");\n	    $end=&tag2value ($Q{$b\
}{open}, \"end\");\n	    $startM=&tag2value ($M{$b\
}{open}, \"start\");\n	    $endM=&tag2value ($M{$b\
}{open}, \"end\");\n	    $coverage=(($end-$start)*\
100)/$L;\n	    \n	   # print \"$id: ID: $identity \
COV: $coverage [$start $end]\\n\";\n	    \n	    if\
 ($identity>$maxid || $identity<$minid || $coverag\
e<$mincov){next;}\n	    # print \"KEEP\\n\";\n\n	 \
   \n	    @lr1=(split (//,$qs));\n	    @lr2=(split\
 (//,$ms));\n	    $l=$#lr1+1;\n	    for ($c=0;$c<$\
L;$c++){$p[$nhits][$c]=\"-\";}\n	    for ($d=0,$c=\
0; $c<$l; $c++)\n	      {\n		$r=$lr1[$c];\n		if ( \
$r=~/[A-Za-z]/)\n		  {\n		    \n		    $p[$nhits][$\
d + $start-1]=$lr2[$c];\n		    $d++;\n		  }\n	    \
  }\n	  \n	    $Qseq[$nhits]=$qs;\n	    $Hseq[$nhi\
ts]=$ms;\n	    $QstartL[$nhits]=$start;\n	    $Hst\
artL[$nhits]=$Hstart;\n	    $identityL[$nhits]=$id\
entity;\n	    $endL[$nhits]=$end;\n	    $definitio\
nL[$nhits]=$definition;\n	    $identifyerL[$nhits]\
=$identifyer;\n	    $comment[$nhits]=\"$ldb|$ident\
ifyer [Eval=$expectation][id=$identity%][start=$st\
artM end=$endM]\";\n	    $nhits++;\n	  }\n      }\\
n    \n    $profile{n}=0;\n    $profile{$profile{n\
}}{name}=$name;\n    $profile{$profile{n}}{seq}=$s\
eq;\n    $profile {n}++;\n    \n    for ($a=0; $a<\
$nhits; $a++)\n      {\n	$n=$a+1;\n	$profile{$n}{n\
ame}=\"$name\\_$a\";\n	$profile{$n}{seq}=\"\";\n	$\
profile{$n}{Qseq}=$Qseq[$a];\n	$profile{$n}{Hseq}=\
$Hseq[$a];\n	$profile{$n}{Qstart}=$QstartL[$a];\n	\
$profile{$n}{Hstart}=$HstartL[$a];\n	$profile{$n}{\
identity}=$identityL[$a];\n	$profile{$n}{definitio\
n}=$definitionL[$a];	\n	$profile{$n}{identifyer}=$\
identifyerL[$a];\n	$profile{$n}{comment}=$comment[\
$a];\n\n	for ($b=0; $b<$L; $b++)\n	  {\n	    if ($\
p[$a][$b])\n	      {\n		$profile{$n}{seq}.=$p[$a][\
$b];\n	      }\n	    else\n	      {\n		$profile{$n\
}{seq}.=\"-\";\n	      }\n	  }\n      }\n    $prof\
ile{n}=$nhits+1;\n    \n    return %profile;\n  }\\
nsub output_profile\n  {\n    my ($outfile,$profil\
eR, $trim)=(@_);\n    my ($a);\n    my %profile=%$\
profileR;\n    my $P= new FileHandle;\n    my $tmp\
=vtmpnam();\n    \n    open ($P, \">$tmp\");\n    \
for ($a=0; $a<$profile{n}; $a++)\n      {\n	print \
$P \">$profile{$a}{name} $profile{$a}{comment}\\n$\
profile{$a}{seq}\\n\";\n      }\n    close ($P);\n\
\n    if ( $trim)\n      {\n	&safe_system (\"t_cof\
fee -other_pg seq_reformat -in $tmp -action +trim \
_aln_%%$trim\\_K1 -output fasta_aln -out $outfile\\
");\n      }\n    else\n      {\n	&safe_system (\"\
mv $tmp $outfile\");\n      }\n    return;\n  }\ns\
ub blast_xml2hit_list\n  {\n    my $string=(@_[0])\
;\n    return &xml2tag_list ($string, \"hit\");\n \
 }\nsub xmltag2value\n  {\n    my ($string_in, $ta\
g)=@_;\n    my %TAG;\n    %TAG=xml2tag_list ($stri\
ng_in, $tag);\n    return $TAG{0}{body};\n  }\n   \
   \nsub xml2tag_list  \n  {\n    my ($string_in,$\
tag)=@_;\n    my $tag_in, $tag_out;\n    my %tag;\\
n    \n    if (-e $string_in)\n      {\n	$string=&\
file2string ($string_in);\n      }\n    else\n    \
  {\n	$string=$string_in;\n      }\n    $tag_in1=\\
"<$tag \";\n    $tag_in2=\"<$tag>\";\n    $tag_out\
=\"/$tag>\";\n    $string=~s/>/>##1/g;\n    $strin\
g=~s/</##2</g;\n    $string=~s/##1/<#/g;\n    $str\
ing=~s/##2/#>/g;\n    @l=($string=~/(\\<[^>]+\\>)/\
g);\n    $tag{n}=0;\n    $in=0;$n=-1;\n  \n \n\n  \
  foreach $t (@l)\n      {\n\n	$t=~s/<#//;\n	$t=~s\
/#>//;\n	\n	if ( $t=~/$tag_in1/ || $t=~/$tag_in2/)\
\n	  {\n	 \n	    $in=1;\n	    $tag{$tag{n}}{open}=\
$t;\n	    $n++;\n	    \n	  }\n	elsif ($t=~/$tag_ou\
t/)\n	  {\n	    \n\n	    $tag{$tag{n}}{close}=$t;\\
n	    $tag{n}++;\n	    $in=0;\n	  }\n	elsif ($in)\\
n	  {\n	   \n	    $tag{$tag{n}}{body}.=$t;\n	  }\n\
      }\n  \n    return %tag;\n  }\n\n\nsub seq2go\
r_prediction \n  {\n    my ($name, $seq,$infile, $\
outfile, $gor_seq, $gor_obs)=(@_);\n    my ($l);\n\
    \n    `gorIV -prd $infile -seq $gor_seq -obs $\
gor_obs > gor_tmp`;\n    open (GR, \">$outfile\");\
\n    open (OG, \"gor_tmp\");\n\n    while (<OG>)\\
n      {\n	\n	$l=$_;\n	if ($l=~/\\>/){print GR \"$\
l\";}\n	elsif ( $l=~/Predicted Sec. Struct./)\n	  \
{\n	    $l=~s/Predicted Sec. Struct\\.//;\n	    pr\
int GR \"$l\";\n	  }\n      }\n    close (GR);\n  \
  close (OG);\n    return;\n  }\nsub seq2msa_tm_pr\
ediction \n  {\n    my ($name, $seq, $db, $infile,\
 $outfile, $arch, $psv)=(@_);\n    my (%p,%gseq,%R\
, $blast_output, %s, $l);\n    my $R2=new FileHand\
le;\n    my $db=\"uniprot\";\n    my $method=\"psi\
tm\";\n    my $SERVER=\"EBI\";\n    \n    $blast_o\
utput=&run_blast ($name,\"blastp\", $db, $infile, \
\"outfile\");\n    \n    if (&cache_file(\"GET\",$\
infile,$name,$method,$db,$outfile,$SERVER))\n     \
 {\n	print \"\\tPSITM: USE Cache\\n\";\n	return $o\
utfile;\n      }\n    else\n      {\n	$CACHE_STATU\
S=\"COMPUTE CACHE\";\n	%p=blast_xml2profile($name,\
$seq,$maxid, $minid,$mincov,$blast_output);\n	\n	\\
n	open (F, \">tm_input\");\n	for (my $a=0; $a<$p{n\
}; $a++)\n	  {\n	    my $s;\n	    \n	    $s=$p{$a}\
{seq};\n	    $s=uc($s);\n	    print F \">$p{$a}{na\
me}\\n$s\\n\";\n	    #print stdout \">$p{$a}{name}\
\\n$s\\n\";\n	  }\n	close (F);\n	print \"\\tPSITM:\
 kept  $p{n} Homologues for Sequence $p{0}{name}\\\
n\";\n	&safe_system (\"t_coffee -other_pg fasta_se\
q2hmmtop_fasta.pl -in=tm_input -out=$outfile -outp\
ut=cons -cov=70 -trim=95 -arch=$arch -psv=$psv\");\
\n	unlink (\"tm_input\");\n	&cache_file(\"SET\",$i\
nfile,$name,$method,$db,$outfile,$SERVER);\n	retur\
n;\n      }\n  }\n\n\nsub seq2msa_gor_prediction \\
n  {\n    my ($name, $seq,$infile, $outfile, $gor_\
seq, $gor_obs)=(@_);\n    my (%p,%gseq,%R, $blast_\
output, %s, $l);\n    my $R2=new FileHandle;\n    \
my $db=\"uniprot\";\n    my $method=\"psigor\";\n \
   my $SERVER=\"EBI\";\n    \n    $blast_output=&r\
un_blast ($name,\"blastp\", \"uniprot\", $infile, \
\"outfile\");\n    \n    if (&cache_file(\"GET\",$\
infile,$name,$method,$db,$outfile,$SERVER))\n     \
 {\n	print \"\\tPSIGOR: USE Cache\\n\";\n	return $\
outfile;\n      }\n    else\n      {\n	$CACHE_STAT\
US=\"COMPUTE CACHE\";\n	%p=blast_xml2profile($name\
,$seq,$maxid, $minid,$mincov,$blast_output);\n	\n	\
\n	open (F, \">gor_input\");\n	for (my $a=0; $a<$p\
{n}; $a++)\n	  {\n	    my $s;\n	    \n	    $s=$p{$\
a}{seq};\n	    $s=uc($s);\n	    print F \">$p{$a}{\
name}\\n$s\\n\";\n	    #print stdout \">$p{$a}{nam\
e}\\n$s\\n\";\n	  }\n	close (F);\n	print \"\\tGORT\
M: kept  $p{n} Homologues for Sequence $p{0}{name}\
\\n\";\n	&safe_system (\"t_coffee -other_pg fasta_\
seq2hmmtop_fasta.pl -in=gor_input -out=$outfile -o\
utput=cons -cov=70 -trim=95 -gor_seq=$gor_seq -gor\
_obs=$gor_obs -mode=gor\");\n	unlink (\"tm_input\"\
);\n	&cache_file(\"SET\",$infile,$name,$method,$db\
,$outfile,$SERVER);\n	return;\n      }\n  }\n\n\n\\
nsub run_blast\n  {\n    my ($name, $method, $db, \
$infile, $outfile, $run)=(@_);\n    if (!$run){$ru\
n=1;}\n    \n    \n    if (&cache_file(\"GET\",$in\
file,$name,$method,$db,$outfile,$SERVER) && is_val\
id_blast_xml ($outfile))\n      {return $outfile;}\
\n    else\n      {\n	$CACHE_STATUS=\"COMPUTE CACH\
E\";\n	if ( $SERVER eq \"EBI_SOAP\")\n	  {\n	    &\
check_configuration (\"EMAIL\",\"SOAP::Light\",\"I\
NTERNET\");\n	    \n	    $cl_method=$method;\n	   \
 if ($cl_method =~/wu/)\n	      {\n		$cl_method=~s\
/wu//;\n		if ( $cl_method eq \"psiblast\")\n		  {\\
n		    add_warning($$,$$,\"PSI BLAST cannot be use\
d with the wuBLAST Client. Use server=EBI Or serve\
r=LOCAL. blastp will be used instead\");\n		    $c\
l_method=\"blastp\";\n		  }\n		\n		$command=\"t_co\
ffee -other_pg wublast.pl --email $EMAIL $infile -\
D $db -p $cl_method --outfile $outfile -o xml>/dev\
/null 2>/dev/null\";\n		&safe_system ( $command);\\
n		if (-e \"$outfile.xml\") {`mv $outfile.xml $out\
file`;}\n	      }\n	    else\n	      {\n		if ($cl_\
method eq \"psiblast\"){$cl_method =\"blastp -j5\"\
;}\n		\n		$command=\"t_coffee -other_pg blastpgp.p\
l --email $EMAIL $infile -d $db --outfile $outfile\
 -p $cl_method --mode PSI-Blast>/dev/null 2>/dev/n\
ull\";\n		&safe_system ( $command);\n		\n		if (-e \
\"$outfile.xml\") {`mv $outfile.xml $outfile`;}\n	\
      }\n	  }\n	elsif ($SERVER eq \"EBI_REST\" || \
$SERVER eq \"EBI\")\n	  {\n	   \n	    $cl_method=$\
method;\n	    &check_configuration(\"EMAIL\",\"XML\
::Simple\", \"INTERNET\");\n	    if ($db eq \"unip\
rot\"){$db1=\"uniprotkb\";}\n	    else {$db1=$db;}\
\n	    \n\n	    if ($cl_method =~/wu/)\n	      {\n\
		$cl_method=~s/wu//;\n		if ( $cl_method eq \"psib\
last\"){$cl_method=\"blastp\";}\n		\n		$command=\"\
t_coffee -other_pg wublast_lwp.pl --email $EMAIL -\
D $db1 -p $cl_method --outfile $outfile --align 7 \
--stype protein $infile>/dev/null 2>/dev/null\";\n\
		\n	      }\n	    else\n	      {\n		if ( $cl_meth\
od =~/psiblast/){$cl_method =\"blastp -j5\";}\n		$\
command=\"t_coffee -other_pg ncbiblast_lwp.pl --em\
ail $EMAIL -D $db1 -p $cl_method --outfile $outfil\
e --align 7 --stype protein $infile>/dev/null 2>/d\
ev/null\";\n	      }\n	    &safe_system ( $command\
,5);\n	    if (-e \"$outfile.out.xml\") {`mv $outf\
ile.out.xml $outfile`;}\n	    elsif (-e \"$outfile\
.xml.xml\"){`mv $outfile.xml.xml $outfile`;}\n	   \
 elsif (-e \"$outfile.out..xml\") {`mv $outfile.ou\
t..xml $outfile`;}\n	    elsif (-e \"$outfile.xml.\
.xml\"){`mv $outfile.xml..xml $outfile`;}\n	  }\n	\
elsif ($SERVER eq \"NCBI\")\n	  {\n	    &check_con\
figuration (\"blastcl3\",\"INTERNET\");\n	    if (\
$db eq \"uniprot\"){$cl_db=\"nr\";}\n	    else {$c\
l_db=$db;}\n	    \n	    if ( $method eq \"psiblast\
\")\n	      {\n		add_warning($$,$$,\"PSI BLAST can\
not be used with the NCBI BLAST Client. Use server\
=EBI Or server=LOCAL. blastp will be used instead\\
");\n		$cl_method=\"blastp\";\n	      }\n	    else\
\n	      {\n		$cl_method=$method;\n	      }\n	    \
$command=\"blastcl3 -p $cl_method -d $cl_db -i $in\
file -o $outfile -m 7\";\n	    &safe_system ($comm\
and);\n	  }\n	elsif ($SERVER =~/CLIENT_(.*)/)\n	  \
{\n	    my $client=$1;\n	    $command=\"$client -p\
 $method -d $db -i $infile -o $outfile -m 7\";\n	 \
   &safe_system ($command);\n	  }\n	elsif ( $SERVE\
R eq \"LOCAL_blastall\")\n	  {\n	    &check_config\
uration (\"blastall\");\n	    if ($method eq \"bla\
stp\")\n	      {\n		$command=\"blastall -d $db -i \
$infile -o $outfile -m7 -p blastp\";\n	      }\n	 \
   &safe_system ($command);\n	  }\n	elsif ( $SERVE\
R eq \"LOCAL\")\n	  {\n	    \n	    if ($ENV{\"BLAS\
T_DB_DIR\"})\n	      {\n		$x=$ENV{\"BLAST_DB_DIR\"\
};\n		$cl_db=\"$x$db\";\n	      }\n	    else\n	   \
   {\n		$cl_db=$db;\n	      }\n	    \n	    if ($me\
thod eq \"blastp\")\n	      {\n		&check_configurat\
ion(\"blastpgp\");\n		$command=\"blastpgp -d $cl_d\
b -i $infile -o $outfile -m7 -j1\";\n	      }\n	  \
  elsif ($method eq \"psiblast\")\n	      {\n		&ch\
eck_configuration(\"blastpgp\");\n		$command=\"bla\
stpgp -d $cl_db -i $infile -o $outfile -m7 -j5\";\\
n	      }\n	    elsif ($method eq \"blastn\")\n	  \
    {\n		&check_configuration(\"blastall\");\n		$c\
ommand=\"blastall -p blastn -d $cl_db -i $infile -\
o $outfile -m7 -W6\";\n	      }	\n	    &safe_syste\
m ($command);\n	  }\n	else\n	  {\n	    \n	    myex\
it(add_error (EXIT_FAILURE,$$,$$,getppid(), \"BLAS\
T_FAILURE::UnknownServer\",$CL));\n	  }\n	\n	if ( \
!-e $outfile || !&is_valid_blast_xml($outfile))\n	\
  {\n	    \n	    if ( -e $outfile)\n	      {\n		ad\
d_warning ($$,$$,\"Corrupted Blast Output (Run $ru\
n)\");\n		unlink($outfile);\n	      }\n	    \n	   \
 if ( $run==$BLAST_MAX_NRUNS)\n	      {\n	\n		myex\
it(add_error (EXIT_FAILURE,$$,$$,getppid(), \"BLAS\
T_FAILURE::UnknownReason\", \"$command\"));\n	    \
  }\n	    else\n	      {\n		my $out;\n		if ($SERVE\
R eq \"NCBI\") {$SERVER=\"EBI\"; }\n		elsif ($SERV\
ER eq \"EBI\"){$SERVER=\"NCBI\";}\n		add_warning (\
$$,$$,\"Blast for $name failed (Run: $run out of $\
BLAST_MAX_NRUNS. Use $SERVER)\");\n		$out=&run_bla\
st ($name, $method, $db,$infile, $outfile, $run+1)\
;\n		if ($SERVER eq \"NCBI\") {$SERVER=\"EBI\"; }\\
n		elsif ($SERVER eq \"EBI\"){$SERVER=\"NCBI\";}\n\
		return $out;\n	      }\n	  }\n	\n	&cache_file(\"\
SET\",$infile,$name,$method,$db,$outfile,$SERVER);\
\n	#system (\"cp $outfile ~/Dropbox/tmp/cedric.out\
\");\n	#die;\n	return $outfile;\n      }\n  }\n\ns\
ub cache_file\n  {\n    my ($cache_mode,$infile,$n\
ame,$method,$db, $outfile,$server)=(@_);\n    my $\
cache_file;\n    #Protect names so that they can b\
e turned into legal filenames\n    $name=&clean_fi\
le_name ($name);\n\n    if ($db=~/\\//)\n      {\n\
	$db=~/([^\\/]+)$/;\n	$db=$1;\n      }\n    $cache\
_file_sh=\"$name.$method.$db.$server.tmp\";\n    $\
cache_file=\"$CACHE/$name.$method.$db.$server.tmp\\
";\n    \n    if ($infile ne \"\")\n      {\n	$cac\
he_file_infile_sh=\"$name.$method.$db.$server.infi\
le.tmp\";\n	$cache_file_infile=\"$CACHE/$name.$met\
hod.$db.$server.infile.tmp\";\n      }\n    \n    \
if ($cache_mode eq \"GET\")\n      {\n	if ($CACHE \
eq \"\" || $CACHE eq \"no\" || $CACHE eq \"ignore\\
"  || $CACHE eq \"local\" || $CACHE eq \"update\")\
{return 0;}\n	elsif ( !-d $CACHE)\n	  {\n	    prin\
t STDERR \"ERROR: Cache Dir: $CACHE Does not Exist\
\";\n	    return 0;\n	  }\n	else\n	  {\n	    if ( \
-e $cache_file && &fasta_file1_eq_fasta_file2($inf\
ile,$cache_file_infile)==1)\n	      {\n		`cp $cach\
e_file $outfile`;\n		$CACHE_STATUS=\"READ CACHE\";\
\n		return 1;\n	      }\n	  }\n      }\n    elsif \
($cache_mode eq \"SET\")\n      {\n	if ($CACHE eq \
\"\" || $CACHE eq \"no\" || $CACHE eq \"ignore\"  \
|| $CACHE eq \"local\" || $CACHE eq \"update\"){re\
turn 0;}\n	elsif ( !-d $CACHE)\n	  {\n	    print S\
TDERR \"ERROR: Cache Dir: $CACHE Does not Exist\";\
\n	    return 0;\n	  }\n	elsif (-e $outfile)\n	  {\
\n	    `cp $outfile $cache_file`;\n	    if ($cache\
_file_infile ne \"\"){ `cp $infile $cache_file_inf\
ile`;}\n\n	    #functions for updating the cache\n\
	    #`t_coffee -other_pg clean_cache.pl -file $ca\
che_file_sh -dir $CACHE`;\n	    #`t_coffee -other_\
pg clean_cache.pl -file $cache_file_infile_sh -dir\
 $CACHE`;\n	    return 1;\n	  }\n      }\n    $CAC\
HE_STATUS=\"COMPUTE CACHE\";\n    return 0;\n  }\n\
sub file1_eq_file2\n  {\n    my ($f1, $f2)=@_;\n  \
  if ( $f1 eq \"\"){return 1;}\n    elsif ( $f2 eq\
 \"\"){return 1;}\n    elsif ( !-e $f1){return 0;}\
\n    elsif ( !-e $f2){return 0;}\n    elsif ($f1 \
eq \"\" || $f2 eq \"\" || `diff $f1 $f2` eq \"\"){\
return 1;}\n    \n    return 0;\n  }\nsub clean_fi\
le_name \n  {\n    my $name=@_[0];\n    \n    $nam\
e=~s/[^A-Za-z1-9.-]/_/g;\n    return $name;\n  }\n\
sub url2file\n  {\n    my ($address, $out)=(@_);\n\
    \n    if (&pg_is_installed (\"wget\"))\n	{\n	 \
 return &safe_system (\"wget $address -O$out >/dev\
/null 2>/dev/null\");\n	}\n    elsif (&pg_is_insta\
lled (\"curl\"))\n      {\n	return &safe_system (\\
"curl $address -o$out >/dev/null 2>/dev/null\");\n\
      }\n    else\n      {\n	myexit(flus_error(\"n\
either curl nor wget are installed. Imnpossible to\
 fectch remote file\"));\n	exit ($EXIT_FAILURE);\n\
      }\n  }\nsub fasta_file1_eq_fasta_file2\n  {\\
n    my ($f1, $f2)=@_;\n    my (%s1, %s2);\n    my\
 @names;\n    %s1=read_fasta_seq ($f1);\n    %s2=r\
ead_fasta_seq ($f2);\n\n    @names=(keys (%s1));\n\
    \n    foreach $n (keys(%s1))\n      {\n	if ($s\
1{$n}{seq} ne $s2{$n}{seq}){return 0;}\n      } \n\
    \n    foreach $n (keys(%s2))\n      {\n	if ($s\
1{$n}{seq} ne $s2{$n}{seq}){return 0;}\n      }\n \
   return 1;\n  }\n	\n\n\nsub read_template_file\n\
{\n	my $pdb_templates = @_[0];\n	open (TEMP, \"<$p\
db_templates\");\n	my %temp_h;\n	while (<TEMP>)\n{\
\n		$line = $_;\n 		$line =~/(\\S+)\\s(\\S+)/;\n 	\
	$temp_h{$1}= $2;\n}\n	close(TEMP);\n	return %temp\
_h;\n}\n\nsub calc_rna_template\n{\n	my ($mode, $i\
nfile, $pdbfile, $outfile)=@_;\n	my %s, %h ;\n	my \
$result;\n	my (@profiles);\n	&set_temporary_dir (\\
"set\",$infile,\"seq.pep\");\n	%s=read_fasta_seq (\
\"seq.pep\");\n	\n	%pdb_template_h = &read_templat\
e_file($pdbfile);\n	my $pdb_chain;\n	open (R, \">r\
esult.aln\");\n\n\n	#print stdout \"\\n\";\n	forea\
ch $seq (keys(%s))\n	{\n		if ($pdb_template_h{$seq\
} eq \"\")\n		{\n			next;\n		}\n		open (F, \">seqf\
ile\");\n		print (F \">$s{$seq}{name}\\n$s{$seq}{s\
eq}\\n\");\n		close (F);\n		$pdb_chain = $pdb_temp\
late_h{$seq};\n		$lib_name=\"$s{$seq}{name}.rfold\\
";\n		$lib_name=&clean_file_name ($lib_name);\n		\\
n 		safe_system (\"secondary_struc.py seqfile $CAC\
HE$pdb_chain  $lib_name\");\n		\n		if ( !-e $lib_n\
ame)\n		{\n		myexit(flush_error(\"RNAplfold failed\
 to compute the secondary structure of $s{$seq}{na\
me}\"));\n			myexit ($EXIT_FAILURE);\n		}\n		else\\
n		{\n			print stdout \"\\tProcess: >$s{$seq}{name\
} _F_ $lib_name\\n\";\n			print R \">$s{$seq}{name\
} _F_ $lib_name\\n\";\n		}\n		unshift (@profiles, \
$lib_name);\n	}\n	close (R);\n	&set_temporary_dir \
(\"unset\",$mode, $method,\"result.aln\",$outfile,\
 @profiles);\n}\n\n\n\nsub seq2rna_pair{\n	my ($mo\
de, $pdbfile1, $pdbfile2, $method, $param, $outfil\
e)=@_;\n	\n	if ($method eq \"runsara.py\")\n	{\n		\
open(TMP,\"<$pdbfile1\");\n		my $count = 0;\n		my \
$line;\n		while (<TMP>)\n		{\n			$line = $_;\n			i\
f ($count ==1)\n			{\n				last;\n			}\n			$count +\
= 1;\n		}\n\n		\n		$chain1 = substr($line,length($\
line)-3,1);\n\n		close TMP;\n		open(TMP,\"<$pdbfil\
e2\");\n		my $count = 0;\n		while (<TMP>)\n		{\n		\
	$line = $_;\n			if ($count ==1)\n			{\n				last;\\
n			}\n			$count += 1;\n		}\n		$chain2 = substr($l\
ine,length($line)-3,1);\n		close TMP;\n\n		$tmp_fi\
le=&vtmpnam();\n	\n		safe_system(\"runsara.py $pdb\
file1 $chain1 $pdbfile2 $chain2 -s -o $tmp_file --\
limitation 5000 > /dev/null 2> /dev/null\") == 0 o\
r die \"sara did not work $!\\n\";\n		open(TMP,\"<\
$tmp_file\") or die \"cannot open the sara tmp fil\
e:$!\\n\";\n		open(OUT,\">$outfile\") or die \"can\
not open the $outfile file:$!\\n\";\n\n		my $switc\
h = 0;\n		my $seqNum = 0;\n		foreach my $line (<TM\
P>)\n		{\n			next unless ($line=~/SARAALI/);\n			i\
f ($line=~/>/)\n			{\n				$switch =0;\n				print O\
UT \">seq$seqNum\\n\";\n				$seqNum++;				\n			}\n\
			if ($switch < 2){\n				$switch++;\n				next;\n	\
		}\n	\n			if ($line =~/REMARK\\s+SARAALI\\s+([^\\\
*]+)\\*/)\n			{\n				my $string = $1;\n				print O\
UT \"$string\\n\";\n			}\n		}\n		close TMP; \n		cl\
ose OUT;\n		unlink($tmp_file);\n	}\n}\n\nsub seq2t\
blastx_lib\n  {\n    my ($mode, $infile, $outfile)\
=@_;\n    my (%s, $method,$nseq);\n\n    $method=$\
mode;\n    &set_temporary_dir (\"set\",$infile,\"i\
nfile\");\n    %s=read_fasta_seq(\"infile\");\n   \
 \n    \n    foreach $seq (keys(%s))\n      {\n	$s\
list[$s{$seq}{order}]=$s{$seq}{seq};\n	$sname[$s{$\
seq}{order}]=$s{$seq}{name};\n	$slen[$s{$seq}{orde\
r}]=length ($s{$seq}{seq});\n      }\n    $nseq=$#\
sname+1;\n    open (F, \">outfile\");\n    print F\
 \"! TC_LIB_FORMAT_01\\n\";\n    print F \"$nseq\\\
n\";\n    for ($a=0; $a<$nseq;$a++)\n      {\n	pri\
nt F \"$sname[$a] $slen[$a]  $slist[$a]\\n\"\n    \
  }\n    close (F);\n    &safe_system (\"formatdb \
-i infile -p F\");\n    &safe_system (\"blastall -\
p tblastx -i infile -d infile -m 7 -S1>blast.outpu\
t\");\n    \n    ncbi_tblastx_xml2lib_file (\"outf\
ile\", file2string (\"blast.output\"));\n    &set_\
temporary_dir (\"unset\",$mode, $method, \"outfile\
\",$outfile);\n    myexit ($EXIT_SUCCESS);\n    }\\
nsub seq2tblastpx_lib\n  {\n    my ($mode, $infile\
, $outfile)=@_;\n    my (%s, $method,$nseq);\n    \
$method=$mode;\n    &set_temporary_dir (\"set\",$i\
nfile,\"infile\");\n    %s=read_fasta_seq(\"infile\
\");\n    \n    foreach $seq (keys(%s))\n      {\n\
	$slist[$s{$seq}{order}]=$s{$seq}{seq};\n	$sname[$\
s{$seq}{order}]=$s{$seq}{name};\n	$slen[$s{$seq}{o\
rder}]=length ($s{$seq}{seq});\n      }\n    $nseq\
=$#sname+1;\n    open (F, \">outfile\");\n    prin\
t F \"! TC_LIB_FORMAT_01\\n\";\n    print F \"$nse\
q\\n\";\n    for ($a=0; $a<$nseq;$a++)\n      {\n	\
print F \"$sname[$a] $slen[$a]  $slist[$a]\\n\"\n \
     }\n    close (F);\n    &safe_system(\"t_coffe\
e -other_pg seq_reformat -in infile -output tblast\
x_db1 > tblastxdb\");\n    &safe_system (\"formatd\
b -i tblastxdb -p T\");\n    &safe_system (\"blast\
all -p blastp -i tblastxdb -d tblastxdb -m7 >blast\
.output\");\n    ncbi_tblastpx_xml2lib_file (\"out\
file\", file2string (\"blast.output\"), %s);\n    \
&set_temporary_dir (\"unset\",$mode, $method, \"ou\
tfile\",$outfile);\n    myexit ($EXIT_SUCCESS);\n \
   }\n\n\n    \n\n\n\nsub file2head\n      {\n	my \
$file = shift;\n	my $size = shift;\n	my $f= new Fi\
leHandle;\n	my $line;\n	open ($f,$file);\n	read ($\
f,$line, $size);\n	close ($f);\n	return $line;\n  \
    }\nsub file2tail\n      {\n	my $file = shift;\\
n	my $size = shift;\n	my $f= new FileHandle;\n	my \
$line;\n	\n	open ($f,$file);\n	seek ($f,$size*-1, \
2);\n	read ($f,$line, $size);\n	close ($f);\n	retu\
rn $line;\n      }\n\n\nsub vtmpnam\n      {\n	my \
$r=rand(100000);\n	my $f=\"file.$r.$$\";\n	while (\
-e $f)\n	  {\n	    $f=vtmpnam();\n	  }\n	push (@TM\
PFILE_LIST, $f);\n	return $f;\n      }\n\nsub myex\
it\n  {\n    my $code=@_[0];\n    if ($CLEAN_EXIT_\
STARTED==1){return;}\n    else {$CLEAN_EXIT_STARTE\
D=1;}\n    ### ONLY BARE EXIT\n    exit ($code);\n\
  }\nsub set_error_lock\n    {\n      my $name = s\
hift;\n      my $pid=$$;\n\n      \n      &lock4tc\
 ($$,\"LERROR\", \"LSET\", \"$$ -- ERROR: $name $P\
ROGRAM\\n\");\n      return;\n    }\nsub set_lock\\
n  {\n    my $pid=shift;\n    my $msg= shift;\n   \
 my $p=getppid();\n    &lock4tc ($pid,\"LLOCK\",\"\
LRESET\",\"$p$msg\\n\");\n  }\nsub unset_lock\n   \
{\n     \n    my $pid=shift;\n    &lock4tc ($pid,\\
"LLOCK\",\"LRELEASE\",\"\");\n  }\nsub shift_lock\\
n  {\n    my $from=shift;\n    my $to=shift;\n    \
my $from_type=shift;\n    my $to_type=shift;\n    \
my $action=shift;\n    my $msg;\n    \n    if (!&l\
ock4tc($from, $from_type, \"LCHECK\", \"\")){retur\
n 0;}\n    $msg=&lock4tc ($from, $from_type, \"LRE\
AD\", \"\");\n    &lock4tc ($from, $from_type,\"LR\
ELEASE\", $msg);\n    &lock4tc ($to, $to_type, $ac\
tion, $msg);\n    return;\n  }\nsub isshellpid\n  \
{\n    my $p=shift;\n    if (!lock4tc ($p, \"LLOCK\
\", \"LCHECK\")){return 0;}\n    else\n      {\n	m\
y $c=lock4tc($p, \"LLOCK\", \"LREAD\");\n	if ( $c=\
~/-SHELL-/){return 1;}\n      }\n    return 0;\n  \
}\nsub isrootpid\n  {\n    if(lock4tc (getppid(), \
\"LLOCK\", \"LCHECK\")){return 0;}\n    else {retu\
rn 1;}\n  }\nsub lock4tc\n	{\n	  my ($pid,$type,$a\
ction,$value)=@_;\n	  my $fname;\n	  my $host=host\
name;\n	  \n	  if ($type eq \"LLOCK\"){$fname=\"$L\
OCKDIR/.$pid.$host.lock4tcoffee\";}\n	  elsif ( $t\
ype eq \"LERROR\"){ $fname=\"$LOCKDIR/.$pid.$host.\
error4tcoffee\";}\n	  elsif ( $type eq \"LWARNING\\
"){ $fname=\"$LOCKDIR/.$pid.$host.warning4tcoffee\\
";}\n	  \n	  if ($debug_lock)\n	    {\n	      prin\
t STDERR \"\\n\\t---lock4tc(tcg): $action => $fnam\
e =>$value (RD: $LOCKDIR)\\n\";\n	    }\n\n	  if  \
  ($action eq \"LCHECK\") {return -e $fname;}\n	  \
elsif ($action eq \"LREAD\"){return file2string($f\
name);}\n	  elsif ($action eq \"LSET\") {return st\
ring2file ($value, $fname, \">>\");}\n	  elsif ($a\
ction eq \"LRESET\") {return string2file ($value, \
$fname, \">\");}\n	  elsif ($action eq \"LRELEASE\\
") \n	    {\n	      if ( $debug_lock)\n		{\n		  my\
 $g=new FileHandle;\n		  open ($g, \">>$fname\");\\
n		  print $g \"\\nDestroyed by $$\\n\";\n		  clos\
e ($g);\n		  safe_system (\"mv $fname $fname.old\"\
);\n		}\n	      else\n		{\n		  unlink ($fname);\n	\
	}\n	    }\n	  return \"\";\n	}\n	\nsub file2strin\
g\n	{\n	  my $file=@_[0];\n	  my $f=new FileHandle\
;\n	  my $r;\n	  open ($f, \"$file\");\n	  while (\
<$f>){$r.=$_;}\n	  close ($f);\n	  return $r;\n	}\\
nsub string2file \n    {\n    my ($s,$file,$mode)=\
@_;\n    my $f=new FileHandle;\n    \n    open ($f\
, \"$mode$file\");\n    print $f  \"$s\";\n    clo\
se ($f);\n  }\n\nBEGIN\n    {\n      srand;\n    \\
n      $SIG{'SIGUP'}='signal_cleanup';\n      $SIG\
{'SIGINT'}='signal_cleanup';\n      $SIG{'SIGQUIT'\
}='signal_cleanup';\n      $SIG{'SIGILL'}='signal_\
cleanup';\n      $SIG{'SIGTRAP'}='signal_cleanup';\
\n      $SIG{'SIGABRT'}='signal_cleanup';\n      $\
SIG{'SIGEMT'}='signal_cleanup';\n      $SIG{'SIGFP\
E'}='signal_cleanup';\n      \n      $SIG{'SIGKILL\
'}='signal_cleanup';\n      $SIG{'SIGPIPE'}='signa\
l_cleanup';\n      $SIG{'SIGSTOP'}='signal_cleanup\
';\n      $SIG{'SIGTTIN'}='signal_cleanup';\n     \
 $SIG{'SIGXFSZ'}='signal_cleanup';\n      $SIG{'SI\
GINFO'}='signal_cleanup';\n      \n      $SIG{'SIG\
BUS'}='signal_cleanup';\n      $SIG{'SIGALRM'}='si\
gnal_cleanup';\n      $SIG{'SIGTSTP'}='signal_clea\
nup';\n      $SIG{'SIGTTOU'}='signal_cleanup';\n  \
    $SIG{'SIGVTALRM'}='signal_cleanup';\n      $SI\
G{'SIGUSR1'}='signal_cleanup';\n\n\n      $SIG{'SI\
GSEGV'}='signal_cleanup';\n      $SIG{'SIGTERM'}='\
signal_cleanup';\n      $SIG{'SIGCONT'}='signal_cl\
eanup';\n      $SIG{'SIGIO'}='signal_cleanup';\n  \
    $SIG{'SIGPROF'}='signal_cleanup';\n      $SIG{\
'SIGUSR2'}='signal_cleanup';\n\n      $SIG{'SIGSYS\
'}='signal_cleanup';\n      $SIG{'SIGURG'}='signal\
_cleanup';\n      $SIG{'SIGCHLD'}='signal_cleanup'\
;\n      $SIG{'SIGXCPU'}='signal_cleanup';\n      \
$SIG{'SIGWINCH'}='signal_cleanup';\n      \n      \
$SIG{'INT'}='signal_cleanup';\n      $SIG{'TERM'}=\
'signal_cleanup';\n      $SIG{'KILL'}='signal_clea\
nup';\n      $SIG{'QUIT'}='signal_cleanup';\n     \
 \n      our $debug_lock=$ENV{\"DEBUG_LOCK\"};\n  \
    \n      \n      \n      \n      foreach my $a \
(@ARGV){$CL.=\" $a\";}\n      if ( $debug_lock ){p\
rint STDERR \"\\n\\n\\n********** START PG: $PROGR\
AM *************\\n\";}\n      if ( $debug_lock ){\
print STDERR \"\\n\\n\\n**********(tcg) LOCKDIR: $\
LOCKDIR $$ *************\\n\";}\n      if ( $debug\
_lock ){print STDERR \"\\n --- $$ -- $CL\\n\";}\n \
     \n	     \n      \n      \n    }\nsub flush_er\
ror\n  {\n    my $msg=shift;\n    return add_error\
 ($EXIT_FAILURE,$$, $$,getppid(), $msg, $CL);\n  }\
\nsub add_error \n  {\n    my $code=shift;\n    my\
 $rpid=shift;\n    my $pid=shift;\n    my $ppid=sh\
ift;\n    my $type=shift;\n    my $com=shift;\n   \
 \n    $ERROR_DONE=1;\n    lock4tc ($rpid, \"LERRO\
R\",\"LSET\",\"$pid -- ERROR: $type\\n\");\n    lo\
ck4tc ($$, \"LERROR\",\"LSET\", \"$pid -- COM: $co\
m\\n\");\n    lock4tc ($$, \"LERROR\",\"LSET\", \"\
$pid -- STACK: $ppid -> $pid\\n\");\n   \n    retu\
rn $code;\n  }\nsub add_warning \n  {\n    my $rpi\
d=shift;\n    my $pid =shift;\n    my $command=shi\
ft;\n    my $msg=\"$$ -- WARNING: $command\\n\";\n\
    print STDERR \"$msg\";\n    lock4tc ($$, \"LWA\
RNING\", \"LSET\", $msg);\n  }\n\nsub signal_clean\
up\n  {\n    print dtderr \"\\n**** $$ (tcg) was k\
illed\\n\";\n    &cleanup;\n    exit ($EXIT_FAILUR\
E);\n  }\nsub clean_dir\n  {\n    my $dir=@_[0];\n\
    if ( !-d $dir){return ;}\n    elsif (!($dir=~/\
tmp/)){return ;}#safety check 1\n    elsif (($dir=\
~/\\*/)){return ;}#safety check 2\n    else\n     \
 {\n	`rm -rf $dir`;\n      }\n    return;\n  }\nsu\
b cleanup\n  {\n    #print stderr \"\\n----tc: $$ \
Kills $PIDCHILD\\n\";\n    #kill (SIGTERM,$PIDCHIL\
D);\n    my $p=getppid();\n    $CLEAN_EXIT_STARTED\
=1;\n    \n    \n    \n    if (&lock4tc($$,\"LERRO\
R\", \"LCHECK\", \"\"))\n      {\n	my $ppid=getppi\
d();\n	if (!$ERROR_DONE) \n	  {\n	    &lock4tc($$,\
\"LERROR\", \"LSET\", \"$$ -- STACK: $p -> $$\\n\"\
);\n	    &lock4tc($$,\"LERROR\", \"LSET\", \"$$ --\
 COM: $CL\\n\");\n	  }\n      }\n    my $warning=&\
lock4tc($$, \"LWARNING\", \"LREAD\", \"\");\n    m\
y $error=&lock4tc($$,  \"LERROR\", \"LREAD\", \"\"\
);\n    #release error and warning lock if root\n \
   \n    if (isrootpid() && ($warning || $error) )\
\n      {\n	\n	print STDERR \"**************** Sum\
mary *************\\n$error\\n$warning\\n\";\n\n	&\
lock4tc($$,\"LERROR\",\"RELEASE\",\"\");\n	&lock4t\
c($$,\"LWARNING\",\"RELEASE\",\"\");\n      } \n  \
  \n    \n    foreach my $f (@TMPFILE_LIST)\n     \
 {\n	if (-e $f){unlink ($f);} \n      }\n    forea\
ch my $d (@TMPDIR_LIST)\n      {\n	clean_dir ($d);\
\n      }\n    #No More Lock Release\n    #&lock4t\
c($$,\"LLOCK\",\"LRELEASE\",\"\"); #release lock \\
n\n    if ( $debug_lock ){print STDERR \"\\n\\n\\n\
********** END PG: $PROGRAM ($$) *************\\n\\
";}\n    if ( $debug_lock ){print STDERR \"\\n\\n\\
\n**********(tcg) LOCKDIR: $LOCKDIR $$ ***********\
**\\n\";}\n  }\nEND \n  {\n    \n    &cleanup();\n\
  }\n   \nsub blast_com2new_blast_com\n    {\n    \
  my $com=shift;\n      if ($ENV{\"NCBI_BLAST_4_TC\
OFFEE\"} eq \"OLD\"){return $com;}\n      elsif (!\
&pg_is_installed(\"legacy_blast.pl\")){return $com\
;}\n      else \n	{\n	  if ($com=~/formatdb/)\n	  \
  {\n	      $com=~s/formatdb/makeblastdb/;\n	     \
 $com=~s/\\-i/\\-in/;\n	      if ($com =~/pF/){$co\
m=~s/\\-pF/\\-dbtype nucl/;}\n	      if ($com =~/p\
 F/){$com=~s/\\-p F/\\-dbtype nucl/;}\n	      $com\
=\"$com -logfile /dev/null\";\n	      return $com;\
\n	    }\n	  elsif (&is_blast_package($com))\n	   \
 {\n	      my $path;\n	      \n	      if ( $ENV{\"\
NCBI_BIN_4_TCOFFEE\"}){$path=$ENV{\"NCBI_BLAST_4_T\
COFFEE\"};}\n	      else\n		{\n		  $path=`which le\
gacy_blast.pl`;\n		  $path=~s/\\/legacy_blast\\.pl\
//;\n		  chomp ($path);\n		}\n	      $path=\"--pat\
h $path\";\n	      if ( $com=~/\\>\\>/){$com=~s/\\\
>\\>/ $path \\>\\>/;}\n	      elsif ( $com=~/\\>/)\
{$com=~s/\\>/ $path \\>/;}\n	      else {$com.=\" \
$path\";}\n	      $com=\"legacy_blast.pl $com\";\n\
	      \n	      return $com;\n	    }\n	}\n    }\ns\
ub safe_system \n{\n  my $com=shift;\n  my $ntry=s\
hift;\n  my $ctry=shift;\n  my $pid;\n  my $status\
;\n  my $ppid=getppid();\n  if ($com eq \"\"){retu\
rn 1;}\n  \n  if ( ($com=~/^blast/) ||($com=~/^for\
matdb/)){$com=&blast_com2new_blast_com($com);} \n\\
n  if (($pid = fork ()) < 0){return (-1);}\n  if (\
$pid == 0)\n    {\n      set_lock($$, \" -SHELL- $\
com (tcg)\");\n      exec ($com);\n    }\n  else\n\
    {\n      lock4tc ($$, \"LLOCK\", \"LSET\", \"$\
pid\\n\");#update parent\n      $PIDCHILD=$pid;\n \
   }\n  if ($debug_lock){printf STDERR \"\\n\\t ..\
.. safe_system (fasta_seq2hmm)  p: $$ c: $pid COM:\
 $com\\n\";}\n\n  waitpid ($pid,WTERMSIG);\n\n  sh\
ift_lock ($pid,$$, \"LWARNING\",\"LWARNING\", \"LS\
ET\");\n\n  if ($? == $EXIT_FAILURE || lock4tc($pi\
d, \"LERROR\", \"LCHECK\", \"\"))\n    {\n      if\
 ($ntry && $ctry <$ntry)\n	{\n\n	  add_warning ($$\
,$$,\"$com failed [retry: $ctry out of $ntry]\");\\
n	  lock4tc ($pid, \"LRELEASE\", \"LERROR\", \"\")\
;\n	  #if ($com=~/EBI/){$com=~s/EBI/NCBI/;}\n	  #e\
lsif ($com=~/NCBI/){$com=~s/NCBI/EBI/;}\n	  \n	  r\
eturn safe_system ($com, $ntry, ++$ctry);\n	}\n   \
   elsif ($ntry == -1)\n	{\n	  if (!shift_lock ($p\
id, $$, \"LERROR\", \"LWARNING\", \"LSET\"))\n	   \
 {\n	      add_warning ($$,$$,\"$com failed\");\n	\
    }\n	  else\n	    {\n	      lock4tc ($pid, \"LR\
ELEASE\", \"LERROR\", \"\");\n	    }\n	  return $?\
;}\n      else\n	{\n	  if (!shift_lock ($pid,$$, \\
"LERROR\",\"LERROR\", \"LSET\"))\n	    {\n	      m\
yexit(add_error ($EXIT_FAILURE,$$,$pid,getppid(), \
\"UNSPECIFIED system\", $com));\n	    }\n	}\n    }\
\n  return $?;\n}\n\nsub check_configuration \n   \
 {\n      my @l=@_;\n      my $v;\n      foreach m\
y $p (@l)\n	{\n	  \n	  if   ( $p eq \"EMAIL\")\n	 \
   { \n	      if ( !($EMAIL=~/@/))\n		{\n		add_war\
ning($$,$$,\"Could Not Use EMAIL\");\n		myexit(add\
_error ($EXIT_FAILURE,$$,$$,getppid(),\"EMAIL\",\"\
$CL\"));\n	      }\n	    }\n	  elsif( $p eq \"INTE\
RNET\")\n	    {\n	      if ( !&check_internet_conn\
ection())\n		{\n		  myexit(add_error ($EXIT_FAILUR\
E,$$,$$,getppid(),\"INTERNET\",\"$CL\"));\n		}\n	 \
   }\n	  elsif( $p eq \"wget\")\n	    {\n	      if\
 (!&pg_is_installed (\"wget\") && !&pg_is_installe\
d (\"curl\"))\n		{\n		  myexit(add_error ($EXIT_FA\
ILURE,$$,$$,getppid(),\"PG_NOT_INSTALLED:wget\",\"\
$CL\"));\n		}\n	    }\n	  elsif( !(&pg_is_installe\
d ($p)))\n	    {\n	      myexit(add_error ($EXIT_F\
AILURE,$$,$$,getppid(),\"PG_NOT_INSTALLED:$p\",\"$\
CL\"));\n	    }\n	}\n      return 1;\n    }\n\n$pr\
ogram=\"T-COFFEE (Version_8.98)\";\n\n","*TC_METHO\
D_FORMAT_01\n******************generic_method.tc_m\
ethod*************\n*\n*       Incorporating new m\
ethods in T-Coffee\n*       Cedric Notredame 26/08\
/08\n*\n******************************************\
*************\n*This file is a method file\n*Copy \
it and adapt it to your need so that the method \n\
*you want to use can be incorporated within T-Coff\
ee\n**********************************************\
*********\n*                  USAGE               \
               *\n********************************\
***********************\n*This file is passed to t\
_coffee via -in:\n*\n*	t_coffee -in Mgeneric_metho\
d.method\n*\n*	The method is passed to the shell u\
sing the following\n*call:\n*<EXECUTABLE><PARAM1><\
IN_FLAG><seq_file><PARAM2><OUT_FLAG><outname><PARA\
M>\n*\n*Conventions:\n*<FLAG_NAME> 	<TYPE>		<VALUE\
>\n*<VALUE>:	no_name 	<=> Replaced with a space\n*\
<VALUE>:	&nbsp	<=> Replaced with a space\n*\n*****\
**************************************************\
\n*                  ALN_MODE                     \
      *\n*****************************************\
**************\n*pairwise   ->all Vs all (no self \
)[(n2-n)/2aln]\n*m_pairwise ->all Vs all (no self)\
[n^2-n]^2\n*s_pairwise ->all Vs all (self): [n^2-n\
]/2 + n\n*multiple   ->All the sequences in one go\
\n*\nALN_MODE		pairwise\n*\n**********************\
*********************************\n*              \
    OUT_MODE                           *\n********\
***********************************************\n*\
 mode for the output:\n*External methods: \n* aln \
-> alignmnent File (Fasta or ClustalW Format)\n* l\
ib-> Lib file (TC_LIB_FORMAT_01)\n*Internal Method\
s:\n* fL -> Internal Function returning a List (Li\
brairie)\n* fA -> Internal Function returning an A\
lignmnent\n*\nOUT_MODE		aln\n*********************\
**********************************\n*             \
     SEQ_TYPE                           *\n*******\
************************************************\n\
*G: Genomic, S: Sequence, P: PDB, R: Profile\n*Exa\
mples:\n*SEQTYPE	S	sequences against sequences (de\
fault)\n*SEQTYPE	S_P	sequence against structure\n*\
SEQTYPE	P_P	structure against structure\n*SEQTYPE	\
PS	mix of sequences and structure	\n*\nSEQ_TYPE	S\\
n*\n\n********************************************\
***********\n*                COMMAND LINE        \
                 *\n*EXECUTABLE PARAM1 IN_FLAG OUT\
_FLAG PARAM             *\n***********************\
********************************\n****************\
***************************************\n*        \
          EXECUTABLE                         *\n**\
**************************************************\
***\n*name of the executable\n*passed to the shell\
: executable\n*	\nEXECUTABLE	tc_generic_method.pl\\
n*\n**********************************************\
*********\n*                  IN_FLAG             \
                *\n*******************************\
************************\n*IN_FLAG\n*flag indicati\
ng the name of the in coming sequences\n*IN_FLAG S\
 no_name ->no flag\n*IN_FLAG S &bnsp-in&bnsp -> \"\
 -in \"\n*\nIN_FLAG		-infile=\n*\n****************\
***************************************\n*        \
          OUT_FLAG                           *\n**\
**************************************************\
***\n*OUT_FLAG\n*flag indicating the name of the o\
ut-coming data\n*same conventions as IN_FLAG\n*OUT\
_FLAG	S no_name ->no flag\n*if you want to redirec\
t, pass the parameters via PARAM1\n*set OUT_FLAG t\
o >\n*\nOUT_FLAG		-outfile=\n*\n******************\
*************************************\n*          \
        PARAM_1                              *\n**\
**************************************************\
***\n*<EXECUTABLE><PARAM1><IN_FLAG><seq_file><PARA\
M2><OUT_FLAG><outname><PARAM>\n*Parameters sent to\
 the EXECUTABLE and specified *before* IN_FLAG \n*\
If there is more than 1 PARAM line, the lines are\\
n*concatenated\n*Command_line: @EP@PARAM@-gapopen%\
e10%s-gapext%e20\n*	%s white space\n*	%e equal sig\
n\n*\n*PARAM1	\n*\n*\n*\n*************************\
******************************\n*                 \
 PARAM_2                              *\n*********\
**********************************************\n*<\
EXECUTABLE><PARAM1><IN_FLAG><seq_file><PARAM2><OUT\
_FLAG><outname><PARAM>\n*Parameters sent to the EX\
ECUTABLE and specified \n*after* IN_FLAG and *befo\
re* OUT_FLAG\n*If there is more than 1 PARAM line,\
 the lines are\n*concatenated\n*\n*PARAM1	\n*\n*\n\
**************************************************\
*****\n*                  PARAM                   \
           *\n************************************\
*******************\n*<EXECUTABLE><PARAM1><IN_FLAG\
><seq_file><PARAM2><OUT_FLAG><outname><PARAM>\n*Pa\
rameters sent to the EXECUTABLE and specified *aft\
er* OUT_FLAG\n*If there is more than 1 PARAM line,\
 the lines are\n*concatenated\n*\nPARAM	-mode=seq_\
msa -method=clustalw\nPARAM   -OUTORDER=INPUT -NEW\
TREE=core -align -gapopen=-15\n*\n****************\
***************************************\n*        \
          END                                *\n**\
**************************************************\
***\n","*TC_METHOD_FORMAT_01\n***************clust\
alw_method.tc_method*********\nEXECUTABLE	clustalw\
\nALN_MODE		pairwise\nIN_FLAG		-INFILE=\nOUT_FLAG	\
	-OUTFILE=\nOUT_MODE		aln\nPARAM		-gapopen=-10\nSE\
Q_TYPE		S\n***************************************\
**********\n","$VersionTag =                      \
                                                  \
                                                  \
         2.43;\nuse Env;\nuse FileHandle;\nuse Cwd\
;\nuse File::Path;\nuse Sys::Hostname;\nour $PIDCH\
ILD;\nour $ERROR_DONE;\nour @TMPFILE_LIST;\nour $E\
XIT_FAILURE=1;\nour $EXIT_SUCCESS=0;\n\nour $REFDI\
R=getcwd;\nour $EXIT_SUCCESS=0;\nour $EXIT_FAILURE\
=1;\n\nour $PROGRAM=\"extract_from_pdb\";\nour $CL\
=$PROGRAM;\n\nour $CLEAN_EXIT_STARTED;\nour $debug\
_lock=$ENV{\"DEBUG_LOCK\"};\nour $LOCKDIR=$ENV{\"L\
OCKDIR_4_TCOFFEE\"};\nif (!$LOCKDIR){$LOCKDIR=getc\
wd();}\nour $ERRORDIR=$ENV{\"ERRORDIR_4_TCOFFEE\"}\
;\nour $ERRORFILE=$ENV{\"ERRORFILE_4_TCOFFEE\"};\n\
&set_lock ($$);\nif (isshellpid(getppid())){lock4t\
c(getppid(), \"LLOCK\", \"LSET\", \"$$\\n\");}\n  \
    \nour $SILENT=\" >/dev/null 2>/dev/null\";\nou\
r $INTERNET=-1;\n\n\n\n\n\n\n\nour $BLAST_MAX_NRUN\
S=2;\nour $EXIT_SUCCESS=0;\nour $EXIT_FAILURE=1;\n\
our $CONFIGURATION=-1;\nour $REF_EMAIL=\"\";\nour \
$PROGRAM=\"extract_from_pdb\";\n\n\nmy %onelett_pr\
ot=&fill_onelett_prot();\nmy %threelett_prot=&fill\
_threelett_prot();\nmy %onelett_RNA=&fill_onelett_\
RNA();\nmy %threelett_RNA=&fill_threelett_RNA();\n\
my %onelett_DNA=&fill_onelett_DNA();\nmy %threelet\
t_DNA=&fill_threelett_DNA();\n\n\n\n\n\nmy %onelet\
t = (\n'P' => \\%onelett_prot,\n'D' => \\%onelett_\
DNA,\n'R' => \\%onelett_RNA\n);\n\n\nmy %threelett\
 = (\n'P' => \\%threelett_prot,\n'D' => \\%threele\
tt_DNA,\n'R' => \\%threelett_RNA\n);\n\n\n\n\n\n\n\
\nif($ARGV[0]=~/help/ ||$ARGV[0]=~/man/ || $ARGV[0\
]=~/HELP/ || $ARGV[0]=~/Man/ || $ARGV[0] eq \"-h\"\
  || $ARGV[0] eq \"-H\"  )\n{die \"SYNTAX: extract\
_from_pdb Version $VersionTag	\n	Minimum:         \
    [extract_from_pdb file] \n			   OR \n			     [\
... | extract_from_pdb]\n 	Flags (Default setting \
on the first line)\n	   -version..................\
.[Returns the Version Number]\n           -force..\
...................[Forces the file to be treated \
like a PDB file]\n                                \
      [Regenerates the header and SEQRES fields]\n\
           -force_name................[Forces the \
file to be named after name]]\n           -infile.\
....file...........[Flag can be omited]\n			      \
        [File must be pdb or fro pgm]\n           \
                           [File can also be compr\
essed Z or gz]\n                                  \
    [In the case of a compressed file, you can omi\
t the gz|Z extension]\n           -netfile........\
...........[File will be fetch from the net using \
wget]\n                                      [wget\
 or curl must be installed]\n                     \
                 [ftp://ftp.gnu.org/pub/gnu/wget/]\
\n                                      [http://cu\
rl.haxx.se/]\n                                    \
  [Must also be used to retrieve the file from a l\
ocal pdb copy (cf netaddress)]\n           -netadd\
ress................[Address used for the retrievi\
ng the netfile]\n                                 \
     [http://www.rcsb.org/pdb/cgi/export.cgi/%%.pd\
b.gz?format=PDB&pdbId=%%&compression=gz]\n        \
                              [http://www.expasy.c\
h/cgi-bin/get-pdb-entry.pl?%%]\n                  \
                    [local -> will get the file fr\
om pdb_dir (see pdb_dir)]\n           -netcompress\
ion............[Extension if the netfile comes com\
pressed]\n                                      [g\
z]\n           -pdb_dir...................[address\
 of the repertory where the pdb is installed]\n   \
                                   [Supports stand\
ard ftp style installation OR every stru in DIR]\n\
                                      [Give the ..\
../pdb/structure/ dir]\n                          \
            [If value omitted, the pg gets it from\
 the env variable PDB_DIR]\n           -netcompres\
sion_pg.........[gunzip]\n           -is_pdb_name.\
.........name.[Returns 1 if the name is a PDB ID, \
0 otherwise]\n           -model_type...........nam\
e.[Returns the model type if valid PDB name]\n    \
       -is_released_pdb_name name.[Returns 1 if th\
e name corresponds to a released PDB file]\n      \
     -get_pdb_chains.....name...[Returns the list \
of chains corresponding to the entry]\n           \
-get_pdb_id.........name...[Returns the PDB id wit\
hin the provided pdb file]\n           -get_fugue_\
name.....name...[Turns a name into a name valid fo\
r fugue]\n                                      [U\
ses the netaddress to do so]\n	   -chain......FIRS\
T..........[Extract the first chain only]\n		     \
  A B C..........[Extract Several chains if needed\
]\n		       ALL............[Extract all the chains\
]	\n           -ligand.....ALL............[Extract\
 the ligands in the chain (HETATM)]\n             \
          <name1>,<name2>[Extract All the named li\
gnds]\n	   -ligand_only...............[Extract onl\
y the ligands]\n           -ligand_list...........\
....[Extract the list of ligands]\n	   -coor......\
.<start>..<end>.[Coordinates of the fragment to ex\
tract]\n			              [Omit end to include the \
Cter]\n           -num........absolute.......[abso\
lute: relative to the seq] \n                     \
  file...........[file: relative to file]\n       \
    -num_out....new............[new: start 1->L]\n\
                       old............[old: keep t\
he file coordinates]\n           -delete.....<star\
t>..<end>.[Delete from residue start to residue en\
d]\n	   -atom.......CA.............[Atoms to inclu\
de, ALL for all of them]\n		       CA O N.........\
[Indicate several atoms if needed]\n	   -code.....\
..3..............[Use the 1 letter code or the 3 l\
etters code]\n	   -mode.......raw............[Outp\
ut original pdb file]\n                       pdb.\
...........[Output something that looks like pdb]\\
n		       fasta..........[Output the sequences in \
fasta format]\n		       simple.........[Output a f\
ormat easy to parse in C ]\n            -seq_field\
..ATOM...........[Field used to extract the sequen\
ce]\n		       SEQRES.........[Use the complete seq\
uence]\n	   -seq.......................[Equivalent\
 to  -mode fasta]\n	   -model......1..............\
[Chosen Model in an NMR file]\n           -nodiagn\
ostic..............[Switches Error Messages off]\n\
           -debug.....................[Sets the DE\
BUG ON]\n           -no_remote_pdb_dir.........[Do\
 not look for a remote file]\n           -cache_pd\
b.................[Cache Value, default is $HOME/.\
t_coffee/cache, other values: NO<=> No cache]\n\n \
     Environement Variables\n           These vari\
ables can be set from the environement\n          \
 Command line values with the corresponding flag s\
uperseed evironement value\n           NO_REMOTE_P\
DB_DIR..........[Prevents the program from searchi\
ng remote file: faster]\n           PDB_DIR.......\
.............[Indicates where PDB file must be fet\
ched (localy)]\n\n	 PROBLEMS: please contact cedri\
c.notredame\\@europe.com\\n\";\n	 exit ($EXIT_SUCC\
ESS);\n}\n\n$np=0;\n$n_para=$#ARGV;\n$model=1;\n$p\
db_dir=$ENV{'PDB_DIR'};if ($pdb_dir){$pdb_dir.=\"/\
\";}\n$debug=$ENV{'DEBUG_EXTRACT_FROM_PDB'};\n\n$n\
o_remote_pdb_dir=$ENV{NO_REMOTE_PDB_DIR};\n$HOME=$\
ENV{'HOME'};\nif ( $ENV{CACHE_4_TCOFFEE})\n{$cache\
=$ENV{CACHE_4_TCOFFEE};}\nelse\n{\n    $cache=\"$H\
OME/.t_coffee/cache/\";\n}\n\n   \n$netaddress=\"h\
ttp://www.rcsb.org/pdb/files/%%.pdb.gz\";\n$netcom\
pression_pg=\"gunzip\";\n$netcompression=\"gz\";\n\
\nforeach ($np=0; $np<=$n_para; $np++)\n  {       \
 \n    $value=$ARGV[$np];\n   \n    if  ($np==0 &&\
 !($value=~/^-.*/))\n      { \n	$pdb_file= $ARGV[$\
np];\n      }\n    elsif ( !($value=~/^-.*/))\n   \
   {\n	print \"@ARGV\";\n	die;\n      } \n    \n  \
  elsif ($value eq \"-nodiagnostic\"){$nodiagnosti\
c=1;}\n    elsif ($value eq \"-force\")\n      {\n\
	$force_pdb=1;\n      }\n    elsif ($value eq \"-f\
orce_name\")\n      {\n	$force_name=$ARGV[++$np];\\
n	$force_pdb=1;\n      }\n    \n    elsif ($value \
eq \"-is_pdb_name\")\n      {\n	$pdb_file= $ARGV[+\
+$np];	\n	$is_pdb_name=1;	\n      } \n    elsif ($\
value eq \"-is_released_pdb_name\")\n      {\n	$pd\
b_file= $ARGV[++$np];	\n	$is_released_pdb_name=1;\\
n      }\n    elsif ($value eq \"-model_type\")\n \
     {\n	$pdb_file= $ARGV[++$np];	\n	$model_type=1\
;\n      }\n    elsif ($value eq \"-debug\")\n{\n	\
$debug=1;\n}\n    elsif ($value eq \"-get_pdb_chai\
ns\")\n{\n	$pdb_file= $ARGV[++$np];\n	$get_pdb_cha\
ins=1;\n}\n    elsif ($value eq \"-get_pdb_ligands\
\")\n{\n	$get_pdb_ligands=1;\n}\n    \n    elsif (\
$value eq \"-get_pdb_id\")\n{\n	$pdb_file= $ARGV[+\
+$np];\n	$get_pdb_id=1;\n	\n}\n    \n    elsif ( $\
value eq \"-get_fugue_name\")\n{\n	$pdb_file= $ARG\
V[++$np];\n	$get_fugue_name=1;\n}\n    elsif ( $va\
lue eq \"-infile\")\n{\n       $pdb_file= $ARGV[++\
$np];\n} \n    elsif ($value eq \"-netfile\")\n{\n\
	$netfile=1;\n	if ( !($ARGV[$np+1]=~/^-.*/)){$pdb_\
file= $ARGV[++$np];}\n}\n    elsif (  $value eq \"\
-num\")\n{\n       $numbering= $ARGV[++$np];\n}\n \
   elsif (  $value eq \"-num_out\")\n{\n       $nu\
mbering_out= $ARGV[++$np];\n}\n    elsif ( $value \
eq \"-netaddress\")\n{\n	$netadress=$ARGV[++$np];\\
n}\n     \n    elsif ( $value eq \"-netcompression\
\")\n{\n	 $netcompression=$ARGV[++$np];\n}\n    el\
sif ( $value eq \"-pdb_dir\")\n{\n	 if ( !($ARGV[$\
np+1]=~/^-.*/)){$pdb_dir= \"$ARGV[++$np]/\";}\n}\n\
     elsif ( $value eq \"-no_remote_pdb_dir\")\n{\\
n	$no_remote_pdb_dir=1;\n	if ( !($ARGV[$np+1]=~/^-\
.*/)){$pdb_dir= \"$ARGV[++$np]/\";}\n}\n    elsif \
( $value eq \"-cache\")\n{\n	$cache=$ARGV[++$np];\\
n}\n    \n    elsif ($value eq \"-netcompression_p\
g\")\n{\n	  $netcompression_pg=$ARGV[++$np];\n}\n \
    elsif ($value eq \"-mode\")\n{\n       $MODE=$\
ARGV[++$np];\n}\n\n    elsif ( $value eq \"-model\\
")\n{\n       $model= $ARGV[++$np];\n}\n    elsif \
($value eq \"-seq_field\" )\n{\n       $seq_field=\
 $ARGV[++$np];\n}   \n    elsif ($value eq \"-coor\
\" )\n{\n       $start= $ARGV[++$np];\n  \n       \
if (($ARGV[$np+1] eq \"\") ||($ARGV[$np+1]=~/^-.*/\
)){$end=\"*\";} \n       else {$end=   $ARGV[++$np\
];}     \n       $coor_set=1;\n}\n    elsif ($valu\
e eq \"-delete\" )\n{\n       $delete_start= $ARGV\
[++$np];\n       $delete_end= $ARGV[++$np];\n     \
  $delete_set=1;\n}\n    elsif  ($value eq \"-code\
\")\n{\n       $code= $ARGV[++$np];\n}\n    elsif \
 ($value eq \"-no_hetatm\")\n{\n       $no_hetatm=\
1;\n}\n    elsif ($value eq \"-chain\")\n{\n      \
 while (!($ARGV[$np+1] eq \"\") &&!($ARGV[$np+1]=~\
/^-.*/))\n{\n	      ++$np;\n	      @c_chain=(@chai\
n,  $ARGV[$np]);\n	      $hc_chain{$ARGV[$np]}=$#c\
_chain+1;\n}           \n}\n    elsif ($value eq \\
"-atom\")\n{\n\n       while (!($ARGV[$np+1] eq \"\
\") && !($ARGV[$np+1]=~/^-.*/))\n{\n	      ++$np;\\
n	      $atom[$n_atom++]=  $ARGV[$np];\n	      $at\
om_list{$ARGV[$np]}=1;	      \n} \n       \n}\n   \
 elsif ( $value eq \"-unfold\")\n{\n	$unfold=1;\n}\
\n    elsif ($value eq \"-seq\" ||$value eq \"-fas\
ta\" )\n{\n       $MODE=\"fasta\";\n}\n    elsif (\
 $value eq \"-version\")\n{\n	print STDERR  \"\\ne\
xtract_from_pdb: Version $VersionTag\\n\";\n	&myex\
it ($EXIT_SUCCESS);\n}\n    elsif ( $value eq \"-l\
igand\")\n{\n	while (!($ARGV[$np+1] eq \"\") && !(\
$ARGV[$np+1]=~/^-.*/))\n{\n	    ++$np;\n	    $liga\
nd=1;\n	    $ligand_list{$ARGV[$np]}=1;	      \n} \
\n	$hc_chain{'LIGAND'}=1;\n}\n    elsif ( $value e\
q \"-ligand_only\")\n{\n	$ligand_only=1;\n}\n}\nif\
 ( $debug)\n{\n    print STDERR \"\\n[DEBUG:extrac\
t_from_pdb] NO_REMOTE_PDB_DIR: $no_remote_pdb_dir\\
\n\";\n    print STDERR \"\\n[DEBUG:extract_from_p\
db] PDB_DIR: $pdb_dir\\n\";\n}\n\n\nif ( $is_pdb_n\
ame)\n  {\n    if (&remote_is_pdb_name($pdb_file))\
\n      {\n	print \"1\";\n      }\n    else\n     \
 {\n	print \"0\";\n      }\n    exit ($EXIT_SUCCES\
S);\n  }\n\nif ( $is_released_pdb_name)\n  {\n    \
\n    if (&is_released($pdb_file))\n      {\n	prin\
t \"1\";\n      }\n    else\n      {\n	print \"0\"\
;\n      }\n    exit ($EXIT_SUCCESS);\n  }\nif ($m\
odel_type)\n  {\n   \n    printf \"%s\", &pdb2mode\
l_type($pdb_file);\n    exit ($EXIT_SUCCESS);\n   \
 \n  }\n    \n\nif (!$force_name)\n{\n    $pdb_fil\
e=~/([^\\/]*)$/;\n    $force_name=$1;\n}\n\n$local\
_pdb_file=$pdb_file;\n\nif ( $debug){print STDERR \
\"\\n[DEBUG: extract_from_pdb] Scan For $local_pdb\
_file\\n\";}\n\n$mem=$no_remote_pdb_dir;\n$no_remo\
te_pdb_dir=1;\n$tmp_pdb_file=get_pdb_file ($local_\
pdb_file);\n\nif ( !-e $tmp_pdb_file || $tmp_pdb_f\
ile eq \"\")\n  {\n    $local_pdb_file=$pdb_file;\\
n    ($local_pdb_file, $suffix_chain)=&pdb_name2na\
me_and_chain($local_pdb_file);\n\n    if ($local_p\
db_file)\n      {\n	if ( $debug){print STDERR \"\\\
nSplit $pdb_file into $local_pdb_file and $suffix_\
chain \\n\";}\n	$tmp_pdb_file=get_pdb_file ($local\
_pdb_file);\n	if ( $tmp_pdb_file ne \"\")\n	  {\n	\
    @c_chain=();\n	    @c_chain=($suffix_chain);\n\
	    %hc_chain=();\n	    $hc_chain{$suffix_chain}=\
1;\n	  }\n      }\n  }\n\n$no_remote_pdb_dir=$mem;\
\nif ($no_remote_pdb_dir==0)\n  {\n    \n    if ( \
!-e $tmp_pdb_file || $tmp_pdb_file eq \"\")\n     \
 {\n	\n	$local_pdb_file=$pdb_file;\n	($local_pdb_f\
ile, $suffix_chain)=&pdb_name2name_and_chain($loca\
l_pdb_file);\n	if ($local_pdb_file)\n	  {\n	    \n\
	    if ( $debug){print STDERR \"\\nSplit $pdb_fil\
e into $local_pdb_file and $suffix_chain \\n\";}\n\
	    \n	    $tmp_pdb_file=get_pdb_file ($local_pdb\
_file);    \n	    \n	    if ( $tmp_pdb_file ne \"\\
")\n	      {\n		@c_chain=();\n		@c_chain=($suffix_\
chain);\n		%hc_chain=();\n		$hc_chain{$suffix_chai\
n}=1;\n	      }\n	  }\n      }\n  }\n\nif ( $debug\
){print STDERR \"\\n$pdb_file copied into ##$tmp_p\
db_file##\\n\";}\n\nif ( !-e $tmp_pdb_file || $tmp\
_pdb_file eq \"\")\n{\n	if ($is_pdb_name)\n{\n	   \
 print \"0\\n\"; exit ($EXIT_SUCCESS);\n}\n	else\n\
{\n  \n	    print \"\\nEXTRACT_FROM_PDB: NO RESULT\
 for $pdb_file\\n\";\n	    &myexit ($EXIT_SUCCESS)\
;	\n}\n}\n\n\n\n\n%molecule_type=&pdbfile2chaintyp\
e($tmp_pdb_file);\nif ( $debug){print STDERR \"\\n\
\\tSequence Type determined\\n\";}\n\n$pdb_id=&get\
_pdb_id ($tmp_pdb_file);\nif ( $debug){print STDER\
R \"\\n\\tID: $pdb_id (for $tmp_pdb_file)\\n\";}\n\
\nif ( $pdb_id eq \"\"){$pdb_id=$force_name;}\n\n@\
f_chain=&get_chain_list ($tmp_pdb_file);\nif ( $de\
bug){print STDERR \"\\n\\tChain_list:@f_chain\\n\"\
;}\n\nif ( $get_pdb_chains)\n{\n    print \"@f_cha\
in\\n\";\n    &myexit ($EXIT_SUCCESS);\n}\nif ( $g\
et_pdb_ligands)\n{\n    %complete_ligand_list=&get\
_ligand_list ($tmp_pdb_file);\n    print $complete\
_ligand_list{\"result\"};\n    &myexit ($EXIT_SUCC\
ESS);\n}\n\nelsif ( $get_pdb_id ||$get_fugue_name \
)\n{\n    if    (@c_chain && $c_chain[0] eq \"FIRS\
T\"){$pdb_id=$pdb_id.$f_chain[0];}\n    elsif (@c_\
chain && $c_chain[0] ne \" \"){$pdb_id=$pdb_id.$c_\
chain[0];}\n    \n    print \"$pdb_id\\n\";\n    &\
myexit ($EXIT_SUCCESS);\n    \n}\nelsif ( $is_pdb_\
name)\n{\n    printf \"1\\n\";\n    &myexit ($EXIT\
_SUCCESS);\n}\n\n\n\n$structure_file=vtmpnam();\n\\
nif ( $debug){print STDERR \"\\n\\tCheck_point #1:\
 $tmp_pdb_file  $structure_file\\n\";}\n\n$INFILE=\
vfopen (\"$tmp_pdb_file\", \"r\"); \n$TMP=vfopen (\
\"$structure_file\", \"w\");\n\n$print_model=1;\n$\
in_model=0;\n\nif ( $debug){print STDERR \"\\n\\tC\
heck_point #2\\n\";}\nwhile ( <$INFILE>)\n{\n  my \
$first_model=0;\n  $line=$_;\n\n  if ( !$first_mod\
el && ($line =~/^MODEL\\s*(\\d*)/))\n    {\n      \
$first_model=$1;\n      if ($model==1){$model=$fir\
st_model;}\n    }\n  \n  if (($line =~/^MODEL\\s*(\
\\d*)/))\n    {\n      if ($1==$model)\n	{\n	  $in\
_model=1;\n	  $print_model=1;\n	  $is_nmr=1;\n	}\n\
      elsif ( $in_model==0)\n	{\n	  $print_model=0\
;\n	}\n      elsif ( $in_model==1)\n	{\n	  last;\n\
	}\n    }\n  if ($print_model){print $TMP $line;} \
 \n}\nclose ($TMP);\nclose ($INFILE);\n\nif ( $deb\
ug){print STDERR \"\\n\\tCheck_point #3\\n\";}	\n\\
n  if ($numbering eq \"\"){$numbering=\"absolute\"\
;}\n  if ($numbering_out eq \"\"){$numbering_out=\\
"new\";}\n\n  if ( $delete_set && $coor_set) {die \
\"-delete and -coor are mutually exclusive, sorry\\
\n\";}\n  if ( $n_atom==0){$atom_list[$n_atom++]=\\
"ALL\";$atom_list{$atom_list[0]}=1;}\n  if ( $seq_\
field eq \"\"){$seq_field=\"ATOM\";}\n  \n  if ( $\
MODE eq \"\"){$MODE=\"pdb\";}\n  elsif ( $MODE eq \
\"simple\" && $code==0){$code=1;}\n\n  if ( $code=\
=0){$code=3;}\n\n\nif ($f_chain[0] eq \" \"){$hc_c\
hain{' '}=1;$c_chain[0]=\" \";}\nelsif (!@c_chain)\
{$hc_chain{FIRST}=1;$c_chain[0]=\"FIRST\";}#make s\
ure the first chain is taken by default\n\nif    (\
$hc_chain{ALL}) \n{\n      @c_chain=@f_chain;\n   \
   foreach $e (@c_chain){$hc_chain{$e}=1;}\n}\nels\
if($hc_chain{FIRST})\n{\n	@c_chain=($f_chain[0]);\\
n	$hc_chain{$f_chain[0]}=1;\n}\n\n\n$MAIN_HOM_CODE\
=&get_main_hom_code ($structure_file);\n$INFILE=vf\
open ($structure_file, \"r\");\n\n\nif ( $MODE eq \
\"raw_pdb\" || $MODE eq \"raw\")\n{\n    while (<$\
INFILE>)\n{	print \"$_\";}\n    close ( $INFILE);\\
n    &myexit($EXIT_SUCCESS);\n}    \nif ( $MODE eq\
 \"raw4fugue\" )\n{\n    while (<$INFILE>)\n{	\n	$\
l=$_;\n	if ($l=~/^SEQRES/)\n{\n	    \n	    $c= sub\
str($l,11,1);\n	    if ($hc_chain {$c}){print \"$l\
\";}\n}\n	elsif ( $l=~/^ATOM/)\n{\n	    $c=substr(\
$l,21,1);\n	    if ($hc_chain {$c}){print \"$l\";}\
\n}\n}\n    close ( $INFILE);\n    &myexit($EXIT_S\
UCCESS);\n}    \n\nif ( $MODE eq \"pdb\")\n{\n\n  \
  $read_header=0;\n    while (<$INFILE>) \n{\n	   \
 $line=$_;\n	    if    ($line =~ /^HEADER/){print \
\"$line\";$read_header=1;}\n}\n    close ($INFILE)\
;\n\n    if (!$read_header)\n{\n	print \"HEADER   \
 UNKNOWN                                 00-JAN-00\
   $force_name\\n\";\n}\n\n    $INFILE=vfopen ($st\
ructure_file, \"r\");\n    \n    print \"COMPND   \
1 CHAIN:\";\n    $last=pop(@c_chain);\n    foreach\
 $c ( @c_chain){ print \" $c,\";}\n    if ( $last \
eq \" \"){print \" NULL;\\n\";}\n    else \n{\n   \
   print \" $last;\\n\";\n}\n    @c_chain=(@c_chai\
n, $last);\n    \n    print \"REMARK Output of the\
 program extract_from_pdb (Version $VersionTag)\\n\
\";\n    print \"REMARK Legal PDB format not Guara\
nteed\\n\";\n    print \"REMARK This format is not\
 meant to be used in place of the PDB format\\n\";\
\n    print \"REMARK The header refers to the orig\
inal entry\\n\";\n    print \"REMARK The sequence \
from the original file has been taken in the field\
: $seq_field\\n\";\n    print \"REMARK extract_fro\
m_pdb, 2001, 2002, 2003, 2004, 2005 2006 (c) CNRS \
and Cedric Notredame\\n\";   \n    if ( $coor_set)\
\n{\n       print \"REMARK Partial chain: Start $s\
tart End $end\\n\";\n}\n    if ( $is_nmr)\n{\n    \
   print \"REMARK NMR structure: MODEL $model\\n\"\
;\n}\n   \n    if ( $n_atom!=0)\n{\n       print \\
"REMARK Contains Coordinates of: \";\n       forea\
ch $a (@atom){print \"$a \";}\n       print \"\\n\\
";\n}  \n}\n\n\n\n\nmy $residue_index = -999;\nmy \
$old_c = \"TemporaryChain\";\n\nwhile (<$INFILE>) \
\n{\n	$line=$_;\n\n\n	if ($line =~ /^SEQRES/)\n{\n\
\n		@field=/(\\S*)\\s*/g;\n\n		$c= substr($_,11,1)\
;\n\n		\n		$l=$#field;\n		for ($a=4; $a<$#field ;)\
\n{\n			if (!$onelett{$molecule_type{$c}}->{$field\
[$a]})\n{\n				splice @field, $a, 1;\n}\n			else \\
n{\n				$a++;\n}\n}\n	\n		if ( $c ne $in_chain)\n{\
\n			$pdb_chain_list[$n_pdb_chains]=$c;\n			$pdb_c\
hain_len [$n_pdb_chains]=$len;\n			$in_chain=$c;\n\
			$n_pdb_chains++;\n}\n	\n		for ( $a=4; $a<$#fiel\
d;$a++)\n{\n			@{$complete_seq{$c}}->[$complete_se\
q_len{$c}++]=$field[$a];\n}\n}\n    elsif ( $line=\
~/^ATOM/ || ($line=~/^HETATM/ && &is_aa(substr($li\
ne,17,3),substr($line,21,1)) && !$no_hetatm))\n{\n\
\n	 \n    $RAW_AT_ID=$AT_ID=substr($line,12,4);\n	\
$RES_ID=&is_aa(substr($line,17,3),substr($line,21,\
1));\n	$CHAIN=substr($line,21,1);\n\n    $RES_NO=s\
ubstr($line,22,4);\n	$HOM_CODE=substr ($line, 26, \
1);\n	$TEMP=substr($line,60,6);\n	\n	$TEMP=~s/\\s/\
/g;\n        $AT_ID=~s/\\s//g;\n	$RES_ID=~s/\\s//g\
;\n        $RES_NO=~s/\\s//g;\n		\n	if ( $HOM_CODE\
 ne $MAIN_HOM_CODE){next;}\n	elsif ( $already_read\
2{$CHAIN}{$RES_ID}{$AT_ID}{$RES_NO}){next;}\n	else\
{$already_read2{$CHAIN}{$RES_ID}{$AT_ID}{$RES_NO}=\
1;}\n	\n	\n	if ($coor_set && $numbering eq \"file\\
" && $residue_index ne $RES_NO)\n{\n	    \n	    if\
 ( $RES_NO<=$start){$real_start{$CHAIN}++;}\n	    \
if ( $RES_NO<=$end){$real_end{$CHAIN}++;}\n}\n	els\
if ($numbering eq \"absolute\")\n{\n	    $real_sta\
rt{$CHAIN}=$start;\n	    $real_end{$CHAIN}=$end;\n\
}\n\n        $KEY=\"ALL\";\n        if ( $CHAIN ne\
 $in_atom_chain)\n{\n	    \n	  $pdb_atom_chain_lis\
t[$n_pdb_atom_chains]=$c;\n	  $pdb_atom_chain_len \
[$n_pdb_atom_chains]=$len;\n	  $in_atom_chain=$c;\\
n	  $n_pdb_atom_chains++;\n}\n	\n	if ( $residue_in\
dex ne $RES_NO)\n{\n	     $residue_index = $RES_NO\
;\n	     @{$atom_seq{$CHAIN}}->[$atom_seq_len{$CHA\
IN}++]=$RES_ID;;\n}\n}\n\n}\nclose ($INFILE);\n\n\\
n\n\n\n\n$INFILE=vfopen ($structure_file, \"r\");\\
nforeach $c (@c_chain)\n{\n\n	if    ( $seq_field e\
q \"SEQRES\"){@pdb_seq=@{$complete_seq{$c}};}\n	el\
sif ( $seq_field eq \"ATOM\")  {@pdb_seq=@{$atom_s\
eq{$c}};}\n	\n\n	$full_length=$l=$#pdb_seq+1;\n		\\
n	if ( $real_end{$c}==\"*\"){$real_end{$c}=$full_l\
ength;}\n	if ( $coor_set)\n{	   \n\n	   if ( $real\
_end{$c} < $l){splice @pdb_seq, $real_end{$c}, $l;\
}\n	   if ( $real_start{$c} < $l){splice @pdb_seq,\
 0, $real_start{$c}-1;}	  	   \n	   $l=$#pdb_seq;\\
n}\n\n	elsif ( $delete_set)\n{\n	   splice @pdb_se\
q, $delete_start, $delete_end-$delete_start+1;\n	 \
  $l=$#pdb_seq;\n}\n	\n	$new_fasta_name=\"$pdb_id$\
c\";\n	if ( $coor_set)\n{\n	   if ( $n_pdb_chains=\
=0){$new_fasta_name=\"$new_fasta_name$c\";}\n	   $\
new_fasta_name= $new_fasta_name.\"\\_$start\\_$end\
\";\n}\n	   \n	if ( $MODE eq \"pdb\")\n{\n	   $nl=\
1;\n	   $n=0;\n	   \n	   foreach $res ( @pdb_seq)\\
n		{\n		if ( !$n)\n		{\n		\n		 printf \"SEQRES %3d\
 %1s %4d  \", $nl,$c, $l;\n		 $nl++;\n	}\n	     $r\
es=~s/\\s//g;\n	     \n	     if ($code==1){ printf\
 \"%3s \",$onelett{$molecule_type{$c}}->{$res};}\n\
	     elsif  ($code==3){ printf \"%3s \",$res};\n	\
     \n	     $n++;		  \n	     if ( $n==13){$n=0;pr\
int \"\\n\";}\n}\n	  if ( $n!=0){print \"\\n\"; $n\
=0;}\n}\n	elsif ( $MODE eq \"simple\")\n{\n	  prin\
t \"# SIMPLE_PDB_FORMAT\\n\";\n	  if ( $new_fasta_\
name eq \" \"){$new_fasta_name=\"dummy_name\";}\n	\
  print \">$new_fasta_name\\n\";\n\n	  foreach $re\
s ( @pdb_seq)\n{\n	      print \"$onelett{$molecul\
e_type{$c}}->{$res}\";\n}\n	  print \"\\n\";\n}\n	\
elsif ( $MODE eq \"fasta\")\n{\n	  $n=0;\n	  print\
 \">$new_fasta_name\\n\";\n	  \n	  foreach $res ( \
@pdb_seq)\n{\n\n	      print \"$onelett{$molecule_\
type{$c}}->{$res}\";\n              $n++;\n	      \
if ( $n==60){print \"\\n\"; $n=0;}\n}\n	  print \"\
\\n\"; \n}\n}\n\nif ( $MODE eq \"fasta\")\n{\n    \
 &myexit($EXIT_SUCCESS);\n  \n}\n\n  \n  $charcoun\
t=0;\n  $inchain=\"BEGIN\";\n  $n=0;\n  while (<$I\
NFILE>) \n{\n    $line=$_;\n     \n    if ($line =\
~/^ATOM/  ||  ($line=~/^HETATM/))\n{\n	$line_heade\
r=\"UNKNWN\";\n	$RES_ID=substr($line,17,3);\n	$cha\
in = substr($line,21,1);\n\n	if ($line =~/^ATOM/)\\
n{\n	    $line_header=\"ATOM\";\n	    $RES_ID=(&is\
_aa($RES_ID,$chain))?&is_aa($RES_ID,$chain):$RES_I\
D;\n}\n	elsif ($line=~/^HETATM/ && ($ligand_list {\
$RES_ID} ||$ligand_list {'ALL'} || !&is_aa($RES_ID\
,$chain)))\n{\n	    $line_header=\"HETATM\";\n}\n	\
elsif ($line=~/^HETATM/ && (&is_aa($RES_ID,$chain)\
 && !$no_hetatm))\n{\n	    $line_header=\"ATOM\";\\
n	    $RES_ID=&is_aa($RES_ID,$chain);\n}\n	else\n{\
\n	    next;\n}\n\n	\n\n	$X=substr($line,30,8);   \
  \n	$Y=substr($line,38,8);\n	$Z=substr($line,46,8\
);\n	$TEMP=substr($line,60,6);\n	\n	$RAW_AT_ID=$AT\
_ID=substr($line,12,4);\n	$CHAIN=substr($line,21,1\
);\n	$RES_NO=substr($line,22,4);\n	$HOM_CODE=subst\
r ($line, 26, 1);\n	\n	$X=~s/\\s//g;\n	$Y=~s/\\s//\
g;\n	$Z=~s/\\s//g;\n	$TEMP=~s/\\s//g;\n	\n	$AT_ID=\
~s/\\s//g;\n	$RES_ID=~s/\\s//g;\n	$RES_NO=~s/\\s//\
g;\n\n	\n	if ( $HOM_CODE ne $MAIN_HOM_CODE){next;}\
\n	elsif ( $already_read{$CHAIN}{$RES_ID}{$AT_ID}{\
$RES_NO}){next;}\n	else{$already_read{$CHAIN}{$RES\
_ID}{$AT_ID}{$RES_NO}=1;}\n	\n	$KEY=\"ALL\";\n\n  \
    	if ( $RES_NO ==0){$start_at_zero=1;}\n\n	$RES\
_NO+=$start_at_zero;    \n	\n	if ( $current_chain \
ne $CHAIN)\n{\n	    $current_chain=$CHAIN;\n	    $\
pos=$current_residue=0;\n	    $offset=($coor_set)?\
($real_start{$CHAIN}-1):0;\n	    if    ( $seq_fiel\
d eq \"SEQRES\"){@ref_seq=@{$complete_seq{$CHAIN}}\
;}\n	    elsif ( $seq_field eq \"ATOM\")  {@ref_se\
q=@{$atom_seq{$CHAIN}};}\n}\n	\n	if ($current_resi\
due != $RES_NO)\n{\n	    $current_residue=$RES_NO;\
\n	    if    ( $seq_field eq \"SEQRES\"){$pos=$cur\
rent_residue;}\n	    elsif ( $seq_field eq \"ATOM\\
"){$pos++;}\n}\n	\n	\n	if ($n_atom==0 || $atom_lis\
t{$AT_ID}==1 || $atom_list{$KEY}==1)\n{ 	\n	    \n\
	    $do_it=(!@c_chain || $hc_chain{$CHAIN} ||$hc_\
chain{'LIGAND'} );\n	    \n	    $do_it= ($do_it==1\
) && ($coor_set==0 ||($pos>=$real_start{$CHAIN} &&\
 $pos<=$real_end{$CHAIN}));\n	    $do_it= ($do_it=\
=1) && ($delete_set==0 || $pos<$delete_start ||$po\
s>$delete_end );\n	    if ($ligand==0 && $line_hea\
der eq \"HETATM\" ){$do_it=0;}\n	    if ($ligand_o\
nly==1 && $line_header eq \"ATOM\" ){$do_it=0;}\n	\
    if ($ligand==1 && $line_header eq \"HETATM\" &\
& $ligand_list{$RES_ID}==0 && $ligand_list{\"ALL\"\
}==0){$do_it=0;} \n	    \n	    \n	    if ( $do_it)\
\n{\n		$n++;\n		$out_pos=$pos;\n		\n	       if ( $\
delete_set)\n{\n		  if ( $out_pos< $delete_start){\
;}\n		  else {$offset=$delete_end-$delete_start;}\\
n}       \n	       \n	       if ( $numbering_out e\
q \"new\"){$out_pos-=$offset;}\n	       elsif ( $n\
umbering_out eq \"old\"){$out_pos=$RES_NO;}\n	    \
   \n       \n	       \n	       if ( $code==1){$RE\
S_ID=$onelett{$molecule_type{$c}}->{$RES_ID};}\n	 \
   \n	       if ($unfold)\n{\n		   $unfolded_x+=5;\
\n		   $X=$unfolded_x;\n		   $Y=0;\n		   $Z=0;\n		\
   $float=1;\n}\n	       else\n{\n		   $float=3;\n\
}\n\n	       if ( $MODE eq \"pdb\")\n{\n		   print\
f \"%-6s%5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f  1.0\
0 %5.2f\\n\",$line_header, $n, $RAW_AT_ID,$RES_ID,\
$CHAIN,$out_pos, $X, $Y, $Z,$TEMP;		  \n}\n	      \
 elsif ( $MODE eq \"simple\")\n{\n		    if ( $RES_\
ID eq \"\"){$RES_ID=\"X\";}\n		  printf \"%-6s %5s\
 %s %2s %4d    %8.3f %8.3f %8.3f\\n\",$line_header\
, $AT_ID, $RES_ID,($CHAIN eq\"\" || $CHAIN eq \" \\
")?\"A\":$CHAIN,$out_pos, $X, $Y, $Z,$TEMP;\n}\n\n\
}\n}\n}\n}\nprint \"\\n\";\nclose($INFILE);\n\n\ni\
f ( $error ne \"\") \n{$error=$error.\"\\nDiagnost\
ic:    SEQRES and the residues in ATOM are probabl\
y Incompatible\\n\";\n    $error=$error.  \"Recome\
ndation: Rerun with '-fix 1' in order to ignore th\
e SEQRES sequences\\n\";\n}\nif (!$nodiagnostic){p\
rint STDERR $error;}\n&myexit ( $EXIT_SUCCESS);\n\\
nsub is_released \n  {\n    my ($r);\n    my $in=@\
_[0];\n    my $name=&remote_is_pdb_name ($in);\n  \
  my $hold=&remote_is_on_hold($in);\n    \n    $r=\
($name && !$hold)?1:0;\n    return $r;\n  }\nsub r\
emote_is_pdb_name \n{\n    my $in=@_[0];\n    my (\
$ref_file, $pdb);\n    my ($value,$value1,$value2)\
;\n\n    if ( $in=~/[^\\w\\d\\:\\_]/){return 0;}\n\
    $ref_file=\"$cache/pdb_entry_type.txt\";\n    \
\n    if ( !-e $ref_file || (-M $ref_file)>2 || -z\
 $ref_file)\n      {\n	&url2file(\"ftp://ftp.wwpdb\
.org/pub/pdb/derived_data/pdb_entry_type.txt\", $r\
ef_file);\n      }\n    $pdb=substr ($in,0, 4);\n \
   chomp(($value1=`grep -c $pdb $ref_file`));\n   \
 $pdb=lc($pdb);\n    chomp(($value2=`grep -c $pdb \
$ref_file`));\n    $value=($value1 || $value2)?1:0\
;\n    $value=($value>0)?1:0;\n    \n    return $v\
alue;\n  }\n\nsub pdb2model_type\n{\n    my $in=@_\
[0];\n    my ($ref_file, $pdb);\n    my ($value, $\
ret);\n\n    if ( $in=~/[^\\w\\d\\:\\_]/){return 0\
;}\n    $ref_file=\"$cache/pdb_entry_type.txt\";\n\
    \n    if ( !-e $ref_file || (-M $ref_file)>2 |\
| -z $ref_file)\n      {\n	&url2file(\"ftp://ftp.w\
wpdb.org/pub/pdb/derived_data/pdb_entry_type.txt\"\
, $ref_file);\n      }\n    $pdb=substr ($in,0, 4)\
;\n    $pdb=lc($pdb);\n    \n    chomp(($value=`gr\
ep $pdb $ref_file`));\n    \n    $value=~/^\\S+\\s\
+\\S+\\s+(\\S+)/;\n    $ret=$1;\n    if ( $ret eq\\
"\"){return \"UNKNOWN\";}\n    \n    return $ret;\\
n  }\nsub remote_is_on_hold\n  {\n    my $in=@_[0]\
;\n    my ($ref_file, $pdb);\n    my ($value1, $va\
lue2,$value);\n    \n    if ( $in=~/[^\\w\\d\\:\\_\
]/){return 0;}\n    $ref_file=\"$cache/unreleased.\
xml\";\n    \n    if ( !-e $ref_file || (-M $ref_f\
ile)>2 || -z $ref_file)\n      {\n	&url2file(\"htt\
p://www.rcsb.org/pdb/rest/getUnreleased\",$ref_fil\
e);\n      }\n    \n    $pdb=substr ($in,0, 4);\n \
   chomp(($value1=`grep -c $pdb $ref_file`));\n   \
 $pdb=lc($pdb);\n    chomp(($value2=`grep -c $pdb \
$ref_file`));\n    $value=($value1 || $value2)?1:0\
;\n    $value=($value>0)?1:0;\n    return $value;\\
n  }\nsub is_pdb_file\n{\n    my @arg=@_;\n\n    i\
f ( !-e $arg[0]){return 0;}\n    \n    $F=vfopen (\
$arg[0], \"r\");\n    while ( <$F>)\n{\n	if (/^HEA\
DER/)\n{\n	    close $F;\n	    return 1;\n}\n	elsi\
f ( /^SEQRES/)\n{\n	    close $F;\n	    return 1;\\
n}\n	elsif ( /^ATOM/)\n{\n	    close $F;\n	    ret\
urn 1;\n}\n}\n    return 0;\n}\nsub get_pdb_id\n{\\
n    my $header_file=@_[0];\n    my $id;\n    my $\
F= new FileHandle;\n    \n    \n    $F=vfopen (\"$\
header_file\", \"r\");\n\n    while ( <$F>)\n     \
 {\n	if ( /HEADER/)\n	  {\n	    if ($debug){print \
\"$_\";}\n	    $id=substr($_,62,4 );\n	    return \
$id;\n	  }\n      }\n    close ($F);\n    \n    re\
turn \"\";\n}\n\nsub get_ligand_list\n{\n    my $p\
db_file=@_[0];\n    my $chain;\n    my $ligand;\n \
   my %complete_ligand_list;\n    \n\n    $F=vfope\
n ($pdb_file, \"r\");\n    while ( <$F>)\n{\n	if (\
 /^HETATM/)\n{\n	    $line=$_;\n	    $chain=substr\
($line,21,1);\n	    $ligand=substr($line,17,3);\n	\
    \n	    if (!$complete_ligand_list{$chain}{$lig\
and})\n{\n		\n		$complete_ligand_list{\"result\"}.\
=\"CHAIN $chain LIGAND $ligand\\n\";\n		$complete_\
ligand_list{$chain}{$ligand}=1;\n}\n}\n}\n    clos\
e ($F);\n    return %complete_ligand_list;\n}\n\ns\
ub get_chain_list \n{\n    my $header_file;\n    m\
y @chain_list;\n    my @list;\n    my $n_chains;\n\
    my %chain_hasch;\n    my $pdb_file=@_[0];\n   \
 my $c;\n    my %hasch;\n    my $chain;\n  \n    \\
n    $F=vfopen ($pdb_file, \"r\");\n    while ( <$\
F>)\n{\n\n\n	if (/SEQRES\\s+\\d+\\s+(\\S+)/)\n	  {\
\n	    $chain = substr($_,11,1);$chain=~s/\\s//g;i\
f ( $chain eq \"\"){$chain=\" \";}\n	    if (!$has\
ch{$chain}){$hasch{$chain}=1;push @chain_list, $ch\
ain;}\n	  }\n	if (/^ATOM/ || /^HETATM/)\n	  {\n	  \
  $chain = substr($_,21,1); $chain=~s/\\s//g;if ( \
$chain eq \"\"){$chain=\" \";}\n	    if (!$hasch{$\
chain}){$hasch{$chain}=1;push @chain_list, $chain;\
}\n	  }\n      }\n\n\nclose ($F);\nif (!@chain_lis\
t)\n  {\n    @chain_list=(\"A\");\n  }\n\n\nreturn\
 @chain_list;\n}\n\nsub token_is_in_list\n{\n\n   \
 my @list=@_;\n    my $a;\n    \n    for ($a=1; $a\
<=$#list; $a++)\n{\n	if ( $list[$a] eq $list[0]){r\
eturn $a;}\n}\n}\n\nsub pdb_name2name_and_chain \n\
{\n    my $pdb_file=@_[0];\n    my $pdb_file_in;\n\
    my @array;\n    my $chain;\n    my $c;\n\n    \
$pdb_file_in=$pdb_file;\n\n    $pdb_file=~/^(.{4})\
/;$pdb_id=$1;\n    @array=($pdb_file=~/([\\w])/g);\
\n  \n  \n    $chain=uc ($array[4]);\n    $chain=(\
$chain eq \"\")?\"FIRST\":$chain;\n    \n    retur\
n ( $pdb_id, $chain);\n\n    if ( $#array==3){retu\
rn ($pdb_id, \"FIRST\");}\n    elsif ( $#array<4){\
 return ($pdb_id, \"\");}\n    else {return ( $pdb\
_id, $chain);}\n      \n    \n    \n}\nsub get_mai\
n_hom_code \n{\n    my $pdb_file=@_[0];\n    my %h\
om, $n, $best, $best_h;\n    open (F, $pdb_file);\\
n    while (<F>)\n{\n	if ( $_=~/^ATOM/)\n{\n	    $\
h=substr ($_,26, 1);\n	    $n=++$hom{$h};\n	    if\
 ($n>$best)\n{\n		$best=$n;\n		$best_h=$h;\n}\n}\n\
}\n    close (F);\n    return $best_h;\n}\n\n\nsub\
 get_pdb_file \n{\n    my ($pdb_file_in)=(@_);\n  \
  my $result;\n    my @letter;\n    my @chain;\n  \
  my $v;\n    my $pdb_file=$pdb_file_in;\n\n    $p\
db_file=($pdb_file_in=~/\\S+_S_(\\S+)/)?$1:$pdb_fi\
le_in;\n    \n    if ($no_remote_pdb_dir==0)\n    \
  {\n	$no_remote_pdb_dir=1;\n	$result=get_pdb_file\
3 ($pdb_file);\n	$no_remote_pdb_dir=0;\n	if ( $res\
ult){return $result;}\n	else\n	  {\n	    \n	    lc\
 ($pdb_file);\n	    $result=get_pdb_file3($pdb_fil\
e);\n	    return  $result;\n	  }\n      }\n    els\
e\n      {\n	return get_pdb_file3 ($pdb_file);\n  \
    }\n    \n  }\n\nsub get_pdb_file3 \n{\n    my \
$pdb_file_in=@_[0];\n    my $result;\n    my @lett\
er;\n    my @chain;\n    my $lcfile;\n    my $ucfi\
le;\n    my $pdb_file=$pdb_file_in;\n    \n    $lc\
file=lc $pdb_file;\n    $ucfile=uc $pdb_file;\n\n \
   if ( ($result=get_pdb_file2 ($pdb_file))){retur\
n $result;}\n    \n\n    if ($lcfile ne $pdb_file \
&& ($result=get_pdb_file2 ($lcfile))){return $resu\
lt;}\n    if ($ucfile ne $pdb_file && ($result=get\
_pdb_file2 ($ucfile))){return $result;}\n    \n   \
\n    \n    return \"\";\n}\nsub get_pdb_file2\n{\\
n    my $pdb_file=@_[0];\n    my $return_value;\n \
   \n    $return_value=\"\";\n    \n    if ( ($res\
ult=get_pdb_file1 ($pdb_file))){$return_value=$res\
ult;}\n    elsif ( !($pdb_file=~/\\.pdb/) && !($pd\
b_file=~/\\.PDB/))\n{\n	if ( ($result=get_pdb_file\
1 (\"$pdb_file.pdb\"))){$return_value=$result;}\n	\
elsif ( ($result=get_pdb_file1 (\"$pdb_file.PDB\")\
)){$return_value=$result;}\n\n	elsif ( ($result=ge\
t_pdb_file1 (\"pdb$pdb_file.pdb\"))){$return_value\
=$result;}	\n	elsif ( ($result=get_pdb_file1 (\"pd\
b$pdb_file.PDB\"))){$return_value=$result;}\n	elsi\
f ( ($result=get_pdb_file1 (\"PDB$pdb_file.PDB\"))\
){$return_value=$result;}\n	elsif ( ($result=get_p\
db_file1 (\"PDB$pdb_file.pdb\"))){$return_value=$r\
esult;}\n	\n	\n	elsif ( ($result=get_pdb_file1 (\"\
$pdb_file.ent\"))){$return_value=$result;}\n	elsif\
 ( ($result=get_pdb_file1 (\"pdb$pdb_file.ent\")))\
{$return_value=$result;}\n	elsif ( ($result=get_pd\
b_file1 (\"PDB$pdb_file.ent\"))){$return_value=$re\
sult;}\n\n	elsif ( ($result=get_pdb_file1 (\"$pdb_\
file.ENT\"))){$return_value=$result;}\n	elsif ( ($\
result=get_pdb_file1 (\"pdb$pdb_file.ENT\"))){$ret\
urn_value=$result;}\n	elsif ( ($result=get_pdb_fil\
e1 (\"PDB$pdb_file.ENT\"))){$return_value=$result;\
}\n	\n	\n	\n}\n    return $return_value;\n}\n    \\
nsub get_pdb_file1\n{\n    my ($pdb_file)=(@_);\n \
   my $return_value;\n    \n\n    $return_value=\"\
\";\n    if ( ($result=get_pdb_file0 ($pdb_file)))\
{$return_value=$result;}\n    elsif ( ($result=get\
_pdb_file0 (\"$pdb_file.Z\"))){$return_value=$resu\
lt;}\n    elsif ( ($result=get_pdb_file0 (\"$pdb_f\
ile.gz\"))){$return_value=$result;}\n    elsif ( (\
$result=get_pdb_file0 (\"$pdb_file.GZ\"))){$return\
_value=$result;}\n    return $return_value;\n}\nsu\
b get_pdb_file0 \n{ \n    my ($pdb_file, $attempt)\
=(@_);\n    my $pdb_file=@_[0];\n    my $tmp_pdb_f\
ile;    \n    my $return_value;\n\n    if ( !$atte\
mpt){$attempt=1;}\n    \n    $local_pdb_file=\"$pd\
b_file\";\n    if ( $local_pdb_file eq \"\")\n{\n	\
$tmp_pdb_file=vtmpnam();\n	open F, \">$tmp_pdb_fil\
e\";\n	\n	while (<STDIN>){print F \"$_\";}\n	close\
 (F);\n	\n	if (-e $tmp_pdb_file && &is_pdb_file ( \
$local_pdb_file))\n{return $tmp_pdb_file;}\n}\n\n \
   $local_pdb_file=\"$pdb_file\";\n    &debug_prin\
t (\"\\nTry access local file: $local_pdb_file\");\
\n    \n    $local_pdb_file=&check_pdb_file4compre\
ssion ($local_pdb_file);\n    if ( -e $local_pdb_f\
ile && (&is_pdb_file ($local_pdb_file) || $force_p\
db))\n{\n	&debug_print ( \"\\n\\tIs in Current Dir\
\");\n	$tmp_pdb_file=vtmpnam();\n	`cp $local_pdb_f\
ile $tmp_pdb_file`;\n	return $tmp_pdb_file;\n}\n  \
  else\n{\n	&debug_print (\"\\n\\tFile Not in Curr\
ent Dir\");\n}\n\n    if ($pdb_file=~/^pdb/||$pdb_\
file=~/^PDB/){$pdb_div=substr ($pdb_file, 4, 2);}\\
n    else\n{\n	  $pdb_div=substr ($pdb_file, 1, 2)\
;\n}\n    $local_pdb_file=\"$pdb_dir/$pdb_div/$pdb\
_file\";\n    $local_pdb_file=&check_pdb_file4comp\
ression ( $local_pdb_file);\n    &debug_print (\"\\
\nTry access file From PDB_DIR: $local_pdb_file\")\
;\n    if ($pdb_dir && -e $local_pdb_file && &is_p\
db_file ($local_pdb_file))\n{\n	&debug_print ( \"\\
\n\\tIs in Local PDB DIR\");\n	$tmp_pdb_file=vtmpn\
am();\n	`cp $local_pdb_file $tmp_pdb_file`;\n	retu\
rn $tmp_pdb_file;\n}\n\n    $local_pdb_file=\"$pdb\
_dir/$pdb_file\";\n    $local_pdb_file=&check_pdb_\
file4compression ( $local_pdb_file);\n    &debug_p\
rint (\"\\nTry access file From PDB_DIR: local_pdb\
_file\");\n    if ($pdb_dir && -e $local_pdb_file \
&& &is_pdb_file ($local_pdb_file))\n{\n	&debug_pri\
nt ( \"\\n\\tIs in Local PDB DIR\");\n	$tmp_pdb_fi\
le=vtmpnam();\n	`cp $local_pdb_file $tmp_pdb_file`\
;\n	return $tmp_pdb_file;\n}\n\n    $local_pdb_fil\
e=\"$pdb_dir$pdb_file\";\n    $local_pdb_file=&che\
ck_pdb_file4compression ( $local_pdb_file);\n    &\
debug_print (\"\\nTry access file From PDB_DIR: $l\
ocal_pdb_file\");\n    if ($pdb_dir && -e $local_p\
db_file && &is_pdb_file ($local_pdb_file))\n{\n	&d\
ebug_print ( \"\\n\\tIs in Local PDB DIR\");\n	$tm\
p_pdb_file=vtmpnam();\n	`cp $local_pdb_file $tmp_p\
db_file`;\n	return $tmp_pdb_file;\n}\n    else\n{&\
debug_print ( \"\\n\\tNot In Local Pdb Dir\");}\n\\
n    if ($cache ne \"NO\" && $cache ne \"no\")\n{\\
n\n	$local_pdb_file=\"$cache/$pdb_file\";\n	$local\
_pdb_file=&check_pdb_file4compression ( $local_pdb\
_file);\n	&debug_print(\"\\nTry access file From C\
ache: $local_pdb_file\");\n	if (-e $local_pdb_file\
 && &is_pdb_file ($local_pdb_file))\n{\n	    &debu\
g_print ( \"\\n\\tIs in T-Coffee Cache\");\n	    $\
tmp_pdb_file=vtmpnam();\n	    `cp $local_pdb_file \
$tmp_pdb_file`;\n	    return $tmp_pdb_file;\n}\n	e\
lse{&debug_print ( \"\\n\\tNot in Cache Dir\");}\n\
}\n\nif (!$no_remote_pdb_dir) \n  {\n    my $value\
=&is_released ($pdb_file);\n    my $return_value=\\
"\";\n    if ($value==1)\n      {\n	\n	&debug_prin\
t (\"\\n******************************************\
***********\\nTry Remote Access for $pdb_file\");\\
n	$tmp_pdb_file=vtmpnam();\n	$netcommand=$netaddre\
ss;\n	$netcommand=~s/%%/$pdb_file/g;\n	&url2file(\\
"$netcommand\", \"$tmp_pdb_file.$netcompression\")\
;\n	&debug_print(\"\\nREMOTE: $netcommand\\n\");\n\
	\n	$compressed_tmp_file_name=\"$tmp_pdb_file.$net\
compression\";\n	\n	if ($netcompression && -B $com\
pressed_tmp_file_name)\n	  {\n	    my $r;\n	    &d\
ebug_print (\"\\n\\tFile Found Remotely\");\n	    \
if (($r=safe_system ( \"$netcompression_pg $compre\
ssed_tmp_file_name\")!=$EXIT_SUCCESS) && $attempts\
<5)\n	      {\n		&debug_print (\"\\n\\tProper Down\
load Failed Try again\");\n		unlink $compressed_tm\
p_file_name;\n		print \"\\nFailed to Download $com\
pressed_tmp_file_name. New Attempt $attempt/5\\n\"\
;\n		return &get_pdb_file0($pdb_file, $attempt+1);\
\n	      }\n	    elsif ($r== $EXIT_SUCCESS)\n	    \
  {\n		&debug_print (\"\\n\\tProper Download Succe\
eded \");\n		$return_value=$tmp_pdb_file;\n	      \
}\n	    else\n	      {\n		&debug_print (\"\\n\\tPr\
oper Download Failed \");\n		&debug_print (\"\\nFi\
le Not Found Remotely\");\n		unlink $compressed_tm\
p_file_name;\n	      }\n	  }\n	else\n	  {\n\n	    \
&debug_print (\"\\nFile Not Found Remotely\");\n	 \
   unlink $compressed_tmp_file_name;\n	  }\n	#Upda\
te cache if required\n	if ($cache ne \"no\" && $ca\
che ne \"update\" && -e $return_value)\n	  {\n	   \
 `cp $return_value $cache/$pdb_file.pdb`;\n	    #`\
t_coffee -other_pg clean_cache.pl -file $pdb_file.\
pdb -dir $cache`;\n	  }\n      }\n    &debug_print\
 (\"\\nRemote Download Finished\");\n    return $r\
eturn_value;\n  }\nreturn \"\";\n}\n\nsub check_pd\
b_file4compression \n{\n    my $file=@_[0];\n    m\
y $tmp;\n    my $r;\n    \n    $tmp=&vtmpnam();\n \
   if (-e $tmp){unlink $tmp;}\n    \n    $file=~s/\
\\/\\//\\//g;\n    if    (-B $file && ($file=~/\\.\
Z/)) {`cp $file $tmp.Z`;`rm $tmp`;`gunzip $tmp.Z $\
SILENT`;$r=$tmp;}\n    elsif (-B $file && ($file=~\
/\\.gz/)){`cp $file $tmp.gz`;`gunzip $tmp.gz $SILE\
NT`;return $r=$tmp;}\n    elsif (-B $file ){`cp $f\
ile $tmp.gz`;`gunzip $tmp.gz $SILENT`;$r=$tmp;}\n \
   elsif ( -e $file ) {$r= $file;}\n    elsif ( -e\
 \"$file.gz\" ){ `cp $file.gz $tmp.gz`;`gunzip    \
 $tmp.gz $SILENT`;$r=$tmp;}    \n    elsif ( -e \"\
$file.Z\") {`cp $file.Z  $tmp.Z`; `gunzip $tmp.Z $\
SILENT`;$r=$tmp;}\n    else  {$r= $file;}\n\n    i\
f ( -e \"$tmp.Z\"){unlink \"$tmp.Z\";}\n    if ( -\
e \"$tmp.gz\"){unlink \"$tmp.gz\";}\n    \n    ret\
urn $r;\n    \n}\n\n\n\n\n\n    \n\n\n\n\n\n\n\nsu\
b vfopen \n{\n    my $file=@_[0];\n    my $mode=@_\
[1];\n    my $tmp;\n    my $F = new FileHandle;\n \
   \n    \n    $tmp=$file;\n	\n    \n    if ( $mod\
e eq \"r\" && !-e $file){ myexit(flush_error (\"Ca\
nnot open file $file\"));}\n    elsif ($mode eq \"\
w\"){$tmp=\">$file\";}\n    elsif ($mode eq \"a\")\
{$tmp=\">>$file\";}\n    \n    \n    open ($F,$tmp\
);\n    return $F;\n}\nsub debug_print\n{\n    my \
$message =@_[0];\n    if ($debug){print STDERR \"N\
O_REMOTE_PDB_DIR: $no_remote_pdb_dir - $message [D\
EBUG:extract_from_pdb]\";}\n    return;\n}\nsub is\
_aa \n{\n    my ($aa, $chain) =@_;\n\n    my $one;\
\n    my $trhee;\n    \n    if ( $onelett{$molecul\
e_type{$chain}}->{$aa} eq 'X' || !$onelett{$molecu\
le_type{$chain}}->{$aa} ){return '';}\n    else\n \
     {\n	$one=$onelett{$molecule_type{$chain}}->{$\
aa};\n\n	$three=$threelett{$molecule_type{$chain}}\
->{$one};\n	\n\n	return $three;\n      }\n  }\n\n\\
n\n\n\nsub url2file\n{\n    my ($address, $out, $w\
get_arg, $curl_arg)=(@_);\n    my ($pg, $flag, $r,\
 $arg, $count);\n    \n    if (!$CONFIGURATION){&c\
heck_configuration (\"wget\", \"INTERNET\", \"gzip\
\");$CONFIGURATION=1;}\n    \n    if (&pg_is_insta\
lled (\"wget\"))   {$pg=\"wget\"; $flag=\"-O\";$ar\
g=$wget_arg;}\n    elsif (&pg_is_installed (\"curl\
\")){$pg=\"curl\"; $flag=\"-o\";$arg=$curl_arg;}\n\
    return safe_system (\"$pg $flag$out $address >\
/dev/null 2>/dev/null\");\n\n}\n\n\n\n\nsub pdbfil\
e2chaintype\n  {\n    my $file=@_[0];\n    my %ct;\
\n    my $F;\n    \n    $F=vfopen ($file, \"r\");\\
n    while (<$F>)\n      {\n	my $line=$_;\n	if ($l\
ine =~/^ATOM/)\n	  {\n	    my $C=substr($line,21,1\
);\n	    if (!$ct{$C})\n	      {	\n		my $r=substr(\
$line,17,3);\n		$r=~s/\\s+//;\n		if (length ($r)==\
1){$ct{$C}=\"R\";}\n		elsif (length ($r)==2){$ct{$\
C}=\"D\";}\n		elsif (length ($r)==3){$ct{$C}=\"P\"\
;}\n		else \n		  {\n		    myexit(flush_error(\"ERR\
OR: Could not read RES_ID field in file $file\"));\
\n		  }\n	      }\n	  }\n      }\n    close ($F);\\
n    return %ct;\n  }\n   \n    \n\n\n\nsub fill_t\
hreelett_RNA\n{\n\n	my %threelett=(\n	'A', '  A',\\
n	'T', '  T',\n	'U', '  U',\n	'C', '  C',\n	'G', '\
  G',\n	'I', '  I', #Inosine\n	);\n	\n	return %thr\
eelett;\n\n}\n\n\nsub fill_onelett_RNA\n{\n	my   %\
onelett=(\n	'  A' => 'A',\n	'  T' => 'T',\n	'  U' \
=> 'U',\n	'  C' => 'C',\n	'  G' => 'G',\n	'CSL' =>\
 'X',\n	'UMS' => 'X',\n	'  I' => 'I',\n	'A' => 'A'\
,\n	'T' => 'T',\n	'U' => 'U',\n	'C' => 'C',\n	'G' \
=> 'G',\n	'I' => 'I',\n	);\n\n	return %onelett;\n\\
n}\n\n\nsub fill_onelett_DNA\n{\n	my   %onelett=(\\
n	' DA', 'A',\n	' DT', 'T',\n	' DC', 'C',\n	' DG',\
 'G',\n	'DA', 'A',\n	'DT', 'T',\n	'DC', 'C',\n	'DG\
', 'G',\n	);\n\n	return %onelett;\n\n}\n\nsub fill\
_threelett_DNA\n{\n\n	my %threelett=(\n	'A', ' DA'\
,\n	'T', ' DT',\n	'C', ' DC',\n	'G', ' DG',\n	);\n\
\n	return %threelett;\n\n}\n\n\n\n\nsub fill_three\
lett_prot\n{  \n  my %threelett;\n\n  %threelett=(\
\n'A', 'ALA',\n'C', 'CYS',\n'D', 'ASP',\n'E', 'GLU\
',\n'F', 'PHE',\n'G', 'GLY',\n'H', 'HIS',\n'I', 'I\
LE',\n'K', 'LYS',\n'L', 'LEU',\n'N', 'ASN',\n'M', \
'MET',\n'P', 'PRO',\n'Q', 'GLN',\n'R', 'ARG',\n'S'\
, 'SER',\n'T', 'THR',\n'V', 'VAL',\n'W', 'TRP',\n'\
Y', 'TYR',\n);\n\nreturn %threelett;\n\n\n}\n\nsub\
 fill_onelett_prot\n{\n    my %onelett;\n    \n   \
 %onelett=(\n\n'10A', 'X',\n'11O', 'X',\n'12A', 'X\
',\n'13P', 'X',\n'13R', 'X',\n'13S', 'X',\n'14W', \
'X',\n'15P', 'X',\n'16A', 'X',\n'16G', 'X',\n'1AN'\
, 'X',\n'1AP', 'X',\n'1AR', 'X',\n'1BH', 'X',\n'1B\
O', 'X',\n'1C5', 'X',\n'1CU', 'X',\n'1DA', 'X',\n'\
1GL', 'X',\n'1GN', 'X',\n'1IN', 'X',\n'1LU', 'L',\\
n'1MA', 'X',\n'1MC', 'X',\n'1MG', 'X',\n'1MZ', 'X'\
,\n'1NA', 'X',\n'1NB', 'X',\n'1NI', 'X',\n'1PA', '\
A',\n'1PC', 'X',\n'1PE', 'X',\n'1PG', 'X',\n'1PI',\
 'A',\n'1PM', 'X',\n'1PN', 'X',\n'1PU', 'X',\n'1PY\
', 'X',\n'1UN', 'X',\n'24T', 'X',\n'25T', 'X',\n'2\
6P', 'X',\n'2AB', 'X',\n'2AM', 'X',\n'2AN', 'X',\n\
'2AP', 'X',\n'2AR', 'X',\n'2AS', 'D',\n'2BL', 'X',\
\n'2BM', 'X',\n'2CP', 'X',\n'2DA', 'X',\n'2DG', 'X\
',\n'2DP', 'X',\n'2DT', 'X',\n'2EP', 'X',\n'2EZ', \
'X',\n'2FG', 'X',\n'2FL', 'X',\n'2FP', 'X',\n'2FU'\
, 'X',\n'2GL', 'X',\n'2GP', 'X',\n'2HP', 'X',\n'2I\
B', 'X',\n'2IP', 'X',\n'2LU', 'L',\n'2MA', 'X',\n'\
2MD', 'X',\n'2ME', 'X',\n'2MG', 'X',\n'2ML', 'L',\\
n'2MO', 'X',\n'2MR', 'R',\n'2MU', 'X',\n'2MZ', 'X'\
,\n'2NO', 'X',\n'2NP', 'X',\n'2OG', 'X',\n'2PA', '\
X',\n'2PC', 'X',\n'2PE', 'X',\n'2PG', 'X',\n'2PH',\
 'X',\n'2PI', 'X',\n'2PL', 'X',\n'2PP', 'X',\n'2PU\
', 'X',\n'2SI', 'X',\n'2TB', 'X',\n'34C', 'X',\n'3\
5G', 'X',\n'3AA', 'X',\n'3AD', 'X',\n'3AH', 'H',\n\
'3AN', 'X',\n'3AP', 'X',\n'3AT', 'X',\n'3BT', 'X',\
\n'3CH', 'X',\n'3CN', 'X',\n'3CO', 'X',\n'3CP', 'X\
',\n'3DR', 'X',\n'3EP', 'X',\n'3FM', 'X',\n'3GA', \
'X',\n'3GP', 'X',\n'3HB', 'X',\n'3HC', 'X',\n'3HP'\
, 'X',\n'3IB', 'X',\n'3ID', 'X',\n'3IN', 'X',\n'3M\
A', 'X',\n'3MB', 'X',\n'3MC', 'X',\n'3MD', 'D',\n'\
3MF', 'X',\n'3MP', 'X',\n'3MT', 'X',\n'3OL', 'X',\\
n'3PA', 'X',\n'3PG', 'X',\n'3PO', 'X',\n'3PP', 'X'\
,\n'3PY', 'X',\n'49A', 'X',\n'4AB', 'X',\n'4AM', '\
X',\n'4AN', 'X',\n'4AP', 'X',\n'4BA', 'X',\n'4BT',\
 'X',\n'4CA', 'X',\n'4CO', 'X',\n'4HP', 'X',\n'4IP\
', 'X',\n'4MO', 'X',\n'4MV', 'X',\n'4MZ', 'X',\n'4\
NC', 'X',\n'4NP', 'X',\n'4OX', 'X',\n'4PB', 'X',\n\
'4PN', 'X',\n'4PP', 'X',\n'4SC', 'X',\n'4SU', 'X',\
\n'4TB', 'X',\n'55C', 'X',\n'5AD', 'X',\n'5AN', 'X\
',\n'5AT', 'X',\n'5CM', 'X',\n'5GP', 'X',\n'5HP', \
'E',\n'5HT', 'X',\n'5IT', 'X',\n'5IU', 'X',\n'5MB'\
, 'X',\n'5MC', 'X',\n'5MD', 'X',\n'5MP', 'X',\n'5M\
U', 'X',\n'5NC', 'X',\n'5OB', 'X',\n'5PA', 'X',\n'\
5PV', 'X',\n'6AB', 'X',\n'6CT', 'X',\n'6HA', 'X',\\
n'6HC', 'X',\n'6HG', 'X',\n'6HT', 'X',\n'6IN', 'X'\
,\n'6MO', 'X',\n'6MP', 'X',\n'6PG', 'X',\n'6WO', '\
X',\n'70U', 'X',\n'7DG', 'X',\n'7HP', 'X',\n'7I2',\
 'X',\n'7MG', 'X',\n'7MQ', 'X',\n'7NI', 'X',\n'87Y\
', 'X',\n'8AD', 'X',\n'8BR', 'X',\n'8IG', 'X',\n'8\
IN', 'X',\n'8OG', 'X',\n'95A', 'X',\n'9AD', 'X',\n\
'9AM', 'X',\n'9AP', 'X',\n'9DG', 'X',\n'9DI', 'X',\
\n'9HX', 'X',\n'9OH', 'X',\n'9TA', 'X',\n'A12', 'X\
',\n'A15', 'X',\n'A23', 'X',\n'A24', 'X',\n'A26', \
'X',\n'A2G', 'X',\n'A2P', 'X',\n'A32', 'X',\n'A3P'\
, 'X',\n'A4P', 'X',\n'A5P', 'X',\n'A70', 'X',\n'A7\
6', 'X',\n'A77', 'X',\n'A78', 'X',\n'A79', 'X',\n'\
A80', 'X',\n'A85', 'X',\n'A88', 'X',\n'A9A', 'X',\\
n'AA3', 'X',\n'AA4', 'X',\n'AA6', 'X',\n'AAA', 'X'\
,\n'AAB', 'X',\n'AAC', 'X',\n'AAE', 'X',\n'AAG', '\
R',\n'AAH', 'X',\n'AAM', 'X',\n'AAN', 'X',\n'AAP',\
 'X',\n'AAR', 'R',\n'AAS', 'X',\n'AAT', 'X',\n'ABA\
', 'X',\n'ABC', 'X',\n'ABD', 'X',\n'ABE', 'X',\n'A\
BH', 'X',\n'ABI', 'X',\n'ABK', 'X',\n'ABM', 'X',\n\
'ABN', 'X',\n'ABP', 'X',\n'ABR', 'X',\n'ABS', 'X',\
\n'ABU', 'X',\n'AC1', 'X',\n'AC2', 'X',\n'ACA', 'X\
',\n'ACB', 'D',\n'ACC', 'C',\n'ACD', 'X',\n'ACE', \
'X',\n'ACH', 'X',\n'ACI', 'X',\n'ACL', 'R',\n'ACM'\
, 'X',\n'ACN', 'X',\n'ACO', 'X',\n'ACP', 'X',\n'AC\
Q', 'X',\n'ACR', 'X',\n'ACS', 'X',\n'ACT', 'X',\n'\
ACV', 'V',\n'ACX', 'X',\n'ACY', 'X',\n'AD2', 'X',\\
n'AD3', 'X',\n'ADC', 'X',\n'ADD', 'X',\n'ADE', 'X'\
,\n'ADH', 'X',\n'ADI', 'X',\n'ADM', 'X',\n'ADN', '\
X',\n'ADP', 'X',\n'ADQ', 'X',\n'ADR', 'X',\n'ADS',\
 'X',\n'ADT', 'X',\n'ADU', 'X',\n'ADW', 'X',\n'ADX\
', 'X',\n'AE2', 'X',\n'AEA', 'X',\n'AEB', 'X',\n'A\
EI', 'D',\n'AEN', 'X',\n'AET', 'T',\n'AF1', 'X',\n\
'AF3', 'X',\n'AFA', 'D',\n'AFP', 'X',\n'AG7', 'X',\
\n'AGB', 'X',\n'AGF', 'X',\n'AGL', 'X',\n'AGM', 'R\
',\n'AGN', 'X',\n'AGP', 'X',\n'AGS', 'X',\n'AGU', \
'X',\n'AH0', 'X',\n'AH1', 'X',\n'AHA', 'X',\n'AHB'\
, 'D',\n'AHC', 'X',\n'AHF', 'X',\n'AHG', 'X',\n'AH\
H', 'X',\n'AHM', 'X',\n'AHO', 'X',\n'AHP', 'X',\n'\
AHS', 'X',\n'AHT', 'Y',\n'AHU', 'X',\n'AHX', 'X',\\
n'AI1', 'X',\n'AI2', 'X',\n'AIB', 'X',\n'AIC', 'X'\
,\n'AIM', 'X',\n'AIP', 'X',\n'AIQ', 'X',\n'AIR', '\
X',\n'AJ3', 'X',\n'AKB', 'X',\n'AKG', 'X',\n'AKR',\
 'X',\n'AL1', 'X',\n'AL2', 'X',\n'AL3', 'X',\n'AL4\
', 'X',\n'AL5', 'X',\n'AL6', 'X',\n'AL7', 'X',\n'A\
L8', 'X',\n'AL9', 'X',\n'ALA', 'A',\n'ALB', 'X',\n\
'ALC', 'X',\n'ALD', 'L',\n'ALE', 'X',\n'ALF', 'X',\
\n'ALG', 'X',\n'ALL', 'X',\n'ALM', 'A',\n'ALN', 'A\
',\n'ALO', 'T',\n'ALP', 'X',\n'ALQ', 'X',\n'ALR', \
'X',\n'ALS', 'X',\n'ALT', 'A',\n'ALY', 'K',\n'ALZ'\
, 'X',\n'AMA', 'X',\n'AMB', 'X',\n'AMC', 'X',\n'AM\
D', 'X',\n'AMG', 'X',\n'AMH', 'X',\n'AMI', 'X',\n'\
AML', 'X',\n'AMN', 'X',\n'AMO', 'X',\n'AMP', 'X',\\
n'AMQ', 'X',\n'AMR', 'X',\n'AMS', 'X',\n'AMT', 'X'\
,\n'AMU', 'X',\n'AMW', 'X',\n'AMX', 'X',\n'AMY', '\
X',\n'ANA', 'X',\n'ANB', 'X',\n'ANC', 'X',\n'AND',\
 'X',\n'ANE', 'X',\n'ANI', 'X',\n'ANL', 'X',\n'ANO\
', 'X',\n'ANP', 'X',\n'ANS', 'X',\n'ANT', 'X',\n'A\
OE', 'X',\n'AOP', 'X',\n'AP1', 'X',\n'AP2', 'X',\n\
'AP3', 'X',\n'AP4', 'X',\n'AP5', 'X',\n'AP6', 'X',\
\n'APA', 'X',\n'APB', 'X',\n'APC', 'X',\n'APE', 'F\
',\n'APF', 'X',\n'APG', 'X',\n'APH', 'A',\n'API', \
'X',\n'APL', 'X',\n'APM', 'X',\n'APN', 'G',\n'APP'\
, 'X',\n'APQ', 'X',\n'APR', 'X',\n'APS', 'X',\n'AP\
T', 'X',\n'APU', 'X',\n'APX', 'X',\n'APY', 'X',\n'\
APZ', 'X',\n'AQS', 'X',\n'AR1', 'X',\n'AR2', 'X',\\
n'ARA', 'X',\n'ARB', 'X',\n'ARC', 'X',\n'ARD', 'X'\
,\n'ARG', 'R',\n'ARH', 'X',\n'ARI', 'X',\n'ARM', '\
R',\n'ARN', 'X',\n'ARO', 'R',\n'ARP', 'X',\n'ARQ',\
 'X',\n'ARS', 'X',\n'AS1', 'R',\n'AS2', 'X',\n'ASA\
', 'D',\n'ASB', 'D',\n'ASC', 'X',\n'ASD', 'X',\n'A\
SE', 'X',\n'ASF', 'X',\n'ASI', 'X',\n'ASK', 'D',\n\
'ASL', 'X',\n'ASM', 'N',\n'ASO', 'X',\n'ASP', 'D',\
\n'ASQ', 'X',\n'ASU', 'X',\n'ATA', 'X',\n'ATC', 'X\
',\n'ATD', 'X',\n'ATF', 'X',\n'ATG', 'X',\n'ATH', \
'X',\n'ATM', 'X',\n'ATO', 'X',\n'ATP', 'X',\n'ATQ'\
, 'X',\n'ATR', 'X',\n'ATT', 'X',\n'ATY', 'X',\n'AT\
Z', 'X',\n'AUC', 'X',\n'AUR', 'X',\n'AVG', 'X',\n'\
AXP', 'X',\n'AYA', 'A',\n'AZ2', 'X',\n'AZA', 'X',\\
n'AZC', 'X',\n'AZD', 'X',\n'AZE', 'X',\n'AZI', 'X'\
,\n'AZL', 'X',\n'AZM', 'X',\n'AZR', 'X',\n'AZT', '\
X',\n'B12', 'X',\n'B1F', 'F',\n'B2A', 'A',\n'B2F',\
 'F',\n'B2I', 'I',\n'B2V', 'V',\n'B3I', 'X',\n'B3P\
', 'X',\n'B7G', 'X',\n'B96', 'X',\n'B9A', 'X',\n'B\
A1', 'X',\n'BAA', 'X',\n'BAB', 'X',\n'BAC', 'X',\n\
'BAF', 'X',\n'BAH', 'X',\n'BAI', 'X',\n'BAK', 'X',\
\n'BAL', 'A',\n'BAM', 'X',\n'BAO', 'X',\n'BAP', 'X\
',\n'BAR', 'X',\n'BAS', 'X',\n'BAT', 'F',\n'BAY', \
'X',\n'BAZ', 'X',\n'BB1', 'X',\n'BB2', 'X',\n'BBA'\
, 'X',\n'BBH', 'X',\n'BBS', 'X',\n'BBT', 'X',\n'BB\
Z', 'X',\n'BCA', 'X',\n'BCB', 'X',\n'BCC', 'X',\n'\
BCD', 'X',\n'BCL', 'X',\n'BCN', 'X',\n'BCR', 'X',\\
n'BCS', 'C',\n'BCT', 'X',\n'BCY', 'X',\n'BCZ', 'X'\
,\n'BDA', 'X',\n'BDG', 'X',\n'BDK', 'X',\n'BDM', '\
X',\n'BDN', 'X',\n'BDS', 'X',\n'BE1', 'X',\n'BE2',\
 'X',\n'BEA', 'X',\n'BEF', 'X',\n'BEN', 'X',\n'BEO\
', 'X',\n'BEP', 'X',\n'BER', 'X',\n'BES', 'X',\n'B\
ET', 'X',\n'BEZ', 'X',\n'BF2', 'X',\n'BFA', 'X',\n\
'BFD', 'X',\n'BFP', 'X',\n'BFS', 'X',\n'BFU', 'X',\
\n'BG6', 'X',\n'BGF', 'X',\n'BGG', 'X',\n'BGL', 'X\
',\n'BGN', 'X',\n'BGP', 'X',\n'BGX', 'X',\n'BH4', \
'X',\n'BHA', 'X',\n'BHC', 'X',\n'BHD', 'D',\n'BHO'\
, 'X',\n'BHS', 'X',\n'BIC', 'X',\n'BIN', 'X',\n'BI\
O', 'X',\n'BIP', 'X',\n'BIS', 'X',\n'BIZ', 'X',\n'\
BJH', 'X',\n'BJI', 'X',\n'BJP', 'X',\n'BLA', 'X',\\
n'BLB', 'X',\n'BLE', 'L',\n'BLG', 'P',\n'BLI', 'X'\
,\n'BLM', 'X',\n'BLV', 'X',\n'BLY', 'K',\n'BM1', '\
X',\n'BM2', 'X',\n'BM5', 'X',\n'BM9', 'X',\n'BMA',\
 'X',\n'BMD', 'X',\n'BME', 'X',\n'BMP', 'X',\n'BMQ\
', 'X',\n'BMS', 'X',\n'BMT', 'T',\n'BMU', 'X',\n'B\
MY', 'X',\n'BMZ', 'X',\n'BNA', 'X',\n'BNG', 'X',\n\
'BNI', 'X',\n'BNN', 'F',\n'BNO', 'L',\n'BNS', 'X',\
\n'BNZ', 'X',\n'BO3', 'X',\n'BO4', 'X',\n'BOC', 'X\
',\n'BOG', 'X',\n'BOM', 'X',\n'BOT', 'X',\n'BOX', \
'X',\n'BOZ', 'X',\n'BPA', 'X',\n'BPB', 'X',\n'BPD'\
, 'X',\n'BPG', 'X',\n'BPH', 'X',\n'BPI', 'X',\n'BP\
J', 'X',\n'BPM', 'X',\n'BPN', 'X',\n'BPO', 'X',\n'\
BPP', 'X',\n'BPT', 'X',\n'BPY', 'X',\n'BRB', 'X',\\
n'BRC', 'X',\n'BRE', 'X',\n'BRI', 'X',\n'BRL', 'X'\
,\n'BRM', 'X',\n'BRN', 'X',\n'BRO', 'X',\n'BRS', '\
X',\n'BRU', 'X',\n'BRZ', 'X',\n'BSB', 'X',\n'BSI',\
 'X',\n'BSP', 'X',\n'BT1', 'X',\n'BT2', 'X',\n'BT3\
', 'X',\n'BTA', 'L',\n'BTB', 'X',\n'BTC', 'C',\n'B\
TD', 'X',\n'BTN', 'X',\n'BTP', 'X',\n'BTR', 'W',\n\
'BU1', 'X',\n'BUA', 'X',\n'BUB', 'X',\n'BUC', 'X',\
\n'BUG', 'X',\n'BUL', 'X',\n'BUM', 'X',\n'BUQ', 'X\
',\n'BUT', 'X',\n'BVD', 'X',\n'BX3', 'X',\n'BYS', \
'X',\n'BZ1', 'X',\n'BZA', 'X',\n'BZB', 'X',\n'BZC'\
, 'X',\n'BZD', 'X',\n'BZF', 'X',\n'BZI', 'X',\n'BZ\
M', 'X',\n'BZO', 'X',\n'BZP', 'X',\n'BZQ', 'X',\n'\
BZS', 'X',\n'BZT', 'X',\n'C02', 'X',\n'C11', 'X',\\
n'C1O', 'X',\n'C20', 'X',\n'C24', 'X',\n'C2F', 'X'\
,\n'C2O', 'X',\n'C2P', 'X',\n'C3M', 'X',\n'C3P', '\
X',\n'C3X', 'X',\n'C48', 'X',\n'C4M', 'X',\n'C4X',\
 'X',\n'C5C', 'X',\n'C5M', 'X',\n'C5P', 'X',\n'C5X\
', 'X',\n'C60', 'X',\n'C6C', 'X',\n'C6M', 'X',\n'C\
78', 'X',\n'C8E', 'X',\n'CA3', 'X',\n'CA5', 'X',\n\
'CAA', 'X',\n'CAB', 'X',\n'CAC', 'X',\n'CAD', 'X',\
\n'CAF', 'C',\n'CAG', 'X',\n'CAH', 'X',\n'CAL', 'X\
',\n'CAM', 'X',\n'CAN', 'X',\n'CAO', 'X',\n'CAP', \
'X',\n'CAQ', 'X',\n'CAR', 'X',\n'CAS', 'C',\n'CAT'\
, 'X',\n'CAV', 'X',\n'CAY', 'C',\n'CAZ', 'X',\n'CB\
3', 'X',\n'CB4', 'X',\n'CBA', 'X',\n'CBD', 'X',\n'\
CBG', 'X',\n'CBI', 'X',\n'CBL', 'X',\n'CBM', 'X',\\
n'CBN', 'X',\n'CBO', 'X',\n'CBP', 'X',\n'CBS', 'X'\
,\n'CBX', 'X',\n'CBZ', 'X',\n'CC0', 'X',\n'CC1', '\
X',\n'CCC', 'X',\n'CCH', 'X',\n'CCI', 'X',\n'CCM',\
 'X',\n'CCN', 'X',\n'CCO', 'X',\n'CCP', 'X',\n'CCR\
', 'X',\n'CCS', 'C',\n'CCV', 'X',\n'CCY', 'X',\n'C\
D1', 'X',\n'CDC', 'X',\n'CDE', 'X',\n'CDF', 'X',\n\
'CDI', 'X',\n'CDL', 'X',\n'CDM', 'X',\n'CDP', 'X',\
\n'CDR', 'X',\n'CDU', 'X',\n'CE1', 'X',\n'CEA', 'C\
',\n'CEB', 'X',\n'CEC', 'X',\n'CED', 'X',\n'CEF', \
'X',\n'CEH', 'X',\n'CEM', 'X',\n'CEO', 'X',\n'CEP'\
, 'X',\n'CEQ', 'X',\n'CER', 'X',\n'CES', 'G',\n'CE\
T', 'X',\n'CFC', 'X',\n'CFF', 'X',\n'CFM', 'X',\n'\
CFO', 'X',\n'CFP', 'X',\n'CFS', 'X',\n'CFX', 'X',\\
n'CGN', 'X',\n'CGP', 'X',\n'CGS', 'X',\n'CGU', 'E'\
,\n'CH2', 'X',\n'CH3', 'X',\n'CHA', 'X',\n'CHB', '\
X',\n'CHD', 'X',\n'CHF', 'X',\n'CHG', 'G',\n'CHI',\
 'X',\n'CHN', 'X',\n'CHO', 'X',\n'CHP', 'G',\n'CHR\
', 'X',\n'CHS', 'F',\n'CHT', 'X',\n'CHX', 'X',\n'C\
IC', 'X',\n'CIN', 'X',\n'CIP', 'X',\n'CIR', 'X',\n\
'CIT', 'X',\n'CIU', 'X',\n'CKI', 'X',\n'CL1', 'X',\
\n'CL2', 'X',\n'CLA', 'X',\n'CLB', 'A',\n'CLC', 'S\
',\n'CLD', 'A',\n'CLE', 'L',\n'CLF', 'X',\n'CLK', \
'S',\n'CLL', 'X',\n'CLM', 'X',\n'CLN', 'X',\n'CLO'\
, 'X',\n'CLP', 'X',\n'CLQ', 'X',\n'CLR', 'X',\n'CL\
S', 'X',\n'CLT', 'X',\n'CLX', 'X',\n'CLY', 'X',\n'\
CMA', 'R',\n'CMC', 'X',\n'CMD', 'X',\n'CME', 'C',\\
n'CMG', 'X',\n'CMK', 'X',\n'CMN', 'X',\n'CMO', 'X'\
,\n'CMP', 'X',\n'CMR', 'X',\n'CMS', 'X',\n'CMT', '\
C',\n'CMX', 'X',\n'CNA', 'X',\n'CNC', 'X',\n'CND',\
 'X',\n'CNH', 'X',\n'CNM', 'X',\n'CNN', 'X',\n'CNO\
', 'X',\n'CNP', 'X',\n'CO2', 'X',\n'CO3', 'X',\n'C\
O5', 'X',\n'CO8', 'X',\n'COA', 'X',\n'COB', 'X',\n\
'COC', 'X',\n'COD', 'X',\n'COE', 'X',\n'COF', 'X',\
\n'COH', 'X',\n'COI', 'X',\n'COJ', 'X',\n'COL', 'X\
',\n'COM', 'X',\n'CON', 'X',\n'COP', 'X',\n'COR', \
'X',\n'COS', 'X',\n'COT', 'X',\n'COY', 'X',\n'CP1'\
, 'G',\n'CP2', 'X',\n'CP4', 'X',\n'CPA', 'X',\n'CP\
B', 'X',\n'CPC', 'X',\n'CPD', 'X',\n'CPG', 'X',\n'\
CPH', 'X',\n'CPI', 'X',\n'CPM', 'X',\n'CPN', 'G',\\
n'CPO', 'X',\n'CPP', 'X',\n'CPQ', 'X',\n'CPR', 'X'\
,\n'CPS', 'X',\n'CPT', 'X',\n'CPU', 'X',\n'CPV', '\
X',\n'CPY', 'X',\n'CR1', 'X',\n'CR6', 'X',\n'CRA',\
 'X',\n'CRB', 'X',\n'CRC', 'X',\n'CRG', 'X',\n'CRH\
', 'X',\n'CRO', 'T',\n'CRP', 'X',\n'CRQ', 'X',\n'C\
RS', 'X',\n'CRT', 'X',\n'CRY', 'X',\n'CSA', 'C',\n\
'CSB', 'X',\n'CSD', 'C',\n'CSE', 'C',\n'CSH', 'X',\
\n'CSI', 'X',\n'CSN', 'X',\n'CSO', 'C',\n'CSP', 'C\
',\n'CSR', 'C',\n'CSS', 'C',\n'CST', 'X',\n'CSW', \
'C',\n'CSX', 'C',\n'CSY', 'X',\n'CSZ', 'C',\n'CT3'\
, 'X',\n'CTA', 'X',\n'CTB', 'X',\n'CTC', 'X',\n'CT\
D', 'X',\n'CTH', 'T',\n'CTO', 'X',\n'CTP', 'X',\n'\
CTR', 'X',\n'CTS', 'X',\n'CTT', 'X',\n'CTY', 'X',\\
n'CTZ', 'X',\n'CU1', 'X',\n'CUA', 'X',\n'CUC', 'X'\
,\n'CUL', 'X',\n'CUO', 'X',\n'CUZ', 'X',\n'CVI', '\
X',\n'CXF', 'X',\n'CXL', 'X',\n'CXM', 'M',\n'CXN',\
 'X',\n'CXP', 'X',\n'CXS', 'X',\n'CY1', 'C',\n'CY3\
', 'X',\n'CYB', 'X',\n'CYC', 'X',\n'CYF', 'C',\n'C\
YG', 'C',\n'CYH', 'X',\n'CYL', 'X',\n'CYM', 'C',\n\
'CYN', 'X',\n'CYO', 'X',\n'CYP', 'X',\n'CYQ', 'C',\
\n'CYS', 'C',\n'CYU', 'X',\n'CYY', 'X',\n'CYZ', 'X\
',\n'CZH', 'X',\n'CZZ', 'C',\n'D12', 'X',\n'D13', \
'X',\n'D16', 'X',\n'D18', 'X',\n'D19', 'X',\n'D1P'\
, 'X',\n'D24', 'X',\n'D34', 'X',\n'D35', 'X',\n'D4\
D', 'X',\n'D4T', 'X',\n'D6G', 'X',\n'DA2', 'R',\n'\
DA3', 'X',\n'DA6', 'X',\n'DA7', 'X',\n'DAA', 'X',\\
n'DAB', 'X',\n'DAC', 'X',\n'DAD', 'X',\n'DAE', 'X'\
,\n'DAF', 'X',\n'DAG', 'X',\n'DAH', 'A',\n'DAJ', '\
X',\n'DAK', 'X',\n'DAL', 'A',\n'DAM', 'A',\n'DAN',\
 'X',\n'DAO', 'X',\n'DAP', 'X',\n'DAQ', 'X',\n'DAR\
', 'R',\n'DAS', 'D',\n'DAT', 'X',\n'DAU', 'X',\n'D\
AV', 'X',\n'DBA', 'X',\n'DBD', 'X',\n'DBF', 'X',\n\
'DBG', 'X',\n'DBI', 'X',\n'DBV', 'X',\n'DBY', 'Y',\
\n'DCA', 'X',\n'DCB', 'X',\n'DCE', 'X',\n'DCF', 'X\
',\n'DCG', 'X',\n'DCH', 'X',\n'DCI', 'I',\n'DCL', \
'X',\n'DCM', 'X',\n'DCP', 'X',\n'DCS', 'X',\n'DCT'\
, 'X',\n'DCY', 'C',\n'DCZ', 'X',\n'DDA', 'X',\n'DD\
B', 'X',\n'DDC', 'X',\n'DDF', 'X',\n'DDG', 'X',\n'\
DDH', 'X',\n'DDL', 'X',\n'DDM', 'X',\n'DDO', 'L',\\
n'DDP', 'X',\n'DDQ', 'X',\n'DDT', 'Y',\n'DDU', 'X'\
,\n'DEA', 'X',\n'DEB', 'X',\n'DEC', 'X',\n'DEF', '\
X',\n'DEL', 'X',\n'DEM', 'X',\n'DEN', 'X',\n'DEP',\
 'X',\n'DEQ', 'X',\n'DES', 'X',\n'DET', 'X',\n'DFC\
', 'X',\n'DFG', 'X',\n'DFI', 'X',\n'DFL', 'X',\n'D\
FO', 'X',\n'DFP', 'X',\n'DFR', 'X',\n'DFT', 'X',\n\
'DFV', 'X',\n'DFX', 'X',\n'DG2', 'X',\n'DG3', 'X',\
\n'DG6', 'X',\n'DGA', 'X',\n'DGD', 'X',\n'DGG', 'X\
',\n'DGL', 'E',\n'DGN', 'Q',\n'DGP', 'X',\n'DGT', \
'X',\n'DGX', 'X',\n'DH2', 'X',\n'DHA', 'A',\n'DHB'\
, 'X',\n'DHC', 'X',\n'DHD', 'X',\n'DHE', 'X',\n'DH\
F', 'X',\n'DHG', 'X',\n'DHI', 'H',\n'DHL', 'X',\n'\
DHM', 'X',\n'DHN', 'V',\n'DHP', 'X',\n'DHQ', 'X',\\
n'DHR', 'X',\n'DHS', 'X',\n'DHT', 'X',\n'DHU', 'X'\
,\n'DHY', 'X',\n'DHZ', 'X',\n'DI2', 'X',\n'DI3', '\
G',\n'DI4', 'X',\n'DI5', 'X',\n'DIA', 'X',\n'DIC',\
 'X',\n'DIF', 'X',\n'DIG', 'X',\n'DII', 'X',\n'DIL\
', 'I',\n'DIM', 'X',\n'DIO', 'X',\n'DIP', 'X',\n'D\
IQ', 'X',\n'DIS', 'X',\n'DIT', 'X',\n'DIV', 'V',\n\
'DIX', 'X',\n'DIY', 'X',\n'DKA', 'X',\n'DLA', 'X',\
\n'DLE', 'L',\n'DLF', 'X',\n'DLS', 'K',\n'DLY', 'K\
',\n'DM1', 'X',\n'DM2', 'X',\n'DM3', 'X',\n'DM4', \
'X',\n'DM5', 'X',\n'DM6', 'X',\n'DM7', 'X',\n'DM8'\
, 'X',\n'DM9', 'X',\n'DMA', 'X',\n'DMB', 'X',\n'DM\
C', 'X',\n'DMD', 'X',\n'DME', 'X',\n'DMF', 'X',\n'\
DMG', 'G',\n'DMH', 'N',\n'DMI', 'X',\n'DMJ', 'X',\\
n'DML', 'X',\n'DMM', 'X',\n'DMN', 'X',\n'DMO', 'X'\
,\n'DMP', 'X',\n'DMQ', 'X',\n'DMR', 'X',\n'DMS', '\
X',\n'DMT', 'X',\n'DMV', 'X',\n'DMY', 'X',\n'DNC',\
 'X',\n'DND', 'X',\n'DNH', 'X',\n'DNJ', 'X',\n'DNN\
', 'X',\n'DNP', 'X',\n'DNQ', 'X',\n'DNR', 'X',\n'D\
O2', 'X',\n'DO3', 'X',\n'DOA', 'X',\n'DOB', 'X',\n\
'DOC', 'X',\n'DOH', 'D',\n'DOM', 'X',\n'DOS', 'X',\
\n'DOX', 'X',\n'DP5', 'X',\n'DP7', 'X',\n'DPA', 'X\
',\n'DPC', 'X',\n'DPD', 'X',\n'DPE', 'X',\n'DPG', \
'X',\n'DPH', 'F',\n'DPM', 'X',\n'DPN', 'F',\n'DPO'\
, 'X',\n'DPP', 'X',\n'DPR', 'P',\n'DPS', 'X',\n'DP\
T', 'X',\n'DPX', 'X',\n'DPY', 'X',\n'DPZ', 'X',\n'\
DQH', 'X',\n'DQN', 'X',\n'DR1', 'X',\n'DRB', 'X',\\
n'DRC', 'X',\n'DRI', 'X',\n'DRP', 'X',\n'DRT', 'X'\
,\n'DRU', 'X',\n'DSA', 'X',\n'DSB', 'X',\n'DSC', '\
X',\n'DSD', 'X',\n'DSE', 'S',\n'DSI', 'X',\n'DSN',\
 'S',\n'DSP', 'D',\n'DSR', 'X',\n'DSS', 'X',\n'DSX\
', 'X',\n'DSY', 'X',\n'DTB', 'X',\n'DTD', 'X',\n'D\
TH', 'T',\n'DTN', 'X',\n'DTO', 'X',\n'DTP', 'X',\n\
'DTQ', 'X',\n'DTR', 'W',\n'DTT', 'X',\n'DTY', 'Y',\
\n'DUD', 'X',\n'DUO', 'X',\n'DUR', 'X',\n'DUT', 'X\
',\n'DVA', 'V',\n'DVR', 'X',\n'DX9', 'X',\n'DXA', \
'X',\n'DXB', 'X',\n'DXC', 'X',\n'DXG', 'X',\n'DXX'\
, 'X',\n'DZF', 'X',\n'E09', 'X',\n'E20', 'X',\n'E2\
P', 'X',\n'E3G', 'X',\n'E4N', 'X',\n'E4P', 'X',\n'\
E64', 'X',\n'E6C', 'X',\n'E96', 'X',\n'E97', 'X',\\
n'EA2', 'X',\n'EAA', 'X',\n'EAP', 'X',\n'EBP', 'X'\
,\n'EBW', 'X',\n'ECO', 'X',\n'EDA', 'X',\n'EDC', '\
X',\n'EDE', 'X',\n'EDO', 'X',\n'EDR', 'X',\n'EEB',\
 'X',\n'EEE', 'X',\n'EFC', 'X',\n'EFZ', 'X',\n'EG1\
', 'X',\n'EG2', 'X',\n'EG3', 'X',\n'EGC', 'X',\n'E\
GL', 'X',\n'EHP', 'A',\n'EIC', 'X',\n'EJT', 'X',\n\
'ELA', 'X',\n'EMB', 'X',\n'EMC', 'X',\n'EMD', 'X',\
\n'EMM', 'X',\n'EMO', 'X',\n'EMP', 'X',\n'EMR', 'X\
',\n'ENA', 'X',\n'ENC', 'X',\n'ENH', 'X',\n'ENO', \
'X',\n'ENP', 'X',\n'EOA', 'X',\n'EOH', 'X',\n'EOT'\
, 'X',\n'EOX', 'X',\n'EPA', 'X',\n'EPE', 'X',\n'EP\
H', 'X',\n'EPI', 'X',\n'EPN', 'X',\n'EPO', 'X',\n'\
EPT', 'X',\n'EPU', 'X',\n'EPX', 'X',\n'EPY', 'X',\\
n'EQI', 'X',\n'EQP', 'X',\n'EQU', 'X',\n'ERG', 'X'\
,\n'ERI', 'X',\n'ERY', 'X',\n'ESC', 'X',\n'ESD', '\
X',\n'ESI', 'X',\n'ESO', 'X',\n'ESP', 'X',\n'EST',\
 'X',\n'ESX', 'X',\n'ETA', 'X',\n'ETC', 'X',\n'ETD\
', 'X',\n'ETF', 'X',\n'ETH', 'X',\n'ETI', 'X',\n'E\
TN', 'X',\n'ETO', 'X',\n'ETP', 'X',\n'ETR', 'X',\n\
'ETS', 'X',\n'ETY', 'X',\n'EU3', 'X',\n'EUG', 'X',\
\n'EYS', 'C',\n'F09', 'X',\n'F2B', 'X',\n'F3S', 'X\
',\n'F42', 'X',\n'F43', 'X',\n'F4S', 'X',\n'F6B', \
'X',\n'F6P', 'X',\n'F89', 'X',\n'FA1', 'X',\n'FA5'\
, 'F',\n'FAA', 'X',\n'FAB', 'X',\n'FAC', 'X',\n'FA\
D', 'X',\n'FAF', 'X',\n'FAG', 'X',\n'FAM', 'X',\n'\
FAR', 'X',\n'FAS', 'X',\n'FAT', 'X',\n'FBA', 'X',\\
n'FBE', 'X',\n'FBI', 'X',\n'FBP', 'X',\n'FBQ', 'X'\
,\n'FBS', 'X',\n'FBT', 'X',\n'FBU', 'X',\n'FCA', '\
X',\n'FCB', 'X',\n'FCI', 'X',\n'FCN', 'X',\n'FCO',\
 'X',\n'FCR', 'X',\n'FCT', 'X',\n'FCX', 'X',\n'FCY\
', 'C',\n'FD1', 'F',\n'FD2', 'F',\n'FD3', 'F',\n'F\
D4', 'F',\n'FDA', 'X',\n'FDC', 'X',\n'FDI', 'X',\n\
'FDP', 'X',\n'FDS', 'X',\n'FE2', 'X',\n'FEA', 'X',\
\n'FEL', 'X',\n'FEM', 'X',\n'FEN', 'X',\n'FEO', 'X\
',\n'FEP', 'X',\n'FER', 'X',\n'FES', 'X',\n'FFB', \
'X',\n'FFC', 'X',\n'FFF', 'X',\n'FFO', 'X',\n'FGL'\
, 'G',\n'FHB', 'X',\n'FHC', 'X',\n'FHP', 'X',\n'FH\
U', 'X',\n'FID', 'X',\n'FII', 'X',\n'FIP', 'X',\n'\
FK5', 'X',\n'FKA', 'X',\n'FKI', 'X',\n'FKP', 'X',\\
n'FL2', 'X',\n'FL9', 'X',\n'FLA', 'A',\n'FLC', 'X'\
,\n'FLD', 'X',\n'FLE', 'L',\n'FLF', 'X',\n'FLO', '\
X',\n'FLP', 'X',\n'FLT', 'Y',\n'FLU', 'X',\n'FLX',\
 'X',\n'FM1', 'X',\n'FM2', 'X',\n'FMA', 'X',\n'FMB\
', 'X',\n'FMC', 'X',\n'FME', 'M',\n'FMN', 'X',\n'F\
MP', 'X',\n'FMR', 'X',\n'FMS', 'X',\n'FMT', 'X',\n\
'FNE', 'X',\n'FNP', 'X',\n'FNS', 'X',\n'FOC', 'X',\
\n'FOE', 'X',\n'FOG', 'F',\n'FOH', 'X',\n'FOK', 'X\
',\n'FOL', 'X',\n'FON', 'X',\n'FOP', 'X',\n'FOR', \
'X',\n'FOS', 'X',\n'FPA', 'X',\n'FPC', 'X',\n'FPI'\
, 'X',\n'FPO', 'X',\n'FPP', 'X',\n'FPT', 'X',\n'FQ\
P', 'X',\n'FRA', 'X',\n'FRD', 'F',\n'FRU', 'X',\n'\
FS3', 'X',\n'FS4', 'X',\n'FSB', 'X',\n'FSO', 'X',\\
n'FSX', 'X',\n'FTC', 'X',\n'FTP', 'X',\n'FTR', 'W'\
,\n'FTT', 'X',\n'FTY', 'Y',\n'FUA', 'X',\n'FUC', '\
X',\n'FUM', 'X',\n'FUP', 'X',\n'FVF', 'X',\n'FXP',\
 'X',\n'FXV', 'X',\n'FYA', 'F',\n'G16', 'X',\n'G1P\
', 'X',\n'G20', 'X',\n'G21', 'X',\n'G23', 'X',\n'G\
26', 'X',\n'G28', 'X',\n'G2F', 'X',\n'G37', 'X',\n\
'G39', 'X',\n'G3H', 'X',\n'G3P', 'X',\n'G4D', 'X',\
\n'G6D', 'X',\n'G6P', 'X',\n'G6Q', 'X',\n'G7M', 'X\
',\n'GA2', 'X',\n'GAA', 'X',\n'GAB', 'X',\n'GAC', \
'X',\n'GAI', 'X',\n'GAL', 'X',\n'GAM', 'X',\n'GAN'\
, 'X',\n'GAO', 'X',\n'GAP', 'X',\n'GAR', 'G',\n'GA\
S', 'X',\n'GAT', 'X',\n'GBC', 'X',\n'GBI', 'X',\n'\
GBP', 'X',\n'GBS', 'X',\n'GBX', 'X',\n'GC4', 'X',\\
n'GCA', 'X',\n'GCD', 'X',\n'GCG', 'G',\n'GCH', 'G'\
,\n'GCK', 'X',\n'GCL', 'X',\n'GCM', 'X',\n'GCN', '\
X',\n'GCO', 'X',\n'GCP', 'X',\n'GCR', 'X',\n'GCS',\
 'X',\n'GCU', 'X',\n'GD3', 'X',\n'GDB', 'X',\n'GDM\
', 'X',\n'GDN', 'X',\n'GDP', 'X',\n'GDS', 'X',\n'G\
DU', 'X',\n'GE1', 'X',\n'GE2', 'X',\n'GE3', 'X',\n\
'GEA', 'X',\n'GEL', 'X',\n'GEM', 'X',\n'GEN', 'X',\
\n'GEP', 'X',\n'GER', 'X',\n'GFP', 'X',\n'GGB', 'X\
',\n'GGL', 'E',\n'GGP', 'X',\n'GHP', 'G',\n'GIP', \
'X',\n'GIS', 'X',\n'GKR', 'X',\n'GL2', 'X',\n'GL3'\
, 'G',\n'GL4', 'X',\n'GL5', 'X',\n'GL7', 'X',\n'GL\
9', 'X',\n'GLA', 'X',\n'GLB', 'X',\n'GLC', 'X',\n'\
GLD', 'X',\n'GLE', 'X',\n'GLF', 'X',\n'GLG', 'X',\\
n'GLH', 'Q',\n'GLI', 'X',\n'GLL', 'X',\n'GLM', 'G'\
,\n'GLN', 'Q',\n'GLO', 'X',\n'GLP', 'X',\n'GLR', '\
X',\n'GLS', 'X',\n'GLT', 'X',\n'GLU', 'E',\n'GLV',\
 'X',\n'GLW', 'X',\n'GLY', 'G',\n'GLZ', 'X',\n'GM1\
', 'X',\n'GMA', 'X',\n'GMC', 'X',\n'GMH', 'X',\n'G\
MP', 'X',\n'GMY', 'X',\n'GN7', 'X',\n'GNA', 'X',\n\
'GNB', 'X',\n'GNH', 'X',\n'GNP', 'X',\n'GNT', 'X',\
\n'GOA', 'X',\n'GOL', 'X',\n'GOX', 'X',\n'GP1', 'X\
',\n'GP3', 'X',\n'GP4', 'X',\n'GP6', 'X',\n'GP8', \
'X',\n'GPB', 'E',\n'GPC', 'X',\n'GPE', 'X',\n'GPG'\
, 'X',\n'GPI', 'X',\n'GPJ', 'X',\n'GPL', 'K',\n'GP\
M', 'X',\n'GPN', 'G',\n'GPP', 'X',\n'GPR', 'X',\n'\
GPS', 'X',\n'GPX', 'X',\n'GR1', 'X',\n'GR3', 'X',\\
n'GR4', 'X',\n'GSA', 'X',\n'GSB', 'X',\n'GSC', 'G'\
,\n'GSE', 'S',\n'GSH', 'X',\n'GSP', 'X',\n'GSR', '\
X',\n'GSS', 'X',\n'GT9', 'C',\n'GTA', 'X',\n'GTB',\
 'X',\n'GTD', 'X',\n'GTE', 'X',\n'GTH', 'T',\n'GTN\
', 'X',\n'GTO', 'X',\n'GTP', 'X',\n'GTR', 'X',\n'G\
TS', 'X',\n'GTT', 'X',\n'GTX', 'X',\n'GTZ', 'X',\n\
'GU7', 'X',\n'GUA', 'X',\n'GUD', 'X',\n'GUM', 'X',\
\n'GUN', 'X',\n'GUP', 'X',\n'GUR', 'X',\n'GW3', 'X\
',\n'GZZ', 'X',\n'H2B', 'X',\n'H2P', 'H',\n'H2S', \
'X',\n'H2U', 'X',\n'H4B', 'X',\n'H5M', 'P',\n'H5P'\
, 'X',\n'HAA', 'X',\n'HAB', 'X',\n'HAC', 'A',\n'HA\
D', 'X',\n'HAE', 'X',\n'HAG', 'X',\n'HAI', 'X',\n'\
HAM', 'X',\n'HAP', 'X',\n'HAQ', 'X',\n'HAR', 'R',\\
n'HAS', 'X',\n'HAV', 'V',\n'HAX', 'X',\n'HAZ', 'X'\
,\n'HBA', 'X',\n'HBC', 'X',\n'HBD', 'X',\n'HBI', '\
X',\n'HBO', 'X',\n'HBU', 'X',\n'HBY', 'X',\n'HC0',\
 'X',\n'HC1', 'X',\n'HC4', 'X',\n'HCA', 'X',\n'HCC\
', 'X',\n'HCI', 'X',\n'HCS', 'X',\n'HDA', 'X',\n'H\
DD', 'X',\n'HDF', 'X',\n'HDN', 'X',\n'HDS', 'X',\n\
'HDZ', 'X',\n'HE1', 'X',\n'HE6', 'X',\n'HEA', 'X',\
\n'HEB', 'X',\n'HEC', 'X',\n'HED', 'X',\n'HEE', 'X\
',\n'HEF', 'X',\n'HEG', 'X',\n'HEM', 'X',\n'HEN', \
'X',\n'HEO', 'X',\n'HEP', 'X',\n'HEU', 'X',\n'HEV'\
, 'X',\n'HEX', 'X',\n'HEZ', 'X',\n'HF1', 'X',\n'HF\
A', 'X',\n'HFP', 'X',\n'HGA', 'Q',\n'HGB', 'X',\n'\
HGC', 'X',\n'HGI', 'X',\n'HGU', 'X',\n'HHO', 'X',\\
n'HHP', 'X',\n'HIB', 'X',\n'HIC', 'H',\n'HII', 'X'\
,\n'HIN', 'X',\n'HIO', 'X',\n'HIP', 'H',\n'HIS', '\
H',\n'HLE', 'X',\n'HLT', 'X',\n'HMA', 'A',\n'HMB',\
 'X',\n'HMC', 'X',\n'HMD', 'X',\n'HMF', 'A',\n'HMG\
', 'X',\n'HMH', 'X',\n'HMI', 'L',\n'HMM', 'X',\n'H\
MN', 'X',\n'HMO', 'X',\n'HMP', 'X',\n'HMR', 'R',\n\
'HNI', 'X',\n'HNP', 'X',\n'HOA', 'X',\n'HOE', 'X',\
\n'HOH', 'X',\n'HOM', 'X',\n'HOP', 'X',\n'HOQ', 'X\
',\n'HP1', 'A',\n'HP2', 'A',\n'HP3', 'X',\n'HPA', \
'X',\n'HPB', 'X',\n'HPC', 'X',\n'HPD', 'X',\n'HPE'\
, 'A',\n'HPG', 'X',\n'HPH', 'F',\n'HPP', 'X',\n'HP\
Q', 'F',\n'HPR', 'X',\n'HPT', 'X',\n'HPY', 'X',\n'\
HQO', 'X',\n'HQQ', 'X',\n'HQU', 'X',\n'HRG', 'R',\\
n'HRI', 'X',\n'HSA', 'X',\n'HSE', 'S',\n'HSF', 'X'\
,\n'HSM', 'X',\n'HSO', 'H',\n'HSP', 'X',\n'HT1', '\
X',\n'HT2', 'X',\n'HTA', 'X',\n'HTL', 'X',\n'HTO',\
 'X',\n'HTP', 'X',\n'HTR', 'W',\n'HUP', 'X',\n'HUX\
', 'X',\n'HV5', 'A',\n'HV7', 'X',\n'HV8', 'X',\n'H\
XA', 'X',\n'HXC', 'X',\n'HXP', 'X',\n'HY1', 'X',\n\
'HYA', 'X',\n'HYB', 'X',\n'HYD', 'X',\n'HYG', 'X',\
\n'HYP', 'P',\n'I06', 'X',\n'I10', 'X',\n'I11', 'X\
',\n'I17', 'X',\n'I2P', 'X',\n'I3N', 'X',\n'I3P', \
'X',\n'I40', 'X',\n'I48', 'X',\n'I4B', 'X',\n'I52'\
, 'X',\n'I5P', 'X',\n'I84', 'G',\n'IAG', 'G',\n'IA\
S', 'X',\n'IB2', 'X',\n'IBB', 'X',\n'IBP', 'X',\n'\
IBR', 'X',\n'IBS', 'X',\n'IBZ', 'X',\n'IC1', 'X',\\
n'ICA', 'X',\n'ICI', 'X',\n'ICL', 'X',\n'ICP', 'X'\
,\n'ICT', 'X',\n'ICU', 'X',\n'ID2', 'X',\n'IDC', '\
X',\n'IDG', 'X',\n'IDH', 'X',\n'IDM', 'X',\n'IDO',\
 'X',\n'IDP', 'X',\n'IDR', 'X',\n'IDS', 'X',\n'IDT\
', 'X',\n'IDU', 'X',\n'IFG', 'X',\n'IFP', 'X',\n'I\
GL', 'X',\n'IGN', 'X',\n'IGP', 'X',\n'IGU', 'X',\n\
'IH1', 'X',\n'IH2', 'X',\n'IH3', 'X',\n'IHB', 'X',\
\n'IHN', 'X',\n'IHP', 'X',\n'IIC', 'X',\n'IIL', 'I\
',\n'IIP', 'X',\n'IK2', 'X',\n'IKT', 'X',\n'ILA', \
'I',\n'ILE', 'I',\n'ILG', 'X',\n'ILO', 'X',\n'ILX'\
, 'I',\n'IM1', 'X',\n'IM2', 'X',\n'IMC', 'X',\n'IM\
D', 'X',\n'IME', 'X',\n'IMF', 'X',\n'IMG', 'X',\n'\
IMH', 'X',\n'IMI', 'X',\n'IML', 'I',\n'IMM', 'X',\\
n'IMN', 'X',\n'IMO', 'X',\n'IMP', 'X',\n'IMR', 'X'\
,\n'IMU', 'X',\n'IN0', 'D',\n'IN1', 'R',\n'IN2', '\
K',\n'IN3', 'L',\n'IN4', 'X',\n'IN5', 'A',\n'IN6',\
 'L',\n'IN7', 'X',\n'IN8', 'X',\n'IN9', 'X',\n'INA\
', 'L',\n'INB', 'X',\n'INC', 'X',\n'IND', 'X',\n'I\
NE', 'X',\n'INF', 'F',\n'ING', 'F',\n'INH', 'R',\n\
'INI', 'X',\n'INJ', 'X',\n'INK', 'X',\n'INL', 'X',\
\n'INM', 'X',\n'INN', 'A',\n'INO', 'X',\n'INP', 'X\
',\n'INQ', 'X',\n'INR', 'X',\n'INS', 'X',\n'INT', \
'V',\n'INU', 'X',\n'INV', 'X',\n'INW', 'X',\n'INX'\
, 'X',\n'INY', 'X',\n'INZ', 'X',\n'IOA', 'X',\n'IO\
B', 'X',\n'IOC', 'X',\n'IOD', 'X',\n'IOE', 'X',\n'\
IOF', 'X',\n'IOH', 'X',\n'IOL', 'X',\n'IOP', 'X',\\
n'IP1', 'X',\n'IP2', 'X',\n'IP3', 'X',\n'IP4', 'X'\
,\n'IPA', 'X',\n'IPB', 'X',\n'IPD', 'X',\n'IPG', '\
G',\n'IPH', 'X',\n'IPL', 'X',\n'IPM', 'X',\n'IPN',\
 'X',\n'IPO', 'F',\n'IPP', 'X',\n'IPS', 'X',\n'IPT\
', 'X',\n'IPU', 'X',\n'IPY', 'A',\n'IQB', 'X',\n'I\
QP', 'X',\n'IQS', 'X',\n'IR3', 'X',\n'IRI', 'X',\n\
'IRP', 'X',\n'ISA', 'X',\n'ISF', 'X',\n'ISO', 'X',\
\n'ISP', 'X',\n'ISQ', 'X',\n'ISU', 'X',\n'ITM', 'X\
',\n'ITP', 'X',\n'ITR', 'W',\n'ITS', 'X',\n'ITU', \
'X',\n'IU5', 'X',\n'IUM', 'X',\n'IUR', 'X',\n'IVA'\
, 'X',\n'IYG', 'G',\n'IYR', 'Y',\n'J77', 'X',\n'J7\
8', 'X',\n'J80', 'X',\n'JE2', 'X',\n'JEN', 'X',\n'\
JST', 'X',\n'K21', 'X',\n'KAH', 'X',\n'KAI', 'X',\\
n'KAM', 'X',\n'KAN', 'X',\n'KAP', 'X',\n'KCP', 'X'\
,\n'KCX', 'K',\n'KDO', 'X',\n'KEF', 'X',\n'KET', '\
X',\n'KGR', 'X',\n'KH1', 'X',\n'KIF', 'X',\n'KIV',\
 'V',\n'KNI', 'X',\n'KPH', 'K',\n'KTH', 'X',\n'KTN\
', 'X',\n'KTP', 'X',\n'KWT', 'X',\n'L04', 'X',\n'L\
1P', 'X',\n'L24', 'E',\n'L2P', 'X',\n'L34', 'E',\n\
'L37', 'E',\n'L3P', 'X',\n'L4P', 'X',\n'L75', 'X',\
\n'LAC', 'X',\n'LAD', 'X',\n'LAK', 'X',\n'LAM', 'X\
',\n'LAR', 'X',\n'LAT', 'X',\n'LAX', 'X',\n'LCO', \
'X',\n'LCP', 'X',\n'LCS', 'X',\n'LDA', 'X',\n'LDO'\
, 'L',\n'LDP', 'X',\n'LEA', 'X',\n'LEO', 'X',\n'LE\
U', 'L',\n'LG2', 'X',\n'LG6', 'X',\n'LGC', 'X',\n'\
LGP', 'X',\n'LHG', 'X',\n'LHY', 'F',\n'LI1', 'X',\\
n'LIG', 'X',\n'LIL', 'X',\n'LIM', 'X',\n'LIN', 'X'\
,\n'LIO', 'X',\n'LIP', 'X',\n'LLA', 'X',\n'LLP', '\
K',\n'LLY', 'K',\n'LMG', 'X',\n'LML', 'X',\n'LMT',\
 'X',\n'LMU', 'X',\n'LMZ', 'X',\n'LNK', 'X',\n'LNL\
', 'X',\n'LNO', 'X',\n'LOF', 'X',\n'LOL', 'L',\n'L\
OM', 'X',\n'LOR', 'X',\n'LOS', 'X',\n'LOV', 'L',\n\
'LOX', 'X',\n'LP1', 'X',\n'LP2', 'R',\n'LPA', 'X',\
\n'LPC', 'X',\n'LPF', 'X',\n'LPL', 'X',\n'LPM', 'X\
',\n'LPP', 'X',\n'LRB', 'X',\n'LRU', 'X',\n'LS1', \
'X',\n'LS2', 'X',\n'LS3', 'X',\n'LS4', 'X',\n'LS5'\
, 'X',\n'LTA', 'X',\n'LTL', 'X',\n'LTR', 'W',\n'LU\
M', 'X',\n'LVS', 'L',\n'LXC', 'X',\n'LY2', 'X',\n'\
LY3', 'X',\n'LYA', 'X',\n'LYB', 'X',\n'LYC', 'X',\\
n'LYD', 'X',\n'LYM', 'K',\n'LYN', 'X',\n'LYS', 'K'\
,\n'LYT', 'X',\n'LYW', 'X',\n'LYZ', 'K',\n'M1A', '\
X',\n'M1G', 'X',\n'M2G', 'X',\n'M3L', 'K',\n'M6P',\
 'X',\n'M6T', 'X',\n'M7G', 'X',\n'MA1', 'X',\n'MA2\
', 'X',\n'MA3', 'X',\n'MA4', 'X',\n'MA6', 'X',\n'M\
AA', 'A',\n'MAB', 'X',\n'MAC', 'X',\n'MAE', 'X',\n\
'MAG', 'X',\n'MAH', 'X',\n'MAI', 'R',\n'MAK', 'X',\
\n'MAL', 'X',\n'MAM', 'X',\n'MAN', 'X',\n'MAO', 'X\
',\n'MAP', 'X',\n'MAR', 'X',\n'MAS', 'X',\n'MAT', \
'X',\n'MAU', 'X',\n'MAZ', 'X',\n'MBA', 'X',\n'MBD'\
, 'X',\n'MBG', 'X',\n'MBH', 'X',\n'MBN', 'X',\n'MB\
O', 'X',\n'MBR', 'X',\n'MBS', 'X',\n'MBV', 'X',\n'\
MBZ', 'X',\n'MCA', 'X',\n'MCD', 'X',\n'MCE', 'X',\\
n'MCG', 'G',\n'MCI', 'X',\n'MCN', 'X',\n'MCP', 'X'\
,\n'MCT', 'X',\n'MCY', 'X',\n'MD2', 'X',\n'MDA', '\
X',\n'MDC', 'X',\n'MDG', 'X',\n'MDH', 'X',\n'MDL',\
 'X',\n'MDM', 'X',\n'MDN', 'X',\n'MDP', 'X',\n'ME6\
', 'X',\n'MEB', 'X',\n'MEC', 'X',\n'MEL', 'X',\n'M\
EN', 'N',\n'MEP', 'X',\n'MER', 'X',\n'MES', 'X',\n\
'MET', 'M',\n'MEV', 'X',\n'MF2', 'X',\n'MF3', 'M',\
\n'MFB', 'X',\n'MFD', 'X',\n'MFU', 'X',\n'MG7', 'X\
',\n'MGA', 'X',\n'MGB', 'X',\n'MGD', 'X',\n'MGG', \
'R',\n'MGL', 'X',\n'MGN', 'Q',\n'MGO', 'X',\n'MGP'\
, 'X',\n'MGR', 'X',\n'MGS', 'X',\n'MGT', 'X',\n'MG\
U', 'X',\n'MGY', 'G',\n'MHB', 'X',\n'MHF', 'X',\n'\
MHL', 'L',\n'MHM', 'X',\n'MHO', 'M',\n'MHS', 'H',\\
n'MHZ', 'X',\n'MIA', 'X',\n'MIC', 'X',\n'MID', 'X'\
,\n'MIL', 'X',\n'MIM', 'X',\n'MIN', 'G',\n'MIP', '\
X',\n'MIS', 'S',\n'MIT', 'X',\n'MJI', 'X',\n'MK1',\
 'X',\n'MKC', 'X',\n'MLA', 'X',\n'MLC', 'X',\n'MLE\
', 'L',\n'MLN', 'X',\n'MLT', 'X',\n'MLY', 'K',\n'M\
LZ', 'K',\n'MM3', 'X',\n'MM4', 'X',\n'MMA', 'X',\n\
'MMC', 'X',\n'MME', 'M',\n'MMO', 'R',\n'MMP', 'X',\
\n'MMQ', 'X',\n'MMT', 'X',\n'MN1', 'X',\n'MN2', 'X\
',\n'MN3', 'X',\n'MN5', 'X',\n'MN7', 'X',\n'MN8', \
'X',\n'MNA', 'X',\n'MNB', 'X',\n'MNC', 'X',\n'MNG'\
, 'X',\n'MNL', 'L',\n'MNO', 'X',\n'MNP', 'X',\n'MN\
Q', 'X',\n'MNS', 'X',\n'MNT', 'X',\n'MNV', 'V',\n'\
MO1', 'X',\n'MO2', 'X',\n'MO3', 'X',\n'MO4', 'X',\\
n'MO5', 'X',\n'MO6', 'X',\n'MOA', 'X',\n'MOB', 'X'\
,\n'MOC', 'X',\n'MOE', 'X',\n'MOG', 'X',\n'MOH', '\
X',\n'MOL', 'X',\n'MOO', 'X',\n'MOP', 'X',\n'MOR',\
 'X',\n'MOS', 'X',\n'MOT', 'X',\n'MOX', 'X',\n'MP1\
', 'X',\n'MP3', 'X',\n'MPA', 'X',\n'MPB', 'X',\n'M\
PC', 'X',\n'MPD', 'X',\n'MPG', 'X',\n'MPH', 'M',\n\
'MPI', 'X',\n'MPJ', 'M',\n'MPL', 'X',\n'MPN', 'X',\
\n'MPO', 'X',\n'MPP', 'X',\n'MPQ', 'G',\n'MPR', 'X\
',\n'MPS', 'X',\n'MQ0', 'X',\n'MQ7', 'X',\n'MQ8', \
'X',\n'MQ9', 'X',\n'MQI', 'X',\n'MR2', 'X',\n'MRC'\
, 'X',\n'MRM', 'X',\n'MRP', 'X',\n'MS2', 'X',\n'MS\
A', 'X',\n'MSB', 'X',\n'MSD', 'X',\n'MSE', 'M',\n'\
MSF', 'X',\n'MSI', 'X',\n'MSO', 'M',\n'MSQ', 'X',\\
n'MST', 'X',\n'MSU', 'X',\n'MTA', 'X',\n'MTB', 'X'\
,\n'MTC', 'X',\n'MTD', 'X',\n'MTE', 'X',\n'MTF', '\
X',\n'MTG', 'X',\n'MTO', 'X',\n'MTS', 'X',\n'MTT',\
 'X',\n'MTX', 'X',\n'MTY', 'Y',\n'MUG', 'X',\n'MUP\
', 'X',\n'MUR', 'X',\n'MVA', 'V',\n'MW1', 'X',\n'M\
W2', 'X',\n'MXA', 'X',\n'MXY', 'X',\n'MYA', 'X',\n\
'MYC', 'X',\n'MYG', 'X',\n'MYR', 'X',\n'MYS', 'X',\
\n'MYT', 'X',\n'MZM', 'X',\n'N1T', 'X',\n'N25', 'X\
',\n'N2B', 'X',\n'N3T', 'X',\n'N4B', 'X',\n'NA2', \
'X',\n'NA5', 'X',\n'NA6', 'X',\n'NAA', 'X',\n'NAB'\
, 'X',\n'NAC', 'X',\n'NAD', 'X',\n'NAE', 'X',\n'NA\
F', 'X',\n'NAG', 'X',\n'NAH', 'X',\n'NAI', 'X',\n'\
NAL', 'A',\n'NAM', 'A',\n'NAN', 'X',\n'NAO', 'X',\\
n'NAP', 'X',\n'NAQ', 'X',\n'NAR', 'X',\n'NAS', 'X'\
,\n'NAU', 'X',\n'NAV', 'X',\n'NAW', 'X',\n'NAX', '\
X',\n'NAY', 'X',\n'NBA', 'X',\n'NBD', 'X',\n'NBE',\
 'X',\n'NBG', 'X',\n'NBN', 'X',\n'NBP', 'X',\n'NBS\
', 'X',\n'NBU', 'X',\n'NCA', 'X',\n'NCB', 'A',\n'N\
CD', 'X',\n'NCH', 'X',\n'NCM', 'X',\n'NCN', 'X',\n\
'NCO', 'X',\n'NCR', 'X',\n'NCS', 'X',\n'ND4', 'X',\
\n'NDA', 'X',\n'NDC', 'X',\n'NDD', 'X',\n'NDO', 'X\
',\n'NDP', 'X',\n'NDT', 'X',\n'NEA', 'X',\n'NEB', \
'X',\n'NED', 'X',\n'NEM', 'H',\n'NEN', 'X',\n'NEO'\
, 'X',\n'NEP', 'H',\n'NEQ', 'X',\n'NES', 'X',\n'NE\
T', 'X',\n'NEV', 'X',\n'NFA', 'F',\n'NFE', 'X',\n'\
NFG', 'X',\n'NFP', 'X',\n'NFS', 'X',\n'NG6', 'X',\\
n'NGA', 'X',\n'NGL', 'X',\n'NGM', 'X',\n'NGO', 'X'\
,\n'NGP', 'X',\n'NGT', 'X',\n'NGU', 'X',\n'NH2', '\
X',\n'NH3', 'X',\n'NH4', 'X',\n'NHD', 'X',\n'NHE',\
 'X',\n'NHM', 'X',\n'NHP', 'X',\n'NHR', 'X',\n'NHS\
', 'X',\n'NI1', 'X',\n'NI2', 'X',\n'NIC', 'X',\n'N\
ID', 'X',\n'NIK', 'X',\n'NIO', 'X',\n'NIP', 'X',\n\
'NIT', 'X',\n'NIU', 'X',\n'NIY', 'Y',\n'NLA', 'X',\
\n'NLE', 'L',\n'NLG', 'X',\n'NLN', 'L',\n'NLP', 'L\
',\n'NM1', 'X',\n'NMA', 'A',\n'NMB', 'X',\n'NMC', \
'G',\n'NMD', 'X',\n'NME', 'X',\n'NMN', 'X',\n'NMO'\
, 'X',\n'NMQ', 'X',\n'NMX', 'X',\n'NMY', 'X',\n'NN\
H', 'R',\n'NNO', 'X',\n'NO2', 'X',\n'NO3', 'X',\n'\
NOA', 'X',\n'NOD', 'X',\n'NOJ', 'X',\n'NON', 'X',\\
n'NOP', 'X',\n'NOR', 'X',\n'NOS', 'X',\n'NOV', 'X'\
,\n'NOX', 'X',\n'NP3', 'X',\n'NPA', 'X',\n'NPC', '\
X',\n'NPD', 'X',\n'NPE', 'X',\n'NPF', 'X',\n'NPH',\
 'C',\n'NPI', 'X',\n'NPL', 'X',\n'NPN', 'X',\n'NPO\
', 'X',\n'NPP', 'X',\n'NPT', 'X',\n'NPY', 'X',\n'N\
RG', 'R',\n'NRI', 'X',\n'NS1', 'X',\n'NS5', 'X',\n\
'NSP', 'X',\n'NTA', 'X',\n'NTB', 'X',\n'NTC', 'X',\
\n'NTH', 'X',\n'NTM', 'X',\n'NTP', 'X',\n'NTS', 'X\
',\n'NTU', 'X',\n'NTZ', 'X',\n'NU1', 'X',\n'NVA', \
'V',\n'NVI', 'X',\n'NVP', 'X',\n'NW1', 'X',\n'NYP'\
, 'X',\n'O4M', 'X',\n'OAA', 'X',\n'OAI', 'X',\n'OA\
P', 'X',\n'OAR', 'X',\n'OAS', 'S',\n'OBA', 'X',\n'\
OBN', 'X',\n'OC1', 'X',\n'OC2', 'X',\n'OC3', 'X',\\
n'OC4', 'X',\n'OC5', 'X',\n'OC6', 'X',\n'OC7', 'X'\
,\n'OCL', 'X',\n'OCM', 'X',\n'OCN', 'X',\n'OCO', '\
X',\n'OCP', 'X',\n'OCS', 'C',\n'OCT', 'X',\n'OCV',\
 'K',\n'OCY', 'C',\n'ODA', 'X',\n'ODS', 'X',\n'OES\
', 'X',\n'OET', 'X',\n'OF1', 'X',\n'OF2', 'X',\n'O\
F3', 'X',\n'OFL', 'X',\n'OFO', 'X',\n'OHE', 'X',\n\
'OHO', 'X',\n'OHT', 'X',\n'OIC', 'X',\n'OIP', 'X',\
\n'OKA', 'X',\n'OLA', 'X',\n'OLE', 'X',\n'OLI', 'X\
',\n'OLO', 'X',\n'OMB', 'X',\n'OMC', 'X',\n'OMD', \
'X',\n'OME', 'X',\n'OMG', 'X',\n'OMP', 'X',\n'OMT'\
, 'M',\n'OMU', 'X',\n'ONE', 'X',\n'ONL', 'L',\n'ON\
P', 'X',\n'OPA', 'X',\n'OPD', 'X',\n'OPE', 'X',\n'\
OPG', 'X',\n'OPH', 'X',\n'OPN', 'X',\n'OPP', 'X',\\
n'OPR', 'R',\n'ORN', 'X',\n'ORO', 'X',\n'ORP', 'X'\
,\n'OSB', 'X',\n'OSS', 'X',\n'OTA', 'X',\n'OTB', '\
X',\n'OTE', 'X',\n'OTG', 'X',\n'OUT', 'X',\n'OVA',\
 'X',\n'OWQ', 'X',\n'OXA', 'X',\n'OXE', 'X',\n'OXI\
', 'X',\n'OXL', 'X',\n'OXM', 'X',\n'OXN', 'X',\n'O\
XO', 'X',\n'OXP', 'X',\n'OXS', 'X',\n'OXY', 'X',\n\
'P11', 'A',\n'P24', 'X',\n'P28', 'X',\n'P2P', 'X',\
\n'P2U', 'X',\n'P3M', 'X',\n'P4C', 'X',\n'P4P', 'X\
',\n'P5P', 'X',\n'P6G', 'X',\n'PA1', 'X',\n'PA2', \
'X',\n'PA3', 'X',\n'PA4', 'X',\n'PA5', 'X',\n'PAA'\
, 'X',\n'PAB', 'X',\n'PAC', 'X',\n'PAD', 'X',\n'PA\
E', 'X',\n'PAG', 'X',\n'PAH', 'X',\n'PAI', 'X',\n'\
PAL', 'D',\n'PAM', 'X',\n'PAN', 'X',\n'PAO', 'X',\\
n'PAP', 'A',\n'PAQ', 'F',\n'PAR', 'X',\n'PAS', 'X'\
,\n'PAT', 'W',\n'PBA', 'X',\n'PBB', 'X',\n'PBC', '\
X',\n'PBF', 'F',\n'PBG', 'X',\n'PBI', 'X',\n'PBM',\
 'X',\n'PBN', 'X',\n'PBP', 'X',\n'PBR', 'X',\n'PBZ\
', 'X',\n'PC2', 'X',\n'PCA', 'E',\n'PCB', 'X',\n'P\
CD', 'X',\n'PCE', 'X',\n'PCG', 'X',\n'PCH', 'X',\n\
'PCL', 'X',\n'PCM', 'X',\n'PCP', 'X',\n'PCR', 'X',\
\n'PCS', 'X',\n'PCU', 'X',\n'PCV', 'X',\n'PCY', 'X\
',\n'PD1', 'X',\n'PDA', 'X',\n'PDC', 'X',\n'PDD', \
'A',\n'PDE', 'A',\n'PDI', 'X',\n'PDL', 'A',\n'PDN'\
, 'X',\n'PDO', 'X',\n'PDP', 'X',\n'PDT', 'X',\n'PD\
U', 'X',\n'PE2', 'X',\n'PE6', 'X',\n'PEA', 'X',\n'\
PEB', 'X',\n'PEC', 'X',\n'PED', 'X',\n'PEE', 'X',\\
n'PEF', 'X',\n'PEG', 'X',\n'PEL', 'X',\n'PEO', 'X'\
,\n'PEP', 'X',\n'PEQ', 'X',\n'PER', 'X',\n'PET', '\
X',\n'PFB', 'X',\n'PFC', 'X',\n'PFG', 'X',\n'PFL',\
 'X',\n'PFM', 'X',\n'PFZ', 'X',\n'PG4', 'X',\n'PG5\
', 'X',\n'PG6', 'X',\n'PGA', 'X',\n'PGC', 'X',\n'P\
GD', 'X',\n'PGE', 'X',\n'PGG', 'G',\n'PGH', 'X',\n\
'PGL', 'X',\n'PGO', 'X',\n'PGP', 'X',\n'PGQ', 'X',\
\n'PGR', 'X',\n'PGS', 'X',\n'PGU', 'X',\n'PGX', 'X\
',\n'PGY', 'G',\n'PH1', 'X',\n'PH2', 'X',\n'PH3', \
'X',\n'PHA', 'F',\n'PHB', 'X',\n'PHC', 'X',\n'PHD'\
, 'X',\n'PHE', 'F',\n'PHG', 'X',\n'PHH', 'X',\n'PH\
I', 'F',\n'PHL', 'F',\n'PHM', 'X',\n'PHN', 'X',\n'\
PHO', 'X',\n'PHP', 'X',\n'PHQ', 'X',\n'PHS', 'H',\\
n'PHT', 'X',\n'PHW', 'P',\n'PHY', 'X',\n'PI1', 'X'\
,\n'PI2', 'X',\n'PI3', 'X',\n'PI4', 'X',\n'PI5', '\
X',\n'PI6', 'X',\n'PI7', 'X',\n'PI8', 'X',\n'PI9',\
 'X',\n'PIA', 'X',\n'PIB', 'X',\n'PIC', 'X',\n'PID\
', 'X',\n'PIG', 'X',\n'PIH', 'X',\n'PIM', 'X',\n'P\
IN', 'X',\n'PIO', 'X',\n'PIP', 'X',\n'PIQ', 'X',\n\
'PIR', 'X',\n'PIV', 'X',\n'PKF', 'X',\n'PL1', 'X',\
\n'PL9', 'X',\n'PLA', 'D',\n'PLC', 'X',\n'PLE', 'L\
',\n'PLG', 'G',\n'PLH', 'X',\n'PLM', 'X',\n'PLP', \
'X',\n'PLS', 'S',\n'PLT', 'W',\n'PLU', 'L',\n'PLY'\
, 'X',\n'PMA', 'X',\n'PMB', 'X',\n'PMC', 'X',\n'PM\
E', 'F',\n'PML', 'X',\n'PMM', 'X',\n'PMO', 'X',\n'\
PMP', 'X',\n'PMS', 'X',\n'PMY', 'X',\n'PN2', 'X',\\
n'PNA', 'X',\n'PNB', 'X',\n'PNC', 'G',\n'PND', 'X'\
,\n'PNE', 'A',\n'PNF', 'X',\n'PNG', 'X',\n'PNI', '\
X',\n'PNL', 'X',\n'PNM', 'X',\n'PNN', 'X',\n'PNO',\
 'X',\n'PNP', 'X',\n'PNQ', 'X',\n'PNS', 'X',\n'PNT\
', 'X',\n'PNU', 'X',\n'PO2', 'X',\n'PO4', 'X',\n'P\
OB', 'X',\n'POC', 'X',\n'POL', 'X',\n'POM', 'P',\n\
'PON', 'X',\n'POP', 'X',\n'POR', 'X',\n'POS', 'X',\
\n'PP1', 'X',\n'PP2', 'X',\n'PP3', 'A',\n'PP4', 'X\
',\n'PP5', 'X',\n'PP6', 'X',\n'PP7', 'X',\n'PP8', \
'N',\n'PP9', 'X',\n'PPB', 'X',\n'PPC', 'X',\n'PPD'\
, 'X',\n'PPE', 'E',\n'PPG', 'X',\n'PPH', 'F',\n'PP\
I', 'X',\n'PPJ', 'V',\n'PPL', 'X',\n'PPM', 'X',\n'\
PPN', 'A',\n'PPO', 'X',\n'PPP', 'X',\n'PPQ', 'X',\\
n'PPR', 'X',\n'PPS', 'X',\n'PPT', 'X',\n'PPU', 'X'\
,\n'PPX', 'F',\n'PPY', 'X',\n'PPZ', 'X',\n'PQ0', '\
X',\n'PQN', 'X',\n'PQQ', 'X',\n'PR1', 'X',\n'PR2',\
 'X',\n'PR3', 'X',\n'PRA', 'X',\n'PRB', 'X',\n'PRC\
', 'X',\n'PRD', 'X',\n'PRE', 'X',\n'PRF', 'X',\n'P\
RH', 'X',\n'PRI', 'P',\n'PRL', 'X',\n'PRN', 'X',\n\
'PRO', 'P',\n'PRP', 'X',\n'PRR', 'A',\n'PRS', 'P',\
\n'PRZ', 'X',\n'PS0', 'X',\n'PSA', 'X',\n'PSD', 'X\
',\n'PSE', 'X',\n'PSF', 'S',\n'PSG', 'X',\n'PSI', \
'X',\n'PSO', 'X',\n'PSQ', 'X',\n'PSS', 'X',\n'PST'\
, 'X',\n'PSU', 'X',\n'PT1', 'X',\n'PT3', 'X',\n'PT\
A', 'X',\n'PTC', 'X',\n'PTD', 'X',\n'PTE', 'X',\n'\
PTH', 'Y',\n'PTL', 'X',\n'PTM', 'Y',\n'PTN', 'X',\\
n'PTO', 'X',\n'PTP', 'X',\n'PTR', 'Y',\n'PTS', 'X'\
,\n'PTT', 'X',\n'PTU', 'X',\n'PTY', 'X',\n'PUA', '\
X',\n'PUB', 'X',\n'PUR', 'X',\n'PUT', 'X',\n'PVA',\
 'X',\n'PVB', 'X',\n'PVH', 'H',\n'PVL', 'X',\n'PXA\
', 'X',\n'PXF', 'X',\n'PXG', 'X',\n'PXP', 'X',\n'P\
XY', 'X',\n'PXZ', 'X',\n'PY2', 'X',\n'PY4', 'X',\n\
'PY5', 'X',\n'PY6', 'X',\n'PYA', 'A',\n'PYC', 'X',\
\n'PYD', 'X',\n'PYE', 'X',\n'PYL', 'X',\n'PYM', 'X\
',\n'PYO', 'X',\n'PYP', 'X',\n'PYQ', 'X',\n'PYR', \
'X',\n'PYS', 'X',\n'PYT', 'X',\n'PYX', 'X',\n'PYY'\
, 'X',\n'PYZ', 'X',\n'PZQ', 'X',\n'Q82', 'X',\n'QN\
C', 'X',\n'QND', 'X',\n'QSI', 'Q',\n'QTR', 'X',\n'\
QUA', 'X',\n'QUE', 'X',\n'QUI', 'X',\n'QUO', 'X',\\
n'R11', 'X',\n'R12', 'X',\n'R13', 'X',\n'R18', 'X'\
,\n'R1P', 'X',\n'R56', 'X',\n'R5P', 'X',\n'RA2', '\
X',\n'RAD', 'X',\n'RAI', 'X',\n'RAL', 'X',\n'RAM',\
 'X',\n'RAN', 'X',\n'RAP', 'X',\n'RBF', 'X',\n'RBU\
', 'X',\n'RCA', 'X',\n'RCL', 'X',\n'RCO', 'X',\n'R\
DC', 'X',\n'RDF', 'W',\n'RE9', 'X',\n'REA', 'X',\n\
'RED', 'K',\n'REO', 'X',\n'REP', 'X',\n'RET', 'X',\
\n'RFA', 'X',\n'RFB', 'X',\n'RFL', 'X',\n'RFP', 'X\
',\n'RG1', 'X',\n'RGS', 'X',\n'RH1', 'X',\n'RHA', \
'X',\n'RHC', 'X',\n'RHD', 'X',\n'RHM', 'X',\n'RHO'\
, 'X',\n'RHQ', 'X',\n'RHS', 'X',\n'RIA', 'X',\n'RI\
B', 'X',\n'RIC', 'X',\n'RIF', 'X',\n'RIN', 'X',\n'\
RIP', 'X',\n'RIT', 'X',\n'RMB', 'X',\n'RMN', 'X',\\
n'RMP', 'X',\n'RNG', 'X',\n'RNS', 'X',\n'RNT', 'X'\
,\n'RO2', 'X',\n'RO4', 'X',\n'ROC', 'N',\n'ROI', '\
X',\n'ROM', 'X',\n'RON', 'V',\n'ROP', 'X',\n'ROS',\
 'X',\n'ROX', 'X',\n'RPA', 'X',\n'RPD', 'X',\n'RPH\
', 'X',\n'RPL', 'X',\n'RPP', 'X',\n'RPR', 'X',\n'R\
PX', 'X',\n'RQ3', 'X',\n'RR1', 'X',\n'RR6', 'X',\n\
'RRS', 'X',\n'RS1', 'X',\n'RS2', 'X',\n'RS7', 'X',\
\n'RSS', 'X',\n'RTA', 'X',\n'RTB', 'X',\n'RTC', 'X\
',\n'RTL', 'X',\n'RUB', 'X',\n'RUN', 'X',\n'RWJ', \
'X',\n'RXP', 'X',\n'S02', 'X',\n'S11', 'X',\n'S1H'\
, 'S',\n'S27', 'X',\n'S2C', 'C',\n'S3P', 'X',\n'S4\
U', 'X',\n'S57', 'X',\n'S58', 'X',\n'S5H', 'X',\n'\
S6G', 'X',\n'S80', 'X',\n'SAA', 'X',\n'SAB', 'X',\\
n'SAC', 'S',\n'SAD', 'X',\n'SAE', 'X',\n'SAF', 'X'\
,\n'SAH', 'C',\n'SAI', 'C',\n'SAL', 'X',\n'SAM', '\
M',\n'SAN', 'X',\n'SAP', 'X',\n'SAR', 'X',\n'SAS',\
 'X',\n'SB1', 'X',\n'SB2', 'X',\n'SB3', 'X',\n'SB4\
', 'X',\n'SB5', 'X',\n'SB6', 'X',\n'SBA', 'L',\n'S\
BB', 'X',\n'SBD', 'A',\n'SBI', 'X',\n'SBL', 'A',\n\
'SBN', 'X',\n'SBO', 'X',\n'SBR', 'X',\n'SBS', 'X',\
\n'SBT', 'X',\n'SBU', 'X',\n'SBX', 'X',\n'SC4', 'X\
',\n'SCA', 'X',\n'SCC', 'X',\n'SCD', 'X',\n'SCH', \
'C',\n'SCI', 'X',\n'SCL', 'X',\n'SCM', 'X',\n'SCN'\
, 'X',\n'SCO', 'X',\n'SCP', 'S',\n'SCR', 'X',\n'SC\
S', 'X',\n'SCV', 'C',\n'SCY', 'C',\n'SD8', 'X',\n'\
SDK', 'X',\n'SDZ', 'X',\n'SE4', 'X',\n'SEA', 'X',\\
n'SEB', 'S',\n'SEC', 'X',\n'SEG', 'A',\n'SEI', 'X'\
,\n'SEL', 'S',\n'SEM', 'X',\n'SEO', 'X',\n'SEP', '\
S',\n'SER', 'S',\n'SES', 'X',\n'SET', 'S',\n'SEU',\
 'X',\n'SF4', 'X',\n'SFG', 'X',\n'SFN', 'X',\n'SFO\
', 'X',\n'SGA', 'X',\n'SGC', 'X',\n'SGL', 'X',\n'S\
GM', 'X',\n'SGN', 'X',\n'SGP', 'X',\n'SHA', 'X',\n\
'SHC', 'X',\n'SHF', 'X',\n'SHH', 'X',\n'SHP', 'G',\
\n'SHR', 'E',\n'SHT', 'T',\n'SHU', 'X',\n'SI2', 'X\
',\n'SIA', 'X',\n'SIF', 'X',\n'SIG', 'X',\n'SIH', \
'X',\n'SIM', 'X',\n'SIN', 'X',\n'SKD', 'X',\n'SKF'\
, 'X',\n'SLB', 'X',\n'SLE', 'X',\n'SLZ', 'K',\n'SM\
A', 'X',\n'SMC', 'C',\n'SME', 'M',\n'SML', 'X',\n'\
SMM', 'M',\n'SMN', 'X',\n'SMP', 'X',\n'SMS', 'X',\\
n'SN1', 'X',\n'SN6', 'X',\n'SN7', 'X',\n'SNC', 'C'\
,\n'SNN', 'X',\n'SNP', 'X',\n'SO1', 'X',\n'SO2', '\
X',\n'SO3', 'X',\n'SO4', 'X',\n'SOA', 'X',\n'SOC',\
 'C',\n'SOM', 'X',\n'SOR', 'X',\n'SOT', 'X',\n'SOX\
', 'X',\n'SPA', 'X',\n'SPB', 'X',\n'SPC', 'X',\n'S\
PD', 'X',\n'SPE', 'X',\n'SPG', 'X',\n'SPH', 'X',\n\
'SPI', 'X',\n'SPK', 'X',\n'SPM', 'X',\n'SPN', 'X',\
\n'SPO', 'X',\n'SPP', 'X',\n'SPS', 'X',\n'SPY', 'X\
',\n'SQU', 'X',\n'SRA', 'X',\n'SRB', 'X',\n'SRD', \
'X',\n'SRL', 'X',\n'SRM', 'X',\n'SRS', 'X',\n'SRY'\
, 'X',\n'SSA', 'X',\n'SSB', 'X',\n'SSG', 'X',\n'SS\
P', 'X',\n'ST1', 'X',\n'ST2', 'X',\n'ST3', 'X',\n'\
ST4', 'X',\n'ST5', 'X',\n'ST6', 'X',\n'STA', 'X',\\
n'STB', 'X',\n'STE', 'X',\n'STG', 'X',\n'STI', 'X'\
,\n'STL', 'X',\n'STN', 'X',\n'STO', 'X',\n'STP', '\
X',\n'STR', 'X',\n'STU', 'X',\n'STY', 'Y',\n'SU1',\
 'X',\n'SU2', 'X',\n'SUC', 'X',\n'SUI', 'X',\n'SUL\
', 'X',\n'SUR', 'X',\n'SVA', 'S',\n'SWA', 'X',\n'T\
16', 'X',\n'T19', 'X',\n'T23', 'X',\n'T29', 'X',\n\
'T33', 'X',\n'T3P', 'X',\n'T42', 'A',\n'T44', 'X',\
\n'T5A', 'X',\n'T6A', 'T',\n'T6P', 'X',\n'T80', 'X\
',\n'T87', 'X',\n'TA1', 'X',\n'TAA', 'X',\n'TAB', \
'X',\n'TAC', 'X',\n'TAD', 'X',\n'TAF', 'X',\n'TAM'\
, 'X',\n'TAP', 'X',\n'TAR', 'X',\n'TAS', 'X',\n'TA\
U', 'X',\n'TAX', 'X',\n'TAZ', 'X',\n'TB9', 'X',\n'\
TBA', 'X',\n'TBD', 'X',\n'TBG', 'G',\n'TBH', 'X',\\
n'TBM', 'T',\n'TBO', 'X',\n'TBP', 'X',\n'TBR', 'X'\
,\n'TBS', 'X',\n'TBT', 'X',\n'TBU', 'X',\n'TBZ', '\
X',\n'TC4', 'X',\n'TCA', 'X',\n'TCB', 'X',\n'TCH',\
 'X',\n'TCK', 'X',\n'TCL', 'X',\n'TCM', 'X',\n'TCN\
', 'X',\n'TCP', 'X',\n'TCR', 'W',\n'TCS', 'X',\n'T\
CZ', 'X',\n'TDA', 'X',\n'TDB', 'X',\n'TDG', 'X',\n\
'TDP', 'X',\n'TDR', 'X',\n'TDX', 'X',\n'TEA', 'X',\
\n'TEM', 'X',\n'TEN', 'X',\n'TEO', 'X',\n'TEP', 'X\
',\n'TER', 'X',\n'TES', 'X',\n'TET', 'X',\n'TFA', \
'X',\n'TFB', 'X',\n'TFH', 'X',\n'TFI', 'X',\n'TFK'\
, 'X',\n'TFP', 'X',\n'THA', 'X',\n'THB', 'X',\n'TH\
C', 'T',\n'THD', 'X',\n'THE', 'X',\n'THF', 'X',\n'\
THJ', 'X',\n'THK', 'X',\n'THM', 'X',\n'THN', 'X',\\
n'THO', 'T',\n'THP', 'X',\n'THQ', 'X',\n'THR', 'T'\
,\n'THS', 'X',\n'THT', 'X',\n'THU', 'X',\n'THX', '\
X',\n'THZ', 'X',\n'TI1', 'X',\n'TI2', 'X',\n'TI3',\
 'P',\n'TIA', 'X',\n'TIH', 'A',\n'TK4', 'X',\n'TLA\
', 'X',\n'TLC', 'X',\n'TLM', 'X',\n'TLN', 'X',\n'T\
LX', 'X',\n'TM5', 'X',\n'TM6', 'X',\n'TMA', 'X',\n\
'TMB', 'T',\n'TMC', 'X',\n'TMD', 'T',\n'TME', 'X',\
\n'TMF', 'X',\n'TML', 'K',\n'TMM', 'X',\n'TMN', 'X\
',\n'TMP', 'X',\n'TMQ', 'X',\n'TMR', 'X',\n'TMT', \
'X',\n'TMZ', 'X',\n'TNB', 'C',\n'TND', 'X',\n'TNK'\
, 'X',\n'TNP', 'X',\n'TNT', 'X',\n'TOA', 'X',\n'TO\
B', 'X',\n'TOC', 'X',\n'TOL', 'X',\n'TOP', 'X',\n'\
TOS', 'X',\n'TOT', 'X',\n'TP1', 'G',\n'TP2', 'P',\\
n'TP3', 'E',\n'TP4', 'E',\n'TP7', 'T',\n'TPA', 'X'\
,\n'TPE', 'X',\n'TPF', 'X',\n'TPI', 'X',\n'TPL', '\
W',\n'TPM', 'X',\n'TPN', 'G',\n'TPO', 'T',\n'TPP',\
 'X',\n'TPQ', 'A',\n'TPR', 'P',\n'TPS', 'X',\n'TPT\
', 'X',\n'TPV', 'X',\n'TPX', 'X',\n'TPY', 'X',\n'T\
Q3', 'X',\n'TQ4', 'X',\n'TQ5', 'X',\n'TQ6', 'X',\n\
'TR1', 'X',\n'TRA', 'X',\n'TRB', 'X',\n'TRC', 'X',\
\n'TRD', 'X',\n'TRE', 'X',\n'TRF', 'W',\n'TRG', 'K\
',\n'TRH', 'X',\n'TRI', 'X',\n'TRJ', 'X',\n'TRM', \
'X',\n'TRN', 'W',\n'TRO', 'W',\n'TRP', 'W',\n'TRQ'\
, 'X',\n'TRS', 'X',\n'TRX', 'W',\n'TRZ', 'X',\n'TS\
2', 'X',\n'TS3', 'X',\n'TS4', 'X',\n'TS5', 'X',\n'\
TSA', 'X',\n'TSB', 'X',\n'TSI', 'X',\n'TSM', 'X',\\
n'TSN', 'X',\n'TSP', 'X',\n'TSU', 'X',\n'TTA', 'X'\
,\n'TTE', 'X',\n'TTN', 'X',\n'TTO', 'X',\n'TTP', '\
X',\n'TTX', 'X',\n'TXL', 'X',\n'TYA', 'Y',\n'TYB',\
 'Y',\n'TYD', 'X',\n'TYI', 'Y',\n'TYL', 'X',\n'TYM\
', 'W',\n'TYN', 'Y',\n'TYQ', 'Y',\n'TYR', 'Y',\n'T\
YS', 'Y',\n'TYV', 'X',\n'TYY', 'A',\n'TZB', 'X',\n\
'TZC', 'X',\n'TZE', 'X',\n'TZL', 'X',\n'TZO', 'X',\
\n'TZP', 'X',\n'U01', 'X',\n'U02', 'X',\n'U03', 'X\
',\n'U04', 'X',\n'U05', 'X',\n'U0E', 'X',\n'U10', \
'X',\n'U18', 'X',\n'U2G', 'X',\n'U3P', 'X',\n'U49'\
, 'X',\n'U55', 'X',\n'U5P', 'X',\n'U66', 'X',\n'U8\
9', 'X',\n'U8U', 'X',\n'UAA', 'X',\n'UAG', 'A',\n'\
UAP', 'X',\n'UAR', 'X',\n'UC1', 'X',\n'UC2', 'X',\\
n'UC3', 'X',\n'UC4', 'X',\n'UD1', 'X',\n'UD2', 'X'\
,\n'UDP', 'X',\n'UDX', 'X',\n'UFG', 'X',\n'UFM', '\
X',\n'UFP', 'X',\n'UGA', 'X',\n'UIN', 'X',\n'UKP',\
 'A',\n'UM3', 'X',\n'UMA', 'A',\n'UMG', 'X',\n'UMP\
', 'X',\n'UNA', 'X',\n'UND', 'X',\n'UNI', 'X',\n'U\
NK', 'X',\n'UNN', 'X',\n'UNX', 'X',\n'UP5', 'X',\n\
'UP6', 'X',\n'UPA', 'X',\n'UPF', 'X',\n'UPG', 'X',\
\n'UPP', 'X',\n'UQ1', 'X',\n'UQ2', 'X',\n'UQ6', 'X\
',\n'UR2', 'X',\n'URA', 'X',\n'URE', 'X',\n'URF', \
'X',\n'URI', 'X',\n'URS', 'X',\n'UTP', 'X',\n'UVC'\
, 'X',\n'UVW', 'X',\n'V35', 'X',\n'V36', 'X',\n'V4\
O', 'X',\n'V7O', 'X',\n'VAA', 'V',\n'VAC', 'X',\n'\
VAD', 'V',\n'VAF', 'V',\n'VAG', 'X',\n'VAL', 'V',\\
n'VAN', 'X',\n'VAS', 'X',\n'VAX', 'X',\n'VDX', 'X'\
,\n'VDY', 'X',\n'VG1', 'X',\n'VIB', 'X',\n'VIR', '\
X',\n'VIT', 'X',\n'VK3', 'X',\n'VO3', 'X',\n'VO4',\
 'X',\n'VS1', 'F',\n'VS2', 'F',\n'VS3', 'F',\n'VS4\
', 'F',\n'VXA', 'X',\n'W01', 'X',\n'W02', 'X',\n'W\
03', 'X',\n'W11', 'X',\n'W33', 'X',\n'W35', 'X',\n\
'W42', 'X',\n'W43', 'X',\n'W54', 'X',\n'W56', 'X',\
\n'W59', 'X',\n'W71', 'X',\n'W84', 'X',\n'W8R', 'X\
',\n'W91', 'X',\n'WAY', 'X',\n'WCC', 'X',\n'WO2', \
'X',\n'WO4', 'X',\n'WRB', 'X',\n'WRR', 'X',\n'WRS'\
, 'X',\n'WW7', 'X',\n'X2F', 'X',\n'X7O', 'X',\n'XA\
A', 'X',\n'XAN', 'X',\n'XAO', 'X',\n'XBB', 'X',\n'\
XBP', 'X',\n'XDN', 'X',\n'XDP', 'X',\n'XIF', 'X',\\
n'XIM', 'X',\n'XK2', 'X',\n'XL1', 'X',\n'XLS', 'X'\
,\n'XMP', 'X',\n'XN1', 'X',\n'XN2', 'X',\n'XN3', '\
X',\n'XUL', 'X',\n'XV6', 'X',\n'XYD', 'X',\n'XYH',\
 'X',\n'XYL', 'X',\n'XYP', 'X',\n'XYS', 'X',\n'YOF\
', 'Y',\n'YRR', 'X',\n'YT3', 'X',\n'YZ9', 'X',\n'Z\
34', 'G',\n'Z5A', 'X',\n'ZAF', 'X',\n'ZAP', 'X',\n\
'ZEB', 'X',\n'ZEN', 'X',\n'ZES', 'X',\n'ZID', 'X',\
\n'ZMR', 'X',\n'ZN3', 'X',\n'ZNH', 'X',\n'ZNO', 'X\
',\n'ZO3', 'X',\n'ZPR', 'P',\n'ZRA', 'A',\n'ZST', \
'X',\n'ZYA', 'A',\n\n\n'ASN','N');\n} \n\n\nsub fi\
le2head\n      {\n	my $file = shift;\n	my $size = \
shift;\n	my $f= new FileHandle;\n	my $line;\n	open\
 ($f,$file);\n	read ($f,$line, $size);\n	close ($f\
);\n	return $line;\n      }\nsub file2tail\n      \
{\n	my $file = shift;\n	my $size = shift;\n	my $f=\
 new FileHandle;\n	my $line;\n	\n	open ($f,$file);\
\n	seek ($f,$size*-1, 2);\n	read ($f,$line, $size)\
;\n	close ($f);\n	return $line;\n      }\n\n\nsub \
vtmpnam\n      {\n	my $r=rand(100000);\n	my $f=\"f\
ile.$r.$$\";\n	while (-e $f)\n	  {\n	    $f=vtmpna\
m();\n	  }\n	push (@TMPFILE_LIST, $f);\n	return $f\
;\n      }\n\nsub myexit\n  {\n    my $code=@_[0];\
\n    if ($CLEAN_EXIT_STARTED==1){return;}\n    el\
se {$CLEAN_EXIT_STARTED=1;}\n    ### ONLY BARE EXI\
T\n    exit ($code);\n  }\nsub set_error_lock\n   \
 {\n      my $name = shift;\n      my $pid=$$;\n\n\
      \n      &lock4tc ($$,\"LERROR\", \"LSET\", \\
"$$ -- ERROR: $name $PROGRAM\\n\");\n      return;\
\n    }\nsub set_lock\n  {\n    my $pid=shift;\n  \
  my $msg= shift;\n    my $p=getppid();\n    &lock\
4tc ($pid,\"LLOCK\",\"LRESET\",\"$p$msg\\n\");\n  \
}\nsub unset_lock\n   {\n     \n    my $pid=shift;\
\n    &lock4tc ($pid,\"LLOCK\",\"LRELEASE\",\"\");\
\n  }\nsub shift_lock\n  {\n    my $from=shift;\n \
   my $to=shift;\n    my $from_type=shift;\n    my\
 $to_type=shift;\n    my $action=shift;\n    my $m\
sg;\n    \n    if (!&lock4tc($from, $from_type, \"\
LCHECK\", \"\")){return 0;}\n    $msg=&lock4tc ($f\
rom, $from_type, \"LREAD\", \"\");\n    &lock4tc (\
$from, $from_type,\"LRELEASE\", $msg);\n    &lock4\
tc ($to, $to_type, $action, $msg);\n    return;\n \
 }\nsub isshellpid\n  {\n    my $p=shift;\n    if \
(!lock4tc ($p, \"LLOCK\", \"LCHECK\")){return 0;}\\
n    else\n      {\n	my $c=lock4tc($p, \"LLOCK\", \
\"LREAD\");\n	if ( $c=~/-SHELL-/){return 1;}\n    \
  }\n    return 0;\n  }\nsub isrootpid\n  {\n    i\
f(lock4tc (getppid(), \"LLOCK\", \"LCHECK\")){retu\
rn 0;}\n    else {return 1;}\n  }\nsub lock4tc\n	{\
\n	  my ($pid,$type,$action,$value)=@_;\n	  my $fn\
ame;\n	  my $host=hostname;\n	  \n	  if ($type eq \
\"LLOCK\"){$fname=\"$LOCKDIR/.$pid.$host.lock4tcof\
fee\";}\n	  elsif ( $type eq \"LERROR\"){ $fname=\\
"$LOCKDIR/.$pid.$host.error4tcoffee\";}\n	  elsif \
( $type eq \"LWARNING\"){ $fname=\"$LOCKDIR/.$pid.\
$host.warning4tcoffee\";}\n	  \n	  if ($debug_lock\
)\n	    {\n	      print STDERR \"\\n\\t---lock4tc(\
tcg): $action => $fname =>$value (RD: $LOCKDIR)\\n\
\";\n	    }\n\n	  if    ($action eq \"LCHECK\") {r\
eturn -e $fname;}\n	  elsif ($action eq \"LREAD\")\
{return file2string($fname);}\n	  elsif ($action e\
q \"LSET\") {return string2file ($value, $fname, \\
">>\");}\n	  elsif ($action eq \"LRESET\") {return\
 string2file ($value, $fname, \">\");}\n	  elsif (\
$action eq \"LRELEASE\") \n	    {\n	      if ( $de\
bug_lock)\n		{\n		  my $g=new FileHandle;\n		  ope\
n ($g, \">>$fname\");\n		  print $g \"\\nDestroyed\
 by $$\\n\";\n		  close ($g);\n		  safe_system (\"\
mv $fname $fname.old\");\n		}\n	      else\n		{\n	\
	  unlink ($fname);\n		}\n	    }\n	  return \"\";\\
n	}\n	\nsub file2string\n	{\n	  my $file=@_[0];\n	\
  my $f=new FileHandle;\n	  my $r;\n	  open ($f, \\
"$file\");\n	  while (<$f>){$r.=$_;}\n	  close ($f\
);\n	  return $r;\n	}\nsub string2file \n    {\n  \
  my ($s,$file,$mode)=@_;\n    my $f=new FileHandl\
e;\n    \n    open ($f, \"$mode$file\");\n    prin\
t $f  \"$s\";\n    close ($f);\n  }\n\nBEGIN\n    \
{\n      srand;\n    \n      $SIG{'SIGUP'}='signal\
_cleanup';\n      $SIG{'SIGINT'}='signal_cleanup';\
\n      $SIG{'SIGQUIT'}='signal_cleanup';\n      $\
SIG{'SIGILL'}='signal_cleanup';\n      $SIG{'SIGTR\
AP'}='signal_cleanup';\n      $SIG{'SIGABRT'}='sig\
nal_cleanup';\n      $SIG{'SIGEMT'}='signal_cleanu\
p';\n      $SIG{'SIGFPE'}='signal_cleanup';\n     \
 \n      $SIG{'SIGKILL'}='signal_cleanup';\n      \
$SIG{'SIGPIPE'}='signal_cleanup';\n      $SIG{'SIG\
STOP'}='signal_cleanup';\n      $SIG{'SIGTTIN'}='s\
ignal_cleanup';\n      $SIG{'SIGXFSZ'}='signal_cle\
anup';\n      $SIG{'SIGINFO'}='signal_cleanup';\n \
     \n      $SIG{'SIGBUS'}='signal_cleanup';\n   \
   $SIG{'SIGALRM'}='signal_cleanup';\n      $SIG{'\
SIGTSTP'}='signal_cleanup';\n      $SIG{'SIGTTOU'}\
='signal_cleanup';\n      $SIG{'SIGVTALRM'}='signa\
l_cleanup';\n      $SIG{'SIGUSR1'}='signal_cleanup\
';\n\n\n      $SIG{'SIGSEGV'}='signal_cleanup';\n \
     $SIG{'SIGTERM'}='signal_cleanup';\n      $SIG\
{'SIGCONT'}='signal_cleanup';\n      $SIG{'SIGIO'}\
='signal_cleanup';\n      $SIG{'SIGPROF'}='signal_\
cleanup';\n      $SIG{'SIGUSR2'}='signal_cleanup';\
\n\n      $SIG{'SIGSYS'}='signal_cleanup';\n      \
$SIG{'SIGURG'}='signal_cleanup';\n      $SIG{'SIGC\
HLD'}='signal_cleanup';\n      $SIG{'SIGXCPU'}='si\
gnal_cleanup';\n      $SIG{'SIGWINCH'}='signal_cle\
anup';\n      \n      $SIG{'INT'}='signal_cleanup'\
;\n      $SIG{'TERM'}='signal_cleanup';\n      $SI\
G{'KILL'}='signal_cleanup';\n      $SIG{'QUIT'}='s\
ignal_cleanup';\n      \n      our $debug_lock=$EN\
V{\"DEBUG_LOCK\"};\n      \n      \n      \n      \
\n      foreach my $a (@ARGV){$CL.=\" $a\";}\n    \
  if ( $debug_lock ){print STDERR \"\\n\\n\\n*****\
***** START PG: $PROGRAM *************\\n\";}\n   \
   if ( $debug_lock ){print STDERR \"\\n\\n\\n****\
******(tcg) LOCKDIR: $LOCKDIR $$ *************\\n\\
";}\n      if ( $debug_lock ){print STDERR \"\\n -\
-- $$ -- $CL\\n\";}\n      \n	     \n      \n     \
 \n    }\nsub flush_error\n  {\n    my $msg=shift;\
\n    return add_error ($EXIT_FAILURE,$$, $$,getpp\
id(), $msg, $CL);\n  }\nsub add_error \n  {\n    m\
y $code=shift;\n    my $rpid=shift;\n    my $pid=s\
hift;\n    my $ppid=shift;\n    my $type=shift;\n \
   my $com=shift;\n    \n    $ERROR_DONE=1;\n    l\
ock4tc ($rpid, \"LERROR\",\"LSET\",\"$pid -- ERROR\
: $type\\n\");\n    lock4tc ($$, \"LERROR\",\"LSET\
\", \"$pid -- COM: $com\\n\");\n    lock4tc ($$, \\
"LERROR\",\"LSET\", \"$pid -- STACK: $ppid -> $pid\
\\n\");\n   \n    return $code;\n  }\nsub add_warn\
ing \n  {\n    my $rpid=shift;\n    my $pid =shift\
;\n    my $command=shift;\n    my $msg=\"$$ -- WAR\
NING: $command\\n\";\n    print STDERR \"$msg\";\n\
    lock4tc ($$, \"LWARNING\", \"LSET\", $msg);\n \
 }\n\nsub signal_cleanup\n  {\n    print dtderr \"\
\\n**** $$ (tcg) was killed\\n\";\n    &cleanup;\n\
    exit ($EXIT_FAILURE);\n  }\nsub clean_dir\n  {\
\n    my $dir=@_[0];\n    if ( !-d $dir){return ;}\
\n    elsif (!($dir=~/tmp/)){return ;}#safety chec\
k 1\n    elsif (($dir=~/\\*/)){return ;}#safety ch\
eck 2\n    else\n      {\n	`rm -rf $dir`;\n      }\
\n    return;\n  }\nsub cleanup\n  {\n    #print s\
tderr \"\\n----tc: $$ Kills $PIDCHILD\\n\";\n    #\
kill (SIGTERM,$PIDCHILD);\n    my $p=getppid();\n \
   $CLEAN_EXIT_STARTED=1;\n    \n    \n    \n    i\
f (&lock4tc($$,\"LERROR\", \"LCHECK\", \"\"))\n   \
   {\n	my $ppid=getppid();\n	if (!$ERROR_DONE) \n	\
  {\n	    &lock4tc($$,\"LERROR\", \"LSET\", \"$$ -\
- STACK: $p -> $$\\n\");\n	    &lock4tc($$,\"LERRO\
R\", \"LSET\", \"$$ -- COM: $CL\\n\");\n	  }\n    \
  }\n    my $warning=&lock4tc($$, \"LWARNING\", \"\
LREAD\", \"\");\n    my $error=&lock4tc($$,  \"LER\
ROR\", \"LREAD\", \"\");\n    #release error and w\
arning lock if root\n    \n    if (isrootpid() && \
($warning || $error) )\n      {\n	\n	print STDERR \
\"**************** Summary *************\\n$error\\
\n$warning\\n\";\n\n	&lock4tc($$,\"LERROR\",\"RELE\
ASE\",\"\");\n	&lock4tc($$,\"LWARNING\",\"RELEASE\\
",\"\");\n      } \n    \n    \n    foreach my $f \
(@TMPFILE_LIST)\n      {\n	if (-e $f){unlink ($f);\
} \n      }\n    foreach my $d (@TMPDIR_LIST)\n   \
   {\n	clean_dir ($d);\n      }\n    #No More Lock\
 Release\n    #&lock4tc($$,\"LLOCK\",\"LRELEASE\",\
\"\"); #release lock \n\n    if ( $debug_lock ){pr\
int STDERR \"\\n\\n\\n********** END PG: $PROGRAM \
($$) *************\\n\";}\n    if ( $debug_lock ){\
print STDERR \"\\n\\n\\n**********(tcg) LOCKDIR: $\
LOCKDIR $$ *************\\n\";}\n  }\nEND \n  {\n \
   \n    &cleanup();\n  }\n   \n\nsub safe_system \
\n{\n  my $com=shift;\n  my $ntry=shift;\n  my $ct\
ry=shift;\n  my $pid;\n  my $status;\n  my $ppid=g\
etppid();\n  if ($com eq \"\"){return 1;}\n  \n  \\
n\n  if (($pid = fork ()) < 0){return (-1);}\n  if\
 ($pid == 0)\n    {\n      set_lock($$, \" -SHELL-\
 $com (tcg)\");\n      exec ($com);\n    }\n  else\
\n    {\n      lock4tc ($$, \"LLOCK\", \"LSET\", \\
"$pid\\n\");#update parent\n      $PIDCHILD=$pid;\\
n    }\n  if ($debug_lock){printf STDERR \"\\n\\t \
.... safe_system (fasta_seq2hmm)  p: $$ c: $pid CO\
M: $com\\n\";}\n\n  waitpid ($pid,WTERMSIG);\n\n  \
shift_lock ($pid,$$, \"LWARNING\",\"LWARNING\", \"\
LSET\");\n\n  if ($? == $EXIT_FAILURE || lock4tc($\
pid, \"LERROR\", \"LCHECK\", \"\"))\n    {\n      \
if ($ntry && $ctry <$ntry)\n	{\n	  add_warning ($$\
,$$,\"$com failed [retry: $ctry]\");\n	  lock4tc (\
$pid, \"LRELEASE\", \"LERROR\", \"\");\n	  return \
safe_system ($com, $ntry, ++$ctry);\n	}\n      els\
if ($ntry == -1)\n	{\n	  if (!shift_lock ($pid, $$\
, \"LERROR\", \"LWARNING\", \"LSET\"))\n	    {\n	 \
     add_warning ($$,$$,\"$com failed\");\n	    }\\
n	  else\n	    {\n	      lock4tc ($pid, \"LRELEASE\
\", \"LERROR\", \"\");\n	    }\n	  return $?;}\n  \
    else\n	{\n	  if (!shift_lock ($pid,$$, \"LERRO\
R\",\"LERROR\", \"LSET\"))\n	    {\n	      myexit(\
add_error ($EXIT_FAILURE,$$,$pid,getppid(), \"UNSP\
ECIFIED system\", $com));\n	    }\n	}\n    }\n  re\
turn $?;\n}\n\nsub check_configuration \n    {\n  \
    my @l=@_;\n      my $v;\n      foreach my $p (\
@l)\n	{\n	  \n	  if   ( $p eq \"EMAIL\")\n	    { \\
n	      if ( !($EMAIL=~/@/))\n		{\n		add_warning($\
$,$$,\"Could Not Use EMAIL\");\n		myexit(add_error\
 ($EXIT_FAILURE,$$,$$,getppid(),\"EMAIL\",\"$CL\")\
);\n	      }\n	    }\n	  elsif( $p eq \"INTERNET\"\
)\n	    {\n	      if ( !&check_internet_connection\
())\n		{\n		  myexit(add_error ($EXIT_FAILURE,$$,$\
$,getppid(),\"INTERNET\",\"$CL\"));\n		}\n	    }\n\
	  elsif( $p eq \"wget\")\n	    {\n	      if (!&pg\
_is_installed (\"wget\") && !&pg_is_installed (\"c\
url\"))\n		{\n		  myexit(add_error ($EXIT_FAILURE,\
$$,$$,getppid(),\"PG_NOT_INSTALLED:wget\",\"$CL\")\
);\n		}\n	    }\n	  elsif( !(&pg_is_installed ($p)\
))\n	    {\n	      myexit(add_error ($EXIT_FAILURE\
,$$,$$,getppid(),\"PG_NOT_INSTALLED:$p\",\"$CL\"))\
;\n	    }\n	}\n      return 1;\n    }\nsub pg_is_i\
nstalled\n  {\n    my @ml=@_;\n    my $r, $p, $m;\\
n    my $supported=0;\n    \n    my $p=shift (@ml)\
;\n    if ($p=~/::/)\n      {\n	if (safe_system (\\
"perl -M$p -e 1\")==$EXIT_SUCCESS){return 1;}\n	el\
se {return 0;}\n      }\n    else\n      {\n	$r=`w\
hich $p 2>/dev/null`;\n	if ($r eq \"\"){return 0;}\
\n	else {return 1;}\n      }\n  }\n\n\n\nsub check\
_internet_connection\n  {\n    my $internet;\n    \
my $tmp;\n    &check_configuration ( \"wget\"); \n\
    \n    $tmp=&vtmpnam ();\n    \n    if     (&pg\
_is_installed    (\"wget\")){`wget www.google.com \
-O$tmp >/dev/null 2>/dev/null`;}\n    elsif  (&pg_\
is_installed    (\"curl\")){`curl www.google.com -\
o$tmp >/dev/null 2>/dev/null`;}\n    \n    if ( !-\
e $tmp || -s $tmp < 10){$internet=0;}\n    else {$\
internet=1;}\n    if (-e $tmp){unlink $tmp;}\n\n  \
  return $internet;\n  }\nsub check_pg_is_installe\
d\n  {\n    my @ml=@_;\n    my $r=&pg_is_installed\
 (@ml);\n    if (!$r && $p=~/::/)\n      {\n	print\
 STDERR \"\\nYou Must Install the perl package $p \
on your system.\\nRUN:\\n\\tsudo perl -MCPAN -e 'i\
nstall $pg'\\n\";\n      }\n    elsif (!$r)\n     \
 {\n	myexit(flush_error(\"\\nProgram $p Supported \
but Not Installed on your system\"));\n      }\n  \
  else\n      {\n	return 1;\n      }\n  }\n\n\n\n"\
,"use Cwd;\nuse Env;\nuse File::Path;\nuse FileHan\
dle;\nuse strict;\n\n\nour (%MODE, %PG, %ENV_SET, \
%SUPPORTED_OS);\n\n\nour $EXIT_SUCCESS=0;\nour $EX\
IT_FAILURE=1;\nour $INTERNET=0;\n\nour $CP=\"cp \"\
; #was causing a crash on MacOSX\nour $SILENT=\">/\
dev/null 2>/dev/null\";\nour $WEB_BASE=\"http://ww\
w.tcoffee.org\";\nour $TCLINKDB_ADDRESS=\"$WEB_BAS\
E/Resources/tclinkdb.txt\";\nour $OS=get_os();\nou\
r $ROOT=&get_root();\nour $CD=cwd();\nour $CDIR=$C\
D;\nour $HOME=$ENV{'HOME'};\n\nour $OSNAME=$ENV{'O\
SNAME'};\nour $OSARCH=$ENV{'OSARCH'};\nour $REPO_R\
OOT=\"\";\n\nour $TCDIR;\nour $TCCACHE;\nour $TCTM\
P;\nour $TCM;\nour $TCMETHODS;\nour $TCPLUGINS;\no\
ur $PLUGINS_DIR=\"\";\nour $INSTALL_DIR=\"\";\n\no\
ur $CXX=\"g++\";\nour $CXXFLAGS=\"\";\n\nour $CPP=\
\"g++\";\nour $CPPFLAGS=\"\";\n\nour $CC=\"gcc\";\\
nour $CFLAGS=\"\";\n\nour $FC=\"f77\";\nour $FFLAG\
S=\"\";\n\nmy $install=\"all\";\nmy $default_updat\
e_action=\"no_update\";\nmy @required_applications\
=(\"wget_OR_curl\");\nmy @smode=(\"all\", \"clean\\
", \"install\");\n\n&initialize_PG();\n\nmy $cl=jo\
in( \" \", @ARGV);\nif ($#ARGV==-1 || ($cl=~/-h/) \
||($cl=~/-H/) )\n  {\n     print \"\\n!!!!!!! ./in\
stall  t_coffee             --> installs t_coffee \
only\";\n     print \"\\n!!!!!!! ./install  all   \
               --> installs all the modes [mcoffee\
, expresso, psicoffee,rcoffee..]\";\n     print \"\
\\n!!!!!!! ./install  [mcoffee|rcoffee|..] --> ins\
talls the specified mode\";\n     print \"\\n!!!!!\
!! ./install  -h                   --> print usage\
\\n\\n\";\n     if ( $#ARGV==-1){exit ($EXIT_FAILU\
RE);}\n   }\n     \nif (($cl=~/-h/) ||($cl=~/-H/) \
)\n  {\n    my $m;\n    print \"\\n\\n!!!!!!! adva\
nced mode\\n\";\n    foreach $m ((keys (%MODE)),@s\
mode)\n      {\n	print \"!!!!!!!       ./install $\
m\\n\";\n      }\n    \n    print \"!!!!!!! ./inst\
all [target:package|mode|] [-update|-force|-exec=d\
ir|-dis=dir|-root|-tclinkdb=file|-] [CC=|FCC=|CXX=\
|CFLAGS=|CXXFLAGS=]\\n\";\n    print \"!!!!!!! ./i\
nstall clean    [removes all executables]\\n\";\n \
   print \"!!!!!!! ./install [optional:target] -up\
date               [updates package already instal\
led]\\n\";\n    print \"!!!!!!! ./install [optiona\
l:target] -force                [Forces recompilat\
ion over everything]\\n\";\n    \n    print \"!!!!\
!!! ./install [optional:target] -root             \
    [You are running as root]\\n\";\n    print \"!\
!!!!!! ./install [optional:target] -exec=/foo/bar/\
       [address for the T-Coffee executable]\\n\";\
\n    print \"!!!!!!! ./install [optional:target] \
-dis=/foo/bar/        [Address where distributions\
 should be stored]\\n\";\n    print \"!!!!!!! ./in\
stall [optional:target] -tclinkdb=foo|update  [fil\
e containing all the packages to be installed]\\n\\
";\n    print \"!!!!!!! ./install [optional:target\
] -clean                [clean everything]\\n\";\n\
    print \"!!!!!!! ./install [optional:target] -p\
lugins              [plugins directory]\\n\";\n   \
 print \"!!!!!!! ./install [optional:target] -tcdi\
r=/foor/bar      [base path where T-Coffee will be\
 installed]\\n\";\n    print \"!!!!!!! ./install [\
optional:target] -repo=/path/to/repo   [binaries r\
epository root directory]\\n\";\n    print \"!!!!!\
!! mode:\";\n    foreach $m (keys(%MODE)){print \"\
$m \";}\n    print \"\\n\";\n    print \"!!!!!!! P\
ackages:\";\n    foreach $m (keys (%PG)){print \"$\
m \";}\n    print \"\\n\";\n    \n    print \"\\n\\
\n\";\n    exit ($EXIT_FAILURE);\n  }\n\n\n\nmy (@\
argl)=($cl=~/(\\S+=[^=]+)\\s\\w+=/g);\npush (@argl\
, ($cl=~/(\\S+=[^=]+\\S)\\s*$/g));\n\nforeach $a (\
@argl)\n  {\n    if ( ($cl=~/CXX=(.*)/)){$CXX=$1;}\
\n    if ( ($cl=~/-CC=(.*)/    )){$CC=$1;}\n    if\
 ( ($cl=~/-FC=(.*)/    )){$FC=$1;}\n    if ( ($cl=\
~/-CFLAGS=(.*)/)){$CFLAGS=$1;}\n    if ( ($cl=~/-C\
XXFLAGS=(.*)/)){$CXXFLAGS=$1;}\n  }\nour ($ROOT_IN\
STALL, $NO_QUESTION, $default_update_action,$BINAR\
IES_ONLY,$force, $default_update_action, $INSTALL_\
DIR, $PLUGINS_DIR, $DISTRIBUTIONS,$tclinkdb, $prox\
y, $clean);\nif ( ($cl=~/-root/)){$ROOT_INSTALL=1;\
}\nif ( ($cl=~/-no_question/)){$NO_QUESTION=1;}\ni\
f ( ($cl=~/-update/)){$default_update_action=\"upd\
ate\";}\nif ( ($cl=~/-binaries/)){$BINARIES_ONLY=1\
;}\nif ( ($cl=~/-force/)){$force=1;$default_update\
_action=\"update\"}\nif ( ($cl=~/-exec=\\s*(\\S+)/\
)){$INSTALL_DIR=$1;}\nif ( ($cl=~/-plugins=\\s*(\\\
S+)/)){$PLUGINS_DIR=$1;}\nif ( ($cl=~/-dis=\\s*(\\\
S+)/)){$DISTRIBUTIONS=$1;}\n\nif ( ($cl=~/-tclinkd\
b=\\s*(\\S+)/)){$tclinkdb=$1;}\nif ( ($cl=~/-proxy\
=\\s*(\\S+)/)){$proxy=$1;}\nif ( ($cl=~/-clean/)){\
$clean=1;}\nif ( ($cl=~/-repo=\\s*(\\S+)/)){ $REPO\
_ROOT=$1; }\nif ( ($cl=~/-tcdir=\\s*(\\S+)/)){ $TC\
DIR=$1; }\nif ($tclinkdb){&update_tclinkdb ($tclin\
kdb);}\n\n\nif( $REPO_ROOT ne \"\" ) {\n	if( $OSNA\
ME eq \"\" ) { print \"You have specified the repo\
sitory folder but the required \\\"OSNAME\\\" envi\
roment variable is missing. \\n\"; exit 1; } \n	if\
( $OSARCH eq \"\" ) { print \"You have specified t\
he repository folder but the required \\\"OSARCH\\\
\" enviroment variable is missing. \\n\"; exit 1; \
} \n}\n\n\nif(!$TCDIR) { $TCDIR=\"$HOME/.t_coffee\\
"; }\n&add_dir ($TCDIR);\n&add_dir ($TCCACHE=\"$TC\
DIR/cache\");\n&add_dir ($TCTMP=\"$CDIR/tmp\");\n&\
add_dir ($TCM=\"$TCDIR/mcoffee\");\n&add_dir ($TCM\
ETHODS=\"$TCDIR/methods\");\n&add_dir ($TCPLUGINS=\
\"$TCDIR/plugins/$OS\");\n\n\nour $BASE=\"$CD/bin\\
";\nour $BIN=\"$BASE/binaries/$OS\";\nour $DOWNLOA\
D_DIR=\"$BASE/download\";\nour $DOWNLOAD_FILE=\"$D\
OWNLOAD_DIR/files\";\nour $TMP=\"$BASE/tmp\";\n\n&\
add_dir($BASE);\n&add_dir($BIN);\n&add_dir($DOWNLO\
AD_DIR);\n&add_dir($DOWNLOAD_FILE);\nif (!$DISTRIB\
UTIONS){$DISTRIBUTIONS=\"$DOWNLOAD_DIR/distributio\
ns\";}\n&add_dir ($DISTRIBUTIONS);\n&add_dir ($TMP\
);\n\n\nif    (!$PLUGINS_DIR && !$ROOT_INSTALL){$P\
LUGINS_DIR=$TCPLUGINS;}\nelsif (!$PLUGINS_DIR &&  \
$ROOT_INSTALL){$PLUGINS_DIR=\"/usr/local/bin/\";}\\
n\nif    (!$INSTALL_DIR && !$ROOT_INSTALL){$INSTAL\
L_DIR=\"$HOME/bin/\";mkpath ($INSTALL_DIR);}\nelsi\
f (!$INSTALL_DIR &&  $ROOT_INSTALL){$INSTALL_DIR=\\
"/usr/local/bin/\";}\n\nif (-d \"mcoffee\"){`cp mc\
offee/* $TCM`;}\n\n\nour $ENV_FILE=\"$TCDIR/t_coff\
ee_env\";\n&env_file2putenv ($ENV_FILE);\n&set_pro\
xy($proxy);\nmy ($target, $p, $r);\n$target=$p;\n\\
nforeach $p (  ((keys (%PG)),(keys(%MODE)),(@smode\
)) )\n  {\n    if ($ARGV[0] eq $p && $target eq \"\
\"){$target=$p;}\n  }\nif ($target eq \"\"){exit (\
$EXIT_FAILURE);}\n\n\nforeach $r (@required_applic\
ations)\n  {\n    my @app_list;\n    my $i;\n    $\
i=0;\n    \n    @app_list=split (/_OR_/, $r);\n   \
 foreach my $pg (@app_list)\n      {\n	$i+=&pg_is_\
installed ($pg);\n      }\n    if ($i==0)\n      {\
\n      print \"One of the following packages must\
 be installed to proceed: \";\n      foreach my $p\
g (@app_list)\n	{\n	  print (\"$pg \");\n	}\n     \
 die;\n    }\n  }\n\n\n\n\n\n\n&sign_license_ni();\
\n\n\n$PG{C}{compiler}=get_C_compiler($CC);\n$PG{F\
ortran}{compiler}=get_F_compiler($FC);\n$PG{CXX}{c\
ompiler}=$PG{CPP}{compiler}=$PG{GPP}{compiler}=get\
_CXX_compiler($CXX);\nif ($CXXFLAGS){$PG{CPP}{opti\
ons}=$PG{GPP}{options}=$PG{CXX}{options}=$CXXFLAGS\
;}\nif ($CFLAGS){$PG{C}{options}=$CFLAGS;}\nforeac\
h my $c (keys(%PG))\n  {\n    my $arguments;\n    \
if ($PG{$c}{compiler})\n      {\n	$arguments=\"$PG\
{$c}{compiler_flag}=$PG{$c}{compiler} \";\n	if ($P\
G{$c}{options})\n	  {\n	    $arguments.=\"$PG{$c}{\
options_flag}=$PG{$c}{options} \";\n	  }\n	$PG{$c}\
{arguments}=$arguments;\n      }\n  }\n\nif ($PG{$\
target}){$PG{$target}{install}=1;}\nelse\n  {\n   \
 foreach my $pg (keys(%PG))\n      {\n	if ( $targe\
t eq \"all\" || ($PG{$pg}{mode}=~/$target/))\n	  {\
\n	    $PG{$pg} {install}=1;\n	  }\n      }\n  }\n\
\nforeach my $pg (keys(%PG))\n  {\n    if (!$PG{$p\
g}{update_action}){$PG{$pg}{update_action}=$defaul\
t_update_action;}\n    elsif ($PG{$pg}{update_acti\
on} eq \"never\"){$PG{$pg}{install}=0;}\n    if ( \
$force && $PG{$pg}{install})\n      {\n	`rm $BIN/$\
pg $BIN/$pg.exe $SILENT`;\n      }\n    if ($PG{$p\
g}{update_action} eq \"update\" && $PG{$pg}{instal\
l}){$PG{$pg}{update}=1;}\n  }\n\nif (($target=~/cl\
ean/))\n  {\n    print \"------- cleaning executab\
les -----\\n\";\n    `rm bin/* $SILENT`;\n    exit\
 ($EXIT_SUCCESS);\n  }\n\nif ( !$PG{$target}){prin\
t \"------- Installing T-Coffee Modes\\n\";}\n\nfo\
reach my $m (keys(%MODE))\n  {\n    if ( $target e\
q \"all\" || $target eq $m)\n      {\n	print \"\\n\
------- The installer will now install the $m comp\
onents $MODE{$m}{description}\\n\";\n	foreach my $\
pg (keys(%PG))\n	  {\n	    if ( $PG{$pg}{mode} =~/\
$m/ && $PG{$pg}{install})\n	      {\n		if ($PG{$pg\
}{touched}){print \"------- $PG{$pg}{dname}: alrea\
dy processed\\n\";}\n		else {$PG{$pg}{success}=&in\
stall_pg($pg);$PG{$pg}{touched}=1;}\n	      }\n	  \
}\n      }\n  }\n\nif ( $PG{$target}){print \"----\
--- Installing Individual Package\\n\";}\nforeach \
my $pg (keys (%PG))\n  {\n    \n    if ( $PG{$pg}{\
install} && !$PG{$pg}{touched})\n      {\n	print \\
"\\n------- Install $pg\\n\";\n	$PG{$pg}{success}=\
&install_pg($pg);$PG{$pg}{touched}=1;\n      }\n  \
}\nprint \"------- Finishing The installation\\n\"\
;\nmy $final_report=&install ($INSTALL_DIR);\n\npr\
int \"\\n\";\nprint \"****************************\
*****************************************\\n\";\np\
rint \"********              INSTALLATION SUMMARY \
         *****************\\n\";\nprint \"********\
**************************************************\
***********\\n\";\nprint \"------- SUMMARY package\
 Installation:\\n\";\nprint \"-------   Executable\
 Installed in: $PLUGINS_DIR\\n\";\n\nforeach my $p\
g (keys(%PG))\n  {\n    if ( $PG{$pg}{install})\n \
     {\n	my $bin_status=($PG{$pg}{from_binary} && \
$PG{$pg}{success})?\"[from binary]\":\"\";\n	if   \
  ( $PG{$pg}{new} && !$PG{$pg}{old})              \
       {print \"*------        $PG{$pg}{dname}: in\
stalled $bin_status\\n\"; $PG{$pg}{status}=1;}\n	e\
lsif  ( $PG{$pg}{new} &&  $PG{$pg}{old})          \
           {print \"*------        $PG{$pg}{dname}\
: updated $bin_status\\n\"  ; $PG{$pg}{status}=1;}\
 \n	elsif  (!$PG{$pg}{new} &&  $PG{$pg}{old} && !$\
PG{$pg}{update}){print \"*------        $PG{$pg}{d\
name}: previous\\n\" ; $PG{$pg}{status}=1;}\n	elsi\
f  (!$PG{$pg}{new} &&  $PG{$pg}{old} &&  $PG{$pg}{\
update}){print \"*------        $PG{$pg}{dname}: f\
ailed update (previous installation available)\\n\\
";$PG{$pg}{status}=0;}\n	else                     \
                                     {print \"*---\
---        $PG{$pg}{dname}: failed installation\\n\
\";$PG{$pg}{status}=0;}\n      }\n  }\nmy $failure\
;\n\nif ( !$PG{$target}){print \"*------ SUMMARY m\
ode Installation:\\n\";}\nforeach my $m (keys(%MOD\
E))\n  {\n  \n    if ( $target eq \"all\" || $targ\
et eq $m)\n      {\n	my $succesful=1;\n	foreach my\
 $pg (keys(%PG))\n	  {\n	    if (($PG{$pg}{mode}=~\
/$m/) && $PG{$pg}{install} && $PG{$pg}{status}==0)\
\n	      {\n		$succesful=0;\n		print \"*!!!!!!    \
   $PG{$pg}{dname}: Missing\\n\";\n	      }\n	  }\\
n	if ( $succesful)\n	  {\n	    $MODE{$m}{status}=1\
;\n	    print \"*------       MODE $MODE{$m}{dname\
} SUCCESSFULLY installed\\n\";\n	  }\n	else\n	  {\\
n	    $failure++;\n	    $MODE{$m}{status}=0;\n	   \
 print \"*!!!!!!       MODE $MODE{$m}{dname} UNSUC\
CESSFULLY installed\\n\";\n	  }\n      }\n  }\n\n \
   \n      \nif ($clean==1 && ($BASE=~/install4tco\
ffee/) ){print \"*------ Clean Installation Direct\
ory: $BASE\\n\";`rm -rf $BASE`;}\nforeach my $pg (\
keys(%PG)){if ($PG{$pg}{install} && $PG{$pg}{statu\
s}==0){exit ($EXIT_FAILURE);}}\n\nif ($failure)\n \
 {\n    print \"**********************************\
***********************************\\n\";\n    pri\
nt \"********     SOME PACKAGES FAILED TO INSTALL \
       *****************\\n\";\n    print \"******\
**************************************************\
*************\\n\";\n    print \"\\nSome of the re\
ported failures may be due to connectivity problem\
s\";\n    print \"\\nRerun the installation and th\
e installer will specifically try to install the m\
issing packages\";\n    print \"\\nIf this Fails, \
go to the original website and install the package\
 manually\";\n  }\n\nprint \"*********************\
************************************************\\\
n\";\nprint \"********              FINALIZE YOUR \
INSTALLATION    *****************\\n\";\nprint \"*\
**************************************************\
******************\\n\";\nprint \"------- Your exe\
cutables are in:\\n\"; \nprint \"-------       $PL\
UGINS_DIR:\\n\";\nprint \"------- Add this directo\
ry to your path with the following command:\\n\";\\
nprint \"-------       export PATH=$PLUGINS_DIR:\\\
$PATH\\n\";\nprint \"------- Make this permanent b\
y adding this line to the file:\\n\";\nprint \"---\
----       $HOME/.bashrc\\n\";\nexit ($EXIT_SUCCES\
S);  \n  \nsub get_CXX_compiler\n  {\n    my $c=@_\
[0];\n    my (@clist)=(\"g++\");\n    \n    return\
 get_compil ($c, @clist);\n }\nsub get_C_compiler\\
n  {\n    my $c=@_[0];\n    my (@clist)=(\"gcc\", \
\"cc\", \"icc\");\n    \n    return get_compil ($c\
, @clist);\n }\n\nsub get_F_compiler\n  {\n    my \
($c)=@_[0];\n    my @clist=(\"f77\", \"g77\",\"g95\
\", \"gfortran\", \"ifort\");\n    return get_comp\
il ($c, @clist);\n  } \n       \nsub get_compil\n \
 {\n    my ($fav,@clist)=(@_);\n    \n    #return \
the first compiler found installed in the system. \
Check first the favorite\n    foreach my $c ($fav,\
@clist)\n      {\n	if  (&pg_is_installed ($c)){ret\
urn $c;}\n      }\n    return \"\";\n  }\nsub exit\
_if_pg_not_installed\n  {\n    my (@arg)=(@_);\n  \
  \n    foreach my $p (@arg)\n      {\n	if ( !&pg_\
is_installed ($p))\n	  {\n	    print \"!!!!!!!! Th\
e $p utility must be installed for this installati\
on to proceed [FATAL]\\n\";\n	    die;\n	  }\n    \
  }\n    return 1;\n  }\nsub set_proxy\n  {\n    m\
y ($proxy)=(@_);\n    my (@list,$p);\n    \n    @l\
ist= (\"HTTP_proxy\", \"http_proxy\", \"HTTP_PROXY\
\", \"ALL_proxy\", \"all_proxy\",\"HTTP_proxy_4_TC\
OFFEE\",\"http_proxy_4_TCOFFEE\");\n    \n    if (\
!$proxy)\n      {\n	foreach my $p (@list)\n	  {\n	\
    if ( ($ENV_SET{$p}) || $ENV{$p}){$proxy=$ENV{$\
p};}\n	  }\n      }\n    foreach my $p(@list){$ENV\
{$p}=$proxy;}\n  }\n	\nsub check_internet_connecti\
on\n  {\n    my $internet;\n    \n    if ( -e \"x\\
"){unlink (\"x\");}\n    if     (&pg_is_installed \
   (\"wget\")){`wget www.google.com -Ox >/dev/null\
 2>/dev/null`;}\n    elsif  (&pg_is_installed    (\
\"curl\")){`curl www.google.com -ox >/dev/null 2>/\
dev/null`;}\n    else\n      {\n	printf stderr \"\\
\nERROR: No pg for remote file fetching [wget or c\
url][FATAL]\\n\";\n	exit ($EXIT_FAILURE);\n      }\
\n    \n    if ( !-e \"x\" || -s \"x\" < 10){$inte\
rnet=0;}\n    else {$internet=1;}\n    if (-e \"x\\
"){unlink \"x\";}\n    return $internet;\n  }\nsub\
 url2file\n  {\n    my ($cmd, $file,$wget_arg, $cu\
rl_arg)=(@_);\n    my ($exit,$flag, $pg, $arg);\n \
   \n    if ($INTERNET || check_internet_connectio\
n ()){$INTERNET=1;}\n    else\n      {\n	print STD\
ERR \"ERROR: No Internet Connection [FATAL:install\
.pl]\\n\";\n	exit ($EXIT_FAILURE);\n      }\n    \\
n    if     (&pg_is_installed    (\"wget\")){$pg=\\
"wget\"; $flag=\"-O\";$arg=\"--tries=2 --connect-t\
imeout=10 $wget_arg\";}\n    elsif  (&pg_is_instal\
led    (\"curl\")){$pg=\"curl\"; $flag=\"-o\";$arg\
=$curl_arg;}\n    else\n      {\n	printf stderr \"\
\\nERROR: No pg for remote file fetching [wget or \
curl][FATAL]\\n\";\n	exit ($EXIT_FAILURE);\n      \
}\n    \n    \n    if (-e $file){unlink($file);}\n\
    $exit=system \"$pg $cmd $flag$file $arg\";\n  \
  return $exit;\n  }\n\nsub pg_is_installed\n  {\n\
    my ($p, $dir)=(@_);\n    my ($r,$m, $ret);\n  \
  my ($supported, $language, $compil);\n    \n  \n\
    if ( $PG{$p})\n      {\n	$language=$PG{$p}{lan\
guage2};\n	$compil=$PG{$language}{compiler};\n    \
  }\n    \n    if ( $compil eq \"CPAN\")\n      {\\
n	if ( system (\"perl -M$p -e 1\")==$EXIT_SUCCESS)\
{$ret=1;}\n	else {$ret=0;}\n      }\n    elsif ($d\
ir)\n      {\n	if (-e \"$dir/$p\" || -e \"$dir/$p\\
\.exe\"){$ret=1;}\n	else {$ret=0;}\n      }\n    e\
lsif (-e \"$PLUGINS_DIR/$p\" || -e \"$PLUGINS_DIR/\
$p.exe\"){$ret=1;}\n    else\n      {\n	$r=`which \
$p 2>/dev/null`;\n	if ($r eq \"\"){$ret=0;}\n	else\
 {$ret=1;}\n      }\n   \n    return $ret;\n  }\ns\
ub install\n  {\n    my ($new_bin)=(@_);\n    my (\
$copied, $report);\n\n    \n    if (!$ROOT_INSTALL\
)\n      {\n	\n	if (-e \"$BIN/t_coffee\"){`$CP $BI\
N/t_coffee $INSTALL_DIR`};\n	`cp $BIN/* $PLUGINS_D\
IR`;\n	$copied=1;\n      }\n    else\n      {\n	$c\
opied=&root_run (\"You must be root to finalize th\
e installation\", \"$CP $BIN/* $INSTALL_DIR $SILEN\
T\");\n      }\n    \n     \n  if ( !$copied)\n   \
 {\n      $report=\"*!!!!!! Installation unsuccesf\
ul. The executables have been left in $BASE/bin\\n\
\";\n    }\n  elsif ( $copied && $ROOT)\n    {\n  \
    $report=\"*------ Installation succesful. Your\
 executables have been copied in $new_bin and are \
on your PATH\\n\";\n    }\n  elsif ( $copied && !$\
ROOT)\n    {\n      $report= \"*!!!!!! T-Coffee an\
d associated packages have been copied in: $new_bi\
n\\n\";\n      $report.=\"*!!!!!! This address is \
NOT in your PATH sytem variable\\n\";\n      $repo\
rt.=\"*!!!!!! You can do so by adding the followin\
g line in your ~/.bashrc file:\\n\";\n      $repor\
t.=\"*!!!!!! export PATH=$new_bin:\\$PATH\\n\";\n \
   }\n  return $report;\n}\n\nsub sign_license_ni\\
n  {\n    my $F=new FileHandle;\n    open ($F, \"l\
icense.txt\");\n    while (<$F>)\n      {\n	print \
\"$_\";\n      }\n    close ($F);\n    \n    retur\
n;\n  }\n\nsub install_pg\n  {\n    my ($pg)=(@_);\
\n    my ($report, $previous, $language, $compiler\
, $return);\n    \n    if (!$PG{$pg}{install}){ret\
urn 1;}\n    \n    $previous=&pg_is_installed ($pg\
);\n    \n    if ($PG{$pg}{update_action} eq \"no_\
update\" && $previous)\n      {\n	$PG{$pg}{old}=1;\
\n	$PG{$pg}{new}=0;\n	$return=1;\n      }\n    els\
e\n      {\n	$PG{$pg}{old}=$previous;\n	\n	if ($PG\
{$pg} {language2} eq \"Perl\"){&install_perl_packa\
ge ($pg);}\n	elsif ($BINARIES_ONLY && &install_bin\
ary_package ($pg)){$PG{$pg}{from_binary}=1;}\n	els\
if (&install_source_package ($pg)){;}\n	else \n	  \
{\n	    \n	    if (!&supported_os($OS))\n	      {\\
n		print \"!!!!!!!! $pg compilation failed, binary\
 unsupported for $OS\\n\"; \n	      }\n	    elsif \
(!($PG{$pg}{from_binary}=&install_binary_package (\
$pg)))\n	      {\n		print \"!!!!!!!! $pg compilati\
on and  binary installation failed\\n\";\n	      }\
\n	  }\n	$PG{$pg}{new}=$return=&pg_is_installed ($\
pg,$BIN);\n      }\n\n    \n    return $return;\n \
 }\nsub install_perl_package\n  {\n    my ($pg)=(@\
_);\n    my ($report, $language, $compiler);\n    \
\n    $language=$PG{$pg} {language2};\n    $compil\
er=$PG{$language}{compiler};\n    \n    if (!&pg_i\
s_installed ($pg))\n      {\n	if ( $OS eq \"window\
s\"){`perl -M$compiler -e 'install $pg'`;}\n	elsif\
 ( $ROOT eq \"sudo\"){system (\"sudo perl -M$compi\
ler -e 'install $pg'\");}\n	else {system (\"su roo\
t -c perl -M$compiler -e 'install $pg'\");}\n     \
 }\n    return &pg_is_installed ($pg);\n  }\n\n\n\\
nsub install_source_package\n  {\n    my ($pg)=(@_\
);\n    my ($report, $download, $arguments, $langu\
age, $address, $name, $ext, $main_dir, $distrib);\\
n    my $wget_tmp=\"$TMP/wget.tmp\";\n    my (@fl)\
;\n    if ( -e \"$BIN/$pg\" || -e \"$BIN/$pg.exe\"\
){return 1;}\n    \n    #\n    # check if the modu\
le exists in the repository cache \n    #\n	if( re\
po_load($pg) ) {\n		return 1;\n	}\n    \n    if ($\
pg eq \"t_coffee\")  {return   &install_t_coffee (\
$pg);}\n    elsif ($pg eq \"TMalign\"){return   &i\
nstall_TMalign ($pg);}\n    \n    chdir $DISTRIBUT\
IONS;\n    \n    $download=$PG{$pg}{source};\n    \
\n    if (($download =~/tgz/))\n      {\n	($addres\
s,$name,$ext)=($download=~/(.+\\/)([^\\/]+)(\\.tgz\
).*/);\n      }\n    elsif (($download=~/tar\\.gz/\
))\n      {\n	($address,$name,$ext)=($download=~/(\
.+\\/)([^\\/]+)(\\.tar\\.gz).*/);\n      }\n    el\
sif (($download=~/tar/))\n      {\n	($address,$nam\
e,$ext)=($download=~/(.+\\/)([^\\/]+)(\\.tar).*/);\
\n      }\n    else\n      {\n	($address,$name)=($\
download=~/(.+\\/)([^\\/]+)/);\n	$ext=\"\";\n     \
 }\n    $distrib=\"$name$ext\";\n    \n    if ( !-\
d $pg){mkdir $pg;}\n    chdir $pg;\n   \n    #get \
the distribution if available\n    if ( -e \"$DOWN\
LOAD_DIR/$distrib\")\n      {\n	`$CP $DOWNLOAD_DIR\
/$distrib .`;\n      }\n    #UNTAR and Prepare eve\
rything\n    if (!-e \"$name.tar\" && !-e \"$name\\
")\n      {\n	&check_rm ($wget_tmp);\n	print \"\\n\
------- Downloading/Installing $pg\\n\";\n	\n	if (\
!-e $distrib && &url2file (\"$download\", \"$wget_\
tmp\")==$EXIT_SUCCESS)\n	  {\n	    \n	    `mv $wge\
t_tmp $distrib`;\n	    `$CP $distrib $DOWNLOAD_DIR\
/`;\n	  }\n\n	if (!-e $distrib)\n	  {\n	    print \
\"!!!!!!! Download of $pg distribution failed\\n\"\
;\n	    print \"!!!!!!! Check Address: $PG{$pg}{so\
urce}\\n\";\n	    return 0;\n	  }\n	print \"\\n---\
---- unzipping/untaring $name\\n\";\n	if (($ext =~\
/z/))\n	  { \n	    &flush_command (\"gunzip $name$\
ext\");\n	    \n	  }\n	if (($ext =~/tar/) || ($ext\
 =~/tgz/))\n	  {\n	    &flush_command(\"tar -xvf $\
name.tar\");\n	  }\n      }\n    #Guess and enter \
the distribution directory\n    @fl=ls($p);\n    f\
oreach my $f (@fl)\n      {\n	if (-d $f)\n	  {\n	 \
   $main_dir=$f;\n	  }\n      }\n    if (-d $main_\
dir)\n	  \n      {\n	chdir $main_dir;}\n    else\n\
      {\n	print \"Error: $main_dir does not exist\\
";\n      }\n    print \"\\n------- Compiling/Inst\
alling $pg\\n\";\n    `make clean $SILENT`;\n    \\
n    \n    #\n    # SAP module\n    #\n    if ($pg\
 eq \"sap\")\n      {\n	if (-e \"./configure\")\n	\
  {\n	    #new sap distribution\n	    if ($OS eq \\
"macosx\")\n	      {\n		&replace_line_in_file (\".\
/src/galloc.h\", \"malloc.h\",  \"\");\n		&replace\
_line_in_file (\"./src/pdbprot.h\", \"malloc.h\", \
\"\");\n		&replace_line_in_file (\"./src/pdbprot.c\
\", \"malloc.h\", \"\");\n	      }\n	    \n	    &f\
lush_command (\"./configure\");\n	    &flush_comma\
nd (\"make clean\");\n	    &flush_command (\"make\\
");\n	    &check_cp (\"./src/$pg\", \"$BIN\");\n	 \
   repo_store(\"./src/$pg\");\n	  }\n	else\n	  {\n\
	    #old style distribution\n	    `rm *.o sap  sa\
p.exe ./util/aa/*.o  ./util/wt/.o $SILENT`;\n	    \
&flush_command (\"make $arguments sap\");\n	    &c\
heck_cp ($pg, \"$BIN\");\n	    repo_store($pg);\n	\
  }\n      }\n    \n    #\n    # CLUSTALW2 module\\
n    #\n    elsif ($pg eq \"clustalw2\")\n      {\\
n	&flush_command(\"./configure\");\n	&flush_comman\
d(\"make $arguments\");\n	&check_cp (\"./src/$pg\"\
, \"$BIN\");\n	repo_store(\"./src/$pg\");\n      }\
\n    \n    #\n    # FSA module\n    # \n    elsif\
 ($pg eq \"fsa\")\n      {\n	&flush_command(\"./co\
nfigure --prefix=$BIN\");\n	&flush_command(\"make \
$arguments\");\n	&flush_command (\"make install\")\
;\n\n	repo_store(\"fsa\", \"$BIN/bin\");\n	`mv $BI\
N/bin/* $BIN`;\n	`rmdir $BIN/bin`;\n      }\n    \\
n    #\n    # CLUSTALW module\n    #\n    elsif ($\
pg eq \"clustalw\")\n      {\n	&flush_command(\"ma\
ke $arguments clustalw\");\n	`$CP $pg $BIN $SILENT\
`;\n	repo_store($pg);\n      }\n    \n    #\n    #\
 MAFFT module\n    #\n    elsif ($pg eq \"mafft\")\
\n      {\n	my $base=cwd();\n	my $c;\n	\n	#compile\
 core\n	mkpath (\"./mafft/bin\");\n	mkpath (\"./ma\
fft/lib\");\n	chdir \"$base/core\";\n	`make clean \
$SILENT`;\n	&flush_command (\"make $arguments\");\\
n	&flush_command (\"make install LIBDIR=../mafft/l\
ib BINDIR=../mafft/bin\");\n	\n	#compile extension\
\n	chdir \"$base/extensions\";\n	`make clean $SILE\
NT`;\n	&flush_command (\"make $arguments\");\n	&fl\
ush_command (\"make install LIBDIR=../mafft/lib BI\
NDIR=../mafft/bin\");\n	\n	#put everything in maff\
t and copy the compiled stuff in bin\n	chdir \"$ba\
se\";\n	if ($ROOT_INSTALL)\n	  {\n	    &root_run (\
\"You Must be Root to Install MAFFT\\n\", \"mkdir \
/usr/local/mafft/;$CP mafft/lib/* /usr/local/mafft\
;$CP mafft/lib/mafft* /usr/local/bin ;$CP mafft/bi\
n/mafft /usr/local/bin/; \");\n	  }\n	else\n	  {\n\
	    `$CP mafft/lib/*  $BIN`;\n	    `$CP mafft/bin\
/mafft  $BIN`;\n	  }\n	`tar -cvf mafft.tar mafft`;\
\n	`gzip mafft.tar`;\n	`mv mafft.tar.gz $BIN`;\n	\\
n	repo_store(\"mafft/bin/mafft\", \"mafft/lib/\", \
\"$BIN/mafft.tar.gz\");\n      }\n      \n    #\n \
   # DIALIGN-TX module\n    #\n    elsif ( $pg eq \
\"dialign-tx\" )\n      {\n	my $f;\n	my $base=cwd(\
);\n\n	chdir \"./source\";\n	if ($OS eq \"macosx\"\
){&flush_command (\"cp makefile.MAC_OS makefile\")\
;}\n\n	&flush_command (\" make CPPFLAGS='-O3 -funr\
oll-loops' all\");\n	\n	chdir \"..\";\n	&check_cp \
(\"./source/$pg\", \"$BIN\");\n	repo_store(\"./sou\
rce/$pg\");\n      }\n      \n    #\n    # DIALIGN\
-T module \n    # (is the same as dialign-tx, but \
it is mantained for backward name compatibility wi\
th tcoffee)\n    #\n    elsif ( $pg eq \"dialign-t\
\" )\n      {\n	my $f;\n	my $base=cwd();\n\n	chdir\
 \"./source\";\n	if ($OS eq \"macosx\"){&flush_com\
mand (\"cp makefile.MAC_OS makefile\");}\n\n	&flus\
h_command (\" make CPPFLAGS='-O3 -funroll-loops' a\
ll\");\n	\n	chdir \"..\";\n	&check_cp (\"./source/\
dialign-tx\", \"$BIN/dialign-t\");\n	repo_store(\"\
$BIN/dialign-t\");	\n      }      \n      \n    #\\
n    # POA module\n    #\n    elsif ($pg eq \"poa\\
")\n      {\n	&flush_command (\"make $arguments po\
a\");\n	&check_cp (\"$pg\", \"$BIN\");\n	repo_stor\
e(\"$pg\");\n      }\n     \n     \n    #\n    # P\
ROBCONS module\n    #\n    elsif ( $pg eq \"probco\
ns\")\n      {\n	&add_C_libraries(\"./Probabilisti\
cModel.h\", \"list\", \"cstring\");\n	\n	`rm *.exe\
 $SILENT`;\n	&flush_command (\"make $arguments pro\
bcons\");\n	&check_cp(\"$pg\", \"$BIN/$pg\");\n	re\
po_store(\"$pg\");\n      }\n      \n    #\n    # \
PROBCONS RNA module\n    #\n    elsif ( $pg eq \"p\
robconsRNA\")\n      {\n	&add_C_libraries(\"./Prob\
abilisticModel.h\", \"list\", \"cstring\");\n	&add\
_C_libraries(\"./Main.cc\", \"iomanip\", \"cstring\
\",\"climits\");\n	`rm *.exe $SILENT`;\n	&flush_co\
mmand (\"make $arguments probcons\");\n	&check_cp(\
\"probcons\", \"$BIN/$pg\");\n	repo_store(\"$BIN/$\
pg\");\n      }\n\n	#\n	# MUSCLE module\n	#\n    e\
lsif (  $pg eq \"muscle\")\n      {	\n	`rm *.o mus\
cle muscle.exe $SILENT`;\n	if ($OS eq \"macosx\" |\
| $OS eq \"linux\")\n	  {\n	    &replace_line_in_f\
ile (\"./Makefile\", \"LDLIBS = -lm -static\",  \"\
LDLIBS = -lm\");\n	  }\n	elsif ($OS eq \"windows\"\
)\n	  {\n	    &replace_line_in_file (\"./intmath.c\
pp\",  \"double log2e\",      \"double cedric_log\\
");\n	    &replace_line_in_file (\"./intmath.cpp\"\
,  \"double log2\",       \"double log_notuse\");\\
n	    &replace_line_in_file (\"./intmath.cpp\",  \\
"double cedric_log\", \"double log2e\");\n	  }\n	&\
flush_command (\"make $arguments all\");\n	&check_\
cp(\"$pg\", \"$BIN\");\n	repo_store(\"$pg\");	\n  \
    }\n      \n     #\n     # MUS4 module\n     #\\
n     elsif (  $pg eq \"mus4\")\n      {\n	`rm *.o\
 muscle muscle.exe $SILENT`;\n	&flush_command (\".\
/mk\");\n	&check_cp(\"$pg\", \"$BIN\");\n	repo_sto\
re(\"$pg\");	\n      }\n      \n    #\n    # PCMA \
module\n    #\n    elsif ( $pg eq \"pcma\")\n     \
 {\n	if ($OS eq \"macosx\")\n	  {\n	    &replace_l\
ine_in_file (\"./alcomp2.c\", \"malloc.h\",  \"\")\
;\n	  }\n	&flush_command (\"make $arguments pcma\"\
);\n	&check_cp(\"$pg\", \"$BIN\");\n	repo_store(\"\
$pg\");	\n      }\n      \n    #\n    # KALIGN mod\
ule\n    #\n    elsif ($pg eq \"kalign\")\n      {\
\n	&flush_command (\"./configure\");\n	&flush_comm\
and(\"make $arguments\");\n	&check_cp (\"$pg\",$BI\
N);\n	repo_store(\"$pg\");	\n      }\n      \n    \
#\n    # AMAP module\n    #\n    elsif ( $pg eq \"\
amap\")\n      {\n	&add_C_libraries(\"./Amap.cc\",\
 \"iomanip\", \"cstring\",\"climits\");	\n	`make c\
lean $SILENT`;\n	&flush_command (\"make $arguments\
 all\");\n	&check_cp (\"$pg\", $BIN);\n	repo_store\
(\"$pg\");	\n      }\n      \n    #\n    # PRODA m\
odule\n    #\n    elsif ( $pg eq \"proda\")\n     \
 {\n	&add_C_libraries(\"AlignedFragment.h\", \"vec\
tor\", \"iostream\", \"cstring\",\"cstdlib\");\n	&\
add_C_libraries(\"Main.cc\", \"vector\", \"climits\
\");	\n	&add_C_libraries(\"Sequence.cc\", \"stdlib\
.h\", \"cstdio\");	\n	&flush_command (\"make $argu\
ments all\");\n	&check_cp (\"$pg\", $BIN);\n	repo_\
store(\"$pg\");	\n      }\n      \n    #\n    # PR\
ANK module\n    #\n    elsif ( $pg eq \"prank\")\n\
      {\n	&flush_command (\"make $arguments all\")\
;\n	&check_cp (\"$pg\", $BIN);\n	repo_store(\"$pg\\
");	\n      }\n      \n    #\n    # !!!! MUSTANG m\
odule\n    #\n     elsif ( $pg eq \"mustang\")\n  \
    {\n	&flush_command (\"rm ./bin/*\");\n	&flush_\
command (\"make $arguments all\");\n\n	if ( $OS=~/\
windows/){&flush_command(\"cp ./bin/* $BIN/mustang\
.exe\");}\n	else {&flush_command(\"cp ./bin/* $BIN\
/mustang\");}\n	\n	repo_store(\"$BIN/mustang\");\n\
      }\n\n	#\n	# RNAplfold module\n	#\n    elsif \
( $pg eq \"RNAplfold\")\n      {\n	&flush_command(\
\"./configure\");\n	&flush_command (\"make $argume\
nts all\");\n	&check_cp(\"./Progs/RNAplfold\", \"$\
BIN\");\n	&check_cp(\"./Progs/RNAalifold\", \"$BIN\
\");\n	&check_cp(\"./Progs/RNAfold\", \"$BIN\");\n\
	\n	repo_store(\"./Progs/RNAplfold\", \"./Progs/RN\
Aalifold\", \"./Progs/RNAfold\");\n      }\n      \
\n    #\n    # !!! RETREE module\n    #\n    elsif\
 ( $pg eq \"retree\")\n      {\n	chdir \"src\";\n	\
&flush_command (\"make $arguments all\");\n	&flush\
_command (\"make put\");\n	system \"cp ../exe/* $B\
IN\";\n	\n	repo_store(\"retree\", \"../exe\");\n  \
    }\n	\n    chdir $CDIR;\n    return &pg_is_inst\
alled ($pg, $BIN);\n  }\n\nsub install_t_coffee\n \
 {\n    my ($pg)=(@_);\n    my ($report,$cflags, $\
arguments, $language, $compiler) ;\n    #1-Install\
 T-Coffee\n    chdir \"t_coffee_source\";\n    &fl\
ush_command (\"make clean\");\n    print \"\\n----\
--- Compiling T-Coffee\\n\";\n    $language=$PG{$p\
g} {language2};\n    $arguments=$PG{$language}{arg\
uments};\n    if (!($arguments =~/CFLAGS/)){$argum\
ents .= \" CFLAGS=-O2 \";}\n\n    if ( $CC ne \"\"\
){&flush_command (\"make -i $arguments t_coffee\")\
;}\n    &check_cp ($pg, $BIN);\n    \n    chdir $C\
DIR;\n    return &pg_is_installed ($pg, $BIN);\n  \
}\nsub install_TMalign\n  {\n    my ($pg)=(@_);\n \
   my $report;\n    chdir \"t_coffee_source\";\n  \
  print \"\\n------- Compiling TMalign\\n\";\n    \
`rm TMalign TMalign.exe $SILENT`;\n    if ( $FC ne\
 \"\"){&flush_command (\"make -i $PG{Fortran}{argu\
ments} TMalign\");}\n    &check_cp ($pg, $BIN);\n \
   repo_store($pg);\n\n    if ( !-e \"$BIN/$pg\" &\
& pg_has_binary_distrib ($pg))\n      {\n	print \"\
!!!!!!! Compilation of $pg impossible. Will try to\
 install from binary\\n\";\n	return &install_binar\
y_package ($pg);\n      }\n    chdir $CDIR;\n    r\
eturn &pg_is_installed ($pg, $BIN);\n  }\n\nsub pg\
_has_binary_distrib\n  {\n    my ($pg)=(@_);\n    \
if ($PG{$pg}{windows}){return 1;}\n    elsif ($PG{\
$pg}{osx}){return 1;}\n    elsif ($PG{$pg}{linux})\
{return 1;}\n    return 0;\n  }\nsub install_binar\
y_package\n  {\n    my ($pg)=(@_);\n    my ($base,\
$report,$name, $download, $arguments, $language, $\
dir);\n    my $isdir;\n    &input_os();\n    \n   \
 if (!&supported_os($OS)){return 0;}\n    if ( $PG\
{$pg}{binary}){$name=$PG{$pg}{binary};}\n    else \
\n      {\n	$name=$pg;\n	if ( $OS eq \"windows\"){\
$name.=\".exe\";}\n      }\n    \n    $download=\"\
$WEB_BASE/Packages/Binaries/$OS/$name\";\n    \n  \
  $base=cwd();\n    chdir $TMP;\n    \n    if (!-e\
 $name)\n      {\n	`rm x $SILENT`;\n	if ( url2file\
(\"$download\",\"x\")==$EXIT_SUCCESS)\n	  {\n	    \
`mv x $name`;\n	  }\n      }\n    \n    if (!-e $n\
ame)\n      {\n	print \"!!!!!!! $PG{$pg}{dname}: D\
ownload of $pg binary failed\\n\";\n	print \"!!!!!\
!! $PG{$pg}{dname}: Check Address: $download\\n\";\
\n	return 0;\n      }\n    print \"\\n------- Inst\
alling $pg\\n\";\n    \n    if ($name =~/tar\\.gz/\
)\n      {\n	`gunzip  $name`;\n	`tar -xvf $pg.tar`\
;\n	chdir $pg;\n	if ( $pg eq \"mafft\")\n	  {\n	  \
  if ($ROOT_INSTALL)\n	      {\n		&root_run (\"You\
 Must be Roor to Install MAFFT\\n\", \"$CP mafft/b\
in/* /usr/local/mafft;mkdir /usr/local/mafft/; $CP\
 mafft/lib/* /usr/local/bin/\");\n	      }\n	    e\
lse\n	      {\n		`$CP $TMP/$pg/bin/* $BIN $SILENT`\
;\n		`$CP $TMP/$pg/lib/* $BIN $SILENT`;\n	      }\\
n	  }\n	else\n	  {\n	    if (-e \"$TMP/$pg/data\")\
{`$CP $TMP/$pg/data/* $TCM $SILENT`;}\n	    if (!(\
$pg=~/\\*/)){`rm -rf $pg`;}\n	  }\n      }\n    el\
se\n      {\n	&check_cp (\"$pg\", \"$BIN\");\n	`ch\
mod u+x $BIN/$pg`; \n	unlink ($pg);\n      }\n    \
chdir $base;\n    $PG{$pg}{from_binary}=1;\n    re\
turn &pg_is_installed ($pg, $BIN);\n  }\n\nsub add\
_dir \n  {\n    my $dir=@_[0];\n    \n    if (!-e \
$dir && !-d $dir)\n      {\n	my @l;\n	umask (0000)\
;\n	@l=mkpath ($dir,{mode => 0777});\n	\n      }\n\
    else\n      {\n	return 0;\n      }\n  }\nsub c\
heck_rm \n  {\n    my ($file)=(@_);\n    \n    if \
( -e $file)\n      {\n	return unlink($file);\n    \
  }\n    return 0;\n  }\nsub check_cp\n  {\n    my\
 ($from, $to)=(@_);\n    if ( !-e $from && -e \"$f\
rom\\.exe\"){$from=\"$from\\.exe\";}\n    if ( !-e\
 $from){return 0;}\n        \n    `$CP $from $to`;\
\n    return 1;\n  }\n\nsub repo_store \n{\n   # c\
heck that all required data are available\n   if( \
$REPO_ROOT eq \"\" ) { return; }\n\n\n    # extrac\
t the package name from the specified path\n    my\
 $pg =`basename $_[0]`;\n    chomp($pg);\n	\n    m\
y $VER = $PG{$pg}{version};\n    my $CACHE = \"$RE\
PO_ROOT/$pg/$VER/$OSNAME-$OSARCH\"; \n    \n    pr\
int \"-------- Storing package: \\\"$pg\\\" to pat\
h: $CACHE\\n\";\n    \n    # clean the cache path \
if exists and create it again\n    `rm -rf $CACHE`\
;\n    `mkdir -p $CACHE`;\n    \n 	for my $path (@\
_) {\n\n	    # check if it is a single file \n	 	i\
f( -f $path ) {\n	    	`cp $path $CACHE`;\n		}\n		\
# .. or a directory, in this case copy all the con\
tent \n		elsif( -d $path ) {\n			opendir(IMD, $pat\
h);\n			my @thefiles= readdir(IMD);\n			closedir(I\
MD);\n			\n			for my $_file (@thefiles) {\n				if(\
 $_file ne \".\" && $_file ne \"..\") {\n	    			`\
cp $path/$_file $CACHE`;\n				}\n			}\n		} \n	}	  \
 \n    \n	\n}   \n\nsub repo_load \n{\n    my ($pg\
)=(@_);\n\n    # check that all required data are \
available\n    if( $REPO_ROOT eq \"\" ) { return 0\
; }\n\n    my $VER = $PG{$pg}{version};\n    my $C\
ACHE = \"$REPO_ROOT/$pg/$VER/$OSNAME-$OSARCH\"; \n\
    if( !-e \"$CACHE/$pg\" ) {\n   	 	print \"----\
---- Module \\\"$pg\\\" NOT found on repository ca\
che.\\n\";\n    	return 0;\n    }\n    \n    print\
 \"-------- Module \\\"$pg\\\" found on repository\
 cache. Using copy on path: $CACHE\\n\";\n    `cp \
$CACHE/* $BIN`;\n    return 1;\n}\n\nsub check_fil\
e_list_exists \n  {\n    my ($base, @flist)=(@_);\\
n    my $f;\n\n    foreach $f (@flist)\n      {\n	\
if ( !-e \"$base/$f\"){return 0;}\n      }\n    re\
turn 1;\n  }\nsub ls\n  {\n    my $f=@_[0];\n    m\
y @fl;\n    chomp(@fl=`ls -1 $f`);\n    return @fl\
;\n  }\nsub flush_command\n  {\n    my $command=@_\
[0];\n    my $F=new FileHandle;\n    open ($F, \"$\
command|\");\n    while (<$F>){print \"    --- $_\\
";}\n    close ($F);\n  }    \n\nsub input_install\
ation_directory\n  {\n    my $dir=@_[0];\n    my $\
new;\n    \n    print \"------- The current instal\
lation directory is: [$dir]\\n\";\n    print \"???\
???? Return to keep the default or new value:\";\n\
   \n    if ($NO_QUESTION==0)\n      {\n	chomp ($n\
ew=<stdin>);\n	while ( $new ne \"\" && !input_yes \
(\"You have entered $new. Is this correct? ([y]/n)\
:\"))\n	  {\n	    print \"???????New installation \
directory:\";\n	    chomp ($new=<stdin>);\n	  }\n	\
$dir=($new eq \"\")?$dir:$new;\n	$dir=~s/\\/$//;\n\
      }\n    \n    if ( -d $dir){return $dir;}\n  \
  elsif (&root_run (\"You must be root to create $\
dir\",\"mkdir $dir\")==$EXIT_SUCCESS){return $dir;\
}\n    else\n      {\n	print \"!!!!!!! $dir could \
not be created\\n\";\n	if ( $NO_QUESTION)\n	  {\n	\
    return \"\";\n	  }\n	elsif ( &input_yes (\"???\
???? Do you want to provide a new directory([y]/n)\
?:\"))\n	  {\n	    return input_installation_direc\
tory ($dir);\n	  }\n	else\n	  {\n	    return \"\";\
\n	  }\n      }\n    \n  }\nsub input_yes\n  {\n  \
  my $question =@_[0];\n    my $answer;\n\n    if \
($NO_QUESTION==1){return 1;}\n    \n    if ($quest\
ion eq \"\"){$question=\"??????? Do you wish to pr\
oceed ([y]/n)?:\";}\n    print $question;\n    cho\
mp($answer=lc(<STDIN>));\n    if (($answer=~/^y/) \
|| $answer eq \"\"){return 1;}\n    elsif ( ($answ\
er=~/^n/)){return 0;}\n    else\n      {\n	return \
input_yes($question);\n      }\n  }\nsub root_run\\
n  {\n    my ($txt, $cmd)=(@_);\n    \n    if ( sy\
stem ($cmd)==$EXIT_SUCCESS){return $EXIT_SUCCESS;}\
\n    else \n      {\n	print \"------- $txt\\n\";\\
n	if ( $ROOT eq \"sudo\"){return system (\"sudo $c\
md\");}\n	else {return system (\"su root -c \\\"$c\
md\\\"\");}\n      }\n  }\nsub get_root\n  {\n    \
if (&pg_is_installed (\"sudo\")){return \"sudo\";}\
\n    else {return \"su\";}\n  }\n\nsub get_os\n  \
{\n    my $raw_os=`uname`;\n    my $os;\n\n    $ra\
w_os=lc ($raw_os);\n    \n    if ($raw_os =~/cygwi\
n/){$os=\"windows\";}\n    elsif ($raw_os =~/linux\
/){$os=\"linux\";}\n    elsif ($raw_os =~/osx/){$o\
s=\"macosx\";}\n    elsif ($raw_os =~/darwin/){$os\
=\"macosx\";}\n    else\n      {\n	$os=$raw_os;\n \
     }\n    return $os;\n  }\nsub input_os\n  {\n \
   my $answer;\n    if ($OS) {return $OS;}\n    \n\
    print \"??????? which os do you use: [w]indows\
, [l]inux, [m]acosx:?\";\n    $answer=lc(<STDIN>);\
\n\n    if (($answer=~/^m/)){$OS=\"macosx\";}\n   \
 elsif ( ($answer=~/^w/)){$OS=\"windows\";}\n    e\
lsif ( ($answer=~/^linux/)){$OS=\"linux\";}\n    \\
n    else\n      {\n	return &input_os();\n      }\\
n    return $OS;\n  }\n\nsub supported_os\n  {\n  \
  my ($os)=(@_[0]);\n    return $SUPPORTED_OS{$os}\
;\n  }\n    \n    \n\n\nsub update_tclinkdb \n  {\\
n    my $file =@_[0];\n    my $name;\n    my $F=ne\
w FileHandle;\n    my ($download, $address, $name,\
 $l, $db);\n    \n    if ( $file eq \"update\"){$f\
ile=$TCLINKDB_ADDRESS;}\n    \n    if ( $file =~/h\
ttp:\\/\\// || $file =~/ftp:\\/\\//)\n      {\n	($\
address, $name)=($download=~/(.*)\\/([^\\/]+)$/);\\
n	`rm x $SILENT`;\n	if (&url2file ($file,\"x\")==$\
EXIT_SUCCESS)\n	  {\n	    print \"------- Susscess\
ful upload of $name\";\n	    `mv x $name`;\n	    $\
file=$name;\n	  }\n      }\n    open ($F, \"$file\\
");\n    while (<$F>)\n      {\n	my $l=$_;\n	if ((\
$l =~/^\\/\\//) || ($db=~/^#/)){;}\n	elsif ( !($l \
=~/\\w/)){;}\n	else\n	  {\n	    my @v=split (/\\s+\
/, $l);\n	    if ( $l=~/^MODE/)\n	      {\n		$MODE\
{$v[1]}{$v[2]}=$v[3];\n	      }\n	    elsif ($l=~/\
^PG/)\n	      {\n		$PG{$v[1]}{$v[2]}=$v[3];\n	    \
  }\n	  }\n      }\n    close ($F);\n    &post_pro\
cess_PG();\n    return;\n  }\n\n\n\nsub initialize\
_PG\n  {\n\n$PG{\"t_coffee\"}{\"4_TCOFFEE\"}=\"TCO\
FFEE\";\n$PG{\"t_coffee\"}{\"type\"}=\"sequence_mu\
ltiple_aligner\";\n$PG{\"t_coffee\"}{\"ADDRESS\"}=\
\"http://www.tcoffee.org\";\n$PG{\"t_coffee\"}{\"l\
anguage\"}=\"C\";\n$PG{\"t_coffee\"}{\"language2\"\
}=\"C\";\n$PG{\"t_coffee\"}{\"source\"}=\"http://w\
ww.tcoffee.org/Packages/T-COFFEE_distribution.tar.\
gz\";\n$PG{\"t_coffee\"}{\"update_action\"}=\"alwa\
ys\";\n$PG{\"t_coffee\"}{\"mode\"}=\"tcoffee,mcoff\
ee,rcoffee,expresso,3dcoffee\";\n$PG{\"clustalw2\"\
}{\"4_TCOFFEE\"}=\"CLUSTALW2\";\n$PG{\"clustalw2\"\
}{\"type\"}=\"sequence_multiple_aligner\";\n$PG{\"\
clustalw2\"}{\"ADDRESS\"}=\"http://www.clustal.org\
\";\n$PG{\"clustalw2\"}{\"language\"}=\"C++\";\n$P\
G{\"clustalw2\"}{\"language2\"}=\"CXX\";\n$PG{\"cl\
ustalw2\"}{\"source\"}=\"http://www.clustal.org/do\
wnload/2.0.10/clustalw-2.0.10-src.tar.gz\";\n$PG{\\
"clustalw2\"}{\"mode\"}=\"mcoffee,rcoffee\";\n$PG{\
\"clustalw2\"}{\"version\"}=\"2.0.10\";\n$PG{\"clu\
stalw\"}{\"4_TCOFFEE\"}=\"CLUSTALW\";\n$PG{\"clust\
alw\"}{\"type\"}=\"sequence_multiple_aligner\";\n$\
PG{\"clustalw\"}{\"ADDRESS\"}=\"http://www.clustal\
.org\";\n$PG{\"clustalw\"}{\"language\"}=\"C\";\n$\
PG{\"clustalw\"}{\"language2\"}=\"C\";\n$PG{\"clus\
talw\"}{\"source\"}=\"http://www.clustal.org/downl\
oad/1.X/ftp-igbmc.u-strasbg.fr/pub/ClustalW/clusta\
lw1.82.UNIX.tar.gz\";\n$PG{\"clustalw\"}{\"mode\"}\
=\"mcoffee,rcoffee\";\n$PG{\"clustalw\"}{\"version\
\"}=\"1.82\";\n$PG{\"dialign-t\"}{\"4_TCOFFEE\"}=\\
"DIALIGNT\";\n$PG{\"dialign-t\"}{\"type\"}=\"seque\
nce_multiple_aligner\";\n$PG{\"dialign-t\"}{\"ADDR\
ESS\"}=\"http://dialign-tx.gobics.de/\";\n$PG{\"di\
align-t\"}{\"DIR\"}=\"/usr/share/dialign-tx/\";\n$\
PG{\"dialign-t\"}{\"language\"}=\"C\";\n$PG{\"dial\
ign-t\"}{\"language2\"}=\"C\";\n$PG{\"dialign-t\"}\
{\"source\"}=\"http://dialign-tx.gobics.de/DIALIGN\
-TX_1.0.2.tar.gz\";\n$PG{\"dialign-t\"}{\"mode\"}=\
\"mcoffee\";\n$PG{\"dialign-t\"}{\"binary\"}=\"dia\
lign-t\";\n$PG{\"dialign-t\"}{\"version\"}=\"1.0.2\
\";\n$PG{\"dialign-tx\"}{\"4_TCOFFEE\"}=\"DIALIGNT\
X\";\n$PG{\"dialign-tx\"}{\"type\"}=\"sequence_mul\
tiple_aligner\";\n$PG{\"dialign-tx\"}{\"ADDRESS\"}\
=\"http://dialign-tx.gobics.de/\";\n$PG{\"dialign-\
tx\"}{\"DIR\"}=\"/usr/share/dialign-tx/\";\n$PG{\"\
dialign-tx\"}{\"language\"}=\"C\";\n$PG{\"dialign-\
tx\"}{\"language2\"}=\"C\";\n$PG{\"dialign-tx\"}{\\
"source\"}=\"http://dialign-tx.gobics.de/DIALIGN-T\
X_1.0.2.tar.gz\";\n$PG{\"dialign-tx\"}{\"mode\"}=\\
"mcoffee\";\n$PG{\"dialign-tx\"}{\"binary\"}=\"dia\
lign-tx\";\n$PG{\"dialign-tx\"}{\"version\"}=\"1.0\
.2\";\n$PG{\"poa\"}{\"4_TCOFFEE\"}=\"POA\";\n$PG{\\
"poa\"}{\"type\"}=\"sequence_multiple_aligner\";\n\
$PG{\"poa\"}{\"ADDRESS\"}=\"http://www.bioinformat\
ics.ucla.edu/poa/\";\n$PG{\"poa\"}{\"language\"}=\\
"C\";\n$PG{\"poa\"}{\"language2\"}=\"C\";\n$PG{\"p\
oa\"}{\"source\"}=\"http://downloads.sourceforge.n\
et/poamsa/poaV2.tar.gz\";\n$PG{\"poa\"}{\"DIR\"}=\\
"/usr/share/\";\n$PG{\"poa\"}{\"FILE1\"}=\"blosum8\
0.mat\";\n$PG{\"poa\"}{\"mode\"}=\"mcoffee\";\n$PG\
{\"poa\"}{\"binary\"}=\"poa\";\n$PG{\"poa\"}{\"ver\
sion\"}=\"2.0\";\n$PG{\"probcons\"}{\"4_TCOFFEE\"}\
=\"PROBCONS\";\n$PG{\"probcons\"}{\"type\"}=\"sequ\
ence_multiple_aligner\";\n$PG{\"probcons\"}{\"ADDR\
ESS\"}=\"http://probcons.stanford.edu/\";\n$PG{\"p\
robcons\"}{\"language2\"}=\"CXX\";\n$PG{\"probcons\
\"}{\"language\"}=\"C++\";\n$PG{\"probcons\"}{\"so\
urce\"}=\"http://probcons.stanford.edu/probcons_v1\
_12.tar.gz\";\n$PG{\"probcons\"}{\"mode\"}=\"mcoff\
ee\";\n$PG{\"probcons\"}{\"binary\"}=\"probcons\";\
\n$PG{\"probcons\"}{\"version\"}=\"1.12\";\n$PG{\"\
mafft\"}{\"4_TCOFFEE\"}=\"MAFFT\";\n$PG{\"mafft\"}\
{\"type\"}=\"sequence_multiple_aligner\";\n$PG{\"m\
afft\"}{\"ADDRESS\"}=\"http://align.bmr.kyushu-u.a\
c.jp/mafft/online/server/\";\n$PG{\"mafft\"}{\"lan\
guage\"}=\"C\";\n$PG{\"mafft\"}{\"language\"}=\"C\\
";\n$PG{\"mafft\"}{\"source\"}=\"http://align.bmr.\
kyushu-u.ac.jp/mafft/software/mafft-6.603-with-ext\
ensions-src.tgz\";\n$PG{\"mafft\"}{\"windows\"}=\"\
http://align.bmr.kyushu-u.ac.jp/mafft/software/maf\
ft-6.603-mingw.tar\";\n$PG{\"mafft\"}{\"mode\"}=\"\
mcoffee,rcoffee\";\n$PG{\"mafft\"}{\"binary\"}=\"m\
afft.tar.gz\";\n$PG{\"mafft\"}{\"version\"}=\"6.60\
3\";\n$PG{\"muscle\"}{\"4_TCOFFEE\"}=\"MUSCLE\";\n\
$PG{\"muscle\"}{\"type\"}=\"sequence_multiple_alig\
ner\";\n$PG{\"muscle\"}{\"ADDRESS\"}=\"http://www.\
drive5.com/muscle/\";\n$PG{\"muscle\"}{\"language\\
"}=\"C++\";\n$PG{\"muscle\"}{\"language2\"}=\"GPP\\
";\n$PG{\"muscle\"}{\"source\"}=\"http://www.drive\
5.com/muscle/downloads3.7/muscle3.7_src.tar.gz\";\\
n$PG{\"muscle\"}{\"windows\"}=\"http://www.drive5.\
com/muscle/downloads3.7/muscle3.7_win32.zip\";\n$P\
G{\"muscle\"}{\"linux\"}=\"http://www.drive5.com/m\
uscle/downloads3.7/muscle3.7_linux_ia32.tar.gz\";\\
n$PG{\"muscle\"}{\"mode\"}=\"mcoffee,rcoffee\";\n$\
PG{\"muscle\"}{\"version\"}=\"3.7\";\n$PG{\"mus4\"\
}{\"4_TCOFFEE\"}=\"MUS4\";\n$PG{\"mus4\"}{\"type\"\
}=\"sequence_multiple_aligner\";\n$PG{\"mus4\"}{\"\
ADDRESS\"}=\"http://www.drive5.com/muscle/\";\n$PG\
{\"mus4\"}{\"language\"}=\"C++\";\n$PG{\"mus4\"}{\\
"language2\"}=\"GPP\";\n$PG{\"mus4\"}{\"source\"}=\
\"http://www.drive5.com/muscle/muscle4.0_src.tar.g\
z\";\n$PG{\"mus4\"}{\"mode\"}=\"mcoffee,rcoffee\";\
\n$PG{\"mus4\"}{\"version\"}=\"4.0\";\n$PG{\"pcma\\
"}{\"4_TCOFFEE\"}=\"PCMA\";\n$PG{\"pcma\"}{\"type\\
"}=\"sequence_multiple_aligner\";\n$PG{\"pcma\"}{\\
"ADDRESS\"}=\"ftp://iole.swmed.edu/pub/PCMA/\";\n$\
PG{\"pcma\"}{\"language\"}=\"C\";\n$PG{\"pcma\"}{\\
"language2\"}=\"C\";\n$PG{\"pcma\"}{\"source\"}=\"\
ftp://iole.swmed.edu/pub/PCMA/pcma.tar.gz\";\n$PG{\
\"pcma\"}{\"mode\"}=\"mcoffee\";\n$PG{\"pcma\"}{\"\
version\"}=\"1.0\";\n$PG{\"kalign\"}{\"4_TCOFFEE\"\
}=\"KALIGN\";\n$PG{\"kalign\"}{\"type\"}=\"sequenc\
e_multiple_aligner\";\n$PG{\"kalign\"}{\"ADDRESS\"\
}=\"http://msa.cgb.ki.se\";\n$PG{\"kalign\"}{\"lan\
guage\"}=\"C\";\n$PG{\"kalign\"}{\"language2\"}=\"\
C\";\n$PG{\"kalign\"}{\"source\"}=\"http://msa.cgb\
.ki.se/downloads/kalign/current.tar.gz\";\n$PG{\"k\
align\"}{\"mode\"}=\"mcoffee\";\n$PG{\"kalign\"}{\\
"version\"}=\"1.0\";\n$PG{\"amap\"}{\"4_TCOFFEE\"}\
=\"AMAP\";\n$PG{\"amap\"}{\"type\"}=\"sequence_mul\
tiple_aligner\";\n$PG{\"amap\"}{\"ADDRESS\"}=\"htt\
p://bio.math.berkeley.edu/amap/\";\n$PG{\"amap\"}{\
\"language\"}=\"C++\";\n$PG{\"amap\"}{\"language2\\
"}=\"CXX\";\n$PG{\"amap\"}{\"source\"}=\"http://am\
ap-align.googlecode.com/files/amap.2.0.tar.gz\";\n\
$PG{\"amap\"}{\"mode\"}=\"mcoffee\";\n$PG{\"amap\"\
}{\"version\"}=\"2.0\";\n$PG{\"proda\"}{\"4_TCOFFE\
E\"}=\"PRODA\";\n$PG{\"proda\"}{\"type\"}=\"sequen\
ce_multiple_aligner\";\n$PG{\"proda\"}{\"ADDRESS\"\
}=\"http://proda.stanford.edu\";\n$PG{\"proda\"}{\\
"language\"}=\"C++\";\n$PG{\"proda\"}{\"language2\\
"}=\"CXX\";\n$PG{\"proda\"}{\"source\"}=\"http://p\
roda.stanford.edu/proda_1_0.tar.gz\";\n$PG{\"proda\
\"}{\"mode\"}=\"mcoffee\";\n$PG{\"proda\"}{\"versi\
on\"}=\"1.0\";\n$PG{\"fsa\"}{\"4_TCOFFEE\"}=\"FSA\\
";\n$PG{\"fsa\"}{\"type\"}=\"sequence_multiple_ali\
gner\";\n$PG{\"fsa\"}{\"ADDRESS\"}=\"http://fsa.so\
urceforge.net/\";\n$PG{\"fsa\"}{\"language\"}=\"C+\
+\";\n$PG{\"fsa\"}{\"language2\"}=\"CXX\";\n$PG{\"\
fsa\"}{\"source\"}=\"http://sourceforge.net/projec\
ts/fsa/files/fsa-1.15.3.tar.gz/download/\";\n$PG{\\
"fsa\"}{\"mode\"}=\"mcoffee\";\n$PG{\"fsa\"}{\"ver\
sion\"}=\"1.15.3\";\n$PG{\"prank\"}{\"4_TCOFFEE\"}\
=\"PRANK\";\n$PG{\"prank\"}{\"type\"}=\"sequence_m\
ultiple_aligner\";\n$PG{\"prank\"}{\"ADDRESS\"}=\"\
http://www.ebi.ac.uk/goldman-srv/prank/\";\n$PG{\"\
prank\"}{\"language\"}=\"C++\";\n$PG{\"prank\"}{\"\
language2\"}=\"CXX\";\n$PG{\"prank\"}{\"source\"}=\
\"http://www.ebi.ac.uk/goldman-srv/prank/src/prank\
/prank.src.100303.tgz\";\n$PG{\"prank\"}{\"mode\"}\
=\"mcoffee\";\n$PG{\"prank\"}{\"version\"}=\"10030\
3\";\n$PG{\"sap\"}{\"4_TCOFFEE\"}=\"SAP\";\n$PG{\"\
sap\"}{\"type\"}=\"structure_pairwise_aligner\";\n\
$PG{\"sap\"}{\"ADDRESS\"}=\"http://mathbio.nimr.mr\
c.ac.uk/wiki/Software\";\n$PG{\"sap\"}{\"language\\
"}=\"C\";\n$PG{\"sap\"}{\"language2\"}=\"C\";\n$PG\
{\"sap\"}{\"source\"}=\"http://mathbio.nimr.mrc.ac\
.uk/download/sap-1.1.1.tar.gz\";\n$PG{\"sap\"}{\"m\
ode\"}=\"expresso,3dcoffee\";\n$PG{\"sap\"}{\"vers\
ion\"}=\"1.1.1\";\n$PG{\"TMalign\"}{\"4_TCOFFEE\"}\
=\"TMALIGN\";\n$PG{\"TMalign\"}{\"type\"}=\"struct\
ure_pairwise_aligner\";\n$PG{\"TMalign\"}{\"ADDRES\
S\"}=\"http://zhang.bioinformatics.ku.edu/TM-align\
/TMalign.f\";\n$PG{\"TMalign\"}{\"language\"}=\"Fo\
rtran\";\n$PG{\"TMalign\"}{\"language2\"}=\"Fortra\
n\";\n$PG{\"TMalign\"}{\"source\"}=\"http://zhang.\
bioinformatics.ku.edu/TM-align/TMalign.f\";\n$PG{\\
"TMalign\"}{\"linux\"}=\"http://zhang.bioinformati\
cs.ku.edu/TM-align/TMalign_32.gz\";\n$PG{\"TMalign\
\"}{\"mode\"}=\"expresso,3dcoffee\";\n$PG{\"TMalig\
n\"}{\"version\"}=\"1.0\";\n$PG{\"mustang\"}{\"4_T\
COFFEE\"}=\"MUSTANG\";\n$PG{\"mustang\"}{\"type\"}\
=\"structure_pairwise_aligner\";\n$PG{\"mustang\"}\
{\"ADDRESS\"}=\"http://www.cs.mu.oz.au/~arun/musta\
ng\";\n$PG{\"mustang\"}{\"language\"}=\"C++\";\n$P\
G{\"mustang\"}{\"language2\"}=\"CXX\";\n$PG{\"must\
ang\"}{\"source\"}=\"http://ww2.cs.mu.oz.au/~arun/\
mustang/mustang_v3.2.1.tgz\";\n$PG{\"mustang\"}{\"\
mode\"}=\"expresso,3dcoffee\";\n$PG{\"mustang\"}{\\
"version\"}=\"3.2.1\";\n$PG{\"lsqman\"}{\"4_TCOFFE\
E\"}=\"LSQMAN\";\n$PG{\"lsqman\"}{\"type\"}=\"stru\
cture_pairwise_aligner\";\n$PG{\"lsqman\"}{\"ADDRE\
SS\"}=\"empty\";\n$PG{\"lsqman\"}{\"language\"}=\"\
empty\";\n$PG{\"lsqman\"}{\"language2\"}=\"empty\"\
;\n$PG{\"lsqman\"}{\"source\"}=\"empty\";\n$PG{\"l\
sqman\"}{\"update_action\"}=\"never\";\n$PG{\"lsqm\
an\"}{\"mode\"}=\"expresso,3dcoffee\";\n$PG{\"alig\
n_pdb\"}{\"4_TCOFFEE\"}=\"ALIGN_PDB\";\n$PG{\"alig\
n_pdb\"}{\"type\"}=\"structure_pairwise_aligner\";\
\n$PG{\"align_pdb\"}{\"ADDRESS\"}=\"empty\";\n$PG{\
\"align_pdb\"}{\"language\"}=\"empty\";\n$PG{\"ali\
gn_pdb\"}{\"language2\"}=\"empty\";\n$PG{\"align_p\
db\"}{\"source\"}=\"empty\";\n$PG{\"align_pdb\"}{\\
"update_action\"}=\"never\";\n$PG{\"align_pdb\"}{\\
"mode\"}=\"expresso,3dcoffee\";\n$PG{\"fugueali\"}\
{\"4_TCOFFEE\"}=\"FUGUE\";\n$PG{\"fugueali\"}{\"ty\
pe\"}=\"structure_pairwise_aligner\";\n$PG{\"fugue\
ali\"}{\"ADDRESS\"}=\"http://www-cryst.bioc.cam.ac\
.uk/fugue/download.html\";\n$PG{\"fugueali\"}{\"la\
nguage\"}=\"empty\";\n$PG{\"fugueali\"}{\"language\
2\"}=\"empty\";\n$PG{\"fugueali\"}{\"source\"}=\"e\
mpty\";\n$PG{\"fugueali\"}{\"update_action\"}=\"ne\
ver\";\n$PG{\"fugueali\"}{\"mode\"}=\"expresso,3dc\
offee\";\n$PG{\"dalilite.pl\"}{\"4_TCOFFEE\"}=\"DA\
LILITEc\";\n$PG{\"dalilite.pl\"}{\"type\"}=\"struc\
ture_pairwise_aligner\";\n$PG{\"dalilite.pl\"}{\"A\
DDRESS\"}=\"built_in\";\n$PG{\"dalilite.pl\"}{\"AD\
DRESS2\"}=\"http://www.ebi.ac.uk/Tools/webservices\
/services/dalilite\";\n$PG{\"dalilite.pl\"}{\"lang\
uage\"}=\"Perl\";\n$PG{\"dalilite.pl\"}{\"language\
2\"}=\"Perl\";\n$PG{\"dalilite.pl\"}{\"source\"}=\\
"empty\";\n$PG{\"dalilite.pl\"}{\"update_action\"}\
=\"never\";\n$PG{\"dalilite.pl\"}{\"mode\"}=\"expr\
esso,3dcoffee\";\n$PG{\"probconsRNA\"}{\"4_TCOFFEE\
\"}=\"PROBCONSRNA\";\n$PG{\"probconsRNA\"}{\"type\\
"}=\"RNA_multiple_aligner\";\n$PG{\"probconsRNA\"}\
{\"ADDRESS\"}=\"http://probcons.stanford.edu/\";\n\
$PG{\"probconsRNA\"}{\"language\"}=\"C++\";\n$PG{\\
"probconsRNA\"}{\"language2\"}=\"CXX\";\n$PG{\"pro\
bconsRNA\"}{\"source\"}=\"http://probcons.stanford\
.edu/probconsRNA.tar.gz\";\n$PG{\"probconsRNA\"}{\\
"mode\"}=\"mcoffee,rcoffee\";\n$PG{\"probconsRNA\"\
}{\"version\"}=\"1.0\";\n$PG{\"sfold\"}{\"4_TCOFFE\
E\"}=\"CONSAN\";\n$PG{\"sfold\"}{\"type\"}=\"RNA_p\
airwise_aligner\";\n$PG{\"sfold\"}{\"ADDRESS\"}=\"\
http://selab.janelia.org/software/consan/\";\n$PG{\
\"sfold\"}{\"language\"}=\"empty\";\n$PG{\"sfold\"\
}{\"language2\"}=\"empty\";\n$PG{\"sfold\"}{\"sour\
ce\"}=\"empty\";\n$PG{\"sfold\"}{\"update_action\"\
}=\"never\";\n$PG{\"sfold\"}{\"mode\"}=\"rcoffee\"\
;\n$PG{\"RNAplfold\"}{\"4_TCOFFEE\"}=\"RNAPLFOLD\"\
;\n$PG{\"RNAplfold\"}{\"type\"}=\"RNA_secondarystr\
ucture_predictor\";\n$PG{\"RNAplfold\"}{\"ADDRESS\\
"}=\"http://www.tbi.univie.ac.at/~ivo/RNA/\";\n$PG\
{\"RNAplfold\"}{\"language\"}=\"C\";\n$PG{\"RNAplf\
old\"}{\"language2\"}=\"C\";\n$PG{\"RNAplfold\"}{\\
"source\"}=\"http://www.tbi.univie.ac.at/~ivo/RNA/\
ViennaRNA-1.7.2.tar.gz\";\n$PG{\"RNAplfold\"}{\"mo\
de\"}=\"rcoffee,\";\n$PG{\"RNAplfold\"}{\"version\\
"}=\"1.7.2\";\n$PG{\"retree\"}{\"4_TCOFFEE\"}=\"PH\
YLIP\";\n$PG{\"retree\"}{\"type\"}=\"RNA_secondary\
structure_predictor\";\n$PG{\"retree\"}{\"ADDRESS\\
"}=\"http://evolution.gs.washington.edu/phylip/\";\
\n$PG{\"retree\"}{\"language\"}=\"C\";\n$PG{\"retr\
ee\"}{\"language2\"}=\"C\";\n$PG{\"retree\"}{\"sou\
rce\"}=\"http://evolution.gs.washington.edu/phylip\
/download/phylip-3.69.tar.gz\";\n$PG{\"retree\"}{\\
"mode\"}=\"trmsd,\";\n$PG{\"retree\"}{\"version\"}\
=\"3.69\";\n$PG{\"hmmtop\"}{\"4_TCOFFEE\"}=\"HMMTO\
P\";\n$PG{\"hmmtop\"}{\"type\"}=\"protein_secondar\
ystructure_predictor\";\n$PG{\"hmmtop\"}{\"ADDRESS\
\"}=\"www.enzim.hu/hmmtop/\";\n$PG{\"hmmtop\"}{\"l\
anguage\"}=\"C\";\n$PG{\"hmmtop\"}{\"language2\"}=\
\"C\";\n$PG{\"hmmtop\"}{\"source\"}=\"empty\";\n$P\
G{\"hmmtop\"}{\"update_action\"}=\"never\";\n$PG{\\
"hmmtop\"}{\"mode\"}=\"tcoffee\";\n$PG{\"gorIV\"}{\
\"4_TCOFFEE\"}=\"GOR4\";\n$PG{\"gorIV\"}{\"type\"}\
=\"protein_secondarystructure_predictor\";\n$PG{\"\
gorIV\"}{\"ADDRESS\"}=\"http://mig.jouy.inra.fr/lo\
giciels/gorIV/\";\n$PG{\"gorIV\"}{\"language\"}=\"\
C\";\n$PG{\"gorIV\"}{\"language2\"}=\"C\";\n$PG{\"\
gorIV\"}{\"source\"}=\"http://mig.jouy.inra.fr/log\
iciels/gorIV/GOR_IV.tar.gz\";\n$PG{\"gorIV\"}{\"up\
date_action\"}=\"never\";\n$PG{\"gorIV\"}{\"mode\"\
}=\"tcoffee\";\n$PG{\"wublast.pl\"}{\"4_TCOFFEE\"}\
=\"EBIWUBLASTc\";\n$PG{\"wublast.pl\"}{\"type\"}=\\
"protein_homology_predictor\";\n$PG{\"wublast.pl\"\
}{\"ADDRESS\"}=\"built_in\";\n$PG{\"wublast.pl\"}{\
\"ADDRESS2\"}=\"http://www.ebi.ac.uk/Tools/webserv\
ices/services/wublast\";\n$PG{\"wublast.pl\"}{\"la\
nguage\"}=\"Perl\";\n$PG{\"wublast.pl\"}{\"languag\
e2\"}=\"Perl\";\n$PG{\"wublast.pl\"}{\"source\"}=\\
"empty\";\n$PG{\"wublast.pl\"}{\"update_action\"}=\
\"never\";\n$PG{\"wublast.pl\"}{\"mode\"}=\"psicof\
fee,expresso,accurate\";\n$PG{\"blastpgp.pl\"}{\"4\
_TCOFFEE\"}=\"EBIBLASTPGPc\";\n$PG{\"blastpgp.pl\"\
}{\"type\"}=\"protein_homology_predictor\";\n$PG{\\
"blastpgp.pl\"}{\"ADDRESS\"}=\"built_in\";\n$PG{\"\
blastpgp.pl\"}{\"ADDRESS2\"}=\"http://www.ebi.ac.u\
k/Tools/webservices/services/blastpgp\";\n$PG{\"bl\
astpgp.pl\"}{\"language\"}=\"Perl\";\n$PG{\"blastp\
gp.pl\"}{\"language2\"}=\"Perl\";\n$PG{\"blastpgp.\
pl\"}{\"source\"}=\"empty\";\n$PG{\"blastpgp.pl\"}\
{\"update_action\"}=\"never\";\n$PG{\"blastpgp.pl\\
"}{\"mode\"}=\"psicoffee,expresso,accurate\";\n$PG\
{\"blastcl3\"}{\"4_TCOFFEE\"}=\"NCBIWEBBLAST\";\n$\
PG{\"blastcl3\"}{\"type\"}=\"protein_homology_pred\
ictor\";\n$PG{\"blastcl3\"}{\"ADDRESS\"}=\"ftp://f\
tp.ncbi.nih.gov/blast/executables/LATEST\";\n$PG{\\
"blastcl3\"}{\"language\"}=\"C\";\n$PG{\"blastcl3\\
"}{\"language2\"}=\"C\";\n$PG{\"blastcl3\"}{\"sour\
ce\"}=\"empty\";\n$PG{\"blastcl3\"}{\"update_actio\
n\"}=\"never\";\n$PG{\"blastcl3\"}{\"mode\"}=\"psi\
coffee,expresso,3dcoffee\";\n$PG{\"blastpgp\"}{\"4\
_TCOFFEE\"}=\"NCBIBLAST\";\n$PG{\"blastpgp\"}{\"ty\
pe\"}=\"protein_homology_predictor\";\n$PG{\"blast\
pgp\"}{\"ADDRESS\"}=\"ftp://ftp.ncbi.nih.gov/blast\
/executables/LATEST\";\n$PG{\"blastpgp\"}{\"langua\
ge\"}=\"C\";\n$PG{\"blastpgp\"}{\"language2\"}=\"C\
\";\n$PG{\"blastpgp\"}{\"source\"}=\"empty\";\n$PG\
{\"blastpgp\"}{\"update_action\"}=\"never\";\n$PG{\
\"blastpgp\"}{\"mode\"}=\"psicoffee,expresso,3dcof\
fee\";\n$PG{\"SOAP::Lite\"}{\"4_TCOFFEE\"}=\"SOAPL\
ITE\";\n$PG{\"SOAP::Lite\"}{\"type\"}=\"library\";\
\n$PG{\"SOAP::Lite\"}{\"ADDRESS\"}=\"http://cpanse\
arch.perl.org/src/MKUTTER/SOAP-Lite-0.710.08/Makef\
ile.PL\";\n$PG{\"SOAP::Lite\"}{\"language\"}=\"Per\
l\";\n$PG{\"SOAP::Lite\"}{\"language2\"}=\"Perl\";\
\n$PG{\"SOAP::Lite\"}{\"source\"}=\"empty\";\n$PG{\
\"blastpgp\"}{\"update_action\"}=\"never\";\n$PG{\\
"SOAP::Lite\"}{\"mode\"}=\"none\";\n$PG{\"XML::Sim\
ple\"}{\"4_TCOFFEE\"}=\"XMLSIMPLE\";\n$PG{\"XML::S\
imple\"}{\"type\"}=\"library\";\n$PG{\"XML::Simple\
\"}{\"ADDRESS\"}=\"http://search.cpan.org/~grantm/\
XML-Simple-2.18/lib/XML/Simple.pm\";\n$PG{\"XML::S\
imple\"}{\"language\"}=\"Perl\";\n$PG{\"XML::Simpl\
e\"}{\"language2\"}=\"Perl\";\n$PG{\"XML::Simple\"\
}{\"source\"}=\"empty\";\n$PG{\"XML::Simple\"}{\"m\
ode\"}=\"psicoffee,expresso,accurate\";\n$MODE{\"t\
coffee\"}{\"name\"}=\"tcoffee\";\n$MODE{\"rcoffee\\
"}{\"name\"}=\"rcoffee\";\n$MODE{\"3dcoffee\"}{\"n\
ame\"}=\"3dcoffee\";\n$MODE{\"mcoffee\"}{\"name\"}\
=\"mcoffee\";\n$MODE{\"expresso\"}{\"name\"}=\"exp\
resso\";\n$MODE{\"trmsd\"}{\"name\"}=\"trmsd\";\n$\
MODE{\"accurate\"}{\"name\"}=\"accurate\";\n$MODE{\
\"seq_reformat\"}{\"name\"}=\"seq_reformat\";\n\n\\
n$PG{C}{compiler}=\"gcc\";\n$PG{C}{compiler_flag}=\
\"CC\";\n$PG{C}{options}=\"\";\n$PG{C}{options_fla\
g}=\"CFLAGS\";\n$PG{C}{type}=\"compiler\";\n\n$PG{\
\"CXX\"}{compiler}=\"g++\";\n$PG{\"CXX\"}{compiler\
_flag}=\"CXX\";\n$PG{\"CXX\"}{options}=\"\";\n$PG{\
\"CXX\"}{options_flag}=\"CXXFLAGS\";\n$PG{CXX}{typ\
e}=\"compiler\";\n\n$PG{\"CPP\"}{compiler}=\"g++\"\
;\n$PG{\"CPP\"}{compiler_flag}=\"CPP\";\n$PG{\"CPP\
\"}{options}=\"\";\n$PG{\"CPP\"}{options_flag}=\"C\
PPFLAGS\";\n$PG{CPP}{type}=\"compiler\";\n\n$PG{\"\
GPP\"}{compiler}=\"g++\";\n$PG{\"GPP\"}{compiler_f\
lag}=\"GPP\";\n$PG{\"GPP\"}{options}=\"\";\n$PG{\"\
GPP\"}{options_flag}=\"CFLAGS\";\n$PG{GPP}{type}=\\
"compiler\";\n\n$PG{Fortran}{compiler}=\"g77\";\n$\
PG{Fortran}{compiler_flag}=\"FCC\";\n$PG{Fortran}{\
type}=\"compiler\";\n\n$PG{Perl}{compiler}=\"CPAN\\
";\n$PG{Perl}{type}=\"compiler\";\n\n$SUPPORTED_OS\
{macox}=\"Macintosh\";\n$SUPPORTED_OS{linux}=\"Lin\
ux\";\n$SUPPORTED_OS{windows}=\"Cygwin\";\n\n\n\n$\
MODE{t_coffee}{description}=\" for regular multipl\
e sequence alignments\";\n$MODE{rcoffee} {descript\
ion}=\" for RNA multiple sequence alignments\";\n\\
n$MODE{psicoffee} {description}=\" for Homology Ex\
tended multiple sequence alignments\";\n$MODE{expr\
esso}{description}=\" for very accurate structure \
based multiple sequence alignments\";\n$MODE{\"3dc\
offee\"}{description}=\" for multiple structure al\
ignments\";\n$MODE{mcoffee} {description}=\" for c\
ombining alternative multiple sequence alignment p\
ackages\\n------- into a unique meta-package. The \
installer will upload several MSA packages and com\
pile them\\n\n\";\n\n\n&post_process_PG();\nreturn\
;\n}\n\nsub post_process_PG\n  {\n    my $p;\n    \
\n    %PG=&name2dname (%PG);\n    %MODE=&name2dnam\
e(%MODE);\n    foreach $p (keys(%PG)){if ( $PG{$p}\
{type} eq \"compiler\"){$PG{$p}{update_action}=\"n\
ever\";}}\n    \n  }\n\nsub name2dname\n  {\n    m\
y (%L)=(@_);\n    my ($l, $ml);\n    \n    foreach\
 my $pg (keys(%L))\n      {\n	$l=length ($pg);\n	i\
f ( $l>$ml){$ml=$l;}\n      }\n    $ml+=1;\n    fo\
reach my $pg (keys(%L))\n      {\n	my $name;\n	$l=\
$ml-length ($pg);\n	$name=$pg;\n	for ( $b=0; $b<$l\
; $b++)\n	  {\n	    $name .=\" \";\n	  }\n	$L{$pg}\
{dname}=$name;\n      }\n    return %L;\n  }\n\nsu\
b env_file2putenv\n  {\n    my $f=@_[0];\n    my $\
F=new FileHandle;\n    my $n;\n    \n    open ($F,\
 \"$f\");\n    while (<$F>)\n      {\n	my $line=$_\
;\n	my($var, $value)=($_=~/(\\S+)\\=(\\S*)/);\n	$E\
NV{$var}=$value;\n	$ENV_SET{$var}=1;\n	$n++;\n    \
  }\n    close ($F);\n    return $n;\n  }\n\nsub r\
eplace_line_in_file\n  {\n    my ($file, $wordin, \
$wordout)=@_;\n    my $O=new FileHandle;\n    my $\
I=new FileHandle;\n    my $l;\n    if (!-e $file){\
return;}\n    \n    system (\"mv $file $file.old\"\
);\n    open ($O, \">$file\");\n    open ($I, \"$f\
ile.old\");\n    while (<$I>)\n      {\n	$l=$_;\n	\
if (!($l=~/$wordin/)){print $O \"$l\";}\n	elsif ( \
$wordout ne \"\"){$l=~s/$wordin/$wordout/g;print $\
O \"$l\";}\n      }\n    close ($O);\n    close ($\
I);\n    return;\n  }\n\nsub add_C_libraries\n  {\\
n   my ($file,$first,@list)=@_;\n   \n    my $O=ne\
w FileHandle;\n    my $I=new FileHandle;\n    my (\
$l,$anchor);\n    if (!-e $file){return;}\n   \n  \
  $anchor=\"#include <$first>\";\n	 \n    system (\
\"mv $file $file.old\");\n    open ($O, \">$file\"\
);\n    open ($I, \"$file.old\");\n    while (<$I>\
)\n      {\n	$l=$_;\n	print $O \"$l\";\n	if (!($l=\
~/$anchor/))\n	   {\n	    \n	    foreach my $lib (\
@list)\n	       {\n                  print $O \"#i\
nclude <$lib>\\n\";\n	       }\n           }\n    \
  }\n    close ($O);\n    close ($I);\n    return;\
\n    }\n","use Env;\nuse Cwd;\n@suffix=(\"tmp\", \
\"temp\", \"cache\", \"t_coffee\", \"core\", \"tco\
ffee\");\n\nif ($#ARGV==-1)\n  {\n    print \"clea\
n_cache.pl -file <file to add in -dir> -dir=<dir> \
-size=<value in Mb>\\n0: unlimited -1 always.\\nWi\
ll only clean directories matching:[\";\n    forea\
ch $k(@suffix){print \"*$k* \";}\n    print \"]\\n\
\";\n    exit (EXIT_FAILURE);\n  }\n\n$cl=join (\"\
 \",@ARGV);\nif (($cl=~/\\-no_action/))\n  {\n    \
exit (EXIT_SUCCESS);\n  }\n\nif (($cl=~/\\-debug/)\
)\n  {\n    $DEBUG=1;\n  }\nelse\n  {\n    $DEBUG=\
0;\n  }\n\nif (($cl=~/\\-dir=(\\S+)/))\n  {\n    $\
dir=$1;\n  }\nelse\n  {\n    $dir=\"./\";\n  }\n\n\
if ($cl=~/\\-file=(\\S+)/)\n  {\n    $file=$1;\n  \
}\nelse\n  {\n    $file=0;\n  }\n\nif ($cl=~/\\-si\
ze=(\\S+)/)\n  {\n    $max_size=$1;\n  }\nelse\n  \
{\n    $max_size=0;#unlimited\n  }\nif ($cl=~/\\-f\
orce/)\n  {\n    $force=1;\n  }\nelse\n  {\n    $f\
orce=0;\n  }\n\nif ($cl=~/\\-age=(\\S+)/)\n  {\n  \
  $max_age=$1;\n  }\nelse\n  {\n    $max_age=0;#un\
limited\n  }\n\n$max_size*=1000000;\nif ( ! -d $di\
r)\n  {\n    print STDERR \"\\nCannot process $dir\
: does not exist \\n\";\n    exit (EXIT_FAILURE);\\
n  }\n\nif ( !($dir=~/^\\//))\n  {\n    $base=cwd(\
);\n    $dir=\"$base/$dir\";\n  }\n\n$proceed=0;\n\
foreach $s (@suffix)\n  {\n    \n    if (($dir=~/$\
s/)){$proceed=1;}\n    $s=uc ($s);\n    if (($dir=\
~/$s/)){$proceed=1;}\n  }\nif ( $proceed==0)\n  {\\
n    print STDERR \"Clean_cache.pl can only clean \
directories whose absolute path name contains the \
following strings:\";\n    foreach $w (@suffix) {p\
rint STDERR \"$w \";$w=lc($w); print STDERR \"$w \\
";}\n    print STDERR \"\\nCannot process $dir\\n\\
";\n    exit (EXIT_FAILURE);\n  }\n\n$name_file=\"\
$dir/name_file.txt\";\n$size_file=\"$dir/size_file\
.txt\";\nif ( $force){&create_ref_file ($dir,$name\
_file,$size_file);}\nif ($file){&add_file ($dir, $\
name_file, $size_file, $file);}\n&clean_dir ($dir,\
 $name_file, $size_file, $max_size,$max_age);\nexi\
t (EXIT_SUCCESS);\n\nsub clean_dir \n  {\n    my (\
$dir, $name_file, $size_file, $max_size, $max_age)\
=@_;\n    my ($tot_size, $size, $f, $s);\n\n  \n  \
  $tot_size=&get_tot_size ($dir, $name_file, $size\
_file);\n\n    if ( $tot_size<=$max_size){return ;\
}\n    else {$max_size/=2;}\n    \n    #recreate t\
he name file in case some temprary files have not \
been properly registered\n    &create_ref_file ($d\
ir, $name_file, $size_file, $max_age);\n  \n    $n\
ew_name_file=&vtmpnam();\n    open (R, \"$name_fil\
e\");\n    open (W, \">$new_name_file\");\n    whi\
le (<R>)\n      {\n	my $line=$_;\n	\n	($f, $s)=($l\
ine=~/(\\S+) (\\S+)/);\n	if ( !($f=~/\\S/)){next;}\
\n	\n	elsif ($max_size && $tot_size>=$max_size && \
!($f=~/name_file/))\n	  {\n	    remove ( \"$dir/$f\
\");\n	    $tot_size-=$s;\n	  }\n	elsif ( $max_age\
 && -M(\"$dir/$f\")>=$max_age)\n	  {\n	    remove \
( \"$dir/$f\");\n	    $tot_size-=$s;\n	  }\n	else\\
n	  {\n	    print W \"$f $s\\n\";\n	  }\n      }\n\
    close (R);\n    close (W);\n    open (F, \">$s\
ize_file\");\n    print F \"$tot_size\";\n    if (\
 -e $new_name_file){`mv $new_name_file $name_file`\
;}\n    close (F);\n  }\nsub get_tot_size\n  {\n  \
  my ($dir, $name_file, $size_file)=@_;\n    my $s\
ize;\n    \n    if ( !-d $dir){return 0;}\n    if \
( !-e $name_file)\n      {\n	\n	&create_ref_file (\
$dir, $name_file, $size_file);\n      }\n    open \
(F, \"$size_file\");\n    $size=<F>;\n    close (F\
);\n    chomp ($size);\n    return $size;\n  }\nsu\
b size \n  {\n    my $f=@_[0];\n\n    if ( !-d $f)\
{return -s($f);}\n    else {return &dir2size($f);}\
\n  }\nsub dir2size\n  {\n    my $d=@_[0];\n    my\
 ($s, $f);\n    \n    if ( !-d $d) {return 0;}\n  \
  \n    foreach $f (&dir2list ($d))\n      {\n	if \
( -d $f){$s+=&dir2size (\"$d/$f\");}\n	else {$s+= \
-s \"$dir/$f\";}\n      }\n    return $s;\n  }\n\n\
sub remove \n  {\n    my $file=@_[0];\n    my ($f)\
;\n    \n    debug_print( \"--- $file ---\\n\");\n\
    if (($file eq \".\") || ($file eq \"..\") || (\
$file=~/\\*/)){return EXIT_FAILURE;}\n    elsif ( \
!-d $file)\n      {\n	debug_print (\"unlink $file\\
\n\");\n	if (-e $file){unlink ($file);}\n      }\n\
    elsif ( -d $file)\n      {\n	debug_print (\"++\
++++++ $file +++++++\\n\");\n	foreach $f (&dir2lis\
t($file))\n	  {\n	    &remove (\"$file/$f\");\n	  \
}\n	debug_print (\"rmdir $file\\n\");\n	rmdir $fil\
e;\n      }\n    else\n      {\n	debug_print (\"??\
??????? $file ????????\\n\");\n      }\n    return\
 EXIT_SUCCESS;\n  }\n\nsub dir2list\n  {\n    my $\
dir=@_[0];\n    my (@list1, @list2,@list3, $l);\n\\
n    opendir (DIR,$dir);\n    @list1=readdir (DIR)\
;\n    closedir (DIR);\n    \n    foreach $l (@lis\
t1)\n      {\n	if ( $l ne \".\" && $l ne \"..\"){@\
list2=(@list2, $l);}\n      }\n    @list3 = sort {\
 (-M \"$dir/$list2[$b]\") <=> (-M \"$dir/$list2[$a\
]\")} @list2;\n    return @list3;\n    \n  }\n\nsu\
b debug_print\n  {\n    \n    if ($DEBUG==1){print\
 @_;}\n    \n  }\nsub create_ref_file\n  {\n    my\
 ($dir,$name_file,$size_file)=@_;\n    my ($f, $s,\
 $tot_size, @l);\n    \n    if ( !-d $dir){return;\
}\n    \n    @l=&dir2list ($dir);\n    open (F, \"\
>$name_file\");\n    foreach $f (@l)\n      {\n	$s\
=&size(\"$dir/$f\");\n	$tot_size+=$s;\n	print F \"\
$f $s\\n\";\n      }\n    &myecho ($tot_size, \">$\
size_file\");\n    close (F);\n  }\nsub add_file \\
n  {\n    my ($dir,$name_file,$size_file,$file)=@_\
;\n    my ($s, $tot_size);\n    \n    if ( !-d $di\
r)   {return;}\n    if ( !-e \"$dir/$file\" ) {ret\
urn;}\n    if ( !-e $name_file){&create_ref_file (\
$dir,$name_file,$size_file);}\n					    \n    $s=&\
size(\"$dir/$file\");\n    open (F, \">>$name_file\
\");\n    print F \"$file\\n\";\n    close (F);\n\\
n    $tot_size=&get_tot_size ($dir,$name_file,$siz\
e_file);\n    $tot_size+=$s;\n    &myecho ($tot_si\
ze, \">$size_file\");\n    \n  }\n	\nsub myecho\n \
 {\n    my ($string, $file)=@_;\n    open (ECHO, $\
file) || die;\n    print ECHO \"$string\";\n    cl\
ose (ECHO);\n  }\n    \n		\n	\nsub vtmpnam\n  {\n \
   my $tmp_file_name;\n    $tmp_name_counter++;\n \
   $tmp_file_name=\"tmp_file_for_clean_cache_pdb$$\
.$tmp_name_counter\";\n    $tmp_file_list[$ntmp_fi\
le++]=$tmp_file_name;\n    if ( -e $tmp_file_name)\
 {return &vtmpnam ();}\n    else {return $tmp_file\
_name;}\n  }\n","\nmy $address=\"http://www.tcoffe\
e.org/Projects/Datasets/NatureMethodsDataset.tar.g\
z\";\n&url2file ($address);\n\nif ( -e \"NatureMet\
hodsDataset.tar.gz\" )\n  {\n    \n    system (\"g\
unzip NatureMethodsDataset.tar.gz\");\n    system \
(\"tar -xvf NatureMethodsDataset.tar\");\n    \n  \
  print \"Your Data Set is in the Folder NatureMet\
hodsDataset\\n\";\n  }\nelse \n  {\n    print \"Co\
uld not Download Dataset --- Web site may be down \
-- Try again later\\n\";\n  }\n\n\n\n\nsub url2fil\
e\n{\n    my ($address, $out, $wget_arg, $curl_arg\
)=(@_);\n    my ($pg, $flag, $r, $arg, $count);\n \
   \n    if (!$CONFIGURATION){&check_configuration\
 (\"wget\", \"INTERNET\", \"gzip\");$CONFIGURATION\
=1;}\n    \n    if (&pg_is_installed (\"wget\"))  \
 {$pg=\"wget\"; $flag=\"-O\";$arg=$wget_arg;}\n   \
 elsif (&pg_is_installed (\"curl\")){$pg=\"curl\";\
 $flag=\"-o\";$arg=$curl_arg;}\n    return system \
(\"$pg $address >/dev/null 2>/dev/null\");\n\n}\n\\
nsub pg_is_installed\n  {\n    my @ml=@_;\n    my \
$r, $p, $m;\n    my $supported=0;\n    \n    my $p\
=shift (@ml);\n    if ($p=~/::/)\n      {\n	if (sy\
stem (\"perl -M$p -e 1\")==$EXIT_SUCCESS){return 1\
;}\n	else {return 0;}\n      }\n    else\n      {\\
n	$r=`which $p 2>/dev/null`;\n	if ($r eq \"\"){ret\
urn 0;}\n	else {return 1;}\n      }\n  }\nsub chec\
k_configuration \n    {\n      my @l=@_;\n      my\
 $v;\n      foreach my $p (@l)\n	{\n	  \n	  if   (\
 $p eq \"EMAIL\")\n	    { \n	      if ( !($EMAIL=~\
/@/))\n		{\n		  exit (EXIT_FAILURE);\n		}\n	    }\\
n	  elsif( $p eq \"INTERNET\")\n	    {\n	      if \
( !&check_internet_connection())\n		{\n		  exit (E\
XIT_FAILURE);\n		}\n	    }\n	  elsif( $p eq \"wget\
\")\n	    {\n	      if (!&pg_is_installed (\"wget\\
") && !&pg_is_installed (\"curl\"))\n		{\n		  exit\
 (EXIT_FAILURE);\n		}\n	    }\n	  elsif( !(&pg_is_\
installed ($p)))\n	    {\n	      exit (EXIT_FAILUR\
E);\n	    }\n	}\n      return 1;\n    }\nsub check\
_internet_connection\n  {\n    my $internet;\n    \
my $tmp;\n    &check_configuration ( \"wget\"); \n\
    \n    $tmp=&vtmpnam ();\n    \n    if     (&pg\
_is_installed    (\"wget\")){`wget www.google.com \
-O$tmp >/dev/null 2>/dev/null`;}\n    elsif  (&pg_\
is_installed    (\"curl\")){`curl www.google.com -\
o$tmp >/dev/null 2>/dev/null`;}\n    \n    if ( !-\
e $tmp || -s $tmp < 10){$internet=0;}\n    else {$\
internet=1;}\n    if (-e $tmp){unlink $tmp;}\n\n  \
  return $internet;\n  }\n\nsub vtmpnam\n      {\n\
	my $r=rand(100000);\n	my $f=\"file.$r.$$\";\n	whi\
le (-e $f)\n	  {\n	    $f=vtmpnam();\n	  }\n	push \
(@TMPFILE_LIST, $f);\n	return $f;\n      }\n","\n$\
t_coffee=\"t_coffee\";\n\nforeach $value ( @ARGV)\\
n  {\n    $seq_file=$seq_file.\" \".$value;\n  }\n\
\n$name=$ARGV[0];\n$name=~s/\\.[^\\.]*$//;\n$lib_n\
ame=\"$name.mocca_lib\";\n$type=`t_coffee $seq_fil\
e -get_type -quiet`;\nchop ($type);\n\nif ( $type \
eq \"PROTEIN\"){$lib_mode=\"lalign_rs_s_pair -lali\
gn_n_top 20\";}\nelsif ( $type eq\"DNA\"){$lib_mod\
e=\"lalign_rs_s_dna_pair -lalign_n_top 40\";}\n\ni\
f ( !(-e $lib_name))\n  {\n	  \n  $command=\"$t_co\
ffee -mocca -seq_weight=no -cosmetic_penalty=0 -mo\
cca_interactive -in $lib_mode -out_lib $lib_name -\
infile $seq_file\";\n  \n  }\nelsif ( (-e $lib_nam\
e))\n  {\n  $command=\"$t_coffee -mocca -seq_weigh\
t=no -cosmetic_penalty=0 -mocca_interactive -in $l\
ib_name -infile $seq_file\";\n  \n  }\n\nsystem ($\
command);\n\nexit;\n\n","my $WSDL = 'http://www.eb\
i.ac.uk/Tools/webservices/wsdl/WSDaliLite.wsdl';\n\
\nuse SOAP::Lite;\nuse Data::Dumper;\nuse Getopt::\
Long qw(:config no_ignore_case bundling);\nuse Fil\
e::Basename;\n\nmy $checkInterval = 5;\n\nmy %para\
ms=(\n	    'async' => '1', # Use async mode and si\
mulate sync mode in client\n	    );\nGetOptions(\n\
    'pdb1=s'     => \\$params{'sequence1'},\n    '\
chainid1=s' => \\$params{'chainid1'},\n    'pdb2=s\
'     => \\$params{'sequence2'},\n    'chainid2=s'\
 => \\$params{'chainid2'},\n    \"help|h\"	 => \\$\
help, # Usage info\n    \"async|a\"	 => \\$async, \
# Asynchronous submission\n    \"polljob\"	 => \\$\
polljob, # Get results\n    \"status\"	 => \\$stat\
us, # Get status\n    \"jobid|j=s\"  => \\$jobid, \
# JobId\n    \"email|S=s\"  => \\$params{email}, #\
 E-mail address\n    \"trace\"      => \\$trace, #\
 SOAP messages\n    \"sequence=s\" => \\$sequence,\
 # Input PDB\n    );\n\nmy $scriptName = basename(\
$0, ());\nif($help) {\n    &usage();\n    exit(0);\
\n}\n\nif($trace) {\n    print \"Tracing active\\n\
\";\n    SOAP::Lite->import(+trace => 'debug');\n}\
\n\nmy $soap = SOAP::Lite\n    ->service($WSDL)\n \
   ->on_fault(sub {\n        my $soap = shift;\n  \
      my $res = shift;\n        # Throw an excepti\
on for all faults\n        if(ref($res) eq '') {\n\
            die($res);\n        } else {\n        \
    die($res->faultstring);\n        }\n        re\
turn new SOAP::SOM;\n    }\n               );\n\ni\
f( !($polljob || $status) &&\n    !( defined($para\
ms{'sequence1'}) && defined($params{'sequence2'}) \
)\n    ) {\n    print STDERR 'Error: bad option co\
mbination', \"\\n\";\n    &usage();\n    exit(1);\\
n}\nelsif($polljob && defined($jobid)) {\n    prin\
t \"Getting results for job $jobid\\n\";\n    getR\
esults($jobid);\n}\nelsif($status && defined($jobi\
d)) {\n    print STDERR \"Getting status for job $\
jobid\\n\";\n    my $result = $soap->checkStatus($\
jobid);\n    print STDOUT \"$result\", \"\\n\";\n \
   if($result eq 'DONE') {\n	print STDERR \"To get\
 results: $scriptName --polljob --jobid $jobid\\n\\
";\n    }\n}\nelse {\n    if(-f $params{'sequence1\
'}) {\n	$params{'sequence1'} = read_file($params{'\
sequence1'});\n    }\n    if(-f $params{'sequence2\
'}) {\n	$params{'sequence2'} = read_file($params{'\
sequence2'});\n    }\n\n    my $jobid;\n    my $pa\
ramsData = SOAP::Data->name('params')->type(map=>\\
\%params);\n    # For SOAP::Lite 0.60 and earlier \
parameters are passed directly\n    if($SOAP::Lite\
::VERSION eq '0.60' || $SOAP::Lite::VERSION =~ /0\\
\.[1-5]/) {\n        $jobid = $soap->runDaliLite($\
paramsData);\n    }\n    # For SOAP::Lite 0.69 and\
 later parameter handling is different, so pass\n \
   # undef's for templated params, and then pass t\
he formatted args.\n    else {\n        $jobid = $\
soap->runDaliLite(undef,\n				     $paramsData);\n\
    }\n\n    if (defined($async)) {\n	print STDOUT\
 $jobid, \"\\n\";\n        print STDERR \"To check\
 status: $scriptName --status --jobid $jobid\\n\";\
\n    } else { # Synchronous mode\n        print S\
TDERR \"JobId: $jobid\\n\";\n        sleep 1;\n   \
     getResults($jobid);\n    }\n}\n\nsub clientPo\
ll($) {\n    my $jobid = shift;\n    my $result = \
'PENDING';\n    # Check status and wait if not fin\
ished\n    #print STDERR \"Checking status: $jobid\
\\n\";\n    while($result eq 'RUNNING' || $result \
eq 'PENDING') {\n        $result = $soap->checkSta\
tus($jobid);\n        print STDERR \"$result\\n\";\
\n        if($result eq 'RUNNING' || $result eq 'P\
ENDING') {\n            # Wait before polling agai\
n.\n            sleep $checkInterval;\n        }\n\
    }\n}\n\nsub getResults($) {\n    $jobid = shif\
t;\n    # Check status, and wait if not finished\n\
    clientPoll($jobid);\n    # Use JobId if output\
 file name is not defined\n    unless(defined($out\
file)) {\n        $outfile=$jobid;\n    }\n    # G\
et list of data types\n    my $resultTypes = $soap\
->getResults($jobid);\n    # Get the data and writ\
e it to a file\n    if(defined($outformat)) { # Sp\
ecified data type\n        my $selResultType;\n   \
     foreach my $resultType (@$resultTypes) {\n   \
         if($resultType->{type} eq $outformat) {\n\
                $selResultType = $resultType;\n   \
         }\n        }\n        $res=$soap->poll($j\
obid, $selResultType->{type});\n        write_file\
($outfile.'.'.$selResultType->{ext}, $res);\n    }\
 else { # Data types available\n        # Write a \
file for each output type\n        for my $resultT\
ype (@$resultTypes){\n            #print \"Getting\
 $resultType->{type}\\n\";\n            $res=$soap\
->poll($jobid, $resultType->{type});\n            \
write_file($outfile.'.'.$resultType->{ext}, $res);\
\n        }\n    }\n}\n\nsub read_file($) {\n    m\
y $filename = shift;\n    open(FILE, $filename);\n\
    my $content;\n    my $buffer;\n    while(sysre\
ad(FILE, $buffer, 1024)) {\n	$content.= $buffer;\n\
    }\n    close(FILE);\n    return $content;\n}\n\
\nsub write_file($$) {\n    my ($tmp,$entity) = @_\
;\n    print STDERR \"Creating result file: \".$tm\
p.\"\\n\";\n    unless(open (FILE, \">$tmp\")) {\n\
	return 0;\n    }\n    syswrite(FILE, $entity);\n \
   close (FILE);\n    return 1;\n}\n\nsub usage {\\
n    print STDERR <<EOF\nDaliLite\n========\n\nPai\
rwise comparison of protein structures\n\n[Require\
d]\n\n  --pdb1                : str  : PDB ID for \
structure 1\n  --pdb2                : str  : PDB \
ID for structure 2\n\n[Optional]\n\n  --chain1    \
          : str  : Chain identifer in structure 1\\
n  --chain2              : str  : Chain identifer \
in structure 2\n\n[General]\n\n  -h, --help       \
     :      : prints this help text\n  -S, --email\
           : str  : user email address\n  -a, --as\
ync           :      : asynchronous submission\n  \
    --status          :      : poll for the status\
 of a job\n      --polljob         :      : poll f\
or the results of a job\n  -j, --jobid           :\
 str  : jobid for an asynchronous job\n  -O, --out\
file         : str  : file name for results (defau\
lt is jobid)\n      --trace	        :      : show \
SOAP messages being interchanged \n\nSynchronous j\
ob:\n\n  The results/errors are returned as soon a\
s the job is finished.\n  Usage: $scriptName --ema\
il <your\\@email> [options] pdbFile [--outfile str\
ing]\n  Returns: saves the results to disk\n\nAsyn\
chronous job:\n\n  Use this if you want to retriev\
e the results at a later time. The results \n  are\
 stored for up to 24 hours. \n  The asynchronous s\
ubmission mode is recommended when users are submi\
tting \n  batch jobs or large database searches	\n\
  Usage: $scriptName --email <your\\@email> --asyn\
c [options] pdbFile\n  Returns: jobid\n\n  Use the\
 jobid to query for the status of the job. \n  Usa\
ge: $scriptName --status --jobid <jobId>\n  Return\
s: string indicating the status of the job:\n    D\
ONE - job has finished\n    RUNNING - job is runni\
ng\n    NOT_FOUND - job cannot be found\n    ERROR\
 - the jobs has encountered an error\n\n  When don\
e, use the jobid to retrieve the status of the job\
. \n  Usage: $scriptName --polljob --jobid <jobId>\
 [--outfile string]\n\n[Help]\n\n  For more detail\
ed help information refer to\n  http://www.ebi.ac.\
uk/DaliLite/\nEOF\n;\n}\n","my $WSDL = 'http://www\
.ebi.ac.uk/Tools/webservices/wsdl/WSWUBlast.wsdl';\
\n\nuse strict;\nuse SOAP::Lite;\nuse Getopt::Long\
 qw(:config no_ignore_case bundling);\nuse File::B\
asename;\n\nmy $checkInterval = 15;\n\nmy $numOpts\
 = scalar(@ARGV);\nmy ($outfile, $outformat, $help\
, $async, $polljob, $status, $ids, $jobid, $trace,\
 $sequence);\nmy %params= ( # Defaults\n	      'as\
ync' => 1, # Force into async mode\n	      'exp' =\
> 10.0, # E-value threshold\n	      'numal' => 50,\
 # Maximum number of alignments\n	      'scores' =\
> 100, # Maximum number of scores\n            );\\
nGetOptions( # Map the options into variables\n   \
 \"program|p=s\"     => \\$params{program}, # BLAS\
T program\n    \"database|D=s\"    => \\$params{da\
tabase}, # Search database\n    \"matrix|m=s\"    \
  => \\$params{matrix}, # Scoring matrix\n    \"ex\
p|E=f\"         => \\$params{exp}, # E-value thres\
hold\n    \"echofilter|e\"    => \\$params{echofil\
ter}, # Display filtered sequence\n    \"filter|f=\
s\"      => \\$params{filter}, # Low complexity fi\
lter name\n    \"alignments|b=i\"  => \\$params{nu\
mal}, # Number of alignments\n    \"scores|s=i\"  \
    => \\$params{scores}, # Number of scores\n    \
\"sensitivity|S=s\" => \\$params{sensitivity}, # S\
earch sensitivity\n    \"sort|t=s\"	      => \\$pa\
rams{sort}, # Sort hits by...\n    \"stats|T=s\"  \
     => \\$params{stats}, # Scoring statistic to u\
se\n    \"strand|d=s\"      => \\$params{strand}, \
# Strand to use in DNA vs. DNA search\n    \"topco\
mbon|c=i\"   => \\$params{topcombon}, # Consistent\
 sets of HSPs\n    \"outfile=s\"       => \\$outfi\
le, # Output file\n    \"outformat|o=s\"   => \\$o\
utformat, # Output format\n    \"help|h\"	      =>\
 \\$help, # Usage info\n    \"async|a\"	      => \\
\$async, # Asynchronous mode\n    \"polljob\"	    \
  => \\$polljob, # Get results\n    \"status\"	   \
   => \\$status, # Get job status\n    \"ids\"    \
         => \\$ids, # Get ids from result\n    \"j\
obid|j=s\"       => \\$jobid, # JobId\n    \"email\
=s\"         => \\$params{email}, # E-mail address\
\n    \"trace\"           => \\$trace, # SOAP trac\
e\n    \"sequence=s\"      => \\$sequence, # Query\
 sequence\n    );\n\nmy $scriptName = basename($0,\
 ());\nif($help || $numOpts == 0) {\n    &usage();\
\n    exit(0);\n}\n\nif($trace){\n    print STDERR\
 \"Tracing active\\n\";\n    SOAP::Lite->import(+t\
race => 'debug');\n}\n\nmy $soap = SOAP::Lite\n   \
 ->service($WSDL)\n    ->proxy('http://localhost/'\
,\n    #proxy => ['http' => 'http://your.proxy.ser\
ver/'], # HTTP proxy\n    timeout => 600, # HTTP c\
onnection timeout\n    )\n    ->on_fault(sub { # S\
OAP fault handler\n        my $soap = shift;\n    \
    my $res = shift;\n        # Throw an exception\
 for all faults\n        if(ref($res) eq '') {\n  \
          die($res);\n        } else {\n          \
  die($res->faultstring);\n        }\n        retu\
rn new SOAP::SOM;\n    }\n               );\n\nif(\
 !($polljob || $status || $ids) &&\n    !( defined\
($ARGV[0]) || defined($sequence) )\n    ) {\n    p\
rint STDERR 'Error: bad option combination', \"\\n\
\";\n    &usage();\n    exit(1);\n}\nelsif($polljo\
b && defined($jobid)) {\n    print \"Getting resul\
ts for job $jobid\\n\";\n    getResults($jobid);\n\
}\nelsif($status && defined($jobid)) {\n    print \
STDERR \"Getting status for job $jobid\\n\";\n    \
my $result = $soap->checkStatus($jobid);\n    prin\
t STDOUT \"$result\\n\";\n    if($result eq 'DONE'\
) {\n	print STDERR \"To get results: $scriptName -\
-polljob --jobid $jobid\\n\";\n    }\n}  \nelsif($\
ids && defined($jobid)) {\n    print STDERR \"Gett\
ing ids from job $jobid\\n\";\n    getIds($jobid);\
\n}\nelse {\n    # Prepare input data\n    my $con\
tent;\n    my (@contents) = ();\n    if(-f $ARGV[0\
] || $ARGV[0] eq '-') {	\n	$content={type=>'sequen\
ce',content=>read_file($ARGV[0])};	\n    }\n    if\
($sequence) {	\n	if(-f $sequence || $sequence eq '\
-') {	\n	    $content={type=>'sequence',content=>r\
ead_file($ARGV[0])};	\n	} else {\n	    $content={t\
ype=>'sequence',content=>$sequence};\n	}\n    }\n \
   push @contents, $content;\n\n    # Submit the j\
ob\n    my $paramsData = SOAP::Data->name('params'\
)->type(map=>\\%params);\n    my $contentData = SO\
AP::Data->name('content')->value(\\@contents);\n  \
  # For SOAP::Lite 0.60 and earlier parameters are\
 passed directly\n    if($SOAP::Lite::VERSION eq '\
0.60' || $SOAP::Lite::VERSION =~ /0\\.[1-5]/) {\n \
       $jobid = $soap->runWUBlast($paramsData, $co\
ntentData);\n    }\n    # For SOAP::Lite 0.69 and \
later parameter handling is different, so pass\n  \
  # undef's for templated params, and then pass th\
e formatted args.\n    else {\n        $jobid = $s\
oap->runWUBlast(undef, undef,\n				   $paramsData,\
 $contentData);\n    }\n\n    # Asynchronous mode:\
 output jobid and exit.\n    if (defined($async)) \
{\n	print STDOUT $jobid, \"\\n\";\n        print S\
TDERR \"To check status: $scriptName --status --jo\
bid $jobid\\n\";\n    }\n    # Synchronous mode: t\
ry to get results\n    else {\n        print STDER\
R \"JobId: $jobid\\n\";\n        sleep 1;\n       \
 getResults($jobid);\n    }\n}\n\nsub getIds($) {\\
n    my $jobid = shift;\n    my $results = $soap->\
getIds($jobid);\n    for my $result (@$results){\n\
	print \"$result\\n\";\n    }\n}\n\nsub clientPoll\
($) {\n    my $jobid = shift;\n    my $result = 'P\
ENDING';\n    # Check status and wait if not finis\
hed\n    while($result eq 'RUNNING' || $result eq \
'PENDING') {\n        $result = $soap->checkStatus\
($jobid);\n        print STDERR \"$result\\n\";\n \
       if($result eq 'RUNNING' || $result eq 'PEND\
ING') {\n            # Wait before polling again.\\
n            sleep $checkInterval;\n        }\n   \
 }\n}\n\nsub getResults($) {\n    my $jobid = shif\
t;\n    my $res;\n    # Check status, and wait if \
not finished\n    clientPoll($jobid);\n    # Use J\
obId if output file name is not defined\n    unles\
s(defined($outfile)) {\n        $outfile=$jobid;\n\
    }\n    # Get list of data types\n    my $resul\
tTypes = $soap->getResults($jobid);\n    # Get the\
 data and write it to a file\n    if(defined($outf\
ormat)) { # Specified data type\n	if($outformat eq\
 'xml') {$outformat = 'toolxml';}\n	if($outformat \
eq 'txt') {$outformat = 'tooloutput';}\n        my\
 $selResultType;\n        foreach my $resultType (\
@$resultTypes) {\n            if($resultType->{typ\
e} eq $outformat) {\n                $selResultTyp\
e = $resultType;\n            }\n        }\n      \
  $res=$soap->poll($jobid, $selResultType->{type})\
;\n	if($outfile eq '-') {\n	     write_file($outfi\
le, $res);\n	} else {\n	    write_file($outfile.'.\
'.$selResultType->{ext}, $res);\n	}\n    } else { \
# Data types available\n        # Write a file for\
 each output type\n        for my $resultType (@$r\
esultTypes){\n            #print STDERR \"Getting \
$resultType->{type}\\n\";\n            $res=$soap-\
>poll($jobid, $resultType->{type});\n	    if($outf\
ile eq '-') {\n		write_file($outfile, $res);\n	   \
 } else {\n		write_file($outfile.'.'.$resultType->\
{ext}, $res);\n	    }\n        }\n    }\n}\n\nsub \
read_file($) {\n    my $filename = shift;\n    my \
($content, $buffer);\n    if($filename eq '-') {\n\
	while(sysread(STDIN, $buffer, 1024)) {\n	    $con\
tent .= $buffer;\n	}\n    }\n    else { # File\n	o\
pen(FILE, $filename) or die \"Error: unable to ope\
n input file\";\n	while(sysread(FILE, $buffer, 102\
4)) {\n	    $content .= $buffer;\n	}\n	close(FILE)\
;\n    }\n    return $content;\n}\n\nsub write_fil\
e($$) {\n    my ($filename, $data) = @_;\n    prin\
t STDERR 'Creating result file: ' . $filename . \"\
\\n\";\n    if($filename eq '-') {\n	print STDOUT \
$data;\n    }\n    else {\n	open(FILE, \">$filenam\
e\") or die \"Error: unable to open output file\";\
\n	syswrite(FILE, $data);\n	close(FILE);\n    }\n}\
\n\nsub usage {\n    print STDERR <<EOF\nWU-BLAST\\
n========\n\nRapid sequence database search progra\
ms utilizing the BLAST algorithm.\n   \n[Required]\
\n\n      --email       : str  : user email addres\
s \n  -p, --program	    : str  : BLAST program to \
use: blastn, blastp, blastx, \n                   \
          tblastn or tblastx\n  -D, --database    \
: str  : database to search\n  seqFile           :\
 file : query sequence data file (\"-\" for STDIN)\
\n\n[Optional]\n\n  -m, --matrix	    : str  : scor\
ing matrix\n  -E, --exp	    : real : 0<E<= 1000. S\
tatistical significance threshold\n               \
              for reporting database sequence matc\
hes.\n  -e, --echofilter  :      : display the fil\
tered query sequence in the output\n  -f, --filter\
	    : str  : activates filtering of the query seq\
uence\n  -b, --alignments  : int  : number of alig\
nments to be reported\n  -s, --scores	    : int  :\
 number of scores to be reported\n  -S, --sensitiv\
ity : str  :\n  -t, --sort	    : str  :\n  -T, --s\
tats       : str  :\n  -d, --strand      : str  : \
DNA strand to search with in DNA vs. DNA searches \
\n  -c, --topcombon   :      :\n\n[General]	\n\n  \
-h, --help       :      : prints this help text\n \
 -a, --async      :      : forces to make an async\
hronous query\n      --status     :      : poll fo\
r the status of a job\n      --polljob    :      :\
 poll for the results of a job\n  -j, --jobid     \
 : str  : jobid that was returned when an asynchro\
nous job \n                            was submitt\
ed.\n  -O, --outfile    : str  : name of the file \
results should be written to \n                   \
         (default is based on the jobid; \"-\" for\
 STDOUT)\n  -o, --outformat  : str  : txt or xml o\
utput (no file is written)\n      --trace	   :    \
  : show SOAP messages being interchanged \n\nSync\
hronous job:\n\n  The results/errors are returned \
as soon as the job is finished.\n  Usage: $scriptN\
ame --email <your\\@email> [options...] seqFile\n \
 Returns: saves the results to disk\n\nAsynchronou\
s job:\n\n  Use this if you want to retrieve the r\
esults at a later time. The results \n  are stored\
 for up to 24 hours. \n  The asynchronous submissi\
on mode is recommended when users are submitting \\
n  batch jobs or large database searches	\n  Usage\
: $scriptName --async --email <your\\@email> [opti\
ons...] seqFile\n  Returns : jobid\n\n  Use the jo\
bid to query for the status of the job. \n  Usage:\
 $scriptName --status --jobid <jobId>\n  Returns :\
 string indicating the status of the job:\n    DON\
E - job has finished\n    RUNNING - job is running\
\n    NOT_FOUND - job cannot be found\n    ERROR -\
 the jobs has encountered an error\n\n  When done,\
 use the jobid to retrieve the status of the job. \
\n  Usage: $scriptName --polljob --jobid <jobId> [\
--outfile string]\n  Returns: saves the results to\
 disk\n\n[Help]\n\nFor more detailed help informat\
ion refer to \nhttp://www.ebi.ac.uk/blast2/WU-Blas\
t2_Help_frame.html\n \nEOF\n;\n}\n","\nmy $WSDL = \
'http://www.ebi.ac.uk/Tools/webservices/wsdl/WSBla\
stpgp.wsdl';\n\nuse SOAP::Lite;\nuse Getopt::Long \
qw(:config no_ignore_case bundling);\nuse File::Ba\
sename;\n\nmy $checkInterval = 15;\n\nmy %params=(\
\n	    'async' => '1', # Use async mode and simula\
te sync mode in client\n	    );\nGetOptions(\n    \
\"mode=s\"           => \\$params{mode}, # Search \
mode: PSI-Blast or PHI-Blast\n    \"database|d=s\"\
     => \\$params{database}, # Database to search\\
n    \"matrix|M=s\"       => \\$params{matrix},# S\
coring maxtrix\n    \"exp|e=f\"          => \\$par\
ams{exp}, # E-value\n    \"expmulti|h=f\"     => \\
\$params{expmulti}, # E-value\n    \"filter|F=s\" \
      => \\$params{filter}, # Low complexity filte\
r\n    \"dropoff|X=i\"      => \\$params{dropoff},\
 # Dropoff score\n    \"finaldropoff|Z=i\" => \\$p\
arams{finaldropoff}, # Final dropoff score\n    \"\
scores|v=i\"       => \\$params{scores}, # Max num\
ber of scores\n    \"align=i\"          => \\$para\
ms{align}, # Alignment view\n    \"startregion|S=i\
\"  => \\$params{startregion}, # Start of region i\
n query\n    \"endregion|H=i\"    => \\$params{end\
region}, # End of region in query\n    \"maxpasses\
|j=i\"    => \\$params{maxpasses}, # Number of PSI\
 iterations\n    \"opengap|G=i\"      => \\$params\
{opengap}, # Gap open penalty\n    \"extendgap|E=i\
\"    => \\$params{extendgap}, # Gap extension pen\
alty\n    \"pattern=s\"        => \\$params{patter\
n}, # PHI-BLAST pattern\n    \"usagemode|p=s\"    \
=> \\$params{usagemode}, # PHI-BLAST program\n    \
\"appxml=s\"         => \\$params{appxml}, # Appli\
cation XML\n    \"sequence=s\"       => \\$sequenc\
e, # Query sequence\n    \"help\"	       => \\$hel\
p, # Usage info\n    \"polljob\"	       => \\$poll\
job, # Get results\n    \"status\"	       => \\$st\
atus, # Get status\n    \"ids\"      	       => \\\
$ids, # Get ids from result\n    \"jobid=s\"      \
    => \\$jobid, # JobId\n    \"outfile=s\"       \
 => \\$outfile, # Output filename\n    \"outformat\
|o=s\"    => \\$outformat, # Output file format\n \
   \"async|a\"	       => \\$async, # Async submiss\
ion\n    \"email=s\"          => \\$params{email},\
 # User e-mail address\n    \"trace\"            =\
> \\$trace, # Show SOAP messages\n    );\n\nmy $sc\
riptName = basename($0, ());\nif($help) {\n    &us\
age();\n    exit(0);\n}\n\nif ($trace){\n    print\
 \"Tracing active\\n\";\n    SOAP::Lite->import(+t\
race => 'debug');\n}\n\nmy $soap = SOAP::Lite\n   \
 ->service($WSDL)\n    ->on_fault(sub {\n        m\
y $soap = shift;\n        my $res = shift;\n      \
  # Throw an exception for all faults\n        if(\
ref($res) eq '') {\n            die($res);\n      \
  } else {\n            die($res->faultstring);\n \
       }\n        return new SOAP::SOM;\n    }\n  \
             );\n\nif( !($polljob || $status || $i\
ds) &&\n    !( (defined($ARGV[0]) && -f $ARGV[0]) \
|| defined($sequence) )\n    ) {\n    print STDERR\
 'Error: bad option combination', \"\\n\";\n    &u\
sage();\n    exit(1);\n}\nelsif($polljob && define\
d($jobid)) {\n    print \"Getting results for job \
$jobid\\n\";\n    getResults($jobid);\n}\nelsif($s\
tatus && defined($jobid)) {\n    print STDERR \"Ge\
tting status for job $jobid\\n\";\n    my $result \
= $soap->checkStatus($jobid);\n    print STDOUT $r\
esult, \"\\n\";\n    if($result eq 'DONE') {\n	pri\
nt STDERR \"To get results: $scriptName --polljob \
--jobid $jobid\\n\";\n    }\n}  \nelsif($ids && de\
fined($jobid)) {\n    print STDERR \"Getting ids f\
rom job $jobid\\n\";\n    getIds($jobid);\n}\nelse\
 {\n    if(-f $ARGV[0]) {	\n	$content={type=>'sequ\
ence', content=>read_file($ARGV[0])};	\n    }\n   \
 if($sequence) {	\n	if(-f $sequence) {\n	    $cont\
ent={type=>'sequence', content=>read_file($sequenc\
e)};	\n	} else {\n	    $content={type=>'sequence',\
 content=>$sequence};\n	}\n    }\n    push @conten\
t, $content;\n\n    my $jobid;\n    my $paramsData\
 = SOAP::Data->name('params')->type(map=>\\%params\
);\n    my $contentData = SOAP::Data->name('conten\
t')->value(\\@content);\n    # For SOAP::Lite 0.60\
 and earlier parameters are passed directly\n    i\
f($SOAP::Lite::VERSION eq '0.60' || $SOAP::Lite::V\
ERSION =~ /0\\.[1-5]/) {\n        $jobid = $soap->\
runBlastpgp($paramsData, $contentData);\n    }\n  \
  # For SOAP::Lite 0.69 and later parameter handli\
ng is different, so pass\n    # undef's for templa\
ted params, and then pass the formatted args.\n   \
 else {\n        $jobid = $soap->runBlastpgp(undef\
, undef,\n				    $paramsData, $contentData);\n   \
 }\n\n    if (defined($async)) {\n	print STDOUT $j\
obid, \"\\n\";\n        print STDERR \"To check st\
atus: $scriptName --status --jobid $jobid\\n\";\n \
   } else { # Synchronous mode\n        print STDE\
RR \"JobId: $jobid\\n\";\n        sleep 1;\n      \
  getResults($jobid);\n    }\n}\n\nsub getIds($) {\
\n    $jobid = shift;\n    my $results = $soap->ge\
tIds($jobid);\n    for $result (@$results){\n	prin\
t \"$result\\n\";\n    }\n}\n\nsub clientPoll($) {\
\n    my $jobid = shift;\n    my $result = 'PENDIN\
G';\n    # Check status and wait if not finished\n\
    #print STDERR \"Checking status: $jobid\\n\";\\
n    while($result eq 'RUNNING' || $result eq 'PEN\
DING') {\n        $result = $soap->checkStatus($jo\
bid);\n        print STDERR \"$result\\n\";\n     \
   if($result eq 'RUNNING' || $result eq 'PENDING'\
) {\n            # Wait before polling again.\n   \
         sleep $checkInterval;\n        }\n    }\n\
}\n\nsub getResults($) {\n    $jobid = shift;\n   \
 # Check status, and wait if not finished\n    cli\
entPoll($jobid);\n    # Use JobId if output file n\
ame is not defined\n    unless(defined($outfile)) \
{\n        $outfile=$jobid;\n    }\n    # Get list\
 of data types\n    my $resultTypes = $soap->getRe\
sults($jobid);\n    # Get the data and write it to\
 a file\n    if(defined($outformat)) { # Specified\
 data type\n        my $selResultType;\n        fo\
reach my $resultType (@$resultTypes) {\n          \
  if($resultType->{type} eq $outformat) {\n       \
         $selResultType = $resultType;\n          \
  }\n        }\n        $res=$soap->poll($jobid, $\
selResultType->{type});\n        write_file($outfi\
le.'.'.$selResultType->{ext}, $res);\n    } else {\
 # Data types available\n        # Write a file fo\
r each output type\n        for my $resultType (@$\
resultTypes){\n            #print \"Getting $resul\
tType->{type}\\n\";\n            $res=$soap->poll(\
$jobid, $resultType->{type});\n            write_f\
ile($outfile.'.'.$resultType->{ext}, $res);\n     \
   }\n    }\n}\n\nsub read_file($) {\n    my $file\
name = shift;\n    open(FILE, $filename);\n    my \
$content;\n    my $buffer;\n    while(sysread(FILE\
, $buffer, 1024)) {\n	$content.= $buffer;\n    }\n\
    close(FILE);  \n    return $content;\n}\n\nsub\
 write_file($$) {\n    my ($tmp,$entity) = @_;\n  \
  print STDERR \"Creating result file: \".$tmp.\"\\
\n\";\n    unless(open (FILE, \">$tmp\")) {\n	retu\
rn 0;\n    }\n    syswrite(FILE, $entity);\n    cl\
ose (FILE);\n    return 1;\n}\n\nsub usage {\n    \
print STDERR <<EOF\nBlastpgp\n========\n   \nThe b\
lastpgp program implements the PSI-BLAST and PHI-B\
LAST variations\nof NCBI BLAST.\n\nFor more detail\
ed help information refer to\nhttp://www.ebi.ac.uk\
/blastpgp/blastpsi_help_frame.html\n \nBlastpgp sp\
ecific options:\n\n[Required]\n\n      --mode     \
       : str  : search mode to use: PSI-Blast or P\
HI-Blast\n  -d, --database        : str  : protein\
 database to search\n  seqFile               : fil\
e : query sequence\n\n[Optional]\n\n  -M, --matrix\
          : str  : scoring matrix\n  -e, --exp    \
         : real : Expectation value\n  -h, --expmu\
lti        : real : threshold (multipass model)\n \
 -F, --filter          : str  : filter query seque\
nce with SEG [T,F]\n  -m, --align           : int \
 : alignment view option:\n                       \
          0 - pairwise, 1 - M/S identities,\n     \
                            2 - M/S non-identities\
, 3 - Flat identities,\n                          \
       4 - Flat non-identities\n  -G, --opengap   \
      : int  : cost to open a gap\n  -E, --extendg\
ap       : int  : cost to extend a gap\n  -g, --ga\
palign        : str  : Gapped [T,F]\n  -v, --score\
s          : int  : number of scores to be reporte\
d\n  -j, --maxpasses       : int  : number of iter\
ations\n  -X, --dropoff         : int  : Dropoff s\
core\n  -Z, --finaldropoff    : int  : Dropoff for\
 final alignment\n  -S, --startregion     : int  :\
 Start of required region in query\n  -H, --endreg\
ion       : int  : End of required region in query\
\n  -k, --pattern         : str  : Hit File (PHI-B\
LAST only)\n  -p, --usagemode       : str  : Progr\
am option (PHI-BLAST only):\n                     \
            blastpgp, patseedp, seedp\n\n[General]\
\n\n      --help            :      : prints this h\
elp text\n  -a, --async           :      : forces \
to make an asynchronous query\n      --status     \
     :      : poll for the status of a job\n      \
--polljob         :      : poll for the results of\
 a job\n      --jobid           : str  : jobid of \
an asynchronous job\n      --ids             :    \
  : get hit identifiers for result \n  -O, --outfi\
le         : str  : name of the file results shoul\
d be written to\n                                 \
(default is based on the jobid)\n  -o, --outformat\
       : str  : txt or xml output (no file is writ\
ten)\n      --trace           :      : show SOAP m\
essages being interchanged\n\nSynchronous job:\n\n\
  The results/errors are returned as soon as the j\
ob is finished.\n  Usage: blastpgp.pl --email <you\
r@email> [options...] seqfile\n  Returns: saves th\
e results to disk\n\nAsynchronous job:\n\n  Use th\
is if you want to retrieve the results at a later \
time. The results\n  are stored for up to 24 hours\
.\n  The asynchronous submission mode is recommend\
ed when users are submitting\n  batch jobs or larg\
e database searches\n  Usage: blastpgp.pl --email \
<your@email> --async [options...] seqFile\n  Retur\
ns: jobid\n\n  Use the jobid to query for the stat\
us of the job.\n  Usage: blastpgp.pl --status --jo\
bid <jobId>\n  Returns: string indicating the stat\
us of the job\n    DONE - job has finished\n    RU\
NNING - job is running\n    NOT_FOUND - job cannot\
 be found\n    ERROR - the jobs has encountered an\
 error\n\n  When done, use the jobid to retrieve t\
he results of the job.\n  Usage: blastpgp.pl --pol\
ljob --jobid <jobId> [--outfile <fileName>]\n  Ret\
urns: saves the results to disk\nEOF\n;\n}\n","\n=\
head1 NAME\n\nncbiblast_lwp.pl\n\n=head1 DESCRIPTI\
ON\n\nNCBI BLAST REST web service Perl client usin\
g L<LWP>.\n\nTested with:\n\n=over\n\n=item *\nL<L\
WP> 5.79, L<XML::Simple> 2.12 and Perl 5.8.3\n\n=i\
tem *\nL<LWP> 5.805, L<XML::Simple> 2.14 and Perl \
5.8.7\n\n=item *\nL<LWP> 5.820, L<XML::Simple> 2.1\
8 and Perl 5.10.0 (Ubuntu 9.04)\n\n=back\n\nFor fu\
rther information see:\n\n=over\n\n=item *\nL<http\
://www.ebi.ac.uk/Tools/webservices/services/sss/nc\
bi_blast_rest>\n\n=item *\nL<http://www.ebi.ac.uk/\
Tools/webservices/tutorials/perl>\n\n=back\n\n=hea\
d1 VERSION\n\n$Id: ncbiblast_lwp.pl 1317 2009-09-0\
3 15:44:11Z hpm $\n\n=cut\n\nuse strict;\nuse warn\
ings;\n\nuse English;\nuse LWP;\nuse XML::Simple;\\
nuse Getopt::Long qw(:config no_ignore_case bundli\
ng);\nuse File::Basename;\nuse Data::Dumper;\n\nmy\
 $baseUrl = 'http://www.ebi.ac.uk/Tools/services/r\
est/ncbiblast';\n\nmy $checkInterval = 3;\n\nmy $o\
utputLevel = 1;\n\nmy $numOpts = scalar(@ARGV);\nm\
y %params = ( 'debugLevel' => 0 );\n\nmy %tool_par\
ams = ();\nGetOptions(\n\n	# Tool specific options\
\n	'program|p=s'  => \\$tool_params{'program'},   \
# blastp, blastn, blastx, etc.\n	'database|D=s' =>\
 \\$params{'database'},       # Database(s) to sea\
rch\n	'matrix|m=s'   => \\$tool_params{'matrix'}, \
   # Scoring martix to use\n	'exp|E=f'      => \\$\
tool_params{'exp'},       # E-value threshold\n	'f\
ilter|f=s'   => \\$tool_params{'filter'},    # Low\
 complexity filter\n	'align|A=i'    => \\$tool_par\
ams{'align'},     # Pairwise alignment format\n	's\
cores|s=i'   => \\$tool_params{'scores'},    # Num\
ber of scores\n	'alignments|n=i' => \\$tool_params\
{'alignments'},   # Number of alignments\n	'dropof\
f|d=i'    => \\$tool_params{'dropoff'},      # Dro\
poff score\n	'match_scores=s' => \\$tool_params{'m\
atch_scores'}, # Match/missmatch scores\n	'match|u\
=i'      => \\$params{'match'},             # Matc\
h score\n	'mismatch|v=i'   => \\$params{'mismatch'\
},          # Mismatch score\n	'gapopen|o=i'    =>\
 \\$tool_params{'gapopen'},      # Open gap penalt\
y\n	'gapext|x=i'     => \\$tool_params{'gapext'}, \
      # Gap extension penality\n	'gapalign|g'     \
=> \\$tool_params{'gapalign'},     # Optimise gap \
alignments\n	'stype=s' => \\$tool_params{'stype'},\
    # Sequence type\n	'seqrange=s' => \\$tool_para\
ms{'seqrange'},    # Query subsequence\n	'sequence\
=s' => \\$params{'sequence'},         # Query sequ\
ence\n	'multifasta' => \\$params{'multifasta'},   \
    # Multiple fasta input\n\n	# Compatability opt\
ions, old command-line\n	'numal|n=i'     => \\$par\
ams{'numal'},        # Number of alignments\n	'ope\
ngap|o=i'   => \\$params{'opengap'},      # Open g\
ap penalty\n	'extendgap|x=i' => \\$params{'extendg\
ap'},    # Gap extension penality\n	\n	# Generic o\
ptions\n	'email=s'       => \\$params{'email'},   \
       # User e-mail address\n	'title=s'       => \
\\$params{'title'},          # Job title\n	'outfil\
e=s'     => \\$params{'outfile'},        # Output \
file name\n	'outformat=s'   => \\$params{'outforma\
t'},      # Output file type\n	'jobid=s'       => \
\\$params{'jobid'},          # JobId\n	'help|h'   \
     => \\$params{'help'},           # Usage help\\
n	'async'         => \\$params{'async'},          \
# Asynchronous submission\n	'polljob'       => \\$\
params{'polljob'},        # Get results\n	'resultT\
ypes'   => \\$params{'resultTypes'},    # Get resu\
lt types\n	'status'        => \\$params{'status'},\
         # Get status\n	'params'        => \\$para\
ms{'params'},         # List input parameters\n	'p\
aramDetail=s' => \\$params{'paramDetail'},    # Ge\
t details for parameter\n	'quiet'         => \\$pa\
rams{'quiet'},          # Decrease output level\n	\
'verbose'       => \\$params{'verbose'},        # \
Increase output level\n	'debugLevel=i'  => \\$para\
ms{'debugLevel'},     # Debug output level\n	'base\
Url=s'     => \\$baseUrl,                  # Base \
URL for service.\n);\nif ( $params{'verbose'} ) { \
$outputLevel++ }\nif ( $params{'$quiet'} )  { $out\
putLevel-- }\n\n&print_debug_message( 'MAIN', 'LWP\
::VERSION: ' . $LWP::VERSION,\n	1 );\n\n&print_deb\
ug_message( 'MAIN', \"params:\\n\" . Dumper( \\%pa\
rams ),           11 );\n&print_debug_message( 'MA\
IN', \"tool_params:\\n\" . Dumper( \\%tool_params \
), 11 );\n\nmy $scriptName = basename( $0, () );\n\
\nif ( $params{'help'} || $numOpts == 0 ) {\n	&usa\
ge();\n	exit(0);\n}\n\n&print_debug_message( 'MAIN\
', 'baseUrl: ' . $baseUrl, 1 );\n\nif (\n	!(\n		  \
 $params{'polljob'}\n		|| $params{'resultTypes'}\n\
		|| $params{'status'}\n		|| $params{'params'}\n		\
|| $params{'paramDetail'}\n	)\n	&& !( defined( $AR\
GV[0] ) || defined( $params{'sequence'} ) )\n  )\n\
{\n\n	# Bad argument combination, so print error m\
essage and usage\n	print STDERR 'Error: bad option\
 combination', \"\\n\";\n	&usage();\n	exit(1);\n}\\
n\nelsif ( $params{'params'} ) {\n	&print_tool_par\
ams();\n}\n\nelsif ( $params{'paramDetail'} ) {\n	\
&print_param_details( $params{'paramDetail'} );\n}\
\n\nelsif ( $params{'status'} && defined( $params{\
'jobid'} ) ) {\n	&print_job_status( $params{'jobid\
'} );\n}\n\nelsif ( $params{'resultTypes'} && defi\
ned( $params{'jobid'} ) ) {\n	&print_result_types(\
 $params{'jobid'} );\n}\n\nelsif ( $params{'polljo\
b'} && defined( $params{'jobid'} ) ) {\n	&get_resu\
lts( $params{'jobid'} );\n}\n\nelse {\n\n	# Multip\
le input sequence mode, assume fasta format.\n	if \
( $params{'multifasta'} ) {\n		&multi_submit_job()\
;\n	}\n\n	# Entry identifier list file.\n	elsif ((\
 defined( $params{'sequence'} ) && $params{'sequen\
ce'} =~ m/^\\@/ )\n		|| ( defined( $ARGV[0] ) && $\
ARGV[0] =~ m/^\\@/ ) )\n	{\n		my $list_filename = \
$params{'sequence'} || $ARGV[0];\n		$list_filename\
 =~ s/^\\@//;\n		&list_file_submit_job($list_filen\
ame);\n	}\n\n	# Default: single sequence/identifie\
r.\n	else {\n\n		# Load the sequence data and subm\
it.\n		&submit_job( &load_data() );\n	}\n}\n\n=hea\
d1 FUNCTIONS\n\n=cut\n\n\n=head2 rest_request()\n\\
nPerform a REST request.\n\n  my $response_str = &\
rest_request($url);\n\n=cut\n\nsub rest_request {\\
n	print_debug_message( 'rest_request', 'Begin', 11\
 );\n	my $requestUrl = shift;\n	print_debug_messag\
e( 'rest_request', 'URL: ' . $requestUrl, 11 );\n\\
n	# Create a user agent\n	my $ua = LWP::UserAgent-\
>new();\n	'$Revision: 1317 $' =~ m/(\\d+)/;\n	$ua-\
>agent(\"EBI-Sample-Client/$1 ($scriptName; $OSNAM\
E) \" . $ua->agent());\n	$ua->env_proxy;\n\n	# Per\
form the request\n	my $response = $ua->get($reques\
tUrl);\n	print_debug_message( 'rest_request', 'HTT\
P status: ' . $response->code,\n		11 );\n\n	# Chec\
k for HTTP error codes\n	if ( $response->is_error \
) {\n		$response->content() =~ m/<h1>([^<]+)<\\/h1\
>/;\n		die 'http status: ' . $response->code . ' '\
 . $response->message . '  ' . $1;\n	}\n	print_deb\
ug_message( 'rest_request', 'End', 11 );\n\n	# Ret\
urn the response data\n	return $response->content(\
);\n}\n\n=head2 rest_get_parameters()\n\nGet list \
of tool parameter names.\n\n  my (@param_list) = &\
rest_get_parameters();\n\n=cut\n\nsub rest_get_par\
ameters {\n	print_debug_message( 'rest_get_paramet\
ers', 'Begin', 1 );\n	my $url                = $ba\
seUrl . '/parameters/';\n	my $param_list_xml_str =\
 rest_request($url);\n	my $param_list_xml     = XM\
Lin($param_list_xml_str);\n	my (@param_list)      \
 = @{ $param_list_xml->{'id'} };\n	print_debug_mes\
sage( 'rest_get_parameters', 'End', 1 );\n	return \
(@param_list);\n}\n\n=head2 rest_get_parameter_det\
ails()\n\nGet details of a tool parameter.\n\n  my\
 $paramDetail = &rest_get_parameter_details($param\
_name);\n\n=cut\n\nsub rest_get_parameter_details \
{\n	print_debug_message( 'rest_get_parameter_detai\
ls', 'Begin', 1 );\n	my $parameterId = shift;\n	pr\
int_debug_message( 'rest_get_parameter_details',\n\
		'parameterId: ' . $parameterId, 1 );\n	my $url  \
                = $baseUrl . '/parameterdetails/' \
. $parameterId;\n	my $param_detail_xml_str = rest_\
request($url);\n	my $param_detail_xml     = XMLin(\
$param_detail_xml_str);\n	print_debug_message( 're\
st_get_parameter_details', 'End', 1 );\n	return ($\
param_detail_xml);\n}\n\n=head2 rest_run()\n\nSubm\
it a job.\n\n  my $job_id = &rest_run($email, $tit\
le, \\%params );\n\n=cut\n\nsub rest_run {\n	print\
_debug_message( 'rest_run', 'Begin', 1 );\n	my $em\
ail  = shift;\n	my $title  = shift;\n	my $params =\
 shift;\n	print_debug_message( 'rest_run', 'email:\
 ' . $email, 1 );\n	if ( defined($title) ) {\n		pr\
int_debug_message( 'rest_run', 'title: ' . $title,\
 1 );\n	}\n	print_debug_message( 'rest_run', 'para\
ms: ' . Dumper($params), 1 );\n\n	# User agent to \
perform http requests\n	my $ua = LWP::UserAgent->n\
ew();\n	$ua->env_proxy;\n\n	# Clean up parameters\\
n	my (%tmp_params) = %{$params};\n	$tmp_params{'em\
ail'} = $email;\n	$tmp_params{'title'} = $title;\n\
	foreach my $param_name ( keys(%tmp_params) ) {\n	\
	if ( !defined( $tmp_params{$param_name} ) ) {\n		\
	delete $tmp_params{$param_name};\n		}\n	}\n\n	# S\
ubmit the job as a POST\n	my $url = $baseUrl . '/r\
un';\n	my $response = $ua->post( $url, \\%tmp_para\
ms );\n	print_debug_message( 'rest_run', 'HTTP sta\
tus: ' . $response->code, 11 );\n	print_debug_mess\
age( 'rest_run',\n		'request: ' . $response->reque\
st()->content(), 11 );\n\n	# Check for HTTP error \
codes\n	if ( $response->is_error ) {\n		$response-\
>content() =~ m/<h1>([^<]+)<\\/h1>/;\n		die 'http \
status: ' . $response->code . ' ' . $response->mes\
sage . '  ' . $1;\n	}\n\n	# The job id is returned\
\n	my $job_id = $response->content();\n	print_debu\
g_message( 'rest_run', 'End', 1 );\n	return $job_i\
d;\n}\n\n=head2 rest_get_status()\n\nCheck the sta\
tus of a job.\n\n  my $status = &rest_get_status($\
job_id);\n\n=cut\n\nsub rest_get_status {\n	print_\
debug_message( 'rest_get_status', 'Begin', 1 );\n	\
my $job_id = shift;\n	print_debug_message( 'rest_g\
et_status', 'jobid: ' . $job_id, 2 );\n	my $status\
_str = 'UNKNOWN';\n	my $url        = $baseUrl . '/\
status/' . $job_id;\n	$status_str = &rest_request(\
$url);\n	print_debug_message( 'rest_get_status', '\
status_str: ' . $status_str, 2 );\n	print_debug_me\
ssage( 'rest_get_status', 'End', 1 );\n	return $st\
atus_str;\n}\n\n=head2 rest_get_result_types()\n\n\
Get list of result types for finished job.\n\n  my\
 (@result_types) = &rest_get_result_types($job_id)\
;\n\n=cut\n\nsub rest_get_result_types {\n	print_d\
ebug_message( 'rest_get_result_types', 'Begin', 1 \
);\n	my $job_id = shift;\n	print_debug_message( 'r\
est_get_result_types', 'jobid: ' . $job_id, 2 );\n\
	my (@resultTypes);\n	my $url                     \
 = $baseUrl . '/resulttypes/' . $job_id;\n	my $res\
ult_type_list_xml_str = &rest_request($url);\n	my \
$result_type_list_xml     = XMLin($result_type_lis\
t_xml_str);\n	(@resultTypes) = @{ $result_type_lis\
t_xml->{'type'} };\n	print_debug_message( 'rest_ge\
t_result_types',\n		scalar(@resultTypes) . ' resul\
t types', 2 );\n	print_debug_message( 'rest_get_re\
sult_types', 'End', 1 );\n	return (@resultTypes);\\
n}\n\n=head2 rest_get_result()\n\nGet result data \
of a specified type for a finished job.\n\n  my $r\
esult = rest_get_result($job_id, $result_type);\n\\
n=cut\n\nsub rest_get_result {\n	print_debug_messa\
ge( 'rest_get_result', 'Begin', 1 );\n	my $job_id \
= shift;\n	my $type   = shift;\n	print_debug_messa\
ge( 'rest_get_result', 'jobid: ' . $job_id, 1 );\n\
	print_debug_message( 'rest_get_result', 'type: ' \
. $type,    1 );\n	my $url    = $baseUrl . '/resul\
t/' . $job_id . '/' . $type;\n	my $result = &rest_\
request($url);\n	print_debug_message( 'rest_get_re\
sult', length($result) . ' characters',\n		1 );\n	\
print_debug_message( 'rest_get_result', 'End', 1 )\
;\n	return $result;\n}\n\n\n=head2 print_debug_mes\
sage()\n\nPrint debug message at specified debug l\
evel.\n\n  &print_debug_message($method_name, $mes\
sage, $level);\n\n=cut\n\nsub print_debug_message \
{\n	my $function_name = shift;\n	my $message      \
 = shift;\n	my $level         = shift;\n	if ( $lev\
el <= $params{'debugLevel'} ) {\n		print STDERR '[\
', $function_name, '()] ', $message, \"\\n\";\n	}\\
n}\n\n=head2 print_tool_params()\n\nPrint list of \
tool parameters.\n\n  &print_tool_params();\n\n=cu\
t\n\nsub print_tool_params {\n	print_debug_message\
( 'print_tool_params', 'Begin', 1 );\n	my (@param_\
list) = &rest_get_parameters();\n	foreach my $para\
m ( sort(@param_list) ) {\n		print $param, \"\\n\"\
;\n	}\n	print_debug_message( 'print_tool_params', \
'End', 1 );\n}\n\n=head2 print_param_details()\n\n\
Print details of a tool parameter.\n\n  &print_par\
am_details($param_name);\n\n=cut\n\nsub print_para\
m_details {\n	print_debug_message( 'print_param_de\
tails', 'Begin', 1 );\n	my $paramName = shift;\n	p\
rint_debug_message( 'print_param_details', 'paramN\
ame: ' . $paramName, 2 );\n	my $paramDetail = &res\
t_get_parameter_details($paramName);\n	print $para\
mDetail->{'name'}, \"\\t\", $paramDetail->{'type'}\
, \"\\n\";\n	print $paramDetail->{'description'}, \
\"\\n\";\n	foreach my $value ( @{ $paramDetail->{'\
values'}->{'value'} } ) {\n		print $value->{'value\
'};\n		if ( $value->{'defaultValue'} eq 'true' ) {\
\n			print \"\\t\", 'default';\n		}\n		print \"\\n\
\";\n		print \"\\t\", $value->{'label'}, \"\\n\";\\
n	}\n	print_debug_message( 'print_param_details', \
'End', 1 );\n}\n\n=head2 print_job_status()\n\nPri\
nt status of a job.\n\n  &print_job_status($job_id\
);\n\n=cut\n\nsub print_job_status {\n	print_debug\
_message( 'print_job_status', 'Begin', 1 );\n	my $\
jobid = shift;\n	print_debug_message( 'print_job_s\
tatus', 'jobid: ' . $jobid, 1 );\n	if ( $outputLev\
el > 0 ) {\n		print STDERR 'Getting status for job\
 ', $jobid, \"\\n\";\n	}\n	my $result = &rest_get_\
status($jobid);\n	print \"$result\\n\";\n	if ( $re\
sult eq 'FINISHED' && $outputLevel > 0 ) {\n		prin\
t STDERR \"To get results: $scriptName --polljob -\
-jobid \" . $jobid\n		  . \"\\n\";\n	}\n	print_deb\
ug_message( 'print_job_status', 'End', 1 );\n}\n\n\
=head2 print_result_types()\n\nPrint available res\
ult types for a job.\n\n  &print_result_types($job\
_id);\n\n=cut\n\nsub print_result_types {\n	print_\
debug_message( 'result_types', 'Begin', 1 );\n	my \
$jobid = shift;\n	print_debug_message( 'result_typ\
es', 'jobid: ' . $jobid, 1 );\n	if ( $outputLevel \
> 0 ) {\n		print STDERR 'Getting result types for \
job ', $jobid, \"\\n\";\n	}\n	my $status = &rest_g\
et_status($jobid);\n	if ( $status eq 'PENDING' || \
$status eq 'RUNNING' ) {\n		print STDERR 'Error: J\
ob status is ', $status,\n		  '. To get result typ\
es the job must be finished.', \"\\n\";\n	}\n	else\
 {\n		my (@resultTypes) = &rest_get_result_types($\
jobid);\n		if ( $outputLevel > 0 ) {\n			print STD\
OUT 'Available result types:', \"\\n\";\n		}\n		fo\
reach my $resultType (@resultTypes) {\n			print ST\
DOUT $resultType->{'identifier'}, \"\\n\";\n			if \
( defined( $resultType->{'label'} ) ) {\n				print\
 STDOUT \"\\t\", $resultType->{'label'}, \"\\n\";\\
n			}\n			if ( defined( $resultType->{'description\
'} ) ) {\n				print STDOUT \"\\t\", $resultType->{\
'description'}, \"\\n\";\n			}\n			if ( defined( $\
resultType->{'mediaType'} ) ) {\n				print STDOUT \
\"\\t\", $resultType->{'mediaType'}, \"\\n\";\n			\
}\n			if ( defined( $resultType->{'fileSuffix'} ) \
) {\n				print STDOUT \"\\t\", $resultType->{'file\
Suffix'}, \"\\n\";\n			}\n		}\n		if ( $status eq '\
FINISHED' && $outputLevel > 0 ) {\n			print STDERR\
 \"\\n\", 'To get results:', \"\\n\",\n			  \"  $s\
criptName --polljob --jobid \" . $params{'jobid'} \
. \"\\n\",\n			  \"  $scriptName --polljob --outfo\
rmat <type> --jobid \"\n			  . $params{'jobid'} . \
\"\\n\";\n		}\n	}\n	print_debug_message( 'result_t\
ypes', 'End', 1 );\n}\n\n=head2 submit_job()\n\nSu\
bmit a job to the service.\n\n  &submit_job($seq);\
\n\n=cut\n\nsub submit_job {\n	print_debug_message\
( 'submit_job', 'Begin', 1 );\n\n	# Set input sequ\
ence\n	$tool_params{'sequence'} = shift;\n\n	# Loa\
d parameters\n	&load_params();\n\n	# Submit the jo\
b\n	my $jobid = &rest_run( $params{'email'}, $para\
ms{'title'}, \\%tool_params );\n\n	# Simulate sync\
/async mode\n	if ( defined( $params{'async'} ) ) {\
\n		print STDOUT $jobid, \"\\n\";\n		if ( $outputL\
evel > 0 ) {\n			print STDERR\n			  \"To check sta\
tus: $scriptName --status --jobid $jobid\\n\";\n		\
}\n	}\n	else {\n		if ( $outputLevel > 0 ) {\n			pr\
int STDERR \"JobId: $jobid\\n\";\n		}\n		sleep 1;\\
n		&get_results($jobid);\n	}\n	print_debug_message\
( 'submit_job', 'End', 1 );\n}\n\n=head2 multi_sub\
mit_job()\n\nSubmit multiple jobs assuming input i\
s a collection of fasta formatted sequences.\n\n  \
&multi_submit_job();\n\n=cut\n\nsub multi_submit_j\
ob {\n	print_debug_message( 'multi_submit_job', 'B\
egin', 1 );\n	my $jobIdForFilename = 1;\n	$jobIdFo\
rFilename = 0 if ( defined( $params{'outfile'} ) )\
;\n	my (@filename_list) = ();\n\n	# Query sequence\
\n	if ( defined( $ARGV[0] ) ) {    # Bare option\n\
		if ( -f $ARGV[0] || $ARGV[0] eq '-' ) {    # Fil\
e\n			push( @filename_list, $ARGV[0] );\n		}\n	}\n\
	if ( $params{'sequence'} ) {                   # \
Via --sequence\n		if ( -f $params{'sequence'} || $\
params{'sequence'} eq '-' ) {    # File\n			push( \
@filename_list, $params{'sequence'} );\n		}\n	}\n\\
n	$/ = '>';\n	foreach my $filename (@filename_list\
) {\n		open( my $INFILE, '<', $filename )\n		  or \
die \"Error: unable to open file $filename ($!)\";\
\n		while (<$INFILE>) {\n			my $seq = $_;\n			$seq\
 =~ s/>$//;\n			if ( $seq =~ m/(\\S+)/ ) {\n				pr\
int STDERR \"Submitting job for: $1\\n\"\n				  if\
 ( $outputLevel > 0 );\n				$seq = '>' . $seq;\n		\
		&print_debug_message( 'multi_submit_job', $seq, \
11 );\n				&submit_job($seq);\n				$params{'outfil\
e'} = undef if ( $jobIdForFilename == 1 );\n			}\n\
		}\n		close $INFILE;\n	}\n	print_debug_message( '\
multi_submit_job', 'End', 1 );\n}\n\n=head2 list_f\
ile_submit_job()\n\nSubmit multiple jobs using a f\
ile containing a list of entry identifiers as \nin\
put.\n\n  &list_file_submit_job($list_filename)\n\\
n=cut\n\nsub list_file_submit_job {\n	my $filename\
         = shift;\n	my $jobIdForFilename = 1;\n	$j\
obIdForFilename = 0 if ( defined( $params{'outfile\
'} ) );\n\n	# Iterate over identifiers, submitting\
 each job\n	open( my $LISTFILE, '<', $filename )\n\
	  or die 'Error: unable to open file ' . $filenam\
e . ' (' . $! . ')';\n	while (<$LISTFILE>) {\n		my\
 $line = $_;\n		chomp($line);\n		if ( $line ne '' \
) {\n			&print_debug_message( 'list_file_submit_jo\
b', 'line: ' . $line, 2 );\n			if ( $line =~ m/\\w\
:\\w/ ) {    # Check this is an identifier\n				pr\
int STDERR \"Submitting job for: $line\\n\"\n				 \
 if ( $outputLevel > 0 );\n				&submit_job($line);\
\n			}\n			else {\n				print STDERR\n\"Warning: li\
ne \\\"$line\\\" is not recognised as an identifie\
r\\n\";\n			}\n		}\n		$params{'outfile'} = undef i\
f ( $jobIdForFilename == 1 );\n	}\n	close $LISTFIL\
E;\n}\n\n=head2 load_data()\n\nLoad sequence data \
from file or option specified on the command-line.\
\n\n  &load_data();\n\n=cut\n\nsub load_data {\n	p\
rint_debug_message( 'load_data', 'Begin', 1 );\n	m\
y $retSeq;\n\n	# Query sequence\n	if ( defined( $A\
RGV[0] ) ) {    # Bare option\n		if ( -f $ARGV[0] \
|| $ARGV[0] eq '-' ) {    # File\n			$retSeq = &re\
ad_file( $ARGV[0] );\n		}\n		else {               \
                      # DB:ID or sequence\n			$ret\
Seq = $ARGV[0];\n		}\n	}\n	if ( $params{'sequence'\
} ) {                   # Via --sequence\n		if ( -\
f $params{'sequence'} || $params{'sequence'} eq '-\
' ) {    # File\n			$retSeq = &read_file( $params{\
'sequence'} );\n		}\n		else {    # DB:ID or sequen\
ce\n			$retSeq = $params{'sequence'};\n		}\n	}\n	p\
rint_debug_message( 'load_data', 'End', 1 );\n	ret\
urn $retSeq;\n}\n\n=head2 load_params()\n\nLoad jo\
b parameters from command-line options.\n\n  &load\
_params();\n\n=cut\n\nsub load_params {\n	print_de\
bug_message( 'load_params', 'Begin', 1 );\n\n	# Da\
tabase(s) to search\n	my (@dbList) = split /[ ,]/,\
 $params{'database'};\n	$tool_params{'database'} =\
 \\@dbList;\n\n	# Match/missmatch\n	if ( $params{'\
match'} && $params{'missmatch'} ) {\n		$tool_param\
s{'match_scores'} =\n		  $params{'match'} . ',' . \
$params{'missmatch'};\n	}\n	\n	# Compatability opt\
ions, old command-line\n	if(!$tool_params{'alignme\
nts'} && $params{'numal'}) {\n		$tool_params{'alig\
nments'} = $params{'numal'};\n	}\n	if(!$tool_param\
s{'gapopen'} && $params{'opengap'}) {\n		$tool_par\
ams{'gapopen'} = $params{'opengap'};\n	}\n	if(!$to\
ol_params{'gapext'} && $params{'extendgap'}) {\n		\
$tool_params{'gapext'} = $params{'extendgap'};\n	}\
\n\n	print_debug_message( 'load_params', 'End', 1 \
);\n}\n\n=head2 client_poll()\n\nClient-side job p\
olling.\n\n  &client_poll($job_id);\n\n=cut\n\nsub\
 client_poll {\n	print_debug_message( 'client_poll\
', 'Begin', 1 );\n	my $jobid  = shift;\n	my $statu\
s = 'PENDING';\n\n	my $errorCount = 0;\n	while ($s\
tatus eq 'RUNNING'\n		|| $status eq 'PENDING'\n		|\
| ( $status eq 'ERROR' && $errorCount < 2 ) )\n	{\\
n		$status = rest_get_status($jobid);\n		print STD\
ERR \"$status\\n\" if ( $outputLevel > 0 );\n		if \
( $status eq 'ERROR' ) {\n			$errorCount++;\n		}\n\
		elsif ( $errorCount > 0 ) {\n			$errorCount--;\n\
		}\n		if (   $status eq 'RUNNING'\n			|| $status \
eq 'PENDING'\n			|| $status eq 'ERROR' )\n		{\n\n	\
		# Wait before polling again.\n			sleep $checkInt\
erval;\n		}\n	}\n	print_debug_message( 'client_pol\
l', 'End', 1 );\n	return $status;\n}\n\n=head2 get\
_results()\n\nGet the results for a job identifier\
.\n\n  &get_results($job_id);\n\n=cut\n\nsub get_r\
esults {\n	print_debug_message( 'get_results', 'Be\
gin', 1 );\n	my $jobid = shift;\n	print_debug_mess\
age( 'get_results', 'jobid: ' . $jobid, 1 );\n\n	#\
 Verbose\n	if ( $outputLevel > 1 ) {\n		print 'Get\
ting results for job ', $jobid, \"\\n\";\n	}\n\n	#\
 Check status, and wait if not finished\n	client_p\
oll($jobid);\n\n	# Use JobId if output file name i\
s not defined\n	unless ( defined( $params{'outfile\
'} ) ) {\n		$params{'outfile'} = $jobid;\n	}\n\n	#\
 Get list of data types\n	my (@resultTypes) = rest\
_get_result_types($jobid);\n\n	# Get the data and \
write it to a file\n	if ( defined( $params{'outfor\
mat'} ) ) {    # Specified data type\n		my $selRes\
ultType;\n		foreach my $resultType (@resultTypes) \
{\n			if ( $resultType->{'identifier'} eq $params{\
'outformat'} ) {\n				$selResultType = $resultType\
;\n			}\n		}\n		if ( defined($selResultType) ) {\n\
			my $result =\n			  rest_get_result( $jobid, $se\
lResultType->{'identifier'} );\n			if ( $params{'o\
utfile'} eq '-' ) {\n				write_file( $params{'outf\
ile'}, $result );\n			}\n			else {\n				write_file\
(\n					$params{'outfile'} . '.'\n					  . $selRes\
ultType->{'identifier'} . '.'\n					  . $selResult\
Type->{'fileSuffix'},\n					$result\n				);\n			}\\
n		}\n		else {\n			die 'Error: unknown result form\
at \"' . $params{'outformat'} . '\"';\n		}\n	}\n	e\
lse {    # Data types available\n		      # Write a\
 file for each output type\n		for my $resultType (\
@resultTypes) {\n			if ( $outputLevel > 1 ) {\n			\
	print STDERR 'Getting ', $resultType->{'identifie\
r'}, \"\\n\";\n			}\n			my $result = rest_get_resu\
lt( $jobid, $resultType->{'identifier'} );\n			if \
( $params{'outfile'} eq '-' ) {\n				write_file( $\
params{'outfile'}, $result );\n			}\n			else {\n		\
		write_file(\n					$params{'outfile'} . '.'\n				\
	  . $resultType->{'identifier'} . '.'\n					  . $\
resultType->{'fileSuffix'},\n					$result\n				);\\
n			}\n		}\n	}\n	print_debug_message( 'get_results\
', 'End', 1 );\n}\n\n=head2 read_file()\n\nRead a \
file into a scalar. The special filename '-' can b\
e used to read from \nstandard input (STDIN).\n\n \
 my $data = &read_file($filename);\n\n=cut\n\nsub \
read_file {\n	print_debug_message( 'read_file', 'B\
egin', 1 );\n	my $filename = shift;\n	print_debug_\
message( 'read_file', 'filename: ' . $filename, 2 \
);\n	my ( $content, $buffer );\n	if ( $filename eq\
 '-' ) {\n		while ( sysread( STDIN, $buffer, 1024 \
) ) {\n			$content .= $buffer;\n		}\n	}\n	else {  \
  # File\n		open( my $FILE, '<', $filename )\n		  \
or die \"Error: unable to open input file $filenam\
e ($!)\";\n		while ( sysread( $FILE, $buffer, 1024\
 ) ) {\n			$content .= $buffer;\n		}\n		close($FIL\
E);\n	}\n	print_debug_message( 'read_file', 'End',\
 1 );\n	return $content;\n}\n\n=head2 write_file()\
\n\nWrite data to a file. The special filename '-'\
 can be used to write to \nstandard output (STDOUT\
).\n\n  &write_file($filename, $data);\n\n=cut\n\n\
sub write_file {\n	print_debug_message( 'write_fil\
e', 'Begin', 1 );\n	my ( $filename, $data ) = @_;\\
n	print_debug_message( 'write_file', 'filename: ' \
. $filename, 2 );\n	if ( $outputLevel > 0 ) {\n		p\
rint STDERR 'Creating result file: ' . $filename .\
 \"\\n\";\n	}\n	if ( $filename eq '-' ) {\n		print\
 STDOUT $data;\n	}\n	else {\n		open( my $FILE, '>'\
, $filename )\n		  or die \"Error: unable to open \
output file $filename ($!)\";\n		syswrite( $FILE, \
$data );\n		close($FILE);\n	}\n	print_debug_messag\
e( 'write_file', 'End', 1 );\n}\n\n=head2 usage()\\
n\nPrint program usage message.\n\n  &usage();\n\n\
=cut\n\nsub usage {\n	print STDERR <<EOF\nNCBI BLA\
ST\n==========\n   \nRapid sequence database searc\
h programs utilizing the BLAST algorithm\n    \n[R\
equired]\n\n  -p, --program      : str  : BLAST pr\
ogram to use, see --paramDetail program\n  -D, --d\
atabase     : str  : database(s) to search, space \
separated. See\n                              --pa\
ramDetail database\n      --stype        : str  : \
query sequence type, see --paramDetail stype\n  se\
qFile            : file : query sequence (\"-\" fo\
r STDIN, \\@filename for\n                        \
      identifier list file)\n\n[Optional]\n\n  -m,\
 --matrix       : str  : scoring matrix, see --par\
amDetail matrix\n  -e, --exp          : real : 0<E\
<= 1000. Statistical significance threshold \n    \
                          for reporting database s\
equence matches.\n  -f, --filter       :      : fi\
lter the query sequence for low complexity \n     \
                         regions, see --paramDetai\
l filter\n  -A, --align        : int  : pairwise a\
lignment format, see --paramDetail align\n  -s, --\
scores       : int  : number of scores to be repor\
ted\n  -n, --alignments   : int  : number of align\
ments to report\n  -u, --match        : int  : Mat\
ch score (BLASTN only)\n  -v, --mismatch     : int\
  : Mismatch score (BLASTN only)\n  -o, --gapopen \
     : int  : Gap open penalty\n  -x, --gapext    \
   : int  : Gap extension penalty\n  -d, --dropoff\
      : int  : Drop-off\n  -g, --gapalign     :   \
   : Optimise gapped alignments\n      --seqrange \
    : str  : region within input to use as query\n\
      --multifasta   :      : treat input as a set\
 of fasta formatted sequences\n\n[General]\n\n  -h\
, --help        :      : prints this help text\n  \
    --async       :      : forces to make an async\
hronous query\n      --email       : str  : e-mail\
 address\n      --title       : str  : title for j\
ob\n      --status      :      : get job status\n \
     --resultTypes :      : get available result t\
ypes for job\n      --polljob     :      : poll fo\
r the status of a job\n      --jobid       : str  \
: jobid that was returned when an asynchronous job\
 \n                             was submitted.\n  \
    --outfile     : str  : file name for results (\
default is jobid;\n                             \"\
-\" for STDOUT)\n      --outformat   : str  : resu\
lt format to retrieve\n      --params      :      \
: list input parameters\n      --paramDetail : str\
  : display details for input parameter\n      --q\
uiet       :      : decrease output\n      --verbo\
se     :      : increase output\n      --trace    \
   :      : show SOAP messages being interchanged \
\n   \nSynchronous job:\n\n  The results/errors ar\
e returned as soon as the job is finished.\n  Usag\
e: $scriptName --email <your\\@email> [options...]\
 seqFile\n  Returns: results as an attachment\n\nA\
synchronous job:\n\n  Use this if you want to retr\
ieve the results at a later time. The results \n  \
are stored for up to 24 hours. 	\n  Usage: $script\
Name --async --email <your\\@email> [options...] s\
eqFile\n  Returns: jobid\n\n  Use the jobid to que\
ry for the status of the job. If the job is finish\
ed, \n  it also returns the results/errors.\n  Usa\
ge: $scriptName --polljob --jobid <jobId> [--outfi\
le string]\n  Returns: string indicating the statu\
s of the job and if applicable, results \n  as an \
attachment.\n\nFurther information:\n\n  http://ww\
w.ebi.ac.uk/Tools/webservices/services/sss/ncbi_bl\
ast_rest\n  http://www.ebi.ac.uk/Tools/webservices\
/tutorials/perl\n\nSupport/Feedback:\n\n  http://w\
ww.ebi.ac.uk/support/\nEOF\n}\n\n=head1 FEEDBACK/S\
UPPORT\n\nPlease contact us at L<http://www.ebi.ac\
.uk/support/> if you have any \nfeedback, suggesti\
ons or issues with the service or this client.\n\n\
=cut\n","\n=head1 NAME\n\nwublast_lwp.pl\n\n=head1\
 DESCRIPTION\n\nWU-BLAST REST web service Perl cli\
ent using L<LWP>.\n\nTested with:\n\n=over\n\n=ite\
m *\nL<LWP> 5.79, L<XML::Simple> 2.12 and Perl 5.8\
.3\n\n=item *\nL<LWP> 5.805, L<XML::Simple> 2.14 a\
nd Perl 5.8.7\n\n=item *\nL<LWP> 5.820, L<XML::Sim\
ple> 2.18 and Perl 5.10.0 (Ubuntu 9.04)\n\n=back\n\
\nFor further information see:\n\n=over\n\n=item *\
\nL<http://www.ebi.ac.uk/Tools/webservices/service\
s/sss/wu_blast_rest>\n\n=item *\nL<http://www.ebi.\
ac.uk/Tools/webservices/tutorials/perl>\n\n=back\n\
\n=head1 VERSION\n\n$Id: wublast_lwp.pl 1317 2009-\
09-03 15:44:11Z hpm $\n\n=cut\n\nuse strict;\nuse \
warnings;\n\nuse English;\nuse LWP;\nuse XML::Simp\
le;\nuse Getopt::Long qw(:config no_ignore_case bu\
ndling);\nuse File::Basename;\nuse Data::Dumper;\n\
\nmy $baseUrl = 'http://www.ebi.ac.uk/Tools/servic\
es/rest/wublast';\n\nmy $checkInterval = 3;\n\nmy \
$outputLevel = 1;\n\nmy $numOpts = scalar(@ARGV);\\
nmy %params = ( 'debugLevel' => 0 );\n\nmy %tool_p\
arams = ();\nGetOptions(\n\n	# Tool specific optio\
ns\n	'program|p=s'     => \\$tool_params{'program'\
},      # BLAST program\n	'database|D=s'    => \\$\
params{'database'},     # Search database\n	'matri\
x|m=s'      => \\$tool_params{'matrix'},       # S\
coring matrix\n	'exp|E=f'         => \\$tool_param\
s{'exp'},          # E-value threshold\n	'viewfilt\
er|e'    => \\$tool_params{'viewfilter'},   # Disp\
lay filtered sequence\n	'filter|f=s'      => \\$to\
ol_params{'filter'},       # Low complexity filter\
 name\n	'alignments|n=i'  => \\$tool_params{'align\
ments'},   # Number of alignments\n	'scores|s=i'  \
    => \\$tool_params{'scores'},       # Number of\
 scores\n	'sensitivity|S=s' => \\$tool_params{'sen\
sitivity'},  # Search sensitivity\n	'sort|t=s'    \
    => \\$tool_params{'sort'},         # Sort hits\
 by...\n	'stats|T=s'       => \\$tool_params{'stat\
s'},        # Scoring statistic to use\n	'strand|d\
=s'      => \\$tool_params{'strand'},       # Stra\
nd to use\n	'topcombon|c=i'   => \\$tool_params{'t\
opcombon'},    # Consistent sets of HSPs\n	'align|\
A=i'       => \\$tool_params{'align'},   # Pairwis\
e alignment format\n	'stype=s' => \\$tool_params{'\
stype'},    # Sequence type 'protein' or 'dna'\n	'\
sequence=s' => \\$params{'sequence'},         # Qu\
ery sequence file or DB:ID\n	'multifasta' => \\$pa\
rams{'multifasta'},       # Multiple fasta input\n\
\n	# Compatability options, old command-line.\n	'e\
chofilter|e'    => \\$params{'echofilter'},   # Di\
splay filtered sequence\n	'b=i'  => \\$params{'num\
al'},        # Number of alignments\n	'appxml=s'  \
      => \\$params{'appxml'},       # Application \
XML\n\n	# Generic options\n	'email=s'       => \\$\
params{'email'},          # User e-mail address\n	\
'title=s'       => \\$params{'title'},          # \
Job title\n	'outfile=s'     => \\$params{'outfile'\
},        # Output file name\n	'outformat=s'   => \
\\$params{'outformat'},      # Output file type\n	\
'jobid=s'       => \\$params{'jobid'},          # \
JobId\n	'help|h'        => \\$params{'help'},     \
      # Usage help\n	'async'         => \\$params{\
'async'},          # Asynchronous submission\n	'po\
lljob'       => \\$params{'polljob'},        # Get\
 results\n	'resultTypes'   => \\$params{'resultTyp\
es'},    # Get result types\n	'status'        => \\
\$params{'status'},         # Get status\n	'params\
'        => \\$params{'params'},         # List in\
put parameters\n	'paramDetail=s' => \\$params{'par\
amDetail'},    # Get details for parameter\n	'quie\
t'         => \\$params{'quiet'},          # Decre\
ase output level\n	'verbose'       => \\$params{'v\
erbose'},        # Increase output level\n	'debugL\
evel=i'  => \\$params{'debugLevel'},     # Debug o\
utput level\n	'baseUrl=s'     => \\$baseUrl,      \
            # Base URL for service.\n);\nif ( $par\
ams{'verbose'} ) { $outputLevel++ }\nif ( $params{\
'$quiet'} )  { $outputLevel-- }\n\n&print_debug_me\
ssage( 'MAIN', 'LWP::VERSION: ' . $LWP::VERSION,\n\
	1 );\n\n&print_debug_message( 'MAIN', \"params:\\\
n\" . Dumper( \\%params ),           11 );\n&print\
_debug_message( 'MAIN', \"tool_params:\\n\" . Dump\
er( \\%tool_params ), 11 );\n\nmy $scriptName = ba\
sename( $0, () );\n\nif ( $params{'help'} || $numO\
pts == 0 ) {\n	&usage();\n	exit(0);\n}\n\n&print_d\
ebug_message( 'MAIN', 'baseUrl: ' . $baseUrl, 1 );\
\n\nif (\n	!(\n		   $params{'polljob'}\n		|| $para\
ms{'resultTypes'}\n		|| $params{'status'}\n		|| $p\
arams{'params'}\n		|| $params{'paramDetail'}\n	)\n\
	&& !( defined( $ARGV[0] ) || defined( $params{'se\
quence'} ) )\n  )\n{\n\n	# Bad argument combinatio\
n, so print error message and usage\n	print STDERR\
 'Error: bad option combination', \"\\n\";\n	&usag\
e();\n	exit(1);\n}\n\nelsif ( $params{'params'} ) \
{\n	&print_tool_params();\n}\n\nelsif ( $params{'p\
aramDetail'} ) {\n	&print_param_details( $params{'\
paramDetail'} );\n}\n\nelsif ( $params{'status'} &\
& defined( $params{'jobid'} ) ) {\n	&print_job_sta\
tus( $params{'jobid'} );\n}\n\nelsif ( $params{'re\
sultTypes'} && defined( $params{'jobid'} ) ) {\n	&\
print_result_types( $params{'jobid'} );\n}\n\nelsi\
f ( $params{'polljob'} && defined( $params{'jobid'\
} ) ) {\n	&get_results( $params{'jobid'} );\n}\n\n\
else {\n\n	# Multiple input sequence mode, assume \
fasta format.\n	if ( $params{'multifasta'} ) {\n		\
&multi_submit_job();\n	}\n\n	# Entry identifier li\
st file.\n	elsif (( defined( $params{'sequence'} )\
 && $params{'sequence'} =~ m/^\\@/ )\n		|| ( defin\
ed( $ARGV[0] ) && $ARGV[0] =~ m/^\\@/ ) )\n	{\n		m\
y $list_filename = $params{'sequence'} || $ARGV[0]\
;\n		$list_filename =~ s/^\\@//;\n		&list_file_sub\
mit_job($list_filename);\n	}\n\n	# Default: single\
 sequence/identifier.\n	else {\n\n		# Load the seq\
uence data and submit.\n		&submit_job( &load_data(\
) );\n	}\n}\n\n=head1 FUNCTIONS\n\n=cut\n\n\n=head\
2 rest_request()\n\nPerform a REST request.\n\n  m\
y $response_str = &rest_request($url);\n\n=cut\n\n\
sub rest_request {\n	print_debug_message( 'rest_re\
quest', 'Begin', 11 );\n	my $requestUrl = shift;\n\
	print_debug_message( 'rest_request', 'URL: ' . $r\
equestUrl, 11 );\n\n	# Create a user agent\n	my $u\
a = LWP::UserAgent->new();\n	'$Revision: 1317 $' =\
~ m/(\\d+)/;\n	$ua->agent(\"EBI-Sample-Client/$1 (\
$scriptName; $OSNAME) \" . $ua->agent());\n	$ua->e\
nv_proxy;\n\n	# Perform the request\n	my $response\
 = $ua->get($requestUrl);\n	print_debug_message( '\
rest_request', 'HTTP status: ' . $response->code,\\
n		11 );\n\n	# Check for HTTP error codes\n	if ( $\
response->is_error ) {\n		$response->content() =~ \
m/<h1>([^<]+)<\\/h1>/;\n		die 'http status: ' . $r\
esponse->code . ' ' . $response->message . '  ' . \
$1;\n	}\n	print_debug_message( 'rest_request', 'En\
d', 11 );\n\n	# Return the response data\n	return \
$response->content();\n}\n\n=head2 rest_get_parame\
ters()\n\nGet list of tool parameter names.\n\n  m\
y (@param_list) = &rest_get_parameters();\n\n=cut\\
n\nsub rest_get_parameters {\n	print_debug_message\
( 'rest_get_parameters', 'Begin', 1 );\n	my $url  \
              = $baseUrl . '/parameters/';\n	my $p\
aram_list_xml_str = rest_request($url);\n	my $para\
m_list_xml     = XMLin($param_list_xml_str);\n	my \
(@param_list)       = @{ $param_list_xml->{'id'} }\
;\n	print_debug_message( 'rest_get_parameters', 'E\
nd', 1 );\n	return (@param_list);\n}\n\n=head2 res\
t_get_parameter_details()\n\nGet details of a tool\
 parameter.\n\n  my $paramDetail = &rest_get_param\
eter_details($param_name);\n\n=cut\n\nsub rest_get\
_parameter_details {\n	print_debug_message( 'rest_\
get_parameter_details', 'Begin', 1 );\n	my $parame\
terId = shift;\n	print_debug_message( 'rest_get_pa\
rameter_details',\n		'parameterId: ' . $parameterI\
d, 1 );\n	my $url                  = $baseUrl . '/\
parameterdetails/' . $parameterId;\n	my $param_det\
ail_xml_str = rest_request($url);\n	my $param_deta\
il_xml     = XMLin($param_detail_xml_str);\n	print\
_debug_message( 'rest_get_parameter_details', 'End\
', 1 );\n	return ($param_detail_xml);\n}\n\n=head2\
 rest_run()\n\nSubmit a job.\n\n  my $job_id = &re\
st_run($email, $title, \\%params );\n\n=cut\n\nsub\
 rest_run {\n	print_debug_message( 'rest_run', 'Be\
gin', 1 );\n	my $email  = shift;\n	my $title  = sh\
ift;\n	my $params = shift;\n	print_debug_message( \
'rest_run', 'email: ' . $email, 1 );\n	if ( define\
d($title) ) {\n		print_debug_message( 'rest_run', \
'title: ' . $title, 1 );\n	}\n	print_debug_message\
( 'rest_run', 'params: ' . Dumper($params), 1 );\n\
\n	# User agent to perform http requests\n	my $ua \
= LWP::UserAgent->new();\n	$ua->env_proxy;\n\n	# C\
lean up parameters\n	my (%tmp_params) = %{$params}\
;\n	$tmp_params{'email'} = $email;\n	$tmp_params{'\
title'} = $title;\n	foreach my $param_name ( keys(\
%tmp_params) ) {\n		if ( !defined( $tmp_params{$pa\
ram_name} ) ) {\n			delete $tmp_params{$param_name\
};\n		}\n	}\n\n	# Submit the job as a POST\n	my $u\
rl = $baseUrl . '/run';\n	my $response = $ua->post\
( $url, \\%tmp_params );\n	print_debug_message( 'r\
est_run', 'HTTP status: ' . $response->code, 11 );\
\n	print_debug_message( 'rest_run',\n		'request: '\
 . $response->request()->content(), 11 );\n\n	# Ch\
eck for HTTP error codes\n	if ( $response->is_erro\
r ) {\n		$response->content() =~ m/<h1>([^<]+)<\\/\
h1>/;\n		die 'http status: ' . $response->code . '\
 ' . $response->message . '  ' . $1;\n	}\n\n	# The\
 job id is returned\n	my $job_id = $response->cont\
ent();\n	print_debug_message( 'rest_run', 'End', 1\
 );\n	return $job_id;\n}\n\n=head2 rest_get_status\
()\n\nCheck the status of a job.\n\n  my $status =\
 &rest_get_status($job_id);\n\n=cut\n\nsub rest_ge\
t_status {\n	print_debug_message( 'rest_get_status\
', 'Begin', 1 );\n	my $job_id = shift;\n	print_deb\
ug_message( 'rest_get_status', 'jobid: ' . $job_id\
, 2 );\n	my $status_str = 'UNKNOWN';\n	my $url    \
    = $baseUrl . '/status/' . $job_id;\n	$status_s\
tr = &rest_request($url);\n	print_debug_message( '\
rest_get_status', 'status_str: ' . $status_str, 2 \
);\n	print_debug_message( 'rest_get_status', 'End'\
, 1 );\n	return $status_str;\n}\n\n=head2 rest_get\
_result_types()\n\nGet list of result types for fi\
nished job.\n\n  my (@result_types) = &rest_get_re\
sult_types($job_id);\n\n=cut\n\nsub rest_get_resul\
t_types {\n	print_debug_message( 'rest_get_result_\
types', 'Begin', 1 );\n	my $job_id = shift;\n	prin\
t_debug_message( 'rest_get_result_types', 'jobid: \
' . $job_id, 2 );\n	my (@resultTypes);\n	my $url  \
                    = $baseUrl . '/resulttypes/' .\
 $job_id;\n	my $result_type_list_xml_str = &rest_r\
equest($url);\n	my $result_type_list_xml     = XML\
in($result_type_list_xml_str);\n	(@resultTypes) = \
@{ $result_type_list_xml->{'type'} };\n	print_debu\
g_message( 'rest_get_result_types',\n		scalar(@res\
ultTypes) . ' result types', 2 );\n	print_debug_me\
ssage( 'rest_get_result_types', 'End', 1 );\n	retu\
rn (@resultTypes);\n}\n\n=head2 rest_get_result()\\
n\nGet result data of a specified type for a finis\
hed job.\n\n  my $result = rest_get_result($job_id\
, $result_type);\n\n=cut\n\nsub rest_get_result {\\
n	print_debug_message( 'rest_get_result', 'Begin',\
 1 );\n	my $job_id = shift;\n	my $type   = shift;\\
n	print_debug_message( 'rest_get_result', 'jobid: \
' . $job_id, 1 );\n	print_debug_message( 'rest_get\
_result', 'type: ' . $type,    1 );\n	my $url    =\
 $baseUrl . '/result/' . $job_id . '/' . $type;\n	\
my $result = &rest_request($url);\n	print_debug_me\
ssage( 'rest_get_result', length($result) . ' char\
acters',\n		1 );\n	print_debug_message( 'rest_get_\
result', 'End', 1 );\n	return $result;\n}\n\n\n=he\
ad2 print_debug_message()\n\nPrint debug message a\
t specified debug level.\n\n  &print_debug_message\
($method_name, $message, $level);\n\n=cut\n\nsub p\
rint_debug_message {\n	my $function_name = shift;\\
n	my $message       = shift;\n	my $level         =\
 shift;\n	if ( $level <= $params{'debugLevel'} ) {\
\n		print STDERR '[', $function_name, '()] ', $mes\
sage, \"\\n\";\n	}\n}\n\n=head2 print_tool_params(\
)\n\nPrint list of tool parameters.\n\n  &print_to\
ol_params();\n\n=cut\n\nsub print_tool_params {\n	\
print_debug_message( 'print_tool_params', 'Begin',\
 1 );\n	my (@param_list) = &rest_get_parameters();\
\n	foreach my $param ( sort(@param_list) ) {\n		pr\
int $param, \"\\n\";\n	}\n	print_debug_message( 'p\
rint_tool_params', 'End', 1 );\n}\n\n=head2 print_\
param_details()\n\nPrint details of a tool paramet\
er.\n\n  &print_param_details($param_name);\n\n=cu\
t\n\nsub print_param_details {\n	print_debug_messa\
ge( 'print_param_details', 'Begin', 1 );\n	my $par\
amName = shift;\n	print_debug_message( 'print_para\
m_details', 'paramName: ' . $paramName, 2 );\n	my \
$paramDetail = &rest_get_parameter_details($paramN\
ame);\n	print $paramDetail->{'name'}, \"\\t\", $pa\
ramDetail->{'type'}, \"\\n\";\n	print $paramDetail\
->{'description'}, \"\\n\";\n	foreach my $value ( \
@{ $paramDetail->{'values'}->{'value'} } ) {\n		pr\
int $value->{'value'};\n		if ( $value->{'defaultVa\
lue'} eq 'true' ) {\n			print \"\\t\", 'default';\\
n		}\n		print \"\\n\";\n		print \"\\t\", $value->{\
'label'}, \"\\n\";\n	}\n	print_debug_message( 'pri\
nt_param_details', 'End', 1 );\n}\n\n=head2 print_\
job_status()\n\nPrint status of a job.\n\n  &print\
_job_status($job_id);\n\n=cut\n\nsub print_job_sta\
tus {\n	print_debug_message( 'print_job_status', '\
Begin', 1 );\n	my $jobid = shift;\n	print_debug_me\
ssage( 'print_job_status', 'jobid: ' . $jobid, 1 )\
;\n	if ( $outputLevel > 0 ) {\n		print STDERR 'Get\
ting status for job ', $jobid, \"\\n\";\n	}\n	my $\
result = &rest_get_status($jobid);\n	print \"$resu\
lt\\n\";\n	if ( $result eq 'FINISHED' && $outputLe\
vel > 0 ) {\n		print STDERR \"To get results: $scr\
iptName --polljob --jobid \" . $jobid\n		  . \"\\n\
\";\n	}\n	print_debug_message( 'print_job_status',\
 'End', 1 );\n}\n\n=head2 print_result_types()\n\n\
Print available result types for a job.\n\n  &prin\
t_result_types($job_id);\n\n=cut\n\nsub print_resu\
lt_types {\n	print_debug_message( 'result_types', \
'Begin', 1 );\n	my $jobid = shift;\n	print_debug_m\
essage( 'result_types', 'jobid: ' . $jobid, 1 );\n\
	if ( $outputLevel > 0 ) {\n		print STDERR 'Gettin\
g result types for job ', $jobid, \"\\n\";\n	}\n	m\
y $status = &rest_get_status($jobid);\n	if ( $stat\
us eq 'PENDING' || $status eq 'RUNNING' ) {\n		pri\
nt STDERR 'Error: Job status is ', $status,\n		  '\
. To get result types the job must be finished.', \
\"\\n\";\n	}\n	else {\n		my (@resultTypes) = &rest\
_get_result_types($jobid);\n		if ( $outputLevel > \
0 ) {\n			print STDOUT 'Available result types:', \
\"\\n\";\n		}\n		foreach my $resultType (@resultTy\
pes) {\n			print STDOUT $resultType->{'identifier'\
}, \"\\n\";\n			if ( defined( $resultType->{'label\
'} ) ) {\n				print STDOUT \"\\t\", $resultType->{\
'label'}, \"\\n\";\n			}\n			if ( defined( $result\
Type->{'description'} ) ) {\n				print STDOUT \"\\\
t\", $resultType->{'description'}, \"\\n\";\n			}\\
n			if ( defined( $resultType->{'mediaType'} ) ) {\
\n				print STDOUT \"\\t\", $resultType->{'mediaTy\
pe'}, \"\\n\";\n			}\n			if ( defined( $resultType\
->{'fileSuffix'} ) ) {\n				print STDOUT \"\\t\", \
$resultType->{'fileSuffix'}, \"\\n\";\n			}\n		}\n\
		if ( $status eq 'FINISHED' && $outputLevel > 0 )\
 {\n			print STDERR \"\\n\", 'To get results:', \"\
\\n\",\n			  \"  $scriptName --polljob --jobid \" \
. $params{'jobid'} . \"\\n\",\n			  \"  $scriptNam\
e --polljob --outformat <type> --jobid \"\n			  . \
$params{'jobid'} . \"\\n\";\n		}\n	}\n	print_debug\
_message( 'result_types', 'End', 1 );\n}\n\n=head2\
 submit_job()\n\nSubmit a job to the service.\n\n \
 &submit_job($seq);\n\n=cut\n\nsub submit_job {\n	\
print_debug_message( 'submit_job', 'Begin', 1 );\n\
\n	# Set input sequence\n	$tool_params{'sequence'}\
 = shift;\n\n	# Load parameters\n	&load_params();\\
n\n	# Submit the job\n	my $jobid = &rest_run( $par\
ams{'email'}, $params{'title'}, \\%tool_params );\\
n\n	# Simulate sync/async mode\n	if ( defined( $pa\
rams{'async'} ) ) {\n		print STDOUT $jobid, \"\\n\\
";\n		if ( $outputLevel > 0 ) {\n			print STDERR\n\
			  \"To check status: $scriptName --status --job\
id $jobid\\n\";\n		}\n	}\n	else {\n		if ( $outputL\
evel > 0 ) {\n			print STDERR \"JobId: $jobid\\n\"\
;\n		}\n		sleep 1;\n		&get_results($jobid);\n	}\n	\
print_debug_message( 'submit_job', 'End', 1 );\n}\\
n\n=head2 multi_submit_job()\n\nSubmit multiple jo\
bs assuming input is a collection of fasta formatt\
ed sequences.\n\n  &multi_submit_job();\n\n=cut\n\\
nsub multi_submit_job {\n	print_debug_message( 'mu\
lti_submit_job', 'Begin', 1 );\n	my $jobIdForFilen\
ame = 1;\n	$jobIdForFilename = 0 if ( defined( $pa\
rams{'outfile'} ) );\n	my (@filename_list) = ();\n\
\n	# Query sequence\n	if ( defined( $ARGV[0] ) ) {\
    # Bare option\n		if ( -f $ARGV[0] || $ARGV[0] \
eq '-' ) {    # File\n			push( @filename_list, $AR\
GV[0] );\n		}\n	}\n	if ( $params{'sequence'} ) {  \
                 # Via --sequence\n		if ( -f $para\
ms{'sequence'} || $params{'sequence'} eq '-' ) {  \
  # File\n			push( @filename_list, $params{'sequen\
ce'} );\n		}\n	}\n\n	$/ = '>';\n	foreach my $filen\
ame (@filename_list) {\n		open( my $INFILE, '<', $\
filename )\n		  or die \"Error: unable to open fil\
e $filename ($!)\";\n		while (<$INFILE>) {\n			my \
$seq = $_;\n			$seq =~ s/>$//;\n			if ( $seq =~ m/\
(\\S+)/ ) {\n				print STDERR \"Submitting job for\
: $1\\n\"\n				  if ( $outputLevel > 0 );\n				$se\
q = '>' . $seq;\n				&print_debug_message( 'multi_\
submit_job', $seq, 11 );\n				&submit_job($seq);\n\
				$params{'outfile'} = undef if ( $jobIdForFilen\
ame == 1 );\n			}\n		}\n		close $INFILE;\n	}\n	pri\
nt_debug_message( 'multi_submit_job', 'End', 1 );\\
n}\n\n=head2 list_file_submit_job()\n\nSubmit mult\
iple jobs using a file containing a list of entry \
identifiers as \ninput.\n\n  &list_file_submit_job\
($list_filename)\n\n=cut\n\nsub list_file_submit_j\
ob {\n	my $filename         = shift;\n	my $jobIdFo\
rFilename = 1;\n	$jobIdForFilename = 0 if ( define\
d( $params{'outfile'} ) );\n\n	# Iterate over iden\
tifiers, submitting each job\n	open( my $LISTFILE,\
 '<', $filename )\n	  or die 'Error: unable to ope\
n file ' . $filename . ' (' . $! . ')';\n	while (<\
$LISTFILE>) {\n		my $line = $_;\n		chomp($line);\n\
		if ( $line ne '' ) {\n			&print_debug_message( '\
list_file_submit_job', 'line: ' . $line, 2 );\n			\
if ( $line =~ m/\\w:\\w/ ) {    # Check this is an\
 identifier\n				print STDERR \"Submitting job for\
: $line\\n\"\n				  if ( $outputLevel > 0 );\n				\
&submit_job($line);\n			}\n			else {\n				print ST\
DERR\n\"Warning: line \\\"$line\\\" is not recogni\
sed as an identifier\\n\";\n			}\n		}\n		$params{'\
outfile'} = undef if ( $jobIdForFilename == 1 );\n\
	}\n	close $LISTFILE;\n}\n\n=head2 load_data()\n\n\
Load sequence data from file or option specified o\
n the command-line.\n\n  &load_data();\n\n=cut\n\n\
sub load_data {\n	print_debug_message( 'load_data'\
, 'Begin', 1 );\n	my $retSeq;\n\n	# Query sequence\
\n	if ( defined( $ARGV[0] ) ) {    # Bare option\n\
		if ( -f $ARGV[0] || $ARGV[0] eq '-' ) {    # Fil\
e\n			$retSeq = &read_file( $ARGV[0] );\n		}\n		el\
se {                                     # DB:ID o\
r sequence\n			$retSeq = $ARGV[0];\n		}\n	}\n	if (\
 $params{'sequence'} ) {                   # Via -\
-sequence\n		if ( -f $params{'sequence'} || $param\
s{'sequence'} eq '-' ) {    # File\n			$retSeq = &\
read_file( $params{'sequence'} );\n		}\n		else {  \
  # DB:ID or sequence\n			$retSeq = $params{'seque\
nce'};\n		}\n	}\n	print_debug_message( 'load_data'\
, 'End', 1 );\n	return $retSeq;\n}\n\n=head2 load_\
params()\n\nLoad job parameters from command-line \
options.\n\n  &load_params();\n\n=cut\n\nsub load_\
params {\n	print_debug_message( 'load_params', 'Be\
gin', 1 );\n\n	# Database(s) to search\n	my (@dbLi\
st) = split /[ ,]/, $params{'database'};\n	$tool_p\
arams{'database'} = \\@dbList;\n\n	# Compatability\
 options, old command-line.\n	if(!$tool_params{'vi\
ewfilter'} && $params{'echofilter'}) {\n		$tool_pa\
rams{'viewfilter'} = 'true';\n	}\n	if(!$tool_param\
s{'alignments'} && $params{'numal'}) {\n		$tool_pa\
rams{'alignments'} = $params{'numal'};\n	}\n	# TOD\
O: set alignment format option to get NCBI BLAST X\
ML.\n	if($params{'appxml'}) {\n		$tool_params{'ali\
gn'} = '';\n	}\n\n	print_debug_message( 'load_para\
ms', 'End', 1 );\n}\n\n=head2 client_poll()\n\nCli\
ent-side job polling.\n\n  &client_poll($job_id);\\
n\n=cut\n\nsub client_poll {\n	print_debug_message\
( 'client_poll', 'Begin', 1 );\n	my $jobid  = shif\
t;\n	my $status = 'PENDING';\n\n	my $errorCount = \
0;\n	while ($status eq 'RUNNING'\n		|| $status eq \
'PENDING'\n		|| ( $status eq 'ERROR' && $errorCoun\
t < 2 ) )\n	{\n		$status = rest_get_status($jobid)\
;\n		print STDERR \"$status\\n\" if ( $outputLevel\
 > 0 );\n		if ( $status eq 'ERROR' ) {\n			$errorC\
ount++;\n		}\n		elsif ( $errorCount > 0 ) {\n			$e\
rrorCount--;\n		}\n		if (   $status eq 'RUNNING'\n\
			|| $status eq 'PENDING'\n			|| $status eq 'ERRO\
R' )\n		{\n\n			# Wait before polling again.\n			s\
leep $checkInterval;\n		}\n	}\n	print_debug_messag\
e( 'client_poll', 'End', 1 );\n	return $status;\n}\
\n\n=head2 get_results()\n\nGet the results for a \
job identifier.\n\n  &get_results($job_id);\n\n=cu\
t\n\nsub get_results {\n	print_debug_message( 'get\
_results', 'Begin', 1 );\n	my $jobid = shift;\n	pr\
int_debug_message( 'get_results', 'jobid: ' . $job\
id, 1 );\n\n	# Verbose\n	if ( $outputLevel > 1 ) {\
\n		print 'Getting results for job ', $jobid, \"\\\
n\";\n	}\n\n	# Check status, and wait if not finis\
hed\n	client_poll($jobid);\n\n	# Use JobId if outp\
ut file name is not defined\n	unless ( defined( $p\
arams{'outfile'} ) ) {\n		$params{'outfile'} = $jo\
bid;\n	}\n\n	# Get list of data types\n	my (@resul\
tTypes) = rest_get_result_types($jobid);\n\n	# Get\
 the data and write it to a file\n	if ( defined( $\
params{'outformat'} ) ) {    # Specified data type\
\n		my $selResultType;\n		foreach my $resultType (\
@resultTypes) {\n			if ( $resultType->{'identifier\
'} eq $params{'outformat'} ) {\n				$selResultType\
 = $resultType;\n			}\n		}\n		if ( defined($selRes\
ultType) ) {\n			my $result =\n			  rest_get_resul\
t( $jobid, $selResultType->{'identifier'} );\n			i\
f ( $params{'outfile'} eq '-' ) {\n				write_file(\
 $params{'outfile'}, $result );\n			}\n			else {\n\
				write_file(\n					$params{'outfile'} . '.'\n		\
			  . $selResultType->{'identifier'} . '.'\n					\
  . $selResultType->{'fileSuffix'},\n					$result\\
n				);\n			}\n		}\n		else {\n			die 'Error: unkno\
wn result format \"' . $params{'outformat'} . '\"'\
;\n		}\n	}\n	else {    # Data types available\n		 \
     # Write a file for each output type\n		for my\
 $resultType (@resultTypes) {\n			if ( $outputLeve\
l > 1 ) {\n				print STDERR 'Getting ', $resultTyp\
e->{'identifier'}, \"\\n\";\n			}\n			my $result =\
 rest_get_result( $jobid, $resultType->{'identifie\
r'} );\n			if ( $params{'outfile'} eq '-' ) {\n			\
	write_file( $params{'outfile'}, $result );\n			}\\
n			else {\n				write_file(\n					$params{'outfile\
'} . '.'\n					  . $resultType->{'identifier'} . '\
.'\n					  . $resultType->{'fileSuffix'},\n					$r\
esult\n				);\n			}\n		}\n	}\n	print_debug_message\
( 'get_results', 'End', 1 );\n}\n\n=head2 read_fil\
e()\n\nRead a file into a scalar. The special file\
name '-' can be used to read from \nstandard input\
 (STDIN).\n\n  my $data = &read_file($filename);\n\
\n=cut\n\nsub read_file {\n	print_debug_message( '\
read_file', 'Begin', 1 );\n	my $filename = shift;\\
n	print_debug_message( 'read_file', 'filename: ' .\
 $filename, 2 );\n	my ( $content, $buffer );\n	if \
( $filename eq '-' ) {\n		while ( sysread( STDIN, \
$buffer, 1024 ) ) {\n			$content .= $buffer;\n		}\\
n	}\n	else {    # File\n		open( my $FILE, '<', $fi\
lename )\n		  or die \"Error: unable to open input\
 file $filename ($!)\";\n		while ( sysread( $FILE,\
 $buffer, 1024 ) ) {\n			$content .= $buffer;\n		}\
\n		close($FILE);\n	}\n	print_debug_message( 'read\
_file', 'End', 1 );\n	return $content;\n}\n\n=head\
2 write_file()\n\nWrite data to a file. The specia\
l filename '-' can be used to write to \nstandard \
output (STDOUT).\n\n  &write_file($filename, $data\
);\n\n=cut\n\nsub write_file {\n	print_debug_messa\
ge( 'write_file', 'Begin', 1 );\n	my ( $filename, \
$data ) = @_;\n	print_debug_message( 'write_file',\
 'filename: ' . $filename, 2 );\n	if ( $outputLeve\
l > 0 ) {\n		print STDERR 'Creating result file: '\
 . $filename . \"\\n\";\n	}\n	if ( $filename eq '-\
' ) {\n		print STDOUT $data;\n	}\n	else {\n		open(\
 my $FILE, '>', $filename )\n		  or die \"Error: u\
nable to open output file $filename ($!)\";\n		sys\
write( $FILE, $data );\n		close($FILE);\n	}\n	prin\
t_debug_message( 'write_file', 'End', 1 );\n}\n\n=\
head2 usage()\n\nPrint program usage message.\n\n \
 &usage();\n\n=cut\n\nsub usage {\n	print STDERR <\
<EOF\nWU-BLAST\n========\n   \nRapid sequence data\
base search programs utilizing the BLAST algorithm\
\n    \n[Required]\n\n  -p, --program      : str  \
: BLAST program to use, see --paramDetail program\\
n  -D, --database     : str  : database(s) to sear\
ch, space separated. See\n                        \
      --paramDetail database\n      --stype       \
 : str  : query sequence type, see --paramDetail s\
type\n  seqFile            : file : query sequence\
 (\"-\" for STDIN, \\@filename for\n              \
                identifier list file)\n\n[Optional\
]\n\n  -m, --matrix       : str  : scoring matrix,\
 see --paramDetail matrix\n  -e, --exp          : \
real : 0<E<= 1000. Statistical significance thresh\
old \n                              for reporting \
database sequence matches.\n  -e, --viewfilter   :\
      : display the filtered query sequence\n  -f,\
 --filter       : str  : filter the query sequence\
 for low complexity \n                            \
  regions, see --paramDetail filter\n  -A, --align\
        : int  : pairwise alignment format, see --\
paramDetail align\n  -s, --scores       : int  : n\
umber of scores to be reported\n  -b, --alignments\
   : int  : number of alignments to report\n  -S, \
--sensitivity  : str  : sensitivity of the search,\
 \n                              see --paramDetail\
 sensitivity\n  -t, --sort	     : str  : sort orde\
r for hits, see --paramDetail sort\n  -T, --stats \
       : str  : statistical model, see --paramDeta\
il stats\n  -d, --strand       : str  : DNA strand\
 to search with,\n                              se\
e --paramDetail strand\n  -c, --topcombon    : str\
  : consistent sets of HSPs\n      --multifasta   \
:      : treat input as a set of fasta formatted s\
equences\n\n[General]\n\n  -h, --help        :    \
  : prints this help text\n      --async       :  \
    : forces to make an asynchronous query\n      \
--email       : str  : e-mail address\n      --tit\
le       : str  : title for job\n      --status   \
   :      : get job status\n      --resultTypes : \
     : get available result types for job\n      -\
-polljob     :      : poll for the status of a job\
\n      --jobid       : str  : jobid that was retu\
rned when an asynchronous job \n                  \
           was submitted.\n      --outfile     : s\
tr  : file name for results (default is jobid;\n  \
                           \"-\" for STDOUT)\n    \
  --outformat   : str  : result format to retrieve\
\n      --params      :      : list input paramete\
rs\n      --paramDetail : str  : display details f\
or input parameter\n      --quiet       :      : d\
ecrease output\n      --verbose     :      : incre\
ase output\n      --trace       :      : show SOAP\
 messages being interchanged \n   \nSynchronous jo\
b:\n\n  The results/errors are returned as soon as\
 the job is finished.\n  Usage: $scriptName --emai\
l <your\\@email> [options...] seqFile\n  Returns: \
results as an attachment\n\nAsynchronous job:\n\n \
 Use this if you want to retrieve the results at a\
 later time. The results \n  are stored for up to \
24 hours. 	\n  Usage: $scriptName --async --email \
<your\\@email> [options...] seqFile\n  Returns: jo\
bid\n\n  Use the jobid to query for the status of \
the job. If the job is finished, \n  it also retur\
ns the results/errors.\n  Usage: $scriptName --pol\
ljob --jobid <jobId> [--outfile string]\n  Returns\
: string indicating the status of the job and if a\
pplicable, results \n  as an attachment.\n\nFurthe\
r information:\n\n  http://www.ebi.ac.uk/Tools/web\
services/services/sss/wu_blast_rest\n  http://www.\
ebi.ac.uk/Tools/webservices/tutorials/perl\n\nSupp\
ort/Feedback:\n\n  http://www.ebi.ac.uk/support/\n\
EOF\n}\n\n=head1 FEEDBACK/SUPPORT\n\nPlease contac\
t us at L<http://www.ebi.ac.uk/support/> if you ha\
ve any \nfeedback, suggestions or issues with the \
service or this client.\n\n=cut\n","\n\n\nmy $PROB\
TRESH = 0.3;# base pairs below this prob threshold\
 will be ignored\nmy $WEIGHT = 100.0; # float!!\nm\
y $NUCALPH = \"ACGTUNRYMKSWHBVD\";\nuse vars qw($N\
UCALPH $WEIGHT);\n\nmy $myname = basename($0);\n\n\
use strict;\nuse warnings;\n\nuse File::Basename;\\
nuse Getopt::Long;\nuse File::Glob ':glob';\nuse F\
ile::Spec;\nuse File::Temp qw/ tempfile tempdir /;\
\n\n\n\n\nsub tcoffeelib_header($;$)\n{\n    my ($\
nseq, $fd) = @_;\n    if (! defined($fd)) {\n     \
   $fd = *STDOUT;\n    }\n    printf $fd \"! TC_LI\
B_FORMAT_01\\n\";\n    printf $fd \"%d\\n\", $nseq\
;\n}\n\n\nsub tcoffeelib_header_addseq($$;$)\n{\n \
   my ($id, $seq, $fd) = @_;\n    if (! defined($f\
d)) {\n        $fd = *STDOUT;\n    }\n    printf $\
fd \"%s %d %s\\n\", $id, length($seq), $seq;\n}\n\\
n\nsub tcoffeelib_comment($;$)\n{\n    my ($commen\
t, $fd) = @_;\n    if (! defined($fd)) {\n        \
$fd = *STDOUT;\n    }\n    printf $fd \"!\" . $com\
ment . \"\\n\";\n}\n\n\nsub tcoffeelib_struct($$$;\
$)\n{\n    my ($nseq, $len, $bpm, $fd) = @_;\n\n  \
  if (! defined($fd)) {\n        $fd = *STDOUT;\n \
   }\n\n    # output basepair indices with fixed w\
eight\n    printf $fd \"#%d %d\\n\", $nseq, $nseq;\
\n    # output basepairs (only once) and with unit\
-offset\n    for (my $i=0; $i<$len; $i++) {\n     \
   for (my $j=$i+1; $j<$len; $j++) {\n            \
if (! defined($bpm->[$i][$j])) {\n                \
print STDERR \"ERROR: \\$bpm->[$i][$j] undefined\\\
n\";\n            }\n            if ($bpm->[$i][$j\
]>0) {\n                print $fd $i+1;\n         \
       print $fd \" \";\n                print $fd\
 $j+1;\n                print $fd \" \" . $bpm->[$\
i][$j] . \"\\n\";\n            }\n        }\n    }\
\n}\n\n\nsub tcoffeelib_footer(;$)\n{\n    my ($fd\
) = @_;\n    if (! defined($fd)) {\n        $fd = \
*STDOUT;\n    }\n    print $fd \"! SEQ_1_TO_N\\n\"\
;\n}\n\n\n    \nsub plfold($$$)\n{    \n    my ($i\
d, $seq, $probtresh) = @_;\n    my (@struct);# ret\
urn\n    my ($templ, $fhtmp, $fnametmp, $cmd, $ctr\
, $window_size);\n    our $ntemp++;\n    \n    $te\
mpl = $myname . \".pid-\" . $$ .$ntemp .\".XXXXXX\\
";\n    ($fhtmp, $fnametmp) = tempfile($templ, UNL\
INK => 1); \n    print $fhtmp \">$id\\n$seq\\n\";\\
n\n    # --- init basepair array\n    #\n    for (\
my $i=0; $i<length($seq); $i++) {\n        for (my\
 $j=$i+1; $j<length($seq); $j++) {\n            $s\
truct[$i][$j]=0;\n        }\n    }\n\n\n    # --- \
call rnaplfold and drop a readme\n    #\n    $wind\
ow_size=(length($seq)<70)?length($seq):70;\n    $c\
md = \"RNAplfold -W $window_size < $fnametmp >/dev\
/null\";\n    system($cmd);\n    \n    if ($? != 0\
) {\n        printf STDERR \"ERROR: RNAplfold ($cm\
d) exited with error status %d\\n\", $? >> 8;\n   \
     return;\n    }\n    #unlink($fnametmp);\n    \
my $fps = sprintf(\"%s_dp.ps\", $id); # check long\
 name\n    \n    if (! -s $fps) {\n      {\n\n	$fp\
s = sprintf(\"%s_dp.ps\", substr($id,0,12)); # che\
ck short name\n 	if (! -s $fps)\n	  {\n	    die(\"\
couldn't find expected file $fps\\n\");\n	    retu\
rn;\n	  }\n      }\n    }\n\n    \n    # --- read \
base pairs from created postscript\n    #\n    ope\
n(FH, $fps);\n    while (my $line = <FH>) {\n     \
   my ($nti, $ntj, $prob);\n        chomp($line); \
       \n        # line: bp bp sqrt-prob ubox\n   \
     my @match = ($line =~ m/^([0-9]+) +([0-9]+) +\
([0-9\\.]+) +ubox$/);\n        if (scalar(@match))\
 {\n            $nti=$1;\n            $ntj=$2;\n  \
          $prob=$3*$3;# prob stored as square root\
\n\n            if ($prob>$probtresh) {\n         \
       #printf STDERR \"\\$struct[$nti][$ntj] sqrt\
prob=$3 prob=$prob > $probtresh\\n\";\n           \
     $struct[$nti-1][$ntj-1] = $WEIGHT\n          \
  }\n            # store with zero-offset\n       \
 }\n    }\n    close(FH);\n\n    # remove or gzi p\
ostscript\n    #\n    unlink($fps);\n    #\n    # \
or gzip\n    #$cmd = \"gzip -qf $fps\";\n    #syst\
em($cmd);\n    #if ($? != 0) {\n    #    printf ST\
DERR \"ERROR: gzip ($cmd) exited with error status\
 %d\\n\", $? >> 8;\n    #}\n\n    return \\@struct\
;\n}\n\n\n\n\n\nsub rnaseqfmt($)\n{\n    my ($seq)\
 = @_;\n    # remove gaps\n    $seq =~ s/-//g;\n  \
  # uppercase RNA\n    $seq = uc($seq);\n    # T -\
> U\n    $seq =~ s/T/U/g;\n    # check for invalid\
 charaters\n    $_ = $seq;\n    s/[^$NUCALPH]//g;\\
n    return $_;\n}\n\n\n\n\nsub usage(;$)\n{    \n\
    my ($errmsg) = @_;\n    if ($errmsg) {\n      \
  print STDERR \"ERROR: $errmsg\\n\";\n    }\n    \
print STDERR << \"EOF\";\n$myname:\n Creates a T-C\
offee RNA structure library from RNAplfold predict\
ion.\n See FIXME:citation\nUsage:\n $myname -in se\
q_file -out tcoffee_lib\nEOF\n    exit(1);\n}\n\ns\
ub read_fasta_seq \n  {\n    my $f=$_[0];\n    my \
%hseq;\n    my (@seq, @com, @name);\n    my ($a, $\
s,$nseq);\n\n    open (F, $f);\n    while (<F>)\n \
     {\n	$s.=$_;\n      }\n    close (F);\n\n    \\
n    @name=($s=~/>(\\S*).*\\n[^>]*/g);\n    \n    \
@seq =($s=~/>.*.*\\n([^>]*)/g);\n    @com =($s=~/>\
(\\S*)(.*)\\n([^>]*)/g);\n\n\n    $nseq=$#name+1;\\
n  \n    for ($a=0; $a<$nseq; $a++)\n      {\n	my \
$n=$name[$a];\n	my $s;\n	$hseq{$n}{name}=$n;\n	$s=\
$seq[$a];$s=~s/\\s//g;\n	\n	$hseq{$n}{seq}=$s;\n	$\
hseq{$n}{com}=$com[$a];\n      }\n    return %hseq\
;\n  }\n\n\n\n\n\n\n\nmy $fmsq = \"\";\nmy $flib =\
 \"\";\nmy %OPTS;\nmy %seq;\nmy ($id, $nseq, $i);\\
nmy @nl;\n\nGetOptions(\"in=s\" => \\$fmsq, \"out=\
s\" => \\$flib);\n\nif (! -s $fmsq) {\n    usage(\\
"empty or non-existant file \\\"$fmsq\\\"\")\n}\ni\
f (length($flib)==0) {\n    usage(\"empty out-file\
name\")\n}\n\n\n\n\n\n\n%seq=read_fasta_seq($fmsq)\
;\n\n\n@nl=keys(%seq);\n\n$nseq=$#nl+1;\nopen FD_L\
IB, \">$flib\" or die \"can't open $flib!\";\ntcof\
feelib_header($nseq, *FD_LIB);\nforeach $id (keys \
(%seq))\n  {\n    my ($seq, $fmtseq);\n    \n    $\
seq = $seq{$id}{seq};\n    \n    $fmtseq = rnaseqf\
mt($seq);# check here, formatting for folding impo\
rtant later\n    if (length($seq)!=length($fmtseq)\
) {\n        print STDERR \"ERROR: invalid sequenc\
e $id is not an RNA sequence. read seq is: $seq\\n\
\";\n        exit\n      }\n   \n    tcoffeelib_he\
ader_addseq($id, uc($seq), *FD_LIB);\n  }\ntcoffee\
lib_comment(\"generated by $myname on \" . localti\
me(), *FD_LIB);\n\n\n\n$i=0;\nforeach $id (keys (%\
seq))\n  {\n    my ($cleanid, $seq, $bpm);\n    $s\
eq=$seq{$id}{seq};\n    $cleanid = $id;\n    $clea\
nid =~ s,[/ ],_,g;# needed for rnaplfold\n    $seq\
 = rnaseqfmt($seq);\n    \n    $bpm = plfold($clea\
nid, rnaseqfmt($seq), $PROBTRESH);       \n    \n \
   tcoffeelib_struct($i+1, length($seq), $bpm, *FD\
_LIB);\n    $i++;\n}\n\n\ntcoffeelib_footer(*FD_LI\
B);\nclose FD_LIB;\nexit (0);\n\n","\n\n\n\n\n$cmd\
=join ' ', @ARGV;\nif ($cmd=~/-infile=(\\S+)/){ $s\
eqfile=$1;}\nif ($cmd=~/-outfile=(\\S+)/){ $libfil\
e=$1;}\n\n\n\n%s=read_fasta_seq ($seqfile);\n\nope\
n (F, \">$libfile\");\nforeach $name (keys (%s))\n\
  {\n    my $tclib=\"$name.RNAplfold_tclib\";\n   \
 print (F \">$name _F_ $tclib\\n\");\n    seq2RNAp\
lfold2tclib ($name, $s{$name}{seq}, $tclib);\n  }\\
nclose (F);\nexit (EXIT_SUCCESS);\n\nsub seq2RNApl\
fold2tclib\n  {\n    my ($name, $seq, $tclib)=@_;\\
n    my ($tmp);\n    $n++;\n    $tmp=\"tmp4seq2RNA\
plfold_tclib.$$.$n.pep\";\n    open (RF, \">$tmp\"\
);\n    print (RF \">$name\\n$seq\\n\");\n    clos\
e (RF);\n    \n    system \"t_coffee -other_pg RNA\
plfold2tclib.pl -in=$tmp -out=$tclib\";\n    \n   \
 unlink ($tmp);\n    return $tclib;\n  }\n    \n  \
  \nsub read_fasta_seq \n  {\n    my $f=@_[0];\n  \
  my %hseq;\n    my (@seq, @com, @name);\n    my (\
$a, $s,$nseq);\n\n    open (F, $f);\n    while (<F\
>)\n      {\n	$s.=$_;\n      }\n    close (F);\n\n\
    \n    @name=($s=~/>(\\S*).*\\n[^>]*/g);\n    \\
n    @seq =($s=~/>.*.*\\n([^>]*)/g);\n    @com =($\
s=~/>\\S*(.*)\\n([^>]*)/g);\n\n    \n    $nseq=$#n\
ame+1;\n    \n    for ($a=0; $a<$nseq; $a++)\n    \
  {\n	my $n=$name[$a];\n	$hseq{$n}{name}=$n;\n	$hs\
eq{$n}{seq}=$seq[$a];\n	$hseq{$n}{com}=$com[$a];\n\
      }\n    return %hseq;\n  }\n","use Getopt::Lo\
ng;\nuse File::Path;\nuse Env;\nuse FileHandle;\nu\
se Cwd;\nuse Sys::Hostname;\nour $PIDCHILD;\nour $\
ERROR_DONE;\nour @TMPFILE_LIST;\nour $EXIT_FAILURE\
=1;\nour $EXIT_SUCCESS=0;\n\nour $REFDIR=getcwd;\n\
our $EXIT_SUCCESS=0;\nour $EXIT_FAILURE=1;\n\nour \
$PROGRAM=\"tc_generic_method.pl\";\nour $CL=$PROGR\
AM;\n\nour $CLEAN_EXIT_STARTED;\nour $debug_lock=$\
ENV{\"DEBUG_LOCK\"};\nour $LOCKDIR=$ENV{\"LOCKDIR_\
4_TCOFFEE\"};\nif (!$LOCKDIR){$LOCKDIR=getcwd();}\\
nour $ERRORDIR=$ENV{\"ERRORDIR_4_TCOFFEE\"};\nour \
$ERRORFILE=$ENV{\"ERRORFILE_4_TCOFFEE\"};\n&set_lo\
ck ($$);\nif (isshellpid(getppid())){lock4tc(getpp\
id(), \"LLOCK\", \"LSET\", \"$$\\n\");}\n      \no\
ur $print;\nmy ($fmsq1, $fmsq2, $output, $outfile,\
 $arch, $psv, $hmmtop_home, $trim, $cov, $sample, \
$mode, $gor_home, $gor_seq, $gor_obs);\n\nGetOptio\
ns(\"-in=s\" => \\$fmsq1,\"-output=s\" =>\\$output\
 ,\"-out=s\" => \\$outfile, \"-arch=s\" => \\$arch\
,\"-psv=s\" => \\$psv, \"-hmmtop_home=s\", \\$hmmt\
op_home,\"-trim=s\" =>\\$trim ,\"-print=s\" =>\\$p\
rint,\"-cov=s\" =>\\$cov , \"-sample=s\" =>\\$samp\
le, \"-mode=s\" =>\\$mode, \"-gor_home=s\"=>\\$gor\
_home, \"-gor_seq=s\"=>\\$gor_seq,\"-gor_obs=s\"=>\
\\$gor_obs);\n\n\nif (!$mode){$mode = \"hmmtop\"}\\
nelsif ($mode eq \"hmmtop\"){;}\nelsif ($mode eq \\
"gor\"){;}\nelse {myexit(flush_error (\"-mode=$mod\
e is unknown\"));}\n\n\nour $HOME=$ENV{\"HOME\"};\\
nour $MCOFFEE=($ENV{\"MCOFFEE_4_TCOFFEE\"})?$ENV{\\
"MCOFFEE_4_TCOFFEE\"}:\"$HOME/.t_coffee/mcoffee\";\
\n\nif ($mode eq \"hmmtop\")\n  {\n    check_confi\
guration (\"hmmtop\");\n    if (-e $arch){$ENV{'HM\
MTOP_ARCH'}=$arch;}\n    elsif (-e $ENV{HMMTOP_ARC\
H}){$arch=$ENV{HMMTOP_ARCH};}\n    elsif (-e \"$MC\
OFFEE/hmmtop.arch\"){$arch=$ENV{'HMMTOP_ARCH'}=\"$\
MCOFFEE/hmmtop.arch\";}\n    elsif (-e \"$hmmtop_h\
ome/hmmtop.arc\"){$arch=$ENV{'HMMTOP_ARCH'}=\"$hmm\
top_home/hmmtop.arc\";}\n    else {myexit(flush_er\
ror ( \"Could not find ARCH file for hmmtop\"));}\\
n    \n    \n    if (-e $psv){$ENV{'HMMTOP_PSV'}=$\
psv;}\n    elsif (-e $ENV{HMMTOP_PSV}){$psv=$ENV{H\
MMTOP_PSV};}\n    elsif (-e \"$MCOFFEE/hmmtop.psv\\
"){$psv=$ENV{'HMMTOP_PSV'}=\"$MCOFFEE/hmmtop.psv\"\
;}\n    elsif (-e \"$hmmtop_home/hmmtop.psv\"){$ps\
v=$ENV{'HMMTOP_PSV'}=\"$hmmtop_home/hmmtop.psv\";}\
\n    else {myexit(flush_error ( \"Could not find \
PSV file for hmmtop\"));}\n  }\nelsif ($mode eq \"\
gor\")\n  {\n    our $GOR_SEQ;\n    our $GOR_OBS;\\
n    \n    check_configuration (\"gorIV\");\n    i\
f (-e $gor_seq){$GOR_SEQ=$gor_seq;}\n    elsif (-e\
 $ENV{GOR_SEQ}){$GOR_SEQ=$ENV{GOR_SEQ};}\n    elsi\
f (-e \"$MCOFFEE/New_KS.267.seq\"){$GOR_SEQ=\"$MCO\
FFEE/New_KS.267.seq\";}\n    elsif (-e \"$gor_home\
/New_KS.267.seq\"){$GOR_SEQ=\"$gor_home/New_KS.267\
.seq\";}\n    else {myexit(flush_error ( \"Could n\
ot find SEQ file for gor\"));}\n\n    if (-e $gor_\
obs){$GOR_OBS=$gor_obs;}\n    elsif (-e $ENV{GOR_O\
BS}){$GOR_OBS=$ENV{GOR_OBS};}\n    elsif (-e \"$MC\
OFFEE/New_KS.267.obs\"){$GOR_OBS=\"$MCOFFEE/New_KS\
.267.obs\";}\n    elsif (-e \"$gor_home/New_KS.267\
.obs\"){$GOR_OBS=\"$gor_home/New_KS.267.obs\";}\n \
   else {myexit(flush_error ( \"Could not find OBS\
 file for gor\"));}\n  }\n\n\nif ( ! -e $fmsq1){my\
exit(flush_error (\"Could Not Read Input file $fms\
q1\"));}\n\n\nmy $fmsq2=vtmpnam();\nmy $fmsq3=vtmp\
nam();\nmy $tmpfile=vtmpnam();\nmy $predfile=vtmpn\
am();\n\nif ($trim){$trim_action=\" +trim _aln_%%$\
trim\\_K1 \";}\nif ($cov) {$cov_action= \" +sim_fi\
lter _aln_c$cov \";}\n&safe_system(\"t_coffee -oth\
er_pg seq_reformat -in $fmsq1 -action +convert 'BO\
UJXZ-' $cov_action $trim_action -output fasta_aln \
-out $fmsq2\");\nmy (%pred, %seq, %predA);\n\n\n%s\
eq=read_fasta_seq($fmsq2);\n%seq=fasta2sample(\\%s\
eq, $sample);\n\nif (1==2 && $mode eq \"hmmtop\" &\
& $output eq \"cons\")\n  {\n    fasta2hmmtop_cons\
($outfile,\\%seq);\n  }\nelse\n  {\n    %pred=fast\
a2pred(\\%seq, $mode);\n    %predA=pred2aln (\\%pr\
ed, \\%seq);\n    \n    \n    if (!$output || $out\
put eq \"prediction\"){output_fasta_seq (\\%predA,\
 $outfile);}\n    elsif ($output eq \"color_html\"\
){pred2color (\\%pred,\\%seq, $outfile);}\n    els\
if ($output eq \"cons\"){pred2cons($outfile,\\%pre\
dA);}\n    else {flush_error (\"$output is an unkn\
own output mode\");}\n  }\n\nsub fasta2sample\n  {\
\n    my $SR=shift;\n    my $it=shift;\n    my %S=\
%$SR;\n    \n    my $seq=index2seq_name (\\%S, 1);\
\n    my $l=length($S{$seq}{seq});\n    my @sl=key\
s(%S);\n    my $nseq=$#sl+1;\n    my $index=$nseq;\
\n  \n    if (!$sample) {return %S;}\n    for (my \
$a=0; $a<$it; $a++)\n      {\n	my $newseq=\"\";\n	\
my $nname=\"$seq\\_sampled_$index\";\n	for (my $p=\
0; $p<$l; $p++)\n	  {\n	    my $i=int(rand($nseq))\
;\n	    \n	    my $name = $sl[$i];\n	    my $seq=$\
S{$name}{seq};\n	    my $r=substr ($seq, $p, 1);\n\
	    $newseq.=$r;\n	  }\n	$S{$nname}{name}=$nname;\
\n	$S{$nname}{seq}=$newseq;\n	$S{$nname}{com}=\"sa\
mpled\";\n	$S{$nname}{index}=++$index;\n      }\n \
   return %S;\n  }\n	      \nsub fasta2pred\n  {\n\
    my $s=shift;\n    my $mode=shift;\n\n    if ( \
$mode eq \"hmmtop\"){return fasta2hmmtop_pred($s);\
}\n    elsif ($mode eq \"gor\"){return fasta2gor_p\
red ($s);}\n  }\nsub fasta2hmmtop_cons\n  {\n    m\
y $outfile=shift;\n    my $SR=shift;\n    \n    my\
 $o = new FileHandle;\n    my $i = new FileHandle;\
\n    my $tmp_in =vtmpnam();\n    my $tmp_out=vtmp\
nam();\n    my %seq=%$SR;\n    my %pred;\n    my $\
N=keys(%seq);\n    \n    output_fasta_seq (\\%seq,\
$tmp_in, \"seq\");\n    `hmmtop -pi=mpred -if=$tmp\
_in -sf=FAS -pl 2>/dev/null >$tmp_out`;\n    open \
($o, \">$outfile\");\n    open ($i, \"$tmp_out\");\
\n    while (<$i>)\n      {\n	my $l=$_;\n	if (($l=\
~/>HP\\:\\s+(\\d+)\\s+(.*)/)){my $line=\">$2 NSEQ:\
 $N\\n\";print $o \"$line\";}\n	elsif ( ($l=~/.*pr\
ed(.*)/))  {my $line=\"$1\\n\";print $o \"$line\";\
}\n      }\n    close ($o);\n    close ($i);\n    \
return read_fasta_seq($tmp);\n  }\nsub fasta2hmmto\
p_pred\n  {\n    my $SR=shift;\n    my $o = new Fi\
leHandle;\n    my $i = new FileHandle;\n    my $tm\
p    =vtmpnam();\n    my $tmp_in =vtmpnam();\n    \
my $tmp_out=vtmpnam();\n    my %seq=%$SR;\n    my \
%pred;\n    \n\n    output_fasta_seq (\\%seq,$tmp_\
in, \"seq\");\n    `hmmtop -if=$tmp_in -sf=FAS -pl\
 2>/dev/null >$tmp_out`;\n    open ($o, \">$tmp\")\
;\n    open ($i, \"$tmp_out\");\n    while (<$i>)\\
n      {\n	my $l=$_;\n	if (($l=~/>HP\\:\\s+(\\d+)\\
\s+(.*)/)){my $line=\">$2\\n\";print $o \"$line\";\
}\n	elsif ( ($l=~/.*pred(.*)/))  {my $line=\"$1\\n\
\";print $o \"$line\";}\n      }\n    close ($o);\\
n    close ($i);\n    return read_fasta_seq($tmp);\
\n  }\n    \n	\n	\n	    \n	\n	\n\n	\nsub fasta2gor\
_pred\n  {\n    my $SR=shift;\n    my $o = new Fil\
eHandle;\n    my $i = new FileHandle;\n    my $tmp\
    =vtmpnam();\n    my $tmp_in =vtmpnam();\n    m\
y $tmp_out=vtmpnam();\n    my %seq=%$SR;\n    my %\
pred;\n    \n\n    output_fasta_seq (\\%seq,$tmp_i\
n, \"seq\");\n    `gorIV -prd $tmp_in -seq $GOR_SE\
Q -obs $GOR_OBS >$tmp_out`;\n    open ($o, \">$tmp\
\");\n    open ($i, \"$tmp_out\");\n    while (<$i\
>)\n      {\n	my $l=$_;\n\n	\n	if ( $l=~/>/){print\
 $o \"$l\";}\n	elsif ( $l=~/Predicted Sec. Struct.\
/){$l=~s/Predicted Sec. Struct\\.//;print $o \"$l\\
";}\n      }\n    close ($o);\n    close ($i);\n  \
  return read_fasta_seq($tmp);\n  }\n			\n			     \
\nsub index2seq_name\n  {\n    \n    my $SR=shift;\
\n    my $index=shift;\n    \n    \n    my %S=%$SR\
;\n    \n    foreach my $s (%S)\n      {\n	if ( $S\
{$s}{index}==$index){return $s;}\n      }\n    ret\
urn \"\";\n  }\n\nsub pred2cons\n  {\n    my $outf\
ile=shift;\n    my $predR=shift;\n    my $seq=shif\
t;\n    my %P=%$predR;\n    my %C;\n    my ($s,@r,\
$nseq);\n    my $f= new FileHandle;\n\n    open ($\
f, \">$outfile\");\n\n    if (!$seq){$seq=index2se\
q_name(\\%P,1);}\n    foreach my $s (keys(%P))\n  \
    {\n	$nseq++;\n	$string= $P{$s}{seq};\n	$string\
 = uc $string;\n	my @r=split (//,$string);\n	for (\
my $a=0; $a<=$#r; $a++)\n	  {\n	    if (($r[$a]=~/\
[OHICE]/)){$C{$a}{$r[$a]}++;}\n	  }\n      }\n    \
@l=keys(%C);\n    \n    \n    $s=$P{$seq}{seq};\n \
   print $f \">$seq pred based on $nseq\\n\";\n   \
 @r=split (//,$s);\n    \n    for (my $x=0; $x<=$#\
r; $x++)\n      {\n	if ($r[$x] ne \"-\")\n	  {\n	 \
   my $h=$C{$x}{H};\n	    my $i=$C{$x}{I};\n	    m\
y $o=$C{$x}{O};\n	    my $c=$C{$x}{C};\n	    my $e\
=$C{$x}{E};\n	    my $l=$i+$o;\n	    \n	    if ($h\
>=$i && $h>=$o && $h>=$c && $h>=$e){$r[$x]='H';}\n\
	    elsif ($i>=$o && $i>=$c && $i>=$e){$r[$x]='I'\
;}\n	    elsif ($o>=$c && $o>=$e){$r[$x]='O';}\n	 \
   elsif ($c>=$e){$r[$x]='C';}\n	    else {$r[$x]=\
'E';}\n	  }\n      }\n    $j=join ('', @r);\n    p\
rint $f \"$j\\n\";\n    close ($f);\n    return $j\
;\n  }\n\nsub pred2aln\n  {\n    my $PR=shift;\n  \
  my $AR=shift;\n    \n    my $f=new FileHandle;\n\
    my %P=%$PR;\n    my %A=%$AR;\n    my %PA;\n   \
 my $tmp=vtmpnam();\n    my $f= new FileHandle;\n \
   \n    open ($f, \">$tmp\");\n    foreach my $s \
(sort{$A{$a}{index}<=>$A{$b}{index}}(keys (%A)))\n\
      {\n	my (@list, $seq, @plist, @pseq, $L, $PL,\
 $c, $w);\n	my $seq;\n	my $seq=$A{$s}{seq};\n	my $\
pred=$P{$s}{seq};\n	$seq=pred2alnS($P{$s}{seq},$A{\
$s}{seq});\n	print $f \">$s\\n$seq\\n\";\n      }\\
n    close ($f);\n    return read_fasta_seq ($tmp)\
;\n  }\nsub pred2alnS\n  {\n    my $pred=shift;\n \
   my $aln= shift;\n    my ($j,$a,$b);\n    my @P=\
split (//, $pred);\n    my @A=split (//, $aln);\n \
   for ($a=$b=0;$a<=$#A; $a++)\n      {\n	if ($A[$\
a] ne \"-\"){$A[$a]=$P[$b++];}\n      }\n    if ($\
b!= ($#P+1)){add_warning (\"Could not thread seque\
nce: $b $#P\");}\n    \n    $j= join ('', @A);\n  \
  return $j;\n  }\nsub pred2color\n  {\n    my $pr\
edP=shift;\n    my $alnP=shift;\n    my $out=shift\
;\n    my $F=new FileHandle;\n    my $struc=vtmpna\
m();\n    my $aln=vtmpnam();\n    \n\n    output_f\
asta_seq ($alnP, $aln);\n    my %p=%$predP;\n    \\
n    open ($F, \">$struc\");\n    \n    \n    fore\
ach my $s (keys(%p))\n      {\n	\n	print $F \">$s\\
\n\";\n	my $s=uc($p{$s}{seq});\n	\n	$s=~s/[Oo]/0/g\
;\n	$s=~s/[Ee]/0/g;\n	\n	$s=~s/[Ii]/5/g;\n	$s=~s/[\
Cc]/5/g;\n	\n	$s=~s/[Hh]/9/g;\n	\n	print $F \"$s\\\
n\";\n      }\n    close ($F);\n    \n    \n    \n\
    safe_system ( \"t_coffee -other_pg seq_reforma\
t -in $aln -struc_in $struc -struc_in_f number_fas\
ta -output color_html -out $out\");\n    return;\n\
  }\n	  \n    \nsub display_fasta_seq\n  {\n    my\
 $SR=shift;\n    my %S=%$SR;\n    \n    foreach my\
 $s (sort{$S{$a}{index}<=>$S{$b}{index}}(keys (%S)\
))\n      {\n	print STDERR \">$s\\n$S{$s}{seq}\\n\\
";\n      }\n    close ($f);\n  }\nsub output_fast\
a_seq\n  {\n    my $SR=shift;\n    my $outfile=shi\
ft;\n    my $mode =shift;\n    my $f= new FileHand\
le;\n    my %S=%$SR;\n    \n    \n    open ($f, \"\
>$outfile\");\n    foreach my $s (sort{$S{$a}{inde\
x}<=>$S{$b}{index}}(keys (%S)))\n      {\n	my $seq\
=$S{$s}{seq};\n	if ( $mode eq \"seq\"){$seq=~s/\\-\
//g;}\n	print $f \">$s\\n$seq\\n\";\n      }\n    \
close ($f);\n  }\n      \nsub read_fasta_seq \n  {\
\n    my $f=$_[0];\n    my %hseq;\n    my (@seq, @\
com, @name);\n    my ($a, $s,$nseq);\n    my $inde\
x;\n    open (F, $f);\n    while (<F>)\n      {\n	\
$s.=$_;\n      }\n    close (F);\n\n    \n    @nam\
e=($s=~/>(\\S*).*\\n[^>]*/g);\n    \n    @seq =($s\
=~/>.*.*\\n([^>]*)/g);\n    @com =($s=~/>.*(.*)\\n\
([^>]*)/g);\n\n\n    $nseq=$#name+1;\n    \n  \n  \
  for ($a=0; $a<$nseq; $a++)\n      {\n	my $n=$nam\
e[$a];\n	my $s;\n	$hseq{$n}{name}=$n;\n	$s=$seq[$a\
];$s=~s/\\s//g;\n	$hseq{$n}{index}=++$index;\n	$hs\
eq{$n}{seq}=$s;\n	$hseq{$n}{com}=$com[$a];\n      \
}\n    return %hseq;\n  }\n\n\nsub file2head\n    \
  {\n	my $file = shift;\n	my $size = shift;\n	my $\
f= new FileHandle;\n	my $line;\n	open ($f,$file);\\
n	read ($f,$line, $size);\n	close ($f);\n	return $\
line;\n      }\nsub file2tail\n      {\n	my $file \
= shift;\n	my $size = shift;\n	my $f= new FileHand\
le;\n	my $line;\n	\n	open ($f,$file);\n	seek ($f,$\
size*-1, 2);\n	read ($f,$line, $size);\n	close ($f\
);\n	return $line;\n      }\n\n\nsub vtmpnam\n    \
  {\n	my $r=rand(100000);\n	my $f=\"file.$r.$$\";\\
n	while (-e $f)\n	  {\n	    $f=vtmpnam();\n	  }\n	\
push (@TMPFILE_LIST, $f);\n	return $f;\n      }\n\\
nsub myexit\n  {\n    my $code=@_[0];\n    if ($CL\
EAN_EXIT_STARTED==1){return;}\n    else {$CLEAN_EX\
IT_STARTED=1;}\n    ### ONLY BARE EXIT\n    exit (\
$code);\n  }\nsub set_error_lock\n    {\n      my \
$name = shift;\n      my $pid=$$;\n\n      \n     \
 &lock4tc ($$,\"LERROR\", \"LSET\", \"$$ -- ERROR:\
 $name $PROGRAM\\n\");\n      return;\n    }\nsub \
set_lock\n  {\n    my $pid=shift;\n    my $msg= sh\
ift;\n    my $p=getppid();\n    &lock4tc ($pid,\"L\
LOCK\",\"LRESET\",\"$p$msg\\n\");\n  }\nsub unset_\
lock\n   {\n     \n    my $pid=shift;\n    &lock4t\
c ($pid,\"LLOCK\",\"LRELEASE\",\"\");\n  }\nsub sh\
ift_lock\n  {\n    my $from=shift;\n    my $to=shi\
ft;\n    my $from_type=shift;\n    my $to_type=shi\
ft;\n    my $action=shift;\n    my $msg;\n    \n  \
  if (!&lock4tc($from, $from_type, \"LCHECK\", \"\\
")){return 0;}\n    $msg=&lock4tc ($from, $from_ty\
pe, \"LREAD\", \"\");\n    &lock4tc ($from, $from_\
type,\"LRELEASE\", $msg);\n    &lock4tc ($to, $to_\
type, $action, $msg);\n    return;\n  }\nsub isshe\
llpid\n  {\n    my $p=shift;\n    if (!lock4tc ($p\
, \"LLOCK\", \"LCHECK\")){return 0;}\n    else\n  \
    {\n	my $c=lock4tc($p, \"LLOCK\", \"LREAD\");\n\
	if ( $c=~/-SHELL-/){return 1;}\n      }\n    retu\
rn 0;\n  }\nsub isrootpid\n  {\n    if(lock4tc (ge\
tppid(), \"LLOCK\", \"LCHECK\")){return 0;}\n    e\
lse {return 1;}\n  }\nsub lock4tc\n	{\n	  my ($pid\
,$type,$action,$value)=@_;\n	  my $fname;\n	  my $\
host=hostname;\n	  \n	  if ($type eq \"LLOCK\"){$f\
name=\"$LOCKDIR/.$pid.$host.lock4tcoffee\";}\n	  e\
lsif ( $type eq \"LERROR\"){ $fname=\"$LOCKDIR/.$p\
id.$host.error4tcoffee\";}\n	  elsif ( $type eq \"\
LWARNING\"){ $fname=\"$LOCKDIR/.$pid.$host.warning\
4tcoffee\";}\n	  \n	  if ($debug_lock)\n	    {\n	 \
     print STDERR \"\\n\\t---lock4tc(tcg): $action\
 => $fname =>$value (RD: $LOCKDIR)\\n\";\n	    }\n\
\n	  if    ($action eq \"LCHECK\") {return -e $fna\
me;}\n	  elsif ($action eq \"LREAD\"){return file2\
string($fname);}\n	  elsif ($action eq \"LSET\") {\
return string2file ($value, $fname, \">>\");}\n	  \
elsif ($action eq \"LRESET\") {return string2file \
($value, $fname, \">\");}\n	  elsif ($action eq \"\
LRELEASE\") \n	    {\n	      if ( $debug_lock)\n		\
{\n		  my $g=new FileHandle;\n		  open ($g, \">>$f\
name\");\n		  print $g \"\\nDestroyed by $$\\n\";\\
n		  close ($g);\n		  safe_system (\"mv $fname $fn\
ame.old\");\n		}\n	      else\n		{\n		  unlink ($f\
name);\n		}\n	    }\n	  return \"\";\n	}\n	\nsub f\
ile2string\n	{\n	  my $file=@_[0];\n	  my $f=new F\
ileHandle;\n	  my $r;\n	  open ($f, \"$file\");\n	\
  while (<$f>){$r.=$_;}\n	  close ($f);\n	  return\
 $r;\n	}\nsub string2file \n    {\n    my ($s,$fil\
e,$mode)=@_;\n    my $f=new FileHandle;\n    \n   \
 open ($f, \"$mode$file\");\n    print $f  \"$s\";\
\n    close ($f);\n  }\n\nBEGIN\n    {\n      sran\
d;\n    \n      $SIG{'SIGUP'}='signal_cleanup';\n \
     $SIG{'SIGINT'}='signal_cleanup';\n      $SIG{\
'SIGQUIT'}='signal_cleanup';\n      $SIG{'SIGILL'}\
='signal_cleanup';\n      $SIG{'SIGTRAP'}='signal_\
cleanup';\n      $SIG{'SIGABRT'}='signal_cleanup';\
\n      $SIG{'SIGEMT'}='signal_cleanup';\n      $S\
IG{'SIGFPE'}='signal_cleanup';\n      \n      $SIG\
{'SIGKILL'}='signal_cleanup';\n      $SIG{'SIGPIPE\
'}='signal_cleanup';\n      $SIG{'SIGSTOP'}='signa\
l_cleanup';\n      $SIG{'SIGTTIN'}='signal_cleanup\
';\n      $SIG{'SIGXFSZ'}='signal_cleanup';\n     \
 $SIG{'SIGINFO'}='signal_cleanup';\n      \n      \
$SIG{'SIGBUS'}='signal_cleanup';\n      $SIG{'SIGA\
LRM'}='signal_cleanup';\n      $SIG{'SIGTSTP'}='si\
gnal_cleanup';\n      $SIG{'SIGTTOU'}='signal_clea\
nup';\n      $SIG{'SIGVTALRM'}='signal_cleanup';\n\
      $SIG{'SIGUSR1'}='signal_cleanup';\n\n\n     \
 $SIG{'SIGSEGV'}='signal_cleanup';\n      $SIG{'SI\
GTERM'}='signal_cleanup';\n      $SIG{'SIGCONT'}='\
signal_cleanup';\n      $SIG{'SIGIO'}='signal_clea\
nup';\n      $SIG{'SIGPROF'}='signal_cleanup';\n  \
    $SIG{'SIGUSR2'}='signal_cleanup';\n\n      $SI\
G{'SIGSYS'}='signal_cleanup';\n      $SIG{'SIGURG'\
}='signal_cleanup';\n      $SIG{'SIGCHLD'}='signal\
_cleanup';\n      $SIG{'SIGXCPU'}='signal_cleanup'\
;\n      $SIG{'SIGWINCH'}='signal_cleanup';\n     \
 \n      $SIG{'INT'}='signal_cleanup';\n      $SIG\
{'TERM'}='signal_cleanup';\n      $SIG{'KILL'}='si\
gnal_cleanup';\n      $SIG{'QUIT'}='signal_cleanup\
';\n      \n      our $debug_lock=$ENV{\"DEBUG_LOC\
K\"};\n      \n      \n      \n      \n      forea\
ch my $a (@ARGV){$CL.=\" $a\";}\n      if ( $debug\
_lock ){print STDERR \"\\n\\n\\n********** START P\
G: $PROGRAM *************\\n\";}\n      if ( $debu\
g_lock ){print STDERR \"\\n\\n\\n**********(tcg) L\
OCKDIR: $LOCKDIR $$ *************\\n\";}\n      if\
 ( $debug_lock ){print STDERR \"\\n --- $$ -- $CL\\
\n\";}\n      \n	     \n      \n      \n    }\nsub\
 flush_error\n  {\n    my $msg=shift;\n    return \
add_error ($EXIT_FAILURE,$$, $$,getppid(), $msg, $\
CL);\n  }\nsub add_error \n  {\n    my $code=shift\
;\n    my $rpid=shift;\n    my $pid=shift;\n    my\
 $ppid=shift;\n    my $type=shift;\n    my $com=sh\
ift;\n    \n    $ERROR_DONE=1;\n    lock4tc ($rpid\
, \"LERROR\",\"LSET\",\"$pid -- ERROR: $type\\n\")\
;\n    lock4tc ($$, \"LERROR\",\"LSET\", \"$pid --\
 COM: $com\\n\");\n    lock4tc ($$, \"LERROR\",\"L\
SET\", \"$pid -- STACK: $ppid -> $pid\\n\");\n   \\
n    return $code;\n  }\nsub add_warning \n  {\n  \
  my $rpid=shift;\n    my $pid =shift;\n    my $co\
mmand=shift;\n    my $msg=\"$$ -- WARNING: $comman\
d\\n\";\n    print STDERR \"$msg\";\n    lock4tc (\
$$, \"LWARNING\", \"LSET\", $msg);\n  }\n\nsub sig\
nal_cleanup\n  {\n    print dtderr \"\\n**** $$ (t\
cg) was killed\\n\";\n    &cleanup;\n    exit ($EX\
IT_FAILURE);\n  }\nsub clean_dir\n  {\n    my $dir\
=@_[0];\n    if ( !-d $dir){return ;}\n    elsif (\
!($dir=~/tmp/)){return ;}#safety check 1\n    elsi\
f (($dir=~/\\*/)){return ;}#safety check 2\n    el\
se\n      {\n	`rm -rf $dir`;\n      }\n    return;\
\n  }\nsub cleanup\n  {\n    #print stderr \"\\n--\
--tc: $$ Kills $PIDCHILD\\n\";\n    #kill (SIGTERM\
,$PIDCHILD);\n    my $p=getppid();\n    $CLEAN_EXI\
T_STARTED=1;\n    \n    \n    \n    if (&lock4tc($\
$,\"LERROR\", \"LCHECK\", \"\"))\n      {\n	my $pp\
id=getppid();\n	if (!$ERROR_DONE) \n	  {\n	    &lo\
ck4tc($$,\"LERROR\", \"LSET\", \"$$ -- STACK: $p -\
> $$\\n\");\n	    &lock4tc($$,\"LERROR\", \"LSET\"\
, \"$$ -- COM: $CL\\n\");\n	  }\n      }\n    my $\
warning=&lock4tc($$, \"LWARNING\", \"LREAD\", \"\"\
);\n    my $error=&lock4tc($$,  \"LERROR\", \"LREA\
D\", \"\");\n    #release error and warning lock i\
f root\n    \n    if (isrootpid() && ($warning || \
$error) )\n      {\n	\n	print STDERR \"***********\
***** Summary *************\\n$error\\n$warning\\n\
\";\n\n	&lock4tc($$,\"LERROR\",\"RELEASE\",\"\");\\
n	&lock4tc($$,\"LWARNING\",\"RELEASE\",\"\");\n   \
   } \n    \n    \n    foreach my $f (@TMPFILE_LIS\
T)\n      {\n	if (-e $f){unlink ($f);} \n      }\n\
    foreach my $d (@TMPDIR_LIST)\n      {\n	clean_\
dir ($d);\n      }\n    #No More Lock Release\n   \
 #&lock4tc($$,\"LLOCK\",\"LRELEASE\",\"\"); #relea\
se lock \n\n    if ( $debug_lock ){print STDERR \"\
\\n\\n\\n********** END PG: $PROGRAM ($$) ********\
*****\\n\";}\n    if ( $debug_lock ){print STDERR \
\"\\n\\n\\n**********(tcg) LOCKDIR: $LOCKDIR $$ **\
***********\\n\";}\n  }\nEND \n  {\n    \n    &cle\
anup();\n  }\n   \n\nsub safe_system \n{\n  my $co\
m=shift;\n  my $ntry=shift;\n  my $ctry=shift;\n  \
my $pid;\n  my $status;\n  my $ppid=getppid();\n  \
if ($com eq \"\"){return 1;}\n  \n  \n\n  if (($pi\
d = fork ()) < 0){return (-1);}\n  if ($pid == 0)\\
n    {\n      set_lock($$, \" -SHELL- $com (tcg)\"\
);\n      exec ($com);\n    }\n  else\n    {\n    \
  lock4tc ($$, \"LLOCK\", \"LSET\", \"$pid\\n\");#\
update parent\n      $PIDCHILD=$pid;\n    }\n  if \
($debug_lock){printf STDERR \"\\n\\t .... safe_sys\
tem (fasta_seq2hmm)  p: $$ c: $pid COM: $com\\n\";\
}\n\n  waitpid ($pid,WTERMSIG);\n\n  shift_lock ($\
pid,$$, \"LWARNING\",\"LWARNING\", \"LSET\");\n\n \
 if ($? == $EXIT_FAILURE || lock4tc($pid, \"LERROR\
\", \"LCHECK\", \"\"))\n    {\n      if ($ntry && \
$ctry <$ntry)\n	{\n	  add_warning ($$,$$,\"$com fa\
iled [retry: $ctry]\");\n	  lock4tc ($pid, \"LRELE\
ASE\", \"LERROR\", \"\");\n	  return safe_system (\
$com, $ntry, ++$ctry);\n	}\n      elsif ($ntry == \
-1)\n	{\n	  if (!shift_lock ($pid, $$, \"LERROR\",\
 \"LWARNING\", \"LSET\"))\n	    {\n	      add_warn\
ing ($$,$$,\"$com failed\");\n	    }\n	  else\n	  \
  {\n	      lock4tc ($pid, \"LRELEASE\", \"LERROR\\
", \"\");\n	    }\n	  return $?;}\n      else\n	{\\
n	  if (!shift_lock ($pid,$$, \"LERROR\",\"LERROR\\
", \"LSET\"))\n	    {\n	      myexit(add_error ($E\
XIT_FAILURE,$$,$pid,getppid(), \"UNSPECIFIED syste\
m\", $com));\n	    }\n	}\n    }\n  return $?;\n}\n\
\nsub check_configuration \n    {\n      my @l=@_;\
\n      my $v;\n      foreach my $p (@l)\n	{\n	  \\
n	  if   ( $p eq \"EMAIL\")\n	    { \n	      if ( \
!($EMAIL=~/@/))\n		{\n		add_warning($$,$$,\"Could \
Not Use EMAIL\");\n		myexit(add_error ($EXIT_FAILU\
RE,$$,$$,getppid(),\"EMAIL\",\"$CL\"));\n	      }\\
n	    }\n	  elsif( $p eq \"INTERNET\")\n	    {\n	 \
     if ( !&check_internet_connection())\n		{\n		 \
 myexit(add_error ($EXIT_FAILURE,$$,$$,getppid(),\\
"INTERNET\",\"$CL\"));\n		}\n	    }\n	  elsif( $p \
eq \"wget\")\n	    {\n	      if (!&pg_is_installed\
 (\"wget\") && !&pg_is_installed (\"curl\"))\n		{\\
n		  myexit(add_error ($EXIT_FAILURE,$$,$$,getppid\
(),\"PG_NOT_INSTALLED:wget\",\"$CL\"));\n		}\n	   \
 }\n	  elsif( !(&pg_is_installed ($p)))\n	    {\n	\
      myexit(add_error ($EXIT_FAILURE,$$,$$,getppi\
d(),\"PG_NOT_INSTALLED:$p\",\"$CL\"));\n	    }\n	}\
\n      return 1;\n    }\nsub pg_is_installed\n  {\
\n    my @ml=@_;\n    my $r, $p, $m;\n    my $supp\
orted=0;\n    \n    my $p=shift (@ml);\n    if ($p\
=~/::/)\n      {\n	if (safe_system (\"perl -M$p -e\
 1\")==$EXIT_SUCCESS){return 1;}\n	else {return 0;\
}\n      }\n    else\n      {\n	$r=`which $p 2>/de\
v/null`;\n	if ($r eq \"\"){return 0;}\n	else {retu\
rn 1;}\n      }\n  }\n\n\n\nsub check_internet_con\
nection\n  {\n    my $internet;\n    my $tmp;\n   \
 &check_configuration ( \"wget\"); \n    \n    $tm\
p=&vtmpnam ();\n    \n    if     (&pg_is_installed\
    (\"wget\")){`wget www.google.com -O$tmp >/dev/\
null 2>/dev/null`;}\n    elsif  (&pg_is_installed \
   (\"curl\")){`curl www.google.com -o$tmp >/dev/n\
ull 2>/dev/null`;}\n    \n    if ( !-e $tmp || -s \
$tmp < 10){$internet=0;}\n    else {$internet=1;}\\
n    if (-e $tmp){unlink $tmp;}\n\n    return $int\
ernet;\n  }\nsub check_pg_is_installed\n  {\n    m\
y @ml=@_;\n    my $r=&pg_is_installed (@ml);\n    \
if (!$r && $p=~/::/)\n      {\n	print STDERR \"\\n\
You Must Install the perl package $p on your syste\
m.\\nRUN:\\n\\tsudo perl -MCPAN -e 'install $pg'\\\
n\";\n      }\n    elsif (!$r)\n      {\n	myexit(f\
lush_error(\"\\nProgram $p Supported but Not Insta\
lled on your system\"));\n      }\n    else\n     \
 {\n	return 1;\n      }\n  }\n\n\n\n","\n\n\n\n\nm\
y $FMODEL =\"\"; \nmy $TMPDIR = \"/tmp\";\n\n\n\n\\
nmy $NUCALPH = \"ACGTUNRYMKSWHBVD\";\nmy $PRIMNUCA\
LPH = \"ACGTUN\";\nuse vars qw($NUCALPH $PRIMNUCAL\
PH $TMPDIR);\n\n\nmy $errmsg;\nuse vars qw($errmsg\
);\n\n\n\nuse Getopt::Long;\nuse Cwd;\nuse File::B\
asename;\nuse File::Temp qw/ tempfile tempdir /;\n\
use File::Copy;\nuse File::Path;\n\n\n\nsub usage(\
;$)\n{\n    my ($errmsg) = @_;\n    my $myname = b\
asename($0);\n\n    if ($errmsg) {\n        print \
STDERR \"ERROR: $errmsg\\n\";\n    }\n\n    print \
STDERR << \"EOF\";\n    \n$myname: align two seque\
nces by means of consan\\'s sfold\nUsage:\n $mynam\
e -i file -o file -d path\nOptions:\n -i|--in : pa\
irwise input sequence file\n -o|--out: output alig\
nment\n -d|--directory containing data\n\nEOF\n}\n\
\nsub read_stk_aln \n  {\n    my $f=$_[0];\n    my\
 ($seq, $id);\n    \n    my %hseq;\n\n    open (ST\
K, \"$f\");\n    while (<STK>)\n      {\n	if ( /^#\
/ || /^\\/\\// || /^\\s*$/){;}\n	else\n	  {\n	    \
($id,$seq)=/(\\S+)\\s+(\\S+)/;\n	    $hseq{$id}{'s\
eq'}.=$seq;\n	  }\n      }\n    close (STK);\n    \
return %hseq;\n  }\nsub read_fasta_seq \n  {\n    \
my $f=$_[0];\n    my %hseq;\n    my (@seq, @com, @\
name);\n    my ($a, $s,$nseq);\n\n    open (F, $f)\
;\n    while (<F>)\n      {\n	$s.=$_;\n      }\n  \
  close (F);\n\n    \n    @name=($s=~/>(.*).*\\n[^\
>]*/g);\n    \n    @seq =($s=~/>.*.*\\n([^>]*)/g);\
\n    @com =($s=~/>.*(.*)\\n([^>]*)/g);\n\n    \n \
   $nseq=$#name+1;\n    \n    for ($a=0; $a<$nseq;\
 $a++)\n      {\n	my $n=$name[$a];\n	$hseq{$n}{nam\
e}=$n;\n	$hseq{$n}{seq}=$seq[$a];\n	$hseq{$n}{com}\
=$com[$a];\n      }\n    return %hseq;\n  }\n\n\n\\
nsub sfold_parseoutput($$)\n{\n    my ($frawout, $\
foutfa) = @_;\n    my %haln;\n    my ($fstk, $cmd,\
 $id);\n    open FOUTFA, \">$foutfa\";\n    \n    \
$fstk = $frawout . \".stk\";\n    \n    # first li\
ne of raw out contains info\n    # remaining stuff\
 is stockholm formatted\n    $cmd = \"sed -e '1d' \
$frawout\";\n    system(\"$cmd > $fstk\");\n    if\
 ($? != 0) {\n        $errmsg = \"command failed w\
ith exit status $?.\";\n        $errmsg .=  \"Comm\
and was \\\"$cmd\\\"\";\n        return -1;\n    }\
\n\n    # this gives an error message. just ignore\
 it...\n    %haln=read_stk_aln ( $fstk);\n    fore\
ach $i (keys (%haln))\n      {\n	my $s;\n	$s=$haln\
{$i}{'seq'};\n	$s =~ s/\\./-/g;\n	print FOUTFA \">\
$i\\n$s\\n\";\n      }\n    close FOUTFA;\n    ret\
urn 0;\n}\n\n\n\n\nsub sfold_wrapper($$$$)\n{\n   \
 \n    my ($fs1, $fs2, $fmodel, $foutfa) = @_;\n  \
  \n\n    my ($cmd, $frawout, $ferrlog, $freadme, \
$ftimelog, $fstk);\n\n    # add  basename($fmsqin)\
 (unknown here!)\n    $frawout = \"sfold.log\";\n \
   $ferrlog = \"sfold.err\";\n    $ftimelog = \"sf\
old.time\";\n    $freadme =  \"sfold.README\";\n  \
  $fstk = \"sfold.stk\";\n    \n    # prepare exec\
ution...\n    #\n    # ./tmp is essential for dswp\
align\n    # otherwise you'll get a segfault\n    \
mkdir \"./tmp\";\n    \n    $cmd = \"sfold -m $fmo\
del $fs1 $fs2\";\n    open(FREADME,\">$freadme\");\
\n    print FREADME \"$cmd\\n\"; \n    close(FREAD\
ME);\n\n    # and go\n    #\n    system(\"/usr/bin\
/time -p -o $ftimelog $cmd >$frawout 2>$ferrlog\")\
;\n    if ($? != 0) {\n        $errmsg = \"command\
 failed with exit status $?\";\n        $errmsg .=\
 \"command was \\\"$cmd\\\". See \" . getcwd . \"\\
\n\";\n        return -1;\n    }\n\n    return sfo\
ld_parseoutput($frawout, $foutfa);\n}\n\n\n\n\n\n\\
n\nmy ($help, $fmsqin, $fmsaout);\nGetOptions(\"he\
lp\"  => \\$help,\n           \"in=s\" => \\$fmsqi\
n,\n           \"out=s\" => \\$fmsaout,\n	   \"dat\
a=s\" => \\$ref_dir);\n\n\n\nif ($help) {\n    usa\
ge();\n    exit(0);\n}\nif (! defined($fmsqin)) {\\
n    usage('missing input filename');\n    exit(1)\
;\n}\nif (! defined($fmsaout)) {\n    usage('missi\
ng output filename');\n    exit(1);\n\n}\nif (scal\
ar(@ARGV)) {\n    usage('Unknown remaining args');\
\n    exit(1);\n}\n\n$FMODEL = \"$ref_dir/mix80.mo\
d\";\nif (! -e \"$FMODEL\") {\n    die(\"couldn't \
find sfold grammar model file. Expected $FMODEL\\n\
\");\n}\n\n\nmy %hseq=read_fasta_seq ($fmsqin);\nm\
y $id;\n\nforeach $id (keys(%hseq))\n  {\n    push\
(@seq_array, $hseq{$id});\n  }\n\nif ( scalar(@seq\
_array) != 2 ) {\n    die(\"Need *exactly* two seq\
uences as input (pairwise alignment!).\")\n}\n\n\n\
\nmy ($sec, $min, $hour, $mday, $mon, $year, $wday\
, $yday, $isdst) = localtime(time);\nmy $datei = s\
printf(\"%4d-%02d-%02d\", $year+1900, $mon+1, $mda\
y);\nmy $templ = basename($0) . \".\" . $datei . \\
".pid-\" . $$ . \".XXXXXX\";\nmy $wd = tempdir ( $\
templ, DIR => $TMPDIR);\n\ncopy($fmsqin, \"$wd/\" \
. basename($fmsqin) . \".org\"); # for reproductio\
n\ncopy($FMODEL, \"$wd\");\nmy $fmodel = basename(\
$FMODEL);\nmy $orgwd = getcwd;\nchdir $wd;\n\n\n\n\
my @sepseqfiles;\nforeach $id (keys(%hseq)) {\n   \
 my ($seq, $orgseq, $fname, $sout);\n    $seq=$hse\
q{$id}{'seq'};\n    \n    $fname = basename($fmsqi\
n) . \"_$id.fa\";\n    # replace funnies in file/i\
d name (e.g. \"/\" \" \" etc)\n    $fname =~ s,[/ \
],_,g;\n    open (PF, \">$fname\");\n    print (PF\
 \">$id\\n$seq\\n\");\n    close (PF);\n\n    push\
(@sepseqfiles, $fname);\n}\n\nmy ($f1, $f2, $fout)\
;\n$f1 = $sepseqfiles[0];\n$f2 = $sepseqfiles[1];\\
n$fout = $wd . basename($fmsqin) . \".out.fa\";\ni\
f (sfold_wrapper($f1, $f2, $fmodel, \"$fout\") != \
0) {\n    printf STDERR \"ERROR: See logs in $wd\\\
n\";\n    exit(1);\n} else {\n    chdir $orgwd;\n \
   copy($fout, $fmsaout);\n    rmtree($wd);\n   ex\
it(0);\n}\n","\nuse Env qw(HOST);\nuse Env qw(HOME\
);\nuse Env qw(USER);\n\n\n$tmp=clean_cr ($ARGV[0]\
);\nopen (F, $tmp);\n\nwhile ( <F>)\n  {\n    my $\
l=$_;\n    if ( $l=~/^# STOCKHOLM/){$stockholm=1;}\
\n    elsif ( $stockholm && $l=~/^#/)\n      {\n	$\
l=~/^#(\\S+)\\s+(\\S+)\\s+(\\S*)/g;\n	$l=\"_stockh\
olmhasch_$1\\_stockholmspace_$2 $3\\n\";\n      }\\
n    $file.=$l;\n  }\nclose (F);\nunlink($tmp);\n$\
file1=$file;\n\n$file=~s/\\#/_hash_symbol_/g;\n$fi\
le=~s/\\@/_arobase_symbol_/g;\n\n\n$file=~s/\\n[\\\
.:*\\s]+\\n/\\n\\n/g;\n\n$file=~s/\\n[ \\t\\r\\f]+\
(\\b)/\\n\\1/g;\n\n\n$file=~s/(\\n\\S+)(\\s+)(\\S)\
/\\1_blank_\\3/g;\n\n$file=~s/[ ]//g;\n$file=~s/_b\
lank_/ /g;\n\n\n\n$file =~s/\\n\\s*\\n/#/g;\n\n$fi\
le.=\"#\";\n$file =~s/\\n/@/g;\n\n\n\n\n@blocks=sp\
lit /\\#/, $file;\nshift (@blocks);\n@s=split /\\@\
/, $blocks[0];\n$nseq=$#s+1;\n\n\n\n$file=join '@'\
, @blocks;\n@lines=split /\\@/,$file;\n\n$c=0;\n\n\
foreach $l (@lines)\n  {\n    if (!($l=~/\\S/)){ne\
xt;}\n    elsif ($stockholm && ($l=~/^\\/\\// || $\
l=~/STOCKHOLM/)){next;}#get read of STOCHOLM Termi\
nator\n   \n    $l=~/(\\S+)\\s+(\\S*)/g;\n    $n=$\
1; $s=$2;\n    \n    $seq[$c].=$s;\n    $name[$c]=\
$n;\n    $c++;\n    \n    if ( $c==$nseq){$c=0;}\n\
    \n  } \n\nif ( $c!=0)\n      {\n	print STDERR \
\"ERROR: $ARGV[0] is NOT an MSA in Clustalw format\
: make sure there is no blank line within a block \
[ERROR]\\n\";\n	exit (EXIT_FAILURE);\n      }\n\nf\
or ($a=0; $a< $nseq; $a++)\n  {\n    $name[$a]=cle\
anstring ($name[$a]);\n    $seq[$a]=cleanstring ($\
seq[$a]);\n    $seq[$a]=breakstring($seq[$a], 60);\
\n    \n    $line=\">$name[$a]\\n$seq[$a]\\n\";\n \
   \n    print \"$line\";\n  }\nexit (EXIT_SUCCESS\
);\n\nsub cleanstring\n  {\n    my $s=@_[0];\n    \
$s=~s/_hash_symbol_/\\#/g;\n    $s=~s/_arobase_sym\
bol_/\\@/g;\n    $s=~s/[ \\t]//g;\n    return $s;\\
n  }\nsub breakstring\n  {\n    my $s=@_[0];\n    \
my $size=@_[1];\n    my @list;\n    my $n,$ns, $sy\
mbol;\n    \n    @list=split //,$s;\n    $n=0;$ns=\
\"\";\n    foreach $symbol (@list)\n      {\n	if (\
 $n==$size)\n	  {\n	    $ns.=\"\\n\";\n	    $n=0;\\
n	  }\n	$ns.=$symbol;\n	$n++;\n      }\n    return\
 $ns;\n    }\n\nsub clean_cr\n  {\n    my $f=@_[0]\
;\n    my $file;\n    \n    $tmp=\"f$.$$\";\n    \\
n    \n    open (IN, $f);\n    open (OUT, \">$tmp\\
");\n    \n    while ( <IN>)\n      {\n	$file=$_;\\
n	$file=~s/\\r\\n/\\n/g;\n	$file=~s/\\n\\r/\\n/g;\\
n	$file=~s/\\r\\r/\\n/g;\n	$file=~s/\\r/\\n/g;\n	p\
rint OUT \"$file\";\n      }\n    \n    close (IN)\
;\n    close (OUT);\n    return $tmp;\n  }\n","use\
 Env qw(HOST);\nuse Env qw(HOME);\nuse Env qw(USER\
);\n\n\n$query_start=-1;\n$query_end=-1;\n\nwhile \
(<>)\n  {\n    if ( /\\/\\//){$in_aln=1;}\n    els\
if ( $in_aln && /(\\S+)\\s+(.*)/)\n      {\n\n\n	$\
name=$1;\n	\n\n	$seq=$2;\n	$seq=~s/\\s//g;\n      \
  $seq=~s/\\~/\\-/g;\n	$seq=~s/\\./\\-/g;\n	if ( $\
list{$n}{'name'} && $list{$n}{'name'} ne $name)\n	\
  {\n	    print \"$list{$n}{'name'} Vs $name\";\n	\
    \n	    exit (EXIT_FAILURE);\n	  }\n	else\n	  {\
\n	    $list{$n}{'name'}= $name;\n	  }\n\n	$list{$\
n}{'seq'}=$list{$n}{'seq'}.$seq;\n	\n	$nseq=++$n;\\
n	\n      }\n    else\n      {$n=0;}\n  }\n\n\nfor\
 ($a=0; $a<$nseq; $a++)\n  {\n    print \">$list{$\
a}{'name'}\\n$list{$a}{'seq'}\\n\";\n  }\n      \n\
","\nuse Env qw(HOST);\nuse Env qw(HOME);\nuse Env\
 qw(USER);\n\n                                    \
                    \nuse strict;                 \
                            \nuse warnings;\nuse d\
iagnostics;\n\nmy $in_hit_list, my $in_aln=0, my(%\
name_list)=(),my (%list)=(),my $n_seq=0; my $test=\
0;\nmy($j)=0, my $n=0, my $nom, my $lg_query, my %\
vu=();\n\nopen (F, \">tmp\");\n\n$/=\"\\n\";\nwhil\
e (<>)\n{\n    print F $_;\n    if($_ =~ /Query=\\\
s*(.+?)\\s/i) { $nom=$1;}\n\n    if ( /Sequences p\
roducing significant alignments/){$in_hit_list=1;}\
\n    \n    if ($_=~ /^pdb\\|/i) { $_=~ s/pdb\\|//\
g; }\n    if ($_=~ /^(1_\\d+)\\s+\\d+/) { $_=~ s/$\
1/QUERY/;}\n      \n    if ( /^(\\S+).+?\\s+[\\d.]\
+\\s+([\\de.-]+)\\s+$/ && $in_hit_list)	\n    {\n	\
my($id)=$1; # \n	$id=~ s/\\|/_/g; #\n	if ($id =~ /\
.+_$/) { chop($id) }; #\n	$name_list{$n_seq++}=$id\
;\n	$name_list{$n_seq-1}=~ s/.*\\|//g;     \n    }\
\n  \n    if (/query/i) {$in_aln=1;}\n    if ( /^(\
\\S+)\\s+(\\d+)\\s+([a-zA-Z-]+)\\s+(\\d+)/ || /^(\\
\S+)(\\s+)(\\-+)(\\s+)/ && ($in_aln == 1))\n    {\\
n	my $name=$1;\n	my $start=$2;\n	my $seq=$3;\n	my \
$end=$4;\n		\n	if ($name =~ /QUERY/i) { $lg_query=\
length($seq); }\n\n	unless ($test > $n) #m\n	{\n	 \
   my(@seqq)= split('',$seq);\n	    my($gap_missin\
g)= scalar(@seqq);\n	    \n	    while ($gap_missin\
g != $lg_query)  { unshift (@seqq,\"-\"); $gap_mis\
sing= scalar(@seqq); }\n	    $seq=join('',@seqq); \
 #m\n	}\n	\n	if ($name =~ /QUERY/i)\n	{\n	    $n=0\
; %vu=(); $j=0;\n	    $list{$n}{'real_name'}=\"$no\
m\";\n	}	\n	else\n	{\n	    unless (exists $vu{$nam\
e}) { ++$j;}	\n	    $list{$n}{'real_name'}=$name_l\
ist{$j-1};\n	}\n		\n	$list{$n}{'name'}=$name;\n\n	\
$seq=~tr/a-z/A-Z/;\n	$list{$n}{'seq'}=$list{$n}{'s\
eq'};\n	$list{$n}{'seq'}.=$seq;\n\n	$n++;\n	$vu{$n\
ame}++;\n	$test++;\n   } \n    \n}\n\nmy @numero=(\
);\n\nfor (my $a=0; $a<$n; $a++) #m\n{\n    my $lo\
ng=length($list{0}{'seq'});  \n    my $long1= leng\
th($list{$a}{'seq'});\n  \n    while ($long1 ne $l\
ong)\n    {\n	$list{$a}{'seq'}.=\"-\";\n	$long1= l\
ength ($list{$a}{'seq'});\n    } \n \n    push (@n\
umero,\"$list{$a}{'name'} $list{$a}{'real_name'}\\\
n\");\n}\n\nmy %dejavu=();\n\n\nfor (my $i=0; $i<=\
$#numero; $i++)\n{\n    my $s=\">$list{$i}{'real_n\
ame'}\\n$list{$i}{'seq'}\\n\";\n    my $k=0;\n    \
\n    if (exists $dejavu{$numero[$i]}) {next;}\n  \
  else\n    {	\n	for ($j=0; $j<$n ; $j++)\n	{\n	  \
  if (\"$numero[$i]\" eq \"$numero[$j]\" && $j != \
$i )\n	    {\n		++$k;\n		$s .=\">$list{$j}{'real_n\
ame'}\\n$list{$j}{'seq'}\\n\";\n	    }\n	}	\n    }\
\n    \n    if ($k>0) \n    {\n	my $cons;\n	open (\
SOR,\">tempo_aln2cons\"); print SOR $s;  close SOR\
 ;\n	open (COM,\"t_coffee -other_pg seq_reformat -\
in tempo_aln2cons -action +aln2cons +upper |\") ; \
\n     	while (<COM>)\n	{	\n	    if (/^>/) { $cons\
 =\">$list{$i}{'real_name'}\\n\"; next;}\n	    $_=\
~ s/\\n//g;\n	    $cons .=$_;\n	}\n	close COM; unl\
ink (\"tempo_aln2cons\");\n	print $cons,\"\\n\"; p\
rint F $cons,\"\\n\";\n    }	\n    else  { print $\
s;  print F $s; }\n    \n    $dejavu{$numero[$i]}+\
+;\n} #m\n\nexit;\n\n\n\n\n\n\n\n\n\n\n\n","use En\
v;\n\n\n$tmp_dir=\"\";\n$init_dir=\"\";\n$program=\
\"tc_generic_method.pl\";\n\n$blast=@ARGV[0];\n\n$\
name=\"query\";$seq=\"\";\n%p=blast_xml2profile($n\
ame,$seq,100, 0, 0, $blast);\n&output_profile (%p)\
;\n\n\nsub output_profile\n  {\n    my (%profile)=\
(@_);\n    my ($a);\n    for ($a=0; $a<$profile{n}\
; $a++)\n      {\n	\n	print \">$profile{$a}{name} \
$profile{$a}{comment}\\n$profile{$a}{seq}\\n\";\n \
     }\n    return;\n  }\nsub file_contains \n  {\\
n    my ($file, $tag, $max)=(@_);\n    my ($n);\n \
   $n=0;\n    \n    if ( !-e $file && ($file =~/$t\
ag/)) {return 1;}\n    elsif ( !-e $file){return 0\
;}\n    else \n      {\n	open (FC, \"$file\");\n	w\
hile ( <FC>)\n	  {\n	    if ( ($_=~/$tag/))\n	    \
  {\n		close (FC);\n		return 1;\n	      }\n	    el\
sif ($max && $n>$max)\n	      {\n		close (FC);\n		\
return 0;\n	      }\n	    $n++;\n	  }\n      }\n  \
  close (FC);\n    return 0;\n  }\n	    \n	  \nsub\
 file2string\n  {\n    my $f=@_[0];\n    my $strin\
g, $l;\n    open (F,\"$f\");\n    while (<F>)\n   \
   {\n\n	$l=$_;\n	#chomp ($l);\n	$string.=$l;\n   \
   }\n    close (F);\n    $string=~s/\\r\\n//g;\n \
   $string=~s/\\n//g;\n    return $string;\n  }\n\\
n\n\nsub tag2value \n  {\n    \n    my $tag=(@_[0]\
);\n    my $word=(@_[1]);\n    my $return;\n    \n\
    $tag=~/$word=\"([^\"]+)\"/;\n    $return=$1;\n\
    return $return;\n  }\n      \nsub hit_tag2pdbi\
d\n  {\n    my $tag=(@_[0]);\n    my $pdbid;\n    \
   \n    $tag=~/id=\"(\\S+)\"/;\n    $pdbid=$1;\n \
   $pdbid=~s/_//;\n    return $pdbid;\n  }\nsub id\
2pdbid \n  {\n    my $id=@_[0];\n  \n    if ($id =\
~/pdb/)\n      {\n	$id=~/pdb(.*)/;\n	$id=$1;\n    \
  }\n    $id=~s/[|�_]//g;\n    return $id;\n  }\ns\
ub set_blast_type \n  {\n    my $file =@_[0];\n   \
 if (&file_contains ($file,\"EBIApplicationResult\\
",100)){$BLAST_TYPE=\"EBI\";}\n    elsif (&file_co\
ntains ($file,\"NCBI_BlastOutput\",100)) {$BLAST_T\
YPE=\"NCBI\";}\n    else\n      {\n	$BLAST_TYPE=\"\
\";\n      }\n    return $BLAST_TYPE;\n  }\nsub bl\
ast_xml2profile \n  {\n    my ($name,$seq,$maxid, \
$minid, $mincov, $file)=(@_);\n    my (%p, $a, $st\
ring, $n);\n    \n\n\n    if ($BLAST_TYPE eq \"EBI\
\" || &file_contains ($file,\"EBIApplicationResult\
\",100)){%p=ebi_blast_xml2profile(@_);}\n    elsif\
 ($BLAST_TYPE eq \"NCBI\" || &file_contains ($file\
,\"NCBI_BlastOutput\",100)){%p=ncbi_blast_xml2prof\
ile(@_);}\n    else \n      {\n	print \"**********\
** ERROR: Blast Returned an unknown XML Format ***\
*******************\";\n	die;\n      }\n    for ($\
a=0; $a<$p{n}; $a++)\n      {\n	my $name=$p{$a}{na\
me};\n	$p{$name}{seq}=$p{$a}{seq};\n      }\n    r\
eturn %p;\n  }\nsub ncbi_blast_xml2profile \n  {\n\
    my ($name,$seq,$maxid, $minid, $mincov, $strin\
g)=(@_);\n    my ($L,$l, $a,$b,$c,$d,$nhits,@ident\
ifyerL);\n    \n    \n    $seq=~s/[^a-zA-Z]//g;\n \
   $L=length ($seq);\n    \n    %hit=&xml2tag_list\
 ($string, \"Hit\");\n    \n    \n    for ($nhits=\
0,$a=0; $a<$hit{n}; $a++)\n      {\n	my ($ldb,$id,\
 $identity, $expectation, $start, $end, $coverage,\
 $r);\n	my (%ID,%DE,%HSP);\n	\n	$ldb=\"\";\n\n	%ID\
=&xml2tag_list ($hit{$a}{body}, \"Hit_id\");\n	$id\
entifyer=$ID{0}{body};\n	\n	%DE=&xml2tag_list ($hi\
t{$a}{body}, \"Hit_def\");\n	$definition=$DE{0}{bo\
dy};\n	\n	%HSP=&xml2tag_list ($hit{$a}{body}, \"Hs\
p\");\n	for ($b=0; $b<$HSP{n}; $b++)\n	  {\n	    m\
y (%START,%END,%E,%I,%Q,%M);\n\n	 \n	    %START=&x\
ml2tag_list ($HSP{$b}{body}, \"Hsp_query-from\");\\
n	    %HSTART=&xml2tag_list ($HSP{$b}{body}, \"Hsp\
_hit-from\");\n	    \n	    %LEN=  &xml2tag_list ($\
HSP{$b}{body}, \"Hsp_align-len\");\n	    %END=  &x\
ml2tag_list ($HSP{$b}{body}, \"Hsp_query-to\");\n	\
    %HEND=  &xml2tag_list ($HSP{$b}{body}, \"Hsp_h\
it-to\");\n	    %E=&xml2tag_list     ($HSP{$b}{bod\
y}, \"Hsp_evalue\");\n	    %I=&xml2tag_list     ($\
HSP{$b}{body}, \"Hsp_identity\");\n	    %Q=&xml2ta\
g_list     ($HSP{$b}{body}, \"Hsp_qseq\");\n	    %\
M=&xml2tag_list     ($HSP{$b}{body}, \"Hsp_hseq\")\
;\n	    \n	    for ($e=0; $e<$Q{n}; $e++)\n\n	    \
  {\n		$qs=$Q{$e}{body};\n		$ms=$M{$e}{body};\n		i\
f ($seq eq\"\"){$seq=$qs;$L=length($seq);}\n		\n		\
$expectation=$E{$e}{body};\n		$identity=($LEN{$e}{\
body}==0)?0:$I{$e}{body}/$LEN{$e}{body}*100;\n		$s\
tart=$START{$e}{body};\n		$end=$END{$e}{body};\n		\
$Hstart=$HSTART{$e}{body};\n		$Hend=$HEND{$e}{body\
};\n	\n		$coverage=(($end-$start)*100)/$L;\n\n	\n	\
	if ($identity>$maxid || $identity<$minid || $cove\
rage<$mincov){next;}\n		@lr1=(split (//,$qs));\n		\
@lr2=(split (//,$ms));\n		$l=$#lr1+1;\n		for ($c=0\
;$c<$L;$c++){$p[$nhits][$c]=\"-\";}\n		for ($d=0,$\
c=0; $c<$l; $c++)\n		  {\n		    $r=$lr1[$c];\n		  \
  if ( $r=~/[A-Za-z]/)\n		      {\n			\n			$p[$nhi\
ts][$d + $start-1]=$lr2[$c];\n			$d++;\n		      }\\
n		  }\n		$Qseq[$nhits]=$qs;\n		$Hseq[$nhits]=$ms;\
\n		$QstartL[$nhits]=$start;\n		$HstartL[$nhits]=$\
Hstart;\n		$identityL[$nhits]=$identity;\n		$endL[\
$nhits]=$end;\n		$definitionL[$nhits]=$definition;\
\n		$identifyerL[$nhits]=$identifyer;\n		$comment[\
$nhits]=\"$ldb|$identifyer [Eval=$expectation][id=\
$identity%][start=$Hstart end=$Hend]\";\n		$nhits+\
+;\n	      }\n	  }\n      }\n    \n    $profile{n}\
=0;\n    $profile{$profile{n}}{name}=$name;\n    $\
profile{$profile{n}}{seq}=$seq;\n    $profile {n}+\
+;\n    \n    for ($a=0; $a<$nhits; $a++)\n      {\
\n	$n=$a+1;\n	\n	$profile{$n}{name}=\"$name\\_$a\"\
;\n	$profile{$n}{seq}=\"\";\n	$profile{$n}{Qseq}=$\
Qseq[$a];\n	$profile{$n}{Hseq}=$Hseq[$a];\n	$profi\
le{$n}{Qstart}=$QstartL[$a];\n	$profile{$n}{Hstart\
}=$HstartL[$a];\n	$profile{$n}{identity}=$identity\
L[$a];\n	$profile{$n}{definition}=$definitionL[$a]\
;\n	$profile{$n}{identifyer}=$identifyerL[$a];\n	$\
profile{$n}{comment}=$comment[$a];\n	for ($b=0; $b\
<$L; $b++)\n	  {\n	    if ($p[$a][$b])\n	      {\n\
		$profile{$n}{seq}.=$p[$a][$b];\n	      }\n	    e\
lse\n	      {\n		$profile{$n}{seq}.=\"-\";\n	     \
 }\n	  }\n      }\n    \n    $profile{n}=$nhits+1;\
\n    return %profile;\n  }\nsub ebi_blast_xml2pro\
file \n  {\n    my ($name,$seq,$maxid, $minid, $mi\
ncov, $string)=(@_);\n    my ($L,$l, $a,$b,$c,$d,$\
nhits,@identifyerL,$identifyer);\n    \n\n    \n  \
  $seq=~s/[^a-zA-Z]//g;\n    $L=length ($seq);\n  \
  %hit=&xml2tag_list ($string, \"hit\");\n    \n  \
  for ($nhits=0,$a=0; $a<$hit{n}; $a++)\n      {\n\
	my ($ldb,$id, $identity, $expectation, $start, $e\
nd, $coverage, $r);\n	my (%Q,%M,%E,%I);\n	\n	$ldb=\
&tag2value ($hit{$a}{open}, \"database\");\n	$iden\
tifyer=&tag2value ($hit{$a}{open}, \"id\");\n\n	$d\
escription=&tag2value ($hit{$a}{open}, \"descripti\
on\");\n	\n	%Q=&xml2tag_list ($hit{$a}{body}, \"qu\
erySeq\");\n	%M=&xml2tag_list ($hit{$a}{body}, \"m\
atchSeq\");\n	%E=&xml2tag_list ($hit{$a}{body}, \"\
expectation\");\n	%I=&xml2tag_list ($hit{$a}{body}\
, \"identity\");\n	\n\n	for ($b=0; $b<$Q{n}; $b++)\
\n	  {\n	    \n	    \n	    $qs=$Q{$b}{body};\n	   \
 $ms=$M{$b}{body};\n	    if ($seq eq\"\"){$seq=$qs\
;$L=length($seq);}\n\n	    $expectation=$E{$b}{bod\
y};\n	    $identity=$I{$b}{body};\n	    \n	    	  \
  \n	    $start=&tag2value ($Q{$b}{open}, \"start\\
");\n	    $end=&tag2value ($Q{$b}{open}, \"end\");\
\n	    $startM=&tag2value ($M{$b}{open}, \"start\"\
);\n	    $endM=&tag2value ($M{$b}{open}, \"end\");\
\n	    $coverage=(($end-$start)*100)/$L;\n	    \n	\
   # print \"$id: ID: $identity COV: $coverage [$s\
tart $end]\\n\";\n	    \n	    \n	    if ($identity\
>$maxid || $identity<$minid || $coverage<$mincov){\
next;}\n	    # print \"KEEP\\n\";\n\n	    \n	    @\
lr1=(split (//,$qs));\n	    @lr2=(split (//,$ms));\
\n	    $l=$#lr1+1;\n	    for ($c=0;$c<$L;$c++){$p[\
$nhits][$c]=\"-\";}\n	    for ($d=0,$c=0; $c<$l; $\
c++)\n	      {\n		$r=$lr1[$c];\n		if ( $r=~/[A-Za-\
z]/)\n		  {\n		    \n		    $p[$nhits][$d + $start-\
1]=$lr2[$c];\n		    $d++;\n		  }\n	      }\n	  \n	\
    \n	    $identifyerL[$nhits]=$identifyer;\n	   \
 $comment[$nhits]=\"$ldb|$identifyer [Eval=$expect\
ation][id=$identity%][start=$startM end=$endM]\";\\
n	    $nhits++;\n	  }\n      }\n    \n    $profile\
{n}=0;\n    $profile{$profile{n}}{name}=$name;\n  \
  $profile{$profile{n}}{seq}=$seq;\n    $profile {\
n}++;\n    \n    for ($a=0; $a<$nhits; $a++)\n    \
  {\n	$n=$a+1;\n	$profile{$n}{name}=\"$name\\_$a\"\
;\n	$profile{$n}{seq}=\"\";\n	$profile{$n}{identif\
yer}=$identifyerL[$a];\n	\n	$profile{$n}{comment}=\
$comment[$a];\n	for ($b=0; $b<$L; $b++)\n	  {\n	  \
  if ($p[$a][$b])\n	      {\n		$profile{$n}{seq}.=\
$p[$a][$b];\n	      }\n	    else\n	      {\n		$pro\
file{$n}{seq}.=\"-\";\n	      }\n	  }\n      }\n  \
  $profile{n}=$nhits+1;\n    \n    return %profile\
;\n  }\n\nsub blast_xml2hit_list\n  {\n    my $str\
ing=(@_[0]);\n    return &xml2tag_list ($string, \\
"hit\");\n  }\nsub xml2tag_list  \n  {\n    my ($s\
tring_in,$tag)=@_;\n    my $tag_in, $tag_out;\n   \
 my %tag;\n    \n    if (-e $string_in)\n      {\n\
	$string=&file2string ($string_in);\n      }\n    \
else\n      {\n	$string=$string_in;\n      }\n    \
$tag_in1=\"<$tag \";\n    $tag_in2=\"<$tag>\";\n  \
  $tag_out=\"/$tag>\";\n    $string=~s/>/>##1/g;\n\
    $string=~s/</##2</g;\n    $string=~s/##1/<#/g;\
\n    $string=~s/##2/#>/g;\n    @l=($string=~/(\\<\
[^>]+\\>)/g);\n    $tag{n}=0;\n    $in=0;$n=-1;\n \
 \n \n\n    foreach $t (@l)\n      {\n\n	$t=~s/<#/\
/;\n	$t=~s/#>//;\n	\n	if ( $t=~/$tag_in1/ || $t=~/\
$tag_in2/)\n	  {\n	 \n	    $in=1;\n	    $tag{$tag{\
n}}{open}=$t;\n	    $n++;\n	    \n	  }\n	elsif ($t\
=~/$tag_out/)\n	  {\n	    \n\n	    $tag{$tag{n}}{c\
lose}=$t;\n	    $tag{n}++;\n	    $in=0;\n	  }\n	el\
sif ($in)\n	  {\n	   \n	    $tag{$tag{n}}{body}.=$\
t;\n	  }\n      }\n  \n    return %tag;\n  }\n\n\n\
\n\n","use Env qw(HOST);\nuse Env qw(HOME);\nuse E\
nv qw(USER);\nwhile (<>)\n  {\n    if ( /^>(\\S+)/\
)\n      {\n	if ($list{$1})\n	  {\n	    print \">$\
1_$list{$1}\\n\";\n	    $list{$1}++;\n	  }\n	else\\
n	  {\n	    print $_;\n	    $list{$1}=1;\n	  }\n  \
    }\n    else\n      {\n	print $_;\n      }\n  }\
\n      \n","\n\n\nuse Env qw(HOST);\nuse Env qw(H\
OME);\nuse Env qw(USER);\n\n\nopen (F,$ARGV[0]);\n\
while ( <>)\n  {\n    @x=/([^:,;\\)\\(\\s]+):[^:,;\
\\)\\(]*/g;\n    @list=(@list,@x);\n  }\n$n=$#list\
+1;\nforeach $n(@list){print \">$n\\nsequence\\n\"\
;}\n\n\nclose (F);\n","\nopen (F, $ARGV[0]);\n\nwh\
ile ( <F>)\n  {\n    @l=($_=~/(\\S+)/g);\n    \n  \
  $name=shift @l;\n    \n    print STDOUT \"\\n>$n\
ame\\n\";\n    foreach $e (@l){$e=($e eq \"0\")?\"\
O\":\"I\";print \"$e\";}\n  }\nclose (F);\n\n		   \
    \n    \n","use Env qw(HOST);\nuse Env qw(HOME)\
;\nuse Env qw(USER);\n\n$tmp=\"$ARGV[0].$$\";\nope\
n (IN, $ARGV[0]);\nopen (OUT, \">$tmp\");\n\nwhile\
 ( <IN>)\n  {\n    $file=$_;\n    $file=~s/\\r\\n/\
\\n/g;\n    $file=~s/\\n\\r/\\n/g;\n    $file=~s/\\
\r\\r/\\n/g;\n    $file=~s/\\r/\\n/g;\n    print O\
UT \"$file\";\n  }\nclose (IN);\nclose (OUT);\n\no\
pen (OUT, \">$ARGV[0]\");\nopen (IN, \"$tmp\");\n\\
nwhile ( <IN>)\n{\n  print OUT \"$_\";\n}\nclose (\
IN);\nclose (OUT);\nunlink ($tmp);\n\n"};
