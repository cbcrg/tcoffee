char *PerlScriptName[]={"rec_sum.pl","count.pl","p\
rocess_list.pl","make_license.pl","CCsed.script","\
msa2bootstrap.pl","t_coffee_dpa","t_coffee_dpa2","\
tc_generic_method.pl","generic_method.tc_method","\
clustalw_method.tc_method","extract_from_pdb","ins\
tall.pl","clean_cache.pl","mocca","dalilite.pl","w\
ublast.pl","blastpgp.pl","ncbiblast_lwp.pl","wubla\
st_lwp.pl","RNAplfold2tclib.pl","fasta_seq2RNAplfo\
ld_templatefile.pl","fasta_seq2hmmtop_fasta.pl","f\
asta_seq2consan_aln.pl","clustalw_aln2fasta_aln.pl\
","msf_aln2fasta_aln.pl","blast_aln2fasta_aln.pl",\
"blast_xml2fasta_aln.pl","fasta_aln2fasta_aln_uniq\
ue_name.pl","newick2name_list.pl","excel2fasta.pl"\
,"any_file2unix_file.pl","EndList"};char *PerlScri\
ptFile[]={"use File::Copy;\nuse Env qw(HOST);\nuse\
 Env qw(HOME);\nuse Env qw(USER);\n$x_field=0;\n$y\
_field=1;\n$interval=0;\n$file=\"stdin\";\n$print_\
avg=1;\n$print_sd=0;\n$print_sum=0;\n$print_n=0;\n\
foreach $value ( @ARGV)\n    {\n	if ($value ne $AR\
GV[$np]) \n	    {\n	    ;\n	    }\n	elsif($value e\
q \"-print_all\")\n	    {\n	    $print_sd=$print_a\
vg=$print_n=$print_sum=1;\n	    $np++;\n	    }\n	e\
lsif($value eq \"-print_sum\")\n	    {\n	    $prin\
t_sum=1;\n	    $print_avg=0;\n	    $np++;\n	    }\\
n	elsif($value eq \"-print_n\")\n	    {\n	    $pri\
nt_n=1;\n	    $print_avg=0;\n	    $np++;\n	    }\n\
	elsif($value eq \"-print_avg\")\n	    {\n	    $pr\
int_avg=1;\n	    $print_avg=0;\n	    $np++;\n	    \
}\n	elsif($value eq \"-sd\")\n	    {\n	    $print_\
sd=1;\n	    $print_avg=0;\n	    $np++;\n	    }\n	e\
lsif($value eq \"-h\")\n	    {\n	    $header=1;\n	\
    $np++;\n	    }\n	elsif ($value eq \"-i\")\n	  \
  {\n	    $interval= $ARGV[++$np];\n	    $np++;\n \
   	    }\n	elsif ($value eq \"-r\")\n	    {\n	   \
 $min= $ARGV[++$np];\n	    $max= $ARGV[++$np];\n	 \
   $np++;\n    	    }\n	\n	elsif ($value eq \"-x\"\
)\n	    {\n	    $x_field= $ARGV[++$np]-1;\n	    $n\
p++;\n    	    }\n	elsif ($value eq \"-y\")\n	    \
{\n	      \n	    while ($ARGV[$np+1] && !($ARGV[$n\
p+1]=~/\\-/))\n	      {\n		$y_field[$nyf++]=$ARGV[\
++$np]-1;\n		$y_field_set=1;\n	      }\n\n	    $np\
++;\n    	    }\n	elsif ($value eq \"-file\")\n	  \
  {\n	    $file= $ARGV[++$np];\n	    $file_set=1;\\
n	    $np++;\n    	    }       \n	elsif ( $value e\
q \"h\" ||  $value eq \"-h\" || $value eq \"-H\" |\
| $value eq \"-help\" || $value eq \"help\")\n	  {\
\n	    print STDOUT \"data_analyse: Analyse and di\
scretization of data\\n\";\n	    print STDOUT \"  \
     -file:    <file containing the data to analyz\
e>,.<def=STDIN>\\n\";\n	    print STDOUT \"       \
-x: <field containing the X>,...............<Def=0\
>\\n\";\n	    print STDOUT \"       -y: <field con\
taining the Y>,...............<Def=1>\\n\";\n	    \
print STDOUT \"       -i:<Interval size on the X>,\
...............<Def=0>\\n\";\n	    print STDOUT \"\
       -i:<0:only one interval>\\n\";\n	    print \
STDOUT \"       -r:<Range of the X>\\n\";\n	    pr\
int STDOUT \"       -sd: print standard deviation \
on the Y\";\n	    print STDOUT \"       -h  : prin\
t column header \\n\";\n	    exit (0);\n	  }\n	els\
if ($value=~/-/)\n	  {\n	    print \"$value is not\
 a valid FLAG[FATAL]\\n\";\n	    exit (0);\n	   } \
\n	elsif ($list eq \"\") \n	    {\n	    $file=$ARG\
V[$np];\n	    $np++;\n	    }\n	\n	\n      }\n\n\n\\
n\n\nif ($file eq \"stdin\")\n	{\n	$remove_file=1;\
\n	$file=\"tmp$$\";\n	open (F, \">$file\");\n	whil\
e (<STDIN>)\n		{\n		print F $_;\n		}\n	close (F);\\
n	 \n	;}\n\n\nopen(F,$file);\n\nif ($interval)\n  \
{\n    $interval_size=($max-$min)/$interval;\n  }\\
nwhile (<F>)\n  {\n    $line=$_;\n    if (!/\\S/){\
next;}\n    @list=($line=~/(\\S+)/g);\n    \n    i\
f ($interval==0){$bin=0;}\n    else{$bin=int (($li\
st[$x_field]-$min)/($interval_size));}\n\n    \n  \
  if ($bin && $bin==$interval){$bin--;}\n    for (\
 $a=0; $a<$nyf; $a++)\n      {\n	$sum{$a}{$bin}+=$\
list[$y_field[$a]];\n	$sum2{$a}{$bin}+=$list[$y_fi\
eld[$a]]*$list[$y_field[$a]];\n	$n{$a}{$bin}++;\n \
     }\n  }\n\nif (!$interval){$interval=1;}\nfor \
( $a=0; $a<$interval; $a++)\n  {\n    printf ( \"%\
3d %3d \", $interval_size*$a, $interval_size*($a+1\
));\n    for ( $b=0; $b<$nyf; $b++)	\n      {\n	$i\
=$interval*$a;\n	if ( $n{$b}{$a}==0)\n	  {\n	    $\
avg=0;\n	    $sd=0;\n	  }\n	else\n	  {\n	    $avg=\
$sum{$b}{$a}/$n{$b}{$a};\n	    $sd=sqrt($sum2{$b}{\
$a}*$n{$b}{$a}-$sum{$b}{$a}*$sum{$b}{$a})/($n{$b}{\
$a}*$n{$b}{$a});\n	  }\n	if ($print_n) {printf \"%\
10.4f \", $n{$b}{$a};}\n	if ($print_sum){printf \"\
%10.4f \", $sum{$b}{$a};}\n	if ($print_avg){printf\
 \"%10.4f \", $avg}\n	if ($print_sd) {printf \"%10\
.4f \", $sd;}\n      }\n    printf (\"\\n\");\n  }\
\n\n\nif ( $remove_file){unlink $file;}\n","use Fi\
le::Copy;\nuse Env qw(HOST);\nuse Env qw(HOME);\nu\
se Env qw(USER);\n\nforeach $v (@ARGV){$cl.=$v;}\n\
\n\nif ( $cl=~/-k(\\d+)/){$k=$1;}\nelse {$k=1;}\ni\
f ( $cl=~/-w(\\d+)/){$w=$1;}\nelse {$w=-1;}\nif ( \
$cl=~/-p(\\d+)/){$p=$1;}\nelse {$p=-1;}\n\nwhile (\
<STDIN>)\n  {\n    @l=($_=~/(\\S+)/g);\n    $v=$l[\
$k-1];\n    if ( !$h{$v}){@ll=($v, @ll);}\n    \n \
   if ( $w==-1)\n      {$h{$v}++;}\n    else\n    \
  {$h{$v}+=$l[$w-1];}\n\n    if ($p!=-1){$print{$v\
}=$l[$p-1];}\n\n  }\nforeach $v (@ll)\n  {\n    pr\
int \"$v $print{$v} $h{$v}\\n\";\n  }\n","\nuse En\
v qw(HOST);\nuse Env qw(HOME);\nuse Env qw(USER);\\
n$random_tag=int (rand 10000)+1;\n$unique_prefix=\\
"$$.$HOST.$random_tag\";\n$queue=\"distillery.and.\
mid\";\n$monitor=0;\n$stderr_file=\"/dev/null\";\n\
$stdio_file=\"/dev/null\";\n$log_file=\"/dev/null\\
";\n$pause_time=0;\n$max_sub_jobs=60;\n$min_sub_jo\
bs=30;\n$output_all=0;\n$var='\\$';\n\nforeach $va\
lue ( @ARGV)\n    {\n	if ($value ne $ARGV[$np]) \n\
	    {\n	    ;\n	    }\n	elsif ($value eq \"-max_s\
ub_jobs\")\n	    {\n	    $max_sub_jobs= $ARGV[++$n\
p];\n	    $np++;\n    	    }	\n	elsif ($value eq \\
"-min_sub_jobs\" )\n	    {\n	    $min_sub_jobs= $A\
RGV[++$np];\n	    $np++;\n    	    }\n	elsif ($val\
ue eq \"-para\")\n	    {\n	    $para=1;\n	    $mon\
itor=1;\n	    $np++;\n    	    }\n	elsif ($value e\
q \"-monitor\") \n	    {\n	    $monitor=1;\n	    $\
np++;\n	    }\n	elsif ($value eq \"-no_monitor\") \
\n	    {\n	    $monitor=0;\n	    $np++;\n	    }\n	\
elsif ($value eq \"-queue\")\n	    {\n	    $queue=\
$ARGV[++$np];\n	    $np++;\n	    }	\n	elsif ($valu\
e eq \"-stderr_file\")\n	    {\n	    $stderr_file=\
$ARGV[++$np];\n	    $np++;\n	    }\n	elsif ($value\
 eq \"-stdio_file\")\n	    {\n	    $stdio_file=$AR\
GV[++$np];\n	    $np++;\n	    }\n	elsif ($value eq\
 \"-output_all\")\n	    {\n	    $output_all=1;\n	 \
   $np++;\n	    }\n	elsif ($value eq \"-pause\") \\
n	    {\n	    $pause_time=$ARGV[++$np];\n	    $np+\
+;\n	    }\n	elsif ($value eq \"-log\")\n	      {\\
n	       $log=1;\n	       \n	       if ($ARGV[$np+\
1]=~/\\-\\S+/) \n	          {\n		  $log_file=\"std\
err\";\n	          }\n	       else \n	          {\\
n		  $log_file=$ARGV[++$np]; \n		  ++$np;\n		 \n	 \
         }\n	      }\n	elsif ( $value eq \"-com\")\
\n	    {\n		\n		if (!$ARGV[$np+1]=~/^\\'/) { $com=\
$ARGV[++$np];}\n		else {$com=$ARGV[++$np];}\n\n	  \
   $np++;\n	    }\n	elsif ( $value eq \"-check\")\\
n	  {\n	    \n	    if (!$ARGV[$np+1]=~/^\\'/) { $c\
heck=$ARGV[++$np];}\n	    else {$check=$ARGV[++$np\
];}\n	    $np++;\n	  }\n	elsif ($com eq \"\") \n	 \
   {\n	    $com_set=1;\n	    $com=$ARGV[$np];\n	  \
  \n	    $np++;\n	    }\n	elsif ($list eq \"\") \n\
	    {\n	    $list_set=1;\n	    $list=$ARGV[$np];\\
n	    $np++;\n	    }\n	elsif ( $var_set eq \"\")\n\
	    {\n	    $var_set=1;\n	    $var=$ARGV[$np];\n	\
    $np++;\n	    }\n	}\n\n\n\n\nif ( $com eq \"\")\
{print \"You Need to Provide a Command [FATAL]\\n\\
";\n	      die;\n	     }\n\n\n\nif ($list_set==0) \
\n    {\n    $x= int (rand 100000)+1;\n    $tmp_fi\
le_name=\"tmp_file_$x\";\n    open ( TMP, \">$tmp_\
file_name\");\n    while (<STDIN>)\n      {\n	prin\
t TMP $_;\n      }\n    close (TMP);\n    open (F,\
 $tmp_file_name);\n    }\nelse \n    {\n    open (\
F, $list);\n    }\n\nif ($para==0) \n    {\n\n    \
 @tc_list= <F>;\n     close (F); \n     \n     for\
each $val(@tc_list) \n	    {\n	      \n	      \n	 \
     \n	      $loc_com=$com;\n	      if ($check){$\
loc_check=$check;}\n	      \n	      @i_val=($val=~\
/([^\\s]+)/g);\n	      \n	      if ( $#i_val==0)\n\
		{\n		  if ($check){$loc_check=~s/$var/$i_val[0]/\
g;}\n		  $loc_com=~s/$var/$i_val[0]/g;\n		}\n	    \
  else\n		{\n		  for ($n=1; $n<=$#i_val+1;$n++ )\n\
		    {\n		      \n		      $sub=\"$var$n\";\n		   \
   \n		      $loc_com=~s/$sub/$i_val[$n-1]/g;\n		 \
     if ($check){$loc_check=~s/$var/$i_val[0]/g;}\\
n		    }\n		}\n	      if ( $check && -e $loc_check\
)\n		{\n		  print STDERR \"skipping $loc_com...\\n\
\";\n		  }\n	      else\n		{\n		  system \"$loc_co\
m\";\n		}\n	    }\n    exit;\n    }\n\nelsif ($par\
a==1) \n    {\n    print STDERR \"do parallel exec\
ution of: \\\"$com $list\\\"\\n\";\n    \n    if (\
$log==1) \n	{\n	if ($log_file eq \"stdout\" || $lo\
g_file eq \"stderr\" ) \n		{\n		$log_file=\"\";\n	\
        }\n\n        else \n		{\n		system \"echo L\
OG FILE> $log_file\";\n		\n	        }\n	}\n    els\
e	\n	{\n	open ( OUT, \">/dev/null\");\n	}\n	\n    \
\n    $id=0;\n    $n_sub=0;\n    while ($val=<F>) \
\n	    {	    	    \n	    $job_log[$id]=\"$HOME/tmp\
/$unique_prefix.$id.log_file\";\n	    \n	    $job=\
$unique_prefix.\"_$id\";\n	    open (JOB, \">$job\\
");\n	    \n	    $loc_com=$com;\n	    chop $val;\n\
\n	    $loc_com=~s/\\$/$val/g;\n	 \n	    print JOB\
 \"#!/bin/csh\\n\";\n	    print JOB \"#\\$ -cwd\\n\
\";\n	    print JOB \"#\\$ -N $unique_prefix\\n\";\
\n	    if ($queue && !($queue eq \" \")) {print JO\
B \"#\\$ -l $queue\\n\";}\n	    print JOB \"#\\n\"\
;	    \n            print JOB \"$loc_com\\n\";\n	 \
   print JOB \"echo FINISHED  >> $job_log[$id]\\n\\
";\n	    print JOB \"pwd\\n\";\n	    \n	    close \
(JOB);\n	    if ( $output_all==1)\n		{\n		system \\
"qsub $job >  $unique_prefix\";		\n	        }\n	  \
  else\n		{system \"qsub $job -e $stderr_file -o $\
stdio_file >$unique_prefix\";	        \n	        }\
 \n\n\n\n	    print STDERR \"$id: $output_all\\n\"\
;\n	    $n_sub++;\n	    if ( $max_sub_jobs && $n_s\
ub==$max_sub_jobs) \n		{\n		$n_sub=monitor_process\
($min_sub_jobs,@job_log); 		 \n		\n	        }	\n	 \
  \n            unlink $unique_prefix;\n	    sleep\
 $pause_time;\n	    $id++;\n	    }\n\n    close (O\
UT);\n    close (F);\n\n    print STDERR \"Your $i\
d Jobs Have Been Submited (NAME=$unique_prefix)\\n\
\";\n    monitor_process (0, @job_log);\n    forea\
ch $file(@job_log) {if (-e $file) {unlink($file);}\
}\n    \n    }\n\nsub monitor_process ( @job_list)\
\n    {\n    my (@job_list)=@_;\n    my $min_sub_j\
obs=shift (@job_list);\n    my $n_sub_jobs;\n    m\
y $finished;\n    my $n=0;\n\n    $n_sub_jobs=-1;\\
n    $finished=0;\n    print STDERR \"\\nMonitor B\
atch: [$min_sub_jobs]\";\n       \n    while (!$fi\
nished && (($n_sub_jobs>$min_sub_jobs)|| $n_sub_jo\
bs==-1) ) \n	{\n	$finished=1;\n	$n_sub_jobs=0;\n	$\
n=0;\n	foreach $file (@job_list)\n	        {\n	\n	\
	if (-e $file){;}\n		else \n		    {\n		    $finish\
ed=0; $n_sub_jobs++;\n	            }\n	        }\n\
	system \"sleep 1\";\n        }\n    \n    return \
$n_sub_jobs;\n    }\n    \n    \nif ($tmp_file_nam\
e){unlink($tmp_file_name);}\n","\n\nforeach ($np=0\
; $np<=$#ARGV; $np++)\n    {\n    $value=$ARGV[$np\
];\n\n    if ($value eq \"-file\")\n      {\n     \
 $file= $ARGV[++$np];\n      }\n    elsif ($value \
eq \"-type\")\n      {\n        $type= $ARGV[++$np\
];\n      }\n    elsif ($value eq \"-institute\")\\
n      {\n        $institute= $ARGV[++$np];\n     \
 }\n    elsif ($value eq \"-author\")\n      {\n  \
      $author= $ARGV[++$np];\n      }\n    elsif (\
$value eq \"-date\")\n      {\n        $date= $ARG\
V[++$np];\n      }\n     elsif ($value eq \"-progr\
am\")\n      {\n        $program= $ARGV[++$np];\n \
     }\n    elsif ($value eq \"-email\")\n      {\\
n        $email= $ARGV[++$np];\n      }\n    else\\
n      {\n	print \"$value is an unkown argument[FA\
TAL]\\n\";\n	exit (1);\n      }\n  }\n\n\n\nopen F\
, $file || die;\nprint $INSTITUTE;\nif ( $type eq \
\"c\"){print \"/*********************************C\
OPYRIGHT NOTICE**********************************/\
\\n\";}\nif ( $type eq \"perl\"){print \"#########\
########################COPYRIGHT NOTICE##########\
#######################/\\n\";}\nif ( $type eq \"t\
xt\"){print \"----------------------------------CO\
PYRIGHT NOTICE---------------------------------/\\\
n\";}\n\n\nwhile (<F>)\n  {\n  s/\\$INSTITUTE/$ins\
titute/g;\n  s/\\$AUTHOR/$author/g;\n  s/\\$DATE/$\
date/g;\n  s/\\$PROGRAM/$program/g;  \n  s/\\$EMAI\
L/$email/g;  \n  if ( $type eq \"txt\"){print $_;}\
\n  elsif ($type eq \"c\"){chop $_; print \"\\/*$_\
*\\/\\n\";}\n  elsif ($type eq \"perl\"){print \"\\
\#$_\";}\n}\nclose (F);\nif ( $type eq \"c\"){prin\
t \"/*********************************COPYRIGHT NO\
TICE**********************************/\\n\";}\nif\
 ( $type eq \"perl\"){print \"####################\
#############COPYRIGHT NOTICE#####################\
############/\\n\";}\nif ( $type eq \"txt\"){print\
 \"----------------------------------COPYRIGHT NOT\
ICE---------------------------------/\\n\";}\n\n",\
"\nwhile (<>)	\n	{\n	s/\\=cc/123456789/g;\n	s/\\bc\
c/\\$\\(CC\\)/g;\n	s/123456789/\\=cc/g;\n	print $_\
;\n	}\n\n","$version=\"1.00\";\n$rseed= int(rand(1\
00000))+1;\n\n\nif ( $#ARGV==-1)\n  {\n    print \\
"msa2bootstrap -i <input_file> -input <seq|msa|mat\
rix|tree> -n <N-Boostrap> -o <outtree> -tmode <nj|\
upgma|parsimony|ml> -dmode <kimura> -alignpg <t_co\
ffee | muscle | clustalw> -rtree <file> -stype <pr\
ot|cdna|dna> -recompute -system <cygwin|unix>\";\n\
    print \"\\n\\t-i: input file, can be sequneces\
, msa, matrix, trees, type is specified via -input\
\";\n    print \"\\n\\t-input: Type of input data\\
";\n    print \"\\n\\t\\tmsa: msa in fasta format\\
";\n    print \"\\n\\t\\tseq: compute an msa with \
-alignpg\";\n    print \"\\n\\t\\tmatrix: phylipp \
distance matrix fed directly to method -tmode [cav\
eat: tmode=nj or upgma]\";\n    print \"\\n\\t\\tt\
ree: list of newick trees directly fed to consence\
 in order to generate a bootstraped tree\";\n    \\
n    print \"\\n\\t-n: number of bootstrap replica\
tes\";\n    print \"\\n\\t-o: name of the output t\
ree. Files are not overwritten. Use -recompute to \
overwrite existing file\";\n    print \"\\n\\t-tmo\
de: tree mode: nj|upgma|parsimony|ml\";\n    print\
 \"\\n\\t-dmode: distance mode\";\n    print \"\\n\
\\t-alignpg: program for aligning sequences (t_cof\
fee=default)\";\n    print \"\\n\\t-rtree: replica\
te tree file (default: no file)\";\n    print \"\\\
n\\t-rmsa: replicate msa file (default: no file)\"\
;\n    print \"\\n\\t-rmat: replicate matrix file \
(default: no file)\";\n    print \"\\n\\t-stype: s\
equence type: protein, dna or cdna\";\n    print \\
"\\n\\t-recompute: force files to be overwritten\"\
;\n    print \"\\n\\t-system: cygwin|unix\";\n    \
  \n\n    \n    &my_exit (EXIT_FAILURE);\n  }\nfor\
each $arg (@ARGV){$command.=\"$arg \";}\n\nprint \\
"CLINE: $command\\n\";\n$threshold=100;\n$trim_msa\
=0;\n$stype=\"prot\";\nprint \"msa2bootstrap \";\n\
\n$system=\"cygwin\";\nif(($command=~/\\-system (\\
\S+)/))\n  {\n    $system=$1;\n    if ( $system eq\
 \"cygwin\")\n      {\n	$exec_extension=\".exe\";\\
n      }\n    elsif ( $system eq \"unix\")\n      \
{\n	$exec_extension=\"\";\n	print \"system=Unix\";\
die;\n      }\n    else\n      {\n	print \"msa2boo\
strap: -system=$system is an unknown mode [FATAL]\\
\n\"; die;\n      }\n    \n    print \"-system $sy\
stem \";\n  }\nif(($command=~/\\-stype (\\S+)/))\n\
  {\n    $stype=$1;\n  }\nprint \"-stype=$stype \"\
;\n\n\n\nif(($command=~/\\-i (\\S+)/))\n  {\n    $\
msa=$1;\n    print \"-i $msa \";\n  }\n\nif(($comm\
and=~/\\-rtree (\\S+)/))\n  {\n    $rtree=$1;\n   \
 print \"-rtree=$rtree \";\n  }\n\nif(($command=~/\
\\-rmsa (\\S+)/))\n  {\n    $rmsa=$1;\n  }\nif(($c\
ommand=~/\\-rmat (\\S+)/))\n  {\n    $rmat=$1;\n  \
}\n$input=\"seq\";\nif(($command=~/\\-input (\\S+)\
/))\n  {\n    $input=$1;\n  }\nprint \"-input=$inp\
ut \";\n\n$dmode=\"kimura\";\nif(($command=~/\\-dm\
ode (\\S+)/))\n  {\n    $dmode=$1;\n  }\nprint \"-\
dmode=$dmode \";\n$alignpg=\"muscle\";\nif(($comma\
nd=~/\\-alignpg (\\S+)/))\n  {\n    $alignpg=$1;\n\
  }\nprint \"-alignpg=$dmode \";\n\n$tmode=\"nj\";\
\nif(($command=~/\\-tmode (\\S+)/))\n  {\n    $tmo\
de=$1;\n  }\nprint \"-tmode=$tmode \";\n$recompute\
=0;\nif(($command=~/\\-recompute/))\n  {\n    $rec\
ompute=1;\n    print \"-recompute \";\n  }\n\n$out\
=$msa;\n$out=~s/\\..*//;\n$out.=\".bph\";\nif(($co\
mmand=~/\\-o (\\S+)/))\n  {\n    $out=$1;\n    \n \
 }\nprint \"-out=$out \";\nif (-e $out && !$recomp\
ute)\n  {\n    print \"\\nNo Computation Required \
$out already exists\\n\";\n    &my_exit (EXIT_SUCC\
ESS);\n    \n  }\n\n$n=100;\nif(($command=~/\\-n (\
\\d+)/))\n  {\n    $n=$1;\n  }\nprint \"-n=$n \";\\
n$seed=3;\nif(($command=~/\\-s (\\d+)/))\n  {\n   \
 $seed=$1;\n  }\nprint \"-s=$seed \";\n\nif(($comm\
and=~/\\-run_name (\\d+)/))\n  {\n    $suffix=$1;\\
n  }\nelse\n  {\n    $msa=~/([^.]+)/;\n    $suffix\
=$1;\n  }\nprint \"-run_name=$suffix\\n\";\n\n\nif\
 ( $input eq \"seq\")\n  {\n    $seq=$msa;\n    $m\
sa=\"$suffix.prot_msa\";\n    \n    if ($stype eq \
\"cdna\")\n      {\n	$cdna_seq=$seq;\n	$clean_cdna\
_seq=&vtmpnam();\n	$seq=&vtmpnam();\n	`t_coffee -o\
ther_pg seq_reformat -in $cdna_seq -action +clean_\
cdna >$clean_cdna_seq`;\n	`t_coffee -other_pg seq_\
reformat -in $clean_cdna_seq -action +translate >$\
seq`;\n	\n      }\n\n    if (!-e $msa || $recomput\
e)\n      {\n	print \"\\n#####   Compute an MSA Wi\
th $alignpg\\n\";\n	\n	if ( $alignpg eq \"t_coffee\
\")\n	  {`$alignpg $seq -outfile=$msa >/dev/null 2\
>/dev/null`;}\n	elsif ( $alignpg eq \"muscle\")\n	\
  {\n	    `$alignpg -in $seq > $msa 2>/dev/null`;\\
n	  }\n	elsif ( $alignpg eq \"clustalw\")\n	  {\n	\
    `$alignpg -infile=$seq -outfile=$msa -quicktre\
e >/dev/null 2>/dev/null`;\n	  }\n	elsif ( $align \
eq \"mafft\")\n	  {\n	    `$alignpg $seq > $msa >/\
dev/null 2>/dev/null`;\n	  }\n	else\n	  {\n	    `$\
alignpg -in=$seq -outfile=$msa`;\n	  }\n      }\n \
   if (!-e $msa)\n      {\n	print \"\\nError: $ali\
gnpg Could Not produce the MSA $msa [FATAL]\\n\";\\
n      }\n\n    if ($stype eq \"cdna\")\n      {\n\
	$msa2=\"$suffix.cdna_msa\";\n	`t_coffee -other_pg\
 seq_reformat -in $clean_cdna_seq -in2 $msa -actio\
n +thread_dna_on_prot_aln -output fasta_aln  >$msa\
2`;\n	$msa=$msa2;\n      }\n    \n    $input=\"msa\
\";\n  }\n\n\n\n$seqboot_o=&vtmpnam();\n$seqboot_c\
=&vtmpnam();\n\n$protdist_o=&vtmpnam();\n$protdist\
_c=&vtmpnam();\nif ( $input eq \"msa\")\n  {\n    \
if ($tmode eq \"nj\" || $tmode eq \"upgma\"){$inpu\
t=\"matrix\";}\n    \n    $lmsa= &vtmpnam ();\n   \
 `t_coffee -other_pg seq_reformat -in $msa -output\
 phylip_aln > $lmsa`;\n    \n    if ( -e \"outfile\
\"){unlink (\"outfile\");}\n    # run seqboot\n  \\
n    if ( $n>1)\n      {\n	print \"Run SeqBoot ...\
..\";\n	open (F, \">$seqboot_c\");\n	print F \"$lm\
sa\\nR\\n$n\\nY\\n$seed\\n\";\n	close (F);\n	`seqb\
oot$exec_extension  < $seqboot_c`;\n	if ( -e \"out\
file\"){ print \"[OK]\\n\";}\n	else { print \"[FAI\
LED]\\n\";&my_exit (EXIT_FAILURE);}\n	`mv outfile \
$seqboot_o`;\n      }\n    else\n      {\n	`cp $lm\
sa $seqboot_o`;\n      }\n\n    if ($rmsa){`cp $se\
qboot_o $rmsa`;}\n    \n    if ($tmode eq \"nj\" |\
| $tmode eq \"upgma\")\n      {\n	if ( $stype eq \\
"prot\")\n	  {\n	    # run protdist\n	    print \"\
Run Protdist [dmode=$dmode]\";\n	    if ($dmode eq\
 \"kimura\")\n	      {\n		$dmode=\"P\\nP\\nP\";\n	\
      }\n	    else\n	      {\n		print \"\\n$dmode \
is an unknown mode for Protdist [FATAL:msa2bootstr\
ap.pl]\\n\";\n		&my_exit (EXIT_FAILURE);\n	      }\
\n	    open (F, \">$protdist_c\");\n	    if ($n>1)\
{print F \"$seqboot_o\\n$dmode\\nM\\nD\\n$n\\nY\\n\
\";}\n	    else {printf F \"$seqboot_o\\n$dmode\\n\
Y\\n\";}\n	    close (F);\n	    `protdist$exec_ext\
ension  < $protdist_c`;\n	    if ( -e \"outfile\")\
{ print \"[OK]\\n\";}\n	    else { print \"[FAILED\
]\\n\";&my_exit (EXIT_FAILURE);}\n	    `mv outfile\
 $protdist_o`;\n	 \n	  }\n	elsif ( $stype eq \"cdn\
a\" || $stype eq \"dna\")\n	  {\n	    print \"Run \
dnadist [dmode=default\";\n	    open (F, \">$protd\
ist_c\");\n	    if ($n>1){print F \"$seqboot_o\\nM\
\\nD\\n$n\\nY\\n\";}\n	    else {printf F \"$seqbo\
ot_o\\nY\\n\";}\n	    close (F);\n	    `protdist$e\
xec_extension  < $protdist_c`;\n	    if ( -e \"out\
file\"){ print \"[OK]\\n\";}\n	    else { print \"\
[FAILED]\\n\";&my_exit (EXIT_FAILURE);}\n	    `mv \
outfile $protdist_o`;\n	  }\n      }\n  }\nelsif (\
 $input eq \"matrix\")\n  {\n    $protdist_o=&vtmp\
nam();\n    print \"MSA: $msa\\n\";\n    `cp $msa \
$protdist_o`;\n    $n=1;\n  }\n\n\n\n\n\n$nb_o=&vt\
mpnam();\n$nb_c=&vtmpnam();\nif ($input eq \"matri\
x\" && $tmode ne \"parsimony\" && $tmode ne \"ml\"\
)\n  {\n    print \"Run neighbor [tmode=$tmode]\";\
\n\n    if ($tmode eq \"nj\")\n      {\n	$tmode=\"\
\\nN\\nN\";\n      }\n    elsif ( $tmode eq \"upgm\
a\")\n      {\n	$tmode = \"\\nN\";\n      }\n    e\
lse\n      {\n	print \"\\n ERROR: $tmode is an unk\
nown tree computation mode\\n\";\n	&my_exit (EXIT_\
FAILURE);\n      }\n\n    open (F, \">$nb_c\");\n \
   if ($n>1){print F \"$protdist_o$tmode\\nM\\n$n\\
\n$seed\\nY\\n\";}\n    else {print F \"$protdist_\
o$tmode\\nY\\n\";}\n    close (F);\n\n    `neighbo\
r$exec_extension  < $nb_c`;\n    if ( -e \"outtree\
\"){ print \"[Neighbor OK]\\n\";}\n    else { prin\
t \"[FAILED]\\n\";&my_exit (EXIT_FAILURE);}\n    `\
mv outtree $nb_o`;\n    unlink (\"outfile\");\n  }\
\nelsif ($input eq \"msa\" && $tmode eq \"parsimon\
y\")\n  {\n    if ( -e \"outfile\"){unlink (\"outf\
ile\");}\n    if ( -e \"outtree\"){unlink (\"outtr\
ee\");}\n    \n    if ($stype eq \"prot\")\n      \
{\n	print \"Run protpars [tmode=$tmode]\";\n	open \
(F, \">$nb_c\");\n	if ($n>1){print F \"$seqboot_o\\
\nM\\nD\\n$n\\n$seed\\n10\\nY\\n\";}\n	else {print\
 F \"$seqboot_o\\nY\\n\";}\n	close (F);\n	`protpar\
s$exec_extension  < $nb_c`;\n      }\n    elsif ( \
$stype eq \"dna\" || $stype eq \"cdna\")\n      {\\
n	print \"Run dnapars [tmode=$tmode]\";\n	open (F,\
 \">$nb_c\");\n	if ($n>1){print F \"$seqboot_o\\nM\
\\nD\\n$n\\n$seed\\n10\\nY\\n\";}\n	else {print F \
\"$seqboot_o\\nY\\n\";}\n	close (F);\n	`dnapars$ex\
ec_extension  < $nb_c`;\n      }\n    if ( -e \"ou\
ttree\"){ print \"[OK]\\n\";}\n    else { print \"\
[FAILED]\\n\";&my_exit (EXIT_FAILURE);}\n    `mv o\
uttree $nb_o`;\n   unlink (\"outfile\");\n  }\nels\
if ($input eq \"msa\" && $tmode eq \"ml\")\n  {\n \
   if ( -e \"outfile\"){unlink (\"outfile\");}\n  \
  if ( -e \"outtree\"){unlink (\"outtree\");}\n   \
 \n    if ($stype eq \"prot\")\n      {\n	print \"\
Error: ML impossible with Protein Sequences [ERROR\
]\";\n	&my_exit (EXIT_FAILURE);\n      }\n    elsi\
f ( $stype eq \"dna\" || $stype eq \"cdna\")\n    \
  {\n	print \"Run dnaml [tmode=$tmode]\";\n	open (\
F, \">$nb_c\");\n	if ($n>1){print F \"$seqboot_o\\\
nM\\nD\\n$n\\n$seed\\n10\\nY\\n\";}\n	else {print \
F \"$seqboot_o\\nY\\n\";}\n	close (F);\n	`dnaml$ex\
ec_extension  < $nb_c`;\n      }\n    if ( -e \"ou\
ttree\"){ print \"[OK]\\n\";}\n    else { print \"\
[FAILED]\\n\";&my_exit (EXIT_FAILURE);}\n    `mv o\
uttree $nb_o`;\n   unlink (\"outfile\");\n  }\n\n\\
nelse\n  {\n    `cp $msa $nb_o`;\n    $n=2;\n  }\n\
\nif ($rmsa && -e $seqboot_o){print \"\\nOutput Li\
st of $n Replicate MSA: $rmsa\\n\";`cp $seqboot_o \
$rmsa`;}\nif ($rmat && -e $protdist_o){print \"\\n\
Output List of $n Replicate MATRICES: $rmat\\n\";`\
cp $protdist_o $rmat`;}\nif ($rtree && -e $nb_o){p\
rint \"\\nOutput List of $n Replicate TREES: $rtre\
e\\n\";`cp $nb_o $rtree`;}\n\n\n\n$con_o=&vtmpnam(\
);\n$con_c=&vtmpnam();\nif ($n >1)\n  {\n    print\
 \"Run Consense.....\";\n    open (F, \">$con_c\")\
;\n    print F \"$nb_o\\nY\\n\";\n    close (F);\n\
    `consense$exec_extension  < $con_c`;\n    if (\
 -s \"outtree\"  > 0) { print \"[OK]\\n\";}\n    e\
lse { print \"[FAILED]\\n\";&my_exit (EXIT_FAILURE\
);}\n    `mv outtree $con_o`;\n    unlink (\"outfi\
le\");\n  }\nelse\n  {\n    `cp $nb_o $con_o`;\n  \
}\n\n\n`cp $con_o $out`;\nif ( !-e $out)\n  {\n   \
 print \"Tree Computation failed [FAILED]\\n\";\n \
   &my_exit (EXIT_FAILURE);\n  }\nelsif ($n>1)\n  \
{\n    print \"\\nOutput Bootstrapped Tree: $out\\\
n\";\n    $avg=`t_coffee -other_pg seq_reformat -i\
n $out -action +avg_bootstrap`;\n    $avg=~s/\\n//\
g;\n    print \"$avg\\n\";\n  }\nelse\n  {\n    pr\
int \"\\nOutput Tree: $out\\n\";\n  }\n\nopen (F, \
\"$out\");\nwhile (<F>)\n  {\n    \n    $tree.=$_;\
\n  }\nclose (F);\n$tree=~s/\\n//g;\nprint \"BPH: \
$tree\\n\";\n\n\n&my_exit (EXIT_SUCCESS);\n\nsub m\
y_exit \n  {\n    my $m=@_[0];\n    &clean_vtmpnam\
();\n    exit ($m);\n  }\nsub vtmpnam \n  {\n    m\
y $file;\n\n\n    $ntmp++;\n    $file=\"tmp4msa2bo\
otstrap.$rseed.$$.$ntmp\";\n    \n    push (@tmpfi\
le, $file);\n    return $file;\n  }\nsub clean_vtm\
pnam \n  {\n    my $t;\n    foreach $t (@tmpfile)\\
n      {\n	if ( -e $t){unlink ($t)};\n      }\n  }\
\n","use Env;\n$seq_reformat=\"t_coffee -other_pg \
seq_reformat \";\n$VersionTag=\"1.00\";\n$step=1;\\
n$unset=\"\";\n$scoreT1=$scoreT2=$nseqT=$dp_limit=\
$unset;\n@tl=();\nchomp($tc_version=`t_coffee -ver\
sion`);$tc_version=~s/PROGRAM: //;\n\n\nprint STDE\
RR \"\\n******************************************\
***********************\";\nprint STDERR \"\\n*   \
        HIGH LEVEL PROGRAM: T-COFFEE_DPA Version $\
VersionTag\";\nprint STDERR \"\\n*           LOW  \
LEVEL PROGRAM: $tc_version \";\nprint STDERR \"\\n\
**************************************************\
***************\";\n\nif (!@ARGV)\n  {\n    print \
\"t_coffee_dpa accepts every t_coffee_flag.\\nType\
 t_coffee to obtain a list\\n\";\n    print \"Requ\
ires $TC_VERSION\\n\";\n    print \"Requires \";\n\
    print \"t_coffee_dpa specific flags:\\n\";\n  \
  print \"\\t-dpa_master_aln....................Ma\
ster alignment: provided OR computed\\n\";\n    pr\
int \"\\t-dpa_master_aln....................By def\
ault, Computed with t_coffee -very_fast\\n\";\n   \
 print \"\\t-dpa_master_aln=<file>.............Use\
 file, (must be an aln in Fasta or ClustalW\\n\";\\
n    print \"\\t-dpa_master_aln=<program>.........\
.Compute aln with pg -in seq -out aln`\\n\";\n    \
print \"\\t-dpa_maxnseq.......................Maxi\
mum number of sequences in subgroups\\n\";\n    pr\
int \"\\t-dpa_min_score1....................Minimu\
m Id for two sequences to be grouped in ref_aln\\n\
\";\n    print \"\\t-dpa_min_score2...............\
.....Minimum Id within a subgroup\\n\";\n    print\
 \"\\t-dpa_debug.........................Keep Tmp \
File (for debug purpose)\\n\\n\";\n    \n    exit \
(0);\n  }\nforeach $arg (@ARGV)\n  {\n    $arg_lis\
t.=\" $arg\";\n  }\n$arg_list=~s/[=,;]/ /g;\n\n\n(\
$seq0, $arg_list)=&extract_val_from_arg_list(\"^\"\
,$arg_list, \"SPLICE\",\"unset\");\n($seq1, $arg_l\
ist)=&extract_val_from_arg_list(\"-seq\",$arg_list\
, \"SPLICE\",\"unset\");\n($seq2, $arg_list)=&extr\
act_val_from_arg_list(\"-in\",$arg_list, \"KEEP\",\
\"unset\");\n($seq3, $arg_list)=&extract_val_from_\
arg_list(\"-infile\",$arg_list, \"SPLICE\",\"unset\
\");\n($prf,  $arg_list)=&extract_val_from_arg_lis\
t(\"-profile\",$arg_list, \"SPLICE\",\"unset\");\n\
\n$gl{'Seq'}=$seq=&vtmpnam();#file containing all \
the sequences\n\n   #1-remove sequences from -in\n\
if ( $arg_list =~/\\-in\\b/)\n  {\n    my $save, $\
name;\n    while($arg_list=~/\\-in\\b[^-]+(\\bS[\\\
w.]+)/)\n      {\n	$name=$1;$name=~s/^.//;\n	if ( \
!-e $name){$save.=\" S$name \";}\n\n	$arg_list=~s/\
S$name/ /;\n      }\n    $arg_list=~s/\\-in\\b/\\-\
in $save /;\n  }\n   #2-prepare \n\nif (!($arg_lis\
t=~/\\-outorder/))\n  {\n    \n    $output_cl .=\"\
 -outorder=$seq\";\n  }\n@output_flag=(\"-output\"\
,\"-outfile\", \"-run_name\", \"-outorder\"); \nfo\
reach $v1 (@output_flag)\n  {\n    ($v2, $arg_list\
)=&extract_val_from_arg_list($v1,$arg_list, \"SPLI\
CE\",\"unset\");\n    if ($v2 ne \"\")\n      {\n\\
n	if ($v1 eq \"-run_name\"){$run_name=$v2;$output_\
cl .=\" $v1 $v2 \";}\n	elsif ( $v1 eq \"-outorder\\
")\n	  {\n	    if ( $v2 eq \"input\"){$v2=$seq;}\n\
	    $outorder=$v2;$output_cl .=\" $v1 $v2 \";\n	 \
 }\n	else\n	  {\n	    $output_cl .=\" $v1 $v2 \";\\
n	  }\n      }\n }\n\n\n($dpa_master_aln, $arg_lis\
t)  =&extract_val_from_arg_list(\"-dpa_master_aln\\
",$arg_list, \"SPLICE\", \"t_coffee\");\n$dpa_mast\
er_aln=~s/\\s//g;\n($nseqT, $arg_list)           =\
&extract_val_from_arg_list(\"-dpa_maxnseq\",$arg_l\
ist, \"SPLICE\", 30);\n($scoreT1, $arg_list)      \
   =&extract_val_from_arg_list(\"-dpa_min_score1\"\
,$arg_list, \"SPLICE\", 80);\n($scoreT2, $arg_list\
)         =&extract_val_from_arg_list(\"-dpa_min_s\
core2\"    ,$arg_list, \"SPLICE\", 30);\n($dpa_lim\
it, $arg_list)       =&extract_val_from_arg_list(\\
"-dpa_limit\"        ,$arg_list, \"SPLICE\", 0);\n\
($dpa_delta_id, $arg_list)    =&extract_val_from_a\
rg_list(\"-dpa_delta_id\"        ,$arg_list, \"SPL\
ICE\", 1);\n($dpa_debug, $arg_list)       =&extrac\
t_val_from_arg_list(\"-dpa_debug\"           ,$arg\
_list, \"SPLICE\", 0);\n\n\n$in_seq=$seq0.\" \".$s\
eq1.\" \".$seq2.\" \".$seq3;\n$in_prf=(($prf ne $u\
nset)?\"$prf \":\"\");\n&exit_dpa (($in_seq eq \"\\
" && $in_prf eq \"\")?1:0, \"ERROR: You did not Pr\
ovide any sequences. Use the -seq flag [FATAL: t_c\
offee_dpa]\\n\", EXIT_FAILURE);\n\n\nprint STDERR \
\"\\nSTART DPA COMPUTATION\";\n\n\n\nif ($in_seq=~\
/\\S+/)\n  {\n    \n    print STDERR \"\\n Step $s\
tep: Gather all the sequences into the tmp file: [\
$seq]\";$step++;	\n    &my_system (\"t_coffee $in_\
seq -convert -quiet -output fasta_seq -outfile=$se\
q -maxnseq 0\");\n  }\n\nif ( !-e $seq){$seq=\"\";\
}\n\nif ($in_prf=~/\\S+/)\n  {\n    $seq_in_type=\\
"profile\"; \n    $seq.= $in_prf; \n  }\nif ($seq \
eq \"\"){ &exit_dpa (1, \"\\nERROR: No Sequence FO\
und. Provide Sequences with the -seq flag [FATAL: \
t_coffee_dpa]\", EXIT_FAILURE);}\n\n \n\nif ( $run\
_name)\n  {\n    $suffix=$run_name;\n  }\nelsif ($\
in_seq=~/\\b(S[\\w.]+\\b)/)\n  {\n    my $suffix1,\
 $sufffix2;\n    $suffix1=$suffix2=$1;\n    $suffi\
x2=~s/^S//;\n    if ( -e $suffix1){$suffix=$suffix\
1;}\n    elsif ( -e $suffix2){$suffix=$suffix2;}\n\
    else\n      {\n	$suffix=&vtmpnam();	\n      }\\
n    $suffix=~s/\\.\\w+//;\n  }\n\nelse\n  {\n    \
$suffix=&vtmpnam();\n  }\n\n\nif (!$run_name){$out\
put_cl.=\" -run_name $suffix \";}\n\n\n$gl{'Tree'}\
=&seq2dpa_tree ($seq, \"$suffix.dpadnd\");\n\nprin\
t STDERR \"\\n Step $step: Prepare guide tree: $gl\
{'Tree'}\";$step++;\n\nprint STDERR \"\\n Step $st\
ep: Identify and Align Closely Related Groups\";$s\
tep++;\n%gl=&make_one_pass (0, $scoreT1,\"Align\",\
%gl);\n\nprint STDERR \"\\n Step $step: Make Multi\
ple Group Alignment\";$step++;\nwhile (!%gl ||$gl{\
'Ng'}>$nseqT)\n  {\n    %gl=&make_one_pass ($nseqT\
, $scoreT2,\"t_coffee\",%gl);\n    if ( $gl{'Newgr\
oups'}==0){$scoreT2--;}    \n  }\nprint STDERR \"\\
\n Step $step: Make The Final Alignment\";$step++;\
\n\n\n$arg_list .=$output_cl;\n\n\n%gl=&tree2group\
 (0,0, %gl);\n$gl{$gl{'0'}{'File'}}{'Output'}=\"\"\
;\n$a=0;\n&align_groups (\"t_coffee\",'0', $arg_li\
st, \" \", %gl);\n\n\n\nif ( !$dpa_keep_tmpfile){&\
clean_tmp_file (@tl);}\n\n\n\nsub seq2dpa_tree \n \
 {\n    my $seq=@_[0];\n    my $newtree=@_[1];\n  \
  my $aln=&vtmpnam ();\n\n    &my_system (\"t_coff\
ee -special_mode quickaln -in $seq -outfile $aln -\
quiet\");\n    &my_system (\"$seq_reformat -in $al\
n -action +aln2tree +tree2dpatree -output newick >\
$newtree\");\n    return $newtree;\n  }	\nsub seq2\
dpa_tree_old \n  {\n    my $aln=@_[0];\n    my $ne\
wtree=@_[1];\n    \n    \n    &my_system(\"$seq_re\
format -in $aln -action +seq2dpatree -output newic\
k > $newtree\");\n    return $newtree;\n  }\nsub a\
ln2dpa_tree \n  {\n    my $aln=@_[0];\n    my $new\
tree=&vtmpnam();\n    \n    &my_system(\"$seq_refo\
rmat -in $aln -action +aln2tree +tree2dpatree -out\
put newick > $newtree\");\n    return $newtree;\n \
 }\nsub group_file2ngroups\n  {\n    my $file=@_[0\
];\n    my $n;\n    \n    open ( F, $file);\n    w\
hile (<F>)\n      {\n	$n+=/\\>/;\n      }\n    clo\
se (F);\n    return $n;\n  }\n\nsub make_one_pass\\
n  {\n    my ($N, $ID,$pg, %gl)=@_;\n    my $a;\n\\
n    %gl=&tree2group ($N,$ID,%gl);\n    if (!$gl{'\
Newgroups'}){return %gl;}\n    else\n      {\n	for\
 ( $a=0; $a< $ng; $a++)\n	  {\n	    if ($gl{$gl{$a\
}{'File'}}{'Ng'}>1){&display_group($a, %gl);}\n	  \
  &align_groups ($pg, $a, $arg_list, \" -quiet=qui\
et \", %gl);\n	  }\n	return %gl;\n      }\n  }\n\n\
sub tree2group \n  {\n    my ($N, $ID, %gl)=@_;\n \
   my $prefix=&vtmpnam();\n    my $group_file=&vtm\
pnam();\n    my $file;\n    my $oldtree=&vtmpnam()\
;\n    my $n;\n    my $tree;\n\n\n    if ( $gl{'Ng\
'}==1){return %gl;}\n    $tree=$gl{'Tree'}; \n    \
\n    #1 extract the groups\n    &my_system (\"$se\
q_reformat -in $tree -action +tree2group $N $ID $p\
refix > $group_file\");\n    $n=group_file2ngroups\
($group_file);\n    \n    \n    $gl{'Newgroups'}=1\
;\n    if ( $n==$gl{'Ng'})\n      {\n	$gl{'Newgrou\
ps'}=0;\n	return %gl;\n      }\n    $gl{'Iteration\
'}++;\n    $gl{'MaxNseq'}=$N;$gl{'MinID'}=$ID;\n  \
  $gl{'GroupFile'}=$group_file;$gl{'Ng'}=$ng=0;\n \
   #2 Process the group list into the hash\n    op\
en (F, $group_file);\n    while (<F>)\n      {\n	$\
gl{'File'}.=$_;\n	if (/\\>/)\n	  {\n	    $line=$_;\
\n	    $line=~s/\\>//;\n	    @list=($line=~/(\\S+)\
/g);\n	    $file=$gl{$ng}{'File'}=shift @list;\n	 \
   $gl{$file}{'Output'}=$file;\n	    \n	    $gl{$f\
ile}{'Ng'}=$#list+1;\n	    if ($gl{$file}{'Ng'}>1)\
{ $gl{$file}{'Tlist'}=$gl{$file}{'Alist'}=\"(\";}\\
n	    foreach $l (@list)\n	      {\n	\n		$gl{$file\
}{'List'}.=\" $l \";\n		\n		if (!$gl{$l}{'Tlist'})\
\n		  {\n		    $gl{$l}{'Tlist'}=\"$l\";\n		    $gl\
{$l}{'Alist'}=\"$l\";\n		    $gl{$l}{'Nseq'}=1;\n	\
	    $gl{$l}{'Ng'}=1;\n		  }\n		$gl{$file}{'Tlist'\
}.=\"$gl{$l}{'Tlist'},\";\n		$gl{$file}{'Alist'}.=\
\"$gl{$l}{'Tlist'}|\";\n		$gl{$file}{'Nseq'}+=$gl{\
$l}{'Nseq'};\n	      }\n	    \n\n	    chop($gl{$fi\
le}{'Tlist'});chop($gl{$file}{'Alist'});\n	    if \
($gl{$file}{'Ng'}>1){$gl{$file}{'Tlist'}.=\")\"; $\
gl{$file}{'Alist'}.=\");\";}\n	    $ng++;\n	  }	\n\
      }\n    $gl{'Ng'}=$ng;\n    close (F);\n    \\
n    #3 Update the old tree with the new groups\n \
   $gl{'Tree'}=&vtmpnam();\n    &my_system (\"$seq\
_reformat -in $tree -action +collapse_tree $group_\
file -output newick > $gl{'Tree'}\");\n    \n    r\
eturn %gl;\n  }\n\nsub display_group \n  {\n    my\
 ($g,%gl)=@_;\n    my $f;\n    \n    if ( $g==-1)\\
n      {\n	print STDERR \"\\nIteration $gl{'Iterat\
ion'} [MaxN=$gl{'MaxNseq'}][MinID=$gl{'MinID'}]\";\
\n      }\n    else\n      {\n\n	$f=$gl{$g}{'File'\
};\n	$action=($gl{$f}{'Ng'}==1 || $gl{'Iteration'}\
==1)?\"KEEP  \":\"ALIGN \";\n        print STDERR \
\"\\n\\t[$action][MaxN=$gl{'MaxNseq'}][MinID=$gl{'\
MinID'}][File $f][Nseq=$gl{$f}{'Nseq'}][Ngroups=$g\
l{$f}{'Ng'}][$gl{$f}{'Alist'}]\";\n      }\n  }\n \
     \n\n\nsub align_groups\n  {\n    my ($pg, $g,\
 $arg, $extra_arg,%gl)=@_;\n    my $f;\n    my $Ou\
tput,$Outflag;\n    \n    \n    $f=$gl{$g}{'File'}\
;\n    $Output=($gl{$f}{'Output'});\n    \n    if \
( $pg eq \"Align\")\n      {\n	if ( !-e $f)\n	  {\\
n	    $command=\"$seq_reformat -in $gl{'Seq'}  -ac\
tion +extract_aln $gl{'GroupFile'}\";\n	    if ($g\
l{$f}{'Ng'}>1)\n	      {\n		&my_system ($command);\
\n		$command=\"t_coffee -special_mode quick_aln  S\
$f -outfile=$Output -quiet\";\n	      }\n	  }\n	el\
se \n	  {$command=\"\";}\n      }\n    elsif ( -e \
$f)\n      {	\n	$Outflag=($Output)?\"-outfile=$Out\
put\":\"\";\n	$command=\"$pg -infile $f $Outflag -\
quiet stdout $arg $extra_arg -maxnseq 0 -convert -\
quiet stdout\";\n      }\n    elsif ( $gl{$f}{'Ng'\
}==1)\n      {\n	$action=($dpa_debug)?\"cp\":\"mv\\
";\n	$command=\"$action $gl{$f}{'List'} $Output\";\
\n      }\n    else\n      {\n	$Outflag=($Output)?\
\"-outfile=$Output\":\"\";\n	$command=\"$pg -profi\
le $gl{$f}{'List'} $Outflag $arg $extra_arg -maxns\
eq 0\";\n      }\n    \n    &my_system ($command);\
\n    return $outfile;\n  }\n    \nsub my_system \\
n  {\n    my $command=@_[0];\n    my $force=@_[1];\
\n    my $status;\n\n    if ( $dpa_debug) {print S\
TDERR \"\\nCOMMAND: $command\";}\n    $status=syst\
em ($command);\n\n    if (!$force)\n       {\n	 &e\
xit_dpa (($status==1), \"Failed in Command:\\n$com\
mand\\n[FATAL: t_coffee_dpa]\\n\", EXIT_FAILURE);\\
n       }\n    \n    return $status;\n  }\n\nsub v\
tmpnam\n  {\n    my $prefix=@_[0];\n    my $tmp_fi\
le_name;\n\n    $tmp_prefix=($prefix)?$prefix:\"dp\
a_tmp_file_$$\";\n   \n    $tmp_count++;\n    $tmp\
_file_name=\"$tmp_prefix\".\"$tmp_count\";\n    $t\
l[$#tl+1]=$tmp_file_name;\n    return $tmp_file_na\
me;\n  }\n\nsub clean_tmp_file\n  {\n\n    my $lis\
t;\n    my $file;\n    \n    if ($dpa_debug){retur\
n;}\n    $list=vtmpnam();\n    `ls -1 | grep $tmp_\
prefix>$list`;\n    \n    open (F,$list);\n    whi\
le ( <F>)\n      {\n	$file=$_;\n	chop $file;\n	if \
( -e $file){unlink $file;}\n      }\n    close (F)\
;\n    unlink $list;\n  }\n\n\nsub exit_dpa\n  {\n\
  my $condition=@_[0];\n  my $error_msg=@_[1];\n  \
my $exit_value=@_[2];\n  if ( $condition)\n    {\n\
      print \"$error_msg\\n\";\n      exit ($exit_\
value);\n    }\n  else\n    {\n      return;\n    \
}\n  \n}\nsub extract_val_from_arg_list\n  {\n    \
my $arg=@_[0];\n    my $arg_list=@_[1];\n    my $k\
eep_flag=@_[2];\n    my $default_value=@_[3];\n   \
 my $val=\"\";\n    \n    #protect\n    $arg_list=\
~s/\\s-/ \\@/g;\n    $arg=~s/-/\\@/g;\n    \n    #\
search\n    if ($arg eq \"^\")\n      {\n	$arg_lis\
t=~/^([^@]*)/;\n	$val=$1;\n      }\n    else\n    \
  {$arg_list=~/$arg ([^@]*)/;$val=$1;}\n    \n    \
#remove trailing spaces\n    $val=~s/\\s*$//;\n   \
 \n    #remove the parsed sequence if needed\n    \
if (($val ne \"\") && $keep_flag ne \"KEEP\")\n   \
   {\n	if ( $arg eq \"^\"){$arg_list=~s/$val/ /;}\\
n	else {$arg_list=~s/($arg [^@]*)/ /;}\n      }\n	\
\n    #unprotect\n    $arg_list=~s/\\@/-/g;\n    $\
arg=~s/\\@/-/g;\n    \n    if (($val eq \"\") && $\
default_value ne \"unset\"){$val=$default_value;}\\
n    \n    return $val, $arg_list;\n  }\n$program=\
\"T-COFFEE (Version_8.79)\";\n\n","\n$DEBUG=1;\n$d\
pa_nseq=10;\n$dpa_sim=0;\nif (!@ARGV)\n  {\n    `t\
_coffee`;\n    exit (0);\n  }\nforeach $arg (@ARGV\
)\n  {\n    $arg_list.=\" $arg\";\n  }\n$max_nseq=\
10;\n($seq0, $arg_list)=&extract_val_from_arg_list\
(\"^\",$arg_list);\n($seq1, $arg_list)=&extract_va\
l_from_arg_list(\"-seq\",$arg_list);\n($seq2, $arg\
_list)=&extract_val_from_arg_list(\"-in\",$arg_lis\
t, \"KEEP\");\n($seq3, $arg_list)=&extract_val_fro\
m_arg_list(\"-infile\",$arg_list);\n$in_seq=$seq0.\
\" \".$seq1.\" \".$seq2.\" \".$seq3;\n\n$seq=vtmpn\
am();\n`t_coffee $in_seq -convert -output fasta_se\
q -outfile=$seq`;\n\n\n($dpa_nseq, $arg_list)=&ext\
ract_val_from_arg_list(\"-dpa_nseq\",$arg_list);\n\
($master_aln, $arg_list)=&extract_val_from_arg_lis\
t(\"-master_aln\",$arg_list);\n($sim_matrix, $arg_\
list)=&extract_val_from_arg_list(\"-sim_matrix\",$\
arg_list);\n($core_seq, $arg_list)=&extract_val_fr\
om_arg_list(\"-core_seq\",$arg_list);\n($dpa_sim, \
$arg_list)=&extract_val_from_arg_list(\"-dpa_sim\"\
,$arg_list);\n($run_name, $arg_list)=&extract_val_\
from_arg_list(\"-run_name\",$arg_list);\n($output,\
 $arg_list)=&extract_val_from_arg_list(\"-output\"\
,$arg_list);\n\n\n\nif (!$sim_mat && !$master_aln)\
#Compute the fast alignment\n  {\n    $ref_aln=vtm\
pnam();\n    `t_coffee -seq=$seq -very_fast -outfi\
le=$ref_aln -quiet`;\n    \n  }\n\nif (!$sim_mat)\\
n  {\n    $sim_mat=vtmpnam();\n    `seq_reformat -\
in $ref_aln -output sim > $sim_mat`;\n  }\n\nif ( \
!$core_seq)\n  {\n    $core_seq=vtmpnam();\n    `s\
eq_reformat -in $ref_aln -action +trimTC N$max_nse\
q -output fasta_seq > $core_seq`;\n  }\n@core_name\
=`seq_reformat -in $core_seq -output name `; \n\n@\
tot_name=`seq_reformat -in $seq -output name `;\n\\
nforeach $s (@core_name){$s=~s/\\s//g;$hcore{$s}=1\
;}\nforeach $s (@tot_name){$s=~s/\\s//g;}\nprint S\
TDERR \"T-Coffee_dpa:\\n\";\nprint STDERR \"\\tTOT\
AL  SEQ: @tot_name\\n\";\nprint STDERR \"\\tCHOSEN\
 SEQ: @core_name\\n\";\n\n\n\nopen (F, $sim_mat);\\
nwhile ( <F>)\n  {\n    @l=($_=~/(\\b[\\S]+\\b)/g)\
;\n    if (($l[0] eq \"TOP\" || $l[0] eq \"BOT\"))\
\n      {\n	$s1=$l[1];$s2=$l[2];$v=$l[3];\n	if ($h\
core{$s1} && !$hcore{$s2})\n	  {\n	    if (!$hseq{\
$s2}{\"sim\"} || $v>$hseq{$s2}{\"sim\"})\n	      {\
\n		$hseq{$s2}{\"sim\"}=$v;$hseq{$s2}{\"seq\"}=$s1\
;\n	      }\n	  }\n      }\n  }\nclose (F);\nforea\
ch $s (@tot_name)\n  {\n\n    if ( !$hseq{$s}{\"se\
q\"}){;}\n    else\n      {\n	$s2=$hseq{$s}{\"seq\\
"};\n	$v=$hseq{$s}{\"sim\"};\n		\n	if ($v>$dpa_sim\
)\n	  {\n	    $hseq{$s}{'used'}=1;\n	    $seq_list\
{$s2}{$seq_list{$s2}{'nseq'}++}=$s;\n	  }\n      }\
\n  }\nforeach $s (@core_name){$seq_list{$s}{$seq_\
list{$s}{'nseq'}++}=$s;$hseq{$s}{'used'}=1;}\nfore\
ach $s (@tot_name){if (!$hseq{$s}{'used'}){$seq_li\
st{'unused'}{$seq_list{'unused'}{'nseq'}++}=$s;}}\\
n\n\n$n=0;\nforeach $s (@core_name)\n  {\n    $ng+\
+;\n    $n=$seq_list{$s}{'nseq'};\n    for (@g_lis\
t=(), $a=0; $a<$n; $a++){@g_list=(@g_list,$seq_lis\
t{$s}{$a});}\n\n    $g_seq=vtmpnam();\n    $g_aln=\
vtmpnam();\n    \n    print STDERR \"Group $ng: $#\
g_list Seq: @g_list: \";\n    \n    \n    `seq_ref\
ormat -in $seq -action +lower +keep_name +extract_\
seq  @g_list -output fasta_seq > $g_seq`;\n    \n \
   \n    if ( $#g_list==0)\n      {\n	print STDERR\
 \"[No aln]\\n\";\n	$g_aln=$g_seq;\n      }\n    e\
lsif ($#g_list<$max_nseq) \n      {\n	print STDERR\
 \"[t_coffee]\\n\";\n	`t_coffee $g_seq -outfile=$g\
_aln -quiet $arg_list`;\n      }\n    else\n      \
{\n	print STDERR \"[t_coffee_dpa]\\n\";\n	`t_coffe\
e_dpa2 $g_seq -outfile=$g_aln $arg_list -sim_matri\
x $sim_matrix -dpa_nseq $dpa_nseq`;\n      }\n    \
@profile_list=(@profile_list, $g_aln);\n  }\n\n\np\
rint \"UNUSED $seq_list{'unused'}{'nseq'}\";\n\nif\
 ($seq_list{'unused'}{'nseq'})\n    {\n      $prf=\
vtmpnam();\n      \n      `t_coffee -profile @prof\
ile_list $arg_list -outfile=$prf -quiet`;\n      $\
n=$seq_list{\"unused\"}{'nseq'};\n      $new_seq=v\
tmpnam();\n      $new_prf=vtmpnam();\n      for ($\
a=0; $a<$n-1; $a++)\n	{\n	  $s=$seq_list{\"unused\\
"}{$a};\n	  print STDERR \"\\nADD Sequence $s\";\n\
	  \n	  `seq_reformat -in $seq -action +lower +kee\
p_name +extract_seq $s  -output fasta_seq > $new_s\
eq`;\n	  `t_coffee -profile $prf $new_seq $arg_lis\
t -outfile=$new_prf`;\n	  `cp $new_prf $prf`;\n	}\\
n      $s=$seq_list{\"unused\"}{$a};\n      `seq_r\
eformat -in $seq -action +lower +keep_name +extrac\
t_seq $s  -output fasta_seq > $new_seq`;\n      @p\
rofile_list=($prf, $new_seq);\n    }\n    \n      \
\nif ($run_name){$arg_list.=\" -run_name $run_name\
\";}\nelse \n  {\n    $in_seq=~/([\\w-]+)/;\n    $\
arg_list.=\" -run_name $1\";\n  }\nif ( $output){$\
arg_list.=\" -output $output \";}\n\n`t_coffee -pr\
ofile @profile_list $arg_list`;\n\n\n&clean (@tmp_\
file_list);\n\n\nsub vtmpnam\n  {\n    my $tmp_fil\
e_name;\n    $tmp_name_counter++;\n    $tmp_file_n\
ame=\"tmp_file_$tmp_name_counter\\_Pid$$\";\n    $\
tmp_file_list[$ntmp_file++]=$tmp_file_name;\n    r\
eturn $tmp_file_name;\n  }\nsub clean\n  {\n  my @\
fl=@_;\n  my $file;\n  return;\n\n  foreach $file \
( @fl)\n    {\n      if ( -e $file){unlink($file);\
}\n    }\n}\nsub extract_val_from_arg_list\n  {\n \
   my $arg=@_[0];\n    my $arg_list=@_[1];\n    my\
 $keep_flag=@_[2];\n    #protect\n    $arg_list=~s\
/\\s-/ \\@/g;\n    $arg=~s/-/\\@/g;\n    \n    #se\
arch\n    if ($arg eq \"^\")\n      {\n	$arg_list=\
~/^([^@]*)/;\n	$val=$1;\n      }\n    else\n      \
{$arg_list=~/$arg ([^@]*)/;$val=$1;}\n    \n    #r\
emove the parsed sequence if needed\n    if ($val \
&& $keep_flag ne \"KEEP\")\n      {\n	if ( $arg eq\
 \"^\"){$arg_list=~s/$val/ /;}\n	else {$arg_list=~\
s/($arg [^@]*)/ /;}\n      }\n	\n    #unprotect\n \
   $arg_list=~s/\\@/-/g;\n    $arg=~s/\\@/-/g;\n  \
  \n    return $val, $arg_list;\n  }\n\n","use Env\
;\nuse FileHandle;\nuse Cwd;\nuse File::Path;\nuse\
 Sys::Hostname;\n\nour $PIDCHILD;\nour $ERROR_DONE\
;\nour @TMPFILE_LIST;\nour $EXIT_FAILURE=1;\nour $\
EXIT_SUCCESS=0;\n\nour $REFDIR=getcwd;\nour $EXIT_\
SUCCESS=0;\nour $EXIT_FAILURE=1;\n\nour $PROGRAM=\\
"tc_generic_method.pl\";\nour $CL=$PROGRAM;\n\nour\
 $CLEAN_EXIT_STARTED;\nour $debug_lock=$ENV{\"DEBU\
G_LOCK\"};\nour $LOCKDIR=$ENV{\"LOCKDIR_4_TCOFFEE\\
"};\nif (!$LOCKDIR){$LOCKDIR=getcwd();}\nour $ERRO\
RDIR=$ENV{\"ERRORDIR_4_TCOFFEE\"};\nour $ERRORFILE\
=$ENV{\"ERRORFILE_4_TCOFFEE\"};\n&set_lock ($$);\n\
if (isshellpid(getppid())){lock4tc(getppid(), \"LL\
OCK\", \"LSET\", \"$$\\n\");}\n      \n\n\n\n\nour\
 $BLAST_MAX_NRUNS=2;\nour $COMMAND;\nour $PIDCHILD\
;\n\n$REF_EMAIL=\"\";\n$tmp_dir=\"\";\n$init_dir=\\
"\";\n\n\n$test=0;\nif ($test==1)\n  {\n    $SERVE\
R=\"NCBI\";\n    $query=$ARGV[0];\n    $hitf=$ARGV\
[1];\n    %s=read_fasta_seq($query);\n    @sl=keys\
(%s);\n    &blast_xml2profile (\"xx\", $s{$sl[0]}{\
seq},$maxid,$minid,$mincov, $hitf);\n    myexit ($\
EXIT_FAILURE);\n  }\n\nforeach $v(@ARGV){$cl.=\"$v\
 \";}\n$COMMAND=$cl;\n($mode)=&my_get_opt ( $cl, \\
"-mode=\",1,0);\n\n($A)=(&my_get_opt ( $cl, \"-nam\
e1=\",0,0));\n($B)=(&my_get_opt ( $cl, \"-name2=\"\
,0,0));\n($TMPDIR)=(&my_get_opt ( $cl, \"-tmpdir=\\
",0,0));\n($CACHE)=(&my_get_opt ( $cl, \"-cache=\"\
,0,0));\n($SERVER)=((&my_get_opt ( $cl, \"-server=\
\",0,0)));\n($EMAIL)=((&my_get_opt ( $cl, \"-email\
=\",0,0)));\n\nif (!$A){$A=\"A\";}\nif (!$B){$B=\"\
B\";}\n\n\nif (!$TMPDIR)\n  {\n    $HOME=$ENV{HOME\
};\n    if ($ENV{TMP_4_TCOFFEE}){$TMPDIR=$ENV{TMP_\
4_TCOFFEE};}\n    else{$TMPDIR=\"$HOME/.t_coffee/t\
mp/\";}\n  }\nif ( ! -d $TMPDIR)\n  {\n    mkdir $\
TMPDIR;\n  }\nif ( ! -d $TMPDIR)\n  {\n    print \\
"ERROR: Could not create temporary dir: $TMPDIR\\n\
\";\n    myexit ($EXIT_FAILURE);\n  }\n\n$EMAIL=~s\
/XEMAILX/\\@/g;\nif (!$EMAIL)\n  {\n    if ($ENV{E\
MAIL_4_TCOFFEE}){$EMAIL=$ENV{EMAIL_4_TCOFFEE};}\n \
   elsif ($ENV{EMAIL}){$EMAIL=$ENV{EMAIL};}\n    e\
lse {$EMAIL=$REF_EMAIL;}\n  }\n\n($maxid,$minid,$m\
incov)=(&my_get_opt ( $cl, \"-maxid=\",0,0, \"-min\
id=\",0,0,\"-mincov=\",0,0));\nif (!$cl=~/\\-maxid\
\\=/){$maxid=95;}\nif (!$cl=~/\\-minid\\=/){$minid\
=35;}\nif (!$cl=~/\\-mincov\\=/){$mincov=80;}\n\n\\
n\nif ($mode eq \"seq_msa\")\n  {\n    &seq2msa($m\
ode,&my_get_opt ( $cl, \"-infile=\",1,1, \"-method\
=\",1,2, \"-param=\",0,0, \"-outfile=\",1,0));\n  \
}\nelsif ( $mode eq \"tblastx_msa\")\n  {\n    &se\
q2tblastx_lib ($mode,&my_get_opt ( $cl, \"-infile=\
\",1,1, \"-outfile=\",1,0));\n  }\nelsif ( $mode e\
q \"tblastpx_msa\")\n  {\n    &seq2tblastpx_lib ($\
mode,&my_get_opt ( $cl, \"-infile=\",1,1, \"-outfi\
le=\",1,0));\n  }\nelsif ( $mode eq \"thread_pair\\
")\n  {\n    &seq2thread_pair($mode,&my_get_opt ( \
$cl, \"-infile=\",1,1, \"-pdbfile1=\",1,1, \"-meth\
od=\",1,2,\"-param=\",0,0, \"-outfile=\",1,0, ));\\
n  }\nelsif ( $mode eq \"pdbid_pair\")\n  {\n    &\
seq2pdbid_pair($mode,&my_get_opt ( $cl, \"-pdbfile\
1=\",1,0, \"-pdbfile2=\",1,0, \"-method=\",1,2,\"-\
param=\",0,0, \"-outfile=\",1,0, ));\n  }\nelsif (\
 $mode eq \"pdb_pair\")\n  {\n    &seq2pdb_pair($m\
ode,&my_get_opt ( $cl, \"-pdbfile1=\",1,1, \"-pdbf\
ile2=\",1,1, \"-method=\",1,2,\"-param=\",0,0, \"-\
outfile=\",1,0, ));\n  }\nelsif ( $mode eq \"profi\
le_pair\")\n  {\n     &seq2profile_pair($mode,&my_\
get_opt ( $cl, \"-profile1=\",1,1, \"-profile2=\",\
1,1, \"-method=\",1,2,\"-param=\",0,0, \"-outfile=\
\",1,0, ));\n  }\nelsif ( $mode eq \"pdb_template\\
")\n  {\n    &blast2pdb_template ($mode,&my_get_op\
t ( $cl, \"-infile=\",1,1, \"-database=\",1,0, \"-\
method=\",1,0, \"-outfile=\",1,0));\n  }\nelsif ( \
$mode eq \"profile_template\")\n  {\n    \n    &ps\
iblast2profile_template ($mode,&my_get_opt ( $cl, \
\"-infile=\",1,1, \"-database=\",1,0, \"-method=\"\
,1,0, \"-outfile=\",1,0));\n  }\nelsif ( $mode eq \
\"psiprofile_template\")\n  {\n    &psiblast2profi\
le_template ($mode,&my_get_opt ( $cl, \"-infile=\"\
,1,1, \"-database=\",1,0, \"-method=\",1,0, \"-out\
file=\",1,0));\n  }\nelsif ( $mode eq \"RNA_templa\
te\")\n  {\n    &seq2RNA_template ($mode,&my_get_o\
pt ( $cl, \"-infile=\",1,1, \"-outfile=\",1,0));\n\
  }\nelsif ( $mode eq \"tm_template\")\n  {\n    &\
seq2tm_template ($mode, \"\", &my_get_opt ( $cl, \\
"-infile=\",1,1,\"-arch=\",1,1,\"-psv=\",1,1, \"-o\
utfile=\",1,0,));\n  }\nelsif ( $mode eq \"psitm_t\
emplate\")\n  {\n    &seq2tm_template ($mode,&my_g\
et_opt ( $cl, \"-database=\",1,0, \"-infile=\",1,1\
, \"-arch=\",1,1,\"-psv=\",1,1, \"-outfile=\",1,0,\
));\n  }\nelsif ( $mode eq \"ssp_template\")\n  {\\
n    &seq2ssp_template ($mode,&my_get_opt ( $cl, \\
"-infile=\",1,1,\"-seq=\",1,1,\"-obs=\",1,1, \"-ou\
tfile=\",1,0));\n  }\nelsif ( $mode eq \"psissp_te\
mplate\")\n  {\n    &seq2ssp_template ($mode,&my_g\
et_opt ( $cl, \"-infile=\",1,1,\"-seq=\",1,1,\"-ob\
s=\",1,1, \"-outfile=\",1,0));\n  }\n\nelsif ( $mo\
de eq \"rna_pair\")\n{\n    &seq2rna_pair($mode,&m\
y_get_opt ( $cl, \"-pdbfile1=\",1,1, \"-pdbfile2=\\
",1,1, \"-method=\",1,2,\"-param=\",0,0, \"-outfil\
e=\",1,0, ));\n}\nelsif ( $mode eq \"calc_rna_temp\
late\")\n{\n    &calc_rna_template($mode,&my_get_o\
pt ( $cl, \"-infile=\",1,1,\"-pdbfile=\",1,1, \"-o\
utfile=\",1,0));\n}\nelse\n  {\n    myexit(flush_e\
rror( \"$mode is an unknown mode of tc_generic_met\
hod.pl\"));\n  }\nmyexit ($EXIT_SUCCESS);\n\n\nsub\
 seq2ssp_template\n  {\n  my ($mode, $infile,$gor_\
seq,$gor_obs,$outfile)=@_;\n  my %s, %h;\n  my $re\
sult;\n  my (@profiles);\n  &set_temporary_dir (\"\
set\",$infile,\"seq.pep\");\n  %s=read_fasta_seq (\
\"seq.pep\");\n\n  \n  open (R, \">result.aln\");\\
n  \n  #print stdout \"\\n\";\n  foreach $seq (key\
s(%s))\n    {\n      \n      open (F, \">seqfile\"\
);\n      $s{$seq}{seq}=uc$s{$seq}{seq};\n      pr\
int (F \">$s{$seq}{name}\\n$s{$seq}{seq}\\n\");\n \
     close (F);\n      $lib_name=\"$s{$seq}{name}.\
ssp\";\n      $lib_name=&clean_file_name ($lib_nam\
e);\n      \n      if ($mode eq \"ssp_template\"){\
&seq2gor_prediction ($s{$seq}{name},$s{$seq}{seq},\
 \"seqfile\", $lib_name,$gor_seq, $gor_obs);}\n   \
   elsif ($mode eq \"psissp_template\")\n	{\n	  &s\
eq2msa_gor_prediction ($s{$seq}{name},$s{$seq}{seq\
},\"seqfile\", $lib_name,$gor_seq, $gor_obs);\n	}\\
n    \n      if ( !-e $lib_name)\n	{\n	  myexit(fl\
ush_error(\"GORIV failed to compute the secondary \
structure of $s{$seq}{name}\"));\n	  myexit ($EXIT\
_FAILURE);\n	}\n      else\n	{\n	  print stdout \"\
\\tProcess: >$s{$seq}{name} _E_ $lib_name \\n\";\n\
	  print R \">$s{$seq}{name} _E_ $lib_name\\n\";\n\
	}\n      unshift (@profiles, $lib_name);\n    }\n\
  close (R);\n  &set_temporary_dir (\"unset\",$mod\
e, $method,\"result.aln\",$outfile, @profiles);\n}\
\n\nsub seq2tm_template\n  {\n  my ($mode, $db, $i\
nfile,$arch,$psv,$outfile)=@_;\n  my %s, %h;\n  my\
 $result;\n  my (@profiles);\n  &set_temporary_dir\
 (\"set\",$infile,\"seq.pep\");\n  %s=read_fasta_s\
eq (\"seq.pep\");\n\n  \n  open (R, \">result.aln\\
");\n  \n  #print stdout \"\\n\";\n  foreach $seq \
(keys(%s))\n    {\n      open (F, \">seqfile\");\n\
      print (F \">$s{$seq}{name}\\n$s{$seq}{seq}\\\
n\");\n      close (F);\n      $lib_name=\"$s{$seq\
}{name}.tmp\";\n      $lib_name=&clean_file_name (\
$lib_name);\n\n      if ($mode eq \"tm_template\")\
\n	{\n	  &safe_system (\"t_coffee -other_pg fasta_\
seq2hmmtop_fasta.pl -in=seqfile -out=$lib_name -ar\
ch=$arch -psv=$psv\");\n	}\n      elsif ( $mode eq\
 \"psitm_template\")\n	{\n	  &seq2msa_tm_predictio\
n ($s{$seq}{name},$s{$seq}{seq}, $db, \"seqfile\",\
 $lib_name,$arch, $psv);\n	}\n      if ( !-e $lib_\
name)\n	{\n	  myexit(flush_error(\"RNAplfold faile\
d to compute the secondary structure of $s{$seq}{n\
ame}\"));\n	  myexit ($EXIT_FAILURE);\n	}\n      e\
lse\n	{\n	  print stdout \"\\tProcess: >$s{$seq}{n\
ame} _T_ $lib_name\\n\";\n	  print R \">$s{$seq}{n\
ame} _T_ $lib_name\\n\";\n	}\n      unshift (@prof\
iles, $lib_name);\n    }\n  close (R);\n  &set_tem\
porary_dir (\"unset\",$mode, $method,\"result.aln\\
",$outfile, @profiles);\n}\n\nsub seq2RNA_template\
\n  {\n  my ($mode, $infile,$outfile)=@_;\n  my %s\
, %h, ;\n  my $result;\n  my (@profiles);\n  &set_\
temporary_dir (\"set\",$infile,\"seq.pep\");\n  %s\
=read_fasta_seq (\"seq.pep\");\n\n  \n  open (R, \\
">result.aln\");\n  \n  #print stdout \"\\n\";\n  \
foreach $seq (keys(%s))\n    {\n      open (F, \">\
seqfile\");\n      print (F \">$s{$seq}{name}\\n$s\
{$seq}{seq}\\n\");\n      close (F);\n      $lib_n\
ame=\"$s{$seq}{name}.rfold\";\n      $lib_name=&cl\
ean_file_name ($lib_name);\n      &safe_system (\"\
t_coffee -other_pg RNAplfold2tclib.pl -in=seqfile \
-out=$lib_name\");\n      \n      if ( !-e $lib_na\
me)\n	{\n	 myexit(flush_error(\"RNAplfold failed t\
o compute the secondary structure of $s{$seq}{name\
}\"));\n	  myexit ($EXIT_FAILURE);\n	}\n      else\
\n	{\n	  print stdout \"\\tProcess: >$s{$seq}{name\
} _F_ $lib_name\\n\";\n	  print R \">$s{$seq}{name\
} _F_ $lib_name\\n\";\n	}\n      unshift (@profile\
s, $lib_name);\n    }\n  close (R);\n  &set_tempor\
ary_dir (\"unset\",$mode, $method,\"result.aln\",$\
outfile, @profiles);\n}\n\nsub psiblast2profile_te\
mplate \n  {\n  my ($mode, $infile, $db, $method, \
$outfile)=@_;\n  my %s, %h, ;\n  my ($result,$psib\
last_output,$profile_name,@profiles);\n  my $trim=\
0;\n  &set_temporary_dir (\"set\",$infile,\"seq.pe\
p\");\n  %s=read_fasta_seq (\"seq.pep\");\n  open \
(R, \">result.aln\");\n  \n  #print stdout \"\\n\"\
;\n  foreach $seq (keys(%s))\n    {\n      open (F\
, \">seqfile\");\n      print (F \">$A\\n$s{$seq}{\
seq}\\n\");\n      close (F);\n      $psiblast_out\
put=&run_blast ($s{$seq}{name},$method, $db, \"seq\
file\",\"outfile\");\n      if ( -e $psiblast_outp\
ut)\n	{\n	  %profile=blast_xml2profile($s{$seq}{na\
me}, $s{$seq}{seq},$maxid, $minid,$mincov,$psiblas\
t_output);\n	  unlink ($psiblast_output);\n	  \n	 \
 $profile_name=\"$s{$seq}{name}.prf\";\n	  $profil\
e_name=&clean_file_name ($profile_name);\n	  unshi\
ft (@profiles, $profile_name);\n	  output_profile \
($profile_name, \\%profile, $trim);\n	  print stdo\
ut \"\\tProcess: >$s{$seq}{name} _R_ $profile_name\
 [$profile{n} Seq.] [$SERVER/blast/$db][$CACHE_STA\
TUS]\\n\";\n	  print R \">$s{$seq}{name} _R_ $prof\
ile_name\\n\";\n	}\n    }\n  close (R);\n  &set_te\
mporary_dir (\"unset\",$mode, $method,\"result.aln\
\",$outfile, @profiles);\n}\n\nsub blast2pdb_templ\
ate \n  {\n  my ($mode, $infile, $db, $method, $ou\
tfile)=@_;\n  my %s, %h, ;\n  my ($result,$blast_o\
utput);\n  &set_temporary_dir (\"set\",$infile,\"s\
eq.pep\");\n  %s=read_fasta_seq (\"seq.pep\");\n  \
open (R, \">result.aln\");\n  \n \n  #print stdout\
 \"\\n\";\n  foreach $seq (keys(%s))\n    {\n     \
 open (F, \">seqfile\");\n      print (F \">$A\\n$\
s{$seq}{seq}\\n\");\n      close (F);\n     \n    \
  $blast_output=&run_blast ($s{$seq}{name},$method\
, $db, \"seqfile\",\"outfile\");\n     \n      %p=\
blast_xml2profile($s{$seq}{name}, $s{$seq}{seq},$m\
axid, $minid,$mincov,$blast_output);\n      unlink\
 ($blast_output);\n      if ($p{n}>1)\n	{\n	  $pdb\
id=id2pdbid($p{1}{identifyer});\n	  if ( length ($\
pdbid)>5){$pdbid=id2pdbid($p{1}{definition});}\n	 \
 \n	  print R \">$s{$seq}{name} _P_ $pdbid\\n\";\n\
	  print stdout \"\\tProcess: >$s{$seq}{name} _P_ \
$pdbid [$SERVER/blast/$db][$CACHE_STATUS]\\n\";\n	\
}\n      else\n	{\n	  print R \">$s{$seq}{name}\\n\
\";\n	  print stdout \"\\tProcess: >$s{$seq}{name}\
 _P_ No Template Found [$SERVER/blast/$db][$CACHE_\
STATUS]\\n\";\n	}\n    }\n  close (R);\n  &set_tem\
porary_dir (\"unset\",$mode, $method,\"result.aln\\
",$outfile);\n}\nsub blast_msa_o2a\n  {\n    my ($\
infile,$outfile)=@_;\n    my ($a, %seq);\n    my $\
seqfile;\n    my $SEQ=new FileHandle;\n    my $seq\
file=\"seqfile\";\n   \n    open (SEQ, \">$seqfile\
\");\n    %s1=&read_fasta_seq ($infile);\n    fore\
ach $s (keys (%s1))\n      {\n	$i=$s1{$s}{order};\\
n	if ( $i==0){print SEQ \">$s1{$s}{name}\\n$s1{$s}\
{seq}\\n\";}\n	\n	$s{$i}{name}=$s;\n	$s{$i}{seq}=$\
s1{$s}{seq};\n	$s{$i}{len}=length( $s{$i}{seq});\n\
	$s{n}++;\n      }\n    close (SEQ);\n    \n    `f\
ormatdb -i $infile`;\n    `blastpgp -i $seqfile -d\
 $infile -m7 -j4 > io`;\n    &set_blast_type (\"io\
\");\n    \n    %FB=&xml2tag_list (\"io\", \"Blast\
Output\");\n    \n    open (F, \">$outfile\");\n  \
  print F \"! TC_LIB_FORMAT_01\\n\";\n    print F \
\"$s{n}\\n\";\n    for ( $a=0; $a<$s{n}; $a++)\n  \
    {\n	print F \"$s{$a}{name} $s{$a}{len} $s{$a}{\
seq}\\n\";\n      }\n\n\n    for ( $a=0; $a<$FB{n}\
; $a++)\n      {\n	%p=blast_xml2profile ($s{$a}{na\
me}, $s{$a}{seq},100, 0, 0, $FB{$a}{body});\n	for \
($b=1; $b<$p{n}; $b++)\n	  {\n	    my $l=length ($\
p{$b}{Qseq});\n	    my $hit=$p{$b}{definition};\n	\
    my $Qstart=$p{$b}{Qstart};\n	    my $Hstart=$p\
{$b}{Hstart};\n	    my $identity=$p{$b}{identity};\
\n	    my @lrQ=split (//,$p{$b}{Qseq});\n	    my @\
lrH=split (//,$p{$b}{Hseq});\n	    my $i= $s1{$s{$\
a}{name}}{order}+1;\n	    my $j= $s1{$hit}{order}+\
1;\n	    #if ( $j==$i){next;}\n	    printf F \"# %\
d %d\\n\", $i, $j;\n	    #  print  F \"\\n$p{$b}{Q\
seq} ($Qstart)\\n$p{$b}{Hseq} ($Hstart)\";\n	    f\
or ($c=0; $c<$l; $c++)\n	      {\n		my $rQ=$lrQ[$c\
];\n		my $rH=$lrH[$c];\n		my $n=0;\n		\n		if ($rQ \
ne \"-\"){$n++, $Qstart++;}\n		if ($rH ne \"-\"){$\
n++; $Hstart++;}\n		\n		if ( $n==2)\n		  {\n		    \
printf F \"\\t%d %d %d\\n\", $Qstart-1, $Hstart-1,\
$identity;\n		  }\n	      }\n	  }\n      }\n    pr\
int F \"! SEQ_1_TO_N\\n\";\n    close (F);\n    re\
turn $output;\n  \n  }\n\nsub blast_msa\n  {\n    \
my ($infile,$outfile)=@_;\n    my ($a, %seq);\n   \
 %s1=&read_fasta_seq ($infile);\n    foreach $s (k\
eys (%s1))\n      {\n	$i=$s1{$s}{order};\n	$s{$i}{\
name}=$s;\n	$s{$i}{seq}=$s1{$s}{seq};\n	$s{$i}{len\
}=length( $s{$i}{seq});\n	$s{n}++;\n      }\n    `\
formatdb -i $infile`;\n    `blastpgp -i $infile -d\
 $infile -m7 -j4 > io`;\n    &set_blast_type (\"io\
\");\n    \n    %FB=&xml2tag_list (\"io\", \"Blast\
Output\");\n    \n    open (F, \">$outfile\");\n  \
  print F \"! TC_LIB_FORMAT_01\\n\";\n    print F \
\"$s{n}\\n\";\n    for ( $a=0; $a<$s{n}; $a++)\n  \
    {\n	print F \"$s{$a}{name} $s{$a}{len} $s{$a}{\
seq}\\n\";\n      }\n    for ( $a=0; $a<$FB{n}; $a\
++)\n      {\n	%p=blast_xml2profile ($s{$a}{name},\
 $s{$a}{seq},100, 0, 0, $FB{$a}{body});\n	for ($b=\
1; $b<$p{n}; $b++)\n	  {\n	    my $l=length ($p{$b\
}{Qseq});\n	    my $hit=$p{$b}{definition};\n	    \
my $Qstart=$p{$b}{Qstart};\n	    my $Hstart=$p{$b}\
{Hstart};\n	    my $identity=$p{$b}{identity};\n	 \
   my @lrQ=split (//,$p{$b}{Qseq});\n	    my @lrH=\
split (//,$p{$b}{Hseq});\n	    my $i= $s1{$s{$a}{n\
ame}}{order}+1;\n	    my $j= $s1{$hit}{order}+1;\n\
	    #if ( $j==$i){next;}\n	    printf F \"# %d %d\
\\n\", $i, $j;\n	    #  print  F \"\\n$p{$b}{Qseq}\
 ($Qstart)\\n$p{$b}{Hseq} ($Hstart)\";\n	    for (\
$c=0; $c<$l; $c++)\n	      {\n		my $rQ=$lrQ[$c];\n\
		my $rH=$lrH[$c];\n		my $n=0;\n		\n		if ($rQ ne \\
"-\"){$n++, $Qstart++;}\n		if ($rH ne \"-\"){$n++;\
 $Hstart++;}\n		\n		if ( $n==2)\n		  {\n		    prin\
tf F \"\\t%d %d %d\\n\", $Qstart-1, $Hstart-1,$ide\
ntity;\n		  }\n	      }\n	  }\n      }\n    print \
F \"! SEQ_1_TO_N\\n\";\n    close (F);\n    return\
 $output;\n  \n  }\n\nsub seq2msa\n  {\n    my ($m\
ode, $infile, $method, $param, $outfile)=@_;\n    \
&set_temporary_dir (\"set\",$infile,\"seq.pep\");\\
n    $param.=\" >/dev/null 2>&1 \";\n    \n    #ma\
ke sure test.pep is in FASTA\n    &safe_system (\"\
t_coffee -other_pg seq_reformat -in seq.pep -outpu\
t fasta_seq > x\");\n    `mv x seq.pep`;\n    \n  \
  if ( $method eq \"blastpgp\")\n      {\n	if ( $p\
aram =~/o2a/){&blast_msa_o2a (\"seq.pep\", \"resul\
t.aln\");}\n	else {&blast_msa (\"seq.pep\", \"resu\
lt.aln\");}\n      }\n    elsif ( $method eq \"mus\
cle\")\n      {\n	`muscle -in seq.pep -out result.\
aln $param`;\n      }\n    elsif ( $method eq \"pr\
obcons\")\n      {\n	`probcons seq.pep >result.aln\
 2>/dev/null`;\n      }\n    elsif ( $method eq \"\
mafft\")\n      {\n	`mafft --quiet --localpair --m\
axiterate 1000 seq.pep> result.aln  2>/dev/null`\n\
      }\n    elsif ( $method=~/prank/)\n      {\n	\
`$method -d=seq.pep -o=result.aln -quiet 2>/dev/nu\
ll`;\n	`mv result.aln.1.fas result.aln`;\n      }\\
n    else\n      {\n	`$method -infile=seq.pep -out\
file=result.aln`;\n      }\n    \n    &set_tempora\
ry_dir (\"unset\",$mode, $method,\"result.aln\",$o\
utfile);\n    myexit ($EXIT_SUCCESS);\n  }\n\nsub \
seq2thread_pair\n  {\n    my ($mode, $infile, $pdb\
file1, $method, $param, $outfile)=@_;\n    &set_te\
mporary_dir (\"set\",$infile,\"seq.pep\",$pdbfile1\
,\"struc.pdb\");\n    if ($method eq \"fugueali\")\
\n      {\n	#Env Variable that need to be defined \
for Fugue\n	if (!$ENV{FUGUE_LIB_LIST}){$ENV{FUGUE_\
LIB_LIST}=\"DUMMY\";}\n	if (!$ENV{HOMSTRAD_PATH}) \
 {$ENV{HOMSTRAD_PATH}=\"DUMMY\";}\n	if (!$ENV{HOMS\
_PATH}){$ENV{HOMS_PATH}=\"DUMMY\";}\n	\n	`joy stru\
c.pdb >x 2>x`;\n	&check_file(\"struc.tem\", \"Joy \
failed [FATAL:$PROGRAM/$method]\");\n	`melody -t s\
truc.tem >x 2>x`;\n	&check_file(\"struc.tem\", \"M\
elody failed [FATAL:$PROGRAM/$method]\");\n	`fugue\
ali -seq seq.pep -prf struc.fug -print > tmp_resul\
t.aln`;\n	\n	&check_file(\"tmp_result.aln\", \"Fug\
ue failed [FATAL:$PROGRAM/$method]\");\n	&safe_sys\
tem (\"t_coffee -other_pg seq_reformat -in tmp_res\
ult.aln -output fasta_aln >result.aln\");\n      }\
\n    elsif ( $method eq \"t_coffee\")\n      {\n	\
&safe_system (\"t_coffee -in Pstruc.pdb Sseq.pep M\
slow_pair -outfile result.aln -quiet\");\n      }\\
n    else\n      {\n	&safe_system (\"$method -infi\
le=seq.pep -pdbfile1=struc.pdb -outfile=result.aln\
 $param>x 2>x\");\n      }\n    &set_temporary_dir\
 (\"unset\",$mode,$method,\"result.aln\",$outfile)\
;\n    myexit ($EXIT_SUCCESS);\n  }\nsub seq2pdbid\
_pair\n  {\n    my ($mode, $pdbfile1, $pdbfile2, $\
method, $param, $outfile)=@_;\n    my ($name);\n\n\
    \n    &set_temporary_dir (\"set\");\n    $name\
=$pdbfile1.\" \".$pdbfile2;\n\n    if (    &cache_\
file(\"GET\",\"\",\"$name\",\"$method\",\"dali\",$\
outfile,\"EBI\"))\n      {return $outfile;}\n    e\
lse\n      {\n	if ($method eq \"daliweb\")\n	  {\n\
	    $pdbfile1=~/(....)(.)/;\n	    $id1=$1; $c1=$2\
;\n	    \n	    $pdbfile2=~/(....)(.)/;\n	    $id2=\
$1; $c2=$2;\n	    \n	    $command=\"t_coffee -othe\
r_pg dalilite.pl --pdb1 $id1 --chainid1 $c1 --pdb2\
 $id2 --chainid2 $c2 --email=$EMAIL  >dali_stderr \
2>dali_stderr\";\n	    $dali=`$command`;\n	    \n	\
    open (F, \"dali_stderr\");\n	    while (<F>)\n\
	      {\n		if ( /JobId: dalilite-(\\S+)/)\n		{\n	\
	  $jobid=$1;\n		}\n	      }\n	    close (F);\n	  \
  unlink (\"dali_stderr\");\n	    \n	    $output1=\
\"dalilite-$jobid.txt\";\n	    if ( -e $output1)\n\
	      {\n		unlink ($output1);\n		&url2file (\"htt\
p://www.ebi.ac.uk/Tools/es/cgi-bin/jobresults.cgi/\
dalilite/dalilite-$jobid/aln.html\", \"output2\");\
\n		\n		if ( -e \"output2\")\n		  {\n		    my ($se\
q1, $seq2);\n		    $seq1=$seq2=\"\";\n		    \n		  \
  open (F, \"output2\");\n		    while (<F>)\n		   \
   {\n			$l=$_;\n			if ( $l=~/Query\\s+(\\S+)/)\n	\
		  {\n			    $seq1.=$1;\n			  }\n			elsif ( $l=~/\
Sbjct\\s+(\\S+)/)\n			  {\n			    $seq2.=$1;\n			 \
 }\n		      }\n		    close (F);\n		    unlink (\"o\
utput2\");\n		    if ($seq1 ne \"\" && $seq2 ne \"\
\")\n		      {\n			$output3=\">$A\\n$seq1\\n>$B\\n\
$seq2\\n\";\n			$output3=~s/\\./-/g;\n			open (F, \
\">result.aln\");\n			print F \"$output3\";\n			cl\
ose (F);\n		      }\n		  }\n	      }\n	  }\n      \
}\n    &cache_file(\"SET\",\"\",\"$name\",\"$metho\
d\",\"dali\",\"result.aln\",\"EBI\");\n    &set_te\
mporary_dir (\"unset\",$mode, $method, \"result.al\
n\",$outfile);\n    myexit ($EXIT_SUCCESS);\n  }\n\
sub seq2pdb_pair\n  {\n    my ($mode, $pdbfile1, $\
pdbfile2, $method, $param, $outfile)=@_;\n    \n  \
  &set_temporary_dir (\"set\",$pdbfile1,\"pdb1.pdb\
\",$pdbfile2,\"pdb2.pdb\");\n    if ($method eq \"\
t_coffee\")\n      {\n	&safe_system (\"t_coffee -i\
n Ppdb1.pdb Ppdb2.pdb -quiet -outfile=result.aln\"\
);\n      }\n    elsif ( $method eq \"DaliLite\")\\
n      {\n	if ( &safe_system (\"DaliLite -pairwise\
 pdb1.pdb pdb2.pdb >tmp1\")==$EXIT_SUCCESS)\n	  {\\
n	     my ($seq1, $seq2);\n	     $seq1=$seq2=\"\";\
\n		    \n	     open (F, \"tmp1\");\n	     while (\
<F>)\n	       {\n		 $l=$_;\n		 if ( $l=~/Query\\s+\
(\\S+)/)\n		   {\n		     $seq1.=$1;\n		   }\n		 el\
sif ( $l=~/Sbjct\\s+(\\S+)/)\n		   {\n		     $seq2\
.=$1;\n		   }\n	       }\n	     close (F);\n	     \
unlink (\"tmp1\");\n	     if ($seq1 ne \"\" && $se\
q2 ne \"\")\n	       {\n		 my $output3=\">$A\\n$se\
q1\\n>$B\\n$seq2\\n\";\n		 $output3=~s/\\./-/g;\n	\
	 open (F, \">result.aln\");\n		 print F \"$output\
3\";\n		 close (F);\n	       }\n	   }\n	else\n	  {\
\n	    print \"ERROR: DalLite failed to align the \
considered structures[tc_generic_method.pl]\\n\";\\
n	  }    \n      }\n    elsif ( $method eq \"TMali\
gn\")\n      {\n	if ( &safe_system (\"TMalign pdb1\
.pdb pdb2.pdb >tmp1\")==$EXIT_SUCCESS)\n	  {\n	   \
 `tail -4 tmp1 > tmp2`;\n	    \n	    open (F, \"tm\
p2\");\n	    while (<F>)\n	      {\n		unshift(@l, \
$_);\n	      }\n	    close (F);\n	    open (F, \">\
result.aln\");\n	    $l[3]=~s/[^a-zA-Z0-9-]/\\-/g;\
\n	    $l[1]=~s/[^a-zA-Z0-9-]/\\-/g;\n	    print F\
 \">$A\\n$l[3]\\n>$B\\n$l[1]\\n\";\n	    close (F)\
;\n	  }\n	else\n	  {\n	    print \"ERROR: TMalign \
failed to align the considered structures[tc_gener\
ic_method.pl]\\n\";\n	    `rm result.aln >/dev/nul\
l 2>/dev/null`;\n	  }\n      }\n    elsif ( $metho\
d eq \"mustang\")\n      {\n	if ( &safe_system (\"\
mustang -i pdb1.pdb pdb2.pdb -F fasta >/dev/null 2\
>/dev/null\")==$EXIT_SUCCESS)\n	  {\n	    `mv resu\
lts.afasta result.aln`;\n	  }\n	else\n	  {\n	    p\
rint \"ERROR: mustang failed to align the consider\
ed structures[tc_generic_method.pl]\\n\";\n	    `r\
m result.aln >/dev/null 2>/dev/null`;\n	  }\n     \
 }\n    else\n      {\n	if ( &safe_system (\"$meth\
od -pdbfile1=pdb1.pep -pdbfile2=pdb2.pdb -outfile=\
result.aln $param>x 2>x\")==$EXIT_SUCCESS)\n	  {\n\
	    `mv results.afasta result.aln`;\n	  }\n	else\\
n	  {\n	    print \"ERROR: $method failed to align\
 the considered structures[tc_generic_method.pl]\\\
n\";\n	    `rm result.aln >/dev/null 2>/dev/null`;\
\n	  }\n      }\n    &set_temporary_dir (\"unset\"\
,$mode, $method, \"result.aln\",$outfile);\n    my\
exit ($EXIT_SUCCESS);\n  }\n\nsub seq2profile_pair\
\n  {\n    my ($mode, $profile1, $profile2, $metho\
d, $param, $outfile)=@_;\n    \n    \n    if ($met\
hod eq \"clustalw\")\n      {\n	&set_temporary_dir\
 (\"set\",$profile1,\"prf1.aln\",$profile2,\"prf2.\
aln\");\n	`clustalw -profile1=prf1.aln -profile2=p\
rf2.aln -outfile=result.aln`;\n	&set_temporary_dir\
 (\"unset\",$mode, $method, \"result.aln\",$outfil\
e);\n      }\n    elsif ( $method eq \"hhalign\")\\
n      {\n	hhalign ( $profile1,$profile2,$outfile,\
$param);\n      }\n    else\n      {\n	\n	`$method\
 -profile1=prf1.aln -profile2=prf2.aln -outfile=re\
sult.aln $param>x 2>x`;\n      }\n   \n    myexit \
($EXIT_SUCCESS);\n  }\n\nsub pg_is_installed\n  {\\
n    my @ml=@_;\n    my $r, $p, $m;\n    my $suppo\
rted=0;\n    \n    my $p=shift (@ml);\n    if ($p=\
~/::/)\n      {\n	if (safe_system (\"perl -M$p -e \
1\")==$EXIT_SUCCESS){return 1;}\n	else {return 0;}\
\n      }\n    else\n      {\n	$r=`which $p 2>/dev\
/null`;\n	if ($r eq \"\"){return 0;}\n	else {retur\
n 1;}\n      }\n  }\n\n\n\nsub check_internet_conn\
ection\n  {\n    my $internet;\n    my $tmp;\n    \
&check_configuration ( \"wget\"); \n    \n    $tmp\
=&vtmpnam ();\n    \n    if     (&pg_is_installed \
   (\"wget\")){`wget www.google.com -O$tmp >/dev/n\
ull 2>/dev/null`;}\n    elsif  (&pg_is_installed  \
  (\"curl\")){`curl www.google.com -o$tmp >/dev/nu\
ll 2>/dev/null`;}\n    \n    if ( !-e $tmp || -s $\
tmp < 10){$internet=0;}\n    else {$internet=1;}\n\
    if (-e $tmp){unlink $tmp;}\n\n    return $inte\
rnet;\n  }\nsub check_pg_is_installed\n  {\n    my\
 @ml=@_;\n    my $r=&pg_is_installed (@ml);\n    i\
f (!$r && $p=~/::/)\n      {\n	print STDERR \"\\nY\
ou Must Install the perl package $p on your system\
.\\nRUN:\\n\\tsudo perl -MCPAN -e 'install $pg'\\n\
\";\n      }\n    elsif (!$r)\n      {\n	myexit(fl\
ush_error(\"\\nProgram $p Supported but Not Instal\
led on your system\"));\n      }\n    else\n      \
{\n	return 1;\n      }\n  }\nsub set_temporary_dir\
\n  {\n    my @list=@_;\n    my $dir_mode, $a, $mo\
de, $method;\n  \n    $dir_mode=shift (@list);\n\n\
    \n    if ( $dir_mode eq \"set\")\n      {\n	$i\
nitial_dir=cwd();\n	if ( !$tmp_dir)\n	  {\n	    $r\
and=rand (100000);\n	    $tmp_dir=\"$TMPDIR/tmp4tc\
offee_profile_pair_dir_$$\\_P_$rand\";\n	  }\n	if \
( !-d $tmp_dir)\n	  {\n	    push (@TMPDIR_LIST, $t\
mp_dir);\n	    `mkdir $tmp_dir`;\n	  }\n	\n	for ( \
$a=0; $a<=$#list; $a+=2)\n	      {\n		`cp $list[$a\
] $tmp_dir/$list[$a+1]`;\n	      }\n	chdir $tmp_di\
r;\n      }\n    elsif ( $dir_mode eq \"unset\")\n\
      {\n	$mode=shift (@list);\n	$method=shift (@l\
ist);\n	\n	if (!-e $list[0])\n	  {\n	   myexit(flu\
sh_error(\"Program $method failed to produce $list\
[1]\" ));\n	    myexit ($EXIT_FAILURE);\n	  }\n	el\
se\n	  {\n	    chdir $initial_dir;\n	    # `t_coff\
ee -other_pg seq_reformat -in $tmp_dir/$list[0] -o\
utput fasta_aln -out $tmp_dir/result2.aln`;\n	    \
`cp $tmp_dir/$list[0] $tmp_dir/result2.aln`;\n	   \
 if ( $list[1] eq \"stdout\")\n	      {\n		open (F\
, \"$tmp_dir/result2.aln\");\n		while (<F>){print \
$_;}close(F);\n	      }\n	    else\n	      {\n		`m\
v $tmp_dir/result2.aln $list[1]`;\n	      }\n	    \
shift (@list); shift (@list);\n	    foreach $f (@l\
ist)\n	      {\n		`mv $tmp_dir/$f .`;\n	      }\n	\
  }\n      }\n  }\n\n\n\n\nsub my_get_opt\n  {\n  \
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
=1 && $argv eq \"\")\n	  {\n	    myexit(flush_erro\
r( \"ERROR: Option $option must be set\"));\n	    \
myexit ($EXIT_FAILURE);\n	  }\n	if ($status==0){;}\
\n	elsif ($status ==1 && $argv ne \"\" && !-e $arg\
v)\n	  {\n	    myexit(flush_error( \"File $argv mu\
st exist\"));\n	    myexit ($EXIT_FAILURE);\n	  }\\
n	elsif ( $status==2 && $argv ne \"\" && &check_pg\
_is_installed ($argv)==0)\n	  {\n	    myexit(flush\
_error( \" $argv is not installed\"));\n	    myexi\
t ($EXIT_FAILURE);\n	  }\n      }\n\n    return @a\
rgl;\n    }\n\nsub check_file \n  {\n    my ($file\
, $msg)=@_;\n\n    if ( !-e $file)\n      {\n	myex\
it(flush_error(\"$msg\"));\n      }\n    }\nsub hh\
align\n  {\n    my ($aln1, $aln2, $outfile, $param\
)=@_;\n    my $h1, $h2;\n    \n    $h{0}{index}=0;\
\n    $h{1}{index}=1;\n    \n    $h{0}{aln}=$aln1;\
\n    $h{1}{aln}=$aln2;\n\n   \n\n    %{$h{0}}=aln\
2psi_profile (%{$h{0}});\n    %{$h{1}}=aln2psi_pro\
file (%{$h{1}});\n\n    $param=~s/#S/ /g;\n    $pa\
ram=~s/#M/\\-/g;\n    $param=~s/#E/\\=/g;\n    \n\\
n    \n    $command=\"hhalign -i $h{0}{a3m} -t $h{\
1}{a3m} -tc $outfile.tmp -rank 1 -mapt 0 $param\";\
\n    `$command`;\n    \n  #  `hhalign -i $h{0}{a3\
m} -t $h{1}{a3m} -tc $outfile.tmp -rank 1 -mapt 0 \
-gapf 0.8 -gapg 0.8`;\n    \n\n    # To run global\
 use the following\n    \n    open (I, \"$outfile.\
tmp\");\n    open (O, \">$outfile\");\n    $h{0}{c\
ons}=s/\\./x/g;\n    $h{1}{cons}=s/\\./x/g;\n\n   \
 print O \"! TC_LIB_FORMAT_01\\n2\\n$h{0}{name} $h\
{0}{len} $h{0}{seq}\\n$h{1}{name} $h{1}{len} $h{1}\
{seq}\\n#1 2\\n\";\n    \n    while (<I>)\n      {\
\n	if (/(\\d+)\\s+(\\d+)\\s+(\\d+)/)\n	  {\n	    p\
rint O \"\\t$h{0}{$1}\\t$h{1}{$2}\\t$3\\n\";\n	  }\
\n      }\n    print O \"! SEQ_1_TO_N\\n\";\n\n   \
 close (O);\n    close (I);\n  }\n\nsub aln2psi_pr\
ofile\n  {\n    my (%h)=@_;\n    my ($aln,$i,$hv, \
$a, @c, $n);\n   \n    $i=$h{index};\n    $aln=$h{\
aln};\n\n    `cp $aln $$.hhh_aln`;\n    $command=\\
"t_coffee -other_pg seq_reformat -in $aln -output \
hasch\";\n    $hv=`$command`;chomp ($hv);\n    \n \
   $h{a2m}=\"$tmp/$hv.tmp4hhpred.a2m\";\n    $h{a3\
m}=\"$tmp/$hv.tmp4hhpred.a3m\";\n    if ( -e $h{a3\
m}){;}\n    else\n      {\n	`hhconsensus  -M 50 -i\
 $h{aln} -oa2m $h{a2m}`;\n	if (!-e $h{a2m})\n	  {\\
n	    print STDERR \"Program tc_generic_method.pl \
FAILED to run:\\n\\thhconsensus  -M 50 -i $h{aln} \
-oa2m $h{a2m}\";\n	    myexit ($EXIT_FAILURE);\n	 \
 }\n	\n	`hhconsensus  -M 50 -i $h{aln} -oa3m $h{a3\
m}`;\n	if (!-e $h{a3m})\n	  {\n	    print STDERR \\
"Program tc_generic_method.pl FAILED to run:\\n\\t\
hhconsensus  -M 50 -i $h{aln} -oa3m $h{a3m}\";\n	 \
   myexit ($EXIT_FAILURE);\n	  }\n       `buildali\
.pl $h{a3m} -n 1`;\n      }\n    \n    \n    $h{a2\
m_seq}=`head -n 2 $h{a2m} | grep -v \">\"`;chomp (\
$h{a2m_seq});\n    $h{a3m_seq}=`head -n 2 $h{a3m} \
| grep -v \">\"`;chomp ($h{a3m_seq});\n    $h{cons\
}=$h{a2m_seq};\n    $h{seq}=`head -n 2 $h{aln} | g\
rep -v \">\"`;chomp ($h{seq});\n    \n    \n\n    \
@c=split (//, $h{cons});\n    $h{len}=$#c+1;\n    \
for ($n=0,$a=0, $b=0; $a<$h{len};$a++)\n      {\n	\
if ( $c[$a]=~/[A-Z]/)\n	  {\n	    $h{++$n}=++$b;\n\
\n	  }\n	elsif ( $c[$a]=~/[a-z\\.]/)\n	  {\n	    +\
+$b;\n	  }\n      }\n    \n    $name=`head -n 2 $h\
{aln} | grep \">\"`;\n    $name=~/\\>(\\S+)/;\n   \
 $h{name}=$1;\n    \n    `cp $h{a2m} $i.a2m`;\n   \
 `cp $h{a3m} $i.a3m`;\n    `cp $h{aln} $i.hh_aln`;\
\n    \n    return %h;\n  }\n\nsub read_fasta_seq \
\n  {\n    my $f=@_[0];\n    my %hseq;\n    my (@s\
eq, @com, @name);\n    my ($a, $s,$nseq);\n\n    o\
pen (F, $f);\n    while (<F>)\n      {\n	$s.=$_;\n\
      }\n    close (F);\n\n    \n    @name=($s=~/>\
(\\S*).*\\n[^>]*/g);\n    \n    @seq =($s=~/>.*.*\\
\n([^>]*)/g);\n    @com =($s=~/>\\S*(.*)\\n([^>]*)\
/g);\n\n    \n    $nseq=$#name+1;\n    \n    for (\
$a=0; $a<$nseq; $a++)\n      {\n	my $s;\n	my $n=$n\
ame[$a];\n	$hseq{$n}{name}=$n;\n	$seq[$a]=~s/[^A-Z\
a-z]//g;\n	$hseq{$n}{order}=$a;\n	$hseq{$n}{seq}=$\
seq[$a];\n	$hseq{$n}{com}=$com[$a];\n	\n      }\n \
   return %hseq;\n  }\n\nsub file_contains \n  {\n\
    my ($file, $tag, $max)=(@_);\n    my ($n);\n  \
  $n=0;\n    \n    if ( !-e $file && ($file =~/$ta\
g/)) {return 1;}\n    elsif ( !-e $file){return 0;\
}\n    else \n      {\n	open (FC, \"$file\");\n	wh\
ile ( <FC>)\n	  {\n	    if ( ($_=~/$tag/))\n	     \
 {\n		close (FC);\n		return 1;\n	      }\n	    els\
if ($max && $n>$max)\n	      {\n		close (FC);\n		r\
eturn 0;\n	      }\n	    $n++;\n	  }\n      }\n   \
 close (FC);\n    return 0;\n  }\n	    \n	  \nsub \
file2string\n  {\n    my $f=@_[0];\n    my $string\
, $l;\n    open (F,\"$f\");\n    while (<F>)\n    \
  {\n\n	$l=$_;\n	#chomp ($l);\n	$string.=$l;\n    \
  }\n    close (F);\n    $string=~s/\\r\\n//g;\n  \
  $string=~s/\\n//g;\n    return $string;\n  }\n\n\
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
{\n\n	    myexit(flush_error(\"Option $option must\
 be set\"));\n	   \n	  }\n	if ($status==0){;}\n	el\
sif ($status ==1 && $argv ne \"\" && !-e $argv)\n	\
  {\n	     myexit(flush_error(\"File $argv must ex\
ist\"));\n	   \n	  }\n	elsif ( $status==2 && $argv\
 ne \"\" && &check_pg_is_installed ($argv)==0)\n	 \
 {\n	    myexit(flush_error(\"$argv is not install\
ed\"));\n	   \n	  }\n      }\n\n    return @argl;\\
n    }\n\nsub tag2value \n  {\n    \n    my $tag=(\
@_[0]);\n    my $word=(@_[1]);\n    my $return;\n \
   \n    $tag=~/$word=\"([^\"]+)\"/;\n    $return=\
$1;\n    return $return;\n  }\n      \nsub hit_tag\
2pdbid\n  {\n    my $tag=(@_[0]);\n    my $pdbid;\\
n       \n    $tag=~/id=\"(\\S+)\"/;\n    $pdbid=$\
1;\n    $pdbid=~s/_//;\n    return $pdbid;\n  }\ns\
ub id2pdbid \n  {\n    my $in=@_[0];\n    my $id;\\
n    \n    $in=~/(\\S+)/;\n    $id=$in;\n    \n   \
 if ($id =~/pdb/)\n      {\n	$id=~/pdb(.*)/;\n	$id\
=$1;\n      }\n    $id=~s/[|��_]//g;\n    retu\
rn $id;\n  }\nsub set_blast_type \n  {\n    my $fi\
le =@_[0];\n    if (&file_contains ($file,\"EBIApp\
licationResult\",100)){$BLAST_TYPE=\"EBI\";}\n    \
elsif (&file_contains ($file,\"NCBI_BlastOutput\",\
100)) {$BLAST_TYPE=\"NCBI\";}\n    else\n      {\n\
	$BLAST_TYPE=\"\";\n      }\n    return $BLAST_TYP\
E;\n  }\nsub is_valid_blast_xml\n    {\n      my $\
file=shift;\n      my $line;\n      \n      \n    \
  if ( !-e $file) {return 0;}\n      $line=&file2t\
ail ($file,100);\n      if ( $line=~/<\\/EBIApplic\
ationResult/ || $line=~/<\\/NCBI_BlastOutput/ || $\
line=~/<\\/BlastOutput/ ){return 1;}\n      return\
 0;\n    }\nsub file2blast_flavor\n      {\n	my $f\
ile=shift;\n	if (&file_contains ($file,\"EBIApplic\
ationResult\",100)){return \"EBI\";}\n	elsif (&fil\
e_contains ($file,\"EBIApplicationResult\",100)){r\
eturn \"NCBI\";}\n	else {return \"UNKNOWN\";}\n   \
   }\nsub blast_xml2profile \n  {\n    my ($name,$\
seq,$maxid, $minid, $mincov, $file)=(@_);\n    my \
(%p, $a, $string, $n);\n    \n\n\n    if ($BLAST_T\
YPE eq \"EBI\" || &file_contains ($file,\"EBIAppli\
cationResult\",100)){%p=ebi_blast_xml2profile(@_);\
}\n    elsif ($BLAST_TYPE eq \"NCBI\" || &file_con\
tains ($file,\"NCBI_BlastOutput\",100)){%p=ncbi_bl\
ast_xml2profile(@_);}\n    else \n      {\n	myexit\
(add_error ( $$,$$,getppid(), \"BLAST_FAILURE::unk\
own XML\",$CL));\n      }\n    for ($a=0; $a<$p{n}\
; $a++)\n      {\n	my $name=$p{$a}{name};\n	$p{$na\
me}{seq}=$p{$a}{seq};\n	$p{$name}{index}=$a;\n    \
  }\n    return %p;\n  }\nsub ncbi_tblastx_xml2lib\
_file \n  {\n    my  ($outlib,$string)=(@_);\n    \
my ($L,$l, $a,$b,$c,$d,$i,$nhits,@identifyerL);\n \
   my (%ITERATION);\n      \n    open (F, \">>$out\
lib\");\n    \n    $seq=~s/[^a-zA-Z]//g;\n    $L=l\
ength ($seq);\n    \n    %ITERATION=xml2tag_list (\
$string, \"Iteration\");\n    for ($i=0; $i<$ITERA\
TION{n};$i++)\n      {\n	my ($qindex, $qlen, %hit,\
 $string);\n	$string=$ITERATION{$i}{body};\n\n	$qi\
ndex=xmltag2value($string,\"Iteration_iter-num\");\
\n	$qlen  =xmltag2value($string,\"Iteration_query-\
len\");\n	%hit=&xml2tag_list  ($string, \"Hit\");\\
n\n	for ($a=0; $a<$hit{n}; $a++)\n	  {\n	    my ($\
string);\n	    $string=$hit{$a}{body};\n	 \n	    $\
hindex=xmltag2value($string,\"Hit_accession\")+1;\\
n	    if ($hindex<=$qindex){next;}\n	    else  {pr\
int F  \"# $qindex $hindex\\n\";}\n		   \n	   \n	 \
   $hlen=xmltag2value  ($string,\"Hit_len\");\n	  \
  %HSP=&xml2tag_list  ($string, \"Hsp\");\n	   \n	\
    for ($b=0; $b<$HSP{n}; $b++)\n	      {\n		my (\
$string, $qs,$qe,$qf,$hs,$he,$hf,$s, $d, $e);\n		$\
string=$HSP{$b}{body};\n	\n		$qs=xmltag2value  ($s\
tring,\"Hsp_query-from\");\n		$qe=xmltag2value  ($\
string,\"Hsp_query-to\");\n		$qf=xmltag2value  ($s\
tring,\"Hsp_query-frame\");\n\n		$hs=xmltag2value \
 ($string,\"Hsp_hit-from\");\n		$he=xmltag2value  \
($string,\"Hsp_hit-to\");\n		$hf=xmltag2value  ($s\
tring,\"Hsp_hit-frame\");\n		\n		$s=xmltag2value  \
($string,\"Hsp_identity\");\n		$l=xmltag2value  ($\
string,\"Hsp_align-len\");\n		$s=int(($s*100)/$l);\
\n		\n		if ($qf>0)\n		  {$rqs=$qs; $rqe=$qe;}\n		e\
lse\n		  {\n		    $rqe=($qlen-$qs)+1;\n		    $rqs=\
($qlen-$qe)+1;\n		  }\n		\n		if ($hf>0)\n		  {$rhs\
=$hs; $rhe=$he;}\n		else\n		  {\n		    $rhe=($hlen\
-$hs)+1;\n		    $rhs=($hlen-$he)+1;\n		  }\n		for \
($d=0,$e=$rqs; $e<$rqe; $e++,$d++)\n		  {\n		    m\
y ($r1,$r2);\n		    $r1=$e;\n		    $r2=$rhs+$d;\n	\
	    print F \" $r1 $r2 $s 0\\n\";\n		  }\n	      \
}\n	  }\n      }\n    print F \"! SEQ_1_TO_N\\n\";\
\n    \n    close (F);\n    return %lib;\n  }\n\ns\
ub ncbi_tblastpx_xml2lib_file \n  {\n    my  ($out\
lib,$string,%s)=(@_);\n    my ($L,$l, $a,$b,$c,$d,\
$i,$nhits,@identifyerL);\n    my (%ITERATION,%hdes\
, %qdes);\n      \n    open (F, \">>$outlib\");\n \
   \n    $seq=~s/[^a-zA-Z]//g;\n    $L=length ($se\
q);\n    \n    %ITERATION=xml2tag_list ($string, \\
"Iteration\");\n    for ($i=0; $i<$ITERATION{n};$i\
++)\n      {\n	my ($qindex, $qlen, %hit, $string);\
\n	$string=$ITERATION{$i}{body};\n\n	$qdef=xmltag2\
value($string,\"Iteration_query-def\");\n	%qdes=&t\
blastpx_name2description($qdef,%s);\n	$qlen  =xmlt\
ag2value($string,\"Iteration_query-len\");\n	%hit=\
&xml2tag_list  ($string, \"Hit\");\n\n	for ($a=0; \
$a<$hit{n}; $a++)\n	  {\n	    my ($string);\n	    \
$string=$hit{$a}{body};\n	    $hdef=xmltag2value($\
string,\"Hit_def\");\n	    %hdes=&tblastpx_name2de\
scription($hdef,%s);\n	    if ($hdes{index}<=$qdes\
{index}){next;}\n	    else  {print F  \"# $qdes{in\
dex} $hdes{index}\\n\";}\n		   \n	   \n	    $hlen=\
xmltag2value  ($string,\"Hit_len\");\n	    %HSP=&x\
ml2tag_list  ($string, \"Hsp\");\n	   \n	    for (\
$b=0; $b<$HSP{n}; $b++)\n	      {\n		my ($string, \
$l,$qs,$qe,$qf,$hs,$he,$hf,$s, $d, $e, @s1, @s2);\\
n		$string=$HSP{$b}{body};\n	\n		$qs=xmltag2value \
 ($string,\"Hsp_query-from\");\n		$qe=xmltag2value\
  ($string,\"Hsp_query-to\");\n		$qf=$qdes{frame};\
\n		$qseq=xmltag2value  ($string,\"Hsp_qseq\");\n	\
	\n		$hs=xmltag2value  ($string,\"Hsp_hit-from\");\
\n		$he=xmltag2value  ($string,\"Hsp_hit-to\");\n	\
	$hf=$hdes{frame};\n		$hseq=xmltag2value  ($string\
,\"Hsp_hseq\");\n		\n		$s=xmltag2value  ($string,\\
"Hsp_identity\");\n		$l=xmltag2value  ($string,\"H\
sp_align-len\");\n		$s=int(($s*100)/$l);\n		@s1=tb\
lastpx_hsp2coordinates($qseq,$qs,$qe,%qdes);\n		@s\
2=tblastpx_hsp2coordinates($hseq,$hs,$he,%hdes);\n\
		\n		\n		for ($f=0; $f<=$#s1; $f++)\n		  {\n		   \
 if ($s1[$f]==-1 || $s2[$f]==-1){next;}\n		    els\
e \n		      {\n			print F \" $s1[$f] $s2[$f] $s 0\\
\n\";\n		      }\n		  }\n	      }\n	  }\n      }\n\
    print F \"! SEQ_1_TO_N\\n\";\n    \n    close \
(F);\n    return %lib;\n  }\nsub tblastpx_hsp2coor\
dinates\n  {\n    my ($seq, $s, $e, %des)=@_;\n   \
 my @list;\n    my @sa;\n    my @gap=(-1,-1,-1);\n\
    \n    $s=$des{start}+3*($s-1);\n  \n    if ($d\
es{strand} eq \"d\"){;}\n    else {$s=($des{length\
}-$s)+1;}\n    \n    foreach $c (split (//,$seq))\\
n      {\n	if ( $c eq '-'){push (@list,@gap);}\n	e\
lsif ($des{strand} eq \"d\")\n	  {\n	    push(@lis\
t,$s++,$s++,$s++);\n	  }\n	else\n	  {\n	    push(@\
list, $s--,$s--,$s--);\n	  }\n      }\n    return \
@list;\n  }\n\nsub tblastpx_name2description\n  {\\
n    my ($name, %s)=@_;\n    my @at=split(\"__\", \
$name);\n    my %des;\n\n    $des{name}=$at[0];\n \
   $des{strand}=$at[1];\n    \n    $des{start}=$at\
[2];\n    $des{end}=$at[3];\n    $des{length}=$at[\
4];\n    $des{index}=$s{$at[0]}{order}+1;\n    ret\
urn %des;\n  }  \nsub ncbi_blast_xml2profile \n  {\
\n    my ($name,$seq,$maxid, $minid, $mincov, $str\
ing)=(@_);\n    my ($L,$l, $a,$b,$c,$d,$nhits,@ide\
ntifyerL);\n    \n    \n    $seq=~s/[^a-zA-Z]//g;\\
n    $L=length ($seq);\n    \n    %hit=&xml2tag_li\
st ($string, \"Hit\");\n    \n    \n    for ($nhit\
s=0,$a=0; $a<$hit{n}; $a++)\n      {\n	my ($ldb,$i\
d, $identity, $expectation, $start, $end, $coverag\
e, $r);\n	my (%ID,%DE,%HSP);\n	\n	$ldb=\"\";\n\n	%\
ID=&xml2tag_list ($hit{$a}{body}, \"Hit_id\");\n	$\
identifyer=$ID{0}{body};\n	\n	%DE=&xml2tag_list ($\
hit{$a}{body}, \"Hit_def\");\n	$definition=$DE{0}{\
body};\n	\n	%HSP=&xml2tag_list ($hit{$a}{body}, \"\
Hsp\");\n	for ($b=0; $b<$HSP{n}; $b++)\n	  {\n	   \
 my (%START,%END,%E,%I,%Q,%M);\n\n	 \n	    %START=\
&xml2tag_list ($HSP{$b}{body}, \"Hsp_query-from\")\
;\n	    %HSTART=&xml2tag_list ($HSP{$b}{body}, \"H\
sp_hit-from\");\n	    \n	    %LEN=  &xml2tag_list \
($HSP{$b}{body}, \"Hsp_align-len\");\n	    %END=  \
&xml2tag_list ($HSP{$b}{body}, \"Hsp_query-to\");\\
n	    %HEND=  &xml2tag_list ($HSP{$b}{body}, \"Hsp\
_hit-to\");\n	    %E=&xml2tag_list     ($HSP{$b}{b\
ody}, \"Hsp_evalue\");\n	    %I=&xml2tag_list     \
($HSP{$b}{body}, \"Hsp_identity\");\n	    %Q=&xml2\
tag_list     ($HSP{$b}{body}, \"Hsp_qseq\");\n	   \
 %M=&xml2tag_list     ($HSP{$b}{body}, \"Hsp_hseq\\
");\n	    \n	    for ($e=0; $e<$Q{n}; $e++)\n\n	  \
    {\n		$qs=$Q{$e}{body};\n		$ms=$M{$e}{body};\n	\
	\n		$expectation=$E{$e}{body};\n		$identity=($LEN\
{$e}{body}==0)?0:$I{$e}{body}/$LEN{$e}{body}*100;\\
n		$start=$START{$e}{body};\n		$end=$END{$e}{body}\
;\n		$Hstart=$HSTART{$e}{body};\n		$Hend=$HEND{$e}\
{body};\n	\n		$coverage=(($end-$start)*100)/$L;\n	\
\n		if ($identity>$maxid || $identity<$minid || $c\
overage<$mincov){next;}\n		@lr1=(split (//,$qs));\\
n		@lr2=(split (//,$ms));\n		$l=$#lr1+1;\n		for ($\
c=0;$c<$L;$c++){$p[$nhits][$c]=\"-\";}\n		for ($d=\
0,$c=0; $c<$l; $c++)\n		  {\n		    $r=$lr1[$c];\n	\
	    if ( $r=~/[A-Za-z]/)\n		      {\n			\n			$p[$\
nhits][$d + $start-1]=$lr2[$c];\n			$d++;\n		     \
 }\n		  }\n		$Qseq[$nhits]=$qs;\n		$Hseq[$nhits]=$\
ms;\n		$QstartL[$nhits]=$start;\n		$HstartL[$nhits\
]=$Hstart;\n		$identityL[$nhits]=$identity;\n		$en\
dL[$nhits]=$end;\n		$definitionL[$nhits]=$definiti\
on;\n		$identifyerL[$nhits]=$identifyer;\n		$comme\
nt[$nhits]=\"$ldb|$identifyer [Eval=$expectation][\
id=$identity%][start=$Hstart end=$Hend]\";\n		$nhi\
ts++;\n	      }\n	  }\n      }\n    \n    $profile\
{n}=0;\n    $profile{$profile{n}}{name}=$name;\n  \
  $profile{$profile{n}}{seq}=$seq;\n    $profile {\
n}++;\n    \n    for ($a=0; $a<$nhits; $a++)\n    \
  {\n	$n=$a+1;\n	\n	$profile{$n}{name}=\"$name\\_$\
a\";\n	$profile{$n}{seq}=\"\";\n	$profile{$n}{Qseq\
}=$Qseq[$a];\n	$profile{$n}{Hseq}=$Hseq[$a];\n	$pr\
ofile{$n}{Qstart}=$QstartL[$a];\n	$profile{$n}{Hst\
art}=$HstartL[$a];\n	$profile{$n}{identity}=$ident\
ityL[$a];\n	$profile{$n}{definition}=$definitionL[\
$a];\n	$profile{$n}{identifyer}=$identifyerL[$a];\\
n	$profile{$n}{comment}=$comment[$a];\n\n	for ($b=\
0; $b<$L; $b++)\n	  {\n	    if ($p[$a][$b])\n	    \
  {\n		$profile{$n}{seq}.=$p[$a][$b];\n	      }\n	\
    else\n	      {\n		$profile{$n}{seq}.=\"-\";\n	\
      }\n	  }\n      }\n    \n    $profile{n}=$nhi\
ts+1;\n    return %profile;\n  }\nsub ebi_blast_xm\
l2profile \n  {\n    my ($name,$seq,$maxid, $minid\
, $mincov, $string)=(@_);\n    my ($L,$l, $a,$b,$c\
,$d,$nhits,@identifyerL,$identifyer);\n    \n\n   \
 \n    $seq=~s/[^a-zA-Z]//g;\n    $L=length ($seq)\
;\n    %hit=&xml2tag_list ($string, \"hit\");\n   \
 \n    for ($nhits=0,$a=0; $a<$hit{n}; $a++)\n    \
  {\n	my ($ldb,$id, $identity, $expectation, $star\
t, $end, $coverage, $r);\n	my (%Q,%M,%E,%I);\n	\n	\
$ldb=&tag2value ($hit{$a}{open}, \"database\");\n	\
$identifyer=&tag2value ($hit{$a}{open}, \"id\");\n\
\n	$description=&tag2value ($hit{$a}{open}, \"desc\
ription\");\n	\n	%Q=&xml2tag_list ($hit{$a}{body},\
 \"querySeq\");\n	%M=&xml2tag_list ($hit{$a}{body}\
, \"matchSeq\");\n	%E=&xml2tag_list ($hit{$a}{body\
}, \"expectation\");\n	%I=&xml2tag_list ($hit{$a}{\
body}, \"identity\");\n	\n\n	for ($b=0; $b<$Q{n}; \
$b++)\n	  {\n\n	    $qs=$Q{$b}{body};\n	    $ms=$M\
{$b}{body};\n	    \n	    $expectation=$E{$b}{body}\
;\n	    $identity=$I{$b}{body};\n	    \n	    	    \
\n	    $start=&tag2value ($Q{$b}{open}, \"start\")\
;\n	    $end=&tag2value ($Q{$b}{open}, \"end\");\n\
	    $startM=&tag2value ($M{$b}{open}, \"start\");\
\n	    $endM=&tag2value ($M{$b}{open}, \"end\");\n\
	    $coverage=(($end-$start)*100)/$L;\n	    \n	  \
 # print \"$id: ID: $identity COV: $coverage [$sta\
rt $end]\\n\";\n	    \n	    \n	    if ($identity>$\
maxid || $identity<$minid || $coverage<$mincov){ne\
xt;}\n	    # print \"KEEP\\n\";\n\n	    \n	    @lr\
1=(split (//,$qs));\n	    @lr2=(split (//,$ms));\n\
	    $l=$#lr1+1;\n	    for ($c=0;$c<$L;$c++){$p[$n\
hits][$c]=\"-\";}\n	    for ($d=0,$c=0; $c<$l; $c+\
+)\n	      {\n		$r=$lr1[$c];\n		if ( $r=~/[A-Za-z]\
/)\n		  {\n		    \n		    $p[$nhits][$d + $start-1]\
=$lr2[$c];\n		    $d++;\n		  }\n	      }\n	  \n	  \
  $Qseq[$nhits]=$qs;\n	    $Hseq[$nhits]=$ms;\n	  \
  $QstartL[$nhits]=$start;\n	    $HstartL[$nhits]=\
$Hstart;\n	    $identityL[$nhits]=$identity;\n	   \
 $endL[$nhits]=$end;\n	    $definitionL[$nhits]=$d\
efinition;\n	    $identifyerL[$nhits]=$identifyer;\
\n	    $comment[$nhits]=\"$ldb|$identifyer [Eval=$\
expectation][id=$identity%][start=$startM end=$end\
M]\";\n	    $nhits++;\n	  }\n      }\n    \n    $p\
rofile{n}=0;\n    $profile{$profile{n}}{name}=$nam\
e;\n    $profile{$profile{n}}{seq}=$seq;\n    $pro\
file {n}++;\n    \n    for ($a=0; $a<$nhits; $a++)\
\n      {\n	$n=$a+1;\n	$profile{$n}{name}=\"$name\\
\_$a\";\n	$profile{$n}{seq}=\"\";\n	$profile{$n}{Q\
seq}=$Qseq[$a];\n	$profile{$n}{Hseq}=$Hseq[$a];\n	\
$profile{$n}{Qstart}=$QstartL[$a];\n	$profile{$n}{\
Hstart}=$HstartL[$a];\n	$profile{$n}{identity}=$id\
entityL[$a];\n	$profile{$n}{definition}=$definitio\
nL[$a];	\n	$profile{$n}{identifyer}=$identifyerL[$\
a];\n	$profile{$n}{comment}=$comment[$a];\n\n	for \
($b=0; $b<$L; $b++)\n	  {\n	    if ($p[$a][$b])\n	\
      {\n		$profile{$n}{seq}.=$p[$a][$b];\n	      \
}\n	    else\n	      {\n		$profile{$n}{seq}.=\"-\"\
;\n	      }\n	  }\n      }\n    $profile{n}=$nhits\
+1;\n    \n    return %profile;\n  }\nsub output_p\
rofile\n  {\n    my ($outfile,$profileR, $trim)=(@\
_);\n    my ($a);\n    my %profile=%$profileR;\n  \
  my $P= new FileHandle;\n    my $tmp=vtmpnam();\n\
    \n    open ($P, \">$tmp\");\n    for ($a=0; $a\
<$profile{n}; $a++)\n      {\n	print $P \">$profil\
e{$a}{name} $profile{$a}{comment}\\n$profile{$a}{s\
eq}\\n\";\n      }\n    close ($P);\n\n    if ( $t\
rim)\n      {\n	&safe_system (\"t_coffee -other_pg\
 seq_reformat -in $tmp -action +trim _aln_%%$trim\\
\_K1 -output fasta_aln -out $outfile\");\n      }\\
n    else\n      {\n	&safe_system (\"mv $tmp $outf\
ile\");\n      }\n    return;\n  }\nsub blast_xml2\
hit_list\n  {\n    my $string=(@_[0]);\n    return\
 &xml2tag_list ($string, \"hit\");\n  }\nsub xmlta\
g2value\n  {\n    my ($string_in, $tag)=@_;\n    m\
y %TAG;\n    %TAG=xml2tag_list ($string_in, $tag);\
\n    return $TAG{0}{body};\n  }\n      \nsub xml2\
tag_list  \n  {\n    my ($string_in,$tag)=@_;\n   \
 my $tag_in, $tag_out;\n    my %tag;\n    \n    if\
 (-e $string_in)\n      {\n	$string=&file2string (\
$string_in);\n      }\n    else\n      {\n	$string\
=$string_in;\n      }\n    $tag_in1=\"<$tag \";\n \
   $tag_in2=\"<$tag>\";\n    $tag_out=\"/$tag>\";\\
n    $string=~s/>/>##1/g;\n    $string=~s/</##2</g\
;\n    $string=~s/##1/<#/g;\n    $string=~s/##2/#>\
/g;\n    @l=($string=~/(\\<[^>]+\\>)/g);\n    $tag\
{n}=0;\n    $in=0;$n=-1;\n  \n \n\n    foreach $t \
(@l)\n      {\n\n	$t=~s/<#//;\n	$t=~s/#>//;\n	\n	i\
f ( $t=~/$tag_in1/ || $t=~/$tag_in2/)\n	  {\n	 \n	\
    $in=1;\n	    $tag{$tag{n}}{open}=$t;\n	    $n+\
+;\n	    \n	  }\n	elsif ($t=~/$tag_out/)\n	  {\n	 \
   \n\n	    $tag{$tag{n}}{close}=$t;\n	    $tag{n}\
++;\n	    $in=0;\n	  }\n	elsif ($in)\n	  {\n	   \n\
	    $tag{$tag{n}}{body}.=$t;\n	  }\n      }\n  \n\
    return %tag;\n  }\n\n\nsub seq2gor_prediction \
\n  {\n    my ($name, $seq,$infile, $outfile, $gor\
_seq, $gor_obs)=(@_);\n    my ($l);\n    \n    `go\
rIV -prd $infile -seq $gor_seq -obs $gor_obs > gor\
_tmp`;\n    open (GR, \">$outfile\");\n    open (O\
G, \"gor_tmp\");\n\n    while (<OG>)\n      {\n	\n\
	$l=$_;\n	if ($l=~/\\>/){print GR \"$l\";}\n	elsif\
 ( $l=~/Predicted Sec. Struct./)\n	  {\n	    $l=~s\
/Predicted Sec. Struct\\.//;\n	    print GR \"$l\"\
;\n	  }\n      }\n    close (GR);\n    close (OG);\
\n    return;\n  }\nsub seq2msa_tm_prediction \n  \
{\n    my ($name, $seq, $db, $infile, $outfile, $a\
rch, $psv)=(@_);\n    my (%p,%gseq,%R, $blast_outp\
ut, %s, $l);\n    my $R2=new FileHandle;\n    my $\
db=\"uniprot\";\n    my $method=\"psitm\";\n    my\
 $SERVER=\"EBI\";\n    \n    $blast_output=&run_bl\
ast ($name,\"blastp\", $db, $infile, \"outfile\");\
\n    \n    if (&cache_file(\"GET\",$infile,$name,\
$method,$db,$outfile,$SERVER))\n      {\n	print \"\
\\tPSITM: USE Cache\\n\";\n	return $outfile;\n    \
  }\n    else\n      {\n	$CACHE_STATUS=\"COMPUTE C\
ACHE\";\n	%p=blast_xml2profile($name,$seq,$maxid, \
$minid,$mincov,$blast_output);\n	\n	\n	open (F, \"\
>tm_input\");\n	for (my $a=0; $a<$p{n}; $a++)\n	  \
{\n	    my $s;\n	    \n	    $s=$p{$a}{seq};\n	    \
$s=uc($s);\n	    print F \">$p{$a}{name}\\n$s\\n\"\
;\n	    #print stdout \">$p{$a}{name}\\n$s\\n\";\n\
	  }\n	close (F);\n	print \"\\tPSITM: kept  $p{n} \
Homologues for Sequence $p{0}{name}\\n\";\n	&safe_\
system (\"t_coffee -other_pg fasta_seq2hmmtop_fast\
a.pl -in=tm_input -out=$outfile -output=cons -cov=\
70 -trim=95 -arch=$arch -psv=$psv\");\n	unlink (\"\
tm_input\");\n	&cache_file(\"SET\",$infile,$name,$\
method,$db,$outfile,$SERVER);\n	return;\n      }\n\
  }\n\n\nsub seq2msa_gor_prediction \n  {\n    my \
($name, $seq,$infile, $outfile, $gor_seq, $gor_obs\
)=(@_);\n    my (%p,%gseq,%R, $blast_output, %s, $\
l);\n    my $R2=new FileHandle;\n    my $db=\"unip\
rot\";\n    my $method=\"psigor\";\n    my $SERVER\
=\"EBI\";\n    \n    $blast_output=&run_blast ($na\
me,\"blastp\", \"uniprot\", $infile, \"outfile\");\
\n    \n    if (&cache_file(\"GET\",$infile,$name,\
$method,$db,$outfile,$SERVER))\n      {\n	print \"\
\\tPSIGOR: USE Cache\\n\";\n	return $outfile;\n   \
   }\n    else\n      {\n	$CACHE_STATUS=\"COMPUTE \
CACHE\";\n	%p=blast_xml2profile($name,$seq,$maxid,\
 $minid,$mincov,$blast_output);\n	\n	\n	open (F, \\
">gor_input\");\n	for (my $a=0; $a<$p{n}; $a++)\n	\
  {\n	    my $s;\n	    \n	    $s=$p{$a}{seq};\n	  \
  $s=uc($s);\n	    print F \">$p{$a}{name}\\n$s\\n\
\";\n	    #print stdout \">$p{$a}{name}\\n$s\\n\";\
\n	  }\n	close (F);\n	print \"\\tGORTM: kept  $p{n\
} Homologues for Sequence $p{0}{name}\\n\";\n	&saf\
e_system (\"t_coffee -other_pg fasta_seq2hmmtop_fa\
sta.pl -in=gor_input -out=$outfile -output=cons -c\
ov=70 -trim=95 -gor_seq=$gor_seq -gor_obs=$gor_obs\
 -mode=gor\");\n	unlink (\"tm_input\");\n	&cache_f\
ile(\"SET\",$infile,$name,$method,$db,$outfile,$SE\
RVER);\n	return;\n      }\n  }\n\n\n\nsub run_blas\
t\n  {\n    my ($name, $method, $db, $infile, $out\
file, $run)=(@_);\n    if (!$run){$run=1;}\n    \n\
    \n    if (&cache_file(\"GET\",$infile,$name,$m\
ethod,$db,$outfile,$SERVER) && is_valid_blast_xml \
($outfile))\n      {return $outfile;}\n    else\n \
     {\n	$CACHE_STATUS=\"COMPUTE CACHE\";\n	if ( $\
SERVER eq \"EBI_SOAP\")\n	  {\n	    &check_configu\
ration (\"EMAIL\",\"SOAP::Light\",\"INTERNET\");\n\
	    \n	    $cl_method=$method;\n	    if ($cl_meth\
od =~/wu/)\n	      {\n		$cl_method=~s/wu//;\n		if \
( $cl_method eq \"psiblast\")\n		  {\n		    add_wa\
rning($$,$$,\"PSI BLAST cannot be used with the wu\
BLAST Client. Use server=EBI Or server=LOCAL. blas\
tp will be used instead\");\n		    $cl_method=\"bl\
astp\";\n		  }\n		\n		$command=\"t_coffee -other_p\
g wublast.pl --email $EMAIL $infile -D $db -p $cl_\
method --outfile $outfile -o xml>/dev/null 2>/dev/\
null\";\n		&safe_system ( $command);\n		if (-e \"$\
outfile.xml\") {`mv $outfile.xml $outfile`;}\n	   \
   }\n	    else\n	      {\n		if ($cl_method eq \"p\
siblast\"){$cl_method =\"blastp -j5\";}\n		\n		$co\
mmand=\"t_coffee -other_pg blastpgp.pl --email $EM\
AIL $infile -d $db --outfile $outfile -p $cl_metho\
d --mode PSI-Blast>/dev/null 2>/dev/null\";\n		&sa\
fe_system ( $command);\n		\n		if (-e \"$outfile.xm\
l\") {`mv $outfile.xml $outfile`;}\n	      }\n	  }\
\n	elsif ($SERVER eq \"EBI_REST\" || $SERVER eq \"\
EBI\")\n	  {\n	   \n	    $cl_method=$method;\n	   \
 &check_configuration(\"EMAIL\",\"XML::Simple\", \\
"INTERNET\");\n	    if ($db eq \"uniprot\"){$db1=\\
"uniprotkb\";}\n	    else {$db1=$db;}\n	    \n\n	 \
   if ($cl_method =~/wu/)\n	      {\n		$cl_method=\
~s/wu//;\n		if ( $cl_method eq \"psiblast\")\n		  \
{\n		    $cl_method=\"blastp\";\n		  }\n		\n		$com\
mand=\"t_coffee -other_pg wublast_lwp.pl --email $\
EMAIL -D $db1 -p $cl_method --outfile $outfile --o\
utformat xml --stype protein $infile>/dev/null 2>/\
dev/null\";\n		&safe_system ( $command,5);\n		if (\
-e \"$outfile.xml\") {`mv $outfile.xml $outfile`;}\
\n		elsif (-e \"$outfile.xml.xml\"){`mv $outfile.x\
ml.xml $outfile`;}\n	      }\n	    else\n	      {\\
n		if ( $cl_method =~/psiblast/){$cl_method =\"bla\
stp -j5\";}\n		$command=\"t_coffee -other_pg ncbib\
last_lwp.pl --email $EMAIL -D $db1 -p $cl_method -\
-outfile $outfile --outformat xml --stype protein \
$infile>/dev/null 2>/dev/null\";\n		#$command=\"t_\
coffee -other_pg ncbiblast_lwp.pl --email $EMAIL -\
D $db1 -p $cl_method --outfile $outfile --outforma\
t xml --stype protein $infile>/dev/null\";\n		&saf\
e_system ( $command,5);\n		if (-e \"$outfile.xml\"\
) {`mv $outfile.xml $outfile`;}\n		elsif (-e \"$ou\
tfile.xml.xml\"){`mv $outfile.xml.xml $outfile`;} \
\n		\n	      }\n	    \n	 }\n	elsif ($SERVER eq \"N\
CBI\")\n	  {\n	    &check_configuration (\"blastcl\
3\",\"INTERNET\");\n	    if ($db eq \"uniprot\"){$\
cl_db=\"nr\";}\n	    else {$cl_db=$db;}\n	    \n	 \
   if ( $method eq \"psiblast\")\n	      {\n		add_\
warning($$,$$,\"PSI BLAST cannot be used with the \
NCBI BLAST Client. Use server=EBI Or server=LOCAL.\
 blastp will be used instead\");\n		$cl_method=\"b\
lastp\";\n	      }\n	    else\n	      {\n		$cl_met\
hod=$method;\n	      }\n	    $command=\"blastcl3 -\
p $cl_method -d $cl_db -i $infile -o $outfile -m 7\
\";\n	    &safe_system ($command);\n	  }\n	elsif (\
$SERVER =~/CLIENT_(.*)/)\n	  {\n	    my $client=$1\
;\n	    $command=\"$client -p $method -d $db -i $i\
nfile -o $outfile -m 7\";\n	    &safe_system ($com\
mand);\n	  }\n	elsif ( $SERVER eq \"LOCAL_blastall\
\")\n	  {\n	    &check_configuration (\"blastall\"\
);\n	    if ($method eq \"blastp\")\n	      {\n		$\
command=\"blastall -d $db -i $infile -o $outfile -\
m7 -p blastp\";\n	      }\n	    &safe_system ($com\
mand);\n	  }\n	elsif ( $SERVER eq \"LOCAL\")\n	  {\
\n	    \n	    if ($ENV{\"BLAST_DB_DIR\"})\n	      \
{\n		$x=$ENV{\"BLAST_DB_DIR\"};\n		$cl_db=\"$x$db\\
";\n	      }\n	    else\n	      {\n		$cl_db=$db;\n\
	      }\n	    \n	    if ($method eq \"blastp\")\n\
	      {\n		&check_configuration(\"blastpgp\");\n	\
	$command=\"blastpgp -d $cl_db -i $infile -o $outf\
ile -m7 -j1\";\n	      }\n	    elsif ($method eq \\
"psiblast\")\n	      {\n		&check_configuration(\"b\
lastpgp\");\n		$command=\"blastpgp -d $cl_db -i $i\
nfile -o $outfile -m7 -j5\";\n	      }\n	    elsif\
 ($method eq \"blastn\")\n	      {\n		&check_confi\
guration(\"blastall\");\n		$command=\"blastall -p \
blastn -d $cl_db -i $infile -o $outfile -m7 -W6\";\
\n	      }	\n	    &safe_system ($command);\n	  }\n\
	else\n	  {\n	    \n	    myexit(add_error (EXIT_FA\
ILURE,$$,$$,getppid(), \"BLAST_FAILURE::UnknownSer\
ver\",$CL));\n	  }\n	\n	if ( !-e $outfile || !&is_\
valid_blast_xml($outfile))\n	  {\n	    \n	    if (\
 -e $outfile)\n	      {\n		add_warning ($$,$$,\"Co\
rrupted Blast Output (Run $run)\");\n		unlink($out\
file);\n	      }\n	    \n	    if ( $run==$BLAST_MA\
X_NRUNS)\n	      {\n	\n		myexit(add_error (EXIT_FA\
ILURE,$$,$$,getppid(), \"BLAST_FAILURE::UnknownRea\
son\", \"$command\"));\n	      }\n	    else\n	    \
  {\n	\n		add_warning ($$,$$,\"Blast for $name fai\
led (Run: $run)\");\n		\n		return run_blast ($name\
, $method, $db,$infile, $outfile, $run+1);\n	     \
 }\n	  }\n	\n	&cache_file(\"SET\",$infile,$name,$m\
ethod,$db,$outfile,$SERVER);\n	return $outfile;\n \
     }\n  }\n\nsub cache_file\n  {\n    my ($cache\
_mode,$infile,$name,$method,$db, $outfile,$server)\
=(@_);\n    my $cache_file;\n    #Protect names so\
 that they can be turned into legal filenames\n   \
 $name=&clean_file_name ($name);\n\n    if ($db=~/\
\\//)\n      {\n	$db=~/([^\\/]+)$/;\n	$db=$1;\n   \
   }\n    $cache_file_sh=\"$name.$method.$db.$serv\
er.tmp\";\n    $cache_file=\"$CACHE/$name.$method.\
$db.$server.tmp\";\n    \n    if ($infile ne \"\")\
\n      {\n	$cache_file_infile_sh=\"$name.$method.\
$db.$server.infile.tmp\";\n	$cache_file_infile=\"$\
CACHE/$name.$method.$db.$server.infile.tmp\";\n   \
   }\n    \n    if ($cache_mode eq \"GET\")\n     \
 {\n	if ($CACHE eq \"\" || $CACHE eq \"no\" || $CA\
CHE eq \"ignore\"  || $CACHE eq \"local\" || $CACH\
E eq \"update\"){return 0;}\n	elsif ( !-d $CACHE)\\
n	  {\n	    print STDERR \"ERROR: Cache Dir: $CACH\
E Does not Exist\";\n	    return 0;\n	  }\n	else\n\
	  {\n	    if ( -e $cache_file && &fasta_file1_eq_\
fasta_file2($infile,$cache_file_infile)==1)\n	    \
  {\n		`cp $cache_file $outfile`;\n		$CACHE_STATUS\
=\"READ CACHE\";\n		return 1;\n	      }\n	  }\n   \
   }\n    elsif ($cache_mode eq \"SET\")\n      {\\
n	if ($CACHE eq \"\" || $CACHE eq \"no\" || $CACHE\
 eq \"ignore\"  || $CACHE eq \"local\" || $CACHE e\
q \"update\"){return 0;}\n	elsif ( !-d $CACHE)\n	 \
 {\n	    print STDERR \"ERROR: Cache Dir: $CACHE D\
oes not Exist\";\n	    return 0;\n	  }\n	elsif (-e\
 $outfile)\n	  {\n	    `cp $outfile $cache_file`;\\
n	    if ($cache_file_infile ne \"\"){ `cp $infile\
 $cache_file_infile`;}\n\n	    #functions for upda\
ting the cache\n	    #`t_coffee -other_pg clean_ca\
che.pl -file $cache_file_sh -dir $CACHE`;\n	    #`\
t_coffee -other_pg clean_cache.pl -file $cache_fil\
e_infile_sh -dir $CACHE`;\n	    return 1;\n	  }\n \
     }\n    $CACHE_STATUS=\"COMPUTE CACHE\";\n    \
return 0;\n  }\nsub file1_eq_file2\n  {\n    my ($\
f1, $f2)=@_;\n    if ( $f1 eq \"\"){return 1;}\n  \
  elsif ( $f2 eq \"\"){return 1;}\n    elsif ( !-e\
 $f1){return 0;}\n    elsif ( !-e $f2){return 0;}\\
n    elsif ($f1 eq \"\" || $f2 eq \"\" || `diff $f\
1 $f2` eq \"\"){return 1;}\n    \n    return 0;\n \
 }\nsub clean_file_name \n  {\n    my $name=@_[0];\
\n    \n    $name=~s/[^A-Za-z1-9.-]/_/g;\n    retu\
rn $name;\n  }\nsub url2file\n  {\n    my ($addres\
s, $out)=(@_);\n    \n    if (&pg_is_installed (\"\
wget\"))\n	{\n	  return &safe_system (\"wget $addr\
ess -O$out >/dev/null 2>/dev/null\");\n	}\n    els\
if (&pg_is_installed (\"curl\"))\n      {\n	return\
 &safe_system (\"curl $address -o$out >/dev/null 2\
>/dev/null\");\n      }\n    else\n      {\n	myexi\
t(flus_error(\"neither curl nor wget are installed\
. Imnpossible to fectch remote file\"));\n	exit ($\
EXIT_FAILURE);\n      }\n  }\nsub fasta_file1_eq_f\
asta_file2\n  {\n    my ($f1, $f2)=@_;\n    my (%s\
1, %s2);\n    my @names;\n    %s1=read_fasta_seq (\
$f1);\n    %s2=read_fasta_seq ($f2);\n\n    @names\
=(keys (%s1));\n    \n    foreach $n (keys(%s1))\n\
      {\n	if ($s1{$n}{seq} ne $s2{$n}{seq}){return\
 0;}\n      } \n    \n    foreach $n (keys(%s2))\n\
      {\n	if ($s1{$n}{seq} ne $s2{$n}{seq}){return\
 0;}\n      }\n    return 1;\n  }\n	\n\n\nsub read\
_template_file\n{\n	my $pdb_templates = @_[0];\n	o\
pen (TEMP, \"<$pdb_templates\");\n	my %temp_h;\n	w\
hile (<TEMP>)\n{\n		$line = $_;\n 		$line =~/(\\S+\
)\\s(\\S+)/;\n 		$temp_h{$1}= $2;\n}\n	close(TEMP)\
;\n	return %temp_h;\n}\n\nsub calc_rna_template\n{\
\n	my ($mode, $infile, $pdbfile, $outfile)=@_;\n	m\
y %s, %h ;\n	my $result;\n	my (@profiles);\n	&set_\
temporary_dir (\"set\",$infile,\"seq.pep\");\n	%s=\
read_fasta_seq (\"seq.pep\");\n	\n	%pdb_template_h\
 = &read_template_file($pdbfile);\n	my $pdb_chain;\
\n	open (R, \">result.aln\");\n\n\n	#print stdout \
\"\\n\";\n	foreach $seq (keys(%s))\n	{\n		if ($pdb\
_template_h{$seq} eq \"\")\n		{\n			next;\n		}\n		\
open (F, \">seqfile\");\n		print (F \">$s{$seq}{na\
me}\\n$s{$seq}{seq}\\n\");\n		close (F);\n		$pdb_c\
hain = $pdb_template_h{$seq};\n		$lib_name=\"$s{$s\
eq}{name}.rfold\";\n		$lib_name=&clean_file_name (\
$lib_name);\n		\n 		safe_system (\"secondary_struc\
.py seqfile $CACHE$pdb_chain  $lib_name\");\n		\n	\
	if ( !-e $lib_name)\n		{\n		myexit(flush_error(\"\
RNAplfold failed to compute the secondary structur\
e of $s{$seq}{name}\"));\n			myexit ($EXIT_FAILURE\
);\n		}\n		else\n		{\n			print stdout \"\\tProcess\
: >$s{$seq}{name} _F_ $lib_name\\n\";\n			print R \
\">$s{$seq}{name} _F_ $lib_name\\n\";\n		}\n		unsh\
ift (@profiles, $lib_name);\n	}\n	close (R);\n	&se\
t_temporary_dir (\"unset\",$mode, $method,\"result\
.aln\",$outfile, @profiles);\n}\n\n\n\nsub seq2rna\
_pair{\n	my ($mode, $pdbfile1, $pdbfile2, $method,\
 $param, $outfile)=@_;\n	\n	if ($method eq \"runsa\
ra.py\")\n	{\n		open(TMP,\"<$pdbfile1\");\n		my $c\
ount = 0;\n		my $line;\n		while (<TMP>)\n		{\n			$\
line = $_;\n			if ($count ==1)\n			{\n				last;\n	\
		}\n			$count += 1;\n		}\n\n		\n		$chain1 = subst\
r($line,length($line)-3,1);\n\n		close TMP;\n		ope\
n(TMP,\"<$pdbfile2\");\n		my $count = 0;\n		while \
(<TMP>)\n		{\n			$line = $_;\n			if ($count ==1)\n\
			{\n				last;\n			}\n			$count += 1;\n		}\n		$ch\
ain2 = substr($line,length($line)-3,1);\n		close T\
MP;\n\n		$tmp_file=&vtmpnam();\n	\n		safe_system(\\
"runsara.py $pdbfile1 $chain1 $pdbfile2 $chain2 -s\
 -o $tmp_file --limitation 5000 > /dev/null 2> /de\
v/null\") == 0 or die \"sara did not work $!\\n\";\
\n		open(TMP,\"<$tmp_file\") or die \"cannot open \
the sara tmp file:$!\\n\";\n		open(OUT,\">$outfile\
\") or die \"cannot open the $outfile file:$!\\n\"\
;\n\n		my $switch = 0;\n		my $seqNum = 0;\n		forea\
ch my $line (<TMP>)\n		{\n			next unless ($line=~/\
SARAALI/);\n			if ($line=~/>/)\n			{\n				$switch \
=0;\n				print OUT \">seq$seqNum\\n\";\n				$seqNu\
m++;				\n			}\n			if ($switch < 2){\n				$switch+\
+;\n				next;\n			}\n	\n			if ($line =~/REMARK\\s+\
SARAALI\\s+([^\\*]+)\\*/)\n			{\n				my $string = \
$1;\n				print OUT \"$string\\n\";\n			}\n		}\n		c\
lose TMP; \n		close OUT;\n		unlink($tmp_file);\n	}\
\n}\n\nsub seq2tblastx_lib\n  {\n    my ($mode, $i\
nfile, $outfile)=@_;\n    my (%s, $method,$nseq);\\
n\n    $method=$mode;\n    &set_temporary_dir (\"s\
et\",$infile,\"infile\");\n    %s=read_fasta_seq(\\
"infile\");\n    \n    \n    foreach $seq (keys(%s\
))\n      {\n	$slist[$s{$seq}{order}]=$s{$seq}{seq\
};\n	$sname[$s{$seq}{order}]=$s{$seq}{name};\n	$sl\
en[$s{$seq}{order}]=length ($s{$seq}{seq});\n     \
 }\n    $nseq=$#sname+1;\n    open (F, \">outfile\\
");\n    print F \"! TC_LIB_FORMAT_01\\n\";\n    p\
rint F \"$nseq\\n\";\n    for ($a=0; $a<$nseq;$a++\
)\n      {\n	print F \"$sname[$a] $slen[$a]  $slis\
t[$a]\\n\"\n      }\n    close (F);\n    `formatdb\
 -i infile -p F`;\n\n    `blastall -p tblastx -i i\
nfile -d infile -m 7 -S1>blast.output`;\n    \n   \
 ncbi_tblastx_xml2lib_file (\"outfile\", file2stri\
ng (\"blast.output\"));\n    &set_temporary_dir (\\
"unset\",$mode, $method, \"outfile\",$outfile);\n \
   myexit ($EXIT_SUCCESS);\n    }\nsub seq2tblastp\
x_lib\n  {\n    my ($mode, $infile, $outfile)=@_;\\
n    my (%s, $method,$nseq);\n    $method=$mode;\n\
    &set_temporary_dir (\"set\",$infile,\"infile\"\
);\n    %s=read_fasta_seq(\"infile\");\n    \n    \
foreach $seq (keys(%s))\n      {\n	$slist[$s{$seq}\
{order}]=$s{$seq}{seq};\n	$sname[$s{$seq}{order}]=\
$s{$seq}{name};\n	$slen[$s{$seq}{order}]=length ($\
s{$seq}{seq});\n      }\n    $nseq=$#sname+1;\n   \
 open (F, \">outfile\");\n    print F \"! TC_LIB_F\
ORMAT_01\\n\";\n    print F \"$nseq\\n\";\n    for\
 ($a=0; $a<$nseq;$a++)\n      {\n	print F \"$sname\
[$a] $slen[$a]  $slist[$a]\\n\"\n      }\n    clos\
e (F);\n    `printenv > /home/notredame/tmp/ce`;\n\
    `t_coffee -other_pg seq_reformat -in infile -o\
utput tblastx_db1 > tblastxdb`;\n    `formatdb -i \
tblastxdb -p T`;\n    `blastall -p blastp -i tblas\
txdb -d tblastxdb -m7 >blast.output`;\n    ncbi_tb\
lastpx_xml2lib_file (\"outfile\", file2string (\"b\
last.output\"), %s);\n    &set_temporary_dir (\"un\
set\",$mode, $method, \"outfile\",$outfile);\n    \
myexit ($EXIT_SUCCESS);\n    }\n\n\n    \n\n\n\nsu\
b file2head\n      {\n	my $file = shift;\n	my $siz\
e = shift;\n	my $f= new FileHandle;\n	my $line;\n	\
open ($f,$file);\n	read ($f,$line, $size);\n	close\
 ($f);\n	return $line;\n      }\nsub file2tail\n  \
    {\n	my $file = shift;\n	my $size = shift;\n	my\
 $f= new FileHandle;\n	my $line;\n	\n	open ($f,$fi\
le);\n	seek ($f,$size*-1, 2);\n	read ($f,$line, $s\
ize);\n	close ($f);\n	return $line;\n      }\n\n\n\
sub vtmpnam\n      {\n	my $r=rand(100000);\n	my $f\
=\"file.$r.$$\";\n	while (-e $f)\n	  {\n	    $f=vt\
mpnam();\n	  }\n	push (@TMPFILE_LIST, $f);\n	retur\
n $f;\n      }\n\nsub myexit\n  {\n    my $code=@_\
[0];\n    if ($CLEAN_EXIT_STARTED==1){return;}\n  \
  else {$CLEAN_EXIT_STARTED=1;}\n    ### ONLY BARE\
 EXIT\n    exit ($code);\n  }\nsub set_error_lock\\
n    {\n      my $name = shift;\n      my $pid=$$;\
\n\n      \n      &lock4tc ($$,\"LERROR\", \"LSET\\
", \"$$ -- ERROR: $name $PROGRAM\\n\");\n      ret\
urn;\n    }\nsub set_lock\n  {\n    my $pid=shift;\
\n    my $msg= shift;\n    my $p=getppid();\n    &\
lock4tc ($pid,\"LLOCK\",\"LRESET\",\"$p$msg\\n\");\
\n  }\nsub unset_lock\n   {\n     \n    my $pid=sh\
ift;\n    &lock4tc ($pid,\"LLOCK\",\"LRELEASE\",\"\
\");\n  }\nsub shift_lock\n  {\n    my $from=shift\
;\n    my $to=shift;\n    my $from_type=shift;\n  \
  my $to_type=shift;\n    my $action=shift;\n    m\
y $msg;\n    \n    if (!&lock4tc($from, $from_type\
, \"LCHECK\", \"\")){return 0;}\n    $msg=&lock4tc\
 ($from, $from_type, \"LREAD\", \"\");\n    &lock4\
tc ($from, $from_type,\"LRELEASE\", $msg);\n    &l\
ock4tc ($to, $to_type, $action, $msg);\n    return\
;\n  }\nsub isshellpid\n  {\n    my $p=shift;\n   \
 if (!lock4tc ($p, \"LLOCK\", \"LCHECK\")){return \
0;}\n    else\n      {\n	my $c=lock4tc($p, \"LLOCK\
\", \"LREAD\");\n	if ( $c=~/-SHELL-/){return 1;}\n\
      }\n    return 0;\n  }\nsub isrootpid\n  {\n \
   if(lock4tc (getppid(), \"LLOCK\", \"LCHECK\")){\
return 0;}\n    else {return 1;}\n  }\nsub lock4tc\
\n	{\n	  my ($pid,$type,$action,$value)=@_;\n	  my\
 $fname;\n	  my $host=hostname;\n	  \n	  if ($type\
 eq \"LLOCK\"){$fname=\"$LOCKDIR/.$pid.$host.lock4\
tcoffee\";}\n	  elsif ( $type eq \"LERROR\"){ $fna\
me=\"$LOCKDIR/.$pid.$host.error4tcoffee\";}\n	  el\
sif ( $type eq \"LWARNING\"){ $fname=\"$LOCKDIR/.$\
pid.$host.warning4tcoffee\";}\n	  \n	  if ($debug_\
lock)\n	    {\n	      print STDERR \"\\n\\t---lock\
4tc(tcg): $action => $fname =>$value (RD: $LOCKDIR\
)\\n\";\n	    }\n\n	  if    ($action eq \"LCHECK\"\
) {return -e $fname;}\n	  elsif ($action eq \"LREA\
D\"){return file2string($fname);}\n	  elsif ($acti\
on eq \"LSET\") {return string2file ($value, $fnam\
e, \">>\");}\n	  elsif ($action eq \"LRESET\") {re\
turn string2file ($value, $fname, \">\");}\n	  els\
if ($action eq \"LRELEASE\") \n	    {\n	      if (\
 $debug_lock)\n		{\n		  my $g=new FileHandle;\n		 \
 open ($g, \">>$fname\");\n		  print $g \"\\nDestr\
oyed by $$\\n\";\n		  close ($g);\n		  safe_system\
 (\"mv $fname $fname.old\");\n		}\n	      else\n		\
{\n		  unlink ($fname);\n		}\n	    }\n	  return \"\
\";\n	}\n	\nsub file2string\n	{\n	  my $file=@_[0]\
;\n	  my $f=new FileHandle;\n	  my $r;\n	  open ($\
f, \"$file\");\n	  while (<$f>){$r.=$_;}\n	  close\
 ($f);\n	  return $r;\n	}\nsub string2file \n    {\
\n    my ($s,$file,$mode)=@_;\n    my $f=new FileH\
andle;\n    \n    open ($f, \"$mode$file\");\n    \
print $f  \"$s\";\n    close ($f);\n  }\n\nBEGIN\n\
    {\n      srand;\n    \n      $SIG{'SIGUP'}='si\
gnal_cleanup';\n      $SIG{'SIGINT'}='signal_clean\
up';\n      $SIG{'SIGQUIT'}='signal_cleanup';\n   \
   $SIG{'SIGILL'}='signal_cleanup';\n      $SIG{'S\
IGTRAP'}='signal_cleanup';\n      $SIG{'SIGABRT'}=\
'signal_cleanup';\n      $SIG{'SIGEMT'}='signal_cl\
eanup';\n      $SIG{'SIGFPE'}='signal_cleanup';\n \
     \n      $SIG{'SIGKILL'}='signal_cleanup';\n  \
    $SIG{'SIGPIPE'}='signal_cleanup';\n      $SIG{\
'SIGSTOP'}='signal_cleanup';\n      $SIG{'SIGTTIN'\
}='signal_cleanup';\n      $SIG{'SIGXFSZ'}='signal\
_cleanup';\n      $SIG{'SIGINFO'}='signal_cleanup'\
;\n      \n      $SIG{'SIGBUS'}='signal_cleanup';\\
n      $SIG{'SIGALRM'}='signal_cleanup';\n      $S\
IG{'SIGTSTP'}='signal_cleanup';\n      $SIG{'SIGTT\
OU'}='signal_cleanup';\n      $SIG{'SIGVTALRM'}='s\
ignal_cleanup';\n      $SIG{'SIGUSR1'}='signal_cle\
anup';\n\n\n      $SIG{'SIGSEGV'}='signal_cleanup'\
;\n      $SIG{'SIGTERM'}='signal_cleanup';\n      \
$SIG{'SIGCONT'}='signal_cleanup';\n      $SIG{'SIG\
IO'}='signal_cleanup';\n      $SIG{'SIGPROF'}='sig\
nal_cleanup';\n      $SIG{'SIGUSR2'}='signal_clean\
up';\n\n      $SIG{'SIGSYS'}='signal_cleanup';\n  \
    $SIG{'SIGURG'}='signal_cleanup';\n      $SIG{'\
SIGCHLD'}='signal_cleanup';\n      $SIG{'SIGXCPU'}\
='signal_cleanup';\n      $SIG{'SIGWINCH'}='signal\
_cleanup';\n      \n      $SIG{'INT'}='signal_clea\
nup';\n      $SIG{'TERM'}='signal_cleanup';\n     \
 $SIG{'KILL'}='signal_cleanup';\n      $SIG{'QUIT'\
}='signal_cleanup';\n      \n      our $debug_lock\
=$ENV{\"DEBUG_LOCK\"};\n      \n      \n      \n  \
    \n      foreach my $a (@ARGV){$CL.=\" $a\";}\n\
      if ( $debug_lock ){print STDERR \"\\n\\n\\n*\
********* START PG: $PROGRAM *************\\n\";}\\
n      if ( $debug_lock ){print STDERR \"\\n\\n\\n\
**********(tcg) LOCKDIR: $LOCKDIR $$ *************\
\\n\";}\n      if ( $debug_lock ){print STDERR \"\\
\n --- $$ -- $CL\\n\";}\n      \n	     \n      \n \
     \n    }\nsub flush_error\n  {\n    my $msg=sh\
ift;\n    return add_error ($EXIT_FAILURE,$$, $$,g\
etppid(), $msg, $CL);\n  }\nsub add_error \n  {\n \
   my $code=shift;\n    my $rpid=shift;\n    my $p\
id=shift;\n    my $ppid=shift;\n    my $type=shift\
;\n    my $com=shift;\n    \n    $ERROR_DONE=1;\n \
   lock4tc ($rpid, \"LERROR\",\"LSET\",\"$pid -- E\
RROR: $type\\n\");\n    lock4tc ($$, \"LERROR\",\"\
LSET\", \"$pid -- COM: $com\\n\");\n    lock4tc ($\
$, \"LERROR\",\"LSET\", \"$pid -- STACK: $ppid -> \
$pid\\n\");\n   \n    return $code;\n  }\nsub add_\
warning \n  {\n    my $rpid=shift;\n    my $pid =s\
hift;\n    my $command=shift;\n    my $msg=\"$$ --\
 WARNING: $command\\n\";\n    print STDERR \"$msg\\
";\n    lock4tc ($$, \"LWARNING\", \"LSET\", $msg)\
;\n  }\n\nsub signal_cleanup\n  {\n    print dtder\
r \"\\n**** $$ (tcg) was killed\\n\";\n    &cleanu\
p;\n    exit ($EXIT_FAILURE);\n  }\nsub clean_dir\\
n  {\n    my $dir=@_[0];\n    if ( !-d $dir){retur\
n ;}\n    elsif (!($dir=~/tmp/)){return ;}#safety \
check 1\n    elsif (($dir=~/\\*/)){return ;}#safet\
y check 2\n    else\n      {\n	`rm -rf $dir`;\n   \
   }\n    return;\n  }\nsub cleanup\n  {\n    #pri\
nt stderr \"\\n----tc: $$ Kills $PIDCHILD\\n\";\n \
   #kill (SIGTERM,$PIDCHILD);\n    my $p=getppid()\
;\n    $CLEAN_EXIT_STARTED=1;\n    \n    \n    \n \
   if (&lock4tc($$,\"LERROR\", \"LCHECK\", \"\"))\\
n      {\n	my $ppid=getppid();\n	if (!$ERROR_DONE)\
 \n	  {\n	    &lock4tc($$,\"LERROR\", \"LSET\", \"\
$$ -- STACK: $p -> $$\\n\");\n	    &lock4tc($$,\"L\
ERROR\", \"LSET\", \"$$ -- COM: $CL\\n\");\n	  }\n\
      }\n    my $warning=&lock4tc($$, \"LWARNING\"\
, \"LREAD\", \"\");\n    my $error=&lock4tc($$,  \\
"LERROR\", \"LREAD\", \"\");\n    #release error a\
nd warning lock if root\n    \n    if (isrootpid()\
 && ($warning || $error) )\n      {\n	\n	print STD\
ERR \"**************** Summary *************\\n$er\
ror\\n$warning\\n\";\n\n	&lock4tc($$,\"LERROR\",\"\
RELEASE\",\"\");\n	&lock4tc($$,\"LWARNING\",\"RELE\
ASE\",\"\");\n      } \n    \n    \n    foreach my\
 $f (@TMPFILE_LIST)\n      {\n	if (-e $f){unlink (\
$f);} \n      }\n    foreach my $d (@TMPDIR_LIST)\\
n      {\n	clean_dir ($d);\n      }\n    #No More \
Lock Release\n    #&lock4tc($$,\"LLOCK\",\"LRELEAS\
E\",\"\"); #release lock \n\n    if ( $debug_lock \
){print STDERR \"\\n\\n\\n********** END PG: $PROG\
RAM ($$) *************\\n\";}\n    if ( $debug_loc\
k ){print STDERR \"\\n\\n\\n**********(tcg) LOCKDI\
R: $LOCKDIR $$ *************\\n\";}\n  }\nEND \n  \
{\n    \n    &cleanup();\n  }\n   \n\nsub safe_sys\
tem \n{\n  my $com=shift;\n  my $ntry=shift;\n  my\
 $ctry=shift;\n  my $pid;\n  my $status;\n  my $pp\
id=getppid();\n  if ($com eq \"\"){return 1;}\n  \\
n  \n\n  if (($pid = fork ()) < 0){return (-1);}\n\
  if ($pid == 0)\n    {\n      set_lock($$, \" -SH\
ELL- $com (tcg)\");\n      exec ($com);\n    }\n  \
else\n    {\n      lock4tc ($$, \"LLOCK\", \"LSET\\
", \"$pid\\n\");#update parent\n      $PIDCHILD=$p\
id;\n    }\n  if ($debug_lock){printf STDERR \"\\n\
\\t .... safe_system (fasta_seq2hmm)  p: $$ c: $pi\
d COM: $com\\n\";}\n\n  waitpid ($pid,WTERMSIG);\n\
\n  shift_lock ($pid,$$, \"LWARNING\",\"LWARNING\"\
, \"LSET\");\n\n  if ($? == $EXIT_FAILURE || lock4\
tc($pid, \"LERROR\", \"LCHECK\", \"\"))\n    {\n  \
    if ($ntry && $ctry <$ntry)\n	{\n	  add_warning\
 ($$,$$,\"$com failed [retry: $ctry]\");\n	  lock4\
tc ($pid, \"LRELEASE\", \"LERROR\", \"\");\n	  ret\
urn safe_system ($com, $ntry, ++$ctry);\n	}\n     \
 elsif ($ntry == -1)\n	{\n	  if (!shift_lock ($pid\
, $$, \"LERROR\", \"LWARNING\", \"LSET\"))\n	    {\
\n	      add_warning ($$,$$,\"$com failed\");\n	  \
  }\n	  else\n	    {\n	      lock4tc ($pid, \"LREL\
EASE\", \"LERROR\", \"\");\n	    }\n	  return $?;}\
\n      else\n	{\n	  if (!shift_lock ($pid,$$, \"L\
ERROR\",\"LERROR\", \"LSET\"))\n	    {\n	      mye\
xit(add_error ($EXIT_FAILURE,$$,$pid,getppid(), \"\
UNSPECIFIED system\", $com));\n	    }\n	}\n    }\n\
  return $?;\n}\n\nsub check_configuration \n    {\
\n      my @l=@_;\n      my $v;\n      foreach my \
$p (@l)\n	{\n	  \n	  if   ( $p eq \"EMAIL\")\n	   \
 { \n	      if ( !($EMAIL=~/@/))\n		{\n		add_warni\
ng($$,$$,\"Could Not Use EMAIL\");\n		myexit(add_e\
rror ($EXIT_FAILURE,$$,$$,getppid(),\"EMAIL\",\"$C\
L\"));\n	      }\n	    }\n	  elsif( $p eq \"INTERN\
ET\")\n	    {\n	      if ( !&check_internet_connec\
tion())\n		{\n		  myexit(add_error ($EXIT_FAILURE,\
$$,$$,getppid(),\"INTERNET\",\"$CL\"));\n		}\n	   \
 }\n	  elsif( $p eq \"wget\")\n	    {\n	      if (\
!&pg_is_installed (\"wget\") && !&pg_is_installed \
(\"curl\"))\n		{\n		  myexit(add_error ($EXIT_FAIL\
URE,$$,$$,getppid(),\"PG_NOT_INSTALLED:wget\",\"$C\
L\"));\n		}\n	    }\n	  elsif( !(&pg_is_installed \
($p)))\n	    {\n	      myexit(add_error ($EXIT_FAI\
LURE,$$,$$,getppid(),\"PG_NOT_INSTALLED:$p\",\"$CL\
\"));\n	    }\n	}\n      return 1;\n    }\nsub pg_\
is_installed\n  {\n    my @ml=@_;\n    my $r, $p, \
$m;\n    my $supported=0;\n    \n    my $p=shift (\
@ml);\n    if ($p=~/::/)\n      {\n	if (safe_syste\
m (\"perl -M$p -e 1\")==$EXIT_SUCCESS){return 1;}\\
n	else {return 0;}\n      }\n    else\n      {\n	$\
r=`which $p 2>/dev/null`;\n	if ($r eq \"\"){return\
 0;}\n	else {return 1;}\n      }\n  }\n\n\n\nsub c\
heck_internet_connection\n  {\n    my $internet;\n\
    my $tmp;\n    &check_configuration ( \"wget\")\
; \n    \n    $tmp=&vtmpnam ();\n    \n    if     \
(&pg_is_installed    (\"wget\")){`wget www.google.\
com -O$tmp >/dev/null 2>/dev/null`;}\n    elsif  (\
&pg_is_installed    (\"curl\")){`curl www.google.c\
om -o$tmp >/dev/null 2>/dev/null`;}\n    \n    if \
( !-e $tmp || -s $tmp < 10){$internet=0;}\n    els\
e {$internet=1;}\n    if (-e $tmp){unlink $tmp;}\n\
\n    return $internet;\n  }\nsub check_pg_is_inst\
alled\n  {\n    my @ml=@_;\n    my $r=&pg_is_insta\
lled (@ml);\n    if (!$r && $p=~/::/)\n      {\n	p\
rint STDERR \"\\nYou Must Install the perl package\
 $p on your system.\\nRUN:\\n\\tsudo perl -MCPAN -\
e 'install $pg'\\n\";\n      }\n    elsif (!$r)\n \
     {\n	myexit(flush_error(\"\\nProgram $p Suppor\
ted but Not Installed on your system\"));\n      }\
\n    else\n      {\n	return 1;\n      }\n  }\n\n$\
program=\"T-COFFEE (Version_8.95)\";\n\n","*TC_MET\
HOD_FORMAT_01\n******************generic_method.tc\
_method*************\n*\n*       Incorporating new\
 methods in T-Coffee\n*       Cedric Notredame 26/\
08/08\n*\n****************************************\
***************\n*This file is a method file\n*Cop\
y it and adapt it to your need so that the method \
\n*you want to use can be incorporated within T-Co\
ffee\n********************************************\
***********\n*                  USAGE             \
                 *\n******************************\
*************************\n*This file is passed to\
 t_coffee via -in:\n*\n*	t_coffee -in Mgeneric_met\
hod.method\n*\n*	The method is passed to the shell\
 using the following\n*call:\n*<EXECUTABLE><PARAM1\
><IN_FLAG><seq_file><PARAM2><OUT_FLAG><outname><PA\
RAM>\n*\n*Conventions:\n*<FLAG_NAME> 	<TYPE>		<VAL\
UE>\n*<VALUE>:	no_name 	<=> Replaced with a space\\
n*<VALUE>:	&nbsp	<=> Replaced with a space\n*\n***\
**************************************************\
**\n*                  ALN_MODE                   \
        *\n***************************************\
****************\n*pairwise   ->all Vs all (no sel\
f )[(n2-n)/2aln]\n*m_pairwise ->all Vs all (no sel\
f)[n^2-n]^2\n*s_pairwise ->all Vs all (self): [n^2\
-n]/2 + n\n*multiple   ->All the sequences in one \
go\n*\nALN_MODE		pairwise\n*\n********************\
***********************************\n*            \
      OUT_MODE                           *\n******\
*************************************************\\
n* mode for the output:\n*External methods: \n* al\
n -> alignmnent File (Fasta or ClustalW Format)\n*\
 lib-> Lib file (TC_LIB_FORMAT_01)\n*Internal Meth\
ods:\n* fL -> Internal Function returning a List (\
Librairie)\n* fA -> Internal Function returning an\
 Alignmnent\n*\nOUT_MODE		aln\n*******************\
************************************\n*           \
       SEQ_TYPE                           *\n*****\
**************************************************\
\n*G: Genomic, S: Sequence, P: PDB, R: Profile\n*E\
xamples:\n*SEQTYPE	S	sequences against sequences (\
default)\n*SEQTYPE	S_P	sequence against structure\\
n*SEQTYPE	P_P	structure against structure\n*SEQTYP\
E	PS	mix of sequences and structure	\n*\nSEQ_TYPE	\
S\n*\n\n******************************************\
*************\n*                COMMAND LINE      \
                   *\n*EXECUTABLE PARAM1 IN_FLAG O\
UT_FLAG PARAM             *\n*********************\
**********************************\n**************\
*****************************************\n*      \
            EXECUTABLE                         *\n\
**************************************************\
*****\n*name of the executable\n*passed to the she\
ll: executable\n*	\nEXECUTABLE	tc_generic_method.p\
l\n*\n********************************************\
***********\n*                  IN_FLAG           \
                  *\n*****************************\
**************************\n*IN_FLAG\n*flag indica\
ting the name of the in coming sequences\n*IN_FLAG\
 S no_name ->no flag\n*IN_FLAG S &bnsp-in&bnsp -> \
\" -in \"\n*\nIN_FLAG		-infile=\n*\n**************\
*****************************************\n*      \
            OUT_FLAG                           *\n\
**************************************************\
*****\n*OUT_FLAG\n*flag indicating the name of the\
 out-coming data\n*same conventions as IN_FLAG\n*O\
UT_FLAG	S no_name ->no flag\n*if you want to redir\
ect, pass the parameters via PARAM1\n*set OUT_FLAG\
 to >\n*\nOUT_FLAG		-outfile=\n*\n****************\
***************************************\n*        \
          PARAM_1                              *\n\
**************************************************\
*****\n*<EXECUTABLE><PARAM1><IN_FLAG><seq_file><PA\
RAM2><OUT_FLAG><outname><PARAM>\n*Parameters sent \
to the EXECUTABLE and specified *before* IN_FLAG \\
n*If there is more than 1 PARAM line, the lines ar\
e\n*concatenated\n*Command_line: @EP@PARAM@-gapope\
n%e10%s-gapext%e20\n*	%s white space\n*	%e equal s\
ign\n*\n*PARAM1	\n*\n*\n*\n***********************\
********************************\n*               \
   PARAM_2                              *\n*******\
************************************************\n\
*<EXECUTABLE><PARAM1><IN_FLAG><seq_file><PARAM2><O\
UT_FLAG><outname><PARAM>\n*Parameters sent to the \
EXECUTABLE and specified \n*after* IN_FLAG and *be\
fore* OUT_FLAG\n*If there is more than 1 PARAM lin\
e, the lines are\n*concatenated\n*\n*PARAM1	\n*\n*\
\n************************************************\
*******\n*                  PARAM                 \
             *\n**********************************\
*********************\n*<EXECUTABLE><PARAM1><IN_FL\
AG><seq_file><PARAM2><OUT_FLAG><outname><PARAM>\n*\
Parameters sent to the EXECUTABLE and specified *a\
fter* OUT_FLAG\n*If there is more than 1 PARAM lin\
e, the lines are\n*concatenated\n*\nPARAM	-mode=se\
q_msa -method=clustalw\nPARAM   -OUTORDER=INPUT -N\
EWTREE=core -align -gapopen=-15\n*\n**************\
*****************************************\n*      \
            END                                *\n\
**************************************************\
*****\n","*TC_METHOD_FORMAT_01\n***************clu\
stalw_method.tc_method*********\nEXECUTABLE	clusta\
lw\nALN_MODE		pairwise\nIN_FLAG		-INFILE=\nOUT_FLA\
G		-OUTFILE=\nOUT_MODE		aln\nPARAM		-gapopen=-10\n\
SEQ_TYPE		S\n*************************************\
************\n","$VersionTag =                    \
                                                  \
                                                  \
           2.43;\nuse Env;\nuse FileHandle;\nuse C\
wd;\nuse File::Path;\nuse Sys::Hostname;\nour $PID\
CHILD;\nour $ERROR_DONE;\nour @TMPFILE_LIST;\nour \
$EXIT_FAILURE=1;\nour $EXIT_SUCCESS=0;\n\nour $REF\
DIR=getcwd;\nour $EXIT_SUCCESS=0;\nour $EXIT_FAILU\
RE=1;\n\nour $PROGRAM=\"extract_from_pdb\";\nour $\
CL=$PROGRAM;\n\nour $CLEAN_EXIT_STARTED;\nour $deb\
ug_lock=$ENV{\"DEBUG_LOCK\"};\nour $LOCKDIR=$ENV{\\
"LOCKDIR_4_TCOFFEE\"};\nif (!$LOCKDIR){$LOCKDIR=ge\
tcwd();}\nour $ERRORDIR=$ENV{\"ERRORDIR_4_TCOFFEE\\
"};\nour $ERRORFILE=$ENV{\"ERRORFILE_4_TCOFFEE\"};\
\n&set_lock ($$);\nif (isshellpid(getppid())){lock\
4tc(getppid(), \"LLOCK\", \"LSET\", \"$$\\n\");}\n\
      \nour $SILENT=\" >/dev/null 2>/dev/null\";\n\
our $INTERNET=-1;\n\n\n\n\n\n\nour $BLAST_MAX_NRUN\
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
.......name...[Returns 1 if the name is a PDB ID, \
0 otherwise]\n           -get_pdb_chains.....name.\
..[Returns the list of chains corresponding to the\
 entry]\n           -get_pdb_id.........name...[Re\
turns the PDB id within the provided pdb file]\n  \
         -get_fugue_name.....name...[Turns a name \
into a name valid for fugue]\n                    \
                  [Uses the netaddress to do so]\n\
	   -chain......FIRST..........[Extract the first \
chain only]\n		       A B C..........[Extract Seve\
ral chains if needed]\n		       ALL............[Ex\
tract all the chains]	\n           -ligand.....ALL\
............[Extract the ligands in the chain (HET\
ATM)]\n                       <name1>,<name2>[Extr\
act All the named lignds]\n	   -ligand_only.......\
........[Extract only the ligands]\n           -li\
gand_list...............[Extract the list of ligan\
ds]\n	   -coor.......<start>..<end>.[Coordinates o\
f the fragment to extract]\n			              [Omit\
 end to include the Cter]\n           -num........\
absolute.......[absolute: relative to the seq] \n \
                      file...........[file: relati\
ve to file]\n           -num_out....new...........\
.[new: start 1->L]\n                       old....\
........[old: keep the file coordinates]\n        \
   -delete.....<start>..<end>.[Delete from residue\
 start to residue end]\n	   -atom.......CA........\
.....[Atoms to include, ALL for all of them]\n		  \
     CA O N.........[Indicate several atoms if nee\
ded]\n	   -code.......3..............[Use the 1 le\
tter code or the 3 letters code]\n	   -mode.......\
raw............[Output original pdb file]\n       \
                pdb............[Output something t\
hat looks like pdb]\n		       fasta..........[Outp\
ut the sequences in fasta format]\n		       simple\
.........[Output a format easy to parse in C ]\n  \
          -seq_field..ATOM...........[Field used t\
o extract the sequence]\n		       SEQRES.........[\
Use the complete sequence]\n	   -seq..............\
.........[Equivalent to  -mode fasta]\n	   -model.\
.....1..............[Chosen Model in an NMR file]\\
n           -nodiagnostic..............[Switches E\
rror Messages off]\n           -debug.............\
........[Sets the DEBUG ON]\n           -no_remote\
_pdb_dir.........[Do not look for a remote file]\n\
           -cache_pdb.................[Cache Value\
, default is $HOME/.t_coffee/cache, other values: \
NO<=> No cache]\n\n      Environement Variables\n \
          These variables can be set from the envi\
ronement\n           Command line values with the \
corresponding flag superseed evironement value\n  \
         NO_REMOTE_PDB_DIR..........[Prevents the \
program from searching remote file: faster]\n     \
      PDB_DIR....................[Indicates where \
PDB file must be fetched (localy)]\n\n	 PROBLEMS: \
please contact cedric.notredame\\@europe.com\\n\";\
\n	 exit ($EXIT_SUCCESS);\n}\n\n$np=0;\n$n_para=$#\
ARGV;\n$model=1;\n$pdb_dir=$ENV{'PDB_DIR'};if ($pd\
b_dir){$pdb_dir.=\"/\";}\n$debug=$ENV{'DEBUG_EXTRA\
CT_FROM_PDB'};\n\n$no_remote_pdb_dir=$ENV{NO_REMOT\
E_PDB_DIR};\n$HOME=$ENV{'HOME'};\nif ( $ENV{CACHE_\
4_TCOFFEE})\n{$cache=$ENV{CACHE_4_TCOFFEE};}\nelse\
\n{\n    $cache=\"$HOME/.t_coffee/cache/\";\n}\n\n\
   \n$netaddress=\"http://www.rcsb.org/pdb/files/%\
%.pdb.gz\";\n$netcompression_pg=\"gunzip\";\n$netc\
ompression=\"gz\";\n\n  foreach ($np=0; $np<=$n_pa\
ra; $np++)\n{        \n    $value=$ARGV[$np];\n   \
 \n    if  ($np==0 && !($value=~/^-.*/))\n{ \n    \
   $pdb_file= $ARGV[$np];\n}\n    elsif ( !($value\
=~/^-.*/))\n{\n	print \"@ARGV\";\n	die;\n} \n    \\
n    elsif ($value eq \"-nodiagnostic\"){$nodiagno\
stic=1;}\n    elsif ($value eq \"-force\")\n{\n	$f\
orce_pdb=1;\n}\n    elsif ($value eq \"-force_name\
\")\n{\n	$force_name=$ARGV[++$np];\n	$force_pdb=1;\
\n}\n    \n    elsif ($value eq \"-is_pdb_name\")\\
n{\n	$pdb_file= $ARGV[++$np];\n	\n	$is_pdb_name=1;\
\n	\n} \n    elsif ($value eq \"-debug\")\n{\n	$de\
bug=1;\n}\n    elsif ($value eq \"-get_pdb_chains\\
")\n{\n	$pdb_file= $ARGV[++$np];\n	$get_pdb_chains\
=1;\n}\n    elsif ($value eq \"-get_pdb_ligands\")\
\n{\n	$get_pdb_ligands=1;\n}\n    \n    elsif ($va\
lue eq \"-get_pdb_id\")\n{\n	$pdb_file= $ARGV[++$n\
p];\n	$get_pdb_id=1;\n	\n}\n    \n    elsif ( $val\
ue eq \"-get_fugue_name\")\n{\n	$pdb_file= $ARGV[+\
+$np];\n	$get_fugue_name=1;\n}\n    elsif ( $value\
 eq \"-infile\")\n{\n       $pdb_file= $ARGV[++$np\
];\n} \n    elsif ($value eq \"-netfile\")\n{\n	$n\
etfile=1;\n	if ( !($ARGV[$np+1]=~/^-.*/)){$pdb_fil\
e= $ARGV[++$np];}\n}\n    elsif (  $value eq \"-nu\
m\")\n{\n       $numbering= $ARGV[++$np];\n}\n    \
elsif (  $value eq \"-num_out\")\n{\n       $numbe\
ring_out= $ARGV[++$np];\n}\n    elsif ( $value eq \
\"-netaddress\")\n{\n	$netadress=$ARGV[++$np];\n}\\
n     \n    elsif ( $value eq \"-netcompression\")\
\n{\n	 $netcompression=$ARGV[++$np];\n}\n    elsif\
 ( $value eq \"-pdb_dir\")\n{\n	 if ( !($ARGV[$np+\
1]=~/^-.*/)){$pdb_dir= \"$ARGV[++$np]/\";}\n}\n   \
  elsif ( $value eq \"-no_remote_pdb_dir\")\n{\n	$\
no_remote_pdb_dir=1;\n	if ( !($ARGV[$np+1]=~/^-.*/\
)){$pdb_dir= \"$ARGV[++$np]/\";}\n}\n    elsif ( $\
value eq \"-cache\")\n{\n	$cache=$ARGV[++$np];\n}\\
n    \n    elsif ($value eq \"-netcompression_pg\"\
)\n{\n	  $netcompression_pg=$ARGV[++$np];\n}\n    \
 elsif ($value eq \"-mode\")\n{\n       $MODE=$ARG\
V[++$np];\n}\n\n    elsif ( $value eq \"-model\")\\
n{\n       $model= $ARGV[++$np];\n}\n    elsif ($v\
alue eq \"-seq_field\" )\n{\n       $seq_field= $A\
RGV[++$np];\n}   \n    elsif ($value eq \"-coor\" \
)\n{\n       $start= $ARGV[++$np];\n  \n       if \
(($ARGV[$np+1] eq \"\") ||($ARGV[$np+1]=~/^-.*/)){\
$end=\"*\";} \n       else {$end=   $ARGV[++$np];}\
     \n       $coor_set=1;\n}\n    elsif ($value e\
q \"-delete\" )\n{\n       $delete_start= $ARGV[++\
$np];\n       $delete_end= $ARGV[++$np];\n       $\
delete_set=1;\n}\n    elsif  ($value eq \"-code\")\
\n{\n       $code= $ARGV[++$np];\n}\n    elsif  ($\
value eq \"-no_hetatm\")\n{\n       $no_hetatm=1;\\
n}\n    elsif ($value eq \"-chain\")\n{\n       wh\
ile (!($ARGV[$np+1] eq \"\") &&!($ARGV[$np+1]=~/^-\
.*/))\n{\n	      ++$np;\n	      @c_chain=(@chain, \
 $ARGV[$np]);\n	      $hc_chain{$ARGV[$np]}=$#c_ch\
ain+1;\n}           \n}\n    elsif ($value eq \"-a\
tom\")\n{\n\n       while (!($ARGV[$np+1] eq \"\")\
 && !($ARGV[$np+1]=~/^-.*/))\n{\n	      ++$np;\n	 \
     $atom[$n_atom++]=  $ARGV[$np];\n	      $atom_\
list{$ARGV[$np]}=1;	      \n} \n       \n}\n    el\
sif ( $value eq \"-unfold\")\n{\n	$unfold=1;\n}\n \
   elsif ($value eq \"-seq\" ||$value eq \"-fasta\\
" )\n{\n       $MODE=\"fasta\";\n}\n    elsif ( $v\
alue eq \"-version\")\n{\n	print STDERR  \"\\nextr\
act_from_pdb: Version $VersionTag\\n\";\n	&myexit \
($EXIT_SUCCESS);\n}\n    elsif ( $value eq \"-liga\
nd\")\n{\n	while (!($ARGV[$np+1] eq \"\") && !($AR\
GV[$np+1]=~/^-.*/))\n{\n	    ++$np;\n	    $ligand=\
1;\n	    $ligand_list{$ARGV[$np]}=1;	      \n} \n	\
$hc_chain{'LIGAND'}=1;\n}\n    elsif ( $value eq \\
"-ligand_only\")\n{\n	$ligand_only=1;\n}\n}\nif ( \
$debug)\n{\n    print STDERR \"\\n[DEBUG:extract_f\
rom_pdb] NO_REMOTE_PDB_DIR: $no_remote_pdb_dir\\n\\
";\n    print STDERR \"\\n[DEBUG:extract_from_pdb]\
 PDB_DIR: $pdb_dir\\n\";\n}\n\nif ( $is_pdb_name)\\
n{\n    if (remote_is_pdb_name($pdb_file, $netaddr\
ess))\n{\n	print \"1\";\n}\n    else\n{\n	print \"\
0\";\n}\n    exit ($EXIT_SUCCESS);\n}\n    \n\nif \
(!$force_name)\n{\n    $pdb_file=~/([^\\/]*)$/;\n \
   $force_name=$1;\n}\n\n$local_pdb_file=$pdb_file\
;\n\nif ( $debug){print STDERR \"\\n[DEBUG: extrac\
t_from_pdb] Scan For $local_pdb_file\\n\";}\n\n$me\
m=$no_remote_pdb_dir;\n$no_remote_pdb_dir=1;\n$tmp\
_pdb_file=get_pdb_file ($local_pdb_file);\n\nif ( \
!-e $tmp_pdb_file || $tmp_pdb_file eq \"\")\n  {\n\
    $local_pdb_file=$pdb_file;\n    ($local_pdb_fi\
le, $suffix_chain)=&pdb_name2name_and_chain($local\
_pdb_file);\n\n    if ($local_pdb_file)\n      {\n\
	if ( $debug){print STDERR \"\\nSplit $pdb_file in\
to $local_pdb_file and $suffix_chain \\n\";}\n	$tm\
p_pdb_file=get_pdb_file ($local_pdb_file);\n	if ( \
$tmp_pdb_file ne \"\")\n	  {\n	    @c_chain=();\n	\
    @c_chain=($suffix_chain);\n	    %hc_chain=();\\
n	    $hc_chain{$suffix_chain}=1;\n	  }\n      }\n\
  }\n\n$no_remote_pdb_dir=$mem;\nif ($no_remote_pd\
b_dir==0)\n{\n    if ( !-e $tmp_pdb_file || $tmp_p\
db_file eq \"\")\n{\n	\n	$local_pdb_file=$pdb_file\
;\n	($local_pdb_file, $suffix_chain)=&pdb_name2nam\
e_and_chain($local_pdb_file);\n	if ($local_pdb_fil\
e)\n{\n	    if ( $debug){print STDERR \"\\nSplit $\
pdb_file into $local_pdb_file and $suffix_chain \\\
n\";}\n	    $tmp_pdb_file=get_pdb_file ($local_pdb\
_file);    \n	    if ( $tmp_pdb_file ne \"\")\n{\n\
		@c_chain=();\n		@c_chain=($suffix_chain);\n		%hc\
_chain=();\n		$hc_chain{$suffix_chain}=1;\n}\n}\n}\
\n}\n\nif ( $debug){print STDERR \"\\n$pdb_file co\
pied into ##$tmp_pdb_file##\\n\";}\n\n\nif ( !-e $\
tmp_pdb_file || $tmp_pdb_file eq \"\")\n{\n	if ($i\
s_pdb_name)\n{\n	    print \"0\\n\"; exit ($EXIT_S\
UCCESS);\n}\n	else\n{\n  \n	    print \"\\nEXTRACT\
_FROM_PDB: NO RESULT for $pdb_file\\n\";\n	    &my\
exit ($EXIT_SUCCESS);	\n}\n}\n\n\n\n\n%molecule_ty\
pe=&pdbfile2chaintype($tmp_pdb_file);\nif ( $debug\
){print STDERR \"\\n\\tSequence Type determined\\n\
\";}\n\n$pdb_id=&get_pdb_id ($tmp_pdb_file);\nif (\
 $debug){print STDERR \"\\n\\tID: $pdb_id (for $tm\
p_pdb_file)\\n\";}\n\nif ( $pdb_id eq \"\"){$pdb_i\
d=$force_name;}\n\n@f_chain=&get_chain_list ($tmp_\
pdb_file);\nif ( $debug){print STDERR \"\\n\\tChai\
n_list:@f_chain\\n\";}\n\nif ( $get_pdb_chains)\n{\
\n    print \"@f_chain\\n\";\n    &myexit ($EXIT_S\
UCCESS);\n}\nif ( $get_pdb_ligands)\n{\n    %compl\
ete_ligand_list=&get_ligand_list ($tmp_pdb_file);\\
n    print $complete_ligand_list{\"result\"};\n   \
 &myexit ($EXIT_SUCCESS);\n}\n\nelsif ( $get_pdb_i\
d ||$get_fugue_name )\n{\n    if    (@c_chain && $\
c_chain[0] eq \"FIRST\"){$pdb_id=$pdb_id.$f_chain[\
0];}\n    elsif (@c_chain && $c_chain[0] ne \" \")\
{$pdb_id=$pdb_id.$c_chain[0];}\n    \n    print \"\
$pdb_id\\n\";\n    &myexit ($EXIT_SUCCESS);\n    \\
n}\nelsif ( $is_pdb_name)\n{\n    printf \"1\\n\";\
\n    &myexit ($EXIT_SUCCESS);\n}\n\n\n\n$structur\
e_file=vtmpnam();\n\nif ( $debug){print STDERR \"\\
\n\\tCheck_point #1: $tmp_pdb_file  $structure_fil\
e\\n\";}\n\n$INFILE=vfopen (\"$tmp_pdb_file\", \"r\
\"); \n$TMP=vfopen (\"$structure_file\", \"w\");\n\
\n$print_model=1;\n$in_model=0;\n\nif ( $debug){pr\
int STDERR \"\\n\\tCheck_point #2\\n\";}\nwhile ( \
<$INFILE>)\n{\n  my $first_model=0;\n  $line=$_;\n\
\n  if ( !$first_model && ($line =~/^MODEL\\s*(\\d\
*)/))\n    {\n      $first_model=$1;\n      if ($m\
odel==1){$model=$first_model;}\n    }\n  \n  if ((\
$line =~/^MODEL\\s*(\\d*)/))\n    {\n      if ($1=\
=$model)\n	{\n	  $in_model=1;\n	  $print_model=1;\\
n	  $is_nmr=1;\n	}\n      elsif ( $in_model==0)\n	\
{\n	  $print_model=0;\n	}\n      elsif ( $in_model\
==1)\n	{\n	  last;\n	}\n    }\n  if ($print_model)\
{print $TMP $line;}  \n}\nclose ($TMP);\nclose ($I\
NFILE);\n\nif ( $debug){print STDERR \"\\n\\tCheck\
_point #3\\n\";}	\n\n  if ($numbering eq \"\"){$nu\
mbering=\"absolute\";}\n  if ($numbering_out eq \"\
\"){$numbering_out=\"new\";}\n\n  if ( $delete_set\
 && $coor_set) {die \"-delete and -coor are mutual\
ly exclusive, sorry\\n\";}\n  if ( $n_atom==0){$at\
om_list[$n_atom++]=\"ALL\";$atom_list{$atom_list[0\
]}=1;}\n  if ( $seq_field eq \"\"){$seq_field=\"AT\
OM\";}\n  \n  if ( $MODE eq \"\"){$MODE=\"pdb\";}\\
n  elsif ( $MODE eq \"simple\" && $code==0){$code=\
1;}\n\n  if ( $code==0){$code=3;}\n\n\nif ($f_chai\
n[0] eq \" \"){$hc_chain{' '}=1;$c_chain[0]=\" \";\
}\nelsif (!@c_chain){$hc_chain{FIRST}=1;$c_chain[0\
]=\"FIRST\";}#make sure the first chain is taken b\
y default\n\nif    ($hc_chain{ALL}) \n{\n      @c_\
chain=@f_chain;\n      foreach $e (@c_chain){$hc_c\
hain{$e}=1;}\n}\nelsif($hc_chain{FIRST})\n{\n	@c_c\
hain=($f_chain[0]);\n	$hc_chain{$f_chain[0]}=1;\n}\
\n\n\n$MAIN_HOM_CODE=&get_main_hom_code ($structur\
e_file);\n$INFILE=vfopen ($structure_file, \"r\");\
\n\n\nif ( $MODE eq \"raw_pdb\" || $MODE eq \"raw\\
")\n{\n    while (<$INFILE>)\n{	print \"$_\";}\n  \
  close ( $INFILE);\n    &myexit($EXIT_SUCCESS);\n\
}    \nif ( $MODE eq \"raw4fugue\" )\n{\n    while\
 (<$INFILE>)\n{	\n	$l=$_;\n	if ($l=~/^SEQRES/)\n{\\
n	    \n	    $c= substr($l,11,1);\n	    if ($hc_ch\
ain {$c}){print \"$l\";}\n}\n	elsif ( $l=~/^ATOM/)\
\n{\n	    $c=substr($l,21,1);\n	    if ($hc_chain \
{$c}){print \"$l\";}\n}\n}\n    close ( $INFILE);\\
n    &myexit($EXIT_SUCCESS);\n}    \n\nif ( $MODE \
eq \"pdb\")\n{\n\n    $read_header=0;\n    while (\
<$INFILE>) \n{\n	    $line=$_;\n	    if    ($line \
=~ /^HEADER/){print \"$line\";$read_header=1;}\n}\\
n    close ($INFILE);\n\n    if (!$read_header)\n{\
\n	print \"HEADER    UNKNOWN                      \
           00-JAN-00   $force_name\\n\";\n}\n\n   \
 $INFILE=vfopen ($structure_file, \"r\");\n    \n \
   print \"COMPND   1 CHAIN:\";\n    $last=pop(@c_\
chain);\n    foreach $c ( @c_chain){ print \" $c,\\
";}\n    if ( $last eq \" \"){print \" NULL;\\n\";\
}\n    else \n{\n      print \" $last;\\n\";\n}\n \
   @c_chain=(@c_chain, $last);\n    \n    print \"\
REMARK Output of the program extract_from_pdb (Ver\
sion $VersionTag)\\n\";\n    print \"REMARK Legal \
PDB format not Guaranteed\\n\";\n    print \"REMAR\
K This format is not meant to be used in place of \
the PDB format\\n\";\n    print \"REMARK The heade\
r refers to the original entry\\n\";\n    print \"\
REMARK The sequence from the original file has bee\
n taken in the field: $seq_field\\n\";\n    print \
\"REMARK extract_from_pdb, 2001, 2002, 2003, 2004,\
 2005 2006 (c) CNRS and Cedric Notredame\\n\";   \\
n    if ( $coor_set)\n{\n       print \"REMARK Par\
tial chain: Start $start End $end\\n\";\n}\n    if\
 ( $is_nmr)\n{\n       print \"REMARK NMR structur\
e: MODEL $model\\n\";\n}\n   \n    if ( $n_atom!=0\
)\n{\n       print \"REMARK Contains Coordinates o\
f: \";\n       foreach $a (@atom){print \"$a \";}\\
n       print \"\\n\";\n}  \n}\n\n\n\n\nmy $residu\
e_index = -999;\nmy $old_c = \"TemporaryChain\";\n\
\nwhile (<$INFILE>) \n{\n	$line=$_;\n\n\n	if ($lin\
e =~ /^SEQRES/)\n{\n\n		@field=/(\\S*)\\s*/g;\n\n	\
	$c= substr($_,11,1);\n\n		\n		$l=$#field;\n		for \
($a=4; $a<$#field ;)\n{\n			if (!$onelett{$molecul\
e_type{$c}}->{$field[$a]})\n{\n				splice @field, \
$a, 1;\n}\n			else \n{\n				$a++;\n}\n}\n	\n		if (\
 $c ne $in_chain)\n{\n			$pdb_chain_list[$n_pdb_ch\
ains]=$c;\n			$pdb_chain_len [$n_pdb_chains]=$len;\
\n			$in_chain=$c;\n			$n_pdb_chains++;\n}\n	\n		f\
or ( $a=4; $a<$#field;$a++)\n{\n			@{$complete_seq\
{$c}}->[$complete_seq_len{$c}++]=$field[$a];\n}\n}\
\n    elsif ( $line=~/^ATOM/ || ($line=~/^HETATM/ \
&& &is_aa(substr($line,17,3),substr($line,21,1)) &\
& !$no_hetatm))\n{\n\n	 \n    $RAW_AT_ID=$AT_ID=su\
bstr($line,12,4);\n	$RES_ID=&is_aa(substr($line,17\
,3),substr($line,21,1));\n	$CHAIN=substr($line,21,\
1);\n\n    $RES_NO=substr($line,22,4);\n	$HOM_CODE\
=substr ($line, 26, 1);\n	$TEMP=substr($line,60,6)\
;\n	\n	$TEMP=~s/\\s//g;\n        $AT_ID=~s/\\s//g;\
\n	$RES_ID=~s/\\s//g;\n        $RES_NO=~s/\\s//g;\\
n		\n	if ( $HOM_CODE ne $MAIN_HOM_CODE){next;}\n	e\
lsif ( $already_read2{$CHAIN}{$RES_ID}{$AT_ID}{$RE\
S_NO}){next;}\n	else{$already_read2{$CHAIN}{$RES_I\
D}{$AT_ID}{$RES_NO}=1;}\n	\n	\n	if ($coor_set && $\
numbering eq \"file\" && $residue_index ne $RES_NO\
)\n{\n	    \n	    if ( $RES_NO<=$start){$real_star\
t{$CHAIN}++;}\n	    if ( $RES_NO<=$end){$real_end{\
$CHAIN}++;}\n}\n	elsif ($numbering eq \"absolute\"\
)\n{\n	    $real_start{$CHAIN}=$start;\n	    $real\
_end{$CHAIN}=$end;\n}\n\n        $KEY=\"ALL\";\n  \
      if ( $CHAIN ne $in_atom_chain)\n{\n	    \n	 \
 $pdb_atom_chain_list[$n_pdb_atom_chains]=$c;\n	  \
$pdb_atom_chain_len [$n_pdb_atom_chains]=$len;\n	 \
 $in_atom_chain=$c;\n	  $n_pdb_atom_chains++;\n}\n\
	\n	if ( $residue_index ne $RES_NO)\n{\n	     $res\
idue_index = $RES_NO;\n	     @{$atom_seq{$CHAIN}}-\
>[$atom_seq_len{$CHAIN}++]=$RES_ID;;\n}\n}\n\n}\nc\
lose ($INFILE);\n\n\n\n\n\n\n$INFILE=vfopen ($stru\
cture_file, \"r\");\nforeach $c (@c_chain)\n{\n\n	\
if    ( $seq_field eq \"SEQRES\"){@pdb_seq=@{$comp\
lete_seq{$c}};}\n	elsif ( $seq_field eq \"ATOM\") \
 {@pdb_seq=@{$atom_seq{$c}};}\n	\n\n	$full_length=\
$l=$#pdb_seq+1;\n		\n	if ( $real_end{$c}==\"*\"){$\
real_end{$c}=$full_length;}\n	if ( $coor_set)\n{	 \
  \n\n	   if ( $real_end{$c} < $l){splice @pdb_seq\
, $real_end{$c}, $l;}\n	   if ( $real_start{$c} < \
$l){splice @pdb_seq, 0, $real_start{$c}-1;}	  	   \
\n	   $l=$#pdb_seq;\n}\n\n	elsif ( $delete_set)\n{\
\n	   splice @pdb_seq, $delete_start, $delete_end-\
$delete_start+1;\n	   $l=$#pdb_seq;\n}\n	\n	$new_f\
asta_name=\"$pdb_id$c\";\n	if ( $coor_set)\n{\n	  \
 if ( $n_pdb_chains==0){$new_fasta_name=\"$new_fas\
ta_name$c\";}\n	   $new_fasta_name= $new_fasta_nam\
e.\"\\_$start\\_$end\";\n}\n	   \n	if ( $MODE eq \\
"pdb\")\n{\n	   $nl=1;\n	   $n=0;\n	   \n	   forea\
ch $res ( @pdb_seq)\n		{\n		if ( !$n)\n		{\n		\n		\
 printf \"SEQRES %3d %1s %4d  \", $nl,$c, $l;\n		 \
$nl++;\n	}\n	     $res=~s/\\s//g;\n	     \n	     i\
f ($code==1){ printf \"%3s \",$onelett{$molecule_t\
ype{$c}}->{$res};}\n	     elsif  ($code==3){ print\
f \"%3s \",$res};\n	     \n	     $n++;		  \n	     \
if ( $n==13){$n=0;print \"\\n\";}\n}\n	  if ( $n!=\
0){print \"\\n\"; $n=0;}\n}\n	elsif ( $MODE eq \"s\
imple\")\n{\n	  print \"# SIMPLE_PDB_FORMAT\\n\";\\
n	  if ( $new_fasta_name eq \" \"){$new_fasta_name\
=\"dummy_name\";}\n	  print \">$new_fasta_name\\n\\
";\n\n	  foreach $res ( @pdb_seq)\n{\n	      print\
 \"$onelett{$molecule_type{$c}}->{$res}\";\n}\n	  \
print \"\\n\";\n}\n	elsif ( $MODE eq \"fasta\")\n{\
\n	  $n=0;\n	  print \">$new_fasta_name\\n\";\n	  \
\n	  foreach $res ( @pdb_seq)\n{\n\n	      print \\
"$onelett{$molecule_type{$c}}->{$res}\";\n        \
      $n++;\n	      if ( $n==60){print \"\\n\"; $n\
=0;}\n}\n	  print \"\\n\"; \n}\n}\n\nif ( $MODE eq\
 \"fasta\")\n{\n     &myexit($EXIT_SUCCESS);\n  \n\
}\n\n  \n  $charcount=0;\n  $inchain=\"BEGIN\";\n \
 $n=0;\n  while (<$INFILE>) \n{\n    $line=$_;\n  \
   \n    if ($line =~/^ATOM/  ||  ($line=~/^HETATM\
/))\n{\n	$line_header=\"UNKNWN\";\n	$RES_ID=substr\
($line,17,3);\n	$chain = substr($line,21,1);\n\n	i\
f ($line =~/^ATOM/)\n{\n	    $line_header=\"ATOM\"\
;\n	    $RES_ID=(&is_aa($RES_ID,$chain))?&is_aa($R\
ES_ID,$chain):$RES_ID;\n}\n	elsif ($line=~/^HETATM\
/ && ($ligand_list {$RES_ID} ||$ligand_list {'ALL'\
} || !&is_aa($RES_ID,$chain)))\n{\n	    $line_head\
er=\"HETATM\";\n}\n	elsif ($line=~/^HETATM/ && (&i\
s_aa($RES_ID,$chain) && !$no_hetatm))\n{\n	    $li\
ne_header=\"ATOM\";\n	    $RES_ID=&is_aa($RES_ID,$\
chain);\n}\n	else\n{\n	    next;\n}\n\n	\n\n	$X=su\
bstr($line,30,8);     \n	$Y=substr($line,38,8);\n	\
$Z=substr($line,46,8);\n	$TEMP=substr($line,60,6);\
\n	\n	$RAW_AT_ID=$AT_ID=substr($line,12,4);\n	$CHA\
IN=substr($line,21,1);\n	$RES_NO=substr($line,22,4\
);\n	$HOM_CODE=substr ($line, 26, 1);\n	\n	$X=~s/\\
\s//g;\n	$Y=~s/\\s//g;\n	$Z=~s/\\s//g;\n	$TEMP=~s/\
\\s//g;\n	\n	$AT_ID=~s/\\s//g;\n	$RES_ID=~s/\\s//g\
;\n	$RES_NO=~s/\\s//g;\n\n	\n	if ( $HOM_CODE ne $M\
AIN_HOM_CODE){next;}\n	elsif ( $already_read{$CHAI\
N}{$RES_ID}{$AT_ID}{$RES_NO}){next;}\n	else{$alrea\
dy_read{$CHAIN}{$RES_ID}{$AT_ID}{$RES_NO}=1;}\n	\n\
	$KEY=\"ALL\";\n\n      	if ( $RES_NO ==0){$start_\
at_zero=1;}\n\n	$RES_NO+=$start_at_zero;    \n	\n	\
if ( $current_chain ne $CHAIN)\n{\n	    $current_c\
hain=$CHAIN;\n	    $pos=$current_residue=0;\n	    \
$offset=($coor_set)?($real_start{$CHAIN}-1):0;\n	 \
   if    ( $seq_field eq \"SEQRES\"){@ref_seq=@{$c\
omplete_seq{$CHAIN}};}\n	    elsif ( $seq_field eq\
 \"ATOM\")  {@ref_seq=@{$atom_seq{$CHAIN}};}\n}\n	\
\n	if ($current_residue != $RES_NO)\n{\n	    $curr\
ent_residue=$RES_NO;\n	    if    ( $seq_field eq \\
"SEQRES\"){$pos=$current_residue;}\n	    elsif ( $\
seq_field eq \"ATOM\"){$pos++;}\n}\n	\n	\n	if ($n_\
atom==0 || $atom_list{$AT_ID}==1 || $atom_list{$KE\
Y}==1)\n{ 	\n	    \n	    $do_it=(!@c_chain || $hc_\
chain{$CHAIN} ||$hc_chain{'LIGAND'} );\n	    \n	  \
  $do_it= ($do_it==1) && ($coor_set==0 ||($pos>=$r\
eal_start{$CHAIN} && $pos<=$real_end{$CHAIN}));\n	\
    $do_it= ($do_it==1) && ($delete_set==0 || $pos\
<$delete_start ||$pos>$delete_end );\n	    if ($li\
gand==0 && $line_header eq \"HETATM\" ){$do_it=0;}\
\n	    if ($ligand_only==1 && $line_header eq \"AT\
OM\" ){$do_it=0;}\n	    if ($ligand==1 && $line_he\
ader eq \"HETATM\" && $ligand_list{$RES_ID}==0 && \
$ligand_list{\"ALL\"}==0){$do_it=0;} \n	    \n	   \
 \n	    if ( $do_it)\n{\n		$n++;\n		$out_pos=$pos;\
\n		\n	       if ( $delete_set)\n{\n		  if ( $out_\
pos< $delete_start){;}\n		  else {$offset=$delete_\
end-$delete_start;}\n}       \n	       \n	       i\
f ( $numbering_out eq \"new\"){$out_pos-=$offset;}\
\n	       elsif ( $numbering_out eq \"old\"){$out_\
pos=$RES_NO;}\n	       \n       \n	       \n	     \
  if ( $code==1){$RES_ID=$onelett{$molecule_type{$\
c}}->{$RES_ID};}\n	    \n	       if ($unfold)\n{\n\
		   $unfolded_x+=5;\n		   $X=$unfolded_x;\n		   $\
Y=0;\n		   $Z=0;\n		   $float=1;\n}\n	       else\\
n{\n		   $float=3;\n}\n\n	       if ( $MODE eq \"p\
db\")\n{\n		   printf \"%-6s%5d %-4s %3s %s%4d    \
%8.3f%8.3f%8.3f  1.00 %5.2f\\n\",$line_header, $n,\
 $RAW_AT_ID,$RES_ID,$CHAIN,$out_pos, $X, $Y, $Z,$T\
EMP;		  \n}\n	       elsif ( $MODE eq \"simple\")\\
n{\n		    if ( $RES_ID eq \"\"){$RES_ID=\"X\";}\n	\
	  printf \"%-6s %5s %s %2s %4d    %8.3f %8.3f %8.\
3f\\n\",$line_header, $AT_ID, $RES_ID,($CHAIN eq\"\
\" || $CHAIN eq \" \")?\"A\":$CHAIN,$out_pos, $X, \
$Y, $Z,$TEMP;\n}\n\n}\n}\n}\n}\nprint \"\\n\";\ncl\
ose($INFILE);\n\n\nif ( $error ne \"\") \n{$error=\
$error.\"\\nDiagnostic:    SEQRES and the residues\
 in ATOM are probably Incompatible\\n\";\n    $err\
or=$error.  \"Recomendation: Rerun with '-fix 1' i\
n order to ignore the SEQRES sequences\\n\";\n}\ni\
f (!$nodiagnostic){print STDERR $error;}\n&myexit \
( $EXIT_SUCCESS);\n\nsub remote_is_pdb_name \n{\n \
   my $in=@_[0];\n    my $netaddress=@_[1];\n    m\
y $ref_file, $pdb;\n    my $value;\n\n    if ( $in\
=~/[^\\w\\d\\:\\_]/){return 0;}\n    \n    \n    $\
ref_file=\"$cache/pdb_entry_type.txt\";\n    \n   \
 if ( !-e $ref_file || (-M $ref_file)>2 || -z $ref\
_file)\n{\n	&url2file(\"ftp://ftp.wwpdb.org/pub/pd\
b/derived_data/pdb_entry_type.txt\", $ref_file);\n\
}\n \n    $pdb=substr ($in,0, 4);\n    \n    \n   \
 $value=`grep -c $pdb $ref_file`;\n    return $val\
ue;\n}\n      \nsub is_pdb_file\n{\n    my @arg=@_\
;\n\n    if ( !-e $arg[0]){return 0;}\n    \n    $\
F=vfopen ($arg[0], \"r\");\n    while ( <$F>)\n{\n\
	if (/^HEADER/)\n{\n	    close $F;\n	    return 1;\
\n}\n	elsif ( /^SEQRES/)\n{\n	    close $F;\n	    \
return 1;\n}\n	elsif ( /^ATOM/)\n{\n	    close $F;\
\n	    return 1;\n}\n}\n    return 0;\n}\nsub get_\
pdb_id\n{\n    my $header_file=@_[0];\n    my $id;\
\n    my $F= new FileHandle;\n    \n    \n    $F=v\
fopen (\"$header_file\", \"r\");\n\n    while ( <$\
F>)\n      {\n	if ( /HEADER/)\n	  {\n	    if ($deb\
ug){print \"$_\";}\n	    $id=substr($_,62,4 );\n	 \
   return $id;\n	  }\n      }\n    close ($F);\n  \
  \n    return \"\";\n}\n\nsub get_ligand_list\n{\\
n    my $pdb_file=@_[0];\n    my $chain;\n    my $\
ligand;\n    my %complete_ligand_list;\n    \n\n  \
  $F=vfopen ($pdb_file, \"r\");\n    while ( <$F>)\
\n{\n	if ( /^HETATM/)\n{\n	    $line=$_;\n	    $ch\
ain=substr($line,21,1);\n	    $ligand=substr($line\
,17,3);\n	    \n	    if (!$complete_ligand_list{$c\
hain}{$ligand})\n{\n		\n		$complete_ligand_list{\"\
result\"}.=\"CHAIN $chain LIGAND $ligand\\n\";\n		\
$complete_ligand_list{$chain}{$ligand}=1;\n}\n}\n}\
\n    close ($F);\n    return %complete_ligand_lis\
t;\n}\n\nsub get_chain_list \n{\n    my $header_fi\
le;\n    my @chain_list;\n    my @list;\n    my $n\
_chains;\n    my %chain_hasch;\n    my $pdb_file=@\
_[0];\n    my $c;\n    my %hasch;\n    my $chain;\\
n  \n    \n    $F=vfopen ($pdb_file, \"r\");\n    \
while ( <$F>)\n{\n\n\n	if (/SEQRES\\s+\\d+\\s+(\\S\
+)/)\n	  {\n	    $chain = substr($_,11,1);$chain=~\
s/\\s//g;if ( $chain eq \"\"){$chain=\" \";}\n	   \
 if (!$hasch{$chain}){$hasch{$chain}=1;push @chain\
_list, $chain;}\n	  }\n	if (/^ATOM/ || /^HETATM/)\\
n	  {\n	    $chain = substr($_,21,1); $chain=~s/\\\
s//g;if ( $chain eq \"\"){$chain=\" \";}\n	    if \
(!$hasch{$chain}){$hasch{$chain}=1;push @chain_lis\
t, $chain;}\n	  }\n      }\n\n\nclose ($F);\nif (!\
@chain_list)\n  {\n    @chain_list=(\"A\");\n  }\n\
\n\nreturn @chain_list;\n}\n\nsub token_is_in_list\
\n{\n\n    my @list=@_;\n    my $a;\n    \n    for\
 ($a=1; $a<=$#list; $a++)\n{\n	if ( $list[$a] eq $\
list[0]){return $a;}\n}\n}\n\nsub pdb_name2name_an\
d_chain \n{\n    my $pdb_file=@_[0];\n    my $pdb_\
file_in;\n    my @array;\n    my $chain;\n    my $\
c;\n\n    $pdb_file_in=$pdb_file;\n\n    $pdb_file\
=~/^(.{4})/;$pdb_id=$1;\n    @array=($pdb_file=~/(\
[\\w])/g);\n  \n  \n    $chain=uc ($array[4]);\n  \
  $chain=($chain eq \"\")?\"FIRST\":$chain;\n    \\
n    return ( $pdb_id, $chain);\n\n    if ( $#arra\
y==3){return ($pdb_id, \"FIRST\");}\n    elsif ( $\
#array<4){ return ($pdb_id, \"\");}\n    else {ret\
urn ( $pdb_id, $chain);}\n      \n    \n    \n}\ns\
ub get_main_hom_code \n{\n    my $pdb_file=@_[0];\\
n    my %hom, $n, $best, $best_h;\n    open (F, $p\
db_file);\n    while (<F>)\n{\n	if ( $_=~/^ATOM/)\\
n{\n	    $h=substr ($_,26, 1);\n	    $n=++$hom{$h}\
;\n	    if ($n>$best)\n{\n		$best=$n;\n		$best_h=$\
h;\n}\n}\n}\n    close (F);\n    return $best_h;\n\
}\n\n\nsub get_pdb_file \n{\n    my ($pdb_file_in)\
=(@_);\n    my $result;\n    my @letter;\n    my @\
chain;\n    my $v;\n    my $pdb_file=$pdb_file_in;\
\n\n    $pdb_file=($pdb_file_in=~/\\S+_S_(\\S+)/)?\
$1:$pdb_file_in;\n\n    if ($no_remote_pdb_dir==0)\
\n{\n	$no_remote_pdb_dir=1;\n	$result=get_pdb_file\
3 ($pdb_file);\n	$no_remote_pdb_dir=0;\n	if ( $res\
ult){return $result;}\n	else\n{\n	    \n	    lc ($\
pdb_file);\n	    $result=get_pdb_file3($pdb_file);\
\n	    return  $result;\n}\n}\n    else\n{\n	retur\
n get_pdb_file3 ($pdb_file);\n}\n    \n}\n\nsub ge\
t_pdb_file3 \n{\n    my $pdb_file_in=@_[0];\n    m\
y $result;\n    my @letter;\n    my @chain;\n    m\
y $lcfile;\n    my $ucfile;\n    my $pdb_file=$pdb\
_file_in;\n    \n    $lcfile=lc $pdb_file;\n    $u\
cfile=uc $pdb_file;\n\n    if ( ($result=get_pdb_f\
ile2 ($pdb_file))){return $result;}\n    \n\n    i\
f ($lcfile ne $pdb_file && ($result=get_pdb_file2 \
($lcfile))){return $result;}\n    if ($ucfile ne $\
pdb_file && ($result=get_pdb_file2 ($ucfile))){ret\
urn $result;}\n    \n   \n    \n    return \"\";\n\
}\nsub get_pdb_file2\n{\n    my $pdb_file=@_[0];\n\
    my $return_value;\n    \n    $return_value=\"\\
";\n    \n    if ( ($result=get_pdb_file1 ($pdb_fi\
le))){$return_value=$result;}\n    elsif ( !($pdb_\
file=~/\\.pdb/) && !($pdb_file=~/\\.PDB/))\n{\n	if\
 ( ($result=get_pdb_file1 (\"$pdb_file.pdb\"))){$r\
eturn_value=$result;}\n	elsif ( ($result=get_pdb_f\
ile1 (\"$pdb_file.PDB\"))){$return_value=$result;}\
\n\n	elsif ( ($result=get_pdb_file1 (\"pdb$pdb_fil\
e.pdb\"))){$return_value=$result;}	\n	elsif ( ($re\
sult=get_pdb_file1 (\"pdb$pdb_file.PDB\"))){$retur\
n_value=$result;}\n	elsif ( ($result=get_pdb_file1\
 (\"PDB$pdb_file.PDB\"))){$return_value=$result;}\\
n	elsif ( ($result=get_pdb_file1 (\"PDB$pdb_file.p\
db\"))){$return_value=$result;}\n	\n	\n	elsif ( ($\
result=get_pdb_file1 (\"$pdb_file.ent\"))){$return\
_value=$result;}\n	elsif ( ($result=get_pdb_file1 \
(\"pdb$pdb_file.ent\"))){$return_value=$result;}\n\
	elsif ( ($result=get_pdb_file1 (\"PDB$pdb_file.en\
t\"))){$return_value=$result;}\n\n	elsif ( ($resul\
t=get_pdb_file1 (\"$pdb_file.ENT\"))){$return_valu\
e=$result;}\n	elsif ( ($result=get_pdb_file1 (\"pd\
b$pdb_file.ENT\"))){$return_value=$result;}\n	elsi\
f ( ($result=get_pdb_file1 (\"PDB$pdb_file.ENT\"))\
){$return_value=$result;}\n	\n	\n	\n}\n    return \
$return_value;\n}\n    \nsub get_pdb_file1\n{\n   \
 my ($pdb_file)=(@_);\n    my $return_value;\n    \
\n\n    $return_value=\"\";\n    if ( ($result=get\
_pdb_file0 ($pdb_file))){$return_value=$result;}\n\
    elsif ( ($result=get_pdb_file0 (\"$pdb_file.Z\\
"))){$return_value=$result;}\n    elsif ( ($result\
=get_pdb_file0 (\"$pdb_file.gz\"))){$return_value=\
$result;}\n    elsif ( ($result=get_pdb_file0 (\"$\
pdb_file.GZ\"))){$return_value=$result;}\n    retu\
rn $return_value;\n}\nsub get_pdb_file0 \n{ \n    \
my ($pdb_file, $attempt)=(@_);\n    my $pdb_file=@\
_[0];\n    my $tmp_pdb_file;    \n    my $return_v\
alue;\n\n    if ( !$attempt){$attempt=1;}\n    \n \
   $local_pdb_file=\"$pdb_file\";\n    if ( $local\
_pdb_file eq \"\")\n{\n	$tmp_pdb_file=vtmpnam();\n\
	open F, \">$tmp_pdb_file\";\n	\n	while (<STDIN>){\
print F \"$_\";}\n	close (F);\n	\n	if (-e $tmp_pdb\
_file && &is_pdb_file ( $local_pdb_file))\n{return\
 $tmp_pdb_file;}\n}\n\n    $local_pdb_file=\"$pdb_\
file\";\n    &debug_print (\"\\nTry access local f\
ile: $local_pdb_file\");\n    \n    $local_pdb_fil\
e=&check_pdb_file4compression ($local_pdb_file);\n\
    if ( -e $local_pdb_file && (&is_pdb_file ($loc\
al_pdb_file) || $force_pdb))\n{\n	&debug_print ( \\
"\\n\\tIs in Current Dir\");\n	$tmp_pdb_file=vtmpn\
am();\n	`cp $local_pdb_file $tmp_pdb_file`;\n	retu\
rn $tmp_pdb_file;\n}\n    else\n{\n	&debug_print (\
\"\\n\\tFile Not in Current Dir\");\n}\n\n    if (\
$pdb_file=~/^pdb/||$pdb_file=~/^PDB/){$pdb_div=sub\
str ($pdb_file, 4, 2);}\n    else\n{\n	  $pdb_div=\
substr ($pdb_file, 1, 2);\n}\n    $local_pdb_file=\
\"$pdb_dir/$pdb_div/$pdb_file\";\n    $local_pdb_f\
ile=&check_pdb_file4compression ( $local_pdb_file)\
;\n    &debug_print (\"\\nTry access file From PDB\
_DIR: $local_pdb_file\");\n    if ($pdb_dir && -e \
$local_pdb_file && &is_pdb_file ($local_pdb_file))\
\n{\n	&debug_print ( \"\\n\\tIs in Local PDB DIR\"\
);\n	$tmp_pdb_file=vtmpnam();\n	`cp $local_pdb_fil\
e $tmp_pdb_file`;\n	return $tmp_pdb_file;\n}\n\n  \
  $local_pdb_file=\"$pdb_dir/$pdb_file\";\n    $lo\
cal_pdb_file=&check_pdb_file4compression ( $local_\
pdb_file);\n    &debug_print (\"\\nTry access file\
 From PDB_DIR: local_pdb_file\");\n    if ($pdb_di\
r && -e $local_pdb_file && &is_pdb_file ($local_pd\
b_file))\n{\n	&debug_print ( \"\\n\\tIs in Local P\
DB DIR\");\n	$tmp_pdb_file=vtmpnam();\n	`cp $local\
_pdb_file $tmp_pdb_file`;\n	return $tmp_pdb_file;\\
n}\n\n    $local_pdb_file=\"$pdb_dir$pdb_file\";\n\
    $local_pdb_file=&check_pdb_file4compression ( \
$local_pdb_file);\n    &debug_print (\"\\nTry acce\
ss file From PDB_DIR: $local_pdb_file\");\n    if \
($pdb_dir && -e $local_pdb_file && &is_pdb_file ($\
local_pdb_file))\n{\n	&debug_print ( \"\\n\\tIs in\
 Local PDB DIR\");\n	$tmp_pdb_file=vtmpnam();\n	`c\
p $local_pdb_file $tmp_pdb_file`;\n	return $tmp_pd\
b_file;\n}\n    else\n{&debug_print ( \"\\n\\tNot \
In Local Pdb Dir\");}\n\n    if ($cache ne \"NO\" \
&& $cache ne \"no\")\n{\n\n	$local_pdb_file=\"$cac\
he/$pdb_file\";\n	$local_pdb_file=&check_pdb_file4\
compression ( $local_pdb_file);\n	&debug_print(\"\\
\nTry access file From Cache: $local_pdb_file\");\\
n	if (-e $local_pdb_file && &is_pdb_file ($local_p\
db_file))\n{\n	    &debug_print ( \"\\n\\tIs in T-\
Coffee Cache\");\n	    $tmp_pdb_file=vtmpnam();\n	\
    `cp $local_pdb_file $tmp_pdb_file`;\n	    retu\
rn $tmp_pdb_file;\n}\n	else{&debug_print ( \"\\n\\\
tNot in Cache Dir\");}\n}\n\nif (!$no_remote_pdb_d\
ir) \n  {\n    my $value=remote_is_pdb_name ($pdb_\
file, $netaddress);\n    \n    my $return_value=\"\
\";\n    if ( &remote_is_pdb_name ($pdb_file, $net\
address)==1)\n      {\n	\n	&debug_print (\"\\n****\
**************************************************\
*\\nTry Remote Access for $pdb_file\");\n	$tmp_pdb\
_file=vtmpnam();\n	$netcommand=$netaddress;\n	$net\
command=~s/%%/$pdb_file/g;\n	&url2file(\"$netcomma\
nd\", \"$tmp_pdb_file.$netcompression\");\n	&debug\
_print(\"\\nREMOTE: $netcommand\\n\");\n	\n	$compr\
essed_tmp_file_name=\"$tmp_pdb_file.$netcompressio\
n\";\n	\n	if ($netcompression && -B $compressed_tm\
p_file_name)\n	  {\n	    my $r;\n	    &debug_print\
 (\"\\n\\tFile Found Remotely\");\n	    if (($r=sa\
fe_system ( \"$netcompression_pg $compressed_tmp_f\
ile_name\")!=$EXIT_SUCCESS) && $attempts<5)\n	    \
  {\n		&debug_print (\"\\n\\tProper Download Faile\
d Try again\");\n		unlink $compressed_tmp_file_nam\
e;\n		print \"\\nFailed to Download $compressed_tm\
p_file_name. New Attempt $attempt/5\\n\";\n		retur\
n &get_pdb_file0($pdb_file, $attempt+1);\n	      }\
\n	    elsif ($r== $EXIT_SUCCESS)\n	      {\n		&de\
bug_print (\"\\n\\tProper Download Succeeded \");\\
n		$return_value=$tmp_pdb_file;\n	      }\n	    el\
se\n	      {\n		&debug_print (\"\\n\\tProper Downl\
oad Failed \");\n		&debug_print (\"\\nFile Not Fou\
nd Remotely\");\n		unlink $compressed_tmp_file_nam\
e;\n	      }\n	  }\n	else\n	  {\n\n	    &debug_pri\
nt (\"\\nFile Not Found Remotely\");\n	    unlink \
$compressed_tmp_file_name;\n	  }\n	#Update cache i\
f required\n	if ($cache ne \"no\" && $cache ne \"u\
pdate\" && -e $return_value)\n	  {\n	    `cp $retu\
rn_value $cache/$pdb_file.pdb`;\n	    #`t_coffee -\
other_pg clean_cache.pl -file $pdb_file.pdb -dir $\
cache`;\n	  }\n      }\n    &debug_print (\"\\nRem\
ote Download Finished\");\n    return $return_valu\
e;\n  }\nreturn \"\";\n}\n\nsub check_pdb_file4com\
pression \n{\n    my $file=@_[0];\n    my $tmp;\n \
   my $r;\n    \n    $tmp=&vtmpnam();\n    if (-e \
$tmp){unlink $tmp;}\n    \n    $file=~s/\\/\\//\\/\
/g;\n    if    (-B $file && ($file=~/\\.Z/)) {`cp \
$file $tmp.Z`;`rm $tmp`;`gunzip $tmp.Z $SILENT`;$r\
=$tmp;}\n    elsif (-B $file && ($file=~/\\.gz/)){\
`cp $file $tmp.gz`;`gunzip $tmp.gz $SILENT`;return\
 $r=$tmp;}\n    elsif (-B $file ){`cp $file $tmp.g\
z`;`gunzip $tmp.gz $SILENT`;$r=$tmp;}\n    elsif (\
 -e $file ) {$r= $file;}\n    elsif ( -e \"$file.g\
z\" ){ `cp $file.gz $tmp.gz`;`gunzip     $tmp.gz $\
SILENT`;$r=$tmp;}    \n    elsif ( -e \"$file.Z\")\
 {`cp $file.Z  $tmp.Z`; `gunzip $tmp.Z $SILENT`;$r\
=$tmp;}\n    else  {$r= $file;}\n\n    if ( -e \"$\
tmp.Z\"){unlink \"$tmp.Z\";}\n    if ( -e \"$tmp.g\
z\"){unlink \"$tmp.gz\";}\n    \n    return $r;\n \
   \n}\n\n\n\n\n\n    \n\n\n\n\n\n\n\nsub vfopen \\
n{\n    my $file=@_[0];\n    my $mode=@_[1];\n    \
my $tmp;\n    my $F = new FileHandle;\n    \n    \\
n    $tmp=$file;\n	\n    \n    if ( $mode eq \"r\"\
 && !-e $file){ myexit(flush_error (\"Cannot open \
file $file\"));}\n    elsif ($mode eq \"w\"){$tmp=\
\">$file\";}\n    elsif ($mode eq \"a\"){$tmp=\">>\
$file\";}\n    \n    \n    open ($F,$tmp);\n    re\
turn $F;\n}\nsub debug_print\n{\n    my $message =\
@_[0];\n    if ($debug){print STDERR \"NO_REMOTE_P\
DB_DIR: $no_remote_pdb_dir - $message [DEBUG:extra\
ct_from_pdb]\";}\n    return;\n}\nsub is_aa \n{\n \
   my ($aa, $chain) =@_;\n\n    my $one;\n    my $\
trhee;\n    \n    if ( $onelett{$molecule_type{$ch\
ain}}->{$aa} eq 'X' || !$onelett{$molecule_type{$c\
hain}}->{$aa} ){return '';}\n    else\n      {\n	$\
one=$onelett{$molecule_type{$chain}}->{$aa};\n\n	$\
three=$threelett{$molecule_type{$chain}}->{$one};\\
n	\n\n	return $three;\n      }\n  }\n\n\n\n\n\nsub\
 url2file\n{\n    my ($address, $out, $wget_arg, $\
curl_arg)=(@_);\n    my ($pg, $flag, $r, $arg, $co\
unt);\n    \n    if (!$CONFIGURATION){&check_confi\
guration (\"wget\", \"INTERNET\", \"gzip\");$CONFI\
GURATION=1;}\n    \n    if (&pg_is_installed (\"wg\
et\"))   {$pg=\"wget\"; $flag=\"-O\";$arg=$wget_ar\
g;}\n    elsif (&pg_is_installed (\"curl\")){$pg=\\
"curl\"; $flag=\"-o\";$arg=$curl_arg;}\n    return\
 safe_system (\"$pg $flag$out $address >/dev/null \
2>/dev/null\");\n\n}\n\n\n\n\nsub pdbfile2chaintyp\
e\n  {\n    my $file=@_[0];\n    my %ct;\n    my $\
F;\n    \n    $F=vfopen ($file, \"r\");\n    while\
 (<$F>)\n      {\n	my $line=$_;\n	if ($line =~/^AT\
OM/)\n	  {\n	    my $C=substr($line,21,1);\n	    i\
f (!$ct{$C})\n	      {	\n		my $r=substr($line,17,3\
);\n		$r=~s/\\s+//;\n		if (length ($r)==1){$ct{$C}\
=\"R\";}\n		elsif (length ($r)==2){$ct{$C}=\"D\";}\
\n		elsif (length ($r)==3){$ct{$C}=\"P\";}\n		else\
 \n		  {\n		    myexit(flush_error(\"ERROR: Could \
not read RES_ID field in file $file\"));\n		  }\n	\
      }\n	  }\n      }\n    close ($F);\n    retur\
n %ct;\n  }\n   \n    \n\n\n\nsub fill_threelett_R\
NA\n{\n\n	my %threelett=(\n	'A', '  A',\n	'T', '  \
T',\n	'U', '  U',\n	'C', '  C',\n	'G', '  G',\n	'I\
', '  I', #Inosine\n	);\n	\n	return %threelett;\n\\
n}\n\n\nsub fill_onelett_RNA\n{\n	my   %onelett=(\\
n	'  A' => 'A',\n	'  T' => 'T',\n	'  U' => 'U',\n	\
'  C' => 'C',\n	'  G' => 'G',\n	'CSL' => 'X',\n	'U\
MS' => 'X',\n	'  I' => 'I',\n	'A' => 'A',\n	'T' =>\
 'T',\n	'U' => 'U',\n	'C' => 'C',\n	'G' => 'G',\n	\
'I' => 'I',\n	);\n\n	return %onelett;\n\n}\n\n\nsu\
b fill_onelett_DNA\n{\n	my   %onelett=(\n	' DA', '\
A',\n	' DT', 'T',\n	' DC', 'C',\n	' DG', 'G',\n	'D\
A', 'A',\n	'DT', 'T',\n	'DC', 'C',\n	'DG', 'G',\n	\
);\n\n	return %onelett;\n\n}\n\nsub fill_threelett\
_DNA\n{\n\n	my %threelett=(\n	'A', ' DA',\n	'T', '\
 DT',\n	'C', ' DC',\n	'G', ' DG',\n	);\n\n	return \
%threelett;\n\n}\n\n\n\n\nsub fill_threelett_prot\\
n{  \n  my %threelett;\n\n  %threelett=(\n'A', 'AL\
A',\n'C', 'CYS',\n'D', 'ASP',\n'E', 'GLU',\n'F', '\
PHE',\n'G', 'GLY',\n'H', 'HIS',\n'I', 'ILE',\n'K',\
 'LYS',\n'L', 'LEU',\n'N', 'ASN',\n'M', 'MET',\n'P\
', 'PRO',\n'Q', 'GLN',\n'R', 'ARG',\n'S', 'SER',\n\
'T', 'THR',\n'V', 'VAL',\n'W', 'TRP',\n'Y', 'TYR',\
\n);\n\nreturn %threelett;\n\n\n}\n\nsub fill_onel\
ett_prot\n{\n    my %onelett;\n    \n    %onelett=\
(\n\n'10A', 'X',\n'11O', 'X',\n'12A', 'X',\n'13P',\
 'X',\n'13R', 'X',\n'13S', 'X',\n'14W', 'X',\n'15P\
', 'X',\n'16A', 'X',\n'16G', 'X',\n'1AN', 'X',\n'1\
AP', 'X',\n'1AR', 'X',\n'1BH', 'X',\n'1BO', 'X',\n\
'1C5', 'X',\n'1CU', 'X',\n'1DA', 'X',\n'1GL', 'X',\
\n'1GN', 'X',\n'1IN', 'X',\n'1LU', 'L',\n'1MA', 'X\
',\n'1MC', 'X',\n'1MG', 'X',\n'1MZ', 'X',\n'1NA', \
'X',\n'1NB', 'X',\n'1NI', 'X',\n'1PA', 'A',\n'1PC'\
, 'X',\n'1PE', 'X',\n'1PG', 'X',\n'1PI', 'A',\n'1P\
M', 'X',\n'1PN', 'X',\n'1PU', 'X',\n'1PY', 'X',\n'\
1UN', 'X',\n'24T', 'X',\n'25T', 'X',\n'26P', 'X',\\
n'2AB', 'X',\n'2AM', 'X',\n'2AN', 'X',\n'2AP', 'X'\
,\n'2AR', 'X',\n'2AS', 'D',\n'2BL', 'X',\n'2BM', '\
X',\n'2CP', 'X',\n'2DA', 'X',\n'2DG', 'X',\n'2DP',\
 'X',\n'2DT', 'X',\n'2EP', 'X',\n'2EZ', 'X',\n'2FG\
', 'X',\n'2FL', 'X',\n'2FP', 'X',\n'2FU', 'X',\n'2\
GL', 'X',\n'2GP', 'X',\n'2HP', 'X',\n'2IB', 'X',\n\
'2IP', 'X',\n'2LU', 'L',\n'2MA', 'X',\n'2MD', 'X',\
\n'2ME', 'X',\n'2MG', 'X',\n'2ML', 'L',\n'2MO', 'X\
',\n'2MR', 'R',\n'2MU', 'X',\n'2MZ', 'X',\n'2NO', \
'X',\n'2NP', 'X',\n'2OG', 'X',\n'2PA', 'X',\n'2PC'\
, 'X',\n'2PE', 'X',\n'2PG', 'X',\n'2PH', 'X',\n'2P\
I', 'X',\n'2PL', 'X',\n'2PP', 'X',\n'2PU', 'X',\n'\
2SI', 'X',\n'2TB', 'X',\n'34C', 'X',\n'35G', 'X',\\
n'3AA', 'X',\n'3AD', 'X',\n'3AH', 'H',\n'3AN', 'X'\
,\n'3AP', 'X',\n'3AT', 'X',\n'3BT', 'X',\n'3CH', '\
X',\n'3CN', 'X',\n'3CO', 'X',\n'3CP', 'X',\n'3DR',\
 'X',\n'3EP', 'X',\n'3FM', 'X',\n'3GA', 'X',\n'3GP\
', 'X',\n'3HB', 'X',\n'3HC', 'X',\n'3HP', 'X',\n'3\
IB', 'X',\n'3ID', 'X',\n'3IN', 'X',\n'3MA', 'X',\n\
'3MB', 'X',\n'3MC', 'X',\n'3MD', 'D',\n'3MF', 'X',\
\n'3MP', 'X',\n'3MT', 'X',\n'3OL', 'X',\n'3PA', 'X\
',\n'3PG', 'X',\n'3PO', 'X',\n'3PP', 'X',\n'3PY', \
'X',\n'49A', 'X',\n'4AB', 'X',\n'4AM', 'X',\n'4AN'\
, 'X',\n'4AP', 'X',\n'4BA', 'X',\n'4BT', 'X',\n'4C\
A', 'X',\n'4CO', 'X',\n'4HP', 'X',\n'4IP', 'X',\n'\
4MO', 'X',\n'4MV', 'X',\n'4MZ', 'X',\n'4NC', 'X',\\
n'4NP', 'X',\n'4OX', 'X',\n'4PB', 'X',\n'4PN', 'X'\
,\n'4PP', 'X',\n'4SC', 'X',\n'4SU', 'X',\n'4TB', '\
X',\n'55C', 'X',\n'5AD', 'X',\n'5AN', 'X',\n'5AT',\
 'X',\n'5CM', 'X',\n'5GP', 'X',\n'5HP', 'E',\n'5HT\
', 'X',\n'5IT', 'X',\n'5IU', 'X',\n'5MB', 'X',\n'5\
MC', 'X',\n'5MD', 'X',\n'5MP', 'X',\n'5MU', 'X',\n\
'5NC', 'X',\n'5OB', 'X',\n'5PA', 'X',\n'5PV', 'X',\
\n'6AB', 'X',\n'6CT', 'X',\n'6HA', 'X',\n'6HC', 'X\
',\n'6HG', 'X',\n'6HT', 'X',\n'6IN', 'X',\n'6MO', \
'X',\n'6MP', 'X',\n'6PG', 'X',\n'6WO', 'X',\n'70U'\
, 'X',\n'7DG', 'X',\n'7HP', 'X',\n'7I2', 'X',\n'7M\
G', 'X',\n'7MQ', 'X',\n'7NI', 'X',\n'87Y', 'X',\n'\
8AD', 'X',\n'8BR', 'X',\n'8IG', 'X',\n'8IN', 'X',\\
n'8OG', 'X',\n'95A', 'X',\n'9AD', 'X',\n'9AM', 'X'\
,\n'9AP', 'X',\n'9DG', 'X',\n'9DI', 'X',\n'9HX', '\
X',\n'9OH', 'X',\n'9TA', 'X',\n'A12', 'X',\n'A15',\
 'X',\n'A23', 'X',\n'A24', 'X',\n'A26', 'X',\n'A2G\
', 'X',\n'A2P', 'X',\n'A32', 'X',\n'A3P', 'X',\n'A\
4P', 'X',\n'A5P', 'X',\n'A70', 'X',\n'A76', 'X',\n\
'A77', 'X',\n'A78', 'X',\n'A79', 'X',\n'A80', 'X',\
\n'A85', 'X',\n'A88', 'X',\n'A9A', 'X',\n'AA3', 'X\
',\n'AA4', 'X',\n'AA6', 'X',\n'AAA', 'X',\n'AAB', \
'X',\n'AAC', 'X',\n'AAE', 'X',\n'AAG', 'R',\n'AAH'\
, 'X',\n'AAM', 'X',\n'AAN', 'X',\n'AAP', 'X',\n'AA\
R', 'R',\n'AAS', 'X',\n'AAT', 'X',\n'ABA', 'X',\n'\
ABC', 'X',\n'ABD', 'X',\n'ABE', 'X',\n'ABH', 'X',\\
n'ABI', 'X',\n'ABK', 'X',\n'ABM', 'X',\n'ABN', 'X'\
,\n'ABP', 'X',\n'ABR', 'X',\n'ABS', 'X',\n'ABU', '\
X',\n'AC1', 'X',\n'AC2', 'X',\n'ACA', 'X',\n'ACB',\
 'D',\n'ACC', 'C',\n'ACD', 'X',\n'ACE', 'X',\n'ACH\
', 'X',\n'ACI', 'X',\n'ACL', 'R',\n'ACM', 'X',\n'A\
CN', 'X',\n'ACO', 'X',\n'ACP', 'X',\n'ACQ', 'X',\n\
'ACR', 'X',\n'ACS', 'X',\n'ACT', 'X',\n'ACV', 'V',\
\n'ACX', 'X',\n'ACY', 'X',\n'AD2', 'X',\n'AD3', 'X\
',\n'ADC', 'X',\n'ADD', 'X',\n'ADE', 'X',\n'ADH', \
'X',\n'ADI', 'X',\n'ADM', 'X',\n'ADN', 'X',\n'ADP'\
, 'X',\n'ADQ', 'X',\n'ADR', 'X',\n'ADS', 'X',\n'AD\
T', 'X',\n'ADU', 'X',\n'ADW', 'X',\n'ADX', 'X',\n'\
AE2', 'X',\n'AEA', 'X',\n'AEB', 'X',\n'AEI', 'D',\\
n'AEN', 'X',\n'AET', 'T',\n'AF1', 'X',\n'AF3', 'X'\
,\n'AFA', 'D',\n'AFP', 'X',\n'AG7', 'X',\n'AGB', '\
X',\n'AGF', 'X',\n'AGL', 'X',\n'AGM', 'R',\n'AGN',\
 'X',\n'AGP', 'X',\n'AGS', 'X',\n'AGU', 'X',\n'AH0\
', 'X',\n'AH1', 'X',\n'AHA', 'X',\n'AHB', 'D',\n'A\
HC', 'X',\n'AHF', 'X',\n'AHG', 'X',\n'AHH', 'X',\n\
'AHM', 'X',\n'AHO', 'X',\n'AHP', 'X',\n'AHS', 'X',\
\n'AHT', 'Y',\n'AHU', 'X',\n'AHX', 'X',\n'AI1', 'X\
',\n'AI2', 'X',\n'AIB', 'X',\n'AIC', 'X',\n'AIM', \
'X',\n'AIP', 'X',\n'AIQ', 'X',\n'AIR', 'X',\n'AJ3'\
, 'X',\n'AKB', 'X',\n'AKG', 'X',\n'AKR', 'X',\n'AL\
1', 'X',\n'AL2', 'X',\n'AL3', 'X',\n'AL4', 'X',\n'\
AL5', 'X',\n'AL6', 'X',\n'AL7', 'X',\n'AL8', 'X',\\
n'AL9', 'X',\n'ALA', 'A',\n'ALB', 'X',\n'ALC', 'X'\
,\n'ALD', 'L',\n'ALE', 'X',\n'ALF', 'X',\n'ALG', '\
X',\n'ALL', 'X',\n'ALM', 'A',\n'ALN', 'A',\n'ALO',\
 'T',\n'ALP', 'X',\n'ALQ', 'X',\n'ALR', 'X',\n'ALS\
', 'X',\n'ALT', 'A',\n'ALY', 'K',\n'ALZ', 'X',\n'A\
MA', 'X',\n'AMB', 'X',\n'AMC', 'X',\n'AMD', 'X',\n\
'AMG', 'X',\n'AMH', 'X',\n'AMI', 'X',\n'AML', 'X',\
\n'AMN', 'X',\n'AMO', 'X',\n'AMP', 'X',\n'AMQ', 'X\
',\n'AMR', 'X',\n'AMS', 'X',\n'AMT', 'X',\n'AMU', \
'X',\n'AMW', 'X',\n'AMX', 'X',\n'AMY', 'X',\n'ANA'\
, 'X',\n'ANB', 'X',\n'ANC', 'X',\n'AND', 'X',\n'AN\
E', 'X',\n'ANI', 'X',\n'ANL', 'X',\n'ANO', 'X',\n'\
ANP', 'X',\n'ANS', 'X',\n'ANT', 'X',\n'AOE', 'X',\\
n'AOP', 'X',\n'AP1', 'X',\n'AP2', 'X',\n'AP3', 'X'\
,\n'AP4', 'X',\n'AP5', 'X',\n'AP6', 'X',\n'APA', '\
X',\n'APB', 'X',\n'APC', 'X',\n'APE', 'F',\n'APF',\
 'X',\n'APG', 'X',\n'APH', 'A',\n'API', 'X',\n'APL\
', 'X',\n'APM', 'X',\n'APN', 'G',\n'APP', 'X',\n'A\
PQ', 'X',\n'APR', 'X',\n'APS', 'X',\n'APT', 'X',\n\
'APU', 'X',\n'APX', 'X',\n'APY', 'X',\n'APZ', 'X',\
\n'AQS', 'X',\n'AR1', 'X',\n'AR2', 'X',\n'ARA', 'X\
',\n'ARB', 'X',\n'ARC', 'X',\n'ARD', 'X',\n'ARG', \
'R',\n'ARH', 'X',\n'ARI', 'X',\n'ARM', 'R',\n'ARN'\
, 'X',\n'ARO', 'R',\n'ARP', 'X',\n'ARQ', 'X',\n'AR\
S', 'X',\n'AS1', 'R',\n'AS2', 'X',\n'ASA', 'D',\n'\
ASB', 'D',\n'ASC', 'X',\n'ASD', 'X',\n'ASE', 'X',\\
n'ASF', 'X',\n'ASI', 'X',\n'ASK', 'D',\n'ASL', 'X'\
,\n'ASM', 'N',\n'ASO', 'X',\n'ASP', 'D',\n'ASQ', '\
X',\n'ASU', 'X',\n'ATA', 'X',\n'ATC', 'X',\n'ATD',\
 'X',\n'ATF', 'X',\n'ATG', 'X',\n'ATH', 'X',\n'ATM\
', 'X',\n'ATO', 'X',\n'ATP', 'X',\n'ATQ', 'X',\n'A\
TR', 'X',\n'ATT', 'X',\n'ATY', 'X',\n'ATZ', 'X',\n\
'AUC', 'X',\n'AUR', 'X',\n'AVG', 'X',\n'AXP', 'X',\
\n'AYA', 'A',\n'AZ2', 'X',\n'AZA', 'X',\n'AZC', 'X\
',\n'AZD', 'X',\n'AZE', 'X',\n'AZI', 'X',\n'AZL', \
'X',\n'AZM', 'X',\n'AZR', 'X',\n'AZT', 'X',\n'B12'\
, 'X',\n'B1F', 'F',\n'B2A', 'A',\n'B2F', 'F',\n'B2\
I', 'I',\n'B2V', 'V',\n'B3I', 'X',\n'B3P', 'X',\n'\
B7G', 'X',\n'B96', 'X',\n'B9A', 'X',\n'BA1', 'X',\\
n'BAA', 'X',\n'BAB', 'X',\n'BAC', 'X',\n'BAF', 'X'\
,\n'BAH', 'X',\n'BAI', 'X',\n'BAK', 'X',\n'BAL', '\
A',\n'BAM', 'X',\n'BAO', 'X',\n'BAP', 'X',\n'BAR',\
 'X',\n'BAS', 'X',\n'BAT', 'F',\n'BAY', 'X',\n'BAZ\
', 'X',\n'BB1', 'X',\n'BB2', 'X',\n'BBA', 'X',\n'B\
BH', 'X',\n'BBS', 'X',\n'BBT', 'X',\n'BBZ', 'X',\n\
'BCA', 'X',\n'BCB', 'X',\n'BCC', 'X',\n'BCD', 'X',\
\n'BCL', 'X',\n'BCN', 'X',\n'BCR', 'X',\n'BCS', 'C\
',\n'BCT', 'X',\n'BCY', 'X',\n'BCZ', 'X',\n'BDA', \
'X',\n'BDG', 'X',\n'BDK', 'X',\n'BDM', 'X',\n'BDN'\
, 'X',\n'BDS', 'X',\n'BE1', 'X',\n'BE2', 'X',\n'BE\
A', 'X',\n'BEF', 'X',\n'BEN', 'X',\n'BEO', 'X',\n'\
BEP', 'X',\n'BER', 'X',\n'BES', 'X',\n'BET', 'X',\\
n'BEZ', 'X',\n'BF2', 'X',\n'BFA', 'X',\n'BFD', 'X'\
,\n'BFP', 'X',\n'BFS', 'X',\n'BFU', 'X',\n'BG6', '\
X',\n'BGF', 'X',\n'BGG', 'X',\n'BGL', 'X',\n'BGN',\
 'X',\n'BGP', 'X',\n'BGX', 'X',\n'BH4', 'X',\n'BHA\
', 'X',\n'BHC', 'X',\n'BHD', 'D',\n'BHO', 'X',\n'B\
HS', 'X',\n'BIC', 'X',\n'BIN', 'X',\n'BIO', 'X',\n\
'BIP', 'X',\n'BIS', 'X',\n'BIZ', 'X',\n'BJH', 'X',\
\n'BJI', 'X',\n'BJP', 'X',\n'BLA', 'X',\n'BLB', 'X\
',\n'BLE', 'L',\n'BLG', 'P',\n'BLI', 'X',\n'BLM', \
'X',\n'BLV', 'X',\n'BLY', 'K',\n'BM1', 'X',\n'BM2'\
, 'X',\n'BM5', 'X',\n'BM9', 'X',\n'BMA', 'X',\n'BM\
D', 'X',\n'BME', 'X',\n'BMP', 'X',\n'BMQ', 'X',\n'\
BMS', 'X',\n'BMT', 'T',\n'BMU', 'X',\n'BMY', 'X',\\
n'BMZ', 'X',\n'BNA', 'X',\n'BNG', 'X',\n'BNI', 'X'\
,\n'BNN', 'F',\n'BNO', 'L',\n'BNS', 'X',\n'BNZ', '\
X',\n'BO3', 'X',\n'BO4', 'X',\n'BOC', 'X',\n'BOG',\
 'X',\n'BOM', 'X',\n'BOT', 'X',\n'BOX', 'X',\n'BOZ\
', 'X',\n'BPA', 'X',\n'BPB', 'X',\n'BPD', 'X',\n'B\
PG', 'X',\n'BPH', 'X',\n'BPI', 'X',\n'BPJ', 'X',\n\
'BPM', 'X',\n'BPN', 'X',\n'BPO', 'X',\n'BPP', 'X',\
\n'BPT', 'X',\n'BPY', 'X',\n'BRB', 'X',\n'BRC', 'X\
',\n'BRE', 'X',\n'BRI', 'X',\n'BRL', 'X',\n'BRM', \
'X',\n'BRN', 'X',\n'BRO', 'X',\n'BRS', 'X',\n'BRU'\
, 'X',\n'BRZ', 'X',\n'BSB', 'X',\n'BSI', 'X',\n'BS\
P', 'X',\n'BT1', 'X',\n'BT2', 'X',\n'BT3', 'X',\n'\
BTA', 'L',\n'BTB', 'X',\n'BTC', 'C',\n'BTD', 'X',\\
n'BTN', 'X',\n'BTP', 'X',\n'BTR', 'W',\n'BU1', 'X'\
,\n'BUA', 'X',\n'BUB', 'X',\n'BUC', 'X',\n'BUG', '\
X',\n'BUL', 'X',\n'BUM', 'X',\n'BUQ', 'X',\n'BUT',\
 'X',\n'BVD', 'X',\n'BX3', 'X',\n'BYS', 'X',\n'BZ1\
', 'X',\n'BZA', 'X',\n'BZB', 'X',\n'BZC', 'X',\n'B\
ZD', 'X',\n'BZF', 'X',\n'BZI', 'X',\n'BZM', 'X',\n\
'BZO', 'X',\n'BZP', 'X',\n'BZQ', 'X',\n'BZS', 'X',\
\n'BZT', 'X',\n'C02', 'X',\n'C11', 'X',\n'C1O', 'X\
',\n'C20', 'X',\n'C24', 'X',\n'C2F', 'X',\n'C2O', \
'X',\n'C2P', 'X',\n'C3M', 'X',\n'C3P', 'X',\n'C3X'\
, 'X',\n'C48', 'X',\n'C4M', 'X',\n'C4X', 'X',\n'C5\
C', 'X',\n'C5M', 'X',\n'C5P', 'X',\n'C5X', 'X',\n'\
C60', 'X',\n'C6C', 'X',\n'C6M', 'X',\n'C78', 'X',\\
n'C8E', 'X',\n'CA3', 'X',\n'CA5', 'X',\n'CAA', 'X'\
,\n'CAB', 'X',\n'CAC', 'X',\n'CAD', 'X',\n'CAF', '\
C',\n'CAG', 'X',\n'CAH', 'X',\n'CAL', 'X',\n'CAM',\
 'X',\n'CAN', 'X',\n'CAO', 'X',\n'CAP', 'X',\n'CAQ\
', 'X',\n'CAR', 'X',\n'CAS', 'C',\n'CAT', 'X',\n'C\
AV', 'X',\n'CAY', 'C',\n'CAZ', 'X',\n'CB3', 'X',\n\
'CB4', 'X',\n'CBA', 'X',\n'CBD', 'X',\n'CBG', 'X',\
\n'CBI', 'X',\n'CBL', 'X',\n'CBM', 'X',\n'CBN', 'X\
',\n'CBO', 'X',\n'CBP', 'X',\n'CBS', 'X',\n'CBX', \
'X',\n'CBZ', 'X',\n'CC0', 'X',\n'CC1', 'X',\n'CCC'\
, 'X',\n'CCH', 'X',\n'CCI', 'X',\n'CCM', 'X',\n'CC\
N', 'X',\n'CCO', 'X',\n'CCP', 'X',\n'CCR', 'X',\n'\
CCS', 'C',\n'CCV', 'X',\n'CCY', 'X',\n'CD1', 'X',\\
n'CDC', 'X',\n'CDE', 'X',\n'CDF', 'X',\n'CDI', 'X'\
,\n'CDL', 'X',\n'CDM', 'X',\n'CDP', 'X',\n'CDR', '\
X',\n'CDU', 'X',\n'CE1', 'X',\n'CEA', 'C',\n'CEB',\
 'X',\n'CEC', 'X',\n'CED', 'X',\n'CEF', 'X',\n'CEH\
', 'X',\n'CEM', 'X',\n'CEO', 'X',\n'CEP', 'X',\n'C\
EQ', 'X',\n'CER', 'X',\n'CES', 'G',\n'CET', 'X',\n\
'CFC', 'X',\n'CFF', 'X',\n'CFM', 'X',\n'CFO', 'X',\
\n'CFP', 'X',\n'CFS', 'X',\n'CFX', 'X',\n'CGN', 'X\
',\n'CGP', 'X',\n'CGS', 'X',\n'CGU', 'E',\n'CH2', \
'X',\n'CH3', 'X',\n'CHA', 'X',\n'CHB', 'X',\n'CHD'\
, 'X',\n'CHF', 'X',\n'CHG', 'G',\n'CHI', 'X',\n'CH\
N', 'X',\n'CHO', 'X',\n'CHP', 'G',\n'CHR', 'X',\n'\
CHS', 'F',\n'CHT', 'X',\n'CHX', 'X',\n'CIC', 'X',\\
n'CIN', 'X',\n'CIP', 'X',\n'CIR', 'X',\n'CIT', 'X'\
,\n'CIU', 'X',\n'CKI', 'X',\n'CL1', 'X',\n'CL2', '\
X',\n'CLA', 'X',\n'CLB', 'A',\n'CLC', 'S',\n'CLD',\
 'A',\n'CLE', 'L',\n'CLF', 'X',\n'CLK', 'S',\n'CLL\
', 'X',\n'CLM', 'X',\n'CLN', 'X',\n'CLO', 'X',\n'C\
LP', 'X',\n'CLQ', 'X',\n'CLR', 'X',\n'CLS', 'X',\n\
'CLT', 'X',\n'CLX', 'X',\n'CLY', 'X',\n'CMA', 'R',\
\n'CMC', 'X',\n'CMD', 'X',\n'CME', 'C',\n'CMG', 'X\
',\n'CMK', 'X',\n'CMN', 'X',\n'CMO', 'X',\n'CMP', \
'X',\n'CMR', 'X',\n'CMS', 'X',\n'CMT', 'C',\n'CMX'\
, 'X',\n'CNA', 'X',\n'CNC', 'X',\n'CND', 'X',\n'CN\
H', 'X',\n'CNM', 'X',\n'CNN', 'X',\n'CNO', 'X',\n'\
CNP', 'X',\n'CO2', 'X',\n'CO3', 'X',\n'CO5', 'X',\\
n'CO8', 'X',\n'COA', 'X',\n'COB', 'X',\n'COC', 'X'\
,\n'COD', 'X',\n'COE', 'X',\n'COF', 'X',\n'COH', '\
X',\n'COI', 'X',\n'COJ', 'X',\n'COL', 'X',\n'COM',\
 'X',\n'CON', 'X',\n'COP', 'X',\n'COR', 'X',\n'COS\
', 'X',\n'COT', 'X',\n'COY', 'X',\n'CP1', 'G',\n'C\
P2', 'X',\n'CP4', 'X',\n'CPA', 'X',\n'CPB', 'X',\n\
'CPC', 'X',\n'CPD', 'X',\n'CPG', 'X',\n'CPH', 'X',\
\n'CPI', 'X',\n'CPM', 'X',\n'CPN', 'G',\n'CPO', 'X\
',\n'CPP', 'X',\n'CPQ', 'X',\n'CPR', 'X',\n'CPS', \
'X',\n'CPT', 'X',\n'CPU', 'X',\n'CPV', 'X',\n'CPY'\
, 'X',\n'CR1', 'X',\n'CR6', 'X',\n'CRA', 'X',\n'CR\
B', 'X',\n'CRC', 'X',\n'CRG', 'X',\n'CRH', 'X',\n'\
CRO', 'T',\n'CRP', 'X',\n'CRQ', 'X',\n'CRS', 'X',\\
n'CRT', 'X',\n'CRY', 'X',\n'CSA', 'C',\n'CSB', 'X'\
,\n'CSD', 'C',\n'CSE', 'C',\n'CSH', 'X',\n'CSI', '\
X',\n'CSN', 'X',\n'CSO', 'C',\n'CSP', 'C',\n'CSR',\
 'C',\n'CSS', 'C',\n'CST', 'X',\n'CSW', 'C',\n'CSX\
', 'C',\n'CSY', 'X',\n'CSZ', 'C',\n'CT3', 'X',\n'C\
TA', 'X',\n'CTB', 'X',\n'CTC', 'X',\n'CTD', 'X',\n\
'CTH', 'T',\n'CTO', 'X',\n'CTP', 'X',\n'CTR', 'X',\
\n'CTS', 'X',\n'CTT', 'X',\n'CTY', 'X',\n'CTZ', 'X\
',\n'CU1', 'X',\n'CUA', 'X',\n'CUC', 'X',\n'CUL', \
'X',\n'CUO', 'X',\n'CUZ', 'X',\n'CVI', 'X',\n'CXF'\
, 'X',\n'CXL', 'X',\n'CXM', 'M',\n'CXN', 'X',\n'CX\
P', 'X',\n'CXS', 'X',\n'CY1', 'C',\n'CY3', 'X',\n'\
CYB', 'X',\n'CYC', 'X',\n'CYF', 'C',\n'CYG', 'C',\\
n'CYH', 'X',\n'CYL', 'X',\n'CYM', 'C',\n'CYN', 'X'\
,\n'CYO', 'X',\n'CYP', 'X',\n'CYQ', 'C',\n'CYS', '\
C',\n'CYU', 'X',\n'CYY', 'X',\n'CYZ', 'X',\n'CZH',\
 'X',\n'CZZ', 'C',\n'D12', 'X',\n'D13', 'X',\n'D16\
', 'X',\n'D18', 'X',\n'D19', 'X',\n'D1P', 'X',\n'D\
24', 'X',\n'D34', 'X',\n'D35', 'X',\n'D4D', 'X',\n\
'D4T', 'X',\n'D6G', 'X',\n'DA2', 'R',\n'DA3', 'X',\
\n'DA6', 'X',\n'DA7', 'X',\n'DAA', 'X',\n'DAB', 'X\
',\n'DAC', 'X',\n'DAD', 'X',\n'DAE', 'X',\n'DAF', \
'X',\n'DAG', 'X',\n'DAH', 'A',\n'DAJ', 'X',\n'DAK'\
, 'X',\n'DAL', 'A',\n'DAM', 'A',\n'DAN', 'X',\n'DA\
O', 'X',\n'DAP', 'X',\n'DAQ', 'X',\n'DAR', 'R',\n'\
DAS', 'D',\n'DAT', 'X',\n'DAU', 'X',\n'DAV', 'X',\\
n'DBA', 'X',\n'DBD', 'X',\n'DBF', 'X',\n'DBG', 'X'\
,\n'DBI', 'X',\n'DBV', 'X',\n'DBY', 'Y',\n'DCA', '\
X',\n'DCB', 'X',\n'DCE', 'X',\n'DCF', 'X',\n'DCG',\
 'X',\n'DCH', 'X',\n'DCI', 'I',\n'DCL', 'X',\n'DCM\
', 'X',\n'DCP', 'X',\n'DCS', 'X',\n'DCT', 'X',\n'D\
CY', 'C',\n'DCZ', 'X',\n'DDA', 'X',\n'DDB', 'X',\n\
'DDC', 'X',\n'DDF', 'X',\n'DDG', 'X',\n'DDH', 'X',\
\n'DDL', 'X',\n'DDM', 'X',\n'DDO', 'L',\n'DDP', 'X\
',\n'DDQ', 'X',\n'DDT', 'Y',\n'DDU', 'X',\n'DEA', \
'X',\n'DEB', 'X',\n'DEC', 'X',\n'DEF', 'X',\n'DEL'\
, 'X',\n'DEM', 'X',\n'DEN', 'X',\n'DEP', 'X',\n'DE\
Q', 'X',\n'DES', 'X',\n'DET', 'X',\n'DFC', 'X',\n'\
DFG', 'X',\n'DFI', 'X',\n'DFL', 'X',\n'DFO', 'X',\\
n'DFP', 'X',\n'DFR', 'X',\n'DFT', 'X',\n'DFV', 'X'\
,\n'DFX', 'X',\n'DG2', 'X',\n'DG3', 'X',\n'DG6', '\
X',\n'DGA', 'X',\n'DGD', 'X',\n'DGG', 'X',\n'DGL',\
 'E',\n'DGN', 'Q',\n'DGP', 'X',\n'DGT', 'X',\n'DGX\
', 'X',\n'DH2', 'X',\n'DHA', 'A',\n'DHB', 'X',\n'D\
HC', 'X',\n'DHD', 'X',\n'DHE', 'X',\n'DHF', 'X',\n\
'DHG', 'X',\n'DHI', 'H',\n'DHL', 'X',\n'DHM', 'X',\
\n'DHN', 'V',\n'DHP', 'X',\n'DHQ', 'X',\n'DHR', 'X\
',\n'DHS', 'X',\n'DHT', 'X',\n'DHU', 'X',\n'DHY', \
'X',\n'DHZ', 'X',\n'DI2', 'X',\n'DI3', 'G',\n'DI4'\
, 'X',\n'DI5', 'X',\n'DIA', 'X',\n'DIC', 'X',\n'DI\
F', 'X',\n'DIG', 'X',\n'DII', 'X',\n'DIL', 'I',\n'\
DIM', 'X',\n'DIO', 'X',\n'DIP', 'X',\n'DIQ', 'X',\\
n'DIS', 'X',\n'DIT', 'X',\n'DIV', 'V',\n'DIX', 'X'\
,\n'DIY', 'X',\n'DKA', 'X',\n'DLA', 'X',\n'DLE', '\
L',\n'DLF', 'X',\n'DLS', 'K',\n'DLY', 'K',\n'DM1',\
 'X',\n'DM2', 'X',\n'DM3', 'X',\n'DM4', 'X',\n'DM5\
', 'X',\n'DM6', 'X',\n'DM7', 'X',\n'DM8', 'X',\n'D\
M9', 'X',\n'DMA', 'X',\n'DMB', 'X',\n'DMC', 'X',\n\
'DMD', 'X',\n'DME', 'X',\n'DMF', 'X',\n'DMG', 'G',\
\n'DMH', 'N',\n'DMI', 'X',\n'DMJ', 'X',\n'DML', 'X\
',\n'DMM', 'X',\n'DMN', 'X',\n'DMO', 'X',\n'DMP', \
'X',\n'DMQ', 'X',\n'DMR', 'X',\n'DMS', 'X',\n'DMT'\
, 'X',\n'DMV', 'X',\n'DMY', 'X',\n'DNC', 'X',\n'DN\
D', 'X',\n'DNH', 'X',\n'DNJ', 'X',\n'DNN', 'X',\n'\
DNP', 'X',\n'DNQ', 'X',\n'DNR', 'X',\n'DO2', 'X',\\
n'DO3', 'X',\n'DOA', 'X',\n'DOB', 'X',\n'DOC', 'X'\
,\n'DOH', 'D',\n'DOM', 'X',\n'DOS', 'X',\n'DOX', '\
X',\n'DP5', 'X',\n'DP7', 'X',\n'DPA', 'X',\n'DPC',\
 'X',\n'DPD', 'X',\n'DPE', 'X',\n'DPG', 'X',\n'DPH\
', 'F',\n'DPM', 'X',\n'DPN', 'F',\n'DPO', 'X',\n'D\
PP', 'X',\n'DPR', 'P',\n'DPS', 'X',\n'DPT', 'X',\n\
'DPX', 'X',\n'DPY', 'X',\n'DPZ', 'X',\n'DQH', 'X',\
\n'DQN', 'X',\n'DR1', 'X',\n'DRB', 'X',\n'DRC', 'X\
',\n'DRI', 'X',\n'DRP', 'X',\n'DRT', 'X',\n'DRU', \
'X',\n'DSA', 'X',\n'DSB', 'X',\n'DSC', 'X',\n'DSD'\
, 'X',\n'DSE', 'S',\n'DSI', 'X',\n'DSN', 'S',\n'DS\
P', 'D',\n'DSR', 'X',\n'DSS', 'X',\n'DSX', 'X',\n'\
DSY', 'X',\n'DTB', 'X',\n'DTD', 'X',\n'DTH', 'T',\\
n'DTN', 'X',\n'DTO', 'X',\n'DTP', 'X',\n'DTQ', 'X'\
,\n'DTR', 'W',\n'DTT', 'X',\n'DTY', 'Y',\n'DUD', '\
X',\n'DUO', 'X',\n'DUR', 'X',\n'DUT', 'X',\n'DVA',\
 'V',\n'DVR', 'X',\n'DX9', 'X',\n'DXA', 'X',\n'DXB\
', 'X',\n'DXC', 'X',\n'DXG', 'X',\n'DXX', 'X',\n'D\
ZF', 'X',\n'E09', 'X',\n'E20', 'X',\n'E2P', 'X',\n\
'E3G', 'X',\n'E4N', 'X',\n'E4P', 'X',\n'E64', 'X',\
\n'E6C', 'X',\n'E96', 'X',\n'E97', 'X',\n'EA2', 'X\
',\n'EAA', 'X',\n'EAP', 'X',\n'EBP', 'X',\n'EBW', \
'X',\n'ECO', 'X',\n'EDA', 'X',\n'EDC', 'X',\n'EDE'\
, 'X',\n'EDO', 'X',\n'EDR', 'X',\n'EEB', 'X',\n'EE\
E', 'X',\n'EFC', 'X',\n'EFZ', 'X',\n'EG1', 'X',\n'\
EG2', 'X',\n'EG3', 'X',\n'EGC', 'X',\n'EGL', 'X',\\
n'EHP', 'A',\n'EIC', 'X',\n'EJT', 'X',\n'ELA', 'X'\
,\n'EMB', 'X',\n'EMC', 'X',\n'EMD', 'X',\n'EMM', '\
X',\n'EMO', 'X',\n'EMP', 'X',\n'EMR', 'X',\n'ENA',\
 'X',\n'ENC', 'X',\n'ENH', 'X',\n'ENO', 'X',\n'ENP\
', 'X',\n'EOA', 'X',\n'EOH', 'X',\n'EOT', 'X',\n'E\
OX', 'X',\n'EPA', 'X',\n'EPE', 'X',\n'EPH', 'X',\n\
'EPI', 'X',\n'EPN', 'X',\n'EPO', 'X',\n'EPT', 'X',\
\n'EPU', 'X',\n'EPX', 'X',\n'EPY', 'X',\n'EQI', 'X\
',\n'EQP', 'X',\n'EQU', 'X',\n'ERG', 'X',\n'ERI', \
'X',\n'ERY', 'X',\n'ESC', 'X',\n'ESD', 'X',\n'ESI'\
, 'X',\n'ESO', 'X',\n'ESP', 'X',\n'EST', 'X',\n'ES\
X', 'X',\n'ETA', 'X',\n'ETC', 'X',\n'ETD', 'X',\n'\
ETF', 'X',\n'ETH', 'X',\n'ETI', 'X',\n'ETN', 'X',\\
n'ETO', 'X',\n'ETP', 'X',\n'ETR', 'X',\n'ETS', 'X'\
,\n'ETY', 'X',\n'EU3', 'X',\n'EUG', 'X',\n'EYS', '\
C',\n'F09', 'X',\n'F2B', 'X',\n'F3S', 'X',\n'F42',\
 'X',\n'F43', 'X',\n'F4S', 'X',\n'F6B', 'X',\n'F6P\
', 'X',\n'F89', 'X',\n'FA1', 'X',\n'FA5', 'F',\n'F\
AA', 'X',\n'FAB', 'X',\n'FAC', 'X',\n'FAD', 'X',\n\
'FAF', 'X',\n'FAG', 'X',\n'FAM', 'X',\n'FAR', 'X',\
\n'FAS', 'X',\n'FAT', 'X',\n'FBA', 'X',\n'FBE', 'X\
',\n'FBI', 'X',\n'FBP', 'X',\n'FBQ', 'X',\n'FBS', \
'X',\n'FBT', 'X',\n'FBU', 'X',\n'FCA', 'X',\n'FCB'\
, 'X',\n'FCI', 'X',\n'FCN', 'X',\n'FCO', 'X',\n'FC\
R', 'X',\n'FCT', 'X',\n'FCX', 'X',\n'FCY', 'C',\n'\
FD1', 'F',\n'FD2', 'F',\n'FD3', 'F',\n'FD4', 'F',\\
n'FDA', 'X',\n'FDC', 'X',\n'FDI', 'X',\n'FDP', 'X'\
,\n'FDS', 'X',\n'FE2', 'X',\n'FEA', 'X',\n'FEL', '\
X',\n'FEM', 'X',\n'FEN', 'X',\n'FEO', 'X',\n'FEP',\
 'X',\n'FER', 'X',\n'FES', 'X',\n'FFB', 'X',\n'FFC\
', 'X',\n'FFF', 'X',\n'FFO', 'X',\n'FGL', 'G',\n'F\
HB', 'X',\n'FHC', 'X',\n'FHP', 'X',\n'FHU', 'X',\n\
'FID', 'X',\n'FII', 'X',\n'FIP', 'X',\n'FK5', 'X',\
\n'FKA', 'X',\n'FKI', 'X',\n'FKP', 'X',\n'FL2', 'X\
',\n'FL9', 'X',\n'FLA', 'A',\n'FLC', 'X',\n'FLD', \
'X',\n'FLE', 'L',\n'FLF', 'X',\n'FLO', 'X',\n'FLP'\
, 'X',\n'FLT', 'Y',\n'FLU', 'X',\n'FLX', 'X',\n'FM\
1', 'X',\n'FM2', 'X',\n'FMA', 'X',\n'FMB', 'X',\n'\
FMC', 'X',\n'FME', 'M',\n'FMN', 'X',\n'FMP', 'X',\\
n'FMR', 'X',\n'FMS', 'X',\n'FMT', 'X',\n'FNE', 'X'\
,\n'FNP', 'X',\n'FNS', 'X',\n'FOC', 'X',\n'FOE', '\
X',\n'FOG', 'F',\n'FOH', 'X',\n'FOK', 'X',\n'FOL',\
 'X',\n'FON', 'X',\n'FOP', 'X',\n'FOR', 'X',\n'FOS\
', 'X',\n'FPA', 'X',\n'FPC', 'X',\n'FPI', 'X',\n'F\
PO', 'X',\n'FPP', 'X',\n'FPT', 'X',\n'FQP', 'X',\n\
'FRA', 'X',\n'FRD', 'F',\n'FRU', 'X',\n'FS3', 'X',\
\n'FS4', 'X',\n'FSB', 'X',\n'FSO', 'X',\n'FSX', 'X\
',\n'FTC', 'X',\n'FTP', 'X',\n'FTR', 'W',\n'FTT', \
'X',\n'FTY', 'Y',\n'FUA', 'X',\n'FUC', 'X',\n'FUM'\
, 'X',\n'FUP', 'X',\n'FVF', 'X',\n'FXP', 'X',\n'FX\
V', 'X',\n'FYA', 'F',\n'G16', 'X',\n'G1P', 'X',\n'\
G20', 'X',\n'G21', 'X',\n'G23', 'X',\n'G26', 'X',\\
n'G28', 'X',\n'G2F', 'X',\n'G37', 'X',\n'G39', 'X'\
,\n'G3H', 'X',\n'G3P', 'X',\n'G4D', 'X',\n'G6D', '\
X',\n'G6P', 'X',\n'G6Q', 'X',\n'G7M', 'X',\n'GA2',\
 'X',\n'GAA', 'X',\n'GAB', 'X',\n'GAC', 'X',\n'GAI\
', 'X',\n'GAL', 'X',\n'GAM', 'X',\n'GAN', 'X',\n'G\
AO', 'X',\n'GAP', 'X',\n'GAR', 'G',\n'GAS', 'X',\n\
'GAT', 'X',\n'GBC', 'X',\n'GBI', 'X',\n'GBP', 'X',\
\n'GBS', 'X',\n'GBX', 'X',\n'GC4', 'X',\n'GCA', 'X\
',\n'GCD', 'X',\n'GCG', 'G',\n'GCH', 'G',\n'GCK', \
'X',\n'GCL', 'X',\n'GCM', 'X',\n'GCN', 'X',\n'GCO'\
, 'X',\n'GCP', 'X',\n'GCR', 'X',\n'GCS', 'X',\n'GC\
U', 'X',\n'GD3', 'X',\n'GDB', 'X',\n'GDM', 'X',\n'\
GDN', 'X',\n'GDP', 'X',\n'GDS', 'X',\n'GDU', 'X',\\
n'GE1', 'X',\n'GE2', 'X',\n'GE3', 'X',\n'GEA', 'X'\
,\n'GEL', 'X',\n'GEM', 'X',\n'GEN', 'X',\n'GEP', '\
X',\n'GER', 'X',\n'GFP', 'X',\n'GGB', 'X',\n'GGL',\
 'E',\n'GGP', 'X',\n'GHP', 'G',\n'GIP', 'X',\n'GIS\
', 'X',\n'GKR', 'X',\n'GL2', 'X',\n'GL3', 'G',\n'G\
L4', 'X',\n'GL5', 'X',\n'GL7', 'X',\n'GL9', 'X',\n\
'GLA', 'X',\n'GLB', 'X',\n'GLC', 'X',\n'GLD', 'X',\
\n'GLE', 'X',\n'GLF', 'X',\n'GLG', 'X',\n'GLH', 'Q\
',\n'GLI', 'X',\n'GLL', 'X',\n'GLM', 'G',\n'GLN', \
'Q',\n'GLO', 'X',\n'GLP', 'X',\n'GLR', 'X',\n'GLS'\
, 'X',\n'GLT', 'X',\n'GLU', 'E',\n'GLV', 'X',\n'GL\
W', 'X',\n'GLY', 'G',\n'GLZ', 'X',\n'GM1', 'X',\n'\
GMA', 'X',\n'GMC', 'X',\n'GMH', 'X',\n'GMP', 'X',\\
n'GMY', 'X',\n'GN7', 'X',\n'GNA', 'X',\n'GNB', 'X'\
,\n'GNH', 'X',\n'GNP', 'X',\n'GNT', 'X',\n'GOA', '\
X',\n'GOL', 'X',\n'GOX', 'X',\n'GP1', 'X',\n'GP3',\
 'X',\n'GP4', 'X',\n'GP6', 'X',\n'GP8', 'X',\n'GPB\
', 'E',\n'GPC', 'X',\n'GPE', 'X',\n'GPG', 'X',\n'G\
PI', 'X',\n'GPJ', 'X',\n'GPL', 'K',\n'GPM', 'X',\n\
'GPN', 'G',\n'GPP', 'X',\n'GPR', 'X',\n'GPS', 'X',\
\n'GPX', 'X',\n'GR1', 'X',\n'GR3', 'X',\n'GR4', 'X\
',\n'GSA', 'X',\n'GSB', 'X',\n'GSC', 'G',\n'GSE', \
'S',\n'GSH', 'X',\n'GSP', 'X',\n'GSR', 'X',\n'GSS'\
, 'X',\n'GT9', 'C',\n'GTA', 'X',\n'GTB', 'X',\n'GT\
D', 'X',\n'GTE', 'X',\n'GTH', 'T',\n'GTN', 'X',\n'\
GTO', 'X',\n'GTP', 'X',\n'GTR', 'X',\n'GTS', 'X',\\
n'GTT', 'X',\n'GTX', 'X',\n'GTZ', 'X',\n'GU7', 'X'\
,\n'GUA', 'X',\n'GUD', 'X',\n'GUM', 'X',\n'GUN', '\
X',\n'GUP', 'X',\n'GUR', 'X',\n'GW3', 'X',\n'GZZ',\
 'X',\n'H2B', 'X',\n'H2P', 'H',\n'H2S', 'X',\n'H2U\
', 'X',\n'H4B', 'X',\n'H5M', 'P',\n'H5P', 'X',\n'H\
AA', 'X',\n'HAB', 'X',\n'HAC', 'A',\n'HAD', 'X',\n\
'HAE', 'X',\n'HAG', 'X',\n'HAI', 'X',\n'HAM', 'X',\
\n'HAP', 'X',\n'HAQ', 'X',\n'HAR', 'R',\n'HAS', 'X\
',\n'HAV', 'V',\n'HAX', 'X',\n'HAZ', 'X',\n'HBA', \
'X',\n'HBC', 'X',\n'HBD', 'X',\n'HBI', 'X',\n'HBO'\
, 'X',\n'HBU', 'X',\n'HBY', 'X',\n'HC0', 'X',\n'HC\
1', 'X',\n'HC4', 'X',\n'HCA', 'X',\n'HCC', 'X',\n'\
HCI', 'X',\n'HCS', 'X',\n'HDA', 'X',\n'HDD', 'X',\\
n'HDF', 'X',\n'HDN', 'X',\n'HDS', 'X',\n'HDZ', 'X'\
,\n'HE1', 'X',\n'HE6', 'X',\n'HEA', 'X',\n'HEB', '\
X',\n'HEC', 'X',\n'HED', 'X',\n'HEE', 'X',\n'HEF',\
 'X',\n'HEG', 'X',\n'HEM', 'X',\n'HEN', 'X',\n'HEO\
', 'X',\n'HEP', 'X',\n'HEU', 'X',\n'HEV', 'X',\n'H\
EX', 'X',\n'HEZ', 'X',\n'HF1', 'X',\n'HFA', 'X',\n\
'HFP', 'X',\n'HGA', 'Q',\n'HGB', 'X',\n'HGC', 'X',\
\n'HGI', 'X',\n'HGU', 'X',\n'HHO', 'X',\n'HHP', 'X\
',\n'HIB', 'X',\n'HIC', 'H',\n'HII', 'X',\n'HIN', \
'X',\n'HIO', 'X',\n'HIP', 'H',\n'HIS', 'H',\n'HLE'\
, 'X',\n'HLT', 'X',\n'HMA', 'A',\n'HMB', 'X',\n'HM\
C', 'X',\n'HMD', 'X',\n'HMF', 'A',\n'HMG', 'X',\n'\
HMH', 'X',\n'HMI', 'L',\n'HMM', 'X',\n'HMN', 'X',\\
n'HMO', 'X',\n'HMP', 'X',\n'HMR', 'R',\n'HNI', 'X'\
,\n'HNP', 'X',\n'HOA', 'X',\n'HOE', 'X',\n'HOH', '\
X',\n'HOM', 'X',\n'HOP', 'X',\n'HOQ', 'X',\n'HP1',\
 'A',\n'HP2', 'A',\n'HP3', 'X',\n'HPA', 'X',\n'HPB\
', 'X',\n'HPC', 'X',\n'HPD', 'X',\n'HPE', 'A',\n'H\
PG', 'X',\n'HPH', 'F',\n'HPP', 'X',\n'HPQ', 'F',\n\
'HPR', 'X',\n'HPT', 'X',\n'HPY', 'X',\n'HQO', 'X',\
\n'HQQ', 'X',\n'HQU', 'X',\n'HRG', 'R',\n'HRI', 'X\
',\n'HSA', 'X',\n'HSE', 'S',\n'HSF', 'X',\n'HSM', \
'X',\n'HSO', 'H',\n'HSP', 'X',\n'HT1', 'X',\n'HT2'\
, 'X',\n'HTA', 'X',\n'HTL', 'X',\n'HTO', 'X',\n'HT\
P', 'X',\n'HTR', 'W',\n'HUP', 'X',\n'HUX', 'X',\n'\
HV5', 'A',\n'HV7', 'X',\n'HV8', 'X',\n'HXA', 'X',\\
n'HXC', 'X',\n'HXP', 'X',\n'HY1', 'X',\n'HYA', 'X'\
,\n'HYB', 'X',\n'HYD', 'X',\n'HYG', 'X',\n'HYP', '\
P',\n'I06', 'X',\n'I10', 'X',\n'I11', 'X',\n'I17',\
 'X',\n'I2P', 'X',\n'I3N', 'X',\n'I3P', 'X',\n'I40\
', 'X',\n'I48', 'X',\n'I4B', 'X',\n'I52', 'X',\n'I\
5P', 'X',\n'I84', 'G',\n'IAG', 'G',\n'IAS', 'X',\n\
'IB2', 'X',\n'IBB', 'X',\n'IBP', 'X',\n'IBR', 'X',\
\n'IBS', 'X',\n'IBZ', 'X',\n'IC1', 'X',\n'ICA', 'X\
',\n'ICI', 'X',\n'ICL', 'X',\n'ICP', 'X',\n'ICT', \
'X',\n'ICU', 'X',\n'ID2', 'X',\n'IDC', 'X',\n'IDG'\
, 'X',\n'IDH', 'X',\n'IDM', 'X',\n'IDO', 'X',\n'ID\
P', 'X',\n'IDR', 'X',\n'IDS', 'X',\n'IDT', 'X',\n'\
IDU', 'X',\n'IFG', 'X',\n'IFP', 'X',\n'IGL', 'X',\\
n'IGN', 'X',\n'IGP', 'X',\n'IGU', 'X',\n'IH1', 'X'\
,\n'IH2', 'X',\n'IH3', 'X',\n'IHB', 'X',\n'IHN', '\
X',\n'IHP', 'X',\n'IIC', 'X',\n'IIL', 'I',\n'IIP',\
 'X',\n'IK2', 'X',\n'IKT', 'X',\n'ILA', 'I',\n'ILE\
', 'I',\n'ILG', 'X',\n'ILO', 'X',\n'ILX', 'I',\n'I\
M1', 'X',\n'IM2', 'X',\n'IMC', 'X',\n'IMD', 'X',\n\
'IME', 'X',\n'IMF', 'X',\n'IMG', 'X',\n'IMH', 'X',\
\n'IMI', 'X',\n'IML', 'I',\n'IMM', 'X',\n'IMN', 'X\
',\n'IMO', 'X',\n'IMP', 'X',\n'IMR', 'X',\n'IMU', \
'X',\n'IN0', 'D',\n'IN1', 'R',\n'IN2', 'K',\n'IN3'\
, 'L',\n'IN4', 'X',\n'IN5', 'A',\n'IN6', 'L',\n'IN\
7', 'X',\n'IN8', 'X',\n'IN9', 'X',\n'INA', 'L',\n'\
INB', 'X',\n'INC', 'X',\n'IND', 'X',\n'INE', 'X',\\
n'INF', 'F',\n'ING', 'F',\n'INH', 'R',\n'INI', 'X'\
,\n'INJ', 'X',\n'INK', 'X',\n'INL', 'X',\n'INM', '\
X',\n'INN', 'A',\n'INO', 'X',\n'INP', 'X',\n'INQ',\
 'X',\n'INR', 'X',\n'INS', 'X',\n'INT', 'V',\n'INU\
', 'X',\n'INV', 'X',\n'INW', 'X',\n'INX', 'X',\n'I\
NY', 'X',\n'INZ', 'X',\n'IOA', 'X',\n'IOB', 'X',\n\
'IOC', 'X',\n'IOD', 'X',\n'IOE', 'X',\n'IOF', 'X',\
\n'IOH', 'X',\n'IOL', 'X',\n'IOP', 'X',\n'IP1', 'X\
',\n'IP2', 'X',\n'IP3', 'X',\n'IP4', 'X',\n'IPA', \
'X',\n'IPB', 'X',\n'IPD', 'X',\n'IPG', 'G',\n'IPH'\
, 'X',\n'IPL', 'X',\n'IPM', 'X',\n'IPN', 'X',\n'IP\
O', 'F',\n'IPP', 'X',\n'IPS', 'X',\n'IPT', 'X',\n'\
IPU', 'X',\n'IPY', 'A',\n'IQB', 'X',\n'IQP', 'X',\\
n'IQS', 'X',\n'IR3', 'X',\n'IRI', 'X',\n'IRP', 'X'\
,\n'ISA', 'X',\n'ISF', 'X',\n'ISO', 'X',\n'ISP', '\
X',\n'ISQ', 'X',\n'ISU', 'X',\n'ITM', 'X',\n'ITP',\
 'X',\n'ITR', 'W',\n'ITS', 'X',\n'ITU', 'X',\n'IU5\
', 'X',\n'IUM', 'X',\n'IUR', 'X',\n'IVA', 'X',\n'I\
YG', 'G',\n'IYR', 'Y',\n'J77', 'X',\n'J78', 'X',\n\
'J80', 'X',\n'JE2', 'X',\n'JEN', 'X',\n'JST', 'X',\
\n'K21', 'X',\n'KAH', 'X',\n'KAI', 'X',\n'KAM', 'X\
',\n'KAN', 'X',\n'KAP', 'X',\n'KCP', 'X',\n'KCX', \
'K',\n'KDO', 'X',\n'KEF', 'X',\n'KET', 'X',\n'KGR'\
, 'X',\n'KH1', 'X',\n'KIF', 'X',\n'KIV', 'V',\n'KN\
I', 'X',\n'KPH', 'K',\n'KTH', 'X',\n'KTN', 'X',\n'\
KTP', 'X',\n'KWT', 'X',\n'L04', 'X',\n'L1P', 'X',\\
n'L24', 'E',\n'L2P', 'X',\n'L34', 'E',\n'L37', 'E'\
,\n'L3P', 'X',\n'L4P', 'X',\n'L75', 'X',\n'LAC', '\
X',\n'LAD', 'X',\n'LAK', 'X',\n'LAM', 'X',\n'LAR',\
 'X',\n'LAT', 'X',\n'LAX', 'X',\n'LCO', 'X',\n'LCP\
', 'X',\n'LCS', 'X',\n'LDA', 'X',\n'LDO', 'L',\n'L\
DP', 'X',\n'LEA', 'X',\n'LEO', 'X',\n'LEU', 'L',\n\
'LG2', 'X',\n'LG6', 'X',\n'LGC', 'X',\n'LGP', 'X',\
\n'LHG', 'X',\n'LHY', 'F',\n'LI1', 'X',\n'LIG', 'X\
',\n'LIL', 'X',\n'LIM', 'X',\n'LIN', 'X',\n'LIO', \
'X',\n'LIP', 'X',\n'LLA', 'X',\n'LLP', 'K',\n'LLY'\
, 'K',\n'LMG', 'X',\n'LML', 'X',\n'LMT', 'X',\n'LM\
U', 'X',\n'LMZ', 'X',\n'LNK', 'X',\n'LNL', 'X',\n'\
LNO', 'X',\n'LOF', 'X',\n'LOL', 'L',\n'LOM', 'X',\\
n'LOR', 'X',\n'LOS', 'X',\n'LOV', 'L',\n'LOX', 'X'\
,\n'LP1', 'X',\n'LP2', 'R',\n'LPA', 'X',\n'LPC', '\
X',\n'LPF', 'X',\n'LPL', 'X',\n'LPM', 'X',\n'LPP',\
 'X',\n'LRB', 'X',\n'LRU', 'X',\n'LS1', 'X',\n'LS2\
', 'X',\n'LS3', 'X',\n'LS4', 'X',\n'LS5', 'X',\n'L\
TA', 'X',\n'LTL', 'X',\n'LTR', 'W',\n'LUM', 'X',\n\
'LVS', 'L',\n'LXC', 'X',\n'LY2', 'X',\n'LY3', 'X',\
\n'LYA', 'X',\n'LYB', 'X',\n'LYC', 'X',\n'LYD', 'X\
',\n'LYM', 'K',\n'LYN', 'X',\n'LYS', 'K',\n'LYT', \
'X',\n'LYW', 'X',\n'LYZ', 'K',\n'M1A', 'X',\n'M1G'\
, 'X',\n'M2G', 'X',\n'M3L', 'K',\n'M6P', 'X',\n'M6\
T', 'X',\n'M7G', 'X',\n'MA1', 'X',\n'MA2', 'X',\n'\
MA3', 'X',\n'MA4', 'X',\n'MA6', 'X',\n'MAA', 'A',\\
n'MAB', 'X',\n'MAC', 'X',\n'MAE', 'X',\n'MAG', 'X'\
,\n'MAH', 'X',\n'MAI', 'R',\n'MAK', 'X',\n'MAL', '\
X',\n'MAM', 'X',\n'MAN', 'X',\n'MAO', 'X',\n'MAP',\
 'X',\n'MAR', 'X',\n'MAS', 'X',\n'MAT', 'X',\n'MAU\
', 'X',\n'MAZ', 'X',\n'MBA', 'X',\n'MBD', 'X',\n'M\
BG', 'X',\n'MBH', 'X',\n'MBN', 'X',\n'MBO', 'X',\n\
'MBR', 'X',\n'MBS', 'X',\n'MBV', 'X',\n'MBZ', 'X',\
\n'MCA', 'X',\n'MCD', 'X',\n'MCE', 'X',\n'MCG', 'G\
',\n'MCI', 'X',\n'MCN', 'X',\n'MCP', 'X',\n'MCT', \
'X',\n'MCY', 'X',\n'MD2', 'X',\n'MDA', 'X',\n'MDC'\
, 'X',\n'MDG', 'X',\n'MDH', 'X',\n'MDL', 'X',\n'MD\
M', 'X',\n'MDN', 'X',\n'MDP', 'X',\n'ME6', 'X',\n'\
MEB', 'X',\n'MEC', 'X',\n'MEL', 'X',\n'MEN', 'N',\\
n'MEP', 'X',\n'MER', 'X',\n'MES', 'X',\n'MET', 'M'\
,\n'MEV', 'X',\n'MF2', 'X',\n'MF3', 'M',\n'MFB', '\
X',\n'MFD', 'X',\n'MFU', 'X',\n'MG7', 'X',\n'MGA',\
 'X',\n'MGB', 'X',\n'MGD', 'X',\n'MGG', 'R',\n'MGL\
', 'X',\n'MGN', 'Q',\n'MGO', 'X',\n'MGP', 'X',\n'M\
GR', 'X',\n'MGS', 'X',\n'MGT', 'X',\n'MGU', 'X',\n\
'MGY', 'G',\n'MHB', 'X',\n'MHF', 'X',\n'MHL', 'L',\
\n'MHM', 'X',\n'MHO', 'M',\n'MHS', 'H',\n'MHZ', 'X\
',\n'MIA', 'X',\n'MIC', 'X',\n'MID', 'X',\n'MIL', \
'X',\n'MIM', 'X',\n'MIN', 'G',\n'MIP', 'X',\n'MIS'\
, 'S',\n'MIT', 'X',\n'MJI', 'X',\n'MK1', 'X',\n'MK\
C', 'X',\n'MLA', 'X',\n'MLC', 'X',\n'MLE', 'L',\n'\
MLN', 'X',\n'MLT', 'X',\n'MLY', 'K',\n'MLZ', 'K',\\
n'MM3', 'X',\n'MM4', 'X',\n'MMA', 'X',\n'MMC', 'X'\
,\n'MME', 'M',\n'MMO', 'R',\n'MMP', 'X',\n'MMQ', '\
X',\n'MMT', 'X',\n'MN1', 'X',\n'MN2', 'X',\n'MN3',\
 'X',\n'MN5', 'X',\n'MN7', 'X',\n'MN8', 'X',\n'MNA\
', 'X',\n'MNB', 'X',\n'MNC', 'X',\n'MNG', 'X',\n'M\
NL', 'L',\n'MNO', 'X',\n'MNP', 'X',\n'MNQ', 'X',\n\
'MNS', 'X',\n'MNT', 'X',\n'MNV', 'V',\n'MO1', 'X',\
\n'MO2', 'X',\n'MO3', 'X',\n'MO4', 'X',\n'MO5', 'X\
',\n'MO6', 'X',\n'MOA', 'X',\n'MOB', 'X',\n'MOC', \
'X',\n'MOE', 'X',\n'MOG', 'X',\n'MOH', 'X',\n'MOL'\
, 'X',\n'MOO', 'X',\n'MOP', 'X',\n'MOR', 'X',\n'MO\
S', 'X',\n'MOT', 'X',\n'MOX', 'X',\n'MP1', 'X',\n'\
MP3', 'X',\n'MPA', 'X',\n'MPB', 'X',\n'MPC', 'X',\\
n'MPD', 'X',\n'MPG', 'X',\n'MPH', 'M',\n'MPI', 'X'\
,\n'MPJ', 'M',\n'MPL', 'X',\n'MPN', 'X',\n'MPO', '\
X',\n'MPP', 'X',\n'MPQ', 'G',\n'MPR', 'X',\n'MPS',\
 'X',\n'MQ0', 'X',\n'MQ7', 'X',\n'MQ8', 'X',\n'MQ9\
', 'X',\n'MQI', 'X',\n'MR2', 'X',\n'MRC', 'X',\n'M\
RM', 'X',\n'MRP', 'X',\n'MS2', 'X',\n'MSA', 'X',\n\
'MSB', 'X',\n'MSD', 'X',\n'MSE', 'M',\n'MSF', 'X',\
\n'MSI', 'X',\n'MSO', 'M',\n'MSQ', 'X',\n'MST', 'X\
',\n'MSU', 'X',\n'MTA', 'X',\n'MTB', 'X',\n'MTC', \
'X',\n'MTD', 'X',\n'MTE', 'X',\n'MTF', 'X',\n'MTG'\
, 'X',\n'MTO', 'X',\n'MTS', 'X',\n'MTT', 'X',\n'MT\
X', 'X',\n'MTY', 'Y',\n'MUG', 'X',\n'MUP', 'X',\n'\
MUR', 'X',\n'MVA', 'V',\n'MW1', 'X',\n'MW2', 'X',\\
n'MXA', 'X',\n'MXY', 'X',\n'MYA', 'X',\n'MYC', 'X'\
,\n'MYG', 'X',\n'MYR', 'X',\n'MYS', 'X',\n'MYT', '\
X',\n'MZM', 'X',\n'N1T', 'X',\n'N25', 'X',\n'N2B',\
 'X',\n'N3T', 'X',\n'N4B', 'X',\n'NA2', 'X',\n'NA5\
', 'X',\n'NA6', 'X',\n'NAA', 'X',\n'NAB', 'X',\n'N\
AC', 'X',\n'NAD', 'X',\n'NAE', 'X',\n'NAF', 'X',\n\
'NAG', 'X',\n'NAH', 'X',\n'NAI', 'X',\n'NAL', 'A',\
\n'NAM', 'A',\n'NAN', 'X',\n'NAO', 'X',\n'NAP', 'X\
',\n'NAQ', 'X',\n'NAR', 'X',\n'NAS', 'X',\n'NAU', \
'X',\n'NAV', 'X',\n'NAW', 'X',\n'NAX', 'X',\n'NAY'\
, 'X',\n'NBA', 'X',\n'NBD', 'X',\n'NBE', 'X',\n'NB\
G', 'X',\n'NBN', 'X',\n'NBP', 'X',\n'NBS', 'X',\n'\
NBU', 'X',\n'NCA', 'X',\n'NCB', 'A',\n'NCD', 'X',\\
n'NCH', 'X',\n'NCM', 'X',\n'NCN', 'X',\n'NCO', 'X'\
,\n'NCR', 'X',\n'NCS', 'X',\n'ND4', 'X',\n'NDA', '\
X',\n'NDC', 'X',\n'NDD', 'X',\n'NDO', 'X',\n'NDP',\
 'X',\n'NDT', 'X',\n'NEA', 'X',\n'NEB', 'X',\n'NED\
', 'X',\n'NEM', 'H',\n'NEN', 'X',\n'NEO', 'X',\n'N\
EP', 'H',\n'NEQ', 'X',\n'NES', 'X',\n'NET', 'X',\n\
'NEV', 'X',\n'NFA', 'F',\n'NFE', 'X',\n'NFG', 'X',\
\n'NFP', 'X',\n'NFS', 'X',\n'NG6', 'X',\n'NGA', 'X\
',\n'NGL', 'X',\n'NGM', 'X',\n'NGO', 'X',\n'NGP', \
'X',\n'NGT', 'X',\n'NGU', 'X',\n'NH2', 'X',\n'NH3'\
, 'X',\n'NH4', 'X',\n'NHD', 'X',\n'NHE', 'X',\n'NH\
M', 'X',\n'NHP', 'X',\n'NHR', 'X',\n'NHS', 'X',\n'\
NI1', 'X',\n'NI2', 'X',\n'NIC', 'X',\n'NID', 'X',\\
n'NIK', 'X',\n'NIO', 'X',\n'NIP', 'X',\n'NIT', 'X'\
,\n'NIU', 'X',\n'NIY', 'Y',\n'NLA', 'X',\n'NLE', '\
L',\n'NLG', 'X',\n'NLN', 'L',\n'NLP', 'L',\n'NM1',\
 'X',\n'NMA', 'A',\n'NMB', 'X',\n'NMC', 'G',\n'NMD\
', 'X',\n'NME', 'X',\n'NMN', 'X',\n'NMO', 'X',\n'N\
MQ', 'X',\n'NMX', 'X',\n'NMY', 'X',\n'NNH', 'R',\n\
'NNO', 'X',\n'NO2', 'X',\n'NO3', 'X',\n'NOA', 'X',\
\n'NOD', 'X',\n'NOJ', 'X',\n'NON', 'X',\n'NOP', 'X\
',\n'NOR', 'X',\n'NOS', 'X',\n'NOV', 'X',\n'NOX', \
'X',\n'NP3', 'X',\n'NPA', 'X',\n'NPC', 'X',\n'NPD'\
, 'X',\n'NPE', 'X',\n'NPF', 'X',\n'NPH', 'C',\n'NP\
I', 'X',\n'NPL', 'X',\n'NPN', 'X',\n'NPO', 'X',\n'\
NPP', 'X',\n'NPT', 'X',\n'NPY', 'X',\n'NRG', 'R',\\
n'NRI', 'X',\n'NS1', 'X',\n'NS5', 'X',\n'NSP', 'X'\
,\n'NTA', 'X',\n'NTB', 'X',\n'NTC', 'X',\n'NTH', '\
X',\n'NTM', 'X',\n'NTP', 'X',\n'NTS', 'X',\n'NTU',\
 'X',\n'NTZ', 'X',\n'NU1', 'X',\n'NVA', 'V',\n'NVI\
', 'X',\n'NVP', 'X',\n'NW1', 'X',\n'NYP', 'X',\n'O\
4M', 'X',\n'OAA', 'X',\n'OAI', 'X',\n'OAP', 'X',\n\
'OAR', 'X',\n'OAS', 'S',\n'OBA', 'X',\n'OBN', 'X',\
\n'OC1', 'X',\n'OC2', 'X',\n'OC3', 'X',\n'OC4', 'X\
',\n'OC5', 'X',\n'OC6', 'X',\n'OC7', 'X',\n'OCL', \
'X',\n'OCM', 'X',\n'OCN', 'X',\n'OCO', 'X',\n'OCP'\
, 'X',\n'OCS', 'C',\n'OCT', 'X',\n'OCV', 'K',\n'OC\
Y', 'C',\n'ODA', 'X',\n'ODS', 'X',\n'OES', 'X',\n'\
OET', 'X',\n'OF1', 'X',\n'OF2', 'X',\n'OF3', 'X',\\
n'OFL', 'X',\n'OFO', 'X',\n'OHE', 'X',\n'OHO', 'X'\
,\n'OHT', 'X',\n'OIC', 'X',\n'OIP', 'X',\n'OKA', '\
X',\n'OLA', 'X',\n'OLE', 'X',\n'OLI', 'X',\n'OLO',\
 'X',\n'OMB', 'X',\n'OMC', 'X',\n'OMD', 'X',\n'OME\
', 'X',\n'OMG', 'X',\n'OMP', 'X',\n'OMT', 'M',\n'O\
MU', 'X',\n'ONE', 'X',\n'ONL', 'L',\n'ONP', 'X',\n\
'OPA', 'X',\n'OPD', 'X',\n'OPE', 'X',\n'OPG', 'X',\
\n'OPH', 'X',\n'OPN', 'X',\n'OPP', 'X',\n'OPR', 'R\
',\n'ORN', 'X',\n'ORO', 'X',\n'ORP', 'X',\n'OSB', \
'X',\n'OSS', 'X',\n'OTA', 'X',\n'OTB', 'X',\n'OTE'\
, 'X',\n'OTG', 'X',\n'OUT', 'X',\n'OVA', 'X',\n'OW\
Q', 'X',\n'OXA', 'X',\n'OXE', 'X',\n'OXI', 'X',\n'\
OXL', 'X',\n'OXM', 'X',\n'OXN', 'X',\n'OXO', 'X',\\
n'OXP', 'X',\n'OXS', 'X',\n'OXY', 'X',\n'P11', 'A'\
,\n'P24', 'X',\n'P28', 'X',\n'P2P', 'X',\n'P2U', '\
X',\n'P3M', 'X',\n'P4C', 'X',\n'P4P', 'X',\n'P5P',\
 'X',\n'P6G', 'X',\n'PA1', 'X',\n'PA2', 'X',\n'PA3\
', 'X',\n'PA4', 'X',\n'PA5', 'X',\n'PAA', 'X',\n'P\
AB', 'X',\n'PAC', 'X',\n'PAD', 'X',\n'PAE', 'X',\n\
'PAG', 'X',\n'PAH', 'X',\n'PAI', 'X',\n'PAL', 'D',\
\n'PAM', 'X',\n'PAN', 'X',\n'PAO', 'X',\n'PAP', 'A\
',\n'PAQ', 'F',\n'PAR', 'X',\n'PAS', 'X',\n'PAT', \
'W',\n'PBA', 'X',\n'PBB', 'X',\n'PBC', 'X',\n'PBF'\
, 'F',\n'PBG', 'X',\n'PBI', 'X',\n'PBM', 'X',\n'PB\
N', 'X',\n'PBP', 'X',\n'PBR', 'X',\n'PBZ', 'X',\n'\
PC2', 'X',\n'PCA', 'E',\n'PCB', 'X',\n'PCD', 'X',\\
n'PCE', 'X',\n'PCG', 'X',\n'PCH', 'X',\n'PCL', 'X'\
,\n'PCM', 'X',\n'PCP', 'X',\n'PCR', 'X',\n'PCS', '\
X',\n'PCU', 'X',\n'PCV', 'X',\n'PCY', 'X',\n'PD1',\
 'X',\n'PDA', 'X',\n'PDC', 'X',\n'PDD', 'A',\n'PDE\
', 'A',\n'PDI', 'X',\n'PDL', 'A',\n'PDN', 'X',\n'P\
DO', 'X',\n'PDP', 'X',\n'PDT', 'X',\n'PDU', 'X',\n\
'PE2', 'X',\n'PE6', 'X',\n'PEA', 'X',\n'PEB', 'X',\
\n'PEC', 'X',\n'PED', 'X',\n'PEE', 'X',\n'PEF', 'X\
',\n'PEG', 'X',\n'PEL', 'X',\n'PEO', 'X',\n'PEP', \
'X',\n'PEQ', 'X',\n'PER', 'X',\n'PET', 'X',\n'PFB'\
, 'X',\n'PFC', 'X',\n'PFG', 'X',\n'PFL', 'X',\n'PF\
M', 'X',\n'PFZ', 'X',\n'PG4', 'X',\n'PG5', 'X',\n'\
PG6', 'X',\n'PGA', 'X',\n'PGC', 'X',\n'PGD', 'X',\\
n'PGE', 'X',\n'PGG', 'G',\n'PGH', 'X',\n'PGL', 'X'\
,\n'PGO', 'X',\n'PGP', 'X',\n'PGQ', 'X',\n'PGR', '\
X',\n'PGS', 'X',\n'PGU', 'X',\n'PGX', 'X',\n'PGY',\
 'G',\n'PH1', 'X',\n'PH2', 'X',\n'PH3', 'X',\n'PHA\
', 'F',\n'PHB', 'X',\n'PHC', 'X',\n'PHD', 'X',\n'P\
HE', 'F',\n'PHG', 'X',\n'PHH', 'X',\n'PHI', 'F',\n\
'PHL', 'F',\n'PHM', 'X',\n'PHN', 'X',\n'PHO', 'X',\
\n'PHP', 'X',\n'PHQ', 'X',\n'PHS', 'H',\n'PHT', 'X\
',\n'PHW', 'P',\n'PHY', 'X',\n'PI1', 'X',\n'PI2', \
'X',\n'PI3', 'X',\n'PI4', 'X',\n'PI5', 'X',\n'PI6'\
, 'X',\n'PI7', 'X',\n'PI8', 'X',\n'PI9', 'X',\n'PI\
A', 'X',\n'PIB', 'X',\n'PIC', 'X',\n'PID', 'X',\n'\
PIG', 'X',\n'PIH', 'X',\n'PIM', 'X',\n'PIN', 'X',\\
n'PIO', 'X',\n'PIP', 'X',\n'PIQ', 'X',\n'PIR', 'X'\
,\n'PIV', 'X',\n'PKF', 'X',\n'PL1', 'X',\n'PL9', '\
X',\n'PLA', 'D',\n'PLC', 'X',\n'PLE', 'L',\n'PLG',\
 'G',\n'PLH', 'X',\n'PLM', 'X',\n'PLP', 'X',\n'PLS\
', 'S',\n'PLT', 'W',\n'PLU', 'L',\n'PLY', 'X',\n'P\
MA', 'X',\n'PMB', 'X',\n'PMC', 'X',\n'PME', 'F',\n\
'PML', 'X',\n'PMM', 'X',\n'PMO', 'X',\n'PMP', 'X',\
\n'PMS', 'X',\n'PMY', 'X',\n'PN2', 'X',\n'PNA', 'X\
',\n'PNB', 'X',\n'PNC', 'G',\n'PND', 'X',\n'PNE', \
'A',\n'PNF', 'X',\n'PNG', 'X',\n'PNI', 'X',\n'PNL'\
, 'X',\n'PNM', 'X',\n'PNN', 'X',\n'PNO', 'X',\n'PN\
P', 'X',\n'PNQ', 'X',\n'PNS', 'X',\n'PNT', 'X',\n'\
PNU', 'X',\n'PO2', 'X',\n'PO4', 'X',\n'POB', 'X',\\
n'POC', 'X',\n'POL', 'X',\n'POM', 'P',\n'PON', 'X'\
,\n'POP', 'X',\n'POR', 'X',\n'POS', 'X',\n'PP1', '\
X',\n'PP2', 'X',\n'PP3', 'A',\n'PP4', 'X',\n'PP5',\
 'X',\n'PP6', 'X',\n'PP7', 'X',\n'PP8', 'N',\n'PP9\
', 'X',\n'PPB', 'X',\n'PPC', 'X',\n'PPD', 'X',\n'P\
PE', 'E',\n'PPG', 'X',\n'PPH', 'F',\n'PPI', 'X',\n\
'PPJ', 'V',\n'PPL', 'X',\n'PPM', 'X',\n'PPN', 'A',\
\n'PPO', 'X',\n'PPP', 'X',\n'PPQ', 'X',\n'PPR', 'X\
',\n'PPS', 'X',\n'PPT', 'X',\n'PPU', 'X',\n'PPX', \
'F',\n'PPY', 'X',\n'PPZ', 'X',\n'PQ0', 'X',\n'PQN'\
, 'X',\n'PQQ', 'X',\n'PR1', 'X',\n'PR2', 'X',\n'PR\
3', 'X',\n'PRA', 'X',\n'PRB', 'X',\n'PRC', 'X',\n'\
PRD', 'X',\n'PRE', 'X',\n'PRF', 'X',\n'PRH', 'X',\\
n'PRI', 'P',\n'PRL', 'X',\n'PRN', 'X',\n'PRO', 'P'\
,\n'PRP', 'X',\n'PRR', 'A',\n'PRS', 'P',\n'PRZ', '\
X',\n'PS0', 'X',\n'PSA', 'X',\n'PSD', 'X',\n'PSE',\
 'X',\n'PSF', 'S',\n'PSG', 'X',\n'PSI', 'X',\n'PSO\
', 'X',\n'PSQ', 'X',\n'PSS', 'X',\n'PST', 'X',\n'P\
SU', 'X',\n'PT1', 'X',\n'PT3', 'X',\n'PTA', 'X',\n\
'PTC', 'X',\n'PTD', 'X',\n'PTE', 'X',\n'PTH', 'Y',\
\n'PTL', 'X',\n'PTM', 'Y',\n'PTN', 'X',\n'PTO', 'X\
',\n'PTP', 'X',\n'PTR', 'Y',\n'PTS', 'X',\n'PTT', \
'X',\n'PTU', 'X',\n'PTY', 'X',\n'PUA', 'X',\n'PUB'\
, 'X',\n'PUR', 'X',\n'PUT', 'X',\n'PVA', 'X',\n'PV\
B', 'X',\n'PVH', 'H',\n'PVL', 'X',\n'PXA', 'X',\n'\
PXF', 'X',\n'PXG', 'X',\n'PXP', 'X',\n'PXY', 'X',\\
n'PXZ', 'X',\n'PY2', 'X',\n'PY4', 'X',\n'PY5', 'X'\
,\n'PY6', 'X',\n'PYA', 'A',\n'PYC', 'X',\n'PYD', '\
X',\n'PYE', 'X',\n'PYL', 'X',\n'PYM', 'X',\n'PYO',\
 'X',\n'PYP', 'X',\n'PYQ', 'X',\n'PYR', 'X',\n'PYS\
', 'X',\n'PYT', 'X',\n'PYX', 'X',\n'PYY', 'X',\n'P\
YZ', 'X',\n'PZQ', 'X',\n'Q82', 'X',\n'QNC', 'X',\n\
'QND', 'X',\n'QSI', 'Q',\n'QTR', 'X',\n'QUA', 'X',\
\n'QUE', 'X',\n'QUI', 'X',\n'QUO', 'X',\n'R11', 'X\
',\n'R12', 'X',\n'R13', 'X',\n'R18', 'X',\n'R1P', \
'X',\n'R56', 'X',\n'R5P', 'X',\n'RA2', 'X',\n'RAD'\
, 'X',\n'RAI', 'X',\n'RAL', 'X',\n'RAM', 'X',\n'RA\
N', 'X',\n'RAP', 'X',\n'RBF', 'X',\n'RBU', 'X',\n'\
RCA', 'X',\n'RCL', 'X',\n'RCO', 'X',\n'RDC', 'X',\\
n'RDF', 'W',\n'RE9', 'X',\n'REA', 'X',\n'RED', 'K'\
,\n'REO', 'X',\n'REP', 'X',\n'RET', 'X',\n'RFA', '\
X',\n'RFB', 'X',\n'RFL', 'X',\n'RFP', 'X',\n'RG1',\
 'X',\n'RGS', 'X',\n'RH1', 'X',\n'RHA', 'X',\n'RHC\
', 'X',\n'RHD', 'X',\n'RHM', 'X',\n'RHO', 'X',\n'R\
HQ', 'X',\n'RHS', 'X',\n'RIA', 'X',\n'RIB', 'X',\n\
'RIC', 'X',\n'RIF', 'X',\n'RIN', 'X',\n'RIP', 'X',\
\n'RIT', 'X',\n'RMB', 'X',\n'RMN', 'X',\n'RMP', 'X\
',\n'RNG', 'X',\n'RNS', 'X',\n'RNT', 'X',\n'RO2', \
'X',\n'RO4', 'X',\n'ROC', 'N',\n'ROI', 'X',\n'ROM'\
, 'X',\n'RON', 'V',\n'ROP', 'X',\n'ROS', 'X',\n'RO\
X', 'X',\n'RPA', 'X',\n'RPD', 'X',\n'RPH', 'X',\n'\
RPL', 'X',\n'RPP', 'X',\n'RPR', 'X',\n'RPX', 'X',\\
n'RQ3', 'X',\n'RR1', 'X',\n'RR6', 'X',\n'RRS', 'X'\
,\n'RS1', 'X',\n'RS2', 'X',\n'RS7', 'X',\n'RSS', '\
X',\n'RTA', 'X',\n'RTB', 'X',\n'RTC', 'X',\n'RTL',\
 'X',\n'RUB', 'X',\n'RUN', 'X',\n'RWJ', 'X',\n'RXP\
', 'X',\n'S02', 'X',\n'S11', 'X',\n'S1H', 'S',\n'S\
27', 'X',\n'S2C', 'C',\n'S3P', 'X',\n'S4U', 'X',\n\
'S57', 'X',\n'S58', 'X',\n'S5H', 'X',\n'S6G', 'X',\
\n'S80', 'X',\n'SAA', 'X',\n'SAB', 'X',\n'SAC', 'S\
',\n'SAD', 'X',\n'SAE', 'X',\n'SAF', 'X',\n'SAH', \
'C',\n'SAI', 'C',\n'SAL', 'X',\n'SAM', 'M',\n'SAN'\
, 'X',\n'SAP', 'X',\n'SAR', 'X',\n'SAS', 'X',\n'SB\
1', 'X',\n'SB2', 'X',\n'SB3', 'X',\n'SB4', 'X',\n'\
SB5', 'X',\n'SB6', 'X',\n'SBA', 'L',\n'SBB', 'X',\\
n'SBD', 'A',\n'SBI', 'X',\n'SBL', 'A',\n'SBN', 'X'\
,\n'SBO', 'X',\n'SBR', 'X',\n'SBS', 'X',\n'SBT', '\
X',\n'SBU', 'X',\n'SBX', 'X',\n'SC4', 'X',\n'SCA',\
 'X',\n'SCC', 'X',\n'SCD', 'X',\n'SCH', 'C',\n'SCI\
', 'X',\n'SCL', 'X',\n'SCM', 'X',\n'SCN', 'X',\n'S\
CO', 'X',\n'SCP', 'S',\n'SCR', 'X',\n'SCS', 'X',\n\
'SCV', 'C',\n'SCY', 'C',\n'SD8', 'X',\n'SDK', 'X',\
\n'SDZ', 'X',\n'SE4', 'X',\n'SEA', 'X',\n'SEB', 'S\
',\n'SEC', 'X',\n'SEG', 'A',\n'SEI', 'X',\n'SEL', \
'S',\n'SEM', 'X',\n'SEO', 'X',\n'SEP', 'S',\n'SER'\
, 'S',\n'SES', 'X',\n'SET', 'S',\n'SEU', 'X',\n'SF\
4', 'X',\n'SFG', 'X',\n'SFN', 'X',\n'SFO', 'X',\n'\
SGA', 'X',\n'SGC', 'X',\n'SGL', 'X',\n'SGM', 'X',\\
n'SGN', 'X',\n'SGP', 'X',\n'SHA', 'X',\n'SHC', 'X'\
,\n'SHF', 'X',\n'SHH', 'X',\n'SHP', 'G',\n'SHR', '\
E',\n'SHT', 'T',\n'SHU', 'X',\n'SI2', 'X',\n'SIA',\
 'X',\n'SIF', 'X',\n'SIG', 'X',\n'SIH', 'X',\n'SIM\
', 'X',\n'SIN', 'X',\n'SKD', 'X',\n'SKF', 'X',\n'S\
LB', 'X',\n'SLE', 'X',\n'SLZ', 'K',\n'SMA', 'X',\n\
'SMC', 'C',\n'SME', 'M',\n'SML', 'X',\n'SMM', 'M',\
\n'SMN', 'X',\n'SMP', 'X',\n'SMS', 'X',\n'SN1', 'X\
',\n'SN6', 'X',\n'SN7', 'X',\n'SNC', 'C',\n'SNN', \
'X',\n'SNP', 'X',\n'SO1', 'X',\n'SO2', 'X',\n'SO3'\
, 'X',\n'SO4', 'X',\n'SOA', 'X',\n'SOC', 'C',\n'SO\
M', 'X',\n'SOR', 'X',\n'SOT', 'X',\n'SOX', 'X',\n'\
SPA', 'X',\n'SPB', 'X',\n'SPC', 'X',\n'SPD', 'X',\\
n'SPE', 'X',\n'SPG', 'X',\n'SPH', 'X',\n'SPI', 'X'\
,\n'SPK', 'X',\n'SPM', 'X',\n'SPN', 'X',\n'SPO', '\
X',\n'SPP', 'X',\n'SPS', 'X',\n'SPY', 'X',\n'SQU',\
 'X',\n'SRA', 'X',\n'SRB', 'X',\n'SRD', 'X',\n'SRL\
', 'X',\n'SRM', 'X',\n'SRS', 'X',\n'SRY', 'X',\n'S\
SA', 'X',\n'SSB', 'X',\n'SSG', 'X',\n'SSP', 'X',\n\
'ST1', 'X',\n'ST2', 'X',\n'ST3', 'X',\n'ST4', 'X',\
\n'ST5', 'X',\n'ST6', 'X',\n'STA', 'X',\n'STB', 'X\
',\n'STE', 'X',\n'STG', 'X',\n'STI', 'X',\n'STL', \
'X',\n'STN', 'X',\n'STO', 'X',\n'STP', 'X',\n'STR'\
, 'X',\n'STU', 'X',\n'STY', 'Y',\n'SU1', 'X',\n'SU\
2', 'X',\n'SUC', 'X',\n'SUI', 'X',\n'SUL', 'X',\n'\
SUR', 'X',\n'SVA', 'S',\n'SWA', 'X',\n'T16', 'X',\\
n'T19', 'X',\n'T23', 'X',\n'T29', 'X',\n'T33', 'X'\
,\n'T3P', 'X',\n'T42', 'A',\n'T44', 'X',\n'T5A', '\
X',\n'T6A', 'T',\n'T6P', 'X',\n'T80', 'X',\n'T87',\
 'X',\n'TA1', 'X',\n'TAA', 'X',\n'TAB', 'X',\n'TAC\
', 'X',\n'TAD', 'X',\n'TAF', 'X',\n'TAM', 'X',\n'T\
AP', 'X',\n'TAR', 'X',\n'TAS', 'X',\n'TAU', 'X',\n\
'TAX', 'X',\n'TAZ', 'X',\n'TB9', 'X',\n'TBA', 'X',\
\n'TBD', 'X',\n'TBG', 'G',\n'TBH', 'X',\n'TBM', 'T\
',\n'TBO', 'X',\n'TBP', 'X',\n'TBR', 'X',\n'TBS', \
'X',\n'TBT', 'X',\n'TBU', 'X',\n'TBZ', 'X',\n'TC4'\
, 'X',\n'TCA', 'X',\n'TCB', 'X',\n'TCH', 'X',\n'TC\
K', 'X',\n'TCL', 'X',\n'TCM', 'X',\n'TCN', 'X',\n'\
TCP', 'X',\n'TCR', 'W',\n'TCS', 'X',\n'TCZ', 'X',\\
n'TDA', 'X',\n'TDB', 'X',\n'TDG', 'X',\n'TDP', 'X'\
,\n'TDR', 'X',\n'TDX', 'X',\n'TEA', 'X',\n'TEM', '\
X',\n'TEN', 'X',\n'TEO', 'X',\n'TEP', 'X',\n'TER',\
 'X',\n'TES', 'X',\n'TET', 'X',\n'TFA', 'X',\n'TFB\
', 'X',\n'TFH', 'X',\n'TFI', 'X',\n'TFK', 'X',\n'T\
FP', 'X',\n'THA', 'X',\n'THB', 'X',\n'THC', 'T',\n\
'THD', 'X',\n'THE', 'X',\n'THF', 'X',\n'THJ', 'X',\
\n'THK', 'X',\n'THM', 'X',\n'THN', 'X',\n'THO', 'T\
',\n'THP', 'X',\n'THQ', 'X',\n'THR', 'T',\n'THS', \
'X',\n'THT', 'X',\n'THU', 'X',\n'THX', 'X',\n'THZ'\
, 'X',\n'TI1', 'X',\n'TI2', 'X',\n'TI3', 'P',\n'TI\
A', 'X',\n'TIH', 'A',\n'TK4', 'X',\n'TLA', 'X',\n'\
TLC', 'X',\n'TLM', 'X',\n'TLN', 'X',\n'TLX', 'X',\\
n'TM5', 'X',\n'TM6', 'X',\n'TMA', 'X',\n'TMB', 'T'\
,\n'TMC', 'X',\n'TMD', 'T',\n'TME', 'X',\n'TMF', '\
X',\n'TML', 'K',\n'TMM', 'X',\n'TMN', 'X',\n'TMP',\
 'X',\n'TMQ', 'X',\n'TMR', 'X',\n'TMT', 'X',\n'TMZ\
', 'X',\n'TNB', 'C',\n'TND', 'X',\n'TNK', 'X',\n'T\
NP', 'X',\n'TNT', 'X',\n'TOA', 'X',\n'TOB', 'X',\n\
'TOC', 'X',\n'TOL', 'X',\n'TOP', 'X',\n'TOS', 'X',\
\n'TOT', 'X',\n'TP1', 'G',\n'TP2', 'P',\n'TP3', 'E\
',\n'TP4', 'E',\n'TP7', 'T',\n'TPA', 'X',\n'TPE', \
'X',\n'TPF', 'X',\n'TPI', 'X',\n'TPL', 'W',\n'TPM'\
, 'X',\n'TPN', 'G',\n'TPO', 'T',\n'TPP', 'X',\n'TP\
Q', 'A',\n'TPR', 'P',\n'TPS', 'X',\n'TPT', 'X',\n'\
TPV', 'X',\n'TPX', 'X',\n'TPY', 'X',\n'TQ3', 'X',\\
n'TQ4', 'X',\n'TQ5', 'X',\n'TQ6', 'X',\n'TR1', 'X'\
,\n'TRA', 'X',\n'TRB', 'X',\n'TRC', 'X',\n'TRD', '\
X',\n'TRE', 'X',\n'TRF', 'W',\n'TRG', 'K',\n'TRH',\
 'X',\n'TRI', 'X',\n'TRJ', 'X',\n'TRM', 'X',\n'TRN\
', 'W',\n'TRO', 'W',\n'TRP', 'W',\n'TRQ', 'X',\n'T\
RS', 'X',\n'TRX', 'W',\n'TRZ', 'X',\n'TS2', 'X',\n\
'TS3', 'X',\n'TS4', 'X',\n'TS5', 'X',\n'TSA', 'X',\
\n'TSB', 'X',\n'TSI', 'X',\n'TSM', 'X',\n'TSN', 'X\
',\n'TSP', 'X',\n'TSU', 'X',\n'TTA', 'X',\n'TTE', \
'X',\n'TTN', 'X',\n'TTO', 'X',\n'TTP', 'X',\n'TTX'\
, 'X',\n'TXL', 'X',\n'TYA', 'Y',\n'TYB', 'Y',\n'TY\
D', 'X',\n'TYI', 'Y',\n'TYL', 'X',\n'TYM', 'W',\n'\
TYN', 'Y',\n'TYQ', 'Y',\n'TYR', 'Y',\n'TYS', 'Y',\\
n'TYV', 'X',\n'TYY', 'A',\n'TZB', 'X',\n'TZC', 'X'\
,\n'TZE', 'X',\n'TZL', 'X',\n'TZO', 'X',\n'TZP', '\
X',\n'U01', 'X',\n'U02', 'X',\n'U03', 'X',\n'U04',\
 'X',\n'U05', 'X',\n'U0E', 'X',\n'U10', 'X',\n'U18\
', 'X',\n'U2G', 'X',\n'U3P', 'X',\n'U49', 'X',\n'U\
55', 'X',\n'U5P', 'X',\n'U66', 'X',\n'U89', 'X',\n\
'U8U', 'X',\n'UAA', 'X',\n'UAG', 'A',\n'UAP', 'X',\
\n'UAR', 'X',\n'UC1', 'X',\n'UC2', 'X',\n'UC3', 'X\
',\n'UC4', 'X',\n'UD1', 'X',\n'UD2', 'X',\n'UDP', \
'X',\n'UDX', 'X',\n'UFG', 'X',\n'UFM', 'X',\n'UFP'\
, 'X',\n'UGA', 'X',\n'UIN', 'X',\n'UKP', 'A',\n'UM\
3', 'X',\n'UMA', 'A',\n'UMG', 'X',\n'UMP', 'X',\n'\
UNA', 'X',\n'UND', 'X',\n'UNI', 'X',\n'UNK', 'X',\\
n'UNN', 'X',\n'UNX', 'X',\n'UP5', 'X',\n'UP6', 'X'\
,\n'UPA', 'X',\n'UPF', 'X',\n'UPG', 'X',\n'UPP', '\
X',\n'UQ1', 'X',\n'UQ2', 'X',\n'UQ6', 'X',\n'UR2',\
 'X',\n'URA', 'X',\n'URE', 'X',\n'URF', 'X',\n'URI\
', 'X',\n'URS', 'X',\n'UTP', 'X',\n'UVC', 'X',\n'U\
VW', 'X',\n'V35', 'X',\n'V36', 'X',\n'V4O', 'X',\n\
'V7O', 'X',\n'VAA', 'V',\n'VAC', 'X',\n'VAD', 'V',\
\n'VAF', 'V',\n'VAG', 'X',\n'VAL', 'V',\n'VAN', 'X\
',\n'VAS', 'X',\n'VAX', 'X',\n'VDX', 'X',\n'VDY', \
'X',\n'VG1', 'X',\n'VIB', 'X',\n'VIR', 'X',\n'VIT'\
, 'X',\n'VK3', 'X',\n'VO3', 'X',\n'VO4', 'X',\n'VS\
1', 'F',\n'VS2', 'F',\n'VS3', 'F',\n'VS4', 'F',\n'\
VXA', 'X',\n'W01', 'X',\n'W02', 'X',\n'W03', 'X',\\
n'W11', 'X',\n'W33', 'X',\n'W35', 'X',\n'W42', 'X'\
,\n'W43', 'X',\n'W54', 'X',\n'W56', 'X',\n'W59', '\
X',\n'W71', 'X',\n'W84', 'X',\n'W8R', 'X',\n'W91',\
 'X',\n'WAY', 'X',\n'WCC', 'X',\n'WO2', 'X',\n'WO4\
', 'X',\n'WRB', 'X',\n'WRR', 'X',\n'WRS', 'X',\n'W\
W7', 'X',\n'X2F', 'X',\n'X7O', 'X',\n'XAA', 'X',\n\
'XAN', 'X',\n'XAO', 'X',\n'XBB', 'X',\n'XBP', 'X',\
\n'XDN', 'X',\n'XDP', 'X',\n'XIF', 'X',\n'XIM', 'X\
',\n'XK2', 'X',\n'XL1', 'X',\n'XLS', 'X',\n'XMP', \
'X',\n'XN1', 'X',\n'XN2', 'X',\n'XN3', 'X',\n'XUL'\
, 'X',\n'XV6', 'X',\n'XYD', 'X',\n'XYH', 'X',\n'XY\
L', 'X',\n'XYP', 'X',\n'XYS', 'X',\n'YOF', 'Y',\n'\
YRR', 'X',\n'YT3', 'X',\n'YZ9', 'X',\n'Z34', 'G',\\
n'Z5A', 'X',\n'ZAF', 'X',\n'ZAP', 'X',\n'ZEB', 'X'\
,\n'ZEN', 'X',\n'ZES', 'X',\n'ZID', 'X',\n'ZMR', '\
X',\n'ZN3', 'X',\n'ZNH', 'X',\n'ZNO', 'X',\n'ZO3',\
 'X',\n'ZPR', 'P',\n'ZRA', 'A',\n'ZST', 'X',\n'ZYA\
', 'A',\n\n\n'ASN','N');\n} \n\n\nsub file2head\n \
     {\n	my $file = shift;\n	my $size = shift;\n	m\
y $f= new FileHandle;\n	my $line;\n	open ($f,$file\
);\n	read ($f,$line, $size);\n	close ($f);\n	retur\
n $line;\n      }\nsub file2tail\n      {\n	my $fi\
le = shift;\n	my $size = shift;\n	my $f= new FileH\
andle;\n	my $line;\n	\n	open ($f,$file);\n	seek ($\
f,$size*-1, 2);\n	read ($f,$line, $size);\n	close \
($f);\n	return $line;\n      }\n\n\nsub vtmpnam\n \
     {\n	my $r=rand(100000);\n	my $f=\"file.$r.$$\\
";\n	while (-e $f)\n	  {\n	    $f=vtmpnam();\n	  }\
\n	push (@TMPFILE_LIST, $f);\n	return $f;\n      }\
\n\nsub myexit\n  {\n    my $code=@_[0];\n    if (\
$CLEAN_EXIT_STARTED==1){return;}\n    else {$CLEAN\
_EXIT_STARTED=1;}\n    ### ONLY BARE EXIT\n    exi\
t ($code);\n  }\nsub set_error_lock\n    {\n      \
my $name = shift;\n      my $pid=$$;\n\n      \n  \
    &lock4tc ($$,\"LERROR\", \"LSET\", \"$$ -- ERR\
OR: $name $PROGRAM\\n\");\n      return;\n    }\ns\
ub set_lock\n  {\n    my $pid=shift;\n    my $msg=\
 shift;\n    my $p=getppid();\n    &lock4tc ($pid,\
\"LLOCK\",\"LRESET\",\"$p$msg\\n\");\n  }\nsub uns\
et_lock\n   {\n     \n    my $pid=shift;\n    &loc\
k4tc ($pid,\"LLOCK\",\"LRELEASE\",\"\");\n  }\nsub\
 shift_lock\n  {\n    my $from=shift;\n    my $to=\
shift;\n    my $from_type=shift;\n    my $to_type=\
shift;\n    my $action=shift;\n    my $msg;\n    \\
n    if (!&lock4tc($from, $from_type, \"LCHECK\", \
\"\")){return 0;}\n    $msg=&lock4tc ($from, $from\
_type, \"LREAD\", \"\");\n    &lock4tc ($from, $fr\
om_type,\"LRELEASE\", $msg);\n    &lock4tc ($to, $\
to_type, $action, $msg);\n    return;\n  }\nsub is\
shellpid\n  {\n    my $p=shift;\n    if (!lock4tc \
($p, \"LLOCK\", \"LCHECK\")){return 0;}\n    else\\
n      {\n	my $c=lock4tc($p, \"LLOCK\", \"LREAD\")\
;\n	if ( $c=~/-SHELL-/){return 1;}\n      }\n    r\
eturn 0;\n  }\nsub isrootpid\n  {\n    if(lock4tc \
(getppid(), \"LLOCK\", \"LCHECK\")){return 0;}\n  \
  else {return 1;}\n  }\nsub lock4tc\n	{\n	  my ($\
pid,$type,$action,$value)=@_;\n	  my $fname;\n	  m\
y $host=hostname;\n	  \n	  if ($type eq \"LLOCK\")\
{$fname=\"$LOCKDIR/.$pid.$host.lock4tcoffee\";}\n	\
  elsif ( $type eq \"LERROR\"){ $fname=\"$LOCKDIR/\
.$pid.$host.error4tcoffee\";}\n	  elsif ( $type eq\
 \"LWARNING\"){ $fname=\"$LOCKDIR/.$pid.$host.warn\
ing4tcoffee\";}\n	  \n	  if ($debug_lock)\n	    {\\
n	      print STDERR \"\\n\\t---lock4tc(tcg): $act\
ion => $fname =>$value (RD: $LOCKDIR)\\n\";\n	    \
}\n\n	  if    ($action eq \"LCHECK\") {return -e $\
fname;}\n	  elsif ($action eq \"LREAD\"){return fi\
le2string($fname);}\n	  elsif ($action eq \"LSET\"\
) {return string2file ($value, $fname, \">>\");}\n\
	  elsif ($action eq \"LRESET\") {return string2fi\
le ($value, $fname, \">\");}\n	  elsif ($action eq\
 \"LRELEASE\") \n	    {\n	      if ( $debug_lock)\\
n		{\n		  my $g=new FileHandle;\n		  open ($g, \">\
>$fname\");\n		  print $g \"\\nDestroyed by $$\\n\\
";\n		  close ($g);\n		  safe_system (\"mv $fname \
$fname.old\");\n		}\n	      else\n		{\n		  unlink \
($fname);\n		}\n	    }\n	  return \"\";\n	}\n	\nsu\
b file2string\n	{\n	  my $file=@_[0];\n	  my $f=ne\
w FileHandle;\n	  my $r;\n	  open ($f, \"$file\");\
\n	  while (<$f>){$r.=$_;}\n	  close ($f);\n	  ret\
urn $r;\n	}\nsub string2file \n    {\n    my ($s,$\
file,$mode)=@_;\n    my $f=new FileHandle;\n    \n\
    open ($f, \"$mode$file\");\n    print $f  \"$s\
\";\n    close ($f);\n  }\n\nBEGIN\n    {\n      s\
rand;\n    \n      $SIG{'SIGUP'}='signal_cleanup';\
\n      $SIG{'SIGINT'}='signal_cleanup';\n      $S\
IG{'SIGQUIT'}='signal_cleanup';\n      $SIG{'SIGIL\
L'}='signal_cleanup';\n      $SIG{'SIGTRAP'}='sign\
al_cleanup';\n      $SIG{'SIGABRT'}='signal_cleanu\
p';\n      $SIG{'SIGEMT'}='signal_cleanup';\n     \
 $SIG{'SIGFPE'}='signal_cleanup';\n      \n      $\
SIG{'SIGKILL'}='signal_cleanup';\n      $SIG{'SIGP\
IPE'}='signal_cleanup';\n      $SIG{'SIGSTOP'}='si\
gnal_cleanup';\n      $SIG{'SIGTTIN'}='signal_clea\
nup';\n      $SIG{'SIGXFSZ'}='signal_cleanup';\n  \
    $SIG{'SIGINFO'}='signal_cleanup';\n      \n   \
   $SIG{'SIGBUS'}='signal_cleanup';\n      $SIG{'S\
IGALRM'}='signal_cleanup';\n      $SIG{'SIGTSTP'}=\
'signal_cleanup';\n      $SIG{'SIGTTOU'}='signal_c\
leanup';\n      $SIG{'SIGVTALRM'}='signal_cleanup'\
;\n      $SIG{'SIGUSR1'}='signal_cleanup';\n\n\n  \
    $SIG{'SIGSEGV'}='signal_cleanup';\n      $SIG{\
'SIGTERM'}='signal_cleanup';\n      $SIG{'SIGCONT'\
}='signal_cleanup';\n      $SIG{'SIGIO'}='signal_c\
leanup';\n      $SIG{'SIGPROF'}='signal_cleanup';\\
n      $SIG{'SIGUSR2'}='signal_cleanup';\n\n      \
$SIG{'SIGSYS'}='signal_cleanup';\n      $SIG{'SIGU\
RG'}='signal_cleanup';\n      $SIG{'SIGCHLD'}='sig\
nal_cleanup';\n      $SIG{'SIGXCPU'}='signal_clean\
up';\n      $SIG{'SIGWINCH'}='signal_cleanup';\n  \
    \n      $SIG{'INT'}='signal_cleanup';\n      $\
SIG{'TERM'}='signal_cleanup';\n      $SIG{'KILL'}=\
'signal_cleanup';\n      $SIG{'QUIT'}='signal_clea\
nup';\n      \n      our $debug_lock=$ENV{\"DEBUG_\
LOCK\"};\n      \n      \n      \n      \n      fo\
reach my $a (@ARGV){$CL.=\" $a\";}\n      if ( $de\
bug_lock ){print STDERR \"\\n\\n\\n********** STAR\
T PG: $PROGRAM *************\\n\";}\n      if ( $d\
ebug_lock ){print STDERR \"\\n\\n\\n**********(tcg\
) LOCKDIR: $LOCKDIR $$ *************\\n\";}\n     \
 if ( $debug_lock ){print STDERR \"\\n --- $$ -- $\
CL\\n\";}\n      \n	     \n      \n      \n    }\n\
sub flush_error\n  {\n    my $msg=shift;\n    retu\
rn add_error ($EXIT_FAILURE,$$, $$,getppid(), $msg\
, $CL);\n  }\nsub add_error \n  {\n    my $code=sh\
ift;\n    my $rpid=shift;\n    my $pid=shift;\n   \
 my $ppid=shift;\n    my $type=shift;\n    my $com\
=shift;\n    \n    $ERROR_DONE=1;\n    lock4tc ($r\
pid, \"LERROR\",\"LSET\",\"$pid -- ERROR: $type\\n\
\");\n    lock4tc ($$, \"LERROR\",\"LSET\", \"$pid\
 -- COM: $com\\n\");\n    lock4tc ($$, \"LERROR\",\
\"LSET\", \"$pid -- STACK: $ppid -> $pid\\n\");\n \
  \n    return $code;\n  }\nsub add_warning \n  {\\
n    my $rpid=shift;\n    my $pid =shift;\n    my \
$command=shift;\n    my $msg=\"$$ -- WARNING: $com\
mand\\n\";\n    print STDERR \"$msg\";\n    lock4t\
c ($$, \"LWARNING\", \"LSET\", $msg);\n  }\n\nsub \
signal_cleanup\n  {\n    print dtderr \"\\n**** $$\
 (tcg) was killed\\n\";\n    &cleanup;\n    exit (\
$EXIT_FAILURE);\n  }\nsub clean_dir\n  {\n    my $\
dir=@_[0];\n    if ( !-d $dir){return ;}\n    elsi\
f (!($dir=~/tmp/)){return ;}#safety check 1\n    e\
lsif (($dir=~/\\*/)){return ;}#safety check 2\n   \
 else\n      {\n	`rm -rf $dir`;\n      }\n    retu\
rn;\n  }\nsub cleanup\n  {\n    #print stderr \"\\\
n----tc: $$ Kills $PIDCHILD\\n\";\n    #kill (SIGT\
ERM,$PIDCHILD);\n    my $p=getppid();\n    $CLEAN_\
EXIT_STARTED=1;\n    \n    \n    \n    if (&lock4t\
c($$,\"LERROR\", \"LCHECK\", \"\"))\n      {\n	my \
$ppid=getppid();\n	if (!$ERROR_DONE) \n	  {\n	    \
&lock4tc($$,\"LERROR\", \"LSET\", \"$$ -- STACK: $\
p -> $$\\n\");\n	    &lock4tc($$,\"LERROR\", \"LSE\
T\", \"$$ -- COM: $CL\\n\");\n	  }\n      }\n    m\
y $warning=&lock4tc($$, \"LWARNING\", \"LREAD\", \\
"\");\n    my $error=&lock4tc($$,  \"LERROR\", \"L\
READ\", \"\");\n    #release error and warning loc\
k if root\n    \n    if (isrootpid() && ($warning \
|| $error) )\n      {\n	\n	print STDERR \"********\
******** Summary *************\\n$error\\n$warning\
\\n\";\n\n	&lock4tc($$,\"LERROR\",\"RELEASE\",\"\"\
);\n	&lock4tc($$,\"LWARNING\",\"RELEASE\",\"\");\n\
      } \n    \n    \n    foreach my $f (@TMPFILE_\
LIST)\n      {\n	if (-e $f){unlink ($f);} \n      \
}\n    foreach my $d (@TMPDIR_LIST)\n      {\n	cle\
an_dir ($d);\n      }\n    #No More Lock Release\n\
    #&lock4tc($$,\"LLOCK\",\"LRELEASE\",\"\"); #re\
lease lock \n\n    if ( $debug_lock ){print STDERR\
 \"\\n\\n\\n********** END PG: $PROGRAM ($$) *****\
********\\n\";}\n    if ( $debug_lock ){print STDE\
RR \"\\n\\n\\n**********(tcg) LOCKDIR: $LOCKDIR $$\
 *************\\n\";}\n  }\nEND \n  {\n    \n    &\
cleanup();\n  }\n   \n\nsub safe_system \n{\n  my \
$com=shift;\n  my $ntry=shift;\n  my $ctry=shift;\\
n  my $pid;\n  my $status;\n  my $ppid=getppid();\\
n  if ($com eq \"\"){return 1;}\n  \n  \n\n  if ((\
$pid = fork ()) < 0){return (-1);}\n  if ($pid == \
0)\n    {\n      set_lock($$, \" -SHELL- $com (tcg\
)\");\n      exec ($com);\n    }\n  else\n    {\n \
     lock4tc ($$, \"LLOCK\", \"LSET\", \"$pid\\n\"\
);#update parent\n      $PIDCHILD=$pid;\n    }\n  \
if ($debug_lock){printf STDERR \"\\n\\t .... safe_\
system (fasta_seq2hmm)  p: $$ c: $pid COM: $com\\n\
\";}\n\n  waitpid ($pid,WTERMSIG);\n\n  shift_lock\
 ($pid,$$, \"LWARNING\",\"LWARNING\", \"LSET\");\n\
\n  if ($? == $EXIT_FAILURE || lock4tc($pid, \"LER\
ROR\", \"LCHECK\", \"\"))\n    {\n      if ($ntry \
&& $ctry <$ntry)\n	{\n	  add_warning ($$,$$,\"$com\
 failed [retry: $ctry]\");\n	  lock4tc ($pid, \"LR\
ELEASE\", \"LERROR\", \"\");\n	  return safe_syste\
m ($com, $ntry, ++$ctry);\n	}\n      elsif ($ntry \
== -1)\n	{\n	  if (!shift_lock ($pid, $$, \"LERROR\
\", \"LWARNING\", \"LSET\"))\n	    {\n	      add_w\
arning ($$,$$,\"$com failed\");\n	    }\n	  else\n\
	    {\n	      lock4tc ($pid, \"LRELEASE\", \"LERR\
OR\", \"\");\n	    }\n	  return $?;}\n      else\n\
	{\n	  if (!shift_lock ($pid,$$, \"LERROR\",\"LERR\
OR\", \"LSET\"))\n	    {\n	      myexit(add_error \
($EXIT_FAILURE,$$,$pid,getppid(), \"UNSPECIFIED sy\
stem\", $com));\n	    }\n	}\n    }\n  return $?;\n\
}\n\nsub check_configuration \n    {\n      my @l=\
@_;\n      my $v;\n      foreach my $p (@l)\n	{\n	\
  \n	  if   ( $p eq \"EMAIL\")\n	    { \n	      if\
 ( !($EMAIL=~/@/))\n		{\n		add_warning($$,$$,\"Cou\
ld Not Use EMAIL\");\n		myexit(add_error ($EXIT_FA\
ILURE,$$,$$,getppid(),\"EMAIL\",\"$CL\"));\n	     \
 }\n	    }\n	  elsif( $p eq \"INTERNET\")\n	    {\\
n	      if ( !&check_internet_connection())\n		{\n\
		  myexit(add_error ($EXIT_FAILURE,$$,$$,getppid(\
),\"INTERNET\",\"$CL\"));\n		}\n	    }\n	  elsif( \
$p eq \"wget\")\n	    {\n	      if (!&pg_is_instal\
led (\"wget\") && !&pg_is_installed (\"curl\"))\n	\
	{\n		  myexit(add_error ($EXIT_FAILURE,$$,$$,getp\
pid(),\"PG_NOT_INSTALLED:wget\",\"$CL\"));\n		}\n	\
    }\n	  elsif( !(&pg_is_installed ($p)))\n	    {\
\n	      myexit(add_error ($EXIT_FAILURE,$$,$$,get\
ppid(),\"PG_NOT_INSTALLED:$p\",\"$CL\"));\n	    }\\
n	}\n      return 1;\n    }\nsub pg_is_installed\n\
  {\n    my @ml=@_;\n    my $r, $p, $m;\n    my $s\
upported=0;\n    \n    my $p=shift (@ml);\n    if \
($p=~/::/)\n      {\n	if (safe_system (\"perl -M$p\
 -e 1\")==$EXIT_SUCCESS){return 1;}\n	else {return\
 0;}\n      }\n    else\n      {\n	$r=`which $p 2>\
/dev/null`;\n	if ($r eq \"\"){return 0;}\n	else {r\
eturn 1;}\n      }\n  }\n\n\n\nsub check_internet_\
connection\n  {\n    my $internet;\n    my $tmp;\n\
    &check_configuration ( \"wget\"); \n    \n    \
$tmp=&vtmpnam ();\n    \n    if     (&pg_is_instal\
led    (\"wget\")){`wget www.google.com -O$tmp >/d\
ev/null 2>/dev/null`;}\n    elsif  (&pg_is_install\
ed    (\"curl\")){`curl www.google.com -o$tmp >/de\
v/null 2>/dev/null`;}\n    \n    if ( !-e $tmp || \
-s $tmp < 10){$internet=0;}\n    else {$internet=1\
;}\n    if (-e $tmp){unlink $tmp;}\n\n    return $\
internet;\n  }\nsub check_pg_is_installed\n  {\n  \
  my @ml=@_;\n    my $r=&pg_is_installed (@ml);\n \
   if (!$r && $p=~/::/)\n      {\n	print STDERR \"\
\\nYou Must Install the perl package $p on your sy\
stem.\\nRUN:\\n\\tsudo perl -MCPAN -e 'install $pg\
'\\n\";\n      }\n    elsif (!$r)\n      {\n	myexi\
t(flush_error(\"\\nProgram $p Supported but Not In\
stalled on your system\"));\n      }\n    else\n  \
    {\n	return 1;\n      }\n  }\n\n\n","use Cwd;\n\
use Env;\nuse File::Path;\nuse FileHandle;\nuse st\
rict;\n\n\nour (%MODE, %PG, %ENV_SET, %SUPPORTED_O\
S);\n\n\nour $EXIT_SUCCESS=0;\nour $EXIT_FAILURE=1\
;\nour $INTERNET=0;\n\nour $CP=\"cp \"; #was causi\
ng a crash on MacOSX\nour $SILENT=\">/dev/null 2>/\
dev/null\";\nour $WEB_BASE=\"http://www.tcoffee.or\
g\";\nour $TCLINKDB_ADDRESS=\"$WEB_BASE/Resources/\
tclinkdb.txt\";\nour $OS=get_os();\nour $ROOT=&get\
_root();\nour $CD=cwd();\nour $CDIR=$CD;\nour $HOM\
E=$ENV{'HOME'};\n\nour $OSNAME=$ENV{'OSNAME'};\nou\
r $OSARCH=$ENV{'OSARCH'};\nour $REPO_ROOT=\"\";\n\\
nour $TCDIR;\nour $TCCACHE;\nour $TCTMP;\nour $TCM\
;\nour $TCMETHODS;\nour $TCPLUGINS;\nour $PLUGINS_\
DIR=\"\";\nour $INSTALL_DIR=\"\";\n\nour $CXX=\"g+\
+\";\nour $CXXFLAGS=\"\";\n\nour $CPP=\"g++\";\nou\
r $CPPFLAGS=\"\";\n\nour $CC=\"gcc\";\nour $CFLAGS\
=\"\";\n\nour $FC=\"f77\";\nour $FFLAGS=\"\";\n\nm\
y $install=\"all\";\nmy $default_update_action=\"n\
o_update\";\nmy @required_applications=(\"wget_OR_\
curl\");\nmy @smode=(\"all\", \"clean\", \"install\
\");\n\n&initialize_PG();\n\nmy $cl=join( \" \", @\
ARGV);\nif ($#ARGV==-1 || ($cl=~/-h/) ||($cl=~/-H/\
) )\n  {\n     print \"\\n!!!!!!! ./install  t_cof\
fee             --> installs t_coffee only\";\n   \
  print \"\\n!!!!!!! ./install  all               \
   --> installs all the modes [mcoffee, expresso, \
psicoffee,rcoffee..]\";\n     print \"\\n!!!!!!! .\
/install  [mcoffee|rcoffee|..] --> installs the sp\
ecified mode\";\n     print \"\\n!!!!!!! ./install\
  -h                   --> print usage\\n\\n\";\n \
    if ( $#ARGV==-1){exit ($EXIT_FAILURE);}\n   }\\
n     \nif (($cl=~/-h/) ||($cl=~/-H/) )\n  {\n    \
my $m;\n    print \"\\n\\n!!!!!!! advanced mode\\n\
\";\n    foreach $m ((keys (%MODE)),@smode)\n     \
 {\n	print \"!!!!!!!       ./install $m\\n\";\n   \
   }\n    \n    print \"!!!!!!! ./install [target:\
package|mode|] [-update|-force|-exec=dir|-dis=dir|\
-root|-tclinkdb=file|-] [CC=|FCC=|CXX=|CFLAGS=|CXX\
FLAGS=]\\n\";\n    print \"!!!!!!! ./install clean\
    [removes all executables]\\n\";\n    print \"!\
!!!!!! ./install [optional:target] -update        \
       [updates package already installed]\\n\";\n\
    print \"!!!!!!! ./install [optional:target] -f\
orce                [Forces recompilation over eve\
rything]\\n\";\n    \n    print \"!!!!!!! ./instal\
l [optional:target] -root                 [You are\
 running as root]\\n\";\n    print \"!!!!!!! ./ins\
tall [optional:target] -exec=/foo/bar/       [addr\
ess for the T-Coffee executable]\\n\";\n    print \
\"!!!!!!! ./install [optional:target] -dis=/foo/ba\
r/        [Address where distributions should be s\
tored]\\n\";\n    print \"!!!!!!! ./install [optio\
nal:target] -tclinkdb=foo|update  [file containing\
 all the packages to be installed]\\n\";\n    prin\
t \"!!!!!!! ./install [optional:target] -clean    \
            [clean everything]\\n\";\n    print \"\
!!!!!!! ./install [optional:target] -plugins      \
        [plugins directory]\\n\";\n    print \"!!!\
!!!! ./install [optional:target] -tcdir=/foor/bar \
     [base path where T-Coffee will be installed]\\
\n\";\n    print \"!!!!!!! ./install [optional:tar\
get] -repo=/path/to/repo   [binaries repository ro\
ot directory]\\n\";\n    print \"!!!!!!! mode:\";\\
n    foreach $m (keys(%MODE)){print \"$m \";}\n   \
 print \"\\n\";\n    print \"!!!!!!! Packages:\";\\
n    foreach $m (keys (%PG)){print \"$m \";}\n    \
print \"\\n\";\n    \n    print \"\\n\\n\";\n    e\
xit ($EXIT_FAILURE);\n  }\n\n\n\nmy (@argl)=($cl=~\
/(\\S+=[^=]+)\\s\\w+=/g);\npush (@argl, ($cl=~/(\\\
S+=[^=]+\\S)\\s*$/g));\n\nforeach $a (@argl)\n  {\\
n    if ( ($cl=~/CXX=(.*)/)){$CXX=$1;}\n    if ( (\
$cl=~/-CC=(.*)/    )){$CC=$1;}\n    if ( ($cl=~/-F\
C=(.*)/    )){$FC=$1;}\n    if ( ($cl=~/-CFLAGS=(.\
*)/)){$CFLAGS=$1;}\n    if ( ($cl=~/-CXXFLAGS=(.*)\
/)){$CXXFLAGS=$1;}\n  }\nour ($ROOT_INSTALL, $NO_Q\
UESTION, $default_update_action,$BINARIES_ONLY,$fo\
rce, $default_update_action, $INSTALL_DIR, $PLUGIN\
S_DIR, $DISTRIBUTIONS,$tclinkdb, $proxy, $clean);\\
nif ( ($cl=~/-root/)){$ROOT_INSTALL=1;}\nif ( ($cl\
=~/-no_question/)){$NO_QUESTION=1;}\nif ( ($cl=~/-\
update/)){$default_update_action=\"update\";}\nif \
( ($cl=~/-binaries/)){$BINARIES_ONLY=1;}\nif ( ($c\
l=~/-force/)){$force=1;$default_update_action=\"up\
date\"}\nif ( ($cl=~/-exec=\\s*(\\S+)/)){$INSTALL_\
DIR=$1;}\nif ( ($cl=~/-plugins=\\s*(\\S+)/)){$PLUG\
INS_DIR=$1;}\nif ( ($cl=~/-dis=\\s*(\\S+)/)){$DIST\
RIBUTIONS=$1;}\n\nif ( ($cl=~/-tclinkdb=\\s*(\\S+)\
/)){$tclinkdb=$1;}\nif ( ($cl=~/-proxy=\\s*(\\S+)/\
)){$proxy=$1;}\nif ( ($cl=~/-clean/)){$clean=1;}\n\
if ( ($cl=~/-repo=\\s*(\\S+)/)){ $REPO_ROOT=$1; }\\
nif ( ($cl=~/-tcdir=\\s*(\\S+)/)){ $TCDIR=$1; }\ni\
f ($tclinkdb){&update_tclinkdb ($tclinkdb);}\n\n\n\
if( $REPO_ROOT ne \"\" ) {\n	if( $OSNAME eq \"\" )\
 { print \"You have specified the repository folde\
r but the required \\\"OSNAME\\\" enviroment varia\
ble is missing. \\n\"; exit 1; } \n	if( $OSARCH eq\
 \"\" ) { print \"You have specified the repositor\
y folder but the required \\\"OSARCH\\\" enviromen\
t variable is missing. \\n\"; exit 1; } \n}\n\n\ni\
f(!$TCDIR) { $TCDIR=\"$HOME/.t_coffee\"; }\n&add_d\
ir ($TCDIR);\n&add_dir ($TCCACHE=\"$TCDIR/cache\")\
;\n&add_dir ($TCTMP=\"$CDIR/tmp\");\n&add_dir ($TC\
M=\"$TCDIR/mcoffee\");\n&add_dir ($TCMETHODS=\"$TC\
DIR/methods\");\n&add_dir ($TCPLUGINS=\"$TCDIR/plu\
gins/$OS\");\n\n\nour $BASE=\"$CD/bin\";\nour $BIN\
=\"$BASE/binaries/$OS\";\nour $DOWNLOAD_DIR=\"$BAS\
E/download\";\nour $DOWNLOAD_FILE=\"$DOWNLOAD_DIR/\
files\";\nour $TMP=\"$BASE/tmp\";\n\n&add_dir($BAS\
E);\n&add_dir($BIN);\n&add_dir($DOWNLOAD_DIR);\n&a\
dd_dir($DOWNLOAD_FILE);\nif (!$DISTRIBUTIONS){$DIS\
TRIBUTIONS=\"$DOWNLOAD_DIR/distributions\";}\n&add\
_dir ($DISTRIBUTIONS);\n&add_dir ($TMP);\n\n\nif  \
  (!$PLUGINS_DIR && !$ROOT_INSTALL){$PLUGINS_DIR=$\
TCPLUGINS;}\nelsif (!$PLUGINS_DIR &&  $ROOT_INSTAL\
L){$PLUGINS_DIR=\"/usr/local/bin/\";}\n\nif    (!$\
INSTALL_DIR && !$ROOT_INSTALL){$INSTALL_DIR=\"$HOM\
E/bin/\";mkpath ($INSTALL_DIR);}\nelsif (!$INSTALL\
_DIR &&  $ROOT_INSTALL){$INSTALL_DIR=\"/usr/local/\
bin/\";}\n\nif (-d \"mcoffee\"){`cp mcoffee/* $TCM\
`;}\n\n\nour $ENV_FILE=\"$TCDIR/t_coffee_env\";\n&\
env_file2putenv ($ENV_FILE);\n&set_proxy($proxy);\\
nmy ($target, $p, $r);\n$target=$p;\n\nforeach $p \
(  ((keys (%PG)),(keys(%MODE)),(@smode)) )\n  {\n \
   if ($ARGV[0] eq $p && $target eq \"\"){$target=\
$p;}\n  }\nif ($target eq \"\"){exit ($EXIT_FAILUR\
E);}\n\n\nforeach $r (@required_applications)\n  {\
\n    my @app_list;\n    my $i;\n    $i=0;\n    \n\
    @app_list=split (/_OR_/, $r);\n    foreach my \
$pg (@app_list)\n      {\n	$i+=&pg_is_installed ($\
pg);\n      }\n    if ($i==0)\n      {\n      prin\
t \"One of the following packages must be installe\
d to proceed: \";\n      foreach my $pg (@app_list\
)\n	{\n	  print (\"$pg \");\n	}\n      die;\n    }\
\n  }\n\n\n\n\n\n\n&sign_license_ni();\n\n\n$PG{C}\
{compiler}=get_C_compiler($CC);\n$PG{Fortran}{comp\
iler}=get_F_compiler($FC);\n$PG{CXX}{compiler}=$PG\
{CPP}{compiler}=$PG{GPP}{compiler}=get_CXX_compile\
r($CXX);\nif ($CXXFLAGS){$PG{CPP}{options}=$PG{GPP\
}{options}=$PG{CXX}{options}=$CXXFLAGS;}\nif ($CFL\
AGS){$PG{C}{options}=$CFLAGS;}\nforeach my $c (key\
s(%PG))\n  {\n    my $arguments;\n    if ($PG{$c}{\
compiler})\n      {\n	$arguments=\"$PG{$c}{compile\
r_flag}=$PG{$c}{compiler} \";\n	if ($PG{$c}{option\
s})\n	  {\n	    $arguments.=\"$PG{$c}{options_flag\
}=$PG{$c}{options} \";\n	  }\n	$PG{$c}{arguments}=\
$arguments;\n      }\n  }\n\nif ($PG{$target}){$PG\
{$target}{install}=1;}\nelse\n  {\n    foreach my \
$pg (keys(%PG))\n      {\n	if ( $target eq \"all\"\
 || ($PG{$pg}{mode}=~/$target/))\n	  {\n	    $PG{$\
pg} {install}=1;\n	  }\n      }\n  }\n\nforeach my\
 $pg (keys(%PG))\n  {\n    if (!$PG{$pg}{update_ac\
tion}){$PG{$pg}{update_action}=$default_update_act\
ion;}\n    elsif ($PG{$pg}{update_action} eq \"nev\
er\"){$PG{$pg}{install}=0;}\n    if ( $force && $P\
G{$pg}{install})\n      {\n	`rm $BIN/$pg $BIN/$pg.\
exe $SILENT`;\n      }\n    if ($PG{$pg}{update_ac\
tion} eq \"update\" && $PG{$pg}{install}){$PG{$pg}\
{update}=1;}\n  }\n\nif (($target=~/clean/))\n  {\\
n    print \"------- cleaning executables -----\\n\
\";\n    `rm bin/* $SILENT`;\n    exit ($EXIT_SUCC\
ESS);\n  }\n\nif ( !$PG{$target}){print \"------- \
Installing T-Coffee Modes\\n\";}\n\nforeach my $m \
(keys(%MODE))\n  {\n    if ( $target eq \"all\" ||\
 $target eq $m)\n      {\n	print \"\\n------- The \
installer will now install the $m components $MODE\
{$m}{description}\\n\";\n	foreach my $pg (keys(%PG\
))\n	  {\n	    if ( $PG{$pg}{mode} =~/$m/ && $PG{$\
pg}{install})\n	      {\n		if ($PG{$pg}{touched}){\
print \"------- $PG{$pg}{dname}: already processed\
\\n\";}\n		else {$PG{$pg}{success}=&install_pg($pg\
);$PG{$pg}{touched}=1;}\n	      }\n	  }\n      }\n\
  }\n\nif ( $PG{$target}){print \"------- Installi\
ng Individual Package\\n\";}\nforeach my $pg (keys\
 (%PG))\n  {\n    \n    if ( $PG{$pg}{install} && \
!$PG{$pg}{touched})\n      {\n	print \"\\n------- \
Install $pg\\n\";\n	$PG{$pg}{success}=&install_pg(\
$pg);$PG{$pg}{touched}=1;\n      }\n  }\nprint \"-\
------ Finishing The installation\\n\";\nmy $final\
_report=&install ($INSTALL_DIR);\n\nprint \"\\n\";\
\nprint \"****************************************\
*****************************\\n\";\nprint \"*****\
***              INSTALLATION SUMMARY          ***\
**************\\n\";\nprint \"********************\
*************************************************\\
\n\";\nprint \"------- SUMMARY package Installatio\
n:\\n\";\nprint \"-------   Executable Installed i\
n: $PLUGINS_DIR\\n\";\n\nforeach my $pg (keys(%PG)\
)\n  {\n    if ( $PG{$pg}{install})\n      {\n	my \
$bin_status=($PG{$pg}{from_binary} && $PG{$pg}{suc\
cess})?\"[from binary]\":\"\";\n	if     ( $PG{$pg}\
{new} && !$PG{$pg}{old})                     {prin\
t \"*------        $PG{$pg}{dname}: installed $bin\
_status\\n\"; $PG{$pg}{status}=1;}\n	elsif  ( $PG{\
$pg}{new} &&  $PG{$pg}{old})                     {\
print \"*------        $PG{$pg}{dname}: updated $b\
in_status\\n\"  ; $PG{$pg}{status}=1;} \n	elsif  (\
!$PG{$pg}{new} &&  $PG{$pg}{old} && !$PG{$pg}{upda\
te}){print \"*------        $PG{$pg}{dname}: previ\
ous\\n\" ; $PG{$pg}{status}=1;}\n	elsif  (!$PG{$pg\
}{new} &&  $PG{$pg}{old} &&  $PG{$pg}{update}){pri\
nt \"*------        $PG{$pg}{dname}: failed update\
 (previous installation available)\\n\";$PG{$pg}{s\
tatus}=0;}\n	else                                 \
                         {print \"*------        $\
PG{$pg}{dname}: failed installation\\n\";$PG{$pg}{\
status}=0;}\n      }\n  }\nmy $failure;\n\nif ( !$\
PG{$target}){print \"*------ SUMMARY mode Installa\
tion:\\n\";}\nforeach my $m (keys(%MODE))\n  {\n  \
\n    if ( $target eq \"all\" || $target eq $m)\n \
     {\n	my $succesful=1;\n	foreach my $pg (keys(%\
PG))\n	  {\n	    if (($PG{$pg}{mode}=~/$m/) && $PG\
{$pg}{install} && $PG{$pg}{status}==0)\n	      {\n\
		$succesful=0;\n		print \"*!!!!!!       $PG{$pg}{\
dname}: Missing\\n\";\n	      }\n	  }\n	if ( $succ\
esful)\n	  {\n	    $MODE{$m}{status}=1;\n	    prin\
t \"*------       MODE $MODE{$m}{dname} SUCCESSFUL\
LY installed\\n\";\n	  }\n	else\n	  {\n	    $failu\
re++;\n	    $MODE{$m}{status}=0;\n	    print \"*!!\
!!!!       MODE $MODE{$m}{dname} UNSUCCESSFULLY in\
stalled\\n\";\n	  }\n      }\n  }\n\n    \n      \\
nif ($clean==1 && ($BASE=~/install4tcoffee/) ){pri\
nt \"*------ Clean Installation Directory: $BASE\\\
n\";`rm -rf $BASE`;}\nforeach my $pg (keys(%PG)){i\
f ($PG{$pg}{install} && $PG{$pg}{status}==0){exit \
($EXIT_FAILURE);}}\n\nif ($failure)\n  {\n    prin\
t \"**********************************************\
***********************\\n\";\n    print \"*******\
*     SOME PACKAGES FAILED TO INSTALL        *****\
************\\n\";\n    print \"******************\
**************************************************\
*\\n\";\n    print \"\\nSome of the reported failu\
res may be due to connectivity problems\";\n    pr\
int \"\\nRerun the installation and the installer \
will specifically try to install the missing packa\
ges\";\n    print \"\\nIf this Fails, go to the or\
iginal website and install the package manually\";\
\n  }\n\nprint \"*********************************\
************************************\\n\";\nprint \
\"********              FINALIZE YOUR INSTALLATION\
    *****************\\n\";\nprint \"*************\
**************************************************\
******\\n\";\nprint \"------- Your executables are\
 in:\\n\"; \nprint \"-------       $PLUGINS_DIR:\\\
n\";\nprint \"------- Add this directory to your p\
ath with the following command:\\n\";\nprint \"---\
----       export PATH=$PLUGINS_DIR:\\$PATH\\n\";\\
nprint \"------- Make this permanent by adding thi\
s line to the file:\\n\";\nprint \"-------       $\
HOME/.bashrc\\n\";\nexit ($EXIT_SUCCESS);  \n  \ns\
ub get_CXX_compiler\n  {\n    my $c=@_[0];\n    my\
 (@clist)=(\"g++\");\n    \n    return get_compil \
($c, @clist);\n }\nsub get_C_compiler\n  {\n    my\
 $c=@_[0];\n    my (@clist)=(\"gcc\", \"cc\", \"ic\
c\");\n    \n    return get_compil ($c, @clist);\n\
 }\n\nsub get_F_compiler\n  {\n    my ($c)=@_[0];\\
n    my @clist=(\"f77\", \"g77\",\"g95\", \"gfortr\
an\", \"ifort\");\n    return get_compil ($c, @cli\
st);\n  } \n       \nsub get_compil\n  {\n    my (\
$fav,@clist)=(@_);\n    \n    #return the first co\
mpiler found installed in the system. Check first \
the favorite\n    foreach my $c ($fav,@clist)\n   \
   {\n	if  (&pg_is_installed ($c)){return $c;}\n  \
    }\n    return \"\";\n  }\nsub exit_if_pg_not_i\
nstalled\n  {\n    my (@arg)=(@_);\n    \n    fore\
ach my $p (@arg)\n      {\n	if ( !&pg_is_installed\
 ($p))\n	  {\n	    print \"!!!!!!!! The $p utility\
 must be installed for this installation to procee\
d [FATAL]\\n\";\n	    die;\n	  }\n      }\n    ret\
urn 1;\n  }\nsub set_proxy\n  {\n    my ($proxy)=(\
@_);\n    my (@list,$p);\n    \n    @list= (\"HTTP\
_proxy\", \"http_proxy\", \"HTTP_PROXY\", \"ALL_pr\
oxy\", \"all_proxy\",\"HTTP_proxy_4_TCOFFEE\",\"ht\
tp_proxy_4_TCOFFEE\");\n    \n    if (!$proxy)\n  \
    {\n	foreach my $p (@list)\n	  {\n	    if ( ($E\
NV_SET{$p}) || $ENV{$p}){$proxy=$ENV{$p};}\n	  }\n\
      }\n    foreach my $p(@list){$ENV{$p}=$proxy;\
}\n  }\n	\nsub check_internet_connection\n  {\n   \
 my $internet;\n    \n    if ( -e \"x\"){unlink (\\
"x\");}\n    if     (&pg_is_installed    (\"wget\"\
)){`wget www.google.com -Ox >/dev/null 2>/dev/null\
`;}\n    elsif  (&pg_is_installed    (\"curl\")){`\
curl www.google.com -ox >/dev/null 2>/dev/null`;}\\
n    else\n      {\n	printf stderr \"\\nERROR: No \
pg for remote file fetching [wget or curl][FATAL]\\
\n\";\n	exit ($EXIT_FAILURE);\n      }\n    \n    \
if ( !-e \"x\" || -s \"x\" < 10){$internet=0;}\n  \
  else {$internet=1;}\n    if (-e \"x\"){unlink \"\
x\";}\n    return $internet;\n  }\nsub url2file\n \
 {\n    my ($cmd, $file,$wget_arg, $curl_arg)=(@_)\
;\n    my ($exit,$flag, $pg, $arg);\n    \n    if \
($INTERNET || check_internet_connection ()){$INTER\
NET=1;}\n    else\n      {\n	print STDERR \"ERROR:\
 No Internet Connection [FATAL:install.pl]\\n\";\n\
	exit ($EXIT_FAILURE);\n      }\n    \n    if     \
(&pg_is_installed    (\"wget\")){$pg=\"wget\"; $fl\
ag=\"-O\";$arg=\"--tries=2 --connect-timeout=10 $w\
get_arg\";}\n    elsif  (&pg_is_installed    (\"cu\
rl\")){$pg=\"curl\"; $flag=\"-o\";$arg=$curl_arg;}\
\n    else\n      {\n	printf stderr \"\\nERROR: No\
 pg for remote file fetching [wget or curl][FATAL]\
\\n\";\n	exit ($EXIT_FAILURE);\n      }\n    \n   \
 \n    if (-e $file){unlink($file);}\n    $exit=sy\
stem \"$pg $cmd $flag$file $arg\";\n    return $ex\
it;\n  }\n\nsub pg_is_installed\n  {\n    my ($p, \
$dir)=(@_);\n    my ($r,$m, $ret);\n    my ($suppo\
rted, $language, $compil);\n    \n  \n    if ( $PG\
{$p})\n      {\n	$language=$PG{$p}{language2};\n	$\
compil=$PG{$language}{compiler};\n      }\n    \n \
   if ( $compil eq \"CPAN\")\n      {\n	if ( syste\
m (\"perl -M$p -e 1\")==$EXIT_SUCCESS){$ret=1;}\n	\
else {$ret=0;}\n      }\n    elsif ($dir)\n      {\
\n	if (-e \"$dir/$p\" || -e \"$dir/$p\\.exe\"){$re\
t=1;}\n	else {$ret=0;}\n      }\n    elsif (-e \"$\
PLUGINS_DIR/$p\" || -e \"$PLUGINS_DIR/$p.exe\"){$r\
et=1;}\n    else\n      {\n	$r=`which $p 2>/dev/nu\
ll`;\n	if ($r eq \"\"){$ret=0;}\n	else {$ret=1;}\n\
      }\n   \n    return $ret;\n  }\nsub install\n\
  {\n    my ($new_bin)=(@_);\n    my ($copied, $re\
port);\n\n    \n    if (!$ROOT_INSTALL)\n      {\n\
	\n	if (-e \"$BIN/t_coffee\"){`$CP $BIN/t_coffee $\
INSTALL_DIR`};\n	`cp $BIN/* $PLUGINS_DIR`;\n	$copi\
ed=1;\n      }\n    else\n      {\n	$copied=&root_\
run (\"You must be root to finalize the installati\
on\", \"$CP $BIN/* $INSTALL_DIR $SILENT\");\n     \
 }\n    \n     \n  if ( !$copied)\n    {\n      $r\
eport=\"*!!!!!! Installation unsuccesful. The exec\
utables have been left in $BASE/bin\\n\";\n    }\n\
  elsif ( $copied && $ROOT)\n    {\n      $report=\
\"*------ Installation succesful. Your executables\
 have been copied in $new_bin and are on your PATH\
\\n\";\n    }\n  elsif ( $copied && !$ROOT)\n    {\
\n      $report= \"*!!!!!! T-Coffee and associated\
 packages have been copied in: $new_bin\\n\";\n   \
   $report.=\"*!!!!!! This address is NOT in your \
PATH sytem variable\\n\";\n      $report.=\"*!!!!!\
! You can do so by adding the following line in yo\
ur ~/.bashrc file:\\n\";\n      $report.=\"*!!!!!!\
 export PATH=$new_bin:\\$PATH\\n\";\n    }\n  retu\
rn $report;\n}\n\nsub sign_license_ni\n  {\n    my\
 $F=new FileHandle;\n    open ($F, \"license.txt\"\
);\n    while (<$F>)\n      {\n	print \"$_\";\n   \
   }\n    close ($F);\n    \n    return;\n  }\n\ns\
ub install_pg\n  {\n    my ($pg)=(@_);\n    my ($r\
eport, $previous, $language, $compiler, $return);\\
n    \n    if (!$PG{$pg}{install}){return 1;}\n   \
 \n    $previous=&pg_is_installed ($pg);\n    \n  \
  if ($PG{$pg}{update_action} eq \"no_update\" && \
$previous)\n      {\n	$PG{$pg}{old}=1;\n	$PG{$pg}{\
new}=0;\n	$return=1;\n      }\n    else\n      {\n\
	$PG{$pg}{old}=$previous;\n	\n	if ($PG{$pg} {langu\
age2} eq \"Perl\"){&install_perl_package ($pg);}\n\
	elsif ($BINARIES_ONLY && &install_binary_package \
($pg)){$PG{$pg}{from_binary}=1;}\n	elsif (&install\
_source_package ($pg)){;}\n	else \n	  {\n	    \n	 \
   if (!&supported_os($OS))\n	      {\n		print \"!\
!!!!!!! $pg compilation failed, binary unsupported\
 for $OS\\n\"; \n	      }\n	    elsif (!($PG{$pg}{\
from_binary}=&install_binary_package ($pg)))\n	   \
   {\n		print \"!!!!!!!! $pg compilation and  bina\
ry installation failed\\n\";\n	      }\n	  }\n	$PG\
{$pg}{new}=$return=&pg_is_installed ($pg,$BIN);\n \
     }\n\n    \n    return $return;\n  }\nsub inst\
all_perl_package\n  {\n    my ($pg)=(@_);\n    my \
($report, $language, $compiler);\n    \n    $langu\
age=$PG{$pg} {language2};\n    $compiler=$PG{$lang\
uage}{compiler};\n    \n    if (!&pg_is_installed \
($pg))\n      {\n	if ( $OS eq \"windows\"){`perl -\
M$compiler -e 'install $pg'`;}\n	elsif ( $ROOT eq \
\"sudo\"){system (\"sudo perl -M$compiler -e 'inst\
all $pg'\");}\n	else {system (\"su root -c perl -M\
$compiler -e 'install $pg'\");}\n      }\n    retu\
rn &pg_is_installed ($pg);\n  }\n\n\n\nsub install\
_source_package\n  {\n    my ($pg)=(@_);\n    my (\
$report, $download, $arguments, $language, $addres\
s, $name, $ext, $main_dir, $distrib);\n    my $wge\
t_tmp=\"$TMP/wget.tmp\";\n    my (@fl);\n    if ( \
-e \"$BIN/$pg\" || -e \"$BIN/$pg.exe\"){return 1;}\
\n    \n    #\n    # check if the module exists in\
 the repository cache \n    #\n	if( repo_load($pg)\
 ) {\n		return 1;\n	}\n    \n    if ($pg eq \"t_co\
ffee\")  {return   &install_t_coffee ($pg);}\n    \
elsif ($pg eq \"TMalign\"){return   &install_TMali\
gn ($pg);}\n    \n    chdir $DISTRIBUTIONS;\n    \\
n    $download=$PG{$pg}{source};\n    \n    if (($\
download =~/tgz/))\n      {\n	($address,$name,$ext\
)=($download=~/(.+\\/)([^\\/]+)(\\.tgz).*/);\n    \
  }\n    elsif (($download=~/tar\\.gz/))\n      {\\
n	($address,$name,$ext)=($download=~/(.+\\/)([^\\/\
]+)(\\.tar\\.gz).*/);\n      }\n    elsif (($downl\
oad=~/tar/))\n      {\n	($address,$name,$ext)=($do\
wnload=~/(.+\\/)([^\\/]+)(\\.tar).*/);\n      }\n \
   else\n      {\n	($address,$name)=($download=~/(\
.+\\/)([^\\/]+)/);\n	$ext=\"\";\n      }\n    $dis\
trib=\"$name$ext\";\n    \n    if ( !-d $pg){mkdir\
 $pg;}\n    chdir $pg;\n   \n    #get the distribu\
tion if available\n    if ( -e \"$DOWNLOAD_DIR/$di\
strib\")\n      {\n	`$CP $DOWNLOAD_DIR/$distrib .`\
;\n      }\n    #UNTAR and Prepare everything\n   \
 if (!-e \"$name.tar\" && !-e \"$name\")\n      {\\
n	&check_rm ($wget_tmp);\n	print \"\\n------- Down\
loading/Installing $pg\\n\";\n	\n	if (!-e $distrib\
 && &url2file (\"$download\", \"$wget_tmp\")==$EXI\
T_SUCCESS)\n	  {\n	    \n	    `mv $wget_tmp $distr\
ib`;\n	    `$CP $distrib $DOWNLOAD_DIR/`;\n	  }\n\\
n	if (!-e $distrib)\n	  {\n	    print \"!!!!!!! Do\
wnload of $pg distribution failed\\n\";\n	    prin\
t \"!!!!!!! Check Address: $PG{$pg}{source}\\n\";\\
n	    return 0;\n	  }\n	print \"\\n------- unzippi\
ng/untaring $name\\n\";\n	if (($ext =~/z/))\n	  { \
\n	    &flush_command (\"gunzip $name$ext\");\n	  \
  \n	  }\n	if (($ext =~/tar/) || ($ext =~/tgz/))\n\
	  {\n	    &flush_command(\"tar -xvf $name.tar\");\
\n	  }\n      }\n    #Guess and enter the distribu\
tion directory\n    @fl=ls($p);\n    foreach my $f\
 (@fl)\n      {\n	if (-d $f)\n	  {\n	    $main_dir\
=$f;\n	  }\n      }\n    if (-d $main_dir)\n	  \n \
     {\n	chdir $main_dir;}\n    else\n      {\n	pr\
int \"Error: $main_dir does not exist\";\n      }\\
n    print \"\\n------- Compiling/Installing $pg\\\
n\";\n    `make clean $SILENT`;\n    \n    \n    #\
\n    # SAP module\n    #\n    if ($pg eq \"sap\")\
\n      {\n	if (-e \"./configure\")\n	  {\n	    #n\
ew sap distribution\n	    if ($OS eq \"macosx\")\n\
	      {\n		&replace_line_in_file (\"./src/galloc.\
h\", \"malloc.h\",  \"\");\n		&replace_line_in_fil\
e (\"./src/pdbprot.h\", \"malloc.h\", \"\");\n		&r\
eplace_line_in_file (\"./src/pdbprot.c\", \"malloc\
.h\", \"\");\n	      }\n	    \n	    &flush_command\
 (\"./configure\");\n	    &flush_command (\"make c\
lean\");\n	    &flush_command (\"make\");\n	    &c\
heck_cp (\"./src/$pg\", \"$BIN\");\n	    repo_stor\
e(\"./src/$pg\");\n	  }\n	else\n	  {\n	    #old st\
yle distribution\n	    `rm *.o sap  sap.exe ./util\
/aa/*.o  ./util/wt/.o $SILENT`;\n	    &flush_comma\
nd (\"make $arguments sap\");\n	    &check_cp ($pg\
, \"$BIN\");\n	    repo_store($pg);\n	  }\n      }\
\n    \n    #\n    # CLUSTALW2 module\n    #\n    \
elsif ($pg eq \"clustalw2\")\n      {\n	&flush_com\
mand(\"./configure\");\n	&flush_command(\"make $ar\
guments\");\n	&check_cp (\"./src/$pg\", \"$BIN\");\
\n	repo_store(\"./src/$pg\");\n      }\n    \n    \
#\n    # FSA module\n    # \n    elsif ($pg eq \"f\
sa\")\n      {\n	&flush_command(\"./configure --pr\
efix=$BIN\");\n	&flush_command(\"make $arguments\"\
);\n	&flush_command (\"make install\");\n\n	repo_s\
tore(\"fsa\", \"$BIN/bin\");\n	`mv $BIN/bin/* $BIN\
`;\n	`rmdir $BIN/bin`;\n      }\n    \n    #\n    \
# CLUSTALW module\n    #\n    elsif ($pg eq \"clus\
talw\")\n      {\n	&flush_command(\"make $argument\
s clustalw\");\n	`$CP $pg $BIN $SILENT`;\n	repo_st\
ore($pg);\n      }\n    \n    #\n    # MAFFT modul\
e\n    #\n    elsif ($pg eq \"mafft\")\n      {\n	\
my $base=cwd();\n	my $c;\n	\n	#compile core\n	mkpa\
th (\"./mafft/bin\");\n	mkpath (\"./mafft/lib\");\\
n	chdir \"$base/core\";\n	`make clean $SILENT`;\n	\
&flush_command (\"make $arguments\");\n	&flush_com\
mand (\"make install LIBDIR=../mafft/lib BINDIR=..\
/mafft/bin\");\n	\n	#compile extension\n	chdir \"$\
base/extensions\";\n	`make clean $SILENT`;\n	&flus\
h_command (\"make $arguments\");\n	&flush_command \
(\"make install LIBDIR=../mafft/lib BINDIR=../maff\
t/bin\");\n	\n	#put everything in mafft and copy t\
he compiled stuff in bin\n	chdir \"$base\";\n	if (\
$ROOT_INSTALL)\n	  {\n	    &root_run (\"You Must b\
e Root to Install MAFFT\\n\", \"mkdir /usr/local/m\
afft/;$CP mafft/lib/* /usr/local/mafft;$CP mafft/l\
ib/mafft* /usr/local/bin ;$CP mafft/bin/mafft /usr\
/local/bin/; \");\n	  }\n	else\n	  {\n	    `$CP ma\
fft/lib/*  $BIN`;\n	    `$CP mafft/bin/mafft  $BIN\
`;\n	  }\n	`tar -cvf mafft.tar mafft`;\n	`gzip maf\
ft.tar`;\n	`mv mafft.tar.gz $BIN`;\n	\n	repo_store\
(\"mafft/bin/mafft\", \"mafft/lib/\", \"$BIN/mafft\
.tar.gz\");\n      }\n      \n    #\n    # DIALIGN\
-TX module\n    #\n    elsif ( $pg eq \"dialign-tx\
\" )\n      {\n	my $f;\n	my $base=cwd();\n\n	chdir\
 \"./source\";\n	if ($OS eq \"macosx\"){&flush_com\
mand (\"cp makefile.MAC_OS makefile\");}\n\n	&flus\
h_command (\" make CPPFLAGS='-O3 -funroll-loops' a\
ll\");\n	\n	chdir \"..\";\n	&check_cp (\"./source/\
$pg\", \"$BIN\");\n	repo_store(\"./source/$pg\");\\
n      }\n      \n    #\n    # DIALIGN-T module \n\
    # (is the same as dialign-tx, but it is mantai\
ned for backward name compatibility with tcoffee)\\
n    #\n    elsif ( $pg eq \"dialign-t\" )\n      \
{\n	my $f;\n	my $base=cwd();\n\n	chdir \"./source\\
";\n	if ($OS eq \"macosx\"){&flush_command (\"cp m\
akefile.MAC_OS makefile\");}\n\n	&flush_command (\\
" make CPPFLAGS='-O3 -funroll-loops' all\");\n	\n	\
chdir \"..\";\n	&check_cp (\"./source/dialign-tx\"\
, \"$BIN/dialign-t\");\n	repo_store(\"$BIN/dialign\
-t\");	\n      }      \n      \n    #\n    # POA m\
odule\n    #\n    elsif ($pg eq \"poa\")\n      {\\
n	&flush_command (\"make $arguments poa\");\n	&che\
ck_cp (\"$pg\", \"$BIN\");\n	repo_store(\"$pg\");\\
n      }\n     \n     \n    #\n    # PROBCONS modu\
le\n    #\n    elsif ( $pg eq \"probcons\")\n     \
 {\n	&add_C_libraries(\"./ProbabilisticModel.h\", \
\"list\", \"cstring\");\n	\n	`rm *.exe $SILENT`;\n\
	&flush_command (\"make $arguments probcons\");\n	\
&check_cp(\"$pg\", \"$BIN/$pg\");\n	repo_store(\"$\
pg\");\n      }\n      \n    #\n    # PROBCONS RNA\
 module\n    #\n    elsif ( $pg eq \"probconsRNA\"\
)\n      {\n	&add_C_libraries(\"./ProbabilisticMod\
el.h\", \"list\", \"cstring\");\n	&add_C_libraries\
(\"./Main.cc\", \"iomanip\", \"cstring\",\"climits\
\");\n	`rm *.exe $SILENT`;\n	&flush_command (\"mak\
e $arguments probcons\");\n	&check_cp(\"probcons\"\
, \"$BIN/$pg\");\n	repo_store(\"$BIN/$pg\");\n    \
  }\n\n	#\n	# MUSCLE module\n	#\n    elsif (  $pg \
eq \"muscle\")\n      {	\n	`rm *.o muscle muscle.e\
xe $SILENT`;\n	if ($OS eq \"macosx\" || $OS eq \"l\
inux\")\n	  {\n	    &replace_line_in_file (\"./Mak\
efile\", \"LDLIBS = -lm -static\",  \"LDLIBS = -lm\
\");\n	  }\n	elsif ($OS eq \"windows\")\n	  {\n	  \
  &replace_line_in_file (\"./intmath.cpp\",  \"dou\
ble log2e\",      \"double cedric_log\");\n	    &r\
eplace_line_in_file (\"./intmath.cpp\",  \"double \
log2\",       \"double log_notuse\");\n	    &repla\
ce_line_in_file (\"./intmath.cpp\",  \"double cedr\
ic_log\", \"double log2e\");\n	  }\n	&flush_comman\
d (\"make $arguments all\");\n	&check_cp(\"$pg\", \
\"$BIN\");\n	repo_store(\"$pg\");	\n      }\n     \
 \n     #\n     # MUS4 module\n     #\n     elsif \
(  $pg eq \"mus4\")\n      {\n	`rm *.o muscle musc\
le.exe $SILENT`;\n	&flush_command (\"./mk\");\n	&c\
heck_cp(\"$pg\", \"$BIN\");\n	repo_store(\"$pg\");\
	\n      }\n      \n    #\n    # PCMA module\n    \
#\n    elsif ( $pg eq \"pcma\")\n      {\n	if ($OS\
 eq \"macosx\")\n	  {\n	    &replace_line_in_file \
(\"./alcomp2.c\", \"malloc.h\",  \"\");\n	  }\n	&f\
lush_command (\"make $arguments pcma\");\n	&check_\
cp(\"$pg\", \"$BIN\");\n	repo_store(\"$pg\");	\n  \
    }\n      \n    #\n    # KALIGN module\n    #\n\
    elsif ($pg eq \"kalign\")\n      {\n	&flush_co\
mmand (\"./configure\");\n	&flush_command(\"make $\
arguments\");\n	&check_cp (\"$pg\",$BIN);\n	repo_s\
tore(\"$pg\");	\n      }\n      \n    #\n    # AMA\
P module\n    #\n    elsif ( $pg eq \"amap\")\n   \
   {\n	&add_C_libraries(\"./Amap.cc\", \"iomanip\"\
, \"cstring\",\"climits\");	\n	`make clean $SILENT\
`;\n	&flush_command (\"make $arguments all\");\n	&\
check_cp (\"$pg\", $BIN);\n	repo_store(\"$pg\");	\\
n      }\n      \n    #\n    # PRODA module\n    #\
\n    elsif ( $pg eq \"proda\")\n      {\n	&add_C_\
libraries(\"AlignedFragment.h\", \"vector\", \"ios\
tream\", \"cstring\",\"cstdlib\");\n	&add_C_librar\
ies(\"Main.cc\", \"vector\", \"climits\");	\n	&add\
_C_libraries(\"Sequence.cc\", \"stdlib.h\", \"cstd\
io\");	\n	&flush_command (\"make $arguments all\")\
;\n	&check_cp (\"$pg\", $BIN);\n	repo_store(\"$pg\\
");	\n      }\n      \n    #\n    # PRANK module\n\
    #\n    elsif ( $pg eq \"prank\")\n      {\n	&f\
lush_command (\"make $arguments all\");\n	&check_c\
p (\"$pg\", $BIN);\n	repo_store(\"$pg\");	\n      \
}\n      \n    #\n    # !!!! MUSTANG module\n    #\
\n     elsif ( $pg eq \"mustang\")\n      {\n	&flu\
sh_command (\"rm ./bin/*\");\n	&flush_command (\"m\
ake $arguments all\");\n\n	if ( $OS=~/windows/){&f\
lush_command(\"cp ./bin/* $BIN/mustang.exe\");}\n	\
else {&flush_command(\"cp ./bin/* $BIN/mustang\");\
}\n	\n	repo_store(\"$BIN/mustang\");\n      }\n\n	\
#\n	# RNAplfold module\n	#\n    elsif ( $pg eq \"R\
NAplfold\")\n      {\n	&flush_command(\"./configur\
e\");\n	&flush_command (\"make $arguments all\");\\
n	&check_cp(\"./Progs/RNAplfold\", \"$BIN\");\n	&c\
heck_cp(\"./Progs/RNAalifold\", \"$BIN\");\n	&chec\
k_cp(\"./Progs/RNAfold\", \"$BIN\");\n	\n	repo_sto\
re(\"./Progs/RNAplfold\", \"./Progs/RNAalifold\", \
\"./Progs/RNAfold\");\n      }\n      \n    #\n   \
 # !!! RETREE module\n    #\n    elsif ( $pg eq \"\
retree\")\n      {\n	chdir \"src\";\n	&flush_comma\
nd (\"make $arguments all\");\n	&flush_command (\"\
make put\");\n	system \"cp ../exe/* $BIN\";\n	\n	r\
epo_store(\"retree\", \"../exe\");\n      }\n	\n  \
  chdir $CDIR;\n    return &pg_is_installed ($pg, \
$BIN);\n  }\n\nsub install_t_coffee\n  {\n    my (\
$pg)=(@_);\n    my ($report,$cflags, $arguments, $\
language, $compiler) ;\n    #1-Install T-Coffee\n \
   chdir \"t_coffee_source\";\n    &flush_command \
(\"make clean\");\n    print \"\\n------- Compilin\
g T-Coffee\\n\";\n    $language=$PG{$pg} {language\
2};\n    $arguments=$PG{$language}{arguments};\n  \
  if (!($arguments =~/CFLAGS/)){$arguments .= \" C\
FLAGS=-O2 \";}\n\n    if ( $CC ne \"\"){&flush_com\
mand (\"make -i $arguments t_coffee\");}\n    &che\
ck_cp ($pg, $BIN);\n    \n    chdir $CDIR;\n    re\
turn &pg_is_installed ($pg, $BIN);\n  }\nsub insta\
ll_TMalign\n  {\n    my ($pg)=(@_);\n    my $repor\
t;\n    chdir \"t_coffee_source\";\n    print \"\\\
n------- Compiling TMalign\\n\";\n    `rm TMalign \
TMalign.exe $SILENT`;\n    if ( $FC ne \"\"){&flus\
h_command (\"make -i $PG{Fortran}{arguments} TMali\
gn\");}\n    &check_cp ($pg, $BIN);\n    repo_stor\
e($pg);\n\n    if ( !-e \"$BIN/$pg\" && pg_has_bin\
ary_distrib ($pg))\n      {\n	print \"!!!!!!! Comp\
ilation of $pg impossible. Will try to install fro\
m binary\\n\";\n	return &install_binary_package ($\
pg);\n      }\n    chdir $CDIR;\n    return &pg_is\
_installed ($pg, $BIN);\n  }\n\nsub pg_has_binary_\
distrib\n  {\n    my ($pg)=(@_);\n    if ($PG{$pg}\
{windows}){return 1;}\n    elsif ($PG{$pg}{osx}){r\
eturn 1;}\n    elsif ($PG{$pg}{linux}){return 1;}\\
n    return 0;\n  }\nsub install_binary_package\n \
 {\n    my ($pg)=(@_);\n    my ($base,$report,$nam\
e, $download, $arguments, $language, $dir);\n    m\
y $isdir;\n    &input_os();\n    \n    if (!&suppo\
rted_os($OS)){return 0;}\n    if ( $PG{$pg}{binary\
}){$name=$PG{$pg}{binary};}\n    else \n      {\n	\
$name=$pg;\n	if ( $OS eq \"windows\"){$name.=\".ex\
e\";}\n      }\n    \n    $download=\"$WEB_BASE/Pa\
ckages/Binaries/$OS/$name\";\n    \n    $base=cwd(\
);\n    chdir $TMP;\n    \n    if (!-e $name)\n   \
   {\n	`rm x $SILENT`;\n	if ( url2file(\"$download\
\",\"x\")==$EXIT_SUCCESS)\n	  {\n	    `mv x $name`\
;\n	  }\n      }\n    \n    if (!-e $name)\n      \
{\n	print \"!!!!!!! $PG{$pg}{dname}: Download of $\
pg binary failed\\n\";\n	print \"!!!!!!! $PG{$pg}{\
dname}: Check Address: $download\\n\";\n	return 0;\
\n      }\n    print \"\\n------- Installing $pg\\\
n\";\n    \n    if ($name =~/tar\\.gz/)\n      {\n\
	`gunzip  $name`;\n	`tar -xvf $pg.tar`;\n	chdir $p\
g;\n	if ( $pg eq \"mafft\")\n	  {\n	    if ($ROOT_\
INSTALL)\n	      {\n		&root_run (\"You Must be Roo\
r to Install MAFFT\\n\", \"$CP mafft/bin/* /usr/lo\
cal/mafft;mkdir /usr/local/mafft/; $CP mafft/lib/*\
 /usr/local/bin/\");\n	      }\n	    else\n	      \
{\n		`$CP $TMP/$pg/bin/* $BIN $SILENT`;\n		`$CP $T\
MP/$pg/lib/* $BIN $SILENT`;\n	      }\n	  }\n	else\
\n	  {\n	    if (-e \"$TMP/$pg/data\"){`$CP $TMP/$\
pg/data/* $TCM $SILENT`;}\n	    if (!($pg=~/\\*/))\
{`rm -rf $pg`;}\n	  }\n      }\n    else\n      {\\
n	&check_cp (\"$pg\", \"$BIN\");\n	`chmod u+x $BIN\
/$pg`; \n	unlink ($pg);\n      }\n    chdir $base;\
\n    $PG{$pg}{from_binary}=1;\n    return &pg_is_\
installed ($pg, $BIN);\n  }\n\nsub add_dir \n  {\n\
    my $dir=@_[0];\n    \n    if (!-e $dir && !-d \
$dir)\n      {\n	my @l;\n	umask (0000);\n	@l=mkpat\
h ($dir,{mode => 0777});\n	\n      }\n    else\n  \
    {\n	return 0;\n      }\n  }\nsub check_rm \n  \
{\n    my ($file)=(@_);\n    \n    if ( -e $file)\\
n      {\n	return unlink($file);\n      }\n    ret\
urn 0;\n  }\nsub check_cp\n  {\n    my ($from, $to\
)=(@_);\n    if ( !-e $from && -e \"$from\\.exe\")\
{$from=\"$from\\.exe\";}\n    if ( !-e $from){retu\
rn 0;}\n        \n    `$CP $from $to`;\n    return\
 1;\n  }\n\nsub repo_store \n{\n   # check that al\
l required data are available\n   if( $REPO_ROOT e\
q \"\" ) { return; }\n\n\n    # extract the packag\
e name from the specified path\n    my $pg =`basen\
ame $_[0]`;\n    chomp($pg);\n	\n    my $VER = $PG\
{$pg}{version};\n    my $CACHE = \"$REPO_ROOT/$pg/\
$VER/$OSNAME-$OSARCH\"; \n    \n    print \"------\
-- Storing package: \\\"$pg\\\" to path: $CACHE\\n\
\";\n    \n    # clean the cache path if exists an\
d create it again\n    `rm -rf $CACHE`;\n    `mkdi\
r -p $CACHE`;\n    \n 	for my $path (@_) {\n\n	   \
 # check if it is a single file \n	 	if( -f $path \
) {\n	    	`cp $path $CACHE`;\n		}\n		# .. or a di\
rectory, in this case copy all the content \n		els\
if( -d $path ) {\n			opendir(IMD, $path);\n			my @\
thefiles= readdir(IMD);\n			closedir(IMD);\n			\n	\
		for my $_file (@thefiles) {\n				if( $_file ne \\
".\" && $_file ne \"..\") {\n	    			`cp $path/$_f\
ile $CACHE`;\n				}\n			}\n		} \n	}	   \n    \n	\n\
}   \n\nsub repo_load \n{\n    my ($pg)=(@_);\n\n \
   # check that all required data are available\n \
   if( $REPO_ROOT eq \"\" ) { return 0; }\n\n    m\
y $VER = $PG{$pg}{version};\n    my $CACHE = \"$RE\
PO_ROOT/$pg/$VER/$OSNAME-$OSARCH\"; \n    if( !-e \
\"$CACHE/$pg\" ) {\n   	 	print \"-------- Module \
\\\"$pg\\\" NOT found on repository cache.\\n\";\n\
    	return 0;\n    }\n    \n    print \"-------- \
Module \\\"$pg\\\" found on repository cache. Usin\
g copy on path: $CACHE\\n\";\n    `cp $CACHE/* $BI\
N`;\n    return 1;\n}\n\nsub check_file_list_exist\
s \n  {\n    my ($base, @flist)=(@_);\n    my $f;\\
n\n    foreach $f (@flist)\n      {\n	if ( !-e \"$\
base/$f\"){return 0;}\n      }\n    return 1;\n  }\
\nsub ls\n  {\n    my $f=@_[0];\n    my @fl;\n    \
chomp(@fl=`ls -1 $f`);\n    return @fl;\n  }\nsub \
flush_command\n  {\n    my $command=@_[0];\n    my\
 $F=new FileHandle;\n    open ($F, \"$command|\");\
\n    while (<$F>){print \"    --- $_\";}\n    clo\
se ($F);\n  }    \n\nsub input_installation_direct\
ory\n  {\n    my $dir=@_[0];\n    my $new;\n    \n\
    print \"------- The current installation direc\
tory is: [$dir]\\n\";\n    print \"??????? Return \
to keep the default or new value:\";\n   \n    if \
($NO_QUESTION==0)\n      {\n	chomp ($new=<stdin>);\
\n	while ( $new ne \"\" && !input_yes (\"You have \
entered $new. Is this correct? ([y]/n):\"))\n	  {\\
n	    print \"???????New installation directory:\"\
;\n	    chomp ($new=<stdin>);\n	  }\n	$dir=($new e\
q \"\")?$dir:$new;\n	$dir=~s/\\/$//;\n      }\n   \
 \n    if ( -d $dir){return $dir;}\n    elsif (&ro\
ot_run (\"You must be root to create $dir\",\"mkdi\
r $dir\")==$EXIT_SUCCESS){return $dir;}\n    else\\
n      {\n	print \"!!!!!!! $dir could not be creat\
ed\\n\";\n	if ( $NO_QUESTION)\n	  {\n	    return \\
"\";\n	  }\n	elsif ( &input_yes (\"??????? Do you \
want to provide a new directory([y]/n)?:\"))\n	  {\
\n	    return input_installation_directory ($dir);\
\n	  }\n	else\n	  {\n	    return \"\";\n	  }\n    \
  }\n    \n  }\nsub input_yes\n  {\n    my $questi\
on =@_[0];\n    my $answer;\n\n    if ($NO_QUESTIO\
N==1){return 1;}\n    \n    if ($question eq \"\")\
{$question=\"??????? Do you wish to proceed ([y]/n\
)?:\";}\n    print $question;\n    chomp($answer=l\
c(<STDIN>));\n    if (($answer=~/^y/) || $answer e\
q \"\"){return 1;}\n    elsif ( ($answer=~/^n/)){r\
eturn 0;}\n    else\n      {\n	return input_yes($q\
uestion);\n      }\n  }\nsub root_run\n  {\n    my\
 ($txt, $cmd)=(@_);\n    \n    if ( system ($cmd)=\
=$EXIT_SUCCESS){return $EXIT_SUCCESS;}\n    else \\
n      {\n	print \"------- $txt\\n\";\n	if ( $ROOT\
 eq \"sudo\"){return system (\"sudo $cmd\");}\n	el\
se {return system (\"su root -c \\\"$cmd\\\"\");}\\
n      }\n  }\nsub get_root\n  {\n    if (&pg_is_i\
nstalled (\"sudo\")){return \"sudo\";}\n    else {\
return \"su\";}\n  }\n\nsub get_os\n  {\n    my $r\
aw_os=`uname`;\n    my $os;\n\n    $raw_os=lc ($ra\
w_os);\n    \n    if ($raw_os =~/cygwin/){$os=\"wi\
ndows\";}\n    elsif ($raw_os =~/linux/){$os=\"lin\
ux\";}\n    elsif ($raw_os =~/osx/){$os=\"macosx\"\
;}\n    elsif ($raw_os =~/darwin/){$os=\"macosx\";\
}\n    else\n      {\n	$os=$raw_os;\n      }\n    \
return $os;\n  }\nsub input_os\n  {\n    my $answe\
r;\n    if ($OS) {return $OS;}\n    \n    print \"\
??????? which os do you use: [w]indows, [l]inux, [\
m]acosx:?\";\n    $answer=lc(<STDIN>);\n\n    if (\
($answer=~/^m/)){$OS=\"macosx\";}\n    elsif ( ($a\
nswer=~/^w/)){$OS=\"windows\";}\n    elsif ( ($ans\
wer=~/^linux/)){$OS=\"linux\";}\n    \n    else\n \
     {\n	return &input_os();\n      }\n    return \
$OS;\n  }\n\nsub supported_os\n  {\n    my ($os)=(\
@_[0]);\n    return $SUPPORTED_OS{$os};\n  }\n    \
\n    \n\n\nsub update_tclinkdb \n  {\n    my $fil\
e =@_[0];\n    my $name;\n    my $F=new FileHandle\
;\n    my ($download, $address, $name, $l, $db);\n\
    \n    if ( $file eq \"update\"){$file=$TCLINKD\
B_ADDRESS;}\n    \n    if ( $file =~/http:\\/\\// \
|| $file =~/ftp:\\/\\//)\n      {\n	($address, $na\
me)=($download=~/(.*)\\/([^\\/]+)$/);\n	`rm x $SIL\
ENT`;\n	if (&url2file ($file,\"x\")==$EXIT_SUCCESS\
)\n	  {\n	    print \"------- Susscessful upload o\
f $name\";\n	    `mv x $name`;\n	    $file=$name;\\
n	  }\n      }\n    open ($F, \"$file\");\n    whi\
le (<$F>)\n      {\n	my $l=$_;\n	if (($l =~/^\\/\\\
//) || ($db=~/^#/)){;}\n	elsif ( !($l =~/\\w/)){;}\
\n	else\n	  {\n	    my @v=split (/\\s+/, $l);\n	  \
  if ( $l=~/^MODE/)\n	      {\n		$MODE{$v[1]}{$v[2\
]}=$v[3];\n	      }\n	    elsif ($l=~/^PG/)\n	    \
  {\n		$PG{$v[1]}{$v[2]}=$v[3];\n	      }\n	  }\n \
     }\n    close ($F);\n    &post_process_PG();\n\
    return;\n  }\n\n\n\nsub initialize_PG\n  {\n\n\
$PG{\"t_coffee\"}{\"4_TCOFFEE\"}=\"TCOFFEE\";\n$PG\
{\"t_coffee\"}{\"type\"}=\"sequence_multiple_align\
er\";\n$PG{\"t_coffee\"}{\"ADDRESS\"}=\"http://www\
.tcoffee.org\";\n$PG{\"t_coffee\"}{\"language\"}=\\
"C\";\n$PG{\"t_coffee\"}{\"language2\"}=\"C\";\n$P\
G{\"t_coffee\"}{\"source\"}=\"http://www.tcoffee.o\
rg/Packages/T-COFFEE_distribution.tar.gz\";\n$PG{\\
"t_coffee\"}{\"update_action\"}=\"always\";\n$PG{\\
"t_coffee\"}{\"mode\"}=\"tcoffee,mcoffee,rcoffee,e\
xpresso,3dcoffee\";\n$PG{\"clustalw2\"}{\"4_TCOFFE\
E\"}=\"CLUSTALW2\";\n$PG{\"clustalw2\"}{\"type\"}=\
\"sequence_multiple_aligner\";\n$PG{\"clustalw2\"}\
{\"ADDRESS\"}=\"http://www.clustal.org\";\n$PG{\"c\
lustalw2\"}{\"language\"}=\"C++\";\n$PG{\"clustalw\
2\"}{\"language2\"}=\"CXX\";\n$PG{\"clustalw2\"}{\\
"source\"}=\"http://www.clustal.org/download/2.0.1\
0/clustalw-2.0.10-src.tar.gz\";\n$PG{\"clustalw2\"\
}{\"mode\"}=\"mcoffee,rcoffee\";\n$PG{\"clustalw2\\
"}{\"version\"}=\"2.0.10\";\n$PG{\"clustalw\"}{\"4\
_TCOFFEE\"}=\"CLUSTALW\";\n$PG{\"clustalw\"}{\"typ\
e\"}=\"sequence_multiple_aligner\";\n$PG{\"clustal\
w\"}{\"ADDRESS\"}=\"http://www.clustal.org\";\n$PG\
{\"clustalw\"}{\"language\"}=\"C\";\n$PG{\"clustal\
w\"}{\"language2\"}=\"C\";\n$PG{\"clustalw\"}{\"so\
urce\"}=\"http://www.clustal.org/download/1.X/ftp-\
igbmc.u-strasbg.fr/pub/ClustalW/clustalw1.82.UNIX.\
tar.gz\";\n$PG{\"clustalw\"}{\"mode\"}=\"mcoffee,r\
coffee\";\n$PG{\"clustalw\"}{\"version\"}=\"1.82\"\
;\n$PG{\"dialign-t\"}{\"4_TCOFFEE\"}=\"DIALIGNT\";\
\n$PG{\"dialign-t\"}{\"type\"}=\"sequence_multiple\
_aligner\";\n$PG{\"dialign-t\"}{\"ADDRESS\"}=\"htt\
p://dialign-tx.gobics.de/\";\n$PG{\"dialign-t\"}{\\
"DIR\"}=\"/usr/share/dialign-tx/\";\n$PG{\"dialign\
-t\"}{\"language\"}=\"C\";\n$PG{\"dialign-t\"}{\"l\
anguage2\"}=\"C\";\n$PG{\"dialign-t\"}{\"source\"}\
=\"http://dialign-tx.gobics.de/DIALIGN-TX_1.0.2.ta\
r.gz\";\n$PG{\"dialign-t\"}{\"mode\"}=\"mcoffee\";\
\n$PG{\"dialign-t\"}{\"binary\"}=\"dialign-t\";\n$\
PG{\"dialign-t\"}{\"version\"}=\"1.0.2\";\n$PG{\"d\
ialign-tx\"}{\"4_TCOFFEE\"}=\"DIALIGNTX\";\n$PG{\"\
dialign-tx\"}{\"type\"}=\"sequence_multiple_aligne\
r\";\n$PG{\"dialign-tx\"}{\"ADDRESS\"}=\"http://di\
align-tx.gobics.de/\";\n$PG{\"dialign-tx\"}{\"DIR\\
"}=\"/usr/share/dialign-tx/\";\n$PG{\"dialign-tx\"\
}{\"language\"}=\"C\";\n$PG{\"dialign-tx\"}{\"lang\
uage2\"}=\"C\";\n$PG{\"dialign-tx\"}{\"source\"}=\\
"http://dialign-tx.gobics.de/DIALIGN-TX_1.0.2.tar.\
gz\";\n$PG{\"dialign-tx\"}{\"mode\"}=\"mcoffee\";\\
n$PG{\"dialign-tx\"}{\"binary\"}=\"dialign-tx\";\n\
$PG{\"dialign-tx\"}{\"version\"}=\"1.0.2\";\n$PG{\\
"poa\"}{\"4_TCOFFEE\"}=\"POA\";\n$PG{\"poa\"}{\"ty\
pe\"}=\"sequence_multiple_aligner\";\n$PG{\"poa\"}\
{\"ADDRESS\"}=\"http://www.bioinformatics.ucla.edu\
/poa/\";\n$PG{\"poa\"}{\"language\"}=\"C\";\n$PG{\\
"poa\"}{\"language2\"}=\"C\";\n$PG{\"poa\"}{\"sour\
ce\"}=\"http://downloads.sourceforge.net/poamsa/po\
aV2.tar.gz\";\n$PG{\"poa\"}{\"DIR\"}=\"/usr/share/\
\";\n$PG{\"poa\"}{\"FILE1\"}=\"blosum80.mat\";\n$P\
G{\"poa\"}{\"mode\"}=\"mcoffee\";\n$PG{\"poa\"}{\"\
binary\"}=\"poa\";\n$PG{\"poa\"}{\"version\"}=\"2.\
0\";\n$PG{\"probcons\"}{\"4_TCOFFEE\"}=\"PROBCONS\\
";\n$PG{\"probcons\"}{\"type\"}=\"sequence_multipl\
e_aligner\";\n$PG{\"probcons\"}{\"ADDRESS\"}=\"htt\
p://probcons.stanford.edu/\";\n$PG{\"probcons\"}{\\
"language2\"}=\"CXX\";\n$PG{\"probcons\"}{\"langua\
ge\"}=\"C++\";\n$PG{\"probcons\"}{\"source\"}=\"ht\
tp://probcons.stanford.edu/probcons_v1_12.tar.gz\"\
;\n$PG{\"probcons\"}{\"mode\"}=\"mcoffee\";\n$PG{\\
"probcons\"}{\"binary\"}=\"probcons\";\n$PG{\"prob\
cons\"}{\"version\"}=\"1.12\";\n$PG{\"mafft\"}{\"4\
_TCOFFEE\"}=\"MAFFT\";\n$PG{\"mafft\"}{\"type\"}=\\
"sequence_multiple_aligner\";\n$PG{\"mafft\"}{\"AD\
DRESS\"}=\"http://align.bmr.kyushu-u.ac.jp/mafft/o\
nline/server/\";\n$PG{\"mafft\"}{\"language\"}=\"C\
\";\n$PG{\"mafft\"}{\"language\"}=\"C\";\n$PG{\"ma\
fft\"}{\"source\"}=\"http://align.bmr.kyushu-u.ac.\
jp/mafft/software/mafft-6.603-with-extensions-src.\
tgz\";\n$PG{\"mafft\"}{\"windows\"}=\"http://align\
.bmr.kyushu-u.ac.jp/mafft/software/mafft-6.603-min\
gw.tar\";\n$PG{\"mafft\"}{\"mode\"}=\"mcoffee,rcof\
fee\";\n$PG{\"mafft\"}{\"binary\"}=\"mafft.tar.gz\\
";\n$PG{\"mafft\"}{\"version\"}=\"6.603\";\n$PG{\"\
muscle\"}{\"4_TCOFFEE\"}=\"MUSCLE\";\n$PG{\"muscle\
\"}{\"type\"}=\"sequence_multiple_aligner\";\n$PG{\
\"muscle\"}{\"ADDRESS\"}=\"http://www.drive5.com/m\
uscle/\";\n$PG{\"muscle\"}{\"language\"}=\"C++\";\\
n$PG{\"muscle\"}{\"language2\"}=\"GPP\";\n$PG{\"mu\
scle\"}{\"source\"}=\"http://www.drive5.com/muscle\
/downloads3.7/muscle3.7_src.tar.gz\";\n$PG{\"muscl\
e\"}{\"windows\"}=\"http://www.drive5.com/muscle/d\
ownloads3.7/muscle3.7_win32.zip\";\n$PG{\"muscle\"\
}{\"linux\"}=\"http://www.drive5.com/muscle/downlo\
ads3.7/muscle3.7_linux_ia32.tar.gz\";\n$PG{\"muscl\
e\"}{\"mode\"}=\"mcoffee,rcoffee\";\n$PG{\"muscle\\
"}{\"version\"}=\"3.7\";\n$PG{\"mus4\"}{\"4_TCOFFE\
E\"}=\"MUS4\";\n$PG{\"mus4\"}{\"type\"}=\"sequence\
_multiple_aligner\";\n$PG{\"mus4\"}{\"ADDRESS\"}=\\
"http://www.drive5.com/muscle/\";\n$PG{\"mus4\"}{\\
"language\"}=\"C++\";\n$PG{\"mus4\"}{\"language2\"\
}=\"GPP\";\n$PG{\"mus4\"}{\"source\"}=\"http://www\
.drive5.com/muscle/muscle4.0_src.tar.gz\";\n$PG{\"\
mus4\"}{\"mode\"}=\"mcoffee,rcoffee\";\n$PG{\"mus4\
\"}{\"version\"}=\"4.0\";\n$PG{\"pcma\"}{\"4_TCOFF\
EE\"}=\"PCMA\";\n$PG{\"pcma\"}{\"type\"}=\"sequenc\
e_multiple_aligner\";\n$PG{\"pcma\"}{\"ADDRESS\"}=\
\"ftp://iole.swmed.edu/pub/PCMA/\";\n$PG{\"pcma\"}\
{\"language\"}=\"C\";\n$PG{\"pcma\"}{\"language2\"\
}=\"C\";\n$PG{\"pcma\"}{\"source\"}=\"ftp://iole.s\
wmed.edu/pub/PCMA/pcma.tar.gz\";\n$PG{\"pcma\"}{\"\
mode\"}=\"mcoffee\";\n$PG{\"pcma\"}{\"version\"}=\\
"1.0\";\n$PG{\"kalign\"}{\"4_TCOFFEE\"}=\"KALIGN\"\
;\n$PG{\"kalign\"}{\"type\"}=\"sequence_multiple_a\
ligner\";\n$PG{\"kalign\"}{\"ADDRESS\"}=\"http://m\
sa.cgb.ki.se\";\n$PG{\"kalign\"}{\"language\"}=\"C\
\";\n$PG{\"kalign\"}{\"language2\"}=\"C\";\n$PG{\"\
kalign\"}{\"source\"}=\"http://msa.cgb.ki.se/downl\
oads/kalign/current.tar.gz\";\n$PG{\"kalign\"}{\"m\
ode\"}=\"mcoffee\";\n$PG{\"kalign\"}{\"version\"}=\
\"1.0\";\n$PG{\"amap\"}{\"4_TCOFFEE\"}=\"AMAP\";\n\
$PG{\"amap\"}{\"type\"}=\"sequence_multiple_aligne\
r\";\n$PG{\"amap\"}{\"ADDRESS\"}=\"http://bio.math\
.berkeley.edu/amap/\";\n$PG{\"amap\"}{\"language\"\
}=\"C++\";\n$PG{\"amap\"}{\"language2\"}=\"CXX\";\\
n$PG{\"amap\"}{\"source\"}=\"http://amap-align.goo\
glecode.com/files/amap.2.0.tar.gz\";\n$PG{\"amap\"\
}{\"mode\"}=\"mcoffee\";\n$PG{\"amap\"}{\"version\\
"}=\"2.0\";\n$PG{\"proda\"}{\"4_TCOFFEE\"}=\"PRODA\
\";\n$PG{\"proda\"}{\"type\"}=\"sequence_multiple_\
aligner\";\n$PG{\"proda\"}{\"ADDRESS\"}=\"http://p\
roda.stanford.edu\";\n$PG{\"proda\"}{\"language\"}\
=\"C++\";\n$PG{\"proda\"}{\"language2\"}=\"CXX\";\\
n$PG{\"proda\"}{\"source\"}=\"http://proda.stanfor\
d.edu/proda_1_0.tar.gz\";\n$PG{\"proda\"}{\"mode\"\
}=\"mcoffee\";\n$PG{\"proda\"}{\"version\"}=\"1.0\\
";\n$PG{\"fsa\"}{\"4_TCOFFEE\"}=\"FSA\";\n$PG{\"fs\
a\"}{\"type\"}=\"sequence_multiple_aligner\";\n$PG\
{\"fsa\"}{\"ADDRESS\"}=\"http://fsa.sourceforge.ne\
t/\";\n$PG{\"fsa\"}{\"language\"}=\"C++\";\n$PG{\"\
fsa\"}{\"language2\"}=\"CXX\";\n$PG{\"fsa\"}{\"sou\
rce\"}=\"http://sourceforge.net/projects/fsa/files\
/fsa-1.15.3.tar.gz/download/\";\n$PG{\"fsa\"}{\"mo\
de\"}=\"mcoffee\";\n$PG{\"fsa\"}{\"version\"}=\"1.\
15.3\";\n$PG{\"prank\"}{\"4_TCOFFEE\"}=\"PRANK\";\\
n$PG{\"prank\"}{\"type\"}=\"sequence_multiple_alig\
ner\";\n$PG{\"prank\"}{\"ADDRESS\"}=\"http://www.e\
bi.ac.uk/goldman-srv/prank/\";\n$PG{\"prank\"}{\"l\
anguage\"}=\"C++\";\n$PG{\"prank\"}{\"language2\"}\
=\"CXX\";\n$PG{\"prank\"}{\"source\"}=\"http://www\
.ebi.ac.uk/goldman-srv/prank/src/prank/prank.src.1\
00303.tgz\";\n$PG{\"prank\"}{\"mode\"}=\"mcoffee\"\
;\n$PG{\"prank\"}{\"version\"}=\"100303\";\n$PG{\"\
sap\"}{\"4_TCOFFEE\"}=\"SAP\";\n$PG{\"sap\"}{\"typ\
e\"}=\"structure_pairwise_aligner\";\n$PG{\"sap\"}\
{\"ADDRESS\"}=\"http://mathbio.nimr.mrc.ac.uk/wiki\
/Software\";\n$PG{\"sap\"}{\"language\"}=\"C\";\n$\
PG{\"sap\"}{\"language2\"}=\"C\";\n$PG{\"sap\"}{\"\
source\"}=\"http://mathbio.nimr.mrc.ac.uk/download\
/sap-1.1.1.tar.gz\";\n$PG{\"sap\"}{\"mode\"}=\"exp\
resso,3dcoffee\";\n$PG{\"sap\"}{\"version\"}=\"1.1\
.1\";\n$PG{\"TMalign\"}{\"4_TCOFFEE\"}=\"TMALIGN\"\
;\n$PG{\"TMalign\"}{\"type\"}=\"structure_pairwise\
_aligner\";\n$PG{\"TMalign\"}{\"ADDRESS\"}=\"http:\
//zhang.bioinformatics.ku.edu/TM-align/TMalign.f\"\
;\n$PG{\"TMalign\"}{\"language\"}=\"Fortran\";\n$P\
G{\"TMalign\"}{\"language2\"}=\"Fortran\";\n$PG{\"\
TMalign\"}{\"source\"}=\"http://zhang.bioinformati\
cs.ku.edu/TM-align/TMalign.f\";\n$PG{\"TMalign\"}{\
\"linux\"}=\"http://zhang.bioinformatics.ku.edu/TM\
-align/TMalign_32.gz\";\n$PG{\"TMalign\"}{\"mode\"\
}=\"expresso,3dcoffee\";\n$PG{\"TMalign\"}{\"versi\
on\"}=\"1.0\";\n$PG{\"mustang\"}{\"4_TCOFFEE\"}=\"\
MUSTANG\";\n$PG{\"mustang\"}{\"type\"}=\"structure\
_pairwise_aligner\";\n$PG{\"mustang\"}{\"ADDRESS\"\
}=\"http://www.cs.mu.oz.au/~arun/mustang\";\n$PG{\\
"mustang\"}{\"language\"}=\"C++\";\n$PG{\"mustang\\
"}{\"language2\"}=\"CXX\";\n$PG{\"mustang\"}{\"sou\
rce\"}=\"http://ww2.cs.mu.oz.au/~arun/mustang/must\
ang_v3.2.1.tgz\";\n$PG{\"mustang\"}{\"mode\"}=\"ex\
presso,3dcoffee\";\n$PG{\"mustang\"}{\"version\"}=\
\"3.2.1\";\n$PG{\"lsqman\"}{\"4_TCOFFEE\"}=\"LSQMA\
N\";\n$PG{\"lsqman\"}{\"type\"}=\"structure_pairwi\
se_aligner\";\n$PG{\"lsqman\"}{\"ADDRESS\"}=\"empt\
y\";\n$PG{\"lsqman\"}{\"language\"}=\"empty\";\n$P\
G{\"lsqman\"}{\"language2\"}=\"empty\";\n$PG{\"lsq\
man\"}{\"source\"}=\"empty\";\n$PG{\"lsqman\"}{\"u\
pdate_action\"}=\"never\";\n$PG{\"lsqman\"}{\"mode\
\"}=\"expresso,3dcoffee\";\n$PG{\"align_pdb\"}{\"4\
_TCOFFEE\"}=\"ALIGN_PDB\";\n$PG{\"align_pdb\"}{\"t\
ype\"}=\"structure_pairwise_aligner\";\n$PG{\"alig\
n_pdb\"}{\"ADDRESS\"}=\"empty\";\n$PG{\"align_pdb\\
"}{\"language\"}=\"empty\";\n$PG{\"align_pdb\"}{\"\
language2\"}=\"empty\";\n$PG{\"align_pdb\"}{\"sour\
ce\"}=\"empty\";\n$PG{\"align_pdb\"}{\"update_acti\
on\"}=\"never\";\n$PG{\"align_pdb\"}{\"mode\"}=\"e\
xpresso,3dcoffee\";\n$PG{\"fugueali\"}{\"4_TCOFFEE\
\"}=\"FUGUE\";\n$PG{\"fugueali\"}{\"type\"}=\"stru\
cture_pairwise_aligner\";\n$PG{\"fugueali\"}{\"ADD\
RESS\"}=\"http://www-cryst.bioc.cam.ac.uk/fugue/do\
wnload.html\";\n$PG{\"fugueali\"}{\"language\"}=\"\
empty\";\n$PG{\"fugueali\"}{\"language2\"}=\"empty\
\";\n$PG{\"fugueali\"}{\"source\"}=\"empty\";\n$PG\
{\"fugueali\"}{\"update_action\"}=\"never\";\n$PG{\
\"fugueali\"}{\"mode\"}=\"expresso,3dcoffee\";\n$P\
G{\"dalilite.pl\"}{\"4_TCOFFEE\"}=\"DALILITEc\";\n\
$PG{\"dalilite.pl\"}{\"type\"}=\"structure_pairwis\
e_aligner\";\n$PG{\"dalilite.pl\"}{\"ADDRESS\"}=\"\
built_in\";\n$PG{\"dalilite.pl\"}{\"ADDRESS2\"}=\"\
http://www.ebi.ac.uk/Tools/webservices/services/da\
lilite\";\n$PG{\"dalilite.pl\"}{\"language\"}=\"Pe\
rl\";\n$PG{\"dalilite.pl\"}{\"language2\"}=\"Perl\\
";\n$PG{\"dalilite.pl\"}{\"source\"}=\"empty\";\n$\
PG{\"dalilite.pl\"}{\"update_action\"}=\"never\";\\
n$PG{\"dalilite.pl\"}{\"mode\"}=\"expresso,3dcoffe\
e\";\n$PG{\"probconsRNA\"}{\"4_TCOFFEE\"}=\"PROBCO\
NSRNA\";\n$PG{\"probconsRNA\"}{\"type\"}=\"RNA_mul\
tiple_aligner\";\n$PG{\"probconsRNA\"}{\"ADDRESS\"\
}=\"http://probcons.stanford.edu/\";\n$PG{\"probco\
nsRNA\"}{\"language\"}=\"C++\";\n$PG{\"probconsRNA\
\"}{\"language2\"}=\"CXX\";\n$PG{\"probconsRNA\"}{\
\"source\"}=\"http://probcons.stanford.edu/probcon\
sRNA.tar.gz\";\n$PG{\"probconsRNA\"}{\"mode\"}=\"m\
coffee,rcoffee\";\n$PG{\"probconsRNA\"}{\"version\\
"}=\"1.0\";\n$PG{\"sfold\"}{\"4_TCOFFEE\"}=\"CONSA\
N\";\n$PG{\"sfold\"}{\"type\"}=\"RNA_pairwise_alig\
ner\";\n$PG{\"sfold\"}{\"ADDRESS\"}=\"http://selab\
.janelia.org/software/consan/\";\n$PG{\"sfold\"}{\\
"language\"}=\"empty\";\n$PG{\"sfold\"}{\"language\
2\"}=\"empty\";\n$PG{\"sfold\"}{\"source\"}=\"empt\
y\";\n$PG{\"sfold\"}{\"update_action\"}=\"never\";\
\n$PG{\"sfold\"}{\"mode\"}=\"rcoffee\";\n$PG{\"RNA\
plfold\"}{\"4_TCOFFEE\"}=\"RNAPLFOLD\";\n$PG{\"RNA\
plfold\"}{\"type\"}=\"RNA_secondarystructure_predi\
ctor\";\n$PG{\"RNAplfold\"}{\"ADDRESS\"}=\"http://\
www.tbi.univie.ac.at/~ivo/RNA/\";\n$PG{\"RNAplfold\
\"}{\"language\"}=\"C\";\n$PG{\"RNAplfold\"}{\"lan\
guage2\"}=\"C\";\n$PG{\"RNAplfold\"}{\"source\"}=\\
"http://www.tbi.univie.ac.at/~ivo/RNA/ViennaRNA-1.\
7.2.tar.gz\";\n$PG{\"RNAplfold\"}{\"mode\"}=\"rcof\
fee,\";\n$PG{\"RNAplfold\"}{\"version\"}=\"1.7.2\"\
;\n$PG{\"retree\"}{\"4_TCOFFEE\"}=\"PHYLIP\";\n$PG\
{\"retree\"}{\"type\"}=\"RNA_secondarystructure_pr\
edictor\";\n$PG{\"retree\"}{\"ADDRESS\"}=\"http://\
evolution.gs.washington.edu/phylip/\";\n$PG{\"retr\
ee\"}{\"language\"}=\"C\";\n$PG{\"retree\"}{\"lang\
uage2\"}=\"C\";\n$PG{\"retree\"}{\"source\"}=\"htt\
p://evolution.gs.washington.edu/phylip/download/ph\
ylip-3.69.tar.gz\";\n$PG{\"retree\"}{\"mode\"}=\"t\
rmsd,\";\n$PG{\"retree\"}{\"version\"}=\"3.69\";\n\
$PG{\"hmmtop\"}{\"4_TCOFFEE\"}=\"HMMTOP\";\n$PG{\"\
hmmtop\"}{\"type\"}=\"protein_secondarystructure_p\
redictor\";\n$PG{\"hmmtop\"}{\"ADDRESS\"}=\"www.en\
zim.hu/hmmtop/\";\n$PG{\"hmmtop\"}{\"language\"}=\\
"C\";\n$PG{\"hmmtop\"}{\"language2\"}=\"C\";\n$PG{\
\"hmmtop\"}{\"source\"}=\"empty\";\n$PG{\"hmmtop\"\
}{\"update_action\"}=\"never\";\n$PG{\"hmmtop\"}{\\
"mode\"}=\"tcoffee\";\n$PG{\"gorIV\"}{\"4_TCOFFEE\\
"}=\"GOR4\";\n$PG{\"gorIV\"}{\"type\"}=\"protein_s\
econdarystructure_predictor\";\n$PG{\"gorIV\"}{\"A\
DDRESS\"}=\"http://mig.jouy.inra.fr/logiciels/gorI\
V/\";\n$PG{\"gorIV\"}{\"language\"}=\"C\";\n$PG{\"\
gorIV\"}{\"language2\"}=\"C\";\n$PG{\"gorIV\"}{\"s\
ource\"}=\"http://mig.jouy.inra.fr/logiciels/gorIV\
/GOR_IV.tar.gz\";\n$PG{\"gorIV\"}{\"update_action\\
"}=\"never\";\n$PG{\"gorIV\"}{\"mode\"}=\"tcoffee\\
";\n$PG{\"wublast.pl\"}{\"4_TCOFFEE\"}=\"EBIWUBLAS\
Tc\";\n$PG{\"wublast.pl\"}{\"type\"}=\"protein_hom\
ology_predictor\";\n$PG{\"wublast.pl\"}{\"ADDRESS\\
"}=\"built_in\";\n$PG{\"wublast.pl\"}{\"ADDRESS2\"\
}=\"http://www.ebi.ac.uk/Tools/webservices/service\
s/wublast\";\n$PG{\"wublast.pl\"}{\"language\"}=\"\
Perl\";\n$PG{\"wublast.pl\"}{\"language2\"}=\"Perl\
\";\n$PG{\"wublast.pl\"}{\"source\"}=\"empty\";\n$\
PG{\"wublast.pl\"}{\"update_action\"}=\"never\";\n\
$PG{\"wublast.pl\"}{\"mode\"}=\"psicoffee,expresso\
,accurate\";\n$PG{\"blastpgp.pl\"}{\"4_TCOFFEE\"}=\
\"EBIBLASTPGPc\";\n$PG{\"blastpgp.pl\"}{\"type\"}=\
\"protein_homology_predictor\";\n$PG{\"blastpgp.pl\
\"}{\"ADDRESS\"}=\"built_in\";\n$PG{\"blastpgp.pl\\
"}{\"ADDRESS2\"}=\"http://www.ebi.ac.uk/Tools/webs\
ervices/services/blastpgp\";\n$PG{\"blastpgp.pl\"}\
{\"language\"}=\"Perl\";\n$PG{\"blastpgp.pl\"}{\"l\
anguage2\"}=\"Perl\";\n$PG{\"blastpgp.pl\"}{\"sour\
ce\"}=\"empty\";\n$PG{\"blastpgp.pl\"}{\"update_ac\
tion\"}=\"never\";\n$PG{\"blastpgp.pl\"}{\"mode\"}\
=\"psicoffee,expresso,accurate\";\n$PG{\"blastcl3\\
"}{\"4_TCOFFEE\"}=\"NCBIWEBBLAST\";\n$PG{\"blastcl\
3\"}{\"type\"}=\"protein_homology_predictor\";\n$P\
G{\"blastcl3\"}{\"ADDRESS\"}=\"ftp://ftp.ncbi.nih.\
gov/blast/executables/LATEST\";\n$PG{\"blastcl3\"}\
{\"language\"}=\"C\";\n$PG{\"blastcl3\"}{\"languag\
e2\"}=\"C\";\n$PG{\"blastcl3\"}{\"source\"}=\"empt\
y\";\n$PG{\"blastcl3\"}{\"update_action\"}=\"never\
\";\n$PG{\"blastcl3\"}{\"mode\"}=\"psicoffee,expre\
sso,3dcoffee\";\n$PG{\"blastpgp\"}{\"4_TCOFFEE\"}=\
\"NCBIBLAST\";\n$PG{\"blastpgp\"}{\"type\"}=\"prot\
ein_homology_predictor\";\n$PG{\"blastpgp\"}{\"ADD\
RESS\"}=\"ftp://ftp.ncbi.nih.gov/blast/executables\
/LATEST\";\n$PG{\"blastpgp\"}{\"language\"}=\"C\";\
\n$PG{\"blastpgp\"}{\"language2\"}=\"C\";\n$PG{\"b\
lastpgp\"}{\"source\"}=\"empty\";\n$PG{\"blastpgp\\
"}{\"update_action\"}=\"never\";\n$PG{\"blastpgp\"\
}{\"mode\"}=\"psicoffee,expresso,3dcoffee\";\n$PG{\
\"SOAP::Lite\"}{\"4_TCOFFEE\"}=\"SOAPLITE\";\n$PG{\
\"SOAP::Lite\"}{\"type\"}=\"library\";\n$PG{\"SOAP\
::Lite\"}{\"ADDRESS\"}=\"http://cpansearch.perl.or\
g/src/MKUTTER/SOAP-Lite-0.710.08/Makefile.PL\";\n$\
PG{\"SOAP::Lite\"}{\"language\"}=\"Perl\";\n$PG{\"\
SOAP::Lite\"}{\"language2\"}=\"Perl\";\n$PG{\"SOAP\
::Lite\"}{\"source\"}=\"empty\";\n$PG{\"blastpgp\"\
}{\"update_action\"}=\"never\";\n$PG{\"SOAP::Lite\\
"}{\"mode\"}=\"none\";\n$PG{\"XML::Simple\"}{\"4_T\
COFFEE\"}=\"XMLSIMPLE\";\n$PG{\"XML::Simple\"}{\"t\
ype\"}=\"library\";\n$PG{\"XML::Simple\"}{\"ADDRES\
S\"}=\"http://search.cpan.org/~grantm/XML-Simple-2\
.18/lib/XML/Simple.pm\";\n$PG{\"XML::Simple\"}{\"l\
anguage\"}=\"Perl\";\n$PG{\"XML::Simple\"}{\"langu\
age2\"}=\"Perl\";\n$PG{\"XML::Simple\"}{\"source\"\
}=\"empty\";\n$PG{\"XML::Simple\"}{\"mode\"}=\"psi\
coffee,expresso,accurate\";\n$MODE{\"tcoffee\"}{\"\
name\"}=\"tcoffee\";\n$MODE{\"rcoffee\"}{\"name\"}\
=\"rcoffee\";\n$MODE{\"3dcoffee\"}{\"name\"}=\"3dc\
offee\";\n$MODE{\"mcoffee\"}{\"name\"}=\"mcoffee\"\
;\n$MODE{\"expresso\"}{\"name\"}=\"expresso\";\n$M\
ODE{\"trmsd\"}{\"name\"}=\"trmsd\";\n$MODE{\"accur\
ate\"}{\"name\"}=\"accurate\";\n$MODE{\"seq_reform\
at\"}{\"name\"}=\"seq_reformat\";\n\n\n$PG{C}{comp\
iler}=\"gcc\";\n$PG{C}{compiler_flag}=\"CC\";\n$PG\
{C}{options}=\"\";\n$PG{C}{options_flag}=\"CFLAGS\\
";\n$PG{C}{type}=\"compiler\";\n\n$PG{\"CXX\"}{com\
piler}=\"g++\";\n$PG{\"CXX\"}{compiler_flag}=\"CXX\
\";\n$PG{\"CXX\"}{options}=\"\";\n$PG{\"CXX\"}{opt\
ions_flag}=\"CXXFLAGS\";\n$PG{CXX}{type}=\"compile\
r\";\n\n$PG{\"CPP\"}{compiler}=\"g++\";\n$PG{\"CPP\
\"}{compiler_flag}=\"CPP\";\n$PG{\"CPP\"}{options}\
=\"\";\n$PG{\"CPP\"}{options_flag}=\"CPPFLAGS\";\n\
$PG{CPP}{type}=\"compiler\";\n\n$PG{\"GPP\"}{compi\
ler}=\"g++\";\n$PG{\"GPP\"}{compiler_flag}=\"GPP\"\
;\n$PG{\"GPP\"}{options}=\"\";\n$PG{\"GPP\"}{optio\
ns_flag}=\"CFLAGS\";\n$PG{GPP}{type}=\"compiler\";\
\n\n$PG{Fortran}{compiler}=\"g77\";\n$PG{Fortran}{\
compiler_flag}=\"FCC\";\n$PG{Fortran}{type}=\"comp\
iler\";\n\n$PG{Perl}{compiler}=\"CPAN\";\n$PG{Perl\
}{type}=\"compiler\";\n\n$SUPPORTED_OS{macox}=\"Ma\
cintosh\";\n$SUPPORTED_OS{linux}=\"Linux\";\n$SUPP\
ORTED_OS{windows}=\"Cygwin\";\n\n\n\n$MODE{t_coffe\
e}{description}=\" for regular multiple sequence a\
lignments\";\n$MODE{rcoffee} {description}=\" for \
RNA multiple sequence alignments\";\n\n$MODE{psico\
ffee} {description}=\" for Homology Extended multi\
ple sequence alignments\";\n$MODE{expresso}{descri\
ption}=\" for very accurate structure based multip\
le sequence alignments\";\n$MODE{\"3dcoffee\"}{des\
cription}=\" for multiple structure alignments\";\\
n$MODE{mcoffee} {description}=\" for combining alt\
ernative multiple sequence alignment packages\\n--\
----- into a unique meta-package. The installer wi\
ll upload several MSA packages and compile them\\n\
\n\";\n\n\n&post_process_PG();\nreturn;\n}\n\nsub \
post_process_PG\n  {\n    my $p;\n    \n    %PG=&n\
ame2dname (%PG);\n    %MODE=&name2dname(%MODE);\n \
   foreach $p (keys(%PG)){if ( $PG{$p}{type} eq \"\
compiler\"){$PG{$p}{update_action}=\"never\";}}\n \
   \n  }\n\nsub name2dname\n  {\n    my (%L)=(@_);\
\n    my ($l, $ml);\n    \n    foreach my $pg (key\
s(%L))\n      {\n	$l=length ($pg);\n	if ( $l>$ml){\
$ml=$l;}\n      }\n    $ml+=1;\n    foreach my $pg\
 (keys(%L))\n      {\n	my $name;\n	$l=$ml-length (\
$pg);\n	$name=$pg;\n	for ( $b=0; $b<$l; $b++)\n	  \
{\n	    $name .=\" \";\n	  }\n	$L{$pg}{dname}=$nam\
e;\n      }\n    return %L;\n  }\n\nsub env_file2p\
utenv\n  {\n    my $f=@_[0];\n    my $F=new FileHa\
ndle;\n    my $n;\n    \n    open ($F, \"$f\");\n \
   while (<$F>)\n      {\n	my $line=$_;\n	my($var,\
 $value)=($_=~/(\\S+)\\=(\\S*)/);\n	$ENV{$var}=$va\
lue;\n	$ENV_SET{$var}=1;\n	$n++;\n      }\n    clo\
se ($F);\n    return $n;\n  }\n\nsub replace_line_\
in_file\n  {\n    my ($file, $wordin, $wordout)=@_\
;\n    my $O=new FileHandle;\n    my $I=new FileHa\
ndle;\n    my $l;\n    if (!-e $file){return;}\n  \
  \n    system (\"mv $file $file.old\");\n    open\
 ($O, \">$file\");\n    open ($I, \"$file.old\");\\
n    while (<$I>)\n      {\n	$l=$_;\n	if (!($l=~/$\
wordin/)){print $O \"$l\";}\n	elsif ( $wordout ne \
\"\"){$l=~s/$wordin/$wordout/g;print $O \"$l\";}\n\
      }\n    close ($O);\n    close ($I);\n    ret\
urn;\n  }\n\nsub add_C_libraries\n  {\n   my ($fil\
e,$first,@list)=@_;\n   \n    my $O=new FileHandle\
;\n    my $I=new FileHandle;\n    my ($l,$anchor);\
\n    if (!-e $file){return;}\n   \n    $anchor=\"\
#include <$first>\";\n	 \n    system (\"mv $file $\
file.old\");\n    open ($O, \">$file\");\n    open\
 ($I, \"$file.old\");\n    while (<$I>)\n      {\n\
	$l=$_;\n	print $O \"$l\";\n	if (!($l=~/$anchor/))\
\n	   {\n	    \n	    foreach my $lib (@list)\n	   \
    {\n                  print $O \"#include <$lib\
>\\n\";\n	       }\n           }\n      }\n    clo\
se ($O);\n    close ($I);\n    return;\n    }\n","\
use Env;\nuse Cwd;\n@suffix=(\"tmp\", \"temp\", \"\
cache\", \"t_coffee\", \"core\", \"tcoffee\");\n\n\
if ($#ARGV==-1)\n  {\n    print \"clean_cache.pl -\
file <file to add in -dir> -dir=<dir> -size=<value\
 in Mb>\\n0: unlimited -1 always.\\nWill only clea\
n directories matching:[\";\n    foreach $k(@suffi\
x){print \"*$k* \";}\n    print \"]\\n\";\n    exi\
t (EXIT_FAILURE);\n  }\n\n$cl=join (\" \",@ARGV);\\
nif (($cl=~/\\-no_action/))\n  {\n    exit (EXIT_S\
UCCESS);\n  }\n\nif (($cl=~/\\-debug/))\n  {\n    \
$DEBUG=1;\n  }\nelse\n  {\n    $DEBUG=0;\n  }\n\ni\
f (($cl=~/\\-dir=(\\S+)/))\n  {\n    $dir=$1;\n  }\
\nelse\n  {\n    $dir=\"./\";\n  }\n\nif ($cl=~/\\\
-file=(\\S+)/)\n  {\n    $file=$1;\n  }\nelse\n  {\
\n    $file=0;\n  }\n\nif ($cl=~/\\-size=(\\S+)/)\\
n  {\n    $max_size=$1;\n  }\nelse\n  {\n    $max_\
size=0;#unlimited\n  }\nif ($cl=~/\\-force/)\n  {\\
n    $force=1;\n  }\nelse\n  {\n    $force=0;\n  }\
\n\nif ($cl=~/\\-age=(\\S+)/)\n  {\n    $max_age=$\
1;\n  }\nelse\n  {\n    $max_age=0;#unlimited\n  }\
\n\n$max_size*=1000000;\nif ( ! -d $dir)\n  {\n   \
 print STDERR \"\\nCannot process $dir: does not e\
xist \\n\";\n    exit (EXIT_FAILURE);\n  }\n\nif (\
 !($dir=~/^\\//))\n  {\n    $base=cwd();\n    $dir\
=\"$base/$dir\";\n  }\n\n$proceed=0;\nforeach $s (\
@suffix)\n  {\n    \n    if (($dir=~/$s/)){$procee\
d=1;}\n    $s=uc ($s);\n    if (($dir=~/$s/)){$pro\
ceed=1;}\n  }\nif ( $proceed==0)\n  {\n    print S\
TDERR \"Clean_cache.pl can only clean directories \
whose absolute path name contains the following st\
rings:\";\n    foreach $w (@suffix) {print STDERR \
\"$w \";$w=lc($w); print STDERR \"$w \";}\n    pri\
nt STDERR \"\\nCannot process $dir\\n\";\n    exit\
 (EXIT_FAILURE);\n  }\n\n$name_file=\"$dir/name_fi\
le.txt\";\n$size_file=\"$dir/size_file.txt\";\nif \
( $force){&create_ref_file ($dir,$name_file,$size_\
file);}\nif ($file){&add_file ($dir, $name_file, $\
size_file, $file);}\n&clean_dir ($dir, $name_file,\
 $size_file, $max_size,$max_age);\nexit (EXIT_SUCC\
ESS);\n\nsub clean_dir \n  {\n    my ($dir, $name_\
file, $size_file, $max_size, $max_age)=@_;\n    my\
 ($tot_size, $size, $f, $s);\n\n  \n    $tot_size=\
&get_tot_size ($dir, $name_file, $size_file);\n\n \
   if ( $tot_size<=$max_size){return ;}\n    else \
{$max_size/=2;}\n    \n    #recreate the name file\
 in case some temprary files have not been properl\
y registered\n    &create_ref_file ($dir, $name_fi\
le, $size_file, $max_age);\n  \n    $new_name_file\
=&vtmpnam();\n    open (R, \"$name_file\");\n    o\
pen (W, \">$new_name_file\");\n    while (<R>)\n  \
    {\n	my $line=$_;\n	\n	($f, $s)=($line=~/(\\S+)\
 (\\S+)/);\n	if ( !($f=~/\\S/)){next;}\n	\n	elsif \
($max_size && $tot_size>=$max_size && !($f=~/name_\
file/))\n	  {\n	    remove ( \"$dir/$f\");\n	    $\
tot_size-=$s;\n	  }\n	elsif ( $max_age && -M(\"$di\
r/$f\")>=$max_age)\n	  {\n	    remove ( \"$dir/$f\\
");\n	    $tot_size-=$s;\n	  }\n	else\n	  {\n	    \
print W \"$f $s\\n\";\n	  }\n      }\n    close (R\
);\n    close (W);\n    open (F, \">$size_file\");\
\n    print F \"$tot_size\";\n    if ( -e $new_nam\
e_file){`mv $new_name_file $name_file`;}\n    clos\
e (F);\n  }\nsub get_tot_size\n  {\n    my ($dir, \
$name_file, $size_file)=@_;\n    my $size;\n    \n\
    if ( !-d $dir){return 0;}\n    if ( !-e $name_\
file)\n      {\n	\n	&create_ref_file ($dir, $name_\
file, $size_file);\n      }\n    open (F, \"$size_\
file\");\n    $size=<F>;\n    close (F);\n    chom\
p ($size);\n    return $size;\n  }\nsub size \n  {\
\n    my $f=@_[0];\n\n    if ( !-d $f){return -s($\
f);}\n    else {return &dir2size($f);}\n  }\nsub d\
ir2size\n  {\n    my $d=@_[0];\n    my ($s, $f);\n\
    \n    if ( !-d $d) {return 0;}\n    \n    fore\
ach $f (&dir2list ($d))\n      {\n	if ( -d $f){$s+\
=&dir2size (\"$d/$f\");}\n	else {$s+= -s \"$dir/$f\
\";}\n      }\n    return $s;\n  }\n\nsub remove \\
n  {\n    my $file=@_[0];\n    my ($f);\n    \n   \
 debug_print( \"--- $file ---\\n\");\n    if (($fi\
le eq \".\") || ($file eq \"..\") || ($file=~/\\*/\
)){return EXIT_FAILURE;}\n    elsif ( !-d $file)\n\
      {\n	debug_print (\"unlink $file\\n\");\n	if \
(-e $file){unlink ($file);}\n      }\n    elsif ( \
-d $file)\n      {\n	debug_print (\"++++++++ $file\
 +++++++\\n\");\n	foreach $f (&dir2list($file))\n	\
  {\n	    &remove (\"$file/$f\");\n	  }\n	debug_pr\
int (\"rmdir $file\\n\");\n	rmdir $file;\n      }\\
n    else\n      {\n	debug_print (\"????????? $fil\
e ????????\\n\");\n      }\n    return EXIT_SUCCES\
S;\n  }\n\nsub dir2list\n  {\n    my $dir=@_[0];\n\
    my (@list1, @list2,@list3, $l);\n\n    opendir\
 (DIR,$dir);\n    @list1=readdir (DIR);\n    close\
dir (DIR);\n    \n    foreach $l (@list1)\n      {\
\n	if ( $l ne \".\" && $l ne \"..\"){@list2=(@list\
2, $l);}\n      }\n    @list3 = sort { (-M \"$dir/\
$list2[$b]\") <=> (-M \"$dir/$list2[$a]\")} @list2\
;\n    return @list3;\n    \n  }\n\nsub debug_prin\
t\n  {\n    \n    if ($DEBUG==1){print @_;}\n    \\
n  }\nsub create_ref_file\n  {\n    my ($dir,$name\
_file,$size_file)=@_;\n    my ($f, $s, $tot_size, \
@l);\n    \n    if ( !-d $dir){return;}\n    \n   \
 @l=&dir2list ($dir);\n    open (F, \">$name_file\\
");\n    foreach $f (@l)\n      {\n	$s=&size(\"$di\
r/$f\");\n	$tot_size+=$s;\n	print F \"$f $s\\n\";\\
n      }\n    &myecho ($tot_size, \">$size_file\")\
;\n    close (F);\n  }\nsub add_file \n  {\n    my\
 ($dir,$name_file,$size_file,$file)=@_;\n    my ($\
s, $tot_size);\n    \n    if ( !-d $dir)   {return\
;}\n    if ( !-e \"$dir/$file\" ) {return;}\n    i\
f ( !-e $name_file){&create_ref_file ($dir,$name_f\
ile,$size_file);}\n					    \n    $s=&size(\"$dir/\
$file\");\n    open (F, \">>$name_file\");\n    pr\
int F \"$file\\n\";\n    close (F);\n\n    $tot_si\
ze=&get_tot_size ($dir,$name_file,$size_file);\n  \
  $tot_size+=$s;\n    &myecho ($tot_size, \">$size\
_file\");\n    \n  }\n	\nsub myecho\n  {\n    my (\
$string, $file)=@_;\n    open (ECHO, $file) || die\
;\n    print ECHO \"$string\";\n    close (ECHO);\\
n  }\n    \n		\n	\nsub vtmpnam\n  {\n    my $tmp_f\
ile_name;\n    $tmp_name_counter++;\n    $tmp_file\
_name=\"tmp_file_for_clean_cache_pdb$$.$tmp_name_c\
ounter\";\n    $tmp_file_list[$ntmp_file++]=$tmp_f\
ile_name;\n    if ( -e $tmp_file_name) {return &vt\
mpnam ();}\n    else {return $tmp_file_name;}\n  }\
\n","\n$t_coffee=\"t_coffee\";\n\nforeach $value (\
 @ARGV)\n  {\n    $seq_file=$seq_file.\" \".$value\
;\n  }\n\n$name=$ARGV[0];\n$name=~s/\\.[^\\.]*$//;\
\n$lib_name=\"$name.mocca_lib\";\n$type=`t_coffee \
$seq_file -get_type -quiet`;\nchop ($type);\n\nif \
( $type eq \"PROTEIN\"){$lib_mode=\"lalign_rs_s_pa\
ir -lalign_n_top 20\";}\nelsif ( $type eq\"DNA\"){\
$lib_mode=\"lalign_rs_s_dna_pair -lalign_n_top 40\\
";}\n\nif ( !(-e $lib_name))\n  {\n	  \n  $command\
=\"$t_coffee -mocca -seq_weight=no -cosmetic_penal\
ty=0 -mocca_interactive -in $lib_mode -out_lib $li\
b_name -infile $seq_file\";\n  \n  }\nelsif ( (-e \
$lib_name))\n  {\n  $command=\"$t_coffee -mocca -s\
eq_weight=no -cosmetic_penalty=0 -mocca_interactiv\
e -in $lib_name -infile $seq_file\";\n  \n  }\n\ns\
ystem ($command);\n\nexit;\n\n","my $WSDL = 'http:\
//www.ebi.ac.uk/Tools/webservices/wsdl/WSDaliLite.\
wsdl';\n\nuse SOAP::Lite;\nuse Data::Dumper;\nuse \
Getopt::Long qw(:config no_ignore_case bundling);\\
nuse File::Basename;\n\nmy $checkInterval = 5;\n\n\
my %params=(\n	    'async' => '1', # Use async mod\
e and simulate sync mode in client\n	    );\nGetOp\
tions(\n    'pdb1=s'     => \\$params{'sequence1'}\
,\n    'chainid1=s' => \\$params{'chainid1'},\n   \
 'pdb2=s'     => \\$params{'sequence2'},\n    'cha\
inid2=s' => \\$params{'chainid2'},\n    \"help|h\"\
	 => \\$help, # Usage info\n    \"async|a\"	 => \\\
$async, # Asynchronous submission\n    \"polljob\"\
	 => \\$polljob, # Get results\n    \"status\"	 =>\
 \\$status, # Get status\n    \"jobid|j=s\"  => \\\
$jobid, # JobId\n    \"email|S=s\"  => \\$params{e\
mail}, # E-mail address\n    \"trace\"      => \\$\
trace, # SOAP messages\n    \"sequence=s\" => \\$s\
equence, # Input PDB\n    );\n\nmy $scriptName = b\
asename($0, ());\nif($help) {\n    &usage();\n    \
exit(0);\n}\n\nif($trace) {\n    print \"Tracing a\
ctive\\n\";\n    SOAP::Lite->import(+trace => 'deb\
ug');\n}\n\nmy $soap = SOAP::Lite\n    ->service($\
WSDL)\n    ->on_fault(sub {\n        my $soap = sh\
ift;\n        my $res = shift;\n        # Throw an\
 exception for all faults\n        if(ref($res) eq\
 '') {\n            die($res);\n        } else {\n\
            die($res->faultstring);\n        }\n  \
      return new SOAP::SOM;\n    }\n              \
 );\n\nif( !($polljob || $status) &&\n    !( defin\
ed($params{'sequence1'}) && defined($params{'seque\
nce2'}) )\n    ) {\n    print STDERR 'Error: bad o\
ption combination', \"\\n\";\n    &usage();\n    e\
xit(1);\n}\nelsif($polljob && defined($jobid)) {\n\
    print \"Getting results for job $jobid\\n\";\n\
    getResults($jobid);\n}\nelsif($status && defin\
ed($jobid)) {\n    print STDERR \"Getting status f\
or job $jobid\\n\";\n    my $result = $soap->check\
Status($jobid);\n    print STDOUT \"$result\", \"\\
\n\";\n    if($result eq 'DONE') {\n	print STDERR \
\"To get results: $scriptName --polljob --jobid $j\
obid\\n\";\n    }\n}\nelse {\n    if(-f $params{'s\
equence1'}) {\n	$params{'sequence1'} = read_file($\
params{'sequence1'});\n    }\n    if(-f $params{'s\
equence2'}) {\n	$params{'sequence2'} = read_file($\
params{'sequence2'});\n    }\n\n    my $jobid;\n  \
  my $paramsData = SOAP::Data->name('params')->typ\
e(map=>\\%params);\n    # For SOAP::Lite 0.60 and \
earlier parameters are passed directly\n    if($SO\
AP::Lite::VERSION eq '0.60' || $SOAP::Lite::VERSIO\
N =~ /0\\.[1-5]/) {\n        $jobid = $soap->runDa\
liLite($paramsData);\n    }\n    # For SOAP::Lite \
0.69 and later parameter handling is different, so\
 pass\n    # undef's for templated params, and the\
n pass the formatted args.\n    else {\n        $j\
obid = $soap->runDaliLite(undef,\n				     $params\
Data);\n    }\n\n    if (defined($async)) {\n	prin\
t STDOUT $jobid, \"\\n\";\n        print STDERR \"\
To check status: $scriptName --status --jobid $job\
id\\n\";\n    } else { # Synchronous mode\n       \
 print STDERR \"JobId: $jobid\\n\";\n        sleep\
 1;\n        getResults($jobid);\n    }\n}\n\nsub \
clientPoll($) {\n    my $jobid = shift;\n    my $r\
esult = 'PENDING';\n    # Check status and wait if\
 not finished\n    #print STDERR \"Checking status\
: $jobid\\n\";\n    while($result eq 'RUNNING' || \
$result eq 'PENDING') {\n        $result = $soap->\
checkStatus($jobid);\n        print STDERR \"$resu\
lt\\n\";\n        if($result eq 'RUNNING' || $resu\
lt eq 'PENDING') {\n            # Wait before poll\
ing again.\n            sleep $checkInterval;\n   \
     }\n    }\n}\n\nsub getResults($) {\n    $jobi\
d = shift;\n    # Check status, and wait if not fi\
nished\n    clientPoll($jobid);\n    # Use JobId i\
f output file name is not defined\n    unless(defi\
ned($outfile)) {\n        $outfile=$jobid;\n    }\\
n    # Get list of data types\n    my $resultTypes\
 = $soap->getResults($jobid);\n    # Get the data \
and write it to a file\n    if(defined($outformat)\
) { # Specified data type\n        my $selResultTy\
pe;\n        foreach my $resultType (@$resultTypes\
) {\n            if($resultType->{type} eq $outfor\
mat) {\n                $selResultType = $resultTy\
pe;\n            }\n        }\n        $res=$soap-\
>poll($jobid, $selResultType->{type});\n        wr\
ite_file($outfile.'.'.$selResultType->{ext}, $res)\
;\n    } else { # Data types available\n        # \
Write a file for each output type\n        for my \
$resultType (@$resultTypes){\n            #print \\
"Getting $resultType->{type}\\n\";\n            $r\
es=$soap->poll($jobid, $resultType->{type});\n    \
        write_file($outfile.'.'.$resultType->{ext}\
, $res);\n        }\n    }\n}\n\nsub read_file($) \
{\n    my $filename = shift;\n    open(FILE, $file\
name);\n    my $content;\n    my $buffer;\n    whi\
le(sysread(FILE, $buffer, 1024)) {\n	$content.= $b\
uffer;\n    }\n    close(FILE);\n    return $conte\
nt;\n}\n\nsub write_file($$) {\n    my ($tmp,$enti\
ty) = @_;\n    print STDERR \"Creating result file\
: \".$tmp.\"\\n\";\n    unless(open (FILE, \">$tmp\
\")) {\n	return 0;\n    }\n    syswrite(FILE, $ent\
ity);\n    close (FILE);\n    return 1;\n}\n\nsub \
usage {\n    print STDERR <<EOF\nDaliLite\n=======\
=\n\nPairwise comparison of protein structures\n\n\
[Required]\n\n  --pdb1                : str  : PDB\
 ID for structure 1\n  --pdb2                : str\
  : PDB ID for structure 2\n\n[Optional]\n\n  --ch\
ain1              : str  : Chain identifer in stru\
cture 1\n  --chain2              : str  : Chain id\
entifer in structure 2\n\n[General]\n\n  -h, --hel\
p            :      : prints this help text\n  -S,\
 --email           : str  : user email address\n  \
-a, --async           :      : asynchronous submis\
sion\n      --status          :      : poll for th\
e status of a job\n      --polljob         :      \
: poll for the results of a job\n  -j, --jobid    \
       : str  : jobid for an asynchronous job\n  -\
O, --outfile         : str  : file name for result\
s (default is jobid)\n      --trace	        :     \
 : show SOAP messages being interchanged \n\nSynch\
ronous job:\n\n  The results/errors are returned a\
s soon as the job is finished.\n  Usage: $scriptNa\
me --email <your\\@email> [options] pdbFile [--out\
file string]\n  Returns: saves the results to disk\
\n\nAsynchronous job:\n\n  Use this if you want to\
 retrieve the results at a later time. The results\
 \n  are stored for up to 24 hours. \n  The asynch\
ronous submission mode is recommended when users a\
re submitting \n  batch jobs or large database sea\
rches	\n  Usage: $scriptName --email <your\\@email\
> --async [options] pdbFile\n  Returns: jobid\n\n \
 Use the jobid to query for the status of the job.\
 \n  Usage: $scriptName --status --jobid <jobId>\n\
  Returns: string indicating the status of the job\
:\n    DONE - job has finished\n    RUNNING - job \
is running\n    NOT_FOUND - job cannot be found\n \
   ERROR - the jobs has encountered an error\n\n  \
When done, use the jobid to retrieve the status of\
 the job. \n  Usage: $scriptName --polljob --jobid\
 <jobId> [--outfile string]\n\n[Help]\n\n  For mor\
e detailed help information refer to\n  http://www\
.ebi.ac.uk/DaliLite/\nEOF\n;\n}\n","my $WSDL = 'ht\
tp://www.ebi.ac.uk/Tools/webservices/wsdl/WSWUBlas\
t.wsdl';\n\nuse strict;\nuse SOAP::Lite;\nuse Geto\
pt::Long qw(:config no_ignore_case bundling);\nuse\
 File::Basename;\n\nmy $checkInterval = 15;\n\nmy \
$numOpts = scalar(@ARGV);\nmy ($outfile, $outforma\
t, $help, $async, $polljob, $status, $ids, $jobid,\
 $trace, $sequence);\nmy %params= ( # Defaults\n	 \
     'async' => 1, # Force into async mode\n	     \
 'exp' => 10.0, # E-value threshold\n	      'numal\
' => 50, # Maximum number of alignments\n	      's\
cores' => 100, # Maximum number of scores\n       \
     );\nGetOptions( # Map the options into variab\
les\n    \"program|p=s\"     => \\$params{program}\
, # BLAST program\n    \"database|D=s\"    => \\$p\
arams{database}, # Search database\n    \"matrix|m\
=s\"      => \\$params{matrix}, # Scoring matrix\n\
    \"exp|E=f\"         => \\$params{exp}, # E-val\
ue threshold\n    \"echofilter|e\"    => \\$params\
{echofilter}, # Display filtered sequence\n    \"f\
ilter|f=s\"      => \\$params{filter}, # Low compl\
exity filter name\n    \"alignments|b=i\"  => \\$p\
arams{numal}, # Number of alignments\n    \"scores\
|s=i\"      => \\$params{scores}, # Number of scor\
es\n    \"sensitivity|S=s\" => \\$params{sensitivi\
ty}, # Search sensitivity\n    \"sort|t=s\"	      \
=> \\$params{sort}, # Sort hits by...\n    \"stats\
|T=s\"       => \\$params{stats}, # Scoring statis\
tic to use\n    \"strand|d=s\"      => \\$params{s\
trand}, # Strand to use in DNA vs. DNA search\n   \
 \"topcombon|c=i\"   => \\$params{topcombon}, # Co\
nsistent sets of HSPs\n    \"outfile=s\"       => \
\\$outfile, # Output file\n    \"outformat|o=s\"  \
 => \\$outformat, # Output format\n    \"help|h\"	\
      => \\$help, # Usage info\n    \"async|a\"	  \
    => \\$async, # Asynchronous mode\n    \"polljo\
b\"	      => \\$polljob, # Get results\n    \"stat\
us\"	      => \\$status, # Get job status\n    \"i\
ds\"             => \\$ids, # Get ids from result\\
n    \"jobid|j=s\"       => \\$jobid, # JobId\n   \
 \"email=s\"         => \\$params{email}, # E-mail\
 address\n    \"trace\"           => \\$trace, # S\
OAP trace\n    \"sequence=s\"      => \\$sequence,\
 # Query sequence\n    );\n\nmy $scriptName = base\
name($0, ());\nif($help || $numOpts == 0) {\n    &\
usage();\n    exit(0);\n}\n\nif($trace){\n    prin\
t STDERR \"Tracing active\\n\";\n    SOAP::Lite->i\
mport(+trace => 'debug');\n}\n\nmy $soap = SOAP::L\
ite\n    ->service($WSDL)\n    ->proxy('http://loc\
alhost/',\n    #proxy => ['http' => 'http://your.p\
roxy.server/'], # HTTP proxy\n    timeout => 600, \
# HTTP connection timeout\n    )\n    ->on_fault(s\
ub { # SOAP fault handler\n        my $soap = shif\
t;\n        my $res = shift;\n        # Throw an e\
xception for all faults\n        if(ref($res) eq '\
') {\n            die($res);\n        } else {\n  \
          die($res->faultstring);\n        }\n    \
    return new SOAP::SOM;\n    }\n               )\
;\n\nif( !($polljob || $status || $ids) &&\n    !(\
 defined($ARGV[0]) || defined($sequence) )\n    ) \
{\n    print STDERR 'Error: bad option combination\
', \"\\n\";\n    &usage();\n    exit(1);\n}\nelsif\
($polljob && defined($jobid)) {\n    print \"Getti\
ng results for job $jobid\\n\";\n    getResults($j\
obid);\n}\nelsif($status && defined($jobid)) {\n  \
  print STDERR \"Getting status for job $jobid\\n\\
";\n    my $result = $soap->checkStatus($jobid);\n\
    print STDOUT \"$result\\n\";\n    if($result e\
q 'DONE') {\n	print STDERR \"To get results: $scri\
ptName --polljob --jobid $jobid\\n\";\n    }\n}  \\
nelsif($ids && defined($jobid)) {\n    print STDER\
R \"Getting ids from job $jobid\\n\";\n    getIds(\
$jobid);\n}\nelse {\n    # Prepare input data\n   \
 my $content;\n    my (@contents) = ();\n    if(-f\
 $ARGV[0] || $ARGV[0] eq '-') {	\n	$content={type=\
>'sequence',content=>read_file($ARGV[0])};	\n    }\
\n    if($sequence) {	\n	if(-f $sequence || $seque\
nce eq '-') {	\n	    $content={type=>'sequence',co\
ntent=>read_file($ARGV[0])};	\n	} else {\n	    $co\
ntent={type=>'sequence',content=>$sequence};\n	}\n\
    }\n    push @contents, $content;\n\n    # Subm\
it the job\n    my $paramsData = SOAP::Data->name(\
'params')->type(map=>\\%params);\n    my $contentD\
ata = SOAP::Data->name('content')->value(\\@conten\
ts);\n    # For SOAP::Lite 0.60 and earlier parame\
ters are passed directly\n    if($SOAP::Lite::VERS\
ION eq '0.60' || $SOAP::Lite::VERSION =~ /0\\.[1-5\
]/) {\n        $jobid = $soap->runWUBlast($paramsD\
ata, $contentData);\n    }\n    # For SOAP::Lite 0\
.69 and later parameter handling is different, so \
pass\n    # undef's for templated params, and then\
 pass the formatted args.\n    else {\n        $jo\
bid = $soap->runWUBlast(undef, undef,\n				   $par\
amsData, $contentData);\n    }\n\n    # Asynchrono\
us mode: output jobid and exit.\n    if (defined($\
async)) {\n	print STDOUT $jobid, \"\\n\";\n       \
 print STDERR \"To check status: $scriptName --sta\
tus --jobid $jobid\\n\";\n    }\n    # Synchronous\
 mode: try to get results\n    else {\n        pri\
nt STDERR \"JobId: $jobid\\n\";\n        sleep 1;\\
n        getResults($jobid);\n    }\n}\n\nsub getI\
ds($) {\n    my $jobid = shift;\n    my $results =\
 $soap->getIds($jobid);\n    for my $result (@$res\
ults){\n	print \"$result\\n\";\n    }\n}\n\nsub cl\
ientPoll($) {\n    my $jobid = shift;\n    my $res\
ult = 'PENDING';\n    # Check status and wait if n\
ot finished\n    while($result eq 'RUNNING' || $re\
sult eq 'PENDING') {\n        $result = $soap->che\
ckStatus($jobid);\n        print STDERR \"$result\\
\n\";\n        if($result eq 'RUNNING' || $result \
eq 'PENDING') {\n            # Wait before polling\
 again.\n            sleep $checkInterval;\n      \
  }\n    }\n}\n\nsub getResults($) {\n    my $jobi\
d = shift;\n    my $res;\n    # Check status, and \
wait if not finished\n    clientPoll($jobid);\n   \
 # Use JobId if output file name is not defined\n \
   unless(defined($outfile)) {\n        $outfile=$\
jobid;\n    }\n    # Get list of data types\n    m\
y $resultTypes = $soap->getResults($jobid);\n    #\
 Get the data and write it to a file\n    if(defin\
ed($outformat)) { # Specified data type\n	if($outf\
ormat eq 'xml') {$outformat = 'toolxml';}\n	if($ou\
tformat eq 'txt') {$outformat = 'tooloutput';}\n  \
      my $selResultType;\n        foreach my $resu\
ltType (@$resultTypes) {\n            if($resultTy\
pe->{type} eq $outformat) {\n                $selR\
esultType = $resultType;\n            }\n        }\
\n        $res=$soap->poll($jobid, $selResultType-\
>{type});\n	if($outfile eq '-') {\n	     write_fil\
e($outfile, $res);\n	} else {\n	    write_file($ou\
tfile.'.'.$selResultType->{ext}, $res);\n	}\n    }\
 else { # Data types available\n        # Write a \
file for each output type\n        for my $resultT\
ype (@$resultTypes){\n            #print STDERR \"\
Getting $resultType->{type}\\n\";\n            $re\
s=$soap->poll($jobid, $resultType->{type});\n	    \
if($outfile eq '-') {\n		write_file($outfile, $res\
);\n	    } else {\n		write_file($outfile.'.'.$resu\
ltType->{ext}, $res);\n	    }\n        }\n    }\n}\
\n\nsub read_file($) {\n    my $filename = shift;\\
n    my ($content, $buffer);\n    if($filename eq \
'-') {\n	while(sysread(STDIN, $buffer, 1024)) {\n	\
    $content .= $buffer;\n	}\n    }\n    else { # \
File\n	open(FILE, $filename) or die \"Error: unabl\
e to open input file\";\n	while(sysread(FILE, $buf\
fer, 1024)) {\n	    $content .= $buffer;\n	}\n	clo\
se(FILE);\n    }\n    return $content;\n}\n\nsub w\
rite_file($$) {\n    my ($filename, $data) = @_;\n\
    print STDERR 'Creating result file: ' . $filen\
ame . \"\\n\";\n    if($filename eq '-') {\n	print\
 STDOUT $data;\n    }\n    else {\n	open(FILE, \">\
$filename\") or die \"Error: unable to open output\
 file\";\n	syswrite(FILE, $data);\n	close(FILE);\n\
    }\n}\n\nsub usage {\n    print STDERR <<EOF\nW\
U-BLAST\n========\n\nRapid sequence database searc\
h programs utilizing the BLAST algorithm.\n   \n[R\
equired]\n\n      --email       : str  : user emai\
l address \n  -p, --program	    : str  : BLAST pro\
gram to use: blastn, blastp, blastx, \n           \
                  tblastn or tblastx\n  -D, --data\
base    : str  : database to search\n  seqFile    \
       : file : query sequence data file (\"-\" fo\
r STDIN)\n\n[Optional]\n\n  -m, --matrix	    : str\
  : scoring matrix\n  -E, --exp	    : real : 0<E<=\
 1000. Statistical significance threshold\n       \
                      for reporting database seque\
nce matches.\n  -e, --echofilter  :      : display\
 the filtered query sequence in the output\n  -f, \
--filter	    : str  : activates filtering of the q\
uery sequence\n  -b, --alignments  : int  : number\
 of alignments to be reported\n  -s, --scores	    \
: int  : number of scores to be reported\n  -S, --\
sensitivity : str  :\n  -t, --sort	    : str  :\n \
 -T, --stats       : str  :\n  -d, --strand      :\
 str  : DNA strand to search with in DNA vs. DNA s\
earches \n  -c, --topcombon   :      :\n\n[General\
]	\n\n  -h, --help       :      : prints this help\
 text\n  -a, --async      :      : forces to make \
an asynchronous query\n      --status     :      :\
 poll for the status of a job\n      --polljob    \
:      : poll for the results of a job\n  -j, --jo\
bid      : str  : jobid that was returned when an \
asynchronous job \n                            was\
 submitted.\n  -O, --outfile    : str  : name of t\
he file results should be written to \n           \
                 (default is based on the jobid; \\
"-\" for STDOUT)\n  -o, --outformat  : str  : txt \
or xml output (no file is written)\n      --trace	\
   :      : show SOAP messages being interchanged \
\n\nSynchronous job:\n\n  The results/errors are r\
eturned as soon as the job is finished.\n  Usage: \
$scriptName --email <your\\@email> [options...] se\
qFile\n  Returns: saves the results to disk\n\nAsy\
nchronous job:\n\n  Use this if you want to retrie\
ve the results at a later time. The results \n  ar\
e stored for up to 24 hours. \n  The asynchronous \
submission mode is recommended when users are subm\
itting \n  batch jobs or large database searches	\\
n  Usage: $scriptName --async --email <your\\@emai\
l> [options...] seqFile\n  Returns : jobid\n\n  Us\
e the jobid to query for the status of the job. \n\
  Usage: $scriptName --status --jobid <jobId>\n  R\
eturns : string indicating the status of the job:\\
n    DONE - job has finished\n    RUNNING - job is\
 running\n    NOT_FOUND - job cannot be found\n   \
 ERROR - the jobs has encountered an error\n\n  Wh\
en done, use the jobid to retrieve the status of t\
he job. \n  Usage: $scriptName --polljob --jobid <\
jobId> [--outfile string]\n  Returns: saves the re\
sults to disk\n\n[Help]\n\nFor more detailed help \
information refer to \nhttp://www.ebi.ac.uk/blast2\
/WU-Blast2_Help_frame.html\n \nEOF\n;\n}\n","\nmy \
$WSDL = 'http://www.ebi.ac.uk/Tools/webservices/ws\
dl/WSBlastpgp.wsdl';\n\nuse SOAP::Lite;\nuse Getop\
t::Long qw(:config no_ignore_case bundling);\nuse \
File::Basename;\n\nmy $checkInterval = 15;\n\nmy %\
params=(\n	    'async' => '1', # Use async mode an\
d simulate sync mode in client\n	    );\nGetOption\
s(\n    \"mode=s\"           => \\$params{mode}, #\
 Search mode: PSI-Blast or PHI-Blast\n    \"databa\
se|d=s\"     => \\$params{database}, # Database to\
 search\n    \"matrix|M=s\"       => \\$params{mat\
rix},# Scoring maxtrix\n    \"exp|e=f\"          =\
> \\$params{exp}, # E-value\n    \"expmulti|h=f\" \
    => \\$params{expmulti}, # E-value\n    \"filte\
r|F=s\"       => \\$params{filter}, # Low complexi\
ty filter\n    \"dropoff|X=i\"      => \\$params{d\
ropoff}, # Dropoff score\n    \"finaldropoff|Z=i\"\
 => \\$params{finaldropoff}, # Final dropoff score\
\n    \"scores|v=i\"       => \\$params{scores}, #\
 Max number of scores\n    \"align=i\"          =>\
 \\$params{align}, # Alignment view\n    \"startre\
gion|S=i\"  => \\$params{startregion}, # Start of \
region in query\n    \"endregion|H=i\"    => \\$pa\
rams{endregion}, # End of region in query\n    \"m\
axpasses|j=i\"    => \\$params{maxpasses}, # Numbe\
r of PSI iterations\n    \"opengap|G=i\"      => \\
\$params{opengap}, # Gap open penalty\n    \"exten\
dgap|E=i\"    => \\$params{extendgap}, # Gap exten\
sion penalty\n    \"pattern=s\"        => \\$param\
s{pattern}, # PHI-BLAST pattern\n    \"usagemode|p\
=s\"    => \\$params{usagemode}, # PHI-BLAST progr\
am\n    \"appxml=s\"         => \\$params{appxml},\
 # Application XML\n    \"sequence=s\"       => \\\
$sequence, # Query sequence\n    \"help\"	       =\
> \\$help, # Usage info\n    \"polljob\"	       =>\
 \\$polljob, # Get results\n    \"status\"	       \
=> \\$status, # Get status\n    \"ids\"      	    \
   => \\$ids, # Get ids from result\n    \"jobid=s\
\"          => \\$jobid, # JobId\n    \"outfile=s\\
"        => \\$outfile, # Output filename\n    \"o\
utformat|o=s\"    => \\$outformat, # Output file f\
ormat\n    \"async|a\"	       => \\$async, # Async\
 submission\n    \"email=s\"          => \\$params\
{email}, # User e-mail address\n    \"trace\"     \
       => \\$trace, # Show SOAP messages\n    );\n\
\nmy $scriptName = basename($0, ());\nif($help) {\\
n    &usage();\n    exit(0);\n}\n\nif ($trace){\n \
   print \"Tracing active\\n\";\n    SOAP::Lite->i\
mport(+trace => 'debug');\n}\n\nmy $soap = SOAP::L\
ite\n    ->service($WSDL)\n    ->on_fault(sub {\n \
       my $soap = shift;\n        my $res = shift;\
\n        # Throw an exception for all faults\n   \
     if(ref($res) eq '') {\n            die($res);\
\n        } else {\n            die($res->faultstr\
ing);\n        }\n        return new SOAP::SOM;\n \
   }\n               );\n\nif( !($polljob || $stat\
us || $ids) &&\n    !( (defined($ARGV[0]) && -f $A\
RGV[0]) || defined($sequence) )\n    ) {\n    prin\
t STDERR 'Error: bad option combination', \"\\n\";\
\n    &usage();\n    exit(1);\n}\nelsif($polljob &\
& defined($jobid)) {\n    print \"Getting results \
for job $jobid\\n\";\n    getResults($jobid);\n}\n\
elsif($status && defined($jobid)) {\n    print STD\
ERR \"Getting status for job $jobid\\n\";\n    my \
$result = $soap->checkStatus($jobid);\n    print S\
TDOUT $result, \"\\n\";\n    if($result eq 'DONE')\
 {\n	print STDERR \"To get results: $scriptName --\
polljob --jobid $jobid\\n\";\n    }\n}  \nelsif($i\
ds && defined($jobid)) {\n    print STDERR \"Getti\
ng ids from job $jobid\\n\";\n    getIds($jobid);\\
n}\nelse {\n    if(-f $ARGV[0]) {	\n	$content={typ\
e=>'sequence', content=>read_file($ARGV[0])};	\n  \
  }\n    if($sequence) {	\n	if(-f $sequence) {\n	 \
   $content={type=>'sequence', content=>read_file(\
$sequence)};	\n	} else {\n	    $content={type=>'se\
quence', content=>$sequence};\n	}\n    }\n    push\
 @content, $content;\n\n    my $jobid;\n    my $pa\
ramsData = SOAP::Data->name('params')->type(map=>\\
\%params);\n    my $contentData = SOAP::Data->name\
('content')->value(\\@content);\n    # For SOAP::L\
ite 0.60 and earlier parameters are passed directl\
y\n    if($SOAP::Lite::VERSION eq '0.60' || $SOAP:\
:Lite::VERSION =~ /0\\.[1-5]/) {\n        $jobid =\
 $soap->runBlastpgp($paramsData, $contentData);\n \
   }\n    # For SOAP::Lite 0.69 and later paramete\
r handling is different, so pass\n    # undef's fo\
r templated params, and then pass the formatted ar\
gs.\n    else {\n        $jobid = $soap->runBlastp\
gp(undef, undef,\n				    $paramsData, $contentDat\
a);\n    }\n\n    if (defined($async)) {\n	print S\
TDOUT $jobid, \"\\n\";\n        print STDERR \"To \
check status: $scriptName --status --jobid $jobid\\
\n\";\n    } else { # Synchronous mode\n        pr\
int STDERR \"JobId: $jobid\\n\";\n        sleep 1;\
\n        getResults($jobid);\n    }\n}\n\nsub get\
Ids($) {\n    $jobid = shift;\n    my $results = $\
soap->getIds($jobid);\n    for $result (@$results)\
{\n	print \"$result\\n\";\n    }\n}\n\nsub clientP\
oll($) {\n    my $jobid = shift;\n    my $result =\
 'PENDING';\n    # Check status and wait if not fi\
nished\n    #print STDERR \"Checking status: $jobi\
d\\n\";\n    while($result eq 'RUNNING' || $result\
 eq 'PENDING') {\n        $result = $soap->checkSt\
atus($jobid);\n        print STDERR \"$result\\n\"\
;\n        if($result eq 'RUNNING' || $result eq '\
PENDING') {\n            # Wait before polling aga\
in.\n            sleep $checkInterval;\n        }\\
n    }\n}\n\nsub getResults($) {\n    $jobid = shi\
ft;\n    # Check status, and wait if not finished\\
n    clientPoll($jobid);\n    # Use JobId if outpu\
t file name is not defined\n    unless(defined($ou\
tfile)) {\n        $outfile=$jobid;\n    }\n    # \
Get list of data types\n    my $resultTypes = $soa\
p->getResults($jobid);\n    # Get the data and wri\
te it to a file\n    if(defined($outformat)) { # S\
pecified data type\n        my $selResultType;\n  \
      foreach my $resultType (@$resultTypes) {\n  \
          if($resultType->{type} eq $outformat) {\\
n                $selResultType = $resultType;\n  \
          }\n        }\n        $res=$soap->poll($\
jobid, $selResultType->{type});\n        write_fil\
e($outfile.'.'.$selResultType->{ext}, $res);\n    \
} else { # Data types available\n        # Write a\
 file for each output type\n        for my $result\
Type (@$resultTypes){\n            #print \"Gettin\
g $resultType->{type}\\n\";\n            $res=$soa\
p->poll($jobid, $resultType->{type});\n           \
 write_file($outfile.'.'.$resultType->{ext}, $res)\
;\n        }\n    }\n}\n\nsub read_file($) {\n    \
my $filename = shift;\n    open(FILE, $filename);\\
n    my $content;\n    my $buffer;\n    while(sysr\
ead(FILE, $buffer, 1024)) {\n	$content.= $buffer;\\
n    }\n    close(FILE);  \n    return $content;\n\
}\n\nsub write_file($$) {\n    my ($tmp,$entity) =\
 @_;\n    print STDERR \"Creating result file: \".\
$tmp.\"\\n\";\n    unless(open (FILE, \">$tmp\")) \
{\n	return 0;\n    }\n    syswrite(FILE, $entity);\
\n    close (FILE);\n    return 1;\n}\n\nsub usage\
 {\n    print STDERR <<EOF\nBlastpgp\n========\n  \
 \nThe blastpgp program implements the PSI-BLAST a\
nd PHI-BLAST variations\nof NCBI BLAST.\n\nFor mor\
e detailed help information refer to\nhttp://www.e\
bi.ac.uk/blastpgp/blastpsi_help_frame.html\n \nBla\
stpgp specific options:\n\n[Required]\n\n      --m\
ode            : str  : search mode to use: PSI-Bl\
ast or PHI-Blast\n  -d, --database        : str  :\
 protein database to search\n  seqFile            \
   : file : query sequence\n\n[Optional]\n\n  -M, \
--matrix          : str  : scoring matrix\n  -e, -\
-exp             : real : Expectation value\n  -h,\
 --expmulti        : real : threshold (multipass m\
odel)\n  -F, --filter          : str  : filter que\
ry sequence with SEG [T,F]\n  -m, --align         \
  : int  : alignment view option:\n               \
                  0 - pairwise, 1 - M/S identities\
,\n                                 2 - M/S non-id\
entities, 3 - Flat identities,\n                  \
               4 - Flat non-identities\n  -G, --op\
engap         : int  : cost to open a gap\n  -E, -\
-extendgap       : int  : cost to extend a gap\n  \
-g, --gapalign        : str  : Gapped [T,F]\n  -v,\
 --scores          : int  : number of scores to be\
 reported\n  -j, --maxpasses       : int  : number\
 of iterations\n  -X, --dropoff         : int  : D\
ropoff score\n  -Z, --finaldropoff    : int  : Dro\
poff for final alignment\n  -S, --startregion     \
: int  : Start of required region in query\n  -H, \
--endregion       : int  : End of required region \
in query\n  -k, --pattern         : str  : Hit Fil\
e (PHI-BLAST only)\n  -p, --usagemode       : str \
 : Program option (PHI-BLAST only):\n             \
                    blastpgp, patseedp, seedp\n\n[\
General]\n\n      --help            :      : print\
s this help text\n  -a, --async           :      :\
 forces to make an asynchronous query\n      --sta\
tus          :      : poll for the status of a job\
\n      --polljob         :      : poll for the re\
sults of a job\n      --jobid           : str  : j\
obid of an asynchronous job\n      --ids          \
   :      : get hit identifiers for result \n  -O,\
 --outfile         : str  : name of the file resul\
ts should be written to\n                         \
        (default is based on the jobid)\n  -o, --o\
utformat       : str  : txt or xml output (no file\
 is written)\n      --trace           :      : sho\
w SOAP messages being interchanged\n\nSynchronous \
job:\n\n  The results/errors are returned as soon \
as the job is finished.\n  Usage: blastpgp.pl --em\
ail <your@email> [options...] seqfile\n  Returns: \
saves the results to disk\n\nAsynchronous job:\n\n\
  Use this if you want to retrieve the results at \
a later time. The results\n  are stored for up to \
24 hours.\n  The asynchronous submission mode is r\
ecommended when users are submitting\n  batch jobs\
 or large database searches\n  Usage: blastpgp.pl \
--email <your@email> --async [options...] seqFile\\
n  Returns: jobid\n\n  Use the jobid to query for \
the status of the job.\n  Usage: blastpgp.pl --sta\
tus --jobid <jobId>\n  Returns: string indicating \
the status of the job\n    DONE - job has finished\
\n    RUNNING - job is running\n    NOT_FOUND - jo\
b cannot be found\n    ERROR - the jobs has encoun\
tered an error\n\n  When done, use the jobid to re\
trieve the results of the job.\n  Usage: blastpgp.\
pl --polljob --jobid <jobId> [--outfile <fileName>\
]\n  Returns: saves the results to disk\nEOF\n;\n}\
\n","\n=head1 NAME\n\nncbiblast_lwp.pl\n\n=head1 D\
ESCRIPTION\n\nNCBI BLAST REST web service Perl cli\
ent using L<LWP>.\n\nTested with:\n\n=over\n\n=ite\
m *\nL<LWP> 5.79, L<XML::Simple> 2.12 and Perl 5.8\
.3\n\n=item *\nL<LWP> 5.805, L<XML::Simple> 2.14 a\
nd Perl 5.8.7\n\n=item *\nL<LWP> 5.820, L<XML::Sim\
ple> 2.18 and Perl 5.10.0 (Ubuntu 9.04)\n\n=back\n\
\nFor further information see:\n\n=over\n\n=item *\
\nL<http://www.ebi.ac.uk/Tools/webservices/service\
s/sss/ncbi_blast_rest>\n\n=item *\nL<http://www.eb\
i.ac.uk/Tools/webservices/tutorials/perl>\n\n=back\
\n\n=head1 VERSION\n\n$Id: ncbiblast_lwp.pl 1317 2\
009-09-03 15:44:11Z hpm $\n\n=cut\n\nuse strict;\n\
use warnings;\n\nuse English;\nuse LWP;\nuse XML::\
Simple;\nuse Getopt::Long qw(:config no_ignore_cas\
e bundling);\nuse File::Basename;\nuse Data::Dumpe\
r;\n\nmy $baseUrl = 'http://www.ebi.ac.uk/Tools/se\
rvices/rest/ncbiblast';\n\nmy $checkInterval = 3;\\
n\nmy $outputLevel = 1;\n\nmy $numOpts = scalar(@A\
RGV);\nmy %params = ( 'debugLevel' => 0 );\n\nmy %\
tool_params = ();\nGetOptions(\n\n	# Tool specific\
 options\n	'program|p=s'  => \\$tool_params{'progr\
am'},   # blastp, blastn, blastx, etc.\n	'database\
|D=s' => \\$params{'database'},       # Database(s\
) to search\n	'matrix|m=s'   => \\$tool_params{'ma\
trix'},    # Scoring martix to use\n	'exp|E=f'    \
  => \\$tool_params{'exp'},       # E-value thresh\
old\n	'filter|f=s'   => \\$tool_params{'filter'}, \
   # Low complexity filter\n	'align|A=i'    => \\$\
tool_params{'align'},     # Pairwise alignment for\
mat\n	'scores|s=i'   => \\$tool_params{'scores'}, \
   # Number of scores\n	'alignments|n=i' => \\$too\
l_params{'alignments'},   # Number of alignments\n\
	'dropoff|d=i'    => \\$tool_params{'dropoff'},   \
   # Dropoff score\n	'match_scores=s' => \\$tool_p\
arams{'match_scores'}, # Match/missmatch scores\n	\
'match|u=i'      => \\$params{'match'},           \
  # Match score\n	'mismatch|v=i'   => \\$params{'m\
ismatch'},          # Mismatch score\n	'gapopen|o=\
i'    => \\$tool_params{'gapopen'},      # Open ga\
p penalty\n	'gapext|x=i'     => \\$tool_params{'ga\
pext'},       # Gap extension penality\n	'gapalign\
|g'     => \\$tool_params{'gapalign'},     # Optim\
ise gap alignments\n	'stype=s' => \\$tool_params{'\
stype'},    # Sequence type\n	'seqrange=s' => \\$t\
ool_params{'seqrange'},    # Query subsequence\n	'\
sequence=s' => \\$params{'sequence'},         # Qu\
ery sequence\n	'multifasta' => \\$params{'multifas\
ta'},       # Multiple fasta input\n\n	# Compatabi\
lity options, old command-line\n	'numal|n=i'     =\
> \\$params{'numal'},        # Number of alignment\
s\n	'opengap|o=i'   => \\$params{'opengap'},      \
# Open gap penalty\n	'extendgap|x=i' => \\$params{\
'extendgap'},    # Gap extension penality\n	\n	# G\
eneric options\n	'email=s'       => \\$params{'ema\
il'},          # User e-mail address\n	'title=s'  \
     => \\$params{'title'},          # Job title\n\
	'outfile=s'     => \\$params{'outfile'},        #\
 Output file name\n	'outformat=s'   => \\$params{'\
outformat'},      # Output file type\n	'jobid=s'  \
     => \\$params{'jobid'},          # JobId\n	'he\
lp|h'        => \\$params{'help'},           # Usa\
ge help\n	'async'         => \\$params{'async'},  \
        # Asynchronous submission\n	'polljob'     \
  => \\$params{'polljob'},        # Get results\n	\
'resultTypes'   => \\$params{'resultTypes'},    # \
Get result types\n	'status'        => \\$params{'s\
tatus'},         # Get status\n	'params'        =>\
 \\$params{'params'},         # List input paramet\
ers\n	'paramDetail=s' => \\$params{'paramDetail'},\
    # Get details for parameter\n	'quiet'         \
=> \\$params{'quiet'},          # Decrease output \
level\n	'verbose'       => \\$params{'verbose'},  \
      # Increase output level\n	'debugLevel=i'  =>\
 \\$params{'debugLevel'},     # Debug output level\
\n	'baseUrl=s'     => \\$baseUrl,                 \
 # Base URL for service.\n);\nif ( $params{'verbos\
e'} ) { $outputLevel++ }\nif ( $params{'$quiet'} )\
  { $outputLevel-- }\n\n&print_debug_message( 'MAI\
N', 'LWP::VERSION: ' . $LWP::VERSION,\n	1 );\n\n&p\
rint_debug_message( 'MAIN', \"params:\\n\" . Dumpe\
r( \\%params ),           11 );\n&print_debug_mess\
age( 'MAIN', \"tool_params:\\n\" . Dumper( \\%tool\
_params ), 11 );\n\nmy $scriptName = basename( $0,\
 () );\n\nif ( $params{'help'} || $numOpts == 0 ) \
{\n	&usage();\n	exit(0);\n}\n\n&print_debug_messag\
e( 'MAIN', 'baseUrl: ' . $baseUrl, 1 );\n\nif (\n	\
!(\n		   $params{'polljob'}\n		|| $params{'resultT\
ypes'}\n		|| $params{'status'}\n		|| $params{'para\
ms'}\n		|| $params{'paramDetail'}\n	)\n	&& !( defi\
ned( $ARGV[0] ) || defined( $params{'sequence'} ) \
)\n  )\n{\n\n	# Bad argument combination, so print\
 error message and usage\n	print STDERR 'Error: ba\
d option combination', \"\\n\";\n	&usage();\n	exit\
(1);\n}\n\nelsif ( $params{'params'} ) {\n	&print_\
tool_params();\n}\n\nelsif ( $params{'paramDetail'\
} ) {\n	&print_param_details( $params{'paramDetail\
'} );\n}\n\nelsif ( $params{'status'} && defined( \
$params{'jobid'} ) ) {\n	&print_job_status( $param\
s{'jobid'} );\n}\n\nelsif ( $params{'resultTypes'}\
 && defined( $params{'jobid'} ) ) {\n	&print_resul\
t_types( $params{'jobid'} );\n}\n\nelsif ( $params\
{'polljob'} && defined( $params{'jobid'} ) ) {\n	&\
get_results( $params{'jobid'} );\n}\n\nelse {\n\n	\
# Multiple input sequence mode, assume fasta forma\
t.\n	if ( $params{'multifasta'} ) {\n		&multi_subm\
it_job();\n	}\n\n	# Entry identifier list file.\n	\
elsif (( defined( $params{'sequence'} ) && $params\
{'sequence'} =~ m/^\\@/ )\n		|| ( defined( $ARGV[0\
] ) && $ARGV[0] =~ m/^\\@/ ) )\n	{\n		my $list_fil\
ename = $params{'sequence'} || $ARGV[0];\n		$list_\
filename =~ s/^\\@//;\n		&list_file_submit_job($li\
st_filename);\n	}\n\n	# Default: single sequence/i\
dentifier.\n	else {\n\n		# Load the sequence data \
and submit.\n		&submit_job( &load_data() );\n	}\n}\
\n\n=head1 FUNCTIONS\n\n=cut\n\n\n=head2 rest_requ\
est()\n\nPerform a REST request.\n\n  my $response\
_str = &rest_request($url);\n\n=cut\n\nsub rest_re\
quest {\n	print_debug_message( 'rest_request', 'Be\
gin', 11 );\n	my $requestUrl = shift;\n	print_debu\
g_message( 'rest_request', 'URL: ' . $requestUrl, \
11 );\n\n	# Create a user agent\n	my $ua = LWP::Us\
erAgent->new();\n	'$Revision: 1317 $' =~ m/(\\d+)/\
;\n	$ua->agent(\"EBI-Sample-Client/$1 ($scriptName\
; $OSNAME) \" . $ua->agent());\n	$ua->env_proxy;\n\
\n	# Perform the request\n	my $response = $ua->get\
($requestUrl);\n	print_debug_message( 'rest_reques\
t', 'HTTP status: ' . $response->code,\n		11 );\n\\
n	# Check for HTTP error codes\n	if ( $response->i\
s_error ) {\n		$response->content() =~ m/<h1>([^<]\
+)<\\/h1>/;\n		die 'http status: ' . $response->co\
de . ' ' . $response->message . '  ' . $1;\n	}\n	p\
rint_debug_message( 'rest_request', 'End', 11 );\n\
\n	# Return the response data\n	return $response->\
content();\n}\n\n=head2 rest_get_parameters()\n\nG\
et list of tool parameter names.\n\n  my (@param_l\
ist) = &rest_get_parameters();\n\n=cut\n\nsub rest\
_get_parameters {\n	print_debug_message( 'rest_get\
_parameters', 'Begin', 1 );\n	my $url             \
   = $baseUrl . '/parameters/';\n	my $param_list_x\
ml_str = rest_request($url);\n	my $param_list_xml \
    = XMLin($param_list_xml_str);\n	my (@param_lis\
t)       = @{ $param_list_xml->{'id'} };\n	print_d\
ebug_message( 'rest_get_parameters', 'End', 1 );\n\
	return (@param_list);\n}\n\n=head2 rest_get_param\
eter_details()\n\nGet details of a tool parameter.\
\n\n  my $paramDetail = &rest_get_parameter_detail\
s($param_name);\n\n=cut\n\nsub rest_get_parameter_\
details {\n	print_debug_message( 'rest_get_paramet\
er_details', 'Begin', 1 );\n	my $parameterId = shi\
ft;\n	print_debug_message( 'rest_get_parameter_det\
ails',\n		'parameterId: ' . $parameterId, 1 );\n	m\
y $url                  = $baseUrl . '/parameterde\
tails/' . $parameterId;\n	my $param_detail_xml_str\
 = rest_request($url);\n	my $param_detail_xml     \
= XMLin($param_detail_xml_str);\n	print_debug_mess\
age( 'rest_get_parameter_details', 'End', 1 );\n	r\
eturn ($param_detail_xml);\n}\n\n=head2 rest_run()\
\n\nSubmit a job.\n\n  my $job_id = &rest_run($ema\
il, $title, \\%params );\n\n=cut\n\nsub rest_run {\
\n	print_debug_message( 'rest_run', 'Begin', 1 );\\
n	my $email  = shift;\n	my $title  = shift;\n	my $\
params = shift;\n	print_debug_message( 'rest_run',\
 'email: ' . $email, 1 );\n	if ( defined($title) )\
 {\n		print_debug_message( 'rest_run', 'title: ' .\
 $title, 1 );\n	}\n	print_debug_message( 'rest_run\
', 'params: ' . Dumper($params), 1 );\n\n	# User a\
gent to perform http requests\n	my $ua = LWP::User\
Agent->new();\n	$ua->env_proxy;\n\n	# Clean up par\
ameters\n	my (%tmp_params) = %{$params};\n	$tmp_pa\
rams{'email'} = $email;\n	$tmp_params{'title'} = $\
title;\n	foreach my $param_name ( keys(%tmp_params\
) ) {\n		if ( !defined( $tmp_params{$param_name} )\
 ) {\n			delete $tmp_params{$param_name};\n		}\n	}\
\n\n	# Submit the job as a POST\n	my $url = $baseU\
rl . '/run';\n	my $response = $ua->post( $url, \\%\
tmp_params );\n	print_debug_message( 'rest_run', '\
HTTP status: ' . $response->code, 11 );\n	print_de\
bug_message( 'rest_run',\n		'request: ' . $respons\
e->request()->content(), 11 );\n\n	# Check for HTT\
P error codes\n	if ( $response->is_error ) {\n		$r\
esponse->content() =~ m/<h1>([^<]+)<\\/h1>/;\n		di\
e 'http status: ' . $response->code . ' ' . $respo\
nse->message . '  ' . $1;\n	}\n\n	# The job id is \
returned\n	my $job_id = $response->content();\n	pr\
int_debug_message( 'rest_run', 'End', 1 );\n	retur\
n $job_id;\n}\n\n=head2 rest_get_status()\n\nCheck\
 the status of a job.\n\n  my $status = &rest_get_\
status($job_id);\n\n=cut\n\nsub rest_get_status {\\
n	print_debug_message( 'rest_get_status', 'Begin',\
 1 );\n	my $job_id = shift;\n	print_debug_message(\
 'rest_get_status', 'jobid: ' . $job_id, 2 );\n	my\
 $status_str = 'UNKNOWN';\n	my $url        = $base\
Url . '/status/' . $job_id;\n	$status_str = &rest_\
request($url);\n	print_debug_message( 'rest_get_st\
atus', 'status_str: ' . $status_str, 2 );\n	print_\
debug_message( 'rest_get_status', 'End', 1 );\n	re\
turn $status_str;\n}\n\n=head2 rest_get_result_typ\
es()\n\nGet list of result types for finished job.\
\n\n  my (@result_types) = &rest_get_result_types(\
$job_id);\n\n=cut\n\nsub rest_get_result_types {\n\
	print_debug_message( 'rest_get_result_types', 'Be\
gin', 1 );\n	my $job_id = shift;\n	print_debug_mes\
sage( 'rest_get_result_types', 'jobid: ' . $job_id\
, 2 );\n	my (@resultTypes);\n	my $url             \
         = $baseUrl . '/resulttypes/' . $job_id;\n\
	my $result_type_list_xml_str = &rest_request($url\
);\n	my $result_type_list_xml     = XMLin($result_\
type_list_xml_str);\n	(@resultTypes) = @{ $result_\
type_list_xml->{'type'} };\n	print_debug_message( \
'rest_get_result_types',\n		scalar(@resultTypes) .\
 ' result types', 2 );\n	print_debug_message( 'res\
t_get_result_types', 'End', 1 );\n	return (@result\
Types);\n}\n\n=head2 rest_get_result()\n\nGet resu\
lt data of a specified type for a finished job.\n\\
n  my $result = rest_get_result($job_id, $result_t\
ype);\n\n=cut\n\nsub rest_get_result {\n	print_deb\
ug_message( 'rest_get_result', 'Begin', 1 );\n	my \
$job_id = shift;\n	my $type   = shift;\n	print_deb\
ug_message( 'rest_get_result', 'jobid: ' . $job_id\
, 1 );\n	print_debug_message( 'rest_get_result', '\
type: ' . $type,    1 );\n	my $url    = $baseUrl .\
 '/result/' . $job_id . '/' . $type;\n	my $result \
= &rest_request($url);\n	print_debug_message( 'res\
t_get_result', length($result) . ' characters',\n	\
	1 );\n	print_debug_message( 'rest_get_result', 'E\
nd', 1 );\n	return $result;\n}\n\n\n=head2 print_d\
ebug_message()\n\nPrint debug message at specified\
 debug level.\n\n  &print_debug_message($method_na\
me, $message, $level);\n\n=cut\n\nsub print_debug_\
message {\n	my $function_name = shift;\n	my $messa\
ge       = shift;\n	my $level         = shift;\n	i\
f ( $level <= $params{'debugLevel'} ) {\n		print S\
TDERR '[', $function_name, '()] ', $message, \"\\n\
\";\n	}\n}\n\n=head2 print_tool_params()\n\nPrint \
list of tool parameters.\n\n  &print_tool_params()\
;\n\n=cut\n\nsub print_tool_params {\n	print_debug\
_message( 'print_tool_params', 'Begin', 1 );\n	my \
(@param_list) = &rest_get_parameters();\n	foreach \
my $param ( sort(@param_list) ) {\n		print $param,\
 \"\\n\";\n	}\n	print_debug_message( 'print_tool_p\
arams', 'End', 1 );\n}\n\n=head2 print_param_detai\
ls()\n\nPrint details of a tool parameter.\n\n  &p\
rint_param_details($param_name);\n\n=cut\n\nsub pr\
int_param_details {\n	print_debug_message( 'print_\
param_details', 'Begin', 1 );\n	my $paramName = sh\
ift;\n	print_debug_message( 'print_param_details',\
 'paramName: ' . $paramName, 2 );\n	my $paramDetai\
l = &rest_get_parameter_details($paramName);\n	pri\
nt $paramDetail->{'name'}, \"\\t\", $paramDetail->\
{'type'}, \"\\n\";\n	print $paramDetail->{'descrip\
tion'}, \"\\n\";\n	foreach my $value ( @{ $paramDe\
tail->{'values'}->{'value'} } ) {\n		print $value-\
>{'value'};\n		if ( $value->{'defaultValue'} eq 't\
rue' ) {\n			print \"\\t\", 'default';\n		}\n		pri\
nt \"\\n\";\n		print \"\\t\", $value->{'label'}, \\
"\\n\";\n	}\n	print_debug_message( 'print_param_de\
tails', 'End', 1 );\n}\n\n=head2 print_job_status(\
)\n\nPrint status of a job.\n\n  &print_job_status\
($job_id);\n\n=cut\n\nsub print_job_status {\n	pri\
nt_debug_message( 'print_job_status', 'Begin', 1 )\
;\n	my $jobid = shift;\n	print_debug_message( 'pri\
nt_job_status', 'jobid: ' . $jobid, 1 );\n	if ( $o\
utputLevel > 0 ) {\n		print STDERR 'Getting status\
 for job ', $jobid, \"\\n\";\n	}\n	my $result = &r\
est_get_status($jobid);\n	print \"$result\\n\";\n	\
if ( $result eq 'FINISHED' && $outputLevel > 0 ) {\
\n		print STDERR \"To get results: $scriptName --p\
olljob --jobid \" . $jobid\n		  . \"\\n\";\n	}\n	p\
rint_debug_message( 'print_job_status', 'End', 1 )\
;\n}\n\n=head2 print_result_types()\n\nPrint avail\
able result types for a job.\n\n  &print_result_ty\
pes($job_id);\n\n=cut\n\nsub print_result_types {\\
n	print_debug_message( 'result_types', 'Begin', 1 \
);\n	my $jobid = shift;\n	print_debug_message( 're\
sult_types', 'jobid: ' . $jobid, 1 );\n	if ( $outp\
utLevel > 0 ) {\n		print STDERR 'Getting result ty\
pes for job ', $jobid, \"\\n\";\n	}\n	my $status =\
 &rest_get_status($jobid);\n	if ( $status eq 'PEND\
ING' || $status eq 'RUNNING' ) {\n		print STDERR '\
Error: Job status is ', $status,\n		  '. To get re\
sult types the job must be finished.', \"\\n\";\n	\
}\n	else {\n		my (@resultTypes) = &rest_get_result\
_types($jobid);\n		if ( $outputLevel > 0 ) {\n			p\
rint STDOUT 'Available result types:', \"\\n\";\n	\
	}\n		foreach my $resultType (@resultTypes) {\n			\
print STDOUT $resultType->{'identifier'}, \"\\n\";\
\n			if ( defined( $resultType->{'label'} ) ) {\n	\
			print STDOUT \"\\t\", $resultType->{'label'}, \\
"\\n\";\n			}\n			if ( defined( $resultType->{'des\
cription'} ) ) {\n				print STDOUT \"\\t\", $resul\
tType->{'description'}, \"\\n\";\n			}\n			if ( de\
fined( $resultType->{'mediaType'} ) ) {\n				print\
 STDOUT \"\\t\", $resultType->{'mediaType'}, \"\\n\
\";\n			}\n			if ( defined( $resultType->{'fileSuf\
fix'} ) ) {\n				print STDOUT \"\\t\", $resultType\
->{'fileSuffix'}, \"\\n\";\n			}\n		}\n		if ( $sta\
tus eq 'FINISHED' && $outputLevel > 0 ) {\n			prin\
t STDERR \"\\n\", 'To get results:', \"\\n\",\n			\
  \"  $scriptName --polljob --jobid \" . $params{'\
jobid'} . \"\\n\",\n			  \"  $scriptName --polljob\
 --outformat <type> --jobid \"\n			  . $params{'jo\
bid'} . \"\\n\";\n		}\n	}\n	print_debug_message( '\
result_types', 'End', 1 );\n}\n\n=head2 submit_job\
()\n\nSubmit a job to the service.\n\n  &submit_jo\
b($seq);\n\n=cut\n\nsub submit_job {\n	print_debug\
_message( 'submit_job', 'Begin', 1 );\n\n	# Set in\
put sequence\n	$tool_params{'sequence'} = shift;\n\
\n	# Load parameters\n	&load_params();\n\n	# Submi\
t the job\n	my $jobid = &rest_run( $params{'email'\
}, $params{'title'}, \\%tool_params );\n\n	# Simul\
ate sync/async mode\n	if ( defined( $params{'async\
'} ) ) {\n		print STDOUT $jobid, \"\\n\";\n		if ( \
$outputLevel > 0 ) {\n			print STDERR\n			  \"To c\
heck status: $scriptName --status --jobid $jobid\\\
n\";\n		}\n	}\n	else {\n		if ( $outputLevel > 0 ) \
{\n			print STDERR \"JobId: $jobid\\n\";\n		}\n		s\
leep 1;\n		&get_results($jobid);\n	}\n	print_debug\
_message( 'submit_job', 'End', 1 );\n}\n\n=head2 m\
ulti_submit_job()\n\nSubmit multiple jobs assuming\
 input is a collection of fasta formatted sequence\
s.\n\n  &multi_submit_job();\n\n=cut\n\nsub multi_\
submit_job {\n	print_debug_message( 'multi_submit_\
job', 'Begin', 1 );\n	my $jobIdForFilename = 1;\n	\
$jobIdForFilename = 0 if ( defined( $params{'outfi\
le'} ) );\n	my (@filename_list) = ();\n\n	# Query \
sequence\n	if ( defined( $ARGV[0] ) ) {    # Bare \
option\n		if ( -f $ARGV[0] || $ARGV[0] eq '-' ) { \
   # File\n			push( @filename_list, $ARGV[0] );\n	\
	}\n	}\n	if ( $params{'sequence'} ) {             \
      # Via --sequence\n		if ( -f $params{'sequenc\
e'} || $params{'sequence'} eq '-' ) {    # File\n	\
		push( @filename_list, $params{'sequence'} );\n		\
}\n	}\n\n	$/ = '>';\n	foreach my $filename (@filen\
ame_list) {\n		open( my $INFILE, '<', $filename )\\
n		  or die \"Error: unable to open file $filename\
 ($!)\";\n		while (<$INFILE>) {\n			my $seq = $_;\\
n			$seq =~ s/>$//;\n			if ( $seq =~ m/(\\S+)/ ) {\
\n				print STDERR \"Submitting job for: $1\\n\"\n\
				  if ( $outputLevel > 0 );\n				$seq = '>' . $\
seq;\n				&print_debug_message( 'multi_submit_job'\
, $seq, 11 );\n				&submit_job($seq);\n				$params\
{'outfile'} = undef if ( $jobIdForFilename == 1 );\
\n			}\n		}\n		close $INFILE;\n	}\n	print_debug_me\
ssage( 'multi_submit_job', 'End', 1 );\n}\n\n=head\
2 list_file_submit_job()\n\nSubmit multiple jobs u\
sing a file containing a list of entry identifiers\
 as \ninput.\n\n  &list_file_submit_job($list_file\
name)\n\n=cut\n\nsub list_file_submit_job {\n	my $\
filename         = shift;\n	my $jobIdForFilename =\
 1;\n	$jobIdForFilename = 0 if ( defined( $params{\
'outfile'} ) );\n\n	# Iterate over identifiers, su\
bmitting each job\n	open( my $LISTFILE, '<', $file\
name )\n	  or die 'Error: unable to open file ' . \
$filename . ' (' . $! . ')';\n	while (<$LISTFILE>)\
 {\n		my $line = $_;\n		chomp($line);\n		if ( $lin\
e ne '' ) {\n			&print_debug_message( 'list_file_s\
ubmit_job', 'line: ' . $line, 2 );\n			if ( $line \
=~ m/\\w:\\w/ ) {    # Check this is an identifier\
\n				print STDERR \"Submitting job for: $line\\n\\
"\n				  if ( $outputLevel > 0 );\n				&submit_job\
($line);\n			}\n			else {\n				print STDERR\n\"War\
ning: line \\\"$line\\\" is not recognised as an i\
dentifier\\n\";\n			}\n		}\n		$params{'outfile'} =\
 undef if ( $jobIdForFilename == 1 );\n	}\n	close \
$LISTFILE;\n}\n\n=head2 load_data()\n\nLoad sequen\
ce data from file or option specified on the comma\
nd-line.\n\n  &load_data();\n\n=cut\n\nsub load_da\
ta {\n	print_debug_message( 'load_data', 'Begin', \
1 );\n	my $retSeq;\n\n	# Query sequence\n	if ( def\
ined( $ARGV[0] ) ) {    # Bare option\n		if ( -f $\
ARGV[0] || $ARGV[0] eq '-' ) {    # File\n			$retS\
eq = &read_file( $ARGV[0] );\n		}\n		else {       \
                              # DB:ID or sequence\\
n			$retSeq = $ARGV[0];\n		}\n	}\n	if ( $params{'s\
equence'} ) {                   # Via --sequence\n\
		if ( -f $params{'sequence'} || $params{'sequence\
'} eq '-' ) {    # File\n			$retSeq = &read_file( \
$params{'sequence'} );\n		}\n		else {    # DB:ID o\
r sequence\n			$retSeq = $params{'sequence'};\n		}\
\n	}\n	print_debug_message( 'load_data', 'End', 1 \
);\n	return $retSeq;\n}\n\n=head2 load_params()\n\\
nLoad job parameters from command-line options.\n\\
n  &load_params();\n\n=cut\n\nsub load_params {\n	\
print_debug_message( 'load_params', 'Begin', 1 );\\
n\n	# Database(s) to search\n	my (@dbList) = split\
 /[ ,]/, $params{'database'};\n	$tool_params{'data\
base'} = \\@dbList;\n\n	# Match/missmatch\n	if ( $\
params{'match'} && $params{'missmatch'} ) {\n		$to\
ol_params{'match_scores'} =\n		  $params{'match'} \
. ',' . $params{'missmatch'};\n	}\n	\n	# Compatabi\
lity options, old command-line\n	if(!$tool_params{\
'alignments'} && $params{'numal'}) {\n		$tool_para\
ms{'alignments'} = $params{'numal'};\n	}\n	if(!$to\
ol_params{'gapopen'} && $params{'opengap'}) {\n		$\
tool_params{'gapopen'} = $params{'opengap'};\n	}\n\
	if(!$tool_params{'gapext'} && $params{'extendgap'\
}) {\n		$tool_params{'gapext'} = $params{'extendga\
p'};\n	}\n\n	print_debug_message( 'load_params', '\
End', 1 );\n}\n\n=head2 client_poll()\n\nClient-si\
de job polling.\n\n  &client_poll($job_id);\n\n=cu\
t\n\nsub client_poll {\n	print_debug_message( 'cli\
ent_poll', 'Begin', 1 );\n	my $jobid  = shift;\n	m\
y $status = 'PENDING';\n\n	my $errorCount = 0;\n	w\
hile ($status eq 'RUNNING'\n		|| $status eq 'PENDI\
NG'\n		|| ( $status eq 'ERROR' && $errorCount < 2 \
) )\n	{\n		$status = rest_get_status($jobid);\n		p\
rint STDERR \"$status\\n\" if ( $outputLevel > 0 )\
;\n		if ( $status eq 'ERROR' ) {\n			$errorCount++\
;\n		}\n		elsif ( $errorCount > 0 ) {\n			$errorCo\
unt--;\n		}\n		if (   $status eq 'RUNNING'\n			|| \
$status eq 'PENDING'\n			|| $status eq 'ERROR' )\n\
		{\n\n			# Wait before polling again.\n			sleep $\
checkInterval;\n		}\n	}\n	print_debug_message( 'cl\
ient_poll', 'End', 1 );\n	return $status;\n}\n\n=h\
ead2 get_results()\n\nGet the results for a job id\
entifier.\n\n  &get_results($job_id);\n\n=cut\n\ns\
ub get_results {\n	print_debug_message( 'get_resul\
ts', 'Begin', 1 );\n	my $jobid = shift;\n	print_de\
bug_message( 'get_results', 'jobid: ' . $jobid, 1 \
);\n\n	# Verbose\n	if ( $outputLevel > 1 ) {\n		pr\
int 'Getting results for job ', $jobid, \"\\n\";\n\
	}\n\n	# Check status, and wait if not finished\n	\
client_poll($jobid);\n\n	# Use JobId if output fil\
e name is not defined\n	unless ( defined( $params{\
'outfile'} ) ) {\n		$params{'outfile'} = $jobid;\n\
	}\n\n	# Get list of data types\n	my (@resultTypes\
) = rest_get_result_types($jobid);\n\n	# Get the d\
ata and write it to a file\n	if ( defined( $params\
{'outformat'} ) ) {    # Specified data type\n		my\
 $selResultType;\n		foreach my $resultType (@resul\
tTypes) {\n			if ( $resultType->{'identifier'} eq \
$params{'outformat'} ) {\n				$selResultType = $re\
sultType;\n			}\n		}\n		if ( defined($selResultTyp\
e) ) {\n			my $result =\n			  rest_get_result( $jo\
bid, $selResultType->{'identifier'} );\n			if ( $p\
arams{'outfile'} eq '-' ) {\n				write_file( $para\
ms{'outfile'}, $result );\n			}\n			else {\n				wr\
ite_file(\n					$params{'outfile'} . '.'\n					  .\
 $selResultType->{'identifier'} . '.'\n					  . $s\
elResultType->{'fileSuffix'},\n					$result\n				)\
;\n			}\n		}\n		else {\n			die 'Error: unknown res\
ult format \"' . $params{'outformat'} . '\"';\n		}\
\n	}\n	else {    # Data types available\n		      #\
 Write a file for each output type\n		for my $resu\
ltType (@resultTypes) {\n			if ( $outputLevel > 1 \
) {\n				print STDERR 'Getting ', $resultType->{'i\
dentifier'}, \"\\n\";\n			}\n			my $result = rest_\
get_result( $jobid, $resultType->{'identifier'} );\
\n			if ( $params{'outfile'} eq '-' ) {\n				write\
_file( $params{'outfile'}, $result );\n			}\n			el\
se {\n				write_file(\n					$params{'outfile'} . '\
.'\n					  . $resultType->{'identifier'} . '.'\n		\
			  . $resultType->{'fileSuffix'},\n					$result\\
n				);\n			}\n		}\n	}\n	print_debug_message( 'get\
_results', 'End', 1 );\n}\n\n=head2 read_file()\n\\
nRead a file into a scalar. The special filename '\
-' can be used to read from \nstandard input (STDI\
N).\n\n  my $data = &read_file($filename);\n\n=cut\
\n\nsub read_file {\n	print_debug_message( 'read_f\
ile', 'Begin', 1 );\n	my $filename = shift;\n	prin\
t_debug_message( 'read_file', 'filename: ' . $file\
name, 2 );\n	my ( $content, $buffer );\n	if ( $fil\
ename eq '-' ) {\n		while ( sysread( STDIN, $buffe\
r, 1024 ) ) {\n			$content .= $buffer;\n		}\n	}\n	\
else {    # File\n		open( my $FILE, '<', $filename\
 )\n		  or die \"Error: unable to open input file \
$filename ($!)\";\n		while ( sysread( $FILE, $buff\
er, 1024 ) ) {\n			$content .= $buffer;\n		}\n		cl\
ose($FILE);\n	}\n	print_debug_message( 'read_file'\
, 'End', 1 );\n	return $content;\n}\n\n=head2 writ\
e_file()\n\nWrite data to a file. The special file\
name '-' can be used to write to \nstandard output\
 (STDOUT).\n\n  &write_file($filename, $data);\n\n\
=cut\n\nsub write_file {\n	print_debug_message( 'w\
rite_file', 'Begin', 1 );\n	my ( $filename, $data \
) = @_;\n	print_debug_message( 'write_file', 'file\
name: ' . $filename, 2 );\n	if ( $outputLevel > 0 \
) {\n		print STDERR 'Creating result file: ' . $fi\
lename . \"\\n\";\n	}\n	if ( $filename eq '-' ) {\\
n		print STDOUT $data;\n	}\n	else {\n		open( my $F\
ILE, '>', $filename )\n		  or die \"Error: unable \
to open output file $filename ($!)\";\n		syswrite(\
 $FILE, $data );\n		close($FILE);\n	}\n	print_debu\
g_message( 'write_file', 'End', 1 );\n}\n\n=head2 \
usage()\n\nPrint program usage message.\n\n  &usag\
e();\n\n=cut\n\nsub usage {\n	print STDERR <<EOF\n\
NCBI BLAST\n==========\n   \nRapid sequence databa\
se search programs utilizing the BLAST algorithm\n\
    \n[Required]\n\n  -p, --program      : str  : \
BLAST program to use, see --paramDetail program\n \
 -D, --database     : str  : database(s) to search\
, space separated. See\n                          \
    --paramDetail database\n      --stype        :\
 str  : query sequence type, see --paramDetail sty\
pe\n  seqFile            : file : query sequence (\
\"-\" for STDIN, \\@filename for\n                \
              identifier list file)\n\n[Optional]\\
n\n  -m, --matrix       : str  : scoring matrix, s\
ee --paramDetail matrix\n  -e, --exp          : re\
al : 0<E<= 1000. Statistical significance threshol\
d \n                              for reporting da\
tabase sequence matches.\n  -f, --filter       :  \
    : filter the query sequence for low complexity\
 \n                              regions, see --pa\
ramDetail filter\n  -A, --align        : int  : pa\
irwise alignment format, see --paramDetail align\n\
  -s, --scores       : int  : number of scores to \
be reported\n  -n, --alignments   : int  : number \
of alignments to report\n  -u, --match        : in\
t  : Match score (BLASTN only)\n  -v, --mismatch  \
   : int  : Mismatch score (BLASTN only)\n  -o, --\
gapopen      : int  : Gap open penalty\n  -x, --ga\
pext       : int  : Gap extension penalty\n  -d, -\
-dropoff      : int  : Drop-off\n  -g, --gapalign \
    :      : Optimise gapped alignments\n      --s\
eqrange     : str  : region within input to use as\
 query\n      --multifasta   :      : treat input \
as a set of fasta formatted sequences\n\n[General]\
\n\n  -h, --help        :      : prints this help \
text\n      --async       :      : forces to make \
an asynchronous query\n      --email       : str  \
: e-mail address\n      --title       : str  : tit\
le for job\n      --status      :      : get job s\
tatus\n      --resultTypes :      : get available \
result types for job\n      --polljob     :      :\
 poll for the status of a job\n      --jobid      \
 : str  : jobid that was returned when an asynchro\
nous job \n                             was submit\
ted.\n      --outfile     : str  : file name for r\
esults (default is jobid;\n                       \
      \"-\" for STDOUT)\n      --outformat   : str\
  : result format to retrieve\n      --params     \
 :      : list input parameters\n      --paramDeta\
il : str  : display details for input parameter\n \
     --quiet       :      : decrease output\n     \
 --verbose     :      : increase output\n      --t\
race       :      : show SOAP messages being inter\
changed \n   \nSynchronous job:\n\n  The results/e\
rrors are returned as soon as the job is finished.\
\n  Usage: $scriptName --email <your\\@email> [opt\
ions...] seqFile\n  Returns: results as an attachm\
ent\n\nAsynchronous job:\n\n  Use this if you want\
 to retrieve the results at a later time. The resu\
lts \n  are stored for up to 24 hours. 	\n  Usage:\
 $scriptName --async --email <your\\@email> [optio\
ns...] seqFile\n  Returns: jobid\n\n  Use the jobi\
d to query for the status of the job. If the job i\
s finished, \n  it also returns the results/errors\
.\n  Usage: $scriptName --polljob --jobid <jobId> \
[--outfile string]\n  Returns: string indicating t\
he status of the job and if applicable, results \n\
  as an attachment.\n\nFurther information:\n\n  h\
ttp://www.ebi.ac.uk/Tools/webservices/services/sss\
/ncbi_blast_rest\n  http://www.ebi.ac.uk/Tools/web\
services/tutorials/perl\n\nSupport/Feedback:\n\n  \
http://www.ebi.ac.uk/support/\nEOF\n}\n\n=head1 FE\
EDBACK/SUPPORT\n\nPlease contact us at L<http://ww\
w.ebi.ac.uk/support/> if you have any \nfeedback, \
suggestions or issues with the service or this cli\
ent.\n\n=cut\n","\n=head1 NAME\n\nwublast_lwp.pl\n\
\n=head1 DESCRIPTION\n\nWU-BLAST REST web service \
Perl client using L<LWP>.\n\nTested with:\n\n=over\
\n\n=item *\nL<LWP> 5.79, L<XML::Simple> 2.12 and \
Perl 5.8.3\n\n=item *\nL<LWP> 5.805, L<XML::Simple\
> 2.14 and Perl 5.8.7\n\n=item *\nL<LWP> 5.820, L<\
XML::Simple> 2.18 and Perl 5.10.0 (Ubuntu 9.04)\n\\
n=back\n\nFor further information see:\n\n=over\n\\
n=item *\nL<http://www.ebi.ac.uk/Tools/webservices\
/services/sss/wu_blast_rest>\n\n=item *\nL<http://\
www.ebi.ac.uk/Tools/webservices/tutorials/perl>\n\\
n=back\n\n=head1 VERSION\n\n$Id: wublast_lwp.pl 13\
17 2009-09-03 15:44:11Z hpm $\n\n=cut\n\nuse stric\
t;\nuse warnings;\n\nuse English;\nuse LWP;\nuse X\
ML::Simple;\nuse Getopt::Long qw(:config no_ignore\
_case bundling);\nuse File::Basename;\nuse Data::D\
umper;\n\nmy $baseUrl = 'http://www.ebi.ac.uk/Tool\
s/services/rest/wublast';\n\nmy $checkInterval = 3\
;\n\nmy $outputLevel = 1;\n\nmy $numOpts = scalar(\
@ARGV);\nmy %params = ( 'debugLevel' => 0 );\n\nmy\
 %tool_params = ();\nGetOptions(\n\n	# Tool specif\
ic options\n	'program|p=s'     => \\$tool_params{'\
program'},      # BLAST program\n	'database|D=s'  \
  => \\$params{'database'},     # Search database\\
n	'matrix|m=s'      => \\$tool_params{'matrix'},  \
     # Scoring matrix\n	'exp|E=f'         => \\$to\
ol_params{'exp'},          # E-value threshold\n	'\
viewfilter|e'    => \\$tool_params{'viewfilter'}, \
  # Display filtered sequence\n	'filter|f=s'      \
=> \\$tool_params{'filter'},       # Low complexit\
y filter name\n	'alignments|n=i'  => \\$tool_param\
s{'alignments'},   # Number of alignments\n	'score\
s|s=i'      => \\$tool_params{'scores'},       # N\
umber of scores\n	'sensitivity|S=s' => \\$tool_par\
ams{'sensitivity'},  # Search sensitivity\n	'sort|\
t=s'        => \\$tool_params{'sort'},         # S\
ort hits by...\n	'stats|T=s'       => \\$tool_para\
ms{'stats'},        # Scoring statistic to use\n	'\
strand|d=s'      => \\$tool_params{'strand'},     \
  # Strand to use\n	'topcombon|c=i'   => \\$tool_p\
arams{'topcombon'},    # Consistent sets of HSPs\n\
	'align|A=i'       => \\$tool_params{'align'},   #\
 Pairwise alignment format\n	'stype=s' => \\$tool_\
params{'stype'},    # Sequence type 'protein' or '\
dna'\n	'sequence=s' => \\$params{'sequence'},     \
    # Query sequence file or DB:ID\n	'multifasta' \
=> \\$params{'multifasta'},       # Multiple fasta\
 input\n\n	# Compatability options, old command-li\
ne.\n	'echofilter|e'    => \\$params{'echofilter'}\
,   # Display filtered sequence\n	'b=i'  => \\$par\
ams{'numal'},        # Number of alignments\n	'app\
xml=s'        => \\$params{'appxml'},       # Appl\
ication XML\n\n	# Generic options\n	'email=s'     \
  => \\$params{'email'},          # User e-mail ad\
dress\n	'title=s'       => \\$params{'title'},    \
      # Job title\n	'outfile=s'     => \\$params{'\
outfile'},        # Output file name\n	'outformat=\
s'   => \\$params{'outformat'},      # Output file\
 type\n	'jobid=s'       => \\$params{'jobid'},    \
      # JobId\n	'help|h'        => \\$params{'help\
'},           # Usage help\n	'async'         => \\\
$params{'async'},          # Asynchronous submissi\
on\n	'polljob'       => \\$params{'polljob'},     \
   # Get results\n	'resultTypes'   => \\$params{'r\
esultTypes'},    # Get result types\n	'status'    \
    => \\$params{'status'},         # Get status\n\
	'params'        => \\$params{'params'},         #\
 List input parameters\n	'paramDetail=s' => \\$par\
ams{'paramDetail'},    # Get details for parameter\
\n	'quiet'         => \\$params{'quiet'},         \
 # Decrease output level\n	'verbose'       => \\$p\
arams{'verbose'},        # Increase output level\n\
	'debugLevel=i'  => \\$params{'debugLevel'},     #\
 Debug output level\n	'baseUrl=s'     => \\$baseUr\
l,                  # Base URL for service.\n);\ni\
f ( $params{'verbose'} ) { $outputLevel++ }\nif ( \
$params{'$quiet'} )  { $outputLevel-- }\n\n&print_\
debug_message( 'MAIN', 'LWP::VERSION: ' . $LWP::VE\
RSION,\n	1 );\n\n&print_debug_message( 'MAIN', \"p\
arams:\\n\" . Dumper( \\%params ),           11 );\
\n&print_debug_message( 'MAIN', \"tool_params:\\n\\
" . Dumper( \\%tool_params ), 11 );\n\nmy $scriptN\
ame = basename( $0, () );\n\nif ( $params{'help'} \
|| $numOpts == 0 ) {\n	&usage();\n	exit(0);\n}\n\n\
&print_debug_message( 'MAIN', 'baseUrl: ' . $baseU\
rl, 1 );\n\nif (\n	!(\n		   $params{'polljob'}\n		\
|| $params{'resultTypes'}\n		|| $params{'status'}\\
n		|| $params{'params'}\n		|| $params{'paramDetail\
'}\n	)\n	&& !( defined( $ARGV[0] ) || defined( $pa\
rams{'sequence'} ) )\n  )\n{\n\n	# Bad argument co\
mbination, so print error message and usage\n	prin\
t STDERR 'Error: bad option combination', \"\\n\";\
\n	&usage();\n	exit(1);\n}\n\nelsif ( $params{'par\
ams'} ) {\n	&print_tool_params();\n}\n\nelsif ( $p\
arams{'paramDetail'} ) {\n	&print_param_details( $\
params{'paramDetail'} );\n}\n\nelsif ( $params{'st\
atus'} && defined( $params{'jobid'} ) ) {\n	&print\
_job_status( $params{'jobid'} );\n}\n\nelsif ( $pa\
rams{'resultTypes'} && defined( $params{'jobid'} )\
 ) {\n	&print_result_types( $params{'jobid'} );\n}\
\n\nelsif ( $params{'polljob'} && defined( $params\
{'jobid'} ) ) {\n	&get_results( $params{'jobid'} )\
;\n}\n\nelse {\n\n	# Multiple input sequence mode,\
 assume fasta format.\n	if ( $params{'multifasta'}\
 ) {\n		&multi_submit_job();\n	}\n\n	# Entry ident\
ifier list file.\n	elsif (( defined( $params{'sequ\
ence'} ) && $params{'sequence'} =~ m/^\\@/ )\n		||\
 ( defined( $ARGV[0] ) && $ARGV[0] =~ m/^\\@/ ) )\\
n	{\n		my $list_filename = $params{'sequence'} || \
$ARGV[0];\n		$list_filename =~ s/^\\@//;\n		&list_\
file_submit_job($list_filename);\n	}\n\n	# Default\
: single sequence/identifier.\n	else {\n\n		# Load\
 the sequence data and submit.\n		&submit_job( &lo\
ad_data() );\n	}\n}\n\n=head1 FUNCTIONS\n\n=cut\n\\
n\n=head2 rest_request()\n\nPerform a REST request\
.\n\n  my $response_str = &rest_request($url);\n\n\
=cut\n\nsub rest_request {\n	print_debug_message( \
'rest_request', 'Begin', 11 );\n	my $requestUrl = \
shift;\n	print_debug_message( 'rest_request', 'URL\
: ' . $requestUrl, 11 );\n\n	# Create a user agent\
\n	my $ua = LWP::UserAgent->new();\n	'$Revision: 1\
317 $' =~ m/(\\d+)/;\n	$ua->agent(\"EBI-Sample-Cli\
ent/$1 ($scriptName; $OSNAME) \" . $ua->agent());\\
n	$ua->env_proxy;\n\n	# Perform the request\n	my $\
response = $ua->get($requestUrl);\n	print_debug_me\
ssage( 'rest_request', 'HTTP status: ' . $response\
->code,\n		11 );\n\n	# Check for HTTP error codes\\
n	if ( $response->is_error ) {\n		$response->conte\
nt() =~ m/<h1>([^<]+)<\\/h1>/;\n		die 'http status\
: ' . $response->code . ' ' . $response->message .\
 '  ' . $1;\n	}\n	print_debug_message( 'rest_reque\
st', 'End', 11 );\n\n	# Return the response data\n\
	return $response->content();\n}\n\n=head2 rest_ge\
t_parameters()\n\nGet list of tool parameter names\
.\n\n  my (@param_list) = &rest_get_parameters();\\
n\n=cut\n\nsub rest_get_parameters {\n	print_debug\
_message( 'rest_get_parameters', 'Begin', 1 );\n	m\
y $url                = $baseUrl . '/parameters/';\
\n	my $param_list_xml_str = rest_request($url);\n	\
my $param_list_xml     = XMLin($param_list_xml_str\
);\n	my (@param_list)       = @{ $param_list_xml->\
{'id'} };\n	print_debug_message( 'rest_get_paramet\
ers', 'End', 1 );\n	return (@param_list);\n}\n\n=h\
ead2 rest_get_parameter_details()\n\nGet details o\
f a tool parameter.\n\n  my $paramDetail = &rest_g\
et_parameter_details($param_name);\n\n=cut\n\nsub \
rest_get_parameter_details {\n	print_debug_message\
( 'rest_get_parameter_details', 'Begin', 1 );\n	my\
 $parameterId = shift;\n	print_debug_message( 'res\
t_get_parameter_details',\n		'parameterId: ' . $pa\
rameterId, 1 );\n	my $url                  = $base\
Url . '/parameterdetails/' . $parameterId;\n	my $p\
aram_detail_xml_str = rest_request($url);\n	my $pa\
ram_detail_xml     = XMLin($param_detail_xml_str);\
\n	print_debug_message( 'rest_get_parameter_detail\
s', 'End', 1 );\n	return ($param_detail_xml);\n}\n\
\n=head2 rest_run()\n\nSubmit a job.\n\n  my $job_\
id = &rest_run($email, $title, \\%params );\n\n=cu\
t\n\nsub rest_run {\n	print_debug_message( 'rest_r\
un', 'Begin', 1 );\n	my $email  = shift;\n	my $tit\
le  = shift;\n	my $params = shift;\n	print_debug_m\
essage( 'rest_run', 'email: ' . $email, 1 );\n	if \
( defined($title) ) {\n		print_debug_message( 'res\
t_run', 'title: ' . $title, 1 );\n	}\n	print_debug\
_message( 'rest_run', 'params: ' . Dumper($params)\
, 1 );\n\n	# User agent to perform http requests\n\
	my $ua = LWP::UserAgent->new();\n	$ua->env_proxy;\
\n\n	# Clean up parameters\n	my (%tmp_params) = %{\
$params};\n	$tmp_params{'email'} = $email;\n	$tmp_\
params{'title'} = $title;\n	foreach my $param_name\
 ( keys(%tmp_params) ) {\n		if ( !defined( $tmp_pa\
rams{$param_name} ) ) {\n			delete $tmp_params{$pa\
ram_name};\n		}\n	}\n\n	# Submit the job as a POST\
\n	my $url = $baseUrl . '/run';\n	my $response = $\
ua->post( $url, \\%tmp_params );\n	print_debug_mes\
sage( 'rest_run', 'HTTP status: ' . $response->cod\
e, 11 );\n	print_debug_message( 'rest_run',\n		're\
quest: ' . $response->request()->content(), 11 );\\
n\n	# Check for HTTP error codes\n	if ( $response-\
>is_error ) {\n		$response->content() =~ m/<h1>([^\
<]+)<\\/h1>/;\n		die 'http status: ' . $response->\
code . ' ' . $response->message . '  ' . $1;\n	}\n\
\n	# The job id is returned\n	my $job_id = $respon\
se->content();\n	print_debug_message( 'rest_run', \
'End', 1 );\n	return $job_id;\n}\n\n=head2 rest_ge\
t_status()\n\nCheck the status of a job.\n\n  my $\
status = &rest_get_status($job_id);\n\n=cut\n\nsub\
 rest_get_status {\n	print_debug_message( 'rest_ge\
t_status', 'Begin', 1 );\n	my $job_id = shift;\n	p\
rint_debug_message( 'rest_get_status', 'jobid: ' .\
 $job_id, 2 );\n	my $status_str = 'UNKNOWN';\n	my \
$url        = $baseUrl . '/status/' . $job_id;\n	$\
status_str = &rest_request($url);\n	print_debug_me\
ssage( 'rest_get_status', 'status_str: ' . $status\
_str, 2 );\n	print_debug_message( 'rest_get_status\
', 'End', 1 );\n	return $status_str;\n}\n\n=head2 \
rest_get_result_types()\n\nGet list of result type\
s for finished job.\n\n  my (@result_types) = &res\
t_get_result_types($job_id);\n\n=cut\n\nsub rest_g\
et_result_types {\n	print_debug_message( 'rest_get\
_result_types', 'Begin', 1 );\n	my $job_id = shift\
;\n	print_debug_message( 'rest_get_result_types', \
'jobid: ' . $job_id, 2 );\n	my (@resultTypes);\n	m\
y $url                      = $baseUrl . '/resultt\
ypes/' . $job_id;\n	my $result_type_list_xml_str =\
 &rest_request($url);\n	my $result_type_list_xml  \
   = XMLin($result_type_list_xml_str);\n	(@resultT\
ypes) = @{ $result_type_list_xml->{'type'} };\n	pr\
int_debug_message( 'rest_get_result_types',\n		sca\
lar(@resultTypes) . ' result types', 2 );\n	print_\
debug_message( 'rest_get_result_types', 'End', 1 )\
;\n	return (@resultTypes);\n}\n\n=head2 rest_get_r\
esult()\n\nGet result data of a specified type for\
 a finished job.\n\n  my $result = rest_get_result\
($job_id, $result_type);\n\n=cut\n\nsub rest_get_r\
esult {\n	print_debug_message( 'rest_get_result', \
'Begin', 1 );\n	my $job_id = shift;\n	my $type   =\
 shift;\n	print_debug_message( 'rest_get_result', \
'jobid: ' . $job_id, 1 );\n	print_debug_message( '\
rest_get_result', 'type: ' . $type,    1 );\n	my $\
url    = $baseUrl . '/result/' . $job_id . '/' . $\
type;\n	my $result = &rest_request($url);\n	print_\
debug_message( 'rest_get_result', length($result) \
. ' characters',\n		1 );\n	print_debug_message( 'r\
est_get_result', 'End', 1 );\n	return $result;\n}\\
n\n\n=head2 print_debug_message()\n\nPrint debug m\
essage at specified debug level.\n\n  &print_debug\
_message($method_name, $message, $level);\n\n=cut\\
n\nsub print_debug_message {\n	my $function_name =\
 shift;\n	my $message       = shift;\n	my $level  \
       = shift;\n	if ( $level <= $params{'debugLev\
el'} ) {\n		print STDERR '[', $function_name, '()]\
 ', $message, \"\\n\";\n	}\n}\n\n=head2 print_tool\
_params()\n\nPrint list of tool parameters.\n\n  &\
print_tool_params();\n\n=cut\n\nsub print_tool_par\
ams {\n	print_debug_message( 'print_tool_params', \
'Begin', 1 );\n	my (@param_list) = &rest_get_param\
eters();\n	foreach my $param ( sort(@param_list) )\
 {\n		print $param, \"\\n\";\n	}\n	print_debug_mes\
sage( 'print_tool_params', 'End', 1 );\n}\n\n=head\
2 print_param_details()\n\nPrint details of a tool\
 parameter.\n\n  &print_param_details($param_name)\
;\n\n=cut\n\nsub print_param_details {\n	print_deb\
ug_message( 'print_param_details', 'Begin', 1 );\n\
	my $paramName = shift;\n	print_debug_message( 'pr\
int_param_details', 'paramName: ' . $paramName, 2 \
);\n	my $paramDetail = &rest_get_parameter_details\
($paramName);\n	print $paramDetail->{'name'}, \"\\\
t\", $paramDetail->{'type'}, \"\\n\";\n	print $par\
amDetail->{'description'}, \"\\n\";\n	foreach my $\
value ( @{ $paramDetail->{'values'}->{'value'} } )\
 {\n		print $value->{'value'};\n		if ( $value->{'d\
efaultValue'} eq 'true' ) {\n			print \"\\t\", 'de\
fault';\n		}\n		print \"\\n\";\n		print \"\\t\", $\
value->{'label'}, \"\\n\";\n	}\n	print_debug_messa\
ge( 'print_param_details', 'End', 1 );\n}\n\n=head\
2 print_job_status()\n\nPrint status of a job.\n\n\
  &print_job_status($job_id);\n\n=cut\n\nsub print\
_job_status {\n	print_debug_message( 'print_job_st\
atus', 'Begin', 1 );\n	my $jobid = shift;\n	print_\
debug_message( 'print_job_status', 'jobid: ' . $jo\
bid, 1 );\n	if ( $outputLevel > 0 ) {\n		print STD\
ERR 'Getting status for job ', $jobid, \"\\n\";\n	\
}\n	my $result = &rest_get_status($jobid);\n	print\
 \"$result\\n\";\n	if ( $result eq 'FINISHED' && $\
outputLevel > 0 ) {\n		print STDERR \"To get resul\
ts: $scriptName --polljob --jobid \" . $jobid\n		 \
 . \"\\n\";\n	}\n	print_debug_message( 'print_job_\
status', 'End', 1 );\n}\n\n=head2 print_result_typ\
es()\n\nPrint available result types for a job.\n\\
n  &print_result_types($job_id);\n\n=cut\n\nsub pr\
int_result_types {\n	print_debug_message( 'result_\
types', 'Begin', 1 );\n	my $jobid = shift;\n	print\
_debug_message( 'result_types', 'jobid: ' . $jobid\
, 1 );\n	if ( $outputLevel > 0 ) {\n		print STDERR\
 'Getting result types for job ', $jobid, \"\\n\";\
\n	}\n	my $status = &rest_get_status($jobid);\n	if\
 ( $status eq 'PENDING' || $status eq 'RUNNING' ) \
{\n		print STDERR 'Error: Job status is ', $status\
,\n		  '. To get result types the job must be fini\
shed.', \"\\n\";\n	}\n	else {\n		my (@resultTypes)\
 = &rest_get_result_types($jobid);\n		if ( $output\
Level > 0 ) {\n			print STDOUT 'Available result t\
ypes:', \"\\n\";\n		}\n		foreach my $resultType (@\
resultTypes) {\n			print STDOUT $resultType->{'ide\
ntifier'}, \"\\n\";\n			if ( defined( $resultType-\
>{'label'} ) ) {\n				print STDOUT \"\\t\", $resul\
tType->{'label'}, \"\\n\";\n			}\n			if ( defined(\
 $resultType->{'description'} ) ) {\n				print STD\
OUT \"\\t\", $resultType->{'description'}, \"\\n\"\
;\n			}\n			if ( defined( $resultType->{'mediaType\
'} ) ) {\n				print STDOUT \"\\t\", $resultType->{\
'mediaType'}, \"\\n\";\n			}\n			if ( defined( $re\
sultType->{'fileSuffix'} ) ) {\n				print STDOUT \\
"\\t\", $resultType->{'fileSuffix'}, \"\\n\";\n			\
}\n		}\n		if ( $status eq 'FINISHED' && $outputLev\
el > 0 ) {\n			print STDERR \"\\n\", 'To get resul\
ts:', \"\\n\",\n			  \"  $scriptName --polljob --j\
obid \" . $params{'jobid'} . \"\\n\",\n			  \"  $s\
criptName --polljob --outformat <type> --jobid \"\\
n			  . $params{'jobid'} . \"\\n\";\n		}\n	}\n	pri\
nt_debug_message( 'result_types', 'End', 1 );\n}\n\
\n=head2 submit_job()\n\nSubmit a job to the servi\
ce.\n\n  &submit_job($seq);\n\n=cut\n\nsub submit_\
job {\n	print_debug_message( 'submit_job', 'Begin'\
, 1 );\n\n	# Set input sequence\n	$tool_params{'se\
quence'} = shift;\n\n	# Load parameters\n	&load_pa\
rams();\n\n	# Submit the job\n	my $jobid = &rest_r\
un( $params{'email'}, $params{'title'}, \\%tool_pa\
rams );\n\n	# Simulate sync/async mode\n	if ( defi\
ned( $params{'async'} ) ) {\n		print STDOUT $jobid\
, \"\\n\";\n		if ( $outputLevel > 0 ) {\n			print \
STDERR\n			  \"To check status: $scriptName --stat\
us --jobid $jobid\\n\";\n		}\n	}\n	else {\n		if ( \
$outputLevel > 0 ) {\n			print STDERR \"JobId: $jo\
bid\\n\";\n		}\n		sleep 1;\n		&get_results($jobid)\
;\n	}\n	print_debug_message( 'submit_job', 'End', \
1 );\n}\n\n=head2 multi_submit_job()\n\nSubmit mul\
tiple jobs assuming input is a collection of fasta\
 formatted sequences.\n\n  &multi_submit_job();\n\\
n=cut\n\nsub multi_submit_job {\n	print_debug_mess\
age( 'multi_submit_job', 'Begin', 1 );\n	my $jobId\
ForFilename = 1;\n	$jobIdForFilename = 0 if ( defi\
ned( $params{'outfile'} ) );\n	my (@filename_list)\
 = ();\n\n	# Query sequence\n	if ( defined( $ARGV[\
0] ) ) {    # Bare option\n		if ( -f $ARGV[0] || $\
ARGV[0] eq '-' ) {    # File\n			push( @filename_l\
ist, $ARGV[0] );\n		}\n	}\n	if ( $params{'sequence\
'} ) {                   # Via --sequence\n		if ( \
-f $params{'sequence'} || $params{'sequence'} eq '\
-' ) {    # File\n			push( @filename_list, $params\
{'sequence'} );\n		}\n	}\n\n	$/ = '>';\n	foreach m\
y $filename (@filename_list) {\n		open( my $INFILE\
, '<', $filename )\n		  or die \"Error: unable to \
open file $filename ($!)\";\n		while (<$INFILE>) {\
\n			my $seq = $_;\n			$seq =~ s/>$//;\n			if ( $s\
eq =~ m/(\\S+)/ ) {\n				print STDERR \"Submitting\
 job for: $1\\n\"\n				  if ( $outputLevel > 0 );\\
n				$seq = '>' . $seq;\n				&print_debug_message(\
 'multi_submit_job', $seq, 11 );\n				&submit_job(\
$seq);\n				$params{'outfile'} = undef if ( $jobId\
ForFilename == 1 );\n			}\n		}\n		close $INFILE;\n\
	}\n	print_debug_message( 'multi_submit_job', 'End\
', 1 );\n}\n\n=head2 list_file_submit_job()\n\nSub\
mit multiple jobs using a file containing a list o\
f entry identifiers as \ninput.\n\n  &list_file_su\
bmit_job($list_filename)\n\n=cut\n\nsub list_file_\
submit_job {\n	my $filename         = shift;\n	my \
$jobIdForFilename = 1;\n	$jobIdForFilename = 0 if \
( defined( $params{'outfile'} ) );\n\n	# Iterate o\
ver identifiers, submitting each job\n	open( my $L\
ISTFILE, '<', $filename )\n	  or die 'Error: unabl\
e to open file ' . $filename . ' (' . $! . ')';\n	\
while (<$LISTFILE>) {\n		my $line = $_;\n		chomp($\
line);\n		if ( $line ne '' ) {\n			&print_debug_me\
ssage( 'list_file_submit_job', 'line: ' . $line, 2\
 );\n			if ( $line =~ m/\\w:\\w/ ) {    # Check th\
is is an identifier\n				print STDERR \"Submitting\
 job for: $line\\n\"\n				  if ( $outputLevel > 0 \
);\n				&submit_job($line);\n			}\n			else {\n				\
print STDERR\n\"Warning: line \\\"$line\\\" is not\
 recognised as an identifier\\n\";\n			}\n		}\n		$\
params{'outfile'} = undef if ( $jobIdForFilename =\
= 1 );\n	}\n	close $LISTFILE;\n}\n\n=head2 load_da\
ta()\n\nLoad sequence data from file or option spe\
cified on the command-line.\n\n  &load_data();\n\n\
=cut\n\nsub load_data {\n	print_debug_message( 'lo\
ad_data', 'Begin', 1 );\n	my $retSeq;\n\n	# Query \
sequence\n	if ( defined( $ARGV[0] ) ) {    # Bare \
option\n		if ( -f $ARGV[0] || $ARGV[0] eq '-' ) { \
   # File\n			$retSeq = &read_file( $ARGV[0] );\n	\
	}\n		else {                                     #\
 DB:ID or sequence\n			$retSeq = $ARGV[0];\n		}\n	\
}\n	if ( $params{'sequence'} ) {                  \
 # Via --sequence\n		if ( -f $params{'sequence'} |\
| $params{'sequence'} eq '-' ) {    # File\n			$re\
tSeq = &read_file( $params{'sequence'} );\n		}\n		\
else {    # DB:ID or sequence\n			$retSeq = $param\
s{'sequence'};\n		}\n	}\n	print_debug_message( 'lo\
ad_data', 'End', 1 );\n	return $retSeq;\n}\n\n=hea\
d2 load_params()\n\nLoad job parameters from comma\
nd-line options.\n\n  &load_params();\n\n=cut\n\ns\
ub load_params {\n	print_debug_message( 'load_para\
ms', 'Begin', 1 );\n\n	# Database(s) to search\n	m\
y (@dbList) = split /[ ,]/, $params{'database'};\n\
	$tool_params{'database'} = \\@dbList;\n\n	# Compa\
tability options, old command-line.\n	if(!$tool_pa\
rams{'viewfilter'} && $params{'echofilter'}) {\n		\
$tool_params{'viewfilter'} = 'true';\n	}\n	if(!$to\
ol_params{'alignments'} && $params{'numal'}) {\n		\
$tool_params{'alignments'} = $params{'numal'};\n	}\
\n	# TODO: set alignment format option to get NCBI\
 BLAST XML.\n	if($params{'appxml'}) {\n		$tool_par\
ams{'align'} = '';\n	}\n\n	print_debug_message( 'l\
oad_params', 'End', 1 );\n}\n\n=head2 client_poll(\
)\n\nClient-side job polling.\n\n  &client_poll($j\
ob_id);\n\n=cut\n\nsub client_poll {\n	print_debug\
_message( 'client_poll', 'Begin', 1 );\n	my $jobid\
  = shift;\n	my $status = 'PENDING';\n\n	my $error\
Count = 0;\n	while ($status eq 'RUNNING'\n		|| $st\
atus eq 'PENDING'\n		|| ( $status eq 'ERROR' && $e\
rrorCount < 2 ) )\n	{\n		$status = rest_get_status\
($jobid);\n		print STDERR \"$status\\n\" if ( $out\
putLevel > 0 );\n		if ( $status eq 'ERROR' ) {\n		\
	$errorCount++;\n		}\n		elsif ( $errorCount > 0 ) \
{\n			$errorCount--;\n		}\n		if (   $status eq 'RU\
NNING'\n			|| $status eq 'PENDING'\n			|| $status \
eq 'ERROR' )\n		{\n\n			# Wait before polling agai\
n.\n			sleep $checkInterval;\n		}\n	}\n	print_debu\
g_message( 'client_poll', 'End', 1 );\n	return $st\
atus;\n}\n\n=head2 get_results()\n\nGet the result\
s for a job identifier.\n\n  &get_results($job_id)\
;\n\n=cut\n\nsub get_results {\n	print_debug_messa\
ge( 'get_results', 'Begin', 1 );\n	my $jobid = shi\
ft;\n	print_debug_message( 'get_results', 'jobid: \
' . $jobid, 1 );\n\n	# Verbose\n	if ( $outputLevel\
 > 1 ) {\n		print 'Getting results for job ', $job\
id, \"\\n\";\n	}\n\n	# Check status, and wait if n\
ot finished\n	client_poll($jobid);\n\n	# Use JobId\
 if output file name is not defined\n	unless ( def\
ined( $params{'outfile'} ) ) {\n		$params{'outfile\
'} = $jobid;\n	}\n\n	# Get list of data types\n	my\
 (@resultTypes) = rest_get_result_types($jobid);\n\
\n	# Get the data and write it to a file\n	if ( de\
fined( $params{'outformat'} ) ) {    # Specified d\
ata type\n		my $selResultType;\n		foreach my $resu\
ltType (@resultTypes) {\n			if ( $resultType->{'id\
entifier'} eq $params{'outformat'} ) {\n				$selRe\
sultType = $resultType;\n			}\n		}\n		if ( defined\
($selResultType) ) {\n			my $result =\n			  rest_g\
et_result( $jobid, $selResultType->{'identifier'} \
);\n			if ( $params{'outfile'} eq '-' ) {\n				wri\
te_file( $params{'outfile'}, $result );\n			}\n			\
else {\n				write_file(\n					$params{'outfile'} .\
 '.'\n					  . $selResultType->{'identifier'} . '.\
'\n					  . $selResultType->{'fileSuffix'},\n					\
$result\n				);\n			}\n		}\n		else {\n			die 'Erro\
r: unknown result format \"' . $params{'outformat'\
} . '\"';\n		}\n	}\n	else {    # Data types availa\
ble\n		      # Write a file for each output type\n\
		for my $resultType (@resultTypes) {\n			if ( $ou\
tputLevel > 1 ) {\n				print STDERR 'Getting ', $r\
esultType->{'identifier'}, \"\\n\";\n			}\n			my $\
result = rest_get_result( $jobid, $resultType->{'i\
dentifier'} );\n			if ( $params{'outfile'} eq '-' \
) {\n				write_file( $params{'outfile'}, $result )\
;\n			}\n			else {\n				write_file(\n					$params{\
'outfile'} . '.'\n					  . $resultType->{'identifi\
er'} . '.'\n					  . $resultType->{'fileSuffix'},\\
n					$result\n				);\n			}\n		}\n	}\n	print_debug\
_message( 'get_results', 'End', 1 );\n}\n\n=head2 \
read_file()\n\nRead a file into a scalar. The spec\
ial filename '-' can be used to read from \nstanda\
rd input (STDIN).\n\n  my $data = &read_file($file\
name);\n\n=cut\n\nsub read_file {\n	print_debug_me\
ssage( 'read_file', 'Begin', 1 );\n	my $filename =\
 shift;\n	print_debug_message( 'read_file', 'filen\
ame: ' . $filename, 2 );\n	my ( $content, $buffer \
);\n	if ( $filename eq '-' ) {\n		while ( sysread(\
 STDIN, $buffer, 1024 ) ) {\n			$content .= $buffe\
r;\n		}\n	}\n	else {    # File\n		open( my $FILE, \
'<', $filename )\n		  or die \"Error: unable to op\
en input file $filename ($!)\";\n		while ( sysread\
( $FILE, $buffer, 1024 ) ) {\n			$content .= $buff\
er;\n		}\n		close($FILE);\n	}\n	print_debug_messag\
e( 'read_file', 'End', 1 );\n	return $content;\n}\\
n\n=head2 write_file()\n\nWrite data to a file. Th\
e special filename '-' can be used to write to \ns\
tandard output (STDOUT).\n\n  &write_file($filenam\
e, $data);\n\n=cut\n\nsub write_file {\n	print_deb\
ug_message( 'write_file', 'Begin', 1 );\n	my ( $fi\
lename, $data ) = @_;\n	print_debug_message( 'writ\
e_file', 'filename: ' . $filename, 2 );\n	if ( $ou\
tputLevel > 0 ) {\n		print STDERR 'Creating result\
 file: ' . $filename . \"\\n\";\n	}\n	if ( $filena\
me eq '-' ) {\n		print STDOUT $data;\n	}\n	else {\\
n		open( my $FILE, '>', $filename )\n		  or die \"\
Error: unable to open output file $filename ($!)\"\
;\n		syswrite( $FILE, $data );\n		close($FILE);\n	\
}\n	print_debug_message( 'write_file', 'End', 1 );\
\n}\n\n=head2 usage()\n\nPrint program usage messa\
ge.\n\n  &usage();\n\n=cut\n\nsub usage {\n	print \
STDERR <<EOF\nWU-BLAST\n========\n   \nRapid seque\
nce database search programs utilizing the BLAST a\
lgorithm\n    \n[Required]\n\n  -p, --program     \
 : str  : BLAST program to use, see --paramDetail \
program\n  -D, --database     : str  : database(s)\
 to search, space separated. See\n                \
              --paramDetail database\n      --styp\
e        : str  : query sequence type, see --param\
Detail stype\n  seqFile            : file : query \
sequence (\"-\" for STDIN, \\@filename for\n      \
                        identifier list file)\n\n[\
Optional]\n\n  -m, --matrix       : str  : scoring\
 matrix, see --paramDetail matrix\n  -e, --exp    \
      : real : 0<E<= 1000. Statistical significanc\
e threshold \n                              for re\
porting database sequence matches.\n  -e, --viewfi\
lter   :      : display the filtered query sequenc\
e\n  -f, --filter       : str  : filter the query \
sequence for low complexity \n                    \
          regions, see --paramDetail filter\n  -A,\
 --align        : int  : pairwise alignment format\
, see --paramDetail align\n  -s, --scores       : \
int  : number of scores to be reported\n  -b, --al\
ignments   : int  : number of alignments to report\
\n  -S, --sensitivity  : str  : sensitivity of the\
 search, \n                              see --par\
amDetail sensitivity\n  -t, --sort	     : str  : s\
ort order for hits, see --paramDetail sort\n  -T, \
--stats        : str  : statistical model, see --p\
aramDetail stats\n  -d, --strand       : str  : DN\
A strand to search with,\n                        \
      see --paramDetail strand\n  -c, --topcombon \
   : str  : consistent sets of HSPs\n      --multi\
fasta   :      : treat input as a set of fasta for\
matted sequences\n\n[General]\n\n  -h, --help     \
   :      : prints this help text\n      --async  \
     :      : forces to make an asynchronous query\
\n      --email       : str  : e-mail address\n   \
   --title       : str  : title for job\n      --s\
tatus      :      : get job status\n      --result\
Types :      : get available result types for job\\
n      --polljob     :      : poll for the status \
of a job\n      --jobid       : str  : jobid that \
was returned when an asynchronous job \n          \
                   was submitted.\n      --outfile\
     : str  : file name for results (default is jo\
bid;\n                             \"-\" for STDOU\
T)\n      --outformat   : str  : result format to \
retrieve\n      --params      :      : list input \
parameters\n      --paramDetail : str  : display d\
etails for input parameter\n      --quiet       : \
     : decrease output\n      --verbose     :     \
 : increase output\n      --trace       :      : s\
how SOAP messages being interchanged \n   \nSynchr\
onous job:\n\n  The results/errors are returned as\
 soon as the job is finished.\n  Usage: $scriptNam\
e --email <your\\@email> [options...] seqFile\n  R\
eturns: results as an attachment\n\nAsynchronous j\
ob:\n\n  Use this if you want to retrieve the resu\
lts at a later time. The results \n  are stored fo\
r up to 24 hours. 	\n  Usage: $scriptName --async \
--email <your\\@email> [options...] seqFile\n  Ret\
urns: jobid\n\n  Use the jobid to query for the st\
atus of the job. If the job is finished, \n  it al\
so returns the results/errors.\n  Usage: $scriptNa\
me --polljob --jobid <jobId> [--outfile string]\n \
 Returns: string indicating the status of the job \
and if applicable, results \n  as an attachment.\n\
\nFurther information:\n\n  http://www.ebi.ac.uk/T\
ools/webservices/services/sss/wu_blast_rest\n  htt\
p://www.ebi.ac.uk/Tools/webservices/tutorials/perl\
\n\nSupport/Feedback:\n\n  http://www.ebi.ac.uk/su\
pport/\nEOF\n}\n\n=head1 FEEDBACK/SUPPORT\n\nPleas\
e contact us at L<http://www.ebi.ac.uk/support/> i\
f you have any \nfeedback, suggestions or issues w\
ith the service or this client.\n\n=cut\n","\n\n\n\
my $PROBTRESH = 0.3;# base pairs below this prob t\
hreshold will be ignored\nmy $WEIGHT = 100.0; # fl\
oat!!\nmy $NUCALPH = \"ACGTUNRYMKSWHBVD\";\nuse va\
rs qw($NUCALPH $WEIGHT);\n\nmy $myname = basename(\
$0);\n\nuse strict;\nuse warnings;\n\nuse File::Ba\
sename;\nuse Getopt::Long;\nuse File::Glob ':glob'\
;\nuse File::Spec;\nuse File::Temp qw/ tempfile te\
mpdir /;\n\n\n\n\nsub tcoffeelib_header($;$)\n{\n \
   my ($nseq, $fd) = @_;\n    if (! defined($fd)) \
{\n        $fd = *STDOUT;\n    }\n    printf $fd \\
"! TC_LIB_FORMAT_01\\n\";\n    printf $fd \"%d\\n\\
", $nseq;\n}\n\n\nsub tcoffeelib_header_addseq($$;\
$)\n{\n    my ($id, $seq, $fd) = @_;\n    if (! de\
fined($fd)) {\n        $fd = *STDOUT;\n    }\n    \
printf $fd \"%s %d %s\\n\", $id, length($seq), $se\
q;\n}\n\n\nsub tcoffeelib_comment($;$)\n{\n    my \
($comment, $fd) = @_;\n    if (! defined($fd)) {\n\
        $fd = *STDOUT;\n    }\n    printf $fd \"!\\
" . $comment . \"\\n\";\n}\n\n\nsub tcoffeelib_str\
uct($$$;$)\n{\n    my ($nseq, $len, $bpm, $fd) = @\
_;\n\n    if (! defined($fd)) {\n        $fd = *ST\
DOUT;\n    }\n\n    # output basepair indices with\
 fixed weight\n    printf $fd \"#%d %d\\n\", $nseq\
, $nseq;\n    # output basepairs (only once) and w\
ith unit-offset\n    for (my $i=0; $i<$len; $i++) \
{\n        for (my $j=$i+1; $j<$len; $j++) {\n    \
        if (! defined($bpm->[$i][$j])) {\n        \
        print STDERR \"ERROR: \\$bpm->[$i][$j] und\
efined\\n\";\n            }\n            if ($bpm-\
>[$i][$j]>0) {\n                print $fd $i+1;\n \
               print $fd \" \";\n                p\
rint $fd $j+1;\n                print $fd \" \" . \
$bpm->[$i][$j] . \"\\n\";\n            }\n        \
}\n    }\n}\n\n\nsub tcoffeelib_footer(;$)\n{\n   \
 my ($fd) = @_;\n    if (! defined($fd)) {\n      \
  $fd = *STDOUT;\n    }\n    print $fd \"! SEQ_1_T\
O_N\\n\";\n}\n\n\n    \nsub plfold($$$)\n{    \n  \
  my ($id, $seq, $probtresh) = @_;\n    my (@struc\
t);# return\n    my ($templ, $fhtmp, $fnametmp, $c\
md, $ctr, $window_size);\n    our $ntemp++;\n    \\
n    $templ = $myname . \".pid-\" . $$ .$ntemp .\"\
.XXXXXX\";\n    ($fhtmp, $fnametmp) = tempfile($te\
mpl, UNLINK => 1); \n    print $fhtmp \">$id\\n$se\
q\\n\";\n\n    # --- init basepair array\n    #\n \
   for (my $i=0; $i<length($seq); $i++) {\n       \
 for (my $j=$i+1; $j<length($seq); $j++) {\n      \
      $struct[$i][$j]=0;\n        }\n    }\n\n\n  \
  # --- call rnaplfold and drop a readme\n    #\n \
   $window_size=(length($seq)<70)?length($seq):70;\
\n    $cmd = \"RNAplfold -W $window_size < $fnamet\
mp >/dev/null\";\n    system($cmd);\n    \n    if \
($? != 0) {\n        printf STDERR \"ERROR: RNAplf\
old ($cmd) exited with error status %d\\n\", $? >>\
 8;\n        return;\n    }\n    #unlink($fnametmp\
);\n    my $fps = sprintf(\"%s_dp.ps\", $id); # ch\
eck long name\n    \n    if (! -s $fps) {\n      {\
\n\n	$fps = sprintf(\"%s_dp.ps\", substr($id,0,12)\
); # check short name\n 	if (! -s $fps)\n	  {\n	  \
  die(\"couldn't find expected file $fps\\n\");\n	\
    return;\n	  }\n      }\n    }\n\n    \n    # -\
-- read base pairs from created postscript\n    #\\
n    open(FH, $fps);\n    while (my $line = <FH>) \
{\n        my ($nti, $ntj, $prob);\n        chomp(\
$line);        \n        # line: bp bp sqrt-prob u\
box\n        my @match = ($line =~ m/^([0-9]+) +([\
0-9]+) +([0-9\\.]+) +ubox$/);\n        if (scalar(\
@match)) {\n            $nti=$1;\n            $ntj\
=$2;\n            $prob=$3*$3;# prob stored as squ\
are root\n\n            if ($prob>$probtresh) {\n \
               #printf STDERR \"\\$struct[$nti][$n\
tj] sqrtprob=$3 prob=$prob > $probtresh\\n\";\n   \
             $struct[$nti-1][$ntj-1] = $WEIGHT\n  \
          }\n            # store with zero-offset\\
n        }\n    }\n    close(FH);\n\n    # remove \
or gzi postscript\n    #\n    unlink($fps);\n    #\
\n    # or gzip\n    #$cmd = \"gzip -qf $fps\";\n \
   #system($cmd);\n    #if ($? != 0) {\n    #    p\
rintf STDERR \"ERROR: gzip ($cmd) exited with erro\
r status %d\\n\", $? >> 8;\n    #}\n\n    return \\
\@struct;\n}\n\n\n\n\n\nsub rnaseqfmt($)\n{\n    m\
y ($seq) = @_;\n    # remove gaps\n    $seq =~ s/-\
//g;\n    # uppercase RNA\n    $seq = uc($seq);\n \
   # T -> U\n    $seq =~ s/T/U/g;\n    # check for\
 invalid charaters\n    $_ = $seq;\n    s/[^$NUCAL\
PH]//g;\n    return $_;\n}\n\n\n\n\nsub usage(;$)\\
n{    \n    my ($errmsg) = @_;\n    if ($errmsg) {\
\n        print STDERR \"ERROR: $errmsg\\n\";\n   \
 }\n    print STDERR << \"EOF\";\n$myname:\n Creat\
es a T-Coffee RNA structure library from RNAplfold\
 prediction.\n See FIXME:citation\nUsage:\n $mynam\
e -in seq_file -out tcoffee_lib\nEOF\n    exit(1);\
\n}\n\nsub read_fasta_seq \n  {\n    my $f=$_[0];\\
n    my %hseq;\n    my (@seq, @com, @name);\n    m\
y ($a, $s,$nseq);\n\n    open (F, $f);\n    while \
(<F>)\n      {\n	$s.=$_;\n      }\n    close (F);\\
n\n    \n    @name=($s=~/>(\\S*).*\\n[^>]*/g);\n  \
  \n    @seq =($s=~/>.*.*\\n([^>]*)/g);\n    @com \
=($s=~/>(\\S*)(.*)\\n([^>]*)/g);\n\n\n    $nseq=$#\
name+1;\n  \n    for ($a=0; $a<$nseq; $a++)\n     \
 {\n	my $n=$name[$a];\n	my $s;\n	$hseq{$n}{name}=$\
n;\n	$s=$seq[$a];$s=~s/\\s//g;\n	\n	$hseq{$n}{seq}\
=$s;\n	$hseq{$n}{com}=$com[$a];\n      }\n    retu\
rn %hseq;\n  }\n\n\n\n\n\n\n\nmy $fmsq = \"\";\nmy\
 $flib = \"\";\nmy %OPTS;\nmy %seq;\nmy ($id, $nse\
q, $i);\nmy @nl;\n\nGetOptions(\"in=s\" => \\$fmsq\
, \"out=s\" => \\$flib);\n\nif (! -s $fmsq) {\n   \
 usage(\"empty or non-existant file \\\"$fmsq\\\"\\
")\n}\nif (length($flib)==0) {\n    usage(\"empty \
out-filename\")\n}\n\n\n\n\n\n\n%seq=read_fasta_se\
q($fmsq);\n\n\n@nl=keys(%seq);\n\n$nseq=$#nl+1;\no\
pen FD_LIB, \">$flib\" or die \"can't open $flib!\\
";\ntcoffeelib_header($nseq, *FD_LIB);\nforeach $i\
d (keys (%seq))\n  {\n    my ($seq, $fmtseq);\n   \
 \n    $seq = $seq{$id}{seq};\n    \n    $fmtseq =\
 rnaseqfmt($seq);# check here, formatting for fold\
ing important later\n    if (length($seq)!=length(\
$fmtseq)) {\n        print STDERR \"ERROR: invalid\
 sequence $id is not an RNA sequence. read seq is:\
 $seq\\n\";\n        exit\n      }\n   \n    tcoff\
eelib_header_addseq($id, uc($seq), *FD_LIB);\n  }\\
ntcoffeelib_comment(\"generated by $myname on \" .\
 localtime(), *FD_LIB);\n\n\n\n$i=0;\nforeach $id \
(keys (%seq))\n  {\n    my ($cleanid, $seq, $bpm);\
\n    $seq=$seq{$id}{seq};\n    $cleanid = $id;\n \
   $cleanid =~ s,[/ ],_,g;# needed for rnaplfold\n\
    $seq = rnaseqfmt($seq);\n    \n    $bpm = plfo\
ld($cleanid, rnaseqfmt($seq), $PROBTRESH);       \\
n    \n    tcoffeelib_struct($i+1, length($seq), $\
bpm, *FD_LIB);\n    $i++;\n}\n\n\ntcoffeelib_foote\
r(*FD_LIB);\nclose FD_LIB;\nexit (0);\n\n","\n\n\n\
\n\n$cmd=join ' ', @ARGV;\nif ($cmd=~/-infile=(\\S\
+)/){ $seqfile=$1;}\nif ($cmd=~/-outfile=(\\S+)/){\
 $libfile=$1;}\n\n\n\n%s=read_fasta_seq ($seqfile)\
;\n\nopen (F, \">$libfile\");\nforeach $name (keys\
 (%s))\n  {\n    my $tclib=\"$name.RNAplfold_tclib\
\";\n    print (F \">$name _F_ $tclib\\n\");\n    \
seq2RNAplfold2tclib ($name, $s{$name}{seq}, $tclib\
);\n  }\nclose (F);\nexit (EXIT_SUCCESS);\n\nsub s\
eq2RNAplfold2tclib\n  {\n    my ($name, $seq, $tcl\
ib)=@_;\n    my ($tmp);\n    $n++;\n    $tmp=\"tmp\
4seq2RNAplfold_tclib.$$.$n.pep\";\n    open (RF, \\
">$tmp\");\n    print (RF \">$name\\n$seq\\n\");\n\
    close (RF);\n    \n    system \"t_coffee -othe\
r_pg RNAplfold2tclib.pl -in=$tmp -out=$tclib\";\n \
   \n    unlink ($tmp);\n    return $tclib;\n  }\n\
    \n    \nsub read_fasta_seq \n  {\n    my $f=@_\
[0];\n    my %hseq;\n    my (@seq, @com, @name);\n\
    my ($a, $s,$nseq);\n\n    open (F, $f);\n    w\
hile (<F>)\n      {\n	$s.=$_;\n      }\n    close \
(F);\n\n    \n    @name=($s=~/>(\\S*).*\\n[^>]*/g)\
;\n    \n    @seq =($s=~/>.*.*\\n([^>]*)/g);\n    \
@com =($s=~/>\\S*(.*)\\n([^>]*)/g);\n\n    \n    $\
nseq=$#name+1;\n    \n    for ($a=0; $a<$nseq; $a+\
+)\n      {\n	my $n=$name[$a];\n	$hseq{$n}{name}=$\
n;\n	$hseq{$n}{seq}=$seq[$a];\n	$hseq{$n}{com}=$co\
m[$a];\n      }\n    return %hseq;\n  }\n","use Ge\
topt::Long;\nuse File::Path;\nuse Env;\nuse FileHa\
ndle;\nuse Cwd;\nuse Sys::Hostname;\nour $PIDCHILD\
;\nour $ERROR_DONE;\nour @TMPFILE_LIST;\nour $EXIT\
_FAILURE=1;\nour $EXIT_SUCCESS=0;\n\nour $REFDIR=g\
etcwd;\nour $EXIT_SUCCESS=0;\nour $EXIT_FAILURE=1;\
\n\nour $PROGRAM=\"tc_generic_method.pl\";\nour $C\
L=$PROGRAM;\n\nour $CLEAN_EXIT_STARTED;\nour $debu\
g_lock=$ENV{\"DEBUG_LOCK\"};\nour $LOCKDIR=$ENV{\"\
LOCKDIR_4_TCOFFEE\"};\nif (!$LOCKDIR){$LOCKDIR=get\
cwd();}\nour $ERRORDIR=$ENV{\"ERRORDIR_4_TCOFFEE\"\
};\nour $ERRORFILE=$ENV{\"ERRORFILE_4_TCOFFEE\"};\\
n&set_lock ($$);\nif (isshellpid(getppid())){lock4\
tc(getppid(), \"LLOCK\", \"LSET\", \"$$\\n\");}\n \
     \nour $print;\nmy ($fmsq1, $fmsq2, $output, $\
outfile, $arch, $psv, $hmmtop_home, $trim, $cov, $\
sample, $mode, $gor_home, $gor_seq, $gor_obs);\n\n\
GetOptions(\"-in=s\" => \\$fmsq1,\"-output=s\" =>\\
\$output ,\"-out=s\" => \\$outfile, \"-arch=s\" =>\
 \\$arch,\"-psv=s\" => \\$psv, \"-hmmtop_home=s\",\
 \\$hmmtop_home,\"-trim=s\" =>\\$trim ,\"-print=s\\
" =>\\$print,\"-cov=s\" =>\\$cov , \"-sample=s\" =\
>\\$sample, \"-mode=s\" =>\\$mode, \"-gor_home=s\"\
=>\\$gor_home, \"-gor_seq=s\"=>\\$gor_seq,\"-gor_o\
bs=s\"=>\\$gor_obs);\n\n\nif (!$mode){$mode = \"hm\
mtop\"}\nelsif ($mode eq \"hmmtop\"){;}\nelsif ($m\
ode eq \"gor\"){;}\nelse {myexit(flush_error (\"-m\
ode=$mode is unknown\"));}\n\n\nour $HOME=$ENV{\"H\
OME\"};\nour $MCOFFEE=($ENV{\"MCOFFEE_4_TCOFFEE\"}\
)?$ENV{\"MCOFFEE_4_TCOFFEE\"}:\"$HOME/.t_coffee/mc\
offee\";\n\nif ($mode eq \"hmmtop\")\n  {\n    che\
ck_configuration (\"hmmtop\");\n    if (-e $arch){\
$ENV{'HMMTOP_ARCH'}=$arch;}\n    elsif (-e $ENV{HM\
MTOP_ARCH}){$arch=$ENV{HMMTOP_ARCH};}\n    elsif (\
-e \"$MCOFFEE/hmmtop.arch\"){$arch=$ENV{'HMMTOP_AR\
CH'}=\"$MCOFFEE/hmmtop.arch\";}\n    elsif (-e \"$\
hmmtop_home/hmmtop.arc\"){$arch=$ENV{'HMMTOP_ARCH'\
}=\"$hmmtop_home/hmmtop.arc\";}\n    else {myexit(\
flush_error ( \"Could not find ARCH file for hmmto\
p\"));}\n    \n    \n    if (-e $psv){$ENV{'HMMTOP\
_PSV'}=$psv;}\n    elsif (-e $ENV{HMMTOP_PSV}){$ps\
v=$ENV{HMMTOP_PSV};}\n    elsif (-e \"$MCOFFEE/hmm\
top.psv\"){$psv=$ENV{'HMMTOP_PSV'}=\"$MCOFFEE/hmmt\
op.psv\";}\n    elsif (-e \"$hmmtop_home/hmmtop.ps\
v\"){$psv=$ENV{'HMMTOP_PSV'}=\"$hmmtop_home/hmmtop\
.psv\";}\n    else {myexit(flush_error ( \"Could n\
ot find PSV file for hmmtop\"));}\n  }\nelsif ($mo\
de eq \"gor\")\n  {\n    our $GOR_SEQ;\n    our $G\
OR_OBS;\n    \n    check_configuration (\"gorIV\")\
;\n    if (-e $gor_seq){$GOR_SEQ=$gor_seq;}\n    e\
lsif (-e $ENV{GOR_SEQ}){$GOR_SEQ=$ENV{GOR_SEQ};}\n\
    elsif (-e \"$MCOFFEE/New_KS.267.seq\"){$GOR_SE\
Q=\"$MCOFFEE/New_KS.267.seq\";}\n    elsif (-e \"$\
gor_home/New_KS.267.seq\"){$GOR_SEQ=\"$gor_home/Ne\
w_KS.267.seq\";}\n    else {myexit(flush_error ( \\
"Could not find SEQ file for gor\"));}\n\n    if (\
-e $gor_obs){$GOR_OBS=$gor_obs;}\n    elsif (-e $E\
NV{GOR_OBS}){$GOR_OBS=$ENV{GOR_OBS};}\n    elsif (\
-e \"$MCOFFEE/New_KS.267.obs\"){$GOR_OBS=\"$MCOFFE\
E/New_KS.267.obs\";}\n    elsif (-e \"$gor_home/Ne\
w_KS.267.obs\"){$GOR_OBS=\"$gor_home/New_KS.267.ob\
s\";}\n    else {myexit(flush_error ( \"Could not \
find OBS file for gor\"));}\n  }\n\n\nif ( ! -e $f\
msq1){myexit(flush_error (\"Could Not Read Input f\
ile $fmsq1\"));}\n\n\nmy $fmsq2=vtmpnam();\nmy $fm\
sq3=vtmpnam();\nmy $tmpfile=vtmpnam();\nmy $predfi\
le=vtmpnam();\n\nif ($trim){$trim_action=\" +trim \
_aln_%%$trim\\_K1 \";}\nif ($cov) {$cov_action= \"\
 +sim_filter _aln_c$cov \";}\n&safe_system(\"t_cof\
fee -other_pg seq_reformat -in $fmsq1 -action +con\
vert 'BOUJXZ-' $cov_action $trim_action -output fa\
sta_aln -out $fmsq2\");\nmy (%pred, %seq, %predA);\
\n\n\n%seq=read_fasta_seq($fmsq2);\n%seq=fasta2sam\
ple(\\%seq, $sample);\n\nif (1==2 && $mode eq \"hm\
mtop\" && $output eq \"cons\")\n  {\n    fasta2hmm\
top_cons($outfile,\\%seq);\n  }\nelse\n  {\n    %p\
red=fasta2pred(\\%seq, $mode);\n    %predA=pred2al\
n (\\%pred, \\%seq);\n    \n    \n    if (!$output\
 || $output eq \"prediction\"){output_fasta_seq (\\
\%predA, $outfile);}\n    elsif ($output eq \"colo\
r_html\"){pred2color (\\%pred,\\%seq, $outfile);}\\
n    elsif ($output eq \"cons\"){pred2cons($outfil\
e,\\%predA);}\n    else {flush_error (\"$output is\
 an unknown output mode\");}\n  }\n\nsub fasta2sam\
ple\n  {\n    my $SR=shift;\n    my $it=shift;\n  \
  my %S=%$SR;\n    \n    my $seq=index2seq_name (\\
\%S, 1);\n    my $l=length($S{$seq}{seq});\n    my\
 @sl=keys(%S);\n    my $nseq=$#sl+1;\n    my $inde\
x=$nseq;\n  \n    if (!$sample) {return %S;}\n    \
for (my $a=0; $a<$it; $a++)\n      {\n	my $newseq=\
\"\";\n	my $nname=\"$seq\\_sampled_$index\";\n	for\
 (my $p=0; $p<$l; $p++)\n	  {\n	    my $i=int(rand\
($nseq));\n	    \n	    my $name = $sl[$i];\n	    m\
y $seq=$S{$name}{seq};\n	    my $r=substr ($seq, $\
p, 1);\n	    $newseq.=$r;\n	  }\n	$S{$nname}{name}\
=$nname;\n	$S{$nname}{seq}=$newseq;\n	$S{$nname}{c\
om}=\"sampled\";\n	$S{$nname}{index}=++$index;\n  \
    }\n    return %S;\n  }\n	      \nsub fasta2pre\
d\n  {\n    my $s=shift;\n    my $mode=shift;\n\n \
   if ( $mode eq \"hmmtop\"){return fasta2hmmtop_p\
red($s);}\n    elsif ($mode eq \"gor\"){return fas\
ta2gor_pred ($s);}\n  }\nsub fasta2hmmtop_cons\n  \
{\n    my $outfile=shift;\n    my $SR=shift;\n    \
\n    my $o = new FileHandle;\n    my $i = new Fil\
eHandle;\n    my $tmp_in =vtmpnam();\n    my $tmp_\
out=vtmpnam();\n    my %seq=%$SR;\n    my %pred;\n\
    my $N=keys(%seq);\n    \n    output_fasta_seq \
(\\%seq,$tmp_in, \"seq\");\n    `hmmtop -pi=mpred \
-if=$tmp_in -sf=FAS -pl 2>/dev/null >$tmp_out`;\n \
   open ($o, \">$outfile\");\n    open ($i, \"$tmp\
_out\");\n    while (<$i>)\n      {\n	my $l=$_;\n	\
if (($l=~/>HP\\:\\s+(\\d+)\\s+(.*)/)){my $line=\">\
$2 NSEQ: $N\\n\";print $o \"$line\";}\n	elsif ( ($\
l=~/.*pred(.*)/))  {my $line=\"$1\\n\";print $o \"\
$line\";}\n      }\n    close ($o);\n    close ($i\
);\n    return read_fasta_seq($tmp);\n  }\nsub fas\
ta2hmmtop_pred\n  {\n    my $SR=shift;\n    my $o \
= new FileHandle;\n    my $i = new FileHandle;\n  \
  my $tmp    =vtmpnam();\n    my $tmp_in =vtmpnam(\
);\n    my $tmp_out=vtmpnam();\n    my %seq=%$SR;\\
n    my %pred;\n    \n\n    output_fasta_seq (\\%s\
eq,$tmp_in, \"seq\");\n    `hmmtop -if=$tmp_in -sf\
=FAS -pl 2>/dev/null >$tmp_out`;\n    open ($o, \"\
>$tmp\");\n    open ($i, \"$tmp_out\");\n    while\
 (<$i>)\n      {\n	my $l=$_;\n	if (($l=~/>HP\\:\\s\
+(\\d+)\\s+(.*)/)){my $line=\">$2\\n\";print $o \"\
$line\";}\n	elsif ( ($l=~/.*pred(.*)/))  {my $line\
=\"$1\\n\";print $o \"$line\";}\n      }\n    clos\
e ($o);\n    close ($i);\n    return read_fasta_se\
q($tmp);\n  }\n    \n	\n	\n	    \n	\n	\n\n	\nsub f\
asta2gor_pred\n  {\n    my $SR=shift;\n    my $o =\
 new FileHandle;\n    my $i = new FileHandle;\n   \
 my $tmp    =vtmpnam();\n    my $tmp_in =vtmpnam()\
;\n    my $tmp_out=vtmpnam();\n    my %seq=%$SR;\n\
    my %pred;\n    \n\n    output_fasta_seq (\\%se\
q,$tmp_in, \"seq\");\n    `gorIV -prd $tmp_in -seq\
 $GOR_SEQ -obs $GOR_OBS >$tmp_out`;\n    open ($o,\
 \">$tmp\");\n    open ($i, \"$tmp_out\");\n    wh\
ile (<$i>)\n      {\n	my $l=$_;\n\n	\n	if ( $l=~/>\
/){print $o \"$l\";}\n	elsif ( $l=~/Predicted Sec.\
 Struct./){$l=~s/Predicted Sec. Struct\\.//;print \
$o \"$l\";}\n      }\n    close ($o);\n    close (\
$i);\n    return read_fasta_seq($tmp);\n  }\n			\n\
			     \nsub index2seq_name\n  {\n    \n    my $S\
R=shift;\n    my $index=shift;\n    \n    \n    my\
 %S=%$SR;\n    \n    foreach my $s (%S)\n      {\n\
	if ( $S{$s}{index}==$index){return $s;}\n      }\\
n    return \"\";\n  }\n\nsub pred2cons\n  {\n    \
my $outfile=shift;\n    my $predR=shift;\n    my $\
seq=shift;\n    my %P=%$predR;\n    my %C;\n    my\
 ($s,@r,$nseq);\n    my $f= new FileHandle;\n\n   \
 open ($f, \">$outfile\");\n\n    if (!$seq){$seq=\
index2seq_name(\\%P,1);}\n    foreach my $s (keys(\
%P))\n      {\n	$nseq++;\n	$string= $P{$s}{seq};\n\
	$string = uc $string;\n	my @r=split (//,$string);\
\n	for (my $a=0; $a<=$#r; $a++)\n	  {\n	    if (($\
r[$a]=~/[OHICE]/)){$C{$a}{$r[$a]}++;}\n	  }\n     \
 }\n    @l=keys(%C);\n    \n    \n    $s=$P{$seq}{\
seq};\n    print $f \">$seq pred based on $nseq\\n\
\";\n    @r=split (//,$s);\n    \n    for (my $x=0\
; $x<=$#r; $x++)\n      {\n	if ($r[$x] ne \"-\")\n\
	  {\n	    my $h=$C{$x}{H};\n	    my $i=$C{$x}{I};\
\n	    my $o=$C{$x}{O};\n	    my $c=$C{$x}{C};\n	 \
   my $e=$C{$x}{E};\n	    my $l=$i+$o;\n	    \n	  \
  if ($h>=$i && $h>=$o && $h>=$c && $h>=$e){$r[$x]\
='H';}\n	    elsif ($i>=$o && $i>=$c && $i>=$e){$r\
[$x]='I';}\n	    elsif ($o>=$c && $o>=$e){$r[$x]='\
O';}\n	    elsif ($c>=$e){$r[$x]='C';}\n	    else \
{$r[$x]='E';}\n	  }\n      }\n    $j=join ('', @r)\
;\n    print $f \"$j\\n\";\n    close ($f);\n    r\
eturn $j;\n  }\n\nsub pred2aln\n  {\n    my $PR=sh\
ift;\n    my $AR=shift;\n    \n    my $f=new FileH\
andle;\n    my %P=%$PR;\n    my %A=%$AR;\n    my %\
PA;\n    my $tmp=vtmpnam();\n    my $f= new FileHa\
ndle;\n    \n    open ($f, \">$tmp\");\n    foreac\
h my $s (sort{$A{$a}{index}<=>$A{$b}{index}}(keys \
(%A)))\n      {\n	my (@list, $seq, @plist, @pseq, \
$L, $PL, $c, $w);\n	my $seq;\n	my $seq=$A{$s}{seq}\
;\n	my $pred=$P{$s}{seq};\n	$seq=pred2alnS($P{$s}{\
seq},$A{$s}{seq});\n	print $f \">$s\\n$seq\\n\";\n\
      }\n    close ($f);\n    return read_fasta_se\
q ($tmp);\n  }\nsub pred2alnS\n  {\n    my $pred=s\
hift;\n    my $aln= shift;\n    my ($j,$a,$b);\n  \
  my @P=split (//, $pred);\n    my @A=split (//, $\
aln);\n    for ($a=$b=0;$a<=$#A; $a++)\n      {\n	\
if ($A[$a] ne \"-\"){$A[$a]=$P[$b++];}\n      }\n \
   if ($b!= ($#P+1)){add_warning (\"Could not thre\
ad sequence: $b $#P\");}\n    \n    $j= join ('', \
@A);\n    return $j;\n  }\nsub pred2color\n  {\n  \
  my $predP=shift;\n    my $alnP=shift;\n    my $o\
ut=shift;\n    my $F=new FileHandle;\n    my $stru\
c=vtmpnam();\n    my $aln=vtmpnam();\n    \n\n    \
output_fasta_seq ($alnP, $aln);\n    my %p=%$predP\
;\n    \n    open ($F, \">$struc\");\n    \n    \n\
    foreach my $s (keys(%p))\n      {\n	\n	print $\
F \">$s\\n\";\n	my $s=uc($p{$s}{seq});\n	\n	$s=~s/\
[Oo]/0/g;\n	$s=~s/[Ee]/0/g;\n	\n	$s=~s/[Ii]/5/g;\n\
	$s=~s/[Cc]/5/g;\n	\n	$s=~s/[Hh]/9/g;\n	\n	print $\
F \"$s\\n\";\n      }\n    close ($F);\n    \n    \
\n    \n    safe_system ( \"t_coffee -other_pg seq\
_reformat -in $aln -struc_in $struc -struc_in_f nu\
mber_fasta -output color_html -out $out\");\n    r\
eturn;\n  }\n	  \n    \nsub display_fasta_seq\n  {\
\n    my $SR=shift;\n    my %S=%$SR;\n    \n    fo\
reach my $s (sort{$S{$a}{index}<=>$S{$b}{index}}(k\
eys (%S)))\n      {\n	print STDERR \">$s\\n$S{$s}{\
seq}\\n\";\n      }\n    close ($f);\n  }\nsub out\
put_fasta_seq\n  {\n    my $SR=shift;\n    my $out\
file=shift;\n    my $mode =shift;\n    my $f= new \
FileHandle;\n    my %S=%$SR;\n    \n    \n    open\
 ($f, \">$outfile\");\n    foreach my $s (sort{$S{\
$a}{index}<=>$S{$b}{index}}(keys (%S)))\n      {\n\
	my $seq=$S{$s}{seq};\n	if ( $mode eq \"seq\"){$se\
q=~s/\\-//g;}\n	print $f \">$s\\n$seq\\n\";\n     \
 }\n    close ($f);\n  }\n      \nsub read_fasta_s\
eq \n  {\n    my $f=$_[0];\n    my %hseq;\n    my \
(@seq, @com, @name);\n    my ($a, $s,$nseq);\n    \
my $index;\n    open (F, $f);\n    while (<F>)\n  \
    {\n	$s.=$_;\n      }\n    close (F);\n\n    \n\
    @name=($s=~/>(\\S*).*\\n[^>]*/g);\n    \n    @\
seq =($s=~/>.*.*\\n([^>]*)/g);\n    @com =($s=~/>.\
*(.*)\\n([^>]*)/g);\n\n\n    $nseq=$#name+1;\n    \
\n  \n    for ($a=0; $a<$nseq; $a++)\n      {\n	my\
 $n=$name[$a];\n	my $s;\n	$hseq{$n}{name}=$n;\n	$s\
=$seq[$a];$s=~s/\\s//g;\n	$hseq{$n}{index}=++$inde\
x;\n	$hseq{$n}{seq}=$s;\n	$hseq{$n}{com}=$com[$a];\
\n      }\n    return %hseq;\n  }\n\n\nsub file2he\
ad\n      {\n	my $file = shift;\n	my $size = shift\
;\n	my $f= new FileHandle;\n	my $line;\n	open ($f,\
$file);\n	read ($f,$line, $size);\n	close ($f);\n	\
return $line;\n      }\nsub file2tail\n      {\n	m\
y $file = shift;\n	my $size = shift;\n	my $f= new \
FileHandle;\n	my $line;\n	\n	open ($f,$file);\n	se\
ek ($f,$size*-1, 2);\n	read ($f,$line, $size);\n	c\
lose ($f);\n	return $line;\n      }\n\n\nsub vtmpn\
am\n      {\n	my $r=rand(100000);\n	my $f=\"file.$\
r.$$\";\n	while (-e $f)\n	  {\n	    $f=vtmpnam();\\
n	  }\n	push (@TMPFILE_LIST, $f);\n	return $f;\n  \
    }\n\nsub myexit\n  {\n    my $code=@_[0];\n   \
 if ($CLEAN_EXIT_STARTED==1){return;}\n    else {$\
CLEAN_EXIT_STARTED=1;}\n    ### ONLY BARE EXIT\n  \
  exit ($code);\n  }\nsub set_error_lock\n    {\n \
     my $name = shift;\n      my $pid=$$;\n\n     \
 \n      &lock4tc ($$,\"LERROR\", \"LSET\", \"$$ -\
- ERROR: $name $PROGRAM\\n\");\n      return;\n   \
 }\nsub set_lock\n  {\n    my $pid=shift;\n    my \
$msg= shift;\n    my $p=getppid();\n    &lock4tc (\
$pid,\"LLOCK\",\"LRESET\",\"$p$msg\\n\");\n  }\nsu\
b unset_lock\n   {\n     \n    my $pid=shift;\n   \
 &lock4tc ($pid,\"LLOCK\",\"LRELEASE\",\"\");\n  }\
\nsub shift_lock\n  {\n    my $from=shift;\n    my\
 $to=shift;\n    my $from_type=shift;\n    my $to_\
type=shift;\n    my $action=shift;\n    my $msg;\n\
    \n    if (!&lock4tc($from, $from_type, \"LCHEC\
K\", \"\")){return 0;}\n    $msg=&lock4tc ($from, \
$from_type, \"LREAD\", \"\");\n    &lock4tc ($from\
, $from_type,\"LRELEASE\", $msg);\n    &lock4tc ($\
to, $to_type, $action, $msg);\n    return;\n  }\ns\
ub isshellpid\n  {\n    my $p=shift;\n    if (!loc\
k4tc ($p, \"LLOCK\", \"LCHECK\")){return 0;}\n    \
else\n      {\n	my $c=lock4tc($p, \"LLOCK\", \"LRE\
AD\");\n	if ( $c=~/-SHELL-/){return 1;}\n      }\n\
    return 0;\n  }\nsub isrootpid\n  {\n    if(loc\
k4tc (getppid(), \"LLOCK\", \"LCHECK\")){return 0;\
}\n    else {return 1;}\n  }\nsub lock4tc\n	{\n	  \
my ($pid,$type,$action,$value)=@_;\n	  my $fname;\\
n	  my $host=hostname;\n	  \n	  if ($type eq \"LLO\
CK\"){$fname=\"$LOCKDIR/.$pid.$host.lock4tcoffee\"\
;}\n	  elsif ( $type eq \"LERROR\"){ $fname=\"$LOC\
KDIR/.$pid.$host.error4tcoffee\";}\n	  elsif ( $ty\
pe eq \"LWARNING\"){ $fname=\"$LOCKDIR/.$pid.$host\
.warning4tcoffee\";}\n	  \n	  if ($debug_lock)\n	 \
   {\n	      print STDERR \"\\n\\t---lock4tc(tcg):\
 $action => $fname =>$value (RD: $LOCKDIR)\\n\";\n\
	    }\n\n	  if    ($action eq \"LCHECK\") {return\
 -e $fname;}\n	  elsif ($action eq \"LREAD\"){retu\
rn file2string($fname);}\n	  elsif ($action eq \"L\
SET\") {return string2file ($value, $fname, \">>\"\
);}\n	  elsif ($action eq \"LRESET\") {return stri\
ng2file ($value, $fname, \">\");}\n	  elsif ($acti\
on eq \"LRELEASE\") \n	    {\n	      if ( $debug_l\
ock)\n		{\n		  my $g=new FileHandle;\n		  open ($g\
, \">>$fname\");\n		  print $g \"\\nDestroyed by $\
$\\n\";\n		  close ($g);\n		  safe_system (\"mv $f\
name $fname.old\");\n		}\n	      else\n		{\n		  un\
link ($fname);\n		}\n	    }\n	  return \"\";\n	}\n\
	\nsub file2string\n	{\n	  my $file=@_[0];\n	  my \
$f=new FileHandle;\n	  my $r;\n	  open ($f, \"$fil\
e\");\n	  while (<$f>){$r.=$_;}\n	  close ($f);\n	\
  return $r;\n	}\nsub string2file \n    {\n    my \
($s,$file,$mode)=@_;\n    my $f=new FileHandle;\n \
   \n    open ($f, \"$mode$file\");\n    print $f \
 \"$s\";\n    close ($f);\n  }\n\nBEGIN\n    {\n  \
    srand;\n    \n      $SIG{'SIGUP'}='signal_clea\
nup';\n      $SIG{'SIGINT'}='signal_cleanup';\n   \
   $SIG{'SIGQUIT'}='signal_cleanup';\n      $SIG{'\
SIGILL'}='signal_cleanup';\n      $SIG{'SIGTRAP'}=\
'signal_cleanup';\n      $SIG{'SIGABRT'}='signal_c\
leanup';\n      $SIG{'SIGEMT'}='signal_cleanup';\n\
      $SIG{'SIGFPE'}='signal_cleanup';\n      \n  \
    $SIG{'SIGKILL'}='signal_cleanup';\n      $SIG{\
'SIGPIPE'}='signal_cleanup';\n      $SIG{'SIGSTOP'\
}='signal_cleanup';\n      $SIG{'SIGTTIN'}='signal\
_cleanup';\n      $SIG{'SIGXFSZ'}='signal_cleanup'\
;\n      $SIG{'SIGINFO'}='signal_cleanup';\n      \
\n      $SIG{'SIGBUS'}='signal_cleanup';\n      $S\
IG{'SIGALRM'}='signal_cleanup';\n      $SIG{'SIGTS\
TP'}='signal_cleanup';\n      $SIG{'SIGTTOU'}='sig\
nal_cleanup';\n      $SIG{'SIGVTALRM'}='signal_cle\
anup';\n      $SIG{'SIGUSR1'}='signal_cleanup';\n\\
n\n      $SIG{'SIGSEGV'}='signal_cleanup';\n      \
$SIG{'SIGTERM'}='signal_cleanup';\n      $SIG{'SIG\
CONT'}='signal_cleanup';\n      $SIG{'SIGIO'}='sig\
nal_cleanup';\n      $SIG{'SIGPROF'}='signal_clean\
up';\n      $SIG{'SIGUSR2'}='signal_cleanup';\n\n \
     $SIG{'SIGSYS'}='signal_cleanup';\n      $SIG{\
'SIGURG'}='signal_cleanup';\n      $SIG{'SIGCHLD'}\
='signal_cleanup';\n      $SIG{'SIGXCPU'}='signal_\
cleanup';\n      $SIG{'SIGWINCH'}='signal_cleanup'\
;\n      \n      $SIG{'INT'}='signal_cleanup';\n  \
    $SIG{'TERM'}='signal_cleanup';\n      $SIG{'KI\
LL'}='signal_cleanup';\n      $SIG{'QUIT'}='signal\
_cleanup';\n      \n      our $debug_lock=$ENV{\"D\
EBUG_LOCK\"};\n      \n      \n      \n      \n   \
   foreach my $a (@ARGV){$CL.=\" $a\";}\n      if \
( $debug_lock ){print STDERR \"\\n\\n\\n**********\
 START PG: $PROGRAM *************\\n\";}\n      if\
 ( $debug_lock ){print STDERR \"\\n\\n\\n*********\
*(tcg) LOCKDIR: $LOCKDIR $$ *************\\n\";}\n\
      if ( $debug_lock ){print STDERR \"\\n --- $$\
 -- $CL\\n\";}\n      \n	     \n      \n      \n  \
  }\nsub flush_error\n  {\n    my $msg=shift;\n   \
 return add_error ($EXIT_FAILURE,$$, $$,getppid(),\
 $msg, $CL);\n  }\nsub add_error \n  {\n    my $co\
de=shift;\n    my $rpid=shift;\n    my $pid=shift;\
\n    my $ppid=shift;\n    my $type=shift;\n    my\
 $com=shift;\n    \n    $ERROR_DONE=1;\n    lock4t\
c ($rpid, \"LERROR\",\"LSET\",\"$pid -- ERROR: $ty\
pe\\n\");\n    lock4tc ($$, \"LERROR\",\"LSET\", \\
"$pid -- COM: $com\\n\");\n    lock4tc ($$, \"LERR\
OR\",\"LSET\", \"$pid -- STACK: $ppid -> $pid\\n\"\
);\n   \n    return $code;\n  }\nsub add_warning \\
n  {\n    my $rpid=shift;\n    my $pid =shift;\n  \
  my $command=shift;\n    my $msg=\"$$ -- WARNING:\
 $command\\n\";\n    print STDERR \"$msg\";\n    l\
ock4tc ($$, \"LWARNING\", \"LSET\", $msg);\n  }\n\\
nsub signal_cleanup\n  {\n    print dtderr \"\\n**\
** $$ (tcg) was killed\\n\";\n    &cleanup;\n    e\
xit ($EXIT_FAILURE);\n  }\nsub clean_dir\n  {\n   \
 my $dir=@_[0];\n    if ( !-d $dir){return ;}\n   \
 elsif (!($dir=~/tmp/)){return ;}#safety check 1\n\
    elsif (($dir=~/\\*/)){return ;}#safety check 2\
\n    else\n      {\n	`rm -rf $dir`;\n      }\n   \
 return;\n  }\nsub cleanup\n  {\n    #print stderr\
 \"\\n----tc: $$ Kills $PIDCHILD\\n\";\n    #kill \
(SIGTERM,$PIDCHILD);\n    my $p=getppid();\n    $C\
LEAN_EXIT_STARTED=1;\n    \n    \n    \n    if (&l\
ock4tc($$,\"LERROR\", \"LCHECK\", \"\"))\n      {\\
n	my $ppid=getppid();\n	if (!$ERROR_DONE) \n	  {\n\
	    &lock4tc($$,\"LERROR\", \"LSET\", \"$$ -- STA\
CK: $p -> $$\\n\");\n	    &lock4tc($$,\"LERROR\", \
\"LSET\", \"$$ -- COM: $CL\\n\");\n	  }\n      }\n\
    my $warning=&lock4tc($$, \"LWARNING\", \"LREAD\
\", \"\");\n    my $error=&lock4tc($$,  \"LERROR\"\
, \"LREAD\", \"\");\n    #release error and warnin\
g lock if root\n    \n    if (isrootpid() && ($war\
ning || $error) )\n      {\n	\n	print STDERR \"***\
************* Summary *************\\n$error\\n$wa\
rning\\n\";\n\n	&lock4tc($$,\"LERROR\",\"RELEASE\"\
,\"\");\n	&lock4tc($$,\"LWARNING\",\"RELEASE\",\"\\
");\n      } \n    \n    \n    foreach my $f (@TMP\
FILE_LIST)\n      {\n	if (-e $f){unlink ($f);} \n \
     }\n    foreach my $d (@TMPDIR_LIST)\n      {\\
n	clean_dir ($d);\n      }\n    #No More Lock Rele\
ase\n    #&lock4tc($$,\"LLOCK\",\"LRELEASE\",\"\")\
; #release lock \n\n    if ( $debug_lock ){print S\
TDERR \"\\n\\n\\n********** END PG: $PROGRAM ($$) \
*************\\n\";}\n    if ( $debug_lock ){print\
 STDERR \"\\n\\n\\n**********(tcg) LOCKDIR: $LOCKD\
IR $$ *************\\n\";}\n  }\nEND \n  {\n    \n\
    &cleanup();\n  }\n   \n\nsub safe_system \n{\n\
  my $com=shift;\n  my $ntry=shift;\n  my $ctry=sh\
ift;\n  my $pid;\n  my $status;\n  my $ppid=getppi\
d();\n  if ($com eq \"\"){return 1;}\n  \n  \n\n  \
if (($pid = fork ()) < 0){return (-1);}\n  if ($pi\
d == 0)\n    {\n      set_lock($$, \" -SHELL- $com\
 (tcg)\");\n      exec ($com);\n    }\n  else\n   \
 {\n      lock4tc ($$, \"LLOCK\", \"LSET\", \"$pid\
\\n\");#update parent\n      $PIDCHILD=$pid;\n    \
}\n  if ($debug_lock){printf STDERR \"\\n\\t .... \
safe_system (fasta_seq2hmm)  p: $$ c: $pid COM: $c\
om\\n\";}\n\n  waitpid ($pid,WTERMSIG);\n\n  shift\
_lock ($pid,$$, \"LWARNING\",\"LWARNING\", \"LSET\\
");\n\n  if ($? == $EXIT_FAILURE || lock4tc($pid, \
\"LERROR\", \"LCHECK\", \"\"))\n    {\n      if ($\
ntry && $ctry <$ntry)\n	{\n	  add_warning ($$,$$,\\
"$com failed [retry: $ctry]\");\n	  lock4tc ($pid,\
 \"LRELEASE\", \"LERROR\", \"\");\n	  return safe_\
system ($com, $ntry, ++$ctry);\n	}\n      elsif ($\
ntry == -1)\n	{\n	  if (!shift_lock ($pid, $$, \"L\
ERROR\", \"LWARNING\", \"LSET\"))\n	    {\n	      \
add_warning ($$,$$,\"$com failed\");\n	    }\n	  e\
lse\n	    {\n	      lock4tc ($pid, \"LRELEASE\", \\
"LERROR\", \"\");\n	    }\n	  return $?;}\n      e\
lse\n	{\n	  if (!shift_lock ($pid,$$, \"LERROR\",\\
"LERROR\", \"LSET\"))\n	    {\n	      myexit(add_e\
rror ($EXIT_FAILURE,$$,$pid,getppid(), \"UNSPECIFI\
ED system\", $com));\n	    }\n	}\n    }\n  return \
$?;\n}\n\nsub check_configuration \n    {\n      m\
y @l=@_;\n      my $v;\n      foreach my $p (@l)\n\
	{\n	  \n	  if   ( $p eq \"EMAIL\")\n	    { \n	   \
   if ( !($EMAIL=~/@/))\n		{\n		add_warning($$,$$,\
\"Could Not Use EMAIL\");\n		myexit(add_error ($EX\
IT_FAILURE,$$,$$,getppid(),\"EMAIL\",\"$CL\"));\n	\
      }\n	    }\n	  elsif( $p eq \"INTERNET\")\n	 \
   {\n	      if ( !&check_internet_connection())\n\
		{\n		  myexit(add_error ($EXIT_FAILURE,$$,$$,get\
ppid(),\"INTERNET\",\"$CL\"));\n		}\n	    }\n	  el\
sif( $p eq \"wget\")\n	    {\n	      if (!&pg_is_i\
nstalled (\"wget\") && !&pg_is_installed (\"curl\"\
))\n		{\n		  myexit(add_error ($EXIT_FAILURE,$$,$$\
,getppid(),\"PG_NOT_INSTALLED:wget\",\"$CL\"));\n	\
	}\n	    }\n	  elsif( !(&pg_is_installed ($p)))\n	\
    {\n	      myexit(add_error ($EXIT_FAILURE,$$,$\
$,getppid(),\"PG_NOT_INSTALLED:$p\",\"$CL\"));\n	 \
   }\n	}\n      return 1;\n    }\nsub pg_is_instal\
led\n  {\n    my @ml=@_;\n    my $r, $p, $m;\n    \
my $supported=0;\n    \n    my $p=shift (@ml);\n  \
  if ($p=~/::/)\n      {\n	if (safe_system (\"perl\
 -M$p -e 1\")==$EXIT_SUCCESS){return 1;}\n	else {r\
eturn 0;}\n      }\n    else\n      {\n	$r=`which \
$p 2>/dev/null`;\n	if ($r eq \"\"){return 0;}\n	el\
se {return 1;}\n      }\n  }\n\n\n\nsub check_inte\
rnet_connection\n  {\n    my $internet;\n    my $t\
mp;\n    &check_configuration ( \"wget\"); \n    \\
n    $tmp=&vtmpnam ();\n    \n    if     (&pg_is_i\
nstalled    (\"wget\")){`wget www.google.com -O$tm\
p >/dev/null 2>/dev/null`;}\n    elsif  (&pg_is_in\
stalled    (\"curl\")){`curl www.google.com -o$tmp\
 >/dev/null 2>/dev/null`;}\n    \n    if ( !-e $tm\
p || -s $tmp < 10){$internet=0;}\n    else {$inter\
net=1;}\n    if (-e $tmp){unlink $tmp;}\n\n    ret\
urn $internet;\n  }\nsub check_pg_is_installed\n  \
{\n    my @ml=@_;\n    my $r=&pg_is_installed (@ml\
);\n    if (!$r && $p=~/::/)\n      {\n	print STDE\
RR \"\\nYou Must Install the perl package $p on yo\
ur system.\\nRUN:\\n\\tsudo perl -MCPAN -e 'instal\
l $pg'\\n\";\n      }\n    elsif (!$r)\n      {\n	\
myexit(flush_error(\"\\nProgram $p Supported but N\
ot Installed on your system\"));\n      }\n    els\
e\n      {\n	return 1;\n      }\n  }\n\n\n\n","\n\\
n\n\n\nmy $FMODEL =\"\"; \nmy $TMPDIR = \"/tmp\";\\
n\n\n\n\nmy $NUCALPH = \"ACGTUNRYMKSWHBVD\";\nmy $\
PRIMNUCALPH = \"ACGTUN\";\nuse vars qw($NUCALPH $P\
RIMNUCALPH $TMPDIR);\n\n\nmy $errmsg;\nuse vars qw\
($errmsg);\n\n\n\nuse Getopt::Long;\nuse Cwd;\nuse\
 File::Basename;\nuse File::Temp qw/ tempfile temp\
dir /;\nuse File::Copy;\nuse File::Path;\n\n\n\nsu\
b usage(;$)\n{\n    my ($errmsg) = @_;\n    my $my\
name = basename($0);\n\n    if ($errmsg) {\n      \
  print STDERR \"ERROR: $errmsg\\n\";\n    }\n\n  \
  print STDERR << \"EOF\";\n    \n$myname: align t\
wo sequences by means of consan\\'s sfold\nUsage:\\
n $myname -i file -o file -d path\nOptions:\n -i|-\
-in : pairwise input sequence file\n -o|--out: out\
put alignment\n -d|--directory containing data\n\n\
EOF\n}\n\nsub read_stk_aln \n  {\n    my $f=$_[0];\
\n    my ($seq, $id);\n    \n    my %hseq;\n\n    \
open (STK, \"$f\");\n    while (<STK>)\n      {\n	\
if ( /^#/ || /^\\/\\// || /^\\s*$/){;}\n	else\n	  \
{\n	    ($id,$seq)=/(\\S+)\\s+(\\S+)/;\n	    $hseq\
{$id}{'seq'}.=$seq;\n	  }\n      }\n    close (STK\
);\n    return %hseq;\n  }\nsub read_fasta_seq \n \
 {\n    my $f=$_[0];\n    my %hseq;\n    my (@seq,\
 @com, @name);\n    my ($a, $s,$nseq);\n\n    open\
 (F, $f);\n    while (<F>)\n      {\n	$s.=$_;\n   \
   }\n    close (F);\n\n    \n    @name=($s=~/>(.*\
).*\\n[^>]*/g);\n    \n    @seq =($s=~/>.*.*\\n([^\
>]*)/g);\n    @com =($s=~/>.*(.*)\\n([^>]*)/g);\n\\
n    \n    $nseq=$#name+1;\n    \n    for ($a=0; $\
a<$nseq; $a++)\n      {\n	my $n=$name[$a];\n	$hseq\
{$n}{name}=$n;\n	$hseq{$n}{seq}=$seq[$a];\n	$hseq{\
$n}{com}=$com[$a];\n      }\n    return %hseq;\n  \
}\n\n\n\nsub sfold_parseoutput($$)\n{\n    my ($fr\
awout, $foutfa) = @_;\n    my %haln;\n    my ($fst\
k, $cmd, $id);\n    open FOUTFA, \">$foutfa\";\n  \
  \n    $fstk = $frawout . \".stk\";\n    \n    # \
first line of raw out contains info\n    # remaini\
ng stuff is stockholm formatted\n    $cmd = \"sed \
-e '1d' $frawout\";\n    system(\"$cmd > $fstk\");\
\n    if ($? != 0) {\n        $errmsg = \"command \
failed with exit status $?.\";\n        $errmsg .=\
  \"Command was \\\"$cmd\\\"\";\n        return -1\
;\n    }\n\n    # this gives an error message. jus\
t ignore it...\n    %haln=read_stk_aln ( $fstk);\n\
    foreach $i (keys (%haln))\n      {\n	my $s;\n	\
$s=$haln{$i}{'seq'};\n	$s =~ s/\\./-/g;\n	print FO\
UTFA \">$i\\n$s\\n\";\n      }\n    close FOUTFA;\\
n    return 0;\n}\n\n\n\n\nsub sfold_wrapper($$$$)\
\n{\n    \n    my ($fs1, $fs2, $fmodel, $foutfa) =\
 @_;\n    \n\n    my ($cmd, $frawout, $ferrlog, $f\
readme, $ftimelog, $fstk);\n\n    # add  basename(\
$fmsqin) (unknown here!)\n    $frawout = \"sfold.l\
og\";\n    $ferrlog = \"sfold.err\";\n    $ftimelo\
g = \"sfold.time\";\n    $freadme =  \"sfold.READM\
E\";\n    $fstk = \"sfold.stk\";\n    \n    # prep\
are execution...\n    #\n    # ./tmp is essential \
for dswpalign\n    # otherwise you'll get a segfau\
lt\n    mkdir \"./tmp\";\n    \n    $cmd = \"sfold\
 -m $fmodel $fs1 $fs2\";\n    open(FREADME,\">$fre\
adme\");\n    print FREADME \"$cmd\\n\"; \n    clo\
se(FREADME);\n\n    # and go\n    #\n    system(\"\
/usr/bin/time -p -o $ftimelog $cmd >$frawout 2>$fe\
rrlog\");\n    if ($? != 0) {\n        $errmsg = \\
"command failed with exit status $?\";\n        $e\
rrmsg .= \"command was \\\"$cmd\\\". See \" . getc\
wd . \"\\n\";\n        return -1;\n    }\n\n    re\
turn sfold_parseoutput($frawout, $foutfa);\n}\n\n\\
n\n\n\n\n\nmy ($help, $fmsqin, $fmsaout);\nGetOpti\
ons(\"help\"  => \\$help,\n           \"in=s\" => \
\\$fmsqin,\n           \"out=s\" => \\$fmsaout,\n	\
   \"data=s\" => \\$ref_dir);\n\n\n\nif ($help) {\\
n    usage();\n    exit(0);\n}\nif (! defined($fms\
qin)) {\n    usage('missing input filename');\n   \
 exit(1);\n}\nif (! defined($fmsaout)) {\n    usag\
e('missing output filename');\n    exit(1);\n\n}\n\
if (scalar(@ARGV)) {\n    usage('Unknown remaining\
 args');\n    exit(1);\n}\n\n$FMODEL = \"$ref_dir/\
mix80.mod\";\nif (! -e \"$FMODEL\") {\n    die(\"c\
ouldn't find sfold grammar model file. Expected $F\
MODEL\\n\");\n}\n\n\nmy %hseq=read_fasta_seq ($fms\
qin);\nmy $id;\n\nforeach $id (keys(%hseq))\n  {\n\
    push(@seq_array, $hseq{$id});\n  }\n\nif ( sca\
lar(@seq_array) != 2 ) {\n    die(\"Need *exactly*\
 two sequences as input (pairwise alignment!).\")\\
n}\n\n\n\nmy ($sec, $min, $hour, $mday, $mon, $yea\
r, $wday, $yday, $isdst) = localtime(time);\nmy $d\
atei = sprintf(\"%4d-%02d-%02d\", $year+1900, $mon\
+1, $mday);\nmy $templ = basename($0) . \".\" . $d\
atei . \".pid-\" . $$ . \".XXXXXX\";\nmy $wd = tem\
pdir ( $templ, DIR => $TMPDIR);\n\ncopy($fmsqin, \\
"$wd/\" . basename($fmsqin) . \".org\"); # for rep\
roduction\ncopy($FMODEL, \"$wd\");\nmy $fmodel = b\
asename($FMODEL);\nmy $orgwd = getcwd;\nchdir $wd;\
\n\n\n\nmy @sepseqfiles;\nforeach $id (keys(%hseq)\
) {\n    my ($seq, $orgseq, $fname, $sout);\n    $\
seq=$hseq{$id}{'seq'};\n    \n    $fname = basenam\
e($fmsqin) . \"_$id.fa\";\n    # replace funnies i\
n file/id name (e.g. \"/\" \" \" etc)\n    $fname \
=~ s,[/ ],_,g;\n    open (PF, \">$fname\");\n    p\
rint (PF \">$id\\n$seq\\n\");\n    close (PF);\n\n\
    push(@sepseqfiles, $fname);\n}\n\nmy ($f1, $f2\
, $fout);\n$f1 = $sepseqfiles[0];\n$f2 = $sepseqfi\
les[1];\n$fout = $wd . basename($fmsqin) . \".out.\
fa\";\nif (sfold_wrapper($f1, $f2, $fmodel, \"$fou\
t\") != 0) {\n    printf STDERR \"ERROR: See logs \
in $wd\\n\";\n    exit(1);\n} else {\n    chdir $o\
rgwd;\n    copy($fout, $fmsaout);\n    rmtree($wd)\
;\n   exit(0);\n}\n","\nuse Env qw(HOST);\nuse Env\
 qw(HOME);\nuse Env qw(USER);\n\n\n$tmp=clean_cr (\
$ARGV[0]);\nopen (F, $tmp);\n\nwhile ( <F>)\n  {\n\
    my $l=$_;\n    if ( $l=~/^# STOCKHOLM/){$stock\
holm=1;}\n    elsif ( $stockholm && $l=~/^#/)\n   \
   {\n	$l=~/^#(\\S+)\\s+(\\S+)\\s+(\\S*)/g;\n	$l=\\
"_stockholmhasch_$1\\_stockholmspace_$2 $3\\n\";\n\
      }\n    $file.=$l;\n  }\nclose (F);\nunlink($\
tmp);\n$file1=$file;\n\n$file=~s/\\#/_hash_symbol_\
/g;\n$file=~s/\\@/_arobase_symbol_/g;\n\n\n$file=~\
s/\\n[\\.:*\\s]+\\n/\\n\\n/g;\n\n$file=~s/\\n[ \\t\
\\r\\f]+(\\b)/\\n\\1/g;\n\n\n$file=~s/(\\n\\S+)(\\\
s+)(\\S)/\\1_blank_\\3/g;\n\n$file=~s/[ ]//g;\n$fi\
le=~s/_blank_/ /g;\n\n\n\n$file =~s/\\n\\s*\\n/#/g\
;\n\n$file.=\"#\";\n$file =~s/\\n/@/g;\n\n\n\n\n@b\
locks=split /\\#/, $file;\nshift (@blocks);\n@s=sp\
lit /\\@/, $blocks[0];\n$nseq=$#s+1;\n\n\n\n$file=\
join '@', @blocks;\n@lines=split /\\@/,$file;\n\n$\
c=0;\n\nforeach $l (@lines)\n  {\n    if (!($l=~/\\
\S/)){next;}\n    elsif ($stockholm && ($l=~/^\\/\\
\// || $l=~/STOCKHOLM/)){next;}#get read of STOCHO\
LM Terminator\n   \n    $l=~/(\\S+)\\s+(\\S*)/g;\n\
    $n=$1; $s=$2;\n    \n    $seq[$c].=$s;\n    $n\
ame[$c]=$n;\n    $c++;\n    \n    if ( $c==$nseq){\
$c=0;}\n    \n  } \n\nif ( $c!=0)\n      {\n	print\
 STDERR \"ERROR: $ARGV[0] is NOT an MSA in Clustal\
w format: make sure there is no blank line within \
a block [ERROR]\\n\";\n	exit (EXIT_FAILURE);\n    \
  }\n\nfor ($a=0; $a< $nseq; $a++)\n  {\n    $name\
[$a]=cleanstring ($name[$a]);\n    $seq[$a]=cleans\
tring ($seq[$a]);\n    $seq[$a]=breakstring($seq[$\
a], 60);\n    \n    $line=\">$name[$a]\\n$seq[$a]\\
\n\";\n    \n    print \"$line\";\n  }\nexit (EXIT\
_SUCCESS);\n\nsub cleanstring\n  {\n    my $s=@_[0\
];\n    $s=~s/_hash_symbol_/\\#/g;\n    $s=~s/_aro\
base_symbol_/\\@/g;\n    $s=~s/[ \\t]//g;\n    ret\
urn $s;\n  }\nsub breakstring\n  {\n    my $s=@_[0\
];\n    my $size=@_[1];\n    my @list;\n    my $n,\
$ns, $symbol;\n    \n    @list=split //,$s;\n    $\
n=0;$ns=\"\";\n    foreach $symbol (@list)\n      \
{\n	if ( $n==$size)\n	  {\n	    $ns.=\"\\n\";\n	  \
  $n=0;\n	  }\n	$ns.=$symbol;\n	$n++;\n      }\n  \
  return $ns;\n    }\n\nsub clean_cr\n  {\n    my \
$f=@_[0];\n    my $file;\n    \n    $tmp=\"f$.$$\"\
;\n    \n    \n    open (IN, $f);\n    open (OUT, \
\">$tmp\");\n    \n    while ( <IN>)\n      {\n	$f\
ile=$_;\n	$file=~s/\\r\\n/\\n/g;\n	$file=~s/\\n\\r\
/\\n/g;\n	$file=~s/\\r\\r/\\n/g;\n	$file=~s/\\r/\\\
n/g;\n	print OUT \"$file\";\n      }\n    \n    cl\
ose (IN);\n    close (OUT);\n    return $tmp;\n  }\
\n","use Env qw(HOST);\nuse Env qw(HOME);\nuse Env\
 qw(USER);\n\n\n$query_start=-1;\n$query_end=-1;\n\
\nwhile (<>)\n  {\n    if ( /\\/\\//){$in_aln=1;}\\
n    elsif ( $in_aln && /(\\S+)\\s+(.*)/)\n      {\
\n\n\n	$name=$1;\n	\n\n	$seq=$2;\n	$seq=~s/\\s//g;\
\n        $seq=~s/\\~/\\-/g;\n	$seq=~s/\\./\\-/g;\\
n	if ( $list{$n}{'name'} && $list{$n}{'name'} ne $\
name)\n	  {\n	    print \"$list{$n}{'name'} Vs $na\
me\";\n	    \n	    exit (EXIT_FAILURE);\n	  }\n	el\
se\n	  {\n	    $list{$n}{'name'}= $name;\n	  }\n\n\
	$list{$n}{'seq'}=$list{$n}{'seq'}.$seq;\n	\n	$nse\
q=++$n;\n	\n      }\n    else\n      {$n=0;}\n  }\\
n\n\nfor ($a=0; $a<$nseq; $a++)\n  {\n    print \"\
>$list{$a}{'name'}\\n$list{$a}{'seq'}\\n\";\n  }\n\
      \n","\nuse Env qw(HOST);\nuse Env qw(HOME);\\
nuse Env qw(USER);\n\n                            \
                            \nuse strict;         \
                                    \nuse warnings\
;\nuse diagnostics;\n\nmy $in_hit_list, my $in_aln\
=0, my(%name_list)=(),my (%list)=(),my $n_seq=0; m\
y $test=0;\nmy($j)=0, my $n=0, my $nom, my $lg_que\
ry, my %vu=();\n\nopen (F, \">tmp\");\n\n$/=\"\\n\\
";\nwhile (<>)\n{\n    print F $_;\n    if($_ =~ /\
Query=\\s*(.+?)\\s/i) { $nom=$1;}\n\n    if ( /Seq\
uences producing significant alignments/){$in_hit_\
list=1;}\n    \n    if ($_=~ /^pdb\\|/i) { $_=~ s/\
pdb\\|//g; }\n    if ($_=~ /^(1_\\d+)\\s+\\d+/) { \
$_=~ s/$1/QUERY/;}\n      \n    if ( /^(\\S+).+?\\\
s+[\\d.]+\\s+([\\de.-]+)\\s+$/ && $in_hit_list)	\n\
    {\n	my($id)=$1; # \n	$id=~ s/\\|/_/g; #\n	if (\
$id =~ /.+_$/) { chop($id) }; #\n	$name_list{$n_se\
q++}=$id;\n	$name_list{$n_seq-1}=~ s/.*\\|//g;    \
 \n    }\n  \n    if (/query/i) {$in_aln=1;}\n    \
if ( /^(\\S+)\\s+(\\d+)\\s+([a-zA-Z-]+)\\s+(\\d+)/\
 || /^(\\S+)(\\s+)(\\-+)(\\s+)/ && ($in_aln == 1))\
\n    {\n	my $name=$1;\n	my $start=$2;\n	my $seq=$\
3;\n	my $end=$4;\n		\n	if ($name =~ /QUERY/i) { $l\
g_query=length($seq); }\n\n	unless ($test > $n) #m\
\n	{\n	    my(@seqq)= split('',$seq);\n	    my($ga\
p_missing)= scalar(@seqq);\n	    \n	    while ($ga\
p_missing != $lg_query)  { unshift (@seqq,\"-\"); \
$gap_missing= scalar(@seqq); }\n	    $seq=join('',\
@seqq);  #m\n	}\n	\n	if ($name =~ /QUERY/i)\n	{\n	\
    $n=0; %vu=(); $j=0;\n	    $list{$n}{'real_name\
'}=\"$nom\";\n	}	\n	else\n	{\n	    unless (exists \
$vu{$name}) { ++$j;}	\n	    $list{$n}{'real_name'}\
=$name_list{$j-1};\n	}\n		\n	$list{$n}{'name'}=$na\
me;\n\n	$seq=~tr/a-z/A-Z/;\n	$list{$n}{'seq'}=$lis\
t{$n}{'seq'};\n	$list{$n}{'seq'}.=$seq;\n\n	$n++;\\
n	$vu{$name}++;\n	$test++;\n   } \n    \n}\n\nmy @\
numero=();\n\nfor (my $a=0; $a<$n; $a++) #m\n{\n  \
  my $long=length($list{0}{'seq'});  \n    my $lon\
g1= length($list{$a}{'seq'});\n  \n    while ($lon\
g1 ne $long)\n    {\n	$list{$a}{'seq'}.=\"-\";\n	$\
long1= length ($list{$a}{'seq'});\n    } \n \n    \
push (@numero,\"$list{$a}{'name'} $list{$a}{'real_\
name'}\\n\");\n}\n\nmy %dejavu=();\n\n\nfor (my $i\
=0; $i<=$#numero; $i++)\n{\n    my $s=\">$list{$i}\
{'real_name'}\\n$list{$i}{'seq'}\\n\";\n    my $k=\
0;\n    \n    if (exists $dejavu{$numero[$i]}) {ne\
xt;}\n    else\n    {	\n	for ($j=0; $j<$n ; $j++)\\
n	{\n	    if (\"$numero[$i]\" eq \"$numero[$j]\" &\
& $j != $i )\n	    {\n		++$k;\n		$s .=\">$list{$j}\
{'real_name'}\\n$list{$j}{'seq'}\\n\";\n	    }\n	}\
	\n    }\n    \n    if ($k>0) \n    {\n	my $cons;\\
n	open (SOR,\">tempo_aln2cons\"); print SOR $s;  c\
lose SOR ;\n	open (COM,\"t_coffee -other_pg seq_re\
format -in tempo_aln2cons -action +aln2cons +upper\
 |\") ; \n     	while (<COM>)\n	{	\n	    if (/^>/)\
 { $cons =\">$list{$i}{'real_name'}\\n\"; next;}\n\
	    $_=~ s/\\n//g;\n	    $cons .=$_;\n	}\n	close \
COM; unlink (\"tempo_aln2cons\");\n	print $cons,\"\
\\n\"; print F $cons,\"\\n\";\n    }	\n    else  {\
 print $s;  print F $s; }\n    \n    $dejavu{$nume\
ro[$i]}++;\n} #m\n\nexit;\n\n\n\n\n\n\n\n\n\n\n\n"\
,"use Env;\n\n\n$tmp_dir=\"\";\n$init_dir=\"\";\n$\
program=\"tc_generic_method.pl\";\n\n$blast=@ARGV[\
0];\n\n$name=\"query\";$seq=\"\";\n%p=blast_xml2pr\
ofile($name,$seq,100, 0, 0, $blast);\n&output_prof\
ile (%p);\n\n\nsub output_profile\n  {\n    my (%p\
rofile)=(@_);\n    my ($a);\n    for ($a=0; $a<$pr\
ofile{n}; $a++)\n      {\n	\n	print \">$profile{$a\
}{name} $profile{$a}{comment}\\n$profile{$a}{seq}\\
\n\";\n      }\n    return;\n  }\nsub file_contain\
s \n  {\n    my ($file, $tag, $max)=(@_);\n    my \
($n);\n    $n=0;\n    \n    if ( !-e $file && ($fi\
le =~/$tag/)) {return 1;}\n    elsif ( !-e $file){\
return 0;}\n    else \n      {\n	open (FC, \"$file\
\");\n	while ( <FC>)\n	  {\n	    if ( ($_=~/$tag/)\
)\n	      {\n		close (FC);\n		return 1;\n	      }\\
n	    elsif ($max && $n>$max)\n	      {\n		close (\
FC);\n		return 0;\n	      }\n	    $n++;\n	  }\n   \
   }\n    close (FC);\n    return 0;\n  }\n	    \n\
	  \nsub file2string\n  {\n    my $f=@_[0];\n    m\
y $string, $l;\n    open (F,\"$f\");\n    while (<\
F>)\n      {\n\n	$l=$_;\n	#chomp ($l);\n	$string.=\
$l;\n      }\n    close (F);\n    $string=~s/\\r\\\
n//g;\n    $string=~s/\\n//g;\n    return $string;\
\n  }\n\n\n\nsub tag2value \n  {\n    \n    my $ta\
g=(@_[0]);\n    my $word=(@_[1]);\n    my $return;\
\n    \n    $tag=~/$word=\"([^\"]+)\"/;\n    $retu\
rn=$1;\n    return $return;\n  }\n      \nsub hit_\
tag2pdbid\n  {\n    my $tag=(@_[0]);\n    my $pdbi\
d;\n       \n    $tag=~/id=\"(\\S+)\"/;\n    $pdbi\
d=$1;\n    $pdbid=~s/_//;\n    return $pdbid;\n  }\
\nsub id2pdbid \n  {\n    my $id=@_[0];\n  \n    i\
f ($id =~/pdb/)\n      {\n	$id=~/pdb(.*)/;\n	$id=$\
1;\n      }\n    $id=~s/[|�_]//g;\n    return $id;\
\n  }\nsub set_blast_type \n  {\n    my $file =@_[\
0];\n    if (&file_contains ($file,\"EBIApplicatio\
nResult\",100)){$BLAST_TYPE=\"EBI\";}\n    elsif (\
&file_contains ($file,\"NCBI_BlastOutput\",100)) {\
$BLAST_TYPE=\"NCBI\";}\n    else\n      {\n	$BLAST\
_TYPE=\"\";\n      }\n    return $BLAST_TYPE;\n  }\
\nsub blast_xml2profile \n  {\n    my ($name,$seq,\
$maxid, $minid, $mincov, $file)=(@_);\n    my (%p,\
 $a, $string, $n);\n    \n\n\n    if ($BLAST_TYPE \
eq \"EBI\" || &file_contains ($file,\"EBIApplicati\
onResult\",100)){%p=ebi_blast_xml2profile(@_);}\n \
   elsif ($BLAST_TYPE eq \"NCBI\" || &file_contain\
s ($file,\"NCBI_BlastOutput\",100)){%p=ncbi_blast_\
xml2profile(@_);}\n    else \n      {\n	print \"**\
********** ERROR: Blast Returned an unknown XML Fo\
rmat **********************\";\n	die;\n      }\n  \
  for ($a=0; $a<$p{n}; $a++)\n      {\n	my $name=$\
p{$a}{name};\n	$p{$name}{seq}=$p{$a}{seq};\n      \
}\n    return %p;\n  }\nsub ncbi_blast_xml2profile\
 \n  {\n    my ($name,$seq,$maxid, $minid, $mincov\
, $string)=(@_);\n    my ($L,$l, $a,$b,$c,$d,$nhit\
s,@identifyerL);\n    \n    \n    $seq=~s/[^a-zA-Z\
]//g;\n    $L=length ($seq);\n    \n    %hit=&xml2\
tag_list ($string, \"Hit\");\n    \n    \n    for \
($nhits=0,$a=0; $a<$hit{n}; $a++)\n      {\n	my ($\
ldb,$id, $identity, $expectation, $start, $end, $c\
overage, $r);\n	my (%ID,%DE,%HSP);\n	\n	$ldb=\"\";\
\n\n	%ID=&xml2tag_list ($hit{$a}{body}, \"Hit_id\"\
);\n	$identifyer=$ID{0}{body};\n	\n	%DE=&xml2tag_l\
ist ($hit{$a}{body}, \"Hit_def\");\n	$definition=$\
DE{0}{body};\n	\n	%HSP=&xml2tag_list ($hit{$a}{bod\
y}, \"Hsp\");\n	for ($b=0; $b<$HSP{n}; $b++)\n	  {\
\n	    my (%START,%END,%E,%I,%Q,%M);\n\n	 \n	    %\
START=&xml2tag_list ($HSP{$b}{body}, \"Hsp_query-f\
rom\");\n	    %HSTART=&xml2tag_list ($HSP{$b}{body\
}, \"Hsp_hit-from\");\n	    \n	    %LEN=  &xml2tag\
_list ($HSP{$b}{body}, \"Hsp_align-len\");\n	    %\
END=  &xml2tag_list ($HSP{$b}{body}, \"Hsp_query-t\
o\");\n	    %HEND=  &xml2tag_list ($HSP{$b}{body},\
 \"Hsp_hit-to\");\n	    %E=&xml2tag_list     ($HSP\
{$b}{body}, \"Hsp_evalue\");\n	    %I=&xml2tag_lis\
t     ($HSP{$b}{body}, \"Hsp_identity\");\n	    %Q\
=&xml2tag_list     ($HSP{$b}{body}, \"Hsp_qseq\");\
\n	    %M=&xml2tag_list     ($HSP{$b}{body}, \"Hsp\
_hseq\");\n	    \n	    for ($e=0; $e<$Q{n}; $e++)\\
n\n	      {\n		$qs=$Q{$e}{body};\n		$ms=$M{$e}{bod\
y};\n		if ($seq eq\"\"){$seq=$qs;$L=length($seq);}\
\n		\n		$expectation=$E{$e}{body};\n		$identity=($\
LEN{$e}{body}==0)?0:$I{$e}{body}/$LEN{$e}{body}*10\
0;\n		$start=$START{$e}{body};\n		$end=$END{$e}{bo\
dy};\n		$Hstart=$HSTART{$e}{body};\n		$Hend=$HEND{\
$e}{body};\n	\n		$coverage=(($end-$start)*100)/$L;\
\n\n	\n		if ($identity>$maxid || $identity<$minid \
|| $coverage<$mincov){next;}\n		@lr1=(split (//,$q\
s));\n		@lr2=(split (//,$ms));\n		$l=$#lr1+1;\n		f\
or ($c=0;$c<$L;$c++){$p[$nhits][$c]=\"-\";}\n		for\
 ($d=0,$c=0; $c<$l; $c++)\n		  {\n		    $r=$lr1[$c\
];\n		    if ( $r=~/[A-Za-z]/)\n		      {\n			\n		\
	$p[$nhits][$d + $start-1]=$lr2[$c];\n			$d++;\n		\
      }\n		  }\n		$Qseq[$nhits]=$qs;\n		$Hseq[$nhi\
ts]=$ms;\n		$QstartL[$nhits]=$start;\n		$HstartL[$\
nhits]=$Hstart;\n		$identityL[$nhits]=$identity;\n\
		$endL[$nhits]=$end;\n		$definitionL[$nhits]=$def\
inition;\n		$identifyerL[$nhits]=$identifyer;\n		$\
comment[$nhits]=\"$ldb|$identifyer [Eval=$expectat\
ion][id=$identity%][start=$Hstart end=$Hend]\";\n	\
	$nhits++;\n	      }\n	  }\n      }\n    \n    $pr\
ofile{n}=0;\n    $profile{$profile{n}}{name}=$name\
;\n    $profile{$profile{n}}{seq}=$seq;\n    $prof\
ile {n}++;\n    \n    for ($a=0; $a<$nhits; $a++)\\
n      {\n	$n=$a+1;\n	\n	$profile{$n}{name}=\"$nam\
e\\_$a\";\n	$profile{$n}{seq}=\"\";\n	$profile{$n}\
{Qseq}=$Qseq[$a];\n	$profile{$n}{Hseq}=$Hseq[$a];\\
n	$profile{$n}{Qstart}=$QstartL[$a];\n	$profile{$n\
}{Hstart}=$HstartL[$a];\n	$profile{$n}{identity}=$\
identityL[$a];\n	$profile{$n}{definition}=$definit\
ionL[$a];\n	$profile{$n}{identifyer}=$identifyerL[\
$a];\n	$profile{$n}{comment}=$comment[$a];\n	for (\
$b=0; $b<$L; $b++)\n	  {\n	    if ($p[$a][$b])\n	 \
     {\n		$profile{$n}{seq}.=$p[$a][$b];\n	      }\
\n	    else\n	      {\n		$profile{$n}{seq}.=\"-\";\
\n	      }\n	  }\n      }\n    \n    $profile{n}=$\
nhits+1;\n    return %profile;\n  }\nsub ebi_blast\
_xml2profile \n  {\n    my ($name,$seq,$maxid, $mi\
nid, $mincov, $string)=(@_);\n    my ($L,$l, $a,$b\
,$c,$d,$nhits,@identifyerL,$identifyer);\n    \n\n\
    \n    $seq=~s/[^a-zA-Z]//g;\n    $L=length ($s\
eq);\n    %hit=&xml2tag_list ($string, \"hit\");\n\
    \n    for ($nhits=0,$a=0; $a<$hit{n}; $a++)\n \
     {\n	my ($ldb,$id, $identity, $expectation, $s\
tart, $end, $coverage, $r);\n	my (%Q,%M,%E,%I);\n	\
\n	$ldb=&tag2value ($hit{$a}{open}, \"database\");\
\n	$identifyer=&tag2value ($hit{$a}{open}, \"id\")\
;\n\n	$description=&tag2value ($hit{$a}{open}, \"d\
escription\");\n	\n	%Q=&xml2tag_list ($hit{$a}{bod\
y}, \"querySeq\");\n	%M=&xml2tag_list ($hit{$a}{bo\
dy}, \"matchSeq\");\n	%E=&xml2tag_list ($hit{$a}{b\
ody}, \"expectation\");\n	%I=&xml2tag_list ($hit{$\
a}{body}, \"identity\");\n	\n\n	for ($b=0; $b<$Q{n\
}; $b++)\n	  {\n	    \n	    \n	    $qs=$Q{$b}{body\
};\n	    $ms=$M{$b}{body};\n	    if ($seq eq\"\"){\
$seq=$qs;$L=length($seq);}\n\n	    $expectation=$E\
{$b}{body};\n	    $identity=$I{$b}{body};\n	    \n\
	    	    \n	    $start=&tag2value ($Q{$b}{open}, \
\"start\");\n	    $end=&tag2value ($Q{$b}{open}, \\
"end\");\n	    $startM=&tag2value ($M{$b}{open}, \\
"start\");\n	    $endM=&tag2value ($M{$b}{open}, \\
"end\");\n	    $coverage=(($end-$start)*100)/$L;\n\
	    \n	   # print \"$id: ID: $identity COV: $cove\
rage [$start $end]\\n\";\n	    \n	    \n	    if ($\
identity>$maxid || $identity<$minid || $coverage<$\
mincov){next;}\n	    # print \"KEEP\\n\";\n\n	    \
\n	    @lr1=(split (//,$qs));\n	    @lr2=(split (/\
/,$ms));\n	    $l=$#lr1+1;\n	    for ($c=0;$c<$L;$\
c++){$p[$nhits][$c]=\"-\";}\n	    for ($d=0,$c=0; \
$c<$l; $c++)\n	      {\n		$r=$lr1[$c];\n		if ( $r=\
~/[A-Za-z]/)\n		  {\n		    \n		    $p[$nhits][$d +\
 $start-1]=$lr2[$c];\n		    $d++;\n		  }\n	      }\
\n	  \n	    \n	    $identifyerL[$nhits]=$identifye\
r;\n	    $comment[$nhits]=\"$ldb|$identifyer [Eval\
=$expectation][id=$identity%][start=$startM end=$e\
ndM]\";\n	    $nhits++;\n	  }\n      }\n    \n    \
$profile{n}=0;\n    $profile{$profile{n}}{name}=$n\
ame;\n    $profile{$profile{n}}{seq}=$seq;\n    $p\
rofile {n}++;\n    \n    for ($a=0; $a<$nhits; $a+\
+)\n      {\n	$n=$a+1;\n	$profile{$n}{name}=\"$nam\
e\\_$a\";\n	$profile{$n}{seq}=\"\";\n	$profile{$n}\
{identifyer}=$identifyerL[$a];\n	\n	$profile{$n}{c\
omment}=$comment[$a];\n	for ($b=0; $b<$L; $b++)\n	\
  {\n	    if ($p[$a][$b])\n	      {\n		$profile{$n\
}{seq}.=$p[$a][$b];\n	      }\n	    else\n	      {\
\n		$profile{$n}{seq}.=\"-\";\n	      }\n	  }\n   \
   }\n    $profile{n}=$nhits+1;\n    \n    return \
%profile;\n  }\n\nsub blast_xml2hit_list\n  {\n   \
 my $string=(@_[0]);\n    return &xml2tag_list ($s\
tring, \"hit\");\n  }\nsub xml2tag_list  \n  {\n  \
  my ($string_in,$tag)=@_;\n    my $tag_in, $tag_o\
ut;\n    my %tag;\n    \n    if (-e $string_in)\n \
     {\n	$string=&file2string ($string_in);\n     \
 }\n    else\n      {\n	$string=$string_in;\n     \
 }\n    $tag_in1=\"<$tag \";\n    $tag_in2=\"<$tag\
>\";\n    $tag_out=\"/$tag>\";\n    $string=~s/>/>\
##1/g;\n    $string=~s/</##2</g;\n    $string=~s/#\
#1/<#/g;\n    $string=~s/##2/#>/g;\n    @l=($strin\
g=~/(\\<[^>]+\\>)/g);\n    $tag{n}=0;\n    $in=0;$\
n=-1;\n  \n \n\n    foreach $t (@l)\n      {\n\n	$\
t=~s/<#//;\n	$t=~s/#>//;\n	\n	if ( $t=~/$tag_in1/ \
|| $t=~/$tag_in2/)\n	  {\n	 \n	    $in=1;\n	    $t\
ag{$tag{n}}{open}=$t;\n	    $n++;\n	    \n	  }\n	e\
lsif ($t=~/$tag_out/)\n	  {\n	    \n\n	    $tag{$t\
ag{n}}{close}=$t;\n	    $tag{n}++;\n	    $in=0;\n	\
  }\n	elsif ($in)\n	  {\n	   \n	    $tag{$tag{n}}{\
body}.=$t;\n	  }\n      }\n  \n    return %tag;\n \
 }\n\n\n\n\n","use Env qw(HOST);\nuse Env qw(HOME)\
;\nuse Env qw(USER);\nwhile (<>)\n  {\n    if ( /^\
>(\\S+)/)\n      {\n	if ($list{$1})\n	  {\n	    pr\
int \">$1_$list{$1}\\n\";\n	    $list{$1}++;\n	  }\
\n	else\n	  {\n	    print $_;\n	    $list{$1}=1;\n\
	  }\n      }\n    else\n      {\n	print $_;\n    \
  }\n  }\n      \n","\n\n\nuse Env qw(HOST);\nuse \
Env qw(HOME);\nuse Env qw(USER);\n\n\nopen (F,$ARG\
V[0]);\nwhile ( <>)\n  {\n    @x=/([^:,;\\)\\(\\s]\
+):[^:,;\\)\\(]*/g;\n    @list=(@list,@x);\n  }\n$\
n=$#list+1;\nforeach $n(@list){print \">$n\\nseque\
nce\\n\";}\n\n\nclose (F);\n","\nopen (F, $ARGV[0]\
);\n\nwhile ( <F>)\n  {\n    @l=($_=~/(\\S+)/g);\n\
    \n    $name=shift @l;\n    \n    print STDOUT \
\"\\n>$name\\n\";\n    foreach $e (@l){$e=($e eq \\
"0\")?\"O\":\"I\";print \"$e\";}\n  }\nclose (F);\\
n\n		       \n    \n","use Env qw(HOST);\nuse Env \
qw(HOME);\nuse Env qw(USER);\n\n$tmp=\"$ARGV[0].$$\
\";\nopen (IN, $ARGV[0]);\nopen (OUT, \">$tmp\");\\
n\nwhile ( <IN>)\n  {\n    $file=$_;\n    $file=~s\
/\\r\\n/\\n/g;\n    $file=~s/\\n\\r/\\n/g;\n    $f\
ile=~s/\\r\\r/\\n/g;\n    $file=~s/\\r/\\n/g;\n   \
 print OUT \"$file\";\n  }\nclose (IN);\nclose (OU\
T);\n\nopen (OUT, \">$ARGV[0]\");\nopen (IN, \"$tm\
p\");\n\nwhile ( <IN>)\n{\n  print OUT \"$_\";\n}\\
nclose (IN);\nclose (OUT);\nunlink ($tmp);\n\n"};
