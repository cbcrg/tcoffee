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
",$outfile);\n}\nsub blast_msa\n  {\n    my ($infi\
le,$outfile)=@_;\n    my ($a, %seq);\n    %s1=&rea\
d_fasta_seq ($infile);\n    foreach $s (keys (%s1)\
)\n      {\n	$i=$s1{$s}{order};\n	$s{$i}{name}=$s;\
\n	$s{$i}{seq}=$s1{$s}{seq};\n	$s{$i}{len}=length(\
 $s{$i}{seq});\n	$s{n}++;\n      }\n    `formatdb \
-i $infile`;\n    `blastpgp -i $infile -d $infile \
-m7 -j4 > io`;\n    &set_blast_type (\"io\");\n   \
 \n    %FB=&xml2tag_list (\"io\", \"BlastOutput\")\
;\n    \n    open (F, \">$outfile\");\n    print F\
 \"! TC_LIB_FORMAT_01\\n\";\n    print F \"$s{n}\\\
n\";\n    for ( $a=0; $a<$s{n}; $a++)\n      {\n	p\
rint F \"$s{$a}{name} $s{$a}{len} $s{$a}{seq}\\n\"\
;\n      }\n    for ( $a=0; $a<$FB{n}; $a++)\n    \
  {\n	%p=blast_xml2profile ($s{$a}{name}, $s{$a}{s\
eq},100, 0, 0, $FB{$a}{body});\n	for ($b=1; $b<$p{\
n}; $b++)\n	  {\n	    my $l=length ($p{$b}{Qseq});\
\n	    my $hit=$p{$b}{definition};\n	    my $Qstar\
t=$p{$b}{Qstart};\n	    my $Hstart=$p{$b}{Hstart};\
\n	    my $identity=$p{$b}{identity};\n	    my @lr\
Q=split (//,$p{$b}{Qseq});\n	    my @lrH=split (//\
,$p{$b}{Hseq});\n	    my $i= $s1{$s{$a}{name}}{ord\
er}+1;\n	    my $j= $s1{$hit}{order}+1;\n	    #if \
( $j==$i){next;}\n	    printf F \"# %d %d\\n\", $i\
, $j;\n	    #  print  F \"\\n$p{$b}{Qseq} ($Qstart\
)\\n$p{$b}{Hseq} ($Hstart)\";\n	    for ($c=0; $c<\
$l; $c++)\n	      {\n		my $rQ=$lrQ[$c];\n		my $rH=\
$lrH[$c];\n		my $n=0;\n		\n		if ($rQ ne \"-\"){$n+\
+, $Qstart++;}\n		if ($rH ne \"-\"){$n++; $Hstart+\
+;}\n		\n		if ( $n==2)\n		  {\n		    printf F \"\\\
t%d %d %d\\n\", $Qstart-1, $Hstart-1,$identity;\n	\
	  }\n	      }\n	  }\n      }\n    print F \"! SEQ\
_1_TO_N\\n\";\n    close (F);\n    return $output;\
\n  \n  }\n\nsub seq2msa\n  {\n    my ($mode, $inf\
ile, $method, $param, $outfile)=@_;\n    &set_temp\
orary_dir (\"set\",$infile,\"seq.pep\");\n    $par\
am.=\" >/dev/null 2>&1 \";\n    \n    #make sure t\
est.pep is in FASTA\n    &safe_system (\"t_coffee \
-other_pg seq_reformat -in seq.pep -output fasta_s\
eq > x\");\n    `mv x seq.pep`;\n    \n    if ( $m\
ethod eq \"blastpgp\")\n      {\n	&blast_msa (\"se\
q.pep\", \"result.aln\");\n      }\n    elsif ( $m\
ethod eq \"muscle\")\n      {\n	`muscle -in seq.pe\
p -out result.aln $param`;\n      }\n    elsif ( $\
method eq \"probcons\")\n      {\n	`probcons seq.p\
ep >result.aln 2>/dev/null`;\n      }\n    elsif (\
 $method eq \"mafft\")\n      {\n	`mafft --quiet -\
-localpair --maxiterate 1000 seq.pep> result.aln  \
2>/dev/null`\n      }\n    elsif ( $method=~/prank\
/)\n      {\n	`$method -d=seq.pep -o=result.aln -q\
uiet 2>/dev/null`;\n	`mv result.aln.1.fas result.a\
ln`;\n      }\n    else\n      {\n	`$method -infil\
e=seq.pep -outfile=result.aln`;\n      }\n    \n  \
  &set_temporary_dir (\"unset\",$mode, $method,\"r\
esult.aln\",$outfile);\n    myexit ($EXIT_SUCCESS)\
;\n  }\n\nsub seq2thread_pair\n  {\n    my ($mode,\
 $infile, $pdbfile1, $method, $param, $outfile)=@_\
;\n    &set_temporary_dir (\"set\",$infile,\"seq.p\
ep\",$pdbfile1,\"struc.pdb\");\n    if ($method eq\
 \"fugueali\")\n      {\n	#Env Variable that need \
to be defined for Fugue\n	if (!$ENV{FUGUE_LIB_LIST\
}){$ENV{FUGUE_LIB_LIST}=\"DUMMY\";}\n	if (!$ENV{HO\
MSTRAD_PATH})  {$ENV{HOMSTRAD_PATH}=\"DUMMY\";}\n	\
if (!$ENV{HOMS_PATH}){$ENV{HOMS_PATH}=\"DUMMY\";}\\
n	\n	`joy struc.pdb >x 2>x`;\n	&check_file(\"struc\
.tem\", \"Joy failed [FATAL:$PROGRAM/$method]\");\\
n	`melody -t struc.tem >x 2>x`;\n	&check_file(\"st\
ruc.tem\", \"Melody failed [FATAL:$PROGRAM/$method\
]\");\n	`fugueali -seq seq.pep -prf struc.fug -pri\
nt > tmp_result.aln`;\n	\n	&check_file(\"tmp_resul\
t.aln\", \"Fugue failed [FATAL:$PROGRAM/$method]\"\
);\n	&safe_system (\"t_coffee -other_pg seq_reform\
at -in tmp_result.aln -output fasta_aln >result.al\
n\");\n      }\n    elsif ( $method eq \"t_coffee\\
")\n      {\n	&safe_system (\"t_coffee -in Pstruc.\
pdb Sseq.pep Mslow_pair -outfile result.aln -quiet\
\");\n      }\n    else\n      {\n	&safe_system (\\
"$method -infile=seq.pep -pdbfile1=struc.pdb -outf\
ile=result.aln $param>x 2>x\");\n      }\n    &set\
_temporary_dir (\"unset\",$mode,$method,\"result.a\
ln\",$outfile);\n    myexit ($EXIT_SUCCESS);\n  }\\
nsub seq2pdbid_pair\n  {\n    my ($mode, $pdbfile1\
, $pdbfile2, $method, $param, $outfile)=@_;\n    m\
y ($name);\n\n    \n    &set_temporary_dir (\"set\\
");\n    $name=$pdbfile1.\" \".$pdbfile2;\n\n    i\
f (    &cache_file(\"GET\",\"\",\"$name\",\"$metho\
d\",\"dali\",$outfile,\"EBI\"))\n      {return $ou\
tfile;}\n    else\n      {\n	if ($method eq \"dali\
web\")\n	  {\n	    $pdbfile1=~/(....)(.)/;\n	    $\
id1=$1; $c1=$2;\n	    \n	    $pdbfile2=~/(....)(.)\
/;\n	    $id2=$1; $c2=$2;\n	    \n	    $command=\"\
t_coffee -other_pg dalilite.pl --pdb1 $id1 --chain\
id1 $c1 --pdb2 $id2 --chainid2 $c2 --email=$EMAIL \
 >dali_stderr 2>dali_stderr\";\n	    $dali=`$comma\
nd`;\n	    \n	    open (F, \"dali_stderr\");\n	   \
 while (<F>)\n	      {\n		if ( /JobId: dalilite-(\\
\S+)/)\n		{\n		  $jobid=$1;\n		}\n	      }\n	    c\
lose (F);\n	    unlink (\"dali_stderr\");\n	    \n\
	    $output1=\"dalilite-$jobid.txt\";\n	    if ( \
-e $output1)\n	      {\n		unlink ($output1);\n		&u\
rl2file (\"http://www.ebi.ac.uk/Tools/es/cgi-bin/j\
obresults.cgi/dalilite/dalilite-$jobid/aln.html\",\
 \"output2\");\n		\n		if ( -e \"output2\")\n		  {\\
n		    my ($seq1, $seq2);\n		    $seq1=$seq2=\"\";\
\n		    \n		    open (F, \"output2\");\n		    whil\
e (<F>)\n		      {\n			$l=$_;\n			if ( $l=~/Query\\
\s+(\\S+)/)\n			  {\n			    $seq1.=$1;\n			  }\n		\
	elsif ( $l=~/Sbjct\\s+(\\S+)/)\n			  {\n			    $s\
eq2.=$1;\n			  }\n		      }\n		    close (F);\n		 \
   unlink (\"output2\");\n		    if ($seq1 ne \"\" \
&& $seq2 ne \"\")\n		      {\n			$output3=\">$A\\n\
$seq1\\n>$B\\n$seq2\\n\";\n			$output3=~s/\\./-/g;\
\n			open (F, \">result.aln\");\n			print F \"$out\
put3\";\n			close (F);\n		      }\n		  }\n	      }\
\n	  }\n      }\n    &cache_file(\"SET\",\"\",\"$n\
ame\",\"$method\",\"dali\",\"result.aln\",\"EBI\")\
;\n    &set_temporary_dir (\"unset\",$mode, $metho\
d, \"result.aln\",$outfile);\n    myexit ($EXIT_SU\
CCESS);\n  }\nsub seq2pdb_pair\n  {\n    my ($mode\
, $pdbfile1, $pdbfile2, $method, $param, $outfile)\
=@_;\n    \n    &set_temporary_dir (\"set\",$pdbfi\
le1,\"pdb1.pdb\",$pdbfile2,\"pdb2.pdb\");\n    if \
($method eq \"t_coffee\")\n      {\n	&safe_system \
(\"t_coffee -in Ppdb1.pdb Ppdb2.pdb -quiet -outfil\
e=result.aln\");\n      }\n    elsif ( $method eq \
\"DaliLite\")\n      {\n	if ( &safe_system (\"Dali\
Lite -pairwise pdb1.pdb pdb2.pdb >tmp1\")==$EXIT_S\
UCCESS)\n	  {\n	     my ($seq1, $seq2);\n	     $se\
q1=$seq2=\"\";\n		    \n	     open (F, \"tmp1\");\\
n	     while (<F>)\n	       {\n		 $l=$_;\n		 if ( \
$l=~/Query\\s+(\\S+)/)\n		   {\n		     $seq1.=$1;\\
n		   }\n		 elsif ( $l=~/Sbjct\\s+(\\S+)/)\n		   {\
\n		     $seq2.=$1;\n		   }\n	       }\n	     clos\
e (F);\n	     unlink (\"tmp1\");\n	     if ($seq1 \
ne \"\" && $seq2 ne \"\")\n	       {\n		 my $outpu\
t3=\">$A\\n$seq1\\n>$B\\n$seq2\\n\";\n		 $output3=\
~s/\\./-/g;\n		 open (F, \">result.aln\");\n		 pri\
nt F \"$output3\";\n		 close (F);\n	       }\n	   \
}\n	else\n	  {\n	    print \"ERROR: DalLite failed\
 to align the considered structures[tc_generic_met\
hod.pl]\\n\";\n	  }    \n      }\n    elsif ( $met\
hod eq \"TMalign\")\n      {\n	if ( &safe_system (\
\"TMalign pdb1.pdb pdb2.pdb >tmp1\")==$EXIT_SUCCES\
S)\n	  {\n	    `tail -4 tmp1 > tmp2`;\n	    \n	   \
 open (F, \"tmp2\");\n	    while (<F>)\n	      {\n\
		unshift(@l, $_);\n	      }\n	    close (F);\n	  \
  open (F, \">result.aln\");\n	    $l[3]=~s/[^a-zA\
-Z0-9-]/\\-/g;\n	    $l[1]=~s/[^a-zA-Z0-9-]/\\-/g;\
\n	    print F \">$A\\n$l[3]\\n>$B\\n$l[1]\\n\";\n\
	    close (F);\n	  }\n	else\n	  {\n	    print \"E\
RROR: TMalign failed to align the considered struc\
tures[tc_generic_method.pl]\\n\";\n	    `rm result\
.aln >/dev/null 2>/dev/null`;\n	  }\n      }\n    \
elsif ( $method eq \"mustang\")\n      {\n	if ( &s\
afe_system (\"mustang -i pdb1.pdb pdb2.pdb -F fast\
a >/dev/null 2>/dev/null\")==$EXIT_SUCCESS)\n	  {\\
n	    `mv results.afasta result.aln`;\n	  }\n	else\
\n	  {\n	    print \"ERROR: mustang failed to alig\
n the considered structures[tc_generic_method.pl]\\
\n\";\n	    `rm result.aln >/dev/null 2>/dev/null`\
;\n	  }\n      }\n    else\n      {\n	if ( &safe_s\
ystem (\"$method -pdbfile1=pdb1.pep -pdbfile2=pdb2\
.pdb -outfile=result.aln $param>x 2>x\")==$EXIT_SU\
CCESS)\n	  {\n	    `mv results.afasta result.aln`;\
\n	  }\n	else\n	  {\n	    print \"ERROR: $method f\
ailed to align the considered structures[tc_generi\
c_method.pl]\\n\";\n	    `rm result.aln >/dev/null\
 2>/dev/null`;\n	  }\n      }\n    &set_temporary_\
dir (\"unset\",$mode, $method, \"result.aln\",$out\
file);\n    myexit ($EXIT_SUCCESS);\n  }\n\nsub se\
q2profile_pair\n  {\n    my ($mode, $profile1, $pr\
ofile2, $method, $param, $outfile)=@_;\n    \n    \
\n    if ($method eq \"clustalw\")\n      {\n	&set\
_temporary_dir (\"set\",$profile1,\"prf1.aln\",$pr\
ofile2,\"prf2.aln\");\n	`clustalw -profile1=prf1.a\
ln -profile2=prf2.aln -outfile=result.aln`;\n	&set\
_temporary_dir (\"unset\",$mode, $method, \"result\
.aln\",$outfile);\n      }\n    elsif ( $method eq\
 \"hhalign\")\n      {\n	hhalign ( $profile1,$prof\
ile2,$outfile,$param);\n      }\n    else\n      {\
\n	\n	`$method -profile1=prf1.aln -profile2=prf2.a\
ln -outfile=result.aln $param>x 2>x`;\n      }\n  \
 \n    myexit ($EXIT_SUCCESS);\n  }\n\nsub pg_is_i\
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
  else\n      {\n	return 1;\n      }\n  }\nsub set\
_temporary_dir\n  {\n    my @list=@_;\n    my $dir\
_mode, $a, $mode, $method;\n  \n    $dir_mode=shif\
t (@list);\n\n    \n    if ( $dir_mode eq \"set\")\
\n      {\n	$initial_dir=cwd();\n	if ( !$tmp_dir)\\
n	  {\n	    $rand=rand (100000);\n	    $tmp_dir=\"\
$TMPDIR/tmp4tcoffee_profile_pair_dir_$$\\_P_$rand\\
";\n	  }\n	if ( !-d $tmp_dir)\n	  {\n	    push (@T\
MPDIR_LIST, $tmp_dir);\n	    `mkdir $tmp_dir`;\n	 \
 }\n	\n	for ( $a=0; $a<=$#list; $a+=2)\n	      {\n\
		`cp $list[$a] $tmp_dir/$list[$a+1]`;\n	      }\n\
	chdir $tmp_dir;\n      }\n    elsif ( $dir_mode e\
q \"unset\")\n      {\n	$mode=shift (@list);\n	$me\
thod=shift (@list);\n	\n	if (!-e $list[0])\n	  {\n\
	   myexit(flush_error(\"Program $method failed to\
 produce $list[1]\" ));\n	    myexit ($EXIT_FAILUR\
E);\n	  }\n	else\n	  {\n	    chdir $initial_dir;\n\
	    # `t_coffee -other_pg seq_reformat -in $tmp_d\
ir/$list[0] -output fasta_aln -out $tmp_dir/result\
2.aln`;\n	    `cp $tmp_dir/$list[0] $tmp_dir/resul\
t2.aln`;\n	    if ( $list[1] eq \"stdout\")\n	    \
  {\n		open (F, \"$tmp_dir/result2.aln\");\n		whil\
e (<F>){print $_;}close(F);\n	      }\n	    else\n\
	      {\n		`mv $tmp_dir/result2.aln $list[1]`;\n	\
      }\n	    shift (@list); shift (@list);\n	    \
foreach $f (@list)\n	      {\n		`mv $tmp_dir/$f .`\
;\n	      }\n	  }\n      }\n  }\n\n\n\n\nsub my_ge\
t_opt\n  {\n    my @list=@_;\n    my $cl, $a, $arg\
v, @argl;\n    \n    @argl=();\n    $cl=shift @lis\
t;\n    for ( $a=0; $a<=$#list; $a+=3)\n      {\n	\
$option=$list[$a];\n	$optional=$list[$a+1];\n	$sta\
tus=$list[$a+2];\n	$argv=\"\";\n	if ($cl=~/$option\
(\\S+)/){$argv=$1;}\n	@argl=(@argl,$argv);\n	\n	\n\
	#$optional:0=>optional\n	#$optional:1=>must be se\
t\n	#$status: 0=>no requirement\n	#$status: 1=>mus\
t be an existing file\n	#$status: 2=>must be an in\
stalled package\n	\n\n	if ($optional==0){;}\n	elsi\
f ( $optional==1 && $argv eq \"\")\n	  {\n	    mye\
xit(flush_error( \"ERROR: Option $option must be s\
et\"));\n	    myexit ($EXIT_FAILURE);\n	  }\n	if (\
$status==0){;}\n	elsif ($status ==1 && $argv ne \"\
\" && !-e $argv)\n	  {\n	    myexit(flush_error( \\
"File $argv must exist\"));\n	    myexit ($EXIT_FA\
ILURE);\n	  }\n	elsif ( $status==2 && $argv ne \"\\
" && &check_pg_is_installed ($argv)==0)\n	  {\n	  \
  myexit(flush_error( \" $argv is not installed\")\
);\n	    myexit ($EXIT_FAILURE);\n	  }\n      }\n\\
n    return @argl;\n    }\n\nsub check_file \n  {\\
n    my ($file, $msg)=@_;\n\n    if ( !-e $file)\n\
      {\n	myexit(flush_error(\"$msg\"));\n      }\\
n    }\nsub hhalign\n  {\n    my ($aln1, $aln2, $o\
utfile, $param)=@_;\n    my $h1, $h2;\n    \n    $\
h{0}{index}=0;\n    $h{1}{index}=1;\n    \n    $h{\
0}{aln}=$aln1;\n    $h{1}{aln}=$aln2;\n\n   \n\n  \
  %{$h{0}}=aln2psi_profile (%{$h{0}});\n    %{$h{1\
}}=aln2psi_profile (%{$h{1}});\n\n    $param=~s/#S\
/ /g;\n    $param=~s/#M/\\-/g;\n    $param=~s/#E/\\
\=/g;\n    \n\n    \n    $command=\"hhalign -i $h{\
0}{a3m} -t $h{1}{a3m} -tc $outfile.tmp -rank 1 -ma\
pt 0 $param\";\n    `$command`;\n    \n  #  `hhali\
gn -i $h{0}{a3m} -t $h{1}{a3m} -tc $outfile.tmp -r\
ank 1 -mapt 0 -gapf 0.8 -gapg 0.8`;\n    \n\n    #\
 To run global use the following\n    \n    open (\
I, \"$outfile.tmp\");\n    open (O, \">$outfile\")\
;\n    $h{0}{cons}=s/\\./x/g;\n    $h{1}{cons}=s/\\
\./x/g;\n\n    print O \"! TC_LIB_FORMAT_01\\n2\\n\
$h{0}{name} $h{0}{len} $h{0}{seq}\\n$h{1}{name} $h\
{1}{len} $h{1}{seq}\\n#1 2\\n\";\n    \n    while \
(<I>)\n      {\n	if (/(\\d+)\\s+(\\d+)\\s+(\\d+)/)\
\n	  {\n	    print O \"\\t$h{0}{$1}\\t$h{1}{$2}\\t\
$3\\n\";\n	  }\n      }\n    print O \"! SEQ_1_TO_\
N\\n\";\n\n    close (O);\n    close (I);\n  }\n\n\
sub aln2psi_profile\n  {\n    my (%h)=@_;\n    my \
($aln,$i,$hv, $a, @c, $n);\n   \n    $i=$h{index};\
\n    $aln=$h{aln};\n\n    `cp $aln $$.hhh_aln`;\n\
    $command=\"t_coffee -other_pg seq_reformat -in\
 $aln -output hasch\";\n    $hv=`$command`;chomp (\
$hv);\n    \n    $h{a2m}=\"$tmp/$hv.tmp4hhpred.a2m\
\";\n    $h{a3m}=\"$tmp/$hv.tmp4hhpred.a3m\";\n   \
 if ( -e $h{a3m}){;}\n    else\n      {\n	`hhconse\
nsus  -M 50 -i $h{aln} -oa2m $h{a2m}`;\n	if (!-e $\
h{a2m})\n	  {\n	    print STDERR \"Program tc_gene\
ric_method.pl FAILED to run:\\n\\thhconsensus  -M \
50 -i $h{aln} -oa2m $h{a2m}\";\n	    myexit ($EXIT\
_FAILURE);\n	  }\n	\n	`hhconsensus  -M 50 -i $h{al\
n} -oa3m $h{a3m}`;\n	if (!-e $h{a3m})\n	  {\n	    \
print STDERR \"Program tc_generic_method.pl FAILED\
 to run:\\n\\thhconsensus  -M 50 -i $h{aln} -oa3m \
$h{a3m}\";\n	    myexit ($EXIT_FAILURE);\n	  }\n  \
     `buildali.pl $h{a3m} -n 1`;\n      }\n    \n \
   \n    $h{a2m_seq}=`head -n 2 $h{a2m} | grep -v \
\">\"`;chomp ($h{a2m_seq});\n    $h{a3m_seq}=`head\
 -n 2 $h{a3m} | grep -v \">\"`;chomp ($h{a3m_seq})\
;\n    $h{cons}=$h{a2m_seq};\n    $h{seq}=`head -n\
 2 $h{aln} | grep -v \">\"`;chomp ($h{seq});\n    \
\n    \n\n    @c=split (//, $h{cons});\n    $h{len\
}=$#c+1;\n    for ($n=0,$a=0, $b=0; $a<$h{len};$a+\
+)\n      {\n	if ( $c[$a]=~/[A-Z]/)\n	  {\n	    $h\
{++$n}=++$b;\n\n	  }\n	elsif ( $c[$a]=~/[a-z\\.]/)\
\n	  {\n	    ++$b;\n	  }\n      }\n    \n    $name\
=`head -n 2 $h{aln} | grep \">\"`;\n    $name=~/\\\
>(\\S+)/;\n    $h{name}=$1;\n    \n    `cp $h{a2m}\
 $i.a2m`;\n    `cp $h{a3m} $i.a3m`;\n    `cp $h{al\
n} $i.hh_aln`;\n    \n    return %h;\n  }\n\nsub r\
ead_fasta_seq \n  {\n    my $f=@_[0];\n    my %hse\
q;\n    my (@seq, @com, @name);\n    my ($a, $s,$n\
seq);\n\n    open (F, $f);\n    while (<F>)\n     \
 {\n	$s.=$_;\n      }\n    close (F);\n\n    \n   \
 @name=($s=~/>(\\S*).*\\n[^>]*/g);\n    \n    @seq\
 =($s=~/>.*.*\\n([^>]*)/g);\n    @com =($s=~/>\\S*\
(.*)\\n([^>]*)/g);\n\n    \n    $nseq=$#name+1;\n \
   \n    for ($a=0; $a<$nseq; $a++)\n      {\n	my \
$s;\n	my $n=$name[$a];\n	$hseq{$n}{name}=$n;\n	$se\
q[$a]=~s/[^A-Za-z]//g;\n	$hseq{$n}{order}=$a;\n	$h\
seq{$n}{seq}=$seq[$a];\n	$hseq{$n}{com}=$com[$a];\\
n	\n      }\n    return %hseq;\n  }\n\nsub file_co\
ntains \n  {\n    my ($file, $tag, $max)=(@_);\n  \
  my ($n);\n    $n=0;\n    \n    if ( !-e $file &&\
 ($file =~/$tag/)) {return 1;}\n    elsif ( !-e $f\
ile){return 0;}\n    else \n      {\n	open (FC, \"\
$file\");\n	while ( <FC>)\n	  {\n	    if ( ($_=~/$\
tag/))\n	      {\n		close (FC);\n		return 1;\n	   \
   }\n	    elsif ($max && $n>$max)\n	      {\n		cl\
ose (FC);\n		return 0;\n	      }\n	    $n++;\n	  }\
\n      }\n    close (FC);\n    return 0;\n  }\n	 \
   \n	  \nsub file2string\n  {\n    my $f=@_[0];\n\
    my $string, $l;\n    open (F,\"$f\");\n    whi\
le (<F>)\n      {\n\n	$l=$_;\n	#chomp ($l);\n	$str\
ing.=$l;\n      }\n    close (F);\n    $string=~s/\
\\r\\n//g;\n    $string=~s/\\n//g;\n    return $st\
ring;\n  }\n\n\nsub my_get_opt\n  {\n    my @list=\
@_;\n    my $cl, $a, $argv, @argl;\n    \n    @arg\
l=();\n    $cl=shift @list;\n    for ( $a=0; $a<=$\
#list; $a+=3)\n      {\n	$option=$list[$a];\n	$opt\
ional=$list[$a+1];\n	$status=$list[$a+2];\n	$argv=\
\"\";\n	if ($cl=~/$option(\\S+)/){$argv=$1;}\n	@ar\
gl=(@argl,$argv);\n	\n	\n	#$optional:0=>optional\n\
	#$optional:1=>must be set\n	#$status: 0=>no requi\
rement\n	#$status: 1=>must be an existing file\n	#\
$status: 2=>must be an installed package\n	\n\n	if\
 ($optional==0){;}\n	elsif ( $optional==1 && $argv\
 eq \"\")\n	  {\n\n	    myexit(flush_error(\"Optio\
n $option must be set\"));\n	   \n	  }\n	if ($stat\
us==0){;}\n	elsif ($status ==1 && $argv ne \"\" &&\
 !-e $argv)\n	  {\n	     myexit(flush_error(\"File\
 $argv must exist\"));\n	   \n	  }\n	elsif ( $stat\
us==2 && $argv ne \"\" && &check_pg_is_installed (\
$argv)==0)\n	  {\n	    myexit(flush_error(\"$argv \
is not installed\"));\n	   \n	  }\n      }\n\n    \
return @argl;\n    }\n\nsub tag2value \n  {\n    \\
n    my $tag=(@_[0]);\n    my $word=(@_[1]);\n    \
my $return;\n    \n    $tag=~/$word=\"([^\"]+)\"/;\
\n    $return=$1;\n    return $return;\n  }\n     \
 \nsub hit_tag2pdbid\n  {\n    my $tag=(@_[0]);\n \
   my $pdbid;\n       \n    $tag=~/id=\"(\\S+)\"/;\
\n    $pdbid=$1;\n    $pdbid=~s/_//;\n    return $\
pdbid;\n  }\nsub id2pdbid \n  {\n    my $in=@_[0];\
\n    my $id;\n    \n    $in=~/(\\S+)/;\n    $id=$\
in;\n    \n    if ($id =~/pdb/)\n      {\n	$id=~/p\
db(.*)/;\n	$id=$1;\n      }\n    $id=~s/[|��_]\
//g;\n    return $id;\n  }\nsub set_blast_type \n \
 {\n    my $file =@_[0];\n    if (&file_contains (\
$file,\"EBIApplicationResult\",100)){$BLAST_TYPE=\\
"EBI\";}\n    elsif (&file_contains ($file,\"NCBI_\
BlastOutput\",100)) {$BLAST_TYPE=\"NCBI\";}\n    e\
lse\n      {\n	$BLAST_TYPE=\"\";\n      }\n    ret\
urn $BLAST_TYPE;\n  }\nsub is_valid_blast_xml\n   \
 {\n      my $file=shift;\n      my $line;\n      \
\n      \n      if ( !-e $file) {return 0;}\n     \
 $line=&file2tail ($file,100);\n      if ( $line=~\
/<\\/EBIApplicationResult/ || $line=~/<\\/NCBI_Bla\
stOutput/){return 1;}\n      return 0;\n    }\nsub\
 file2blast_flavor\n      {\n	my $file=shift;\n	if\
 (&file_contains ($file,\"EBIApplicationResult\",1\
00)){return \"EBI\";}\n	elsif (&file_contains ($fi\
le,\"EBIApplicationResult\",100)){return \"NCBI\";\
}\n	else {return \"UNKNOWN\";}\n      }\nsub blast\
_xml2profile \n  {\n    my ($name,$seq,$maxid, $mi\
nid, $mincov, $file)=(@_);\n    my (%p, $a, $strin\
g, $n);\n    \n\n\n    if ($BLAST_TYPE eq \"EBI\" \
|| &file_contains ($file,\"EBIApplicationResult\",\
100)){%p=ebi_blast_xml2profile(@_);}\n    elsif ($\
BLAST_TYPE eq \"NCBI\" || &file_contains ($file,\"\
NCBI_BlastOutput\",100)){%p=ncbi_blast_xml2profile\
(@_);}\n    else \n      {\n	myexit(add_error ( $$\
,$$,getppid(), \"BLAST_FAILURE::unkown XML\",$CL))\
;\n      }\n    for ($a=0; $a<$p{n}; $a++)\n      \
{\n	my $name=$p{$a}{name};\n	$p{$name}{seq}=$p{$a}\
{seq};\n	$p{$name}{index}=$a;\n      }\n    return\
 %p;\n  }\nsub ncbi_tblastx_xml2lib_file \n  {\n  \
  my  ($outlib,$string)=(@_);\n    my ($L,$l, $a,$\
b,$c,$d,$i,$nhits,@identifyerL);\n    my (%ITERATI\
ON);\n      \n    open (F, \">>$outlib\");\n    \n\
    $seq=~s/[^a-zA-Z]//g;\n    $L=length ($seq);\n\
    \n    %ITERATION=xml2tag_list ($string, \"Iter\
ation\");\n    for ($i=0; $i<$ITERATION{n};$i++)\n\
      {\n	my ($qindex, $qlen, %hit, $string);\n	$s\
tring=$ITERATION{$i}{body};\n\n	$qindex=xmltag2val\
ue($string,\"Iteration_iter-num\");\n	$qlen  =xmlt\
ag2value($string,\"Iteration_query-len\");\n	%hit=\
&xml2tag_list  ($string, \"Hit\");\n\n	for ($a=0; \
$a<$hit{n}; $a++)\n	  {\n	    my ($string);\n	    \
$string=$hit{$a}{body};\n	 \n	    $hindex=xmltag2v\
alue($string,\"Hit_accession\")+1;\n	    if ($hind\
ex<=$qindex){next;}\n	    else  {print F  \"# $qin\
dex $hindex\\n\";}\n		   \n	   \n	    $hlen=xmltag\
2value  ($string,\"Hit_len\");\n	    %HSP=&xml2tag\
_list  ($string, \"Hsp\");\n	   \n	    for ($b=0; \
$b<$HSP{n}; $b++)\n	      {\n		my ($string, $qs,$q\
e,$qf,$hs,$he,$hf,$s, $d, $e);\n		$string=$HSP{$b}\
{body};\n	\n		$qs=xmltag2value  ($string,\"Hsp_que\
ry-from\");\n		$qe=xmltag2value  ($string,\"Hsp_qu\
ery-to\");\n		$qf=xmltag2value  ($string,\"Hsp_que\
ry-frame\");\n\n		$hs=xmltag2value  ($string,\"Hsp\
_hit-from\");\n		$he=xmltag2value  ($string,\"Hsp_\
hit-to\");\n		$hf=xmltag2value  ($string,\"Hsp_hit\
-frame\");\n		\n		$s=xmltag2value  ($string,\"Hsp_\
identity\");\n		$l=xmltag2value  ($string,\"Hsp_al\
ign-len\");\n		$s=int(($s*100)/$l);\n		\n		if ($qf\
>0)\n		  {$rqs=$qs; $rqe=$qe;}\n		else\n		  {\n		 \
   $rqe=($qlen-$qs)+1;\n		    $rqs=($qlen-$qe)+1;\\
n		  }\n		\n		if ($hf>0)\n		  {$rhs=$hs; $rhe=$he;\
}\n		else\n		  {\n		    $rhe=($hlen-$hs)+1;\n		   \
 $rhs=($hlen-$he)+1;\n		  }\n		for ($d=0,$e=$rqs; \
$e<$rqe; $e++,$d++)\n		  {\n		    my ($r1,$r2);\n	\
	    $r1=$e;\n		    $r2=$rhs+$d;\n		    print F \"\
 $r1 $r2 $s 0\\n\";\n		  }\n	      }\n	  }\n      \
}\n    print F \"! SEQ_1_TO_N\\n\";\n    \n    clo\
se (F);\n    return %lib;\n  }\n\nsub ncbi_tblastp\
x_xml2lib_file \n  {\n    my  ($outlib,$string,%s)\
=(@_);\n    my ($L,$l, $a,$b,$c,$d,$i,$nhits,@iden\
tifyerL);\n    my (%ITERATION,%hdes, %qdes);\n    \
  \n    open (F, \">>$outlib\");\n    \n    $seq=~\
s/[^a-zA-Z]//g;\n    $L=length ($seq);\n    \n    \
%ITERATION=xml2tag_list ($string, \"Iteration\");\\
n    for ($i=0; $i<$ITERATION{n};$i++)\n      {\n	\
my ($qindex, $qlen, %hit, $string);\n	$string=$ITE\
RATION{$i}{body};\n\n	$qdef=xmltag2value($string,\\
"Iteration_query-def\");\n	%qdes=&tblastpx_name2de\
scription($qdef,%s);\n	$qlen  =xmltag2value($strin\
g,\"Iteration_query-len\");\n	%hit=&xml2tag_list  \
($string, \"Hit\");\n\n	for ($a=0; $a<$hit{n}; $a+\
+)\n	  {\n	    my ($string);\n	    $string=$hit{$a\
}{body};\n	    $hdef=xmltag2value($string,\"Hit_de\
f\");\n	    %hdes=&tblastpx_name2description($hdef\
,%s);\n	    if ($hdes{index}<=$qdes{index}){next;}\
\n	    else  {print F  \"# $qdes{index} $hdes{inde\
x}\\n\";}\n		   \n	   \n	    $hlen=xmltag2value  (\
$string,\"Hit_len\");\n	    %HSP=&xml2tag_list  ($\
string, \"Hsp\");\n	   \n	    for ($b=0; $b<$HSP{n\
}; $b++)\n	      {\n		my ($string, $l,$qs,$qe,$qf,\
$hs,$he,$hf,$s, $d, $e, @s1, @s2);\n		$string=$HSP\
{$b}{body};\n	\n		$qs=xmltag2value  ($string,\"Hsp\
_query-from\");\n		$qe=xmltag2value  ($string,\"Hs\
p_query-to\");\n		$qf=$qdes{frame};\n		$qseq=xmlta\
g2value  ($string,\"Hsp_qseq\");\n		\n		$hs=xmltag\
2value  ($string,\"Hsp_hit-from\");\n		$he=xmltag2\
value  ($string,\"Hsp_hit-to\");\n		$hf=$hdes{fram\
e};\n		$hseq=xmltag2value  ($string,\"Hsp_hseq\");\
\n		\n		$s=xmltag2value  ($string,\"Hsp_identity\"\
);\n		$l=xmltag2value  ($string,\"Hsp_align-len\")\
;\n		$s=int(($s*100)/$l);\n		@s1=tblastpx_hsp2coor\
dinates($qseq,$qs,$qe,%qdes);\n		@s2=tblastpx_hsp2\
coordinates($hseq,$hs,$he,%hdes);\n		\n		\n		for (\
$f=0; $f<=$#s1; $f++)\n		  {\n		    if ($s1[$f]==-\
1 || $s2[$f]==-1){next;}\n		    else \n		      {\n\
			print F \" $s1[$f] $s2[$f] $s 0\\n\";\n		      \
}\n		  }\n	      }\n	  }\n      }\n    print F \"!\
 SEQ_1_TO_N\\n\";\n    \n    close (F);\n    retur\
n %lib;\n  }\nsub tblastpx_hsp2coordinates\n  {\n \
   my ($seq, $s, $e, %des)=@_;\n    my @list;\n   \
 my @sa;\n    my @gap=(-1,-1,-1);\n    \n    $s=$d\
es{start}+3*($s-1);\n  \n    if ($des{strand} eq \\
"d\"){;}\n    else {$s=($des{length}-$s)+1;}\n    \
\n    foreach $c (split (//,$seq))\n      {\n	if (\
 $c eq '-'){push (@list,@gap);}\n	elsif ($des{stra\
nd} eq \"d\")\n	  {\n	    push(@list,$s++,$s++,$s+\
+);\n	  }\n	else\n	  {\n	    push(@list, $s--,$s--\
,$s--);\n	  }\n      }\n    return @list;\n  }\n\n\
sub tblastpx_name2description\n  {\n    my ($name,\
 %s)=@_;\n    my @at=split(\"__\", $name);\n    my\
 %des;\n\n    $des{name}=$at[0];\n    $des{strand}\
=$at[1];\n    \n    $des{start}=$at[2];\n    $des{\
end}=$at[3];\n    $des{length}=$at[4];\n    $des{i\
ndex}=$s{$at[0]}{order}+1;\n    return %des;\n  } \
 \nsub ncbi_blast_xml2profile \n  {\n    my ($name\
,$seq,$maxid, $minid, $mincov, $string)=(@_);\n   \
 my ($L,$l, $a,$b,$c,$d,$nhits,@identifyerL);\n   \
 \n    \n    $seq=~s/[^a-zA-Z]//g;\n    $L=length \
($seq);\n    \n    %hit=&xml2tag_list ($string, \"\
Hit\");\n    \n    \n    for ($nhits=0,$a=0; $a<$h\
it{n}; $a++)\n      {\n	my ($ldb,$id, $identity, $\
expectation, $start, $end, $coverage, $r);\n	my (%\
ID,%DE,%HSP);\n	\n	$ldb=\"\";\n\n	%ID=&xml2tag_lis\
t ($hit{$a}{body}, \"Hit_id\");\n	$identifyer=$ID{\
0}{body};\n	\n	%DE=&xml2tag_list ($hit{$a}{body}, \
\"Hit_def\");\n	$definition=$DE{0}{body};\n	\n	%HS\
P=&xml2tag_list ($hit{$a}{body}, \"Hsp\");\n	for (\
$b=0; $b<$HSP{n}; $b++)\n	  {\n	    my (%START,%EN\
D,%E,%I,%Q,%M);\n\n	 \n	    %START=&xml2tag_list (\
$HSP{$b}{body}, \"Hsp_query-from\");\n	    %HSTART\
=&xml2tag_list ($HSP{$b}{body}, \"Hsp_hit-from\");\
\n	    \n	    %LEN=  &xml2tag_list ($HSP{$b}{body}\
, \"Hsp_align-len\");\n	    %END=  &xml2tag_list (\
$HSP{$b}{body}, \"Hsp_query-to\");\n	    %HEND=  &\
xml2tag_list ($HSP{$b}{body}, \"Hsp_hit-to\");\n	 \
   %E=&xml2tag_list     ($HSP{$b}{body}, \"Hsp_eva\
lue\");\n	    %I=&xml2tag_list     ($HSP{$b}{body}\
, \"Hsp_identity\");\n	    %Q=&xml2tag_list     ($\
HSP{$b}{body}, \"Hsp_qseq\");\n	    %M=&xml2tag_li\
st     ($HSP{$b}{body}, \"Hsp_hseq\");\n	    \n	  \
  for ($e=0; $e<$Q{n}; $e++)\n\n	      {\n		$qs=$Q\
{$e}{body};\n		$ms=$M{$e}{body};\n		\n		$expectati\
on=$E{$e}{body};\n		$identity=($LEN{$e}{body}==0)?\
0:$I{$e}{body}/$LEN{$e}{body}*100;\n		$start=$STAR\
T{$e}{body};\n		$end=$END{$e}{body};\n		$Hstart=$H\
START{$e}{body};\n		$Hend=$HEND{$e}{body};\n	\n		$\
coverage=(($end-$start)*100)/$L;\n	\n		if ($identi\
ty>$maxid || $identity<$minid || $coverage<$mincov\
){next;}\n		@lr1=(split (//,$qs));\n		@lr2=(split \
(//,$ms));\n		$l=$#lr1+1;\n		for ($c=0;$c<$L;$c++)\
{$p[$nhits][$c]=\"-\";}\n		for ($d=0,$c=0; $c<$l; \
$c++)\n		  {\n		    $r=$lr1[$c];\n		    if ( $r=~/\
[A-Za-z]/)\n		      {\n			\n			$p[$nhits][$d + $st\
art-1]=$lr2[$c];\n			$d++;\n		      }\n		  }\n		$Q\
seq[$nhits]=$qs;\n		$Hseq[$nhits]=$ms;\n		$QstartL\
[$nhits]=$start;\n		$HstartL[$nhits]=$Hstart;\n		$\
identityL[$nhits]=$identity;\n		$endL[$nhits]=$end\
;\n		$definitionL[$nhits]=$definition;\n		$identif\
yerL[$nhits]=$identifyer;\n		$comment[$nhits]=\"$l\
db|$identifyer [Eval=$expectation][id=$identity%][\
start=$Hstart end=$Hend]\";\n		$nhits++;\n	      }\
\n	  }\n      }\n    \n    $profile{n}=0;\n    $pr\
ofile{$profile{n}}{name}=$name;\n    $profile{$pro\
file{n}}{seq}=$seq;\n    $profile {n}++;\n    \n  \
  for ($a=0; $a<$nhits; $a++)\n      {\n	$n=$a+1;\\
n	\n	$profile{$n}{name}=\"$name\\_$a\";\n	$profile\
{$n}{seq}=\"\";\n	$profile{$n}{Qseq}=$Qseq[$a];\n	\
$profile{$n}{Hseq}=$Hseq[$a];\n	$profile{$n}{Qstar\
t}=$QstartL[$a];\n	$profile{$n}{Hstart}=$HstartL[$\
a];\n	$profile{$n}{identity}=$identityL[$a];\n	$pr\
ofile{$n}{definition}=$definitionL[$a];\n	$profile\
{$n}{identifyer}=$identifyerL[$a];\n	$profile{$n}{\
comment}=$comment[$a];\n\n	for ($b=0; $b<$L; $b++)\
\n	  {\n	    if ($p[$a][$b])\n	      {\n		$profile\
{$n}{seq}.=$p[$a][$b];\n	      }\n	    else\n	    \
  {\n		$profile{$n}{seq}.=\"-\";\n	      }\n	  }\n\
      }\n    \n    $profile{n}=$nhits+1;\n    retu\
rn %profile;\n  }\nsub ebi_blast_xml2profile \n  {\
\n    my ($name,$seq,$maxid, $minid, $mincov, $str\
ing)=(@_);\n    my ($L,$l, $a,$b,$c,$d,$nhits,@ide\
ntifyerL,$identifyer);\n    \n\n    \n    $seq=~s/\
[^a-zA-Z]//g;\n    $L=length ($seq);\n    %hit=&xm\
l2tag_list ($string, \"hit\");\n    \n    for ($nh\
its=0,$a=0; $a<$hit{n}; $a++)\n      {\n	my ($ldb,\
$id, $identity, $expectation, $start, $end, $cover\
age, $r);\n	my (%Q,%M,%E,%I);\n	\n	$ldb=&tag2value\
 ($hit{$a}{open}, \"database\");\n	$identifyer=&ta\
g2value ($hit{$a}{open}, \"id\");\n\n	$description\
=&tag2value ($hit{$a}{open}, \"description\");\n	\\
n	%Q=&xml2tag_list ($hit{$a}{body}, \"querySeq\");\
\n	%M=&xml2tag_list ($hit{$a}{body}, \"matchSeq\")\
;\n	%E=&xml2tag_list ($hit{$a}{body}, \"expectatio\
n\");\n	%I=&xml2tag_list ($hit{$a}{body}, \"identi\
ty\");\n	\n\n	for ($b=0; $b<$Q{n}; $b++)\n	  {\n\n\
	    $qs=$Q{$b}{body};\n	    $ms=$M{$b}{body};\n	 \
   \n	    $expectation=$E{$b}{body};\n	    $identi\
ty=$I{$b}{body};\n	    \n	    	    \n	    $start=&\
tag2value ($Q{$b}{open}, \"start\");\n	    $end=&t\
ag2value ($Q{$b}{open}, \"end\");\n	    $startM=&t\
ag2value ($M{$b}{open}, \"start\");\n	    $endM=&t\
ag2value ($M{$b}{open}, \"end\");\n	    $coverage=\
(($end-$start)*100)/$L;\n	    \n	   # print \"$id:\
 ID: $identity COV: $coverage [$start $end]\\n\";\\
n	    \n	    \n	    if ($identity>$maxid || $ident\
ity<$minid || $coverage<$mincov){next;}\n	    # pr\
int \"KEEP\\n\";\n\n	    \n	    @lr1=(split (//,$q\
s));\n	    @lr2=(split (//,$ms));\n	    $l=$#lr1+1\
;\n	    for ($c=0;$c<$L;$c++){$p[$nhits][$c]=\"-\"\
;}\n	    for ($d=0,$c=0; $c<$l; $c++)\n	      {\n	\
	$r=$lr1[$c];\n		if ( $r=~/[A-Za-z]/)\n		  {\n		  \
  \n		    $p[$nhits][$d + $start-1]=$lr2[$c];\n		 \
   $d++;\n		  }\n	      }\n	  \n	    $Qseq[$nhits]\
=$qs;\n	    $Hseq[$nhits]=$ms;\n	    $QstartL[$nhi\
ts]=$start;\n	    $HstartL[$nhits]=$Hstart;\n	    \
$identityL[$nhits]=$identity;\n	    $endL[$nhits]=\
$end;\n	    $definitionL[$nhits]=$definition;\n	  \
  $identifyerL[$nhits]=$identifyer;\n	    $comment\
[$nhits]=\"$ldb|$identifyer [Eval=$expectation][id\
=$identity%][start=$startM end=$endM]\";\n	    $nh\
its++;\n	  }\n      }\n    \n    $profile{n}=0;\n \
   $profile{$profile{n}}{name}=$name;\n    $profil\
e{$profile{n}}{seq}=$seq;\n    $profile {n}++;\n  \
  \n    for ($a=0; $a<$nhits; $a++)\n      {\n	$n=\
$a+1;\n	$profile{$n}{name}=\"$name\\_$a\";\n	$prof\
ile{$n}{seq}=\"\";\n	$profile{$n}{Qseq}=$Qseq[$a];\
\n	$profile{$n}{Hseq}=$Hseq[$a];\n	$profile{$n}{Qs\
tart}=$QstartL[$a];\n	$profile{$n}{Hstart}=$Hstart\
L[$a];\n	$profile{$n}{identity}=$identityL[$a];\n	\
$profile{$n}{definition}=$definitionL[$a];	\n	$pro\
file{$n}{identifyer}=$identifyerL[$a];\n	$profile{\
$n}{comment}=$comment[$a];\n\n	for ($b=0; $b<$L; $\
b++)\n	  {\n	    if ($p[$a][$b])\n	      {\n		$pro\
file{$n}{seq}.=$p[$a][$b];\n	      }\n	    else\n	\
      {\n		$profile{$n}{seq}.=\"-\";\n	      }\n	 \
 }\n      }\n    $profile{n}=$nhits+1;\n    \n    \
return %profile;\n  }\nsub output_profile\n  {\n  \
  my ($outfile,$profileR, $trim)=(@_);\n    my ($a\
);\n    my %profile=%$profileR;\n    my $P= new Fi\
leHandle;\n    my $tmp=vtmpnam();\n    \n    open \
($P, \">$tmp\");\n    for ($a=0; $a<$profile{n}; $\
a++)\n      {\n	print $P \">$profile{$a}{name} $pr\
ofile{$a}{comment}\\n$profile{$a}{seq}\\n\";\n    \
  }\n    close ($P);\n\n    if ( $trim)\n      {\n\
	&safe_system (\"t_coffee -other_pg seq_reformat -\
in $tmp -action +trim _aln_%%$trim\\_K1 -output fa\
sta_aln -out $outfile\");\n      }\n    else\n    \
  {\n	&safe_system (\"mv $tmp $outfile\");\n      \
}\n    return;\n  }\nsub blast_xml2hit_list\n  {\n\
    my $string=(@_[0]);\n    return &xml2tag_list \
($string, \"hit\");\n  }\nsub xmltag2value\n  {\n \
   my ($string_in, $tag)=@_;\n    my %TAG;\n    %T\
AG=xml2tag_list ($string_in, $tag);\n    return $T\
AG{0}{body};\n  }\n      \nsub xml2tag_list  \n  {\
\n    my ($string_in,$tag)=@_;\n    my $tag_in, $t\
ag_out;\n    my %tag;\n    \n    if (-e $string_in\
)\n      {\n	$string=&file2string ($string_in);\n \
     }\n    else\n      {\n	$string=$string_in;\n \
     }\n    $tag_in1=\"<$tag \";\n    $tag_in2=\"<\
$tag>\";\n    $tag_out=\"/$tag>\";\n    $string=~s\
/>/>##1/g;\n    $string=~s/</##2</g;\n    $string=\
~s/##1/<#/g;\n    $string=~s/##2/#>/g;\n    @l=($s\
tring=~/(\\<[^>]+\\>)/g);\n    $tag{n}=0;\n    $in\
=0;$n=-1;\n  \n \n\n    foreach $t (@l)\n      {\n\
\n	$t=~s/<#//;\n	$t=~s/#>//;\n	\n	if ( $t=~/$tag_i\
n1/ || $t=~/$tag_in2/)\n	  {\n	 \n	    $in=1;\n	  \
  $tag{$tag{n}}{open}=$t;\n	    $n++;\n	    \n	  }\
\n	elsif ($t=~/$tag_out/)\n	  {\n	    \n\n	    $ta\
g{$tag{n}}{close}=$t;\n	    $tag{n}++;\n	    $in=0\
;\n	  }\n	elsif ($in)\n	  {\n	   \n	    $tag{$tag{\
n}}{body}.=$t;\n	  }\n      }\n  \n    return %tag\
;\n  }\n\n\nsub seq2gor_prediction \n  {\n    my (\
$name, $seq,$infile, $outfile, $gor_seq, $gor_obs)\
=(@_);\n    my ($l);\n    \n    `gorIV -prd $infil\
e -seq $gor_seq -obs $gor_obs > gor_tmp`;\n    ope\
n (GR, \">$outfile\");\n    open (OG, \"gor_tmp\")\
;\n\n    while (<OG>)\n      {\n	\n	$l=$_;\n	if ($\
l=~/\\>/){print GR \"$l\";}\n	elsif ( $l=~/Predict\
ed Sec. Struct./)\n	  {\n	    $l=~s/Predicted Sec.\
 Struct\\.//;\n	    print GR \"$l\";\n	  }\n      \
}\n    close (GR);\n    close (OG);\n    return;\n\
  }\nsub seq2msa_tm_prediction \n  {\n    my ($nam\
e, $seq, $db, $infile, $outfile, $arch, $psv)=(@_)\
;\n    my (%p,%gseq,%R, $blast_output, %s, $l);\n \
   my $R2=new FileHandle;\n    my $db=\"uniprot\";\
\n    my $method=\"psitm\";\n    my $SERVER=\"EBI\\
";\n    \n    $blast_output=&run_blast ($name,\"bl\
astp\", $db, $infile, \"outfile\");\n    \n    if \
(&cache_file(\"GET\",$infile,$name,$method,$db,$ou\
tfile,$SERVER))\n      {\n	print \"\\tPSITM: USE C\
ache\\n\";\n	return $outfile;\n      }\n    else\n\
      {\n	$CACHE_STATUS=\"COMPUTE CACHE\";\n	%p=bl\
ast_xml2profile($name,$seq,$maxid, $minid,$mincov,\
$blast_output);\n	\n	\n	open (F, \">tm_input\");\n\
	for (my $a=0; $a<$p{n}; $a++)\n	  {\n	    my $s;\\
n	    \n	    $s=$p{$a}{seq};\n	    $s=uc($s);\n	  \
  print F \">$p{$a}{name}\\n$s\\n\";\n	    #print \
stdout \">$p{$a}{name}\\n$s\\n\";\n	  }\n	close (F\
);\n	print \"\\tPSITM: kept  $p{n} Homologues for \
Sequence $p{0}{name}\\n\";\n	&safe_system (\"t_cof\
fee -other_pg fasta_seq2hmmtop_fasta.pl -in=tm_inp\
ut -out=$outfile -output=cons -cov=70 -trim=95 -ar\
ch=$arch -psv=$psv\");\n	unlink (\"tm_input\");\n	\
&cache_file(\"SET\",$infile,$name,$method,$db,$out\
file,$SERVER);\n	return;\n      }\n  }\n\n\nsub se\
q2msa_gor_prediction \n  {\n    my ($name, $seq,$i\
nfile, $outfile, $gor_seq, $gor_obs)=(@_);\n    my\
 (%p,%gseq,%R, $blast_output, %s, $l);\n    my $R2\
=new FileHandle;\n    my $db=\"uniprot\";\n    my \
$method=\"psigor\";\n    my $SERVER=\"EBI\";\n    \
\n    $blast_output=&run_blast ($name,\"blastp\", \
\"uniprot\", $infile, \"outfile\");\n    \n    if \
(&cache_file(\"GET\",$infile,$name,$method,$db,$ou\
tfile,$SERVER))\n      {\n	print \"\\tPSIGOR: USE \
Cache\\n\";\n	return $outfile;\n      }\n    else\\
n      {\n	$CACHE_STATUS=\"COMPUTE CACHE\";\n	%p=b\
last_xml2profile($name,$seq,$maxid, $minid,$mincov\
,$blast_output);\n	\n	\n	open (F, \">gor_input\");\
\n	for (my $a=0; $a<$p{n}; $a++)\n	  {\n	    my $s\
;\n	    \n	    $s=$p{$a}{seq};\n	    $s=uc($s);\n	\
    print F \">$p{$a}{name}\\n$s\\n\";\n	    #prin\
t stdout \">$p{$a}{name}\\n$s\\n\";\n	  }\n	close \
(F);\n	print \"\\tGORTM: kept  $p{n} Homologues fo\
r Sequence $p{0}{name}\\n\";\n	&safe_system (\"t_c\
offee -other_pg fasta_seq2hmmtop_fasta.pl -in=gor_\
input -out=$outfile -output=cons -cov=70 -trim=95 \
-gor_seq=$gor_seq -gor_obs=$gor_obs -mode=gor\");\\
n	unlink (\"tm_input\");\n	&cache_file(\"SET\",$in\
file,$name,$method,$db,$outfile,$SERVER);\n	return\
;\n      }\n  }\n\n\n\nsub run_blast\n  {\n    my \
($name, $method, $db, $infile, $outfile, $run)=(@_\
);\n    if (!$run){$run=1;}\n    \n    \n    if (&\
cache_file(\"GET\",$infile,$name,$method,$db,$outf\
ile,$SERVER) && is_valid_blast_xml ($outfile))\n  \
    {return $outfile;}\n    else\n      {\n	$CACHE\
_STATUS=\"COMPUTE CACHE\";\n	if ( $SERVER eq \"EBI\
_SOAP\")\n	  {\n	    &check_configuration (\"EMAIL\
\",\"SOAP::Light\",\"INTERNET\");\n	    \n	    $cl\
_method=$method;\n	    if ($cl_method =~/wu/)\n	  \
    {\n		$cl_method=~s/wu//;\n		if ( $cl_method eq\
 \"psiblast\")\n		  {\n		    add_warning($$,$$,\"P\
SI BLAST cannot be used with the wuBLAST Client. U\
se server=EBI Or server=LOCAL. blastp will be used\
 instead\");\n		    $cl_method=\"blastp\";\n		  }\\
n		\n		$command=\"t_coffee -other_pg wublast.pl --\
email $EMAIL $infile -D $db -p $cl_method --outfil\
e $outfile -o xml>/dev/null 2>/dev/null\";\n		&saf\
e_system ( $command);\n		if (-e \"$outfile.xml\") \
{`mv $outfile.xml $outfile`;}\n	      }\n	    else\
\n	      {\n		if ($cl_method eq \"psiblast\"){$cl_\
method =\"blastp -j5\";}\n		\n		$command=\"t_coffe\
e -other_pg blastpgp.pl --email $EMAIL $infile -d \
$db --outfile $outfile -p $cl_method --mode PSI-Bl\
ast>/dev/null 2>/dev/null\";\n		&safe_system ( $co\
mmand);\n		\n		if (-e \"$outfile.xml\") {`mv $outf\
ile.xml $outfile`;}\n	      }\n	  }\n	elsif ($SERV\
ER eq \"EBI_REST\" || $SERVER eq \"EBI\")\n	  {\n	\
   \n	    $cl_method=$method;\n	    &check_configu\
ration(\"EMAIL\",\"XML::Simple\", \"INTERNET\");\n\
	    if ($db eq \"uniprot\"){$db1=\"uniprotkb\";}\\
n	    else {$db1=$db;}\n	    \n\n	    if ($cl_meth\
od =~/wu/)\n	      {\n		$cl_method=~s/wu//;\n		if \
( $cl_method eq \"psiblast\")\n		  {\n		    $cl_me\
thod=\"blastp\";\n		  }\n		\n		$command=\"t_coffee\
 -other_pg wublast_lwp.pl --email $EMAIL -D $db1 -\
p $cl_method --outfile $outfile --outformat xml --\
stype protein $infile>/dev/null 2>/dev/null\";\n		\
&safe_system ( $command,5);\n		if (-e \"$outfile.x\
ml\") {`mv $outfile.xml $outfile`;}\n		elsif (-e \\
"$outfile.xml.xml\"){`mv $outfile.xml.xml $outfile\
`;}\n	      }\n	    else\n	      {\n		if ( $cl_met\
hod =~/psiblast/){$cl_method =\"blastp -j5\";}\n		\
$command=\"t_coffee -other_pg ncbiblast_lwp.pl --e\
mail $EMAIL -D $db1 -p $cl_method --outfile $outfi\
le --outformat xml --stype protein $infile>/dev/nu\
ll 2>/dev/null\";\n		#$command=\"t_coffee -other_p\
g ncbiblast_lwp.pl --email $EMAIL -D $db1 -p $cl_m\
ethod --outfile $outfile --outformat xml --stype p\
rotein $infile>/dev/null\";\n		&safe_system ( $com\
mand,5);\n		if (-e \"$outfile.xml\") {`mv $outfile\
.xml $outfile`;}\n		elsif (-e \"$outfile.xml.xml\"\
){`mv $outfile.xml.xml $outfile`;} \n		\n	      }\\
n	    \n	 }\n	elsif ($SERVER eq \"NCBI\")\n	  {\n	\
    &check_configuration (\"blastcl3\",\"INTERNET\\
");\n	    if ($db eq \"uniprot\"){$cl_db=\"nr\";}\\
n	    else {$cl_db=$db;}\n	    \n	    if ( $method\
 eq \"psiblast\")\n	      {\n		add_warning($$,$$,\\
"PSI BLAST cannot be used with the NCBI BLAST Clie\
nt. Use server=EBI Or server=LOCAL. blastp will be\
 used instead\");\n		$cl_method=\"blastp\";\n	    \
  }\n	    else\n	      {\n		$cl_method=$method;\n	\
      }\n	    $command=\"blastcl3 -p $cl_method -d\
 $cl_db -i $infile -o $outfile -m 7\";\n	    &safe\
_system ($command);\n	  }\n	elsif ($SERVER =~/CLIE\
NT_(.*)/)\n	  {\n	    my $client=$1;\n	    $comman\
d=\"$client -p $method -d $db -i $infile -o $outfi\
le -m 7\";\n	    &safe_system ($command);\n	  }\n	\
elsif ( $SERVER eq \"LOCAL_blastall\")\n	  {\n	   \
 &check_configuration (\"blastall\");\n	    if ($m\
ethod eq \"blastp\")\n	      {\n		$command=\"blast\
all -d $db -i $infile -o $outfile -m7 -p blastp\";\
\n	      }\n	    &safe_system ($command);\n	  }\n	\
elsif ( $SERVER eq \"LOCAL\")\n	  {\n	    \n	    i\
f ($ENV{\"BLAST_DB_DIR\"})\n	      {\n		$x=$ENV{\"\
BLAST_DB_DIR\"};\n		$cl_db=\"$x$db\";\n	      }\n	\
    else\n	      {\n		$cl_db=$db;\n	      }\n	    \
\n	    if ($method eq \"blastp\")\n	      {\n		&ch\
eck_configuration(\"blastpgp\");\n		$command=\"bla\
stpgp -d $cl_db -i $infile -o $outfile -m7 -j1\";\\
n	      }\n	    elsif ($method eq \"psiblast\")\n	\
      {\n		&check_configuration(\"blastpgp\");\n		\
$command=\"blastpgp -d $cl_db -i $infile -o $outfi\
le -m7 -j5\";\n	      }\n	    elsif ($method eq \"\
blastn\")\n	      {\n		&check_configuration(\"blas\
tall\");\n		$command=\"blastall -p blastn -d $cl_d\
b -i $infile -o $outfile -m7 -W6\";\n	      }	\n	 \
   &safe_system ($command);\n	  }\n	else\n	  {\n	 \
   \n	    myexit(add_error (EXIT_FAILURE,$$,$$,get\
ppid(), \"BLAST_FAILURE::UnknownServer\",$CL));\n	\
  }\n	\n	if ( !-e $outfile || !&is_valid_blast_xml\
($outfile))\n	  {\n	    \n	    if ( -e $outfile)\n\
	      {\n		add_warning ($$,$$,\"Corrupted Blast O\
utput (Run $run)\");\n		unlink($outfile);\n	      \
}\n	    \n	    if ( $run==$BLAST_MAX_NRUNS)\n	    \
  {\n	\n		myexit(add_error (EXIT_FAILURE,$$,$$,get\
ppid(), \"BLAST_FAILURE::UnknownReason\", \"$comma\
nd\"));\n	      }\n	    else\n	      {\n	\n		add_w\
arning ($$,$$,\"Blast for $name failed (Run: $run)\
\");\n		\n		return run_blast ($name, $method, $db,\
$infile, $outfile, $run+1);\n	      }\n	  }\n	\n	&\
cache_file(\"SET\",$infile,$name,$method,$db,$outf\
ile,$SERVER);\n	return $outfile;\n      }\n  }\n\n\
sub cache_file\n  {\n    my ($cache_mode,$infile,$\
name,$method,$db, $outfile,$server)=(@_);\n    my \
$cache_file;\n    #Protect names so that they can \
be turned into legal filenames\n    $name=&clean_f\
ile_name ($name);\n\n    if ($db=~/\\//)\n      {\\
n	$db=~/([^\\/]+)$/;\n	$db=$1;\n      }\n    $cach\
e_file_sh=\"$name.$method.$db.$server.tmp\";\n    \
$cache_file=\"$CACHE/$name.$method.$db.$server.tmp\
\";\n    \n    if ($infile ne \"\")\n      {\n	$ca\
che_file_infile_sh=\"$name.$method.$db.$server.inf\
ile.tmp\";\n	$cache_file_infile=\"$CACHE/$name.$me\
thod.$db.$server.infile.tmp\";\n      }\n    \n   \
 if ($cache_mode eq \"GET\")\n      {\n	if ($CACHE\
 eq \"\" || $CACHE eq \"no\" || $CACHE eq \"ignore\
\"  || $CACHE eq \"local\" || $CACHE eq \"update\"\
){return 0;}\n	elsif ( !-d $CACHE)\n	  {\n	    pri\
nt STDERR \"ERROR: Cache Dir: $CACHE Does not Exis\
t\";\n	    return 0;\n	  }\n	else\n	  {\n	    if (\
 -e $cache_file && &fasta_file1_eq_fasta_file2($in\
file,$cache_file_infile)==1)\n	      {\n		`cp $cac\
he_file $outfile`;\n		$CACHE_STATUS=\"READ CACHE\"\
;\n		return 1;\n	      }\n	  }\n      }\n    elsif\
 ($cache_mode eq \"SET\")\n      {\n	if ($CACHE eq\
 \"\" || $CACHE eq \"no\" || $CACHE eq \"ignore\" \
 || $CACHE eq \"local\" || $CACHE eq \"update\"){r\
eturn 0;}\n	elsif ( !-d $CACHE)\n	  {\n	    print \
STDERR \"ERROR: Cache Dir: $CACHE Does not Exist\"\
;\n	    return 0;\n	  }\n	elsif (-e $outfile)\n	  \
{\n	    `cp $outfile $cache_file`;\n	    if ($cach\
e_file_infile ne \"\"){ `cp $infile $cache_file_in\
file`;}\n\n	    #functions for updating the cache\\
n	    #`t_coffee -other_pg clean_cache.pl -file $c\
ache_file_sh -dir $CACHE`;\n	    #`t_coffee -other\
_pg clean_cache.pl -file $cache_file_infile_sh -di\
r $CACHE`;\n	    return 1;\n	  }\n      }\n    $CA\
CHE_STATUS=\"COMPUTE CACHE\";\n    return 0;\n  }\\
nsub file1_eq_file2\n  {\n    my ($f1, $f2)=@_;\n \
   if ( $f1 eq \"\"){return 1;}\n    elsif ( $f2 e\
q \"\"){return 1;}\n    elsif ( !-e $f1){return 0;\
}\n    elsif ( !-e $f2){return 0;}\n    elsif ($f1\
 eq \"\" || $f2 eq \"\" || `diff $f1 $f2` eq \"\")\
{return 1;}\n    \n    return 0;\n  }\nsub clean_f\
ile_name \n  {\n    my $name=@_[0];\n    \n    $na\
me=~s/[^A-Za-z1-9.-]/_/g;\n    return $name;\n  }\\
nsub url2file\n  {\n    my ($address, $out)=(@_);\\
n    \n    if (&pg_is_installed (\"wget\"))\n	{\n	\
  return &safe_system (\"wget $address -O$out >/de\
v/null 2>/dev/null\");\n	}\n    elsif (&pg_is_inst\
alled (\"curl\"))\n      {\n	return &safe_system (\
\"curl $address -o$out >/dev/null 2>/dev/null\");\\
n      }\n    else\n      {\n	myexit(flus_error(\"\
neither curl nor wget are installed. Imnpossible t\
o fectch remote file\"));\n	exit ($EXIT_FAILURE);\\
n      }\n  }\nsub fasta_file1_eq_fasta_file2\n  {\
\n    my ($f1, $f2)=@_;\n    my (%s1, %s2);\n    m\
y @names;\n    %s1=read_fasta_seq ($f1);\n    %s2=\
read_fasta_seq ($f2);\n\n    @names=(keys (%s1));\\
n    \n    foreach $n (keys(%s1))\n      {\n	if ($\
s1{$n}{seq} ne $s2{$n}{seq}){return 0;}\n      } \\
n    \n    foreach $n (keys(%s2))\n      {\n	if ($\
s1{$n}{seq} ne $s2{$n}{seq}){return 0;}\n      }\n\
    return 1;\n  }\n	\n\n\nsub read_template_file\\
n{\n	my $pdb_templates = @_[0];\n	open (TEMP, \"<$\
pdb_templates\");\n	my %temp_h;\n	while (<TEMP>)\n\
{\n		$line = $_;\n 		$line =~/(\\S+)\\s(\\S+)/;\n \
		$temp_h{$1}= $2;\n}\n	close(TEMP);\n	return %tem\
p_h;\n}\n\nsub calc_rna_template\n{\n	my ($mode, $\
infile, $pdbfile, $outfile)=@_;\n	my %s, %h ;\n	my\
 $result;\n	my (@profiles);\n	&set_temporary_dir (\
\"set\",$infile,\"seq.pep\");\n	%s=read_fasta_seq \
(\"seq.pep\");\n	\n	%pdb_template_h = &read_templa\
te_file($pdbfile);\n	my $pdb_chain;\n	open (R, \">\
result.aln\");\n\n\n	#print stdout \"\\n\";\n	fore\
ach $seq (keys(%s))\n	{\n		if ($pdb_template_h{$se\
q} eq \"\")\n		{\n			next;\n		}\n		open (F, \">seq\
file\");\n		print (F \">$s{$seq}{name}\\n$s{$seq}{\
seq}\\n\");\n		close (F);\n		$pdb_chain = $pdb_tem\
plate_h{$seq};\n		$lib_name=\"$s{$seq}{name}.rfold\
\";\n		$lib_name=&clean_file_name ($lib_name);\n		\
\n 		safe_system (\"secondary_struc.py seqfile $CA\
CHE$pdb_chain  $lib_name\");\n		\n		if ( !-e $lib_\
name)\n		{\n		myexit(flush_error(\"RNAplfold faile\
d to compute the secondary structure of $s{$seq}{n\
ame}\"));\n			myexit ($EXIT_FAILURE);\n		}\n		else\
\n		{\n			print stdout \"\\tProcess: >$s{$seq}{nam\
e} _F_ $lib_name\\n\";\n			print R \">$s{$seq}{nam\
e} _F_ $lib_name\\n\";\n		}\n		unshift (@profiles,\
 $lib_name);\n	}\n	close (R);\n	&set_temporary_dir\
 (\"unset\",$mode, $method,\"result.aln\",$outfile\
, @profiles);\n}\n\n\n\nsub seq2rna_pair{\n	my ($m\
ode, $pdbfile1, $pdbfile2, $method, $param, $outfi\
le)=@_;\n	\n	if ($method eq \"runsara.py\")\n	{\n	\
	open(TMP,\"<$pdbfile1\");\n		my $count = 0;\n		my\
 $line;\n		while (<TMP>)\n		{\n			$line = $_;\n			\
if ($count ==1)\n			{\n				last;\n			}\n			$count \
+= 1;\n		}\n\n		\n		$chain1 = substr($line,length(\
$line)-3,1);\n\n		close TMP;\n		open(TMP,\"<$pdbfi\
le2\");\n		my $count = 0;\n		while (<TMP>)\n		{\n	\
		$line = $_;\n			if ($count ==1)\n			{\n				last;\
\n			}\n			$count += 1;\n		}\n		$chain2 = substr($\
line,length($line)-3,1);\n		close TMP;\n\n		$tmp_f\
ile=&vtmpnam();\n	\n		safe_system(\"runsara.py $pd\
bfile1 $chain1 $pdbfile2 $chain2 -s -o $tmp_file -\
-limitation 5000 > /dev/null 2> /dev/null\") == 0 \
or die \"sara did not work $!\\n\";\n		open(TMP,\"\
<$tmp_file\") or die \"cannot open the sara tmp fi\
le:$!\\n\";\n		open(OUT,\">$outfile\") or die \"ca\
nnot open the $outfile file:$!\\n\";\n\n		my $swit\
ch = 0;\n		my $seqNum = 0;\n		foreach my $line (<T\
MP>)\n		{\n			next unless ($line=~/SARAALI/);\n			\
if ($line=~/>/)\n			{\n				$switch =0;\n				print \
OUT \">seq$seqNum\\n\";\n				$seqNum++;				\n			}\\
n			if ($switch < 2){\n				$switch++;\n				next;\n\
			}\n	\n			if ($line =~/REMARK\\s+SARAALI\\s+([^\\
\*]+)\\*/)\n			{\n				my $string = $1;\n				print \
OUT \"$string\\n\";\n			}\n		}\n		close TMP; \n		c\
lose OUT;\n		unlink($tmp_file);\n	}\n}\n\nsub seq2\
tblastx_lib\n  {\n    my ($mode, $infile, $outfile\
)=@_;\n    my (%s, $method,$nseq);\n\n    $method=\
$mode;\n    &set_temporary_dir (\"set\",$infile,\"\
infile\");\n    %s=read_fasta_seq(\"infile\");\n  \
  \n    \n    foreach $seq (keys(%s))\n      {\n	$\
slist[$s{$seq}{order}]=$s{$seq}{seq};\n	$sname[$s{\
$seq}{order}]=$s{$seq}{name};\n	$slen[$s{$seq}{ord\
er}]=length ($s{$seq}{seq});\n      }\n    $nseq=$\
#sname+1;\n    open (F, \">outfile\");\n    print \
F \"! TC_LIB_FORMAT_01\\n\";\n    print F \"$nseq\\
\n\";\n    for ($a=0; $a<$nseq;$a++)\n      {\n	pr\
int F \"$sname[$a] $slen[$a]  $slist[$a]\\n\"\n   \
   }\n    close (F);\n    `formatdb -i infile -p F\
`;\n\n    `blastall -p tblastx -i infile -d infile\
 -m 7 -S1>blast.output`;\n    \n    ncbi_tblastx_x\
ml2lib_file (\"outfile\", file2string (\"blast.out\
put\"));\n    &set_temporary_dir (\"unset\",$mode,\
 $method, \"outfile\",$outfile);\n    myexit ($EXI\
T_SUCCESS);\n    }\nsub seq2tblastpx_lib\n  {\n   \
 my ($mode, $infile, $outfile)=@_;\n    my (%s, $m\
ethod,$nseq);\n    $method=$mode;\n    &set_tempor\
ary_dir (\"set\",$infile,\"infile\");\n    %s=read\
_fasta_seq(\"infile\");\n    \n    foreach $seq (k\
eys(%s))\n      {\n	$slist[$s{$seq}{order}]=$s{$se\
q}{seq};\n	$sname[$s{$seq}{order}]=$s{$seq}{name};\
\n	$slen[$s{$seq}{order}]=length ($s{$seq}{seq});\\
n      }\n    $nseq=$#sname+1;\n    open (F, \">ou\
tfile\");\n    print F \"! TC_LIB_FORMAT_01\\n\";\\
n    print F \"$nseq\\n\";\n    for ($a=0; $a<$nse\
q;$a++)\n      {\n	print F \"$sname[$a] $slen[$a] \
 $slist[$a]\\n\"\n      }\n    close (F);\n    `pr\
intenv > /home/notredame/tmp/ce`;\n    `t_coffee -\
other_pg seq_reformat -in infile -output tblastx_d\
b1 > tblastxdb`;\n    `formatdb -i tblastxdb -p T`\
;\n    `blastall -p blastp -i tblastxdb -d tblastx\
db -m7 >blast.output`;\n    ncbi_tblastpx_xml2lib_\
file (\"outfile\", file2string (\"blast.output\"),\
 %s);\n    &set_temporary_dir (\"unset\",$mode, $m\
ethod, \"outfile\",$outfile);\n    myexit ($EXIT_S\
UCCESS);\n    }\n\n\n    \n\n\n\nsub file2head\n  \
    {\n	my $file = shift;\n	my $size = shift;\n	my\
 $f= new FileHandle;\n	my $line;\n	open ($f,$file)\
;\n	read ($f,$line, $size);\n	close ($f);\n	return\
 $line;\n      }\nsub file2tail\n      {\n	my $fil\
e = shift;\n	my $size = shift;\n	my $f= new FileHa\
ndle;\n	my $line;\n	\n	open ($f,$file);\n	seek ($f\
,$size*-1, 2);\n	read ($f,$line, $size);\n	close (\
$f);\n	return $line;\n      }\n\n\nsub vtmpnam\n  \
    {\n	my $r=rand(100000);\n	my $f=\"file.$r.$$\"\
;\n	while (-e $f)\n	  {\n	    $f=vtmpnam();\n	  }\\
n	push (@TMPFILE_LIST, $f);\n	return $f;\n      }\\
n\nsub myexit\n  {\n    my $code=@_[0];\n    if ($\
CLEAN_EXIT_STARTED==1){return;}\n    else {$CLEAN_\
EXIT_STARTED=1;}\n    ### ONLY BARE EXIT\n    exit\
 ($code);\n  }\nsub set_error_lock\n    {\n      m\
y $name = shift;\n      my $pid=$$;\n\n      \n   \
   &lock4tc ($$,\"LERROR\", \"LSET\", \"$$ -- ERRO\
R: $name $PROGRAM\\n\");\n      return;\n    }\nsu\
b set_lock\n  {\n    my $pid=shift;\n    my $msg= \
shift;\n    my $p=getppid();\n    &lock4tc ($pid,\\
"LLOCK\",\"LRESET\",\"$p$msg\\n\");\n  }\nsub unse\
t_lock\n   {\n     \n    my $pid=shift;\n    &lock\
4tc ($pid,\"LLOCK\",\"LRELEASE\",\"\");\n  }\nsub \
shift_lock\n  {\n    my $from=shift;\n    my $to=s\
hift;\n    my $from_type=shift;\n    my $to_type=s\
hift;\n    my $action=shift;\n    my $msg;\n    \n\
    if (!&lock4tc($from, $from_type, \"LCHECK\", \\
"\")){return 0;}\n    $msg=&lock4tc ($from, $from_\
type, \"LREAD\", \"\");\n    &lock4tc ($from, $fro\
m_type,\"LRELEASE\", $msg);\n    &lock4tc ($to, $t\
o_type, $action, $msg);\n    return;\n  }\nsub iss\
hellpid\n  {\n    my $p=shift;\n    if (!lock4tc (\
$p, \"LLOCK\", \"LCHECK\")){return 0;}\n    else\n\
      {\n	my $c=lock4tc($p, \"LLOCK\", \"LREAD\");\
\n	if ( $c=~/-SHELL-/){return 1;}\n      }\n    re\
turn 0;\n  }\nsub isrootpid\n  {\n    if(lock4tc (\
getppid(), \"LLOCK\", \"LCHECK\")){return 0;}\n   \
 else {return 1;}\n  }\nsub lock4tc\n	{\n	  my ($p\
id,$type,$action,$value)=@_;\n	  my $fname;\n	  my\
 $host=hostname;\n	  \n	  if ($type eq \"LLOCK\"){\
$fname=\"$LOCKDIR/.$pid.$host.lock4tcoffee\";}\n	 \
 elsif ( $type eq \"LERROR\"){ $fname=\"$LOCKDIR/.\
$pid.$host.error4tcoffee\";}\n	  elsif ( $type eq \
\"LWARNING\"){ $fname=\"$LOCKDIR/.$pid.$host.warni\
ng4tcoffee\";}\n	  \n	  if ($debug_lock)\n	    {\n\
	      print STDERR \"\\n\\t---lock4tc(tcg): $acti\
on => $fname =>$value (RD: $LOCKDIR)\\n\";\n	    }\
\n\n	  if    ($action eq \"LCHECK\") {return -e $f\
name;}\n	  elsif ($action eq \"LREAD\"){return fil\
e2string($fname);}\n	  elsif ($action eq \"LSET\")\
 {return string2file ($value, $fname, \">>\");}\n	\
  elsif ($action eq \"LRESET\") {return string2fil\
e ($value, $fname, \">\");}\n	  elsif ($action eq \
\"LRELEASE\") \n	    {\n	      if ( $debug_lock)\n\
		{\n		  my $g=new FileHandle;\n		  open ($g, \">>\
$fname\");\n		  print $g \"\\nDestroyed by $$\\n\"\
;\n		  close ($g);\n		  safe_system (\"mv $fname $\
fname.old\");\n		}\n	      else\n		{\n		  unlink (\
$fname);\n		}\n	    }\n	  return \"\";\n	}\n	\nsub\
 file2string\n	{\n	  my $file=@_[0];\n	  my $f=new\
 FileHandle;\n	  my $r;\n	  open ($f, \"$file\");\\
n	  while (<$f>){$r.=$_;}\n	  close ($f);\n	  retu\
rn $r;\n	}\nsub string2file \n    {\n    my ($s,$f\
ile,$mode)=@_;\n    my $f=new FileHandle;\n    \n \
   open ($f, \"$mode$file\");\n    print $f  \"$s\\
";\n    close ($f);\n  }\n\nBEGIN\n    {\n      sr\
and;\n    \n      $SIG{'SIGUP'}='signal_cleanup';\\
n      $SIG{'SIGINT'}='signal_cleanup';\n      $SI\
G{'SIGQUIT'}='signal_cleanup';\n      $SIG{'SIGILL\
'}='signal_cleanup';\n      $SIG{'SIGTRAP'}='signa\
l_cleanup';\n      $SIG{'SIGABRT'}='signal_cleanup\
';\n      $SIG{'SIGEMT'}='signal_cleanup';\n      \
$SIG{'SIGFPE'}='signal_cleanup';\n      \n      $S\
IG{'SIGKILL'}='signal_cleanup';\n      $SIG{'SIGPI\
PE'}='signal_cleanup';\n      $SIG{'SIGSTOP'}='sig\
nal_cleanup';\n      $SIG{'SIGTTIN'}='signal_clean\
up';\n      $SIG{'SIGXFSZ'}='signal_cleanup';\n   \
   $SIG{'SIGINFO'}='signal_cleanup';\n      \n    \
  $SIG{'SIGBUS'}='signal_cleanup';\n      $SIG{'SI\
GALRM'}='signal_cleanup';\n      $SIG{'SIGTSTP'}='\
signal_cleanup';\n      $SIG{'SIGTTOU'}='signal_cl\
eanup';\n      $SIG{'SIGVTALRM'}='signal_cleanup';\
\n      $SIG{'SIGUSR1'}='signal_cleanup';\n\n\n   \
   $SIG{'SIGSEGV'}='signal_cleanup';\n      $SIG{'\
SIGTERM'}='signal_cleanup';\n      $SIG{'SIGCONT'}\
='signal_cleanup';\n      $SIG{'SIGIO'}='signal_cl\
eanup';\n      $SIG{'SIGPROF'}='signal_cleanup';\n\
      $SIG{'SIGUSR2'}='signal_cleanup';\n\n      $\
SIG{'SIGSYS'}='signal_cleanup';\n      $SIG{'SIGUR\
G'}='signal_cleanup';\n      $SIG{'SIGCHLD'}='sign\
al_cleanup';\n      $SIG{'SIGXCPU'}='signal_cleanu\
p';\n      $SIG{'SIGWINCH'}='signal_cleanup';\n   \
   \n      $SIG{'INT'}='signal_cleanup';\n      $S\
IG{'TERM'}='signal_cleanup';\n      $SIG{'KILL'}='\
signal_cleanup';\n      $SIG{'QUIT'}='signal_clean\
up';\n      \n      our $debug_lock=$ENV{\"DEBUG_L\
OCK\"};\n      \n      \n      \n      \n      for\
each my $a (@ARGV){$CL.=\" $a\";}\n      if ( $deb\
ug_lock ){print STDERR \"\\n\\n\\n********** START\
 PG: $PROGRAM *************\\n\";}\n      if ( $de\
bug_lock ){print STDERR \"\\n\\n\\n**********(tcg)\
 LOCKDIR: $LOCKDIR $$ *************\\n\";}\n      \
if ( $debug_lock ){print STDERR \"\\n --- $$ -- $C\
L\\n\";}\n      \n	     \n      \n      \n    }\ns\
ub flush_error\n  {\n    my $msg=shift;\n    retur\
n add_error ($EXIT_FAILURE,$$, $$,getppid(), $msg,\
 $CL);\n  }\nsub add_error \n  {\n    my $code=shi\
ft;\n    my $rpid=shift;\n    my $pid=shift;\n    \
my $ppid=shift;\n    my $type=shift;\n    my $com=\
shift;\n    \n    $ERROR_DONE=1;\n    lock4tc ($rp\
id, \"LERROR\",\"LSET\",\"$pid -- ERROR: $type\\n\\
");\n    lock4tc ($$, \"LERROR\",\"LSET\", \"$pid \
-- COM: $com\\n\");\n    lock4tc ($$, \"LERROR\",\\
"LSET\", \"$pid -- STACK: $ppid -> $pid\\n\");\n  \
 \n    return $code;\n  }\nsub add_warning \n  {\n\
    my $rpid=shift;\n    my $pid =shift;\n    my $\
command=shift;\n    my $msg=\"$$ -- WARNING: $comm\
and\\n\";\n    print STDERR \"$msg\";\n    lock4tc\
 ($$, \"LWARNING\", \"LSET\", $msg);\n  }\n\nsub s\
ignal_cleanup\n  {\n    print dtderr \"\\n**** $$ \
(tcg) was killed\\n\";\n    &cleanup;\n    exit ($\
EXIT_FAILURE);\n  }\nsub clean_dir\n  {\n    my $d\
ir=@_[0];\n    if ( !-d $dir){return ;}\n    elsif\
 (!($dir=~/tmp/)){return ;}#safety check 1\n    el\
sif (($dir=~/\\*/)){return ;}#safety check 2\n    \
else\n      {\n	`rm -rf $dir`;\n      }\n    retur\
n;\n  }\nsub cleanup\n  {\n    #print stderr \"\\n\
----tc: $$ Kills $PIDCHILD\\n\";\n    #kill (SIGTE\
RM,$PIDCHILD);\n    my $p=getppid();\n    $CLEAN_E\
XIT_STARTED=1;\n    \n    \n    \n    if (&lock4tc\
($$,\"LERROR\", \"LCHECK\", \"\"))\n      {\n	my $\
ppid=getppid();\n	if (!$ERROR_DONE) \n	  {\n	    &\
lock4tc($$,\"LERROR\", \"LSET\", \"$$ -- STACK: $p\
 -> $$\\n\");\n	    &lock4tc($$,\"LERROR\", \"LSET\
\", \"$$ -- COM: $CL\\n\");\n	  }\n      }\n    my\
 $warning=&lock4tc($$, \"LWARNING\", \"LREAD\", \"\
\");\n    my $error=&lock4tc($$,  \"LERROR\", \"LR\
EAD\", \"\");\n    #release error and warning lock\
 if root\n    \n    if (isrootpid() && ($warning |\
| $error) )\n      {\n	\n	print STDERR \"*********\
******* Summary *************\\n$error\\n$warning\\
\n\";\n\n	&lock4tc($$,\"LERROR\",\"RELEASE\",\"\")\
;\n	&lock4tc($$,\"LWARNING\",\"RELEASE\",\"\");\n \
     } \n    \n    \n    foreach my $f (@TMPFILE_L\
IST)\n      {\n	if (-e $f){unlink ($f);} \n      }\
\n    foreach my $d (@TMPDIR_LIST)\n      {\n	clea\
n_dir ($d);\n      }\n    #No More Lock Release\n \
   #&lock4tc($$,\"LLOCK\",\"LRELEASE\",\"\"); #rel\
ease lock \n\n    if ( $debug_lock ){print STDERR \
\"\\n\\n\\n********** END PG: $PROGRAM ($$) ******\
*******\\n\";}\n    if ( $debug_lock ){print STDER\
R \"\\n\\n\\n**********(tcg) LOCKDIR: $LOCKDIR $$ \
*************\\n\";}\n  }\nEND \n  {\n    \n    &c\
leanup();\n  }\n   \n\nsub safe_system \n{\n  my $\
com=shift;\n  my $ntry=shift;\n  my $ctry=shift;\n\
  my $pid;\n  my $status;\n  my $ppid=getppid();\n\
  if ($com eq \"\"){return 1;}\n  \n  \n\n  if (($\
pid = fork ()) < 0){return (-1);}\n  if ($pid == 0\
)\n    {\n      set_lock($$, \" -SHELL- $com (tcg)\
\");\n      exec ($com);\n    }\n  else\n    {\n  \
    lock4tc ($$, \"LLOCK\", \"LSET\", \"$pid\\n\")\
;#update parent\n      $PIDCHILD=$pid;\n    }\n  i\
f ($debug_lock){printf STDERR \"\\n\\t .... safe_s\
ystem (fasta_seq2hmm)  p: $$ c: $pid COM: $com\\n\\
";}\n\n  waitpid ($pid,WTERMSIG);\n\n  shift_lock \
($pid,$$, \"LWARNING\",\"LWARNING\", \"LSET\");\n\\
n  if ($? == $EXIT_FAILURE || lock4tc($pid, \"LERR\
OR\", \"LCHECK\", \"\"))\n    {\n      if ($ntry &\
& $ctry <$ntry)\n	{\n	  add_warning ($$,$$,\"$com \
failed [retry: $ctry]\");\n	  lock4tc ($pid, \"LRE\
LEASE\", \"LERROR\", \"\");\n	  return safe_system\
 ($com, $ntry, ++$ctry);\n	}\n      elsif ($ntry =\
= -1)\n	{\n	  if (!shift_lock ($pid, $$, \"LERROR\\
", \"LWARNING\", \"LSET\"))\n	    {\n	      add_wa\
rning ($$,$$,\"$com failed\");\n	    }\n	  else\n	\
    {\n	      lock4tc ($pid, \"LRELEASE\", \"LERRO\
R\", \"\");\n	    }\n	  return $?;}\n      else\n	\
{\n	  if (!shift_lock ($pid,$$, \"LERROR\",\"LERRO\
R\", \"LSET\"))\n	    {\n	      myexit(add_error (\
$EXIT_FAILURE,$$,$pid,getppid(), \"UNSPECIFIED sys\
tem\", $com));\n	    }\n	}\n    }\n  return $?;\n}\
\n\nsub check_configuration \n    {\n      my @l=@\
_;\n      my $v;\n      foreach my $p (@l)\n	{\n	 \
 \n	  if   ( $p eq \"EMAIL\")\n	    { \n	      if \
( !($EMAIL=~/@/))\n		{\n		add_warning($$,$$,\"Coul\
d Not Use EMAIL\");\n		myexit(add_error ($EXIT_FAI\
LURE,$$,$$,getppid(),\"EMAIL\",\"$CL\"));\n	      \
}\n	    }\n	  elsif( $p eq \"INTERNET\")\n	    {\n\
	      if ( !&check_internet_connection())\n		{\n	\
	  myexit(add_error ($EXIT_FAILURE,$$,$$,getppid()\
,\"INTERNET\",\"$CL\"));\n		}\n	    }\n	  elsif( $\
p eq \"wget\")\n	    {\n	      if (!&pg_is_install\
ed (\"wget\") && !&pg_is_installed (\"curl\"))\n		\
{\n		  myexit(add_error ($EXIT_FAILURE,$$,$$,getpp\
id(),\"PG_NOT_INSTALLED:wget\",\"$CL\"));\n		}\n	 \
   }\n	  elsif( !(&pg_is_installed ($p)))\n	    {\\
n	      myexit(add_error ($EXIT_FAILURE,$$,$$,getp\
pid(),\"PG_NOT_INSTALLED:$p\",\"$CL\"));\n	    }\n\
	}\n      return 1;\n    }\nsub pg_is_installed\n \
 {\n    my @ml=@_;\n    my $r, $p, $m;\n    my $su\
pported=0;\n    \n    my $p=shift (@ml);\n    if (\
$p=~/::/)\n      {\n	if (safe_system (\"perl -M$p \
-e 1\")==$EXIT_SUCCESS){return 1;}\n	else {return \
0;}\n      }\n    else\n      {\n	$r=`which $p 2>/\
dev/null`;\n	if ($r eq \"\"){return 0;}\n	else {re\
turn 1;}\n      }\n  }\n\n\n\nsub check_internet_c\
onnection\n  {\n    my $internet;\n    my $tmp;\n \
   &check_configuration ( \"wget\"); \n    \n    $\
tmp=&vtmpnam ();\n    \n    if     (&pg_is_install\
ed    (\"wget\")){`wget www.google.com -O$tmp >/de\
v/null 2>/dev/null`;}\n    elsif  (&pg_is_installe\
d    (\"curl\")){`curl www.google.com -o$tmp >/dev\
/null 2>/dev/null`;}\n    \n    if ( !-e $tmp || -\
s $tmp < 10){$internet=0;}\n    else {$internet=1;\
}\n    if (-e $tmp){unlink $tmp;}\n\n    return $i\
nternet;\n  }\nsub check_pg_is_installed\n  {\n   \
 my @ml=@_;\n    my $r=&pg_is_installed (@ml);\n  \
  if (!$r && $p=~/::/)\n      {\n	print STDERR \"\\
\nYou Must Install the perl package $p on your sys\
tem.\\nRUN:\\n\\tsudo perl -MCPAN -e 'install $pg'\
\\n\";\n      }\n    elsif (!$r)\n      {\n	myexit\
(flush_error(\"\\nProgram $p Supported but Not Ins\
talled on your system\"));\n      }\n    else\n   \
   {\n	return 1;\n      }\n  }\n\n$program=\"T-COF\
FEE (Version_8.90)\";\n\n","*TC_METHOD_FORMAT_01\n\
******************generic_method.tc_method********\
*****\n*\n*       Incorporating new methods in T-C\
offee\n*       Cedric Notredame 26/08/08\n*\n*****\
**************************************************\
\n*This file is a method file\n*Copy it and adapt \
it to your need so that the method \n*you want to \
use can be incorporated within T-Coffee\n*********\
**********************************************\n* \
                 USAGE                            \
  *\n*********************************************\
**********\n*This file is passed to t_coffee via -\
in:\n*\n*	t_coffee -in Mgeneric_method.method\n*\n\
*	The method is passed to the shell using the foll\
owing\n*call:\n*<EXECUTABLE><PARAM1><IN_FLAG><seq_\
file><PARAM2><OUT_FLAG><outname><PARAM>\n*\n*Conve\
ntions:\n*<FLAG_NAME> 	<TYPE>		<VALUE>\n*<VALUE>:	\
no_name 	<=> Replaced with a space\n*<VALUE>:	&nbs\
p	<=> Replaced with a space\n*\n******************\
*************************************\n*          \
        ALN_MODE                           *\n****\
**************************************************\
*\n*pairwise   ->all Vs all (no self )[(n2-n)/2aln\
]\n*m_pairwise ->all Vs all (no self)[n^2-n]^2\n*s\
_pairwise ->all Vs all (self): [n^2-n]/2 + n\n*mul\
tiple   ->All the sequences in one go\n*\nALN_MODE\
		pairwise\n*\n***********************************\
********************\n*                  OUT_MODE \
                          *\n*********************\
**********************************\n* mode for the\
 output:\n*External methods: \n* aln -> alignmnent\
 File (Fasta or ClustalW Format)\n* lib-> Lib file\
 (TC_LIB_FORMAT_01)\n*Internal Methods:\n* fL -> I\
nternal Function returning a List (Librairie)\n* f\
A -> Internal Function returning an Alignmnent\n*\\
nOUT_MODE		aln\n**********************************\
*********************\n*                  SEQ_TYPE\
                           *\n********************\
***********************************\n*G: Genomic, \
S: Sequence, P: PDB, R: Profile\n*Examples:\n*SEQT\
YPE	S	sequences against sequences (default)\n*SEQT\
YPE	S_P	sequence against structure\n*SEQTYPE	P_P	s\
tructure against structure\n*SEQTYPE	PS	mix of seq\
uences and structure	\n*\nSEQ_TYPE	S\n*\n\n*******\
************************************************\n\
*                COMMAND LINE                     \
    *\n*EXECUTABLE PARAM1 IN_FLAG OUT_FLAG PARAM  \
           *\n************************************\
*******************\n*****************************\
**************************\n*                  EXE\
CUTABLE                         *\n***************\
****************************************\n*name of\
 the executable\n*passed to the shell: executable\\
n*	\nEXECUTABLE	tc_generic_method.pl\n*\n*********\
**********************************************\n* \
                 IN_FLAG                          \
   *\n********************************************\
***********\n*IN_FLAG\n*flag indicating the name o\
f the in coming sequences\n*IN_FLAG S no_name ->no\
 flag\n*IN_FLAG S &bnsp-in&bnsp -> \" -in \"\n*\nI\
N_FLAG		-infile=\n*\n*****************************\
**************************\n*                  OUT\
_FLAG                           *\n***************\
****************************************\n*OUT_FLA\
G\n*flag indicating the name of the out-coming dat\
a\n*same conventions as IN_FLAG\n*OUT_FLAG	S no_na\
me ->no flag\n*if you want to redirect, pass the p\
arameters via PARAM1\n*set OUT_FLAG to >\n*\nOUT_F\
LAG		-outfile=\n*\n*******************************\
************************\n*                  PARAM\
_1                              *\n***************\
****************************************\n*<EXECUT\
ABLE><PARAM1><IN_FLAG><seq_file><PARAM2><OUT_FLAG>\
<outname><PARAM>\n*Parameters sent to the EXECUTAB\
LE and specified *before* IN_FLAG \n*If there is m\
ore than 1 PARAM line, the lines are\n*concatenate\
d\n*Command_line: @EP@PARAM@-gapopen%e10%s-gapext%\
e20\n*	%s white space\n*	%e equal sign\n*\n*PARAM1\
	\n*\n*\n*\n**************************************\
*****************\n*                  PARAM_2     \
                         *\n**********************\
*********************************\n*<EXECUTABLE><P\
ARAM1><IN_FLAG><seq_file><PARAM2><OUT_FLAG><outnam\
e><PARAM>\n*Parameters sent to the EXECUTABLE and \
specified \n*after* IN_FLAG and *before* OUT_FLAG\\
n*If there is more than 1 PARAM line, the lines ar\
e\n*concatenated\n*\n*PARAM1	\n*\n*\n*************\
******************************************\n*     \
             PARAM                              *\\
n*************************************************\
******\n*<EXECUTABLE><PARAM1><IN_FLAG><seq_file><P\
ARAM2><OUT_FLAG><outname><PARAM>\n*Parameters sent\
 to the EXECUTABLE and specified *after* OUT_FLAG\\
n*If there is more than 1 PARAM line, the lines ar\
e\n*concatenated\n*\nPARAM	-mode=seq_msa -method=c\
lustalw\nPARAM   -OUTORDER=INPUT -NEWTREE=core -al\
ign -gapopen=-15\n*\n*****************************\
**************************\n*                  END\
                                *\n***************\
****************************************\n","*TC_M\
ETHOD_FORMAT_01\n***************clustalw_method.tc\
_method*********\nEXECUTABLE	clustalw\nALN_MODE		p\
airwise\nIN_FLAG		-INFILE=\nOUT_FLAG		-OUTFILE=\nO\
UT_MODE		aln\nPARAM		-gapopen=-10\nSEQ_TYPE		S\n**\
***********************************************\n"\
,"$VersionTag =                                   \
                                                  \
                                              2.43\
;\nuse Env;\nuse FileHandle;\nuse Cwd;\nuse File::\
Path;\nuse Sys::Hostname;\nour $PIDCHILD;\nour $ER\
ROR_DONE;\nour @TMPFILE_LIST;\nour $EXIT_FAILURE=1\
;\nour $EXIT_SUCCESS=0;\n\nour $REFDIR=getcwd;\nou\
r $EXIT_SUCCESS=0;\nour $EXIT_FAILURE=1;\n\nour $P\
ROGRAM=\"extract_from_pdb\";\nour $CL=$PROGRAM;\n\\
nour $CLEAN_EXIT_STARTED;\nour $debug_lock=$ENV{\"\
DEBUG_LOCK\"};\nour $LOCKDIR=$ENV{\"LOCKDIR_4_TCOF\
FEE\"};\nif (!$LOCKDIR){$LOCKDIR=getcwd();}\nour $\
ERRORDIR=$ENV{\"ERRORDIR_4_TCOFFEE\"};\nour $ERROR\
FILE=$ENV{\"ERRORFILE_4_TCOFFEE\"};\n&set_lock ($$\
);\nif (isshellpid(getppid())){lock4tc(getppid(), \
\"LLOCK\", \"LSET\", \"$$\\n\");}\n      \nour $SI\
LENT=\" >/dev/null 2>/dev/null\";\nour $INTERNET=-\
1;\n\n\n\n\n\n\nour $BLAST_MAX_NRUNS=2;\nour $EXIT\
_SUCCESS=0;\nour $EXIT_FAILURE=1;\nour $CONFIGURAT\
ION=-1;\nour $REF_EMAIL=\"\";\nour $PROGRAM=\"extr\
act_from_pdb\";\n\n\nmy %onelett_prot=&fill_onelet\
t_prot();\nmy %threelett_prot=&fill_threelett_prot\
();\nmy %onelett_RNA=&fill_onelett_RNA();\nmy %thr\
eelett_RNA=&fill_threelett_RNA();\nmy %onelett_DNA\
=&fill_onelett_DNA();\nmy %threelett_DNA=&fill_thr\
eelett_DNA();\n\n\n\n\n\nmy %onelett = (\n'P' => \\
\%onelett_prot,\n'D' => \\%onelett_DNA,\n'R' => \\\
%onelett_RNA\n);\n\n\nmy %threelett = (\n'P' => \\\
%threelett_prot,\n'D' => \\%threelett_DNA,\n'R' =>\
 \\%threelett_RNA\n);\n\n\n\n\n\n\n\nif($ARGV[0]=~\
/help/ ||$ARGV[0]=~/man/ || $ARGV[0]=~/HELP/ || $A\
RGV[0]=~/Man/ || $ARGV[0] eq \"-h\"  || $ARGV[0] e\
q \"-H\"  )\n{die \"SYNTAX: extract_from_pdb Versi\
on $VersionTag	\n	Minimum:             [extract_fr\
om_pdb file] \n			   OR \n			     [... | extract_f\
rom_pdb]\n 	Flags (Default setting on the first li\
ne)\n	   -version...................[Returns the V\
ersion Number]\n           -force.................\
....[Forces the file to be treated like a PDB file\
]\n                                      [Regenera\
tes the header and SEQRES fields]\n           -for\
ce_name................[Forces the file to be name\
d after name]]\n           -infile.....file.......\
....[Flag can be omited]\n			              [File m\
ust be pdb or fro pgm]\n                          \
            [File can also be compressed Z or gz]\\
n                                      [In the cas\
e of a compressed file, you can omit the gz|Z exte\
nsion]\n           -netfile...................[Fil\
e will be fetch from the net using wget]\n        \
                              [wget or curl must b\
e installed]\n                                    \
  [ftp://ftp.gnu.org/pub/gnu/wget/]\n             \
                         [http://curl.haxx.se/]\n \
                                     [Must also be\
 used to retrieve the file from a local pdb copy (\
cf netaddress)]\n           -netaddress...........\
.....[Address used for the retrieving the netfile]\
\n                                      [http://ww\
w.rcsb.org/pdb/cgi/export.cgi/%%.pdb.gz?format=PDB\
&pdbId=%%&compression=gz]\n                       \
               [http://www.expasy.ch/cgi-bin/get-p\
db-entry.pl?%%]\n                                 \
     [local -> will get the file from pdb_dir (see\
 pdb_dir)]\n           -netcompression............\
[Extension if the netfile comes compressed]\n     \
                                 [gz]\n           \
-pdb_dir...................[address of the reperto\
ry where the pdb is installed]\n                  \
                    [Supports standard ftp style i\
nstallation OR every stru in DIR]\n               \
                       [Give the ..../pdb/structur\
e/ dir]\n                                      [If\
 value omitted, the pg gets it from the env variab\
le PDB_DIR]\n           -netcompression_pg........\
.[gunzip]\n           -is_pdb_name........name...[\
Returns 1 if the name is a PDB ID, 0 otherwise]\n \
          -get_pdb_chains.....name...[Returns the \
list of chains corresponding to the entry]\n      \
     -get_pdb_id.........name...[Returns the PDB i\
d within the provided pdb file]\n           -get_f\
ugue_name.....name...[Turns a name into a name val\
id for fugue]\n                                   \
   [Uses the netaddress to do so]\n	   -chain.....\
.FIRST..........[Extract the first chain only]\n		\
       A B C..........[Extract Several chains if n\
eeded]\n		       ALL............[Extract all the c\
hains]	\n           -ligand.....ALL............[Ex\
tract the ligands in the chain (HETATM)]\n        \
               <name1>,<name2>[Extract All the nam\
ed lignds]\n	   -ligand_only...............[Extrac\
t only the ligands]\n           -ligand_list......\
.........[Extract the list of ligands]\n	   -coor.\
......<start>..<end>.[Coordinates of the fragment \
to extract]\n			              [Omit end to include\
 the Cter]\n           -num........absolute.......\
[absolute: relative to the seq] \n                \
       file...........[file: relative to file]\n  \
         -num_out....new............[new: start 1-\
>L]\n                       old............[old: k\
eep the file coordinates]\n           -delete.....\
<start>..<end>.[Delete from residue start to resid\
ue end]\n	   -atom.......CA.............[Atoms to \
include, ALL for all of them]\n		       CA O N....\
.....[Indicate several atoms if needed]\n	   -code\
.......3..............[Use the 1 letter code or th\
e 3 letters code]\n	   -mode.......raw............\
[Output original pdb file]\n                      \
 pdb............[Output something that looks like \
pdb]\n		       fasta..........[Output the sequence\
s in fasta format]\n		       simple.........[Outpu\
t a format easy to parse in C ]\n            -seq_\
field..ATOM...........[Field used to extract the s\
equence]\n		       SEQRES.........[Use the complet\
e sequence]\n	   -seq.......................[Equiv\
alent to  -mode fasta]\n	   -model......1.........\
.....[Chosen Model in an NMR file]\n           -no\
diagnostic..............[Switches Error Messages o\
ff]\n           -debug.....................[Sets t\
he DEBUG ON]\n           -no_remote_pdb_dir.......\
..[Do not look for a remote file]\n           -cac\
he_pdb.................[Cache Value, default is $H\
OME/.t_coffee/cache, other values: NO<=> No cache]\
\n\n      Environement Variables\n           These\
 variables can be set from the environement\n     \
      Command line values with the corresponding f\
lag superseed evironement value\n           NO_REM\
OTE_PDB_DIR..........[Prevents the program from se\
arching remote file: faster]\n           PDB_DIR..\
..................[Indicates where PDB file must b\
e fetched (localy)]\n\n	 PROBLEMS: please contact \
cedric.notredame\\@europe.com\\n\";\n	 exit ($EXIT\
_SUCCESS);\n}\n\n$np=0;\n$n_para=$#ARGV;\n$model=1\
;\n$pdb_dir=$ENV{'PDB_DIR'};if ($pdb_dir){$pdb_dir\
.=\"/\";}\n$debug=$ENV{'DEBUG_EXTRACT_FROM_PDB'};\\
n\n$no_remote_pdb_dir=$ENV{NO_REMOTE_PDB_DIR};\n$H\
OME=$ENV{'HOME'};\nif ( $ENV{CACHE_4_TCOFFEE})\n{$\
cache=$ENV{CACHE_4_TCOFFEE};}\nelse\n{\n    $cache\
=\"$HOME/.t_coffee/cache/\";\n}\n\n   \n$netaddres\
s=\"http://www.rcsb.org/pdb/files/%%.pdb.gz\";\n$n\
etcompression_pg=\"gunzip\";\n$netcompression=\"gz\
\";\n\n  foreach ($np=0; $np<=$n_para; $np++)\n{  \
      \n    $value=$ARGV[$np];\n    \n    if  ($np\
==0 && !($value=~/^-.*/))\n{ \n       $pdb_file= $\
ARGV[$np];\n}\n    elsif ( !($value=~/^-.*/))\n{\n\
	print \"@ARGV\";\n	die;\n} \n    \n    elsif ($va\
lue eq \"-nodiagnostic\"){$nodiagnostic=1;}\n    e\
lsif ($value eq \"-force\")\n{\n	$force_pdb=1;\n}\\
n    elsif ($value eq \"-force_name\")\n{\n	$force\
_name=$ARGV[++$np];\n	$force_pdb=1;\n}\n    \n    \
elsif ($value eq \"-is_pdb_name\")\n{\n	$pdb_file=\
 $ARGV[++$np];\n	\n	$is_pdb_name=1;\n	\n} \n    el\
sif ($value eq \"-debug\")\n{\n	$debug=1;\n}\n    \
elsif ($value eq \"-get_pdb_chains\")\n{\n	$pdb_fi\
le= $ARGV[++$np];\n	$get_pdb_chains=1;\n}\n    els\
if ($value eq \"-get_pdb_ligands\")\n{\n	$get_pdb_\
ligands=1;\n}\n    \n    elsif ($value eq \"-get_p\
db_id\")\n{\n	$pdb_file= $ARGV[++$np];\n	$get_pdb_\
id=1;\n	\n}\n    \n    elsif ( $value eq \"-get_fu\
gue_name\")\n{\n	$pdb_file= $ARGV[++$np];\n	$get_f\
ugue_name=1;\n}\n    elsif ( $value eq \"-infile\"\
)\n{\n       $pdb_file= $ARGV[++$np];\n} \n    els\
if ($value eq \"-netfile\")\n{\n	$netfile=1;\n	if \
( !($ARGV[$np+1]=~/^-.*/)){$pdb_file= $ARGV[++$np]\
;}\n}\n    elsif (  $value eq \"-num\")\n{\n      \
 $numbering= $ARGV[++$np];\n}\n    elsif (  $value\
 eq \"-num_out\")\n{\n       $numbering_out= $ARGV\
[++$np];\n}\n    elsif ( $value eq \"-netaddress\"\
)\n{\n	$netadress=$ARGV[++$np];\n}\n     \n    els\
if ( $value eq \"-netcompression\")\n{\n	 $netcomp\
ression=$ARGV[++$np];\n}\n    elsif ( $value eq \"\
-pdb_dir\")\n{\n	 if ( !($ARGV[$np+1]=~/^-.*/)){$p\
db_dir= \"$ARGV[++$np]/\";}\n}\n     elsif ( $valu\
e eq \"-no_remote_pdb_dir\")\n{\n	$no_remote_pdb_d\
ir=1;\n	if ( !($ARGV[$np+1]=~/^-.*/)){$pdb_dir= \"\
$ARGV[++$np]/\";}\n}\n    elsif ( $value eq \"-cac\
he\")\n{\n	$cache=$ARGV[++$np];\n}\n    \n    elsi\
f ($value eq \"-netcompression_pg\")\n{\n	  $netco\
mpression_pg=$ARGV[++$np];\n}\n     elsif ($value \
eq \"-mode\")\n{\n       $MODE=$ARGV[++$np];\n}\n\\
n    elsif ( $value eq \"-model\")\n{\n       $mod\
el= $ARGV[++$np];\n}\n    elsif ($value eq \"-seq_\
field\" )\n{\n       $seq_field= $ARGV[++$np];\n} \
  \n    elsif ($value eq \"-coor\" )\n{\n       $s\
tart= $ARGV[++$np];\n  \n       if (($ARGV[$np+1] \
eq \"\") ||($ARGV[$np+1]=~/^-.*/)){$end=\"*\";} \n\
       else {$end=   $ARGV[++$np];}     \n       $\
coor_set=1;\n}\n    elsif ($value eq \"-delete\" )\
\n{\n       $delete_start= $ARGV[++$np];\n       $\
delete_end= $ARGV[++$np];\n       $delete_set=1;\n\
}\n    elsif  ($value eq \"-code\")\n{\n       $co\
de= $ARGV[++$np];\n}\n    elsif  ($value eq \"-no_\
hetatm\")\n{\n       $no_hetatm=1;\n}\n    elsif (\
$value eq \"-chain\")\n{\n       while (!($ARGV[$n\
p+1] eq \"\") &&!($ARGV[$np+1]=~/^-.*/))\n{\n	    \
  ++$np;\n	      @c_chain=(@chain,  $ARGV[$np]);\n\
	      $hc_chain{$ARGV[$np]}=$#c_chain+1;\n}      \
     \n}\n    elsif ($value eq \"-atom\")\n{\n\n  \
     while (!($ARGV[$np+1] eq \"\") && !($ARGV[$np\
+1]=~/^-.*/))\n{\n	      ++$np;\n	      $atom[$n_a\
tom++]=  $ARGV[$np];\n	      $atom_list{$ARGV[$np]\
}=1;	      \n} \n       \n}\n    elsif ( $value eq\
 \"-unfold\")\n{\n	$unfold=1;\n}\n    elsif ($valu\
e eq \"-seq\" ||$value eq \"-fasta\" )\n{\n       \
$MODE=\"fasta\";\n}\n    elsif ( $value eq \"-vers\
ion\")\n{\n	print STDERR  \"\\nextract_from_pdb: V\
ersion $VersionTag\\n\";\n	&myexit ($EXIT_SUCCESS)\
;\n}\n    elsif ( $value eq \"-ligand\")\n{\n	whil\
e (!($ARGV[$np+1] eq \"\") && !($ARGV[$np+1]=~/^-.\
*/))\n{\n	    ++$np;\n	    $ligand=1;\n	    $ligan\
d_list{$ARGV[$np]}=1;	      \n} \n	$hc_chain{'LIGA\
ND'}=1;\n}\n    elsif ( $value eq \"-ligand_only\"\
)\n{\n	$ligand_only=1;\n}\n}\nif ( $debug)\n{\n   \
 print STDERR \"\\n[DEBUG:extract_from_pdb] NO_REM\
OTE_PDB_DIR: $no_remote_pdb_dir\\n\";\n    print S\
TDERR \"\\n[DEBUG:extract_from_pdb] PDB_DIR: $pdb_\
dir\\n\";\n}\n\nif ( $is_pdb_name)\n{\n    if (rem\
ote_is_pdb_name($pdb_file, $netaddress))\n{\n	prin\
t \"1\";\n}\n    else\n{\n	print \"0\";\n}\n    ex\
it ($EXIT_SUCCESS);\n}\n    \n\nif (!$force_name)\\
n{\n    $pdb_file=~/([^\\/]*)$/;\n    $force_name=\
$1;\n}\n\n$local_pdb_file=$pdb_file;\n\nif ( $debu\
g){print STDERR \"\\n[DEBUG: extract_from_pdb] Sca\
n For $local_pdb_file\\n\";}\n\n$mem=$no_remote_pd\
b_dir;\n$no_remote_pdb_dir=1;\n$tmp_pdb_file=get_p\
db_file ($local_pdb_file);\n\nif ( !-e $tmp_pdb_fi\
le || $tmp_pdb_file eq \"\")\n  {\n    $local_pdb_\
file=$pdb_file;\n    ($local_pdb_file, $suffix_cha\
in)=&pdb_name2name_and_chain($local_pdb_file);\n\n\
    if ($local_pdb_file)\n      {\n	if ( $debug){p\
rint STDERR \"\\nSplit $pdb_file into $local_pdb_f\
ile and $suffix_chain \\n\";}\n	$tmp_pdb_file=get_\
pdb_file ($local_pdb_file);\n	if ( $tmp_pdb_file n\
e \"\")\n	  {\n	    @c_chain=();\n	    @c_chain=($\
suffix_chain);\n	    %hc_chain=();\n	    $hc_chain\
{$suffix_chain}=1;\n	  }\n      }\n  }\n\n$no_remo\
te_pdb_dir=$mem;\nif ($no_remote_pdb_dir==0)\n{\n \
   if ( !-e $tmp_pdb_file || $tmp_pdb_file eq \"\"\
)\n{\n	\n	$local_pdb_file=$pdb_file;\n	($local_pdb\
_file, $suffix_chain)=&pdb_name2name_and_chain($lo\
cal_pdb_file);\n	if ($local_pdb_file)\n{\n	    if \
( $debug){print STDERR \"\\nSplit $pdb_file into $\
local_pdb_file and $suffix_chain \\n\";}\n	    $tm\
p_pdb_file=get_pdb_file ($local_pdb_file);    \n	 \
   if ( $tmp_pdb_file ne \"\")\n{\n		@c_chain=();\\
n		@c_chain=($suffix_chain);\n		%hc_chain=();\n		$\
hc_chain{$suffix_chain}=1;\n}\n}\n}\n}\n\nif ( $de\
bug){print STDERR \"\\n$pdb_file copied into ##$tm\
p_pdb_file##\\n\";}\n\n\nif ( !-e $tmp_pdb_file ||\
 $tmp_pdb_file eq \"\")\n{\n	if ($is_pdb_name)\n{\\
n	    print \"0\\n\"; exit ($EXIT_SUCCESS);\n}\n	e\
lse\n{\n  \n	    print \"\\nEXTRACT_FROM_PDB: NO R\
ESULT for $pdb_file\\n\";\n	    &myexit ($EXIT_SUC\
CESS);	\n}\n}\n\n\n\n\n%molecule_type=&pdbfile2cha\
intype($tmp_pdb_file);\nif ( $debug){print STDERR \
\"\\n\\tSequence Type determined\\n\";}\n\n$pdb_id\
=&get_pdb_id ($tmp_pdb_file);\nif ( $debug){print \
STDERR \"\\n\\tID: $pdb_id (for $tmp_pdb_file)\\n\\
";}\n\nif ( $pdb_id eq \"\"){$pdb_id=$force_name;}\
\n\n@f_chain=&get_chain_list ($tmp_pdb_file);\nif \
( $debug){print STDERR \"\\n\\tChain_list:@f_chain\
\\n\";}\n\nif ( $get_pdb_chains)\n{\n    print \"@\
f_chain\\n\";\n    &myexit ($EXIT_SUCCESS);\n}\nif\
 ( $get_pdb_ligands)\n{\n    %complete_ligand_list\
=&get_ligand_list ($tmp_pdb_file);\n    print $com\
plete_ligand_list{\"result\"};\n    &myexit ($EXIT\
_SUCCESS);\n}\n\nelsif ( $get_pdb_id ||$get_fugue_\
name )\n{\n    if    (@c_chain && $c_chain[0] eq \\
"FIRST\"){$pdb_id=$pdb_id.$f_chain[0];}\n    elsif\
 (@c_chain && $c_chain[0] ne \" \"){$pdb_id=$pdb_i\
d.$c_chain[0];}\n    \n    print \"$pdb_id\\n\";\n\
    &myexit ($EXIT_SUCCESS);\n    \n}\nelsif ( $is\
_pdb_name)\n{\n    printf \"1\\n\";\n    &myexit (\
$EXIT_SUCCESS);\n}\n\n\n\n$structure_file=vtmpnam(\
);\n\nif ( $debug){print STDERR \"\\n\\tCheck_poin\
t #1: $tmp_pdb_file  $structure_file\\n\";}\n\n$IN\
FILE=vfopen (\"$tmp_pdb_file\", \"r\"); \n$TMP=vfo\
pen (\"$structure_file\", \"w\");\n\n$print_model=\
1;\n$in_model=0;\n\nif ( $debug){print STDERR \"\\\
n\\tCheck_point #2\\n\";}\nwhile ( <$INFILE>)\n{\n\
  my $first_model=0;\n  $line=$_;\n\n  if ( !$firs\
t_model && ($line =~/^MODEL\\s*(\\d*)/))\n    {\n \
     $first_model=$1;\n      if ($model==1){$model\
=$first_model;}\n    }\n  \n  if (($line =~/^MODEL\
\\s*(\\d*)/))\n    {\n      if ($1==$model)\n	{\n	\
  $in_model=1;\n	  $print_model=1;\n	  $is_nmr=1;\\
n	}\n      elsif ( $in_model==0)\n	{\n	  $print_mo\
del=0;\n	}\n      elsif ( $in_model==1)\n	{\n	  la\
st;\n	}\n    }\n  if ($print_model){print $TMP $li\
ne;}  \n}\nclose ($TMP);\nclose ($INFILE);\n\nif (\
 $debug){print STDERR \"\\n\\tCheck_point #3\\n\";\
}	\n\n  if ($numbering eq \"\"){$numbering=\"absol\
ute\";}\n  if ($numbering_out eq \"\"){$numbering_\
out=\"new\";}\n\n  if ( $delete_set && $coor_set) \
{die \"-delete and -coor are mutually exclusive, s\
orry\\n\";}\n  if ( $n_atom==0){$atom_list[$n_atom\
++]=\"ALL\";$atom_list{$atom_list[0]}=1;}\n  if ( \
$seq_field eq \"\"){$seq_field=\"ATOM\";}\n  \n  i\
f ( $MODE eq \"\"){$MODE=\"pdb\";}\n  elsif ( $MOD\
E eq \"simple\" && $code==0){$code=1;}\n\n  if ( $\
code==0){$code=3;}\n\n\nif ($f_chain[0] eq \" \"){\
$hc_chain{' '}=1;$c_chain[0]=\" \";}\nelsif (!@c_c\
hain){$hc_chain{FIRST}=1;$c_chain[0]=\"FIRST\";}#m\
ake sure the first chain is taken by default\n\nif\
    ($hc_chain{ALL}) \n{\n      @c_chain=@f_chain;\
\n      foreach $e (@c_chain){$hc_chain{$e}=1;}\n}\
\nelsif($hc_chain{FIRST})\n{\n	@c_chain=($f_chain[\
0]);\n	$hc_chain{$f_chain[0]}=1;\n}\n\n\n$MAIN_HOM\
_CODE=&get_main_hom_code ($structure_file);\n$INFI\
LE=vfopen ($structure_file, \"r\");\n\n\nif ( $MOD\
E eq \"raw_pdb\" || $MODE eq \"raw\")\n{\n    whil\
e (<$INFILE>)\n{	print \"$_\";}\n    close ( $INFI\
LE);\n    &myexit($EXIT_SUCCESS);\n}    \nif ( $MO\
DE eq \"raw4fugue\" )\n{\n    while (<$INFILE>)\n{\
	\n	$l=$_;\n	if ($l=~/^SEQRES/)\n{\n	    \n	    $c\
= substr($l,11,1);\n	    if ($hc_chain {$c}){print\
 \"$l\";}\n}\n	elsif ( $l=~/^ATOM/)\n{\n	    $c=su\
bstr($l,21,1);\n	    if ($hc_chain {$c}){print \"$\
l\";}\n}\n}\n    close ( $INFILE);\n    &myexit($E\
XIT_SUCCESS);\n}    \n\nif ( $MODE eq \"pdb\")\n{\\
n\n    $read_header=0;\n    while (<$INFILE>) \n{\\
n	    $line=$_;\n	    if    ($line =~ /^HEADER/){p\
rint \"$line\";$read_header=1;}\n}\n    close ($IN\
FILE);\n\n    if (!$read_header)\n{\n	print \"HEAD\
ER    UNKNOWN                                 00-J\
AN-00   $force_name\\n\";\n}\n\n    $INFILE=vfopen\
 ($structure_file, \"r\");\n    \n    print \"COMP\
ND   1 CHAIN:\";\n    $last=pop(@c_chain);\n    fo\
reach $c ( @c_chain){ print \" $c,\";}\n    if ( $\
last eq \" \"){print \" NULL;\\n\";}\n    else \n{\
\n      print \" $last;\\n\";\n}\n    @c_chain=(@c\
_chain, $last);\n    \n    print \"REMARK Output o\
f the program extract_from_pdb (Version $VersionTa\
g)\\n\";\n    print \"REMARK Legal PDB format not \
Guaranteed\\n\";\n    print \"REMARK This format i\
s not meant to be used in place of the PDB format\\
\n\";\n    print \"REMARK The header refers to the\
 original entry\\n\";\n    print \"REMARK The sequ\
ence from the original file has been taken in the \
field: $seq_field\\n\";\n    print \"REMARK extrac\
t_from_pdb, 2001, 2002, 2003, 2004, 2005 2006 (c) \
CNRS and Cedric Notredame\\n\";   \n    if ( $coor\
_set)\n{\n       print \"REMARK Partial chain: Sta\
rt $start End $end\\n\";\n}\n    if ( $is_nmr)\n{\\
n       print \"REMARK NMR structure: MODEL $model\
\\n\";\n}\n   \n    if ( $n_atom!=0)\n{\n       pr\
int \"REMARK Contains Coordinates of: \";\n       \
foreach $a (@atom){print \"$a \";}\n       print \\
"\\n\";\n}  \n}\n\n\n\n\nmy $residue_index = -999;\
\nmy $old_c = \"TemporaryChain\";\n\nwhile (<$INFI\
LE>) \n{\n	$line=$_;\n\n\n	if ($line =~ /^SEQRES/)\
\n{\n\n		@field=/(\\S*)\\s*/g;\n\n		$c= substr($_,\
11,1);\n\n		\n		$l=$#field;\n		for ($a=4; $a<$#fie\
ld ;)\n{\n			if (!$onelett{$molecule_type{$c}}->{$\
field[$a]})\n{\n				splice @field, $a, 1;\n}\n			e\
lse \n{\n				$a++;\n}\n}\n	\n		if ( $c ne $in_chai\
n)\n{\n			$pdb_chain_list[$n_pdb_chains]=$c;\n			$\
pdb_chain_len [$n_pdb_chains]=$len;\n			$in_chain=\
$c;\n			$n_pdb_chains++;\n}\n	\n		for ( $a=4; $a<$\
#field;$a++)\n{\n			@{$complete_seq{$c}}->[$comple\
te_seq_len{$c}++]=$field[$a];\n}\n}\n    elsif ( $\
line=~/^ATOM/ || ($line=~/^HETATM/ && &is_aa(subst\
r($line,17,3),substr($line,21,1)) && !$no_hetatm))\
\n{\n\n	 \n    $RAW_AT_ID=$AT_ID=substr($line,12,4\
);\n	$RES_ID=&is_aa(substr($line,17,3),substr($lin\
e,21,1));\n	$CHAIN=substr($line,21,1);\n\n    $RES\
_NO=substr($line,22,4);\n	$HOM_CODE=substr ($line,\
 26, 1);\n	$TEMP=substr($line,60,6);\n	\n	$TEMP=~s\
/\\s//g;\n        $AT_ID=~s/\\s//g;\n	$RES_ID=~s/\\
\s//g;\n        $RES_NO=~s/\\s//g;\n		\n	if ( $HOM\
_CODE ne $MAIN_HOM_CODE){next;}\n	elsif ( $already\
_read2{$CHAIN}{$RES_ID}{$AT_ID}{$RES_NO}){next;}\n\
	else{$already_read2{$CHAIN}{$RES_ID}{$AT_ID}{$RES\
_NO}=1;}\n	\n	\n	if ($coor_set && $numbering eq \"\
file\" && $residue_index ne $RES_NO)\n{\n	    \n	 \
   if ( $RES_NO<=$start){$real_start{$CHAIN}++;}\n\
	    if ( $RES_NO<=$end){$real_end{$CHAIN}++;}\n}\\
n	elsif ($numbering eq \"absolute\")\n{\n	    $rea\
l_start{$CHAIN}=$start;\n	    $real_end{$CHAIN}=$e\
nd;\n}\n\n        $KEY=\"ALL\";\n        if ( $CHA\
IN ne $in_atom_chain)\n{\n	    \n	  $pdb_atom_chai\
n_list[$n_pdb_atom_chains]=$c;\n	  $pdb_atom_chain\
_len [$n_pdb_atom_chains]=$len;\n	  $in_atom_chain\
=$c;\n	  $n_pdb_atom_chains++;\n}\n	\n	if ( $resid\
ue_index ne $RES_NO)\n{\n	     $residue_index = $R\
ES_NO;\n	     @{$atom_seq{$CHAIN}}->[$atom_seq_len\
{$CHAIN}++]=$RES_ID;;\n}\n}\n\n}\nclose ($INFILE);\
\n\n\n\n\n\n\n$INFILE=vfopen ($structure_file, \"r\
\");\nforeach $c (@c_chain)\n{\n\n	if    ( $seq_fi\
eld eq \"SEQRES\"){@pdb_seq=@{$complete_seq{$c}};}\
\n	elsif ( $seq_field eq \"ATOM\")  {@pdb_seq=@{$a\
tom_seq{$c}};}\n	\n\n	$full_length=$l=$#pdb_seq+1;\
\n		\n	if ( $real_end{$c}==\"*\"){$real_end{$c}=$f\
ull_length;}\n	if ( $coor_set)\n{	   \n\n	   if ( \
$real_end{$c} < $l){splice @pdb_seq, $real_end{$c}\
, $l;}\n	   if ( $real_start{$c} < $l){splice @pdb\
_seq, 0, $real_start{$c}-1;}	  	   \n	   $l=$#pdb_\
seq;\n}\n\n	elsif ( $delete_set)\n{\n	   splice @p\
db_seq, $delete_start, $delete_end-$delete_start+1\
;\n	   $l=$#pdb_seq;\n}\n	\n	$new_fasta_name=\"$pd\
b_id$c\";\n	if ( $coor_set)\n{\n	   if ( $n_pdb_ch\
ains==0){$new_fasta_name=\"$new_fasta_name$c\";}\n\
	   $new_fasta_name= $new_fasta_name.\"\\_$start\\\
_$end\";\n}\n	   \n	if ( $MODE eq \"pdb\")\n{\n	  \
 $nl=1;\n	   $n=0;\n	   \n	   foreach $res ( @pdb_\
seq)\n		{\n		if ( !$n)\n		{\n		\n		 printf \"SEQRE\
S %3d %1s %4d  \", $nl,$c, $l;\n		 $nl++;\n	}\n	  \
   $res=~s/\\s//g;\n	     \n	     if ($code==1){ p\
rintf \"%3s \",$onelett{$molecule_type{$c}}->{$res\
};}\n	     elsif  ($code==3){ printf \"%3s \",$res\
};\n	     \n	     $n++;		  \n	     if ( $n==13){$n\
=0;print \"\\n\";}\n}\n	  if ( $n!=0){print \"\\n\\
"; $n=0;}\n}\n	elsif ( $MODE eq \"simple\")\n{\n	 \
 print \"# SIMPLE_PDB_FORMAT\\n\";\n	  if ( $new_f\
asta_name eq \" \"){$new_fasta_name=\"dummy_name\"\
;}\n	  print \">$new_fasta_name\\n\";\n\n	  foreac\
h $res ( @pdb_seq)\n{\n	      print \"$onelett{$mo\
lecule_type{$c}}->{$res}\";\n}\n	  print \"\\n\";\\
n}\n	elsif ( $MODE eq \"fasta\")\n{\n	  $n=0;\n	  \
print \">$new_fasta_name\\n\";\n	  \n	  foreach $r\
es ( @pdb_seq)\n{\n\n	      print \"$onelett{$mole\
cule_type{$c}}->{$res}\";\n              $n++;\n	 \
     if ( $n==60){print \"\\n\"; $n=0;}\n}\n	  pri\
nt \"\\n\"; \n}\n}\n\nif ( $MODE eq \"fasta\")\n{\\
n     &myexit($EXIT_SUCCESS);\n  \n}\n\n  \n  $cha\
rcount=0;\n  $inchain=\"BEGIN\";\n  $n=0;\n  while\
 (<$INFILE>) \n{\n    $line=$_;\n     \n    if ($l\
ine =~/^ATOM/  ||  ($line=~/^HETATM/))\n{\n	$line_\
header=\"UNKNWN\";\n	$RES_ID=substr($line,17,3);\n\
	$chain = substr($line,21,1);\n\n	if ($line =~/^AT\
OM/)\n{\n	    $line_header=\"ATOM\";\n	    $RES_ID\
=(&is_aa($RES_ID,$chain))?&is_aa($RES_ID,$chain):$\
RES_ID;\n}\n	elsif ($line=~/^HETATM/ && ($ligand_l\
ist {$RES_ID} ||$ligand_list {'ALL'} || !&is_aa($R\
ES_ID,$chain)))\n{\n	    $line_header=\"HETATM\";\\
n}\n	elsif ($line=~/^HETATM/ && (&is_aa($RES_ID,$c\
hain) && !$no_hetatm))\n{\n	    $line_header=\"ATO\
M\";\n	    $RES_ID=&is_aa($RES_ID,$chain);\n}\n	el\
se\n{\n	    next;\n}\n\n	\n\n	$X=substr($line,30,8\
);     \n	$Y=substr($line,38,8);\n	$Z=substr($line\
,46,8);\n	$TEMP=substr($line,60,6);\n	\n	$RAW_AT_I\
D=$AT_ID=substr($line,12,4);\n	$CHAIN=substr($line\
,21,1);\n	$RES_NO=substr($line,22,4);\n	$HOM_CODE=\
substr ($line, 26, 1);\n	\n	$X=~s/\\s//g;\n	$Y=~s/\
\\s//g;\n	$Z=~s/\\s//g;\n	$TEMP=~s/\\s//g;\n	\n	$A\
T_ID=~s/\\s//g;\n	$RES_ID=~s/\\s//g;\n	$RES_NO=~s/\
\\s//g;\n\n	\n	if ( $HOM_CODE ne $MAIN_HOM_CODE){n\
ext;}\n	elsif ( $already_read{$CHAIN}{$RES_ID}{$AT\
_ID}{$RES_NO}){next;}\n	else{$already_read{$CHAIN}\
{$RES_ID}{$AT_ID}{$RES_NO}=1;}\n	\n	$KEY=\"ALL\";\\
n\n      	if ( $RES_NO ==0){$start_at_zero=1;}\n\n\
	$RES_NO+=$start_at_zero;    \n	\n	if ( $current_c\
hain ne $CHAIN)\n{\n	    $current_chain=$CHAIN;\n	\
    $pos=$current_residue=0;\n	    $offset=($coor_\
set)?($real_start{$CHAIN}-1):0;\n	    if    ( $seq\
_field eq \"SEQRES\"){@ref_seq=@{$complete_seq{$CH\
AIN}};}\n	    elsif ( $seq_field eq \"ATOM\")  {@r\
ef_seq=@{$atom_seq{$CHAIN}};}\n}\n	\n	if ($current\
_residue != $RES_NO)\n{\n	    $current_residue=$RE\
S_NO;\n	    if    ( $seq_field eq \"SEQRES\"){$pos\
=$current_residue;}\n	    elsif ( $seq_field eq \"\
ATOM\"){$pos++;}\n}\n	\n	\n	if ($n_atom==0 || $ato\
m_list{$AT_ID}==1 || $atom_list{$KEY}==1)\n{ 	\n	 \
   \n	    $do_it=(!@c_chain || $hc_chain{$CHAIN} |\
|$hc_chain{'LIGAND'} );\n	    \n	    $do_it= ($do_\
it==1) && ($coor_set==0 ||($pos>=$real_start{$CHAI\
N} && $pos<=$real_end{$CHAIN}));\n	    $do_it= ($d\
o_it==1) && ($delete_set==0 || $pos<$delete_start \
||$pos>$delete_end );\n	    if ($ligand==0 && $lin\
e_header eq \"HETATM\" ){$do_it=0;}\n	    if ($lig\
and_only==1 && $line_header eq \"ATOM\" ){$do_it=0\
;}\n	    if ($ligand==1 && $line_header eq \"HETAT\
M\" && $ligand_list{$RES_ID}==0 && $ligand_list{\"\
ALL\"}==0){$do_it=0;} \n	    \n	    \n	    if ( $d\
o_it)\n{\n		$n++;\n		$out_pos=$pos;\n		\n	       i\
f ( $delete_set)\n{\n		  if ( $out_pos< $delete_st\
art){;}\n		  else {$offset=$delete_end-$delete_sta\
rt;}\n}       \n	       \n	       if ( $numbering_\
out eq \"new\"){$out_pos-=$offset;}\n	       elsif\
 ( $numbering_out eq \"old\"){$out_pos=$RES_NO;}\n\
	       \n       \n	       \n	       if ( $code==1\
){$RES_ID=$onelett{$molecule_type{$c}}->{$RES_ID};\
}\n	    \n	       if ($unfold)\n{\n		   $unfolded_\
x+=5;\n		   $X=$unfolded_x;\n		   $Y=0;\n		   $Z=0\
;\n		   $float=1;\n}\n	       else\n{\n		   $float\
=3;\n}\n\n	       if ( $MODE eq \"pdb\")\n{\n		   \
printf \"%-6s%5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f\
  1.00 %5.2f\\n\",$line_header, $n, $RAW_AT_ID,$RE\
S_ID,$CHAIN,$out_pos, $X, $Y, $Z,$TEMP;		  \n}\n	 \
      elsif ( $MODE eq \"simple\")\n{\n		    if ( \
$RES_ID eq \"\"){$RES_ID=\"X\";}\n		  printf \"%-6\
s %5s %s %2s %4d    %8.3f %8.3f %8.3f\\n\",$line_h\
eader, $AT_ID, $RES_ID,($CHAIN eq\"\" || $CHAIN eq\
 \" \")?\"A\":$CHAIN,$out_pos, $X, $Y, $Z,$TEMP;\n\
}\n\n}\n}\n}\n}\nprint \"\\n\";\nclose($INFILE);\n\
\n\nif ( $error ne \"\") \n{$error=$error.\"\\nDia\
gnostic:    SEQRES and the residues in ATOM are pr\
obably Incompatible\\n\";\n    $error=$error.  \"R\
ecomendation: Rerun with '-fix 1' in order to igno\
re the SEQRES sequences\\n\";\n}\nif (!$nodiagnost\
ic){print STDERR $error;}\n&myexit ( $EXIT_SUCCESS\
);\n\nsub remote_is_pdb_name \n{\n    my $in=@_[0]\
;\n    my $netaddress=@_[1];\n    my $ref_file, $p\
db;\n    my $value;\n\n    if ( $in=~/[^\\w\\d\\:\\
\_]/){return 0;}\n    \n    \n    $ref_file=\"$cac\
he/pdb_entry_type.txt\";\n    \n    if ( !-e $ref_\
file || (-M $ref_file)>2 || -z $ref_file)\n{\n	&ur\
l2file(\"ftp://ftp.wwpdb.org/pub/pdb/derived_data/\
pdb_entry_type.txt\", $ref_file);\n}\n \n    $pdb=\
substr ($in,0, 4);\n    \n    \n    $value=`grep -\
c $pdb $ref_file`;\n    return $value;\n}\n      \\
nsub is_pdb_file\n{\n    my @arg=@_;\n\n    if ( !\
-e $arg[0]){return 0;}\n    \n    $F=vfopen ($arg[\
0], \"r\");\n    while ( <$F>)\n{\n	if (/^HEADER/)\
\n{\n	    close $F;\n	    return 1;\n}\n	elsif ( /\
^SEQRES/)\n{\n	    close $F;\n	    return 1;\n}\n	\
elsif ( /^ATOM/)\n{\n	    close $F;\n	    return 1\
;\n}\n}\n    return 0;\n}\nsub get_pdb_id\n{\n    \
my $header_file=@_[0];\n    my $id;\n    my $F= ne\
w FileHandle;\n    \n    \n    $F=vfopen (\"$heade\
r_file\", \"r\");\n\n    while ( <$F>)\n      {\n	\
if ( /HEADER/)\n	  {\n	    if ($debug){print \"$_\\
";}\n	    $id=substr($_,62,4 );\n	    return $id;\\
n	  }\n      }\n    close ($F);\n    \n    return \
\"\";\n}\n\nsub get_ligand_list\n{\n    my $pdb_fi\
le=@_[0];\n    my $chain;\n    my $ligand;\n    my\
 %complete_ligand_list;\n    \n\n    $F=vfopen ($p\
db_file, \"r\");\n    while ( <$F>)\n{\n	if ( /^HE\
TATM/)\n{\n	    $line=$_;\n	    $chain=substr($lin\
e,21,1);\n	    $ligand=substr($line,17,3);\n	    \\
n	    if (!$complete_ligand_list{$chain}{$ligand})\
\n{\n		\n		$complete_ligand_list{\"result\"}.=\"CH\
AIN $chain LIGAND $ligand\\n\";\n		$complete_ligan\
d_list{$chain}{$ligand}=1;\n}\n}\n}\n    close ($F\
);\n    return %complete_ligand_list;\n}\n\nsub ge\
t_chain_list \n{\n    my $header_file;\n    my @ch\
ain_list;\n    my @list;\n    my $n_chains;\n    m\
y %chain_hasch;\n    my $pdb_file=@_[0];\n    my $\
c;\n    my %hasch;\n    my $chain;\n  \n    \n    \
$F=vfopen ($pdb_file, \"r\");\n    while ( <$F>)\n\
{\n\n\n	if (/SEQRES\\s+\\d+\\s+(\\S+)/)\n	  {\n	  \
  $chain = substr($_,11,1);$chain=~s/\\s//g;if ( $\
chain eq \"\"){$chain=\" \";}\n	    if (!$hasch{$c\
hain}){$hasch{$chain}=1;push @chain_list, $chain;}\
\n	  }\n	if (/^ATOM/ || /^HETATM/)\n	  {\n	    $ch\
ain = substr($_,21,1); $chain=~s/\\s//g;if ( $chai\
n eq \"\"){$chain=\" \";}\n	    if (!$hasch{$chain\
}){$hasch{$chain}=1;push @chain_list, $chain;}\n	 \
 }\n      }\n\n\nclose ($F);\nif (!@chain_list)\n \
 {\n    @chain_list=(\"A\");\n  }\n\n\nreturn @cha\
in_list;\n}\n\nsub token_is_in_list\n{\n\n    my @\
list=@_;\n    my $a;\n    \n    for ($a=1; $a<=$#l\
ist; $a++)\n{\n	if ( $list[$a] eq $list[0]){return\
 $a;}\n}\n}\n\nsub pdb_name2name_and_chain \n{\n  \
  my $pdb_file=@_[0];\n    my $pdb_file_in;\n    m\
y @array;\n    my $chain;\n    my $c;\n\n    $pdb_\
file_in=$pdb_file;\n\n    $pdb_file=~/^(.{4})/;$pd\
b_id=$1;\n    @array=($pdb_file=~/([\\w])/g);\n  \\
n  \n    $chain=uc ($array[4]);\n    $chain=($chai\
n eq \"\")?\"FIRST\":$chain;\n    \n    return ( $\
pdb_id, $chain);\n\n    if ( $#array==3){return ($\
pdb_id, \"FIRST\");}\n    elsif ( $#array<4){ retu\
rn ($pdb_id, \"\");}\n    else {return ( $pdb_id, \
$chain);}\n      \n    \n    \n}\nsub get_main_hom\
_code \n{\n    my $pdb_file=@_[0];\n    my %hom, $\
n, $best, $best_h;\n    open (F, $pdb_file);\n    \
while (<F>)\n{\n	if ( $_=~/^ATOM/)\n{\n	    $h=sub\
str ($_,26, 1);\n	    $n=++$hom{$h};\n	    if ($n>\
$best)\n{\n		$best=$n;\n		$best_h=$h;\n}\n}\n}\n  \
  close (F);\n    return $best_h;\n}\n\n\nsub get_\
pdb_file \n{\n    my ($pdb_file_in)=(@_);\n    my \
$result;\n    my @letter;\n    my @chain;\n    my \
$v;\n    my $pdb_file=$pdb_file_in;\n\n    $pdb_fi\
le=($pdb_file_in=~/\\S+_S_(\\S+)/)?$1:$pdb_file_in\
;\n\n    if ($no_remote_pdb_dir==0)\n{\n	$no_remot\
e_pdb_dir=1;\n	$result=get_pdb_file3 ($pdb_file);\\
n	$no_remote_pdb_dir=0;\n	if ( $result){return $re\
sult;}\n	else\n{\n	    \n	    lc ($pdb_file);\n	  \
  $result=get_pdb_file3($pdb_file);\n	    return  \
$result;\n}\n}\n    else\n{\n	return get_pdb_file3\
 ($pdb_file);\n}\n    \n}\n\nsub get_pdb_file3 \n{\
\n    my $pdb_file_in=@_[0];\n    my $result;\n   \
 my @letter;\n    my @chain;\n    my $lcfile;\n   \
 my $ucfile;\n    my $pdb_file=$pdb_file_in;\n    \
\n    $lcfile=lc $pdb_file;\n    $ucfile=uc $pdb_f\
ile;\n\n    if ( ($result=get_pdb_file2 ($pdb_file\
))){return $result;}\n    \n\n    if ($lcfile ne $\
pdb_file && ($result=get_pdb_file2 ($lcfile))){ret\
urn $result;}\n    if ($ucfile ne $pdb_file && ($r\
esult=get_pdb_file2 ($ucfile))){return $result;}\n\
    \n   \n    \n    return \"\";\n}\nsub get_pdb_\
file2\n{\n    my $pdb_file=@_[0];\n    my $return_\
value;\n    \n    $return_value=\"\";\n    \n    i\
f ( ($result=get_pdb_file1 ($pdb_file))){$return_v\
alue=$result;}\n    elsif ( !($pdb_file=~/\\.pdb/)\
 && !($pdb_file=~/\\.PDB/))\n{\n	if ( ($result=get\
_pdb_file1 (\"$pdb_file.pdb\"))){$return_value=$re\
sult;}\n	elsif ( ($result=get_pdb_file1 (\"$pdb_fi\
le.PDB\"))){$return_value=$result;}\n\n	elsif ( ($\
result=get_pdb_file1 (\"pdb$pdb_file.pdb\"))){$ret\
urn_value=$result;}	\n	elsif ( ($result=get_pdb_fi\
le1 (\"pdb$pdb_file.PDB\"))){$return_value=$result\
;}\n	elsif ( ($result=get_pdb_file1 (\"PDB$pdb_fil\
e.PDB\"))){$return_value=$result;}\n	elsif ( ($res\
ult=get_pdb_file1 (\"PDB$pdb_file.pdb\"))){$return\
_value=$result;}\n	\n	\n	elsif ( ($result=get_pdb_\
file1 (\"$pdb_file.ent\"))){$return_value=$result;\
}\n	elsif ( ($result=get_pdb_file1 (\"pdb$pdb_file\
.ent\"))){$return_value=$result;}\n	elsif ( ($resu\
lt=get_pdb_file1 (\"PDB$pdb_file.ent\"))){$return_\
value=$result;}\n\n	elsif ( ($result=get_pdb_file1\
 (\"$pdb_file.ENT\"))){$return_value=$result;}\n	e\
lsif ( ($result=get_pdb_file1 (\"pdb$pdb_file.ENT\\
"))){$return_value=$result;}\n	elsif ( ($result=ge\
t_pdb_file1 (\"PDB$pdb_file.ENT\"))){$return_value\
=$result;}\n	\n	\n	\n}\n    return $return_value;\\
n}\n    \nsub get_pdb_file1\n{\n    my ($pdb_file)\
=(@_);\n    my $return_value;\n    \n\n    $return\
_value=\"\";\n    if ( ($result=get_pdb_file0 ($pd\
b_file))){$return_value=$result;}\n    elsif ( ($r\
esult=get_pdb_file0 (\"$pdb_file.Z\"))){$return_va\
lue=$result;}\n    elsif ( ($result=get_pdb_file0 \
(\"$pdb_file.gz\"))){$return_value=$result;}\n    \
elsif ( ($result=get_pdb_file0 (\"$pdb_file.GZ\"))\
){$return_value=$result;}\n    return $return_valu\
e;\n}\nsub get_pdb_file0 \n{ \n    my ($pdb_file, \
$attempt)=(@_);\n    my $pdb_file=@_[0];\n    my $\
tmp_pdb_file;    \n    my $return_value;\n\n    if\
 ( !$attempt){$attempt=1;}\n    \n    $local_pdb_f\
ile=\"$pdb_file\";\n    if ( $local_pdb_file eq \"\
\")\n{\n	$tmp_pdb_file=vtmpnam();\n	open F, \">$tm\
p_pdb_file\";\n	\n	while (<STDIN>){print F \"$_\";\
}\n	close (F);\n	\n	if (-e $tmp_pdb_file && &is_pd\
b_file ( $local_pdb_file))\n{return $tmp_pdb_file;\
}\n}\n\n    $local_pdb_file=\"$pdb_file\";\n    &d\
ebug_print (\"\\nTry access local file: $local_pdb\
_file\");\n    \n    $local_pdb_file=&check_pdb_fi\
le4compression ($local_pdb_file);\n    if ( -e $lo\
cal_pdb_file && (&is_pdb_file ($local_pdb_file) ||\
 $force_pdb))\n{\n	&debug_print ( \"\\n\\tIs in Cu\
rrent Dir\");\n	$tmp_pdb_file=vtmpnam();\n	`cp $lo\
cal_pdb_file $tmp_pdb_file`;\n	return $tmp_pdb_fil\
e;\n}\n    else\n{\n	&debug_print (\"\\n\\tFile No\
t in Current Dir\");\n}\n\n    if ($pdb_file=~/^pd\
b/||$pdb_file=~/^PDB/){$pdb_div=substr ($pdb_file,\
 4, 2);}\n    else\n{\n	  $pdb_div=substr ($pdb_fi\
le, 1, 2);\n}\n    $local_pdb_file=\"$pdb_dir/$pdb\
_div/$pdb_file\";\n    $local_pdb_file=&check_pdb_\
file4compression ( $local_pdb_file);\n    &debug_p\
rint (\"\\nTry access file From PDB_DIR: $local_pd\
b_file\");\n    if ($pdb_dir && -e $local_pdb_file\
 && &is_pdb_file ($local_pdb_file))\n{\n	&debug_pr\
int ( \"\\n\\tIs in Local PDB DIR\");\n	$tmp_pdb_f\
ile=vtmpnam();\n	`cp $local_pdb_file $tmp_pdb_file\
`;\n	return $tmp_pdb_file;\n}\n\n    $local_pdb_fi\
le=\"$pdb_dir/$pdb_file\";\n    $local_pdb_file=&c\
heck_pdb_file4compression ( $local_pdb_file);\n   \
 &debug_print (\"\\nTry access file From PDB_DIR: \
local_pdb_file\");\n    if ($pdb_dir && -e $local_\
pdb_file && &is_pdb_file ($local_pdb_file))\n{\n	&\
debug_print ( \"\\n\\tIs in Local PDB DIR\");\n	$t\
mp_pdb_file=vtmpnam();\n	`cp $local_pdb_file $tmp_\
pdb_file`;\n	return $tmp_pdb_file;\n}\n\n    $loca\
l_pdb_file=\"$pdb_dir$pdb_file\";\n    $local_pdb_\
file=&check_pdb_file4compression ( $local_pdb_file\
);\n    &debug_print (\"\\nTry access file From PD\
B_DIR: $local_pdb_file\");\n    if ($pdb_dir && -e\
 $local_pdb_file && &is_pdb_file ($local_pdb_file)\
)\n{\n	&debug_print ( \"\\n\\tIs in Local PDB DIR\\
");\n	$tmp_pdb_file=vtmpnam();\n	`cp $local_pdb_fi\
le $tmp_pdb_file`;\n	return $tmp_pdb_file;\n}\n   \
 else\n{&debug_print ( \"\\n\\tNot In Local Pdb Di\
r\");}\n\n    if ($cache ne \"NO\" && $cache ne \"\
no\")\n{\n\n	$local_pdb_file=\"$cache/$pdb_file\";\
\n	$local_pdb_file=&check_pdb_file4compression ( $\
local_pdb_file);\n	&debug_print(\"\\nTry access fi\
le From Cache: $local_pdb_file\");\n	if (-e $local\
_pdb_file && &is_pdb_file ($local_pdb_file))\n{\n	\
    &debug_print ( \"\\n\\tIs in T-Coffee Cache\")\
;\n	    $tmp_pdb_file=vtmpnam();\n	    `cp $local_\
pdb_file $tmp_pdb_file`;\n	    return $tmp_pdb_fil\
e;\n}\n	else{&debug_print ( \"\\n\\tNot in Cache D\
ir\");}\n}\n\nif (!$no_remote_pdb_dir) \n  {\n    \
my $value=remote_is_pdb_name ($pdb_file, $netaddre\
ss);\n    \n    my $return_value=\"\";\n    if ( &\
remote_is_pdb_name ($pdb_file, $netaddress)==1)\n \
     {\n	\n	&debug_print (\"\\n*******************\
************************************\\nTry Remote \
Access for $pdb_file\");\n	$tmp_pdb_file=vtmpnam()\
;\n	$netcommand=$netaddress;\n	$netcommand=~s/%%/$\
pdb_file/g;\n	&url2file(\"$netcommand\", \"$tmp_pd\
b_file.$netcompression\");\n	&debug_print(\"\\nREM\
OTE: $netcommand\\n\");\n	\n	$compressed_tmp_file_\
name=\"$tmp_pdb_file.$netcompression\";\n	\n	if ($\
netcompression && -B $compressed_tmp_file_name)\n	\
  {\n	    my $r;\n	    &debug_print (\"\\n\\tFile \
Found Remotely\");\n	    if (($r=safe_system ( \"$\
netcompression_pg $compressed_tmp_file_name\")!=$E\
XIT_SUCCESS) && $attempts<5)\n	      {\n		&debug_p\
rint (\"\\n\\tProper Download Failed Try again\");\
\n		unlink $compressed_tmp_file_name;\n		print \"\\
\nFailed to Download $compressed_tmp_file_name. Ne\
w Attempt $attempt/5\\n\";\n		return &get_pdb_file\
0($pdb_file, $attempt+1);\n	      }\n	    elsif ($\
r== $EXIT_SUCCESS)\n	      {\n		&debug_print (\"\\\
n\\tProper Download Succeeded \");\n		$return_valu\
e=$tmp_pdb_file;\n	      }\n	    else\n	      {\n	\
	&debug_print (\"\\n\\tProper Download Failed \");\
\n		&debug_print (\"\\nFile Not Found Remotely\");\
\n		unlink $compressed_tmp_file_name;\n	      }\n	\
  }\n	else\n	  {\n\n	    &debug_print (\"\\nFile N\
ot Found Remotely\");\n	    unlink $compressed_tmp\
_file_name;\n	  }\n	#Update cache if required\n	if\
 ($cache ne \"no\" && $cache ne \"update\" && -e $\
return_value)\n	  {\n	    `cp $return_value $cache\
/$pdb_file.pdb`;\n	    #`t_coffee -other_pg clean_\
cache.pl -file $pdb_file.pdb -dir $cache`;\n	  }\n\
      }\n    &debug_print (\"\\nRemote Download Fi\
nished\");\n    return $return_value;\n  }\nreturn\
 \"\";\n}\n\nsub check_pdb_file4compression \n{\n \
   my $file=@_[0];\n    my $tmp;\n    my $r;\n    \
\n    $tmp=&vtmpnam();\n    if (-e $tmp){unlink $t\
mp;}\n    \n    $file=~s/\\/\\//\\//g;\n    if    \
(-B $file && ($file=~/\\.Z/)) {`cp $file $tmp.Z`;`\
rm $tmp`;`gunzip $tmp.Z $SILENT`;$r=$tmp;}\n    el\
sif (-B $file && ($file=~/\\.gz/)){`cp $file $tmp.\
gz`;`gunzip $tmp.gz $SILENT`;return $r=$tmp;}\n   \
 elsif (-B $file ){`cp $file $tmp.gz`;`gunzip $tmp\
.gz $SILENT`;$r=$tmp;}\n    elsif ( -e $file ) {$r\
= $file;}\n    elsif ( -e \"$file.gz\" ){ `cp $fil\
e.gz $tmp.gz`;`gunzip     $tmp.gz $SILENT`;$r=$tmp\
;}    \n    elsif ( -e \"$file.Z\") {`cp $file.Z  \
$tmp.Z`; `gunzip $tmp.Z $SILENT`;$r=$tmp;}\n    el\
se  {$r= $file;}\n\n    if ( -e \"$tmp.Z\"){unlink\
 \"$tmp.Z\";}\n    if ( -e \"$tmp.gz\"){unlink \"$\
tmp.gz\";}\n    \n    return $r;\n    \n}\n\n\n\n\\
n\n    \n\n\n\n\n\n\n\nsub vfopen \n{\n    my $fil\
e=@_[0];\n    my $mode=@_[1];\n    my $tmp;\n    m\
y $F = new FileHandle;\n    \n    \n    $tmp=$file\
;\n	\n    \n    if ( $mode eq \"r\" && !-e $file){\
 myexit(flush_error (\"Cannot open file $file\"));\
}\n    elsif ($mode eq \"w\"){$tmp=\">$file\";}\n \
   elsif ($mode eq \"a\"){$tmp=\">>$file\";}\n    \
\n    \n    open ($F,$tmp);\n    return $F;\n}\nsu\
b debug_print\n{\n    my $message =@_[0];\n    if \
($debug){print STDERR \"NO_REMOTE_PDB_DIR: $no_rem\
ote_pdb_dir - $message [DEBUG:extract_from_pdb]\";\
}\n    return;\n}\nsub is_aa \n{\n    my ($aa, $ch\
ain) =@_;\n\n    my $one;\n    my $trhee;\n    \n \
   if ( $onelett{$molecule_type{$chain}}->{$aa} eq\
 'X' || !$onelett{$molecule_type{$chain}}->{$aa} )\
{return '';}\n    else\n      {\n	$one=$onelett{$m\
olecule_type{$chain}}->{$aa};\n\n	$three=$threelet\
t{$molecule_type{$chain}}->{$one};\n	\n\n	return $\
three;\n      }\n  }\n\n\n\n\n\nsub url2file\n{\n \
   my ($address, $out, $wget_arg, $curl_arg)=(@_);\
\n    my ($pg, $flag, $r, $arg, $count);\n    \n  \
  if (!$CONFIGURATION){&check_configuration (\"wge\
t\", \"INTERNET\", \"gzip\");$CONFIGURATION=1;}\n \
   \n    if (&pg_is_installed (\"wget\"))   {$pg=\\
"wget\"; $flag=\"-O\";$arg=$wget_arg;}\n    elsif \
(&pg_is_installed (\"curl\")){$pg=\"curl\"; $flag=\
\"-o\";$arg=$curl_arg;}\n    return safe_system (\\
"$pg $flag$out $address >/dev/null 2>/dev/null\");\
\n\n}\n\n\n\n\nsub pdbfile2chaintype\n  {\n    my \
$file=@_[0];\n    my %ct;\n    my $F;\n    \n    $\
F=vfopen ($file, \"r\");\n    while (<$F>)\n      \
{\n	my $line=$_;\n	if ($line =~/^ATOM/)\n	  {\n	  \
  my $C=substr($line,21,1);\n	    if (!$ct{$C})\n	\
      {	\n		my $r=substr($line,17,3);\n		$r=~s/\\s\
+//;\n		if (length ($r)==1){$ct{$C}=\"R\";}\n		els\
if (length ($r)==2){$ct{$C}=\"D\";}\n		elsif (leng\
th ($r)==3){$ct{$C}=\"P\";}\n		else \n		  {\n		   \
 myexit(flush_error(\"ERROR: Could not read RES_ID\
 field in file $file\"));\n		  }\n	      }\n	  }\n\
      }\n    close ($F);\n    return %ct;\n  }\n  \
 \n    \n\n\n\nsub fill_threelett_RNA\n{\n\n	my %t\
hreelett=(\n	'A', '  A',\n	'T', '  T',\n	'U', '  U\
',\n	'C', '  C',\n	'G', '  G',\n	'I', '  I', #Inos\
ine\n	);\n	\n	return %threelett;\n\n}\n\n\nsub fil\
l_onelett_RNA\n{\n	my   %onelett=(\n	'  A' => 'A',\
\n	'  T' => 'T',\n	'  U' => 'U',\n	'  C' => 'C',\n\
	'  G' => 'G',\n	'CSL' => 'X',\n	'UMS' => 'X',\n	'\
  I' => 'I',\n	'A' => 'A',\n	'T' => 'T',\n	'U' => \
'U',\n	'C' => 'C',\n	'G' => 'G',\n	'I' => 'I',\n	)\
;\n\n	return %onelett;\n\n}\n\n\nsub fill_onelett_\
DNA\n{\n	my   %onelett=(\n	' DA', 'A',\n	' DT', 'T\
',\n	' DC', 'C',\n	' DG', 'G',\n	'DA', 'A',\n	'DT'\
, 'T',\n	'DC', 'C',\n	'DG', 'G',\n	);\n\n	return %\
onelett;\n\n}\n\nsub fill_threelett_DNA\n{\n\n	my \
%threelett=(\n	'A', ' DA',\n	'T', ' DT',\n	'C', ' \
DC',\n	'G', ' DG',\n	);\n\n	return %threelett;\n\n\
}\n\n\n\n\nsub fill_threelett_prot\n{  \n  my %thr\
eelett;\n\n  %threelett=(\n'A', 'ALA',\n'C', 'CYS'\
,\n'D', 'ASP',\n'E', 'GLU',\n'F', 'PHE',\n'G', 'GL\
Y',\n'H', 'HIS',\n'I', 'ILE',\n'K', 'LYS',\n'L', '\
LEU',\n'N', 'ASN',\n'M', 'MET',\n'P', 'PRO',\n'Q',\
 'GLN',\n'R', 'ARG',\n'S', 'SER',\n'T', 'THR',\n'V\
', 'VAL',\n'W', 'TRP',\n'Y', 'TYR',\n);\n\nreturn \
%threelett;\n\n\n}\n\nsub fill_onelett_prot\n{\n  \
  my %onelett;\n    \n    %onelett=(\n\n'10A', 'X'\
,\n'11O', 'X',\n'12A', 'X',\n'13P', 'X',\n'13R', '\
X',\n'13S', 'X',\n'14W', 'X',\n'15P', 'X',\n'16A',\
 'X',\n'16G', 'X',\n'1AN', 'X',\n'1AP', 'X',\n'1AR\
', 'X',\n'1BH', 'X',\n'1BO', 'X',\n'1C5', 'X',\n'1\
CU', 'X',\n'1DA', 'X',\n'1GL', 'X',\n'1GN', 'X',\n\
'1IN', 'X',\n'1LU', 'L',\n'1MA', 'X',\n'1MC', 'X',\
\n'1MG', 'X',\n'1MZ', 'X',\n'1NA', 'X',\n'1NB', 'X\
',\n'1NI', 'X',\n'1PA', 'A',\n'1PC', 'X',\n'1PE', \
'X',\n'1PG', 'X',\n'1PI', 'A',\n'1PM', 'X',\n'1PN'\
, 'X',\n'1PU', 'X',\n'1PY', 'X',\n'1UN', 'X',\n'24\
T', 'X',\n'25T', 'X',\n'26P', 'X',\n'2AB', 'X',\n'\
2AM', 'X',\n'2AN', 'X',\n'2AP', 'X',\n'2AR', 'X',\\
n'2AS', 'D',\n'2BL', 'X',\n'2BM', 'X',\n'2CP', 'X'\
,\n'2DA', 'X',\n'2DG', 'X',\n'2DP', 'X',\n'2DT', '\
X',\n'2EP', 'X',\n'2EZ', 'X',\n'2FG', 'X',\n'2FL',\
 'X',\n'2FP', 'X',\n'2FU', 'X',\n'2GL', 'X',\n'2GP\
', 'X',\n'2HP', 'X',\n'2IB', 'X',\n'2IP', 'X',\n'2\
LU', 'L',\n'2MA', 'X',\n'2MD', 'X',\n'2ME', 'X',\n\
'2MG', 'X',\n'2ML', 'L',\n'2MO', 'X',\n'2MR', 'R',\
\n'2MU', 'X',\n'2MZ', 'X',\n'2NO', 'X',\n'2NP', 'X\
',\n'2OG', 'X',\n'2PA', 'X',\n'2PC', 'X',\n'2PE', \
'X',\n'2PG', 'X',\n'2PH', 'X',\n'2PI', 'X',\n'2PL'\
, 'X',\n'2PP', 'X',\n'2PU', 'X',\n'2SI', 'X',\n'2T\
B', 'X',\n'34C', 'X',\n'35G', 'X',\n'3AA', 'X',\n'\
3AD', 'X',\n'3AH', 'H',\n'3AN', 'X',\n'3AP', 'X',\\
n'3AT', 'X',\n'3BT', 'X',\n'3CH', 'X',\n'3CN', 'X'\
,\n'3CO', 'X',\n'3CP', 'X',\n'3DR', 'X',\n'3EP', '\
X',\n'3FM', 'X',\n'3GA', 'X',\n'3GP', 'X',\n'3HB',\
 'X',\n'3HC', 'X',\n'3HP', 'X',\n'3IB', 'X',\n'3ID\
', 'X',\n'3IN', 'X',\n'3MA', 'X',\n'3MB', 'X',\n'3\
MC', 'X',\n'3MD', 'D',\n'3MF', 'X',\n'3MP', 'X',\n\
'3MT', 'X',\n'3OL', 'X',\n'3PA', 'X',\n'3PG', 'X',\
\n'3PO', 'X',\n'3PP', 'X',\n'3PY', 'X',\n'49A', 'X\
',\n'4AB', 'X',\n'4AM', 'X',\n'4AN', 'X',\n'4AP', \
'X',\n'4BA', 'X',\n'4BT', 'X',\n'4CA', 'X',\n'4CO'\
, 'X',\n'4HP', 'X',\n'4IP', 'X',\n'4MO', 'X',\n'4M\
V', 'X',\n'4MZ', 'X',\n'4NC', 'X',\n'4NP', 'X',\n'\
4OX', 'X',\n'4PB', 'X',\n'4PN', 'X',\n'4PP', 'X',\\
n'4SC', 'X',\n'4SU', 'X',\n'4TB', 'X',\n'55C', 'X'\
,\n'5AD', 'X',\n'5AN', 'X',\n'5AT', 'X',\n'5CM', '\
X',\n'5GP', 'X',\n'5HP', 'E',\n'5HT', 'X',\n'5IT',\
 'X',\n'5IU', 'X',\n'5MB', 'X',\n'5MC', 'X',\n'5MD\
', 'X',\n'5MP', 'X',\n'5MU', 'X',\n'5NC', 'X',\n'5\
OB', 'X',\n'5PA', 'X',\n'5PV', 'X',\n'6AB', 'X',\n\
'6CT', 'X',\n'6HA', 'X',\n'6HC', 'X',\n'6HG', 'X',\
\n'6HT', 'X',\n'6IN', 'X',\n'6MO', 'X',\n'6MP', 'X\
',\n'6PG', 'X',\n'6WO', 'X',\n'70U', 'X',\n'7DG', \
'X',\n'7HP', 'X',\n'7I2', 'X',\n'7MG', 'X',\n'7MQ'\
, 'X',\n'7NI', 'X',\n'87Y', 'X',\n'8AD', 'X',\n'8B\
R', 'X',\n'8IG', 'X',\n'8IN', 'X',\n'8OG', 'X',\n'\
95A', 'X',\n'9AD', 'X',\n'9AM', 'X',\n'9AP', 'X',\\
n'9DG', 'X',\n'9DI', 'X',\n'9HX', 'X',\n'9OH', 'X'\
,\n'9TA', 'X',\n'A12', 'X',\n'A15', 'X',\n'A23', '\
X',\n'A24', 'X',\n'A26', 'X',\n'A2G', 'X',\n'A2P',\
 'X',\n'A32', 'X',\n'A3P', 'X',\n'A4P', 'X',\n'A5P\
', 'X',\n'A70', 'X',\n'A76', 'X',\n'A77', 'X',\n'A\
78', 'X',\n'A79', 'X',\n'A80', 'X',\n'A85', 'X',\n\
'A88', 'X',\n'A9A', 'X',\n'AA3', 'X',\n'AA4', 'X',\
\n'AA6', 'X',\n'AAA', 'X',\n'AAB', 'X',\n'AAC', 'X\
',\n'AAE', 'X',\n'AAG', 'R',\n'AAH', 'X',\n'AAM', \
'X',\n'AAN', 'X',\n'AAP', 'X',\n'AAR', 'R',\n'AAS'\
, 'X',\n'AAT', 'X',\n'ABA', 'X',\n'ABC', 'X',\n'AB\
D', 'X',\n'ABE', 'X',\n'ABH', 'X',\n'ABI', 'X',\n'\
ABK', 'X',\n'ABM', 'X',\n'ABN', 'X',\n'ABP', 'X',\\
n'ABR', 'X',\n'ABS', 'X',\n'ABU', 'X',\n'AC1', 'X'\
,\n'AC2', 'X',\n'ACA', 'X',\n'ACB', 'D',\n'ACC', '\
C',\n'ACD', 'X',\n'ACE', 'X',\n'ACH', 'X',\n'ACI',\
 'X',\n'ACL', 'R',\n'ACM', 'X',\n'ACN', 'X',\n'ACO\
', 'X',\n'ACP', 'X',\n'ACQ', 'X',\n'ACR', 'X',\n'A\
CS', 'X',\n'ACT', 'X',\n'ACV', 'V',\n'ACX', 'X',\n\
'ACY', 'X',\n'AD2', 'X',\n'AD3', 'X',\n'ADC', 'X',\
\n'ADD', 'X',\n'ADE', 'X',\n'ADH', 'X',\n'ADI', 'X\
',\n'ADM', 'X',\n'ADN', 'X',\n'ADP', 'X',\n'ADQ', \
'X',\n'ADR', 'X',\n'ADS', 'X',\n'ADT', 'X',\n'ADU'\
, 'X',\n'ADW', 'X',\n'ADX', 'X',\n'AE2', 'X',\n'AE\
A', 'X',\n'AEB', 'X',\n'AEI', 'D',\n'AEN', 'X',\n'\
AET', 'T',\n'AF1', 'X',\n'AF3', 'X',\n'AFA', 'D',\\
n'AFP', 'X',\n'AG7', 'X',\n'AGB', 'X',\n'AGF', 'X'\
,\n'AGL', 'X',\n'AGM', 'R',\n'AGN', 'X',\n'AGP', '\
X',\n'AGS', 'X',\n'AGU', 'X',\n'AH0', 'X',\n'AH1',\
 'X',\n'AHA', 'X',\n'AHB', 'D',\n'AHC', 'X',\n'AHF\
', 'X',\n'AHG', 'X',\n'AHH', 'X',\n'AHM', 'X',\n'A\
HO', 'X',\n'AHP', 'X',\n'AHS', 'X',\n'AHT', 'Y',\n\
'AHU', 'X',\n'AHX', 'X',\n'AI1', 'X',\n'AI2', 'X',\
\n'AIB', 'X',\n'AIC', 'X',\n'AIM', 'X',\n'AIP', 'X\
',\n'AIQ', 'X',\n'AIR', 'X',\n'AJ3', 'X',\n'AKB', \
'X',\n'AKG', 'X',\n'AKR', 'X',\n'AL1', 'X',\n'AL2'\
, 'X',\n'AL3', 'X',\n'AL4', 'X',\n'AL5', 'X',\n'AL\
6', 'X',\n'AL7', 'X',\n'AL8', 'X',\n'AL9', 'X',\n'\
ALA', 'A',\n'ALB', 'X',\n'ALC', 'X',\n'ALD', 'L',\\
n'ALE', 'X',\n'ALF', 'X',\n'ALG', 'X',\n'ALL', 'X'\
,\n'ALM', 'A',\n'ALN', 'A',\n'ALO', 'T',\n'ALP', '\
X',\n'ALQ', 'X',\n'ALR', 'X',\n'ALS', 'X',\n'ALT',\
 'A',\n'ALY', 'K',\n'ALZ', 'X',\n'AMA', 'X',\n'AMB\
', 'X',\n'AMC', 'X',\n'AMD', 'X',\n'AMG', 'X',\n'A\
MH', 'X',\n'AMI', 'X',\n'AML', 'X',\n'AMN', 'X',\n\
'AMO', 'X',\n'AMP', 'X',\n'AMQ', 'X',\n'AMR', 'X',\
\n'AMS', 'X',\n'AMT', 'X',\n'AMU', 'X',\n'AMW', 'X\
',\n'AMX', 'X',\n'AMY', 'X',\n'ANA', 'X',\n'ANB', \
'X',\n'ANC', 'X',\n'AND', 'X',\n'ANE', 'X',\n'ANI'\
, 'X',\n'ANL', 'X',\n'ANO', 'X',\n'ANP', 'X',\n'AN\
S', 'X',\n'ANT', 'X',\n'AOE', 'X',\n'AOP', 'X',\n'\
AP1', 'X',\n'AP2', 'X',\n'AP3', 'X',\n'AP4', 'X',\\
n'AP5', 'X',\n'AP6', 'X',\n'APA', 'X',\n'APB', 'X'\
,\n'APC', 'X',\n'APE', 'F',\n'APF', 'X',\n'APG', '\
X',\n'APH', 'A',\n'API', 'X',\n'APL', 'X',\n'APM',\
 'X',\n'APN', 'G',\n'APP', 'X',\n'APQ', 'X',\n'APR\
', 'X',\n'APS', 'X',\n'APT', 'X',\n'APU', 'X',\n'A\
PX', 'X',\n'APY', 'X',\n'APZ', 'X',\n'AQS', 'X',\n\
'AR1', 'X',\n'AR2', 'X',\n'ARA', 'X',\n'ARB', 'X',\
\n'ARC', 'X',\n'ARD', 'X',\n'ARG', 'R',\n'ARH', 'X\
',\n'ARI', 'X',\n'ARM', 'R',\n'ARN', 'X',\n'ARO', \
'R',\n'ARP', 'X',\n'ARQ', 'X',\n'ARS', 'X',\n'AS1'\
, 'R',\n'AS2', 'X',\n'ASA', 'D',\n'ASB', 'D',\n'AS\
C', 'X',\n'ASD', 'X',\n'ASE', 'X',\n'ASF', 'X',\n'\
ASI', 'X',\n'ASK', 'D',\n'ASL', 'X',\n'ASM', 'N',\\
n'ASO', 'X',\n'ASP', 'D',\n'ASQ', 'X',\n'ASU', 'X'\
,\n'ATA', 'X',\n'ATC', 'X',\n'ATD', 'X',\n'ATF', '\
X',\n'ATG', 'X',\n'ATH', 'X',\n'ATM', 'X',\n'ATO',\
 'X',\n'ATP', 'X',\n'ATQ', 'X',\n'ATR', 'X',\n'ATT\
', 'X',\n'ATY', 'X',\n'ATZ', 'X',\n'AUC', 'X',\n'A\
UR', 'X',\n'AVG', 'X',\n'AXP', 'X',\n'AYA', 'A',\n\
'AZ2', 'X',\n'AZA', 'X',\n'AZC', 'X',\n'AZD', 'X',\
\n'AZE', 'X',\n'AZI', 'X',\n'AZL', 'X',\n'AZM', 'X\
',\n'AZR', 'X',\n'AZT', 'X',\n'B12', 'X',\n'B1F', \
'F',\n'B2A', 'A',\n'B2F', 'F',\n'B2I', 'I',\n'B2V'\
, 'V',\n'B3I', 'X',\n'B3P', 'X',\n'B7G', 'X',\n'B9\
6', 'X',\n'B9A', 'X',\n'BA1', 'X',\n'BAA', 'X',\n'\
BAB', 'X',\n'BAC', 'X',\n'BAF', 'X',\n'BAH', 'X',\\
n'BAI', 'X',\n'BAK', 'X',\n'BAL', 'A',\n'BAM', 'X'\
,\n'BAO', 'X',\n'BAP', 'X',\n'BAR', 'X',\n'BAS', '\
X',\n'BAT', 'F',\n'BAY', 'X',\n'BAZ', 'X',\n'BB1',\
 'X',\n'BB2', 'X',\n'BBA', 'X',\n'BBH', 'X',\n'BBS\
', 'X',\n'BBT', 'X',\n'BBZ', 'X',\n'BCA', 'X',\n'B\
CB', 'X',\n'BCC', 'X',\n'BCD', 'X',\n'BCL', 'X',\n\
'BCN', 'X',\n'BCR', 'X',\n'BCS', 'C',\n'BCT', 'X',\
\n'BCY', 'X',\n'BCZ', 'X',\n'BDA', 'X',\n'BDG', 'X\
',\n'BDK', 'X',\n'BDM', 'X',\n'BDN', 'X',\n'BDS', \
'X',\n'BE1', 'X',\n'BE2', 'X',\n'BEA', 'X',\n'BEF'\
, 'X',\n'BEN', 'X',\n'BEO', 'X',\n'BEP', 'X',\n'BE\
R', 'X',\n'BES', 'X',\n'BET', 'X',\n'BEZ', 'X',\n'\
BF2', 'X',\n'BFA', 'X',\n'BFD', 'X',\n'BFP', 'X',\\
n'BFS', 'X',\n'BFU', 'X',\n'BG6', 'X',\n'BGF', 'X'\
,\n'BGG', 'X',\n'BGL', 'X',\n'BGN', 'X',\n'BGP', '\
X',\n'BGX', 'X',\n'BH4', 'X',\n'BHA', 'X',\n'BHC',\
 'X',\n'BHD', 'D',\n'BHO', 'X',\n'BHS', 'X',\n'BIC\
', 'X',\n'BIN', 'X',\n'BIO', 'X',\n'BIP', 'X',\n'B\
IS', 'X',\n'BIZ', 'X',\n'BJH', 'X',\n'BJI', 'X',\n\
'BJP', 'X',\n'BLA', 'X',\n'BLB', 'X',\n'BLE', 'L',\
\n'BLG', 'P',\n'BLI', 'X',\n'BLM', 'X',\n'BLV', 'X\
',\n'BLY', 'K',\n'BM1', 'X',\n'BM2', 'X',\n'BM5', \
'X',\n'BM9', 'X',\n'BMA', 'X',\n'BMD', 'X',\n'BME'\
, 'X',\n'BMP', 'X',\n'BMQ', 'X',\n'BMS', 'X',\n'BM\
T', 'T',\n'BMU', 'X',\n'BMY', 'X',\n'BMZ', 'X',\n'\
BNA', 'X',\n'BNG', 'X',\n'BNI', 'X',\n'BNN', 'F',\\
n'BNO', 'L',\n'BNS', 'X',\n'BNZ', 'X',\n'BO3', 'X'\
,\n'BO4', 'X',\n'BOC', 'X',\n'BOG', 'X',\n'BOM', '\
X',\n'BOT', 'X',\n'BOX', 'X',\n'BOZ', 'X',\n'BPA',\
 'X',\n'BPB', 'X',\n'BPD', 'X',\n'BPG', 'X',\n'BPH\
', 'X',\n'BPI', 'X',\n'BPJ', 'X',\n'BPM', 'X',\n'B\
PN', 'X',\n'BPO', 'X',\n'BPP', 'X',\n'BPT', 'X',\n\
'BPY', 'X',\n'BRB', 'X',\n'BRC', 'X',\n'BRE', 'X',\
\n'BRI', 'X',\n'BRL', 'X',\n'BRM', 'X',\n'BRN', 'X\
',\n'BRO', 'X',\n'BRS', 'X',\n'BRU', 'X',\n'BRZ', \
'X',\n'BSB', 'X',\n'BSI', 'X',\n'BSP', 'X',\n'BT1'\
, 'X',\n'BT2', 'X',\n'BT3', 'X',\n'BTA', 'L',\n'BT\
B', 'X',\n'BTC', 'C',\n'BTD', 'X',\n'BTN', 'X',\n'\
BTP', 'X',\n'BTR', 'W',\n'BU1', 'X',\n'BUA', 'X',\\
n'BUB', 'X',\n'BUC', 'X',\n'BUG', 'X',\n'BUL', 'X'\
,\n'BUM', 'X',\n'BUQ', 'X',\n'BUT', 'X',\n'BVD', '\
X',\n'BX3', 'X',\n'BYS', 'X',\n'BZ1', 'X',\n'BZA',\
 'X',\n'BZB', 'X',\n'BZC', 'X',\n'BZD', 'X',\n'BZF\
', 'X',\n'BZI', 'X',\n'BZM', 'X',\n'BZO', 'X',\n'B\
ZP', 'X',\n'BZQ', 'X',\n'BZS', 'X',\n'BZT', 'X',\n\
'C02', 'X',\n'C11', 'X',\n'C1O', 'X',\n'C20', 'X',\
\n'C24', 'X',\n'C2F', 'X',\n'C2O', 'X',\n'C2P', 'X\
',\n'C3M', 'X',\n'C3P', 'X',\n'C3X', 'X',\n'C48', \
'X',\n'C4M', 'X',\n'C4X', 'X',\n'C5C', 'X',\n'C5M'\
, 'X',\n'C5P', 'X',\n'C5X', 'X',\n'C60', 'X',\n'C6\
C', 'X',\n'C6M', 'X',\n'C78', 'X',\n'C8E', 'X',\n'\
CA3', 'X',\n'CA5', 'X',\n'CAA', 'X',\n'CAB', 'X',\\
n'CAC', 'X',\n'CAD', 'X',\n'CAF', 'C',\n'CAG', 'X'\
,\n'CAH', 'X',\n'CAL', 'X',\n'CAM', 'X',\n'CAN', '\
X',\n'CAO', 'X',\n'CAP', 'X',\n'CAQ', 'X',\n'CAR',\
 'X',\n'CAS', 'C',\n'CAT', 'X',\n'CAV', 'X',\n'CAY\
', 'C',\n'CAZ', 'X',\n'CB3', 'X',\n'CB4', 'X',\n'C\
BA', 'X',\n'CBD', 'X',\n'CBG', 'X',\n'CBI', 'X',\n\
'CBL', 'X',\n'CBM', 'X',\n'CBN', 'X',\n'CBO', 'X',\
\n'CBP', 'X',\n'CBS', 'X',\n'CBX', 'X',\n'CBZ', 'X\
',\n'CC0', 'X',\n'CC1', 'X',\n'CCC', 'X',\n'CCH', \
'X',\n'CCI', 'X',\n'CCM', 'X',\n'CCN', 'X',\n'CCO'\
, 'X',\n'CCP', 'X',\n'CCR', 'X',\n'CCS', 'C',\n'CC\
V', 'X',\n'CCY', 'X',\n'CD1', 'X',\n'CDC', 'X',\n'\
CDE', 'X',\n'CDF', 'X',\n'CDI', 'X',\n'CDL', 'X',\\
n'CDM', 'X',\n'CDP', 'X',\n'CDR', 'X',\n'CDU', 'X'\
,\n'CE1', 'X',\n'CEA', 'C',\n'CEB', 'X',\n'CEC', '\
X',\n'CED', 'X',\n'CEF', 'X',\n'CEH', 'X',\n'CEM',\
 'X',\n'CEO', 'X',\n'CEP', 'X',\n'CEQ', 'X',\n'CER\
', 'X',\n'CES', 'G',\n'CET', 'X',\n'CFC', 'X',\n'C\
FF', 'X',\n'CFM', 'X',\n'CFO', 'X',\n'CFP', 'X',\n\
'CFS', 'X',\n'CFX', 'X',\n'CGN', 'X',\n'CGP', 'X',\
\n'CGS', 'X',\n'CGU', 'E',\n'CH2', 'X',\n'CH3', 'X\
',\n'CHA', 'X',\n'CHB', 'X',\n'CHD', 'X',\n'CHF', \
'X',\n'CHG', 'G',\n'CHI', 'X',\n'CHN', 'X',\n'CHO'\
, 'X',\n'CHP', 'G',\n'CHR', 'X',\n'CHS', 'F',\n'CH\
T', 'X',\n'CHX', 'X',\n'CIC', 'X',\n'CIN', 'X',\n'\
CIP', 'X',\n'CIR', 'X',\n'CIT', 'X',\n'CIU', 'X',\\
n'CKI', 'X',\n'CL1', 'X',\n'CL2', 'X',\n'CLA', 'X'\
,\n'CLB', 'A',\n'CLC', 'S',\n'CLD', 'A',\n'CLE', '\
L',\n'CLF', 'X',\n'CLK', 'S',\n'CLL', 'X',\n'CLM',\
 'X',\n'CLN', 'X',\n'CLO', 'X',\n'CLP', 'X',\n'CLQ\
', 'X',\n'CLR', 'X',\n'CLS', 'X',\n'CLT', 'X',\n'C\
LX', 'X',\n'CLY', 'X',\n'CMA', 'R',\n'CMC', 'X',\n\
'CMD', 'X',\n'CME', 'C',\n'CMG', 'X',\n'CMK', 'X',\
\n'CMN', 'X',\n'CMO', 'X',\n'CMP', 'X',\n'CMR', 'X\
',\n'CMS', 'X',\n'CMT', 'C',\n'CMX', 'X',\n'CNA', \
'X',\n'CNC', 'X',\n'CND', 'X',\n'CNH', 'X',\n'CNM'\
, 'X',\n'CNN', 'X',\n'CNO', 'X',\n'CNP', 'X',\n'CO\
2', 'X',\n'CO3', 'X',\n'CO5', 'X',\n'CO8', 'X',\n'\
COA', 'X',\n'COB', 'X',\n'COC', 'X',\n'COD', 'X',\\
n'COE', 'X',\n'COF', 'X',\n'COH', 'X',\n'COI', 'X'\
,\n'COJ', 'X',\n'COL', 'X',\n'COM', 'X',\n'CON', '\
X',\n'COP', 'X',\n'COR', 'X',\n'COS', 'X',\n'COT',\
 'X',\n'COY', 'X',\n'CP1', 'G',\n'CP2', 'X',\n'CP4\
', 'X',\n'CPA', 'X',\n'CPB', 'X',\n'CPC', 'X',\n'C\
PD', 'X',\n'CPG', 'X',\n'CPH', 'X',\n'CPI', 'X',\n\
'CPM', 'X',\n'CPN', 'G',\n'CPO', 'X',\n'CPP', 'X',\
\n'CPQ', 'X',\n'CPR', 'X',\n'CPS', 'X',\n'CPT', 'X\
',\n'CPU', 'X',\n'CPV', 'X',\n'CPY', 'X',\n'CR1', \
'X',\n'CR6', 'X',\n'CRA', 'X',\n'CRB', 'X',\n'CRC'\
, 'X',\n'CRG', 'X',\n'CRH', 'X',\n'CRO', 'T',\n'CR\
P', 'X',\n'CRQ', 'X',\n'CRS', 'X',\n'CRT', 'X',\n'\
CRY', 'X',\n'CSA', 'C',\n'CSB', 'X',\n'CSD', 'C',\\
n'CSE', 'C',\n'CSH', 'X',\n'CSI', 'X',\n'CSN', 'X'\
,\n'CSO', 'C',\n'CSP', 'C',\n'CSR', 'C',\n'CSS', '\
C',\n'CST', 'X',\n'CSW', 'C',\n'CSX', 'C',\n'CSY',\
 'X',\n'CSZ', 'C',\n'CT3', 'X',\n'CTA', 'X',\n'CTB\
', 'X',\n'CTC', 'X',\n'CTD', 'X',\n'CTH', 'T',\n'C\
TO', 'X',\n'CTP', 'X',\n'CTR', 'X',\n'CTS', 'X',\n\
'CTT', 'X',\n'CTY', 'X',\n'CTZ', 'X',\n'CU1', 'X',\
\n'CUA', 'X',\n'CUC', 'X',\n'CUL', 'X',\n'CUO', 'X\
',\n'CUZ', 'X',\n'CVI', 'X',\n'CXF', 'X',\n'CXL', \
'X',\n'CXM', 'M',\n'CXN', 'X',\n'CXP', 'X',\n'CXS'\
, 'X',\n'CY1', 'C',\n'CY3', 'X',\n'CYB', 'X',\n'CY\
C', 'X',\n'CYF', 'C',\n'CYG', 'C',\n'CYH', 'X',\n'\
CYL', 'X',\n'CYM', 'C',\n'CYN', 'X',\n'CYO', 'X',\\
n'CYP', 'X',\n'CYQ', 'C',\n'CYS', 'C',\n'CYU', 'X'\
,\n'CYY', 'X',\n'CYZ', 'X',\n'CZH', 'X',\n'CZZ', '\
C',\n'D12', 'X',\n'D13', 'X',\n'D16', 'X',\n'D18',\
 'X',\n'D19', 'X',\n'D1P', 'X',\n'D24', 'X',\n'D34\
', 'X',\n'D35', 'X',\n'D4D', 'X',\n'D4T', 'X',\n'D\
6G', 'X',\n'DA2', 'R',\n'DA3', 'X',\n'DA6', 'X',\n\
'DA7', 'X',\n'DAA', 'X',\n'DAB', 'X',\n'DAC', 'X',\
\n'DAD', 'X',\n'DAE', 'X',\n'DAF', 'X',\n'DAG', 'X\
',\n'DAH', 'A',\n'DAJ', 'X',\n'DAK', 'X',\n'DAL', \
'A',\n'DAM', 'A',\n'DAN', 'X',\n'DAO', 'X',\n'DAP'\
, 'X',\n'DAQ', 'X',\n'DAR', 'R',\n'DAS', 'D',\n'DA\
T', 'X',\n'DAU', 'X',\n'DAV', 'X',\n'DBA', 'X',\n'\
DBD', 'X',\n'DBF', 'X',\n'DBG', 'X',\n'DBI', 'X',\\
n'DBV', 'X',\n'DBY', 'Y',\n'DCA', 'X',\n'DCB', 'X'\
,\n'DCE', 'X',\n'DCF', 'X',\n'DCG', 'X',\n'DCH', '\
X',\n'DCI', 'I',\n'DCL', 'X',\n'DCM', 'X',\n'DCP',\
 'X',\n'DCS', 'X',\n'DCT', 'X',\n'DCY', 'C',\n'DCZ\
', 'X',\n'DDA', 'X',\n'DDB', 'X',\n'DDC', 'X',\n'D\
DF', 'X',\n'DDG', 'X',\n'DDH', 'X',\n'DDL', 'X',\n\
'DDM', 'X',\n'DDO', 'L',\n'DDP', 'X',\n'DDQ', 'X',\
\n'DDT', 'Y',\n'DDU', 'X',\n'DEA', 'X',\n'DEB', 'X\
',\n'DEC', 'X',\n'DEF', 'X',\n'DEL', 'X',\n'DEM', \
'X',\n'DEN', 'X',\n'DEP', 'X',\n'DEQ', 'X',\n'DES'\
, 'X',\n'DET', 'X',\n'DFC', 'X',\n'DFG', 'X',\n'DF\
I', 'X',\n'DFL', 'X',\n'DFO', 'X',\n'DFP', 'X',\n'\
DFR', 'X',\n'DFT', 'X',\n'DFV', 'X',\n'DFX', 'X',\\
n'DG2', 'X',\n'DG3', 'X',\n'DG6', 'X',\n'DGA', 'X'\
,\n'DGD', 'X',\n'DGG', 'X',\n'DGL', 'E',\n'DGN', '\
Q',\n'DGP', 'X',\n'DGT', 'X',\n'DGX', 'X',\n'DH2',\
 'X',\n'DHA', 'A',\n'DHB', 'X',\n'DHC', 'X',\n'DHD\
', 'X',\n'DHE', 'X',\n'DHF', 'X',\n'DHG', 'X',\n'D\
HI', 'H',\n'DHL', 'X',\n'DHM', 'X',\n'DHN', 'V',\n\
'DHP', 'X',\n'DHQ', 'X',\n'DHR', 'X',\n'DHS', 'X',\
\n'DHT', 'X',\n'DHU', 'X',\n'DHY', 'X',\n'DHZ', 'X\
',\n'DI2', 'X',\n'DI3', 'G',\n'DI4', 'X',\n'DI5', \
'X',\n'DIA', 'X',\n'DIC', 'X',\n'DIF', 'X',\n'DIG'\
, 'X',\n'DII', 'X',\n'DIL', 'I',\n'DIM', 'X',\n'DI\
O', 'X',\n'DIP', 'X',\n'DIQ', 'X',\n'DIS', 'X',\n'\
DIT', 'X',\n'DIV', 'V',\n'DIX', 'X',\n'DIY', 'X',\\
n'DKA', 'X',\n'DLA', 'X',\n'DLE', 'L',\n'DLF', 'X'\
,\n'DLS', 'K',\n'DLY', 'K',\n'DM1', 'X',\n'DM2', '\
X',\n'DM3', 'X',\n'DM4', 'X',\n'DM5', 'X',\n'DM6',\
 'X',\n'DM7', 'X',\n'DM8', 'X',\n'DM9', 'X',\n'DMA\
', 'X',\n'DMB', 'X',\n'DMC', 'X',\n'DMD', 'X',\n'D\
ME', 'X',\n'DMF', 'X',\n'DMG', 'G',\n'DMH', 'N',\n\
'DMI', 'X',\n'DMJ', 'X',\n'DML', 'X',\n'DMM', 'X',\
\n'DMN', 'X',\n'DMO', 'X',\n'DMP', 'X',\n'DMQ', 'X\
',\n'DMR', 'X',\n'DMS', 'X',\n'DMT', 'X',\n'DMV', \
'X',\n'DMY', 'X',\n'DNC', 'X',\n'DND', 'X',\n'DNH'\
, 'X',\n'DNJ', 'X',\n'DNN', 'X',\n'DNP', 'X',\n'DN\
Q', 'X',\n'DNR', 'X',\n'DO2', 'X',\n'DO3', 'X',\n'\
DOA', 'X',\n'DOB', 'X',\n'DOC', 'X',\n'DOH', 'D',\\
n'DOM', 'X',\n'DOS', 'X',\n'DOX', 'X',\n'DP5', 'X'\
,\n'DP7', 'X',\n'DPA', 'X',\n'DPC', 'X',\n'DPD', '\
X',\n'DPE', 'X',\n'DPG', 'X',\n'DPH', 'F',\n'DPM',\
 'X',\n'DPN', 'F',\n'DPO', 'X',\n'DPP', 'X',\n'DPR\
', 'P',\n'DPS', 'X',\n'DPT', 'X',\n'DPX', 'X',\n'D\
PY', 'X',\n'DPZ', 'X',\n'DQH', 'X',\n'DQN', 'X',\n\
'DR1', 'X',\n'DRB', 'X',\n'DRC', 'X',\n'DRI', 'X',\
\n'DRP', 'X',\n'DRT', 'X',\n'DRU', 'X',\n'DSA', 'X\
',\n'DSB', 'X',\n'DSC', 'X',\n'DSD', 'X',\n'DSE', \
'S',\n'DSI', 'X',\n'DSN', 'S',\n'DSP', 'D',\n'DSR'\
, 'X',\n'DSS', 'X',\n'DSX', 'X',\n'DSY', 'X',\n'DT\
B', 'X',\n'DTD', 'X',\n'DTH', 'T',\n'DTN', 'X',\n'\
DTO', 'X',\n'DTP', 'X',\n'DTQ', 'X',\n'DTR', 'W',\\
n'DTT', 'X',\n'DTY', 'Y',\n'DUD', 'X',\n'DUO', 'X'\
,\n'DUR', 'X',\n'DUT', 'X',\n'DVA', 'V',\n'DVR', '\
X',\n'DX9', 'X',\n'DXA', 'X',\n'DXB', 'X',\n'DXC',\
 'X',\n'DXG', 'X',\n'DXX', 'X',\n'DZF', 'X',\n'E09\
', 'X',\n'E20', 'X',\n'E2P', 'X',\n'E3G', 'X',\n'E\
4N', 'X',\n'E4P', 'X',\n'E64', 'X',\n'E6C', 'X',\n\
'E96', 'X',\n'E97', 'X',\n'EA2', 'X',\n'EAA', 'X',\
\n'EAP', 'X',\n'EBP', 'X',\n'EBW', 'X',\n'ECO', 'X\
',\n'EDA', 'X',\n'EDC', 'X',\n'EDE', 'X',\n'EDO', \
'X',\n'EDR', 'X',\n'EEB', 'X',\n'EEE', 'X',\n'EFC'\
, 'X',\n'EFZ', 'X',\n'EG1', 'X',\n'EG2', 'X',\n'EG\
3', 'X',\n'EGC', 'X',\n'EGL', 'X',\n'EHP', 'A',\n'\
EIC', 'X',\n'EJT', 'X',\n'ELA', 'X',\n'EMB', 'X',\\
n'EMC', 'X',\n'EMD', 'X',\n'EMM', 'X',\n'EMO', 'X'\
,\n'EMP', 'X',\n'EMR', 'X',\n'ENA', 'X',\n'ENC', '\
X',\n'ENH', 'X',\n'ENO', 'X',\n'ENP', 'X',\n'EOA',\
 'X',\n'EOH', 'X',\n'EOT', 'X',\n'EOX', 'X',\n'EPA\
', 'X',\n'EPE', 'X',\n'EPH', 'X',\n'EPI', 'X',\n'E\
PN', 'X',\n'EPO', 'X',\n'EPT', 'X',\n'EPU', 'X',\n\
'EPX', 'X',\n'EPY', 'X',\n'EQI', 'X',\n'EQP', 'X',\
\n'EQU', 'X',\n'ERG', 'X',\n'ERI', 'X',\n'ERY', 'X\
',\n'ESC', 'X',\n'ESD', 'X',\n'ESI', 'X',\n'ESO', \
'X',\n'ESP', 'X',\n'EST', 'X',\n'ESX', 'X',\n'ETA'\
, 'X',\n'ETC', 'X',\n'ETD', 'X',\n'ETF', 'X',\n'ET\
H', 'X',\n'ETI', 'X',\n'ETN', 'X',\n'ETO', 'X',\n'\
ETP', 'X',\n'ETR', 'X',\n'ETS', 'X',\n'ETY', 'X',\\
n'EU3', 'X',\n'EUG', 'X',\n'EYS', 'C',\n'F09', 'X'\
,\n'F2B', 'X',\n'F3S', 'X',\n'F42', 'X',\n'F43', '\
X',\n'F4S', 'X',\n'F6B', 'X',\n'F6P', 'X',\n'F89',\
 'X',\n'FA1', 'X',\n'FA5', 'F',\n'FAA', 'X',\n'FAB\
', 'X',\n'FAC', 'X',\n'FAD', 'X',\n'FAF', 'X',\n'F\
AG', 'X',\n'FAM', 'X',\n'FAR', 'X',\n'FAS', 'X',\n\
'FAT', 'X',\n'FBA', 'X',\n'FBE', 'X',\n'FBI', 'X',\
\n'FBP', 'X',\n'FBQ', 'X',\n'FBS', 'X',\n'FBT', 'X\
',\n'FBU', 'X',\n'FCA', 'X',\n'FCB', 'X',\n'FCI', \
'X',\n'FCN', 'X',\n'FCO', 'X',\n'FCR', 'X',\n'FCT'\
, 'X',\n'FCX', 'X',\n'FCY', 'C',\n'FD1', 'F',\n'FD\
2', 'F',\n'FD3', 'F',\n'FD4', 'F',\n'FDA', 'X',\n'\
FDC', 'X',\n'FDI', 'X',\n'FDP', 'X',\n'FDS', 'X',\\
n'FE2', 'X',\n'FEA', 'X',\n'FEL', 'X',\n'FEM', 'X'\
,\n'FEN', 'X',\n'FEO', 'X',\n'FEP', 'X',\n'FER', '\
X',\n'FES', 'X',\n'FFB', 'X',\n'FFC', 'X',\n'FFF',\
 'X',\n'FFO', 'X',\n'FGL', 'G',\n'FHB', 'X',\n'FHC\
', 'X',\n'FHP', 'X',\n'FHU', 'X',\n'FID', 'X',\n'F\
II', 'X',\n'FIP', 'X',\n'FK5', 'X',\n'FKA', 'X',\n\
'FKI', 'X',\n'FKP', 'X',\n'FL2', 'X',\n'FL9', 'X',\
\n'FLA', 'A',\n'FLC', 'X',\n'FLD', 'X',\n'FLE', 'L\
',\n'FLF', 'X',\n'FLO', 'X',\n'FLP', 'X',\n'FLT', \
'Y',\n'FLU', 'X',\n'FLX', 'X',\n'FM1', 'X',\n'FM2'\
, 'X',\n'FMA', 'X',\n'FMB', 'X',\n'FMC', 'X',\n'FM\
E', 'M',\n'FMN', 'X',\n'FMP', 'X',\n'FMR', 'X',\n'\
FMS', 'X',\n'FMT', 'X',\n'FNE', 'X',\n'FNP', 'X',\\
n'FNS', 'X',\n'FOC', 'X',\n'FOE', 'X',\n'FOG', 'F'\
,\n'FOH', 'X',\n'FOK', 'X',\n'FOL', 'X',\n'FON', '\
X',\n'FOP', 'X',\n'FOR', 'X',\n'FOS', 'X',\n'FPA',\
 'X',\n'FPC', 'X',\n'FPI', 'X',\n'FPO', 'X',\n'FPP\
', 'X',\n'FPT', 'X',\n'FQP', 'X',\n'FRA', 'X',\n'F\
RD', 'F',\n'FRU', 'X',\n'FS3', 'X',\n'FS4', 'X',\n\
'FSB', 'X',\n'FSO', 'X',\n'FSX', 'X',\n'FTC', 'X',\
\n'FTP', 'X',\n'FTR', 'W',\n'FTT', 'X',\n'FTY', 'Y\
',\n'FUA', 'X',\n'FUC', 'X',\n'FUM', 'X',\n'FUP', \
'X',\n'FVF', 'X',\n'FXP', 'X',\n'FXV', 'X',\n'FYA'\
, 'F',\n'G16', 'X',\n'G1P', 'X',\n'G20', 'X',\n'G2\
1', 'X',\n'G23', 'X',\n'G26', 'X',\n'G28', 'X',\n'\
G2F', 'X',\n'G37', 'X',\n'G39', 'X',\n'G3H', 'X',\\
n'G3P', 'X',\n'G4D', 'X',\n'G6D', 'X',\n'G6P', 'X'\
,\n'G6Q', 'X',\n'G7M', 'X',\n'GA2', 'X',\n'GAA', '\
X',\n'GAB', 'X',\n'GAC', 'X',\n'GAI', 'X',\n'GAL',\
 'X',\n'GAM', 'X',\n'GAN', 'X',\n'GAO', 'X',\n'GAP\
', 'X',\n'GAR', 'G',\n'GAS', 'X',\n'GAT', 'X',\n'G\
BC', 'X',\n'GBI', 'X',\n'GBP', 'X',\n'GBS', 'X',\n\
'GBX', 'X',\n'GC4', 'X',\n'GCA', 'X',\n'GCD', 'X',\
\n'GCG', 'G',\n'GCH', 'G',\n'GCK', 'X',\n'GCL', 'X\
',\n'GCM', 'X',\n'GCN', 'X',\n'GCO', 'X',\n'GCP', \
'X',\n'GCR', 'X',\n'GCS', 'X',\n'GCU', 'X',\n'GD3'\
, 'X',\n'GDB', 'X',\n'GDM', 'X',\n'GDN', 'X',\n'GD\
P', 'X',\n'GDS', 'X',\n'GDU', 'X',\n'GE1', 'X',\n'\
GE2', 'X',\n'GE3', 'X',\n'GEA', 'X',\n'GEL', 'X',\\
n'GEM', 'X',\n'GEN', 'X',\n'GEP', 'X',\n'GER', 'X'\
,\n'GFP', 'X',\n'GGB', 'X',\n'GGL', 'E',\n'GGP', '\
X',\n'GHP', 'G',\n'GIP', 'X',\n'GIS', 'X',\n'GKR',\
 'X',\n'GL2', 'X',\n'GL3', 'G',\n'GL4', 'X',\n'GL5\
', 'X',\n'GL7', 'X',\n'GL9', 'X',\n'GLA', 'X',\n'G\
LB', 'X',\n'GLC', 'X',\n'GLD', 'X',\n'GLE', 'X',\n\
'GLF', 'X',\n'GLG', 'X',\n'GLH', 'Q',\n'GLI', 'X',\
\n'GLL', 'X',\n'GLM', 'G',\n'GLN', 'Q',\n'GLO', 'X\
',\n'GLP', 'X',\n'GLR', 'X',\n'GLS', 'X',\n'GLT', \
'X',\n'GLU', 'E',\n'GLV', 'X',\n'GLW', 'X',\n'GLY'\
, 'G',\n'GLZ', 'X',\n'GM1', 'X',\n'GMA', 'X',\n'GM\
C', 'X',\n'GMH', 'X',\n'GMP', 'X',\n'GMY', 'X',\n'\
GN7', 'X',\n'GNA', 'X',\n'GNB', 'X',\n'GNH', 'X',\\
n'GNP', 'X',\n'GNT', 'X',\n'GOA', 'X',\n'GOL', 'X'\
,\n'GOX', 'X',\n'GP1', 'X',\n'GP3', 'X',\n'GP4', '\
X',\n'GP6', 'X',\n'GP8', 'X',\n'GPB', 'E',\n'GPC',\
 'X',\n'GPE', 'X',\n'GPG', 'X',\n'GPI', 'X',\n'GPJ\
', 'X',\n'GPL', 'K',\n'GPM', 'X',\n'GPN', 'G',\n'G\
PP', 'X',\n'GPR', 'X',\n'GPS', 'X',\n'GPX', 'X',\n\
'GR1', 'X',\n'GR3', 'X',\n'GR4', 'X',\n'GSA', 'X',\
\n'GSB', 'X',\n'GSC', 'G',\n'GSE', 'S',\n'GSH', 'X\
',\n'GSP', 'X',\n'GSR', 'X',\n'GSS', 'X',\n'GT9', \
'C',\n'GTA', 'X',\n'GTB', 'X',\n'GTD', 'X',\n'GTE'\
, 'X',\n'GTH', 'T',\n'GTN', 'X',\n'GTO', 'X',\n'GT\
P', 'X',\n'GTR', 'X',\n'GTS', 'X',\n'GTT', 'X',\n'\
GTX', 'X',\n'GTZ', 'X',\n'GU7', 'X',\n'GUA', 'X',\\
n'GUD', 'X',\n'GUM', 'X',\n'GUN', 'X',\n'GUP', 'X'\
,\n'GUR', 'X',\n'GW3', 'X',\n'GZZ', 'X',\n'H2B', '\
X',\n'H2P', 'H',\n'H2S', 'X',\n'H2U', 'X',\n'H4B',\
 'X',\n'H5M', 'P',\n'H5P', 'X',\n'HAA', 'X',\n'HAB\
', 'X',\n'HAC', 'A',\n'HAD', 'X',\n'HAE', 'X',\n'H\
AG', 'X',\n'HAI', 'X',\n'HAM', 'X',\n'HAP', 'X',\n\
'HAQ', 'X',\n'HAR', 'R',\n'HAS', 'X',\n'HAV', 'V',\
\n'HAX', 'X',\n'HAZ', 'X',\n'HBA', 'X',\n'HBC', 'X\
',\n'HBD', 'X',\n'HBI', 'X',\n'HBO', 'X',\n'HBU', \
'X',\n'HBY', 'X',\n'HC0', 'X',\n'HC1', 'X',\n'HC4'\
, 'X',\n'HCA', 'X',\n'HCC', 'X',\n'HCI', 'X',\n'HC\
S', 'X',\n'HDA', 'X',\n'HDD', 'X',\n'HDF', 'X',\n'\
HDN', 'X',\n'HDS', 'X',\n'HDZ', 'X',\n'HE1', 'X',\\
n'HE6', 'X',\n'HEA', 'X',\n'HEB', 'X',\n'HEC', 'X'\
,\n'HED', 'X',\n'HEE', 'X',\n'HEF', 'X',\n'HEG', '\
X',\n'HEM', 'X',\n'HEN', 'X',\n'HEO', 'X',\n'HEP',\
 'X',\n'HEU', 'X',\n'HEV', 'X',\n'HEX', 'X',\n'HEZ\
', 'X',\n'HF1', 'X',\n'HFA', 'X',\n'HFP', 'X',\n'H\
GA', 'Q',\n'HGB', 'X',\n'HGC', 'X',\n'HGI', 'X',\n\
'HGU', 'X',\n'HHO', 'X',\n'HHP', 'X',\n'HIB', 'X',\
\n'HIC', 'H',\n'HII', 'X',\n'HIN', 'X',\n'HIO', 'X\
',\n'HIP', 'H',\n'HIS', 'H',\n'HLE', 'X',\n'HLT', \
'X',\n'HMA', 'A',\n'HMB', 'X',\n'HMC', 'X',\n'HMD'\
, 'X',\n'HMF', 'A',\n'HMG', 'X',\n'HMH', 'X',\n'HM\
I', 'L',\n'HMM', 'X',\n'HMN', 'X',\n'HMO', 'X',\n'\
HMP', 'X',\n'HMR', 'R',\n'HNI', 'X',\n'HNP', 'X',\\
n'HOA', 'X',\n'HOE', 'X',\n'HOH', 'X',\n'HOM', 'X'\
,\n'HOP', 'X',\n'HOQ', 'X',\n'HP1', 'A',\n'HP2', '\
A',\n'HP3', 'X',\n'HPA', 'X',\n'HPB', 'X',\n'HPC',\
 'X',\n'HPD', 'X',\n'HPE', 'A',\n'HPG', 'X',\n'HPH\
', 'F',\n'HPP', 'X',\n'HPQ', 'F',\n'HPR', 'X',\n'H\
PT', 'X',\n'HPY', 'X',\n'HQO', 'X',\n'HQQ', 'X',\n\
'HQU', 'X',\n'HRG', 'R',\n'HRI', 'X',\n'HSA', 'X',\
\n'HSE', 'S',\n'HSF', 'X',\n'HSM', 'X',\n'HSO', 'H\
',\n'HSP', 'X',\n'HT1', 'X',\n'HT2', 'X',\n'HTA', \
'X',\n'HTL', 'X',\n'HTO', 'X',\n'HTP', 'X',\n'HTR'\
, 'W',\n'HUP', 'X',\n'HUX', 'X',\n'HV5', 'A',\n'HV\
7', 'X',\n'HV8', 'X',\n'HXA', 'X',\n'HXC', 'X',\n'\
HXP', 'X',\n'HY1', 'X',\n'HYA', 'X',\n'HYB', 'X',\\
n'HYD', 'X',\n'HYG', 'X',\n'HYP', 'P',\n'I06', 'X'\
,\n'I10', 'X',\n'I11', 'X',\n'I17', 'X',\n'I2P', '\
X',\n'I3N', 'X',\n'I3P', 'X',\n'I40', 'X',\n'I48',\
 'X',\n'I4B', 'X',\n'I52', 'X',\n'I5P', 'X',\n'I84\
', 'G',\n'IAG', 'G',\n'IAS', 'X',\n'IB2', 'X',\n'I\
BB', 'X',\n'IBP', 'X',\n'IBR', 'X',\n'IBS', 'X',\n\
'IBZ', 'X',\n'IC1', 'X',\n'ICA', 'X',\n'ICI', 'X',\
\n'ICL', 'X',\n'ICP', 'X',\n'ICT', 'X',\n'ICU', 'X\
',\n'ID2', 'X',\n'IDC', 'X',\n'IDG', 'X',\n'IDH', \
'X',\n'IDM', 'X',\n'IDO', 'X',\n'IDP', 'X',\n'IDR'\
, 'X',\n'IDS', 'X',\n'IDT', 'X',\n'IDU', 'X',\n'IF\
G', 'X',\n'IFP', 'X',\n'IGL', 'X',\n'IGN', 'X',\n'\
IGP', 'X',\n'IGU', 'X',\n'IH1', 'X',\n'IH2', 'X',\\
n'IH3', 'X',\n'IHB', 'X',\n'IHN', 'X',\n'IHP', 'X'\
,\n'IIC', 'X',\n'IIL', 'I',\n'IIP', 'X',\n'IK2', '\
X',\n'IKT', 'X',\n'ILA', 'I',\n'ILE', 'I',\n'ILG',\
 'X',\n'ILO', 'X',\n'ILX', 'I',\n'IM1', 'X',\n'IM2\
', 'X',\n'IMC', 'X',\n'IMD', 'X',\n'IME', 'X',\n'I\
MF', 'X',\n'IMG', 'X',\n'IMH', 'X',\n'IMI', 'X',\n\
'IML', 'I',\n'IMM', 'X',\n'IMN', 'X',\n'IMO', 'X',\
\n'IMP', 'X',\n'IMR', 'X',\n'IMU', 'X',\n'IN0', 'D\
',\n'IN1', 'R',\n'IN2', 'K',\n'IN3', 'L',\n'IN4', \
'X',\n'IN5', 'A',\n'IN6', 'L',\n'IN7', 'X',\n'IN8'\
, 'X',\n'IN9', 'X',\n'INA', 'L',\n'INB', 'X',\n'IN\
C', 'X',\n'IND', 'X',\n'INE', 'X',\n'INF', 'F',\n'\
ING', 'F',\n'INH', 'R',\n'INI', 'X',\n'INJ', 'X',\\
n'INK', 'X',\n'INL', 'X',\n'INM', 'X',\n'INN', 'A'\
,\n'INO', 'X',\n'INP', 'X',\n'INQ', 'X',\n'INR', '\
X',\n'INS', 'X',\n'INT', 'V',\n'INU', 'X',\n'INV',\
 'X',\n'INW', 'X',\n'INX', 'X',\n'INY', 'X',\n'INZ\
', 'X',\n'IOA', 'X',\n'IOB', 'X',\n'IOC', 'X',\n'I\
OD', 'X',\n'IOE', 'X',\n'IOF', 'X',\n'IOH', 'X',\n\
'IOL', 'X',\n'IOP', 'X',\n'IP1', 'X',\n'IP2', 'X',\
\n'IP3', 'X',\n'IP4', 'X',\n'IPA', 'X',\n'IPB', 'X\
',\n'IPD', 'X',\n'IPG', 'G',\n'IPH', 'X',\n'IPL', \
'X',\n'IPM', 'X',\n'IPN', 'X',\n'IPO', 'F',\n'IPP'\
, 'X',\n'IPS', 'X',\n'IPT', 'X',\n'IPU', 'X',\n'IP\
Y', 'A',\n'IQB', 'X',\n'IQP', 'X',\n'IQS', 'X',\n'\
IR3', 'X',\n'IRI', 'X',\n'IRP', 'X',\n'ISA', 'X',\\
n'ISF', 'X',\n'ISO', 'X',\n'ISP', 'X',\n'ISQ', 'X'\
,\n'ISU', 'X',\n'ITM', 'X',\n'ITP', 'X',\n'ITR', '\
W',\n'ITS', 'X',\n'ITU', 'X',\n'IU5', 'X',\n'IUM',\
 'X',\n'IUR', 'X',\n'IVA', 'X',\n'IYG', 'G',\n'IYR\
', 'Y',\n'J77', 'X',\n'J78', 'X',\n'J80', 'X',\n'J\
E2', 'X',\n'JEN', 'X',\n'JST', 'X',\n'K21', 'X',\n\
'KAH', 'X',\n'KAI', 'X',\n'KAM', 'X',\n'KAN', 'X',\
\n'KAP', 'X',\n'KCP', 'X',\n'KCX', 'K',\n'KDO', 'X\
',\n'KEF', 'X',\n'KET', 'X',\n'KGR', 'X',\n'KH1', \
'X',\n'KIF', 'X',\n'KIV', 'V',\n'KNI', 'X',\n'KPH'\
, 'K',\n'KTH', 'X',\n'KTN', 'X',\n'KTP', 'X',\n'KW\
T', 'X',\n'L04', 'X',\n'L1P', 'X',\n'L24', 'E',\n'\
L2P', 'X',\n'L34', 'E',\n'L37', 'E',\n'L3P', 'X',\\
n'L4P', 'X',\n'L75', 'X',\n'LAC', 'X',\n'LAD', 'X'\
,\n'LAK', 'X',\n'LAM', 'X',\n'LAR', 'X',\n'LAT', '\
X',\n'LAX', 'X',\n'LCO', 'X',\n'LCP', 'X',\n'LCS',\
 'X',\n'LDA', 'X',\n'LDO', 'L',\n'LDP', 'X',\n'LEA\
', 'X',\n'LEO', 'X',\n'LEU', 'L',\n'LG2', 'X',\n'L\
G6', 'X',\n'LGC', 'X',\n'LGP', 'X',\n'LHG', 'X',\n\
'LHY', 'F',\n'LI1', 'X',\n'LIG', 'X',\n'LIL', 'X',\
\n'LIM', 'X',\n'LIN', 'X',\n'LIO', 'X',\n'LIP', 'X\
',\n'LLA', 'X',\n'LLP', 'K',\n'LLY', 'K',\n'LMG', \
'X',\n'LML', 'X',\n'LMT', 'X',\n'LMU', 'X',\n'LMZ'\
, 'X',\n'LNK', 'X',\n'LNL', 'X',\n'LNO', 'X',\n'LO\
F', 'X',\n'LOL', 'L',\n'LOM', 'X',\n'LOR', 'X',\n'\
LOS', 'X',\n'LOV', 'L',\n'LOX', 'X',\n'LP1', 'X',\\
n'LP2', 'R',\n'LPA', 'X',\n'LPC', 'X',\n'LPF', 'X'\
,\n'LPL', 'X',\n'LPM', 'X',\n'LPP', 'X',\n'LRB', '\
X',\n'LRU', 'X',\n'LS1', 'X',\n'LS2', 'X',\n'LS3',\
 'X',\n'LS4', 'X',\n'LS5', 'X',\n'LTA', 'X',\n'LTL\
', 'X',\n'LTR', 'W',\n'LUM', 'X',\n'LVS', 'L',\n'L\
XC', 'X',\n'LY2', 'X',\n'LY3', 'X',\n'LYA', 'X',\n\
'LYB', 'X',\n'LYC', 'X',\n'LYD', 'X',\n'LYM', 'K',\
\n'LYN', 'X',\n'LYS', 'K',\n'LYT', 'X',\n'LYW', 'X\
',\n'LYZ', 'K',\n'M1A', 'X',\n'M1G', 'X',\n'M2G', \
'X',\n'M3L', 'K',\n'M6P', 'X',\n'M6T', 'X',\n'M7G'\
, 'X',\n'MA1', 'X',\n'MA2', 'X',\n'MA3', 'X',\n'MA\
4', 'X',\n'MA6', 'X',\n'MAA', 'A',\n'MAB', 'X',\n'\
MAC', 'X',\n'MAE', 'X',\n'MAG', 'X',\n'MAH', 'X',\\
n'MAI', 'R',\n'MAK', 'X',\n'MAL', 'X',\n'MAM', 'X'\
,\n'MAN', 'X',\n'MAO', 'X',\n'MAP', 'X',\n'MAR', '\
X',\n'MAS', 'X',\n'MAT', 'X',\n'MAU', 'X',\n'MAZ',\
 'X',\n'MBA', 'X',\n'MBD', 'X',\n'MBG', 'X',\n'MBH\
', 'X',\n'MBN', 'X',\n'MBO', 'X',\n'MBR', 'X',\n'M\
BS', 'X',\n'MBV', 'X',\n'MBZ', 'X',\n'MCA', 'X',\n\
'MCD', 'X',\n'MCE', 'X',\n'MCG', 'G',\n'MCI', 'X',\
\n'MCN', 'X',\n'MCP', 'X',\n'MCT', 'X',\n'MCY', 'X\
',\n'MD2', 'X',\n'MDA', 'X',\n'MDC', 'X',\n'MDG', \
'X',\n'MDH', 'X',\n'MDL', 'X',\n'MDM', 'X',\n'MDN'\
, 'X',\n'MDP', 'X',\n'ME6', 'X',\n'MEB', 'X',\n'ME\
C', 'X',\n'MEL', 'X',\n'MEN', 'N',\n'MEP', 'X',\n'\
MER', 'X',\n'MES', 'X',\n'MET', 'M',\n'MEV', 'X',\\
n'MF2', 'X',\n'MF3', 'M',\n'MFB', 'X',\n'MFD', 'X'\
,\n'MFU', 'X',\n'MG7', 'X',\n'MGA', 'X',\n'MGB', '\
X',\n'MGD', 'X',\n'MGG', 'R',\n'MGL', 'X',\n'MGN',\
 'Q',\n'MGO', 'X',\n'MGP', 'X',\n'MGR', 'X',\n'MGS\
', 'X',\n'MGT', 'X',\n'MGU', 'X',\n'MGY', 'G',\n'M\
HB', 'X',\n'MHF', 'X',\n'MHL', 'L',\n'MHM', 'X',\n\
'MHO', 'M',\n'MHS', 'H',\n'MHZ', 'X',\n'MIA', 'X',\
\n'MIC', 'X',\n'MID', 'X',\n'MIL', 'X',\n'MIM', 'X\
',\n'MIN', 'G',\n'MIP', 'X',\n'MIS', 'S',\n'MIT', \
'X',\n'MJI', 'X',\n'MK1', 'X',\n'MKC', 'X',\n'MLA'\
, 'X',\n'MLC', 'X',\n'MLE', 'L',\n'MLN', 'X',\n'ML\
T', 'X',\n'MLY', 'K',\n'MLZ', 'K',\n'MM3', 'X',\n'\
MM4', 'X',\n'MMA', 'X',\n'MMC', 'X',\n'MME', 'M',\\
n'MMO', 'R',\n'MMP', 'X',\n'MMQ', 'X',\n'MMT', 'X'\
,\n'MN1', 'X',\n'MN2', 'X',\n'MN3', 'X',\n'MN5', '\
X',\n'MN7', 'X',\n'MN8', 'X',\n'MNA', 'X',\n'MNB',\
 'X',\n'MNC', 'X',\n'MNG', 'X',\n'MNL', 'L',\n'MNO\
', 'X',\n'MNP', 'X',\n'MNQ', 'X',\n'MNS', 'X',\n'M\
NT', 'X',\n'MNV', 'V',\n'MO1', 'X',\n'MO2', 'X',\n\
'MO3', 'X',\n'MO4', 'X',\n'MO5', 'X',\n'MO6', 'X',\
\n'MOA', 'X',\n'MOB', 'X',\n'MOC', 'X',\n'MOE', 'X\
',\n'MOG', 'X',\n'MOH', 'X',\n'MOL', 'X',\n'MOO', \
'X',\n'MOP', 'X',\n'MOR', 'X',\n'MOS', 'X',\n'MOT'\
, 'X',\n'MOX', 'X',\n'MP1', 'X',\n'MP3', 'X',\n'MP\
A', 'X',\n'MPB', 'X',\n'MPC', 'X',\n'MPD', 'X',\n'\
MPG', 'X',\n'MPH', 'M',\n'MPI', 'X',\n'MPJ', 'M',\\
n'MPL', 'X',\n'MPN', 'X',\n'MPO', 'X',\n'MPP', 'X'\
,\n'MPQ', 'G',\n'MPR', 'X',\n'MPS', 'X',\n'MQ0', '\
X',\n'MQ7', 'X',\n'MQ8', 'X',\n'MQ9', 'X',\n'MQI',\
 'X',\n'MR2', 'X',\n'MRC', 'X',\n'MRM', 'X',\n'MRP\
', 'X',\n'MS2', 'X',\n'MSA', 'X',\n'MSB', 'X',\n'M\
SD', 'X',\n'MSE', 'M',\n'MSF', 'X',\n'MSI', 'X',\n\
'MSO', 'M',\n'MSQ', 'X',\n'MST', 'X',\n'MSU', 'X',\
\n'MTA', 'X',\n'MTB', 'X',\n'MTC', 'X',\n'MTD', 'X\
',\n'MTE', 'X',\n'MTF', 'X',\n'MTG', 'X',\n'MTO', \
'X',\n'MTS', 'X',\n'MTT', 'X',\n'MTX', 'X',\n'MTY'\
, 'Y',\n'MUG', 'X',\n'MUP', 'X',\n'MUR', 'X',\n'MV\
A', 'V',\n'MW1', 'X',\n'MW2', 'X',\n'MXA', 'X',\n'\
MXY', 'X',\n'MYA', 'X',\n'MYC', 'X',\n'MYG', 'X',\\
n'MYR', 'X',\n'MYS', 'X',\n'MYT', 'X',\n'MZM', 'X'\
,\n'N1T', 'X',\n'N25', 'X',\n'N2B', 'X',\n'N3T', '\
X',\n'N4B', 'X',\n'NA2', 'X',\n'NA5', 'X',\n'NA6',\
 'X',\n'NAA', 'X',\n'NAB', 'X',\n'NAC', 'X',\n'NAD\
', 'X',\n'NAE', 'X',\n'NAF', 'X',\n'NAG', 'X',\n'N\
AH', 'X',\n'NAI', 'X',\n'NAL', 'A',\n'NAM', 'A',\n\
'NAN', 'X',\n'NAO', 'X',\n'NAP', 'X',\n'NAQ', 'X',\
\n'NAR', 'X',\n'NAS', 'X',\n'NAU', 'X',\n'NAV', 'X\
',\n'NAW', 'X',\n'NAX', 'X',\n'NAY', 'X',\n'NBA', \
'X',\n'NBD', 'X',\n'NBE', 'X',\n'NBG', 'X',\n'NBN'\
, 'X',\n'NBP', 'X',\n'NBS', 'X',\n'NBU', 'X',\n'NC\
A', 'X',\n'NCB', 'A',\n'NCD', 'X',\n'NCH', 'X',\n'\
NCM', 'X',\n'NCN', 'X',\n'NCO', 'X',\n'NCR', 'X',\\
n'NCS', 'X',\n'ND4', 'X',\n'NDA', 'X',\n'NDC', 'X'\
,\n'NDD', 'X',\n'NDO', 'X',\n'NDP', 'X',\n'NDT', '\
X',\n'NEA', 'X',\n'NEB', 'X',\n'NED', 'X',\n'NEM',\
 'H',\n'NEN', 'X',\n'NEO', 'X',\n'NEP', 'H',\n'NEQ\
', 'X',\n'NES', 'X',\n'NET', 'X',\n'NEV', 'X',\n'N\
FA', 'F',\n'NFE', 'X',\n'NFG', 'X',\n'NFP', 'X',\n\
'NFS', 'X',\n'NG6', 'X',\n'NGA', 'X',\n'NGL', 'X',\
\n'NGM', 'X',\n'NGO', 'X',\n'NGP', 'X',\n'NGT', 'X\
',\n'NGU', 'X',\n'NH2', 'X',\n'NH3', 'X',\n'NH4', \
'X',\n'NHD', 'X',\n'NHE', 'X',\n'NHM', 'X',\n'NHP'\
, 'X',\n'NHR', 'X',\n'NHS', 'X',\n'NI1', 'X',\n'NI\
2', 'X',\n'NIC', 'X',\n'NID', 'X',\n'NIK', 'X',\n'\
NIO', 'X',\n'NIP', 'X',\n'NIT', 'X',\n'NIU', 'X',\\
n'NIY', 'Y',\n'NLA', 'X',\n'NLE', 'L',\n'NLG', 'X'\
,\n'NLN', 'L',\n'NLP', 'L',\n'NM1', 'X',\n'NMA', '\
A',\n'NMB', 'X',\n'NMC', 'G',\n'NMD', 'X',\n'NME',\
 'X',\n'NMN', 'X',\n'NMO', 'X',\n'NMQ', 'X',\n'NMX\
', 'X',\n'NMY', 'X',\n'NNH', 'R',\n'NNO', 'X',\n'N\
O2', 'X',\n'NO3', 'X',\n'NOA', 'X',\n'NOD', 'X',\n\
'NOJ', 'X',\n'NON', 'X',\n'NOP', 'X',\n'NOR', 'X',\
\n'NOS', 'X',\n'NOV', 'X',\n'NOX', 'X',\n'NP3', 'X\
',\n'NPA', 'X',\n'NPC', 'X',\n'NPD', 'X',\n'NPE', \
'X',\n'NPF', 'X',\n'NPH', 'C',\n'NPI', 'X',\n'NPL'\
, 'X',\n'NPN', 'X',\n'NPO', 'X',\n'NPP', 'X',\n'NP\
T', 'X',\n'NPY', 'X',\n'NRG', 'R',\n'NRI', 'X',\n'\
NS1', 'X',\n'NS5', 'X',\n'NSP', 'X',\n'NTA', 'X',\\
n'NTB', 'X',\n'NTC', 'X',\n'NTH', 'X',\n'NTM', 'X'\
,\n'NTP', 'X',\n'NTS', 'X',\n'NTU', 'X',\n'NTZ', '\
X',\n'NU1', 'X',\n'NVA', 'V',\n'NVI', 'X',\n'NVP',\
 'X',\n'NW1', 'X',\n'NYP', 'X',\n'O4M', 'X',\n'OAA\
', 'X',\n'OAI', 'X',\n'OAP', 'X',\n'OAR', 'X',\n'O\
AS', 'S',\n'OBA', 'X',\n'OBN', 'X',\n'OC1', 'X',\n\
'OC2', 'X',\n'OC3', 'X',\n'OC4', 'X',\n'OC5', 'X',\
\n'OC6', 'X',\n'OC7', 'X',\n'OCL', 'X',\n'OCM', 'X\
',\n'OCN', 'X',\n'OCO', 'X',\n'OCP', 'X',\n'OCS', \
'C',\n'OCT', 'X',\n'OCV', 'K',\n'OCY', 'C',\n'ODA'\
, 'X',\n'ODS', 'X',\n'OES', 'X',\n'OET', 'X',\n'OF\
1', 'X',\n'OF2', 'X',\n'OF3', 'X',\n'OFL', 'X',\n'\
OFO', 'X',\n'OHE', 'X',\n'OHO', 'X',\n'OHT', 'X',\\
n'OIC', 'X',\n'OIP', 'X',\n'OKA', 'X',\n'OLA', 'X'\
,\n'OLE', 'X',\n'OLI', 'X',\n'OLO', 'X',\n'OMB', '\
X',\n'OMC', 'X',\n'OMD', 'X',\n'OME', 'X',\n'OMG',\
 'X',\n'OMP', 'X',\n'OMT', 'M',\n'OMU', 'X',\n'ONE\
', 'X',\n'ONL', 'L',\n'ONP', 'X',\n'OPA', 'X',\n'O\
PD', 'X',\n'OPE', 'X',\n'OPG', 'X',\n'OPH', 'X',\n\
'OPN', 'X',\n'OPP', 'X',\n'OPR', 'R',\n'ORN', 'X',\
\n'ORO', 'X',\n'ORP', 'X',\n'OSB', 'X',\n'OSS', 'X\
',\n'OTA', 'X',\n'OTB', 'X',\n'OTE', 'X',\n'OTG', \
'X',\n'OUT', 'X',\n'OVA', 'X',\n'OWQ', 'X',\n'OXA'\
, 'X',\n'OXE', 'X',\n'OXI', 'X',\n'OXL', 'X',\n'OX\
M', 'X',\n'OXN', 'X',\n'OXO', 'X',\n'OXP', 'X',\n'\
OXS', 'X',\n'OXY', 'X',\n'P11', 'A',\n'P24', 'X',\\
n'P28', 'X',\n'P2P', 'X',\n'P2U', 'X',\n'P3M', 'X'\
,\n'P4C', 'X',\n'P4P', 'X',\n'P5P', 'X',\n'P6G', '\
X',\n'PA1', 'X',\n'PA2', 'X',\n'PA3', 'X',\n'PA4',\
 'X',\n'PA5', 'X',\n'PAA', 'X',\n'PAB', 'X',\n'PAC\
', 'X',\n'PAD', 'X',\n'PAE', 'X',\n'PAG', 'X',\n'P\
AH', 'X',\n'PAI', 'X',\n'PAL', 'D',\n'PAM', 'X',\n\
'PAN', 'X',\n'PAO', 'X',\n'PAP', 'A',\n'PAQ', 'F',\
\n'PAR', 'X',\n'PAS', 'X',\n'PAT', 'W',\n'PBA', 'X\
',\n'PBB', 'X',\n'PBC', 'X',\n'PBF', 'F',\n'PBG', \
'X',\n'PBI', 'X',\n'PBM', 'X',\n'PBN', 'X',\n'PBP'\
, 'X',\n'PBR', 'X',\n'PBZ', 'X',\n'PC2', 'X',\n'PC\
A', 'E',\n'PCB', 'X',\n'PCD', 'X',\n'PCE', 'X',\n'\
PCG', 'X',\n'PCH', 'X',\n'PCL', 'X',\n'PCM', 'X',\\
n'PCP', 'X',\n'PCR', 'X',\n'PCS', 'X',\n'PCU', 'X'\
,\n'PCV', 'X',\n'PCY', 'X',\n'PD1', 'X',\n'PDA', '\
X',\n'PDC', 'X',\n'PDD', 'A',\n'PDE', 'A',\n'PDI',\
 'X',\n'PDL', 'A',\n'PDN', 'X',\n'PDO', 'X',\n'PDP\
', 'X',\n'PDT', 'X',\n'PDU', 'X',\n'PE2', 'X',\n'P\
E6', 'X',\n'PEA', 'X',\n'PEB', 'X',\n'PEC', 'X',\n\
'PED', 'X',\n'PEE', 'X',\n'PEF', 'X',\n'PEG', 'X',\
\n'PEL', 'X',\n'PEO', 'X',\n'PEP', 'X',\n'PEQ', 'X\
',\n'PER', 'X',\n'PET', 'X',\n'PFB', 'X',\n'PFC', \
'X',\n'PFG', 'X',\n'PFL', 'X',\n'PFM', 'X',\n'PFZ'\
, 'X',\n'PG4', 'X',\n'PG5', 'X',\n'PG6', 'X',\n'PG\
A', 'X',\n'PGC', 'X',\n'PGD', 'X',\n'PGE', 'X',\n'\
PGG', 'G',\n'PGH', 'X',\n'PGL', 'X',\n'PGO', 'X',\\
n'PGP', 'X',\n'PGQ', 'X',\n'PGR', 'X',\n'PGS', 'X'\
,\n'PGU', 'X',\n'PGX', 'X',\n'PGY', 'G',\n'PH1', '\
X',\n'PH2', 'X',\n'PH3', 'X',\n'PHA', 'F',\n'PHB',\
 'X',\n'PHC', 'X',\n'PHD', 'X',\n'PHE', 'F',\n'PHG\
', 'X',\n'PHH', 'X',\n'PHI', 'F',\n'PHL', 'F',\n'P\
HM', 'X',\n'PHN', 'X',\n'PHO', 'X',\n'PHP', 'X',\n\
'PHQ', 'X',\n'PHS', 'H',\n'PHT', 'X',\n'PHW', 'P',\
\n'PHY', 'X',\n'PI1', 'X',\n'PI2', 'X',\n'PI3', 'X\
',\n'PI4', 'X',\n'PI5', 'X',\n'PI6', 'X',\n'PI7', \
'X',\n'PI8', 'X',\n'PI9', 'X',\n'PIA', 'X',\n'PIB'\
, 'X',\n'PIC', 'X',\n'PID', 'X',\n'PIG', 'X',\n'PI\
H', 'X',\n'PIM', 'X',\n'PIN', 'X',\n'PIO', 'X',\n'\
PIP', 'X',\n'PIQ', 'X',\n'PIR', 'X',\n'PIV', 'X',\\
n'PKF', 'X',\n'PL1', 'X',\n'PL9', 'X',\n'PLA', 'D'\
,\n'PLC', 'X',\n'PLE', 'L',\n'PLG', 'G',\n'PLH', '\
X',\n'PLM', 'X',\n'PLP', 'X',\n'PLS', 'S',\n'PLT',\
 'W',\n'PLU', 'L',\n'PLY', 'X',\n'PMA', 'X',\n'PMB\
', 'X',\n'PMC', 'X',\n'PME', 'F',\n'PML', 'X',\n'P\
MM', 'X',\n'PMO', 'X',\n'PMP', 'X',\n'PMS', 'X',\n\
'PMY', 'X',\n'PN2', 'X',\n'PNA', 'X',\n'PNB', 'X',\
\n'PNC', 'G',\n'PND', 'X',\n'PNE', 'A',\n'PNF', 'X\
',\n'PNG', 'X',\n'PNI', 'X',\n'PNL', 'X',\n'PNM', \
'X',\n'PNN', 'X',\n'PNO', 'X',\n'PNP', 'X',\n'PNQ'\
, 'X',\n'PNS', 'X',\n'PNT', 'X',\n'PNU', 'X',\n'PO\
2', 'X',\n'PO4', 'X',\n'POB', 'X',\n'POC', 'X',\n'\
POL', 'X',\n'POM', 'P',\n'PON', 'X',\n'POP', 'X',\\
n'POR', 'X',\n'POS', 'X',\n'PP1', 'X',\n'PP2', 'X'\
,\n'PP3', 'A',\n'PP4', 'X',\n'PP5', 'X',\n'PP6', '\
X',\n'PP7', 'X',\n'PP8', 'N',\n'PP9', 'X',\n'PPB',\
 'X',\n'PPC', 'X',\n'PPD', 'X',\n'PPE', 'E',\n'PPG\
', 'X',\n'PPH', 'F',\n'PPI', 'X',\n'PPJ', 'V',\n'P\
PL', 'X',\n'PPM', 'X',\n'PPN', 'A',\n'PPO', 'X',\n\
'PPP', 'X',\n'PPQ', 'X',\n'PPR', 'X',\n'PPS', 'X',\
\n'PPT', 'X',\n'PPU', 'X',\n'PPX', 'F',\n'PPY', 'X\
',\n'PPZ', 'X',\n'PQ0', 'X',\n'PQN', 'X',\n'PQQ', \
'X',\n'PR1', 'X',\n'PR2', 'X',\n'PR3', 'X',\n'PRA'\
, 'X',\n'PRB', 'X',\n'PRC', 'X',\n'PRD', 'X',\n'PR\
E', 'X',\n'PRF', 'X',\n'PRH', 'X',\n'PRI', 'P',\n'\
PRL', 'X',\n'PRN', 'X',\n'PRO', 'P',\n'PRP', 'X',\\
n'PRR', 'A',\n'PRS', 'P',\n'PRZ', 'X',\n'PS0', 'X'\
,\n'PSA', 'X',\n'PSD', 'X',\n'PSE', 'X',\n'PSF', '\
S',\n'PSG', 'X',\n'PSI', 'X',\n'PSO', 'X',\n'PSQ',\
 'X',\n'PSS', 'X',\n'PST', 'X',\n'PSU', 'X',\n'PT1\
', 'X',\n'PT3', 'X',\n'PTA', 'X',\n'PTC', 'X',\n'P\
TD', 'X',\n'PTE', 'X',\n'PTH', 'Y',\n'PTL', 'X',\n\
'PTM', 'Y',\n'PTN', 'X',\n'PTO', 'X',\n'PTP', 'X',\
\n'PTR', 'Y',\n'PTS', 'X',\n'PTT', 'X',\n'PTU', 'X\
',\n'PTY', 'X',\n'PUA', 'X',\n'PUB', 'X',\n'PUR', \
'X',\n'PUT', 'X',\n'PVA', 'X',\n'PVB', 'X',\n'PVH'\
, 'H',\n'PVL', 'X',\n'PXA', 'X',\n'PXF', 'X',\n'PX\
G', 'X',\n'PXP', 'X',\n'PXY', 'X',\n'PXZ', 'X',\n'\
PY2', 'X',\n'PY4', 'X',\n'PY5', 'X',\n'PY6', 'X',\\
n'PYA', 'A',\n'PYC', 'X',\n'PYD', 'X',\n'PYE', 'X'\
,\n'PYL', 'X',\n'PYM', 'X',\n'PYO', 'X',\n'PYP', '\
X',\n'PYQ', 'X',\n'PYR', 'X',\n'PYS', 'X',\n'PYT',\
 'X',\n'PYX', 'X',\n'PYY', 'X',\n'PYZ', 'X',\n'PZQ\
', 'X',\n'Q82', 'X',\n'QNC', 'X',\n'QND', 'X',\n'Q\
SI', 'Q',\n'QTR', 'X',\n'QUA', 'X',\n'QUE', 'X',\n\
'QUI', 'X',\n'QUO', 'X',\n'R11', 'X',\n'R12', 'X',\
\n'R13', 'X',\n'R18', 'X',\n'R1P', 'X',\n'R56', 'X\
',\n'R5P', 'X',\n'RA2', 'X',\n'RAD', 'X',\n'RAI', \
'X',\n'RAL', 'X',\n'RAM', 'X',\n'RAN', 'X',\n'RAP'\
, 'X',\n'RBF', 'X',\n'RBU', 'X',\n'RCA', 'X',\n'RC\
L', 'X',\n'RCO', 'X',\n'RDC', 'X',\n'RDF', 'W',\n'\
RE9', 'X',\n'REA', 'X',\n'RED', 'K',\n'REO', 'X',\\
n'REP', 'X',\n'RET', 'X',\n'RFA', 'X',\n'RFB', 'X'\
,\n'RFL', 'X',\n'RFP', 'X',\n'RG1', 'X',\n'RGS', '\
X',\n'RH1', 'X',\n'RHA', 'X',\n'RHC', 'X',\n'RHD',\
 'X',\n'RHM', 'X',\n'RHO', 'X',\n'RHQ', 'X',\n'RHS\
', 'X',\n'RIA', 'X',\n'RIB', 'X',\n'RIC', 'X',\n'R\
IF', 'X',\n'RIN', 'X',\n'RIP', 'X',\n'RIT', 'X',\n\
'RMB', 'X',\n'RMN', 'X',\n'RMP', 'X',\n'RNG', 'X',\
\n'RNS', 'X',\n'RNT', 'X',\n'RO2', 'X',\n'RO4', 'X\
',\n'ROC', 'N',\n'ROI', 'X',\n'ROM', 'X',\n'RON', \
'V',\n'ROP', 'X',\n'ROS', 'X',\n'ROX', 'X',\n'RPA'\
, 'X',\n'RPD', 'X',\n'RPH', 'X',\n'RPL', 'X',\n'RP\
P', 'X',\n'RPR', 'X',\n'RPX', 'X',\n'RQ3', 'X',\n'\
RR1', 'X',\n'RR6', 'X',\n'RRS', 'X',\n'RS1', 'X',\\
n'RS2', 'X',\n'RS7', 'X',\n'RSS', 'X',\n'RTA', 'X'\
,\n'RTB', 'X',\n'RTC', 'X',\n'RTL', 'X',\n'RUB', '\
X',\n'RUN', 'X',\n'RWJ', 'X',\n'RXP', 'X',\n'S02',\
 'X',\n'S11', 'X',\n'S1H', 'S',\n'S27', 'X',\n'S2C\
', 'C',\n'S3P', 'X',\n'S4U', 'X',\n'S57', 'X',\n'S\
58', 'X',\n'S5H', 'X',\n'S6G', 'X',\n'S80', 'X',\n\
'SAA', 'X',\n'SAB', 'X',\n'SAC', 'S',\n'SAD', 'X',\
\n'SAE', 'X',\n'SAF', 'X',\n'SAH', 'C',\n'SAI', 'C\
',\n'SAL', 'X',\n'SAM', 'M',\n'SAN', 'X',\n'SAP', \
'X',\n'SAR', 'X',\n'SAS', 'X',\n'SB1', 'X',\n'SB2'\
, 'X',\n'SB3', 'X',\n'SB4', 'X',\n'SB5', 'X',\n'SB\
6', 'X',\n'SBA', 'L',\n'SBB', 'X',\n'SBD', 'A',\n'\
SBI', 'X',\n'SBL', 'A',\n'SBN', 'X',\n'SBO', 'X',\\
n'SBR', 'X',\n'SBS', 'X',\n'SBT', 'X',\n'SBU', 'X'\
,\n'SBX', 'X',\n'SC4', 'X',\n'SCA', 'X',\n'SCC', '\
X',\n'SCD', 'X',\n'SCH', 'C',\n'SCI', 'X',\n'SCL',\
 'X',\n'SCM', 'X',\n'SCN', 'X',\n'SCO', 'X',\n'SCP\
', 'S',\n'SCR', 'X',\n'SCS', 'X',\n'SCV', 'C',\n'S\
CY', 'C',\n'SD8', 'X',\n'SDK', 'X',\n'SDZ', 'X',\n\
'SE4', 'X',\n'SEA', 'X',\n'SEB', 'S',\n'SEC', 'X',\
\n'SEG', 'A',\n'SEI', 'X',\n'SEL', 'S',\n'SEM', 'X\
',\n'SEO', 'X',\n'SEP', 'S',\n'SER', 'S',\n'SES', \
'X',\n'SET', 'S',\n'SEU', 'X',\n'SF4', 'X',\n'SFG'\
, 'X',\n'SFN', 'X',\n'SFO', 'X',\n'SGA', 'X',\n'SG\
C', 'X',\n'SGL', 'X',\n'SGM', 'X',\n'SGN', 'X',\n'\
SGP', 'X',\n'SHA', 'X',\n'SHC', 'X',\n'SHF', 'X',\\
n'SHH', 'X',\n'SHP', 'G',\n'SHR', 'E',\n'SHT', 'T'\
,\n'SHU', 'X',\n'SI2', 'X',\n'SIA', 'X',\n'SIF', '\
X',\n'SIG', 'X',\n'SIH', 'X',\n'SIM', 'X',\n'SIN',\
 'X',\n'SKD', 'X',\n'SKF', 'X',\n'SLB', 'X',\n'SLE\
', 'X',\n'SLZ', 'K',\n'SMA', 'X',\n'SMC', 'C',\n'S\
ME', 'M',\n'SML', 'X',\n'SMM', 'M',\n'SMN', 'X',\n\
'SMP', 'X',\n'SMS', 'X',\n'SN1', 'X',\n'SN6', 'X',\
\n'SN7', 'X',\n'SNC', 'C',\n'SNN', 'X',\n'SNP', 'X\
',\n'SO1', 'X',\n'SO2', 'X',\n'SO3', 'X',\n'SO4', \
'X',\n'SOA', 'X',\n'SOC', 'C',\n'SOM', 'X',\n'SOR'\
, 'X',\n'SOT', 'X',\n'SOX', 'X',\n'SPA', 'X',\n'SP\
B', 'X',\n'SPC', 'X',\n'SPD', 'X',\n'SPE', 'X',\n'\
SPG', 'X',\n'SPH', 'X',\n'SPI', 'X',\n'SPK', 'X',\\
n'SPM', 'X',\n'SPN', 'X',\n'SPO', 'X',\n'SPP', 'X'\
,\n'SPS', 'X',\n'SPY', 'X',\n'SQU', 'X',\n'SRA', '\
X',\n'SRB', 'X',\n'SRD', 'X',\n'SRL', 'X',\n'SRM',\
 'X',\n'SRS', 'X',\n'SRY', 'X',\n'SSA', 'X',\n'SSB\
', 'X',\n'SSG', 'X',\n'SSP', 'X',\n'ST1', 'X',\n'S\
T2', 'X',\n'ST3', 'X',\n'ST4', 'X',\n'ST5', 'X',\n\
'ST6', 'X',\n'STA', 'X',\n'STB', 'X',\n'STE', 'X',\
\n'STG', 'X',\n'STI', 'X',\n'STL', 'X',\n'STN', 'X\
',\n'STO', 'X',\n'STP', 'X',\n'STR', 'X',\n'STU', \
'X',\n'STY', 'Y',\n'SU1', 'X',\n'SU2', 'X',\n'SUC'\
, 'X',\n'SUI', 'X',\n'SUL', 'X',\n'SUR', 'X',\n'SV\
A', 'S',\n'SWA', 'X',\n'T16', 'X',\n'T19', 'X',\n'\
T23', 'X',\n'T29', 'X',\n'T33', 'X',\n'T3P', 'X',\\
n'T42', 'A',\n'T44', 'X',\n'T5A', 'X',\n'T6A', 'T'\
,\n'T6P', 'X',\n'T80', 'X',\n'T87', 'X',\n'TA1', '\
X',\n'TAA', 'X',\n'TAB', 'X',\n'TAC', 'X',\n'TAD',\
 'X',\n'TAF', 'X',\n'TAM', 'X',\n'TAP', 'X',\n'TAR\
', 'X',\n'TAS', 'X',\n'TAU', 'X',\n'TAX', 'X',\n'T\
AZ', 'X',\n'TB9', 'X',\n'TBA', 'X',\n'TBD', 'X',\n\
'TBG', 'G',\n'TBH', 'X',\n'TBM', 'T',\n'TBO', 'X',\
\n'TBP', 'X',\n'TBR', 'X',\n'TBS', 'X',\n'TBT', 'X\
',\n'TBU', 'X',\n'TBZ', 'X',\n'TC4', 'X',\n'TCA', \
'X',\n'TCB', 'X',\n'TCH', 'X',\n'TCK', 'X',\n'TCL'\
, 'X',\n'TCM', 'X',\n'TCN', 'X',\n'TCP', 'X',\n'TC\
R', 'W',\n'TCS', 'X',\n'TCZ', 'X',\n'TDA', 'X',\n'\
TDB', 'X',\n'TDG', 'X',\n'TDP', 'X',\n'TDR', 'X',\\
n'TDX', 'X',\n'TEA', 'X',\n'TEM', 'X',\n'TEN', 'X'\
,\n'TEO', 'X',\n'TEP', 'X',\n'TER', 'X',\n'TES', '\
X',\n'TET', 'X',\n'TFA', 'X',\n'TFB', 'X',\n'TFH',\
 'X',\n'TFI', 'X',\n'TFK', 'X',\n'TFP', 'X',\n'THA\
', 'X',\n'THB', 'X',\n'THC', 'T',\n'THD', 'X',\n'T\
HE', 'X',\n'THF', 'X',\n'THJ', 'X',\n'THK', 'X',\n\
'THM', 'X',\n'THN', 'X',\n'THO', 'T',\n'THP', 'X',\
\n'THQ', 'X',\n'THR', 'T',\n'THS', 'X',\n'THT', 'X\
',\n'THU', 'X',\n'THX', 'X',\n'THZ', 'X',\n'TI1', \
'X',\n'TI2', 'X',\n'TI3', 'P',\n'TIA', 'X',\n'TIH'\
, 'A',\n'TK4', 'X',\n'TLA', 'X',\n'TLC', 'X',\n'TL\
M', 'X',\n'TLN', 'X',\n'TLX', 'X',\n'TM5', 'X',\n'\
TM6', 'X',\n'TMA', 'X',\n'TMB', 'T',\n'TMC', 'X',\\
n'TMD', 'T',\n'TME', 'X',\n'TMF', 'X',\n'TML', 'K'\
,\n'TMM', 'X',\n'TMN', 'X',\n'TMP', 'X',\n'TMQ', '\
X',\n'TMR', 'X',\n'TMT', 'X',\n'TMZ', 'X',\n'TNB',\
 'C',\n'TND', 'X',\n'TNK', 'X',\n'TNP', 'X',\n'TNT\
', 'X',\n'TOA', 'X',\n'TOB', 'X',\n'TOC', 'X',\n'T\
OL', 'X',\n'TOP', 'X',\n'TOS', 'X',\n'TOT', 'X',\n\
'TP1', 'G',\n'TP2', 'P',\n'TP3', 'E',\n'TP4', 'E',\
\n'TP7', 'T',\n'TPA', 'X',\n'TPE', 'X',\n'TPF', 'X\
',\n'TPI', 'X',\n'TPL', 'W',\n'TPM', 'X',\n'TPN', \
'G',\n'TPO', 'T',\n'TPP', 'X',\n'TPQ', 'A',\n'TPR'\
, 'P',\n'TPS', 'X',\n'TPT', 'X',\n'TPV', 'X',\n'TP\
X', 'X',\n'TPY', 'X',\n'TQ3', 'X',\n'TQ4', 'X',\n'\
TQ5', 'X',\n'TQ6', 'X',\n'TR1', 'X',\n'TRA', 'X',\\
n'TRB', 'X',\n'TRC', 'X',\n'TRD', 'X',\n'TRE', 'X'\
,\n'TRF', 'W',\n'TRG', 'K',\n'TRH', 'X',\n'TRI', '\
X',\n'TRJ', 'X',\n'TRM', 'X',\n'TRN', 'W',\n'TRO',\
 'W',\n'TRP', 'W',\n'TRQ', 'X',\n'TRS', 'X',\n'TRX\
', 'W',\n'TRZ', 'X',\n'TS2', 'X',\n'TS3', 'X',\n'T\
S4', 'X',\n'TS5', 'X',\n'TSA', 'X',\n'TSB', 'X',\n\
'TSI', 'X',\n'TSM', 'X',\n'TSN', 'X',\n'TSP', 'X',\
\n'TSU', 'X',\n'TTA', 'X',\n'TTE', 'X',\n'TTN', 'X\
',\n'TTO', 'X',\n'TTP', 'X',\n'TTX', 'X',\n'TXL', \
'X',\n'TYA', 'Y',\n'TYB', 'Y',\n'TYD', 'X',\n'TYI'\
, 'Y',\n'TYL', 'X',\n'TYM', 'W',\n'TYN', 'Y',\n'TY\
Q', 'Y',\n'TYR', 'Y',\n'TYS', 'Y',\n'TYV', 'X',\n'\
TYY', 'A',\n'TZB', 'X',\n'TZC', 'X',\n'TZE', 'X',\\
n'TZL', 'X',\n'TZO', 'X',\n'TZP', 'X',\n'U01', 'X'\
,\n'U02', 'X',\n'U03', 'X',\n'U04', 'X',\n'U05', '\
X',\n'U0E', 'X',\n'U10', 'X',\n'U18', 'X',\n'U2G',\
 'X',\n'U3P', 'X',\n'U49', 'X',\n'U55', 'X',\n'U5P\
', 'X',\n'U66', 'X',\n'U89', 'X',\n'U8U', 'X',\n'U\
AA', 'X',\n'UAG', 'A',\n'UAP', 'X',\n'UAR', 'X',\n\
'UC1', 'X',\n'UC2', 'X',\n'UC3', 'X',\n'UC4', 'X',\
\n'UD1', 'X',\n'UD2', 'X',\n'UDP', 'X',\n'UDX', 'X\
',\n'UFG', 'X',\n'UFM', 'X',\n'UFP', 'X',\n'UGA', \
'X',\n'UIN', 'X',\n'UKP', 'A',\n'UM3', 'X',\n'UMA'\
, 'A',\n'UMG', 'X',\n'UMP', 'X',\n'UNA', 'X',\n'UN\
D', 'X',\n'UNI', 'X',\n'UNK', 'X',\n'UNN', 'X',\n'\
UNX', 'X',\n'UP5', 'X',\n'UP6', 'X',\n'UPA', 'X',\\
n'UPF', 'X',\n'UPG', 'X',\n'UPP', 'X',\n'UQ1', 'X'\
,\n'UQ2', 'X',\n'UQ6', 'X',\n'UR2', 'X',\n'URA', '\
X',\n'URE', 'X',\n'URF', 'X',\n'URI', 'X',\n'URS',\
 'X',\n'UTP', 'X',\n'UVC', 'X',\n'UVW', 'X',\n'V35\
', 'X',\n'V36', 'X',\n'V4O', 'X',\n'V7O', 'X',\n'V\
AA', 'V',\n'VAC', 'X',\n'VAD', 'V',\n'VAF', 'V',\n\
'VAG', 'X',\n'VAL', 'V',\n'VAN', 'X',\n'VAS', 'X',\
\n'VAX', 'X',\n'VDX', 'X',\n'VDY', 'X',\n'VG1', 'X\
',\n'VIB', 'X',\n'VIR', 'X',\n'VIT', 'X',\n'VK3', \
'X',\n'VO3', 'X',\n'VO4', 'X',\n'VS1', 'F',\n'VS2'\
, 'F',\n'VS3', 'F',\n'VS4', 'F',\n'VXA', 'X',\n'W0\
1', 'X',\n'W02', 'X',\n'W03', 'X',\n'W11', 'X',\n'\
W33', 'X',\n'W35', 'X',\n'W42', 'X',\n'W43', 'X',\\
n'W54', 'X',\n'W56', 'X',\n'W59', 'X',\n'W71', 'X'\
,\n'W84', 'X',\n'W8R', 'X',\n'W91', 'X',\n'WAY', '\
X',\n'WCC', 'X',\n'WO2', 'X',\n'WO4', 'X',\n'WRB',\
 'X',\n'WRR', 'X',\n'WRS', 'X',\n'WW7', 'X',\n'X2F\
', 'X',\n'X7O', 'X',\n'XAA', 'X',\n'XAN', 'X',\n'X\
AO', 'X',\n'XBB', 'X',\n'XBP', 'X',\n'XDN', 'X',\n\
'XDP', 'X',\n'XIF', 'X',\n'XIM', 'X',\n'XK2', 'X',\
\n'XL1', 'X',\n'XLS', 'X',\n'XMP', 'X',\n'XN1', 'X\
',\n'XN2', 'X',\n'XN3', 'X',\n'XUL', 'X',\n'XV6', \
'X',\n'XYD', 'X',\n'XYH', 'X',\n'XYL', 'X',\n'XYP'\
, 'X',\n'XYS', 'X',\n'YOF', 'Y',\n'YRR', 'X',\n'YT\
3', 'X',\n'YZ9', 'X',\n'Z34', 'G',\n'Z5A', 'X',\n'\
ZAF', 'X',\n'ZAP', 'X',\n'ZEB', 'X',\n'ZEN', 'X',\\
n'ZES', 'X',\n'ZID', 'X',\n'ZMR', 'X',\n'ZN3', 'X'\
,\n'ZNH', 'X',\n'ZNO', 'X',\n'ZO3', 'X',\n'ZPR', '\
P',\n'ZRA', 'A',\n'ZST', 'X',\n'ZYA', 'A',\n\n\n'A\
SN','N');\n} \n\n\nsub file2head\n      {\n	my $fi\
le = shift;\n	my $size = shift;\n	my $f= new FileH\
andle;\n	my $line;\n	open ($f,$file);\n	read ($f,$\
line, $size);\n	close ($f);\n	return $line;\n     \
 }\nsub file2tail\n      {\n	my $file = shift;\n	m\
y $size = shift;\n	my $f= new FileHandle;\n	my $li\
ne;\n	\n	open ($f,$file);\n	seek ($f,$size*-1, 2);\
\n	read ($f,$line, $size);\n	close ($f);\n	return \
$line;\n      }\n\n\nsub vtmpnam\n      {\n	my $r=\
rand(100000);\n	my $f=\"file.$r.$$\";\n	while (-e \
$f)\n	  {\n	    $f=vtmpnam();\n	  }\n	push (@TMPFI\
LE_LIST, $f);\n	return $f;\n      }\n\nsub myexit\\
n  {\n    my $code=@_[0];\n    if ($CLEAN_EXIT_STA\
RTED==1){return;}\n    else {$CLEAN_EXIT_STARTED=1\
;}\n    ### ONLY BARE EXIT\n    exit ($code);\n  }\
\nsub set_error_lock\n    {\n      my $name = shif\
t;\n      my $pid=$$;\n\n      \n      &lock4tc ($\
$,\"LERROR\", \"LSET\", \"$$ -- ERROR: $name $PROG\
RAM\\n\");\n      return;\n    }\nsub set_lock\n  \
{\n    my $pid=shift;\n    my $msg= shift;\n    my\
 $p=getppid();\n    &lock4tc ($pid,\"LLOCK\",\"LRE\
SET\",\"$p$msg\\n\");\n  }\nsub unset_lock\n   {\n\
     \n    my $pid=shift;\n    &lock4tc ($pid,\"LL\
OCK\",\"LRELEASE\",\"\");\n  }\nsub shift_lock\n  \
{\n    my $from=shift;\n    my $to=shift;\n    my \
$from_type=shift;\n    my $to_type=shift;\n    my \
$action=shift;\n    my $msg;\n    \n    if (!&lock\
4tc($from, $from_type, \"LCHECK\", \"\")){return 0\
;}\n    $msg=&lock4tc ($from, $from_type, \"LREAD\\
", \"\");\n    &lock4tc ($from, $from_type,\"LRELE\
ASE\", $msg);\n    &lock4tc ($to, $to_type, $actio\
n, $msg);\n    return;\n  }\nsub isshellpid\n  {\n\
    my $p=shift;\n    if (!lock4tc ($p, \"LLOCK\",\
 \"LCHECK\")){return 0;}\n    else\n      {\n	my $\
c=lock4tc($p, \"LLOCK\", \"LREAD\");\n	if ( $c=~/-\
SHELL-/){return 1;}\n      }\n    return 0;\n  }\n\
sub isrootpid\n  {\n    if(lock4tc (getppid(), \"L\
LOCK\", \"LCHECK\")){return 0;}\n    else {return \
1;}\n  }\nsub lock4tc\n	{\n	  my ($pid,$type,$acti\
on,$value)=@_;\n	  my $fname;\n	  my $host=hostnam\
e;\n	  \n	  if ($type eq \"LLOCK\"){$fname=\"$LOCK\
DIR/.$pid.$host.lock4tcoffee\";}\n	  elsif ( $type\
 eq \"LERROR\"){ $fname=\"$LOCKDIR/.$pid.$host.err\
or4tcoffee\";}\n	  elsif ( $type eq \"LWARNING\"){\
 $fname=\"$LOCKDIR/.$pid.$host.warning4tcoffee\";}\
\n	  \n	  if ($debug_lock)\n	    {\n	      print S\
TDERR \"\\n\\t---lock4tc(tcg): $action => $fname =\
>$value (RD: $LOCKDIR)\\n\";\n	    }\n\n	  if    (\
$action eq \"LCHECK\") {return -e $fname;}\n	  els\
if ($action eq \"LREAD\"){return file2string($fnam\
e);}\n	  elsif ($action eq \"LSET\") {return strin\
g2file ($value, $fname, \">>\");}\n	  elsif ($acti\
on eq \"LRESET\") {return string2file ($value, $fn\
ame, \">\");}\n	  elsif ($action eq \"LRELEASE\") \
\n	    {\n	      if ( $debug_lock)\n		{\n		  my $g\
=new FileHandle;\n		  open ($g, \">>$fname\");\n		\
  print $g \"\\nDestroyed by $$\\n\";\n		  close (\
$g);\n		  safe_system (\"mv $fname $fname.old\");\\
n		}\n	      else\n		{\n		  unlink ($fname);\n		}\\
n	    }\n	  return \"\";\n	}\n	\nsub file2string\n\
	{\n	  my $file=@_[0];\n	  my $f=new FileHandle;\n\
	  my $r;\n	  open ($f, \"$file\");\n	  while (<$f\
>){$r.=$_;}\n	  close ($f);\n	  return $r;\n	}\nsu\
b string2file \n    {\n    my ($s,$file,$mode)=@_;\
\n    my $f=new FileHandle;\n    \n    open ($f, \\
"$mode$file\");\n    print $f  \"$s\";\n    close \
($f);\n  }\n\nBEGIN\n    {\n      srand;\n    \n  \
    $SIG{'SIGUP'}='signal_cleanup';\n      $SIG{'S\
IGINT'}='signal_cleanup';\n      $SIG{'SIGQUIT'}='\
signal_cleanup';\n      $SIG{'SIGILL'}='signal_cle\
anup';\n      $SIG{'SIGTRAP'}='signal_cleanup';\n \
     $SIG{'SIGABRT'}='signal_cleanup';\n      $SIG\
{'SIGEMT'}='signal_cleanup';\n      $SIG{'SIGFPE'}\
='signal_cleanup';\n      \n      $SIG{'SIGKILL'}=\
'signal_cleanup';\n      $SIG{'SIGPIPE'}='signal_c\
leanup';\n      $SIG{'SIGSTOP'}='signal_cleanup';\\
n      $SIG{'SIGTTIN'}='signal_cleanup';\n      $S\
IG{'SIGXFSZ'}='signal_cleanup';\n      $SIG{'SIGIN\
FO'}='signal_cleanup';\n      \n      $SIG{'SIGBUS\
'}='signal_cleanup';\n      $SIG{'SIGALRM'}='signa\
l_cleanup';\n      $SIG{'SIGTSTP'}='signal_cleanup\
';\n      $SIG{'SIGTTOU'}='signal_cleanup';\n     \
 $SIG{'SIGVTALRM'}='signal_cleanup';\n      $SIG{'\
SIGUSR1'}='signal_cleanup';\n\n\n      $SIG{'SIGSE\
GV'}='signal_cleanup';\n      $SIG{'SIGTERM'}='sig\
nal_cleanup';\n      $SIG{'SIGCONT'}='signal_clean\
up';\n      $SIG{'SIGIO'}='signal_cleanup';\n     \
 $SIG{'SIGPROF'}='signal_cleanup';\n      $SIG{'SI\
GUSR2'}='signal_cleanup';\n\n      $SIG{'SIGSYS'}=\
'signal_cleanup';\n      $SIG{'SIGURG'}='signal_cl\
eanup';\n      $SIG{'SIGCHLD'}='signal_cleanup';\n\
      $SIG{'SIGXCPU'}='signal_cleanup';\n      $SI\
G{'SIGWINCH'}='signal_cleanup';\n      \n      $SI\
G{'INT'}='signal_cleanup';\n      $SIG{'TERM'}='si\
gnal_cleanup';\n      $SIG{'KILL'}='signal_cleanup\
';\n      $SIG{'QUIT'}='signal_cleanup';\n      \n\
      our $debug_lock=$ENV{\"DEBUG_LOCK\"};\n     \
 \n      \n      \n      \n      foreach my $a (@A\
RGV){$CL.=\" $a\";}\n      if ( $debug_lock ){prin\
t STDERR \"\\n\\n\\n********** START PG: $PROGRAM \
*************\\n\";}\n      if ( $debug_lock ){pri\
nt STDERR \"\\n\\n\\n**********(tcg) LOCKDIR: $LOC\
KDIR $$ *************\\n\";}\n      if ( $debug_lo\
ck ){print STDERR \"\\n --- $$ -- $CL\\n\";}\n    \
  \n	     \n      \n      \n    }\nsub flush_error\
\n  {\n    my $msg=shift;\n    return add_error ($\
EXIT_FAILURE,$$, $$,getppid(), $msg, $CL);\n  }\ns\
ub add_error \n  {\n    my $code=shift;\n    my $r\
pid=shift;\n    my $pid=shift;\n    my $ppid=shift\
;\n    my $type=shift;\n    my $com=shift;\n    \n\
    $ERROR_DONE=1;\n    lock4tc ($rpid, \"LERROR\"\
,\"LSET\",\"$pid -- ERROR: $type\\n\");\n    lock4\
tc ($$, \"LERROR\",\"LSET\", \"$pid -- COM: $com\\\
n\");\n    lock4tc ($$, \"LERROR\",\"LSET\", \"$pi\
d -- STACK: $ppid -> $pid\\n\");\n   \n    return \
$code;\n  }\nsub add_warning \n  {\n    my $rpid=s\
hift;\n    my $pid =shift;\n    my $command=shift;\
\n    my $msg=\"$$ -- WARNING: $command\\n\";\n   \
 print STDERR \"$msg\";\n    lock4tc ($$, \"LWARNI\
NG\", \"LSET\", $msg);\n  }\n\nsub signal_cleanup\\
n  {\n    print dtderr \"\\n**** $$ (tcg) was kill\
ed\\n\";\n    &cleanup;\n    exit ($EXIT_FAILURE);\
\n  }\nsub clean_dir\n  {\n    my $dir=@_[0];\n   \
 if ( !-d $dir){return ;}\n    elsif (!($dir=~/tmp\
/)){return ;}#safety check 1\n    elsif (($dir=~/\\
\*/)){return ;}#safety check 2\n    else\n      {\\
n	`rm -rf $dir`;\n      }\n    return;\n  }\nsub c\
leanup\n  {\n    #print stderr \"\\n----tc: $$ Kil\
ls $PIDCHILD\\n\";\n    #kill (SIGTERM,$PIDCHILD);\
\n    my $p=getppid();\n    $CLEAN_EXIT_STARTED=1;\
\n    \n    \n    \n    if (&lock4tc($$,\"LERROR\"\
, \"LCHECK\", \"\"))\n      {\n	my $ppid=getppid()\
;\n	if (!$ERROR_DONE) \n	  {\n	    &lock4tc($$,\"L\
ERROR\", \"LSET\", \"$$ -- STACK: $p -> $$\\n\");\\
n	    &lock4tc($$,\"LERROR\", \"LSET\", \"$$ -- CO\
M: $CL\\n\");\n	  }\n      }\n    my $warning=&loc\
k4tc($$, \"LWARNING\", \"LREAD\", \"\");\n    my $\
error=&lock4tc($$,  \"LERROR\", \"LREAD\", \"\");\\
n    #release error and warning lock if root\n    \
\n    if (isrootpid() && ($warning || $error) )\n \
     {\n	\n	print STDERR \"**************** Summar\
y *************\\n$error\\n$warning\\n\";\n\n	&loc\
k4tc($$,\"LERROR\",\"RELEASE\",\"\");\n	&lock4tc($\
$,\"LWARNING\",\"RELEASE\",\"\");\n      } \n    \\
n    \n    foreach my $f (@TMPFILE_LIST)\n      {\\
n	if (-e $f){unlink ($f);} \n      }\n    foreach \
my $d (@TMPDIR_LIST)\n      {\n	clean_dir ($d);\n \
     }\n    #No More Lock Release\n    #&lock4tc($\
$,\"LLOCK\",\"LRELEASE\",\"\"); #release lock \n\n\
    if ( $debug_lock ){print STDERR \"\\n\\n\\n***\
******* END PG: $PROGRAM ($$) *************\\n\";}\
\n    if ( $debug_lock ){print STDERR \"\\n\\n\\n*\
*********(tcg) LOCKDIR: $LOCKDIR $$ *************\\
\n\";}\n  }\nEND \n  {\n    \n    &cleanup();\n  }\
\n   \n\nsub safe_system \n{\n  my $com=shift;\n  \
my $ntry=shift;\n  my $ctry=shift;\n  my $pid;\n  \
my $status;\n  my $ppid=getppid();\n  if ($com eq \
\"\"){return 1;}\n  \n  \n\n  if (($pid = fork ())\
 < 0){return (-1);}\n  if ($pid == 0)\n    {\n    \
  set_lock($$, \" -SHELL- $com (tcg)\");\n      ex\
ec ($com);\n    }\n  else\n    {\n      lock4tc ($\
$, \"LLOCK\", \"LSET\", \"$pid\\n\");#update paren\
t\n      $PIDCHILD=$pid;\n    }\n  if ($debug_lock\
){printf STDERR \"\\n\\t .... safe_system (fasta_s\
eq2hmm)  p: $$ c: $pid COM: $com\\n\";}\n\n  waitp\
id ($pid,WTERMSIG);\n\n  shift_lock ($pid,$$, \"LW\
ARNING\",\"LWARNING\", \"LSET\");\n\n  if ($? == $\
EXIT_FAILURE || lock4tc($pid, \"LERROR\", \"LCHECK\
\", \"\"))\n    {\n      if ($ntry && $ctry <$ntry\
)\n	{\n	  add_warning ($$,$$,\"$com failed [retry:\
 $ctry]\");\n	  lock4tc ($pid, \"LRELEASE\", \"LER\
ROR\", \"\");\n	  return safe_system ($com, $ntry,\
 ++$ctry);\n	}\n      elsif ($ntry == -1)\n	{\n	  \
if (!shift_lock ($pid, $$, \"LERROR\", \"LWARNING\\
", \"LSET\"))\n	    {\n	      add_warning ($$,$$,\\
"$com failed\");\n	    }\n	  else\n	    {\n	      \
lock4tc ($pid, \"LRELEASE\", \"LERROR\", \"\");\n	\
    }\n	  return $?;}\n      else\n	{\n	  if (!shi\
ft_lock ($pid,$$, \"LERROR\",\"LERROR\", \"LSET\")\
)\n	    {\n	      myexit(add_error ($EXIT_FAILURE,\
$$,$pid,getppid(), \"UNSPECIFIED system\", $com));\
\n	    }\n	}\n    }\n  return $?;\n}\n\nsub check_\
configuration \n    {\n      my @l=@_;\n      my $\
v;\n      foreach my $p (@l)\n	{\n	  \n	  if   ( $\
p eq \"EMAIL\")\n	    { \n	      if ( !($EMAIL=~/@\
/))\n		{\n		add_warning($$,$$,\"Could Not Use EMAI\
L\");\n		myexit(add_error ($EXIT_FAILURE,$$,$$,get\
ppid(),\"EMAIL\",\"$CL\"));\n	      }\n	    }\n	  \
elsif( $p eq \"INTERNET\")\n	    {\n	      if ( !&\
check_internet_connection())\n		{\n		  myexit(add_\
error ($EXIT_FAILURE,$$,$$,getppid(),\"INTERNET\",\
\"$CL\"));\n		}\n	    }\n	  elsif( $p eq \"wget\")\
\n	    {\n	      if (!&pg_is_installed (\"wget\") \
&& !&pg_is_installed (\"curl\"))\n		{\n		  myexit(\
add_error ($EXIT_FAILURE,$$,$$,getppid(),\"PG_NOT_\
INSTALLED:wget\",\"$CL\"));\n		}\n	    }\n	  elsif\
( !(&pg_is_installed ($p)))\n	    {\n	      myexit\
(add_error ($EXIT_FAILURE,$$,$$,getppid(),\"PG_NOT\
_INSTALLED:$p\",\"$CL\"));\n	    }\n	}\n      retu\
rn 1;\n    }\nsub pg_is_installed\n  {\n    my @ml\
=@_;\n    my $r, $p, $m;\n    my $supported=0;\n  \
  \n    my $p=shift (@ml);\n    if ($p=~/::/)\n   \
   {\n	if (safe_system (\"perl -M$p -e 1\")==$EXIT\
_SUCCESS){return 1;}\n	else {return 0;}\n      }\n\
    else\n      {\n	$r=`which $p 2>/dev/null`;\n	i\
f ($r eq \"\"){return 0;}\n	else {return 1;}\n    \
  }\n  }\n\n\n\nsub check_internet_connection\n  {\
\n    my $internet;\n    my $tmp;\n    &check_conf\
iguration ( \"wget\"); \n    \n    $tmp=&vtmpnam (\
);\n    \n    if     (&pg_is_installed    (\"wget\\
")){`wget www.google.com -O$tmp >/dev/null 2>/dev/\
null`;}\n    elsif  (&pg_is_installed    (\"curl\"\
)){`curl www.google.com -o$tmp >/dev/null 2>/dev/n\
ull`;}\n    \n    if ( !-e $tmp || -s $tmp < 10){$\
internet=0;}\n    else {$internet=1;}\n    if (-e \
$tmp){unlink $tmp;}\n\n    return $internet;\n  }\\
nsub check_pg_is_installed\n  {\n    my @ml=@_;\n \
   my $r=&pg_is_installed (@ml);\n    if (!$r && $\
p=~/::/)\n      {\n	print STDERR \"\\nYou Must Ins\
tall the perl package $p on your system.\\nRUN:\\n\
\\tsudo perl -MCPAN -e 'install $pg'\\n\";\n      \
}\n    elsif (!$r)\n      {\n	myexit(flush_error(\\
"\\nProgram $p Supported but Not Installed on your\
 system\"));\n      }\n    else\n      {\n	return \
1;\n      }\n  }\n\n\n","use Cwd;\nuse Env;\nuse F\
ile::Path;\nuse FileHandle;\nuse strict;\n\n\nour \
(%MODE, %PG, %ENV_SET, %SUPPORTED_OS);\n\n\nour $E\
XIT_SUCCESS=0;\nour $EXIT_FAILURE=1;\nour $INTERNE\
T=0;\n\nour $CP=\"cp \"; #was causing a crash on M\
acOSX\nour $SILENT=\">/dev/null 2>/dev/null\";\nou\
r $WEB_BASE=\"http://www.tcoffee.org\";\nour $TCLI\
NKDB_ADDRESS=\"$WEB_BASE/Resources/tclinkdb.txt\";\
\nour $OS=get_os();\nour $ROOT=&get_root();\nour $\
CD=cwd();\nour $CDIR=$CD;\nour $HOME=$ENV{'HOME'};\
\n\nour $OSNAME=$ENV{'OSNAME'};\nour $OSARCH=$ENV{\
'OSARCH'};\nour $REPO_ROOT=\"\";\n\nour $CXX=\"g++\
\";\nour $CXXFLAGS=\"\";\n\nour $CPP=\"g++\";\nour\
 $CPPFLAGS=\"\";\n\nour $CC=\"gcc\";\nour $CFLAGS=\
\"\";\n\nour $FC=\"f77\";\nour $FFLAGS=\"\";\n\nmy\
 $install=\"all\";\nmy $default_update_action=\"no\
_update\";\nmy @required_applications=(\"wget_OR_c\
url\");\nmy @smode=(\"all\", \"clean\", \"install\\
");\n\n&initialize_PG();\n\nmy $cl=join( \" \", @A\
RGV);\nif ($#ARGV==-1 || ($cl=~/-h/) ||($cl=~/-H/)\
 )\n  {\n     print \"\\n!!!!!!! ./install  t_coff\
ee             --> installs t_coffee only\";\n    \
 print \"\\n!!!!!!! ./install  all                \
  --> installs all the modes [mcoffee, expresso, p\
sicoffee,rcoffee..]\";\n     print \"\\n!!!!!!! ./\
install  [mcoffee|rcoffee|..] --> installs the spe\
cified mode\";\n     print \"\\n!!!!!!! ./install \
 -h                   --> print usage\\n\\n\";\n  \
   if ( $#ARGV==-1){exit ($EXIT_FAILURE);}\n   }\n\
     \nif (($cl=~/-h/) ||($cl=~/-H/) )\n  {\n    m\
y $m;\n    print \"\\n\\n!!!!!!! advanced mode\\n\\
";\n    foreach $m ((keys (%MODE)),@smode)\n      \
{\n	print \"!!!!!!!       ./install $m\\n\";\n    \
  }\n    \n    print \"!!!!!!! ./install [target:p\
ackage|mode|] [-update|-force|-exec=dir|-dis=dir|-\
root|-tclinkdb=file|-] [CC=|FCC=|CXX=|CFLAGS=|CXXF\
LAGS=]\\n\";\n    print \"!!!!!!! ./install clean \
   [removes all executables]\\n\";\n    print \"!!\
!!!!! ./install [optional:target] -update         \
      [updates package already installed]\\n\";\n \
   print \"!!!!!!! ./install [optional:target] -fo\
rce                [Forces recompilation over ever\
ything]\\n\";\n    \n    print \"!!!!!!! ./install\
 [optional:target] -root                 [You are \
running as root]\\n\";\n    print \"!!!!!!! ./inst\
all [optional:target] -exec=/foo/bar/       [addre\
ss for the T-Coffee executable]\\n\";\n    print \\
"!!!!!!! ./install [optional:target] -dis=/foo/bar\
/        [Address where distributions should be st\
ored]\\n\";\n    print \"!!!!!!! ./install [option\
al:target] -tclinkdb=foo|update  [file containing \
all the packages to be installed]\\n\";\n    print\
 \"!!!!!!! ./install [optional:target] -clean     \
           [clean everything]\\n\";\n    print \"!\
!!!!!! ./install [optional:target] -plugins       \
       [plugins directory]\\n\";\n    print \"!!!!\
!!! ./install [optional:target] -repo=/path/to/rep\
o   [binaries repository root directory]\\n\";\n  \
  print \"!!!!!!! mode:\";\n    foreach $m (keys(%\
MODE)){print \"$m \";}\n    print \"\\n\";\n    pr\
int \"!!!!!!! Packages:\";\n    foreach $m (keys (\
%PG)){print \"$m \";}\n    print \"\\n\";\n    \n \
   print \"\\n\\n\";\n    exit ($EXIT_FAILURE);\n \
 }\n\n\n\nmy (@argl)=($cl=~/(\\S+=[^=]+)\\s\\w+=/g\
);\npush (@argl, ($cl=~/(\\S+=[^=]+\\S)\\s*$/g));\\
n\nforeach $a (@argl)\n  {\n    if ( ($cl=~/CXX=(.\
*)/)){$CXX=$1;}\n    if ( ($cl=~/-CC=(.*)/    )){$\
CC=$1;}\n    if ( ($cl=~/-FC=(.*)/    )){$FC=$1;}\\
n    if ( ($cl=~/-CFLAGS=(.*)/)){$CFLAGS=$1;}\n   \
 if ( ($cl=~/-CXXFLAGS=(.*)/)){$CXXFLAGS=$1;}\n  }\
\nour ($ROOT_INSTALL, $NO_QUESTION, $default_updat\
e_action,$BINARIES_ONLY,$force, $default_update_ac\
tion, $INSTALL_DIR, $PLUGINS_DIR, $DISTRIBUTIONS,$\
tclinkdb, $proxy, $clean);\nif ( ($cl=~/-root/)){$\
ROOT_INSTALL=1;}\nif ( ($cl=~/-no_question/)){$NO_\
QUESTION=1;}\nif ( ($cl=~/-update/)){$default_upda\
te_action=\"update\";}\nif ( ($cl=~/-binaries/)){$\
BINARIES_ONLY=1;}\nif ( ($cl=~/-force/)){$force=1;\
$default_update_action=\"update\"}\nif ( ($cl=~/-e\
xec=\\s*(\\S+)/)){$INSTALL_DIR=$1;}\nif ( ($cl=~/-\
plugins=\\s*(\\S+)/)){$PLUGINS_DIR=$1;}\nif ( ($cl\
=~/-dis=\\s*(\\S+)/)){$DISTRIBUTIONS=$1;}\n\nif ( \
($cl=~/-tclinkdb=\\s*(\\S+)/)){$tclinkdb=$1;}\nif \
( ($cl=~/-proxy=\\s*(\\S+)/)){$proxy=$1;}\nif ( ($\
cl=~/-clean/)){$clean=1;}\nif ( ($cl=~/-repo=\\s*(\
\\S+)/)){ $REPO_ROOT=$1; }\nif ($tclinkdb){&update\
_tclinkdb ($tclinkdb);}\n\n\nif( $REPO_ROOT ne \"\\
" ) {\n	if( $OSNAME eq \"\" ) { print \"You have s\
pecified the repository folder but the required \\\
\"OSNAME\\\" enviroment variable is missing. \\n\"\
; exit 1; } \n	if( $OSARCH eq \"\" ) { print \"You\
 have specified the repository folder but the requ\
ired \\\"OSARCH\\\" enviroment variable is missing\
. \\n\"; exit 1; } \n}\n\nour $TCDIR=$ENV{DIR_4_TC\
OFFEE};\nour $TCCACHE=$ENV{CACHE_4_TCOFFEE};\nour \
$TCTMP=$ENV{CACHE_4_TCOFFEE};\nour $TCM=$ENV{MCOFF\
EE_4_TCOFFEE};\nour $TCMETHODS=$ENV{METHODS_4_TCOF\
FEE};\nour $TCPLUGINS=$ENV{PLUGINS_4_TCOFFEE};\nou\
r $PLUGINS_DIR=\"\";\nour $INSTALL_DIR=\"\";\n\n&a\
dd_dir ($TCDIR=\"$HOME/.t_coffee\");\n&add_dir ($T\
CCACHE=\"$TCDIR/cache\");\n&add_dir ($TCTMP=\"$CDI\
R/tmp\");\n&add_dir ($TCM=\"$TCDIR/mcoffee\");\n&a\
dd_dir ($TCMETHODS=\"$TCDIR/methods\");\n&add_dir \
($TCPLUGINS=\"$TCDIR/plugins/$OS\");\n\n\nour $BAS\
E=\"$CD/bin\";\nour $BIN=\"$BASE/binaries/$OS\";\n\
our $DOWNLOAD_DIR=\"$BASE/download\";\nour $DOWNLO\
AD_FILE=\"$DOWNLOAD_DIR/files\";\nour $TMP=\"$BASE\
/tmp\";\n\n&add_dir($BASE);\n&add_dir($BIN);\n&add\
_dir($DOWNLOAD_DIR);\n&add_dir($DOWNLOAD_FILE);\ni\
f (!$DISTRIBUTIONS){$DISTRIBUTIONS=\"$DOWNLOAD_DIR\
/distributions\";}\n&add_dir ($DISTRIBUTIONS);\n&a\
dd_dir ($TMP);\n\n\nif    (!$PLUGINS_DIR && !$ROOT\
_INSTALL){$PLUGINS_DIR=$TCPLUGINS;}\nelsif (!$PLUG\
INS_DIR &&  $ROOT_INSTALL){$PLUGINS_DIR=\"/usr/loc\
al/bin/\";}\n\nif    (!$INSTALL_DIR && !$ROOT_INST\
ALL){$INSTALL_DIR=\"$HOME/bin/\";mkpath ($INSTALL_\
DIR);}\nelsif (!$INSTALL_DIR &&  $ROOT_INSTALL){$I\
NSTALL_DIR=\"/usr/local/bin/\";}\n\nif (-d \"mcoff\
ee\"){`cp mcoffee/* $TCM`;}\n\n\nour $ENV_FILE=\"$\
TCDIR/t_coffee_env\";\n&env_file2putenv ($ENV_FILE\
);\n&set_proxy($proxy);\nmy ($target, $p, $r);\n$t\
arget=$p;\n\nforeach $p (  ((keys (%PG)),(keys(%MO\
DE)),(@smode)) )\n  {\n    if ($ARGV[0] eq $p && $\
target eq \"\"){$target=$p;}\n  }\nif ($target eq \
\"\"){exit ($EXIT_FAILURE);}\n\n\nforeach $r (@req\
uired_applications)\n  {\n    my @app_list;\n    m\
y $i;\n    $i=0;\n    \n    @app_list=split (/_OR_\
/, $r);\n    foreach my $pg (@app_list)\n      {\n\
	$i+=&pg_is_installed ($pg);\n      }\n    if ($i=\
=0)\n      {\n      print \"One of the following p\
ackages must be installed to proceed: \";\n      f\
oreach my $pg (@app_list)\n	{\n	  print (\"$pg \")\
;\n	}\n      die;\n    }\n  }\n\n\n\n\n\n\n&sign_l\
icense_ni();\n\n\n$PG{C}{compiler}=get_C_compiler(\
$CC);\n$PG{Fortran}{compiler}=get_F_compiler($FC);\
\n$PG{CXX}{compiler}=$PG{CPP}{compiler}=$PG{GPP}{c\
ompiler}=get_CXX_compiler($CXX);\nif ($CXXFLAGS){$\
PG{CPP}{options}=$PG{GPP}{options}=$PG{CXX}{option\
s}=$CXXFLAGS;}\nif ($CFLAGS){$PG{C}{options}=$CFLA\
GS;}\nforeach my $c (keys(%PG))\n  {\n    my $argu\
ments;\n    if ($PG{$c}{compiler})\n      {\n	$arg\
uments=\"$PG{$c}{compiler_flag}=$PG{$c}{compiler} \
\";\n	if ($PG{$c}{options})\n	  {\n	    $arguments\
.=\"$PG{$c}{options_flag}=$PG{$c}{options} \";\n	 \
 }\n	$PG{$c}{arguments}=$arguments;\n      }\n  }\\
n\nif ($PG{$target}){$PG{$target}{install}=1;}\nel\
se\n  {\n    foreach my $pg (keys(%PG))\n      {\n\
	if ( $target eq \"all\" || ($PG{$pg}{mode}=~/$tar\
get/))\n	  {\n	    $PG{$pg} {install}=1;\n	  }\n  \
    }\n  }\n\nforeach my $pg (keys(%PG))\n  {\n   \
 if (!$PG{$pg}{update_action}){$PG{$pg}{update_act\
ion}=$default_update_action;}\n    elsif ($PG{$pg}\
{update_action} eq \"never\"){$PG{$pg}{install}=0;\
}\n    if ( $force && $PG{$pg}{install})\n      {\\
n	`rm $BIN/$pg $BIN/$pg.exe $SILENT`;\n      }\n  \
  if ($PG{$pg}{update_action} eq \"update\" && $PG\
{$pg}{install}){$PG{$pg}{update}=1;}\n  }\n\nif ((\
$target=~/clean/))\n  {\n    print \"------- clean\
ing executables -----\\n\";\n    `rm bin/* $SILENT\
`;\n    exit ($EXIT_SUCCESS);\n  }\n\nif ( !$PG{$t\
arget}){print \"------- Installing T-Coffee Modes\\
\n\";}\n\nforeach my $m (keys(%MODE))\n  {\n    if\
 ( $target eq \"all\" || $target eq $m)\n      {\n\
	print \"\\n------- The installer will now install\
 the $m components $MODE{$m}{description}\\n\";\n	\
foreach my $pg (keys(%PG))\n	  {\n	    if ( $PG{$p\
g}{mode} =~/$m/ && $PG{$pg}{install})\n	      {\n	\
	if ($PG{$pg}{touched}){print \"------- $PG{$pg}{d\
name}: already processed\\n\";}\n		else {$PG{$pg}{\
success}=&install_pg($pg);$PG{$pg}{touched}=1;}\n	\
      }\n	  }\n      }\n  }\n\nif ( $PG{$target}){\
print \"------- Installing Individual Package\\n\"\
;}\nforeach my $pg (keys (%PG))\n  {\n    \n    if\
 ( $PG{$pg}{install} && !$PG{$pg}{touched})\n     \
 {\n	print \"\\n------- Install $pg\\n\";\n	$PG{$p\
g}{success}=&install_pg($pg);$PG{$pg}{touched}=1;\\
n      }\n  }\nprint \"------- Finishing The insta\
llation\\n\";\nmy $final_report=&install ($INSTALL\
_DIR);\n\nprint \"\\n\";\nprint \"****************\
**************************************************\
***\\n\";\nprint \"********              INSTALLAT\
ION SUMMARY          *****************\\n\";\nprin\
t \"**********************************************\
***********************\\n\";\nprint \"------- SUM\
MARY package Installation:\\n\";\nprint \"------- \
  Executable Installed in: $PLUGINS_DIR\\n\";\n\nf\
oreach my $pg (keys(%PG))\n  {\n    if ( $PG{$pg}{\
install})\n      {\n	my $bin_status=($PG{$pg}{from\
_binary} && $PG{$pg}{success})?\"[from binary]\":\\
"\";\n	if     ( $PG{$pg}{new} && !$PG{$pg}{old})  \
                   {print \"*------        $PG{$pg\
}{dname}: installed $bin_status\\n\"; $PG{$pg}{sta\
tus}=1;}\n	elsif  ( $PG{$pg}{new} &&  $PG{$pg}{old\
})                     {print \"*------        $PG\
{$pg}{dname}: updated $bin_status\\n\"  ; $PG{$pg}\
{status}=1;} \n	elsif  (!$PG{$pg}{new} &&  $PG{$pg\
}{old} && !$PG{$pg}{update}){print \"*------      \
  $PG{$pg}{dname}: previous\\n\" ; $PG{$pg}{status\
}=1;}\n	elsif  (!$PG{$pg}{new} &&  $PG{$pg}{old} &\
&  $PG{$pg}{update}){print \"*------        $PG{$p\
g}{dname}: failed update (previous installation av\
ailable)\\n\";$PG{$pg}{status}=0;}\n	else         \
                                                 {\
print \"*------        $PG{$pg}{dname}: failed ins\
tallation\\n\";$PG{$pg}{status}=0;}\n      }\n  }\\
nmy $failure;\n\nif ( !$PG{$target}){print \"*----\
-- SUMMARY mode Installation:\\n\";}\nforeach my $\
m (keys(%MODE))\n  {\n  \n    if ( $target eq \"al\
l\" || $target eq $m)\n      {\n	my $succesful=1;\\
n	foreach my $pg (keys(%PG))\n	  {\n	    if (($PG{\
$pg}{mode}=~/$m/) && $PG{$pg}{install} && $PG{$pg}\
{status}==0)\n	      {\n		$succesful=0;\n		print \\
"*!!!!!!       $PG{$pg}{dname}: Missing\\n\";\n	  \
    }\n	  }\n	if ( $succesful)\n	  {\n	    $MODE{$\
m}{status}=1;\n	    print \"*------       MODE $MO\
DE{$m}{dname} SUCCESFULY installed\\n\";\n	  }\n	e\
lse\n	  {\n	    $failure++;\n	    $MODE{$m}{status\
}=0;\n	    print \"*!!!!!!       MODE $MODE{$m}{dn\
ame} UNSUCCESFULY installed\\n\";\n	  }\n      }\n\
  }\n\n    \n      \nif ($clean==1 && ($BASE=~/ins\
tall4tcoffee/) ){print \"*------ Clean Installatio\
n Directory: $BASE\\n\";`rm -rf $BASE`;}\nforeach \
my $pg (keys(%PG)){if ($PG{$pg}{install} && $PG{$p\
g}{status}==0){exit ($EXIT_FAILURE);}}\n\nif ($fai\
lure)\n  {\n    print \"**************************\
*******************************************\\n\";\\
n    print \"********     SOME PACKAGES FAILED TO \
INSTALL        *****************\\n\";\n    print \
\"************************************************\
*********************\\n\";\n    print \"\\nSome o\
f the reported failures may be due to connectivity\
 problems\";\n    print \"\\nRerun the installatio\
n and the installer will specifically try to insta\
ll the missing packages\";\n    print \"\\nIf this\
 Fails, go to the original website and install the\
 package manually\";\n  }\n\nprint \"*************\
**************************************************\
******\\n\";\nprint \"********              FINALI\
ZE YOUR INSTALLATION    *****************\\n\";\np\
rint \"*******************************************\
**************************\\n\";\nprint \"------- \
Your executables are in:\\n\"; \nprint \"-------  \
     $PLUGINS_DIR:\\n\";\nprint \"------- Add this\
 directory to your path with the following command\
:\\n\";\nprint \"-------       export PATH=$PLUGIN\
S_DIR:\\$PATH\\n\";\nprint \"------- Make this per\
manent by adding this line to the file:\\n\";\npri\
nt \"-------       $HOME/.bashrc\\n\";\nexit ($EXI\
T_SUCCESS);  \n  \nsub get_CXX_compiler\n  {\n    \
my $c=@_[0];\n    my (@clist)=(\"g++\");\n    \n  \
  return get_compil ($c, @clist);\n }\nsub get_C_c\
ompiler\n  {\n    my $c=@_[0];\n    my (@clist)=(\\
"gcc\", \"cc\", \"icc\");\n    \n    return get_co\
mpil ($c, @clist);\n }\n\nsub get_F_compiler\n  {\\
n    my ($c)=@_[0];\n    my @clist=(\"f77\", \"g77\
\",\"g95\", \"gfortran\", \"ifort\");\n    return \
get_compil ($c, @clist);\n  } \n       \nsub get_c\
ompil\n  {\n    my ($fav,@clist)=(@_);\n    \n    \
#return the first compiler found installed in the \
system. Check first the favorite\n    foreach my $\
c ($fav,@clist)\n      {\n	if  (&pg_is_installed (\
$c)){return $c;}\n      }\n    return \"\";\n  }\n\
sub exit_if_pg_not_installed\n  {\n    my (@arg)=(\
@_);\n    \n    foreach my $p (@arg)\n      {\n	if\
 ( !&pg_is_installed ($p))\n	  {\n	    print \"!!!\
!!!!! The $p utility must be installed for this in\
stallation to proceed [FATAL]\\n\";\n	    die;\n	 \
 }\n      }\n    return 1;\n  }\nsub set_proxy\n  \
{\n    my ($proxy)=(@_);\n    my (@list,$p);\n    \
\n    @list= (\"HTTP_proxy\", \"http_proxy\", \"HT\
TP_PROXY\", \"ALL_proxy\", \"all_proxy\",\"HTTP_pr\
oxy_4_TCOFFEE\",\"http_proxy_4_TCOFFEE\");\n    \n\
    if (!$proxy)\n      {\n	foreach my $p (@list)\\
n	  {\n	    if ( ($ENV_SET{$p}) || $ENV{$p}){$prox\
y=$ENV{$p};}\n	  }\n      }\n    foreach my $p(@li\
st){$ENV{$p}=$proxy;}\n  }\n	\nsub check_internet_\
connection\n  {\n    my $internet;\n    \n    if (\
 -e \"x\"){unlink (\"x\");}\n    if     (&pg_is_in\
stalled    (\"wget\")){`wget www.google.com -Ox >/\
dev/null 2>/dev/null`;}\n    elsif  (&pg_is_instal\
led    (\"curl\")){`curl www.google.com -ox >/dev/\
null 2>/dev/null`;}\n    else\n      {\n	printf st\
derr \"\\nERROR: No pg for remote file fetching [w\
get or curl][FATAL]\\n\";\n	exit ($EXIT_FAILURE);\\
n      }\n    \n    if ( !-e \"x\" || -s \"x\" < 1\
0){$internet=0;}\n    else {$internet=1;}\n    if \
(-e \"x\"){unlink \"x\";}\n    return $internet;\n\
  }\nsub url2file\n  {\n    my ($cmd, $file,$wget_\
arg, $curl_arg)=(@_);\n    my ($exit,$flag, $pg, $\
arg);\n    \n    if ($INTERNET || check_internet_c\
onnection ()){$INTERNET=1;}\n    else\n      {\n	p\
rint STDERR \"ERROR: No Internet Connection [FATAL\
:install.pl]\\n\";\n	exit ($EXIT_FAILURE);\n      \
}\n    \n    if     (&pg_is_installed    (\"wget\"\
)){$pg=\"wget\"; $flag=\"-O\";$arg=\"--tries=2 --c\
onnect-timeout=10 $wget_arg\";}\n    elsif  (&pg_i\
s_installed    (\"curl\")){$pg=\"curl\"; $flag=\"-\
o\";$arg=$curl_arg;}\n    else\n      {\n	printf s\
tderr \"\\nERROR: No pg for remote file fetching [\
wget or curl][FATAL]\\n\";\n	exit ($EXIT_FAILURE);\
\n      }\n    \n    \n    if (-e $file){unlink($f\
ile);}\n    $exit=system \"$pg $cmd $flag$file $ar\
g\";\n    return $exit;\n  }\n\nsub pg_is_installe\
d\n  {\n    my ($p, $dir)=(@_);\n    my ($r,$m, $r\
et);\n    my ($supported, $language, $compil);\n  \
  \n  \n    if ( $PG{$p})\n      {\n	$language=$PG\
{$p}{language2};\n	$compil=$PG{$language}{compiler\
};\n      }\n    \n    if ( $compil eq \"CPAN\")\n\
      {\n	if ( system (\"perl -M$p -e 1\")==$EXIT_\
SUCCESS){$ret=1;}\n	else {$ret=0;}\n      }\n    e\
lsif ($dir)\n      {\n	if (-e \"$dir/$p\" || -e \"\
$dir/$p\\.exe\"){$ret=1;}\n	else {$ret=0;}\n      \
}\n    elsif (-e \"$PLUGINS_DIR/$p\" || -e \"$PLUG\
INS_DIR/$p.exe\"){$ret=1;}\n    else\n      {\n	$r\
=`which $p 2>/dev/null`;\n	if ($r eq \"\"){$ret=0;\
}\n	else {$ret=1;}\n      }\n   \n    return $ret;\
\n  }\nsub install\n  {\n    my ($new_bin)=(@_);\n\
    my ($copied, $report);\n\n    \n    if (!$ROOT\
_INSTALL)\n      {\n	\n	if (-e \"$BIN/t_coffee\"){\
`$CP $BIN/t_coffee $INSTALL_DIR`};\n	`cp $BIN/* $P\
LUGINS_DIR`;\n	$copied=1;\n      }\n    else\n    \
  {\n	$copied=&root_run (\"You must be root to fin\
alize the installation\", \"$CP $BIN/* $INSTALL_DI\
R $SILENT\");\n      }\n    \n     \n  if ( !$copi\
ed)\n    {\n      $report=\"*!!!!!! Installation u\
nsuccesful. The executables have been left in $BAS\
E/bin\\n\";\n    }\n  elsif ( $copied && $ROOT)\n \
   {\n      $report=\"*------ Installation succesf\
ul. Your executables have been copied in $new_bin \
and are on your PATH\\n\";\n    }\n  elsif ( $copi\
ed && !$ROOT)\n    {\n      $report= \"*!!!!!! T-C\
offee and associated packages have been copied in:\
 $new_bin\\n\";\n      $report.=\"*!!!!!! This add\
ress is NOT in your PATH sytem variable\\n\";\n   \
   $report.=\"*!!!!!! You can do so by adding the \
following line in your ~/.bashrc file:\\n\";\n    \
  $report.=\"*!!!!!! export PATH=$new_bin:\\$PATH\\
\n\";\n    }\n  return $report;\n}\n\nsub sign_lic\
ense_ni\n  {\n    my $F=new FileHandle;\n    open \
($F, \"license.txt\");\n    while (<$F>)\n      {\\
n	print \"$_\";\n      }\n    close ($F);\n    \n \
   return;\n  }\n\nsub install_pg\n  {\n    my ($p\
g)=(@_);\n    my ($report, $previous, $language, $\
compiler, $return);\n    \n    if (!$PG{$pg}{insta\
ll}){return 1;}\n    \n    $previous=&pg_is_instal\
led ($pg);\n    \n    if ($PG{$pg}{update_action} \
eq \"no_update\" && $previous)\n      {\n	$PG{$pg}\
{old}=1;\n	$PG{$pg}{new}=0;\n	$return=1;\n      }\\
n    else\n      {\n	$PG{$pg}{old}=$previous;\n	\n\
	if ($PG{$pg} {language2} eq \"Perl\"){&install_pe\
rl_package ($pg);}\n	elsif ($BINARIES_ONLY && &ins\
tall_binary_package ($pg)){$PG{$pg}{from_binary}=1\
;}\n	elsif (&install_source_package ($pg)){;}\n	el\
se \n	  {\n	    \n	    if (!&supported_os($OS))\n	\
      {\n		print \"!!!!!!!! $pg compilation failed\
, binary unsupported for $OS\\n\"; \n	      }\n	  \
  elsif (!($PG{$pg}{from_binary}=&install_binary_p\
ackage ($pg)))\n	      {\n		print \"!!!!!!!! $pg c\
ompilation and  binary installation failed\\n\";\n\
	      }\n	  }\n	$PG{$pg}{new}=$return=&pg_is_inst\
alled ($pg,$BIN);\n      }\n\n    \n    return $re\
turn;\n  }\nsub install_perl_package\n  {\n    my \
($pg)=(@_);\n    my ($report, $language, $compiler\
);\n    \n    $language=$PG{$pg} {language2};\n   \
 $compiler=$PG{$language}{compiler};\n    \n    if\
 (!&pg_is_installed ($pg))\n      {\n	if ( $OS eq \
\"windows\"){`perl -M$compiler -e 'install $pg'`;}\
\n	elsif ( $ROOT eq \"sudo\"){system (\"sudo perl \
-M$compiler -e 'install $pg'\");}\n	else {system (\
\"su root -c perl -M$compiler -e 'install $pg'\");\
}\n      }\n    return &pg_is_installed ($pg);\n  \
}\n\n\n\nsub install_source_package\n  {\n    my (\
$pg)=(@_);\n    my ($report, $download, $arguments\
, $language, $address, $name, $ext, $main_dir, $di\
strib);\n    my $wget_tmp=\"$TMP/wget.tmp\";\n    \
my (@fl);\n    if ( -e \"$BIN/$pg\" || -e \"$BIN/$\
pg.exe\"){return 1;}\n    \n    #\n    # check if \
the module exists in the repository cache \n    #\\
n	if( repo_load($pg) ) {\n		return 1;\n	}\n    \n \
   if ($pg eq \"t_coffee\")  {return   &install_t_\
coffee ($pg);}\n    elsif ($pg eq \"TMalign\"){ret\
urn   &install_TMalign ($pg);}\n    \n    chdir $D\
ISTRIBUTIONS;\n    \n    $download=$PG{$pg}{source\
};\n    \n    if (($download =~/tgz/))\n      {\n	\
($address,$name,$ext)=($download=~/(.+\\/)([^\\/]+\
)(\\.tgz).*/);\n      }\n    elsif (($download=~/t\
ar\\.gz/))\n      {\n	($address,$name,$ext)=($down\
load=~/(.+\\/)([^\\/]+)(\\.tar\\.gz).*/);\n      }\
\n    elsif (($download=~/tar/))\n      {\n	($addr\
ess,$name,$ext)=($download=~/(.+\\/)([^\\/]+)(\\.t\
ar).*/);\n      }\n    else\n      {\n	($address,$\
name)=($download=~/(.+\\/)([^\\/]+)/);\n	$ext=\"\"\
;\n      }\n    $distrib=\"$name$ext\";\n    \n   \
 if ( !-d $pg){mkdir $pg;}\n    chdir $pg;\n   \n \
   #get the distribution if available\n    if ( -e\
 \"$DOWNLOAD_DIR/$distrib\")\n      {\n	`$CP $DOWN\
LOAD_DIR/$distrib .`;\n      }\n    #UNTAR and Pre\
pare everything\n    if (!-e \"$name.tar\" && !-e \
\"$name\")\n      {\n	&check_rm ($wget_tmp);\n	pri\
nt \"\\n------- Downloading/Installing $pg\\n\";\n\
	\n	if (!-e $distrib && &url2file (\"$download\", \
\"$wget_tmp\")==$EXIT_SUCCESS)\n	  {\n	    \n	    \
`mv $wget_tmp $distrib`;\n	    `$CP $distrib $DOWN\
LOAD_DIR/`;\n	  }\n\n	if (!-e $distrib)\n	  {\n	  \
  print \"!!!!!!! Download of $pg distribution fai\
led\\n\";\n	    print \"!!!!!!! Check Address: $PG\
{$pg}{source}\\n\";\n	    return 0;\n	  }\n	print \
\"\\n------- unzipping/untaring $name\\n\";\n	if (\
($ext =~/z/))\n	  { \n	    &flush_command (\"gunzi\
p $name$ext\");\n	    \n	  }\n	if (($ext =~/tar/) \
|| ($ext =~/tgz/))\n	  {\n	    &flush_command(\"ta\
r -xvf $name.tar\");\n	  }\n      }\n    #Guess an\
d enter the distribution directory\n    @fl=ls($p)\
;\n    foreach my $f (@fl)\n      {\n	if (-d $f)\n\
	  {\n	    $main_dir=$f;\n	  }\n      }\n    if (-\
d $main_dir)\n	  \n      {\n	chdir $main_dir;}\n  \
  else\n      {\n	print \"Error: $main_dir does no\
t exist\";\n      }\n    print \"\\n------- Compil\
ing/Installing $pg\\n\";\n    `make clean $SILENT`\
;\n    \n    \n    #\n    # SAP module\n    #\n   \
 if ($pg eq \"sap\")\n      {\n	if (-e \"./configu\
re\")\n	  {\n	    #new sap distribution\n	    if (\
$OS eq \"macosx\")\n	      {\n		&replace_line_in_f\
ile (\"./src/galloc.h\", \"malloc.h\",  \"\");\n		\
&replace_line_in_file (\"./src/pdbprot.h\", \"mall\
oc.h\", \"\");\n		&replace_line_in_file (\"./src/p\
dbprot.c\", \"malloc.h\", \"\");\n	      }\n	    \\
n	    &flush_command (\"./configure\");\n	    &flu\
sh_command (\"make clean\");\n	    &flush_command \
(\"make\");\n	    &check_cp (\"./src/$pg\", \"$BIN\
\");\n	    repo_store(\"./src/$pg\");\n	  }\n	else\
\n	  {\n	    #old style distribution\n	    `rm *.o\
 sap  sap.exe ./util/aa/*.o  ./util/wt/.o $SILENT`\
;\n	    &flush_command (\"make $arguments sap\");\\
n	    &check_cp ($pg, \"$BIN\");\n	    repo_store(\
$pg);\n	  }\n      }\n    \n    #\n    # CLUSTALW2\
 module\n    #\n    elsif ($pg eq \"clustalw2\")\n\
      {\n	&flush_command(\"./configure\");\n	&flus\
h_command(\"make $arguments\");\n	&check_cp (\"./s\
rc/$pg\", \"$BIN\");\n	repo_store(\"./src/$pg\");\\
n      }\n    \n    #\n    # FSA module\n    # \n \
   elsif ($pg eq \"fsa\")\n      {\n	&flush_comman\
d(\"./configure --prefix=$BIN\");\n	&flush_command\
(\"make $arguments\");\n	&flush_command (\"make in\
stall\");\n\n	repo_store(\"fsa\", \"$BIN/bin\");\n\
	`mv $BIN/bin/* $BIN`;\n	`rmdir $BIN/bin`;\n      \
}\n    \n    #\n    # CLUSTALW module\n    #\n    \
elsif ($pg eq \"clustalw\")\n      {\n	&flush_comm\
and(\"make $arguments clustalw\");\n	`$CP $pg $BIN\
 $SILENT`;\n	repo_store($pg);\n      }\n    \n    \
#\n    # MAFFT module\n    #\n    elsif ($pg eq \"\
mafft\")\n      {\n	my $base=cwd();\n	my $c;\n	\n	\
#compile core\n	mkpath (\"./mafft/bin\");\n	mkpath\
 (\"./mafft/lib\");\n	chdir \"$base/core\";\n	`mak\
e clean $SILENT`;\n	&flush_command (\"make $argume\
nts\");\n	&flush_command (\"make install LIBDIR=..\
/mafft/lib BINDIR=../mafft/bin\");\n	\n	#compile e\
xtension\n	chdir \"$base/extensions\";\n	`make cle\
an $SILENT`;\n	&flush_command (\"make $arguments\"\
);\n	&flush_command (\"make install LIBDIR=../maff\
t/lib BINDIR=../mafft/bin\");\n	\n	#put everything\
 in mafft and copy the compiled stuff in bin\n	chd\
ir \"$base\";\n	if ($ROOT_INSTALL)\n	  {\n	    &ro\
ot_run (\"You Must be Root to Install MAFFT\\n\", \
\"mkdir /usr/local/mafft/;$CP mafft/lib/* /usr/loc\
al/mafft;$CP mafft/lib/mafft* /usr/local/bin ;$CP \
mafft/bin/mafft /usr/local/bin/; \");\n	  }\n	else\
\n	  {\n	    `$CP mafft/lib/*  $BIN`;\n	    `$CP m\
afft/bin/mafft  $BIN`;\n	  }\n	`tar -cvf mafft.tar\
 mafft`;\n	`gzip mafft.tar`;\n	`mv mafft.tar.gz $B\
IN`;\n	\n	repo_store(\"mafft/bin/mafft\", \"mafft/\
lib/\", \"$BIN/mafft.tar.gz\");\n      }\n      \n\
    #\n    # DIALIGN-TX module\n    #\n    elsif (\
 $pg eq \"dialign-tx\" )\n      {\n	my $f;\n	my $b\
ase=cwd();\n\n	chdir \"./source\";\n	if ($OS eq \"\
macosx\"){&flush_command (\"cp makefile.MAC_OS mak\
efile\");}\n\n	&flush_command (\" make CPPFLAGS='-\
O3 -funroll-loops' all\");\n	\n	chdir \"..\";\n	&c\
heck_cp (\"./source/$pg\", \"$BIN\");\n	repo_store\
(\"./source/$pg\");\n      }\n      \n    #\n    #\
 DIALIGN-T module \n    # (is the same as dialign-\
tx, but it is mantained for backward name compatib\
ility with tcoffee)\n    #\n    elsif ( $pg eq \"d\
ialign-t\" )\n      {\n	my $f;\n	my $base=cwd();\n\
\n	chdir \"./source\";\n	if ($OS eq \"macosx\"){&f\
lush_command (\"cp makefile.MAC_OS makefile\");}\n\
\n	&flush_command (\" make CPPFLAGS='-O3 -funroll-\
loops' all\");\n	\n	chdir \"..\";\n	&check_cp (\".\
/source/dialign-tx\", \"$BIN/dialign-t\");\n	repo_\
store(\"$BIN/dialign-t\");	\n      }      \n      \
\n    #\n    # POA module\n    #\n    elsif ($pg e\
q \"poa\")\n      {\n	&flush_command (\"make $argu\
ments poa\");\n	&check_cp (\"$pg\", \"$BIN\");\n	r\
epo_store(\"$pg\");\n      }\n     \n     \n    #\\
n    # PROBCONS module\n    #\n    elsif ( $pg eq \
\"probcons\")\n      {\n	&add_C_libraries(\"./Prob\
abilisticModel.h\", \"list\", \"cstring\");\n	\n	`\
rm *.exe $SILENT`;\n	&flush_command (\"make $argum\
ents probcons\");\n	&check_cp(\"$pg\", \"$BIN/$pg\\
");\n	repo_store(\"$pg\");\n      }\n      \n    #\
\n    # PROBCONS RNA module\n    #\n    elsif ( $p\
g eq \"probconsRNA\")\n      {\n	&add_C_libraries(\
\"./ProbabilisticModel.h\", \"list\", \"cstring\")\
;\n	&add_C_libraries(\"./Main.cc\", \"iomanip\", \\
"cstring\",\"climits\");\n	`rm *.exe $SILENT`;\n	&\
flush_command (\"make $arguments probcons\");\n	&c\
heck_cp(\"probcons\", \"$BIN/$pg\");\n	repo_store(\
\"$BIN/$pg\");\n      }\n\n	#\n	# MUSCLE module\n	\
#\n    elsif (  $pg eq \"muscle\")\n      {	\n	`rm\
 *.o muscle muscle.exe $SILENT`;\n	if ($OS eq \"ma\
cosx\" || $OS eq \"linux\")\n	  {\n	    &replace_l\
ine_in_file (\"./Makefile\", \"LDLIBS = -lm -stati\
c\",  \"LDLIBS = -lm\");\n	  }\n	elsif ($OS eq \"w\
indows\")\n	  {\n	    &replace_line_in_file (\"./i\
ntmath.cpp\",  \"double log2e\",      \"double ced\
ric_log\");\n	    &replace_line_in_file (\"./intma\
th.cpp\",  \"double log2\",       \"double log_not\
use\");\n	    &replace_line_in_file (\"./intmath.c\
pp\",  \"double cedric_log\", \"double log2e\");\n\
	  }\n	&flush_command (\"make $arguments all\");\n\
	&check_cp(\"$pg\", \"$BIN\");\n	repo_store(\"$pg\\
");	\n      }\n      \n     #\n     # MUS4 module\\
n     #\n     elsif (  $pg eq \"mus4\")\n      {\n\
	`rm *.o muscle muscle.exe $SILENT`;\n	&flush_comm\
and (\"./mk\");\n	&check_cp(\"$pg\", \"$BIN\");\n	\
repo_store(\"$pg\");	\n      }\n      \n    #\n   \
 # PCMA module\n    #\n    elsif ( $pg eq \"pcma\"\
)\n      {\n	if ($OS eq \"macosx\")\n	  {\n	    &r\
eplace_line_in_file (\"./alcomp2.c\", \"malloc.h\"\
,  \"\");\n	  }\n	&flush_command (\"make $argument\
s pcma\");\n	&check_cp(\"$pg\", \"$BIN\");\n	repo_\
store(\"$pg\");	\n      }\n      \n    #\n    # KA\
LIGN module\n    #\n    elsif ($pg eq \"kalign\")\\
n      {\n	&flush_command (\"./configure\");\n	&fl\
ush_command(\"make $arguments\");\n	&check_cp (\"$\
pg\",$BIN);\n	repo_store(\"$pg\");	\n      }\n    \
  \n    #\n    # AMAP module\n    #\n    elsif ( $\
pg eq \"amap\")\n      {\n	&add_C_libraries(\"./Am\
ap.cc\", \"iomanip\", \"cstring\",\"climits\");	\n\
	`make clean $SILENT`;\n	&flush_command (\"make $a\
rguments all\");\n	&check_cp (\"$pg\", $BIN);\n	re\
po_store(\"$pg\");	\n      }\n      \n    #\n    #\
 PRODA module\n    #\n    elsif ( $pg eq \"proda\"\
)\n      {\n	&add_C_libraries(\"AlignedFragment.h\\
", \"vector\", \"iostream\", \"cstring\",\"cstdlib\
\");\n	&add_C_libraries(\"Main.cc\", \"vector\", \\
"climits\");	\n	&add_C_libraries(\"Sequence.cc\", \
\"stdlib.h\", \"cstdio\");	\n	&flush_command (\"ma\
ke $arguments all\");\n	&check_cp (\"$pg\", $BIN);\
\n	repo_store(\"$pg\");	\n      }\n      \n    #\n\
    # PRANK module\n    #\n    elsif ( $pg eq \"pr\
ank\")\n      {\n	&flush_command (\"make $argument\
s all\");\n	&check_cp (\"$pg\", $BIN);\n	repo_stor\
e(\"$pg\");	\n      }\n      \n    #\n    # !!!! M\
USTANG module\n    #\n     elsif ( $pg eq \"mustan\
g\")\n      {\n	&flush_command (\"rm ./bin/*\");\n\
	&flush_command (\"make $arguments all\");\n\n	if \
( $OS=~/windows/){&flush_command(\"cp ./bin/* $BIN\
/mustang.exe\");}\n	else {&flush_command(\"cp ./bi\
n/* $BIN/mustang\");}\n	\n	repo_store(\"$BIN/musta\
ng\");\n      }\n\n	#\n	# RNAplfold module\n	#\n  \
  elsif ( $pg eq \"RNAplfold\")\n      {\n	&flush_\
command(\"./configure\");\n	&flush_command (\"make\
 $arguments all\");\n	&check_cp(\"./Progs/RNAplfol\
d\", \"$BIN\");\n	&check_cp(\"./Progs/RNAalifold\"\
, \"$BIN\");\n	&check_cp(\"./Progs/RNAfold\", \"$B\
IN\");\n	\n	repo_store(\"./Progs/RNAplfold\", \"./\
Progs/RNAalifold\", \"./Progs/RNAfold\");\n      }\
\n      \n    #\n    # !!! RETREE module\n    #\n \
   elsif ( $pg eq \"retree\")\n      {\n	chdir \"s\
rc\";\n	&flush_command (\"make $arguments all\");\\
n	&flush_command (\"make put\");\n	system \"cp ../\
exe/* $BIN\";\n	\n	repo_store(\"retree\", \"../exe\
\");\n      }\n	\n    chdir $CDIR;\n    return &pg\
_is_installed ($pg, $BIN);\n  }\n\nsub install_t_c\
offee\n  {\n    my ($pg)=(@_);\n    my ($report,$c\
flags, $arguments, $language, $compiler) ;\n    #1\
-Install T-Coffee\n    chdir \"t_coffee_source\";\\
n    &flush_command (\"make clean\");\n    print \\
"\\n------- Compiling T-Coffee\\n\";\n    $languag\
e=$PG{$pg} {language2};\n    $arguments=$PG{$langu\
age}{arguments};\n    if (!($arguments =~/CFLAGS/)\
){$arguments .= \" CFLAGS=-O2 \";}\n\n    if ( $CC\
 ne \"\"){&flush_command (\"make -i $arguments t_c\
offee\");}\n    &check_cp ($pg, $BIN);\n    \n    \
chdir $CDIR;\n    return &pg_is_installed ($pg, $B\
IN);\n  }\nsub install_TMalign\n  {\n    my ($pg)=\
(@_);\n    my $report;\n    chdir \"t_coffee_sourc\
e\";\n    print \"\\n------- Compiling TMalign\\n\\
";\n    `rm TMalign TMalign.exe $SILENT`;\n    if \
( $FC ne \"\"){&flush_command (\"make -i $PG{Fortr\
an}{arguments} TMalign\");}\n    &check_cp ($pg, $\
BIN);\n    repo_store($pg);\n\n    if ( !-e \"$BIN\
/$pg\" && pg_has_binary_distrib ($pg))\n      {\n	\
print \"!!!!!!! Compilation of $pg impossible. Wil\
l try to install from binary\\n\";\n	return &insta\
ll_binary_package ($pg);\n      }\n    chdir $CDIR\
;\n    return &pg_is_installed ($pg, $BIN);\n  }\n\
\nsub pg_has_binary_distrib\n  {\n    my ($pg)=(@_\
);\n    if ($PG{$pg}{windows}){return 1;}\n    els\
if ($PG{$pg}{osx}){return 1;}\n    elsif ($PG{$pg}\
{linux}){return 1;}\n    return 0;\n  }\nsub insta\
ll_binary_package\n  {\n    my ($pg)=(@_);\n    my\
 ($base,$report,$name, $download, $arguments, $lan\
guage, $dir);\n    my $isdir;\n    &input_os();\n \
   \n    if (!&supported_os($OS)){return 0;}\n    \
if ( $PG{$pg}{binary}){$name=$PG{$pg}{binary};}\n \
   else \n      {\n	$name=$pg;\n	if ( $OS eq \"win\
dows\"){$name.=\".exe\";}\n      }\n    \n    $dow\
nload=\"$WEB_BASE/Packages/Binaries/$OS/$name\";\n\
    \n    $base=cwd();\n    chdir $TMP;\n    \n   \
 if (!-e $name)\n      {\n	`rm x $SILENT`;\n	if ( \
url2file(\"$download\",\"x\")==$EXIT_SUCCESS)\n	  \
{\n	    `mv x $name`;\n	  }\n      }\n    \n    if\
 (!-e $name)\n      {\n	print \"!!!!!!! $PG{$pg}{d\
name}: Download of $pg binary failed\\n\";\n	print\
 \"!!!!!!! $PG{$pg}{dname}: Check Address: $downlo\
ad\\n\";\n	return 0;\n      }\n    print \"\\n----\
--- Installing $pg\\n\";\n    \n    if ($name =~/t\
ar\\.gz/)\n      {\n	`gunzip  $name`;\n	`tar -xvf \
$pg.tar`;\n	chdir $pg;\n	if ( $pg eq \"mafft\")\n	\
  {\n	    if ($ROOT_INSTALL)\n	      {\n		&root_ru\
n (\"You Must be Roor to Install MAFFT\\n\", \"$CP\
 mafft/bin/* /usr/local/mafft;mkdir /usr/local/maf\
ft/; $CP mafft/lib/* /usr/local/bin/\");\n	      }\
\n	    else\n	      {\n		`$CP $TMP/$pg/bin/* $BIN \
$SILENT`;\n		`$CP $TMP/$pg/lib/* $BIN $SILENT`;\n	\
      }\n	  }\n	else\n	  {\n	    if (-e \"$TMP/$pg\
/data\"){`$CP $TMP/$pg/data/* $TCM $SILENT`;}\n	  \
  if (!($pg=~/\\*/)){`rm -rf $pg`;}\n	  }\n      }\
\n    else\n      {\n	&check_cp (\"$pg\", \"$BIN\"\
);\n	`chmod u+x $BIN/$pg`; \n	unlink ($pg);\n     \
 }\n    chdir $base;\n    $PG{$pg}{from_binary}=1;\
\n    return &pg_is_installed ($pg, $BIN);\n  }\n\\
nsub add_dir \n  {\n    my $dir=@_[0];\n    \n    \
if (!-e $dir && !-d $dir)\n      {\n	my @l;\n	umas\
k (0000);\n	@l=mkpath ($dir,{mode => 0777});\n	\n \
     }\n    else\n      {\n	return 0;\n      }\n  \
}\nsub check_rm \n  {\n    my ($file)=(@_);\n    \\
n    if ( -e $file)\n      {\n	return unlink($file\
);\n      }\n    return 0;\n  }\nsub check_cp\n  {\
\n    my ($from, $to)=(@_);\n    if ( !-e $from &&\
 -e \"$from\\.exe\"){$from=\"$from\\.exe\";}\n    \
if ( !-e $from){return 0;}\n        \n    `$CP $fr\
om $to`;\n    return 1;\n  }\n\nsub repo_store \n{\
\n   # check that all required data are available\\
n   if( $REPO_ROOT eq \"\" ) { return; }\n\n\n    \
# extract the package name from the specified path\
\n    my $pg =`basename $_[0]`;\n    chomp($pg);\n\
	\n    my $VER = $PG{$pg}{version};\n    my $CACHE\
 = \"$REPO_ROOT/$pg/$VER/$OSNAME-$OSARCH\"; \n    \
\n    print \"~~~~~~~ Storing package: \\\"$pg\\\"\
 to path: $CACHE\\n\";\n    \n    # clean the cach\
e path if exists and create it again\n    `rm -rf \
$CACHE`;\n    `mkdir -p $CACHE`;\n    \n 	for my $\
path (@_) {\n\n	    # check if it is a single file\
 \n	 	if( -f $path ) {\n	    	`cp $path $CACHE`;\n\
		}\n		# .. or a directory, in this case copy all \
the content \n		elsif( -d $path ) {\n			opendir(IM\
D, $path);\n			my @thefiles= readdir(IMD);\n			clo\
sedir(IMD);\n			\n			for my $_file (@thefiles) {\n\
				if( $_file ne \".\" && $_file ne \"..\") {\n	 \
   			`cp $path/$_file $CACHE`;\n				}\n			}\n		} \
\n	}	   \n    \n	\n}   \n\nsub repo_load \n{\n    \
my ($pg)=(@_);\n\n    # check that all required da\
ta are available\n    if( $REPO_ROOT eq \"\" ) { r\
eturn 0; }\n\n    my $VER = $PG{$pg}{version};\n  \
  my $CACHE = \"$REPO_ROOT/$pg/$VER/$OSNAME-$OSARC\
H\"; \n    if( !-e \"$CACHE/$pg\" ) {\n   	 	print\
 \"~~~~~~~ Module \\\"$pg\\\" NOT found on reposit\
ory cache.\\n\";\n    	return 0;\n    }\n    \n   \
 print \"~~~~~~~ Module \\\"$pg\\\" found on repos\
itory cache. Using copy on path: $CACHE\\n\";\n   \
 `cp $CACHE/* $BIN`;\n    return 1;\n}\n\nsub chec\
k_file_list_exists \n  {\n    my ($base, @flist)=(\
@_);\n    my $f;\n\n    foreach $f (@flist)\n     \
 {\n	if ( !-e \"$base/$f\"){return 0;}\n      }\n \
   return 1;\n  }\nsub ls\n  {\n    my $f=@_[0];\n\
    my @fl;\n    chomp(@fl=`ls -1 $f`);\n    retur\
n @fl;\n  }\nsub flush_command\n  {\n    my $comma\
nd=@_[0];\n    my $F=new FileHandle;\n    open ($F\
, \"$command|\");\n    while (<$F>){print \"    --\
- $_\";}\n    close ($F);\n  }    \n\nsub input_in\
stallation_directory\n  {\n    my $dir=@_[0];\n   \
 my $new;\n    \n    print \"------- The current i\
nstallation directory is: [$dir]\\n\";\n    print \
\"??????? Return to keep the default or new value:\
\";\n   \n    if ($NO_QUESTION==0)\n      {\n	chom\
p ($new=<stdin>);\n	while ( $new ne \"\" && !input\
_yes (\"You have entered $new. Is this correct? ([\
y]/n):\"))\n	  {\n	    print \"???????New installa\
tion directory:\";\n	    chomp ($new=<stdin>);\n	 \
 }\n	$dir=($new eq \"\")?$dir:$new;\n	$dir=~s/\\/$\
//;\n      }\n    \n    if ( -d $dir){return $dir;\
}\n    elsif (&root_run (\"You must be root to cre\
ate $dir\",\"mkdir $dir\")==$EXIT_SUCCESS){return \
$dir;}\n    else\n      {\n	print \"!!!!!!! $dir c\
ould not be created\\n\";\n	if ( $NO_QUESTION)\n	 \
 {\n	    return \"\";\n	  }\n	elsif ( &input_yes (\
\"??????? Do you want to provide a new directory([\
y]/n)?:\"))\n	  {\n	    return input_installation_\
directory ($dir);\n	  }\n	else\n	  {\n	    return \
\"\";\n	  }\n      }\n    \n  }\nsub input_yes\n  \
{\n    my $question =@_[0];\n    my $answer;\n\n  \
  if ($NO_QUESTION==1){return 1;}\n    \n    if ($\
question eq \"\"){$question=\"??????? Do you wish \
to proceed ([y]/n)?:\";}\n    print $question;\n  \
  chomp($answer=lc(<STDIN>));\n    if (($answer=~/\
^y/) || $answer eq \"\"){return 1;}\n    elsif ( (\
$answer=~/^n/)){return 0;}\n    else\n      {\n	re\
turn input_yes($question);\n      }\n  }\nsub root\
_run\n  {\n    my ($txt, $cmd)=(@_);\n    \n    if\
 ( system ($cmd)==$EXIT_SUCCESS){return $EXIT_SUCC\
ESS;}\n    else \n      {\n	print \"------- $txt\\\
n\";\n	if ( $ROOT eq \"sudo\"){return system (\"su\
do $cmd\");}\n	else {return system (\"su root -c \\
\\"$cmd\\\"\");}\n      }\n  }\nsub get_root\n  {\\
n    if (&pg_is_installed (\"sudo\")){return \"sud\
o\";}\n    else {return \"su\";}\n  }\n\nsub get_o\
s\n  {\n    my $raw_os=`uname`;\n    my $os;\n\n  \
  $raw_os=lc ($raw_os);\n    \n    if ($raw_os =~/\
cygwin/){$os=\"windows\";}\n    elsif ($raw_os =~/\
linux/){$os=\"linux\";}\n    elsif ($raw_os =~/osx\
/){$os=\"macosx\";}\n    elsif ($raw_os =~/darwin/\
){$os=\"macosx\";}\n    else\n      {\n	$os=$raw_o\
s;\n      }\n    return $os;\n  }\nsub input_os\n \
 {\n    my $answer;\n    if ($OS) {return $OS;}\n \
   \n    print \"??????? which os do you use: [w]i\
ndows, [l]inux, [m]acosx:?\";\n    $answer=lc(<STD\
IN>);\n\n    if (($answer=~/^m/)){$OS=\"macosx\";}\
\n    elsif ( ($answer=~/^w/)){$OS=\"windows\";}\n\
    elsif ( ($answer=~/^linux/)){$OS=\"linux\";}\n\
    \n    else\n      {\n	return &input_os();\n   \
   }\n    return $OS;\n  }\n\nsub supported_os\n  \
{\n    my ($os)=(@_[0]);\n    return $SUPPORTED_OS\
{$os};\n  }\n    \n    \n\n\nsub update_tclinkdb \\
n  {\n    my $file =@_[0];\n    my $name;\n    my \
$F=new FileHandle;\n    my ($download, $address, $\
name, $l, $db);\n    \n    if ( $file eq \"update\\
"){$file=$TCLINKDB_ADDRESS;}\n    \n    if ( $file\
 =~/http:\\/\\// || $file =~/ftp:\\/\\//)\n      {\
\n	($address, $name)=($download=~/(.*)\\/([^\\/]+)\
$/);\n	`rm x $SILENT`;\n	if (&url2file ($file,\"x\\
")==$EXIT_SUCCESS)\n	  {\n	    print \"------- Sus\
scessful upload of $name\";\n	    `mv x $name`;\n	\
    $file=$name;\n	  }\n      }\n    open ($F, \"$\
file\");\n    while (<$F>)\n      {\n	my $l=$_;\n	\
if (($l =~/^\\/\\//) || ($db=~/^#/)){;}\n	elsif ( \
!($l =~/\\w/)){;}\n	else\n	  {\n	    my @v=split (\
/\\s+/, $l);\n	    if ( $l=~/^MODE/)\n	      {\n		\
$MODE{$v[1]}{$v[2]}=$v[3];\n	      }\n	    elsif (\
$l=~/^PG/)\n	      {\n		$PG{$v[1]}{$v[2]}=$v[3];\n\
	      }\n	  }\n      }\n    close ($F);\n    &pos\
t_process_PG();\n    return;\n  }\n\n\n\nsub initi\
alize_PG\n  {\n\n$PG{\"t_coffee\"}{\"4_TCOFFEE\"}=\
\"TCOFFEE\";\n$PG{\"t_coffee\"}{\"type\"}=\"sequen\
ce_multiple_aligner\";\n$PG{\"t_coffee\"}{\"ADDRES\
S\"}=\"http://www.tcoffee.org\";\n$PG{\"t_coffee\"\
}{\"language\"}=\"C\";\n$PG{\"t_coffee\"}{\"langua\
ge2\"}=\"C\";\n$PG{\"t_coffee\"}{\"source\"}=\"htt\
p://www.tcoffee.org/Packages/T-COFFEE_distribution\
.tar.gz\";\n$PG{\"t_coffee\"}{\"update_action\"}=\\
"always\";\n$PG{\"t_coffee\"}{\"mode\"}=\"tcoffee,\
mcoffee,rcoffee,expresso,3dcoffee\";\n$PG{\"clusta\
lw2\"}{\"4_TCOFFEE\"}=\"CLUSTALW2\";\n$PG{\"clusta\
lw2\"}{\"type\"}=\"sequence_multiple_aligner\";\n$\
PG{\"clustalw2\"}{\"ADDRESS\"}=\"http://www.clusta\
l.org\";\n$PG{\"clustalw2\"}{\"language\"}=\"C++\"\
;\n$PG{\"clustalw2\"}{\"language2\"}=\"CXX\";\n$PG\
{\"clustalw2\"}{\"source\"}=\"http://www.clustal.o\
rg/download/2.0.10/clustalw-2.0.10-src.tar.gz\";\n\
$PG{\"clustalw2\"}{\"mode\"}=\"mcoffee,rcoffee\";\\
n$PG{\"clustalw2\"}{\"version\"}=\"2.0.10\";\n$PG{\
\"clustalw\"}{\"4_TCOFFEE\"}=\"CLUSTALW\";\n$PG{\"\
clustalw\"}{\"type\"}=\"sequence_multiple_aligner\\
";\n$PG{\"clustalw\"}{\"ADDRESS\"}=\"http://www.cl\
ustal.org\";\n$PG{\"clustalw\"}{\"language\"}=\"C\\
";\n$PG{\"clustalw\"}{\"language2\"}=\"C\";\n$PG{\\
"clustalw\"}{\"source\"}=\"http://www.clustal.org/\
download/1.X/ftp-igbmc.u-strasbg.fr/pub/ClustalW/c\
lustalw1.82.UNIX.tar.gz\";\n$PG{\"clustalw\"}{\"mo\
de\"}=\"mcoffee,rcoffee\";\n$PG{\"clustalw\"}{\"ve\
rsion\"}=\"1.82\";\n$PG{\"dialign-t\"}{\"4_TCOFFEE\
\"}=\"DIALIGNT\";\n$PG{\"dialign-t\"}{\"type\"}=\"\
sequence_multiple_aligner\";\n$PG{\"dialign-t\"}{\\
"ADDRESS\"}=\"http://dialign-tx.gobics.de/\";\n$PG\
{\"dialign-t\"}{\"DIR\"}=\"/usr/share/dialign-tx/\\
";\n$PG{\"dialign-t\"}{\"language\"}=\"C\";\n$PG{\\
"dialign-t\"}{\"language2\"}=\"C\";\n$PG{\"dialign\
-t\"}{\"source\"}=\"http://dialign-tx.gobics.de/DI\
ALIGN-TX_1.0.2.tar.gz\";\n$PG{\"dialign-t\"}{\"mod\
e\"}=\"mcoffee\";\n$PG{\"dialign-t\"}{\"binary\"}=\
\"dialign-t\";\n$PG{\"dialign-t\"}{\"version\"}=\"\
1.0.2\";\n$PG{\"dialign-tx\"}{\"4_TCOFFEE\"}=\"DIA\
LIGNTX\";\n$PG{\"dialign-tx\"}{\"type\"}=\"sequenc\
e_multiple_aligner\";\n$PG{\"dialign-tx\"}{\"ADDRE\
SS\"}=\"http://dialign-tx.gobics.de/\";\n$PG{\"dia\
lign-tx\"}{\"DIR\"}=\"/usr/share/dialign-tx/\";\n$\
PG{\"dialign-tx\"}{\"language\"}=\"C\";\n$PG{\"dia\
lign-tx\"}{\"language2\"}=\"C\";\n$PG{\"dialign-tx\
\"}{\"source\"}=\"http://dialign-tx.gobics.de/DIAL\
IGN-TX_1.0.2.tar.gz\";\n$PG{\"dialign-tx\"}{\"mode\
\"}=\"mcoffee\";\n$PG{\"dialign-tx\"}{\"binary\"}=\
\"dialign-tx\";\n$PG{\"dialign-tx\"}{\"version\"}=\
\"1.0.2\";\n$PG{\"poa\"}{\"4_TCOFFEE\"}=\"POA\";\n\
$PG{\"poa\"}{\"type\"}=\"sequence_multiple_aligner\
\";\n$PG{\"poa\"}{\"ADDRESS\"}=\"http://www.bioinf\
ormatics.ucla.edu/poa/\";\n$PG{\"poa\"}{\"language\
\"}=\"C\";\n$PG{\"poa\"}{\"language2\"}=\"C\";\n$P\
G{\"poa\"}{\"source\"}=\"http://downloads.sourcefo\
rge.net/poamsa/poaV2.tar.gz\";\n$PG{\"poa\"}{\"DIR\
\"}=\"/usr/share/\";\n$PG{\"poa\"}{\"FILE1\"}=\"bl\
osum80.mat\";\n$PG{\"poa\"}{\"mode\"}=\"mcoffee\";\
\n$PG{\"poa\"}{\"binary\"}=\"poa\";\n$PG{\"poa\"}{\
\"version\"}=\"2.0\";\n$PG{\"probcons\"}{\"4_TCOFF\
EE\"}=\"PROBCONS\";\n$PG{\"probcons\"}{\"type\"}=\\
"sequence_multiple_aligner\";\n$PG{\"probcons\"}{\\
"ADDRESS\"}=\"http://probcons.stanford.edu/\";\n$P\
G{\"probcons\"}{\"language2\"}=\"CXX\";\n$PG{\"pro\
bcons\"}{\"language\"}=\"C++\";\n$PG{\"probcons\"}\
{\"source\"}=\"http://probcons.stanford.edu/probco\
ns_v1_12.tar.gz\";\n$PG{\"probcons\"}{\"mode\"}=\"\
mcoffee\";\n$PG{\"probcons\"}{\"binary\"}=\"probco\
ns\";\n$PG{\"probcons\"}{\"version\"}=\"1.12\";\n$\
PG{\"mafft\"}{\"4_TCOFFEE\"}=\"MAFFT\";\n$PG{\"maf\
ft\"}{\"type\"}=\"sequence_multiple_aligner\";\n$P\
G{\"mafft\"}{\"ADDRESS\"}=\"http://align.bmr.kyush\
u-u.ac.jp/mafft/online/server/\";\n$PG{\"mafft\"}{\
\"language\"}=\"C\";\n$PG{\"mafft\"}{\"language\"}\
=\"C\";\n$PG{\"mafft\"}{\"source\"}=\"http://align\
.bmr.kyushu-u.ac.jp/mafft/software/mafft-6.603-wit\
h-extensions-src.tgz\";\n$PG{\"mafft\"}{\"windows\\
"}=\"http://align.bmr.kyushu-u.ac.jp/mafft/softwar\
e/mafft-6.603-mingw.tar\";\n$PG{\"mafft\"}{\"mode\\
"}=\"mcoffee,rcoffee\";\n$PG{\"mafft\"}{\"binary\"\
}=\"mafft.tar.gz\";\n$PG{\"mafft\"}{\"version\"}=\\
"6.603\";\n$PG{\"muscle\"}{\"4_TCOFFEE\"}=\"MUSCLE\
\";\n$PG{\"muscle\"}{\"type\"}=\"sequence_multiple\
_aligner\";\n$PG{\"muscle\"}{\"ADDRESS\"}=\"http:/\
/www.drive5.com/muscle/\";\n$PG{\"muscle\"}{\"lang\
uage\"}=\"C++\";\n$PG{\"muscle\"}{\"language2\"}=\\
"GPP\";\n$PG{\"muscle\"}{\"source\"}=\"http://www.\
drive5.com/muscle/downloads3.7/muscle3.7_src.tar.g\
z\";\n$PG{\"muscle\"}{\"windows\"}=\"http://www.dr\
ive5.com/muscle/downloads3.7/muscle3.7_win32.zip\"\
;\n$PG{\"muscle\"}{\"linux\"}=\"http://www.drive5.\
com/muscle/downloads3.7/muscle3.7_linux_ia32.tar.g\
z\";\n$PG{\"muscle\"}{\"mode\"}=\"mcoffee,rcoffee\\
";\n$PG{\"muscle\"}{\"version\"}=\"3.7\";\n$PG{\"m\
us4\"}{\"4_TCOFFEE\"}=\"MUS4\";\n$PG{\"mus4\"}{\"t\
ype\"}=\"sequence_multiple_aligner\";\n$PG{\"mus4\\
"}{\"ADDRESS\"}=\"http://www.drive5.com/muscle/\";\
\n$PG{\"mus4\"}{\"language\"}=\"C++\";\n$PG{\"mus4\
\"}{\"language2\"}=\"GPP\";\n$PG{\"mus4\"}{\"sourc\
e\"}=\"http://www.drive5.com/muscle/muscle4.0_src.\
tar.gz\";\n$PG{\"mus4\"}{\"mode\"}=\"mcoffee,rcoff\
ee\";\n$PG{\"mus4\"}{\"version\"}=\"4.0\";\n$PG{\"\
pcma\"}{\"4_TCOFFEE\"}=\"PCMA\";\n$PG{\"pcma\"}{\"\
type\"}=\"sequence_multiple_aligner\";\n$PG{\"pcma\
\"}{\"ADDRESS\"}=\"ftp://iole.swmed.edu/pub/PCMA/\\
";\n$PG{\"pcma\"}{\"language\"}=\"C\";\n$PG{\"pcma\
\"}{\"language2\"}=\"C\";\n$PG{\"pcma\"}{\"source\\
"}=\"ftp://iole.swmed.edu/pub/PCMA/pcma.tar.gz\";\\
n$PG{\"pcma\"}{\"mode\"}=\"mcoffee\";\n$PG{\"pcma\\
"}{\"version\"}=\"1.0\";\n$PG{\"kalign\"}{\"4_TCOF\
FEE\"}=\"KALIGN\";\n$PG{\"kalign\"}{\"type\"}=\"se\
quence_multiple_aligner\";\n$PG{\"kalign\"}{\"ADDR\
ESS\"}=\"http://msa.cgb.ki.se\";\n$PG{\"kalign\"}{\
\"language\"}=\"C\";\n$PG{\"kalign\"}{\"language2\\
"}=\"C\";\n$PG{\"kalign\"}{\"source\"}=\"http://ms\
a.cgb.ki.se/downloads/kalign/current.tar.gz\";\n$P\
G{\"kalign\"}{\"mode\"}=\"mcoffee\";\n$PG{\"kalign\
\"}{\"version\"}=\"1.0\";\n$PG{\"amap\"}{\"4_TCOFF\
EE\"}=\"AMAP\";\n$PG{\"amap\"}{\"type\"}=\"sequenc\
e_multiple_aligner\";\n$PG{\"amap\"}{\"ADDRESS\"}=\
\"http://bio.math.berkeley.edu/amap/\";\n$PG{\"ama\
p\"}{\"language\"}=\"C++\";\n$PG{\"amap\"}{\"langu\
age2\"}=\"CXX\";\n$PG{\"amap\"}{\"source\"}=\"http\
://amap-align.googlecode.com/files/amap.2.0.tar.gz\
\";\n$PG{\"amap\"}{\"mode\"}=\"mcoffee\";\n$PG{\"a\
map\"}{\"version\"}=\"2.0\";\n$PG{\"proda\"}{\"4_T\
COFFEE\"}=\"PRODA\";\n$PG{\"proda\"}{\"type\"}=\"s\
equence_multiple_aligner\";\n$PG{\"proda\"}{\"ADDR\
ESS\"}=\"http://proda.stanford.edu\";\n$PG{\"proda\
\"}{\"language\"}=\"C++\";\n$PG{\"proda\"}{\"langu\
age2\"}=\"CXX\";\n$PG{\"proda\"}{\"source\"}=\"htt\
p://proda.stanford.edu/proda_1_0.tar.gz\";\n$PG{\"\
proda\"}{\"mode\"}=\"mcoffee\";\n$PG{\"proda\"}{\"\
version\"}=\"1.0\";\n$PG{\"fsa\"}{\"4_TCOFFEE\"}=\\
"FSA\";\n$PG{\"fsa\"}{\"type\"}=\"sequence_multipl\
e_aligner\";\n$PG{\"fsa\"}{\"ADDRESS\"}=\"http://f\
sa.sourceforge.net/\";\n$PG{\"fsa\"}{\"language\"}\
=\"C++\";\n$PG{\"fsa\"}{\"language2\"}=\"CXX\";\n$\
PG{\"fsa\"}{\"source\"}=\"http://sourceforge.net/p\
rojects/fsa/files/fsa-1.15.3.tar.gz/download/\";\n\
$PG{\"fsa\"}{\"mode\"}=\"mcoffee\";\n$PG{\"fsa\"}{\
\"version\"}=\"1.15.3\";\n$PG{\"prank\"}{\"4_TCOFF\
EE\"}=\"PRANK\";\n$PG{\"prank\"}{\"type\"}=\"seque\
nce_multiple_aligner\";\n$PG{\"prank\"}{\"ADDRESS\\
"}=\"http://www.ebi.ac.uk/goldman-srv/prank/\";\n$\
PG{\"prank\"}{\"language\"}=\"C++\";\n$PG{\"prank\\
"}{\"language2\"}=\"CXX\";\n$PG{\"prank\"}{\"sourc\
e\"}=\"http://www.ebi.ac.uk/goldman-srv/prank/src/\
prank/prank.src.100303.tgz\";\n$PG{\"prank\"}{\"mo\
de\"}=\"mcoffee\";\n$PG{\"prank\"}{\"version\"}=\"\
100303\";\n$PG{\"sap\"}{\"4_TCOFFEE\"}=\"SAP\";\n$\
PG{\"sap\"}{\"type\"}=\"structure_pairwise_aligner\
\";\n$PG{\"sap\"}{\"ADDRESS\"}=\"http://mathbio.ni\
mr.mrc.ac.uk/wiki/Software\";\n$PG{\"sap\"}{\"lang\
uage\"}=\"C\";\n$PG{\"sap\"}{\"language2\"}=\"C\";\
\n$PG{\"sap\"}{\"source\"}=\"http://mathbio.nimr.m\
rc.ac.uk/download/sap-1.1.1.tar.gz\";\n$PG{\"sap\"\
}{\"mode\"}=\"expresso,3dcoffee\";\n$PG{\"sap\"}{\\
"version\"}=\"1.1.1\";\n$PG{\"TMalign\"}{\"4_TCOFF\
EE\"}=\"TMALIGN\";\n$PG{\"TMalign\"}{\"type\"}=\"s\
tructure_pairwise_aligner\";\n$PG{\"TMalign\"}{\"A\
DDRESS\"}=\"http://zhang.bioinformatics.ku.edu/TM-\
align/TMalign.f\";\n$PG{\"TMalign\"}{\"language\"}\
=\"Fortran\";\n$PG{\"TMalign\"}{\"language2\"}=\"F\
ortran\";\n$PG{\"TMalign\"}{\"source\"}=\"http://z\
hang.bioinformatics.ku.edu/TM-align/TMalign.f\";\n\
$PG{\"TMalign\"}{\"linux\"}=\"http://zhang.bioinfo\
rmatics.ku.edu/TM-align/TMalign_32.gz\";\n$PG{\"TM\
align\"}{\"mode\"}=\"expresso,3dcoffee\";\n$PG{\"T\
Malign\"}{\"version\"}=\"1.0\";\n$PG{\"mustang\"}{\
\"4_TCOFFEE\"}=\"MUSTANG\";\n$PG{\"mustang\"}{\"ty\
pe\"}=\"structure_pairwise_aligner\";\n$PG{\"musta\
ng\"}{\"ADDRESS\"}=\"http://www.cs.mu.oz.au/~arun/\
mustang\";\n$PG{\"mustang\"}{\"language\"}=\"C++\"\
;\n$PG{\"mustang\"}{\"language2\"}=\"CXX\";\n$PG{\\
"mustang\"}{\"source\"}=\"http://ww2.cs.mu.oz.au/~\
arun/mustang/mustang_v3.2.1.tgz\";\n$PG{\"mustang\\
"}{\"mode\"}=\"expresso,3dcoffee\";\n$PG{\"mustang\
\"}{\"version\"}=\"3.2.1\";\n$PG{\"lsqman\"}{\"4_T\
COFFEE\"}=\"LSQMAN\";\n$PG{\"lsqman\"}{\"type\"}=\\
"structure_pairwise_aligner\";\n$PG{\"lsqman\"}{\"\
ADDRESS\"}=\"empty\";\n$PG{\"lsqman\"}{\"language\\
"}=\"empty\";\n$PG{\"lsqman\"}{\"language2\"}=\"em\
pty\";\n$PG{\"lsqman\"}{\"source\"}=\"empty\";\n$P\
G{\"lsqman\"}{\"update_action\"}=\"never\";\n$PG{\\
"lsqman\"}{\"mode\"}=\"expresso,3dcoffee\";\n$PG{\\
"align_pdb\"}{\"4_TCOFFEE\"}=\"ALIGN_PDB\";\n$PG{\\
"align_pdb\"}{\"type\"}=\"structure_pairwise_align\
er\";\n$PG{\"align_pdb\"}{\"ADDRESS\"}=\"empty\";\\
n$PG{\"align_pdb\"}{\"language\"}=\"empty\";\n$PG{\
\"align_pdb\"}{\"language2\"}=\"empty\";\n$PG{\"al\
ign_pdb\"}{\"source\"}=\"empty\";\n$PG{\"align_pdb\
\"}{\"update_action\"}=\"never\";\n$PG{\"align_pdb\
\"}{\"mode\"}=\"expresso,3dcoffee\";\n$PG{\"fuguea\
li\"}{\"4_TCOFFEE\"}=\"FUGUE\";\n$PG{\"fugueali\"}\
{\"type\"}=\"structure_pairwise_aligner\";\n$PG{\"\
fugueali\"}{\"ADDRESS\"}=\"http://www-cryst.bioc.c\
am.ac.uk/fugue/download.html\";\n$PG{\"fugueali\"}\
{\"language\"}=\"empty\";\n$PG{\"fugueali\"}{\"lan\
guage2\"}=\"empty\";\n$PG{\"fugueali\"}{\"source\"\
}=\"empty\";\n$PG{\"fugueali\"}{\"update_action\"}\
=\"never\";\n$PG{\"fugueali\"}{\"mode\"}=\"express\
o,3dcoffee\";\n$PG{\"dalilite.pl\"}{\"4_TCOFFEE\"}\
=\"DALILITEc\";\n$PG{\"dalilite.pl\"}{\"type\"}=\"\
structure_pairwise_aligner\";\n$PG{\"dalilite.pl\"\
}{\"ADDRESS\"}=\"built_in\";\n$PG{\"dalilite.pl\"}\
{\"ADDRESS2\"}=\"http://www.ebi.ac.uk/Tools/webser\
vices/services/dalilite\";\n$PG{\"dalilite.pl\"}{\\
"language\"}=\"Perl\";\n$PG{\"dalilite.pl\"}{\"lan\
guage2\"}=\"Perl\";\n$PG{\"dalilite.pl\"}{\"source\
\"}=\"empty\";\n$PG{\"dalilite.pl\"}{\"update_acti\
on\"}=\"never\";\n$PG{\"dalilite.pl\"}{\"mode\"}=\\
"expresso,3dcoffee\";\n$PG{\"probconsRNA\"}{\"4_TC\
OFFEE\"}=\"PROBCONSRNA\";\n$PG{\"probconsRNA\"}{\"\
type\"}=\"RNA_multiple_aligner\";\n$PG{\"probconsR\
NA\"}{\"ADDRESS\"}=\"http://probcons.stanford.edu/\
\";\n$PG{\"probconsRNA\"}{\"language\"}=\"C++\";\n\
$PG{\"probconsRNA\"}{\"language2\"}=\"CXX\";\n$PG{\
\"probconsRNA\"}{\"source\"}=\"http://probcons.sta\
nford.edu/probconsRNA.tar.gz\";\n$PG{\"probconsRNA\
\"}{\"mode\"}=\"mcoffee,rcoffee\";\n$PG{\"probcons\
RNA\"}{\"version\"}=\"1.0\";\n$PG{\"sfold\"}{\"4_T\
COFFEE\"}=\"CONSAN\";\n$PG{\"sfold\"}{\"type\"}=\"\
RNA_pairwise_aligner\";\n$PG{\"sfold\"}{\"ADDRESS\\
"}=\"http://selab.janelia.org/software/consan/\";\\
n$PG{\"sfold\"}{\"language\"}=\"empty\";\n$PG{\"sf\
old\"}{\"language2\"}=\"empty\";\n$PG{\"sfold\"}{\\
"source\"}=\"empty\";\n$PG{\"sfold\"}{\"update_act\
ion\"}=\"never\";\n$PG{\"sfold\"}{\"mode\"}=\"rcof\
fee\";\n$PG{\"RNAplfold\"}{\"4_TCOFFEE\"}=\"RNAPLF\
OLD\";\n$PG{\"RNAplfold\"}{\"type\"}=\"RNA_seconda\
rystructure_predictor\";\n$PG{\"RNAplfold\"}{\"ADD\
RESS\"}=\"http://www.tbi.univie.ac.at/~ivo/RNA/\";\
\n$PG{\"RNAplfold\"}{\"language\"}=\"C\";\n$PG{\"R\
NAplfold\"}{\"language2\"}=\"C\";\n$PG{\"RNAplfold\
\"}{\"source\"}=\"http://www.tbi.univie.ac.at/~ivo\
/RNA/ViennaRNA-1.7.2.tar.gz\";\n$PG{\"RNAplfold\"}\
{\"mode\"}=\"rcoffee,\";\n$PG{\"RNAplfold\"}{\"ver\
sion\"}=\"1.7.2\";\n$PG{\"retree\"}{\"4_TCOFFEE\"}\
=\"PHYLIP\";\n$PG{\"retree\"}{\"type\"}=\"RNA_seco\
ndarystructure_predictor\";\n$PG{\"retree\"}{\"ADD\
RESS\"}=\"http://evolution.gs.washington.edu/phyli\
p/\";\n$PG{\"retree\"}{\"language\"}=\"C\";\n$PG{\\
"retree\"}{\"language2\"}=\"C\";\n$PG{\"retree\"}{\
\"source\"}=\"http://evolution.gs.washington.edu/p\
hylip/download/phylip-3.69.tar.gz\";\n$PG{\"retree\
\"}{\"mode\"}=\"trmsd,\";\n$PG{\"retree\"}{\"versi\
on\"}=\"3.69\";\n$PG{\"hmmtop\"}{\"4_TCOFFEE\"}=\"\
HMMTOP\";\n$PG{\"hmmtop\"}{\"type\"}=\"protein_sec\
ondarystructure_predictor\";\n$PG{\"hmmtop\"}{\"AD\
DRESS\"}=\"www.enzim.hu/hmmtop/\";\n$PG{\"hmmtop\"\
}{\"language\"}=\"C\";\n$PG{\"hmmtop\"}{\"language\
2\"}=\"C\";\n$PG{\"hmmtop\"}{\"source\"}=\"empty\"\
;\n$PG{\"hmmtop\"}{\"update_action\"}=\"never\";\n\
$PG{\"hmmtop\"}{\"mode\"}=\"tcoffee\";\n$PG{\"gorI\
V\"}{\"4_TCOFFEE\"}=\"GOR4\";\n$PG{\"gorIV\"}{\"ty\
pe\"}=\"protein_secondarystructure_predictor\";\n$\
PG{\"gorIV\"}{\"ADDRESS\"}=\"http://mig.jouy.inra.\
fr/logiciels/gorIV/\";\n$PG{\"gorIV\"}{\"language\\
"}=\"C\";\n$PG{\"gorIV\"}{\"language2\"}=\"C\";\n$\
PG{\"gorIV\"}{\"source\"}=\"http://mig.jouy.inra.f\
r/logiciels/gorIV/GOR_IV.tar.gz\";\n$PG{\"gorIV\"}\
{\"update_action\"}=\"never\";\n$PG{\"gorIV\"}{\"m\
ode\"}=\"tcoffee\";\n$PG{\"wublast.pl\"}{\"4_TCOFF\
EE\"}=\"EBIWUBLASTc\";\n$PG{\"wublast.pl\"}{\"type\
\"}=\"protein_homology_predictor\";\n$PG{\"wublast\
.pl\"}{\"ADDRESS\"}=\"built_in\";\n$PG{\"wublast.p\
l\"}{\"ADDRESS2\"}=\"http://www.ebi.ac.uk/Tools/we\
bservices/services/wublast\";\n$PG{\"wublast.pl\"}\
{\"language\"}=\"Perl\";\n$PG{\"wublast.pl\"}{\"la\
nguage2\"}=\"Perl\";\n$PG{\"wublast.pl\"}{\"source\
\"}=\"empty\";\n$PG{\"wublast.pl\"}{\"update_actio\
n\"}=\"never\";\n$PG{\"wublast.pl\"}{\"mode\"}=\"p\
sicoffee,expresso,accurate\";\n$PG{\"blastpgp.pl\"\
}{\"4_TCOFFEE\"}=\"EBIBLASTPGPc\";\n$PG{\"blastpgp\
.pl\"}{\"type\"}=\"protein_homology_predictor\";\n\
$PG{\"blastpgp.pl\"}{\"ADDRESS\"}=\"built_in\";\n$\
PG{\"blastpgp.pl\"}{\"ADDRESS2\"}=\"http://www.ebi\
.ac.uk/Tools/webservices/services/blastpgp\";\n$PG\
{\"blastpgp.pl\"}{\"language\"}=\"Perl\";\n$PG{\"b\
lastpgp.pl\"}{\"language2\"}=\"Perl\";\n$PG{\"blas\
tpgp.pl\"}{\"source\"}=\"empty\";\n$PG{\"blastpgp.\
pl\"}{\"update_action\"}=\"never\";\n$PG{\"blastpg\
p.pl\"}{\"mode\"}=\"psicoffee,expresso,accurate\";\
\n$PG{\"blastcl3\"}{\"4_TCOFFEE\"}=\"NCBIWEBBLAST\\
";\n$PG{\"blastcl3\"}{\"type\"}=\"protein_homology\
_predictor\";\n$PG{\"blastcl3\"}{\"ADDRESS\"}=\"ft\
p://ftp.ncbi.nih.gov/blast/executables/LATEST\";\n\
$PG{\"blastcl3\"}{\"language\"}=\"C\";\n$PG{\"blas\
tcl3\"}{\"language2\"}=\"C\";\n$PG{\"blastcl3\"}{\\
"source\"}=\"empty\";\n$PG{\"blastcl3\"}{\"update_\
action\"}=\"never\";\n$PG{\"blastcl3\"}{\"mode\"}=\
\"psicoffee,expresso,3dcoffee\";\n$PG{\"blastpgp\"\
}{\"4_TCOFFEE\"}=\"NCBIBLAST\";\n$PG{\"blastpgp\"}\
{\"type\"}=\"protein_homology_predictor\";\n$PG{\"\
blastpgp\"}{\"ADDRESS\"}=\"ftp://ftp.ncbi.nih.gov/\
blast/executables/LATEST\";\n$PG{\"blastpgp\"}{\"l\
anguage\"}=\"C\";\n$PG{\"blastpgp\"}{\"language2\"\
}=\"C\";\n$PG{\"blastpgp\"}{\"source\"}=\"empty\";\
\n$PG{\"blastpgp\"}{\"update_action\"}=\"never\";\\
n$PG{\"blastpgp\"}{\"mode\"}=\"psicoffee,expresso,\
3dcoffee\";\n$PG{\"SOAP::Lite\"}{\"4_TCOFFEE\"}=\"\
SOAPLITE\";\n$PG{\"SOAP::Lite\"}{\"type\"}=\"libra\
ry\";\n$PG{\"SOAP::Lite\"}{\"ADDRESS\"}=\"http://c\
pansearch.perl.org/src/MKUTTER/SOAP-Lite-0.710.08/\
Makefile.PL\";\n$PG{\"SOAP::Lite\"}{\"language\"}=\
\"Perl\";\n$PG{\"SOAP::Lite\"}{\"language2\"}=\"Pe\
rl\";\n$PG{\"SOAP::Lite\"}{\"source\"}=\"empty\";\\
n$PG{\"blastpgp\"}{\"update_action\"}=\"never\";\n\
$PG{\"SOAP::Lite\"}{\"mode\"}=\"none\";\n$PG{\"XML\
::Simple\"}{\"4_TCOFFEE\"}=\"XMLSIMPLE\";\n$PG{\"X\
ML::Simple\"}{\"type\"}=\"library\";\n$PG{\"XML::S\
imple\"}{\"ADDRESS\"}=\"http://search.cpan.org/~gr\
antm/XML-Simple-2.18/lib/XML/Simple.pm\";\n$PG{\"X\
ML::Simple\"}{\"language\"}=\"Perl\";\n$PG{\"XML::\
Simple\"}{\"language2\"}=\"Perl\";\n$PG{\"XML::Sim\
ple\"}{\"source\"}=\"empty\";\n$PG{\"XML::Simple\"\
}{\"mode\"}=\"psicoffee,expresso,accurate\";\n$MOD\
E{\"tcoffee\"}{\"name\"}=\"tcoffee\";\n$MODE{\"rco\
ffee\"}{\"name\"}=\"rcoffee\";\n$MODE{\"3dcoffee\"\
}{\"name\"}=\"3dcoffee\";\n$MODE{\"mcoffee\"}{\"na\
me\"}=\"mcoffee\";\n$MODE{\"expresso\"}{\"name\"}=\
\"expresso\";\n$MODE{\"trmsd\"}{\"name\"}=\"trmsd\\
";\n$MODE{\"accurate\"}{\"name\"}=\"accurate\";\n$\
MODE{\"seq_reformat\"}{\"name\"}=\"seq_reformat\";\
\n\n\n$PG{C}{compiler}=\"gcc\";\n$PG{C}{compiler_f\
lag}=\"CC\";\n$PG{C}{options}=\"\";\n$PG{C}{option\
s_flag}=\"CFLAGS\";\n$PG{C}{type}=\"compiler\";\n\\
n$PG{\"CXX\"}{compiler}=\"g++\";\n$PG{\"CXX\"}{com\
piler_flag}=\"CXX\";\n$PG{\"CXX\"}{options}=\"\";\\
n$PG{\"CXX\"}{options_flag}=\"CXXFLAGS\";\n$PG{CXX\
}{type}=\"compiler\";\n\n$PG{\"CPP\"}{compiler}=\"\
g++\";\n$PG{\"CPP\"}{compiler_flag}=\"CPP\";\n$PG{\
\"CPP\"}{options}=\"\";\n$PG{\"CPP\"}{options_flag\
}=\"CPPFLAGS\";\n$PG{CPP}{type}=\"compiler\";\n\n$\
PG{\"GPP\"}{compiler}=\"g++\";\n$PG{\"GPP\"}{compi\
ler_flag}=\"GPP\";\n$PG{\"GPP\"}{options}=\"\";\n$\
PG{\"GPP\"}{options_flag}=\"CFLAGS\";\n$PG{GPP}{ty\
pe}=\"compiler\";\n\n$PG{Fortran}{compiler}=\"g77\\
";\n$PG{Fortran}{compiler_flag}=\"FCC\";\n$PG{Fort\
ran}{type}=\"compiler\";\n\n$PG{Perl}{compiler}=\"\
CPAN\";\n$PG{Perl}{type}=\"compiler\";\n\n$SUPPORT\
ED_OS{macox}=\"Macintosh\";\n$SUPPORTED_OS{linux}=\
\"Linux\";\n$SUPPORTED_OS{windows}=\"Cygwin\";\n\n\
\n\n$MODE{t_coffee}{description}=\" for regular mu\
ltiple sequence alignments\";\n$MODE{rcoffee} {des\
cription}=\" for RNA multiple sequence alignments\\
";\n\n$MODE{psicoffee} {description}=\" for Homolo\
gy Extended multiple sequence alignments\";\n$MODE\
{expresso}{description}=\" for very accurate struc\
ture based multiple sequence alignments\";\n$MODE{\
\"3dcoffee\"}{description}=\" for multiple structu\
re alignments\";\n$MODE{mcoffee} {description}=\" \
for combining alternative multiple sequence alignm\
ent packages\\n------- into a unique meta-package.\
 The installer will upload several MSA packages an\
d compile them\\n\n\";\n\n\n&post_process_PG();\nr\
eturn;\n}\n\nsub post_process_PG\n  {\n    my $p;\\
n    \n    %PG=&name2dname (%PG);\n    %MODE=&name\
2dname(%MODE);\n    foreach $p (keys(%PG)){if ( $P\
G{$p}{type} eq \"compiler\"){$PG{$p}{update_action\
}=\"never\";}}\n    \n  }\n\nsub name2dname\n  {\n\
    my (%L)=(@_);\n    my ($l, $ml);\n    \n    fo\
reach my $pg (keys(%L))\n      {\n	$l=length ($pg)\
;\n	if ( $l>$ml){$ml=$l;}\n      }\n    $ml+=1;\n \
   foreach my $pg (keys(%L))\n      {\n	my $name;\\
n	$l=$ml-length ($pg);\n	$name=$pg;\n	for ( $b=0; \
$b<$l; $b++)\n	  {\n	    $name .=\" \";\n	  }\n	$L\
{$pg}{dname}=$name;\n      }\n    return %L;\n  }\\
n\nsub env_file2putenv\n  {\n    my $f=@_[0];\n   \
 my $F=new FileHandle;\n    my $n;\n    \n    open\
 ($F, \"$f\");\n    while (<$F>)\n      {\n	my $li\
ne=$_;\n	my($var, $value)=($_=~/(\\S+)\\=(\\S*)/);\
\n	$ENV{$var}=$value;\n	$ENV_SET{$var}=1;\n	$n++;\\
n      }\n    close ($F);\n    return $n;\n  }\n\n\
sub replace_line_in_file\n  {\n    my ($file, $wor\
din, $wordout)=@_;\n    my $O=new FileHandle;\n   \
 my $I=new FileHandle;\n    my $l;\n    if (!-e $f\
ile){return;}\n    \n    system (\"mv $file $file.\
old\");\n    open ($O, \">$file\");\n    open ($I,\
 \"$file.old\");\n    while (<$I>)\n      {\n	$l=$\
_;\n	if (!($l=~/$wordin/)){print $O \"$l\";}\n	els\
if ( $wordout ne \"\"){$l=~s/$wordin/$wordout/g;pr\
int $O \"$l\";}\n      }\n    close ($O);\n    clo\
se ($I);\n    return;\n  }\n\nsub add_C_libraries\\
n  {\n   my ($file,$first,@list)=@_;\n   \n    my \
$O=new FileHandle;\n    my $I=new FileHandle;\n   \
 my ($l,$anchor);\n    if (!-e $file){return;}\n  \
 \n    $anchor=\"#include <$first>\";\n	 \n    sys\
tem (\"mv $file $file.old\");\n    open ($O, \">$f\
ile\");\n    open ($I, \"$file.old\");\n    while \
(<$I>)\n      {\n	$l=$_;\n	print $O \"$l\";\n	if (\
!($l=~/$anchor/))\n	   {\n	    \n	    foreach my $\
lib (@list)\n	       {\n                  print $O\
 \"#include <$lib>\\n\";\n	       }\n           }\\
n      }\n    close ($O);\n    close ($I);\n    re\
turn;\n    }\n","use Env;\nuse Cwd;\n@suffix=(\"tm\
p\", \"temp\", \"cache\", \"t_coffee\", \"core\", \
\"tcoffee\");\n\nif ($#ARGV==-1)\n  {\n    print \\
"clean_cache.pl -file <file to add in -dir> -dir=<\
dir> -size=<value in Mb>\\n0: unlimited -1 always.\
\\nWill only clean directories matching:[\";\n    \
foreach $k(@suffix){print \"*$k* \";}\n    print \\
"]\\n\";\n    exit (EXIT_FAILURE);\n  }\n\n$cl=joi\
n (\" \",@ARGV);\nif (($cl=~/\\-no_action/))\n  {\\
n    exit (EXIT_SUCCESS);\n  }\n\nif (($cl=~/\\-de\
bug/))\n  {\n    $DEBUG=1;\n  }\nelse\n  {\n    $D\
EBUG=0;\n  }\n\nif (($cl=~/\\-dir=(\\S+)/))\n  {\n\
    $dir=$1;\n  }\nelse\n  {\n    $dir=\"./\";\n  \
}\n\nif ($cl=~/\\-file=(\\S+)/)\n  {\n    $file=$1\
;\n  }\nelse\n  {\n    $file=0;\n  }\n\nif ($cl=~/\
\\-size=(\\S+)/)\n  {\n    $max_size=$1;\n  }\nels\
e\n  {\n    $max_size=0;#unlimited\n  }\nif ($cl=~\
/\\-force/)\n  {\n    $force=1;\n  }\nelse\n  {\n \
   $force=0;\n  }\n\nif ($cl=~/\\-age=(\\S+)/)\n  \
{\n    $max_age=$1;\n  }\nelse\n  {\n    $max_age=\
0;#unlimited\n  }\n\n$max_size*=1000000;\nif ( ! -\
d $dir)\n  {\n    print STDERR \"\\nCannot process\
 $dir: does not exist \\n\";\n    exit (EXIT_FAILU\
RE);\n  }\n\nif ( !($dir=~/^\\//))\n  {\n    $base\
=cwd();\n    $dir=\"$base/$dir\";\n  }\n\n$proceed\
=0;\nforeach $s (@suffix)\n  {\n    \n    if (($di\
r=~/$s/)){$proceed=1;}\n    $s=uc ($s);\n    if ((\
$dir=~/$s/)){$proceed=1;}\n  }\nif ( $proceed==0)\\
n  {\n    print STDERR \"Clean_cache.pl can only c\
lean directories whose absolute path name contains\
 the following strings:\";\n    foreach $w (@suffi\
x) {print STDERR \"$w \";$w=lc($w); print STDERR \\
"$w \";}\n    print STDERR \"\\nCannot process $di\
r\\n\";\n    exit (EXIT_FAILURE);\n  }\n\n$name_fi\
le=\"$dir/name_file.txt\";\n$size_file=\"$dir/size\
_file.txt\";\nif ( $force){&create_ref_file ($dir,\
$name_file,$size_file);}\nif ($file){&add_file ($d\
ir, $name_file, $size_file, $file);}\n&clean_dir (\
$dir, $name_file, $size_file, $max_size,$max_age);\
\nexit (EXIT_SUCCESS);\n\nsub clean_dir \n  {\n   \
 my ($dir, $name_file, $size_file, $max_size, $max\
_age)=@_;\n    my ($tot_size, $size, $f, $s);\n\n \
 \n    $tot_size=&get_tot_size ($dir, $name_file, \
$size_file);\n\n    if ( $tot_size<=$max_size){ret\
urn ;}\n    else {$max_size/=2;}\n    \n    #recre\
ate the name file in case some temprary files have\
 not been properly registered\n    &create_ref_fil\
e ($dir, $name_file, $size_file, $max_age);\n  \n \
   $new_name_file=&vtmpnam();\n    open (R, \"$nam\
e_file\");\n    open (W, \">$new_name_file\");\n  \
  while (<R>)\n      {\n	my $line=$_;\n	\n	($f, $s\
)=($line=~/(\\S+) (\\S+)/);\n	if ( !($f=~/\\S/)){n\
ext;}\n	\n	elsif ($max_size && $tot_size>=$max_siz\
e && !($f=~/name_file/))\n	  {\n	    remove ( \"$d\
ir/$f\");\n	    $tot_size-=$s;\n	  }\n	elsif ( $ma\
x_age && -M(\"$dir/$f\")>=$max_age)\n	  {\n	    re\
move ( \"$dir/$f\");\n	    $tot_size-=$s;\n	  }\n	\
else\n	  {\n	    print W \"$f $s\\n\";\n	  }\n    \
  }\n    close (R);\n    close (W);\n    open (F, \
\">$size_file\");\n    print F \"$tot_size\";\n   \
 if ( -e $new_name_file){`mv $new_name_file $name_\
file`;}\n    close (F);\n  }\nsub get_tot_size\n  \
{\n    my ($dir, $name_file, $size_file)=@_;\n    \
my $size;\n    \n    if ( !-d $dir){return 0;}\n  \
  if ( !-e $name_file)\n      {\n	\n	&create_ref_f\
ile ($dir, $name_file, $size_file);\n      }\n    \
open (F, \"$size_file\");\n    $size=<F>;\n    clo\
se (F);\n    chomp ($size);\n    return $size;\n  \
}\nsub size \n  {\n    my $f=@_[0];\n\n    if ( !-\
d $f){return -s($f);}\n    else {return &dir2size(\
$f);}\n  }\nsub dir2size\n  {\n    my $d=@_[0];\n \
   my ($s, $f);\n    \n    if ( !-d $d) {return 0;\
}\n    \n    foreach $f (&dir2list ($d))\n      {\\
n	if ( -d $f){$s+=&dir2size (\"$d/$f\");}\n	else {\
$s+= -s \"$dir/$f\";}\n      }\n    return $s;\n  \
}\n\nsub remove \n  {\n    my $file=@_[0];\n    my\
 ($f);\n    \n    debug_print( \"--- $file ---\\n\\
");\n    if (($file eq \".\") || ($file eq \"..\")\
 || ($file=~/\\*/)){return EXIT_FAILURE;}\n    els\
if ( !-d $file)\n      {\n	debug_print (\"unlink $\
file\\n\");\n	if (-e $file){unlink ($file);}\n    \
  }\n    elsif ( -d $file)\n      {\n	debug_print \
(\"++++++++ $file +++++++\\n\");\n	foreach $f (&di\
r2list($file))\n	  {\n	    &remove (\"$file/$f\");\
\n	  }\n	debug_print (\"rmdir $file\\n\");\n	rmdir\
 $file;\n      }\n    else\n      {\n	debug_print \
(\"????????? $file ????????\\n\");\n      }\n    r\
eturn EXIT_SUCCESS;\n  }\n\nsub dir2list\n  {\n   \
 my $dir=@_[0];\n    my (@list1, @list2,@list3, $l\
);\n\n    opendir (DIR,$dir);\n    @list1=readdir \
(DIR);\n    closedir (DIR);\n    \n    foreach $l \
(@list1)\n      {\n	if ( $l ne \".\" && $l ne \"..\
\"){@list2=(@list2, $l);}\n      }\n    @list3 = s\
ort { (-M \"$dir/$list2[$b]\") <=> (-M \"$dir/$lis\
t2[$a]\")} @list2;\n    return @list3;\n    \n  }\\
n\nsub debug_print\n  {\n    \n    if ($DEBUG==1){\
print @_;}\n    \n  }\nsub create_ref_file\n  {\n \
   my ($dir,$name_file,$size_file)=@_;\n    my ($f\
, $s, $tot_size, @l);\n    \n    if ( !-d $dir){re\
turn;}\n    \n    @l=&dir2list ($dir);\n    open (\
F, \">$name_file\");\n    foreach $f (@l)\n      {\
\n	$s=&size(\"$dir/$f\");\n	$tot_size+=$s;\n	print\
 F \"$f $s\\n\";\n      }\n    &myecho ($tot_size,\
 \">$size_file\");\n    close (F);\n  }\nsub add_f\
ile \n  {\n    my ($dir,$name_file,$size_file,$fil\
e)=@_;\n    my ($s, $tot_size);\n    \n    if ( !-\
d $dir)   {return;}\n    if ( !-e \"$dir/$file\" )\
 {return;}\n    if ( !-e $name_file){&create_ref_f\
ile ($dir,$name_file,$size_file);}\n					    \n   \
 $s=&size(\"$dir/$file\");\n    open (F, \">>$name\
_file\");\n    print F \"$file\\n\";\n    close (F\
);\n\n    $tot_size=&get_tot_size ($dir,$name_file\
,$size_file);\n    $tot_size+=$s;\n    &myecho ($t\
ot_size, \">$size_file\");\n    \n  }\n	\nsub myec\
ho\n  {\n    my ($string, $file)=@_;\n    open (EC\
HO, $file) || die;\n    print ECHO \"$string\";\n \
   close (ECHO);\n  }\n    \n		\n	\nsub vtmpnam\n \
 {\n    my $tmp_file_name;\n    $tmp_name_counter+\
+;\n    $tmp_file_name=\"tmp_file_for_clean_cache_\
pdb$$.$tmp_name_counter\";\n    $tmp_file_list[$nt\
mp_file++]=$tmp_file_name;\n    if ( -e $tmp_file_\
name) {return &vtmpnam ();}\n    else {return $tmp\
_file_name;}\n  }\n","\n$t_coffee=\"t_coffee\";\n\\
nforeach $value ( @ARGV)\n  {\n    $seq_file=$seq_\
file.\" \".$value;\n  }\n\n$name=$ARGV[0];\n$name=\
~s/\\.[^\\.]*$//;\n$lib_name=\"$name.mocca_lib\";\\
n$type=`t_coffee $seq_file -get_type -quiet`;\ncho\
p ($type);\n\nif ( $type eq \"PROTEIN\"){$lib_mode\
=\"lalign_rs_s_pair -lalign_n_top 20\";}\nelsif ( \
$type eq\"DNA\"){$lib_mode=\"lalign_rs_s_dna_pair \
-lalign_n_top 40\";}\n\nif ( !(-e $lib_name))\n  {\
\n	  \n  $command=\"$t_coffee -mocca -seq_weight=n\
o -cosmetic_penalty=0 -mocca_interactive -in $lib_\
mode -out_lib $lib_name -infile $seq_file\";\n  \n\
  }\nelsif ( (-e $lib_name))\n  {\n  $command=\"$t\
_coffee -mocca -seq_weight=no -cosmetic_penalty=0 \
-mocca_interactive -in $lib_name -infile $seq_file\
\";\n  \n  }\n\nsystem ($command);\n\nexit;\n\n","\
my $WSDL = 'http://www.ebi.ac.uk/Tools/webservices\
/wsdl/WSDaliLite.wsdl';\n\nuse SOAP::Lite;\nuse Da\
ta::Dumper;\nuse Getopt::Long qw(:config no_ignore\
_case bundling);\nuse File::Basename;\n\nmy $check\
Interval = 5;\n\nmy %params=(\n	    'async' => '1'\
, # Use async mode and simulate sync mode in clien\
t\n	    );\nGetOptions(\n    'pdb1=s'     => \\$pa\
rams{'sequence1'},\n    'chainid1=s' => \\$params{\
'chainid1'},\n    'pdb2=s'     => \\$params{'seque\
nce2'},\n    'chainid2=s' => \\$params{'chainid2'}\
,\n    \"help|h\"	 => \\$help, # Usage info\n    \\
"async|a\"	 => \\$async, # Asynchronous submission\
\n    \"polljob\"	 => \\$polljob, # Get results\n \
   \"status\"	 => \\$status, # Get status\n    \"j\
obid|j=s\"  => \\$jobid, # JobId\n    \"email|S=s\\
"  => \\$params{email}, # E-mail address\n    \"tr\
ace\"      => \\$trace, # SOAP messages\n    \"seq\
uence=s\" => \\$sequence, # Input PDB\n    );\n\nm\
y $scriptName = basename($0, ());\nif($help) {\n  \
  &usage();\n    exit(0);\n}\n\nif($trace) {\n    \
print \"Tracing active\\n\";\n    SOAP::Lite->impo\
rt(+trace => 'debug');\n}\n\nmy $soap = SOAP::Lite\
\n    ->service($WSDL)\n    ->on_fault(sub {\n    \
    my $soap = shift;\n        my $res = shift;\n \
       # Throw an exception for all faults\n      \
  if(ref($res) eq '') {\n            die($res);\n \
       } else {\n            die($res->faultstring\
);\n        }\n        return new SOAP::SOM;\n    \
}\n               );\n\nif( !($polljob || $status)\
 &&\n    !( defined($params{'sequence1'}) && defin\
ed($params{'sequence2'}) )\n    ) {\n    print STD\
ERR 'Error: bad option combination', \"\\n\";\n   \
 &usage();\n    exit(1);\n}\nelsif($polljob && def\
ined($jobid)) {\n    print \"Getting results for j\
ob $jobid\\n\";\n    getResults($jobid);\n}\nelsif\
($status && defined($jobid)) {\n    print STDERR \\
"Getting status for job $jobid\\n\";\n    my $resu\
lt = $soap->checkStatus($jobid);\n    print STDOUT\
 \"$result\", \"\\n\";\n    if($result eq 'DONE') \
{\n	print STDERR \"To get results: $scriptName --p\
olljob --jobid $jobid\\n\";\n    }\n}\nelse {\n   \
 if(-f $params{'sequence1'}) {\n	$params{'sequence\
1'} = read_file($params{'sequence1'});\n    }\n   \
 if(-f $params{'sequence2'}) {\n	$params{'sequence\
2'} = read_file($params{'sequence2'});\n    }\n\n \
   my $jobid;\n    my $paramsData = SOAP::Data->na\
me('params')->type(map=>\\%params);\n    # For SOA\
P::Lite 0.60 and earlier parameters are passed dir\
ectly\n    if($SOAP::Lite::VERSION eq '0.60' || $S\
OAP::Lite::VERSION =~ /0\\.[1-5]/) {\n        $job\
id = $soap->runDaliLite($paramsData);\n    }\n    \
# For SOAP::Lite 0.69 and later parameter handling\
 is different, so pass\n    # undef's for template\
d params, and then pass the formatted args.\n    e\
lse {\n        $jobid = $soap->runDaliLite(undef,\\
n				     $paramsData);\n    }\n\n    if (defined(\
$async)) {\n	print STDOUT $jobid, \"\\n\";\n      \
  print STDERR \"To check status: $scriptName --st\
atus --jobid $jobid\\n\";\n    } else { # Synchron\
ous mode\n        print STDERR \"JobId: $jobid\\n\\
";\n        sleep 1;\n        getResults($jobid);\\
n    }\n}\n\nsub clientPoll($) {\n    my $jobid = \
shift;\n    my $result = 'PENDING';\n    # Check s\
tatus and wait if not finished\n    #print STDERR \
\"Checking status: $jobid\\n\";\n    while($result\
 eq 'RUNNING' || $result eq 'PENDING') {\n        \
$result = $soap->checkStatus($jobid);\n        pri\
nt STDERR \"$result\\n\";\n        if($result eq '\
RUNNING' || $result eq 'PENDING') {\n            #\
 Wait before polling again.\n            sleep $ch\
eckInterval;\n        }\n    }\n}\n\nsub getResult\
s($) {\n    $jobid = shift;\n    # Check status, a\
nd wait if not finished\n    clientPoll($jobid);\n\
    # Use JobId if output file name is not defined\
\n    unless(defined($outfile)) {\n        $outfil\
e=$jobid;\n    }\n    # Get list of data types\n  \
  my $resultTypes = $soap->getResults($jobid);\n  \
  # Get the data and write it to a file\n    if(de\
fined($outformat)) { # Specified data type\n      \
  my $selResultType;\n        foreach my $resultTy\
pe (@$resultTypes) {\n            if($resultType->\
{type} eq $outformat) {\n                $selResul\
tType = $resultType;\n            }\n        }\n  \
      $res=$soap->poll($jobid, $selResultType->{ty\
pe});\n        write_file($outfile.'.'.$selResultT\
ype->{ext}, $res);\n    } else { # Data types avai\
lable\n        # Write a file for each output type\
\n        for my $resultType (@$resultTypes){\n   \
         #print \"Getting $resultType->{type}\\n\"\
;\n            $res=$soap->poll($jobid, $resultTyp\
e->{type});\n            write_file($outfile.'.'.$\
resultType->{ext}, $res);\n        }\n    }\n}\n\n\
sub read_file($) {\n    my $filename = shift;\n   \
 open(FILE, $filename);\n    my $content;\n    my \
$buffer;\n    while(sysread(FILE, $buffer, 1024)) \
{\n	$content.= $buffer;\n    }\n    close(FILE);\n\
    return $content;\n}\n\nsub write_file($$) {\n \
   my ($tmp,$entity) = @_;\n    print STDERR \"Cre\
ating result file: \".$tmp.\"\\n\";\n    unless(op\
en (FILE, \">$tmp\")) {\n	return 0;\n    }\n    sy\
swrite(FILE, $entity);\n    close (FILE);\n    ret\
urn 1;\n}\n\nsub usage {\n    print STDERR <<EOF\n\
DaliLite\n========\n\nPairwise comparison of prote\
in structures\n\n[Required]\n\n  --pdb1           \
     : str  : PDB ID for structure 1\n  --pdb2    \
            : str  : PDB ID for structure 2\n\n[Op\
tional]\n\n  --chain1              : str  : Chain \
identifer in structure 1\n  --chain2              \
: str  : Chain identifer in structure 2\n\n[Genera\
l]\n\n  -h, --help            :      : prints this\
 help text\n  -S, --email           : str  : user \
email address\n  -a, --async           :      : as\
ynchronous submission\n      --status          :  \
    : poll for the status of a job\n      --polljo\
b         :      : poll for the results of a job\n\
  -j, --jobid           : str  : jobid for an asyn\
chronous job\n  -O, --outfile         : str  : fil\
e name for results (default is jobid)\n      --tra\
ce	        :      : show SOAP messages being inter\
changed \n\nSynchronous job:\n\n  The results/erro\
rs are returned as soon as the job is finished.\n \
 Usage: $scriptName --email <your\\@email> [option\
s] pdbFile [--outfile string]\n  Returns: saves th\
e results to disk\n\nAsynchronous job:\n\n  Use th\
is if you want to retrieve the results at a later \
time. The results \n  are stored for up to 24 hour\
s. \n  The asynchronous submission mode is recomme\
nded when users are submitting \n  batch jobs or l\
arge database searches	\n  Usage: $scriptName --em\
ail <your\\@email> --async [options] pdbFile\n  Re\
turns: jobid\n\n  Use the jobid to query for the s\
tatus of the job. \n  Usage: $scriptName --status \
--jobid <jobId>\n  Returns: string indicating the \
status of the job:\n    DONE - job has finished\n \
   RUNNING - job is running\n    NOT_FOUND - job c\
annot be found\n    ERROR - the jobs has encounter\
ed an error\n\n  When done, use the jobid to retri\
eve the status of the job. \n  Usage: $scriptName \
--polljob --jobid <jobId> [--outfile string]\n\n[H\
elp]\n\n  For more detailed help information refer\
 to\n  http://www.ebi.ac.uk/DaliLite/\nEOF\n;\n}\n\
","my $WSDL = 'http://www.ebi.ac.uk/Tools/webservi\
ces/wsdl/WSWUBlast.wsdl';\n\nuse strict;\nuse SOAP\
::Lite;\nuse Getopt::Long qw(:config no_ignore_cas\
e bundling);\nuse File::Basename;\n\nmy $checkInte\
rval = 15;\n\nmy $numOpts = scalar(@ARGV);\nmy ($o\
utfile, $outformat, $help, $async, $polljob, $stat\
us, $ids, $jobid, $trace, $sequence);\nmy %params=\
 ( # Defaults\n	      'async' => 1, # Force into a\
sync mode\n	      'exp' => 10.0, # E-value thresho\
ld\n	      'numal' => 50, # Maximum number of alig\
nments\n	      'scores' => 100, # Maximum number o\
f scores\n            );\nGetOptions( # Map the op\
tions into variables\n    \"program|p=s\"     => \\
\$params{program}, # BLAST program\n    \"database\
|D=s\"    => \\$params{database}, # Search databas\
e\n    \"matrix|m=s\"      => \\$params{matrix}, #\
 Scoring matrix\n    \"exp|E=f\"         => \\$par\
ams{exp}, # E-value threshold\n    \"echofilter|e\\
"    => \\$params{echofilter}, # Display filtered \
sequence\n    \"filter|f=s\"      => \\$params{fil\
ter}, # Low complexity filter name\n    \"alignmen\
ts|b=i\"  => \\$params{numal}, # Number of alignme\
nts\n    \"scores|s=i\"      => \\$params{scores},\
 # Number of scores\n    \"sensitivity|S=s\" => \\\
$params{sensitivity}, # Search sensitivity\n    \"\
sort|t=s\"	      => \\$params{sort}, # Sort hits b\
y...\n    \"stats|T=s\"       => \\$params{stats},\
 # Scoring statistic to use\n    \"strand|d=s\"   \
   => \\$params{strand}, # Strand to use in DNA vs\
. DNA search\n    \"topcombon|c=i\"   => \\$params\
{topcombon}, # Consistent sets of HSPs\n    \"outf\
ile=s\"       => \\$outfile, # Output file\n    \"\
outformat|o=s\"   => \\$outformat, # Output format\
\n    \"help|h\"	      => \\$help, # Usage info\n \
   \"async|a\"	      => \\$async, # Asynchronous m\
ode\n    \"polljob\"	      => \\$polljob, # Get re\
sults\n    \"status\"	      => \\$status, # Get jo\
b status\n    \"ids\"             => \\$ids, # Get\
 ids from result\n    \"jobid|j=s\"       => \\$jo\
bid, # JobId\n    \"email=s\"         => \\$params\
{email}, # E-mail address\n    \"trace\"          \
 => \\$trace, # SOAP trace\n    \"sequence=s\"    \
  => \\$sequence, # Query sequence\n    );\n\nmy $\
scriptName = basename($0, ());\nif($help || $numOp\
ts == 0) {\n    &usage();\n    exit(0);\n}\n\nif($\
trace){\n    print STDERR \"Tracing active\\n\";\n\
    SOAP::Lite->import(+trace => 'debug');\n}\n\nm\
y $soap = SOAP::Lite\n    ->service($WSDL)\n    ->\
proxy('http://localhost/',\n    #proxy => ['http' \
=> 'http://your.proxy.server/'], # HTTP proxy\n   \
 timeout => 600, # HTTP connection timeout\n    )\\
n    ->on_fault(sub { # SOAP fault handler\n      \
  my $soap = shift;\n        my $res = shift;\n   \
     # Throw an exception for all faults\n        \
if(ref($res) eq '') {\n            die($res);\n   \
     } else {\n            die($res->faultstring);\
\n        }\n        return new SOAP::SOM;\n    }\\
n               );\n\nif( !($polljob || $status ||\
 $ids) &&\n    !( defined($ARGV[0]) || defined($se\
quence) )\n    ) {\n    print STDERR 'Error: bad o\
ption combination', \"\\n\";\n    &usage();\n    e\
xit(1);\n}\nelsif($polljob && defined($jobid)) {\n\
    print \"Getting results for job $jobid\\n\";\n\
    getResults($jobid);\n}\nelsif($status && defin\
ed($jobid)) {\n    print STDERR \"Getting status f\
or job $jobid\\n\";\n    my $result = $soap->check\
Status($jobid);\n    print STDOUT \"$result\\n\";\\
n    if($result eq 'DONE') {\n	print STDERR \"To g\
et results: $scriptName --polljob --jobid $jobid\\\
n\";\n    }\n}  \nelsif($ids && defined($jobid)) {\
\n    print STDERR \"Getting ids from job $jobid\\\
n\";\n    getIds($jobid);\n}\nelse {\n    # Prepar\
e input data\n    my $content;\n    my (@contents)\
 = ();\n    if(-f $ARGV[0] || $ARGV[0] eq '-') {	\\
n	$content={type=>'sequence',content=>read_file($A\
RGV[0])};	\n    }\n    if($sequence) {	\n	if(-f $s\
equence || $sequence eq '-') {	\n	    $content={ty\
pe=>'sequence',content=>read_file($ARGV[0])};	\n	}\
 else {\n	    $content={type=>'sequence',content=>\
$sequence};\n	}\n    }\n    push @contents, $conte\
nt;\n\n    # Submit the job\n    my $paramsData = \
SOAP::Data->name('params')->type(map=>\\%params);\\
n    my $contentData = SOAP::Data->name('content')\
->value(\\@contents);\n    # For SOAP::Lite 0.60 a\
nd earlier parameters are passed directly\n    if(\
$SOAP::Lite::VERSION eq '0.60' || $SOAP::Lite::VER\
SION =~ /0\\.[1-5]/) {\n        $jobid = $soap->ru\
nWUBlast($paramsData, $contentData);\n    }\n    #\
 For SOAP::Lite 0.69 and later parameter handling \
is different, so pass\n    # undef's for templated\
 params, and then pass the formatted args.\n    el\
se {\n        $jobid = $soap->runWUBlast(undef, un\
def,\n				   $paramsData, $contentData);\n    }\n\\
n    # Asynchronous mode: output jobid and exit.\n\
    if (defined($async)) {\n	print STDOUT $jobid, \
\"\\n\";\n        print STDERR \"To check status: \
$scriptName --status --jobid $jobid\\n\";\n    }\n\
    # Synchronous mode: try to get results\n    el\
se {\n        print STDERR \"JobId: $jobid\\n\";\n\
        sleep 1;\n        getResults($jobid);\n   \
 }\n}\n\nsub getIds($) {\n    my $jobid = shift;\n\
    my $results = $soap->getIds($jobid);\n    for \
my $result (@$results){\n	print \"$result\\n\";\n \
   }\n}\n\nsub clientPoll($) {\n    my $jobid = sh\
ift;\n    my $result = 'PENDING';\n    # Check sta\
tus and wait if not finished\n    while($result eq\
 'RUNNING' || $result eq 'PENDING') {\n        $re\
sult = $soap->checkStatus($jobid);\n        print \
STDERR \"$result\\n\";\n        if($result eq 'RUN\
NING' || $result eq 'PENDING') {\n            # Wa\
it before polling again.\n            sleep $check\
Interval;\n        }\n    }\n}\n\nsub getResults($\
) {\n    my $jobid = shift;\n    my $res;\n    # C\
heck status, and wait if not finished\n    clientP\
oll($jobid);\n    # Use JobId if output file name \
is not defined\n    unless(defined($outfile)) {\n \
       $outfile=$jobid;\n    }\n    # Get list of \
data types\n    my $resultTypes = $soap->getResult\
s($jobid);\n    # Get the data and write it to a f\
ile\n    if(defined($outformat)) { # Specified dat\
a type\n	if($outformat eq 'xml') {$outformat = 'to\
olxml';}\n	if($outformat eq 'txt') {$outformat = '\
tooloutput';}\n        my $selResultType;\n       \
 foreach my $resultType (@$resultTypes) {\n       \
     if($resultType->{type} eq $outformat) {\n    \
            $selResultType = $resultType;\n       \
     }\n        }\n        $res=$soap->poll($jobid\
, $selResultType->{type});\n	if($outfile eq '-') {\
\n	     write_file($outfile, $res);\n	} else {\n	 \
   write_file($outfile.'.'.$selResultType->{ext}, \
$res);\n	}\n    } else { # Data types available\n \
       # Write a file for each output type\n      \
  for my $resultType (@$resultTypes){\n           \
 #print STDERR \"Getting $resultType->{type}\\n\";\
\n            $res=$soap->poll($jobid, $resultType\
->{type});\n	    if($outfile eq '-') {\n		write_fi\
le($outfile, $res);\n	    } else {\n		write_file($\
outfile.'.'.$resultType->{ext}, $res);\n	    }\n  \
      }\n    }\n}\n\nsub read_file($) {\n    my $f\
ilename = shift;\n    my ($content, $buffer);\n   \
 if($filename eq '-') {\n	while(sysread(STDIN, $bu\
ffer, 1024)) {\n	    $content .= $buffer;\n	}\n   \
 }\n    else { # File\n	open(FILE, $filename) or d\
ie \"Error: unable to open input file\";\n	while(s\
ysread(FILE, $buffer, 1024)) {\n	    $content .= $\
buffer;\n	}\n	close(FILE);\n    }\n    return $con\
tent;\n}\n\nsub write_file($$) {\n    my ($filenam\
e, $data) = @_;\n    print STDERR 'Creating result\
 file: ' . $filename . \"\\n\";\n    if($filename \
eq '-') {\n	print STDOUT $data;\n    }\n    else {\
\n	open(FILE, \">$filename\") or die \"Error: unab\
le to open output file\";\n	syswrite(FILE, $data);\
\n	close(FILE);\n    }\n}\n\nsub usage {\n    prin\
t STDERR <<EOF\nWU-BLAST\n========\n\nRapid sequen\
ce database search programs utilizing the BLAST al\
gorithm.\n   \n[Required]\n\n      --email       :\
 str  : user email address \n  -p, --program	    :\
 str  : BLAST program to use: blastn, blastp, blas\
tx, \n                             tblastn or tbla\
stx\n  -D, --database    : str  : database to sear\
ch\n  seqFile           : file : query sequence da\
ta file (\"-\" for STDIN)\n\n[Optional]\n\n  -m, -\
-matrix	    : str  : scoring matrix\n  -E, --exp	 \
   : real : 0<E<= 1000. Statistical significance t\
hreshold\n                             for reporti\
ng database sequence matches.\n  -e, --echofilter \
 :      : display the filtered query sequence in t\
he output\n  -f, --filter	    : str  : activates f\
iltering of the query sequence\n  -b, --alignments\
  : int  : number of alignments to be reported\n  \
-s, --scores	    : int  : number of scores to be r\
eported\n  -S, --sensitivity : str  :\n  -t, --sor\
t	    : str  :\n  -T, --stats       : str  :\n  -d\
, --strand      : str  : DNA strand to search with\
 in DNA vs. DNA searches \n  -c, --topcombon   :  \
    :\n\n[General]	\n\n  -h, --help       :      :\
 prints this help text\n  -a, --async      :      \
: forces to make an asynchronous query\n      --st\
atus     :      : poll for the status of a job\n  \
    --polljob    :      : poll for the results of \
a job\n  -j, --jobid      : str  : jobid that was \
returned when an asynchronous job \n              \
              was submitted.\n  -O, --outfile    :\
 str  : name of the file results should be written\
 to \n                            (default is base\
d on the jobid; \"-\" for STDOUT)\n  -o, --outform\
at  : str  : txt or xml output (no file is written\
)\n      --trace	   :      : show SOAP messages be\
ing interchanged \n\nSynchronous job:\n\n  The res\
ults/errors are returned as soon as the job is fin\
ished.\n  Usage: $scriptName --email <your\\@email\
> [options...] seqFile\n  Returns: saves the resul\
ts to disk\n\nAsynchronous job:\n\n  Use this if y\
ou want to retrieve the results at a later time. T\
he results \n  are stored for up to 24 hours. \n  \
The asynchronous submission mode is recommended wh\
en users are submitting \n  batch jobs or large da\
tabase searches	\n  Usage: $scriptName --async --e\
mail <your\\@email> [options...] seqFile\n  Return\
s : jobid\n\n  Use the jobid to query for the stat\
us of the job. \n  Usage: $scriptName --status --j\
obid <jobId>\n  Returns : string indicating the st\
atus of the job:\n    DONE - job has finished\n   \
 RUNNING - job is running\n    NOT_FOUND - job can\
not be found\n    ERROR - the jobs has encountered\
 an error\n\n  When done, use the jobid to retriev\
e the status of the job. \n  Usage: $scriptName --\
polljob --jobid <jobId> [--outfile string]\n  Retu\
rns: saves the results to disk\n\n[Help]\n\nFor mo\
re detailed help information refer to \nhttp://www\
.ebi.ac.uk/blast2/WU-Blast2_Help_frame.html\n \nEO\
F\n;\n}\n","\nmy $WSDL = 'http://www.ebi.ac.uk/Too\
ls/webservices/wsdl/WSBlastpgp.wsdl';\n\nuse SOAP:\
:Lite;\nuse Getopt::Long qw(:config no_ignore_case\
 bundling);\nuse File::Basename;\n\nmy $checkInter\
val = 15;\n\nmy %params=(\n	    'async' => '1', # \
Use async mode and simulate sync mode in client\n	\
    );\nGetOptions(\n    \"mode=s\"           => \\
\$params{mode}, # Search mode: PSI-Blast or PHI-Bl\
ast\n    \"database|d=s\"     => \\$params{databas\
e}, # Database to search\n    \"matrix|M=s\"      \
 => \\$params{matrix},# Scoring maxtrix\n    \"exp\
|e=f\"          => \\$params{exp}, # E-value\n    \
\"expmulti|h=f\"     => \\$params{expmulti}, # E-v\
alue\n    \"filter|F=s\"       => \\$params{filter\
}, # Low complexity filter\n    \"dropoff|X=i\"   \
   => \\$params{dropoff}, # Dropoff score\n    \"f\
inaldropoff|Z=i\" => \\$params{finaldropoff}, # Fi\
nal dropoff score\n    \"scores|v=i\"       => \\$\
params{scores}, # Max number of scores\n    \"alig\
n=i\"          => \\$params{align}, # Alignment vi\
ew\n    \"startregion|S=i\"  => \\$params{startreg\
ion}, # Start of region in query\n    \"endregion|\
H=i\"    => \\$params{endregion}, # End of region \
in query\n    \"maxpasses|j=i\"    => \\$params{ma\
xpasses}, # Number of PSI iterations\n    \"openga\
p|G=i\"      => \\$params{opengap}, # Gap open pen\
alty\n    \"extendgap|E=i\"    => \\$params{extend\
gap}, # Gap extension penalty\n    \"pattern=s\"  \
      => \\$params{pattern}, # PHI-BLAST pattern\n\
    \"usagemode|p=s\"    => \\$params{usagemode}, \
# PHI-BLAST program\n    \"appxml=s\"         => \\
\$params{appxml}, # Application XML\n    \"sequenc\
e=s\"       => \\$sequence, # Query sequence\n    \
\"help\"	       => \\$help, # Usage info\n    \"po\
lljob\"	       => \\$polljob, # Get results\n    \\
"status\"	       => \\$status, # Get status\n    \\
"ids\"      	       => \\$ids, # Get ids from resu\
lt\n    \"jobid=s\"          => \\$jobid, # JobId\\
n    \"outfile=s\"        => \\$outfile, # Output \
filename\n    \"outformat|o=s\"    => \\$outformat\
, # Output file format\n    \"async|a\"	       => \
\\$async, # Async submission\n    \"email=s\"     \
     => \\$params{email}, # User e-mail address\n \
   \"trace\"            => \\$trace, # Show SOAP m\
essages\n    );\n\nmy $scriptName = basename($0, (\
));\nif($help) {\n    &usage();\n    exit(0);\n}\n\
\nif ($trace){\n    print \"Tracing active\\n\";\n\
    SOAP::Lite->import(+trace => 'debug');\n}\n\nm\
y $soap = SOAP::Lite\n    ->service($WSDL)\n    ->\
on_fault(sub {\n        my $soap = shift;\n       \
 my $res = shift;\n        # Throw an exception fo\
r all faults\n        if(ref($res) eq '') {\n     \
       die($res);\n        } else {\n            d\
ie($res->faultstring);\n        }\n        return \
new SOAP::SOM;\n    }\n               );\n\nif( !(\
$polljob || $status || $ids) &&\n    !( (defined($\
ARGV[0]) && -f $ARGV[0]) || defined($sequence) )\n\
    ) {\n    print STDERR 'Error: bad option combi\
nation', \"\\n\";\n    &usage();\n    exit(1);\n}\\
nelsif($polljob && defined($jobid)) {\n    print \\
"Getting results for job $jobid\\n\";\n    getResu\
lts($jobid);\n}\nelsif($status && defined($jobid))\
 {\n    print STDERR \"Getting status for job $job\
id\\n\";\n    my $result = $soap->checkStatus($job\
id);\n    print STDOUT $result, \"\\n\";\n    if($\
result eq 'DONE') {\n	print STDERR \"To get result\
s: $scriptName --polljob --jobid $jobid\\n\";\n   \
 }\n}  \nelsif($ids && defined($jobid)) {\n    pri\
nt STDERR \"Getting ids from job $jobid\\n\";\n   \
 getIds($jobid);\n}\nelse {\n    if(-f $ARGV[0]) {\
	\n	$content={type=>'sequence', content=>read_file\
($ARGV[0])};	\n    }\n    if($sequence) {	\n	if(-f\
 $sequence) {\n	    $content={type=>'sequence', co\
ntent=>read_file($sequence)};	\n	} else {\n	    $c\
ontent={type=>'sequence', content=>$sequence};\n	}\
\n    }\n    push @content, $content;\n\n    my $j\
obid;\n    my $paramsData = SOAP::Data->name('para\
ms')->type(map=>\\%params);\n    my $contentData =\
 SOAP::Data->name('content')->value(\\@content);\n\
    # For SOAP::Lite 0.60 and earlier parameters a\
re passed directly\n    if($SOAP::Lite::VERSION eq\
 '0.60' || $SOAP::Lite::VERSION =~ /0\\.[1-5]/) {\\
n        $jobid = $soap->runBlastpgp($paramsData, \
$contentData);\n    }\n    # For SOAP::Lite 0.69 a\
nd later parameter handling is different, so pass\\
n    # undef's for templated params, and then pass\
 the formatted args.\n    else {\n        $jobid =\
 $soap->runBlastpgp(undef, undef,\n				    $params\
Data, $contentData);\n    }\n\n    if (defined($as\
ync)) {\n	print STDOUT $jobid, \"\\n\";\n        p\
rint STDERR \"To check status: $scriptName --statu\
s --jobid $jobid\\n\";\n    } else { # Synchronous\
 mode\n        print STDERR \"JobId: $jobid\\n\";\\
n        sleep 1;\n        getResults($jobid);\n  \
  }\n}\n\nsub getIds($) {\n    $jobid = shift;\n  \
  my $results = $soap->getIds($jobid);\n    for $r\
esult (@$results){\n	print \"$result\\n\";\n    }\\
n}\n\nsub clientPoll($) {\n    my $jobid = shift;\\
n    my $result = 'PENDING';\n    # Check status a\
nd wait if not finished\n    #print STDERR \"Check\
ing status: $jobid\\n\";\n    while($result eq 'RU\
NNING' || $result eq 'PENDING') {\n        $result\
 = $soap->checkStatus($jobid);\n        print STDE\
RR \"$result\\n\";\n        if($result eq 'RUNNING\
' || $result eq 'PENDING') {\n            # Wait b\
efore polling again.\n            sleep $checkInte\
rval;\n        }\n    }\n}\n\nsub getResults($) {\\
n    $jobid = shift;\n    # Check status, and wait\
 if not finished\n    clientPoll($jobid);\n    # U\
se JobId if output file name is not defined\n    u\
nless(defined($outfile)) {\n        $outfile=$jobi\
d;\n    }\n    # Get list of data types\n    my $r\
esultTypes = $soap->getResults($jobid);\n    # Get\
 the data and write it to a file\n    if(defined($\
outformat)) { # Specified data type\n        my $s\
elResultType;\n        foreach my $resultType (@$r\
esultTypes) {\n            if($resultType->{type} \
eq $outformat) {\n                $selResultType =\
 $resultType;\n            }\n        }\n        $\
res=$soap->poll($jobid, $selResultType->{type});\n\
        write_file($outfile.'.'.$selResultType->{e\
xt}, $res);\n    } else { # Data types available\n\
        # Write a file for each output type\n     \
   for my $resultType (@$resultTypes){\n          \
  #print \"Getting $resultType->{type}\\n\";\n    \
        $res=$soap->poll($jobid, $resultType->{typ\
e});\n            write_file($outfile.'.'.$resultT\
ype->{ext}, $res);\n        }\n    }\n}\n\nsub rea\
d_file($) {\n    my $filename = shift;\n    open(F\
ILE, $filename);\n    my $content;\n    my $buffer\
;\n    while(sysread(FILE, $buffer, 1024)) {\n	$co\
ntent.= $buffer;\n    }\n    close(FILE);  \n    r\
eturn $content;\n}\n\nsub write_file($$) {\n    my\
 ($tmp,$entity) = @_;\n    print STDERR \"Creating\
 result file: \".$tmp.\"\\n\";\n    unless(open (F\
ILE, \">$tmp\")) {\n	return 0;\n    }\n    syswrit\
e(FILE, $entity);\n    close (FILE);\n    return 1\
;\n}\n\nsub usage {\n    print STDERR <<EOF\nBlast\
pgp\n========\n   \nThe blastpgp program implement\
s the PSI-BLAST and PHI-BLAST variations\nof NCBI \
BLAST.\n\nFor more detailed help information refer\
 to\nhttp://www.ebi.ac.uk/blastpgp/blastpsi_help_f\
rame.html\n \nBlastpgp specific options:\n\n[Requi\
red]\n\n      --mode            : str  : search mo\
de to use: PSI-Blast or PHI-Blast\n  -d, --databas\
e        : str  : protein database to search\n  se\
qFile               : file : query sequence\n\n[Op\
tional]\n\n  -M, --matrix          : str  : scorin\
g matrix\n  -e, --exp             : real : Expecta\
tion value\n  -h, --expmulti        : real : thres\
hold (multipass model)\n  -F, --filter          : \
str  : filter query sequence with SEG [T,F]\n  -m,\
 --align           : int  : alignment view option:\
\n                                 0 - pairwise, 1\
 - M/S identities,\n                              \
   2 - M/S non-identities, 3 - Flat identities,\n \
                                4 - Flat non-ident\
ities\n  -G, --opengap         : int  : cost to op\
en a gap\n  -E, --extendgap       : int  : cost to\
 extend a gap\n  -g, --gapalign        : str  : Ga\
pped [T,F]\n  -v, --scores          : int  : numbe\
r of scores to be reported\n  -j, --maxpasses     \
  : int  : number of iterations\n  -X, --dropoff  \
       : int  : Dropoff score\n  -Z, --finaldropof\
f    : int  : Dropoff for final alignment\n  -S, -\
-startregion     : int  : Start of required region\
 in query\n  -H, --endregion       : int  : End of\
 required region in query\n  -k, --pattern        \
 : str  : Hit File (PHI-BLAST only)\n  -p, --usage\
mode       : str  : Program option (PHI-BLAST only\
):\n                                 blastpgp, pat\
seedp, seedp\n\n[General]\n\n      --help         \
   :      : prints this help text\n  -a, --async  \
         :      : forces to make an asynchronous q\
uery\n      --status          :      : poll for th\
e status of a job\n      --polljob         :      \
: poll for the results of a job\n      --jobid    \
       : str  : jobid of an asynchronous job\n    \
  --ids             :      : get hit identifiers f\
or result \n  -O, --outfile         : str  : name \
of the file results should be written to\n        \
                         (default is based on the \
jobid)\n  -o, --outformat       : str  : txt or xm\
l output (no file is written)\n      --trace      \
     :      : show SOAP messages being interchange\
d\n\nSynchronous job:\n\n  The results/errors are \
returned as soon as the job is finished.\n  Usage:\
 blastpgp.pl --email <your@email> [options...] seq\
file\n  Returns: saves the results to disk\n\nAsyn\
chronous job:\n\n  Use this if you want to retriev\
e the results at a later time. The results\n  are \
stored for up to 24 hours.\n  The asynchronous sub\
mission mode is recommended when users are submitt\
ing\n  batch jobs or large database searches\n  Us\
age: blastpgp.pl --email <your@email> --async [opt\
ions...] seqFile\n  Returns: jobid\n\n  Use the jo\
bid to query for the status of the job.\n  Usage: \
blastpgp.pl --status --jobid <jobId>\n  Returns: s\
tring indicating the status of the job\n    DONE -\
 job has finished\n    RUNNING - job is running\n \
   NOT_FOUND - job cannot be found\n    ERROR - th\
e jobs has encountered an error\n\n  When done, us\
e the jobid to retrieve the results of the job.\n \
 Usage: blastpgp.pl --polljob --jobid <jobId> [--o\
utfile <fileName>]\n  Returns: saves the results t\
o disk\nEOF\n;\n}\n","\n=head1 NAME\n\nncbiblast_l\
wp.pl\n\n=head1 DESCRIPTION\n\nNCBI BLAST REST web\
 service Perl client using L<LWP>.\n\nTested with:\
\n\n=over\n\n=item *\nL<LWP> 5.79, L<XML::Simple> \
2.12 and Perl 5.8.3\n\n=item *\nL<LWP> 5.805, L<XM\
L::Simple> 2.14 and Perl 5.8.7\n\n=item *\nL<LWP> \
5.820, L<XML::Simple> 2.18 and Perl 5.10.0 (Ubuntu\
 9.04)\n\n=back\n\nFor further information see:\n\\
n=over\n\n=item *\nL<http://www.ebi.ac.uk/Tools/we\
bservices/services/sss/ncbi_blast_rest>\n\n=item *\
\nL<http://www.ebi.ac.uk/Tools/webservices/tutoria\
ls/perl>\n\n=back\n\n=head1 VERSION\n\n$Id: ncbibl\
ast_lwp.pl 1317 2009-09-03 15:44:11Z hpm $\n\n=cut\
\n\nuse strict;\nuse warnings;\n\nuse English;\nus\
e LWP;\nuse XML::Simple;\nuse Getopt::Long qw(:con\
fig no_ignore_case bundling);\nuse File::Basename;\
\nuse Data::Dumper;\n\nmy $baseUrl = 'http://www.e\
bi.ac.uk/Tools/services/rest/ncbiblast';\n\nmy $ch\
eckInterval = 3;\n\nmy $outputLevel = 1;\n\nmy $nu\
mOpts = scalar(@ARGV);\nmy %params = ( 'debugLevel\
' => 0 );\n\nmy %tool_params = ();\nGetOptions(\n\\
n	# Tool specific options\n	'program|p=s'  => \\$t\
ool_params{'program'},   # blastp, blastn, blastx,\
 etc.\n	'database|D=s' => \\$params{'database'},  \
     # Database(s) to search\n	'matrix|m=s'   => \\
\$tool_params{'matrix'},    # Scoring martix to us\
e\n	'exp|E=f'      => \\$tool_params{'exp'},      \
 # E-value threshold\n	'filter|f=s'   => \\$tool_p\
arams{'filter'},    # Low complexity filter\n	'ali\
gn|A=i'    => \\$tool_params{'align'},     # Pairw\
ise alignment format\n	'scores|s=i'   => \\$tool_p\
arams{'scores'},    # Number of scores\n	'alignmen\
ts|n=i' => \\$tool_params{'alignments'},   # Numbe\
r of alignments\n	'dropoff|d=i'    => \\$tool_para\
ms{'dropoff'},      # Dropoff score\n	'match_score\
s=s' => \\$tool_params{'match_scores'}, # Match/mi\
ssmatch scores\n	'match|u=i'      => \\$params{'ma\
tch'},             # Match score\n	'mismatch|v=i' \
  => \\$params{'mismatch'},          # Mismatch sc\
ore\n	'gapopen|o=i'    => \\$tool_params{'gapopen'\
},      # Open gap penalty\n	'gapext|x=i'     => \\
\$tool_params{'gapext'},       # Gap extension pen\
ality\n	'gapalign|g'     => \\$tool_params{'gapali\
gn'},     # Optimise gap alignments\n	'stype=s' =>\
 \\$tool_params{'stype'},    # Sequence type\n	'se\
qrange=s' => \\$tool_params{'seqrange'},    # Quer\
y subsequence\n	'sequence=s' => \\$params{'sequenc\
e'},         # Query sequence\n	'multifasta' => \\\
$params{'multifasta'},       # Multiple fasta inpu\
t\n\n	# Compatability options, old command-line\n	\
'numal|n=i'     => \\$params{'numal'},        # Nu\
mber of alignments\n	'opengap|o=i'   => \\$params{\
'opengap'},      # Open gap penalty\n	'extendgap|x\
=i' => \\$params{'extendgap'},    # Gap extension \
penality\n	\n	# Generic options\n	'email=s'       \
=> \\$params{'email'},          # User e-mail addr\
ess\n	'title=s'       => \\$params{'title'},      \
    # Job title\n	'outfile=s'     => \\$params{'ou\
tfile'},        # Output file name\n	'outformat=s'\
   => \\$params{'outformat'},      # Output file t\
ype\n	'jobid=s'       => \\$params{'jobid'},      \
    # JobId\n	'help|h'        => \\$params{'help'}\
,           # Usage help\n	'async'         => \\$p\
arams{'async'},          # Asynchronous submission\
\n	'polljob'       => \\$params{'polljob'},       \
 # Get results\n	'resultTypes'   => \\$params{'res\
ultTypes'},    # Get result types\n	'status'      \
  => \\$params{'status'},         # Get status\n	'\
params'        => \\$params{'params'},         # L\
ist input parameters\n	'paramDetail=s' => \\$param\
s{'paramDetail'},    # Get details for parameter\n\
	'quiet'         => \\$params{'quiet'},          #\
 Decrease output level\n	'verbose'       => \\$par\
ams{'verbose'},        # Increase output level\n	'\
debugLevel=i'  => \\$params{'debugLevel'},     # D\
ebug output level\n	'baseUrl=s'     => \\$baseUrl,\
                  # Base URL for service.\n);\nif \
( $params{'verbose'} ) { $outputLevel++ }\nif ( $p\
arams{'$quiet'} )  { $outputLevel-- }\n\n&print_de\
bug_message( 'MAIN', 'LWP::VERSION: ' . $LWP::VERS\
ION,\n	1 );\n\n&print_debug_message( 'MAIN', \"par\
ams:\\n\" . Dumper( \\%params ),           11 );\n\
&print_debug_message( 'MAIN', \"tool_params:\\n\" \
. Dumper( \\%tool_params ), 11 );\n\nmy $scriptNam\
e = basename( $0, () );\n\nif ( $params{'help'} ||\
 $numOpts == 0 ) {\n	&usage();\n	exit(0);\n}\n\n&p\
rint_debug_message( 'MAIN', 'baseUrl: ' . $baseUrl\
, 1 );\n\nif (\n	!(\n		   $params{'polljob'}\n		||\
 $params{'resultTypes'}\n		|| $params{'status'}\n	\
	|| $params{'params'}\n		|| $params{'paramDetail'}\
\n	)\n	&& !( defined( $ARGV[0] ) || defined( $para\
ms{'sequence'} ) )\n  )\n{\n\n	# Bad argument comb\
ination, so print error message and usage\n	print \
STDERR 'Error: bad option combination', \"\\n\";\n\
	&usage();\n	exit(1);\n}\n\nelsif ( $params{'param\
s'} ) {\n	&print_tool_params();\n}\n\nelsif ( $par\
ams{'paramDetail'} ) {\n	&print_param_details( $pa\
rams{'paramDetail'} );\n}\n\nelsif ( $params{'stat\
us'} && defined( $params{'jobid'} ) ) {\n	&print_j\
ob_status( $params{'jobid'} );\n}\n\nelsif ( $para\
ms{'resultTypes'} && defined( $params{'jobid'} ) )\
 {\n	&print_result_types( $params{'jobid'} );\n}\n\
\nelsif ( $params{'polljob'} && defined( $params{'\
jobid'} ) ) {\n	&get_results( $params{'jobid'} );\\
n}\n\nelse {\n\n	# Multiple input sequence mode, a\
ssume fasta format.\n	if ( $params{'multifasta'} )\
 {\n		&multi_submit_job();\n	}\n\n	# Entry identif\
ier list file.\n	elsif (( defined( $params{'sequen\
ce'} ) && $params{'sequence'} =~ m/^\\@/ )\n		|| (\
 defined( $ARGV[0] ) && $ARGV[0] =~ m/^\\@/ ) )\n	\
{\n		my $list_filename = $params{'sequence'} || $A\
RGV[0];\n		$list_filename =~ s/^\\@//;\n		&list_fi\
le_submit_job($list_filename);\n	}\n\n	# Default: \
single sequence/identifier.\n	else {\n\n		# Load t\
he sequence data and submit.\n		&submit_job( &load\
_data() );\n	}\n}\n\n=head1 FUNCTIONS\n\n=cut\n\n\\
n=head2 rest_request()\n\nPerform a REST request.\\
n\n  my $response_str = &rest_request($url);\n\n=c\
ut\n\nsub rest_request {\n	print_debug_message( 'r\
est_request', 'Begin', 11 );\n	my $requestUrl = sh\
ift;\n	print_debug_message( 'rest_request', 'URL: \
' . $requestUrl, 11 );\n\n	# Create a user agent\n\
	my $ua = LWP::UserAgent->new();\n	'$Revision: 131\
7 $' =~ m/(\\d+)/;\n	$ua->agent(\"EBI-Sample-Clien\
t/$1 ($scriptName; $OSNAME) \" . $ua->agent());\n	\
$ua->env_proxy;\n\n	# Perform the request\n	my $re\
sponse = $ua->get($requestUrl);\n	print_debug_mess\
age( 'rest_request', 'HTTP status: ' . $response->\
code,\n		11 );\n\n	# Check for HTTP error codes\n	\
if ( $response->is_error ) {\n		$response->content\
() =~ m/<h1>([^<]+)<\\/h1>/;\n		die 'http status: \
' . $response->code . ' ' . $response->message . '\
  ' . $1;\n	}\n	print_debug_message( 'rest_request\
', 'End', 11 );\n\n	# Return the response data\n	r\
eturn $response->content();\n}\n\n=head2 rest_get_\
parameters()\n\nGet list of tool parameter names.\\
n\n  my (@param_list) = &rest_get_parameters();\n\\
n=cut\n\nsub rest_get_parameters {\n	print_debug_m\
essage( 'rest_get_parameters', 'Begin', 1 );\n	my \
$url                = $baseUrl . '/parameters/';\n\
	my $param_list_xml_str = rest_request($url);\n	my\
 $param_list_xml     = XMLin($param_list_xml_str);\
\n	my (@param_list)       = @{ $param_list_xml->{'\
id'} };\n	print_debug_message( 'rest_get_parameter\
s', 'End', 1 );\n	return (@param_list);\n}\n\n=hea\
d2 rest_get_parameter_details()\n\nGet details of \
a tool parameter.\n\n  my $paramDetail = &rest_get\
_parameter_details($param_name);\n\n=cut\n\nsub re\
st_get_parameter_details {\n	print_debug_message( \
'rest_get_parameter_details', 'Begin', 1 );\n	my $\
parameterId = shift;\n	print_debug_message( 'rest_\
get_parameter_details',\n		'parameterId: ' . $para\
meterId, 1 );\n	my $url                  = $baseUr\
l . '/parameterdetails/' . $parameterId;\n	my $par\
am_detail_xml_str = rest_request($url);\n	my $para\
m_detail_xml     = XMLin($param_detail_xml_str);\n\
	print_debug_message( 'rest_get_parameter_details'\
, 'End', 1 );\n	return ($param_detail_xml);\n}\n\n\
=head2 rest_run()\n\nSubmit a job.\n\n  my $job_id\
 = &rest_run($email, $title, \\%params );\n\n=cut\\
n\nsub rest_run {\n	print_debug_message( 'rest_run\
', 'Begin', 1 );\n	my $email  = shift;\n	my $title\
  = shift;\n	my $params = shift;\n	print_debug_mes\
sage( 'rest_run', 'email: ' . $email, 1 );\n	if ( \
defined($title) ) {\n		print_debug_message( 'rest_\
run', 'title: ' . $title, 1 );\n	}\n	print_debug_m\
essage( 'rest_run', 'params: ' . Dumper($params), \
1 );\n\n	# User agent to perform http requests\n	m\
y $ua = LWP::UserAgent->new();\n	$ua->env_proxy;\n\
\n	# Clean up parameters\n	my (%tmp_params) = %{$p\
arams};\n	$tmp_params{'email'} = $email;\n	$tmp_pa\
rams{'title'} = $title;\n	foreach my $param_name (\
 keys(%tmp_params) ) {\n		if ( !defined( $tmp_para\
ms{$param_name} ) ) {\n			delete $tmp_params{$para\
m_name};\n		}\n	}\n\n	# Submit the job as a POST\n\
	my $url = $baseUrl . '/run';\n	my $response = $ua\
->post( $url, \\%tmp_params );\n	print_debug_messa\
ge( 'rest_run', 'HTTP status: ' . $response->code,\
 11 );\n	print_debug_message( 'rest_run',\n		'requ\
est: ' . $response->request()->content(), 11 );\n\\
n	# Check for HTTP error codes\n	if ( $response->i\
s_error ) {\n		$response->content() =~ m/<h1>([^<]\
+)<\\/h1>/;\n		die 'http status: ' . $response->co\
de . ' ' . $response->message . '  ' . $1;\n	}\n\n\
	# The job id is returned\n	my $job_id = $response\
->content();\n	print_debug_message( 'rest_run', 'E\
nd', 1 );\n	return $job_id;\n}\n\n=head2 rest_get_\
status()\n\nCheck the status of a job.\n\n  my $st\
atus = &rest_get_status($job_id);\n\n=cut\n\nsub r\
est_get_status {\n	print_debug_message( 'rest_get_\
status', 'Begin', 1 );\n	my $job_id = shift;\n	pri\
nt_debug_message( 'rest_get_status', 'jobid: ' . $\
job_id, 2 );\n	my $status_str = 'UNKNOWN';\n	my $u\
rl        = $baseUrl . '/status/' . $job_id;\n	$st\
atus_str = &rest_request($url);\n	print_debug_mess\
age( 'rest_get_status', 'status_str: ' . $status_s\
tr, 2 );\n	print_debug_message( 'rest_get_status',\
 'End', 1 );\n	return $status_str;\n}\n\n=head2 re\
st_get_result_types()\n\nGet list of result types \
for finished job.\n\n  my (@result_types) = &rest_\
get_result_types($job_id);\n\n=cut\n\nsub rest_get\
_result_types {\n	print_debug_message( 'rest_get_r\
esult_types', 'Begin', 1 );\n	my $job_id = shift;\\
n	print_debug_message( 'rest_get_result_types', 'j\
obid: ' . $job_id, 2 );\n	my (@resultTypes);\n	my \
$url                      = $baseUrl . '/resulttyp\
es/' . $job_id;\n	my $result_type_list_xml_str = &\
rest_request($url);\n	my $result_type_list_xml    \
 = XMLin($result_type_list_xml_str);\n	(@resultTyp\
es) = @{ $result_type_list_xml->{'type'} };\n	prin\
t_debug_message( 'rest_get_result_types',\n		scala\
r(@resultTypes) . ' result types', 2 );\n	print_de\
bug_message( 'rest_get_result_types', 'End', 1 );\\
n	return (@resultTypes);\n}\n\n=head2 rest_get_res\
ult()\n\nGet result data of a specified type for a\
 finished job.\n\n  my $result = rest_get_result($\
job_id, $result_type);\n\n=cut\n\nsub rest_get_res\
ult {\n	print_debug_message( 'rest_get_result', 'B\
egin', 1 );\n	my $job_id = shift;\n	my $type   = s\
hift;\n	print_debug_message( 'rest_get_result', 'j\
obid: ' . $job_id, 1 );\n	print_debug_message( 're\
st_get_result', 'type: ' . $type,    1 );\n	my $ur\
l    = $baseUrl . '/result/' . $job_id . '/' . $ty\
pe;\n	my $result = &rest_request($url);\n	print_de\
bug_message( 'rest_get_result', length($result) . \
' characters',\n		1 );\n	print_debug_message( 'res\
t_get_result', 'End', 1 );\n	return $result;\n}\n\\
n\n=head2 print_debug_message()\n\nPrint debug mes\
sage at specified debug level.\n\n  &print_debug_m\
essage($method_name, $message, $level);\n\n=cut\n\\
nsub print_debug_message {\n	my $function_name = s\
hift;\n	my $message       = shift;\n	my $level    \
     = shift;\n	if ( $level <= $params{'debugLevel\
'} ) {\n		print STDERR '[', $function_name, '()] '\
, $message, \"\\n\";\n	}\n}\n\n=head2 print_tool_p\
arams()\n\nPrint list of tool parameters.\n\n  &pr\
int_tool_params();\n\n=cut\n\nsub print_tool_param\
s {\n	print_debug_message( 'print_tool_params', 'B\
egin', 1 );\n	my (@param_list) = &rest_get_paramet\
ers();\n	foreach my $param ( sort(@param_list) ) {\
\n		print $param, \"\\n\";\n	}\n	print_debug_messa\
ge( 'print_tool_params', 'End', 1 );\n}\n\n=head2 \
print_param_details()\n\nPrint details of a tool p\
arameter.\n\n  &print_param_details($param_name);\\
n\n=cut\n\nsub print_param_details {\n	print_debug\
_message( 'print_param_details', 'Begin', 1 );\n	m\
y $paramName = shift;\n	print_debug_message( 'prin\
t_param_details', 'paramName: ' . $paramName, 2 );\
\n	my $paramDetail = &rest_get_parameter_details($\
paramName);\n	print $paramDetail->{'name'}, \"\\t\\
", $paramDetail->{'type'}, \"\\n\";\n	print $param\
Detail->{'description'}, \"\\n\";\n	foreach my $va\
lue ( @{ $paramDetail->{'values'}->{'value'} } ) {\
\n		print $value->{'value'};\n		if ( $value->{'def\
aultValue'} eq 'true' ) {\n			print \"\\t\", 'defa\
ult';\n		}\n		print \"\\n\";\n		print \"\\t\", $va\
lue->{'label'}, \"\\n\";\n	}\n	print_debug_message\
( 'print_param_details', 'End', 1 );\n}\n\n=head2 \
print_job_status()\n\nPrint status of a job.\n\n  \
&print_job_status($job_id);\n\n=cut\n\nsub print_j\
ob_status {\n	print_debug_message( 'print_job_stat\
us', 'Begin', 1 );\n	my $jobid = shift;\n	print_de\
bug_message( 'print_job_status', 'jobid: ' . $jobi\
d, 1 );\n	if ( $outputLevel > 0 ) {\n		print STDER\
R 'Getting status for job ', $jobid, \"\\n\";\n	}\\
n	my $result = &rest_get_status($jobid);\n	print \\
"$result\\n\";\n	if ( $result eq 'FINISHED' && $ou\
tputLevel > 0 ) {\n		print STDERR \"To get results\
: $scriptName --polljob --jobid \" . $jobid\n		  .\
 \"\\n\";\n	}\n	print_debug_message( 'print_job_st\
atus', 'End', 1 );\n}\n\n=head2 print_result_types\
()\n\nPrint available result types for a job.\n\n \
 &print_result_types($job_id);\n\n=cut\n\nsub prin\
t_result_types {\n	print_debug_message( 'result_ty\
pes', 'Begin', 1 );\n	my $jobid = shift;\n	print_d\
ebug_message( 'result_types', 'jobid: ' . $jobid, \
1 );\n	if ( $outputLevel > 0 ) {\n		print STDERR '\
Getting result types for job ', $jobid, \"\\n\";\n\
	}\n	my $status = &rest_get_status($jobid);\n	if (\
 $status eq 'PENDING' || $status eq 'RUNNING' ) {\\
n		print STDERR 'Error: Job status is ', $status,\\
n		  '. To get result types the job must be finish\
ed.', \"\\n\";\n	}\n	else {\n		my (@resultTypes) =\
 &rest_get_result_types($jobid);\n		if ( $outputLe\
vel > 0 ) {\n			print STDOUT 'Available result typ\
es:', \"\\n\";\n		}\n		foreach my $resultType (@re\
sultTypes) {\n			print STDOUT $resultType->{'ident\
ifier'}, \"\\n\";\n			if ( defined( $resultType->{\
'label'} ) ) {\n				print STDOUT \"\\t\", $resultT\
ype->{'label'}, \"\\n\";\n			}\n			if ( defined( $\
resultType->{'description'} ) ) {\n				print STDOU\
T \"\\t\", $resultType->{'description'}, \"\\n\";\\
n			}\n			if ( defined( $resultType->{'mediaType'}\
 ) ) {\n				print STDOUT \"\\t\", $resultType->{'m\
ediaType'}, \"\\n\";\n			}\n			if ( defined( $resu\
ltType->{'fileSuffix'} ) ) {\n				print STDOUT \"\\
\t\", $resultType->{'fileSuffix'}, \"\\n\";\n			}\\
n		}\n		if ( $status eq 'FINISHED' && $outputLevel\
 > 0 ) {\n			print STDERR \"\\n\", 'To get results\
:', \"\\n\",\n			  \"  $scriptName --polljob --job\
id \" . $params{'jobid'} . \"\\n\",\n			  \"  $scr\
iptName --polljob --outformat <type> --jobid \"\n	\
		  . $params{'jobid'} . \"\\n\";\n		}\n	}\n	print\
_debug_message( 'result_types', 'End', 1 );\n}\n\n\
=head2 submit_job()\n\nSubmit a job to the service\
.\n\n  &submit_job($seq);\n\n=cut\n\nsub submit_jo\
b {\n	print_debug_message( 'submit_job', 'Begin', \
1 );\n\n	# Set input sequence\n	$tool_params{'sequ\
ence'} = shift;\n\n	# Load parameters\n	&load_para\
ms();\n\n	# Submit the job\n	my $jobid = &rest_run\
( $params{'email'}, $params{'title'}, \\%tool_para\
ms );\n\n	# Simulate sync/async mode\n	if ( define\
d( $params{'async'} ) ) {\n		print STDOUT $jobid, \
\"\\n\";\n		if ( $outputLevel > 0 ) {\n			print ST\
DERR\n			  \"To check status: $scriptName --status\
 --jobid $jobid\\n\";\n		}\n	}\n	else {\n		if ( $o\
utputLevel > 0 ) {\n			print STDERR \"JobId: $jobi\
d\\n\";\n		}\n		sleep 1;\n		&get_results($jobid);\\
n	}\n	print_debug_message( 'submit_job', 'End', 1 \
);\n}\n\n=head2 multi_submit_job()\n\nSubmit multi\
ple jobs assuming input is a collection of fasta f\
ormatted sequences.\n\n  &multi_submit_job();\n\n=\
cut\n\nsub multi_submit_job {\n	print_debug_messag\
e( 'multi_submit_job', 'Begin', 1 );\n	my $jobIdFo\
rFilename = 1;\n	$jobIdForFilename = 0 if ( define\
d( $params{'outfile'} ) );\n	my (@filename_list) =\
 ();\n\n	# Query sequence\n	if ( defined( $ARGV[0]\
 ) ) {    # Bare option\n		if ( -f $ARGV[0] || $AR\
GV[0] eq '-' ) {    # File\n			push( @filename_lis\
t, $ARGV[0] );\n		}\n	}\n	if ( $params{'sequence'}\
 ) {                   # Via --sequence\n		if ( -f\
 $params{'sequence'} || $params{'sequence'} eq '-'\
 ) {    # File\n			push( @filename_list, $params{'\
sequence'} );\n		}\n	}\n\n	$/ = '>';\n	foreach my \
$filename (@filename_list) {\n		open( my $INFILE, \
'<', $filename )\n		  or die \"Error: unable to op\
en file $filename ($!)\";\n		while (<$INFILE>) {\n\
			my $seq = $_;\n			$seq =~ s/>$//;\n			if ( $seq\
 =~ m/(\\S+)/ ) {\n				print STDERR \"Submitting j\
ob for: $1\\n\"\n				  if ( $outputLevel > 0 );\n	\
			$seq = '>' . $seq;\n				&print_debug_message( '\
multi_submit_job', $seq, 11 );\n				&submit_job($s\
eq);\n				$params{'outfile'} = undef if ( $jobIdFo\
rFilename == 1 );\n			}\n		}\n		close $INFILE;\n	}\
\n	print_debug_message( 'multi_submit_job', 'End',\
 1 );\n}\n\n=head2 list_file_submit_job()\n\nSubmi\
t multiple jobs using a file containing a list of \
entry identifiers as \ninput.\n\n  &list_file_subm\
it_job($list_filename)\n\n=cut\n\nsub list_file_su\
bmit_job {\n	my $filename         = shift;\n	my $j\
obIdForFilename = 1;\n	$jobIdForFilename = 0 if ( \
defined( $params{'outfile'} ) );\n\n	# Iterate ove\
r identifiers, submitting each job\n	open( my $LIS\
TFILE, '<', $filename )\n	  or die 'Error: unable \
to open file ' . $filename . ' (' . $! . ')';\n	wh\
ile (<$LISTFILE>) {\n		my $line = $_;\n		chomp($li\
ne);\n		if ( $line ne '' ) {\n			&print_debug_mess\
age( 'list_file_submit_job', 'line: ' . $line, 2 )\
;\n			if ( $line =~ m/\\w:\\w/ ) {    # Check this\
 is an identifier\n				print STDERR \"Submitting j\
ob for: $line\\n\"\n				  if ( $outputLevel > 0 );\
\n				&submit_job($line);\n			}\n			else {\n				pr\
int STDERR\n\"Warning: line \\\"$line\\\" is not r\
ecognised as an identifier\\n\";\n			}\n		}\n		$pa\
rams{'outfile'} = undef if ( $jobIdForFilename == \
1 );\n	}\n	close $LISTFILE;\n}\n\n=head2 load_data\
()\n\nLoad sequence data from file or option speci\
fied on the command-line.\n\n  &load_data();\n\n=c\
ut\n\nsub load_data {\n	print_debug_message( 'load\
_data', 'Begin', 1 );\n	my $retSeq;\n\n	# Query se\
quence\n	if ( defined( $ARGV[0] ) ) {    # Bare op\
tion\n		if ( -f $ARGV[0] || $ARGV[0] eq '-' ) {   \
 # File\n			$retSeq = &read_file( $ARGV[0] );\n		}\
\n		else {                                     # D\
B:ID or sequence\n			$retSeq = $ARGV[0];\n		}\n	}\\
n	if ( $params{'sequence'} ) {                   #\
 Via --sequence\n		if ( -f $params{'sequence'} || \
$params{'sequence'} eq '-' ) {    # File\n			$retS\
eq = &read_file( $params{'sequence'} );\n		}\n		el\
se {    # DB:ID or sequence\n			$retSeq = $params{\
'sequence'};\n		}\n	}\n	print_debug_message( 'load\
_data', 'End', 1 );\n	return $retSeq;\n}\n\n=head2\
 load_params()\n\nLoad job parameters from command\
-line options.\n\n  &load_params();\n\n=cut\n\nsub\
 load_params {\n	print_debug_message( 'load_params\
', 'Begin', 1 );\n\n	# Database(s) to search\n	my \
(@dbList) = split /[ ,]/, $params{'database'};\n	$\
tool_params{'database'} = \\@dbList;\n\n	# Match/m\
issmatch\n	if ( $params{'match'} && $params{'missm\
atch'} ) {\n		$tool_params{'match_scores'} =\n		  \
$params{'match'} . ',' . $params{'missmatch'};\n	}\
\n	\n	# Compatability options, old command-line\n	\
if(!$tool_params{'alignments'} && $params{'numal'}\
) {\n		$tool_params{'alignments'} = $params{'numal\
'};\n	}\n	if(!$tool_params{'gapopen'} && $params{'\
opengap'}) {\n		$tool_params{'gapopen'} = $params{\
'opengap'};\n	}\n	if(!$tool_params{'gapext'} && $p\
arams{'extendgap'}) {\n		$tool_params{'gapext'} = \
$params{'extendgap'};\n	}\n\n	print_debug_message(\
 'load_params', 'End', 1 );\n}\n\n=head2 client_po\
ll()\n\nClient-side job polling.\n\n  &client_poll\
($job_id);\n\n=cut\n\nsub client_poll {\n	print_de\
bug_message( 'client_poll', 'Begin', 1 );\n	my $jo\
bid  = shift;\n	my $status = 'PENDING';\n\n	my $er\
rorCount = 0;\n	while ($status eq 'RUNNING'\n		|| \
$status eq 'PENDING'\n		|| ( $status eq 'ERROR' &&\
 $errorCount < 2 ) )\n	{\n		$status = rest_get_sta\
tus($jobid);\n		print STDERR \"$status\\n\" if ( $\
outputLevel > 0 );\n		if ( $status eq 'ERROR' ) {\\
n			$errorCount++;\n		}\n		elsif ( $errorCount > 0\
 ) {\n			$errorCount--;\n		}\n		if (   $status eq \
'RUNNING'\n			|| $status eq 'PENDING'\n			|| $stat\
us eq 'ERROR' )\n		{\n\n			# Wait before polling a\
gain.\n			sleep $checkInterval;\n		}\n	}\n	print_d\
ebug_message( 'client_poll', 'End', 1 );\n	return \
$status;\n}\n\n=head2 get_results()\n\nGet the res\
ults for a job identifier.\n\n  &get_results($job_\
id);\n\n=cut\n\nsub get_results {\n	print_debug_me\
ssage( 'get_results', 'Begin', 1 );\n	my $jobid = \
shift;\n	print_debug_message( 'get_results', 'jobi\
d: ' . $jobid, 1 );\n\n	# Verbose\n	if ( $outputLe\
vel > 1 ) {\n		print 'Getting results for job ', $\
jobid, \"\\n\";\n	}\n\n	# Check status, and wait i\
f not finished\n	client_poll($jobid);\n\n	# Use Jo\
bId if output file name is not defined\n	unless ( \
defined( $params{'outfile'} ) ) {\n		$params{'outf\
ile'} = $jobid;\n	}\n\n	# Get list of data types\n\
	my (@resultTypes) = rest_get_result_types($jobid)\
;\n\n	# Get the data and write it to a file\n	if (\
 defined( $params{'outformat'} ) ) {    # Specifie\
d data type\n		my $selResultType;\n		foreach my $r\
esultType (@resultTypes) {\n			if ( $resultType->{\
'identifier'} eq $params{'outformat'} ) {\n				$se\
lResultType = $resultType;\n			}\n		}\n		if ( defi\
ned($selResultType) ) {\n			my $result =\n			  res\
t_get_result( $jobid, $selResultType->{'identifier\
'} );\n			if ( $params{'outfile'} eq '-' ) {\n				\
write_file( $params{'outfile'}, $result );\n			}\n\
			else {\n				write_file(\n					$params{'outfile'\
} . '.'\n					  . $selResultType->{'identifier'} .\
 '.'\n					  . $selResultType->{'fileSuffix'},\n		\
			$result\n				);\n			}\n		}\n		else {\n			die 'E\
rror: unknown result format \"' . $params{'outform\
at'} . '\"';\n		}\n	}\n	else {    # Data types ava\
ilable\n		      # Write a file for each output typ\
e\n		for my $resultType (@resultTypes) {\n			if ( \
$outputLevel > 1 ) {\n				print STDERR 'Getting ',\
 $resultType->{'identifier'}, \"\\n\";\n			}\n			m\
y $result = rest_get_result( $jobid, $resultType->\
{'identifier'} );\n			if ( $params{'outfile'} eq '\
-' ) {\n				write_file( $params{'outfile'}, $resul\
t );\n			}\n			else {\n				write_file(\n					$para\
ms{'outfile'} . '.'\n					  . $resultType->{'ident\
ifier'} . '.'\n					  . $resultType->{'fileSuffix'\
},\n					$result\n				);\n			}\n		}\n	}\n	print_de\
bug_message( 'get_results', 'End', 1 );\n}\n\n=hea\
d2 read_file()\n\nRead a file into a scalar. The s\
pecial filename '-' can be used to read from \nsta\
ndard input (STDIN).\n\n  my $data = &read_file($f\
ilename);\n\n=cut\n\nsub read_file {\n	print_debug\
_message( 'read_file', 'Begin', 1 );\n	my $filenam\
e = shift;\n	print_debug_message( 'read_file', 'fi\
lename: ' . $filename, 2 );\n	my ( $content, $buff\
er );\n	if ( $filename eq '-' ) {\n		while ( sysre\
ad( STDIN, $buffer, 1024 ) ) {\n			$content .= $bu\
ffer;\n		}\n	}\n	else {    # File\n		open( my $FIL\
E, '<', $filename )\n		  or die \"Error: unable to\
 open input file $filename ($!)\";\n		while ( sysr\
ead( $FILE, $buffer, 1024 ) ) {\n			$content .= $b\
uffer;\n		}\n		close($FILE);\n	}\n	print_debug_mes\
sage( 'read_file', 'End', 1 );\n	return $content;\\
n}\n\n=head2 write_file()\n\nWrite data to a file.\
 The special filename '-' can be used to write to \
\nstandard output (STDOUT).\n\n  &write_file($file\
name, $data);\n\n=cut\n\nsub write_file {\n	print_\
debug_message( 'write_file', 'Begin', 1 );\n	my ( \
$filename, $data ) = @_;\n	print_debug_message( 'w\
rite_file', 'filename: ' . $filename, 2 );\n	if ( \
$outputLevel > 0 ) {\n		print STDERR 'Creating res\
ult file: ' . $filename . \"\\n\";\n	}\n	if ( $fil\
ename eq '-' ) {\n		print STDOUT $data;\n	}\n	else\
 {\n		open( my $FILE, '>', $filename )\n		  or die\
 \"Error: unable to open output file $filename ($!\
)\";\n		syswrite( $FILE, $data );\n		close($FILE);\
\n	}\n	print_debug_message( 'write_file', 'End', 1\
 );\n}\n\n=head2 usage()\n\nPrint program usage me\
ssage.\n\n  &usage();\n\n=cut\n\nsub usage {\n	pri\
nt STDERR <<EOF\nNCBI BLAST\n==========\n   \nRapi\
d sequence database search programs utilizing the \
BLAST algorithm\n    \n[Required]\n\n  -p, --progr\
am      : str  : BLAST program to use, see --param\
Detail program\n  -D, --database     : str  : data\
base(s) to search, space separated. See\n         \
                     --paramDetail database\n     \
 --stype        : str  : query sequence type, see \
--paramDetail stype\n  seqFile            : file :\
 query sequence (\"-\" for STDIN, \\@filename for\\
n                              identifier list fil\
e)\n\n[Optional]\n\n  -m, --matrix       : str  : \
scoring matrix, see --paramDetail matrix\n  -e, --\
exp          : real : 0<E<= 1000. Statistical sign\
ificance threshold \n                             \
 for reporting database sequence matches.\n  -f, -\
-filter       :      : filter the query sequence f\
or low complexity \n                              \
regions, see --paramDetail filter\n  -A, --align  \
      : int  : pairwise alignment format, see --pa\
ramDetail align\n  -s, --scores       : int  : num\
ber of scores to be reported\n  -n, --alignments  \
 : int  : number of alignments to report\n  -u, --\
match        : int  : Match score (BLASTN only)\n \
 -v, --mismatch     : int  : Mismatch score (BLAST\
N only)\n  -o, --gapopen      : int  : Gap open pe\
nalty\n  -x, --gapext       : int  : Gap extension\
 penalty\n  -d, --dropoff      : int  : Drop-off\n\
  -g, --gapalign     :      : Optimise gapped alig\
nments\n      --seqrange     : str  : region withi\
n input to use as query\n      --multifasta   :   \
   : treat input as a set of fasta formatted seque\
nces\n\n[General]\n\n  -h, --help        :      : \
prints this help text\n      --async       :      \
: forces to make an asynchronous query\n      --em\
ail       : str  : e-mail address\n      --title  \
     : str  : title for job\n      --status      :\
      : get job status\n      --resultTypes :     \
 : get available result types for job\n      --pol\
ljob     :      : poll for the status of a job\n  \
    --jobid       : str  : jobid that was returned\
 when an asynchronous job \n                      \
       was submitted.\n      --outfile     : str  \
: file name for results (default is jobid;\n      \
                       \"-\" for STDOUT)\n      --\
outformat   : str  : result format to retrieve\n  \
    --params      :      : list input parameters\n\
      --paramDetail : str  : display details for i\
nput parameter\n      --quiet       :      : decre\
ase output\n      --verbose     :      : increase \
output\n      --trace       :      : show SOAP mes\
sages being interchanged \n   \nSynchronous job:\n\
\n  The results/errors are returned as soon as the\
 job is finished.\n  Usage: $scriptName --email <y\
our\\@email> [options...] seqFile\n  Returns: resu\
lts as an attachment\n\nAsynchronous job:\n\n  Use\
 this if you want to retrieve the results at a lat\
er time. The results \n  are stored for up to 24 h\
ours. 	\n  Usage: $scriptName --async --email <you\
r\\@email> [options...] seqFile\n  Returns: jobid\\
n\n  Use the jobid to query for the status of the \
job. If the job is finished, \n  it also returns t\
he results/errors.\n  Usage: $scriptName --polljob\
 --jobid <jobId> [--outfile string]\n  Returns: st\
ring indicating the status of the job and if appli\
cable, results \n  as an attachment.\n\nFurther in\
formation:\n\n  http://www.ebi.ac.uk/Tools/webserv\
ices/services/sss/ncbi_blast_rest\n  http://www.eb\
i.ac.uk/Tools/webservices/tutorials/perl\n\nSuppor\
t/Feedback:\n\n  http://www.ebi.ac.uk/support/\nEO\
F\n}\n\n=head1 FEEDBACK/SUPPORT\n\nPlease contact \
us at L<http://www.ebi.ac.uk/support/> if you have\
 any \nfeedback, suggestions or issues with the se\
rvice or this client.\n\n=cut\n","\n=head1 NAME\n\\
nwublast_lwp.pl\n\n=head1 DESCRIPTION\n\nWU-BLAST \
REST web service Perl client using L<LWP>.\n\nTest\
ed with:\n\n=over\n\n=item *\nL<LWP> 5.79, L<XML::\
Simple> 2.12 and Perl 5.8.3\n\n=item *\nL<LWP> 5.8\
05, L<XML::Simple> 2.14 and Perl 5.8.7\n\n=item *\\
nL<LWP> 5.820, L<XML::Simple> 2.18 and Perl 5.10.0\
 (Ubuntu 9.04)\n\n=back\n\nFor further information\
 see:\n\n=over\n\n=item *\nL<http://www.ebi.ac.uk/\
Tools/webservices/services/sss/wu_blast_rest>\n\n=\
item *\nL<http://www.ebi.ac.uk/Tools/webservices/t\
utorials/perl>\n\n=back\n\n=head1 VERSION\n\n$Id: \
wublast_lwp.pl 1317 2009-09-03 15:44:11Z hpm $\n\n\
=cut\n\nuse strict;\nuse warnings;\n\nuse English;\
\nuse LWP;\nuse XML::Simple;\nuse Getopt::Long qw(\
:config no_ignore_case bundling);\nuse File::Basen\
ame;\nuse Data::Dumper;\n\nmy $baseUrl = 'http://w\
ww.ebi.ac.uk/Tools/services/rest/wublast';\n\nmy $\
checkInterval = 3;\n\nmy $outputLevel = 1;\n\nmy $\
numOpts = scalar(@ARGV);\nmy %params = ( 'debugLev\
el' => 0 );\n\nmy %tool_params = ();\nGetOptions(\\
n\n	# Tool specific options\n	'program|p=s'     =>\
 \\$tool_params{'program'},      # BLAST program\n\
	'database|D=s'    => \\$params{'database'},     #\
 Search database\n	'matrix|m=s'      => \\$tool_pa\
rams{'matrix'},       # Scoring matrix\n	'exp|E=f'\
         => \\$tool_params{'exp'},          # E-va\
lue threshold\n	'viewfilter|e'    => \\$tool_param\
s{'viewfilter'},   # Display filtered sequence\n	'\
filter|f=s'      => \\$tool_params{'filter'},     \
  # Low complexity filter name\n	'alignments|n=i' \
 => \\$tool_params{'alignments'},   # Number of al\
ignments\n	'scores|s=i'      => \\$tool_params{'sc\
ores'},       # Number of scores\n	'sensitivity|S=\
s' => \\$tool_params{'sensitivity'},  # Search sen\
sitivity\n	'sort|t=s'        => \\$tool_params{'so\
rt'},         # Sort hits by...\n	'stats|T=s'     \
  => \\$tool_params{'stats'},        # Scoring sta\
tistic to use\n	'strand|d=s'      => \\$tool_param\
s{'strand'},       # Strand to use\n	'topcombon|c=\
i'   => \\$tool_params{'topcombon'},    # Consiste\
nt sets of HSPs\n	'align|A=i'       => \\$tool_par\
ams{'align'},   # Pairwise alignment format\n	'sty\
pe=s' => \\$tool_params{'stype'},    # Sequence ty\
pe 'protein' or 'dna'\n	'sequence=s' => \\$params{\
'sequence'},         # Query sequence file or DB:I\
D\n	'multifasta' => \\$params{'multifasta'},      \
 # Multiple fasta input\n\n	# Compatability option\
s, old command-line.\n	'echofilter|e'    => \\$par\
ams{'echofilter'},   # Display filtered sequence\n\
	'b=i'  => \\$params{'numal'},        # Number of \
alignments\n	'appxml=s'        => \\$params{'appxm\
l'},       # Application XML\n\n	# Generic options\
\n	'email=s'       => \\$params{'email'},         \
 # User e-mail address\n	'title=s'       => \\$par\
ams{'title'},          # Job title\n	'outfile=s'  \
   => \\$params{'outfile'},        # Output file n\
ame\n	'outformat=s'   => \\$params{'outformat'},  \
    # Output file type\n	'jobid=s'       => \\$par\
ams{'jobid'},          # JobId\n	'help|h'        =\
> \\$params{'help'},           # Usage help\n	'asy\
nc'         => \\$params{'async'},          # Asyn\
chronous submission\n	'polljob'       => \\$params\
{'polljob'},        # Get results\n	'resultTypes' \
  => \\$params{'resultTypes'},    # Get result typ\
es\n	'status'        => \\$params{'status'},      \
   # Get status\n	'params'        => \\$params{'pa\
rams'},         # List input parameters\n	'paramDe\
tail=s' => \\$params{'paramDetail'},    # Get deta\
ils for parameter\n	'quiet'         => \\$params{'\
quiet'},          # Decrease output level\n	'verbo\
se'       => \\$params{'verbose'},        # Increa\
se output level\n	'debugLevel=i'  => \\$params{'de\
bugLevel'},     # Debug output level\n	'baseUrl=s'\
     => \\$baseUrl,                  # Base URL fo\
r service.\n);\nif ( $params{'verbose'} ) { $outpu\
tLevel++ }\nif ( $params{'$quiet'} )  { $outputLev\
el-- }\n\n&print_debug_message( 'MAIN', 'LWP::VERS\
ION: ' . $LWP::VERSION,\n	1 );\n\n&print_debug_mes\
sage( 'MAIN', \"params:\\n\" . Dumper( \\%params )\
,           11 );\n&print_debug_message( 'MAIN', \\
"tool_params:\\n\" . Dumper( \\%tool_params ), 11 \
);\n\nmy $scriptName = basename( $0, () );\n\nif (\
 $params{'help'} || $numOpts == 0 ) {\n	&usage();\\
n	exit(0);\n}\n\n&print_debug_message( 'MAIN', 'ba\
seUrl: ' . $baseUrl, 1 );\n\nif (\n	!(\n		   $para\
ms{'polljob'}\n		|| $params{'resultTypes'}\n		|| $\
params{'status'}\n		|| $params{'params'}\n		|| $pa\
rams{'paramDetail'}\n	)\n	&& !( defined( $ARGV[0] \
) || defined( $params{'sequence'} ) )\n  )\n{\n\n	\
# Bad argument combination, so print error message\
 and usage\n	print STDERR 'Error: bad option combi\
nation', \"\\n\";\n	&usage();\n	exit(1);\n}\n\nels\
if ( $params{'params'} ) {\n	&print_tool_params();\
\n}\n\nelsif ( $params{'paramDetail'} ) {\n	&print\
_param_details( $params{'paramDetail'} );\n}\n\nel\
sif ( $params{'status'} && defined( $params{'jobid\
'} ) ) {\n	&print_job_status( $params{'jobid'} );\\
n}\n\nelsif ( $params{'resultTypes'} && defined( $\
params{'jobid'} ) ) {\n	&print_result_types( $para\
ms{'jobid'} );\n}\n\nelsif ( $params{'polljob'} &&\
 defined( $params{'jobid'} ) ) {\n	&get_results( $\
params{'jobid'} );\n}\n\nelse {\n\n	# Multiple inp\
ut sequence mode, assume fasta format.\n	if ( $par\
ams{'multifasta'} ) {\n		&multi_submit_job();\n	}\\
n\n	# Entry identifier list file.\n	elsif (( defin\
ed( $params{'sequence'} ) && $params{'sequence'} =\
~ m/^\\@/ )\n		|| ( defined( $ARGV[0] ) && $ARGV[0\
] =~ m/^\\@/ ) )\n	{\n		my $list_filename = $param\
s{'sequence'} || $ARGV[0];\n		$list_filename =~ s/\
^\\@//;\n		&list_file_submit_job($list_filename);\\
n	}\n\n	# Default: single sequence/identifier.\n	e\
lse {\n\n		# Load the sequence data and submit.\n	\
	&submit_job( &load_data() );\n	}\n}\n\n=head1 FUN\
CTIONS\n\n=cut\n\n\n=head2 rest_request()\n\nPerfo\
rm a REST request.\n\n  my $response_str = &rest_r\
equest($url);\n\n=cut\n\nsub rest_request {\n	prin\
t_debug_message( 'rest_request', 'Begin', 11 );\n	\
my $requestUrl = shift;\n	print_debug_message( 're\
st_request', 'URL: ' . $requestUrl, 11 );\n\n	# Cr\
eate a user agent\n	my $ua = LWP::UserAgent->new()\
;\n	'$Revision: 1317 $' =~ m/(\\d+)/;\n	$ua->agent\
(\"EBI-Sample-Client/$1 ($scriptName; $OSNAME) \" \
. $ua->agent());\n	$ua->env_proxy;\n\n	# Perform t\
he request\n	my $response = $ua->get($requestUrl);\
\n	print_debug_message( 'rest_request', 'HTTP stat\
us: ' . $response->code,\n		11 );\n\n	# Check for \
HTTP error codes\n	if ( $response->is_error ) {\n	\
	$response->content() =~ m/<h1>([^<]+)<\\/h1>/;\n	\
	die 'http status: ' . $response->code . ' ' . $re\
sponse->message . '  ' . $1;\n	}\n	print_debug_mes\
sage( 'rest_request', 'End', 11 );\n\n	# Return th\
e response data\n	return $response->content();\n}\\
n\n=head2 rest_get_parameters()\n\nGet list of too\
l parameter names.\n\n  my (@param_list) = &rest_g\
et_parameters();\n\n=cut\n\nsub rest_get_parameter\
s {\n	print_debug_message( 'rest_get_parameters', \
'Begin', 1 );\n	my $url                = $baseUrl \
. '/parameters/';\n	my $param_list_xml_str = rest_\
request($url);\n	my $param_list_xml     = XMLin($p\
aram_list_xml_str);\n	my (@param_list)       = @{ \
$param_list_xml->{'id'} };\n	print_debug_message( \
'rest_get_parameters', 'End', 1 );\n	return (@para\
m_list);\n}\n\n=head2 rest_get_parameter_details()\
\n\nGet details of a tool parameter.\n\n  my $para\
mDetail = &rest_get_parameter_details($param_name)\
;\n\n=cut\n\nsub rest_get_parameter_details {\n	pr\
int_debug_message( 'rest_get_parameter_details', '\
Begin', 1 );\n	my $parameterId = shift;\n	print_de\
bug_message( 'rest_get_parameter_details',\n		'par\
ameterId: ' . $parameterId, 1 );\n	my $url        \
          = $baseUrl . '/parameterdetails/' . $par\
ameterId;\n	my $param_detail_xml_str = rest_reques\
t($url);\n	my $param_detail_xml     = XMLin($param\
_detail_xml_str);\n	print_debug_message( 'rest_get\
_parameter_details', 'End', 1 );\n	return ($param_\
detail_xml);\n}\n\n=head2 rest_run()\n\nSubmit a j\
ob.\n\n  my $job_id = &rest_run($email, $title, \\\
%params );\n\n=cut\n\nsub rest_run {\n	print_debug\
_message( 'rest_run', 'Begin', 1 );\n	my $email  =\
 shift;\n	my $title  = shift;\n	my $params = shift\
;\n	print_debug_message( 'rest_run', 'email: ' . $\
email, 1 );\n	if ( defined($title) ) {\n		print_de\
bug_message( 'rest_run', 'title: ' . $title, 1 );\\
n	}\n	print_debug_message( 'rest_run', 'params: ' \
. Dumper($params), 1 );\n\n	# User agent to perfor\
m http requests\n	my $ua = LWP::UserAgent->new();\\
n	$ua->env_proxy;\n\n	# Clean up parameters\n	my (\
%tmp_params) = %{$params};\n	$tmp_params{'email'} \
= $email;\n	$tmp_params{'title'} = $title;\n	forea\
ch my $param_name ( keys(%tmp_params) ) {\n		if ( \
!defined( $tmp_params{$param_name} ) ) {\n			delet\
e $tmp_params{$param_name};\n		}\n	}\n\n	# Submit \
the job as a POST\n	my $url = $baseUrl . '/run';\n\
	my $response = $ua->post( $url, \\%tmp_params );\\
n	print_debug_message( 'rest_run', 'HTTP status: '\
 . $response->code, 11 );\n	print_debug_message( '\
rest_run',\n		'request: ' . $response->request()->\
content(), 11 );\n\n	# Check for HTTP error codes\\
n	if ( $response->is_error ) {\n		$response->conte\
nt() =~ m/<h1>([^<]+)<\\/h1>/;\n		die 'http status\
: ' . $response->code . ' ' . $response->message .\
 '  ' . $1;\n	}\n\n	# The job id is returned\n	my \
$job_id = $response->content();\n	print_debug_mess\
age( 'rest_run', 'End', 1 );\n	return $job_id;\n}\\
n\n=head2 rest_get_status()\n\nCheck the status of\
 a job.\n\n  my $status = &rest_get_status($job_id\
);\n\n=cut\n\nsub rest_get_status {\n	print_debug_\
message( 'rest_get_status', 'Begin', 1 );\n	my $jo\
b_id = shift;\n	print_debug_message( 'rest_get_sta\
tus', 'jobid: ' . $job_id, 2 );\n	my $status_str =\
 'UNKNOWN';\n	my $url        = $baseUrl . '/status\
/' . $job_id;\n	$status_str = &rest_request($url);\
\n	print_debug_message( 'rest_get_status', 'status\
_str: ' . $status_str, 2 );\n	print_debug_message(\
 'rest_get_status', 'End', 1 );\n	return $status_s\
tr;\n}\n\n=head2 rest_get_result_types()\n\nGet li\
st of result types for finished job.\n\n  my (@res\
ult_types) = &rest_get_result_types($job_id);\n\n=\
cut\n\nsub rest_get_result_types {\n	print_debug_m\
essage( 'rest_get_result_types', 'Begin', 1 );\n	m\
y $job_id = shift;\n	print_debug_message( 'rest_ge\
t_result_types', 'jobid: ' . $job_id, 2 );\n	my (@\
resultTypes);\n	my $url                      = $ba\
seUrl . '/resulttypes/' . $job_id;\n	my $result_ty\
pe_list_xml_str = &rest_request($url);\n	my $resul\
t_type_list_xml     = XMLin($result_type_list_xml_\
str);\n	(@resultTypes) = @{ $result_type_list_xml-\
>{'type'} };\n	print_debug_message( 'rest_get_resu\
lt_types',\n		scalar(@resultTypes) . ' result type\
s', 2 );\n	print_debug_message( 'rest_get_result_t\
ypes', 'End', 1 );\n	return (@resultTypes);\n}\n\n\
=head2 rest_get_result()\n\nGet result data of a s\
pecified type for a finished job.\n\n  my $result \
= rest_get_result($job_id, $result_type);\n\n=cut\\
n\nsub rest_get_result {\n	print_debug_message( 'r\
est_get_result', 'Begin', 1 );\n	my $job_id = shif\
t;\n	my $type   = shift;\n	print_debug_message( 'r\
est_get_result', 'jobid: ' . $job_id, 1 );\n	print\
_debug_message( 'rest_get_result', 'type: ' . $typ\
e,    1 );\n	my $url    = $baseUrl . '/result/' . \
$job_id . '/' . $type;\n	my $result = &rest_reques\
t($url);\n	print_debug_message( 'rest_get_result',\
 length($result) . ' characters',\n		1 );\n	print_\
debug_message( 'rest_get_result', 'End', 1 );\n	re\
turn $result;\n}\n\n\n=head2 print_debug_message()\
\n\nPrint debug message at specified debug level.\\
n\n  &print_debug_message($method_name, $message, \
$level);\n\n=cut\n\nsub print_debug_message {\n	my\
 $function_name = shift;\n	my $message       = shi\
ft;\n	my $level         = shift;\n	if ( $level <= \
$params{'debugLevel'} ) {\n		print STDERR '[', $fu\
nction_name, '()] ', $message, \"\\n\";\n	}\n}\n\n\
=head2 print_tool_params()\n\nPrint list of tool p\
arameters.\n\n  &print_tool_params();\n\n=cut\n\ns\
ub print_tool_params {\n	print_debug_message( 'pri\
nt_tool_params', 'Begin', 1 );\n	my (@param_list) \
= &rest_get_parameters();\n	foreach my $param ( so\
rt(@param_list) ) {\n		print $param, \"\\n\";\n	}\\
n	print_debug_message( 'print_tool_params', 'End',\
 1 );\n}\n\n=head2 print_param_details()\n\nPrint \
details of a tool parameter.\n\n  &print_param_det\
ails($param_name);\n\n=cut\n\nsub print_param_deta\
ils {\n	print_debug_message( 'print_param_details'\
, 'Begin', 1 );\n	my $paramName = shift;\n	print_d\
ebug_message( 'print_param_details', 'paramName: '\
 . $paramName, 2 );\n	my $paramDetail = &rest_get_\
parameter_details($paramName);\n	print $paramDetai\
l->{'name'}, \"\\t\", $paramDetail->{'type'}, \"\\\
n\";\n	print $paramDetail->{'description'}, \"\\n\\
";\n	foreach my $value ( @{ $paramDetail->{'values\
'}->{'value'} } ) {\n		print $value->{'value'};\n	\
	if ( $value->{'defaultValue'} eq 'true' ) {\n			p\
rint \"\\t\", 'default';\n		}\n		print \"\\n\";\n	\
	print \"\\t\", $value->{'label'}, \"\\n\";\n	}\n	\
print_debug_message( 'print_param_details', 'End',\
 1 );\n}\n\n=head2 print_job_status()\n\nPrint sta\
tus of a job.\n\n  &print_job_status($job_id);\n\n\
=cut\n\nsub print_job_status {\n	print_debug_messa\
ge( 'print_job_status', 'Begin', 1 );\n	my $jobid \
= shift;\n	print_debug_message( 'print_job_status'\
, 'jobid: ' . $jobid, 1 );\n	if ( $outputLevel > 0\
 ) {\n		print STDERR 'Getting status for job ', $j\
obid, \"\\n\";\n	}\n	my $result = &rest_get_status\
($jobid);\n	print \"$result\\n\";\n	if ( $result e\
q 'FINISHED' && $outputLevel > 0 ) {\n		print STDE\
RR \"To get results: $scriptName --polljob --jobid\
 \" . $jobid\n		  . \"\\n\";\n	}\n	print_debug_mes\
sage( 'print_job_status', 'End', 1 );\n}\n\n=head2\
 print_result_types()\n\nPrint available result ty\
pes for a job.\n\n  &print_result_types($job_id);\\
n\n=cut\n\nsub print_result_types {\n	print_debug_\
message( 'result_types', 'Begin', 1 );\n	my $jobid\
 = shift;\n	print_debug_message( 'result_types', '\
jobid: ' . $jobid, 1 );\n	if ( $outputLevel > 0 ) \
{\n		print STDERR 'Getting result types for job ',\
 $jobid, \"\\n\";\n	}\n	my $status = &rest_get_sta\
tus($jobid);\n	if ( $status eq 'PENDING' || $statu\
s eq 'RUNNING' ) {\n		print STDERR 'Error: Job sta\
tus is ', $status,\n		  '. To get result types the\
 job must be finished.', \"\\n\";\n	}\n	else {\n		\
my (@resultTypes) = &rest_get_result_types($jobid)\
;\n		if ( $outputLevel > 0 ) {\n			print STDOUT 'A\
vailable result types:', \"\\n\";\n		}\n		foreach \
my $resultType (@resultTypes) {\n			print STDOUT $\
resultType->{'identifier'}, \"\\n\";\n			if ( defi\
ned( $resultType->{'label'} ) ) {\n				print STDOU\
T \"\\t\", $resultType->{'label'}, \"\\n\";\n			}\\
n			if ( defined( $resultType->{'description'} ) )\
 {\n				print STDOUT \"\\t\", $resultType->{'descr\
iption'}, \"\\n\";\n			}\n			if ( defined( $result\
Type->{'mediaType'} ) ) {\n				print STDOUT \"\\t\\
", $resultType->{'mediaType'}, \"\\n\";\n			}\n			\
if ( defined( $resultType->{'fileSuffix'} ) ) {\n	\
			print STDOUT \"\\t\", $resultType->{'fileSuffix\
'}, \"\\n\";\n			}\n		}\n		if ( $status eq 'FINISH\
ED' && $outputLevel > 0 ) {\n			print STDERR \"\\n\
\", 'To get results:', \"\\n\",\n			  \"  $scriptN\
ame --polljob --jobid \" . $params{'jobid'} . \"\\\
n\",\n			  \"  $scriptName --polljob --outformat <\
type> --jobid \"\n			  . $params{'jobid'} . \"\\n\\
";\n		}\n	}\n	print_debug_message( 'result_types',\
 'End', 1 );\n}\n\n=head2 submit_job()\n\nSubmit a\
 job to the service.\n\n  &submit_job($seq);\n\n=c\
ut\n\nsub submit_job {\n	print_debug_message( 'sub\
mit_job', 'Begin', 1 );\n\n	# Set input sequence\n\
	$tool_params{'sequence'} = shift;\n\n	# Load para\
meters\n	&load_params();\n\n	# Submit the job\n	my\
 $jobid = &rest_run( $params{'email'}, $params{'ti\
tle'}, \\%tool_params );\n\n	# Simulate sync/async\
 mode\n	if ( defined( $params{'async'} ) ) {\n		pr\
int STDOUT $jobid, \"\\n\";\n		if ( $outputLevel >\
 0 ) {\n			print STDERR\n			  \"To check status: $\
scriptName --status --jobid $jobid\\n\";\n		}\n	}\\
n	else {\n		if ( $outputLevel > 0 ) {\n			print ST\
DERR \"JobId: $jobid\\n\";\n		}\n		sleep 1;\n		&ge\
t_results($jobid);\n	}\n	print_debug_message( 'sub\
mit_job', 'End', 1 );\n}\n\n=head2 multi_submit_jo\
b()\n\nSubmit multiple jobs assuming input is a co\
llection of fasta formatted sequences.\n\n  &multi\
_submit_job();\n\n=cut\n\nsub multi_submit_job {\n\
	print_debug_message( 'multi_submit_job', 'Begin',\
 1 );\n	my $jobIdForFilename = 1;\n	$jobIdForFilen\
ame = 0 if ( defined( $params{'outfile'} ) );\n	my\
 (@filename_list) = ();\n\n	# Query sequence\n	if \
( defined( $ARGV[0] ) ) {    # Bare option\n		if (\
 -f $ARGV[0] || $ARGV[0] eq '-' ) {    # File\n			\
push( @filename_list, $ARGV[0] );\n		}\n	}\n	if ( \
$params{'sequence'} ) {                   # Via --\
sequence\n		if ( -f $params{'sequence'} || $params\
{'sequence'} eq '-' ) {    # File\n			push( @filen\
ame_list, $params{'sequence'} );\n		}\n	}\n\n	$/ =\
 '>';\n	foreach my $filename (@filename_list) {\n	\
	open( my $INFILE, '<', $filename )\n		  or die \"\
Error: unable to open file $filename ($!)\";\n		wh\
ile (<$INFILE>) {\n			my $seq = $_;\n			$seq =~ s/\
>$//;\n			if ( $seq =~ m/(\\S+)/ ) {\n				print ST\
DERR \"Submitting job for: $1\\n\"\n				  if ( $ou\
tputLevel > 0 );\n				$seq = '>' . $seq;\n				&pri\
nt_debug_message( 'multi_submit_job', $seq, 11 );\\
n				&submit_job($seq);\n				$params{'outfile'} = \
undef if ( $jobIdForFilename == 1 );\n			}\n		}\n	\
	close $INFILE;\n	}\n	print_debug_message( 'multi_\
submit_job', 'End', 1 );\n}\n\n=head2 list_file_su\
bmit_job()\n\nSubmit multiple jobs using a file co\
ntaining a list of entry identifiers as \ninput.\n\
\n  &list_file_submit_job($list_filename)\n\n=cut\\
n\nsub list_file_submit_job {\n	my $filename      \
   = shift;\n	my $jobIdForFilename = 1;\n	$jobIdFo\
rFilename = 0 if ( defined( $params{'outfile'} ) )\
;\n\n	# Iterate over identifiers, submitting each \
job\n	open( my $LISTFILE, '<', $filename )\n	  or \
die 'Error: unable to open file ' . $filename . ' \
(' . $! . ')';\n	while (<$LISTFILE>) {\n		my $line\
 = $_;\n		chomp($line);\n		if ( $line ne '' ) {\n	\
		&print_debug_message( 'list_file_submit_job', 'l\
ine: ' . $line, 2 );\n			if ( $line =~ m/\\w:\\w/ \
) {    # Check this is an identifier\n				print ST\
DERR \"Submitting job for: $line\\n\"\n				  if ( \
$outputLevel > 0 );\n				&submit_job($line);\n			}\
\n			else {\n				print STDERR\n\"Warning: line \\\\
"$line\\\" is not recognised as an identifier\\n\"\
;\n			}\n		}\n		$params{'outfile'} = undef if ( $j\
obIdForFilename == 1 );\n	}\n	close $LISTFILE;\n}\\
n\n=head2 load_data()\n\nLoad sequence data from f\
ile or option specified on the command-line.\n\n  \
&load_data();\n\n=cut\n\nsub load_data {\n	print_d\
ebug_message( 'load_data', 'Begin', 1 );\n	my $ret\
Seq;\n\n	# Query sequence\n	if ( defined( $ARGV[0]\
 ) ) {    # Bare option\n		if ( -f $ARGV[0] || $AR\
GV[0] eq '-' ) {    # File\n			$retSeq = &read_fil\
e( $ARGV[0] );\n		}\n		else {                     \
                # DB:ID or sequence\n			$retSeq = \
$ARGV[0];\n		}\n	}\n	if ( $params{'sequence'} ) { \
                  # Via --sequence\n		if ( -f $par\
ams{'sequence'} || $params{'sequence'} eq '-' ) { \
   # File\n			$retSeq = &read_file( $params{'seque\
nce'} );\n		}\n		else {    # DB:ID or sequence\n		\
	$retSeq = $params{'sequence'};\n		}\n	}\n	print_d\
ebug_message( 'load_data', 'End', 1 );\n	return $r\
etSeq;\n}\n\n=head2 load_params()\n\nLoad job para\
meters from command-line options.\n\n  &load_param\
s();\n\n=cut\n\nsub load_params {\n	print_debug_me\
ssage( 'load_params', 'Begin', 1 );\n\n	# Database\
(s) to search\n	my (@dbList) = split /[ ,]/, $para\
ms{'database'};\n	$tool_params{'database'} = \\@db\
List;\n\n	# Compatability options, old command-lin\
e.\n	if(!$tool_params{'viewfilter'} && $params{'ec\
hofilter'}) {\n		$tool_params{'viewfilter'} = 'tru\
e';\n	}\n	if(!$tool_params{'alignments'} && $param\
s{'numal'}) {\n		$tool_params{'alignments'} = $par\
ams{'numal'};\n	}\n	# TODO: set alignment format o\
ption to get NCBI BLAST XML.\n	if($params{'appxml'\
}) {\n		$tool_params{'align'} = '';\n	}\n\n	print_\
debug_message( 'load_params', 'End', 1 );\n}\n\n=h\
ead2 client_poll()\n\nClient-side job polling.\n\n\
  &client_poll($job_id);\n\n=cut\n\nsub client_pol\
l {\n	print_debug_message( 'client_poll', 'Begin',\
 1 );\n	my $jobid  = shift;\n	my $status = 'PENDIN\
G';\n\n	my $errorCount = 0;\n	while ($status eq 'R\
UNNING'\n		|| $status eq 'PENDING'\n		|| ( $status\
 eq 'ERROR' && $errorCount < 2 ) )\n	{\n		$status \
= rest_get_status($jobid);\n		print STDERR \"$stat\
us\\n\" if ( $outputLevel > 0 );\n		if ( $status e\
q 'ERROR' ) {\n			$errorCount++;\n		}\n		elsif ( $\
errorCount > 0 ) {\n			$errorCount--;\n		}\n		if (\
   $status eq 'RUNNING'\n			|| $status eq 'PENDING\
'\n			|| $status eq 'ERROR' )\n		{\n\n			# Wait be\
fore polling again.\n			sleep $checkInterval;\n		}\
\n	}\n	print_debug_message( 'client_poll', 'End', \
1 );\n	return $status;\n}\n\n=head2 get_results()\\
n\nGet the results for a job identifier.\n\n  &get\
_results($job_id);\n\n=cut\n\nsub get_results {\n	\
print_debug_message( 'get_results', 'Begin', 1 );\\
n	my $jobid = shift;\n	print_debug_message( 'get_r\
esults', 'jobid: ' . $jobid, 1 );\n\n	# Verbose\n	\
if ( $outputLevel > 1 ) {\n		print 'Getting result\
s for job ', $jobid, \"\\n\";\n	}\n\n	# Check stat\
us, and wait if not finished\n	client_poll($jobid)\
;\n\n	# Use JobId if output file name is not defin\
ed\n	unless ( defined( $params{'outfile'} ) ) {\n	\
	$params{'outfile'} = $jobid;\n	}\n\n	# Get list o\
f data types\n	my (@resultTypes) = rest_get_result\
_types($jobid);\n\n	# Get the data and write it to\
 a file\n	if ( defined( $params{'outformat'} ) ) {\
    # Specified data type\n		my $selResultType;\n	\
	foreach my $resultType (@resultTypes) {\n			if ( \
$resultType->{'identifier'} eq $params{'outformat'\
} ) {\n				$selResultType = $resultType;\n			}\n		\
}\n		if ( defined($selResultType) ) {\n			my $resu\
lt =\n			  rest_get_result( $jobid, $selResultType\
->{'identifier'} );\n			if ( $params{'outfile'} eq\
 '-' ) {\n				write_file( $params{'outfile'}, $res\
ult );\n			}\n			else {\n				write_file(\n					$pa\
rams{'outfile'} . '.'\n					  . $selResultType->{'\
identifier'} . '.'\n					  . $selResultType->{'fil\
eSuffix'},\n					$result\n				);\n			}\n		}\n		els\
e {\n			die 'Error: unknown result format \"' . $p\
arams{'outformat'} . '\"';\n		}\n	}\n	else {    # \
Data types available\n		      # Write a file for e\
ach output type\n		for my $resultType (@resultType\
s) {\n			if ( $outputLevel > 1 ) {\n				print STDE\
RR 'Getting ', $resultType->{'identifier'}, \"\\n\\
";\n			}\n			my $result = rest_get_result( $jobid,\
 $resultType->{'identifier'} );\n			if ( $params{'\
outfile'} eq '-' ) {\n				write_file( $params{'out\
file'}, $result );\n			}\n			else {\n				write_fil\
e(\n					$params{'outfile'} . '.'\n					  . $resul\
tType->{'identifier'} . '.'\n					  . $resultType-\
>{'fileSuffix'},\n					$result\n				);\n			}\n		}\\
n	}\n	print_debug_message( 'get_results', 'End', 1\
 );\n}\n\n=head2 read_file()\n\nRead a file into a\
 scalar. The special filename '-' can be used to r\
ead from \nstandard input (STDIN).\n\n  my $data =\
 &read_file($filename);\n\n=cut\n\nsub read_file {\
\n	print_debug_message( 'read_file', 'Begin', 1 );\
\n	my $filename = shift;\n	print_debug_message( 'r\
ead_file', 'filename: ' . $filename, 2 );\n	my ( $\
content, $buffer );\n	if ( $filename eq '-' ) {\n	\
	while ( sysread( STDIN, $buffer, 1024 ) ) {\n			$\
content .= $buffer;\n		}\n	}\n	else {    # File\n	\
	open( my $FILE, '<', $filename )\n		  or die \"Er\
ror: unable to open input file $filename ($!)\";\n\
		while ( sysread( $FILE, $buffer, 1024 ) ) {\n			\
$content .= $buffer;\n		}\n		close($FILE);\n	}\n	p\
rint_debug_message( 'read_file', 'End', 1 );\n	ret\
urn $content;\n}\n\n=head2 write_file()\n\nWrite d\
ata to a file. The special filename '-' can be use\
d to write to \nstandard output (STDOUT).\n\n  &wr\
ite_file($filename, $data);\n\n=cut\n\nsub write_f\
ile {\n	print_debug_message( 'write_file', 'Begin'\
, 1 );\n	my ( $filename, $data ) = @_;\n	print_deb\
ug_message( 'write_file', 'filename: ' . $filename\
, 2 );\n	if ( $outputLevel > 0 ) {\n		print STDERR\
 'Creating result file: ' . $filename . \"\\n\";\n\
	}\n	if ( $filename eq '-' ) {\n		print STDOUT $da\
ta;\n	}\n	else {\n		open( my $FILE, '>', $filename\
 )\n		  or die \"Error: unable to open output file\
 $filename ($!)\";\n		syswrite( $FILE, $data );\n	\
	close($FILE);\n	}\n	print_debug_message( 'write_f\
ile', 'End', 1 );\n}\n\n=head2 usage()\n\nPrint pr\
ogram usage message.\n\n  &usage();\n\n=cut\n\nsub\
 usage {\n	print STDERR <<EOF\nWU-BLAST\n========\\
n   \nRapid sequence database search programs util\
izing the BLAST algorithm\n    \n[Required]\n\n  -\
p, --program      : str  : BLAST program to use, s\
ee --paramDetail program\n  -D, --database     : s\
tr  : database(s) to search, space separated. See\\
n                              --paramDetail datab\
ase\n      --stype        : str  : query sequence \
type, see --paramDetail stype\n  seqFile          \
  : file : query sequence (\"-\" for STDIN, \\@fil\
ename for\n                              identifie\
r list file)\n\n[Optional]\n\n  -m, --matrix      \
 : str  : scoring matrix, see --paramDetail matrix\
\n  -e, --exp          : real : 0<E<= 1000. Statis\
tical significance threshold \n                   \
           for reporting database sequence matches\
.\n  -e, --viewfilter   :      : display the filte\
red query sequence\n  -f, --filter       : str  : \
filter the query sequence for low complexity \n   \
                           regions, see --paramDet\
ail filter\n  -A, --align        : int  : pairwise\
 alignment format, see --paramDetail align\n  -s, \
--scores       : int  : number of scores to be rep\
orted\n  -b, --alignments   : int  : number of ali\
gnments to report\n  -S, --sensitivity  : str  : s\
ensitivity of the search, \n                      \
        see --paramDetail sensitivity\n  -t, --sor\
t	     : str  : sort order for hits, see --paramDe\
tail sort\n  -T, --stats        : str  : statistic\
al model, see --paramDetail stats\n  -d, --strand \
      : str  : DNA strand to search with,\n       \
                       see --paramDetail strand\n \
 -c, --topcombon    : str  : consistent sets of HS\
Ps\n      --multifasta   :      : treat input as a\
 set of fasta formatted sequences\n\n[General]\n\n\
  -h, --help        :      : prints this help text\
\n      --async       :      : forces to make an a\
synchronous query\n      --email       : str  : e-\
mail address\n      --title       : str  : title f\
or job\n      --status      :      : get job statu\
s\n      --resultTypes :      : get available resu\
lt types for job\n      --polljob     :      : pol\
l for the status of a job\n      --jobid       : s\
tr  : jobid that was returned when an asynchronous\
 job \n                             was submitted.\
\n      --outfile     : str  : file name for resul\
ts (default is jobid;\n                           \
  \"-\" for STDOUT)\n      --outformat   : str  : \
result format to retrieve\n      --params      :  \
    : list input parameters\n      --paramDetail :\
 str  : display details for input parameter\n     \
 --quiet       :      : decrease output\n      --v\
erbose     :      : increase output\n      --trace\
       :      : show SOAP messages being interchan\
ged \n   \nSynchronous job:\n\n  The results/error\
s are returned as soon as the job is finished.\n  \
Usage: $scriptName --email <your\\@email> [options\
...] seqFile\n  Returns: results as an attachment\\
n\nAsynchronous job:\n\n  Use this if you want to \
retrieve the results at a later time. The results \
\n  are stored for up to 24 hours. 	\n  Usage: $sc\
riptName --async --email <your\\@email> [options..\
.] seqFile\n  Returns: jobid\n\n  Use the jobid to\
 query for the status of the job. If the job is fi\
nished, \n  it also returns the results/errors.\n \
 Usage: $scriptName --polljob --jobid <jobId> [--o\
utfile string]\n  Returns: string indicating the s\
tatus of the job and if applicable, results \n  as\
 an attachment.\n\nFurther information:\n\n  http:\
//www.ebi.ac.uk/Tools/webservices/services/sss/wu_\
blast_rest\n  http://www.ebi.ac.uk/Tools/webservic\
es/tutorials/perl\n\nSupport/Feedback:\n\n  http:/\
/www.ebi.ac.uk/support/\nEOF\n}\n\n=head1 FEEDBACK\
/SUPPORT\n\nPlease contact us at L<http://www.ebi.\
ac.uk/support/> if you have any \nfeedback, sugges\
tions or issues with the service or this client.\n\
\n=cut\n","\n\n\nmy $PROBTRESH = 0.3;# base pairs \
below this prob threshold will be ignored\nmy $WEI\
GHT = 100.0; # float!!\nmy $NUCALPH = \"ACGTUNRYMK\
SWHBVD\";\nuse vars qw($NUCALPH $WEIGHT);\n\nmy $m\
yname = basename($0);\n\nuse strict;\nuse warnings\
;\n\nuse File::Basename;\nuse Getopt::Long;\nuse F\
ile::Glob ':glob';\nuse File::Spec;\nuse File::Tem\
p qw/ tempfile tempdir /;\n\n\n\n\nsub tcoffeelib_\
header($;$)\n{\n    my ($nseq, $fd) = @_;\n    if \
(! defined($fd)) {\n        $fd = *STDOUT;\n    }\\
n    printf $fd \"! TC_LIB_FORMAT_01\\n\";\n    pr\
intf $fd \"%d\\n\", $nseq;\n}\n\n\nsub tcoffeelib_\
header_addseq($$;$)\n{\n    my ($id, $seq, $fd) = \
@_;\n    if (! defined($fd)) {\n        $fd = *STD\
OUT;\n    }\n    printf $fd \"%s %d %s\\n\", $id, \
length($seq), $seq;\n}\n\n\nsub tcoffeelib_comment\
($;$)\n{\n    my ($comment, $fd) = @_;\n    if (! \
defined($fd)) {\n        $fd = *STDOUT;\n    }\n  \
  printf $fd \"!\" . $comment . \"\\n\";\n}\n\n\ns\
ub tcoffeelib_struct($$$;$)\n{\n    my ($nseq, $le\
n, $bpm, $fd) = @_;\n\n    if (! defined($fd)) {\n\
        $fd = *STDOUT;\n    }\n\n    # output base\
pair indices with fixed weight\n    printf $fd \"#\
%d %d\\n\", $nseq, $nseq;\n    # output basepairs \
(only once) and with unit-offset\n    for (my $i=0\
; $i<$len; $i++) {\n        for (my $j=$i+1; $j<$l\
en; $j++) {\n            if (! defined($bpm->[$i][\
$j])) {\n                print STDERR \"ERROR: \\$\
bpm->[$i][$j] undefined\\n\";\n            }\n    \
        if ($bpm->[$i][$j]>0) {\n                p\
rint $fd $i+1;\n                print $fd \" \";\n\
                print $fd $j+1;\n                p\
rint $fd \" \" . $bpm->[$i][$j] . \"\\n\";\n      \
      }\n        }\n    }\n}\n\n\nsub tcoffeelib_f\
ooter(;$)\n{\n    my ($fd) = @_;\n    if (! define\
d($fd)) {\n        $fd = *STDOUT;\n    }\n    prin\
t $fd \"! SEQ_1_TO_N\\n\";\n}\n\n\n    \nsub plfol\
d($$$)\n{    \n    my ($id, $seq, $probtresh) = @_\
;\n    my (@struct);# return\n    my ($templ, $fht\
mp, $fnametmp, $cmd, $ctr, $window_size);\n    our\
 $ntemp++;\n    \n    $templ = $myname . \".pid-\"\
 . $$ .$ntemp .\".XXXXXX\";\n    ($fhtmp, $fnametm\
p) = tempfile($templ, UNLINK => 1); \n    print $f\
htmp \">$id\\n$seq\\n\";\n\n    # --- init basepai\
r array\n    #\n    for (my $i=0; $i<length($seq);\
 $i++) {\n        for (my $j=$i+1; $j<length($seq)\
; $j++) {\n            $struct[$i][$j]=0;\n       \
 }\n    }\n\n\n    # --- call rnaplfold and drop a\
 readme\n    #\n    $window_size=(length($seq)<70)\
?length($seq):70;\n    $cmd = \"RNAplfold -W $wind\
ow_size < $fnametmp >/dev/null\";\n    system($cmd\
);\n    \n    if ($? != 0) {\n        printf STDER\
R \"ERROR: RNAplfold ($cmd) exited with error stat\
us %d\\n\", $? >> 8;\n        return;\n    }\n    \
#unlink($fnametmp);\n    my $fps = sprintf(\"%s_dp\
.ps\", $id); # check long name\n    \n    if (! -s\
 $fps) {\n      {\n\n	$fps = sprintf(\"%s_dp.ps\",\
 substr($id,0,12)); # check short name\n 	if (! -s\
 $fps)\n	  {\n	    die(\"couldn't find expected fi\
le $fps\\n\");\n	    return;\n	  }\n      }\n    }\
\n\n    \n    # --- read base pairs from created p\
ostscript\n    #\n    open(FH, $fps);\n    while (\
my $line = <FH>) {\n        my ($nti, $ntj, $prob)\
;\n        chomp($line);        \n        # line: \
bp bp sqrt-prob ubox\n        my @match = ($line =\
~ m/^([0-9]+) +([0-9]+) +([0-9\\.]+) +ubox$/);\n  \
      if (scalar(@match)) {\n            $nti=$1;\\
n            $ntj=$2;\n            $prob=$3*$3;# p\
rob stored as square root\n\n            if ($prob\
>$probtresh) {\n                #printf STDERR \"\\
\$struct[$nti][$ntj] sqrtprob=$3 prob=$prob > $pro\
btresh\\n\";\n                $struct[$nti-1][$ntj\
-1] = $WEIGHT\n            }\n            # store \
with zero-offset\n        }\n    }\n    close(FH);\
\n\n    # remove or gzi postscript\n    #\n    unl\
ink($fps);\n    #\n    # or gzip\n    #$cmd = \"gz\
ip -qf $fps\";\n    #system($cmd);\n    #if ($? !=\
 0) {\n    #    printf STDERR \"ERROR: gzip ($cmd)\
 exited with error status %d\\n\", $? >> 8;\n    #\
}\n\n    return \\@struct;\n}\n\n\n\n\n\nsub rnase\
qfmt($)\n{\n    my ($seq) = @_;\n    # remove gaps\
\n    $seq =~ s/-//g;\n    # uppercase RNA\n    $s\
eq = uc($seq);\n    # T -> U\n    $seq =~ s/T/U/g;\
\n    # check for invalid charaters\n    $_ = $seq\
;\n    s/[^$NUCALPH]//g;\n    return $_;\n}\n\n\n\\
n\nsub usage(;$)\n{    \n    my ($errmsg) = @_;\n \
   if ($errmsg) {\n        print STDERR \"ERROR: $\
errmsg\\n\";\n    }\n    print STDERR << \"EOF\";\\
n$myname:\n Creates a T-Coffee RNA structure libra\
ry from RNAplfold prediction.\n See FIXME:citation\
\nUsage:\n $myname -in seq_file -out tcoffee_lib\n\
EOF\n    exit(1);\n}\n\nsub read_fasta_seq \n  {\n\
    my $f=$_[0];\n    my %hseq;\n    my (@seq, @co\
m, @name);\n    my ($a, $s,$nseq);\n\n    open (F,\
 $f);\n    while (<F>)\n      {\n	$s.=$_;\n      }\
\n    close (F);\n\n    \n    @name=($s=~/>(\\S*).\
*\\n[^>]*/g);\n    \n    @seq =($s=~/>.*.*\\n([^>]\
*)/g);\n    @com =($s=~/>(\\S*)(.*)\\n([^>]*)/g);\\
n\n\n    $nseq=$#name+1;\n  \n    for ($a=0; $a<$n\
seq; $a++)\n      {\n	my $n=$name[$a];\n	my $s;\n	\
$hseq{$n}{name}=$n;\n	$s=$seq[$a];$s=~s/\\s//g;\n	\
\n	$hseq{$n}{seq}=$s;\n	$hseq{$n}{com}=$com[$a];\n\
      }\n    return %hseq;\n  }\n\n\n\n\n\n\n\nmy \
$fmsq = \"\";\nmy $flib = \"\";\nmy %OPTS;\nmy %se\
q;\nmy ($id, $nseq, $i);\nmy @nl;\n\nGetOptions(\"\
in=s\" => \\$fmsq, \"out=s\" => \\$flib);\n\nif (!\
 -s $fmsq) {\n    usage(\"empty or non-existant fi\
le \\\"$fmsq\\\"\")\n}\nif (length($flib)==0) {\n \
   usage(\"empty out-filename\")\n}\n\n\n\n\n\n\n%\
seq=read_fasta_seq($fmsq);\n\n\n@nl=keys(%seq);\n\\
n$nseq=$#nl+1;\nopen FD_LIB, \">$flib\" or die \"c\
an't open $flib!\";\ntcoffeelib_header($nseq, *FD_\
LIB);\nforeach $id (keys (%seq))\n  {\n    my ($se\
q, $fmtseq);\n    \n    $seq = $seq{$id}{seq};\n  \
  \n    $fmtseq = rnaseqfmt($seq);# check here, fo\
rmatting for folding important later\n    if (leng\
th($seq)!=length($fmtseq)) {\n        print STDERR\
 \"ERROR: invalid sequence $id is not an RNA seque\
nce. read seq is: $seq\\n\";\n        exit\n      \
}\n   \n    tcoffeelib_header_addseq($id, uc($seq)\
, *FD_LIB);\n  }\ntcoffeelib_comment(\"generated b\
y $myname on \" . localtime(), *FD_LIB);\n\n\n\n$i\
=0;\nforeach $id (keys (%seq))\n  {\n    my ($clea\
nid, $seq, $bpm);\n    $seq=$seq{$id}{seq};\n    $\
cleanid = $id;\n    $cleanid =~ s,[/ ],_,g;# neede\
d for rnaplfold\n    $seq = rnaseqfmt($seq);\n    \
\n    $bpm = plfold($cleanid, rnaseqfmt($seq), $PR\
OBTRESH);       \n    \n    tcoffeelib_struct($i+1\
, length($seq), $bpm, *FD_LIB);\n    $i++;\n}\n\n\\
ntcoffeelib_footer(*FD_LIB);\nclose FD_LIB;\nexit \
(0);\n\n","\n\n\n\n\n$cmd=join ' ', @ARGV;\nif ($c\
md=~/-infile=(\\S+)/){ $seqfile=$1;}\nif ($cmd=~/-\
outfile=(\\S+)/){ $libfile=$1;}\n\n\n\n%s=read_fas\
ta_seq ($seqfile);\n\nopen (F, \">$libfile\");\nfo\
reach $name (keys (%s))\n  {\n    my $tclib=\"$nam\
e.RNAplfold_tclib\";\n    print (F \">$name _F_ $t\
clib\\n\");\n    seq2RNAplfold2tclib ($name, $s{$n\
ame}{seq}, $tclib);\n  }\nclose (F);\nexit (EXIT_S\
UCCESS);\n\nsub seq2RNAplfold2tclib\n  {\n    my (\
$name, $seq, $tclib)=@_;\n    my ($tmp);\n    $n++\
;\n    $tmp=\"tmp4seq2RNAplfold_tclib.$$.$n.pep\";\
\n    open (RF, \">$tmp\");\n    print (RF \">$nam\
e\\n$seq\\n\");\n    close (RF);\n    \n    system\
 \"t_coffee -other_pg RNAplfold2tclib.pl -in=$tmp \
-out=$tclib\";\n    \n    unlink ($tmp);\n    retu\
rn $tclib;\n  }\n    \n    \nsub read_fasta_seq \n\
  {\n    my $f=@_[0];\n    my %hseq;\n    my (@seq\
, @com, @name);\n    my ($a, $s,$nseq);\n\n    ope\
n (F, $f);\n    while (<F>)\n      {\n	$s.=$_;\n  \
    }\n    close (F);\n\n    \n    @name=($s=~/>(\\
\S*).*\\n[^>]*/g);\n    \n    @seq =($s=~/>.*.*\\n\
([^>]*)/g);\n    @com =($s=~/>\\S*(.*)\\n([^>]*)/g\
);\n\n    \n    $nseq=$#name+1;\n    \n    for ($a\
=0; $a<$nseq; $a++)\n      {\n	my $n=$name[$a];\n	\
$hseq{$n}{name}=$n;\n	$hseq{$n}{seq}=$seq[$a];\n	$\
hseq{$n}{com}=$com[$a];\n      }\n    return %hseq\
;\n  }\n","use Getopt::Long;\nuse File::Path;\nuse\
 Env;\nuse FileHandle;\nuse Cwd;\nuse Sys::Hostnam\
e;\nour $PIDCHILD;\nour $ERROR_DONE;\nour @TMPFILE\
_LIST;\nour $EXIT_FAILURE=1;\nour $EXIT_SUCCESS=0;\
\n\nour $REFDIR=getcwd;\nour $EXIT_SUCCESS=0;\nour\
 $EXIT_FAILURE=1;\n\nour $PROGRAM=\"tc_generic_met\
hod.pl\";\nour $CL=$PROGRAM;\n\nour $CLEAN_EXIT_ST\
ARTED;\nour $debug_lock=$ENV{\"DEBUG_LOCK\"};\nour\
 $LOCKDIR=$ENV{\"LOCKDIR_4_TCOFFEE\"};\nif (!$LOCK\
DIR){$LOCKDIR=getcwd();}\nour $ERRORDIR=$ENV{\"ERR\
ORDIR_4_TCOFFEE\"};\nour $ERRORFILE=$ENV{\"ERRORFI\
LE_4_TCOFFEE\"};\n&set_lock ($$);\nif (isshellpid(\
getppid())){lock4tc(getppid(), \"LLOCK\", \"LSET\"\
, \"$$\\n\");}\n      \nour $print;\nmy ($fmsq1, $\
fmsq2, $output, $outfile, $arch, $psv, $hmmtop_hom\
e, $trim, $cov, $sample, $mode, $gor_home, $gor_se\
q, $gor_obs);\n\nGetOptions(\"-in=s\" => \\$fmsq1,\
\"-output=s\" =>\\$output ,\"-out=s\" => \\$outfil\
e, \"-arch=s\" => \\$arch,\"-psv=s\" => \\$psv, \"\
-hmmtop_home=s\", \\$hmmtop_home,\"-trim=s\" =>\\$\
trim ,\"-print=s\" =>\\$print,\"-cov=s\" =>\\$cov \
, \"-sample=s\" =>\\$sample, \"-mode=s\" =>\\$mode\
, \"-gor_home=s\"=>\\$gor_home, \"-gor_seq=s\"=>\\\
$gor_seq,\"-gor_obs=s\"=>\\$gor_obs);\n\n\nif (!$m\
ode){$mode = \"hmmtop\"}\nelsif ($mode eq \"hmmtop\
\"){;}\nelsif ($mode eq \"gor\"){;}\nelse {myexit(\
flush_error (\"-mode=$mode is unknown\"));}\n\n\no\
ur $HOME=$ENV{\"HOME\"};\nour $MCOFFEE=($ENV{\"MCO\
FFEE_4_TCOFFEE\"})?$ENV{\"MCOFFEE_4_TCOFFEE\"}:\"$\
HOME/.t_coffee/mcoffee\";\n\nif ($mode eq \"hmmtop\
\")\n  {\n    check_configuration (\"hmmtop\");\n \
   if (-e $arch){$ENV{'HMMTOP_ARCH'}=$arch;}\n    \
elsif (-e $ENV{HMMTOP_ARCH}){$arch=$ENV{HMMTOP_ARC\
H};}\n    elsif (-e \"$MCOFFEE/hmmtop.arch\"){$arc\
h=$ENV{'HMMTOP_ARCH'}=\"$MCOFFEE/hmmtop.arch\";}\n\
    elsif (-e \"$hmmtop_home/hmmtop.arc\"){$arch=$\
ENV{'HMMTOP_ARCH'}=\"$hmmtop_home/hmmtop.arc\";}\n\
    else {myexit(flush_error ( \"Could not find AR\
CH file for hmmtop\"));}\n    \n    \n    if (-e $\
psv){$ENV{'HMMTOP_PSV'}=$psv;}\n    elsif (-e $ENV\
{HMMTOP_PSV}){$psv=$ENV{HMMTOP_PSV};}\n    elsif (\
-e \"$MCOFFEE/hmmtop.psv\"){$psv=$ENV{'HMMTOP_PSV'\
}=\"$MCOFFEE/hmmtop.psv\";}\n    elsif (-e \"$hmmt\
op_home/hmmtop.psv\"){$psv=$ENV{'HMMTOP_PSV'}=\"$h\
mmtop_home/hmmtop.psv\";}\n    else {myexit(flush_\
error ( \"Could not find PSV file for hmmtop\"));}\
\n  }\nelsif ($mode eq \"gor\")\n  {\n    our $GOR\
_SEQ;\n    our $GOR_OBS;\n    \n    check_configur\
ation (\"gorIV\");\n    if (-e $gor_seq){$GOR_SEQ=\
$gor_seq;}\n    elsif (-e $ENV{GOR_SEQ}){$GOR_SEQ=\
$ENV{GOR_SEQ};}\n    elsif (-e \"$MCOFFEE/New_KS.2\
67.seq\"){$GOR_SEQ=\"$MCOFFEE/New_KS.267.seq\";}\n\
    elsif (-e \"$gor_home/New_KS.267.seq\"){$GOR_S\
EQ=\"$gor_home/New_KS.267.seq\";}\n    else {myexi\
t(flush_error ( \"Could not find SEQ file for gor\\
"));}\n\n    if (-e $gor_obs){$GOR_OBS=$gor_obs;}\\
n    elsif (-e $ENV{GOR_OBS}){$GOR_OBS=$ENV{GOR_OB\
S};}\n    elsif (-e \"$MCOFFEE/New_KS.267.obs\"){$\
GOR_OBS=\"$MCOFFEE/New_KS.267.obs\";}\n    elsif (\
-e \"$gor_home/New_KS.267.obs\"){$GOR_OBS=\"$gor_h\
ome/New_KS.267.obs\";}\n    else {myexit(flush_err\
or ( \"Could not find OBS file for gor\"));}\n  }\\
n\n\nif ( ! -e $fmsq1){myexit(flush_error (\"Could\
 Not Read Input file $fmsq1\"));}\n\n\nmy $fmsq2=v\
tmpnam();\nmy $fmsq3=vtmpnam();\nmy $tmpfile=vtmpn\
am();\nmy $predfile=vtmpnam();\n\nif ($trim){$trim\
_action=\" +trim _aln_%%$trim\\_K1 \";}\nif ($cov)\
 {$cov_action= \" +sim_filter _aln_c$cov \";}\n&sa\
fe_system(\"t_coffee -other_pg seq_reformat -in $f\
msq1 -action +convert 'BOUJXZ-' $cov_action $trim_\
action -output fasta_aln -out $fmsq2\");\nmy (%pre\
d, %seq, %predA);\n\n\n%seq=read_fasta_seq($fmsq2)\
;\n%seq=fasta2sample(\\%seq, $sample);\n\nif (1==2\
 && $mode eq \"hmmtop\" && $output eq \"cons\")\n \
 {\n    fasta2hmmtop_cons($outfile,\\%seq);\n  }\n\
else\n  {\n    %pred=fasta2pred(\\%seq, $mode);\n \
   %predA=pred2aln (\\%pred, \\%seq);\n    \n    \\
n    if (!$output || $output eq \"prediction\"){ou\
tput_fasta_seq (\\%predA, $outfile);}\n    elsif (\
$output eq \"color_html\"){pred2color (\\%pred,\\%\
seq, $outfile);}\n    elsif ($output eq \"cons\"){\
pred2cons($outfile,\\%predA);}\n    else {flush_er\
ror (\"$output is an unknown output mode\");}\n  }\
\n\nsub fasta2sample\n  {\n    my $SR=shift;\n    \
my $it=shift;\n    my %S=%$SR;\n    \n    my $seq=\
index2seq_name (\\%S, 1);\n    my $l=length($S{$se\
q}{seq});\n    my @sl=keys(%S);\n    my $nseq=$#sl\
+1;\n    my $index=$nseq;\n  \n    if (!$sample) {\
return %S;}\n    for (my $a=0; $a<$it; $a++)\n    \
  {\n	my $newseq=\"\";\n	my $nname=\"$seq\\_sample\
d_$index\";\n	for (my $p=0; $p<$l; $p++)\n	  {\n	 \
   my $i=int(rand($nseq));\n	    \n	    my $name =\
 $sl[$i];\n	    my $seq=$S{$name}{seq};\n	    my $\
r=substr ($seq, $p, 1);\n	    $newseq.=$r;\n	  }\n\
	$S{$nname}{name}=$nname;\n	$S{$nname}{seq}=$newse\
q;\n	$S{$nname}{com}=\"sampled\";\n	$S{$nname}{ind\
ex}=++$index;\n      }\n    return %S;\n  }\n	    \
  \nsub fasta2pred\n  {\n    my $s=shift;\n    my \
$mode=shift;\n\n    if ( $mode eq \"hmmtop\"){retu\
rn fasta2hmmtop_pred($s);}\n    elsif ($mode eq \"\
gor\"){return fasta2gor_pred ($s);}\n  }\nsub fast\
a2hmmtop_cons\n  {\n    my $outfile=shift;\n    my\
 $SR=shift;\n    \n    my $o = new FileHandle;\n  \
  my $i = new FileHandle;\n    my $tmp_in =vtmpnam\
();\n    my $tmp_out=vtmpnam();\n    my %seq=%$SR;\
\n    my %pred;\n    my $N=keys(%seq);\n    \n    \
output_fasta_seq (\\%seq,$tmp_in, \"seq\");\n    `\
hmmtop -pi=mpred -if=$tmp_in -sf=FAS -pl 2>/dev/nu\
ll >$tmp_out`;\n    open ($o, \">$outfile\");\n   \
 open ($i, \"$tmp_out\");\n    while (<$i>)\n     \
 {\n	my $l=$_;\n	if (($l=~/>HP\\:\\s+(\\d+)\\s+(.*\
)/)){my $line=\">$2 NSEQ: $N\\n\";print $o \"$line\
\";}\n	elsif ( ($l=~/.*pred(.*)/))  {my $line=\"$1\
\\n\";print $o \"$line\";}\n      }\n    close ($o\
);\n    close ($i);\n    return read_fasta_seq($tm\
p);\n  }\nsub fasta2hmmtop_pred\n  {\n    my $SR=s\
hift;\n    my $o = new FileHandle;\n    my $i = ne\
w FileHandle;\n    my $tmp    =vtmpnam();\n    my \
$tmp_in =vtmpnam();\n    my $tmp_out=vtmpnam();\n \
   my %seq=%$SR;\n    my %pred;\n    \n\n    outpu\
t_fasta_seq (\\%seq,$tmp_in, \"seq\");\n    `hmmto\
p -if=$tmp_in -sf=FAS -pl 2>/dev/null >$tmp_out`;\\
n    open ($o, \">$tmp\");\n    open ($i, \"$tmp_o\
ut\");\n    while (<$i>)\n      {\n	my $l=$_;\n	if\
 (($l=~/>HP\\:\\s+(\\d+)\\s+(.*)/)){my $line=\">$2\
\\n\";print $o \"$line\";}\n	elsif ( ($l=~/.*pred(\
.*)/))  {my $line=\"$1\\n\";print $o \"$line\";}\n\
      }\n    close ($o);\n    close ($i);\n    ret\
urn read_fasta_seq($tmp);\n  }\n    \n	\n	\n	    \\
n	\n	\n\n	\nsub fasta2gor_pred\n  {\n    my $SR=sh\
ift;\n    my $o = new FileHandle;\n    my $i = new\
 FileHandle;\n    my $tmp    =vtmpnam();\n    my $\
tmp_in =vtmpnam();\n    my $tmp_out=vtmpnam();\n  \
  my %seq=%$SR;\n    my %pred;\n    \n\n    output\
_fasta_seq (\\%seq,$tmp_in, \"seq\");\n    `gorIV \
-prd $tmp_in -seq $GOR_SEQ -obs $GOR_OBS >$tmp_out\
`;\n    open ($o, \">$tmp\");\n    open ($i, \"$tm\
p_out\");\n    while (<$i>)\n      {\n	my $l=$_;\n\
\n	\n	if ( $l=~/>/){print $o \"$l\";}\n	elsif ( $l\
=~/Predicted Sec. Struct./){$l=~s/Predicted Sec. S\
truct\\.//;print $o \"$l\";}\n      }\n    close (\
$o);\n    close ($i);\n    return read_fasta_seq($\
tmp);\n  }\n			\n			     \nsub index2seq_name\n  {\
\n    \n    my $SR=shift;\n    my $index=shift;\n \
   \n    \n    my %S=%$SR;\n    \n    foreach my $\
s (%S)\n      {\n	if ( $S{$s}{index}==$index){retu\
rn $s;}\n      }\n    return \"\";\n  }\n\nsub pre\
d2cons\n  {\n    my $outfile=shift;\n    my $predR\
=shift;\n    my $seq=shift;\n    my %P=%$predR;\n \
   my %C;\n    my ($s,@r,$nseq);\n    my $f= new F\
ileHandle;\n\n    open ($f, \">$outfile\");\n\n   \
 if (!$seq){$seq=index2seq_name(\\%P,1);}\n    for\
each my $s (keys(%P))\n      {\n	$nseq++;\n	$strin\
g= $P{$s}{seq};\n	$string = uc $string;\n	my @r=sp\
lit (//,$string);\n	for (my $a=0; $a<=$#r; $a++)\n\
	  {\n	    if (($r[$a]=~/[OHICE]/)){$C{$a}{$r[$a]}\
++;}\n	  }\n      }\n    @l=keys(%C);\n    \n    \\
n    $s=$P{$seq}{seq};\n    print $f \">$seq pred \
based on $nseq\\n\";\n    @r=split (//,$s);\n    \\
n    for (my $x=0; $x<=$#r; $x++)\n      {\n	if ($\
r[$x] ne \"-\")\n	  {\n	    my $h=$C{$x}{H};\n	   \
 my $i=$C{$x}{I};\n	    my $o=$C{$x}{O};\n	    my \
$c=$C{$x}{C};\n	    my $e=$C{$x}{E};\n	    my $l=$\
i+$o;\n	    \n	    if ($h>=$i && $h>=$o && $h>=$c \
&& $h>=$e){$r[$x]='H';}\n	    elsif ($i>=$o && $i>\
=$c && $i>=$e){$r[$x]='I';}\n	    elsif ($o>=$c &&\
 $o>=$e){$r[$x]='O';}\n	    elsif ($c>=$e){$r[$x]=\
'C';}\n	    else {$r[$x]='E';}\n	  }\n      }\n   \
 $j=join ('', @r);\n    print $f \"$j\\n\";\n    c\
lose ($f);\n    return $j;\n  }\n\nsub pred2aln\n \
 {\n    my $PR=shift;\n    my $AR=shift;\n    \n  \
  my $f=new FileHandle;\n    my %P=%$PR;\n    my %\
A=%$AR;\n    my %PA;\n    my $tmp=vtmpnam();\n    \
my $f= new FileHandle;\n    \n    open ($f, \">$tm\
p\");\n    foreach my $s (sort{$A{$a}{index}<=>$A{\
$b}{index}}(keys (%A)))\n      {\n	my (@list, $seq\
, @plist, @pseq, $L, $PL, $c, $w);\n	my $seq;\n	my\
 $seq=$A{$s}{seq};\n	my $pred=$P{$s}{seq};\n	$seq=\
pred2alnS($P{$s}{seq},$A{$s}{seq});\n	print $f \">\
$s\\n$seq\\n\";\n      }\n    close ($f);\n    ret\
urn read_fasta_seq ($tmp);\n  }\nsub pred2alnS\n  \
{\n    my $pred=shift;\n    my $aln= shift;\n    m\
y ($j,$a,$b);\n    my @P=split (//, $pred);\n    m\
y @A=split (//, $aln);\n    for ($a=$b=0;$a<=$#A; \
$a++)\n      {\n	if ($A[$a] ne \"-\"){$A[$a]=$P[$b\
++];}\n      }\n    if ($b!= ($#P+1)){add_warning \
(\"Could not thread sequence: $b $#P\");}\n    \n \
   $j= join ('', @A);\n    return $j;\n  }\nsub pr\
ed2color\n  {\n    my $predP=shift;\n    my $alnP=\
shift;\n    my $out=shift;\n    my $F=new FileHand\
le;\n    my $struc=vtmpnam();\n    my $aln=vtmpnam\
();\n    \n\n    output_fasta_seq ($alnP, $aln);\n\
    my %p=%$predP;\n    \n    open ($F, \">$struc\\
");\n    \n    \n    foreach my $s (keys(%p))\n   \
   {\n	\n	print $F \">$s\\n\";\n	my $s=uc($p{$s}{s\
eq});\n	\n	$s=~s/[Oo]/0/g;\n	$s=~s/[Ee]/0/g;\n	\n	\
$s=~s/[Ii]/5/g;\n	$s=~s/[Cc]/5/g;\n	\n	$s=~s/[Hh]/\
9/g;\n	\n	print $F \"$s\\n\";\n      }\n    close \
($F);\n    \n    \n    \n    safe_system ( \"t_cof\
fee -other_pg seq_reformat -in $aln -struc_in $str\
uc -struc_in_f number_fasta -output color_html -ou\
t $out\");\n    return;\n  }\n	  \n    \nsub displ\
ay_fasta_seq\n  {\n    my $SR=shift;\n    my %S=%$\
SR;\n    \n    foreach my $s (sort{$S{$a}{index}<=\
>$S{$b}{index}}(keys (%S)))\n      {\n	print STDER\
R \">$s\\n$S{$s}{seq}\\n\";\n      }\n    close ($\
f);\n  }\nsub output_fasta_seq\n  {\n    my $SR=sh\
ift;\n    my $outfile=shift;\n    my $mode =shift;\
\n    my $f= new FileHandle;\n    my %S=%$SR;\n   \
 \n    \n    open ($f, \">$outfile\");\n    foreac\
h my $s (sort{$S{$a}{index}<=>$S{$b}{index}}(keys \
(%S)))\n      {\n	my $seq=$S{$s}{seq};\n	if ( $mod\
e eq \"seq\"){$seq=~s/\\-//g;}\n	print $f \">$s\\n\
$seq\\n\";\n      }\n    close ($f);\n  }\n      \\
nsub read_fasta_seq \n  {\n    my $f=$_[0];\n    m\
y %hseq;\n    my (@seq, @com, @name);\n    my ($a,\
 $s,$nseq);\n    my $index;\n    open (F, $f);\n  \
  while (<F>)\n      {\n	$s.=$_;\n      }\n    clo\
se (F);\n\n    \n    @name=($s=~/>(\\S*).*\\n[^>]*\
/g);\n    \n    @seq =($s=~/>.*.*\\n([^>]*)/g);\n \
   @com =($s=~/>.*(.*)\\n([^>]*)/g);\n\n\n    $nse\
q=$#name+1;\n    \n  \n    for ($a=0; $a<$nseq; $a\
++)\n      {\n	my $n=$name[$a];\n	my $s;\n	$hseq{$\
n}{name}=$n;\n	$s=$seq[$a];$s=~s/\\s//g;\n	$hseq{$\
n}{index}=++$index;\n	$hseq{$n}{seq}=$s;\n	$hseq{$\
n}{com}=$com[$a];\n      }\n    return %hseq;\n  }\
\n\n\nsub file2head\n      {\n	my $file = shift;\n\
	my $size = shift;\n	my $f= new FileHandle;\n	my $\
line;\n	open ($f,$file);\n	read ($f,$line, $size);\
\n	close ($f);\n	return $line;\n      }\nsub file2\
tail\n      {\n	my $file = shift;\n	my $size = shi\
ft;\n	my $f= new FileHandle;\n	my $line;\n	\n	open\
 ($f,$file);\n	seek ($f,$size*-1, 2);\n	read ($f,$\
line, $size);\n	close ($f);\n	return $line;\n     \
 }\n\n\nsub vtmpnam\n      {\n	my $r=rand(100000);\
\n	my $f=\"file.$r.$$\";\n	while (-e $f)\n	  {\n	 \
   $f=vtmpnam();\n	  }\n	push (@TMPFILE_LIST, $f);\
\n	return $f;\n      }\n\nsub myexit\n  {\n    my \
$code=@_[0];\n    if ($CLEAN_EXIT_STARTED==1){retu\
rn;}\n    else {$CLEAN_EXIT_STARTED=1;}\n    ### O\
NLY BARE EXIT\n    exit ($code);\n  }\nsub set_err\
or_lock\n    {\n      my $name = shift;\n      my \
$pid=$$;\n\n      \n      &lock4tc ($$,\"LERROR\",\
 \"LSET\", \"$$ -- ERROR: $name $PROGRAM\\n\");\n \
     return;\n    }\nsub set_lock\n  {\n    my $pi\
d=shift;\n    my $msg= shift;\n    my $p=getppid()\
;\n    &lock4tc ($pid,\"LLOCK\",\"LRESET\",\"$p$ms\
g\\n\");\n  }\nsub unset_lock\n   {\n     \n    my\
 $pid=shift;\n    &lock4tc ($pid,\"LLOCK\",\"LRELE\
ASE\",\"\");\n  }\nsub shift_lock\n  {\n    my $fr\
om=shift;\n    my $to=shift;\n    my $from_type=sh\
ift;\n    my $to_type=shift;\n    my $action=shift\
;\n    my $msg;\n    \n    if (!&lock4tc($from, $f\
rom_type, \"LCHECK\", \"\")){return 0;}\n    $msg=\
&lock4tc ($from, $from_type, \"LREAD\", \"\");\n  \
  &lock4tc ($from, $from_type,\"LRELEASE\", $msg);\
\n    &lock4tc ($to, $to_type, $action, $msg);\n  \
  return;\n  }\nsub isshellpid\n  {\n    my $p=shi\
ft;\n    if (!lock4tc ($p, \"LLOCK\", \"LCHECK\"))\
{return 0;}\n    else\n      {\n	my $c=lock4tc($p,\
 \"LLOCK\", \"LREAD\");\n	if ( $c=~/-SHELL-/){retu\
rn 1;}\n      }\n    return 0;\n  }\nsub isrootpid\
\n  {\n    if(lock4tc (getppid(), \"LLOCK\", \"LCH\
ECK\")){return 0;}\n    else {return 1;}\n  }\nsub\
 lock4tc\n	{\n	  my ($pid,$type,$action,$value)=@_\
;\n	  my $fname;\n	  my $host=hostname;\n	  \n	  i\
f ($type eq \"LLOCK\"){$fname=\"$LOCKDIR/.$pid.$ho\
st.lock4tcoffee\";}\n	  elsif ( $type eq \"LERROR\\
"){ $fname=\"$LOCKDIR/.$pid.$host.error4tcoffee\";\
}\n	  elsif ( $type eq \"LWARNING\"){ $fname=\"$LO\
CKDIR/.$pid.$host.warning4tcoffee\";}\n	  \n	  if \
($debug_lock)\n	    {\n	      print STDERR \"\\n\\\
t---lock4tc(tcg): $action => $fname =>$value (RD: \
$LOCKDIR)\\n\";\n	    }\n\n	  if    ($action eq \"\
LCHECK\") {return -e $fname;}\n	  elsif ($action e\
q \"LREAD\"){return file2string($fname);}\n	  elsi\
f ($action eq \"LSET\") {return string2file ($valu\
e, $fname, \">>\");}\n	  elsif ($action eq \"LRESE\
T\") {return string2file ($value, $fname, \">\");}\
\n	  elsif ($action eq \"LRELEASE\") \n	    {\n	  \
    if ( $debug_lock)\n		{\n		  my $g=new FileHand\
le;\n		  open ($g, \">>$fname\");\n		  print $g \"\
\\nDestroyed by $$\\n\";\n		  close ($g);\n		  saf\
e_system (\"mv $fname $fname.old\");\n		}\n	      \
else\n		{\n		  unlink ($fname);\n		}\n	    }\n	  r\
eturn \"\";\n	}\n	\nsub file2string\n	{\n	  my $fi\
le=@_[0];\n	  my $f=new FileHandle;\n	  my $r;\n	 \
 open ($f, \"$file\");\n	  while (<$f>){$r.=$_;}\n\
	  close ($f);\n	  return $r;\n	}\nsub string2file\
 \n    {\n    my ($s,$file,$mode)=@_;\n    my $f=n\
ew FileHandle;\n    \n    open ($f, \"$mode$file\"\
);\n    print $f  \"$s\";\n    close ($f);\n  }\n\\
nBEGIN\n    {\n      srand;\n    \n      $SIG{'SIG\
UP'}='signal_cleanup';\n      $SIG{'SIGINT'}='sign\
al_cleanup';\n      $SIG{'SIGQUIT'}='signal_cleanu\
p';\n      $SIG{'SIGILL'}='signal_cleanup';\n     \
 $SIG{'SIGTRAP'}='signal_cleanup';\n      $SIG{'SI\
GABRT'}='signal_cleanup';\n      $SIG{'SIGEMT'}='s\
ignal_cleanup';\n      $SIG{'SIGFPE'}='signal_clea\
nup';\n      \n      $SIG{'SIGKILL'}='signal_clean\
up';\n      $SIG{'SIGPIPE'}='signal_cleanup';\n   \
   $SIG{'SIGSTOP'}='signal_cleanup';\n      $SIG{'\
SIGTTIN'}='signal_cleanup';\n      $SIG{'SIGXFSZ'}\
='signal_cleanup';\n      $SIG{'SIGINFO'}='signal_\
cleanup';\n      \n      $SIG{'SIGBUS'}='signal_cl\
eanup';\n      $SIG{'SIGALRM'}='signal_cleanup';\n\
      $SIG{'SIGTSTP'}='signal_cleanup';\n      $SI\
G{'SIGTTOU'}='signal_cleanup';\n      $SIG{'SIGVTA\
LRM'}='signal_cleanup';\n      $SIG{'SIGUSR1'}='si\
gnal_cleanup';\n\n\n      $SIG{'SIGSEGV'}='signal_\
cleanup';\n      $SIG{'SIGTERM'}='signal_cleanup';\
\n      $SIG{'SIGCONT'}='signal_cleanup';\n      $\
SIG{'SIGIO'}='signal_cleanup';\n      $SIG{'SIGPRO\
F'}='signal_cleanup';\n      $SIG{'SIGUSR2'}='sign\
al_cleanup';\n\n      $SIG{'SIGSYS'}='signal_clean\
up';\n      $SIG{'SIGURG'}='signal_cleanup';\n    \
  $SIG{'SIGCHLD'}='signal_cleanup';\n      $SIG{'S\
IGXCPU'}='signal_cleanup';\n      $SIG{'SIGWINCH'}\
='signal_cleanup';\n      \n      $SIG{'INT'}='sig\
nal_cleanup';\n      $SIG{'TERM'}='signal_cleanup'\
;\n      $SIG{'KILL'}='signal_cleanup';\n      $SI\
G{'QUIT'}='signal_cleanup';\n      \n      our $de\
bug_lock=$ENV{\"DEBUG_LOCK\"};\n      \n      \n  \
    \n      \n      foreach my $a (@ARGV){$CL.=\" \
$a\";}\n      if ( $debug_lock ){print STDERR \"\\\
n\\n\\n********** START PG: $PROGRAM *************\
\\n\";}\n      if ( $debug_lock ){print STDERR \"\\
\n\\n\\n**********(tcg) LOCKDIR: $LOCKDIR $$ *****\
********\\n\";}\n      if ( $debug_lock ){print ST\
DERR \"\\n --- $$ -- $CL\\n\";}\n      \n	     \n \
     \n      \n    }\nsub flush_error\n  {\n    my\
 $msg=shift;\n    return add_error ($EXIT_FAILURE,\
$$, $$,getppid(), $msg, $CL);\n  }\nsub add_error \
\n  {\n    my $code=shift;\n    my $rpid=shift;\n \
   my $pid=shift;\n    my $ppid=shift;\n    my $ty\
pe=shift;\n    my $com=shift;\n    \n    $ERROR_DO\
NE=1;\n    lock4tc ($rpid, \"LERROR\",\"LSET\",\"$\
pid -- ERROR: $type\\n\");\n    lock4tc ($$, \"LER\
ROR\",\"LSET\", \"$pid -- COM: $com\\n\");\n    lo\
ck4tc ($$, \"LERROR\",\"LSET\", \"$pid -- STACK: $\
ppid -> $pid\\n\");\n   \n    return $code;\n  }\n\
sub add_warning \n  {\n    my $rpid=shift;\n    my\
 $pid =shift;\n    my $command=shift;\n    my $msg\
=\"$$ -- WARNING: $command\\n\";\n    print STDERR\
 \"$msg\";\n    lock4tc ($$, \"LWARNING\", \"LSET\\
", $msg);\n  }\n\nsub signal_cleanup\n  {\n    pri\
nt dtderr \"\\n**** $$ (tcg) was killed\\n\";\n   \
 &cleanup;\n    exit ($EXIT_FAILURE);\n  }\nsub cl\
ean_dir\n  {\n    my $dir=@_[0];\n    if ( !-d $di\
r){return ;}\n    elsif (!($dir=~/tmp/)){return ;}\
#safety check 1\n    elsif (($dir=~/\\*/)){return \
;}#safety check 2\n    else\n      {\n	`rm -rf $di\
r`;\n      }\n    return;\n  }\nsub cleanup\n  {\n\
    #print stderr \"\\n----tc: $$ Kills $PIDCHILD\\
\n\";\n    #kill (SIGTERM,$PIDCHILD);\n    my $p=g\
etppid();\n    $CLEAN_EXIT_STARTED=1;\n    \n    \\
n    \n    if (&lock4tc($$,\"LERROR\", \"LCHECK\",\
 \"\"))\n      {\n	my $ppid=getppid();\n	if (!$ERR\
OR_DONE) \n	  {\n	    &lock4tc($$,\"LERROR\", \"LS\
ET\", \"$$ -- STACK: $p -> $$\\n\");\n	    &lock4t\
c($$,\"LERROR\", \"LSET\", \"$$ -- COM: $CL\\n\");\
\n	  }\n      }\n    my $warning=&lock4tc($$, \"LW\
ARNING\", \"LREAD\", \"\");\n    my $error=&lock4t\
c($$,  \"LERROR\", \"LREAD\", \"\");\n    #release\
 error and warning lock if root\n    \n    if (isr\
ootpid() && ($warning || $error) )\n      {\n	\n	p\
rint STDERR \"**************** Summary ***********\
**\\n$error\\n$warning\\n\";\n\n	&lock4tc($$,\"LER\
ROR\",\"RELEASE\",\"\");\n	&lock4tc($$,\"LWARNING\\
",\"RELEASE\",\"\");\n      } \n    \n    \n    fo\
reach my $f (@TMPFILE_LIST)\n      {\n	if (-e $f){\
unlink ($f);} \n      }\n    foreach my $d (@TMPDI\
R_LIST)\n      {\n	clean_dir ($d);\n      }\n    #\
No More Lock Release\n    #&lock4tc($$,\"LLOCK\",\\
"LRELEASE\",\"\"); #release lock \n\n    if ( $deb\
ug_lock ){print STDERR \"\\n\\n\\n********** END P\
G: $PROGRAM ($$) *************\\n\";}\n    if ( $d\
ebug_lock ){print STDERR \"\\n\\n\\n**********(tcg\
) LOCKDIR: $LOCKDIR $$ *************\\n\";}\n  }\n\
END \n  {\n    \n    &cleanup();\n  }\n   \n\nsub \
safe_system \n{\n  my $com=shift;\n  my $ntry=shif\
t;\n  my $ctry=shift;\n  my $pid;\n  my $status;\n\
  my $ppid=getppid();\n  if ($com eq \"\"){return \
1;}\n  \n  \n\n  if (($pid = fork ()) < 0){return \
(-1);}\n  if ($pid == 0)\n    {\n      set_lock($$\
, \" -SHELL- $com (tcg)\");\n      exec ($com);\n \
   }\n  else\n    {\n      lock4tc ($$, \"LLOCK\",\
 \"LSET\", \"$pid\\n\");#update parent\n      $PID\
CHILD=$pid;\n    }\n  if ($debug_lock){printf STDE\
RR \"\\n\\t .... safe_system (fasta_seq2hmm)  p: $\
$ c: $pid COM: $com\\n\";}\n\n  waitpid ($pid,WTER\
MSIG);\n\n  shift_lock ($pid,$$, \"LWARNING\",\"LW\
ARNING\", \"LSET\");\n\n  if ($? == $EXIT_FAILURE \
|| lock4tc($pid, \"LERROR\", \"LCHECK\", \"\"))\n \
   {\n      if ($ntry && $ctry <$ntry)\n	{\n	  add\
_warning ($$,$$,\"$com failed [retry: $ctry]\");\n\
	  lock4tc ($pid, \"LRELEASE\", \"LERROR\", \"\");\
\n	  return safe_system ($com, $ntry, ++$ctry);\n	\
}\n      elsif ($ntry == -1)\n	{\n	  if (!shift_lo\
ck ($pid, $$, \"LERROR\", \"LWARNING\", \"LSET\"))\
\n	    {\n	      add_warning ($$,$$,\"$com failed\\
");\n	    }\n	  else\n	    {\n	      lock4tc ($pid\
, \"LRELEASE\", \"LERROR\", \"\");\n	    }\n	  ret\
urn $?;}\n      else\n	{\n	  if (!shift_lock ($pid\
,$$, \"LERROR\",\"LERROR\", \"LSET\"))\n	    {\n	 \
     myexit(add_error ($EXIT_FAILURE,$$,$pid,getpp\
id(), \"UNSPECIFIED system\", $com));\n	    }\n	}\\
n    }\n  return $?;\n}\n\nsub check_configuration\
 \n    {\n      my @l=@_;\n      my $v;\n      for\
each my $p (@l)\n	{\n	  \n	  if   ( $p eq \"EMAIL\\
")\n	    { \n	      if ( !($EMAIL=~/@/))\n		{\n		a\
dd_warning($$,$$,\"Could Not Use EMAIL\");\n		myex\
it(add_error ($EXIT_FAILURE,$$,$$,getppid(),\"EMAI\
L\",\"$CL\"));\n	      }\n	    }\n	  elsif( $p eq \
\"INTERNET\")\n	    {\n	      if ( !&check_interne\
t_connection())\n		{\n		  myexit(add_error ($EXIT_\
FAILURE,$$,$$,getppid(),\"INTERNET\",\"$CL\"));\n	\
	}\n	    }\n	  elsif( $p eq \"wget\")\n	    {\n	  \
    if (!&pg_is_installed (\"wget\") && !&pg_is_in\
stalled (\"curl\"))\n		{\n		  myexit(add_error ($E\
XIT_FAILURE,$$,$$,getppid(),\"PG_NOT_INSTALLED:wge\
t\",\"$CL\"));\n		}\n	    }\n	  elsif( !(&pg_is_in\
stalled ($p)))\n	    {\n	      myexit(add_error ($\
EXIT_FAILURE,$$,$$,getppid(),\"PG_NOT_INSTALLED:$p\
\",\"$CL\"));\n	    }\n	}\n      return 1;\n    }\\
nsub pg_is_installed\n  {\n    my @ml=@_;\n    my \
$r, $p, $m;\n    my $supported=0;\n    \n    my $p\
=shift (@ml);\n    if ($p=~/::/)\n      {\n	if (sa\
fe_system (\"perl -M$p -e 1\")==$EXIT_SUCCESS){ret\
urn 1;}\n	else {return 0;}\n      }\n    else\n   \
   {\n	$r=`which $p 2>/dev/null`;\n	if ($r eq \"\"\
){return 0;}\n	else {return 1;}\n      }\n  }\n\n\\
n\nsub check_internet_connection\n  {\n    my $int\
ernet;\n    my $tmp;\n    &check_configuration ( \\
"wget\"); \n    \n    $tmp=&vtmpnam ();\n    \n   \
 if     (&pg_is_installed    (\"wget\")){`wget www\
.google.com -O$tmp >/dev/null 2>/dev/null`;}\n    \
elsif  (&pg_is_installed    (\"curl\")){`curl www.\
google.com -o$tmp >/dev/null 2>/dev/null`;}\n    \\
n    if ( !-e $tmp || -s $tmp < 10){$internet=0;}\\
n    else {$internet=1;}\n    if (-e $tmp){unlink \
$tmp;}\n\n    return $internet;\n  }\nsub check_pg\
_is_installed\n  {\n    my @ml=@_;\n    my $r=&pg_\
is_installed (@ml);\n    if (!$r && $p=~/::/)\n   \
   {\n	print STDERR \"\\nYou Must Install the perl\
 package $p on your system.\\nRUN:\\n\\tsudo perl \
-MCPAN -e 'install $pg'\\n\";\n      }\n    elsif \
(!$r)\n      {\n	myexit(flush_error(\"\\nProgram $\
p Supported but Not Installed on your system\"));\\
n      }\n    else\n      {\n	return 1;\n      }\n\
  }\n\n\n\n","\n\n\n\n\nmy $FMODEL =\"\"; \nmy $TM\
PDIR = \"/tmp\";\n\n\n\n\nmy $NUCALPH = \"ACGTUNRY\
MKSWHBVD\";\nmy $PRIMNUCALPH = \"ACGTUN\";\nuse va\
rs qw($NUCALPH $PRIMNUCALPH $TMPDIR);\n\n\nmy $err\
msg;\nuse vars qw($errmsg);\n\n\n\nuse Getopt::Lon\
g;\nuse Cwd;\nuse File::Basename;\nuse File::Temp \
qw/ tempfile tempdir /;\nuse File::Copy;\nuse File\
::Path;\n\n\n\nsub usage(;$)\n{\n    my ($errmsg) \
= @_;\n    my $myname = basename($0);\n\n    if ($\
errmsg) {\n        print STDERR \"ERROR: $errmsg\\\
n\";\n    }\n\n    print STDERR << \"EOF\";\n    \\
n$myname: align two sequences by means of consan\\\
's sfold\nUsage:\n $myname -i file -o file -d path\
\nOptions:\n -i|--in : pairwise input sequence fil\
e\n -o|--out: output alignment\n -d|--directory co\
ntaining data\n\nEOF\n}\n\nsub read_stk_aln \n  {\\
n    my $f=$_[0];\n    my ($seq, $id);\n    \n    \
my %hseq;\n\n    open (STK, \"$f\");\n    while (<\
STK>)\n      {\n	if ( /^#/ || /^\\/\\// || /^\\s*$\
/){;}\n	else\n	  {\n	    ($id,$seq)=/(\\S+)\\s+(\\\
S+)/;\n	    $hseq{$id}{'seq'}.=$seq;\n	  }\n      \
}\n    close (STK);\n    return %hseq;\n  }\nsub r\
ead_fasta_seq \n  {\n    my $f=$_[0];\n    my %hse\
q;\n    my (@seq, @com, @name);\n    my ($a, $s,$n\
seq);\n\n    open (F, $f);\n    while (<F>)\n     \
 {\n	$s.=$_;\n      }\n    close (F);\n\n    \n   \
 @name=($s=~/>(.*).*\\n[^>]*/g);\n    \n    @seq =\
($s=~/>.*.*\\n([^>]*)/g);\n    @com =($s=~/>.*(.*)\
\\n([^>]*)/g);\n\n    \n    $nseq=$#name+1;\n    \\
n    for ($a=0; $a<$nseq; $a++)\n      {\n	my $n=$\
name[$a];\n	$hseq{$n}{name}=$n;\n	$hseq{$n}{seq}=$\
seq[$a];\n	$hseq{$n}{com}=$com[$a];\n      }\n    \
return %hseq;\n  }\n\n\n\nsub sfold_parseoutput($$\
)\n{\n    my ($frawout, $foutfa) = @_;\n    my %ha\
ln;\n    my ($fstk, $cmd, $id);\n    open FOUTFA, \
\">$foutfa\";\n    \n    $fstk = $frawout . \".stk\
\";\n    \n    # first line of raw out contains in\
fo\n    # remaining stuff is stockholm formatted\n\
    $cmd = \"sed -e '1d' $frawout\";\n    system(\\
"$cmd > $fstk\");\n    if ($? != 0) {\n        $er\
rmsg = \"command failed with exit status $?.\";\n \
       $errmsg .=  \"Command was \\\"$cmd\\\"\";\n\
        return -1;\n    }\n\n    # this gives an e\
rror message. just ignore it...\n    %haln=read_st\
k_aln ( $fstk);\n    foreach $i (keys (%haln))\n  \
    {\n	my $s;\n	$s=$haln{$i}{'seq'};\n	$s =~ s/\\\
./-/g;\n	print FOUTFA \">$i\\n$s\\n\";\n      }\n \
   close FOUTFA;\n    return 0;\n}\n\n\n\n\nsub sf\
old_wrapper($$$$)\n{\n    \n    my ($fs1, $fs2, $f\
model, $foutfa) = @_;\n    \n\n    my ($cmd, $fraw\
out, $ferrlog, $freadme, $ftimelog, $fstk);\n\n   \
 # add  basename($fmsqin) (unknown here!)\n    $fr\
awout = \"sfold.log\";\n    $ferrlog = \"sfold.err\
\";\n    $ftimelog = \"sfold.time\";\n    $freadme\
 =  \"sfold.README\";\n    $fstk = \"sfold.stk\";\\
n    \n    # prepare execution...\n    #\n    # ./\
tmp is essential for dswpalign\n    # otherwise yo\
u'll get a segfault\n    mkdir \"./tmp\";\n    \n \
   $cmd = \"sfold -m $fmodel $fs1 $fs2\";\n    ope\
n(FREADME,\">$freadme\");\n    print FREADME \"$cm\
d\\n\"; \n    close(FREADME);\n\n    # and go\n   \
 #\n    system(\"/usr/bin/time -p -o $ftimelog $cm\
d >$frawout 2>$ferrlog\");\n    if ($? != 0) {\n  \
      $errmsg = \"command failed with exit status \
$?\";\n        $errmsg .= \"command was \\\"$cmd\\\
\". See \" . getcwd . \"\\n\";\n        return -1;\
\n    }\n\n    return sfold_parseoutput($frawout, \
$foutfa);\n}\n\n\n\n\n\n\n\nmy ($help, $fmsqin, $f\
msaout);\nGetOptions(\"help\"  => \\$help,\n      \
     \"in=s\" => \\$fmsqin,\n           \"out=s\" \
=> \\$fmsaout,\n	   \"data=s\" => \\$ref_dir);\n\n\
\n\nif ($help) {\n    usage();\n    exit(0);\n}\ni\
f (! defined($fmsqin)) {\n    usage('missing input\
 filename');\n    exit(1);\n}\nif (! defined($fmsa\
out)) {\n    usage('missing output filename');\n  \
  exit(1);\n\n}\nif (scalar(@ARGV)) {\n    usage('\
Unknown remaining args');\n    exit(1);\n}\n\n$FMO\
DEL = \"$ref_dir/mix80.mod\";\nif (! -e \"$FMODEL\\
") {\n    die(\"couldn't find sfold grammar model \
file. Expected $FMODEL\\n\");\n}\n\n\nmy %hseq=rea\
d_fasta_seq ($fmsqin);\nmy $id;\n\nforeach $id (ke\
ys(%hseq))\n  {\n    push(@seq_array, $hseq{$id});\
\n  }\n\nif ( scalar(@seq_array) != 2 ) {\n    die\
(\"Need *exactly* two sequences as input (pairwise\
 alignment!).\")\n}\n\n\n\nmy ($sec, $min, $hour, \
$mday, $mon, $year, $wday, $yday, $isdst) = localt\
ime(time);\nmy $datei = sprintf(\"%4d-%02d-%02d\",\
 $year+1900, $mon+1, $mday);\nmy $templ = basename\
($0) . \".\" . $datei . \".pid-\" . $$ . \".XXXXXX\
\";\nmy $wd = tempdir ( $templ, DIR => $TMPDIR);\n\
\ncopy($fmsqin, \"$wd/\" . basename($fmsqin) . \".\
org\"); # for reproduction\ncopy($FMODEL, \"$wd\")\
;\nmy $fmodel = basename($FMODEL);\nmy $orgwd = ge\
tcwd;\nchdir $wd;\n\n\n\nmy @sepseqfiles;\nforeach\
 $id (keys(%hseq)) {\n    my ($seq, $orgseq, $fnam\
e, $sout);\n    $seq=$hseq{$id}{'seq'};\n    \n   \
 $fname = basename($fmsqin) . \"_$id.fa\";\n    # \
replace funnies in file/id name (e.g. \"/\" \" \" \
etc)\n    $fname =~ s,[/ ],_,g;\n    open (PF, \">\
$fname\");\n    print (PF \">$id\\n$seq\\n\");\n  \
  close (PF);\n\n    push(@sepseqfiles, $fname);\n\
}\n\nmy ($f1, $f2, $fout);\n$f1 = $sepseqfiles[0];\
\n$f2 = $sepseqfiles[1];\n$fout = $wd . basename($\
fmsqin) . \".out.fa\";\nif (sfold_wrapper($f1, $f2\
, $fmodel, \"$fout\") != 0) {\n    printf STDERR \\
"ERROR: See logs in $wd\\n\";\n    exit(1);\n} els\
e {\n    chdir $orgwd;\n    copy($fout, $fmsaout);\
\n    rmtree($wd);\n   exit(0);\n}\n","\nuse Env q\
w(HOST);\nuse Env qw(HOME);\nuse Env qw(USER);\n\n\
\n$tmp=clean_cr ($ARGV[0]);\nopen (F, $tmp);\n\nwh\
ile ( <F>)\n  {\n    my $l=$_;\n    if ( $l=~/^# S\
TOCKHOLM/){$stockholm=1;}\n    elsif ( $stockholm \
&& $l=~/^#/)\n      {\n	$l=~/^#(\\S+)\\s+(\\S+)\\s\
+(\\S*)/g;\n	$l=\"_stockholmhasch_$1\\_stockholmsp\
ace_$2 $3\\n\";\n      }\n    $file.=$l;\n  }\nclo\
se (F);\nunlink($tmp);\n$file1=$file;\n\n$file=~s/\
\\#/_hash_symbol_/g;\n$file=~s/\\@/_arobase_symbol\
_/g;\n\n\n$file=~s/\\n[\\.:*\\s]+\\n/\\n\\n/g;\n\n\
$file=~s/\\n[ \\t\\r\\f]+(\\b)/\\n\\1/g;\n\n\n$fil\
e=~s/(\\n\\S+)(\\s+)(\\S)/\\1_blank_\\3/g;\n\n$fil\
e=~s/[ ]//g;\n$file=~s/_blank_/ /g;\n\n\n\n$file =\
~s/\\n\\s*\\n/#/g;\n\n$file.=\"#\";\n$file =~s/\\n\
/@/g;\n\n\n\n\n@blocks=split /\\#/, $file;\nshift \
(@blocks);\n@s=split /\\@/, $blocks[0];\n$nseq=$#s\
+1;\n\n\n\n$file=join '@', @blocks;\n@lines=split \
/\\@/,$file;\n\n$c=0;\n\nforeach $l (@lines)\n  {\\
n    if (!($l=~/\\S/)){next;}\n    elsif ($stockho\
lm && ($l=~/^\\/\\// || $l=~/STOCKHOLM/)){next;}#g\
et read of STOCHOLM Terminator\n   \n    $l=~/(\\S\
+)\\s+(\\S*)/g;\n    $n=$1; $s=$2;\n    \n    $seq\
[$c].=$s;\n    $name[$c]=$n;\n    $c++;\n    \n   \
 if ( $c==$nseq){$c=0;}\n    \n  } \n\nif ( $c!=0)\
\n      {\n	print STDERR \"ERROR: $ARGV[0] is NOT \
an MSA in Clustalw format: make sure there is no b\
lank line within a block [ERROR]\\n\";\n	exit (EXI\
T_FAILURE);\n      }\n\nfor ($a=0; $a< $nseq; $a++\
)\n  {\n    $name[$a]=cleanstring ($name[$a]);\n  \
  $seq[$a]=cleanstring ($seq[$a]);\n    $seq[$a]=b\
reakstring($seq[$a], 60);\n    \n    $line=\">$nam\
e[$a]\\n$seq[$a]\\n\";\n    \n    print \"$line\";\
\n  }\nexit (EXIT_SUCCESS);\n\nsub cleanstring\n  \
{\n    my $s=@_[0];\n    $s=~s/_hash_symbol_/\\#/g\
;\n    $s=~s/_arobase_symbol_/\\@/g;\n    $s=~s/[ \
\\t]//g;\n    return $s;\n  }\nsub breakstring\n  \
{\n    my $s=@_[0];\n    my $size=@_[1];\n    my @\
list;\n    my $n,$ns, $symbol;\n    \n    @list=sp\
lit //,$s;\n    $n=0;$ns=\"\";\n    foreach $symbo\
l (@list)\n      {\n	if ( $n==$size)\n	  {\n	    $\
ns.=\"\\n\";\n	    $n=0;\n	  }\n	$ns.=$symbol;\n	$\
n++;\n      }\n    return $ns;\n    }\n\nsub clean\
_cr\n  {\n    my $f=@_[0];\n    my $file;\n    \n \
   $tmp=\"f$.$$\";\n    \n    \n    open (IN, $f);\
\n    open (OUT, \">$tmp\");\n    \n    while ( <I\
N>)\n      {\n	$file=$_;\n	$file=~s/\\r\\n/\\n/g;\\
n	$file=~s/\\n\\r/\\n/g;\n	$file=~s/\\r\\r/\\n/g;\\
n	$file=~s/\\r/\\n/g;\n	print OUT \"$file\";\n    \
  }\n    \n    close (IN);\n    close (OUT);\n    \
return $tmp;\n  }\n","use Env qw(HOST);\nuse Env q\
w(HOME);\nuse Env qw(USER);\n\n\n$query_start=-1;\\
n$query_end=-1;\n\nwhile (<>)\n  {\n    if ( /\\/\\
\//){$in_aln=1;}\n    elsif ( $in_aln && /(\\S+)\\\
s+(.*)/)\n      {\n\n\n	$name=$1;\n	\n\n	$seq=$2;\\
n	$seq=~s/\\s//g;\n        $seq=~s/\\~/\\-/g;\n	$s\
eq=~s/\\./\\-/g;\n	if ( $list{$n}{'name'} && $list\
{$n}{'name'} ne $name)\n	  {\n	    print \"$list{$\
n}{'name'} Vs $name\";\n	    \n	    exit (EXIT_FAI\
LURE);\n	  }\n	else\n	  {\n	    $list{$n}{'name'}=\
 $name;\n	  }\n\n	$list{$n}{'seq'}=$list{$n}{'seq'\
}.$seq;\n	\n	$nseq=++$n;\n	\n      }\n    else\n  \
    {$n=0;}\n  }\n\n\nfor ($a=0; $a<$nseq; $a++)\n\
  {\n    print \">$list{$a}{'name'}\\n$list{$a}{'s\
eq'}\\n\";\n  }\n      \n","\nuse Env qw(HOST);\nu\
se Env qw(HOME);\nuse Env qw(USER);\n\n           \
                                             \nuse\
 strict;                                          \
   \nuse warnings;\nuse diagnostics;\n\nmy $in_hit\
_list, my $in_aln=0, my(%name_list)=(),my (%list)=\
(),my $n_seq=0; my $test=0;\nmy($j)=0, my $n=0, my\
 $nom, my $lg_query, my %vu=();\n\nopen (F, \">tmp\
\");\n\n$/=\"\\n\";\nwhile (<>)\n{\n    print F $_\
;\n    if($_ =~ /Query=\\s*(.+?)\\s/i) { $nom=$1;}\
\n\n    if ( /Sequences producing significant alig\
nments/){$in_hit_list=1;}\n    \n    if ($_=~ /^pd\
b\\|/i) { $_=~ s/pdb\\|//g; }\n    if ($_=~ /^(1_\\
\d+)\\s+\\d+/) { $_=~ s/$1/QUERY/;}\n      \n    i\
f ( /^(\\S+).+?\\s+[\\d.]+\\s+([\\de.-]+)\\s+$/ &&\
 $in_hit_list)	\n    {\n	my($id)=$1; # \n	$id=~ s/\
\\|/_/g; #\n	if ($id =~ /.+_$/) { chop($id) }; #\n\
	$name_list{$n_seq++}=$id;\n	$name_list{$n_seq-1}=\
~ s/.*\\|//g;     \n    }\n  \n    if (/query/i) {\
$in_aln=1;}\n    if ( /^(\\S+)\\s+(\\d+)\\s+([a-zA\
-Z-]+)\\s+(\\d+)/ || /^(\\S+)(\\s+)(\\-+)(\\s+)/ &\
& ($in_aln == 1))\n    {\n	my $name=$1;\n	my $star\
t=$2;\n	my $seq=$3;\n	my $end=$4;\n		\n	if ($name \
=~ /QUERY/i) { $lg_query=length($seq); }\n\n	unles\
s ($test > $n) #m\n	{\n	    my(@seqq)= split('',$s\
eq);\n	    my($gap_missing)= scalar(@seqq);\n	    \
\n	    while ($gap_missing != $lg_query)  { unshif\
t (@seqq,\"-\"); $gap_missing= scalar(@seqq); }\n	\
    $seq=join('',@seqq);  #m\n	}\n	\n	if ($name =~\
 /QUERY/i)\n	{\n	    $n=0; %vu=(); $j=0;\n	    $li\
st{$n}{'real_name'}=\"$nom\";\n	}	\n	else\n	{\n	  \
  unless (exists $vu{$name}) { ++$j;}	\n	    $list\
{$n}{'real_name'}=$name_list{$j-1};\n	}\n		\n	$lis\
t{$n}{'name'}=$name;\n\n	$seq=~tr/a-z/A-Z/;\n	$lis\
t{$n}{'seq'}=$list{$n}{'seq'};\n	$list{$n}{'seq'}.\
=$seq;\n\n	$n++;\n	$vu{$name}++;\n	$test++;\n   } \
\n    \n}\n\nmy @numero=();\n\nfor (my $a=0; $a<$n\
; $a++) #m\n{\n    my $long=length($list{0}{'seq'}\
);  \n    my $long1= length($list{$a}{'seq'});\n  \
\n    while ($long1 ne $long)\n    {\n	$list{$a}{'\
seq'}.=\"-\";\n	$long1= length ($list{$a}{'seq'});\
\n    } \n \n    push (@numero,\"$list{$a}{'name'}\
 $list{$a}{'real_name'}\\n\");\n}\n\nmy %dejavu=()\
;\n\n\nfor (my $i=0; $i<=$#numero; $i++)\n{\n    m\
y $s=\">$list{$i}{'real_name'}\\n$list{$i}{'seq'}\\
\n\";\n    my $k=0;\n    \n    if (exists $dejavu{\
$numero[$i]}) {next;}\n    else\n    {	\n	for ($j=\
0; $j<$n ; $j++)\n	{\n	    if (\"$numero[$i]\" eq \
\"$numero[$j]\" && $j != $i )\n	    {\n		++$k;\n		\
$s .=\">$list{$j}{'real_name'}\\n$list{$j}{'seq'}\\
\n\";\n	    }\n	}	\n    }\n    \n    if ($k>0) \n \
   {\n	my $cons;\n	open (SOR,\">tempo_aln2cons\");\
 print SOR $s;  close SOR ;\n	open (COM,\"t_coffee\
 -other_pg seq_reformat -in tempo_aln2cons -action\
 +aln2cons +upper |\") ; \n     	while (<COM>)\n	{\
	\n	    if (/^>/) { $cons =\">$list{$i}{'real_name\
'}\\n\"; next;}\n	    $_=~ s/\\n//g;\n	    $cons .\
=$_;\n	}\n	close COM; unlink (\"tempo_aln2cons\");\
\n	print $cons,\"\\n\"; print F $cons,\"\\n\";\n  \
  }	\n    else  { print $s;  print F $s; }\n    \n\
    $dejavu{$numero[$i]}++;\n} #m\n\nexit;\n\n\n\n\
\n\n\n\n\n\n\n\n","use Env;\n\n\n$tmp_dir=\"\";\n$\
init_dir=\"\";\n$program=\"tc_generic_method.pl\";\
\n\n$blast=@ARGV[0];\n\n$name=\"query\";$seq=\"\";\
\n%p=blast_xml2profile($name,$seq,100, 0, 0, $blas\
t);\n&output_profile (%p);\n\n\nsub output_profile\
\n  {\n    my (%profile)=(@_);\n    my ($a);\n    \
for ($a=0; $a<$profile{n}; $a++)\n      {\n	\n	pri\
nt \">$profile{$a}{name} $profile{$a}{comment}\\n$\
profile{$a}{seq}\\n\";\n      }\n    return;\n  }\\
nsub file_contains \n  {\n    my ($file, $tag, $ma\
x)=(@_);\n    my ($n);\n    $n=0;\n    \n    if ( \
!-e $file && ($file =~/$tag/)) {return 1;}\n    el\
sif ( !-e $file){return 0;}\n    else \n      {\n	\
open (FC, \"$file\");\n	while ( <FC>)\n	  {\n	    \
if ( ($_=~/$tag/))\n	      {\n		close (FC);\n		ret\
urn 1;\n	      }\n	    elsif ($max && $n>$max)\n	 \
     {\n		close (FC);\n		return 0;\n	      }\n	   \
 $n++;\n	  }\n      }\n    close (FC);\n    return\
 0;\n  }\n	    \n	  \nsub file2string\n  {\n    my\
 $f=@_[0];\n    my $string, $l;\n    open (F,\"$f\\
");\n    while (<F>)\n      {\n\n	$l=$_;\n	#chomp \
($l);\n	$string.=$l;\n      }\n    close (F);\n   \
 $string=~s/\\r\\n//g;\n    $string=~s/\\n//g;\n  \
  return $string;\n  }\n\n\n\nsub tag2value \n  {\\
n    \n    my $tag=(@_[0]);\n    my $word=(@_[1]);\
\n    my $return;\n    \n    $tag=~/$word=\"([^\"]\
+)\"/;\n    $return=$1;\n    return $return;\n  }\\
n      \nsub hit_tag2pdbid\n  {\n    my $tag=(@_[0\
]);\n    my $pdbid;\n       \n    $tag=~/id=\"(\\S\
+)\"/;\n    $pdbid=$1;\n    $pdbid=~s/_//;\n    re\
turn $pdbid;\n  }\nsub id2pdbid \n  {\n    my $id=\
@_[0];\n  \n    if ($id =~/pdb/)\n      {\n	$id=~/\
pdb(.*)/;\n	$id=$1;\n      }\n    $id=~s/[|�_]//g;\
\n    return $id;\n  }\nsub set_blast_type \n  {\n\
    my $file =@_[0];\n    if (&file_contains ($fil\
e,\"EBIApplicationResult\",100)){$BLAST_TYPE=\"EBI\
\";}\n    elsif (&file_contains ($file,\"NCBI_Blas\
tOutput\",100)) {$BLAST_TYPE=\"NCBI\";}\n    else\\
n      {\n	$BLAST_TYPE=\"\";\n      }\n    return \
$BLAST_TYPE;\n  }\nsub blast_xml2profile \n  {\n  \
  my ($name,$seq,$maxid, $minid, $mincov, $file)=(\
@_);\n    my (%p, $a, $string, $n);\n    \n\n\n   \
 if ($BLAST_TYPE eq \"EBI\" || &file_contains ($fi\
le,\"EBIApplicationResult\",100)){%p=ebi_blast_xml\
2profile(@_);}\n    elsif ($BLAST_TYPE eq \"NCBI\"\
 || &file_contains ($file,\"NCBI_BlastOutput\",100\
)){%p=ncbi_blast_xml2profile(@_);}\n    else \n   \
   {\n	print \"************ ERROR: Blast Returned \
an unknown XML Format **********************\";\n	\
die;\n      }\n    for ($a=0; $a<$p{n}; $a++)\n   \
   {\n	my $name=$p{$a}{name};\n	$p{$name}{seq}=$p{\
$a}{seq};\n      }\n    return %p;\n  }\nsub ncbi_\
blast_xml2profile \n  {\n    my ($name,$seq,$maxid\
, $minid, $mincov, $string)=(@_);\n    my ($L,$l, \
$a,$b,$c,$d,$nhits,@identifyerL);\n    \n    \n   \
 $seq=~s/[^a-zA-Z]//g;\n    $L=length ($seq);\n   \
 \n    %hit=&xml2tag_list ($string, \"Hit\");\n   \
 \n    \n    for ($nhits=0,$a=0; $a<$hit{n}; $a++)\
\n      {\n	my ($ldb,$id, $identity, $expectation,\
 $start, $end, $coverage, $r);\n	my (%ID,%DE,%HSP)\
;\n	\n	$ldb=\"\";\n\n	%ID=&xml2tag_list ($hit{$a}{\
body}, \"Hit_id\");\n	$identifyer=$ID{0}{body};\n	\
\n	%DE=&xml2tag_list ($hit{$a}{body}, \"Hit_def\")\
;\n	$definition=$DE{0}{body};\n	\n	%HSP=&xml2tag_l\
ist ($hit{$a}{body}, \"Hsp\");\n	for ($b=0; $b<$HS\
P{n}; $b++)\n	  {\n	    my (%START,%END,%E,%I,%Q,%\
M);\n\n	 \n	    %START=&xml2tag_list ($HSP{$b}{bod\
y}, \"Hsp_query-from\");\n	    %HSTART=&xml2tag_li\
st ($HSP{$b}{body}, \"Hsp_hit-from\");\n	    \n	  \
  %LEN=  &xml2tag_list ($HSP{$b}{body}, \"Hsp_alig\
n-len\");\n	    %END=  &xml2tag_list ($HSP{$b}{bod\
y}, \"Hsp_query-to\");\n	    %HEND=  &xml2tag_list\
 ($HSP{$b}{body}, \"Hsp_hit-to\");\n	    %E=&xml2t\
ag_list     ($HSP{$b}{body}, \"Hsp_evalue\");\n	  \
  %I=&xml2tag_list     ($HSP{$b}{body}, \"Hsp_iden\
tity\");\n	    %Q=&xml2tag_list     ($HSP{$b}{body\
}, \"Hsp_qseq\");\n	    %M=&xml2tag_list     ($HSP\
{$b}{body}, \"Hsp_hseq\");\n	    \n	    for ($e=0;\
 $e<$Q{n}; $e++)\n\n	      {\n		$qs=$Q{$e}{body};\\
n		$ms=$M{$e}{body};\n		if ($seq eq\"\"){$seq=$qs;\
$L=length($seq);}\n		\n		$expectation=$E{$e}{body}\
;\n		$identity=($LEN{$e}{body}==0)?0:$I{$e}{body}/\
$LEN{$e}{body}*100;\n		$start=$START{$e}{body};\n	\
	$end=$END{$e}{body};\n		$Hstart=$HSTART{$e}{body}\
;\n		$Hend=$HEND{$e}{body};\n	\n		$coverage=(($end\
-$start)*100)/$L;\n\n	\n		if ($identity>$maxid || \
$identity<$minid || $coverage<$mincov){next;}\n		@\
lr1=(split (//,$qs));\n		@lr2=(split (//,$ms));\n	\
	$l=$#lr1+1;\n		for ($c=0;$c<$L;$c++){$p[$nhits][$\
c]=\"-\";}\n		for ($d=0,$c=0; $c<$l; $c++)\n		  {\\
n		    $r=$lr1[$c];\n		    if ( $r=~/[A-Za-z]/)\n	\
	      {\n			\n			$p[$nhits][$d + $start-1]=$lr2[$\
c];\n			$d++;\n		      }\n		  }\n		$Qseq[$nhits]=$\
qs;\n		$Hseq[$nhits]=$ms;\n		$QstartL[$nhits]=$sta\
rt;\n		$HstartL[$nhits]=$Hstart;\n		$identityL[$nh\
its]=$identity;\n		$endL[$nhits]=$end;\n		$definit\
ionL[$nhits]=$definition;\n		$identifyerL[$nhits]=\
$identifyer;\n		$comment[$nhits]=\"$ldb|$identifye\
r [Eval=$expectation][id=$identity%][start=$Hstart\
 end=$Hend]\";\n		$nhits++;\n	      }\n	  }\n     \
 }\n    \n    $profile{n}=0;\n    $profile{$profil\
e{n}}{name}=$name;\n    $profile{$profile{n}}{seq}\
=$seq;\n    $profile {n}++;\n    \n    for ($a=0; \
$a<$nhits; $a++)\n      {\n	$n=$a+1;\n	\n	$profile\
{$n}{name}=\"$name\\_$a\";\n	$profile{$n}{seq}=\"\\
";\n	$profile{$n}{Qseq}=$Qseq[$a];\n	$profile{$n}{\
Hseq}=$Hseq[$a];\n	$profile{$n}{Qstart}=$QstartL[$\
a];\n	$profile{$n}{Hstart}=$HstartL[$a];\n	$profil\
e{$n}{identity}=$identityL[$a];\n	$profile{$n}{def\
inition}=$definitionL[$a];\n	$profile{$n}{identify\
er}=$identifyerL[$a];\n	$profile{$n}{comment}=$com\
ment[$a];\n	for ($b=0; $b<$L; $b++)\n	  {\n	    if\
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
or ($b=0; $b<$Q{n}; $b++)\n	  {\n	    \n	    \n	  \
  $qs=$Q{$b}{body};\n	    $ms=$M{$b}{body};\n	    \
if ($seq eq\"\"){$seq=$qs;$L=length($seq);}\n\n	  \
  $expectation=$E{$b}{body};\n	    $identity=$I{$b\
}{body};\n	    \n	    	    \n	    $start=&tag2valu\
e ($Q{$b}{open}, \"start\");\n	    $end=&tag2value\
 ($Q{$b}{open}, \"end\");\n	    $startM=&tag2value\
 ($M{$b}{open}, \"start\");\n	    $endM=&tag2value\
 ($M{$b}{open}, \"end\");\n	    $coverage=(($end-$\
start)*100)/$L;\n	    \n	   # print \"$id: ID: $id\
entity COV: $coverage [$start $end]\\n\";\n	    \n\
	    \n	    if ($identity>$maxid || $identity<$min\
id || $coverage<$mincov){next;}\n	    # print \"KE\
EP\\n\";\n\n	    \n	    @lr1=(split (//,$qs));\n	 \
   @lr2=(split (//,$ms));\n	    $l=$#lr1+1;\n	    \
for ($c=0;$c<$L;$c++){$p[$nhits][$c]=\"-\";}\n	   \
 for ($d=0,$c=0; $c<$l; $c++)\n	      {\n		$r=$lr1\
[$c];\n		if ( $r=~/[A-Za-z]/)\n		  {\n		    \n		  \
  $p[$nhits][$d + $start-1]=$lr2[$c];\n		    $d++;\
\n		  }\n	      }\n	  \n	    \n	    $identifyerL[$\
nhits]=$identifyer;\n	    $comment[$nhits]=\"$ldb|\
$identifyer [Eval=$expectation][id=$identity%][sta\
rt=$startM end=$endM]\";\n	    $nhits++;\n	  }\n  \
    }\n    \n    $profile{n}=0;\n    $profile{$pro\
file{n}}{name}=$name;\n    $profile{$profile{n}}{s\
eq}=$seq;\n    $profile {n}++;\n    \n    for ($a=\
0; $a<$nhits; $a++)\n      {\n	$n=$a+1;\n	$profile\
{$n}{name}=\"$name\\_$a\";\n	$profile{$n}{seq}=\"\\
";\n	$profile{$n}{identifyer}=$identifyerL[$a];\n	\
\n	$profile{$n}{comment}=$comment[$a];\n	for ($b=0\
; $b<$L; $b++)\n	  {\n	    if ($p[$a][$b])\n	     \
 {\n		$profile{$n}{seq}.=$p[$a][$b];\n	      }\n	 \
   else\n	      {\n		$profile{$n}{seq}.=\"-\";\n	 \
     }\n	  }\n      }\n    $profile{n}=$nhits+1;\n\
    \n    return %profile;\n  }\n\nsub blast_xml2h\
it_list\n  {\n    my $string=(@_[0]);\n    return \
&xml2tag_list ($string, \"hit\");\n  }\nsub xml2ta\
g_list  \n  {\n    my ($string_in,$tag)=@_;\n    m\
y $tag_in, $tag_out;\n    my %tag;\n    \n    if (\
-e $string_in)\n      {\n	$string=&file2string ($s\
tring_in);\n      }\n    else\n      {\n	$string=$\
string_in;\n      }\n    $tag_in1=\"<$tag \";\n   \
 $tag_in2=\"<$tag>\";\n    $tag_out=\"/$tag>\";\n \
   $string=~s/>/>##1/g;\n    $string=~s/</##2</g;\\
n    $string=~s/##1/<#/g;\n    $string=~s/##2/#>/g\
;\n    @l=($string=~/(\\<[^>]+\\>)/g);\n    $tag{n\
}=0;\n    $in=0;$n=-1;\n  \n \n\n    foreach $t (@\
l)\n      {\n\n	$t=~s/<#//;\n	$t=~s/#>//;\n	\n	if \
( $t=~/$tag_in1/ || $t=~/$tag_in2/)\n	  {\n	 \n	  \
  $in=1;\n	    $tag{$tag{n}}{open}=$t;\n	    $n++;\
\n	    \n	  }\n	elsif ($t=~/$tag_out/)\n	  {\n	   \
 \n\n	    $tag{$tag{n}}{close}=$t;\n	    $tag{n}++\
;\n	    $in=0;\n	  }\n	elsif ($in)\n	  {\n	   \n	 \
   $tag{$tag{n}}{body}.=$t;\n	  }\n      }\n  \n  \
  return %tag;\n  }\n\n\n\n\n","use Env qw(HOST);\\
nuse Env qw(HOME);\nuse Env qw(USER);\nwhile (<>)\\
n  {\n    if ( /^>(\\S+)/)\n      {\n	if ($list{$1\
})\n	  {\n	    print \">$1_$list{$1}\\n\";\n	    $\
list{$1}++;\n	  }\n	else\n	  {\n	    print $_;\n	 \
   $list{$1}=1;\n	  }\n      }\n    else\n      {\\
n	print $_;\n      }\n  }\n      \n","\n\n\nuse En\
v qw(HOST);\nuse Env qw(HOME);\nuse Env qw(USER);\\
n\n\nopen (F,$ARGV[0]);\nwhile ( <>)\n  {\n    @x=\
/([^:,;\\)\\(\\s]+):[^:,;\\)\\(]*/g;\n    @list=(@\
list,@x);\n  }\n$n=$#list+1;\nforeach $n(@list){pr\
int \">$n\\nsequence\\n\";}\n\n\nclose (F);\n","\n\
open (F, $ARGV[0]);\n\nwhile ( <F>)\n  {\n    @l=(\
$_=~/(\\S+)/g);\n    \n    $name=shift @l;\n    \n\
    print STDOUT \"\\n>$name\\n\";\n    foreach $e\
 (@l){$e=($e eq \"0\")?\"O\":\"I\";print \"$e\";}\\
n  }\nclose (F);\n\n		       \n    \n","use Env qw\
(HOST);\nuse Env qw(HOME);\nuse Env qw(USER);\n\n$\
tmp=\"$ARGV[0].$$\";\nopen (IN, $ARGV[0]);\nopen (\
OUT, \">$tmp\");\n\nwhile ( <IN>)\n  {\n    $file=\
$_;\n    $file=~s/\\r\\n/\\n/g;\n    $file=~s/\\n\\
\r/\\n/g;\n    $file=~s/\\r\\r/\\n/g;\n    $file=~\
s/\\r/\\n/g;\n    print OUT \"$file\";\n  }\nclose\
 (IN);\nclose (OUT);\n\nopen (OUT, \">$ARGV[0]\");\
\nopen (IN, \"$tmp\");\n\nwhile ( <IN>)\n{\n  prin\
t OUT \"$_\";\n}\nclose (IN);\nclose (OUT);\nunlin\
k ($tmp);\n\n"};
