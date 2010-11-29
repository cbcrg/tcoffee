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
n\n\nif ($mode eq \"seq_msa\")\n  {\n    &seq2msa(\
$mode,&my_get_opt ( $cl, \"-infile=\",1,1, \"-meth\
od=\",1,2, \"-param=\",0,0,\"-outfile=\",1,0, \"-d\
atabase=\",0,0));\n  }\nelsif ( $mode eq \"tblastx\
_msa\")\n  {\n    &seq2tblastx_lib ($mode,&my_get_\
opt ( $cl, \"-infile=\",1,1, \"-outfile=\",1,0));\\
n  }\nelsif ( $mode eq \"tblastpx_msa\")\n  {\n   \
 &seq2tblastpx_lib ($mode,&my_get_opt ( $cl, \"-in\
file=\",1,1, \"-outfile=\",1,0));\n  }\nelsif ( $m\
ode eq \"thread_pair\")\n  {\n    &seq2thread_pair\
($mode,&my_get_opt ( $cl, \"-infile=\",1,1, \"-pdb\
file1=\",1,1, \"-method=\",1,2,\"-param=\",0,0, \"\
-outfile=\",1,0, ));\n  }\nelsif ( $mode eq \"pdbi\
d_pair\")\n  {\n    &seq2pdbid_pair($mode,&my_get_\
opt ( $cl, \"-pdbfile1=\",1,0, \"-pdbfile2=\",1,0,\
 \"-method=\",1,2,\"-param=\",0,0, \"-outfile=\",1\
,0, ));\n  }\nelsif ( $mode eq \"pdb_pair\")\n  {\\
n    &seq2pdb_pair($mode,&my_get_opt ( $cl, \"-pdb\
file1=\",1,1, \"-pdbfile2=\",1,1, \"-method=\",1,2\
,\"-param=\",0,0, \"-outfile=\",1,0, ));\n  }\nels\
if ( $mode eq \"profile_pair\")\n  {\n     &seq2pr\
ofile_pair($mode,&my_get_opt ( $cl, \"-profile1=\"\
,1,1, \"-profile2=\",1,1, \"-method=\",1,2,\"-para\
m=\",0,0, \"-outfile=\",1,0, ));\n  }\nelsif ( $mo\
de eq \"pdb_template\")\n  {\n    &blast2pdb_templ\
ate ($mode,&my_get_opt ( $cl, \"-infile=\",1,1, \"\
-database=\",1,0, \"-method=\",1,0, \"-outfile=\",\
1,0,\"-pdb_type=\",1,0,));\n  }\nelsif ( $mode eq \
\"profile_template\")\n  {\n    \n    &psiblast2pr\
ofile_template ($mode,&my_get_opt ( $cl, \"-infile\
=\",1,1, \"-database=\",1,0, \"-method=\",1,0, \"-\
outfile=\",1,0));\n  }\nelsif ( $mode eq \"psiprof\
ile_template\")\n  {\n    &psiblast2profile_templa\
te ($mode,&my_get_opt ( $cl, \"-infile=\",1,1, \"-\
database=\",1,0, \"-method=\",1,0, \"-outfile=\",1\
,0));\n  }\nelsif ( $mode eq \"RNA_template\")\n  \
{\n    &seq2RNA_template ($mode,&my_get_opt ( $cl,\
 \"-infile=\",1,1, \"-outfile=\",1,0));\n  }\nelsi\
f ( $mode eq \"tm_template\")\n  {\n    &seq2tm_te\
mplate ($mode, \"\", &my_get_opt ( $cl, \"-infile=\
\",1,1,\"-arch=\",1,1,\"-psv=\",1,1, \"-outfile=\"\
,1,0,));\n  }\nelsif ( $mode eq \"psitm_template\"\
)\n  {\n    &seq2tm_template ($mode,&my_get_opt ( \
$cl, \"-database=\",1,0, \"-infile=\",1,1, \"-arch\
=\",1,1,\"-psv=\",1,1, \"-outfile=\",1,0,));\n  }\\
nelsif ( $mode eq \"ssp_template\")\n  {\n    &seq\
2ssp_template ($mode,&my_get_opt ( $cl, \"-infile=\
\",1,1,\"-seq=\",1,1,\"-obs=\",1,1, \"-outfile=\",\
1,0));\n  }\nelsif ( $mode eq \"psissp_template\")\
\n  {\n    &seq2ssp_template ($mode,&my_get_opt ( \
$cl, \"-infile=\",1,1,\"-seq=\",1,1,\"-obs=\",1,1,\
 \"-outfile=\",1,0));\n  }\n\nelsif ( $mode eq \"r\
na_pair\")\n{\n    &seq2rna_pair($mode,&my_get_opt\
 ( $cl, \"-pdbfile1=\",1,1, \"-pdbfile2=\",1,1, \"\
-method=\",1,2,\"-param=\",0,0, \"-outfile=\",1,0,\
 ));\n}\nelsif ( $mode eq \"calc_rna_template\")\n\
{\n    &calc_rna_template($mode,&my_get_opt ( $cl,\
 \"-infile=\",1,1,\"-pdbfile=\",1,1, \"-outfile=\"\
,1,0));\n}\nelse\n  {\n    myexit(flush_error( \"$\
mode is an unknown mode of tc_generic_method.pl\")\
);\n  }\nmyexit ($EXIT_SUCCESS);\n\n\nsub seq2ssp_\
template\n  {\n  my ($mode, $infile,$gor_seq,$gor_\
obs,$outfile)=@_;\n  my %s, %h;\n  my $result;\n  \
my (@profiles);\n  &set_temporary_dir (\"set\",$in\
file,\"seq.pep\");\n  %s=read_fasta_seq (\"seq.pep\
\");\n\n  \n  open (R, \">result.aln\");\n  \n  #p\
rint stdout \"\\n\";\n  foreach $seq (keys(%s))\n \
   {\n      \n      open (F, \">seqfile\");\n     \
 $s{$seq}{seq}=uc$s{$seq}{seq};\n      print (F \"\
>$s{$seq}{name}\\n$s{$seq}{seq}\\n\");\n      clos\
e (F);\n      $lib_name=\"$s{$seq}{name}.ssp\";\n \
     $lib_name=&clean_file_name ($lib_name);\n    \
  \n      if ($mode eq \"ssp_template\"){&seq2gor_\
prediction ($s{$seq}{name},$s{$seq}{seq}, \"seqfil\
e\", $lib_name,$gor_seq, $gor_obs);}\n      elsif \
($mode eq \"psissp_template\")\n	{\n	  &seq2msa_go\
r_prediction ($s{$seq}{name},$s{$seq}{seq},\"seqfi\
le\", $lib_name,$gor_seq, $gor_obs);\n	}\n    \n  \
    if ( !-e $lib_name)\n	{\n	  myexit(flush_error\
(\"GORIV failed to compute the secondary structure\
 of $s{$seq}{name}\"));\n	  myexit ($EXIT_FAILURE)\
;\n	}\n      else\n	{\n	  print stdout \"\\tProces\
s: >$s{$seq}{name} _E_ $lib_name \\n\";\n	  print \
R \">$s{$seq}{name} _E_ $lib_name\\n\";\n	}\n     \
 unshift (@profiles, $lib_name);\n    }\n  close (\
R);\n  &set_temporary_dir (\"unset\",$mode, $metho\
d,\"result.aln\",$outfile, @profiles);\n}\n\nsub s\
eq2tm_template\n  {\n  my ($mode, $db, $infile,$ar\
ch,$psv,$outfile)=@_;\n  my %s, %h;\n  my $result;\
\n  my (@profiles);\n  &set_temporary_dir (\"set\"\
,$infile,\"seq.pep\");\n  %s=read_fasta_seq (\"seq\
.pep\");\n\n  \n  open (R, \">result.aln\");\n  \n\
  #print stdout \"\\n\";\n  foreach $seq (keys(%s)\
)\n    {\n      open (F, \">seqfile\");\n      pri\
nt (F \">$s{$seq}{name}\\n$s{$seq}{seq}\\n\");\n  \
    close (F);\n      $lib_name=\"$s{$seq}{name}.t\
mp\";\n      $lib_name=&clean_file_name ($lib_name\
);\n\n      if ($mode eq \"tm_template\")\n	{\n	  \
&safe_system (\"t_coffee -other_pg fasta_seq2hmmto\
p_fasta.pl -in=seqfile -out=$lib_name -arch=$arch \
-psv=$psv\");\n	}\n      elsif ( $mode eq \"psitm_\
template\")\n	{\n	  &seq2msa_tm_prediction ($s{$se\
q}{name},$s{$seq}{seq}, $db, \"seqfile\", $lib_nam\
e,$arch, $psv);\n	}\n      if ( !-e $lib_name)\n	{\
\n	  myexit(flush_error(\"RNAplfold failed to comp\
ute the secondary structure of $s{$seq}{name}\"));\
\n	  myexit ($EXIT_FAILURE);\n	}\n      else\n	{\n\
	  print stdout \"\\tProcess: >$s{$seq}{name} _T_ \
$lib_name\\n\";\n	  print R \">$s{$seq}{name} _T_ \
$lib_name\\n\";\n	}\n      unshift (@profiles, $li\
b_name);\n    }\n  close (R);\n  &set_temporary_di\
r (\"unset\",$mode, $method,\"result.aln\",$outfil\
e, @profiles);\n}\n\nsub seq2RNA_template\n  {\n  \
my ($mode, $infile,$outfile)=@_;\n  my %s, %h, ;\n\
  my $result;\n  my (@profiles);\n  &set_temporary\
_dir (\"set\",$infile,\"seq.pep\");\n  %s=read_fas\
ta_seq (\"seq.pep\");\n\n  \n  open (R, \">result.\
aln\");\n  \n  #print stdout \"\\n\";\n  foreach $\
seq (keys(%s))\n    {\n      open (F, \">seqfile\"\
);\n      print (F \">$s{$seq}{name}\\n$s{$seq}{se\
q}\\n\");\n      close (F);\n      $lib_name=\"$s{\
$seq}{name}.rfold\";\n      $lib_name=&clean_file_\
name ($lib_name);\n      &safe_system (\"t_coffee \
-other_pg RNAplfold2tclib.pl -in=seqfile -out=$lib\
_name\");\n      \n      if ( !-e $lib_name)\n	{\n\
	 myexit(flush_error(\"RNAplfold failed to compute\
 the secondary structure of $s{$seq}{name}\"));\n	\
  myexit ($EXIT_FAILURE);\n	}\n      else\n	{\n	  \
print stdout \"\\tProcess: >$s{$seq}{name} _F_ $li\
b_name\\n\";\n	  print R \">$s{$seq}{name} _F_ $li\
b_name\\n\";\n	}\n      unshift (@profiles, $lib_n\
ame);\n    }\n  close (R);\n  &set_temporary_dir (\
\"unset\",$mode, $method,\"result.aln\",$outfile, \
@profiles);\n}\n\nsub psiblast2profile_template \n\
  {\n  my ($mode, $infile, $db, $method, $outfile)\
=@_;\n  my %s, %h, ;\n  my ($result,$psiblast_outp\
ut,$profile_name,@profiles);\n  my $trim=0;\n  &se\
t_temporary_dir (\"set\",$infile,\"seq.pep\");\n  \
%s=read_fasta_seq (\"seq.pep\");\n  open (R, \">re\
sult.aln\");\n  \n  #print stdout \"\\n\";\n  fore\
ach $seq (keys(%s))\n    {\n      open (F, \">seqf\
ile\");\n      print (F \">$A\\n$s{$seq}{seq}\\n\"\
);\n      close (F);\n      $psiblast_output=&run_\
blast ($s{$seq}{name},$method, $db, \"seqfile\",\"\
outfile\");\n      if ( -e $psiblast_output)\n	{\n\
	  %profile=blast_xml2profile($s{$seq}{name}, $s{$\
seq}{seq},$maxid, $minid,$mincov,$psiblast_output)\
;\n	  unlink ($psiblast_output);\n	  \n	  $profile\
_name=\"$s{$seq}{name}.prf\";\n	  $profile_name=&c\
lean_file_name ($profile_name);\n	  unshift (@prof\
iles, $profile_name);\n	  output_profile ($profile\
_name, \\%profile, $trim);\n	  print stdout \"\\tP\
rocess: >$s{$seq}{name} _R_ $profile_name [$profil\
e{n} Seq.] [$SERVER/blast/$db][$CACHE_STATUS]\\n\"\
;\n	  print R \">$s{$seq}{name} _R_ $profile_name\\
\n\";\n	}\n    }\n  close (R);\n  &set_temporary_d\
ir (\"unset\",$mode, $method,\"result.aln\",$outfi\
le, @profiles);\n}\n\nsub blast2pdb_template \n  {\
\n  my ($mode, $infile, $db, $method, $outfile,$ty\
pe)=@_;\n  my %s, %h, ;\n  my ($result,$blast_outp\
ut);\n  &set_temporary_dir (\"set\",$infile,\"seq.\
pep\");\n  %s=read_fasta_seq (\"seq.pep\");\n  ope\
n (R, \">result.aln\");\n  \n \n  #print stdout \"\
\\n\";\n  foreach $seq (keys(%s))\n    {\n      my\
 $c;\n      my $found;\n      \n      open (F, \">\
seqfile\");\n      print (F \">$A\\n$s{$seq}{seq}\\
\n\");\n      close (F);\n     \n      $blast_outp\
ut=&run_blast ($s{$seq}{name},$method, $db, \"seqf\
ile\",\"outfile\");\n     \n      %p=blast_xml2pro\
file($s{$seq}{name}, $s{$seq}{seq},$maxid, $minid,\
$mincov,$blast_output);\n      unlink ($blast_outp\
ut);\n      \n      $c=1;\n      print stdout \"\\\
tProcess: >$s{$seq}{name} [$SERVER/blast/$db][$CAC\
HE_STATUS]\\n\";\n      while (!$found && $c<$p{n}\
)\n	{\n	  $pdbid=&id2pdbid($p{$c}{identifyer});\n	\
  if ( length ($pdbid)>5){$pdbid=id2pdbid($p{$c}{d\
efinition});}\n	  \n	  if (!&pdb_is_released($pdbi\
d))\n	    {\n	      print stdout \"\\t\\t**$pdbid \
[PDB NOT RELEASED or WITHDRAWN]\\n\";\n	      $c++\
;\n	    }\n	  elsif (!&pdb_has_right_type ($pdbid,\
$type))\n	    {\n	      print stdout \"\\t\\t**$pd\
bid [PDB with Invalid Type ($type)]\\n\";\n	      \
$c++;\n	    }\n	  else\n	    {\n	      $found=1;\n\
	    }\n	}\n\n      if ($found)\n	{\n	  print R \"\
>$s{$seq}{name} _P_ $pdbid\\n\";\n	  print stdout \
\"\\t\\t >$s{$seq}{name} _P_ $pdbid\\n\";\n	}\n   \
   else\n	{\n	  print R \">$s{$seq}{name}\\n\";\n	\
  print stdout \"\\t\\t >$s{$seq}{name} No Templat\
e Selected\\n\";\n	}\n    }\n  close (R);\n  &set_\
temporary_dir (\"unset\",$mode, $method,\"result.a\
ln\",$outfile);\n}\nsub pdb_has_right_type\n  {\n \
   my $pdb=shift;\n    my $type=shift;\n    my $f=\
vtmpnam();\n    \n    $value= &safe_system (\"t_co\
ffee -other_pg extract_from_pdb -model_type $pdb >\
 $f\");\n    my $r=&file2string ($f);\n    chomp($\
r);\n\n        \n    if ( $r eq \"NMR\" && $type=~\
/n/){return 1;}\n    elsif ( $r eq \"diffraction\"\
 && $type=~/d/){return 1;}\n    elsif ( $r eq \"mo\
del\" && $type=~/m/){return 1;}\n    else {return \
0;}\n  }\nsub pdb_is_released\n  {\n    my $pdb=sh\
ift;\n    my $f=vtmpnam();\n    \n    $value= &saf\
e_system (\"t_coffee -other_pg extract_from_pdb -i\
s_released_pdb_name $pdb > $f\");\n    my $r=&file\
2string ($f);\n    chomp($r);\n    return $r;\n  }\
\nsub blast_msa\n  {\n    my ($infile,$db,$outfile\
)=@_;\n    my ($a, %seq);\n    my $seqfile;\n    m\
y $SEQ=new FileHandle;\n    my $seqfile=\"seqfile\\
";\n    my @txt;\n    \n    \n    %s1=&read_fasta_\
seq ($db);\n    \n    foreach $s (keys (%s1))\n   \
   {\n	$i=$s1{$s}{order};\n	$s{$i}{name}=$s;\n	$s{\
$i}{seq}=$s1{$s}{seq};\n	$s{$i}{len}=length( $s{$i\
}{seq});\n	$s{n}++;\n      }\n    \n    &safe_syst\
em (\"formatdb -i $db\");\n    &safe_system  (\"bl\
astall -i $infile -d $db -m7 -p blastp -o io\");\n\
    &set_blast_type (\"io\");\n    \n    %FB=&xml2\
tag_list (\"io\", \"Iteration\");\n    open (F, \"\
>$outfile\");\n    print F \"! TC_LIB_FORMAT_01\\n\
\";\n    print F \"$s{n}\\n\";\n    for ( $a=0; $a\
<$s{n}; $a++)\n      {\n	print F \"$s{$a}{name} $s\
{$a}{len} $s{$a}{seq}\\n\";\n      }\n\n\n    for \
( $a=0; $a<$FB{n}; $a++)\n      {\n	%p=blast_xml2p\
rofile ($s{$a}{name}, $s{$a}{seq},100, 0, 0, $FB{$\
a}{body});\n	my $query=$p{0}{name};\n	my $i= $s1{$\
query}{order}+1;\n	for ($b=1; $b<$p{n}; $b++)\n	  \
{\n	    my $l=length ($p{$b}{Qseq});\n	    my $hit\
=$p{$b}{definition};\n	    my $Qstart=$p{$b}{Qstar\
t};\n	    my $Hstart=$p{$b}{Hstart};\n	    my $ide\
ntity=$p{$b}{identity};\n	    my @lrQ=split (//,$p\
{$b}{Qseq});\n	    my @lrH=split (//,$p{$b}{Hseq})\
;\n	    \n	    my $j= $s1{$hit}{order}+1;\n	    #i\
f ( $j==$i){next;}\n	    printf F \"# %d %d\\n\", \
$i, $j;\n	    #  print  F \"\\n$p{$b}{Qseq} ($Qsta\
rt)\\n$p{$b}{Hseq} ($Hstart)\";\n	    for ($c=0; $\
c<$l; $c++)\n	      {\n		my $rQ=$lrQ[$c];\n		my $r\
H=$lrH[$c];\n		my $n=0;\n		\n		if ($rQ ne \"-\"){$\
n++, $Qstart++;}\n		if ($rH ne \"-\"){$n++; $Hstar\
t++;}\n		\n		if ( $n==2)\n		  {\n		    printf F \"\
\\t%d %d %d\\n\", $Qstart-1, $Hstart-1,$identity;\\
n		  }\n	      }\n	  }\n      }\n    print F \"! S\
EQ_1_TO_N\\n\";\n    close (F);\n    return $outpu\
t;\n  \n  }\n\nsub blast_msa_old\n  {\n    my ($in\
file,$outfile)=@_;\n    my ($a, %seq);\n    %s1=&r\
ead_fasta_seq ($infile);\n    foreach $s (keys (%s\
1))\n      {\n	$i=$s1{$s}{order};\n	$s{$i}{name}=$\
s;\n	$s{$i}{seq}=$s1{$s}{seq};\n	$s{$i}{len}=lengt\
h( $s{$i}{seq});\n	$s{n}++;\n      }\n    &safe_sy\
stem (\"formatdb -i $infile\");\n    &safe_system \
(\"blastall -i $infile -d $infile -m7 -o io\");\n \
   &set_blast_type (\"io\");\n    \n    %FB=&xml2t\
ag_list (\"io\", \"Iteration\");\n    \n    open (\
F, \">$outfile\");\n    print F \"! TC_LIB_FORMAT_\
01\\n\";\n    print F \"$s{n}\\n\";\n    for ( $a=\
0; $a<$s{n}; $a++)\n      {\n	print F \"$s{$a}{nam\
e} $s{$a}{len} $s{$a}{seq}\\n\";\n      }\n    for\
 ( $a=0; $a<$FB{n}; $a++)\n      {\n	%p=blast_xml2\
profile ($s{$a}{name}, $s{$a}{seq},100, 0, 0, $FB{\
$a}{body});\n	for ($b=1; $b<$p{n}; $b++)\n	  {\n	 \
   my $l=length ($p{$b}{Qseq});\n	    my $hit=$p{$\
b}{definition};\n	    my $Qstart=$p{$b}{Qstart};\n\
	    my $Hstart=$p{$b}{Hstart};\n	    my $identity\
=$p{$b}{identity};\n	    my @lrQ=split (//,$p{$b}{\
Qseq});\n	    my @lrH=split (//,$p{$b}{Hseq});\n	 \
   my $i= $s1{$s{$a}{name}}{order}+1;\n	    my $j=\
 $s1{$hit}{order}+1;\n	    #if ( $j==$i){next;}\n	\
    printf F \"# %d %d\\n\", $i, $j;\n	    #  prin\
t  F \"\\n$p{$b}{Qseq} ($Qstart)\\n$p{$b}{Hseq} ($\
Hstart)\";\n	    for ($c=0; $c<$l; $c++)\n	      {\
\n		my $rQ=$lrQ[$c];\n		my $rH=$lrH[$c];\n		my $n=\
0;\n		\n		if ($rQ ne \"-\"){$n++, $Qstart++;}\n		i\
f ($rH ne \"-\"){$n++; $Hstart++;}\n		\n		if ( $n=\
=2)\n		  {\n		    printf F \"\\t%d %d %d\\n\", $Qs\
tart-1, $Hstart-1,$identity;\n		  }\n	      }\n	  \
}\n      }\n    print F \"! SEQ_1_TO_N\\n\";\n    \
close (F);\n    return $output;\n  \n  }\n\nsub se\
q2msa\n  {\n    my ($mode, $infile, $method, $para\
m, $outfile,$database)=@_;\n    &set_temporary_dir\
 (\"set\",$infile,\"seq.pep\", $database, \"db.pep\
\");\n    $param.=\" >/dev/null 2>&1 \";\n    \n  \
  \n    #make sure test.pep is in FASTA\n    &safe\
_system (\"t_coffee -other_pg seq_reformat -in seq\
.pep -output fasta_seq > x\");\n    `mv x seq.pep`\
;\n    \n    if ( $method eq \"blastp\")\n      {\\
n	&blast_msa (\"seq.pep\", \"db.pep\",\"result.aln\
\");\n      }\n    elsif ( $method eq \"muscle\")\\
n      {\n	`muscle -in seq.pep -out result.aln $pa\
ram`;\n      }\n    elsif ( $method eq \"probcons\\
")\n      {\n	`probcons seq.pep >result.aln 2>/dev\
/null`;\n      }\n    elsif ( $method eq \"mafft\"\
)\n      {\n	`mafft --quiet --localpair --maxitera\
te 1000 seq.pep> result.aln  2>/dev/null`\n      }\
\n    elsif ( $method=~/prank/)\n      {\n	`$metho\
d -d=seq.pep -o=result.aln -quiet 2>/dev/null`;\n	\
`mv result.aln.1.fas result.aln`;\n      }\n    el\
se\n      {\n	`$method -infile=seq.pep -outfile=re\
sult.aln`;\n      }\n    \n    &set_temporary_dir \
(\"unset\",$mode, $method,\"result.aln\",$outfile)\
;\n    myexit ($EXIT_SUCCESS);\n  }\n\nsub seq2thr\
ead_pair\n  {\n    my ($mode, $infile, $pdbfile1, \
$method, $param, $outfile)=@_;\n    &set_temporary\
_dir (\"set\",$infile,\"seq.pep\",$pdbfile1,\"stru\
c.pdb\");\n    if ($method eq \"fugueali\")\n     \
 {\n	#Env Variable that need to be defined for Fug\
ue\n	if (!$ENV{FUGUE_LIB_LIST}){$ENV{FUGUE_LIB_LIS\
T}=\"DUMMY\";}\n	if (!$ENV{HOMSTRAD_PATH})  {$ENV{\
HOMSTRAD_PATH}=\"DUMMY\";}\n	if (!$ENV{HOMS_PATH})\
{$ENV{HOMS_PATH}=\"DUMMY\";}\n	\n	`joy struc.pdb >\
x 2>x`;\n	&check_file(\"struc.tem\", \"Joy failed \
[FATAL:$PROGRAM/$method]\");\n	`melody -t struc.te\
m >x 2>x`;\n	&check_file(\"struc.tem\", \"Melody f\
ailed [FATAL:$PROGRAM/$method]\");\n	`fugueali -se\
q seq.pep -prf struc.fug -print > tmp_result.aln`;\
\n	\n	&check_file(\"tmp_result.aln\", \"Fugue fail\
ed [FATAL:$PROGRAM/$method]\");\n	&safe_system (\"\
t_coffee -other_pg seq_reformat -in tmp_result.aln\
 -output fasta_aln >result.aln\");\n      }\n    e\
lsif ( $method eq \"t_coffee\")\n      {\n	&safe_s\
ystem (\"t_coffee -in Pstruc.pdb Sseq.pep Mslow_pa\
ir -outfile result.aln -quiet\");\n      }\n    el\
se\n      {\n	&safe_system (\"$method -infile=seq.\
pep -pdbfile1=struc.pdb -outfile=result.aln $param\
>x 2>x\");\n      }\n    &set_temporary_dir (\"uns\
et\",$mode,$method,\"result.aln\",$outfile);\n    \
myexit ($EXIT_SUCCESS);\n  }\nsub seq2pdbid_pair\n\
  {\n    my ($mode, $pdbfile1, $pdbfile2, $method,\
 $param, $outfile)=@_;\n    my ($name);\n\n    \n \
   &set_temporary_dir (\"set\");\n    $name=$pdbfi\
le1.\" \".$pdbfile2;\n\n    if (    &cache_file(\"\
GET\",\"\",\"$name\",\"$method\",\"dali\",$outfile\
,\"EBI\"))\n      {return $outfile;}\n    else\n  \
    {\n	if ($method eq \"daliweb\")\n	  {\n	    $p\
dbfile1=~/(....)(.)/;\n	    $id1=$1; $c1=$2;\n	   \
 \n	    $pdbfile2=~/(....)(.)/;\n	    $id2=$1; $c2\
=$2;\n	    \n	    $command=\"t_coffee -other_pg da\
lilite.pl --pdb1 $id1 --chainid1 $c1 --pdb2 $id2 -\
-chainid2 $c2 --email=$EMAIL  >dali_stderr 2>dali_\
stderr\";\n	    $dali=`$command`;\n	    \n	    ope\
n (F, \"dali_stderr\");\n	    while (<F>)\n	      \
{\n		if ( /JobId: dalilite-(\\S+)/)\n		{\n		  $job\
id=$1;\n		}\n	      }\n	    close (F);\n	    unlin\
k (\"dali_stderr\");\n	    \n	    $output1=\"dalil\
ite-$jobid.txt\";\n	    if ( -e $output1)\n	      \
{\n		unlink ($output1);\n		&url2file (\"http://www\
.ebi.ac.uk/Tools/es/cgi-bin/jobresults.cgi/dalilit\
e/dalilite-$jobid/aln.html\", \"output2\");\n		\n	\
	if ( -e \"output2\")\n		  {\n		    my ($seq1, $se\
q2);\n		    $seq1=$seq2=\"\";\n		    \n		    open \
(F, \"output2\");\n		    while (<F>)\n		      {\n	\
		$l=$_;\n			if ( $l=~/Query\\s+(\\S+)/)\n			  {\n\
			    $seq1.=$1;\n			  }\n			elsif ( $l=~/Sbjct\\\
s+(\\S+)/)\n			  {\n			    $seq2.=$1;\n			  }\n		 \
     }\n		    close (F);\n		    unlink (\"output2\\
");\n		    if ($seq1 ne \"\" && $seq2 ne \"\")\n		\
      {\n			$output3=\">$A\\n$seq1\\n>$B\\n$seq2\\\
n\";\n			$output3=~s/\\./-/g;\n			open (F, \">resu\
lt.aln\");\n			print F \"$output3\";\n			close (F)\
;\n		      }\n		  }\n	      }\n	  }\n      }\n    \
&cache_file(\"SET\",\"\",\"$name\",\"$method\",\"d\
ali\",\"result.aln\",\"EBI\");\n    &set_temporary\
_dir (\"unset\",$mode, $method, \"result.aln\",$ou\
tfile);\n    myexit ($EXIT_SUCCESS);\n  }\nsub seq\
2pdb_pair\n  {\n    my ($mode, $pdbfile1, $pdbfile\
2, $method, $param, $outfile)=@_;\n    \n    &set_\
temporary_dir (\"set\",$pdbfile1,\"pdb1.pdb\",$pdb\
file2,\"pdb2.pdb\");\n    if ($method eq \"t_coffe\
e\")\n      {\n	&safe_system (\"t_coffee -in Ppdb1\
.pdb Ppdb2.pdb -quiet -outfile=result.aln\");\n   \
   }\n    elsif ( $method eq \"DaliLite\")\n      \
{\n	if ( &safe_system (\"DaliLite -pairwise pdb1.p\
db pdb2.pdb >tmp1\")==$EXIT_SUCCESS)\n	  {\n	     \
my ($seq1, $seq2);\n	     $seq1=$seq2=\"\";\n		   \
 \n	     open (F, \"tmp1\");\n	     while (<F>)\n	\
       {\n		 $l=$_;\n		 if ( $l=~/Query\\s+(\\S+)/\
)\n		   {\n		     $seq1.=$1;\n		   }\n		 elsif ( $\
l=~/Sbjct\\s+(\\S+)/)\n		   {\n		     $seq2.=$1;\n\
		   }\n	       }\n	     close (F);\n	     unlink \
(\"tmp1\");\n	     if ($seq1 ne \"\" && $seq2 ne \\
"\")\n	       {\n		 my $output3=\">$A\\n$seq1\\n>$\
B\\n$seq2\\n\";\n		 $output3=~s/\\./-/g;\n		 open \
(F, \">result.aln\");\n		 print F \"$output3\";\n	\
	 close (F);\n	       }\n	   }\n	else\n	  {\n	    \
print \"ERROR: DalLite failed to align the conside\
red structures[tc_generic_method.pl]\\n\";\n	  }  \
  \n      }\n    elsif ( $method eq \"TMalign\")\n\
      {\n	if ( &safe_system (\"TMalign pdb1.pdb pd\
b2.pdb >tmp1\")==$EXIT_SUCCESS)\n	  {\n	    `tail \
-4 tmp1 > tmp2`;\n	    \n	    open (F, \"tmp2\");\\
n	    while (<F>)\n	      {\n		unshift(@l, $_);\n	\
      }\n	    close (F);\n	    open (F, \">result.\
aln\");\n	    $l[3]=~s/[^a-zA-Z0-9-]/\\-/g;\n	    \
$l[1]=~s/[^a-zA-Z0-9-]/\\-/g;\n	    print F \">$A\\
\n$l[3]\\n>$B\\n$l[1]\\n\";\n	    close (F);\n	  }\
\n	else\n	  {\n	    print \"ERROR: TMalign failed \
to align the considered structures[tc_generic_meth\
od.pl]\\n\";\n	    `rm result.aln >/dev/null 2>/de\
v/null`;\n	  }\n      }\n    elsif ( $method eq \"\
mustang\")\n      {\n	if ( &safe_system (\"mustang\
 -i pdb1.pdb pdb2.pdb -F fasta >/dev/null 2>/dev/n\
ull\")==$EXIT_SUCCESS)\n	  {\n	    `mv results.afa\
sta result.aln`;\n	  }\n	else\n	  {\n	    print \"\
ERROR: mustang failed to align the considered stru\
ctures[tc_generic_method.pl]\\n\";\n	    `rm resul\
t.aln >/dev/null 2>/dev/null`;\n	  }\n      }\n   \
 else\n      {\n	if ( &safe_system (\"$method -pdb\
file1=pdb1.pep -pdbfile2=pdb2.pdb -outfile=result.\
aln $param>x 2>x\")==$EXIT_SUCCESS)\n	  {\n	    `m\
v results.afasta result.aln`;\n	  }\n	else\n	  {\n\
	    print \"ERROR: $method failed to align the co\
nsidered structures[tc_generic_method.pl]\\n\";\n	\
    `rm result.aln >/dev/null 2>/dev/null`;\n	  }\\
n      }\n    &set_temporary_dir (\"unset\",$mode,\
 $method, \"result.aln\",$outfile);\n    myexit ($\
EXIT_SUCCESS);\n  }\n\nsub seq2profile_pair\n  {\n\
    my ($mode, $profile1, $profile2, $method, $par\
am, $outfile)=@_;\n    \n    \n    if ($method eq \
\"clustalw\")\n      {\n	&set_temporary_dir (\"set\
\",$profile1,\"prf1.aln\",$profile2,\"prf2.aln\");\
\n	`clustalw -profile1=prf1.aln -profile2=prf2.aln\
 -outfile=result.aln`;\n	&set_temporary_dir (\"uns\
et\",$mode, $method, \"result.aln\",$outfile);\n  \
    }\n    elsif ( $method eq \"hhalign\")\n      \
{\n	hhalign ( $profile1,$profile2,$outfile,$param)\
;\n      }\n    else\n      {\n	\n	`$method -profi\
le1=prf1.aln -profile2=prf2.aln -outfile=result.al\
n $param>x 2>x`;\n      }\n   \n    myexit ($EXIT_\
SUCCESS);\n  }\n\nsub pg_is_installed\n  {\n    my\
 @ml=@_;\n    my ($r, $p, $m);\n    my $supported=\
0;\n    \n    my $p=shift (@ml);\n    if ($p=~/::/\
)\n      {\n	if (safe_system (\"perl -M$p -e 1\")=\
=$EXIT_SUCCESS){return 1;}\n	else {return 0;}\n   \
   }\n    else\n      {\n	$r=`which $p 2>/dev/null\
`;\n	if ($r eq \"\"){$r=0;}\n	else {$r1;}\n\n	if (\
$r==0 && is_blast_package ($p)){return pg_is_insta\
lled (\"legacy_blast.pl\");}\n	else {return $r;}\n\
      }\n  }\n\nsub is_blast_package\n  {\n    my \
$p=shift;\n    if ( $p=~/blastp/){return 1;}\n    \
elsif ($p=~/blastall/){return 1;}\n    elsif ($p=~\
/blastn/){return 1;}\n    elsif ($p=~/blastx/){ret\
urn 1;}\n    elsif ($p=~/formatdb/){return 1;}\n  \
  else {return 0;}\n  }\n    \nsub check_internet_\
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
    {\n	return 1;\n      }\n  }\nsub set_temporary\
_dir\n  {\n    my @list=@_;\n    my $dir_mode, $a,\
 $mode, $method;\n  \n    $dir_mode=shift (@list);\
\n\n    \n    if ( $dir_mode eq \"set\")\n      {\\
n	$initial_dir=cwd();\n	if ( !$tmp_dir)\n	  {\n	  \
  $rand=rand (100000);\n	    $tmp_dir=\"$TMPDIR/tm\
p4tcoffee_profile_pair_dir_$$\\_P_$rand\";\n	  }\n\
	if ( !-d $tmp_dir)\n	  {\n	    push (@TMPDIR_LIST\
, $tmp_dir);\n	    `mkdir $tmp_dir`;\n	  }\n	\n	fo\
r ( $a=0; $a<=$#list; $a+=2)\n	      {\n		if (-e $\
list[$a]){ `cp $list[$a] $tmp_dir/$list[$a+1]`;}\n\
	      }\n	chdir $tmp_dir;\n      }\n    elsif ( $\
dir_mode eq \"unset\")\n      {\n	$mode=shift (@li\
st);\n	$method=shift (@list);\n	\n	if (!-e $list[0\
])\n	  {\n	   myexit(flush_error(\"Program $method\
 failed to produce $list[1]\" ));\n	    myexit ($E\
XIT_FAILURE);\n	  }\n	else\n	  {\n	    chdir $init\
ial_dir;\n	    # `t_coffee -other_pg seq_reformat \
-in $tmp_dir/$list[0] -output fasta_aln -out $tmp_\
dir/result2.aln`;\n	    `cp $tmp_dir/$list[0] $tmp\
_dir/result2.aln`;\n	    if ( $list[1] eq \"stdout\
\")\n	      {\n		open (F, \"$tmp_dir/result2.aln\"\
);\n		while (<F>){print $_;}close(F);\n	      }\n	\
    else\n	      {\n		`mv $tmp_dir/result2.aln $li\
st[1]`;\n	      }\n	    shift (@list); shift (@lis\
t);\n	    foreach $f (@list)\n	      {\n		if (-e (\
\"$tmp_dir/$f\")){`mv $tmp_dir/$f .`;}\n	      }\n\
	  }\n      }\n  }\n\n\n\n\nsub my_get_opt\n  {\n \
   my @list=@_;\n    my $cl, $a, $argv, @argl;\n  \
  \n    @argl=();\n    $cl=shift @list;\n    for (\
 $a=0; $a<=$#list; $a+=3)\n      {\n	$option=$list\
[$a];\n	$optional=$list[$a+1];\n	$status=$list[$a+\
2];\n	$argv=\"\";\n	if ($cl=~/$option(\\S+)/){$arg\
v=$1;}\n	@argl=(@argl,$argv);\n	\n	\n	#$optional:0\
=>optional\n	#$optional:1=>must be set\n	#$status:\
 0=>no requirement\n	#$status: 1=>must be an exist\
ing file\n	#$status: 2=>must be an installed packa\
ge\n	\n\n	if ($optional==0){;}\n	elsif ( $optional\
==1 && $argv eq \"\")\n	  {\n	    myexit(flush_err\
or( \"ERROR: Option $option must be set\"));\n	   \
 myexit ($EXIT_FAILURE);\n	  }\n	if ($status==0){;\
}\n	elsif ($status ==1 && $argv ne \"\" && !-e $ar\
gv)\n	  {\n	    myexit(flush_error( \"File $argv m\
ust exist\"));\n	    myexit ($EXIT_FAILURE);\n	  }\
\n	elsif ( $status==2 && $argv ne \"\" && &check_p\
g_is_installed ($argv)==0)\n	  {\n	    myexit(flus\
h_error( \" $argv is not installed\"));\n	    myex\
it ($EXIT_FAILURE);\n	  }\n      }\n\n    return @\
argl;\n    }\n\nsub check_file \n  {\n    my ($fil\
e, $msg)=@_;\n\n    if ( !-e $file)\n      {\n	mye\
xit(flush_error(\"$msg\"));\n      }\n    }\nsub h\
halign\n  {\n    my ($aln1, $aln2, $outfile, $para\
m)=@_;\n    my $h1, $h2;\n    \n    $h{0}{index}=0\
;\n    $h{1}{index}=1;\n    \n    $h{0}{aln}=$aln1\
;\n    $h{1}{aln}=$aln2;\n\n   \n\n    %{$h{0}}=al\
n2psi_profile (%{$h{0}});\n    %{$h{1}}=aln2psi_pr\
ofile (%{$h{1}});\n\n    $param=~s/#S/ /g;\n    $p\
aram=~s/#M/\\-/g;\n    $param=~s/#E/\\=/g;\n    \n\
\n    \n    $command=\"hhalign -i $h{0}{a3m} -t $h\
{1}{a3m} -tc $outfile.tmp -rank 1 -mapt 0 $param\"\
;\n    `$command`;\n    \n  #  `hhalign -i $h{0}{a\
3m} -t $h{1}{a3m} -tc $outfile.tmp -rank 1 -mapt 0\
 -gapf 0.8 -gapg 0.8`;\n    \n\n    # To run globa\
l use the following\n    \n    open (I, \"$outfile\
.tmp\");\n    open (O, \">$outfile\");\n    $h{0}{\
cons}=s/\\./x/g;\n    $h{1}{cons}=s/\\./x/g;\n\n  \
  print O \"! TC_LIB_FORMAT_01\\n2\\n$h{0}{name} $\
h{0}{len} $h{0}{seq}\\n$h{1}{name} $h{1}{len} $h{1\
}{seq}\\n#1 2\\n\";\n    \n    while (<I>)\n      \
{\n	if (/(\\d+)\\s+(\\d+)\\s+(\\d+)/)\n	  {\n	    \
print O \"\\t$h{0}{$1}\\t$h{1}{$2}\\t$3\\n\";\n	  \
}\n      }\n    print O \"! SEQ_1_TO_N\\n\";\n\n  \
  close (O);\n    close (I);\n  }\n\nsub aln2psi_p\
rofile\n  {\n    my (%h)=@_;\n    my ($aln,$i,$hv,\
 $a, @c, $n);\n   \n    $i=$h{index};\n    $aln=$h\
{aln};\n\n    `cp $aln $$.hhh_aln`;\n    $command=\
\"t_coffee -other_pg seq_reformat -in $aln -output\
 hasch\";\n    $hv=`$command`;chomp ($hv);\n    \n\
    $h{a2m}=\"$tmp/$hv.tmp4hhpred.a2m\";\n    $h{a\
3m}=\"$tmp/$hv.tmp4hhpred.a3m\";\n    if ( -e $h{a\
3m}){;}\n    else\n      {\n	`hhconsensus  -M 50 -\
i $h{aln} -oa2m $h{a2m}`;\n	if (!-e $h{a2m})\n	  {\
\n	    print STDERR \"Program tc_generic_method.pl\
 FAILED to run:\\n\\thhconsensus  -M 50 -i $h{aln}\
 -oa2m $h{a2m}\";\n	    myexit ($EXIT_FAILURE);\n	\
  }\n	\n	`hhconsensus  -M 50 -i $h{aln} -oa3m $h{a\
3m}`;\n	if (!-e $h{a3m})\n	  {\n	    print STDERR \
\"Program tc_generic_method.pl FAILED to run:\\n\\\
thhconsensus  -M 50 -i $h{aln} -oa3m $h{a3m}\";\n	\
    myexit ($EXIT_FAILURE);\n	  }\n       `buildal\
i.pl $h{a3m} -n 1`;\n      }\n    \n    \n    $h{a\
2m_seq}=`head -n 2 $h{a2m} | grep -v \">\"`;chomp \
($h{a2m_seq});\n    $h{a3m_seq}=`head -n 2 $h{a3m}\
 | grep -v \">\"`;chomp ($h{a3m_seq});\n    $h{con\
s}=$h{a2m_seq};\n    $h{seq}=`head -n 2 $h{aln} | \
grep -v \">\"`;chomp ($h{seq});\n    \n    \n\n   \
 @c=split (//, $h{cons});\n    $h{len}=$#c+1;\n   \
 for ($n=0,$a=0, $b=0; $a<$h{len};$a++)\n      {\n\
	if ( $c[$a]=~/[A-Z]/)\n	  {\n	    $h{++$n}=++$b;\\
n\n	  }\n	elsif ( $c[$a]=~/[a-z\\.]/)\n	  {\n	    \
++$b;\n	  }\n      }\n    \n    $name=`head -n 2 $\
h{aln} | grep \">\"`;\n    $name=~/\\>(\\S+)/;\n  \
  $h{name}=$1;\n    \n    `cp $h{a2m} $i.a2m`;\n  \
  `cp $h{a3m} $i.a3m`;\n    `cp $h{aln} $i.hh_aln`\
;\n    \n    return %h;\n  }\n\nsub read_fasta_seq\
 \n  {\n    my $f=@_[0];\n    my %hseq;\n    my (@\
seq, @com, @name);\n    my ($a, $s,$nseq);\n\n    \
open (F, $f);\n    while (<F>)\n      {\n	$s.=$_;\\
n      }\n    close (F);\n\n    \n    @name=($s=~/\
>(\\S*).*\\n[^>]*/g);\n    \n    @seq =($s=~/>.*.*\
\\n([^>]*)/g);\n    @com =($s=~/>\\S*(.*)\\n([^>]*\
)/g);\n\n    \n    $nseq=$#name+1;\n    \n    for \
($a=0; $a<$nseq; $a++)\n      {\n	my $s;\n	my $n=$\
name[$a];\n	$hseq{$n}{name}=$n;\n	$seq[$a]=~s/[^A-\
Za-z]//g;\n	$hseq{$n}{order}=$a;\n	$hseq{$n}{seq}=\
$seq[$a];\n	$hseq{$n}{com}=$com[$a];\n	\n      }\n\
    return %hseq;\n  }\n\nsub file_contains \n  {\\
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
n\nsub my_get_opt\n  {\n    my @list=@_;\n    my $\
cl, $a, $argv, @argl;\n    \n    @argl=();\n    $c\
l=shift @list;\n    for ( $a=0; $a<=$#list; $a+=3)\
\n      {\n	$option=$list[$a];\n	$optional=$list[$\
a+1];\n	$status=$list[$a+2];\n	$argv=\"\";\n	if ($\
cl=~/$option(\\S+)/){$argv=$1;}\n	@argl=(@argl,$ar\
gv);\n	\n	\n	#$optional:0=>optional\n	#$optional:1\
=>must be set\n	#$status: 0=>no requirement\n	#$st\
atus: 1=>must be an existing file\n	#$status: 2=>m\
ust be an installed package\n	\n\n	if ($optional==\
0){;}\n	elsif ( $optional==1 && $argv eq \"\")\n	 \
 {\n\n	    myexit(flush_error(\"Option $option mus\
t be set\"));\n	   \n	  }\n	if ($status==0){;}\n	e\
lsif ($status ==1 && $argv ne \"\" && !-e $argv)\n\
	  {\n	     myexit(flush_error(\"File $argv must e\
xist\"));\n	   \n	  }\n	elsif ( $status==2 && $arg\
v ne \"\" && &check_pg_is_installed ($argv)==0)\n	\
  {\n	    myexit(flush_error(\"$argv is not instal\
led\"));\n	   \n	  }\n      }\n\n    return @argl;\
\n    }\n\nsub tag2value \n  {\n    \n    my $tag=\
(@_[0]);\n    my $word=(@_[1]);\n    my $return;\n\
    \n    $tag=~/$word=\"([^\"]+)\"/;\n    $return\
=$1;\n    return $return;\n  }\n      \nsub hit_ta\
g2pdbid\n  {\n    my $tag=(@_[0]);\n    my $pdbid;\
\n       \n    $tag=~/id=\"(\\S+)\"/;\n    $pdbid=\
$1;\n    $pdbid=~s/_//;\n    return $pdbid;\n  }\n\
sub id2pdbid\n  {\n    my $in=@_[0];\n    my $id;\\
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
n    $L=length ($seq);\n   \n    %query=&xml2tag_l\
ist ($string, \"Iteration_query-def\");\n    $name\
=$query{0}{body};\n    \n    %hit=&xml2tag_list ($\
string, \"Hit\");\n    \n    \n    for ($nhits=0,$\
a=0; $a<$hit{n}; $a++)\n      {\n	my ($ldb,$id, $i\
dentity, $expectation, $start, $end, $coverage, $r\
);\n	my (%ID,%DE,%HSP);\n	\n	$ldb=\"\";\n\n	%ID=&x\
ml2tag_list ($hit{$a}{body}, \"Hit_id\");\n	$ident\
ifyer=$ID{0}{body};\n	\n	%DE=&xml2tag_list ($hit{$\
a}{body}, \"Hit_def\");\n	$definition=$DE{0}{body}\
;\n	\n	%HSP=&xml2tag_list ($hit{$a}{body}, \"Hsp\"\
);\n	for ($b=0; $b<$HSP{n}; $b++)\n	  {\n	    my (\
%START,%END,%E,%I,%Q,%M);\n\n	 \n	    %START=&xml2\
tag_list ($HSP{$b}{body}, \"Hsp_query-from\");\n	 \
   %HSTART=&xml2tag_list ($HSP{$b}{body}, \"Hsp_hi\
t-from\");\n	    \n	    %LEN=  &xml2tag_list ($HSP\
{$b}{body}, \"Hsp_align-len\");\n	    %END=  &xml2\
tag_list ($HSP{$b}{body}, \"Hsp_query-to\");\n	   \
 %HEND=  &xml2tag_list ($HSP{$b}{body}, \"Hsp_hit-\
to\");\n	    %E=&xml2tag_list     ($HSP{$b}{body},\
 \"Hsp_evalue\");\n	    %I=&xml2tag_list     ($HSP\
{$b}{body}, \"Hsp_identity\");\n	    %Q=&xml2tag_l\
ist     ($HSP{$b}{body}, \"Hsp_qseq\");\n	    %M=&\
xml2tag_list     ($HSP{$b}{body}, \"Hsp_hseq\");\n\
	    \n	    for ($e=0; $e<$Q{n}; $e++)\n\n	      {\
\n		$qs=$Q{$e}{body};\n		$ms=$M{$e}{body};\n		\n		\
$expectation=$E{$e}{body};\n		$identity=($LEN{$e}{\
body}==0)?0:$I{$e}{body}/$LEN{$e}{body}*100;\n		$s\
tart=$START{$e}{body};\n		$end=$END{$e}{body};\n		\
$Hstart=$HSTART{$e}{body};\n		$Hend=$HEND{$e}{body\
};\n	\n		$coverage=(($end-$start)*100)/$L;\n	\n		i\
f ($identity>$maxid || $identity<$minid || $covera\
ge<$mincov){next;}\n		@lr1=(split (//,$qs));\n		@l\
r2=(split (//,$ms));\n		$l=$#lr1+1;\n		for ($c=0;$\
c<$L;$c++){$p[$nhits][$c]=\"-\";}\n		for ($d=0,$c=\
0; $c<$l; $c++)\n		  {\n		    $r=$lr1[$c];\n		    \
if ( $r=~/[A-Za-z]/)\n		      {\n			\n			$p[$nhits\
][$d + $start-1]=$lr2[$c];\n			$d++;\n		      }\n	\
	  }\n		$Qseq[$nhits]=$qs;\n		$Hseq[$nhits]=$ms;\n\
		$QstartL[$nhits]=$start;\n		$HstartL[$nhits]=$Hs\
tart;\n		$identityL[$nhits]=$identity;\n		$endL[$n\
hits]=$end;\n		$definitionL[$nhits]=$definition;\n\
		$identifyerL[$nhits]=$identifyer;\n		$comment[$n\
hits]=\"$ldb|$identifyer [Eval=$expectation][id=$i\
dentity%][start=$Hstart end=$Hend]\";\n		$nhits++;\
\n	      }\n	  }\n      }\n    \n    \n    $profil\
e{n}=0;\n    $profile{$profile{n}}{name}=$name;\n \
   $profile{$profile{n}}{seq}=$seq;\n    $profile \
{n}++;\n    \n    for ($a=0; $a<$nhits; $a++)\n   \
   {\n	$n=$a+1;\n	\n	$profile{$n}{name}=\"$name\\_\
$a\";\n	$profile{$n}{seq}=\"\";\n	$profile{$n}{Qse\
q}=$Qseq[$a];\n	$profile{$n}{Hseq}=$Hseq[$a];\n	$p\
rofile{$n}{Qstart}=$QstartL[$a];\n	$profile{$n}{Hs\
tart}=$HstartL[$a];\n	$profile{$n}{identity}=$iden\
tityL[$a];\n	$profile{$n}{definition}=$definitionL\
[$a];\n	$profile{$n}{identifyer}=$identifyerL[$a];\
\n	$profile{$n}{comment}=$comment[$a];\n\n	for ($b\
=0; $b<$L; $b++)\n	  {\n	    if ($p[$a][$b])\n	   \
   {\n		$profile{$n}{seq}.=$p[$a][$b];\n	      }\n\
	    else\n	      {\n		$profile{$n}{seq}.=\"-\";\n\
	      }\n	  }\n      }\n    \n    $profile{n}=$nh\
its+1;\n    return %profile;\n  }\nsub ebi_blast_x\
ml2profile \n  {\n    my ($name,$seq,$maxid, $mini\
d, $mincov, $string)=(@_);\n    my ($L,$l, $a,$b,$\
c,$d,$nhits,@identifyerL,$identifyer);\n    \n\n  \
  \n    $seq=~s/[^a-zA-Z]//g;\n    $L=length ($seq\
);\n    %hit=&xml2tag_list ($string, \"hit\");\n  \
  \n    for ($nhits=0,$a=0; $a<$hit{n}; $a++)\n   \
   {\n	my ($ldb,$id, $identity, $expectation, $sta\
rt, $end, $coverage, $r);\n	my (%Q,%M,%E,%I);\n	\n\
	$ldb=&tag2value ($hit{$a}{open}, \"database\");\n\
	$identifyer=&tag2value ($hit{$a}{open}, \"id\");\\
n\n	$description=&tag2value ($hit{$a}{open}, \"des\
cription\");\n	\n	%Q=&xml2tag_list ($hit{$a}{body}\
, \"querySeq\");\n	%M=&xml2tag_list ($hit{$a}{body\
}, \"matchSeq\");\n	%E=&xml2tag_list ($hit{$a}{bod\
y}, \"expectation\");\n	%I=&xml2tag_list ($hit{$a}\
{body}, \"identity\");\n	\n\n	for ($b=0; $b<$Q{n};\
 $b++)\n	  {\n\n	    $qs=$Q{$b}{body};\n	    $ms=$\
M{$b}{body};\n	    \n	    $expectation=$E{$b}{body\
};\n	    $identity=$I{$b}{body};\n	    \n	    	   \
 \n	    $start=&tag2value ($Q{$b}{open}, \"start\"\
);\n	    $end=&tag2value ($Q{$b}{open}, \"end\");\\
n	    $startM=&tag2value ($M{$b}{open}, \"start\")\
;\n	    $endM=&tag2value ($M{$b}{open}, \"end\");\\
n	    $coverage=(($end-$start)*100)/$L;\n	    \n	 \
  # print \"$id: ID: $identity COV: $coverage [$st\
art $end]\\n\";\n	    \n	    \n	    if ($identity>\
$maxid || $identity<$minid || $coverage<$mincov){n\
ext;}\n	    # print \"KEEP\\n\";\n\n	    \n	    @l\
r1=(split (//,$qs));\n	    @lr2=(split (//,$ms));\\
n	    $l=$#lr1+1;\n	    for ($c=0;$c<$L;$c++){$p[$\
nhits][$c]=\"-\";}\n	    for ($d=0,$c=0; $c<$l; $c\
++)\n	      {\n		$r=$lr1[$c];\n		if ( $r=~/[A-Za-z\
]/)\n		  {\n		    \n		    $p[$nhits][$d + $start-1\
]=$lr2[$c];\n		    $d++;\n		  }\n	      }\n	  \n	 \
   $Qseq[$nhits]=$qs;\n	    $Hseq[$nhits]=$ms;\n	 \
   $QstartL[$nhits]=$start;\n	    $HstartL[$nhits]\
=$Hstart;\n	    $identityL[$nhits]=$identity;\n	  \
  $endL[$nhits]=$end;\n	    $definitionL[$nhits]=$\
definition;\n	    $identifyerL[$nhits]=$identifyer\
;\n	    $comment[$nhits]=\"$ldb|$identifyer [Eval=\
$expectation][id=$identity%][start=$startM end=$en\
dM]\";\n	    $nhits++;\n	  }\n      }\n    \n    $\
profile{n}=0;\n    $profile{$profile{n}}{name}=$na\
me;\n    $profile{$profile{n}}{seq}=$seq;\n    $pr\
ofile {n}++;\n    \n    for ($a=0; $a<$nhits; $a++\
)\n      {\n	$n=$a+1;\n	$profile{$n}{name}=\"$name\
\\_$a\";\n	$profile{$n}{seq}=\"\";\n	$profile{$n}{\
Qseq}=$Qseq[$a];\n	$profile{$n}{Hseq}=$Hseq[$a];\n\
	$profile{$n}{Qstart}=$QstartL[$a];\n	$profile{$n}\
{Hstart}=$HstartL[$a];\n	$profile{$n}{identity}=$i\
dentityL[$a];\n	$profile{$n}{definition}=$definiti\
onL[$a];	\n	$profile{$n}{identifyer}=$identifyerL[\
$a];\n	$profile{$n}{comment}=$comment[$a];\n\n	for\
 ($b=0; $b<$L; $b++)\n	  {\n	    if ($p[$a][$b])\n\
	      {\n		$profile{$n}{seq}.=$p[$a][$b];\n	     \
 }\n	    else\n	      {\n		$profile{$n}{seq}.=\"-\\
";\n	      }\n	  }\n      }\n    $profile{n}=$nhit\
s+1;\n    \n    return %profile;\n  }\nsub output_\
profile\n  {\n    my ($outfile,$profileR, $trim)=(\
@_);\n    my ($a);\n    my %profile=%$profileR;\n \
   my $P= new FileHandle;\n    my $tmp=vtmpnam();\\
n    \n    open ($P, \">$tmp\");\n    for ($a=0; $\
a<$profile{n}; $a++)\n      {\n	print $P \">$profi\
le{$a}{name} $profile{$a}{comment}\\n$profile{$a}{\
seq}\\n\";\n      }\n    close ($P);\n\n    if ( $\
trim)\n      {\n	&safe_system (\"t_coffee -other_p\
g seq_reformat -in $tmp -action +trim _aln_%%$trim\
\\_K1 -output fasta_aln -out $outfile\");\n      }\
\n    else\n      {\n	&safe_system (\"mv $tmp $out\
file\");\n      }\n    return;\n  }\nsub blast_xml\
2hit_list\n  {\n    my $string=(@_[0]);\n    retur\
n &xml2tag_list ($string, \"hit\");\n  }\nsub xmlt\
ag2value\n  {\n    my ($string_in, $tag)=@_;\n    \
my %TAG;\n    %TAG=xml2tag_list ($string_in, $tag)\
;\n    return $TAG{0}{body};\n  }\n      \nsub xml\
2tag_list  \n  {\n    my ($string_in,$tag)=@_;\n  \
  my $tag_in, $tag_out;\n    my %tag;\n    \n    i\
f (-e $string_in)\n      {\n	$string=&file2string \
($string_in);\n      }\n    else\n      {\n	$strin\
g=$string_in;\n      }\n    $tag_in1=\"<$tag \";\n\
    $tag_in2=\"<$tag>\";\n    $tag_out=\"/$tag>\";\
\n    $string=~s/>/>##1/g;\n    $string=~s/</##2</\
g;\n    $string=~s/##1/<#/g;\n    $string=~s/##2/#\
>/g;\n    @l=($string=~/(\\<[^>]+\\>)/g);\n    $ta\
g{n}=0;\n    $in=0;$n=-1;\n  \n \n\n    foreach $t\
 (@l)\n      {\n\n	$t=~s/<#//;\n	$t=~s/#>//;\n	\n	\
if ( $t=~/$tag_in1/ || $t=~/$tag_in2/)\n	  {\n	 \n\
	    $in=1;\n	    $tag{$tag{n}}{open}=$t;\n	    $n\
++;\n	    \n	  }\n	elsif ($t=~/$tag_out/)\n	  {\n	\
    \n\n	    $tag{$tag{n}}{close}=$t;\n	    $tag{n\
}++;\n	    $in=0;\n	  }\n	elsif ($in)\n	  {\n	   \\
n	    $tag{$tag{n}}{body}.=$t;\n	  }\n      }\n  \\
n    return %tag;\n  }\n\n\nsub seq2gor_prediction\
 \n  {\n    my ($name, $seq,$infile, $outfile, $go\
r_seq, $gor_obs)=(@_);\n    my ($l);\n    \n    `g\
orIV -prd $infile -seq $gor_seq -obs $gor_obs > go\
r_tmp`;\n    open (GR, \">$outfile\");\n    open (\
OG, \"gor_tmp\");\n\n    while (<OG>)\n      {\n	\\
n	$l=$_;\n	if ($l=~/\\>/){print GR \"$l\";}\n	elsi\
f ( $l=~/Predicted Sec. Struct./)\n	  {\n	    $l=~\
s/Predicted Sec. Struct\\.//;\n	    print GR \"$l\\
";\n	  }\n      }\n    close (GR);\n    close (OG)\
;\n    return;\n  }\nsub seq2msa_tm_prediction \n \
 {\n    my ($name, $seq, $db, $infile, $outfile, $\
arch, $psv)=(@_);\n    my (%p,%gseq,%R, $blast_out\
put, %s, $l);\n    my $R2=new FileHandle;\n    my \
$db=\"uniprot\";\n    my $method=\"psitm\";\n    m\
y $SERVER=\"EBI\";\n    \n    $blast_output=&run_b\
last ($name,\"blastp\", $db, $infile, \"outfile\")\
;\n    \n    if (&cache_file(\"GET\",$infile,$name\
,$method,$db,$outfile,$SERVER))\n      {\n	print \\
"\\tPSITM: USE Cache\\n\";\n	return $outfile;\n   \
   }\n    else\n      {\n	$CACHE_STATUS=\"COMPUTE \
CACHE\";\n	%p=blast_xml2profile($name,$seq,$maxid,\
 $minid,$mincov,$blast_output);\n	\n	\n	open (F, \\
">tm_input\");\n	for (my $a=0; $a<$p{n}; $a++)\n	 \
 {\n	    my $s;\n	    \n	    $s=$p{$a}{seq};\n	   \
 $s=uc($s);\n	    print F \">$p{$a}{name}\\n$s\\n\\
";\n	    #print stdout \">$p{$a}{name}\\n$s\\n\";\\
n	  }\n	close (F);\n	print \"\\tPSITM: kept  $p{n}\
 Homologues for Sequence $p{0}{name}\\n\";\n	&safe\
_system (\"t_coffee -other_pg fasta_seq2hmmtop_fas\
ta.pl -in=tm_input -out=$outfile -output=cons -cov\
=70 -trim=95 -arch=$arch -psv=$psv\");\n	unlink (\\
"tm_input\");\n	&cache_file(\"SET\",$infile,$name,\
$method,$db,$outfile,$SERVER);\n	return;\n      }\\
n  }\n\n\nsub seq2msa_gor_prediction \n  {\n    my\
 ($name, $seq,$infile, $outfile, $gor_seq, $gor_ob\
s)=(@_);\n    my (%p,%gseq,%R, $blast_output, %s, \
$l);\n    my $R2=new FileHandle;\n    my $db=\"uni\
prot\";\n    my $method=\"psigor\";\n    my $SERVE\
R=\"EBI\";\n    \n    $blast_output=&run_blast ($n\
ame,\"blastp\", \"uniprot\", $infile, \"outfile\")\
;\n    \n    if (&cache_file(\"GET\",$infile,$name\
,$method,$db,$outfile,$SERVER))\n      {\n	print \\
"\\tPSIGOR: USE Cache\\n\";\n	return $outfile;\n  \
    }\n    else\n      {\n	$CACHE_STATUS=\"COMPUTE\
 CACHE\";\n	%p=blast_xml2profile($name,$seq,$maxid\
, $minid,$mincov,$blast_output);\n	\n	\n	open (F, \
\">gor_input\");\n	for (my $a=0; $a<$p{n}; $a++)\n\
	  {\n	    my $s;\n	    \n	    $s=$p{$a}{seq};\n	 \
   $s=uc($s);\n	    print F \">$p{$a}{name}\\n$s\\\
n\";\n	    #print stdout \">$p{$a}{name}\\n$s\\n\"\
;\n	  }\n	close (F);\n	print \"\\tGORTM: kept  $p{\
n} Homologues for Sequence $p{0}{name}\\n\";\n	&sa\
fe_system (\"t_coffee -other_pg fasta_seq2hmmtop_f\
asta.pl -in=gor_input -out=$outfile -output=cons -\
cov=70 -trim=95 -gor_seq=$gor_seq -gor_obs=$gor_ob\
s -mode=gor\");\n	unlink (\"tm_input\");\n	&cache_\
file(\"SET\",$infile,$name,$method,$db,$outfile,$S\
ERVER);\n	return;\n      }\n  }\n\n\n\nsub run_bla\
st\n  {\n    my ($name, $method, $db, $infile, $ou\
tfile, $run)=(@_);\n    if (!$run){$run=1;}\n    \\
n    \n    if (&cache_file(\"GET\",$infile,$name,$\
method,$db,$outfile,$SERVER) && is_valid_blast_xml\
 ($outfile))\n      {return $outfile;}\n    else\n\
      {\n	$CACHE_STATUS=\"COMPUTE CACHE\";\n	if ( \
$SERVER eq \"EBI_SOAP\")\n	  {\n	    &check_config\
uration (\"EMAIL\",\"SOAP::Light\",\"INTERNET\");\\
n	    \n	    $cl_method=$method;\n	    if ($cl_met\
hod =~/wu/)\n	      {\n		$cl_method=~s/wu//;\n		if\
 ( $cl_method eq \"psiblast\")\n		  {\n		    add_w\
arning($$,$$,\"PSI BLAST cannot be used with the w\
uBLAST Client. Use server=EBI Or server=LOCAL. bla\
stp will be used instead\");\n		    $cl_method=\"b\
lastp\";\n		  }\n		\n		$command=\"t_coffee -other_\
pg wublast.pl --email $EMAIL $infile -D $db -p $cl\
_method --outfile $outfile -o xml>/dev/null 2>/dev\
/null\";\n		&safe_system ( $command);\n		if (-e \"\
$outfile.xml\") {`mv $outfile.xml $outfile`;}\n	  \
    }\n	    else\n	      {\n		if ($cl_method eq \"\
psiblast\"){$cl_method =\"blastp -j5\";}\n		\n		$c\
ommand=\"t_coffee -other_pg blastpgp.pl --email $E\
MAIL $infile -d $db --outfile $outfile -p $cl_meth\
od --mode PSI-Blast>/dev/null 2>/dev/null\";\n		&s\
afe_system ( $command);\n		\n		if (-e \"$outfile.x\
ml\") {`mv $outfile.xml $outfile`;}\n	      }\n	  \
}\n	elsif ($SERVER eq \"EBI_REST\" || $SERVER eq \\
"EBI\")\n	  {\n	   \n	    $cl_method=$method;\n	  \
  &check_configuration(\"EMAIL\",\"XML::Simple\", \
\"INTERNET\");\n	    if ($db eq \"uniprot\"){$db1=\
\"uniprotkb\";}\n	    else {$db1=$db;}\n	    \n\n	\
    if ($cl_method =~/wu/)\n	      {\n		$cl_method\
=~s/wu//;\n		if ( $cl_method eq \"psiblast\")\n		 \
 {\n		    $cl_method=\"blastp\";\n		  }\n		\n		$co\
mmand=\"t_coffee -other_pg wublast_lwp.pl --email \
$EMAIL -D $db1 -p $cl_method --outfile $outfile --\
outformat xml --stype protein $infile>/dev/null 2>\
/dev/null\";\n		&safe_system ( $command,5);\n		if \
(-e \"$outfile.xml\") {`mv $outfile.xml $outfile`;\
}\n		elsif (-e \"$outfile.xml.xml\"){`mv $outfile.\
xml.xml $outfile`;}\n	      }\n	    else\n	      {\
\n		if ( $cl_method =~/psiblast/){$cl_method =\"bl\
astp -j5\";}\n		$command=\"t_coffee -other_pg ncbi\
blast_lwp.pl --email $EMAIL -D $db1 -p $cl_method \
--outfile $outfile --outformat xml --stype protein\
 $infile>/dev/null 2>/dev/null\";\n		#$command=\"t\
_coffee -other_pg ncbiblast_lwp.pl --email $EMAIL \
-D $db1 -p $cl_method --outfile $outfile --outform\
at xml --stype protein $infile>/dev/null\";\n		&sa\
fe_system ( $command,5);\n		if (-e \"$outfile.xml\\
") {`mv $outfile.xml $outfile`;}\n		elsif (-e \"$o\
utfile.xml.xml\"){`mv $outfile.xml.xml $outfile`;}\
 \n		\n	      }\n	    \n	 }\n	elsif ($SERVER eq \"\
NCBI\")\n	  {\n	    &check_configuration (\"blastc\
l3\",\"INTERNET\");\n	    if ($db eq \"uniprot\"){\
$cl_db=\"nr\";}\n	    else {$cl_db=$db;}\n	    \n	\
    if ( $method eq \"psiblast\")\n	      {\n		add\
_warning($$,$$,\"PSI BLAST cannot be used with the\
 NCBI BLAST Client. Use server=EBI Or server=LOCAL\
. blastp will be used instead\");\n		$cl_method=\"\
blastp\";\n	      }\n	    else\n	      {\n		$cl_me\
thod=$method;\n	      }\n	    $command=\"blastcl3 \
-p $cl_method -d $cl_db -i $infile -o $outfile -m \
7\";\n	    &safe_system ($command);\n	  }\n	elsif \
($SERVER =~/CLIENT_(.*)/)\n	  {\n	    my $client=$\
1;\n	    $command=\"$client -p $method -d $db -i $\
infile -o $outfile -m 7\";\n	    &safe_system ($co\
mmand);\n	  }\n	elsif ( $SERVER eq \"LOCAL_blastal\
l\")\n	  {\n	    &check_configuration (\"blastall\\
");\n	    if ($method eq \"blastp\")\n	      {\n		\
$command=\"blastall -d $db -i $infile -o $outfile \
-m7 -p blastp\";\n	      }\n	    &safe_system ($co\
mmand);\n	  }\n	elsif ( $SERVER eq \"LOCAL\")\n	  \
{\n	    \n	    if ($ENV{\"BLAST_DB_DIR\"})\n	     \
 {\n		$x=$ENV{\"BLAST_DB_DIR\"};\n		$cl_db=\"$x$db\
\";\n	      }\n	    else\n	      {\n		$cl_db=$db;\\
n	      }\n	    \n	    if ($method eq \"blastp\")\\
n	      {\n		&check_configuration(\"blastpgp\");\n\
		$command=\"blastpgp -d $cl_db -i $infile -o $out\
file -m7 -j1\";\n	      }\n	    elsif ($method eq \
\"psiblast\")\n	      {\n		&check_configuration(\"\
blastpgp\");\n		$command=\"blastpgp -d $cl_db -i $\
infile -o $outfile -m7 -j5\";\n	      }\n	    elsi\
f ($method eq \"blastn\")\n	      {\n		&check_conf\
iguration(\"blastall\");\n		$command=\"blastall -p\
 blastn -d $cl_db -i $infile -o $outfile -m7 -W6\"\
;\n	      }	\n	    &safe_system ($command);\n	  }\\
n	else\n	  {\n	    \n	    myexit(add_error (EXIT_F\
AILURE,$$,$$,getppid(), \"BLAST_FAILURE::UnknownSe\
rver\",$CL));\n	  }\n	\n	if ( !-e $outfile || !&is\
_valid_blast_xml($outfile))\n	  {\n	    \n	    if \
( -e $outfile)\n	      {\n		add_warning ($$,$$,\"C\
orrupted Blast Output (Run $run)\");\n		unlink($ou\
tfile);\n	      }\n	    \n	    if ( $run==$BLAST_M\
AX_NRUNS)\n	      {\n	\n		myexit(add_error (EXIT_F\
AILURE,$$,$$,getppid(), \"BLAST_FAILURE::UnknownRe\
ason\", \"$command\"));\n	      }\n	    else\n	   \
   {\n	\n		add_warning ($$,$$,\"Blast for $name fa\
iled (Run: $run)\");\n		\n		return run_blast ($nam\
e, $method, $db,$infile, $outfile, $run+1);\n	    \
  }\n	  }\n	\n	&cache_file(\"SET\",$infile,$name,$\
method,$db,$outfile,$SERVER);\n	return $outfile;\n\
      }\n  }\n\nsub cache_file\n  {\n    my ($cach\
e_mode,$infile,$name,$method,$db, $outfile,$server\
)=(@_);\n    my $cache_file;\n    #Protect names s\
o that they can be turned into legal filenames\n  \
  $name=&clean_file_name ($name);\n\n    if ($db=~\
/\\//)\n      {\n	$db=~/([^\\/]+)$/;\n	$db=$1;\n  \
    }\n    $cache_file_sh=\"$name.$method.$db.$ser\
ver.tmp\";\n    $cache_file=\"$CACHE/$name.$method\
.$db.$server.tmp\";\n    \n    if ($infile ne \"\"\
)\n      {\n	$cache_file_infile_sh=\"$name.$method\
.$db.$server.infile.tmp\";\n	$cache_file_infile=\"\
$CACHE/$name.$method.$db.$server.infile.tmp\";\n  \
    }\n    \n    if ($cache_mode eq \"GET\")\n    \
  {\n	if ($CACHE eq \"\" || $CACHE eq \"no\" || $C\
ACHE eq \"ignore\"  || $CACHE eq \"local\" || $CAC\
HE eq \"update\"){return 0;}\n	elsif ( !-d $CACHE)\
\n	  {\n	    print STDERR \"ERROR: Cache Dir: $CAC\
HE Does not Exist\";\n	    return 0;\n	  }\n	else\\
n	  {\n	    if ( -e $cache_file && &fasta_file1_eq\
_fasta_file2($infile,$cache_file_infile)==1)\n	   \
   {\n		`cp $cache_file $outfile`;\n		$CACHE_STATU\
S=\"READ CACHE\";\n		return 1;\n	      }\n	  }\n  \
    }\n    elsif ($cache_mode eq \"SET\")\n      {\
\n	if ($CACHE eq \"\" || $CACHE eq \"no\" || $CACH\
E eq \"ignore\"  || $CACHE eq \"local\" || $CACHE \
eq \"update\"){return 0;}\n	elsif ( !-d $CACHE)\n	\
  {\n	    print STDERR \"ERROR: Cache Dir: $CACHE \
Does not Exist\";\n	    return 0;\n	  }\n	elsif (-\
e $outfile)\n	  {\n	    `cp $outfile $cache_file`;\
\n	    if ($cache_file_infile ne \"\"){ `cp $infil\
e $cache_file_infile`;}\n\n	    #functions for upd\
ating the cache\n	    #`t_coffee -other_pg clean_c\
ache.pl -file $cache_file_sh -dir $CACHE`;\n	    #\
`t_coffee -other_pg clean_cache.pl -file $cache_fi\
le_infile_sh -dir $CACHE`;\n	    return 1;\n	  }\n\
      }\n    $CACHE_STATUS=\"COMPUTE CACHE\";\n   \
 return 0;\n  }\nsub file1_eq_file2\n  {\n    my (\
$f1, $f2)=@_;\n    if ( $f1 eq \"\"){return 1;}\n \
   elsif ( $f2 eq \"\"){return 1;}\n    elsif ( !-\
e $f1){return 0;}\n    elsif ( !-e $f2){return 0;}\
\n    elsif ($f1 eq \"\" || $f2 eq \"\" || `diff $\
f1 $f2` eq \"\"){return 1;}\n    \n    return 0;\n\
  }\nsub clean_file_name \n  {\n    my $name=@_[0]\
;\n    \n    $name=~s/[^A-Za-z1-9.-]/_/g;\n    ret\
urn $name;\n  }\nsub url2file\n  {\n    my ($addre\
ss, $out)=(@_);\n    \n    if (&pg_is_installed (\\
"wget\"))\n	{\n	  return &safe_system (\"wget $add\
ress -O$out >/dev/null 2>/dev/null\");\n	}\n    el\
sif (&pg_is_installed (\"curl\"))\n      {\n	retur\
n &safe_system (\"curl $address -o$out >/dev/null \
2>/dev/null\");\n      }\n    else\n      {\n	myex\
it(flus_error(\"neither curl nor wget are installe\
d. Imnpossible to fectch remote file\"));\n	exit (\
$EXIT_FAILURE);\n      }\n  }\nsub fasta_file1_eq_\
fasta_file2\n  {\n    my ($f1, $f2)=@_;\n    my (%\
s1, %s2);\n    my @names;\n    %s1=read_fasta_seq \
($f1);\n    %s2=read_fasta_seq ($f2);\n\n    @name\
s=(keys (%s1));\n    \n    foreach $n (keys(%s1))\\
n      {\n	if ($s1{$n}{seq} ne $s2{$n}{seq}){retur\
n 0;}\n      } \n    \n    foreach $n (keys(%s2))\\
n      {\n	if ($s1{$n}{seq} ne $s2{$n}{seq}){retur\
n 0;}\n      }\n    return 1;\n  }\n	\n\n\nsub rea\
d_template_file\n{\n	my $pdb_templates = @_[0];\n	\
open (TEMP, \"<$pdb_templates\");\n	my %temp_h;\n	\
while (<TEMP>)\n{\n		$line = $_;\n 		$line =~/(\\S\
+)\\s(\\S+)/;\n 		$temp_h{$1}= $2;\n}\n	close(TEMP\
);\n	return %temp_h;\n}\n\nsub calc_rna_template\n\
{\n	my ($mode, $infile, $pdbfile, $outfile)=@_;\n	\
my %s, %h ;\n	my $result;\n	my (@profiles);\n	&set\
_temporary_dir (\"set\",$infile,\"seq.pep\");\n	%s\
=read_fasta_seq (\"seq.pep\");\n	\n	%pdb_template_\
h = &read_template_file($pdbfile);\n	my $pdb_chain\
;\n	open (R, \">result.aln\");\n\n\n	#print stdout\
 \"\\n\";\n	foreach $seq (keys(%s))\n	{\n		if ($pd\
b_template_h{$seq} eq \"\")\n		{\n			next;\n		}\n	\
	open (F, \">seqfile\");\n		print (F \">$s{$seq}{n\
ame}\\n$s{$seq}{seq}\\n\");\n		close (F);\n		$pdb_\
chain = $pdb_template_h{$seq};\n		$lib_name=\"$s{$\
seq}{name}.rfold\";\n		$lib_name=&clean_file_name \
($lib_name);\n		\n 		safe_system (\"secondary_stru\
c.py seqfile $CACHE$pdb_chain  $lib_name\");\n		\n\
		if ( !-e $lib_name)\n		{\n		myexit(flush_error(\\
"RNAplfold failed to compute the secondary structu\
re of $s{$seq}{name}\"));\n			myexit ($EXIT_FAILUR\
E);\n		}\n		else\n		{\n			print stdout \"\\tProces\
s: >$s{$seq}{name} _F_ $lib_name\\n\";\n			print R\
 \">$s{$seq}{name} _F_ $lib_name\\n\";\n		}\n		uns\
hift (@profiles, $lib_name);\n	}\n	close (R);\n	&s\
et_temporary_dir (\"unset\",$mode, $method,\"resul\
t.aln\",$outfile, @profiles);\n}\n\n\n\nsub seq2rn\
a_pair{\n	my ($mode, $pdbfile1, $pdbfile2, $method\
, $param, $outfile)=@_;\n	\n	if ($method eq \"runs\
ara.py\")\n	{\n		open(TMP,\"<$pdbfile1\");\n		my $\
count = 0;\n		my $line;\n		while (<TMP>)\n		{\n			\
$line = $_;\n			if ($count ==1)\n			{\n				last;\n\
			}\n			$count += 1;\n		}\n\n		\n		$chain1 = subs\
tr($line,length($line)-3,1);\n\n		close TMP;\n		op\
en(TMP,\"<$pdbfile2\");\n		my $count = 0;\n		while\
 (<TMP>)\n		{\n			$line = $_;\n			if ($count ==1)\\
n			{\n				last;\n			}\n			$count += 1;\n		}\n		$c\
hain2 = substr($line,length($line)-3,1);\n		close \
TMP;\n\n		$tmp_file=&vtmpnam();\n	\n		safe_system(\
\"runsara.py $pdbfile1 $chain1 $pdbfile2 $chain2 -\
s -o $tmp_file --limitation 5000 > /dev/null 2> /d\
ev/null\") == 0 or die \"sara did not work $!\\n\"\
;\n		open(TMP,\"<$tmp_file\") or die \"cannot open\
 the sara tmp file:$!\\n\";\n		open(OUT,\">$outfil\
e\") or die \"cannot open the $outfile file:$!\\n\\
";\n\n		my $switch = 0;\n		my $seqNum = 0;\n		fore\
ach my $line (<TMP>)\n		{\n			next unless ($line=~\
/SARAALI/);\n			if ($line=~/>/)\n			{\n				$switch\
 =0;\n				print OUT \">seq$seqNum\\n\";\n				$seqN\
um++;				\n			}\n			if ($switch < 2){\n				$switch\
++;\n				next;\n			}\n	\n			if ($line =~/REMARK\\s\
+SARAALI\\s+([^\\*]+)\\*/)\n			{\n				my $string =\
 $1;\n				print OUT \"$string\\n\";\n			}\n		}\n		\
close TMP; \n		close OUT;\n		unlink($tmp_file);\n	\
}\n}\n\nsub seq2tblastx_lib\n  {\n    my ($mode, $\
infile, $outfile)=@_;\n    my (%s, $method,$nseq);\
\n\n    $method=$mode;\n    &set_temporary_dir (\"\
set\",$infile,\"infile\");\n    %s=read_fasta_seq(\
\"infile\");\n    \n    \n    foreach $seq (keys(%\
s))\n      {\n	$slist[$s{$seq}{order}]=$s{$seq}{se\
q};\n	$sname[$s{$seq}{order}]=$s{$seq}{name};\n	$s\
len[$s{$seq}{order}]=length ($s{$seq}{seq});\n    \
  }\n    $nseq=$#sname+1;\n    open (F, \">outfile\
\");\n    print F \"! TC_LIB_FORMAT_01\\n\";\n    \
print F \"$nseq\\n\";\n    for ($a=0; $a<$nseq;$a+\
+)\n      {\n	print F \"$sname[$a] $slen[$a]  $sli\
st[$a]\\n\"\n      }\n    close (F);\n    &safe_sy\
stem (\"formatdb -i infile -p F\");\n    &safe_sys\
tem (\"blastall -p tblastx -i infile -d infile -m \
7 -S1>blast.output\");\n    \n    ncbi_tblastx_xml\
2lib_file (\"outfile\", file2string (\"blast.outpu\
t\"));\n    &set_temporary_dir (\"unset\",$mode, $\
method, \"outfile\",$outfile);\n    myexit ($EXIT_\
SUCCESS);\n    }\nsub seq2tblastpx_lib\n  {\n    m\
y ($mode, $infile, $outfile)=@_;\n    my (%s, $met\
hod,$nseq);\n    $method=$mode;\n    &set_temporar\
y_dir (\"set\",$infile,\"infile\");\n    %s=read_f\
asta_seq(\"infile\");\n    \n    foreach $seq (key\
s(%s))\n      {\n	$slist[$s{$seq}{order}]=$s{$seq}\
{seq};\n	$sname[$s{$seq}{order}]=$s{$seq}{name};\n\
	$slen[$s{$seq}{order}]=length ($s{$seq}{seq});\n \
     }\n    $nseq=$#sname+1;\n    open (F, \">outf\
ile\");\n    print F \"! TC_LIB_FORMAT_01\\n\";\n \
   print F \"$nseq\\n\";\n    for ($a=0; $a<$nseq;\
$a++)\n      {\n	print F \"$sname[$a] $slen[$a]  $\
slist[$a]\\n\"\n      }\n    close (F);\n    &safe\
_system(\"t_coffee -other_pg seq_reformat -in infi\
le -output tblastx_db1 > tblastxdb\");\n    &safe_\
system (\"formatdb -i tblastxdb -p T\");\n    &saf\
e_system (\"blastall -p blastp -i tblastxdb -d tbl\
astxdb -m7 >blast.output\");\n    ncbi_tblastpx_xm\
l2lib_file (\"outfile\", file2string (\"blast.outp\
ut\"), %s);\n    &set_temporary_dir (\"unset\",$mo\
de, $method, \"outfile\",$outfile);\n    myexit ($\
EXIT_SUCCESS);\n    }\n\n\n    \n\n\n\nsub file2he\
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
    &cleanup();\n  }\n   \nsub blast_com2new_blast\
_com\n    {\n      my $com=shift;\n      if ($ENV{\
\"NCBI_BLAST_4_TCOFFEE\"} eq \"OLD\"){return $com;\
}\n      elsif (!&pg_is_installed(\"legacy_blast.p\
l\")){return $com;}\n      else \n	{\n	  if ($com=\
~/formatdb/)\n	    {\n	      $com=~s/formatdb/make\
blastdb/;\n	      $com=~s/\\-i/\\-in/;\n	      if \
($com =~/pF/){$com=~s/\\-pF/\\-dbtype nucl/;}\n	  \
    if ($com =~/p F/){$com=~s/\\-p F/\\-dbtype nuc\
l/;}\n	      $com=\"$com -logfile /dev/null\";\n	 \
     return $com;\n	    }\n	  elsif (&is_blast_pac\
kage($com))\n	    {\n	      my $path;\n	      \n	 \
     if ( $ENV{\"NCBI_BIN_4_TCOFFEE\"}){$path=$ENV\
{\"NCBI_BLAST_4_TCOFFEE\"};}\n	      else\n		{\n		\
  $path=`which legacy_blast.pl`;\n		  $path=~s/\\/\
legacy_blast\\.pl//;\n		  chomp ($path);\n		}\n	  \
    $path=\"--path $path\";\n	      if ( $com=~/\\\
>\\>/){$com=~s/\\>\\>/ $path \\>\\>/;}\n	      els\
if ( $com=~/\\>/){$com=~s/\\>/ $path \\>/;}\n	    \
  else {$com.=\" $path\";}\n	      $com=\"legacy_b\
last.pl $com\";\n	      \n	      return $com;\n	  \
  }\n	}\n    }\nsub safe_system \n{\n  my $com=shi\
ft;\n  my $ntry=shift;\n  my $ctry=shift;\n  my $p\
id;\n  my $status;\n  my $ppid=getppid();\n  if ($\
com eq \"\"){return 1;}\n  \n  if ( ($com=~/^blast\
/) ||($com=~/^formatdb/)){$com=&blast_com2new_blas\
t_com($com);} \n\n  if (($pid = fork ()) < 0){retu\
rn (-1);}\n  if ($pid == 0)\n    {\n      set_lock\
($$, \" -SHELL- $com (tcg)\");\n      exec ($com);\
\n    }\n  else\n    {\n      lock4tc ($$, \"LLOCK\
\", \"LSET\", \"$pid\\n\");#update parent\n      $\
PIDCHILD=$pid;\n    }\n  if ($debug_lock){printf S\
TDERR \"\\n\\t .... safe_system (fasta_seq2hmm)  p\
: $$ c: $pid COM: $com\\n\";}\n\n  waitpid ($pid,W\
TERMSIG);\n\n  shift_lock ($pid,$$, \"LWARNING\",\\
"LWARNING\", \"LSET\");\n\n  if ($? == $EXIT_FAILU\
RE || lock4tc($pid, \"LERROR\", \"LCHECK\", \"\"))\
\n    {\n      if ($ntry && $ctry <$ntry)\n	{\n	  \
add_warning ($$,$$,\"$com failed [retry: $ctry]\")\
;\n	  lock4tc ($pid, \"LRELEASE\", \"LERROR\", \"\\
");\n	  return safe_system ($com, $ntry, ++$ctry);\
\n	}\n      elsif ($ntry == -1)\n	{\n	  if (!shift\
_lock ($pid, $$, \"LERROR\", \"LWARNING\", \"LSET\\
"))\n	    {\n	      add_warning ($$,$$,\"$com fail\
ed\");\n	    }\n	  else\n	    {\n	      lock4tc ($\
pid, \"LRELEASE\", \"LERROR\", \"\");\n	    }\n	  \
return $?;}\n      else\n	{\n	  if (!shift_lock ($\
pid,$$, \"LERROR\",\"LERROR\", \"LSET\"))\n	    {\\
n	      myexit(add_error ($EXIT_FAILURE,$$,$pid,ge\
tppid(), \"UNSPECIFIED system\", $com));\n	    }\n\
	}\n    }\n  return $?;\n}\n\nsub check_configurat\
ion \n    {\n      my @l=@_;\n      my $v;\n      \
foreach my $p (@l)\n	{\n	  \n	  if   ( $p eq \"EMA\
IL\")\n	    { \n	      if ( !($EMAIL=~/@/))\n		{\n\
		add_warning($$,$$,\"Could Not Use EMAIL\");\n		m\
yexit(add_error ($EXIT_FAILURE,$$,$$,getppid(),\"E\
MAIL\",\"$CL\"));\n	      }\n	    }\n	  elsif( $p \
eq \"INTERNET\")\n	    {\n	      if ( !&check_inte\
rnet_connection())\n		{\n		  myexit(add_error ($EX\
IT_FAILURE,$$,$$,getppid(),\"INTERNET\",\"$CL\"));\
\n		}\n	    }\n	  elsif( $p eq \"wget\")\n	    {\n\
	      if (!&pg_is_installed (\"wget\") && !&pg_is\
_installed (\"curl\"))\n		{\n		  myexit(add_error \
($EXIT_FAILURE,$$,$$,getppid(),\"PG_NOT_INSTALLED:\
wget\",\"$CL\"));\n		}\n	    }\n	  elsif( !(&pg_is\
_installed ($p)))\n	    {\n	      myexit(add_error\
 ($EXIT_FAILURE,$$,$$,getppid(),\"PG_NOT_INSTALLED\
:$p\",\"$CL\"));\n	    }\n	}\n      return 1;\n   \
 }\n\n$program=\"T-COFFEE (Version_8.97)\";\n\n","\
*TC_METHOD_FORMAT_01\n******************generic_me\
thod.tc_method*************\n*\n*       Incorporat\
ing new methods in T-Coffee\n*       Cedric Notred\
ame 26/08/08\n*\n*********************************\
**********************\n*This file is a method fil\
e\n*Copy it and adapt it to your need so that the \
method \n*you want to use can be incorporated with\
in T-Coffee\n*************************************\
******************\n*                  USAGE      \
                        *\n***********************\
********************************\n*This file is pa\
ssed to t_coffee via -in:\n*\n*	t_coffee -in Mgene\
ric_method.method\n*\n*	The method is passed to th\
e shell using the following\n*call:\n*<EXECUTABLE>\
<PARAM1><IN_FLAG><seq_file><PARAM2><OUT_FLAG><outn\
ame><PARAM>\n*\n*Conventions:\n*<FLAG_NAME> 	<TYPE\
>		<VALUE>\n*<VALUE>:	no_name 	<=> Replaced with a\
 space\n*<VALUE>:	&nbsp	<=> Replaced with a space\\
n*\n**********************************************\
*********\n*                  ALN_MODE            \
               *\n********************************\
***********************\n*pairwise   ->all Vs all \
(no self )[(n2-n)/2aln]\n*m_pairwise ->all Vs all \
(no self)[n^2-n]^2\n*s_pairwise ->all Vs all (self\
): [n^2-n]/2 + n\n*multiple   ->All the sequences \
in one go\n*\nALN_MODE		pairwise\n*\n*************\
******************************************\n*     \
             OUT_MODE                           *\\
n*************************************************\
******\n* mode for the output:\n*External methods:\
 \n* aln -> alignmnent File (Fasta or ClustalW For\
mat)\n* lib-> Lib file (TC_LIB_FORMAT_01)\n*Intern\
al Methods:\n* fL -> Internal Function returning a\
 List (Librairie)\n* fA -> Internal Function retur\
ning an Alignmnent\n*\nOUT_MODE		aln\n************\
*******************************************\n*    \
              SEQ_TYPE                           *\
\n************************************************\
*******\n*G: Genomic, S: Sequence, P: PDB, R: Prof\
ile\n*Examples:\n*SEQTYPE	S	sequences against sequ\
ences (default)\n*SEQTYPE	S_P	sequence against str\
ucture\n*SEQTYPE	P_P	structure against structure\n\
*SEQTYPE	PS	mix of sequences and structure	\n*\nSE\
Q_TYPE	S\n*\n\n***********************************\
********************\n*                COMMAND LIN\
E                         *\n*EXECUTABLE PARAM1 IN\
_FLAG OUT_FLAG PARAM             *\n**************\
*****************************************\n*******\
************************************************\n\
*                  EXECUTABLE                     \
    *\n*******************************************\
************\n*name of the executable\n*passed to \
the shell: executable\n*	\nEXECUTABLE	tc_generic_m\
ethod.pl\n*\n*************************************\
******************\n*                  IN_FLAG    \
                         *\n**********************\
*********************************\n*IN_FLAG\n*flag\
 indicating the name of the in coming sequences\n*\
IN_FLAG S no_name ->no flag\n*IN_FLAG S &bnsp-in&b\
nsp -> \" -in \"\n*\nIN_FLAG		-infile=\n*\n*******\
************************************************\n\
*                  OUT_FLAG                       \
    *\n*******************************************\
************\n*OUT_FLAG\n*flag indicating the name\
 of the out-coming data\n*same conventions as IN_F\
LAG\n*OUT_FLAG	S no_name ->no flag\n*if you want t\
o redirect, pass the parameters via PARAM1\n*set O\
UT_FLAG to >\n*\nOUT_FLAG		-outfile=\n*\n*********\
**********************************************\n* \
                 PARAM_1                          \
    *\n*******************************************\
************\n*<EXECUTABLE><PARAM1><IN_FLAG><seq_f\
ile><PARAM2><OUT_FLAG><outname><PARAM>\n*Parameter\
s sent to the EXECUTABLE and specified *before* IN\
_FLAG \n*If there is more than 1 PARAM line, the l\
ines are\n*concatenated\n*Command_line: @EP@PARAM@\
-gapopen%e10%s-gapext%e20\n*	%s white space\n*	%e \
equal sign\n*\n*PARAM1	\n*\n*\n*\n****************\
***************************************\n*        \
          PARAM_2                              *\n\
**************************************************\
*****\n*<EXECUTABLE><PARAM1><IN_FLAG><seq_file><PA\
RAM2><OUT_FLAG><outname><PARAM>\n*Parameters sent \
to the EXECUTABLE and specified \n*after* IN_FLAG \
and *before* OUT_FLAG\n*If there is more than 1 PA\
RAM line, the lines are\n*concatenated\n*\n*PARAM1\
	\n*\n*\n*****************************************\
**************\n*                  PARAM          \
                    *\n***************************\
****************************\n*<EXECUTABLE><PARAM1\
><IN_FLAG><seq_file><PARAM2><OUT_FLAG><outname><PA\
RAM>\n*Parameters sent to the EXECUTABLE and speci\
fied *after* OUT_FLAG\n*If there is more than 1 PA\
RAM line, the lines are\n*concatenated\n*\nPARAM	-\
mode=seq_msa -method=clustalw\nPARAM   -OUTORDER=I\
NPUT -NEWTREE=core -align -gapopen=-15\n*\n*******\
************************************************\n\
*                  END                            \
    *\n*******************************************\
************\n","*TC_METHOD_FORMAT_01\n***********\
****clustalw_method.tc_method*********\nEXECUTABLE\
	clustalw\nALN_MODE		pairwise\nIN_FLAG		-INFILE=\n\
OUT_FLAG		-OUTFILE=\nOUT_MODE		aln\nPARAM		-gapope\
n=-10\nSEQ_TYPE		S\n******************************\
*******************\n","$VersionTag =             \
                                                  \
                                                  \
                  2.43;\nuse Env;\nuse FileHandle;\
\nuse Cwd;\nuse File::Path;\nuse Sys::Hostname;\no\
ur $PIDCHILD;\nour $ERROR_DONE;\nour @TMPFILE_LIST\
;\nour $EXIT_FAILURE=1;\nour $EXIT_SUCCESS=0;\n\no\
ur $REFDIR=getcwd;\nour $EXIT_SUCCESS=0;\nour $EXI\
T_FAILURE=1;\n\nour $PROGRAM=\"extract_from_pdb\";\
\nour $CL=$PROGRAM;\n\nour $CLEAN_EXIT_STARTED;\no\
ur $debug_lock=$ENV{\"DEBUG_LOCK\"};\nour $LOCKDIR\
=$ENV{\"LOCKDIR_4_TCOFFEE\"};\nif (!$LOCKDIR){$LOC\
KDIR=getcwd();}\nour $ERRORDIR=$ENV{\"ERRORDIR_4_T\
COFFEE\"};\nour $ERRORFILE=$ENV{\"ERRORFILE_4_TCOF\
FEE\"};\n&set_lock ($$);\nif (isshellpid(getppid()\
)){lock4tc(getppid(), \"LLOCK\", \"LSET\", \"$$\\n\
\");}\n      \nour $SILENT=\" >/dev/null 2>/dev/nu\
ll\";\nour $INTERNET=-1;\n\n\n\n\n\n\n\nour $BLAST\
_MAX_NRUNS=2;\nour $EXIT_SUCCESS=0;\nour $EXIT_FAI\
LURE=1;\nour $CONFIGURATION=-1;\nour $REF_EMAIL=\"\
\";\nour $PROGRAM=\"extract_from_pdb\";\n\n\nmy %o\
nelett_prot=&fill_onelett_prot();\nmy %threelett_p\
rot=&fill_threelett_prot();\nmy %onelett_RNA=&fill\
_onelett_RNA();\nmy %threelett_RNA=&fill_threelett\
_RNA();\nmy %onelett_DNA=&fill_onelett_DNA();\nmy \
%threelett_DNA=&fill_threelett_DNA();\n\n\n\n\n\nm\
y %onelett = (\n'P' => \\%onelett_prot,\n'D' => \\\
%onelett_DNA,\n'R' => \\%onelett_RNA\n);\n\n\nmy %\
threelett = (\n'P' => \\%threelett_prot,\n'D' => \\
\%threelett_DNA,\n'R' => \\%threelett_RNA\n);\n\n\\
n\n\n\n\n\nif($ARGV[0]=~/help/ ||$ARGV[0]=~/man/ |\
| $ARGV[0]=~/HELP/ || $ARGV[0]=~/Man/ || $ARGV[0] \
eq \"-h\"  || $ARGV[0] eq \"-H\"  )\n{die \"SYNTAX\
: extract_from_pdb Version $VersionTag	\n	Minimum:\
             [extract_from_pdb file] \n			   OR \n\
			     [... | extract_from_pdb]\n 	Flags (Default\
 setting on the first line)\n	   -version.........\
..........[Returns the Version Number]\n          \
 -force.....................[Forces the file to be\
 treated like a PDB file]\n                       \
               [Regenerates the header and SEQRES \
fields]\n           -force_name................[Fo\
rces the file to be named after name]]\n          \
 -infile.....file...........[Flag can be omited]\n\
			              [File must be pdb or fro pgm]\n  \
                                    [File can also\
 be compressed Z or gz]\n                         \
             [In the case of a compressed file, yo\
u can omit the gz|Z extension]\n           -netfil\
e...................[File will be fetch from the n\
et using wget]\n                                  \
    [wget or curl must be installed]\n            \
                          [ftp://ftp.gnu.org/pub/g\
nu/wget/]\n                                      [\
http://curl.haxx.se/]\n                           \
           [Must also be used to retrieve the file\
 from a local pdb copy (cf netaddress)]\n         \
  -netaddress................[Address used for the\
 retrieving the netfile]\n                        \
              [http://www.rcsb.org/pdb/cgi/export.\
cgi/%%.pdb.gz?format=PDB&pdbId=%%&compression=gz]\\
n                                      [http://www\
.expasy.ch/cgi-bin/get-pdb-entry.pl?%%]\n         \
                             [local -> will get th\
e file from pdb_dir (see pdb_dir)]\n           -ne\
tcompression............[Extension if the netfile \
comes compressed]\n                               \
       [gz]\n           -pdb_dir..................\
.[address of the repertory where the pdb is instal\
led]\n                                      [Suppo\
rts standard ftp style installation OR every stru \
in DIR]\n                                      [Gi\
ve the ..../pdb/structure/ dir]\n                 \
                     [If value omitted, the pg get\
s it from the env variable PDB_DIR]\n           -n\
etcompression_pg.........[gunzip]\n           -is_\
pdb_name..........name.[Returns 1 if the name is a\
 PDB ID, 0 otherwise]\n           -model_type.....\
......name.[Returns the model type if valid PDB na\
me]\n           -is_released_pdb_name name.[Return\
s 1 if the name corresponds to a released PDB file\
]\n           -get_pdb_chains.....name...[Returns \
the list of chains corresponding to the entry]\n  \
         -get_pdb_id.........name...[Returns the P\
DB id within the provided pdb file]\n           -g\
et_fugue_name.....name...[Turns a name into a name\
 valid for fugue]\n                               \
       [Uses the netaddress to do so]\n	   -chain.\
.....FIRST..........[Extract the first chain only]\
\n		       A B C..........[Extract Several chains \
if needed]\n		       ALL............[Extract all t\
he chains]	\n           -ligand.....ALL...........\
.[Extract the ligands in the chain (HETATM)]\n    \
                   <name1>,<name2>[Extract All the\
 named lignds]\n	   -ligand_only...............[Ex\
tract only the ligands]\n           -ligand_list..\
.............[Extract the list of ligands]\n	   -c\
oor.......<start>..<end>.[Coordinates of the fragm\
ent to extract]\n			              [Omit end to inc\
lude the Cter]\n           -num........absolute...\
....[absolute: relative to the seq] \n            \
           file...........[file: relative to file]\
\n           -num_out....new............[new: star\
t 1->L]\n                       old............[ol\
d: keep the file coordinates]\n           -delete.\
....<start>..<end>.[Delete from residue start to r\
esidue end]\n	   -atom.......CA.............[Atoms\
 to include, ALL for all of them]\n		       CA O N\
.........[Indicate several atoms if needed]\n	   -\
code.......3..............[Use the 1 letter code o\
r the 3 letters code]\n	   -mode.......raw........\
....[Output original pdb file]\n                  \
     pdb............[Output something that looks l\
ike pdb]\n		       fasta..........[Output the sequ\
ences in fasta format]\n		       simple.........[O\
utput a format easy to parse in C ]\n            -\
seq_field..ATOM...........[Field used to extract t\
he sequence]\n		       SEQRES.........[Use the com\
plete sequence]\n	   -seq.......................[E\
quivalent to  -mode fasta]\n	   -model......1.....\
.........[Chosen Model in an NMR file]\n          \
 -nodiagnostic..............[Switches Error Messag\
es off]\n           -debug.....................[Se\
ts the DEBUG ON]\n           -no_remote_pdb_dir...\
......[Do not look for a remote file]\n           \
-cache_pdb.................[Cache Value, default i\
s $HOME/.t_coffee/cache, other values: NO<=> No ca\
che]\n\n      Environement Variables\n           T\
hese variables can be set from the environement\n \
          Command line values with the correspondi\
ng flag superseed evironement value\n           NO\
_REMOTE_PDB_DIR..........[Prevents the program fro\
m searching remote file: faster]\n           PDB_D\
IR....................[Indicates where PDB file mu\
st be fetched (localy)]\n\n	 PROBLEMS: please cont\
act cedric.notredame\\@europe.com\\n\";\n	 exit ($\
EXIT_SUCCESS);\n}\n\n$np=0;\n$n_para=$#ARGV;\n$mod\
el=1;\n$pdb_dir=$ENV{'PDB_DIR'};if ($pdb_dir){$pdb\
_dir.=\"/\";}\n$debug=$ENV{'DEBUG_EXTRACT_FROM_PDB\
'};\n\n$no_remote_pdb_dir=$ENV{NO_REMOTE_PDB_DIR};\
\n$HOME=$ENV{'HOME'};\nif ( $ENV{CACHE_4_TCOFFEE})\
\n{$cache=$ENV{CACHE_4_TCOFFEE};}\nelse\n{\n    $c\
ache=\"$HOME/.t_coffee/cache/\";\n}\n\n   \n$netad\
dress=\"http://www.rcsb.org/pdb/files/%%.pdb.gz\";\
\n$netcompression_pg=\"gunzip\";\n$netcompression=\
\"gz\";\n\nforeach ($np=0; $np<=$n_para; $np++)\n \
 {        \n    $value=$ARGV[$np];\n   \n    if  (\
$np==0 && !($value=~/^-.*/))\n      { \n	$pdb_file\
= $ARGV[$np];\n      }\n    elsif ( !($value=~/^-.\
*/))\n      {\n	print \"@ARGV\";\n	die;\n      } \\
n    \n    elsif ($value eq \"-nodiagnostic\"){$no\
diagnostic=1;}\n    elsif ($value eq \"-force\")\n\
      {\n	$force_pdb=1;\n      }\n    elsif ($valu\
e eq \"-force_name\")\n      {\n	$force_name=$ARGV\
[++$np];\n	$force_pdb=1;\n      }\n    \n    elsif\
 ($value eq \"-is_pdb_name\")\n      {\n	$pdb_file\
= $ARGV[++$np];	\n	$is_pdb_name=1;	\n      } \n   \
 elsif ($value eq \"-is_released_pdb_name\")\n    \
  {\n	$pdb_file= $ARGV[++$np];	\n	$is_released_pdb\
_name=1;\n      }\n    elsif ($value eq \"-model_t\
ype\")\n      {\n	$pdb_file= $ARGV[++$np];	\n	$mod\
el_type=1;\n      }\n    elsif ($value eq \"-debug\
\")\n{\n	$debug=1;\n}\n    elsif ($value eq \"-get\
_pdb_chains\")\n{\n	$pdb_file= $ARGV[++$np];\n	$ge\
t_pdb_chains=1;\n}\n    elsif ($value eq \"-get_pd\
b_ligands\")\n{\n	$get_pdb_ligands=1;\n}\n    \n  \
  elsif ($value eq \"-get_pdb_id\")\n{\n	$pdb_file\
= $ARGV[++$np];\n	$get_pdb_id=1;\n	\n}\n    \n    \
elsif ( $value eq \"-get_fugue_name\")\n{\n	$pdb_f\
ile= $ARGV[++$np];\n	$get_fugue_name=1;\n}\n    el\
sif ( $value eq \"-infile\")\n{\n       $pdb_file=\
 $ARGV[++$np];\n} \n    elsif ($value eq \"-netfil\
e\")\n{\n	$netfile=1;\n	if ( !($ARGV[$np+1]=~/^-.*\
/)){$pdb_file= $ARGV[++$np];}\n}\n    elsif (  $va\
lue eq \"-num\")\n{\n       $numbering= $ARGV[++$n\
p];\n}\n    elsif (  $value eq \"-num_out\")\n{\n \
      $numbering_out= $ARGV[++$np];\n}\n    elsif \
( $value eq \"-netaddress\")\n{\n	$netadress=$ARGV\
[++$np];\n}\n     \n    elsif ( $value eq \"-netco\
mpression\")\n{\n	 $netcompression=$ARGV[++$np];\n\
}\n    elsif ( $value eq \"-pdb_dir\")\n{\n	 if ( \
!($ARGV[$np+1]=~/^-.*/)){$pdb_dir= \"$ARGV[++$np]/\
\";}\n}\n     elsif ( $value eq \"-no_remote_pdb_d\
ir\")\n{\n	$no_remote_pdb_dir=1;\n	if ( !($ARGV[$n\
p+1]=~/^-.*/)){$pdb_dir= \"$ARGV[++$np]/\";}\n}\n \
   elsif ( $value eq \"-cache\")\n{\n	$cache=$ARGV\
[++$np];\n}\n    \n    elsif ($value eq \"-netcomp\
ression_pg\")\n{\n	  $netcompression_pg=$ARGV[++$n\
p];\n}\n     elsif ($value eq \"-mode\")\n{\n     \
  $MODE=$ARGV[++$np];\n}\n\n    elsif ( $value eq \
\"-model\")\n{\n       $model= $ARGV[++$np];\n}\n \
   elsif ($value eq \"-seq_field\" )\n{\n       $s\
eq_field= $ARGV[++$np];\n}   \n    elsif ($value e\
q \"-coor\" )\n{\n       $start= $ARGV[++$np];\n  \
\n       if (($ARGV[$np+1] eq \"\") ||($ARGV[$np+1\
]=~/^-.*/)){$end=\"*\";} \n       else {$end=   $A\
RGV[++$np];}     \n       $coor_set=1;\n}\n    els\
if ($value eq \"-delete\" )\n{\n       $delete_sta\
rt= $ARGV[++$np];\n       $delete_end= $ARGV[++$np\
];\n       $delete_set=1;\n}\n    elsif  ($value e\
q \"-code\")\n{\n       $code= $ARGV[++$np];\n}\n \
   elsif  ($value eq \"-no_hetatm\")\n{\n       $n\
o_hetatm=1;\n}\n    elsif ($value eq \"-chain\")\n\
{\n       while (!($ARGV[$np+1] eq \"\") &&!($ARGV\
[$np+1]=~/^-.*/))\n{\n	      ++$np;\n	      @c_cha\
in=(@chain,  $ARGV[$np]);\n	      $hc_chain{$ARGV[\
$np]}=$#c_chain+1;\n}           \n}\n    elsif ($v\
alue eq \"-atom\")\n{\n\n       while (!($ARGV[$np\
+1] eq \"\") && !($ARGV[$np+1]=~/^-.*/))\n{\n	    \
  ++$np;\n	      $atom[$n_atom++]=  $ARGV[$np];\n	\
      $atom_list{$ARGV[$np]}=1;	      \n} \n      \
 \n}\n    elsif ( $value eq \"-unfold\")\n{\n	$unf\
old=1;\n}\n    elsif ($value eq \"-seq\" ||$value \
eq \"-fasta\" )\n{\n       $MODE=\"fasta\";\n}\n  \
  elsif ( $value eq \"-version\")\n{\n	print STDER\
R  \"\\nextract_from_pdb: Version $VersionTag\\n\"\
;\n	&myexit ($EXIT_SUCCESS);\n}\n    elsif ( $valu\
e eq \"-ligand\")\n{\n	while (!($ARGV[$np+1] eq \"\
\") && !($ARGV[$np+1]=~/^-.*/))\n{\n	    ++$np;\n	\
    $ligand=1;\n	    $ligand_list{$ARGV[$np]}=1;	 \
     \n} \n	$hc_chain{'LIGAND'}=1;\n}\n    elsif (\
 $value eq \"-ligand_only\")\n{\n	$ligand_only=1;\\
n}\n}\nif ( $debug)\n{\n    print STDERR \"\\n[DEB\
UG:extract_from_pdb] NO_REMOTE_PDB_DIR: $no_remote\
_pdb_dir\\n\";\n    print STDERR \"\\n[DEBUG:extra\
ct_from_pdb] PDB_DIR: $pdb_dir\\n\";\n}\n\n\nif ( \
$is_pdb_name)\n  {\n    if (&remote_is_pdb_name($p\
db_file))\n      {\n	print \"1\";\n      }\n    el\
se\n      {\n	print \"0\";\n      }\n    exit ($EX\
IT_SUCCESS);\n  }\n\nif ( $is_released_pdb_name)\n\
  {\n    \n    if (&is_released($pdb_file))\n     \
 {\n	print \"1\";\n      }\n    else\n      {\n	pr\
int \"0\";\n      }\n    exit ($EXIT_SUCCESS);\n  \
}\nif ($model_type)\n  {\n   \n    printf \"%s\", \
&pdb2model_type($pdb_file);\n    exit ($EXIT_SUCCE\
SS);\n    \n  }\n    \n\nif (!$force_name)\n{\n   \
 $pdb_file=~/([^\\/]*)$/;\n    $force_name=$1;\n}\\
n\n$local_pdb_file=$pdb_file;\n\nif ( $debug){prin\
t STDERR \"\\n[DEBUG: extract_from_pdb] Scan For $\
local_pdb_file\\n\";}\n\n$mem=$no_remote_pdb_dir;\\
n$no_remote_pdb_dir=1;\n$tmp_pdb_file=get_pdb_file\
 ($local_pdb_file);\n\nif ( !-e $tmp_pdb_file || $\
tmp_pdb_file eq \"\")\n  {\n    $local_pdb_file=$p\
db_file;\n    ($local_pdb_file, $suffix_chain)=&pd\
b_name2name_and_chain($local_pdb_file);\n\n    if \
($local_pdb_file)\n      {\n	if ( $debug){print ST\
DERR \"\\nSplit $pdb_file into $local_pdb_file and\
 $suffix_chain \\n\";}\n	$tmp_pdb_file=get_pdb_fil\
e ($local_pdb_file);\n	if ( $tmp_pdb_file ne \"\")\
\n	  {\n	    @c_chain=();\n	    @c_chain=($suffix_\
chain);\n	    %hc_chain=();\n	    $hc_chain{$suffi\
x_chain}=1;\n	  }\n      }\n  }\n\n$no_remote_pdb_\
dir=$mem;\nif ($no_remote_pdb_dir==0)\n  {\n    \n\
    if ( !-e $tmp_pdb_file || $tmp_pdb_file eq \"\\
")\n      {\n	\n	$local_pdb_file=$pdb_file;\n	($lo\
cal_pdb_file, $suffix_chain)=&pdb_name2name_and_ch\
ain($local_pdb_file);\n	if ($local_pdb_file)\n	  {\
\n	    \n	    if ( $debug){print STDERR \"\\nSplit\
 $pdb_file into $local_pdb_file and $suffix_chain \
\\n\";}\n	    \n	    $tmp_pdb_file=get_pdb_file ($\
local_pdb_file);    \n	    \n	    if ( $tmp_pdb_fi\
le ne \"\")\n	      {\n		@c_chain=();\n		@c_chain=\
($suffix_chain);\n		%hc_chain=();\n		$hc_chain{$su\
ffix_chain}=1;\n	      }\n	  }\n      }\n  }\n\nif\
 ( $debug){print STDERR \"\\n$pdb_file copied into\
 ##$tmp_pdb_file##\\n\";}\n\nif ( !-e $tmp_pdb_fil\
e || $tmp_pdb_file eq \"\")\n{\n	if ($is_pdb_name)\
\n{\n	    print \"0\\n\"; exit ($EXIT_SUCCESS);\n}\
\n	else\n{\n  \n	    print \"\\nEXTRACT_FROM_PDB: \
NO RESULT for $pdb_file\\n\";\n	    &myexit ($EXIT\
_SUCCESS);	\n}\n}\n\n\n\n\n%molecule_type=&pdbfile\
2chaintype($tmp_pdb_file);\nif ( $debug){print STD\
ERR \"\\n\\tSequence Type determined\\n\";}\n\n$pd\
b_id=&get_pdb_id ($tmp_pdb_file);\nif ( $debug){pr\
int STDERR \"\\n\\tID: $pdb_id (for $tmp_pdb_file)\
\\n\";}\n\nif ( $pdb_id eq \"\"){$pdb_id=$force_na\
me;}\n\n@f_chain=&get_chain_list ($tmp_pdb_file);\\
nif ( $debug){print STDERR \"\\n\\tChain_list:@f_c\
hain\\n\";}\n\nif ( $get_pdb_chains)\n{\n    print\
 \"@f_chain\\n\";\n    &myexit ($EXIT_SUCCESS);\n}\
\nif ( $get_pdb_ligands)\n{\n    %complete_ligand_\
list=&get_ligand_list ($tmp_pdb_file);\n    print \
$complete_ligand_list{\"result\"};\n    &myexit ($\
EXIT_SUCCESS);\n}\n\nelsif ( $get_pdb_id ||$get_fu\
gue_name )\n{\n    if    (@c_chain && $c_chain[0] \
eq \"FIRST\"){$pdb_id=$pdb_id.$f_chain[0];}\n    e\
lsif (@c_chain && $c_chain[0] ne \" \"){$pdb_id=$p\
db_id.$c_chain[0];}\n    \n    print \"$pdb_id\\n\\
";\n    &myexit ($EXIT_SUCCESS);\n    \n}\nelsif (\
 $is_pdb_name)\n{\n    printf \"1\\n\";\n    &myex\
it ($EXIT_SUCCESS);\n}\n\n\n\n$structure_file=vtmp\
nam();\n\nif ( $debug){print STDERR \"\\n\\tCheck_\
point #1: $tmp_pdb_file  $structure_file\\n\";}\n\\
n$INFILE=vfopen (\"$tmp_pdb_file\", \"r\"); \n$TMP\
=vfopen (\"$structure_file\", \"w\");\n\n$print_mo\
del=1;\n$in_model=0;\n\nif ( $debug){print STDERR \
\"\\n\\tCheck_point #2\\n\";}\nwhile ( <$INFILE>)\\
n{\n  my $first_model=0;\n  $line=$_;\n\n  if ( !$\
first_model && ($line =~/^MODEL\\s*(\\d*)/))\n    \
{\n      $first_model=$1;\n      if ($model==1){$m\
odel=$first_model;}\n    }\n  \n  if (($line =~/^M\
ODEL\\s*(\\d*)/))\n    {\n      if ($1==$model)\n	\
{\n	  $in_model=1;\n	  $print_model=1;\n	  $is_nmr\
=1;\n	}\n      elsif ( $in_model==0)\n	{\n	  $prin\
t_model=0;\n	}\n      elsif ( $in_model==1)\n	{\n	\
  last;\n	}\n    }\n  if ($print_model){print $TMP\
 $line;}  \n}\nclose ($TMP);\nclose ($INFILE);\n\n\
if ( $debug){print STDERR \"\\n\\tCheck_point #3\\\
n\";}	\n\n  if ($numbering eq \"\"){$numbering=\"a\
bsolute\";}\n  if ($numbering_out eq \"\"){$number\
ing_out=\"new\";}\n\n  if ( $delete_set && $coor_s\
et) {die \"-delete and -coor are mutually exclusiv\
e, sorry\\n\";}\n  if ( $n_atom==0){$atom_list[$n_\
atom++]=\"ALL\";$atom_list{$atom_list[0]}=1;}\n  i\
f ( $seq_field eq \"\"){$seq_field=\"ATOM\";}\n  \\
n  if ( $MODE eq \"\"){$MODE=\"pdb\";}\n  elsif ( \
$MODE eq \"simple\" && $code==0){$code=1;}\n\n  if\
 ( $code==0){$code=3;}\n\n\nif ($f_chain[0] eq \" \
\"){$hc_chain{' '}=1;$c_chain[0]=\" \";}\nelsif (!\
@c_chain){$hc_chain{FIRST}=1;$c_chain[0]=\"FIRST\"\
;}#make sure the first chain is taken by default\n\
\nif    ($hc_chain{ALL}) \n{\n      @c_chain=@f_ch\
ain;\n      foreach $e (@c_chain){$hc_chain{$e}=1;\
}\n}\nelsif($hc_chain{FIRST})\n{\n	@c_chain=($f_ch\
ain[0]);\n	$hc_chain{$f_chain[0]}=1;\n}\n\n\n$MAIN\
_HOM_CODE=&get_main_hom_code ($structure_file);\n$\
INFILE=vfopen ($structure_file, \"r\");\n\n\nif ( \
$MODE eq \"raw_pdb\" || $MODE eq \"raw\")\n{\n    \
while (<$INFILE>)\n{	print \"$_\";}\n    close ( $\
INFILE);\n    &myexit($EXIT_SUCCESS);\n}    \nif (\
 $MODE eq \"raw4fugue\" )\n{\n    while (<$INFILE>\
)\n{	\n	$l=$_;\n	if ($l=~/^SEQRES/)\n{\n	    \n	  \
  $c= substr($l,11,1);\n	    if ($hc_chain {$c}){p\
rint \"$l\";}\n}\n	elsif ( $l=~/^ATOM/)\n{\n	    $\
c=substr($l,21,1);\n	    if ($hc_chain {$c}){print\
 \"$l\";}\n}\n}\n    close ( $INFILE);\n    &myexi\
t($EXIT_SUCCESS);\n}    \n\nif ( $MODE eq \"pdb\")\
\n{\n\n    $read_header=0;\n    while (<$INFILE>) \
\n{\n	    $line=$_;\n	    if    ($line =~ /^HEADER\
/){print \"$line\";$read_header=1;}\n}\n    close \
($INFILE);\n\n    if (!$read_header)\n{\n	print \"\
HEADER    UNKNOWN                                 \
00-JAN-00   $force_name\\n\";\n}\n\n    $INFILE=vf\
open ($structure_file, \"r\");\n    \n    print \"\
COMPND   1 CHAIN:\";\n    $last=pop(@c_chain);\n  \
  foreach $c ( @c_chain){ print \" $c,\";}\n    if\
 ( $last eq \" \"){print \" NULL;\\n\";}\n    else\
 \n{\n      print \" $last;\\n\";\n}\n    @c_chain\
=(@c_chain, $last);\n    \n    print \"REMARK Outp\
ut of the program extract_from_pdb (Version $Versi\
onTag)\\n\";\n    print \"REMARK Legal PDB format \
not Guaranteed\\n\";\n    print \"REMARK This form\
at is not meant to be used in place of the PDB for\
mat\\n\";\n    print \"REMARK The header refers to\
 the original entry\\n\";\n    print \"REMARK The \
sequence from the original file has been taken in \
the field: $seq_field\\n\";\n    print \"REMARK ex\
tract_from_pdb, 2001, 2002, 2003, 2004, 2005 2006 \
(c) CNRS and Cedric Notredame\\n\";   \n    if ( $\
coor_set)\n{\n       print \"REMARK Partial chain:\
 Start $start End $end\\n\";\n}\n    if ( $is_nmr)\
\n{\n       print \"REMARK NMR structure: MODEL $m\
odel\\n\";\n}\n   \n    if ( $n_atom!=0)\n{\n     \
  print \"REMARK Contains Coordinates of: \";\n   \
    foreach $a (@atom){print \"$a \";}\n       pri\
nt \"\\n\";\n}  \n}\n\n\n\n\nmy $residue_index = -\
999;\nmy $old_c = \"TemporaryChain\";\n\nwhile (<$\
INFILE>) \n{\n	$line=$_;\n\n\n	if ($line =~ /^SEQR\
ES/)\n{\n\n		@field=/(\\S*)\\s*/g;\n\n		$c= substr\
($_,11,1);\n\n		\n		$l=$#field;\n		for ($a=4; $a<$\
#field ;)\n{\n			if (!$onelett{$molecule_type{$c}}\
->{$field[$a]})\n{\n				splice @field, $a, 1;\n}\n\
			else \n{\n				$a++;\n}\n}\n	\n		if ( $c ne $in_\
chain)\n{\n			$pdb_chain_list[$n_pdb_chains]=$c;\n\
			$pdb_chain_len [$n_pdb_chains]=$len;\n			$in_ch\
ain=$c;\n			$n_pdb_chains++;\n}\n	\n		for ( $a=4; \
$a<$#field;$a++)\n{\n			@{$complete_seq{$c}}->[$co\
mplete_seq_len{$c}++]=$field[$a];\n}\n}\n    elsif\
 ( $line=~/^ATOM/ || ($line=~/^HETATM/ && &is_aa(s\
ubstr($line,17,3),substr($line,21,1)) && !$no_heta\
tm))\n{\n\n	 \n    $RAW_AT_ID=$AT_ID=substr($line,\
12,4);\n	$RES_ID=&is_aa(substr($line,17,3),substr(\
$line,21,1));\n	$CHAIN=substr($line,21,1);\n\n    \
$RES_NO=substr($line,22,4);\n	$HOM_CODE=substr ($l\
ine, 26, 1);\n	$TEMP=substr($line,60,6);\n	\n	$TEM\
P=~s/\\s//g;\n        $AT_ID=~s/\\s//g;\n	$RES_ID=\
~s/\\s//g;\n        $RES_NO=~s/\\s//g;\n		\n	if ( \
$HOM_CODE ne $MAIN_HOM_CODE){next;}\n	elsif ( $alr\
eady_read2{$CHAIN}{$RES_ID}{$AT_ID}{$RES_NO}){next\
;}\n	else{$already_read2{$CHAIN}{$RES_ID}{$AT_ID}{\
$RES_NO}=1;}\n	\n	\n	if ($coor_set && $numbering e\
q \"file\" && $residue_index ne $RES_NO)\n{\n	    \
\n	    if ( $RES_NO<=$start){$real_start{$CHAIN}++\
;}\n	    if ( $RES_NO<=$end){$real_end{$CHAIN}++;}\
\n}\n	elsif ($numbering eq \"absolute\")\n{\n	    \
$real_start{$CHAIN}=$start;\n	    $real_end{$CHAIN\
}=$end;\n}\n\n        $KEY=\"ALL\";\n        if ( \
$CHAIN ne $in_atom_chain)\n{\n	    \n	  $pdb_atom_\
chain_list[$n_pdb_atom_chains]=$c;\n	  $pdb_atom_c\
hain_len [$n_pdb_atom_chains]=$len;\n	  $in_atom_c\
hain=$c;\n	  $n_pdb_atom_chains++;\n}\n	\n	if ( $r\
esidue_index ne $RES_NO)\n{\n	     $residue_index \
= $RES_NO;\n	     @{$atom_seq{$CHAIN}}->[$atom_seq\
_len{$CHAIN}++]=$RES_ID;;\n}\n}\n\n}\nclose ($INFI\
LE);\n\n\n\n\n\n\n$INFILE=vfopen ($structure_file,\
 \"r\");\nforeach $c (@c_chain)\n{\n\n	if    ( $se\
q_field eq \"SEQRES\"){@pdb_seq=@{$complete_seq{$c\
}};}\n	elsif ( $seq_field eq \"ATOM\")  {@pdb_seq=\
@{$atom_seq{$c}};}\n	\n\n	$full_length=$l=$#pdb_se\
q+1;\n		\n	if ( $real_end{$c}==\"*\"){$real_end{$c\
}=$full_length;}\n	if ( $coor_set)\n{	   \n\n	   i\
f ( $real_end{$c} < $l){splice @pdb_seq, $real_end\
{$c}, $l;}\n	   if ( $real_start{$c} < $l){splice \
@pdb_seq, 0, $real_start{$c}-1;}	  	   \n	   $l=$#\
pdb_seq;\n}\n\n	elsif ( $delete_set)\n{\n	   splic\
e @pdb_seq, $delete_start, $delete_end-$delete_sta\
rt+1;\n	   $l=$#pdb_seq;\n}\n	\n	$new_fasta_name=\\
"$pdb_id$c\";\n	if ( $coor_set)\n{\n	   if ( $n_pd\
b_chains==0){$new_fasta_name=\"$new_fasta_name$c\"\
;}\n	   $new_fasta_name= $new_fasta_name.\"\\_$sta\
rt\\_$end\";\n}\n	   \n	if ( $MODE eq \"pdb\")\n{\\
n	   $nl=1;\n	   $n=0;\n	   \n	   foreach $res ( @\
pdb_seq)\n		{\n		if ( !$n)\n		{\n		\n		 printf \"S\
EQRES %3d %1s %4d  \", $nl,$c, $l;\n		 $nl++;\n	}\\
n	     $res=~s/\\s//g;\n	     \n	     if ($code==1\
){ printf \"%3s \",$onelett{$molecule_type{$c}}->{\
$res};}\n	     elsif  ($code==3){ printf \"%3s \",\
$res};\n	     \n	     $n++;		  \n	     if ( $n==13\
){$n=0;print \"\\n\";}\n}\n	  if ( $n!=0){print \"\
\\n\"; $n=0;}\n}\n	elsif ( $MODE eq \"simple\")\n{\
\n	  print \"# SIMPLE_PDB_FORMAT\\n\";\n	  if ( $n\
ew_fasta_name eq \" \"){$new_fasta_name=\"dummy_na\
me\";}\n	  print \">$new_fasta_name\\n\";\n\n	  fo\
reach $res ( @pdb_seq)\n{\n	      print \"$onelett\
{$molecule_type{$c}}->{$res}\";\n}\n	  print \"\\n\
\";\n}\n	elsif ( $MODE eq \"fasta\")\n{\n	  $n=0;\\
n	  print \">$new_fasta_name\\n\";\n	  \n	  foreac\
h $res ( @pdb_seq)\n{\n\n	      print \"$onelett{$\
molecule_type{$c}}->{$res}\";\n              $n++;\
\n	      if ( $n==60){print \"\\n\"; $n=0;}\n}\n	 \
 print \"\\n\"; \n}\n}\n\nif ( $MODE eq \"fasta\")\
\n{\n     &myexit($EXIT_SUCCESS);\n  \n}\n\n  \n  \
$charcount=0;\n  $inchain=\"BEGIN\";\n  $n=0;\n  w\
hile (<$INFILE>) \n{\n    $line=$_;\n     \n    if\
 ($line =~/^ATOM/  ||  ($line=~/^HETATM/))\n{\n	$l\
ine_header=\"UNKNWN\";\n	$RES_ID=substr($line,17,3\
);\n	$chain = substr($line,21,1);\n\n	if ($line =~\
/^ATOM/)\n{\n	    $line_header=\"ATOM\";\n	    $RE\
S_ID=(&is_aa($RES_ID,$chain))?&is_aa($RES_ID,$chai\
n):$RES_ID;\n}\n	elsif ($line=~/^HETATM/ && ($liga\
nd_list {$RES_ID} ||$ligand_list {'ALL'} || !&is_a\
a($RES_ID,$chain)))\n{\n	    $line_header=\"HETATM\
\";\n}\n	elsif ($line=~/^HETATM/ && (&is_aa($RES_I\
D,$chain) && !$no_hetatm))\n{\n	    $line_header=\\
"ATOM\";\n	    $RES_ID=&is_aa($RES_ID,$chain);\n}\\
n	else\n{\n	    next;\n}\n\n	\n\n	$X=substr($line,\
30,8);     \n	$Y=substr($line,38,8);\n	$Z=substr($\
line,46,8);\n	$TEMP=substr($line,60,6);\n	\n	$RAW_\
AT_ID=$AT_ID=substr($line,12,4);\n	$CHAIN=substr($\
line,21,1);\n	$RES_NO=substr($line,22,4);\n	$HOM_C\
ODE=substr ($line, 26, 1);\n	\n	$X=~s/\\s//g;\n	$Y\
=~s/\\s//g;\n	$Z=~s/\\s//g;\n	$TEMP=~s/\\s//g;\n	\\
n	$AT_ID=~s/\\s//g;\n	$RES_ID=~s/\\s//g;\n	$RES_NO\
=~s/\\s//g;\n\n	\n	if ( $HOM_CODE ne $MAIN_HOM_COD\
E){next;}\n	elsif ( $already_read{$CHAIN}{$RES_ID}\
{$AT_ID}{$RES_NO}){next;}\n	else{$already_read{$CH\
AIN}{$RES_ID}{$AT_ID}{$RES_NO}=1;}\n	\n	$KEY=\"ALL\
\";\n\n      	if ( $RES_NO ==0){$start_at_zero=1;}\
\n\n	$RES_NO+=$start_at_zero;    \n	\n	if ( $curre\
nt_chain ne $CHAIN)\n{\n	    $current_chain=$CHAIN\
;\n	    $pos=$current_residue=0;\n	    $offset=($c\
oor_set)?($real_start{$CHAIN}-1):0;\n	    if    ( \
$seq_field eq \"SEQRES\"){@ref_seq=@{$complete_seq\
{$CHAIN}};}\n	    elsif ( $seq_field eq \"ATOM\") \
 {@ref_seq=@{$atom_seq{$CHAIN}};}\n}\n	\n	if ($cur\
rent_residue != $RES_NO)\n{\n	    $current_residue\
=$RES_NO;\n	    if    ( $seq_field eq \"SEQRES\"){\
$pos=$current_residue;}\n	    elsif ( $seq_field e\
q \"ATOM\"){$pos++;}\n}\n	\n	\n	if ($n_atom==0 || \
$atom_list{$AT_ID}==1 || $atom_list{$KEY}==1)\n{ 	\
\n	    \n	    $do_it=(!@c_chain || $hc_chain{$CHAI\
N} ||$hc_chain{'LIGAND'} );\n	    \n	    $do_it= (\
$do_it==1) && ($coor_set==0 ||($pos>=$real_start{$\
CHAIN} && $pos<=$real_end{$CHAIN}));\n	    $do_it=\
 ($do_it==1) && ($delete_set==0 || $pos<$delete_st\
art ||$pos>$delete_end );\n	    if ($ligand==0 && \
$line_header eq \"HETATM\" ){$do_it=0;}\n	    if (\
$ligand_only==1 && $line_header eq \"ATOM\" ){$do_\
it=0;}\n	    if ($ligand==1 && $line_header eq \"H\
ETATM\" && $ligand_list{$RES_ID}==0 && $ligand_lis\
t{\"ALL\"}==0){$do_it=0;} \n	    \n	    \n	    if \
( $do_it)\n{\n		$n++;\n		$out_pos=$pos;\n		\n	    \
   if ( $delete_set)\n{\n		  if ( $out_pos< $delet\
e_start){;}\n		  else {$offset=$delete_end-$delete\
_start;}\n}       \n	       \n	       if ( $number\
ing_out eq \"new\"){$out_pos-=$offset;}\n	       e\
lsif ( $numbering_out eq \"old\"){$out_pos=$RES_NO\
;}\n	       \n       \n	       \n	       if ( $cod\
e==1){$RES_ID=$onelett{$molecule_type{$c}}->{$RES_\
ID};}\n	    \n	       if ($unfold)\n{\n		   $unfol\
ded_x+=5;\n		   $X=$unfolded_x;\n		   $Y=0;\n		   \
$Z=0;\n		   $float=1;\n}\n	       else\n{\n		   $f\
loat=3;\n}\n\n	       if ( $MODE eq \"pdb\")\n{\n	\
	   printf \"%-6s%5d %-4s %3s %s%4d    %8.3f%8.3f%\
8.3f  1.00 %5.2f\\n\",$line_header, $n, $RAW_AT_ID\
,$RES_ID,$CHAIN,$out_pos, $X, $Y, $Z,$TEMP;		  \n}\
\n	       elsif ( $MODE eq \"simple\")\n{\n		    i\
f ( $RES_ID eq \"\"){$RES_ID=\"X\";}\n		  printf \\
"%-6s %5s %s %2s %4d    %8.3f %8.3f %8.3f\\n\",$li\
ne_header, $AT_ID, $RES_ID,($CHAIN eq\"\" || $CHAI\
N eq \" \")?\"A\":$CHAIN,$out_pos, $X, $Y, $Z,$TEM\
P;\n}\n\n}\n}\n}\n}\nprint \"\\n\";\nclose($INFILE\
);\n\n\nif ( $error ne \"\") \n{$error=$error.\"\\\
nDiagnostic:    SEQRES and the residues in ATOM ar\
e probably Incompatible\\n\";\n    $error=$error. \
 \"Recomendation: Rerun with '-fix 1' in order to \
ignore the SEQRES sequences\\n\";\n}\nif (!$nodiag\
nostic){print STDERR $error;}\n&myexit ( $EXIT_SUC\
CESS);\n\nsub is_released \n  {\n    my ($r);\n   \
 my $in=@_[0];\n    my $name=&remote_is_pdb_name (\
$in);\n    my $hold=&remote_is_on_hold($in);\n    \
\n    $r=($name && !$hold)?1:0;\n    return $r;\n \
 }\nsub remote_is_pdb_name \n{\n    my $in=@_[0];\\
n    my ($ref_file, $pdb);\n    my ($value,$value1\
,$value2);\n\n    if ( $in=~/[^\\w\\d\\:\\_]/){ret\
urn 0;}\n    $ref_file=\"$cache/pdb_entry_type.txt\
\";\n    \n    if ( !-e $ref_file || (-M $ref_file\
)>2 || -z $ref_file)\n      {\n	&url2file(\"ftp://\
ftp.wwpdb.org/pub/pdb/derived_data/pdb_entry_type.\
txt\", $ref_file);\n      }\n    $pdb=substr ($in,\
0, 4);\n    chomp(($value1=`grep -c $pdb $ref_file\
`));\n    $pdb=lc($pdb);\n    chomp(($value2=`grep\
 -c $pdb $ref_file`));\n    $value=($value1 || $va\
lue2)?1:0;\n    $value=($value>0)?1:0;\n    \n    \
return $value;\n  }\n\nsub pdb2model_type\n{\n    \
my $in=@_[0];\n    my ($ref_file, $pdb);\n    my (\
$value, $ret);\n\n    if ( $in=~/[^\\w\\d\\:\\_]/)\
{return 0;}\n    $ref_file=\"$cache/pdb_entry_type\
.txt\";\n    \n    if ( !-e $ref_file || (-M $ref_\
file)>2 || -z $ref_file)\n      {\n	&url2file(\"ft\
p://ftp.wwpdb.org/pub/pdb/derived_data/pdb_entry_t\
ype.txt\", $ref_file);\n      }\n    $pdb=substr (\
$in,0, 4);\n    $pdb=lc($pdb);\n    \n    chomp(($\
value=`grep $pdb $ref_file`));\n    \n    $value=~\
/^\\S+\\s+\\S+\\s+(\\S+)/;\n    $ret=$1;\n    if (\
 $ret eq\"\"){return \"UNKNOWN\";}\n    \n    retu\
rn $ret;\n  }\nsub remote_is_on_hold\n  {\n    my \
$in=@_[0];\n    my ($ref_file, $pdb);\n    my ($va\
lue1, $value2,$value);\n    \n    if ( $in=~/[^\\w\
\\d\\:\\_]/){return 0;}\n    $ref_file=\"$cache/un\
released.xml\";\n    \n    if ( !-e $ref_file || (\
-M $ref_file)>2 || -z $ref_file)\n      {\n	&url2f\
ile(\"http://www.rcsb.org/pdb/rest/getUnreleased\"\
,$ref_file);\n      }\n    \n    $pdb=substr ($in,\
0, 4);\n    chomp(($value1=`grep -c $pdb $ref_file\
`));\n    $pdb=lc($pdb);\n    chomp(($value2=`grep\
 -c $pdb $ref_file`));\n    $value=($value1 || $va\
lue2)?1:0;\n    $value=($value>0)?1:0;\n    return\
 $value;\n  }\nsub is_pdb_file\n{\n    my @arg=@_;\
\n\n    if ( !-e $arg[0]){return 0;}\n    \n    $F\
=vfopen ($arg[0], \"r\");\n    while ( <$F>)\n{\n	\
if (/^HEADER/)\n{\n	    close $F;\n	    return 1;\\
n}\n	elsif ( /^SEQRES/)\n{\n	    close $F;\n	    r\
eturn 1;\n}\n	elsif ( /^ATOM/)\n{\n	    close $F;\\
n	    return 1;\n}\n}\n    return 0;\n}\nsub get_p\
db_id\n{\n    my $header_file=@_[0];\n    my $id;\\
n    my $F= new FileHandle;\n    \n    \n    $F=vf\
open (\"$header_file\", \"r\");\n\n    while ( <$F\
>)\n      {\n	if ( /HEADER/)\n	  {\n	    if ($debu\
g){print \"$_\";}\n	    $id=substr($_,62,4 );\n	  \
  return $id;\n	  }\n      }\n    close ($F);\n   \
 \n    return \"\";\n}\n\nsub get_ligand_list\n{\n\
    my $pdb_file=@_[0];\n    my $chain;\n    my $l\
igand;\n    my %complete_ligand_list;\n    \n\n   \
 $F=vfopen ($pdb_file, \"r\");\n    while ( <$F>)\\
n{\n	if ( /^HETATM/)\n{\n	    $line=$_;\n	    $cha\
in=substr($line,21,1);\n	    $ligand=substr($line,\
17,3);\n	    \n	    if (!$complete_ligand_list{$ch\
ain}{$ligand})\n{\n		\n		$complete_ligand_list{\"r\
esult\"}.=\"CHAIN $chain LIGAND $ligand\\n\";\n		$\
complete_ligand_list{$chain}{$ligand}=1;\n}\n}\n}\\
n    close ($F);\n    return %complete_ligand_list\
;\n}\n\nsub get_chain_list \n{\n    my $header_fil\
e;\n    my @chain_list;\n    my @list;\n    my $n_\
chains;\n    my %chain_hasch;\n    my $pdb_file=@_\
[0];\n    my $c;\n    my %hasch;\n    my $chain;\n\
  \n    \n    $F=vfopen ($pdb_file, \"r\");\n    w\
hile ( <$F>)\n{\n\n\n	if (/SEQRES\\s+\\d+\\s+(\\S+\
)/)\n	  {\n	    $chain = substr($_,11,1);$chain=~s\
/\\s//g;if ( $chain eq \"\"){$chain=\" \";}\n	    \
if (!$hasch{$chain}){$hasch{$chain}=1;push @chain_\
list, $chain;}\n	  }\n	if (/^ATOM/ || /^HETATM/)\n\
	  {\n	    $chain = substr($_,21,1); $chain=~s/\\s\
//g;if ( $chain eq \"\"){$chain=\" \";}\n	    if (\
!$hasch{$chain}){$hasch{$chain}=1;push @chain_list\
, $chain;}\n	  }\n      }\n\n\nclose ($F);\nif (!@\
chain_list)\n  {\n    @chain_list=(\"A\");\n  }\n\\
n\nreturn @chain_list;\n}\n\nsub token_is_in_list\\
n{\n\n    my @list=@_;\n    my $a;\n    \n    for \
($a=1; $a<=$#list; $a++)\n{\n	if ( $list[$a] eq $l\
ist[0]){return $a;}\n}\n}\n\nsub pdb_name2name_and\
_chain \n{\n    my $pdb_file=@_[0];\n    my $pdb_f\
ile_in;\n    my @array;\n    my $chain;\n    my $c\
;\n\n    $pdb_file_in=$pdb_file;\n\n    $pdb_file=\
~/^(.{4})/;$pdb_id=$1;\n    @array=($pdb_file=~/([\
\\w])/g);\n  \n  \n    $chain=uc ($array[4]);\n   \
 $chain=($chain eq \"\")?\"FIRST\":$chain;\n    \n\
    return ( $pdb_id, $chain);\n\n    if ( $#array\
==3){return ($pdb_id, \"FIRST\");}\n    elsif ( $#\
array<4){ return ($pdb_id, \"\");}\n    else {retu\
rn ( $pdb_id, $chain);}\n      \n    \n    \n}\nsu\
b get_main_hom_code \n{\n    my $pdb_file=@_[0];\n\
    my %hom, $n, $best, $best_h;\n    open (F, $pd\
b_file);\n    while (<F>)\n{\n	if ( $_=~/^ATOM/)\n\
{\n	    $h=substr ($_,26, 1);\n	    $n=++$hom{$h};\
\n	    if ($n>$best)\n{\n		$best=$n;\n		$best_h=$h\
;\n}\n}\n}\n    close (F);\n    return $best_h;\n}\
\n\n\nsub get_pdb_file \n{\n    my ($pdb_file_in)=\
(@_);\n    my $result;\n    my @letter;\n    my @c\
hain;\n    my $v;\n    my $pdb_file=$pdb_file_in;\\
n\n    $pdb_file=($pdb_file_in=~/\\S+_S_(\\S+)/)?$\
1:$pdb_file_in;\n    \n    if ($no_remote_pdb_dir=\
=0)\n      {\n	$no_remote_pdb_dir=1;\n	$result=get\
_pdb_file3 ($pdb_file);\n	$no_remote_pdb_dir=0;\n	\
if ( $result){return $result;}\n	else\n	  {\n	    \
\n	    lc ($pdb_file);\n	    $result=get_pdb_file3\
($pdb_file);\n	    return  $result;\n	  }\n      }\
\n    else\n      {\n	return get_pdb_file3 ($pdb_f\
ile);\n      }\n    \n  }\n\nsub get_pdb_file3 \n{\
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
my $value=&is_released ($pdb_file);\n    my $retur\
n_value=\"\";\n    if ($value==1)\n      {\n	\n	&d\
ebug_print (\"\\n*********************************\
********************\\nTry Remote Access for $pdb_\
file\");\n	$tmp_pdb_file=vtmpnam();\n	$netcommand=\
$netaddress;\n	$netcommand=~s/%%/$pdb_file/g;\n	&u\
rl2file(\"$netcommand\", \"$tmp_pdb_file.$netcompr\
ession\");\n	&debug_print(\"\\nREMOTE: $netcommand\
\\n\");\n	\n	$compressed_tmp_file_name=\"$tmp_pdb_\
file.$netcompression\";\n	\n	if ($netcompression &\
& -B $compressed_tmp_file_name)\n	  {\n	    my $r;\
\n	    &debug_print (\"\\n\\tFile Found Remotely\"\
);\n	    if (($r=safe_system ( \"$netcompression_p\
g $compressed_tmp_file_name\")!=$EXIT_SUCCESS) && \
$attempts<5)\n	      {\n		&debug_print (\"\\n\\tPr\
oper Download Failed Try again\");\n		unlink $comp\
ressed_tmp_file_name;\n		print \"\\nFailed to Down\
load $compressed_tmp_file_name. New Attempt $attem\
pt/5\\n\";\n		return &get_pdb_file0($pdb_file, $at\
tempt+1);\n	      }\n	    elsif ($r== $EXIT_SUCCES\
S)\n	      {\n		&debug_print (\"\\n\\tProper Downl\
oad Succeeded \");\n		$return_value=$tmp_pdb_file;\
\n	      }\n	    else\n	      {\n		&debug_print (\\
"\\n\\tProper Download Failed \");\n		&debug_print\
 (\"\\nFile Not Found Remotely\");\n		unlink $comp\
ressed_tmp_file_name;\n	      }\n	  }\n	else\n	  {\
\n\n	    &debug_print (\"\\nFile Not Found Remotel\
y\");\n	    unlink $compressed_tmp_file_name;\n	  \
}\n	#Update cache if required\n	if ($cache ne \"no\
\" && $cache ne \"update\" && -e $return_value)\n	\
  {\n	    `cp $return_value $cache/$pdb_file.pdb`;\
\n	    #`t_coffee -other_pg clean_cache.pl -file $\
pdb_file.pdb -dir $cache`;\n	  }\n      }\n    &de\
bug_print (\"\\nRemote Download Finished\");\n    \
return $return_value;\n  }\nreturn \"\";\n}\n\nsub\
 check_pdb_file4compression \n{\n    my $file=@_[0\
];\n    my $tmp;\n    my $r;\n    \n    $tmp=&vtmp\
nam();\n    if (-e $tmp){unlink $tmp;}\n    \n    \
$file=~s/\\/\\//\\//g;\n    if    (-B $file && ($f\
ile=~/\\.Z/)) {`cp $file $tmp.Z`;`rm $tmp`;`gunzip\
 $tmp.Z $SILENT`;$r=$tmp;}\n    elsif (-B $file &&\
 ($file=~/\\.gz/)){`cp $file $tmp.gz`;`gunzip $tmp\
.gz $SILENT`;return $r=$tmp;}\n    elsif (-B $file\
 ){`cp $file $tmp.gz`;`gunzip $tmp.gz $SILENT`;$r=\
$tmp;}\n    elsif ( -e $file ) {$r= $file;}\n    e\
lsif ( -e \"$file.gz\" ){ `cp $file.gz $tmp.gz`;`g\
unzip     $tmp.gz $SILENT`;$r=$tmp;}    \n    elsi\
f ( -e \"$file.Z\") {`cp $file.Z  $tmp.Z`; `gunzip\
 $tmp.Z $SILENT`;$r=$tmp;}\n    else  {$r= $file;}\
\n\n    if ( -e \"$tmp.Z\"){unlink \"$tmp.Z\";}\n \
   if ( -e \"$tmp.gz\"){unlink \"$tmp.gz\";}\n    \
\n    return $r;\n    \n}\n\n\n\n\n\n    \n\n\n\n\\
n\n\n\nsub vfopen \n{\n    my $file=@_[0];\n    my\
 $mode=@_[1];\n    my $tmp;\n    my $F = new FileH\
andle;\n    \n    \n    $tmp=$file;\n	\n    \n    \
if ( $mode eq \"r\" && !-e $file){ myexit(flush_er\
ror (\"Cannot open file $file\"));}\n    elsif ($m\
ode eq \"w\"){$tmp=\">$file\";}\n    elsif ($mode \
eq \"a\"){$tmp=\">>$file\";}\n    \n    \n    open\
 ($F,$tmp);\n    return $F;\n}\nsub debug_print\n{\
\n    my $message =@_[0];\n    if ($debug){print S\
TDERR \"NO_REMOTE_PDB_DIR: $no_remote_pdb_dir - $m\
essage [DEBUG:extract_from_pdb]\";}\n    return;\n\
}\nsub is_aa \n{\n    my ($aa, $chain) =@_;\n\n   \
 my $one;\n    my $trhee;\n    \n    if ( $onelett\
{$molecule_type{$chain}}->{$aa} eq 'X' || !$onelet\
t{$molecule_type{$chain}}->{$aa} ){return '';}\n  \
  else\n      {\n	$one=$onelett{$molecule_type{$ch\
ain}}->{$aa};\n\n	$three=$threelett{$molecule_type\
{$chain}}->{$one};\n	\n\n	return $three;\n      }\\
n  }\n\n\n\n\n\nsub url2file\n{\n    my ($address,\
 $out, $wget_arg, $curl_arg)=(@_);\n    my ($pg, $\
flag, $r, $arg, $count);\n    \n    if (!$CONFIGUR\
ATION){&check_configuration (\"wget\", \"INTERNET\\
", \"gzip\");$CONFIGURATION=1;}\n    \n    if (&pg\
_is_installed (\"wget\"))   {$pg=\"wget\"; $flag=\\
"-O\";$arg=$wget_arg;}\n    elsif (&pg_is_installe\
d (\"curl\")){$pg=\"curl\"; $flag=\"-o\";$arg=$cur\
l_arg;}\n    return safe_system (\"$pg $flag$out $\
address >/dev/null 2>/dev/null\");\n\n}\n\n\n\n\ns\
ub pdbfile2chaintype\n  {\n    my $file=@_[0];\n  \
  my %ct;\n    my $F;\n    \n    $F=vfopen ($file,\
 \"r\");\n    while (<$F>)\n      {\n	my $line=$_;\
\n	if ($line =~/^ATOM/)\n	  {\n	    my $C=substr($\
line,21,1);\n	    if (!$ct{$C})\n	      {	\n		my $\
r=substr($line,17,3);\n		$r=~s/\\s+//;\n		if (leng\
th ($r)==1){$ct{$C}=\"R\";}\n		elsif (length ($r)=\
=2){$ct{$C}=\"D\";}\n		elsif (length ($r)==3){$ct{\
$C}=\"P\";}\n		else \n		  {\n		    myexit(flush_er\
ror(\"ERROR: Could not read RES_ID field in file $\
file\"));\n		  }\n	      }\n	  }\n      }\n    clo\
se ($F);\n    return %ct;\n  }\n   \n    \n\n\n\ns\
ub fill_threelett_RNA\n{\n\n	my %threelett=(\n	'A'\
, '  A',\n	'T', '  T',\n	'U', '  U',\n	'C', '  C',\
\n	'G', '  G',\n	'I', '  I', #Inosine\n	);\n	\n	re\
turn %threelett;\n\n}\n\n\nsub fill_onelett_RNA\n{\
\n	my   %onelett=(\n	'  A' => 'A',\n	'  T' => 'T',\
\n	'  U' => 'U',\n	'  C' => 'C',\n	'  G' => 'G',\n\
	'CSL' => 'X',\n	'UMS' => 'X',\n	'  I' => 'I',\n	'\
A' => 'A',\n	'T' => 'T',\n	'U' => 'U',\n	'C' => 'C\
',\n	'G' => 'G',\n	'I' => 'I',\n	);\n\n	return %on\
elett;\n\n}\n\n\nsub fill_onelett_DNA\n{\n	my   %o\
nelett=(\n	' DA', 'A',\n	' DT', 'T',\n	' DC', 'C',\
\n	' DG', 'G',\n	'DA', 'A',\n	'DT', 'T',\n	'DC', '\
C',\n	'DG', 'G',\n	);\n\n	return %onelett;\n\n}\n\\
nsub fill_threelett_DNA\n{\n\n	my %threelett=(\n	'\
A', ' DA',\n	'T', ' DT',\n	'C', ' DC',\n	'G', ' DG\
',\n	);\n\n	return %threelett;\n\n}\n\n\n\n\nsub f\
ill_threelett_prot\n{  \n  my %threelett;\n\n  %th\
reelett=(\n'A', 'ALA',\n'C', 'CYS',\n'D', 'ASP',\n\
'E', 'GLU',\n'F', 'PHE',\n'G', 'GLY',\n'H', 'HIS',\
\n'I', 'ILE',\n'K', 'LYS',\n'L', 'LEU',\n'N', 'ASN\
',\n'M', 'MET',\n'P', 'PRO',\n'Q', 'GLN',\n'R', 'A\
RG',\n'S', 'SER',\n'T', 'THR',\n'V', 'VAL',\n'W', \
'TRP',\n'Y', 'TYR',\n);\n\nreturn %threelett;\n\n\\
n}\n\nsub fill_onelett_prot\n{\n    my %onelett;\n\
    \n    %onelett=(\n\n'10A', 'X',\n'11O', 'X',\n\
'12A', 'X',\n'13P', 'X',\n'13R', 'X',\n'13S', 'X',\
\n'14W', 'X',\n'15P', 'X',\n'16A', 'X',\n'16G', 'X\
',\n'1AN', 'X',\n'1AP', 'X',\n'1AR', 'X',\n'1BH', \
'X',\n'1BO', 'X',\n'1C5', 'X',\n'1CU', 'X',\n'1DA'\
, 'X',\n'1GL', 'X',\n'1GN', 'X',\n'1IN', 'X',\n'1L\
U', 'L',\n'1MA', 'X',\n'1MC', 'X',\n'1MG', 'X',\n'\
1MZ', 'X',\n'1NA', 'X',\n'1NB', 'X',\n'1NI', 'X',\\
n'1PA', 'A',\n'1PC', 'X',\n'1PE', 'X',\n'1PG', 'X'\
,\n'1PI', 'A',\n'1PM', 'X',\n'1PN', 'X',\n'1PU', '\
X',\n'1PY', 'X',\n'1UN', 'X',\n'24T', 'X',\n'25T',\
 'X',\n'26P', 'X',\n'2AB', 'X',\n'2AM', 'X',\n'2AN\
', 'X',\n'2AP', 'X',\n'2AR', 'X',\n'2AS', 'D',\n'2\
BL', 'X',\n'2BM', 'X',\n'2CP', 'X',\n'2DA', 'X',\n\
'2DG', 'X',\n'2DP', 'X',\n'2DT', 'X',\n'2EP', 'X',\
\n'2EZ', 'X',\n'2FG', 'X',\n'2FL', 'X',\n'2FP', 'X\
',\n'2FU', 'X',\n'2GL', 'X',\n'2GP', 'X',\n'2HP', \
'X',\n'2IB', 'X',\n'2IP', 'X',\n'2LU', 'L',\n'2MA'\
, 'X',\n'2MD', 'X',\n'2ME', 'X',\n'2MG', 'X',\n'2M\
L', 'L',\n'2MO', 'X',\n'2MR', 'R',\n'2MU', 'X',\n'\
2MZ', 'X',\n'2NO', 'X',\n'2NP', 'X',\n'2OG', 'X',\\
n'2PA', 'X',\n'2PC', 'X',\n'2PE', 'X',\n'2PG', 'X'\
,\n'2PH', 'X',\n'2PI', 'X',\n'2PL', 'X',\n'2PP', '\
X',\n'2PU', 'X',\n'2SI', 'X',\n'2TB', 'X',\n'34C',\
 'X',\n'35G', 'X',\n'3AA', 'X',\n'3AD', 'X',\n'3AH\
', 'H',\n'3AN', 'X',\n'3AP', 'X',\n'3AT', 'X',\n'3\
BT', 'X',\n'3CH', 'X',\n'3CN', 'X',\n'3CO', 'X',\n\
'3CP', 'X',\n'3DR', 'X',\n'3EP', 'X',\n'3FM', 'X',\
\n'3GA', 'X',\n'3GP', 'X',\n'3HB', 'X',\n'3HC', 'X\
',\n'3HP', 'X',\n'3IB', 'X',\n'3ID', 'X',\n'3IN', \
'X',\n'3MA', 'X',\n'3MB', 'X',\n'3MC', 'X',\n'3MD'\
, 'D',\n'3MF', 'X',\n'3MP', 'X',\n'3MT', 'X',\n'3O\
L', 'X',\n'3PA', 'X',\n'3PG', 'X',\n'3PO', 'X',\n'\
3PP', 'X',\n'3PY', 'X',\n'49A', 'X',\n'4AB', 'X',\\
n'4AM', 'X',\n'4AN', 'X',\n'4AP', 'X',\n'4BA', 'X'\
,\n'4BT', 'X',\n'4CA', 'X',\n'4CO', 'X',\n'4HP', '\
X',\n'4IP', 'X',\n'4MO', 'X',\n'4MV', 'X',\n'4MZ',\
 'X',\n'4NC', 'X',\n'4NP', 'X',\n'4OX', 'X',\n'4PB\
', 'X',\n'4PN', 'X',\n'4PP', 'X',\n'4SC', 'X',\n'4\
SU', 'X',\n'4TB', 'X',\n'55C', 'X',\n'5AD', 'X',\n\
'5AN', 'X',\n'5AT', 'X',\n'5CM', 'X',\n'5GP', 'X',\
\n'5HP', 'E',\n'5HT', 'X',\n'5IT', 'X',\n'5IU', 'X\
',\n'5MB', 'X',\n'5MC', 'X',\n'5MD', 'X',\n'5MP', \
'X',\n'5MU', 'X',\n'5NC', 'X',\n'5OB', 'X',\n'5PA'\
, 'X',\n'5PV', 'X',\n'6AB', 'X',\n'6CT', 'X',\n'6H\
A', 'X',\n'6HC', 'X',\n'6HG', 'X',\n'6HT', 'X',\n'\
6IN', 'X',\n'6MO', 'X',\n'6MP', 'X',\n'6PG', 'X',\\
n'6WO', 'X',\n'70U', 'X',\n'7DG', 'X',\n'7HP', 'X'\
,\n'7I2', 'X',\n'7MG', 'X',\n'7MQ', 'X',\n'7NI', '\
X',\n'87Y', 'X',\n'8AD', 'X',\n'8BR', 'X',\n'8IG',\
 'X',\n'8IN', 'X',\n'8OG', 'X',\n'95A', 'X',\n'9AD\
', 'X',\n'9AM', 'X',\n'9AP', 'X',\n'9DG', 'X',\n'9\
DI', 'X',\n'9HX', 'X',\n'9OH', 'X',\n'9TA', 'X',\n\
'A12', 'X',\n'A15', 'X',\n'A23', 'X',\n'A24', 'X',\
\n'A26', 'X',\n'A2G', 'X',\n'A2P', 'X',\n'A32', 'X\
',\n'A3P', 'X',\n'A4P', 'X',\n'A5P', 'X',\n'A70', \
'X',\n'A76', 'X',\n'A77', 'X',\n'A78', 'X',\n'A79'\
, 'X',\n'A80', 'X',\n'A85', 'X',\n'A88', 'X',\n'A9\
A', 'X',\n'AA3', 'X',\n'AA4', 'X',\n'AA6', 'X',\n'\
AAA', 'X',\n'AAB', 'X',\n'AAC', 'X',\n'AAE', 'X',\\
n'AAG', 'R',\n'AAH', 'X',\n'AAM', 'X',\n'AAN', 'X'\
,\n'AAP', 'X',\n'AAR', 'R',\n'AAS', 'X',\n'AAT', '\
X',\n'ABA', 'X',\n'ABC', 'X',\n'ABD', 'X',\n'ABE',\
 'X',\n'ABH', 'X',\n'ABI', 'X',\n'ABK', 'X',\n'ABM\
', 'X',\n'ABN', 'X',\n'ABP', 'X',\n'ABR', 'X',\n'A\
BS', 'X',\n'ABU', 'X',\n'AC1', 'X',\n'AC2', 'X',\n\
'ACA', 'X',\n'ACB', 'D',\n'ACC', 'C',\n'ACD', 'X',\
\n'ACE', 'X',\n'ACH', 'X',\n'ACI', 'X',\n'ACL', 'R\
',\n'ACM', 'X',\n'ACN', 'X',\n'ACO', 'X',\n'ACP', \
'X',\n'ACQ', 'X',\n'ACR', 'X',\n'ACS', 'X',\n'ACT'\
, 'X',\n'ACV', 'V',\n'ACX', 'X',\n'ACY', 'X',\n'AD\
2', 'X',\n'AD3', 'X',\n'ADC', 'X',\n'ADD', 'X',\n'\
ADE', 'X',\n'ADH', 'X',\n'ADI', 'X',\n'ADM', 'X',\\
n'ADN', 'X',\n'ADP', 'X',\n'ADQ', 'X',\n'ADR', 'X'\
,\n'ADS', 'X',\n'ADT', 'X',\n'ADU', 'X',\n'ADW', '\
X',\n'ADX', 'X',\n'AE2', 'X',\n'AEA', 'X',\n'AEB',\
 'X',\n'AEI', 'D',\n'AEN', 'X',\n'AET', 'T',\n'AF1\
', 'X',\n'AF3', 'X',\n'AFA', 'D',\n'AFP', 'X',\n'A\
G7', 'X',\n'AGB', 'X',\n'AGF', 'X',\n'AGL', 'X',\n\
'AGM', 'R',\n'AGN', 'X',\n'AGP', 'X',\n'AGS', 'X',\
\n'AGU', 'X',\n'AH0', 'X',\n'AH1', 'X',\n'AHA', 'X\
',\n'AHB', 'D',\n'AHC', 'X',\n'AHF', 'X',\n'AHG', \
'X',\n'AHH', 'X',\n'AHM', 'X',\n'AHO', 'X',\n'AHP'\
, 'X',\n'AHS', 'X',\n'AHT', 'Y',\n'AHU', 'X',\n'AH\
X', 'X',\n'AI1', 'X',\n'AI2', 'X',\n'AIB', 'X',\n'\
AIC', 'X',\n'AIM', 'X',\n'AIP', 'X',\n'AIQ', 'X',\\
n'AIR', 'X',\n'AJ3', 'X',\n'AKB', 'X',\n'AKG', 'X'\
,\n'AKR', 'X',\n'AL1', 'X',\n'AL2', 'X',\n'AL3', '\
X',\n'AL4', 'X',\n'AL5', 'X',\n'AL6', 'X',\n'AL7',\
 'X',\n'AL8', 'X',\n'AL9', 'X',\n'ALA', 'A',\n'ALB\
', 'X',\n'ALC', 'X',\n'ALD', 'L',\n'ALE', 'X',\n'A\
LF', 'X',\n'ALG', 'X',\n'ALL', 'X',\n'ALM', 'A',\n\
'ALN', 'A',\n'ALO', 'T',\n'ALP', 'X',\n'ALQ', 'X',\
\n'ALR', 'X',\n'ALS', 'X',\n'ALT', 'A',\n'ALY', 'K\
',\n'ALZ', 'X',\n'AMA', 'X',\n'AMB', 'X',\n'AMC', \
'X',\n'AMD', 'X',\n'AMG', 'X',\n'AMH', 'X',\n'AMI'\
, 'X',\n'AML', 'X',\n'AMN', 'X',\n'AMO', 'X',\n'AM\
P', 'X',\n'AMQ', 'X',\n'AMR', 'X',\n'AMS', 'X',\n'\
AMT', 'X',\n'AMU', 'X',\n'AMW', 'X',\n'AMX', 'X',\\
n'AMY', 'X',\n'ANA', 'X',\n'ANB', 'X',\n'ANC', 'X'\
,\n'AND', 'X',\n'ANE', 'X',\n'ANI', 'X',\n'ANL', '\
X',\n'ANO', 'X',\n'ANP', 'X',\n'ANS', 'X',\n'ANT',\
 'X',\n'AOE', 'X',\n'AOP', 'X',\n'AP1', 'X',\n'AP2\
', 'X',\n'AP3', 'X',\n'AP4', 'X',\n'AP5', 'X',\n'A\
P6', 'X',\n'APA', 'X',\n'APB', 'X',\n'APC', 'X',\n\
'APE', 'F',\n'APF', 'X',\n'APG', 'X',\n'APH', 'A',\
\n'API', 'X',\n'APL', 'X',\n'APM', 'X',\n'APN', 'G\
',\n'APP', 'X',\n'APQ', 'X',\n'APR', 'X',\n'APS', \
'X',\n'APT', 'X',\n'APU', 'X',\n'APX', 'X',\n'APY'\
, 'X',\n'APZ', 'X',\n'AQS', 'X',\n'AR1', 'X',\n'AR\
2', 'X',\n'ARA', 'X',\n'ARB', 'X',\n'ARC', 'X',\n'\
ARD', 'X',\n'ARG', 'R',\n'ARH', 'X',\n'ARI', 'X',\\
n'ARM', 'R',\n'ARN', 'X',\n'ARO', 'R',\n'ARP', 'X'\
,\n'ARQ', 'X',\n'ARS', 'X',\n'AS1', 'R',\n'AS2', '\
X',\n'ASA', 'D',\n'ASB', 'D',\n'ASC', 'X',\n'ASD',\
 'X',\n'ASE', 'X',\n'ASF', 'X',\n'ASI', 'X',\n'ASK\
', 'D',\n'ASL', 'X',\n'ASM', 'N',\n'ASO', 'X',\n'A\
SP', 'D',\n'ASQ', 'X',\n'ASU', 'X',\n'ATA', 'X',\n\
'ATC', 'X',\n'ATD', 'X',\n'ATF', 'X',\n'ATG', 'X',\
\n'ATH', 'X',\n'ATM', 'X',\n'ATO', 'X',\n'ATP', 'X\
',\n'ATQ', 'X',\n'ATR', 'X',\n'ATT', 'X',\n'ATY', \
'X',\n'ATZ', 'X',\n'AUC', 'X',\n'AUR', 'X',\n'AVG'\
, 'X',\n'AXP', 'X',\n'AYA', 'A',\n'AZ2', 'X',\n'AZ\
A', 'X',\n'AZC', 'X',\n'AZD', 'X',\n'AZE', 'X',\n'\
AZI', 'X',\n'AZL', 'X',\n'AZM', 'X',\n'AZR', 'X',\\
n'AZT', 'X',\n'B12', 'X',\n'B1F', 'F',\n'B2A', 'A'\
,\n'B2F', 'F',\n'B2I', 'I',\n'B2V', 'V',\n'B3I', '\
X',\n'B3P', 'X',\n'B7G', 'X',\n'B96', 'X',\n'B9A',\
 'X',\n'BA1', 'X',\n'BAA', 'X',\n'BAB', 'X',\n'BAC\
', 'X',\n'BAF', 'X',\n'BAH', 'X',\n'BAI', 'X',\n'B\
AK', 'X',\n'BAL', 'A',\n'BAM', 'X',\n'BAO', 'X',\n\
'BAP', 'X',\n'BAR', 'X',\n'BAS', 'X',\n'BAT', 'F',\
\n'BAY', 'X',\n'BAZ', 'X',\n'BB1', 'X',\n'BB2', 'X\
',\n'BBA', 'X',\n'BBH', 'X',\n'BBS', 'X',\n'BBT', \
'X',\n'BBZ', 'X',\n'BCA', 'X',\n'BCB', 'X',\n'BCC'\
, 'X',\n'BCD', 'X',\n'BCL', 'X',\n'BCN', 'X',\n'BC\
R', 'X',\n'BCS', 'C',\n'BCT', 'X',\n'BCY', 'X',\n'\
BCZ', 'X',\n'BDA', 'X',\n'BDG', 'X',\n'BDK', 'X',\\
n'BDM', 'X',\n'BDN', 'X',\n'BDS', 'X',\n'BE1', 'X'\
,\n'BE2', 'X',\n'BEA', 'X',\n'BEF', 'X',\n'BEN', '\
X',\n'BEO', 'X',\n'BEP', 'X',\n'BER', 'X',\n'BES',\
 'X',\n'BET', 'X',\n'BEZ', 'X',\n'BF2', 'X',\n'BFA\
', 'X',\n'BFD', 'X',\n'BFP', 'X',\n'BFS', 'X',\n'B\
FU', 'X',\n'BG6', 'X',\n'BGF', 'X',\n'BGG', 'X',\n\
'BGL', 'X',\n'BGN', 'X',\n'BGP', 'X',\n'BGX', 'X',\
\n'BH4', 'X',\n'BHA', 'X',\n'BHC', 'X',\n'BHD', 'D\
',\n'BHO', 'X',\n'BHS', 'X',\n'BIC', 'X',\n'BIN', \
'X',\n'BIO', 'X',\n'BIP', 'X',\n'BIS', 'X',\n'BIZ'\
, 'X',\n'BJH', 'X',\n'BJI', 'X',\n'BJP', 'X',\n'BL\
A', 'X',\n'BLB', 'X',\n'BLE', 'L',\n'BLG', 'P',\n'\
BLI', 'X',\n'BLM', 'X',\n'BLV', 'X',\n'BLY', 'K',\\
n'BM1', 'X',\n'BM2', 'X',\n'BM5', 'X',\n'BM9', 'X'\
,\n'BMA', 'X',\n'BMD', 'X',\n'BME', 'X',\n'BMP', '\
X',\n'BMQ', 'X',\n'BMS', 'X',\n'BMT', 'T',\n'BMU',\
 'X',\n'BMY', 'X',\n'BMZ', 'X',\n'BNA', 'X',\n'BNG\
', 'X',\n'BNI', 'X',\n'BNN', 'F',\n'BNO', 'L',\n'B\
NS', 'X',\n'BNZ', 'X',\n'BO3', 'X',\n'BO4', 'X',\n\
'BOC', 'X',\n'BOG', 'X',\n'BOM', 'X',\n'BOT', 'X',\
\n'BOX', 'X',\n'BOZ', 'X',\n'BPA', 'X',\n'BPB', 'X\
',\n'BPD', 'X',\n'BPG', 'X',\n'BPH', 'X',\n'BPI', \
'X',\n'BPJ', 'X',\n'BPM', 'X',\n'BPN', 'X',\n'BPO'\
, 'X',\n'BPP', 'X',\n'BPT', 'X',\n'BPY', 'X',\n'BR\
B', 'X',\n'BRC', 'X',\n'BRE', 'X',\n'BRI', 'X',\n'\
BRL', 'X',\n'BRM', 'X',\n'BRN', 'X',\n'BRO', 'X',\\
n'BRS', 'X',\n'BRU', 'X',\n'BRZ', 'X',\n'BSB', 'X'\
,\n'BSI', 'X',\n'BSP', 'X',\n'BT1', 'X',\n'BT2', '\
X',\n'BT3', 'X',\n'BTA', 'L',\n'BTB', 'X',\n'BTC',\
 'C',\n'BTD', 'X',\n'BTN', 'X',\n'BTP', 'X',\n'BTR\
', 'W',\n'BU1', 'X',\n'BUA', 'X',\n'BUB', 'X',\n'B\
UC', 'X',\n'BUG', 'X',\n'BUL', 'X',\n'BUM', 'X',\n\
'BUQ', 'X',\n'BUT', 'X',\n'BVD', 'X',\n'BX3', 'X',\
\n'BYS', 'X',\n'BZ1', 'X',\n'BZA', 'X',\n'BZB', 'X\
',\n'BZC', 'X',\n'BZD', 'X',\n'BZF', 'X',\n'BZI', \
'X',\n'BZM', 'X',\n'BZO', 'X',\n'BZP', 'X',\n'BZQ'\
, 'X',\n'BZS', 'X',\n'BZT', 'X',\n'C02', 'X',\n'C1\
1', 'X',\n'C1O', 'X',\n'C20', 'X',\n'C24', 'X',\n'\
C2F', 'X',\n'C2O', 'X',\n'C2P', 'X',\n'C3M', 'X',\\
n'C3P', 'X',\n'C3X', 'X',\n'C48', 'X',\n'C4M', 'X'\
,\n'C4X', 'X',\n'C5C', 'X',\n'C5M', 'X',\n'C5P', '\
X',\n'C5X', 'X',\n'C60', 'X',\n'C6C', 'X',\n'C6M',\
 'X',\n'C78', 'X',\n'C8E', 'X',\n'CA3', 'X',\n'CA5\
', 'X',\n'CAA', 'X',\n'CAB', 'X',\n'CAC', 'X',\n'C\
AD', 'X',\n'CAF', 'C',\n'CAG', 'X',\n'CAH', 'X',\n\
'CAL', 'X',\n'CAM', 'X',\n'CAN', 'X',\n'CAO', 'X',\
\n'CAP', 'X',\n'CAQ', 'X',\n'CAR', 'X',\n'CAS', 'C\
',\n'CAT', 'X',\n'CAV', 'X',\n'CAY', 'C',\n'CAZ', \
'X',\n'CB3', 'X',\n'CB4', 'X',\n'CBA', 'X',\n'CBD'\
, 'X',\n'CBG', 'X',\n'CBI', 'X',\n'CBL', 'X',\n'CB\
M', 'X',\n'CBN', 'X',\n'CBO', 'X',\n'CBP', 'X',\n'\
CBS', 'X',\n'CBX', 'X',\n'CBZ', 'X',\n'CC0', 'X',\\
n'CC1', 'X',\n'CCC', 'X',\n'CCH', 'X',\n'CCI', 'X'\
,\n'CCM', 'X',\n'CCN', 'X',\n'CCO', 'X',\n'CCP', '\
X',\n'CCR', 'X',\n'CCS', 'C',\n'CCV', 'X',\n'CCY',\
 'X',\n'CD1', 'X',\n'CDC', 'X',\n'CDE', 'X',\n'CDF\
', 'X',\n'CDI', 'X',\n'CDL', 'X',\n'CDM', 'X',\n'C\
DP', 'X',\n'CDR', 'X',\n'CDU', 'X',\n'CE1', 'X',\n\
'CEA', 'C',\n'CEB', 'X',\n'CEC', 'X',\n'CED', 'X',\
\n'CEF', 'X',\n'CEH', 'X',\n'CEM', 'X',\n'CEO', 'X\
',\n'CEP', 'X',\n'CEQ', 'X',\n'CER', 'X',\n'CES', \
'G',\n'CET', 'X',\n'CFC', 'X',\n'CFF', 'X',\n'CFM'\
, 'X',\n'CFO', 'X',\n'CFP', 'X',\n'CFS', 'X',\n'CF\
X', 'X',\n'CGN', 'X',\n'CGP', 'X',\n'CGS', 'X',\n'\
CGU', 'E',\n'CH2', 'X',\n'CH3', 'X',\n'CHA', 'X',\\
n'CHB', 'X',\n'CHD', 'X',\n'CHF', 'X',\n'CHG', 'G'\
,\n'CHI', 'X',\n'CHN', 'X',\n'CHO', 'X',\n'CHP', '\
G',\n'CHR', 'X',\n'CHS', 'F',\n'CHT', 'X',\n'CHX',\
 'X',\n'CIC', 'X',\n'CIN', 'X',\n'CIP', 'X',\n'CIR\
', 'X',\n'CIT', 'X',\n'CIU', 'X',\n'CKI', 'X',\n'C\
L1', 'X',\n'CL2', 'X',\n'CLA', 'X',\n'CLB', 'A',\n\
'CLC', 'S',\n'CLD', 'A',\n'CLE', 'L',\n'CLF', 'X',\
\n'CLK', 'S',\n'CLL', 'X',\n'CLM', 'X',\n'CLN', 'X\
',\n'CLO', 'X',\n'CLP', 'X',\n'CLQ', 'X',\n'CLR', \
'X',\n'CLS', 'X',\n'CLT', 'X',\n'CLX', 'X',\n'CLY'\
, 'X',\n'CMA', 'R',\n'CMC', 'X',\n'CMD', 'X',\n'CM\
E', 'C',\n'CMG', 'X',\n'CMK', 'X',\n'CMN', 'X',\n'\
CMO', 'X',\n'CMP', 'X',\n'CMR', 'X',\n'CMS', 'X',\\
n'CMT', 'C',\n'CMX', 'X',\n'CNA', 'X',\n'CNC', 'X'\
,\n'CND', 'X',\n'CNH', 'X',\n'CNM', 'X',\n'CNN', '\
X',\n'CNO', 'X',\n'CNP', 'X',\n'CO2', 'X',\n'CO3',\
 'X',\n'CO5', 'X',\n'CO8', 'X',\n'COA', 'X',\n'COB\
', 'X',\n'COC', 'X',\n'COD', 'X',\n'COE', 'X',\n'C\
OF', 'X',\n'COH', 'X',\n'COI', 'X',\n'COJ', 'X',\n\
'COL', 'X',\n'COM', 'X',\n'CON', 'X',\n'COP', 'X',\
\n'COR', 'X',\n'COS', 'X',\n'COT', 'X',\n'COY', 'X\
',\n'CP1', 'G',\n'CP2', 'X',\n'CP4', 'X',\n'CPA', \
'X',\n'CPB', 'X',\n'CPC', 'X',\n'CPD', 'X',\n'CPG'\
, 'X',\n'CPH', 'X',\n'CPI', 'X',\n'CPM', 'X',\n'CP\
N', 'G',\n'CPO', 'X',\n'CPP', 'X',\n'CPQ', 'X',\n'\
CPR', 'X',\n'CPS', 'X',\n'CPT', 'X',\n'CPU', 'X',\\
n'CPV', 'X',\n'CPY', 'X',\n'CR1', 'X',\n'CR6', 'X'\
,\n'CRA', 'X',\n'CRB', 'X',\n'CRC', 'X',\n'CRG', '\
X',\n'CRH', 'X',\n'CRO', 'T',\n'CRP', 'X',\n'CRQ',\
 'X',\n'CRS', 'X',\n'CRT', 'X',\n'CRY', 'X',\n'CSA\
', 'C',\n'CSB', 'X',\n'CSD', 'C',\n'CSE', 'C',\n'C\
SH', 'X',\n'CSI', 'X',\n'CSN', 'X',\n'CSO', 'C',\n\
'CSP', 'C',\n'CSR', 'C',\n'CSS', 'C',\n'CST', 'X',\
\n'CSW', 'C',\n'CSX', 'C',\n'CSY', 'X',\n'CSZ', 'C\
',\n'CT3', 'X',\n'CTA', 'X',\n'CTB', 'X',\n'CTC', \
'X',\n'CTD', 'X',\n'CTH', 'T',\n'CTO', 'X',\n'CTP'\
, 'X',\n'CTR', 'X',\n'CTS', 'X',\n'CTT', 'X',\n'CT\
Y', 'X',\n'CTZ', 'X',\n'CU1', 'X',\n'CUA', 'X',\n'\
CUC', 'X',\n'CUL', 'X',\n'CUO', 'X',\n'CUZ', 'X',\\
n'CVI', 'X',\n'CXF', 'X',\n'CXL', 'X',\n'CXM', 'M'\
,\n'CXN', 'X',\n'CXP', 'X',\n'CXS', 'X',\n'CY1', '\
C',\n'CY3', 'X',\n'CYB', 'X',\n'CYC', 'X',\n'CYF',\
 'C',\n'CYG', 'C',\n'CYH', 'X',\n'CYL', 'X',\n'CYM\
', 'C',\n'CYN', 'X',\n'CYO', 'X',\n'CYP', 'X',\n'C\
YQ', 'C',\n'CYS', 'C',\n'CYU', 'X',\n'CYY', 'X',\n\
'CYZ', 'X',\n'CZH', 'X',\n'CZZ', 'C',\n'D12', 'X',\
\n'D13', 'X',\n'D16', 'X',\n'D18', 'X',\n'D19', 'X\
',\n'D1P', 'X',\n'D24', 'X',\n'D34', 'X',\n'D35', \
'X',\n'D4D', 'X',\n'D4T', 'X',\n'D6G', 'X',\n'DA2'\
, 'R',\n'DA3', 'X',\n'DA6', 'X',\n'DA7', 'X',\n'DA\
A', 'X',\n'DAB', 'X',\n'DAC', 'X',\n'DAD', 'X',\n'\
DAE', 'X',\n'DAF', 'X',\n'DAG', 'X',\n'DAH', 'A',\\
n'DAJ', 'X',\n'DAK', 'X',\n'DAL', 'A',\n'DAM', 'A'\
,\n'DAN', 'X',\n'DAO', 'X',\n'DAP', 'X',\n'DAQ', '\
X',\n'DAR', 'R',\n'DAS', 'D',\n'DAT', 'X',\n'DAU',\
 'X',\n'DAV', 'X',\n'DBA', 'X',\n'DBD', 'X',\n'DBF\
', 'X',\n'DBG', 'X',\n'DBI', 'X',\n'DBV', 'X',\n'D\
BY', 'Y',\n'DCA', 'X',\n'DCB', 'X',\n'DCE', 'X',\n\
'DCF', 'X',\n'DCG', 'X',\n'DCH', 'X',\n'DCI', 'I',\
\n'DCL', 'X',\n'DCM', 'X',\n'DCP', 'X',\n'DCS', 'X\
',\n'DCT', 'X',\n'DCY', 'C',\n'DCZ', 'X',\n'DDA', \
'X',\n'DDB', 'X',\n'DDC', 'X',\n'DDF', 'X',\n'DDG'\
, 'X',\n'DDH', 'X',\n'DDL', 'X',\n'DDM', 'X',\n'DD\
O', 'L',\n'DDP', 'X',\n'DDQ', 'X',\n'DDT', 'Y',\n'\
DDU', 'X',\n'DEA', 'X',\n'DEB', 'X',\n'DEC', 'X',\\
n'DEF', 'X',\n'DEL', 'X',\n'DEM', 'X',\n'DEN', 'X'\
,\n'DEP', 'X',\n'DEQ', 'X',\n'DES', 'X',\n'DET', '\
X',\n'DFC', 'X',\n'DFG', 'X',\n'DFI', 'X',\n'DFL',\
 'X',\n'DFO', 'X',\n'DFP', 'X',\n'DFR', 'X',\n'DFT\
', 'X',\n'DFV', 'X',\n'DFX', 'X',\n'DG2', 'X',\n'D\
G3', 'X',\n'DG6', 'X',\n'DGA', 'X',\n'DGD', 'X',\n\
'DGG', 'X',\n'DGL', 'E',\n'DGN', 'Q',\n'DGP', 'X',\
\n'DGT', 'X',\n'DGX', 'X',\n'DH2', 'X',\n'DHA', 'A\
',\n'DHB', 'X',\n'DHC', 'X',\n'DHD', 'X',\n'DHE', \
'X',\n'DHF', 'X',\n'DHG', 'X',\n'DHI', 'H',\n'DHL'\
, 'X',\n'DHM', 'X',\n'DHN', 'V',\n'DHP', 'X',\n'DH\
Q', 'X',\n'DHR', 'X',\n'DHS', 'X',\n'DHT', 'X',\n'\
DHU', 'X',\n'DHY', 'X',\n'DHZ', 'X',\n'DI2', 'X',\\
n'DI3', 'G',\n'DI4', 'X',\n'DI5', 'X',\n'DIA', 'X'\
,\n'DIC', 'X',\n'DIF', 'X',\n'DIG', 'X',\n'DII', '\
X',\n'DIL', 'I',\n'DIM', 'X',\n'DIO', 'X',\n'DIP',\
 'X',\n'DIQ', 'X',\n'DIS', 'X',\n'DIT', 'X',\n'DIV\
', 'V',\n'DIX', 'X',\n'DIY', 'X',\n'DKA', 'X',\n'D\
LA', 'X',\n'DLE', 'L',\n'DLF', 'X',\n'DLS', 'K',\n\
'DLY', 'K',\n'DM1', 'X',\n'DM2', 'X',\n'DM3', 'X',\
\n'DM4', 'X',\n'DM5', 'X',\n'DM6', 'X',\n'DM7', 'X\
',\n'DM8', 'X',\n'DM9', 'X',\n'DMA', 'X',\n'DMB', \
'X',\n'DMC', 'X',\n'DMD', 'X',\n'DME', 'X',\n'DMF'\
, 'X',\n'DMG', 'G',\n'DMH', 'N',\n'DMI', 'X',\n'DM\
J', 'X',\n'DML', 'X',\n'DMM', 'X',\n'DMN', 'X',\n'\
DMO', 'X',\n'DMP', 'X',\n'DMQ', 'X',\n'DMR', 'X',\\
n'DMS', 'X',\n'DMT', 'X',\n'DMV', 'X',\n'DMY', 'X'\
,\n'DNC', 'X',\n'DND', 'X',\n'DNH', 'X',\n'DNJ', '\
X',\n'DNN', 'X',\n'DNP', 'X',\n'DNQ', 'X',\n'DNR',\
 'X',\n'DO2', 'X',\n'DO3', 'X',\n'DOA', 'X',\n'DOB\
', 'X',\n'DOC', 'X',\n'DOH', 'D',\n'DOM', 'X',\n'D\
OS', 'X',\n'DOX', 'X',\n'DP5', 'X',\n'DP7', 'X',\n\
'DPA', 'X',\n'DPC', 'X',\n'DPD', 'X',\n'DPE', 'X',\
\n'DPG', 'X',\n'DPH', 'F',\n'DPM', 'X',\n'DPN', 'F\
',\n'DPO', 'X',\n'DPP', 'X',\n'DPR', 'P',\n'DPS', \
'X',\n'DPT', 'X',\n'DPX', 'X',\n'DPY', 'X',\n'DPZ'\
, 'X',\n'DQH', 'X',\n'DQN', 'X',\n'DR1', 'X',\n'DR\
B', 'X',\n'DRC', 'X',\n'DRI', 'X',\n'DRP', 'X',\n'\
DRT', 'X',\n'DRU', 'X',\n'DSA', 'X',\n'DSB', 'X',\\
n'DSC', 'X',\n'DSD', 'X',\n'DSE', 'S',\n'DSI', 'X'\
,\n'DSN', 'S',\n'DSP', 'D',\n'DSR', 'X',\n'DSS', '\
X',\n'DSX', 'X',\n'DSY', 'X',\n'DTB', 'X',\n'DTD',\
 'X',\n'DTH', 'T',\n'DTN', 'X',\n'DTO', 'X',\n'DTP\
', 'X',\n'DTQ', 'X',\n'DTR', 'W',\n'DTT', 'X',\n'D\
TY', 'Y',\n'DUD', 'X',\n'DUO', 'X',\n'DUR', 'X',\n\
'DUT', 'X',\n'DVA', 'V',\n'DVR', 'X',\n'DX9', 'X',\
\n'DXA', 'X',\n'DXB', 'X',\n'DXC', 'X',\n'DXG', 'X\
',\n'DXX', 'X',\n'DZF', 'X',\n'E09', 'X',\n'E20', \
'X',\n'E2P', 'X',\n'E3G', 'X',\n'E4N', 'X',\n'E4P'\
, 'X',\n'E64', 'X',\n'E6C', 'X',\n'E96', 'X',\n'E9\
7', 'X',\n'EA2', 'X',\n'EAA', 'X',\n'EAP', 'X',\n'\
EBP', 'X',\n'EBW', 'X',\n'ECO', 'X',\n'EDA', 'X',\\
n'EDC', 'X',\n'EDE', 'X',\n'EDO', 'X',\n'EDR', 'X'\
,\n'EEB', 'X',\n'EEE', 'X',\n'EFC', 'X',\n'EFZ', '\
X',\n'EG1', 'X',\n'EG2', 'X',\n'EG3', 'X',\n'EGC',\
 'X',\n'EGL', 'X',\n'EHP', 'A',\n'EIC', 'X',\n'EJT\
', 'X',\n'ELA', 'X',\n'EMB', 'X',\n'EMC', 'X',\n'E\
MD', 'X',\n'EMM', 'X',\n'EMO', 'X',\n'EMP', 'X',\n\
'EMR', 'X',\n'ENA', 'X',\n'ENC', 'X',\n'ENH', 'X',\
\n'ENO', 'X',\n'ENP', 'X',\n'EOA', 'X',\n'EOH', 'X\
',\n'EOT', 'X',\n'EOX', 'X',\n'EPA', 'X',\n'EPE', \
'X',\n'EPH', 'X',\n'EPI', 'X',\n'EPN', 'X',\n'EPO'\
, 'X',\n'EPT', 'X',\n'EPU', 'X',\n'EPX', 'X',\n'EP\
Y', 'X',\n'EQI', 'X',\n'EQP', 'X',\n'EQU', 'X',\n'\
ERG', 'X',\n'ERI', 'X',\n'ERY', 'X',\n'ESC', 'X',\\
n'ESD', 'X',\n'ESI', 'X',\n'ESO', 'X',\n'ESP', 'X'\
,\n'EST', 'X',\n'ESX', 'X',\n'ETA', 'X',\n'ETC', '\
X',\n'ETD', 'X',\n'ETF', 'X',\n'ETH', 'X',\n'ETI',\
 'X',\n'ETN', 'X',\n'ETO', 'X',\n'ETP', 'X',\n'ETR\
', 'X',\n'ETS', 'X',\n'ETY', 'X',\n'EU3', 'X',\n'E\
UG', 'X',\n'EYS', 'C',\n'F09', 'X',\n'F2B', 'X',\n\
'F3S', 'X',\n'F42', 'X',\n'F43', 'X',\n'F4S', 'X',\
\n'F6B', 'X',\n'F6P', 'X',\n'F89', 'X',\n'FA1', 'X\
',\n'FA5', 'F',\n'FAA', 'X',\n'FAB', 'X',\n'FAC', \
'X',\n'FAD', 'X',\n'FAF', 'X',\n'FAG', 'X',\n'FAM'\
, 'X',\n'FAR', 'X',\n'FAS', 'X',\n'FAT', 'X',\n'FB\
A', 'X',\n'FBE', 'X',\n'FBI', 'X',\n'FBP', 'X',\n'\
FBQ', 'X',\n'FBS', 'X',\n'FBT', 'X',\n'FBU', 'X',\\
n'FCA', 'X',\n'FCB', 'X',\n'FCI', 'X',\n'FCN', 'X'\
,\n'FCO', 'X',\n'FCR', 'X',\n'FCT', 'X',\n'FCX', '\
X',\n'FCY', 'C',\n'FD1', 'F',\n'FD2', 'F',\n'FD3',\
 'F',\n'FD4', 'F',\n'FDA', 'X',\n'FDC', 'X',\n'FDI\
', 'X',\n'FDP', 'X',\n'FDS', 'X',\n'FE2', 'X',\n'F\
EA', 'X',\n'FEL', 'X',\n'FEM', 'X',\n'FEN', 'X',\n\
'FEO', 'X',\n'FEP', 'X',\n'FER', 'X',\n'FES', 'X',\
\n'FFB', 'X',\n'FFC', 'X',\n'FFF', 'X',\n'FFO', 'X\
',\n'FGL', 'G',\n'FHB', 'X',\n'FHC', 'X',\n'FHP', \
'X',\n'FHU', 'X',\n'FID', 'X',\n'FII', 'X',\n'FIP'\
, 'X',\n'FK5', 'X',\n'FKA', 'X',\n'FKI', 'X',\n'FK\
P', 'X',\n'FL2', 'X',\n'FL9', 'X',\n'FLA', 'A',\n'\
FLC', 'X',\n'FLD', 'X',\n'FLE', 'L',\n'FLF', 'X',\\
n'FLO', 'X',\n'FLP', 'X',\n'FLT', 'Y',\n'FLU', 'X'\
,\n'FLX', 'X',\n'FM1', 'X',\n'FM2', 'X',\n'FMA', '\
X',\n'FMB', 'X',\n'FMC', 'X',\n'FME', 'M',\n'FMN',\
 'X',\n'FMP', 'X',\n'FMR', 'X',\n'FMS', 'X',\n'FMT\
', 'X',\n'FNE', 'X',\n'FNP', 'X',\n'FNS', 'X',\n'F\
OC', 'X',\n'FOE', 'X',\n'FOG', 'F',\n'FOH', 'X',\n\
'FOK', 'X',\n'FOL', 'X',\n'FON', 'X',\n'FOP', 'X',\
\n'FOR', 'X',\n'FOS', 'X',\n'FPA', 'X',\n'FPC', 'X\
',\n'FPI', 'X',\n'FPO', 'X',\n'FPP', 'X',\n'FPT', \
'X',\n'FQP', 'X',\n'FRA', 'X',\n'FRD', 'F',\n'FRU'\
, 'X',\n'FS3', 'X',\n'FS4', 'X',\n'FSB', 'X',\n'FS\
O', 'X',\n'FSX', 'X',\n'FTC', 'X',\n'FTP', 'X',\n'\
FTR', 'W',\n'FTT', 'X',\n'FTY', 'Y',\n'FUA', 'X',\\
n'FUC', 'X',\n'FUM', 'X',\n'FUP', 'X',\n'FVF', 'X'\
,\n'FXP', 'X',\n'FXV', 'X',\n'FYA', 'F',\n'G16', '\
X',\n'G1P', 'X',\n'G20', 'X',\n'G21', 'X',\n'G23',\
 'X',\n'G26', 'X',\n'G28', 'X',\n'G2F', 'X',\n'G37\
', 'X',\n'G39', 'X',\n'G3H', 'X',\n'G3P', 'X',\n'G\
4D', 'X',\n'G6D', 'X',\n'G6P', 'X',\n'G6Q', 'X',\n\
'G7M', 'X',\n'GA2', 'X',\n'GAA', 'X',\n'GAB', 'X',\
\n'GAC', 'X',\n'GAI', 'X',\n'GAL', 'X',\n'GAM', 'X\
',\n'GAN', 'X',\n'GAO', 'X',\n'GAP', 'X',\n'GAR', \
'G',\n'GAS', 'X',\n'GAT', 'X',\n'GBC', 'X',\n'GBI'\
, 'X',\n'GBP', 'X',\n'GBS', 'X',\n'GBX', 'X',\n'GC\
4', 'X',\n'GCA', 'X',\n'GCD', 'X',\n'GCG', 'G',\n'\
GCH', 'G',\n'GCK', 'X',\n'GCL', 'X',\n'GCM', 'X',\\
n'GCN', 'X',\n'GCO', 'X',\n'GCP', 'X',\n'GCR', 'X'\
,\n'GCS', 'X',\n'GCU', 'X',\n'GD3', 'X',\n'GDB', '\
X',\n'GDM', 'X',\n'GDN', 'X',\n'GDP', 'X',\n'GDS',\
 'X',\n'GDU', 'X',\n'GE1', 'X',\n'GE2', 'X',\n'GE3\
', 'X',\n'GEA', 'X',\n'GEL', 'X',\n'GEM', 'X',\n'G\
EN', 'X',\n'GEP', 'X',\n'GER', 'X',\n'GFP', 'X',\n\
'GGB', 'X',\n'GGL', 'E',\n'GGP', 'X',\n'GHP', 'G',\
\n'GIP', 'X',\n'GIS', 'X',\n'GKR', 'X',\n'GL2', 'X\
',\n'GL3', 'G',\n'GL4', 'X',\n'GL5', 'X',\n'GL7', \
'X',\n'GL9', 'X',\n'GLA', 'X',\n'GLB', 'X',\n'GLC'\
, 'X',\n'GLD', 'X',\n'GLE', 'X',\n'GLF', 'X',\n'GL\
G', 'X',\n'GLH', 'Q',\n'GLI', 'X',\n'GLL', 'X',\n'\
GLM', 'G',\n'GLN', 'Q',\n'GLO', 'X',\n'GLP', 'X',\\
n'GLR', 'X',\n'GLS', 'X',\n'GLT', 'X',\n'GLU', 'E'\
,\n'GLV', 'X',\n'GLW', 'X',\n'GLY', 'G',\n'GLZ', '\
X',\n'GM1', 'X',\n'GMA', 'X',\n'GMC', 'X',\n'GMH',\
 'X',\n'GMP', 'X',\n'GMY', 'X',\n'GN7', 'X',\n'GNA\
', 'X',\n'GNB', 'X',\n'GNH', 'X',\n'GNP', 'X',\n'G\
NT', 'X',\n'GOA', 'X',\n'GOL', 'X',\n'GOX', 'X',\n\
'GP1', 'X',\n'GP3', 'X',\n'GP4', 'X',\n'GP6', 'X',\
\n'GP8', 'X',\n'GPB', 'E',\n'GPC', 'X',\n'GPE', 'X\
',\n'GPG', 'X',\n'GPI', 'X',\n'GPJ', 'X',\n'GPL', \
'K',\n'GPM', 'X',\n'GPN', 'G',\n'GPP', 'X',\n'GPR'\
, 'X',\n'GPS', 'X',\n'GPX', 'X',\n'GR1', 'X',\n'GR\
3', 'X',\n'GR4', 'X',\n'GSA', 'X',\n'GSB', 'X',\n'\
GSC', 'G',\n'GSE', 'S',\n'GSH', 'X',\n'GSP', 'X',\\
n'GSR', 'X',\n'GSS', 'X',\n'GT9', 'C',\n'GTA', 'X'\
,\n'GTB', 'X',\n'GTD', 'X',\n'GTE', 'X',\n'GTH', '\
T',\n'GTN', 'X',\n'GTO', 'X',\n'GTP', 'X',\n'GTR',\
 'X',\n'GTS', 'X',\n'GTT', 'X',\n'GTX', 'X',\n'GTZ\
', 'X',\n'GU7', 'X',\n'GUA', 'X',\n'GUD', 'X',\n'G\
UM', 'X',\n'GUN', 'X',\n'GUP', 'X',\n'GUR', 'X',\n\
'GW3', 'X',\n'GZZ', 'X',\n'H2B', 'X',\n'H2P', 'H',\
\n'H2S', 'X',\n'H2U', 'X',\n'H4B', 'X',\n'H5M', 'P\
',\n'H5P', 'X',\n'HAA', 'X',\n'HAB', 'X',\n'HAC', \
'A',\n'HAD', 'X',\n'HAE', 'X',\n'HAG', 'X',\n'HAI'\
, 'X',\n'HAM', 'X',\n'HAP', 'X',\n'HAQ', 'X',\n'HA\
R', 'R',\n'HAS', 'X',\n'HAV', 'V',\n'HAX', 'X',\n'\
HAZ', 'X',\n'HBA', 'X',\n'HBC', 'X',\n'HBD', 'X',\\
n'HBI', 'X',\n'HBO', 'X',\n'HBU', 'X',\n'HBY', 'X'\
,\n'HC0', 'X',\n'HC1', 'X',\n'HC4', 'X',\n'HCA', '\
X',\n'HCC', 'X',\n'HCI', 'X',\n'HCS', 'X',\n'HDA',\
 'X',\n'HDD', 'X',\n'HDF', 'X',\n'HDN', 'X',\n'HDS\
', 'X',\n'HDZ', 'X',\n'HE1', 'X',\n'HE6', 'X',\n'H\
EA', 'X',\n'HEB', 'X',\n'HEC', 'X',\n'HED', 'X',\n\
'HEE', 'X',\n'HEF', 'X',\n'HEG', 'X',\n'HEM', 'X',\
\n'HEN', 'X',\n'HEO', 'X',\n'HEP', 'X',\n'HEU', 'X\
',\n'HEV', 'X',\n'HEX', 'X',\n'HEZ', 'X',\n'HF1', \
'X',\n'HFA', 'X',\n'HFP', 'X',\n'HGA', 'Q',\n'HGB'\
, 'X',\n'HGC', 'X',\n'HGI', 'X',\n'HGU', 'X',\n'HH\
O', 'X',\n'HHP', 'X',\n'HIB', 'X',\n'HIC', 'H',\n'\
HII', 'X',\n'HIN', 'X',\n'HIO', 'X',\n'HIP', 'H',\\
n'HIS', 'H',\n'HLE', 'X',\n'HLT', 'X',\n'HMA', 'A'\
,\n'HMB', 'X',\n'HMC', 'X',\n'HMD', 'X',\n'HMF', '\
A',\n'HMG', 'X',\n'HMH', 'X',\n'HMI', 'L',\n'HMM',\
 'X',\n'HMN', 'X',\n'HMO', 'X',\n'HMP', 'X',\n'HMR\
', 'R',\n'HNI', 'X',\n'HNP', 'X',\n'HOA', 'X',\n'H\
OE', 'X',\n'HOH', 'X',\n'HOM', 'X',\n'HOP', 'X',\n\
'HOQ', 'X',\n'HP1', 'A',\n'HP2', 'A',\n'HP3', 'X',\
\n'HPA', 'X',\n'HPB', 'X',\n'HPC', 'X',\n'HPD', 'X\
',\n'HPE', 'A',\n'HPG', 'X',\n'HPH', 'F',\n'HPP', \
'X',\n'HPQ', 'F',\n'HPR', 'X',\n'HPT', 'X',\n'HPY'\
, 'X',\n'HQO', 'X',\n'HQQ', 'X',\n'HQU', 'X',\n'HR\
G', 'R',\n'HRI', 'X',\n'HSA', 'X',\n'HSE', 'S',\n'\
HSF', 'X',\n'HSM', 'X',\n'HSO', 'H',\n'HSP', 'X',\\
n'HT1', 'X',\n'HT2', 'X',\n'HTA', 'X',\n'HTL', 'X'\
,\n'HTO', 'X',\n'HTP', 'X',\n'HTR', 'W',\n'HUP', '\
X',\n'HUX', 'X',\n'HV5', 'A',\n'HV7', 'X',\n'HV8',\
 'X',\n'HXA', 'X',\n'HXC', 'X',\n'HXP', 'X',\n'HY1\
', 'X',\n'HYA', 'X',\n'HYB', 'X',\n'HYD', 'X',\n'H\
YG', 'X',\n'HYP', 'P',\n'I06', 'X',\n'I10', 'X',\n\
'I11', 'X',\n'I17', 'X',\n'I2P', 'X',\n'I3N', 'X',\
\n'I3P', 'X',\n'I40', 'X',\n'I48', 'X',\n'I4B', 'X\
',\n'I52', 'X',\n'I5P', 'X',\n'I84', 'G',\n'IAG', \
'G',\n'IAS', 'X',\n'IB2', 'X',\n'IBB', 'X',\n'IBP'\
, 'X',\n'IBR', 'X',\n'IBS', 'X',\n'IBZ', 'X',\n'IC\
1', 'X',\n'ICA', 'X',\n'ICI', 'X',\n'ICL', 'X',\n'\
ICP', 'X',\n'ICT', 'X',\n'ICU', 'X',\n'ID2', 'X',\\
n'IDC', 'X',\n'IDG', 'X',\n'IDH', 'X',\n'IDM', 'X'\
,\n'IDO', 'X',\n'IDP', 'X',\n'IDR', 'X',\n'IDS', '\
X',\n'IDT', 'X',\n'IDU', 'X',\n'IFG', 'X',\n'IFP',\
 'X',\n'IGL', 'X',\n'IGN', 'X',\n'IGP', 'X',\n'IGU\
', 'X',\n'IH1', 'X',\n'IH2', 'X',\n'IH3', 'X',\n'I\
HB', 'X',\n'IHN', 'X',\n'IHP', 'X',\n'IIC', 'X',\n\
'IIL', 'I',\n'IIP', 'X',\n'IK2', 'X',\n'IKT', 'X',\
\n'ILA', 'I',\n'ILE', 'I',\n'ILG', 'X',\n'ILO', 'X\
',\n'ILX', 'I',\n'IM1', 'X',\n'IM2', 'X',\n'IMC', \
'X',\n'IMD', 'X',\n'IME', 'X',\n'IMF', 'X',\n'IMG'\
, 'X',\n'IMH', 'X',\n'IMI', 'X',\n'IML', 'I',\n'IM\
M', 'X',\n'IMN', 'X',\n'IMO', 'X',\n'IMP', 'X',\n'\
IMR', 'X',\n'IMU', 'X',\n'IN0', 'D',\n'IN1', 'R',\\
n'IN2', 'K',\n'IN3', 'L',\n'IN4', 'X',\n'IN5', 'A'\
,\n'IN6', 'L',\n'IN7', 'X',\n'IN8', 'X',\n'IN9', '\
X',\n'INA', 'L',\n'INB', 'X',\n'INC', 'X',\n'IND',\
 'X',\n'INE', 'X',\n'INF', 'F',\n'ING', 'F',\n'INH\
', 'R',\n'INI', 'X',\n'INJ', 'X',\n'INK', 'X',\n'I\
NL', 'X',\n'INM', 'X',\n'INN', 'A',\n'INO', 'X',\n\
'INP', 'X',\n'INQ', 'X',\n'INR', 'X',\n'INS', 'X',\
\n'INT', 'V',\n'INU', 'X',\n'INV', 'X',\n'INW', 'X\
',\n'INX', 'X',\n'INY', 'X',\n'INZ', 'X',\n'IOA', \
'X',\n'IOB', 'X',\n'IOC', 'X',\n'IOD', 'X',\n'IOE'\
, 'X',\n'IOF', 'X',\n'IOH', 'X',\n'IOL', 'X',\n'IO\
P', 'X',\n'IP1', 'X',\n'IP2', 'X',\n'IP3', 'X',\n'\
IP4', 'X',\n'IPA', 'X',\n'IPB', 'X',\n'IPD', 'X',\\
n'IPG', 'G',\n'IPH', 'X',\n'IPL', 'X',\n'IPM', 'X'\
,\n'IPN', 'X',\n'IPO', 'F',\n'IPP', 'X',\n'IPS', '\
X',\n'IPT', 'X',\n'IPU', 'X',\n'IPY', 'A',\n'IQB',\
 'X',\n'IQP', 'X',\n'IQS', 'X',\n'IR3', 'X',\n'IRI\
', 'X',\n'IRP', 'X',\n'ISA', 'X',\n'ISF', 'X',\n'I\
SO', 'X',\n'ISP', 'X',\n'ISQ', 'X',\n'ISU', 'X',\n\
'ITM', 'X',\n'ITP', 'X',\n'ITR', 'W',\n'ITS', 'X',\
\n'ITU', 'X',\n'IU5', 'X',\n'IUM', 'X',\n'IUR', 'X\
',\n'IVA', 'X',\n'IYG', 'G',\n'IYR', 'Y',\n'J77', \
'X',\n'J78', 'X',\n'J80', 'X',\n'JE2', 'X',\n'JEN'\
, 'X',\n'JST', 'X',\n'K21', 'X',\n'KAH', 'X',\n'KA\
I', 'X',\n'KAM', 'X',\n'KAN', 'X',\n'KAP', 'X',\n'\
KCP', 'X',\n'KCX', 'K',\n'KDO', 'X',\n'KEF', 'X',\\
n'KET', 'X',\n'KGR', 'X',\n'KH1', 'X',\n'KIF', 'X'\
,\n'KIV', 'V',\n'KNI', 'X',\n'KPH', 'K',\n'KTH', '\
X',\n'KTN', 'X',\n'KTP', 'X',\n'KWT', 'X',\n'L04',\
 'X',\n'L1P', 'X',\n'L24', 'E',\n'L2P', 'X',\n'L34\
', 'E',\n'L37', 'E',\n'L3P', 'X',\n'L4P', 'X',\n'L\
75', 'X',\n'LAC', 'X',\n'LAD', 'X',\n'LAK', 'X',\n\
'LAM', 'X',\n'LAR', 'X',\n'LAT', 'X',\n'LAX', 'X',\
\n'LCO', 'X',\n'LCP', 'X',\n'LCS', 'X',\n'LDA', 'X\
',\n'LDO', 'L',\n'LDP', 'X',\n'LEA', 'X',\n'LEO', \
'X',\n'LEU', 'L',\n'LG2', 'X',\n'LG6', 'X',\n'LGC'\
, 'X',\n'LGP', 'X',\n'LHG', 'X',\n'LHY', 'F',\n'LI\
1', 'X',\n'LIG', 'X',\n'LIL', 'X',\n'LIM', 'X',\n'\
LIN', 'X',\n'LIO', 'X',\n'LIP', 'X',\n'LLA', 'X',\\
n'LLP', 'K',\n'LLY', 'K',\n'LMG', 'X',\n'LML', 'X'\
,\n'LMT', 'X',\n'LMU', 'X',\n'LMZ', 'X',\n'LNK', '\
X',\n'LNL', 'X',\n'LNO', 'X',\n'LOF', 'X',\n'LOL',\
 'L',\n'LOM', 'X',\n'LOR', 'X',\n'LOS', 'X',\n'LOV\
', 'L',\n'LOX', 'X',\n'LP1', 'X',\n'LP2', 'R',\n'L\
PA', 'X',\n'LPC', 'X',\n'LPF', 'X',\n'LPL', 'X',\n\
'LPM', 'X',\n'LPP', 'X',\n'LRB', 'X',\n'LRU', 'X',\
\n'LS1', 'X',\n'LS2', 'X',\n'LS3', 'X',\n'LS4', 'X\
',\n'LS5', 'X',\n'LTA', 'X',\n'LTL', 'X',\n'LTR', \
'W',\n'LUM', 'X',\n'LVS', 'L',\n'LXC', 'X',\n'LY2'\
, 'X',\n'LY3', 'X',\n'LYA', 'X',\n'LYB', 'X',\n'LY\
C', 'X',\n'LYD', 'X',\n'LYM', 'K',\n'LYN', 'X',\n'\
LYS', 'K',\n'LYT', 'X',\n'LYW', 'X',\n'LYZ', 'K',\\
n'M1A', 'X',\n'M1G', 'X',\n'M2G', 'X',\n'M3L', 'K'\
,\n'M6P', 'X',\n'M6T', 'X',\n'M7G', 'X',\n'MA1', '\
X',\n'MA2', 'X',\n'MA3', 'X',\n'MA4', 'X',\n'MA6',\
 'X',\n'MAA', 'A',\n'MAB', 'X',\n'MAC', 'X',\n'MAE\
', 'X',\n'MAG', 'X',\n'MAH', 'X',\n'MAI', 'R',\n'M\
AK', 'X',\n'MAL', 'X',\n'MAM', 'X',\n'MAN', 'X',\n\
'MAO', 'X',\n'MAP', 'X',\n'MAR', 'X',\n'MAS', 'X',\
\n'MAT', 'X',\n'MAU', 'X',\n'MAZ', 'X',\n'MBA', 'X\
',\n'MBD', 'X',\n'MBG', 'X',\n'MBH', 'X',\n'MBN', \
'X',\n'MBO', 'X',\n'MBR', 'X',\n'MBS', 'X',\n'MBV'\
, 'X',\n'MBZ', 'X',\n'MCA', 'X',\n'MCD', 'X',\n'MC\
E', 'X',\n'MCG', 'G',\n'MCI', 'X',\n'MCN', 'X',\n'\
MCP', 'X',\n'MCT', 'X',\n'MCY', 'X',\n'MD2', 'X',\\
n'MDA', 'X',\n'MDC', 'X',\n'MDG', 'X',\n'MDH', 'X'\
,\n'MDL', 'X',\n'MDM', 'X',\n'MDN', 'X',\n'MDP', '\
X',\n'ME6', 'X',\n'MEB', 'X',\n'MEC', 'X',\n'MEL',\
 'X',\n'MEN', 'N',\n'MEP', 'X',\n'MER', 'X',\n'MES\
', 'X',\n'MET', 'M',\n'MEV', 'X',\n'MF2', 'X',\n'M\
F3', 'M',\n'MFB', 'X',\n'MFD', 'X',\n'MFU', 'X',\n\
'MG7', 'X',\n'MGA', 'X',\n'MGB', 'X',\n'MGD', 'X',\
\n'MGG', 'R',\n'MGL', 'X',\n'MGN', 'Q',\n'MGO', 'X\
',\n'MGP', 'X',\n'MGR', 'X',\n'MGS', 'X',\n'MGT', \
'X',\n'MGU', 'X',\n'MGY', 'G',\n'MHB', 'X',\n'MHF'\
, 'X',\n'MHL', 'L',\n'MHM', 'X',\n'MHO', 'M',\n'MH\
S', 'H',\n'MHZ', 'X',\n'MIA', 'X',\n'MIC', 'X',\n'\
MID', 'X',\n'MIL', 'X',\n'MIM', 'X',\n'MIN', 'G',\\
n'MIP', 'X',\n'MIS', 'S',\n'MIT', 'X',\n'MJI', 'X'\
,\n'MK1', 'X',\n'MKC', 'X',\n'MLA', 'X',\n'MLC', '\
X',\n'MLE', 'L',\n'MLN', 'X',\n'MLT', 'X',\n'MLY',\
 'K',\n'MLZ', 'K',\n'MM3', 'X',\n'MM4', 'X',\n'MMA\
', 'X',\n'MMC', 'X',\n'MME', 'M',\n'MMO', 'R',\n'M\
MP', 'X',\n'MMQ', 'X',\n'MMT', 'X',\n'MN1', 'X',\n\
'MN2', 'X',\n'MN3', 'X',\n'MN5', 'X',\n'MN7', 'X',\
\n'MN8', 'X',\n'MNA', 'X',\n'MNB', 'X',\n'MNC', 'X\
',\n'MNG', 'X',\n'MNL', 'L',\n'MNO', 'X',\n'MNP', \
'X',\n'MNQ', 'X',\n'MNS', 'X',\n'MNT', 'X',\n'MNV'\
, 'V',\n'MO1', 'X',\n'MO2', 'X',\n'MO3', 'X',\n'MO\
4', 'X',\n'MO5', 'X',\n'MO6', 'X',\n'MOA', 'X',\n'\
MOB', 'X',\n'MOC', 'X',\n'MOE', 'X',\n'MOG', 'X',\\
n'MOH', 'X',\n'MOL', 'X',\n'MOO', 'X',\n'MOP', 'X'\
,\n'MOR', 'X',\n'MOS', 'X',\n'MOT', 'X',\n'MOX', '\
X',\n'MP1', 'X',\n'MP3', 'X',\n'MPA', 'X',\n'MPB',\
 'X',\n'MPC', 'X',\n'MPD', 'X',\n'MPG', 'X',\n'MPH\
', 'M',\n'MPI', 'X',\n'MPJ', 'M',\n'MPL', 'X',\n'M\
PN', 'X',\n'MPO', 'X',\n'MPP', 'X',\n'MPQ', 'G',\n\
'MPR', 'X',\n'MPS', 'X',\n'MQ0', 'X',\n'MQ7', 'X',\
\n'MQ8', 'X',\n'MQ9', 'X',\n'MQI', 'X',\n'MR2', 'X\
',\n'MRC', 'X',\n'MRM', 'X',\n'MRP', 'X',\n'MS2', \
'X',\n'MSA', 'X',\n'MSB', 'X',\n'MSD', 'X',\n'MSE'\
, 'M',\n'MSF', 'X',\n'MSI', 'X',\n'MSO', 'M',\n'MS\
Q', 'X',\n'MST', 'X',\n'MSU', 'X',\n'MTA', 'X',\n'\
MTB', 'X',\n'MTC', 'X',\n'MTD', 'X',\n'MTE', 'X',\\
n'MTF', 'X',\n'MTG', 'X',\n'MTO', 'X',\n'MTS', 'X'\
,\n'MTT', 'X',\n'MTX', 'X',\n'MTY', 'Y',\n'MUG', '\
X',\n'MUP', 'X',\n'MUR', 'X',\n'MVA', 'V',\n'MW1',\
 'X',\n'MW2', 'X',\n'MXA', 'X',\n'MXY', 'X',\n'MYA\
', 'X',\n'MYC', 'X',\n'MYG', 'X',\n'MYR', 'X',\n'M\
YS', 'X',\n'MYT', 'X',\n'MZM', 'X',\n'N1T', 'X',\n\
'N25', 'X',\n'N2B', 'X',\n'N3T', 'X',\n'N4B', 'X',\
\n'NA2', 'X',\n'NA5', 'X',\n'NA6', 'X',\n'NAA', 'X\
',\n'NAB', 'X',\n'NAC', 'X',\n'NAD', 'X',\n'NAE', \
'X',\n'NAF', 'X',\n'NAG', 'X',\n'NAH', 'X',\n'NAI'\
, 'X',\n'NAL', 'A',\n'NAM', 'A',\n'NAN', 'X',\n'NA\
O', 'X',\n'NAP', 'X',\n'NAQ', 'X',\n'NAR', 'X',\n'\
NAS', 'X',\n'NAU', 'X',\n'NAV', 'X',\n'NAW', 'X',\\
n'NAX', 'X',\n'NAY', 'X',\n'NBA', 'X',\n'NBD', 'X'\
,\n'NBE', 'X',\n'NBG', 'X',\n'NBN', 'X',\n'NBP', '\
X',\n'NBS', 'X',\n'NBU', 'X',\n'NCA', 'X',\n'NCB',\
 'A',\n'NCD', 'X',\n'NCH', 'X',\n'NCM', 'X',\n'NCN\
', 'X',\n'NCO', 'X',\n'NCR', 'X',\n'NCS', 'X',\n'N\
D4', 'X',\n'NDA', 'X',\n'NDC', 'X',\n'NDD', 'X',\n\
'NDO', 'X',\n'NDP', 'X',\n'NDT', 'X',\n'NEA', 'X',\
\n'NEB', 'X',\n'NED', 'X',\n'NEM', 'H',\n'NEN', 'X\
',\n'NEO', 'X',\n'NEP', 'H',\n'NEQ', 'X',\n'NES', \
'X',\n'NET', 'X',\n'NEV', 'X',\n'NFA', 'F',\n'NFE'\
, 'X',\n'NFG', 'X',\n'NFP', 'X',\n'NFS', 'X',\n'NG\
6', 'X',\n'NGA', 'X',\n'NGL', 'X',\n'NGM', 'X',\n'\
NGO', 'X',\n'NGP', 'X',\n'NGT', 'X',\n'NGU', 'X',\\
n'NH2', 'X',\n'NH3', 'X',\n'NH4', 'X',\n'NHD', 'X'\
,\n'NHE', 'X',\n'NHM', 'X',\n'NHP', 'X',\n'NHR', '\
X',\n'NHS', 'X',\n'NI1', 'X',\n'NI2', 'X',\n'NIC',\
 'X',\n'NID', 'X',\n'NIK', 'X',\n'NIO', 'X',\n'NIP\
', 'X',\n'NIT', 'X',\n'NIU', 'X',\n'NIY', 'Y',\n'N\
LA', 'X',\n'NLE', 'L',\n'NLG', 'X',\n'NLN', 'L',\n\
'NLP', 'L',\n'NM1', 'X',\n'NMA', 'A',\n'NMB', 'X',\
\n'NMC', 'G',\n'NMD', 'X',\n'NME', 'X',\n'NMN', 'X\
',\n'NMO', 'X',\n'NMQ', 'X',\n'NMX', 'X',\n'NMY', \
'X',\n'NNH', 'R',\n'NNO', 'X',\n'NO2', 'X',\n'NO3'\
, 'X',\n'NOA', 'X',\n'NOD', 'X',\n'NOJ', 'X',\n'NO\
N', 'X',\n'NOP', 'X',\n'NOR', 'X',\n'NOS', 'X',\n'\
NOV', 'X',\n'NOX', 'X',\n'NP3', 'X',\n'NPA', 'X',\\
n'NPC', 'X',\n'NPD', 'X',\n'NPE', 'X',\n'NPF', 'X'\
,\n'NPH', 'C',\n'NPI', 'X',\n'NPL', 'X',\n'NPN', '\
X',\n'NPO', 'X',\n'NPP', 'X',\n'NPT', 'X',\n'NPY',\
 'X',\n'NRG', 'R',\n'NRI', 'X',\n'NS1', 'X',\n'NS5\
', 'X',\n'NSP', 'X',\n'NTA', 'X',\n'NTB', 'X',\n'N\
TC', 'X',\n'NTH', 'X',\n'NTM', 'X',\n'NTP', 'X',\n\
'NTS', 'X',\n'NTU', 'X',\n'NTZ', 'X',\n'NU1', 'X',\
\n'NVA', 'V',\n'NVI', 'X',\n'NVP', 'X',\n'NW1', 'X\
',\n'NYP', 'X',\n'O4M', 'X',\n'OAA', 'X',\n'OAI', \
'X',\n'OAP', 'X',\n'OAR', 'X',\n'OAS', 'S',\n'OBA'\
, 'X',\n'OBN', 'X',\n'OC1', 'X',\n'OC2', 'X',\n'OC\
3', 'X',\n'OC4', 'X',\n'OC5', 'X',\n'OC6', 'X',\n'\
OC7', 'X',\n'OCL', 'X',\n'OCM', 'X',\n'OCN', 'X',\\
n'OCO', 'X',\n'OCP', 'X',\n'OCS', 'C',\n'OCT', 'X'\
,\n'OCV', 'K',\n'OCY', 'C',\n'ODA', 'X',\n'ODS', '\
X',\n'OES', 'X',\n'OET', 'X',\n'OF1', 'X',\n'OF2',\
 'X',\n'OF3', 'X',\n'OFL', 'X',\n'OFO', 'X',\n'OHE\
', 'X',\n'OHO', 'X',\n'OHT', 'X',\n'OIC', 'X',\n'O\
IP', 'X',\n'OKA', 'X',\n'OLA', 'X',\n'OLE', 'X',\n\
'OLI', 'X',\n'OLO', 'X',\n'OMB', 'X',\n'OMC', 'X',\
\n'OMD', 'X',\n'OME', 'X',\n'OMG', 'X',\n'OMP', 'X\
',\n'OMT', 'M',\n'OMU', 'X',\n'ONE', 'X',\n'ONL', \
'L',\n'ONP', 'X',\n'OPA', 'X',\n'OPD', 'X',\n'OPE'\
, 'X',\n'OPG', 'X',\n'OPH', 'X',\n'OPN', 'X',\n'OP\
P', 'X',\n'OPR', 'R',\n'ORN', 'X',\n'ORO', 'X',\n'\
ORP', 'X',\n'OSB', 'X',\n'OSS', 'X',\n'OTA', 'X',\\
n'OTB', 'X',\n'OTE', 'X',\n'OTG', 'X',\n'OUT', 'X'\
,\n'OVA', 'X',\n'OWQ', 'X',\n'OXA', 'X',\n'OXE', '\
X',\n'OXI', 'X',\n'OXL', 'X',\n'OXM', 'X',\n'OXN',\
 'X',\n'OXO', 'X',\n'OXP', 'X',\n'OXS', 'X',\n'OXY\
', 'X',\n'P11', 'A',\n'P24', 'X',\n'P28', 'X',\n'P\
2P', 'X',\n'P2U', 'X',\n'P3M', 'X',\n'P4C', 'X',\n\
'P4P', 'X',\n'P5P', 'X',\n'P6G', 'X',\n'PA1', 'X',\
\n'PA2', 'X',\n'PA3', 'X',\n'PA4', 'X',\n'PA5', 'X\
',\n'PAA', 'X',\n'PAB', 'X',\n'PAC', 'X',\n'PAD', \
'X',\n'PAE', 'X',\n'PAG', 'X',\n'PAH', 'X',\n'PAI'\
, 'X',\n'PAL', 'D',\n'PAM', 'X',\n'PAN', 'X',\n'PA\
O', 'X',\n'PAP', 'A',\n'PAQ', 'F',\n'PAR', 'X',\n'\
PAS', 'X',\n'PAT', 'W',\n'PBA', 'X',\n'PBB', 'X',\\
n'PBC', 'X',\n'PBF', 'F',\n'PBG', 'X',\n'PBI', 'X'\
,\n'PBM', 'X',\n'PBN', 'X',\n'PBP', 'X',\n'PBR', '\
X',\n'PBZ', 'X',\n'PC2', 'X',\n'PCA', 'E',\n'PCB',\
 'X',\n'PCD', 'X',\n'PCE', 'X',\n'PCG', 'X',\n'PCH\
', 'X',\n'PCL', 'X',\n'PCM', 'X',\n'PCP', 'X',\n'P\
CR', 'X',\n'PCS', 'X',\n'PCU', 'X',\n'PCV', 'X',\n\
'PCY', 'X',\n'PD1', 'X',\n'PDA', 'X',\n'PDC', 'X',\
\n'PDD', 'A',\n'PDE', 'A',\n'PDI', 'X',\n'PDL', 'A\
',\n'PDN', 'X',\n'PDO', 'X',\n'PDP', 'X',\n'PDT', \
'X',\n'PDU', 'X',\n'PE2', 'X',\n'PE6', 'X',\n'PEA'\
, 'X',\n'PEB', 'X',\n'PEC', 'X',\n'PED', 'X',\n'PE\
E', 'X',\n'PEF', 'X',\n'PEG', 'X',\n'PEL', 'X',\n'\
PEO', 'X',\n'PEP', 'X',\n'PEQ', 'X',\n'PER', 'X',\\
n'PET', 'X',\n'PFB', 'X',\n'PFC', 'X',\n'PFG', 'X'\
,\n'PFL', 'X',\n'PFM', 'X',\n'PFZ', 'X',\n'PG4', '\
X',\n'PG5', 'X',\n'PG6', 'X',\n'PGA', 'X',\n'PGC',\
 'X',\n'PGD', 'X',\n'PGE', 'X',\n'PGG', 'G',\n'PGH\
', 'X',\n'PGL', 'X',\n'PGO', 'X',\n'PGP', 'X',\n'P\
GQ', 'X',\n'PGR', 'X',\n'PGS', 'X',\n'PGU', 'X',\n\
'PGX', 'X',\n'PGY', 'G',\n'PH1', 'X',\n'PH2', 'X',\
\n'PH3', 'X',\n'PHA', 'F',\n'PHB', 'X',\n'PHC', 'X\
',\n'PHD', 'X',\n'PHE', 'F',\n'PHG', 'X',\n'PHH', \
'X',\n'PHI', 'F',\n'PHL', 'F',\n'PHM', 'X',\n'PHN'\
, 'X',\n'PHO', 'X',\n'PHP', 'X',\n'PHQ', 'X',\n'PH\
S', 'H',\n'PHT', 'X',\n'PHW', 'P',\n'PHY', 'X',\n'\
PI1', 'X',\n'PI2', 'X',\n'PI3', 'X',\n'PI4', 'X',\\
n'PI5', 'X',\n'PI6', 'X',\n'PI7', 'X',\n'PI8', 'X'\
,\n'PI9', 'X',\n'PIA', 'X',\n'PIB', 'X',\n'PIC', '\
X',\n'PID', 'X',\n'PIG', 'X',\n'PIH', 'X',\n'PIM',\
 'X',\n'PIN', 'X',\n'PIO', 'X',\n'PIP', 'X',\n'PIQ\
', 'X',\n'PIR', 'X',\n'PIV', 'X',\n'PKF', 'X',\n'P\
L1', 'X',\n'PL9', 'X',\n'PLA', 'D',\n'PLC', 'X',\n\
'PLE', 'L',\n'PLG', 'G',\n'PLH', 'X',\n'PLM', 'X',\
\n'PLP', 'X',\n'PLS', 'S',\n'PLT', 'W',\n'PLU', 'L\
',\n'PLY', 'X',\n'PMA', 'X',\n'PMB', 'X',\n'PMC', \
'X',\n'PME', 'F',\n'PML', 'X',\n'PMM', 'X',\n'PMO'\
, 'X',\n'PMP', 'X',\n'PMS', 'X',\n'PMY', 'X',\n'PN\
2', 'X',\n'PNA', 'X',\n'PNB', 'X',\n'PNC', 'G',\n'\
PND', 'X',\n'PNE', 'A',\n'PNF', 'X',\n'PNG', 'X',\\
n'PNI', 'X',\n'PNL', 'X',\n'PNM', 'X',\n'PNN', 'X'\
,\n'PNO', 'X',\n'PNP', 'X',\n'PNQ', 'X',\n'PNS', '\
X',\n'PNT', 'X',\n'PNU', 'X',\n'PO2', 'X',\n'PO4',\
 'X',\n'POB', 'X',\n'POC', 'X',\n'POL', 'X',\n'POM\
', 'P',\n'PON', 'X',\n'POP', 'X',\n'POR', 'X',\n'P\
OS', 'X',\n'PP1', 'X',\n'PP2', 'X',\n'PP3', 'A',\n\
'PP4', 'X',\n'PP5', 'X',\n'PP6', 'X',\n'PP7', 'X',\
\n'PP8', 'N',\n'PP9', 'X',\n'PPB', 'X',\n'PPC', 'X\
',\n'PPD', 'X',\n'PPE', 'E',\n'PPG', 'X',\n'PPH', \
'F',\n'PPI', 'X',\n'PPJ', 'V',\n'PPL', 'X',\n'PPM'\
, 'X',\n'PPN', 'A',\n'PPO', 'X',\n'PPP', 'X',\n'PP\
Q', 'X',\n'PPR', 'X',\n'PPS', 'X',\n'PPT', 'X',\n'\
PPU', 'X',\n'PPX', 'F',\n'PPY', 'X',\n'PPZ', 'X',\\
n'PQ0', 'X',\n'PQN', 'X',\n'PQQ', 'X',\n'PR1', 'X'\
,\n'PR2', 'X',\n'PR3', 'X',\n'PRA', 'X',\n'PRB', '\
X',\n'PRC', 'X',\n'PRD', 'X',\n'PRE', 'X',\n'PRF',\
 'X',\n'PRH', 'X',\n'PRI', 'P',\n'PRL', 'X',\n'PRN\
', 'X',\n'PRO', 'P',\n'PRP', 'X',\n'PRR', 'A',\n'P\
RS', 'P',\n'PRZ', 'X',\n'PS0', 'X',\n'PSA', 'X',\n\
'PSD', 'X',\n'PSE', 'X',\n'PSF', 'S',\n'PSG', 'X',\
\n'PSI', 'X',\n'PSO', 'X',\n'PSQ', 'X',\n'PSS', 'X\
',\n'PST', 'X',\n'PSU', 'X',\n'PT1', 'X',\n'PT3', \
'X',\n'PTA', 'X',\n'PTC', 'X',\n'PTD', 'X',\n'PTE'\
, 'X',\n'PTH', 'Y',\n'PTL', 'X',\n'PTM', 'Y',\n'PT\
N', 'X',\n'PTO', 'X',\n'PTP', 'X',\n'PTR', 'Y',\n'\
PTS', 'X',\n'PTT', 'X',\n'PTU', 'X',\n'PTY', 'X',\\
n'PUA', 'X',\n'PUB', 'X',\n'PUR', 'X',\n'PUT', 'X'\
,\n'PVA', 'X',\n'PVB', 'X',\n'PVH', 'H',\n'PVL', '\
X',\n'PXA', 'X',\n'PXF', 'X',\n'PXG', 'X',\n'PXP',\
 'X',\n'PXY', 'X',\n'PXZ', 'X',\n'PY2', 'X',\n'PY4\
', 'X',\n'PY5', 'X',\n'PY6', 'X',\n'PYA', 'A',\n'P\
YC', 'X',\n'PYD', 'X',\n'PYE', 'X',\n'PYL', 'X',\n\
'PYM', 'X',\n'PYO', 'X',\n'PYP', 'X',\n'PYQ', 'X',\
\n'PYR', 'X',\n'PYS', 'X',\n'PYT', 'X',\n'PYX', 'X\
',\n'PYY', 'X',\n'PYZ', 'X',\n'PZQ', 'X',\n'Q82', \
'X',\n'QNC', 'X',\n'QND', 'X',\n'QSI', 'Q',\n'QTR'\
, 'X',\n'QUA', 'X',\n'QUE', 'X',\n'QUI', 'X',\n'QU\
O', 'X',\n'R11', 'X',\n'R12', 'X',\n'R13', 'X',\n'\
R18', 'X',\n'R1P', 'X',\n'R56', 'X',\n'R5P', 'X',\\
n'RA2', 'X',\n'RAD', 'X',\n'RAI', 'X',\n'RAL', 'X'\
,\n'RAM', 'X',\n'RAN', 'X',\n'RAP', 'X',\n'RBF', '\
X',\n'RBU', 'X',\n'RCA', 'X',\n'RCL', 'X',\n'RCO',\
 'X',\n'RDC', 'X',\n'RDF', 'W',\n'RE9', 'X',\n'REA\
', 'X',\n'RED', 'K',\n'REO', 'X',\n'REP', 'X',\n'R\
ET', 'X',\n'RFA', 'X',\n'RFB', 'X',\n'RFL', 'X',\n\
'RFP', 'X',\n'RG1', 'X',\n'RGS', 'X',\n'RH1', 'X',\
\n'RHA', 'X',\n'RHC', 'X',\n'RHD', 'X',\n'RHM', 'X\
',\n'RHO', 'X',\n'RHQ', 'X',\n'RHS', 'X',\n'RIA', \
'X',\n'RIB', 'X',\n'RIC', 'X',\n'RIF', 'X',\n'RIN'\
, 'X',\n'RIP', 'X',\n'RIT', 'X',\n'RMB', 'X',\n'RM\
N', 'X',\n'RMP', 'X',\n'RNG', 'X',\n'RNS', 'X',\n'\
RNT', 'X',\n'RO2', 'X',\n'RO4', 'X',\n'ROC', 'N',\\
n'ROI', 'X',\n'ROM', 'X',\n'RON', 'V',\n'ROP', 'X'\
,\n'ROS', 'X',\n'ROX', 'X',\n'RPA', 'X',\n'RPD', '\
X',\n'RPH', 'X',\n'RPL', 'X',\n'RPP', 'X',\n'RPR',\
 'X',\n'RPX', 'X',\n'RQ3', 'X',\n'RR1', 'X',\n'RR6\
', 'X',\n'RRS', 'X',\n'RS1', 'X',\n'RS2', 'X',\n'R\
S7', 'X',\n'RSS', 'X',\n'RTA', 'X',\n'RTB', 'X',\n\
'RTC', 'X',\n'RTL', 'X',\n'RUB', 'X',\n'RUN', 'X',\
\n'RWJ', 'X',\n'RXP', 'X',\n'S02', 'X',\n'S11', 'X\
',\n'S1H', 'S',\n'S27', 'X',\n'S2C', 'C',\n'S3P', \
'X',\n'S4U', 'X',\n'S57', 'X',\n'S58', 'X',\n'S5H'\
, 'X',\n'S6G', 'X',\n'S80', 'X',\n'SAA', 'X',\n'SA\
B', 'X',\n'SAC', 'S',\n'SAD', 'X',\n'SAE', 'X',\n'\
SAF', 'X',\n'SAH', 'C',\n'SAI', 'C',\n'SAL', 'X',\\
n'SAM', 'M',\n'SAN', 'X',\n'SAP', 'X',\n'SAR', 'X'\
,\n'SAS', 'X',\n'SB1', 'X',\n'SB2', 'X',\n'SB3', '\
X',\n'SB4', 'X',\n'SB5', 'X',\n'SB6', 'X',\n'SBA',\
 'L',\n'SBB', 'X',\n'SBD', 'A',\n'SBI', 'X',\n'SBL\
', 'A',\n'SBN', 'X',\n'SBO', 'X',\n'SBR', 'X',\n'S\
BS', 'X',\n'SBT', 'X',\n'SBU', 'X',\n'SBX', 'X',\n\
'SC4', 'X',\n'SCA', 'X',\n'SCC', 'X',\n'SCD', 'X',\
\n'SCH', 'C',\n'SCI', 'X',\n'SCL', 'X',\n'SCM', 'X\
',\n'SCN', 'X',\n'SCO', 'X',\n'SCP', 'S',\n'SCR', \
'X',\n'SCS', 'X',\n'SCV', 'C',\n'SCY', 'C',\n'SD8'\
, 'X',\n'SDK', 'X',\n'SDZ', 'X',\n'SE4', 'X',\n'SE\
A', 'X',\n'SEB', 'S',\n'SEC', 'X',\n'SEG', 'A',\n'\
SEI', 'X',\n'SEL', 'S',\n'SEM', 'X',\n'SEO', 'X',\\
n'SEP', 'S',\n'SER', 'S',\n'SES', 'X',\n'SET', 'S'\
,\n'SEU', 'X',\n'SF4', 'X',\n'SFG', 'X',\n'SFN', '\
X',\n'SFO', 'X',\n'SGA', 'X',\n'SGC', 'X',\n'SGL',\
 'X',\n'SGM', 'X',\n'SGN', 'X',\n'SGP', 'X',\n'SHA\
', 'X',\n'SHC', 'X',\n'SHF', 'X',\n'SHH', 'X',\n'S\
HP', 'G',\n'SHR', 'E',\n'SHT', 'T',\n'SHU', 'X',\n\
'SI2', 'X',\n'SIA', 'X',\n'SIF', 'X',\n'SIG', 'X',\
\n'SIH', 'X',\n'SIM', 'X',\n'SIN', 'X',\n'SKD', 'X\
',\n'SKF', 'X',\n'SLB', 'X',\n'SLE', 'X',\n'SLZ', \
'K',\n'SMA', 'X',\n'SMC', 'C',\n'SME', 'M',\n'SML'\
, 'X',\n'SMM', 'M',\n'SMN', 'X',\n'SMP', 'X',\n'SM\
S', 'X',\n'SN1', 'X',\n'SN6', 'X',\n'SN7', 'X',\n'\
SNC', 'C',\n'SNN', 'X',\n'SNP', 'X',\n'SO1', 'X',\\
n'SO2', 'X',\n'SO3', 'X',\n'SO4', 'X',\n'SOA', 'X'\
,\n'SOC', 'C',\n'SOM', 'X',\n'SOR', 'X',\n'SOT', '\
X',\n'SOX', 'X',\n'SPA', 'X',\n'SPB', 'X',\n'SPC',\
 'X',\n'SPD', 'X',\n'SPE', 'X',\n'SPG', 'X',\n'SPH\
', 'X',\n'SPI', 'X',\n'SPK', 'X',\n'SPM', 'X',\n'S\
PN', 'X',\n'SPO', 'X',\n'SPP', 'X',\n'SPS', 'X',\n\
'SPY', 'X',\n'SQU', 'X',\n'SRA', 'X',\n'SRB', 'X',\
\n'SRD', 'X',\n'SRL', 'X',\n'SRM', 'X',\n'SRS', 'X\
',\n'SRY', 'X',\n'SSA', 'X',\n'SSB', 'X',\n'SSG', \
'X',\n'SSP', 'X',\n'ST1', 'X',\n'ST2', 'X',\n'ST3'\
, 'X',\n'ST4', 'X',\n'ST5', 'X',\n'ST6', 'X',\n'ST\
A', 'X',\n'STB', 'X',\n'STE', 'X',\n'STG', 'X',\n'\
STI', 'X',\n'STL', 'X',\n'STN', 'X',\n'STO', 'X',\\
n'STP', 'X',\n'STR', 'X',\n'STU', 'X',\n'STY', 'Y'\
,\n'SU1', 'X',\n'SU2', 'X',\n'SUC', 'X',\n'SUI', '\
X',\n'SUL', 'X',\n'SUR', 'X',\n'SVA', 'S',\n'SWA',\
 'X',\n'T16', 'X',\n'T19', 'X',\n'T23', 'X',\n'T29\
', 'X',\n'T33', 'X',\n'T3P', 'X',\n'T42', 'A',\n'T\
44', 'X',\n'T5A', 'X',\n'T6A', 'T',\n'T6P', 'X',\n\
'T80', 'X',\n'T87', 'X',\n'TA1', 'X',\n'TAA', 'X',\
\n'TAB', 'X',\n'TAC', 'X',\n'TAD', 'X',\n'TAF', 'X\
',\n'TAM', 'X',\n'TAP', 'X',\n'TAR', 'X',\n'TAS', \
'X',\n'TAU', 'X',\n'TAX', 'X',\n'TAZ', 'X',\n'TB9'\
, 'X',\n'TBA', 'X',\n'TBD', 'X',\n'TBG', 'G',\n'TB\
H', 'X',\n'TBM', 'T',\n'TBO', 'X',\n'TBP', 'X',\n'\
TBR', 'X',\n'TBS', 'X',\n'TBT', 'X',\n'TBU', 'X',\\
n'TBZ', 'X',\n'TC4', 'X',\n'TCA', 'X',\n'TCB', 'X'\
,\n'TCH', 'X',\n'TCK', 'X',\n'TCL', 'X',\n'TCM', '\
X',\n'TCN', 'X',\n'TCP', 'X',\n'TCR', 'W',\n'TCS',\
 'X',\n'TCZ', 'X',\n'TDA', 'X',\n'TDB', 'X',\n'TDG\
', 'X',\n'TDP', 'X',\n'TDR', 'X',\n'TDX', 'X',\n'T\
EA', 'X',\n'TEM', 'X',\n'TEN', 'X',\n'TEO', 'X',\n\
'TEP', 'X',\n'TER', 'X',\n'TES', 'X',\n'TET', 'X',\
\n'TFA', 'X',\n'TFB', 'X',\n'TFH', 'X',\n'TFI', 'X\
',\n'TFK', 'X',\n'TFP', 'X',\n'THA', 'X',\n'THB', \
'X',\n'THC', 'T',\n'THD', 'X',\n'THE', 'X',\n'THF'\
, 'X',\n'THJ', 'X',\n'THK', 'X',\n'THM', 'X',\n'TH\
N', 'X',\n'THO', 'T',\n'THP', 'X',\n'THQ', 'X',\n'\
THR', 'T',\n'THS', 'X',\n'THT', 'X',\n'THU', 'X',\\
n'THX', 'X',\n'THZ', 'X',\n'TI1', 'X',\n'TI2', 'X'\
,\n'TI3', 'P',\n'TIA', 'X',\n'TIH', 'A',\n'TK4', '\
X',\n'TLA', 'X',\n'TLC', 'X',\n'TLM', 'X',\n'TLN',\
 'X',\n'TLX', 'X',\n'TM5', 'X',\n'TM6', 'X',\n'TMA\
', 'X',\n'TMB', 'T',\n'TMC', 'X',\n'TMD', 'T',\n'T\
ME', 'X',\n'TMF', 'X',\n'TML', 'K',\n'TMM', 'X',\n\
'TMN', 'X',\n'TMP', 'X',\n'TMQ', 'X',\n'TMR', 'X',\
\n'TMT', 'X',\n'TMZ', 'X',\n'TNB', 'C',\n'TND', 'X\
',\n'TNK', 'X',\n'TNP', 'X',\n'TNT', 'X',\n'TOA', \
'X',\n'TOB', 'X',\n'TOC', 'X',\n'TOL', 'X',\n'TOP'\
, 'X',\n'TOS', 'X',\n'TOT', 'X',\n'TP1', 'G',\n'TP\
2', 'P',\n'TP3', 'E',\n'TP4', 'E',\n'TP7', 'T',\n'\
TPA', 'X',\n'TPE', 'X',\n'TPF', 'X',\n'TPI', 'X',\\
n'TPL', 'W',\n'TPM', 'X',\n'TPN', 'G',\n'TPO', 'T'\
,\n'TPP', 'X',\n'TPQ', 'A',\n'TPR', 'P',\n'TPS', '\
X',\n'TPT', 'X',\n'TPV', 'X',\n'TPX', 'X',\n'TPY',\
 'X',\n'TQ3', 'X',\n'TQ4', 'X',\n'TQ5', 'X',\n'TQ6\
', 'X',\n'TR1', 'X',\n'TRA', 'X',\n'TRB', 'X',\n'T\
RC', 'X',\n'TRD', 'X',\n'TRE', 'X',\n'TRF', 'W',\n\
'TRG', 'K',\n'TRH', 'X',\n'TRI', 'X',\n'TRJ', 'X',\
\n'TRM', 'X',\n'TRN', 'W',\n'TRO', 'W',\n'TRP', 'W\
',\n'TRQ', 'X',\n'TRS', 'X',\n'TRX', 'W',\n'TRZ', \
'X',\n'TS2', 'X',\n'TS3', 'X',\n'TS4', 'X',\n'TS5'\
, 'X',\n'TSA', 'X',\n'TSB', 'X',\n'TSI', 'X',\n'TS\
M', 'X',\n'TSN', 'X',\n'TSP', 'X',\n'TSU', 'X',\n'\
TTA', 'X',\n'TTE', 'X',\n'TTN', 'X',\n'TTO', 'X',\\
n'TTP', 'X',\n'TTX', 'X',\n'TXL', 'X',\n'TYA', 'Y'\
,\n'TYB', 'Y',\n'TYD', 'X',\n'TYI', 'Y',\n'TYL', '\
X',\n'TYM', 'W',\n'TYN', 'Y',\n'TYQ', 'Y',\n'TYR',\
 'Y',\n'TYS', 'Y',\n'TYV', 'X',\n'TYY', 'A',\n'TZB\
', 'X',\n'TZC', 'X',\n'TZE', 'X',\n'TZL', 'X',\n'T\
ZO', 'X',\n'TZP', 'X',\n'U01', 'X',\n'U02', 'X',\n\
'U03', 'X',\n'U04', 'X',\n'U05', 'X',\n'U0E', 'X',\
\n'U10', 'X',\n'U18', 'X',\n'U2G', 'X',\n'U3P', 'X\
',\n'U49', 'X',\n'U55', 'X',\n'U5P', 'X',\n'U66', \
'X',\n'U89', 'X',\n'U8U', 'X',\n'UAA', 'X',\n'UAG'\
, 'A',\n'UAP', 'X',\n'UAR', 'X',\n'UC1', 'X',\n'UC\
2', 'X',\n'UC3', 'X',\n'UC4', 'X',\n'UD1', 'X',\n'\
UD2', 'X',\n'UDP', 'X',\n'UDX', 'X',\n'UFG', 'X',\\
n'UFM', 'X',\n'UFP', 'X',\n'UGA', 'X',\n'UIN', 'X'\
,\n'UKP', 'A',\n'UM3', 'X',\n'UMA', 'A',\n'UMG', '\
X',\n'UMP', 'X',\n'UNA', 'X',\n'UND', 'X',\n'UNI',\
 'X',\n'UNK', 'X',\n'UNN', 'X',\n'UNX', 'X',\n'UP5\
', 'X',\n'UP6', 'X',\n'UPA', 'X',\n'UPF', 'X',\n'U\
PG', 'X',\n'UPP', 'X',\n'UQ1', 'X',\n'UQ2', 'X',\n\
'UQ6', 'X',\n'UR2', 'X',\n'URA', 'X',\n'URE', 'X',\
\n'URF', 'X',\n'URI', 'X',\n'URS', 'X',\n'UTP', 'X\
',\n'UVC', 'X',\n'UVW', 'X',\n'V35', 'X',\n'V36', \
'X',\n'V4O', 'X',\n'V7O', 'X',\n'VAA', 'V',\n'VAC'\
, 'X',\n'VAD', 'V',\n'VAF', 'V',\n'VAG', 'X',\n'VA\
L', 'V',\n'VAN', 'X',\n'VAS', 'X',\n'VAX', 'X',\n'\
VDX', 'X',\n'VDY', 'X',\n'VG1', 'X',\n'VIB', 'X',\\
n'VIR', 'X',\n'VIT', 'X',\n'VK3', 'X',\n'VO3', 'X'\
,\n'VO4', 'X',\n'VS1', 'F',\n'VS2', 'F',\n'VS3', '\
F',\n'VS4', 'F',\n'VXA', 'X',\n'W01', 'X',\n'W02',\
 'X',\n'W03', 'X',\n'W11', 'X',\n'W33', 'X',\n'W35\
', 'X',\n'W42', 'X',\n'W43', 'X',\n'W54', 'X',\n'W\
56', 'X',\n'W59', 'X',\n'W71', 'X',\n'W84', 'X',\n\
'W8R', 'X',\n'W91', 'X',\n'WAY', 'X',\n'WCC', 'X',\
\n'WO2', 'X',\n'WO4', 'X',\n'WRB', 'X',\n'WRR', 'X\
',\n'WRS', 'X',\n'WW7', 'X',\n'X2F', 'X',\n'X7O', \
'X',\n'XAA', 'X',\n'XAN', 'X',\n'XAO', 'X',\n'XBB'\
, 'X',\n'XBP', 'X',\n'XDN', 'X',\n'XDP', 'X',\n'XI\
F', 'X',\n'XIM', 'X',\n'XK2', 'X',\n'XL1', 'X',\n'\
XLS', 'X',\n'XMP', 'X',\n'XN1', 'X',\n'XN2', 'X',\\
n'XN3', 'X',\n'XUL', 'X',\n'XV6', 'X',\n'XYD', 'X'\
,\n'XYH', 'X',\n'XYL', 'X',\n'XYP', 'X',\n'XYS', '\
X',\n'YOF', 'Y',\n'YRR', 'X',\n'YT3', 'X',\n'YZ9',\
 'X',\n'Z34', 'G',\n'Z5A', 'X',\n'ZAF', 'X',\n'ZAP\
', 'X',\n'ZEB', 'X',\n'ZEN', 'X',\n'ZES', 'X',\n'Z\
ID', 'X',\n'ZMR', 'X',\n'ZN3', 'X',\n'ZNH', 'X',\n\
'ZNO', 'X',\n'ZO3', 'X',\n'ZPR', 'P',\n'ZRA', 'A',\
\n'ZST', 'X',\n'ZYA', 'A',\n\n\n'ASN','N');\n} \n\\
n\nsub file2head\n      {\n	my $file = shift;\n	my\
 $size = shift;\n	my $f= new FileHandle;\n	my $lin\
e;\n	open ($f,$file);\n	read ($f,$line, $size);\n	\
close ($f);\n	return $line;\n      }\nsub file2tai\
l\n      {\n	my $file = shift;\n	my $size = shift;\
\n	my $f= new FileHandle;\n	my $line;\n	\n	open ($\
f,$file);\n	seek ($f,$size*-1, 2);\n	read ($f,$lin\
e, $size);\n	close ($f);\n	return $line;\n      }\\
n\n\nsub vtmpnam\n      {\n	my $r=rand(100000);\n	\
my $f=\"file.$r.$$\";\n	while (-e $f)\n	  {\n	    \
$f=vtmpnam();\n	  }\n	push (@TMPFILE_LIST, $f);\n	\
return $f;\n      }\n\nsub myexit\n  {\n    my $co\
de=@_[0];\n    if ($CLEAN_EXIT_STARTED==1){return;\
}\n    else {$CLEAN_EXIT_STARTED=1;}\n    ### ONLY\
 BARE EXIT\n    exit ($code);\n  }\nsub set_error_\
lock\n    {\n      my $name = shift;\n      my $pi\
d=$$;\n\n      \n      &lock4tc ($$,\"LERROR\", \"\
LSET\", \"$$ -- ERROR: $name $PROGRAM\\n\");\n    \
  return;\n    }\nsub set_lock\n  {\n    my $pid=s\
hift;\n    my $msg= shift;\n    my $p=getppid();\n\
    &lock4tc ($pid,\"LLOCK\",\"LRESET\",\"$p$msg\\\
n\");\n  }\nsub unset_lock\n   {\n     \n    my $p\
id=shift;\n    &lock4tc ($pid,\"LLOCK\",\"LRELEASE\
\",\"\");\n  }\nsub shift_lock\n  {\n    my $from=\
shift;\n    my $to=shift;\n    my $from_type=shift\
;\n    my $to_type=shift;\n    my $action=shift;\n\
    my $msg;\n    \n    if (!&lock4tc($from, $from\
_type, \"LCHECK\", \"\")){return 0;}\n    $msg=&lo\
ck4tc ($from, $from_type, \"LREAD\", \"\");\n    &\
lock4tc ($from, $from_type,\"LRELEASE\", $msg);\n \
   &lock4tc ($to, $to_type, $action, $msg);\n    r\
eturn;\n  }\nsub isshellpid\n  {\n    my $p=shift;\
\n    if (!lock4tc ($p, \"LLOCK\", \"LCHECK\")){re\
turn 0;}\n    else\n      {\n	my $c=lock4tc($p, \"\
LLOCK\", \"LREAD\");\n	if ( $c=~/-SHELL-/){return \
1;}\n      }\n    return 0;\n  }\nsub isrootpid\n \
 {\n    if(lock4tc (getppid(), \"LLOCK\", \"LCHECK\
\")){return 0;}\n    else {return 1;}\n  }\nsub lo\
ck4tc\n	{\n	  my ($pid,$type,$action,$value)=@_;\n\
	  my $fname;\n	  my $host=hostname;\n	  \n	  if (\
$type eq \"LLOCK\"){$fname=\"$LOCKDIR/.$pid.$host.\
lock4tcoffee\";}\n	  elsif ( $type eq \"LERROR\"){\
 $fname=\"$LOCKDIR/.$pid.$host.error4tcoffee\";}\n\
	  elsif ( $type eq \"LWARNING\"){ $fname=\"$LOCKD\
IR/.$pid.$host.warning4tcoffee\";}\n	  \n	  if ($d\
ebug_lock)\n	    {\n	      print STDERR \"\\n\\t--\
-lock4tc(tcg): $action => $fname =>$value (RD: $LO\
CKDIR)\\n\";\n	    }\n\n	  if    ($action eq \"LCH\
ECK\") {return -e $fname;}\n	  elsif ($action eq \\
"LREAD\"){return file2string($fname);}\n	  elsif (\
$action eq \"LSET\") {return string2file ($value, \
$fname, \">>\");}\n	  elsif ($action eq \"LRESET\"\
) {return string2file ($value, $fname, \">\");}\n	\
  elsif ($action eq \"LRELEASE\") \n	    {\n	     \
 if ( $debug_lock)\n		{\n		  my $g=new FileHandle;\
\n		  open ($g, \">>$fname\");\n		  print $g \"\\n\
Destroyed by $$\\n\";\n		  close ($g);\n		  safe_s\
ystem (\"mv $fname $fname.old\");\n		}\n	      els\
e\n		{\n		  unlink ($fname);\n		}\n	    }\n	  retu\
rn \"\";\n	}\n	\nsub file2string\n	{\n	  my $file=\
@_[0];\n	  my $f=new FileHandle;\n	  my $r;\n	  op\
en ($f, \"$file\");\n	  while (<$f>){$r.=$_;}\n	  \
close ($f);\n	  return $r;\n	}\nsub string2file \n\
    {\n    my ($s,$file,$mode)=@_;\n    my $f=new \
FileHandle;\n    \n    open ($f, \"$mode$file\");\\
n    print $f  \"$s\";\n    close ($f);\n  }\n\nBE\
GIN\n    {\n      srand;\n    \n      $SIG{'SIGUP'\
}='signal_cleanup';\n      $SIG{'SIGINT'}='signal_\
cleanup';\n      $SIG{'SIGQUIT'}='signal_cleanup';\
\n      $SIG{'SIGILL'}='signal_cleanup';\n      $S\
IG{'SIGTRAP'}='signal_cleanup';\n      $SIG{'SIGAB\
RT'}='signal_cleanup';\n      $SIG{'SIGEMT'}='sign\
al_cleanup';\n      $SIG{'SIGFPE'}='signal_cleanup\
';\n      \n      $SIG{'SIGKILL'}='signal_cleanup'\
;\n      $SIG{'SIGPIPE'}='signal_cleanup';\n      \
$SIG{'SIGSTOP'}='signal_cleanup';\n      $SIG{'SIG\
TTIN'}='signal_cleanup';\n      $SIG{'SIGXFSZ'}='s\
ignal_cleanup';\n      $SIG{'SIGINFO'}='signal_cle\
anup';\n      \n      $SIG{'SIGBUS'}='signal_clean\
up';\n      $SIG{'SIGALRM'}='signal_cleanup';\n   \
   $SIG{'SIGTSTP'}='signal_cleanup';\n      $SIG{'\
SIGTTOU'}='signal_cleanup';\n      $SIG{'SIGVTALRM\
'}='signal_cleanup';\n      $SIG{'SIGUSR1'}='signa\
l_cleanup';\n\n\n      $SIG{'SIGSEGV'}='signal_cle\
anup';\n      $SIG{'SIGTERM'}='signal_cleanup';\n \
     $SIG{'SIGCONT'}='signal_cleanup';\n      $SIG\
{'SIGIO'}='signal_cleanup';\n      $SIG{'SIGPROF'}\
='signal_cleanup';\n      $SIG{'SIGUSR2'}='signal_\
cleanup';\n\n      $SIG{'SIGSYS'}='signal_cleanup'\
;\n      $SIG{'SIGURG'}='signal_cleanup';\n      $\
SIG{'SIGCHLD'}='signal_cleanup';\n      $SIG{'SIGX\
CPU'}='signal_cleanup';\n      $SIG{'SIGWINCH'}='s\
ignal_cleanup';\n      \n      $SIG{'INT'}='signal\
_cleanup';\n      $SIG{'TERM'}='signal_cleanup';\n\
      $SIG{'KILL'}='signal_cleanup';\n      $SIG{'\
QUIT'}='signal_cleanup';\n      \n      our $debug\
_lock=$ENV{\"DEBUG_LOCK\"};\n      \n      \n     \
 \n      \n      foreach my $a (@ARGV){$CL.=\" $a\\
";}\n      if ( $debug_lock ){print STDERR \"\\n\\\
n\\n********** START PG: $PROGRAM *************\\n\
\";}\n      if ( $debug_lock ){print STDERR \"\\n\\
\n\\n**********(tcg) LOCKDIR: $LOCKDIR $$ ********\
*****\\n\";}\n      if ( $debug_lock ){print STDER\
R \"\\n --- $$ -- $CL\\n\";}\n      \n	     \n    \
  \n      \n    }\nsub flush_error\n  {\n    my $m\
sg=shift;\n    return add_error ($EXIT_FAILURE,$$,\
 $$,getppid(), $msg, $CL);\n  }\nsub add_error \n \
 {\n    my $code=shift;\n    my $rpid=shift;\n    \
my $pid=shift;\n    my $ppid=shift;\n    my $type=\
shift;\n    my $com=shift;\n    \n    $ERROR_DONE=\
1;\n    lock4tc ($rpid, \"LERROR\",\"LSET\",\"$pid\
 -- ERROR: $type\\n\");\n    lock4tc ($$, \"LERROR\
\",\"LSET\", \"$pid -- COM: $com\\n\");\n    lock4\
tc ($$, \"LERROR\",\"LSET\", \"$pid -- STACK: $ppi\
d -> $pid\\n\");\n   \n    return $code;\n  }\nsub\
 add_warning \n  {\n    my $rpid=shift;\n    my $p\
id =shift;\n    my $command=shift;\n    my $msg=\"\
$$ -- WARNING: $command\\n\";\n    print STDERR \"\
$msg\";\n    lock4tc ($$, \"LWARNING\", \"LSET\", \
$msg);\n  }\n\nsub signal_cleanup\n  {\n    print \
dtderr \"\\n**** $$ (tcg) was killed\\n\";\n    &c\
leanup;\n    exit ($EXIT_FAILURE);\n  }\nsub clean\
_dir\n  {\n    my $dir=@_[0];\n    if ( !-d $dir){\
return ;}\n    elsif (!($dir=~/tmp/)){return ;}#sa\
fety check 1\n    elsif (($dir=~/\\*/)){return ;}#\
safety check 2\n    else\n      {\n	`rm -rf $dir`;\
\n      }\n    return;\n  }\nsub cleanup\n  {\n   \
 #print stderr \"\\n----tc: $$ Kills $PIDCHILD\\n\\
";\n    #kill (SIGTERM,$PIDCHILD);\n    my $p=getp\
pid();\n    $CLEAN_EXIT_STARTED=1;\n    \n    \n  \
  \n    if (&lock4tc($$,\"LERROR\", \"LCHECK\", \"\
\"))\n      {\n	my $ppid=getppid();\n	if (!$ERROR_\
DONE) \n	  {\n	    &lock4tc($$,\"LERROR\", \"LSET\\
", \"$$ -- STACK: $p -> $$\\n\");\n	    &lock4tc($\
$,\"LERROR\", \"LSET\", \"$$ -- COM: $CL\\n\");\n	\
  }\n      }\n    my $warning=&lock4tc($$, \"LWARN\
ING\", \"LREAD\", \"\");\n    my $error=&lock4tc($\
$,  \"LERROR\", \"LREAD\", \"\");\n    #release er\
ror and warning lock if root\n    \n    if (isroot\
pid() && ($warning || $error) )\n      {\n	\n	prin\
t STDERR \"**************** Summary *************\\
\n$error\\n$warning\\n\";\n\n	&lock4tc($$,\"LERROR\
\",\"RELEASE\",\"\");\n	&lock4tc($$,\"LWARNING\",\\
"RELEASE\",\"\");\n      } \n    \n    \n    forea\
ch my $f (@TMPFILE_LIST)\n      {\n	if (-e $f){unl\
ink ($f);} \n      }\n    foreach my $d (@TMPDIR_L\
IST)\n      {\n	clean_dir ($d);\n      }\n    #No \
More Lock Release\n    #&lock4tc($$,\"LLOCK\",\"LR\
ELEASE\",\"\"); #release lock \n\n    if ( $debug_\
lock ){print STDERR \"\\n\\n\\n********** END PG: \
$PROGRAM ($$) *************\\n\";}\n    if ( $debu\
g_lock ){print STDERR \"\\n\\n\\n**********(tcg) L\
OCKDIR: $LOCKDIR $$ *************\\n\";}\n  }\nEND\
 \n  {\n    \n    &cleanup();\n  }\n   \n\nsub saf\
e_system \n{\n  my $com=shift;\n  my $ntry=shift;\\
n  my $ctry=shift;\n  my $pid;\n  my $status;\n  m\
y $ppid=getppid();\n  if ($com eq \"\"){return 1;}\
\n  \n  \n\n  if (($pid = fork ()) < 0){return (-1\
);}\n  if ($pid == 0)\n    {\n      set_lock($$, \\
" -SHELL- $com (tcg)\");\n      exec ($com);\n    \
}\n  else\n    {\n      lock4tc ($$, \"LLOCK\", \"\
LSET\", \"$pid\\n\");#update parent\n      $PIDCHI\
LD=$pid;\n    }\n  if ($debug_lock){printf STDERR \
\"\\n\\t .... safe_system (fasta_seq2hmm)  p: $$ c\
: $pid COM: $com\\n\";}\n\n  waitpid ($pid,WTERMSI\
G);\n\n  shift_lock ($pid,$$, \"LWARNING\",\"LWARN\
ING\", \"LSET\");\n\n  if ($? == $EXIT_FAILURE || \
lock4tc($pid, \"LERROR\", \"LCHECK\", \"\"))\n    \
{\n      if ($ntry && $ctry <$ntry)\n	{\n	  add_wa\
rning ($$,$$,\"$com failed [retry: $ctry]\");\n	  \
lock4tc ($pid, \"LRELEASE\", \"LERROR\", \"\");\n	\
  return safe_system ($com, $ntry, ++$ctry);\n	}\n\
      elsif ($ntry == -1)\n	{\n	  if (!shift_lock \
($pid, $$, \"LERROR\", \"LWARNING\", \"LSET\"))\n	\
    {\n	      add_warning ($$,$$,\"$com failed\");\
\n	    }\n	  else\n	    {\n	      lock4tc ($pid, \\
"LRELEASE\", \"LERROR\", \"\");\n	    }\n	  return\
 $?;}\n      else\n	{\n	  if (!shift_lock ($pid,$$\
, \"LERROR\",\"LERROR\", \"LSET\"))\n	    {\n	    \
  myexit(add_error ($EXIT_FAILURE,$$,$pid,getppid(\
), \"UNSPECIFIED system\", $com));\n	    }\n	}\n  \
  }\n  return $?;\n}\n\nsub check_configuration \n\
    {\n      my @l=@_;\n      my $v;\n      foreac\
h my $p (@l)\n	{\n	  \n	  if   ( $p eq \"EMAIL\")\\
n	    { \n	      if ( !($EMAIL=~/@/))\n		{\n		add_\
warning($$,$$,\"Could Not Use EMAIL\");\n		myexit(\
add_error ($EXIT_FAILURE,$$,$$,getppid(),\"EMAIL\"\
,\"$CL\"));\n	      }\n	    }\n	  elsif( $p eq \"I\
NTERNET\")\n	    {\n	      if ( !&check_internet_c\
onnection())\n		{\n		  myexit(add_error ($EXIT_FAI\
LURE,$$,$$,getppid(),\"INTERNET\",\"$CL\"));\n		}\\
n	    }\n	  elsif( $p eq \"wget\")\n	    {\n	     \
 if (!&pg_is_installed (\"wget\") && !&pg_is_insta\
lled (\"curl\"))\n		{\n		  myexit(add_error ($EXIT\
_FAILURE,$$,$$,getppid(),\"PG_NOT_INSTALLED:wget\"\
,\"$CL\"));\n		}\n	    }\n	  elsif( !(&pg_is_insta\
lled ($p)))\n	    {\n	      myexit(add_error ($EXI\
T_FAILURE,$$,$$,getppid(),\"PG_NOT_INSTALLED:$p\",\
\"$CL\"));\n	    }\n	}\n      return 1;\n    }\nsu\
b pg_is_installed\n  {\n    my @ml=@_;\n    my $r,\
 $p, $m;\n    my $supported=0;\n    \n    my $p=sh\
ift (@ml);\n    if ($p=~/::/)\n      {\n	if (safe_\
system (\"perl -M$p -e 1\")==$EXIT_SUCCESS){return\
 1;}\n	else {return 0;}\n      }\n    else\n      \
{\n	$r=`which $p 2>/dev/null`;\n	if ($r eq \"\"){r\
eturn 0;}\n	else {return 1;}\n      }\n  }\n\n\n\n\
sub check_internet_connection\n  {\n    my $intern\
et;\n    my $tmp;\n    &check_configuration ( \"wg\
et\"); \n    \n    $tmp=&vtmpnam ();\n    \n    if\
     (&pg_is_installed    (\"wget\")){`wget www.go\
ogle.com -O$tmp >/dev/null 2>/dev/null`;}\n    els\
if  (&pg_is_installed    (\"curl\")){`curl www.goo\
gle.com -o$tmp >/dev/null 2>/dev/null`;}\n    \n  \
  if ( !-e $tmp || -s $tmp < 10){$internet=0;}\n  \
  else {$internet=1;}\n    if (-e $tmp){unlink $tm\
p;}\n\n    return $internet;\n  }\nsub check_pg_is\
_installed\n  {\n    my @ml=@_;\n    my $r=&pg_is_\
installed (@ml);\n    if (!$r && $p=~/::/)\n      \
{\n	print STDERR \"\\nYou Must Install the perl pa\
ckage $p on your system.\\nRUN:\\n\\tsudo perl -MC\
PAN -e 'install $pg'\\n\";\n      }\n    elsif (!$\
r)\n      {\n	myexit(flush_error(\"\\nProgram $p S\
upported but Not Installed on your system\"));\n  \
    }\n    else\n      {\n	return 1;\n      }\n  }\
\n\n\n","use Cwd;\nuse Env;\nuse File::Path;\nuse \
FileHandle;\nuse strict;\n\n\nour (%MODE, %PG, %EN\
V_SET, %SUPPORTED_OS);\n\n\nour $EXIT_SUCCESS=0;\n\
our $EXIT_FAILURE=1;\nour $INTERNET=0;\n\nour $CP=\
\"cp \"; #was causing a crash on MacOSX\nour $SILE\
NT=\">/dev/null 2>/dev/null\";\nour $WEB_BASE=\"ht\
tp://www.tcoffee.org\";\nour $TCLINKDB_ADDRESS=\"$\
WEB_BASE/Resources/tclinkdb.txt\";\nour $OS=get_os\
();\nour $ROOT=&get_root();\nour $CD=cwd();\nour $\
CDIR=$CD;\nour $HOME=$ENV{'HOME'};\n\nour $OSNAME=\
$ENV{'OSNAME'};\nour $OSARCH=$ENV{'OSARCH'};\nour \
$REPO_ROOT=\"\";\n\nour $TCDIR;\nour $TCCACHE;\nou\
r $TCTMP;\nour $TCM;\nour $TCMETHODS;\nour $TCPLUG\
INS;\nour $PLUGINS_DIR=\"\";\nour $INSTALL_DIR=\"\\
";\n\nour $CXX=\"g++\";\nour $CXXFLAGS=\"\";\n\nou\
r $CPP=\"g++\";\nour $CPPFLAGS=\"\";\n\nour $CC=\"\
gcc\";\nour $CFLAGS=\"\";\n\nour $FC=\"f77\";\nour\
 $FFLAGS=\"\";\n\nmy $install=\"all\";\nmy $defaul\
t_update_action=\"no_update\";\nmy @required_appli\
cations=(\"wget_OR_curl\");\nmy @smode=(\"all\", \\
"clean\", \"install\");\n\n&initialize_PG();\n\nmy\
 $cl=join( \" \", @ARGV);\nif ($#ARGV==-1 || ($cl=\
~/-h/) ||($cl=~/-H/) )\n  {\n     print \"\\n!!!!!\
!! ./install  t_coffee             --> installs t_\
coffee only\";\n     print \"\\n!!!!!!! ./install \
 all                  --> installs all the modes [\
mcoffee, expresso, psicoffee,rcoffee..]\";\n     p\
rint \"\\n!!!!!!! ./install  [mcoffee|rcoffee|..] \
--> installs the specified mode\";\n     print \"\\
\n!!!!!!! ./install  -h                   --> prin\
t usage\\n\\n\";\n     if ( $#ARGV==-1){exit ($EXI\
T_FAILURE);}\n   }\n     \nif (($cl=~/-h/) ||($cl=\
~/-H/) )\n  {\n    my $m;\n    print \"\\n\\n!!!!!\
!! advanced mode\\n\";\n    foreach $m ((keys (%MO\
DE)),@smode)\n      {\n	print \"!!!!!!!       ./in\
stall $m\\n\";\n      }\n    \n    print \"!!!!!!!\
 ./install [target:package|mode|] [-update|-force|\
-exec=dir|-dis=dir|-root|-tclinkdb=file|-] [CC=|FC\
C=|CXX=|CFLAGS=|CXXFLAGS=]\\n\";\n    print \"!!!!\
!!! ./install clean    [removes all executables]\\\
n\";\n    print \"!!!!!!! ./install [optional:targ\
et] -update               [updates package already\
 installed]\\n\";\n    print \"!!!!!!! ./install [\
optional:target] -force                [Forces rec\
ompilation over everything]\\n\";\n    \n    print\
 \"!!!!!!! ./install [optional:target] -root      \
           [You are running as root]\\n\";\n    pr\
int \"!!!!!!! ./install [optional:target] -exec=/f\
oo/bar/       [address for the T-Coffee executable\
]\\n\";\n    print \"!!!!!!! ./install [optional:t\
arget] -dis=/foo/bar/        [Address where distri\
butions should be stored]\\n\";\n    print \"!!!!!\
!! ./install [optional:target] -tclinkdb=foo|updat\
e  [file containing all the packages to be install\
ed]\\n\";\n    print \"!!!!!!! ./install [optional\
:target] -clean                [clean everything]\\
\n\";\n    print \"!!!!!!! ./install [optional:tar\
get] -plugins              [plugins directory]\\n\\
";\n    print \"!!!!!!! ./install [optional:target\
] -tcdir=/foor/bar      [base path where T-Coffee \
will be installed]\\n\";\n    print \"!!!!!!! ./in\
stall [optional:target] -repo=/path/to/repo   [bin\
aries repository root directory]\\n\";\n    print \
\"!!!!!!! mode:\";\n    foreach $m (keys(%MODE)){p\
rint \"$m \";}\n    print \"\\n\";\n    print \"!!\
!!!!! Packages:\";\n    foreach $m (keys (%PG)){pr\
int \"$m \";}\n    print \"\\n\";\n    \n    print\
 \"\\n\\n\";\n    exit ($EXIT_FAILURE);\n  }\n\n\n\
\nmy (@argl)=($cl=~/(\\S+=[^=]+)\\s\\w+=/g);\npush\
 (@argl, ($cl=~/(\\S+=[^=]+\\S)\\s*$/g));\n\nforea\
ch $a (@argl)\n  {\n    if ( ($cl=~/CXX=(.*)/)){$C\
XX=$1;}\n    if ( ($cl=~/-CC=(.*)/    )){$CC=$1;}\\
n    if ( ($cl=~/-FC=(.*)/    )){$FC=$1;}\n    if \
( ($cl=~/-CFLAGS=(.*)/)){$CFLAGS=$1;}\n    if ( ($\
cl=~/-CXXFLAGS=(.*)/)){$CXXFLAGS=$1;}\n  }\nour ($\
ROOT_INSTALL, $NO_QUESTION, $default_update_action\
,$BINARIES_ONLY,$force, $default_update_action, $I\
NSTALL_DIR, $PLUGINS_DIR, $DISTRIBUTIONS,$tclinkdb\
, $proxy, $clean);\nif ( ($cl=~/-root/)){$ROOT_INS\
TALL=1;}\nif ( ($cl=~/-no_question/)){$NO_QUESTION\
=1;}\nif ( ($cl=~/-update/)){$default_update_actio\
n=\"update\";}\nif ( ($cl=~/-binaries/)){$BINARIES\
_ONLY=1;}\nif ( ($cl=~/-force/)){$force=1;$default\
_update_action=\"update\"}\nif ( ($cl=~/-exec=\\s*\
(\\S+)/)){$INSTALL_DIR=$1;}\nif ( ($cl=~/-plugins=\
\\s*(\\S+)/)){$PLUGINS_DIR=$1;}\nif ( ($cl=~/-dis=\
\\s*(\\S+)/)){$DISTRIBUTIONS=$1;}\n\nif ( ($cl=~/-\
tclinkdb=\\s*(\\S+)/)){$tclinkdb=$1;}\nif ( ($cl=~\
/-proxy=\\s*(\\S+)/)){$proxy=$1;}\nif ( ($cl=~/-cl\
ean/)){$clean=1;}\nif ( ($cl=~/-repo=\\s*(\\S+)/))\
{ $REPO_ROOT=$1; }\nif ( ($cl=~/-tcdir=\\s*(\\S+)/\
)){ $TCDIR=$1; }\nif ($tclinkdb){&update_tclinkdb \
($tclinkdb);}\n\n\nif( $REPO_ROOT ne \"\" ) {\n	if\
( $OSNAME eq \"\" ) { print \"You have specified t\
he repository folder but the required \\\"OSNAME\\\
\" enviroment variable is missing. \\n\"; exit 1; \
} \n	if( $OSARCH eq \"\" ) { print \"You have spec\
ified the repository folder but the required \\\"O\
SARCH\\\" enviroment variable is missing. \\n\"; e\
xit 1; } \n}\n\n\nif(!$TCDIR) { $TCDIR=\"$HOME/.t_\
coffee\"; }\n&add_dir ($TCDIR);\n&add_dir ($TCCACH\
E=\"$TCDIR/cache\");\n&add_dir ($TCTMP=\"$CDIR/tmp\
\");\n&add_dir ($TCM=\"$TCDIR/mcoffee\");\n&add_di\
r ($TCMETHODS=\"$TCDIR/methods\");\n&add_dir ($TCP\
LUGINS=\"$TCDIR/plugins/$OS\");\n\n\nour $BASE=\"$\
CD/bin\";\nour $BIN=\"$BASE/binaries/$OS\";\nour $\
DOWNLOAD_DIR=\"$BASE/download\";\nour $DOWNLOAD_FI\
LE=\"$DOWNLOAD_DIR/files\";\nour $TMP=\"$BASE/tmp\\
";\n\n&add_dir($BASE);\n&add_dir($BIN);\n&add_dir(\
$DOWNLOAD_DIR);\n&add_dir($DOWNLOAD_FILE);\nif (!$\
DISTRIBUTIONS){$DISTRIBUTIONS=\"$DOWNLOAD_DIR/dist\
ributions\";}\n&add_dir ($DISTRIBUTIONS);\n&add_di\
r ($TMP);\n\n\nif    (!$PLUGINS_DIR && !$ROOT_INST\
ALL){$PLUGINS_DIR=$TCPLUGINS;}\nelsif (!$PLUGINS_D\
IR &&  $ROOT_INSTALL){$PLUGINS_DIR=\"/usr/local/bi\
n/\";}\n\nif    (!$INSTALL_DIR && !$ROOT_INSTALL){\
$INSTALL_DIR=\"$HOME/bin/\";mkpath ($INSTALL_DIR);\
}\nelsif (!$INSTALL_DIR &&  $ROOT_INSTALL){$INSTAL\
L_DIR=\"/usr/local/bin/\";}\n\nif (-d \"mcoffee\")\
{`cp mcoffee/* $TCM`;}\n\n\nour $ENV_FILE=\"$TCDIR\
/t_coffee_env\";\n&env_file2putenv ($ENV_FILE);\n&\
set_proxy($proxy);\nmy ($target, $p, $r);\n$target\
=$p;\n\nforeach $p (  ((keys (%PG)),(keys(%MODE)),\
(@smode)) )\n  {\n    if ($ARGV[0] eq $p && $targe\
t eq \"\"){$target=$p;}\n  }\nif ($target eq \"\")\
{exit ($EXIT_FAILURE);}\n\n\nforeach $r (@required\
_applications)\n  {\n    my @app_list;\n    my $i;\
\n    $i=0;\n    \n    @app_list=split (/_OR_/, $r\
);\n    foreach my $pg (@app_list)\n      {\n	$i+=\
&pg_is_installed ($pg);\n      }\n    if ($i==0)\n\
      {\n      print \"One of the following packag\
es must be installed to proceed: \";\n      foreac\
h my $pg (@app_list)\n	{\n	  print (\"$pg \");\n	}\
\n      die;\n    }\n  }\n\n\n\n\n\n\n&sign_licens\
e_ni();\n\n\n$PG{C}{compiler}=get_C_compiler($CC);\
\n$PG{Fortran}{compiler}=get_F_compiler($FC);\n$PG\
{CXX}{compiler}=$PG{CPP}{compiler}=$PG{GPP}{compil\
er}=get_CXX_compiler($CXX);\nif ($CXXFLAGS){$PG{CP\
P}{options}=$PG{GPP}{options}=$PG{CXX}{options}=$C\
XXFLAGS;}\nif ($CFLAGS){$PG{C}{options}=$CFLAGS;}\\
nforeach my $c (keys(%PG))\n  {\n    my $arguments\
;\n    if ($PG{$c}{compiler})\n      {\n	$argument\
s=\"$PG{$c}{compiler_flag}=$PG{$c}{compiler} \";\n\
	if ($PG{$c}{options})\n	  {\n	    $arguments.=\"$\
PG{$c}{options_flag}=$PG{$c}{options} \";\n	  }\n	\
$PG{$c}{arguments}=$arguments;\n      }\n  }\n\nif\
 ($PG{$target}){$PG{$target}{install}=1;}\nelse\n \
 {\n    foreach my $pg (keys(%PG))\n      {\n	if (\
 $target eq \"all\" || ($PG{$pg}{mode}=~/$target/)\
)\n	  {\n	    $PG{$pg} {install}=1;\n	  }\n      }\
\n  }\n\nforeach my $pg (keys(%PG))\n  {\n    if (\
!$PG{$pg}{update_action}){$PG{$pg}{update_action}=\
$default_update_action;}\n    elsif ($PG{$pg}{upda\
te_action} eq \"never\"){$PG{$pg}{install}=0;}\n  \
  if ( $force && $PG{$pg}{install})\n      {\n	`rm\
 $BIN/$pg $BIN/$pg.exe $SILENT`;\n      }\n    if \
($PG{$pg}{update_action} eq \"update\" && $PG{$pg}\
{install}){$PG{$pg}{update}=1;}\n  }\n\nif (($targ\
et=~/clean/))\n  {\n    print \"------- cleaning e\
xecutables -----\\n\";\n    `rm bin/* $SILENT`;\n \
   exit ($EXIT_SUCCESS);\n  }\n\nif ( !$PG{$target\
}){print \"------- Installing T-Coffee Modes\\n\";\
}\n\nforeach my $m (keys(%MODE))\n  {\n    if ( $t\
arget eq \"all\" || $target eq $m)\n      {\n	prin\
t \"\\n------- The installer will now install the \
$m components $MODE{$m}{description}\\n\";\n	forea\
ch my $pg (keys(%PG))\n	  {\n	    if ( $PG{$pg}{mo\
de} =~/$m/ && $PG{$pg}{install})\n	      {\n		if (\
$PG{$pg}{touched}){print \"------- $PG{$pg}{dname}\
: already processed\\n\";}\n		else {$PG{$pg}{succe\
ss}=&install_pg($pg);$PG{$pg}{touched}=1;}\n	     \
 }\n	  }\n      }\n  }\n\nif ( $PG{$target}){print\
 \"------- Installing Individual Package\\n\";}\nf\
oreach my $pg (keys (%PG))\n  {\n    \n    if ( $P\
G{$pg}{install} && !$PG{$pg}{touched})\n      {\n	\
print \"\\n------- Install $pg\\n\";\n	$PG{$pg}{su\
ccess}=&install_pg($pg);$PG{$pg}{touched}=1;\n    \
  }\n  }\nprint \"------- Finishing The installati\
on\\n\";\nmy $final_report=&install ($INSTALL_DIR)\
;\n\nprint \"\\n\";\nprint \"*********************\
************************************************\\\
n\";\nprint \"********              INSTALLATION S\
UMMARY          *****************\\n\";\nprint \"*\
**************************************************\
******************\\n\";\nprint \"------- SUMMARY \
package Installation:\\n\";\nprint \"-------   Exe\
cutable Installed in: $PLUGINS_DIR\\n\";\n\nforeac\
h my $pg (keys(%PG))\n  {\n    if ( $PG{$pg}{insta\
ll})\n      {\n	my $bin_status=($PG{$pg}{from_bina\
ry} && $PG{$pg}{success})?\"[from binary]\":\"\";\\
n	if     ( $PG{$pg}{new} && !$PG{$pg}{old})       \
              {print \"*------        $PG{$pg}{dna\
me}: installed $bin_status\\n\"; $PG{$pg}{status}=\
1;}\n	elsif  ( $PG{$pg}{new} &&  $PG{$pg}{old})   \
                  {print \"*------        $PG{$pg}\
{dname}: updated $bin_status\\n\"  ; $PG{$pg}{stat\
us}=1;} \n	elsif  (!$PG{$pg}{new} &&  $PG{$pg}{old\
} && !$PG{$pg}{update}){print \"*------        $PG\
{$pg}{dname}: previous\\n\" ; $PG{$pg}{status}=1;}\
\n	elsif  (!$PG{$pg}{new} &&  $PG{$pg}{old} &&  $P\
G{$pg}{update}){print \"*------        $PG{$pg}{dn\
ame}: failed update (previous installation availab\
le)\\n\";$PG{$pg}{status}=0;}\n	else              \
                                            {print\
 \"*------        $PG{$pg}{dname}: failed installa\
tion\\n\";$PG{$pg}{status}=0;}\n      }\n  }\nmy $\
failure;\n\nif ( !$PG{$target}){print \"*------ SU\
MMARY mode Installation:\\n\";}\nforeach my $m (ke\
ys(%MODE))\n  {\n  \n    if ( $target eq \"all\" |\
| $target eq $m)\n      {\n	my $succesful=1;\n	for\
each my $pg (keys(%PG))\n	  {\n	    if (($PG{$pg}{\
mode}=~/$m/) && $PG{$pg}{install} && $PG{$pg}{stat\
us}==0)\n	      {\n		$succesful=0;\n		print \"*!!!\
!!!       $PG{$pg}{dname}: Missing\\n\";\n	      }\
\n	  }\n	if ( $succesful)\n	  {\n	    $MODE{$m}{st\
atus}=1;\n	    print \"*------       MODE $MODE{$m\
}{dname} SUCCESSFULLY installed\\n\";\n	  }\n	else\
\n	  {\n	    $failure++;\n	    $MODE{$m}{status}=0\
;\n	    print \"*!!!!!!       MODE $MODE{$m}{dname\
} UNSUCCESSFULLY installed\\n\";\n	  }\n      }\n \
 }\n\n    \n      \nif ($clean==1 && ($BASE=~/inst\
all4tcoffee/) ){print \"*------ Clean Installation\
 Directory: $BASE\\n\";`rm -rf $BASE`;}\nforeach m\
y $pg (keys(%PG)){if ($PG{$pg}{install} && $PG{$pg\
}{status}==0){exit ($EXIT_FAILURE);}}\n\nif ($fail\
ure)\n  {\n    print \"***************************\
******************************************\\n\";\n\
    print \"********     SOME PACKAGES FAILED TO I\
NSTALL        *****************\\n\";\n    print \\
"*************************************************\
********************\\n\";\n    print \"\\nSome of\
 the reported failures may be due to connectivity \
problems\";\n    print \"\\nRerun the installation\
 and the installer will specifically try to instal\
l the missing packages\";\n    print \"\\nIf this \
Fails, go to the original website and install the \
package manually\";\n  }\n\nprint \"**************\
**************************************************\
*****\\n\";\nprint \"********              FINALIZ\
E YOUR INSTALLATION    *****************\\n\";\npr\
int \"********************************************\
*************************\\n\";\nprint \"------- Y\
our executables are in:\\n\"; \nprint \"-------   \
    $PLUGINS_DIR:\\n\";\nprint \"------- Add this \
directory to your path with the following command:\
\\n\";\nprint \"-------       export PATH=$PLUGINS\
_DIR:\\$PATH\\n\";\nprint \"------- Make this perm\
anent by adding this line to the file:\\n\";\nprin\
t \"-------       $HOME/.bashrc\\n\";\nexit ($EXIT\
_SUCCESS);  \n  \nsub get_CXX_compiler\n  {\n    m\
y $c=@_[0];\n    my (@clist)=(\"g++\");\n    \n   \
 return get_compil ($c, @clist);\n }\nsub get_C_co\
mpiler\n  {\n    my $c=@_[0];\n    my (@clist)=(\"\
gcc\", \"cc\", \"icc\");\n    \n    return get_com\
pil ($c, @clist);\n }\n\nsub get_F_compiler\n  {\n\
    my ($c)=@_[0];\n    my @clist=(\"f77\", \"g77\\
",\"g95\", \"gfortran\", \"ifort\");\n    return g\
et_compil ($c, @clist);\n  } \n       \nsub get_co\
mpil\n  {\n    my ($fav,@clist)=(@_);\n    \n    #\
return the first compiler found installed in the s\
ystem. Check first the favorite\n    foreach my $c\
 ($fav,@clist)\n      {\n	if  (&pg_is_installed ($\
c)){return $c;}\n      }\n    return \"\";\n  }\ns\
ub exit_if_pg_not_installed\n  {\n    my (@arg)=(@\
_);\n    \n    foreach my $p (@arg)\n      {\n	if \
( !&pg_is_installed ($p))\n	  {\n	    print \"!!!!\
!!!! The $p utility must be installed for this ins\
tallation to proceed [FATAL]\\n\";\n	    die;\n	  \
}\n      }\n    return 1;\n  }\nsub set_proxy\n  {\
\n    my ($proxy)=(@_);\n    my (@list,$p);\n    \\
n    @list= (\"HTTP_proxy\", \"http_proxy\", \"HTT\
P_PROXY\", \"ALL_proxy\", \"all_proxy\",\"HTTP_pro\
xy_4_TCOFFEE\",\"http_proxy_4_TCOFFEE\");\n    \n \
   if (!$proxy)\n      {\n	foreach my $p (@list)\n\
	  {\n	    if ( ($ENV_SET{$p}) || $ENV{$p}){$proxy\
=$ENV{$p};}\n	  }\n      }\n    foreach my $p(@lis\
t){$ENV{$p}=$proxy;}\n  }\n	\nsub check_internet_c\
onnection\n  {\n    my $internet;\n    \n    if ( \
-e \"x\"){unlink (\"x\");}\n    if     (&pg_is_ins\
talled    (\"wget\")){`wget www.google.com -Ox >/d\
ev/null 2>/dev/null`;}\n    elsif  (&pg_is_install\
ed    (\"curl\")){`curl www.google.com -ox >/dev/n\
ull 2>/dev/null`;}\n    else\n      {\n	printf std\
err \"\\nERROR: No pg for remote file fetching [wg\
et or curl][FATAL]\\n\";\n	exit ($EXIT_FAILURE);\n\
      }\n    \n    if ( !-e \"x\" || -s \"x\" < 10\
){$internet=0;}\n    else {$internet=1;}\n    if (\
-e \"x\"){unlink \"x\";}\n    return $internet;\n \
 }\nsub url2file\n  {\n    my ($cmd, $file,$wget_a\
rg, $curl_arg)=(@_);\n    my ($exit,$flag, $pg, $a\
rg);\n    \n    if ($INTERNET || check_internet_co\
nnection ()){$INTERNET=1;}\n    else\n      {\n	pr\
int STDERR \"ERROR: No Internet Connection [FATAL:\
install.pl]\\n\";\n	exit ($EXIT_FAILURE);\n      }\
\n    \n    if     (&pg_is_installed    (\"wget\")\
){$pg=\"wget\"; $flag=\"-O\";$arg=\"--tries=2 --co\
nnect-timeout=10 $wget_arg\";}\n    elsif  (&pg_is\
_installed    (\"curl\")){$pg=\"curl\"; $flag=\"-o\
\";$arg=$curl_arg;}\n    else\n      {\n	printf st\
derr \"\\nERROR: No pg for remote file fetching [w\
get or curl][FATAL]\\n\";\n	exit ($EXIT_FAILURE);\\
n      }\n    \n    \n    if (-e $file){unlink($fi\
le);}\n    $exit=system \"$pg $cmd $flag$file $arg\
\";\n    return $exit;\n  }\n\nsub pg_is_installed\
\n  {\n    my ($p, $dir)=(@_);\n    my ($r,$m, $re\
t);\n    my ($supported, $language, $compil);\n   \
 \n  \n    if ( $PG{$p})\n      {\n	$language=$PG{\
$p}{language2};\n	$compil=$PG{$language}{compiler}\
;\n      }\n    \n    if ( $compil eq \"CPAN\")\n \
     {\n	if ( system (\"perl -M$p -e 1\")==$EXIT_S\
UCCESS){$ret=1;}\n	else {$ret=0;}\n      }\n    el\
sif ($dir)\n      {\n	if (-e \"$dir/$p\" || -e \"$\
dir/$p\\.exe\"){$ret=1;}\n	else {$ret=0;}\n      }\
\n    elsif (-e \"$PLUGINS_DIR/$p\" || -e \"$PLUGI\
NS_DIR/$p.exe\"){$ret=1;}\n    else\n      {\n	$r=\
`which $p 2>/dev/null`;\n	if ($r eq \"\"){$ret=0;}\
\n	else {$ret=1;}\n      }\n   \n    return $ret;\\
n  }\nsub install\n  {\n    my ($new_bin)=(@_);\n \
   my ($copied, $report);\n\n    \n    if (!$ROOT_\
INSTALL)\n      {\n	\n	if (-e \"$BIN/t_coffee\"){`\
$CP $BIN/t_coffee $INSTALL_DIR`};\n	`cp $BIN/* $PL\
UGINS_DIR`;\n	$copied=1;\n      }\n    else\n     \
 {\n	$copied=&root_run (\"You must be root to fina\
lize the installation\", \"$CP $BIN/* $INSTALL_DIR\
 $SILENT\");\n      }\n    \n     \n  if ( !$copie\
d)\n    {\n      $report=\"*!!!!!! Installation un\
succesful. The executables have been left in $BASE\
/bin\\n\";\n    }\n  elsif ( $copied && $ROOT)\n  \
  {\n      $report=\"*------ Installation succesfu\
l. Your executables have been copied in $new_bin a\
nd are on your PATH\\n\";\n    }\n  elsif ( $copie\
d && !$ROOT)\n    {\n      $report= \"*!!!!!! T-Co\
ffee and associated packages have been copied in: \
$new_bin\\n\";\n      $report.=\"*!!!!!! This addr\
ess is NOT in your PATH sytem variable\\n\";\n    \
  $report.=\"*!!!!!! You can do so by adding the f\
ollowing line in your ~/.bashrc file:\\n\";\n     \
 $report.=\"*!!!!!! export PATH=$new_bin:\\$PATH\\\
n\";\n    }\n  return $report;\n}\n\nsub sign_lice\
nse_ni\n  {\n    my $F=new FileHandle;\n    open (\
$F, \"license.txt\");\n    while (<$F>)\n      {\n\
	print \"$_\";\n      }\n    close ($F);\n    \n  \
  return;\n  }\n\nsub install_pg\n  {\n    my ($pg\
)=(@_);\n    my ($report, $previous, $language, $c\
ompiler, $return);\n    \n    if (!$PG{$pg}{instal\
l}){return 1;}\n    \n    $previous=&pg_is_install\
ed ($pg);\n    \n    if ($PG{$pg}{update_action} e\
q \"no_update\" && $previous)\n      {\n	$PG{$pg}{\
old}=1;\n	$PG{$pg}{new}=0;\n	$return=1;\n      }\n\
    else\n      {\n	$PG{$pg}{old}=$previous;\n	\n	\
if ($PG{$pg} {language2} eq \"Perl\"){&install_per\
l_package ($pg);}\n	elsif ($BINARIES_ONLY && &inst\
all_binary_package ($pg)){$PG{$pg}{from_binary}=1;\
}\n	elsif (&install_source_package ($pg)){;}\n	els\
e \n	  {\n	    \n	    if (!&supported_os($OS))\n	 \
     {\n		print \"!!!!!!!! $pg compilation failed,\
 binary unsupported for $OS\\n\"; \n	      }\n	   \
 elsif (!($PG{$pg}{from_binary}=&install_binary_pa\
ckage ($pg)))\n	      {\n		print \"!!!!!!!! $pg co\
mpilation and  binary installation failed\\n\";\n	\
      }\n	  }\n	$PG{$pg}{new}=$return=&pg_is_insta\
lled ($pg,$BIN);\n      }\n\n    \n    return $ret\
urn;\n  }\nsub install_perl_package\n  {\n    my (\
$pg)=(@_);\n    my ($report, $language, $compiler)\
;\n    \n    $language=$PG{$pg} {language2};\n    \
$compiler=$PG{$language}{compiler};\n    \n    if \
(!&pg_is_installed ($pg))\n      {\n	if ( $OS eq \\
"windows\"){`perl -M$compiler -e 'install $pg'`;}\\
n	elsif ( $ROOT eq \"sudo\"){system (\"sudo perl -\
M$compiler -e 'install $pg'\");}\n	else {system (\\
"su root -c perl -M$compiler -e 'install $pg'\");}\
\n      }\n    return &pg_is_installed ($pg);\n  }\
\n\n\n\nsub install_source_package\n  {\n    my ($\
pg)=(@_);\n    my ($report, $download, $arguments,\
 $language, $address, $name, $ext, $main_dir, $dis\
trib);\n    my $wget_tmp=\"$TMP/wget.tmp\";\n    m\
y (@fl);\n    if ( -e \"$BIN/$pg\" || -e \"$BIN/$p\
g.exe\"){return 1;}\n    \n    #\n    # check if t\
he module exists in the repository cache \n    #\n\
	if( repo_load($pg) ) {\n		return 1;\n	}\n    \n  \
  if ($pg eq \"t_coffee\")  {return   &install_t_c\
offee ($pg);}\n    elsif ($pg eq \"TMalign\"){retu\
rn   &install_TMalign ($pg);}\n    \n    chdir $DI\
STRIBUTIONS;\n    \n    $download=$PG{$pg}{source}\
;\n    \n    if (($download =~/tgz/))\n      {\n	(\
$address,$name,$ext)=($download=~/(.+\\/)([^\\/]+)\
(\\.tgz).*/);\n      }\n    elsif (($download=~/ta\
r\\.gz/))\n      {\n	($address,$name,$ext)=($downl\
oad=~/(.+\\/)([^\\/]+)(\\.tar\\.gz).*/);\n      }\\
n    elsif (($download=~/tar/))\n      {\n	($addre\
ss,$name,$ext)=($download=~/(.+\\/)([^\\/]+)(\\.ta\
r).*/);\n      }\n    else\n      {\n	($address,$n\
ame)=($download=~/(.+\\/)([^\\/]+)/);\n	$ext=\"\";\
\n      }\n    $distrib=\"$name$ext\";\n    \n    \
if ( !-d $pg){mkdir $pg;}\n    chdir $pg;\n   \n  \
  #get the distribution if available\n    if ( -e \
\"$DOWNLOAD_DIR/$distrib\")\n      {\n	`$CP $DOWNL\
OAD_DIR/$distrib .`;\n      }\n    #UNTAR and Prep\
are everything\n    if (!-e \"$name.tar\" && !-e \\
"$name\")\n      {\n	&check_rm ($wget_tmp);\n	prin\
t \"\\n------- Downloading/Installing $pg\\n\";\n	\
\n	if (!-e $distrib && &url2file (\"$download\", \\
"$wget_tmp\")==$EXIT_SUCCESS)\n	  {\n	    \n	    `\
mv $wget_tmp $distrib`;\n	    `$CP $distrib $DOWNL\
OAD_DIR/`;\n	  }\n\n	if (!-e $distrib)\n	  {\n	   \
 print \"!!!!!!! Download of $pg distribution fail\
ed\\n\";\n	    print \"!!!!!!! Check Address: $PG{\
$pg}{source}\\n\";\n	    return 0;\n	  }\n	print \\
"\\n------- unzipping/untaring $name\\n\";\n	if ((\
$ext =~/z/))\n	  { \n	    &flush_command (\"gunzip\
 $name$ext\");\n	    \n	  }\n	if (($ext =~/tar/) |\
| ($ext =~/tgz/))\n	  {\n	    &flush_command(\"tar\
 -xvf $name.tar\");\n	  }\n      }\n    #Guess and\
 enter the distribution directory\n    @fl=ls($p);\
\n    foreach my $f (@fl)\n      {\n	if (-d $f)\n	\
  {\n	    $main_dir=$f;\n	  }\n      }\n    if (-d\
 $main_dir)\n	  \n      {\n	chdir $main_dir;}\n   \
 else\n      {\n	print \"Error: $main_dir does not\
 exist\";\n      }\n    print \"\\n------- Compili\
ng/Installing $pg\\n\";\n    `make clean $SILENT`;\
\n    \n    \n    #\n    # SAP module\n    #\n    \
if ($pg eq \"sap\")\n      {\n	if (-e \"./configur\
e\")\n	  {\n	    #new sap distribution\n	    if ($\
OS eq \"macosx\")\n	      {\n		&replace_line_in_fi\
le (\"./src/galloc.h\", \"malloc.h\",  \"\");\n		&\
replace_line_in_file (\"./src/pdbprot.h\", \"mallo\
c.h\", \"\");\n		&replace_line_in_file (\"./src/pd\
bprot.c\", \"malloc.h\", \"\");\n	      }\n	    \n\
	    &flush_command (\"./configure\");\n	    &flus\
h_command (\"make clean\");\n	    &flush_command (\
\"make\");\n	    &check_cp (\"./src/$pg\", \"$BIN\\
");\n	    repo_store(\"./src/$pg\");\n	  }\n	else\\
n	  {\n	    #old style distribution\n	    `rm *.o \
sap  sap.exe ./util/aa/*.o  ./util/wt/.o $SILENT`;\
\n	    &flush_command (\"make $arguments sap\");\n\
	    &check_cp ($pg, \"$BIN\");\n	    repo_store($\
pg);\n	  }\n      }\n    \n    #\n    # CLUSTALW2 \
module\n    #\n    elsif ($pg eq \"clustalw2\")\n \
     {\n	&flush_command(\"./configure\");\n	&flush\
_command(\"make $arguments\");\n	&check_cp (\"./sr\
c/$pg\", \"$BIN\");\n	repo_store(\"./src/$pg\");\n\
      }\n    \n    #\n    # FSA module\n    # \n  \
  elsif ($pg eq \"fsa\")\n      {\n	&flush_command\
(\"./configure --prefix=$BIN\");\n	&flush_command(\
\"make $arguments\");\n	&flush_command (\"make ins\
tall\");\n\n	repo_store(\"fsa\", \"$BIN/bin\");\n	\
`mv $BIN/bin/* $BIN`;\n	`rmdir $BIN/bin`;\n      }\
\n    \n    #\n    # CLUSTALW module\n    #\n    e\
lsif ($pg eq \"clustalw\")\n      {\n	&flush_comma\
nd(\"make $arguments clustalw\");\n	`$CP $pg $BIN \
$SILENT`;\n	repo_store($pg);\n      }\n    \n    #\
\n    # MAFFT module\n    #\n    elsif ($pg eq \"m\
afft\")\n      {\n	my $base=cwd();\n	my $c;\n	\n	#\
compile core\n	mkpath (\"./mafft/bin\");\n	mkpath \
(\"./mafft/lib\");\n	chdir \"$base/core\";\n	`make\
 clean $SILENT`;\n	&flush_command (\"make $argumen\
ts\");\n	&flush_command (\"make install LIBDIR=../\
mafft/lib BINDIR=../mafft/bin\");\n	\n	#compile ex\
tension\n	chdir \"$base/extensions\";\n	`make clea\
n $SILENT`;\n	&flush_command (\"make $arguments\")\
;\n	&flush_command (\"make install LIBDIR=../mafft\
/lib BINDIR=../mafft/bin\");\n	\n	#put everything \
in mafft and copy the compiled stuff in bin\n	chdi\
r \"$base\";\n	if ($ROOT_INSTALL)\n	  {\n	    &roo\
t_run (\"You Must be Root to Install MAFFT\\n\", \\
"mkdir /usr/local/mafft/;$CP mafft/lib/* /usr/loca\
l/mafft;$CP mafft/lib/mafft* /usr/local/bin ;$CP m\
afft/bin/mafft /usr/local/bin/; \");\n	  }\n	else\\
n	  {\n	    `$CP mafft/lib/*  $BIN`;\n	    `$CP ma\
fft/bin/mafft  $BIN`;\n	  }\n	`tar -cvf mafft.tar \
mafft`;\n	`gzip mafft.tar`;\n	`mv mafft.tar.gz $BI\
N`;\n	\n	repo_store(\"mafft/bin/mafft\", \"mafft/l\
ib/\", \"$BIN/mafft.tar.gz\");\n      }\n      \n \
   #\n    # DIALIGN-TX module\n    #\n    elsif ( \
$pg eq \"dialign-tx\" )\n      {\n	my $f;\n	my $ba\
se=cwd();\n\n	chdir \"./source\";\n	if ($OS eq \"m\
acosx\"){&flush_command (\"cp makefile.MAC_OS make\
file\");}\n\n	&flush_command (\" make CPPFLAGS='-O\
3 -funroll-loops' all\");\n	\n	chdir \"..\";\n	&ch\
eck_cp (\"./source/$pg\", \"$BIN\");\n	repo_store(\
\"./source/$pg\");\n      }\n      \n    #\n    # \
DIALIGN-T module \n    # (is the same as dialign-t\
x, but it is mantained for backward name compatibi\
lity with tcoffee)\n    #\n    elsif ( $pg eq \"di\
align-t\" )\n      {\n	my $f;\n	my $base=cwd();\n\\
n	chdir \"./source\";\n	if ($OS eq \"macosx\"){&fl\
ush_command (\"cp makefile.MAC_OS makefile\");}\n\\
n	&flush_command (\" make CPPFLAGS='-O3 -funroll-l\
oops' all\");\n	\n	chdir \"..\";\n	&check_cp (\"./\
source/dialign-tx\", \"$BIN/dialign-t\");\n	repo_s\
tore(\"$BIN/dialign-t\");	\n      }      \n      \\
n    #\n    # POA module\n    #\n    elsif ($pg eq\
 \"poa\")\n      {\n	&flush_command (\"make $argum\
ents poa\");\n	&check_cp (\"$pg\", \"$BIN\");\n	re\
po_store(\"$pg\");\n      }\n     \n     \n    #\n\
    # PROBCONS module\n    #\n    elsif ( $pg eq \\
"probcons\")\n      {\n	&add_C_libraries(\"./Proba\
bilisticModel.h\", \"list\", \"cstring\");\n	\n	`r\
m *.exe $SILENT`;\n	&flush_command (\"make $argume\
nts probcons\");\n	&check_cp(\"$pg\", \"$BIN/$pg\"\
);\n	repo_store(\"$pg\");\n      }\n      \n    #\\
n    # PROBCONS RNA module\n    #\n    elsif ( $pg\
 eq \"probconsRNA\")\n      {\n	&add_C_libraries(\\
"./ProbabilisticModel.h\", \"list\", \"cstring\");\
\n	&add_C_libraries(\"./Main.cc\", \"iomanip\", \"\
cstring\",\"climits\");\n	`rm *.exe $SILENT`;\n	&f\
lush_command (\"make $arguments probcons\");\n	&ch\
eck_cp(\"probcons\", \"$BIN/$pg\");\n	repo_store(\\
"$BIN/$pg\");\n      }\n\n	#\n	# MUSCLE module\n	#\
\n    elsif (  $pg eq \"muscle\")\n      {	\n	`rm \
*.o muscle muscle.exe $SILENT`;\n	if ($OS eq \"mac\
osx\" || $OS eq \"linux\")\n	  {\n	    &replace_li\
ne_in_file (\"./Makefile\", \"LDLIBS = -lm -static\
\",  \"LDLIBS = -lm\");\n	  }\n	elsif ($OS eq \"wi\
ndows\")\n	  {\n	    &replace_line_in_file (\"./in\
tmath.cpp\",  \"double log2e\",      \"double cedr\
ic_log\");\n	    &replace_line_in_file (\"./intmat\
h.cpp\",  \"double log2\",       \"double log_notu\
se\");\n	    &replace_line_in_file (\"./intmath.cp\
p\",  \"double cedric_log\", \"double log2e\");\n	\
  }\n	&flush_command (\"make $arguments all\");\n	\
&check_cp(\"$pg\", \"$BIN\");\n	repo_store(\"$pg\"\
);	\n      }\n      \n     #\n     # MUS4 module\n\
     #\n     elsif (  $pg eq \"mus4\")\n      {\n	\
`rm *.o muscle muscle.exe $SILENT`;\n	&flush_comma\
nd (\"./mk\");\n	&check_cp(\"$pg\", \"$BIN\");\n	r\
epo_store(\"$pg\");	\n      }\n      \n    #\n    \
# PCMA module\n    #\n    elsif ( $pg eq \"pcma\")\
\n      {\n	if ($OS eq \"macosx\")\n	  {\n	    &re\
place_line_in_file (\"./alcomp2.c\", \"malloc.h\",\
  \"\");\n	  }\n	&flush_command (\"make $arguments\
 pcma\");\n	&check_cp(\"$pg\", \"$BIN\");\n	repo_s\
tore(\"$pg\");	\n      }\n      \n    #\n    # KAL\
IGN module\n    #\n    elsif ($pg eq \"kalign\")\n\
      {\n	&flush_command (\"./configure\");\n	&flu\
sh_command(\"make $arguments\");\n	&check_cp (\"$p\
g\",$BIN);\n	repo_store(\"$pg\");	\n      }\n     \
 \n    #\n    # AMAP module\n    #\n    elsif ( $p\
g eq \"amap\")\n      {\n	&add_C_libraries(\"./Ama\
p.cc\", \"iomanip\", \"cstring\",\"climits\");	\n	\
`make clean $SILENT`;\n	&flush_command (\"make $ar\
guments all\");\n	&check_cp (\"$pg\", $BIN);\n	rep\
o_store(\"$pg\");	\n      }\n      \n    #\n    # \
PRODA module\n    #\n    elsif ( $pg eq \"proda\")\
\n      {\n	&add_C_libraries(\"AlignedFragment.h\"\
, \"vector\", \"iostream\", \"cstring\",\"cstdlib\\
");\n	&add_C_libraries(\"Main.cc\", \"vector\", \"\
climits\");	\n	&add_C_libraries(\"Sequence.cc\", \\
"stdlib.h\", \"cstdio\");	\n	&flush_command (\"mak\
e $arguments all\");\n	&check_cp (\"$pg\", $BIN);\\
n	repo_store(\"$pg\");	\n      }\n      \n    #\n \
   # PRANK module\n    #\n    elsif ( $pg eq \"pra\
nk\")\n      {\n	&flush_command (\"make $arguments\
 all\");\n	&check_cp (\"$pg\", $BIN);\n	repo_store\
(\"$pg\");	\n      }\n      \n    #\n    # !!!! MU\
STANG module\n    #\n     elsif ( $pg eq \"mustang\
\")\n      {\n	&flush_command (\"rm ./bin/*\");\n	\
&flush_command (\"make $arguments all\");\n\n	if (\
 $OS=~/windows/){&flush_command(\"cp ./bin/* $BIN/\
mustang.exe\");}\n	else {&flush_command(\"cp ./bin\
/* $BIN/mustang\");}\n	\n	repo_store(\"$BIN/mustan\
g\");\n      }\n\n	#\n	# RNAplfold module\n	#\n   \
 elsif ( $pg eq \"RNAplfold\")\n      {\n	&flush_c\
ommand(\"./configure\");\n	&flush_command (\"make \
$arguments all\");\n	&check_cp(\"./Progs/RNAplfold\
\", \"$BIN\");\n	&check_cp(\"./Progs/RNAalifold\",\
 \"$BIN\");\n	&check_cp(\"./Progs/RNAfold\", \"$BI\
N\");\n	\n	repo_store(\"./Progs/RNAplfold\", \"./P\
rogs/RNAalifold\", \"./Progs/RNAfold\");\n      }\\
n      \n    #\n    # !!! RETREE module\n    #\n  \
  elsif ( $pg eq \"retree\")\n      {\n	chdir \"sr\
c\";\n	&flush_command (\"make $arguments all\");\n\
	&flush_command (\"make put\");\n	system \"cp ../e\
xe/* $BIN\";\n	\n	repo_store(\"retree\", \"../exe\\
");\n      }\n	\n    chdir $CDIR;\n    return &pg_\
is_installed ($pg, $BIN);\n  }\n\nsub install_t_co\
ffee\n  {\n    my ($pg)=(@_);\n    my ($report,$cf\
lags, $arguments, $language, $compiler) ;\n    #1-\
Install T-Coffee\n    chdir \"t_coffee_source\";\n\
    &flush_command (\"make clean\");\n    print \"\
\\n------- Compiling T-Coffee\\n\";\n    $language\
=$PG{$pg} {language2};\n    $arguments=$PG{$langua\
ge}{arguments};\n    if (!($arguments =~/CFLAGS/))\
{$arguments .= \" CFLAGS=-O2 \";}\n\n    if ( $CC \
ne \"\"){&flush_command (\"make -i $arguments t_co\
ffee\");}\n    &check_cp ($pg, $BIN);\n    \n    c\
hdir $CDIR;\n    return &pg_is_installed ($pg, $BI\
N);\n  }\nsub install_TMalign\n  {\n    my ($pg)=(\
@_);\n    my $report;\n    chdir \"t_coffee_source\
\";\n    print \"\\n------- Compiling TMalign\\n\"\
;\n    `rm TMalign TMalign.exe $SILENT`;\n    if (\
 $FC ne \"\"){&flush_command (\"make -i $PG{Fortra\
n}{arguments} TMalign\");}\n    &check_cp ($pg, $B\
IN);\n    repo_store($pg);\n\n    if ( !-e \"$BIN/\
$pg\" && pg_has_binary_distrib ($pg))\n      {\n	p\
rint \"!!!!!!! Compilation of $pg impossible. Will\
 try to install from binary\\n\";\n	return &instal\
l_binary_package ($pg);\n      }\n    chdir $CDIR;\
\n    return &pg_is_installed ($pg, $BIN);\n  }\n\\
nsub pg_has_binary_distrib\n  {\n    my ($pg)=(@_)\
;\n    if ($PG{$pg}{windows}){return 1;}\n    elsi\
f ($PG{$pg}{osx}){return 1;}\n    elsif ($PG{$pg}{\
linux}){return 1;}\n    return 0;\n  }\nsub instal\
l_binary_package\n  {\n    my ($pg)=(@_);\n    my \
($base,$report,$name, $download, $arguments, $lang\
uage, $dir);\n    my $isdir;\n    &input_os();\n  \
  \n    if (!&supported_os($OS)){return 0;}\n    i\
f ( $PG{$pg}{binary}){$name=$PG{$pg}{binary};}\n  \
  else \n      {\n	$name=$pg;\n	if ( $OS eq \"wind\
ows\"){$name.=\".exe\";}\n      }\n    \n    $down\
load=\"$WEB_BASE/Packages/Binaries/$OS/$name\";\n \
   \n    $base=cwd();\n    chdir $TMP;\n    \n    \
if (!-e $name)\n      {\n	`rm x $SILENT`;\n	if ( u\
rl2file(\"$download\",\"x\")==$EXIT_SUCCESS)\n	  {\
\n	    `mv x $name`;\n	  }\n      }\n    \n    if \
(!-e $name)\n      {\n	print \"!!!!!!! $PG{$pg}{dn\
ame}: Download of $pg binary failed\\n\";\n	print \
\"!!!!!!! $PG{$pg}{dname}: Check Address: $downloa\
d\\n\";\n	return 0;\n      }\n    print \"\\n-----\
-- Installing $pg\\n\";\n    \n    if ($name =~/ta\
r\\.gz/)\n      {\n	`gunzip  $name`;\n	`tar -xvf $\
pg.tar`;\n	chdir $pg;\n	if ( $pg eq \"mafft\")\n	 \
 {\n	    if ($ROOT_INSTALL)\n	      {\n		&root_run\
 (\"You Must be Roor to Install MAFFT\\n\", \"$CP \
mafft/bin/* /usr/local/mafft;mkdir /usr/local/maff\
t/; $CP mafft/lib/* /usr/local/bin/\");\n	      }\\
n	    else\n	      {\n		`$CP $TMP/$pg/bin/* $BIN $\
SILENT`;\n		`$CP $TMP/$pg/lib/* $BIN $SILENT`;\n	 \
     }\n	  }\n	else\n	  {\n	    if (-e \"$TMP/$pg/\
data\"){`$CP $TMP/$pg/data/* $TCM $SILENT`;}\n	   \
 if (!($pg=~/\\*/)){`rm -rf $pg`;}\n	  }\n      }\\
n    else\n      {\n	&check_cp (\"$pg\", \"$BIN\")\
;\n	`chmod u+x $BIN/$pg`; \n	unlink ($pg);\n      \
}\n    chdir $base;\n    $PG{$pg}{from_binary}=1;\\
n    return &pg_is_installed ($pg, $BIN);\n  }\n\n\
sub add_dir \n  {\n    my $dir=@_[0];\n    \n    i\
f (!-e $dir && !-d $dir)\n      {\n	my @l;\n	umask\
 (0000);\n	@l=mkpath ($dir,{mode => 0777});\n	\n  \
    }\n    else\n      {\n	return 0;\n      }\n  }\
\nsub check_rm \n  {\n    my ($file)=(@_);\n    \n\
    if ( -e $file)\n      {\n	return unlink($file)\
;\n      }\n    return 0;\n  }\nsub check_cp\n  {\\
n    my ($from, $to)=(@_);\n    if ( !-e $from && \
-e \"$from\\.exe\"){$from=\"$from\\.exe\";}\n    i\
f ( !-e $from){return 0;}\n        \n    `$CP $fro\
m $to`;\n    return 1;\n  }\n\nsub repo_store \n{\\
n   # check that all required data are available\n\
   if( $REPO_ROOT eq \"\" ) { return; }\n\n\n    #\
 extract the package name from the specified path\\
n    my $pg =`basename $_[0]`;\n    chomp($pg);\n	\
\n    my $VER = $PG{$pg}{version};\n    my $CACHE \
= \"$REPO_ROOT/$pg/$VER/$OSNAME-$OSARCH\"; \n    \\
n    print \"-------- Storing package: \\\"$pg\\\"\
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
 \"-------- Module \\\"$pg\\\" NOT found on reposi\
tory cache.\\n\";\n    	return 0;\n    }\n    \n  \
  print \"-------- Module \\\"$pg\\\" found on rep\
ository cache. Using copy on path: $CACHE\\n\";\n \
   `cp $CACHE/* $BIN`;\n    return 1;\n}\n\nsub ch\
eck_file_list_exists \n  {\n    my ($base, @flist)\
=(@_);\n    my $f;\n\n    foreach $f (@flist)\n   \
   {\n	if ( !-e \"$base/$f\"){return 0;}\n      }\\
n    return 1;\n  }\nsub ls\n  {\n    my $f=@_[0];\
\n    my @fl;\n    chomp(@fl=`ls -1 $f`);\n    ret\
urn @fl;\n  }\nsub flush_command\n  {\n    my $com\
mand=@_[0];\n    my $F=new FileHandle;\n    open (\
$F, \"$command|\");\n    while (<$F>){print \"    \
--- $_\";}\n    close ($F);\n  }    \n\nsub input_\
installation_directory\n  {\n    my $dir=@_[0];\n \
   my $new;\n    \n    print \"------- The current\
 installation directory is: [$dir]\\n\";\n    prin\
t \"??????? Return to keep the default or new valu\
e:\";\n   \n    if ($NO_QUESTION==0)\n      {\n	ch\
omp ($new=<stdin>);\n	while ( $new ne \"\" && !inp\
ut_yes (\"You have entered $new. Is this correct? \
([y]/n):\"))\n	  {\n	    print \"???????New instal\
lation directory:\";\n	    chomp ($new=<stdin>);\n\
	  }\n	$dir=($new eq \"\")?$dir:$new;\n	$dir=~s/\\\
/$//;\n      }\n    \n    if ( -d $dir){return $di\
r;}\n    elsif (&root_run (\"You must be root to c\
reate $dir\",\"mkdir $dir\")==$EXIT_SUCCESS){retur\
n $dir;}\n    else\n      {\n	print \"!!!!!!! $dir\
 could not be created\\n\";\n	if ( $NO_QUESTION)\n\
	  {\n	    return \"\";\n	  }\n	elsif ( &input_yes\
 (\"??????? Do you want to provide a new directory\
([y]/n)?:\"))\n	  {\n	    return input_installatio\
n_directory ($dir);\n	  }\n	else\n	  {\n	    retur\
n \"\";\n	  }\n      }\n    \n  }\nsub input_yes\n\
  {\n    my $question =@_[0];\n    my $answer;\n\n\
    if ($NO_QUESTION==1){return 1;}\n    \n    if \
($question eq \"\"){$question=\"??????? Do you wis\
h to proceed ([y]/n)?:\";}\n    print $question;\n\
    chomp($answer=lc(<STDIN>));\n    if (($answer=\
~/^y/) || $answer eq \"\"){return 1;}\n    elsif (\
 ($answer=~/^n/)){return 0;}\n    else\n      {\n	\
return input_yes($question);\n      }\n  }\nsub ro\
ot_run\n  {\n    my ($txt, $cmd)=(@_);\n    \n    \
if ( system ($cmd)==$EXIT_SUCCESS){return $EXIT_SU\
CCESS;}\n    else \n      {\n	print \"------- $txt\
\\n\";\n	if ( $ROOT eq \"sudo\"){return system (\"\
sudo $cmd\");}\n	else {return system (\"su root -c\
 \\\"$cmd\\\"\");}\n      }\n  }\nsub get_root\n  \
{\n    if (&pg_is_installed (\"sudo\")){return \"s\
udo\";}\n    else {return \"su\";}\n  }\n\nsub get\
_os\n  {\n    my $raw_os=`uname`;\n    my $os;\n\n\
    $raw_os=lc ($raw_os);\n    \n    if ($raw_os =\
~/cygwin/){$os=\"windows\";}\n    elsif ($raw_os =\
~/linux/){$os=\"linux\";}\n    elsif ($raw_os =~/o\
sx/){$os=\"macosx\";}\n    elsif ($raw_os =~/darwi\
n/){$os=\"macosx\";}\n    else\n      {\n	$os=$raw\
_os;\n      }\n    return $os;\n  }\nsub input_os\\
n  {\n    my $answer;\n    if ($OS) {return $OS;}\\
n    \n    print \"??????? which os do you use: [w\
]indows, [l]inux, [m]acosx:?\";\n    $answer=lc(<S\
TDIN>);\n\n    if (($answer=~/^m/)){$OS=\"macosx\"\
;}\n    elsif ( ($answer=~/^w/)){$OS=\"windows\";}\
\n    elsif ( ($answer=~/^linux/)){$OS=\"linux\";}\
\n    \n    else\n      {\n	return &input_os();\n \
     }\n    return $OS;\n  }\n\nsub supported_os\n\
  {\n    my ($os)=(@_[0]);\n    return $SUPPORTED_\
OS{$os};\n  }\n    \n    \n\n\nsub update_tclinkdb\
 \n  {\n    my $file =@_[0];\n    my $name;\n    m\
y $F=new FileHandle;\n    my ($download, $address,\
 $name, $l, $db);\n    \n    if ( $file eq \"updat\
e\"){$file=$TCLINKDB_ADDRESS;}\n    \n    if ( $fi\
le =~/http:\\/\\// || $file =~/ftp:\\/\\//)\n     \
 {\n	($address, $name)=($download=~/(.*)\\/([^\\/]\
+)$/);\n	`rm x $SILENT`;\n	if (&url2file ($file,\"\
x\")==$EXIT_SUCCESS)\n	  {\n	    print \"------- S\
usscessful upload of $name\";\n	    `mv x $name`;\\
n	    $file=$name;\n	  }\n      }\n    open ($F, \\
"$file\");\n    while (<$F>)\n      {\n	my $l=$_;\\
n	if (($l =~/^\\/\\//) || ($db=~/^#/)){;}\n	elsif \
( !($l =~/\\w/)){;}\n	else\n	  {\n	    my @v=split\
 (/\\s+/, $l);\n	    if ( $l=~/^MODE/)\n	      {\n\
		$MODE{$v[1]}{$v[2]}=$v[3];\n	      }\n	    elsif\
 ($l=~/^PG/)\n	      {\n		$PG{$v[1]}{$v[2]}=$v[3];\
\n	      }\n	  }\n      }\n    close ($F);\n    &p\
ost_process_PG();\n    return;\n  }\n\n\n\nsub ini\
tialize_PG\n  {\n\n$PG{\"t_coffee\"}{\"4_TCOFFEE\"\
}=\"TCOFFEE\";\n$PG{\"t_coffee\"}{\"type\"}=\"sequ\
ence_multiple_aligner\";\n$PG{\"t_coffee\"}{\"ADDR\
ESS\"}=\"http://www.tcoffee.org\";\n$PG{\"t_coffee\
\"}{\"language\"}=\"C\";\n$PG{\"t_coffee\"}{\"lang\
uage2\"}=\"C\";\n$PG{\"t_coffee\"}{\"source\"}=\"h\
ttp://www.tcoffee.org/Packages/T-COFFEE_distributi\
on.tar.gz\";\n$PG{\"t_coffee\"}{\"update_action\"}\
=\"always\";\n$PG{\"t_coffee\"}{\"mode\"}=\"tcoffe\
e,mcoffee,rcoffee,expresso,3dcoffee\";\n$PG{\"clus\
talw2\"}{\"4_TCOFFEE\"}=\"CLUSTALW2\";\n$PG{\"clus\
talw2\"}{\"type\"}=\"sequence_multiple_aligner\";\\
n$PG{\"clustalw2\"}{\"ADDRESS\"}=\"http://www.clus\
tal.org\";\n$PG{\"clustalw2\"}{\"language\"}=\"C++\
\";\n$PG{\"clustalw2\"}{\"language2\"}=\"CXX\";\n$\
PG{\"clustalw2\"}{\"source\"}=\"http://www.clustal\
.org/download/2.0.10/clustalw-2.0.10-src.tar.gz\";\
\n$PG{\"clustalw2\"}{\"mode\"}=\"mcoffee,rcoffee\"\
;\n$PG{\"clustalw2\"}{\"version\"}=\"2.0.10\";\n$P\
G{\"clustalw\"}{\"4_TCOFFEE\"}=\"CLUSTALW\";\n$PG{\
\"clustalw\"}{\"type\"}=\"sequence_multiple_aligne\
r\";\n$PG{\"clustalw\"}{\"ADDRESS\"}=\"http://www.\
clustal.org\";\n$PG{\"clustalw\"}{\"language\"}=\"\
C\";\n$PG{\"clustalw\"}{\"language2\"}=\"C\";\n$PG\
{\"clustalw\"}{\"source\"}=\"http://www.clustal.or\
g/download/1.X/ftp-igbmc.u-strasbg.fr/pub/ClustalW\
/clustalw1.82.UNIX.tar.gz\";\n$PG{\"clustalw\"}{\"\
mode\"}=\"mcoffee,rcoffee\";\n$PG{\"clustalw\"}{\"\
version\"}=\"1.82\";\n$PG{\"dialign-t\"}{\"4_TCOFF\
EE\"}=\"DIALIGNT\";\n$PG{\"dialign-t\"}{\"type\"}=\
\"sequence_multiple_aligner\";\n$PG{\"dialign-t\"}\
{\"ADDRESS\"}=\"http://dialign-tx.gobics.de/\";\n$\
PG{\"dialign-t\"}{\"DIR\"}=\"/usr/share/dialign-tx\
/\";\n$PG{\"dialign-t\"}{\"language\"}=\"C\";\n$PG\
{\"dialign-t\"}{\"language2\"}=\"C\";\n$PG{\"diali\
gn-t\"}{\"source\"}=\"http://dialign-tx.gobics.de/\
DIALIGN-TX_1.0.2.tar.gz\";\n$PG{\"dialign-t\"}{\"m\
ode\"}=\"mcoffee\";\n$PG{\"dialign-t\"}{\"binary\"\
}=\"dialign-t\";\n$PG{\"dialign-t\"}{\"version\"}=\
\"1.0.2\";\n$PG{\"dialign-tx\"}{\"4_TCOFFEE\"}=\"D\
IALIGNTX\";\n$PG{\"dialign-tx\"}{\"type\"}=\"seque\
nce_multiple_aligner\";\n$PG{\"dialign-tx\"}{\"ADD\
RESS\"}=\"http://dialign-tx.gobics.de/\";\n$PG{\"d\
ialign-tx\"}{\"DIR\"}=\"/usr/share/dialign-tx/\";\\
n$PG{\"dialign-tx\"}{\"language\"}=\"C\";\n$PG{\"d\
ialign-tx\"}{\"language2\"}=\"C\";\n$PG{\"dialign-\
tx\"}{\"source\"}=\"http://dialign-tx.gobics.de/DI\
ALIGN-TX_1.0.2.tar.gz\";\n$PG{\"dialign-tx\"}{\"mo\
de\"}=\"mcoffee\";\n$PG{\"dialign-tx\"}{\"binary\"\
}=\"dialign-tx\";\n$PG{\"dialign-tx\"}{\"version\"\
}=\"1.0.2\";\n$PG{\"poa\"}{\"4_TCOFFEE\"}=\"POA\";\
\n$PG{\"poa\"}{\"type\"}=\"sequence_multiple_align\
er\";\n$PG{\"poa\"}{\"ADDRESS\"}=\"http://www.bioi\
nformatics.ucla.edu/poa/\";\n$PG{\"poa\"}{\"langua\
ge\"}=\"C\";\n$PG{\"poa\"}{\"language2\"}=\"C\";\n\
$PG{\"poa\"}{\"source\"}=\"http://downloads.source\
forge.net/poamsa/poaV2.tar.gz\";\n$PG{\"poa\"}{\"D\
IR\"}=\"/usr/share/\";\n$PG{\"poa\"}{\"FILE1\"}=\"\
blosum80.mat\";\n$PG{\"poa\"}{\"mode\"}=\"mcoffee\\
";\n$PG{\"poa\"}{\"binary\"}=\"poa\";\n$PG{\"poa\"\
}{\"version\"}=\"2.0\";\n$PG{\"probcons\"}{\"4_TCO\
FFEE\"}=\"PROBCONS\";\n$PG{\"probcons\"}{\"type\"}\
=\"sequence_multiple_aligner\";\n$PG{\"probcons\"}\
{\"ADDRESS\"}=\"http://probcons.stanford.edu/\";\n\
$PG{\"probcons\"}{\"language2\"}=\"CXX\";\n$PG{\"p\
robcons\"}{\"language\"}=\"C++\";\n$PG{\"probcons\\
"}{\"source\"}=\"http://probcons.stanford.edu/prob\
cons_v1_12.tar.gz\";\n$PG{\"probcons\"}{\"mode\"}=\
\"mcoffee\";\n$PG{\"probcons\"}{\"binary\"}=\"prob\
cons\";\n$PG{\"probcons\"}{\"version\"}=\"1.12\";\\
n$PG{\"mafft\"}{\"4_TCOFFEE\"}=\"MAFFT\";\n$PG{\"m\
afft\"}{\"type\"}=\"sequence_multiple_aligner\";\n\
$PG{\"mafft\"}{\"ADDRESS\"}=\"http://align.bmr.kyu\
shu-u.ac.jp/mafft/online/server/\";\n$PG{\"mafft\"\
}{\"language\"}=\"C\";\n$PG{\"mafft\"}{\"language\\
"}=\"C\";\n$PG{\"mafft\"}{\"source\"}=\"http://ali\
gn.bmr.kyushu-u.ac.jp/mafft/software/mafft-6.603-w\
ith-extensions-src.tgz\";\n$PG{\"mafft\"}{\"window\
s\"}=\"http://align.bmr.kyushu-u.ac.jp/mafft/softw\
are/mafft-6.603-mingw.tar\";\n$PG{\"mafft\"}{\"mod\
e\"}=\"mcoffee,rcoffee\";\n$PG{\"mafft\"}{\"binary\
\"}=\"mafft.tar.gz\";\n$PG{\"mafft\"}{\"version\"}\
=\"6.603\";\n$PG{\"muscle\"}{\"4_TCOFFEE\"}=\"MUSC\
LE\";\n$PG{\"muscle\"}{\"type\"}=\"sequence_multip\
le_aligner\";\n$PG{\"muscle\"}{\"ADDRESS\"}=\"http\
://www.drive5.com/muscle/\";\n$PG{\"muscle\"}{\"la\
nguage\"}=\"C++\";\n$PG{\"muscle\"}{\"language2\"}\
=\"GPP\";\n$PG{\"muscle\"}{\"source\"}=\"http://ww\
w.drive5.com/muscle/downloads3.7/muscle3.7_src.tar\
.gz\";\n$PG{\"muscle\"}{\"windows\"}=\"http://www.\
drive5.com/muscle/downloads3.7/muscle3.7_win32.zip\
\";\n$PG{\"muscle\"}{\"linux\"}=\"http://www.drive\
5.com/muscle/downloads3.7/muscle3.7_linux_ia32.tar\
.gz\";\n$PG{\"muscle\"}{\"mode\"}=\"mcoffee,rcoffe\
e\";\n$PG{\"muscle\"}{\"version\"}=\"3.7\";\n$PG{\\
"mus4\"}{\"4_TCOFFEE\"}=\"MUS4\";\n$PG{\"mus4\"}{\\
"type\"}=\"sequence_multiple_aligner\";\n$PG{\"mus\
4\"}{\"ADDRESS\"}=\"http://www.drive5.com/muscle/\\
";\n$PG{\"mus4\"}{\"language\"}=\"C++\";\n$PG{\"mu\
s4\"}{\"language2\"}=\"GPP\";\n$PG{\"mus4\"}{\"sou\
rce\"}=\"http://www.drive5.com/muscle/muscle4.0_sr\
c.tar.gz\";\n$PG{\"mus4\"}{\"mode\"}=\"mcoffee,rco\
ffee\";\n$PG{\"mus4\"}{\"version\"}=\"4.0\";\n$PG{\
\"pcma\"}{\"4_TCOFFEE\"}=\"PCMA\";\n$PG{\"pcma\"}{\
\"type\"}=\"sequence_multiple_aligner\";\n$PG{\"pc\
ma\"}{\"ADDRESS\"}=\"ftp://iole.swmed.edu/pub/PCMA\
/\";\n$PG{\"pcma\"}{\"language\"}=\"C\";\n$PG{\"pc\
ma\"}{\"language2\"}=\"C\";\n$PG{\"pcma\"}{\"sourc\
e\"}=\"ftp://iole.swmed.edu/pub/PCMA/pcma.tar.gz\"\
;\n$PG{\"pcma\"}{\"mode\"}=\"mcoffee\";\n$PG{\"pcm\
a\"}{\"version\"}=\"1.0\";\n$PG{\"kalign\"}{\"4_TC\
OFFEE\"}=\"KALIGN\";\n$PG{\"kalign\"}{\"type\"}=\"\
sequence_multiple_aligner\";\n$PG{\"kalign\"}{\"AD\
DRESS\"}=\"http://msa.cgb.ki.se\";\n$PG{\"kalign\"\
}{\"language\"}=\"C\";\n$PG{\"kalign\"}{\"language\
2\"}=\"C\";\n$PG{\"kalign\"}{\"source\"}=\"http://\
msa.cgb.ki.se/downloads/kalign/current.tar.gz\";\n\
$PG{\"kalign\"}{\"mode\"}=\"mcoffee\";\n$PG{\"kali\
gn\"}{\"version\"}=\"1.0\";\n$PG{\"amap\"}{\"4_TCO\
FFEE\"}=\"AMAP\";\n$PG{\"amap\"}{\"type\"}=\"seque\
nce_multiple_aligner\";\n$PG{\"amap\"}{\"ADDRESS\"\
}=\"http://bio.math.berkeley.edu/amap/\";\n$PG{\"a\
map\"}{\"language\"}=\"C++\";\n$PG{\"amap\"}{\"lan\
guage2\"}=\"CXX\";\n$PG{\"amap\"}{\"source\"}=\"ht\
tp://amap-align.googlecode.com/files/amap.2.0.tar.\
gz\";\n$PG{\"amap\"}{\"mode\"}=\"mcoffee\";\n$PG{\\
"amap\"}{\"version\"}=\"2.0\";\n$PG{\"proda\"}{\"4\
_TCOFFEE\"}=\"PRODA\";\n$PG{\"proda\"}{\"type\"}=\\
"sequence_multiple_aligner\";\n$PG{\"proda\"}{\"AD\
DRESS\"}=\"http://proda.stanford.edu\";\n$PG{\"pro\
da\"}{\"language\"}=\"C++\";\n$PG{\"proda\"}{\"lan\
guage2\"}=\"CXX\";\n$PG{\"proda\"}{\"source\"}=\"h\
ttp://proda.stanford.edu/proda_1_0.tar.gz\";\n$PG{\
\"proda\"}{\"mode\"}=\"mcoffee\";\n$PG{\"proda\"}{\
\"version\"}=\"1.0\";\n$PG{\"fsa\"}{\"4_TCOFFEE\"}\
=\"FSA\";\n$PG{\"fsa\"}{\"type\"}=\"sequence_multi\
ple_aligner\";\n$PG{\"fsa\"}{\"ADDRESS\"}=\"http:/\
/fsa.sourceforge.net/\";\n$PG{\"fsa\"}{\"language\\
"}=\"C++\";\n$PG{\"fsa\"}{\"language2\"}=\"CXX\";\\
n$PG{\"fsa\"}{\"source\"}=\"http://sourceforge.net\
/projects/fsa/files/fsa-1.15.3.tar.gz/download/\";\
\n$PG{\"fsa\"}{\"mode\"}=\"mcoffee\";\n$PG{\"fsa\"\
}{\"version\"}=\"1.15.3\";\n$PG{\"prank\"}{\"4_TCO\
FFEE\"}=\"PRANK\";\n$PG{\"prank\"}{\"type\"}=\"seq\
uence_multiple_aligner\";\n$PG{\"prank\"}{\"ADDRES\
S\"}=\"http://www.ebi.ac.uk/goldman-srv/prank/\";\\
n$PG{\"prank\"}{\"language\"}=\"C++\";\n$PG{\"pran\
k\"}{\"language2\"}=\"CXX\";\n$PG{\"prank\"}{\"sou\
rce\"}=\"http://www.ebi.ac.uk/goldman-srv/prank/sr\
c/prank/prank.src.100303.tgz\";\n$PG{\"prank\"}{\"\
mode\"}=\"mcoffee\";\n$PG{\"prank\"}{\"version\"}=\
\"100303\";\n$PG{\"sap\"}{\"4_TCOFFEE\"}=\"SAP\";\\
n$PG{\"sap\"}{\"type\"}=\"structure_pairwise_align\
er\";\n$PG{\"sap\"}{\"ADDRESS\"}=\"http://mathbio.\
nimr.mrc.ac.uk/wiki/Software\";\n$PG{\"sap\"}{\"la\
nguage\"}=\"C\";\n$PG{\"sap\"}{\"language2\"}=\"C\\
";\n$PG{\"sap\"}{\"source\"}=\"http://mathbio.nimr\
.mrc.ac.uk/download/sap-1.1.1.tar.gz\";\n$PG{\"sap\
\"}{\"mode\"}=\"expresso,3dcoffee\";\n$PG{\"sap\"}\
{\"version\"}=\"1.1.1\";\n$PG{\"TMalign\"}{\"4_TCO\
FFEE\"}=\"TMALIGN\";\n$PG{\"TMalign\"}{\"type\"}=\\
"structure_pairwise_aligner\";\n$PG{\"TMalign\"}{\\
"ADDRESS\"}=\"http://zhang.bioinformatics.ku.edu/T\
M-align/TMalign.f\";\n$PG{\"TMalign\"}{\"language\\
"}=\"Fortran\";\n$PG{\"TMalign\"}{\"language2\"}=\\
"Fortran\";\n$PG{\"TMalign\"}{\"source\"}=\"http:/\
/zhang.bioinformatics.ku.edu/TM-align/TMalign.f\";\
\n$PG{\"TMalign\"}{\"linux\"}=\"http://zhang.bioin\
formatics.ku.edu/TM-align/TMalign_32.gz\";\n$PG{\"\
TMalign\"}{\"mode\"}=\"expresso,3dcoffee\";\n$PG{\\
"TMalign\"}{\"version\"}=\"1.0\";\n$PG{\"mustang\"\
}{\"4_TCOFFEE\"}=\"MUSTANG\";\n$PG{\"mustang\"}{\"\
type\"}=\"structure_pairwise_aligner\";\n$PG{\"mus\
tang\"}{\"ADDRESS\"}=\"http://www.cs.mu.oz.au/~aru\
n/mustang\";\n$PG{\"mustang\"}{\"language\"}=\"C++\
\";\n$PG{\"mustang\"}{\"language2\"}=\"CXX\";\n$PG\
{\"mustang\"}{\"source\"}=\"http://ww2.cs.mu.oz.au\
/~arun/mustang/mustang_v3.2.1.tgz\";\n$PG{\"mustan\
g\"}{\"mode\"}=\"expresso,3dcoffee\";\n$PG{\"musta\
ng\"}{\"version\"}=\"3.2.1\";\n$PG{\"lsqman\"}{\"4\
_TCOFFEE\"}=\"LSQMAN\";\n$PG{\"lsqman\"}{\"type\"}\
=\"structure_pairwise_aligner\";\n$PG{\"lsqman\"}{\
\"ADDRESS\"}=\"empty\";\n$PG{\"lsqman\"}{\"languag\
e\"}=\"empty\";\n$PG{\"lsqman\"}{\"language2\"}=\"\
empty\";\n$PG{\"lsqman\"}{\"source\"}=\"empty\";\n\
$PG{\"lsqman\"}{\"update_action\"}=\"never\";\n$PG\
{\"lsqman\"}{\"mode\"}=\"expresso,3dcoffee\";\n$PG\
{\"align_pdb\"}{\"4_TCOFFEE\"}=\"ALIGN_PDB\";\n$PG\
{\"align_pdb\"}{\"type\"}=\"structure_pairwise_ali\
gner\";\n$PG{\"align_pdb\"}{\"ADDRESS\"}=\"empty\"\
;\n$PG{\"align_pdb\"}{\"language\"}=\"empty\";\n$P\
G{\"align_pdb\"}{\"language2\"}=\"empty\";\n$PG{\"\
align_pdb\"}{\"source\"}=\"empty\";\n$PG{\"align_p\
db\"}{\"update_action\"}=\"never\";\n$PG{\"align_p\
db\"}{\"mode\"}=\"expresso,3dcoffee\";\n$PG{\"fugu\
eali\"}{\"4_TCOFFEE\"}=\"FUGUE\";\n$PG{\"fugueali\\
"}{\"type\"}=\"structure_pairwise_aligner\";\n$PG{\
\"fugueali\"}{\"ADDRESS\"}=\"http://www-cryst.bioc\
.cam.ac.uk/fugue/download.html\";\n$PG{\"fugueali\\
"}{\"language\"}=\"empty\";\n$PG{\"fugueali\"}{\"l\
anguage2\"}=\"empty\";\n$PG{\"fugueali\"}{\"source\
\"}=\"empty\";\n$PG{\"fugueali\"}{\"update_action\\
"}=\"never\";\n$PG{\"fugueali\"}{\"mode\"}=\"expre\
sso,3dcoffee\";\n$PG{\"dalilite.pl\"}{\"4_TCOFFEE\\
"}=\"DALILITEc\";\n$PG{\"dalilite.pl\"}{\"type\"}=\
\"structure_pairwise_aligner\";\n$PG{\"dalilite.pl\
\"}{\"ADDRESS\"}=\"built_in\";\n$PG{\"dalilite.pl\\
"}{\"ADDRESS2\"}=\"http://www.ebi.ac.uk/Tools/webs\
ervices/services/dalilite\";\n$PG{\"dalilite.pl\"}\
{\"language\"}=\"Perl\";\n$PG{\"dalilite.pl\"}{\"l\
anguage2\"}=\"Perl\";\n$PG{\"dalilite.pl\"}{\"sour\
ce\"}=\"empty\";\n$PG{\"dalilite.pl\"}{\"update_ac\
tion\"}=\"never\";\n$PG{\"dalilite.pl\"}{\"mode\"}\
=\"expresso,3dcoffee\";\n$PG{\"probconsRNA\"}{\"4_\
TCOFFEE\"}=\"PROBCONSRNA\";\n$PG{\"probconsRNA\"}{\
\"type\"}=\"RNA_multiple_aligner\";\n$PG{\"probcon\
sRNA\"}{\"ADDRESS\"}=\"http://probcons.stanford.ed\
u/\";\n$PG{\"probconsRNA\"}{\"language\"}=\"C++\";\
\n$PG{\"probconsRNA\"}{\"language2\"}=\"CXX\";\n$P\
G{\"probconsRNA\"}{\"source\"}=\"http://probcons.s\
tanford.edu/probconsRNA.tar.gz\";\n$PG{\"probconsR\
NA\"}{\"mode\"}=\"mcoffee,rcoffee\";\n$PG{\"probco\
nsRNA\"}{\"version\"}=\"1.0\";\n$PG{\"sfold\"}{\"4\
_TCOFFEE\"}=\"CONSAN\";\n$PG{\"sfold\"}{\"type\"}=\
\"RNA_pairwise_aligner\";\n$PG{\"sfold\"}{\"ADDRES\
S\"}=\"http://selab.janelia.org/software/consan/\"\
;\n$PG{\"sfold\"}{\"language\"}=\"empty\";\n$PG{\"\
sfold\"}{\"language2\"}=\"empty\";\n$PG{\"sfold\"}\
{\"source\"}=\"empty\";\n$PG{\"sfold\"}{\"update_a\
ction\"}=\"never\";\n$PG{\"sfold\"}{\"mode\"}=\"rc\
offee\";\n$PG{\"RNAplfold\"}{\"4_TCOFFEE\"}=\"RNAP\
LFOLD\";\n$PG{\"RNAplfold\"}{\"type\"}=\"RNA_secon\
darystructure_predictor\";\n$PG{\"RNAplfold\"}{\"A\
DDRESS\"}=\"http://www.tbi.univie.ac.at/~ivo/RNA/\\
";\n$PG{\"RNAplfold\"}{\"language\"}=\"C\";\n$PG{\\
"RNAplfold\"}{\"language2\"}=\"C\";\n$PG{\"RNAplfo\
ld\"}{\"source\"}=\"http://www.tbi.univie.ac.at/~i\
vo/RNA/ViennaRNA-1.7.2.tar.gz\";\n$PG{\"RNAplfold\\
"}{\"mode\"}=\"rcoffee,\";\n$PG{\"RNAplfold\"}{\"v\
ersion\"}=\"1.7.2\";\n$PG{\"retree\"}{\"4_TCOFFEE\\
"}=\"PHYLIP\";\n$PG{\"retree\"}{\"type\"}=\"RNA_se\
condarystructure_predictor\";\n$PG{\"retree\"}{\"A\
DDRESS\"}=\"http://evolution.gs.washington.edu/phy\
lip/\";\n$PG{\"retree\"}{\"language\"}=\"C\";\n$PG\
{\"retree\"}{\"language2\"}=\"C\";\n$PG{\"retree\"\
}{\"source\"}=\"http://evolution.gs.washington.edu\
/phylip/download/phylip-3.69.tar.gz\";\n$PG{\"retr\
ee\"}{\"mode\"}=\"trmsd,\";\n$PG{\"retree\"}{\"ver\
sion\"}=\"3.69\";\n$PG{\"hmmtop\"}{\"4_TCOFFEE\"}=\
\"HMMTOP\";\n$PG{\"hmmtop\"}{\"type\"}=\"protein_s\
econdarystructure_predictor\";\n$PG{\"hmmtop\"}{\"\
ADDRESS\"}=\"www.enzim.hu/hmmtop/\";\n$PG{\"hmmtop\
\"}{\"language\"}=\"C\";\n$PG{\"hmmtop\"}{\"langua\
ge2\"}=\"C\";\n$PG{\"hmmtop\"}{\"source\"}=\"empty\
\";\n$PG{\"hmmtop\"}{\"update_action\"}=\"never\";\
\n$PG{\"hmmtop\"}{\"mode\"}=\"tcoffee\";\n$PG{\"go\
rIV\"}{\"4_TCOFFEE\"}=\"GOR4\";\n$PG{\"gorIV\"}{\"\
type\"}=\"protein_secondarystructure_predictor\";\\
n$PG{\"gorIV\"}{\"ADDRESS\"}=\"http://mig.jouy.inr\
a.fr/logiciels/gorIV/\";\n$PG{\"gorIV\"}{\"languag\
e\"}=\"C\";\n$PG{\"gorIV\"}{\"language2\"}=\"C\";\\
n$PG{\"gorIV\"}{\"source\"}=\"http://mig.jouy.inra\
.fr/logiciels/gorIV/GOR_IV.tar.gz\";\n$PG{\"gorIV\\
"}{\"update_action\"}=\"never\";\n$PG{\"gorIV\"}{\\
"mode\"}=\"tcoffee\";\n$PG{\"wublast.pl\"}{\"4_TCO\
FFEE\"}=\"EBIWUBLASTc\";\n$PG{\"wublast.pl\"}{\"ty\
pe\"}=\"protein_homology_predictor\";\n$PG{\"wubla\
st.pl\"}{\"ADDRESS\"}=\"built_in\";\n$PG{\"wublast\
.pl\"}{\"ADDRESS2\"}=\"http://www.ebi.ac.uk/Tools/\
webservices/services/wublast\";\n$PG{\"wublast.pl\\
"}{\"language\"}=\"Perl\";\n$PG{\"wublast.pl\"}{\"\
language2\"}=\"Perl\";\n$PG{\"wublast.pl\"}{\"sour\
ce\"}=\"empty\";\n$PG{\"wublast.pl\"}{\"update_act\
ion\"}=\"never\";\n$PG{\"wublast.pl\"}{\"mode\"}=\\
"psicoffee,expresso,accurate\";\n$PG{\"blastpgp.pl\
\"}{\"4_TCOFFEE\"}=\"EBIBLASTPGPc\";\n$PG{\"blastp\
gp.pl\"}{\"type\"}=\"protein_homology_predictor\";\
\n$PG{\"blastpgp.pl\"}{\"ADDRESS\"}=\"built_in\";\\
n$PG{\"blastpgp.pl\"}{\"ADDRESS2\"}=\"http://www.e\
bi.ac.uk/Tools/webservices/services/blastpgp\";\n$\
PG{\"blastpgp.pl\"}{\"language\"}=\"Perl\";\n$PG{\\
"blastpgp.pl\"}{\"language2\"}=\"Perl\";\n$PG{\"bl\
astpgp.pl\"}{\"source\"}=\"empty\";\n$PG{\"blastpg\
p.pl\"}{\"update_action\"}=\"never\";\n$PG{\"blast\
pgp.pl\"}{\"mode\"}=\"psicoffee,expresso,accurate\\
";\n$PG{\"blastcl3\"}{\"4_TCOFFEE\"}=\"NCBIWEBBLAS\
T\";\n$PG{\"blastcl3\"}{\"type\"}=\"protein_homolo\
gy_predictor\";\n$PG{\"blastcl3\"}{\"ADDRESS\"}=\"\
ftp://ftp.ncbi.nih.gov/blast/executables/LATEST\";\
\n$PG{\"blastcl3\"}{\"language\"}=\"C\";\n$PG{\"bl\
astcl3\"}{\"language2\"}=\"C\";\n$PG{\"blastcl3\"}\
{\"source\"}=\"empty\";\n$PG{\"blastcl3\"}{\"updat\
e_action\"}=\"never\";\n$PG{\"blastcl3\"}{\"mode\"\
}=\"psicoffee,expresso,3dcoffee\";\n$PG{\"blastpgp\
\"}{\"4_TCOFFEE\"}=\"NCBIBLAST\";\n$PG{\"blastpgp\\
"}{\"type\"}=\"protein_homology_predictor\";\n$PG{\
\"blastpgp\"}{\"ADDRESS\"}=\"ftp://ftp.ncbi.nih.go\
v/blast/executables/LATEST\";\n$PG{\"blastpgp\"}{\\
"language\"}=\"C\";\n$PG{\"blastpgp\"}{\"language2\
\"}=\"C\";\n$PG{\"blastpgp\"}{\"source\"}=\"empty\\
";\n$PG{\"blastpgp\"}{\"update_action\"}=\"never\"\
;\n$PG{\"blastpgp\"}{\"mode\"}=\"psicoffee,express\
o,3dcoffee\";\n$PG{\"SOAP::Lite\"}{\"4_TCOFFEE\"}=\
\"SOAPLITE\";\n$PG{\"SOAP::Lite\"}{\"type\"}=\"lib\
rary\";\n$PG{\"SOAP::Lite\"}{\"ADDRESS\"}=\"http:/\
/cpansearch.perl.org/src/MKUTTER/SOAP-Lite-0.710.0\
8/Makefile.PL\";\n$PG{\"SOAP::Lite\"}{\"language\"\
}=\"Perl\";\n$PG{\"SOAP::Lite\"}{\"language2\"}=\"\
Perl\";\n$PG{\"SOAP::Lite\"}{\"source\"}=\"empty\"\
;\n$PG{\"blastpgp\"}{\"update_action\"}=\"never\";\
\n$PG{\"SOAP::Lite\"}{\"mode\"}=\"none\";\n$PG{\"X\
ML::Simple\"}{\"4_TCOFFEE\"}=\"XMLSIMPLE\";\n$PG{\\
"XML::Simple\"}{\"type\"}=\"library\";\n$PG{\"XML:\
:Simple\"}{\"ADDRESS\"}=\"http://search.cpan.org/~\
grantm/XML-Simple-2.18/lib/XML/Simple.pm\";\n$PG{\\
"XML::Simple\"}{\"language\"}=\"Perl\";\n$PG{\"XML\
::Simple\"}{\"language2\"}=\"Perl\";\n$PG{\"XML::S\
imple\"}{\"source\"}=\"empty\";\n$PG{\"XML::Simple\
\"}{\"mode\"}=\"psicoffee,expresso,accurate\";\n$M\
ODE{\"tcoffee\"}{\"name\"}=\"tcoffee\";\n$MODE{\"r\
coffee\"}{\"name\"}=\"rcoffee\";\n$MODE{\"3dcoffee\
\"}{\"name\"}=\"3dcoffee\";\n$MODE{\"mcoffee\"}{\"\
name\"}=\"mcoffee\";\n$MODE{\"expresso\"}{\"name\"\
}=\"expresso\";\n$MODE{\"trmsd\"}{\"name\"}=\"trms\
d\";\n$MODE{\"accurate\"}{\"name\"}=\"accurate\";\\
n$MODE{\"seq_reformat\"}{\"name\"}=\"seq_reformat\\
";\n\n\n$PG{C}{compiler}=\"gcc\";\n$PG{C}{compiler\
_flag}=\"CC\";\n$PG{C}{options}=\"\";\n$PG{C}{opti\
ons_flag}=\"CFLAGS\";\n$PG{C}{type}=\"compiler\";\\
n\n$PG{\"CXX\"}{compiler}=\"g++\";\n$PG{\"CXX\"}{c\
ompiler_flag}=\"CXX\";\n$PG{\"CXX\"}{options}=\"\"\
;\n$PG{\"CXX\"}{options_flag}=\"CXXFLAGS\";\n$PG{C\
XX}{type}=\"compiler\";\n\n$PG{\"CPP\"}{compiler}=\
\"g++\";\n$PG{\"CPP\"}{compiler_flag}=\"CPP\";\n$P\
G{\"CPP\"}{options}=\"\";\n$PG{\"CPP\"}{options_fl\
ag}=\"CPPFLAGS\";\n$PG{CPP}{type}=\"compiler\";\n\\
n$PG{\"GPP\"}{compiler}=\"g++\";\n$PG{\"GPP\"}{com\
piler_flag}=\"GPP\";\n$PG{\"GPP\"}{options}=\"\";\\
n$PG{\"GPP\"}{options_flag}=\"CFLAGS\";\n$PG{GPP}{\
type}=\"compiler\";\n\n$PG{Fortran}{compiler}=\"g7\
7\";\n$PG{Fortran}{compiler_flag}=\"FCC\";\n$PG{Fo\
rtran}{type}=\"compiler\";\n\n$PG{Perl}{compiler}=\
\"CPAN\";\n$PG{Perl}{type}=\"compiler\";\n\n$SUPPO\
RTED_OS{macox}=\"Macintosh\";\n$SUPPORTED_OS{linux\
}=\"Linux\";\n$SUPPORTED_OS{windows}=\"Cygwin\";\n\
\n\n\n$MODE{t_coffee}{description}=\" for regular \
multiple sequence alignments\";\n$MODE{rcoffee} {d\
escription}=\" for RNA multiple sequence alignment\
s\";\n\n$MODE{psicoffee} {description}=\" for Homo\
logy Extended multiple sequence alignments\";\n$MO\
DE{expresso}{description}=\" for very accurate str\
ucture based multiple sequence alignments\";\n$MOD\
E{\"3dcoffee\"}{description}=\" for multiple struc\
ture alignments\";\n$MODE{mcoffee} {description}=\\
" for combining alternative multiple sequence alig\
nment packages\\n------- into a unique meta-packag\
e. The installer will upload several MSA packages \
and compile them\\n\n\";\n\n\n&post_process_PG();\\
nreturn;\n}\n\nsub post_process_PG\n  {\n    my $p\
;\n    \n    %PG=&name2dname (%PG);\n    %MODE=&na\
me2dname(%MODE);\n    foreach $p (keys(%PG)){if ( \
$PG{$p}{type} eq \"compiler\"){$PG{$p}{update_acti\
on}=\"never\";}}\n    \n  }\n\nsub name2dname\n  {\
\n    my (%L)=(@_);\n    my ($l, $ml);\n    \n    \
foreach my $pg (keys(%L))\n      {\n	$l=length ($p\
g);\n	if ( $l>$ml){$ml=$l;}\n      }\n    $ml+=1;\\
n    foreach my $pg (keys(%L))\n      {\n	my $name\
;\n	$l=$ml-length ($pg);\n	$name=$pg;\n	for ( $b=0\
; $b<$l; $b++)\n	  {\n	    $name .=\" \";\n	  }\n	\
$L{$pg}{dname}=$name;\n      }\n    return %L;\n  \
}\n\nsub env_file2putenv\n  {\n    my $f=@_[0];\n \
   my $F=new FileHandle;\n    my $n;\n    \n    op\
en ($F, \"$f\");\n    while (<$F>)\n      {\n	my $\
line=$_;\n	my($var, $value)=($_=~/(\\S+)\\=(\\S*)/\
);\n	$ENV{$var}=$value;\n	$ENV_SET{$var}=1;\n	$n++\
;\n      }\n    close ($F);\n    return $n;\n  }\n\
\nsub replace_line_in_file\n  {\n    my ($file, $w\
ordin, $wordout)=@_;\n    my $O=new FileHandle;\n \
   my $I=new FileHandle;\n    my $l;\n    if (!-e \
$file){return;}\n    \n    system (\"mv $file $fil\
e.old\");\n    open ($O, \">$file\");\n    open ($\
I, \"$file.old\");\n    while (<$I>)\n      {\n	$l\
=$_;\n	if (!($l=~/$wordin/)){print $O \"$l\";}\n	e\
lsif ( $wordout ne \"\"){$l=~s/$wordin/$wordout/g;\
print $O \"$l\";}\n      }\n    close ($O);\n    c\
lose ($I);\n    return;\n  }\n\nsub add_C_librarie\
s\n  {\n   my ($file,$first,@list)=@_;\n   \n    m\
y $O=new FileHandle;\n    my $I=new FileHandle;\n \
   my ($l,$anchor);\n    if (!-e $file){return;}\n\
   \n    $anchor=\"#include <$first>\";\n	 \n    s\
ystem (\"mv $file $file.old\");\n    open ($O, \">\
$file\");\n    open ($I, \"$file.old\");\n    whil\
e (<$I>)\n      {\n	$l=$_;\n	print $O \"$l\";\n	if\
 (!($l=~/$anchor/))\n	   {\n	    \n	    foreach my\
 $lib (@list)\n	       {\n                  print \
$O \"#include <$lib>\\n\";\n	       }\n           \
}\n      }\n    close ($O);\n    close ($I);\n    \
return;\n    }\n","use Env;\nuse Cwd;\n@suffix=(\"\
tmp\", \"temp\", \"cache\", \"t_coffee\", \"core\"\
, \"tcoffee\");\n\nif ($#ARGV==-1)\n  {\n    print\
 \"clean_cache.pl -file <file to add in -dir> -dir\
=<dir> -size=<value in Mb>\\n0: unlimited -1 alway\
s.\\nWill only clean directories matching:[\";\n  \
  foreach $k(@suffix){print \"*$k* \";}\n    print\
 \"]\\n\";\n    exit (EXIT_FAILURE);\n  }\n\n$cl=j\
oin (\" \",@ARGV);\nif (($cl=~/\\-no_action/))\n  \
{\n    exit (EXIT_SUCCESS);\n  }\n\nif (($cl=~/\\-\
debug/))\n  {\n    $DEBUG=1;\n  }\nelse\n  {\n    \
$DEBUG=0;\n  }\n\nif (($cl=~/\\-dir=(\\S+)/))\n  {\
\n    $dir=$1;\n  }\nelse\n  {\n    $dir=\"./\";\n\
  }\n\nif ($cl=~/\\-file=(\\S+)/)\n  {\n    $file=\
$1;\n  }\nelse\n  {\n    $file=0;\n  }\n\nif ($cl=\
~/\\-size=(\\S+)/)\n  {\n    $max_size=$1;\n  }\ne\
lse\n  {\n    $max_size=0;#unlimited\n  }\nif ($cl\
=~/\\-force/)\n  {\n    $force=1;\n  }\nelse\n  {\\
n    $force=0;\n  }\n\nif ($cl=~/\\-age=(\\S+)/)\n\
  {\n    $max_age=$1;\n  }\nelse\n  {\n    $max_ag\
e=0;#unlimited\n  }\n\n$max_size*=1000000;\nif ( !\
 -d $dir)\n  {\n    print STDERR \"\\nCannot proce\
ss $dir: does not exist \\n\";\n    exit (EXIT_FAI\
LURE);\n  }\n\nif ( !($dir=~/^\\//))\n  {\n    $ba\
se=cwd();\n    $dir=\"$base/$dir\";\n  }\n\n$proce\
ed=0;\nforeach $s (@suffix)\n  {\n    \n    if (($\
dir=~/$s/)){$proceed=1;}\n    $s=uc ($s);\n    if \
(($dir=~/$s/)){$proceed=1;}\n  }\nif ( $proceed==0\
)\n  {\n    print STDERR \"Clean_cache.pl can only\
 clean directories whose absolute path name contai\
ns the following strings:\";\n    foreach $w (@suf\
fix) {print STDERR \"$w \";$w=lc($w); print STDERR\
 \"$w \";}\n    print STDERR \"\\nCannot process $\
dir\\n\";\n    exit (EXIT_FAILURE);\n  }\n\n$name_\
file=\"$dir/name_file.txt\";\n$size_file=\"$dir/si\
ze_file.txt\";\nif ( $force){&create_ref_file ($di\
r,$name_file,$size_file);}\nif ($file){&add_file (\
$dir, $name_file, $size_file, $file);}\n&clean_dir\
 ($dir, $name_file, $size_file, $max_size,$max_age\
);\nexit (EXIT_SUCCESS);\n\nsub clean_dir \n  {\n \
   my ($dir, $name_file, $size_file, $max_size, $m\
ax_age)=@_;\n    my ($tot_size, $size, $f, $s);\n\\
n  \n    $tot_size=&get_tot_size ($dir, $name_file\
, $size_file);\n\n    if ( $tot_size<=$max_size){r\
eturn ;}\n    else {$max_size/=2;}\n    \n    #rec\
reate the name file in case some temprary files ha\
ve not been properly registered\n    &create_ref_f\
ile ($dir, $name_file, $size_file, $max_age);\n  \\
n    $new_name_file=&vtmpnam();\n    open (R, \"$n\
ame_file\");\n    open (W, \">$new_name_file\");\n\
    while (<R>)\n      {\n	my $line=$_;\n	\n	($f, \
$s)=($line=~/(\\S+) (\\S+)/);\n	if ( !($f=~/\\S/))\
{next;}\n	\n	elsif ($max_size && $tot_size>=$max_s\
ize && !($f=~/name_file/))\n	  {\n	    remove ( \"\
$dir/$f\");\n	    $tot_size-=$s;\n	  }\n	elsif ( $\
max_age && -M(\"$dir/$f\")>=$max_age)\n	  {\n	    \
remove ( \"$dir/$f\");\n	    $tot_size-=$s;\n	  }\\
n	else\n	  {\n	    print W \"$f $s\\n\";\n	  }\n  \
    }\n    close (R);\n    close (W);\n    open (F\
, \">$size_file\");\n    print F \"$tot_size\";\n \
   if ( -e $new_name_file){`mv $new_name_file $nam\
e_file`;}\n    close (F);\n  }\nsub get_tot_size\n\
  {\n    my ($dir, $name_file, $size_file)=@_;\n  \
  my $size;\n    \n    if ( !-d $dir){return 0;}\n\
    if ( !-e $name_file)\n      {\n	\n	&create_ref\
_file ($dir, $name_file, $size_file);\n      }\n  \
  open (F, \"$size_file\");\n    $size=<F>;\n    c\
lose (F);\n    chomp ($size);\n    return $size;\n\
  }\nsub size \n  {\n    my $f=@_[0];\n\n    if ( \
!-d $f){return -s($f);}\n    else {return &dir2siz\
e($f);}\n  }\nsub dir2size\n  {\n    my $d=@_[0];\\
n    my ($s, $f);\n    \n    if ( !-d $d) {return \
0;}\n    \n    foreach $f (&dir2list ($d))\n      \
{\n	if ( -d $f){$s+=&dir2size (\"$d/$f\");}\n	else\
 {$s+= -s \"$dir/$f\";}\n      }\n    return $s;\n\
  }\n\nsub remove \n  {\n    my $file=@_[0];\n    \
my ($f);\n    \n    debug_print( \"--- $file ---\\\
n\");\n    if (($file eq \".\") || ($file eq \"..\\
") || ($file=~/\\*/)){return EXIT_FAILURE;}\n    e\
lsif ( !-d $file)\n      {\n	debug_print (\"unlink\
 $file\\n\");\n	if (-e $file){unlink ($file);}\n  \
    }\n    elsif ( -d $file)\n      {\n	debug_prin\
t (\"++++++++ $file +++++++\\n\");\n	foreach $f (&\
dir2list($file))\n	  {\n	    &remove (\"$file/$f\"\
);\n	  }\n	debug_print (\"rmdir $file\\n\");\n	rmd\
ir $file;\n      }\n    else\n      {\n	debug_prin\
t (\"????????? $file ????????\\n\");\n      }\n   \
 return EXIT_SUCCESS;\n  }\n\nsub dir2list\n  {\n \
   my $dir=@_[0];\n    my (@list1, @list2,@list3, \
$l);\n\n    opendir (DIR,$dir);\n    @list1=readdi\
r (DIR);\n    closedir (DIR);\n    \n    foreach $\
l (@list1)\n      {\n	if ( $l ne \".\" && $l ne \"\
..\"){@list2=(@list2, $l);}\n      }\n    @list3 =\
 sort { (-M \"$dir/$list2[$b]\") <=> (-M \"$dir/$l\
ist2[$a]\")} @list2;\n    return @list3;\n    \n  \
}\n\nsub debug_print\n  {\n    \n    if ($DEBUG==1\
){print @_;}\n    \n  }\nsub create_ref_file\n  {\\
n    my ($dir,$name_file,$size_file)=@_;\n    my (\
$f, $s, $tot_size, @l);\n    \n    if ( !-d $dir){\
return;}\n    \n    @l=&dir2list ($dir);\n    open\
 (F, \">$name_file\");\n    foreach $f (@l)\n     \
 {\n	$s=&size(\"$dir/$f\");\n	$tot_size+=$s;\n	pri\
nt F \"$f $s\\n\";\n      }\n    &myecho ($tot_siz\
e, \">$size_file\");\n    close (F);\n  }\nsub add\
_file \n  {\n    my ($dir,$name_file,$size_file,$f\
ile)=@_;\n    my ($s, $tot_size);\n    \n    if ( \
!-d $dir)   {return;}\n    if ( !-e \"$dir/$file\"\
 ) {return;}\n    if ( !-e $name_file){&create_ref\
_file ($dir,$name_file,$size_file);}\n					    \n \
   $s=&size(\"$dir/$file\");\n    open (F, \">>$na\
me_file\");\n    print F \"$file\\n\";\n    close \
(F);\n\n    $tot_size=&get_tot_size ($dir,$name_fi\
le,$size_file);\n    $tot_size+=$s;\n    &myecho (\
$tot_size, \">$size_file\");\n    \n  }\n	\nsub my\
echo\n  {\n    my ($string, $file)=@_;\n    open (\
ECHO, $file) || die;\n    print ECHO \"$string\";\\
n    close (ECHO);\n  }\n    \n		\n	\nsub vtmpnam\\
n  {\n    my $tmp_file_name;\n    $tmp_name_counte\
r++;\n    $tmp_file_name=\"tmp_file_for_clean_cach\
e_pdb$$.$tmp_name_counter\";\n    $tmp_file_list[$\
ntmp_file++]=$tmp_file_name;\n    if ( -e $tmp_fil\
e_name) {return &vtmpnam ();}\n    else {return $t\
mp_file_name;}\n  }\n","\n$t_coffee=\"t_coffee\";\\
n\nforeach $value ( @ARGV)\n  {\n    $seq_file=$se\
q_file.\" \".$value;\n  }\n\n$name=$ARGV[0];\n$nam\
e=~s/\\.[^\\.]*$//;\n$lib_name=\"$name.mocca_lib\"\
;\n$type=`t_coffee $seq_file -get_type -quiet`;\nc\
hop ($type);\n\nif ( $type eq \"PROTEIN\"){$lib_mo\
de=\"lalign_rs_s_pair -lalign_n_top 20\";}\nelsif \
( $type eq\"DNA\"){$lib_mode=\"lalign_rs_s_dna_pai\
r -lalign_n_top 40\";}\n\nif ( !(-e $lib_name))\n \
 {\n	  \n  $command=\"$t_coffee -mocca -seq_weight\
=no -cosmetic_penalty=0 -mocca_interactive -in $li\
b_mode -out_lib $lib_name -infile $seq_file\";\n  \
\n  }\nelsif ( (-e $lib_name))\n  {\n  $command=\"\
$t_coffee -mocca -seq_weight=no -cosmetic_penalty=\
0 -mocca_interactive -in $lib_name -infile $seq_fi\
le\";\n  \n  }\n\nsystem ($command);\n\nexit;\n\n"\
,"my $WSDL = 'http://www.ebi.ac.uk/Tools/webservic\
es/wsdl/WSDaliLite.wsdl';\n\nuse SOAP::Lite;\nuse \
Data::Dumper;\nuse Getopt::Long qw(:config no_igno\
re_case bundling);\nuse File::Basename;\n\nmy $che\
ckInterval = 5;\n\nmy %params=(\n	    'async' => '\
1', # Use async mode and simulate sync mode in cli\
ent\n	    );\nGetOptions(\n    'pdb1=s'     => \\$\
params{'sequence1'},\n    'chainid1=s' => \\$param\
s{'chainid1'},\n    'pdb2=s'     => \\$params{'seq\
uence2'},\n    'chainid2=s' => \\$params{'chainid2\
'},\n    \"help|h\"	 => \\$help, # Usage info\n   \
 \"async|a\"	 => \\$async, # Asynchronous submissi\
on\n    \"polljob\"	 => \\$polljob, # Get results\\
n    \"status\"	 => \\$status, # Get status\n    \\
"jobid|j=s\"  => \\$jobid, # JobId\n    \"email|S=\
s\"  => \\$params{email}, # E-mail address\n    \"\
trace\"      => \\$trace, # SOAP messages\n    \"s\
equence=s\" => \\$sequence, # Input PDB\n    );\n\\
nmy $scriptName = basename($0, ());\nif($help) {\n\
    &usage();\n    exit(0);\n}\n\nif($trace) {\n  \
  print \"Tracing active\\n\";\n    SOAP::Lite->im\
port(+trace => 'debug');\n}\n\nmy $soap = SOAP::Li\
te\n    ->service($WSDL)\n    ->on_fault(sub {\n  \
      my $soap = shift;\n        my $res = shift;\\
n        # Throw an exception for all faults\n    \
    if(ref($res) eq '') {\n            die($res);\\
n        } else {\n            die($res->faultstri\
ng);\n        }\n        return new SOAP::SOM;\n  \
  }\n               );\n\nif( !($polljob || $statu\
s) &&\n    !( defined($params{'sequence1'}) && def\
ined($params{'sequence2'}) )\n    ) {\n    print S\
TDERR 'Error: bad option combination', \"\\n\";\n \
   &usage();\n    exit(1);\n}\nelsif($polljob && d\
efined($jobid)) {\n    print \"Getting results for\
 job $jobid\\n\";\n    getResults($jobid);\n}\nels\
if($status && defined($jobid)) {\n    print STDERR\
 \"Getting status for job $jobid\\n\";\n    my $re\
sult = $soap->checkStatus($jobid);\n    print STDO\
UT \"$result\", \"\\n\";\n    if($result eq 'DONE'\
) {\n	print STDERR \"To get results: $scriptName -\
-polljob --jobid $jobid\\n\";\n    }\n}\nelse {\n \
   if(-f $params{'sequence1'}) {\n	$params{'sequen\
ce1'} = read_file($params{'sequence1'});\n    }\n \
   if(-f $params{'sequence2'}) {\n	$params{'sequen\
ce2'} = read_file($params{'sequence2'});\n    }\n\\
n    my $jobid;\n    my $paramsData = SOAP::Data->\
name('params')->type(map=>\\%params);\n    # For S\
OAP::Lite 0.60 and earlier parameters are passed d\
irectly\n    if($SOAP::Lite::VERSION eq '0.60' || \
$SOAP::Lite::VERSION =~ /0\\.[1-5]/) {\n        $j\
obid = $soap->runDaliLite($paramsData);\n    }\n  \
  # For SOAP::Lite 0.69 and later parameter handli\
ng is different, so pass\n    # undef's for templa\
ted params, and then pass the formatted args.\n   \
 else {\n        $jobid = $soap->runDaliLite(undef\
,\n				     $paramsData);\n    }\n\n    if (define\
d($async)) {\n	print STDOUT $jobid, \"\\n\";\n    \
    print STDERR \"To check status: $scriptName --\
status --jobid $jobid\\n\";\n    } else { # Synchr\
onous mode\n        print STDERR \"JobId: $jobid\\\
n\";\n        sleep 1;\n        getResults($jobid)\
;\n    }\n}\n\nsub clientPoll($) {\n    my $jobid \
= shift;\n    my $result = 'PENDING';\n    # Check\
 status and wait if not finished\n    #print STDER\
R \"Checking status: $jobid\\n\";\n    while($resu\
lt eq 'RUNNING' || $result eq 'PENDING') {\n      \
  $result = $soap->checkStatus($jobid);\n        p\
rint STDERR \"$result\\n\";\n        if($result eq\
 'RUNNING' || $result eq 'PENDING') {\n           \
 # Wait before polling again.\n            sleep $\
checkInterval;\n        }\n    }\n}\n\nsub getResu\
lts($) {\n    $jobid = shift;\n    # Check status,\
 and wait if not finished\n    clientPoll($jobid);\
\n    # Use JobId if output file name is not defin\
ed\n    unless(defined($outfile)) {\n        $outf\
ile=$jobid;\n    }\n    # Get list of data types\n\
    my $resultTypes = $soap->getResults($jobid);\n\
    # Get the data and write it to a file\n    if(\
defined($outformat)) { # Specified data type\n    \
    my $selResultType;\n        foreach my $result\
Type (@$resultTypes) {\n            if($resultType\
->{type} eq $outformat) {\n                $selRes\
ultType = $resultType;\n            }\n        }\n\
        $res=$soap->poll($jobid, $selResultType->{\
type});\n        write_file($outfile.'.'.$selResul\
tType->{ext}, $res);\n    } else { # Data types av\
ailable\n        # Write a file for each output ty\
pe\n        for my $resultType (@$resultTypes){\n \
           #print \"Getting $resultType->{type}\\n\
\";\n            $res=$soap->poll($jobid, $resultT\
ype->{type});\n            write_file($outfile.'.'\
.$resultType->{ext}, $res);\n        }\n    }\n}\n\
\nsub read_file($) {\n    my $filename = shift;\n \
   open(FILE, $filename);\n    my $content;\n    m\
y $buffer;\n    while(sysread(FILE, $buffer, 1024)\
) {\n	$content.= $buffer;\n    }\n    close(FILE);\
\n    return $content;\n}\n\nsub write_file($$) {\\
n    my ($tmp,$entity) = @_;\n    print STDERR \"C\
reating result file: \".$tmp.\"\\n\";\n    unless(\
open (FILE, \">$tmp\")) {\n	return 0;\n    }\n    \
syswrite(FILE, $entity);\n    close (FILE);\n    r\
eturn 1;\n}\n\nsub usage {\n    print STDERR <<EOF\
\nDaliLite\n========\n\nPairwise comparison of pro\
tein structures\n\n[Required]\n\n  --pdb1         \
       : str  : PDB ID for structure 1\n  --pdb2  \
              : str  : PDB ID for structure 2\n\n[\
Optional]\n\n  --chain1              : str  : Chai\
n identifer in structure 1\n  --chain2            \
  : str  : Chain identifer in structure 2\n\n[Gene\
ral]\n\n  -h, --help            :      : prints th\
is help text\n  -S, --email           : str  : use\
r email address\n  -a, --async           :      : \
asynchronous submission\n      --status          :\
      : poll for the status of a job\n      --poll\
job         :      : poll for the results of a job\
\n  -j, --jobid           : str  : jobid for an as\
ynchronous job\n  -O, --outfile         : str  : f\
ile name for results (default is jobid)\n      --t\
race	        :      : show SOAP messages being int\
erchanged \n\nSynchronous job:\n\n  The results/er\
rors are returned as soon as the job is finished.\\
n  Usage: $scriptName --email <your\\@email> [opti\
ons] pdbFile [--outfile string]\n  Returns: saves \
the results to disk\n\nAsynchronous job:\n\n  Use \
this if you want to retrieve the results at a late\
r time. The results \n  are stored for up to 24 ho\
urs. \n  The asynchronous submission mode is recom\
mended when users are submitting \n  batch jobs or\
 large database searches	\n  Usage: $scriptName --\
email <your\\@email> --async [options] pdbFile\n  \
Returns: jobid\n\n  Use the jobid to query for the\
 status of the job. \n  Usage: $scriptName --statu\
s --jobid <jobId>\n  Returns: string indicating th\
e status of the job:\n    DONE - job has finished\\
n    RUNNING - job is running\n    NOT_FOUND - job\
 cannot be found\n    ERROR - the jobs has encount\
ered an error\n\n  When done, use the jobid to ret\
rieve the status of the job. \n  Usage: $scriptNam\
e --polljob --jobid <jobId> [--outfile string]\n\n\
[Help]\n\n  For more detailed help information ref\
er to\n  http://www.ebi.ac.uk/DaliLite/\nEOF\n;\n}\
\n","my $WSDL = 'http://www.ebi.ac.uk/Tools/webser\
vices/wsdl/WSWUBlast.wsdl';\n\nuse strict;\nuse SO\
AP::Lite;\nuse Getopt::Long qw(:config no_ignore_c\
ase bundling);\nuse File::Basename;\n\nmy $checkIn\
terval = 15;\n\nmy $numOpts = scalar(@ARGV);\nmy (\
$outfile, $outformat, $help, $async, $polljob, $st\
atus, $ids, $jobid, $trace, $sequence);\nmy %param\
s= ( # Defaults\n	      'async' => 1, # Force into\
 async mode\n	      'exp' => 10.0, # E-value thres\
hold\n	      'numal' => 50, # Maximum number of al\
ignments\n	      'scores' => 100, # Maximum number\
 of scores\n            );\nGetOptions( # Map the \
options into variables\n    \"program|p=s\"     =>\
 \\$params{program}, # BLAST program\n    \"databa\
se|D=s\"    => \\$params{database}, # Search datab\
ase\n    \"matrix|m=s\"      => \\$params{matrix},\
 # Scoring matrix\n    \"exp|E=f\"         => \\$p\
arams{exp}, # E-value threshold\n    \"echofilter|\
e\"    => \\$params{echofilter}, # Display filtere\
d sequence\n    \"filter|f=s\"      => \\$params{f\
ilter}, # Low complexity filter name\n    \"alignm\
ents|b=i\"  => \\$params{numal}, # Number of align\
ments\n    \"scores|s=i\"      => \\$params{scores\
}, # Number of scores\n    \"sensitivity|S=s\" => \
\\$params{sensitivity}, # Search sensitivity\n    \
\"sort|t=s\"	      => \\$params{sort}, # Sort hits\
 by...\n    \"stats|T=s\"       => \\$params{stats\
}, # Scoring statistic to use\n    \"strand|d=s\" \
     => \\$params{strand}, # Strand to use in DNA \
vs. DNA search\n    \"topcombon|c=i\"   => \\$para\
ms{topcombon}, # Consistent sets of HSPs\n    \"ou\
tfile=s\"       => \\$outfile, # Output file\n    \
\"outformat|o=s\"   => \\$outformat, # Output form\
at\n    \"help|h\"	      => \\$help, # Usage info\\
n    \"async|a\"	      => \\$async, # Asynchronous\
 mode\n    \"polljob\"	      => \\$polljob, # Get \
results\n    \"status\"	      => \\$status, # Get \
job status\n    \"ids\"             => \\$ids, # G\
et ids from result\n    \"jobid|j=s\"       => \\$\
jobid, # JobId\n    \"email=s\"         => \\$para\
ms{email}, # E-mail address\n    \"trace\"        \
   => \\$trace, # SOAP trace\n    \"sequence=s\"  \
    => \\$sequence, # Query sequence\n    );\n\nmy\
 $scriptName = basename($0, ());\nif($help || $num\
Opts == 0) {\n    &usage();\n    exit(0);\n}\n\nif\
($trace){\n    print STDERR \"Tracing active\\n\";\
\n    SOAP::Lite->import(+trace => 'debug');\n}\n\\
nmy $soap = SOAP::Lite\n    ->service($WSDL)\n    \
->proxy('http://localhost/',\n    #proxy => ['http\
' => 'http://your.proxy.server/'], # HTTP proxy\n \
   timeout => 600, # HTTP connection timeout\n    \
)\n    ->on_fault(sub { # SOAP fault handler\n    \
    my $soap = shift;\n        my $res = shift;\n \
       # Throw an exception for all faults\n      \
  if(ref($res) eq '') {\n            die($res);\n \
       } else {\n            die($res->faultstring\
);\n        }\n        return new SOAP::SOM;\n    \
}\n               );\n\nif( !($polljob || $status \
|| $ids) &&\n    !( defined($ARGV[0]) || defined($\
sequence) )\n    ) {\n    print STDERR 'Error: bad\
 option combination', \"\\n\";\n    &usage();\n   \
 exit(1);\n}\nelsif($polljob && defined($jobid)) {\
\n    print \"Getting results for job $jobid\\n\";\
\n    getResults($jobid);\n}\nelsif($status && def\
ined($jobid)) {\n    print STDERR \"Getting status\
 for job $jobid\\n\";\n    my $result = $soap->che\
ckStatus($jobid);\n    print STDOUT \"$result\\n\"\
;\n    if($result eq 'DONE') {\n	print STDERR \"To\
 get results: $scriptName --polljob --jobid $jobid\
\\n\";\n    }\n}  \nelsif($ids && defined($jobid))\
 {\n    print STDERR \"Getting ids from job $jobid\
\\n\";\n    getIds($jobid);\n}\nelse {\n    # Prep\
are input data\n    my $content;\n    my (@content\
s) = ();\n    if(-f $ARGV[0] || $ARGV[0] eq '-') {\
	\n	$content={type=>'sequence',content=>read_file(\
$ARGV[0])};	\n    }\n    if($sequence) {	\n	if(-f \
$sequence || $sequence eq '-') {	\n	    $content={\
type=>'sequence',content=>read_file($ARGV[0])};	\n\
	} else {\n	    $content={type=>'sequence',content\
=>$sequence};\n	}\n    }\n    push @contents, $con\
tent;\n\n    # Submit the job\n    my $paramsData \
= SOAP::Data->name('params')->type(map=>\\%params)\
;\n    my $contentData = SOAP::Data->name('content\
')->value(\\@contents);\n    # For SOAP::Lite 0.60\
 and earlier parameters are passed directly\n    i\
f($SOAP::Lite::VERSION eq '0.60' || $SOAP::Lite::V\
ERSION =~ /0\\.[1-5]/) {\n        $jobid = $soap->\
runWUBlast($paramsData, $contentData);\n    }\n   \
 # For SOAP::Lite 0.69 and later parameter handlin\
g is different, so pass\n    # undef's for templat\
ed params, and then pass the formatted args.\n    \
else {\n        $jobid = $soap->runWUBlast(undef, \
undef,\n				   $paramsData, $contentData);\n    }\\
n\n    # Asynchronous mode: output jobid and exit.\
\n    if (defined($async)) {\n	print STDOUT $jobid\
, \"\\n\";\n        print STDERR \"To check status\
: $scriptName --status --jobid $jobid\\n\";\n    }\
\n    # Synchronous mode: try to get results\n    \
else {\n        print STDERR \"JobId: $jobid\\n\";\
\n        sleep 1;\n        getResults($jobid);\n \
   }\n}\n\nsub getIds($) {\n    my $jobid = shift;\
\n    my $results = $soap->getIds($jobid);\n    fo\
r my $result (@$results){\n	print \"$result\\n\";\\
n    }\n}\n\nsub clientPoll($) {\n    my $jobid = \
shift;\n    my $result = 'PENDING';\n    # Check s\
tatus and wait if not finished\n    while($result \
eq 'RUNNING' || $result eq 'PENDING') {\n        $\
result = $soap->checkStatus($jobid);\n        prin\
t STDERR \"$result\\n\";\n        if($result eq 'R\
UNNING' || $result eq 'PENDING') {\n            # \
Wait before polling again.\n            sleep $che\
ckInterval;\n        }\n    }\n}\n\nsub getResults\
($) {\n    my $jobid = shift;\n    my $res;\n    #\
 Check status, and wait if not finished\n    clien\
tPoll($jobid);\n    # Use JobId if output file nam\
e is not defined\n    unless(defined($outfile)) {\\
n        $outfile=$jobid;\n    }\n    # Get list o\
f data types\n    my $resultTypes = $soap->getResu\
lts($jobid);\n    # Get the data and write it to a\
 file\n    if(defined($outformat)) { # Specified d\
ata type\n	if($outformat eq 'xml') {$outformat = '\
toolxml';}\n	if($outformat eq 'txt') {$outformat =\
 'tooloutput';}\n        my $selResultType;\n     \
   foreach my $resultType (@$resultTypes) {\n     \
       if($resultType->{type} eq $outformat) {\n  \
              $selResultType = $resultType;\n     \
       }\n        }\n        $res=$soap->poll($job\
id, $selResultType->{type});\n	if($outfile eq '-')\
 {\n	     write_file($outfile, $res);\n	} else {\n\
	    write_file($outfile.'.'.$selResultType->{ext}\
, $res);\n	}\n    } else { # Data types available\\
n        # Write a file for each output type\n    \
    for my $resultType (@$resultTypes){\n         \
   #print STDERR \"Getting $resultType->{type}\\n\\
";\n            $res=$soap->poll($jobid, $resultTy\
pe->{type});\n	    if($outfile eq '-') {\n		write_\
file($outfile, $res);\n	    } else {\n		write_file\
($outfile.'.'.$resultType->{ext}, $res);\n	    }\n\
        }\n    }\n}\n\nsub read_file($) {\n    my \
$filename = shift;\n    my ($content, $buffer);\n \
   if($filename eq '-') {\n	while(sysread(STDIN, $\
buffer, 1024)) {\n	    $content .= $buffer;\n	}\n \
   }\n    else { # File\n	open(FILE, $filename) or\
 die \"Error: unable to open input file\";\n	while\
(sysread(FILE, $buffer, 1024)) {\n	    $content .=\
 $buffer;\n	}\n	close(FILE);\n    }\n    return $c\
ontent;\n}\n\nsub write_file($$) {\n    my ($filen\
ame, $data) = @_;\n    print STDERR 'Creating resu\
lt file: ' . $filename . \"\\n\";\n    if($filenam\
e eq '-') {\n	print STDOUT $data;\n    }\n    else\
 {\n	open(FILE, \">$filename\") or die \"Error: un\
able to open output file\";\n	syswrite(FILE, $data\
);\n	close(FILE);\n    }\n}\n\nsub usage {\n    pr\
int STDERR <<EOF\nWU-BLAST\n========\n\nRapid sequ\
ence database search programs utilizing the BLAST \
algorithm.\n   \n[Required]\n\n      --email      \
 : str  : user email address \n  -p, --program	   \
 : str  : BLAST program to use: blastn, blastp, bl\
astx, \n                             tblastn or tb\
lastx\n  -D, --database    : str  : database to se\
arch\n  seqFile           : file : query sequence \
data file (\"-\" for STDIN)\n\n[Optional]\n\n  -m,\
 --matrix	    : str  : scoring matrix\n  -E, --exp\
	    : real : 0<E<= 1000. Statistical significance\
 threshold\n                             for repor\
ting database sequence matches.\n  -e, --echofilte\
r  :      : display the filtered query sequence in\
 the output\n  -f, --filter	    : str  : activates\
 filtering of the query sequence\n  -b, --alignmen\
ts  : int  : number of alignments to be reported\n\
  -s, --scores	    : int  : number of scores to be\
 reported\n  -S, --sensitivity : str  :\n  -t, --s\
ort	    : str  :\n  -T, --stats       : str  :\n  \
-d, --strand      : str  : DNA strand to search wi\
th in DNA vs. DNA searches \n  -c, --topcombon   :\
      :\n\n[General]	\n\n  -h, --help       :     \
 : prints this help text\n  -a, --async      :    \
  : forces to make an asynchronous query\n      --\
status     :      : poll for the status of a job\n\
      --polljob    :      : poll for the results o\
f a job\n  -j, --jobid      : str  : jobid that wa\
s returned when an asynchronous job \n            \
                was submitted.\n  -O, --outfile   \
 : str  : name of the file results should be writt\
en to \n                            (default is ba\
sed on the jobid; \"-\" for STDOUT)\n  -o, --outfo\
rmat  : str  : txt or xml output (no file is writt\
en)\n      --trace	   :      : show SOAP messages \
being interchanged \n\nSynchronous job:\n\n  The r\
esults/errors are returned as soon as the job is f\
inished.\n  Usage: $scriptName --email <your\\@ema\
il> [options...] seqFile\n  Returns: saves the res\
ults to disk\n\nAsynchronous job:\n\n  Use this if\
 you want to retrieve the results at a later time.\
 The results \n  are stored for up to 24 hours. \n\
  The asynchronous submission mode is recommended \
when users are submitting \n  batch jobs or large \
database searches	\n  Usage: $scriptName --async -\
-email <your\\@email> [options...] seqFile\n  Retu\
rns : jobid\n\n  Use the jobid to query for the st\
atus of the job. \n  Usage: $scriptName --status -\
-jobid <jobId>\n  Returns : string indicating the \
status of the job:\n    DONE - job has finished\n \
   RUNNING - job is running\n    NOT_FOUND - job c\
annot be found\n    ERROR - the jobs has encounter\
ed an error\n\n  When done, use the jobid to retri\
eve the status of the job. \n  Usage: $scriptName \
--polljob --jobid <jobId> [--outfile string]\n  Re\
turns: saves the results to disk\n\n[Help]\n\nFor \
more detailed help information refer to \nhttp://w\
ww.ebi.ac.uk/blast2/WU-Blast2_Help_frame.html\n \n\
EOF\n;\n}\n","\nmy $WSDL = 'http://www.ebi.ac.uk/T\
ools/webservices/wsdl/WSBlastpgp.wsdl';\n\nuse SOA\
P::Lite;\nuse Getopt::Long qw(:config no_ignore_ca\
se bundling);\nuse File::Basename;\n\nmy $checkInt\
erval = 15;\n\nmy %params=(\n	    'async' => '1', \
# Use async mode and simulate sync mode in client\\
n	    );\nGetOptions(\n    \"mode=s\"           =>\
 \\$params{mode}, # Search mode: PSI-Blast or PHI-\
Blast\n    \"database|d=s\"     => \\$params{datab\
ase}, # Database to search\n    \"matrix|M=s\"    \
   => \\$params{matrix},# Scoring maxtrix\n    \"e\
xp|e=f\"          => \\$params{exp}, # E-value\n  \
  \"expmulti|h=f\"     => \\$params{expmulti}, # E\
-value\n    \"filter|F=s\"       => \\$params{filt\
er}, # Low complexity filter\n    \"dropoff|X=i\" \
     => \\$params{dropoff}, # Dropoff score\n    \\
"finaldropoff|Z=i\" => \\$params{finaldropoff}, # \
Final dropoff score\n    \"scores|v=i\"       => \\
\$params{scores}, # Max number of scores\n    \"al\
ign=i\"          => \\$params{align}, # Alignment \
view\n    \"startregion|S=i\"  => \\$params{startr\
egion}, # Start of region in query\n    \"endregio\
n|H=i\"    => \\$params{endregion}, # End of regio\
n in query\n    \"maxpasses|j=i\"    => \\$params{\
maxpasses}, # Number of PSI iterations\n    \"open\
gap|G=i\"      => \\$params{opengap}, # Gap open p\
enalty\n    \"extendgap|E=i\"    => \\$params{exte\
ndgap}, # Gap extension penalty\n    \"pattern=s\"\
        => \\$params{pattern}, # PHI-BLAST pattern\
\n    \"usagemode|p=s\"    => \\$params{usagemode}\
, # PHI-BLAST program\n    \"appxml=s\"         =>\
 \\$params{appxml}, # Application XML\n    \"seque\
nce=s\"       => \\$sequence, # Query sequence\n  \
  \"help\"	       => \\$help, # Usage info\n    \"\
polljob\"	       => \\$polljob, # Get results\n   \
 \"status\"	       => \\$status, # Get status\n   \
 \"ids\"      	       => \\$ids, # Get ids from re\
sult\n    \"jobid=s\"          => \\$jobid, # JobI\
d\n    \"outfile=s\"        => \\$outfile, # Outpu\
t filename\n    \"outformat|o=s\"    => \\$outform\
at, # Output file format\n    \"async|a\"	       =\
> \\$async, # Async submission\n    \"email=s\"   \
       => \\$params{email}, # User e-mail address\\
n    \"trace\"            => \\$trace, # Show SOAP\
 messages\n    );\n\nmy $scriptName = basename($0,\
 ());\nif($help) {\n    &usage();\n    exit(0);\n}\
\n\nif ($trace){\n    print \"Tracing active\\n\";\
\n    SOAP::Lite->import(+trace => 'debug');\n}\n\\
nmy $soap = SOAP::Lite\n    ->service($WSDL)\n    \
->on_fault(sub {\n        my $soap = shift;\n     \
   my $res = shift;\n        # Throw an exception \
for all faults\n        if(ref($res) eq '') {\n   \
         die($res);\n        } else {\n           \
 die($res->faultstring);\n        }\n        retur\
n new SOAP::SOM;\n    }\n               );\n\nif( \
!($polljob || $status || $ids) &&\n    !( (defined\
($ARGV[0]) && -f $ARGV[0]) || defined($sequence) )\
\n    ) {\n    print STDERR 'Error: bad option com\
bination', \"\\n\";\n    &usage();\n    exit(1);\n\
}\nelsif($polljob && defined($jobid)) {\n    print\
 \"Getting results for job $jobid\\n\";\n    getRe\
sults($jobid);\n}\nelsif($status && defined($jobid\
)) {\n    print STDERR \"Getting status for job $j\
obid\\n\";\n    my $result = $soap->checkStatus($j\
obid);\n    print STDOUT $result, \"\\n\";\n    if\
($result eq 'DONE') {\n	print STDERR \"To get resu\
lts: $scriptName --polljob --jobid $jobid\\n\";\n \
   }\n}  \nelsif($ids && defined($jobid)) {\n    p\
rint STDERR \"Getting ids from job $jobid\\n\";\n \
   getIds($jobid);\n}\nelse {\n    if(-f $ARGV[0])\
 {	\n	$content={type=>'sequence', content=>read_fi\
le($ARGV[0])};	\n    }\n    if($sequence) {	\n	if(\
-f $sequence) {\n	    $content={type=>'sequence', \
content=>read_file($sequence)};	\n	} else {\n	    \
$content={type=>'sequence', content=>$sequence};\n\
	}\n    }\n    push @content, $content;\n\n    my \
$jobid;\n    my $paramsData = SOAP::Data->name('pa\
rams')->type(map=>\\%params);\n    my $contentData\
 = SOAP::Data->name('content')->value(\\@content);\
\n    # For SOAP::Lite 0.60 and earlier parameters\
 are passed directly\n    if($SOAP::Lite::VERSION \
eq '0.60' || $SOAP::Lite::VERSION =~ /0\\.[1-5]/) \
{\n        $jobid = $soap->runBlastpgp($paramsData\
, $contentData);\n    }\n    # For SOAP::Lite 0.69\
 and later parameter handling is different, so pas\
s\n    # undef's for templated params, and then pa\
ss the formatted args.\n    else {\n        $jobid\
 = $soap->runBlastpgp(undef, undef,\n				    $para\
msData, $contentData);\n    }\n\n    if (defined($\
async)) {\n	print STDOUT $jobid, \"\\n\";\n       \
 print STDERR \"To check status: $scriptName --sta\
tus --jobid $jobid\\n\";\n    } else { # Synchrono\
us mode\n        print STDERR \"JobId: $jobid\\n\"\
;\n        sleep 1;\n        getResults($jobid);\n\
    }\n}\n\nsub getIds($) {\n    $jobid = shift;\n\
    my $results = $soap->getIds($jobid);\n    for \
$result (@$results){\n	print \"$result\\n\";\n    \
}\n}\n\nsub clientPoll($) {\n    my $jobid = shift\
;\n    my $result = 'PENDING';\n    # Check status\
 and wait if not finished\n    #print STDERR \"Che\
cking status: $jobid\\n\";\n    while($result eq '\
RUNNING' || $result eq 'PENDING') {\n        $resu\
lt = $soap->checkStatus($jobid);\n        print ST\
DERR \"$result\\n\";\n        if($result eq 'RUNNI\
NG' || $result eq 'PENDING') {\n            # Wait\
 before polling again.\n            sleep $checkIn\
terval;\n        }\n    }\n}\n\nsub getResults($) \
{\n    $jobid = shift;\n    # Check status, and wa\
it if not finished\n    clientPoll($jobid);\n    #\
 Use JobId if output file name is not defined\n   \
 unless(defined($outfile)) {\n        $outfile=$jo\
bid;\n    }\n    # Get list of data types\n    my \
$resultTypes = $soap->getResults($jobid);\n    # G\
et the data and write it to a file\n    if(defined\
($outformat)) { # Specified data type\n        my \
$selResultType;\n        foreach my $resultType (@\
$resultTypes) {\n            if($resultType->{type\
} eq $outformat) {\n                $selResultType\
 = $resultType;\n            }\n        }\n       \
 $res=$soap->poll($jobid, $selResultType->{type});\
\n        write_file($outfile.'.'.$selResultType->\
{ext}, $res);\n    } else { # Data types available\
\n        # Write a file for each output type\n   \
     for my $resultType (@$resultTypes){\n        \
    #print \"Getting $resultType->{type}\\n\";\n  \
          $res=$soap->poll($jobid, $resultType->{t\
ype});\n            write_file($outfile.'.'.$resul\
tType->{ext}, $res);\n        }\n    }\n}\n\nsub r\
ead_file($) {\n    my $filename = shift;\n    open\
(FILE, $filename);\n    my $content;\n    my $buff\
er;\n    while(sysread(FILE, $buffer, 1024)) {\n	$\
content.= $buffer;\n    }\n    close(FILE);  \n   \
 return $content;\n}\n\nsub write_file($$) {\n    \
my ($tmp,$entity) = @_;\n    print STDERR \"Creati\
ng result file: \".$tmp.\"\\n\";\n    unless(open \
(FILE, \">$tmp\")) {\n	return 0;\n    }\n    syswr\
ite(FILE, $entity);\n    close (FILE);\n    return\
 1;\n}\n\nsub usage {\n    print STDERR <<EOF\nBla\
stpgp\n========\n   \nThe blastpgp program impleme\
nts the PSI-BLAST and PHI-BLAST variations\nof NCB\
I BLAST.\n\nFor more detailed help information ref\
er to\nhttp://www.ebi.ac.uk/blastpgp/blastpsi_help\
_frame.html\n \nBlastpgp specific options:\n\n[Req\
uired]\n\n      --mode            : str  : search \
mode to use: PSI-Blast or PHI-Blast\n  -d, --datab\
ase        : str  : protein database to search\n  \
seqFile               : file : query sequence\n\n[\
Optional]\n\n  -M, --matrix          : str  : scor\
ing matrix\n  -e, --exp             : real : Expec\
tation value\n  -h, --expmulti        : real : thr\
eshold (multipass model)\n  -F, --filter          \
: str  : filter query sequence with SEG [T,F]\n  -\
m, --align           : int  : alignment view optio\
n:\n                                 0 - pairwise,\
 1 - M/S identities,\n                            \
     2 - M/S non-identities, 3 - Flat identities,\\
n                                 4 - Flat non-ide\
ntities\n  -G, --opengap         : int  : cost to \
open a gap\n  -E, --extendgap       : int  : cost \
to extend a gap\n  -g, --gapalign        : str  : \
Gapped [T,F]\n  -v, --scores          : int  : num\
ber of scores to be reported\n  -j, --maxpasses   \
    : int  : number of iterations\n  -X, --dropoff\
         : int  : Dropoff score\n  -Z, --finaldrop\
off    : int  : Dropoff for final alignment\n  -S,\
 --startregion     : int  : Start of required regi\
on in query\n  -H, --endregion       : int  : End \
of required region in query\n  -k, --pattern      \
   : str  : Hit File (PHI-BLAST only)\n  -p, --usa\
gemode       : str  : Program option (PHI-BLAST on\
ly):\n                                 blastpgp, p\
atseedp, seedp\n\n[General]\n\n      --help       \
     :      : prints this help text\n  -a, --async\
           :      : forces to make an asynchronous\
 query\n      --status          :      : poll for \
the status of a job\n      --polljob         :    \
  : poll for the results of a job\n      --jobid  \
         : str  : jobid of an asynchronous job\n  \
    --ids             :      : get hit identifiers\
 for result \n  -O, --outfile         : str  : nam\
e of the file results should be written to\n      \
                           (default is based on th\
e jobid)\n  -o, --outformat       : str  : txt or \
xml output (no file is written)\n      --trace    \
       :      : show SOAP messages being interchan\
ged\n\nSynchronous job:\n\n  The results/errors ar\
e returned as soon as the job is finished.\n  Usag\
e: blastpgp.pl --email <your@email> [options...] s\
eqfile\n  Returns: saves the results to disk\n\nAs\
ynchronous job:\n\n  Use this if you want to retri\
eve the results at a later time. The results\n  ar\
e stored for up to 24 hours.\n  The asynchronous s\
ubmission mode is recommended when users are submi\
tting\n  batch jobs or large database searches\n  \
Usage: blastpgp.pl --email <your@email> --async [o\
ptions...] seqFile\n  Returns: jobid\n\n  Use the \
jobid to query for the status of the job.\n  Usage\
: blastpgp.pl --status --jobid <jobId>\n  Returns:\
 string indicating the status of the job\n    DONE\
 - job has finished\n    RUNNING - job is running\\
n    NOT_FOUND - job cannot be found\n    ERROR - \
the jobs has encountered an error\n\n  When done, \
use the jobid to retrieve the results of the job.\\
n  Usage: blastpgp.pl --polljob --jobid <jobId> [-\
-outfile <fileName>]\n  Returns: saves the results\
 to disk\nEOF\n;\n}\n","\n=head1 NAME\n\nncbiblast\
_lwp.pl\n\n=head1 DESCRIPTION\n\nNCBI BLAST REST w\
eb service Perl client using L<LWP>.\n\nTested wit\
h:\n\n=over\n\n=item *\nL<LWP> 5.79, L<XML::Simple\
> 2.12 and Perl 5.8.3\n\n=item *\nL<LWP> 5.805, L<\
XML::Simple> 2.14 and Perl 5.8.7\n\n=item *\nL<LWP\
> 5.820, L<XML::Simple> 2.18 and Perl 5.10.0 (Ubun\
tu 9.04)\n\n=back\n\nFor further information see:\\
n\n=over\n\n=item *\nL<http://www.ebi.ac.uk/Tools/\
webservices/services/sss/ncbi_blast_rest>\n\n=item\
 *\nL<http://www.ebi.ac.uk/Tools/webservices/tutor\
ials/perl>\n\n=back\n\n=head1 VERSION\n\n$Id: ncbi\
blast_lwp.pl 1317 2009-09-03 15:44:11Z hpm $\n\n=c\
ut\n\nuse strict;\nuse warnings;\n\nuse English;\n\
use LWP;\nuse XML::Simple;\nuse Getopt::Long qw(:c\
onfig no_ignore_case bundling);\nuse File::Basenam\
e;\nuse Data::Dumper;\n\nmy $baseUrl = 'http://www\
.ebi.ac.uk/Tools/services/rest/ncbiblast';\n\nmy $\
checkInterval = 3;\n\nmy $outputLevel = 1;\n\nmy $\
numOpts = scalar(@ARGV);\nmy %params = ( 'debugLev\
el' => 0 );\n\nmy %tool_params = ();\nGetOptions(\\
n\n	# Tool specific options\n	'program|p=s'  => \\\
$tool_params{'program'},   # blastp, blastn, blast\
x, etc.\n	'database|D=s' => \\$params{'database'},\
       # Database(s) to search\n	'matrix|m=s'   =>\
 \\$tool_params{'matrix'},    # Scoring martix to \
use\n	'exp|E=f'      => \\$tool_params{'exp'},    \
   # E-value threshold\n	'filter|f=s'   => \\$tool\
_params{'filter'},    # Low complexity filter\n	'a\
lign|A=i'    => \\$tool_params{'align'},     # Pai\
rwise alignment format\n	'scores|s=i'   => \\$tool\
_params{'scores'},    # Number of scores\n	'alignm\
ents|n=i' => \\$tool_params{'alignments'},   # Num\
ber of alignments\n	'dropoff|d=i'    => \\$tool_pa\
rams{'dropoff'},      # Dropoff score\n	'match_sco\
res=s' => \\$tool_params{'match_scores'}, # Match/\
missmatch scores\n	'match|u=i'      => \\$params{'\
match'},             # Match score\n	'mismatch|v=i\
'   => \\$params{'mismatch'},          # Mismatch \
score\n	'gapopen|o=i'    => \\$tool_params{'gapope\
n'},      # Open gap penalty\n	'gapext|x=i'     =>\
 \\$tool_params{'gapext'},       # Gap extension p\
enality\n	'gapalign|g'     => \\$tool_params{'gapa\
lign'},     # Optimise gap alignments\n	'stype=s' \
=> \\$tool_params{'stype'},    # Sequence type\n	'\
seqrange=s' => \\$tool_params{'seqrange'},    # Qu\
ery subsequence\n	'sequence=s' => \\$params{'seque\
nce'},         # Query sequence\n	'multifasta' => \
\\$params{'multifasta'},       # Multiple fasta in\
put\n\n	# Compatability options, old command-line\\
n	'numal|n=i'     => \\$params{'numal'},        # \
Number of alignments\n	'opengap|o=i'   => \\$param\
s{'opengap'},      # Open gap penalty\n	'extendgap\
|x=i' => \\$params{'extendgap'},    # Gap extensio\
n penality\n	\n	# Generic options\n	'email=s'     \
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
	$tool_params{'database'} = \\@dbList;\n\n	# Match\
/missmatch\n	if ( $params{'match'} && $params{'mis\
smatch'} ) {\n		$tool_params{'match_scores'} =\n		\
  $params{'match'} . ',' . $params{'missmatch'};\n\
	}\n	\n	# Compatability options, old command-line\\
n	if(!$tool_params{'alignments'} && $params{'numal\
'}) {\n		$tool_params{'alignments'} = $params{'num\
al'};\n	}\n	if(!$tool_params{'gapopen'} && $params\
{'opengap'}) {\n		$tool_params{'gapopen'} = $param\
s{'opengap'};\n	}\n	if(!$tool_params{'gapext'} && \
$params{'extendgap'}) {\n		$tool_params{'gapext'} \
= $params{'extendgap'};\n	}\n\n	print_debug_messag\
e( 'load_params', 'End', 1 );\n}\n\n=head2 client_\
poll()\n\nClient-side job polling.\n\n  &client_po\
ll($job_id);\n\n=cut\n\nsub client_poll {\n	print_\
debug_message( 'client_poll', 'Begin', 1 );\n	my $\
jobid  = shift;\n	my $status = 'PENDING';\n\n	my $\
errorCount = 0;\n	while ($status eq 'RUNNING'\n		|\
| $status eq 'PENDING'\n		|| ( $status eq 'ERROR' \
&& $errorCount < 2 ) )\n	{\n		$status = rest_get_s\
tatus($jobid);\n		print STDERR \"$status\\n\" if (\
 $outputLevel > 0 );\n		if ( $status eq 'ERROR' ) \
{\n			$errorCount++;\n		}\n		elsif ( $errorCount >\
 0 ) {\n			$errorCount--;\n		}\n		if (   $status e\
q 'RUNNING'\n			|| $status eq 'PENDING'\n			|| $st\
atus eq 'ERROR' )\n		{\n\n			# Wait before polling\
 again.\n			sleep $checkInterval;\n		}\n	}\n	print\
_debug_message( 'client_poll', 'End', 1 );\n	retur\
n $status;\n}\n\n=head2 get_results()\n\nGet the r\
esults for a job identifier.\n\n  &get_results($jo\
b_id);\n\n=cut\n\nsub get_results {\n	print_debug_\
message( 'get_results', 'Begin', 1 );\n	my $jobid \
= shift;\n	print_debug_message( 'get_results', 'jo\
bid: ' . $jobid, 1 );\n\n	# Verbose\n	if ( $output\
Level > 1 ) {\n		print 'Getting results for job ',\
 $jobid, \"\\n\";\n	}\n\n	# Check status, and wait\
 if not finished\n	client_poll($jobid);\n\n	# Use \
JobId if output file name is not defined\n	unless \
( defined( $params{'outfile'} ) ) {\n		$params{'ou\
tfile'} = $jobid;\n	}\n\n	# Get list of data types\
\n	my (@resultTypes) = rest_get_result_types($jobi\
d);\n\n	# Get the data and write it to a file\n	if\
 ( defined( $params{'outformat'} ) ) {    # Specif\
ied data type\n		my $selResultType;\n		foreach my \
$resultType (@resultTypes) {\n			if ( $resultType-\
>{'identifier'} eq $params{'outformat'} ) {\n				$\
selResultType = $resultType;\n			}\n		}\n		if ( de\
fined($selResultType) ) {\n			my $result =\n			  r\
est_get_result( $jobid, $selResultType->{'identifi\
er'} );\n			if ( $params{'outfile'} eq '-' ) {\n		\
		write_file( $params{'outfile'}, $result );\n			}\
\n			else {\n				write_file(\n					$params{'outfil\
e'} . '.'\n					  . $selResultType->{'identifier'}\
 . '.'\n					  . $selResultType->{'fileSuffix'},\n\
					$result\n				);\n			}\n		}\n		else {\n			die \
'Error: unknown result format \"' . $params{'outfo\
rmat'} . '\"';\n		}\n	}\n	else {    # Data types a\
vailable\n		      # Write a file for each output t\
ype\n		for my $resultType (@resultTypes) {\n			if \
( $outputLevel > 1 ) {\n				print STDERR 'Getting \
', $resultType->{'identifier'}, \"\\n\";\n			}\n		\
	my $result = rest_get_result( $jobid, $resultType\
->{'identifier'} );\n			if ( $params{'outfile'} eq\
 '-' ) {\n				write_file( $params{'outfile'}, $res\
ult );\n			}\n			else {\n				write_file(\n					$pa\
rams{'outfile'} . '.'\n					  . $resultType->{'ide\
ntifier'} . '.'\n					  . $resultType->{'fileSuffi\
x'},\n					$result\n				);\n			}\n		}\n	}\n	print_\
debug_message( 'get_results', 'End', 1 );\n}\n\n=h\
ead2 read_file()\n\nRead a file into a scalar. The\
 special filename '-' can be used to read from \ns\
tandard input (STDIN).\n\n  my $data = &read_file(\
$filename);\n\n=cut\n\nsub read_file {\n	print_deb\
ug_message( 'read_file', 'Begin', 1 );\n	my $filen\
ame = shift;\n	print_debug_message( 'read_file', '\
filename: ' . $filename, 2 );\n	my ( $content, $bu\
ffer );\n	if ( $filename eq '-' ) {\n		while ( sys\
read( STDIN, $buffer, 1024 ) ) {\n			$content .= $\
buffer;\n		}\n	}\n	else {    # File\n		open( my $F\
ILE, '<', $filename )\n		  or die \"Error: unable \
to open input file $filename ($!)\";\n		while ( sy\
sread( $FILE, $buffer, 1024 ) ) {\n			$content .= \
$buffer;\n		}\n		close($FILE);\n	}\n	print_debug_m\
essage( 'read_file', 'End', 1 );\n	return $content\
;\n}\n\n=head2 write_file()\n\nWrite data to a fil\
e. The special filename '-' can be used to write t\
o \nstandard output (STDOUT).\n\n  &write_file($fi\
lename, $data);\n\n=cut\n\nsub write_file {\n	prin\
t_debug_message( 'write_file', 'Begin', 1 );\n	my \
( $filename, $data ) = @_;\n	print_debug_message( \
'write_file', 'filename: ' . $filename, 2 );\n	if \
( $outputLevel > 0 ) {\n		print STDERR 'Creating r\
esult file: ' . $filename . \"\\n\";\n	}\n	if ( $f\
ilename eq '-' ) {\n		print STDOUT $data;\n	}\n	el\
se {\n		open( my $FILE, '>', $filename )\n		  or d\
ie \"Error: unable to open output file $filename (\
$!)\";\n		syswrite( $FILE, $data );\n		close($FILE\
);\n	}\n	print_debug_message( 'write_file', 'End',\
 1 );\n}\n\n=head2 usage()\n\nPrint program usage \
message.\n\n  &usage();\n\n=cut\n\nsub usage {\n	p\
rint STDERR <<EOF\nNCBI BLAST\n==========\n   \nRa\
pid sequence database search programs utilizing th\
e BLAST algorithm\n    \n[Required]\n\n  -p, --pro\
gram      : str  : BLAST program to use, see --par\
amDetail program\n  -D, --database     : str  : da\
tabase(s) to search, space separated. See\n       \
                       --paramDetail database\n   \
   --stype        : str  : query sequence type, se\
e --paramDetail stype\n  seqFile            : file\
 : query sequence (\"-\" for STDIN, \\@filename fo\
r\n                              identifier list f\
ile)\n\n[Optional]\n\n  -m, --matrix       : str  \
: scoring matrix, see --paramDetail matrix\n  -e, \
--exp          : real : 0<E<= 1000. Statistical si\
gnificance threshold \n                           \
   for reporting database sequence matches.\n  -f,\
 --filter       :      : filter the query sequence\
 for low complexity \n                            \
  regions, see --paramDetail filter\n  -A, --align\
        : int  : pairwise alignment format, see --\
paramDetail align\n  -s, --scores       : int  : n\
umber of scores to be reported\n  -n, --alignments\
   : int  : number of alignments to report\n  -u, \
--match        : int  : Match score (BLASTN only)\\
n  -v, --mismatch     : int  : Mismatch score (BLA\
STN only)\n  -o, --gapopen      : int  : Gap open \
penalty\n  -x, --gapext       : int  : Gap extensi\
on penalty\n  -d, --dropoff      : int  : Drop-off\
\n  -g, --gapalign     :      : Optimise gapped al\
ignments\n      --seqrange     : str  : region wit\
hin input to use as query\n      --multifasta   : \
     : treat input as a set of fasta formatted seq\
uences\n\n[General]\n\n  -h, --help        :      \
: prints this help text\n      --async       :    \
  : forces to make an asynchronous query\n      --\
email       : str  : e-mail address\n      --title\
       : str  : title for job\n      --status     \
 :      : get job status\n      --resultTypes :   \
   : get available result types for job\n      --p\
olljob     :      : poll for the status of a job\n\
      --jobid       : str  : jobid that was return\
ed when an asynchronous job \n                    \
         was submitted.\n      --outfile     : str\
  : file name for results (default is jobid;\n    \
                         \"-\" for STDOUT)\n      \
--outformat   : str  : result format to retrieve\n\
      --params      :      : list input parameters\
\n      --paramDetail : str  : display details for\
 input parameter\n      --quiet       :      : dec\
rease output\n      --verbose     :      : increas\
e output\n      --trace       :      : show SOAP m\
essages being interchanged \n   \nSynchronous job:\
\n\n  The results/errors are returned as soon as t\
he job is finished.\n  Usage: $scriptName --email \
<your\\@email> [options...] seqFile\n  Returns: re\
sults as an attachment\n\nAsynchronous job:\n\n  U\
se this if you want to retrieve the results at a l\
ater time. The results \n  are stored for up to 24\
 hours. 	\n  Usage: $scriptName --async --email <y\
our\\@email> [options...] seqFile\n  Returns: jobi\
d\n\n  Use the jobid to query for the status of th\
e job. If the job is finished, \n  it also returns\
 the results/errors.\n  Usage: $scriptName --pollj\
ob --jobid <jobId> [--outfile string]\n  Returns: \
string indicating the status of the job and if app\
licable, results \n  as an attachment.\n\nFurther \
information:\n\n  http://www.ebi.ac.uk/Tools/webse\
rvices/services/sss/ncbi_blast_rest\n  http://www.\
ebi.ac.uk/Tools/webservices/tutorials/perl\n\nSupp\
ort/Feedback:\n\n  http://www.ebi.ac.uk/support/\n\
EOF\n}\n\n=head1 FEEDBACK/SUPPORT\n\nPlease contac\
t us at L<http://www.ebi.ac.uk/support/> if you ha\
ve any \nfeedback, suggestions or issues with the \
service or this client.\n\n=cut\n","\n=head1 NAME\\
n\nwublast_lwp.pl\n\n=head1 DESCRIPTION\n\nWU-BLAS\
T REST web service Perl client using L<LWP>.\n\nTe\
sted with:\n\n=over\n\n=item *\nL<LWP> 5.79, L<XML\
::Simple> 2.12 and Perl 5.8.3\n\n=item *\nL<LWP> 5\
.805, L<XML::Simple> 2.14 and Perl 5.8.7\n\n=item \
*\nL<LWP> 5.820, L<XML::Simple> 2.18 and Perl 5.10\
.0 (Ubuntu 9.04)\n\n=back\n\nFor further informati\
on see:\n\n=over\n\n=item *\nL<http://www.ebi.ac.u\
k/Tools/webservices/services/sss/wu_blast_rest>\n\\
n=item *\nL<http://www.ebi.ac.uk/Tools/webservices\
/tutorials/perl>\n\n=back\n\n=head1 VERSION\n\n$Id\
: wublast_lwp.pl 1317 2009-09-03 15:44:11Z hpm $\n\
\n=cut\n\nuse strict;\nuse warnings;\n\nuse Englis\
h;\nuse LWP;\nuse XML::Simple;\nuse Getopt::Long q\
w(:config no_ignore_case bundling);\nuse File::Bas\
ename;\nuse Data::Dumper;\n\nmy $baseUrl = 'http:/\
/www.ebi.ac.uk/Tools/services/rest/wublast';\n\nmy\
 $checkInterval = 3;\n\nmy $outputLevel = 1;\n\nmy\
 $numOpts = scalar(@ARGV);\nmy %params = ( 'debugL\
evel' => 0 );\n\nmy %tool_params = ();\nGetOptions\
(\n\n	# Tool specific options\n	'program|p=s'     \
=> \\$tool_params{'program'},      # BLAST program\
\n	'database|D=s'    => \\$params{'database'},    \
 # Search database\n	'matrix|m=s'      => \\$tool_\
params{'matrix'},       # Scoring matrix\n	'exp|E=\
f'         => \\$tool_params{'exp'},          # E-\
value threshold\n	'viewfilter|e'    => \\$tool_par\
ams{'viewfilter'},   # Display filtered sequence\n\
	'filter|f=s'      => \\$tool_params{'filter'},   \
    # Low complexity filter name\n	'alignments|n=i\
'  => \\$tool_params{'alignments'},   # Number of \
alignments\n	'scores|s=i'      => \\$tool_params{'\
scores'},       # Number of scores\n	'sensitivity|\
S=s' => \\$tool_params{'sensitivity'},  # Search s\
ensitivity\n	'sort|t=s'        => \\$tool_params{'\
sort'},         # Sort hits by...\n	'stats|T=s'   \
    => \\$tool_params{'stats'},        # Scoring s\
tatistic to use\n	'strand|d=s'      => \\$tool_par\
ams{'strand'},       # Strand to use\n	'topcombon|\
c=i'   => \\$tool_params{'topcombon'},    # Consis\
tent sets of HSPs\n	'align|A=i'       => \\$tool_p\
arams{'align'},   # Pairwise alignment format\n	's\
type=s' => \\$tool_params{'stype'},    # Sequence \
type 'protein' or 'dna'\n	'sequence=s' => \\$param\
s{'sequence'},         # Query sequence file or DB\
:ID\n	'multifasta' => \\$params{'multifasta'},    \
   # Multiple fasta input\n\n	# Compatability opti\
ons, old command-line.\n	'echofilter|e'    => \\$p\
arams{'echofilter'},   # Display filtered sequence\
\n	'b=i'  => \\$params{'numal'},        # Number o\
f alignments\n	'appxml=s'        => \\$params{'app\
xml'},       # Application XML\n\n	# Generic optio\
ns\n	'email=s'       => \\$params{'email'},       \
   # User e-mail address\n	'title=s'       => \\$p\
arams{'title'},          # Job title\n	'outfile=s'\
     => \\$params{'outfile'},        # Output file\
 name\n	'outformat=s'   => \\$params{'outformat'},\
      # Output file type\n	'jobid=s'       => \\$p\
arams{'jobid'},          # JobId\n	'help|h'       \
 => \\$params{'help'},           # Usage help\n	'a\
sync'         => \\$params{'async'},          # As\
ynchronous submission\n	'polljob'       => \\$para\
ms{'polljob'},        # Get results\n	'resultTypes\
'   => \\$params{'resultTypes'},    # Get result t\
ypes\n	'status'        => \\$params{'status'},    \
     # Get status\n	'params'        => \\$params{'\
params'},         # List input parameters\n	'param\
Detail=s' => \\$params{'paramDetail'},    # Get de\
tails for parameter\n	'quiet'         => \\$params\
{'quiet'},          # Decrease output level\n	'ver\
bose'       => \\$params{'verbose'},        # Incr\
ease output level\n	'debugLevel=i'  => \\$params{'\
debugLevel'},     # Debug output level\n	'baseUrl=\
s'     => \\$baseUrl,                  # Base URL \
for service.\n);\nif ( $params{'verbose'} ) { $out\
putLevel++ }\nif ( $params{'$quiet'} )  { $outputL\
evel-- }\n\n&print_debug_message( 'MAIN', 'LWP::VE\
RSION: ' . $LWP::VERSION,\n	1 );\n\n&print_debug_m\
essage( 'MAIN', \"params:\\n\" . Dumper( \\%params\
 ),           11 );\n&print_debug_message( 'MAIN',\
 \"tool_params:\\n\" . Dumper( \\%tool_params ), 1\
1 );\n\nmy $scriptName = basename( $0, () );\n\nif\
 ( $params{'help'} || $numOpts == 0 ) {\n	&usage()\
;\n	exit(0);\n}\n\n&print_debug_message( 'MAIN', '\
baseUrl: ' . $baseUrl, 1 );\n\nif (\n	!(\n		   $pa\
rams{'polljob'}\n		|| $params{'resultTypes'}\n		||\
 $params{'status'}\n		|| $params{'params'}\n		|| $\
params{'paramDetail'}\n	)\n	&& !( defined( $ARGV[0\
] ) || defined( $params{'sequence'} ) )\n  )\n{\n\\
n	# Bad argument combination, so print error messa\
ge and usage\n	print STDERR 'Error: bad option com\
bination', \"\\n\";\n	&usage();\n	exit(1);\n}\n\ne\
lsif ( $params{'params'} ) {\n	&print_tool_params(\
);\n}\n\nelsif ( $params{'paramDetail'} ) {\n	&pri\
nt_param_details( $params{'paramDetail'} );\n}\n\n\
elsif ( $params{'status'} && defined( $params{'job\
id'} ) ) {\n	&print_job_status( $params{'jobid'} )\
;\n}\n\nelsif ( $params{'resultTypes'} && defined(\
 $params{'jobid'} ) ) {\n	&print_result_types( $pa\
rams{'jobid'} );\n}\n\nelsif ( $params{'polljob'} \
&& defined( $params{'jobid'} ) ) {\n	&get_results(\
 $params{'jobid'} );\n}\n\nelse {\n\n	# Multiple i\
nput sequence mode, assume fasta format.\n	if ( $p\
arams{'multifasta'} ) {\n		&multi_submit_job();\n	\
}\n\n	# Entry identifier list file.\n	elsif (( def\
ined( $params{'sequence'} ) && $params{'sequence'}\
 =~ m/^\\@/ )\n		|| ( defined( $ARGV[0] ) && $ARGV\
[0] =~ m/^\\@/ ) )\n	{\n		my $list_filename = $par\
ams{'sequence'} || $ARGV[0];\n		$list_filename =~ \
s/^\\@//;\n		&list_file_submit_job($list_filename)\
;\n	}\n\n	# Default: single sequence/identifier.\n\
	else {\n\n		# Load the sequence data and submit.\\
n		&submit_job( &load_data() );\n	}\n}\n\n=head1 F\
UNCTIONS\n\n=cut\n\n\n=head2 rest_request()\n\nPer\
form a REST request.\n\n  my $response_str = &rest\
_request($url);\n\n=cut\n\nsub rest_request {\n	pr\
int_debug_message( 'rest_request', 'Begin', 11 );\\
n	my $requestUrl = shift;\n	print_debug_message( '\
rest_request', 'URL: ' . $requestUrl, 11 );\n\n	# \
Create a user agent\n	my $ua = LWP::UserAgent->new\
();\n	'$Revision: 1317 $' =~ m/(\\d+)/;\n	$ua->age\
nt(\"EBI-Sample-Client/$1 ($scriptName; $OSNAME) \\
" . $ua->agent());\n	$ua->env_proxy;\n\n	# Perform\
 the request\n	my $response = $ua->get($requestUrl\
);\n	print_debug_message( 'rest_request', 'HTTP st\
atus: ' . $response->code,\n		11 );\n\n	# Check fo\
r HTTP error codes\n	if ( $response->is_error ) {\\
n		$response->content() =~ m/<h1>([^<]+)<\\/h1>/;\\
n		die 'http status: ' . $response->code . ' ' . $\
response->message . '  ' . $1;\n	}\n	print_debug_m\
essage( 'rest_request', 'End', 11 );\n\n	# Return \
the response data\n	return $response->content();\n\
}\n\n=head2 rest_get_parameters()\n\nGet list of t\
ool parameter names.\n\n  my (@param_list) = &rest\
_get_parameters();\n\n=cut\n\nsub rest_get_paramet\
ers {\n	print_debug_message( 'rest_get_parameters'\
, 'Begin', 1 );\n	my $url                = $baseUr\
l . '/parameters/';\n	my $param_list_xml_str = res\
t_request($url);\n	my $param_list_xml     = XMLin(\
$param_list_xml_str);\n	my (@param_list)       = @\
{ $param_list_xml->{'id'} };\n	print_debug_message\
( 'rest_get_parameters', 'End', 1 );\n	return (@pa\
ram_list);\n}\n\n=head2 rest_get_parameter_details\
()\n\nGet details of a tool parameter.\n\n  my $pa\
ramDetail = &rest_get_parameter_details($param_nam\
e);\n\n=cut\n\nsub rest_get_parameter_details {\n	\
print_debug_message( 'rest_get_parameter_details',\
 'Begin', 1 );\n	my $parameterId = shift;\n	print_\
debug_message( 'rest_get_parameter_details',\n		'p\
arameterId: ' . $parameterId, 1 );\n	my $url      \
            = $baseUrl . '/parameterdetails/' . $p\
arameterId;\n	my $param_detail_xml_str = rest_requ\
est($url);\n	my $param_detail_xml     = XMLin($par\
am_detail_xml_str);\n	print_debug_message( 'rest_g\
et_parameter_details', 'End', 1 );\n	return ($para\
m_detail_xml);\n}\n\n=head2 rest_run()\n\nSubmit a\
 job.\n\n  my $job_id = &rest_run($email, $title, \
\\%params );\n\n=cut\n\nsub rest_run {\n	print_deb\
ug_message( 'rest_run', 'Begin', 1 );\n	my $email \
 = shift;\n	my $title  = shift;\n	my $params = shi\
ft;\n	print_debug_message( 'rest_run', 'email: ' .\
 $email, 1 );\n	if ( defined($title) ) {\n		print_\
debug_message( 'rest_run', 'title: ' . $title, 1 )\
;\n	}\n	print_debug_message( 'rest_run', 'params: \
' . Dumper($params), 1 );\n\n	# User agent to perf\
orm http requests\n	my $ua = LWP::UserAgent->new()\
;\n	$ua->env_proxy;\n\n	# Clean up parameters\n	my\
 (%tmp_params) = %{$params};\n	$tmp_params{'email'\
} = $email;\n	$tmp_params{'title'} = $title;\n	for\
each my $param_name ( keys(%tmp_params) ) {\n		if \
( !defined( $tmp_params{$param_name} ) ) {\n			del\
ete $tmp_params{$param_name};\n		}\n	}\n\n	# Submi\
t the job as a POST\n	my $url = $baseUrl . '/run';\
\n	my $response = $ua->post( $url, \\%tmp_params )\
;\n	print_debug_message( 'rest_run', 'HTTP status:\
 ' . $response->code, 11 );\n	print_debug_message(\
 'rest_run',\n		'request: ' . $response->request()\
->content(), 11 );\n\n	# Check for HTTP error code\
s\n	if ( $response->is_error ) {\n		$response->con\
tent() =~ m/<h1>([^<]+)<\\/h1>/;\n		die 'http stat\
us: ' . $response->code . ' ' . $response->message\
 . '  ' . $1;\n	}\n\n	# The job id is returned\n	m\
y $job_id = $response->content();\n	print_debug_me\
ssage( 'rest_run', 'End', 1 );\n	return $job_id;\n\
}\n\n=head2 rest_get_status()\n\nCheck the status \
of a job.\n\n  my $status = &rest_get_status($job_\
id);\n\n=cut\n\nsub rest_get_status {\n	print_debu\
g_message( 'rest_get_status', 'Begin', 1 );\n	my $\
job_id = shift;\n	print_debug_message( 'rest_get_s\
tatus', 'jobid: ' . $job_id, 2 );\n	my $status_str\
 = 'UNKNOWN';\n	my $url        = $baseUrl . '/stat\
us/' . $job_id;\n	$status_str = &rest_request($url\
);\n	print_debug_message( 'rest_get_status', 'stat\
us_str: ' . $status_str, 2 );\n	print_debug_messag\
e( 'rest_get_status', 'End', 1 );\n	return $status\
_str;\n}\n\n=head2 rest_get_result_types()\n\nGet \
list of result types for finished job.\n\n  my (@r\
esult_types) = &rest_get_result_types($job_id);\n\\
n=cut\n\nsub rest_get_result_types {\n	print_debug\
_message( 'rest_get_result_types', 'Begin', 1 );\n\
	my $job_id = shift;\n	print_debug_message( 'rest_\
get_result_types', 'jobid: ' . $job_id, 2 );\n	my \
(@resultTypes);\n	my $url                      = $\
baseUrl . '/resulttypes/' . $job_id;\n	my $result_\
type_list_xml_str = &rest_request($url);\n	my $res\
ult_type_list_xml     = XMLin($result_type_list_xm\
l_str);\n	(@resultTypes) = @{ $result_type_list_xm\
l->{'type'} };\n	print_debug_message( 'rest_get_re\
sult_types',\n		scalar(@resultTypes) . ' result ty\
pes', 2 );\n	print_debug_message( 'rest_get_result\
_types', 'End', 1 );\n	return (@resultTypes);\n}\n\
\n=head2 rest_get_result()\n\nGet result data of a\
 specified type for a finished job.\n\n  my $resul\
t = rest_get_result($job_id, $result_type);\n\n=cu\
t\n\nsub rest_get_result {\n	print_debug_message( \
'rest_get_result', 'Begin', 1 );\n	my $job_id = sh\
ift;\n	my $type   = shift;\n	print_debug_message( \
'rest_get_result', 'jobid: ' . $job_id, 1 );\n	pri\
nt_debug_message( 'rest_get_result', 'type: ' . $t\
ype,    1 );\n	my $url    = $baseUrl . '/result/' \
. $job_id . '/' . $type;\n	my $result = &rest_requ\
est($url);\n	print_debug_message( 'rest_get_result\
', length($result) . ' characters',\n		1 );\n	prin\
t_debug_message( 'rest_get_result', 'End', 1 );\n	\
return $result;\n}\n\n\n=head2 print_debug_message\
()\n\nPrint debug message at specified debug level\
.\n\n  &print_debug_message($method_name, $message\
, $level);\n\n=cut\n\nsub print_debug_message {\n	\
my $function_name = shift;\n	my $message       = s\
hift;\n	my $level         = shift;\n	if ( $level <\
= $params{'debugLevel'} ) {\n		print STDERR '[', $\
function_name, '()] ', $message, \"\\n\";\n	}\n}\n\
\n=head2 print_tool_params()\n\nPrint list of tool\
 parameters.\n\n  &print_tool_params();\n\n=cut\n\\
nsub print_tool_params {\n	print_debug_message( 'p\
rint_tool_params', 'Begin', 1 );\n	my (@param_list\
) = &rest_get_parameters();\n	foreach my $param ( \
sort(@param_list) ) {\n		print $param, \"\\n\";\n	\
}\n	print_debug_message( 'print_tool_params', 'End\
', 1 );\n}\n\n=head2 print_param_details()\n\nPrin\
t details of a tool parameter.\n\n  &print_param_d\
etails($param_name);\n\n=cut\n\nsub print_param_de\
tails {\n	print_debug_message( 'print_param_detail\
s', 'Begin', 1 );\n	my $paramName = shift;\n	print\
_debug_message( 'print_param_details', 'paramName:\
 ' . $paramName, 2 );\n	my $paramDetail = &rest_ge\
t_parameter_details($paramName);\n	print $paramDet\
ail->{'name'}, \"\\t\", $paramDetail->{'type'}, \"\
\\n\";\n	print $paramDetail->{'description'}, \"\\\
n\";\n	foreach my $value ( @{ $paramDetail->{'valu\
es'}->{'value'} } ) {\n		print $value->{'value'};\\
n		if ( $value->{'defaultValue'} eq 'true' ) {\n		\
	print \"\\t\", 'default';\n		}\n		print \"\\n\";\\
n		print \"\\t\", $value->{'label'}, \"\\n\";\n	}\\
n	print_debug_message( 'print_param_details', 'End\
', 1 );\n}\n\n=head2 print_job_status()\n\nPrint s\
tatus of a job.\n\n  &print_job_status($job_id);\n\
\n=cut\n\nsub print_job_status {\n	print_debug_mes\
sage( 'print_job_status', 'Begin', 1 );\n	my $jobi\
d = shift;\n	print_debug_message( 'print_job_statu\
s', 'jobid: ' . $jobid, 1 );\n	if ( $outputLevel >\
 0 ) {\n		print STDERR 'Getting status for job ', \
$jobid, \"\\n\";\n	}\n	my $result = &rest_get_stat\
us($jobid);\n	print \"$result\\n\";\n	if ( $result\
 eq 'FINISHED' && $outputLevel > 0 ) {\n		print ST\
DERR \"To get results: $scriptName --polljob --job\
id \" . $jobid\n		  . \"\\n\";\n	}\n	print_debug_m\
essage( 'print_job_status', 'End', 1 );\n}\n\n=hea\
d2 print_result_types()\n\nPrint available result \
types for a job.\n\n  &print_result_types($job_id)\
;\n\n=cut\n\nsub print_result_types {\n	print_debu\
g_message( 'result_types', 'Begin', 1 );\n	my $job\
id = shift;\n	print_debug_message( 'result_types',\
 'jobid: ' . $jobid, 1 );\n	if ( $outputLevel > 0 \
) {\n		print STDERR 'Getting result types for job \
', $jobid, \"\\n\";\n	}\n	my $status = &rest_get_s\
tatus($jobid);\n	if ( $status eq 'PENDING' || $sta\
tus eq 'RUNNING' ) {\n		print STDERR 'Error: Job s\
tatus is ', $status,\n		  '. To get result types t\
he job must be finished.', \"\\n\";\n	}\n	else {\n\
		my (@resultTypes) = &rest_get_result_types($jobi\
d);\n		if ( $outputLevel > 0 ) {\n			print STDOUT \
'Available result types:', \"\\n\";\n		}\n		foreac\
h my $resultType (@resultTypes) {\n			print STDOUT\
 $resultType->{'identifier'}, \"\\n\";\n			if ( de\
fined( $resultType->{'label'} ) ) {\n				print STD\
OUT \"\\t\", $resultType->{'label'}, \"\\n\";\n			\
}\n			if ( defined( $resultType->{'description'} )\
 ) {\n				print STDOUT \"\\t\", $resultType->{'des\
cription'}, \"\\n\";\n			}\n			if ( defined( $resu\
ltType->{'mediaType'} ) ) {\n				print STDOUT \"\\\
t\", $resultType->{'mediaType'}, \"\\n\";\n			}\n	\
		if ( defined( $resultType->{'fileSuffix'} ) ) {\\
n				print STDOUT \"\\t\", $resultType->{'fileSuff\
ix'}, \"\\n\";\n			}\n		}\n		if ( $status eq 'FINI\
SHED' && $outputLevel > 0 ) {\n			print STDERR \"\\
\n\", 'To get results:', \"\\n\",\n			  \"  $scrip\
tName --polljob --jobid \" . $params{'jobid'} . \"\
\\n\",\n			  \"  $scriptName --polljob --outformat\
 <type> --jobid \"\n			  . $params{'jobid'} . \"\\\
n\";\n		}\n	}\n	print_debug_message( 'result_types\
', 'End', 1 );\n}\n\n=head2 submit_job()\n\nSubmit\
 a job to the service.\n\n  &submit_job($seq);\n\n\
=cut\n\nsub submit_job {\n	print_debug_message( 's\
ubmit_job', 'Begin', 1 );\n\n	# Set input sequence\
\n	$tool_params{'sequence'} = shift;\n\n	# Load pa\
rameters\n	&load_params();\n\n	# Submit the job\n	\
my $jobid = &rest_run( $params{'email'}, $params{'\
title'}, \\%tool_params );\n\n	# Simulate sync/asy\
nc mode\n	if ( defined( $params{'async'} ) ) {\n		\
print STDOUT $jobid, \"\\n\";\n		if ( $outputLevel\
 > 0 ) {\n			print STDERR\n			  \"To check status:\
 $scriptName --status --jobid $jobid\\n\";\n		}\n	\
}\n	else {\n		if ( $outputLevel > 0 ) {\n			print \
STDERR \"JobId: $jobid\\n\";\n		}\n		sleep 1;\n		&\
get_results($jobid);\n	}\n	print_debug_message( 's\
ubmit_job', 'End', 1 );\n}\n\n=head2 multi_submit_\
job()\n\nSubmit multiple jobs assuming input is a \
collection of fasta formatted sequences.\n\n  &mul\
ti_submit_job();\n\n=cut\n\nsub multi_submit_job {\
\n	print_debug_message( 'multi_submit_job', 'Begin\
', 1 );\n	my $jobIdForFilename = 1;\n	$jobIdForFil\
ename = 0 if ( defined( $params{'outfile'} ) );\n	\
my (@filename_list) = ();\n\n	# Query sequence\n	i\
f ( defined( $ARGV[0] ) ) {    # Bare option\n		if\
 ( -f $ARGV[0] || $ARGV[0] eq '-' ) {    # File\n	\
		push( @filename_list, $ARGV[0] );\n		}\n	}\n	if \
( $params{'sequence'} ) {                   # Via \
--sequence\n		if ( -f $params{'sequence'} || $para\
ms{'sequence'} eq '-' ) {    # File\n			push( @fil\
ename_list, $params{'sequence'} );\n		}\n	}\n\n	$/\
 = '>';\n	foreach my $filename (@filename_list) {\\
n		open( my $INFILE, '<', $filename )\n		  or die \
\"Error: unable to open file $filename ($!)\";\n		\
while (<$INFILE>) {\n			my $seq = $_;\n			$seq =~ \
s/>$//;\n			if ( $seq =~ m/(\\S+)/ ) {\n				print \
STDERR \"Submitting job for: $1\\n\"\n				  if ( $\
outputLevel > 0 );\n				$seq = '>' . $seq;\n				&p\
rint_debug_message( 'multi_submit_job', $seq, 11 )\
;\n				&submit_job($seq);\n				$params{'outfile'} \
= undef if ( $jobIdForFilename == 1 );\n			}\n		}\\
n		close $INFILE;\n	}\n	print_debug_message( 'mult\
i_submit_job', 'End', 1 );\n}\n\n=head2 list_file_\
submit_job()\n\nSubmit multiple jobs using a file \
containing a list of entry identifiers as \ninput.\
\n\n  &list_file_submit_job($list_filename)\n\n=cu\
t\n\nsub list_file_submit_job {\n	my $filename    \
     = shift;\n	my $jobIdForFilename = 1;\n	$jobId\
ForFilename = 0 if ( defined( $params{'outfile'} )\
 );\n\n	# Iterate over identifiers, submitting eac\
h job\n	open( my $LISTFILE, '<', $filename )\n	  o\
r die 'Error: unable to open file ' . $filename . \
' (' . $! . ')';\n	while (<$LISTFILE>) {\n		my $li\
ne = $_;\n		chomp($line);\n		if ( $line ne '' ) {\\
n			&print_debug_message( 'list_file_submit_job', \
'line: ' . $line, 2 );\n			if ( $line =~ m/\\w:\\w\
/ ) {    # Check this is an identifier\n				print \
STDERR \"Submitting job for: $line\\n\"\n				  if \
( $outputLevel > 0 );\n				&submit_job($line);\n		\
	}\n			else {\n				print STDERR\n\"Warning: line \\
\\"$line\\\" is not recognised as an identifier\\n\
\";\n			}\n		}\n		$params{'outfile'} = undef if ( \
$jobIdForFilename == 1 );\n	}\n	close $LISTFILE;\n\
}\n\n=head2 load_data()\n\nLoad sequence data from\
 file or option specified on the command-line.\n\n\
  &load_data();\n\n=cut\n\nsub load_data {\n	print\
_debug_message( 'load_data', 'Begin', 1 );\n	my $r\
etSeq;\n\n	# Query sequence\n	if ( defined( $ARGV[\
0] ) ) {    # Bare option\n		if ( -f $ARGV[0] || $\
ARGV[0] eq '-' ) {    # File\n			$retSeq = &read_f\
ile( $ARGV[0] );\n		}\n		else {                   \
                  # DB:ID or sequence\n			$retSeq \
= $ARGV[0];\n		}\n	}\n	if ( $params{'sequence'} ) \
{                   # Via --sequence\n		if ( -f $p\
arams{'sequence'} || $params{'sequence'} eq '-' ) \
{    # File\n			$retSeq = &read_file( $params{'seq\
uence'} );\n		}\n		else {    # DB:ID or sequence\n\
			$retSeq = $params{'sequence'};\n		}\n	}\n	print\
_debug_message( 'load_data', 'End', 1 );\n	return \
$retSeq;\n}\n\n=head2 load_params()\n\nLoad job pa\
rameters from command-line options.\n\n  &load_par\
ams();\n\n=cut\n\nsub load_params {\n	print_debug_\
message( 'load_params', 'Begin', 1 );\n\n	# Databa\
se(s) to search\n	my (@dbList) = split /[ ,]/, $pa\
rams{'database'};\n	$tool_params{'database'} = \\@\
dbList;\n\n	# Compatability options, old command-l\
ine.\n	if(!$tool_params{'viewfilter'} && $params{'\
echofilter'}) {\n		$tool_params{'viewfilter'} = 't\
rue';\n	}\n	if(!$tool_params{'alignments'} && $par\
ams{'numal'}) {\n		$tool_params{'alignments'} = $p\
arams{'numal'};\n	}\n	# TODO: set alignment format\
 option to get NCBI BLAST XML.\n	if($params{'appxm\
l'}) {\n		$tool_params{'align'} = '';\n	}\n\n	prin\
t_debug_message( 'load_params', 'End', 1 );\n}\n\n\
=head2 client_poll()\n\nClient-side job polling.\n\
\n  &client_poll($job_id);\n\n=cut\n\nsub client_p\
oll {\n	print_debug_message( 'client_poll', 'Begin\
', 1 );\n	my $jobid  = shift;\n	my $status = 'PEND\
ING';\n\n	my $errorCount = 0;\n	while ($status eq \
'RUNNING'\n		|| $status eq 'PENDING'\n		|| ( $stat\
us eq 'ERROR' && $errorCount < 2 ) )\n	{\n		$statu\
s = rest_get_status($jobid);\n		print STDERR \"$st\
atus\\n\" if ( $outputLevel > 0 );\n		if ( $status\
 eq 'ERROR' ) {\n			$errorCount++;\n		}\n		elsif (\
 $errorCount > 0 ) {\n			$errorCount--;\n		}\n		if\
 (   $status eq 'RUNNING'\n			|| $status eq 'PENDI\
NG'\n			|| $status eq 'ERROR' )\n		{\n\n			# Wait \
before polling again.\n			sleep $checkInterval;\n	\
	}\n	}\n	print_debug_message( 'client_poll', 'End'\
, 1 );\n	return $status;\n}\n\n=head2 get_results(\
)\n\nGet the results for a job identifier.\n\n  &g\
et_results($job_id);\n\n=cut\n\nsub get_results {\\
n	print_debug_message( 'get_results', 'Begin', 1 )\
;\n	my $jobid = shift;\n	print_debug_message( 'get\
_results', 'jobid: ' . $jobid, 1 );\n\n	# Verbose\\
n	if ( $outputLevel > 1 ) {\n		print 'Getting resu\
lts for job ', $jobid, \"\\n\";\n	}\n\n	# Check st\
atus, and wait if not finished\n	client_poll($jobi\
d);\n\n	# Use JobId if output file name is not def\
ined\n	unless ( defined( $params{'outfile'} ) ) {\\
n		$params{'outfile'} = $jobid;\n	}\n\n	# Get list\
 of data types\n	my (@resultTypes) = rest_get_resu\
lt_types($jobid);\n\n	# Get the data and write it \
to a file\n	if ( defined( $params{'outformat'} ) )\
 {    # Specified data type\n		my $selResultType;\\
n		foreach my $resultType (@resultTypes) {\n			if \
( $resultType->{'identifier'} eq $params{'outforma\
t'} ) {\n				$selResultType = $resultType;\n			}\n\
		}\n		if ( defined($selResultType) ) {\n			my $re\
sult =\n			  rest_get_result( $jobid, $selResultTy\
pe->{'identifier'} );\n			if ( $params{'outfile'} \
eq '-' ) {\n				write_file( $params{'outfile'}, $r\
esult );\n			}\n			else {\n				write_file(\n					$\
params{'outfile'} . '.'\n					  . $selResultType->\
{'identifier'} . '.'\n					  . $selResultType->{'f\
ileSuffix'},\n					$result\n				);\n			}\n		}\n		e\
lse {\n			die 'Error: unknown result format \"' . \
$params{'outformat'} . '\"';\n		}\n	}\n	else {    \
# Data types available\n		      # Write a file for\
 each output type\n		for my $resultType (@resultTy\
pes) {\n			if ( $outputLevel > 1 ) {\n				print ST\
DERR 'Getting ', $resultType->{'identifier'}, \"\\\
n\";\n			}\n			my $result = rest_get_result( $jobi\
d, $resultType->{'identifier'} );\n			if ( $params\
{'outfile'} eq '-' ) {\n				write_file( $params{'o\
utfile'}, $result );\n			}\n			else {\n				write_f\
ile(\n					$params{'outfile'} . '.'\n					  . $res\
ultType->{'identifier'} . '.'\n					  . $resultTyp\
e->{'fileSuffix'},\n					$result\n				);\n			}\n		\
}\n	}\n	print_debug_message( 'get_results', 'End',\
 1 );\n}\n\n=head2 read_file()\n\nRead a file into\
 a scalar. The special filename '-' can be used to\
 read from \nstandard input (STDIN).\n\n  my $data\
 = &read_file($filename);\n\n=cut\n\nsub read_file\
 {\n	print_debug_message( 'read_file', 'Begin', 1 \
);\n	my $filename = shift;\n	print_debug_message( \
'read_file', 'filename: ' . $filename, 2 );\n	my (\
 $content, $buffer );\n	if ( $filename eq '-' ) {\\
n		while ( sysread( STDIN, $buffer, 1024 ) ) {\n		\
	$content .= $buffer;\n		}\n	}\n	else {    # File\\
n		open( my $FILE, '<', $filename )\n		  or die \"\
Error: unable to open input file $filename ($!)\";\
\n		while ( sysread( $FILE, $buffer, 1024 ) ) {\n	\
		$content .= $buffer;\n		}\n		close($FILE);\n	}\n\
	print_debug_message( 'read_file', 'End', 1 );\n	r\
eturn $content;\n}\n\n=head2 write_file()\n\nWrite\
 data to a file. The special filename '-' can be u\
sed to write to \nstandard output (STDOUT).\n\n  &\
write_file($filename, $data);\n\n=cut\n\nsub write\
_file {\n	print_debug_message( 'write_file', 'Begi\
n', 1 );\n	my ( $filename, $data ) = @_;\n	print_d\
ebug_message( 'write_file', 'filename: ' . $filena\
me, 2 );\n	if ( $outputLevel > 0 ) {\n		print STDE\
RR 'Creating result file: ' . $filename . \"\\n\";\
\n	}\n	if ( $filename eq '-' ) {\n		print STDOUT $\
data;\n	}\n	else {\n		open( my $FILE, '>', $filena\
me )\n		  or die \"Error: unable to open output fi\
le $filename ($!)\";\n		syswrite( $FILE, $data );\\
n		close($FILE);\n	}\n	print_debug_message( 'write\
_file', 'End', 1 );\n}\n\n=head2 usage()\n\nPrint \
program usage message.\n\n  &usage();\n\n=cut\n\ns\
ub usage {\n	print STDERR <<EOF\nWU-BLAST\n=======\
=\n   \nRapid sequence database search programs ut\
ilizing the BLAST algorithm\n    \n[Required]\n\n \
 -p, --program      : str  : BLAST program to use,\
 see --paramDetail program\n  -D, --database     :\
 str  : database(s) to search, space separated. Se\
e\n                              --paramDetail dat\
abase\n      --stype        : str  : query sequenc\
e type, see --paramDetail stype\n  seqFile        \
    : file : query sequence (\"-\" for STDIN, \\@f\
ilename for\n                              identif\
ier list file)\n\n[Optional]\n\n  -m, --matrix    \
   : str  : scoring matrix, see --paramDetail matr\
ix\n  -e, --exp          : real : 0<E<= 1000. Stat\
istical significance threshold \n                 \
             for reporting database sequence match\
es.\n  -e, --viewfilter   :      : display the fil\
tered query sequence\n  -f, --filter       : str  \
: filter the query sequence for low complexity \n \
                             regions, see --paramD\
etail filter\n  -A, --align        : int  : pairwi\
se alignment format, see --paramDetail align\n  -s\
, --scores       : int  : number of scores to be r\
eported\n  -b, --alignments   : int  : number of a\
lignments to report\n  -S, --sensitivity  : str  :\
 sensitivity of the search, \n                    \
          see --paramDetail sensitivity\n  -t, --s\
ort	     : str  : sort order for hits, see --param\
Detail sort\n  -T, --stats        : str  : statist\
ical model, see --paramDetail stats\n  -d, --stran\
d       : str  : DNA strand to search with,\n     \
                         see --paramDetail strand\\
n  -c, --topcombon    : str  : consistent sets of \
HSPs\n      --multifasta   :      : treat input as\
 a set of fasta formatted sequences\n\n[General]\n\
\n  -h, --help        :      : prints this help te\
xt\n      --async       :      : forces to make an\
 asynchronous query\n      --email       : str  : \
e-mail address\n      --title       : str  : title\
 for job\n      --status      :      : get job sta\
tus\n      --resultTypes :      : get available re\
sult types for job\n      --polljob     :      : p\
oll for the status of a job\n      --jobid       :\
 str  : jobid that was returned when an asynchrono\
us job \n                             was submitte\
d.\n      --outfile     : str  : file name for res\
ults (default is jobid;\n                         \
    \"-\" for STDOUT)\n      --outformat   : str  \
: result format to retrieve\n      --params      :\
      : list input parameters\n      --paramDetail\
 : str  : display details for input parameter\n   \
   --quiet       :      : decrease output\n      -\
-verbose     :      : increase output\n      --tra\
ce       :      : show SOAP messages being interch\
anged \n   \nSynchronous job:\n\n  The results/err\
ors are returned as soon as the job is finished.\n\
  Usage: $scriptName --email <your\\@email> [optio\
ns...] seqFile\n  Returns: results as an attachmen\
t\n\nAsynchronous job:\n\n  Use this if you want t\
o retrieve the results at a later time. The result\
s \n  are stored for up to 24 hours. 	\n  Usage: $\
scriptName --async --email <your\\@email> [options\
...] seqFile\n  Returns: jobid\n\n  Use the jobid \
to query for the status of the job. If the job is \
finished, \n  it also returns the results/errors.\\
n  Usage: $scriptName --polljob --jobid <jobId> [-\
-outfile string]\n  Returns: string indicating the\
 status of the job and if applicable, results \n  \
as an attachment.\n\nFurther information:\n\n  htt\
p://www.ebi.ac.uk/Tools/webservices/services/sss/w\
u_blast_rest\n  http://www.ebi.ac.uk/Tools/webserv\
ices/tutorials/perl\n\nSupport/Feedback:\n\n  http\
://www.ebi.ac.uk/support/\nEOF\n}\n\n=head1 FEEDBA\
CK/SUPPORT\n\nPlease contact us at L<http://www.eb\
i.ac.uk/support/> if you have any \nfeedback, sugg\
estions or issues with the service or this client.\
\n\n=cut\n","\n\n\nmy $PROBTRESH = 0.3;# base pair\
s below this prob threshold will be ignored\nmy $W\
EIGHT = 100.0; # float!!\nmy $NUCALPH = \"ACGTUNRY\
MKSWHBVD\";\nuse vars qw($NUCALPH $WEIGHT);\n\nmy \
$myname = basename($0);\n\nuse strict;\nuse warnin\
gs;\n\nuse File::Basename;\nuse Getopt::Long;\nuse\
 File::Glob ':glob';\nuse File::Spec;\nuse File::T\
emp qw/ tempfile tempdir /;\n\n\n\n\nsub tcoffeeli\
b_header($;$)\n{\n    my ($nseq, $fd) = @_;\n    i\
f (! defined($fd)) {\n        $fd = *STDOUT;\n    \
}\n    printf $fd \"! TC_LIB_FORMAT_01\\n\";\n    \
printf $fd \"%d\\n\", $nseq;\n}\n\n\nsub tcoffeeli\
b_header_addseq($$;$)\n{\n    my ($id, $seq, $fd) \
= @_;\n    if (! defined($fd)) {\n        $fd = *S\
TDOUT;\n    }\n    printf $fd \"%s %d %s\\n\", $id\
, length($seq), $seq;\n}\n\n\nsub tcoffeelib_comme\
nt($;$)\n{\n    my ($comment, $fd) = @_;\n    if (\
! defined($fd)) {\n        $fd = *STDOUT;\n    }\n\
    printf $fd \"!\" . $comment . \"\\n\";\n}\n\n\\
nsub tcoffeelib_struct($$$;$)\n{\n    my ($nseq, $\
len, $bpm, $fd) = @_;\n\n    if (! defined($fd)) {\
\n        $fd = *STDOUT;\n    }\n\n    # output ba\
sepair indices with fixed weight\n    printf $fd \\
"#%d %d\\n\", $nseq, $nseq;\n    # output basepair\
s (only once) and with unit-offset\n    for (my $i\
=0; $i<$len; $i++) {\n        for (my $j=$i+1; $j<\
$len; $j++) {\n            if (! defined($bpm->[$i\
][$j])) {\n                print STDERR \"ERROR: \\
\$bpm->[$i][$j] undefined\\n\";\n            }\n  \
          if ($bpm->[$i][$j]>0) {\n               \
 print $fd $i+1;\n                print $fd \" \";\
\n                print $fd $j+1;\n               \
 print $fd \" \" . $bpm->[$i][$j] . \"\\n\";\n    \
        }\n        }\n    }\n}\n\n\nsub tcoffeelib\
_footer(;$)\n{\n    my ($fd) = @_;\n    if (! defi\
ned($fd)) {\n        $fd = *STDOUT;\n    }\n    pr\
int $fd \"! SEQ_1_TO_N\\n\";\n}\n\n\n    \nsub plf\
old($$$)\n{    \n    my ($id, $seq, $probtresh) = \
@_;\n    my (@struct);# return\n    my ($templ, $f\
htmp, $fnametmp, $cmd, $ctr, $window_size);\n    o\
ur $ntemp++;\n    \n    $templ = $myname . \".pid-\
\" . $$ .$ntemp .\".XXXXXX\";\n    ($fhtmp, $fname\
tmp) = tempfile($templ, UNLINK => 1); \n    print \
$fhtmp \">$id\\n$seq\\n\";\n\n    # --- init basep\
air array\n    #\n    for (my $i=0; $i<length($seq\
); $i++) {\n        for (my $j=$i+1; $j<length($se\
q); $j++) {\n            $struct[$i][$j]=0;\n     \
   }\n    }\n\n\n    # --- call rnaplfold and drop\
 a readme\n    #\n    $window_size=(length($seq)<7\
0)?length($seq):70;\n    $cmd = \"RNAplfold -W $wi\
ndow_size < $fnametmp >/dev/null\";\n    system($c\
md);\n    \n    if ($? != 0) {\n        printf STD\
ERR \"ERROR: RNAplfold ($cmd) exited with error st\
atus %d\\n\", $? >> 8;\n        return;\n    }\n  \
  #unlink($fnametmp);\n    my $fps = sprintf(\"%s_\
dp.ps\", $id); # check long name\n    \n    if (! \
-s $fps) {\n      {\n\n	$fps = sprintf(\"%s_dp.ps\\
", substr($id,0,12)); # check short name\n 	if (! \
-s $fps)\n	  {\n	    die(\"couldn't find expected \
file $fps\\n\");\n	    return;\n	  }\n      }\n   \
 }\n\n    \n    # --- read base pairs from created\
 postscript\n    #\n    open(FH, $fps);\n    while\
 (my $line = <FH>) {\n        my ($nti, $ntj, $pro\
b);\n        chomp($line);        \n        # line\
: bp bp sqrt-prob ubox\n        my @match = ($line\
 =~ m/^([0-9]+) +([0-9]+) +([0-9\\.]+) +ubox$/);\n\
        if (scalar(@match)) {\n            $nti=$1\
;\n            $ntj=$2;\n            $prob=$3*$3;#\
 prob stored as square root\n\n            if ($pr\
ob>$probtresh) {\n                #printf STDERR \\
"\\$struct[$nti][$ntj] sqrtprob=$3 prob=$prob > $p\
robtresh\\n\";\n                $struct[$nti-1][$n\
tj-1] = $WEIGHT\n            }\n            # stor\
e with zero-offset\n        }\n    }\n    close(FH\
);\n\n    # remove or gzi postscript\n    #\n    u\
nlink($fps);\n    #\n    # or gzip\n    #$cmd = \"\
gzip -qf $fps\";\n    #system($cmd);\n    #if ($? \
!= 0) {\n    #    printf STDERR \"ERROR: gzip ($cm\
d) exited with error status %d\\n\", $? >> 8;\n   \
 #}\n\n    return \\@struct;\n}\n\n\n\n\n\nsub rna\
seqfmt($)\n{\n    my ($seq) = @_;\n    # remove ga\
ps\n    $seq =~ s/-//g;\n    # uppercase RNA\n    \
$seq = uc($seq);\n    # T -> U\n    $seq =~ s/T/U/\
g;\n    # check for invalid charaters\n    $_ = $s\
eq;\n    s/[^$NUCALPH]//g;\n    return $_;\n}\n\n\\
n\n\nsub usage(;$)\n{    \n    my ($errmsg) = @_;\\
n    if ($errmsg) {\n        print STDERR \"ERROR:\
 $errmsg\\n\";\n    }\n    print STDERR << \"EOF\"\
;\n$myname:\n Creates a T-Coffee RNA structure lib\
rary from RNAplfold prediction.\n See FIXME:citati\
on\nUsage:\n $myname -in seq_file -out tcoffee_lib\
\nEOF\n    exit(1);\n}\n\nsub read_fasta_seq \n  {\
\n    my $f=$_[0];\n    my %hseq;\n    my (@seq, @\
com, @name);\n    my ($a, $s,$nseq);\n\n    open (\
F, $f);\n    while (<F>)\n      {\n	$s.=$_;\n     \
 }\n    close (F);\n\n    \n    @name=($s=~/>(\\S*\
).*\\n[^>]*/g);\n    \n    @seq =($s=~/>.*.*\\n([^\
>]*)/g);\n    @com =($s=~/>(\\S*)(.*)\\n([^>]*)/g)\
;\n\n\n    $nseq=$#name+1;\n  \n    for ($a=0; $a<\
$nseq; $a++)\n      {\n	my $n=$name[$a];\n	my $s;\\
n	$hseq{$n}{name}=$n;\n	$s=$seq[$a];$s=~s/\\s//g;\\
n	\n	$hseq{$n}{seq}=$s;\n	$hseq{$n}{com}=$com[$a];\
\n      }\n    return %hseq;\n  }\n\n\n\n\n\n\n\nm\
y $fmsq = \"\";\nmy $flib = \"\";\nmy %OPTS;\nmy %\
seq;\nmy ($id, $nseq, $i);\nmy @nl;\n\nGetOptions(\
\"in=s\" => \\$fmsq, \"out=s\" => \\$flib);\n\nif \
(! -s $fmsq) {\n    usage(\"empty or non-existant \
file \\\"$fmsq\\\"\")\n}\nif (length($flib)==0) {\\
n    usage(\"empty out-filename\")\n}\n\n\n\n\n\n\\
n%seq=read_fasta_seq($fmsq);\n\n\n@nl=keys(%seq);\\
n\n$nseq=$#nl+1;\nopen FD_LIB, \">$flib\" or die \\
"can't open $flib!\";\ntcoffeelib_header($nseq, *F\
D_LIB);\nforeach $id (keys (%seq))\n  {\n    my ($\
seq, $fmtseq);\n    \n    $seq = $seq{$id}{seq};\n\
    \n    $fmtseq = rnaseqfmt($seq);# check here, \
formatting for folding important later\n    if (le\
ngth($seq)!=length($fmtseq)) {\n        print STDE\
RR \"ERROR: invalid sequence $id is not an RNA seq\
uence. read seq is: $seq\\n\";\n        exit\n    \
  }\n   \n    tcoffeelib_header_addseq($id, uc($se\
q), *FD_LIB);\n  }\ntcoffeelib_comment(\"generated\
 by $myname on \" . localtime(), *FD_LIB);\n\n\n\n\
$i=0;\nforeach $id (keys (%seq))\n  {\n    my ($cl\
eanid, $seq, $bpm);\n    $seq=$seq{$id}{seq};\n   \
 $cleanid = $id;\n    $cleanid =~ s,[/ ],_,g;# nee\
ded for rnaplfold\n    $seq = rnaseqfmt($seq);\n  \
  \n    $bpm = plfold($cleanid, rnaseqfmt($seq), $\
PROBTRESH);       \n    \n    tcoffeelib_struct($i\
+1, length($seq), $bpm, *FD_LIB);\n    $i++;\n}\n\\
n\ntcoffeelib_footer(*FD_LIB);\nclose FD_LIB;\nexi\
t (0);\n\n","\n\n\n\n\n$cmd=join ' ', @ARGV;\nif (\
$cmd=~/-infile=(\\S+)/){ $seqfile=$1;}\nif ($cmd=~\
/-outfile=(\\S+)/){ $libfile=$1;}\n\n\n\n%s=read_f\
asta_seq ($seqfile);\n\nopen (F, \">$libfile\");\n\
foreach $name (keys (%s))\n  {\n    my $tclib=\"$n\
ame.RNAplfold_tclib\";\n    print (F \">$name _F_ \
$tclib\\n\");\n    seq2RNAplfold2tclib ($name, $s{\
$name}{seq}, $tclib);\n  }\nclose (F);\nexit (EXIT\
_SUCCESS);\n\nsub seq2RNAplfold2tclib\n  {\n    my\
 ($name, $seq, $tclib)=@_;\n    my ($tmp);\n    $n\
++;\n    $tmp=\"tmp4seq2RNAplfold_tclib.$$.$n.pep\\
";\n    open (RF, \">$tmp\");\n    print (RF \">$n\
ame\\n$seq\\n\");\n    close (RF);\n    \n    syst\
em \"t_coffee -other_pg RNAplfold2tclib.pl -in=$tm\
p -out=$tclib\";\n    \n    unlink ($tmp);\n    re\
turn $tclib;\n  }\n    \n    \nsub read_fasta_seq \
\n  {\n    my $f=@_[0];\n    my %hseq;\n    my (@s\
eq, @com, @name);\n    my ($a, $s,$nseq);\n\n    o\
pen (F, $f);\n    while (<F>)\n      {\n	$s.=$_;\n\
      }\n    close (F);\n\n    \n    @name=($s=~/>\
(\\S*).*\\n[^>]*/g);\n    \n    @seq =($s=~/>.*.*\\
\n([^>]*)/g);\n    @com =($s=~/>\\S*(.*)\\n([^>]*)\
/g);\n\n    \n    $nseq=$#name+1;\n    \n    for (\
$a=0; $a<$nseq; $a++)\n      {\n	my $n=$name[$a];\\
n	$hseq{$n}{name}=$n;\n	$hseq{$n}{seq}=$seq[$a];\n\
	$hseq{$n}{com}=$com[$a];\n      }\n    return %hs\
eq;\n  }\n","use Getopt::Long;\nuse File::Path;\nu\
se Env;\nuse FileHandle;\nuse Cwd;\nuse Sys::Hostn\
ame;\nour $PIDCHILD;\nour $ERROR_DONE;\nour @TMPFI\
LE_LIST;\nour $EXIT_FAILURE=1;\nour $EXIT_SUCCESS=\
0;\n\nour $REFDIR=getcwd;\nour $EXIT_SUCCESS=0;\no\
ur $EXIT_FAILURE=1;\n\nour $PROGRAM=\"tc_generic_m\
ethod.pl\";\nour $CL=$PROGRAM;\n\nour $CLEAN_EXIT_\
STARTED;\nour $debug_lock=$ENV{\"DEBUG_LOCK\"};\no\
ur $LOCKDIR=$ENV{\"LOCKDIR_4_TCOFFEE\"};\nif (!$LO\
CKDIR){$LOCKDIR=getcwd();}\nour $ERRORDIR=$ENV{\"E\
RRORDIR_4_TCOFFEE\"};\nour $ERRORFILE=$ENV{\"ERROR\
FILE_4_TCOFFEE\"};\n&set_lock ($$);\nif (isshellpi\
d(getppid())){lock4tc(getppid(), \"LLOCK\", \"LSET\
\", \"$$\\n\");}\n      \nour $print;\nmy ($fmsq1,\
 $fmsq2, $output, $outfile, $arch, $psv, $hmmtop_h\
ome, $trim, $cov, $sample, $mode, $gor_home, $gor_\
seq, $gor_obs);\n\nGetOptions(\"-in=s\" => \\$fmsq\
1,\"-output=s\" =>\\$output ,\"-out=s\" => \\$outf\
ile, \"-arch=s\" => \\$arch,\"-psv=s\" => \\$psv, \
\"-hmmtop_home=s\", \\$hmmtop_home,\"-trim=s\" =>\\
\$trim ,\"-print=s\" =>\\$print,\"-cov=s\" =>\\$co\
v , \"-sample=s\" =>\\$sample, \"-mode=s\" =>\\$mo\
de, \"-gor_home=s\"=>\\$gor_home, \"-gor_seq=s\"=>\
\\$gor_seq,\"-gor_obs=s\"=>\\$gor_obs);\n\n\nif (!\
$mode){$mode = \"hmmtop\"}\nelsif ($mode eq \"hmmt\
op\"){;}\nelsif ($mode eq \"gor\"){;}\nelse {myexi\
t(flush_error (\"-mode=$mode is unknown\"));}\n\n\\
nour $HOME=$ENV{\"HOME\"};\nour $MCOFFEE=($ENV{\"M\
COFFEE_4_TCOFFEE\"})?$ENV{\"MCOFFEE_4_TCOFFEE\"}:\\
"$HOME/.t_coffee/mcoffee\";\n\nif ($mode eq \"hmmt\
op\")\n  {\n    check_configuration (\"hmmtop\");\\
n    if (-e $arch){$ENV{'HMMTOP_ARCH'}=$arch;}\n  \
  elsif (-e $ENV{HMMTOP_ARCH}){$arch=$ENV{HMMTOP_A\
RCH};}\n    elsif (-e \"$MCOFFEE/hmmtop.arch\"){$a\
rch=$ENV{'HMMTOP_ARCH'}=\"$MCOFFEE/hmmtop.arch\";}\
\n    elsif (-e \"$hmmtop_home/hmmtop.arc\"){$arch\
=$ENV{'HMMTOP_ARCH'}=\"$hmmtop_home/hmmtop.arc\";}\
\n    else {myexit(flush_error ( \"Could not find \
ARCH file for hmmtop\"));}\n    \n    \n    if (-e\
 $psv){$ENV{'HMMTOP_PSV'}=$psv;}\n    elsif (-e $E\
NV{HMMTOP_PSV}){$psv=$ENV{HMMTOP_PSV};}\n    elsif\
 (-e \"$MCOFFEE/hmmtop.psv\"){$psv=$ENV{'HMMTOP_PS\
V'}=\"$MCOFFEE/hmmtop.psv\";}\n    elsif (-e \"$hm\
mtop_home/hmmtop.psv\"){$psv=$ENV{'HMMTOP_PSV'}=\"\
$hmmtop_home/hmmtop.psv\";}\n    else {myexit(flus\
h_error ( \"Could not find PSV file for hmmtop\"))\
;}\n  }\nelsif ($mode eq \"gor\")\n  {\n    our $G\
OR_SEQ;\n    our $GOR_OBS;\n    \n    check_config\
uration (\"gorIV\");\n    if (-e $gor_seq){$GOR_SE\
Q=$gor_seq;}\n    elsif (-e $ENV{GOR_SEQ}){$GOR_SE\
Q=$ENV{GOR_SEQ};}\n    elsif (-e \"$MCOFFEE/New_KS\
.267.seq\"){$GOR_SEQ=\"$MCOFFEE/New_KS.267.seq\";}\
\n    elsif (-e \"$gor_home/New_KS.267.seq\"){$GOR\
_SEQ=\"$gor_home/New_KS.267.seq\";}\n    else {mye\
xit(flush_error ( \"Could not find SEQ file for go\
r\"));}\n\n    if (-e $gor_obs){$GOR_OBS=$gor_obs;\
}\n    elsif (-e $ENV{GOR_OBS}){$GOR_OBS=$ENV{GOR_\
OBS};}\n    elsif (-e \"$MCOFFEE/New_KS.267.obs\")\
{$GOR_OBS=\"$MCOFFEE/New_KS.267.obs\";}\n    elsif\
 (-e \"$gor_home/New_KS.267.obs\"){$GOR_OBS=\"$gor\
_home/New_KS.267.obs\";}\n    else {myexit(flush_e\
rror ( \"Could not find OBS file for gor\"));}\n  \
}\n\n\nif ( ! -e $fmsq1){myexit(flush_error (\"Cou\
ld Not Read Input file $fmsq1\"));}\n\n\nmy $fmsq2\
=vtmpnam();\nmy $fmsq3=vtmpnam();\nmy $tmpfile=vtm\
pnam();\nmy $predfile=vtmpnam();\n\nif ($trim){$tr\
im_action=\" +trim _aln_%%$trim\\_K1 \";}\nif ($co\
v) {$cov_action= \" +sim_filter _aln_c$cov \";}\n&\
safe_system(\"t_coffee -other_pg seq_reformat -in \
$fmsq1 -action +convert 'BOUJXZ-' $cov_action $tri\
m_action -output fasta_aln -out $fmsq2\");\nmy (%p\
red, %seq, %predA);\n\n\n%seq=read_fasta_seq($fmsq\
2);\n%seq=fasta2sample(\\%seq, $sample);\n\nif (1=\
=2 && $mode eq \"hmmtop\" && $output eq \"cons\")\\
n  {\n    fasta2hmmtop_cons($outfile,\\%seq);\n  }\
\nelse\n  {\n    %pred=fasta2pred(\\%seq, $mode);\\
n    %predA=pred2aln (\\%pred, \\%seq);\n    \n   \
 \n    if (!$output || $output eq \"prediction\"){\
output_fasta_seq (\\%predA, $outfile);}\n    elsif\
 ($output eq \"color_html\"){pred2color (\\%pred,\\
\%seq, $outfile);}\n    elsif ($output eq \"cons\"\
){pred2cons($outfile,\\%predA);}\n    else {flush_\
error (\"$output is an unknown output mode\");}\n \
 }\n\nsub fasta2sample\n  {\n    my $SR=shift;\n  \
  my $it=shift;\n    my %S=%$SR;\n    \n    my $se\
q=index2seq_name (\\%S, 1);\n    my $l=length($S{$\
seq}{seq});\n    my @sl=keys(%S);\n    my $nseq=$#\
sl+1;\n    my $index=$nseq;\n  \n    if (!$sample)\
 {return %S;}\n    for (my $a=0; $a<$it; $a++)\n  \
    {\n	my $newseq=\"\";\n	my $nname=\"$seq\\_samp\
led_$index\";\n	for (my $p=0; $p<$l; $p++)\n	  {\n\
	    my $i=int(rand($nseq));\n	    \n	    my $name\
 = $sl[$i];\n	    my $seq=$S{$name}{seq};\n	    my\
 $r=substr ($seq, $p, 1);\n	    $newseq.=$r;\n	  }\
\n	$S{$nname}{name}=$nname;\n	$S{$nname}{seq}=$new\
seq;\n	$S{$nname}{com}=\"sampled\";\n	$S{$nname}{i\
ndex}=++$index;\n      }\n    return %S;\n  }\n	  \
    \nsub fasta2pred\n  {\n    my $s=shift;\n    m\
y $mode=shift;\n\n    if ( $mode eq \"hmmtop\"){re\
turn fasta2hmmtop_pred($s);}\n    elsif ($mode eq \
\"gor\"){return fasta2gor_pred ($s);}\n  }\nsub fa\
sta2hmmtop_cons\n  {\n    my $outfile=shift;\n    \
my $SR=shift;\n    \n    my $o = new FileHandle;\n\
    my $i = new FileHandle;\n    my $tmp_in =vtmpn\
am();\n    my $tmp_out=vtmpnam();\n    my %seq=%$S\
R;\n    my %pred;\n    my $N=keys(%seq);\n    \n  \
  output_fasta_seq (\\%seq,$tmp_in, \"seq\");\n   \
 `hmmtop -pi=mpred -if=$tmp_in -sf=FAS -pl 2>/dev/\
null >$tmp_out`;\n    open ($o, \">$outfile\");\n \
   open ($i, \"$tmp_out\");\n    while (<$i>)\n   \
   {\n	my $l=$_;\n	if (($l=~/>HP\\:\\s+(\\d+)\\s+(\
.*)/)){my $line=\">$2 NSEQ: $N\\n\";print $o \"$li\
ne\";}\n	elsif ( ($l=~/.*pred(.*)/))  {my $line=\"\
$1\\n\";print $o \"$line\";}\n      }\n    close (\
$o);\n    close ($i);\n    return read_fasta_seq($\
tmp);\n  }\nsub fasta2hmmtop_pred\n  {\n    my $SR\
=shift;\n    my $o = new FileHandle;\n    my $i = \
new FileHandle;\n    my $tmp    =vtmpnam();\n    m\
y $tmp_in =vtmpnam();\n    my $tmp_out=vtmpnam();\\
n    my %seq=%$SR;\n    my %pred;\n    \n\n    out\
put_fasta_seq (\\%seq,$tmp_in, \"seq\");\n    `hmm\
top -if=$tmp_in -sf=FAS -pl 2>/dev/null >$tmp_out`\
;\n    open ($o, \">$tmp\");\n    open ($i, \"$tmp\
_out\");\n    while (<$i>)\n      {\n	my $l=$_;\n	\
if (($l=~/>HP\\:\\s+(\\d+)\\s+(.*)/)){my $line=\">\
$2\\n\";print $o \"$line\";}\n	elsif ( ($l=~/.*pre\
d(.*)/))  {my $line=\"$1\\n\";print $o \"$line\";}\
\n      }\n    close ($o);\n    close ($i);\n    r\
eturn read_fasta_seq($tmp);\n  }\n    \n	\n	\n	   \
 \n	\n	\n\n	\nsub fasta2gor_pred\n  {\n    my $SR=\
shift;\n    my $o = new FileHandle;\n    my $i = n\
ew FileHandle;\n    my $tmp    =vtmpnam();\n    my\
 $tmp_in =vtmpnam();\n    my $tmp_out=vtmpnam();\n\
    my %seq=%$SR;\n    my %pred;\n    \n\n    outp\
ut_fasta_seq (\\%seq,$tmp_in, \"seq\");\n    `gorI\
V -prd $tmp_in -seq $GOR_SEQ -obs $GOR_OBS >$tmp_o\
ut`;\n    open ($o, \">$tmp\");\n    open ($i, \"$\
tmp_out\");\n    while (<$i>)\n      {\n	my $l=$_;\
\n\n	\n	if ( $l=~/>/){print $o \"$l\";}\n	elsif ( \
$l=~/Predicted Sec. Struct./){$l=~s/Predicted Sec.\
 Struct\\.//;print $o \"$l\";}\n      }\n    close\
 ($o);\n    close ($i);\n    return read_fasta_seq\
($tmp);\n  }\n			\n			     \nsub index2seq_name\n \
 {\n    \n    my $SR=shift;\n    my $index=shift;\\
n    \n    \n    my %S=%$SR;\n    \n    foreach my\
 $s (%S)\n      {\n	if ( $S{$s}{index}==$index){re\
turn $s;}\n      }\n    return \"\";\n  }\n\nsub p\
red2cons\n  {\n    my $outfile=shift;\n    my $pre\
dR=shift;\n    my $seq=shift;\n    my %P=%$predR;\\
n    my %C;\n    my ($s,@r,$nseq);\n    my $f= new\
 FileHandle;\n\n    open ($f, \">$outfile\");\n\n \
   if (!$seq){$seq=index2seq_name(\\%P,1);}\n    f\
oreach my $s (keys(%P))\n      {\n	$nseq++;\n	$str\
ing= $P{$s}{seq};\n	$string = uc $string;\n	my @r=\
split (//,$string);\n	for (my $a=0; $a<=$#r; $a++)\
\n	  {\n	    if (($r[$a]=~/[OHICE]/)){$C{$a}{$r[$a\
]}++;}\n	  }\n      }\n    @l=keys(%C);\n    \n   \
 \n    $s=$P{$seq}{seq};\n    print $f \">$seq pre\
d based on $nseq\\n\";\n    @r=split (//,$s);\n   \
 \n    for (my $x=0; $x<=$#r; $x++)\n      {\n	if \
($r[$x] ne \"-\")\n	  {\n	    my $h=$C{$x}{H};\n	 \
   my $i=$C{$x}{I};\n	    my $o=$C{$x}{O};\n	    m\
y $c=$C{$x}{C};\n	    my $e=$C{$x}{E};\n	    my $l\
=$i+$o;\n	    \n	    if ($h>=$i && $h>=$o && $h>=$\
c && $h>=$e){$r[$x]='H';}\n	    elsif ($i>=$o && $\
i>=$c && $i>=$e){$r[$x]='I';}\n	    elsif ($o>=$c \
&& $o>=$e){$r[$x]='O';}\n	    elsif ($c>=$e){$r[$x\
]='C';}\n	    else {$r[$x]='E';}\n	  }\n      }\n \
   $j=join ('', @r);\n    print $f \"$j\\n\";\n   \
 close ($f);\n    return $j;\n  }\n\nsub pred2aln\\
n  {\n    my $PR=shift;\n    my $AR=shift;\n    \n\
    my $f=new FileHandle;\n    my %P=%$PR;\n    my\
 %A=%$AR;\n    my %PA;\n    my $tmp=vtmpnam();\n  \
  my $f= new FileHandle;\n    \n    open ($f, \">$\
tmp\");\n    foreach my $s (sort{$A{$a}{index}<=>$\
A{$b}{index}}(keys (%A)))\n      {\n	my (@list, $s\
eq, @plist, @pseq, $L, $PL, $c, $w);\n	my $seq;\n	\
my $seq=$A{$s}{seq};\n	my $pred=$P{$s}{seq};\n	$se\
q=pred2alnS($P{$s}{seq},$A{$s}{seq});\n	print $f \\
">$s\\n$seq\\n\";\n      }\n    close ($f);\n    r\
eturn read_fasta_seq ($tmp);\n  }\nsub pred2alnS\n\
  {\n    my $pred=shift;\n    my $aln= shift;\n   \
 my ($j,$a,$b);\n    my @P=split (//, $pred);\n   \
 my @A=split (//, $aln);\n    for ($a=$b=0;$a<=$#A\
; $a++)\n      {\n	if ($A[$a] ne \"-\"){$A[$a]=$P[\
$b++];}\n      }\n    if ($b!= ($#P+1)){add_warnin\
g (\"Could not thread sequence: $b $#P\");}\n    \\
n    $j= join ('', @A);\n    return $j;\n  }\nsub \
pred2color\n  {\n    my $predP=shift;\n    my $aln\
P=shift;\n    my $out=shift;\n    my $F=new FileHa\
ndle;\n    my $struc=vtmpnam();\n    my $aln=vtmpn\
am();\n    \n\n    output_fasta_seq ($alnP, $aln);\
\n    my %p=%$predP;\n    \n    open ($F, \">$stru\
c\");\n    \n    \n    foreach my $s (keys(%p))\n \
     {\n	\n	print $F \">$s\\n\";\n	my $s=uc($p{$s}\
{seq});\n	\n	$s=~s/[Oo]/0/g;\n	$s=~s/[Ee]/0/g;\n	\\
n	$s=~s/[Ii]/5/g;\n	$s=~s/[Cc]/5/g;\n	\n	$s=~s/[Hh\
]/9/g;\n	\n	print $F \"$s\\n\";\n      }\n    clos\
e ($F);\n    \n    \n    \n    safe_system ( \"t_c\
offee -other_pg seq_reformat -in $aln -struc_in $s\
truc -struc_in_f number_fasta -output color_html -\
out $out\");\n    return;\n  }\n	  \n    \nsub dis\
play_fasta_seq\n  {\n    my $SR=shift;\n    my %S=\
%$SR;\n    \n    foreach my $s (sort{$S{$a}{index}\
<=>$S{$b}{index}}(keys (%S)))\n      {\n	print STD\
ERR \">$s\\n$S{$s}{seq}\\n\";\n      }\n    close \
($f);\n  }\nsub output_fasta_seq\n  {\n    my $SR=\
shift;\n    my $outfile=shift;\n    my $mode =shif\
t;\n    my $f= new FileHandle;\n    my %S=%$SR;\n \
   \n    \n    open ($f, \">$outfile\");\n    fore\
ach my $s (sort{$S{$a}{index}<=>$S{$b}{index}}(key\
s (%S)))\n      {\n	my $seq=$S{$s}{seq};\n	if ( $m\
ode eq \"seq\"){$seq=~s/\\-//g;}\n	print $f \">$s\\
\n$seq\\n\";\n      }\n    close ($f);\n  }\n     \
 \nsub read_fasta_seq \n  {\n    my $f=$_[0];\n   \
 my %hseq;\n    my (@seq, @com, @name);\n    my ($\
a, $s,$nseq);\n    my $index;\n    open (F, $f);\n\
    while (<F>)\n      {\n	$s.=$_;\n      }\n    c\
lose (F);\n\n    \n    @name=($s=~/>(\\S*).*\\n[^>\
]*/g);\n    \n    @seq =($s=~/>.*.*\\n([^>]*)/g);\\
n    @com =($s=~/>.*(.*)\\n([^>]*)/g);\n\n\n    $n\
seq=$#name+1;\n    \n  \n    for ($a=0; $a<$nseq; \
$a++)\n      {\n	my $n=$name[$a];\n	my $s;\n	$hseq\
{$n}{name}=$n;\n	$s=$seq[$a];$s=~s/\\s//g;\n	$hseq\
{$n}{index}=++$index;\n	$hseq{$n}{seq}=$s;\n	$hseq\
{$n}{com}=$com[$a];\n      }\n    return %hseq;\n \
 }\n\n\nsub file2head\n      {\n	my $file = shift;\
\n	my $size = shift;\n	my $f= new FileHandle;\n	my\
 $line;\n	open ($f,$file);\n	read ($f,$line, $size\
);\n	close ($f);\n	return $line;\n      }\nsub fil\
e2tail\n      {\n	my $file = shift;\n	my $size = s\
hift;\n	my $f= new FileHandle;\n	my $line;\n	\n	op\
en ($f,$file);\n	seek ($f,$size*-1, 2);\n	read ($f\
,$line, $size);\n	close ($f);\n	return $line;\n   \
   }\n\n\nsub vtmpnam\n      {\n	my $r=rand(100000\
);\n	my $f=\"file.$r.$$\";\n	while (-e $f)\n	  {\n\
	    $f=vtmpnam();\n	  }\n	push (@TMPFILE_LIST, $f\
);\n	return $f;\n      }\n\nsub myexit\n  {\n    m\
y $code=@_[0];\n    if ($CLEAN_EXIT_STARTED==1){re\
turn;}\n    else {$CLEAN_EXIT_STARTED=1;}\n    ###\
 ONLY BARE EXIT\n    exit ($code);\n  }\nsub set_e\
rror_lock\n    {\n      my $name = shift;\n      m\
y $pid=$$;\n\n      \n      &lock4tc ($$,\"LERROR\\
", \"LSET\", \"$$ -- ERROR: $name $PROGRAM\\n\");\\
n      return;\n    }\nsub set_lock\n  {\n    my $\
pid=shift;\n    my $msg= shift;\n    my $p=getppid\
();\n    &lock4tc ($pid,\"LLOCK\",\"LRESET\",\"$p$\
msg\\n\");\n  }\nsub unset_lock\n   {\n     \n    \
my $pid=shift;\n    &lock4tc ($pid,\"LLOCK\",\"LRE\
LEASE\",\"\");\n  }\nsub shift_lock\n  {\n    my $\
from=shift;\n    my $to=shift;\n    my $from_type=\
shift;\n    my $to_type=shift;\n    my $action=shi\
ft;\n    my $msg;\n    \n    if (!&lock4tc($from, \
$from_type, \"LCHECK\", \"\")){return 0;}\n    $ms\
g=&lock4tc ($from, $from_type, \"LREAD\", \"\");\n\
    &lock4tc ($from, $from_type,\"LRELEASE\", $msg\
);\n    &lock4tc ($to, $to_type, $action, $msg);\n\
    return;\n  }\nsub isshellpid\n  {\n    my $p=s\
hift;\n    if (!lock4tc ($p, \"LLOCK\", \"LCHECK\"\
)){return 0;}\n    else\n      {\n	my $c=lock4tc($\
p, \"LLOCK\", \"LREAD\");\n	if ( $c=~/-SHELL-/){re\
turn 1;}\n      }\n    return 0;\n  }\nsub isrootp\
id\n  {\n    if(lock4tc (getppid(), \"LLOCK\", \"L\
CHECK\")){return 0;}\n    else {return 1;}\n  }\ns\
ub lock4tc\n	{\n	  my ($pid,$type,$action,$value)=\
@_;\n	  my $fname;\n	  my $host=hostname;\n	  \n	 \
 if ($type eq \"LLOCK\"){$fname=\"$LOCKDIR/.$pid.$\
host.lock4tcoffee\";}\n	  elsif ( $type eq \"LERRO\
R\"){ $fname=\"$LOCKDIR/.$pid.$host.error4tcoffee\\
";}\n	  elsif ( $type eq \"LWARNING\"){ $fname=\"$\
LOCKDIR/.$pid.$host.warning4tcoffee\";}\n	  \n	  i\
f ($debug_lock)\n	    {\n	      print STDERR \"\\n\
\\t---lock4tc(tcg): $action => $fname =>$value (RD\
: $LOCKDIR)\\n\";\n	    }\n\n	  if    ($action eq \
\"LCHECK\") {return -e $fname;}\n	  elsif ($action\
 eq \"LREAD\"){return file2string($fname);}\n	  el\
sif ($action eq \"LSET\") {return string2file ($va\
lue, $fname, \">>\");}\n	  elsif ($action eq \"LRE\
SET\") {return string2file ($value, $fname, \">\")\
;}\n	  elsif ($action eq \"LRELEASE\") \n	    {\n	\
      if ( $debug_lock)\n		{\n		  my $g=new FileHa\
ndle;\n		  open ($g, \">>$fname\");\n		  print $g \
\"\\nDestroyed by $$\\n\";\n		  close ($g);\n		  s\
afe_system (\"mv $fname $fname.old\");\n		}\n	    \
  else\n		{\n		  unlink ($fname);\n		}\n	    }\n	 \
 return \"\";\n	}\n	\nsub file2string\n	{\n	  my $\
file=@_[0];\n	  my $f=new FileHandle;\n	  my $r;\n\
	  open ($f, \"$file\");\n	  while (<$f>){$r.=$_;}\
\n	  close ($f);\n	  return $r;\n	}\nsub string2fi\
le \n    {\n    my ($s,$file,$mode)=@_;\n    my $f\
=new FileHandle;\n    \n    open ($f, \"$mode$file\
\");\n    print $f  \"$s\";\n    close ($f);\n  }\\
n\nBEGIN\n    {\n      srand;\n    \n      $SIG{'S\
IGUP'}='signal_cleanup';\n      $SIG{'SIGINT'}='si\
gnal_cleanup';\n      $SIG{'SIGQUIT'}='signal_clea\
nup';\n      $SIG{'SIGILL'}='signal_cleanup';\n   \
   $SIG{'SIGTRAP'}='signal_cleanup';\n      $SIG{'\
SIGABRT'}='signal_cleanup';\n      $SIG{'SIGEMT'}=\
'signal_cleanup';\n      $SIG{'SIGFPE'}='signal_cl\
eanup';\n      \n      $SIG{'SIGKILL'}='signal_cle\
anup';\n      $SIG{'SIGPIPE'}='signal_cleanup';\n \
     $SIG{'SIGSTOP'}='signal_cleanup';\n      $SIG\
{'SIGTTIN'}='signal_cleanup';\n      $SIG{'SIGXFSZ\
'}='signal_cleanup';\n      $SIG{'SIGINFO'}='signa\
l_cleanup';\n      \n      $SIG{'SIGBUS'}='signal_\
cleanup';\n      $SIG{'SIGALRM'}='signal_cleanup';\
\n      $SIG{'SIGTSTP'}='signal_cleanup';\n      $\
SIG{'SIGTTOU'}='signal_cleanup';\n      $SIG{'SIGV\
TALRM'}='signal_cleanup';\n      $SIG{'SIGUSR1'}='\
signal_cleanup';\n\n\n      $SIG{'SIGSEGV'}='signa\
l_cleanup';\n      $SIG{'SIGTERM'}='signal_cleanup\
';\n      $SIG{'SIGCONT'}='signal_cleanup';\n     \
 $SIG{'SIGIO'}='signal_cleanup';\n      $SIG{'SIGP\
ROF'}='signal_cleanup';\n      $SIG{'SIGUSR2'}='si\
gnal_cleanup';\n\n      $SIG{'SIGSYS'}='signal_cle\
anup';\n      $SIG{'SIGURG'}='signal_cleanup';\n  \
    $SIG{'SIGCHLD'}='signal_cleanup';\n      $SIG{\
'SIGXCPU'}='signal_cleanup';\n      $SIG{'SIGWINCH\
'}='signal_cleanup';\n      \n      $SIG{'INT'}='s\
ignal_cleanup';\n      $SIG{'TERM'}='signal_cleanu\
p';\n      $SIG{'KILL'}='signal_cleanup';\n      $\
SIG{'QUIT'}='signal_cleanup';\n      \n      our $\
debug_lock=$ENV{\"DEBUG_LOCK\"};\n      \n      \n\
      \n      \n      foreach my $a (@ARGV){$CL.=\\
" $a\";}\n      if ( $debug_lock ){print STDERR \"\
\\n\\n\\n********** START PG: $PROGRAM ***********\
**\\n\";}\n      if ( $debug_lock ){print STDERR \\
"\\n\\n\\n**********(tcg) LOCKDIR: $LOCKDIR $$ ***\
**********\\n\";}\n      if ( $debug_lock ){print \
STDERR \"\\n --- $$ -- $CL\\n\";}\n      \n	     \\
n      \n      \n    }\nsub flush_error\n  {\n    \
my $msg=shift;\n    return add_error ($EXIT_FAILUR\
E,$$, $$,getppid(), $msg, $CL);\n  }\nsub add_erro\
r \n  {\n    my $code=shift;\n    my $rpid=shift;\\
n    my $pid=shift;\n    my $ppid=shift;\n    my $\
type=shift;\n    my $com=shift;\n    \n    $ERROR_\
DONE=1;\n    lock4tc ($rpid, \"LERROR\",\"LSET\",\\
"$pid -- ERROR: $type\\n\");\n    lock4tc ($$, \"L\
ERROR\",\"LSET\", \"$pid -- COM: $com\\n\");\n    \
lock4tc ($$, \"LERROR\",\"LSET\", \"$pid -- STACK:\
 $ppid -> $pid\\n\");\n   \n    return $code;\n  }\
\nsub add_warning \n  {\n    my $rpid=shift;\n    \
my $pid =shift;\n    my $command=shift;\n    my $m\
sg=\"$$ -- WARNING: $command\\n\";\n    print STDE\
RR \"$msg\";\n    lock4tc ($$, \"LWARNING\", \"LSE\
T\", $msg);\n  }\n\nsub signal_cleanup\n  {\n    p\
rint dtderr \"\\n**** $$ (tcg) was killed\\n\";\n \
   &cleanup;\n    exit ($EXIT_FAILURE);\n  }\nsub \
clean_dir\n  {\n    my $dir=@_[0];\n    if ( !-d $\
dir){return ;}\n    elsif (!($dir=~/tmp/)){return \
;}#safety check 1\n    elsif (($dir=~/\\*/)){retur\
n ;}#safety check 2\n    else\n      {\n	`rm -rf $\
dir`;\n      }\n    return;\n  }\nsub cleanup\n  {\
\n    #print stderr \"\\n----tc: $$ Kills $PIDCHIL\
D\\n\";\n    #kill (SIGTERM,$PIDCHILD);\n    my $p\
=getppid();\n    $CLEAN_EXIT_STARTED=1;\n    \n   \
 \n    \n    if (&lock4tc($$,\"LERROR\", \"LCHECK\\
", \"\"))\n      {\n	my $ppid=getppid();\n	if (!$E\
RROR_DONE) \n	  {\n	    &lock4tc($$,\"LERROR\", \"\
LSET\", \"$$ -- STACK: $p -> $$\\n\");\n	    &lock\
4tc($$,\"LERROR\", \"LSET\", \"$$ -- COM: $CL\\n\"\
);\n	  }\n      }\n    my $warning=&lock4tc($$, \"\
LWARNING\", \"LREAD\", \"\");\n    my $error=&lock\
4tc($$,  \"LERROR\", \"LREAD\", \"\");\n    #relea\
se error and warning lock if root\n    \n    if (i\
srootpid() && ($warning || $error) )\n      {\n	\n\
	print STDERR \"**************** Summary *********\
****\\n$error\\n$warning\\n\";\n\n	&lock4tc($$,\"L\
ERROR\",\"RELEASE\",\"\");\n	&lock4tc($$,\"LWARNIN\
G\",\"RELEASE\",\"\");\n      } \n    \n    \n    \
foreach my $f (@TMPFILE_LIST)\n      {\n	if (-e $f\
){unlink ($f);} \n      }\n    foreach my $d (@TMP\
DIR_LIST)\n      {\n	clean_dir ($d);\n      }\n   \
 #No More Lock Release\n    #&lock4tc($$,\"LLOCK\"\
,\"LRELEASE\",\"\"); #release lock \n\n    if ( $d\
ebug_lock ){print STDERR \"\\n\\n\\n********** END\
 PG: $PROGRAM ($$) *************\\n\";}\n    if ( \
$debug_lock ){print STDERR \"\\n\\n\\n**********(t\
cg) LOCKDIR: $LOCKDIR $$ *************\\n\";}\n  }\
\nEND \n  {\n    \n    &cleanup();\n  }\n   \n\nsu\
b safe_system \n{\n  my $com=shift;\n  my $ntry=sh\
ift;\n  my $ctry=shift;\n  my $pid;\n  my $status;\
\n  my $ppid=getppid();\n  if ($com eq \"\"){retur\
n 1;}\n  \n  \n\n  if (($pid = fork ()) < 0){retur\
n (-1);}\n  if ($pid == 0)\n    {\n      set_lock(\
$$, \" -SHELL- $com (tcg)\");\n      exec ($com);\\
n    }\n  else\n    {\n      lock4tc ($$, \"LLOCK\\
", \"LSET\", \"$pid\\n\");#update parent\n      $P\
IDCHILD=$pid;\n    }\n  if ($debug_lock){printf ST\
DERR \"\\n\\t .... safe_system (fasta_seq2hmm)  p:\
 $$ c: $pid COM: $com\\n\";}\n\n  waitpid ($pid,WT\
ERMSIG);\n\n  shift_lock ($pid,$$, \"LWARNING\",\"\
LWARNING\", \"LSET\");\n\n  if ($? == $EXIT_FAILUR\
E || lock4tc($pid, \"LERROR\", \"LCHECK\", \"\"))\\
n    {\n      if ($ntry && $ctry <$ntry)\n	{\n	  a\
dd_warning ($$,$$,\"$com failed [retry: $ctry]\");\
\n	  lock4tc ($pid, \"LRELEASE\", \"LERROR\", \"\"\
);\n	  return safe_system ($com, $ntry, ++$ctry);\\
n	}\n      elsif ($ntry == -1)\n	{\n	  if (!shift_\
lock ($pid, $$, \"LERROR\", \"LWARNING\", \"LSET\"\
))\n	    {\n	      add_warning ($$,$$,\"$com faile\
d\");\n	    }\n	  else\n	    {\n	      lock4tc ($p\
id, \"LRELEASE\", \"LERROR\", \"\");\n	    }\n	  r\
eturn $?;}\n      else\n	{\n	  if (!shift_lock ($p\
id,$$, \"LERROR\",\"LERROR\", \"LSET\"))\n	    {\n\
	      myexit(add_error ($EXIT_FAILURE,$$,$pid,get\
ppid(), \"UNSPECIFIED system\", $com));\n	    }\n	\
}\n    }\n  return $?;\n}\n\nsub check_configurati\
on \n    {\n      my @l=@_;\n      my $v;\n      f\
oreach my $p (@l)\n	{\n	  \n	  if   ( $p eq \"EMAI\
L\")\n	    { \n	      if ( !($EMAIL=~/@/))\n		{\n	\
	add_warning($$,$$,\"Could Not Use EMAIL\");\n		my\
exit(add_error ($EXIT_FAILURE,$$,$$,getppid(),\"EM\
AIL\",\"$CL\"));\n	      }\n	    }\n	  elsif( $p e\
q \"INTERNET\")\n	    {\n	      if ( !&check_inter\
net_connection())\n		{\n		  myexit(add_error ($EXI\
T_FAILURE,$$,$$,getppid(),\"INTERNET\",\"$CL\"));\\
n		}\n	    }\n	  elsif( $p eq \"wget\")\n	    {\n	\
      if (!&pg_is_installed (\"wget\") && !&pg_is_\
installed (\"curl\"))\n		{\n		  myexit(add_error (\
$EXIT_FAILURE,$$,$$,getppid(),\"PG_NOT_INSTALLED:w\
get\",\"$CL\"));\n		}\n	    }\n	  elsif( !(&pg_is_\
installed ($p)))\n	    {\n	      myexit(add_error \
($EXIT_FAILURE,$$,$$,getppid(),\"PG_NOT_INSTALLED:\
$p\",\"$CL\"));\n	    }\n	}\n      return 1;\n    \
}\nsub pg_is_installed\n  {\n    my @ml=@_;\n    m\
y $r, $p, $m;\n    my $supported=0;\n    \n    my \
$p=shift (@ml);\n    if ($p=~/::/)\n      {\n	if (\
safe_system (\"perl -M$p -e 1\")==$EXIT_SUCCESS){r\
eturn 1;}\n	else {return 0;}\n      }\n    else\n \
     {\n	$r=`which $p 2>/dev/null`;\n	if ($r eq \"\
\"){return 0;}\n	else {return 1;}\n      }\n  }\n\\
n\n\nsub check_internet_connection\n  {\n    my $i\
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
\n  }\n\n\n\n","\n\n\n\n\nmy $FMODEL =\"\"; \nmy $\
TMPDIR = \"/tmp\";\n\n\n\n\nmy $NUCALPH = \"ACGTUN\
RYMKSWHBVD\";\nmy $PRIMNUCALPH = \"ACGTUN\";\nuse \
vars qw($NUCALPH $PRIMNUCALPH $TMPDIR);\n\n\nmy $e\
rrmsg;\nuse vars qw($errmsg);\n\n\n\nuse Getopt::L\
ong;\nuse Cwd;\nuse File::Basename;\nuse File::Tem\
p qw/ tempfile tempdir /;\nuse File::Copy;\nuse Fi\
le::Path;\n\n\n\nsub usage(;$)\n{\n    my ($errmsg\
) = @_;\n    my $myname = basename($0);\n\n    if \
($errmsg) {\n        print STDERR \"ERROR: $errmsg\
\\n\";\n    }\n\n    print STDERR << \"EOF\";\n   \
 \n$myname: align two sequences by means of consan\
\\'s sfold\nUsage:\n $myname -i file -o file -d pa\
th\nOptions:\n -i|--in : pairwise input sequence f\
ile\n -o|--out: output alignment\n -d|--directory \
containing data\n\nEOF\n}\n\nsub read_stk_aln \n  \
{\n    my $f=$_[0];\n    my ($seq, $id);\n    \n  \
  my %hseq;\n\n    open (STK, \"$f\");\n    while \
(<STK>)\n      {\n	if ( /^#/ || /^\\/\\// || /^\\s\
*$/){;}\n	else\n	  {\n	    ($id,$seq)=/(\\S+)\\s+(\
\\S+)/;\n	    $hseq{$id}{'seq'}.=$seq;\n	  }\n    \
  }\n    close (STK);\n    return %hseq;\n  }\nsub\
 read_fasta_seq \n  {\n    my $f=$_[0];\n    my %h\
seq;\n    my (@seq, @com, @name);\n    my ($a, $s,\
$nseq);\n\n    open (F, $f);\n    while (<F>)\n   \
   {\n	$s.=$_;\n      }\n    close (F);\n\n    \n \
   @name=($s=~/>(.*).*\\n[^>]*/g);\n    \n    @seq\
 =($s=~/>.*.*\\n([^>]*)/g);\n    @com =($s=~/>.*(.\
*)\\n([^>]*)/g);\n\n    \n    $nseq=$#name+1;\n   \
 \n    for ($a=0; $a<$nseq; $a++)\n      {\n	my $n\
=$name[$a];\n	$hseq{$n}{name}=$n;\n	$hseq{$n}{seq}\
=$seq[$a];\n	$hseq{$n}{com}=$com[$a];\n      }\n  \
  return %hseq;\n  }\n\n\n\nsub sfold_parseoutput(\
$$)\n{\n    my ($frawout, $foutfa) = @_;\n    my %\
haln;\n    my ($fstk, $cmd, $id);\n    open FOUTFA\
, \">$foutfa\";\n    \n    $fstk = $frawout . \".s\
tk\";\n    \n    # first line of raw out contains \
info\n    # remaining stuff is stockholm formatted\
\n    $cmd = \"sed -e '1d' $frawout\";\n    system\
(\"$cmd > $fstk\");\n    if ($? != 0) {\n        $\
errmsg = \"command failed with exit status $?.\";\\
n        $errmsg .=  \"Command was \\\"$cmd\\\"\";\
\n        return -1;\n    }\n\n    # this gives an\
 error message. just ignore it...\n    %haln=read_\
stk_aln ( $fstk);\n    foreach $i (keys (%haln))\n\
      {\n	my $s;\n	$s=$haln{$i}{'seq'};\n	$s =~ s/\
\\./-/g;\n	print FOUTFA \">$i\\n$s\\n\";\n      }\\
n    close FOUTFA;\n    return 0;\n}\n\n\n\n\nsub \
sfold_wrapper($$$$)\n{\n    \n    my ($fs1, $fs2, \
$fmodel, $foutfa) = @_;\n    \n\n    my ($cmd, $fr\
awout, $ferrlog, $freadme, $ftimelog, $fstk);\n\n \
   # add  basename($fmsqin) (unknown here!)\n    $\
frawout = \"sfold.log\";\n    $ferrlog = \"sfold.e\
rr\";\n    $ftimelog = \"sfold.time\";\n    $fread\
me =  \"sfold.README\";\n    $fstk = \"sfold.stk\"\
;\n    \n    # prepare execution...\n    #\n    # \
./tmp is essential for dswpalign\n    # otherwise \
you'll get a segfault\n    mkdir \"./tmp\";\n    \\
n    $cmd = \"sfold -m $fmodel $fs1 $fs2\";\n    o\
pen(FREADME,\">$freadme\");\n    print FREADME \"$\
cmd\\n\"; \n    close(FREADME);\n\n    # and go\n \
   #\n    system(\"/usr/bin/time -p -o $ftimelog $\
cmd >$frawout 2>$ferrlog\");\n    if ($? != 0) {\n\
        $errmsg = \"command failed with exit statu\
s $?\";\n        $errmsg .= \"command was \\\"$cmd\
\\\". See \" . getcwd . \"\\n\";\n        return -\
1;\n    }\n\n    return sfold_parseoutput($frawout\
, $foutfa);\n}\n\n\n\n\n\n\n\nmy ($help, $fmsqin, \
$fmsaout);\nGetOptions(\"help\"  => \\$help,\n    \
       \"in=s\" => \\$fmsqin,\n           \"out=s\\
" => \\$fmsaout,\n	   \"data=s\" => \\$ref_dir);\n\
\n\n\nif ($help) {\n    usage();\n    exit(0);\n}\\
nif (! defined($fmsqin)) {\n    usage('missing inp\
ut filename');\n    exit(1);\n}\nif (! defined($fm\
saout)) {\n    usage('missing output filename');\n\
    exit(1);\n\n}\nif (scalar(@ARGV)) {\n    usage\
('Unknown remaining args');\n    exit(1);\n}\n\n$F\
MODEL = \"$ref_dir/mix80.mod\";\nif (! -e \"$FMODE\
L\") {\n    die(\"couldn't find sfold grammar mode\
l file. Expected $FMODEL\\n\");\n}\n\n\nmy %hseq=r\
ead_fasta_seq ($fmsqin);\nmy $id;\n\nforeach $id (\
keys(%hseq))\n  {\n    push(@seq_array, $hseq{$id}\
);\n  }\n\nif ( scalar(@seq_array) != 2 ) {\n    d\
ie(\"Need *exactly* two sequences as input (pairwi\
se alignment!).\")\n}\n\n\n\nmy ($sec, $min, $hour\
, $mday, $mon, $year, $wday, $yday, $isdst) = loca\
ltime(time);\nmy $datei = sprintf(\"%4d-%02d-%02d\\
", $year+1900, $mon+1, $mday);\nmy $templ = basena\
me($0) . \".\" . $datei . \".pid-\" . $$ . \".XXXX\
XX\";\nmy $wd = tempdir ( $templ, DIR => $TMPDIR);\
\n\ncopy($fmsqin, \"$wd/\" . basename($fmsqin) . \\
".org\"); # for reproduction\ncopy($FMODEL, \"$wd\\
");\nmy $fmodel = basename($FMODEL);\nmy $orgwd = \
getcwd;\nchdir $wd;\n\n\n\nmy @sepseqfiles;\nforea\
ch $id (keys(%hseq)) {\n    my ($seq, $orgseq, $fn\
ame, $sout);\n    $seq=$hseq{$id}{'seq'};\n    \n \
   $fname = basename($fmsqin) . \"_$id.fa\";\n    \
# replace funnies in file/id name (e.g. \"/\" \" \\
" etc)\n    $fname =~ s,[/ ],_,g;\n    open (PF, \\
">$fname\");\n    print (PF \">$id\\n$seq\\n\");\n\
    close (PF);\n\n    push(@sepseqfiles, $fname);\
\n}\n\nmy ($f1, $f2, $fout);\n$f1 = $sepseqfiles[0\
];\n$f2 = $sepseqfiles[1];\n$fout = $wd . basename\
($fmsqin) . \".out.fa\";\nif (sfold_wrapper($f1, $\
f2, $fmodel, \"$fout\") != 0) {\n    printf STDERR\
 \"ERROR: See logs in $wd\\n\";\n    exit(1);\n} e\
lse {\n    chdir $orgwd;\n    copy($fout, $fmsaout\
);\n    rmtree($wd);\n   exit(0);\n}\n","\nuse Env\
 qw(HOST);\nuse Env qw(HOME);\nuse Env qw(USER);\n\
\n\n$tmp=clean_cr ($ARGV[0]);\nopen (F, $tmp);\n\n\
while ( <F>)\n  {\n    my $l=$_;\n    if ( $l=~/^#\
 STOCKHOLM/){$stockholm=1;}\n    elsif ( $stockhol\
m && $l=~/^#/)\n      {\n	$l=~/^#(\\S+)\\s+(\\S+)\\
\s+(\\S*)/g;\n	$l=\"_stockholmhasch_$1\\_stockholm\
space_$2 $3\\n\";\n      }\n    $file.=$l;\n  }\nc\
lose (F);\nunlink($tmp);\n$file1=$file;\n\n$file=~\
s/\\#/_hash_symbol_/g;\n$file=~s/\\@/_arobase_symb\
ol_/g;\n\n\n$file=~s/\\n[\\.:*\\s]+\\n/\\n\\n/g;\n\
\n$file=~s/\\n[ \\t\\r\\f]+(\\b)/\\n\\1/g;\n\n\n$f\
ile=~s/(\\n\\S+)(\\s+)(\\S)/\\1_blank_\\3/g;\n\n$f\
ile=~s/[ ]//g;\n$file=~s/_blank_/ /g;\n\n\n\n$file\
 =~s/\\n\\s*\\n/#/g;\n\n$file.=\"#\";\n$file =~s/\\
\n/@/g;\n\n\n\n\n@blocks=split /\\#/, $file;\nshif\
t (@blocks);\n@s=split /\\@/, $blocks[0];\n$nseq=$\
#s+1;\n\n\n\n$file=join '@', @blocks;\n@lines=spli\
t /\\@/,$file;\n\n$c=0;\n\nforeach $l (@lines)\n  \
{\n    if (!($l=~/\\S/)){next;}\n    elsif ($stock\
holm && ($l=~/^\\/\\// || $l=~/STOCKHOLM/)){next;}\
#get read of STOCHOLM Terminator\n   \n    $l=~/(\\
\S+)\\s+(\\S*)/g;\n    $n=$1; $s=$2;\n    \n    $s\
eq[$c].=$s;\n    $name[$c]=$n;\n    $c++;\n    \n \
   if ( $c==$nseq){$c=0;}\n    \n  } \n\nif ( $c!=\
0)\n      {\n	print STDERR \"ERROR: $ARGV[0] is NO\
T an MSA in Clustalw format: make sure there is no\
 blank line within a block [ERROR]\\n\";\n	exit (E\
XIT_FAILURE);\n      }\n\nfor ($a=0; $a< $nseq; $a\
++)\n  {\n    $name[$a]=cleanstring ($name[$a]);\n\
    $seq[$a]=cleanstring ($seq[$a]);\n    $seq[$a]\
=breakstring($seq[$a], 60);\n    \n    $line=\">$n\
ame[$a]\\n$seq[$a]\\n\";\n    \n    print \"$line\\
";\n  }\nexit (EXIT_SUCCESS);\n\nsub cleanstring\n\
  {\n    my $s=@_[0];\n    $s=~s/_hash_symbol_/\\#\
/g;\n    $s=~s/_arobase_symbol_/\\@/g;\n    $s=~s/\
[ \\t]//g;\n    return $s;\n  }\nsub breakstring\n\
  {\n    my $s=@_[0];\n    my $size=@_[1];\n    my\
 @list;\n    my $n,$ns, $symbol;\n    \n    @list=\
split //,$s;\n    $n=0;$ns=\"\";\n    foreach $sym\
bol (@list)\n      {\n	if ( $n==$size)\n	  {\n	   \
 $ns.=\"\\n\";\n	    $n=0;\n	  }\n	$ns.=$symbol;\n\
	$n++;\n      }\n    return $ns;\n    }\n\nsub cle\
an_cr\n  {\n    my $f=@_[0];\n    my $file;\n    \\
n    $tmp=\"f$.$$\";\n    \n    \n    open (IN, $f\
);\n    open (OUT, \">$tmp\");\n    \n    while ( \
<IN>)\n      {\n	$file=$_;\n	$file=~s/\\r\\n/\\n/g\
;\n	$file=~s/\\n\\r/\\n/g;\n	$file=~s/\\r\\r/\\n/g\
;\n	$file=~s/\\r/\\n/g;\n	print OUT \"$file\";\n  \
    }\n    \n    close (IN);\n    close (OUT);\n  \
  return $tmp;\n  }\n","use Env qw(HOST);\nuse Env\
 qw(HOME);\nuse Env qw(USER);\n\n\n$query_start=-1\
;\n$query_end=-1;\n\nwhile (<>)\n  {\n    if ( /\\\
/\\//){$in_aln=1;}\n    elsif ( $in_aln && /(\\S+)\
\\s+(.*)/)\n      {\n\n\n	$name=$1;\n	\n\n	$seq=$2\
;\n	$seq=~s/\\s//g;\n        $seq=~s/\\~/\\-/g;\n	\
$seq=~s/\\./\\-/g;\n	if ( $list{$n}{'name'} && $li\
st{$n}{'name'} ne $name)\n	  {\n	    print \"$list\
{$n}{'name'} Vs $name\";\n	    \n	    exit (EXIT_F\
AILURE);\n	  }\n	else\n	  {\n	    $list{$n}{'name'\
}= $name;\n	  }\n\n	$list{$n}{'seq'}=$list{$n}{'se\
q'}.$seq;\n	\n	$nseq=++$n;\n	\n      }\n    else\n\
      {$n=0;}\n  }\n\n\nfor ($a=0; $a<$nseq; $a++)\
\n  {\n    print \">$list{$a}{'name'}\\n$list{$a}{\
'seq'}\\n\";\n  }\n      \n","\nuse Env qw(HOST);\\
nuse Env qw(HOME);\nuse Env qw(USER);\n\n         \
                                               \nu\
se strict;                                        \
     \nuse warnings;\nuse diagnostics;\n\nmy $in_h\
it_list, my $in_aln=0, my(%name_list)=(),my (%list\
)=(),my $n_seq=0; my $test=0;\nmy($j)=0, my $n=0, \
my $nom, my $lg_query, my %vu=();\n\nopen (F, \">t\
mp\");\n\n$/=\"\\n\";\nwhile (<>)\n{\n    print F \
$_;\n    if($_ =~ /Query=\\s*(.+?)\\s/i) { $nom=$1\
;}\n\n    if ( /Sequences producing significant al\
ignments/){$in_hit_list=1;}\n    \n    if ($_=~ /^\
pdb\\|/i) { $_=~ s/pdb\\|//g; }\n    if ($_=~ /^(1\
_\\d+)\\s+\\d+/) { $_=~ s/$1/QUERY/;}\n      \n   \
 if ( /^(\\S+).+?\\s+[\\d.]+\\s+([\\de.-]+)\\s+$/ \
&& $in_hit_list)	\n    {\n	my($id)=$1; # \n	$id=~ \
s/\\|/_/g; #\n	if ($id =~ /.+_$/) { chop($id) }; #\
\n	$name_list{$n_seq++}=$id;\n	$name_list{$n_seq-1\
}=~ s/.*\\|//g;     \n    }\n  \n    if (/query/i)\
 {$in_aln=1;}\n    if ( /^(\\S+)\\s+(\\d+)\\s+([a-\
zA-Z-]+)\\s+(\\d+)/ || /^(\\S+)(\\s+)(\\-+)(\\s+)/\
 && ($in_aln == 1))\n    {\n	my $name=$1;\n	my $st\
art=$2;\n	my $seq=$3;\n	my $end=$4;\n		\n	if ($nam\
e =~ /QUERY/i) { $lg_query=length($seq); }\n\n	unl\
ess ($test > $n) #m\n	{\n	    my(@seqq)= split('',\
$seq);\n	    my($gap_missing)= scalar(@seqq);\n	  \
  \n	    while ($gap_missing != $lg_query)  { unsh\
ift (@seqq,\"-\"); $gap_missing= scalar(@seqq); }\\
n	    $seq=join('',@seqq);  #m\n	}\n	\n	if ($name \
=~ /QUERY/i)\n	{\n	    $n=0; %vu=(); $j=0;\n	    $\
list{$n}{'real_name'}=\"$nom\";\n	}	\n	else\n	{\n	\
    unless (exists $vu{$name}) { ++$j;}	\n	    $li\
st{$n}{'real_name'}=$name_list{$j-1};\n	}\n		\n	$l\
ist{$n}{'name'}=$name;\n\n	$seq=~tr/a-z/A-Z/;\n	$l\
ist{$n}{'seq'}=$list{$n}{'seq'};\n	$list{$n}{'seq'\
}.=$seq;\n\n	$n++;\n	$vu{$name}++;\n	$test++;\n   \
} \n    \n}\n\nmy @numero=();\n\nfor (my $a=0; $a<\
$n; $a++) #m\n{\n    my $long=length($list{0}{'seq\
'});  \n    my $long1= length($list{$a}{'seq'});\n\
  \n    while ($long1 ne $long)\n    {\n	$list{$a}\
{'seq'}.=\"-\";\n	$long1= length ($list{$a}{'seq'}\
);\n    } \n \n    push (@numero,\"$list{$a}{'name\
'} $list{$a}{'real_name'}\\n\");\n}\n\nmy %dejavu=\
();\n\n\nfor (my $i=0; $i<=$#numero; $i++)\n{\n   \
 my $s=\">$list{$i}{'real_name'}\\n$list{$i}{'seq'\
}\\n\";\n    my $k=0;\n    \n    if (exists $dejav\
u{$numero[$i]}) {next;}\n    else\n    {	\n	for ($\
j=0; $j<$n ; $j++)\n	{\n	    if (\"$numero[$i]\" e\
q \"$numero[$j]\" && $j != $i )\n	    {\n		++$k;\n\
		$s .=\">$list{$j}{'real_name'}\\n$list{$j}{'seq'\
}\\n\";\n	    }\n	}	\n    }\n    \n    if ($k>0) \\
n    {\n	my $cons;\n	open (SOR,\">tempo_aln2cons\"\
); print SOR $s;  close SOR ;\n	open (COM,\"t_coff\
ee -other_pg seq_reformat -in tempo_aln2cons -acti\
on +aln2cons +upper |\") ; \n     	while (<COM>)\n\
	{	\n	    if (/^>/) { $cons =\">$list{$i}{'real_na\
me'}\\n\"; next;}\n	    $_=~ s/\\n//g;\n	    $cons\
 .=$_;\n	}\n	close COM; unlink (\"tempo_aln2cons\"\
);\n	print $cons,\"\\n\"; print F $cons,\"\\n\";\n\
    }	\n    else  { print $s;  print F $s; }\n    \
\n    $dejavu{$numero[$i]}++;\n} #m\n\nexit;\n\n\n\
\n\n\n\n\n\n\n\n\n","use Env;\n\n\n$tmp_dir=\"\";\\
n$init_dir=\"\";\n$program=\"tc_generic_method.pl\\
";\n\n$blast=@ARGV[0];\n\n$name=\"query\";$seq=\"\\
";\n%p=blast_xml2profile($name,$seq,100, 0, 0, $bl\
ast);\n&output_profile (%p);\n\n\nsub output_profi\
le\n  {\n    my (%profile)=(@_);\n    my ($a);\n  \
  for ($a=0; $a<$profile{n}; $a++)\n      {\n	\n	p\
rint \">$profile{$a}{name} $profile{$a}{comment}\\\
n$profile{$a}{seq}\\n\";\n      }\n    return;\n  \
}\nsub file_contains \n  {\n    my ($file, $tag, $\
max)=(@_);\n    my ($n);\n    $n=0;\n    \n    if \
( !-e $file && ($file =~/$tag/)) {return 1;}\n    \
elsif ( !-e $file){return 0;}\n    else \n      {\\
n	open (FC, \"$file\");\n	while ( <FC>)\n	  {\n	  \
  if ( ($_=~/$tag/))\n	      {\n		close (FC);\n		r\
eturn 1;\n	      }\n	    elsif ($max && $n>$max)\n\
	      {\n		close (FC);\n		return 0;\n	      }\n	 \
   $n++;\n	  }\n      }\n    close (FC);\n    retu\
rn 0;\n  }\n	    \n	  \nsub file2string\n  {\n    \
my $f=@_[0];\n    my $string, $l;\n    open (F,\"$\
f\");\n    while (<F>)\n      {\n\n	$l=$_;\n	#chom\
p ($l);\n	$string.=$l;\n      }\n    close (F);\n \
   $string=~s/\\r\\n//g;\n    $string=~s/\\n//g;\n\
    return $string;\n  }\n\n\n\nsub tag2value \n  \
{\n    \n    my $tag=(@_[0]);\n    my $word=(@_[1]\
);\n    my $return;\n    \n    $tag=~/$word=\"([^\\
"]+)\"/;\n    $return=$1;\n    return $return;\n  \
}\n      \nsub hit_tag2pdbid\n  {\n    my $tag=(@_\
[0]);\n    my $pdbid;\n       \n    $tag=~/id=\"(\\
\S+)\"/;\n    $pdbid=$1;\n    $pdbid=~s/_//;\n    \
return $pdbid;\n  }\nsub id2pdbid \n  {\n    my $i\
d=@_[0];\n  \n    if ($id =~/pdb/)\n      {\n	$id=\
~/pdb(.*)/;\n	$id=$1;\n      }\n    $id=~s/[|�_]//\
g;\n    return $id;\n  }\nsub set_blast_type \n  {\
\n    my $file =@_[0];\n    if (&file_contains ($f\
ile,\"EBIApplicationResult\",100)){$BLAST_TYPE=\"E\
BI\";}\n    elsif (&file_contains ($file,\"NCBI_Bl\
astOutput\",100)) {$BLAST_TYPE=\"NCBI\";}\n    els\
e\n      {\n	$BLAST_TYPE=\"\";\n      }\n    retur\
n $BLAST_TYPE;\n  }\nsub blast_xml2profile \n  {\n\
    my ($name,$seq,$maxid, $minid, $mincov, $file)\
=(@_);\n    my (%p, $a, $string, $n);\n    \n\n\n \
   if ($BLAST_TYPE eq \"EBI\" || &file_contains ($\
file,\"EBIApplicationResult\",100)){%p=ebi_blast_x\
ml2profile(@_);}\n    elsif ($BLAST_TYPE eq \"NCBI\
\" || &file_contains ($file,\"NCBI_BlastOutput\",1\
00)){%p=ncbi_blast_xml2profile(@_);}\n    else \n \
     {\n	print \"************ ERROR: Blast Returne\
d an unknown XML Format **********************\";\\
n	die;\n      }\n    for ($a=0; $a<$p{n}; $a++)\n \
     {\n	my $name=$p{$a}{name};\n	$p{$name}{seq}=$\
p{$a}{seq};\n      }\n    return %p;\n  }\nsub ncb\
i_blast_xml2profile \n  {\n    my ($name,$seq,$max\
id, $minid, $mincov, $string)=(@_);\n    my ($L,$l\
, $a,$b,$c,$d,$nhits,@identifyerL);\n    \n    \n \
   $seq=~s/[^a-zA-Z]//g;\n    $L=length ($seq);\n \
   \n    %hit=&xml2tag_list ($string, \"Hit\");\n \
   \n    \n    for ($nhits=0,$a=0; $a<$hit{n}; $a+\
+)\n      {\n	my ($ldb,$id, $identity, $expectatio\
n, $start, $end, $coverage, $r);\n	my (%ID,%DE,%HS\
P);\n	\n	$ldb=\"\";\n\n	%ID=&xml2tag_list ($hit{$a\
}{body}, \"Hit_id\");\n	$identifyer=$ID{0}{body};\\
n	\n	%DE=&xml2tag_list ($hit{$a}{body}, \"Hit_def\\
");\n	$definition=$DE{0}{body};\n	\n	%HSP=&xml2tag\
_list ($hit{$a}{body}, \"Hsp\");\n	for ($b=0; $b<$\
HSP{n}; $b++)\n	  {\n	    my (%START,%END,%E,%I,%Q\
,%M);\n\n	 \n	    %START=&xml2tag_list ($HSP{$b}{b\
ody}, \"Hsp_query-from\");\n	    %HSTART=&xml2tag_\
list ($HSP{$b}{body}, \"Hsp_hit-from\");\n	    \n	\
    %LEN=  &xml2tag_list ($HSP{$b}{body}, \"Hsp_al\
ign-len\");\n	    %END=  &xml2tag_list ($HSP{$b}{b\
ody}, \"Hsp_query-to\");\n	    %HEND=  &xml2tag_li\
st ($HSP{$b}{body}, \"Hsp_hit-to\");\n	    %E=&xml\
2tag_list     ($HSP{$b}{body}, \"Hsp_evalue\");\n	\
    %I=&xml2tag_list     ($HSP{$b}{body}, \"Hsp_id\
entity\");\n	    %Q=&xml2tag_list     ($HSP{$b}{bo\
dy}, \"Hsp_qseq\");\n	    %M=&xml2tag_list     ($H\
SP{$b}{body}, \"Hsp_hseq\");\n	    \n	    for ($e=\
0; $e<$Q{n}; $e++)\n\n	      {\n		$qs=$Q{$e}{body}\
;\n		$ms=$M{$e}{body};\n		if ($seq eq\"\"){$seq=$q\
s;$L=length($seq);}\n		\n		$expectation=$E{$e}{bod\
y};\n		$identity=($LEN{$e}{body}==0)?0:$I{$e}{body\
}/$LEN{$e}{body}*100;\n		$start=$START{$e}{body};\\
n		$end=$END{$e}{body};\n		$Hstart=$HSTART{$e}{bod\
y};\n		$Hend=$HEND{$e}{body};\n	\n		$coverage=(($e\
nd-$start)*100)/$L;\n\n	\n		if ($identity>$maxid |\
| $identity<$minid || $coverage<$mincov){next;}\n	\
	@lr1=(split (//,$qs));\n		@lr2=(split (//,$ms));\\
n		$l=$#lr1+1;\n		for ($c=0;$c<$L;$c++){$p[$nhits]\
[$c]=\"-\";}\n		for ($d=0,$c=0; $c<$l; $c++)\n		  \
{\n		    $r=$lr1[$c];\n		    if ( $r=~/[A-Za-z]/)\\
n		      {\n			\n			$p[$nhits][$d + $start-1]=$lr2\
[$c];\n			$d++;\n		      }\n		  }\n		$Qseq[$nhits]\
=$qs;\n		$Hseq[$nhits]=$ms;\n		$QstartL[$nhits]=$s\
tart;\n		$HstartL[$nhits]=$Hstart;\n		$identityL[$\
nhits]=$identity;\n		$endL[$nhits]=$end;\n		$defin\
itionL[$nhits]=$definition;\n		$identifyerL[$nhits\
]=$identifyer;\n		$comment[$nhits]=\"$ldb|$identif\
yer [Eval=$expectation][id=$identity%][start=$Hsta\
rt end=$Hend]\";\n		$nhits++;\n	      }\n	  }\n   \
   }\n    \n    $profile{n}=0;\n    $profile{$prof\
ile{n}}{name}=$name;\n    $profile{$profile{n}}{se\
q}=$seq;\n    $profile {n}++;\n    \n    for ($a=0\
; $a<$nhits; $a++)\n      {\n	$n=$a+1;\n	\n	$profi\
le{$n}{name}=\"$name\\_$a\";\n	$profile{$n}{seq}=\\
"\";\n	$profile{$n}{Qseq}=$Qseq[$a];\n	$profile{$n\
}{Hseq}=$Hseq[$a];\n	$profile{$n}{Qstart}=$QstartL\
[$a];\n	$profile{$n}{Hstart}=$HstartL[$a];\n	$prof\
ile{$n}{identity}=$identityL[$a];\n	$profile{$n}{d\
efinition}=$definitionL[$a];\n	$profile{$n}{identi\
fyer}=$identifyerL[$a];\n	$profile{$n}{comment}=$c\
omment[$a];\n	for ($b=0; $b<$L; $b++)\n	  {\n	    \
if ($p[$a][$b])\n	      {\n		$profile{$n}{seq}.=$p\
[$a][$b];\n	      }\n	    else\n	      {\n		$profi\
le{$n}{seq}.=\"-\";\n	      }\n	  }\n      }\n    \
\n    $profile{n}=$nhits+1;\n    return %profile;\\
n  }\nsub ebi_blast_xml2profile \n  {\n    my ($na\
me,$seq,$maxid, $minid, $mincov, $string)=(@_);\n \
   my ($L,$l, $a,$b,$c,$d,$nhits,@identifyerL,$ide\
ntifyer);\n    \n\n    \n    $seq=~s/[^a-zA-Z]//g;\
\n    $L=length ($seq);\n    %hit=&xml2tag_list ($\
string, \"hit\");\n    \n    for ($nhits=0,$a=0; $\
a<$hit{n}; $a++)\n      {\n	my ($ldb,$id, $identit\
y, $expectation, $start, $end, $coverage, $r);\n	m\
y (%Q,%M,%E,%I);\n	\n	$ldb=&tag2value ($hit{$a}{op\
en}, \"database\");\n	$identifyer=&tag2value ($hit\
{$a}{open}, \"id\");\n\n	$description=&tag2value (\
$hit{$a}{open}, \"description\");\n	\n	%Q=&xml2tag\
_list ($hit{$a}{body}, \"querySeq\");\n	%M=&xml2ta\
g_list ($hit{$a}{body}, \"matchSeq\");\n	%E=&xml2t\
ag_list ($hit{$a}{body}, \"expectation\");\n	%I=&x\
ml2tag_list ($hit{$a}{body}, \"identity\");\n	\n\n\
	for ($b=0; $b<$Q{n}; $b++)\n	  {\n	    \n	    \n	\
    $qs=$Q{$b}{body};\n	    $ms=$M{$b}{body};\n	  \
  if ($seq eq\"\"){$seq=$qs;$L=length($seq);}\n\n	\
    $expectation=$E{$b}{body};\n	    $identity=$I{\
$b}{body};\n	    \n	    	    \n	    $start=&tag2va\
lue ($Q{$b}{open}, \"start\");\n	    $end=&tag2val\
ue ($Q{$b}{open}, \"end\");\n	    $startM=&tag2val\
ue ($M{$b}{open}, \"start\");\n	    $endM=&tag2val\
ue ($M{$b}{open}, \"end\");\n	    $coverage=(($end\
-$start)*100)/$L;\n	    \n	   # print \"$id: ID: $\
identity COV: $coverage [$start $end]\\n\";\n	    \
\n	    \n	    if ($identity>$maxid || $identity<$m\
inid || $coverage<$mincov){next;}\n	    # print \"\
KEEP\\n\";\n\n	    \n	    @lr1=(split (//,$qs));\n\
	    @lr2=(split (//,$ms));\n	    $l=$#lr1+1;\n	  \
  for ($c=0;$c<$L;$c++){$p[$nhits][$c]=\"-\";}\n	 \
   for ($d=0,$c=0; $c<$l; $c++)\n	      {\n		$r=$l\
r1[$c];\n		if ( $r=~/[A-Za-z]/)\n		  {\n		    \n		\
    $p[$nhits][$d + $start-1]=$lr2[$c];\n		    $d+\
+;\n		  }\n	      }\n	  \n	    \n	    $identifyerL\
[$nhits]=$identifyer;\n	    $comment[$nhits]=\"$ld\
b|$identifyer [Eval=$expectation][id=$identity%][s\
tart=$startM end=$endM]\";\n	    $nhits++;\n	  }\n\
      }\n    \n    $profile{n}=0;\n    $profile{$p\
rofile{n}}{name}=$name;\n    $profile{$profile{n}}\
{seq}=$seq;\n    $profile {n}++;\n    \n    for ($\
a=0; $a<$nhits; $a++)\n      {\n	$n=$a+1;\n	$profi\
le{$n}{name}=\"$name\\_$a\";\n	$profile{$n}{seq}=\\
"\";\n	$profile{$n}{identifyer}=$identifyerL[$a];\\
n	\n	$profile{$n}{comment}=$comment[$a];\n	for ($b\
=0; $b<$L; $b++)\n	  {\n	    if ($p[$a][$b])\n	   \
   {\n		$profile{$n}{seq}.=$p[$a][$b];\n	      }\n\
	    else\n	      {\n		$profile{$n}{seq}.=\"-\";\n\
	      }\n	  }\n      }\n    $profile{n}=$nhits+1;\
\n    \n    return %profile;\n  }\n\nsub blast_xml\
2hit_list\n  {\n    my $string=(@_[0]);\n    retur\
n &xml2tag_list ($string, \"hit\");\n  }\nsub xml2\
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
    return %tag;\n  }\n\n\n\n\n","use Env qw(HOST)\
;\nuse Env qw(HOME);\nuse Env qw(USER);\nwhile (<>\
)\n  {\n    if ( /^>(\\S+)/)\n      {\n	if ($list{\
$1})\n	  {\n	    print \">$1_$list{$1}\\n\";\n	   \
 $list{$1}++;\n	  }\n	else\n	  {\n	    print $_;\n\
	    $list{$1}=1;\n	  }\n      }\n    else\n      \
{\n	print $_;\n      }\n  }\n      \n","\n\n\nuse \
Env qw(HOST);\nuse Env qw(HOME);\nuse Env qw(USER)\
;\n\n\nopen (F,$ARGV[0]);\nwhile ( <>)\n  {\n    @\
x=/([^:,;\\)\\(\\s]+):[^:,;\\)\\(]*/g;\n    @list=\
(@list,@x);\n  }\n$n=$#list+1;\nforeach $n(@list){\
print \">$n\\nsequence\\n\";}\n\n\nclose (F);\n","\
\nopen (F, $ARGV[0]);\n\nwhile ( <F>)\n  {\n    @l\
=($_=~/(\\S+)/g);\n    \n    $name=shift @l;\n    \
\n    print STDOUT \"\\n>$name\\n\";\n    foreach \
$e (@l){$e=($e eq \"0\")?\"O\":\"I\";print \"$e\";\
}\n  }\nclose (F);\n\n		       \n    \n","use Env \
qw(HOST);\nuse Env qw(HOME);\nuse Env qw(USER);\n\\
n$tmp=\"$ARGV[0].$$\";\nopen (IN, $ARGV[0]);\nopen\
 (OUT, \">$tmp\");\n\nwhile ( <IN>)\n  {\n    $fil\
e=$_;\n    $file=~s/\\r\\n/\\n/g;\n    $file=~s/\\\
n\\r/\\n/g;\n    $file=~s/\\r\\r/\\n/g;\n    $file\
=~s/\\r/\\n/g;\n    print OUT \"$file\";\n  }\nclo\
se (IN);\nclose (OUT);\n\nopen (OUT, \">$ARGV[0]\"\
);\nopen (IN, \"$tmp\");\n\nwhile ( <IN>)\n{\n  pr\
int OUT \"$_\";\n}\nclose (IN);\nclose (OUT);\nunl\
ink ($tmp);\n\n"};
