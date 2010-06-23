#!/usr/bin/env perl
$version="1.00";
$rseed= int(rand(100000))+1;


if ( $#ARGV==-1)
  {
    print "msa2bootstrap -i <input_file> -input <seq|msa|matrix|tree> -n <N-Boostrap> -o <outtree> -tmode <nj|upgma|parsimony|ml> -dmode <kimura> -alignpg <t_coffee | muscle | clustalw> -rtree <file> -stype <prot|cdna|dna> -recompute -system <cygwin|unix>";
    print "\n\t-i: input file, can be sequneces, msa, matrix, trees, type is specified via -input";
    print "\n\t-input: Type of input data";
    print "\n\t\tmsa: msa in fasta format";
    print "\n\t\tseq: compute an msa with -alignpg";
    print "\n\t\tmatrix: phylipp distance matrix fed directly to method -tmode [caveat: tmode=nj or upgma]";
    print "\n\t\ttree: list of newick trees directly fed to consence in order to generate a bootstraped tree";
    
    print "\n\t-n: number of bootstrap replicates";
    print "\n\t-o: name of the output tree. Files are not overwritten. Use -recompute to overwrite existing file";
    print "\n\t-tmode: tree mode: nj|upgma|parsimony|ml";
    print "\n\t-dmode: distance mode";
    print "\n\t-alignpg: program for aligning sequences (t_coffee=default)";
    print "\n\t-rtree: replicate tree file (default: no file)";
    print "\n\t-rmsa: replicate msa file (default: no file)";
    print "\n\t-rmat: replicate matrix file (default: no file)";
    print "\n\t-stype: sequence type: protein, dna or cdna";
    print "\n\t-recompute: force files to be overwritten";
    print "\n\t-system: cygwin|unix";
      

    
    &my_exit (EXIT_FAILURE);
  }
foreach $arg (@ARGV){$command.="$arg ";}

print "CLINE: $command\n";
$threshold=100;
$trim_msa=0;
$stype="prot";
print "msa2bootstrap ";

$system="cygwin";
if(($command=~/\-system (\S+)/))
  {
    $system=$1;
    if ( $system eq "cygwin")
      {
	$exec_extension=".exe";
      }
    elsif ( $system eq "unix")
      {
	$exec_extension="";
	print "system=Unix";die;
      }
    else
      {
	print "msa2boostrap: -system=$system is an unknown mode [FATAL]\n"; die;
      }
    
    print "-system $system ";
  }
if(($command=~/\-stype (\S+)/))
  {
    $stype=$1;
  }
print "-stype=$stype ";



if(($command=~/\-i (\S+)/))
  {
    $msa=$1;
    print "-i $msa ";
  }

if(($command=~/\-rtree (\S+)/))
  {
    $rtree=$1;
    print "-rtree=$rtree ";
  }

if(($command=~/\-rmsa (\S+)/))
  {
    $rmsa=$1;
  }
if(($command=~/\-rmat (\S+)/))
  {
    $rmat=$1;
  }
$input="seq";
if(($command=~/\-input (\S+)/))
  {
    $input=$1;
  }
print "-input=$input ";

$dmode="kimura";
if(($command=~/\-dmode (\S+)/))
  {
    $dmode=$1;
  }
print "-dmode=$dmode ";
$alignpg="muscle";
if(($command=~/\-alignpg (\S+)/))
  {
    $alignpg=$1;
  }
print "-alignpg=$dmode ";

$tmode="nj";
if(($command=~/\-tmode (\S+)/))
  {
    $tmode=$1;
  }
print "-tmode=$tmode ";
$recompute=0;
if(($command=~/\-recompute/))
  {
    $recompute=1;
    print "-recompute ";
  }

$out=$msa;
$out=~s/\..*//;
$out.=".bph";
if(($command=~/\-o (\S+)/))
  {
    $out=$1;
    
  }
print "-out=$out ";
if (-e $out && !$recompute)
  {
    print "\nNo Computation Required $out already exists\n";
    &my_exit (EXIT_SUCCESS);
    
  }

$n=100;
if(($command=~/\-n (\d+)/))
  {
    $n=$1;
  }
print "-n=$n ";
$seed=3;
if(($command=~/\-s (\d+)/))
  {
    $seed=$1;
  }
print "-s=$seed ";

if(($command=~/\-run_name (\d+)/))
  {
    $suffix=$1;
  }
else
  {
    $msa=~/([^.]+)/;
    $suffix=$1;
  }
print "-run_name=$suffix\n";


if ( $input eq "seq")
  {
    $seq=$msa;
    $msa="$suffix.prot_msa";
    
    if ($stype eq "cdna")
      {
	$cdna_seq=$seq;
	$clean_cdna_seq=&vtmpnam();
	$seq=&vtmpnam();
	`t_coffee -other_pg seq_reformat -in $cdna_seq -action +clean_cdna >$clean_cdna_seq`;
	`t_coffee -other_pg seq_reformat -in $clean_cdna_seq -action +translate >$seq`;
	
      }

    if (!-e $msa || $recompute)
      {
	print "\n#####   Compute an MSA With $alignpg\n";
	
	if ( $alignpg eq "t_coffee")
	  {`$alignpg $seq -outfile=$msa >/dev/null 2>/dev/null`;}
	elsif ( $alignpg eq "muscle")
	  {
	    `$alignpg -in $seq > $msa 2>/dev/null`;
	  }
	elsif ( $alignpg eq "clustalw")
	  {
	    `$alignpg -infile=$seq -outfile=$msa -quicktree >/dev/null 2>/dev/null`;
	  }
	elsif ( $align eq "mafft")
	  {
	    `$alignpg $seq > $msa >/dev/null 2>/dev/null`;
	  }
	else
	  {
	    `$alignpg -in=$seq -outfile=$msa`;
	  }
      }
    if (!-e $msa)
      {
	print "\nError: $alignpg Could Not produce the MSA $msa [FATAL]\n";
      }

    if ($stype eq "cdna")
      {
	$msa2="$suffix.cdna_msa";
	`t_coffee -other_pg seq_reformat -in $clean_cdna_seq -in2 $msa -action +thread_dna_on_prot_aln -output fasta_aln  >$msa2`;
	$msa=$msa2;
      }
    
    $input="msa";
  }



$seqboot_o=&vtmpnam();
$seqboot_c=&vtmpnam();

$protdist_o=&vtmpnam();
$protdist_c=&vtmpnam();
if ( $input eq "msa")
  {
    if ($tmode eq "nj" || $tmode eq "upgma"){$input="matrix";}
    
    $lmsa= &vtmpnam ();
    `t_coffee -other_pg seq_reformat -in $msa -output phylip_aln > $lmsa`;
    
    if ( -e "outfile"){unlink ("outfile");}
    # run seqboot
  
    if ( $n>1)
      {
	print "Run SeqBoot .....";
	open (F, ">$seqboot_c");
	print F "$lmsa\nR\n$n\nY\n$seed\n";
	close (F);
	`seqboot$exec_extension  < $seqboot_c`;
	if ( -e "outfile"){ print "[OK]\n";}
	else { print "[FAILED]\n";&my_exit (EXIT_FAILURE);}
	`mv outfile $seqboot_o`;
      }
    else
      {
	`cp $lmsa $seqboot_o`;
      }

    if ($rmsa){`cp $seqboot_o $rmsa`;}
    
    if ($tmode eq "nj" || $tmode eq "upgma")
      {
	if ( $stype eq "prot")
	  {
	    # run protdist
	    print "Run Protdist [dmode=$dmode]";
	    if ($dmode eq "kimura")
	      {
		$dmode="P\nP\nP";
	      }
	    else
	      {
		print "\n$dmode is an unknown mode for Protdist [FATAL:msa2bootstrap.pl]\n";
		&my_exit (EXIT_FAILURE);
	      }
	    open (F, ">$protdist_c");
	    if ($n>1){print F "$seqboot_o\n$dmode\nM\nD\n$n\nY\n";}
	    else {printf F "$seqboot_o\n$dmode\nY\n";}
	    close (F);
	    `protdist$exec_extension  < $protdist_c`;
	    if ( -e "outfile"){ print "[OK]\n";}
	    else { print "[FAILED]\n";&my_exit (EXIT_FAILURE);}
	    `mv outfile $protdist_o`;
	 
	  }
	elsif ( $stype eq "cdna" || $stype eq "dna")
	  {
	    print "Run dnadist [dmode=default";
	    open (F, ">$protdist_c");
	    if ($n>1){print F "$seqboot_o\nM\nD\n$n\nY\n";}
	    else {printf F "$seqboot_o\nY\n";}
	    close (F);
	    `protdist$exec_extension  < $protdist_c`;
	    if ( -e "outfile"){ print "[OK]\n";}
	    else { print "[FAILED]\n";&my_exit (EXIT_FAILURE);}
	    `mv outfile $protdist_o`;
	  }
      }
  }
elsif ( $input eq "matrix")
  {
    $protdist_o=&vtmpnam();
    print "MSA: $msa\n";
    `cp $msa $protdist_o`;
    $n=1;
  }





$nb_o=&vtmpnam();
$nb_c=&vtmpnam();
if ($input eq "matrix" && $tmode ne "parsimony" && $tmode ne "ml")
  {
    print "Run neighbor [tmode=$tmode]";

    if ($tmode eq "nj")
      {
	$tmode="\nN\nN";
      }
    elsif ( $tmode eq "upgma")
      {
	$tmode = "\nN";
      }
    else
      {
	print "\n ERROR: $tmode is an unknown tree computation mode\n";
	&my_exit (EXIT_FAILURE);
      }

    open (F, ">$nb_c");
    if ($n>1){print F "$protdist_o$tmode\nM\n$n\n$seed\nY\n";}
    else {print F "$protdist_o$tmode\nY\n";}
    close (F);

    `neighbor$exec_extension  < $nb_c`;
    if ( -e "outtree"){ print "[Neighbor OK]\n";}
    else { print "[FAILED]\n";&my_exit (EXIT_FAILURE);}
    `mv outtree $nb_o`;
    unlink ("outfile");
  }
elsif ($input eq "msa" && $tmode eq "parsimony")
  {
    if ( -e "outfile"){unlink ("outfile");}
    if ( -e "outtree"){unlink ("outtree");}
    
    if ($stype eq "prot")
      {
	print "Run protpars [tmode=$tmode]";
	open (F, ">$nb_c");
	if ($n>1){print F "$seqboot_o\nM\nD\n$n\n$seed\n10\nY\n";}
	else {print F "$seqboot_o\nY\n";}
	close (F);
	`protpars$exec_extension  < $nb_c`;
      }
    elsif ( $stype eq "dna" || $stype eq "cdna")
      {
	print "Run dnapars [tmode=$tmode]";
	open (F, ">$nb_c");
	if ($n>1){print F "$seqboot_o\nM\nD\n$n\n$seed\n10\nY\n";}
	else {print F "$seqboot_o\nY\n";}
	close (F);
	`dnapars$exec_extension  < $nb_c`;
      }
    if ( -e "outtree"){ print "[OK]\n";}
    else { print "[FAILED]\n";&my_exit (EXIT_FAILURE);}
    `mv outtree $nb_o`;
   unlink ("outfile");
  }
elsif ($input eq "msa" && $tmode eq "ml")
  {
    if ( -e "outfile"){unlink ("outfile");}
    if ( -e "outtree"){unlink ("outtree");}
    
    if ($stype eq "prot")
      {
	print "Error: ML impossible with Protein Sequences [ERROR]";
	&my_exit (EXIT_FAILURE);
      }
    elsif ( $stype eq "dna" || $stype eq "cdna")
      {
	print "Run dnaml [tmode=$tmode]";
	open (F, ">$nb_c");
	if ($n>1){print F "$seqboot_o\nM\nD\n$n\n$seed\n10\nY\n";}
	else {print F "$seqboot_o\nY\n";}
	close (F);
	`dnaml$exec_extension  < $nb_c`;
      }
    if ( -e "outtree"){ print "[OK]\n";}
    else { print "[FAILED]\n";&my_exit (EXIT_FAILURE);}
    `mv outtree $nb_o`;
   unlink ("outfile");
  }


else
  {
    `cp $msa $nb_o`;
    $n=2;
  }

if ($rmsa && -e $seqboot_o){print "\nOutput List of $n Replicate MSA: $rmsa\n";`cp $seqboot_o $rmsa`;}
if ($rmat && -e $protdist_o){print "\nOutput List of $n Replicate MATRICES: $rmat\n";`cp $protdist_o $rmat`;}
if ($rtree && -e $nb_o){print "\nOutput List of $n Replicate TREES: $rtree\n";`cp $nb_o $rtree`;}



$con_o=&vtmpnam();
$con_c=&vtmpnam();
if ($n >1)
  {
    print "Run Consense.....";
    open (F, ">$con_c");
    print F "$nb_o\nY\n";
    close (F);
    `consense$exec_extension  < $con_c`;
    if ( -s "outtree"  > 0) { print "[OK]\n";}
    else { print "[FAILED]\n";&my_exit (EXIT_FAILURE);}
    `mv outtree $con_o`;
    unlink ("outfile");
  }
else
  {
    `cp $nb_o $con_o`;
  }


`cp $con_o $out`;
if ( !-e $out)
  {
    print "Tree Computation failed [FAILED]\n";
    &my_exit (EXIT_FAILURE);
  }
elsif ($n>1)
  {
    print "\nOutput Bootstrapped Tree: $out\n";
    $avg=`t_coffee -other_pg seq_reformat -in $out -action +avg_bootstrap`;
    $avg=~s/\n//g;
    print "$avg\n";
  }
else
  {
    print "\nOutput Tree: $out\n";
  }

open (F, "$out");
while (<F>)
  {
    
    $tree.=$_;
  }
close (F);
$tree=~s/\n//g;
print "BPH: $tree\n";


&my_exit (EXIT_SUCCESS);

sub my_exit 
  {
    my $m=@_[0];
    &clean_vtmpnam();
    exit ($m);
  }
sub vtmpnam 
  {
    my $file;


    $ntmp++;
    $file="tmp4msa2bootstrap.$rseed.$$.$ntmp";
    
    push (@tmpfile, $file);
    return $file;
  }
sub clean_vtmpnam 
  {
    my $t;
    foreach $t (@tmpfile)
      {
	if ( -e $t){unlink ($t)};
      }
  }

