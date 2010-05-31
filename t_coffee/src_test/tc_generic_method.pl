#!/usr/bin/env perl
use Env;

$BLAST_MAX_NRUNS=2;
$EXIT_SUCCESS=0;
$EXIT_FAILURE=1;

use Cwd;
$REF_EMAIL="";


$tmp_dir="";
$init_dir="";
$program="tc_generic_method.pl";

$test=0;
if ($test==1)
  {
    $SERVER="NCBI";
    $query=$ARGV[0];
    $hitf=$ARGV[1];
    %s=read_fasta_seq($query);
    @sl=keys(%s);
    &blast_xml2profile ("xx", $s{$sl[0]}{seq},$maxid,$minid,$mincov, $hitf);
    myexit ($EXIT_FAILURE);
  }

foreach $v(@ARGV){$cl.="$v ";}
($mode)=&my_get_opt ( $cl, "-mode=",1,0);

($A)=(&my_get_opt ( $cl, "-name1=",0,0));
($B)=(&my_get_opt ( $cl, "-name2=",0,0));
($TMPDIR)=(&my_get_opt ( $cl, "-tmpdir=",0,0));
($CACHE)=(&my_get_opt ( $cl, "-cache=",0,0));
($SERVER)=((&my_get_opt ( $cl, "-server=",0,0)));
($EMAIL)=((&my_get_opt ( $cl, "-email=",0,0)));

if (!$A){$A="A";}
if (!$B){$B="B";}


if (!$TMPDIR)
  {
    $HOME=$ENV{HOME};
    if ($ENV{TMP_4_TCOFFEE}){$TMPDIR=$ENV{TMP_4_TCOFFEE};}
    else{$TMPDIR="$HOME/.t_coffee/tmp/";}
  }
if ( ! -d $TMPDIR)
  {
    mkdir $TMPDIR;
  }
if ( ! -d $TMPDIR)
  {
    print "ERROR: Could not create temporary dir: $TMPDIR\n";
    myexit ($EXIT_FAILURE);
  }

$EMAIL=~s/XEMAILX/\@/g;
if (!$EMAIL)
  {
    if ($ENV{EMAIL_4_TCOFFEE}){$EMAIL=$ENV{EMAIL_4_TCOFFEE};}
    elsif ($ENV{EMAIL}){$EMAIL=$ENV{EMAIL};}
    else {$EMAIL=$REF_EMAIL;}
  }

($maxid,$minid,$mincov)=(&my_get_opt ( $cl, "-maxid=",0,0, "-minid=",0,0,"-mincov=",0,0));
if (!$cl=~/\-maxid\=/){$maxid=95;}
if (!$cl=~/\-minid\=/){$minid=35;}
if (!$cl=~/\-mincov\=/){$mincov=80;}



if ($mode eq "seq_msa")
  {
    &seq2msa($mode,&my_get_opt ( $cl, "-infile=",1,1, "-method=",1,2, "-param=",0,0, "-outfile=",1,0));
  }

elsif ( $mode eq "thread_pair")
  {
    &seq2thread_pair($mode,&my_get_opt ( $cl, "-infile=",1,1, "-pdbfile1=",1,1, "-method=",1,2,"-param=",0,0, "-outfile=",1,0, ));
  }
elsif ( $mode eq "pdbid_pair")
  {
    &seq2pdbid_pair($mode,&my_get_opt ( $cl, "-pdbfile1=",1,0, "-pdbfile2=",1,0, "-method=",1,2,"-param=",0,0, "-outfile=",1,0, ));
  }
elsif ( $mode eq "pdb_pair")
  {
    &seq2pdb_pair($mode,&my_get_opt ( $cl, "-pdbfile1=",1,1, "-pdbfile2=",1,1, "-method=",1,2,"-param=",0,0, "-outfile=",1,0, ));
  }
elsif ( $mode eq "profile_pair")
  {
     &seq2profile_pair($mode,&my_get_opt ( $cl, "-profile1=",1,1, "-profile2=",1,1, "-method=",1,2,"-param=",0,0, "-outfile=",1,0, ));
  }
elsif ( $mode eq "pdb_template")
  {
    &blast2pdb_template ($mode,&my_get_opt ( $cl, "-infile=",1,1, "-database=",1,0, "-method=",1,0, "-outfile=",1,0));
  }
elsif ( $mode eq "profile_template")
  {
    &psiblast2profile_template ($mode,&my_get_opt ( $cl, "-infile=",1,1, "-database=",1,0, "-method=",1,0, "-outfile=",1,0));
  }
elsif ( $mode eq "psiprofile_template")
  {
    &psiblast2profile_template ($mode,&my_get_opt ( $cl, "-infile=",1,1, "-database=",1,0, "-method=",1,0, "-outfile=",1,0));
  }
elsif ( $mode eq "RNA_template")
  {
    &seq2RNA_template ($mode,&my_get_opt ( $cl, "-infile=",1,1, "-outfile=",1,0));
  }
elsif ( $mode eq "tm_template")
  {
    &seq2tm_template ($mode,&my_get_opt ( $cl, "-infile=",1,1,"-arch=",1,1,"-psv=",1,1, "-outfile=",1,0,));
  }
elsif ( $mode eq "psitm_template")
  {
    &seq2tm_template ($mode,&my_get_opt ( $cl, "-infile=",1,1,"-arch=",1,1,"-psv=",1,1, "-outfile=",1,0,));
  }
elsif ( $mode eq "ssp_template")
  {
    &seq2ssp_template ($mode,&my_get_opt ( $cl, "-infile=",1,1,"-seq=",1,1,"-obs=",1,1, "-outfile=",1,0));
  }
elsif ( $mode eq "psissp_template")
  {
    &seq2ssp_template ($mode,&my_get_opt ( $cl, "-infile=",1,1,"-seq=",1,1,"-obs=",1,1, "-outfile=",1,0));
  }
elsif ( $mode eq "rna_pair")
{
    &seq2rna_pair($mode,&my_get_opt ( $cl, "-pdbfile1=",1,1, "-pdbfile2=",1,1, "-method=",1,2,"-param=",0,0, "-outfile=",1,0, ));
}elsif ( $mode eq "calc_rna_template")
{
    &calc_rna_template($mode,&my_get_opt ( $cl, "-infile=",1,1,"-pdbfile=",1,1, "-outfile=",1,0));
}
else
  {
    print STDERR "$mode is an unknown mode of tc_generic_method.pl [FATAL]\n";
  }
myexit ($EXIT_SUCCESS);
sub seq2ssp_template
  {
  my ($mode, $infile,$gor_seq,$gor_obs,$outfile)=@_;
  my %s, %h;
  my $result;
  my (@profiles);
  &set_temporary_dir ("set",$infile,"seq.pep");
  %s=read_fasta_seq ("seq.pep");

  
  open (R, ">result.aln");
  
  #print stdout "\n";
  foreach $seq (keys(%s))
    {
      
      open (F, ">seqfile");
      $s{$seq}{seq}=uc$s{$seq}{seq};
      print (F ">$s{$seq}{name}\n$s{$seq}{seq}\n");
      close (F);
      $lib_name="$s{$seq}{name}.ssp";
      $lib_name=&clean_file_name ($lib_name);
      
      if ($mode eq "ssp_template"){&seq2gor_prediction ($s{$seq}{name},$s{$seq}{seq}, "seqfile", $lib_name,$gor_seq, $gor_obs);}
      elsif ($mode eq "psissp_template")
	{
	  &seq2msa_gor_prediction ($s{$seq}{name},$s{$seq}{seq},"seqfile", $lib_name,$gor_seq, $gor_obs);
	}
    
      if ( !-e $lib_name)
	{
	  print STDERR ("GORIV failed to compute the secondary structure of $s{$seq}{name} [FATAL:$mode/$method/$program]\n");
	  myexit ($EXIT_FAILURE);
	}
      else
	{
	  print stdout "\tProcess: >$s{$seq}{name} _E_ $lib_name \n";
	  print R ">$s{$seq}{name} _E_ $lib_name\n";
	}
      unshift (@profiles, $lib_name);
    }
  close (R);
  &set_temporary_dir ("unset",$mode, $method,"result.aln",$outfile, @profiles);
}

sub seq2tm_template
  {
  my ($mode, $infile,$arch,$psv,$outfile)=@_;
  my %s, %h;
  my $result;
  my (@profiles);
  &set_temporary_dir ("set",$infile,"seq.pep");
  %s=read_fasta_seq ("seq.pep");

  
  open (R, ">result.aln");
  
  #print stdout "\n";
  foreach $seq (keys(%s))
    {
      open (F, ">seqfile");
      print (F ">$s{$seq}{name}\n$s{$seq}{seq}\n");
      close (F);
      $lib_name="$s{$seq}{name}.tmp";
      $lib_name=&clean_file_name ($lib_name);

      if ($mode eq "tm_template")
	{
	  &safe_system ("t_coffee -other_pg fasta_seq2hmmtop_fasta.pl -in=seqfile -out=$lib_name -arch=$arch -psv=$psv");
	}
      elsif ( $mode eq "psitm_template")
	{
	  &seq2msa_tm_prediction ($s{$seq}{name},$s{$seq}{seq},"seqfile", $lib_name,$arch, $psv);
	}
      if ( !-e $lib_name)
	{
	  print STDERR ("RNAplfold failed to compute the secondary structure of $s{$seq}{name} [FATAL:$mode/$method/$program]\n");
	  myexit ($EXIT_FAILURE);
	}
      else
	{
	  print stdout "\tProcess: >$s{$seq}{name} _T_ $lib_name\n";
	  print R ">$s{$seq}{name} _T_ $lib_name\n";
	}
      unshift (@profiles, $lib_name);
    }
  close (R);
  &set_temporary_dir ("unset",$mode, $method,"result.aln",$outfile, @profiles);
}

sub seq2RNA_template
  {
  my ($mode, $infile,$outfile)=@_;
  my %s, %h, ;
  my $result;
  my (@profiles);
  &set_temporary_dir ("set",$infile,"seq.pep");
  %s=read_fasta_seq ("seq.pep");

  
  open (R, ">result.aln");
  
  #print stdout "\n";
  foreach $seq (keys(%s))
    {
      open (F, ">seqfile");
      print (F ">$s{$seq}{name}\n$s{$seq}{seq}\n");
      close (F);
      $lib_name="$s{$seq}{name}.rfold";
      $lib_name=&clean_file_name ($lib_name);
      &safe_system ("t_coffee -other_pg RNAplfold2tclib.pl -in=seqfile -out=$lib_name");
      
      if ( !-e $lib_name)
	{
	  print STDERR ("RNAplfold failed to compute the secondary structure of $s{$seq}{name} [FATAL:$mode/$method/$program]\n");
	  myexit ($EXIT_FAILURE);
	}
      else
	{
	  print stdout "\tProcess: >$s{$seq}{name} _F_ $lib_name\n";
	  print R ">$s{$seq}{name} _F_ $lib_name\n";
	}
      unshift (@profiles, $lib_name);
    }
  close (R);
  &set_temporary_dir ("unset",$mode, $method,"result.aln",$outfile, @profiles);
}

sub psiblast2profile_template 
  {
  my ($mode, $infile, $db, $method, $outfile)=@_;
  my %s, %h, ;
  my ($result,$psiblast_output,$profile_name,@profiles);
  
  &set_temporary_dir ("set",$infile,"seq.pep");
  %s=read_fasta_seq ("seq.pep");
  open (R, ">result.aln");
  
  #print stdout "\n";
  foreach $seq (keys(%s))
    {
      open (F, ">seqfile");
      print (F ">$A\n$s{$seq}{seq}\n");
      close (F);
      $psiblast_output=&run_blast ($s{$seq}{name},$method, $db, "seqfile","outfile");
      if ( -e $psiblast_output)
	{
	  %profile=blast_xml2profile($s{$seq}{name}, $s{$seq}{seq},$maxid, $minid,$mincov,$psiblast_output);
	  unlink ($psiblast_output);
	  
	  $profile_name="$s{$seq}{name}.prf";
	  $profile_name=&clean_file_name ($profile_name);
	  unshift (@profiles, $profile_name);
	  output_profile ($profile_name, %profile);
	  print stdout "\tProcess: >$s{$seq}{name} _R_ $profile_name [$profile{n} Seq.] [$SERVER/blast/$db][$CACHE_STATUS]\n";
	  print R ">$s{$seq}{name} _R_ $profile_name\n";
	}
    }
  close (R);
  &set_temporary_dir ("unset",$mode, $method,"result.aln",$outfile, @profiles);
}

sub blast2pdb_template 
  {
  my ($mode, $infile, $db, $method, $outfile)=@_;
  my %s, %h, ;
  my ($result,$blast_output);
  &set_temporary_dir ("set",$infile,"seq.pep");
  %s=read_fasta_seq ("seq.pep");
  open (R, ">result.aln");
  
 
  #print stdout "\n";
  foreach $seq (keys(%s))
    {
      open (F, ">seqfile");
      print (F ">$A\n$s{$seq}{seq}\n");
      close (F);
      
      $blast_output=&run_blast ($s{$seq}{name},$method, $db, "seqfile","outfile");
      %p=blast_xml2profile($s{$seq}{name}, $s{$seq}{seq},$maxid, $minid,$mincov,$blast_output);
      unlink ($blast_output);
      if ($p{n}>1)
	{
	  $pdbid=id2pdbid($p{1}{identifyer});
	  if ( length ($pdbid)>5){$pdbid=id2pdbid($p{1}{definition});}
	  
	  print R ">$s{$seq}{name} _P_ $pdbid\n";
	  print stdout "\tProcess: >$s{$seq}{name} _P_ $pdbid [$SERVER/blast/$db][$CACHE_STATUS]\n";
	}
      else
	{
	  print R ">$s{$seq}{name}\n";
	  print stdout "\tProcess: >$s{$seq}{name} _P_ No Template Found [$SERVER/blast/$db][$CACHE_STATUS]\n";
	}
    }
  close (R);
  &set_temporary_dir ("unset",$mode, $method,"result.aln",$outfile);
}
sub blast_msa
  {
    my ($infile,$outfile)=@_;
    my ($a, %seq);
    %s1=&read_fasta_seq ($infile);
    foreach $s (keys (%s1))
      {
	$i=$s1{$s}{order};
	$s{$i}{name}=$s;
	$s{$i}{seq}=$s1{$s}{seq};
	$s{$i}{len}=length( $s{$i}{seq});
	$s{n}++;
      }
    `formatdb -i $infile`;
    `blastpgp -i $infile -d $infile -m7 -j4 > io`;
    &set_blast_type ("io");
    
    %FB=&xml2tag_list ("io", "BlastOutput");
    
    open (F, ">$outfile");
    print F "! TC_LIB_FORMAT_01\n";
    print F "$s{n}\n";
    for ( $a=0; $a<$s{n}; $a++)
      {
	print F "$s{$a}{name} $s{$a}{len} $s{$a}{seq}\n";
      }
    for ( $a=0; $a<$FB{n}; $a++)
      {
	%p=blast_xml2profile ($s{$a}{name}, $s{$a}{seq},100, 0, 0, $FB{$a}{body});
	for ($b=1; $b<$p{n}; $b++)
	  {
	    my $l=length ($p{$b}{Qseq});
	    my $hit=$p{$b}{definition};
	    my $Qstart=$p{$b}{Qstart};
	    my $Hstart=$p{$b}{Hstart};
	    my $identity=$p{$b}{identity};
	    my @lrQ=split (//,$p{$b}{Qseq});
	    my @lrH=split (//,$p{$b}{Hseq});
	    my $i= $s1{$s{$a}{name}}{order}+1;
	    my $j= $s1{$hit}{order}+1;
	    #if ( $j==$i){next;}
	    printf F "# %d %d\n", $i, $j;
	    #  print  F "\n$p{$b}{Qseq} ($Qstart)\n$p{$b}{Hseq} ($Hstart)";
	    for ($c=0; $c<$l; $c++)
	      {
		my $rQ=$lrQ[$c];
		my $rH=$lrH[$c];
		my $n=0;
		
		if ($rQ ne "-"){$n++, $Qstart++;}
		if ($rH ne "-"){$n++; $Hstart++;}
		
		if ( $n==2)
		  {
		    printf F "\t%d %d %d\n", $Qstart-1, $Hstart-1,$identity;
		  }
	      }
	  }
      }
    print F "! SEQ_1_TO_N\n";
    close (F);
    return $output;
  
  }

sub seq2msa
  {
    my ($mode, $infile, $method, $param, $outfile)=@_;
    &set_temporary_dir ("set",$infile,"seq.pep");
    $param.=" >/dev/null 2>&1 ";
    
    #make sure test.pep is in FASTA
    &safe_system ("t_coffee -other_pg seq_reformat -in seq.pep -output fasta_seq > x");
    `mv x seq.pep`;
    
    if ( $method eq "blastpgp")
      {
	&blast_msa ("seq.pep", "result.aln");
      }
    elsif ( $method eq "muscle")
      {
	`muscle -in seq.pep -out result.aln $param`;
      }
    elsif ( $method eq "probcons")
      {
	`probcons seq.pep >result.aln 2>/dev/null`;
      }
    elsif ( $method eq "mafft")
      {
	`mafft --quiet --localpair --maxiterate 1000 seq.pep> result.aln  2>/dev/null`
      }
    elsif ( $method=~/prank/)
      {
	`$method -d=seq.pep -o=result.aln -quiet 2>/dev/null`;
	`mv result.aln.1.fas result.aln`;
      }
    else
      {
	`$method -infile=seq.pep -outfile=result.aln`;
      }
    
    &set_temporary_dir ("unset",$mode, $method,"result.aln",$outfile);
    myexit ($EXIT_SUCCESS);
  }

sub seq2thread_pair
  {
    my ($mode, $infile, $pdbfile1, $method, $param, $outfile)=@_;
    &set_temporary_dir ("set",$infile,"seq.pep",$pdbfile1,"struc.pdb");
    if ($method eq "fugueali")
      {
	#Env Variable that need to be defined for Fugue
	if (!$ENV{FUGUE_LIB_LIST}){$ENV{FUGUE_LIB_LIST}="DUMMY";}
	if (!$ENV{HOMSTRAD_PATH})  {$ENV{HOMSTRAD_PATH}="DUMMY";}
	if (!$ENV{HOMS_PATH}){$ENV{HOMS_PATH}="DUMMY";}
	
	`joy struc.pdb >x 2>x`;
	&check_file("struc.tem", "Joy failed [FATAL:$program/$method]");
	`melody -t struc.tem >x 2>x`;
	&check_file("struc.tem", "Melody failed [FATAL:$program/$method]");
	`fugueali -seq seq.pep -prf struc.fug -print > tmp_result.aln`;
	
	&check_file("tmp_result.aln", "Fugue failed [FATAL:$program/$method]");
	&safe_system ("t_coffee -other_pg seq_reformat -in tmp_result.aln -output fasta_aln >result.aln");
      }
    elsif ( $method eq "t_coffee")
      {
	&safe_system ("t_coffee -in Pstruc.pdb Sseq.pep Mslow_pair -outfile result.aln -quiet");
      }
    else
      {
	&safe_system ("$method -infile=seq.pep -pdbfile1=struc.pdb -outfile=result.aln $param>x 2>x");
      }
    &set_temporary_dir ("unset",$mode,$method,"result.aln",$outfile);
    myexit ($EXIT_SUCCESS);
  }
sub seq2pdbid_pair
  {
    my ($mode, $pdbfile1, $pdbfile2, $method, $param, $outfile)=@_;
    my ($name);

    
    &set_temporary_dir ("set");
    $name=$pdbfile1." ".$pdbfile2;

    if (    &cache_file("GET","","$name","$method","dali",$outfile,"EBI"))
      {return $outfile;}
    else
      {
	if ($method eq "dalilite")
	  {
	    $pdbfile1=~/(....)(.)/;
	    $id1=$1; $c1=$2;
	    
	    $pdbfile2=~/(....)(.)/;
	    $id2=$1; $c2=$2;
	    
	    $command="t_coffee -other_pg dalilite.pl --pdb1 $id1 --chainid1 $c1 --pdb2 $id2 --chainid2 $c2 --email=$EMAIL  >dali_stderr 2>dali_stderr";
	    $dali=`$command`;
	    
	    open (F, "dali_stderr");
	    while (<F>)
	      {
		if ( /JobId: dalilite-(\S+)/)
		{
		  $jobid=$1;
		}
	      }
	    close (F);
	    unlink ("dali_stderr");
	    
	    $output1="dalilite-$jobid.txt";
	    if ( -e $output1)
	      {
		unlink ($output1);
		&url2file ("http://www.ebi.ac.uk/Tools/es/cgi-bin/jobresults.cgi/dalilite/dalilite-$jobid/aln.html", "output2");
		
		if ( -e "output2")
		  {
		    my ($seq1, $seq2);
		    $seq1=$seq2="";
		    
		    open (F, "output2");
		    while (<F>)
		      {
			$l=$_;
			if ( $l=~/Query\s+(\S+)/)
			  {
			    $seq1.=$1;
			  }
			elsif ( $l=~/Sbjct\s+(\S+)/)
			  {
			    $seq2.=$1;
			  }
		      }
		    close (F);
		    unlink ("output2");
		    if ($seq1 ne "" && $seq2 ne "")
		      {
			$output3=">$A\n$seq1\n>$B\n$seq2\n";
			$output3=~s/\./-/g;
			open (F, ">result.aln");
			print F "$output3";
			close (F);
		      }
		  }
	      }
	  }
      }
    &cache_file("SET","","$name","$method","dali","result.aln","EBI");
    &set_temporary_dir ("unset",$mode, $method, "result.aln",$outfile);
    myexit ($EXIT_SUCCESS);
  }
sub seq2pdb_pair
  {
    my ($mode, $pdbfile1, $pdbfile2, $method, $param, $outfile)=@_;
    
    &set_temporary_dir ("set",$pdbfile1,"pdb1.pdb",$pdbfile2,"pdb2.pdb");
    if ($method eq "t_coffee")
      {
	&safe_system ("t_coffee -in Ppdb1.pdb Ppdb2.pdb -quiet -outfile=result.aln");
      }
    elsif ( $method eq "TMalign")
      {
	if ( &safe_system ("TMalign pdb1.pdb pdb2.pdb >tmp1")==$EXIT_SUCCESS)
	  {
	    `tail -4 tmp1 > tmp2`;
	    
	    open (F, "tmp2");
	    while (<F>)
	      {
		unshift(@l, $_);
	      }
	    close (F);
	    open (F, ">result.aln");
	    $l[3]=~s/[^a-zA-Z0-9-]/\-/g;
	    $l[1]=~s/[^a-zA-Z0-9-]/\-/g;
	    print F ">$A\n$l[3]\n>$B\n$l[1]\n";
	    close (F);
	  }
	else
	  {
	    print "ERROR: TMalign failed to align the considered structures[tc_generic_method.pl]\n";
	    `rm result.aln >/dev/null 2>/dev/null`;
	  }
      }
    elsif ( $method eq "mustang")
      {
	if ( &safe_system ("mustang -i pdb1.pdb pdb2.pdb -F fasta >/dev/null 2>/dev/null")==$EXIT_SUCCESS)
	  {
	    `mv results.afasta result.aln`;
	  }
	else
	  {
	    print "ERROR: mustang failed to align the considered structures[tc_generic_method.pl]\n";
	    `rm result.aln >/dev/null 2>/dev/null`;
	  }
      }
    else
      {
	if ( &safe_system ("$method -pdbfile1=pdb1.pep -pdbfile2=pdb2.pdb -outfile=result.aln $param>x 2>x")==$EXIT_SUCCESS)
	  {
	    `mv results.afasta result.aln`;
	  }
	else
	  {
	    print "ERROR: $method failed to align the considered structures[tc_generic_method.pl]\n";
	    `rm result.aln >/dev/null 2>/dev/null`;
	  }
      }
    &set_temporary_dir ("unset",$mode, $method, "result.aln",$outfile);
    myexit ($EXIT_SUCCESS);
  }

sub seq2profile_pair
  {
    my ($mode, $profile1, $profile2, $method, $param, $outfile)=@_;
    
    
    if ($method eq "clustalw")
      {
	&set_temporary_dir ("set",$profile1,"prf1.aln",$profile2,"prf2.aln");
	`clustalw -profile1=prf1.aln -profile2=prf2.aln -outfile=result.aln`;
	&set_temporary_dir ("unset",$mode, $method, "result.aln",$outfile);
      }
    elsif ( $method eq "hhalign")
      {
	hhalign ( $profile1,$profile2,$outfile,$param);
      }
    else
      {
	
	`$method -profile1=prf1.aln -profile2=prf2.aln -outfile=result.aln $param>x 2>x`;
      }
   
    myexit ($EXIT_SUCCESS);
  }

sub pg_is_installed
  {
    my @ml=@_;
    my $r, $p, $m;
    my $supported=0;
    
    my $p=shift (@ml);
    $r=`which $p 2>/dev/null`;
    if ($r eq ""){return 0;}
    else {return 1;}
  }
sub check_pg_is_installed
  {
    my @ml=@_;
    my $r=&pg_is_installed (@ml);
    if (!$r)
      {
	print STDERR "\nProgram $p Supported but Not Installed on your system [FATAL:tc_generic_method]\n";
	myexit ($EXIT_FAILURE);
      }
    else
      {
	return 1;
      }
  }
sub set_temporary_dir
  {
    my @list=@_;
    my $dir_mode, $a, $mode, $method;

    $dir_mode=shift (@list);

    
    if ( $dir_mode eq "set")
      {
	$initial_dir=cwd();
	if ( !$tmp_dir)
	  {
	    srand;
	    $rand=rand (100000);
	    $tmp_dir="$TMPDIR/tmp4tcoffee_profile_pair_dir_$$_P_$rand";
	  }
	if ( !-d $tm_dir)
	  {
	    `mkdir $tmp_dir`;
	  }
	
	for ( $a=0; $a<=$#list; $a+=2)
	      {
		`cp $list[$a] $tmp_dir/$list[$a+1]`;
	      }
	chdir $tmp_dir;
      }
    elsif ( $dir_mode eq "unset")
      {
	$mode=shift (@list);
	$method=shift (@list);
	
	if (!-e $list[0])
	  {
	    print STDERR ("Program $method failed to produce $list[1] [FATAL:$mode/$method/$program]\n");
	    myexit ($EXIT_FAILURE);
	  }
	else
	  {
	    chdir $initial_dir;
	    # `t_coffee -other_pg seq_reformat -in $tmp_dir/$list[0] -output fasta_aln -out $tmp_dir/result2.aln`;
	    `cp $tmp_dir/$list[0] $tmp_dir/result2.aln`;
	    if ( $list[1] eq "stdout")
	      {
		open (F, "$tmp_dir/result2.aln");
		while (<F>){print $_;}close(F);
	      }
	    else
	      {
		`mv $tmp_dir/result2.aln $list[1]`;
	      }
	    shift (@list); shift (@list);
	    foreach $f (@list)
	      {
		`mv $tmp_dir/$f .`;
	      }
	  }
      }
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

sub myexit
  {
    my $code=@_[0];
    &clean_dir ($tmp_dir);
    exit ($code);
  }

sub my_get_opt
  {
    my @list=@_;
    my $cl, $a, $argv, @argl;
    
    @argl=();
    $cl=shift @list;
    for ( $a=0; $a<=$#list; $a+=3)
      {
	$option=$list[$a];
	$optional=$list[$a+1];
	$status=$list[$a+2];
	$argv="";
	if ($cl=~/$option(\S+)/){$argv=$1;}
	@argl=(@argl,$argv);
	
	
	#$optional:0=>optional
	#$optional:1=>must be set
	#$status: 0=>no requirement
	#$status: 1=>must be an existing file
	#$status: 2=>must be an installed package
	

	if ($optional==0){;}
	elsif ( $optional==1 && $argv eq "")
	  {
	    print STDERR "ERROR: Option $option must be set [FATAL:$program/$mode/$method]\n";
	    myexit ($EXIT_FAILURE);
	  }
	if ($status==0){;}
	elsif ($status ==1 && $argv ne "" && !-e $argv)
	  {
	    print STDERR "ERROR: File $argv must exist [FATAL:$program/$mode/$method]\n";
	    myexit ($EXIT_FAILURE);
	  }
	elsif ( $status==2 && $argv ne "" && &check_pg_is_installed ($argv)==0)
	  {
	    print STDERR "ERROR: $argv is not installed [FATAL:$program/$mode/$method]\n";
	    myexit ($EXIT_FAILURE);
	  }
      }

    return @argl;
    }

sub check_file 
  {
    my ($file, $msg)=@_;

    if ( !-e $file)
      {
	print "\n$msg\n";
	myexit ($EXIT_FAILURE);
      }
    }
sub hhalign
  {
    my ($aln1, $aln2, $outfile, $param)=@_;
    my $h1, $h2;
    
    $h{0}{index}=0;
    $h{1}{index}=1;
    
    $h{0}{aln}=$aln1;
    $h{1}{aln}=$aln2;

   

    %{$h{0}}=aln2psi_profile (%{$h{0}});
    %{$h{1}}=aln2psi_profile (%{$h{1}});

    $param=~s/#S/ /g;
    $param=~s/#M/\-/g;
    $param=~s/#E/\=/g;
    

    
    $command="hhalign -i $h{0}{a3m} -t $h{1}{a3m} -tc $outfile.tmp -rank 1 -mapt 0 $param";
    `$command`;
    
  #  `hhalign -i $h{0}{a3m} -t $h{1}{a3m} -tc $outfile.tmp -rank 1 -mapt 0 -gapf 0.8 -gapg 0.8`;
    

    # To run global use the following
    
    open (I, "$outfile.tmp");
    open (O, ">$outfile");
    $h{0}{cons}=s/\./x/g;
    $h{1}{cons}=s/\./x/g;

    print O "! TC_LIB_FORMAT_01\n2\n$h{0}{name} $h{0}{len} $h{0}{seq}\n$h{1}{name} $h{1}{len} $h{1}{seq}\n#1 2\n";
    
    while (<I>)
      {
	if (/(\d+)\s+(\d+)\s+(\d+)/)
	  {
	    print O "\t$h{0}{$1}\t$h{1}{$2}\t$3\n";
	  }
      }
    print O "! SEQ_1_TO_N\n";

    close (O);
    close (I);
  }

sub aln2psi_profile
  {
    my (%h)=@_;
    my ($aln,$i,$hv, $a, @c, $n);
   
    $i=$h{index};
    $aln=$h{aln};

    `cp $aln $$.hhh_aln`;
    $command="t_coffee -other_pg seq_reformat -in $aln -output hasch";
    $hv=`$command`;chomp ($hv);
    
    $h{a2m}="$tmp/$hv.tmp4hhpred.a2m";
    $h{a3m}="$tmp/$hv.tmp4hhpred.a3m";
    if ( -e $h{a3m}){;}
    else
      {
	`hhconsensus  -M 50 -i $h{aln} -oa2m $h{a2m}`;
	if (!-e $h{a2m})
	  {
	    print STDERR "Program tc_generic_method.pl FAILED to run:\n\thhconsensus  -M 50 -i $h{aln} -oa2m $h{a2m}";
	    myexit ($EXIT_FAILURE);
	  }
	
	`hhconsensus  -M 50 -i $h{aln} -oa3m $h{a3m}`;
	if (!-e $h{a3m})
	  {
	    print STDERR "Program tc_generic_method.pl FAILED to run:\n\thhconsensus  -M 50 -i $h{aln} -oa3m $h{a3m}";
	    myexit ($EXIT_FAILURE);
	  }
       `buildali.pl $h{a3m} -n 1`;
      }
    
    
    $h{a2m_seq}=`head -n 2 $h{a2m} | grep -v ">"`;chomp ($h{a2m_seq});
    $h{a3m_seq}=`head -n 2 $h{a3m} | grep -v ">"`;chomp ($h{a3m_seq});
    $h{cons}=$h{a2m_seq};
    $h{seq}=`head -n 2 $h{aln} | grep -v ">"`;chomp ($h{seq});
    
    

    @c=split (//, $h{cons});
    $h{len}=$#c+1;
    for ($n=0,$a=0, $b=0; $a<$h{len};$a++)
      {
	if ( $c[$a]=~/[A-Z]/)
	  {
	    $h{++$n}=++$b;

	  }
	elsif ( $c[$a]=~/[a-z\.]/)
	  {
	    ++$b;
	  }
      }
    
    $name=`head -n 2 $h{aln} | grep ">"`;
    $name=~/\>(\S+)/;
    $h{name}=$1;
    
    `cp $h{a2m} $i.a2m`;
    `cp $h{a3m} $i.a3m`;
    `cp $h{aln} $i.hh_aln`;
    
    return %h;
  }

sub read_fasta_seq 
  {
    my $f=@_[0];
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
    @com =($s=~/>\S*(.*)\n([^>]*)/g);

    
    $nseq=$#name+1;
    
    for ($a=0; $a<$nseq; $a++)
      {
	my $s;
	my $n=$name[$a];
	$hseq{$n}{name}=$n;
	$seq[$a]=~s/[^A-Za-z]//g;
	$hseq{$n}{order}=$a;
	$hseq{$n}{seq}=$seq[$a];
	$hseq{$n}{com}=$com[$a];
	
      }
    return %hseq;
  }

sub file_contains 
  {
    my ($file, $tag, $max)=(@_);
    my ($n);
    $n=0;
    
    if ( !-e $file && ($file =~/$tag/)) {return 1;}
    elsif ( !-e $file){return 0;}
    else 
      {
	open (FC, "$file");
	while ( <FC>)
	  {
	    if ( ($_=~/$tag/))
	      {
		close (FC);
		return 1;
	      }
	    elsif ($max && $n>$max)
	      {
		close (FC);
		return 0;
	      }
	    $n++;
	  }
      }
    close (FC);
    return 0;
  }
	    
	  
sub file2string
  {
    my $f=@_[0];
    my $string, $l;
    open (F,"$f");
    while (<F>)
      {

	$l=$_;
	#chomp ($l);
	$string.=$l;
      }
    close (F);
    $string=~s/\r\n//g;
    $string=~s/\n//g;
    return $string;
  }


sub my_get_opt
  {
    my @list=@_;
    my $cl, $a, $argv, @argl;
    
    @argl=();
    $cl=shift @list;
    for ( $a=0; $a<=$#list; $a+=3)
      {
	$option=$list[$a];
	$optional=$list[$a+1];
	$status=$list[$a+2];
	$argv="";
	if ($cl=~/$option(\S+)/){$argv=$1;}
	@argl=(@argl,$argv);
	
	
	#$optional:0=>optional
	#$optional:1=>must be set
	#$status: 0=>no requirement
	#$status: 1=>must be an existing file
	#$status: 2=>must be an installed package
	

	if ($optional==0){;}
	elsif ( $optional==1 && $argv eq "")
	  {
	    print STDERR "ERROR: Option $option must be set [FATAL:$program/$mode/$method]\n";
	    myexit ($EXIT_FAILURE);
	  }
	if ($status==0){;}
	elsif ($status ==1 && $argv ne "" && !-e $argv)
	  {
	    print STDERR "ERROR: File $argv must exist [FATAL:$program/$mode/$method]\n";
	    myexit ($EXIT_FAILURE);
	  }
	elsif ( $status==2 && $argv ne "" && &check_pg_is_installed ($argv)==0)
	  {
	    print STDERR "ERROR: $argv is not installed [FATAL:$program/$mode/$method]\n";
	    myexit ($EXIT_FAILURE);
	  }
      }

    return @argl;
    }

sub tag2value 
  {
    
    my $tag=(@_[0]);
    my $word=(@_[1]);
    my $return;
    
    $tag=~/$word="([^"]+)"/;
    $return=$1;
    return $return;
  }
      
sub hit_tag2pdbid
  {
    my $tag=(@_[0]);
    my $pdbid;
       
    $tag=~/id="(\S+)"/;
    $pdbid=$1;
    $pdbid=~s/_//;
    return $pdbid;
  }
sub id2pdbid 
  {
    my $in=@_[0];
    my $id;
    
    $in=~/(\S+)/;
    $id=$in;
    
    if ($id =~/pdb/)
      {
	$id=~/pdb(.*)/;
	$id=$1;
      }
    $id=~s/[|��_]//g;
    return $id;
  }
sub set_blast_type 
  {
    my $file =@_[0];
    if (&file_contains ($file,"EBIApplicationResult",100)){$BLAST_TYPE="EBI";}
    elsif (&file_contains ($file,"NCBI_BlastOutput",100)) {$BLAST_TYPE="NCBI";}
    else
      {
	$BLAST_TYPE="";
      }
    return $BLAST_TYPE;
  }
sub blast_xml2profile 
  {
    my ($name,$seq,$maxid, $minid, $mincov, $file)=(@_);
    my (%p, $a, $string, $n);
    


    if ($BLAST_TYPE eq "EBI" || &file_contains ($file,"EBIApplicationResult",100)){%p=ebi_blast_xml2profile(@_);}
    elsif ($BLAST_TYPE eq "NCBI" || &file_contains ($file,"NCBI_BlastOutput",100)){%p=ncbi_blast_xml2profile(@_);}
    else 
      {
	print "************ ERROR: Blast Returned an unknown XML Format **********************";
	myexit ($EXIT_FAILURE);
      }
    for ($a=0; $a<$p{n}; $a++)
      {
	my $name=$p{$a}{name};
	$p{$name}{seq}=$p{$a}{seq};
      }
    return %p;
  }
sub ncbi_blast_xml2profile 
  {
    my ($name,$seq,$maxid, $minid, $mincov, $string)=(@_);
    my ($L,$l, $a,$b,$c,$d,$nhits,@identifyerL);
    
    
    $seq=~s/[^a-zA-Z]//g;
    $L=length ($seq);
    
    %hit=&xml2tag_list ($string, "Hit");
    
    
    for ($nhits=0,$a=0; $a<$hit{n}; $a++)
      {
	my ($ldb,$id, $identity, $expectation, $start, $end, $coverage, $r);
	my (%ID,%DE,%HSP);
	
	$ldb="";

	%ID=&xml2tag_list ($hit{$a}{body}, "Hit_id");
	$identifyer=$ID{0}{body};
	
	%DE=&xml2tag_list ($hit{$a}{body}, "Hit_def");
	$definition=$DE{0}{body};
	
	%HSP=&xml2tag_list ($hit{$a}{body}, "Hsp");
	for ($b=0; $b<$HSP{n}; $b++)
	  {
	    my (%START,%END,%E,%I,%Q,%M);

	 
	    %START=&xml2tag_list ($HSP{$b}{body}, "Hsp_query-from");
	    %HSTART=&xml2tag_list ($HSP{$b}{body}, "Hsp_hit-from");
	    
	    %LEN=  &xml2tag_list ($HSP{$b}{body}, "Hsp_align-len");
	    %END=  &xml2tag_list ($HSP{$b}{body}, "Hsp_query-to");
	    %HEND=  &xml2tag_list ($HSP{$b}{body}, "Hsp_hit-to");
	    %E=&xml2tag_list     ($HSP{$b}{body}, "Hsp_evalue");
	    %I=&xml2tag_list     ($HSP{$b}{body}, "Hsp_identity");
	    %Q=&xml2tag_list     ($HSP{$b}{body}, "Hsp_qseq");
	    %M=&xml2tag_list     ($HSP{$b}{body}, "Hsp_hseq");
	    
	    for ($e=0; $e<$Q{n}; $e++)

	      {
		$qs=$Q{$e}{body};
		$ms=$M{$e}{body};
		
		$expectation=$E{$e}{body};
		$identity=($LEN{$e}{body}==0)?0:$I{$e}{body}/$LEN{$e}{body}*100;
		$start=$START{$e}{body};
		$end=$END{$e}{body};
		$Hstart=$HSTART{$e}{body};
		$Hend=$HEND{$e}{body};
	
		$coverage=(($end-$start)*100)/$L;

	
		if ($identity>$maxid || $identity<$minid || $coverage<$mincov){next;}
		@lr1=(split (//,$qs));
		@lr2=(split (//,$ms));
		$l=$#lr1+1;
		for ($c=0;$c<$L;$c++){$p[$nhits][$c]="-";}
		for ($d=0,$c=0; $c<$l; $c++)
		  {
		    $r=$lr1[$c];
		    if ( $r=~/[A-Za-z]/)
		      {
			
			$p[$nhits][$d + $start-1]=$lr2[$c];
			$d++;
		      }
		  }
		$Qseq[$nhits]=$qs;
		$Hseq[$nhits]=$ms;
		$QstartL[$nhits]=$start;
		$HstartL[$nhits]=$Hstart;
		$identityL[$nhits]=$identity;
		$endL[$nhits]=$end;
		$definitionL[$nhits]=$definition;
		$identifyerL[$nhits]=$identifyer;
		$comment[$nhits]="$ldb|$identifyer [Eval=$expectation][id=$identity%][start=$Hstart end=$Hend]";
		$nhits++;
	      }
	  }
      }
    
    $profile{n}=0;
    $profile{$profile{n}}{name}=$name;
    $profile{$profile{n}}{seq}=$seq;
    $profile {n}++;
    
    for ($a=0; $a<$nhits; $a++)
      {
	$n=$a+1;
	
	$profile{$n}{name}="$name\_$a";
	$profile{$n}{seq}="";
	$profile{$n}{Qseq}=$Qseq[$a];
	$profile{$n}{Hseq}=$Hseq[$a];
	$profile{$n}{Qstart}=$QstartL[$a];
	$profile{$n}{Hstart}=$HstartL[$a];
	$profile{$n}{identity}=$identityL[$a];
	$profile{$n}{definition}=$definitionL[$a];
	$profile{$n}{identifyer}=$identifyerL[$a];
	$profile{$n}{comment}=$comment[$a];
	for ($b=0; $b<$L; $b++)
	  {
	    if ($p[$a][$b])
	      {
		$profile{$n}{seq}.=$p[$a][$b];
	      }
	    else
	      {
		$profile{$n}{seq}.="-";
	      }
	  }
      }
    
    $profile{n}=$nhits+1;
    return %profile;
  }
sub ebi_blast_xml2profile 
  {
    my ($name,$seq,$maxid, $minid, $mincov, $string)=(@_);
    my ($L,$l, $a,$b,$c,$d,$nhits,@identifyerL,$identifyer);
    

    
    $seq=~s/[^a-zA-Z]//g;
    $L=length ($seq);
    %hit=&xml2tag_list ($string, "hit");
    
    for ($nhits=0,$a=0; $a<$hit{n}; $a++)
      {
	my ($ldb,$id, $identity, $expectation, $start, $end, $coverage, $r);
	my (%Q,%M,%E,%I);
	
	$ldb=&tag2value ($hit{$a}{open}, "database");
	$identifyer=&tag2value ($hit{$a}{open}, "id");

	$description=&tag2value ($hit{$a}{open}, "description");
	
	%Q=&xml2tag_list ($hit{$a}{body}, "querySeq");
	%M=&xml2tag_list ($hit{$a}{body}, "matchSeq");
	%E=&xml2tag_list ($hit{$a}{body}, "expectation");
	%I=&xml2tag_list ($hit{$a}{body}, "identity");
	

	for ($b=0; $b<$Q{n}; $b++)
	  {

	    $qs=$Q{$b}{body};
	    $ms=$M{$b}{body};
	    
	    $expectation=$E{$b}{body};
	    $identity=$I{$b}{body};
	    
	    	    
	    $start=&tag2value ($Q{$b}{open}, "start");
	    $end=&tag2value ($Q{$b}{open}, "end");
	    $startM=&tag2value ($M{$b}{open}, "start");
	    $endM=&tag2value ($M{$b}{open}, "end");
	    $coverage=(($end-$start)*100)/$L;
	    
	   # print "$id: ID: $identity COV: $coverage [$start $end]\n";
	    
	    
	    if ($identity>$maxid || $identity<$minid || $coverage<$mincov){next;}
	    # print "KEEP\n";

	    
	    @lr1=(split (//,$qs));
	    @lr2=(split (//,$ms));
	    $l=$#lr1+1;
	    for ($c=0;$c<$L;$c++){$p[$nhits][$c]="-";}
	    for ($d=0,$c=0; $c<$l; $c++)
	      {
		$r=$lr1[$c];
		if ( $r=~/[A-Za-z]/)
		  {
		    
		    $p[$nhits][$d + $start-1]=$lr2[$c];
		    $d++;
		  }
	      }
	  
	    
	    $identifyerL[$nhits]=$identifyer;
	    $comment[$nhits]="$ldb|$identifyer [Eval=$expectation][id=$identity%][start=$startM end=$endM]";
	    $nhits++;
	  }
      }
    
    $profile{n}=0;
    $profile{$profile{n}}{name}=$name;
    $profile{$profile{n}}{seq}=$seq;
    $profile {n}++;
    
    for ($a=0; $a<$nhits; $a++)
      {
	$n=$a+1;
	$profile{$n}{name}="$name\_$a";
	$profile{$n}{seq}="";
	$profile{$n}{identifyer}=$identifyerL[$a];
	
	$profile{$n}{comment}=$comment[$a];
	for ($b=0; $b<$L; $b++)
	  {
	    if ($p[$a][$b])
	      {
		$profile{$n}{seq}.=$p[$a][$b];
	      }
	    else
	      {
		$profile{$n}{seq}.="-";
	      }
	  }
      }
    $profile{n}=$nhits+1;
    
    return %profile;
  }
sub output_profile
  {
    my ($name,%profile)=(@_);
    my ($a);
    open (P, ">$name");
    for ($a=0; $a<$profile{n}; $a++)
      {
	print P ">$profile{$a}{name} $profile{$a}{comment}\n$profile{$a}{seq}\n";
      }
    close (P);
    return;
  }
sub blast_xml2hit_list
  {
    my $string=(@_[0]);
    return &xml2tag_list ($string, "hit");
  }
sub xml2tag_list  
  {
    my ($string_in,$tag)=@_;
    my $tag_in, $tag_out;
    my %tag;
    
    if (-e $string_in)
      {
	$string=&file2string ($string_in);
      }
    else
      {
	$string=$string_in;
      }
    $tag_in1="<$tag ";
    $tag_in2="<$tag>";
    $tag_out="/$tag>";
    $string=~s/>/>##1/g;
    $string=~s/</##2</g;
    $string=~s/##1/<#/g;
    $string=~s/##2/#>/g;
    @l=($string=~/(\<[^>]+\>)/g);
    $tag{n}=0;
    $in=0;$n=-1;
  
 

    foreach $t (@l)
      {

	$t=~s/<#//;
	$t=~s/#>//;
	
	if ( $t=~/$tag_in1/ || $t=~/$tag_in2/)
	  {
	 
	    $in=1;
	    $tag{$tag{n}}{open}=$t;
	    $n++;
	    
	  }
	elsif ($t=~/$tag_out/)
	  {
	    

	    $tag{$tag{n}}{close}=$t;
	    $tag{n}++;
	    $in=0;
	  }
	elsif ($in)
	  {
	   
	    $tag{$tag{n}}{body}.=$t;
	  }
      }
  
    return %tag;
  }


sub seq2gor_prediction 
  {
    my ($name, $seq,$infile, $outfile, $gor_seq, $gor_obs)=(@_);
    my ($l);
    
    `gorIV -prd $infile -seq $gor_seq -obs $gor_obs > gor_tmp`;
    open (GR, ">$outfile");
    open (OG, "gor_tmp");

    while (<OG>)
      {
	
	$l=$_;
	if ($l=~/\>/){print GR "$l";}
	elsif ( $l=~/Predicted Sec. Struct./)
	  {
	    $l=~s/Predicted Sec. Struct\.//;
	    print GR "$l";
	  }
      }
    close (GR);
    close (OG);
    return;
  }
sub seq2msa_tm_prediction 
  {
    my ($name, $seq,$infile, $outfile, $arch, $psv)=(@_);
    my (%p,%gseq,%R, $blast_output, %s, $l);
    
    $blast_output=&run_blast ($name,"blastp", "uniprot", $infile, "outfile");
    
    
    %p=blast_xml2profile($name,$seq,$maxid, $minid,$mincov,$blast_output);
    
    
    open (F, ">tm_input");
    for ($a=0; $a<$p{n}; $a++)
      {
	my $s;
	
	$s=$p{$a}{seq};
	$s=uc($s);
	print F ">$p{$a}{name}\n$s\n";
	#print stdout ">$p{$a}{name}\n$s\n";
      }
    close (F);
    print "\tPSITM: kept  $p{n} Homologues for Sequence $p{0}{name}\n";
    &safe_system ("t_coffee -other_pg fasta_seq2hmmtop_fasta.pl -in=tm_input -out=tm_output -arch=$arch -psv=$psv");
    unlink ("tm_input");
    %gs=read_fasta_seq("tm_output");
    foreach $s (keys(%gs))
      {
	my (@list, $seq, @plist, @pseq, $L, $PL);
	
	
	#Prediction
	$seq=$gs{$s}{seq};
	$seq=uc($seq);
	$L=length($seq);
	@list=split //, $seq;
	
	#Original Profile Sequence
	$pseq=$p{$s}{seq};
	$pseq=uc($pseq);
	$PL=length($pseq);
	@plist=split //, $pseq;
	
	for ($c=0,$b=0; $b<$PL; $b++)
	  {
	    my $r=$plist[$b];
	    if($r ne "-" && $r ne "X")
	      {
		$r=$plist[$b]=$list[$c++];
	      }
	  }
	
	if ($c!=$L)
	  {
	    print "ERROR: Could Not Thread the Prediction Back [FATAL:tc_generic_method.pl]\n";
	    myexit ($EXIT_FAILURE);
	  }
	for ($b=0;$b<$PL; $b++)
	  {
	    my $r=$plist[$b];
	    if ( $r ne "-" && $r ne "X")
	      {
		$R{$b}{$r}++;
	     }
	  }
      }
    $L=length ($p{0}{seq});
    open (R2, ">$outfile");
    print R2 ">$name\n";
    
    for ($a=0; $a<$L; $a++)
      {
	
	my ($v,$v_max,$r,$r_max, @rl);
	
	$v=$v_max=0;
	@rl=keys (%{$R{$a}});
	foreach $r (@rl)
	  {

	    $v=$R{$a}{$r};
	    if ($v>=$v_max)
	      {
		$v_max=$v;
		$r_max=$r;
	      }
	  }
	print R2 "$r_max";
      }
    print R2 "\n";
    close (R2);
    return;
  }
sub seq2msa_gor_prediction 
  {
    my ($name, $seq,$infile, $outfile, $gor_seq, $gor_obs)=(@_);
    my (%p,%gseq,%R, $blast_output, %s, $l);
    
    
    $blast_output=&run_blast ($name,"blastp", "uniprot", $infile, "outfile");
    %p=blast_xml2profile($name,$seq,$maxid, $minid,$mincov,$blast_output);
    
    open (F, ">gor_input");
    for ($a=0; $a<$p{n}; $a++)
      {
	my $s;
	
	$s=$p{$a}{seq};
	$s=~s/\-//g;
	$s=~s/X//g;
	
	$s=uc($s);
	print F ">$p{$a}{name}\n$s\n";
      }
    close (F);
    print "\tPSIGOR: kept  $p{n} Homologues for Sequence $p{0}{name}\n";
    
    `gorIV -prd gor_input -seq $gor_seq -obs $gor_obs > gor_tmp`;
    unlink ("gor_input");
    
    open (GR, ">gor_output");
    open (OG, "gor_tmp");
    
    while (<OG>)
      {
	
	my $l;
	$l=$_;
	
	if ($l=~/\>/){print GR "$l";}
	elsif ( $l=~/Predicted Sec. Struct./)
	  {
	    $l=~s/Predicted Sec. Struct\.//;
	    print GR "$l";
	  }
      }
    close (GR);
    close (OG);
    

    %gs=read_fasta_seq("gor_output");
     foreach $s (keys(%gs))
      {
	my (@list, $seq, @plist, @pseq, $L, $PL);
	
	
	#Prediction
	$seq=$gs{$s}{seq};
	$seq=uc($seq);
	$L=length($seq);
	@list=split //, $seq;
	
	#Original Profile Sequence
	$pseq=$p{$s}{seq};
	$pseq=uc($pseq);
	$PL=length($pseq);
	@plist=split //, $pseq;
	
	$tseq="";
	for ($c=0,$b=0; $b<$PL; $b++)
	  {
	    my $r=$plist[$b];
	    if($r ne "-" && $r ne "X")
	      {
		$r=$plist[$b]=$list[$c++];
	      }
	    $tseq.=$r;
	  }
	
	if ($c!=$L)
	  {
	    print "ERROR: Could Not Thread the Prediction Back [FATAL:tc_generic_method.pl]\n";
	    print "SEQ:$seq\nPSEQ:$pseq\nTSEQ:$tseq";
	    
	    myexit ($EXIT_FAILURE);
	  }
	for ($b=0;$b<$PL; $b++)
	  {
	    my $r=$plist[$b];
	    if ( $r ne "-" && $r ne "X")
	      {
		$R{$b}{$r}++;
	     }
	  }
      }
   
    $L=length ($p{0}{seq});
    open (R2, ">$outfile");
    print R2 ">$name\n";
    
    for ($a=0; $a<$L; $a++)
      {
	
	my ($v,$v_max,$r,$r_max, @rl);
	
	$v=$v_max=0;
	@rl=keys (%{$R{$a}});
	foreach $r (@rl)
	  {

	    $v=$R{$a}{$r};
	    if ($v>=$v_max)
	      {
		$v_max=$v;
		$r_max=$r;
	      }
	  }
	print R2 "$r_max";
      }
    print R2 "\n";
    close (R2);
    return;
  }

sub run_blast
  {
    my ($name, $method, $db,$infile, $outfile, $run)=(@_);
    if (!$run){$run=1;}
    
    
    if (&cache_file("GET",$infile,$name,$method,$db,$outfile,$SERVER)){return $outfile;}
    else
      {
	
	if ( $SERVER eq "EBI")
	  {
	    $cl_method=$method;
	    if ($cl_method =~/wu/)
	      {
		$cl_method=~s/wu//;
		if ( $cl_method eq "psiblast")
		  {
		    print STDERR "\n***************WARNING: PSI BLAST cannot be used with the NCBI BLAST Client. Use server=EBI Or server=LOCAL. blastp will be used instead***********\n";
		    $cl_method="blastp";
		  }
		
		$command="t_coffee -other_pg wublast.pl --email $EMAIL $infile -D $db -p $cl_method --outfile $outfile -o xml>/dev/null 2>/dev/null";
		&safe_system ( $command);
		if (-e "$outfile.xml") {`mv $outfile.xml $outfile`;}
	      }
	    else
	      {
		if ($cl_method eq "psiblast"){$cl_method ="blastp -j5";}
		
		$command="t_coffee -other_pg blastpgp.pl --email $EMAIL $infile -d $db --outfile $outfile -p $cl_method --mode PSI-Blast>/dev/null 2>/dev/null";
		&safe_system ( $command);
		
		if (-e "$outfile.xml") {`mv $outfile.xml $outfile`;}
	      }
	  }
	elsif ($SERVER eq "NCBI")
	  {
	    if ($db eq "uniprot"){$cl_db="nr";}
	    else {$cl_db=$db;}
	    
	    if ( $method eq "psiblast")
	      {
		print STDERR "\n***************WARNING: PSI BLAST cannot be used with the NCBI BLAST Client. Use server=EBI Or server=LOCAL. blastp will be used instead***********\n";
		$cl_method="blastp";
	      }
	    else
	      {
		$cl_method=$method;
	      }
	    $command="blastcl3 -p $cl_method -d $cl_db -i $infile -o $outfile -m 7";
	    &mysystem ($command);
	  }
	elsif ($SERVER =~/CLIENT_(.*)/)
	  {
	    my $client=$1;
	    $command="$client -p $method -d $db -i $infile -o $outfile -m 7";
	    &mysystem ($command);
	  }
	elsif ( $SERVER eq "LOCAL_blastall")
	  {
	    if ($method eq "blastp")
	      {
		$command="blastall -d $db -i $infile -o $outfile -m7 -p blastp";
	      }
	    &mysystem ($command);
	  }
	elsif ( $SERVER eq "LOCAL")
	  {

	    if ($ENV{"BLAST_DB_DIR"})
	      {
		$x=$ENV{"BLAST_DB_DIR"};
		$cl_db="$x$db";
	      }
	    else
	      {
		$cl_db=$db;
	      }
	    
	    if ($method eq "blastp")
	      {
		$command="blastpgp -d $cl_db -i $infile -o $outfile -m7 -j1";
	      }
	    elsif ($method eq "psiblast")
	      {
		$command="blastpgp -d $cl_db -i $infile -o $outfile -m7 -j5";
	      }
		elsif ($method eq "blastn")
		{
			$command="blastall -p blastn -d $cl_db -i $infile -o $outfile -m7 -W6";
		}	
	    &mysystem ($command);
	  }
	else
	  {
	    print ("*************** ERROR: $SERVER is an Unknown Server***********");
	  }
	
	if ( !-e $outfile)
	  {
	    
	    if ( $run==$BLAST_MAX_NRUNS)
	      {
		print STDERR "COM: $command\n";
		print STDERR ("BLAST failed against $name [FATAL:$mode/$method/$program]\n");
		if  ( $SERVER eq "EBI" && !($method=~/wu/))
		  {
		    print STDERR ("Try WuBlast instead");
		    return run_blast ($name,"wublastp", $db,$infile, $outfile);
		  }
	      }
	    else
	      {
		print STDERR "(Blast for $name failed  [$command][Attempt $run/$BLAST_MAX_NRUNS] [Try again]\n";
		return run_blast ($name, $method, $db,$infile, $outfile, $run+1);
	      }
	  }

	&cache_file("SET",$infile,$name,$method,$db,$outfile,$SERVER);
	return $outfile;
      }
  }
sub mysystem 
  {
    my $command=@_[0];
    my $count=0;
    my $r;
    
    while (($r=&safe_system($command))!=$EXIT_SUCCESS && $count<5)
      {
	print "\nCOMMAND $command Failed. Will try again\n";
	$count++;
      }
    return $r;
  }
sub cache_file
  {
    my ($cache_mode,$infile,$name,$method,$db, $outfile,$server)=(@_);
    my $cache_file;
    #Protect names so that they can be turned into legal filenames
    $name=&clean_file_name ($name);

    if ($db=~/\//)
      {
	$db=~/([^\/]+)$/;
	$db=$1;
      }
    $cache_file_sh="$name.$method.$db.$server.tmp";
    $cache_file="$CACHE/$name.$method.$db.$server.tmp";
    
    if ($infile ne "")
      {
	$cache_file_infile_sh="$name.$method.$db.$server.infile.tmp";
	$cache_file_infile="$CACHE/$name.$method.$db.$server.infile.tmp";
      }
    
    if ($cache_mode eq "GET")
      {
	if ($CACHE eq "" || $CACHE eq "no" || $CACHE eq "ignore"  || $CACHE eq "local" || $CACHE eq "update"){return 0;}
	elsif ( !-d $CACHE)
	  {
	    print STDERR "ERROR: Cache Dir: $CACHE Does not Exist";
	    return 0;
	  }
	else
	  {
	    if ( -e $cache_file && &fasta_file1_eq_fasta_file2($infile,$cache_file_infile)==1)
	      {
		`cp $cache_file $outfile`;
		$CACHE_STATUS="READ CACHE";
		return 1;
	      }
	  }
      }
    elsif ($cache_mode eq "SET")
      {
	if ($CACHE eq "" || $CACHE eq "no" || $CACHE eq "ignore"  || $CACHE eq "local" || $CACHE eq "update"){return 0;}
	elsif ( !-d $CACHE)
	  {
	    print STDERR "ERROR: Cache Dir: $CACHE Does not Exist";
	    return 0;
	  }
	elsif (-e $outfile)
	  {
	    `cp $outfile $cache_file`;
	    if ($cache_file_infile ne ""){ `cp $infile $cache_file_infile`;}

	    #functions for updating the cache
	    #`t_coffee -other_pg clean_cache.pl -file $cache_file_sh -dir $CACHE`;
	    #`t_coffee -other_pg clean_cache.pl -file $cache_file_infile_sh -dir $CACHE`;
	    return 1;
	  }
      }
    $CACHE_STATUS="COMPUTE CACHE";
    return 0;
  }
sub file1_eq_file2
  {
    my ($f1, $f2)=@_;
    if ( $f1 eq ""){return 1;}
    elsif ( $f2 eq ""){return 1;}
    elsif ( !-e $f1){return 0;}
    elsif ( !-e $f2){return 0;}
    elsif ($f1 eq "" || $f2 eq "" || `diff $f1 $f2` eq ""){return 1;}
    
    return 0;
  }
sub clean_file_name 
  {
    my $name=@_[0];
    
    $name=~s/[^A-Za-z1-9.-]/_/g;
    return $name;
  }
sub url2file
  {
    my ($address, $out)=(@_);
    
    if (&pg_is_installed ("wget"))
	{
	  return &safe_system ("wget $address -O$out >/dev/null 2>/dev/null");
	}
    elsif (&pg_is_installed ("curl"))
      {
	return &safe_system ("curl $address -o$out >/dev/null 2>/dev/null");
      }
    else
      {
	print stderr "ERROR: neither curl nor wget are installed. Imnpossible to fectch remote file [FATAL]\n";
	exit ($EXIT_FAILURE);
      }
  }
sub fasta_file1_eq_fasta_file2
  {
    my ($f1, $f2)=@_;
    my (%s1, %s2);
    my @names;
    %s1=read_fasta_seq (%f1);
    %s2=read_fasta_seq (%f2);

    @names=(keys (%s1));
    
    foreach $n (keys(%s1))
      {
	if ($s1{$n}{seq} ne $s2{$n}{seq}){return 0;}
      } 
    
    foreach $n (keys(%s2))
      {
	if ($s1{$n}{seq} ne $s2{$n}{seq}){return 0;}
      }
    return 1;
  }
	
sub safe_system 
{
  my $com=@_[0];
  my $pid;
  my $status;
  if ($com eq ""){return 1;}


  if (($pid = fork ()) < 0){return (-1);}
  if ($pid == 0)
    {
      exec ($com);
    }
  else
      {
	$PIDCHILD=$pid;
      }
  
  waitpid ($pid,WTERMSIG);
  return $?; #contains the status of the exit
}
END {
  kill ($PIDCHILD);
}


sub read_template_file
{
	my $pdb_templates = @_[0];
	open (TEMP, "<$pdb_templates");
	my %temp_h;
	while (<TEMP>)
{
		$line = $_;
 		$line =~/(\S+)\s(\S+)/;
 		$temp_h{$1}= $2;
}
	close(TEMP);
	return %temp_h;
}

sub calc_rna_template
{
	my ($mode, $infile, $pdbfile, $outfile)=@_;
	my %s, %h ;
	my $result;
	my (@profiles);
	&set_temporary_dir ("set",$infile,"seq.pep");
	%s=read_fasta_seq ("seq.pep");
	
	%pdb_template_h = &read_template_file($pdbfile);
	my $pdb_chain;
	open (R, ">result.aln");


	#print stdout "\n";
	foreach $seq (keys(%s))
	{
		if ($pdb_template_h{$seq} eq "")
		{
			next;
		}
		open (F, ">seqfile");
		print (F ">$s{$seq}{name}\n$s{$seq}{seq}\n");
		close (F);
		$pdb_chain = $pdb_template_h{$seq};
		$lib_name="$s{$seq}{name}.rfold";
		$lib_name=&clean_file_name ($lib_name);
		safe_system ("t_coffee -other_pg RNAplfold2tclib.pl -in=seqfile -out=$lib_name");
		
 		safe_system ("secondary_struc.py seqfile $CACHE$pdb_chain  $lib_name");
		
		if ( !-e $lib_name)
		{
			print STDERR ("RNAplfold failed to compute the secondary structure of $s{$seq}{name} [FATAL:$mode/$method/$program]\n");
			myexit ($EXIT_FAILURE);
		}
		else
		{
			print stdout "\tProcess: >$s{$seq}{name} _F_ $lib_name\n";
			print R ">$s{$seq}{name} _F_ $lib_name\n";
		}
		unshift (@profiles, $lib_name);
	}
	close (R);
	&set_temporary_dir ("unset",$mode, $method,"result.aln",$outfile, @profiles);
}



sub seq2rna_pair{
	my ($mode, $pdbfile1, $pdbfile2, $method, $param, $outfile)=@_;

	if ($method eq "runsara.py")
	{
		open(TMP,"<$pdbfile1");
		my $count = 0;
		my $line;
		while (<TMP>)
		{
			$line = $_;
			if ($count ==1)
			{
				last;
			}
			$count += 1;
		}
	
		my $y = length($line);

		$chain1 = substr($line,length($line)-3,1);
		close TMP;
		open(TMP,"<$pdbfile2");
		my $count = 0;
		while (<TMP>)
		{
			$line = $_;
			if ($count ==1)
			{
				last;
			}
			$count += 1;
		}
		$chain2 = substr($line,length($line)-3,1);
		close TMP;
		
		
		system("runsara.py $pdbfile1 $chain1 $pdbfile2 $chain2 -s -o tmp >/dev/null 2>/dev/null");
		open(TMP,"<tmp") or die "cannot open the sara tmp file:$!\n";
		open(OUT,">$outfile") or die "cannot open the $outfile file:$!\n";

		my $switch = 0;
		my $seqNum = 0;
		foreach my $line (<TMP>)
		{
			next unless ($line=~/SARAALI/);
			if ($line=~/>/)
			{
				$switch =0;
				print OUT ">seq$seqNum\n";
				$seqNum++;				
			}
			if ($switch < 2){
				$switch++;
				next;
			}
	
			if ($line =~/REMARK\s+SARAALI\s+([^\*]+)\*/)
			{
				my $string = $1;
				print OUT "$string\n";
			}
		}
		close TMP; 
		close OUT;
	}
}$program="T-COFFEE (Version_8.07)";\n

