#!/usr/bin/perl
#
#
#date : 12/10/
#PArameters handling;
$program="blastp";


foreach $v(@ARGV){$cl.="$v ";}
($mode)=&my_get_opt ( $cl, "-program=",0,0);
($db)=&my_get_opt ( $cl, "-database=",0,0);
($infile)=&my_get_opt ( $cl, "-infile=",1,1);
($outfile)=&my_get_opt ( $cl, "-outfile=",1,0);




 usage: web_blast.pl -infile <fasta file> -method <pdbid/geneid or profile> options []


       -program ...... Program Name (blastp)
                       Default = blastp
       -database ..... Database at NCBI (nr, pdb, swissprot,refseq_protein) or indicate a local fasta file
                       Default = pdb at NCBI
       -infile ....... Query_file = a list of sequences in fasta format
       -outfile ...... Name the outfile to make a template file for t_coffee
                       Default = STDOUT or  default.profile if method is profile
       -evalue ....... Evalue threshold Default = 1;
       -matrix ....... PAM30 PAM70 BLOSUM45 BLOSUM80
                       Default BLOSUM62
       -method ....... geneid,pdbid, profile
       -gigablast..... yes/no FASTER REMOTE BLAST with Gigablaster (Stephane Audic program: http://www.igs.cnrs-mrs.fr/adele/~database/remoteblast.c
                       Default no
       -filter ....... T or F locally, L or R or M or C or V for distant blast
                       Default = Off
       -organism ..... Fungi, Bos_taurus, Vertebrata, Mus_musculus, Bacteria, Gallus_gallus, Caenorhabditis_elegans, Primates, Escherichia_coli, Ara
, Drosophila_melanogaster, Eukaryota, Mammalia, Viruses, Archaea  are available
                       Default is All_organisms
       -identity ..... blast identity threshold = provide a % for view only the results upper or equal to the threshold
                       Default 50
       -cover ........ Cover threshold = provide a % : sequence covering Default: 30
       -hits ........  Number of hits
                       Default = 1
       -processor .... Number of processors to use
                       Default = 1
       -blast_dir .... Indicates where your BLAST directory  is installed localy
       -quiet ........ on : do not display all the default/defined blast parameters
								    Default off
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
	    exit (EXIT_FAILURE);
	  }
	if ($status==0){;}
	elsif ($status ==1 && $argv ne "" && !-e $argv)
	  {
	    print STDERR "ERROR: File $argv must exist [FATAL:$program/$mode/$method]\n";
	    exit (EXIT_FAILURE);
	  }
	elsif ( $status==2 && $argv ne "" && &check_pg_is_installed ($argv)==0)
	  {
	    print STDERR "ERROR: $argv is not installed [FATAL:$program/$mode/$method]\n";
	    exit (EXIT_FAILURE);
	  }
      }

    return @argl;
    }
