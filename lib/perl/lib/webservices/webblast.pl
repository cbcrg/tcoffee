#!/usr/bin/env perl
#
#
#date : 27/09/05
#prog : webblast.pl
#subj : make a BLAST/WU-BLAST (by HTTP request or locally) against a database with a file containing sequences in fasta format
####### method genid, pdbid and profile 
#
#############################################################################################
use Env qw(HOME);
use lib "$HOME/.lib_webblast/";
use LWP::UserAgent; 
use HTML::Parser;                                            # @@@@@@@ #       
use HTTP::Request::Common qw(POST);                         # @/^   ^\@ # 
use URI::Escape;                                           # @/ -   - \@ #
use Getopt::Long;                                         ##  \   ^   /  ##
use strict;                                              ##    |  0  |    ##
use warnings;                                           ####### \ _ / #######
##############################################################################################


############################  EXPRESSO PARAM  ################################
my $database_expresso="pdb_seqres.txt"; #PDB database name                   #
my $blast_dir_expresso="/home/igs/public_html/Tcoffee/Exe/blastall";         #############
my $BLASTMAT="export BLASTMAT=/home/igs/Src/blast/data/"; #MAtrix directory for blastall #
my $BLASTDB="export BLASTDB=/banks/_Proteins/PDB_seqres/";#PDB directoty     #############
##############################################################################
    
my $runblast="runblast.pl";

my(@list_encoded)=(), my(@list_pdb)=(), my(%deja_vu)=(), my(@pdb_list)=(),my($i)=0, my(@names)=(), my $locale=0, my $distant=0, my $database, my $blast_way;
my($ua)= LWP::UserAgent->new;


##-- Variables d'environnements 

my($database_var)= $ENV { 'DATABASE' };
my($blast_var)   = $ENV { 'BLAST_DIRECTORY' };

##-- Recupere Options/parametres du BLAST && controle des options ds OPTIONS_GET
 
my($program,$database_line,$blast_line,$query_file,
   $out_file,$identity_treshold,$cover_tresh,$Eval,
   $align,$matrix,$filter,$method,$orgn,$process,$quiet,$gigablast)= &OPTIONS_GET();

##-- Determination BLAST LOCAL /DISTANT && Controle database/programme

unless (-e $query_file ) { print STDERR "\nfile does not exist!\n";exit;}
unless (-s $query_file ) { print STDERR "\nyour file is empty!\n";exit; }

if ((($database_line || $database_var) && ($blast_line || $blast_var)) ) 
{   
    if ($database_line=~/expressopdb/ && $blast_line=~/blastexpresso/)
    {
	#mode special pour fichier de configuration du serveur Expresso
	$locale=1;
	$database="expressopdb";
	unless ($quiet=~ /on/i) { print STDERR "\nRUN BLAST LOCALY\n"; }
    }
    else
    {
	($database)=$database_line || $database_var;	
	my($blast_tp)=$blast_var || $blast_line; $locale=1;
	$blast_way = &CONTROLE_DB_PG($database,$blast_tp,$program);
	unless ($quiet=~ /on/i) { print STDERR "\nRUN BLAST LOCALY\n"; }
    }
}

else
{
    $database = &NCBI_DATABASE($database_line); $distant=1;    
    if ($gigablast=~ /^yes$/i) { $locale=2; $distant=0; unless ($quiet=~ /on/i) { print STDERR "\nRUN GIGABLASTER\n"; }}   
    else { unless ($quiet=~ /on/i)  {print STDERR "\nRUN BLAST AT THE NCBI\n";}} 
}
 
##-fixation de parametres selon la valeur du flag -method

if ($method=~ /^pdbid$/i)
{   
  
    if    ($gigablast=~ /^yes$/i)          
    {
	if ($database ne "pdb") { print STDERR "\nprovide a valid database name FOR RUN GIGABLASTER: nr,pdb or refseq_protein\n";exit;} else {$database="pdbaa";}
    }
    elsif ($gigablast=~ /^no$/i && $distant==1) { $database="pdb";}
          
}
elsif ($method=~ /^geneid$/i)
{
    if ($distant==1) { unless ($database eq "nr" || $database eq "swissprot" || $database eq "pdb") { $database="refseq_protein";}}   
}
elsif ($method=~/^profile$/i)
{
    if ($distant==1) {	unless ($database eq "nr" || $database eq "swissprot" || $database eq "refseq_protein" )  { $database="pdb"; } }
}
else { die "unknown method\n";}

if (($orgn !~ /All\+organisms/) && ($locale=~/1|2/))
   { print STDERR "-organism option can't be used locally or with -gigablast option!\n";exit;}

##---AFFICHAGE des valeurs des options
unless ($quiet =~ /on/i )
{
    print STDERR "   
              
             Program : $program
             Database : $database
             Method : $method
     
             Query_file : $query_file
             Out_file : $out_file 
             ";
  print STDOUT "      
             Evalue threshold : $Eval
             Matrix : $matrix
             Filter : $filter
             Blast_identity_threshold : $identity_treshold
             Cover threshold : $cover_tresh
              ";
    print STDERR "
             Number of hits :  $align
             Number of processors used : $process 
             ";
    if ($gigablast=~ /^yes$/i) { print STDERR "
             gigablast: yes\n" }
    unless ($locale) { print STDERR "
             Organism : $orgn\n" }
            
                        
print STDERR "
***************************************************************\n\n";
}


#-- LOCAL/DISTANT BLASTP
   
if    ($locale==1 || $locale==2) {@list_pdb= &LOCAL_BLAST ($blast_way,$database,$query_file,$Eval,$align,$method,$matrix,$filter,$process,$gigablast,$database_expresso,$blast_dir_expresso,$runblast);   }
elsif ($distant==1)  { @list_pdb= &WEB_BLAST   ($query_file,$Eval,$program,$database,$matrix,$method,$align,$orgn,$filter); }
else              { die " Report bug to armougom\@igs.cnrs-mrs.fr\n"}; 

#-- PARSE BLAST RESULTS -> MAKE A PDB_ID LIST
if ($method =~ /^pdbid$/i) 
{ 
    my(@result_sort)= &PARSING (\@list_pdb,$locale,$distant,$method,$quiet,$database,$gigablast);
   
                      &AFFICHAGE_PDB_PARSING (\@result_sort,$cover_tresh,$identity_treshold,$out_file);
                      exit;
}

#-- PARSE BLAST RESULT -> MAKE LIST OF REFSEQ ID
elsif ($method =~ /^geneid$/i)
{ 
    my(@result_sort)= &PARSING (\@list_pdb,$locale,$distant,$method,$quiet,$database,$gigablast);
                      &AFFICHAGE_REFSEQ_PARSING (\@result_sort,$cover_tresh,$identity_treshold,$out_file);
                      exit;
}

#-- PARSE BLAST RESULT -> MAKE PROFILE
elsif ($method=~ /^profile$/i) { &PROFILE (\@list_pdb,$out_file,$distant); exit; }
else { die " \nFATAL ERROR :  Method or database error\n" ;}

exit;
                   
                                               ##############
###############################################  FONCTIONS  ####################################################################
                                              ##############
sub CONTROLE_DB_PG
{
    my($database,$blast_dir,$program)=@_;
    
    if (! -e $database) { die "$database file  does not exist\n";}
    if ( -d $database) { die "$database must be a file, not a directory\n"; }    
    if ($blast_dir !~ /\/$/) { $blast_dir.="/"; }    
    my ($blastall) = $blast_dir . "blastall";
    if (! -e $blastall) { die "$blastall program not found \n";} 

    return ($blastall);
}

#-------------------------------------------------------------------------------------------------------------------------------------
sub NCBI_DATABASE
{
    my($ncbi_db)=@_;
    
    my  (%all_db)=
	(
	         'nr'                 =>'1',
	         'pdb'               =>'1',
	         'swissprot'        =>'1',
                 'refseq_protein'  =>'1',
                	       
	);
    
    if (exists $all_db{$ncbi_db}) { return($ncbi_db); }    
    elsif ($ncbi_db eq "")        { return (""); }
    else                          { return (1);  }
} 
#------------------------------------------------------------------------------------------------------------------------
sub HELP
{   
    my($org,@orga)= &LIST_ORGA();
    my ($list_orga)=join(', ',@orga);

    print STDERR "
                      usage: web_blast.pl -infile <fasta file> -method <pdbid/geneid or profile> options []\n

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
                            -gigablast..... yes/no FASTER REMOTE BLAST with Gigablaster (Stephane Audic program: http://www.igs.cnrs-mrs.fr/adele/~database/remoteblast.cgi)  
                                            Default no
                            -filter ....... T or F locally, L or R or M or C or V for distant blast
                                            Default = Off
                            -organism ..... $list_orga  are available 
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

                     Environement Variables
                     These variables can be set from the environement
           DATABASE......................[Indicates where your database file must be fetched (localy)]
           BLAST_DIRECTORY...............[Indicates where your BLAST directory is installed localy]
                           
";
	
    exit;
    
}
#-----------------------------------------------------------------------------------------
sub OPTIONS_GET
{   
    my %opt=();
 
    GetOptions 
              (
	       "infile=s"    =>\$opt{infile},
	       "outfile=s"    =>\$opt{outfile},
	       "program=s"     =>\$opt{program},
	       "database=s"     =>\$opt{database},
	       "blast_dir=s"     =>\$opt{blast_dir},
	       "identity=f"       =>\$opt{treshold},
	       "cover=f"           =>\$opt{cover},
	       "evalue=f"             =>\$opt{evalue},     
	       "hits=i"                 =>\$opt{hits},	   
	       "matrix=s"                 =>\$opt{matrix},
	       "filter=s"                  =>\$opt{filter},
	       "method=s"                   =>\$opt{method},
               "organism=s"                  =>\$opt{organism},
	       "processor=i"                  =>\$opt{processor},
	       "quiet=s"                       =>\$opt{quiet},
	       "gigablast=s"                    =>\$opt{gigablast},
	       );
  
    if ($ARGV[0]) {print "Unprocessed by Getopt::Long\n $ARGV[0]\n"; &HELP();} 
   
  

    my($evalue_tresh)= $opt{'evalue'};       unless ($evalue_tresh) { $evalue_tresh=1;};
    my($cover_tresh) = $opt{'cover'};         unless (defined $cover_tresh)  { $cover_tresh=30;};
    my($query_file)  = $opt{'infile'};         unless ($query_file)   { print STDERR "Flag -infile must be defined\n"; &HELP();};
    my($outfil)      = $opt{'outfile'};         unless ($outfil)       { $outfil="";};
    my($treshold)    = $opt{'treshold'};         unless (defined $treshold)     { $treshold=50;};
    my($blast_dir)   = $opt{'blast_dir'};         unless ($blast_dir)    { $blast_dir="";};
    my($database)    = $opt{'database'};           unless ($database)     { $database="";};
    my($program)     = $opt{'program'};             unless ($program)      { $program="blastp";};  
    my($align)       = $opt{'hits'};                  unless (defined $align)       { $align=1;};  
    my($matrix)      = $opt{'matrix'};                  unless ($matrix)      { $matrix="BLOSUM62";};
    my($filter)       = $opt{'filter'};                  unless ($filter)      { $filter="F";};
    my($method)      = $opt{'method'};                     unless ($method)      {print STDERR "Flag -method must be defined\n"; &HELP();};
    my($organism)    = $opt{'organism'};                    unless ($organism)    { $organism="All organisms"};
    my($process)     = $opt{'processor'};                    unless ($process)     { $process="1" };   
    my($param)       = $opt{'quiet'};                         unless ($param)       { $param= "off" }; 
    my($gigablast)   = $opt{'gigablast'};                      unless ($gigablast)    {$gigablast="no"}; 
    if ($method !~ /(^geneid$|^pdbid$|^profile$)/i)  { print STDERR "unknown method for the flag -method\n";&HELP(); }

	if ($treshold <0 || $treshold >100)               { print STDERR "\nout of range for the option -treshold \n"; &HELP();}  
	if ($cover_tresh <0 || $cover_tresh >100)         { print STDERR "\nout of range for the option -cover \n"; &HELP();} 
	if ($align <0)                                    { print STDERR "\n error with option   align\n"; &HELP();}
	if ($gigablast!~/^yes$|^no$/i)                    { print STDERR "invalid argument for gigaglast option : yes/no\n";exit;};
	if ($filter!~ /^[TFRLMCV]{1}$|^off$/i)                  {print STDERR  "valid values for -filter are T,F,R,L,M,C,or V!\n";exit;}
	if ($matrix!~ /PAM30|PAM70|BLOSUM45|BLOSUM80|BLOSUM62/) { print STDERR "valid values for -matrix  are PAM30,PAM70,BLOSUM45,BLOSUM80 or BLOSUM62\n";exit }
	if ($outfil eq "" && $method=~ /^profile$/i) { $outfil="default_profile.template"}
	
	if ($param!~ /^on$|^off$/i)                       { print STDERR "valid values for -quiet is on or off\n";exit;}
	my($orgn,@all_orgn)= &ORGN($organism);
	return ($program,$database,$blast_dir,$query_file,
		    $outfil,$treshold,$cover_tresh,$evalue_tresh,
		    $align,$matrix,$filter,$method,$orgn,$process,$param,$gigablast);
}

#--------------------------------------------------------------------------------------------------
sub RECOVER
{
    my($pdb_result,$aln_length,$length_query)=@_;
    my $nb_gap=0;
    
    if ($pdb_result=~ /(score.+?\n\n\n).+?score/ism) #cas ou plusieurs HSP, prend que le 1er
    { $pdb_result= $1;}
    
    $length_query  =~ s/,//g;
    my ($requete)  =  join('',($pdb_result=~/^Query(.*)\n/gm));
    $requete       =~ s/[^A-Z-]//g;
    my(@sequence)  =  split('',$requete);
 
    for (my $i=0; $i<=$#sequence; $i++)
    {
	if($sequence[$i] eq "-"){ ++$nb_gap; }
    }
    my($recouvrement)= sprintf("%-3d",(($aln_length-$nb_gap)/$length_query)*100);
    undef(@sequence);
    return ($recouvrement,$nb_gap); 
   
}
#--------------------------------------------------------------------------------------------------
sub WEB_BLAST
{
    open (SOR1,">web_tempo.result") or die;
    my($query_file,$Eval,$program,$database,$matrix,$method,$align,$orgn,$filter)=@_;
    my $aln_view, my $format="Txt";
    my($description)=$align;
    if ($method=~/^profile$/i) { $aln_view ="FlatQueryAnchoredNoIdentities"} else { $aln_view ="Pairwise"}
   
    if ($filter eq "F") { $filter="off";}

    $/=">";
    open(FIC,$query_file) or die "can not open $query_file $!\n";
    my(@sequences)=<FIC>;
    close FIC;
    shift(@sequences);

    foreach my $sequence(@sequences)	
    {
	$sequence=~ s/>//g;
	$sequence=">$sequence";
	
	my($name)=($sequence=~ /^>(.+)\n/);
	push(@names, $name); 
	my($encoded_query)= uri_escape($sequence);
	push (@list_encoded, $encoded_query);    
    } 
    
    undef(@sequences);
    if (scalar (@names != @list_encoded)) { die "error $!";}     
    foreach my $encoded_seq(@list_encoded)
    {    
	my $nb=0;
	print STDERR "BLAST $names[$i]...";
	
#-- BUILD THE REQUEST
		
	my($arguments) = "CMD=Put&ENTREZ_QUERY=$orgn&CDD_SEARCH=off&FILTER=$filter&MATRIX_NAME=$matrix&PROGRAM=$program&DATABASE=$database&QUERY=" . $encoded_seq;
	
	my($req) = new HTTP::Request POST => 'http://www.ncbi.nlm.nih.gov/blast/Blast.cgi';
	$req -> content_type('application/x-www-form-urlencoded');
	$req -> content($arguments);
	
#-- GET THE RESPONSE : PARSE OUT THE REQUEST ID and THE ESTIMATED TIME
	my($response) = $ua -> request($req);
	
	if ($response -> content =~ /Server Error/i) { die "Server Error at NCBI!!Sorry try later\n"; }
	$response -> content =~ /^\s{4}RID = (.*)$/m;   my($rid) = $1;
	$response -> content =~ /^\s{4}RTOE = (.*)$/m;	my($wait)= $1;
	unless ($rid && $wait)             { die "parse error: $!" };
	for (my $j=0; $j<=$wait/2; $j++)   {    print STDERR ".";	sleep 2;   }
	
	my($verif)=0;
	
	while ()
	{  		
		for (my $j=0; $j<=5; $j++)  { print  STDERR ".";  sleep 1; }
		
		$req = new HTTP::Request GET =>
		    "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=$rid";   
		$response = $ua->request($req);	   
		if    ($response->content =~ /Status=WAITING/im) {  next; }	
		elsif ($response->content =~ /Status=FAILED/im)  { print STDERR "Search $rid failed\n"; $verif=1; last; }	    
		elsif ($response->content =~ /Status=UNKNOWN/im) { print STDERR "Search $rid expired\n"; $verif=1; last; }	    
		elsif ($response->content =~ /Status=READY/im) 
		{	       
		    if   ($response->content =~ /ThereAreHits=yes/im){last;}	       
		    else { print STDERR "No hits found.\n";$verif=1;last;  }
		}
		elsif ($response->content =~ /can\'t connect/im)
		{ 
		    print STDERR "\nCan't connect to www.ncbi.nlm.nih.gov:80...new attempt"; 
		    if ($nb <3) { ++$nb; next; } 
		    else { print STDERR "sorry, BLAST $names[$i] failed after 3 attempts!!\n"; $verif=1; last;}
		}
		else { print STDERR "unknown error\n"; $verif=1; last; }
	    } 
	
	if($verif==1){ ++$i; next; }
	
#-- GET RESULT
	
	while ()
	{
	    sleep 3;
	    $req = new HTTP::Request GET => "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=$format&FILTER=off&EXPECT=$Eval&ALIGNMENTS=$align&DESCRIPTIONS=$align&ALIGNMENT_VIEW=$aln_view&RID=$rid";
		$response = $ua -> request($req);
	    
	    if   ($response->content =~ /Altschul/i) {  print STDERR "Search Complete\n"; push(@list_pdb,$response -> content);last; }
	    else { next; }
	}
	print SOR1 (@list_pdb);
	++$i;
    }
    
    undef (@list_encoded);
    
    close SOR1;
    return (@list_pdb);
    
}
#-----------------------------------------------------------------------------------------------------------------------
sub LOCAL_BLAST
{
    my ($blast_dir,$database,$query_file,$Eval,$align,$method,$matrix,$filter,$process,$gigablast,$database_expresso,$blast_dir_expresso,$runblast)=@_;
    my $n=0;
    if ($method=~ /^profile$/i && $gigablast=~ /^no$/i)      
    { 
	open (COM,"$blast_dir -p blastp -d $database -i $query_file -m 6 -M $matrix -v $align -b $align -F $filter -e $Eval -a $process|") or die;
    }
    elsif ($method=~ /^geneid$/i && $gigablast=~ /^no$/i) 
    { 
      
	open (COM,"$blast_dir -p blastp -d $database -i $query_file  -v $align  -b $align -F $filter -M $matrix -e $Eval -a $process|") or die;
    }

    elsif ($method=~ /^geneid$|^pdbid$/i && ($gigablast=~ /^yes$/i))
    {

	unless ($database eq "nr" || $database eq "pdb" || $database eq "refseq_protein" || $database eq "pdbaa" ) { print STDERR "\nsorry invalid database for gigablast\n";exit;} ;
	if ($database eq 'pdb')            { $database='pdbaa';}
	if ($database eq 'refseq_protein') { $database='refprot';}
	if ($database eq '')               { print STDERR "provide a valid database!\n" ;exit;}
	open (COM,"$runblast -d $database -p blastp  -e $Eval -v $align -F F \<$query_file |");
    }

    elsif ($method=~ /^profile$/i && ($gigablast=~ /^yes$/i)) { print STDERR "\nSorry method profile  can't be used with -gigablast option\n";exit;} 
    elsif ($method=~ /^pdbid$/i && $database=~ /expressopdb/)
    {	
	#BLAST pour Expresso  	
	open (COM,"$BLASTMAT; $BLASTDB;$blast_dir_expresso -p blastp -d $database_expresso -i $query_file -F $filter -e $Eval -M $matrix -v $align -b $align |") or die;       
    }

    else
    {  
	open (COM,"$blast_dir -p blastp -d $database -i $query_file -v 1 -b 1  -F $filter -e $Eval -M $matrix   -v $align -b $align -a $process |") or die;
    }
    
    unless ($quiet=~ /on/) { print STDERR "\nrun BLAST..."; } 
    
    my $name_database, my $posted, my $version;

    open (SOR2,">blast_result.txt") or die;
    
    $/="Query=";
    while (<COM>) 
    {
	if ($_=~ /Database: (\S+)/g)      { $name_database=$1;}
	if ($_=~ /Posted date: (.+?)\n/)  { $posted=$1;       }
       	if ($_=~ /BLASTP\s+(\S+)/o)        { $version=$1;}
	print SOR2 $_;
	push (@list_pdb,$_) ; 
	if ($_=~ /\s*(.+?)\s/) { print STDERR "\n$1 done";} 
    }
    close COM;
    close SOR2;
    print STDERR "\n";

    unless ($quiet=~ /on/i) { 
	                print STDOUT "
             Version: BLASTP $version
             Database: $name_database
             Posted date: $posted\n\n";
		           }
    shift (@list_pdb);  
    return (@list_pdb);
}

#-----------------------------------------------------------------------------------------------------------------------------
sub PARSING
{    
    my($list_pdb,$locale,$distant,$method,$quiet,$database,$gigablast)=@_;
    my(@list_pdb)=@$list_pdb; my(@result_not_sort)=();my $n=0;
    open (SOR,">webblast.log") or die;

    if ($gigablast=~ /^yes$/i) { $locale=2;$distant=0;}
    if ($gigablast=~ /^no$/i)  { $locale=1;}
    if ($distant==1)           { $locale=0;}
    
    foreach my $pdb_result(@list_pdb)
    { 	 
	my $query, my $length_query, my($pdb_id), my $comp=0;
       
	if ($pdb_result=~/No hits found/m) {  print SOR $pdb_result; next;}

	$pdb_result=~ s/ALIGNMENTS//;
	local $/=undef;
	my(@intra_res)= split(/(?=\n\n>)/s,$pdb_result);

       
	if ($distant==1) 
	{
	    my $version_d, my $database_d, my $poste_d;
	    undef $/; ($query,$length_query)=($intra_res[0] =~ /Query=\s+(\S+)\s+Length=\s*(\d+)/smo);  
	    $/="\n";
	    open (F3,"web_tempo.result") or die ;
	    while ($_=<F3>) 
	    {
		if ($_=~ /BLASTP\s+(\S+)/o)        { $version_d=$1;}
		if ($_=~ /Database:\s+(.+?)$/o)    { $database_d=$1;}
		if ($_=~ /Posted date:\s*(.+?)$/o) { $poste_d=$1; last;}
	    } 
	    close F3;
	    unless ($quiet=~ /on/i || $n>0) {++$n; 
		         print STDOUT "
             Version: BLASTP $version_d
             Database: $database_d
             Posted date: $poste_d\n\n";
		                       }

	} 
	else { ($query,$length_query)=($intra_res[0] =~ /\s*(.+?)\s.+?\(([\d,]+) letters/smo);}
	
	shift(@intra_res);
  	
	foreach my $intra_res(@intra_res) #look for the different results of the query
	{	    
	    my($aln_length,$identity)  = ($intra_res=~ /^\sIdentities = \d+\/(\d+)\s\((.+?)\)/im);
	    my($recouvrement,$gap)     = &RECOVER($intra_res,$aln_length,$length_query);	
	    my($evalue)                = ($intra_res=~ /Expect = (.+?)\s/im);
	    my($bits  )                = ($intra_res=~ /Score =\s+([\d.]+)\s/im);
	    
	    
	    unless ($method !~ /^geneid$/i) { if ($comp<=$bits) { $comp=$bits;} else { last;} }
	    
	    if ($query  eq "" || $length_query eq "" || $aln_length eq "" || $identity eq "" || $recouvrement eq "" || $gap eq "") 
	    { print SOR " can't parse $pdb_result"; next; }
	       
	    if ($method =~ /^pdbid$/i)
	    {
		if ($locale == 1)  { ($pdb_id) = ($intra_res=~ /^>(.{6})/im); $pdb_id=~ s/_//; $pdb_id=uc($pdb_id);}
		else               { ($pdb_id) = ($intra_res=~ /^>pdb\|(.{6})/im); $pdb_id=~ s/\|//; }	
		($evalue)                      = ($intra_res=~ /Expect = (.+?)\s/im);
	    
		push (@result_not_sort,("$query\t$pdb_id\t$evalue\t$identity\t$recouvrement\t"));
	    }	     
	    elsif ($method =~/^geneid$/i)
	    {	    
		if ($database !~ /pdb/i && $database !~ /swiss/i && ($locale=~/1|2/))
		{
	   
		     while  ($intra_res=~ />.*?(gb|prf|emb|sp|pir|tpe|ref|prf|dbj|ddbj|pdb)[\|]+([A-Za-z0-9_\.]+?)(\s|\|(.{1}))/sg) 
		     { 
			 my $databank =$1;
			 my $last     =$4;
			 my $refseq   =$2;
			 if ($databank eq "pdb") { $refseq.=$last } 
			 $refseq=~ s/\.\d+$//;
			 push (@result_not_sort,"$query\t$refseq\t$identity\t$recouvrement\t$bits\t$evalue\t$databank");
		     }
		   	  
		}  
		elsif ($database=~ /pdb|pdbaa/i && ($locale==1 || $locale==2))
		{	  
		   
		    my $refseq;
		    if($locale==1) {($refseq)  = ($intra_res=~ />(.*?)\s/o);      $refseq=~ s/_//; }
		    else           {($refseq)  = ($intra_res=~ /^>pdb\|(.{6})/im);$refseq=~ s/\|//;}
		    
		    unless ($refseq)  { print SOR $intra_res; next; }		    
		    push (@result_not_sort,("$query\t$refseq\t$identity\t$recouvrement\t$bits\t$evalue\tpdb"));  
		}  	 
		elsif ($distant==1 )
		{	  		    
		    my($resul)=&MULTI_EQUIVALENT($query,$identity,$recouvrement,$bits,$evalue,$intra_res);
		    push (@result_not_sort,"$resul");				  
		}
		elsif ($database=~ /swiss/i)
		{	  
		    my($refseq)       = ($intra_res=~ />.*?sp\|(.+?)\|/o);
		    unless ($refseq)  { print SOR $pdb_result; next; }	
		    $refseq=~ s/\.\d+$//;
		    push (@result_not_sort,("$query\t$refseq\t$identity\t$recouvrement\t$bits\tswiss_prot"));	
		}  
	    }	
		else {die;}		    	
	}
    }	
    close SOR;
    undef (@list_pdb);
    
    if ($method =~/^geneid$/i) { return (@result_not_sort); }
    else                                                                                                                         
    {
	my(@result_sort)= 
	    map {$_->[1]} 
	sort { $b->[0]<=>$a->[0]} 
	map {[/\t([\d.]+)%/,$_]} 
	@result_not_sort;
	
	undef(@result_not_sort);	
	return (@result_sort);	    
    }    
}

#-------------------------------------------------------------------------------------------------------------------

sub MULTI_EQUIVALENT
{
    my($query,$identity,$recouvrement,$bits,$evalue,$intra_res)=@_;

    my @result=();
    while  ($intra_res=~ />.*?(gb|prf|emb|sp|pir|tpe|ref|prf|dbj|ddbj|pdb)[\|]+([A-Za-z0-9_\.]+?)(\s|\|(.{1}))/g) 
    { 
	my $databank =$1;
	my $last     =$4;
	my $refseq   =$2;
	
	if ($databank eq "pdb") { $refseq.=$last } 
	$refseq=~ s/\.\d+$//;
	push (@result,"$query\t$refseq\t$identity\t$recouvrement\t$bits\t$evalue\t$databank");
    }
    return (@result);
}
#--------------------------------------------------------------------------------------------------------------------
sub AFFICHAGE_REFSEQ_PARSING
{
    my($result_sort,$cover_tresh,$identity_treshold,$out_file)=@_;
    my(@result_sort)=@$result_sort, my(@name_gid)=();my@resultats=();my $afficher="";
   
(my($entete)= sprintf("%-40s %-25s %-10s %-12s %-10s %-10s %-10s","Sequence Name","Accession number","Databank","%Identity","%Cover","BITS","Evalue")); 
    
    foreach my $result_sort(@result_sort)
    {     
	my($seq_name,$refseq_name,$identiq,$cover,$bits,$evalue,$bank)= split("\t",$result_sort);  	    
	($identiq)= split(/%/,$identiq);
	
	if ($identiq >= $identity_treshold && $cover >= $cover_tresh)
	{
	    push (@name_gid,">$seq_name\@$bank\_\_$refseq_name\n");
	    (($afficher).=  sprintf("%-40s %-25s %-10s %-12s %-10s %-10s %-10s ",$seq_name,$refseq_name,$bank,$identiq,$cover,$bits,$evalue));
	    $afficher.="\n";
	} 
	else {next;}	
    }

if ($afficher) { print "\n$entete\n\n"; print $afficher; }


if (@name_gid) { print STDOUT "\n**********************************************************************\n\n"; }
if ($out_file) { open (SOR,">$out_file") or die "can not open $out_file"; print SOR @name_gid; }
print STDOUT "\n", @name_gid;
close SOR;
}
#-------------------------------------------------------------------------------------------------------------

sub AFFICHAGE_PDB_PARSING 
{

    my($result_sort,$cover_tresh,$identity_treshold,$out_file)=@_;
    my(@result_sort)=@$result_sort, my @sortie=();
    
    print STDOUT "\n\n",(my($en_tete)= sprintf("%-40s %-10s %-10s %-12s %-10s","Sequence Name","PDB_id","Evalue","Identity(%)","Cover(%)")),"\n\n"; 
    
    foreach my $result_sort(@result_sort)
    {     
	my($seq_name,$pdb_name,$EValue,$identiq,$cover)= split("\t",$result_sort);  	    
	($identiq)= split(/%/,$identiq);
	
	if ($identiq >= $identity_treshold && $cover >= $cover_tresh)
	{
	    push (@pdb_list,$pdb_name);
	    print STDOUT ((my $afficher)= sprintf("%-40s %-10s %-10s %-12s %-10s",$seq_name,$pdb_name,$EValue,$identiq,$cover)),"\n";
	    push (@sortie,">$seq_name _P_ $pdb_name\n");
	} 
	else {next;}		
    }
    undef(@result_sort);
    print STDOUT "\n**********************************************************************\n\n";

#-- OUTFILE /STDOUT
    if   ($out_file) { open (SOR,">$out_file") or die "can not open $out_file"; print SOR @sortie; }
    print STDOUT @sortie;
    close SOR;

}
#-----------------------------------------------------------------------------------------------------------------------------------
sub PROFILE
{
    my($list_pdb,$out_file,$distant)=@_;
    my(@list_pdb)=@$list_pdb,  my(@sortie)=();
    my %names=();   my $i=0;    my($name)="";

    open (SOR1,">$out_file") or die;
    foreach my $pdb_result(@list_pdb)
    {     	
	if ($pdb_result =~ /No hits found/i) { next; }
	else
	{
	    ++$i;
	    if ($distant==1) {($name)   =($pdb_result =~ /Query=\s*(.+?)Length/smoi)  or die "\nparse error in distant profile\n";}
	    else             {($name)   =($pdb_result =~ /\s*(.+?)\(.*?letters/ismo)  or die "\nparse error in profile\n";}
	  
	    my($name1)= ($name=~ /(.+?)\s+$/);
 
	    open (SOR,">tempo_file_profile") or die "can not open tempo_file_profile";
	    print SOR "Query= $pdb_result";
	    close SOR;	    
	    
	    open(COM,"|t_coffee -other_pg seq_reformat -input blast_aln -in tempo_file_profile -output fasta_aln -out ${i}.profile");	     
	    close COM;
	    push (@sortie,">$name1 _R_ ${i}.profile\n");
     
	} 	
    }    
    unlink("tempo_file_profile");
    undef(@list_pdb); 
  
    print STDERR "\n**********************************************************************\n\n";
#-- OUTFILE /STDOUT
    if   ($out_file) { open (SOR1,">$out_file") or die "can not open $out_file"; print SOR1 @sortie; }
    print STDOUT @sortie;
    close SOR1;
         
}

#--------------------------------------------------------------------------------------------------------------------------------
sub ORGN
{   
    my($organism)=@_;
    $organism=~ s/_/ /;

    my(%orgs)= (
		
		'Homo sapiens'         =>'1',
		'Bos taurus'             =>'1',
		'Gallus gallus'         =>'1',
		'Viruses'              =>'1',
		'Bacteria'            =>'1',           
		'Eukaryota'            =>'1',
		'Mammalia'              =>'1',
		'Vertebrata'              =>'1',
		'All organisms'          =>'1',
		'Fungi'                 =>'1',
		'Primates'             =>'1',
		'Archaea'               =>'1',
                'Arabidopsis thaliana'   =>'1',
                'Caenorhabditis elegans'  =>'1',
                'Escherichia coli'        =>'1',
		'Mus musculus'             =>'1',
                'Drosophila melanogaster'   =>'1',
		);

    if (exists $orgs{$organism}) 
    { $organism=~ s/ /+/g; return ($organism); }
    else { print STDERR "organism not valid or syntax error, replace space by \"_\" \n"; &HELP(); }
   
} 
#------------------------------------------------------------------------------------------------------------------------------------

sub LIST_ORGA
{
    my(%orgs)= (
		
		'Homo sapiens'         =>'1',
		'Bos taurus'             =>'1',
		'Gallus gallus'         =>'1',
		'Viruses'              =>'1',
		'Bacteria'            =>'1',           
		'Eukaryota'            =>'1',
		'Mammalia'              =>'1',
		'Vertebrata'              =>'1',
		'All organisms'          =>'1',
		'Fungi'                 =>'1',
		'Primates'             =>'1',
		'Archaea'               =>'1',
                'Arabidopsis thaliana'   =>'1',
                'Caenorhabditis elegans'  =>'1',
                'Escherichia coli'        =>'1',
		'Mus musculus'             =>'1',
                'Drosophila melanogaster'   =>'1',
		);
    
    my (@cle)=keys(%orgs);
    foreach my $cle(@cle){ $cle=~ s/ /_/; }
    return (@cle);
}

