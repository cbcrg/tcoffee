#!/usr/bin/env perl
#
#
#date : 04/04/05
#prog : webblast.pl
#subj : make a BLAST/WU-BLAST (by HTTP request or locally) against a database with a file containing sequences in fasta format
#auth : Fabrice Armougom
####### method genid, pdbid and profile 
#
#############################################################################################

use Env qw(HOME);
use lib "/lib/perl5/site_perl";
use lib "/lib/perl5/5.8.1";
use LWP::UserAgent;                                          # @@@@@@@ #                 
use HTTP::Request::Common qw(POST);                         # @/^   ^\@ # 
use URI::Escape;                                           # @/ -   - \@ #
use Getopt::Long;                                         ##  \   ^   /  ##
use strict;                                              ##    |  0  |    ##
use warnings;                                           ####### \ _ / #######
#########################################################################################


my(@list_encoded)=(), my(@list_pdb)=(), my(%deja_vu)=(), my(@pdb_list)=(),
my($i)=0, my(@names)=(), my $locale, my $distant, my $database, my $blast_way;

my($ua)= LWP::UserAgent->new;

##-- Variables d'environnements 

my($database_var)= $ENV { 'DATABASE' };
my($blast_var)   = $ENV { 'BLAST_DIRECTORY' };

##-- Recupere Options/parametres du BLAST && controle des options ds OPTIONS_GET
 
my($program,$database_line,$blast_line,$query_file,
   $out_file,$identity_treshold,$cover_tresh,$Eval,
   $iter,$i_tresh,$align,$description,$matrix,$filter,$method,$orgn)= &OPTIONS_GET();

##-- Determination BLAST || WU-BLAST, LOCAL || DISTANT && Controle database/programme

if ((($database_line || $database_var) && ($blast_line || $blast_var)) && (my($db_web)=&NCBI_DATABASE ($database_line)) eq "1") 
{   
    ($database)=$database_line || $database_var;
    my($blast_tp)=$blast_var || $blast_line; $locale=1;
    $blast_way = &CONTROLE_DB_PG($database,$blast_tp,$program);
    print STDERR "\nRUN BLAST LOCALY\n";
}
else
{
    $database = &NCBI_DATABASE ($database_line); $distant=1;
    print STDERR "\nRUN BLAST AT THE NCBI\n"; 
}
 
##-- Fixation de certains parametres selon la methode geneid/pdbid/profile

if ($method=~  /pdbid/i && $distant)
{
    $database="pdb";  
    $align="1";  
    $description="1";
}

if ($method=~ /geneid/i && $distant)
{  
    $matrix = "BLOSUM80";
    unless ($database eq "nr" || $database eq "swissprot" || $database eq "pdb")  { $database="refseq_protein";} ; 
    $align="1";  
    $description="1";
}

if ($method=~ /geneid/i && $locale) 
{ 
    $matrix = "BLOSUM80"; 
    $align="1";  
    $description="1"; 
}  

if ($method=~ /profile/i && $distant) 
{  
     unless ($database eq "nr" || $database eq "swissprot")  { $database="pdb"; } ;
}
###

print STDERR "   
                 
             Program : $program
             Database : $database
             Method : $method
             Query_file : $query_file
             Out_file : $out_file
            
             Evalue treshold : $Eval
             Matrix : $matrix
             Filter : $filter
             Organism : $orgn
            
             Blast_identity_threshold : $identity_treshold
             Cover threshold : $cover_tresh

             Number of descriptions : $description
             Number of alignments :  $align
            
***************************************************************\n\n";


#-- LOCAL/DISTANT BLASTP
   
if    ($locale)   { @list_pdb= &LOCAL_BLAST ($blast_way,$database,$query_file,$Eval,$align,$description,$method,$matrix,$filter);   }
elsif ($distant)  { @list_pdb= &WEB_BLAST   ($query_file,$Eval,$program,$database,$matrix,$method,$align,$description,$orgn,$filter); }
else              { die " Report bug to armougom\@igs.cnrs-mrs.fr\n"}; 

#-- PARSE BLAST RESULTS -> MAKE A PDB_ID LIST
if ($method =~ /^pdbid$/i) 
{ 
    my(@result_sort)= &PARSING (\@list_pdb,$locale,$distant,$method);
                      &AFFICHAGE_PDB_PARSING (\@result_sort,$cover_tresh,$identity_treshold,$out_file);
                      exit;
}

#-- PARSE BLAST RESULT -> MAKE LIST OF REFSEQ ID
elsif ($method =~ /^geneid$/i)
{ 
    my(@result_sort)= &PARSING (\@list_pdb,$locale,$distant,$method);    
                      &AFFICHAGE_REFSEQ_PARSING (\@result_sort,$cover_tresh,$identity_treshold,$out_file);
                      exit;
}

#-- PARSE BLAST RESULT -> MAKE PROFILE
elsif ($method=~ /^profile$/i) { &PROFILE (\@list_pdb,$out_file); exit; }


else { die "FATAL ERROR :  Method or database error" ;}


exit;
                                               ##############
###############################################  FONCTIONS  ####################################################################
                                              ##############

sub CONTROLE_DB_PG
{
    my($database,$blast_dir,$program)=@_;
   
    unless ($database =~ /refseq_protein$/) #refseq format different
    {
if (! -e $database) { die "$database file  does not exist\n";}
elsif ( -d $database) { die "$database must be a file, not a directory\n"; }
else 
{    
    open (FIC1,"$database") or die "can not open $database\n";
    while ($_=<FIC1>)  { if ($_=~ /^>.+/) {last;}  else {die "$database is not in fasta format \n"; } }
    close FIC1;
} 
    }
   
    if ($blast_dir !~ /\/$/) { $blast_dir.="/"; }    
    my ($blastall) = $blast_dir . "blastall";
    if (! -e $blastall) { $blastall=$blast_dir."wu-blastall"; }

    if ($program=~ /blastp/i) { if (! -e $blastall) { die "$blastall program not found \n";} }
    
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
    else                          { return (1); }
} 
#------------------------------------------------------------------------------------------------------------------------
sub HELP
{   
    my($org,@orga)= &LIST_ORGA();
    my ($list_orga)=join(', ',@orga);

    print STDERR "
                      usage: web_blast.pl -infile <fasta file> options []\n
                            -program ...... Program Name (blastp)
                                    Default = blastp
                            -database ..... Database at NCBI (nr, pdb, swissprot,refseq_protein) or indicate a local fasta file
                                    Default = pdb at NCBI
                    -infile ....... Query_file = a list of sequences in fasta format
                    -outfile ...... Name the outfile for make a template file for t_coffee
                                    Default = STDOUT
                            -evalue ....... Evalue treshold Default = 1;
                            -matrix ....... PAM30 PAM70 BLOSUM45 BLOSUM80 
                                            Default BLOSUM62
                            -method ....... geneid,pdbid, profile
                                            Default pdbid
                            -filter ....... L or R or m or C or V (for more information see NCBI doc for filter)
                                            Default = Off
                            -organism ..... $list_orga  are available 
                                           Default is All_organisms
                    -identity ..... blast_identity_thresold = provide a % for view only the results upper or equal to the threshold
                                    Default = 0
                            -cover ........ Cover threshold = provide a % : sequence covering 
                                            Default = 0
                            -align ........ Number of alignments to display
                                            Default = 10
                            -description .. display (up to) this number of descriptions
                                            Default = 10
                            -blast_dir .... Indicates where your BLAST directory  is installed localy
                           

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
       "iter=i"               =>\$opt{iter},
       "i_tresh=f"             =>\$opt{i_tresh},
       "align=i"                =>\$opt{align},
       "description=i"           =>\$opt{description},
       "matrix=s"                 =>\$opt{matrix},
       "filter=s"                  =>\$opt{filter},
       "method=s"                   =>\$opt{method},
               "organism=s"                  =>\$opt{organism},
       );
  
    if ($ARGV[0]) {print "Unprocessed by Getopt::Long\n $ARGV[0]\n"; &HELP();} 
   
  

    my($evalue_tresh)= $opt{'evalue'};       unless ($evalue_tresh) { $evalue_tresh=1;};
    my($cover_tresh) = $opt{'cover'};         unless ($cover_tresh)  { $cover_tresh=0;};
    my($query_file)  = $opt{'infile'};         unless ($query_file)   { print STDERR "Flag -infile must be defined\n"; &HELP();};
    my($outfil)     = $opt{'outfile'};          unless ($outfil)       { $outfil="";};
    my($treshold)    = $opt{'treshold'};         unless ($treshold)     { $treshold=0;};
    my($blast_dir)   = $opt{'blast_dir'};         unless ($blast_dir)    { $blast_dir="";};
    my($database)    = $opt{'database'};           unless ($database)     {$database="";};
    my($program)     = $opt{'program'};             unless ($program)      { $program="blastp";};  
    my($iter)        = $opt{'iter'};                  unless ($iter)        {$iter=2;} ;
    my($i_tresh)     = $opt{'i_tresh'};                unless ($i_tresh)     {$i_tresh=0.002;}; 
    my($align)       = $opt{'align'};                   unless ($align)       {$align= 10;};
    my($description) = $opt{'description'};              unless ($description) {$description=10;};
    my($matrix)      = $opt{'matrix'};                    unless ($matrix)      {$matrix="BLOSUM62";};
    my($fiter)       = $opt{'filter'};                     unless ($filter)      {$filter="F";};
    my($method)      = $opt{'method'};                      unless ($method)      {$method="pdbid";};
    my($organism)    = $opt{'organism'};                     unless ($organism)    {$organism="All organisms"};

  
    unless ($method =~ /(^geneid$|^pdbid$|^profile$)/i)  { print STDERR "unknown method for the flag -method\n";&HELP(); }

    if ($treshold <0 || $treshold >100)               { print STDERR "\nout of range for the option -treshold \n"; &HELP();}  
    if ($cover_tresh <0 || $cover_tresh >100)         { print STDERR "\nout of range for the option -cover \n"; &HELP();} 
    if ($description<=0 || $align <=0)                { print STDERR "\n error with option description or align\n"; &HELP();}
   
    my($orgn,@all_orgn)= &ORGN($organism);
  
    return ($program,$database,$blast_dir,$query_file,
    $outfil,$treshold,$cover_tresh,$evalue_tresh,
    $iter,$i_tresh,$align,$description,$matrix,$filter,$method,$orgn);
}

#--------------------------------------------------------------------------------------------------
sub RECOVER
{
    my($pdb_result,$aln_length,$length_query)=@_;
    my $nb_gap=0;

    $length_query  =~ s/,//g;
    my ($requete)  =  join('',($pdb_result=~/^Query:(.*)\n/gm));
        $requete   =~ s/[^A-Z-]//g;
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
    my($query_file,$Eval,$program,$database,$matrix,$method,$align,$description,$orgn,$filter)=@_;
    my $aln_view, my $format="Txt";
    
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
print STDERR "blast $names[$i]...";

#-- BUILD THE REQUEST
my($arguments) = "CMD=Put&ENTREZ_QUERY=$orgn&CDD_SEARCH=off&FILTER=$filter&MATRIX_NAME=$matrix&PROGRAM=$program&DATABASE=$database&QUERY=" . $encoded_seq;

my($req) = new HTTP::Request POST => 'http://www.ncbi.nlm.nih.gov/blast/Blast.cgi';
$req -> content_type('application/x-www-form-urlencoded');
$req -> content($arguments);

#-- GET THE RESPONSE : PARSE OUT THE REQUEST ID and THE ESTIMATED TIME
my($response) = $ua -> request($req);
     
        if ($response -> content =~ /Server Error/i) { die "Server Error at NCBI!!Sorry but try later\n"; }
$response -> content =~ /^\s{4}RID = (.*)$/m;   my($rid) = $1;
$response -> content =~ /^\s{4}RTOE = (.*)$/m; my($wait)= $1;
unless ($rid && $wait) { die "parse error: $!" };
sleep $wait;

my($verif)=0;

while ()
{         
    sleep 3;
    
    $req = new HTTP::Request GET =>
            "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=$rid";   
    $response = $ua->request($req);    
   
    if    ($response->content =~ /Status=WAITING/im) { next; } 
    elsif ($response->content =~ /Status=FAILED/im)  { print STDERR "Search $rid failed\n"; $verif=1; last; }     
    elsif ($response->content =~ /Status=UNKNOWN/im) { print STDERR "Search $rid expired\n"; $verif=1; last; }     
    elsif ($response->content =~ /Status=READY/im) 
    {        
if ($response->content =~ /ThereAreHits=yes/im){ print STDERR "Search complete\n"; last;}        
else { print STDERR "No hits found.\n";$verif=1;last; }
    }
    else { die "unknown error\n"; }
} 

if($verif==1){ ++$i; next; }

#-- GET RESULT

while ()
{
    sleep 4;
    $req = new HTTP::Request GET => "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=$format&FILTER=off&EXPECT=$Eval&ALIGNMENTS=$align&DESCRIPTIONS=$description&ALIGNMENT_VIEW=$aln_view&RID=$rid";
    $response = $ua -> request($req);

    if   ($response->content =~ /Altschul/i) { push(@list_pdb,$response -> content);last; }
    else { next; }
}
++$i;
    }
 
    undef (@list_encoded);
    return (@list_pdb);
    
}
#-----------------------------------------------------------------------------------------------------------------------
sub LOCAL_BLAST
{
    my ($blast_dir,$database,$query_file,$Eval,$align,$description,$method,$matrix,$filter)=@_;
    if ($blast_dir =~ /wu-/) { $matrix= lc($matrix); }
    
    if ($method=~ /^profile$/i)      { `$blast_dir -p blastp -d $database -i $query_file -m 6 -M $matrix -v $align -b $description  -F $filter -e $Eval -o tempo_blast.result`;}
    elsif ($method=~ /^geneid$/i) { print STDERR "run BLAST..."; `$blast_dir -p blastp -d $database -i $query_file  -v $align -b $description  -F $filter -M $matrix -e $Eval -o tempo_blast.result`;}
    else
    {  print STDERR "run BLAST..."; `$blast_dir -p blastp -d $database -i $query_file -v 1 -b 1 -F $filter -e $Eval -M $matrix -o tempo_blast.result`;}
    
    $/="BLASTP";
    open (FIC1,"tempo_blast.result") or die "\ncan not open tempo_blast.result:$!";
    my @list_pdb=<FIC1>;
    close FIC1;
    unlink ("tempo_blast.result");
    shift (@list_pdb);
    print STDERR "\nBLAST end\n\n";
    return (@list_pdb);

}

#-----------------------------------------------------------------------------------------------------------------------------
sub PARSING
{
    
    my($list_pdb,$locale,$distant,$method)=@_;
    my(@list_pdb)=@$list_pdb; my(@result_not_sort)=();
    
    print STDERR "parsing\n";

    foreach my $pdb_result(@list_pdb)
    { 
my($pdb_id), my($evalue);

if ($pdb_result=~/No hits found/m) { next;}
       
my($query,$length_query)   = ($pdb_result =~ /Query= (.+?)\s.+?\(([\d,]+) letters/smo) or die;

if ($locale  && $method =~ /^pdbid$/i) 
{  ($pdb_id)               = ($pdb_result=~ /^>(.{6})/im) or die; $pdb_id=~ s/_//; $pdb_id=uc($pdb_id);} 

if ($distant && $method =~ /^pdbid$/i)
{  ($pdb_id)               = ($pdb_result=~ /^>pdb\|(.{6})/im) or die; $pdb_id=~ s/\|//; } 

if ($method =~ /^pdbid$/i)
{ ($evalue)                = ($pdb_result=~ /Expect = (.+?)\s/im) or die; }

my($aln_length,$identity)  = ($pdb_result=~ /^\sIdentities = \d+\/(\d+)\s\((.+?)\)/im) or die;
my($recouvrement,$gap)     = &RECOVER($pdb_result,$aln_length,$length_query) or die;

if ($method =~/^geneid$/i)
{   
    my($refseq)            = ($pdb_result=~ />.*?ref\|(.+?)\|/o) or die;
    $refseq=~ s/\.\d+$//;

    push (@result_not_sort,("$query\t$refseq\t$identity\t$recouvrement\t")); 
}  
if ($method =~ /^pdbid$/i)
{
    push (@result_not_sort,("$query\t$pdb_id\t$evalue\t$identity\t$recouvrement\t")); 
}

print STDERR "$query done\n";
    }
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

sub AFFICHAGE_REFSEQ_PARSING
{
    my($result_sort,$cover_tresh,$identity_treshold,$out_file)=@_;
    my(@result_sort)=@$result_sort, my(@name_gid)=();
   
    print STDERR "\n\n",(my($en_tete)= sprintf("%-40s %-25s %-15s %-15s ","Sequence Name","Accession number","Identity%","Cover%")),"\n\n"; 
    
    foreach my $result_sort(@result_sort)
    {     
my($seq_name,$refseq_name,$identiq,$cover)= split("\t",$result_sort);      
($identiq)= split(/%/,$identiq);

if ($identiq >= $identity_treshold && $cover >= $cover_tresh)
{
    push (@name_gid,">$seq_name\_\_$refseq_name\n");
  
    print STDERR ((my $afficher)= sprintf("%-40s %-25s %-15s %-15s",$seq_name,$refseq_name,$identiq,$cover)),"\n";
} 
else {next;} 
    }
    
    print STDERR "\n**********************************************************************\n\n";
    if   ($out_file) { open (SOR,">$out_file") or die "can not open $out_file"; print SOR @name_gid; }
    print STDOUT "\n", @name_gid;
    close SOR;
}
#-------------------------------------------------------------------------------------------------------------

sub AFFICHAGE_PDB_PARSING 
{
    my($result_sort,$cover_tresh,$identity_treshold,$out_file)=@_;
    my(@result_sort)=@$result_sort, my @sortie=();

    print STDERR "\n\n",(my($en_tete)= sprintf("%-40s %-10s %-10s %-12s %-10s","Sequence Name","PDB_id","Evalue","Identity(%)","Cover(%)")),"\n\n"; 
    
    foreach my $result_sort(@result_sort)
    {     
my($seq_name,$pdb_name,$EValue,$identiq,$cover)= split("\t",$result_sort);      
($identiq)= split(/%/,$identiq);

if ($identiq >= $identity_treshold && $cover >= $cover_tresh)
{
    push (@pdb_list,$pdb_name);
    print STDERR ((my $afficher)= sprintf("%-40s %-10s %-10s %-12s %-10s",$seq_name,$pdb_name,$EValue,$identiq,$cover)),"\n";
} 
else {next;}

foreach my $nom(@pdb_list)
{ 
    unless ($deja_vu{$nom}) 
    {
push (@sortie, ">$seq_name \_P\_ $nom\n");      
$deja_vu{$nom}=1;
    } 
}
    }
    undef(@result_sort);
    print STDERR "\n**********************************************************************\n\n";

#-- OUTFILE /STDOUT
    if   ($out_file) { open (SOR,">$out_file") or die "can not open $out_file"; print SOR @sortie; }
    print STDOUT @sortie;
    close SOR;

}
#-----------------------------------------------------------------------------------------------------------------------------------
sub PROFILE
{
    my($list_pdb,$out_file)=@_;
    my(@list_pdb)=@$list_pdb,  my(@profile_file)=();
    my($home)=$ENV{'HOME'};

    foreach my $pdb_result(@list_pdb)
    {     
if ($pdb_result =~ /No hits found/i) { next; }
else
{
    my $name_stdout, my $name_bis;
    
    my($name)=($pdb_result =~ /^Query= (.+?)\s/m) or die "\nparse error in profile\n"; 
    ($name_bis=$name)=~ s/\|/\\|/g;

    if (-e  "$home/.t_coffee/cache/$name.profile") 
    {     
push(@profile_file,">$name _R_ $home/.t_coffee/cache/$name.profile\n");
    }
    else
    {
open (SOR,">tempo_file_profile") or die "can not open tempo_file_profile";
print SOR "$pdb_result";
close SOR;
       `t_coffee -other_pg seq_reformat -input blast_aln -in tempo_file_profile -output clustalw_aln -out $home/.t_coffee/cache/$name_bis.profile`;
unlink("tempo_file_profile");

#($name_stdout=$name_bis)=~ s/\\//g;  
push (@profile_file,">$name _R_ $home/.t_coffee/cache/$name.profile\n");
    } 
   
}
      
    }
    undef(@list_pdb);
    if   ($out_file) { open (SOR1,">$out_file") or die "can not open $out_file"; print SOR @profile_file; }
    else { print STDOUT @profile_file;}
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
