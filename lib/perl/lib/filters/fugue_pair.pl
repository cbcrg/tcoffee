#!/usr/bin/env perl
#

use Cwd;

if (!$ENV{FUGUE_LIB_LIST}){$ENV{FUGUE_LIB_LIST}="DUMMY";}
if (!$ENV{HOMSTRAD_PATH})  {$ENV{HOMSTRAD_PATH}="DUMMY";}
if (!$ENV{HOMS_PATH}){$ENV{HOMS_PATH}="DUMMY";}
$HOME=$ENV{HOME};

if ($ENV{TMP_4_TCOFFEE}){$tmp=$ENV{TMP_4_TCOFFEE};}
else
  {
    $tmp="$HOME/.t_coffee/tmp/";
   }


$out="stdout";
foreach $v(@ARGV){$cl.="$v ";}

if ( $cl=~/\-pdbfile1=(\S+)/){$pdb= $1;}
if ( $cl=~/\-infile=(\S+)/){$pep= $1;}
if ( $cl=~/\-outfile=(\S+)/){$out= $1;}

if ( ! -e $pdb)
  {
    print "COULD NOT READ PDB FILE |$pdb| [FATAL:fugue_pair]\n";
    exit (EXIT_FAILURE);
  }
elsif ( !-e $pep)
  {
    print "COULD NOT READ SEQ FILE |$pep| [FATAL:fugue_pair]\n";
    exit (EXIT_FAILURE);
  }


#Prepare the temporary directory
$ini_dir=cwd();
srand; 
$rand=rand(1000000);
$tmp_dir="$tmp/fugue_dir_$$_R_$rand";

`mkdir $tmp_dir`;
`cp $pdb $tmp_dir/struc.pdb`;
`cp $pep $tmp_dir/seq.pep`;
chdir $tmp_dir;
#Compute the profile/pdb ALignment

`joy struc.pdb >x 2>x`;
if ( !-e "struc.tem")
  {
    die "joy $pdb failed [FATAL:fugue_pair]\n";
    exit (EXIT_FAILURE);
  }
`melody -t struc.tem >x 2>x`;
if ( !-e "struc.fug")
  {
    die "melody failed [FATAL:fugue_pair]\n";
    exit (EXIT_FAILURE);
  }

`fugueali -seq seq.pep -prf struc.fug -print > struc.ali`;
if ( !-e "struc.ali")
  {
    die "fugueali failed [FATAL:fugue_pair]\n";
    exit (EXIT_FAILURE);
  }

`t_coffee -other_pg seq_reformat -in struc.ali -output fasta_aln >struc.fasta_aln`;
if ( !-e "struc.fasta_aln")
  {
    die "seq_reformat failed [FATAL:fugue_pair]\n";
    exit (EXIT_FAILURE);
  }
chdir $ini_dir;


if ($out eq "stdout")
  {
    open ( F,"$tmp_dir/struc.fasta_aln");
    while (<F>){print $_;}
    close (F);
  }
else
  {
    `cp $tmp_dir/struc.fasta_aln $out`;
  }
`rm $tmp_dir/*`;
`rmdir $tmp_dir`;


