#!/usr/bin/env perl
#

use Cwd;

$HOME=$ENV{HOME};
if ($ENV{TMP_4_TCOFFEE}){$tmp=$ENV{TMP_4_TCOFFEE};}
else
  {
    $tmp="$HOME/.t_coffee/tmp/";
   }


$out="stdout";
foreach $v(@ARGV){$cl.="$v ";}

if ( $cl=~/\-profile1=(\S+)/){$prf1= $1;}
if ( $cl=~/\-profile2=(\S+)/){$prf2= $1;}
if ( $cl=~/\-outfile=(\S+)/){$out= $1;}

if ( ! -e $prf1)
  {
    print "COULD NOT READ PRF1 FILE |$prf1| [FATAL:profile_pair]\n";
    exit (EXIT_FAILURE);
  }
elsif ( !-e $prf2)
  {
    print "COULD NOT READ PRF2 FILE |$prf2| [FATAL:profile_pair]\n";
    exit (EXIT_FAILURE);
  }


#Prepare the temporary directory
$ini_dir=cwd();
srand; 
$rand=rand(1000000);
$tmp_dir="$tmp/profile_pair_dir_$$_P_$rand";

`mkdir $tmp_dir`;
`cp $prf1 $tmp_dir/prf1.aln`;
`cp $prf2 $tmp_dir/prf2.aln`;
chdir $tmp_dir;
#Compute the profile/pdb ALignment
`clustalw -profile1=prf1.aln -profile2=prf2.aln -outfile=aligned_prf.aln`;

if ( !-e "aligned_prf.aln")
  {
    die "profile_pair failed [FATAL:profile_pair]\n";
    exit (EXIT_FAILURE);
  }
chdir $ini_dir;


if ($out eq "stdout")
  {
    open ( F,"$tmp_dir/aligned_prf.aln");
    while (<F>){print $_;}
    close (F);
  }
else
  {
    `cp $tmp_dir/aligned_prf.aln $out`;
  }
`rm $tmp_dir/*`;
`rmdir $tmp_dir`;


