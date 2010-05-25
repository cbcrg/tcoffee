#!/usr/bin/env perl
#
#Initialization
$program="thread_pair.pl";
use Cwd;
$HOME=$ENV{HOME};
if ($ENV{TMP_4_TCOFFEE}){$tmp=$ENV{TMP_4_TCOFFEE};}
else{$tmp="$HOME/.t_coffee/tmp/";}

foreach $v(@ARGV){$cl.="$v ";}
($infile, $pdbfile1, $thread_package, $out)=&my_get_opt ( $cl, "-infile=",1,1, "-pdbfile1=",1,1, "-thread_package=",0,2, "-outfile=",1,0);
if ( $out eq ""){$out="stdout"};

#Prepare the temporary directory
$ini_dir=cwd();
srand; 
$rand=rand(1000000);
$tmp_dir="$tmp/thread_pair_dir_$$_P_$rand";
`mkdir $tmp_dir`;
`cp $infile $tmp_dir/infile.pep`;
`cp $pdbfile1 $tmp_dir/struc.pdb`;
chdir $tmp_dir;

#Compute the profile/pdb ALignment

if ( $thread_package eq "")
  {
    `t_coffee -in Pstruc.pdb Sinfile.pep Mslow_pair -outfile result.aln -quiet`;
  }

if ( !-e "result.aln")
  {
    die "ERROR: $program failed [FATAL:$program]\n";
    exit (EXIT_FAILURE);
  }
chdir $ini_dir;


if ($out eq "stdout")
  {
    open ( F,"$tmp_dir/result.aln");
    while (<F>){print $_;}
    close (F);
  }
else
  {
    `cp $tmp_dir/result.aln $out`;
  }
`rm $tmp_dir/*`;
`rmdir $tmp_dir`;


sub my_get_opt
  {
    my @list=@_;
    my $cl, $a, $argv, @argl;
    
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
	
	$installed=`which $argv`;
	if ($optional==0){;}
	elsif ( $optional==1 && $argv eq "")
	  {
	    print STDERR "ERROR: Option $option must be set [FATAL:$program]\n";
	    exit (EXIT_FAILURE);
	  }
	if ($status==0){;}
	elsif ($status ==1 && $argv ne "" && !-e $argv)
	  {
	    print STDERR "ERROR: Option $argv must exist [FATAL:$program]\n";
	    exit (EXIT_FAILURE);
	  }
	elsif ( $status==2 && $argv ne "" && $installed=~/not found/)
	  {
	    print STDERR "ERROR: $argv is not installed [FATAL:$program]\n";
	    exit (EXIT_FAILURE);
	  }
      }
    return @argl;
    }
