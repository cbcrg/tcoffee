#! /usr/bin/perl
$version="1.00";

# Read   the command line
if ( $#ARGV==-1)
  {
    print "coglist2conc -i <list_file> -n <N-Boostrap> -o <outfile> -e MSA extension\n";
    print "list_file: 1 cog per line\n";
    print "-i: list of COG MSAs or List of COG names with the -e extension";
    print "-e: MSA extension";
    print "-o: default= coglist2conc.output\n";
    print "-n: default=100";
  }
foreach $arg (@ARGV){$command.="$arg ";}

#extract parameters
$threshold=100;
$trim_msa=0;

if(($command=~/\-e (\S+)/))
  {
    $extension=$1;
  }

if(($command=~/\-i (\S+)/))
  {
    open (F, "$1");
    while (<F>)
      {
	$m=$_;
	chomp($m);
	$m.=".$extension";
	@list=(@list,$m);
      }
    close (F);
  }
$n_rep=100;
if(($command=~/\-n (\d+)/))
  {
    $n_rep=$1;
  }

$output="coglist2conc.output";
if(($command=~/\-o (\S+)/))
  {
    open (O, ">$1");
  }
print O "coglist2conc.pl Version $version\n";
print O "CLINE: $command\n";

foreach $m (@list)
  {
    print "$m\n";
    if ( !-e $m){print "ERROR: $m Could not be found [FATAL]\n";}
  }


foreach $m (@list)
  {
    $r=&make_concatenated_aln ($m);
  }
close (O);
exit (EXIT_SUCCESS);

sub make_concatenated_aln
  {
    my $msa=@_[0];
    
    $pconc="conc$nconc.msa";
    $nconc++;
    $conc="conc$nconc.msa";
    
    if ($len==0)
      {
	`cp $msa $conc`;
      }
    else
      {
	`t_coffee -other_pg seq_reformat -in $pconc -in2 $msa -action +cat_aln -output fasta_aln > $conc`;
      }
    $len++;
    
    @conc_list=(@conc_list, $msa);
    
    $r=`msa2bootstrap.pl -i $conc -n $n_rep`;
    
    $r=~/AVERAGE BOOTSRAP:\ ([\d.]+)/;
    $avg=$1;

    $r=~/Bootstrapped Tree:\ (.+)/;
    $conc_tree=$1;
    
    $list=join (";", @conc_list);
    $r="CONC_ALN: $conc\nCONC_TREE: $conc_tree\nAVERAGE BOOTSTRAP: $nconc $avg\nLIST: $list\n";
    print STDOUT "$r";
    print O "$r";
  }
