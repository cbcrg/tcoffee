#! /usr/bin/perl
$version="5.09";
$pg="clean_tree";

# Read   the command line
if ( $#ARGV==-1)
  {
    print "clean_tree -trees <tree_file> -msa <list_of_msa_files> -add <redundancy file> -min_cogs <number> -min_species <number> -keep_cogs <list,> -keep_species <list,> -threshold <0-1> -trim_msa <0-1> -ttm <filename> -ttph <file name> -catorder <sim|tree> -bcog <avg|cons|xxx> -input <ptree|tree_file|msa_file>\n";
    print "\t-trim_msa <1-0>. 1= remove non common species, 0: keep all species, replace missing ones with a gap.\n";
    print "\t-ttm : filename for the distance matrix of the Tree of Trees\n";
    print "\t-ttph  : filename for newick upgma tree of trees\n";
    print "\t-catorder  : order for concatenating the COGs\n";
    print "\t-bcog: Starting point of the concatenation: avg, cons or any COG of the list\n";
    die;
  }
foreach $arg (@ARGV){$command.="$arg ";}
print "clean_tree Version $version\n";
print "CLINE: $command\n";
#extract parameters
$threshold=100;
$trim_msa=0;



$bcog="cons";
if(($command=~/\-bcog (\S+)/))
  {
    $bcog=$1;
  }


$catorder="tree";
if(($command=~/\-catorder (\S+)/))
  {
    $catorder=$1;
  }
$ttm_file="$pg.ttm";
if(($command=~/\-ttm (\S+)/))
  {
    $ttm_file=$1;
  }
$ttph_file="$pg.ttph";
if(($command=~/\-ttph (\S+)/))
  {
    $ttph_file=$1;
  }

if(($command=~/\-trim_msa/))
  {
    $trim_msa=1;
  }
if(($command=~/\-trees (\S+)/))
  {
    $tree_list=$1;
  }
$input="ptree";
if(($command=~/\-input (\S+)/))
  {
    $input=$1;
  }
if(($command=~/\-add (\S+)/))
  {
    $taxa_list=$1;
  }

if(($command=~/\-min_cogs (\d+)/))
  {
    $min_cogs=$1;
  }
if(($command=~/\-min_species (\d+)/))
  {
    $min_species=$1;
  }
if(($command=~/\-threshold ([\d.]+)/))
  {
    $threshold=$1;
  }

if ( $command=~/\-keep_cogs([^-]+)/)
  {
    @keep_cog_file=split (/[\s,]+/, $1);
    foreach $f (@keep_cog_file)
      {
	if ( !-e $f && $f){push (@keep_cog, $f);}
	else
	  {
	    open (F, "$f");
	    while (<F>)
	      {
		$line=$_;chomp $line;
		@array=split (/[\s,]+/, $line);
		
		push (@keep_cog, $array[2]);
	      }
	    close (F);
	  }
      }
  }

if ( $command=~/\-keep_species([^-]+)/)
  {
    @species_file=split (/[\s,]+/, $1);
    foreach $f (@species_file)
      {
	if ( !-e $f && $f){push (@keep_species, $f);}
	else
	  {
	    open (F, "$f");
	    while (<F>)
	      {
		$line=$_;chomp $line;
		@array=split (/,/, $line);
		
		push (@keep_species, $array[2]);
	      }
	    close (F);
	  }
      }
  }

if ( $input eq "ptree")
  {
    #Read the trees
    open (F, $tree_list);
    while (<F>)
      {
	$file.=$_;
      }
    close (F);
    
    $file=~s/\n//g;
    $file=~s/\r//g;
    
    
    
    
    @cog=($file=~/(COG\d+);/g);
    $file=~s/COG\d+;//g;
    @trees=split (/;/, $file);
    
    
    
    for ($a=0; $a<=@cog; $a++)
      {
	$hcog{$cog[$a]}{'ref_tree'}= $trees[$a];
      }
    
    
    #Clean Trees
    foreach $c (@cog)
      {
	$l=$hcog{$c}{'ref_tree'};
	$nt++;
	if ( !($l=~/NHX/)){;}
	else
	  {
	    $l=~s/GI/gi/g;
	    $l=~s/lcl\|/gi\|/g;
	    $l=~s/\[[^\]]+T\=(\d+)[^\]]+\]/\[$1]/g;
	    
	    
	    $l=~s/\[[^\]]+unknown[^\]]+\]//g;
	    #print "$l\n\n";
	    
	    
	    $l=~s/(gi[^:]+):([^\[]+)\[([^\]]+)\]/TI\@\1\@\3#T:\2/g;
	    @gi_taxid=($l=~/TI@([^#]+)/g);
	    
	    $l=~s/\@[^@]+\@//g;
	    $l=~s/\#//g;
	    
	    
	    $l=~s/\n//g;
	    
	    $l=~s/\[[^\]]+\]//g;
	    
	    
	    $l=~s/\)[^:,;T]+\:/\)\:/g;
	    $l=~s/\:[0123456789\-.]+E[0123456789\-.]+/:0.000/g;
	    $l=~s/\:(\d\.\d{4})\d+/:\1/g;
	    $hcog{$c}{'gi_taxid'}=join (":::", @gi_taxid);
	  }
	
	@seq=($l=~/TI[^T]+T/g);
	push (@clean_tree, $l);
	$hcog{$c}{'clean_tree'}=$l;
	
      }
  }
elsif ( $input eq "tree_file")
  {
    open (F, "$tree_list");
    while (<F>)
      {
	$l=$_;chomp($f);
	if ( $l=~/>/)
	  {
	    $l=~/>(\S+)/;
	    $f=$1;
	  }
	else
	  {
	    $f=$l;
	  }
	@cog=(@cog, $f);
	
	open (G, "$f");
	while (<G>)
	  {
	    $hcog{$f}{'clean_tree'}.=$_;
	  }
	close (G);
      }
    close (F);
  }
elsif ( $input eq "msa_file")
  {
    my $l;
    open (F, "$tree_list");
    while (<F>)
      {
	$l .=$_;
      }
    @msa_list=split (/\s+/, $msa_list);
    foreach $m (@msa_list)
      {
	$m=~/([^.]+)/;
	$c=$1;
	push (@cog, $c);
	$hcog{$c}{'ref_msa'}=$m;
	$ref_tree="$c.ref_tree";
	`msa2bootstrap -i $msa -n 1 -o $ref_tree`;
	open (F, "$ref_tree");
	while (<F>)
	  {
	    $hcog{$c}{'clean_tree'}.=$_;
	  }
	close (F);
      }
  }

#Check that the tree and the Alignments contain the same sequences
foreach $c (@cog)
  {
    my $g,@seq,$msa, $tree;
    $msa="$c.fasta";
    
    $g=$hcog{$c}{'gi_taxid'};
    if ( !-e "$c.fasta")
      {
	print "\nWARNING: Could Not Find $msa";
      }
    open (F, "$msa");
    while (<F>)
      {
	if ( />(\S+)/)
	  {
	    $n=$1;
	    if (!($g=~/$n/))
	      {
		$stop=1;
		print "Could Not Find Sequence $n in $msa\n";
	      }
	  }
      }
    close (F);
  }
if ( $stop)
  {
    die;
  }

#Expand Trees with omitted Taxa
open (F, "$taxa_list");
while (<F>)
  {

    @l=split (/\s+/, $_);
    print "@l\n";
    $c=$l[0];
    $t=$l[1];

    $sup1{$c}.="$t ";
    for ($d=2; $d<=@l;$d++)
      {
	$sup2{$c}{$t}.="$l[$d] ";
      }
  }
close (F);


foreach $c (@cog)
  {
    $hcog{$c}{'exp_tree'}=$hcog{$c}{'clean_tree'};
    if ($sup1{$c})
      {
	foreach $t1 (split (/ /, $sup1{$c}))
	  {
	    print "\nExpand: $c $t1 $sup2{$c}{$t1}";
	    foreach $t2 (split (/ /, $sup2{$c}{$t1}))
		     {
		       $hcog{$c}{'add_taxid'}.="$t1\@$t2\:\:\:";
		       $hcog{$c}{'exp_tree'}=&insert_leaf ($hcog{$c}{'exp_tree'},$t1,$t2);
		     }
	  }
      }
  }
#Filter cogs with Species to keep
@new_cog=();
print "@keep_species";
foreach $c (@cog)
  {
    $l=$hcog{$c}{'exp_tree'};
    
    $keep=1;
    foreach $s (@keep_species)
      {
	if ($s && !($l=~/$s/))
	  {
	    $keep=0;
	  }
      }
    if ( $keep)
      {@new_cog=(@new_cog, $c);}
    else 
      {
	print "Remove_COG: $c\n";
      }
  }

@cog=@new_cog;


#Trim the trees to the minimal subset
$ns=0;
#1-identify species in each cog
$n=0;
foreach $c (@cog)
  {

    $l=$hcog{$c}{'exp_tree'};
    @species=($l=~/TI[^T]+T/g);
    $hcog{$c}{'n_species'}=$#species+1;
    $hcog{$c}{'species'}=join ";", @species;
    foreach $s (@species)
      {
	$species_list{$c}{$s}++;
	$species_list{'full'}{$s}++;
	$species{$s}{$c}=1;

      }
    print "COGLIST $n $c Contains $hcog{$c}{'n_species'} SPECIES: $hcog{$c}{'species'}\n";
    $n++;
  }

#2-display species in each cog
print "SPECIES $species[0] ";
$n=0;
foreach $c (@cog)
  {
    printf "%d", $n/10;
    $n++;
  }
print "\nSPECIES $species[0] ";
$n=0;
foreach $c (@cog)
  {
    printf "%d", $n%10;
    $n++;
  }
print "\n";
#Filter species list with keep_cog list
print "\nKEEP_COG: @keep_cog\n";
if (@keep_cog)
  {
    @keep_species=();
    foreach $s (keys %{$species_list{'full'}})
      {
	$keep=1;
	foreach $c (@keep_cog)
	  {
	    if (!$species_list{$c}{$s})
	      {
		$keep=0;
	      }
	  }
	if (!$keep){$species_list{full}{$s}=0;}
	else
	  {
	    print "\nKEEP $s\n";
	    @keep_species=(@keep_species, $s);
	  }
      }
  $min_species=$#keep_species;
  }


#add species to keep (set them to the maximum
foreach $s (keys %{$species_list{'full'}})
  {
    foreach $s2 (@keep_species)
      {
	if ($s=~/$s2/)
	  {
	    $species_list{full}{$s}=$#cog+1;
	  }
      }
  }
    
foreach $s (sort {$species_list{'full'}{$b}<=>$species_list{'full'}{$a}} keys %{$species_list{'full'}})
  {
    printf "SPECIES %-10s",$s;
    foreach $c (@cog)
      {
	if ($species_list{$c}{$s}){print "*";}
	else {print "-";}
      }
    
    $ns++;
    $list.="$s ";
    $n=&count_trees ($list,%hcog);
    if ($min_cogs!=0 && $n>=$min_cogs)
      {
	$min_list=$list;
	$n_species=$ns;
	$keep_n_cogs=$n;
      }
    elsif ( $min_species!=0 && $ns<=$min_species)
      {
	$min_list=$list;
	$n_species=$ns;
	$keep_n_cogs=$n;

      }
    elsif (!$min_cogs && !$min_species)
      {
	$min_list=$list;
	$n_species=$ns;
	$$keep_n_cogs=$n;
      }
    print " Keep $n Trees with $ns Species\n";
  }
print "\n#KEEP $keep_n_cogs COGs with $n_species SPECIES: $min_list";


#Save All the trees
printf "#Save All the trees:\n";
printf "#\t*.clean_tree: original tree with geneames repalced by taxID\n";
printf "#\t*.exp_tree  : clean_tree with omitted species added\n";
printf "#\t*.trim_tree : exp_tree with non-common species removed (empty trees are not reported)\n";

foreach $c (@cog)
  {
    open (F, ">$c.clean_tree");
    print F "$hcog{$c}{'clean_tree'}";
    close (F);
    
    open (F,">$c.exp_tree");
    print  F "$hcog{$c}{'exp_tree'};";
    close (F);

    
    if ( ($hcog{$c}{'trim_tree'}=&trim_cog_tree ($min_list, $hcog{$c}{'exp_tree'})))
      {
	print "KeepCOG $c\n";
	push (@rcog, $c);
	$hcog{$c}{'file'}="$c.trim_tree";
	open (F,">$hcog{$c}{'file'}");
	print F "$hcog{$c}{'trim_tree'}";
	close (F);
      }
  }

#
#
# add the consensus_tree
#
#

$cname="COG000";
&treelist2ctree ($cname);
push (@rcog, $cname);


@cog=@rcog;
$nc=(@rcog)?$#cog+1:0;


#
#   Produce the distance Matrix between every pair of COG
#
print "#Kept a total of $nc Cog Trees with $n_species Common Species\n";

open (M, ">$ttm_file");


print "\nPMAT $nc";
print M "\n$nc";

for ($a=0; $a<$nc; $a++)
  {

    $c1=$cog[$a];
    printf "\nPMAT %-10s ",$c1;
    printf M "\n%-10s ",$c1;
    
    for ($b=0; $b<$nc; $b++)
      {

	$c2=$cog[$b];


	$t1=$hcog{$c1}{'file'};
	$t2=$hcog{$c2}{'file'};


	$command="t_coffee -other_pg seq_reformat -in $t1 -in2 $t2 -input newick_tree -input2 newick_tree -action +tree_cmp ";
	$r=`$command`;
	
	#print "$r";
	
	$r=~/T\: (\d+)/g;
	$v=$1;
	$v=(100-$v)/100;
	
	printf "%.4f ",$v;
	printf M "%.4f ",$v;
	
	$matrix{$c1}{$c2}=$v;
	$matrix {$c1}{tot}+=$v;
      }
  }
#
#
#   Compute the Tree of Trees
#
#

$ttph_file=&mat2upgma_tree ($ttm_file, $ttph_file);
printf "\nTTPH %s", &display_file ($ttph_file);
print  "\nOutput Tree of Trees: $ttph_file";

#
#
#   Identify The Best COG
#
#




foreach $c (@cog)
      {
	$matrix{$c}{avg}=$matrix{$c}{tot}/$nc;
      }

@scog=@cog;
if ($bcog eq "sim")
  {
    @scog=sort {$matrix{$a}{avg}<=>$matrix{$b}{avg}} @cog;
    $best_cog=$scog[0];
  }
elsif ( $bcog eq "cons")
  {
    $best_cog="COG000";
  }
else
  {
    $best_cog=$bcog;
  }

foreach $c (@scog)
      {
	$matrix{$c}{avg}=$matrix{$c}{tot}/$nc;
	printf "\n\tAVG %-9s %6.2f %6.2f",$c,$matrix{$c}{tot}, $matrix{$c}{avg};
	if ( $best_cog eq $c){ print "\t**** BEST COG ****"};
      }
printf "\n";
print "\n#Filter trees: Identify the tree closest to the centroid.\n";


#
#
# Order the COGs, starting from the best
#
#


if ($catorder eq "sim")
  {
    @scog=sort {$matrix{$best_cog}{$a}<=>$matrix{$best_cog}{$b}} @cog;
  }
elsif ($catorder eq "tree")
  {
    @scog=tree_sort ($ttph_file, $best_cog);
  }
else
  {
    print "ERROR: $catorder is an unknown mode for sorting COGs\n";
    exit (EXIT_FAILURE);
  }

foreach $c (@scog)
      {
	printf "\n\tBEST_COG:%-9s %-9s %6.2f ",$best_cog, $c,$matrix{$best_cog}{$c};
	if ( $matrix{$best_cog}{$c}<=$threshold) { print " -->keep";}
      }
printf "\n";
printf ("\n#Keep all the neighbor trees less than $threshold %% away from it\n");
@cog=();
foreach $c (@scog)
      {
	if ($matrix{$best_cog}{$c}<=$threshold)
	  {
	    printf "\nBEST_TREE:$hcog{$c}{'trim_tree'}";
	    if ($c ne "COG000")
	      {
		@cog=(@cog, $c);
	       print "\nCOG_MSA: $c.clean_msa"; 
	      }
	  }
      }
printf "\n";




#
#
#  Produce The Concatenated MSA
#
#
$t=$hcog{$cog[0]}{'trim_tree'};
if ($trim_msa)
  {
    @order=($t=~/TI([^T]+)T/g);
  }
else
  {
    foreach $t(keys (%{$species_list{'full'}}))
      {
	$t=~/TI([^T]+)T/g;
	push (@order,$1);
      }
  }


$cog_list=join ("__", @cog);
foreach $c (@cog)
  {
    open (G, ">$c.clean_msa");
    %aln=read_aln ($c,$hcog{$c}{'gi_taxid'}, $hcog{$c}{'add_taxid'});
    if (!%aln){next;}
    foreach $s(@order)
      {
	$seq=($aln{$s})?$aln{$s}:$aln{pad};
	$caln{$s}.="COG$seq";
	print G ">$s\n$seq\n";
      }
    close (G);
  }

foreach $s (@order)
  {
    print "\nCONC_ALN:>$s $cog_list";
    print "\nCONC_ALN:$caln{$s}\n";
  }
print "\n";
exit (EXIT_SUCCESS);

#Reads an alignment
#Replaces every sequence with its TaxID
#Replaces every 
sub read_aln
  {
    my ($aln, $gi_taxid, $add_taxid)=@_;
    my ($a,$pad,$g, $t, $t1, $t2, $name,%aln, @gi_taxid, @add_taxid, $gi, $taxid);

   
    $aln.=".fasta";
    if ( !-e $aln)
      {
	print "\nERROR: $aln Does Not Exist.";
	return %aln;
      }
    print "Read, rename and expand $aln\n";
    
    open (F, "$aln");
    while (<F>)
      {
	$s=$_;
	if ( $s=~/\>/)
	  {
	    $s=~s/\>GI/\>gi/;
	    $s=~/\>(\S+)/;
	    $name=$1;
	    $l=0;
	    
	  }
	else
	  {
	    chomp ($s);
	    $aln{$name}.=$s;
	    $l+=length ($s);
	  }
      }
    close (F);
    
    for ($a=0; $a<$l;$a++)
      {$aln{pad}.='-';}
    

    @gi_taxid=split (/:::/, $gi_taxid);
    foreach $g (@gi_taxid)
      {
	($gi, $taxid)=split (/@/, $g);
	if (!$aln{$gi})
	  {
	    $aln{$taxid}=$aln{pad};
	  }
	$aln{$taxid}=$aln{$gi};
	
      }

    @add_taxid=split (/:::/, $add_taxid);
    foreach $t (@add_taxid)
      {

	($t1, $t2)=split (/@/, $t);
	$aln {$t2}=$aln{$t1};
      }
   
    return %aln;
  }
sub insert_leaf
  {
    my $tree=@_[0];
    my $t1=@_[1];
    my $t2=      @_[2];
    my ($nt1,$nt2,$zero);

    
    $nt1="TI$t1\T";
    $nt2="TI$t2\T";
    $zero="0.0000";

    $v=($tree=~s/$nt1/\($nt1\:$zero,$nt2\:$zero\)/);


    return $tree;
  }
sub count_trees
  {
    my ($s, %h)=(@_);
    my (@list, $tree, $keep, $tot);
    
    @list=split (/\s/, $s);
    
    $tot=0;
    foreach $t (keys (%h))
      {

	$tree=$h{$t}{'clean_tree'};
	$keep=1;
	foreach $s (@list)
	  {
	    if (!($tree=~/$s/)){$keep=0;}
	  }
	$tot+=$keep;
      }
    return $tot;
    }



sub trim_cog_tree
  {
    my $species =shift (@_);
    my $tree=shift (@_);
    my ($s, @list, $new_tree);

    @list=split (/\s/, $species);

    foreach $s (@list)
      {
	if (!($tree=~/$s/)){return 0;}
      }
    
    open (F, ">tmp_seq");
    foreach $s (@list)
      {
	print F ">$s\n";
      }
    close (F);
    open (F, ">tmp_tree");
    print F "$tree;\n";
    close (F);

    $new_tree=`t_coffee -other_pg seq_reformat -in tmp_tree -input newick_tree -in2 tmp_seq -action  +tree_prune -output newick_tree`;
    return $new_tree;
  }
sub treelist2ctree 
  {
    my ($nc)=@_[0];
    my ($t,$f, $c, $ct, $command);
    
    $f=&vtmpnam();
    open (C, ">$f");
    foreach $c (keys(%hcog))
      {
	if ( $hcog {$c}{trim_tree})
	  {
	    print "$hcog{$c}{trim_tree}\n";
	    print C "$hcog{$c}{trim_tree}\n";
	  }
      }
    close (C);
    $ct="$nc.trim_tree";
    $command="msa2bootstrap.pl -i $f -input tree -o $nc.trim_tree";
    `$command`;
    print "$command";
    
    $hcog{$nc}{'file'}=$ct;
    $hcog{$nc}{'trim_tree'}=&display_file ($ct);
  }

sub msa2tree
  {
    my $msa=@_[0];
    my $tree=@_[1];
    if (!$tree){$tree=&vtmpnam();}

    `msa2bootstrap.pl -i $msa -n 1 -o $tree`;
    return $tree;
  }
sub mat2upgma_tree 
  {
    my $mat=@_[0];
    my $tree=@_[1];
    if (!$tree){$tree=&vtmpnam();}

   
    $command="msa2bootstrap.pl -i $mat -input matrix -n 1 -o $tree";
    `$command`;
    if ( !-e $tree)
      {
	print "ERROR: Could Not Run $command";
	exit (EXIT_FAILURE);
      }
    else
      {
	return $tree;
      }
  }

#
#
# Tree_sort
# 
#
#

sub tree_sort
  {
    my ($t,$c)=@_;
    my ($content, @sorted);
    
    $content=`t_coffee -other_pg seq_reformat -in $t -action +node_sort $c`;
    @sorted=($content=~/>(\S+)/g);
    return @sorted;
  }

#
#
# Temporary Files 
#
#


sub vtmpnam 
  {
    my $file;

    $ntmp++;
    $file="tmp4_$pg.$$.$ntmp";
    
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
sub display_file
  {
    my $file=@_[0];
    my $content; 

    open (F, "$file");
    while (<F>)
      {
	$content .=$_;
      }
    close (F);
    return $content;
  }
