#!/usr/bin/env perl
use Statistics::Test::WilcoxonRankSum;
use Statistics::ChisqIndep;

foreach my $mode ("avg", "min", "max","geo")
{
    analyze_combo($ARGV[0],$mode);
}

sub analyze_combo
    {
	my $file=shift;
	my $mode=shift;
		
	@combolist=("I__","_E_","__L", "IE_", "I_L", "_EL","IEL");
	%H=();
	%B=();
	%NIN=();
	open (F, $file);

	while (<F>)
	{
	    my $line=$_;
	    
	    my @words = split(/\s+/, $line);

	    #print ("$words[4]\n");

	    if ($words[0] eq "AUC" && $words[4] eq $mode)
	    {
		$family=$words[2];
		$ref    =$words[6];
		$ppmode =$words[8];
		$hfamily{$family}=1;
		
		
		for ( $a=11; $a<94; $a+=4)
		{
		    $H{$family}{$ref}{$ppmode}{$words[$a+1]}{auc}=$words[$a+2];
		    $H{$family}{$ref}{$ppmode}{$words[$a+1]}{isone}=$words[$a+3];
		    $H{$family}{$ref}{$ppmode}{$words[$a+1]}{npp}=$words[10];
		    
		    #print ("$words[$a+1],$words[$a+2] $words[$a+3]\n");
		}
		
	    }
	    elsif ($words[0] eq "PPLIST")
	    {
		$family =$words[2];
		$nseq   =$words[4];
		$ref    =$words[6];
		$ppmode =$words[8];
		$depth  =$words[12];
		$rdepth =$words[14];
		
		$B{$ref}{$ppmode}{npp}++;
		$B{$ref}{$ppmode}{depth}+=$depth;
		$B{$ref}{$ppmode}{rdepth}+=$rdepth;
		$NIN{$family}=$nseq-3;
	    }
	}

	close (F);

	%finalauc=();
	foreach $ref ("I__", "_E_", "__L")
	{
	    
	    foreach $combo1 (@combolist)
	    {
		
		foreach $bst ("000", "080", "100")
		{
		    $ppmode=$combo1.$bst;
		    
		    foreach $ncolumns ("025", "100", "200")
		    {
			my %auc=h2auc($ref,$ppmode,$ncolumns);
			foreach $combo2 (@combolist)
			{
			    my $comboN=$combo2.$ncolumns;
			    #printf "MEASURE: $ref $ppmode $comboN %.2f %.2f outof %d trees and %d branches using $mode for the bs combination\n", $auc{$comboN}{auc}, $auc{$comboN}{isone}, $auc{$comboN}{n}, $auc{$comboN}{npp};
			    $finalauc{$ref.$ppmode.$comboN}{auc}  =$auc{$comboN}{auc};
			    $finalauc{$ref.$ppmode.$comboN}{isone}=$auc{$comboN}{isone};
			    $finalauc{$ref.$ppmode.$comboN}{n}    =$auc{$comboN}{n};
			    $finalauc{$ref.$ppmode.$comboN}{npp}  =$auc{$comboN}{npp};
			    $finalauc{$ref.$ppmode.$comboN}{nin}  =$auc{$comboN}{nin};
			    
			    
			}
		    }
		}
	    }
	}
    

	
	if ($dotats)
	{
	    h2p();
	}
	
	
	foreach my $ncol ("025", "100","200")
	{
	    foreach my $ref ("__L", "_E_", "I__")
	    {
		;
		print_table ($mode,$ref,"TABLEAUC:", "auc",$ncol);
		print_table ($mode,$ref,"TABLEISONE:", "isone", $ncol);
	    }
	}
    }
sub print_table
{
    my $mode=shift;
    my $test_tree=shift;
    my $name=shift;
    my $content=shift;
    my $ncol=shift;
    my %colH;

    printf ("TABLE: ##################################################################################\n");
    printf ("TABLE: #                                                                                #\n");
    printf ("TABLE: #                                                                                #\n");
    printf ("TABLE: #         ref_tree: $test_tree                                                   #\n");
    printf ("TABLE: #         ncolumns: $ncol                                                        #\n");
    printf ("TABLE: #         content : $content                                                     #\n");
    printf ("TABLE: #         mode    : $mode                                                        #\n");
    printf ("TABLE: #                                                                                #\n");
    printf ("TABLE: #                                                                                #\n");
    printf ("TABLE: ##################################################################################\n");
   
	

    printf "$name\t\t";
    foreach $c2 (@combolist){printf "$c2\t";}
    printf ("npp\t%%br.\t#trees\trdepth\tdepth\n");
    foreach $c1 (@combolist)
    {
	printf "$name\t$c1\t";
	foreach $c2 (@combolist)
	{
	    $measure=$test_tree.$c1."080".$c2.$ncol;
	    printf ("%.3f\t", $finalauc{$measure}{$content});
	    push (@{$colH{$c2}},$finalauc{$measure}{$content});
	}
	printf("%d\t"   , $finalauc{$measure}{npp});
	
	printf("%6.2f\t" ,($finalauc{$measure}{nin}==0)?0:100*($finalauc{$measure}{npp}/$finalauc{$measure}{nin}));
	printf("%d\t"  , $finalauc{$measure}{n});
	printf ("%.2f\t"  , ($B{$test_tree}{$c1."080"}{npp}==0)?0:$B{$test_tree}{$c1."080"}{rdepth}/$B{$test_tree}{$c1."080"}{npp});
	printf ("%.2f\t"  , ($B{$test_tree}{$c1."080"}{npp}==0)?0:$B{$test_tree}{$c1."080"}{depth }/$B{$test_tree}{$c1."080"}{npp});
	printf ("\n");
    }
    my %done;
    foreach $c1(@combolist)
    {
	foreach $c2 (@combolist)
	{
	    my @list1;
	    my @list2;
	    
	    for ($a=0; $a<$#combolist; $a++)
	    {
		push (@list1,  $colH{$c1}->[$a]);
		push (@list2,  $colH{$c2}->[$a]);
	    }
	    if ($c2 eq "I_L"){;}
	    elsif (!$done{$c1}{$c2})
	    {
		my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();
		$wilcox_test->load_data(\@list1, \@list2);
		printf ("WILCOXTABLE:\t$content\t$c1\tvs\t$c2\t%.5f\n", $wilcox_test->probability());
		$done{$c1}{$c2}=$done{$c2}{$c1}=1;
	    }
	}
    }

}

sub h2auc
{
    my $ref=shift;
    my $ppmode=shift;
    my $ncolumns=shift;
    
    my %auc;
    
    for my $f (keys (%hfamily))
    {
	foreach my $combo (@combolist)
	{
	    my $comboN=$combo.$ncolumns;
	    #print ("$f $ref $ppmode $comboN --- $H{$f}{$ref}{$ppmode}{$comboN}{auc}\n");
	    if ($H{$f}{$ref}{$ppmode}{$comboN})
	    {
	
		$auc{$comboN}{auc}  +=$H{$f}{$ref}{$ppmode}{$comboN}{auc};
		$auc{$comboN}{isone}+=$H{$f}{$ref}{$ppmode}{$comboN}{isone};
		$auc{$comboN}{npp}  +=$H{$f}{$ref}{$ppmode}{$comboN}{npp};
		$auc{$comboN}{nin}  +=$NIN{$f};
		$auc{$comboN}{n}+=1;
	    }
		
	}
    }
    foreach my $k (keys(%auc))
    {
	#print ("$k --  $auc{$k}{auc} -- $auc{$k}{ione} -- $auc{$k}{n}\n");
	if ($auc{$k}{n}>0)
	{
	    $auc{$k}{auc}/=$auc{$k}{n};
	    $auc{$k}{isone}/=$auc{$k}{n};
	}
    }
    return %auc;
}

  
sub h2p
{
    foreach my $ref ("I__", "_E_", "__L")
    {
	foreach my $combo (@combolist)
	{
	    foreach my $bst ("000", "080", "100")
	    {
		    my $ppmode=$combo.$bst;
		    for my $combo1 (@combolist)
		    {
			for my $combo2 (@combolist)
			{
			    foreach my $ncol ("025", "100", "200")
			    {
				my $wilcox_test1 = Statistics::Test::WilcoxonRankSum->new();
				my $wilcox_test2 = Statistics::Test::WilcoxonRankSum->new();

				my @aauc1;
				my @aauc2;
				my @aisone1;
				my @aisone2;
				my ($nauc, $nisone,$ntest);

			
				foreach my $f (keys(%hfamily))
				{
				    if (exists ($H{$f}) && exists ($H{$f}{$ref}) && exists ($H{$f}{$ref}{$ppmode}) && exists ($H{$f}{$ref}{$ppmode}{$combo1.$ncol}{auc}) && exists ($H{$f}{$ref}{$ppmode}{$combo2.$ncol}{auc}))
					{
					    my $auc1  =$H{$f}{$ref}{$ppmode}{$combo1.$ncol}{auc};
					    my $auc2  =$H{$f}{$ref}{$ppmode}{$combo2.$ncol}{auc};
					    my $isone1=$H{$f}{$ref}{$ppmode}{$combo1.$ncol}{isone};
					    my $isone2=$H{$f}{$ref}{$ppmode}{$combo2.$ncol}{isone};

					    
					    
					    if ($auc1>0 || $auc2>0)  {$nauc++;}
					    if ($isone1>0 || $isone2>0){$nisone++;}
					    $ntest+=2;
					    
					    push (@aauc1,$auc1);   
					    push (@aauc2,$auc2);
					    push (@aisone1, $isone1);
					    push (@aisone2, $isone2);
					}
				}

				
				
				if ($nauc>0 && $nauc<$ntest)
				{
				    my @sauc1= sort {$a <=> $b} @aauc1;
				    my @sauc2= sort {$a <=> $b} @aauc2;
				    
				    $wilcox_test1->load_data(\@sauc1, \@aauc2);
				    my $prob_auc = $wilcox_test1->probability();
				    my $auc1=$finalauc{$ref.$ppmode.$combo1.$ncol}{auc};
				    my $auc2=$finalauc{$ref.$ppmode.$combo2.$ncol}{auc};
				    my $n   =$finalauc{$ref.$ppmode.$combo1.$ncol}{n};
				    				    
				    printf "WILCOXAUC: $ref$ppmode$combo1$ncol $ref$ppmode$combo2$ncol %.2f %.2f %.8f on %d samples\n",$auc1, $auc2, $prob_auc, $n;
				}
				if ($nisone>0 && $nisone<$ntest)
				{
				    my @sisone1= sort {$a <=> $b} @aisone1;
				    my @sisone2= sort {$a <=> $b} @aisone2;
				    $wilcox_test2->load_data(\@sisone1, \@sisone2);				    
				    my $prob_isone = $wilcox_test2->probability();

				    
				    my $isone1=$finalauc{$ref.$ppmode.$combo1.$ncol}{isone};
				    my $isone2=$finalauc{$ref.$ppmode.$combo2.$ncol}{isone};
				    my $n     =$finalauc{$ref.$ppmode.$combo1.$ncol}{n};

				    my $p1=$isone1*$n;
				    my $n1=$n-$p1;

				    my $p2=$isone2*$n;
				    my $n2=$n-$p2;

				    my $obs = new Statistics::ChisqIndep;
				    $obs->load_data([[$p1,$n1],[$p2,$n2],]);
				    $prob_isone=$obs->{p_value};
				    
				    printf "WILCOXISONE: $ref$ppmode$combo1$ncol $ref$ppmode$combo2$ncol %.2f %.2f %.8f on %d samples\n",$isone1, $isone2, $prob_isone, $n;
				}
			    }
			}
		    }
	    }
	}
    }
    return;
}

