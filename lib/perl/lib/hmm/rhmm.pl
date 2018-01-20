#!/usr/bin/env perl
# -filter field <field name> contains|equals|min|max value
# -bin
use HTTP::Date;
use Data::Dumper;
use strict;
use FileHandle;

srand();
my $TAG=0;
my $BIN=0;

my $index;
my $LOG_ZERO=-999999999;
my $LOG_UNDERFLOW=0.000000001;
my $DELTA=0.000001; #minimum difference for two values to be considered equal
my $HEADER;

if ($#ARGV ==-1)
  {
    
    print "****************** Description * ***********************\n";
    print "Trains Hidden Markov Models (Bawm and Welch)\n";
    print "Decode data using pre-computed models\n";
    print "Evaluate the decoding capacity using reference data\n";
    
    print "****************** Command Line  ***********************\n";
    print "rhmm.pl -test -model <model_file>  -constraints <model_file> -data <data_file> -train <trinmode> -nrounds <number of rounds> -niter <number of iterations> -outmodel <file> -outdata <file> -evaluate\n";
    print "****************** Flags      **************************\n";
    print "  -test...........................outputs two sample files\n";
    print "  ................................out.rhmm.model and out.rhmm.data that can be used as input\n";
    print "  ................................rhmm.pl -model out.rhmm.model -data rhmm.out.data -train bw -evaluate\n";
    
    print "  -model       <file>,<mode> .....File: input model from file.\n";
    print "  ................................Mode: 'topology' creates a model with -nst states\n";
    print "  ................................number of emission is set to the number of bin in data or to -nem\n";
    
    print "  -nem         <int>..............Number of emmissions when -model=topology\n";
    print "  -nst         <int>..............Number of states when -model=topology\n";
    print "  -constraints <file>,<mode>'.....File: input a model where every transition/emission become a seed.\n";
    print "  ................................Mode: 1 to forbid every undeclared transition/emission.\n";
    print "  -data  <file1 file2.. >,<mode>..File: input data from file(s).\n";
    print "  ................................Mode: 'model' to generate simulated data from the model -model.\n";
    print "  -action      <mode> ............'bw' to apply baum and welch training.\n";
    print "  ................................'decode' to apply the provided model onto the data\n";
    print "  ................................'model2data' generates -n artificial dataset -len long (with RST)\n";  print "  -single.........................Set model states and emission to single char values (max=62 values)";
    print "  ................................'data2model' estimates model from data (must have RST)\n";
    
    print "  -nrounds     <value>............Value; number of trainned models.\n";
    print "  -nit         <value>............Value: number of iterations in each round of trainning\n";
    print "  -outmodel    <file>.............File:  name of the output file containning the trainned model.\n";
    print "  -outmodel_f  <mode>.............Mode:  zeros or default.\n";
    print "  -output      <mode>.............Mode: data or fasta.\n";
    print "  -outdata     <file>.............File:  name of the output file containning the decoded data.\n";
    print "  -evaluate    <mode>.............Mode: 'posterior' evaluate posterior decoding predictions.\n";
    print "  ................................Mode: 'viterbi'   evaluate viterbi path predictions,\n";
    print "  ................................NOTE: RST field must be set in data.\n";
    print "  ................................NOTE: RST field declares the reference sate of each data item.\n";
    print "  -evaluation  <mode>.............Mode: 'sequence' evaluates the probability of a sequence of observations\n";
    print "****************** Data Format **************************\n";
    print "  Raw Data\n";
    print "  #d;<exp>;<record id>;bin;<emission bin>;\n";
    print "  Decoded Data\n";
    print "  #d;<exp>;<record id>;bin;<emission bin>;posterior;<posterior label>;viterbi;<viterbi label>\n";
    print "  Reference Data (used for the evaluation or the estimation of a model)\n";
    print "  #d;<exp>;<record id>;RST;<reference label>;\n";
    print "  Fasta\n";
    print "  ><seq1> <field name>\n";
    print "  <list of 1 letter symbols>\n";
    
    print "Note: extra fields are ignored\n";
    print "Note: output files can be used as input\n";
    print "Note: run rhmm.pl -test to get sample formats\n";
    print "****************** Model Format **************************\n";
    print "  #set;<state label>;<state label>;proba;\n";
    print "  #set;<state label>;<emission bin>;proba;\n";
    print "  #comment;<free text>;\n";
    print "\n";
    print "Note: to declare constraints one simply needs to declare transitions\n";
    print "******************Contact **************************\n";
    print "contact: cedric.notredame\@gmail.com\n\n";
    
    die;
  }

my $param;

$param=process_param (@ARGV);
if (!exists ($param->{c})){$param->{c}=1;}
my $nrounds=1;
my $nit=100;
my ($Data,$Model,$P);
my $newModel;

#process Parameters
if (!$param->{n}){$param->{n}=1;}
if (!$param->{len}){$param->{len}=3000;}
if ( $param-> {test})
  {
    $param->{model}="test.rhmm.model";
    my $Model=ODHC2model($Model);
    dump_model ($Model,$param->{model},"Reference_test");
    
    my $Data=model2data($Model, $param->{n}, $param->{len});
    $param->{data}="test.rhmm.data";
    dump_data ($Data,$param->{data}, "ReferenceData");
    printf "      Data : $param->{data}\n";
    printf "      Model: $param->{model}\n";
    if (!$param->{action}){exit(0);}
  }

#Read the data, the model and estimate the constraints
$Data =main_undump_data ($Data, $param);

#Read/Create the model
if (-e $param->{model}){$Model=undump_model($param->{model});}
elsif ($param->{model} eq "topology")
  {
    print STDERR "----------------------\n";  
    $Model=topology2model($Data,$param->{nst},$param->{nem});
  }
if (-e $param->{constraints}){$Model=undump_model($param->{constraints});}
elsif ($Model){$Model=model2constraints($Model, $param->{constraints})};

#Reformat the Model
if ($param->{single}){$Model=model2single_symbol_state($Model);}

if ($param->{action} eq "bw")
  {
    check_data ($Data, "data", $Model, "model"); 
    
    if (!$param->{nrounds}){$param->{nrounds}=1;}
    if (!$param->{nit}){$param->{nit}=100;}
    ($Data,$Model,$P)=&baum_welch($Data,$Model,$param);
    if ($param->{sort}){($Data,$Model)=data2recode_on_em($Data, $Model);}
    $Model=modelL2model($Model);
  }
elsif ($param->{action} eq "model2data")
  {
    check_data ($Model, "model"); 
    $Data=model2data($Model, $param->{n}, $param->{len});
    ($Data,$P)=decode ($Model,$Data,$param);
  }
elsif ($param->{action} eq "data2model")
  {
    check_data ($Data, "data"); 
    $Model=data2model ($Data);
    ($Data,$P)=decode ($Model, $Data, $param);
  }
elsif ($param->{action} eq "decode")
  {
    check_data ($Data, "data", $Model, "model");
    ($Data,$P)=decode ($Model, $Data, $param);
  }
elsif ($param->{evaluation} eq "sequence")
  {    
    check_data ($Data, "data", $Model, "model");
    $P = &seq2proba ($Model, $Data, $param);	  
  }
elsif  ($param->{action} eq "convert")
  {
    ;
  }
else
  {
    print STDERR "ERROR: $param->{action} is an unknown action \n";
  }

      

$param=set_output_name ($param, model2nst($Model),model2nem($Model));

    
print "\n";
print "----- rhhm.pl Output:-----\n";
if ($Data && $param->{outdata} ne "no")
    {
      printf "      Data Score  : log(P)=%.3f\n",$P;
      if ($param->{evaluate})
	{
	  data2recode_on_RST ($Data,$Model, "viterbi");
    my $sc1=data2evaluation  ($Data, $param->{evaluate});
    my $sc2=data2evaluation2 ($Data, $Model,$param->{evaluate});
    printf "      Data Pred   : sc1=%.2f sc2=%.2f\n", $sc1,$sc2;
	   	  
	}
      
      if (!$param->{output} || $param->{output} eq "data")
	{
	  dump_data ($Data, $param->{outdata},hash2string($param, "-"));
	}
      elsif ($param->{output} eq "fasta")
	{
	  $Model=long_state2short_state ($Model,$Data);
	  dump_seq ($Data, $param->{outdata});
	}
      printf "      Decoded Data: $param->{outdata}\n";
    }
if ($Model && $param->{outmodel} ne "no")
  {
    printf "      Model       : $param->{outmodel}\n";
    printf "      Constraints : $param->{outconstraints}\n\n";

    if ($param->{outmodel_format} eq "zeros")
      {
	dump_model_with_zero    ($Model, $param->{outmodel_with_zero}, hash2string($param, "-"));
      }
    elsif ($param->{outmodel_format} eq "default" || !$param->{outmodel_format})
      {
	dump_model      ($Model,$param->{outmodel}      ,hash2string($param, "-"));
      }
    dump_constraints($Model,$param->{outconstraints}, hash2string($param, "-"));
  }

if ($P && $param->{outEval} ne "no")
  {
    printf "   Forward eval: log (P)  $param->{outEval}\n";
    dump_evalForward ($P, $param->{outEval});
  }
  
exit (0);


####################################################################
#                                                                  #
#                                                                  #
#                                                                  #
#                    Parameters                                    #
#                                                                  #
#                                                                  #
####################################################################
sub set_output_name
  {
    my $param=shift;
    my $nst=shift;
    my $nem=shift;
    if (!$param->{out})
      {
	$param->{out}=$param->{data};
	$param->{out}=~s/\.[^\.]*$//;
	if (!$param->{out})
	  {
	    $param->{out}=$param->{model};
	    $param->{out}=~s/\.[^\.]*$//;
	  }
	if ($param->{out} eq "in.rhmm"){$param->{out}="out";}
      }
    $param->{out}=~s/\.rhmm//;
    
    if (!$param->{outdata})       {$param->{outdata}       ="$param->{out}.rhmm.$nst\_$nem.decoded";}
    if (!$param->{outmodel})      {$param->{outmodel}      ="$param->{out}.rhmm.$nst\_$nem.model";}
    if (!$param->{outEval})      {$param->{outEval}      ="$param->{out}.rhmm.$nst\_$nem.eval";}
    if (!$param->{outmodel_with_zero})      {$param->{outmodel_with_zero}      ="$param->{out}.rhmm.$nst\_$nem.model_zero";}
    if (!$param->{outconstraints}){$param->{outconstraints}="$param->{out}.rhmm.$nst\_$nem.constraints";}
    
    return $param;
  }
sub check_data
  {
    my @l=@_;

    for (my $a=0; $a<=$#l; $a+=2)
      {
	if ( !$l[$a])
	  {print "\n*****ERROR: could not use $l[$a+1] [FATAL]\n";
	   exit (1);
	 }
      }
    return 1;
  }
		      
sub check_parameters 
  {
    my $p=shift;
    my $rp={};

    $rp->{test}=1;
    $rp->{n}=1;
    $rp->{len}=1;
    $rp->{model}=1;
    $rp->{data}=1;
    $rp->{data2}=1;
    
    $rp->{constraints}=1;
    $rp->{action}=1;
    $rp->{nrounds}=1;
    $rp->{nit}=1;
    $rp->{outdata}=1;
    $rp->{outmodel}=1;
    $rp->{outEval}=1;
    $rp->{outmodel_format}=1;
    $rp->{evaluate}=1;
    $rp->{evaluation}=1;
    $rp->{output}=1;
    $rp->{out}=1;
    $rp->{single}=1;
    $rp->{nst}=1;
    $rp->{nem}=1;
    $rp->{w}=1;
    $rp->{gw}=1;
    
    $rp->{c}=1;
    
    foreach my $k (keys(%$p))
      {
	if (!$rp->{$k})
	  {
	    print "\n****ERROR: $k is an unknown pararmeter[FATAL]***\n";
	    die;
	  }
	else
	  {
	    print "PARAM: -$k ---> [$p->{$k}]\n";
	  }
      }
    return $p;
  }
  
sub process_param
  {
    my @arg=@_;
    my $cl=join(" ", @arg);
    
    my @commands=split (/\s\-+/,$cl);
    my $param={};
    
    
    foreach my $c (@commands)
      {
	if (!($c=~/\S/)){next;}
	$c=~/(\w+)\s*(.*)\s*/;
	my $k=$1;
	if (!$2){$param->{$k}=1;}
	else {$param->{$k}=$2;}
	$param->{$k}=~s/\s*$//;
      }
    
    return check_parameters ($param);
  }

####################################################################
#                                                                  #
#                                                                  #
#                                                                  #
#                    HMM Trainning                                 #
#                                                                  #
#                                                                  #
####################################################################
sub viterbi_trainningL
  {
    #untested
    my $M=shift;
    my $S=shift;
    
    my $A={};
    my $E={};
    my $ns;
    my $VP;
    
    ($S,$VP)=viterbiL($M,$S);
    
    foreach my $j (keys(%$S))
      {
	my $L=keys(%{$S->{$j}});
	for (my $i=1; $i<=$L; $i++)
	  {

	    my $l=$S->{$j}{$i}{viterbi};
	    if ($i>1)
	      {
		my $k=$S->{$j}{$i-1}{viterbi};
		$A->{$k}{$l}++;
	      }
	    $E->{$l}{$S->{$j}{$i}{bin}}++;
	  }
      }
    
    foreach my $k (keys (%$M))
      {
	foreach my $l (keys(%{$M->{$k}}))
	  {
	    if (($l=~/ST::/)){$A->{$k}{$l}+=1;}
	    else {$E->{$k}{$l}+=1;}
	  }
      }
    #update A/Model
    foreach my $k (keys (%$M))
      {
	my $num;
	foreach my $lp (keys(%$M)){$num+=$A->{$k}{$lp};}
	foreach my $l (keys(%$M)){$M->{$k}{$l}{v}=$A->{$k}{$l}/$num;}
      }
    
    # update E/model 
    foreach my $k (keys (%$M))
      {
	my $num;
	foreach my $lp (keys(%{$M->{$k}}))
	  {
	    if (!($lp =~/ST::/)){$num+=$E->{$k}{$lp};}
	  }
 	foreach my $l (keys(%{$M->{$k}}))
	  {
	    if (!($l =~/ST::/)){$M->{$k}{$l}{v}=$E->{$k}{$l}/$num;}
	  }
      }
    
    return (model2modelL($M), $VP);
  }
sub baum_welch
  {
    my $S=shift;
    my $M=shift;
    my $p=shift;
    if ($p->{c}) {return baum_welch_C($S,$M,$p);}
    else         {return baum_welch_P($S,$M,$p);}
  }

sub baum_welch_P
  {
    my $S=shift;
    my $M=shift;
    my $p=shift;
    
    my $BM;
    my $maxidle=5;
    my $score;
    my $delta;
    my ($P, $PP, $BP);
    

    for (my $i=1;$i<=$p->{nrounds}; $i++)
      {
	
	print "---- $i -----\n";
	
	$M=model2modelR ($M);
	#display_model ($M);
	$M=model2modelL($M);
	

	my $idle;
	my $cont=1;
	
	my $it;
	$PP=0;
	for ($it=0; $it<$p->{nit} && $idle<$maxidle && $cont; $it++)
	  {
	    ($M,$P)=baum_welch_iteration ($S,$M);
	    if ($PP)
	      {
		$delta=$P-$PP;
		if ($delta<0){$idle++;}
		elsif ($delta<0.001 || $delta<0){$idle++;}
		else {$idle=0;}
	      }
	    printf "\t$it ==> %.3f\n", $P;
	    $PP=$P;
	  }
	if (!$BM || $P>$BP)
	  {
	    $BM=$M;
	    $BP=$P;
	    dump_model ($BM, "current_best.rhmm.model", "SCORE=$PP");
	  }
      }
    
    
    
    ($S, $score)=viterbiL($BM, $S);
    
    posteriorL($BM, $S);
    
    print "\n----- Best Model: $BP";
    if ($p->{evaluate})
	{
	  
	  ($S,$M)=data2recode_on_RST ($S,$M,   $p->{evaluate});
	  my $sc1 =data2evaluation  ($S,    $p->{evaluate});
	  my $sc2 =data2evaluation2 ($S,$BM,$p->{evaluate});
	  printf " EVALUATION $p->{mode} sc1: %.3f sc2: %.3f",$sc1,$sc2;
	}
    print "\n";
    return ($S,$BM,$BP);
  }
sub baum_welch_C
  {
    my $D=shift;
    my $M=shift;
    my $p=shift;
    my ($P,$VP,$PP);
    my $mfile="$$.model.m4c";
    my $dfile="$$.data.d4c";
    my $cfile="$$.constraint.4dc";
    my $nmfile="$$.new.m4c";
    
    my $I=model2index($M);
    if ($p->{c} eq "1"){$p->{c}="rhmmC";}
    dump_data_C ($D, $I,$dfile);
    dump_model_C($M,$I,$mfile,$cfile);
    
    print "RUN $p->{c}\n";
    print "$p->{c} bw $dfile $mfile $nmfile $cfile $p->{nrounds} $p->{nit}\n";
    
    system ("$p->{c} bw $dfile $mfile $nmfile $cfile $p->{nrounds} $p->{nit}");

    if (!-e "$nmfile.model")
      {
	print "ERROR: $p->{c} is Not properly installed--- Will run the Perl version (100X slower)\n";
	return baum_welch_P($D,$M,$p);
      }
    
    ($M,$P)=undump_model_C  ($M,$I,"$nmfile.model");
    ($D,$PP)=undump_viterbi_C ($D,$I,"$nmfile.viterbi");
    ($D,$VP)=undump_posterior_C ($D,$I,"$nmfile.posterior");

    unlink(($mfile,$dfile,$nmfile,$cfile,"$nmfile.model","$nmfile.viterbi", "$nmfile.posterior"));
    
    return ($D,$M,$P);
  }
  
sub baum_welch_iteration
  {

    my $S=shift;
    my $M=shift;
    my $A={};
    my $E={};
    my $TP;
    my $new_proc=1; #use the new procedure
    
    foreach my $j (keys (%$S))
      {
	my $L=keys (%{$S->{$j}});
	my $F={};
	my $B={};
	my $P;
	
	($B,$P)=backwardL($M, $S->{$j});#log_space
	($F,$P)=forwardL ($M, $S->{$j});
		
	#Update A
	foreach my $k (keys (%$M))
	  {
	    foreach my $l(keys(%$M))
	      {
		$A->{$j}{$k}{$l}=$LOG_ZERO;
		if (!$A->{-1}{$k}{$l}){$A->{-1}{$k}{$l}=$LOG_ZERO;}
		
		for (my $i=1; $i<$L; $i++)
		  {
		    my $symbol=$S->{$j}{$i+1}{'bin'};
		    
		    if (!(exists ($M->{$l}{$symbol}{v}))){next;}
		    my $fo=$F->{$i}{$k};#log space
		    my $ba=$B->{$i+1}{$l};#log_space
		    my $tr=$M->{$k}{$l}{v}; #log_space
		    my $em=$M->{$l}{$symbol}{v};;#log_space
		    $A->{$j}{$k}{$l}=log_add ($A->{$j}{$k}{$l},log_multiply($fo,$tr,$em,$ba));
		  }
		$A->{$j}{$k}{$l}=log_divide($A->{$j}{$k}{$l},$P);
		$A->{-1}{$k}{$l}=log_add ($A->{-1}{$k}{$l},$A->{$j}{$k}{$l});
	      }
	  }
	#update Emissions
	foreach my $k (keys (%$M))
	  {
	    foreach my $b (keys (%{$M->{$k}}))
	      {

		if ($b=~/ST::/){next;}
		$E->{$j}{$k}{$b}=$LOG_ZERO;
		if ( ! $E->{-1}{$k}{$b}){$E->{-1}{$k}{$b}=$LOG_ZERO;}
		
		for (my $i=1; $i<=$L; $i++)
		  {
		    if ($S->{$j}{$i}{bin} eq $b)
		      {
			my $p=$E->{$j}{$k}{$b};
			my $q=log_multiply ($F->{$i}{$k},$B->{$i}{$k});
			$E->{$j}{$k}{$b}=log_add($p,$q);
		      }
		  }
		
		$E->{$j}{$k}{$b}=log_divide($E->{$j}{$k}{$b},$P);
		$E->{-1}{$k}{$b}=log_add($E->{-1}{$k}{$b},$E->{$j}{$k}{$b});
	      }
	  }
	$TP=(!$TP)?$P:log_multiply ($TP, $P);
      }
    
    #turn log counts into counts
    my $pseudocount=0;
    foreach my $k (keys(%$M))
      {
	foreach my $l (keys (%{$M->{$k}}))
	  {
	    if (($l=~/ST::/ )){$A->{-1}{$k}{$l}=exp($A->{-1}{$k}{$l})+$pseudocount;}
	    else {$E->{-1}{$k}{$l}=exp($E->{-1}{$k}{$l})+$pseudocount;}
	  }
      }
    if ($new_proc)
      {
	foreach my $s1(keys(%$M))
	  {
	    my ($ts,$te,$fs,$fe,$v);
	    foreach my $s2 (keys(%{$M->{$s1}}))
	      {
		$v=$M->{$s1}{$s2}{v}=($s2 =~/ST::/)?$A->{-1}{$s1}{$s2}:$E->{-1}{$s1}{$s2};
		if (!exists($M->{$s1}{$s2}{f}))
		  {
		    if ($s2=~/ST::/){$ts+=$v;}
		    else {$te+=$v;}
		  }
		else
		  {
		    if ($s2=~/ST::/){$fs+=$M->{$s1}{$s2}{f};}
		    else {$fe+=$M->{$s1}{$s2}{f};}
		  }
	      }
	    $ts/=($fs==1)?1:(1-$fs);
	    $te/=($fe==1)?1:(1-$fe);
	    
	    foreach my $s2 (keys(%{$M->{$s1}}))
	      {
		if (exists($M->{$s1}{$s2}{f})){$M->{$s1}{$s2}{v}=$M->{$s1}{$s2}{f};}
		elsif ($s2=~/ST::/)
		  {$M->{$s1}{$s2}{v}/=($ts==0)?1:$ts;}
		else
		  {$M->{$s1}{$s2}{v}/=($te==0)?1:$te;}
	      }
	  }
      }
    else
      {
	#update A/Model
	foreach my $k (keys (%$M))
	  {
	    my $num;
	    
	    
	    foreach my $l (keys(%$M)){$num+=$A->{-1}{$k}{$l};}
	    foreach my $l (keys(%$M)){$M->{$k}{$l}{v}=($num==0)?0:$A->{-1}{$k}{$l}/$num;}
	  }
	
	# update E/model 
	foreach my $k (keys (%$M))
	  {
	    my $num;
	    foreach my $l (keys(%{$M->{$k}}))
	      {
		if (!($l =~/ST::/)){$num+=$E->{-1}{$k}{$l};}
	      }
	    foreach my $l (keys(%{$M->{$k}}))
	      {
		if (!($l =~/ST::/)){$M->{$k}{$l}{v}=($num==0)?0:$E->{-1}{$k}{$l}/$num;}
	      }
	  }
      }
    $M=model2modelL($M);
    
    return ($M,$TP);
  }
####################################################################
#                                                                  #
#                                                                  #
#                                                                  #
#                      HMM Decoding                                 #
#                                                                  #
#                                                                  #
####################################################################    
sub seq2proba
  {
    my $M=shift;
    my $D= shift;
    my $p=shift;

    if ($p->{c}) 
      {
    	($P) = forward_C ($M,$D,$p);
      }
    else 
      {
        my $f;
        $P = seq2probaL ($M,$D);
        print "******** Result forward: log(P)= $P";
      } 
    $param->{outdata} = "no";
    return ($P); 
  }
  
sub seq2probaL
  {
    my $M=shift;
    my $S=shift;
    $M=model2modelL($M);
    
    my $TP;

    for my $j (keys(%$S))
      {
	     my ($f,$P)=forwardL($M, $S->{$j});

	     $TP+=$P;
      }
    return $TP;

  }

sub forwardL 
  {    
    my $M=shift;    
    my $S=shift;
    my $f={};
    my $P;
    my $L=keys(%$S);

    foreach my $k (keys(%$M)){$f->{0}{$k}=0;}
      
    for (my $i=1; $i<=$L; $i++)
	   {
	     foreach my $l (keys(%$M))
	       {

	         $f->{$i}{$l}=$LOG_ZERO;
	         my $emit=(!exists($M->{$l}{$S->{$i}{bin}}{v}))?$LOG_ZERO:$M->{$l}{$S->{$i}{bin}}{v};
	         foreach my $k (keys(%$M))
		        {
#		  my $v = log_multiply ($f->{$i-1}{$k},$M->{$k}{$l}{v}); print ";;;;;;$f->{$i-1}{$k};;;;;;;;$M->{$k}{$l}{v};;;;;;;;$v\n";#tag2del
		          $f->{$i}{$l}=log_add($f->{$i}{$l}, log_multiply ($f->{$i-1}{$k},$M->{$k}{$l}{v}));
		        }
	         $f->{$i}{$l}=log_multiply($f->{$i}{$l},$emit);

	       }
     }
 
      $P=$LOG_ZERO;
      foreach my $k (keys (%$M))
	     {
	       $P=log_add ($P, $f->{$L}{$k});
	     }
      
      return ($f,$P);
    }
    
#forward_C ($M,$D,$p)
sub forward_C
  {
    my $M=shift;
    my $D=shift;
    my $p=shift;
    my ($P,$VP,$PP, $pEvalFor);
    my $mfile="$$.model.m4c";
    my $dfile="$$.data.d4c";
    my $cfile="$$.constraint.4dc";
    my $nmfile="$$.new.m4c";
        
    my $I=model2index($M);
    if ($p->{c} eq "1"){$p->{c}="rhmmC";}
    dump_data_C ($D, $I,$dfile);
    dump_model_C($M,$I,$mfile,$cfile);
    
    print "RUN $p->{c}\n";
    print "$p->{c} ev $dfile $mfile $nmfile $cfile\n";
    
    system ("$p->{c} ev $dfile $mfile $nmfile $cfile\n");
    
    if (!-e "$nmfile.eval")
      {
	     print "ERROR: $p->{c} is Not properly installed--- Will run the Perl version (100X slower)\n";
	     $M=model2modelL($M);
	     return forwardL ($D,$M,$p);
      }
    
    $pEvalFor = undump_evalForward_C ("$nmfile.eval");

    unlink (($dfile,$mfile,$cfile, $nmfile,"$nmfile.eval"));
    return ($pEvalFor);
  }
    
sub posteriorL
    {
      my $M=shift;
      my $S=shift;
      my $P;
      
      foreach my $j (keys (%$S))
      {
	
	my $L=keys (%{$S->{$j}});
	my $F={};
	my $B={};
	
	($B,$P)=backwardL($M,$S->{$j});#log_space
	($F,$P)=forwardL ($M, $S->{$j});
	
	for (my $i=1; $i<=$L; $i++)
	  {
	    my $symbol=$S->{$j}{$i}{'bin'};
	    my $bpost_score;
	    my $bpost_k;
	    foreach my $k (keys (%$M))
	      {
		if (!(exists ($M->{$k}{$symbol}{v}))){next;}
		my $p=log_divide (log_multiply($F->{$i}{$k},$B->{$i}{$k}),$P);
		
		if (!$bpost_score || $p>$bpost_score)
		  {
		    $bpost_score=$p;
		    $bpost_k=$k;
		  }
	      }
	    $S->{$j}{$i}{posterior}=$bpost_k;
	    $S->{$j}{$i}{bpost_score}=myexp($bpost_score);
	  }
      }
   return ($S,0);
 }

sub backwardL 
    {
      my $M=shift;
      my $S=shift; 
      
      my $B={};
     
      my $P;
      my $L=keys (%$S);
      

      foreach my $k (keys(%$M))
	{
	  $B->{$L}{$k}=$LOG_ZERO;
	  foreach my $l(keys(%$M))
	    {
	      $B->{$L}{$k}=log_add($B->{$L}{$k}, $M->{$k}{$l}{v});
	    }
	}
      
      for (my $i=$L-1; $i>=1; $i--)
	{
	  foreach my $k (keys(%$M))
	    {
	      $B->{$i}{$k}=$LOG_ZERO;
	      
	      foreach my $l (keys(%$M))
		{
		  if (!exists($M->{$l}{$S->{$i+1}{bin}}{v})){next;}
		  my $x=$M->{$k}{$l}{v};
		  my $y=$M->{$l}{$S->{$i+1}{bin}}{v};
		  my $z=$B->{$i+1}{$l};
		  my $p=$B->{$i}{$k};
		  my $q=$x+$y+$z;
		  $B->{$i}{$k}=log_add ( $B->{$i}{$k}, log_multiply($x,$y,$z));
		}
	    }
	}
      
      $P=$LOG_ZERO;
      foreach my $l (keys(%$M))
	{
	  my $emission=$M->{$l}{$S->{1}{bin}}{v};
	  $P=log_add ($P, log_multiply($B->{1}{$l},$emission));
	}
      return ($B,$P);
    }
sub decode
  {
    my $M=shift;
    my $D= shift;
    my $p=shift;
    my $P;

    
    $Model=model2modelL($Model);
    
    if ( $p->{c}){($D,$P)=decode_C($M,$D,$p);}
    else {($D,$P)=decode_P($M,$D,$p);}
    $Model=modelL2model($Model);

    return ($D,$P);
  }
sub decode_P
  {
    my $M=shift;
    my $D= shift;
    my $p=shift;
    my $P;
    
    
    ($D, $P)=posteriorL($Model, $Data);
    ($D, $P)=viterbiL  ($Model, $Data);
    $Model=modelL2model($Model);
    
    return ($D, $P);
  }
sub decode_C
  {
    my $M=shift;
    my $D= shift;
    my $p=shift;
    my ($VP, $PP);
    
    my $mfile="$$.model.m4c";
    my $dfile="$$.data.d4c";
    my $nmfile="$$.new.m4c";
    
    if ($p->{c} eq "1"){$p->{c}="rhmmC";}

    my $I=model2index($M);
    dump_data_C ($D, $I,$dfile);
    dump_model_C($M, $I,$mfile);
    system ("$p->{c} decode $dfile $mfile $nmfile");
    
    if (!-e "$nmfile.viterbi")
      {
	print "ERROR: $p->{c} is Not properly installed--- Will run the Perl version (100X slower)\n";
	return baum_welch_P($D,$M,$p);
      }
    

    ($D,$PP)=undump_viterbi_C ($D,$I,"$nmfile.viterbi");
    ($D,$VP)=undump_posterior_C ($D,$I,"$nmfile.posterior");
    unlink (($dfile,$mfile,$nmfile,"$nmfile.viterbi","$nmfile.posterior"));
    return ($D,$PP);
  }
sub viterbiL 
  {
    my $M=shift;
    my $S= shift;
    my ($max_k, $ptr_k);
    
    foreach my $j (keys(%$S))
      {
	print "key: $j\n"
	my $L=keys(%{$S->{$j}});
	my $PTR={};
	my $V={};
	my ($path, $ppath);
	
	foreach my $k (keys(%$M)){$V->{0}{$k}=0;}
	
	for (my $i=1; $i<=$L; $i++)
	  {
	    my $symbol=$S->{$j}{$i}{bin};
	    
	    foreach my $l (keys (%$M))
	      {
		if ( !exists ($M->{$l}{$symbol}{v}))
		  {
		    $V->{$i}{$l}=$LOG_ZERO;
		    $PTR->{$i}{$l}=$LOG_ZERO;
		  }
		else
		  {
		    $max_k=$LOG_ZERO;
		    $ptr_k="";
		    foreach my $k (keys (%$M))
		      {
			my $v=log_multiply($V->{$i-1}{$k},$M->{$k}{$l}{v});
			if ($v>$max_k || $max_k==$LOG_ZERO)
			  {
			    $max_k=$v;
			    $ptr_k=$k;
			  }
		      }
		    $V->{$i}{$l}=log_multiply($M->{$l}{$S->{$j}{$i}{bin}}{v},$max_k);
		    $PTR->{$i}{$l}=$ptr_k;
		  }
	      }
	  }
	
	$max_k=$LOG_ZERO;
	$ptr_k="";
	foreach my $k (keys (%$M))
	  {
	    
	    my $vv=$V->{$L}{$k};
	    if ($vv>$max_k  || $max_k==$LOG_ZERO)
	      {
		$max_k=$vv;
		$ptr_k=$k;
	      }
	  }
	for (my $i=$L; $i>=1; $i--)
	  {
	    $S->{$j}{$i}{viterbi}=$ptr_k;
	    $ptr_k=$PTR->{$i}{$ptr_k};
	  }
	
      }
    return ($S,$max_k);
  }
####################################################################
#                                                                  #
#                                                                  #
#                                                                  #
#                      Occasionally Dishonnest Casino              #
#                                                                  #
#                                                                  #
####################################################################
sub ODHC2model
  {
    my $M={};

   
    $M->{'ST::f'} {'ST::f'  }{v}=0.95;
    $M->{'ST::f'} {'ST::u'}{v}=0.05;
    
    $M->{'ST::u'} {'ST::f'  }{v}=0.05;
    $M->{'ST::u'} {'ST::u'}{v}=0.95;
    
   
    $M->{'ST::f'} {'1'}{v}=1/6;
    $M->{'ST::f'} {'2'}{v}=1/6;
    $M->{'ST::f'} {'3'}{v}=1/6;
    $M->{'ST::f'} {'4'}{v}=1/6;
    $M->{'ST::f'} {'5'}{v}=1/6;
    $M->{'ST::f'} {'6'}{v}=1/6;
    
    $M->{'ST::u'} {'1'}{v}=1/10;
    $M->{'ST::u'} {'2'}{v}=1/10;
    $M->{'ST::u'} {'3'}{v}=1/10;
    $M->{'ST::u'} {'4'}{v}=1/10;
    $M->{'ST::u'} {'5'}{v}=1/10;
    $M->{'ST::u'} {'6'}{v}=5/10;
    

    
    $M=set_model_null_entries($M);
    return $M;
  }
####################################################################
#                                                                  #
#                                                                  #
#                                                                  #
#                      Data Conversion                             #
#                                                                  #
#                                                                  #
####################################################################
sub data2size
    {
      my $d=shift;
      my $L;
      
      foreach my $e (keys (%$d))
	{
	  $L+=keys(%{$d->{$e}});
	}
      return $L;
    }
sub data2field_list
  {
    my $d=shift;
    my $list={};
    foreach my $exp (keys(%$d))
	{
	  foreach my $r (keys (%{$d->{$exp}})) 
	    {
	      foreach my $k (keys (%{$d->{$exp}{$r}}))
		{
		  $list->{$k}=1;
		}
	    }
	}
    return $list;
  }
sub data2bin
  {
    my $d=shift;
    my $bin={};
    
    foreach my $exp (keys(%$d))
      {
	my $pst="";
	foreach my $r (keys (%{$d->{$exp}})) 
	  {
	    $bin->{$d->{$exp}{$r}{bin}}=1;
	  }
      }
    return $bin;
  }
sub data2model
    {
      my $d=shift;
      my $M={};
      my $pseudo=1;
      my %tot;
      
      if (!data_has_RST($d))
	{
	  print STDERR "ERROR: data has no RST field. Impossible to estimate model from data\n";
	  die;
	}
      elsif(!data_has_bin($d)) 
	{
	  print STDERR "ERROR: data is not binned. Impossible to estimate model from data\n";
	  die;
	}
      foreach my $exp (keys(%$d))
	{
	  my $pst="";
	  foreach my $r (keys (%{$d->{$exp}})) 
	    {
	      if ($d->{$exp}{$r}{RST})
		{
		  my $cst=$d->{$exp}{$r}{RST};
		  my $em= $d->{$exp}{$r}{bin};
		  
		  $M->{$cst}{$em}{v}++;
		  if ($pst){$M->{$pst}{$cst}{v}++;}
		  $tot{$cst}++;
		  $pst=$cst;
		}
	    }
	}
      foreach my $s1 (keys (%$M))
	{
	  foreach my $s2 (keys (%{$M->{$s1}}))
	    {
	      $M->{$s1}{$s2}{v}+=$pseudo;
	      $tot{$s1}+=$pseudo;
	      $M->{$s1}{$s2}{v}/=$tot{$s1};
	    }
	}
      
      return $M;
    }
sub data2recode_on_em
    {
      #recode data by comparing the viterbi assignments and the RST field.
      #returns recoded data and recoded model
      
      my $d=shift;
      my $M=shift;
      my $best={};
      my $nm={};
      my $convert={};
      
      foreach my $st (%$M)
	{
	  my $bv=$LOG_ZERO;
	  foreach my $em (keys (%{$M->{$st}}))
	    {
	      if ($em=~/ST::/){next;}
	      my $v=$M->{$st}{$em}{v};
	      if ($bv==-1 || $v>$bv)
		{
		  $best->{$st}=$em;
		  $bv=$v;
		}
	    }
	}
      my $n;
      foreach my $st (sort {$best->{$a}<=>$best->{$b}}keys(%$best))
	{
	  $n++;
	  print "\n----- $st: $best->{$st}  ---> ST::$n";
	  $convert->{$st}="ST::$n";
	}
      
      #relabel data
      foreach my $k1 (keys(%$d))
	{
	  foreach my $k2 (keys (%{$d->{$k1}}))
	    {
	      foreach my $k3 (keys (%{$d->{$k1}{$k2}}))
		{
		  my $v1=$d->{$k1}{$k2}{$k3};
		  my $v2=$convert->{$v1};
		  if ($k3 eq "viterbi" || $k3 eq "posterior")
		    {
		      $d->{$k1}{$k2}{$k3}=$v2;
		    }
		}
	    }
	}
    
      #relabel model
      foreach my $k1 (keys(%$M))
	{
	  foreach my $k2 (keys(%{$M->{$k1}}))
	    {
	      my $nk1=$convert->{$k1};
	      my $nk2=$convert->{$k2};
	      
	      $nk1=($nk1)?$nk1:$k1;
	      $nk2=($nk2)?$nk2:$k2;
	      
	      $nm->{$nk1}{$nk2}{v}=$M->{$k1}{$k2}{v};
	    }
	}
      return ($d,$nm);
    }
sub data2recode_on_RST
    {
      #recode data by comparing the viterbi assignments and the RST field.
      #returns recoded data and recoded model
      
      my $d=shift;
      my $m=shift;
      my $field=shift;
      my $nm={};
      my $convert={};
      my $best_score;
      my $r=[];
      my $border;
      
      if (!data_has_RST ($d)){return ($d,$m);}
      
      my @stl=keys(%$m);
      $r=permute(@stl);
      foreach my $order (@$r)
	{
	  my $lconvert={};
	  my $score;
	  
	  for (my $a=0; $a<=$#stl;$a++)
	    {
	      $lconvert->{$stl[$a]}=$order->[$a];
	    }
	  $score=data2evaluation ($d,$field, $lconvert);
	  if (!$best_score || $score>$best_score){$best_score=$score;$border=$order;}
	}
      
      for (my $a=0; $a<=$#stl;$a++)
	    {
	      $convert->{$stl[$a]}=$border->[$a];
	    }
      
      #relabel data
      foreach my $k1 (keys(%$d))
	{
	  foreach my $k2 (keys (%{$d->{$k1}}))
	    {
	      foreach my $k3 (keys (%{$d->{$k1}{$k2}}))
		{
		  my $v1=$d->{$k1}{$k2}{$k3};
		  my $v2=$convert->{$v1};
		  
		  if ($v2 && $k3 ne "RST")
		    {
		      #print "$v1 ---> $v2\n";
		      $d->{$k1}{$k2}{$k3}=$v2;
		    }
		}
	    }
	}
    
      #relabel model
      foreach my $k1 (keys(%$m))
	{
	  foreach my $k2 (keys(%{$m->{$k1}}))
	    {
	      my $nk1=$convert->{$k1};
	      my $nk2=$convert->{$k2};
	      
	      $nk1=($nk1)?$nk1:$k1;
	      $nk2=($nk2)?$nk2:$k2;
	      
	      $nm->{$nk1}{$nk2}{v}=$m->{$k1}{$k2}{v};
	    }
	}
      return ($d,$nm);
    }
sub display_field
  {
    my $d=shift;
    my $field=shift;
    my $sep=shift;

    foreach my $k1 (keys(%$d))
      {
	foreach my $k2 (keys (%{$d->{$k1}}))
	  {
	    print "$d->{$k1}{$k2}{$field}$sep";
	  }
      }
  }
sub data2evaluation
    {
      my $d=shift;
      my $field=shift;
      my $convert=shift;
      
      my ($tot,$score);
      
      if (!data_has_RST ($d)){return (0,0);}
      if ($field eq "1"){$field="viterbi";}
      
      foreach my $k1 (keys(%$d))
	{
	  foreach my $k2 (keys (%{$d->{$k1}}))
	    {
	      my $ref =$d->{$k1}{$k2}{RST};
	      my $pred=$d->{$k1}{$k2}{$field};
	      
	      if ($convert){$pred=$convert->{$pred};}
	      if ($ref eq $pred){$score++;}
	      $tot++;
	    }
	}

      $score/=($tot)?$tot:1;
      
      return $score;
    }
sub data2evaluation2
    {
      my $d=shift;
      my $m=shift;
      my $field=shift;
      my $r={};
      my $score;
      if (!data_has_RST ($d)){return 0;}
      if ($field eq "1"){$field="viterbi";}
      foreach my $st (keys (%$m))
	{
	  my ($tp,$tn,$fp,$fn);
	  foreach my $k1 (keys(%$d))
	    {
	      foreach my $k2 (keys (%{$d->{$k1}}))
		{
		  my $re =$d->{$k1}{$k2}{RST};
		  my $pr =$d->{$k1}{$k2}{$field};
		  
		  if    ($re eq $st && $pr eq $st){$tp++;}
		  elsif ($re eq $st && $pr ne $st){$fn++;}
		  elsif ($re ne $st && $pr ne $st){$tn++;}
		  elsif ($re ne $st && $pr eq $st){$fp++;}
		}
	    }

	  $r->{$st}{sn} =(($tp+$fn)==0)?0:($tp)/($tp+$fn);
	  $r->{$st}{sp} =(($tn+$fp)==0)?0:($tn)/($tn+$fp);
	  $r->{$st}{sn2}=(($tp+$fp)==0)?0:($tp)/($tp+$fp);
	}
      $score=100;
      foreach my $k1(keys(%$r))
	{
	  foreach my $k2 (keys (%{$r->{$k1}}))
	    {
	      my $v=$r->{$k1}{$k2};
	      $score=($v<$score)?$v:$score;
	    }
	}
      return $score;
    }

sub data_has_bin
      {
	my $d=shift;
	
	foreach my $k1 (keys(%$d))
	  {
	  foreach my $k2 (keys (%{$d->{$k1}}))
	    {
	      if (exists($d->{$k1}{$k2}{bin})){return 1;}
	      else {return 0;}
	    }
	 }
	return 1;
      }
sub data_has_RST
      {
	my $d=shift;
	
	foreach my $k1 (keys(%$d))
	  {
	  foreach my $k2 (keys (%{$d->{$k1}}))
	    {
	      if (exists($d->{$k1}{$k2}{RST})){return 1;}
	      else {return 0;}
	    }
	 }
	return 1;
      }
####################################################################
#                                                                  #
#                                                                  #
#                                                                  #
#                      Model Convertion                            #
#                                                                  #
#                                                                  #
####################################################################	
sub convert_state
  {
    my $m=shift;
    my $st1=shift;
    my $st2=split; 
    my $nm={};
    
    foreach my $k (keys (%$m))
      {
	foreach my $l (keys (%{$m->{$k}}))
	  {
	    my $nk=($k eq $st1)?$st2:$k;
	    my $nl=($l eq $st1)?$st2:$l;
	    
	    $nm->{$nk}{$nl}{v}=$m->{$k}{$l}{v};
	    if ( exists ( $m->{$k}{$l}{f})){$nm->{$nk}{$nl}{f}=$m->{$k}{$l}{f};}
	  }
      }
    return $nm;
  }


sub long_state2short_state
  {
    my $m=shift;
    my $d=shift;
    my @symbol=split (//,"1234567890abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ");
    
    my ($a,$b);
    my $b=0;
    my $nm={};
    my $stl={};
    
    if (model_has_short_state ($m)){return $m;}
    
    foreach my $k (keys(%$m))
      {
	if (!exists ($stl->{$k})){$stl->{$k}="ST::$symbol[$a++]";}
    	foreach my $l (keys(%{$m->{$k}}))
	  {
	    if (!($l=~/ST::/) && !exists ($stl->{$l})){$stl->{$l}="$symbol[$b++]";}
	  }
      }
    if ($b>=62 || $a>=62){error ("Too many states or emissions for single");}
    #print STDOUT "\nModel Conversion into Short State (for Fasta output)\n";
    #foreach my $k (keys (%$stl)){if ($k=~/ST::/){	print "$k -----> $stl->{$k}\n";}}
    #foreach my $k (keys (%$stl)){if (!$k=~/ST::/){	print "$k -----> $stl->{$k}\n";}}
    
    foreach my $k(keys(%$m))
      {
	foreach my $l (keys(%{$m->{$k}}))
	  {
	    my $nk=$stl->{$k};
	    my $nl=$stl->{$l};
	    $nm->{$nk}{$nl}{v}=$m->{$k}{$l}{v};
	    if ( exists ($m->{$k}{$l}{f})){$nm->{$nk}{$nl}{f}=$m->{$k}{$l}{f};}
	  }
      }
    if ($d)
      {
	my @fl=("bin","viterbi","RST","posterior");
	foreach my $exp (keys (%$d))
	  {
	    foreach my $rec (keys(%{$d->{$exp}}))
	      {
		foreach my $f (@fl)
		  {
		    if ( exists ($d->{$exp}{$rec}{$f}))
		      {
			my $v=$d->{$exp}{$rec}{$f};
			$v=$stl->{$v};
			$d->{$exp}{$rec}{$f}=$v;
		      }
		  }
	      }
	  }
      }
    

    return $nm;
  }
sub model_has_short_state
  {
    my $m=shift;
    foreach my $k (keys (%$m))
      {
	$k=~/ST::(.*)/;
	if ( length($1)>1){return 0;}
	foreach my $l (keys (%{$m->{$k}}))
	  {
	    if (!$l=~/ST::/)
	      {
		if (length($l)>1){return 0;}
	      }
	  }
      }
    return 1;
  }
sub model_is_Count
  {
    my $m=shift;
    foreach my $s1 (keys (%$m))
      {
	foreach my $s2 (keys (%{$m->{$s1}}))
	  {
	    if ($m->{$s1}{$s2}{v}>1){return 1;}
	  }
      }
    return 0;
  }
sub model_is_Log
  {
    my $m=shift;
    foreach my $s1 (keys (%$m))
      {
	foreach my $s2 (keys (%{$m->{$s1}}))
	  {
	    if ($m->{$s1}{$s2}{v}<0){return 1};
	  }
      }
    return 0;
  }
sub model2constraints
  {
    my $M=shift;
    my $mode=shift;
    
    my $stl=model2state_list($M);
    my $eml=model2emission_list($M);
    if ($mode eq "null" || $mode eq "1" || $mode eq "default")
      {
	foreach my $s1 (keys(%$stl))
	  {
	    foreach my $s2 (keys(%$stl))
	      {
		if (!exists($M->{$s1}{$s2}{v}) || $M->{$s1}{$s2}{v}==0)
		  {
		    $M->{$s1}{$s2}{f}=0;
		  }
	      }
	    foreach my $s2 (keys(%$eml))
	      {
		if (!exists($M->{$s1}{$s2}{v}) || $M->{$s1}{$s2}{v}==0)
		  {
		    $M->{$s1}{$s2}{f}=0;
		  }
	      }
	  }
      }
   
    return $M;
  }

sub set_model_null_entries
  {
    my $M=shift;
    my $sl=model2state_list($M);
    foreach my $s1 (keys(%$sl))
      {
	foreach my $s2 (keys(%$sl))
	  {
	    if (!exists ($M->{$s1}{$s2}{v})){$M->{$s1}{$s2}{v}=0;}
	  }
      }
    return $M;
  }

sub model2data
  {
    my $M=shift;
    my $n=shift; #number of strings
    my $l=shift; #string length
    my $S={};
    my $sl=model2state_list ($M);    
    my $state;
    my $symbol;
    
    my @sll=keys (%$sl);
    $state=shift(@sll);
    
    for (my $j=1; $j<=$n; $j++)
      {
	for (my $i=1; $i<=$l; $i++)
	  {
	    ($state, $symbol)=model2emit($M, $state);
	    
	    $S->{$j}{$i}{bin}=$symbol;
	    $S->{$j}{$i}{RST}=$state;
	   
	  }
      }
    return $S;
  }

sub model2emit
  {
    my $M=shift;
    my $start=shift;
    my ($state, $bin);
    my $r_state=rand(1);
    my $r_bin=rand(1);
    my $p=0;
    
    foreach my $k (keys(%$M))
      {
	$p+=$M->{$start}{$k}{v};
	if ( $r_state<=$p){$state=$k;last;}
      }
    
    $p=0;
    foreach my $bin (keys (%{$M->{$state}}))
      {
	if ( $bin =~/ST::/){;}
	else
	  {
	    $p+=$M->{$state}{$bin}{v};
	    if ($r_bin<=$p)
	      {
		return ($state, $bin);
	      }
	  }
      }
    
    return ($state, $bin);
  }


sub model2nst
  {
    my $M=shift;
    my $l=model2state_list($M);
    my @ll=keys (%$l);
    return $#ll+1;
  }
sub model2nem
  {
    my $M=shift;
    my $l=model2emission_list($M);
    my @ll=keys (%$l);
    return $#ll+1;
  }

sub model2state_list
  {
    my $M=shift;
    my $list={};

    foreach my $k1 (keys(%$M))
      {
	$list->{$k1}=1;
	foreach my $k2 (keys(%{$M->{$k1}}))
	  {
	    if ($k2=~/ST::/){$list->{$k2}=1;}
	  }
      }
    return $list;
  }

sub model2emission_list
  {
    my $M=shift;
    my $list={};

    foreach my $k1 (keys(%$M))
      {
	foreach my $k2 (keys(%{$M->{$k1}}))
	  {
	    if (!($k2=~/ST::/)){$list->{$k2}=1;}
	  }
      }
    return $list;
  }

sub model2index
  {
    my $M=shift;
    my $list={};
    my $n=0;
    my $ls=model2state_list($M);
    my $le=model2emission_list ($M);
    
    foreach my $k (keys(%$ls))
      {
	$list->{i2l}{$n++}=$k;
	$list->{ns}++;
	$list->{t}++;
      }
    
    foreach my $k (keys(%$le))
      {
	$list->{i2l}{$n++}=$k;
	$list->{ne}++;
	$list->{t}++;
      }
    
    for (my $a=0; $a<$n; $a++)
      {
	my $l=$list->{i2l}{$a};
	$list->{l2i}{$l}=$a;
      }
    return $list;
  }

sub undumped_model2model
  {
    my $M=shift;

    if (check_model ($M, "f")==0)
      {
	print STDERR "FATAL::rhmm.pl\n";
	die;
      }
    
    foreach my $s1 (keys (%$M))
      {
	my ($ts,$te,$fs,$fe,$v);
	foreach my $s2 (keys(%{$M->{$s1}}))
	  {
	   
	    $v=$M->{$s1}{$s2}{v};
	    if (!exists($M->{$s1}{$s2}{f}))
	      {
		if ($s2=~/ST::/){$ts+=$v;}
		else {$te+=$v;}
	      }
	    else
	      {
		if ($s2=~/ST::/){$fs+=$M->{$s1}{$s2}{f};}
		else {$fe+=$M->{$s1}{$s2}{f};}
	      }
	  }
	$ts/=($fs==1)?1:(1-$fs);
	$te/=($fe==1)?1:(1-$fe);
	
	foreach my $s2 (keys(%{$M->{$s1}}))
	  {
	    if (exists ($M->{$s1}{$s2}{f})){$M->{$s1}{$s2}{v}=$M->{$s1}{$s2}{f};}
	    elsif ($s2=~/ST::/)
	      {$M->{$s1}{$s2}{v}/=($ts==0)?1:$ts;}
	    else
	      {$M->{$s1}{$s2}{v}/=($te==0)?1:$te;}
	  }
      }
    
    return $M;
  }
sub check_model 
  {
    my $M=shift;
    my $f=shift;
    my $no_error=1;
    foreach my $s1 (keys (%$M))
      {
	my $te;
	my $ts;
	foreach my $s2 (keys(%{$M->{$s1}}))
	  {
	    if ($s2=~/ST::/){$ts+=$M->{$s1}{$s2}{$f};}
	    else {$te+=$M->{$s1}{$s2}{$f};}
	  }
	if ($te>1.0001)
	  {
	    print STDERR "\nERROR: emsissions ($f) proba of $s1 do not sum to 1: $te\n";
	    $no_error=0;
	  }
	if ($ts>1.0001)
	  {
	    print STDERR "\nERROR: transitions ($f) proba of $s1 do not sum to 1: $ts\n";
	    $no_error=0;
	  }
      }
    return $no_error;
  }

sub model2modelR
  {
    #randomize a model
    my $M=shift;
    
    foreach my $s1 (keys (%$M))
      {
	my ($ts,$te,$fs,$fe,$v);
	foreach my $s2 (keys(%{$M->{$s1}}))
	  {
	    $v=rand(10000)+1;
	    $M->{$s1}{$s2}{v}=$v;
	    
	    if (!exists ($M->{$s1}{$s2}{f}))
	      {
		if ($s2=~/ST::/){$ts+=$v;}
		else {$te+=$v;}
	      }
	    else
	      {
		if ($s2=~/ST::/){$fs+=$M->{$s1}{$s2}{f};}
		else {$fe+=$M->{$s1}{$s2}{f};}
	      }
	  }
	$ts/=($fs==1)?1:(1-$fs);
	$te/=($fe==1)?1:(1-$fe);
	
	foreach my $s2 (keys(%{$M->{$s1}}))
	  {
	    if (exists ($M->{$s1}{$s2}{f})){$M->{$s1}{$s2}{v}=$M->{$s1}{$s2}{f};}
	    elsif ($s2=~/ST::/)
	      {$M->{$s1}{$s2}{v}/=($ts==0)?1:$ts;}
	    else
	      {$M->{$s1}{$s2}{v}/=($te==0)?1:$te;}
	  }
      }
    return $M;
  }

sub model2modelL
  {
    #turns a proba model into a Log model
     my $M=shift;
     my $tag=shift;
     
     foreach my $k(keys(%$M))
       {
	 foreach my $l (keys(%$M))
	   {
	     $M->{$k}{$l}{v}=(!$M->{$k}{$l}{v} ||$M->{$k}{$l}{v}<$DELTA)?$LOG_ZERO:mylog($M->{$k}{$l}{v});
	   }
	}
      
      foreach my $k(keys(%$M))
	{
	  foreach my $l (keys(%{$M->{$k}}))
	    {
	      
	      if (!($l=~/ST::/))
		{
		  $M->{$k}{$l}{v}=(!$M->{$k}{$l}{v} ||$M->{$k}{$l}{v}<$DELTA)?$LOG_ZERO:mylog($M->{$k}{$l}{v});
		}
	    }
	}
     return $M;
    }
sub modelL2model
  {
     my $M=shift;
     my $tag=shift;
     
     foreach my $k(keys(%$M))
       {
	 foreach my $l (keys(%$M))
	   {
	     if ($M->{$k}{$l}{v}>0){return $M;}
	   }
	}
     
     foreach my $k(keys(%$M))
       {
	 foreach my $l (keys(%$M))
	   {
	     $M->{$k}{$l}{v}=(!$M->{$k}{$l}{v} ||$M->{$k}{$l}{v}==$LOG_ZERO)?0:exp($M->{$k}{$l}{v});
	   }
	}
     
      foreach my $k(keys(%$M))
	{
	  foreach my $l (keys(%{$M->{$k}}))
	    {
	      if (!($l=~/ST::/))
		{
		  $M->{$k}{$l}{v}=(!$M->{$k}{$l}{v} ||$M->{$k}{$l}{v}==$LOG_ZERO)?0:exp($M->{$k}{$l}{v});
		}
	    }
	}
      return $M;
    }

sub modelC2model
  {
    #turns a count model into a proba model
    my $m=shift;

    foreach my $s1 (keys (%$m))
      {
	my $tot_em;
	my $tot_st;
	foreach my $s2 (keys (%{$m->{$s1}}))
	  {
	    if ($s2=~/ST::/){$tot_st+=$m->{$s1}{$s2}{v};}
	    else {$tot_em+=$m->{$s1}{$s2}{v};}
	  }
	foreach my $s2 (keys (%{$m->{$s1}}))
	  {
	    if ($s2=~/ST::/){$m->{$s1}{$s2}{v}/=($tot_st)?$tot_st:1;}
	    else {$m->{$s1}{$s2}{v}/=($tot_em)?$tot_em:1;}
	  }
      }
    return $m;
  }

sub topology2model
  {
    my $data=shift;
    my $nst=shift;
    my $nem=shift;
    my $bin;
    
    my $M={};
    
    if ($nem)
      {
	$bin={};
	for (my $a=1; $a<=$nem; $a++){$bin->{$a}=1;}
      }
    elsif ($data) 
      {
	$bin=data2bin($data);
      }
    else 
      {
	print STDERR "ERROR: Number of Emission unknown";
	exit (0);
      }
    if (!$nst)
      {
	print STDERR "ERROR: Number of States unknown";
	exit (0);
      }
    for (my $a=1; $a<=$nst; $a++)
      {
	for (my $b=1; $b<=$nst; $b++)
	  {
	    $M->{"ST::$a"}{"ST::$b"}{v}=1;
	  }
	foreach my $b (keys (%$bin))
	  {
	    $M->{"ST::$a"}{$b}{v}=1;
	  }
      }
    $M=modelC2model($M);
    return $M;
    
  }
####################################################################
#                                                                  #
#                                                                  #
#                                                                  #
#                      Math                                        #
#                                                                  #
#                                                                  #
####################################################################
sub myexp
    {
      
      my $x=shift;
      if ( $x==$LOG_ZERO){return 0;}
      return exp($x);
    }
sub mylog
  {
   my $x=shift;
   if ( $x<$LOG_UNDERFLOW){return $LOG_ZERO;}
   return log($x);
 }

sub log_add 
  {
    my ($x,$y)=@_;
    
    
    if ($x==$LOG_ZERO){return $y;}
    elsif ($y==$LOG_ZERO){return $x;}
    elsif ($x>=$y)
      {
	return $x+mylog(1+myexp($y-$x));
      }
    else
      {
	return $y+mylog(1+myexp($x-$y));
      }
  }
sub log_divide 
  {
    my $x=shift;
    my $y=shift;
    
    if ($x==$LOG_ZERO || $y==$LOG_ZERO || $y==0){return $LOG_ZERO;}
    return $x-$y;
  }
sub log_multiply
  {
    my @l=@_;
    my $r;
    
    foreach my $v (@l)
      {
	if ($v==$LOG_ZERO){return $LOG_ZERO;}
	$r+=$v;
      }
    return $r;
  }

####################################################################
#                                                                  #
#                                                                  #
#                                                                  #
#                      Strings Mani                                #
#                                                                  #
#                                                                  #
####################################################################
sub array2hash

  {
    my $arrayR=shift;
    my $A=shift;
    my ($v, $k);
    my @array=@$arrayR;
    
    while (($k=shift(@array)))
      {
	my $v=shift (@array);
	$k=~s/-//g;
	$A->{$k}=$v;
      }
    return $A;
  }

sub hash2string
  {
    my $p=shift;
    my $sep=shift;
    my $cl;
    
    foreach my $k (keys (%$p))
      {
	$cl.="$sep$k $p->{$k} ";
      }
    return $cl;
  }
		   
sub string2hash 
  {
    my $s=shift;
    my $h=shift;

    my @l=split (/\s+/, $s);
    shift @l;
    return array2hash (\@l, $h);
  }
####################################################################
#                                                                  #
#                                                                  #
#                                                                  #
#                      Input/output  Data                          #
#                                                                  #
#                                                                  #
####################################################################
sub main_undump_data
  {
    my $d=shift;
    my $p=shift;
    
    if ($p->{data}) 
      {
	my @fl=split (/\s+/, $p->{data});
	foreach my $ff (@fl)
	  {
	    if ( -e $ff) {$d=generic_undump_data($ff,$d, $p);}
	  }
      }
    
    if ($param->{data2})
      {
	my @fl=split (/\s+/, $param->{data2});
	for (my $a=0; $a<=$#fl; $a+=2)
	  {
	    if ( -e $fl[$a]){$d=generic_undump_data($fl[$a],$d, $p,$fl[$a+1]);}
	  }
      }
    return $d;
  }
sub generic_undump_data
  {
    my $file=shift;
    my $d=shift;
    my $p=shift;
    my $field=shift;
    my $F= new FileHandle;

    vfopen ($F,$file);
    my $l=<$F>;
    close ($F);
    
    

    if    ($l=~/Format: rhmm.data.01/){return undump_data ($file,$d,$p,$field);}
    elsif ($l=~/^>/){return undump_seq ($file,$d,$p,$field);}
    else 
      {
	return undump_data ($file,$d,$p,$field);
	print STDERR "***** ERROR: Format of $file is unknown [FATAL]\n";
	die;
      }
  }
sub display_data
  {
    my $d=shift;
    my $file=shift;
    my $F= new FileHandle;
    
    if (!$file){open ($F, ">-");}
    else {vfopen ($F, ">$file");}
   
    print $F "$HEADER";
    foreach my $c (sort ({$a<=>$b}keys(%$d)))
      {
	foreach my $i (sort {$a<=>$b}keys (%{$d->{$c}}))
	  {
	    print $F "#d;";
	    foreach my $k (sort (keys (%{$d->{$c}{$i}})))
	      {
		print $F "$k;$d->{$c}{$i}{$k};";
	      }
	    print $F "\n";
	  }
      }
    close ($F);
  }
sub undump_data
      {
	my $file=shift;
	my $d=shift;
	my $F= new FileHandle;
	my $internalID;
	
	if (!$d){$d={};}
	vfopen ($F, $file);

	while (<$F>)
	  {
	    my $l=$_;
	   
	    chomp($l);
	    if ( $l=~/#d/)
	      {
		my @v=($l=~/([^;]+)/g);
		shift @v;	
		my $exp=shift (@v);
		my $record=shift(@v);
		$internalID+=1;
		for (my $a=0; $a<=$#v; $a+=2)
		  {
		    $d->{$exp}{$internalID}{$v[$a]}=$v[$a+1];
		  }
		$d->{$exp}{$internalID}{externalID}=$record; 
	      }
	  }
	close ($F);
	return $d;
      }
sub dump_data_C
  {
    my $d=shift;
    my $I=shift;
    my $file=shift;
    my $F=new FileHandle;

    
    vfopen ($F, ">$file");
    printf $F "%d", data2size($d);
    
    foreach my $exp (sort {$a<=>$b}keys(%$d))
      {
	foreach my $rec (sort {$a<=>$b}keys(%{$d->{$exp}}))
	  {
	    my $bin=$d->{$exp}{$rec}{bin};
	    $bin=$I->{l2i}{$bin};
	    print $F " $bin";
	  }
      }
    close ($F);
    return;
  }
sub undump_viterbi_C
  {
    my $d=shift;
    my $I=shift;
    my $file=shift;
    my $F=new FileHandle;
    
    vfopen ($F, "$file");
    my $cl=<$F>;
    close ($F);
    my @l=split (/\s+/, $cl);
    my $L=shift (@l);#get rid of L
    my $P=shift (@l);
    
    foreach my $exp (sort {$a<=>$b}keys(%$d))
      {
	foreach my $rec (sort {$a<=>$b}keys(%{$d->{$exp}}))
	  {
	    $d->{$exp}{$rec}{viterbi}=$I->{i2l}{shift(@l)};
	  }
      }
    return ($d,$P);
  }
sub undump_posterior_C
  {
    my $d=shift;
    my $I=shift;
    my $file=shift;
    my $F=new FileHandle;
    
    
    vfopen ($F, "$file");
    my $cl=<$F>;
    close ($F);
    
    
    my @l=split (/\s+/, $cl);
    my $L=shift (@l);#get rid of L
    my $P=shift (@l);
    
    foreach my $exp (sort {$a<=>$b}keys(%$d))
      {
	foreach my $rec (sort {$a<=>$b}keys(%{$d->{$exp}}))
	  {
	    $d->{$exp}{$rec}{posterior}=$I->{i2l}{shift(@l)};
	    $d->{$exp}{$rec}{bpost_score}=myexp(shift(@l));
	  }
      }
    return ($d,$P);
  }
  
sub undump_evalForward_C
  {

    my $file=shift;
    my $F=new FileHandle;
        
    vfopen ($F, "$file");
    my $cl=<$F>;
    close ($F);
    return ($cl);
  }  

sub dump_data
      {
	my $d=shift;
	my $file =shift;
	my $comment=shift;
	
	my $F= new FileHandle;

	vfopen ($F, ">$file");
	print $F "#comment;Format: rhmm.data.01\n";
	print $F "#comment;$comment;\n";
	foreach my $exp (sort {$a<=>$b}(keys(%$d)))
	  {
	    foreach my $record (sort {$a<=>$b}(keys(%{$d->{$exp}})))
	      {
		my $nrecord;
		if (exists($d->{$exp}{$record}{externalID}))
		    {
		      $nrecord=$d->{$exp}{$record}{externalID};
		    }
		else 
		  {
		    $nrecord=$record;
		  }
		
		
		print $F "#d;$exp;$nrecord;";
		foreach my $k (sort (keys (%{$d->{$exp}{$record}})))
		  {
		    
		    if ($k ne "externalID"){print $F "$k;$d->{$exp}{$record}{$k};";}
		  }
		print $F "\n";
	      }
	  }
	close ($F);
      }

sub undump_seq
  {
    my $file=shift;
    my $d=shift;
    my $p=shift;
    my $f=shift;
    my $F=new FileHandle;
    my ($w,$gw);
    $w=(exists($p->{w}))?$p->{w}:1;
    if ($p->{gw}){$w=$gw=$p->{gw};}
    
    
    foreach my $e ((file2entry($file,">")))
      {
	if ($e eq ">"){next;}
	my (@s)=split (/\n/, $e);
	my $h=shift(@s);
	my $s=join('',@s);
	$s=~s/>//g;
	
	my ($name, $f2)=split (/\s+/,$h);
	if (!$f){$f=$f2;}
	
	if (!$f)
	  {
	    print STDERR "ERROR: untaged data must be provided via -data2\n";
	    die;
	  }
	
	my @rl=split (//,$s);
	
	for (my $b=0;$b<=($#rl-($w-1)); $b++)
	  {
	    my $value;
	    if ($f eq "bin" && !$gw)  {for (my $c=0; $c<$w; $c++){$value.=$rl[$b+$c];}}
	    elsif ($f eq "bin" && $gw){$value=$rl[$b].$rl[$b+$gw];}
	    elsif ($f eq "RST"){$value="ST::$rl[$b]";}
	    
	    $d->{$name}{$b+1}{$f}=$value;
	  }
      }
    return $d;
  }

	
sub dump_seq
  {
    my $d=shift;
    my $f=shift;
    my @list=@_;
    my $F=new FileHandle;
    
    vfopen ($F, ">$f");
    my $fl=data2field_list($d);
    if ($#list==-1){@list=("RST","posterior","viterbi","bin");}
    
        
    foreach my $field (@list)
      {
	if (!$fl->{$field}){next;}
	else
	  {
	    foreach my $s (keys (%$d))
	      {
		print $F ">$s $field\n";
		foreach my $r (sort {$a<=>$b}keys (%{$d->{$s}}))
		  {
		    my $vv;
		    my $v=$d->{$s}{$r}{$field};
		    
		    if ($field eq "bin"){$v=~/^(.).*/;$vv=$1;}
		    else {$v=~/^ST::(.).*/;$vv=$1;}
		    print $F "$vv";
		  }
		print $F "\n";
	      }
	  }
      }
    close ($F);
  }

sub dump_evalForward
  {
	my $P =shift;
	my $file =shift;
	my $F= new FileHandle;

	vfopen ($F, ">$file");

	print $F "$P";

	print $F "\n";
	close ($F);
   }
      
####################################################################
#                                                                  #
#                                                                  #
#                                                                  #
#                      Input/output Model                          #
#                                                                  #
#                                                                  #
####################################################################
sub display_model
  {
    my $M=shift;
    
    print "#### STATES\n";
    foreach my $k(sort ({$a<=>$b}keys(%$M)))
      {
	foreach my $l (sort ({$a<=>$b} keys(%$M)))
	  {
	   if ($M->{$k}{$l}{v}){printf "$k;$l;%7.5f\n",$M->{$k}{$l}{v};}
	 }
      }
    print "#### EMISSIONS\n";
    foreach my $k(sort ({$a<=>$b}keys(%$M)))
      {
	foreach my $l (sort ({$a<=>$b}keys(%{$M->{$k}})))
	  {
	    if (!($l=~/ST::/))
	      {
		printf "$k;$l;%7.5f\n",$M->{$k}{$l}{v};
	      }
	  }
      }
    print "#### END\n";
    return;
  }
sub undump_model
      {

	#return a model with values set as proba
	my $file=shift;
	my $m={};
	my $C;
	my $L;
	my $F= new FileHandle;
	vfopen ($F, $file);
	while (<$F>)
	  {
	    my $l=$_;
	    chomp($l);
	    
	    if ($l=~/^#graph/)     {$m=add_graph($m,$l);}
	    elsif ($l=~/^#set/)    {$m=set_model_value ($m, $l);}
	    elsif ($l=~/^#bset/)   {$m=batch_set_model_value ($m, $l);}
	    elsif ($l=~/^#comment/){;}
	    elsif ($l=~/^#convert/)
	      {
		my @sl=split(/;/,$l);
		shift (@sl);
		$m=convert_state($m,shift(@sl), shift(@sl));
	      }
	    elsif (!$l=~/^#/){;}
	    else
	      {
		print STDERR "\n---ERROR: unknown model definition: \n$l\n";
		exit (0);
	      }
	    
	    
	  }
	close ($F);
	
	if    (model_is_Count($m)){$m=modelC2model($m);}
	elsif (model_is_Log($m)){$m=modelL2model($m);}
	$m=undumped_model2model($m);
	
	
	return $m
      }



sub add_graph
      {
	my $m=shift;
	my $l=shift;

	my @gg=split (/;/, $l);
	shift (@gg);
	my $gname=shift(@gg);
	my $nst=shift(@gg);
	my $nem=shift (@gg);
	my $bconnect=shift(@gg);
	my $sconnect=shift(@gg);
	my $fconnect=shift(@gg);
	my $topo=shift (@gg);
	my @sl;
	
	
	for (my $a=1; $a<=$nst; $a++)
	  {
	    push (@sl,"ST::$gname#$a");
	  }
	for (my $a=0; $a<$nst; $a++)
	  {
	    my $s1=$sl[$a];
	    for (my $b=0; $b<$nem; $b++){$m->{$s1}{$b+1}{v}=1;}
	    for (my $b=0; $b<$nst; $b++)
	      {
		my $s2=$sl[$b];
		if ($bconnect eq "full" && $b<$a)     {$m->{$s1}{$s2}{v}=1;}
		elsif ($bconnect eq "nb" && $b==$a-1) {$m->{$s1}{$s2}{v}=1;}
		
		if ($sconnect eq "full" && $b==$a){$m->{$s1}{$s2}{v}=1;}
		
		if ($fconnect eq "full" && $b>$a)     {$m->{$s1}{$s2}{v}=1;}
		elsif ($fconnect eq "nb" && $b==$a+1) {$m->{$s1}{$s2}{v}=1;}
	      }
	    if ($a==$nst-1 && $topo eq "circular")
	      {
		my $s2=$sl[0];
		
		if ($bconnect eq "full") {$m->{$s2}{$s1}{v}=1;}
		elsif ($bconnect eq "nb"){$m->{$s2}{$s1}{v}=1;}
	      
		if ($fconnect eq "full") {$m->{$s1}{$s2}{v}=1;}
		elsif ($fconnect eq "nb"){$m->{$s1}{$s2}{v}=1;}
	      }
	  }
	return $m;
      }

sub set_model_value
  {
    my $m=shift;
    my $l=shift;
    
    my @g=split (/;/, $l);
    my $type =shift (@g);
    my $st1  =shift(@g);
    my $st2  =shift(@g);
    my $v1   =shift(@g);
    my $mode =shift(@g);
    
    $m->{$st1}{$st2}{v}=$v1;
    if    ($mode eq "fixed"){ $m->{$st1}{$st2}{f}=$v1;}
    elsif ($mode eq "unfix"){ undef($m->{$st1}{$st2}{f});}
    return $m;
  }
  
sub batch_set_model_value
  {
    my $m=shift;
    my $l=shift;
    
    my @g=split (/;/, $l);
    my $type =shift (@g);
    my $st1  =shift(@g);
    my $st2  =shift(@g);
    my $v1   =shift(@g);
    my $mode =shift(@g);
    
    foreach my $k (keys (%$m))
      {
	if ($k=~/^$st1$/)
	  {
	    foreach my $l (keys(%{$m->{$k}}))
	      {
		if ($l=~/ST::/ && !($st2=~/ST::/)){next;}
		if ($l=~/^$st2$/)
		  {
		    $m->{$k}{$l}{v}=$v1;
		    if ($mode eq "fixed"){$m->{$k}{$l}{f}=$v1;}
		    elsif ($mode eq "unfixed"){undef ($m->{$k}{$l}{f});}
		    elsif ($mode)
		      {
			print STDERR "\n---ERROR: unknown bset mode:\n$l\n";
		      }
		  }
	      }
	  }
      }
    return $m;
  }
sub dump_constraints
      {
	my $model=shift;
	my $file =shift;
	my $comment =shift;

	my $F= new FileHandle;
	
	
	if ($file){vfopen ($F, ">$file");}
	else {$F=*STDOUT;}
	
	print $F "#comment;constraints\n";
	
	foreach my $k1 (sort(keys(%$model)))
	  {
	    foreach my $k2 (sort(keys (%{$model->{$k1}})))
	      {
		if ($k2=~/^ST/) 
		  {
		    if (exists($model->{$k1}{$k2}{f}))
		      {printf $F "#st;$k1;$k2;%.4f;\n",$model->{$k1}{$k2}{f};}
		  }
	      }
	  }
	foreach my $k1 (sort(keys(%$model)))
	  {
	    foreach my $k2 (sort{$a<=>$b}(keys (%{$model->{$k1}})))
	      {
		if (!($k2=~/^ST/)) 
		  {
		    if (exists($model->{$k1}{$k2}{f})){printf $F "#set;$k1;$k2;%.4f;\n",$model->{$k1}{$k2}{v};}
		  }
	      }
	  }
	if ($F!=*STDOUT){close ($F);}
      }

sub dump_model_with_zero
      {
	my $m=shift;
	my $file =shift;
	my $comment =shift;
	my $sl=model2state_list ($m);
	my $el=model2emission_list ($m);
	my $F= new FileHandle;
	
	
	if ($file){vfopen ($F, ">$file");}
	else {$F=*STDOUT;}
	
	print $F "#comment;$comment\n";
	
	foreach my $k1 (keys (%$sl))
	  {
	    foreach my $k2 (keys (%$sl))
	      {
		printf $F "#set;$k1;$k2;%.4f;",$m->{$k1}{$k2}{v};
		if (exists ($m->{$k1}{$k2}{f})){print$F "fixed;"}
		print $F "\n";
	      }
	  }
	foreach my $k1 (keys (%$sl))
	  {
	    foreach my $k2 (keys (%$el))
	      {
		printf $F "#set;$k1;$k2;%.4f;",$m->{$k1}{$k2}{v};
		if (exists ($m->{$k1}{$k2}{f})){print$F "fixed;"}
		print $F "\n";
	      }
	  }
	if ($F!=*STDOUT){close ($F);}
      }
sub dump_model
      {
	my $m=shift;
	my $file =shift;
	my $comment =shift;

	my $F= new FileHandle;
	
	
	if ($file){vfopen ($F, ">$file");}
	else {$F=*STDOUT;}
	
	print $F "#comment;$comment\n";
	
	foreach my $k1 (sort(keys(%$m)))
	  {
	    foreach my $k2 (sort(keys (%{$m->{$k1}})))
	      {
		if ($k2=~/^ST/ && $m->{$k1}{$k2}{v}>$DELTA)
		  {
		    printf $F "#set;$k1;$k2;%.4f;",$m->{$k1}{$k2}{v};
		    if (exists ($m->{$k1}{$k2}{f})){print$F "fixed;"}
		    print $F "\n";
		  }
	      }
	  }
	foreach my $k1 (sort(keys(%$m)))
	  {
	    foreach my $k2 (sort{$a<=>$b}(keys (%{$m->{$k1}})))
	      {
		if (!($k2=~/^ST/) && $m->{$k1}{$k2}{v}>$DELTA)
		  {
		    printf $F "#set;$k1;$k2;%.4f;",$m->{$k1}{$k2}{v};
		    if (exists ($m->{$k1}{$k2}{f})){print$F "fixed;"}
		    print $F "\n";
		  }
	      }
	  }
	if ($F!=*STDOUT){close ($F);}
      }
sub dump_model_C
  {
    my $M=shift;
    my $I=shift;
    my $mfile=shift;
    my $cfile=shift;
    
    my $F1=new FileHandle;
    my $F2=new FileHandle;
    
    vfopen ($F1, ">$mfile");
    vfopen ($F2, ">$cfile");
    print $F1 "0 $I->{ns} $I->{ne}";
    print $F2 "0 $I->{ns} $I->{ne}";
    
    for (my $k=0; $k<$I->{ns}; $k++)
      {
	for (my $l=0; $l<$I->{t}; $l++)
	  {
	    my $s1=$I->{i2l}{$k};
	    my $s2=$I->{i2l}{$l};
	    
	    my $v1=(exists ($M->{$s1}{$s2}{v}))?$M->{$s1}{$s2}{v}:$LOG_ZERO;
	    print $F1 " $v1";
	    
	    my $v2=(exists($M->{$s1}{$s2}{f}))?$M->{$s1}{$s2}{f}:$LOG_ZERO;
	    print $F2 " $v2";
	  }
      }
    close ($F1);
    close ($F2);
    
    return $M;
  }
sub undump_model_C
  {
    my $M=shift;
    my $I=shift;
    my $file=shift;
    my $F=new FileHandle;

    vfopen ($F, "$file");
    my $l=<$F>;
    close ($F);
    my @list=split (/\s+/, $l);
    my $P =shift(@list);
    my $ns=shift(@list);
    my $ne=shift(@list);
    
    for (my $k=0; $k<$ns; $k++)
      {
	for (my $l=0; $l<($ne+$ns); $l++)
	  {
	    my $s1=$I->{i2l}{$k};
	    my $s2=$I->{i2l}{$l};
	    $M->{$s1}{$s2}{v}=shift(@list);
	  }
      }
    return ($M,$P);
  }

####################################################################
#                                                                  #
#                                                                  #
#                                                                  #
#                      Misc                                        #
#                                                                  #
#                                                                  #
####################################################################
sub error
  {
    my $msg=shift;
    
    print STDERR "ERROR: $msg [FATAL]\n";
    exit (0);
  }

sub vfopen 
  {
    my $f=shift;
    my $file=shift;

    if (($file =~/^\>/) && !($file =~/^\>\>/ )){open ($f, $file); return $f;}
    elsif (($file =~/^\>\>(.*)/))
      {
	if (!-e $1){	print STDERR "\nERROR: $file does not exist [FATAL]\n";exit(1);}
      }
    elsif (!-e $file){	print STDERR "\nERROR: $file does not exist [FATAL]\n";exit(1);}
   
    open ($f,$file);
    return $f;
  }
	
sub file2entry
  {
    my $f=shift;
    my $sep=shift;
    my $F=new FileHandle;
    
    my $rsep=$/;
    if ($sep){$/=$sep;}
    
    vfopen ($F, "$f");
    my @e=<$F>;
    close ($F);
    return @e;
  }
sub permute
  {
    my @list=@_;
    my $nl=[];
    my $nk={};
    my $r=[];

    foreach my $e1 (@list){$nk->{$e1}=1;}
    
    $r=sub_permute ($#list+1,$nk,$nl,$r);
    return $r;
   }
    
sub sub_permute
  {
    my $max=shift;
    my $k=shift;
    my $l=shift;
    my $r=shift;
    
    foreach my $e1 (keys (%$k))
      {
	my $nk={};
	
	foreach my $e2  (keys (%$k)){if ($e2 ne $e1){$nk->{$e2}=1;}}
	
	@$l=(@$l,$e1);
	
	if (@$l==($max))
	  {
	    
	    push (@$r,[@$l]);
	    my $a=0;
	  }
	else
	  {
	    $r=sub_permute ($max,$nk,$l,$r);
	  }
	pop (@$l);
      }
    return $r;
  }

