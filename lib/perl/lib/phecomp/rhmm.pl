#!/usr/bin/env perl
# -filter field <field name> contains|equals|min|max value
# -bin
use HTTP::Date;
use strict;
use FileHandle;


srand();
my $TAG=0;
my $BIN=0;

my $index;
my $LOG_ZERO=-999999999;
my $LOG_UNDERFLOW=0.0000000000000000000001;
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
    print "  ................................Mode: 'data' estimates model from data if RST field is set.\n";
    print "  -constraints <file>,<mode>'.....File: input a model where every transition/emission become a seed.\n";
    print "  ................................Mode: 1 to forbid every undeclared transition/emission.\n";
    print "  -data: <file1 file2.. >,<mode>..File: input data from file(s).\n";
    print "  ................................Mode: 'model' to generate simulated data from the model -model.\n";
    print "  -train:      <mode> ............Mode: 'bw' to apply aum and welch training.\n";
    print "  -nrounds     <value>............Value; number of trainned models.\n";
    print "  -nit         <value>............Value: number of iterations in each round of trainning\n";
    print "  -outmodel    <file>.............File:  name of the output file containning the trainned model.\n";
    print "  -output      <mode>.............Mode: data or fasta.\n";
    print "  -outdata     <file>.............File:  name of the output file containning the decoded data.\n";
    print "  -evaluate    <mode>.............Mode: 'posterior' evaluate posterior decoding predictions.\n";
    print "  ................................Mode: 'viterbi'   evaluate viterbi path predictions,\n";
    print "  ................................NOTE: RST field must be set in data.\n";
    print "  ................................NOTE: RST field declares the reference sate of each data item.\n";
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
    print "  #st;<state label>;<state label>;proba;\n";
    print "  #em;<state label>;<emission bin>;proba;\n";
    print "  #comment;<free text>;\n";
    print "\n";
    print "Note: to declare constraints one simply needs to declare transitions\n";
    print "******************Contact **************************\n";
    print "contact: cedric.notredame\@gmail.com\n\n";
    
    die;
  }


my $param=process_param (@ARGV);
my $nrounds=1;
my $nit=100;
my ($Data,$Model,$Constraints,$Score);

#process Parameters
if (!$param->{nref}){$param->{nref}=1;}
if (!$param->{lenref}){$param->{lenref}=3000;}
if ($param-> {test})
  {
    $param->{model}="in.rhmm.model";
    my $Model=ODHC2model($Model);
    dump_model ($Model,$param->{model},"Reference_test");
    
    my $Data=model2data($Model, $param->{nref}, $param->{lenref});
    $param->{data}="in.rhmm.data";
    dump_data ($Data,$param->{data}, "ReferenceData");
  }

if ($param->{model} eq "data"){;}
elsif ( $param->{model} eq "topology"){$Model=topology2model($param->{nstate}, $param->{nem});}
elsif (-e $param->{model}){$Model=undump_model($param->{model});}

if ($param->{constraints} eq "1" )
  {
    #makes non specified transitions/emissions forbiden in the seed model
    $Constraints=model2constraints($Model, $param->{constraints});
  }
elsif    (-e $param->{constraints}){$Constraints=undump_model($param->{constraints});}

if ($param->{data} eq "model"){;}
elsif (-e $param->{data}) 
  {
    my @flist=split (/\s+/, $param->{data});
    foreach my $ff (@flist)
      {
	$Data=generic_undump_data($ff,$Data);
      }
  }



if ($param->{model} eq "data" && $param->{data} eq "model")
  {
    print STDERR "Error: you must provide some data or a model!\n";die;
  }
elsif ($param->{data} eq "model")
  {
    $Data=model2data ($param->{nref}, $param->{lenref});
  }
elsif ($param->{model} eq "data")
  {
    $Model=data2model ($Data);
  }


if ($param->{train} eq "bw")
  {
   
    if (!$param->{nrounds}){$param->{nrounds}=1;}
    if (!$param->{nit}){$param->{nit}=100;}
    ($Score, $Model)=&baum_welch($Data,$Model,$Constraints,$param->{nrounds},$param->{nit}, $param->{evaluate});
  }

($Data,$Score)=viterbiL ($Model, $Data);
($Data)=posteriorL($Model, $Data);

$Model=modelL2model($Model);


if (!$param->{outdata}){$param->{outdata}="out.rhmm.data";}
if (!$param->{outmodel}){$param->{outmodel}="out.rhmm.model";}
($Data,$Model)=data2recode ($Data,$Model, "viterbi");

print "\n----- output decoded data: $param->{outdata} [$param->{output}]\n";

if (!$param->{output} || $param->{output} eq "data")
  {dump_data ($Data, $param->{outdata},hash2string($param, "-"));}
elsif ($param->{output} eq "fasta")
  {
    dump_seq ($Data, $param->{outdata});
  }

print "\n----- output best model  : $param->{outmodel}\n";
dump_model($Model,$param->{outmodel}, hash2string($param, "-"));
#display_data ($Data);
if ($param->{evaluate})
  {
    my $sc1=data2evaluation  ($Data, $param->{evaluate});
    my $sc2=data2evaluation2 ($Data, $Model,$param->{evaluate});
    
    print "\nBest Model score: $sc1 $sc2\n";
  }
die;

####################################################################
#                                                                  #
#                                                                  #
#                                                                  #
#                    Parameters                                    #
#                                                                  #
#                                                                  #
####################################################################

sub check_parameters 
  {
    my $p=shift;
    my $rp={};

    $rp->{test}=1;
    $rp->{model}=1;
    $rp->{data}=1;
    $rp->{constraints}=1;
    $rp->{train}=1;
    $rp->{nrounds}=1;
    $rp->{nit}=1;
    $rp->{outdata}=1;
    $rp->{outmodel}=1;
    $rp->{evaluate}=1;
    $rp->{output}=1;
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
    
    my @commands=split (/\-+/,$cl);
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
    my $score;
    
    ($S,$score)=viterbiL($M,$S);
    
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
	foreach my $l (keys(%$M)){$M->{$k}{$l}=$A->{$k}{$l}/$num;}
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
	    if (!($l =~/ST::/)){$M->{$k}{$l}=$E->{$k}{$l}/$num;}
	  }
      }
    
    return model2modelL($M);
  }

sub baum_welch
  {
    my $S=shift;
    my $M=shift;
    my $C=shift;
    
    my $multi=shift;
    my $nit=shift;
    my $evaluate=shift;
    
    my $best_score;
    my $BM;
    my $maxidle=5;
    my $score;
    
    for (my $i=0;$i<$multi; $i++)
      {
	my ($idle, $PP);
	print "---- $i -----\n";
	my $P;
	$M=model2modelR ($M,$C);
	display_model ($M);
	$M=model2modelL($M);
	
	my $cont=1;
	for (my $it=0; $it<$nit && $idle<$maxidle && $cont; $it++)
	  {
	    my ($sc1,$sc2);
	    ($P,$M)=baum_welch_iteration ($S,$M);
	    if ($evaluate)
	      {
		($S, $score)=viterbiL($M, $S);
		$S=posteriorL($M, $S);
		($S,$M)=data2recode ($S,$M, $evaluate);
		$sc1 =data2evaluation  ($S, $evaluate);
		$sc2 =data2evaluation2 ($S,$M,$evaluate);
	      }
	    
	    if ($PP)
	      {
		my $delta=$P-$PP;
		if ($delta<-1){$cont=0;}
		elsif ($delta<0.001){$idle++;}
		else {$idle=0;}
	      }
	    printf "\t$it ==> %.3f", $P;
	    if ($evaluate){printf " sc1: %.3f sc2: %.3f",$sc1,$sc2;}
	    print "\n";
	    
	    $PP=$P;
	  }
	my $score=seq2probaL($M, $S);
	if (!$BM || $score>$best_score)
	  {
	    $BM=$M;
	    $best_score=$score;
	  }
      }
    
    
    ($S, $score)=viterbiL($BM, $S);
    $S=posteriorL($BM, $S);
    
    return ($best_score, $BM);
  }
sub baum_welch_iteration
  {

    my $S=shift;
    my $M=shift;
    my $A={};
    my $E={};
    my $P;

    
    foreach my $j (keys (%$S))
      {
	
	my $L=keys (%{$S->{$j}});
	my $F={};
	my $B={};
	
	my ($P,$B)=backwardL($M, $S->{$j});#log_space
	my ($P,$F)=forwardL ($M, $S->{$j});
	
	
	
	#Update A
	foreach my $k (keys (%$M))
	  {
	    foreach my $l(keys(%$M))
	      {
		$A->{$j}{$k}{$l}=$LOG_ZERO;
		if (!$A->{$k}{$l}){$A->{$k}{$l}=$LOG_ZERO;}
		
		for (my $i=1; $i<$L; $i++)
		  {
		    my $symbol=$S->{$j}{$i+1}{'bin'};
		    if (!(exists ($M->{$l}{$symbol}))){next;}
		    
		    my $fo=$F->{$i}{$k};#log space
		    my $ba=$B->{$i+1}{$l};#log_space
		    my $tr=$M->{$k}{$l}; #log_space
		    my $em=$M->{$l}{$symbol};;#log_space
		    $A->{$j}{$k}{$l}=log_add ($A->{$j}{$k}{$l},log_multiply($fo,$tr,$em,$ba));
		  }
		$A->{$j}{$k}{$l}=log_divide($A->{$j}{$k}{$l},$P);
		$A->{$k}{$l}=log_add ($A->{$k}{$l},$A->{$j}{$k}{$l});
	      }
	  }
	
	#update Emissions
	foreach my $k (keys (%$M))
	  {
	    foreach my $b (keys (%{$M->{$k}}))
	      {
		if ($b=~/ST::/){next;}
		$E->{$j}{$k}{$b}=$LOG_ZERO;
		if ( ! exists ($E->{$k}{$b})){$E->{$k}{$b}=$LOG_ZERO;}
		
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
		$E->{$k}{$b}=log_add($E->{$k}{$b},$E->{$j}{$k}{$b});
	      }
	  }
      }
    
    
    foreach my $k (keys(%$M))
      {
	foreach my $l (keys (%{$M->{$k}}))
	  {
	    if (($l=~/ST::/ )){$A->{$k}{$l}=exp($A->{$k}{$l});}
	    else {$E->{$k}{$l}=exp($E->{$k}{$l});}
	  }
      }
    
     #add pseudo-counts
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
	foreach my $l (keys(%$M)){$M->{$k}{$l}=$A->{$k}{$l}/$num;}
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
	    if (!($l =~/ST::/)){$M->{$k}{$l}=$E->{$k}{$l}/$num;}
	  }
      }
    $M=model2modelL($M);
    $P=seq2probaL($M,$S);
    return ($P,$M);
  }
####################################################################
#                                                                  #
#                                                                  #
#                                                                  #
#                      HMM Decoding                                 #
#                                                                  #
#                                                                  #
####################################################################    
sub seq2probaL
  {
    my $M=shift;
    my $S=shift;
    
    my $TP;

    for my $j (keys(%$S))
      {
	my ($P, $f)=forwardL($M, $S->{$j});
	$TP+=$P;
      }
    return $TP;
  }
sub posteriorL
    {
      my $M=shift;
      my $S=shift;
      
      foreach my $j (keys (%$S))
      {
	
	my $L=keys (%{$S->{$j}});
	my $F={};
	my $B={};
	
	my ($P,$B)=backwardL($M,$S->{$j});#log_space
	my ($P,$F)=forwardL ($M, $S->{$j});
	
	for (my $i=1; $i<=$L; $i++)
	  {
	    my $symbol=$S->{$j}{$i}{'bin'};
	    my $bpost_score;
	    my $bpost_k;
	    foreach my $k (keys (%$M))
	      {
		if (!(exists ($M->{$k}{$symbol}))){next;}
		my $p=log_divide (log_multiply($F->{$i}{$k},$B->{$i}{$k}),$P);
		#$S->{$j}{$i}{posterior}{$k}=$p;
		if (!$bpost_score || $p>$bpost_score)
		  {
		    $bpost_score=$p;
		    $bpost_k=$k;
		  }
	      }
	    $S->{$j}{$i}{posterior}=$bpost_k;
	    $S->{$j}{$i}{bpost_score}=$bpost_score;
	  }
      }
   return $S;
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
	      my $emit=(!exists($M->{$l}{$S->{$i}{bin}}))?$LOG_ZERO:$M->{$l}{$S->{$i}{bin}};
	      
	      foreach my $k (keys(%$M))
		{
		  $f->{$i}{$l}=log_add($f->{$i}{$l}, log_multiply ($f->{$i-1}{$k},$M->{$k}{$l}));
		  my $v1=myexp($f->{$i}{$l});
		  my $v2=myexp($f->{$i-1}{$k});
		  my $v3=myexp($M->{$k}{$l});
		  #print "\t----L: V1: $v1 V2: $v2 V3: $v3\n";
		}
	      my $v1=myexp($f->{$i}{$l});
	      my $v2=myexp($emit);
	      
	      $f->{$i}{$l}=log_multiply($f->{$i}{$l},$emit);
	      my $v3=myexp($f->{$i}{$l});
	      #print "----L: $i V1: $v1 V2: $v2 V3: $v3\n";
	      
	    }
	}

      
      $P=$LOG_ZERO;
      foreach my $k (keys (%$M))
	{
	  $P=log_add ($P, $f->{$L}{$k});
	}
      
      return ($P,$f);
    }

sub backwardL 
    {
      my $M=shift;
      my $S=shift; 
      
      my $B={};
     
      my $P;
      my $L=keys (%$S);
      

      foreach my $k (keys(%$M)){$B->{$L}{$k}=0;}

      for (my $i=$L-1; $i>=1; $i--)
	{
	  foreach my $k (keys(%$M))
	    {
	      $B->{$i}{$k}=$LOG_ZERO;
	      
	      foreach my $l (keys(%$M))
		{
		  if (!exists($M->{$l}{$S->{$i+1}{bin}})){next;}
		  my $x=$M->{$k}{$l};
		  my $y=$M->{$l}{$S->{$i+1}{bin}};
		  my $z=$B->{$i+1}{$l};
		  my $p=$B->{$i}{$k};
		  my $q=$x+$y+$z;
		  
		  $B->{$i}{$k}=log_add ( $B->{$i}{$k}, log_multiply($x,$y,$z));
		  
		}
	    }
	}
      return (0,$B);
    }
sub viterbiL 
  {
    my $M=shift;
    my $S= shift;
    my ($max_k, $ptr_k);
    
    foreach my $j (keys(%$S))
      {
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
		if ( !exists ($M->{$l}{$symbol}))
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
			my $v=log_multiply($V->{$i-1}{$k},$M->{$k}{$l});
			my $v1=myexp($V->{$i-1}{$k});
			my $v2=myexp($M->{$k}{$l});
			my $v3=myexp($v);
			#print "\t$k $l---> V1: $v1 * V2: $v2 = V3: $v3\n";
			if ($v>$max_k || $max_k==$LOG_ZERO)
			  {
			    $max_k=$v;
			    $ptr_k=$k;
			  }
		      }
		    $V->{$i}{$l}=log_multiply($M->{$l}{$S->{$j}{$i}{bin}},$max_k);
		    $PTR->{$i}{$l}=$ptr_k;
		    my $v=exp($max_k);
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

   
    $M->{'ST::fair'} {'ST::fair'  }=0.95;
    $M->{'ST::fair'} {'ST::unfair'}=0.05;
    
    $M->{'ST::unfair'} {'ST::fair'  }=0.05;
    $M->{'ST::unfair'} {'ST::unfair'}=0.95;
    
   
    $M->{'ST::fair'} {'1'}=1/6;
    $M->{'ST::fair'} {'2'}=1/6;
    $M->{'ST::fair'} {'3'}=1/6;
    $M->{'ST::fair'} {'4'}=1/6;
    $M->{'ST::fair'} {'5'}=1/6;
    $M->{'ST::fair'} {'6'}=1/6;
    
    $M->{'ST::unfair'} {'1'}=1/10;
    $M->{'ST::unfair'} {'2'}=1/10;
    $M->{'ST::unfair'} {'3'}=1/10;
    $M->{'ST::unfair'} {'4'}=1/10;
    $M->{'ST::unfair'} {'5'}=1/10;
    $M->{'ST::unfair'} {'6'}=5/10;
    

    
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
sub data2model
    {
      my $d=shift;
      my $M={};
      my $pseudo=1;
      my %tot;
      
      foreach my $exp (keys(%$d))
	{
	  my $pst="";
	  foreach my $r (keys (%{$d->{$exp}})) 
	    {
	      if ($d->{$exp}{$r}{RST})
		{
		  my $cst=$d->{$exp}{$r}{RST};
		  my $em= $d->{$exp}{$r}{bin};
		  $M->{$cst}{$em}++;
		  if ($pst){$M->{$pst}{$cst}++;}
		  $tot{$cst}++;
		  $pst=$cst;
		}
	    }
	}
      foreach my $s1 (keys (%$M))
	{
	  foreach my $s2 (keys (%{$M->{$s1}}))
	    {
	      $M->{$s1}{$s2}+=$pseudo;
	      $tot{$s1}+=$pseudo;
	      $M->{$s1}{$s2}/=$tot{$s1};
	    }
	}
      display_model ($M);die;
      return $M;
    }
sub data2recode
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
	      
	      $nm->{$nk1}{$nk2}=$m->{$k1}{$k2};
	    }
	}
      return ($d,$nm);
    }

sub data2evaluation
    {
      my $d=shift;
      my $field=shift;
      my $convert=shift;
      
      my ($tot,$score);

      if (!data_has_RST ($d)){return (0,0);}
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


sub data_has_RST
      {
	my $d=shift;
	
	foreach my $k1 (keys(%$d))
	  {
	  foreach my $k2 (keys (%{$d->{$k1}}))
	    {
	      if ($d->{$k1}{$k2}{RST}){return 1;}
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
sub model2constraints
  {
    my $M=shift;
    my $mode=shift;
    my $C;
    
    if ($mode eq "null" || $mode eq "1")
      {
	my @stl=model2state_list($M);
	my @eml=model2emission_list($M);

	foreach my $s1 (@stl)
	  {
	    foreach my $s2 (@stl)
	      {
		if (!defined($M->{$s1}{$s2}) || $M->{$s1}{$s2}==0)
		  {
		    $C->{$s1}{$s2}=0;
		  }
	      }
	    foreach my $s2 (@eml)
	      {
		if (!defined($M->{$s1}{$s2}) || $M->{$s1}{$s2}==0)
		  {
		    $C->{$s1}{$s2}=0;
		  }
	      }
	  }
	return $C;
      }
  return $C;
  }
	    
sub set_model_null_entries
  {
    my $M=shift;
    my @sl=model2state_list($M);
    foreach my $s1 (@sl)
      {
	foreach my $s2 (@sl)
	  {
	    if (!$M->{$s1}{$s2}){$M->{$s1}{$s2}=0;}
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
    my @sl=model2state_list ($M);    
    my $state=$sl[0];
    my $symbol;
    
    for (my $j=0; $j<$n; $j++)
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
	$p+=$M->{$start}{$k};
	if ( $r_state<=$p){$state=$k;last;}
      }
    
    $p=0;
    foreach my $bin (keys (%{$M->{$state}}))
      {
	if ( $bin =~/ST::/){;}
	else
	  {
	    $p+=$M->{$state}{$bin};
	    if ($r_bin<=$p)
	      {
		return ($state, $bin);
	      }
	  }
      }
    
    return ($state, $bin);
  }



sub model2state_list
  {
    my $M=shift;
    my %list;

    foreach my $k1 (keys(%$M))
      {
	$list{$k1}=1;
	foreach my $k2 (keys(%{$M->{$k1}}))
	  {
	    if ($k2=~/ST::/){$list{$k2}=1;}
	  }
      }
    return keys (%list);
  }
sub model2emission_list
  {
    my $M=shift;
    my %list;

    foreach my $k1 (keys(%$M))
      {
	foreach my $k2 (keys(%{$M->{$k1}}))
	  {
	    if (!($k2=~/ST::/)){$list{$k2}=1;}
	  }
      }
    return keys (%list);
  }
sub model2modelR
  {
    #randomize a model
    my $M=shift;
    my $C=shift; #constraints imposed on some of the model transitions
    foreach my $s1 (keys (%$M))
      {
	my ($tr,$em,$v);
	foreach my $s2 (keys(%{$M->{$s1}}))
	  {
	    $v=rand(10000)+1;
	    $M->{$s1}{$s2}=$v;
	    
	    if (!defined ($C->{$s1}{$s2}))
	      {
		if ($s2=~/ST::/){$tr+=$v;}
		else {$em+=$v;}
	      }
	  }
	
	foreach my $s2 (keys(%{$M->{$s1}}))
	  {
	    if (defined ($C->{$s1}{$s2})){$M->{$s1}{$s2}=$C->{$s1}{$s2};}
	    elsif ($s2=~/ST::/)
	      {$M->{$s1}{$s2}/=($tr==0)?1:$tr;}
	    else
	      {$M->{$s1}{$s2}/=($em==0)?1:$em;}
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
	     $M->{$k}{$l}=(!$M->{$k}{$l} ||$M->{$k}{$l}<0.00000001)?$LOG_ZERO:mylog($M->{$k}{$l});
	   }
	}
      
      foreach my $k(keys(%$M))
	{
	  foreach my $l (keys(%{$M->{$k}}))
	    {
	      
	      if (!($l=~/ST::/))
		{
		  $M->{$k}{$l}=(!$M->{$k}{$l} ||$M->{$k}{$l}<0.00000001)?$LOG_ZERO:mylog($M->{$k}{$l});
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
	     if ($M->{$k}{$l}>0){return $M;}
	   }
	}
     
     foreach my $k(keys(%$M))
       {
	 foreach my $l (keys(%$M))
	   {
	     $M->{$k}{$l}=(!$M->{$k}{$l} ||$M->{$k}{$l}==$LOG_ZERO)?0:exp($M->{$k}{$l});
	   }
	}
     
      foreach my $k(keys(%$M))
	{
	  foreach my $l (keys(%{$M->{$k}}))
	    {
	      if (!($l=~/ST::/))
		{
		  $M->{$k}{$l}=(!$M->{$k}{$l} ||$M->{$k}{$l}==$LOG_ZERO)?0:exp($M->{$k}{$l});
		}
	    }
	}
      return $M;
    }

sub modelC2model
  {
    #turns a count model into a proba model
    my $m=shift;
    my $tot_em;
    my $tot_st;
    
    foreach my $s1 (keys (%$m))
      {
	foreach my $s2 (keys (%{$m->{$s1}}))
	  {
	    if ($s2=~/ST::/){$tot_st+=$m->{$s1}{$s2};}
	    else {$tot_em+=$m->{$s1}{$s2};}
	  }
	foreach my $s2 (keys (%{$m->{$s1}}))
	  {
	    if ($s2=~/ST::/){$m->{$s1}{$s2}/=($tot_st)?$tot_st:1;}
	    else {$m->{$s1}{$s2}/=($tot_em)?$tot_em:1;}
	  }
      }
    return $m;
  }

sub topology2model
  {
    my ($model,$nst, $nem)=@_;
    my $SL={};
    my $M={};
    

    my @emlist=split(/\s+/,$nem);
    if ($#emlist!=($nst-1))
      {
	for (my$a=0; $a<$nst; $a++){$emlist[$a]=$emlist[0];}
      }
    
    for (my $a=1; $a<=$nst; $a++)
      {
	$SL->{"ST::$a"}=$emlist[$a];
      }
    
    foreach my $st1 (keys (%$SL))
      {
	my $tot;
	
	#set emmissions
	for (my $a=1; $a<=$SL->{$st1}; $a++)
	  {
	    $tot+=$M->{$st1}{$a}=($model eq "rand")?rand (1000):100;
	  }
	for (my $a=1; $a<=$SL->{$st1}; $a++)
	  {
	    $M->{$st1}{$a}/=$tot;
	  }


	#Set Transitions
	$tot=0;
	foreach my $st2 (keys (%$SL))
	  {
	    $tot+=$M->{$st1}{$st2}=($model eq "rand")?rand (1000):100;
	  }
	foreach my $st2 (keys (%$SL))
	  {
	    $M->{$st1}{$st2}/=$tot;
	  }
      }
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
	$x=(($x==$LOG_ZERO) || ($y-$x)>=$LOG_UNDERFLOW)?$y:mylog($y-$x)+$x;
      }
    else{return log_add ($y,$x);}
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
sub generic_undump_data
  {
    my $file=shift;
    my $d=shift;
    my $F= new FileHandle;

    open ($F,$file);
    my $l=<$F>;
    close ($F);
    
    if    ($l=~/Format: rhmm.data.01/){return undump_data ($file,$d);}
    elsif ($l=~/^>/){return undump_seq ($file,$d);}
    else 
      {
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
    else {open ($F, ">$file");}
   
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
	
	if (!$d){$d={};}
	open ($F, $file);
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
		for (my $a=0; $a<=$#v; $a+=2)
		  {
		    $d->{$exp}{$record}{$v[$a]}=$v[$a+1];
		  }
	      }
	  }
	close ($F);
	return $d;
      }
sub dump_data
      {
	my $d=shift;
	my $file =shift;
	my $comment=shift;
	
	my $F= new FileHandle;

	open ($F, ">$file");
	print $F "#comment;Format: rhmm.data.01\n";
	print $F "#comment;$comment;\n";
	foreach my $exp (sort {$a<=>$b}(keys(%$d)))
	  {
	    foreach my $record (sort {$a<=>$b}(keys(%{$d->{$exp}})))
	      {
		print $F "#d;$exp;$record;";
		foreach my $k (sort (keys (%{$d->{$exp}{$record}})))
		  {
		    
		    print $F "$k;$d->{$exp}{$record}{$k};";
		  }
		print $F "\n";
	      }
	  }
	close ($F);
      }

sub undump_seq
  {
    my $f=shift;
    my $d=shift;
    my (@seq, @field, @name,$s);
    my $F=new FileHandle;
    
    open ($F, "$f");
    while (<$F>)
      {
	$s.=$_;
      }
    close ($F);
    
    
    @name=($s=~/>(\S*).*\n[^>]*/g);
    @seq =($s=~/>.*.*\n([^>]*)/g);
    @field =($s=~/>\S*(.*)\n([^>]*)/g);
   
    
    for ($a=0; $a<=$#seq; $a++)
      {
	my @rl=split (//,$seq[$a]);
	my $add=($field[$a] eq "bin")?"":"ST::";
	
	for (my $b=0;$b<=$#rl; $b++)
	  {
	    $d->{$name[$a]}{$b+1}{$field[$a]}="$add$rl[$b]";
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
    
    open ($F, ">$f");
    
    my $fl=data2field_list($d);
    if (!@list){@list=("RST","posterior","viterbi","bin");}
    
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
    foreach my $k(keys(%$M))
      {
	foreach my $l (keys(%$M))
	  {
	    printf "$k;$l;%7.5f\n",$M->{$k}{$l};
	  }
      }
    print "#### EMISSIONS\n";
    foreach my $k(keys(%$M))
      {
	foreach my $l (keys(%{$M->{$k}}))
	  {
	    if (!($l=~/ST::/))
	      {
		printf "$k;$l;%7.5f\n",$M->{$k}{$l};
	      }
	  }
      }
    print "#### END\n";
    return;
  }
sub undump_model
      {
	my $file=shift;
	my $m={};
	my $C;
	my $F= new FileHandle;
	open ($F, $file);
	while (<$F>)
	  {
	    my $l=$_;
	    
	    chomp($l);
	    if ( $l=~/^#st/ || $l=~/^#em/ )
	      {
		my @v=($l=~/([^;]+)/g);
		$m->{$v[1]}{$v[2]}=$v[3];
		if ($v[3]>1){$C=1;}
	      }
	    
	  }
	close ($F);
	if ($C){$m=modelC2model($m);}
	return $m
      }

sub dump_model
      {
	my $model=shift;
	my $file =shift;
	my $comment =shift;

	my $F= new FileHandle;

	

	open ($F, ">$file");
	print $F "#comment;$comment\n";
	
	foreach my $k1 (sort(keys(%$model)))
	  {
	    foreach my $k2 (sort(keys (%{$model->{$k1}})))
	      {
		if ($k2=~/^ST/) 
		  {
		    printf $F "#st;$k1;$k2;%.4f;\n",$model->{$k1}{$k2};
		  }
	      }
	  }
	foreach my $k1 (sort(keys(%$model)))
	  {
	    foreach my $k2 (sort(keys (%{$model->{$k1}})))
	      {
		if (!($k2=~/^ST/)) 
		  {
		    printf $F "#em;$k1;$k2;%.4f;\n",$model->{$k1}{$k2};
		  }
	      }
	  }
	close ($F);
      }

####################################################################
#                                                                  #
#                                                                  #
#                                                                  #
#                      Misc                                        #
#                                                                  #
#                                                                  #
####################################################################
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
####################################################################
#                                                                  #
#                                                                  #
#                                                                  #
#                      Documentation                               #
#                                                                  #
#                                                                  #
####################################################################	
#
# sample model
#
# rhmm -test -model <model_file> -data <data_file> -constraints <model_file> -train <trinmode> -nrounds <number of rounds> -niter <number of iterations> -outmodel <file> -outdata <file>
#
#
###############     -test -train bw -evaluate  #############################
#
# will run the occasionally dishonnest casino and produce two sample files:
#            out.rhmm.model: a model
#            out.rhmm.data : sample data produced by the model
#
#
###############     -model <model file>    #############################
#
#st;<state>;<state>;proba;
#em;<state>;<emission bin>;proba;
# example:
#st;ST::fair;ST::fair;0.9515;
#st;ST::fair;ST::unfair;0.0485;
#st;ST::unfair;ST::fair;0.0486;
#st;ST::unfair;ST::unfair;0.9514;
#em;ST::fair;1;0.1613;
#em;ST::fair;2;0.1948;
#em;ST::fair;3;0.1570;
#em;ST::fair;4;0.1422;
#em;ST::fair;5;0.1620;
#em;ST::fair;6;0.1827;
#em;ST::unfair;1;0.1141;
#em;ST::unfair;2;0.0843;
#em;ST::unfair;3;0.0857;
#em;ST::unfair;4;0.1053;
#em;ST::unfair;5;0.1087;
#em;ST::unfair;6;0.5019;

#
###############     -data <data file>    #############################
#

#d;<experiment><rec identifier>;bin;<emission_bin>
#
#after decoding, the data will also contain extra fields:
#d;<exp>;<rec id>;RST;<reference state(for evaluation);bin;<emission bin>;bpost_score;score of the best posterior decoding;posterior;<best posterior decoding state>;viterbi;<viterbi decoding label>;
#example
#d;0;36;RST;ST::unfair;bin;4;bpost_score;-0.283877115663927;posterior;ST::unfair;viterbi;ST::unfair;
#d;0;37;RST;ST::unfair;bin;6;bpost_score;-0.177523146348904;posterior;ST::unfair;viterbi;ST::unfair;
#d;0;38;RST;ST::unfair;bin;1;bpost_score;-0.14957962835706;posterior;ST::unfair;viterbi;ST::unfair;
#d;0;39;RST;ST::unfair;bin;6;bpost_score;-0.0991868631435864;posterior;ST::unfair;viterbi;ST::unfair;

###############     -constraints <model file>    #############################
#
# every item declred in this file will be used a s aconstraint when trainning the model
# it means that the values declared here will be used to initialize tzhe seed models in the bw trainning
#
#st;<state>;<state>;proba;
#em;<state>;<emission bin>;proba;
# example:
#st;ST::fair;ST::fair;0.3;
#
#


# runs baum and welch on the occasionaly dishonnest casino
# outputs out.rhmm.data and out.rhmm.model
# 

