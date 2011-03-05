#!/usr/bin/env perl
use HTTP::Date;
use strict;
use FileHandle;
my $p=process_param (@ARGV);
my $d1=undump_data ($p->{f1});
my $d2=undump_data ($p->{f2});
data2diff_list ($d1, $d2, $p);






####################################################################
#                                                                  #
#                                                                  #
#                                                                  #
#                    Process_data                                  #
#                                                                  #
#                                                                  #
####################################################################
sub display 
  {
    my $d1=shift;
    
    foreach my $e (keys(%$d1))
      {
	foreach my $r (sort{$a<=>$b}keys (%{$d1->{$e}}))
	  {
	    my $v1=$d1->{$e}{$r}{viterbi};
	    print "$e $r $v1\n";
	  }
      }
  }  
sub data2diff_list
  {
    my $d1=shift;
    my $d2=shift;
    my $p=shift;
    
    my ($start,$len,$end);
    foreach my $e (keys(%$d1))
      {
	foreach my $r (sort{$a<=>$b}keys (%{$d1->{$e}}))
	  {
	    my $v1=$d1->{$e}{$r}{viterbi};
	    my $v2=$d2->{$e}{$r}{viterbi};
	    
	    $v1=~s/ST:://;
	    $v2=~s/ST:://;
	    my $delta=abs(($v1-$v2));
	    if ($delta>=$p->{t})
	      {
		if ($start==-1){$start=$r;};
		$end=$r;
		$len++;
	      }
	    else
	      {
		if ($len>1)
		  {
		    printf "$e\t%10d %10d  %10d\n", $start, $end, $len;
		  }
		$len=0;
		$start=-1;
	      }
	  }
      }
		 
  }	       

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

    $rp->{f1}=1;
    $rp->{f2}=1;
    $rp->{t}=1;
        
    if (!$p)
      {
	foreach my $k (keys (%$rp))
	  {
	    print STDERR "-$k <value> ";
	  }
	print "\n";
	die;
      }
    
    foreach my $k (keys(%$p))
      {
	if (!$rp->{$k})
	  {
	    print STDERR "\n****ERROR: $k is an unknown pararmeter[FATAL]***\n";
	    die;
	  }

      }
    return $p;
  }
sub display_param
  {
    my $p=shift;
    my $F=shift;
    
    print $F "************** PARAM::START *****************\n";
    
    foreach my $v (keys(%$p))
      {
	print $F "PARAM: $v ---> $p->{$v}\n";
      }
    printf $F "Command line:\ndata2bin.pl ";
    foreach my $v (keys(%$p))
      {
	print $F "-$v $p->{$v} ";
      }
    print $F "\n************** PARAM::START *****************\n";
  }
sub process_param
  {
    my @arg=@_;
    my $cl=join(" ", @arg);
    
    my @commands=split (/\s\-+/,$cl);
    my $param={};
    if ($cl)
      {
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
    else
      {
	check_parameters();
      }
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
	return undump_data ($file,$d);
	print STDERR "***** ERROR: Format of $file is unknown [FATAL]\n";
	die;
      }
  }

sub undump_data
      {
	my $file=shift;
	my $d=shift;
	my $F= new FileHandle;
	my $internalID=0;
	
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
		$internalID+=1;
		for (my $a=0; $a<=$#v; $a+=2)
		  {
		    #print "[$v[$a]] -> [$v[$a+1]]\n";
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

    
    open ($F, ">$file");
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
    
    open ($F, "$file");
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
    
    
    open ($F, "$file");
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
	    $d->{$exp}{$rec}{bpost_score}=shift(@l);
	  }
      }
    return ($d,$P);
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
