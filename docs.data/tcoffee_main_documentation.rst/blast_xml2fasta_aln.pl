#!/usr/bin/env perl
use Env;


$tmp_dir="";
$init_dir="";
$program="tc_generic_method.pl";

$blast=@ARGV[0];

$name="query";$seq="";
%p=blast_xml2profile($name,$seq,100, 0, 0, $blast);
&output_profile (%p);


sub output_profile
  {
    my (%profile)=(@_);
    my ($a);
    for ($a=0; $a<$profile{n}; $a++)
      {
	
	print ">$profile{$a}{name} $profile{$a}{comment}\n$profile{$a}{seq}\n";
      }
    return;
  }
sub file_contains 
  {
    my ($file, $tag, $max)=(@_);
    my ($n);
    $n=0;
    
    if ( !-e $file && ($file =~/$tag/)) {return 1;}
    elsif ( !-e $file){return 0;}
    else 
      {
	open (FC, "$file");
	while ( <FC>)
	  {
	    if ( ($_=~/$tag/))
	      {
		close (FC);
		return 1;
	      }
	    elsif ($max && $n>$max)
	      {
		close (FC);
		return 0;
	      }
	    $n++;
	  }
      }
    close (FC);
    return 0;
  }
	    
	  
sub file2string
  {
    my $f=@_[0];
    my $string, $l;
    open (F,"$f");
    while (<F>)
      {

	$l=$_;
	#chomp ($l);
	$string.=$l;
      }
    close (F);
    $string=~s/\r\n//g;
    $string=~s/\n//g;
    return $string;
  }



sub tag2value 
  {
    
    my $tag=(@_[0]);
    my $word=(@_[1]);
    my $return;
    
    $tag=~/$word="([^"]+)"/;
    $return=$1;
    return $return;
  }
      
sub hit_tag2pdbid
  {
    my $tag=(@_[0]);
    my $pdbid;
       
    $tag=~/id="(\S+)"/;
    $pdbid=$1;
    $pdbid=~s/_//;
    return $pdbid;
  }
sub id2pdbid 
  {
    my $id=@_[0];
  
    if ($id =~/pdb/)
      {
	$id=~/pdb(.*)/;
	$id=$1;
      }
    $id=~s/[|¦_]//g;
    return $id;
  }
sub set_blast_type 
  {
    my $file =@_[0];
    if (&file_contains ($file,"EBIApplicationResult",100)){$BLAST_TYPE="EBI";}
    elsif (&file_contains ($file,"NCBI_BlastOutput",100)) {$BLAST_TYPE="NCBI";}
    else
      {
	$BLAST_TYPE="";
      }
    return $BLAST_TYPE;
  }
sub blast_xml2profile 
  {
    my ($name,$seq,$maxid, $minid, $mincov, $file)=(@_);
    my (%p, $a, $string, $n);
    


    if ($BLAST_TYPE eq "EBI" || &file_contains ($file,"EBIApplicationResult",100)){%p=ebi_blast_xml2profile(@_);}
    elsif ($BLAST_TYPE eq "NCBI" || &file_contains ($file,"NCBI_BlastOutput",100)){%p=ncbi_blast_xml2profile(@_);}
    else 
      {
	print "************ ERROR: Blast Returned an unknown XML Format **********************";
	die;
      }
    for ($a=0; $a<$p{n}; $a++)
      {
	my $name=$p{$a}{name};
	$p{$name}{seq}=$p{$a}{seq};
      }
    return %p;
  }
sub ncbi_blast_xml2profile 
  {
    my ($name,$seq,$maxid, $minid, $mincov, $string)=(@_);
    my ($L,$l, $a,$b,$c,$d,$nhits,@identifyerL);
    
    
    $seq=~s/[^a-zA-Z]//g;
    $L=length ($seq);
    
    %hit=&xml2tag_list ($string, "Hit");
    
    
    for ($nhits=0,$a=0; $a<$hit{n}; $a++)
      {
	my ($ldb,$id, $identity, $expectation, $start, $end, $coverage, $r);
	my (%ID,%DE,%HSP);
	
	$ldb="";

	%ID=&xml2tag_list ($hit{$a}{body}, "Hit_id");
	$identifyer=$ID{0}{body};
	
	%DE=&xml2tag_list ($hit{$a}{body}, "Hit_def");
	$definition=$DE{0}{body};
	
	%HSP=&xml2tag_list ($hit{$a}{body}, "Hsp");
	for ($b=0; $b<$HSP{n}; $b++)
	  {
	    my (%START,%END,%E,%I,%Q,%M);

	 
	    %START=&xml2tag_list ($HSP{$b}{body}, "Hsp_query-from");
	    %HSTART=&xml2tag_list ($HSP{$b}{body}, "Hsp_hit-from");
	    
	    %LEN=  &xml2tag_list ($HSP{$b}{body}, "Hsp_align-len");
	    %END=  &xml2tag_list ($HSP{$b}{body}, "Hsp_query-to");
	    %HEND=  &xml2tag_list ($HSP{$b}{body}, "Hsp_hit-to");
	    %E=&xml2tag_list     ($HSP{$b}{body}, "Hsp_evalue");
	    %I=&xml2tag_list     ($HSP{$b}{body}, "Hsp_identity");
	    %Q=&xml2tag_list     ($HSP{$b}{body}, "Hsp_qseq");
	    %M=&xml2tag_list     ($HSP{$b}{body}, "Hsp_hseq");
	    
	    for ($e=0; $e<$Q{n}; $e++)

	      {
		$qs=$Q{$e}{body};
		$ms=$M{$e}{body};
		if ($seq eq""){$seq=$qs;$L=length($seq);}
		
		$expectation=$E{$e}{body};
		$identity=($LEN{$e}{body}==0)?0:$I{$e}{body}/$LEN{$e}{body}*100;
		$start=$START{$e}{body};
		$end=$END{$e}{body};
		$Hstart=$HSTART{$e}{body};
		$Hend=$HEND{$e}{body};
	
		$coverage=(($end-$start)*100)/$L;

	
		if ($identity>$maxid || $identity<$minid || $coverage<$mincov){next;}
		@lr1=(split (//,$qs));
		@lr2=(split (//,$ms));
		$l=$#lr1+1;
		for ($c=0;$c<$L;$c++){$p[$nhits][$c]="-";}
		for ($d=0,$c=0; $c<$l; $c++)
		  {
		    $r=$lr1[$c];
		    if ( $r=~/[A-Za-z]/)
		      {
			
			$p[$nhits][$d + $start-1]=$lr2[$c];
			$d++;
		      }
		  }
		$Qseq[$nhits]=$qs;
		$Hseq[$nhits]=$ms;
		$QstartL[$nhits]=$start;
		$HstartL[$nhits]=$Hstart;
		$identityL[$nhits]=$identity;
		$endL[$nhits]=$end;
		$definitionL[$nhits]=$definition;
		$identifyerL[$nhits]=$identifyer;
		$comment[$nhits]="$ldb|$identifyer [Eval=$expectation][id=$identity%][start=$Hstart end=$Hend]";
		$nhits++;
	      }
	  }
      }
    
    $profile{n}=0;
    $profile{$profile{n}}{name}=$name;
    $profile{$profile{n}}{seq}=$seq;
    $profile {n}++;
    
    for ($a=0; $a<$nhits; $a++)
      {
	$n=$a+1;
	
	$profile{$n}{name}="$name\_$a";
	$profile{$n}{seq}="";
	$profile{$n}{Qseq}=$Qseq[$a];
	$profile{$n}{Hseq}=$Hseq[$a];
	$profile{$n}{Qstart}=$QstartL[$a];
	$profile{$n}{Hstart}=$HstartL[$a];
	$profile{$n}{identity}=$identityL[$a];
	$profile{$n}{definition}=$definitionL[$a];
	$profile{$n}{identifyer}=$identifyerL[$a];
	$profile{$n}{comment}=$comment[$a];
	for ($b=0; $b<$L; $b++)
	  {
	    if ($p[$a][$b])
	      {
		$profile{$n}{seq}.=$p[$a][$b];
	      }
	    else
	      {
		$profile{$n}{seq}.="-";
	      }
	  }
      }
    
    $profile{n}=$nhits+1;
    return %profile;
  }
sub ebi_blast_xml2profile 
  {
    my ($name,$seq,$maxid, $minid, $mincov, $string)=(@_);
    my ($L,$l, $a,$b,$c,$d,$nhits,@identifyerL,$identifyer);
    

    
    $seq=~s/[^a-zA-Z]//g;
    $L=length ($seq);
    %hit=&xml2tag_list ($string, "hit");
    
    for ($nhits=0,$a=0; $a<$hit{n}; $a++)
      {
	my ($ldb,$id, $identity, $expectation, $start, $end, $coverage, $r);
	my (%Q,%M,%E,%I);
	
	$ldb=&tag2value ($hit{$a}{open}, "database");
	$identifyer=&tag2value ($hit{$a}{open}, "id");

	$description=&tag2value ($hit{$a}{open}, "description");
	
	%Q=&xml2tag_list ($hit{$a}{body}, "querySeq");
	%M=&xml2tag_list ($hit{$a}{body}, "matchSeq");
	%E=&xml2tag_list ($hit{$a}{body}, "expectation");
	%I=&xml2tag_list ($hit{$a}{body}, "identity");
	

	for ($b=0; $b<$Q{n}; $b++)
	  {
	    
	    
	    $qs=$Q{$b}{body};
	    $ms=$M{$b}{body};
	    if ($seq eq""){$seq=$qs;$L=length($seq);}

	    $expectation=$E{$b}{body};
	    $identity=$I{$b}{body};
	    
	    	    
	    $start=&tag2value ($Q{$b}{open}, "start");
	    $end=&tag2value ($Q{$b}{open}, "end");
	    $startM=&tag2value ($M{$b}{open}, "start");
	    $endM=&tag2value ($M{$b}{open}, "end");
	    $coverage=(($end-$start)*100)/$L;
	    
	   # print "$id: ID: $identity COV: $coverage [$start $end]\n";
	    
	    
	    if ($identity>$maxid || $identity<$minid || $coverage<$mincov){next;}
	    # print "KEEP\n";

	    
	    @lr1=(split (//,$qs));
	    @lr2=(split (//,$ms));
	    $l=$#lr1+1;
	    for ($c=0;$c<$L;$c++){$p[$nhits][$c]="-";}
	    for ($d=0,$c=0; $c<$l; $c++)
	      {
		$r=$lr1[$c];
		if ( $r=~/[A-Za-z]/)
		  {
		    
		    $p[$nhits][$d + $start-1]=$lr2[$c];
		    $d++;
		  }
	      }
	  
	    
	    $identifyerL[$nhits]=$identifyer;
	    $comment[$nhits]="$ldb|$identifyer [Eval=$expectation][id=$identity%][start=$startM end=$endM]";
	    $nhits++;
	  }
      }
    
    $profile{n}=0;
    $profile{$profile{n}}{name}=$name;
    $profile{$profile{n}}{seq}=$seq;
    $profile {n}++;
    
    for ($a=0; $a<$nhits; $a++)
      {
	$n=$a+1;
	$profile{$n}{name}="$name\_$a";
	$profile{$n}{seq}="";
	$profile{$n}{identifyer}=$identifyerL[$a];
	
	$profile{$n}{comment}=$comment[$a];
	for ($b=0; $b<$L; $b++)
	  {
	    if ($p[$a][$b])
	      {
		$profile{$n}{seq}.=$p[$a][$b];
	      }
	    else
	      {
		$profile{$n}{seq}.="-";
	      }
	  }
      }
    $profile{n}=$nhits+1;
    
    return %profile;
  }

sub blast_xml2hit_list
  {
    my $string=(@_[0]);
    return &xml2tag_list ($string, "hit");
  }
sub xml2tag_list  
  {
    my ($string_in,$tag)=@_;
    my $tag_in, $tag_out;
    my %tag;
    
    if (-e $string_in)
      {
	$string=&file2string ($string_in);
      }
    else
      {
	$string=$string_in;
      }
    $tag_in1="<$tag ";
    $tag_in2="<$tag>";
    $tag_out="/$tag>";
    $string=~s/>/>##1/g;
    $string=~s/</##2</g;
    $string=~s/##1/<#/g;
    $string=~s/##2/#>/g;
    @l=($string=~/(\<[^>]+\>)/g);
    $tag{n}=0;
    $in=0;$n=-1;
  
 

    foreach $t (@l)
      {

	$t=~s/<#//;
	$t=~s/#>//;
	
	if ( $t=~/$tag_in1/ || $t=~/$tag_in2/)
	  {
	 
	    $in=1;
	    $tag{$tag{n}}{open}=$t;
	    $n++;
	    
	  }
	elsif ($t=~/$tag_out/)
	  {
	    

	    $tag{$tag{n}}{close}=$t;
	    $tag{n}++;
	    $in=0;
	  }
	elsif ($in)
	  {
	   
	    $tag{$tag{n}}{body}.=$t;
	  }
      }
  
    return %tag;
  }





