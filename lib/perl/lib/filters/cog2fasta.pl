#!/usr/bin/perl
$whog="whog.txt";
$myva="myva.txt";

if ( $ARGV[0] eq "-cog2fasta")
  {
    %cog=&cog2slist($whog, $ARGV[1]);
    $r=&slist2fasta ($myva, %cog);
    print $r;
    exit (EXIT_SUCCESS);
  }
elsif ( $ARGV[0] eq "-coglist2fasta")
  {
    &coglist2fasta ($whog, $myva,$ARGV[1]);
  }
elsif ( $ARGV[0] eq "-cogfile2coglist")
  {
    @list= &cogfile2coglist ($whog);
    print @list;
  }
elsif ($ARGV[0] eq "-compile_cog")
  {
    print stderr "\nCompile a collection of consistent COGs";
    &compile_cog($whog);
  }
elsif ($ARGV[0] eq "-cog2ccog")
  {
    @list=&cogfile2coglist ($whog);
    foreach $c (@list)
      {
	%cog=&cog2slist ( $whog, $ARGV[1]);
      }
  }
sub coglist2fasta
  {
    my $whog=$_[0];
    my $myva=$_[1];
    my $result;
    my $coglistfile=@_[2];
    my @coglist;
    

    open (F, $coglistfile);
    while (<F>)
      {
	$l.=$_;
      }
    close F;
    
    @coglist=($l=~/(COG0\S+)/g);
    
    foreach $cog (@coglist)
      {
	print "\nWrite $cog.fasta\n";
	%slist=&cog2slist($whog, $cog);
	$result=&slist2fasta ($myva, %slist);
	open ( F, ">$cog.fasta");
	print F "$result";
	close (F);
      }
    return;
    }
sub cogfile2coglist
  {
    my $cogfile=@_[0];
    my @coglist;
    
    open (F,$cogfile);
    while (<F>)
      {
	if (/(COG\S+)/)
	{
	  @coglist=(@coglist, $1);
	}
      }
    close (F);
    return @coglist;
  }

sub compile_cog 
  {
    my $whog=shift (@_);
    my @list, $line, %glist;
    open (F, $whog);
    print STDERR "\nRead COG lists from $whog\n";
    while (<F>)
      {
	
	if ( $_=~/(COG\S+)/ )
	  {
	    $cog=$1;
	    $glist{$cog}=$cog;
	    $glist{'LCOG'}{$glist{'NCOG'}++}=$cog;
	    while (<F> )
	      {
		if ( $_=~/___/)
		  {
		    last;
		  }
		elsif ( !/\:/){;}
		else
		  {
		    $line=$_;
		    $line=~s/:/ /g;
		    @list=($line=~/(\S+)/g);
		    $glist{$cog}{$list[0]}=$#list;
		    $glist{$cog}{$list[1]}{'genome'}=$list[0];
		    $glist{$cog}{$list[1]}{'seq'}=$list[1];
		    $glist{$cog}{$list[1]}{'north'}=$#list;
		    if ( !$glist{$list[0]})
		      {
			$glist{'LG'}{$glist{'NG'}++}=$list[0];
			$glist{$list[0]}=1;
		      }
		  }
	      }
	  }
      }
    close (F);

    print STDERR "\nIdentify Consistent Groups";
    foreach ($a=0; $a<$glist{'NCOG'}; $a++)
      {
	$cog=$glist{'LCOG'}{$a};
	for ( $b=0; $b<$glist{'NG'}; $b++)
	  {
	    $genome=$glist{'LG'}{$b};
	    $v=($glist{$cog}{$genome})?"1":"0";
	    $glist{$cog}{'glist'}.=$v;
	  }
      }

    print STDERR "\nIdentify Identical Groups";
    $ncog=$glist{'NCOG'};
    for ($a=0; $a<$ncog; $a++)
      {
	
	$cog1=$glist{'LCOG'}{$a};
	if ($glist{$cog1}{'used'}){next;}
	
	$string1=$glist{$cog1}{'glist'};
	print STDERR "[$cog1]";
	for ($b=0; $b<$ncog; $b++)
	  {
	    
	    $cog2=$glist {'LCOG'}{$b};
	    $string2=$glist{$cog2}{'glist'};
	    if ( $string1 eq $string2)
	      {
		$glist{$cog1}{'cog_list'}{$glist{$cog1}{'ncog'}++}=$cog2;
		$glist{$cog2}{'used'}=1;
	      }
	  }
      }
    print STDERR "\nPrint Results";
    for ( $a=0; $a<$ncog; $a++)
      {
	
	$cog=$glist{'LCOG'}{$a};
	$string=$glist{$cog}{'glist'};
	$nsub_cog=$glist{$cog}{'ncog'};
	if ( $ncog)
	  {
	    print "\nCOGLIST $nsub_cog ";
	    for ($b=0; $b< $nsub_cog; $b++)
	      {
		print "$glist{$cog}{'cog_list'}{$b} ";
	      }
	    print "\nGENOMELIST $nsub_cog $cog $string";
	  }
      }
    return %glist;
    }
	
sub cog2slist
  {
    my $whog=shift (@_);
    my $cog=shift (@_);
    my %hash_cog;
    my @list;

 
    print "TEST: $cog";
    open (F, $whog);
    while (<F>)
      {
	
	if ( $_=~/$cog/ )
	  {
	    print $_;
	    while (<F> )
	      {
		if ( $_=~/___/)
		  {
		    close(F);
		  }
		else
		  {
		    $line=$_;
		    $line=~s/:/ /g;
		    
		    @list=($line=~/(\S+)/g);
		    $hash_cog {$list[1]}{'genome'}=$list[0];
		    $hash_cog {$list[1]}{'north'}=$#list;
		    $hash_cog {$list[1]}{'seq'}=$list[1];
		    $hash_cog {$cog}{'nseq'}++;
		    $hash_cog {'genomelist'}.="$list[0]#";
		  }
	      }

	  }

      }

    close (F);
    $hash_cog {'cog'}=$cog;
    print "TEST: $cog";
    print "$hash_cog{'cog'} $hash_cog{'genomelist'}\n";
    return %hash_cog;
    }
sub slist2fasta
  {
    my $file=shift (@_);
    my %hash_cog=@_;
    my $result;

    $/=">";
    open (F, $file);
    while (<F>)
      {
	$line=$_;
	$line=~s/\>//g;
	($line=~/(\b\S+\b)/);
	$seq=$1;
	
	if ( $hash_cog{$seq}{'north'})
	  {
	    $line=~s/$seq/\>$hash_cog{$seq}{'genome'} $seq/;
	    $result.=$line;
	  }
      }
    close (F);
    $/="\n";
    return $result;
    }



  

