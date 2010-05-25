#! /usr/bin/perl
$version="5.09";

if ($#ARGV==-1)
  {
    print "coglist2cogfile -species <file> -olist <file> -folist <igs|ncbi>-e <extension>\n";
    print "Turns a list of COGs and Proteomes into a collection of files containing the corresponding sequences, with each gene named according to the genome\n\n";
    print "-e: extension\n";
    print "-species: list of files containing ORFs, can be in the form *.orf\n";
    print "\tgenome1.orf\n";
    print "\tgenome2.orf\n";
    print "\t...\n";
    
    print "-olist  : csv file containint the COGs (NCBI format)\n";
    print "\trig;genome1;genome2;genome3\n";
    print "\tRIG01;ORF1;ORF1;ORF36;ORF4\n";
    print "\t...\n";
    print "generates RIGxx.$extension";
    exit (EXIT_FAILURE);
  }

$command= join (" ", @ARGV);
$command=&protect_dash ($command);


if ( $command=~/\-species([^-]+)/)
  {
    $f=$1;
    $f=unprotect_dash($f);
    print "\nword: $f\n$command\n" ; 
    @species=split (/\s+/, $f);
    print "@species";
   
#    open (F, "$f");
#    while ( <F>)
#      {
#	print "$_";
#	$s.=$_;
#      }
    
#    @species=split (/[,\s]+/, $s);
#    close (F);
  }

$extension="cog_list";
if ( $command=~/\-e([^-]+)/)
  {
    $extension=$1;
    $extension=~s/\s//;
  }

if ( $command=~/\-olist([^-]+)/)
  {
    $olist=$1;
  }

if ( $command=~/\-folist([^-]+)/)
  {
    $folist=$1;
    $folist=~s/\s//g;
  }



if(($command=~/\-o (\S+)/))
  {
    $out=$1
  }
open (O, ">$out");
#1 Read ORF in species files
foreach $f (@species)
  {
    print "$f";
    $f=~/([^-.]+)/;
    $taxid=$1;
    lc($taxid);
    $n=0;
    open (F, "$f");
    $seq{$taxid}{list}=1;
    while (<F>)
      {
	if ( />/)
	  {
	    $x=$_;
	    $x=~/\>(\S+)/;
	    $name=$1;
	    lc($name);
	    $n++;
	  }
	else
	  {
	    $s=$_;
	    $s=~s/\s//g;
	    $seq{$taxid}{$name}.=$s;
	    $seq{$name}{taxid}=$taxid;
	  }
      }
    print "Read $n Sequences in File $f\n";
    close (F);
  }

#Read the COG List

print "FOLIST |$folist|";


if ( $folist eq "igs")
  {
    open (F, "$olist");
    while ( <F>)
      {
	$line=$_;
	$line=~s/[\n\r]+//;
	
	if (!@taxidlist)
	  {
	    lc ($line);
	    @taxidlist=split (/;/, $line);
	    shift (@taxidlist);
	  }
	else
	  {
	    @l=split (/;/, $line);
	    $cog=shift (@l);
	    foreach $s(@taxidlist)
	      {
		$clist{$cog}{$s}=shift (@l);
	      }
	  }
      }
    close (F);
  }
elsif ($folist eq "ncbi")
  {
    print "OPEN $olist";
    open (F, "$olist");
    while (<F>)
      {
	$line=$_;
	$line=~s/[\n\r]+//;
		
	if (($line=~/^(\S+).*(rCOG.*)/))
	  {
	    $name=$1;$cog=$2;
	    lc($name);
	    $taxid=$seq{$name}{taxid};
	    if (!$taxid)
	      {
		print "ERROR: $name could not be found [FATAL]\n"; die;
	      }
	    else
	      {
		$h_taxidlist{$taxid}=1;
		$clist{$cog}{$taxid}=$name;
		$clist{$cog}{$taxid}{'n'}++;
	      }
	  }
      }
    close (F);
  }
@taxidlist=keys(%h_taxidlist);
$tot=0;
$n_species=$#taxidlist+1;
foreach $c (keys (%clist))
  {
    foreach $s (@taxidlist)
      {
	if ($clist{$c}{$s}{'n'}==1)
	  {
	    $cog_stat{$c}{ns}++;
	  }
      }
    print "$c: $cog_stat{$c}{ns}\n";
    if ($cog_stat{$c}{ns}==$n_species)
      {
	$tot++;
      }
  }


foreach $c (keys (%clist))
  {
    
    $n=0;
    if ($cog_stat{$c}{ns}!=$n_species){next;}
    $seq=$clist{$c}{full_cog}="$c.$extension";
    print "Create $seq\n"; 
    open (F, ">$seq");
    foreach $s (@taxidlist)
      {

	if ($seq{$s}{$clist{$c}{$s}})
	  {
	    print F ">$s\n$seq{$s}{$clist{$c}{$s}}\n";
	  }


      }
    close (F);
  }
print "\n$tot COGs processed containing one member only in $n_species species\n";


sub protect_dash 
  {
    my $s=@_[0];
    $s=~s/(\S)\-(\S)/\1\@@@\2/g;
    return $s;
  }
sub unprotect_dash 
  {
    my $s=@_[0];
    $s=~s/(\S)@@@(\S)/\1\-\2/g;
    return $s;
  }
