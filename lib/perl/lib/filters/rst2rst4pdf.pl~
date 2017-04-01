#!/usr/bin/env perl
use strict;

#This script reads an Interleaved MSA and Outputs a Fasta_msa
#if two sequences have the same name, the second one is renamed name_1 and so on
my %Tag;
my @TagL;
my @txt;
my $nl;
my $ptype;

setTag ("<h1", "</h1", "h1");
setTag ("<h2", "</h2", "h2");
setTag ("<h3", "</h3", "h3");
setTag ("<h4", "</h4", "h4");
setTag ("<h5", "</h5", "text");
setTag ("<table class=", "</table", "table");
setTag ("<p class=ChapterTitle", "</p>", "h1");
setTag ("<p class=TitleCover", "</p>", "title");
setTag ("<p class=MsoTitle", "</p>", "title");

setTag ("<p class=MyNorma", "</p>", "text", '"');
setTag ("<p class=Mylist", "</p>", "text", '"');

setTag ("<p class=MsoBodyText", "</p>", "text");
setTag ("<p class=noteboxbody", "</p>", "file");
setTag ("<p class=Note", "</p>", "note");


setTag ("<p class=MyExample", "</p>", "example");
setTag ("<p class=example", "</p>", "example");
setTag ("<p class=Formula", "</p>", "example");

setTag ("<p class=MsoNormal", "</p>", "example");
setTag ("<p class=examplefile", "</p>", "file");
setTag ("<p class=Warning", "</p>", "warning");
setTag ("<p class=Stylenotebox", "</p>", "box");
setTag ("<p class=flag>", "</p>", "flag");
setTag ("<p class=flagusage", "</p>", "flagusage");
setTag ("<p class=instruction", "</p>", "flaginstruction");
setTag ("<p class=flagdescription", "</p>", "flagdescription");
setTag ("<p class=flagheading", "</p>", "flagdescription");

setTag ("<p class=FAQ", "</p>", "faq");
setTag ("<p class=output", "</p>", "text");

@TagL=keys(%Tag);
my $convert_utf=0;
my $check_class=0;

my $infile=$ARGV[0];

if ($convert_utf)
  {
    $infile="tmp.html";
    if (-e $infile){unlink ($infile);}
    system (" cat $ARGV[0]| iconv -f windows-1252 -t utf-8 > $infile");
  }

if (!-e $infile){die;}

if ($check_class)
  {my %clist;
   my @TL=keys(%Tag);
   open (IN, $infile);
   while (<>)
     {
       my $l=$_;
       if ($l=~/\<p class=(\w+)/)
	 {
	   my $c=$1;
	   if (!$clist{$c})
	     {
	       my $no=1;
	       my $current="<p class=$c>";
	       
	       $clist{$c}=1;
	       foreach my $t (@TL)
		 {
		   if ($current=~/$t/){$no=0;last;}
		 }
	       if ($no){print " $c NO\n";}
	     }
	 }
     }
 }




open (IN, $infile);

my $in=0;
my $current="";
my $EndTag;
my $StartTag;
my $n;
while ( <IN>)
  {
    my $l=$_;
    my $parse=0;

    if (!$in)
      {foreach my $ST (@TagL)
	 {
	   if ($l =~/$ST/)
	     {
	       $StartTag=$ST;
	       $EndTag=$Tag{$ST}{end};
	       $parse=1;
	       last;
	     }
	 }
     }
    
    if ($parse)
      {
	if (!($l=~/$EndTag/))
	  {
	    $current.=$l;
	    $in=1;
	  }
	else 
	  {
	    htmlstring2data($l,$StartTag);
	    $n=0;
	  }
      }
    elsif ($in && $l=~/$EndTag/)
      {
	$current .=$l;
	htmlstring2data($current, $StartTag);
	$StartTag="";
	$EndTag="";
	
	$current="";
	$in=0;
	$n=0;
      }
    elsif ($in)
      {
	$current .=$l;
      }
  }

for (my $n=0; $n<$nl; $n++)
  {
    txt2rst($txt[$n][0],$txt[$n][1]);
  }

sub txt2rst 
  {
    my ($st, $tag)=@_;
    my $type=$Tag{$tag}{rst};

    my $len=length $st;


    if ($ptype eq "file" && $type ne "file")
      {
	print "\n\n";
      }
    elsif ($ptype eq "box" && $type ne "box") 
      {
	print "\n\n";
      }
    elsif ($ptype eq "example" && $type ne "example")
      {
	print "\n\n";
      }

    if    ($type eq "title")
      {
 	my $line=generate_string($len, '#');
	print "$line\n$st\n$line\n";
      }
    elsif ($type eq "table")
      {
	
	$st=html_table2rst_table($st);
      }
    elsif ($type eq "h1")
      {
	my $line=generate_string($len, '*');
	print "$line\n$st\n$line\n";
      }
    elsif ($type eq "h2")
      {
	my $line=generate_string($len, '=');
	print "$st\n$line\n";
      }
    elsif ($type eq "h3")
      {
	my $line=generate_string($len, '-');
	print "$st\n$line\n";
      }
    elsif ($type eq "h4")
      {
	my $line=generate_string($len, '^');
	print "$st\n$line\n";
      }
    elsif ($type eq "text")
      {
	print "$st\n\n\n";
      }
    elsif ($type eq "example")
      {
	if ($ptype ne "example"){print "::\n\n";}
	print "  $st\n\n";
      }
    elsif ($type eq "file")
      {
	if ($ptype ne "file"){print "::\n\n";}
	print "  $st\n";
      }
    
    elsif ($type eq "box")
      {
	if ($ptype ne "box"){print "::\n\n";}
	#$st=string2split ($st, 50);
	print "  $st\n";
      }
    elsif ($type eq "warning")
      {
	print ".. warning:: $st\n\n";
      }
    elsif ($type eq "note")
      {
	print ".. note:: $st\n\n";
      }
    elsif ($type eq "flag")
      {
	my $line=generate_string($len, '^');
	print "$st\n$line\n";
      }
    elsif ($type eq "flagusage")
      {
	print "  **$st**\n\n";
      }
    elsif ($type eq "flaginstruction")
      {
	print "   *$st*\n\n";
      }
    elsif ($type eq "flagdescription")
      {
	print "   *$st*\n\n";
      }
    elsif ($type eq "faq")
      {
	print "   **$st**\n\n";
      }
    $ptype=$type;
  }


sub generate_string
  {
    my ($len, $char)=@_;
    my $return;
    for (my $a=0; $a<$len; $a++)
      {
	$return.=$char;
      }
    return $return;
  }
    
sub htmlstring2data 
  {
    my ($string, $tag)=@_;
    my $stringi=$string;
    
    
    if (!($string =~/^<table/)){$string=~s/<[^><]*>//g;}
    $string=~s/\n/ /g;
    $string=~s/&#8211;/\-/g;
    $string=~s/&lt;/\</g;
    $string=~s/&gt;/\>/g;
    $string=~s/&quot;/\'/g;
    $string=~s/&nbsp;/ /g;
    $string=~s/^\s+//g;
    $string=~s/\s+$//g;
    
    #$string=~s/\*/\\\*/g;
    
    my $c='\222';$string=~s/$c/\'/g;
    my $c='\223';$string=~s/$c/\'/g;
    my $c='\224';$string=~s/$c/\'/g;
    my $c='\205';$string=~s/$c/\.\.\./g;
    

    if ($string =~/\S/ || $string =~/Word did not find any entries/)
      {
	$txt[$nl][0]=$string;
	$txt[$nl][1]=$tag;
	$nl++;
      }
  }
sub string2split
  {
    my ($string, $size)=@_;
    my $s1=$string;
    $string="";
    
    my @cl=split (//, $s1);
    my $cc=0;
    my $l=length $s1;
    for (my $c=0; $c<$l; $c++)
      {
	my $ch=$cl[$c];
	
	if ($cc>=60 && $ch eq " "){$string .="\n";$cc=0;}
	else {$string.=$ch;}
	$cc++;
      }
    return $string;
  }
sub setTag
    {
      my ($start,$end,$rst)=@_;
      $Tag{$start}{start}=$start;
      $Tag{$start}{end}=$end;
      $Tag{$start}{rst}=$rst;
      
      return;
    }
    
sub html_table2rst_table 
  {
    my ($string)=(@_);
    
    
    $string=~s/<table[^>]*>/<table>/g;
    $string=~s/<\/table[^>]*>/<\/table>/g;
    $string=~s/<tr[^>]*>/<tr>/g;
    $string=~s/<\/tr[^>]*>/<\/tr>/g;
    $string=~s/<td[^>]*>/<td>/g;
    $string=~s/<\/td[^>]*>/<\/td>/g;
    $string=~s/<p[^>]*>//g;
    $string=~s/\n/ /g;
    
    $string=~m/<table>(.*)<\/table>/;
    my @row=(split (/<\/tr>/, $1));

    my ($cr,$nr, $cc, $nc,@table);
    foreach my $r (@row)
      {
	$cr++;
	$cc=0;
        $r=~s/<tr>//g;
	my @cells=split (/<\/td>/, $r);
	foreach my $c (@cells)
	  {
	    $c=~s/<[^><]*>//g;
	    $table[$cr][$cc++]=$c;
	  }
	if ($cc > $nc){$nc=$cc;}
      }
    $nr=$cr;
    
    #trim entries blank extremities
    my (@rl, @cl);
    for (my $c=0; $c<$nc; $c++)
      {
	for (my $r=0; $r<$nr; $r++)
	  {
	    my $s=$table[$r][$c];
	    $s=~s/^\s+//g;
	    $s=~s/\s+$//g;
	    $table[$r][$c]=$s;
	    $rl[$r]+=length($s);
	    $cl[$c]+=length($s);
	  }
      }
    #get rid pf empty rows
    my (@table2, $cr2,$nr2);
    #trim empty rows
    for (my $r=0; $r<$nr; $r++)
      {
	if ($rl[$r])
	  {
	    for (my $c=0; $c<$nc; $c++)
	      {
		$table2[$cr2][$c]=$table[$r][$c];
	      }
	    $cr2++;
	  }
      }
    my $nr2=$cr2;

    #get rid of empty columns
    my (@table3, $cc2);
    for ( my $c=0; $c<$nc; $c++)
      {
	if ($cl[$c])
	  {
	    for (my $r=0; $r<$nr2; $r++)
	      {
		$table3[$r][$cc2]=$table2[$r][$c];
	      }
	    $cc2++;
	  }
      }
    
    my $nc=$cc2;
    my $nr=$nr2;
    my @table=@table3;
    

    #measure each column max length
    my @cml;
    for (my $r=0; $r<$nr; $r++)
      {
	for (my $c=0; $c<$nc; $c++)
	  {
	    my $l=length ($table[$r][$c]);
	    if ($l>$cml[$c]){$cml[$c]=$l;}
	  }
      }
        
    

    for (my $r=0; $r<$nr; $r++)
      {
	if ($r==0)
	  {
	    for (my $c=0; $c<$nc; $c++)
	      {
		my $pad=generate_string($cml[$c], '=');
		print "$pad ";
	      }
	    print "\n";
	  }
	for (my $c=0; $c<$nc; $c++)
	      {
		my $l=length ($table[$r][$c]);
		print "$table[$r][$c]";
		my $pad=generate_string($cml[$c]-$l, ' ');
		print "$pad ";
	      }
	print "\n";
	if ($r==0 || $r==$nr-1)
	  {
	    for (my $c=0; $c<$nc; $c++)
	      {
		my $pad=generate_string($cml[$c], '=');
		print "$pad ";
	      }
	    print "\n";
	  }
      }
    print "\n";
  }
