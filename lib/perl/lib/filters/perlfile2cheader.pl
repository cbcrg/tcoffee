#!/usr/bin/env perl

foreach $file_name (@ARGV)
  {
    push @file_name_list, $file_name;
    open (F, $file_name);
    
    $file_txt="";
    while (<F>)
      {
	if (!/^\#/){$file_txt.=$_;}
      }
    $file_txt=~s/\\/\\\\/g;
    $file_txt=~s/\r//g;
    $file_txt=~s/\n/\\n/g;
    $file_txt=~s/\"/\\"/g;
    $file_txt=~s/\'/\'/g;
    $file_txt=~s/\`/\`/g;
    
    close (fp);
    push @file_txt_list, $file_txt;
  }

$outfile="";
push @file_name_list, "EndList";
$outfile.="char *PerlScriptName[]={";
foreach $file_name (@file_name_list)
  {
    $outfile.="\"$file_name\",";
  }
chop $outfile;$outfile.="};\n";

$outfile.="char *PerlScriptFile[]={";
foreach $file_txt (@file_txt_list)
  {
    $outfile.="\"$file_txt\",";
  }
chop $outfile;$outfile.="};\n";
@char_list=($outfile=~/(.)/g);
foreach $c (@char_list)
  {
    print "$c";
    $n++;
    if ( $n==50){print "\\";print "\n";$n=0;}
  }
print "\n";





