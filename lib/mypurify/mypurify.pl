#!/usr/bin/env perl


while (<>)
  {
    $line=$_;
    print "$line\n";
#1 protect comments    
    @pstrings=($line=~/(\"[^"]*\")/g);
    $n=0;
    while ($x=($line=~s/\"[^"]*\"/ mypurify_string_$n /)){$n++;}
    print $line;

#2 Protect Structure dereferencing
    $line=~s/\-\>/mypurify_code0/g;

#2 attack []
    $line=~s/\]/\] /g;
    $line=~s/\]\s*\[/\]\[/g;
    $line=~s/\] /\]MP_EOF/g;
    $line=~s/([\w()]+\[)/MP_SOF$1/g;
    print $line;


#unprotect all
    $line=~s/mypurify_code0/->/g;
    foreach ($a=0; $a<$n; $a++)
      {
	$string=" mypurify_string_$a ";
	$line=~s/$string/$pstrings[$a]/;
      }
    print "\n$line";
    
  }
