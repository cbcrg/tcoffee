#!/usr/bin/env perl


while (<>)
  {
    s/-/hyphen/g;
    s/\s/space/g;
    s/\W/\ /g;
    
    
    s/space/\ /g;
    s/hyphen/\-/g;
   
    
    
    tr/A-Z/a-z/;
    @list=/(\w+)/g;
    foreach $word (@list){if ( $hash{$word}==0){$hash{$word}=1;print "\n$word";}}
 

  }
