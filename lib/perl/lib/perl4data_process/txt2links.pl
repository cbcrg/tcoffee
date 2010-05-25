#!/usr/bin/env perl


while (<>)
  {
    @links=/([\w-]+\.[\w-]+\.[\w-]+[\w\/\.-]+)/g;
    foreach $word (@links)
      {if ( $hash{$word}==0)
	 {$hash{$word}=1;
	 
	  @l=($word=~/(.)/g);
	  if ($l[$#l] eq '.'){chop $word;}
	  print "\nhttp://$word";
	}
     }
   @links=/http\:\/\/([\w-]+\.[\w-]+\.[\w-]+[\w\/\.-]+)/g;
    foreach $word (@links)
      {if ( $hash{$word}==0)
	 {$hash{$word}=1;
	 
	  @l=($word=~/(.)/g);
	  if ($l[$#l] eq '.'){chop $word;}
	  print "\nhttp://$word";
	}
     } 
  }
