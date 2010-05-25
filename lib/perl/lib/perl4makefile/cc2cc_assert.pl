#!/usr/bin/env perl

if ( $#ARGV<0)
  {
    `ls -1 *.c > c_file_list`;
    
  }
else
  {
    open (F, ">c_file_list");
    foreach $f (@ARGV)
      {
	print F "$f\n";
      }
    close (F);
    }
open (F, "c_file_list");

while (<F>)
  {
    $file=$_;
    chomp ($file);
    
    `cp $file $file.original`;
    open (IN, "$file.original");
    while ( <IN>) {if ( /\/\*MBC\*\// && /assert/){$clean=1;}}
    close (IN);open (IN, "$file.original");
    
    open (OUT, ">$file");
    if ( $clean) {print (STDERR "Cleanning my_assert in $file\n");}
    else {print (STDERR "Adding my_assert tags in $file\n");}

    if ($clean!=1){print (OUT "#include <assert.h>\n");}
    while (<IN>)	
	{
	  if ( $clean==1 && /assert\.h/)
	    {;}
	  elsif ( /\/\*MBC\*\// && /assert/)
	    {
	      $line=$_;
	      $line=~s/assert.*\)//;
	      print OUT "$line\n";
	    }
	  elsif ( /\/\*MBC\*\//)
	    {

	      $line=protect_string ($_);
	      $instruction=line2instruction ($line);
	      $instruction=unprotect_string ($instruction);
	      $i=$_;chomp ($i);
	      $line="{$i$instruction}";
	      print OUT "$line\n";
	    }
	  else
	    {
	      print OUT "$_";
	    }
	}
 
    close (OUT);
    close (IN);
  }
close (F);

sub line2instruction
  {
  my $line=@_[0];
  my $new_line;
  my @list;
  my $instruction;
  my $full_intruction;
  

  @list=(split (//, $line));

  for ($a=0; $a<=$#list; $a++)
    {
      $c=$list[$a];
      if ( $c eq "["){$index++;$list[$a]="Obracket$index ";}
      elsif ( $c eq "]"){$list[$a]="Cbracket$index ";$index--;}
      if ($index>$max_index){$max_index=$index;}
    }
  
  
  foreach $c (@list)
    {$new_line.=$c;}


  for ($a=1; $a<=$max_index; $a++)
    {
      $b=$a;
      $token_address="#$a"."_";
      $new_line=~s/(\w+Obracket$a)/$token_address$1/g;
    }

  for ($a=1; $a<=$max_index; $a++)
    {
      $buffer_line=$new_line;
      $token_address="#$a"."_";
      $buffer_line=~s/$token_address/\@/g;
      
      $buffer_line=~s/Obracket$a/\[/g;
      $buffer_line=~s/Cbracket$a/\]/g;
      @result=($buffer_line=~/@([^@]*)\[([^]]+)\]/g);
 
	for ($b=0; $b<=$#result; $b+=2)
	{
	  $instruction.="assert (my_assert ($result[$b], $result[$b+1]));";
	}
    }
  return  $instruction;
}

sub protect_string
  {
    my $line=@_[0];


    @string_list=($line=~/(\"[^"]+\")/g);

    $n=0;
    while (($line=~s/(\"[^"]+\")/TOKEN_STRING$n /)){$n++;}

    $line=~s/\-\>/TOKEN_POINTER/g;
    $line=~s/#/TOKEN_HASCH/g;

    $line=~s/\s*\[\s*/\[/g;
    $line=~s/\s*\]\s*/\]/g;
    return $line;
  }

sub unprotect_string
  {
    my $line=@_[0];
    my $a;
    
    for ( $a=0; $a<$n_strings; $a++)
      {
	$line=~s/TOKEN_STRING$a /$string_list[$a] /;
      }
    $line=~s/Obracket\d+/\[/g;
    $line=~s/Cbracket\d+/\]/g;
    
    $line=~s/TOKEN_POINTER/\-\>/g;
    $line=~s/TOKEN_HASCH/#/g;
    $line=~s/#\d+_//g;
    return $line;
    }
