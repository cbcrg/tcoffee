#!/usr/bin/env perl


foreach ($np=0; $np<=$#ARGV; $np++)
    {
    $value=$ARGV[$np];

    if ($value eq "-file")
      {
      $file= $ARGV[++$np];
      }
    elsif ($value eq "-type")
      {
        $type= $ARGV[++$np];
      }
    elsif ($value eq "-institute")
      {
        $institute= $ARGV[++$np];
      }
    elsif ($value eq "-author")
      {
        $author= $ARGV[++$np];
      }
    elsif ($value eq "-date")
      {
        $date= $ARGV[++$np];
      }
     elsif ($value eq "-program")
      {
        $program= $ARGV[++$np];
      }
    elsif ($value eq "-email")
      {
        $email= $ARGV[++$np];
      }
    else
      {
	print "$value is an unkown argument[FATAL]\n";
	exit (1);
      }
  }



open F, $file || die;
print $INSTITUTE;
if ( $type eq "c"){print "/*********************************COPYRIGHT NOTICE**********************************/\n";}
if ( $type eq "perl"){print "#################################COPYRIGHT NOTICE#################################/\n";}
if ( $type eq "txt"){print "----------------------------------COPYRIGHT NOTICE---------------------------------/\n";}


while (<F>)
  {
  s/\$INSTITUTE/$institute/g;
  s/\$AUTHOR/$author/g;
  s/\$DATE/$date/g;
  s/\$PROGRAM/$program/g;  
  s/\$EMAIL/$email/g;  
  if ( $type eq "txt"){print $_;}
  elsif ($type eq "c"){chop $_; print "\/*$_*\/\n";}
  elsif ($type eq "perl"){print "\#$_";}
}
close (F);
if ( $type eq "c"){print "/*********************************COPYRIGHT NOTICE**********************************/\n";}
if ( $type eq "perl"){print "#################################COPYRIGHT NOTICE#################################/\n";}
if ( $type eq "txt"){print "----------------------------------COPYRIGHT NOTICE---------------------------------/\n";}


