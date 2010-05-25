#!/usr/bin/env perl
#Version 1.01 (25/02/03)
$upper_threshold=9;
$lower_threshold=0;
$char='-';
$output="fasta_aln";
$keep=0;

foreach ($np=0; $np<=$#ARGV; $np++)
    {
    $value=$ARGV[$np];

    if ($value eq "-in")
      {
      $infile= $ARGV[++$np];
      }
    elsif ($value eq "-script")
      {
        $script= $ARGV[++$np];
      }
  }

open (F, $script);
while (<F>){$ref_command.=$_;}
chop ($ref_command);chop ($ref_command);chop ($ref_command);


$string=`seq_reformat -in $infile -output name`;
@name_list=($string=~/(\S+)/g);


foreach $name (@name_list)
  {
    $tmp_name=&vtmpnam();
    $hash{$name}{'File'}=$tmp_name;
    
    `seq_reformat -in $infile -action +keep_name +extract_seq $name > $tmp_name`;
  }

foreach $name1 ( @name_list)
  {
    foreach $name2 (@name_list)
      {
	$result=&vtmpnam();
	chomp ($name1);
	chomp ($name2);
	if ( $name1 ne $name2)
	  {
	    $command=$ref_command;
	    $param="$hash{$name1}{'File'} $hash{$name2}{'File'} ";
	    $command=~s/DATAFILE/$param/g;
	    $command.=" \>$result";

	    
	    system ($command);
	    $command="aln_compare -al1 $infile -al2 $result -io_format t -sep '' ";
	    $R=`$command`;chomp($R);chomp($R);
	    print "$R $name1 $name2\n";
	  }
      }
  }



# TMP NAME MANAGEMENT
sub vtmpnam
  {
    my $prefix=@_[0];
    my $tmp_file_name;

    $tmp_prefix=($prefix)?$prefix:"dpa_tmp_file_$$";
   
    $tmp_count++;
    $tmp_file_name="$tmp_prefix"."$tmp_count";
    $tl[$#tl+1]=$tmp_file_name;
    return $tmp_file_name;
  }



sub clean_tmp_file
  {

    my $list;
    my $file;
    if ($noclean){return;}

    $list=vtmpnam();
    `ls -1 | grep $tmp_prefix>$list`;
    
    open (F,$list);
    while ( <F>)
      {
	$file=$_;
	chop $file;
	if ( -e $file){unlink $file;}
      }
    close (F);
    unlink $list;
  }


sub exit_dpa
  {
  my $condition=@_[0];
  my $error_msg=@_[1];
  my $exit_value=@_[2];
  if ( $condition)
    {
      print "$error_msg\n";
      exit ($exit_value);
    }
  else
    {
      return;
    }
  
}


