#!/usr/bin/env perl
use Env;
use Cwd;
@suffix=("tmp", "temp", "cache", "t_coffee", "core", "tcoffee");

if ($#ARGV==-1)
  {
    print "clean_cache.pl -file <file to add in -dir> -dir=<dir> -size=<value in Mb>\n0: unlimited -1 always.\nWill only clean directories matching:[";
    foreach $k(@suffix){print "*$k* ";}
    print "]\n";
    exit (EXIT_FAILURE);
  }

$cl=join (" ",@ARGV);
if (($cl=~/\-no_action/))
  {
    exit (EXIT_SUCCESS);
  }

if (($cl=~/\-debug/))
  {
    $DEBUG=1;
  }
else
  {
    $DEBUG=0;
  }

if (($cl=~/\-dir=(\S+)/))
  {
    $dir=$1;
  }
else
  {
    $dir="./";
  }

if ($cl=~/\-file=(\S+)/)
  {
    $file=$1;
  }
else
  {
    $file=0;
  }

if ($cl=~/\-size=(\S+)/)
  {
    $max_size=$1;
  }
else
  {
    $max_size=0;#unlimited
  }
if ($cl=~/\-force/)
  {
    $force=1;
  }
else
  {
    $force=0;
  }

if ($cl=~/\-age=(\S+)/)
  {
    $max_age=$1;
  }
else
  {
    $max_age=0;#unlimited
  }

$max_size*=1000000;
if ( ! -d $dir)
  {
    print STDERR "\nCannot process $dir: does not exist \n";
    exit (EXIT_FAILURE);
  }

if ( !($dir=~/^\//))
  {
    $base=cwd();
    $dir="$base/$dir";
  }

$proceed=0;
foreach $s (@suffix)
  {
    
    if (($dir=~/$s/)){$proceed=1;}
    $s=uc ($s);
    if (($dir=~/$s/)){$proceed=1;}
  }
if ( $proceed==0)
  {
    print STDERR "Clean_cache.pl can only clean directories whose absolute path name contains the following strings:";
    foreach $w (@suffix) {print STDERR "$w ";$w=lc($w); print STDERR "$w ";}
    print STDERR "\nCannot process $dir\n";
    exit (EXIT_FAILURE);
  }

$name_file="$dir/name_file.txt";
$size_file="$dir/size_file.txt";
if ( $force){&create_ref_file ($dir,$name_file,$size_file);}
if ($file){&add_file ($dir, $name_file, $size_file, $file);}
&clean_dir ($dir, $name_file, $size_file, $max_size,$max_age);
exit (EXIT_SUCCESS);

sub clean_dir 
  {
    my ($dir, $name_file, $size_file, $max_size, $max_age)=@_;
    my ($tot_size, $size, $f, $s);

  
    $tot_size=&get_tot_size ($dir, $name_file, $size_file);

    if ( $tot_size<=$max_size){return ;}
    else {$max_size/=2;}
    
    #recreate the name file in case some temprary files have not been properly registered
    &create_ref_file ($dir, $name_file, $size_file, $max_age);
  
    $new_name_file=&vtmpnam();
    open (R, "$name_file");
    open (W, ">$new_name_file");
    while (<R>)
      {
	my $line=$_;
	
	($f, $s)=($line=~/(\S+) (\S+)/);
	if ( !($f=~/\S/)){next;}
	
	elsif ($max_size && $tot_size>=$max_size && !($f=~/name_file/))
	  {
	    remove ( "$dir/$f");
	    $tot_size-=$s;
	  }
	elsif ( $max_age && -M("$dir/$f")>=$max_age)
	  {
	    remove ( "$dir/$f");
	    $tot_size-=$s;
	  }
	else
	  {
	    print W "$f $s\n";
	  }
      }
    close (R);
    close (W);
    open (F, ">$size_file");
    print F "$tot_size";
    if ( -e $new_name_file){`mv $new_name_file $name_file`;}
    close (F);
  }
sub get_tot_size
  {
    my ($dir, $name_file, $size_file)=@_;
    my $size;
    
    if ( !-d $dir){return 0;}
    if ( !-e $name_file)
      {
	
	&create_ref_file ($dir, $name_file, $size_file);
      }
    open (F, "$size_file");
    $size=<F>;
    close (F);
    chomp ($size);
    return $size;
  }
sub size 
  {
    my $f=@_[0];

    if ( !-d $f){return -s($f);}
    else {return &dir2size($f);}
  }
sub dir2size
  {
    my $d=@_[0];
    my ($s, $f);
    
    if ( !-d $d) {return 0;}
    
    foreach $f (&dir2list ($d))
      {
	if ( -d $f){$s+=&dir2size ("$d/$f");}
	else {$s+= -s "$dir/$f";}
      }
    return $s;
  }

sub remove 
  {
    my $file=@_[0];
    my ($f);
    
    debug_print( "--- $file ---\n");
    if (($file eq ".") || ($file eq "..") || ($file=~/\*/)){return EXIT_FAILURE;}
    elsif ( !-d $file)
      {
	debug_print ("unlink $file\n");
	if (-e $file){unlink ($file);}
      }
    elsif ( -d $file)
      {
	debug_print ("++++++++ $file +++++++\n");
	foreach $f (&dir2list($file))
	  {
	    &remove ("$file/$f");
	  }
	debug_print ("rmdir $file\n");
	rmdir $file;
      }
    else
      {
	debug_print ("????????? $file ????????\n");
      }
    return EXIT_SUCCESS;
  }

sub dir2list
  {
    my $dir=@_[0];
    my (@list1, @list2,@list3, $l);

    opendir (DIR,$dir);
    @list1=readdir (DIR);
    closedir (DIR);
    
    foreach $l (@list1)
      {
	if ( $l ne "." && $l ne ".."){@list2=(@list2, $l);}
      }
    @list3 = sort { (-M "$dir/$list2[$b]") <=> (-M "$dir/$list2[$a]")} @list2;
    return @list3;
    
  }

sub debug_print
  {
    
    if ($DEBUG==1){print @_;}
    
  }
sub create_ref_file
  {
    my ($dir,$name_file,$size_file)=@_;
    my ($f, $s, $tot_size, @l);
    
    if ( !-d $dir){return;}
    
    @l=&dir2list ($dir);
    open (F, ">$name_file");
    foreach $f (@l)
      {
	$s=&size("$dir/$f");
	$tot_size+=$s;
	print F "$f $s\n";
      }
    &myecho ($tot_size, ">$size_file");
    close (F);
  }
sub add_file 
  {
    my ($dir,$name_file,$size_file,$file)=@_;
    my ($s, $tot_size);
    
    if ( !-d $dir)   {return;}
    if ( !-e "$dir/$file" ) {return;}
    if ( !-e $name_file){&create_ref_file ($dir,$name_file,$size_file);}
					    
    $s=&size("$dir/$file");
    open (F, ">>$name_file");
    print F "$file\n";
    close (F);

    $tot_size=&get_tot_size ($dir,$name_file,$size_file);
    $tot_size+=$s;
    &myecho ($tot_size, ">$size_file");
    
  }
	
sub myecho
  {
    my ($string, $file)=@_;
    open (ECHO, $file) || die;
    print ECHO "$string";
    close (ECHO);
  }
    
		
	
sub vtmpnam
  {
    my $tmp_file_name;
    $tmp_name_counter++;
    $tmp_file_name="tmp_file_for_clean_cache_pdb$$.$tmp_name_counter";
    $tmp_file_list[$ntmp_file++]=$tmp_file_name;
    if ( -e $tmp_file_name) {return &vtmpnam ();}
    else {return $tmp_file_name;}
  }

