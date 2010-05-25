#!/usr/bin/env perl
  
#Defaults
  $np=0;
  $n_para=$#ARGV;
  $model=1;


  foreach ($np=0; $np<=$n_para; $np++)
    {        
    $v=$ARGV[$np];
    
    if ( $v eq "-pdb")
       {
       $pdb_name= $ARGV[++$np];
       }
    elsif (  $v eq "-pep")
       {
       $pep_file= $ARGV[++$np];
       }
    }

$pdb_name = lc($pdb_name);

$pdb_name=`export NO_REMOTE_PDB_DIR=1;export PDB_DIR=/Sequences/Pdb/PdbMirror;/usr/local/bin/extract_from_pdb -get_fugue_name $pdb_name`;
chop $pdb_name;

if ( $pdb_name eq non_valid_pdb_name)
  {
    exit (EXIT_FAILURE)
  }
$pdb_file = "$pdb_name.pdb";
$x = `export NO_REMOTE_PDB_DIR=1;export PDB_DIR=/Sequences/Pdb/PdbMirror;extract_from_pdb -netfile $pdb_name -mode raw > $pdb_file`;
if (not -e $pdb_file) {
  $pdb_file = $pdb_name;
  $pdb_file =~ s/.$/.pdb/;
  if (not -e $pdb_file) {
    exit (EXIT_FAILURE)
  }
}

# open (F, $tmp);
open (F, "/home/igs/Tools/T-COFFEE_web/Struct/fugue_align.sh $pep_file $pdb_file|");
while (<F>)
  {
    s/^>P1;/>/;
    s/\*//;
    next if ( /structureX:/ );
    next if ( /sequence:/ );
    next if ( /^ *$/ );
    print $_;
  }
close (F);
unlink $tmp;



	  
    
