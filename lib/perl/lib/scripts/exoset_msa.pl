#!/usr/bin/env perl
#usage exoset_msa.pl <sequences> <coded_sequence>
#sequence : original sequence
#coded_sequence: sequences where exon boundary residues have been repalced with boj


($ARGV[0]=~/^([^.]+)/);
$suffix=$1;
$exon_lib1="$suffix.tc_lib";



`t_coffee $ARGV[1] -in exon_pair -out_lib $suffix.1.tc_lib  -lib_only`;
`t_coffee -other_pg seq_reformat -in $suffix.1.tc_lib -action +weight_lib 1000 +swap_lib_header $ARGV[0] -output tc_lib > $suffix.tc_lib`;
`t_coffee $ARGV[0] -in L$suffix.tc_lib slow_pair lalign_id_pair`;


