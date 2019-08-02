#!/bin/bash          
# Go to the repertory "Post-processing"
# change executable permissions: chmod +x README.sh
# execute the script: ./README.sh
# 
# OR cut and paste these commands in your shell


source ~/.bashrc

 #MSA to fasta
echo -n "Producing sequence file from MSA ... "
t_coffee -other_pg seq_reformat -in sh3.aln -output fasta_aln > sh3_fasta.aln
echo done

 #Extracting blocks
echo -n "Extracting block 10 to 20 (respect consensus) ... "
t_coffee -other_pg seq_reformat -in sh3.aln -action +extract_block cons 10 20
echo done
echo -n "Extracting block 20 to 40 (respect seq 1PHT) ... "
t_coffee -other_pg seq_reformat -in sh3.aln -action +extract_block '1PHT' 20 40
echo done

 #Coloring alignments
echo -n "Coloring aligment ... "
t_coffee -other_pg seq_reformat -in sh3.aln -action +color_residue 1PHT 10 1 -output color_html > color.html
echo done

echo -n "Coloring aligment from file ... "
t_coffee -other_pg seq_reformat -in sh3.aln -action +color_residue color.txt -output color_html > file2color.html
echo done

echo -n "Coloring by conservation ... "
t_coffee -other_pg seq_reformat -in sh3.aln -in3 sh3.aln -action +3evaluate pam250mt -output color_html > cons2color.html
echo done

echo -n "Coloring by boxshade ... "
t_coffee -other_pg seq_reformat -in sh3.aln -in3 sh3.aln -action +3evaluate boxshade -output color_html > boxshade2color.html
echo done

 #Comparison
echo -n "Estimating diversity ... "
t_coffee -other_pg seq_reformat -in sh3.aln -output sim
echo done

echo -n "Comparing alignments ... "
t_coffee -other_pg aln_compare -al1 sh3_ref.aln -al2 sh3.aln -compare_mode column
echo done

echo -n "Comparing trees ... "
t_coffee -other_pg seq_reformat -in crd.struc_tree100 -in2 crd.struc_tree50  -action +tree_cmp 
echo done


