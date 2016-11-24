#!/bin/bash
# Go to the repertory "RNAs"
# change executable permissions: chmod +x README.sh
# execute the script: ./README.sh
# 
# OR cut and paste these commands in your shell


work_dir=$PWD
source ~/.bashrc

 #produce R-Coffee alignment
echo -n "Producing R-Coffee alignment ... "
mkdir -p $work_dir/r_coffee
cd $work_dir/r_coffee
t_coffee -seq $work_dir/RNA_set.fa -mode rcoffee -outfile RNA_rcoffee.aln  >t_coffee.log 2>t_coffee.err
echo done

 #adding structure
echo -n "Predicting structure and projecting it ... "
mkdir -p $work_dir/pred_str_r_coffee
cd $work_dir/pred_str_r_coffee
t_coffee -other_pg seq_reformat -in $work_dir/r_coffee/RNA_rcoffee.aln -action +add_alifold -output stockholm_aln -out RNA_rcoffee_predicted.stk >t_coffee.log 2>t_coffee.err
echo done

echo -n "Projecting known structure ... "
mkdir -p $work_dir/add_str_r_coffee
cd $work_dir/add_str_r_coffee
t_coffee -other_pg seq_reformat -in $work_dir/r_coffee/RNA_rcoffee.aln -in2 $work_dir/1F6X_A_real.struc -input2 alifold -action +add_alifold -output stockholm_aln -out RNA_rcoffee_real.stk >t_coffee.log 2>t_coffee.err
echo done

 #Analysis
echo -n "Producing and analizing T-coffee alignment ... "
mkdir -p $work_dir/analysis_r_coffee
cd $work_dir/analysis_r_coffee
t_coffee -seq $work_dir/RNA_set.fa -outfile RNA_tcoffee.aln >t_coffee.log 2>t_coffee.err 
t_coffee -other_pg seq_reformat -in RNA_tcoffee.aln -action +alifold2analyze color_html  >RNA_tcoffee.comp.html 2>t_coffee.err
echo done

echo -n "Analizing R-coffee alignment ... "
t_coffee -other_pg seq_reformat -in $work_dir/r_coffee/RNA_rcoffee.aln -action +alifold2analyze color_html > RNA_rcoffee.comp.html 2>t_coffee.err
echo done

echo -n "Summarizing results in a table ..."
t_coffee -other_pg seq_reformat -in $work_dir/analysis_r_coffee/RNA_tcoffee.aln -action +alifold2analyze >results_t_coffee.log 2>results_t_coffee.err
t_coffee -other_pg seq_reformat -in $work_dir/r_coffee/RNA_rcoffee.aln -action +alifold2analyze >results_r_coffee.log 2>results_r_coffee.err
echo done

 #Printing table with results
echo -e "t_coffee\tr_coffee" 
Nseq=$(grep "Nseq" results_t_coffee.log | awk '{print $3}')
Length=$(grep "Nseq" results_t_coffee.log | awk '{print $4, $5}')
fold=$(grep "Nseq" results_t_coffee.log | awk '{print $6, $7}')
neut=$(grep "Nseq" results_t_coffee.log | awk '{print $8, $9}')
wat=$(grep "Nseq" results_t_coffee.log | awk '{print $10, $11}')
CompWC=$(grep "Nseq" results_t_coffee.log | awk '{print $12, $13}')
Comp=$(grep "Nseq" results_t_coffee.log | awk '{print $14, $15}')
Inc=$(grep "Nseq" results_t_coffee.log | awk '{print $16, $17}')

r_Nseq=$(grep "Nseq" results_r_coffee.log | awk '{print $3}')
r_Length=$(grep "Nseq" results_r_coffee.log | awk '{print $4, $5}')
r_fold=$(grep "Nseq" results_r_coffee.log | awk '{print $6, $7}')
r_neut=$(grep "Nseq" results_r_coffee.log | awk '{print $8, $9}')
r_wat=$(grep "Nseq" results_r_coffee.log | awk '{print $10, $11}')
r_CompWC=$(grep "Nseq" results_r_coffee.log | awk '{print $12, $13}')
r_Comp=$(grep "Nseq" results_r_coffee.log | awk '{print $14, $15}')
r_Inc=$(grep "Nseq" results_r_coffee.log | awk '{print $16, $17}')

echo -e "$Nseq\t\t$r_Nseq"
echo -e "$Length\t$r_Length"
echo -e "$fold\t\t$r_fold"
echo -e "$neut\t$r_neut"
echo -e "$wat\t$r_wat"
echo -e "$CompWC\t$r_CompWC"
echo -e "$Comp\t\t$r_Comp"
echo -e "$Inc\t$r_Inc"
