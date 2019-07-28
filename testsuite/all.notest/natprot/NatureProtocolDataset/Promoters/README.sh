#!/bin/bash          
# Go to the repertory "Promoters"
# change executable permissions: chmod +x README.sh
# execute the script: ./README.sh
# 
# OR cut and paste these commands in your shell
#
# Correspondence between ENSEMBL codes and figure names:
#
#  ENSG00000177150 --------- HoSa
#  ENSMUSG00000038121 ------ MuMu
#  ENSCAFG00000000173 ------ CaFa
#  ENSBTAG00000009141 ------ BoTa
#  ENSGALG00000013887 ------ GaGa

work_dir=$PWD
source ~/.bashrc

 # pro-coffee
echo -n "Producing Pro-Coffee alignment ... "
mkdir -p $work_dir/pro_coffee
cd $work_dir/pro_coffee
(time t_coffee -seq $work_dir/c18orf19.fa -mode procoffee) >t_coffee.log 2>t_coffee.err
echo Done

 # t_coffee
echo -n "Producing T-coffee alignment ... "
mkdir -p $work_dir/t_coffee
cd $work_dir/t_coffee
(time t_coffee -seq $work_dir/c18orf19.fa) >t_coffee.log 2>t_coffee.err
echo Done

# Alignment extraction 
echo -n "Extraction alignment portion ... "
mkdir -p $work_dir/align_portion
cd $work_dir/align_portion
(time t_coffee -other_pg seq_reformat -in $work_dir/pro_coffee/c18orf19.aln -action +extract_block 'ENSG00000177150' 1852 1983) >c18orf19_chipseq.aln 2>t_coffee.err
echo Done

 #Time stuff
echo -n "Time used ..."
cd $work_dir
for mode in $(echo "t_coffee,pro_coffee,align_portion" | tr ',' '\n');do
timeused=$(grep real $mode/t_coffee.err | cut -f2)
echo -e "$mode\t$timeused"
done;

