#!/bin/bash          
# Go to the repertory "Proteins"
# change executable permissions: chmod +x README.sh
# execute the script: ./README.sh
# 
# OR cut and paste these commands in your shell
#
# Two options are not shown in the protocol:
# -cache_update --> each run of t_coffee perform with 
#                   this option does not use files  
#                   stored in t_coffee cache
# -io_format    --> for aligment comparison, format of 
#                   input/output alignment

work_dir=$PWD
source ~/.bashrc

 # t_coffee
echo -n "Producing T-coffee alignment ... "
mkdir -p $work_dir/t_coffee
cd $work_dir/t_coffee
(time t_coffee -seq $work_dir/sh3.fasta -cache update) >t_coffee.log 2>t_coffee.err
echo Done

 # fmcoffee
echo -n "Producing Fast M-coffee alignment ... "
mkdir -p $work_dir/fmcoffee
cd $work_dir/fmcoffee
(time t_coffee -seq $work_dir/sh3.fasta -mode fmcoffee -cache update) >t_coffee.log 2>t_coffee.err
echo Done

 # fmcoffee with different methods
echo -n "Producing Fast M-coffee alignment with different methods... "
mkdir -p $work_dir/fmcoffee_dif_methods
cd $work_dir/fmcoffee_dif_methods
(time t_coffee -seq $work_dir/sh3.fasta -method muscle_msa probcons_msa clustalw_msa -cache update) >t_coffee.log 2>t_coffee.err
echo Done

 # psicoffee
echo -n "Producing PSI-coffee alignment ... "
mkdir -p $work_dir/psicoffee
cd $work_dir/psicoffee
(time  t_coffee -seq $work_dir/sh3.fasta -mode psicoffee -cache update) >t_coffee.log 2>t_coffee.err
echo Done

 # psicoffee with cache
echo -n "Re-Producing PSI-coffee alignment using the cache ... "
mkdir -p $work_dir/psicoffee_with_cache
cd $work_dir/psicoffee_with_cache
(time t_coffee -seq $work_dir/sh3.fasta -mode psicoffee) >t_coffee.log 2>t_coffee.err
echo Done

 # expresso 
echo -n "Producing Expresso alignment ... "
mkdir -p $work_dir/expresso
cd $work_dir/expresso
(time t_coffee -seq $work_dir/sh3.fasta -mode expresso -pdb_type dn -cache update) >t_coffee.log 2>t_coffee.err
echo Done

 # 3D-Coffee 
echo -n "Producing 3D-coffee alignment ... "
mkdir -p $work_dir/3D_coffee
cd $work_dir/3D_coffee
cp $CACHE_4_TCOFFEE/*.pdb .
(time t_coffee -seq $work_dir/sh3.fasta  -method sap_pair -template_file $work_dir/sh3.template_file -cache update) >t_coffee.log 2>t_coffee.err
rm *.pdb
echo Done

 # alignment comparisons
echo "Comparing all T-coffee alignments to the reference one "
cd $work_dir
for mode in $(echo "t_coffee,fmcoffee,fmcoffee_dif_methods,psicoffee,psicoffee_with_cache,expresso,3D_coffee" | tr ',' '\n');do
result_col=$(t_coffee -other_pg aln_compare -al1 sh3_ref.aln -al2 $mode/sh3.aln -io_format t -compare_mode column | sed s,"  "*," ",g | cut -d ' ' -f4)
timeused=$(grep real $mode/t_coffee.err | cut -f2)
echo -e "$mode\t$result_col\t$timeused"
done;

 # alignment T-RMSD 
echo -n "Producing Expresso structural alignment for T-RMSD ... "
mkdir -p $work_dir/T_RMSD
cd $work_dir/T_RMSD
###Getting PDBs
tar xzf $work_dir/PDBs.tar.gz
mv $work_dir/T_RMSD/PDBs/* . 
(time t_coffee -seq $work_dir/crd.fasta -template_file $work_dir/crd.template_file -method sap_pair TMalign_pair -cache update) >t_coffee.log 2>t_coffee.err
echo Done

 # Structural tree by T-RMSD
echo -n "Producing structural tree using T-RMSD ... "
mkdir -p $work_dir/Tree_T_RMSD
cd $work_dir/Tree_T_RMSD 
mv $work_dir/T_RMSD/*.pdb  $work_dir/Tree_T_RMSD/  
(time t_coffee -other_pg trmsd -aln $work_dir/T_RMSD/crd.aln -template_file $work_dir/crd.template_file -cache update)  >t_coffee.log 2>t_coffee.err
rm *.pdb
echo Done

