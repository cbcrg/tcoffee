###################
T-Coffee Web Server 
###################

.. warning:: This chapter is currently under maintenance. Its aim is to describe what the T-Coffee web server is doing programmatically when you are using it; it provides technical details about the restrictions and commands in use.

In this chapter we describe briefly the available T-Coffee web servers and their usage (command lines and examples files); if you need more advanced options, you should use the T-Coffee package via command lines. The web server is a reduced version of the T-Coffee package, containing all T-Coffee modes for aligning protein, RNA or DNA sequences; it also contains the evaluation and downstream analysis tools (TCS, iRMSD/APDB, STRIKE, T-RMSD). Currently, all reformatting utilities are not available on the web server, however, you can choose some reformatting options related to the ouput format. Here we present briefly the webserver and the T-Coffee commands it contains.

*******
General
*******
Most of the following command lines used by the web server contains lots of different options called with flags; here is a summary of the options common to all comand lines:

  - **-in**: your input file
  - **-output**: specifies the output files you require
  - **-maxnseq**: maximum number of sequences you can run (variable depending on the mode)
  - **-maxlength**: maximum length of your sequences (variable depending on the mode)
  - **-case=upper**: all residues/nucleotids will be upper case
  - **-seqnos=off**: the size of the sequences is not indicated on the MSA
  - **-outorder=input**: orders the sequences in the final MSA as in the input dataset 
  - **-run_name**: name of the job on the cluster
  - **-multi_core=4**: uses only 4 cores on your server when running the job
  - **-quiet=stdout**: standard verbose output
 
*******************
T-Coffee Simple MSA
*******************
The central part of the web server is the T-Coffee aligner, you can use it to align any kind of sequences (Protein, RNA or DNA alike). The command it runs is the following:

::

  $#: t_coffee -in=data_93c5fbb0.in -mode=regular -output=score_html clustalw_aln fasta_aln \
      score_ascii phylip -maxnseq=150 -maxlen=10000 -case=upper -seqnos=off -outorder=input \
      -run_name=result -multi_core=4 -quiet=stdout

 
*****************
Protein Sequences
*****************
Expresso (using 3D structures)
==============================

::

  $#: t_coffee -in=data_93c5fbb0.in -mode=expresso -blast=LOCAL -pdb_db=/db/pdb/derived_data_\
      format/blast/2016-01-01/pdb_seqres.fa -evaluate_mode=t_coffee_slow -output=score_html \
      clustalw_aln fasta_aln score_ascii phylip -maxnseq=150 -maxlen=2500 -case=upper -seqnos= \
      off -outorder=input -run_name=result -multi_core=4 -quiet=stdout


M-Coffee (combining multiple methods)
=====================================

::

  $#: t_coffee -in=data_93c5fbb0.in  Mpcma_msa Mmafft_msa Mclustalw_msa Mdialigntx_msa Mpoa_msa \
      Mmuscle_msa Mprobcons_msa Mt_coffee_msa -output=score_html clustalw_aln fasta_aln \
      score_ascii phylip -tree -maxnseq=150 -maxlen=2500 -case=upper -seqnos=off -outorder=input \
      -run_name=result -multi_core=4 -quiet=stdout
      
    
PSI/TM-Coffee (transmembrane proteins)
======================================
Two options are available (in addition to the choice of the database): without transmembrane prediction or with prediction (displayed color coded on the html file).

::

  $#: tmcoffee.sh -in data_9df741d4.in -mode psicoffee -blast_server LOCAL --search-db 'UniRef50 \
      -- Very Fast/Rough' --search-type '' -prot_min_sim 50 -prot_max_sim 90 -prot_min_cov 70 \
      --search-out 'clustalw_aln fasta_aln score_ascii phylip score_html' -maxnseq 1000 -maxlen \
      =5000 -case upper -seqnos=off -outorder input -run_name result -multi_core 4 -quiet=stdout

  $#: tmcoffee.sh -in data_9df741d4.in -mode psicoffee -blast_server LOCAL --search-db 'UniRef50 \
      -- Very Fast/Rough' --search-type 'transmembrane' -prot_min_sim 50 -prot_max_sim 90 \
      -prot_min_cov 70 --search-out 'clustalw_aln fasta_aln score_ascii phylip score_html' 
      -maxnseq 1000 -maxlen 5000 -case upper -seqnos off -outorder input -run_name result º
      -multi_core 4 -quiet=stdout


PSI-Coffee (homology extension)
===============================

::

  $#: t_coffee -in=data_93c5fbb0.in -mode=psicoffee -blast=LOCAL -protein_db=/db/ncbi/201511/ \
      blast/db/nr.fa -output=score_html clustalw_aln fasta_aln score_ascii phylip -maxnseq=150 \
      -maxlen=2500 -case=upper -seqnos=off -outorder=input -run_name=result -multi_core=4 \
      -quiet=stdout


*************
RNA Sequences
*************
R-Coffee (using 2D prediction)
==============================

::

  $#: t_coffee -in=data_29091222.in -method=mafft_msa muscle_msa probconsRNA_msa -output= \
      score_html clustalw_aln fasta_aln score_ascii phylip -maxnseq=150 -maxlen=2500 -case=upper \
      -seqnos=off -outorder=input -run_name=result -tree -special_mode=rcoffee -method_limits= \
      consan_pair 5 150 -multi_core=4 -quiet=stdout
      
SARA-Coffee (using 3D structures)
==============================
SARA-Coffee is a bit complicated to run, it uses several third party packages and is run through a script and environment variables we set up; here is what it looks like: 

::

  ##: export X3DNA=/data/www-cn/sara_coffee_package/X3DNA; 
  ##: export PDB_DIR=/data/www-cn/sara_coffee_package/PDBdir/; 
  ##: export NO_REMOTE_PDB_DIR=1; 
  ##: unset MAFFT_BINARIES;
  ##: cd $CACHE_4_TCOFFEE
  ##: ln -s /data/www-cn/sara_coffee_package/pdb_entry_type.txt);
  $#: t_coffee -in data_3e6e7aec.in -method sara_pair -template_file \
      /data/www-cn/sara_coffee_package/TEMPLATEFILE,RNA -extend_mode rna2 -relax_lib 0 -transform \
      dna2rna -run_name=result -output score_html clustalw_aln -case=upper -seqnos=off -outorder= \
      input -multi_core=4 -pdb_min_sim 0 -quiet stdout
 
RM-Coffee (combining multiple methods) 
======================================
Not yet available...

*************
DNA Sequences
*************
M-Coffee (combining multiple methods) 
=====================================
For now, M-Coffee by default is the same for DNA, RNA and protein sequences alike. There is no specific M-Coffee for DNA sequences.

Pro-Coffee (homologous promoter regions)
========================================
::

  $#: t_coffee -in=data_476efe5f.in -mode=procoffee -output=score_html clustalw_aln fasta_aln \
      score_ascii phylip -maxnseq=150 -maxlen=10000 -case=upper -seqnos=off -outorder=input \
      -run_name=result -multi_core=4 -quiet=stdout


****************
Evaluation Tools
****************
TCS (Transitive Consistency Score)
==================================

::

  ##: tcs.sh -infile data_a98d61a6.in -in Mproba_pair -score 1 -output clustalw_aln fasta_aln \
      phylip score_ascii tcs_weighted tcs_replicate score_html -maxnseq 1000 -maxlen 8000 \
      -seqnos=off -run_name result -multi_core 4 --filter-type column --filter-min 4 --filter-max \
      9 --filter-gap yes -quiet=stdout

iRMSD/APDB (MSA structural evaluation) (under maintenance...)
======================================

::

  $#: t_coffee -other_pg apdb -aln data_c7151320.in -apdb_outfile default -outfile default \
      -io_format hsg3 -output score_html -maximum_distance 10 -md_threshold 2.0 -similarity_ \
      threshold 70 -template_file EXPRESSO -run_name result -quiet stdout
      
T-RMSD (structural clustering)
==============================

::

  $#: t_coffee -in=data_b89d3438.in -mode=expresso -cache=$PWD -blast=LOCAL -pdb_db=/db/pdb/ \
      derived_data_format/blast/2016-01-01/pdb_seqres.fa -evaluate_mode=t_coffee_slow -output= \
      aln score_html -maxnseq=150 -maxlen=2500 -case=upper -outorder=input -run_name=result \
      -multi_core=4 -quiet=stdout; t_coffee -other_pg trmsd result.aln -template_file \
      result_pdb1.template_list -output color_html 2>&1; [ -e result.struc_tree.consensus ]

STRIKE (MSA evaluation with single structure) 
=============================================

::
  $#: wget 	ftp://ftp.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt
  $#: install   ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST

  $#: t_coffee -other_pg strike <sequence file> -template_file PDB  -pdb_db <pdb_seqres> -blast_server LOCAL







