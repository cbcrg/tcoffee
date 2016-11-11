###################
T-Coffee Web Server 
###################

.. warning:: This chapter is currently under maintenance. It aims is to describe what the T-Coffee web server is doing programmatically when you are using it; it provides technical details about the restrictions and commands in use.

In this chapter we describe briefly the available T-Coffee web servers and their usage (command lines and examples files); if you need more advanced options, you should use the T-Coffee package via command lines. The web server is a reduced version of the T-Coffee package, containing all T-Coffee modes for aligning protein, RNA or DNA sequences; it also contains the evaluation and downstream analysis tools (TCS, iRMSD/APDB, STRIKE, T-RMSD). Currently, all reformatting utilities are not available on the web server, however, you can choose some reformatting options related to the ouput format. Here we present briefly the webserver and the T-Coffee commands it contains.

*******************
T-Coffee Simple MSA
*******************
The central part of the web server is the T-Coffee aligner, you can use it to align any kind of sequences (Protein, RNA or DNA alike). The command it runs is the following:

::

  $#: t_coffee -in=data_93c5fbb0.in -mode=regular -output=score_html clustalw_aln fasta_aln \
      score_ascii phylip -maxnseq=150 -maxlen=10000 -case=upper -seqnos=off -outorder=input \
      -run_name=result -multi_core=4 -quiet=stdout


The different options correspond to:
  - **-mode=regular**: uses the internal **proba_pair** aligner
  - **-output**: specifies the output file format you require
  - **-maxnseq**: restriction on the number of sequences, max=150
  - **-maxlength**: restriction on the length of your sequences, max=1000
  - **-case=upper**: all residues/nucleotids will be upper case
  - **-seqnos=off**: the size of the sequences is not indicated on the MSA
  - **-outorder=input**: orders the sequences in the final MSA as in the input dataset 
  - **-run_name**: name of the job on the cluster
  - **-multi_core=4**: uses only 4 cores on your server when running the job
  - **-quiet=stdout**: standard verbose output
 
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
      Mmuscle_msa Mprobcons_msa Mt_coffee_msa -output=score_html clustalw_aln fasta_aln score_ascii \
      phylip -tree -maxnseq=150 -maxlen=2500 -case=upper -seqnos=off -outorder=input -run_name \
      =result -multi_core=4 -quiet=stdout
      
    
PSI/TM-Coffee (transmembrane proteins)
======================================
Two options are available (in addition to the choice of the database): without transmembrane prediction or with prediction (displayed color coded on the html file).

::

  $#: tmcoffee.sh -in data_9df741d4.in -mode psicoffee -blast_server LOCAL --search-db 'UniRef50 \
      -- Very Fast/Rough' --search-type '' -prot_min_sim 50 -prot_max_sim 90 -prot_min_cov 70 \
      --search-out 'clustalw_aln fasta_aln score_ascii phylip score_html' -maxnseq 1000 -maxlen \
      5000 -case upper -seqnos=off -outorder input -run_name result -multi_core 4 -quiet=stdout

$#: tmcoffee.sh -in data_9df741d4.in -mode psicoffee -blast_server LOCAL --search-db 'UniRef50 \
    -- Very Fast/Rough' --search-type 'transmembrane' -prot_min_sim 50 -prot_max_sim 90 -prot_min_cov \
    70 --search-out 'clustalw_aln fasta_aln score_ascii phylip score_html' -maxnseq 1000 -maxlen 5000 \
    -case upper -seqnos off -outorder input -run_name result -multi_core 4 -quiet=stdout


PSI-Coffee (homology extension)
===============================

::

    $#: t_coffee -in=data_93c5fbb0.in -mode=psicoffee -blast=LOCAL -protein_db=/db/ncbi/201511/blast/ \
        db/nr.fa -output=score_html clustalw_aln fasta_aln score_ascii phylip -maxnseq=150 -maxlen=2500 \
        -case=upper -seqnos=off -outorder=input -run_name=result -multi_core=4 -quiet=stdout


*************
RNA Sequences
*************
R-Coffee (using 2D prediction)
==============================

::

  $#: t_coffee -in=data_29091222.in -method=mafft_msa muscle_msa probconsRNA_msa -output= \
      score_html clustalw_aln fasta_aln score_ascii phylip -maxnseq=150 -maxlen=2500 -case=upper \
      -seqnos=off -outorder=input -run_name=result -tree -special_mode=rcoffee -method_limits=consan_pair \
      5 150 -multi_core=4 -quiet=stdout
      
SARA-Coffee (using 3D structures)
==============================

::

  export X3DNA=/data/www-cn/sara_coffee_package/X3DNA; 
  export PDB_DIR=/data/www-cn/sara_coffee_package/PDBdir/; 
  export NO_REMOTE_PDB_DIR=1; 
  unset MAFFT_BINARIES;
  (cd $CACHE_4_TCOFFEE; ln -s /data/www-cn/sara_coffee_package/pdb_entry_type.txt);
  $#: t_coffee -in data_3e6e7aec.in -method sara_pair -template_file \
      /data/www-cn/sara_coffee_package/TEMPLATEFILE,RNA -extend_mode rna2 -relax_lib 0 -transform \
      dna2rna -run_name=result -output score_html clustalw_aln -case=upper -seqnos=off -outorder= \
      input -multi_core=4 -pdb_min_sim 0 -quiet stdout
 
 
RM-Coffee (combining multiple methods) (uner maintenance...)
======================================



*************
DNA Sequences
*************
M-Coffee (combining multiple methods) (under maintenance...)
=====================================

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


iRMSD/APDB (MSA structural evaluation)
======================================


T-RMSD (structural clustering)
==============================


STRIKE (MSA evaluation with single structure) (under maintenance...)
=============================================








