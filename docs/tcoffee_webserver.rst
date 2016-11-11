###################
T-Coffee Web Server
###################

In this chapter we describe briefly the available T-Coffee web servers and their usage (command lines and examples files); if you need more advanced options, you should use the T-Coffee package via command lines. The web server is a reduced version of the T-Coffee package, containing all T-Coffee modes for aligning protein, RNA or DNA sequences; it also contains the evaluation and downstream analysis tools (TCS, iRMSD/APDB, STRIKE, T-RMSD). Currently, all reformatting utilities are not available on the web server, however, you can choose some reformatting options related to the ouput format. Here we present briefly the webserver and the T-Coffee commands it contains.

*******************
T-Coffee Simple MSA
*******************
The central part of the web server is the T-Coffee aligner, you can use it to align any kind of sequences (Protein, RNA or DNA alike). The command it runs is the following:

::

  $#: t_coffee -in=data_93c5fbb0.in -mode=regular -output=score_html clustalw_aln fasta_aln \
      score_ascii phylip -maxnseq=150 -maxlen=10000 -case=upper -seqnos=off -outorder=input \
      -run_name=result -multi_core=4 -quiet=stdout

 - **-mode=regular**: uses the internal **proba_pair** aligner
 - **-output**: specifies the output file format you require
 - **-maxnseq**: restriction on the number of sequences, max=150
 - **-maxlength**: restriction on the length of your sequences, max=1000
 - **-outorder=input**: orders the sequences in the final MSA as in the input dataset 
 - **-multi_core=4**: uses only 4 cores on your server when running the job
 
*****************
Protein Sequences
*****************
Expresso (using 3D structures)
==============================

::

##: t_coffee -in=data_93c5fbb0.in -mode=expresso -blast=LOCAL -pdb_db=/db/pdb/derived_data_format\
/blast/2016-01-01/pdb_seqres.fa -evaluate_mode=t_coffee_slow -output=score_html clustalw_aln \
fasta_aln score_ascii phylip -maxnseq=150 -maxlen=2500 -case=upper -seqnos=off -outorder=input \
-run_name=result -multi_core=4 -quiet=stdout


M-Coffee (combining multiple methods)
=====================================

::

##: t_coffee -in=data_93c5fbb0.in  Mpcma_msa Mmafft_msa Mclustalw_msa Mdialigntx_msa Mpoa_msa \
Mmuscle_msa Mprobcons_msa Mt_coffee_msa -output=score_html clustalw_aln fasta_aln score_ascii \
phylip -tree -maxnseq=150 -maxlen=2500 -case=upper -seqnos=off -outorder=input -run_name=result \
-multi_core=4 -quiet=stdout
      
    
PSI/TM-Coffee (transmembrane proteins)
======================================


PSI-Coffee (homology extension)
===============================


*************
RNA Sequences
*************
R-Coffee (using 2D prediction)
==============================

SARA-Coffee (using 3D structures)
==============================

RM-Coffee (combining multiple methods)
======================================


*************
DNA Sequences
*************
M-Coffee (combining multiple methods)
=====================================

Pro-Coffee (homologous promoter regions)
========================================


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








