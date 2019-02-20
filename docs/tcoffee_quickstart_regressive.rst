################################
Quick Start Regressive Algorithm
################################

************
Introduction
************

This document introduces the regressive mode of T-Coffee. It is meant to align very large datasets with a high accuracy. In order to ise it, you need a complete installation of T-Coffee with all supported third-party packages. These come by default in the installation package but you will have to install them manualy if you compile T-Coffee from its source code. See the `installation section <https://tcoffee.readthedocs.io/en/latest/tcoffee_installation.html#installation>`_ for more details.

Multiple sequence alignments (MSA) are usualy estimated by progressively aligning all the sequences, starting with the most similar. The regressive alignment is a new procedure that proceeds the other way around and starts by aligning the most distantly related sequences first. Like its progressive sibblin, the regressive algorithm starts with a guide tree. It uses this guide tree to extract the N most diverse sequences. In this first intermediate MSA, each sequence is either a leaf or the representative of a subtree. The algorithm is re-aplied recursively onto every representative sequence until all sequences have been incorporated in an internediate MSA of max size N. The final MSA is then obtained by merging all the intermediate MSAs into the final MSA. The merging is very efficient and does not require further alignment because the intermediate MSAs contain common sequences. 

Aside from its improved accuracy, the main benefit of this protocol is to allow any third party multiple sequence aligner to be applied on the sub-groups (ClustalO, MAfft, etc) and any third party guide-tree protocol. The default T-Coffee supports a large number of `internal and external package <https://tcoffee.readthedocs.io/en/latest/tcoffee_main_documentation.html#internal-external-methods>`_ but if the one you need is not supported, note that it is relatively simple to create a configuration file allowing T-Coffee to interact with the package you need as docimented `here <https://tcoffee.readthedocs.io/en/latest/tcoffee_main_documentation.html#advanced-method-integration>`_. 

When deployed serialy, the regressive algorithm scales linearly with the number of sequences both in time and memory. The amount of memory it needs is roughly the amount of memory required to align N sequences with the selected alignment method. Given TOT sequences, the amount of time is roughly TOT/N*(amount of time required to align N Seq with the selected method). The algorithm implementation is multi-threaded and owing to its embarrasingly parralel nature, it benefits from a very efficient speedup when deployed this way.

In the following document, we list its main options and how these can be deployed from the T-Coffee algorithm. 

.. note:: REG does not support ALL the T-Coffee options. It is especially limited with respect to input and output and will only allow FASTA.

************
Installation
************

::

  ##: get the latest version from https://github.com/cbcrg/tcoffee
  ##: cd tcoffee-master/t_coffee/src
  ##: make t_coffee
  ##: mv t_coffee <Binary Directory>

In order to use T-Coffee you must also have the following pacaked installed

::
  
  ## Mafft:	 	https://mafft.cbrc.jp/alignment/software/
  ## ClustalOmega:      http://www.clustal.org/omega/#Download

********
Exemples
********

Fast and accurate
=================

In the following exemple, the regressive alignment is used to align the sequences in FASTA format. The tree is estimated using the mbed method of Clustal Omega (-reg_tree=mbed), the size of the groups is 100 (-reg_nseq=100) and the method used to align the groups is Clustal Omega:

::

  $$: t_coffee -reg -seq proteases_large.fasta -reg_nseq 100 -reg_tree mbed -reg_method clustalo_msa -outfile proteases_large.aln -outtree proteases_large.mbed

This mode is the default mode and is expected to deliver fast and reasonnably accurate alignments 

Slower and more accurate
========================
This mode is expected to be the most accurate and combines the ginsi mode of MAFFT which is very accurate with the mbed trees of Clustal Omega. Note that it is the regressive deployment that allows ginsi to be used on datasets of unlimited size.

::

  $$: t_coffee -reg -seq proteases_large.fasta -reg_nseq 100 -reg_tree mbed -reg_method mafftginsi_msa -outfile proteases_large.aln -outtree proteases_large.mbed

Very Fast
=========
This mode is expected to be the fastest currently available. Its accuracy is comparable to that of MAFFT-fftnsi running on its own 

::

  $$: t_coffee -reg -seq proteases_large.fasta -reg_nseq 100 -reg_tree parttree -reg_method mafftfftnsi_msa -outfile proteases_large.aln -outtree proteases_large.parttree

****************
Regressive flags
****************

- **-seq** (usage:**-seq=<FASTA sequence file>**)
This flag is mandatory to provide the sequences. The sequences must be in FASTA and must not contain any sequences. Because these sequences are meant to be aligned with third party aligners that may not support the IUPAC extended alphabet, it is advisable to avoid non standard amino-acids symbols. It is also not advisable to provide gapped sequences. 

- **-reg_tree** (usage:**-reg_tree=<method or Newick File>**, default: **mbed** )
This flag defines which method will be used to estimate the tree. The following methods are available

::

  Tree Computation Method:
  - mbed 	: use mBed mode of ClustalO - Default
  - cwdnd 	: use the quicktree mode of ClustalW
  - parttree 	: parttree method of MAFFT - fastest option. Does not support sequences less than 6 AA long	 
  - dpparttree 	: MAFFT fast clustering method
  - fastparttree: MAFFT fast clustering method
  - mafftdnd    : default MAFFT NJ tree - slower than the parttree modes
  - fftns1dnd   : Tree produced after the first iteration MAFFT fftns mode
  - fftns2dnd   : Tree produced after the second iteration MAFFT fftns mode
  - upgma       : upgma tree - warning cubic time computation
  - nj          : Neighbour Joinning tree
  - #<command>  : Runs comamnd <seq> > <tree>. 
  - filename    : Any file in newick format. The seq file and the tree file must match

- **-newtree** (usage:**-newtree=<filename>** , default: <infile>.<reg_tree>)
This flag defines the name of the newly computed ouput tree. Deafult will be filename.reg_tree

- **-outfile**(usage:**-outfile=<filename>** , default: <infile>.aln)
This flag defines the name of the output file containing the multiple sequence alignment


- **-reg_nseq** (usage:**-reg_nseq=N** , default: 1000 for datasets larger than 10,000 and Nseq/10 for smaller datasets )
Sets the maximum size of the subsequence alignments. The recommanded value is 1000. With slow/accurate aligners that do not scale in a linear way, this parameter can have an importnat impact on CPU requirement with small values resulting in faster computation.

- **-reg_thread** (usage:**-reg_nseq=N** , default: 0 to use all available threads)
Sets the maximum number of threads to be used by one instance. Note that if you have scrited your own aligner using a configuration file, you should make sure it does not run in a multi-threaded version as well. 

- **-reg_method**(usage:**-reg_tree=<method or configuration file>** , default: clustalo_msa)
This flag defines which method will be used to estimate the tree. In order to know which methods are available, type he following command line:

::

  $$: t_coffee

All methods the multiple sequence alignment methods xxx_msa are supported.

If you want to use an non-supported method, follow these `guidelines <https://tcoffee.readthedocs.io/en/latest/tcoffee_main_documentation.html#advanced-method-integration>`_. 

