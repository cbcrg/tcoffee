#####################
Regressive Quickstart
#####################

.. note:: This documentation is merely a cheat-sheet that recapitulates the material and the command lines associated with the `_ 

************
Introduction
************

This document introduces the regressive mode of T-Coffee. It is meant to align very large datasets with a high accuracy. In order to ise it, you need a complete installation of T-Coffee with all supported third-party packages. These come by default in the installation package but you will have to install them manualy if you compile T-Coffee from its source code. See the `installation section <https://tcoffee.readthedocs.io/en/latest/tcoffee_installation.html#installation>`_ for more details.

Multiple sequence alignments (MSA) are usualy estimated by progressively aligning all the sequences, starting with the most similar. The regressive alignment is a new procedure that proceeds the other way round and starts by aligning the most distantly related sequences first. Like its progressive sibblin, the regressive algorithm starts with a guide tree. It uses this guide tree to extract the N most diverse sequences. In this first intermediate MSA, each sequence is either a leaf or the representative of a subtree. The algorithm is re-aplied recursively onto very representative sequence until all sequences have been incorporated in an internediate MSA of max size N. The final MSA is then obtained by merging all the intermediate MSAs into one large model. The merging is very efficient and does not require further alignment because the intermediate MSAs contain common sequences. 

Aside from its improved accuracy, the main benefit of this protocol is to allow any third party multiple sequence aligner to be applied on the sub-groups (ClustalO, MAfft, etc) and any third party guide-tree proticol. The default T-Coffee supports a large number of `internal and external package <https://tcoffee.readthedocs.io/en/latest/tcoffee_main_documentation.html#internal-external-methods>`_ but if the one you need is not supported, note that it is relatively simple to create a configuration file allowing T-Coffee to interact with the package you need as docimented `here <https://tcoffee.readthedocs.io/en/latest/tcoffee_main_documentation.html#advanced-method-integration>`_. 

The regressive implmentation of T-Coffee is named DPA. In the following document, we list its main options and how these can be deployed from the T-Coffee algorithm. Note that DPA is a simplified wrapper around T-Coffee.It does not support all the T-Coffee flags and its output/input formats are limited to FASTA for the sequences and Newick for the trees.

.. note:: DPA does not support ALL the T-Coffee options. It is especially limited with respect to input and output and will only allow FASTA.



********
Exemples
********

Fast and accurate
=================

In the following exemple, the regressive alignment is used to align the sequences in FASTA format. The tree is estimated using the mbed method of Clustal Omega (-dpa_tree=mbed), the size of the groups is 100 (-dpa_nseq=100) and the method used to align the groups is Clustal Omega:

  $$: t_coffee -dpa -seq proteases_large.fasta -dpa_nseq 100 -dpa_tree mbed -dpa_method clustalo_msa -outfile proteases_large.aln -outtree proteases_large.mbed

This mode is the default mode and is expected to deliver fast and reasonnably accurate alignments 

Slower and more accurate
========================
This mode is expected to be the most accurate and combines the ginsi mode of MAFFT which is very accurate with the mbed trees of Clustal Omega. Note that it is the regressive deployment that allows ginsi to be used on datasets of unlimited size.

  $$: t_coffee -dpa -seq proteases_large.fasta -dpa_nseq 100 -dpa_tree mbed -dpa_method mafftginsi_msa -outfile proteases_large.aln -outtree proteases_large.mbed

Very Fast
=========
This mode is expected to be the fastest currently available. Its accuracy is comparable to that of MAFFT-fftnsi running on its own 

  $$: t_coffee -dpa -seq proteases_large.fasta -dpa_nseq 100 -dpa_tree parttree -dpa_method mafftfftnsi_msa -outfile proteases_large.aln -outtree proteases_large.parttree

********************
Regressive DPA flags
********************

- **-seq** (usage:**-seq=<FASTA sequence file>**)
This flag is mandatory to provide the sequences. The sequences must be in FASTA and must not contain any sequences. Because these sequences are meant to be aligned with third party aligners that may not support the IUPAC extended alphabet, it is advisable to avoid non standard amino-acids symbols. It is also not advisable to provide gapped sequences. 

- **-dpa_tree**(usage:**-dpa_tree=<method or Newick File>**)
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

- **-newtree**(usage:**-newtree=<filename>**)
This flag defines the name of the newly computed ouput tree. Deafult will be filename.dpa_tree

- **-outfile**(usage:**-outfile=<filename>**)
This flag defines the name of the output file containing the multiple sequence alignment


- **-dpa_nseq** (usage:**-dpa_nseq=N**) [cw]
Sets the maximum size of the subsequence alignments. The recommanded value is 1000. With slow/accurate aligners that do not scale in a linear way, this parameter can have an importnat impact on CPU requirement with small values resulting in faster computation.

- **-dpa_method**(usage:**-dpa_tree=<method or configuration file>**)
This flag defines which method will be used to estimate the tree. In order to know which methods are available, type he following command line:

::

  $$: t_coffee

All methods the multiple sequence alignment methods xxx_msa are supported.

If you want to use an non-supported method, follow these `guidelines <https://tcoffee.readthedocs.io/en/latest/tcoffee_main_documentation.html#advanced-method-integration>`_. 

