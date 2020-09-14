####################
Quick Start Unistrap
####################

************
Introduction
************

This document introduces the unsitsrap strategy that incorpates MSA instability into the estimation of the Felsentein Bootstrap replicates (please cite the `Systematic Biology original publication <https://academic.oup.com/sysbio/article/67/6/997/4948750>`_). 

The principke is very straightorward. It involves shuffling the input order of the sequences so as to generate MSA shuffle replicates, from which columns are then drawn with replacement so as to generate the Bootsrap replicates. As established in the paper, it is also possible to shuffle left and right chidren in the guide tree of any progressive methods. This procedure recapitulates two third of the instability measured when shuffling sequences. Whenever using a progressive algorithm, the tree shuffling procedure should be applied in order to by-pass any arbitrary stabilisation procedures (internal sorting for instance).

The procedure used for the original validation is available on `Github <https://github.com/cbcrg/unistrap>`_). The content of this repository is not especially user friendly and its sole purpose is to allow reproducing the analysis provided in the original paper. For practical usage, we recommand using the re-implementation provided below that can be accessed using the T-Coffee package. This re-implementation supports all the T-Coffee supported aligners when shuffling input sequences (**t_coffee ** to get a list) and a subset of progressive aligners when using the tree shuffle strategy. 


Aside from its improved accuracy, the main benefit of this protocol is to allow any third party multiple sequence aligner to be applied on the sub-groups (ClustalO, MAfft, etc) and any third party guide-tree protocol. The default T-Coffee supports a large number of `internal and external package <https://tcoffee.readthedocs.io/en/latest/tcoffee_main_documentation.html#internal-external-methods>`_ but if the one you need is not supported, note that it is relatively simple to create a configuration file allowing T-Coffee to interact with the package you need as docimented `here <https://tcoffee.readthedocs.io/en/latest/tcoffee_main_documentation.html#advanced-method-integration>`_. 

When deployed serialy, the regressive algorithm scales linearly with the number of sequences both in time and memory. The amount of memory it needs is roughly the amount of memory required to align N sequences with the selected alignment method. Given TOT sequences, the amount of time is roughly TOT/N*(amount of time required to align N Seq with the selected method). The algorithm implementation is multi-threaded and owing to its embarrasingly parralel nature, it benefits from a very efficient speedup when deployed this way.

In the following document, we list its main options and how these can be deployed from the T-Coffee algorithm. 

.. note:: REG does not support ALL the T-Coffee options. It is especially limited with respect to input and output and will only allow FASTA.

***************************
Installation from binaries
**************************

::

  ##: Get the latest stable version from http://www.tcoffee.org/Packages/Stable/Latest/
  ##: Or  the latest Beta   Version from http://www.tcoffee.org/Packages/Beta/Latest/	
  ##: download the *.tar.gz file
  ##: tar -xvf T-COFFEE_distribution_Version_XXXXX.tar.gz
  ##: cd T-COFFFEE_distribution_Version_XXXXX
  ##: ./install all
  ##: add the instructions given at the bottom of the install output to your .profile or your .bashrc 



********
Examples
********

Produce bootsrap replicates by shuffling sequences and guide tree sister nodes
==============================================================================

In the following example, the input order of the sequences will be shuffled 10 times so as to generate as many MSA replicates (-msa). Each time a guide tree will be generated using mbed (-tree), and the sister nodes of this tree will be randomly shuffled (i.e. left and right will be randomléy reassigned to left or right) - note that this shuffling does not change the topology. For each shuffled replicated, an MSA willl be estimated with clustalo (-method) and each MSA this produced will be used to produce 10 replicates. 

.. note:: Note: This mode will re-compute a guide tree for every sequence shuffle. The tree shuffle is not needed but it will increase the probability of each iteration delivering a different MSA. It is much more computationaly efficient to replace the sequence order shuffle with a tree shuffle thgat will generate a comparable amount of MSA diversity at a significantly lower computational cost. 

::

  $$: t_coffee -other_pg unistrap -in sample_seq1.fasta -tree mbed -method clustalo -shuffle_tree -shuffle_seq -msa 10 -bootstrap 10 -outfile replicates.phylip



Produce bootsrap replicates by shuffling the guide tree sister nodes
====================================================================
In the following example, a guide tree will be generated once using mbed (-tree), and the sister nodes of this tree will be randomly shuffled (i.e. left and right will be randomléy reassigned to left or right) 10 times. For each shuffled replicated, an MSA willl be estimated with clustalo (-method) and each MSA this produced will be used to produce 10 replicates. 


::

  $$: t_coffee -other_pg unistrap -in sample_seq1.fasta -tree mbed -method clustalo -shuffle_tree -msa 10 -bootstrap 10 -outfile replicates.phylip

Get list of supported methods
=============================
The following command will list all the guide tree methods and multiple aligners supported by unistrap

::

  $$: t_coffee -other_pg unistrap  -tree list -method list 



****************
unistrap flags
****************

- **-in** (usage:**-seq=<FASTA sequence file>**)
This flag is mandatory to provide the sequences. The sequences must be in FASTA and must not contain any sequences. Because these sequences are meant to be aligned with third party aligners that may not support the IUPAC extended alphabet, it is advisable to avoid non standard amino-acids symbols. It is also not advisable to provide gapped sequences. 

- **-shuffle_tree** (usage:**-shuffle_tree**, default: **on** )
This flag defines the shuffling of sister nodes for each MSA replicate. The guide tree does not need to be re-estimated

- **-shuffle_seq** (usage:**-shuffle_seq**, default: **off** )
This flag defines the shuffling of input sequences for each MSA replicate. The guide tree *must be re-estimated*.

- **-tree** (usage:**-tree=<method,filename,list>** , default: mbed)
This flag defines the method used to estimate the guide tree, it can be a method or an existing tree

- **-method** (usage:**-tree=<method,list>** , default: clustalo)
This flag defines the method used to estimate the MSA.

- **-msa** (usage:**-msa=<integer>** , default: 10)
This flag defines the number of MSA replicates

- **-bootstrap** (usage:**-bootstrap=<integer>** , default: 10)
This flag defines the number of Bootsrap replicates. Note that the final number of replicates will be msa* x bootstrap
	
- **-outfile** (usage:**-outfile=<tree_method>** , default:stdout)
This flag defines the outfile name

