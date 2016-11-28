###########
Quick Start
###########
.. warning:: This chapter has been extensively updated in 11/2016. All T-Coffee modes/tools/commands should be valid for version 9.03 and above, but it is not guaranteed for lower versions of T-Coffee.

******************************
Basic Command Lines (or modes)
******************************
This chapter is a quick overview on how to run T-Coffee alignment using predefined procedures we call "modes". All the files mentioned can be found `here <https://github.com/cbcrg/tcoffee/tree/master/examples>`_. You can also use examples associated with their corresponding command lines from the subsection **Tutorial (practical examples)** published in Nature Protocols (2011). Refer to the section **T-Coffee Manual** for more technical details about T-Coffee usage and tools. Please, use the corresponding citation when using a specific mode or tools of T-Coffee (see **References and Citations**) otherwise cite T-Coffee. 


Protein sequences
=================
::

  ----------------------------------------------------------------------------
  Default              
  ----------------------------------------------------------------------------
  
  $$: t_coffee sample_seq1.fasta
                                              
  Citation: Notredame et al., JMB (2000)                      PMID:10964570   
  ----------------------------------------------------------------------------
  Fast                 
  ----------------------------------------------------------------------------  
  
  $$: t_coffee sample_seq1.fasta -mode quickaln
                       
  Citation: Notredame et al., JMB (2000)                      PMID:10964570
  ----------------------------------------------------------------------------
  Consistent (M-Coffee combines the most common MSA packages)
  ----------------------------------------------------------------------------
  
  $$: t_coffee sample_seq1.fasta -mode mcoffee      

  Citation: Wallace et al., Nucleic Acids Res. (2006)         PMID:16556910
  ----------------------------------------------------------------------------
  Structure (Expresso finds structures homologous to your sequences)         
  ----------------------------------------------------------------------------
  
  $$: t_coffee sample_seq1.fasta -mode expresso

  Citation: Armougom et al. Nucleic Acids Res. (2006)         PMID:16845081
  ----------------------------------------------------------------------------
  Homology (PSI-Coffee enriches your dataset with homologous sequences)
  ----------------------------------------------------------------------------
  
  $$: t_coffee sample_seq1.fasta -mode psicoffee
  
  Citation: Chang et al., BMC Bioinformatics (2012)           PMID:22536955
  ----------------------------------------------------------------------------
  Accurate (combines Structures and Homology)            
  ----------------------------------------------------------------------------  

  $$: t_coffee sample_seq1.fasta -mode accurate
                                             
  Citation: Notredame et al., JMB (2000)                      PMID:10964570


DNA sequences
=============
::

  ----------------------------------------------------------------------------
  Default              
  ----------------------------------------------------------------------------

  $$: t_coffee sample_dnaseq1.fasta                    
                       
  Citation: Notredame et al., JMB (2000)                      PMID:10964570  
  ----------------------------------------------------------------------------
  Functional (Pro-Coffee increases accuracy of functional DNA regions )        
  ----------------------------------------------------------------------------  
  
  $$: t_coffee sample_dnaseq1.fasta -mode procoffee

  Citation: Erb et al., Nucleic Acids Res. (2012)             PMID:22230796


RNA sequences
=============
::

  ----------------------------------------------------------------------------
  Default              
  ----------------------------------------------------------------------------
  
  $$: t_coffee sample_rnaseq1.fasta                    
                       
  Citation: Notredame et al., JMB (2000)                      PMID:10964570  
  ----------------------------------------------------------------------------
  Structure 2D (R-Coffee uses predicted secondary structures)        
  ----------------------------------------------------------------------------
  
  $$: t_coffee sample_rnaseq1.fasta -mode rcoffee
  
  Citation: Wilm et al., Nucleic Acids Res. (2008)            PMID:18420654
  ----------------------------------------------------------------------------
  Structure 3D (R-Coffee combined with Consan structural alignments)
  ----------------------------------------------------------------------------  
  
  $#: t_coffee sample_rnaseq1.fasta -mode rcoffee_consan

  Citation: Wilm et al., Nucleic Acids Res. (2008)            PMID:18420654   
  ----------------------------------------------------------------------------
  Accurate (RM-Coffee use M-Coffee and secondary structure predictions)             
  ----------------------------------------------------------------------------
  
  $$: t_coffee sample_rnaseq1.fasta -mode rmcoffee
                  
  Citation: Wilm et al., Nucleic Acids Res. (2008)            PMID:18420654


********************************
Brief Overview of T-Coffee Tools
********************************
We only give you the very basics here, please go to the **T-Coffee Main Documentation** for a more detailed description. You can also try the **T-Coffee tutorial** for a practical training on T-Coffee alignment and other tools using applied examples on published research data.

Alignment methods
=================
T-Coffee
--------
Write or copy all your sequences (protein, DNA or RNA) in a given text file using one of the following format: Swiss-Prot, FASTA or PIR. Run T-Coffee with the following command line:

::

  $$: t_coffee sample_seq1.fasta


When aligning, T-Coffee will always at least generate three files:

 - ``sample_seq1.aln``  : Multiple Sequence Alignment (ClustalW format by default)
 - ``sample_seq1.dnd``  : guide tree (Newick format) 
 - ``sample_seq1.html`` : colored MSA according to consistency (html format)

In principle, the type of the sequences is automatically detected and the default methods adapted accordingly. Sometimes, however, this may fail either because the sequences are too short or contain too many ambiguity codes. When this happens, you are advised to explicitly set the type of your sequences using the flag **-type**.

::

  $$: t_coffee sample_dnaseq1.fasta -type=dna


.. note:: Please cite: Notredame, C., Higgins, D.G., Heringa, J. T-Coffee: a novel method for fast and accurate multiple sequence alignment. J. Mol. Biol., 302(1):205-217 (2000), PMID:10964570 and/or Magis, C., Taly, J.-F., Bussotti, G., Chang, J.M., Di Tommaso, P., Erb, I., Espinosa-Carrasco, J., Notredame, C. **T-Coffee: tree-based consistency objective function for alignment evaluation**. Methods Mol. Biol., 1079:117-129 (2014), PMID:24170398


M-Coffee
--------
M-Coffee is a meta version of T-Coffee that combines the output of eight aligners (MUSCLE, ProbCons, POA, DIALIGN-T, MAFFT, ClustalW, PCMA and T-Coffee); when installing T-Coffee, all required packages are automatically installed on your computer. To use M-Coffee, write your sequences in a file (format: Swiss-Prot, FASTA or PIR) and run the following command 1. M-Coffee is a predefined combination of different types of aligners; there is a faster version called fm-Coffee (command 2) which combines the fastest aligners (Kalign, MUSCLE and MAFFT). Finally, the user can make its own combination of aligners included in T-Coffee by specifying the list of packages to be combined; here is an example of T-Coffee combining ClustalW, Kalign and ProbCons (command 3).

::

  Command 1: running M-Coffee
  $$: t_coffee sample_seq1.fasta -mode mcoffee

  Command 2: running fm-Coffee
  $$: t_coffee sample_seq1.fasta -mode fmcoffee

  Command 3: user defined multiple methodes
  $$: t_coffee sample_seq1.fasta -method clustalw_pair, kalign_pair, probcons_pair
  

.. warning:: If the program starts complaining one package or the other is missing, this means you will have to go the hard way and install all these packages yourself...

.. note:: Please cite: Wallace, I.M., O'Sullivan, O., Higgins, D.G., Notredame, C. **M-Coffee: combining multiple sequence alignment methods with T-Coffee**. Nucleic Acids Res., 34(6):1692-1699 (2006), PMID:16556910


Expresso
--------
The default installation of T-Coffee provides you with the EBI wublast.pl client required to run Expresso ) command 1). Using this, Expresso will BLAST your sequences against the PDB database, identify the best targets (by default X-RAY structures, minimum 35% identical to your sequences) and use them to align your proteins using a structural aligner. If all the required structural packages for Expresso are not installed or if you want to select another structural aligner, you can select the structural package you want to use, for instance, if can use TM-align rather than SAP (command 2).

::

  Command 1: 
  $$: t_coffee sample_seq1.fasta -mode expresso

  Command 2:
  $$: t_coffee sample_seq1.fasta -template_file PDB -method TMalign_pair


This correspondence between sequences and structures (templates) is declared in a FASTA-like file we call template file. Expresso automatically generates the template file (``<your file name>_pdb1.template_list``) that can be reused for applications, but you can also provide your own with the following format. This template file should have the following format:

::

  > <seq_name> _P_ <PDB structure file or name>

  ******* sample_3Dseq1.template *******
  >TNFR10-2  _P_ 1D4V2.pdb
  >TNFR10-3  _P_ 1D4V3.pdb
  ...
  **************************************
  

.. note:: Please cite: Armougom, F., Moretti, S., Poirot, O., Audic, S., Dumas, P., Schaeli, B., Keduas, V., Notredame. C. **Expresso: automatic incorporation of structural information in multiple sequence alignments using 3D-Coffee**. Nucleic Acids Res., 34:W604-W608 (2006), PMID:16845081

R-Coffee
--------
R-Coffee can be used to align RNA sequences, using their RNApfold predicted secondary structures (command 1). The best results are obtained by using the Consan pairwise method. If you have Consan installed (under maintenance...), you get access to one of the most accurate mode of R-Coffee (command 2). This will only work if your sequences are short enough (less than 200 nucleotides). A good alternative is the rmcoffee mode (command 3) that will run MUSCLE, ProbCons4RNA and MAFFT and then use the secondary structures predicted by RNApfold. Finally, you can also select yourself which methods should be combined by R-Coffee (command 4).

::

  Command 1: R-Coffee
  $$: t_coffee sample_rnaseq1.fasta -mode rcoffee
  
  Command 2: R-Coffee + Consan
  $#: t_coffee sample_rnaseq1.fasta -mode rcoffee_consan

  Command 3: RM-Coffee
  $$: t_coffee sample_rnaseq1.fasta -mode rmcoffee

  Command 4: user defined R-Coffee
  $$: t_coffee sample_rnaseq1.fasta -mode rcoffee -method lalign_id_pair,slow_pair

.. note:: Please cite: Wilm, A., Higgins, D.G., Notredame, C. **R-Coffee: a method for multiple alignment of non-coding RNA**. Nucleic Acids Res., 36(9):e52 (2008), PMID:18420654

Pro-Coffee
----------
Pro-Coffee is a particular mode of T-Coffee designed to align specific functional DNA sequences, in particular regulatory regions. To run Pro-Coffee by default, just use command 1. In order to adjust the quality of the alignment, Pro-Coffee allows you to modify gap penalties (gap-opening and/or gap-extension) with specific flags (command 2).

::

  Command 1: Pro-Coffee default
  $$: t_coffee sample_dnaseq1.fasta -mode procoffee

  Command 2: Pro-Coffee with modified parameters
  $$: t_coffee sample_dnaseq1.fasta -method promo_pair@EP@GOP@-60@GEP@-1

.. note:: Please cite: Erb, I., González-Vallinas, J.R., Bussotti, G., Blanco, E., Eyras, E., Notredame, C. **Use of ChIP-Seq data for the design of a multiple promoter-alignment method**. Nucleic Acids Res., 40(7):e52 (2012), PMID:22230796.


Evaluation tools
================

TCS (MSA evaluation based on consistency)
-----------------------------------------
Transitive Consistency Score (TCS) is an alignment evaluation score that makes it possible to identify the most correct positions in an MSA. It has been shown that these positions are the most likely to be structuraly correct and also the most informative when estimating phylogenetic trees. The TCS evaluation and filtering procedure is implemented in the T-Coffee package and can be used to evaluate and filter any third party MSA (including T-Coffee MSA of course!). 

It's usage is a bit tricky as it comes with a lot of different options, go to the **T-Coffee Main Documentation**, section **Evaluating Your Alignment** to have all the details about TCS.

.. note:: Please cite: Chang, J.-M., Di Tommaso, P., Notredame, C. **TCS: A new multiple sequence alignment reliability measure to estimate alignment accuracy and improve phylogenetic tree reconstruction**. Mol. Biol. Evol., 31(6), 1625–1637 (2014), PMID:24694831 and/or Chang, J.-M., Di Tommaso, P., Lefort, V., Gascuel, O., Notredame, C. **TCS: a web server for multiple sequence alignment evaluation and phylogenetic reconstruction**. Nucleic Acids Res., 43(W1):W3-6 (2015), PMID:25855806

iRMSD/APDB (MSA structural evaluation)
--------------------------------------
iRMSD/APDB is not an alignment tool, it is an evalution tool of a given alignment using structural information. All you need is a file containing the alignment of sequences with a known structure ("template file"; see Expresso). If you don't provide a template file, these sequences must be named according to their PDB ID, followed by the chain index (1aabA for example). In the first example (command 1) names are different therefore it won't deliver any result. In that case, you should declare the correspondence between sequences and structures using your own template file (command 2). All the sequences do not need to have a known structure, but at least two is required otherwise it won't deliver any result. 
::

  Command 1: 
  $$: t_coffee -other_pg irmsd sample_3Dseq1.aln

  Command 2:
  $$: t_coffee -other_pg irmsd sample_3Dseq1.aln -template_file sample_3Dseq1.template


.. note:: Please cite: Armougom, F., Moretti, S., Keduas, V., Notredame, C. **The iRMSD: a local measure of sequence alignment accuracy using structural information**. Bioinformatics, 22(14):e35-e39 (2006), PMID:16873492

STRIKE (single structure MSA evaluation)
----------------------------------------
Under maintenance on the webserver or the T-Coffee package...

T-RMSD (structural clustering)
------------------------------
T-RMSD is a structure based clustering method using the iRMSD to drive the structural clustering of your aligned sequences with an available structure. The T-RMSD supports all the parameters supported by iRMSD or APDB. To run T-RMSD, type:

::

  $$: t_coffee -other_pg trmsd sample_3Dseq1.aln -template_file sample_3Dseq1.template


The program then outputs a series of files:
 - ``sample_3Dseq1.struc_tree.list`` : list of the trees associated with every position.
 - ``sample_3Dseq1.struc_tree.html`` : colored columns supporting the tree.
 - ``sample_3Dseq1.struc_tree.consensus_output`` : schematic display of the results.
 - ``sample_3Dseq1.struc_tree.consensus`` : final consensus structural tree.

.. note:: Please cite: Magis, C., Stricher, F., van der Sloot, A.M., Serrano, L., Notredame, C. **T-RMSD: a fine-grained, structure based classification method and its application to the functional characterization of TNF receptors**. J. Mol. Biol., 400(3):605-617 (2010), PMID:20471393 and/or Magis, C., van der Sloot, A.M., Serrano, L., Notredame, C. **An improved understanding of TNFL/TNFR interactions using structure-based classifications**. Trends Biochem. Sci., 37(9):353-363 (2012), PMID:22789664


*****************************
Tutorial (Practical Examples)
*****************************

.. note:: This documentation is merely a cheat-sheet that recapitulates the material and the command lines associated with the manual. This tutorial itself is adpated from the `T-Coffee Nature Protocols Article <http://www.nature.com/nprot/journal/v6/n11/full/nprot.2011.393.html>`_ that can be followed step by step on the following `website <http://www.tcoffee.org/Projects/tcoffee/workshops/tcoffeetutorials/index.html>`_ 

Introduction
============
T-Coffee is a versatile Multiple Sequence Alignment method suitable for aligning most types of biological sequences. The series of protocols presented here show how the package can be used to multiply align proteins, DNA and RNA sequences. The package is an open source freeware available from `our website <http://www.tcoffee.org>`_.

There are several parts: 1) the protein section presents controlled cases for PSI-Coffee the homology extended mode suitable for remote homologues, Expresso the structure based multiple aligner and M-Coffee, a meta version able to combine several third party aligners into one, 2) we then show how the T-RMSD option can be used to produce a functionally informative structure based clustering, 3) RNA alignment procedures are shown for R-Coffee a mode that produces secondary structure based MSAs, 4) DNA alignments are illustrated with Pro-Coffee, a multiple aligner specific of promoter regions, 5) finally, the last section presents some of the many reformatting utilities bundled with T-Coffee. 

Materials
=========
The list of files (input and output) required by this protocol is available from `here <http://www.tcoffee.org/Packages/NatureProtocols/NatureProtocolDataset.tar.gz>`_. They can be automatically retrieved using the following command:

::

  $$: t_coffee -other_pg nature_protocol.pl    

This will create 4 repertories containing the input sequences necessary for the protocols we report in this section. For each part, all command lines have been collected into the file README.sh.

Procedures
==========
- `Full Tutorial <http://www.tcoffee.org/Projects/tcoffee/workshops/tcoffeetutorials/index.html>`_
- `Installation <http://www.tcoffee.org/Projects/tcoffee/workshops/tcoffeetutorials/installation.html>`_
- `Protein Multiple Sequence Alignments <http://www.tcoffee.org/Projects/tcoffee/workshops/tcoffeetutorials/protein-alignment.html>`_
- `RNA Multiple Sequence Alignments <http://www.tcoffee.org/Projects/tcoffee/workshops/tcoffeetutorials/rna-alignment.html>`_
- `Promoter alignments <http://www.tcoffee.org/Projects/tcoffee/workshops/tcoffeetutorials/promoter-alignment.html>`_
- `Reformat alignments <http://www.tcoffee.org/Projects/tcoffee/workshops/tcoffeetutorials/reformating.html>`_
