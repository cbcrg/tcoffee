###########
Quick Start
###########


******************************
Basic Command Lines (or modes)
******************************
.. Danger:: The T-Coffee documentation is currently under maintenance, the example files can be found here: <https://github.com/cbcrg/tcoffee/tree/master/t_coffee/doc_test/data>. You can also use examples associated with their corresponding command lines from the section **T-Coffee Tutorial** published in Nature Protocols (2011). This is temporary and everything will be working when this message will be removed. Thanks for your understanding.

.. important:: This part is a quick overview on how to run T-Coffee alignment using predefined procedures or "modes". All the files mentioned here can be found in the example directory of the distribution (currently under maintenance). If you use a specific mode of T-Coffee, please use the corresponding citation for publication. For more details about T-Coffee usage, refer to the section **T-Coffee Manual**.


Protein sequences
=================
::

  ----------------------------------------------------------------------------
  Default              t_coffee sample_seq1.fasta
                       use the output.html to visualize the MSA accuracy
                       
  Citation: Notredame et al., JMB (2000)                      PMID:10964570  
  ----------------------------------------------------------------------------
  Fast                 t_coffee sample_seq1.fasta -mode quickaln
                       lower -ndiag if the sequences are very similar

  Citation: Notredame et al., JMB (2000)                      PMID:10964570
  ---------------------------------------------------------------------------- 
  Accurate             t_coffee sample_aln1.fasta -mode accurate
                       combines structures, sequences and profiles
                       
  Citation: Notredame et al., JMB (2000)                      PMID:10964570
  ----------------------------------------------------------------------------
  Consistent           t_coffee sample_aln1.fasta -mode mcoffee
  (M-Coffee)           combines the most common existing MSA packages

  Citation: Wallace et al., Nucleic Acids Res. (2006)         PMID:16556910
  ----------------------------------------------------------------------------
  Structural           t_coffee sample_aln1.fasta -mode expresso
  (Expresso)           finds structures homologous to your sequences

  Citation: Armougom et al. Nucleic Acids Res. (2006)         PMID:16845081
  ----------------------------------------------------------------------------
  Homology             t_coffee sample_aln1.fasta -mode psicoffee
  (PSI-Coffee)         enriches your dataset with homologous sequences
  
  Citation: Chang et al., BMC Bioinformatics (2012)           PMID:22536955
  ----------------------------------------------------------------------------


DNA sequences
=============
::

  ----------------------------------------------------------------------------
  Default              t_coffee three_cdna.fasta
                       use the output.html to visualize the MSA accuracy
                       
  Citation: Notredame et al., JMB (2000)                      PMID:10964570  
  ----------------------------------------------------------------------------
  Functional           t_coffee three_cdna.fasta -mode procoffee
  (Pro-Coffee)         increases accuracy of functional DNA regions
  
  Citation: Erb et al., Nucleic Acids Res. (2012)             PMID:22230796
  ----------------------------------------------------------------------------  


RNA sequences
=============
::

  ----------------------------------------------------------------------------
  Default              t_coffee sample_rnaseq1.fasta
                       use the output.html to visualize the MSA accuracy
                       
  Citation: Notredame et al., JMB (2000)                      PMID:10964570  
  ----------------------------------------------------------------------------
  Structural 2D        t_coffee sample_rnaseq1.fasta -mode rcoffee
  (R-Coffee)           uses the predicted secondary structure
  
  Citation: Wilm et al., Nucleic Acids Res. (2008)            PMID:18420654
  ----------------------------------------------------------------------------
  Structural 3D     t_coffee sample_rnaseq1.fasta -mode 
  (R-Coffee Consan)    uses R-Coffee to combine Consan structural alignments 
  
  Citation: Wilm et al., Nucleic Acids Res. (2008)            PMID:18420654   
  ----------------------------------------------------------------------------
  Accurate             t_coffee sample_rnaseq1.fasta -mode rmcoffee
  (RM-Coffee)          use M-Coffee + secondary structure prediction
                       
  Citation: Wilm et al., Nucleic Acids Res. (2008)            PMID:18420654
  ----------------------------------------------------------------------------



********************************
Brief Overview of T-Coffee Tools
********************************

.. important:: We only give you the very basics here, please go to the T-Coffee manual for a more detailed description and available options for the different tools. You can also try the T-Coffee tutorial for a practical training on T-Coffee alignment and other functions using applied examples on published research data.


T-Coffee
========
Write or copy all your sequences (protein, DNA or RNA) in a given text file using one of the following format: Swiss-Prot, FASTA or PIR; then run T-Coffee with the following command line:


::

  $$: t_coffee sample_seq1.fasta



This will output three files:


::

  sample_seq1.aln  : your multiple sequence alignment (Clustal format by default)

  sample_seq1.dnd  : the guide tree (Newick format)
  
  sample_seq1.html : the color coded MSA according to T-Coffee consistency color scheme (html)


.. tip:: In principle, the type of the sequences is automatically detected and the default methods adapted accordingly. Sometimes, however, this may fail either because the sequences are too short or contain too many ambiguity codes. When this happens, you are advised to explicitly set the type of your sequences using the flag -type.

::

  $$: t_coffee sample_dnaseq1.fasta -type=dna


.. note:: Please cite: Notredame, C., Higgins, D.G., Heringa, J. T-Coffee: a novel method for fast and accurate multiple sequence alignment. J. Mol. Biol., 302(1):205-217 (2000), PMID:10964570 and/or Magis, C., Taly, J.-F., Bussotti, G., Chang, J.M., Di Tommaso, P., Erb, I., Espinosa-Carrasco, J., Notredame, C. T-Coffee: tree-based consistency objective function for alignment evaluation. Methods Mol. Biol., 1079:117-129 (2014), PMID:24170398


M-Coffee
========
M-Coffee is a meta version of T-Coffee that combines the output of eight aligners (MUSCLE, ProbCons, POA, DIALIGN-T, MAFFT, ClustalW, PCMA and T-Coffee); when installing T-Coffee, all required packages are automatically installed on your computer. To use M-Coffee, write your sequences in a file (format: Swiss-Prot, FASTA or PIR) and run the following command line:


::

  $$: t_coffee sample_seq1.fasta -mode mcoffee


M-Coffee is a predefined combination of different types of aligners; there is a faster version called fm-Coffee which combines the fastest aligners (Kalign, MUSCLE and MAFFT):

::

  $$: t_coffee sample_seq1.fasta -mode fmcoffee

Also, the user can make its own combination of aligners included in T-Coffee by specifying the list of packages to be combined; here is an example of T-Coffee combining ClustalW, Kalign and ProbCons:

::

  $$: t_coffee sample_seq1.fasta -method clustalw_pair, kalign_pair, probcons_pair
  
  
If the program starts complaining one package or the other is missing, this means you will have to go the hard way and install all these packages yourself...Proceed to the **T-Coffee Installation** section for more detailed instructions.


.. note:: Please cite: Wallace, I.M., O'Sullivan, O., Higgins, D.G., Notredame, C. M-Coffee: combining multiple sequence alignment methods with T-Coffee. Nucleic Acids Res., 34(6):1692-1699 (2006), PMID:16556910


Expresso
========
The default installation of T-Coffee provides you with the EBI wublast.pl client required to run Expresso. Using this, Expresso will BLAST your sequences against the PDB database, identify the best targets and use them to align your proteins using a structural aligner. Run Expresso with the following command:


::

  $$: t_coffee sample_seq1.fasta -mode expresso



If all the required structural packages for Expresso were not installed or if you want to select another structural aligner, you can select the structural package you want to use. For instance, if can use TM-align rather than SAP:


::

  $$: t_coffee sample_seq1.fasta -template_file expresso -method TMalign_pair


.. note:: Please cite: Armougom, F., Moretti, S., Poirot, O., Audic, S., Dumas, P., Schaeli, B., Keduas, V., Notredame. C. Expresso: automatic incorporation of structural information in multiple sequence alignments using 3D-Coffee. Nucleic Acids Res., 34:W604-W608 (2006), PMID:16845081


MOCCA
=====
MOCCA is a specific tool in T-Coffee designed to deal with highly divergent protein repeats.  Write your sequences in the same file and type:


::

  $$: t_coffee -other_pg mocca sample_seq1.fasta


This command output one files (<your sequences>.mocca_lib) and starts an interactive menu.


.. note:: Please cite: Notredame, C. MOCCA: semi-automatic method for domain hunting. Bioinformatics, 17(4):373-374 (2001), PMID:11301309


Pro-Coffee
==========
Pro-Coffee is a particular mode of T-Coffee designed to align specific functional DNA sequences, in particular regulatory regions. To run Pro-Coffee by default, type:


::

  $$: t_coffee three_cdna.fasta -mode procoffee
  

In order to adjust the quality of the alignment, Pro-Coffee allows you to modify gap penalties (gap-opening and/or gap-extension) using the following command line:


::

  $$: t_coffee three_cdna.fasta -method promo_pair@EP@GOP@-60@GEP@-1


.. note:: Please cite: Erb, I., González-Vallinas, J.R., Bussotti, G., Blanco, E., Eyras, E., Notredame, C. Use of ChIP-Seq data for the design of a multiple promoter-alignment method. Nucleic Acids Res., 40(7):e52 (2012), PMID:22230796.


R-Coffee
========
R-Coffee can be used to align RNA sequences, using their RNApfold predicted secondary structures. The best results are obtained by using the Consan pairwise method. If you have Consan installed, run:


::

  $$: t_coffee sample_rnaseq1.fasta -special_mode rcoffee_consan



This will only work if your sequences are short enough (less than 200 nucleotides). A good alternative is the rmcoffee mode that will run MUSCLE, ProbCons4RNA and MAFFT and then use the secondary structures predicted by RNApfold:


::

  $$: t_coffee sample_rnaseq1.fasta -mode rmcoffee



If you want to select yourself which methods should be combined by R-Coffee, run:


::

  $$: t_coffee sample_rnaseq1.fasta -mode rcoffee -method lalign_id_pair slow_pair


.. note:: Please cite: Wilm, A., Higgins, D.G., Notredame, C. R-Coffee: a method for multiple alignment of non-coding RNA. Nucleic Acids Res., 36(9):e52 (2008), PMID:18420654


iRMSD and APDB
==============
iRMSD/APDB is not an alignment tool, it is an evalution tool of a given alignment using structural information. All you need is a file containing the alignment of sequences with a known structure. These sequences must be named according to their PDB ID, followed by the chain index (1aabA for instance). All the sequences do not need to have a known structure, but at least two is required. Given the alignment, use the following command:


::

  $$: t_coffee -other_pg irmsd -aln 3d_sample5.aln


If the names of the sequences do not correspond to the PDB name, then the user have to declare the correspondence between sequences and structures in a template file (cf. **T-Coffee Manual** section):

::

  $$: t_coffee -other_pg irmsd -aln 3d_sample5.aln -template_file 3d_sample5.template_file


.. note:: Please cite: Armougom, F., Moretti, S., Keduas, V., Notredame, C. The iRMSD: a local measure of sequence alignment accuracy using structural information. Bioinformatics, 22(14):e35-e39 (2006), PMID:16873492


T-RMSD
=====
T-RMSD is a structure based clustering method using the iRMSD to drive the structural clustering of your aligned sequences with an available structure. The T-RMSD supports all the parameters supported by iRMSD or APDB. To run T-RMSD, type:


::

  $$: t_coffee -other_pg trmsd -aln 3d_sample5.aln -template_file 3d_sample5.template_list


3d_sample5.aln is a multiple alignment in which each sequence has a known structure. The file 3d_sample5.template_list is a fasta like file declaring the structure associated with each sequence, in the form:


::

  > <seq_name> _P_ <PDB structure file or name>

  ******* 3d_sample5.template_list ********

  >2UWI-3A _P_ 2UWI-3.pdb

  >2UWI-2A _P_ 2UWI-2.pdb

  ...

  **************************************


The program then outputs a series of files:

3d_sample5.struc_tree.list is a list of the tRMSD tree associated with every position columns
3d_sample5.struc_tree.html is a colored output showing columns accordingg to their support to the tree (red: high, blue: low)
3d_sample5.struc_tree.consensus_output is a schematic representation of the results (it's better to use a tree viewer)
3d_sample5.struc_tree.consensus is the final consensus structural tree 


.. note:: Please cite: Magis, C., Stricher, F., van der Sloot, A.M., Serrano, L., Notredame, C. T-RMSD: a fine-grained, structure based classification method and its application to the functional characterization of TNF receptors. J. Mol. Biol., 400(3):605-617 (2010), PMID:20471393 and/or Magis, C., van der Sloot, A.M., Serrano, L., Notredame, C. An improved understanding of TNFL/TNFR interactions using structure-based classifications. Trends Biochem. Sci., 37(9):353-363 (2012), PMID:22789664


TCS
===

to be done...


STRIKE
======

to be done...


SARA-Coffee
===========

to be done...

