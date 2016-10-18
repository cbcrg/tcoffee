###########
Quick Start
###########


.. warning:: This part is designed to have an quick overview of T-Coffee alignment procedures. All the files mentioned here can be found in the example directory of the distribution. If you use a particular mode of T-Coffee, please use the specified citation for publication.


***************************
Basic Command Lines (-mode)
***************************

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
  Default              t_coffee sample_seq1.fasta
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
  Default              t_coffee sample_seq1.fasta
                       use the output.html to visualize the MSA accuracy
                       
  Citation: Notredame et al., JMB (2000)                      PMID:10964570  
  ----------------------------------------------------------------------------
  Structural 2D        t_coffee sample_rnaseq1.fasta -mode rcoffee
  (R-Coffee)           use the predicted secondary structure
  
  Citation: Wilm et al., Nucleic Acids Res. (2008)            PMID:18420654
  ----------------------------------------------------------------------------
  Structural 3D        t_coffee sample_rnaseq1.fasta -mode rcoffee_consan
  (R-Coffee Consan)    use R-Coffee to combine consan structural alignments 
  
  Citation: Wilm et al., Nucleic Acids Res. (2008)            PMID:18420654
  ----------------------------------------------------------------------------
  Accurate             t_coffee sample_rnaseq1.fasta -mode rmcoffee
  (RM-Coffee)          use M-Coffee + secondary structure prediction
                       
  Citation: Wilm et al., Nucleic Acids Res. (2008)            PMID:18420654
  ----------------------------------------------------------------------------

  

Memory Fix
==========
::

  memory  t_coffee sample_aln1.fasta -mode memory

  "to be defined"


**********************
Starting With T-Coffee
**********************

.. warning:: We only give you the very basics here, please go to the T-Coffee manual for a more detailed description and available options for the different tools. You can also try the T-Coffee tutorial for a practical training on T-Coffee alignment and other functions using applied examples on published research data.


T-Coffee
========
Write your sequences (protein, DNA or RNA) in the same file in one of the following format: Swiss-prot, Fasta or Pir; then run T-Coffee with the following command:


::

  $$: t_coffee sample_seq1.fasta



This will output three files:


::

  sample_seq1.aln  : your multiple sequence alignment (Clustal format by default)

  sample_seq1.dnd  : the guide tree (Newick format)
  
  sample_seq1.html : the color coded MSA according to T-Coffee consistency color scheme (html)


In principle, the type of the sequences should be automatically detected and the default methods should be adapted accordingly. However sometimes this may fail, either because the sequences are too short or contain too many ambiguity codes. When this happens, you are advised to explicitly set the type of your sequences using the flag -type. (Note: the flag -mode=dna is not needed or supported anymore).

::

  $$: t_coffee sample_dnaseq1.fasta -type=dna


Citation: Notredame et al., JMB (2000), PMID:10964570


M-Coffee
========
M-Coffee is a meta version of T-Coffee that makes it possible to combine the output of a combination of eight packages (Muscle, probcons, poa, dialignT, mafft, clustalw, PCMA and T-Coffee).


If all these packages are already installed on your machine. You must:


1) Set the following environment variables:


::

   export POA_DIR=[absolute path of the POA installation dir]

   export DIALIGNT_DIR=[Absolute path of the DIALIGN-T/conf



2) Write your sequences in a file and run T-Coffee using this file (format: Swiss-prot, Fasta or Pir) with:


::

  $$: t_coffee sample_seq1.fasta -mode mcoffee



If the program starts complaining one package or the other is missing, this means you will have to go the hard way and install all these packages yourself... Proceed to the M-Coffee section for more detailed instructions.


Citation: Wallace et al., Nucleic Acids Res. (2006), PMID:16556910


Expresso
========
If you have installed the EBI wublast.pl client, Expresso will BLAST your sequences against the PDB database, identify the best targets and use them to align your proteins using the following command:


::

  $$: t_coffee sample_seq1.fasta -mode expresso



If you did not manage to install all the required structural packages for Expresso you can still run eEpresso by selecting yourself the structural packages you want to use. For instance, if you'd rather use TM-Align than sap, try:



::

  $$: t_coffee sample_seq1.fasta -template_file expresso -method TMalign_pair


Citation: Armougom et al. Nucleic Acids Res. (2006), PMID:16845081


R-Coffee
========
R-Coffee can be used to align RNA sequences, using their RNApfold predicted secondary structures. The best results are obtained by using the consan pairwise method. If you have consan installed, run:


::

  $$: t_coffee sample_rnaseq1.fasta -special_mode rcoffee_consan



This will only work if your sequences are short enough (less than 200 nucleotides). A good alternative is the rmcoffee mode that will run Muscle, Probcons4RNA and Mafft and then use the secondary structures predicted by RNApfold:


::

  $$: t_coffee sample_rnaseq1.fasta -mode rmcoffee



If you want to select yourself which methods should be combined by R-Coffee, run:


::

  $$: t_coffee sample_rnaseq1.fasta -mode rcoffee -method lalign_id_pair slow_pair


Citation: Wilm et al., Nucleic Acids Res. (2008), PMID:18420654


iRMSD and APDB
==============
All you need is a file containing the alignment of sequences with a known structure. These sequences must be named according to their PDB ID, followed by the chain index (1aabA for instance). All the sequences do not need to have a known structure, but at least two is required. Given the alignment, use the following command:


::

  $$: t_coffee -other_pg irmsd -aln 3d_sample4.aln


Citation: Armougom et al., Bioinformatics (2006), PMID:16873492


T-RMSD
=====
T-RMSD is a structure based clustering method using the iRMSD to drive the structural clustering of your sequences with an available structure. The T-RMSD supports all the parameters supported by iRMSD or APDB. To run T-RMSD, type:


::

  $$: t_coffee -other_pg trmsd -aln 3d_sample5.aln -template_file 3d_sample5.template_list


3d_sample5.aln is a multiple alignment in which each sequence has a known structure. The file 3d_sample5.template_list is a fasta like file declaring the structure associated with each sequence, in the form:


::

  > <seq_name> _P_ <PDB structure file or name>

  ******* 3d_sample5.template_list ********

  >2UWI-3A _P_ 2UWI-3.pdb

  >2UWI-2A _P_ 2UWI-2.pdb

  >2UWI-1A _P_ 2UWI-1.pdb

  >2HEY-4R _P_ 2HEY-4.pdb

  ...

  **************************************


The program then outputs a series of files:

3d_sample5.struc_tree.list is a list of the tRMSD tree associated with every position columns
3d_sample5.struc_tree.html is a colored output showing columns accordingg to their support to the tree (red: high, blue: low)
3d_sample5.struc_tree.consensus_output is a schematic representation of the results (it's better to use a tree viewer)
3d_sample5.struc_tree.consensus is the final consensus structural tree 


Citation: Magis et al., JMB (2010), PMID:20471393 and/or Magis et al., Trends Biochem. Sci. (2012), PMID:22789664


MOCCA
=====
MOCCA is a specific tool in T-Coffee designed to deal with highly divergent protein repeats.  Write your sequences in the same file (format: Swiss-prot, Fasta or Pir) and type:


::

  $$: t_coffee -other_pg mocca sample_seq1.fasta


This command output one files (<your sequences>.mocca_lib) and starts an interactive menu.


Citation: Notredame, Bioinformatis (2001), PMID:11301309

