###########
Quick Start
###########


.. warning:: All the files mentionned here (e.g. sample_seq) can be found in the example directory of the distribution.


*******************
Basic Command Lines
*******************

Proteins sequences
==================
::

  Mode Command
  ============================================================================

  ----------------------------------------------------------------------------

  Regular              *t_coffee sample_seq1.fasta*

                       use the output.html to estimate the MSA accuracy
  ----------------------------------------------------------------------------
 
  Very Fast            t_coffee sample_seq1.fasta -mode quickaln

                       lower -ndiag if the sequences are very similar

  ----------------------------------------------------------------------------

  Very Accurate        t_coffee sample_aln1.fasta -mode accurate

                       slow, combines structures, sequences and profiles

  ----------------------------------------------------------------------------

  M-Coffee             t_coffee sample_aln1.fasta -mode mcoffee

                       combines most of the existing MSA packages

  ----------------------------------------------------------------------------

  3D-Coffee            t_coffee sample_aln1.fasta -mode 3dcoffee

                       uses the structure of your sequences named with PDBID

  ----------------------------------------------------------------------------

  Expresso             t_coffee sample_aln1.fasta -mode expresso

                       finds structures homologous to your sequences

  ----------------------------------------------------------------------------

  PSI-Coffee           t_coffee sample_aln1.fasta -mode psicoffee

                       enriches your sequence with profile information

  ----------------------------------------------------------------------------


DNA sequences
=============
::


  Mode Command
  ============================================================================

  ----------------------------------------------------------------------------
  R-Coffee             t_coffee three_cdna.fasta -mode cdna



RNA sequences
=============
::

  Mode Command
  ============================================================================

  ----------------------------------------------------------------------------
  R-Coffee             t_coffee sample_rnaseq1.fasta -mode rcoffee

                       use the predicted secondary structure of your sequences

  ----------------------------------------------------------------------------

  RM-Coffee            t_coffee sample_rnaseq1.fasta -mode rmcoffee

                       use M-Coffee + secondary structure prediction

  ----------------------------------------------------------------------------

  R-Coffee Consan      t_coffee sample_rnaseq1.fasta -mode rcoffee_consan

                       use rcoffee to combine consan alignments (accurate/slow)

  ----------------------------------------------------------------------------
  

Memory Fix
==========
::

  memory  t_coffee sample_aln1.fasta -mode memory





**********************
Starting with T-Coffee
**********************
We only give you the very basics here. Please use the Tutorial for more detailed information on how to use our tools.


T-Coffee
========
Write your sequences in the same file (Swiss-prot, Fasta or Pir) and type.


::

  $$: t_coffee sample_seq1.fasta



This will output two files:


::

  sample_seq1.aln: your Multiple Sequence Alignment

  sample_seq1.dnd: The Guide tree (newick Format)



.. warning:: IMPORTANT: In theory nucleic acids should be automatically detected and the default methods should be adapted appropriately. However, sometimes this may fail, either because the sequences are too short or contain too many ambiguity codes. When this happens, you are advised to explicitly set the type of your sequences using the flag -type. NOTE: the -mode=dna is not needed or supported anymore

::

  $$: t_coffee sample_dnaseq1.fasta -type=dna



M-Coffee
========
M-Coffee is a Meta version of T-Coffee that makes it possible to combine the output of at least eight packages (Muscle, probcons, poa, dialignT, mafft, clustalw, PCMA and T-Coffee).


If all these packages are already installed on your machine. You must:


1) Set the following environment variables:


::

   export POA_DIR=[absolute path of the POA installation dir]

   export DIALIGNT_DIR=[Absolute path of the DIALIGN-T/conf



2) Write your sequences in a file and run T-Coffee using this file (Swiss-prot, Fasta or Pir) with:


::

  $$: t_coffee sample_seq1.fasta -mode mcoffee



If the program starts complaining one package or the other is missing, this means you will have to go the hard way and install all these packages yourself... Proceed to the M-Coffee section for more detailed instructions.


Expresso
========
If you have installed the EBI wublast.pl client, Expresso will BLAST your sequences against the PDB database, identify the best targets and use these to align your proteins using the following commandline:


::

  $$: t_coffee sample_seq1.fasta -mode expresso



If you did not manage to install all the required structural packages for Expresso you can still run eEpresso by selecting yourself the structural packages you want to use. For instance, if you'd rather use TM-Align than sap, try:



::

  $$: t_coffee sample_seq1.fasta -template_file expresso -method TMalign_pair



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



iRMSD and APDB
==============
All you need is a file containing the alignment of sequences with a known structure. These sequences must be named according to their PDB ID, followed by the chain index ( 1aabA for instance). All the sequences do not need to have a known structure, but at least two is required. Given the alignment, use the following command:


::

  $$: t_coffee -other_pg irmsd -aln 3d_sample4.aln



tRMSD
=====
tRMSD is a structure based clustering method using the iRMSD to drive the clustering. The T-RMSD supports all the parameters supported by iRMSD or APDB.


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


MOCCA
=====
Write your sequences in the same file (Swiss-prot, Fasta or Pir) and type:


::

  $$: t_coffee -other_pg mocca sample_seq1.fasta



This command output one files (<your sequences>.mocca_lib) and starts an interactive menu.

