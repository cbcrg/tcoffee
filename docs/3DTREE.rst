**********
3DTree
**********

*Tree estimation procedure based on the comparison of internal distances*

3DTree makes it possible to estimate a tree using either contact conservation or differences in internal distances as a measure of similarity bewtween protein or RNA sequences. The trees thus estimated can be bootsrapped or further analyzed like regular phylogenetic trees. 3DTree also makes it possible to estimate the local support of any structural alignment (i.e. each individual column) for either a full tree or any pre-defined sub-group contained within the dataset. 

Generating a Tree based on distances
=====================================

This option makes it possible to estimate a tree while taking into account the variation of intra-molecular distances within the considered sequences. The following call will generate a 100 replicate nj trees using the difference of distances between pairs of aligned residues, at a maximum cut-off of 15A. Columns with less than 50% residues are ignored


Input:

* ``aln``: Multiple Sequence Alignment in FASTA, MSA or MSF
* ``template``: FASTA name list with templates: ``>name _P_ template``

:: 

  t_coffee -other_pg seq_reformat -in <aln> -in2 <template> -action +tree replicates 100  +evaluate3D distances +tree2bs first -output newick -out tree.dnd


Outputs: 

* ``tree.dnd``: Tree in newick format with bootstrap support   

It is possible to control default parameters using the following extended command line

::

  t_coffee -other_pg seq_reformat -in <aln> -in2 <template> -action +tree replicates 100 gap 0.5 mode nj  +evaluate3D distances 15 +tree2bs first -output newick -out tree.dnd

.. warning: sequences without 3D structure will be excluded from the analysis and from the final output


Generating a Tree based on contact conservation
================================================

This option makes it possible to estimate a tree while taking into account the variation of contact conservation within the considered sequences. This call will generate a 100 replicate nj trees using as a distance metrics the fraction of contacts conserved between pairs of aligned residues, at a maximum cut-off of 1.2 A between VdW radius and ignoring the 3 closests neighbors. Columns with less than 50% residues are ignored. For sequences without 3D information, the strike contact potential is used instead (Watson and crick base pairing propensity for RNA).

:: 

  t_coffee -other_pg seq_reformat -in <seq.aln> -in2 <seq.template> -action +tree replicates 100  +evaluate3D contacts +tree2bs first -output newick -out tree.dnd


Outputs: 

* ``tree.dnd``: Tree in newick format  

It is possible to control default parameters using the following extended command line:

::

  seq_reformat -in <aln> -in2 <template> -action +tree replicates 100 gap 0.5 mode nj  +evaluate3D contacts 1.2 3 +tree2bs first -output newick -out tree.dnd

.. warning: the procedure requires at least 1 sequence with a known 3D structure or with contact information.



Visulizing 3D Conservation
================================================

This same procedure can be used to visualize either intra-molecular distance conservation or contact conservation

::

  seq_reformat -in CRD.aln -in2 CRD.template -action +evaluate3D distances -output score_html 
  seq_reformat -in CRD.aln -in2 CRD.template -action +evaluate3D distances -output score_ascii
  seq_reformat -in CRD.aln -in2 CRD.template -action +evaluate3D distances -output score_raw

Outputs:
* ``score_raw``: Tabulated dump of the numerical values associated with every residue, every sequence and every column of the considered alignment.

Identification of positions 
=============================

If you have a well defined sub-group of sequences (i.e. domains having the same function, same specificty, etc...), it is possible to estimate which columns yield the best support using the following command,

Input:
* ``group.fasta``: A Fasta formatted list of the sequences that form the group whose support you want to analyze

::

 seq_reformat -in <seq.aln> -in2 <seq.template> -action +tree replicates columns  +evaluate3D  distances +evaluateTree <group.fasta> -output score_html -out <aln.html>

Output
* ``aln.score_html`` Colored version of your MSA indicating the sequences that best contribute to your clustering.


Evaluating Clustering capacities
=================================

If you want to check the capacity of an algorithm to bring related sequences within mono-phyletic groups, you should name your sequences according to the group they belong to (XXXX_1, YYYYY_1, ZZZZ_2, KKKK_2, for members of _1 and _2, etc) and use the following evaluation procedure. The output will be the number of monophyletic groups containing sequences belonging to the same group:

The tree can be pre-computed
:: 

  seq_reformat -in <tree> +tree2collapse groups 4 +print nseq -output no

Or it can be computed on the fly
:: 

  seq_reformat -in <aln> -in2 <template> -action +tree replicates 100  +evaluate3D  distances 15 +tree2bs first +tree2collapse groups 4 +print nseq -output no

