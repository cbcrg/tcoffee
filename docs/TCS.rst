Transitive Consistency Score 
=============================

TCS is an alignment evaluation score that makes it possible to identify the most correct positions in an MSA. 
It has been shown that these positions are the most likely to be structuraly correct and also the most informative when estimating phylogenetic trees. The TCS evaluation and filtering procedure is implemented in the T-Coffee package and can be used to evaluate and filter any third party multiple sequence alignment (including T-Coffee of course!)


Evaluate an existing MSA 
------------------------


:: 

  t_coffee -infile prot.aln -evaluate -output score_ascii, aln, score_html

Outputs: 

* `prot.score_ascii`  displays the score of the MSA, the sequences and the resiues. This file can be used to further filter your MSA with seq_reformat. 
* ``prot.score_html`` displays a colored version score of the MSA, the sequences and the residues. 

.. warning:: The color code in the score_html indicates the agreement between the library and the considered alignment. It is important to 
 understand that this score does not only depend on the input MSA, but it also depends on the library.

.. tip:: The TCS is most informative when used to identify low-scoring portions within an MSA. It is also worth noting that the TCS is not informative when aligning less than five sequences.
  
Filter unreliable MSA positions
-------------------------------

:: 

  t_coffee -infile prot.aln -evaluate -output tcs_residue_filter3, tcs_column_filter3, tcs_residue_lower4

Outputs: 

* `prot.tcs_residue_filter3`  All residues with a TCS score lower than 3 are filtered out 
* `prot.tcs_column_filter3`   All columns with a TCS score lower than 3 are filtered out 
* `prot.tcs_residue_lower4`   All residues with a TCS score lower than 3 are lower cased
  
Note that all these outpout functions are also compatible with the default T-Coffee when computing an alignment::

  t_coffee -seq prot.fa -output tcs_residue_filter3, tcs_column_filter3, tcs_residue_lower4

or with ``seq_reformat`` using a T-Coffee `.score_ascii` file:: 

  t_coffee -other_pg seq_reformat -in prot.aln -struc_in prot.score_ascii -struc_in_f number_aln -output tcs_residue_filter3
  


Weight MSA for improved trees
---------------------------------------


:: 

  t_coffee -infile prot.aln -evaluate -output tcs_weighted, tcs_replicate_100

Outputs: 

* `prot.tcs_weighted`       All columns are duplicated according to their TCS score 
* `prot.tcs_replicate_100`  Contains 100 replicates in phylip format with each column drawn with a probability corresponding to its TCS score 


Note that all these output functions are also compatible with the default T-Coffee when computing an alignment::

  t_coffee -seq prot.fa -output tcs_weighted, tcs_replicate_100

or with ``seq_reformat`` using a T-Coffee `.score_ascii` file:: 

  t_coffee -other_pg seq_reformat -in prot.aln -struc_in prot.score_ascii -struc_in_f number_aln -output tcs_weighted



Working with coding DNA
--------------------

When working with DNA, it is advisable to first align the sequences at the protein level and later thread back the DNA onto your aligned proteins.
The filtering must be done in two steps, as shown below. Note that your DNA and protein sequences must have the same name:: 

  t_coffee -infile prot.aln -evaluate -output score_ascii

This first step produces the TCS evaluation file `prot.score_ascii`::
 
  t_coffee -other_pg seq_reformat -in prot.aln -in2 dna.fa -struc_in prot.score_ascii -struc_in_f number_aln -output tcs_replicate_100 -out dna.replicates
  
`dna.replicates` 100 DNA replicates with positions selected according to their AA TCS score::

  t_coffee -other_pg seq_reformat -in prot.aln -in2 dna.fa -struc_in prot.score_ascii -struc_in_f number_aln -output tcs_column_filter5 -out dna.filter  

`dna.filtered` DNA positions filtered according to their TCS column score



Using different libraries
-----------------------------

It is possible to change the way TCS reliability is estimated. 
This can be done by building different T-Coffee libraries. `proba_pair` is the default mode of T-Coffee that runs a pair-HMM to populate the library with residue pairs having the best posterior probabilities.
The following instructions will do this:: 

  t_coffee -infile prot.aln -evaluate -method proba_pair -output score_ascii, aln, score_html

This mode runs a series of fast multiple aligners. It is very fast and used by `ENSEMBL Compara <http://www.ensembl.org/info/genome/compara/index.html>`_:: 

  t_coffee -infile prot.aln -evaluate -method mafft_msa,kalign_msa,muscle_msa -output score_ascii, aln, score_html

This mode runs the orginal default T-Coffee that was combining local and global alignments:: 

  t_coffee -infile prot.aln -evaluate -method clustalw_pair,lalign_id_pair -output score_ascii, aln, score_html


Summary of the output flags
-----------------------------------

============================ 	================
Flags        					Description
============================ 	================
-output=score_ascii	    		outputs a TCS evaluation file
-output=score_html				contains ascii format in html format
-output=score_pdf				will transfer score_html into pdf format
-output=sp_ascii				is a format reporting the TCS score of every aligned pair in the target MSA
-output=tcs_residue_filter_N	Removes all residues with a TCS score lower than `N`
-output=tcs_columns_filter_N	Removes all columns with a TCS score lower than `N`
-output=tcs_weighted	 		Outputs phylip format with duplicated columns according to their TCS score
-output=tcs_replicate_N	 		Generates `N` phylip replicates, with columns drawn according to their TCS score
============================ 	================


Supplement data
---------------

#. Documents

	* A `pre-print <http://www.tcoffee.org/Publications/Pdf/062-MBE-TCS.preprint.pdf>`_
	* The final publication can be accessed from `PubMed <http://www.ncbi.nlm.nih.gov/pubmed/24694831>`_
	* A `powerpoint <http://www.tcoffee.org/Projects/tcs/201403_TCS_MBE.pptx>`_

#. Example 

    * All files used in the above examples can be downloaded from `here <http://www.tcoffee.org/Projects/tcs/files.tar.gz>`_.

#. Structural validation

	* `BAliBASE 3 <http://lbgi.fr/balibase/>`_ [1]_
	* `PREFAB 4 <http://www.drive5.com/bench/>`_ [2]_

#. Phylogenetic validation

	* `yeasts <http://www.ncbi.nlm.nih.gov/pubmed/18218900>`_ [3]_
	* `subset gene list <http://www.tcoffee.org/Projects/tcs/geneList_YeastToL.txt>`_ : at least one aligner yields a phylogeny topology identical to the canonical yeast ToL
	* `tips16 <http://molevol.cmima.csic.es/castresana/data/syst_biol_2007_alignments.tar.Z>`_ from Gblocks [4]_
	* `tips32, tips64 <ftp://trimal.cgenomics.org/trimal/>`_ from trimAl [5]_


.. Note:: Please cite `Chang, J.-M., Tommaso, P. & Notredame, C. TCS: A New Multiple Sequence Alignment Reliability Measure to Estimate Alignment Accuracy and Improve Phylogenetic Tree Reconstruction. Molecular biology and evolution 31, 1625–37 (2014). doi: 10.1093/molbev/msu117 <http://www.ncbi.nlm.nih.gov/pubmed/24694831>`_.


Reference
---------------

.. [1] Thompson J, Koehl P, Ripp R, Poch O 2005. BAliBASE 3.0: Latest developments of the multiple sequence alignment benchmark. Proteins 61: 127-136.
.. [2] Edgar RC 2004. MUSCLE: multiple sequence alignment with high accuracy and high throughput. Nucleic Acids Res 32: 1792-1797.
.. [3] Wong KM, Suchard MA, Huelsenbeck JP 2008. Alignment uncertainty and genomic analysis. Science 319: 473-476.
.. [4] Talavera G, Castresana J 2007. Improvement of phylogenies after removing divergent and ambiguously aligned blocks from protein sequence alignments. Syst Biol 56: 564-577.
.. [5] Capella-Gutierrez S, Silla-Martinez JM, Gabaldon T 2009. trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. BIOINFORMATICS 25: 1972-1973.
