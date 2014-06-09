Transitive consistency score 
================

TCS is an alignment evaluation score that makes it possible to identify in an MSA the most correct positions. 
It has been shown that these positions are the most likely to be structuraly correct and also the most informative when estimating phylogenetic trees. The TCS evaluation and filtering procedure is implemented in the T-Coffee package and can be used to evaluate and filter any third party multiple sequence alignment (including T-Coffee of course!)


Evaluate an existing MSA 
---------------


:: 

  t_coffee -infile test.prot.aln -evaluate -output score_ascii, aln, score_html

Outputs: 

* est.prot.score_ascii: displays the score of the MSA, the sequences and theresiues. This fille can be used to further filter your MSA with seq_reformat 
* test.prot.score_html: displays a colored version score of the MSA, the sequences and the resiues. 


Filter unreliable MSA positions
---------------

:: 

  t_coffee -infile test.prot.aln -evaluate -output tcs_residue_filter3, tcs_column_filter3, tcs_residue_lower4

Outputs: 

* test.prot.tcs_residue_filter3: All residues with a TCS score lower than 3 are filtered out 
* test.prot.tcs_column_filter3: All columns with a TCS score lower than 3 are filtered out 
* test.prot.tcs_residue_lower4: All residues with a TCS score lower than 3 are lower cased
  
Note that all these outpout functions are also compatible with the default T-Coffee when computing an alignment::

  t_coffee -seq test.prot.fa -output tcs_residue_filter3, tcs_column_filter3, tcs_residue_lower4

or with seq_reformat using a T-Coffee .score_ascii file:: 

  t_coffee -other_pg seq_reformat -in test.prot.aln -struc_in test.prot.score_ascii -struc_in_f number_aln -output tcs_residue_filter3
  


Weight MSA positions for improved trees
---------------


:: 

  t_coffee -infile test.prot.aln -evaluate -output tcs_weighted, tcs_replicate_100

Outputs: 

* test.prot.tcs_weighted: All columns are duplicated according to their TCS score 
* test.prot.tcs_replicate_100: Contains 100 replicates in phylip format with each column drawn with a probability corresponding to its TCS score 


Note that all these outpout functions are also compatible with the default T-Coffee when computing an alignment::

  t_coffee -seq test.prot.fa -output tcs_weighted, tcs_replicate_100

or with seq_reformat using a T-Coffee .score_ascii file:: 

  t_coffee -other_pg seq_reformat -in test.prot.aln -struc_in test.prot.score_ascii -struc_in_f number_aln -output tcs_weighted



Work with coding DNA
---------------

When working with DNA, it is advisable to first align the sequences at the protein level and later thread back the DNA onto your aligned proteins.
The filtering must be done in two steps, as shown below. Note that your DNA and protein sequences must have the same name:: 

  t_coffee -infile test.prot.aln -evaluate -output score_ascii

This first step produces the TCS evaluation file test.prot.score_ascii::
 
  t_coffee -other_pg seq_reformat -in test.prot.aln -in2 test.dna.fa -struc_in test.prot.score_ascii -struc_in_f number_aln -output tcs_replicate_100 -out test.dna.replicates
  
test.dna.replicates: 100 DNA replicates with positions selected according to their AA TCS score::

  t_coffee -other_pg seq_reformat -in test.prot.aln -in2 test.dna.fa -struc_in test.prot.score_ascii -struc_in_f number_aln -output tcs_column_filter5 -out test.dna.filter  

test.dna.filtered: DNA positions filtered according to their TCS column score



Using different TCS libraries
---------------

It is possible to change the way TCS reliability is estimated. 
This can be done by building different T-Coffee libraries. proba_pair is the default mode of T-Coffee that runs a pair-HMM to populate the library with residue pairs having the best posterior probabilities.
The following instructions will do this:: 

  t_coffee -infile test.aln -evaluate -method proba_pair -output score_ascii, aln, score_html

This mode runs a series of fast multiple aligners. It is very fast and used by ENSEMBL Compara:: 

  t_coffee -infile test.aln -evaluate -method mafft_msa,kalign_msa,muscle_msa -output score_ascii, aln, score_html

This mode runs the orginal default T-Coffee that was combining local and global alignments:: 

  t_coffee -infile test.aln -evaluate -method clustalw_pair,lalign_id_pair -output score_ascii, aln, score_html


Summary of the various output flags
---------------

============================ 	================
Flags        					Description
============================ 	================
-output=score_ascii	    		outputs a TCS evaluation file
-output=score_html				contains ascii format in html format
-output=score_pdf				will transfer score_html into pdf format
-output=sp_ascii				is a format reporting the TCS score of every aligned pair in the target MSA
-output=tcs_residue_filter_N	Removes all residues with a TCS score lower than N
-output=tcs_columns_filter_N	Removes all columns with a TCS score lower than N
-output=tcs_weighted	 		Duplicates all columns according to their TCS scoreN, The output is in Phylips
-output=tcs_replicate_N	 		Generates N phylips replicates, with columns drawn according to their TCS score
============================ 	================


Reference
---------------

A pre-print is available here and the final publication can be accessed from pubmed.

A powerpoint is available here.

Please Cite:TCS, A new multiple sequence alignment reliability measure to estimate alignment accuracy and improve phylogenetic tree reconstruction. Chang JM, Tommaso PD, Notredame C., MBE (2014), doi: 10.1093/molbev/msu117pubmed


Links
---------------

#. structural validation

* BAliBASE 3
* PREFAB 4

#. phylogenetic validation

* yeasts from Wong et al. Science, 2008
* subset gene list : at least one aligner yields a phylogeny topology identical to the canonical yeast ToL
* tips16 from Gblocks
* tips32, tips64 from trimAl