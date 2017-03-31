###
FAQ
###

*******************************
Abnormal Terminations & Results
*******************************

Q: The program keeps crashing when I give my sequences
------------------------------------------------------
A-1: This may be a format problem. Try to reformat your sequences, we recommend the FASTA format. If the problem persists, contact us.

A-2: Your sequences may not be recognized for what they really are. Normally T-Coffee recognizes the type of your sequences automatically, but if it fails specify the sequence format with **-type** and the type of sequence **DNA, RNA or PROTEIN**.

A-3: Costly computation or data gathered over the net is stored by T-Coffee in a cache directory. Sometimes, some of these files can be corrupted and cause an abnormal termination. You can either empty the cache ( ~/.t_coffee/cache/), request T-Coffee to run without using it (command 1), or update the files corresponding to your data (command 2).


::

  Command 1: No cache
  $$: t_coffee -pdb=struc1.pdb,struc2.pdb,struc3.pdb -method sap_pair -cache=no

  Command 2: Update cache
  $$: t_coffee -pdb=struc1.pdb,struc2.pdb,struc3.pdb -method sap_pair -cache=update


Q: The default alignment is not good enough
--------------------------------------------
A: see next question

Q: The alignment contains obvious mistakes
------------------------------------------
A: This happens with most multiple alignment procedures. However, wrong alignments are sometimes caused by bugs or an implementation mistake. Please report the most unexpected results to the authors.

Q: The program is crashing
--------------------------
A: If you get the message **FAILED TO ALLOCATE REQUIRED MEMORY** refer to the next question. Otherwise, if the program crashes check whether you are using the right syntax. If the problem persists, contact us.

Q: I am running out of memory
-----------------------------
A: You can use a more accurate, slower and less memory hungry dynamic programming mode called myers_miller_pair_wise (command 1). Note that this mode will be much less time efficient than the default, although it may be slightly more accurate. In practice the parameterization associate with special mode turns off every memory expensive heuristic within T-Coffee (command 2). If you keep running out of memory, you may also want to reduce the number of sequences using **-maxnseq**.

::

  Command 1:
  $$: t_coffee sample_seq1.fasta -mode low_memory

  Command 2:
  $$: t_coffee sample_seq1.fasta -method=slow_pair,lalign_id_pair -distance_matrix_mode \
      =idscore -dp_mode=myers_miller_pair_wise


********************
Input/Output control
********************
Q: How many sequences can T_Coffee handle? [Deprecated]
------------------------------------------
A: T-Coffee is limited to a maximum of 50 sequences. Above this number, the program automatically switches to a heuristic mode, named DPA, where DPA stands for Double Progressive Alignment. DPA is still in development and the version currently shipped with T-Coffee is only a beta version.

Q: Can I prevent the output of all the warnings?
-----------------------------------------------
A: Yes, by setting **-no_warning**.

Q: How many ways to pass parameters to t_coffee?
------------------------------------------------
A: Refer to the **T-Coffee Technical Documentation**, Section **T-Coffee Parameters & Flags**. 

Q: How can I change the default output format?
----------------------------------------------
A: See the **-output** option in the **T-Coffee Main Documentation**; nearly all common output formats are recognized by T-Coffee.

Q: My sequences are slightly different between all the alignments.
------------------------------------------------------------------
A: It does not matter. T-Coffee will reconstruct a set of sequences that incorporates all the residues potentially missing in some of the sequences (See flag **-in**).

Q: Is it possible to pipe stuff out of T-Coffee?
------------------------------------------------
A: Specify stderr or stdout as output filename, the output will be redirected accordingly. For instance this instruction will output the tree (in Newick format) and the alignment to stdout.

::

  $$: t_coffee sample_seq1.fasta -outfile=stdout -out_lib=stdout


Q: Is it possible to pipe stuff into T_Coffee?
----------------------------------------------
A: If as a file name you specify stdin, the content of this file will be expected throught pipe:

::

  $$: cat sample_seq1.fasta | t_coffee -infile=stdin

  is equivalent to:

  $$: t_coffee sample_seq1.fasta


If you do not give any argument to T-Coffee, they will be expected to come from pipe:


::

  $$: cat sample_param_file.param | t_coffee -parameters=stdin

  or 
  
  $$: echo -seq=sample_seq1.fasta -method=clustalw_pair | t_coffee -parameters=stdin


Q: Can I read my parameters from a file?
----------------------------------------
A: See the **T-Coffee Technical Documentation**.



Q: I want to decide myself on the name of the output files !!!
--------------------------------------------------------------
A: Use the **-run_name** flag:

::

  $$: t_coffee sample_seq1.fasta -run_name=luke_skywalker


Q: I want to use the sequences in an alignment file
---------------------------------------------------
A: Simply fed your alignment any way you like, but do not forget to append the prefix S for sequence:

::

  $$: t_coffee sample_aln1.aln -in proba_pair

  $$: t_coffee -seq=sample_aln1.aln -method=slow_pair,lalign_id_pair -outfile=outaln


This means that the gaps will be reset and that the alignment you provide will not be considered as an alignment, but as a set of sequences.

Q: I only want to produce a library
-----------------------------------
A: use the **-lib_only** flag; but note that this supersedes the use of the **-convert** flag. Its main advantage is to restrict computation time to the actual library computation.

::

  $$: t_coffee sample_seq1.fasta -out_lib=sample_lib1.tc_lib -lib_only


Q: I want to turn an alignment into a library
---------------------------------------------
A: use the **-lib_only** flag (command 1). It is also possible to control the weight associated with this alignment with the flag **-weight** (command 2).

::

  Command 1:
  $$: t_coffee -in=Asample_aln1.aln -out_lib=sample_lib1.tc_lib -lib_only

  Command 2: 
  $$: t_coffee -aln=sample_aln1.aln -out_lib=sample_lib1.tc_lib -lib_only -weight=1000


Q: I want to concatenate two libraries
--------------------------------------
A: You cannot concatenate these files on their own. You will have to use T-Coffee assuming you want to combine for instance ``tc_lib1.tc_lib`` and ``tc_lib2.tc_lib``:

::

  $$: t_coffee -lib=sample_lib1.tc_lib,sample_lib2.tc_lib -lib_only -out_lib=sample_lib3.tc_lib


Q: What happens to the gaps when an alignment is fed to T-Coffee?
-----------------------------------------------------------------
A: An alignment is **ALWAYS** considered as a library **AND** a set of sequences. If you want your alignment to be considered as a library only, use the S identifier; it will be seen as a sequence file, even if it has an alignment format (gaps will be removed).

::

  $$: t_coffee Ssample_aln1.aln -outfile=outaln


Q: I cannot print the html graphic display!!!
---------------------------------------------
A: This is a problem that has to do with your browser. Instead of requesting the score_html output, request the score_ps output that can be read using ghostview:

::

  Postscript
  $$: t_coffee sample_seq1.fasta -output=score_ps
   
  PDF only if you have ps2pdf installed
  $*: t_coffee sample_seq1.fasta -output=score_pdf


Q: I want to output an html file and a regular file
---------------------------------------------------
A: See the next question.


Q: I would like to output more than one alignment format at the same time
-------------------------------------------------------------------------
A: The flag **-output** accepts more than one parameter. For instance this will output four alignment files in the corresponding formats. Alignments' names will have the format name as an extension.

::

  $$: t_coffee sample_seq1.fasta -output=clustalw,html,score_ps,msf


.. note:: Note: you need to have the converter ps2pdf installed on your system (standard under Linux and Cygwin). The latest versions of Internet Explorer and Netscape now allow the user to print the html display. Do not forget to request background printing.

*********************
Alignment computation
*********************
Q: Is T-Coffee the best? Why not using MUSCLE, MAFFT, or ProbCons???
--------------------------------------------------------------------
A: All these packages are good packages and they sometimes outperform T-Coffee. They also claim to outperform one another... If you have them installed locally, you can have T-Coffee to generate a consensus alignment:

::

  $$: t_coffee sample_seq1.fasta -method muscle_msa,probcons_msa,mafft_msa,lalign_id_pair,slow_pair


Q: Can T_Coffee align nucleic acids ???
---------------------------------------
A: Normally it can, but check in the log that the program recognises the right type. If this fails, you will need to manually set the type using **-type dna**


Q: I do not want to compute the alignment
-----------------------------------------
A: use the **-convert** flag. This command will read the .aln file and turn it into an .msf alignment.

::

  $$: t_coffee sample_aln1.aln -convert -output=gcg


Q: I would like to force some residues to be aligned
----------------------------------------------------
If you want to brutally force some residues to be aligned, you may use as a post processing, the **+force_aln** function of **seq_reformat**. You can either specify single (command 1) or multiple constraints using a TC_LIB_FORMAT_02 file (command 2). When giving more than one constraint, these will be applied one after the other in the order they are provided. This greedy procedure means that the Nth constraint may disrupt the (N-1)th previously imposed constraint, hence the importance of forcing the constraints in the right order, with the most important coming last. We do not recommend imposing hard constraints on an alignment, and it is much more advisable to use the soft constraints provided by standard T-Coffee libraries (cf. **T-Coffee Technical Documentation**, subsection **Creating your own T-Coffee libraries**).

::

  Command 1: single constraint
  $$: t_coffee -other_pg seq_reformat -in sample_aln3.aln -action +force_aln seq1 5 seq2 6
  
  Command 2: multiple constraints
  $$: t_coffee -other_pg seq_reformat -in sample_aln3.aln -action +force_aln sample_lib3.tc_lib02


The TC_LIB_FORMAT_02 is still experimental and unsupported. It can only be used in the context of the force_aln function described here. The tc_lib02 format is as follow:

::

  *TC_LIB_FORMAT_02
  SeqX resY ResY_index  SeqZ ResZ ResZ_index



Q: I would like to use structural alignments
--------------------------------------------
Refer to the **T-Coffee Main Documentation and/or T-Coffee Technical Documentation**.


Q: I want to build my own libraries
-----------------------------------
A: Turn your alignment into a library, forcing the residues to have a very good weight, using structure:

::

  $$: t_coffee -aln=sample_seq1.aln -weight=1000 -out_lib=sample_seq1.tc_lib -lib_only


The value 1000 is simply a high value that should make it more likely for the substitution found in your alignment to reoccur in the final alignment. This will produce the library sample_aln1.tc_lib that you can later use when aligning all the sequences:

::

  $$: t_coffee -seq=sample_seq1.fasta -lib=sample_seq1.tc_lib -outfile sample_seq1.aln


If you only want some of these residues to be aligned, or want to give them individual weights, you will have to edit the library file yourself or use the -force_aln option (cf FAQ: I would like to force some residues to be aligned). A value of N*N * 1000 (N being the number of sequences) usually ensure the respect of a constraint.


Q: I want to use my own tree
----------------------------
A: Use the **-usetree=<your own tree>** flag:

::

  $$: t_coffee sample_seq1.fasta -usetree=sample_seq1_tree_nj.nwk


Q: I want to align coding DNA
-----------------------------
A: Use the **fasta_cdna_pair** method that compares two cDNA using the best reading frame and taking frameshifts into account. Notice that in the resulting alignments (command 1), all the gaps are of modulo3, except one small gap in the first line of sequence hmgl_trybr. This is a frameshift made on purpose. You can realign the same sequences while ignoring their coding potential and treating them like standard DNA (command 2).

::

  Command 1:
  $$: t_coffee three_cdna.fasta -method=cdna_fast_pair

  Command 2:
  $$: t_coffee three_cdna.fasta


.. warning:: This method has not yet been fully tested and is only provided 'as-is' with no warranty. Any feedback will be much appreciated.

Q: I do not want to use all the possible pairs when computing the library
-------------------------------------------------------------------------
See next question.

Q: I only want to use specific pairs to compute the library
-----------------------------------------------------------
A: Simply write in a file the list of sequence groups you want to use. Pairwise methods (slow_pair, proba_pair, <method>_pair...) will only be applied to list of pairs of sequences, while multiple methods (clustalw_msa, mafft_msa, <method_msa...) will be applied to any dataset having more than two sequences.

::

  $$: t_coffee sample_seq1.fasta -method=clustalw_pair,clustalw_msa -lib_list=sample_list1.lib_list

  Format of the list of libraries:
  ***************sample_list1.lib_list****
  2 hmgl_trybr hmgt_mouse
  2 hmgl_trybr hmgb_chite
  2 hmgl_trybr hmgl_wheat
  3 hmgl_trybr hmgl_wheat hmgl_mouse
  ***************sample_list1.lib_list****


Q: There are duplicates or quasi-duplicates in my set.
------------------------------------------------------
A: If you can remove them, this will make the program run faster, otherwise the T-Coffee scoring scheme should be able to avoid overweighting of overrepresented sequences.


*****************************
Using Structures and Profiles
*****************************
Q: Can I align sequences to a profile with T-Coffee?
----------------------------------------------------
A: Yes, you simply need to indicate that your alignment is a profile with the R tag:

::

  $$: t_coffee sample_seq1.fasta -profile=sample_aln2.aln -outfile chewbacca


Q: Can I align sequences two or more profiles? 
----------------------------------------------
A: Yes, you, simply tag your profiles with the letter R and the program will treat them like standard sequences:


::

  $$: t_coffee -profile=sample_aln1.aln,sample_aln2.aln -outfile han_solo



Q: Can I align two profiles according to the structures they contain?
---------------------------------------------------------------------
A: Yes, as long as the structure sequences are named according to their PDB identifier:

::

  $$: t_coffee -profile=sample_profile1.aln,sample_profile2.aln -mode=3dcoffee


Q: T-Coffee becomes very slow when combining sequences and structures
---------------------------------------------------------------------
A: This is true. By default the structures are fetched through the net using RCSB. The problem arises when T-Coffee looks for the structure of sequences WITHOUT structures. One solution is to install the PDB database locally. In that case you will need to set two environment variables:

::

  Variables to set up:
  ##: setenv (or export) PDB_DIR='directory containing the pdb structures' 
  ##: setenv (or export) NO_REMOTE_PDB_DIR=1


Interestingly, the observation that sequences without structures are those that take the most time to be checked is a reminder of the strongest rational argument that I know of against torture: any innocent would require the maximum amount of torture to establish his/her innocence, which sounds...hummmm...strange.. Then again I was never struck by the efficiency of the Bush Jr administration.

Q: Can I use a local installation of PDB?
-----------------------------------------
A: Yes, T-Coffee supports three types of installations:

- *Ad hoc* installation where all your structures are in a directory under the form pdbid.pdb, pdbid.id.Z or pdbid.pdb.gz. In that case, all you need to do is set the environement variables correctly:

::

  Setting up variable
  ##: setenv (or export) PDB_DIR='directory containing the pdb structures' 
  ##: setenv (or export) NO_REMOTE_PDB_DIR=1


- Full standard PDB installation using the all section of PDB. In that case, you must set the variables to:

::

  Setting up variable
  ##: setenv (or export) PDB_DIR='<some absolute path>/data/structures/all/pdb/' 
  ##: setenv (or export) NO_REMOTE_PDB_DIR=1


- Reduced standard PDB installation using the divided section of pdb:

::

  Setting up the PDB:
  ##: setenv (or export) PDB_DIR='<some absolute path>/data/structures/divided/pdb/
  ##: setenv (or export) NO_REMOTE_PDB_DIR=1


If you need to do more clever things, you should know that all the PDB manipulation is made in T-Coffee by a perl script named **extract_from_pdb**. You can then edit the script to suit your needs; T-Coffee will use your edited version if it is in the current directory and issue a warning that it used a local version. If you make extensive modifications, I would appreciate you send me the corrected file so that I can incorporate it in the next distribution. By default, T-Coffee also requires two important PDB files declared using the two following variables. These variables do not need to be set if the considered files are in the cache directory (default behavior): 

::

  Found at: ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_entry_type.txt
  ##: export PDB_ENTRY_TYPE_FILE=<location of the file pdb_entry_type.txt>
 
  Found at: http://www.rcsb.org/pdb/rest/getUnreleased  
  ##: export PDB_UNREALEASED_FILE=<location of the file unrealeased.xml>


.. warning:: Since the file ``unreleased.xml`` is not part of the PDB distribution, T-Coffee will make an attempt to obtain it even when using the **NO_REMOTE_PDB_DIR=1 mode**. You must therefore make sure that the file ``PDB_UNREALEASED_FILE`` is pointing to is read and write.


******************************
Improving/Evaluating Your MSAs
******************************
Q: How can I edit my alignment manually?
----------------------------------------
A: We recommend to use Jalview, a free program for MSA editing that you can find `here <http://www.jalview.org>`_.

Q: Have I improved or not my alignment?
---------------------------------------
A: Using structural information is the only way to establish whether you have improved or not your alignment. The CORE index can also give you some information. Refers to the **T-Coffee Main Documentation**, section **Evaluating Your Alignment**.

Q: How good is my alignment?
----------------------------
A: Refers to the **T-Coffee Main Documentation**, section **Evaluating Your Alignment**. Or just look at the color index ;-)

Q: What is that color index?
----------------------------
A: T-Coffee can provide you with a measure of consistency among all the methods used. An html file is produced by default each time you run an alignment. This html file is a colored version of your MSA that you can visualize with any common browser. As alternatives, you can use **score_ps** (postscript), **score_pdf** (pdf file) or **score_ascii** (text file). If you want more information about the CORE index represented by this color index, have a look at this `chapter <http://www.tcoffee.org/Publications/Pdf/core.pp.pdf>`_.

Q: Can I evaluate alignments NOT produced with T-Coffee?
--------------------------------------------------------
A: Yes !! You may have an alignment produced from any source you like. If you have no library available, the library will be computed on the fly but this can take some time depending on your sample size.

::

  With a library:
  $$: t_coffee -infile=sample_aln1.aln -lib=sample_aln1.tc_lib -special_mode=evaluate

  Without a library:
  $$: t_coffee -infile=sample_aln1.aln -evaluate -method proba_pair


Q: Can I compare two alignments?
--------------------------------
A: Yes. You can treat one of your alignments as a library and compare it with the second alignment. 

::

  $$: t_coffee -infile=sample_aln1_1.aln -aln=sample_aln1_2.aln -special_mode=evaluate


Q: I am aligning sequences with long regions of very good overlap
-----------------------------------------------------------------
A: Increase the ktuple size (up to 4-5 for DNA) and up to 3 for proteins. This will speed up the program. It can be very useful, especially when aligning ESTs.

::

  $$: t_coffee sample_seq1.fasta -ktuple=3



Q: Why is T-Coffee changing the names of my sequences!!!!
---------------------------------------------------------
A: If there is no duplicated name in your sequence set, T-Coffee handles names similarly to Clustalw. If your dataset contains sequences with identical names, these will automatically be renamed by adding an index (integer) to duplicated names even if there are more than 2. Also be careful, if there are spaces in your names, whatever comes after the space is not read.


.. danger:: The behaviour is undefined when this creates two sequence with a similar names.


*************
Release Notes
*************

.. Warning:: This log of modifications is not as thorough and accurate as it should be...but it's a beginning !
- 11.00+: extensive update of the documentation, examples; addition of the latest T-Coffee modes (PSI-Coffee, SARA-Coffee, Pro-Coffee, STRIKE, T-RMSD...); creation of an automated procedure for checking command lines from the documentation **doc2test.pl**.
- 9.86 New data structure for the primary library that results in highly improved running times for mcoffee and significantly decreased memory usage.
- 5.80 Novel assembly algorithm (linked_pair_wise) and the primary library is now made of probcons style pairwise alignments (proba_pair)
- 4.30 and upward: the FAQ has moved into a new tutorial document; **-in** can be replaced by the flags **-profile,-method,-aln,-seq,-pdb**.
- 4.02: **-mode=dna** is still available but not any more needed or supported. Use **-type=protein or dna** if you need to force things
- 3.28: corrected a bug that prevents short sequences from being correctly aligned
- Use of @ as a separator when specifying methods parameters
- The most notable modifications have to do with the structure of the input. From version 2.20, all files must be tagged to indicate their nature (A: alignment, S: Sequence, L: Library...). We are becoming stricter, but that's for your own good... Another important modification has to do with the flag -matrix: it now controls the matrix being used for the computation
 
