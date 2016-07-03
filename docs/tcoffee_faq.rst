################
FAQ for T-Coffee
################
.. tip:: All the files mentioned here (e.g. sample_seq) can be found in the example directory of the distribution.


***************************************
Abnormal Terminations and Wrong Results
***************************************

Q: The program keeps crashing when I give my sequences.
-------------------------------------------------------

A: This may be a format problem. Try to reformat your sequences using any utility (readseq...). We recommend the Fasta format. If the problem persists, contact us.


A: Your sequences may not be recognized for what they really are. Normally T-Coffee recognizes the type of your sequences automatically, but if it fails, use:


::

  $$: t_coffee sample_seq1.fasta -type=PROTEIN



A: Costly computation or data gathered over the net is stored by T-Coffee in a cache directory. Sometimes, some of these files can be corrupted and cause an abnormal termination. You can either empty the cache ( ~/.t_coffee/cache/) or request T-Coffee to run without using it:


::

  $$: t_coffee -pdb=struc1.pdb,struc2.pdb,struc3.pdb -method sap_pair -cache=no



If you do not want to empty your cache, you may also use -cache=update that will only update the files corresponding to your data


::

  $$: t_coffee -pdb=struc1.pdb,struc2.pdb,struc3.pdb -method sap_pair -cache=upd\
 ate



Q: The default alignment is not good enough.
--------------------------------------------

A: see next question


Q: The alignment contains obvious mistakes.
-------------------------------------------

A: This happens with most multiple alignment procedures. However, wrong alignments are sometimes caused by bugs or an implementation mistake. Please report the most unexpected results to the authors.


Q: The program is crashing.
---------------------------

A: If you get the message:

FAILED TO ALLOCATE REQUIRED MEMORY


See the next question.


If the program crashes for some other reason, please check whether you are using the right syntax and if the problem persists get in touch with the authors.


Q: I am running out of memory.
------------------------------

A: You can use a more accurate, slower and less memory hungry dynamic programming mode called myers_miller_pair_wise. Simply indicate the flag:


::

  $$: t_coffee sample_seq1.fasta -special_mode low_memory



Note that this mode will be much less time efficient than the default, although it may be slightly more accurate. In practice the parameterization associate with special mode turns off every memory expensive heuristic within T-Coffee.


::

  $$: t_coffee sample_seq1.fasta -method=slow_pair,lalign_id_pair -distance_mat\
 rix_mode=idscore -dp_mode=myers_miller_pair_wise



If you keep running out of memory, you may also want to lower -maxnseq, to ensure that t_coffee_dpa will be used.

********************
Input/Output control
********************

Q: How many sequences can T_Coffee handle?
------------------------------------------

A: T-Coffee is limited to a maximum of 50 sequences. Above this number, the program automatically switches to a heuristic mode, named DPA, where DPA stands for Double Progressive Alignment.


DPA is still in development and the version currently shipped with T-Coffee is only a beta version.


Q: Can I prevent the output of all the warnings?
------------------------------------------------

A: Yes, by setting -no_warning


Q: How many ways to pass parameters to t_coffee?
------------------------------------------------

A: See the section well behaved parameters


Q: How can I change the default output format?
----------------------------------------------

A: See the -output option, common output formats are:


::

  $$: t_coffee sample_seq1.fasta -output=msf,fasta_aln



Q: My sequences are slightly different between all the alignments.
------------------------------------------------------------------

A: It does not matter. T-Coffee will reconstruct a set of sequences that incorporates all the residues potentially missing in some of the sequences ( see flag -in).


Q: Is it possible to pipe stuff OUT of t_coffee?
------------------------------------------------

A: Specify stderr or stdout as output filename, the output will be redirected accordingly. For instance


::

  $$: t_coffee sample_seq1.fasta -outfile=stdout -out_lib=stdout



This instruction will output the tree (in new hampshire format) and the alignment to stdout.


Q: Is it possible to pipe stuff into T_Coffee?
----------------------------------------------

A: If as a file name, you specify stdin, the content of this file will be expected throught pipe:


::

  $$: cat sample_seq1.fasta | t_coffee -infile=stdin



will be equivalent to:


::

  $$: t_coffee sample_seq1.fasta



If you do not give any argument to t_coffee, they will be expected to come from pipe:


::

  $$: cat sample_param_file.param | t_coffee -parameters=stdin



For instance:


::

  $$: echo -seq=sample_seq1.fasta -method=clustalw_pair | t_coffee -parameters=s\
 tdin



Q: Can I read my parameters from a file?
----------------------------------------

A: See the well behaved parameters section.



Q: I want to decide myself on the name of the output files !!!
--------------------------------------------------------------

A: Use the -run_name flag:


::

  $$: t_coffee sample_seq1.fasta -run_name=guacamole



Q: I want to use the sequences in an alignment file.
----------------------------------------------------

A: Simply fed your alignment, any way you like, but do not forget to append the prefix S for sequence:


::

  $$: t_coffee Ssample_aln1.aln

  $$: t_coffee -infile=Ssample_aln1.aln

  $$: t_coffee -seq=sample_aln1.aln -method=slow_pair,lalign_id_pair -outfile=outaln



This means that the gaps will be reset and that the alignment you provide will not be considered as an alignment, but as a set of sequences.


Q: I only want to produce a library.
------------------------------------

A: use the -lib_only flag:


::

  $$: t_coffee sample_seq1.fasta -out_lib=sample_lib1.tc_lib -lib_only



Please, note that the previous usage supersedes the use of the -convert flag. Its main advantage is to restrict computation time to the actual library computation.


Q: I want to turn an alignment into a library.
----------------------------------------------

A: use the -lib_only flag:


::

  $$: t_coffee -in=Asample_aln1.aln -out_lib=sample_lib1.tc_lib -lib_only



It is also possible to control the weight associated with this alignment (see the -weight section):


::

  $$: t_coffee -aln=sample_aln1.aln -out_lib=sample_lib1.tc_lib -lib_only -weigh\
 t=1000



Q: I want to concatenate two libraries.
---------------------------------------

A: You cannot concatenate these files on their own. You will have to use t_coffee. Assume you want to combine tc_lib1.tc_lib and tc_lib2.tc_lib:


::

  $$: t_coffee -lib=sample_lib1.tc_lib,sample_lib2.tc_lib -lib_only -out_lib=sam\
 ple_lib3.tc_lib



Q: What happens to the gaps when an alignment is fed to T-Coffee?
-----------------------------------------------------------------

A: An alignment is ALWAYS considered as a library AND a set of sequences. If you want your alignment to be considered as a library only, use the S identifier:


::

  $$: t_coffee Ssample_aln1.aln -outfile=outaln



It will be seen as a sequence file, even if it has an alignment format (gaps will be removed).


Q: I cannot print the html graphic display!!!
---------------------------------------------

A: This is a problem that has to do with your browser. Instead of requesting the score_html output, request the score_ps output that can be read using ghostview:


::

  $$: t_coffee sample_seq1.fasta -output=score_ps



or


::

  $$: t_coffee sample_seq2.fasta -output=score_pdf



Q: I want to output an html file and a regular file.
----------------------------------------------------

A: see the next question


Q: I would like to output more than one alignment format at the same time.
--------------------------------------------------------------------------

A: The flag -output accepts more than one parameter. For instance:


::

  $$: t_coffee sample_seq1.fasta -output=clustalw,html,score_ps,msf



This will output founr alignment files in the corresponding formats. Alignments' names will have the format name as an extension.


.. note:: Note: you need to have the converter ps2pdf installed on your system (standard under Linux and cygwin). The latest versions of Internet Explorer and Netscape now allow the user to print the HTML display Do not forget to request Background printing.

Alignment computation
=====================
Q: Is T-Coffee the best? Why not using Muscle, or Mafft, or ProbCons???
-----------------------------------------------------------------------

A: All these packages are good packages and they sometimes outperform T-Coffee. They also claim to outperform one another... If you have them installed locally, you can have T-Coffee to generate a consensus alignment:


::

  $$: t_coffee sample_seq1.fasta -method muscle_msa,probcons_msa, mafft_msa, lal\
 ign_id_pair,slow_pair



Q: Can T_Coffee align nucleic acids ???
---------------------------------------

A: Normally it can, but check in the log that the program recognises the right type (In the INPUT SEQ section, Type: xxxx). If this fails, you will need to manually set the type:


::

  $$: t_coffee sample_dnaseq1.fasta -type dna



Q: I do not want to compute the alignment.
------------------------------------------

A: use the -convert flag:


::

  $$: t_coffee sample_aln1.aln -convert -output=gcg



This command will read the .aln file and turn it into an .msf alignment.


Q: I would like to force some residues to be aligned.
-----------------------------------------------------

If you want to brutally force some residues to be aligned, you may use as a post processing, the force_aln function of seq_reformat:


::

  $$: t_coffee -other_pg seq_reformat -in sample_aln4.aln -action +force_aln seq\
 1 10 seq2 15

  $$: t_coffee -other_pg seq_reformat -in sample_aln4.aln -action +force_aln sample_lib4.tc_lib02



sample_lib4.tc_lib02 is a T-Coffee library using the tc_lib02 format:


::

  *TC_LIB_FORMAT_02

  SeqX resY ResY_index  SeqZ ResZ ResZ_index



.. warning:: The TC_LIB_FORMAT_02 is still experimental and unsupported. It can only be used in the context of the force_aln function described here.

Given more than one constraint, these will be applied one after the other, in the order they are provided. This greedy procedure means that the Nth constraint may disrupt the (N-1)th previously imposed constraint, hence the importance of forcing the constraints in the right order, with the most important coming last.


We do not recommend imposing hard constraints on an alignment, and it is much more advisable to use the soft constraints provided by standard t_coffee libraries (cf. building your own libraries section)


Q: I would like to use structural alignments.
---------------------------------------------

See the section "Using structures in Multiple Sequence Alignments", or see the question I want to build my own libraries.


Q: I want to build my own libraries.
------------------------------------

A: Turn your alignment into a library, forcing the residues to have a very good weight, using structure:


::

  $$: t_coffee -aln=sample_seq1.aln -weight=1000 -out_lib=sample_seq1.tc_lib -li\
 b_only



The value 1000 is simply a high value that should make it more likely for the substitution found in your alignment to reoccur in the final alignment. This will produce the library sample_aln1.tc_lib that you can later use when aligning all the sequences:


::

  $$: t_coffee -seq=sample_seq1.fasta -lib=sample_seq1.tc_lib -outfile sample_se\
 q1.aln



If you only want some of these residues to be aligned, or want to give them individual weights, you will have to edit the library file yourself or use the -force_aln option (cf FAQ: I would like to force some residues to be aligned). A value of N*N * 1000 (N being the number of sequences) usually ensure the respect of a constraint.


Q: I want to use my own tree.
-----------------------------

A: Use the -usetree=<your own tree> flag:


::

  $$: t_coffee sample_seq1.fasta -usetree=sample_tree.dnd



Q: I want to align coding DNA.
------------------------------

A: Use the fasta_cdna_pair method that compares two cDNA using the best reading frame and taking frameshifts into account:


::

  $$: t_coffee three_cdna.fasta -method=cdna_fast_pair



Notice that in the resulting alignments, all the gaps are of modulo3, except one small gap in the first line of sequence hmgl_trybr. This is a framshift, made on purpose. You can realign the same sequences while ignoring their coding potential and treating them like standard DNA:


::

  $$: t_coffee three_cdna.fasta



.. warning:: This method has not yet been fully tested and is only provided 'as-is' with no warranty. Any feedback will be much appreciated.

Q: I do not want to use all the possible pairs when computing the library.
--------------------------------------------------------------------------

Q: I only want to use specific pairs to compute the library.
------------------------------------------------------------

A: Simply write in a file the list of sequence groups you want to use:


::

  $$: t_coffee sample_seq1.fasta -method=clustalw_pair,clustalw_msa -lib_list=sa\
 mple_list1.lib_list -outfile=test

  ***************sample_list1.lib_list****

  2 hmgl_trybr hmgt_mouse

  2 hmgl_trybr hmgb_chite

  2 hmgl_trybr hmgl_wheat

  3 hmgl_trybr hmgl_wheat hmgl_mouse

  ***************sample_list1.lib_list****



.. note:: Note: Pairwise methods (slow_pair...) will only be applied to list of pairs of sequences, while multiple methods (clustalw_aln) will be applied to any dataset having more than two sequences.

Q: There are duplicates or quasi-duplicates in my set.
------------------------------------------------------

A: If you can remove them, this will make the program run faster, otherwise, the t_coffee scoring scheme should be able to avoid over-weighting of over-represented sequences.

*****************************
Using Structures and Profiles
*****************************

Q: Can I align sequences to a profile with T-Coffee?
----------------------------------------------------

A: Yes, you simply need to indicate that your alignment is a profile with the R tag:


::

  $$: t_coffee sample_seq1.fasta -profile=sample_aln2.aln -outfile tacos



Q: Can I align sequences two or more profiles?
----------------------------------------------

A: Yes, you, simply tag your profiles with the letter R and the program will treat them like standard sequences:


::

  $$: t_coffee -profile=sample_aln1.fasta,sample_aln2.aln -outfile tacos



Q: Can I align two profiles according to the structures they contain?
---------------------------------------------------------------------

A: Yes. As long as the structure sequences are named according to their PDB identifier:


::

  $$: t_coffee -profile=sample_profile1.aln,sample_profile2.aln -special_mode=3\
 dcoffee -outfile=aligne_prf.aln



Q: T-Coffee becomes very slow when combining sequences and structures.
----------------------------------------------------------------------

A: This is true. By default the structures are feteched on the net, using RCSB. The problem arises when T-Coffee looks for the structure of sequences WITHOUT structures. One solution is to install PDB locally. In that case you will need to set two environment variables:



::

  setenv (or export) PDB_DIR='directory containing the pdb structures' setenv (\
 or export) NO_REMOTE_PDB_DIR=1




Interestingly, the observation that sequences without structures are those that take the most time to be checked is a reminder of the strongest rational argument that I know of against torture: any innocent would require the maximum amount of torture to establish his/her innocence, which sounds...ahem...strange., and at least inneficient. Then again I was never struck by the efficiency of the Bush administration.


Q: Can I use a local installation of PDB?
-----------------------------------------

A: Yes, T-Coffe supports three types of installations:


 -an add-hoc installation where all your structures are in a directory, under the form pdbid.pdb or pdbid.id.Z or pdbid.pdb.gz. In that case, all you need to do is set the environement variables correctly:


::

  setenv (or export) PDB_DIR='directory containing the pdb structures' setenv (\
 or export) NO_REMOTE_PDB_DIR=1



 -a standard pdb installation using the all section of pdb. In that case, you must set the variables to:


::

  setenv (or export) PDB_DIR='<some absolute path>/data/structures/all/pdb/' se\
 tenv (or export) NO_REMOTE_PDB_DIR=1



 -a standard pdb installation using the divided section of pdb:


::

  setenv (or export) PDB_DIR='<some absolute path>/data/structures/divided/pdb/\
 ' setenv (or export) NO_REMOTE_PDB_DIR=1



If you need to do more clever things, you should know that all the PDB manipulation is made in T-Coffee by a perl script named extract_from_pdb. You can extract this script from T-Coffee:


::

  t_coffee -other_pg unpack_extract_from_pdb  chmod u+x extract_from_pdb



You can then edit the script to suit your needs. T-Coffee will use your edited version if it is in the current directory. It will issue a warning that it used a local version.


If you make extensive modifications, I would appreciate you send me the corrected file so that I can incorporate it in the next distribution.

By default, T-Coffee also requires two important PDB files declared using the two following variables. These variables do not need to be set if the considered files are in the cache directory (default behavior): 

::
  export PDB_ENTRY_TYPE_FILE=<location of the file pdb_entry_type.txt>
  Found at: ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_entry_type.txt

and...

::
  export PDB_UNREALEASED_FILE=<location of the file unrealeased.xml>
  Found at: http://www.rcsb.org/pdb/rest/getUnreleased

.. warning:: Since the file unreleased.xml is not part of the pdb distribution, T-Coffee will make an attempt to obtain it even when using the NO_REMOTE_PDB_DIR=1 mode. You must therefore make sure that the file PDB_UNREALEASED_FILE is pointing to is read and write.

********************
Alignment Evaluation
********************

Q: How good is my alignment?
----------------------------

A: see what is the color index !!!


Q: What is that color index?
----------------------------

A: T-Coffee can provide you with a measure of consistency among all the methods used. You can produce such an output using:


::

  $$: t_coffee sample_seq1.fasta -output=html



This will compute your_seq.score_html that you can view using netscape. An alternative is to use score_ps or score_pdf that can be viewed using ghostview or acroread, score_ascii will give you an alignment that can be parsed as a text file.


A book chapter describing the CORE index is available on:


http://www.tcoffee.org/Publications/Pdf/core.pp.pdf


Q: Can I evaluate alignments NOT produced with T-Coffee?
--------------------------------------------------------

A: Yes. You may have an alignment produced from any source you like. To evaluate it do:


::

  $$: t_coffee -infile=sample_aln1.aln -lib=sample_aln1.tc_lib -special_mode=eva\
 luate



If you have no library available, the library will be computed on the fly using the following command. This can take some time, depending on your sample size. To monitor the progress in a situation where the default library is being built, use:


::

  $$: t_coffee -infile=sample_aln1.aln -special_mode evaluate



Q: Can I compare two alignments?
--------------------------------

A: Yes. You can treat one of your alignments as a library and compare it with the second alignment:


::

  $$: t_coffee -infile=sample_aln1_1.aln -aln=sample_aln1_2.aln -special_mode=ev\
 aluate



If you have no library available, the library will be computed on the fly using the following command. This can take some time, depending on your sample size. To monitor the progress in a situation where the default library is being built, use:


::

  $$: t_coffee -infile=sample_aln1.aln -special_mode evaluate



Q: I am aligning sequences with long regions of very good overlap.
------------------------------------------------------------------

A: Increase the ktuple size ( up to 4 or 5 for DNA) and up to 3 for proteins:


::

  $$: t_coffee sample_seq1.fasta -ktuple=3



This will speed up the program. It can be very useful, especially when aligning ESTs.


Q: Why is T-Coffee changing the names of my sequences!!!!
---------------------------------------------------------

A: If there is no duplicated name in your sequence set, T-Coffee's handling of names is consistent with Clustalw, (Cf Sequence Name Handling in the Format section). If your dataset contains sequences with identical names, these will automatically be renamed to:


::

  ************************

  >seq1

  >seq1

  ************************

  >seq1

  >seq1_1

  ************************



.. warning:: The behaviour is undefined when this creates two sequence with a similar names.

************************
Improving Your Alignment
************************

Q: How can I edit my alignment manually?
----------------------------------------

A: Use jalview, a Java online MSA editor: www.jalview.org


Q: Have I improved or not my alignment?
---------------------------------------

A: Using structural information is the only way to establish whether you have improved or not your alignment. The CORE index can also give you some information.
