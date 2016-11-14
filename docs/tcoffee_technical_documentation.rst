################################
T-Coffee Technical Documentation 
################################

.. Note:: This documentation covers the T-Coffee different flag usage; it describes how to use different T-Coffee options, with their format and restrictions. The most up to date version is available from our `webpage <http://www.tcoffee.org>`_.


*******************
General Information
*******************
T-Coffee general behavior
=========================

.. warning:: T-Coffee is not POSIX compliant (sorry L).

T-Coffee flags
--------------
This documentation gives a list of all the flags that can be used to modify the behavior of T-Coffee; for your convenience, we have grouped them according to their nature. When running T-Coffee, some options and their associated values or parameters will be displayed on screen (command 1). You can also display the list of all the flags (command 2) or a single flag (command 3) used in the version of T-Coffee you are using along with their default value type.

::

  Command 1:
  $$: t_coffee
  
  Command 2: 
  $$: t_coffee -help
  
  Command 3:
  $$: t_coffee -help -<flag>
 
Syntax of T-Coffee commands
---------------------------
You can use any kind of separator you want (i.e. ,; <space>=). The syntax used in this document is meant to be consistent with that of ClustalW. However, in order to take advantage of the automatic filename compleation provided by many shells, you can replace '=' and ',' with a space.

Entering the right parameters
-----------------------------
There are many ways to enter parameters in T-Coffee, see the **-parameters** flag. In general you will not need to use these complicated parameters, yet, if you find yourself typing long command lines on a regular basis, it may be worth reading this section. One may easily feel confused with the various manners in which the parameters can be passed to T-Coffee. The reason for these many mechanisms is that they allow several levels of intervention. For instance, you may install T-Coffee for all the users and decide that the defaults we provide are not the proper ones...In this case, you will need to make your own ``t_coffee_default`` file. Later on, a user may find that he/she needs to keep reusing a specific set of parameters, different from those in t_coffee_default, hence the possibility to write an extra parameter file with the flag **-parameters**. In summary, this means that **-parameters** supersede all the other options, while parameters provided via **-mode** are the weakest:

::

  -parameters > prompt parameters > -t_coffee_defaults > -mode
  

Environment variables
=====================
List of variables
-----------------
It is possible to modify T-Coffee's behavior by setting any of the following environment variables. With the bash shell, use **export VAR='value'**; with the cshell, use **set $VAR='value'**.

 - **http_proxy_4_TCOFFEE**: sets the http_proxy and HTTP_proxy values used by T-Coffee. These values supersede http_proxy and HTTP_proxy while http_proxy_4_TCOFFEE gets superseded by the command line values **-proxy** and **-email**; if you have no proxy, just set this value to an empty string.

 - **email_4_TCOFFEE**: sets the E-mail values provided to web services called upon by T-Coffee; can be overriden by the flag **-email**.

 - **DIR_4_TCOFFEE**: by default this variable is set to $HOME/.t_coffee; this is where T-Coffee expects to find its cache, tmp dir and possibly any temporary data stored by the program.

 - **TMP_4_TCOFFEE**: by default this variable is set to $HOME/.t_coffee/tmp; this is where T-Coffee stores temporary files.

 - **CACHE_4_TCOFFEE**: by default this variable is set to $HOME/.t_coffee/cache; this is where T-Coffee stores any data expensive to obtain: PDB files, structural alignments...

 - **PLUGINS_4_TCOFFEE**: by default all the third party packages are searched in the directory DIR_4_TCOFFEE/plugins/<OS>. This variable overrides the default but can also be overriden by the **-plugins** flag.

 - **NO_ERROR_REPORT_4_TCOFFEE**: by default this variable is no set; set it if you do not want the program to generate a verbose error output file (useful for running a server).

 - **PDB_DIR**:indicate the location of your local PDB installation.

 - **NO_WARNING_4_TCOFFEE**: suppresses all the warnings.

 - **UNIQUE_DIR_4_TCOFFEE**: sets all DIR_4_TCOFFEE, CACHE_4_TCOFFEE, TMP_4_TCOFFEE, PLUGINS_4_TCOFFEE


Setting up the environment variables
------------------------------------
T-Coffee can have its own environment file. This environment is kept in a file named ``$HOME/.t_coffee/t_coffee_env`` and can be edited. The value of any legal variable can be modified through that file. For instance, here are some examples of 1) using a configuration file when not requiring a proxy, 2) setting up any environment variable using the **-setenv** or 3) simply using an export.

::

  1) No proxy
  $: http_proxy_4_TCOFFEE=
  $: EMAIL_4_TCOFFEE=cedric.notredame@gmail.com

  2) Using 'setenv'
  $$: t_coffee ... -setenv ENV_4_TCOFFEE=<location>

  3) Using 'export'
  $: export ENV_4_TCOFFEE=<location>

.. hint:: export > -setenv > -proxy, -email > t_coffee_env > default environment

.. note:: When you use **-setenv** for PATH, the value you provide is concatenated TO THE BEGINNING of the current PATH value. This way you can force T-Coffee to use a specific version of an aligner.

Meta Parameters
===============
Global behavior
---------------
 - **no flag**
*If no flag is provided, your sequence datset must be the first argument; when you do so, the name of your file is used as a name prefix for every output file of the program (changing the extension according to the type of result).*

 - **-mode**
*A T-Coffee mode is a hard coded mode command line calling to specific options predetermined and optimized. By default, they are not used and should be called upon. Here are some examples: expresso, mcoffee, rcoffee, evaluate, accurate, procoffee... These modes have been designed to deliver the best results possible for a specific task; they can work without any parameters but can be controled and modified extensively with extra parameters.*

 - **-parameters**
*Input needs to be a file containing extra parameters for T-Coffee. Parameters read this way behave as if they had been added on the right end of the command line that they either supersede (one value parameter) or complete (list of values). Here is an example of usage that will cause T-Coffee to apply the* **fast_pair** *method onto the sequences contained in* ``sample_seq1.fasta``. *If you wish, you can also pipe these arguments into T-Coffee by naming the parameter file 'stdin' (as a rule, any file named stdin is expected to receive its content via the stdin)*

.. warning:: This parameter file can ONLY contain valid parameters; comments are not allowed. Parameters passed this way will be checked like normal parameters.

::

  $$: t_coffee -parameters=sample_file.param
  or
  $: cat sample_param_file.param | t_coffee -parameters=stdin
  
  **********sample_file.param***********
   -in=Ssample_seq1.fasta,Mfast_pair
   -output=msf_aln
  **************************************

 - **-t_coffee_defaults**
*Input needs to be a file; it will tells the program to use some default parameter file for T-Coffee. The format of that file is the same as the one used with* **-parameters**. *The file used is either:*

 1) <file name> if a name has been specified
 2) ~/.t_coffee_defaults if no file was specified
 3) The file indicated by the environment variable TCOFFEE_DEFAULTS

 - **-evaluate**
*Replaces the former flag* **-score** *which is no longer supported. This flag toggles on the evaluate mode and causes T-Coffee to evaluate a precomputed MSA provided via* **-infile=<MSA>**. *The main purpose of this flag is to let you control every aspect of the evaluation, yet it is advisable to use predefined parameterization* **-mode=evaluate**. *The flag* **-output** *must be set to an appropriate format (refer to the subsection 'Alignments Flags').*

::

  $$: t_coffee -infile=sample_aln1.aln -mode=evaluate

  $$: t_coffee -infile=sample_seq1.aln -in Lsample_lib1.tc_lib -mode=evaluate



 - **-convert [cw]**  
*By default, is turned off. It toggles on the conversion mode and causes T-Coffee to convert the sequences, alignments, libraries or structures provided via the -infile and -in flags. The output format must be set via the -output flag. This flag can also be used if you simply want to compute a library (i.e. you have an alignment and you want to turn it into a library). This flag is ClustalW compliant.*


Misc parameters
---------------
 - **-version**
*Returns the current version number*

 - **-proxy**
*Sets the proxy used by HTTP_proxy AND http_proxy. Setting with the propmpt supersedes ANY other setting. Note that if you use no proxy, you should set* **-proxy**.

 - **-email**
*Sets your email value as provided to web services*

 - **-check_configuration**
*Checks your system to determine whether all the programs T-Coffee can interact with are installed or not.*

 - **-cache**
*By default, t_coffee stores in a cache directory, the results of computationally expensive (structural alignment) or network intensive (BLAST search) operations. The usage is the following:* **-cache=<use, update, ignore, <filename>**.

 - **-update**
*Causes a wget access that checks whether the t_coffee version you are using needs updating.*

 - **-full_log**
*Requires a file name as parameter; it causes t_coffee to output a full log file that contains all the input/output files.*

 - **-plugins**
*The input parameter has to be a directory; the directory is where all third pirty packages used by T-Coffee are kept as an alternative, you can also set the environment variable PLUGINS_4_TCOFFEE. The default directory is ~/.t_coffee/plugins/.*

 - **-other_pg**
*Some rumours claim that Tetris is embedded within T-Coffee and could be ran using some special set of commands. We wish to deny these rumours, although we may admit that several interesting reformatting programs are now embedded in T-Coffee and can be ran through the* **-other_pg** *flag. Among these other programs you can find seq_reformat, aln_compare, extract_from_pdb, irmsd, etc...*

::

  $$: t_coffee -other_pg=seq_reformat

  $$: t_coffee -other_pg=unpack_all

  $$: t_coffee -other_pg=unpack_extract_from_pdb


Input
=====
The most important "-in" flag
---------------------------
The -in flag and its Identifier TAGS <-in> is the real grinder of T-Coffee. Sequences, methods and alignments all pass through so that T-Coffee can turn it all into a single list of constraints (the library). Everything is done automatically with T-Coffee going through each file to extract the sequences it contains. The methods are then applied to the sequences. Pre-compiled constraint list can also be provided. Each file provided via this flag must be preceded with a symbol (Identifier TAG) that indicates its nature to T-Coffee. The TAGs currently supported are the following:

::

  P  PDB structure
  S  for sequences (use it as well to treat an MSA as unaligned sequences)
  M  Methods used to build the library
  L  Pre-computed T-Coffee library
  A  Multiple Alignments that must be turned into a Library
  X  Substitution matrices.
  R  Profiles. 
  
This is a legal multiple alignments that will be treated as single sequences (the sequences it contains will not be realigned). If you do not want to use the TAGS, you will need to use the following flags in replacement of -in. Do not use the TAGS when using these flags:
 - aln:  Alignments  (A)
 - profile: Profiles  (R)
 - method: Method  (M)
 - seq: Sequences  (S)
 - lib: Libraries (L)


The **-in** common usageis: **-in=[<P,S,A,L,M,X><name>]**. *By default: -in=Mlalign_id_pair,Mclustalw_pair*

.. note:: Note: -in can be replaced with the combined usage of -aln, iprofile, .pdb, .lib, -method.

*See the box for an explanation of the -in flag. The following argument passed via -in*

::

  $$: t_coffee -in=Ssample_seq1.fasta,Asample_aln1.aln,Asample_aln2.msf,Mlalign_id_pair, \
      Lsample_lib1.tc_lib -outfile=outaln


This command will trigger the following chain of events:

1) *Gather all the sequences. Sequences within all the provided files are pooled together. Format recognition is automatic. Duplicates are removed (if they have the same name). Duplicates in a single file are only tolerated in FASTA format file, although they will cause sequences to be renamed.In the above case, the total set of sequences will be made of sequences contained in sequences1.seq, alignment1.aln, alignment2.msf and library.lib, plus the sequences initially gathered by -infile.*

2) *Turn alignments into libraries. alignment1.aln and alignment2.msf will be read and turned into libraries. Another library will be produced by applying the method lalign_id_pair to the set of sequences previously obtained (1). The final library used for the alignment will be the combination of all this information.*

*Note as well the following rules:*

1) Order: The order in which sequences, methods, alignments and libraries are fed in is irrelevant.
2) Heterogeneity: There is no need for each element (A, S, L) to contain the same sequences.
3) No Duplicate: Each file should contain only one copy of each sequence. Duplicates are only allowed in FASTA files but will cause the sequences to be renamed.
4) Reconciliation: If two files (for instance two alignments) contain different versions of the same sequence due to an indel, a new sequence will be reconstructed and used instead. This can be useful if you are trying to combine several runs of blast, or structural information where residues may have been deleted. However substitutions are forbidden. If two sequences with the same name cannot be merged, they will cause the program to exit with an information message.

::

  aln 1:      hgab1 AAAAABAAAAA
  aln 2:      hgab1 AAAAAAAAAACCC
  consensus:  hgab1 AAAAABAAAAACCC


5) Methods: The method describer can either be built in (See ### for a list of all the available methods) or be a file describing the method to be used. The exact syntax is provided in part 4 of this manual.
6) Substitution Matrices: If the method is a substitution matrix (X) then no other type of information should be provided. This command results in a progressive alignment carried out on the sequences in seqfile. The procedure does not use any more the T-Coffee concistency based algorithm, but switches to a standard progressive alignment algorithm (like ClustalW or Pileup) much less accurate. In this context, appropriate gap penalties should be provided. The matrices are in the file ``source/matrices.h``. *Ad Hoc* matrices can also be provided by the user (see the matrices format section at the end of this manual).

::

  $$: t_coffee sample_seq1.fasta -in=Xpam250mt -gapopen=-10 -gapext=-1

   
.. warning:: The matrix **X** does not have the same effect as using the **-matrix** flag, which defines the matrix that will be used while compiling the library while the Xmatrix defines the matrix used when assembling the final alignment.

Other sequence input flags
--------------------------
 - **-infile [cw]**
*To remain compatible with ClustalW, it is possible to indicate the sequences with this flag. Common multiple sequence alignments format constitute a valid input format. T-Coffee automatically removes the gaps before doing the alignment. This behaviour is different from that of ClustalW where the gaps are kept.*

::

  $$: t_coffee -infile=sample_seq1.fasta



 - **-get_type**
*Forces t_coffee to identify the sequences type (PROTEIN, DNA).*

 - **-type [cw]**
.. warning:: In case of low complexity or short sequences, it is recommended to set the type manually.

The common usage is **-type=DNA  PROTEIN DNA_PROTEIN**. *Default: -type=<automatically set>*. *This flag sets the type of the sequences. If omitted, the type is guessed automatically. This flag is compatible with ClustalW.*

 - **-seq**

The common usage is **-seq=[<P,S><name>]**. The flag -seq is now the recommended flag to provide your sequences; it behaves mostly like the -in flag.

 - **-seq_source**

The common usage is **-seq_source=<ANY or _LS or LS >**. *You may not want to combine all the provided sequences into a single sequence list. You can do by specifying that you do not want to treat all the -in files as potential sequence sources.The flag -seq_source=_LA indicates that neither sequences provided via the A (Alignment) flag or via the L (Library flag) should be added to the sequence list. The flag -seq_source=S means that only sequences provided via the S tag will be considered. All the other sequences will be ignored.*

.. note:: This flag is mostly designed for interactions between T-Coffee and T-CoffeeDPA (the large scale version of T-Coffee).

Other input flags (structure, tree, profile)
--------------------------------------------
 - **-pdb**
The common usage is **-pdb=<pdbid1>,<pdbid2>...[Max 200]** *It reads or fetch a pdb file. It is possible to specify a chain or even a sub-chain: PDBID(PDB_CHAIN)[opt] (FIRST,LAST)[opt]. It is also possible to input structures via the -in flag. In that case, you will need to use the TAG identifier: -in Ppdb1 Ppdb2...*

 - **-usetree**
*The common usage is* **-usetree=<tree file>**. *Default: No file specified. Format: newick tree format (ClustalW Style). This flag indicates that rather than computing a new dendrogram, t_coffee must use a pre-computed one. The tree files are in phylips format and compatible with ClustalW. In most cases, using a pre-computed tree will halve the computation time required by t_coffee. It is also possible to use trees output by ClustalW, Phylips and any other program.*

 - **-profile**
The common usage is: **-profile=[<name1>,<name2>,...] maximum of 200 profiles.** *This flag causes T-Coffee to treat multiple alignments as a single sequences, thus making it possible to make multiple profile alignments. The profile-profile alignment is controlled by -profile_mode and -profile_comparison. When provided with the -in flag, profiles must be preceded with the letter R. Note that when using -template_file, the program will also look for the templates associated with the profiles, even if the profiles have been provided as templates themselves (however it will not look for the template of the profile templates of the profile templates...)*

::

  $$: t_coffee -profile sample_aln1.aln,sample_aln2.aln -outfile=profile_aln

  $$: t_coffee -in Rsample_aln1.aln,Rsample_aln2.aln,Mslow_pair,Mlalign_id_pair \
      -outfile=profile_aln
    
  - **-profile1 [cw]** & **-profile2 [cw]**
The common usage is: **-profile1=[<prf1>], one name only** and **-profile2=[<prf2>], one name only**. *It is similar to the previous one and was provided for compatibility with ClustalW.*


Alignment Computation
=====================
Library Computation: Methods
----------------------------
 - **-lalign_n_top**

Common usage: **-lalign_n_top=<Integer>**. *Default: -lalign_n_top=10*. *Number of alignment reported by the local method (lalign).*

 - **-align_pdb_param_file** [Unsupported]

 - **-align_pdb_hasch_mode** [Unsupported]


Library Computation: Extension
------------------------------
 - **-lib_list** [Unsupported]

Common usage: **-lib_list=<filename>**. *Default:unset*. *Use this flag if you do not want the library computation to take into account all the possible pairs in your dataset. For instance*

   *Format:*

::

   2 Name1 name2
   2 Name1 name4
   3 Name1 Name2 Name3...
   * (the line 3 would be used by a multiple alignment method).*

 - **-do_normalise**
Common uage: **-do_normalise=<0 or a positive value>**. *Default:-do_normalise=1000*. *Development Only*. *When using a value different from 0, this flag sets the score of the highest scoring pair to 1000.*

 - **-extend**
Common usage: **-extend=<0,1 or a positive value>**. *Default:-extend=1*. *Development Only*. *When turned on, this flag indicates that the library extension should be carried out when performing the multiple alignment. If -extend =0, the extension is not made, if it is set to 1, the extension is made on all the pairs in the library. If the extension is set to another positive value, the extension is only carried out on pairs having a weight value superior to the specified limit.*

 - **-extend_mode**
Common usage: **-extend=<string>**. *Default:-extend=very_fast_triplet*. *Warning: Development Only*. *Controls the algorithm for matrix extension. Available modes include:relative_triplet Unsupported*, *g_coffee Unsupported*, *g_coffee_quadruplets Unsupported*, *fast_triplet Fast triplet extension*, *very_fast_triplet slow triplet extension, limited to the -max_n_pair best sequence pairs when aligning two profiles*, *slow_triplet Exhaustive use of all the triplets*, *mixt Unsupported*, *quadruplet Unsupported*, *test Unsupported*, *matrix Use of the matrix -matrix*, *fast_matrix Use of the matrix -matrix. Profiles are turned into consensus*

 - **-max_n_pair**
Common usage:** -max_n_pair=<integer>**, *Default:-extend=10*, *Development Only*, *Controls the number of pairs considered by the -extend_mode=very_fast_triplet. Setting it to 0 forces all the pairs to be considered (equivalent to -extend_mode=slow_triplet).*

 - **-seq_name_for_quadruplet** [Unsupported]

 - **-compact** [Unsupported]

 - **-clean** [Unsupported]

 - **-maximise** [Unsupported]

 - **-do_self**
*This flag causes the extension to carried out within the sequences (as opposed to between sequences). This is necessary when looking for internal repeats with Mocca.*

 - **-weight**
Common usage: **-weight=<winsimN, sim or sim_<matrix_name or matrix_file> or <integer value>**; *Default: -weight=sim*; *Weight defines the way alignments are weighted when turned into a library. Overweighting can be obtained with the OW<X> weight mode*; *winsimN indicates that the weight assigned to a given pair will be equal to the percent identity within a window of 2N+1 length centered on that pair. For instance winsim10 defines a window of 10 residues around the pair being considered. This gives its own weight to each residue in the output library. In our hands, this type of weighting scheme has not provided any significant improvement over the standard sim value.*

::

  $$: t_coffee sample_seq1.fasta -weight=winsim10 -out_lib=test.tc_lib



*sim indicates that the weight equals the average identity within the sequences containing the matched residues.*
*OW<X> will cause the sim weight to be multiplied by X*
*sim_matrix_name indicates the average identity with two residues regarded as identical when their substitution value is positive. The valid matrices names are in matrices.h (pam250mt) .Matrices not found in this header are considered to be filenames. See the format section for matrices. For instance, -weight=sim_pam250mt indicates that the grouping used for similarity will be the set of classes with positive substitutions.*

::

  $$: t_coffee sample_seq1.fasta -weight=winsim10 -out_lib=test.tc_lib


*Other groups include:
*sim_clustalw_col ( categories of clustalw marked with :)*
*sim_clustalw_dot ( categories of clustalw marked with .)*
*Value indicates that all the pairs found in the alignments must be given the same weight equal to value. This is useful when the alignment one wishes to turn into a library must be given a pre-specified score (for instance if they come from a structure super-imposition program). Value is an integer:*

::

  $$: t_coffee sample_seq1.fasta -weight=1000 -out_lib=test.tc_lib



Tree Computation
----------------
 - **-distance_matrix_mode**
Common usage: **-distance_matrix_mode=<slow, fast, very_fast>** (*Default: very_fast*). *This flag indicates the method used for computing the distance matrix (distance between every pair of sequences) required for the computation of the dendrogram.*
   *Slow  The chosen dp_mode using the extended library,*
   *fast:  The fasta dp_mode using the extended library.*
   *very_fast The fasta dp_mode using blosum62mt.*
   *ktup Ktup matching (Muscle kind)*
   *aln Read the distances on a precomputed MSA*

 - **-quicktree [cw]**
*Description: Causes T-Coffee to compute a fast approximate guide tree*. This flag is kept for compatibility with ClustalW. It indicates that:

::

  $$: t_coffee sample_seq1.fasta -distance_matrix_mode=very_fast

  $$: t_coffee sample_seq1.fasta -quicktree


Pairwise Alignment Computation
------------------------------
Controlling Alignment Computation. Most parameters in this section refer to the alignment mode fasta_pair_wise and cfatsa_pair_wise. When using these alignment modes, things proceed as follow:

1) Sequences are recoded using a degenerated alphabet provided with **-sim_matrix**
2) Recoded sequences are then hashed into ktuples of size <-ktup>
3) Dynamic programming runs on the <-ndiag> best diagonals whose score is higher than **-diag_threshold**, the way diagonals are scored is controlled via **-diag_mode**.
4) The Dynamic computation is made to optimize either the library scoring scheme (as defined by the **-in** flag) or a substitution matrix as provided via the **-matrix** flag. The penalty scheme is defined by **-gapopen** and **-gapext**. If **-gapopen** is undefined, the value defined in **-cosmetic_penalty** is used instead.
5) Terminal gaps are scored according to **-tg_mode**.


 - **-dp_mode**
Common usage: **-dp_mode=<string>** (*Default: -dp_mode=cfasta_fair_wise*). This flag indicates the type of dynamic programming used by the program. Users may find by looking into the code that other modes with fancy names exists (viterby_pair_wise...) Unless mentioned in this documentation, these modes are not supported.
::

  $$: t_coffee sample_seq1.fasta -dp_mode myers_miller_pair_wise


gotoh_pair_wise: implementation of the gotoh algorithm (quadratic in memory and time)
myers_miller_pair_wise: implementation of the Myers and Miller dynamic programming algorithm ( quadratic in time and linear in space). This algorithm is recommended for very long sequences. It is about 2 times slower than gotoh and only accepts tg_mode=1or 2 (i.e. gaps penalized for opening).
fasta_pair_wise: implementation of the fasta algorithm. The sequence is hashed, looking for ktuples words. Dynamic programming is only carried out on the ndiag best scoring diagonals. This is much faster but less accurate than the two previous. This mode is controlled by the parameters -ktuple, -diag_mode and -ndiag
cfasta_pair_wise: c stands for checked. It is the same algorithm. The dynamic programming is made on the ndiag best diagonals, and then on the 2*ndiags, and so on until the scores converge. Complexity will depend on the level of divergence of the sequences, but will usually be L*log(L), with an accuracy comparable to the two first mode ( this was checked on BaliBase). This mode is controlled by the parameters -ktuple, -diag_mode and -ndiag


 - **-ktuple**
Common usage: **-ktuple=<value>** (*Default: -ktuple=1 or 2*). *Indicates the ktuple size for cfasta_pair_wise dp_mode and fasta_pair_wise. It is set to 1 for proteins, and 2 for DNA. The alphabet used for protein can be a degenerated version, set with -sim_matrix..*

 - **-ndiag**
Common usage: **-ndiag=<value>** (*Default: -ndiag=0*). *Indicates the number of diagonals used by the fasta_pair_wise algorithm (cf -dp_mode). When -ndiag=0, n_diag=Log (length of the smallest sequence)+1.* When -ndiag and -diag_threshold are set, diagonals are selected if and only if they fulfill both conditions.*

 - **-diag_mode**
Common usage: **-diag_mode=<value>** (*Default: -diag_mode=0*). *Indicates the manner in which diagonals are scored during the fasta hashing: "0" indicates that the score of a diagonal is equal to the sum of the scores of the exact matches it contains, and "1" indicates that this score is set equal to the score of the best uninterrupted segment (useful when dealing with fragments of sequences).*

 - **-diag_threshold**
Common usage: **-diag_threshold=<value>** (*Default: -diag_threshold=0*). *Sets the value of the threshold when selecting diagonals. A value of 0: indicates that -ndiag should be used to select the diagonals (cf -ndiag section).*

 - **-sim_matrix**
Common usage: **-sim_matrix=<string>** (*Default: -sim_matrix=vasiliky*). *Indicates the manner in which the amino acid alphabet is degenerated when hashing in the fasta_pairwise dynamic programming. Standard ClustalW matrices are all valid. They are used to define groups of amino acids having positive substitution values. In T-Coffee, the default is a 13 letter grouping named Vasiliky, with residues grouped as follows:*

::

  rk, de, qh, vilm, fy (other residues kept alone).


*This alphabet is set with the flag -sim_matrix=vasiliky. In order to keep the alphabet non degenerated, -sim_matrix=idmat can be used to retain the standard alphabet.*

 - **-matrix [cw]**
Common usage: **-matrix=<blosum62mt>** (*Default: -matrix=blosum62mt*). *The usage of this flag has been modified from previous versions, due to frequent mistakes in its usage. This flag sets the matrix that will be used by alignment methods within t_coffee (slow_pair, lalign_id_pair). It does not affect external methods (like clustal_pair, clustal_aln...). Users can also provide their own matrices, using the matrix format described in the appendix.*

 - **-nomatch**
Common usage: **-nomatch=<positive value>** (*Default: -nomatch=0*). *Indicates the penalty to associate with a match. When using a library, all matches are positive or equal to 0. Matches equal to 0 are unsupported by the library but non-penalized. Setting nomatch to a non-negative value makes it possible to penalize these null matches and prevent unrelated sequences from being aligned (this can be useful when the alignments are meant to be used for structural modeling).*

 - **-gapopen**
Common usage: **-gapopen=<negative value>** (*Default: -gapopen=0*). *Indicates the penalty applied for opening a gap. The penalty must be negative. If no value is provided when using a substitution matrix, a value will be automatically computed.*
*Here are some guidelines regarding the tuning of gapopen and gapext. In T-Coffee matches get a score between 0 (match) and 1000 (match perfectly consistent with the library). The default cosmetic penalty is set to -50 (5% of a perfect match). If you want to tune -gapoen and see a strong effect, you should therefore consider values between 0 and -1000.*

 - **-gapext**
Common usage: **-gapext=<negative value>** (*Default: -gapext=0*). *Indicates the penalty applied for extending a gap (cf -gapopen)*

 - **-fgapopen** [Unsupported]

 - **-fgapext** [Unsupported]

 - **-cosmetic_penalty**
Common usage: **-cosmetic_penalty=<negative value>** (*Default: -cosmetic_penalty=-50*). *Indicates the penalty applied for opening a gap. This penalty is set to a very low value. It will only have an influence on the portions of the alignment that are unalignable. It will not make them more correct, but only more pleasing to the eye ( i.e. Avoid stretches of lonely residues). The cosmetic penalty is automatically turned off if a substitution matrix is used rather than a library.*

 - **-tg_mode**
Common usage: -**tg_mode=<0, 1, or 2>** (*Default: -tg_mode=1*).
*0: terminal gaps penalized with -gapopen + -gapext*len*
*1: terminal gaps penalized with a -gapext*len*
*2: terminal gaps unpenalized.*

Weighting Schemes
-----------------
 - **-seq_weight**
Common usage: **-seq_weight=<t_coffee or <file_name>>** (*Default: -seq_weight=t_coffee*). *These are the individual weights assigned to each sequence. The t_coffee weights try to compensate the bias in consistency caused by redundancy in the sequences.*

::

   sim(A,B)=%similarity between A and B, between 0 and 1.
   weight(A)=1/sum(sim(A,X)^3)

*Weights are normalized so that their sum equals the number of sequences. They are applied onto the primary library in the following manner:*

::

   res_score(Ax,By)=Min(weight(A), weight(B))*res_score(Ax, By)


*These are very simple weights. Their main goal is to prevent a single sequence present in many copies to dominate the alignment.*

.. note:: 1) The library output by -out_lib is the un-weighted library. 2) Weights can be output using the -outseqweight flag. 3) You can use your own weights (see the format section).


Multiple Alignment Computation
------------------------------
 - **-msa_mode** [Unsupported]
Common usage: **-msa_mode=<tree,graph,precomputed>** (*Default: -evaluate_mode=tree*).

 - **-one2all**
Common usage: **-one2all=<name>**. *Will generate a one to all library with respect to the specified sequence and will then align all the sequences in turn to that sequence, in a sequence determined by the order in which the sequences were provided.*
*-profile_comparison =profile, the MSAs provided via -profile are vectorized and the function specified by -profile_comparison is used to make profile profile alignments. In that case, the complexity is NL^2*

 - **-profile_comparison**
Common usage: **-profile_mode=<fullN,profile>** (*Default: -profile_mode=full50*). *The profile mode flag controls the multiple profile alignments in T-Coffee. There are two instances where t_coffee can make multiple profile alignments:*
*1-When N, the number of sequences is higher than -maxnseq, the program switches to its multiple profile alignment mode (t_coffee_dpa).*
*2-When MSAs are provided via the -profile flag or via -profile1 and -profile2.*
*In these situations, the -profile_mode value influences the alignment computation, these values are:*
*-profile_comparison =profile, the MSAs provided via -profile are vectorized and the function specified by -profile_comparison is used to make profile profile alignments. In that case, the complexity is NL^2*
*-profile_comparison=fullN, N is an integer value that can omitted. Full indicates that given two profiles, the alignment will be based on a library that includes every possible pair of sequences between the two profiles. If N is set, then the library will be restricted to the N most similar pairs of sequences between the two profiles, as judged from a measure made on a pairwise alignment of these two profiles.*

 - **-profile_mode**
Common usage: **-profile_mode=<cw_profile_profile, muscle_profile_profile, multi_channel>** (*Default: -profile_mode=cw_profile_profile*). *When -profile_comparison=profile, this flag selects a profile scoring function.*

Alignment Post-Processing
-------------------------
 - **-clean_aln**
Common uUsage: **-clean_aln** (*Default:-clean_aln*). *This flag causes T-Coffee to post-process the multiple alignment. Residues that have a reliability score smaller or equal to -clean_threshold (as given by an evaluation that uses -clean_evaluate_mode) are realigned to the rest of the alignment. Residues with a score higher than the threshold constitute a rigid framework that cannot be altered.* *The cleaning algorithm is greedy. It starts from the top left segment of low constituency residues and works its way left to right, top to bottom along the alignment. You can require this operation to be carried out for several cycles using the -clean_iterations flag.* *The rationale behind this operation is mostly cosmetic. In order to ensure a decent looking alignment, the gop is set to -20 and the gep to -1. There is no penalty for terminal gaps, and the matrix is blosum62mt. Gaps are always considered to have a reliability score of 0. The use of the cleaning option can result in memory overflow when aligning large sequences*.

 - **-clean_threshold**
Common usage: **-clean_threshold=<0-9>** (*Default:-clean_aln=1*). See -clean_aln for details.

 - **-clean_iteration**
Common usage: **-clean_iteration=<value between 1 and >** (*Default:-clean_iteration=1*). See -clean_aln for details.

 - **-clean_evaluation_mode**
Common usage: **-clean_iteration=<evaluation_mode >** (*Default:-clean_iteration=t_coffee_non_extended*). *Indicates the mode used for the evaluation that will indicate the segments that should be realigned. See -evaluation_mode for the list of accepted modes.*

 - **-iterate**
Common usage: **-iterate=<integer>** (*Default: -iterate=0*). *Sequences are extracted in turn and realigned to the MSA. If iterate is set to -1, each sequence is realigned, otherwise the number of iterations is set by -iterate.*

Database Searches
=================
BLAST Template Selection Parameters
-----------------------------------
These parameters are used by T-Coffee when running expresso, accurate and psicoffee


-blast_server
^^^^^^^^^^^^^
  **Usage: -blast_server= EBI, NCBI or LOCAL_BLAST**

   *Default: EBI*

   *Defines whih way BLAST will be used*

-prot_min_sim
^^^^^^^^^^^^^
  **Usage: -prot_min_sim= <percent_id>**

   *Default: 40*

   *Minimum id for inclusion of a sequence in a psi-blast profile*

-prot_max_sim
^^^^^^^^^^^^^
  **Usage: -prot_max_sim= <percent_id>**

   *Default: 90*

   *Maximum id for inclusion of a sequence in a psi-blast profile.*

-prot_min_cov
^^^^^^^^^^^^^
  **Usage: -prot_min_cov= <percent>**

   *Default: 40*

   *Minimum coverage for inclusion of a sequence in a psi-blast profile*

-protein_db
^^^^^^^^^^^
  **Usage: -protein_db= <BLAST database>**

   *Default: nr*

   *Database used for construction of psi-blast profiles*

-pdb_min_sim
^^^^^^^^^^^^
  **Usage: -pdb_min_sim= <percent_id>**

   *Default: 35*

   *Minimum id for a PDB template to be selected by expresso*

-pdb_max_sim
^^^^^^^^^^^^
  **Usage: -pdb_max_sim= <percent_id>**

   *Default: 100*

   *Maximum id for a PDB template to be selected by expresso*

-pdb_min_cov
^^^^^^^^^^^^
  **Usage: -pdb_min_cov= <percent>**

   *Default: 50*

   *Minimum coverage for a PDB template to be selected by expresso.*

-pdb_db
^^^^^^^
  **Usage: -protein_db= <BLAST database>**

   *Default: pdb*

   *Database for PDB template to be selected by expresso.*

-pdb_type
^^^^^^^^^
  **Usage: -pdb_type= d,n,m,dnm,dn**

   *Default: d*

   *d: diffraction*

   *n: NMR*

   *m: model*

CPU Control
===========
Multithreading
--------------
-multi_core
^^^^^^^^^^^
  **Usage: -multi_core= templates_jobs_relax_msa**

   *Default: 0*

   *template: fetch the templates in a parallel way*

   *jobs: compute the library*

   *relax: extend the library in a parallel way*

   *msa: compute the msa in a parallel way*

   *Specifies that the steps of T-Coffee that should be multi threaded. by default all relevant steps are parallelized.*

::

  $$: t_coffee sample_seq2.fasta -multi_core jobs



   *In order to prevent the use of the parallel mode it is possible to use:*

::

  $$: t_coffee sample_seq2.fasta -multi_core no



-n_core
^^^^^^^
  **Usage: -n_core= <number of cores>**

   *Default: 0*

   *Default indicates that all cores will be used, as indicated by the environment via:*

::

  $$: t_coffee sample_seq2.fasta -multi_core jobs



Limits
------
-mem_mode
^^^^^^^^^
  **Usage: deprecated**

-ulimit
^^^^^^^
  **Usage: -ulimit=<value>**

   *Default: -ulimit=0*

   *Specifies the upper limit of memory usage (in Megabytes). Processes exceeding this limit will automatically exit. A value 0 indicates that no limit applies.*

-maxlen
^^^^^^^
  **Usage: -maxlen=<value, 0=nolimit>**

   *Default: -maxlen=1000*

   *Indicates the maximum length of the sequences.*

Aligning more than 100 sequences with DPA
-----------------------------------------
-maxnseq
^^^^^^^^
  **Usage: -maxnseq=<value, 0=nolimit>**

   *Default: -maxnseq=50*

   *Indicates the maximum number of sequences before triggering the use of t_coffee_dpa.*

-dpa_master_aln
^^^^^^^^^^^^^^^
  **Usage: -dpa_master_aln=<File, method>**

   *Default: -dpa_master_aln=NO*

   *When using dpa, t_coffee needs a seed alignment that can be computed using any appropriate method. By default, t_coffee computes a fast approximate alignment.*

   *A pre-alignment can be provided through this flag, as well as any program using the following syntax:*

::

  your_script -in <fasta_file> -out <file_name>



-dpa_maxnseq
^^^^^^^^^^^^
  **Usage: -dpa_maxnseq=<integer value>**

   *Default: -dpa_maxnseq=30*

   *Maximum number of sequences aligned simultaneously when DPA is ran. Given the tree computed from the master alignment, a node is sent to computation if it controls more than -dpa_maxnseq OR if it controls a pair of sequences having less than -dpa_min_score2 percent ID.*

-dpa_min_score1
^^^^^^^^^^^^^^^
  **Usage: -dpa_min_score1=<integer value>**

   *Default: -dpa_min_score1=95*

   *Threshold for not realigning the sequences within the master alignment. Given this alignment and the associated tree, sequences below a node are not realigned if none of them has less than -dpa_min_score1 % identity.*

-dpa_min_score2
^^^^^^^^^^^^^^^
  **Usage: -dpa_min_score2**

   *Default: -dpa_min_score2*

   *Maximum number of sequences aligned simultaneously when DPA is ran. Given the tree computed from the master alignment, a node is sent to computation if it controls more than -dpa_maxnseq OR if it controls a pair of sequences having less than -dpa_min_score2 percent ID.*

-dap_tree [NOT IMPLEMENTED]
^^^^^^^^^^^^^^^^^^^^^^^^^^^
  **Usage: -dpa_tree=<filename>**

   *Default: -unset*

   *Guide tree used in DPA. This is a newick tree where the distance associated with each node is set to the minimum pairwise distance among all considered sequences.*

Using Structures
================
Generic
-------
-mode
^^^^^
  **Usage: -mode=3dcoffee**

   *Default: turned off*

   *Runs t_coffee with the 3dcoffee mode (cf next section).*

-check_pdb_status
^^^^^^^^^^^^^^^^^
  **Usage: -check_pdb_status**

   *Default: turned off*

   *Forces t_coffee to run extract_from_pdb to check the pdb status of each sequence. This can considerably slow down the program.*

3D Coffee: Using SAP
--------------------
   *It is possible to use t_coffee to compute multiple structural alignments. To do so, ensure that you have the sap program installed.*

::

  $$: t_coffee -pdb=struc1.pdb,struc2.pdb,struc3.pdb -method sap_pair



   *Will combine the pairwise alignments produced by SAP. There are currently four methods that can be interfaced with t_coffee:*

   *sap_pair: that uses the sap algorithm*

   *align_pdb: uses a t_coffee implementation of sap, not as accurate.*

   *tmaliagn_pair (http://zhang.bioinformatics.ku.edu/TM-align/)*

   *mustang_pair (http://www.cs.mu.oz.au/~arun/mustang)*

   *When providing a PDB file, the computation is only carried out on the first chain of this file. If your original file contains several chain, you should extract the chain you want to work on. You can use t_coffee -other_pg extract_from_pdb or any pdb handling program.*

   *If you are working with public PDB files, you can use the PDB identifier and specify the chain by adding its index to the identifier (i.e. 1pdbC). If your structure is an NMR structure, you are advised to provide the program with one structure only.*

   *If you wish to align only a portion of the structure, you should extract it yourself from the pdb file, using t_coffee -other_pg extract_from_pdb or any pdb handling program.*

   *You can provide t_coffee with a mixture of sequences and structure. In this case, you should use the special mode:*

::

  $$: t_coffee -mode 3dcoffee -seq 3d_sample3.fasta -template_file template_file\
 .template



Using/finding PDB templates for the Sequences
---------------------------------------------
-template_file
^^^^^^^^^^^^^^
  **Usage: -template_file =**

  **<filename,**

  **SCRIPT_scriptame,**

  **SELF_TAG**

  **SEQFILE_TAG_filename,**

  **no>**

   *Default: no*

   *This flag instructs t_coffee on the templates that will be used when combining several types of information. For instance, when using structural information, this file will indicate the structural template that corresponds to your sequences. The identifier T indicates that the file should be a FASTA like file, formatted as follows. There are several ways to pass the templates:*

   *Predefined Modes*

EXPRESSO: will use the EBI server to find _P_ templates


PSIBLAST: will use the EBI sever to find profiles


   *File name*

   *This file contains the sequence/template association it uses a FASTA-like format, as follows:*

::

  ><sequence name> _P_ <pdb template>

  ><sequence name> _G_ <gene template>

  ><sequence name> _R_ <MSA template>

  ><sequence name> _F_ <RNA Secondary Structure>

  ><sequence name> _T_ <Transmembrane Secondary Structure>

  ><sequence name> _E_ <Protein Secondary Structure>



   *Each template will be used in place of the sequence with the appropriate method. For instance, structural templates will be aligned with sap_pair and the information thus generated will be transferred onto the alignment.*

   *Note the following rule:*

   * -Each sequence can have one template of each type (structural, genomics...)*

   * -Each sequence can only have one template of a given type*

   * -Several sequences can share the same template*

   * -All the sequences do not need to have a template*

   *The type of template on which a method works is declared with the SEQ_TYPE parameter in the method configuration file:*

   * SEQ_TYPE S: a method that uses sequences*

   * SEQ_TYPE PS: a pairwise method that aligns sequences and structures*

   * SEQ_TYPE P: a method that aligns structures (sap for instance)*

   *There are 4 tags identifying the template type:*

   *_P_ Structural templates: a pdb identifier OR a pdb file*

   *_G_ Genomic templates: a protein sequence where boundary amino-acid have been recoded with ( o:0, i:1, j:2)*

   *_R_ Profile Templates: a file containing a multiple sequence alignment*

   *_F_ RNA secondary Structures*

   *More than one template file can be provided. There is no need to have one template for every sequence in the dataset.*

   *_P_, _G_, and _R_ are known as template TAGS*

   *2-SCRIPT_<scriptname>*

   *Indicates that filename is a script that will be used to generate a valid template file. The script will run on a file containing all your sequences using the following syntax:*

::

  scriptname -infile=<your sequences> -outfile=<template_file>



   *It is also possible to pass some parameters, use @ as a separator and # in place of the = sign. For instance, if you want to call the a script named blast.pl with the foloowing parameters;*

::

  blast.pl -db=pdb -dir=/local/test



   *Use*

::

  SCRIPT_blast.pl@db#pdb@dir#/local/test



   *Bear in mind that the input output flags will then be concatenated to this command line so that t_coffee ends up calling the program using the following system call:*

::

  blast.pl -db=pdb -dir=/local/test -infile=<some tmp file> -outfile=<another tm\
 p file>



   *3-SELF_TAG*

   *TAG can take the value of any of the known TAGS (_S_, _G_, _P_). SELF indicates that the original name of the sequence will be used to fetch the template:*

::

  $$: t_coffee 3d_sample2.fasta -template_file SELF_P_



   *The previous command will work because the sequences in 3d_sample3 are named*

   *4-SEQFILE_TAG_filename*

   *Use this flag if your templates are in filename, and are named according to the sequences. For instance, if your protein sequences have been recoded with Exon/Intron information, you should have the recoded sequences names according to the original:*

::

  SEQFILE_G_recodedprotein.fasta



-struc_to_use
^^^^^^^^^^^^^
  **Usage: -struc_to_use=<struc1, struc2...>**

   *Default: -struc_to_use=NULL*

   *Restricts the 3Dcoffee to a set of pre-defined structures.*

Domain Analysis
===============
Multiple Local Alignments
-------------------------
It is possible to compute multiple local alignments, using the moca routine. MOCA is a routine that allows extracting all the local alignments that show some similarity with another predefined fragment.


'mocca' is a perl script that calls t-coffee and provides it with the appropriate parameters.


-domain/-mocca
^^^^^^^^^^^^^^
  **Usage: -domain**

   *Default: not set*

   *This flag indicates that t_coffee will run using the domain mode. All the sequences will be concatenated, and the resulting sequence will be compared to itself using lalign_rs_s_pair mode (lalign of the sequence against itself using keeping the lalign raw score). This step is the most computer intensive, and it is advisable to save the resulting file.*

::

  $$: t_coffee -in Ssample_seq1.fasta,Mlalign_rs_s_pair -out_lib=sample_lib1.moc\
 ca_lib -domain -start=100 -len=50



   *This instruction will use the fragment 100-150 on the concatenated sequences, as a template for the extracted repeats. The extraction will only be made once. The library will be placed in the file <lib name>.*

   *If you want, you can test other coordinates for the repeat, such as*

::

  $$: t_coffee -in sample_lib1.mocca_lib -domain -start=100 -len=60



   *This run will use the fragment 100-160, and will be much faster because it does not need to re-compute the lalign library.*

-start
^^^^^^
  **Usage: -start=<int value>**

   *Default: not set*

   *This flag indicates the starting position of the portion of sequence that will be used as a template for the repeat extraction. The value assumes that all the sequences have been concatenated, and is given on the resulting sequence.*

-len
^^^^
  **Usage: -len=<int value>**

   *Default: not set*

   *This flag indicates the length of the portion of sequence that will be used as a template.*

-scale
^^^^^^
  **Usage: -scale=<int value>**

   *Default: -scale=-100*

   *This flag indicates the value of the threshold for extracting the repeats. The actual threshold is equal to:*

   * motif_len*scale*

   *Increase the scale Increase sensitivity  More alignments( i.e. -50).*

-domain_interactive [Examples]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  **Usage: -domain_interactive**

   *Default: unset*

   *Launches an interactive mocca session.*

::

  $$: t_coffee -in Lsample_lib3.tc_lib,Mlalign_rs_s_pair -domain -start=100 -len\
 =60

  TOLB_ECOLI_212_26  211 SKLAYVTFESGR--SALVIQTLANGAVRQV-ASFPRHNGAPAFSPDGSKLAFA

  TOLB_ECOLI_165_218 164 TRIAYVVQTNGGQFPYELRVSDYDGYNQFVVHRSPQPLMSPAWSPDGSKLAYV

  TOLB_ECOLI_256_306 255 SKLAFALSKTGS--LNLYVMDLASGQIRQV-TDGRSNNTEPTWFPDSQNLAFT

  TOLB_ECOLI_307_350 306 -------DQAGR--PQVYKVNINGGAPQRI-TWEGSQNQDADVSSDGKFMVMV

  TOLB_ECOLI_351_393 350 -------SNGGQ--QHIAKQDLATGGV-QV-LSSTFLDETPSLAPNGTMVIYS

   1 * * : . .:. :

   MENU: Type Letter Flag[number] and Return: ex |10

   |x -->Set the START to x

   >x -->Set the LEN to x

   Cx -->Set the sCale to x

   Sname -->Save the Alignment

   Bx -->Save Goes back x it

   return -->Compute the Alignment

   X -->eXit

  [ITERATION 1] [START=211] [LEN= 50] [SCALE=-100] YOUR CHOICE:

  For instance, to set the length of the domain to 40, type:

  [ITERATION 1] [START=211] [LEN= 50] [SCALE=-100] YOUR CHOICE:>40[return]

  [return]

  Which will generate:

  TOLB_ECOLI_212_252 211 SKLAYVTFESGRSALVIQTLANGAVRQVASFPRHNGAPAF 251

  TOLB_ECOLI_256_296 255 SKLAFALSKTGSLNLYVMDLASGQIRQVTDGRSNNTEPTW 295

  TOLB_ECOLI_300_340 299 QNLAFTSDQAGRPQVYKVNINGGAPQRITWEGSQNQDADV 339

  TOLB_ECOLI_344_383 343 KFMVMVSSNGGQQHIAKQDLATGGV-QVLSSTFLDETPSL 382

  TOLB_ECOLI_387_427 386 TMVIYSSSQGMGSVLNLVSTDGRFKARLPATDGQVKFPAW 426

   1 : : : :: . 40

   MENU: Type Letter Flag[number] and Return: ex |10

   |x -->Set the START to x

   >x -->Set the LEN to x

   Cx -->Set the sCale to x

   Sname -->Save the Alignment

   Bx -->Save Goes back x it

   return -->Compute the Alignment

   X -->eXit

  [ITERATION 3] [START=211] [LEN= 40] [SCALE=-100] YOUR CHOICE:



   *If you want to indicate the coordinates, relative to a specific sequence, type:*

::

   |<seq_name>:start



   *Type S<your name> to save the current alignment, and extract a new motif.*

   *Type X when you are done.*

Output Control
==============
Generic
-------
Conventions Regarding Filenames
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
stdout, stderr, stdin, no, /dev/null are valid filenames. They cause the corresponding file to be output in stderr or stdout, for an input file, stdin causes the program to requests the corresponding file through pipe. No causes a suppression of the output, as does /dev/null.


Identifying the Output files automatically
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In the t_coffee output, each output appears in a line:


::

  ##### FILENAME <name> TYPE <Type> FORMAT <Format>



-no_warning
^^^^^^^^^^^
  **Usage: -no_warning**

   *Default: Switched off*

   *Suppresseswarning output.*

Alignments
----------
-outfile
^^^^^^^^
  **Usage: -outfile=<out_aln file,default,no>**

Defau TOC \o '1-1' Word did not find any entries for your table of contents.lt:-outfile=default


   *Indicates the name of the alignment output by t_coffee. If the default is used, the alignment is named <your sequences>.aln*

-output
^^^^^^^
  **Usage: -output=<format1,format2,...>**

   *Default:-output=clustalw*

   *Indicates the format used for outputting the -outfile.*

   *Supported formats are:*

   **

   *clustalw_aln, clustalw : ClustalW format.*

   *gcg, msf_aln  : MSF alignment.*

   *pir_aln : pir alignment.*

   *fasta_aln : fasta alignment.*

   *phylip : Phylip format.*

   *pir_seq : pir sequences (no gap).*

   *fasta_seq : fasta sequences (no gap).*

   **

   *As well as:*

   *score_ascii : causes the output of a reliability flag*

   *score_html : causes the output to be a reliability plot in HTML*

   *score_pdf : idem in PDF (if ps2pdf is installed on your system).*

   *score_ps : idem in postscript.*

   *More than one format can be indicated:*

::

  $$: t_coffee sample_seq1.fasta -output=clustalw,gcg, score_html



   *A publication describing the CORE index is available on:*

http://www.tcoffee.org/Publications/Pdf/core.pp.pdf


-outseqweight
^^^^^^^^^^^^^
  **Usage: -outseqweight=<filename>**

   *Default: not used*

   *Indicates the name of the file in which the sequences weights should be saved..*

-case
^^^^^
  **Usage: -case=<keep,upper,lower>**

   *Default: -case=keep*

Instructs the program on the case to be used in the output file (Clustalw uses upper case). The default keeps the case and makes it possible to maintain a mixture of upper and lower case residues.


If you need to change the case of your file, you can use seq_reformat:


::

  $$: t_coffee -other_pg seq_reformat -in sample_aln1.aln -action +lower -output\
  clustalw



-cpu
^^^^
  **Usage: deprecated**

-outseqweight
^^^^^^^^^^^^^
Usage: -outseqweight=<name of the file containing the weights applied>


Default: -outseqweight=no


Will cause the program to output the weights associated with every sequence in the dataset.


-outorder [cw]
^^^^^^^^^^^^^^
  **Usage: -outorder=<input OR aligned OR filename>**

   *Default:-outorder=input*

   *Sets the order of the sequences in the output alignment: -outorder=input means the sequences are kept in the original order. -outorder=aligned means the sequences come in the order indicated by the tree. This order can be seen as a one-dimensional projection of the tree distances. -outdorder=<filename>Filename is a legal fasta file, whose order will be used in the final alignment.*

-inorder [cw]
^^^^^^^^^^^^^
  **Usage: -inorder=<input OR aligned>**

   *Default:-inorder=aligned*

   *Multiple alignments based on dynamic programming depend slightly on the order in which the incoming sequences are provided. To prevent this effect sequences are arbitrarily sorted at the beginning of the program (-inorder=aligned). However, this affects the sequence order within the library. You can switch this off by ststing -inorder=input.*

-seqnos
^^^^^^^
  **Usage: -seqnos=<on or off>**

   *Default:-seqnos=off*

Causes the output alignment to contain residue numbers at the end of each line:


::

  T-COFFEE

  seq1 aaa---aaaa--------aa 9

  seq2 a-----aa-----------a 4

  seq1 a-----------------a 11

  seq2 aaaaaaaaaaaaaaaaaaa 19



Libraries
---------
Although, it does not necessarily do so explicitly, T-Coffee always end up combining libraries. Libraries are collections of pairs of residues. Given a set of libraries, T-Coffee makes an attempt to assemble the alignment with the highest level of consistence. You can think of the alignment as a timetable. Each library pair would be a request from students or teachers, and the job of T-Coffee would be to assemble the time table that makes as many people as possible happy...


-out_lib
^^^^^^^^
Usage: -out_lib=<name of the library,default,no>


Default:-out_lib=default


   *Sets the name of the library output. Default implies <run_name>.tc_lib*

-lib_only
^^^^^^^^^
  **Usage: -lib_only**

   *Default: unset*

   *Causes the program to stop once the library has been computed. Must be used in conjunction with the flag -out_lib*

Trees
-----
-newtree
^^^^^^^^
  **Usage: -newtree=<tree file>**

   *Default: No file specified*

   *Indicates the name of the file into which the guide tree will be written. The default will be <sequence_name>.dnd, or <run_name.dnd>. The tree is written in the parenthesis format known as newick or New Hampshire and used by Phylips (see the format section).*

.. warning:: Do NOT confuse this guide tree with a phylogenetic tree.

Reliability Estimation
======================
CORE Computation
----------------
The CORE is an index that indicates the consistency between the library of piarwise alignments and the final multiple alignment. Our experiment indicate that the higher this consistency, the more reliable the alignment. A publication describing the CORE index can be found on:


http://www.tcoffee.org/Publications/Pdf/core.pp.pdf


-evaluate_mode
^^^^^^^^^^^^^^
  **Usage: -evaluate_mode=<t_coffee_fast,t_coffee_slow,t_coffee_non_extended >**

   *Default: -evaluate_mode=t_coffee_fast*

   *This flag indicates the mode used to normalize the t_coffee score when computing the reliability score.*

   *t_coffee_fast: Normalization is made using the highest score in the MSA. This evaluation mode was validated and in our hands, pairs of residues with a score of 5 or higher have 90 % chances to be correctly aligned to one another.*

   *t_coffee_slow: Normalization is made using the library. This usually results in lower score and a scoring scheme more sensitive to the number of sequences in the dataset. Note that this scoring scheme is not any more slower, thanks to the implementation of a faster heuristic algorithm.*

   *t_coffee_non_extended: the score of each residue is the ratio between the sum of its non extended scores with the column and the sum of all its possible non extended scores.*

   *These modes will be useful when generating colored version of the output, with the -output flag:*

::

  $$: t_coffee sample_seq1.fasta -evaluate_mode t_coffee_slow -output score_asci\
 i, score_html

  $$: t_coffee sample_seq1.fasta -evaluate_mode t_coffee_fast -output score_ascii, score_html

  $$: t_coffee sample_seq1.fasta -evaluate_mode t_coffee_non_extended -output score_ascii, score_html



Generic Output
==============
Misc
----
-run_name
^^^^^^^^^
  **Usage: -run_name=<your run name>**

   *Default: no default set*

This flag causes the prefix <your sequences> to be replaced by <your run name> when renaming the default output files.


-quiet
^^^^^^
  **Usage: -quiet=<stderr,stdout,file name OR nothing>.**

   *Default:-quiet=stderr*

   *Redirects the standard output to either a file. -quiet on its own redirect the output to /dev/null.*

-align [CW]
^^^^^^^^^^^
This flag indicates that the program must produce the alignment. It is here for compatibility with ClustalW.


Structural Analysis
===================
APDB, iRMSD and tRMSD Parameters
--------------------------------
.. warning:: These flags will only work within the APDB package that can be invoked via the -other_pg parameter of T-Coffee: t_coffee -other_pg apdb -aln <your aln>

-quiet [Same as T-Coffee]
^^^^^^^^^^^^^^^^^^^^^^^^^
-run_name [Same as T-Coffee]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
-aln
^^^^
  **Usage: -aln=<file_name>.**

   *Default:none*

   *Indicates the name of the file containing the sequences that need to be evaluated. The sequences whose structure is meant to be used must be named according to their PDB identifier.*

   *The format can be FASTA, CLUSTAL or any of the formats supported by T-Coffee. APDB only evaluates residues in capital and ignores those in lower case. If your sequences are in lower case, you can upper case them using seq_reformat:*

::

  $$: t_coffee -other_pg seq_reformat -in 3d_sample4.aln -action +upper -output \
 clustalw > 3d_sample4.cw_aln



   *The alignment can then be evaluated using the defaultr of APDB:*

::

  $$: t_coffee -other_pg apdb -aln 3d_sample4.aln



   *The alignment can contain as many structures as you wish.*

-n_excluded_nb
^^^^^^^^^^^^^^
  **Usage: -n_excluded_nb=<integer>.**

   *Default:1*

   *When evaluating the local score of a pair of aligned residues, the residues immediately next to that column should not contribute to the measure. By default the first to the left and first to the right are excluded.*

-maximum_distance
^^^^^^^^^^^^^^^^^
  **Usage: -maximum_distance=<float>.**

   *Default:10*

   *Size of the neighborhood considered around every residue. If .-local_mode is set to sphere, -maximum_distance is the radius of a sphere centered around each residue. If -local_mode is set to window, then -maximum_distance is the size of the half window (i.e. window_size=-maximum_distance*2+1).*

-similarity_threshold
^^^^^^^^^^^^^^^^^^^^^
  **Usage: -similarity_threshold=<integer>.**

   *Default:70*

   *Fraction of the neighborhood that must be supportive for a pair of residue to be considered correct in APDB. The neighborhood is a sphere defined by -maximum_distance, and the support is defined by -md_threshold.*

-local_mode
^^^^^^^^^^^
  **Usage: -local_mode=<sphere,window>.**

   *Default:sphere*

   *Defines the shape of a neighborhood, either as a sphere or as a window.*

-filter
^^^^^^^
  **Usage: -filter=<0.00-1.00>.**

   *Default:1.00*

   *Defines the centiles that should be kept when making the local measure. Foir instance, -filter=0.90 means that the the 10 last centiles will be removed from the evaluation. The filtration is carried out on the iRMSD values.*

-print_rapdb [Unsupported]
^^^^^^^^^^^^^^^^^^^^^^^^^^
  **Usage: -print_rapdb (FLAG)**

   *Default:off*

   *This causes the prints out of the exact neighborhood of every considered pair of residues.*

-outfile [Same as T-Coffee]
^^^^^^^^^^^^^^^^^^^^^^^^^^^
This flag is meant to control the output name of the colored APDB output. This file will either display the local APDB score or the local iRMD, depending on the value of -color_mode. The default format is defined by -ouptut and is score_html.


-color_mode
^^^^^^^^^^^
  **Usage: -color_mode=<apdb, irmsd>**

   *Default:apdb*

This flag is meant to control the colored APDB output (local score). This file will either display the local APDB score or the local iRMD.


*****************
Building a Server
*****************
We maintain a T-Coffee server (www.tcoffee.org). We will be pleased to provide anyone who wants to set up a similar service with the sources


Environment Variables
=====================
T-Coffee stores a lots of information in locations that may be unsuitable when running a server.


By default, T-Coffee will generate and rely on the follwing directory structure:


::

  /home/youraccount/ #HOME_4_TCOFFEE

  HOME_4_TCOFFEE/.t_coffee/  #DIR_4_TCOFFEE

  DIR_4_TCOFFEE/cache #CACHE_4_TCOFFEE

  DIR_4_TCOFFEE/tmp #TMP_4_TCOFFEE

  DIR_4_TCOFFEE/methods #METHOS_4_TCOFFEE

  DIR_4_TCOFFEE/mcoffee #MCOFFEE_4_TCOFFEE



By default, all these directories are automatically created, following the dependencies suggested here.


The first step is the determination of the HOME. By default the program tries to use HOME_4_TCOFFEE, then the HOME variable and TMP or TEMP if HOME is not set on your system or your account. It is your responsibility to make sure that one of these variables is set to some valid location where the T-Coffee process is allowed to read and write.


If no valid location can be found for HOME_4_TCOFFEE, the program exits. If you are running T-Coffee on a server, we recommend to hard set the following locations, where your scratch is a valid location.


::

  HOME_4_TCOFFEE='your scratch'

  TMP_4_TCOFFEE='your scratch'

  DIR_4_TCOFFEE='your scratch'

  CACHE_4_TCOFFEE='your scratch'

  NO_ERROR_REPORT_4_TCOFFEE=1



Note that it is a good idea to have a cron job that cleans up this scratch area, once in a while.


Output of the .dnd file.
========================
A common source of error when running a server: T-Coffee MUST output the .dnd file because it re-reads it to carry out the progressive alignment. By default T-Coffee outputs this file in the directory where the process is running. If the T-Coffee process does not have permission to write in that directory, the computation will abort...


To avoid this, simply specify the name of the output tree:


 -newtree=<writable file (usually in /tmp)>


Chose the name so that two processes may not over-write each other dnd file.


Permissions
===========
The t_coffee process MUST be allowed to write in some scratch area, even when it is ran by Mr nobody... Make sure the /tmp/ partition is not protected.


Other Programs
==============
T-Coffee may call various programs while it runs (lalign2list by defaults). Make sure your process knows where to find these executables.


*******
Formats
*******
Parameter files
===============
Parameter files used with -parameters, -t_coffee_defaults, -dali_defaults... Must contain a valid parameter string where line breaks are allowed. These files cannot contain any comment, the recommended format is one parameter per line:


::

   <parameter name>=<value1>,<value2>....

   <parameter name>=.....



Sequence Name Handling
======================
Sequence name handling is meant to be fully consistent with ClustalW (Version 1.75). This implies that in some cases the names of your sequences may be edited when coming out of the program. Five rules apply:


.. note:: Naming Your Sequences the Right Way

::

  1-No Space
  Names that do contain spaces, for instance:
   >seq1 human_myc
  will be turned into
   >seq1
  It is your responsibility to make sure that the names you provide are not ambi\
 guous after such an editing. This editing is consistent with Clustalw (Version 1\
 .75)
  2-No Strange Character
  Some non alphabetical characters are replaced with underscores. These are: ';:\
 ()'
  Other characters are legal and will be kept unchanged. This editing is meant t\
 o keep in line with Clustalw (Version 1.75).
  3-> is NEVER legal (except as a header token in a FASTA file)
  4-Name length must be below 100 characters, although 15 is recommended for com\
 patibility with other programs.
  5-Duplicated sequences will be renamed (i.e. sequences with the same name in t\
 he same dataset) are allowed but will be renamed according to their original ord\
 er. When sequences come from multiple sources via the -in flag, consistency of t\
 he renaming is not guaranteed. You should avoid duplicated sequences as they wil\
 l cause your input to differ from your output thus making it difficult to track \
 data.


Automatic Format Recognition
============================
Most common formats are automatically recognized by t_coffee. See -in and the next section for more details. If your format is not recognized, use readseq or clustalw to switch to another format. We recommend Fasta.


Structures
==========
PDB format is recognized by T-Coffee. T-Coffee uses extract_from_pdb (cf -other_pg flag). extract_from_pdb is a small embeded module that can be used on its own to extract information from pdb files.


RNA Structures
==============
RNA structures can either be coded as T-Coffee libraries, with each line indicating two paired residues, or as alifold output. The selex format is also partly supported (see the seq_reformat tutorial on RNA sequences handling).


Sequences
=========
Sequences can come in the following formats: fasta, pir, swiss-prot, clustal aln, msf aln and t_coffee aln. These formats are the one automatically recognized. Please replace the '*' sign sometimes used for stop codons with an X.


Alignments
==========
Alignments can come in the following formats: msf, ClustalW, Fasta, Pir and t_coffee. The t_coffee format is very similar to the ClustalW format, but slightly more flexible. Any interleaved format with sequence name on each line will be correctly parsed:


::

  <empy line>  [Facultative]n

  <line of text>  [Required]

  <line of text> [Facultative]n

  <empty line> [Required]

  <empty line> [Facultative]n

  <seq1 name><space><seq1>

  <seq2 name><space><seq2>

  <seq3 name><space><seq3>

  <empty line> [Required]

  <empty line> [Facultative]n

  <seq1 name><space><seq1>

  <seq2 name><space><seq2>

  <seq3 name><space><seq3>

  <empty line> [Required]

  <empty line> [Facultative]n



An empty line is a line that does NOT contain amino-acid. A line that contains the ClustalW annotation (.:\*) is empty.


Spaces are forbidden in the name. When the alignment is being read, non character signs are ignored in the sequence field (such as numbers, annotation...).


.. note:: Note: a different number of lines in the different blocks will cause the program to crash or hang.

Libraries
=========
T-COFFEE_LIB_FORMAT_01
----------------------
This is currently the only supported format.


::

  !<space> TC_LIB_FORMAT_01

  <nseq>

  <seq1 name> <seq1 length> <seq1>

  <seq2 name> <seq2 length> <seq2>

  <seq3 name> <seq3 length> <seq3>

  !Comment

  (!Comment)n

  #Si1 Si2

  Ri1 Ri2 V1 (V2, V3)

  #1 2

  12 13 99 (12/0 vs 13/1, weight 99)

  12 14 70

  15 16 56

  #1 3

  12 13 99

  12 14 70

  15 16 56

  !<space>SEQ_1_TO_N



Si1: index of Sequence 1


Ri1: index of residue 1 in seq1


V1: Integer Value: Weight


V2, V3: optional values


.. note:: Note 1: There is a space between the ! And SEQ_1_TO_N

.. note:: Note 2: The last line (! SEQ_1_TO_N) indicates that:

Sequences and residues are numbered from 1 to N, unless the token SEQ_1_TO_N is omitted, in which case the sequences are numbered from 0 to N-1, and residues are from 1 to N.


Residues do not need to be sorted, and neither do the sequences. The same pair can appear several times in the library. For instance, the following file would be legal:


::

  #1 2

  12 13 99

  #1 2

  15 16 99

  #1 1

  12 14 70



It is also poosible to declare ranges of resdues rather than single pairs. For instance, the following:


::

  #0 1

  +BLOCK+ 10 12 14 99

  +BLOCK+ 15 30 40 99

  #0 2

  15 16 99

  #0 1

  12 14 70



The first statement BLOCK declares a BLOCK of length 10, that starts on position 12 of sequence 1 and position 14 of sequence 2 and where each pair of residues within the block has a score of 99. The second BLOCK starts on residue 30 of 1, residue 40 of 2 and extends for 15 residues.


Blocks can overalp and be incompatible with one another, just like single constraints.





T-COFFEE_LIB_FORMAT_02
----------------------
A simpler format is being developed, however it is not yet fully supported and is only mentioned here for development purpose.


::

  ! TC_LIB_FORMAT_02

  #S1 SEQ1 [OPTIONAL]

  #S2 SEQ2 [OPTIONAL]

  ...

  !comment [OPTIONAL]

  S1 R1 Ri1 S2 R2 Ri2 V1 (V2 V3)

  => N R1 Ri1 S2 R2 Ri2 V1 (V2 V3)

  ...



S1, S2: name of sequence 1 and 2


SEQ1: sequence of S1


Ri1, Ri2: index of the residues in their respective sequence


R1, R2: Residue type


V1, V2, V3: integer Values (V2 and V3 are optional)


Value1, Value 2 and Value3 are optional.


Library List
============
These are lists of pairs of sequences that must be used to compute a library. The format is:


::

  <nseq> <S1> <S2>

  2 hamg2 globav

  3 hamgw hemog singa

  ...



Substitution matrices.
======================
If the required substitution matrix is not available, write your own in a file using the following format:


ClustalW Style [Deprecated]
---------------------------
::

  # CLUSTALW_MATRIX FORMAT

  $

  v1

  v2 v3

  v4 v5 v6

  ...

  $



v1, v2... are integers, possibly negatives.


The order of the amino acids is: ABCDEFGHIKLMNQRSTVWXYZ, which means that v1 is the substitution value for A vs A, v2 for A vs B, v3 for B vs B, v4 for A vs C and so on.


BLAST Format [Recommended]
--------------------------
::

  # BLAST_MATRIX FORMAT

  # ALPHABET=AGCT

  A G C T

  A 0 1 2 3

  G 0 2 3 4

  C 1 1 2 3

  ...



The alphabet can be freely defined


Sequences Weights
=================
Create your own weight file, using the -seq_weight flag:


::

  # SINGLE_SEQ_WEIGHT_FORMAT_01

  seq_name1 v1

  seq_name2 v2

  ...



No duplicate allowed. Sequences not included in the set of sequences provided to t_coffee will be ignored. Order is free. V1 is a float. Un-weighted sequences will see their weight set to 1.


***************
Technical Notes
***************
These notes are only meant for internal development.


Known Problems
==============
 1) Sensitivity to sequence order: it is difficult to implement a MSA algorithm totally insensitive to the order of input of the sequences. In T-Coffee, robustness is increased by sorting the sequences alphabetically before aligning them. Beware that this can result in confusing output where sequences with similar name are unexpectedly close to one another in the final alignment.

 2) Nucleotides sequences with long stretches of Ns will cause problems to lalign, especially when using Mocca. To avoid any problem, filter out these nucleotides before running mocca.

 3) Stop codons are sometimes coded with \* in protein sequences, this will cause the program to crash or hang. Please replace the all \* signs with an X.

 4) Results can differ from one architecture to another, due rounding differences. This is caused by the tree estimation procedcure. If you want to make sure an alignment is reproducible, you should keep the associated dendrogram.


Development
===========
The following examples are only meant for internal development, and are used to insure stability from release to release

profile2list
------------
prf1: profile containing one structure


prf2: profile containing one structure


::

  $$: t_coffee Rsample_profile1.aln,Rsample_profile2.aln -mode=3dcoffee -outfile\
 =aligned_prf.aln

Command Line List
-----------------
These command lines have been checked before every release (along with the other CL in this documentation:

-external methods;

::

  $$: t_coffee sample_seq1.fasta -in=Mclustalw_pair,Mclustalw_msa,Mslow_pair -ou\
 tfile=clustal_text


-fugue_client

::

  $$: t_coffee -in Ssample_seq5.fasta Pstruc4.pdb Mfugue_pair


-A list of command lines kindly provided by James Watson (used to crash the pg before version 3.40)

::

  $$: t_coffee -in Sseq.fas P2PTC Mfugue_pair
  $$: t_coffee -in S2seqs.fas Mfugue_pair -template_file SELF_P_
  $$: t_coffee -mode 3dcoffee -in Sseq.fas P2PTC
  $$: t_coffee -mode 3dcoffee -in S2seqs.fas -template_file SELF_P_


-A list of command lines that crashed the program before 3.81

::

  $$: t_coffee sample_seq6.fasta -in Mfast_pair Msap_pair Mfugue_pair -template_\
 file template_file6.template


 -A command line to read 'relaxed' pdb files...

::

  $$: t_coffee -in Msap_pair Ssample_seq7.fasta -template_file template_file7.te\
 mplate -weight 1001 -out_lib test_lib7.tc_lib -lib_only


 -Parsing of MARNA libraries

::

  $$: t_coffee -in Lmarna.tc_lib -outfile maran.test


 -Parsing of long sequence lines:

::

  $$: t_coffee -in Asample_aln5.aln -outfile test.aln


To do list
==========
Here are some improvement we are planning to do:
 - implement UPGMA tree computation
 - implement seq2dpa_tree
 - debug dpa
 - reconciliate sequences and template when reading the template
 - add the server command lines to the checking procedure


