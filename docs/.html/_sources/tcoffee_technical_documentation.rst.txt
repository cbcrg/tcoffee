################################
T-Coffee Technical Documentation 
################################

.. warning:: This chapter has been extensively updated in 11/2016. The **T-Coffee Technical Documentation** covers the different parameters, flags and options. Everything described here should be working from T-Coffee version 9.03 and higher, yet there is a trade-off between new options present in the most recent versions and old options now deprecated, unsupported or substituted by new ones (don't worry it will specified in the documentation).

.. warning:: T-Coffee is not POSIX compliant (sorry !!!).

***************************
T-Coffee Parameters & Flags
***************************
T-Coffee general information
============================
General syntax 
---------------
About the syntax of T-Coffee command lines, just for you to know that T-Coffee is quite flexible, you can use any kind of separator you want (i.e. , ; <space> =). The syntax used in this document is meant to be consistent with that of ClustalW. However, in order to take advantage of the automatic filename completion provided by many shells, you can replace '=' and ',' with a space.

T-Coffee flags
--------------
This documentation gives a list of all the flags that can be used to modify the behavior of T-Coffee; for your convenience, we have grouped them according to their nature. When running T-Coffee, some options and their associated values or parameters will be displayed on screen (command 1). You can also display the list of all the flags (command 2) or a single flag (command 3) used in the version of T-Coffee you are using along with their default value type.

::

  Command 1:
  $$: t_coffee
  
  Command 2: 
  $$: t_coffee -help
  
  Command 3:
  $#: t_coffee -help -<flag>
 
Setting up the parameters
-------------------------
There are many ways to enter parameters in T-Coffee (see the **-parameters** flag). In general you won't need to use these complicated parameters, yet, if you find yourself typing long command lines on a regular basis it may be worth reading this section. One may easily feel confused with the various manners in which the parameters can be passed to T-Coffee. The reason for these many mechanisms is that they allow several levels of intervention. For instance, you may install T-Coffee for all the users and decide that the defaults we provide are not the proper ones...In this case, you will need to make your own ``t_coffee_default`` file. Later on, a user may find that he/she needs to keep reusing a specific set of parameters, different from those in ``t_coffee_default``, hence the possibility to write an extra parameter file with the flag **-parameters**. In summary, this means that **-parameters** supersede all the other options, while parameters provided via **-mode** are the weakest.


.. hint:: Priorities: "-parameters" > "prompt parameters" > "-t_coffee_defaults" > "-mode"
  
Setting up the variables
------------------------
It is possible to modify T-Coffee's behavior by setting any of the following environment variables. With the bash shell, use **export VAR='value'**; with the cshell, use **set $VAR='value'**. Here is a list of variables in T-Coffee:

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

T-Coffee can have its own environment file, kept in a file named ``t_coffee_env`` in the folder $HOME/.t_coffee/ and can be edited (under maintenance...). The value of any legal variable can be modified through that file. For instance, here are some examples of 1) **using a configuration file** when not requiring a proxy, 2) **setting up any environment variable** using the **-setenv** or 3) simply **using an export**.

::

  1) No proxy
  ##: http_proxy_4_TCOFFEE=
  ##: EMAIL_4_TCOFFEE=cedric.notredame@gmail.com

  2) Using 'setenv'
  ##: t_coffee ... -setenv ENV_4_TCOFFEE=<location>

  3) Using 'export'
  ##: export ENV_4_TCOFFEE=<location>

.. hint:: Priorities: "export" > "-setenv" > "-proxy,-email" > "t_coffee_env" > "default environment"

.. note:: When you use **-setenv** for PATH, the value you provide is concatenated at the beginning of the current PATH value. This way you can force T-Coffee to use a specific version of an aligner.

CPU control
-----------
Multithreading
^^^^^^^^^^^^^^
- **-multi_core** (usage:**-multi_core=[templates,jobs,relax,msa]**/default:none)
Specifies that T-Coffee should be multithreaded or not; by default all relevant steps are parallelized. The different options for the flag are the following: 

  - template: fetch the templates in a parallel way
  - jobs: compute the library
  - relax: extend the library in a parallel way
  - msa: compute the msa in a parallel way
  - no: not parallelized

- **-n_core** (usage:**-n_core= [number of cores+]**/default:none)
Default indicates that all cores will be used as indicated by the environment.

Limits
^^^^^^
- **-maxlen** (usage:**-maxlen=[value,0=nolimit]**/default:**-maxlen=1000**)
Indicates the maximum length of the sequences. 

- **-maxnseq** (usage:**-maxnseq=[value,0=nolimit]**/default:**-maxnseq=??**)
Indicates the maximum number of the sequences. 

- **-ulimit** (usage:**-ulimit=[value]**/default:**-ulimit=0**)
Specifies the upper limit of memory usage (in Megabytes) and processes exceeding this limit will automatically exit. A value 0 indicates that no limit applies.

- **-mem_mode** [Deprecated]


Meta-parameters
---------------
Global parameters
^^^^^^^^^^^^^^^^^
- **no flag**
If no flag is provided, your sequence dataset must be the first argument. When you do so, the name of your file is used as a name prefix for every output file of the program (changing the extension according to the type of result).

- **-mode**
A T-Coffee mode is a hard coded command line calling to specific options predetermined and optimized. By default, they are not used and should be called upon. Here are some examples: **expresso, mcoffee, rcoffee, evaluate, accurate, procoffee**...These modes have been designed to deliver the best results possible for a specific task; they can work without any parameters but can be controlled and modified extensively with extra parameters.

- **-parameters**
The input has to be a file containing extra parameters for T-Coffee. Parameters read this way behave as if they had been added on the right end of the command line that they either supersede one parameter or complete list of parameters. Here is an example (command 1) that will cause T-Coffee to apply the **fast_pair** method onto the sequences contained in ``sample_seq1.fasta``. If you wish, you can also pipe these arguments into T-Coffee (command 2) by naming the parameter file 'stdin' (as a rule, any file named stdin is expected to receive its content via the stdin).

.. warning:: The parameter file can ONLY contain valid parameters; comments are not allowed. Parameters passed this way will be checked like normal parameters.

::

  Command 1: defining parameters for T-coffee
  $$: t_coffee -parameters=sample_param_file.param
  
  Command 2: sending parameters to T-Coffee
  $$: cat sample_param_file.param | t_coffee -parameters=stdin
  
  **********sample_file.param***********
   -in=Ssample_seq1.fasta,Mfast_pair
   -output=msf_aln
  **************************************


- **-t_coffee_defaults**
The input has to be a file; it will tells the program to use some default parameter file for T-Coffee. The format of that file is the same as the one used with **-parameters**. The file used is either:

1) <file name> if a name has been specified
2) ~/.t_coffee_defaults if no file was specified
3) The file indicated by the environment variable **TCOFFEE_DEFAULTS**

- **-evaluate**
Replaces the former flag **-score** which is no longer supported. This flag toggles on the evaluate mode and causes T-Coffee to evaluate a precomputed MSA provided via **-infile=<MSA>**. The main purpose of this flag is to let you control every aspect of the evaluation, yet it is advisable to use predefined parameterization **-mode=evaluate**. The flag **-output** must be set to an appropriate format (refer to the subsection 'Alignments Flags').

::

  $$: t_coffee -infile=sample_aln1.aln -mode=evaluate -method proba_pair

  $$: t_coffee -infile=sample_seq1.aln -in Lsample_seq1_lib1.tc_lib -mode=evaluate


- **-convert** [cw]
By default, is turned off. It toggles on the conversion mode and causes T-Coffee to convert the sequences, alignments, libraries or structures provided via the **-infile** and **-in** flags. The output format must be set via the **-output** flag. This flag can also be used if you simply want to compute a library (i.e. you have an alignment and you want to turn it into a library). This option is ClustalW compliant.

Misc parameters
^^^^^^^^^^^^^^^
- **-version**
Returns the current version number of T-Coffee you are using.

- **-proxy**
Sets the proxy used by **HTTP_proxy** and **http_proxy**. Setting with the propmpt supersedes ANY other setting. Note that if you use no proxy, you should still set **-proxy**.

- **-email**
Sets your email value as provided for web services.

- **-cache** (usage:**-cache=[use, update, ignore, <filename>]**/default:none)
By default, T-Coffee stores in a cache directory the results of computationally expensive (structural alignment for instance) or network intensive operations (BLAST search). The usage is the following:.

- **-update**
Causes a wget access that checks whether the T-Coffee version you are using needs updating.

- **-plugins**
The input parameter has to be the directory, where all third pirty packages used by T-Coffee are kept (~/.t_coffee/plugins/ by default). As an alternative, you can also set the environment variable **PLUGINS_4_TCOFFEE** to your convenience. 

- **-other_pg** (usage:**[seq_reformat,aln_compare,extract_from_pdb,irmsd,trmsd...]**)
Some rumours claim that Tetris is embedded within T-Coffee and could be ran using some special set of commands. We wish to deny these rumours, although we may admit that several interesting reformatting programs are now embedded in T-Coffee and can be ran through the **-other_pg** flag.

::

  $$: t_coffee -other_pg=seq_reformat
  $$: t_coffee -other_pg=unpack_all


- **-check_configuration** [under evaluation]
Checks your system to determine if all the programs T-Coffee can interact with are installed or not.

- **-full_log** [under evaluation]
Requires a file name as parameter; it causes T-Coffee to output a full log file that contains all the input/output files.

Verbose parameters
^^^^^^^^^^^^^^^^^^
- **-quiet** (usage:**-quiet=[stderr,stdout,file name OR nothing]**/default:**-quiet=stderr**)
This control the verbose mode of T-Coffee from the display on the screen or to redirect to a given file; **-quiet** on its own redirect the output to /dev/null.

- **-no_warning** (usage:**-no_warning=[yes,no]**/default:none)
Suppresses all warning output of the verbose mode.


Input(s)
========
The "-in" flag
--------------
The **-in** flag and its identifier TAGs **are the real grinder of T-Coffee**. Sequences, methods, alignments, whatever...all pass through so that T-Coffee can turn them all into a single list of constraints (the library). Everything is done automatically with T-Coffee going through each file to extract the sequences it contains. The methods are then applied to the sequences. Precompiled constraint list can also be provided. Each file provided via this flag must be preceded with a symbol (the identifier TAG) that indicates its nature to T-Coffee. The common usage is **-in=[<P,S,A,L,M,X><name>]**. By default it is set up to **-in=Mlalign_id_pair,Mclustalw_pair**. This is a legal multiple alignments that will be treated as single sequences (the sequences it contains will not be realigned). The TAGs currently supported are the following:

::

  P : PDB structure
  S : Sequences (aligned or unaligned sequences)
  M : Methods used to build the library
  L : Precomputed T-Coffee library
  A : Alignments that must be turned into a Library
  X : Substitution matrices
  R : Profiles
 
If you do not want to use the TAGS, you will need to use the following flags in replacement. Do not use the TAGS when using these flags.
::

 -aln     : Alignments  (A)
 -profile : Profiles    (R)
 -method  : Method      (M)
 -seq     : Sequences   (S)
 -lib     : Libraries   (L)


.. note:: The flag **-in** can be replaced with the combined usage of -aln, -profile, -pdb, -lib, -method depending on what you want.


::

  $$: t_coffee -in=Ssample_seq1.fasta,Asample_seq1_aln2.aln,Asample_seq1_aln2.msf, \
      Mlalign_id_pair,Lsample_seq1_lib1.tc_lib -outfile=outaln


This command will trigger the following chain of events:

1) **Gather all the sequences and pool them together** (format recognition is automatic). Duplicates are removed (if they have the same name). Duplicates in a single file are only tolerated in FASTA format file, although they will cause sequences to be renamed. In the above case, the total set of sequences will be made of sequences contained in ``sample_seq1.fasta``, ``sample_seq1_aln2.aln``, ``sample_seq1_aln2.msf`` and ``sample_seq1_lib1.tc_lib``, plus the sequences initially gathered by **-infile**.

2) **Turn alignment(s) into libraries** (e.g. alignment1.aln and alignment2.msf will be read and turned into libraries). Another library will be produced by applying the method lalign_id_pair to the set of sequences previously obtained (1). The final library used for the alignment will be the combination of all this information.

This procedure follows specific rules within T-Coffee; be carefull with the following rules:

- **Order**: the order in which sequences, methods, alignments and libraries are fed in is irrelevant.
- **Heterogeneity**: there is no need for each element (A, S, L) to contain the same sequences.
- **No Duplicate**: each file should contain only one copy of each sequence. Duplicates are only allowed in FASTA files but will cause the sequences to be renamed.
- **Reconciliation**: if two files (for instance two alignments) contain different versions of the same sequence due to an indel, a new sequence will be reconstructed and used instead. This can be useful if you are trying to combine several runs of blast, or structural information where residues may have been deleted. However substitutions are forbidden. If two sequences with the same name cannot be merged, they will cause the program to exit with an information message.

::

  aln 1:      hgab1 AAAAABAAAAA
  aln 2:      hgab1 AAAAAAAAAACCC
  consensus:  hgab1 AAAAABAAAAACCC

- **Substitution Matrices**: if the method is a substitution matrix (X) then no other type of information should be provided. This command results in a progressive alignment carried out on the sequences in seqfile. The procedure does not use any more the T-Coffee concistency based algorithm, but switches to a standard progressive alignment algorithm (like ClustalW or Pileup) much less accurate. In this context, appropriate gap penalties should be provided. The matrices are in the file ``matrices.h`` in the folder **"$HOME/tcoffee/Version_XX/src/"**. *Ad hoc* matrices can also be provided by the user (see the matrices format section at the end of this manual).

::

  $$: t_coffee sample_seq1.fasta -in=Xpam250mt -gapopen=-10 -gapext=-1

   
.. warning:: The matrix **X** does not have the same effect as using the **-matrix** flag, which defines the matrix that will be used while compiling the library while the Xmatrix defines the matrix used when assembling the final alignment.

- **Methods**: the method describer can either be built-in or be a file describing the method to be used (see chapter **T-Coffee Main Documentation, Internal/External Methods In T-Coffee** for more information). The exact syntax is provided later in this technical documentation.

Sequence input flags
--------------------
- **-infile** (usage:**-infile=<filename>**) [cw]
Common multiple sequence alignments format constitute a valid input format. To remain compatible with ClustalW ([cw]) it is possible to indicate the sequences with this flag. T-Coffee automatically removes the gaps before doing the alignment, and this behaviour is different from that of ClustalW where the gaps are kept.

- **-get_type**
Forces T-Coffee to identify the sequences type (protein, DNA or RNA sequences).

- **-type** (usage:**-type=[DNA,RNA,PROTEIN]**) [cw]
This flag sets the type of the sequences. The. By default, it recognizes the sequence type. If omitted, the type is guessed automatically, but in case of low complexity or short sequences, it is recommended to set the type manually This flag is compatible with ClustalW.

- **-seq** (usage:**-seq=[<P,S><name>]**)
The flag **-seq** is now the recommended flag to provide your sequences; it behaves mostly like the **-in** flag.

- **-seq_source** (usage:**-seq_source=[ANY or _LS or LS]**) [under evaluation]
You may not want to combine all the provided sequences into a single sequence list. You can do by specifying that you do not want to treat all the **-in** files as potential sequence sources.The flag **-seq_source=_LA** indicates that neither sequences provided via the A (Alignment) flag or via the L (Library flag) should be added to the sequence list. The flag **-seq_source=S** means that only sequences provided via the S tag will be considered. All the other sequences will be ignored. This flag was mostly designed for interactions between T-Coffee and T-CoffeeDPA (the large scale version of T-Coffee) which is now deprecated !!!

Other input flags (structure, tree, profile)
--------------------------------------------
- **-pdb** (usage:**-pdb=<pdbid1>,<pdbid2>...**/[max 200])
It reads or fetch a PDB file or even to specify a chain or a sub-chain: PDBID(PDB_CHAIN)[opt] (FIRST,LAST)[opt]. It is also possible to input structures via the **-in** flag but in that case, you will need to use the TAG identifier (Ppdb1 Ppdb2...).

- **-usetree** (usage:**-usetree=<tree file>**) [cw]
This flag indicates that rather than computing a new dendrogram, T-Coffee must use a precomputed one in newick tree format (ClustalW Style). The tree files are in Phylip format and compatible with ClustalW. In most cases, using a precomputed tree will halve the computation time required by T-Coffee. It is also possible to use trees output by ClustalW, Phylip and some other tree generating software. 

- **-profile** (usage:**-profile=[<name1>,<name2>,...]**/[max 200]) 
This flag causes T-Coffee to treat multiple alignments as a single sequences, thus making it possible to make multiple profile alignments. The profile-profile alignment is controlled by **-profile_mode** and **-profile_comparison**. When provided with the **-in** flag, profiles must be preceded with the letter R. Note that when using **-template_file**, the program will also look for the templates associated with the profiles even if the profiles have been provided as templates themselves (however it will not look for the template of the profile templates of the profile templates...).

::

  Turning several MSA into a single profile:
  $$: t_coffee -profile sample_seq1.aln,sample_seq1_aln2.aln -outfile=profile_aln

  Using MSA as profiles:
  $$: t_coffee -in Rsample_seq1.aln,Rsample_seq1_aln2.aln,Mslow_pair,Mlalign_id_pair \
      -outfile=profile_aln


- **-profile1**/**-profile2** (usage:**-profile1=[<prf1>]**/**-profile2=[<prf2>]**/one name only/) [cw]
It is similar to the previous command and was provided for compatibility with ClustalW.It accepts only one name in as parameter.


Output(s)
=========
Those names, stdout, stderr, stdin, no, /dev/null are valid filenames. They cause the corresponding file to be output in stderr or stdout; for an input file, stdin causes the program to requests the corresponding file through pipe. No causes a suppression of the output, as does /dev/null. In the T-Coffee output (displayed on screen), the output results appear in the following format (last lines displayed once the job is done):


::

  OUTPUT RESULTS
  ##### File Type= <type> Format= <format> Name <filename> 
  
  File Type: can be GUIDE TREE, MSA, etc...
  Format   : can be newick, html, aln, etc...
  Name     : prefix is the name of the input, extension corresponds to the type & format.


Output files, format & names
----------------------------
- **-run_name** (usage:**-run_name=<your run name>**)
This flag causes the prefix <your sequences> to be replaced by <your run name> when renaming the default output files.

- **-align** [cw]
This flag indicates that the program MUST produce an alignment. It is here for compatibility with ClustalW.

- **-outfile** (usage:**-outfile=[out_aln,file,default,no]**)
Indicates the name of the alignment output by T-Coffee. If the default is used, the alignment is named <your sequences>.aln.

- **output** (usage:**-output=[format1,format2,...]**/default:**-output=clustalw**)
Indicates the format used for the output of the resulting alignment; more than one format can be indicated. The supported input/output formats are listed below. The scoring output files rely mainly on the T-Coffee CORE index (you can find more information `here <http://www.tcoffee.org/Publications/Pdf/core.pp.pdf>`_).

::

  Sequence format:
  - clustalw_aln, clustalw
  - gcg
  - msf_aln 
  - pir_aln, pir_seq
  - fasta_aln, fasta_seq
  - phylip
  
  Output format:
  - score_ascii : causes the output of a reliability flag
  - score_html  : causes the output to be a reliability plot in HTML
  - score_pdf   : idem in PDF (if ps2pdf is installed on your system)
  - score_ps    : idem in postscript

- **-seqnos** (usage:**-seqnos=[on,off]**/default:**-seqnos=off**)
Causes the output alignment to contain the number of residue at the end of each line.

Evaluation files
----------------
The CORE is an index that indicates the consistency between the library of piarwise alignments and the final multiple alignment. Our experiment indicate that the higher this consistency, the more reliable the alignment. A publication describing the CORE index can be found `here <http://www.tcoffee.org/Publications/Pdf/core.pp.pdf>`_. All these modes will be useful when generating colored version of the output, with the **-output** flag (cf. Chapter **T-Coffee Main Documentation, Section Preparing Your Data**).

- **-evaluate_mode** (usage:**-evaluate_mode=<t_coffee_fast,t_coffee_slow,t_coffee_non_extended>**/default:**-evaluate_mode=t_coffee_fast**)
This flag indicates the mode used to normalize the T-Coffee score when computing the reliability score. Several options are possible: 

  - **t_coffee_fast**: Normalization is made using the highest score in the MSA. This evaluation mode was validated and in our hands, pairs of residues with a score of 5 or higher have 90 % chances to be correctly aligned to one another.
  - **t_coffee_slow**: Normalization is made using the library. This usually results in lower score and a scoring scheme more sensitive to the number of sequences in the dataset. Note that this scoring scheme is not any more slower, thanks to the implementation of a faster heuristic algorithm.
  - **t_coffee_non_extended**: the score of each residue is the ratio between the sum of its non extended scores with the column and the sum of all its possible non extended scores.

::

  $$: t_coffee sample_seq1.fasta -evaluate_mode t_coffee_slow -output score_ascii

  $$: t_coffee sample_seq1.fasta -evaluate_mode t_coffee_fast -output score_ascii

  $$: t_coffee sample_seq1.fasta -evaluate_mode t_coffee_non_extended -output score_ascii
 
 
Alignments
----------
- **-outseqweight** (usage:**-outseqweight=<filename>**/default:none)
Indicates the name of the file in which the sequences weights should be saved...

- **-case** (usage:**-case=[keep,upper,lower]**/default:**-case=keep**)
Instructs the program on the case to be used in the output file (Clustalw uses upper case). The default keeps the case and makes it possible to maintain a mixture of upper and lower case residues. If you need to change the case of your file, refers to the **T-Coffee Main Documentation**.

- **-outorder** (usage:**-outorder=[input,aligned,filename]**/default:**-outorder=input**) [cw]
Sets the order of the sequences in the output alignment; by default, the sequences are kept in the original order of the input dataset. Using **-outorder=aligned** means the sequences come in the order indicated by the tree; this order can be seen as a one dimensional projection of the tree distances. It is also possible to provide a file (a legal FASTA file) whose order will be used in the final MSA via **-outdorder=<filename>**.

- **-inorder** (usage:**-inorder=[input,aligned]**/default:**-inorder=aligned**) [cw]
MSAs based on dynamic programming depend slightly on the order in which the incoming sequences are provided. To prevent this effect sequences are arbitrarily sorted at the beginning of the program (**-inorder=aligned**). However, this affects the sequence order within the library. You can switch this off by stating **-inorder=input**.

- **-cpu** [Deprecated]

Librairies & trees
------------------
- **-out_lib** (usage:**-out_lib=[name of the library,default,no]**/default:**-out_lib=default**)
Sets the name of the library output. Default implies ``<run_name>.tc_lib``.

- **-lib_only**
Causes the program to stop once the library has been computed. Must be used in conjunction with the flag **-out_lib**.

- **-newtree** (usage:**-newtree=<tree file>**)
Indicates the name of the file into which the guide tree will be written. The default will be ``<sequence_name>.dnd``, or ``<run_name.dnd>``. The tree is written in the parenthesis format known as Newick or New Hampshire and used by Phylip.

.. warning:: Do NOT confuse this guide tree with a phylogenetic tree.


Alignment Computation
=====================
Library computation: methods and extension
------------------------------------------
Although it does not necessarily do so explicitly, T-Coffee always end up combining libraries. Libraries are collections of pairs of residues. Given a set of libraries, T-Coffee tries to assemble the alignment with the highest level of consistency. You can think of the alignment as a list of constraints; the job of T-Coffee is to satisfy as many constraints as possible.

- **-lalign_n_top** (usage:**-lalign_n_top=<Integer>**/default:**-lalign_n_top=10**)
Number of alignment reported by the local method (lalign).

- **-align_pdb_param_file** [Unsupported]

- **-align_pdb_hasch_mode** [Unsupported]

- **-do_normalise** (usage:**-do_normalise=[0 or a positive value]**/default:**-do_normalise=1000**)
Development Only. When using a value different from 0, this flag sets the score of the highest scoring pair to 1000.

- **-extend** (usage:**-extend=[0,1 or a positive value]**/default:**-extend=1**)
Development only. When turned on, this flag indicates that the library extension should be carried out when performing the multiple alignment. If **-extend =0**, the extension is not made, if it is set to 1, the extension is made on all the pairs in the library. If the extension is set to another positive value, the extension is only carried out on pairs having a weight value superior to the specified limit.

- **-extend_mode** (usage:**-extend=[see below...]**/default:**-extend=very_fast_triplet**)
Development only. Controls the algorithm for matrix extension. Available SUPPORTED modes include: **fast_triplet**, **very_fast_triplet** (limited to the **-max_n_pair** best sequence pairs when aligning two profiles), **slow_triplet** (exhaustive use of all the triplets), **matrix** (use of the matrix **-matrix**) and **fast_matrix** (use of the matrix **-matrix**). The following are NOT SUPPORTED: relative_triplet, g_coffee, g_coffee_quadruplets, mixt, quadruplet, test. Profiles are turned into consensus. 

- **-max_n_pair** (usage:**-max_n_pair=[integer]**/default:**-extend=10**) 
Development only. Controls the number of pairs considered by the **-extend_mode=very_fast_triplet**. Setting it to 0 forces all the pairs to be considered equivalent to **-extend_mode=slow_triplet**).

- **-weight** (usage:**-weight=[winsimN,sim,sim_<matrix_name,matrix_file>,integer]** / default:**-weight=sim**)
Weight defines the way alignments are weighted when turned into a library. Overweighting can be obtained with the OW<X> weight mode; winsimN indicates that the weight assigned to a given pair will be equal to the percent identity within a window of 2N+1 length centered on that pair. For instance winsim10 defines a window of 10 residues around the pair being considered. This gives its own weight to each residue in the output library. However, in our hands, this type of weighting scheme has not provided any significant improvement over the standard sim value (the value indicates that all the pairs found in the alignments must be given the same weight equal to value). This is useful when the alignment one wishes to turn into a library must be given a prespecified score (for instance if they come from a structure superimposition program). 

::

  $$: t_coffee sample_seq1.fasta -weight=winsim10 -out_lib=test.tc_lib

  $$: t_coffee sample_seq1.fasta -weight=1000 -out_lib=test.tc_lib

  $$: t_coffee sample_seq1.fasta -weight=sim_pam250mt -out_lib=test.tc_lib


Several options are available:

  - **sim**             : indicates that the weight equals the average identity within the sequences containing the matched residues.
  - **OW<X>**           : will cause the sim weight to be multiplied by X.
  - **sim_matrix_name** : indicates the average identity with two residues regarded as identical when their substitution value is positive. The valid matrices names are in ``matrices.h`` (pam250mt). Matrices not found in this header are considered to be filenames (Refer to the next section about matrices). For instance, -weight=sim_pam250mt indicates that the grouping used for similarity will be the set of classes with positive substitutions.
  - **sim_clustalw_col**: categories of clustalw marked with ":".
  - **sim_clustalw_dot**: categories of clustalw marked with ".".


- **-lib_list** (usage:**-lib_list=<filename>**) [Unsupported] 
Use this flag if you do not want the library computation to take into account all the possible pairs in your dataset. 

- **-do_self** [Unsupported]
This flag causes the extension to carried out within the sequences (as opposed to between sequences). This is necessary when looking for internal repeats with Mocca.

- **-seq_name_for_quadruplet** [Unsupported]

- **-compact** [Unsupported]

- **-clean** [Unsupported]

- **-maximise** [Unsupported]


Tree computation
----------------
- **-distance_matrix_mode** (usage:**-distance_matrix_mode=[slow,fast,very_fast]**/default:**very_fast**)
This flag indicates the method used for computing the distance matrix (distance between every pair of sequences) required for the computation of the dendrogram.

:: 

  - slow      : the chosen dp_mode using the extended library
  - fast      : the fasta dp_mode using the extended library
  - very_fast : the fasta dp_mode using blosum62mt
  - ktup      : Ktup matching (MUSCLE-like)
  - aln       : read the distances on a precomputed MSA


- **-quicktree** [cw]
Causes T-Coffee to compute a fast approximate guide tree; this flag is kept for compatibility with ClustalW.

::

  $$: t_coffee sample_seq1.fasta -distance_matrix_mode=very_fast

  $$: t_coffee sample_seq1.fasta -quicktree


Weighting schemes
-----------------
- **-seq_weight** (usage:**-seq_weight=[t_coffee or <file_name>]**/default:**-seq_weight=t_coffee**)
These are the individual weights assigned to each sequence. The t_coffee weights try to compensate the bias in consistency caused by redundancy in the sequences.*

::

   sim(A,B)=%similarity between A and B, between 0 and 1.
   weight(A)=1/sum(sim(A,X)^3)

Weights are normalized so that their sum equals the number of sequences. They are applied onto the primary library in the following manner:

::

   res_score(Ax,By)=Min(weight(A), weight(B))*res_score(Ax, By)


These are very simple weights. Their main goal is to prevent a single sequence present in many copies to dominate the alignment.

.. note:: The library output by **-out_lib** is the unweighted library; weights can be output using the **-outseqweight** flag; you can use your own weights (see **T-Coffee Parameter Files Format** subsection).


Pairwise alignment computation
------------------------------
Most parameters in this section refer to the alignment mode **fasta_pair_wise** and **cfatsa_pair_wise**. When using these alignment modes, things proceed as follow:

1) Sequences are recoded using a degenerated alphabet provided with **-sim_matrix**
2) Recoded sequences are then hashed into ktuples of size **-ktup**
3) Dynamic programming runs on the **-ndiag** best diagonals whose score is higher than **-diag_threshold**, the way diagonals are scored is controlled via **-diag_mode**.
4) The Dynamic computation is made to optimize either the library scoring scheme (as defined by the **-in**) or a substitution matrix as provided via **-matrix**. The penalty scheme is defined by **-gapopen** and **-gapext**. If **-gapopen** is undefined, the value defined via **-cosmetic_penalty** is used instead.
5) Terminal gaps are scored according to **-tg_mode**.


- **-dp_mode** (usage:**-dp_mode=<string>**/default:**-dp_mode=cfasta_fair_wise**). 
This flag indicates the type of dynamic programming used by the program. Users may find by looking into the code that other modes with fancy names exists (viterby_pair_wise...). Unless mentioned in this documentation, these modes are NOT SUPPORTED.

::

  $$: t_coffee sample_seq1.fasta -dp_mode myers_miller_pair_wise

The possible modes are:

  - **gotoh_pair_wise**: implementation of the gotoh algorithm (quadratic in memory and time).
  - **myers_miller_pair_wise**: implementation of the Myers and Miller dynamic programming algorithm ( quadratic in time and linear in space). This algorithm is recommended for very long sequences. It is about 2 times slower than gotoh and only accepts **-tg_mode= 1 or 2** (i.e. gaps penalized for opening).
  - **fasta_pair_wise**: implementation of the fasta algorithm. The sequence is hashed, looking for ktuples words. Dynamic programming is only carried out on the ndiag best scoring diagonals. This is much faster but less accurate than the two previous. This mode is controlled by the parameters **-ktuple, -diag_mode and -ndiag**.
  - **cfasta_pair_wise**: c stands for checked but it is the same algorithm. The dynamic programming is made on the ndiag best diagonals, and then on the 2*ndiags, and so on until the scores converge. Complexity will depend on the level of divergence of the sequences, but will usually be L*log(L), with an accuracy comparable to the two first mode (this was checked on BaliBase). This mode is controlled by the parameters **-ktuple, -diag_mode and -ndiag**.

- **-ktuple** (usage:**-ktuple=[value]**/default:**-ktuple=1 or 2**).
Indicates the ktuple size for cfasta_pair_wise and fasta_pair_wise of **-dp_mode**. It is set to 1 for proteins, and 2 for DNA. The alphabet used for protein can be a degenerated version, set with **-sim_matrix**.

- **-ndiag** (usage:**-ndiag=[value]**/default:**-ndiag=0**)
Indicates the number of diagonals used by the fasta_pair_wise algorithm (**-dp_mode**). When **-ndiag=0**, n_diag=Log (length of the smallest sequence)+1. When **-ndiag & -diag_threshold** are set, diagonals are selected if and only if they fulfill both conditions.

- **-diag_mode** (usage:**-diag_mode=[value]**/default:**-diag_mode=0**)
Indicates the manner in which diagonals are scored during the fasta hashing: "0" indicates that the score of a diagonal is equal to the sum of the scores of the exact matches it contains, and "1" indicates that this score is set equal to the score of the best uninterrupted segment (useful when dealing with fragments of sequences).

- **-diag_threshold** (usage:**-diag_threshold=[value]**/default:**-diag_threshold=0**)
Sets the value of the threshold when selecting diagonals. A value of 0: indicates that **-ndiag** (seen before) should be used to select the diagonals.

- **-sim_matrix** (usage:**-sim_matrix=[string]**/default:**-sim_matrix=vasiliky**)
Indicates the manner in which the aminoacid alphabet is degenerated when hashing in the fasta_pairwise dynamic programming. Standard ClustalW matrices are all valid. They are used to define groups of aminoacids having positive substitution values. In T-Coffee, the default is a 13 letter grouping named Vasiliky, with residues grouped as follows: [RK], [DE], [QH], [VILM], [FY], and all other residues kept alone. To keep the standard alphabet non degenerated, better use **-sim_matrix=idmat**.

- **-matrix** (usage:**-matrix=[blosum62mt,...]**/default:**-matrix=blosum62mt**) [cw]
 The usage of this flag has been modified from previous versions due to frequent mistakes in its usage. This flag sets the matrix that will be used by alignment methods within T-Coffee (slow_pair,lalign_id_pair). It does not affect external methods (like clustal_pair, clustal_aln...). Users can also provide their own matrices, using the matrix format described in the appendix.
 
- **-nomatch** (usage:**-nomatch=[positive value]**/default:**-nomatch=0**)
Indicates the penalty to associate with a match. When using a library, all matches are positive or equal to 0. Matches equal to 0 are unsupported by the library but not penalized. Setting **-nomatch** to a non negative value makes it possible to penalize these null matches and prevent unrelated sequences from being aligned (this can be useful when the alignments are meant to be used for structural modeling).

- **-gapopen** (usage:**-gapopen=[negative value]**/default:**-gapopen=0**)
Indicates the penalty applied for opening a gap. The penalty must be negative. If no value is provided when using a substitution matrix, a value will be automatically computed.

- **-gapext** (usage:**-gapext=[negative value]**/default:**-gapext=0**)
Indicates the penalty applied for extending a gap. The penalty must be negative. If no value is provided when using a substitution matrix, a value will be automatically computed.

.. hint:: Here are some guidelines regarding the tuning of **-gapopen** and **-gapext**. In T-Coffee matches get a score between 0 (match) and 1000 (match perfectly consistent with the library). The default cosmetic penalty is set to -50 (5% of a perfect match). If you want to tune **-gaopen** and see a strong effect, you should therefore consider values between 0 and -1000.

- **-cosmetic_penalty** (usage:**-cosmetic_penalty=[negative value]**/default:**-cosmetic_penalty=-50**)
Indicates the penalty applied for opening a gap. This penalty is set to a very low value, it will only have an influence on the portions of the alignment that are unalignable. It will not make them more correct but only more pleasing to the eye (avoid stretches of lonely residues). The cosmetic penalty is automatically turned off if a substitution matrix is used rather than a library.

- **-tg_mode** (usage:**-tg_mode=[0,1 or 2]**/default:**-tg_mode=1**)
The values indicate a penalty on the terminal gaps: "0" a penalty of -gapopen + -gapext*len; "1" a penalty of -gapext*len; and
"2" for no penalty associated to terminal gaps.

- **-fgapopen** [Unsupported]

- **-fgapext** [Unsupported]


Multiple alignment computation
------------------------------
- **-one2all** (usage:**-one2all=[name]**)
Will generate a one to all library with respect to the specified sequence and will then align all the sequences in turn to that sequence, in a sequence determined by the order in which the sequences were provided.If **-profile_comparison=profile**, the MSAs provided via **-profile** are vectorized and the function specified by **-profile_comparison** is used to make profile-profile alignments. In that case, the complexity is NL^2.

- **-profile_comparison** (usage:**-profile_mode=[fullN,profile]**/default:**-profile_mode=full50**)
The profile mode flag controls the multiple profile alignments in T-Coffee. There are two instances where T-Coffee can make multiple profile alignments:

  1) When N, the number of sequences is higher than **-maxnseq**, the program switches to its multiple profile alignment mode (t_coffee_dpa [Unsupported]).
  2) When MSAs are provided via **-profile, -profile1 or -profile2**. In these situations, the **-profile_mode** value influences the alignment computation, these values are: a) **-profile_comparison=profile**, the MSAs provided via **-profile** are vectorized and the function specified by **-profile_comparison** is used to make profile-profile alignments. In that case, the complexity is NL^2; b) **-profile_comparison=fullN**, N is an integer value that can omitted. Full indicates that given two profiles, the alignment will be based on a library that includes every possible pair of sequences between the two profiles. If N is set, then the library will be restricted to the N most similar pairs of sequences between the two profiles, as judged from a measure made on a pairwise alignment of these two profiles.

- **-profile_mode**(usage:**-profile_mode=[cw_profile_profile, muscle_profile_profile, multi_channel]**/default:**-profile_mode=cw_profile_profile**)
When **-profile_comparison=profile**, this flag selects a profile scoring function.

- **-msa_mode** (usage:**-msa_mode=[tree,graph,precomputed]**/default:**-evaluate_mode=tree**) [Unsupported]


Large scale aligment computation [Unsupported]
----------------------------------------------
- **-dpa** (usage:**-dpa**/default:**unset**) [Unsupported]
 This flag triggers the dpa mode. This mode involves computing the <-dpa_tree> that is then resolved into an alignment on buckets of <-dpa_nseq> sequences using <-dpa_method> as an aligner. If no <-dpa_method> is provided the T-Coffee aligner is used and all the parameters (including modes) are passed to T-Coffee. Note that all the output parameters are not supported. Note tha sequences should be input using the -seq flag.

- **-dpa_nseq** (usage:**-dpa_nseq=[ineger]**/default:**-dpa_nseq=30**) 
 Specifies the maximum bucket size when using the dpa mode. The buckets are provided to the <dpa_method> or to the remainder of the command line

- **-dpa_tree** (usage:**-dpa_tree=[file, method]**/default:**-dpa_tree=kmtree**) 
 Specifies the dpa tree

::

  dpa_tree modes
  ##: catswl       --- caterpilar tree with sequences ordered according to their distance from the average swl vector
  ##:                  the swl vector of a sequence is defined as the SMith and Waterman Length against the -swlN N sequences
  ##: catlong      --- caterpilar tree with sequences ordered according to their length (longuest on root)
  ##: catshort     --- caterpilar tree with sequences ordered according to their length (shortest on root)	
  ##: swldnd       --- Kmeans ran on the <-swlN> dimention vector where each component is the SW length of a sequence against a swlN seed
  ##: kmdnd        --- use kmeans based on triaa vectors (fast)
  ##: codnd        --- clustalo guide tree, obtainned by running clustalo externaly
  ##: cwdnd        --- clustalw guide tree, obtainned by running clustalw externaly
  ##: #<pg>        --- runs pg <seq> > stdout that outputs the tree on the stdou
  ##: parttree     --- MAFFT parttree (external)
  ##: dpparttree   --- MAFFT dpparttree (external)
  ##: fastparttree --- MAFFT fastparttree (external)
 

- **-dpa_swlN** (usage:**-dpa_swlN=[integer]**/default:**-dpa_swlN=10**) 
 When computing sw vetcors, specifies the number of N seed sequences to do an all-against-N. The sequences are selected evenly among the full sequences set sorted by size.

- **-dpa_weight** (usage:**-dpa_weight=[file, method]**/default:**-dpa_weight=longuest**). Weight used to push sequences from buckets to buckets. Sequence with the highest weight gets pushed one lebel up. 
 
::

  dpa_weight modes
  ##: longuest     --- push the longuest sequence
  ##: shortest     --- push the shortest sequence
  ##: name         --- push the sequence with the shorttes name (debug purpose)
  ##: #<pg>        --- runs pg <seq> > stdout that outputs a two column file <seqname> <float weight value> (fasta w/o sequences can be read)
  ##: <file>       --- reads weights from two column file <seqname> <float weight value> (fasta w/o sequences can be read)
  ##: kmeans       --- assign to each sequence the size of its kmeans cluster, as computed using triaa vectors and requesting 100 groups
  ##: swa,iswa     --- assigns to each sequence the average length of its SW alignment (matched residues) against the -dpa_swlN seed sequences, iswa for invert
  ##: swl,iswl     --- assign each sequence a weight equal to its distance from the average sequence using a swlN component vector, iswl for invert
  ##: diaa,idiaa   --- assign each sequence a weight equal to its distance from the average sequence using a diaa component vector, idiaa fro invert
  ##: triaa,itriaa --- assign each sequence a weight equal to its distance from the average sequence using a triaa component vector, itriaa fro invert
  

Local alignments computation [Unsupported]
------------------------------------------
It is possible to compute multiple local alignments, using the moca routine. MOCCA is a routine that allows extracting all the local alignments that show some similarity with another predefined fragment. MOCCA is a perl script that calls T-Coffee and provides the appropriate parameters.

- **-domain/-mocca** (usage:**-domain**) [Unsupported]
This flag indicates that T-Coffee will run using the domain mode. All the sequences will be concatenated, and the resulting sequence will be compared to itself using **lalign_rs_s_pair** mode (lalign of the sequence against itself using keeping the lalign raw score). This step is the most computer intensive, and it is advisable to save the resulting file.

::

  $#: t_coffee -in Ssample_seq1.fasta,Mlalign_rs_s_pair -out_lib=sample_lib1.mocca_lib \
      -domain -start=100 -len=50

This instruction will use the fragment 100-150 on the concatenated sequences, as a template for the extracted repeats. The extraction will only be made once. The library will be placed in the file <lib name>. If you want, you can test other coordinates for the repeat, such as:

::

  $#: t_coffee -in sample_lib1.mocca_lib -domain -start=100 -len=60

This run will use the fragment 100-160, and will be much faster because it does not need to recompute the lalign library.

- **-start** (usage:**-start=[integer]**) [Unsupported]
This flag indicates the starting position of the portion of sequence that will be used as a template for the repeat extraction. The value assumes that all the sequences have been concatenated, and is given on the resulting sequence.

- **-len** (usage:**-len=[integer]**) [Unsupported]
This flag indicates the length of the portion of sequence that will be used as a template.

- **-scale** (usage:**-scale=[integer]**/default:**-scale=-100**) [Unsupported]
This flag indicates the value of the threshold for extracting the repeats. The actual threshold is equal to -motif_len*scale
Increase the scale Increase sensitivity  More alignments( i.e. -50).

- **-domain_interactive** [Unsupported]
Launches an interactive MOCCA session.


Alignment post-processing
-------------------------
- **-clean_aln**
This flag causes T-Coffee to post-process the MSA. Residues that have a reliability score smaller or equal to **-clean_threshold** (as given by an evaluation that uses **-clean_evaluate_mode**) are realigned to the rest of the alignment. Residues with a score higher than the threshold constitute a rigid framework that cannot be altered. The cleaning algorithm is greedy, it starts from the top left segment of low constitency residues and works its way left to right, top to bottom along the alignment. You can require this operation to be carried out for several cycles using the **-clean_iterations** flag. The rationale behind this operation is mostly cosmetic. In order to ensure a decent looking alignment, the GOP is set to -20 and the GEP to -1. There is no penalty for terminal gaps, and the matrix is blosum62mt. Gaps are always considered to have a reliability score of 0. The use of the cleaning option can result in memory overflow when aligning large sequences.

- **-clean_threshold** (usage:**-clean_threshold=[0-9]**/default:**-clean_aln=1**)
See **-clean_aln** for details.

- **-clean_iteration** (usage: **-clean_iteration=[1-X]**/default:**-clean_iteration=1**)
See **-clean_aln** for details.

- **-clean_evaluation_mode** (usage:**-clean_iteration=[evaluation_mode]**/default:**-clean_iteration=t_coffee_non_extended**)
Indicates the mode used for the evaluation that will indicate the segments that should be realigned. See **-evaluation_mode** for the list of accepted modes.

- **-iterate** (usage: **-iterate=[integer]**/default:**-iterate=0**)
Sequences are extracted in turn and realigned to the MSA. If **-iterate** is set to -1, each sequence is realigned; otherwise the number of iterations is set by **-iterate**.


Template based modes
====================
Database searches parameters
----------------------------
These parameters are used when running template based modes of T-Coffee such as EXPRESSO (3D), PSI-Cofee (homology extension), TM/PSI-Coffee (homology extension/reduced databases), or accurate (mixture of profile and 3D templates).

- **-blast_server** (sage:**-blast_server=[EBI,NCBI,LOCAL]**/default:**-blast_server=EBI**)
Defines which way BLAST will be used, either through web services or locally. To have more information about BLAST, refer to the **T-Coffee Installation** chapter.

- **-protein_db** (usage:**-protein_db=<database>**/default:nr database)
Database used for the construction of a PSI-BLAST profile.

- **-prot_min_sim** (usage:**-prot_min_sim= [%identity]**/default:40)
Minimum identity for inclusion of a sequence in a PSI-BLAST profile.

- **-prot_max_sim** (usage:**-prot_max_sim= [%identity]**/default:90)
Maximum identity for inclusion of a sequence in a PSI-BLAST profile.

- **-prot_min_cov** (usage:**-prot_min_cov=[0-100]**/default:40)
Minimum coverage for inclusion of a sequence in a PSI-BLAST profile.

- **-pdb_db** (usage:**-protein_db=[database]**/Default:PDB database)
Database for PDB template to be selected by EXPRESSO. By default it runs on the PDB, but you can specify a version installed locally on your system,

- **-pdb_type** (usage:**-pdb_type=[d,n,m,dnm,dn]**/default:**-pdb_type=d**)
The different types are as follow: "d" stands for XRAY structures (diffraction), "n" stands for NMR structures and "m" stands for models.

- **-pdb_min_sim** (usage:**-pdb_min_sim= [%identity]**/default:35)
Minimum % identity for a PDB template to be selected by Expresso.

- **-pdb_max_sim** (usage:**-pdb_max_sim= <%identity>**/default:100)
Maximum % identity for a PDB template to be selected by Expresso.

- **-pdb_min_cov** (usage:**-pdb_min_cov= <0-100t>**/default:50)
Minimum coverage for a PDB template to be selected by Expresso.


Using structures
----------------
- **-mode** (usage:**-mode=[3dcoffee,expresso]**)
Refer to the **T-Coffee Main Documentation** for the structural modes of T-Coffee.

- **-check_pdb_status**
Forces T-Coffee to run **extract_from_pdb** to check the PDB status of each sequence. This can considerably slow down the program.


Using/finding PDB templates for the sequences
---------------------------------------------
- **-struc_to_use** (usage:**-struc_to_use=[struc1, struc2...]**/ Default:none)
Restricts the 3D-Coffee to a set of predefined structures.

- **-template_file** (usage:**-template_file =[<filename>,<SCRIPT_scriptname>,SELF_TAG_,SEQFILE_TAG_<filename>,PDB,no>**)
This flag instructs T-Coffee on the templates that will be used when combining several types of information. For instance, when using structural information, this file will indicate the structural template that corresponds to your sequences. Each template will be used in place of the sequence with the appropriate method. There are several ways to pass the templates, the format of the template file being as followed "<sequence name> <TAG> <template name>":

**Using a "filename"**:

You provide to T-Coffee an existing template file either created by Expresso or on your own. Different types of templates exist but all with the similar format and rules; the type of template on which a method works is declared with the SEQ_TYPE parameter in the method configuration file. 

::

  Template file format:
  ><sequence name> _P_ <PDB template>
  ><sequence name> _G_ <gene template>
  ><sequence name> _R_ <MSA template>
  ><sequence name> _F_ <RNA Secondary Structure>
  ><sequence name> _T_ <Transmembrane Secondary Structure>
  ><sequence name> _E_ <Protein Secondary Structure>

  Template file rules:
  - Each sequence can have one template of each type (structural, genomics...)
  - Each sequence can only have one template of a given type
  - Several sequences can share the same template
  - All the sequences do not need to have a template

  Template based method types:
 - SEQ_TYPE S: a method that uses sequences
 - SEQ_TYPE PS: a pairwise method that aligns sequences and structures
 - SEQ_TYPE P: a method that aligns structures (SAP for instance)


**Using "SCRIPT_<scriptname>"**:
 
Indicates that filename is a script that will be used to generate a valid template file. The script will run on a file containing all your sequences using the following syntax (See box 1.). It is also possible to pass some parameters, use @ as a separator and # in place of the = sign (See box 2.). For instance, if you want to call the a script named **blast.pl** with the following parameters. Bear in mind that the input/output flags will then be concatenated to this command line so that T-Coffee ends up calling the program using the following system call.

::

  1) Running using a script:
  ##: scriptname -infile=<your sequences> -outfile=<template_file>

  2) Running blast.pl:
  ##: blast.pl -db=pdb -dir=/local/test
  ##: SCRIPT_blast.pl@db#pdb@dir#/local/test
  ##: blast.pl -db=pdb -dir=/local/test -infile=<some tmp file> -outfile=<another tmp file>


**Using "SELF_TAG"**:
  
TAG can take the value of any of the known TAGS (_S_, _G_, _P_). **SELF** indicates that the original name of the sequence will be used to fetch the template.

::

  $$: t_coffee three_pdb.fasta -template_file SELF_P_



**Using "SEQFILE_TAG_filename"**
  
Use this flag if your templates are in filename, and are named according to the sequences. For instance, if your protein sequences have been recoded with Exon/Intron information, you should have the recoded sequences names according to the original ``SEQFILE_G_recodedprotein.fasta``.


Using structure for MSA evaluation
----------------------------------
MSA can be evaluated using structures with T-Coffee special mode APDB/iRMSD (see **T-Coffee Main Documentation**). The requirements to run these modes are similar to running a structure based MSA with Expresso/3D-Coffee. Here we describe the different parameters associated with the APDB mode.

- **-n_excluded_nb** (usage:**-n_excluded_nb=[integer]**/default:1)
When evaluating the local score of a pair of aligned residues, the residues immediately next to that column should not contribute to the measure. By default the first to the left and first to the right are excluded.

- **-maximum_distance** (usage:**-maximum_distance=[float]**/default:10)
Size of the neighborhood considered around every residue. If .-local_mode is set to sphere, -maximum_distance is the radius of a sphere centered around each residue. If -local_mode is set to window, then -maximum_distance is the size of the half window (i.e. window_size=-maximum_distance*2+1).

- **-similarity_threshold** (usage:**-similarity_threshold=[integer]**/default:70)
Fraction of the neighborhood that must be supportive for a pair of residue to be considered correct in APDB. The neighborhood is a sphere defined by -maximum_distance, and the support is defined by **-md_threshold**.

- **-local_mode** (usage:**-local_mode=[sphere,window]**/default:sphere)
Defines the shape of a neighborhood, either as a sphere or as a window.

- **-filter** (usage:**-filter=<0.00-1.00>**/default:1.00)
Defines the centiles that should be kept when making the local measure. Foir instance, -filter=0.90 means that the the 10 last centiles will be removed from the evaluation. The filtration is carried out on the iRMSD values.

- **print_rapdb** (usage:**-print_rapdb=[FLAG]**/default:off) [Unsupported] 
This causes the prints out of the exact neighborhood of every considered pair of residues.

- **color_mode** (usage:**-color_mode=[apdb,irmsd]**/default:apdb)
This flag is meant to control the colored APDB output (local score). This file will either display the local APDB score or the local iRSMD.


*******************************
T-Coffee Parameter Files Format 
*******************************
Sequence name handling
======================
Sequence name handling is meant to be fully consistent with ClustalW (Version 1.75). This implies that in some cases the names of your sequences may be edited when coming out of the program. Five rules apply:

  1) **No space**: It is your responsibility to make sure that the names you provide are not ambiguous after such an editing. This editing is consistent with Clustalw (Version 1..75). For instance, sequences with spaces ">seq1 human_myc"  will be turned into ">seq1". 
  2) **No strange character**: Some non alphabetical characters are replaced with underscores. These are: ';:()'. Other characters are legal and will be kept unchanged. This editing is meant to keep in line with Clustalw (Version 1.75).
  3) **> is NEVER legal** (except as a header token in a FASTA file)
  4) **Name length must be below 100 characters**, although 15 is recommended for compatibility with other programs.
  5) **Duplicated sequences will be renamed** (i.e. sequences with the same name in the same dataset) are allowed but will be renamed according to their original order. When sequences come from multiple sources via the **-in** flag, consistency of the renaming is not guaranteed. You should avoid duplicated sequences as they will cause your input to differ from your output thus making it difficult to track data.


Automatic format recognition
============================
Most common formats are automatically recognized by t_coffee. See -in and the next section for more details. If your format is not recognized, use readseq or clustalw to switch to another format. We recommend FASTA.

- **Sequences**: Sequences can come in the following formats: fasta, pir, swiss-prot, clustal aln, msf aln and t_coffee aln. These formats are the one automatically recognized. Please replace the '*' sign sometimes used for stop codons with an X.

- **Structures**: PDB format is recognized by T-Coffee. T-Coffee uses extract_from_pdb (cf -other_pg flag). extract_from_pdb is a small embeded module that can be used on its own to extract information from pdb files.

- **RNA Structures**: RNA structures can either be coded as T-Coffee libraries, with each line indicating two paired residues, or as alifold output. The selex format is also partly supported (see the seq_reformat tutorial on RNA sequences handling).

- **Alignments**: Alignments can come in the following formats: msf, ClustalW, FASTA, PIR and T-Coffee. The T-Coffee format is very similar to the ClustalW format, but slightly more flexible. Any interleaved format with sequence name on each line will be correctly parsed:

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


An empty line is a line that does NOT contain amino-acid. A line that contains the ClustalW annotation (.:\*) is empty. Spaces are forbidden in the name. When the alignment is being read, non character signs are ignored in the sequence field (such as numbers, annotation...).

.. note:: A different number of lines in the different blocks will cause the program to crash or hang.


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


Note that here is a space between the ! And SEQ_1_TO_N...and the last line (! SEQ_1_TO_N) indicates that sequences and residues are numbered from 1 to N, unless the token SEQ_1_TO_N is omitted, in which case the sequences are numbered from 0 to N-1, and residues are from 1 to N. Residues do not need to be sorted, and neither do the sequences. The same pair can appear several times in the library. For instance, the following file would be legal:

::

  #1 2
  12 13 99
  #1 2
  15 16 99
  #1 1
  12 14 70


It is also possible to declare ranges of residues rather than single pairs. For instance, the following:

::

  #0 1
  +BLOCK+ 10 12 14 99
  +BLOCK+ 15 30 40 99
  #0 2
  15 16 99
  #0 1
  12 14 70


The first statement BLOCK declares a BLOCK of length 10, that starts on position 12 of sequence 1 and position 14 of sequence 2 and where each pair of residues within the block has a score of 99. The second BLOCK starts on residue 30 of 1, residue 40 of 2 and extends for 15 residues. Blocks can overalp and be incompatible with one another, just like single constraints.


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


S1,S2: name of sequence 1 and 2.

SEQ1: sequence of S1.

Ri1, Ri2: index of the residues in their respective sequence.

R1, R2: Residue type.

V1, V2, V3: integer Values (V2 and V3 are optional).

Value1, Value 2 and Value3 are optional.

Library List
------------
These are lists of pairs of sequences that must be used to compute a library. The format is:

::

  <nseq> <S1> <S2>
  2 hamg2 globav
  3 hamgw hemog singa
  ...


Substitution matrices
=====================
If the required substitution matrix is not available, write your own in a file using one of the following format.

BLAST format [Recommended]
--------------------------
The alphabet can be freely defined

::

  # BLAST_MATRIX FORMAT
  # ALPHABET=AGCT
  A G C T
  A 0 1 2 3
  G 0 2 3 4
  C 1 1 2 3
  ...


ClustalW style [Deprecated]
---------------------------
::

  # CLUSTALW_MATRIX FORMAT
  $
  v1
  v2 v3
  v4 v5 v6
  ...
  $


v1, v2... are integers, possibly negatives. The order of the amino acids is: ABCDEFGHIKLMNQRSTVWXYZ, which means that v1 is the substitution value for A vs A, v2 for A vs B, v3 for B vs B, v4 for A vs C and so on.


Sequences weights
=================
Create your own weight file, using the -seq_weight flag:

::

  # SINGLE_SEQ_WEIGHT_FORMAT_01
  seq_name1 v1
  seq_name2 v2
  ...


No duplicate allowed. Sequences not included in the set of sequences provided to T-Coffee will be ignored. Order is free. V1 is a float. Unweighted sequences will see their weight set to 1.


Parameter files
===============
Parameter files used with **-parameters**, **-t_coffee_defaults**, **-dali_defaults**...must contain a valid parameter string where line breaks are allowed. These files cannot contain any comment, the recommended format is one parameter per line:


::

   <parameter name>=<value1>,<value2>....

   <parameter name>=.....



***************
Technical Notes
***************
These notes are only meant for internal development.

Building a T-Coffee server 
==========================
We maintain a T-Coffee server (www.tcoffee.org). We will be pleased to provide anyone who wants to set up a similar service with the sources


Environment Variables
---------------------
T-Coffee stores a lots of information in locations that may be unsuitable when running a server.


By default, T-Coffee will generate and rely on the following directory structure:


::

  Directory structure:
  ##: /home/your account/ #HOME_4_TCOFFEE
  ##: HOME_4_TCOFFEE/.t_coffee/  #DIR_4_TCOFFEE
  ##: DIR_4_TCOFFEE/cache #CACHE_4_TCOFFEE
  ##: DIR_4_TCOFFEE/tmp #TMP_4_TCOFFEE
  ##: DIR_4_TCOFFEE/methods #METHOS_4_TCOFFEE
  ##: DIR_4_TCOFFEE/mcoffee #MCOFFEE_4_TCOFFEE


By default, all these directories are automatically created, following the dependencies suggested here. The first step is the determination of the HOME. By default the program tries to use HOME_4_TCOFFEE, then the HOME variable and TMP or TEMP if HOME is not set on your system or your account. It is your responsibility to make sure that one of these variables is set to some valid location where the T-Coffee process is allowed to read and write. If no valid location can be found for HOME_4_TCOFFEE, the program exits. If you are running T-Coffee on a server, we recommend to hard set the following locations, where your scratch is a valid location.


::

  T-Coffee server configuration:
  
  ##: HOME_4_TCOFFEE='your scratch'
  ##: TMP_4_TCOFFEE='your scratch'
  ##: DIR_4_TCOFFEE='your scratch'
  ##: CACHE_4_TCOFFEE='your scratch'
  ##: NO_ERROR_REPORT_4_TCOFFEE=1


Note that it is a good idea to have a cron job that cleans up this scratch area once in a while.


Output of the .dnd file.
------------------------
A common source of error when running a server: T-Coffee MUST output the .dnd file because it re-reads it to carry out the progressive alignment. By default T-Coffee outputs this file in the directory where the process is running. If the T-Coffee process does not have permission to write in that directory, the computation will abort...To avoid this, simply specify the name of the output tree (refers to the **T-Coffee Parameters & Flags**, in the **Output** subsection about trees). Choose the name so that two processes may not over-write each other dnd file.


Permissions
-----------
The T-Coffee process MUST be allowed to write in some scratch area, even when it is ran by Mr nobody... Make sure the /tmp/ partition is not protected.


Other Programs
--------------
T-Coffee may call various programs while it runs (lalign2list by defaults). Make sure your process knows where to find these executables.


Known Problems
==============
 1) Sensitivity to sequence order: it is difficult to implement a MSA algorithm totally insensitive to the order of input of the sequences. In T-Coffee, robustness is increased by sorting the sequences alphabetically before aligning them. Beware that this can result in confusing output where sequences with similar name are unexpectedly close to one another in the final alignment.

 2) Nucleotides sequences with long stretches of Ns will cause problems to lalign, especially when using Mocca. To avoid any problem, filter out these nucleotides before running mocca.

 3) Stop codons are sometimes coded with \* in protein sequences, this will cause the program to crash or hang. Please replace the all \* signs with an X.

 4) Results can differ from one architecture to another, due rounding differences. This is caused by the tree estimation procedcure. If you want to make sure an alignment is reproducible, you should keep the associated dendrogram.


Development
===========
The following examples are only meant for internal development and are used to insure stability from release to release.

- **profile to list**

prf1: profile containing one structure
prf2: profile containing one structure

::


  $$: t_coffee -profile sample_profile1.aln,sample_profile2.aln -mode=3dcoffee


- **command line list**

These command lines have been checked before every release along with the other command lines in this documentation.

1) External methods:

::

  $$: t_coffee sample_seq1.fasta -in=Mclustalw_pair,Mclustalw_msa,Mslow_pair -outfile=text



2) List of command lines provided by James Watson to crash T-Coffee before version 3.40:

::

  $$: t_coffee -mode 3dcoffee -in Sthree_pdb.fasta P1PPG
  $$: t_coffee -mode 3dcoffee -in Sthree_pdb.fasta -template_file SELF_P_



3) Command line to read 'relaxed' PDB files...

::

  $$: t_coffee -in Msap_pair Ssample_3Dseq1.fasta -template_file sample_3Dseq1.template \
      -weight 1000 -out_lib sample_3Dseq1.tc_lib -lib_only


4) Parsing long sequence lines:

::

  $$: t_coffee -in Aproteases_large.aln -outfile test.aln


*************
T-Coffee Test
*************
Before the realease of a new T-Coffee version, many command lines in the documentation are programmatically tested. The purpose of this document is to explain which command lines are tested, how new command lines can be added to the documentation and how to perform a documentation check before a new release. 

Introducing a new command to test
=================================
In this documentation you will find three types of command lines some of which are up and running on the examples files. Whenever adding a new command, input files must be added to the repository directory ./examples/. New commands can be also be built using the existing files, or they can depend on files newly added to the repository. Here are the three different flag types of command lines you will fine in the documentation "$$", "$#" and "##":

- **$$: t_coffee...**
Corresponds to all command that will be tested systematically in a programmatic manner. Other command lines starting with different symbols are not checked.

- **$#: t_coffee...**
Corresponds to a command that could be tested but will not be, either for the sake of time or because it is currently unstable. When a new release needs to be urgently made available because of a critical fix, it is advisable to comment out this way non critical command lines failing the test.

- **##: t_coffee...**
These commands are never tested, either because they contain system dependant information or non programmatic information.




Checking the documentation
==========================
Our procedure uses **doc2test.pl** to extract all command lines from the .rst in docs and will run these command lines. When a reference output is available in /testsuite/docs/ref/ it will produce an other output and compare it with the reference. The new output will be put in /testsuite/docs/latest/. The script will eventually produce a global report that is to be found in /testsuite/docs/log/validation.log. The output of the failed scripts will be put in /testsuite/docs/failed. All the reference output will be added to git repository. The simplest way to check the documentation is to run the following command:

::

  Check the documentation:
  ##: ./lib/perl/lib/perl4makefile/doc2test.pl -mode update


When changes are made in the list of command line, we use the following mode  whichwill only run the new command lines and report (defined as those not having a reference output in /testsuite/docs/ref/). This mode is the default mode:

::

  Check new command lines:
  ##: ./lib/perl/lib/perl4makefile/doc2test.pl -mode new



This mode will only run the command lines that have previously failed, as indicated by the presence of an output in the failed directory.

::
  
  Check failed command lines:
  ##: ./lib/perl/lib/perl4makefile/doc2test.pl -mode failed



Compiling the documentation
===========================
The documentation can be compiled using the following command; it will cause the test to stop whenever a failed is encountered.

::

  Compile the documentation:
  ##: ./lib/perl/lib/perl4makefile/doc2test.pl -mode update -stop_on_failed


To have a clean start, you can also reset all files. This command will cause all the reference files to be erased:

::

  Reset reference files:
  ##: ./lib/perl/lib/perl4makefile/doc2test.pl -mode reset


To removed unnecessary files, the following command will cause all the files in ./examples/ that are not used for any /docs/*.rst commands to be deleted. They will be deleted from both their current locations and the git repository. The remaining files will be added to the git repository.

::

  Remove unnecessary files:
  ##: ./lib/perl/lib/perl4makefile/doc2test.pl -clean


Making a new release
====================
The Beta and Stable releases can be done with the makefile, from the t_coffee/src/ directory

::

  1-git commit -a -m 'what you did'
  2-make beta_release OR make stable_release -- This will populate the source directories
  3-git push


