###########################
T-Coffee Main Documentation
###########################

.. warning:: This chapter has been extensively updated in 11/2016. We improved its readibility and incorporate all the latest development of the past years (PSI/TM-Coffee for homology extension based MSAs, Pro-Coffee for functional DNA, R-Coffee/SARA-Coffee for RNA, TCS and STRIKE for evaluating MSAs, T-RMSD for structural clustering). If you have any suggestion or correction, don't hesitate to contact us.


*******************
Before You Start...
*******************
Foreword
========
Most of the work presented here emanates from two summer schools that were tentatively called the 'Prosite Workshops' and were held in Marseille, in 2001 and 2002. These workshops were mostly an excuse to go rambling and swimming in the creeks of Cassis (Calanques). Yet, when we got tired of lazing in the sun, we eventually did a bit of work to chill out. Most of our experiments were revolving around the development of sequence analysis tools. Many of the most advanced ideas in T-Coffee were launched during these fruitful sessions. Participants included Phillip Bucher, Laurent Falquet, Marco Pagni, Alexandre Gattiker, Nicolas Hulo, Christian Siegfried, Anne-Lise Veuthey, Virginie Leseau, Lorenzo Ceruti and Cedric Notredame.


Prerequisite for using T-Coffee
===============================
This documentation relies on the assumption that you have installed T-Coffee, version 9.03 or higher (preferably the latest stable version). T-Coffee is a freeware open source running on all Unix-like platforms including Mac OS X, woever, it cannot run on Windows except by using Cygwin, a freeware open source allowing to run a Unix-like command line on Windows (`download <https://www.cygwin.com/>`_). A better option for Windows users would be to install a Unix-like virtualbox on your computer. All the relevant information for installing T-Coffee is contained in the previous chapter **T-Coffee Installation**.


Let's have a try...
===================
This manual is made to help you discover (nearly) all subtleties of T-Coffee, from default applications to sophisticated ones. All along this documentation, we expect you to use Unix-like shell commands running on prepared example files we have created. For now you can find all these files on our github repository for now: `example files <https://github.com/cbcrg/tcoffee/tree/master/examples>`_; however in the future they will be included in the T-Coffee distribution (~/tcoffee/Version.../examples/). All T-Coffee commands indicated with a tag "$$:" are currently up and running; just copy/paste these commands from the documentation (without the tag "$$:") to your terminal provided you previously downloaded the example files. We encourage you to try all these examples which are controlled cases, and also to try these commands with your own data files. All T-Coffee commands indicated with a tag "$#:" are currently under maintenance. Also, be aware that all the values/results given in this manual were obtained with different versions of T-Coffee; depending on the version of T-Coffee you are using, these results may slightly differ.

.. important:: Using T-Coffee package via "command lines" is the best/only way  **if you need to do something sophisticated and/or highly reproducible**, but if you don't want to bother with command line you can use our `web server <http://tcoffee.crg.cat/apps/tcoffee/index.html>`_ or different links on the `Cedric Notredame's lab homepage <http://www.tcoffee.org>`_. We have tried to put as many of T-Coffee functionalities on the webserver and we plan to incorporate even more in the future, but for now some options are not available or limited. 


*******************
What is  T-Coffee ?
*******************
What is T-Coffee?
=================
T-Coffee stands for **"Tree based Consistency Objective Function For alignment Evaluation"**. It is mainly/primarily a Multiple Sequence Alignment (MSA) method but it also provides a collection of useful/powerful tools presented in this manual. Before going deep into the core of the matter, here are a few words to briefly explain the things T-Coffee can do !!!

What does it do?
----------------
T-Coffee has two components allowing you to perform different tasks:
 - mainly it is a Multiple Sequence Alignment method. Given a dataset of sequences previously gathered using database search programs (like BLAST, Ensembl, etc...), T-Coffee will produce a Multiple Sequence Alignment (MSA) (refer to the section **Building Your Multiple Sequence Alignments** for more details).
 - also it is a collection of tools (usually called with the flag **-other_pg "tools"**) and third party software (cf. section **Integrating External Methods in T-Coffee**) able to perform a wide range of operations such as aligning accurately, reformatting data, evaluating results, comparing methods and many more...

What can it align?
------------------
T-Coffee will align nucleic acid (DNA and RNA) and protein sequences alike. T-Coffee is also able to use other type of information such as secondary/tertiary structure information (for protein or RNA sequences with a known/predicted structure), sequence profiles, trees...

.. Hint:: To have a rough idea, on an average computer, T-Coffee can easily align up to a 200 sequences, about 1000 amino acid long in about ~20min. It is better to avoid aligning more than 1000 sequences without using a cluster or T-Coffee fast modes !!!

How can I use it?
-----------------
T-Coffee is not an interactive program. It runs from your Unix or Linux command line and you must provide it with the correct parameters and syntax. If you do not like typing commands, you can still use the latest T-Coffee webserver `here <http://tcoffee.crg.cat/apps/tcoffee/index.html>`_.

.. Tip:: Installing and using T-Coffee requires a minimum acquaintance with the Linux/Unix operating system. If you feel this is beyond your computer skills, we suggest you use one of the available online servers.

Is there an online webserver?
-----------------------------
Yes, at `Cedric Notredame's lab homepage <http://www.tcoffee.org>`_ which contains all the necessary links to our latest  `webserver <http://tcoffee.crg.cat/apps/tcoffee/index.html>`_ !

Is T-Coffee different from ClustalW?
------------------------------------
According to several benchmarks, T-Coffee is on overall much more accurate than ClustalW, but this increase in accuracy comes at a price: **T-Coffee (default mode) is slower than ClustalW** (about N times for N Sequences). Still, if you want to align closely related sequences, **T-Coffee can be used in a fast mode, much faster and about as accurate than ClustalW** (cf. section **Building Multiple Sequence Alignment**). If you are familiar with ClustalW or if you run a ClustalW server, you will find that we have made some efforts to ensure as much compatibility as possible between ClustalW and T-Coffee. Whenever it was relevant, we have kept the flag names and the flag syntax of ClustalW. Yet, you will find that T-Coffee also has many extra possibilities...

Is T-Coffee very accurate?
--------------------------
T-Coffee belongs to the category of consistency based aligners which currently corresponds to the most accurate algorithms available (e.g. ProbCons, MSAprobs...). In addition, T-Coffee can combine (many) methods and therefore be as accurate (and hopefully more) as the methods it combines. For instance, the "accurate" mode of T-Coffee is very slow but also very accurate; on average this mode was shown to be 10 % more accurate than normal aligners on sequences less than 30% similar. If you need a very accurate alignment go to section **Building Multiple Sequence Alignment**.


What T-Coffee can and cannot do for you ...
===========================================
What T-Coffee can't do
----------------------
To be honest, a short answer will be that there is only one thing T-Coffee cannot do for you: **T-Coffee can NOT fetch sequences for you**. You must select the sequences you want to align beforehand and prepare your own dataset. We suggest you use any BLAST server and format your sequences in FASTA so that T-Coffee can use them easily. The  `ExPASy BLAST server <http://www.expasy.ch>`_ provides a nice interface for integrating database searches.

What T-Coffee can do
--------------------
T-Coffee is not only just an aligner program, it comes with multiple tools and third party software increasing the range of its possibilities; here is a non exhaustive list of tasks T-Coffee can perform:

 - **T-Coffee can compute (or at least try to compute!) accurate Multiple Sequence Alignments of DNA, RNA or Protein sequences**. Several modes and options are available and will be presented all along this manual. The default T-Coffee accepts any kind of sequence, although some modes are specific to a given type of sequence.

 - **T-Coffee can help you to reformat, trim, clean, cut, color your input (sequences, structures...) or output (alignments, trees...) data**; meaning that once you have your data and/or results ready, you can always modify them at will.

 - **T-Coffee allows you to combine results obtained with several alignment methods** (see the section **FAQ for T-Coffee** and **Building Multiple Sequence Alignment** for more details). T-Coffee can virtually combine all these MSAs you have to produce a new Multiple Sequence Alignment having the best agreement with all these methods you tried.

 - **One of the most important improvement of T-Coffee is to let you combine sequences and structures**, so that your alignments are of higher quality. You need to have the SAP package installed to fully benefit of this facility (or to use another structural alignment method). 

 - **T-Coffee allows you to extract a collection of repeats from a single sequence or a set of sequences** using MOCCA. In other words, if you know the coordinates of one copy of a repeat, you can extract all the other occurrences. MOCCA needs some time to compute a library and then prompt you with an interactive menu. You just have to follow the instructions.

 - **T-Coffee can be used to measure the reliability of your Multiple Sequence Alignment**. If you want to find out about that, read the section **FAQ for T-Coffee** or the **Technical Documentation** (-output flag). More details will be given anyway in this manual in the section **How Good Is Your Alignment?**.

 - **T-Coffee can be used to compare alternative alignment**; in case you generate several alignments of the same sequences, you can compare these alignments using the most common scores (Sum-of-Pairs or Column Score). In case you have reference alignments, you can directly benchmark your method by comparing your MSAs to your references.

And probably many more options we will discover together all along this manual !

.. warning:: Some options are carried out using the function "wget". If "wget" is not installed on your system, you can get it for free from `wget download <http://www.wget.org>`_. Just type **wget** to be sure it is installed.


How does T-Coffee alignment works?
==================================
If you only want to make a standard Multiple Sequence Alignment, you may skip these explanations. But if you want to do more sophisticated things, these few indications may help before you start reading the documentation and the different articles. 

When you run T-Coffee, the first thing it does is to compute a library. The library is a list of pairs of residues that could be aligned...it is like a christmas list: you can ask anything you fancy, it doesn't imply you will get it. Given a standard library, it is nearly impossible to have all the residues aligned at the same time because all the lines of the library may not agree. For instance:

::

  Line 1 says:
  Residue 1 of seq A with Residue 5 of seq B,
  ...
  Line 100 says:
  Residue 1 of seq A with Residue 29 of seq B,

Each of these constraints comes with a weight and in the end, the T-Coffee algorithm tries to generate the multiple alignment that contains constraints whose sum of weights yields the highest score. In other words, it tries to make happy as many constraints as possible (replace the word constraint with, friends, relatives, collaborators... and you will know exactly what we mean). You can generate this list of constraints the way you like. You may even provide it yourself, forcing important residues to be aligned by giving them high weights (see **FAQ for T-Coffee**). For your convenience, T-Coffee can generate (by default) its own list by making all the possible global pairwise alignments, and the 10 best local alignments associated with each pair of sequences. Each pair of residues observed aligned in these pairwise alignments becomes a line in the library.

.. note:: Be aware that nothing forces you to use a given library and that you could build it using other methods. In protein language, **T-Coffee is synonymous for freedom, the freedom of being aligned however you fancy** (I was probably a Tryptophan in some previous life).


*******************
Preparing Your Data
*******************
.. important:: About the syntax of T-Coffee command lines, just for you to know that it is quite flexible, for instance you can use any kind of separator you want (i.e. , ; <space> =). The syntax used in this document is meant to be consistent with that of ClustalW. However, in order to take advantage of the automatic filename completion provided by many shells, you can replace '=' and ',' with a space. Also, T-Coffee tools/modes are called using different flags to specify input/output files, parameters, modifiers, etc...Some flags are not always mandatory, however, if you use a correct/strict flag usage T-Coffee will always work fine; you just have some degrees of freedom ;-).

The reformatting utility: seq_reformat
======================================
General introduction
--------------------
Nothing is more frustrating than downloading important data and realizing you need to format it before using it. In general, you should avoid manual reformatting: it is by essence inconsistent and will get you into trouble. It will also get you depressed when you realize that you have spent the whole day adding carriage return to each line in your files. T-Coffee comes with several tools to reformat/trim/clean/select your input data but also your output results, especially a very powerful reformatting utility named **seq_reformat**. You can use **seq_reformat** by invoking the t_coffee shell:

::

  $$: t_coffee -other_pg seq_reformat

This will output the online flag usage of seq_reformat meaning a complete list of things seq_reformat can do for you. The seq_reformat is a reformatting utility so it recognizes automatically the most common formats (FASTA, Swiss-Prot,ClustalW, MSF, Phylip...). It reads the input file(s) via the **"-in"** and **"-in2"** flags and outputs in whatever specified format via the **"-output"** flag. In the meantime, you can use the flag **"-action"** to perform a wide range of modification on your data. In this section we give you quite a lot of different examples of you can do with **"-other_pg seq_reformat"**.

.. danger:: After the flag **-other_pg**, the common T-Coffee flags are not recognized anymore; it is like if you were using a different program.

.. tip:: When using T-Coffee seq_reformat, command line may become quite long...a practical way to handle this is to create in your .bashrc an alias to call directly seq_reformat. For example, write in your .bashrc: **alias reformat='t_coffee -other_pg seq_reformat'**. You can now call seq_reformat by tiping **reformat**.

Modification options
--------------------
In order to perform different modifications on your data (residues/sequences/columns...), the **seq_reformat** utility can be followed by the flag **-action** and one or several modifiers listed here (this list is not exhaustive):

:: 

  Options:
  - +upper		: to uppercase your residues
  - +lower		: to lowercase your residues
  - +switchcase		: to selectively toggle the case of your residues
  - +keep		: to only keep the residues within the range
  - +use_cons +keep	: to only keep the columns within the range
  - +remove		: to remove the residues within the range
  - +convert		: to only convert the residues within the range
  - +grep		: to select a given string of character
  - +rm_gap		: to remove columns containing gaps
  - etc...
 
Using a "cache" file
--------------------
What is a cache in T-Coffee?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Several options can be performed easily by using what we call a cache (or cache file). In T-Coffee, a cache is a file containing an alternate version of your alignment where each position of the alignment is replaced by an alternative coding scheme. For instance each residue can be replaced by a score previously evaluated: this score can be the T-Coffee CORE index (cf. section **Evaluating Your Alignment**) or a matrix based evalution (blosum62mt or identity matrix). Then, when performing any modification or reformatting of your alignments, you can just specify the range of positions to be modified according to their respective scores within the cache. We will see some example especially regarding the modification of format of a given alignment; it is not mandatory to use a cache but it is quite practical. To generate a cache before any reformatting using a given evaluation score, you can use one of the following possible option:

::

  Evaluating the T-Coffee CORE index during the alignment procedure:
  $$: t_coffee sample_seq1.fasta -output=score_ascii

  Evaluating the T-Coffee CORE index of a given alignment (-infile is mandatory):
  $$: t_coffee -infile sample_seq1.aln -mode evaluate

  Using an identity matrix:
  $$: t_coffee -other_pg seq_reformat -in sample_seq1.aln -action +evaluate \
      idmat -output=score_ascii

  Using a substitution matrix:
  $$: t_coffee -other_pg seq_reformat -in sample_seq1.aln -action +evaluate \
      blosum62mt -output=score_ascii

Preparing a sequence/alignment cache
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The following command will convert your alignment according to the given parameters: all A/a will be turned into 1, the gaps "-" will be conserved and all the other symbols (#) will be turned into 0. The flag **-action +convert** indicates the actions that must be carried out on the alignment before it is output into cache.

::

  1) Initial alignment:
  CLUSTAL FORMAT
  A CTCCGTgTCTAGGagt-TTACGTggAGT
  B CTGAGA----AGCCGCCTGAGGTCG---
  D CTTCGT----AGTCGT-TTAAGAca---
  C -TTAAGGTCC---AGATTGCGGAGC---


  2) Command line:
  $$: t_coffee -other_pg seq_reformat -in=sample_dnaseq3.aln -output=clustalw_aln -out=cache.aln \
      -action +convert 'Aa1' '.-' '#0'

  3) Cache generated:
  CLUSTAL W (1.83) multiple sequence alignment
  A 0000000000100100-00100000100
  B 000101----100000000100000---
  D 000000----100000-00110101---
  C -001100000---101000000100---


Other alternative are possible, for instance generating a cache for unaligned sequences (**-in** can refer to an alignment or unaligned sequences):

::

  1) Command line:
  $$: t_coffee -other_pg seq_reformat -in=sample_dnaseq3.aln -output=fasta_seq -out=cache.seq \
      -action +convert 'Aa1' '.-' '#0'

  2) Generate the file cache.seq:
  >A
  0000000000100100000100000100
  >B
  0001010000100000000100000000
  >D
  0000000000100000000110101000
  >C
  0001100000000101000000100000


Preparing a library cache
^^^^^^^^^^^^^^^^^^^^^^^^^
The library is a special format used by T-Coffee to declare relationships between pairs of residues. The cache library format can also be used to declare for instance the color of specific residues in an alignment. For instance, the following file
``sample_dnaseq3.tc_lib`` declares which residue of which sequence will receive which color. Note that the sequence number and the residue index are duplicated, owing to the recycling of this format from its original usage. It is also possible to use the BLOCK operator when defining the library (see **Technical Documentation**). The number right after BLOCK indicates the block length (10). The two next numbers (1 1) indicate the position of the first element in the block. The last value is the color.

::

  ! TC_LIB_FORMAT_01
  4
  A 27 CTCCGTgTCTAGGagtTTACGTggAGT
  B 21 CTGAGAAGCCGCCTGAGGTCG
  C 21 TTAAGGTCCAGATTGCGGAGC
  D 20 CTTCGTAGTCGTTTAAGAca
  #1 1
   +BLOCK+ 10 1 1 3
   +BLOCK+ 5 15 15 5
  #3 3
   6 6 1
   9 9 4
  ! CPU 240
  ! SEQ_1_TO_N
      
Modifying the format of your data
=================================
Keeping/Protecting your sequence names
--------------------------------------
Only few programs support long sequence names, and sometimes, when going through some pipeline the names of your sequences can be truncated or modified. To avoid this, **seq_reformat** contains a utility that can automatically rename your sequences into a form that will be machine-friendly, while making it easy to return to the human-friendly form.

1) **Create a code list**: The first thing to do is to generate a list of names that will be used in place of the long original name of the sequences:

::

  $$: t_coffee -other_pg seq_reformat -in proteases_large.fasta -output \
      code_name > proteases_large.code_name

2) **Code your data**: This will create a file where each original name is associated with a coded name (Cxxx). You can then use this file to either code your dataset using the following command:

::

  $$: t_coffee -other_pg seq_reformat -code proteases_large.code_name -in \
      proteases_large.fasta > proteases_large.coded.fasta

3) **Decode your data**: Then you can work with the file sproteases_large.coded.fasta and when you are done, you can decode the names of your sequences with the following command line:

::

  $$: t_coffee -other_pg seq_reformat -decode proteases_large.code_name -in \
      proteases_large.coded.fasta

Changing the sequence format
----------------------------
Sometimes it may be necessary to change from one format to another, for instance when using another software which recognize only a given format. T-Coffee recognizes most common alignment formats and you can find the list of all input or output format recognized by simply typing **t_coffee -other_pg seq_reformat** without any input file(s). It is possible to reformat unaligned or aligned sequences alike although changing the alignment format is probably more interesting in order to use other applications; unaligned sequences format flags are generally preceded by the suffix **"_seq"** and aligned sequences flags by the suffix **"_aln"**. This also allows you to transform any alignment into unaligned sequences by removing the gaps. Here are some examples on how to change the format of your data:

::

  For unaligned sequences (e.g. FASTA to PIR):
  $$: t_coffee -other_pg seq_reformat -in proteases_small.fasta -output pir_seq > \
      proteases_small.pir_seq
  
  For alignements (e.g. ClustalW to MSF):
  $$: t_coffee -other_pg seq_reformat -in proteases_small.aln -output fasta_aln > \
      proteases_small.fasta_aln
      
  From aligned to unaligned sequences:
  $$: t_coffee -other_pg seq_reformat -in proteases_small.aln -output fasta_seq > \
      proteases_small.fasta

.. Warning:: Format recognition is not 100% full proof; occasionally you will have to inform the program about the nature of the file you are trying to reformat with **-input msf_aln -output fasta_aln** for instance.

Changing the case
-----------------
Changing the case of your sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you need to change the case of your sequences, you can use different modifiers embedded in **seq_reformat**. They are accessed via the **-action** flag. For instance, to write your sequences in lower case:

::

  $$: t_coffee -other_pg seq_reformat -in proteases_small.aln -action +lower \
      -output clustalw

.. tip:: No prize for guessing that **+upper** will do exactly the opposite...

Changing the case of specific residues
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you want to change the case of a specific residue, you can use the flag: **+edit_residue <sequence> <residue #> <lower|upper|symbol>**. If you have more than one residue to modify, write all the coordinates in a text file (one coordinate per line) as spans are not yet supported; then give the file to T-Coffee

::

  $$: t_coffee -other_pg seq_reformat -in sample_seq1.aln -action +upper \
      +edit_residue hmgb_chite 10 lower
      

  $$: t_coffee -other_pg seq_reformat -in sample_seq1.aln -action +upper \ 
      +edit_residue sample_seq1.list

.. warning:: If you give a list of coordinates, it has to be a Unix text file (not a word document).

Changing the case with a cache
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you want to change the case depending on the score, you must either evaluate your alignment, or provide a cache. For example, this command line will upper the case of all residue then lower the case of every residue more than 50% identical to other residues in the same column:

::

  Using a cache on the fly:
  $$: t_coffee -other_pg seq_reformat -in sample_dnaseq2.aln -action +upper \
      +evaluate idmat +lower '[5-9]'
      
  Using a cache file previously computed (2 steps):
  $$: t_coffee -other_pg seq_reformat -in sample_dnaseq2.aln -action +evaluate \
      idmat -output=score_ascii > sample_dnaseq2.cache
      
  $$: t_coffee -other_pg seq_reformat -in sample_dnaseq2.aln -struc_in sample_dnaseq2.cache \
     -struc_in_f number_aln -action +lower '[5-9]'
  
Coloring/Editing residues in an alignment
-----------------------------------------
Changing the default colors
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Colors are hard coded in the program, but if you wish, you can change them by simply creating a file named ``seq_reformat.color`` that is used to declare the color values. The name of the file (``seq_reformat.color``) is defined in ``programmes_define.h``, COLOR_FILE and can be changed before compilation. By default, the file is searched in the current directory. For example, the following line written in ``seq_reformat.color`` indicates that the value 0 in the cache corresponds now to #FFAA00 in html, and in RGB 1, 0.2 and 0. 

::

  Format for hard coded colors in T-Coffee:
  0 #FFAA00 1 0.2 0

Coloring specific types of residues/nucleic acids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
You can color all the residues of your sequences on the fly; for instance, the following command line will color all the a's in color 0 (blue):

::

  $$: t_coffee -other_pg seq_reformat -in sample_seq1.aln -action +3convert a0 \
      -output color_html > sample_seq1_color.html

.. warning:: This option is case sensitive so the case of the residues or nucleotides should be the same in the command line (in this command line, only a lower case will be colored). 

Coloring a specific residue of a specific sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you want to color a specific residue/nucleotide, you can use the flag **+color_residue <sequence> <residue #> <color #>**. If you have more than one residue to color, you can put all the coordinates in a file, (one coordinate per line). Spans are not yet supported.

::

  $$: t_coffee -other_pg seq_reformat -in sample_seq1.aln -action +color_residue \
      hmgb_chite 10 1 -output color_html > sample_seq1_single.html


  $$: t_coffee -other_pg seq_reformat -in sample_seq1.aln -action +color_residue \
      sample_seq1_color.list -output color_html > sample_seq1_all.html
      
.. warning:: If you give a list of coordinates, it has to be a Unix text file (not a word document).

Coloring according to the conservation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Use the +evaluate flag if you want to color your alignment according to its conservation level or using the boxshade scoring scheme:

::

  Conservation color scheme:
  $$: t_coffee -other_pg seq_reformat -in sample_seq1.aln -action +3evaluate pam250mt \
      -output color_html > color_cons.html

  Boxshade color scheme:
  $$: t_coffee -other_pg seq_reformat -in sample_aln1.aln -action +3evaluate boxshade \
      -output color_html > color_cons_box.html

Coloring an alignment using a cache
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you have a cache alignment or a cache library, you can use it to color your alignment and either make a post script, html or PDF output. For instance, if you use the file cache.seq:

::

   Produces a html file:
   $$: t_coffee -other_pg seq_reformat -in=sample_dnaseq3.aln -struc_in=sample_dnaseq3.cache \
       -struc_in_f number_fasta -output=color_html -out=color_dnaseq3.html

  Produces a pdf file:
   $*: t_coffee -other_pg seq_reformat -in=sample_dnaseq3.aln -struc_in=sample_dnaseq3.cache \
       -struc_in_f number_fasta -output=color_pdf -out=color_dnaseq3.pdf

  Produces an output using a library:
  $$: t_coffee -other_pg seq_reformat -in=sample_dnaseq3.aln -struc_in=sample_dnaseq3.tc_lib \
      -output=color_html -out=color_dnaseq3_lib.html
      
.. warning:: The script **ps2pdf** must be installed on your system for the pdf options.

Modifying the data itself...
============================
Modifiying sequences in your dataset
------------------------------------
Converting residues
^^^^^^^^^^^^^^^^^^^
It is possible for instance to selectively convert all given characters in a sequence (residues or nucleic acids alike) into another one, for example all G's having a score between 5 and 9 by using the command line:

::

  $$: t_coffee -other_pg seq_reformat -in sample_dnaseq2.aln -struc_in sample_dnaseq2.cache \ 
      -struc_in_f number_aln -action +convert '[5-9]' GX

 
Extracting sequences according to a pattern
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
You can extract any sequence by requesting a specific pattern to be found either in the name (NAME), the comment (COMMENT) or the sequence (SEQ) using the modifier is '+grep'. For instance, if you want to extract all the sequences whose name contain the word HUMAN, the flag NAME/COMMENT/SEQ indicates that the modification is made according to the sequences names, the comment section or the sequence itself, and the flag KEEP/REMOVE means that you will keep/remove all the sequences containing the string HUMAN. Here are some examples:

::

  To keep sequences containing HUMAN in the name:
  $$: t_coffee -other_pg seq_reformat -in proteases_large.aln -action +grep NAME \
      KEEP HUMAN -output clustalw

  To remove sequences containing HUMAN in the name:
  $$: t_coffee -other_pg seq_reformat -in proteases_large.aln -action +grep NAME \
      REMOVE HUMAN -output clustalw

  To keep sequence which contain sapiens in the comment:
  $$: t_coffee -other_pg seq_reformat -in proteases_large.fasta -action +grep COMMENT \
      KEEP sapiens -output clustalw
 
  To remove sequences containing the pattern [ILM]K:
  $$: t_coffee -other_pg seq_reformat -in proteases_large.aln -action +grep SEQ \
      REMOVE '[ILM]K' -output clustalw

.. important:: You should know that the pattern can be any perl legal regular expression, you can visit this  `page <http://www.comp.leeds.ac.uk/Perl/matching.html>`_ for some background on regular expressions. 

.. caution:: This option is case sensitive (Human, HUMAN and hUman will not yield the same results). Be careful !!!

Extracting/Removing specific sequences by names
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you want to extract (command 1) or remove (command 2) several sequences in order to make a subset, you can specify a list of sequences by their full name:

::

  Command 1: keep sequences
  $$: t_coffee -other_pg seq_reformat -in proteases_small.fasta -action +extract_seq_list \
      'sp|P29786|TRY3_AEDAE' 'sp|P35037|TRY3_ANOGA'

  Command 2: remove sequences
  $$: t_coffee -other_pg seq_reformat -in proteases_small.fasta -action +remove_seq \
      'sp|P29786|TRY3_AEDAE' 'sp|P35037|TRY3_ANOGA'

.. note:: Note the single quotes (') are mandatory as they are meant to protect the name of your sequence and prevent the Unix shell to interpret it like an instruction.

Once sequences are extracted or removed, some columns may remain containing only gaps, but it is possible to simply remove empty columns from the resulting dataset (command 3), and even extract specific blocks for the selected sequences either keeping the exact same name (command 4) or the name of the specific blocks extracted (command 5):

::

  Command 3: removing empty columns
  $$: t_coffee -other_pg seq_reformat -in proteases_small.aln -action +extract_seq_list \
      'sp|P29786|TRY3_AEDAE' 'sp|P35037|TRY3_ANOGA' +rm_gap

  Command 4: keeping the initial name after extracting specific blocks and removing empty columns
  $$: t_coffee -other_pg seq_reformat -in proteases_small.aln -action +keep_name \
      +extract_seq 'sp|P29786|TRY3_AEDAE' 20 200 'sp|P35037|TRY3_ANOGA' 10 150 +rm_gap

  Command 5: renaming sequences according to the extracted blocks and removing empty columns
  $$: t_coffee -other_pg seq_reformat -in proteases_small.aln -action +extract_seq \
      'sp|P29786|TRY3_AEDAE' 20 200 'sp|P35037|TRY3_ANOGA' 10 150 +rm_gap 

.. hint:: The tag **+keep_name** must come BEFORE the tag **+extract_seq**.

Extracting the most informative sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Large datasets are problematic because they can be difficult to align and analyze, MSA programs tend to become very slow and inaccurate. In short, the best size for an MSA dataset would be between 20 to 40 sequences to have enough sequences to see the effect of evolution, but in the same time small enough so that you can visualize your alignment and recompute it as many times as needed. More important than its size, a good dataset have to be informative, when each sequence contains information the others do not have. The most informative sequences are the sequences that are as different as possible to one another, within your dataset. You can extract the most informative sequences using flag **+trim** followed by the number of sequences you wish to keep ("n" for a number and "N" for a perrcentage). The following commands will extract the 10 most informative sequences (command 1) or the 20% of most informative sequences (command 2):

::

  Command 1:
  $$: t_coffee -other_pg seq_reformat -in proteases_large.fasta -action +trim _seq_n10 \
      -output fasta_seq
  Command 2:
  $$: t_coffee -other_pg seq_reformat -in proteases_large.fasta -action +trim _seq_N20 \
      -output fasta_seq

.. hint:: The argument to trim include **_seq_**, it means your sequences are provided unaligned. If your sequences are already aligned, you do not need to provide this parameter. It is generally more accurate to use unaligned sequences.

.. note:: For very large dataset, **seq_reformat** will compute the similarity matrix between your sequences once only. It will then store it in its cache to be reused any time you run on the same dataset. In short this means that it will take much longer to run the first time, but be much faster if you need to rerun it.

Extracting/Removing sequences with the % identity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
**Removing too identical sequences (redundant)**

Removing the most similar sequences is often what people have in mind when they talk about removing redundancy. You can do so using the **+trim** option. For instance, you can generate a dataset where no pair of sequences has more than 50% identity either from a dataset of unaligned sequences (command 1) or from any given alignment (command 2). If you start from unaligned sequences, the removal of redundancy can be slow. If your sequences have already been aligned using a fast method, you can take advantage of this by replacing the **"_seq_"** with **"_aln_"**. Just run the following command lines to see the difference un runtime:

::

  Command 1: unaligned sequences
  $$: t_coffee -other_pg seq_reformat -in proteases_large.fasta -action +trim _seq_%%50_

  Command 2: aligned sequences
  $$: t_coffee -other_pg seq_reformat -in proteases_large.fasta -action +trim _aln_%%50_

.. note:: Using aligned sequences results in a fastest trimming, however, it also means that you rely on a more approximate estimation of sequence similarity.

**Removing too different sequences (outliers)**

Sequences that are too distantly related from the rest of the set (called outliers) may have very negative effects on the overall alignment; to prevent this, it is advisable not to use them. The next command line will lead to the removal of all the sequences where no pair of sequences has less than 30% average accuracy with all the other sequences in the dataset (the symbol "_O" stands for Outliers) and more than 80% identity: 

::

  $$: t_coffee -other_pg seq_reformat -in proteases_large.fasta -action +trim _seq_%%80_O30

.. hint:: This particular option is quite powerful as it allows you to decide both inferior and superior tresholds for trimming your dataset based on pairwise identity score, and therefore you can dissect your dataset according to different ranges of identity values. Be careful not to remove too many sequences ;-)

**Forcing specific sequences to be kept**

Sometimes you want to trim based on identity while making sure specific/important sequences remain in your dataset. You can do so by providing a pattern (**"_f"** for field) : it will keep all the sequences whose name contains the given string (**"_fNAME"**, **"_fCOMMENT"** or **"_fSEQ"**). Here are some examples corresponding to the different protected fields while removing all sequences above 50% identity: 

::

  Keep all HUMAN sequences    
  $$: t_coffee -other_pg seq_reformat -in proteases_large.fasta -action +trim \
      _seq_%%50_fNAME HUMAN

  Keep all sequences containing ".apiens"
  $$: t_coffee -other_pg seq_reformat -in proteases_large.fasta -action +trim \
      _seq_%%50_fCOMMENT '.apiens'

  Keep all sequences containing residues
  $$: t_coffee -other_pg seq_reformat -in proteases_large.fasta -action +trim \
      _seq_%%50_fSEQ '[MLV][RK]'

You can also specify the sequences you want to keep by giving another fasta file containing the name of these sequences via the flag **-in2**:

::

  $$: t_coffee -other_pg seq_reformat -in proteases_large.fasta -in2 proteases_small.fasta \
      -action +trim _seq_%%40

Chaining important sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In order to align two distantly related sequences, most multiple sequence alignment packages perform better when provided with many intermediate sequences that make it possible to 'bridge' your two sequences. The modifier **+chain** makes it possible to extract from a dataset a subset of intermediate sequences that chain the sequences you are interested in. For instance, let us consider the two sequences "sp|P21844|MCPT5_MOUSE" and "sp|P29786|TRY3_AEDAE" having 26% identity. This is high enough to make a case for a homology relationship between them, but this is too low to blindly trust any pairwise alignment. With the names of the two sequences written in the file sproteases_pair.fasta, run the following command:

::

  $$: t_coffee -other_pg seq_reformat -in proteases_large.fasta -in2 proteases_pair.fasta \
      -action +chain > proteases_chain.fasta

This will generate a dataset of 21 sequences, with the following chain of similarity between your two sequences. This is probably the best way to generate a high quality alignment of your two sequences when using a progressive method like ClustalW, T-Coffee, MUSCLE or MAFFT.

::

  N: 21 Lower: 40 Sim: 25 DELTA: 15
  #sp|P21844|MCPT5_MOUSE -->93 -->sp|P50339|MCPT3_RAT -->85 -->sp|P50341|MCPT2_M\
  ERUN -->72 -->sp|P52195|MCPT1_PAPHA -->98 -->sp|P56435|MCPT1_MACFA -->97 -->sp\
  |P23946|MCPT1_HUMAN -->81 -->sp|P21842|MCPT1_CANFA -->77 -->sp|P79204|MCPT2_SH\
  EEP -->60 -->sp|P21812|MCPT4_MOUSE -->90 -->sp|P09650|MCPT1_RAT -->83 -->sp|P5\
  0340|MCPT1_MERUN -->73 -->sp|P11034|MCPT1_MOUSE-->76 -->sp|P00770|MCPT2_RAT --\
  >71 -->sp|P97592|MCPT4_RAT -->66 -->sp|Q00356|MCPTX_MOUSE -->97 -->sp|O35164|M\
  CPT9_MOUSE -->61 -->sp|P15119|MCPT2_MOUSE -->50 -->sp|Q06606|GRZ2_RAT -->54 --\
  >sp|P80931|MCT1A_SHEEP -->40 -->sp|Q90629|TRY3_CHICK -->41 -->sp|P29786|TRY3_A\
  EDAE

Modifying columns/blocks in your dataset
----------------------------------------
Removing gapped columns
^^^^^^^^^^^^^^^^^^^^^^^
You can also remove all the columns containing a given proportion of gaps; for instance the following command will delete all the residues occurring in a column that contains 50% or more gaps (use 1 to delete residues from columns having 1 gap or more):

::

  $$: t_coffee -other_pg seq_reformat -in sample_dnaseq3.aln -action +rm_gap 50

Extracting specific columns 
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Extracting portions of a dataset is something very frequently needed. You may need to extract all the sequences that contain the word human in their name, or you may want all the sequences containing a simple motif. We show you here how to do a couple of these things. To do this, you need an evaluation file that may have been generated with T-Coffee, either running a *de novo* alignment (command 1) or evaluating a preexisting alignment (command 2):

::

  Command 1:
  $$: t_coffee sample_seq1.fasta -output score_ascii, aln
  
  Command 2:
  $$: t_coffee -other_pg seq_reformat -in sample_seq1.aln -action +evaluate \
      blosum62mt -output score_ascii > sample_seq1_blosum62.score_ascii

This generates a score_ascii file that you can then use to filter out the bad bits in your alignment considering the individual score of each residue to trigger the filtering (command 3), or according to the whole column score by simply add the **+use_cons** flag (command 4). This last command can also be run on the fly with command 5. The commands 3/4/5 will keep only residues and/or columns having a score between 6 and 9.

::

  Command 3:
  $$: t_coffee -other_pg seq_reformat -in sample_seq1.aln -struc_in sample_seq1.score_ascii \
      -struc_in_f number_aln -action +keep '[6-9]'

  Command 4:
  $$: t_coffee -other_pg seq_reformat -in sample_seq1.aln -struc_in sample_seq1.score_ascii \
      -struc_in_f number_aln -action +use_cons +keep '[6-9]'

  Command 5
  $$: t_coffee -other_pg seq_reformat -in sample_aln1.aln -action +evaluate blosum62mt \
       +use_cons +keep '[6-9]'

.. warning:: Don't forget the simple quotes ('), it's mandatory !!!

Extracting entire blocks
^^^^^^^^^^^^^^^^^^^^^^^^
In case you want to extract a specific block of your alignment, for instance to remove poorly resolved regions, delimit your alignments boundaries or to extract domains, you can do so with the option **+extract_block**. In command 1, the option **cons** indicates that you are counting the positions according to the consensus of the alignment (i.e. the positions correspond to the columns # of the alignment). If you want to extract your block relatively to a specific sequence, you should replace cons with this sequence name (command 2).

::

  Command 1: extract block from MSA
  $$: t_coffee -other_pg seq_reformat -in proteases_small.aln -action +extract_block \
      cons 150 200

  Command 2: extract_block relative to a give sequence of the MSA
  $$: t_coffee -other_pg seq_reformat -in proteases_small.aln -action +extract_block \
      'sp|Q03238|GRAM_RAT' 10 200

.. tip:: It may be sometimes difficult to know where starts the blocks you are interested in except by counting manually the number of column. You can also make some tries by modifying the boundaries until you get the block you want and then redirect the result into the output file name of your choice. 

Concatenating blocks or MSAs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you have extracted several blocks generated using the previous command and you want to glue them together, you can use the **+cat_aln** modifier:

::

  $$: t_coffee -other_pg seq_reformat -in proteases_small.aln -action +extract_block \
      cons 100 120 > block1.aln

  $$: t_coffee -other_pg seq_reformat -in proteases_small.aln -action +extract_block \
      cons 150 200 > block2.aln

  $$: t_coffee -other_pg seq_reformat -in block1.aln -in2 block2.aln -action +cat_aln

.. note:: The alignments do not need to have the same number of sequences and the sequences do not need to come in the same order.


Manipulating DNA sequences
==========================
Translating DNA sequences into protein sequences
------------------------------------------------
If your sequences are DNA coding sequences, it is often safer and more accurate to align them as proteins (as protein sequences are more conserved than their corresponding DNA sequence). The **seq_reformat** options make it easy for you to translate your sequences:

::

  $$: t_coffee -other_pg seq_reformat -in proteases_small_dna.fasta -action \
      +translate -output fasta_seq

Back-translation with the *bona fide* DNA sequences
---------------------------------------------------
Once your sequences have been aligned, you may want to turn your protein alignment back into a DNA alignment, either to do phylogeny, or maybe in order to design PCR probes. To do so, use the following command:

::

  $$: t_coffee -other_pg seq_reformat -in proteases_small_dna.fasta -in2 \
      proteases_small.aln -action +thread_dna_on_prot_aln -output clustalw

Finding the *bona fide* sequences for the back-translation
----------------------------------------------------------
Use the online server `ProtoGen <http://tcoffee.vital-it.ch/apps/tcoffee/do:protogene>`_.


Manipulating RNA Sequences 
==========================
Producing a Stockholm output: adding predicted secondary structures
-------------------------------------------------------------------
Producing/Adding a consensus structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Given an RNA multiple sequence alignment, it is possible to compute (command 1) or add (command 2) the alifold (Vienna package) consensus secondary structure and output in stockholm:

::

  Command 1:
  $$: t_coffee -other_pg seq_reformat -in sample_rnaseq2.aln -action +aln2alifold \
      -output stockholm_aln

  Command 2: 
  $$: t_coffee -other_pg seq_reformat -in sample_rnaseq2.aln -action +add_alifold \
      -output stockholm_aln

Adding a precomputed consensus structure to an alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The file sample_rnaseq2.alifold contains the raw output of the alifold program produced via the RNAalifold `webserver <http://rna.tbi.univie.ac.at/cgi-bin/RNAalifold.cgi>`_ or captured with the command "RNAalifold <sample_rnaseq2.aln > sample_rnaseq2.alifold". It is possible to add this secondary structure to an alignment (command 1) and to stack Stockholm formatted secondary structures (command 2):

::

  Command 1:
  $$: t_coffee -other_pg seq_reformat -in sample_rnaseq2.aln -in2 sample_rnaseq2.alifold \ 
      -input2 alifold -action +add_alifold -output stockholm_aln  

  Command 2:
  $$: t_coffee -other_pg seq_reformat -in sample_rnaseq2.aln -in2 sample_rnaseq2.cons.stk \
      -action +add_alifold -output stockholm_aln

.. warning:: The alifold structure and the alignment MUST be compatible. The function makes no attempt to thread or align the structure, it merely stacks it below the MSA.

Analyzing a RNAalifold secondary structure prediction
-----------------------------------------------------
The following commands can either be applied on a Stockholm or a standard MSA. In the second case (standard MSA) the secondary structure will be automatically recomputed by alifold.

Visualizing compensatory mutations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The following command will output a color coded version of your alignment with matching columns indicated as follows:

 - I: incompatible pair (i.e. at least one pair is not WC)
 - N: pairs are Gus or WC
 - W: all pairs are Watson
 - c: compensatory mutations
 - C: WC compensatory mutations

::

  Standard alignment:
  $$: t_coffee -other_pg seq_reformat -in sample_rnaseq2.aln -action +alifold2analyze aln
  
  Color coded alignment:
  $$: t_coffee -other_pg seq_reformat -in sample_rnaseq2.aln -action +alifold2analyze color_html

.. warning:: Handling gapped columns: by default gapped column are ignored but they can be included by adding the tag **-usegap**.

Analyzing matching columns
^^^^^^^^^^^^^^^^^^^^^^^^^^
The option **+alifold2analyze** will estimate the number of pairs of columns that are perfect Watson and Crick pairings, those that are neutral (including a GU) and those that include correlated mutations (command 1). The WCcomp are the compensated mutations maintaining WC base pairing. Other arguments can given, to display the list of paired positions and their status (compensated, Watson, etc...) use command 2:

::

  Command 1:
  $$: t_coffee -other_pg seq_reformat -in sample_rnaseq2.stk -action +alifold2analyze stat
  
  Command 2:
  $$: t_coffee -other_pg seq_reformat -in sample_rnaseq2.stk -action +alifold2analyze list


Comparing alternative folds
---------------------------
The folds associated with alternative alignments can be compared. This comparison involves counting how many identical pairs of residues are predicted on each sequence in one fold and in the other. The top of the output (@@lines) summarizes the results that are displayed on the input alignment; if the provided alignment do not have a fold, this fold will be estimated with alifold. The folds can be provided as Stockholm alignments:

::

  $$: t_coffee -other_pg seq_reformat -in sample_rnaseq2.cw.stk -in2 sample_rnaseq2.tcoffee.stk \
      -action +RNAfold_cmp


Phylogenetic Trees Manipulation
===============================
Producing phylogenetic trees
----------------------------
The seq_reformat is NOT a phylogeny package, yet over the time it has accumulated a few functions that make it possible to compute simple phylogenetic trees, or similar types of clustering. Given a multiple sequence alignment, it is possible to compute either a UPGMA or an NJ tree. The following commands use an identity matrix to compare your sequences and will output an unrooted NJ tree in Newick format (command 1) or a rooted UPGMA tree (command 2):

::

  Command 1:
  $$: t_coffee -other_pg seq_reformat -in sample_seq1.aln -action +aln2tree -output newick \
      -out sample_seq1_tree_nj.nwk

  Command 2:
  $$: t_coffee -other_pg seq_reformat -in sample_seq1.aln -action +aln2tree _TMODE_upgma  \
      -output newick -out sample_seq1_tree_upgma.nwk

If your data is not data sequence, but a matrix of 1 and Os (i.e. SAR matrix for instance), you can use a different matrix to compute the pairwise distances (command 3), and all these parameters can be concatenated (command 4):

::

  Command 3:
  $#: t_coffee -other_pg seq_reformat -in sample_seq1.aln -action +aln2tree \
      _MATRIX_sarmat -output newick

  Command 4:
  $#: t_coffee -other_pg seq_reformat -in sample_seq1.aln -action +aln2tree \
      _TMODE_upgma_MATRIX_sarmat -output newick

.. warning:: Bootstrap facilities will also be added at some point...We recommend you to use `Phylip <http://evolution.genetics.washington.edu/phylip.html>`_ or any other specific phylogenetic software (PhyML, RAxML, MrBayes, etc...) if you need some serious phylogeny !

Comparing two phylogenetic trees
--------------------------------
A real interesting option is the ability to compare two trees (unrooted) returning some ofthe most common scores used for this including Robinson-Foulds

::

  $$: t_coffee -other_pg seq_reformat -in sample_tree1.dnd -in2 sample_tree2.dnd -action \
      +tree_cmp -output newick

  #tree_cmp|T: 33 W: 20.00 L: 14.88 RF: 2 N: 9 S: 5
  #tree_cmp_def|T: ratio of identical nodes
  #tree_cmp_def|W: ratio of identical nodes weighted with the min Nseq below node
  #tree_cmp_def|L: average branch length similarity
  #tree_cmp_def|RF: Robinson and Foulds
  #tree_cmp_def|N: number of Nodes in T1 [unrooted]
  #tree_cmp_def|S: number of Sequences in T1

About the output scores in more details:

 - T: Fraction of the branches conserved between the two trees. This is obtained by considering the split induced by each branch and by checking whether that split is found in both trees.
 - W: Fraction of the branches conserved between the two trees. Each branch is weighted with MIN, the minimum number of leaf on its left or right.
 - L: Fraction of branch length difference between the two considered trees.
 - RF: is the standard Robinson-Foulds value when comparing trees.

The last line contains a tree where distances have been replaced by the number of leaf under the considered node:

 - Positive values indicate a node common to both trees and correspond to MIN.
 - Negative values indicate a node found in tree1 but not in tree2.
 - The higher this value, the deeper the node.

.. tip:: You can extract this tree for further usage by typing **cat outfile | grep -v 'tree_cmp'**

Scanning phylogenetic trees
---------------------------
It is possible to scan an alignment and locally measure the similarity between an estimated local tree and some reference tree provided from an external source (or computed on the fly) using the following command:

::

  $# : t_coffee -other_pg seq_reformat -in sample_seq1.aln -in2 sample_seq1_tree.nwk -action \
       +tree_scan _MODE_scan__W_10_ > ph_tree_scan.txt

For each position of the alignment, W*2 blocks of size 2*1+1 up to W*2+1 will be extracted, for each of these block a tree will be estimated and the similarity of that tree with the reference tree will be estimated with cmp_tree. For each position, the tree giving the best fit will be reported, along with the size of the block leading to that tree:

::

  P: <position> <block start> <block_end> <block score> <block Length>

Pruning phylogenetic trees
--------------------------
Pruning removes leaves from an existing tree and recomputes distances so that no information is lost. To do this with T-Coffee you need two input files: a tree file in the Newick format and a FASTA-like file where sequences can be omitted, but you can also leave them, at your entire convenience. The second file is merely a list of the sequences to be kept when pruning the tree. The resulting tree will contain only the sequences specified in the list.

::

  Tree file: "sample_3Dseq1.tree"
  ((((TNFR10-2:0.04546,TNFR16-2:0.06640)...,TNFR4-2:0.05255) 45:0.00848); \
  
  FASTA-like file: "sample_3Dseq1.fake"
  >TNFR2-2
  >TNFR4-4
  ...

  Pruning the tree:
  $$: t_coffee -other_pg seq_reformat -in sample_3Dseq1.tree -in2 sample_3Dseq1.fasta -action \
      +tree_prune -output newick

Tree Reformatting for MAFFT
---------------------------
The most common format for phylogenetic trees is newick, but some alternatives exist, including the one used by MAFFT for very large trees, in which the sequence names are replaced by their original index (1..N) in the sequence file. Seq_reformat can be used to cross convert these two formats:

::

  $$: t_coffee -other_pg seq_reformat -in sample_3Dseq1.tree -in2 sample_3Dseq1.fasta -action +newick2mafftnewick

And the reverse operation when processing a mafft guide tree

::
  $$: t_coffee -other_pg seq_reformat -in sample_3Dseq1.mafftnewick -in2 sample_3Dseq1.fasta -action +mafftnewick2newick



Manipulating structure files (PDB)
==================================
Extracting a structure
----------------------
There are many reasons why you may need a structure. T-Coffee contains a powerful utility named **extract_from_pdb** that makes it possible to fetch the PDB coordinates of a structure or its FASTA sequence without requiring a local installation. By default, the option **extract_from_pdb will** start looking for the structure in the current directory; it will then look it up locally (PDB_DIR) and eventually try to fetch it from the web (via a wget to the `PDB <http://www.rcsb.org>`_). All these settings can be customized using environment variables (see next paragraph). For instance if you want to fetch the chain E of the PDB structure 1PPG and/or its sequence in FASTA format, you can use:

::

  Fetch the structure:
  $$: t_coffee -other_pg extract_from_pdb -infile 1PPGE

  Fetch the correpsonding sequence:
  $$: t_coffee -other_pg extract_from_pdb -infile 1PPGE -fasta

Adapting extract_from_pdb to your own environment
-------------------------------------------------
If you have the PDB installed locally, simply set the variable PDB_DIR to the absolute location of the directory in which the PDB is installed. The PDB can either be installed in its divided form or in its full form. If the file you are looking for is neither in the current directory nor in the local PDB version, extract_from_pdb will try to fetch it from rcsb. If you do not want this to happen, you should either set the environment variable NO_REMOTE_PDB_DIR to 1 or use the **-no_remote_pdb_dir** flag:

::

  Setting up the environment:
  ##: export NO_REMOTE_PDB_FILE=1
  
  Running PDB extract:
  $$: t_coffee -other_pg extract_from_pdb -infile 1PPGE -fasta -no_remote_pdb_file

By default, T-Coffee also requires two important PDB files declared using the two following variables. These variables do not need to be set if the considered files are in the cache directory (default behavior):

::

  ##: export PDB_ENTRY_TYPE_FILE=<location of the file pdb_entry_type.txt>
  (Found at: ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_entry_type.txt)
  and
  ##: export PDB_UNREALEASED_FILE=<location of the file unrealeased.xml>
  (Found at: http://www.rcsb.org/pdb/rest/getUnreleased)


.. warning:: Since the file ``unreleased.xml`` is not part of the PDB distribution, T-Coffee will make an attempt to obtain it even when using the NO_REMOTE_PDB_DIR=1 mode. You must therefore make sure that the file PDB_UNREALEASED_FILE is pointing to is read and write.


***********************
Aligning Your Sequences
***********************
General comments on alignments and aligners
===========================================
What is a good alignment?
-------------------------
This is a tricky question, a good answer would be  **"a good alignment is an alignment that makes it possible to do good biology"**. In practice, the alignment community has become used to measuring the accuracy of alignment methods using structures. Structures are relatively easy to align correctly, even when the sequences have diverged quite a lot. The most common usage is therefore to compare structure based alignments with their sequence based counterpart and to evaluate the accuracy of the method using these criterions. Unfortunately it is not easy to establish structure based standards of truth. Several of these exist and they do not necessarily agree. To summarize, the situation is as roughly as follows:

  - **Above 40% identity**, all the reference collections do agree with one another and all the established methods give roughly the same results. These alignments can be trusted blindly.

  - **Below 40% identity**, all the reference collections stop agreeing and the methods do not give consistent results. In this area of similarity it is not necessarily easy to determine who is right and who is wrong, although most studies indicate that consistency based methods (T-Coffee, ProbCons, MAFFT-slow or MSAProbs) have an edge over traditional methods.

When dealing with distantly related sequences, the only way to produce reliable alignments is to use structural information. T-Coffee provides many facilities to do so in a seamless fashion. Several important factors need to be taken into account when selecting an alignment method:

  - **The best methods are not always the best**. Given a difficult dataset, the best method is only more likely to deliver the best alignment, but there is no guaranty it will do so. It is very much like betting on the horse with the best odds.

  - **The difference in accuracy between all the available methods is not incredibly high** (as measured on reference datasets). It is unclear whether this is an artifact caused by the use of 'easy' reference alignments, or whether this is a reality. The only thing that can change dramatically the accuracy of the alignment is the use of structural information.

  - **Keep in mind that these methods have only been evaluated by comparison with reference alignments (benchmarks)**. This is merely one criterion among many. In theory, these methods should be evaluated for their ability to produce alignments that lead to accurate trees, good profiles or good models. Unfortunately, these evaluation procedures do not yet exist.


The main methods and their scope
--------------------------------
.. note:: There are many MSA packages around, the most common ones being ClustalW, MUSCLE, MAFFT, T-Coffee and ProbCons; amongst the latest ones, you can find phylogeny aware aligners (PRANK and SATé) and modifed/improved consistency based aligners (MSAProbs). You can almost forget about the other packages, as there is virtually nothing you could do with them that you will not be able to do with these packages. All these packages offer a complex trade-off between speed, accuracy and versatility.

ClustalW is really everywhere...
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
ClustalW is still the most widely used Multiple Sequence Alignment package. Yet things are changing fast and different tests have consistently shown that ClustalW is neither the most accurate nor the fastest package around. This being said, ClustalW is everywhere and if your sequences are similar enough, it should deliver a fairly reasonable alignment.

MAFFT/MUSCLE to align big datasets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you have many sequences to align MUSCLE or MAFFT are the obvious choice. MAFFT is often described as the fastest and the most efficient. This is not entirely true, in its fast mode (FFT-NS-1), MAFFT is similar to MUSCLE and although it is fairly accurate, about 5 points less accurate than the consistency based packages (ProbCons and T-Coffee). In its most accurate mode (L-INS-i) MAFFT uses local alignments and consistency, however, it becomes much more accurate but also slower, and more sensitive to the number of sequences. More recently, we have seen growing the number of **(ultra) large scale** aligners such as Clustal Omega, PASTA, UPP, and we hope soon the large scale version of T-Coffee (called MEGA-Coffee).

**Suitable for**:
 - Distance-based phylogenetic reconstruction (NJ trees)
 - Secondary structure prediction

**Not suitable for**:
 - Profile construction
 - Structure modeling
 - 3D prediction
 - Function analysis

T-Coffee/ProbCons, slow but accurate !!!
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
T-Coffee works by first assembling a library and then by turning this library into an alignment. The library is a list of potential pairs of residues. All of them are not compatible and the job of the algorithm is to make sure that as many possible constraints as possible find their way into the final alignment: it is very much like trying to choose a meeting date, and each one says something like 'I need my Monday morning', 'I can't come on Thursday afternoon', and so on. In the end you want a schedule that makes everybody happy, if possible. The nice thing about the library is that it can be used as a media to combine as many methods as one wishes. It is just a matter of generating the right constraints with the right method and compile them into the library. ProbCons and MAFFT (L-INS-i) uses a similar algorithm, but with a Bayesian twist in the case of ProbCons. In practice, however, ProbCons and T-Coffee give very similar results and have similar running time. MAFFT is significantly faster.

**Suited for**:
 - Profile reconstruction
 - Structure modeling
 - Function analysis
 - 3D prediction

Choosing the right package (without flipping a coin !)
------------------------------------------------------
Each available package has something to go for it, it is just a matter of knowing what you want to do !! T-Coffee is probably the most versatile, but it comes at a price, its default aligner being currently slower than many alternative packages. In the rest of this tutorial we give some hints on how to carry out each of these applications within the T-Coffee framework.

================= ====== ===== ======== ======== ======== 
Packages          MUSCLE MAFFT ProbCons T-Coffee ClustalW 
================= ====== ===== ======== ======== ======== 
Accuracy          ++     +++   +++      +++      \+        
<100 Seq.         ++     ++    +++      +++      \+        
>100 Seq.         +++    +++   \-       \+       \+        
Remote Homologues ++     +++   +++      +++      \+        
MSA vs Seq.       \-     \-    +++      +++      +++      
MSA vs MSA        \-     \-    \-       +++      +++      
>2 MSAs           \-     \-    \-       +++      \-        
Seq. vs Struc.    \-     \-    \-       +++      \+        
Splicing Var.     \-     +++   \-       +++      \-        
Reformat          \-     \-    \-       +++      ++       
Phylogeny         \-     \-    \-       \+       ++       
Evaluation        \-     \-    \+       \+++     \-        
Speed             +++    +++   \+       \+       ++       
================= ====== ===== ======== ======== ======== 

Table 1. Relative possibilities associated with the main packages. In any of the situations corresponding to each table line, (+++) indicates that the method is the best suited, (++) indicates that the method is not optimal but behaves reasonably well, (+) indicates that it is possible but not recommended (-) indicates that the option is not available.

===================== ====== ===== ======== ======== ======== 
Packages              MUSCLE MAFFT ProbCons T-Coffee ClustalW 
===================== ====== ===== ======== ======== ======== 
Dist Based Phylogeny  +++    +++   ++       ++       ++       
ML or MP Phylogeny    ++     +++   +++      +++      ++       
Profile Construction  ++     +++   +++      +++      ++       
3D Modeling           ++     ++    ++       +++      \+        
2D Predictions        +++    +++   ++       ++       ++       
===================== ====== ===== ======== ======== ======== 

Table 2. Most Suitable Appplications of each package. In any of the situations corresponding to each table line, (+++) indicates that the method is the best suited, (++) indicates that the method is not optimal but behaves reasonably well, (+) indicates that it is possible but not recommended (-) indicates that the option is not available.


Computing simple MSA with T-Coffee 
==================================
General considerations
----------------------
T-Coffee aligner is by default parallelized, meaning that it can use multiple cores when running on a cluster or a computer. By default, T-Coffee will use all available processors to run, but you can parallelize the differents steps and allocate the number of cores you want with the flag **-multi_core** or **n_core**. For more details, refer to the chapter **T-Coffee Technical Documentation**, subsection **CPU control**. 

::

  Using a single core:
  $$: t_coffee -in sample_seq1.fasta -multi_core=msa
  
  Using 12 cores:
  $$: t_coffee -in sample_seq1.fasta -n_core=12
  
Also when running T-Coffee, it displays a lot of information directly on screening while running from general information, options, results, warnings...if you want, you can reduce the display using **-no_warning** to remove all the warnings, or even more strict using **-quiet** removing any display while running T-Coffee.

Computing a simple MSA (default T-Coffee)
-----------------------------------------
T-Coffee default mode will simply compute a Multiple Sequence Alignment of the sequences you provided in input (command 1). It will display the final MSA on the screen and in several files according to the format you asked with command 2 (by default, the MSA is stored in a file .aln in ClustalW format). The headline of the alignment file contains important information such as the version of T-Coffee used, the CPU time, the overall consistency score (normalized to 100 or 1000 depending on the version of T-Coffee) and the total length of the MSA: it is quite practical to have a quick glance at the result. 

::

  Command 1: default MSA
  $$: t_coffee proteases_small.fasta

  Command 2: default MSA, multiple output files
  $$: t_coffee proteases_small.fasta -output=clustalw,fasta_aln,msf
  
Each time you run T-Coffee, 3 files are always generated:

 - ``proteases_small.aln``: the alignment in ClustalW format
 - ``proteases_small.dnd``: the guide tree in Newick format
 - ``proteases_small.html``: the colored MSA in html format

.. warning:: The guide tree is not a phylogenetic tree, it is used in the alignment process for clustering the sequences. 

.. tip:: You can visualize the colored html file with any browser/software you prefer. The display of the sequences should be aligned and formatted; if not, use another browser, it works quite well with Firefox, Safari, etc... If you need to do more sophisticated modifications on your MSA, we recommend to use `Jalview <http://www.jalview.org/>`_ which incorporate the T-Coffee color scheme.

Aligning multiple datasets/Combining multiple MSAs
--------------------------------------------------
If your sequences are spread across several datasets, you can give all the files you want (the limit is 200) via the flag **-seq**, and in any format you want. Just know that 1) if you give an alignment, the gaps will be reset and your alignment will only provide sequences, 2) sequences with the same name between two files are assumed to be the same sequence, 3) ff their sequences differ, they will be aligned and replaced by the consensus of that alignment (process known as sequence reconciliation). To align multiple datasets:

::

  $$: t_coffee -seq=proteases1_small.fasta,proteases2_small.aln -output=clustalw,fasta_aln,msf


You may also have a bunch of alignments (with the same sequences) that you have either precomputed, assembled manually or received from a colleague. You can also combine these alignments. For instance, let us imagine we generated 4 alignments with ClustalW using different gap penalties. To combine them into ONE single alignment, use the **-aln** flag. The final score indicates a high level of consistency (91%) between all these MSAs, meaning that the final MSA is probably correct.

::

  Your 4 different MSAs:
  ##: clustalw -infile=proteases_small.fasta -gapopen=0 -outfile=g0.aln
  ##: clustalw -infile=proteases_small.fasta -gapopen=-5 -outfile=g5.aln
  ##: clustalw -infile=proteases_small.fasta -gapopen=-10 -outfile=g10.aln
  ##: clustalw -infile=proteases_small.fasta -gapopen=-15 -outfile=g15.aln

  Combining multiple MSAs:
  $$: t_coffee proteases_small.fasta -aln g0.aln g5.aln g10.aln g15.aln -output=clustalw,html

Estimating the diversity in your alignment
------------------------------------------
It is easy to measure the level of diversity within your MSA with the **-output** option of **seq_reformat**, it will output all the pairwise identities, as well as the average level of identity between each sequence and the others. There are two possibilities given that your input are unaligned sequences or not: **-output sim_idscore** realign your sequences pairwise so it can accept unaligned or aligned sequences alike; **-output sim** computes the identity using the sequences as they are in your input file so it is only suited for MSAs. You can after redirect, sort and grep the output in order to select the sequences you are interested in.

::

  $$: t_coffee -other_pg seq_reformat -in sample_seq1.aln -output sim


Comparing alternative alignments
--------------------------------
If you change the parameters, you will end up with alternative alignments. It can be interesting to compare them quantitatively. T-Coffee comes along with an alignment comparison module named **aln_compare**. You can use it to estimate the amount of difference between your two alignments either using the Sum-of-Pair score or the column score using the flag **-compare_mode** (sp or column). By default aln_compare returns the SoP score:

::

  $$: t_coffee -other_pg aln_compare -al1 b80.aln -al2 b30.aln -compare_mode sp


This comparison will return the following result:

::

  *****************************************************
  seq1       seq2          Sim   [ALL]           Tot  
  b80           19         33.6    94.2 [100.0]   [79686]


The interpretation of this output is as follow: b80 is the reference MSA, it contains 19 sequences with an average identity of 33.6%, and is 94.2% identical to the second MSA b30.aln (79686 pairs to be precise). Of course, this does not tell you where are the good bits, but you can get this information for instance residues that have lost more than 50% of their pairing partner between the two alignments are in lower case (command 1) or converted in any character you want (command 2).
The sp score is defined as the number of aligned residue pairs (i.e. not gaps) that are common between the reference (-al1) and the test (-al2) divided by the total numbers of pairs in the reference, while excluding gaps. 

Other modes of comparison are available including the TC score

::

  $$: t_coffee -other_pg aln_compare -al1 b80.aln -al2 b30.aln -compare_mode tc

The tc score is defined as the number of columns in the test MSA (-al2)  that are common between the reference and the test, divided by the total number of columns in the reference. This measure, originally reported in Balibase can be problematic as it treats equaly highly and sparcely populated columns. For this reason the 


::

  $$: t_coffee -other_pg aln_compare -al1 b80.aln -al2 b30.aln -compare_mode column

The column score is defined as the total number pairs occuring in columns entirely identical between the reference (-al1) and the test (-al2) divided by the total number of pairs in the reference (excluding gaps). This makes the column score the equivalent of a weighted TC score.

It is also possible to print out the conserved positions between MSAs using the following commands. 

:: 

  Command 1:
  $$: t_coffee -other_pg aln_compare -al1 b30.aln -al2 p350.aln -output_aln \
      -output_aln_threshold 50

  Command 2:
  $$: t_coffee -other_pg aln_compare -al1 b30.aln -al2 p350.aln -output_aln \
      -output_aln_threshold 50 -output_aln_modif x


.. note:: It is important to take into account that the comparison between -al1 and -al2 is not symetrical. It returns the fraction of -al1 residues pairs - by index - that are recovered in -al2. If -al1 is the reference alignment, this number is more or less the equivalent of a sensistivity measure (i.e. how much True positives are recovered). If -al1 is the tartget alignment, this measure becomes the equivalent of a specificity measure (i.e. what is the fraction of residue pairs in the target that are part of the reference). 

.. tip:: aln_compare is particularly interesting if you are modifying the default parameters of T-Coffee and want to monitor the effects of your modifications. 

Comparing subsets of alternative alignments
-------------------------------------------

It is also possible to restrict the comparison to a specific subset of residues within the target alignment (for instance the core residues if using a structural reference, or the catalytic residues if using functional annotation alignment). This requires a cache in which residues are annotated as being part of the core or not. The cache is a plain text ascii file that can be produced in any way convenient to you. We show below how to produce a cache based on the consistency score and how to use this cache for further comparisons.

::

	$$: t_coffee b11.fa -output aln score_ascii

This command generates two output files

 - ``b11.aln``: the alignment in ClustalW format
 - ``b11.score_ascii``: a local reliability score with every residue having a score between 0 and 9

The target file must then be reformatted so as to decide which residues will be part of the evaluation (c=> core) and which ones will be ignored. In this example we will set the threshild at 7 which means that the comparison will ignore all residues that do not have 70% or more support from the T-Coffee library:

::

	$$: t_coffee -other_pg seq_reformat -in b11.score_ascii -input number_aln -action +convert 0123456i 789c -output fasta_seq > b11.cache


One can then produce a target file with another method
::

	$$: t_coffee b11.fa -method mafft_pair -outfile b11.tcmafft

With the cache, it is then possible to compare the original alignment (b11.aln) while only considering pairs of core residues
::

	$$: t_coffee -other_pg aln_compare -al1 b11.aln -al2 b11.tcmafft -st b11.cache pep -io_cat '[c][c]=[core]'

This comparison will return the following result:

::
  *****************************************************
  seq1       seq2          Sim   [core]          Tot  
  b11           11         32.8    99.5 [ 63.6]   [26044]


The output can be interpreted as follows: the core residues - defined as c-c pairs constitute 63.6% of the -al1 alignment. 93.5% of these pairs are recovered in the -al2 alignment. Note that the -io_cat makes it possible to declare as many categories and combinations of categories as one wishes to have. For instance:

::

	$$: t_coffee -other_pg aln_compare -al1 b11.aln -al2 b11.tcmafft -st b11.cache pep\
             -io_cat '[c][c]=[core];[c][i]=[mixed];[ci][ci]=[all1];[*][*]=[all2]'
	


Modifying the default parameters of T-Coffee
--------------------------------------------
.. note:: The main parameters of T-Coffee are similar to those of ClustalW, including a substitution matrix and some gap penalties. In general, T-Coffee's default is adequate. If, however, you are not satisfied with the default parameters, we encourage you to change the following parameters. Interestingly, most of what we say here holds reasonably well for ClustalW.

Can you guess the optimal parameters?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Here is another tricky question...and the general answer is NO. The matrix and the gap penalties are simplistic attempts at modeling evolution. While the matrices do a reasonable job, the penalties are simply inappropriate: they should have a value that depends on the structure of the protein and a uniform value cannot be good enough. Yet, since we do not have better we must use them...In practice, this means that parameter optimality is a very *ad hoc* business. It will change from one dataset to the next and there is no simple way to predict which matrix and which penalty will do better. The problem is also that even after your alignment has been computed, it is not always easy to tell whether your new parameters have improved or degraded your MSA. 

There is no systematic way to evaluate an MSA. In general, people visually evaluate the alignment, count the number of identical columns and consider that one more conserved column is good news. If you are lucky you may know a few functional features that you expect to see aligned. If you are very lucky, you will have one structure and you can check the gaps fall in the loops. If you are extremely lucky, you will have two structures and you can assess the quality of your MSA. An advantage of T-Coffee is the fact that the overall score of the alignment (i.e. the consistency with the library) is correlated with the overall accuracy. In other words, if you alignment score increases, its accuracy probably increases also. All this being said, consistency is merely an empirical way of estimating the change of parameters and it does not have the predictive power of a BLAST E-Value.

Changing the substitution matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
T-Coffee only uses the substitution matrix to make the pairwise alignments that go into the library. These are all the global alignments of every possible pair of sequences, and the ten best local alignments associated with every pair of sequences. 

 - By default, these alignments are computed using a Blosum62 matrix, but you can use any matrix you fancy instead, including: pam120mt, pam160mt, pam250mt, pam350mt, blosum30mt, blosum40mt, blosum45mt, blosum50mt, blosum55mt, blosum62mt, blosum80mt, or even user-provided matrices in the BLAST format (see **T-Coffee Technical Documentation**).

 - PAM matrices: These matrices are allegedly less accurate than the Blosum. The index is correlated to the evolutionary distances, you should therefore use the pam350mt to align very distantly related sequences.

 - Blosum matrices: These matrices are allegedly the most accurate. The index is correlated to the maximum percent identity within the sequences used to estimate the matrix. you should therefore use the blosum30mt to align very distantly related sequences. Blosum matrices are biased toward protein core regions, explaining why these matrices tend to give better alignments, since by design, they can capture the most evolutionary resilient signal contained in proteins.

Unless you have some structural information available, the only way to tell whether your alignment has improved or not is to look at the score. For instance, if you compute the two following alignments:

::

  $$: t_coffee proteases_small.fasta -matrix=blosum30mt -outfile=b30.aln
  $$: t_coffee proteases_small.fasta -matrix=blosum80mt -outfile=b80.aln

You will get two alignments that have roughly the same score but are slightly different. You can still use these two alternative alignments by comparing them to identify regions that have been aligned identically by the two matrices. These regions are usually more trustworthy.

Changing gap penalties
^^^^^^^^^^^^^^^^^^^^^^
.. important:: Gap penalties are the core of the matter when it comes to MSAs. An interesting feature of T-Coffee is that it does not really need such penalties when assembling the MSA, because in theory the penalties have already been applied when computing the library. This is the theory, as in practice penalties can help improve the quality of the alignment.

The penalties can be changed via the flags **-gapopen** for the gap opening penalty and via **-gapext** for the gap extension penalty. The range for gapopen are [-500,-5000], the range for the extension should rather be [-1,-10]. These values do not refer to a substitution matrix, but rather to the values range of the consistency estimation (i.e. ratio) normalized to 10000 for a maximum consistency. The default values are **-gapopen=-50, -gapext=0**. The reasons for these very low values are that they are meant to be cosmetic only, since a trademark of T-Coffee (inherited from Dialign) is not to need explicit penalties. Yet, we know for a fact that alignments with higher gap penalties often look nicer (for publications) and are sometimes more accurate. For instance, you can try:

::

  $$: t_coffee proteases_small.fasta -gapopen -100 -gapext -5

This gap penalty is only applied at the alignment level (i.e. after the library was computed). If you want to change the gap penalties of the methods used to build the library, you will need to go deeper...Two methods are used by default to build the library (command 1). One does global pairwise alignments and is named slow_pair, the other is named lalign_id_pair and produces local alignments. These methods are specified via the **-method** flag. Usually you do not need to write it because it is the default, but if you want to change the default parameters of the constituting methods (command 2), you will need to do so explicitly (the default parameters are for lalign_id_pair **GOP=-10, GEP=-4, MATRIX=blosum50mt** and for slow_pair **GOP=-10, GEP=-1 and MATRIX=blosum62mt**. Using the command 2, the library is now computed using the Blosum62mt with lalign, rather than the Blosum50mt; the good news is that when using this matrix, the score of our alignment increases from 48 to 50. We assume this new alignment is therefore more accurate than the previous one.

::

  Command 1: default T-Coffee
  $$: t_coffee proteases_small.fasta -method=lalign_id_pair,slow_pair

  Command 2: modifiying the parameters
  $$: t_coffee proteases_small.fasta -method lalign_id_pair@EP@MATRIX@blosum62mt, \
      slow_pair -outfile proteases_small.b62_aln

.. warning:: It only makes sense to compare the consistency score of alternative alignments when these alignments have been computed using the same methods (lalign_id_pair and slow_pair for instance).


Aligning (very) large datasets
==============================
Aligning (very) large datasets with MUSCLE
------------------------------------------
To run MUSCLE you can try one of the following command; don't hesitate to MUSCLE tutorial or help to get more information.

::

  Default mode:
  ##: muscle -in proteases_large.fasta > proteases_large.muscle
  
  Fast mode (less accurate):
  ##: muscle -in proteases_large.fasta -maxiters 1 -diags -sv -distance1 kbit20_3 \
  > proteases_large.muscle

Aligning (very) large datasets with MAFFT
-----------------------------------------
MAFFT is can align large datasets by default however it is better to use the fastest mode with MAFFT using the **retree** parameter; don't hesitate to MAFFT tutorial or help to get more information.

::
  
  Default mode:
  ##: mafft input > output
  Fast mode:
  ##: mafft --retree 2 input > output

Aligning (very) large alignments with T-Coffee
----------------------------------------------
T-Coffee is not very well gifted for aligning large datasets (for now), but you can give it a try using a special option that generates approximate fast alignments (command 1). These MSAs should roughly have the same accuracy as ClustalW, and are quite acceptable for sequences more than 40% identical. This mode works by only considering the best diagonals between two sequences, and by default all the diagonals with substitution score >0 are considered. You can lower this value with the flag **-ndiag** to reduce the running time (command 2). This will be very useful if you have long and very similar sequences to align (DNA for instance).

::

  Command 1:
  $$: t_coffee proteases_large.fasta -mode quickaln

  Command 2:
   $$: t_coffee proteases_large.fasta -mode quickaln -ndiag=10

Another alternative to align large datasets is a special mode of T-Coffee, fm-Coffee (command 3), derived from M-Coffee (see next section) and designed to be fast and able to handle large datasets (it is used for example in Ensembl). To do so, T-Coffee used three different fast aligners: MAFFT, MUSCLE and Kalign. 

::

  Command 3:
  $$: t_coffee proteases_large.fasta -mode fmcoffee

.. tip:: Once you have your large MSA, you can always shrink/trim it using reformatting options (see previous section) for instance by extracting the most informative sequences or by defining a % identity cut-off.

.. note:: In the last 10 years, a special effort have been made to improve large scale alignment leading to the development of few new methods among which Clustal Omega, PASTA, UPP and we hope soon a MEGA-Coffee aligner. These methods are not incorporated in T-Coffee so if your datasets are really large (>5000 sequences) don't hesitate to use these methods instead.


Using many methods at once
==========================
One of the most common situation when building MSAs is to have several alignments produced by different alternative methods, and not knowing which one to choose. In this section, we show you how to use M-Coffee to combine many alignments into one single alignment, or how you can specify only the methods you want. M-Coffee is not always the best method, but extensive benchmarks on BAliBASE, PREFAB and HOMSTRAD have shown that it delivers the best alignment 2 times out of 3. If you do not want to use the methods provided by M-Coffee, you can also combine precomputed alignments. 

Using third party aligner via T-Coffee
--------------------------------------
T-Coffee is installed along with many aligners necessary to run M-Coffee for instance, and many more. If you type **t_coffee**, it will display on the screen the different t_coffee options and all the methods included. If you look carefully, you will see that most of the methods exist under two denominations: 1) **<method>_msa** or 2) **<method>_pair**. In the first case, it means that T-Coffee will use the specified method to run your MSA, so you can easily have a ClustalW or a MAFFT alignment using T-Coffee. In the second case, you ask T-Coffee to align every pair of sequence with the specified methods, the final MSA will be computed using the T-Coffee consistency between all the pairs. Go to the **Integrating External Methods in T-Coffee** if you want more information.

Using all the methods at the same time: M-Coffee
------------------------------------------------
To use M-Coffee (M stands for Meta aligner), you will need several packages to be installed (see **T-Coffee Installation** and section **Integrating External Methods in T-Coffee**). If you did a default installation, all the software you need should be there. M-Coffee is a special mode of T-Coffee that you can call using the flag **-mode mcoffee**. It will align your sequence using 8 different aligners: ClustalW, POA, MUSCLE, ProbCons, MAFFT, Dialing-T, PCMA and T-Coffee:

::

  $$: t_coffee proteases_small.fasta -mode mcoffee -output clustalw, html

The final MSA is a combination of all methods. The alignment is colored with the T-Coffee consistency color scheme, but in this case the colors will reflect the consistency between methods: 1) regions in red have a high consistency, so all the methods agree and you can expect them to be fairly accurate, 2) regions in green/blue have the lowest consistency, meaning that all the methods deliver different alignment in these regions and you should not trust them. Overall this alignment has a score of 951 (1000 being the max), which means that it is roughly 95% consistent with the entire collection; this is a fairly high index meaning that you can trust your alignment. 

Using selected methods to compute your MSA
-------------------------------------------
Using the 8 methods predefined in M-Coffee can sometimes be a bit heavy, if you only want to use a subset of your favorite methods, you should know that each of these methods is available via the **-method** flag. You can make all the combination you want !!! For instance, to combine MAFFT, MUSCLE, T-Coffee and ProbCons, you can use:

::

  $$: t_coffee proteases_small.fasta -method=t_coffee_msa,mafft_msa,probcons_msa, \
      muscle_msa -output=html


Aligning profiles 
=================
Sometimes, it is better to prealign a subset of your sequences, and then to use this small alignment as a master for adding sequences (sequence to profile alignment) or even to align several profiles together if your protein family contains distantly related groups. T-Coffee contains most of the facilities available in ClustalW to deal with profiles, and the strategy we outline here can be used to deal with large datasets.

Aligning sequence(s) to profile(s)
----------------------------------
Assuming you have multiple alignment(s) (sproteases_small.aln) or profile(s) here is a simple strategy to align sequence(s) to your profile(s). It can align a variable number of sequences from 1 to N, with a variable number of profiles from 1 ot N: you can mix sequences and profiles in any proportion you like. 

::

  Adding one sequence to your MSA:
  $$: t_coffee proteases_oneseq.fasta -profile proteases_small.aln

  Adding many sequences to many profiles:
  $$: t_coffee -profile=prf1.aln,prf2.aln,prf3.aln -outfile=combined_profiles.aln

.. warning:: You can also use all the methods you want but be aware when using external methods that profiles are nto always supported. When it is not, it is replaced with its consensus sequence which will not be quite as accurate. Methods supporting full profile information are: lalign_id_pair, slow_pair, proba_pair, clustalw_pair and clustalw_msa. All the other methods (internal or external) treat the profile as a consensus (less accurate).

Computing very accurate (but slow) alignments with PSI/TM-Coffee
-----------------------------------------------------------------
PSI-Coffee is currently the most accurate mode of T-Coffee but also the slowest. Its principle is rather simple: it associates every sequence with a profile of homologous sequences gathered using BLAST on a sequence database (nr by default). PSI-Coffee then uses the profiles instead of the initial sequences to makes a multiple profile alignment (command 1). In a last step, your profiles are replaced by their initial query sequence from your initial dataset and returns a MSA of your sequences. PSI-Coffee can also use reduced database instead of nr (installed locally) in order to speed-up the process. A special mode, TM-Coffee, exists using PSI-Coffee but specialized to align transmembrane proteins using a reduced database of TM proteins and also including a prediction of transmembrane domains with the flag **-template_file PSITM** (command 2). It is much faster as the search database is limited to known transmembrane protein, however, it applies in only specific cases unlike PSI-Coffee which is a general method. You can find more information about TM-Coffee `here <http://tcoffee.crg.cat/apps/tcoffee/tutorial_tmcoffee.html>`_. If you want to specify a local BLAST version and a local database of your choice, just add to your command line the flags **-blast_server** and **-protein_db** and the corresponding paths.

::

  Command 1: PSI-Coffee
  $$: t_coffee proteases_small.fasta -mode psicoffee
  
  Command 2: TM-Coffee
  ##: t_coffee proteases_small.fasta -mode psicoffee -template_file PSITM -protein_db <database>
  

.. warning:: PSI/TM-Coffee requires BLAST and a database to search; if you don't have BLAST installed locally, it will use the BLAST default of T-Coffee. More importantly, if you don't specify a reduced database for TM-Coffee, it will run on nr and be equal to PSI-Coffee.

.. hint:: When running PSI/TM-Coffee, T-Coffee will use the BLAST EBI by default; it can happen that the web service is unavailable from time to time, T-Coffee will return a warning asking you to use NCBI BLAST instead. You can either do that or wait and rerun your job later on.


Using protein 2D/3D structural information 
==========================================
Using structural information when aligning sequences is very useful. The reason is that structures diverge slower than sequences. As a consequence, one may still find a discernable homology between two sequences that have been diverging for a long time beyond recognition using their corresponding structure. Yet, when assembling a structure based MSA, you will realize that these sequences contain key conserved residues that a simple alignment procedure was unable to reveal. We show you in this section how to make the best of T-Coffee tools to incorporate structural information in your alignment.

Using 3D structures: Expresso/3D-Coffee
---------------------------------------
Requirements to run Expresso/3D-Coffee
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Expresso needs BLAST you provide with the tag **-blast** and it also requires the PDB database you can specify via the flag **-pdb_db**. The good news is that it is not mandatory to have BLAST or the PDB installed locally as T-Coffee can automatically fetch the structures directly from RCSB (the home of PDB) using EBI BLAST web service. 

How does Expresso/3D-Coffee work?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Expresso/3D-Coffee is probably one of the most important improvement of T-Coffee. The principle is simple: 1) it first runs a BLAST for every sequence in your dataset against the PDB and finds (or not) a structure similar in sequence (35% identity by default) to be used as a template for structurally aligning your sequence, 2) the correspondence between the query sequences and the templates are stored in a template file which is automatically generated by Expresso. Here lies the difference between 3D-Coffee and Expresso: when running Expresso, fetching structures and creating the template file are automated (so you can reuse it for other applications) while using 3D-Coffee is a bit more tricky as it requires the name of the sequences to correspond to the structure file name (and it does not fetch or create anything for you). Of course, you can create and use your own template with the tag **-template_file** with the same format presented here. At the end, whenever there are enough templates (minimum of two obviously), it will align sequences using the structural information, otherwise sequences will be aligned using the standard T-Coffee aligner.  Of course, if your dataset only contains structures, your alignment becomes a structural alignment. 

::

  >sp|P08246|ELNE_HUMAN  _P_ 1PPGE
  >sp|P20160|CAP7_HUMAN  _P_ 1AE5
  ...

.. tip:: The PDB files used as a template should be in the current directory, otherwise you have to declare in the template file the whole path to find your templates.

Running Expresso/3D-Coffee
^^^^^^^^^^^^^^^^^^^^^^^^^^
Both Expresso (command 1) and 3D-Coffee (command 2) are modes of T-Coffee you can call with the flag **-mode**; the correspond to a preestablished configuration with default parameters. The template selection is indicated with the flag **-template_file** followed by either **PDB** (means BLAST will run remotely on the PDB) or **_SELF_P_** (means that the PDB identifier is the name of the sequence so there is no need to run BLAST) or a template file of your choice. As you can see in commands 1 and2, SAP (sap_pair) is the default structural aligner but you can choose alternative aligner(s) installed. You can give any combination of methods with the flag **-method**, but at least one has to be a structural aligner (command 3). You can also specify local installations of BLAST and PDB (command 4) **which is highly recommended if you want reproducible results**.

::

  Command 1: two ways of running Expresso
  $$: t_coffee three_pdb_two_seq.fasta -mode expresso
  $$: t_coffee three_pdb_two_seq.fasta -method sap_pair,slow_pair -template_file PDB

  Command 2: two ways of running 3D-Coffee
  $$: t_coffee three_pdb.fasta -mode 3dcoffee
  $$: t_coffee three_pdb.fasta -method sap_pair,slow_pair -template_file _SELF_P_

  Command 3: Choosing your own templates and methods
  $$: t_coffee three_pdb_two_seq.fasta -method mustang_pair,slow_pair -template_file \
      three_pdb_two_seq_pdb1.template_list
      
  Command 4: Running Expresso using a local BLAST/PDB 
  ##: t_coffee three_pdb_two_seq.fasta -mode expresso -blast=LOCAL -pdb_db=<PDB> \
      -pdb_type d -pdb_min_sim 95 -pdb_min cov 90 -cache $PWD 

.. tip:: By default, structures and structural pairwise alignments are stored in your local ~/.t_coffee/cache/ allowing Expresso to run faster if you reuse similar structures; you can choose to have all these files directly in your working directory by using **-cache=$PWD**. Don't forget to empty your cache directory from time to time otherwise your folder is just getting bigger and bigger (similar comment can be done for any template based mode of T-Coffee). 

Template search paramaters 
^^^^^^^^^^^^^^^^^^^^^^^^^^
To use Expresso, you have different option from an entirely automated procedure to tailored procedure, by selecting either your own structures or by defining different criteria for the template selection. You can have an exhaustive list in the **T-Coffee Technical Documentation** (subsection **Template based T-Coffee modes**) yet the most important parameters for the template selection are the following:

 - **-pdb_type**    : type of structure ("d" for diffraction/XRAY or "n" NMR structures)
 - **-pdb_min_cov** : minimum coverage between query sequence and template (from 0-100%)
 - **-pdb_min_sim** : minimum identity between query sequence and template (from 0-100%)


Aligning sequences and structures
---------------------------------
Mixing sequence profile and structure templates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you want to go further, and be even slower, you can use the accurate mode that will combine profile and structural information. If no structure is available, the template will be a profile (similar to PSI-Coffee, see subsection **Aligning Profiles**). It is probably one of the most accurate way of aligning sequences currently available as it tries to get as much information as possible.

::

  $$: t_coffee proteases_small.fasta -mode accurate

Aligning profile using structural information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you have two profiles to align, an ideal situation is when your profiles each contain one or more structures. These structures will guide the alignment of the profiles, even if they contain very distally related sequences. All you need is a template file that declares which sequences have a known structure. If you only want to align sequences, you can try:

::

  $$: t_coffee -profile=profile1_pdb1.aln, profile2_pdb2.aln -method sap_pair \
      -profile_template_file two_profiles.template_file
     

Using secondary structure predictions
-------------------------------------
T-Coffee can be used to predict secondary structures and transmembrane domains. For secondary structure predictions, the current implementation is only able to run GOR on either single sequences or on a bunch of homologues found by BLAST.

Single sequence prediction
^^^^^^^^^^^^^^^^^^^^^^^^^^
To make a secondary structure prediction T-Coffee used the GOR software. In the command line **-template_file=SSP** is a hard coded mode which prompts the computation of predicted secondary structures. Used this way, the method will produce for each sequence a secondary prediction file (``<sequence_name>.ssp``). GOR run on single sequences with a relatively low accuracy, which can be increased by coupling it with BLAST (command 2). When doing so, the predictions for each sequence are obtained by averaging the GOR predictions on every homologue as reported by a BLAST against nr. By default the BLAST is done remotely at the NCBI using the BLASTPGP web service of the EBI. Transmembrane structures can also be carried out simply, or following the same previous strategy (command 3 and 4).

::

  Command 1: 
  $#: t_coffee sample_seq1.fasta -template_file SSP
  
  Command 2: 
  $#: t_coffee sample_seq1.fasta -template_file PSISSP

  Command 3:
  $$: t_coffee sample_seq1.fasta -template_file TM

  Command 4:
  $$: t_coffee sample_seq1.fasta -template_file PSITM

Incorporation of the prediction in the alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
It is possible to use the secondary prediction (command 1) or the transmembrane domains prediction (command 2) in order to reward the alignment of similar elements:

::

  Command 1:
  $#: t_coffee sample_seq1.fasta -template_file PSISSP -method_evaluate_mode ssp -method \
      lalign_id_pair

  Command 2:
  $$: t_coffee sample_seq1.fasta -template_file PSITM -method_evaluate_mode tm -method \
      lalign_id_pair

The overall effect is very crude and can go up to overweighting by 30% the score obtained when matching two residues in a similar secondary structure state. The net consequence is that residues in similar predicted states tend to be aligned more easily.

.. hint:: In that case, the flag **-method** can only accept one single method

Using other secondary structure predictions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you have your own predictions, you can use them to run T-Coffee providing you give your own template file, where the file containing the secondary structure prediction is declared along with the sequence (with _E_).

::


  Command:
  $#: t_coffee sample_seq1.fasta -template_file sample_seq1_ssp.template -method_evaluate_mode \
      ssp -method lalign_id_pair
      
  Format of the template file:     
  >hmgl_wheat _E_ hmgl_wheat.ssp
  >hmgb_chite _E_ hmgb_chite.ssp
  >hmgl_trybr3 _E_ hmgl_trybr3.ssp
  ...

  Format of the prediction file:
  >hmgl_wheat
  CCCCCCCCCCCCHHHHHHHCCCCCCCCCHHHHHHHHHHHHHHHCCCCHHHHHHHHHHHHHHHCE

Output of the prediction
^^^^^^^^^^^^^^^^^^^^^^^^
You can output a color coded version of your alignment using the secondary predicted structure or transmembrane regions predictions:

::


  Secondary structure prediction:
  $#: t_coffee sample_seq1.fasta -template_file PSISSP -output sec_html

  Transmembrane regions prediction: 
  $$: t_coffee sample_seq1.fasta -template_file PSITM -output tm_html



Aligning RNA sequences 
======================
RNA sequences are very important and almost every-where these days. The main property of RNA sequences is to have a secondary structure that can be used to guide the alignment. While the default T-Coffee has no special RNA alignment method incorporated in, we have developped specific modes and tools for RNA alignment and analysis (see subsection **Manipulating RNA Sequences** for more details).

R-Coffee, many possibilities
----------------------------
Introduction
^^^^^^^^^^^^
R-Coffee is the special mode of T-Coffee developped to handle specifically RNA sequences. It has been proven far more accurate than T-Coffee default because of its specific design. It can be run as a standalone aligner (using secondary structure prediction) or using third party software.

R-Coffee: aligning RNA sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
R-Coffee uses predicted secondary structures via the software RNApfold, in order to improve RNA alignments. Running R-Coffee by default is rather simple (command 1) but as for T-Coffee, you can also specify the methods you prefer (command 2):

::

  Command 1: R-Coffee default
  $$: t_coffee sample_rnaseq1.fasta -mode rcoffee
 
  Command 2: R-Coffee selected methods
  $$: t_coffee sample_rnaseq1.fasta -mode rcoffee -method lalign_id_pair,slow_pair

Improving R-Coffee
^^^^^^^^^^^^^^^^^^
There are two modes we proposed to improve R-Coffee alignments: 1) using the best method for RNA sequences (namely Consan), 2) combining many methods to achieve a better reliability (RM-Coffee).

::

  Using Consan (best):
  $#: t_coffee sample_rnaseq1.fasta -mode rcoffee_consan

 
  Using multiple methods:
  $$: t_coffee sample_rnaseq1.fasta -mode rmcoffee

.. tip:: In order to know if a RNA alignment is better than another one, the best is to visualize the compensatory mutations of the secondary structure: look at the subsection **Manipulating RNA Sequences**.

Using SARA-Coffee
-----------------
SARA-Coffee is a structure based multiple RNA aligner. This is a new algorithm that joins the pairwise RNA structure alignments performed by SARA with the multiple sequence T-Coffee framework. Since setting up the SARA-Coffee dependencies (T-Coffee, `SARA <http://structure.biofold.org/sara/>`_, `X3DNA <http://x3dna.org/>`_, `Numpy <http://www.numpy.org/>`_, `Biopython <http://biopython.org/>`_, Perl, Python 2.7) can be tricky we provide a self-contained Vagrant VM, which downloads and configures all the required pieces of software for you. 

Installing SARA-Coffee VM
^^^^^^^^^^^^^^^^^^^^^^^^^
Follow the procedure:

::

  1) Install or update virtual box from https://www.virtualbox.org/wiki/Downloads

  2) Install or update vagrant from: http://www.vagrantup.com

  3) Clone the sara-coffee virtual machine (this project) in a convenient location
     ##: git clone https://github.com/cbcrg/sara-coffee-vm.git
     
  4) Enter the sara-coffee-vm folder and launch vagrant
     ##: cd sara-coffee-vm/
     ##: vagrant up  

The first time you run it, it will automatically download the virtual machine and all the packages required by SARA-Coffee. It may take some minutes to complete, so be patient. When it boots up and the configuration steps are terminated, login into the VM instance:

::

  1) Login in VM:
     ##: vagrant ssh
  
  2) Go to SARA-Coffee:
     ##: cd sara_coffee_package/
     
  3) Run SARA-Coffee:
     ##: ./sara_coffee.sh <input file> <output file>
     
The folder '/vagrant/' is shared between the Sara-Coffee virtual and your local machine. On your local machine, this folder is the one in which you started vagrant (i.e. sara-coffee-vm). When finished, stop the VM using the command **vagrant halt** or **vagrant destroy**, depending if you want to temporary stop the execution or delete permanently the VM with all its files.   

Docker image for SARA-Coffee
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
SARA-Coffee is also distributed as a Docker container. This will allow you to run it without having to install and configure each single dependency packages. If you have Docker installed simply pull the SARA-Coffee container by using the command 1 and run SARA-Coffee using the command 2:

::

  Command 1: Pull SARA-Coffee container
  ##: docker pull cbcrg/saracoffee**.  
  
  Command 2: Run SARA-Coffee
  ##: docker run -v $PWD:$PWD -w $PWD cbcrg/saracoffee <input> <output>**


.. Note:: this command assumes your input file is located in the working directory. If this is not the case, you will need to mount the input file path accordingly. 


Aligning DNA sequences
======================
Aligning DNA sequences
----------------------
MSA methods are not at their best when aligning DNA sequences. Whenever you can, try using a local MSA package like the Gibbs sampler; yet if you believe your DNA sequence are homologous over their entire length, you can use T-Coffee. In theory, the program automatically recognizes DNA sequences and uses appropriate methods, yet adding the **-type=dna flag** cannot do any harm...

::

  $$: t_coffee sample_dnaseq1.fasta -type=dna

The type declaration (or its automatic detection) triggers the use of the appropriate substitution matrix in most of the methods. However, if you would rather use your own matrix, use the following command where you should replace **idmat** with your own matrix in BLAST format.

::


  $$: t_coffee sample_dnaseq1.fasta -in Mcdna_fast_pair@EP@MATRIX@idmat


Pro-Coffee: Aligning functional DNA regions
-------------------------------------------
Pro-Coffee is a MSA method specifically designed for promoter regions or other orthologous DNA regions containing functional elements (enchancers for instance). Pro-Coffee takes nearest-neighbour nucleotide correlations into account when aligning DNA sequence. For this it first translates sequences into a dinucleotide alphabet and then does the alignment using a specifically designed dinucleotide substitution matrix. This matrix was constructed from binding sites alignments from the Transfac database. A benchmark on multispecies ChIP-seq data shows that validated binding sites will be better aligned than when using off-the-shelf methods. To run Pro-Coffee is easy by defaults (command 1), however the user can change empirically all gap parameters in order to improve the MSA (command 2). By default, parameters were optimized on benchmarks (GOP=-60/GEP=-1) but as every dataset is different you may want to see for yourself. Once your alignment is done, you can always use reformatting options in order to extract the regions of interest (command 3).

::

  Command 1: Pro-Coffee default
  $$: t_coffee -seq c18orf19.fa -mode procoffee

  Command 2: Pro-Coffee parameters
  $$: t_coffee -seq c18orf19.fa -method promo_pair@EP@GOP@-60@GEP@-1

  Command 3: Extracting regions
  $$: t_coffee -other_pg seq_reformat -in c18orf19.aln -action +extract_block \
      'ENSG00000177150' 1852 1983 > c18orf19_chipseq.aln

.. tip:: If you want more details, we suggest you follow the subsection **Quick Start, Tutorial (Practical Examples)** published in Nature Protocols (2011) or refer to the original article.

Splice variants
---------------
Splice variants are especially challenging for most MSA programs, because the splicing variants need very long gaps to be inserted, and most programs attempt to match as many symbols as possible. Standard programs like ClustalW or MUSCLE are not good at dealing with this situation and in our experience, the only programs that can do something with splice variants are those using local information like some flavors of MAFFT and T-Coffee. For instance, if you try muscle on the following dataset with the command 1, you will quickly realize your alignment is not very good and does not show where the alternative splicing co-occurs. On the other hand, if you use T-Coffee (command 2), things become much clearer. The reason why T-Coffee does better than other packages is mostly because it uses local information (lalign_id_pair) and is therefore less sensitive to long gaps. If the default mode does not work for your dataset, you can try to be a bit more aggressive and only use local information to compute your library (command 3). Of course, the most distantly related your sequences, the harder the alignment of splice variants.

::

  Command 1: using MUSCLE
  ##: muscle -in sv.fasta -clw -out sv_muscle.aln
  
  Command 2: using T-Coffee
  $$: t_coffee sv.fasta

  Command 3: using only local information
  $$: t_coffee sv.fasta -method lalign_id_pair

Noisy coding DNA sequences...
-----------------------------
When dealing with coding DNA, the right thing to do is to translate your DNA sequence and thread the DNA onto the protein alignment if you really need some DNA. However, sometimes, your cDNA may not be so clean that you can easily translate it (frameshifts and so on). Whenever this happens, try (no warranty) the following special method. The test case in three_dna_seq.fasta contains the DNA sequences of three proteases with a couple of frameshifts here and there. If you make a regular alignment of these sequences (command 1) you see that many gaps have sizes that are not multiple of 3 (codon size). When using an appropriate alignment method that takes into account all the frames at the same time, we get something much more meaningful (command 2). At the end, you can even recover the corrected protein sequence (command 3) using a special option **+clean cdna**, a small HMM that loops through each sequence and select the frame in order to maximize the similarity within the alignment (part of the **seq_reformat** utility).

::

  Command 1: Default alignment
  $$: t_coffee three_cdna.fasta

  Command 2: Special mode for cDNA
  $$: t_coffee three_cdna.fasta -method cdna_fast_pair

  Command 3: Recovering your protein sequences
  $$: t_coffee -other_pg seq_reformat -in three_cdna.aln -action +clean_cdna +translate



*************************
Evaluating Your Alignment
*************************
There are three possible strategies for evaluating your alignment of protein sequences:

 1) **Sequence based methods** like the CORE index and the TCS if you don't have any structure (quite often). These do pretty well in the core regions of an alignment (which can correspond to protein domains, fold, structural elements, etc...) but can be limited in more variable regions (which can correspond to loops, disordered regions, etc...).
 2) **Structure is a killer**, so if you have at least two structures available for your protein family, you are in an ideal situation and you can use the iRMSD. If you only have one structure available, we developped STRIKE to compare alternative alignment.
 3) **Another killer is the use of functional information**, but it is much less often at hand. If you know some residues MUST be aligned because they are functionally related. As the information is scarce and not standard, there is no automated procedure specifically designed for this kind of analysis but you can still setup an evaluation procedure with T-Coffee.

.. note:: Most of evalution methods are designed for protein sequences (notably structure based methods), however, T-Coffee via the TCS/CORE index offers some possibilities to evaluate also DNA alignments.

Sequence Based Methods
======================
The CORE index
--------------
.. note:: The CORE index is the basis of T-Coffee estimation of the consistency, however for evaluating alignment we recommend to use the TCS procedure describe in the next section. 

Computing the local CORE index
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The CORE index is an estimation of the consistency between your alignment and the computed library. The higher the consistency, the better the alignment. The score reported with every T-Coffee alignment is the concistency score (depending on the version it can be normalized to 100 or 1000). If you want to go further and estimate the local concistency (known as the CORE index), an html file is automatically created each time you run T-Coffee; it is colored version of your alignment where residues are colored according to their consistency score, from blue (low consistency) to red (high consistency). It is not full proof but in principle you can expect positions with a score above 6 to be correctly aligned.

Computing the CORE index of any alignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
You can evaluate any existing alignment with the CORE index. All you need to do is provide that alignment with the -infile flag and specify that you want to evaluate it. For more information on filtering/trimming an alignment using the CORE index score, refer to the subsection **Preparing Your Data: Reformatting, Trimming an more.../Modifying the data itself**.

::

  $$: t_coffee -infile=proteases_small_g10.aln -output=html -score

Transitive Consistency Score (TCS) 
----------------------------------
TCS is an alignment evaluation score that makes it possible to identify the most correct positions in an MSA. 
It has been shown that these positions are the most likely to be structuraly correct and also the most informative when estimating phylogenetic trees. The TCS evaluation and filtering procedure is implemented in the T-Coffee package and can be used to evaluate and filter any third party MSA (including T-Coffee MSA of course!).

.. warning:: TCS has been incorporated to T-Coffee recently, it means that not all distribution have the TCS implemented; you should install one the latest stable version of T-Coffee to have the TCS along with T-Coffee.

Evaluating an existing MSA 
^^^^^^^^^^^^^^^^^^^^^^^^^^
The TCS is most informative when used to identify low-scoring portions within an MSA. It is also worth noting that the TCS is not informative when aligning less than five sequences.

:: 

  $$: t_coffee -infile sample_seq1.aln -evaluate -output=score_ascii,aln,score_html

Output files: 

* ``sample_seq1.score_ascii``  displays the score of the MSA, the sequences and the residues. 
* ``sample_seq1.score_html`` displays a colored version score of the MSA, the sequences and the residues. 

.. warning:: The color code in the score_html indicates the agreement between the library and the considered alignment. It is important to understand that this score does not only depend on the input MSA, but it also depends on the library.

Filtering unreliable MSA positions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
TCS allows you to filter out from your alignment regions that appears unreliable according to the consistency score; the filtering can be made at the residue level or the column level:

:: 

  $$: t_coffee -infile sample_seq1.aln -evaluate -output tcs_residue_filter3, \
      tcs_column_filter3,tcs_residue_lower4

Output files: 

* ``sample_seq1.tcs_residue_filter3``  All residues with a TCS score lower than 3 are filtered out 
* ``sample_seq1.tcs_column_filter3``   All columns with a TCS score lower than 3 are filtered out 
* ``sample_seq1.tcs_residue_lower4``   All residues with a TCS score lower than 3 are in lower case

.. warning:: The TCS will create output files containing your results; if you rerun similar jobs or with the same name, TCS will not overwrite the previous outputs but append the new results in the already existing files.

All these output functions are also compatible with the default T-Coffee (command 1) when computing an alignment or with **seq_reformat** (command 2) using a T-Coffee ``<name>.score_ascii file``.

::

  Command 1:
  $$: t_coffee -seq sample_seq1.fasta -output tcs_residue_filter3, tcs_column_filter3, \
      tcs_residue_lower4

  Command 2:
  $$: t_coffee -other_pg seq_reformat -in sample_seq1.aln -struc_in sample_seq1.score_ascii \
      -struc_in_f number_aln -output tcs_residue_filter3

Weighting MSA for improved trees
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
One cool trick about TCS is the possibility to deliver weigthed MSA, where each columns is multiplied according to its consistency score; this appears to be particularly useful to build phylogenetic tree. Phylogenetic trees are evaluated using a bootstrap score so each column as the same weight, no matter its relevance, in the case of weighted TCS, the more reliable columns are more represented, therefore improving the support of informative and reliable positions of your MSA. 

:: 

  $$: t_coffee -infile sample_seq1.aln -evaluate -output tcs_weighted, tcs_replicate_100

Output files: 

* ``sample_seq1.tcs_weighted``       All columns are duplicated according to their TCS score 
* ``sample_seq1.tcs_replicate_100``  Contains 100 replicates in phylip format with each column drawn with a probability corresponding to its TCS score 

Note that all these output functions are also compatible with the default T-Coffee (command 1) when computing an alignment or with **seq_reformat** (command 2) using a T-Coffee .score_ascii file.

::

  Command 1:
  $$: t_coffee -seq sample_seq1.fasta -output tcs_weighted, tcs_replicate_100

  Command 2:
  $$: t_coffee -other_pg seq_reformat -in sample_seq1.aln -struc_in sample_seq1.score_ascii \
      -struc_in_f number_aln -output tcs_weighted

Using different libraries for TCS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
It is possible to change the way TCS reliability is estimated. This can be done by building different T-Coffee libraries. The proba_pair is the default aligner of T-Coffee that runs a pair-HMM to populate the library with residue pairs having the best posterior probabilities (command 1). You can also combine local and global alignners (command 2). There is also a fast alternative by using a special mode to run a series of fast multiple aligners; it is very fast and used by `ENSEMBL Compara <http://www.ensembl.org/info/genome/compara/index.html>`_ (command 3)

::

  Command 1:
  $$: t_coffee -infile sample_seq1.aln -evaluate -method proba_pair -output score_ascii, \
      aln, score_html

  Command 2:
  $$: t_coffee -infile sample_seq1.aln -evaluate -method proba_pair,lalign_id_pair \
      -output score_ascii,aln, score_html
      
  Command 3:
  $$: t_coffee -infile sample_seq1.aln -evaluate -method mafft_msa,kalign_msa,muscle_msa \
      -output score_ascii, aln, score_html
    
Working with coding DNA
^^^^^^^^^^^^^^^^^^^^^^^
When working with DNA, it is advisable to first align the sequences at the protein level and later thread back the DNA onto your aligned proteins. The filtering must be done in two steps, as shown below. Note that your DNA and protein sequences must have the same name. This first step produces the TCS evaluation file ``sample_prot_thread.score_ascii`` (command 1). Then, the **-out dna.replicates** option produces 100 DNA replicates with positions selected according to their aminoacid TCS score (command 2). Finally, the **-out dna.filtered** option will filter the DNA alignment according to their TCS column score.

::

  Command 1: evaluating the protein MSA
  $$: t_coffee -infile sample_prot_thread.aln -evaluate -output score_ascii

  Command 2: generating replicates for DNA 
  $$: t_coffee -other_pg seq_reformat -in sample_prot_thread.aln -in2 sample_dna_thread.fa \
      -struc_in sample_prot_thread.score_ascii -struc_in_f number_aln -output tcs_replicate_100 \
      -out dna.replicates
  
  Command 3: filtering out according to the TCS score
  $$: t_coffee -other_pg seq_reformat -in sample_prot_thread.aln -in2 sample_dna_thread.fa \
      -struc_in sample_prot_thread.score_ascii -struc_in_f number_aln -output tcs_column_filter5 \
      -out dna.filter  

Summary of the output options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
============================  ================
Flags        		               Description
============================  ================
**score_ascii**	              Outputs a TCS evaluation file
**score_html**	      	        Contains ascii format in html format
**score_pdf**	      	         Transfers score_html into pdf format
**sp_ascii**	      	          Reports the TCS score of every aligned pair in the target MSA
**tcs_residue_filter_N**      Removes all residues with a TCS score lower than `N`
**tcs_columns_filter_N**      Removes all columns with a TCS score lower than `N`
**tcs_weighted**	             Phylip format with duplicated columns according to their TCS score
**tcs_replicate_N**	          Generates `N` replicates of columns drawn according to their TCS score
============================  ================


Structural evaluation of MSAs
=============================
APDB/iRMSD
----------
What is the APDB/iRMSD?
^^^^^^^^^^^^^^^^^^^^^^^
APDB and the iRMSD are two closely related measures meant to evaluate the accuracy of a MSAt without using a structure based reference alignment. The iRMSD is a follow up of the APDB measure and we now recommend using the iRMSD rather than APDB. Although it may seem that the iRMSD was an attempt to get free iPODs from Apple, it is not (or at least we never got the iPODs). The iRMSD is a special RMSD (standing for intramolecular distances based RMSD) where the alignments are evaluated using the structural information of the sequences with known structures.

The strength of the iRMSD is its independence from a specific superposition models. When using the iRMSD to evaluate the score of a MSA, one does not need to superpose the two structures and deduce a sequence alignment that will then be compared with the target alignment. Given two aligned residues (X and Y) the iRMSD score is an attempt to estimate the neighborhood support for the XY alignment. This is done by measuring the difference of distances between X and Y and every other pair of aligned residues within the same sphere (W and Z). The iRMSD is obtained by measuring the average Root Mean Square (RMS) of these differences of distances. The first step of APDB/iRMSD is to measure the distances between the Ca (carbon alpha) of each residue and its neighbors. Neighborhood is defined as a sphere of radius ***-maximum_distance** (10 by default). However, by setting **-local_mode** to 'window', the sphere can be replaced with a window of 1/2 size **-maximum_distance** residues. 

The lower the iRMSD, the better the alignment. Yet, an alignment can obtain a good iRMSD score by simply having few aligned residues. To avoid this, the program also reports a normalized version of the iRMSD, the NiRMSD= MIN(L1,L2)*iRMSD/Number considered columns. It is recommended to use the NiRMSD to compare alternative alignments of different length. From a structural point of view, the NiRMSD has a meaning very similar to the iRMSD and it behaves in a similar fashion from a numerical point of view (similar ranges in Angstroms).

.. note:: APDB is an older measure less robust than the iRMSD; it is an attempt to estimate the fraction of pairs of residues whose alignment seems to be correct form a structural point of view. The higher APDB, the better the alignment, and conversely the lower the NiRMSD, the better the alignment.

How to efficiently use structural information?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
When it comes to evaluating MSAs, nothing is better than structural information. To use the methods we describe here, you will need to have at least two structures, similar enough (>60%) to sequences contained in your dataset. Here an outline of the best way to proceed:

 1) Try to include two structures with distantly related sequences, the other sequences being intermediates.
 2) Align your sequences without using the structural information (i.e. T-Coffee, MUSCLE, MAFFT...).
 3) Evaluate your alignment with iRMSD/NiRMSD (see later in this section); give this alignment the score S1.
 4) Realign your sequences but this time using structural information with Expresso.
 5) Measure the score of that alignment; the score will be S2.

If S1 and S2 are almost similar, it means your distantly related structures were well aligned and you can expect the intermediate sequences to be well aligned as well. If S2 is much better than S1, you can expect the structures to be well aligned in the second alignment, while there is no guarantee that the alignment of the intermediate sequences has improved as well, although in practice it often does.

Evaluating an alignment with the iRMSD package
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Let us evaluate the alignment produced by Expresso using the template file it returns. Evaluating the MSA with iRMSD will deliver a long output (all pairs are compared), the most interesting bit is at the bottom with the global iRMSD score in Angstrom (NiRMSD is the iRMSD score normalized to the length of the MSA). 

::

  Running the iRMSD:
  $$: t_coffee -other_pg irmsd proteases_small.aln -template_file proteases_small.template

  Result of the iRMSD evaluation:
  TOTAL     EVALUATED :  50.17 %
  TOTAL     APDB      :  83.44 %
  TOTAL     iRMSD     :  0.67 Angs
  TOTAL     NiRMSD    :  1.34 Angs

Evaluating alternative alignments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The strength of structure based alignments is that they make it possible to compare alternative alignments. In this case let us consider the following results in the table below (APDB in %, iRMSD/NiRMSD in Angstroms, and the evaluated columns in %). As expected, Expresso delivers the best alignment from a structural point of view. This makes sense, since Expresso explicitely USES structural information. The other figures show us that the structural based alignment is only marginally better than most sequences based alignments, yet with the notable exception of ClustalW.

======== ======== ======== ========= =========
 Method   APDB(%) iRMSD(A) NiRMSD(A) Eval. (%)
======== ======== ======== ========= =========
Expresso   83.44    0.67     1.34      50.17
T-Coffee   83.11    0.68     1.35      50.29
M-Coffee   83.08    0.68     1.36      50.00
ProbCons   83.10    0.68     1.35      50.28
MAFFT      82.99    0.68     1.35      50.25
Kalign     82.42    0.69     1.38      50.02
ClustalW   80.62    0.73     1.47      49.55
======== ======== ======== ========= =========

DEC - Distance Evolutionnary Conservation with msa2distances
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Using a similar approach it is possible to estimate the variation of intra-molecular distances across a multiple structural MSA for every pair of residue of every sequence. The following command also returns the average variation (stdev) for every residue and its neighbours withing <radius>

::

  $$: t_coffee -other_pg seq_reformat  -in proteases_small.aln -in2 proteases_small.template -action +msa2distances  20
  #$: t_coffee -other_pg seq_reformat  -in proteases_small.aln -in2 proteases_small.template -action +msa2distances  <radius: 20A> <threshold: 0.5 A> <min cluster: 4>	

The output comes in two sections that can be grepped with the first keyword:

::

  ##DECPAIR s1:               ##DECPAIR s1:   <seq_name> c1:  <column 1> c2: <column 2> r1:  <residue 1> r2:  <residue 2> pdbr1:  <pos filed on ATOM line> pdbr2:  pos filed on ATOM line> aa1: <amino acid 1> aa2: <amino acid 2> d:  <distance> avg_d:  <average distance in the column> stdev_d:  <standard deviation> Normalized_score: <standard_deviation/average_distance> N: <number of residue pairs> F: <number residue pairs/number sequences> ent1: <shannon entropy col 1> ent2: <shannon entropy col2>

The second fraction summarizes the fate of each residue. Note that the radius parameter potentially makes each residue in a column resulting having different levels of conservation

::

  ##DECRES s1:               <seq name> aa: <amino acid> c1: <msa column> r1:  <residue index 1..N> avg_stdev: <average stdev across radius AA> -Zscore:  <average zscores of the averaged stdev> normZ: <normalized zscore 0..100> Neighborhood:  <# of AA within radius> Radius:  <radius size in Angstrom as specified in the command line, 20 A by default>

This function also outputs colorized versions of the input PDBs using the PDB Bfactor filed. The output files are the following:

====================  ===================================================================
 File                 APDB(%) iRMSD(A) NiRMSD(A) Eval. (%)
====================  ===================================================================
*.bfactor2entropy.pdb replace thee Bfactor with the entropy measured on the MSA
*.bfactor2decres      replace the  Bfactor with the normalized average stdev of 
                      every residue within a sphere of radius <radius>
*.bfactor2deccluster  uses UPGMA to paint in the same shade all groups of residues larger 
                      than <min> with a pairwise stdev lower than <threshold> 
====================  ===================================================================



STRIKE: Contact based evaluations 
---------------------------------
STRIKE uses a contact matrix to evaluate an alignment using one or more structures. When doing so, the contacts are extracted from the sequences with known structures and are tnen projected onto the other sequences. Each sequence is then evaluated for the quality of its predicted contacts. When more than one structure is available, the results are provided using each structure as a template. The following Command line will carry this analysis. It takes as input an MSA containing at least one sequence with one of the templates defiend in the template file. The PDB files must be within the current directory.

::

  $$: t_coffee -other_pg seq_reformat -in sample_3Dseq1.aln -in2 sample_3Dseq1.template -action \
      +evaluate3D strike -output score_ascii -out sample_3D.score_ascii


The first part of the output is the STRIKE RAW SCORE. It reports the score of each sequence when using either one structure as a template, or averaging over all the structures; the highest value in every column is marked with a star. The various fields are as follows:

======== ================= ===========================================
Field    Name              Definition
======== ================= ===========================================
RS       Raw Score         Average STRIKE score
Rn	      Random Score      Average score when scrambling contacts
Bg       Background Score  Average score of all-against-all contacts
======== ================= ===========================================


It is also possible to use three extra parameters. The three parameters must be passed!!! Each one can be replaced with the string **"def"** if the default value is meant to be used. 

========== ======== =================================================
Parameter  Default  Definition
========== ======== =================================================
max        1.2      Max distance between two contacting AA (Angstrom)
end	       3	       Number of excluded neighbours between contacts
matrix     strike   File containing the STRIKE matrix 
========== ======== =================================================

:: 

  Generic common line for STRIKE:
  ##: t_coffee -other_pg seq_reformat -in sample_3Dseq1.aln -in2 sample_3Dseq1.template -action \
      +evaluate3D strike <max> <end> <matrix> -output score_ascii

  Generic common line for colored output:
  ##: t_coffee -other_pg seq_reformat -in sample_3Dseq1.aln -in2 sample_3Dseq1.template -action \
      +evaluate3D strike <max> <end> <matrix> -output score_html -out sample_3Dseq1.html



Evaluating a MSA according to your own criterion
================================================
Any kind of feature can easily be turned into an evaluation grid. For instance, the protease sequences we have been using here have a well characterized binding site. A possible evaluation can be made as follows: let us consider the UniProt annotation of the two distantly related sequences; these two sequences contain the electron relay system of the proteases. We can use it to build an evaluation library: in P29786 (TRY3_AEDAE) the Histidine residue is at position 68 while in P21844 (MCPT5_MOUSE) the functionally equivalent Histidine residue is at position 66. We can therefore build a library that will check whether these two residues are properly aligned in any MSA. The library will look like this:


::

  ! TC_LIB_FORMAT_01
  2
  sp|P21844|MCPT5_MOUSE 247 MHLLTLHLLLLLLGSSTKAGEIIGGTECIPHSRPYMAYLEIVTSENYLSACS\
  GFLIRRNFVLTAAHCAGRSITVLLGAHNKTSKEDTWQKLEVEKQFLHPKYDENLVVHDIMLLKLKEKAKLTLGVGTLP\
  LSANFNFIPPGRMCRAVGWGRTNVNEPASDTLQEVKMRLQEPQACKHFTSFRHNSQLCVGNPKKMQNVYKGDSGGPLL\
  CAGIAQGIASYVHRNAKPPAVFTRISHYRPWINKILREN
  sp|P29786|TRY3_AEDAE 254 MNQFLFVSFCALLDSAKVSAATLSSGRIVGGFQIDIAEVPHQVSLQRSGRHFC\
  GGSIISPRWVLTRAHCTTNTDPAAYTIRAGSTDRTNGGIIVKVKSVIPHPQYNGDTYNYDFSLLELDESIGFSRSIEA\
  IALPDASETVADGAMCTVSGWGDTKNVFEMNTLLRAVNVPSYNQAECAAALVNVVPVTEQMICAGYAAGGKDSCQGDS\
  GGPLVSGDKLVGVVSWGKGCALPNLPGVYARVSTVRQWIREVSEV
  #1 2
  66 68 100
  ! SEQ_1_TO_N


You simply need to cut and paste this library in a file and use this file to measure the concistency between your alignment and the correspondances declared in your library. The following command line also makes it possible to visually display the agreement between your sequences and the library.

::

  $$: t_coffee -infile proteases_small.aln -lib charge_relay_lib.tc_lib -score \
      -output html


*******************
Downstream Analysis
*******************
Contact Evaluations
===================
Contact Evaluation on proteins
------------------------------
This function will report in a library format all the intra-molecular distances of the provided sequences
::

  Command 1:
  $$: t_coffee -other_pg seq_reformat -in sample_3Dseq1.aln -in2 sample_3Dseq1.template -action +pdb2contacts intra  distances 10 -output contact_lib

This command will report all the distances between every pair of residue less than 10 Angstrom appart. In the output every pair of residue will be assocuiated with an integer value equal to 100*distance in Angstrom. Distances are between CA
 
::

Command 1:
  $$: t_coffee -other_pg seq_reformat -in sample_3Dseq1.aln -in2 sample_3Dseq1.template -action +pdb2contacts intra  contacts 1.2  -output contact_lib

This command will report all the pairs of residue featuring two atoms less than 1.2 Angstrom appart 
 



Clustering/Trees based on protein 3D structures
===============================================
This section describes tree estimation procedure based on the comparison of intramolecular distances. These methods are slightly different from the `T-RMSD <http://tcoffee.readthedocs.io/en/latest/tcoffee_main_documentation.html#the-t-rmsd>` that was also originally developped for this purpose and whose usage is decribed in the next section. The methods described in this section are more versatile than the original T-RMSD: they support trees based on contacts or intra-molecular distance variations,alternative tree reconstruction algorithms (NJ, UPGMA), genuine column based generalized boostrap procedures. They also support the STRIKE protocol that involves projecting contacts measured within one or more experimental structures onto the sequences these are aligned with and using the `STRIKE <https://academic.oup.com/bioinformatics/article/27/24/3385/306640/STRIKE-evaluation-of-protein-MSAs-using-a-single>` contact matrice to score these contacts in the sequences without known structures.   

Generating a tree based on 3D intramolecular distances conservation
-------------------------------------------------------------------

This option makes it possible to estimate a tree while taking into account the variation of intramolecular distances within the considered sequences. It requires in input an alignment (FASTA, MSF, ClustalW...) and a template file (structure files associated with the query sequences). The following command (command 1) will generate a 100 replicate NJ trees using the difference of distances between pairs of aligned residues, at a maximum cut-off of 15A. Columns with less than 50% residues are ignored. It is possible to control the defautl parameters (command 2). The ouput will be a tree in the Newick format with bootstrap supports.

::

  Command 1:
  $$: t_coffee -other_pg seq_reformat -in sample_3Dseq1.aln -in2 sample_3Dseq1.template -action \
      +tree replicates 100 +evaluate3D distances +tree2bs first +print_replicates -output newick


This command will print the main tree with bootstrap support values (first line). The 100 replicates (+tree replicates 100) will then be printed in the following lines (+print_replicates flag)  


::

  Command 2:
  $$: t_coffee -other_pg seq_reformat -in sample_3Dseq1.aln -in2 sample_3Dseq1.template -action \
      +tree replicates 100 +evaluate3D distances +tree2bs first +print_replicates -output dm



This command will print the main distance matrix in triangular format (first block). The 100 replicates (+tree replicates 100) will then be printed in the following blocks (+print_replicates flag)  


::

  Command 3:
  $$: t_coffee -other_pg seq_reformat -in sample_3Dseq1.aln -in2 sample_3Dseq1.template -action \
      +tree replicates 100 +evaluate3D distances +tree2bs first +print_replicates -output newick_dm


This command will print the newick trees and their corresponding distance matrix in triangular format (first block). The 100 replicates (+tree replicates 100) will then be printed in the following blocks along with the trees (+print_replicates flag). Note that the lines containing trees can be grepped by the ";" symbol they terminate with.  


:: 
  
  Command 4:
  $$: t_coffee -other_pg seq_reformat -in sample_3Dseq1.aln -in2 sample_3Dseq1.template -action \
      +tree replicates 100 gap 0.5 mode nj +evaluate3D distances 15 +tree2bs first -output newick
      


This command will do the same but pairs of residues will only contribute to the distance computation if they are less than 15 Angstroms appart (i.e. when comparing two aligned pairs, both pairs of residues in each sequence must be less than 15 Angstroms appart)

.. warning:: Sequences without 3D structure will be excluded from the analysis and from the final output.


Generating a tree based on contact conservation
-----------------------------------------------
The following option (command 1) makes it possible to estimate a tree while taking into account the variation of contact conservation within the considered sequences. This call will generate a 100 replicate NJ trees using as a distance metrics the fraction of contacts conserved between pairs of aligned residues, at a maximum cutoff of 1.2 A between Van der Waals radius and ignoring the 3 closest neighbors; columns with less than 50% residues are ignored. For sequences without 3D information, the STRIKE contact potential is used instead (Watson and Crick base pairing propensity for RNA). The output will be a tree in Newick format. It is possible to control default parameters using the following extended command line (command 2).

:: 

  Command 1:
  $$: t_coffee -other_pg seq_reformat -in sample_3Dseq1.aln -in2 sample_3Dseq1.template -action \
      +tree replicates 100 +evaluate3D contacts +tree2bs first -output newick -out tree.dnd

  Command 2:
  $$: t_coffee -other_pg seq_reformat -in sample_3Dseq1.aln -in2 sample_3Dseq1.template -action \
      +tree replicates 100 gap 0.5 mode nj +evaluate3D contacts 1.2 3 +tree2bs first -output \
      newick -out tree.dnd



.. warning:: The procedure requires at least 1 sequence with a known 3D structure or with contact information.

Visualizing contact/distance conservation
-----------------------------------------
This same procedure can be used to visualize either intramolecular distance conservation or contact conservation.  The output option **score_raw** generate a tabulated dump of the numerical values associated with every residue, sequence and column of the considered alignment. The flag **-out** specifies the name of the output file containing the results (you can give the name you want, but if you read the documentation so far, you already know it...). 

::

  $$: t_coffee -other_pg seq_reformat -in sample_3Dseq1.aln -in2 sample_3Dseq1.template \
      -action +evaluate3D distances -output score_html -out out.html
      
  $$: t_coffee -other_pg seq_reformat -in sample_3Dseq1.aln -in2 sample_3Dseq1.template \
      -action +evaluate3D distances -output score_ascii -out out.ascii
      
  $$: t_coffee -other_pg seq_reformat -in sample_3Dseq1.aln -in2 sample_3Dseq1.template \
      -action +evaluate3D distances -output score_raw -out out.tab


Visualizing informative positions 
---------------------------------
If you have a well defined subgroup of sequences (domains having the same function, same specificity, etc...), it is possible to estimate which columns yield the best support using the following command. The input are an alignment of your sequences, a template file containing the list of structure files to use as templates and a FASTA file of the sequences that form the group whose support you want to analyze. The output will be a colored version of your MSA indicating the sequences that best contribute to your clustering.

::

  $$: t_coffee -other_pg seq_reformat -in sample_3Dseq1.aln -in2 sample_3Dseq1.template \
      -action +tree replicates columns +evaluate3D  distances +evaluateTree group_3Dseq1.fasta \
      -output score_html -out out_aln.html


Evaluating clustering capacities (under maintenance...)
--------------------------------
If you want to check the capacity of an algorithm to bring related sequences within monophyletic groups, you should name your sequences according to the group they belong to (XXXX_1 for members of group 1; YYYY_2, for members of group 2; etc...) and use the following evaluation procedure. The output will be the number of monophyletic groups containing sequences belonging to the same group. The tree can be precomputed (command 1) or it can be computed on the fly (command 2).

:: 

  Command 1: precomputed tree
  ##: t_coffee -other_pg seq_reformat -in <tree> +tree2collapse groups 4 +print nseq -output no

  Command 2: computed on the fly
  ##: t_coffee -other_pg seq_reformat -in <aln> -in2 <template> -action +tree replicates 100 \
      +evaluate3D  distances 15 +tree2bs first +tree2collapse groups 4 +print nseq -output no

Structural Tree Parameters
--------------------------
The structural tree parameters come through three distinct flags:
 - **-tree**                       :Used to pass Tree computation parameters 
 - *** replicates <I=integer> ***  : Indicates how many replicates I are to be produced (def=1) when doing the bootsrap
 - *** mode <nj|upgma> ***         : Tree computation mode to be used
 - *** gap  <F=float>***           : ignores columns with a fraction of gaps higher or equal to F (def=0.5)
 
 - **-evaluate3D contacts <MaxD> <Ng>**:Used to estimate the cointact matrix
 - ***<MAxD=float>***              : maximum distance between the closest atoms (including side chains) of two residues <def=1.3>
 - ***<Ng=integer>***              : Number of excluded neighbours (on both side) <def=3>	 
 - ***def***                       : Any value can be left to its default as follows -evaluate3D contacts def 5 will used 1.2 as maxD and 5 as neighbor cutoff

 - **-evaluate3D distances <MaxD> <Ng>**:Used to estimate the cointact matrix
 - ***<MAxD=float>***              : maximum distance between the closest atoms (including side chains) of two residues <def=15>. When set to 0, all against all distances are computed. 
 - ***<Ng=integer>***              : Number of excluded neighbours (on both side) <def=3>	 
 - ***def***                       : Any value can be left to its default as follows -evaluate3D contacts def 5 will used 1.2 as maxD and 5 as neighbor cutoff

 - **-evaluate3D strike <MaxD> <Ng> <mat>**:Used to estimate the cointact matrix
 - ***<MAxD=float>***              : maximum distance between the closest atoms (including side chains) of two residues <def=15>
 - ***<Ng=integer>***              : Number of excluded neighbours (on both side) <def=3>	 
 - ***<mat=string>***              : strike matrix file (def="strike")
 - ***def***                       : Any value can be left to its default as follows -evaluate3D contacts def 5 will used 1.2 as maxD and 5 as neighbor cutoff


Distance based trees measures are based on a normalized difference of distances between pairs of residues:
::
	X AB
	Y CD
	d1=d(A,B)
	d2=d(B,B)
	score(AB, CD)=(Min (d1/d2, d2/d1)^3)*d1
	score (X,Y)=Sum(score(AB,CD))+Sum(score (CD, AB)/Sum (d(AB))+Sum (d(CD)) over all pairs for which d(AB)>=MaxD AND d(CD)>=MaxD


The T-RMSD 
==========
What is the T-RMSD?
-------------------
T-RMSD (Tree based on Root Mean Square Deviation) is a structure based clustering method using the iRMSD to drive the structural clustering of your aligned sequences with available structures. The T-RMSD supports all the parameters supported by iRMSD or APDB. T-RMSD  is a distance RMSD (dRMSD) based method which generate Structural Trees (analogue to phylogenetic tree) to determine fine-grained structural variations associated with a given Multiple Sequence Alignment (generated by the user using the method of is choice). The specificity of T-RMSD compared to other structural comparison methods stems from its capacity to generate a structural tree with values equivalent bootstrap values supporting the structural clustering. Such clustering is achieved by the construction of matrixes of distances, calculated between equivalent residues as defined by the ungapped columns of the given MSA. The resulting matrixes are then combined using the CONSENSE program from the Phylip package to generate a consensus structural tree with equivalent bootstrap values supporting each node (from 0 to 100; 100 indicating that all positions support the clustering). 

Generating a structural tree with support values
------------------------------------------------
T-RMSD is a special modof T-Coffee. To run T-RMSD, you just need a MSA (generated the way you want) and a template file; you need one structure for each sequence, otherwise these sequences will be excluded from the final results. There are several output files: 

 - ``<input name>.struc_tree.list``: list of all individual trees (one per ungapped column)
 - ``<input name>.struc_tree.consense_output``: basic consensus tree rendering
 - ``<input name>.struc_tree.consensus``: resulting structural tree with support (Newick)
 - ``<input name>.struc_tree.html``: colored MSA according to the contribution to the clusterin#
 
::

  $$: t_coffee -other_pg trmsd -aln sample_3Dseq1.aln -template_file sample_3Dseq1.template
  

.. hint:: Another important information displayed only on screen is the number of usable positions used to build the clustering. If your alignment is too gappy, this percentage will be low and the resulting clustering may not be accurate nor informative; in the previous example, 28.57% of the MSA will be used. 

Visualizing the structural clustering
-------------------------------------
T-Coffee have limited tree visualization capacities, yet the T-RMSD delivers two relevant files in this matter: the first one is a basic rendering of the tree with the support values for each node (``<input name>.struc_tree.list``), the second one is a colored html of the MSA where columns are colored according to their individual support to the final consensus tree (``<input name>.struc_tree.html``). If you want to visualize/modify the resulting tree, we recommend to use a dedicated tree viewer such as `PhyloWidget <http://www.phylowidget.org/>`_. By default PhyloWidget is used by the T-RMSD web server, but as the format is standard (Newick format) you can use any viewer you want (iTOL, ETE, FigTree, Phylodendron, etc...)

The TCS
=======
The TCS (Transitive Consistency Score) is a measure of consistency between a multiple sequence alignmment and a library of pairwise alignments of the same sequences. In the original `publication <https://academic.oup.com/mbe/article/31/6/1625/2925802/TCS-A-New-Multiple-Sequence-Alignment-Reliability>` the TCS score was shown to correlate well with local structural accuracy and with the phylogenic reconstruction potential of individual MSA columns. As such the TCS can therefore be used to down-weight or filter out the less reliable positions. The TCS is available as `>web-server http://tcoffee.crg.cat/apps/tcoffee/do:core>` and its command line usage is described in details in an `earlier section <http://tcoffee.readthedocs.io/en/latest/tcoffee_main_documentation.html#transitive-consistency-score-tcs>` of this documentation. 

It can be used to evaluate and compare existing alignmnets 
. The original publication describing it is available
described in the previous subsection as it is mainly an evaluation tool. It has however several option that are relevant for downstream analysis such weighting and/or trimming MSA for phylogenetic reconstruction. It was demonstrated that TCS was able to improve phylogenetic reconstruction, you are welcome to look at the publication for more details.


*************************
Internal/External Methods 
*************************
.. note:: The real power of T-Coffee is its ability to seamlessly combine many methods into one. While we try to integrate as many methods as we can in the default distribution, we do not have the means to be exhaustive and if you desperately need your favorite method to be integrated, you will need to bite the bullet ...


What are the methods already integrated in T-Coffee?
====================================================
Although, it does not necessarily do so explicitly, T-Coffee always end up combining libraries (collections of pairs of residues). Given a set of libraries, T-Coffee makes an attempt to assemble the alignment with the highest level of consistence. You can think of the alignment as a timetable, each library pair being a request from students or teachers, and the job of T-Coffee would be to assemble the time table that makes as many people as possible happy...In T-Coffee, methods replace the students/professors to generate constraints. These methods can be any standard/non standard alignment methods that can be used to generate alignments (pairwise, most of the time). These alignments can be viewed as collections of constraints that must be fit within the final alignment. Of course, the constraints do not have to agree with one another...

This section shows you what are the available method in T-Coffee, and how you can add your own methods either through direct parameterization or via a perl script. There are two kinds of methods: the internal and the external. For the internal methods, you simply need to have T-Coffee up and running. The external methods are generally installed with T-Coffee (when using the installer), but if you have problems with some packages, refer to the **T-Coffee Installation** chapter.

INTERNAL methods (built-in)
---------------------------
Internal methods can be requested using the following names:

- **proba_pair**: adapted from ProbCons, this method (the current default) uses a pair HMM to compute a pairwise alignment with a biphasic gap penalty.

- **fast_pair**: makes a global fasta style pairwise alignment. For proteins: **matrix=blosum62mt, gep=-1, gop=-10, ktup=2**. For DNA, **matrix=idmat (id=10), gep=-1, gop=-20, ktup=5**. Each pair of residue is given a score function of the weighting mode defined by **-weight**.

- **slow_pair**: identical to fast pair, but does a full dynamic programming, using the myers and miller algorithm. This method is recommended if your sequences are distantly related.

- **ifast_pair**: iterative fast_pair.

- **islow_pair**: makes a global fasta alignmnet using the previously computed pairs as a library. `i` stands for iterative. Each pair of residue is given a score function of the weighting mode defined by **-weight**. The library used for the computation is the one computed before the method is used. The result is therefore dependent on the order in methods and library are set via the **-in** flag.

- **align_pdb_pair**: uses the **align_pdb** routine to align two structures. The pairwise scores are those returnes by the **align_pdb** program. If a structure is missing, fast_pair is used instead. Each pair of residue is given a score function defined by align_pdb. [UNSUPPORTED]

- **lalign_id_pair**: uses the ten top non intersecting local alignments, as delivered by lalign. Each alignement is weighted with its average percent identity.

- **lalign_rs_s_pair**: Same as above but does also does self comparison and uses the lalign raw_score (s stands for self). This is needed when extracting repeats. [UNSUPPORTED]

- **Matrix Amy**: matrix can be requested, simply indicate as a method the name of the matrix preceded with an X (i.e. Xpam250mt). If you indicate such a matrix, all the other methods will simply be ignored, and a standard fast progressive alignment will be computed. If you want to change the substitution matrix used by the methods, use the **-matrix** flag.

- **cdna_fast_pair**: this method computes the pairwise alignment of two cDNA sequences. It is a fast_pair alignment that only takes into account the aminoacid similarity with different penalties for insertions and frameshifts. This alignment is turned into a library where matched nucleotides receive a score equal to the average level of identity at the aminoacid level. This mode is intended to clean cDNA obtained from ESTs or to align pseudo-genes. [UNSUPPORTED]

EXTERNAL methods (plug-in)
--------------------------
External methods correspond to packages developed by other groups that you may want to run within T-Coffee. We are very open to extending these options and we welcome any request to add an extra interface. To have a complete list of the methods that can be used as plug-ins, just type **t_coffee** in a terminal, it will we displayed on your screen. Most of these methods can be used as either pairwise (**<method>_pair**) or multiple alignment methods (**<method>_msa**); note that all these methods use Blosum62 as a default. 

One package is a bit different; **fugue_pair** uses a standard FUGUE installation to make a sequence/structure alignment. Installation must be standard but it does not have to include all the FUGUE packages but only:

::

  Required packages:
  1) joy, melody, fugueali, sstruc, hbond
  2) copy fugue/classdef.dat /data/fugue/SUBST/classdef.dat

  Configuration:
  ##: setenv MELODY_CLASSDEF=<location>
  ##: setenv MELODY_SUBST=fugue/allmat.dat


Modifying the parameters of INTERNAL/EXTERNAL methods
=====================================================
Internal methods
----------------
It is possible to modify on the fly the parameters of hard coded methods (**EP** stands for Extra parameters). These parameters will superseed any other parameters.

::

  $$: t_coffee sample_seq1.fasta -method slow_pair@EP@MATRIX@pam250mt@GOP@-10@GEP@-1


External methods
----------------
External methods receive a command line built with the information provided via the parameter file (see next heading). It is possible to produce such a parameter file and to modify it in order to modify the commands passed to the methods. The passed command is built as follows:

::

  ##: <EXECUTABLE><PARAM1><IN_FLAG><seq_file><PARAM2><OUT_FLAG><outname><PARAM>


You should know what is the best place for squizing your extra parameters. It will depend on the application, although PARAM2 is usually a good guess. Now if you want, for instance to modify the gap penalty of clustalw, you can try the following (of course, you must know the command line of the program you are trying to modify (ClustalW in this case):


::

  $$: t_coffee sample_seq1.fasta -method clustalw_msa@EP@PARAM2@-GAPOPEN%e100%s-GAPEXT%e10
  
  <@EP>     : indicates that you will pass an extra parameter
  <@PARAM1> : is the name of this parameter
  <%s>      : replaces spaces
  <%e>      : replaces the equal sign


Integrating EXTERNAL methods
============================
If the method you need is not already included in T-Coffee, you will need to integrate it yourself. We give you here some guidelines on how to do so.

Accessing external methods
--------------------------
A special method exists in T-Coffee that can be used to invoke any existing program (command 1): for instance, ClustalwWis a method that can be ran with the command 2. Here is just an example but ClustalW can be replaced with any method using a similar syntax. If the program you want to use cannot be run this way, you can either write a perl wrapper that fits the bill or write a ``tc_method`` file adapted to your program (cf. next section). This special method (**em** for external method) uses a particular syntax (command 3).

::

  Command 1: T-Coffee using ClustalW
  $*: t_coffee sample_seq1.fasta -method=em@clustalw@pairwise

  Command 2: ClustalW command
  ##: method -infile=<infile> -outfile=<outfile>

  Command 3: Syntax for EXTERNAL methods
  ##: em@<method>@<aln_mode:pairwises_pairwise|multiple>

Customizing an external method 
------------------------------
T-Coffee can run external method using a ``tc_method`` file that can be used in place of an established method. Two such files are incorporated in T-Coffee (``clustalw_method.tc_method`` and ``generic_method.tc_method``). You can dump them and customize them according to your needs. The first file is a very straightforward example on how to run Clustalw via T-Coffee with a set of parameters you may be interested in. Note that **ALN_MODE** instructs T_Coffee to run ClustalW on every pair of sequences. The second file is detailed in the next section.

::

  1) Getting the configuration file:
  $$: t_coffee -other_pg unpack_clustalw_method.tc_method

  2) Format of the configuration file:
  *TC_METHOD_FORMAT_01

  ***************clustalw_method.tc_method*********
  EXECUTABLE clustalw
  ALN_MODE pairwise
  IN_FLAG -INFILE=
  OUT_FLAG -OUTFILE=
  OUT_MODE aln
  PARAM -gapopen=-10
  SEQ_TYPE S
  *************************************************
  
  3) The configuration file will cause T-Coffee to emit the following system call:
  ##: clustalw -INFILE=tmpfile1 -OUTFILE=tmpfile2 -gapopen=-10

  4) Running ClustalW via T-Coffee (in your working DIR):
  $$: t_coffee sample_seq1.fasta -method clustalw_method.tc_method
  
Advanced method integration
---------------------------
It may sometimes be difficult to customize the program you want to use through a ``tc_method file``. In that case, you may rather use an external perl_script to run your external application. This can easily be achieved using the ``generic_method.tc_method`` file which contains many hints on how to customize your new method. The ``tc_method`` files are treated like any standard established method in T-Coffee.

::
 
  1) Using any generic method
  $$: t_coffee -other_pg unpack_generic_method.tc_method

  2) Format of the configuration file:
  *TC_METHOD_FORMAT_01
  ***************generic_method.tc_method*********
  EXECUTABLE tc_generic_method.pl
  ALN_MODE pairwise
  IN_FLAG -infile=
  OUT_FLAG -outfile=
  OUT_MODE aln
  PARAM -method clustalw
  PARAM -gapopen=-10
  SEQ_TYPE S
  *************************************************
  * Note: &amp;bsnp can be used to for white spaces

  3) Run the method: 
  $*: t_coffee sample_seq1.fasta -method generic_method.tc_method


T-Coffee runs the script **tc_generic_method.pl** on your data. It also provides the script with parameters. In the case **-method clustalw** indicates that the script should run ClustalW on your data. Over the time, this script will be the place where novel methods will be integrated and make it possible to run any available method. It will be used to run the script **tc_generic_method.pl**, a perl script automatically generated by T-Coffee.

.. note:: If there is a copy of that script in your local directory, that copy will be used in place of the internal copy of T-Coffee.

The mother of all method files...
---------------------------------
Here is the mother of all configuration file for T-Coffee: 

::

  *TC_METHOD_FORMAT_01
  ******************generic_method.tc_method*************
  *
  * Incorporating new methods in T-Coffee
  * Cedric Notredame 17/04/05
  *
  *******************************************************
  *This file is a method file
  *Copy it and adapt it to your need so that the method
  *you want to use can be incorporated within T-Coffee
  ******************************************************
  * USAGE *
  *******************************************************
  *This file is passed to t_coffee via -in:
  *
  * t_coffee -in Mgeneric_method.method
  *
  *The method is passed to the shell using the following
  *call:
  *<EXECUTABLE><IN_FLAG><seq_file><OUT_FLAG><outname><PARAM>
  *
  *Conventions:
  *<FLAG_NAME>  <TYPE> <VALUE>
  *<VALUE>: no_name  <=> Replaced with a space
  *<VALUE>: &amp;nbsp <=> Replaced with a space
  *
  *******************************************************
  * EXECUTABLE *
  *******************************************************
  *name of the executable
  *passed to the shell: executable
  *
  EXECUTABLE tc_generic_method.pl
  *
  *******************************************************
  * ALN_MODE *
  *******************************************************
  *pairwise ->all Vs all (no self )[(n2-n)/2aln]
  *m_pairwise ->all Vs all (no self)[n^2-n]^2
  *s_pairwise ->all Vs all (self): [n^2-n]/2 + n
  *multiple ->All the sequences in one go
  *
  ALN_MODE pairwise
  *
  *******************************************************
  * OUT_MODE *
  *******************************************************
  *mode for the output:
  *External methods:
  * aln -> Alignmnent file (Fasta or ClustalW Format)
  * lib-> Library file (TC_LIB_FORMAT_01)
  *Internal Methods:
  * fL -> Internal Function returning a Lib (Library)
  * fA -> Internal Function returning an Alignmnent
  *
  OUT_MODE aln
  *
  *******************************************************
  * IN_FLAG *
  *******************************************************
  *IN_FLAG
  *flag indicating the name of the in coming sequences
  *IN_FLAG no_name ->no flag
  *IN_FLAG &amp;nbsp-in&amp;nbsp -> ' -in '
  *
  IN_FLAG -infile=
  *
  *******************************************************
  * OUT_FLAG *
  *******************************************************
  *OUT_FLAG
  *flag indicating the name of the out-coming data
  *same conventions as IN_FLAG
  *OUT_FLAG no_name ->no flag
  *
  OUT_FLAG -outfile=
  *
  *******************************************************
  * SEQ_TYPE *
  *******************************************************
  *G: Genomic, S: Sequence, P: PDB, R: Profile
  *Examples:
  *SEQTYPE S sequences against sequences (default)
  *SEQTYPE S_P sequence against structure
  *SEQTYPE P_P structure against structure
  *SEQTYPE PS mix of sequences and structure
  *
  SEQ_TYPE S
  *
  *******************************************************
  * PARAM *
  *******************************************************
  *Parameters sent to the EXECUTABLE
  *If there is more than 1 PARAM line, the lines are
  *concatenated
  *
  PARAM -method clustalw
  PARAM -OUTORDER=INPUT -NEWTREE=core -align -gapopen=-15
  *
  *******************************************************
  * END *
  *******************************************************

Weighting your method
---------------------
By default, the alignment produced by your method will be weighted according to its percent identity. However, this can be customized via the **WEIGHT** parameter, which supports all the values of the **-weight** flag. The only difference is that the **-weight** value thus declared will only be applied onto your method. If needed you can also modify on the fly the WEIGHT value of your method:

::

  Increases by a factor 2 the weight of slow_pair:
  $$: t_coffee sample_seq1.fasta -method slow_pair@WEIGHT@OW2
 
  Causes every pair of slow_pair to have a weight equal to 250:
  $$: t_coffee sample_seq1.fasta -method slow_pair@WEIGHT@250

Managing a collection of method files
-------------------------------------
It may be convenient to store all the method files in a single location on your system. By default, T-Coffee will go looking into the directory ~/.t_coffee/methods/. You can change this by either modifying the **METHODS_4_TCOFFEE** in ``define_headers.h`` (and recompile) or by modifying the environement variable **METHODS_4_TCOFFEE****.


Creating your own T-Coffee libraries
====================================
If the method you want to use is not integrated or impossible to integrate, you can generate your own libraries, either directly or by turning existing alignments into libraries. You may also want to precompute your libraries, in order to combine them at your convenience.

Using precomputed alignments
-----------------------------
If the method you wish to use is not supported, or if you simply have the alignments, the simplest thing to do is to generate yourself the pairwise/multiple alignments, in FASTA, ClustalW, MSF or PIR format and feed them into T-Coffee. You can even combine multiple MSAs (refers to subsection **Combining multiple MSAs**).

Customizing the weighting scheme
--------------------------------
The previous integration method forces you to use the same weighting scheme for each alignment and the rest of the libraries generated on the fly. This weighting scheme is based on global pairwise sequence identity. If you want to use a more specific weighting scheme with a given method, you should either generate your own library (command 1) or convert your alignment into a library (command 2), or use the **-weight** flag (command 3).

::

  Command 1:
  $$: t_coffee -aln=sample_aln1.aln,sample_aln2.aln -method=fast_pair,lalign_id_pair

  Command 2:
  $$: t_coffee -aln sample_aln1.aln -lib=test_lib.tc_lib 
  
  Command 3:
  $$: t_coffee -aln sample_aln1.aln -out_lib=test_lib.tc_lib -lib_only -weight=sim_pam250mt
      

Generating your own libraries
-----------------------------
This is suitable if you have local alignments, or very detailed information about your potential residue pairs, or if you want to use a very specific weighting scheme. You will need to generate your own libraries using the format described in the last section. You may also want to precompute your libraries in order to save them for further use. For instance, in the following example, we generate the local (command 1) and the global (command 2) libraries and later reuse them for combination into a multiple alignment. Once these libraries have been computed, you can then combine them at your convenience in a single MSA (command 3). Of course you can decide to only use the local or the global library.

::

  Command 1:
  $$: t_coffee sample_seq1.fasta -method lalign_id_pair -out_lib lalign_id_pair_seq1.tc_lib \
      -lib_only

  Command 2:
  $$: t_coffee sample_seq1.fasta -method slow_pair -out_lib slow_pair_seq1.tc_lib -lib_only

  Command 3:
  $$: t_coffee sample_seq1.fasta -lib lalign_id_pair_seq1.tc_lib, slow_pair_seq1.tc_lib


Plug-out: using T-Coffee as a plug-in
=====================================
Just because it enjoys enslaving other methods as plug-in does not mean that T-Coffee does not enjoy being incorporated within other packages. We try to give as much support as possible to anyone who wishes to incorporate T-Coffee in an alignment pipeline. If you want to do so, please work out some way to incorporate T-Coffee in your script. If you need some help along the way, do not hesitate to ask, as we will always be happy to either give assistance or even modify the package so that it accomodates as many needs as possible. Once that procedure is over, set aside a couple of input files with the correct parameterisation and send them to us. These will be included as a distribution test, to insure that any further distribution remains compliant with your application.
