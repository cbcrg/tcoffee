###########################
T-Coffee Main Documentation
###########################
*******************
Before You Start...
*******************

.. Danger:: T-Coffee Manual is currently under heavy maintenance in order to update all the latest T-Coffee developments and tools but also to increase its readibility. For now, the example files can be found here on github: `example files <https://github.com/cbcrg/tcoffee/tree/master/t_coffee/doc_test/data>`_. You can also use examples associated with their corresponding command lines from the section T-Coffee Tutorial published in Nature Protocols (2011). This message will be removed as soon as the maintenance is over. 


Foreword
========
Most of the work presented here emanates from two summer schools that were tentatively called the 'Prosite Workshops' and were held in Marseille, in 2001 and 2002. These workshops were mostly an excuse to go rambling and swimming in the creeks of Cassis (Calanques). Yet, when we got tired of lazing in the sun, we eventually did a bit of work to chill out. Most of our experiments were revolving around the development of sequence analysis tools. Many of the most advanced ideas in T-Coffee were launched during these fruitful sessions. Participants included Phillip Bucher, Laurent Falquet, Marco Pagni, Alexandre Gattiker, Nicolas Hulo, Christian Siegfried, Anne-Lise Veuthey, Virginie Leseau, Lorenzo Ceruti and Cedric Notredame.


.. Note:: This manual is made to help you discover (nearly) all subtleties of T-Coffee and ranges from simple applications to more complicated ones. We have tried to put as many of these functionalities on the webserver but not all are available except by using "command lines", also, **if you need to do something sophisticated and/or highly reproducible, the command line is the only way**. In the course of this manual, we expect you to use a Unix-like command line shell. All the example files can be found `here <https://github.com/cbcrg/tcoffee/tree/master/t_coffee/doc_test/data>`_. If you don't want to bother with command line, try using our online latest `tcoffee webserver <http://tcoffee.crg.cat/apps/tcoffee/index.html>`_ or via the `Cedric Notredame's lab homepage <http://www.tcoffee.org>`_.


Prerequisite for using T-Coffee
===============================
This manual relies on the assumption that you have installed T-Coffee, version 6.18 or higher. All the relevant information for installing T-Coffee is contained in the previous section **T-Coffee Installation**. T-Coffee is a freeware open source running on all Unix-like platforms, including Mac OS X and Cygwin. T-Coffee cannot run on Windows except by using Cygwin, a freeware open source allowing to run a Unix-like command line on Windows (`Cygwin download <https://www.cygwin.com/>`_). Cygwin is free of charge and easy to install, yet, as the first installation requires downloading substantial amounts of data, make sure you have access to a broad-band connection. A better option and a more stable one would be to install a Unix-like virtualbox on your computer using Windows.


Have a try...it's not painful!
==============================
We encourage you to try the following examples with your own sequences/structures. If you want to try with ours, they are available on our github repository : `example files <https://github.com/cbcrg/tcoffee/tree/master/t_coffee/doc_test/data>`_. It contains everything you need to run ALL the examples provided in this manual.



*******************
What is  T-Coffee ?
*******************
What is T-Coffee?
=================
**T-Coffee stands for Tree based Consistency Objective Function For alignment Evaluation**, it is primarily a multiple sequence alignment method but it also provides a collection of useful tools presented in this manual. Before going deep into the core of the matter, here are a few words to quickly explain some of the things T-Coffee will do for you.


What does it do?
----------------
T-Coffee is mainly a Multiple Sequence Alignment program: given a set of sequences previously gathered using database search programs like BLAST, FASTA or Smith and Waterman, T-Coffee will produce a Multiple Sequence Alignment (MSA). To use T-Coffee you must already have your sequences ready. T-Coffee is also bundled with a collection of tools (usually called with the flag **-other_pg "tools"**) and third party software (cf. section **Integrating External Methods in T-Coffee**) allowing you to perform a wide range of different tasks, such as aligning, reformatting, evaluating your data, and more...


What can it align?
------------------
T-Coffee will align nucleic acid (DNA and RNA) and protein sequences alike. T-Coffee is also able to use other type of information such as secondary and tertiary structure information (for protein or RNA sequences with a known/predicted structure) or sequence profiles.

.. Hint:: To have an idea, on an average PC, T-Coffee can easily align up to a 100-200 sequences, about 1000 amino acid long. 


How can I use it?
-----------------
T-Coffee is not an interactive program. It runs from your Unix or Linux command line and you must provide it with the correct parameters and syntax. If you do not like typing commands, you can still use the latest T-Coffee webserver `here <http://tcoffee.crg.cat/apps/tcoffee/index.html>`_.

.. Tip:: Installing and using T-Coffee requires a minimum acquaintance with the Linux/Unix operating system. If you feel this is beyond your computer skills, we suggest you use one of the available online servers.


Is there an online webserver?
-----------------------------
Yes, at `Cedric Notredame's lab homepage <http://www.tcoffee.org>`_ which contains all the necessary links to our latest  `webserver <http://tcoffee.crg.cat/apps/tcoffee/index.html>`_ !


Is T-Coffee different from ClustalW?
------------------------------------
According to several benchmarks, T-Coffee is on overall much more accurate than ClustalW, but this increase in accuracy comes at a price: **T-Coffee (default mode) is slower than ClustalW** (about N times for N Sequences). Still, if you want to align closely related sequences, **T-Coffee can also be used in a fast mode ,much faster than ClustalW, and about as accurate** (cf. section **Building Multiple Sequence Alignment**). If you are familiar with ClustalW or if you run a ClustalW server, you will find that we have made some efforts to ensure as much compatibility as possible between ClustalW and T-Coffee. Whenever it was relevant, we have kept the flag names and the flag syntax of ClustalW. Yet, you will find that T-Coffee also has many extra possibilities...


Is T-Coffee very accurate?
--------------------------
T-Coffee belongs to the category of consistency-based aligners which currently corresponds to the most accurate algorithms available (e.g. ProbCons, MSAprobs...). In addition, T-Coffee can combines (many) methods and therefore be as accurate (and hopefully more) as the methods it combines. For instance, the "accurate" mode of T-Coffee is very slow but also very accurate; on average this mode was shown to be 10 % more accurate than normal aligners on sequences less than 30% similar. If you need a very accurate alignment go to section **Building Multiple Sequence Alignment**.



What T-Coffee can and cannot do for you ...
===========================================

What T-Coffee can't do
----------------------
To be honest, a short answer will be that there is only one thing T-Coffee cannot do for you: **T-Coffee can NOT fetch sequences for you**. You must select the sequences you want to align beforehand and prepare your own dataset. We suggest you use any BLAST server and format your sequences in FASTA so that T-Coffee can use them easily. The  `ExPASy BLAST server <http://www.expasy.ch>`_ provides a nice interface for integrating database searches.


What T-Coffee can do
--------------------
T-Coffee is not only just an aligner program, it comes with multiple tools and third party software increasing the range of its possibilities; here is a non exhaustive list of tasks T-Coffee can perform:

**T-Coffee can compute (or at least try to compute!) accurate Multiple Sequence Alignments of DNA, RNA or Protein sequences**. Several modes and options are available and will be presented all along this manual. The default T-Coffee accepts any kind of sequence, although some modes are specific to a given type of sequence.

**T-Coffee can help you to reformat, trim, clean, cut, color your input (sequences, structures...) or output (alignments, trees...) data**; meaning that once you have your data and/or results ready, you can always modify them at will.

**T-Coffee allows you to combine results obtained with several alignment methods** (see the section **FAQ for T-Coffee** and **Building Multiple Sequence Alignment** for more details). T-Coffee can virtually combine all these MSAs you have to produce a new Multiple Sequence Alignment having the best agreement with all these methods you tried.

**One of the most important improvement of T-Coffee is to let you combine sequences and structures**, so that your alignments are of higher quality. You need to have the SAP package installed to fully benefit of this facility (or to use another structural alignment method). 

.. warning:: In this particular T-Coffee mode, operations are carried out using the function "wget". If "wget" is not installed on your system, you can get it for free from  `wget download <http://www.wget.org>`_. To make sure wget is installed on your system, type **which wget**.

**T-Coffee allows you to extract a serie of repeats from a single sequence or a set of sequences** using MOCCA. In other words, if you know the coordinates of one copy of a repeat, you can extract all the other occurrences. 

.. note:: MOCCA needs some time to compute a library and then prompt you with an interactive menu. You just have to follow the instructions.

**T-Coffee can be used to measure the reliability of your Multiple Sequence Alignment**. If you want to find out about that, read the section **FAQ for T-Coffee** or the **Technical Documentation** (-output flag). More details will be given anyway in this manual in the section **How Good Is Your Alignment?**.

**T-Coffee can be used to compare alternative alignment**; in case you generate several alignments of the same sequences, you can compare these alignments using the most common scores (Sum-of-Pairs or Column Score). In case you have reference alignments, you can directly benchmark your method by comparing your MSAs to your references.


And probably many more options we will discover together all along this manual !


How does T-Coffee alignment works?
==================================
If you only want to make a standard Multiple Sequence Alignment, you may skip these explanations. But if you want to do more sophisticated things, these few indications may help before you start reading the documentation and the different articles. 

When you run T-Coffee, the first thing it does is to compute a library. The library is a list of pairs of residues that could be aligned...it is like a christmas list: you can ask anything you fancy, but it is down to Santa Claus to assemble a collection of toys that won't get him stuck at the airport, while going through the metal detector. 

Given a standard library, it is nearly impossible to have all the residues aligned at the same time because all the lines of the library may not agree. For instance:

::

  Line 1 says:
  Residue 1 of seq A with Residue 5 of seq B,
  ...
  ...
  Line 100 says:
  Residue 1 of seq A with Residue 29 of seq B,

Each of these constraints comes with a weight and in the end, the T-Coffee algorithm tries to generate the multiple alignment that contains constraints whose sum of weights yields the highest score. In other words, it tries to make happy as many constraints as possible (replace the word constraint with, friends, relatives, collaborators... and you will know exactly what we mean).

You can generate this list of constraints the way you like. You may even provide it yourself, forcing important residues to be aligned by giving them high weights (see **FAQ for T-Coffee**). For your convenience, T-Coffee can generate (this is the default) its own list by making all the possible global pairwise alignments, and the 10 best local alignments associated with each pair of sequences. Each pair of residues observed aligned in these pairwise alignments, becomes a line in the library.


.. note:: Be aware that nothing forces you to use a given library and that you could build it using other methods. In protein language, **T-Coffee is synonymous for freedom, the freedom of being aligned however you fancy** (I was a Tryptophan in some previous life).



********************************************************
Preparing Your Data: Reformatting, Trimming, and more... 
********************************************************
The reformatting utility: seq_reformat
======================================
General introduction
--------------------
Nothing is more frustrating than downloading important data and realizing you need to format it before using it. In general, you should avoid manual reformatting: it is by essence inconsistent and will get you into trouble. It will also get you depressed when you realize that you have spend the whole day adding carriage return to each line in your files. T-Coffee comes with several tools to reformat/trim/clean/select your input data but also your output results, especially a very powerful reformatting utility named **seq_reformat**. You can use seq_reformat by invoking the t_coffee shell:

::

  $$: t_coffee -other_pg seq_reformat


This will output the online flag usage of seq_reformat meaning a complete list of things seq_reformat can do for you. The seq_reformat is a reformatting utility so it recognizes automatically the most common formats (FASTA, Swiss-Prot,ClustalW, MSF, Phylip...). It reads the input file(s) via the **-in** and **-in2** flags and outputs in whatever specified format via the **-output** flag. In the meantime, you can use the flag **-action** to perform a wide range of modification on your data. In this section we give you quite a lot of different examples of you can do with **seq_reformat**.

.. danger:: After the flag -other_pg, the common T-Coffee flags are not recognized anymore; it is like if you were using a different program.

Modification options
--------------------
In order to perform different modifications on your data (residues/sequences/columns...), the seq_reformat utility has to be followed by the flag **-action** (within the seq_reformat tool) and one or several modifiers listed here (this list is not exhaustive):

:: 

  Options:
  - +upper          : to uppercase your residues
  - +lower          : to lowercase your residues
  - +switchcase     : to selectively toggle the case of your residues
  - +keep           : to only keep the residues within the range
  - +use_cons +keep : to only keep the columns within the range
  - +remove         : to remove the residues within the range
  - +convert        : to only convert the residues within the range
  - +grep           : to select a given string of character
  - +rm_gap         : to remove columns containing gaps
  - etc...
 
  
Using a "cache" file
--------------------
Several option can be performed easily by using what we call a cache (or cache file). A cache is a file containing an alternate version of your alignment where each position of the alignment is replaced by a score previously evaluated: this score can be the T-Coffee CORE index (cf. section **How Good Is Your Alignment?**) or a matrix-based evalution (blosum62nt or identity matrix). Then, when performing any modification or reformatting of your alignments, you can just specify the range of positions to be modified according to their respective scores within the cache. We will see some example especially regarding the modification of format of a given alignment; it is not mandatory to use a cache but it is rather practical. To generate a cache before any reformatting using a given evaluation score, you can use one of the following possible option:

::

  Evaluating the T-Coffee CORE index during the alignment procedure:
  $$: t_coffee sample_aln1.fasta -output=score_ascii

  Evaluating the T-Coffee CORE index of a given alignment:
  $$: t_coffee -infile sample_aln1.aln -mode evaluate

  Using an identity matrix:
  $$: t_coffee -other_pg seq_reformat -in sample_aln1.aln -action +evaluate \
      idmat -output score_ascii

  Using a substitution matrix:
  $$: t_coffee -other_pg seq_reformat -in sample_aln1.aln -action +evaluate \
      blosum62mt -output score_ascii
      
      
Modifying the format of your data
=================================
Changing the sequence format
----------------------------
Sometimes it may be necessary to change from one format to another, for instance when using another software which recognize only a given format. T-Coffee recognizes most common alignment formats and you can find the list of all input or output format recognized by simply typing:

::

  $$: t_coffee -other_pg seq_reformat 
  
It is possible to reformat unaligned or aligned sequences alike although changing the alignment format is probably more interesting in order to use other applications; unaligned sequences format flags are generally preceded by the suffix "_seq" and aligned sequences flags by the suffix "_aln". This also allows you to transform any alignment into unaligned sequences by removing the gaps. Here are some examples on how to change the format of your data:

::

  For unaligned sequences (e.g. FASTA to PIR):
  $$: t_coffee -other_pg seq_reformat -in sproteases_small.fa -output pir_seq >\
      sproteases_small.fasta_aln
  
  For alignements (e.g. ClustalW to MSF):
  $$: t_coffee -other_pg seq_reformat -in sproteases_small.aln -output fasta_aln >\
      sproteases_small.fasta_aln
      
  From aligned to unaligned sequences:
  $$: t_coffee -other_pg seq_reformat -in sproteases_small.aln -output fasta_seq >\
      sproteases_small.fa

.. Warning:: Format recognition is not 100% full proof; occasionally you will have to inform the program about the nature of the file you are trying to reformat with " -input msf_aln -output fasta_aln" for instance.


Changing the case
-----------------
Changing the case of your sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you need to change the case of your sequences, you can use different modifiers embedded in seq_reformat. They are accessed via the **-action** flag. For instance, to write your sequences in lower case:

::

  $$: t_coffee -other_pg seq_reformat -in sproteases_small.aln -action +lower\
      -output clustalw


.. hint:: No prize for guessing that +upper will do exactly the opposite...


Changing the case of specific residues
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you want to change the case of a specific residue, you can use the flag: +edit_residue <sequence> <residue #> <lower|upper|symbol>. If you have more than one residue to modify, write all the coordinates in a text file (one coordinate per line) as spans are not yet supported; then give the file to T-Coffee

::

  $$: t_coffee -other_pg seq_reformat -in sample_aln1.aln -action +upper \
      +edit_residue hmgb_chite 10 lower
      
  $$: t_coffee -other_pg seq_reformat -in sample_aln1.aln -action +upper \ 
      +edit_residue <your file containing coordinates>

.. warning:: If you give a list of coordinates, it has to be a Unix text file (not a word document).


Changing the case with a cache
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you want to change the case depending on the score, you must either evaluate your alignment, or provide a cache. For example, this command line will upper the case of all residue then lower the case of every residue more than 50% identical to other residues in the same column:

::

  Using a cache on-the-fly:
  $$: t_coffee -other_pg seq_reformat -in sample_aln7.aln -action +upper \
      +evaluate idmat +lower '[5-9]'
      
  Using a cache file previously computed (2 steps):
  $$: t_coffee -other_pg seq_reformat -in sample_aln7.aln -action +evaluate \
      idmat -output score_ascii > sample_aln7.cache
  $$: t_coffee -other_pg seq_reformat -in sample_aln7.aln -struc_in sample_aln7.cache \
      -struc_in number_aln -action +upper '[5-9]'
  

Keeping/Protecting your sequence names
--------------------------------------
Only few programs support long sequence names, and sometimes, when going through some pipeline the names of your sequences can be truncated or modified. To avoid this, seq_reformat contains a utility that can automatically rename your sequences into a form that will be machine-friendly, while making it easy to return to the human-friendly form.

1) **Create a code list**: The first thing to do is to generate a list of names that will be used in place of the long original name of the sequences:

::

  $$: t_coffee -other_pg seq_reformat -in sproteases_large.fasta -output \
      code_name > sproteases_large.code_name

2) **Code your data**: This will create a file where each original name is associated with a coded name (Cxxxx). You can then use this file to either code your dataset using the following command:

::

  $$: t_coffee -other_pg seq_reformat -code sproteases_large.code_name -in \
      sproteases_large.fasta > sproteases_large.coded.fasta

3) **Decode your data**: Then you can work with the file sproteases_large.coded.fasta and when you are done, you can decode the names of your sequences with the following command line:

::

  $$: t_coffee -other_pg seq_reformat -decode sproteases_large.code_name -in \
      sproteases_large.coded.fasta


Colouring/Editing residues in an alignment
------------------------------------------
Changing the default colors
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Colors are hard coded in the program, but if you wish, you can change them by simply creating a file named ``seq_reformat.color`` that is used to declare the color values. The name of the file (seq_reformat.color) is defined in programmes_define.h, COLOR_FILE and can be changed before compilation. By default, the file is searched in the current directory. For example, the following line written in ``seq_reformat.color`` indicates that the value 0 in the cache corresponds now to #FFAA00 in html, and in RGB 1, 0.2 and 0. 

::

  0 #FFAA00 1 0.2 0


Coloring specific types of residues/nucleic acids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
You can color all the residues of your sequences on-the-fly; for instance, the following command line will color all the a's in color 0 (blue):

::

  $$: t_coffee -other_pg seq_reformat -in sample_aln1.aln -action +3convert a0 \
      -output color_html > colored.html

.. warning:: This option is case sensitive so the case of the residues or nucleotides should be the same in the command line (in this command line, only a lower case will be colored). 


Coloring a specific residue of a specific sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you want to color a specific residue/nucleotide, you can use the flag **+color_residue <sequence> <residue #> <color #>**. If you have more than one residue to color, you can put all the coordinates in a file, (one coordinate per line). Spans are not yet supported.

::

  $$: t_coffee -other_pg seq_reformat -in sample_aln1.aln -action +color_residue \
      hmgb_chite 10 1 -output color_html > color.html

.. warning:: If you give a list of coordinates, it has to be a Unix text file (not a word document).

Coloring according to the conservation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Use the +evaluate flag if you want to color your alignment according to its conservation level:

::

  $$: t_coffee -other_pg seq_reformat -in sample_aln1.aln -in3 sample_aln1.aln -action \
      +3evaluate pam250mt- output color_html > color.html



You can also use the boxshade scoring scheme:

::

  $$: t_coffee -other_pg seq_reformat -in sample_aln1.aln -in3 sample_aln1.aln -action \
      +3evaluate boxshade -output color_html > color.html


Coloring/Editing residues in an alignment using a Cache
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In this particular case, the cache is a file where residues to be colored are declared along with the colors. Nine different colors are currently supported. They are set by default but can be modified by the user. The cache can either look like a standard sequence or alignment file (see below) or like a standard T-Coffee library (see next section). In this section we show you how to specifically modify your original sequences to turn them into a cache.


In the cache, the colors of each residue are declared with a number between 0 and 9. Undeclared residues will appear without any color in the final alignment.


Preparing a sequence or alignment cache
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Let us consider the following file:


::

  CLUSTAL FORMAT

  B CTGAGA-AGCCGC---CTGAGG--TCG
  C TTAAGG-TCCAGA---TTGCGG--AGC
  D CTTCGT-AGTCGT---TTAAGA--ca-
  A CTCCGTgTCTAGGagtTTACGTggAGT


The command


::

  $$: t_coffee -other_pg seq_reformat -in=sample_aln6.aln -output=clustalw_aln \
      -out=cache.aln -action +convert 'Aa1' '.--' +convert '#0'



The conversion will proceed as follows:


-conv indicates the filters for character conversion:


 The gaps "-" will remain
 A and a will be turned into 1
 All the other symbols (#) will be turned into 0.

-action +convert, indicates the actions that must be carried out on the alignment before it is output into cache.


This command generates the following alignment (called a cache):

::

  CLUSTAL FORMAT for SEQ_REFORMAT Version 1.00, CPU=0.00 sec, SCORE=0, Nseq=4, Len=27
  B 000101-100000---000100--000
  C 001100-000101---000000--100
  D 000000-100000---001101--01-
  A 000000000010010000100000100

Other alternative are possible. For instance, the following command:

::

  $$: t_coffee -other_pg seq_reformat -in=sample_aln6.aln -output=fasta_seq -out=cache.seq \
      -action +convert 'Aa1' '.--' +convert '#0'


will produce the following file cache_seq

::

  >B
  000101100000000100000
  >C
  001100000101000000100
  >D
  00000010000000110101
  >A
  000000000010010000100000100


where each residue has been replaced with a number according to what was specified by conv. Note that it is not necessary to replace EVERY residue with a code. For instance, the following file would also be suitable as a cache:


::

  $$: t_coffee -other_pg seq_reformat -in=sample_aln6.aln -output=fasta_seq -out=cache \
      -action +convert 'Aa1' '.--'

  >B
  CTG1G11GCCGCCTG1GGTCG
  >C
  TT11GGTCC1G1TTGCGG1GC
  >D
  CTTCGT1GTCGTTT11G1c1
  >A
  CTCCGTgTCT1GG1gtTT1CGTgg1GT


Preparing a library cache
^^^^^^^^^^^^^^^^^^^^^^^^^
The Library is a special format used by T-Coffee to declare special relationships between pairs of residues. The cache library format can also be used to declare the color of specific residues in an alignment. For instance, the following file


::

  ! TC_LIB_FORMAT_01

  4

  A 27 CTCCGTgTCTAGGagtTTACGTggAGT
  B 21 CTGAGAAGCCGCCTGAGGTCG
  C 21 TTAAGGTCCAGATTGCGGAGC
  D 20 CTTCGTAGTCGTTTAAGAca

  #1 1
   1 1 3
   4 4 5
  #3 3
   6 6 1
   9 9 4
   
  ! CPU 240
  ! SEQ_1_TO_N


sample_lib5.tc_lib declares that residue 1 of sequence 3 will be receive color 6, while residue 20 of sequence 4 will receive color 20. Note that the sequence number and the residue index are duplicated, owing to the recycling of this format from its original usage.


It is also possible to use the BLOCK operator when defining the library (c.f. technical doc, library format). For instance:

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


The number right after BLOCK indicates the block length (10). The two next numbers (1 1) indicate the position of the first element in the block. The last value is the color.


Coloring an alignment using a cache
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you have a cache alignment or a cache library, you can use it to color your alignment and either make a post script, html or PDF output. For instance, if you use the file cache.seq:

::

   $$: t_coffee -other_pg seq_reformat -in=sample_aln6.aln -struc_in=sample_aln6.cache \
       -struc_in_f number_fasta -output=color_html -out=x.html

This will produce a colored version readable with any standard web browser, while:

::

   $$: t_coffee -other_pg seq_reformat -in=sample_aln6.aln -struc_in=sample_aln6.cache \
       -struc_in_f number_fasta -output=color_pdf -out=x.pdf


This will produce a colored version readable with acrobat reader.


.. warning:: ps2pdf must be installed on your system

You can also use a cache library like the one shown above (sample_lib5.tc_lib):

::

  $$: t_coffee -other_pg seq_reformat -in=sample_aln6.aln -struc_in=sample_lib5.tc_lib \
      -output=color_html -out=x.html
 

Modifying the data itself...
=============================
Modifiying sequences in your dataset
------------------------------------
Converting residues
^^^^^^^^^^^^^^^^^^^
It is possible for instance to selectively convert all given characters in a sequence (residues or nucleic acids alike) into another one, for example all G's having a score between 1 and 2 by using the command line:

::

  $$: t_coffee -other_pg seq_reformat -in sample_aln7.aln -struc_in sample_aln7.cache_aln \ 
      -struc_in_f number_aln -action +convert '[1-2]' CX
 

Extracting sequences according to a pattern
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
You can extract any sequence by requesting a specific pattern to be found either in the name (NAME), the comment (COMMENT) or the sequence (SEQ) using the modifier is '+grep'. For instance, if you want to extract all the sequences whose name contain the word HUMAN, the flag NAME/COMMENT/SEQ indicates that the modification is made according to the sequences names, the comment section or the sequence itself, and the flag KEEP/REMOVE means that you will keep/remove all the sequences containing the string HUMAN. Here are some examples:

::

  To keep sequences containing HUMAN in the name:
  $$: t_coffee -other_pg seq_reformat -in sproteases_small.aln -action +grep NAME \
      KEEP HUMAN -output clustalw

  To remove sequences containing HUMAN in the name:
  $$: t_coffee -other_pg seq_reformat -in sproteases_small.aln -action +grep NAME \
      REMOVE HUMAN -output clustalw

  To keep sequence which contain sapiens in the comment:
  $$: t_coffee -other_pg seq_reformat -in sproteases_small.aln -action +grep COMMENT \
      KEEP sapiens -output clustalw
 
  To remove sequences containing the pattern [ILM]K:
  $$: t_coffee -other_pg seq_reformat -in sproteases_small.aln -action +grep SEQ \
      REMOVE '[ILM]K' -output clustalw


.. important:: you should know that the pattern can be any perl legal regular expression, you can visit this  `page <http://www.comp.leeds.ac.uk/Perl/matching.html>`_ for some background on regular expressions. 

.. caution:: This option is case sensitive (Human, HUMAN and hUman will not yield the same results). Be careful !!!


Extracting/Removing specific sequences by names
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you want to extract (command 1) or remove (command 2) several sequences in order to make a subset, you can specify a list of sequences by their full name:

::

  Command 1: keep sequences
  $$: t_coffee -other_pg seq_reformat -in sproteases_small.aln -action +extract_seq_list \
      'sp|P29786|TRY3_AEDAE' 'sp|P35037|TRY3_ANOGA'

  Command 2: remove sequences
  $$: t_coffee -other_pg seq_reformat -in sproteases_small.aln -action +remove_seq \
      'sp|P29786|TRY3_AEDAE' 'sp|P35037|TRY3_ANOGA'


.. note:: Note the single quotes (') are mandatory as they are meant to protect the name of your sequence and prevent the Unix shell to interpret it like an instruction.

Once sequences are extracted or removed, some columns may remain containing only gaps, but it is possible to simply remove empty columns from the resulting dataset (command 3), and even extract specific blocks for the selected sequences either keeping the exact same name (command 4) or the name of the specific blocks extracted (command 5):

::

  Command 3: removing empty columns
  $$: t_coffee -other_pg seq_reformat -in sproteases_small.aln -action +extract_seq_list \
      'sp|P29786|TRY3_AEDAE' 'sp|P35037|TRY3_ANOGA' +rm_gap

  Command 4: keeping the initial name after extracting specific blocks and removing empty columns
  $$: t_coffee -other_pg seq_reformat -in sproteases_small.aln -action +keep_name \
      +extract_seq 'sp|P29786|TRY3_AEDAE' 20 200 'sp|P35037|TRY3_ANOGA' 10 150 +rm_gap

  Command 5: renaming sequences according to the extracted blocks and removing empty columns
  $$: t_coffee -other_pg seq_reformat -in sproteases_small.aln -action +extract_seq \
      'sp|P29786|TRY3_AEDAE' 20 200 'sp|P35037|TRY3_ANOGA' 10 150 +rm_gap 


.. hint:: The tag **+keep_name** must come BEFORE the tag **+extract_seq**.


Extracting the most informative sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Large datasets are problematic because they can be difficult to align and analyze, MSA programs tend to become very slow and inaccurate. In short, the best size for an MSA dataset would be between 20 to 40 sequences to have enough sequences to see the effect of evolution, but in the same time small enough so that you can visualize your alignment and recompute it as many times as needed. More important than its size, a good dataset have to be informative, when each sequence contains information the others do not have. The most informative sequences are the sequences that are as different as possible to one another, within your dataset. You can extract the most informative sequences using flag **+trim** followed by the number of sequences you wish to keep ("n" for a number and "N" for a pourcentage). The following commands will extract the 10 most informative sequences (command 1) or the 20% of most informative sequences (command 2):

::

  Command 1:
  $$: t_coffee -other_pg seq_reformat -in sproteases_large.fasta -action +trim _seq_n10 \
      -output fasta_seq
  Command 2:
  $$: t_coffee -other_pg seq_reformat -in sproteases_large.fasta -action +trim _seq_N20 \
      -output fasta_seq


.. hint:: The argument to trim include _seq_, it means your sequences are provided unaligned. If your sequences are already aligned, you do not need to provide this parameter. It is generaly more accurate to use unaligned sequences.

.. note:: For very large dataset, seq_reformat will compute the similarity matrix between your sequences once only. It will then store it in its cache to be reused any time you run on the same dataset. In short this means that it will take much longer to run the first time, but be much faster if you need to rerun it.


Extracting/Removing sequences with the % identity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
**Removing too identical sequences (redundant)**

Removing the most similar sequences is often what people have in mind when they talk about removing redundancy. You can do so using the **+trim** option. For instance, you can generate a dataset where no pair of sequences has more than 50% identity either from a dataset of unaligned sequences (command 1) or from any given alignment (command 2). If you start from unaligned sequences, the removal of redundancy can be slow. If your sequences have already been aligned using a fast method, you can take advantage of this by replacing the "_seq_" with "_aln_". Just run the following command lines to see the difference un runtime:

::

  Command 1: unaligned sequences
  $$: t_coffee -other_pg seq_reformat -in sproteases_large.fasta -action +trim _seq_%%50_

  Command 2: aligned sequences
  $$: t_coffee -other_pg seq_reformat -in sproteases_large.fasta -action +trim _aln_%%50_

.. note:: Using aligned sequences results in a fastest trimming, however, it also means that you rely on a more approximate estimation of sequence similarity.

**Removing too different sequences (outliers)**

Sequences that are too distantly related from the rest of the set (called outliers) may have very negative effects on the overall alignment; to prevent this, it is advisable not to use them. The next command line will lead to the removal of all the sequences where no pair of sequences has less than 30% average accuracy with all the other sequences in the dataset (the symbol "_O" stands for Outliers) and more than 80% identity: 

::

  $$: t_coffee -other_pg seq_reformat -in sproteases_large.fasta -action +trim _seq_%%80_O30

.. hint:: This particular option is quite powerful as it allows you to decide both inferior and superior tresholds for trimming your dataset based on pairwise identity score, and therefore you can dissect your dataset according to different ranges of identity values. Be careful not to remove too many sequences ;-)

**Forcing specific sequences to be kept**

Sometimes you want to trim based on identity while making sure specific/important sequences remain in your dataset. You can do so by providing a pattern ("_f" for field) : it will keep all the sequences whose name contains the given string ("_fNAME", "_fCOMMENT" or "_fSEQ", f standing for field). Here are some examples corresponding to the different protected fields while removing all sequences above 50% identity: 

::

  Keep all HUMAN sequences    
  $$: t_coffee -other_pg seq_reformat -in sproteases_large.fasta -action +trim \
      _seq_%%50_fNAME HUMAN

  Keep all sequences containing ".apiens"
  $$: t_coffee -other_pg seq_reformat -in sproteases_large.fasta -action +trim \
      _seq_%%50_fCOMMENT '.apiens'

  Keep all sequences containing residues
  $$: t_coffee -other_pg seq_reformat -in sproteases_large.fasta -action +trim \
      _seq_%%50_fSEQ '[MLV][RK]'

You can also specify the sequences you want to keep by giving another fasta file containing the name of these sequences via the flag **-in2**:

::

  $$:t_coffee -other_pg seq_reformat -in sproteases_large.fasta -in2 sproteases_small.fasta \
     -action +trim _seq_%%40


Chaining important sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In order to align two distantly related sequences, most multiple sequence alignment packages perform better when provided with many intermediate sequences that make it possible to 'bridge' your two sequences. The modifier **+chain** makes it possible to extract from a dataset a subset of intermediate sequences that chain the sequences you are interested in. For instance, let us consider the two sequences "sp|P21844|MCPT5_MOUSE" and "sp|P29786|TRY3_AEDAE" having 26% identity. This is high enough to make a case for a homology relationship between them, but this is too low to blindly trust any pairwise alignment. With the names of the two sequences written in the file sproteases_pair.fasta, run the following command:

::

  $$: t_coffee -other_pg seq_reformat -in sproteases_large.fasta -in2 sproteases_pair.fasta \
      -action +chain > sproteases_chain.fasta

This will generate a dataset of 21 sequences, with the following chain of similarity between your two sequences:

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

This is probably the best way to generate a high quality alignment of your two sequences when using a progressive method like ClustalW, T-Coffee, MUSCLE or MAFFT.



Modifying columns/blocks in your dataset
----------------------------------------
Removing gapped columns
^^^^^^^^^^^^^^^^^^^^^^^
You can also remove all the columns containing a given proportion of gaps; for instance the following command will delete all the residues occurring in a column that contains 50% or more gaps (use 1 to delete residues from columns having 1 gap or more):

::

  $$: t_coffee -other_pg seq_reformat -in sample_aln7.aln -action +rm_gap 50


Extracting specific columns 
^^^^^^^^^^^^^^^^^^^^^^^^^^^
Extracting portions of a dataset is something very frequently needed. You may need to extract all the sequences that contain the word human in their name, or you may want all the sequences containing a simple motif. We show you here how to do a couple of these things. To do this, you need an evaluation file that may have been generated with T-Coffee, either running a *de novo* alignment (command 1) or evaluating a preexisting alignment (command 2):

::

  Command 1:
  $$: t_coffee sample_seq1.fasta -output score_ascii, aln
  
  Command 2:
  $$: t_coffee -other_pg seq_reformat -in sample_seq1.aln -action +evaluate \
      blosum62mt -output score_ascii

This generates a score_ascii file that you can then use to filter out the bad bits in your alignment considering the individual score of each residue to trigger the filtering (command 3), or according to the whole column score by simply add the '+use_cons' flag (command 4). The commands 3 and 4 will keep only residues and columns having a score between 6 and 9:

::

  Command 3:
  $$: t_coffee -other_pg seq_reformat -in sample_seq1.aln -struc_in sample_seq1.score_ascii \
      -struc_in_f number_aln -action +keep '[6-9]'
  Command 4:
  $$: t_coffee -other_pg seq_reformat -in sample_seq1.aln -struc_in sample_seq1.score_ascii \
      -struc_in_f number_aln -action +use_cons +keep '[6-9]'


It is also possible to use a score_ascii file (as produced in the previous section) in order to extract high scoring portions of an alignment on-the-fly using the following command:

::

   $$: t_coffee -other_pg seq_reformat -in sample_aln1.aln -action +evaluate blosum62mt \
       +use_cons +keep '[5-9]'


.. warning:: Don't forget the simple quotes ('), it's mandatory !!!


Extracting entire blocks
^^^^^^^^^^^^^^^^^^^^^^^^
In case you want to extracting a specific block of your alignment for instance to remove poorly resolved regions, remove your alignments boudnaries or to extract specific domains, you can do so with the modified **+extract_block**. In this command line, the option **cons** (command 1) indicates that you are counting the positions according to the consensus of the alignment (i.e. the positions correspond to the columns # of the alignment). If you want to extract your block relatively to a specific sequence, you should replace cons with this sequence name (command 2).

::

  Command 1: extract block from MSA
  $$: t_coffee -other_pg seq_reformat -in sproteases_small.aln -action +extract_block \
      cons 150 200

  Command 2: extract_block relative to a give sequence of the MSA
  $$: t_coffee -other_pg seq_reformat -in sproteases_small.aln -action +extract_block \
      'sp|Q03238|GRAM_RAT' 10 200


.. tip:: It may be sometimes difficult to know where starts the blocks you are interested in except by counting manually the number of column. You can also make some tries by modifying the boundaries until you get the block you want and then redirect the result into the output file name of your choice. 


Concatenating blocks or MSAs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you have extracted several blocks generated using the previous command and you want to glue them together, you can use the **+cat_aln** modifier:

::

  $$: t_coffee -other_pg seq_reformat -in sproteases_small.aln -action +extract_block \
      cons 100 120 > block1.aln

  $$: t_coffee -other_pg seq_reformat -in sproteases_small.aln -action +extract_block \
      cons 150 200 > block2.aln

  $$: t_coffee -other_pg seq_reformat -in block1.aln -in2 block2.aln -action +cat_aln


.. note:: The alignments do not need to have the same number of sequences and the sequences do not need to come in the same order.


Manipulating DNA sequences
==========================
Translating DNA sequences into protein sequences
------------------------------------------------
If your sequences are DNA coding sequences, it is often safer and more accurate to align them as proteins (as protein sequences are more conserved than their corresponding DNA sequence). The seq_reformat option makes it easy for you to translate your sequences:

::

  $$: t_coffee -other_pg seq_reformat -in sproteases_small_dna.fasta -action \
      +translate -output fasta_seq


Back-translation with the *bona fide* DNA sequences
---------------------------------------------------
Once your sequences have been aligned, you may want to turn your protein alignment back into a DNA alignment, either to do phylogeny, or maybe in order to design PCR probes. To do so, use the following command:

::

  $$: t_coffee -other_pg seq_reformat -in sproteases_small_dna.fasta -in2 \
      sproteases_small.aln -action +thread_dna_on_prot_aln -output clustalw


Finding the *bona fide* sequences for the back-translation
----------------------------------------------------------
Use the online server `ProtoGen <http://tcoffee.vital-it.ch/apps/tcoffee/do:protogene>`_.


Manipulating RNA Sequences 
==========================
Producing a Stockholm output: adding predicted secondary structures
-------------------------------------------------------------------
Producing/Adding a consensus structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Given an RNA multiple sequence alignment, it is possible to compute (command 1) or add (command 2) the alifold (Vienna package) consensus secondary structure and output in in stockholm:

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
  $$: seq_reformat -in sample_rnaseq2.aln -in2 sample_rnaseq2.cons.stk -action +add_alifold \
      -output stockholm_aln


.. warning:: The alifold structure and the alignment MUST be compatible. The function makes no attempt to thread or align the structure, it merely stacks it below the MSA.

Analyzing a RNAalifold secondary structure prediction
-----------------------------------------------------
The following commands can either be applied on a Stockholm or a standard MSA. In the second case (standard MSA) the secondary structure will be automatically recomputed by alifold.


Analyzing matching columns
^^^^^^^^^^^^^^^^^^^^^^^^^^
The option **+alifold2cov_stat** will estimate the number of pairs of columns that are perfect Watson and Crick pairings, those that are neutral (including a GU) and those that include correlated mutations (command 1). The WCcomp are the compensated mutations maintaining WC base pairing. Other arguments can given, to display the list of paired positions and their status (compensated, Watson, etc...) use command 2:

::

  Command 1:
  $$: t_coffee -other_pg seq_reformat -in sample_rnaseq2.stk -action +alifold2analyze stat
  
  Command 2 (display list of options)
  $$: t_coffee -other_pg seq_reformat -in sample_rnaseq2.stk -action +alifold2analyze list


Visualizing compensatory mutations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The following command will output a color coded version of your alignment with matching columns indicated as follows:
I: incompatible pair (i.e. at least one pair is not WC)
N: pairs are Gus or WC
W: all pairs are Watson
c: compensatory mutations
C: WC compensatory mutations

::

  Standard alignment:
  $$: t_coffee -other_pg seq_reformat -in sample_rnaseq2.aln -action +alifold2analyze aln
  
  Color coded alignment:
  $$: t_coffee -other_pg seq_reformat -in sample_rnaseq2.aln -action +alifold2analyze color_html


.. warning:: Handling gapped columns: by default gapped column are ignored but they can be included by adding the tag **-usegap**.


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
The seq_reformat is NOT a phylogeny package, yet over the time it has accumulated a few functions that make it possible to compute simple phylogenetic trees, or similar types of clustering. Given a multiple sequence alignment, it is possible to compute either a UPGM or an NJ tree. The following commands use an identity matrix to compare your sequences and will output an unrooted NJ tree in newick format (command 1) or a rooted UPGMA tree (command 2):

::

  Command 1:
  $$:  t_coffee -other_pg seq_reformat -in <aln> -action +aln2tree -output newick

  Command 2:
  $$: t_coffee -other_pg seq_reformat -in <aln> -action +aln2tree _TMODE_upgma -output newick


If your data is not data sequence, but a matrix of 1 and Os (i.e. SAR matrix for instance), you can use a different matrix to compute the pairwise distances (command 3), and all these parameters can be concatenated (command 4):


::

  Command 3:
  $$: t_coffee -other_pg seq_reformat -in <aln> -action +aln2tree _MATRIX_sarmat -output newick

  Command 4:
  $$: t_coffee -other_pg seq_reformat -in <aln> -action +aln2tree _TMODE_upgma_MATRIX_sarmat \
       -output newick


.. warning:: Bootstrap facilities will also be added at some point...We recommend you to use `Phylip <http://evolution.genetics.washington.edu/phylip.html>`_ or any other specific phylogenetic software (PhyML, RAxML, MrBayes, etc...) if you need some serious phylogeny !


Comparing two phylogenetic trees
--------------------------------
A real interesting option is the ability to compare two trees (unrooted) returning some ofthe most common scores used for this including the famous Robinson-Foulds ;-)

::

  $$: t_coffee -other_pg seq_reformat -in sample_tree2.dnd -in2 sample_tree3.dnd -action \
      +tree_cmp -output newick

  #tree_cmp|T: 33 W: 20.00 L: 14.88 RF: 2 N: 9 S: 5
  #tree_cmp_def|T: ratio of identical nodes
  #tree_cmp_def|W: ratio of identical nodes weighted with the min Nseq below node
  #tree_cmp_def|L: average branch length similarity
  #tree_cmp_def|RF: Robinson and Foulds
  #tree_cmp_def|N: number of Nodes in T1 [unrooted]
  #tree_cmp_def|S: number of Sequences in T1


The output scores in more details:
- T: Fraction of the branches conserved between the two trees. This is obtained by considering the split induced by each branch and by checking whether that split is found in both trees
- W: Fraction of the branches conserved between the two trees. Each branch is weighted with MIN the minimum number of leaf on its left or right (Number leaf left, Number leaf Right)
- L: Fraction of branch length difference between the two considered trees.

The last line contains a tree where distances have been replaced by the number of leaf under the considered node:
- Positive values indicate a node common to both trees and correspond to MIN.
- Negative values indicate a node found in tree1 but not in tree2
- The higher this value, the deeper the node.


.. tip:: You can extract this tree for further usage by typing **cat outfile | grep -v 'tree_cmp'**


Scanning phylogenetic trees
---------------------------
It is possible to scan an alignment and locally measure the similarity between an estimated local tree and some reference tree provided from an external source (or computed on-the-fly) using the following command:

::

  $$ :t_coffee -other_pg seq_reformat -in <aln> -in2 <reftree> -action +tree_scan \
      _MODE_scan__W_10_ > ph_tree_scan.txt

For each position of the alignment, W*2 blocks of size 2*1+1 up to W*2+1 will be extracted, for each of these block a tree will be estimated and the similarity of that tree with the reference tree will be estimated with cmp_tree. For each position, the tree giving the best fit will be reported, along with the size of the block leading to that tree:

::

  P: <position> <block start> <blck_end> <block score> <block Length>

 
Pruning phylogenetic trees
--------------------------
Pruning removes leaves from an existing tree and recomputes distances so that no information is lost. Consider the files sample_tree2.dnd and the file sample_seq8.seq:

::

  sample_tree2.dnd
  (( A:0.50000, C:0.50000):0.00000,( D:0.00500, E:0.00500):0.99000, B:0.50000);
 
  sample_seq8.seq
  >A
  >B
  >C
  >D

.. note:: Sample_seq8 is merely a FASTA file where sequences can be omitted, but you can also leave them, at your entire convenience.

::

  $$: t_coffee -other_pg seq_reformat -in sample_tree2.dnd -in2 sample_seq8.seq -action \
      +tree_prune -output newick

  (( A:0.50000, C:0.50000):0.00000, B:0.50000, D:0.99500);


Manipulating structure files (PDB)
==================================
Extracting a structure
----------------------
There are many reasons why you may need a structure. T-Coffee contains a powerful utility named **extract_from_pdb** that makes it possible to fetch the PDB coordinates of a structure or its FASTA sequence without requiring a local installation. By default, the option **extract_from_pdb will** start looking for the structure in the current directory; it will then look it up locally (PDB_DIR) and eventually try to fetch it from the web (via a wget to www.rcsb.org). All these settings can be customized using environment variables (see next section). For instance if you want to fetch the chain E of the PDB structure 1PPG and/or its sequence in FASTA format, you can use:

::

  Fetch the structure:
  $$: t_coffee -other_pg extract_from_pdb -infile 1PPGE

  Fetch the correpsonding sequence:
  $$: t_coffee -other_pg extract_from_pdb -infile 1PPGE -fasta


Adapting extract_from_pdb to your own environment
-------------------------------------------------
If you have the PDB installed locally, simply set the variable PDB_DIR to the absolute location of the directory in which the PDB is installed. The PDB can either be installed in its divided form or in its full form. If the file you are looking for is neither in the current directory nor in the local PDB version, extract_from_pdb will try to fetch it from rcsb. If you do not want this to happen, you should either set the environment variable NO_REMOTE_PDB_DIR to 1 or use the **-no_remote_pdb_dir** flag:


::

  export NO_REMOTE_PDB_FILE=1
  
  t_coffee -other_pg extract_from_pdb -infile 1PPGE -fasta -no_remote_pdb_file


By default, T-Coffee also requires two important PDB files declared using the two following variables. These variables do not need to be set if the considered files are in the cache directory (default behavior):


::

  export PDB_ENTRY_TYPE_FILE=<location of the file pdb_entry_type.txt>
  (Found at: ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_entry_type.txt)
  
  export PDB_UNREALEASED_FILE=<location of the file unrealeased.xml>
  (Found at: http://www.rcsb.org/pdb/rest/getUnreleased)


.. warning:: Since the file ``unreleased.xml`` is not part of the pdb distribution, T-Coffee will make an attempt to obtain it even when using the NO_REMOTE_PDB_DIR=1 mode. You must therefore make sure that the file PDB_UNREALEASED_FILE is pointing to is read and write.


*************************************
Building Multiple Sequence Alignments
*************************************
General comments on alignments and aligners
===========================================
What is a good alignment?
-------------------------
This is a tricky question, a good answer would be  **"a good alignment is an alignment that makes it possible to do good biology"**. In practice, the alignment community has become used to measuring the accuracy of alignment methods using structures. Structures are relatively easy to align correctly, even when the sequences have diverged quite a lot. The most common usage is therefore to compare structure based alignments with their sequence based counterpart and to evaluate the accuracy of the method using these criterions. Unfortunately it is not easy to establish structure-based standards of truth. Several of these exist and they do not necessarily agree. To summarize, the situation is as roughly as follows:

  - **Above 40% identity** (within the reference dataset), all the reference collections agree with one another and all the established methods give roughly the same results. These alignments can be trusted blindly.

  - **Below 40% identity** (within the reference dataset), allthe reference collections stop agreeing and the methods do not give consistent results. In this area of similarity it is not necessarily easy to determine who is right and who is wrong, although most studies seem to indicate that consistency based methods (T-Coffee, ProbCons, MAFFT-slow or MSAProbs) have an edge over traditional methods.

When dealing with distantly related sequences, the only way to produce reliable alignments is to use structural information. T-Coffee provides many facilities to do so in a seamless fashion. Several important factors need to be taken into account when selecting an alignment method:

  - **The best methods are not always the best**. Given a difficult dataset, the best method is only more likely to deliver the best alignment, but there is no guaranty it will do so. It is very much like betting on the horse with the best odds.

  - **The difference in accuracy (as measured on reference datasets) between all the available methods is not incredibly high**. It is unclear whether this is an artifact caused by the use of 'easy' reference alignments, or whether this is a reality. The only thing that can change dramatically the accuracy of the alignment is the use of structural information.

  - **Keep in mind that these methods have only been evaluated by comparison with reference alignments (benchmarks)**. This is merely one criterion among many. In theory, these methods should be evaluated for their ability to produce alignments that lead to accurate trees, good profiles or good models. Unfortunately, these evaluation procedures do not yet exist.


The main methods and their scope
--------------------------------
There are many MSA packages around, the most common ones being ClustalW, MUSCLE, MAFFT, T-Coffee and ProbCons; amongst the latest ones, you can find phylogeny-aware aligners (PRANK and SATé) and modifed/improved consistency-based aligners (MSAProbs). You can almost forget about the other packages, as there is virtually nothing you could do with them that you will not be able to do with these packages. All these packages offer a complex trade-off between speed, accuracy and versatility.

ClustalW is everywhere you look...
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
ClustalW is still the most widely used Multiple Sequence Alignment package. Yet things are changing fast and different tests have consistently shown that ClustalW is neither the most accurate nor the fastest package around. This being said, ClustalW is everywhere and if your sequences are similar enough, it should deliver a fairly reasonable alignment.

MAFFT/MUSCLE to align big datasets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you have many sequences to align MUSCLE or MAFFT are the obvious choice. MAFFT is often described as the fastest and the most efficient. This is not entirely true, in its fast mode (FFT-NS-1), MAFFT is similar to MUSCLE and although it is fairly accurate, about 5 points less accurate than the consistency-based packages (ProbCons and T-Coffee). In its most accurate mode (L-INS-i) MAFFT uses local alignments and consistency, however, it becomes much more accurate but also slower, and more sensitive to the number of sequences. More recently, we have seen growing the number of **(ultra) large scale** aligners such as ClustalO, PASTA, UPP, and we hope soon the large scale version of T-Coffee (called MEGA-Coffee).

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
T-Coffee works by first assembling a library and then by turning this library into an alignment. The library is a list of potential pairs of residues. All of them are not compatible and the job of the algorithm is to make sure that as many possible constraints as possible find their way into the final alignment. Each library line is a constraint and the purpose is to assemble the alignment that accommodates the more all the constraints.

It is very much like building a high school schedule, where each teachers says something 'I need my Monday morning', 'I can't come on Thursday afternoon', and so on. In the end you want a schedule that makes everybody happy, if possible.The nice thing about the library is that it can be used as a media to combine as many methods as one wishes. It is just a matter of generating the right constraints with the right method and compile them into the library. ProbCons and MAFFT (L-INS-i) uses a similar algorithm, but with a Bayesian twist in the case of Probcons. In practice, however, ProbCons and T-Coffee give very similar results and have similar running time. MAFFT is significantly faster.

**Suited for**:
 - Profile reconstruction
 - Structure modeling
 - Function analysis
 - 3D prediction

Choosing the right package (without flipping a coin !)
------------------------------------------------------
Each available package has something to go for it, it is just a matter of knowing what you want to do !! T-Coffee is probably the most versatile, but it comes at a price and it is currently slower than many alternative packages. In the rest of this tutorial we give some hints on how to carry out each of these applications with T-Coffee.


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
A simple Multiple Sequence Alignment (default)
----------------------------------------------
T-Coffee default mode will simply compute a Multiple Sequence Alignment of the sequences you provided in input (command 1). It will display the final MSA on the screen and in several files according to the format you asked in command 2 (by default, the MSA is stored in a file .aln in ClustalW format). The headline of the alignment file contains important information such as the version of T-Coffee used, the CPU time, the overall consistency score (normalized to 100 or 1000 depending on the version of T-Coffee) and the total length of the MSA: it is quite practical to have a quick glance at the result. 

::

  Command 1: default MSA
  $$: t_coffee sproteases_small.fasta

  Command 2: default MSA, multiple output files
  $$: t_coffee sproteases_small.fasta -output=clustalw,fasta_aln,msf
  
Each time you run T-Coffee, 3 files are always generated:

 - the alignment:	``sproteases_small.aln``
 - the guide tree:	``sproteases_small.dnd``
 - the colored MSA:	``sproteases_small.html``

.. warning:: the guide tree is not a phylogenetic tree, it is used in the alignment process for clustering the sequences. 

.. tip:: you can visualize the colored html file with any browser/software you prefer. The display of the sequences should be aligned and formatted; if not, use another browser, it works quite well with Firefox, Safari, etc... If you need to do more sophisticated modifications on your MSA, we recommend to use `Jalview <http://www.jalview.org/>`_ which incorporate the T-Coffee color scheme.


Aligning multiple datasets/Combining multiple MSAs
--------------------------------------------------
If your sequences are spread across several datasets, you can give all the files you want (the limit is 200) via the flag **-seq**, and in any format you want. Just know that 1) if you give an alignment, the gaps will be reset and your alignment will only provide sequences, 2) sequences with the same name between two files are assumed to be the same sequence, 3) ff their sequences differ, they will be aligned and replaced by the consensus of that alignment (process known as sequence reconciliation). To align multiple datasets:

::

  $$: t_coffee -seq=sprotease1_small.fasta,sprotease2_small.aln -output=clustalw,fasta_aln,msf


You may also have a bunch of alignments (with the same sequences) that you have either precomputed, assembled manually or received from a colleague. You can also combine these alignments. For instance, let us imagine we generated 4 alignments with ClustalW using different gap penalties. To combine them into ONE single alignment, use the **-aln** flag. The final score indicates a high level of consistency (91%) between all these MSAs, meaning that the final MSA is probably correct.

::

  Your 4 different MSAs:
  clustalw -infile=sproteases_small.fasta -gapopen=0 -outfile=g0.aln
  clustalw -infile=sproteases_small.fasta -gapopen=-5 -outfile=g5.aln
  clustalw -infile=sproteases_small.fasta -gapopen=-10 -outfile=g10.aln
  clustalw -infile=sproteases_small.fasta -gapopen=-15 -outfile=g15.aln

  Combining multiple MSAs:
  $$: t_coffee sproteases_small.fasta -aln g0.aln g5.aln g10.aln g15.aln -output\
  clustalw html


Estimating the diversity in your alignment
------------------------------------------
It is easy to measure the level of diversity within your multiple sequence alignment. Will output all the pairwise identities, as well as the average level of identity between each sequence and the others. You can sort and grep in order to select the sequences you are interested in.


::

  $$: t_coffee -other_pg seq_reformat -in sample_aln1.aln -output sim


Comparing alternative alignments
--------------------------------
If you change the parameters, you will end up with alternative alignments. It can be interesting to compare them quantitatively. T-Coffee comes along with an alignment comparison module named **aln_compare**. You can use it to estimate the amount of difference between your two alignments either using the Sum-of-Pair score or the column score using the flag **-compare_mode** (sp or column). By default aln_compare returns the SoP score:

::

  $$: t_coffee -other_pg aln_compare -al1 b80.aln -al2 b30.aln -compare_mode sp


This comparison will return the following result:

::

  *****************************************************
  seq1       seq2          Sim   [ALL]           Tot  
  b80           19         33.5    89.5 [100.0]   [ 8958]

The interpretation of this output is as follow: b80 is the reference MSA, it contains 19 sequences with an average identity of 33.5%, and is 89.5% identical to the second MSA b30.aln (8958 pairs to be precise). Of course, this does not tell you where are the good bits, but you can get this information for instance residues that have lost more than 50% of their pairing partner between the two alignments are now in lower case.

:: 

  $$: t_coffee -other_pg aln_compare -al1 b30.aln -al2 p350.aln -output_aln \
      -output_aln_threshold 50

  $$: t_coffee -other_pg aln_compare -al1 b30.aln -al2 p350.aln -output_aln \
      -output_aln_threshold 50 -output_aln_modif x


.. tip:: This option is particularly interesting if you are modifying the default parameters of T-Coffee and want to monitor the effects of your modifications. 


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

 - By default, these alignments are computed using a Blosum62 matrix, but you can use any matrix you fancy instead, including: pam120mt, pam160mt, pam250mt, pam350mt, blosum30mt, blosum40mt, blosum45mt, blosum50mt, blosum55mt, blosum62mt, blosum80mt, or even user-provided matrices in the BLAST format, as described in the technical manual.

 - PAM matrices: These matrices are allegedly less accurate than the blosum. The index is correlated to the evolutionary distances. You should therefore use the pam350mt to align very distantly related sequences.

 - Blosum matrices: These matrices are allegedly the most accurate. The index is correlated to the maximum percent identity within the sequences used to estimate the matrix. You should therefore use the Blosum30mt to align very distantly related sequences. Blosum matrices are biased toward protein core regions. This may explain why these matrices tend to give better alignments, since by design they can capture the most evolutionary resilient signal contained in proteins.

Unless you have some structural information available, the only way to tell whether your alignment has improved or not is to look at the score. For instance, if you compute the two following alignments:

::

  $$: t_coffee sproteases_small.fasta -matrix=blosum30mt -outfile=b30.aln
  $$: t_coffee sproteases_small.fasta -matrix=blosum80mt -outfile=b80.aln

You will get two alignments that have roughly the same score but are different. You can still use these two alternative alignments by comparing them to identify regions that have been aligned identically by the two matrices. These regions are usually more trustworthy.

Changing gap penalties
^^^^^^^^^^^^^^^^^^^^^^
.. important:: Gap penalties are the core of the matter when it comes to multiple sequence alignments. An interesting feature of T-Coffee is that it does not really need such penalties when assembling the MSA, because in theory the penalties have already been applied when computing the library. This is the theory, as in practice penalties can help improve the quality of the alignment.

The penalties can be changed via the flags **-gapopen** for the gap opening penalty and via **-gapext** for the gap extension penalty. The range for gapopen are [-500,-5000], the range for the extension should rather be [-1,-10]. These values do not refer to a substitution matrix, but rather to the values range of the consistency estimation (i.e. ratio) normalized to 10000 for a maximum consistency. The default values are **-gapopen=-50, -gapext=0**. The reasons for these very low values are that they are meant to be cosmetic only, since a trademark of T-Coffee (inherited from Dialign) is not to need explicit penalties. Yet, we know for a fact that alignments with higher gap penalties often look nicer (for publications) and are sometimes more accurate. For instance, you can try:

::

  $$: t_coffee sproteases_small.fasta -gapopen -100 -gapext -5

This gap penalty is only applied at the alignment level (i.e. after the library was computed). If you want to change the gap penalties of the methods used to build the library, you will need to go deeper...Two methods are used by default to build the library (command 1). One does global pairwise alignments and is named slow_pair, the other is named lalign_id_pair and produces local alignments. These methods are specified via the **-method** flag. Usually you do not need to write it because it is the default, but if you want to change the default parameters of the constituting methods (command 2), you will need to do so explicitly (the default parameters are for lalign_id_pair **GOP=-10, GEP=-4, MATRIX=blosum50mt** and for slow_pair **GOP=-10, GEP=-1 and MATRIX=blosum62mt**. Using the command 2, the library is now computed using the Blosum62mt with lalign, rather than the Blosum50mt; the good news is that when using this matrix, the score of our alignment increases from 48 to 50. We assume this new alignment is therefore more accurate than the previous one.

::

  Command 1: default T-Coffee
  $$: t_coffee sproteases_small.fasta -method=lalign_id_pair,slow_pair

  Command 2: modifiying the parameters
  $$: t_coffee sproteases_small.fasta -method lalign_id_pair@EP@MATRIX@blosum62mt, \
      slow_pair -outfile sproteases_small.b62_aln


.. warning:: It only makes sense to compare the consistency score of alternative alignments when these alignments have been computed using the same methods (lalign_id_pair and slow_pair for instance).


Aligning (very) large datasets
==============================
Aligning (very) large datasets with MUSCLE
------------------------------------------
T-Coffee is not a good choice if you are dealing with very large datasets, use MAFFT or MUSCLE. To align a large dataset with MUSCLE, try:

::

  Default mode
  muscle -infile sproteases_large.fasta > sproteases_large.muscle
  
  Fast mode (less accurate)
  muscle -in sproteases_large.fasta -maxiters 1 -diags -sv -distance1 kbit20_3 \
  > sproteases_large.muscle


Aligning (very) large datasets with MAFFT
-----------------------------------------
The fastest mode with MAFFT can be achieved using: **mafft --retree 2 input > output**.


Aligning (very) large alignments with T-Coffee
----------------------------------------------
T-Coffee is not very well gifted for aligning large datasets (for now), but you can give it a try using a special option that generates approximate fast alignments (command 1). These MSAs should roughly have the same accuracy as ClustalW, and are quite acceptable for sequences more than 40% identical. This mode works by only considering the best diagonals between two sequences, and by default all the diagonals with substitution score >0 are considered, but you can lower this to reduce the running time (command 2). That will only consider the top 10 diagonals. This will be very useful if you have long and very similar sequences to align (DNA for instance).

::

  Command 1:
  $$: t_coffee sproteases_large.fasta -mode quickaln

  Command 2:
   $$: t_coffee sample_seq1.fasta -mode quickaln -ndiag=10

Another alternative to align large datasets is a special mode of T-Coffee, fm-Coffee (command 3), derived from M-Coffee (see next section) and designed to be fast and able to handle large datasets (it is used for example in Ensembl). To do so, T-Coffee used three different fast aligners: MAFFT, MUSCLE and Kalign. 

::

  Command 3:
  $$: t_coffee sproteases_large.fasta -mode fmcoffee

.. tip:: Once you have your large MSA, you can always shrink/trim them using reformatting options (see previous section) for instance by extraction the most informative sequences or by defining %identity cut-off.

.. note:: In the last 10 years, a special effort have been made to improve large scale alignment leading to the development of few new methods among which Clustal Omega, PASTA, UPP and we hope soon a MEGA-Coffee aligner. These methods are not incorporated in T-Coffee so if your datasets are really large (>5000 sequences) don't hesitate to use these methods instead.


Using many methods at once
==========================
One of the most common situation when building MSAs is to have several alignments produced by different alternative methods, and not knowing which one to choose. In this section, we show you how to use M-Coffee to combine many alignments into one single alignment, or how you can specify only the methods you want. M-Coffee is not always the best method, but extensive benchmarks on BaliBase, Prefab and Homstrad have shown that it delivers the best alignment 2 times out of 3. If you do not want to use the methods provided by M-Coffee, you can also combine precomputed alignments.


Using all the methods at the same time: M-Coffee
------------------------------------------------
In M-Coffee, M stands for Meta. To use M-Coffee, you will need several packages to be installed (see **T-Coffee Installation** and section **Integrating External Methods in T-Coffee**). If you did a default installation, all the software you need should be there. M-Coffee is a special mode of T-Coffee that you can call using the flag **-mode mcoffee**. It will align your sequence using 8 different aligners: ClustalW, POA, MUSCLE, ProbCons, MAFFT, Dialing-T, PCMA and T-Coffee:

::

  $$: t_coffee sproteases_small.fasta -mode mcoffee -output clustalw, html


When this is done, all the alignments will be combined into one. If you open the file sproteases_small.html with your favorite web browser, you will see a colored version of your alignment: the alignment is colored according to the consistency of the different methods used. Regions in red have a high consistency, so all the methods agree and you can expect them to be fairly accurate. Regions in green/blue have the lowest consistency, meaning that all the methods deliver different alignment in these regions and you should not trust them. Overall this alignment has a score of 80, which means that it is 80% consistent with the entire collection. This is a fairly high index, which means you can probably trust your alignment (at least where it is red).

Using selected methods to compute your MSA
-------------------------------------------
Using the 8 Methods of M-Coffee can sometimes be a bit heavy, if you only want to use a subset of your favorite methods, you should know that each of these methods is available via the **-method flag**. You can make all the combination you want !!! For instance, to combine MAFFT, MUSCLE, T-Coffee and ProbCons, you can use:

::

  $$: t_coffee sproteases_small.fasta -method=t_coffee_msa,mafft_msa,probcons_msa, \
      muscle_msa -output=html


Aligning profiles 
=================
Sometimes, it is better to prealign a subset of your sequences, and then to use this small alignment as a master for adding sequences (sequence to profile alignment) or even to align several profiles together if your protein family contains distantly related groups. T-Coffee contains most of the facilities available in ClustalW to deal with profiles, and the strategy we outline here can be used to deal with large datasets

Aligning one sequence to a profile
----------------------------------
Assuming you have a multiple alignment (sproteases_small.aln) here is a simple strategy to align one sequence to your profile:

::

  $$: t_coffee sproteases_oneseq.fasta -profile sproteases_small.aln


Aligning many sequences to a profile
------------------------------------
You can align as many sequences as you wish to your profile. Likewise, you can have as many profiles as you want. For instance, the following:


::

  $$: t_coffee sequences.fasta -profile=prf1.aln,prf2.aln,prf3.aln -outfile=comb\
 ined_profiles.aln



Will make a multiple alignment of 3 profiles and 5 sequences. You can mix sequences and profiles in any proportion you like. You can also use all the methods you want although you should be aware that when using external methods (see the external method section in this tutorial), the profile is replaced with its consensus sequence, which will not be quite as accurate. Methods supporting full profile information are: lalign_id_pair, slow_pair and proba_pair, clustalw_pair and clustalw_msa. All the other methods (internal or external) treat the profile as a consensus (less accurate).


Computing very accurate (but slow) alignments with PSI-Coffee
-------------------------------------------------------------
PSI-Coffee is currently the most accurate mode of T-Coffee and also the slowest; its principle is rather simple: it associates every sequence with a profile of homologous sequences gathered using BLAST on a sequence database (nr by default). PSI-Coffee then uses the profiles instead of the initial sequences to makes a multiple profile alignment. In a last step, your profiles are replaced by their initial query sequence from your initial dataset and returns a MSA of your sequences.

::

  $$: t_coffee sproteases_small.fasta -mode psicoffee


Aligning other types of sequences
=================================
Splicing variants
-----------------
Splicing variants are especially challenging for most MSA programs. This is because the splicing variants need very long gaps to be inserted, while most programs attempt to match as many symbols as possible.


Standard programs like ClustalW or Muscle are not good at dealing with this situation and in our experience, the only programs that can do something with splice variants are those using local information like some flavors of Mafft and T-Coffee .


For instance, if you try muscle on the following dataset:


::

  muscle -in sv.fasta -clw



You will quickly realise that your alignment is not very good and does not show where the alternative splicing coocurs. On the other hand, if you use T-Coffee, things become much clearer


::

  $$: t_coffee sv.fasta



The reason why T-Coffee does better than other packages is mostly because it uses local information (lalign_id_pair) and is therefore less sensitive to long gaps. If the default mode does not work for your dataset, you can try to be a bit more aggressive and only use local information to compute your library:


::

  $$: t_coffee sv.fasta -method lalign_id_pair



Of course, the most distantly related your sequences, the harder the alignment of splicing variants


Aligning DNA sequences
----------------------
Multiple Sequence Alignment methods are not at their best when aligning DNA sequences. Whenever you can, try using a local multiple sequence alignment package like the Gibbs sampler. Yet if you believe your DNA sequence are homologous over their entire length, you can use T-Coffee.


In theory, the program automatically recognizes DNA sequences and uses appropriate methods, yet adding the -type=dna flag cannot do any harm...


::

  $$: t_coffee sample_dnaseq1.fasta -type=dna



The type declaration (or its automatic detection) triggers the use of the appropriate substitution matrix in most of the methods. In practice, any time it encounters dna, the program will try to use '4dna' version of the requested methods. These methods have lower penalties and are better suited for dealing with nucleic acid sequences.


However, if you would rather use your own matrix, use:


::

  $$: t_coffee sample_dnaseq1.fasta -in Mlalign_id_pair4dna@EP@MATRIX@idmat



Where you should replace idmat with your own matrix, in BLAST format (see the format section of the Reference Manual).


Aligning RNA sequences
----------------------
RNA sequences are very important and almost every-where these days. The main property of RNA sequences is to have a secondary structure that can be used to guide the alignment. While the default T-Coffee has no special RNA alignment method incorporated in, smart people have thought about this. If you are interested in RNA, check: http://www.bio.inf.uni-jena.de/Software/MARNA/.


Noisy coding DNA sequences...
-----------------------------
When dealing with coding DNA, the right thing to do is to translate your DNA sequence and thread the DNA onto the protein alignment if you really need some DNA. However, sometimes, your cDNA may not be so clean that you can easily translate it (frameshifts and so on). Whenever this happens, try (no warranty) the following special method.


The test case in three_dna_seq.fasta contains the DNA sequences of three proteases with a couple of frameshifts here and there. If you make a regular alignment of these sequences


::

  $$: t_coffee three_cdna.fasta



You can immediately see that many gaps have sizes that are not multiple of 3 (codon size). Most of the information is lost. On the other hand, when using an appropriate alignment method that takes into account all the frames at the same time, we get something much more meaningful:


::

  $$: t_coffee three_cdna.fasta -method cdna_fast_pair



And most importantly, the frameshifts end up at the right place. You can even recover the corrected protein sequence using a special mode of seq_reformat:


::

  $$: t_coffee -other_pg seq_reformat -in three_cdna.aln -action +clean_cdna +tr\
 anslate



+clean cdna is a small HMM that loops through each sequence and select the frame in order to maximize the similarity within the alignment.


**************************************************************************************
Combining Protein 2D and 3D Structural Information In Your Multiple Sequence Alignment
**************************************************************************************
Using structural information when aligning sequences is very useful. The reason is that structures diverge slower than sequences. As a consequence, one may still find a discernable homology between two sequences that have been diverging for so long that their sequences have evolved beyond recognition. Yet, when assembling the correct structure based MSA, you will realize that these sequences contain key conserved residues that a simple alignment procedure was unable to reveal. We show you in this section how to make the best of T-Coffee tools to incorporate structural information in your alignment.


If you are in a hurry: Expresso
===============================
What is Expresso?
-----------------
Expresso is the latest T-Coffee mode. It is not yet available for local installation, but you can run it from the www.tcoffee.org server. The principle of Expresso is simple: the server runs a BLAST between every sequence in your query against the PDB database. If it finds a structure similar enough to a sequence in your dataset (>60% identity), it will use that structure as a template for your sequence.


Template files look something like:


::

  >sp|P08246|ELNE_HUMAN _P_ 1PPGE

  >sp|P20160|CAP7_HUMAN _P_ 1AE5

  >sp|P00757|KLKB4_MOUSE _P_ 1SGFX

  >sp|Q6H321|KLK2_HORSE _P_ 1GVZA

  >sp|P00773|ELA1_RAT _P_ 2D26C

  >sp|Q00871|CTRB1_PENVA _P_ 1AZZB

  >sp|P21844|MCPT5_MOUSE _P_ 1NN6A

  >sp|O35205|GRAK_MOUSE _P_ 1MZDA

  >sp|P07338|CTRB1_RAT _P_ 2CGAB

  >sp|P80015|CAP7_PIG _P_ 1FY3A

  >sp|P03953|CFAD_MOUSE _P_ 1FDPD

  >sp|Q7YRZ7|GRAA_BOVIN _P_ 1OP8F

  >sp|Q06606|GRZ2_RAT _P_ 1EUFA

  >sp|P08884|GRAE_MOUSE _P_ 1FI8B



In a template file, _P_ indicates that the template is of type structure (P for PDB). Template files can be generated manually or automatically by the Expresso server. Whenever possible t_coffee will then align your sequences using the structural information contained in the templates. If it encounters enough structures (as shown here) it will produce a genuine structure based sequence alignment.


Using Expresso
--------------
::

  $$: t_coffee three_pdb_two_seq.fasta -method sap_pair,slow_pair -template_file\
  PDB



Using secondary structure predictions:
======================================
T-Coffee can be used to predict secondary structures and transmembrane domains. For secondary structure predictions, the current implementation is only able to run GOR on either single sequences or on a bunch of homologues found by BLAST.


Single sequence prediction
--------------------------
To make a secondary structure prediction with GOR, run the following. In this command line SSP is a hard coded mode. It prompts the computation of predicted secondary structures.


::

  t_coffee sample_aln.fasta -template_file SSP



The predictions are then displayed in the files:


::

  #### File Type= Template Protein Secondary Structure Format= fasta_seq Name= hmgb_chite.ssp

  #### File Type= Template Protein Secondary Structure Format= fasta_seq Name= hmgl_trybr.ssp

  #### File Type= Template Protein Secondary Structure Format= fasta_seq Name= hmgl_trybr3.ssp

  #### File Type= Template Protein Secondary Structure Format= fasta_seq Name= hmgl_wheat.ssp

  #### File Type= Template Protein Secondary Structure Format= fasta_seq Name= hmgl_wheat2.ssp

  #### File Type= Template Protein Secondary Structure Format= fasta_seq Name= hmgt_mouse.ssp



Transmembrane structures can be carried out with:


::

  $$: t_coffee sample_aln.fasta -template_file TM



Multiple sequences predictions
-----------------------------
Used this way, the method will produce for each sequence a secondary prediction file. GOR is a single sequence with a relatively low accuracy. It is possible to increase the accuracy by coupling BLAST and GOR, this can be achieved with the following command:


::

  $$: t_coffee sample_aln.fasta -template_file PSISSP



When doing so, the predictions for each sequence are obtained by averaging the GOR predictions on every homologue as reported by a BLAST against NR. By default the BLAST is done remotely at the NCBI using the blastpgp web service of the EBI.


A similar output can be obtained for Transmembrane segment predictions:


::

  $$: t_coffee sample_aln.fasta -template_file PSITM



Incorporation of the prediction in the alignment
------------------------------------------------
It is possible to use the secondary prediction in order to reward the alignment of similar elements


::

  $$: t_coffee sample_aln.fasta -template_file PSISSP -method_evaluate_mode ssp -met\
 hod lalign_id_pair slow_pair



Likewise, it is possible to use this information with trans-membrane domains


::

  $$: t_coffee sample_aln.fasta -template_file PSITM -method_evaluate_mode tm -metho\
 d lalign_id_pair slow_pair



The overall effect is very crude and amounts to over-weighting by 30% the score obtained when matching two residues in a similar secondary structure state. The net consequence is that residues in similar predicted states tend to be aligned more easily.


Using other secondary structure predictions
-------------------------------------------
If you have your own predictions, you can use them. All you need is to produce a template file where the file containing the secondary structure prediction is declared along with the sequence:


::

  >hmgl_wheat _E_ hmgl_wheat.ssp

  >hmgb_chite _E_ hmgb_chite.ssp

  >hmgl_trybr3 _E_ hmgl_trybr3.ssp

  >hmgl_wheat2 _E_ hmgl_wheat2.ssp

  >hmgt_mouse _E_ hmgt_mouse.ssp

  >hmgl_trybr _E_ hmgl_trybr.ssp



where each template looks like this:


::

  >hmgl_wheat

  CCCCCCCCCCCCHHHHHHHCCCCCCCCCHHHHHHHHHHHHHHHCCCCHHHHHHHHHHHHHHHCE



You can then run T-Coffee using your own template file


::

  $$: t_coffee sample_aln.fasta -template_file <template_file> -method_evaluate_mode\
  ssp -method lalign_id_pair slow_pair



Output of the prediction
------------------------

You can output a color coded version of your alignment using the predicted structures


::

  $$: t_coffee sample_aln.fasta -template_file PSISSP -output sec_html



A similar result can be obtained with trans-membrane regions:


::

  $$: t_coffee sample_aln.fasta -template_file PSITM -output tm_html



Aligning sequences and structures
=================================
Mixing sequences and structures
-------------------------------
Gather your sequences in the same file. Name your structures according to their PDB identifier. The file three_pdb_two_seq.fasta contains five sequences, three are the sequences of PDB structures and two are regular sequences.


What you want to do is to build a T-Coffee library where sequences with a known structures are aligned with a structure alignment program (like sap) while the other sequences are aligned using regular T-Coffee methods. You can achieve this with the following command:


::

  $$: t_coffee three_pdb_two_seq.fasta -method sap_pair,slow_pair -template_file\
  PDB



The option -template_file is here to tell the program how to find the PDB. In that case. EXPRESSO means that a remote BLAST (the EBI BLAST) will be used to identify the best targets. If your sequences are already named according to their PDB name, you can use:


::

  $$: t_coffee three_pdb_two_seq.fasta -method sap_pair,slow_pair -template_file\
  _SELF_P_



_SELF_ means that the PDB identifier is the name of the sequences, while _P_ is an indication that the template is indeed a PDB. These indications are necessary for T-Coffee to fetch the relevant structures.


The good news is that you do not need to have PDB installed locally as T-Coffee will automatically fetch the structures directly from RCSB (the home of PDB). Of course, if your dataset only contains structures, your alignment becomes a structural alignment.


If you have a fugue license, you can also add the fugue method to your run. Fugue will align the structures with sequences whose structure is unknown (this is called threading).


::

  $$: t_coffee three_pdb_two_seq.fasta -method sap_pair,slow_pair,fugue_pair -te\
 mplate_file _SELF_P_



This can be written more concisely, using one of T-Coffee special_modes:


::

   $$: t_coffee three_pdb_two_seq.fasta -mode 3dcoffee



or


::

   $$: t_coffee three_pdb_two_seq.fasta -mode expresso



Using sequences only
--------------------
What often happens is that you have already built a dataset with sequences that are very similar to PDB sequences but not exactly identical. It may even be the case that the real sequence and the PDB one do not match exactly because of some genetic engineering on the structure. In this case, you have no structure whose sequence is exactly similar to the sequences in your dataset. All you need to do is to declare the equivalence sequences/structures and run T-Coffee, just like Expresso does.


The first step is to fill up a template file that contains an explicit declaration of the structures corresponding to your sequences. The format is very simple and fasta-like. You can use the file: sproteases_small.template_file


::

  >sp|P08246|ELNE_HUMAN _P_ 1PPGE
  >sp|P20160|CAP7_HUMAN _P_ 1AE5
  >sp|P00757|KLKB4_MOUSE _P_ 1SGFX
  >sp|Q6H321|KLK2_HORSE _P_ 1GVZA



In this file, the first line is telling us that sequence sp|P08246|ELNE_HUMAN is associated with the structural template 1PPGE. The sequence and the structure do not need to be identical although we recommend using structural templates more than 60% identical with your actual sequences (i.e. similar enough so that they generate a non ambiguous alignment). If your template file is ready, all you need to do is run the following command.


::

  $$: t_coffee sproteases_small.fasta -method slow_pair, lalign_id_pair, sap_pai\
 r -template_file sproteases_small.template_file



When you run this once, T-Coffee goes and fetches the structures. It will then align them using sap. It takes a lot of time to fetch structures, and it takes even more time to align them with sap. This is why T-Coffee saves these important intermediate results in a special location called the cache. By default, your cache is in ~/.t_coffee/cache, it is a good idea to empty it from time to time...


Aligning profile using structural information
=============================================
If you have two profiles to align, an ideal situation is when your profiles each contain one or more structures. These structures will guide the alignment of the profiles, even if they contain very distally related sequences. We have prepared two such profiles (prf1_pdb1.aln, prf2_pdb2.aln). You have two choices here. All you need is a template file that declares which sequences have a known structure. If you only want to align sequences, you can try:


::

  $$: t_coffee -profile=profile1_pdb1.aln, profile2_pdb2.aln -method sap_pair -p\
 rofile_template_file two_profiles.template_file



****************************
How Good Is Your Alignment ?
****************************
There are three strategies for evaluating your alignment. Structure is a killer so if you have at least two structures available for your protein family, you are in an ideal situation and you can use the iRMSD. If you only have one structure available, we developped STRIKE to compare alternative alignment. If you don't have any structure, you are left with the option of using sequence based methods like the CORE index and the TCS. These do pretty well in the CORE regions, but can be limited in the loops. Another killer, less often at hand, is the use of functional information. If you know some residues MUST be aligned because they are functionally related, you can easily set up an evaluation procedure using T-Coffee.


Transitive Consistency Score - TCS 
==================================
TCS is an alignment evaluation score that makes it possible to identify the most correct positions in an MSA. 
It has been shown that these positions are the most likely to be structuraly correct and also the most informative when estimating phylogenetic trees. The TCS evaluation and filtering procedure is implemented in the T-Coffee package and can be used to evaluate and filter any third party MSA (including T-Coffee MSA of course!)


Evaluate an existing MSA 
------------------------

:: 

  $$: t_coffee -infile prot.aln -evaluate -output score_ascii, aln, score_html

Output files: 

* ``prot.score_ascii``  displays the score of the MSA, the sequences and the residues. This file can be used to further filter your MSA with seq_reformat. 
* ``prot.score_html`` displays a colored version score of the MSA, the sequences and the residues. 

.. warning:: The color code in the score_html indicates the agreement between the library and the considered alignment. It is important to understand that this score does not only depend on the input MSA, but it also depends on the library.

.. tip:: The TCS is most informative when used to identify low-scoring portions within an MSA. It is also worth noting that the TCS is not informative when aligning less than five sequences.
  
Filter unreliable MSA positions
-------------------------------

:: 

  $$: t_coffee -infile prot.aln -evaluate -output tcs_residue_filter3, tcs_column_filter3 \
      , tcs_residue_lower4

Output file: 

* ``prot.tcs_residue_filter3``  All residues with a TCS score lower than 3 are filtered out 
* ``prot.tcs_column_filter3``   All columns with a TCS score lower than 3 are filtered out 
* ``prot.tcs_residue_lower4``   All residues with a TCS score lower than 3 are lower cased
  
Note that all these output functions are also compatible with the default T-Coffee when computing an alignment:

  $$: t_coffee -seq prot.fa -output tcs_residue_filter3, tcs_column_filter3, tcs_residue_lower4

or with **seq_reformat** using a T-Coffee .score_ascii file:: 

  $$: t_coffee -other_pg seq_reformat -in prot.aln -struc_in prot.score_ascii -struc_in_f \
      number_aln -output tcs_residue_filter3
  

Weight MSA for improved trees
-----------------------------

:: 

  $$: t_coffee -infile prot.aln -evaluate -output tcs_weighted, tcs_replicate_100

Output files: 

* ``prot.tcs_weighted``       All columns are duplicated according to their TCS score 
* ``prot.tcs_replicate_100``  Contains 100 replicates in phylip format with each column drawn with a probability corresponding to its TCS score 


Note that all these output functions are also compatible with the default T-Coffee when computing an alignment:

  $$: t_coffee -seq prot.fa -output tcs_weighted, tcs_replicate_100

or with **seq_reformat** using a T-Coffee .score_ascii file::

  $$: t_coffee -other_pg seq_reformat -in prot.aln -struc_in prot.score_ascii -struc_in_f \
      number_aln -output tcs_weighted


Working with coding DNA
-----------------------

When working with DNA, it is advisable to first align the sequences at the protein level and later thread back the DNA onto your aligned proteins. The filtering must be done in two steps, as shown below. Note that your DNA and protein sequences must have the same name::

  $$: t_coffee -infile prot.aln -evaluate -output score_ascii

This first step produces the TCS evaluation file ``prot.score_ascii``:
 
  $$: t_coffee -other_pg seq_reformat -in prot.aln -in2 dna.fa -struc_in prot.score_ascii \
      -struc_in_f number_aln -output tcs_replicate_100 -out dna.replicates
  
The dna.replicates option: 100 DNA replicates with positions selected according to their aminoacid TCS score::

  $$: t_coffee -other_pg seq_reformat -in prot.aln -in2 dna.fa -struc_in prot.score_\
  ascii -struc_in_f number_aln -output tcs_column_filter5 -out dna.filter  

The dna.filtered option: DNA positions filtered according to their TCS column score


Using different libraries
-------------------------

It is possible to change the way TCS reliability is estimated. This can be done by building different T-Coffee libraries. The proba_pair is the default aligner of T-Coffee that runs a pair-HMM to populate the library with residue pairs having the best posterior probabilities. The following instructions will do this: 

  $$: t_coffee -infile prot.aln -evaluate -method proba_pair -output score_ascii, aln, score_html

This mode runs a series of fast multiple aligners; it is very fast and used by `ENSEMBL Compara <http://www.ensembl.org/info/genome/compara/index.html>`_:

  $$: t_coffee -infile prot.aln -evaluate -method mafft_msa,kalign_msa,muscle_msa -output \
      score_ascii, aln, score_html

This mode runs the orginal default T-Coffee that was combining local and global alignments::

  $$: t_coffee -infile prot.aln -evaluate -method clustalw_pair,lalign_id_pair -output \
      score_ascii, aln, score_html


Summary of the output flags
---------------------------

============================  ================
Flags        		      Description
============================  ================
-output=score_ascii	      outputs a TCS evaluation file
-output=score_html	      contains ascii format in html format
-output=score_pdf	      will transfer score_html into pdf format
-output=sp_ascii	      format reporting the TCS score of every aligned pair in the target MSA
-output=tcs_residue_filter_N  removes all residues with a TCS score lower than `N`
-output=tcs_columns_filter_N  removes all columns with a TCS score lower than `N`
-output=tcs_weighted	      phylip format with duplicated columns according to their TCS score
-output=tcs_replicate_N	      generates `N` replicates with columns drawn according to their TCS score
============================  ================


.. Note:: Please cite

	* `Chang, J.-M., Tommaso, P. & Notredame, C. TCS: A New Multiple Sequence Alignment Reliability Measure to Estimate Alignment Accuracy and Improve Phylogenetic Tree Reconstruction. Molecular biology and evolution 31, 1625–37 (2014). <http://www.ncbi.nlm.nih.gov/pubmed/24694831>`_
	* `Chang, J.-M., Di Tommaso, P., Lefort, V., Gascuel, O. & Notredame, C. TCS: a web server for multiple sequence alignment evaluation and phylogenetic reconstruction. Nucleic Acids Res 43(W1):W3-6 (2015). <http://www.ncbi.nlm.nih.gov/pubmed/25855806>`_.


Evaluating alignments with The CORE index
=========================================

.. note:: This procedure is now being deprecated and is meant to be replaced by the TCS procedure outlined in the above section

Computing the local CORE index
------------------------------
The CORE index is an estimation of the consistency between your alignment and the computed library. The higher the consistency, the better the alignment. The score reported with every T-Coffee alignment is the concistency score. However, if you want to go further and estimate the local concistency (known as the CORE index). Simply request one extra output:

::

  $$: t_coffee sproteases_small.fasta -output=clustalw,html


The format html leads to the output of a file named sproteases_small.html. Open this file. It contains a colorized version of your alignment. In this colorized version, positions that have no concistency with the library are in blue, a little in green, better positions in yellow, then orange, then red. You can expect yellow positions to be entirely correct.


Computing the CORE index of any alignment
-------------------------------------------

You can evaluate any existing alignment with the CORE index. All you need to do is provide that alignment with the -infile flag and specify that you want to evaluate it:

::

  $$: t_coffee -infile=sproteases_small.g10.cw_aln -output=html -score


For more information on filtering/trimming an alignment using the CORE index score, refer to the subsection **Preparing Your Data: Reformatting, Trimming an more.../Modifying the data itself**.


Evaluating an alignment using structural information: APDB and iRMSD
====================================================================
What is the iRMSD?
------------------
APDB and the iRMSD are two closely related measures meant to evaluate the accuracy of a sequence alignment without using a structure based reference alignment. The iRMSD is a follow up of the APDB measure and we now recommend using the iRMSD rather than APDB.


Although it may seem that the iRMSD was an attempt to get free iPODs from Apple, it is not (or at least we never got the iPODs). The iRMSD is a special RMSD (standing for intra-molecular distances based RMSD) where the alignments are evaluated using the structural information of the sequences with known structures.


The strength of the iRMSD is its independence from a specific superposition models. When using the iRMSD to evaluate the score of a sequence alignment, one does not need to superpose the two structures and deduce a sequence alignment that will then be compared with the target alignment. In practice, we use a Normalized version of the iRMSD, the NiRMSD that makes it possible to compare alternative alignments of different length. From a structural point of view, the iRMSD has a meaning very similar to the iRMSD and it behaves in a similar fashion from a numerical point of view (similar ranges in Angstroms).


The first step of APDB is to measure the distances between the Ca of each residue and its neighbors. Neighborhood is defined as a sphere of radius -maximum_distance (10 by default). However, by setting -local_mode to 'window', the sphere can be replaced with a window of 1/2 size '-maximum_distance' residues.


Given two aligned residues (X and Y on the Figure) the iRMSD measure is an attempt to estimate the neighborhood support for the XY alignment. This is done by measuring the difference of distances between X and Y and every other pair of aligned residues within the same sphere (W and Z on Figure 1). The iRMSD is obtained by measuring the average Root Mean Square of these differences of distances. The lower the iRMSD, the better the alignment. However, an alignment can obtain a good iRMSD by simply having few aligned residues. To avoid this, the program also reports the NiRMSD= MIN(L1,L2)*iRMSD/Number Considered columns.


How to efficiently use structural information
---------------------------------------------
When it comes to evaluating Multiple Sequence Alignments, nothing beats structural information. To use the methods we describe here, you will need to have at least two structures, similar enough (>60%) to two sequences in your dataset. Here an outline of the best way to proceed:

1) Make sure you include two structures whose sequences are so distantly related that most of the other sequences are intermediates.
2) Align your sequences without using the structural information (i.e. t_coffee, muscle...).
3) Evaluate your alignment with iRMSD (see later in this section); the score will be S1.
4) Realign your sequences, but this time using structural information (Expresso).
5) Measure the score of that alignment; the score wil be S2.

If S1 and S2 are almost similar, it means your distantly related structures were well aligned, and you can expect the intermediate sequences to be well aligned as well. If S2 is much better than S1, you can expect the structures to be well aligned in the second alignment, while there is no guaranty that the alignment of the intermediate sequences has improved as well, although in practice it often does


Evaluating an alignment with the iRMSD package
----------------------------------------------
Let us evaluate the alignment produced by Expresso, using the template_file returned by expresso:


::

  $$: t_coffee -other_pg irmsd sproteases_small.expresso -template_file sproteases_small.template_file

This will deliver a long output. The most interesting bit is at the bottom:

::

  #TOTAL for the Full MSA
   TOTAL EVALUATED: 52.90 %
   TOTAL APDB: 81.59 %  
   TOTAL iRMSD: 0.71 Angs
   TOTAL NiRMSD: 1.33 Angs


APDB is an older measure, less robust than the iRMSD and it is an attempt to estimate the fraction of pairs of residues whose alignment seems to be correct form a structural point of view. The higher APDB, the better the alignment, the lower the NiRMSD, the better the alignment.


Evaluating alternative alignments
---------------------------------
The strength of structure based alignments is that they make it possible to compare alternative alignments. In this case let us consider:


======== ========================= ====== 
Method   File                      NiRMSD 
======== ========================= ====== 
Expresso sproteases_small.expresso 1.33  
T-Coffee sproteases_small.tc_aln   1.35  
ClustalW sproteases_small.cw_aln   1.52  
MAFFT    sproteases_small.mafft    1.36  
MUSCLE   sproteases_small.muscle   1.34  
======== ========================= ====== 

As expected, Expresso delivers the best alignment from a structural point of view. This makes sense, since Expresso explicitely USES structural information. The other figures show us that the structural based alignment is only marginally better than most sequences based alignments. Muscle seems to have a small edge here although the reality is that all these figures are impossible to distinguish with the notable exception of ClustalW


Identifying the most distantly related sequences in your dataset
----------------------------------------------------------------
In order to identify the most distantly related sequences in a dataset, you can use the seq_reformat utility, in order to compare all the sequences two by two and pick up the two having the lowest level of identity:

::

  $$: t_coffee -other_pg seq_reformat sproteases_small.fasta -output sim_idscore\
      | grep TOP |sort -rnk3


This sim_idscore indicates that every pair of sequences will need to be aligned when estimating the similarity. The ouput (below) indicates that the two sequences having the lowest level of identity are AEDAE and MOUSE. It may not be a bad idea to choose these sequences (if possible) for evaluating your MSA.


::

  ...
  TOP 16 10 28.00 sp|P29786|TRY3_AEDAE sp|Q6H321|KLK2_HORSE 28.00
  TOP 16 7 28.00 sp|P29786|TRY3_AEDAE sp|P08246|ELNE_HUMAN 28.00
  TOP 16 1 28.00 sp|P29786|TRY3_AEDAE sp|P08884|GRAE_MOUSE 28.00
  TOP 15 14 27.00 sp|P80015|CAP7_PIG sp|P00757|KLKB4_MOUSE 27.00
  TOP 12 9 27.00 sp|P20160|CAP7_HUMAN sp|Q91VE3|KLK7_MOUSE 27.00
  TOP 9 7 27.00 sp|Q91VE3|KLK7_MOUSE sp|P08246|ELNE_HUMAN 27.00
  TOP 16 2 26.00 sp|P29786|TRY3_AEDAE sp|P21844|MCPT5_MOUSE 26.00



Evaluating an alignment according to your own criterion
=======================================================
Establishing your own criterion
-------------------------------
Any kind of Feature can easily be turned into an evaluation grid. For instance, the protease sequences we have been using here have a well characterized binding site. A possible evaluation can be made as follows. let us consider the Swissprot annotation of the two most distantly related sequences. These two sequences contain the electron relay system of the proteases. We can use it to build an evaluation library: in P29786, the first Histidine is at position 68, while in P21844 this Histidine is on position 66. We can therefore build a library that will check whether these residues are properly aligned in any MSA. The library will look like this:


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


You simply need to cut and paste this library in a file and use this file as a library to measure the concistency between your alignment and the correspondances declared in your library. The following command line also makes it possible to visualy display the agreement between your sequences and the library.


::

  $$: t_coffee -infile sproteases_small.aln -lib charge_relay_lib.tc_lib -score \
 -output html



********************************************
Trees Based on Protein 3D Structures: T-RMSD 
********************************************

*Tree estimation procedure based on the comparison of internal distances*

T-RMSD makes it possible to estimate a tree using either contact conservation or differences in internal distances as a measure of similarity bewtween protein or RNA sequences. The trees thus estimated can be bootsrapped or further analyzed like regular phylogenetic trees. T-RMSD also makes it possible to estimate the local support of any structural alignment (i.e. each individual column) for either a full tree or any pre-defined sub-group contained within the dataset. 

Generating a tree based on structural distances
===============================================
This option makes it possible to estimate a tree while taking into account the variation of intra-molecular distances within the considered sequences. The following call will generate a 100 replicate nj trees using the difference of distances between pairs of aligned residues, at a maximum cut-off of 15A. Columns with less than 50% residues are ignored


Input:

* ``aln``: Multiple Sequence Alignment in FASTA, MSA or MSF
* ``template``: FASTA name list with templates: ``>name _P_ template``

:: 

  t_coffee -other_pg seq_reformat -in <aln> -in2 <template> -action +tree replic\
 ates 100  +evaluate3D distances +tree2bs first -output newick -out tree.dnd


Output: 

* ``tree.dnd``: Tree in newick format with bootstrap support   

It is possible to control default parameters using the following extended command line

::

  t_coffee -other_pg seq_reformat -in <aln> -in2 <template> -action +tree replic\
 ates 100 gap 0.5 mode nj  +evaluate3D distances 15 +tree2bs first -output newick\
  -out tree.dnd

.. warning: sequences without 3D structure will be excluded from the analysis and from the final output


Generating a Tree based on contact conservation
================================================

This option makes it possible to estimate a tree while taking into account the variation of contact conservation within the considered sequences. This call will generate a 100 replicate nj trees using as a distance metrics the fraction of contacts conserved between pairs of aligned residues, at a maximum cut-off of 1.2 A between VdW radius and ignoring the 3 closests neighbors. Columns with less than 50% residues are ignored. For sequences without 3D information, the strike contact potential is used instead (Watson and crick base pairing propensity for RNA).

:: 

  t_coffee -other_pg seq_reformat -in <seq.aln> -in2 <seq.template> -action +tre\
 e replicates 100  +evaluate3D contacts +tree2bs first -output newick -out tree.d\
 nd


Output: 

* ``tree.dnd``: Tree in newick format  

It is possible to control default parameters using the following extended command line:

::

  seq_reformat -in <aln> -in2 <template> -action +tree replicates 100 gap 0.5 mo\
 de nj  +evaluate3D contacts 1.2 3 +tree2bs first -output newick -out tree.dnd

.. warning: the procedure requires at least 1 sequence with a known 3D structure or with contact information.



Visulizing 3D Conservation
==========================

This same procedure can be used to visualize either intramolecular distance conservation or contact conservation:

::

  seq_reformat -in CRD.aln -in2 CRD.template -action +evaluate3D distances -output score_html 
  seq_reformat -in CRD.aln -in2 CRD.template -action +evaluate3D distances -output score_ascii
  seq_reformat -in CRD.aln -in2 CRD.template -action +evaluate3D distances -output score_raw

Output:
* ``score_raw``: Tabulated dump of the numerical values associated with every residue, every sequence and every column of the considered alignment.


Identification of positions 
===========================

If you have a well defined sub-group of sequences (i.e. domains having the same function, same specificty, etc...), it is possible to estimate which columns yield the best support using the following command,

Input:
* ``group.fasta``: A Fasta formatted list of the sequences that form the group whose support you want to analyze

::

 seq_reformat -in <seq.aln> -in2 <seq.template> -action +tree replicates columns\
   +evaluate3D  distances +evaluateTree <group.fasta> -output score_html -out <al\
 n.html>

Output:
* ``aln.score_html`` Colored version of your MSA indicating the sequences that best contribute to your clustering.


Evaluating Clustering capacities
=================================

If you want to check the capacity of an algorithm to bring related sequences within mono-phyletic groups, you should name your sequences according to the group they belong to (XXXX_1, YYYYY_1, ZZZZ_2, KKKK_2, for members of _1 and _2, etc) and use the following evaluation procedure. The output will be the number of monophyletic groups containing sequences belonging to the same group:

The tree can be pre-computed
:: 

  seq_reformat -in <tree> +tree2collapse groups 4 +print nseq -output no

Or it can be computed on the fly
:: 

  seq_reformat -in <aln> -in2 <template> -action +tree replicates 100  +evaluate\
 3D  distances 15 +tree2bs first +tree2collapse groups 4 +print nseq -output no



****************************************
Integrating External Methods In T-Coffee
****************************************
The real power of T-Coffee is its ability to seamlessly combine many methods into one. While we try to integrate as many methods as we can in the default distribution, we do not have the means to be exhaustive and if you desperately need your favourite method to be integrated, you will need to bite the bullet ...


What are the methods already integrated in T-Coffee
===================================================
Although, it does not necessarily do so explicitly, T-Coffee always end up combining libraries. Libraries are collections of pairs of residues. Given a set of libraries, T-Coffee makes an attempt to assemble the alignment with the highest level of consistence. You can think of the alignment as a timetable. Each library pair would be a request from students or teachers, and the job of T-Coffee would be to assemble the time table that makes as many people as possible happy...


In T-Coffee, methods replace the students/professors as constraints generators. These methods can be any standard/non standard alignment methods that can be used to generate alignments (pairwise, most of the time). These alignments can be viewed as collections of constraints that must be fit within the final alignment. Of course, the constraints do not have to agree with one another...


This section shows you what are the vailable method in T-Coffee, and how you can add your own methods, either through direct parameterization or via a perl script. There are two kinds of methods: the internal and the external. For the internal methods, you simply need to have T-Coffee up and running. The external methods will require you to install a package.


List of INTERNAL methods
------------------------
Built in methods can be requested using the following names:

-proba_pair: Adapted from Probcons, this method [the current default] uses a pair HMM to compute a pairwise alignment with a bi-phasic gap penalty.

-fast_pair: Makes a global fasta style pairwise alignment. For proteins, matrix=blosum62mt, gep=-1, gop=-10, ktup=2. For DNA, matrix=idmat (id=10), gep=-1, gop=-20, ktup=5. Each pair of residue is given a score function of the weighting mode defined by -weight.

-slow_pair: Identical to fast pair, but does a full dynamic programming, using the myers and miller algorithm. This method is recommended if your sequences are distantly related.

-ifast_pair

-islow_pair:  Makes a global fasta alignmnet using the previously computed pairs as a library. `i` stands for iterative. Each pair of residue is given a score function of the weighting mode defined by -weight. The Library used for the computation is the one computed before the method is used. The resullt is therefore dependant on the order in methods and library are set via the -in flag.

-align_pdb_pair: Uses the align_pdb routine to align two structures. The pairwise scores are those returnes by the align_pdb program. If a structure is missing, fast_pair is used instead. Each pair of residue is given a score function defined by align_pdb. [UNSUPPORTED]

-lalign_id_pair: Uses the ten top non intersecting local alignments, as delivered by lalign. Each alignement is weighted with its average percent identity.

-lalign_rs_s_pair: Same as above but does also does self comparison and uses the lalign raw_score (s stands for self). This is needed when extracting repeats.

-Matrix Amy: matrix can be requested, simply indicate as a method the name of the matrix preceded with an X (i.e. Xpam250mt). If you indicate such a matrix, all the other methods will simply be ignored, and a standard fast progressive alignment will be computed. If you want to change the substitution matrix used by the methods, use the -matrix flag.

-cdna_fast_pair: This method computes the pairwise alignment of two cDNA sequences. It is a fast_pair alignment that only takes into account the amino-acid similarity and uses different penalties for amino-acid insertions and frameshifts. This alignment is turned into a library where matched nucleotides receive a score equql to the average level of identity at the amino-acid level. This mode is intended to clean cDNA obtained from ESTs, or to align pseudo-genes.


.. warning:: This last method is currently no longer supported.


List of EXTERNAL methods (plug-ins)
-----------------------------------
The following methods are external. They correspond to packages developed by other groups that you may want to run within T-Coffee. We are very open to extending these options and we welcome any request to add an extra interface. The following table lists the methods that can be used as plug-ins:


::

  Package Where From
  =====================================================================
  ClustalW can interact with t_coffee
  ---------------------------------------------------------------------
  POA  http://www.bioinformatics.ucla.edu/poa/
  ---------------------------------------------------------------------
  MUSCLE http://www.bioinformatics.ucla.edu/poa/
  ---------------------------------------------------------------------
  ProbCons http://probcons.stanford.edu/
  ---------------------------------------------------------------------
  MAFFT http://www.biophys.kyoto-u.ac.jp/~katoh/programs/align/mafft/
  ---------------------------------------------------------------------
  Dialign-T http://dialign-t.gobics.de/
  ---------------------------------------------------------------------
  PCMA ftp://iole.swmed.edu/pub/PCMA/
  ---------------------------------------------------------------------
  SAP obtained from W. Taylor, NIMR-MRC
  ---------------------------------------------------------------------
  BLAST  http://www.ncbi.nih.nlm.gov
  ---------------------------------------------------------------------
  Fugue http://www-cryst.bioc.cam.ac.uk/fugue/download.html
  ---------------------------------------------------------------------

Once installed, most of these methods can be used as either pairwise or multiple alignment methods. Note that all these methods use Blosum62 as a default. 
clustalw_pair uses clustalw (default parameters) to align two sequences. Each pair of residue is given a score function of the weighting mode defined by -weight. 
clustalw_msa makes a multiple alignment using ClustalW and adds it to the library. Each pair of residue is given a score function of the weighting mode defined by -weight.

probcons_pair Probcons package: probcons.stanford.edu/.
probcons_msa idem.

muscle_pair Muscle package www.drive5.com/muscle/ .
muscle_msa idem.

mafft_pair www.biophys.kyoto-u.ac.jp/~katoh/programs/align/mafft/ .
mafft_msa idem.

pcma_msa pcma package
pcma_pair pcma package

poa_msa poa package
poa_pair poa package

dialignt_pair dialignt package
dialignt_msa pcma package


sap_pair uses SAP to align two structures. Each pair of residue is given a score function defined by sap. You must have sap installed on your system to use this method.

fugue_pair uses a standard fugue installation to make a sequence /structure alignment. Fugue installation must be standard. It does not have to include all the fugue packages but only:

1- joy, melody, fugueali, sstruc, hbond
2-copy fugue/classdef.dat /data/fugue/SUBST/classdef.dat

OR

Setenv MELODY_CLASSDEF=<location>
Setenv MELODY_SUBST=fugue/allmat.dat


All the configuration files must be in the right location. To request a method, see the -in or the -method flag. For instance, if you wish to request the use of fast_pair and lalign_id_pair (the current default):

::

  $$: t_coffee -seq sample_seq1.fasta -method fast_pair,lalign_id_pair


Modifying the parameters of internal and external Methods
=========================================================
Internal Methods
----------------
It is possible to modify on the fly the parameters of hard coded methods:


::

  $$: t_coffee sample_seq1.fasta -method slow_pair@EP@MATRIX@pam250mt@GOP@-10@GE\
 P@-1



EP stands for Extra parameters. These parameters will superseed any other parameters.


External methods
----------------
External methods receive a command line built with the information provided via the parameter file (see next heading). It is possible to produce such a parameter file and to modify it in order to modify the commands passed to the methods.


The passed command is built as follows:


<EXECUTABLE><PARAM1><IN_FLAG><seq_file><PARAM2><OUT_FLAG><outname><PARAM>


You should know what is the best place for squizing your extra parameters. It will depend on the application, although PARAM2 is usually a good guess. Now if you want, for instance to modify the gap penalty of clustalw, you can try the following:


::

  $$: t_coffee sample_seq1.fasta -method clustalw_msa@EP@PARAM2@-GAPOPEN%e100%s-\
 GAPEXT%e10



@EP is here to indicate that you will pass an extra parameter


@PARAM1 is the name of this parameter


the next filed is the parameter itself, where:


%s replace spaces


%e replaces the equal sign


Of course, you must know the command line of the program you are trying to modify (clustalw in this case).


Integrating external methods
============================
If the method you need is not already included in T-Coffee, you will need to integrate it yourself. We give you here some guidelines on how to do so.


Direct access to external methods
---------------------------------
A special method exists in T-Coffee that can be used to invoke any existing program:


::

  $$: t_coffee sample_seq1.fasta -method=em@clustalw@pairwise



In this context, Clustalw is a method that can be ran with the following command line:


::

  method -infile=<infile> -outfile=<outfile>



Clustalw can be replaced with any method using a similar syntax. If the program you want to use cannot be run this way, you can either write a perl wrapper that fits the bill or write a tc_method file adapted to your program (cf next section).


This special method (em, external method) uses the following syntax:


::

  em@<method>@<aln_mode:pairwises_pairwise|multiple>



Customizing an external method (with parameters) for T-Coffee
-------------------------------------------------------------
T-Coffee can run external methods, using a tc_method file that can be used in place of an established method. Two such files are incorporated in T-Coffee. You can dump them and customize them according to your needs:


For instance if you have ClustalW installed, you can use the following file to run the


::

  $$: t_coffee -other_pg unpack_clustalw_method.tc_method

  $$: t_coffee -other_pg unpack_generic_method.tc_method



The second file (generic_method.tc_method) contains many hints on how to customize your new method. The first file is a very straightforward example on how to have t_coffee to run Clustalw with a set of parameters you may be interested in:


::

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



This configuration file will cause T-Coffee to emit the following system call:


::

  clustalw -INFILE=tmpfile1 -OUTFILE=tmpfile2 -gapopen=-10



Note that ALN_MODE instructs t_coffee to run clustalw on every pair of sequences (cf generic_method.tc_method for more details). The tc_method files are treated like any standard established method in T-Coffee. For instance, if the file clustalw_method.tc_method is in your current directory, run:


::

  $$: t_coffee sample_seq1.fasta -method clustalw_method.tc_method


Managing a collection of method files
-------------------------------------
It may be convenient to store all the method files in a single location on your system. By default, t_coffee will go looking into the directory ~/.t_coffee/methods/. You can change this by either modifying the METHODS_4_TCOFFEE in define_headers.h (and recompile) or by modifying the environement variable METHODS_4_TCOFFEE.


Advanced method integration
===========================
It may sometimes be difficult to customize the program you want to use through a tc_method file. In that case, you may rather use an external perl_script to run your external application. This can easily be achieved using the generic_method.tc_method file.


::

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

When you run this method:


::

  $$: t_coffee -other_pg unpack_generic_method.tc_method
  $$: t_coffee sample_seq1.fasta -method generic_method.tc_method


T-Coffee runs the script tc_generic_method.pl on your data. It also provides the script with parameters. In this case -method clustalw indicates that the script should run clustalw on your data. The script tc_generic_method.pl is incorporated in t_coffee. Over the time, this script will be the place where novel methods will be integrated. It will be used to run the script tc_generic_method.pl. The file tc_generic_method.pl is a perl file, automatically generated by t_coffee. Over the time this file will make it possible to run all available methods. You can dump the script using the following command:

::

  $$: t_coffee -other_pg=unpack_tc_generic_method.pl


.. note:: If there is a copy of that script in your local directory, that copy will be used in place of the internal copy of T-Coffee.


The mother of all method files...
---------------------------------
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
  * The method is passed to the shell using the following
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
  *IN_FLAG S no_name ->no flag
  *IN_FLAG S &amp;nbsp-in&amp;nbsp -> ' -in '
  *
  IN_FLAG -infile=
  *
  *******************************************************
  * OUT_FLAG *
  *******************************************************
  *OUT_FLAG
  *flag indicating the name of the out-coming data
  *same conventions as IN_FLAG
  *OUT_FLAG S no_name ->no flag
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
By default, the alignment produced by your method will be weighted according to its percent identity. However, this can be customized via the WEIGHT parameter. The WEIGHT parameter supports all the values of the -weight flag. The only difference is that the -weight value thus declared will only be applied onto your method. If needed you can also modify on the fly the WEIGHT value of your method:

::

  $$: t_coffee sample_seq1.fasta -method slow_pair@WEIGHT@OW2
  

Will overweight by a factor 2 the weight of slow_pair.

::

  $$: t_coffee sample_seq1.fasta -method slow_pair@WEIGHT@250


Will cause every pair of slow_pair to have a weight equal to 250


Plug-out: using T-Coffee as a plug-in
=====================================
Just because it enjoys enslaving other methods as plug-ins, does not mean that T-Coffee does not enjoy being incorporated within other packages. We try to give as much support as possible to anyone who wishes to incorporate T-Coffee in an alignment pipeline. If you want to do so, please work out some way to incorporate T-Coffee in your script . If you need some help along the ways, do not hesitate to ask, as we will always be happy to either give assistance, or even modify the package so that it accomodates as many needs as possible. Once that procedure is over, set aside a couple of input files with the correct parameterisation and send them to us. These will be included as a distribution test, to insure that any further distribution remains compliant with your application. We currently support:

::

  Package Where From
  ==============================================================
  Marna  http://www.bio.inf.unijena.de/Software/MARNA/download
  --------------------------------------------------------------



Creating your own T-Coffee libraries
====================================
If the method you want to use is not integrated, or impossible to integrate, you can generate your own libraries, either directly or by turning existing alignments into libraries. You may also want to precompute your libraries, in order to combine them at your convenience.


Using pre-computed alignments
-----------------------------
If the method you wish to use is not supported, or if you simply have the alignments, the simplest thing to do is to generate yourself the pairwise/multiple alignments, in FASTA, ClustalW, msf or Pir format and feed them into t_coffee using the -in flag:


::

  $$: t_coffee -aln=sample_aln1_1.aln,sample_aln1_2.aln -outfile=combined_aln.al\
 n



Customizing the weighting scheme
--------------------------------
The previous integration method forces you to use the same weighting scheme for each alignment and the rest of the libraries generated on the fly. This weighting scheme is based on global pairwise sequence identity. If you want to use a more specific weighting scheme with a given method, you should either:


generate your own library (cf next section)


convert your aln into a lib, using the -weight flag:


::

  $$: t_coffee -aln sample_aln1.aln -out_lib=test_lib.tc_lib -lib_only -weight=s\
 im_pam250mt

  $$: t_coffee -aln sample_aln1.aln -lib test_lib.tc_lib -outfile=outaln

  $$: t_coffee -aln=sample_aln1_1.aln,sample_aln1_2.aln -method= fast_pair,lalig\
  n_id_pair -outfile=out_aln



Generating your own libraries
-----------------------------
This is suitable if you have local alignments, or very detailed information about your potential residue pairs, or if you want to use a very specific weighting scheme. You will need to generate your own libraries, using the format described in the last section.


You may also want to pre-compute your libraries in order to save them for further use. For instance, in the following example, we generate the local and the global libraries and later re-use them for combination into a multiple alignment.


::

  $$: t_coffee sample_seq1.fasta -method slow_pair -out_lib slow_pair_seq1.tc_li\
 b -lib_only

  $$: t_coffee sample_seq1.fasta -method lalign_id_pair -out_lib lalign_id_pair_seq1.tc_lib -lib_only



Once these libraries have been computed, you can then combine them at your convenience in a single MSA. Of course you can decide to only use the local or the global library


::

  $$: t_coffee sample_seq1.fasta -lib lalign_id_pair_seq1.tc_lib, slow_pair_seq1\
 .tc_lib









