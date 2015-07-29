.. raw:: html

   <div class="WordSection1">

.. raw:: html

   <div
   style="mso-element:frame;mso-element-wrap:no-wrap-beside;mso-height-rule:
   exactly">

+--------------------------------------------------------------------------+
| .T-coffee Tutorial                                                       |
+--------------------------------------------------------------------------+

.. raw:: html

   </div>

Centre National De LA Recherche scientifique
 CENTRO DE REGULACCIO GENOMICA, Barcelona

.. raw:: html

   <div
   style="mso-element:para-border-div;border:none;border-top:solid windowtext 1.0pt;
   mso-border-top-alt:solid windowtext .75pt;padding:1.0pt 0cm 0cm 0cm">

Cédric Notredame
 www.tcoffee.org

.. raw:: html

   </div>

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid white 1.0pt;mso-border-alt:
   solid white .75pt;padding:31.0pt 31.0pt 31.0pt 31.0pt;background:#E5E5E5;
   mso-shading:windowtext;mso-pattern:gray-10 auto">

T-Coffee:
 Cheat Sheet

Tutorial and FAQ

.. raw:: html

   </div>

.. raw:: html

   <div
   style="mso-element:para-border-div;border:none;border-top:solid windowtext 1.0pt;
   mso-border-top-alt:solid windowtext .75pt;padding:1.0pt 0cm 0cm 0cm">

 

.. raw:: html

   </div>

.. raw:: html

   </div>

.. raw:: html

   <div class="WordSection2">

.. raw:: html

   <div
   style="mso-element:para-border-div;border:none;border-bottom:solid gray 1.0pt;
   mso-border-bottom-alt:solid gray .75pt;padding:0cm 0cm 14.0pt 0cm">

T-Coffee Tutorial
 (Version 6.18, August 2008)
 T-Coffee, PSI-Coffee
 3D-Coffee/Expresso
 M-Coffee
 R-Coffee
 APDB and iRMSD

.. raw:: html

   </div>

ã Cédric Notredame, Centro de Regulaccio Genomica and Centre National de
la Recherche Scientifique, France

.. raw:: html

   </div>

.. raw:: html

   <div class="WordSection3">

`Cheat Sheet: T-Coffee. 6 <#_Toc235985428>`__

`Proteins 6 <#_Toc235985429>`__

`DNA. 6 <#_Toc235985430>`__

`RNA. 6 <#_Toc235985431>`__

`Memory Problems 6 <#_Toc235985432>`__

`Before You Start….. 8 <#_Toc235985433>`__

`Foreword. 8 <#_Toc235985434>`__

`Pre-Requisite 8 <#_Toc235985435>`__

`Getting the Example Files of the Tutorial 9 <#_Toc235985436>`__

`What Is  T-COFFEE ?. 10 <#_Toc235985437>`__

`What is T-Coffee?. 10 <#_Toc235985438>`__

`What does it do?. 10 <#_Toc235985439>`__

`What can it align?. 10 <#_Toc235985440>`__

`How can I use it?. 10 <#_Toc235985441>`__

`Is There an Online Server 11 <#_Toc235985442>`__

`Is T-Coffee different from ClustalW?. 11 <#_Toc235985443>`__

`Is T-Coffee very accurate?. 11 <#_Toc235985444>`__

`What T-Coffee Can and Cannot do for you ….. 12 <#_Toc235985445>`__

`(NOT) Fetching Sequences 12 <#_Toc235985446>`__

`Aligning Sequences 12 <#_Toc235985447>`__

`Combining Alignments 12 <#_Toc235985448>`__

`Evaluating Alignments 12 <#_Toc235985449>`__

`Combining Sequences and Structures 12 <#_Toc235985450>`__

`Identifying Occurrences of a Motif: Mocca. 13 <#_Toc235985451>`__

`How Does T-Coffee works 13 <#_Toc235985452>`__

`Preparing Your Data: Reformatting and Trimming With seq\_reformat
15 <#_Toc235985453>`__

`Seq\_reformat 15 <#_Toc235985454>`__

`Accessing the T-Coffee Reformatting Utility 15 <#_Toc235985455>`__

`An overview of seq\_reformat 16 <#_Toc235985456>`__

`Reformatting your data. 16 <#_Toc235985457>`__

`Changing MSA formats 16 <#_Toc235985458>`__

`Dealing with Non-automatically recognized formats
16 <#_Toc235985459>`__

`Automated Sequence Edition. 16 <#_Toc235985460>`__

`Removing the gaps from an alignment 16 <#_Toc235985461>`__

`Changing the case of your sequences 16 <#_Toc235985462>`__

`Changing the case of specific residues 17 <#_Toc235985463>`__

`Changing the case depending on the score 17 <#_Toc235985464>`__

`Protecting Important Sequence Names 17 <#_Toc235985465>`__

`Colouring/Editing Residues in an Alignment 18 <#_Toc235985466>`__

`Coloring specific types of residues 18 <#_Toc235985467>`__

`Coloring a specific residue of a specific sequence
18 <#_Toc235985468>`__

`Coloring according to the conservation. 18 <#_Toc235985469>`__

`Colouring/Editing residues in an Alignment Using a Cache
19 <#_Toc235985470>`__

`Overview. 19 <#_Toc235985471>`__

`Preparing a Sequence or Alignment Cache 19 <#_Toc235985472>`__

`Preparing a Library Cache 20 <#_Toc235985473>`__

`Coloring an Alignment using a cache 21 <#_Toc235985474>`__

`Changing the default colors 21 <#_Toc235985475>`__

`Evaluating an alignment and producing a cache 22 <#_Toc235985476>`__

`Evaluating an alignment with T-Coffee 22 <#_Toc235985477>`__

`Evaluating the level of conservation with a substitution matrix
22 <#_Toc235985478>`__

`Selective Reformatting. 23 <#_Toc235985479>`__

`Removing gapped columns 23 <#_Toc235985480>`__

`Selectively turn some residues to lower case 23 <#_Toc235985481>`__

`Selectively modifying residues 24 <#_Toc235985482>`__

`Keeping only the best portion of an alignment 24 <#_Toc235985483>`__

`Extracting Portions of Dataset 25 <#_Toc235985484>`__

`Extracting The High Scoring Blocks 25 <#_Toc235985485>`__

`Extracting Sequences According to a Pattern. 26 <#_Toc235985486>`__

`Extracting Sequences by Names 26 <#_Toc235985487>`__

`Removing Sequences by Names 27 <#_Toc235985488>`__

`Extracting Blocks Within Alignment 27 <#_Toc235985489>`__

`Concatenating Alignments 28 <#_Toc235985490>`__

`Analyzing your Multiple Sequence Alignment 28 <#_Toc235985491>`__

`Estimating the diversity in your alignment 28 <#_Toc235985492>`__

`Reducing and improving your dataset 28 <#_Toc235985493>`__

`Extracting the N most informative sequences 29 <#_Toc235985494>`__

`Extracting all the sequences less than X% identical
29 <#_Toc235985495>`__

`Speeding up the process 29 <#_Toc235985496>`__

`Forcing Specific Sequences to be kept 30 <#_Toc235985497>`__

`Identifying and Removing Outlayers 31 <#_Toc235985498>`__

`Chaining Important Sequences 31 <#_Toc235985499>`__

`Manipulating DNA sequences 31 <#_Toc235985500>`__

`Translating DNA sequences into Proteins 31 <#_Toc235985501>`__

`Back-Translation With the Bona-Fide DNA sequences
32 <#_Toc235985502>`__

`Finding the Bona-Fide Sequences for the Back-Translation.
32 <#_Toc235985503>`__

`Guessing Your Back Translation. 32 <#_Toc235985504>`__

`Fetching a Structure 32 <#_Toc235985505>`__

`Fetching a PDB structure 32 <#_Toc235985506>`__

`Fetching The Sequence of a PDB structure 33 <#_Toc235985507>`__

`Adapting extract\_from\_pdb to your own environment
33 <#_Toc235985508>`__

`Manipulating RNA sequences with seq\_reformat 34 <#_Toc235985509>`__

`Producing a Stockholm output: adding predicted secondary structures
34 <#_Toc235985510>`__

`Producing a consensus structure 34 <#_Toc235985511>`__

`Adding a consensus structure to an alignment 34 <#_Toc235985512>`__

`Analyzing an alifold secondary structure prediction.
35 <#_Toc235985513>`__

`Analyzing matching columns 35 <#_Toc235985514>`__

`Visualizing compensatory mutations 36 <#_Toc235985515>`__

`Handling gapped columns 36 <#_Toc235985516>`__

`Comparing alternative folds 36 <#_Toc235985517>`__

`Manipulating Phylogenetic Trees with seq\_reformat
37 <#_Toc235985518>`__

`Producing phylogenetic trees 37 <#_Toc235985519>`__

`Comparing two phylogenetic trees 38 <#_Toc235985520>`__

`Scanning Phylogenetic Trees 38 <#_Toc235985521>`__

`Pruning Phylogenetic Trees 39 <#_Toc235985522>`__

`Building Multiple Sequence Alignments. 40 <#_Toc235985523>`__

`How to generate The Alignment You Need?. 40 <#_Toc235985524>`__

`What is a Good Alignment?. 40 <#_Toc235985525>`__

`The Main Methods and their Scope 41 <#_Toc235985526>`__

`Choosing The Right Package 42 <#_Toc235985527>`__

`Computing Multiple Sequence Alignments With T-Coffee
43 <#_Toc235985528>`__

`Computing Very accurate (but slow) alignments with PSI-Coffee
43 <#_Toc235985529>`__

`A Simple Multiple Sequence Alignment 43 <#_Toc235985530>`__

`Controlling the Output Format 43 <#_Toc235985531>`__

`Computing a Phylogenetic tree 43 <#_Toc235985532>`__

`Using Several Datasets 44 <#_Toc235985533>`__

`How Good is Your Alignment 44 <#_Toc235985534>`__

`Doing it over the WWW.. 44 <#_Toc235985535>`__

`Aligning Many Sequences 45 <#_Toc235985536>`__

`Aligning Very Large Datasets with Muscle 45 <#_Toc235985537>`__

`Aligning Very Large Alignments with Mafft 45 <#_Toc235985538>`__

`Aligning Very Large Alignments with T-Coffee 45 <#_Toc235985539>`__

`Shrinking Large Alignments With T-Coffee 45 <#_Toc235985540>`__

`Modifying the default parameters of T-Coffee 45 <#_Toc235985541>`__

`Changing the Substitution Matrix 46 <#_Toc235985542>`__

`Comparing Two Alternative Alignments 46 <#_Toc235985543>`__

`Changing Gap Penalties 48 <#_Toc235985544>`__

`Can You Guess The Optimal Parameters?. 49 <#_Toc235985545>`__

`Using Many Methods at once 49 <#_Toc235985546>`__

`Using All the Methods at the Same Time: M-Coffee 49 <#_Toc235985547>`__

`Using Selected Methods to Compute your MSA. 50 <#_Toc235985548>`__

`Combining pre-Computed Alignments 51 <#_Toc235985549>`__

`Aligning Profiles 51 <#_Toc235985550>`__

`Using Profiles as templates 51 <#_Toc235985551>`__

`Aligning One sequence to a Profile 51 <#_Toc235985552>`__

`Aligning Many Sequences to a Profile 52 <#_Toc235985553>`__

`Aligning Other Types of Sequences 52 <#_Toc235985554>`__

`Splicing variants 52 <#_Toc235985555>`__

`Aligning DNA sequences 53 <#_Toc235985556>`__

`Aligning RNA sequences 53 <#_Toc235985557>`__

`Noisy Coding DNA Sequences….. 53 <#_Toc235985558>`__

`Using Secondary Structure Predictions: 55 <#_Toc235985559>`__

`Single Sequence prediction. 55 <#_Toc235985560>`__

`Multiple Sequence Predictions 55 <#_Toc235985561>`__

`Incorporation of the prediction in the alignment 56 <#_Toc235985562>`__

`Using other secondary structure predictions 56 <#_Toc235985563>`__

`Output of the prediction. 57 <#_Toc235985564>`__

`Combining Sequences and 3D-Structures. 58 <#_Toc235985565>`__

`If you are in a Hurry: Expresso. 58 <#_Toc235985566>`__

`What is Expresso?. 58 <#_Toc235985567>`__

`Using Expresso. 59 <#_Toc235985568>`__

`Aligning Sequences and Structures 59 <#_Toc235985569>`__

`Mixing Sequences and Structures 59 <#_Toc235985570>`__

`Using Sequences only 60 <#_Toc235985571>`__

`Aligning Profile using Structural Information. 60 <#_Toc235985572>`__

`How Good Is Your Alignment ?. 61 <#_Toc235985573>`__

`Evaluating Alignments with The CORE index. 61 <#_Toc235985574>`__

`Computing the Local CORE Index 61 <#_Toc235985575>`__

`Computing the CORE index of any alignment 61 <#_Toc235985576>`__

`Filtering Bad Residues 61 <#_Toc235985577>`__

`Filtering Gap Columns 62 <#_Toc235985578>`__

`Evaluating an Alignment Using Structural Information: APDB and iRMSD..
63 <#_Toc235985579>`__

`What is the iRMSD?. 63 <#_Toc235985580>`__

`How to Efficiently Use Structural Information. 64 <#_Toc235985581>`__

`Evaluating an Alignment With the iRMSD Package 64 <#_Toc235985582>`__

`Evaluating Alternative Alignments 64 <#_Toc235985583>`__

`Identifying the most distantly related sequences in your dataset
65 <#_Toc235985584>`__

`Evaluating an Alignment according to your own Criterion.
65 <#_Toc235985585>`__

`Establishing Your Own Criterion. 65 <#_Toc235985586>`__

`Integrating External Methods In T-Coffee. 67 <#_Toc235985587>`__

`What Are The Methods Already Integrated in T-Coffee
67 <#_Toc235985588>`__

`List of INTERNAL Methods 67 <#_Toc235985589>`__

`Plug-In: Using Methods Integrated in T-Coffee 68 <#_Toc235985590>`__

`Modifying the parameters of Internal and External Methods
70 <#_Toc235985591>`__

`Internal Methods 70 <#_Toc235985592>`__

`External Methods 70 <#_Toc235985593>`__

`Integrating External Methods 71 <#_Toc235985594>`__

`Direct access to external methods 71 <#_Toc235985595>`__

`Customizing an external method (with parameters) for T-Coffee
71 <#_Toc235985596>`__

`Managing a collection of method files 72 <#_Toc235985597>`__

`Advanced Method Integration. 72 <#_Toc235985598>`__

`The Mother of All method files….. 74 <#_Toc235985599>`__

`Weighting your Method. 75 <#_Toc235985600>`__

`Plug-Out: Using T-Coffee as a Plug-In. 76 <#_Toc235985601>`__

`Creating Your Own T-Coffee Libraries 76 <#_Toc235985602>`__

`Using Pre-Computed Alignments 76 <#_Toc235985603>`__

`Customizing the Weighting Scheme 76 <#_Toc235985604>`__

`Generating Your Own Libraries 77 <#_Toc235985605>`__

`Frequently Asked Questions. 78 <#_Toc235985606>`__

`Abnormal Terminations and Wrong Results 78 <#_Toc235985607>`__

`Q: The program keeps crashing when I give my sequences
78 <#_Toc235985608>`__

`Q: The default alignment is not good enough. 78 <#_Toc235985609>`__

`Q: The alignment contains obvious mistakes 79 <#_Toc235985610>`__

`Q: The program is crashing. 79 <#_Toc235985611>`__

`Q: I am running out of memory 79 <#_Toc235985612>`__

`Input/Output Control 79 <#_Toc235985613>`__

`Q: How many Sequences can t\_coffee handle 79 <#_Toc235985614>`__

`Q: Can I prevent the Output of all the warnings?.
79 <#_Toc235985615>`__

`Q: How many ways to pass parameters to t\_coffee?.
79 <#_Toc235985616>`__

`Q: How can I change the default output format?. 80 <#_Toc235985617>`__

`Q: My sequences are slightly different between all the alignments.
80 <#_Toc235985618>`__

`Q: Is it possible to pipe stuff OUT of t\_coffee?.
80 <#_Toc235985619>`__

`Q: Is it possible to pipe stuff INTO t\_coffee?. 80 <#_Toc235985620>`__

`Q: Can I read my parameters from a file?. 80 <#_Toc235985621>`__

`Q: I want to  decide myself on the name of the output files!!!
81 <#_Toc235985622>`__

`Q: I want to use the sequences in an alignment file
81 <#_Toc235985623>`__

`Q: I only want to produce a library 81 <#_Toc235985624>`__

`Q: I want to turn an alignment into a library 81 <#_Toc235985625>`__

`Q: I want to concatenate two libraries 81 <#_Toc235985626>`__

`Q: What happens to the gaps when an alignment is fed to T-Coffee
82 <#_Toc235985627>`__

`Q: I cannot print the html graphic display!!! 82 <#_Toc235985628>`__

`Q: I want to output an html file and a regular file
82 <#_Toc235985629>`__

`Q: I would like to output more than one alignment format at the same
time 82 <#_Toc235985630>`__

`Alignment Computation. 83 <#_Toc235985631>`__

`Q: Is T-Coffee the best? Why Not Using Muscle, or Mafft, or
ProbCons???. 83 <#_Toc235985632>`__

`Q: Can t\_coffee align Nucleic Acids ???. 83 <#_Toc235985633>`__

`Q: I do not want to compute the alignment. 83 <#_Toc235985634>`__

`Q: I would like to force some residues to be aligned.
83 <#_Toc235985635>`__

`Q: I would like to use structural alignments. 84 <#_Toc235985636>`__

`Q: I want to build my own libraries. 84 <#_Toc235985637>`__

`Q: I want to use my own tree 84 <#_Toc235985638>`__

`Q: I want to align coding DNA. 85 <#_Toc235985639>`__

`Q: I do not want to use all the possible pairs when computing the
library 85 <#_Toc235985640>`__

`Q: I only want to use specific pairs to compute the library
85 <#_Toc235985641>`__

`Q: There are duplicates or quasi-duplicates in my set
85 <#_Toc235985642>`__

`Using Structures and Profiles 86 <#_Toc235985643>`__

`Q: Can I align sequences to a profile with T-Coffee?.
86 <#_Toc235985644>`__

`Q: Can I align sequences Two or More Profiles?. 86 <#_Toc235985645>`__

`Q: Can I align two profiles according to the structures they contain?.
86 <#_Toc235985646>`__

`Q: T-Coffee becomes very slow when combining sequences and structures
86 <#_Toc235985647>`__

`Q: Can I use a local installation of PDB?. 87 <#_Toc235985648>`__

`Alignment Evaluation. 87 <#_Toc235985649>`__

`Q: How good is my alignment?. 87 <#_Toc235985650>`__

`Q: What is that color index?. 87 <#_Toc235985651>`__

`Q: Can I evaluate alignments NOT produced with T-Coffee?.
88 <#_Toc235985652>`__

`Q: Can I Compare Two Alignments?. 88 <#_Toc235985653>`__

`Q: I am aligning sequences with long regions of very good overlap.
88 <#_Toc235985654>`__

`Q: Why is T-Coffee changing the names of my sequences!!!!
89 <#_Toc235985655>`__

`Improving Your Alignment 89 <#_Toc235985656>`__

`Q: How Can I Edit my Alignment Manually?. 89 <#_Toc235985657>`__

`Q: Have I Improved or Not my Alignment?. 89 <#_Toc235985658>`__

`Addresses and Contacts. 90 <#_Toc235985659>`__

`Contributors 90 <#_Toc235985660>`__

`Addresses 90 <#_Toc235985661>`__

`References. 92 <#_Toc235985662>`__

`T-Coffee 92 <#_Toc235985663>`__

`Mocca. 93 <#_Toc235985664>`__

`CORE. 94 <#_Toc235985665>`__

`Other Contributions 94 <#_Toc235985666>`__

`Bug Reports and Feedback. 94 <#_Toc235985667>`__

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 3.0pt;
   padding:1.0pt 4.0pt 1.0pt 4.0pt;background:#CCCCCC;margin-left:0cm;margin-right:
   4.25pt">

 Cheat Sheet: T-Coffee

.. raw:: html

   </div>

.. rubric:: Proteins

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6">

Mode                 Command

============================================================================

Very Fast            t\_coffee sample\_aln1.fasta -mode quickaln

                     lower -ndiag if the sequences are very similar

----------------------------------------------------------------------------

Regular              t\_coffee sample\_aln1.fasta

                     use the output.html to estimate the MSA accuracy

----------------------------------------------------------------------------

Very Accurate       t\_coffee sample\_aln1.fasta -mode accurate

                     slow, combines structures, sequences and profiles

----------------------------------------------------------------------------

M-Coffee             t\_coffee sample\_aln1.fasta -mode mcoffee

                     combines most of the existing MSA packages

----------------------------------------------------------------------------

3D-Coffee            t\_coffee sample\_aln1.fasta -mode 3dcoffee

                     uses the structure of your sequences if named with
PDBID

----------------------------------------------------------------------------

Expresso             t\_coffee sample\_aln1.fasta -mode expresso

                     finds structures homologous to your sequences

----------------------------------------------------------------------------

PSI-Coffee           t\_coffee sample\_aln1.fasta -mode psicoffee

                     enriches your sequence with profile information

----------------------------------------------------------------------------

.. raw:: html

   </div>

.. rubric:: DNA

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6">

R-Coffee             t\_coffee three\_cdna.fasta -mode cdna

.. raw:: html

   </div>

.. rubric:: RNA

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6">

Mode                 Command

============================================================================

R-Coffee             t\_coffee sample\_rnaseq1.fasta -mode rcoffee

                     use the predicted secondary structure of your
sequences

----------------------------------------------------------------------------

RM-Coffee            t\_coffee sample\_rnaseq1.fasta -mode rmcoffee

                     use M-Coffee + secondary structure prediction

----------------------------------------------------------------------------

R-Coffee Consan      t\_coffee sample\_rnaseq1.fasta -mode
rcoffee\_consan

                     use rcoffee to combine consan alignments. Accurate
and Slow

 

.. raw:: html

   </div>

.. rubric:: Memory Problems

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6">

memory               t\_coffee sample\_aln1.fasta -mode memory

.. raw:: html

   </div>

 

 

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 3.0pt;
   padding:1.0pt 4.0pt 1.0pt 4.0pt;background:#CCCCCC;margin-left:0cm;margin-right:
   4.25pt">

Before You Start…

.. raw:: html

   </div>

.. rubric:: Foreword

A lot of the stuff presented here emanates form two summer schools that
were tentatively called the "Prosite Workshops" and were held in
Marseille, in 2001 and 2002. These workshops were mostly an excuse to go
rambling and swimming in the callanques. Yet, when we got tired of
lazing in the sun, we eventually did a bit of work to chill out. Most of
our experiments were revolving around the development of sequence
analysis tools. Many of the most advanced ideas in T-Coffee were
launched during these fruitful sessions. Participants included Phillip
Bucher, Laurent Falquet, Marco Pagni, Alexandre Gattiker, Nicolas Hulo,
Christian Siegfried, Anne-Lise Veuthey, Virginie Leseau, Lorenzo Ceruti
and Cedric Notredame.

This Document contains two main sections. The first one is a tutorial,
where we go from simple things to more complicated and show you how to
use all the subtleties of T-Coffee. We have tried to put as many of
these functionalities on the web (www.tcoffee.org) but if you need to do
something special and highly reproducible, the Command Line is the only
way.  

.. rubric:: Pre-Requisite

This tutorial relies on the assumption that you have installed T-Coffee,
version 6.18 or higher.

T-Coffee is a freeware open source running on all Unix-like platforms,
including MAC-osX and Cygwin. All the relevant information for
installing T-Coffee is contained in the Technical Documentation
(tcoffee\_technical.doc in the doc directory.)

T-Coffee cannot run on the Microsoft Windows shell. If you need to run T
-Coffee on windows, start by installing cygwin (www.cygwin.com). Cygwin
is a freeware open source that makes it possible to run a unix-like
command line on your Microsoft Windows PC without having to reboot.
Cygwin is free of charge and very easy to install. Yet, as the first
installation requires downloading substantial amounts of data, you
should make sure you have access to a broad-band connection.

In the course of this tutorial, we expect you to use a unix-like command
line shell. If you work on Cygwin, this means clicking on the cygwin
icon and typing commands in the window that appears. If you don't want
to bother with command line stuff, try using the online tcoffee
webserver at: **www.tcoffee.org**

.. rubric:: Getting the Example Files of the Tutorial

We encourage you to try all the following examples with your own
sequences/structures. If you want to try with ours, you can get the
material from the example directory of the distribution. If you do not
know where this file leaves or if you do not have access to it, the
simplest thing to do is to:

1-    download T-Coffee's latest version from www.tcoffee.org (Follow
the link to the T-Coffee Home Page)

2-    Download the latest distribution

3-    gunzip <distrib>.tar.gz

4-    tar -xvf <distrib>.tar

5-    go into <distrib>/example

This is all you need to do to run ALL the examples provided in this
tutorial.

 

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 3.0pt;
   padding:1.0pt 4.0pt 1.0pt 4.0pt;background:#CCCCCC;margin-left:0cm;margin-right:
   4.25pt">

What Is
 T-COFFEE
 ?

.. raw:: html

   </div>

.. rubric:: What is T-Coffee?

Before going deep into the core of the matter, here are a few words to
quickly explain some of the things T-Coffee will do for you.

.. rubric:: What does it do?

T-Coffee is a multiple sequence alignment program: given a set of
sequences previously gathered using database search programs like BLAST,
FASTA or Smith and Waterman, T-Coffee will produce a multiple sequence
alignment. **To use T-Coffee you must already have your sequences
ready**.

T-Coffee can also be used to compare alignments, reformat them or
evaluate them using structural information, it is the mode known as
seq\_reformat.

.. rubric:: What can it align?

T-Coffee will align nucleic and protein sequences alike. It will be able
to use structural information for protein sequences with a known
structure or the RNA sequences. On a new PC mid of the range, T-Coffee
will align up to a 100 sequences, about 1000 amino acid long.

.. rubric:: How can I use it?

T-Coffee is not an interactive program. It runs from your UNIX or Linux
command line and you must provide it with the correct parameters. If you
do not like typing commands, here is the simplest available mode where
T-Coffee only needs the name of the sequence file:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

          PROMPT: t\_coffee sample\_seq1.fasta

.. raw:: html

   </div>

Installing and using T-Coffee requires a minimum acquaintance with the
Linux/Unix operating system. If you feel this is beyond your computer
skills, we suggest you use one of the available online servers.

.. rubric:: Is There an Online Server

Yes, at www.tcoffee.org

.. rubric:: Is T-Coffee different from ClustalW?

According to several benchmarks, T-Coffee appears to be more accurate
than ClustalW. Yet, this increased accuracy comes at a price: T-Coffee
is slower than Clustal (about N times fro N Sequences).

If you are familiar with ClustalW, or if you run a ClustalW server, you
will find that we have made some efforts to ensure as much compatibility
as possible between ClustalW and T-COFFEE. Whenever it was relevant, we
have kept the flag names and the flag syntax of ClustalW. Yet, you will
find that T-Coffee also has many extra possibilities…

If you want to align closely related sequences, T-Coffee can also be
used in a fast mode, much faster than ClustalW, and about as accurate:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

          PROMPT: t\_coffee sample\_seq1.fasta -mode quickaln

.. raw:: html

   </div>

This mode works by only considering the best diagonals between two
sequences. By default all the diagonals having a substitution score >0
are considered, but you can lower this by specifying:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

          PROMPT: t\_coffee sample\_seq1.fasta -mode quickaln -ndiag=10

.. raw:: html

   </div>

That will only consider the top 10 diagonals. This will be very useful
if you have long and very similar sequences to align (DNA for instance).

.. rubric:: Is T-Coffee very accurate?

T-Coffee combines methods, and can be made as accurate (and hopefully
more) as the methods it combines. If you need a very accurate alignment
(and you have the full package installed with SOAP);

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

   PROMPT: t\_coffee sample\_seq1.fasta -mode accurate

.. raw:: html

   </div>

If you cannot run this job, go to the first section of the technical
manual (Installing BLAST for T-Coffee). You don't necessary need to
install BLAST locally but you must have access to a remote server (EBI
or NCBI).

This mode is very slow but also very accurate. On average this mode is
about 10 % more accurate than normal aligners on sequences less than 30%
similar. If you want something faster:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

t\_coffee sample\_seq1.fasta

.. raw:: html

   </div>

This is the normal mode. It is one of the most accurate of its kind,
roughly like Probcons.

 

 

        

.. rubric:: What T-Coffee Can and Cannot do for you …

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.5pt;
   padding:0cm 0cm 0cm 0cm;background:#F7CBB7">

IMPORTANT: All the files mentioned here (sample\_seq...) can be found in
the example directory of the distribution.

.. raw:: html

   </div>

.. rubric:: (NOT) Fetching Sequences

T-Coffee will NOT fetch sequences for you: you must select the sequences
you want to align before hand. We suggest you use any BLAST server and
format your sequences in FASTA so that T-COFFEE can use them easily. The
expasy BLAST server (www.expasy.ch) provides a nice interface for
integrating database searches.

Yet, the new modes of

.. rubric:: Aligning Sequences

T-Coffee will compute (or at least try to compute!) accurate multiple
alignments of DNA, RNA or Protein sequences.

.. rubric:: Combining Alignments

T-Coffee allows you to combine results obtained with several alignment
methods. For instance if you have an alignment coming from ClustalW, an
other alignment coming from Dialign, and a structural alignment of some
of your sequences, T-Coffee will combine all that information and
produce a new multiple sequence alignment having the best agreement with
all these methods (see the FAQ for more details)

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee –aln=sproteases\_small.cw\_aln,
sproteases\_small.muscle, sproteases\_small.tc\_aln
–outfile=combined\_aln.aln

.. raw:: html

   </div>

.. rubric:: Evaluating Alignments

You can use T-Coffee to measure the reliability of your Multiple
Sequence alignment. If you want to find out about that, read the FAQ or
the documentation for the -output flag.

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee –infile=sproteases\_small.aln –special\_mode=evaluate

.. raw:: html

   </div>

.. rubric:: Combining Sequences and Structures

One of the latest improvements of T-Coffee is to let you combine
sequences and structures, so that your alignments are of higher quality.
You need to have sap package installed to fully benefit of this
facility. If you have the EBI BLAST client installed (see installation
procedure), you can run the following:

 

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee 3d.fasta –special\_mode=expresso

.. raw:: html

   </div>

BLAST will identify the best PDB target for each sequences, and T-Coffee
will use sap (or any other structural package) to align your structures
and your sequences. If you do not have BLAST installed, or if you want
to specify the templates yourself, you can use 3D-Coffee:

 

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee 3d.fasta –special\_mode=3dcoffee

.. raw:: html

   </div>

In this case, the sequences must be names according to their PDB
targets.  All these network based operations are carried out using wget.
If wget is not installed on your system, you can get it for free from
(www.wget.org). To make sure wget is installed on your system, type

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: which wget

.. raw:: html

   </div>

.. rubric:: Identifying Occurrences of a Motif: Mocca

Mocca is a special mode of T-Coffee that allows you to extract a series
of repeats from a single sequence or a set of sequences. In other words,
if you know the coordinates of one copy of a repeat, you can extract all
the other occurrences. If you want to use Mocca, simply type:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee –other\_pg mocca sample\_seq1.fasta

.. raw:: html

   </div>

The program needs some time to compute a library and it will then prompt
you with an interactive menu. Follow the instructions.

.. rubric:: How Does T-Coffee works

If you only want to make a standard multiple alignments, you may skip
these explanations. But if you want to do more sophisticated things,
these few indications may help before you start reading the doc and the
papers.

When you run T-Coffee, the first thing it does is to compute a library.
The library is a list of pairs of residues that *could* be aligned. It
is like a Xmas list: you can ask anything you fancy, but it is down to
Santa to assemble a collection of Toys that won't get him stuck at the
airport, while going through the metal detector.

Given a standard library, it is not possible to have all the residues
aligned at the same time because all the lines of the library may not
agree. For instance, line 1 may say

Residue 1 of seq A with Residue 5 of seq B,

and line 100 may say

Residue 1 of seq A with Residue 29 of seq B,

Each of these constraints comes with a weight and in the end, the
T-Coffee algorithm tries to generate the multiple alignment that
contains constraints whose sum of weights yields the highest score. In
other words, it tries to make happy as many constraints as possible
(replace the word constraint with, friends, family members,
collaborators… and you will know exactly what we mean).

You can generate this list of constraints however you like. You may even
provide it yourself, forcing important residues to be aligned by giving
them high weights (see the FAQ). For your convenience, T-Coffee can
generate (this is the default) its own list by making all the possible
global pairwise alignments, and the 10 best local alignments associated
with each pair of sequences. Each pair of residues observed aligned in
these pairwise alignments becomes a line in the library.

Yet be aware that nothing forces you to use this library and that you
could build it using other methods (see the FAQ). In protein language,
T-COFEE is synonymous for freedom, the freedom of being aligned however
you fancy ( I was a Tryptophan in some previous life).

 

 

 

 

           

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 3.0pt;
   padding:1.0pt 4.0pt 1.0pt 4.0pt;background:#CCCCCC;margin-left:0cm;margin-right:
   4.25pt">

Preparing Your Data:
 Reformatting and Trimming With seq\_reformat

.. raw:: html

   </div>

Nothing is more frustrating than downloading important data and
realizing you need to format it \*before\* using it. In general, you
should avoid manual reformatting: it is by essence inconsistent and will
get you into trouble. It will also get you depressed when you will
realize that you have spend the whole day adding carriage return to each
line in your files.

.. rubric:: Seq\_reformat

.. rubric:: Accessing the T-Coffee Reformatting Utility

T-Coffee comes along with a very powerful reformatting utility named
seq\_reformat. You can use seq\_reformat by invoking the t\_coffee
shell.

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat

.. raw:: html

   </div>

This will output the online flag usage of seq\_reformat. Seq\_reformat
recognizes automatically the most common formats. You can use it to:

Reformat your sequences.

extract sub-portions of alignments

Extract sequences.

In this section we give you a few examples of things you can do with
seq\_reformat:

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.5pt;
   padding:0cm 0cm 0cm 0cm;background:#F7CBB7">

Warning: after the flag -other\_pg, the T-Coffee flags are not any more
recognized. It is like if you were using a different programme

.. raw:: html

   </div>

 

.. rubric:: An overview of seq\_reformat

seq\_reformat is a reformatting utility. It reads in via the -in and
-in2 flags and outputs in whatever specified format via the -output
flag. In the meantime, you can use the flag '-action' to modify your
data, using any of the flag. If you want a complete list of things
seq\_reformat can do for you, try:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat

.. raw:: html

   </div>

.. rubric:: Reformatting your data

           

.. rubric:: Changing MSA formats

It can be necessary to change from one MSA format to another. If your
sequences are in ClustalW format and you want to turn them into fasta,
while keeping the gaps, try

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sproteases\_small.aln
-output fasta\_aln > sproteases\_small.fasta\_aln

.. raw:: html

   </div>

If you want to turn a clustalw alignment into an alignment having the
pileup format (MSF), try:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sproteases\_small.aln
-output msf > sproteases\_small.msf

.. raw:: html

   </div>

.. rubric:: Dealing with Non-automatically recognized formats

    Format recognition is not 100% full proof. Occasionally you will
have to inform the program about the nature of the file you are trying
to reformat:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

    -in\_f msf\_aln for intance

.. raw:: html

   </div>

.. rubric:: Automated Sequence Edition

.. rubric:: Removing the gaps from an alignment

If you want to recover your sequences from some pre-computed alignment,
you can try:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sproteases\_small.aln
-output fasta\_seq > sproteases\_small.fasta

.. raw:: html

   </div>

This will remove all the gaps.

.. rubric:: Changing the case of your sequences

If you need to change the case of your sequences, you can use more
sophisticated functions embedded in seq\_reformat. We call these
modifiers, and they are accessed via the -action flag. For instance, to
write our sequences in lower case:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sproteases\_small.aln
-action +lower -output clustalw

.. raw:: html

   </div>

No prize for guessing that +upper will do exactly the opposite....

.. rubric:: Changing the case of specific residues

If you want to change the case of a specific residue, you can use the
flag: +edit\_residue <sequence> <residue #> <lower\|upper\|symbol>. If
you have more than one residue to color, you can put all the coordinates
in a file, (one coordinate per line). Spans are not yet supported.

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sample\_aln1.aln -action
+upper +edit\_residue hmgb\_chite 10 lower

.. raw:: html

   </div>

.. rubric:: Changing the case depending on the score

If you want to change the case depending on the score, you must either
evaluate your alignment, or provide cache (see next section for the
cache). If you want to evaluate on the fly, try:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sample\_aln1.aln -in3
sample\_aln1.aln -action +upper +3evaluate idmat +lower '[5-9]'

.. raw:: html

   </div>

Will lower the case of every residue identical to more than 50% of the
residue in its column.

.. rubric:: Protecting Important Sequence Names

Few programs support long sequence names. Sometimes, when going through
some pipeline the names of your sequences can be damaged (truncated or
modified). To avoid this, seq\_reformat contains a utility that can
automatically rename your sequences into a form that will be machine
friendly, while making it easy to return to the human friendly form.

The first thing to do is to generate a list of names that will be used
in place of the long original name of the sequences. For instance:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sproteases\_large.fasta
-output code\_name > sproteases\_large.code\_name

.. raw:: html

   </div>

Will create a file where each original name is associated with a coded
name (Cxxxx). You can then use this file to either code or decode your
dataset. For instance, the following command:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -code
sproteases\_large.code\_name -in sproteases\_large.fasta
>sproteases\_large.coded.fasta

.. raw:: html

   </div>

Will code all the names of the original data. You can work with the file
sproteases\_large.coded.fasta, and when you are done, you can de-code
the names of your sequences using:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -decode
sproteases\_large.code\_name -in sproteases\_large.coded.fasta

.. raw:: html

   </div>

.. rubric:: Colouring/Editing Residues in an Alignment

.. rubric:: Coloring specific types of residues

You can color all the residues of your sequences on the fly. For
instance, the following command:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sample\_aln1.aln -in3
sample\_aln1.aln -action  +3convert a0 -output color\_html >
colored.html

.. raw:: html

   </div>

will color all the As in color 0 (blue).

.. rubric:: Coloring a specific residue of a specific sequence

If you want to color a specific residue, you can use the flag:
+color\_residue <sequence> <residue #> <color #>. If you have more than
one residue to color, you can put all the coordinates in a file, (one
coordinate per line). Spans are not yet supported.

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sample\_aln1.aln -action
+color\_residue hmgb\_chite 10 1 -output color\_html > color.html

.. raw:: html

   </div>

.. rubric:: Coloring according to the conservation

Use the +evaluate flag if you want to color your alignment according to
its conservation level

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sample\_aln1.aln -in3
sample\_aln1.aln -action +3evaluate pam250mt- output color\_html >
color.html

.. raw:: html

   </div>

You can also use the boxshade scoreing scheme:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sample\_aln1.aln -in3
sample\_aln1.aln -action +3evaluate boxshade -output color\_html >
color.html

.. raw:: html

   </div>

 

 

.. rubric:: Colouring/Editing residues in an Alignment Using a Cache

.. rubric:: Overview

To color an alignment, two files are needed: the alignment (aln) and the
cache (cache). The cache is a file where residues to be colored are
declared along with the colors. Nine different colors are currently
supported. They are set by default but can be modified by the user (see
last changing default colors). The cache can either look like a standard
sequence or alignment file (see below) or like a standard T-Coffee
library (see next section). In this section we show you how to
specifically modify your original sequences to turn them into a cache.

In the cache, the colors of each residue are declared with a number
between 0 and 9.  Undeclared residues will appear without any color in
the final alignment.

.. rubric:: Preparing a Sequence or Alignment Cache

Let us consider the following file:

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6;margin-left:0cm;margin-right:1.0cm">

CLUSTAL FORMAT

 

B               CTGAGA-AGCCGC---CTGAGG--TCG

C               TTAAGG-TCCAGA---TTGCGG--AGC

D               CTTCGT-AGTCGT---TTAAGA--ca-

A               CTCCGTgTCTAGGagtTTACGTggAGT

                 \*  \*      \*     \*  \*     

.. raw:: html

   </div>

The command

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in=sample\_aln6.aln
-output=clustalw\_aln -out=cache.aln -action +convert 'Aa1' '.--'
+convert '#0'

.. raw:: html

   </div>

The conversion will proceed as follows:

-conv indicates the filters for character conversion:

                     - will remain -

                     A and a will be turned into 1

                     All the other symbols (#) will be turned into 0.

-action +convert, indicates the actions that must be carried out on the
alignment before it is output into cache.

This command generates the following alignment (called a cache):

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6;margin-left:0cm;margin-right:1.0cm">

CLUSTAL FORMAT for SEQ\_REFORMAT Version 1.00, CPU=0.00 sec, SCORE=0,
Nseq=4, Len=27

 

B               000101-100000---000100--000

C               001100-000101---000000--100

D               000000-100000---001101--01-

A               000000000010010000100000100

.. raw:: html

   </div>

Other alternative are possible. For instance, the following command:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in=sample\_aln6.aln
-output=fasta\_seq -out=cache.seq -action +convert 'Aa1' '.--' +convert
'#0'

.. raw:: html

   </div>

will produce the following file cache\_seq

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6;margin-left:0cm;margin-right:1.0cm">

>B

000101100000000100000

>C

001100000101000000100

>D

00000010000000110101

>A

000000000010010000100000100

.. raw:: html

   </div>

where each residue has been replaced with a number according to what was
specified by conv. Note that it is not necessary to replace EVERY
residue with a code. For instance, the following file would also be
suitable as a cache:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in=sample\_aln6.aln
-output=fasta\_seq -out=cache -action +convert 'Aa1' '.--'

.. raw:: html

   </div>

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6;margin-left:0cm;margin-right:1.0cm">

>B

CTG1G11GCCGCCTG1GGTCG

>C

TT11GGTCC1G1TTGCGG1GC

>D

CTTCGT1GTCGTTT11G1c1

>A

CTCCGTgTCT1GG1gtTT1CGTgg1GT

.. raw:: html

   </div>

.. rubric:: Preparing a Library Cache

The Library is a special format used by T-Coffee to declare special
relationships between pairs of residues. The cache library format can
also be used to declare the color of specific residues in an alignment.
For instance, the following file

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6;margin-left:0cm;margin-right:1.0cm">

! TC\_LIB\_FORMAT\_01

4

A 27 CTCCGTgTCTAGGagtTTACGTggAGT

B 21 CTGAGAAGCCGCCTGAGGTCG

C 21 TTAAGGTCCAGATTGCGGAGC

D 20 CTTCGTAGTCGTTTAAGAca

#1 1

    1     1   3

    4     4   5

#3 3

    6     6   1

    9     9   4

! CPU 240

! SEQ\_1\_TO\_N

.. raw:: html

   </div>

sample\_lib5.tc\_lib declares that residue 1 of sequence 3 will be
receive color 6, while residue 20 of sequence 4 will receive color 20.
Note that the sequence number and the residue index are duplicated,
owing to the recycling of this format from its original usage.

It is also possible to use the BLOCK operator when defining the library
(c.f. technical doc, library format). For instance:

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6;margin-left:0cm;margin-right:1.0cm">

! TC\_LIB\_FORMAT\_01

4

A 27 CTCCGTgTCTAGGagtTTACGTggAGT

B 21 CTGAGAAGCCGCCTGAGGTCG

C 21 TTAAGGTCCAGATTGCGGAGC

D 20 CTTCGTAGTCGTTTAAGAca

#1 1

    +BLOCK+ 10 1     1    3

    +BLOCK+ 5  15    15   5

#3 3

    6     6   1

    9     9   4

! CPU 240

! SEQ\_1\_TO\_N

.. raw:: html

   </div>

The number right after BLOCK indicates the block length (10). The two
next numbers (1 1) indicate the position of the first element in the
block. The last value is the color.

.. rubric:: Coloring an Alignment using a cache

If you have a cache alignment or a cache library, you can use it to
color your alignment and either make a post script, html or PDF output.
For instance, if you use the file cache.seq:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

   PROMPT: t\_coffee -other\_pg seq\_reformat -in=sample\_aln6.aln
-struc\_in=sample\_aln6.cache -struc\_in\_f number\_fasta
-output=color\_html -out=x.html

.. raw:: html

   </div>

This will produce a colored version readable with any standard web
browser, while:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

   PROMPT: t\_coffee -other\_pg seq\_reformat -in=sample\_aln6.aln
-struc\_in=sample\_aln6.cache -struc\_in\_f number\_fasta
-output=color\_pdf -out=x.pdf

.. raw:: html

   </div>

This will produce a colored version readable with acrobat reader.

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.5pt;
   padding:0cm 0cm 0cm 0cm;background:#F7CBB7">

Warning: ps2pdf must be installed on your system

.. raw:: html

   </div>

You can also use a cache library like the one shown above
(sample\_lib5.tc\_lib):

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in=sample\_aln6.aln
-struc\_in=sample\_lib5.tc\_lib -output=color\_html -out=x.html

.. raw:: html

   </div>

.. rubric:: Changing the default colors

Colors are hard coded in the program, but if you wish, you can change
them, simply create a file named:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

   seq\_reformat.color

.. raw:: html

   </div>

That is used to declare the color values:

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6;margin-left:0cm;margin-right:1.0cm">

 

0 #FFAA00 1 0.2 0

 

.. raw:: html

   </div>

This indicates that the value 0 in the cache corresponds now to #FFAA00
in html, and in RGB 1,  0.2 and 0. The name of the file
(seq\_reformat.color) is defined in: programmes\_define.h, COLOR\_FILE.
And can be changed before compilation. **By default, the file is
searched in the current directory**

.. rubric:: Evaluating an alignment and producing a cache

.. rubric:: Evaluating an alignment with T-Coffee

As suggested in a previous section, it is possible to evaluate the
accuracy of any alignment using a T-Coffee library. The simplest way to
do that is to compute a default library and evaluate the target
alignment against this library:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -infile sample\_aln1.aln -mode evaluate

.. raw:: html

   </div>

This command will output a file named sample\_aln1.score\_asccii that
can then be used to either evaluate the local accuracy of the alignment
or automatically filter it using the seq\_reformat utility.

In some circumstances, you may also want to evaluate your alignment
against a pre-computed library. This can be easily achieved:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -infile sample\_aln1.aln -out\_lib
sample\_aln1.tc\_lib -lib\_only

.. raw:: html

   </div>

 

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -infile sample\_aln1.aln -mode evaluate -lib
sample\_aln1.tc\_lib

.. raw:: html

   </div>

When using this last command, the reference library will be the one
provided by the user. The local score thus reported is the CORE index.

.. rubric:: Evaluating the level of conservation with a substitution
   matrix

It is possible to use seq\_reformat in a similar way to infer the local
level of identity, either using an identity matrix or with any regular
matrix, in which case, every residue with a substitution score higher
than 0 is counted as an identity. This can be achieved as follows for
identity measure:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sample\_aln1.aln -action
+evaluate idmat -output score\_ascii

.. raw:: html

   </div>

Or with the following for measuring similarity with a blosum62

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sample\_aln1.aln -action
+evaluate blosum62mt -output score\_ascii

.. raw:: html

   </div>

Finally, it is also possible to display in color the conservation
levels:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

   PROMPT: t\_coffee -other\_pg seq\_reformat -in sample\_aln1.aln
-action +evaluate blosum62mt -output score\_html > x.html

.. raw:: html

   </div>

 

 

.. rubric:: Selective Reformatting

 

.. rubric:: Removing gapped columns

You can remove all the columns containing a certain proportion of gaps.
For instance:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sample\_aln7.aln -action
+rm\_gap 50

.. raw:: html

   </div>

Will delete all the residues occurring in a column that contains 50% or
more gaps (use 1 to delete residues from columns having 1 gap or more).

 

.. rubric:: 

 

.. rubric:: Selectively turn some residues to lower case

Consider the following alignment (sample\_aln7.aln)

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6;margin-left:0cm;margin-right:1.0cm">

CLUSTAL FORMAT for T-COFFEE Version\_4.62 [http://www.tcoffee.org],
CPU=0.04 sec, SCORE=0, Nseq=4, Len=28

 

A               CTCCGTGTCTAGGAGT-TTACGTGGAGT

B               CTGAGA----AGCCGCCTGAGGTCG---

D               CTTCGT----AGTCGT-TTAAGACA---

C               -TTAAGGTCC---AGATTGCGGAGC---

                 \* ..        .\*  \* . \*:

.. raw:: html

   </div>

and the following cache (sample\_aln7.cache\_aln):

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6;margin-left:0cm;margin-right:1.0cm">

CLUSTAL FORMAT for T-COFFEE Version\_4.62 [http://www.tcoffee.org],
CPU=0.04 sec, SCORE=0, Nseq=4, Len=28

 

A               3133212131022021-11032122021

B               312020----023323312022132---

D               311321----021321-11002030---

C               -110022133---020112322023---

.. raw:: html

   </div>

 

You can turn to lower case all the residues having a score between 1 and
2:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sample\_aln7.aln
-struc\_in sample\_aln7.cache\_aln -struc\_in\_f number\_aln -action
+lower '[1-2]'

.. raw:: html

   </div>

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6;margin-left:0cm;margin-right:1.0cm">

CLUSTAL FORMAT for T-COFFEE Version\_4.62 [http://www.tcoffee.org],
CPU=0.05 sec, SCORE=0, Nseq=4, Len=28

 

A               CtCCgtgtCtAggAgt-ttACgtggAgt

B               CtgAgA----AgCCgCCtgAggtCg---

D               CttCgt----AgtCgt-ttAAgACA---

C               -ttAAggtCC---AgAttgCggAgC---

                 \* ..        .\*  \* . \*:

.. raw:: html

   </div>

 

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.5pt;
   padding:1.0pt 4.0pt 1.0pt 4.0pt;background:#FFFF89">

Note: that residues not concerned will keep their original case

.. raw:: html

   </div>

.. rubric:: Selectively modifying residues

The range operator is supported by three other important modifiers:

         \ **-upper:** to uppercase your residues

         \ **-lower:** to lowercase your residues

         \ **-switchcase:** to selectively toggle the case of your
residues

         \ **-keep:** to only keep the residues within the range

         \ **-remove:** to remove the residues within the range

         \ **-convert:** to only convert the residues within the range.

For instance, to selectively turn all the G having a score between 1 and
2, use:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sample\_aln7.aln
-struc\_in sample\_aln7.cache\_aln -struc\_in\_f number\_aln -action
+convert '[1-2]' CX

.. raw:: html

   </div>

.. rubric:: Keeping only the best portion of an alignment

To do this, you need an evaluation file that may have been generated
with T-Coffee, either running a de-novo alignment

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sample\_seq1.fasta -output score\_ascii, aln

.. raw:: html

   </div>

Or evaluating a pre-existing alignment

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sample\_seq1.aln -action
+evaluate blosum62mt -output score\_ascii

.. raw:: html

   </div>

This generates a score\_ascii file that you can then use to filter out
the bad bits in your alignement:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sample\_seq1.aln
-struc\_in sample\_seq1.score\_ascii -struc\_in\_f number\_aln  -action
+keep '[8-9]'

.. raw:: html

   </div>

This command considers the individual score of each residue to trigger
the filtering. It is also possible to do this according to the whole
column. Simply add the "+use\_cons" flag.

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sample\_seq1.aln
-struc\_in sample\_seq1.score\_ascii -struc\_in\_f number\_aln  -action
+use\_cons +keep '[8-9]'

.. raw:: html

   </div>

 

.. rubric:: 

 

.. rubric:: Extracting Portions of Dataset

Extracting portions of a dataset is something very frequently needed.
You may need to extract all the sequences that contain the word human in
their name, or you may want all the sequences containing a simple motif.
We show you here how to do a couple of these things.

.. rubric:: Extracting The High Scoring Blocks

It is possible to use a score\_ascii file ( as produced in the previous
section) in order to extract high scoring portions of an alignment. For
instance, the following command:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

   PROMPT: t\_coffee -other\_pg seq\_reformat -in sample\_aln1.aln
-action +evaluate blosum62mt +use\_cons +keep '[5-9]'

.. raw:: html

   </div>

will keep all the residues having a column conservation score between 5
and 9

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.5pt;
   padding:1.0pt 4.0pt 1.0pt 4.0pt;background:#FFFF89">

Note: Don't forget the simple quotes! (')

.. raw:: html

   </div>

It is also possible to re-use pre-computed score\_ascii files, such as
those obtained when computing a T-Coffee multiple alignment. For
instance, the following series of command will make it possible to
extract the positions having a consistency score between 6 and 9:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

   PROMPT: t\_coffee sample\_aln1.fasta -output score\_ascii -outfile
sample1.score\_ascii

.. raw:: html

   </div>

 

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

   PROMPT: t\_coffee -other\_pg seq\_reformat -in sample\_aln1.aln
-struc\_in sample1.score\_ascii -struc\_in\_f number\_aln -action
+use\_cons +keep '[8-9]'

.. raw:: html

   </div>

 

.. rubric:: Extracting Sequences According to a Pattern

You can extract any sequence by requesting a specific pattern to be
found either in the name, the comment or the sequence. For instance, if
you want to extract all the sequences whose name contain the word HUMAN:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sproteases\_small.aln
-action +grep NAME KEEP HUMAN -output clustalw

.. raw:: html

   </div>

The modifier is "+grep". NAME indicates that the extraction is made
according to the sequences names, and KEEP means that you will keep all
the sequences containing the string HUMAN. If you wanted to remove all
the sequences whose name contains the word HUMAN, you should have typed:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sproteases\_small.aln
-action +grep NAME REMOVE HUMAN -output clustalw

.. raw:: html

   </div>

Note that  HUMAN is case sensitive (Human, HUMAN and hUman will not
yield the same results). You can also select the sequences according to
some pattern found in their COMMENT section or directly in the sequence.
For instance

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sproteases\_small.aln
-action +grep COMMENT KEEP sapiens -output clustalw

.. raw:: html

   </div>

Will keep all the sequences containing the word ***sapiens*** in the
comment section. Last but not least, you should know that the pattern
can be any perl legal regular expression (See
www.comp.leeds.ac.uk/Perl/matching.html for some background on regular
expressions). For instance:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sproteases\_small.aln
-action +grep NAME REMOVE '[ILM]K' -output clustalw

.. raw:: html

   </div>

Will extract all the sequences containing the pattern [ILM]K.

.. rubric:: Extracting Sequences by Names

***Extracting Two Sequences***: If you want to extract several
sequences, in order to make a subset. You can do the following:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sproteases\_small.aln
-action +extract\_seq\_list 'sp\|P29786\|TRY3\_AEDAE'
'sp\|P35037\|TRY3\_ANOGA'

.. raw:: html

   </div>

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.5pt;
   padding:0cm 0cm 0cm 0cm;background:#F7CBB7">

Note the single quotes ('). They are meant to protect the name of your
sequence and prevent the UNIX shell to interpret it like an instruction.

.. raw:: html

   </div>

**Removing Columns of Gaps.** Removing intermediate sequences results in
columns of gaps appearing here and there. Keeping them is convenient if
some features are mapped on your alignment. On the other hand, if you
want to remove these columns you can use:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sproteases\_small.aln
-action +extract\_seq\_list 'sp\|P29786\|TRY3\_AEDAE'
'sp\|P35037\|TRY3\_ANOGA' +rm\_gap

.. raw:: html

   </div>

***Extracting Sub sequences***: You may want to extract portions of your
sequences. This is possible if you specify the coordinates after the
sequences name:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sproteases\_small.aln
-action +extract\_seq 'sp\|P29786\|TRY3\_AEDAE' 20 200
'sp\|P35037\|TRY3\_ANOGA' 10 150 +rm\_gap

.. raw:: html

   </div>

***Keeping the original Sequence Names.*** Note that your sequences are
now renamed according to the extraction coordinates. You can keep the
original names by using the +keep\_name modifier:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sproteases\_small.aln
-action +keep\_name +extract\_seq 'sp\|P29786\|TRY3\_AEDAE' 20 200
'sp\|P35037\|TRY3\_ANOGA' 10 150 +rm\_gap

.. raw:: html

   </div>

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.5pt;
   padding:0cm 0cm 0cm 0cm;background:#F7CBB7">

Note: +keep\_name must come BEFORE +extract\_seq

.. raw:: html

   </div>

.. rubric:: Removing Sequences by Names

***Removing Two Sequences***. If you want to remove several sequences,
use rm\_seq instead of keep\_seq.

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sproteases\_small.aln
-action +remove\_seq 'sp\|P29786\|TRY3\_AEDAE' 'sp\|P35037\|TRY3\_ANOGA'

.. raw:: html

   </div>

.. rubric:: Extracting Blocks Within Alignment

***Extracting a Block.*** If you only want to keep one block in your
alignment, use

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sproteases\_small.aln
-action +extract\_block cons 150 200

.. raw:: html

   </div>

In this command line, ***cons*** indicates that you are counting the
positions according to the consensus of the alignment (i.e. the
positions correspond to the columns # of the alignment). If you want to
extract your block relatively to a specific sequence, you should replace
cons with this sequence name. For instance:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sproteases\_small.aln
-action +extract\_block 'sp\|Q03238\|GRAM\_RAT' 10 200

.. raw:: html

   </div>

.. rubric:: Concatenating Alignments

If you have extracted several blocks and you now want to glue them
together, you can use the cat\_aln function

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sproteases\_small.aln
-action +extract\_block cons 100 120 > block1.aln

PROMPT: t\_coffee -other\_pg seq\_reformat -in sproteases\_small.aln
-action +extract\_block cons 150 200 > block2.aln

PROMPT: t\_coffee -other\_pg seq\_reformat -in block1.aln -in2
block2.aln -action +cat\_aln

 

.. raw:: html

   </div>

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.5pt;
   padding:0cm 0cm 0cm 0cm;background:#F7CBB7">

Note: The alignments do not need to have the same number of sequences
and the sequences do not need to come in the same order.

.. raw:: html

   </div>

.. rubric:: Analyzing your Multiple Sequence Alignment

.. rubric:: Estimating the diversity in your alignment

It is easy to measure the level of diversity within your multiple
sequence alignment. The following command:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sample\_aln1.aln -output
sim

.. raw:: html

   </div>

Will output all the pairwise identities, as well as the average level of
identity between each sequence and the others. You can sort and grep in
order to select the sequences you are interested in.

 

.. rubric:: Reducing and improving your dataset

Large datasets are problematic because they can be difficult to analyze.
The problem is that when there are too many sequences, MSA programs tend
to become very slow and inaccurate. Furthermore, you will find that
large datasets are difficult to display and analyze. In short, the best
size for an MSA dataset is between 20 and 40 sequences. This way you
have enough sequences to ***see*** the effect of evolution, but at the
same time the dataset is small enough so that you can visualize your
alignment and recompute it as many times as needed.

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.5pt;
   padding:0cm 0cm 0cm 0cm;background:#F7CBB7">

Note: If your sequence dataset is very large, seq\_reformat will compute
the similarity matrix between your sequences once only. It will then
keep it in its cache and re-use it any time you re-use that dataset. In
short this means that it will take much longer to run the first time.

.. raw:: html

   </div>

 

.. rubric:: Extracting the N most informative sequences

To be informative, a sequence must contain information the other
sequences do not contain. The N most informative sequences are the N
sequences that are as different as possible to one another, given the
initial dataset.

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sproteases\_large.fasta
-action +trim \_seq\_n10  -output fasta\_seq

.. raw:: html

   </div>

The arguments to trim include ***\_seq\_*** . It means your sequences
are provided unaligned. If your sequences are already aligned, you do
not need to provide this parameter. It is generaly more accurate to use
unaligned sequences.

The argument ***\_n10*** means you want to extract the 10 most
informative sequences. If you would rather extract the 20% most
informative sequences, use

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sproteases\_large.fasta
-action +trim \_seq\_N20  -output fasta\_seq

.. raw:: html

   </div>

 

.. rubric:: Extracting all the sequences less than X% identical

Removing the most similar sequences is often what people have in mind
when they talk about removing redundancy. You can do so using the trim
option. For instance, to generate a dataset where no pair of sequences
has more than 50% identity, use:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sproteases\_large.fasta
-action +trim \_seq\_%%50\_

.. raw:: html

   </div>

.. rubric:: 

 

.. rubric:: Speeding up the process

If you start form unaligned sequences, the removal of redundancy can be
slow. If your sequences have already been aligned using a fast method,
you can take advantage of this by replacing the \_seq\_ with \_aln\_

Note the difference of speed between these two command and the previous
one:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in kinases.aln -action +trim
\_aln\_%%50\_

.. raw:: html

   </div>

 

 

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

t\_coffee -other\_pg seq\_reformat -in kinases.fasta -action +trim
\_seq\_%%50\_

.. raw:: html

   </div>

Of course, using the MSA will mean that you rely on a more approximate
estimation of sequence similarity.

.. rubric:: Forcing Specific Sequences to be kept

Sometimes you want to trim while making sure specific important
sequences remain in your dataset. You can do so by providing trim with a
**string.** Trim will keep all the sequences whose name contains the
string. For instance, if you want to force trim to keep all the
sequences that contain the word HUMAN, no matter how similar they are to
one another, you can run the following command:

 

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sproteases\_large.fasta
-action +trim \_seq\_%%50 HUMAN

.. raw:: html

   </div>

When you give this command, the program will first make sure that all
the HUMAN sequences are kept and it will then assemble your 50% dataset
while keeping the HUMAN sequences. Note that string is a perl regular
expression.

By default, string causes all the sequences whose name it matches to be
kept. You can also make sure that sequences whose COMMENT or SEQUENCE
matches string are kept. For instance, the following line

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sproteases\_large.fasta
-action +trim \_seq\_%%50\_fCOMMENT '.apiens'

.. raw:: html

   </div>

Will cause all the sequences containing the regular expression '.apiens'
in the comment to be kept. The \_f symbol before COMMENT stands for
"\_field" If you want to  make a selection on the sequences:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sproteases\_large.fasta
-action +trim \_seq\_%%50\_fSEQ '[MLV][RK]'

.. raw:: html

   </div>

You can also specify the sequences you want to keep. To do so, give a
fasta file containing the name of these sequences via the -in2 file

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT:t\_coffee -other\_pg seq\_reformat -in sproteases\_large.fasta
-in2 sproteases\_small.fasta -action +trim \_seq\_%%40

.. raw:: html

   </div>

.. rubric:: Identifying and Removing Outlayers

Sequences that are too distantly related from the rest of the set will
sometimes have very negative effects on the overall alignment. To
prevent this, it is advisable not to use them. This can be done when
trimming the sequences. For instance,

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sproteases\_large.fasta 
-action +trim \_seq\_%%50\_O40

.. raw:: html

   </div>

The symbol \_O stands for Outlayers. It will lead to the removal of all
the sequences that have less than 40% average accuracy with all the
other sequences in the dataset.

.. rubric:: Chaining Important Sequences

In order to align two distantly related sequences, most multiple
sequence alignment packages perform better when provided with many
intermediate sequences that make it possible to "bridge" your two
sequences. The modifier ***+chain*** makes it possible to extract from a
dataset a subset of intermediate sequences that chain the sequences you
are interested in.

For instance, le us consider the two sequences:

 sp\|P21844\|MCPT5\_MOUSE      sp\|P29786\|TRY3\_AEDAE

These sequences have 26% identity. This is high enough to make a case
for a homology relationship between them, but this is too low to blindly
trust any pairwise alignment. With the names of the two sequences
written in the file sproteases\_pair.fasta, run the following command:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sproteases\_large.fasta
-in2 sproteases\_pair.fasta -action +chain > sproteases\_chain.fasta

.. raw:: html

   </div>

This will generate a dataset of 21 sequences, whith the following chain
of similarity between your two sequences:

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6;margin-left:0cm;margin-right:1.0cm">

N: 21 Lower: 40 Sim: 25 DELTA: 15

#sp\|P21844\|MCPT5\_MOUSE -->93 -->sp\|P50339\|MCPT3\_RAT -->85
-->sp\|P50341\|MCPT2\_MERUN -->72 -->sp\|P52195\|MCPT1\_PAPHA -->98
-->sp\|P56435\|MCPT1\_MACFA -->97 -->sp\|P23946\|MCPT1\_HUMAN -->8

1 -->sp\|P21842\|MCPT1\_CANFA -->77 -->sp\|P79204\|MCPT2\_SHEEP -->60
-->sp\|P21812\|MCPT4\_MOUSE -->90 -->sp\|P09650\|MCPT1\_RAT -->83
-->sp\|P50340\|MCPT1\_MERUN -->73 -->sp\|P11034\|MCPT1\_MOUSE

-->76 -->sp\|P00770\|MCPT2\_RAT -->71 -->sp\|P97592\|MCPT4\_RAT -->66
-->sp\|Q00356\|MCPTX\_MOUSE -->97 -->sp\|O35164\|MCPT9\_MOUSE -->61
-->sp\|P15119\|MCPT2\_MOUSE -->50 -->sp\|Q06606\|GRZ2\_RAT -

->54 -->sp\|P80931\|MCT1A\_SHEEP -->40 -->sp\|Q90629\|TRY3\_CHICK -->41
-->sp\|P29786\|TRY3\_AEDAE

.. raw:: html

   </div>

This is probably the best way to generate a high quality alignment of
your two sequences when using a progressive method like ClustalW,
T-Coffee, Muscle or Mafft.

.. rubric:: Manipulating DNA sequences

.. rubric:: Translating DNA sequences into Proteins

If your sequences are DNA coding sequences, it is always safer to align
them as proteins. Seq\_reformat makes it easy for you to translate your
sequences:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in
sproteases\_small\_dna.fasta -action +translate -output fasta\_seq

.. raw:: html

   </div>

.. rubric:: Back-Translation With the Bona-Fide DNA sequences

Once your sequences have been aligned, you may want to turn your protein
alignment back into a DNA alignment, either to do phylogeny, or maybe in
order to design PCR probes. To do so, use the following command:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in
sproteases\_small\_dna.fasta -in2 sproteases\_small.aln -action
+thread\_dna\_on\_prot\_aln -output clustalw

.. raw:: html

   </div>

.. rubric:: Finding the Bona-Fide Sequences for the Back-Translation

Use the online server Protogene, available from www.tcoffee.org.

.. rubric:: Guessing Your Back Translation

Back-translating means turning a protein sequence into a DNA sequence.
If you do not have the original DNA sequence, this operation will not be
exact, owing to the fact that the genetic code is degenerated. Yet, if a
random-back translation is fine with you, you can use the following
command.

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in
sproteases\_small\_dna.fasta -in2 sproteases\_small.aln -action
+thread\_dna\_on\_prot\_aln  -output clustalw

.. raw:: html

   </div>

In  this process, codons are chosen randomly. For instance, if an
amino-acid has four codons, the back-translation process will randomly
select one of these. If you need more sophisticated back-translations
that take into account the codon bias, we suggest you use more specific
tools like: alpha.dmi.unict.it/~ctnyu/bbocushelp.html

.. rubric:: Fetching a Structure

There are many reasons why you may need a structure. T-Coffee contains a
powerful utility named ***extract\_from\_pdb*** that makes it possible
to fetch the PDB coordinates of a structure or its FASTA sequence
without requiring a local installation.

By default, extract\_from\_pdb will start looking for the structure in
the current directory; it will then look it up locally (PDB\_DIR) and
eventually try to fetch it from the web (via a wget to www.rcsb.org).
All these settings can be customized using environment variables (see
the last section).

.. rubric:: Fetching a PDB structure

If you want to fetch the chain E of the PDB structure 1PPG, you can use:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg extract\_from\_pdb -infile 1PPGE

.. raw:: html

   </div>

.. rubric:: Fetching The Sequence of a PDB structure

To Fetch the sequence, use:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg extract\_from\_pdb -infile 1PPGE -fasta

.. raw:: html

   </div>

Will fetch the fasta sequence.

 

.. rubric:: Adapting extract\_from\_pdb to your own environment

If you have the PDB installed locally, simply set the variable PDB\_DIR
to the absolute location of the directory in which the PDB is installed.
The PDB can either be installed in its divided form or in its full form.

If the file you are looking for is neither in the current directory nor
in the local PDB version, extract\_from\_pdb will try to fetch it from
rcsb. If you do not want this to happen, you should either set the
environment variable NO\_REMOTE\_PDB\_DIR to 1 or use the
-no\_remote\_pdb\_dir flag:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

export NO\_REMOTE\_PDB\_FILE=1

or

t\_coffee -other\_pg extract\_from\_pdb -infile 1PPGE -fasta
-no\_remote\_pdb\_file

.. raw:: html

   </div>

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 3.0pt;
   padding:1.0pt 4.0pt 1.0pt 4.0pt;background:#CCCCCC;margin-left:0cm;margin-right:
   4.25pt">

Manipulating RNA sequences with seq\_reformat

.. raw:: html

   </div>

.. rubric:: Producing a Stockholm output: adding predicted secondary
   structures

.. rubric:: Producing a consensus structure

 

Given an RNA multiple sequence alignment, it is possible to compute the
alifold (Vienna package) consensus secondary structure and output in in
stockholm:

 

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sample\_rnaseq2.aln
-action +aln2alifold -output stockholm\_aln

.. raw:: html

   </div>

 

.. rubric:: Adding a consensus structure to an alignment

 

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sample\_rnaseq2.aln
-action +add\_alifold -output stockholm\_aln

.. raw:: html

   </div>

.. rubric:: Adding a pre-computed consensus structure to an alignment

The file sample\_rnaseq2.aalifold contains the raw output of the alifold
program captured as follows:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

RNAalifold <sample\_rnaseq2.aln  > sample\_rnaseq2.alifold

.. raw:: html

   </div>

It is possible to add this secondary structure to an alignment using:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sample\_rnaseq2.aln -in2
sample\_rnaseq2.alifold -input2 alifold -action +add\_alifold -output
stockholm\_aln

.. raw:: html

   </div>

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.5pt;
   padding:0cm 0cm 0cm 0cm;background:#F7CBB7">

WARNING:

The alifold structure and the alignment MUST be compatible. The function
makes no attempt to thread or align the structure. It merely stack it
below the MSA.

 

.. raw:: html

   </div>

It is also possible to stack Stockholm formatted secondary structures:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: seq\_reformat -in sample\_rnaseq2.aln -in2
sample\_rnaseq2.cons.stk -action +add\_alifold -output stockholm\_aln

.. raw:: html

   </div>

 

.. rubric:: Analyzing an RNAalifold secondary structure prediction

 

the following commands can either be applied on a Stockholm or a
standard MSA. In the second case (standard MSA) the secondary structure
will be automatically re-computed by alifold.

.. rubric:: Analyzing matching columns

 

***+alifold2cov\_stat*** will estimate the number of pairs of columns
that are perfect Watson and Crick, those that are neutral (including a
GU) and those that include correlated mutations. The WCcomp are the
compensated mutations maintaining WC base pairing

 

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sample\_rnaseq2.stk
-action +alifold2analyze stat

.. raw:: html

   </div>

Other arguments can given, to display the list of paired positions and
their status (compensated, Watson, etc)

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sample\_rnaseq2.stk
-action +alifold2analyze list

.. raw:: html

   </div>

 

.. rubric:: Visualizing compensatory mutations

The following command will output a color coded version of your
alignment with matching columns indicated as follows:

         I: Incompatible pair (i.e. at least one pair is not WC)

N: pairs are Gus or WC

W: All pairs are Watson

c : Compensatory mutations

         C: WC compensatory mutations

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sample\_rnaseq2.aln
-action +alifold2analyze aln

.. raw:: html

   </div>

It is possible to turn this output into a colored one using:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sample\_rnaseq2.aln
-action +alifold2analyze color\_htm

.. raw:: html

   </div>

.. rubric:: Handling gapped columns

by default gapped column are ignored but they can be included by adding
the tag usegap

.. rubric:: Comparing alternative folds

The folds associated with alternative alignments can be compared. This
comparison involves counting how many identical pairs of residues are
predicted on each sequence in one fold and in the other. The folds can
either be provided via Stockholm alignments

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

t\_coffee -other\_pg seq\_reformat -in sample\_rnaseq2.cw.stk -in2
sample\_rnaseq2.tcoffee.stk -action +RNAfold\_cmp

.. raw:: html

   </div>

The top of the output (@@lines) summarizes the results that are
displayed on the -in alignment. If the provided alignments do not have a
fold, this fold will be estimated with alifold.

 

 

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 3.0pt;
   padding:1.0pt 4.0pt 1.0pt 4.0pt;background:#CCCCCC;margin-left:0cm;margin-right:
   4.25pt">

Manipulating Phylogenetic Trees with seq\_reformat

.. raw:: html

   </div>

.. rubric:: Producing phylogenetic trees

Seq\_reformat is NOT a phylogeny package, yet over the time it has
accumulated a few functions that make it possible to compute simple
phylogenetic trees, or similar types of clustering:

Given a multiple sequence alignment, it is possible to compute either a
UPGM or an NJ tree:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

seq\_reformat -in <aln> -action +aln2tree -output newick

.. raw:: html

   </div>

Will use an identity matrix to compare your sequences and will output an
unrooted NJ tree in newick format. If you want to produce a rooted UPGMA
tree:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

seq\_reformat -in <aln> -action +aln2tree \_TMODE\_upgma -output newick

.. raw:: html

   </div>

If your data is not data sequence, but a matrix of 1 and Os (i.e. SAR
matrix for instance), you can use a different matrix to compute the
pairwise distances:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

   seq\_reformat -in <aln> -action +aln2tree \_MATRIX\_sarmat -output
newick

.. raw:: html

   </div>

All these parameters can be concatenated:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

   seq\_reformat -in <aln> -action +aln2tree
\_TMODE\_upgma\_MATRIX\_sarmat -output newick

.. raw:: html

   </div>

Bootstrap facilities will also be added at some point … For now we
recommend you use Phylip if you need some serious phylogeny…

 

 

.. rubric:: Comparing two phylogenetic trees

Consider the following file (sample\_tree1.dnd)

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6;margin-left:0cm;margin-right:1.0cm">

 (( A:0.50000, C:0.50000):0.00000,( D:0.00500, E:0.00500):0.99000,
B:0.50000);

.. raw:: html

   </div>

and the file sample\_tree3.dnd.

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6;margin-left:0cm;margin-right:1.0cm">

 (( E:0.50000, C:0.50000):0.00000,( A:0.00500, B:0.00500):0.99000,
D:0.50000);

.. raw:: html

   </div>

You can compare them using:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:6.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

seq\_reformat -in sample\_tree2.dnd -in2 sample\_tree3.dnd -action
+tree\_cmp -output newick

.. raw:: html

   </div>

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6;margin-left:0cm;margin-right:1.0cm">

tree\_cpm\|T: 75 W: 71.43 L: 50.50

tree\_cpm\|8 Nodes in T1 with 5 Sequences

tree\_cmp\|T: ratio of identical nodes

tree\_cmp\|W: ratio of identical nodes weighted with the min Nseq below
node

tree\_cmp\|L: average branch length similarity

(( A:1.00000, C:1.00000):-2.00000,( D:1.00000, E:1.00000):-2.00000,
B:1.00000);

.. raw:: html

   </div>

Please consider the following aspects when exploiting these results:

-The comparison is made on the unrooted trees

T: Fraction of the branches conserved between the two trees. This is
obtained by considering the split induced by each branch and by checking
whether that split is found in both trees

W: Fraction of the branches conserved between the two trees. Each branch
is weighted with MIN the minimum number of leaf on its left or right
(Number leaf left, Number leaf Right)

L: Fraction of branch length difference between the two considered
trees.

 

The last portion of the output contains a tree where distances have been
replaced by the number of leaf under the considered node

Positive values (i.e. 2, 5) indicate a node common to both trees and
correspond to MIN.

Negative values indicate a node found in tree1 but not in tree2

The higher this value, the deeper the node.

You can extract this tree for further usage by typing:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

   cat outfile \| grep -v "tree\_cmp"

.. raw:: html

   </div>

.. rubric:: Scanning Phylogenetic Trees

It is possible to scan an alignment and locally measure the similarity
between an estimated local tree and some reference tree provided from an
external source (or computed on the fly). The following command:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

seq\_reformat -in <aln> -in2 <reftree> -action +tree\_scan
\_MODE\_scan\_\_W\_10\_ > ph\_tree\_scan.txt

.. raw:: html

   </div>

For each position of the alignment, W\*2 blocks of size 2\*1+1 up to
W\*2+1  will be extracted, for each of these block a tree will be
estimated and the similarity of that tree with the reference tree will
be estimated with cmp\_tree. For each position, the tree giving the best
fit will be reported, along with the size of the block leading to that
tree:

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6;margin-left:0cm;margin-right:1.0cm">

 

P: <position>  <block start> <blck\_end> <block score> <block Length>

 

.. raw:: html

   </div>

 

.. rubric:: Pruning Phylogenetic Trees

 

Pruning removes leaves from an existing tree and recomputes distances so
that no information is lost

Consider the file sample\_tree2.dnd:

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6;margin-left:0cm;margin-right:1.0cm">

 (( A:0.50000, C:0.50000):0.00000,( D:0.00500, E:0.00500):0.99000,
B:0.50000);

.. raw:: html

   </div>

And the file sample\_seq8.seq

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6;margin-left:0cm;margin-right:1.0cm">

>A

>B

>C

>D

.. raw:: html

   </div>

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.5pt;
   padding:1.0pt 4.0pt 1.0pt 4.0pt;background:#FFFF89">

Note: Sample\_seq8 is merely a FASTA file where sequences can be
omitted. Sequences can be omitted, but you can also leave them, at your
entire convenience.

.. raw:: html

   </div>

 

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

seq\_reformat -in sample\_tree2.dnd -in2 sample\_seq8.seq -action
+tree\_prune -output newick

.. raw:: html

   </div>

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6;margin-left:0cm;margin-right:1.0cm">

 (( A:0.50000, C:0.50000):0.00000, B:0.50000, D:0.99500);

.. raw:: html

   </div>

 

.. rubric:: 

 

.. raw:: html

   </div>

.. raw:: html

   <div class="WordSection4">

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 3.0pt;
   padding:1.0pt 4.0pt 1.0pt 4.0pt;background:#CCCCCC;margin-left:0cm;margin-right:
   4.25pt">

Building Multiple Sequence Alignments

.. raw:: html

   </div>

.. rubric:: How to generate The Alignment You Need?

.. rubric:: What is a Good Alignment?

This is a trick question. A good alignment is an alignment that makes it
possible to do good biology. If you want to reconstruct a phylogeny, a
good alignment will be an alignment leading to an accurate
reconstruction.

In practice, the alignment community has become used to measuring the
accuracy of alignment methods using structures. Structures are
relatively easy to align correctly, even when the sequences have
diverged quite a lot. The most common usage is therefore to compare
structure based alignments with their sequence based counterpart and to
evaluate the accuracy of the method using these criterions.

Unfortunately it is not easy to establish structure based standards of
truth. Several of these exist and they do not necessarily agree. To
summarize, the situation is as roughly as follows:

         Above 40% identity (within the reference datasets), all the
reference collections agree with one another and all the established
methods give roughly the same results. These alignments can be trusted
blindly.

         Below 40% accuracy within the reference datasets, the reference
collections stop agreeing and the methods do not give consistent
results. In this area of similarity it is not necessarily easy to
determine who is right and who is wrong, although most studies seem to
indicate that consistency based methods (T-Coffee, Mafft-slow and
ProbCons) have an edge over traditional methods.

When dealing with distantly related sequences, the only way to produce
reliable alignments is to us structural information. T-Coffee provides
many facilities to do so in a seamless fashion. Several important
factors need to be taken into account when selecting an alignment
method:

-The best methods are not ***always*** doing best. Given a difficult
dataset, the best method is only more likely to deliver the best
alignment, but there is no guaranty it will do so. It is very much like
betting on the horse with the best odds.

-Secondly, the difference in accuracy (as measured on reference
datasets)  between all the available methods is not incredibly high. It
is unclear whether this is an artifact caused by the use of "easy"
reference alignments, or whether this is a reality. The only thing that
can change dramatically the accuracy of the alignment is the use of
structural information.

Last, but not least, bear in mind that these methods have only been
evaluated by comparison with reference structure based sequence
alignments. This is merely one criterion among many. In theory, these
methods should be evaluated for their ability to produce alignments that
lead to accurate trees, good profiles or good models. Unfortunately,
these evaluation procedures do not yet exist.

 

.. rubric:: The Main Methods and their Scope

There are many MSA packages around. The main ones are ClustalW, Muscle,
Mafft, T-Coffee and ProbCons. You can almost forget about the other
packages, as there is virtually nothing you could do with them that you
will not be able to do with these packages.

These packages offer a complex trade-off between speed, accuracy and
versatility.

.. rubric:: ClustalW: everywhere you look

ClustalW is still the most widely used multiple sequence alignment
package. Yet things are gradually changing as recent tests have
consistently shown that ClustalW is neither the most accurate nor the
fastest package around. This being said, ClustalW is everywhere and if
your sequences are similar enough, it should deliver a fairly reasonable
alignment.

.. rubric:: Mafft and Muscle: Aligning Many Sequences

If you have many sequences to align Muscle or Mafft are the obvious
choice. Mafft is often described as the fastest and the most efficient.
This is not entirely true. In its fast mode (FFT-NS-1), Mafft is similar
to Muscle and although it is fairly accurate it is about 5 points less
accurate than the consistency based packages (ProbCons and T-Coffee). In
its most accurate mode (L-INS-i) Mafft uses local alignments and
consistency. It becomes much more accurate but also slower, and more
sensitive to the number of sequences.

The alignments generated using the fast modes of these programs will be
very suitable for several important applications such as:

         -Distance based phylogenetic reconstruction (NJ trees)

         -Secondary structure predictions

However they may not be suitable for more refined application such as

         -Profile construction

         -Structure Modeling

         -3D structure prediction

         -Function analysis

In that case you may need to use more accurate methods

.. rubric:: T-Coffee and ProbCons: Slow and Accurate

T-Coffee works by first assembling a library and then by turning this
library into an alignment. The library is a list of potential pairs of
residues. All of them are not compatible and the job of the algorithm is
to make sure that as many possible constraints as possible find their
way into the final alignment. Each library line is a constraint and the
purpose is to assemble the alignment that accommodates the more all the
constraints.

It is very much like building a high school schedule, where each
teachers says something "I need my Monday morning", "I can't come on
Thursday afternoon", and so on. In the end you want a schedule that
makes everybody happy, if possible.The nice thing about the library is
that it can be used as a media to combine as many methods as one wishes.
It is just a matter of generating the right constraints with the right
method and compile them into the library.

ProbCons and Mafft (L-INS-i) uses a similar algorithm, but with a
Bayesian twist in the case of Probcons. In practice, however, probcons
and T-Coffee give very similar results and have similar running time.
Mafft is significantly faster.

All these packages are ideal for the following applications:

         -Profile reconstruction

         -Function analysis

         -3D Prediction

.. rubric:: Choosing The Right Package

Each available package has something to go for it. It is just a matter
of knowing what you want to do. T-Coffee is probably the most versatile,
but it comes at a price and it is currently slower than many alternative
packages.

In the rest of this tutorial we give some hints on how to carry out each
of these applications with T-Coffee.

+--------------+--------------+--------------+--------------+--------------+--------------+
|              | Muscle       | Mafft        | ProbCons     | T-Coffee     | ClustalW     |
|              |              |              |              |              |              |
                                                                                         
+--------------+--------------+--------------+--------------+--------------+--------------+
| Accuracy     | ++           | +++          | +++          | +++          | +            |
+--------------+--------------+--------------+--------------+--------------+--------------+
| <100 Seq.    | ++           | ++           | +++          | +++          | +            |
+--------------+--------------+--------------+--------------+--------------+--------------+
| >100 Seq.    | +++          | +++          | -            | +            | +            |
+--------------+--------------+--------------+--------------+--------------+--------------+
| Remote       | ++           | +++          | +++          | +++          | +            |
| Homologues   |              |              |              |              |              |
+--------------+--------------+--------------+--------------+--------------+--------------+
| MSA vs Seq.  | -            | -            |              | +++          | +++          |
+--------------+--------------+--------------+--------------+--------------+--------------+
| MSA vs MSA   | -            | -            | -            | +++          | +++          |
+--------------+--------------+--------------+--------------+--------------+--------------+
| >2 MSAs      | -            | -            | -            | +++          | -            |
+--------------+--------------+--------------+--------------+--------------+--------------+
| Seq. vs      | -            | -            | -            | +++          | +            |
| Struc.       |              |              |              |              |              |
+--------------+--------------+--------------+--------------+--------------+--------------+
| Splicing     | -            | +++          | -            | +++          | -            |
| Var.         |              |              |              |              |              |
+--------------+--------------+--------------+--------------+--------------+--------------+
| Reformat     | -            | -            | -            | +++          | ++           |
+--------------+--------------+--------------+--------------+--------------+--------------+
| Phylogeny    | -            | -            | -            | +            | ++           |
+--------------+--------------+--------------+--------------+--------------+--------------+
| Evaluation   | -            | -            | +            | +++          | -            |
+--------------+--------------+--------------+--------------+--------------+--------------+
| Speed        | +++          | +++          | +            | +            | ++           |
+--------------+--------------+--------------+--------------+--------------+--------------+

**Table 1.** Relative possibilities associated with the main packages
(T-Coffee Tutorial, C. Notredame, www.tcoffee.org). In any of the
situations corresponding to each table line, (+++) indicates that the
method is the best suited, (++) indicates that the method is not optimal
but behaves reasonably well, (+) indicates that it is possible but not
recommended (-) indicates that the option is not available.

 

+--------------+--------------+--------------+--------------+--------------+--------------+
|              | Muscle       | Mafft        | ProbCons     | T-Coffee     | ClustalW     |
|              |              |              |              |              |              |
                                                                                         
+--------------+--------------+--------------+--------------+--------------+--------------+
| Dist Based   | +++          | +++          | ++           | ++           | ++           |
| Phylogeny    |              |              |              |              |              |
+--------------+--------------+--------------+--------------+--------------+--------------+
| ML or MP     | ++           | +++          | +++          | +++          | ++           |
| Phylogeny    |              |              |              |              |              |
+--------------+--------------+--------------+--------------+--------------+--------------+
| Profile      | ++           | +++          | +++          | +++          | ++           |
| Construction |              |              |              |              |              |
+--------------+--------------+--------------+--------------+--------------+--------------+
| 3D Modeling  | ++           | ++           | ++           | +++          | +            |
+--------------+--------------+--------------+--------------+--------------+--------------+
| Secondary    | +++          | +++          | ++           | ++           | ++           |
| Structure P  |              |              |              |              |              |
+--------------+--------------+--------------+--------------+--------------+--------------+

**Table 2.** Most Suitable Appplications of each package (T-Coffee
Tutorial, C. Notredame, www.tcoffee.org). In any of the situations
corresponding to each table line, (+++) indicates that the method is the
best suited, (++) indicates that the method is not optimal but behaves
reasonably well, (+) indicates that it is possible but not recommended
(-) indicates that the option is not available.

 

.. rubric:: Computing Multiple Sequence Alignments With T-Coffee

.. rubric:: Computing Very accurate (but slow) alignments with
   PSI-Coffee

PSI-Coffee builds a profile associated with each of your input sequence
and then makes a multiple profile alignment. If you do not have any
structure, it is the most accurate mode of T-Coffee.

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sproteases\_small.fasta -mode psicoffee

.. raw:: html

   </div>

If you want to go further, and be even slower, you can use the accurate
mode that will combine profile and structural information

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sproteases\_small.fasta -mode accurate

.. raw:: html

   </div>

It is probably one of the most accurate way of aligning sequences
currently available.

.. rubric:: A Simple Multiple Sequence Alignment

T-Coffee is meant to be run like ClustalW. This means you can use it
like ClustalW for most simple applications. For instance, the following
instruction

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sproteases\_small.fasta

.. raw:: html

   </div>

This instruction will compute a multiple sequence alignment of your
sequences, using the default mode of T-Coffee. It will output the
alignment on the screen and in a file named sproteases\_small.aln. This
file contains your alignment in ClustalW format.

The program will also output a file named sproteases\_small.dnd that
contains the guide tree used to assemble the progressive alignment.

.. rubric:: Controlling the Output Format

If you need to, you can also trigger different ouput formats using the
-output flag:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sproteases\_small.fasta
-output=clustalw,fasta\_aln,msf

.. raw:: html

   </div>

You can specify as many formats as you want.

.. rubric:: Computing a Phylogenetic tree

T-Coffee is not a phylogeny package. Yet, it has some limited abilities
to turn your MSA into a phylogenetic tree. This tree is a Neighbor
Joining Phylogenetic tree, very similar to the one you could compute
using ClustalW.

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sproteases\_small.fasta
-output=clustalw,fasta\_aln,msf

.. raw:: html

   </div>

The phylogenetic tree is the file with the ph extension. Never use the
.dnd tree in place of a genuine phylogenetic tree. The phylogenetic tree
output by T-Coffee is only an indication. You should produce a
bootstrapped phylogenetic tree using packages like Phylip
(bioweb.pasteur.fr/seqanal/phylogeny/phylip-uk.html). You can visualize
your tree using online tree drawing programs like phylodendron
(iubio.bio.indiana.edu/treeapp/treeprint-form.html).

 

.. rubric:: Using Several Datasets

If your sequences are spread across several datasets, you can give all
the files via the -seq flag:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -seq=sprotease1\_small.fasta,sprotease2\_small.aln
-output=clustalw,fasta\_aln,msf

.. raw:: html

   </div>

Note that you can give as many file as you want (the limit is 200) and
that the files can be in any format. If you give an alignment, the gaps
will be reset and your alignment will only provide sequences.

Sequences with the same name between two files are assumed to be the
same sequence. If their sequences differ, they will be aligned and
replaced by the consensus of that alignment. This process is known as
sequence reconciliation.

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.5pt;
   padding:0cm 0cm 0cm 0cm;background:#F7CBB7">

You should make sure that there are no duplicates in your alignment,
especially when providing multiple datasets.

.. raw:: html

   </div>

 

.. rubric:: How Good is Your Alignment

Later in this tutorial we show you how to estimate the accuracy of your
alignment. Before we go into details, you should know that the number
that comes on the first line of the header (in ClustalW format) is the
score of your alignment.

CLUSTAL FORMAT for T-COFFEE Version\_4.32 [http://www.tcoffee.org],
CPU=19.06 sec, **SCORE=37**, Nseq=19, Len=341

You can use this value to compare alternative alignments of the same
sequences. Alignments with a score higher than 40 are usually pretty
good.

.. rubric:: Doing it over the WWW

You can run T-Coffee online at www.tcoffee.org. Use the regular or the
advanced form of the T-Coffee server.

.. rubric:: Aligning Many Sequences

.. rubric:: Aligning Very Large Datasets with Muscle

T-Coffee is not a good choice if you are dealing with very large
datasets, use Mafft or Muscle. To align a large dataset with Muscle,
try:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

muscle -infile sproteases\_large.fasta > sproteases\_large.muscle

.. raw:: html

   </div>

To use the fastest possible mode (less accurate) run:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

muscle -in sproteases\_large.fasta -maxiters 1 -diags -sv -distance1
kbit20\_3 > sproteases\_large.muscle

.. raw:: html

   </div>

.. rubric:: Aligning Very Large Alignments with Mafft

The fastest mode with Mafft can be achieved using:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

mafft --retree 2 input > output

.. raw:: html

   </div>

.. rubric:: Aligning Very Large Alignments with T-Coffee

 

T-Coffee is not very well gifted for aligning large datasets, but you
can give it a try using a special option that generates approximate
alignments. These alignments should roughly have the same accuracy as
ClustalW. They are acceptable for sequences more than 40% identical.

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sproteases\_large.fasta -mode quickaln

.. raw:: html

   </div>

.. rubric:: Shrinking Large Alignments With T-Coffee

Once you have generated your large alignment, you may nedd/want to
shrink it to a smaller one, that will be (hopefuly) as informative and
easier to manipulate. For that purpose, use the trim option (described
in detail in the first section of this document).

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in sproteases\_large.muscle
-action +trim \_n20 -output > sproteases\_large\_muscle\_trim.aln

.. raw:: html

   </div>

.. rubric:: Modifying the default parameters of T-Coffee

The main parameters of T-Coffee are similar to those of ClustalW. They
include the substitution matrix and the gap penalties. In general,
T-Coffee's default is adequate. If, however, you are not satisfied with
the default parameters, we encourage you to change the following
parameters. Interestingly, most of what we say here holds reasonably
well for ClustalW.

.. rubric:: Changing the Substitution Matrix

T-Coffee only uses the substitution matrix to make the pairwise
alignments that go into the library. These are all the global alignments
of every possible pair of sequences, and the ten best local alignments
associated with every pair of sequences.

By default, these alignments are computed using a Blosum62 matrix, but
you can use any matrix you fancy instead, including: pam120mt, pam160mt,
pam250mt, pam350mt, blosum30mt, blosum40mt, blosum45mt, blosum50mt,
blosum55mt, blosum62mt, blosum80mt, or even user-provided matrices in
the BLAST format, as described in the technical manual.

***Pam matrices***: These matrices are allegedly less accurate than the
blosum. The index is correlated to the evolutionary distances. You
should therefore use the pam350mt to align very distantly related
sequences.

***Blosum matrices***: These matrices are allegedly the most accurate.
The index is correlated to the maximum percent identity within the
sequences used to estimate the matrix. You should therefore use the
Blosum30mt to align very distantly related sequences. Blosum matrices
are biased toward protein core regions. This may explain why theses
matrices tend to give better alignments, since by design they can
capture the most evolutionary resilient signal contained in proteins.

Unless you have some structural information available, the only way to
tell whether your alignment has improved or not is to look at the score.
For instance, if you compute the two following alignments:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sproteases\_small.fasta -matrix=blosum30mt
-outfile=b30.aln

PROMPT: t\_coffee sproteases\_small.fasta -matrix=blosum80mt
-outfile=b80.aln

PROMPT: t\_coffee sproteases\_small.fasta -matrix=pam350mt -outfile
p350.aln

.. raw:: html

   </div>

You will get two alignments that have roughly the same score but are
different. You can still use these two alternative alignments by
comparing them to identify regions that have been aligned identically by
the two matrices. These regions are usually more trustworthy.

.. rubric:: Comparing Two Alternative Alignments

If you change the parameters, you will end up with alternative
alignemnts. It can be interesting to compare them quantitatively.
T-Coffee comes along with an alignment comparison module named
aln\_compare. You can use it to estimate the amount of difference
between your two alignments:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg aln\_compare -al1 b30.aln -al2 p350.aln

.. raw:: html

   </div>

This comparison will return the following result:

 

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6">

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

seq1       seq2          Sim   [ALL]           Tot

b30           19         32.6    93.7 [100.0]   [40444]

.. raw:: html

   </div>

Where 93.7 is the percentage of similarity (sums of pairs) between the
two alignments. It means that when considering every pair of aligned
residues in b30 (40444), the program found that 93.7% of these pairs
could be found in the alignment p350.aln.

Of course, this does not tell you where are the good bits, but you can
get this information with the same program:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

t\_coffee -other\_pg aln\_compare -al1 b30.aln -al2 p350.aln
-output\_aln -output\_aln\_threshold 50

.. raw:: html

   </div>

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6">

sp\|O35205\|GRAK\_MOUSE  
M---r----fssw-------ALvslvagvym----------------SSECFHTEIIGGR

sp\|Q7YRZ7\|GRAA\_BOVIN  
M--ni----pfpf--sfppaIClllipgvfp----------------vs---cEGIIGGN

sp\|P08884\|GRAE\_MOUSE  
M--------ppv----------lilltlllp----------------l-GAGAEEIIGGH

sp\|Q06606\|GRZ2\_RAT    
M--------flf----------lfflvailp----------------v-NTEGGEIIWGT

sp\|P21844\|MCPT5\_MOUSE 
M---h----llt----------lhllllllg----------------s-STKAGEIIGGT

sp\|P03953\|CFAD\_MOUSE  
M---h----ssvy-------fvalvilgaav----------------CAAQPRGRILGGQ

sp\|P00773\|ELA1\_RAT    
M---l----rflv--F----ASlvlyghstq----------------DFPETNARVVGGA

sp\|Q00871\|CTRB1\_PENVA 
MIgkl----slll--V----CVavasgnpaagkpwhwKSPKPLVDPRIHVNATPRIVGGV

sp\|P08246\|ELNE\_HUMAN  
M--tlGR--rlac--L----FLacvlpalll----------------GGTALASEIVGGR

sp\|P20160\|CAP7\_HUMAN  
M--t-----rltv--L----ALlagllassr----------------AGSSPLLDIVGGR

sp\|P80015\|CAP7\_PIG    
-------------------------------------------------------IVGGR

sp\|Q03238\|GRAM\_RAT    
l-------------------LLllalktlwa----------------VGNRFEAQIIGGR

sp\|P00757\|KLKB4\_MOUSE 
M-----------w-------flilflalslggid-------------AAPP-----vqsq

sp\|Q6H321\|KLK2\_HORSE  
M-----------w-------flvlcldlslgetg-------------ALPPIQSRIIGGW

sp\|Q91VE3\|KLK7\_MOUSE  
M---------gvw-------llslitvllslale-------------tag-QGERIIDGY

sp\|Q9Y5K2\|KLK4\_HUMAN  
M-ataGN--pwgw-------flgylilgvag-sl-------------vsg-SCSQIINGE

sp\|P29786\|TRY3\_AEDAE  
M-------nqflfVSF---------calldsakvsaa------------tLSSGRIVGGF

sp\|P35037\|TRY3\_ANOGA  
M---iSNKiaillAVLvvav----acaqarvaqqhrsVQALPRFLPRPKYDVGHRIVGGF

sp\|P07338\|CTRB1\_RAT   
M--a------flwlvs---------cfalvgatfgcg---vptiqpv--LTGLSRIVNGE

                 
                                                             : .

 

sp\|O35205\|GRAK\_MOUSE  
EVQPHSRPFMASIQYR----SKHICGGVLIHPQWVLTAAHCYSWFprGHSPTVVLGAHSL

sp\|Q7YRZ7\|GRAA\_BOVIN  
EVAPHTRRYMALIK------GLKLCAGALIKENWVLTAAHCDlk----GNPQVILGAHST

sp\|P08884\|GRAE\_MOUSE 
 VVKPHSRPYMAFVKSVDIEGNRRYCGGFLVQDDFVLTAAHCRN-----RTMTVTLGAHNI

sp\|Q06606\|GRZ2\_RAT    
ESKPHSRPYMAFIKFYDSNSEPHHCGGFLVAKDIVMTAAHCNG-----RNIKVTLGAHNI

sp\|P21844\|MCPT5\_MOUSE 
ECIPHSRPYMAYLEIVTSENYLSACSGFLIRRNFVLTAAHCAG-----RSITVLLGAHNK

sp\|P03953\|CFAD\_MOUSE  
EAAAHARPYMASVQVN----GTHVCGGTLLDEQWVLSAAHCMDGVtdDDSVQVLLGAHSL

sp\|P00773\|ELA1\_RAT    
EARRNSWPSQISLQYLSggswyHTCGGTLIRRNWVMTAAHCVSSQm---TFRVVVGDHNL

sp\|Q00871\|CTRB1\_PENVA 
EATPHSWPHQAALFId----DMYFCGGSLISSEWVLTAAHCMDGAg---FVEVVLGAHNI

sp\|P08246\|ELNE\_HUMAN  
RARPHAWPFMVSLQLr----GGHFCGATLIAPNFVMSAAHCVANVNV-RAVRVVLGAHNL

sp\|P20160\|CAP7\_HUMAN  
KARPRQFPFLASIQNq----GRHFCGGALIHARFVMTAASCFQSQNP-GVSTVVLGAYDL

sp\|P80015\|CAP7\_PIG    
RAQPQEFPFLASIQKq----GRPFCAGALVHPRFVLTAASCFRGKNS-GSASVVLGAYDL

sp\|Q03238\|GRAM\_RAT    
EAVPHSRPYMVSLQNT----KSHMCGGVLVHQKWVLTAAHCLSEP--LQQLKLVFGLHSL

sp\|P00757\|KLKB4\_MOUSE 
vdcENSQPWHVAVYRF----NKYQCGGVLLDRNWVLTAAHCYN-----DKYQVWLGKNNF

sp\|Q6H321\|KLK2\_HORSE  
ECEKHSKPWQVAVYHQ----GHFQCGGVLVHPQWVLTAAHCMS-----DDYQIWLGRHNL

sp\|Q91VE3\|KLK7\_MOUSE  
KCKEGSHPWQVALLKG----NQLHCGGVLVDKYWVLTAAHCKM-----GQYQVQLGSDKI

sp\|Q9Y5K2\|KLK4\_HUMAN  
DCSPHSQPWQAALVME----NELFCSGVLVHPQWVLSAAHCFQ-----NSYTIGLGLHSL

sp\|P29786\|TRY3\_AEDAE  
QIDIAEVPHQVSLQRS----GRHFCGGSIISPRWVLTRAHCTTNTDP-AAYTIRAGStd-

sp\|P35037\|TRY3\_ANOGA  
EIDVSETPYQVSLQYF----NSHRCGGSVLNSKWILTAAHCTVNLQP-SSLAVRLGSsr-

sp\|P07338\|CTRB1\_RAT   
DAIPGSWPWQVSLQDKt---gfHFCGGSLISEDWVVTAAHCGVKT----SDVVVAGEFDQ

                                   :           \*.. ::    ::: \*
\*           :  \*

.. raw:: html

   </div>

This is the alignment al1, but residues that have lost more than 50% of
their pairing partner between the two alignments are now in lower case.
In the section of this tutorial entitled comparing alignments, we show
you more sophisticated ways to do this comparison.

For an even more drastic display, try:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

t\_coffee -other\_pg aln\_compare -al1 b30.aln -al2 p350.aln
-output\_aln -output\_aln\_threshold 50 -output\_aln\_modif x

.. raw:: html

   </div>

 

.. rubric:: Changing Gap Penalties

Gap penalties are the core of the matter when it comes to multiple
sequence alignments. An interesting feature of T-Coffee is that it does
not really need such penalties when assembling the MSA, because in
theory the penalties have already been applied when computing the
library. This is the theory, as in practice penalties can help improve
the quality of the alignment.

The penalties can be changed via the flags -gapopen for the gap opening
penalty and via -gapext for the gap extension penalty. The range for
gapopen are [-500,--5000], the range for the extension should rather be
[-1, -10]. These values do not refer to a substitution matrix, but
rather to the values range of the concistensy estimation (i.e. a ratio)
normalized to 10000 for a maximum consistency.

The default values are -gapopen=-50, -gapext=0. The reasons for these
very low values are that they are meant to be cosmetic only, since a
trademark of T-Coffee (inherited from Dialign) is not to need explicit
penalties. Yet, we know for a fact that alignments with higher gap
penalties often look nicer (for publications) and are sometimes more
accurate. For instance, you can try:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sproteases\_small.fasta -gapopen -100 -gapext -5

.. raw:: html

   </div>

This gap penalty is only applied at the alignment level (i.e. after the
library was computed). If you want to change the gap penalties of the
methods used to build the library, you will need to go deeper into the
core of the matter...

Two methods are used by default to build the library. One does global
pairwise alignments and is named slow\_pair, the other is named
lalign\_id\_pair and and produces local alignments. These methods are
specified via the -method flag. The default of this flag is:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sproteases\_small.fasta
-method=lalign\_id\_pair,slow\_pair

.. raw:: html

   </div>

Usually you do not need to write it because it is the default, but if
you want to change the default parameters of the constituting methods,
you will need to do so explicitely. The default for lalign\_id\_pair is
gop=-10, GEP=-4, MATRIX=blosum50mt. The default for slow\_pair is:
GOP=-10, GEP=-1 and MATRIX=blosum62mt. If you want to change this, try:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sproteases\_small.fasta -method
lalign\_id\_pair@EP@MATRIX@blosum62mt,slow\_pair -outfile
sproteases\_small.b62\_aln

.. raw:: html

   </div>

This means the library is now computed using the Blosum62mt with lalign,
rather than the Blosum50mt. The good news is that when using this
matrix, the score of our alignment increases from 48(default) to 50. We
may assume this new alignment is more accurate than the previous one.

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.5pt;
   padding:0cm 0cm 0cm 0cm;background:#F7CBB7">

WARNING: It only makes sense to compare the consistency score of
alternative alignments when these alignments have been computed using
the same methods (lalign\_id\_pair and slow\_pair for instance).

.. raw:: html

   </div>

.. rubric:: Can You Guess The Optimal Parameters?

It is a trick question, but the general answer is NO. The matrix and the
gap penalties are simplistic attempts at modeling evolution. While the
matrices do a reasonable job, the penalties are simply inappropriate:
they should have a value that depends on the structure of the protein
and a uniform value cannot be good enough. Yet, since we do not have
better we must use them…

In practice, this means that parameter optimality is a very add-hoc
business. It will change from one dataset to the next and there is no
simple way to predict which matrix and which penalty will do better. The
problem is also that even after your alignment has been computed, it is
not always easy to tell whether your new parameters have improved or
degraded your MSA. There is no systematic way to evaluate an MSA.

In general, people visually evaluate the alignment, count the number of
identical columns and consider that one more conserved column is good
news. If you are lucky you may know a few functional features that you
expect to see aligned. If you are very lucky, you will have one
structure and you can check the gaps fall in the loops. If you are
***extremely*** lucky, you will have two structures and you can assess
the quality of your MSA.

An advantage of T-Coffee is the fact that the overall score of the
alignment (i.e. the consistency with the library) is correlated with the
overall accuracy. In other words, if you alignment score increases, its
accuracy probably increases also. All this being said, consistency is
merely an empirical way of estimating the change of parameters and it
does not have the predictive power of a BLAST E-Value.

.. rubric:: Using Many Methods at once

One of the most common situation when building multiple sequence
alignments is to have several alignments produced by several alternative
methods, and not knowing which one to choose. In this section, we show
you that you can use M-Coffee to combine your many alignments into one
single alignment. We show you here that you can either let T-Coffee
compute all the multiple sequence alignments and combine them into one,
or you can specify the methods you want to combine. M-Coffee is not
always the best methods, but extensive benchmarks on BaliBase, Prefab 
and Homstrad have shown that it delivers the best alignment 2 times out
of 3. If you do not want to use the methods provided by M-Coffee, you
can also combine pre-computed alignments.

.. rubric:: Using All the Methods at the Same Time: M-Coffee

In M-Coffee, M stands for Meta. To use M-Coffee, you will need several
packages to be installed (see documentation). The following command:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sproteases\_small.fasta -mode mcoffee -output
clustalw, html

.. raw:: html

   </div>

Will compute a Multiple Sequence Alignment with the following MSA
packages:

clustalw, poa, muscle, probcons, mafft, dialing-T, pcma and T-Coffee.

For those using debian, another mode is available

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sproteases\_small.fasta -mode dmcoffee -output
clustalw, html

.. raw:: html

   </div>

Will compute a Multiple Sequence Alignment with the following MSA
packages:

kalign, poa, muscle, probcons, mafft, dialing-T, and T-Coffee.

 

 

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6">

Package              Where From

==========================================================

ClustalW             can interact with t\_coffee

----------------------------------------------------------

Poa                  http://www.bioinformatics.ucla.edu/poa/

----------------------------------------------------------

Muscle        http://www.bioinformatics.ucla.edu/poa/

----------------------------------------------------------

ProbCons             http://probcons.stanford.edu/

----------------------------------------------------------

MAFFT 
\ `http://www.biophys.kyoto- <http://www.biophys.kyoto-/>`__\ u.ac.jp/~katoh/programs/align/mafft/

----------------------------------------------------------

Dialign-T            \ http://dialign-t.gobics.de/

----------------------------------------------------------

PCMA                
\ `ftp://iole.swmed.edu/pub/PCMA/ <file://localhost/pub/PCMA>`__

----------------------------------------------------------

T-Coffee             www.tcoffee.org

----------------------------------------------------------

Kalign               msa.cgb.ki.se/cgi-bin/msa.cgi

----------------------------------------------------------

amap                 bio.math.berkeley.edu/amap/

----------------------------------------------------------

.. raw:: html

   </div>

When this is done, all the alignments will be combined into one. If you
open the file sproteases\_small.html with your favorite web browser, you
will see a colored version of your alignment.

The alignment is colored according to its consistency with all the MSA
used to compute it. Regions in red have a high consistency and you can
expect them to be fairly accurate. Regions in green/blue have the lowest
consistency and you should not trust them.

Overall this alignment has a score of 80, which means that it is 80%
consistent with the entire collection. This is a fairly high index,
which means you can probably trust your alignment (at least where it is
red).

.. rubric::  Using Selected Methods to Compute your MSA

Using the 8 Methods of M-Coffee8 can sometimes be a bit heavy. If you
only want to use a subset of your favorite methods, you should know that
each of these methods is available via the -method flag. For instance,
to combine MAFFT, Muscle, t\_coffee and ProbCons, you can use:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sproteases\_small.fasta
-method=t\_coffee\_msa,mafft\_msa,probcons\_msa,muscle\_msa -output=html

.. raw:: html

   </div>

This will result in a computation where all the specified methods are
mixed together

.. rubric:: Combining pre-Computed Alignments

You may have a bunch of alignments that you have either pre-computed, or
assembled manually or received from a colleague. You can also combine
these alignments. For instance, let us imagine we generated 4 alignments
with ClustalW using different gap penalties:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

clustalw -infile=sproteases\_small.fasta -gapopen=0 -outfile=g0.aln

clustalw -infile=sproteases\_small.fasta -gapopen=-5 -outfile=g5.aln

clustalw -infile=sproteases\_small.fasta -gapopen=-10 -outfile=g10.aln

clustalw -infile=sproteases\_small.fasta -gapopen=-15 -outfile=g15.aln

 

.. raw:: html

   </div>

To combine them into ONE single alignment, use the -aln flag:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sproteases\_small.fasta -aln g0.aln g5.aln g10.aln
g15.aln -output clustalw html

.. raw:: html

   </div>

As before, the score indicates a high level of consistency (91%) between
all these alignments. This is an indication that the final alignment is
probably correct.

.. rubric:: Aligning Profiles

Sometimes, it is better to pre-align a subset of your sequences, and
then to use this small alignment as a master for adding sequences
(sequence to profile alignment) or even to align several profiles
together if your protein family contains distantly related groups.
T-Coffee contains most of the facilities available in ClustalW to deal
with profiles, and the strategy we outline here can be used to deal with
large datasets

.. rubric:: Using Profiles as templates

 

.. rubric:: Aligning One sequence to a Profile

Assuming you have a multiple alignment (sproteases\_small.aln) here is a
simple strategy to align one sequence to your profile:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sproteases\_oneseq.fasta  -profile
sproteases\_small.aln

.. raw:: html

   </div>

.. rubric:: Aligning Many Sequences to a Profile

You can align as many sequences as you wish to your profile. Likewise,
you can have as many profiles as you want. For instance, the following:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sequences.fasta  -profile=prf1.aln,prf2.aln,prf3.aln
-outfile=combined\_profiles.aln

.. raw:: html

   </div>

Will make a multiple alignment of 3 profiles and 5 sequences. You can
mix sequences and profiles in any proportion you like. You can also use
all the methods you want although you should be aware that when using
external methods (see the external method section in this tutorial), the
profile is replaced with its consensus sequence, which will not be quite
as accurate.

Methods supporting full profile information are: lalign\_id\_pair,
slow\_pair and proba\_pair, clustalw\_pair and clustalw\_msa. All the
other methods (internal or external) treat the profile as a consensus
(less accurate).

 

.. rubric:: Aligning Other Types of Sequences

.. rubric:: Splicing variants

Splicing variants are especially challenging for most MSA programs. This
is because the splicing variants need very long gaps to be inserted,
while most programs attempt to match as many symbols as possible.

Standard programs like ClustalW or Muscle are not good at dealing with
this situation and in our experience, the only programs that can do
something with splice variants are those using local information like
some flavors of Mafft and T-Coffee .

For instance, if you try muscle on the following dataset:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

muscle -in sv.fasta -clw

.. raw:: html

   </div>

You will quickly realise that your alignment is not very good and does
not show where the alternative splicing coocurs. On the other hand, if
you use T-Coffee, things become much clearer

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee  sv.fasta

.. raw:: html

   </div>

The reason why T-Coffee does better than other packages is mostly
because it uses local information (lalign\_id\_pair) and is therefore
less sensitive to long gaps. If the default mode does not work for your
dataset, you can try to be a bit more aggressive and only use local
information to compute your library:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sv.fasta -method lalign\_id\_pair

.. raw:: html

   </div>

Of course, the most distantly related your sequences, the harder the
alignment of splicing variants

.. rubric:: Aligning DNA sequences

Multiple Sequence Alignment methods are not at their best when aligning
DNA sequences. Whenever you can, try using a local multiple sequence
alignment package like the Gibbs sampler. Yet if you believe your DNA
sequence are homologous over their entire length, you can use T-Coffee.

In theory, the program automatically recognizes DNA sequences and uses
appropriate methods, yet adding the -type=dna flag cannot do any harm...

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sample\_dnaseq1.fasta –type=dna

.. raw:: html

   </div>

The type declaration (or its automatic detection) triggers the use of
the appropriate substitution matrix in most of the methods. In practice,
any time it encounters dna, the program will try to use “4dna” version
of the requested methods. These methods have lower penalties and are
better suited for dealing with nucleic acid sequences.

However, if you would rather use your own matrix, use:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sample\_dnaseq1.fasta –in
Mlalign\_id\_pair4dna@EP@MATRIX@idmat

.. raw:: html

   </div>

Where you should replace idmat with your own matrix, in BLAST format
(see the format section of the Reference Manual).

.. rubric:: Aligning RNA sequences

RNA sequences are very important and almost every-where these days. The
main property of RNA sequences is to have a secondary structure that can
be used to guide the alignment. While the default T-Coffee has no
special RNA alignment method incorporated in, smart people have thought
about this. If you are interested in RNA, check:
http://www.bio.inf.uni-jena.de/Software/MARNA/.

.. rubric:: Noisy Coding DNA Sequences…

When dealing with coding DNA, the right thing to do is to translate your
DNA sequence and thread the DNA onto the protein alignment if you really
need some DNA. However, sometimes, your cDNA may not be so clean that
you can easily translate it (frameshifts and so on). Whenever this
happens, try (no warranty) the following special method.

The test case in three\_dna\_seq.fasta contains the DNA sequences of
three proteases with a couple of frameshifts here and there. If you make
a regular alignment of these sequences

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee three\_cdna.fasta

.. raw:: html

   </div>

You can immediately see that many gaps have sizes that are not multiple
of 3 (codon size). Most of the information is lost. On the other hand,
when using an appropriate alignment method that takes into account all
the frames at the same time, we get something much more meaningful:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee three\_cdna.fasta -method cdna\_fast\_pair

.. raw:: html

   </div>

And most importantly, the frameshifts end up at the right place. You can
even recover the corrected protein sequence using a special mode of
seq\_reformat:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in three\_cdna.aln -action
+clean\_cdna +translate

.. raw:: html

   </div>

+clean cdna is a small HMM that loops through each sequence and select
the frame in order to maximize the similarity within the alignment.

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 3.0pt;
   padding:1.0pt 4.0pt 1.0pt 4.0pt;background:#CCCCCC;margin-left:0cm;margin-right:
   4.25pt">

Using Secondary Structure Predictions:

.. raw:: html

   </div>

 

T-Coffee can be used to predict secondary structures and transmembrane
domains. For secondary structure predictions, the current implementation
is only able to run GOR on either single sequences or on a bunch of
homologues found by BLAST.

.. rubric:: Single Sequence prediction

To make a secondary structure prediction with GOR, run the following. In
this command line SSP is a hard coded mode. It prompts the computation
of predicted secondary structures.

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

t\_coffee sample\_aln.fasta -template\_file SSP

.. raw:: html

   </div>

The predictions are then displayed in the files:

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 31.0pt;
   background:#E6E6E6;margin-left:0cm;margin-right:1.0cm">

 

#### File Type= Template Protein Secondary Structure Format=  fasta\_seq
Name= hmgb\_chite.ssp

#### File Type= Template Protein Secondary Structure Format=  fasta\_seq
Name= hmgl\_trybr.ssp

#### File Type= Template Protein Secondary Structure Format=  fasta\_seq
Name= hmgl\_trybr3.ssp

#### File Type= Template Protein Secondary Structure Format=  fasta\_seq
Name= hmgl\_wheat.ssp

#### File Type= Template Protein Secondary Structure Format=  fasta\_seq
Name= hmgl\_wheat2.ssp

#### File Type= Template Protein Secondary Structure Format=  fasta\_seq
Name= hmgt\_mouse.ssp

 

.. raw:: html

   </div>

Transmembrane structures can be carried out with:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

t\_coffee sample\_aln.fasta -template\_file TM

.. raw:: html

   </div>

.. rubric:: Multiple Sequence Predictions

Used this way, the method will produce for each sequence a secondary
prediction file. GOR is a single sequence with a relatively low
accuracy. It is possible to increase the accuracy by coupling BLAST and
GOR, this can be achieved with the following command:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

t\_coffee sample\_aln.fasta -template\_file PSISSP

.. raw:: html

   </div>

When doing so, the predictions for each sequence are obtained by
averaging the GOR predictions on every homologue as reported by a BLAST
against NR. By default the BLAST is done remotely at the NCBI using the
blastpgp web service of the EBI.

A similar output can be obtained for Transmembrane segment predictions:

 

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

t\_coffee sample\_aln.fasta -template\_file PSITM

.. raw:: html

   </div>

 

.. rubric:: Incorporation of the prediction in the alignment

It is possible to use the secondary prediction in order to reward the
alignment of similar elements

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

t\_coffee sample\_aln.fasta -template\_file PSISSP
-method\_evaluate\_mode ssp -method lalign\_id\_pair slow\_pair

.. raw:: html

   </div>

Likewise, it is possible to use this information with trans-membrane
domains

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

t\_coffee sample\_aln.fasta -template\_file PSITM
-method\_evaluate\_mode tm -method lalign\_id\_pair slow\_pair

.. raw:: html

   </div>

The overall effect is very crude and amounts to over-weighting by 30%
the score obtained when matching two residues in a similar secondary
structure state. The net consequence is that residues in similar
predicted states tend to be aligned more easily.

.. rubric:: 

 

.. rubric:: Using other secondary structure predictions

If you have your own predictions, you can use them. All you need is to
produce a template file where the file containing the secondary
structure prediction is declared along with the sequence:

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 31.0pt;
   background:#E6E6E6;margin-left:0cm;margin-right:1.0cm">

 

>hmgl\_wheat \_E\_ hmgl\_wheat.ssp

>hmgb\_chite \_E\_ hmgb\_chite.ssp

>hmgl\_trybr3 \_E\_ hmgl\_trybr3.ssp

>hmgl\_wheat2 \_E\_ hmgl\_wheat2.ssp

>hmgt\_mouse \_E\_ hmgt\_mouse.ssp

>hmgl\_trybr \_E\_ hmgl\_trybr.ssp

 

.. raw:: html

   </div>

where each template looks like this:

 

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 31.0pt;
   background:#E6E6E6;margin-left:0cm;margin-right:1.0cm">

 

>hmgl\_wheat

CCCCCCCCCCCCHHHHHHHCCCCCCCCCHHHHHHHHHHHHHHHCCCCHHHHHHHHHHHHHHHCE

 

.. raw:: html

   </div>

You can then run T-Coffee using your own template file

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

t\_coffee sample\_aln.fasta -template\_file <template\_file>
-method\_evaluate\_mode ssp -method lalign\_id\_pair slow\_pair

.. raw:: html

   </div>

 

.. rubric::     Output of the prediction

You can output a color coded version of your alignment using the
predicted structures

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

t\_coffee sample\_aln.fasta -template\_file PSISSP -output sec\_html

.. raw:: html

   </div>

A similar result can be obtained with trans-membrane regions:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

t\_coffee sample\_aln.fasta -template\_file PSITM -output tm\_html

.. raw:: html

   </div>

 

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 3.0pt;
   padding:1.0pt 4.0pt 1.0pt 4.0pt;background:#CCCCCC;margin-left:0cm;margin-right:
   4.25pt">

Combining Sequences and 3D-Structures

.. raw:: html

   </div>

Using structural information when aligning sequences is very useful. The
reason is that structures diverge slower than sequences. As a
consequence, one may still find a discernable homology between two
sequences that have been diverging for so long that their sequences have
evolved beyond recognition. Yet, when assembling the correct structure
based MSA, you will realize that these sequences contain key conserved
residues that a simple alignment procedure was unable to reveal. We show
you in this section how to make the best of T-Coffee tools to
incorporate structural information in your alignment.

.. rubric:: If you are in a Hurry: Expresso

.. rubric:: What is Expresso?

Expresso is the latest T-Coffee mode. It is not yet available for local
installation, but you can run it from the www.tcoffee.org server. The
principle of Expresso is simple: the server runs a BLAST between every
sequence in your query against the PDB database. If it finds a structure
similar enough to a sequence in your dataset (>60% identity), it will
use that structure as a template for your sequence.

Template files look something like:

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6;margin-left:0cm;margin-right:1.0cm">

>sp\|P08246\|ELNE\_HUMAN \_P\_ 1PPGE

>sp\|P20160\|CAP7\_HUMAN \_P\_ 1AE5

>sp\|P00757\|KLKB4\_MOUSE \_P\_ 1SGFX

>sp\|Q6H321\|KLK2\_HORSE \_P\_ 1GVZA

>sp\|P00773\|ELA1\_RAT \_P\_ 2D26C

>sp\|Q00871\|CTRB1\_PENVA \_P\_ 1AZZB

>sp\|P21844\|MCPT5\_MOUSE \_P\_ 1NN6A

>sp\|O35205\|GRAK\_MOUSE \_P\_ 1MZDA

>sp\|P07338\|CTRB1\_RAT \_P\_ 2CGAB

>sp\|P80015\|CAP7\_PIG \_P\_ 1FY3A

>sp\|P03953\|CFAD\_MOUSE \_P\_ 1FDPD

>sp\|Q7YRZ7\|GRAA\_BOVIN \_P\_ 1OP8F

>sp\|Q06606\|GRZ2\_RAT \_P\_ 1EUFA

>sp\|P08884\|GRAE\_MOUSE \_P\_ 1FI8B

.. raw:: html

   </div>

In a template file,  \ **\_P\_** indicates that the template is of type
structure (P for PDB). Template files can be generated manually or
automatically by the Expresso server. Whenever possible t\_coffee will
then align your sequences using the structural information contained in
the templates. If it encounters enough structures (as shown here) it
will produce a genuine structure based sequence alignment.

.. rubric:: Using Expresso

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee three\_pdb\_two\_seq.fasta -method
sap\_pair,slow\_pair -template\_file PDB

.. raw:: html

   </div>

.. rubric:: Aligning Sequences and Structures

.. rubric:: Mixing Sequences and Structures

Gather your sequences in the same file. Name your structures according
to their PDB identifier. The file three\_pdb\_two\_seq.fasta contains
five sequences, three are the sequences of PDB structures and two are
regular sequences.

What you want to do is to build a T-Coffee library where sequences with
a known structures are aligned with a structure alignment program (like
sap) while the other sequences are aligned using regular T-Coffee
methods. You can achieve this with the following command:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee three\_pdb\_two\_seq.fasta -method
sap\_pair,slow\_pair -template\_file PDB

.. raw:: html

   </div>

The option -template\_file is here to tell the program how to find the
PDB. In that case. EXPRESSO means that a remote BLAST (the EBI BLAST)
will be used to identify the best targets. If your sequences are already
named according to their PDB name, you can use:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

 PROMPT: t\_coffee three\_pdb\_two\_seq.fasta -method
sap\_pair,slow\_pair -template\_file \_SELF\_P\_

.. raw:: html

   </div>

 

\_SELF\_ means that the PDB identifier is the name of the sequences,
while \_P\_ is an indication that the template is indeed a PDB. These
indications are necessary for T-Coffee to fetch the relevant structures.

The good news is that you do not need to have PDB installed locally as
T-Coffee will automatically fetch the structures directly from RCSB (the
home of PDB). Of course, if your dataset only contains structures, your
alignment becomes a structural alignment.

If you have a fugue license, you can also add the fugue method to your
run. Fugue will align the structures with sequences whose structure is
unknown (this is called threading).

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee three\_pdb\_two\_seq.fasta -method
sap\_pair,slow\_pair,fugue\_pair -template\_file \_SELF\_P\_

.. raw:: html

   </div>

This can be written more concisely, using one of T-Coffee
special\_modes:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

   PROMPT: t\_coffee three\_pdb\_two\_seq.fasta -mode 3dcoffee

.. raw:: html

   </div>

or

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

   PROMPT: t\_coffee three\_pdb\_two\_seq.fasta -mode expresso

.. raw:: html

   </div>

.. rubric:: Using Sequences only

What often happens is that you have already built a dataset with
sequences that are very similar to PDB sequences but not exactly
identical. It may even be the case that the real sequence and the PDB
one do not match exactly because of some genetic engineering on the
structure. In this case, you have no structure whose sequence is exactly
similar to the sequences in your dataset. All you need to do is to
declare the equivalence sequences/structures and run T-Coffee, just like
Expresso does.

The first step is to fill up a template file that contains an explicit
declaration of the structures corresponding to your sequences. The
format is very simple and fasta-like. You can use the file:
sproteases\_small.template\_file

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6;margin-left:0cm;margin-right:1.0cm">

>sp\|P08246\|ELNE\_HUMAN \_P\_ 1PPGE

>sp\|P20160\|CAP7\_HUMAN \_P\_ 1AE5

>sp\|P00757\|KLKB4\_MOUSE \_P\_ 1SGFX

>sp\|Q6H321\|KLK2\_HORSE \_P\_ 1GVZA

.. raw:: html

   </div>

In this file, the first line is telling us that sequence
sp\|P08246\|ELNE\_HUMAN is associated with the structural template
1PPGE. The sequence and the structure do not need to be identical
although we recommend using structural templates more than 60% identical
with your actual sequences (i.e. similar enough so that they generate a
non ambiguous alignment). If your template file is ready, all you need
to do is run the following command.

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sproteases\_small.fasta -method slow\_pair,
lalign\_id\_pair, sap\_pair -template\_file
sproteases\_small.template\_file

.. raw:: html

   </div>

When you run this once, T-Coffee goes and fetches the structures. It
will then align them using sap. It takes a lot of time to fetch
structures, and it takes even more time to align them with sap. This is
why T-Coffee saves these important intermediate results in a special
location called the cache. By default, your cache is in
~/.t\_coffee/cache, it is a good idea to empty it from time to time…

.. rubric:: Aligning Profile using Structural Information

If you have two profiles to align, an ideal situation is when your
profiles each contain one or more structures. These structures will
guide the alignment of the profiles, even if they contain very distally
related sequences. We have prepared two such profiles (prf1\_pdb1.aln,
prf2\_pdb2.aln). You have two choices here. All you need is a template
file that declares which sequences have a known structure. If you only
want to align sequences, you can try:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -profile=profile1\_pdb1.aln, profile2\_pdb2.aln
-method sap\_pair -profile\_template\_file two\_profiles.template\_file

.. raw:: html

   </div>

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 3.0pt;
   padding:1.0pt 4.0pt 1.0pt 4.0pt;background:#CCCCCC;margin-left:0cm;margin-right:
   4.25pt">

How Good Is Your Alignment ?

.. raw:: html

   </div>

There are three strategies for evaluating your alignment. Structure is a
killer. If you have two structures available for your protein family,
you are in an ideal situation and you can use the iRMSD. If you don't,
you are left with the option of using sequence based methods like the
CORE index. These do pretty well in the CORE regions, but can be limited
in the loops. Another killer, less often at hand, is the use of
functional information. If you know some residues MUST be aligned
because they are functionally related, you can easily set up an
evaluation procedure using T-Coffee.

.. rubric:: Evaluating Alignments with The CORE index

.. rubric:: Computing the Local CORE Index

The CORE index is an estimation of the consistency between your
alignment and the computed library. The higher the consistency, the
better the alignment. The score reported with every T-Coffee alignment
is the concistency score. However, if you want to go further and
estimate the local concistency (known as the CORE index). Simply request
one extra output:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sproteases\_small.fasta -output=clustalw,html

.. raw:: html

   </div>

The format html leads to the output of a file named
sproteases\_small.html. Open this file. It contains a colorized version
of your alignment. In this colorized version, positions that have no
concistency with the library are in blue, a little in green, better
positions in yellow, then orange, then red. You can expect yellow
positions to be entirely correct.

.. rubric::   Computing the CORE index of any alignment

You can evaluate any existing alignment with the CORE index. All you
need to do is provide that alignment with the -infile flag and specify
that you want to evaluate it:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -infile=sproteases\_small.g10.cw\_aln -output=html
-score

.. raw:: html

   </div>

.. rubric::   Filtering Bad Residues

The local consistency score is between 0 and 9. If you need to build a
profile or identify a signature or do some phylogeny, it may be a good
idea to remove portions that are too unreliable. Here is how you can do
it.

You will first need to produce a score version of your alignment in
machine readable format. This format is called score\_ascii:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -infile=sproteases\_small.g10.cw\_aln
-output=score\_ascii

.. raw:: html

   </div>

You will then need to use seq\_reformat to filter out the portions you
are not interested in, using the file sproteases\_small.g10.score\_ascii
as a cache for filtering. For instance, use the following command to
only keep the residues whose score is between 5 and 9.

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in
sproteases\_small.g10.cw\_aln -struc\_in
sproteases\_small.g10.score\_ascii -struc\_in\_f number\_aln -action 
+keep '[5-9]'

.. raw:: html

   </div>

Not so neat... The reason is that most columns tend to be heterogenous
and contain a couple of unhappy residues. If you would rather generate
nice blocks, filter according to the consencus. This is easy and simply
requires adding the +use\_cons filter

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in
sproteases\_small.g10.cw\_aln -struc\_in
sproteases\_small.g10.score\_ascii -struc\_in\_f number\_aln -action 
+use\_cons +keep '[5-9]'

.. raw:: html

   </div>

Removing columns of gaps is just as easy. You simply need to add the
switch +rm\_gap

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in
sproteases\_small.g10.cw\_aln -struc\_in
sproteases\_small.g10.score\_ascii -struc\_in\_f number\_aln -action 
+use\_cons +keep '[5-9]' +rm\_gap

.. raw:: html

   </div>

.. rubric::   Filtering Gap Columns

You may want to remove columns that contain too many gaps. It is just a
small variation around the rm\_gap switch:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat -in
sproteases\_small.g10.cw\_aln -struc\_in
sproteases\_small.g10.score\_ascii -struc\_in\_f number\_aln -action
+rm\_gap 50

.. raw:: html

   </div>

Will remove all the columns containing 40% or more gaps.

 

.. rubric:: Evaluating an Alignment Using Structural Information: APDB
   and iRMSD

.. rubric:: What is the iRMSD?

APDB and the iRMSD are two closely related measures meant to evaluate
the accuracy of a sequence alignment without using a structure based
reference alignment. The iRMSD is a follow up of the APDB measure and we
now recommend using the iRMSD rather than APDB.

Although it may seem that the iRMSD was an attempt to get free iPODs
from Apple, it is not (or at least we never got the iPODs). The iRMSD is
a special RMSD (It stands for intra-catener) where the alignments are
evaluated using the structural information of the sequences with known
structures.

The strength of the iRMSD is its independence from a specific
superposition models.  When using the iRMSD to evaluate the score of a
sequence alignment, one does not need to superpose the two structures
and deduce a sequence alignment that will then be compared with the
target alignment. In practice, we use a Normalized version of the iRMSD,
the NiRMSD that makes it possible to compare alternative alignments of
different length. From a structural point of view, the iRMSD has a
meaning very similar to the iRMSD and it beahaves in a similar fashion
from a numerical point of view (similar ranges in Angstroms).

The first step of APDB is to measure the distances between the Ca of
each residue and its neighbors. Neighborhood is defined as a sphere of
**radius –maximum\_distance** (10Å by default). However, by setting
**–local\_mode** to “window”, the sphere can be replaced with a window
of 1/2 size '**-maximum\_distance'** residues.

|image0|

Given two aligned residues (X and Y on the Figure) the iRMSD measure is
an attempt to estimate the neighborhood support for the XY alignment.
This is done by measuring the difference of distances between X and Y
and every other pair of aligned residues within the same sphere (W and Z
on Figure 1). The iRMSD is obtained by measuring the average Root Mean
Square of these differences of distances. The lower the iRMSD, the
better the alignment. However, an alignment can obtain a good iRMSD by
simply having few aligned residues. To avoid this, the program also
reports the NiRMSD= MIN(L1,L2)\*iRMSD/Number Considered columns.

.. rubric:: How to Efficiently Use Structural Information

When it comes to evaluating Multiple Sequence Alignments, nothing beats
structural information. To use the methods we describe here, you will
need to have at least two structures, similar enough (>60%) to two
sequences in your dataset.

Here an outline of the best way to proceed:

1-              Make sure you include two structures whose sequences are
so distantly related that most of the other sequences are intermediates.

2-              Align your sequences without using the structural
information (i.e. t\_coffee, muscle...)

3-              Evaluate your alignment with iRMSD (see later  in this
section). The score will be S1

4-              Realign your sequences, but this time using structural
information (Expresso)

5-              Measure the score of that alignment (Score=S2)

If S1 and S2 are almost similar, it means your distantly related
structures were well aligned, and you can expect the intermediate
sequences to be well aligned as well. If S2 is much better than S1, you
can expect the structures to be well aligned in the second alignment,
while there is no guaranty that the alignment of the intermediate
sequences has improved as well, although in practice it often does

.. rubric:: Evaluating an Alignment With the iRMSD Package

Let us evaluate the alignment produced by Expresso, using the
template\_file returned by expresso:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg irmsd sproteases\_small.expresso
-template\_file sproteases\_small.template\_file

.. raw:: html

   </div>

This will deliver a long output. The most interesting bit is at the
bottom:

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6">

#TOTAL for the Full MSA

        TOTAL     EVALUATED:  52.90 %

        TOTAL     APDB:       81.59 %

        TOTAL     iRMSD:       0.71 Angs

        TOTAL     NiRMSD:      1.33 Angs

.. raw:: html

   </div>

APDB is an older measure, less robust than the iRMSD and it is an
attempt to estimate the fraction of pairs of residues whose alignment
seems to be correct form a structural point of view. The higher APDB,
the better the alignment, the lower the NiRMSD, the better the
alignment.

.. rubric:: Evaluating Alternative Alignments

The strength of structure based alignments is that they make it possible
to compare alternative alignments. In this case let us consider:

+--------------------------+--------------------------+--------------------------+
| Method                   | File                     | NiRMSD                   |
+--------------------------+--------------------------+--------------------------+
| Expresso                 | sproteases\_small.expres | 1.33 Å                   |
|                          | so                       |                          |
+--------------------------+--------------------------+--------------------------+
| T-Coffee                 | sproteases\_small.tc\_al | 1.35 Å                   |
|                          | n                        |                          |
+--------------------------+--------------------------+--------------------------+
| ClustalW                 | sproteases\_small.cw\_al | 1.52 Å                   |
|                          | n                        |                          |
+--------------------------+--------------------------+--------------------------+
| Mafft                    | sproteases\_small.mafft  | 1.36 Å                   |
+--------------------------+--------------------------+--------------------------+
| Muscle                   | sproteases\_small.muscle | 1.34 Å                   |
+--------------------------+--------------------------+--------------------------+

 

As expected, Expresso delivers the best alignment from a structural
point of view. This makes sense, since Expresso explicitely USES
structural information. The other figures show us that the structural
based alignment is only marginally better than most sequences based
alignments. Muscle seems to have a small edge here although the reality
is that all these figures are impossible to distinguish with the notable
exception of ClustalW

.. rubric:: Identifying the most distantly related sequences in your
   dataset

In order to identify the most distantly related sequences in a dataset,
you can use the seq\_reformat utility, in order to compare all the
sequences two by two and pick up the two having the lowest level of
identity:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -other\_pg seq\_reformat sproteases\_small.fasta
-output sim\_idscore \| grep TOP \|sort -rnk3

.. raw:: html

   </div>

This sim\_idscore indicates that every pair of sequences will need to be
aligned when estimating the similarity. The ouput (below) indicates that
the two sequences having the lowest level of identity are AEDAE and
MOUSE. It may not be a bad idea to choose these sequences (if possible)
for evaluating your MSA.

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6">

…

TOP        16   10       28.00  sp\|P29786\|TRY3\_AEDAE     
sp\|Q6H321\|KLK2\_HORSE   28.00

TOP        16    7       28.00  sp\|P29786\|TRY3\_AEDAE     
sp\|P08246\|ELNE\_HUMAN   28.00

TOP        16    1       28.00  sp\|P29786\|TRY3\_AEDAE     
sp\|P08884\|GRAE\_MOUSE   28.00

TOP        15   14       27.00    sp\|P80015\|CAP7\_PIG    
sp\|P00757\|KLKB4\_MOUSE   27.00

TOP        12    9       27.00  sp\|P20160\|CAP7\_HUMAN     
sp\|Q91VE3\|KLK7\_MOUSE   27.00

TOP         9    7       27.00  sp\|Q91VE3\|KLK7\_MOUSE     
sp\|P08246\|ELNE\_HUMAN   27.00

TOP        16    2       26.00  sp\|P29786\|TRY3\_AEDAE    
sp\|P21844\|MCPT5\_MOUSE   26.00

.. raw:: html

   </div>

.. rubric:: Evaluating an Alignment according to your own Criterion

.. rubric:: Establishing Your Own Criterion

Any kind of Feature can easily be turned into an evaluation grid. For
instance, the protease sequences we have been using here have a well
characterized binding site. A possible evaluation can be made as
follows. let us consider the Swissprot annotation of the two most
distantly related sequences. These two sequences contain the electron
relay system of the proteases. We can use it to build an evaluation
library: in P29786, the first Histidine is at position 68, while in
P21844 this Histidine is on position 66. We can therefore build a
library that will check whether these residues are properly aligned in
any MSA. The library will look like this:

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6">

! TC\_LIB\_FORMAT\_01

2

sp\|P21844\|MCPT5\_MOUSE 247
MHLLTLHLLLLLLGSSTKAGEIIGGTECIPHSRPYMAYLEIVTSENYLSACSGFLIRRNFVLTAAHCAGRSITVLLGAHNKTSKEDTWQKLEVEKQFLHPKYDENLVVHDIMLLKLKEKAKLTLGVGTLPLSANFNFIPPGRMCRAVGWGRTNV

NEPASDTLQEVKMRLQEPQACKHFTSFRHNSQLCVGNPKKMQNVYKGDSGGPLLCAGIAQGIASYVHRNAKPPAVFTRISHYRPWINKILREN

sp\|P29786\|TRY3\_AEDAE 254
MNQFLFVSFCALLDSAKVSAATLSSGRIVGGFQIDIAEVPHQVSLQRSGRHFCGGSIISPRWVLTRAHCTTNTDPAAYTIRAGSTDRTNGGIIVKVKSVIPHPQYNGDTYNYDFSLLELDESIGFSRSIEAIALPDASETVADGAMCTVSGWGDT

KNVFEMNTLLRAVNVPSYNQAECAAALVNVVPVTEQMICAGYAAGGKDSCQGDSGGPLVSGDKLVGVVSWGKGCALPNLPGVYARVSTVRQWIREVSEV

#1 2

    66  68     100

! SEQ\_1\_TO\_N

.. raw:: html

   </div>

You simply need to cut and paste this library in a file and use this
file as a library to measure the concistency between your alignment and
the correspondances declared in your library. The following command line
also makes it possible to visualy display the agreement between your
sequences and the library.

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -infile sproteases\_small.aln -lib 
charge\_relay\_lib.tc\_lib  -score -output html

.. raw:: html

   </div>

 

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 3.0pt;
   padding:1.0pt 4.0pt 1.0pt 4.0pt;background:#CCCCCC;margin-left:0cm;margin-right:
   4.25pt">

Integrating External Methods In T-Coffee

.. raw:: html

   </div>

The real power of T-Coffee is its ability to seamlessly combine many
methods into one. While we try to integrate as many methods as we can in
the default distribution, we do not have the means to be exhaustive and
if you desperately need your favourite method to be integrated, you will
need to bite the bullet …

.. rubric:: What Are The Methods Already Integrated in T-Coffee

Although, it does not necessarily do so explicitly, T-Coffee always end
up combining libraries. Libraries are collections of pairs of residues.
Given a set of libraries, T-Coffee makes an attempt to assemble the
alignment with the highest level of consistence. You can think of the
alignment as a timetable. Each library pair would be a request from
students or teachers, and the job of T-Coffee would be to assemble the
time table that makes as many people as possible happy…

In T-Coffee, methods replace the students/professors as constraints
generators. These methods can be any standard/non standard alignment
methods that can be used to generate alignments (pairwise, most of the
time). These alignments can be viewed as collections of constraints that
must be fit within the final alignment. Of course, the constraints do
not have to agree with one another…

This section shows you what are the vailable method in T-Coffee, and how
you can add your own methods, either through direct parameterization or
via a perl script. There are two kinds of methods: the internal and the
external. For the internal methods, you simply need to have T-Coffee up
and running. The external methods will require you to instal a package.

.. rubric:: List of INTERNAL Methods

Built in methods methods can be requested using the following names. To

**proba\_pair      **\ Adpated from Probcons, this method **[the current
default]** uses a pair HMM to compute a pairwise alignment with a
bi-phasic gap penalty. ****

**fast\_pair          **\ Makes a global fasta style pairwise alignment.
For proteins, matrix=blosum62mt, gep=-1, gop=-10, ktup=2. For DNA,
matrix=idmat (id=10), gep=-1, gop=-20, ktup=5. Each pair of residue is
given a score function of the weighting mode defined by -weight.

**slow\_pair         **\ Identical to fast pair, but does a full dynamic
programming, using the myers and miller algorithm. This method is
recommended if your sequences are distantly related.

ifast\_pair

**islow\_pair       
** Makes a global fasta alignmnet using the previously computed pairs as
a library. \`i\` stands for iterative. Each pair of residue is given a
score function of the weighting mode defined by -weight. The Library
used for the computation is the one computed before the method is used.
The resullt is therefore dependant on the order in methods and library
are set via the –in flag.

**align\_pdb\_pair**\     Uses the align\_pdb routine to align two
structures. The pairwise scores are those returnes by the align\_pdb
program. If a structure is missing, fast\_pair is used instead. Each
pair of residue is given a score function defined by align\_pdb.
[UNSUPORTED]

**lalign\_id\_pair  **\ Uses the ten top non intersecting local
alignments, as delivered by lalign. Each alignement is weighted with its
average percent identity.

**lalign\_rs\_s\_pair**\    Same as above but does also does self
comparison and uses the lalign raw\_score (s stands for self). This is
needed when extracting repeats.

**Matrix             **\ Amy matrix can be requested. Simply indicate as
a method the name of the matrix preceded with an X (i.e. Xpam250mt). If
you indicate such a matrix, all the other methods will simply be
ignored, and a standard fast progressive alignment will be computed. If
you want to change the substitution matrix used by the methods, use the
–matrix flag.

**cdna\_fast\_pair**\  This method computes the pairwise alignment of
two cDNA sequences. It is a fast\_pair alignment that only takes into
account the amino-acid similarity and uses different penalties for
amino-acid insertions and frameshifts. This alignment is turned into a
library where matched nucleotides receive a score equql to the average
level of identity at the amino-acid level. This mode is intended to
clean cDNA obtained from ESTs, or to align pseudo-genes.

WARNING: This method is currently unsuported.

.. rubric:: PLUG-INs:List OF EXTERNAL METHODS

.. rubric:: Plug-In: Using Methods Integrated in T-Coffee 

 

 

The following methods are external. They correspond to packages
developped by other groups that you may want to run within T-Coffee. We
are very open to extending these options and we welcome any request to
ad an extra interface. The following table lists the methods that can be
used as plug-ins:

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6">

Package              Where From

==========================================================

ClustalW             can interact with t\_coffee

----------------------------------------------------------

Poa                  http://www.bioinformatics.ucla.edu/poa/

----------------------------------------------------------

Muscle        http://www.bioinformatics.ucla.edu/poa/

----------------------------------------------------------

ProbCons             http://probcons.stanford.edu/

----------------------------------------------------------

MAFFT 
\ `http://www.biophys.kyoto- <http://www.biophys.kyoto-/>`__\ u.ac.jp/~katoh/programs/align/mafft/

----------------------------------------------------------

Dialign-T            \ http://dialign-t.gobics.de/

----------------------------------------------------------

PCMA                
\ `ftp://iole.swmed.edu/pub/PCMA/ <file://localhost/pub/PCMA>`__

----------------------------------------------------------

sap                  structure/structure comparisons

(obtain it from W. Taylor, NIMR-MRC).

---------------------------------------------------

Blast               
\ `www.ncbi.nih.nlm.gov <http://www.ncbi.nih.nlm.gov/>`__

---------------------------------------------------

Fugue                protein to structure alignment program

                     http://www-cryst.bioc.cam.ac.uk/fugue/download.html

 

.. raw:: html

   </div>

 

Once installed, most of these methods can be used as either pairwise or
multiple alignment methods. Note that all these methods use Blosum62 as
a default.

**clustalw\_pair   **\ Uses clustalw (default parameters) to align two
sequences. Each pair of residue is given a score function of the
weighting mode defined by -weight.

**clustalw\_msa   **\ Makes a multiple alignment using ClustalW and adds
it to the library. Each pair of residue is given a score function of the
weighting mode defined by -weight.

**probcons\_pair **\ Probcons package: probcons.stanford.edu/.

probcons\_msa idem.

**muscle\_pair     **\ Muscle package www.drive5.com/muscle/ .

muscle\_msa     idem.

**mafft\_pair       
**\ www.biophys.kyoto-u.ac.jp/~katoh/programs/align/mafft/ .

mafft\_msa        idem.

**pcma\_msa       **\ pcma package\ ****

**pcma\_pair       **\ pcma package\ ****

**poa\_msa          **\ poa package\ ****

**poa\_pair          **\ poa package\ ****

**dialignt\_pair    **\ dialignt  package\ ****

**dialignt\_msa**\     pcma package

**sap\_pair          **\ Uses sap to align two structures. Each pair of
residue is given a score function defined by sap. You must have sap
installed on your system to use this method.

**fugue\_pair       **\ Uses a standard fugue installation to make a
sequence /structure alignment. Fugue installation must be standard. It
does not have to include all the fugue packages but only:

1- joy, melody, fugueali, sstruc, hbond

2-copy fugue/classdef.dat            /data/fugue/SUBST/classdef.dat

OR

Setenv MELODY\_CLASSDEF=<location>

Setenv MELODY\_SUBST=fugue/allmat.dat

 All the configuration files must be in the right location.

 

To request a method, see the -in or the -method flag. For instance, if
you wish to request the use of fast\_pair and lalign\_id\_pair (the
current default):

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -seq sample\_seq1.fasta -method
fast\_pair,lalign\_id\_pair

.. raw:: html

   </div>

.. rubric:: Modifying the parameters of Internal and External Methods

.. rubric:: Internal Methods

It is possible to modify on the fly the parameters of hard coded
methods:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sample\_seq1.fasta –method
slow\_pair@EP@MATRIX@pam250mt@GOP@-10@GEP@-1

.. raw:: html

   </div>

EP stands for Extra parameters. These parameters will superseed any
other parameters.

.. rubric:: External Methods

External methods receive a command line built with the information
provided via the parameter file (see next heading). It is possible to
produce such a parameter file and to modify it in order to modify the
commands passed to the methods.

The passed command is built as follows:

<EXECUTABLE><PARAM1><IN\_FLAG><seq\_file><PARAM2><OUT\_FLAG><outname><PARAM>

 

You should know what is the best place for squizing your extra
parameters. It will depend on the application, although PARAM2 is
usually a good guess. Now if you want, for instance to modify the gap
penalty of clustalw, you can try the following:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

 PROMPT: t\_coffee sample\_seq1.fasta -method
clustalw\_msa@EP@PARAM2@-GAPOPEN%e100%s-GAPEXT%e10

.. raw:: html

   </div>

@EP is here to indicate that you will pass an extra parameter

@PARAM1 is the name of this parameter

the next filed is the parameter itself, where:

%s replace spaces

%e replaces the equal sign

Of course, you must know the command line of the program you are trying
to modify (clustalw in this case).

 

.. rubric:: Integrating External Methods

If the method you need is not already included in T-Coffee, you will
need to integrate it yourself. We give you here some guidelines on how
to do so.

.. rubric:: Direct access to external methods

A special method exists in T-Coffee that can be used to invoke any
existing program:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sample\_seq1.fasta –method=em@clustalw@pairwise

.. raw:: html

   </div>

In this context, Clustalw is a method that can be ran with the following
command line:

method –infile=<infile> -outfile=<outfile>

Clustalw can be replaced with any method using a similar syntax. If the
program you want to use cannot be run this way, you can either write a
perl wrapper that fits the bill or write a tc\_method file adapted to
your program (cf next section).

This special method (em, external method) uses the following syntax:

em@<method>@<aln\_mode:pairwise¦s\_pairwise\|multiple>

.. rubric:: Customizing an external method (with parameters) for
   T-Coffee

T-Coffee can run external methods, using a *tc\_method* file that can be
used in place of an established method. Two such files are incorporated
in T-Coffee. You can dump them and customize them according to your
needs:

For instance if you have ClustalW installed, you can use the following
file to run the

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee –other\_pg unpack\_clustalw\_method.tc\_method

PROMPT: t\_coffee –other\_pg unpack\_generic\_method.tc\_method

.. raw:: html

   </div>

The second file (generic\_method.tc\_method) contains many hints on how
to customize your new method. The first file is a very straightforward
example on how to have t\_coffee to run Clustalw with a set of
parameters you may be interested in:

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6">

 

\*TC\_METHOD\_FORMAT\_01

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*clustalw\_method.tc\_method\*\*\*\*\*\*\*\*\*

EXECUTABLE    clustalw

ALN\_MODE             pairwise

IN\_FLAG              -INFILE=

OUT\_FLAG             -OUTFILE=

OUT\_MODE             aln

PARAM         -gapopen=-10

SEQ\_TYPE             S

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

.. raw:: html

   </div>

This configuration file will cause T-Coffee to emit the following system
call:

clustalw –INFILE=tmpfile1 –OUTFILE=tmpfile2 –gapopen=-10

Note that ALN\_MODE instructs t\_coffee to run clustalw on every pair of
sequences (cf generic\_method.tc\_method for more details).

The tc\_method files are treated like any standard established method in
T-Coffee. For instance, if the file *clustalw\_method.tc\_method* is in
your current directory, run:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sample\_seq1.fasta –method clustalw\_method.tc\_method

.. raw:: html

   </div>

.. rubric:: Managing a collection of method files

It may be convenient to store all the method files in a single location
on your system. By default, t\_coffee will go looking into the directory
~/.t\_coffee/methods/. You can change this by either modifying the
METHODS\_4\_TCOFFEE in define\_headers.h (and recompile) or by modifying
the envoronement variable METHODS\_4\_TCOFFEE.

.. rubric:: Advanced Method Integration

It may sometimes be difficult to customize the program you want to use
through a tc\_method file. In that case, you may rather use an external
perl\_script to run your external application. This can easily be
achieved using the generic\_method.tc\_method file.

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6">

\*TC\_METHOD\_FORMAT\_01

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*generic\_method.tc\_method\*\*\*\*\*\*\*\*\*

EXECUTABLE    tc\_generic\_method.pl

ALN\_MODE             pairwise

IN\_FLAG              -infile=

OUT\_FLAG             -outfile=

OUT\_MODE             aln

PARAM         -method clustalw

PARAM         -gapopen=-10

SEQ\_TYPE             S

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

\* Note: &bsnp can be used to for  white spaces

.. raw:: html

   </div>

When you run this method:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee –other\_pg unpack\_generic\_method.tc\_method

PROMPT: t\_coffee sample\_seq1.fasta –method generic\_method.tc\_method

.. raw:: html

   </div>

T-Coffee runs the script tc\_generic\_method.pl on your data. It also
provides the script with parameters. In this case –method clustalw
indicates that the script should run clustalw on your data. The script
tc\_generic\_method.pl is incorporated in t\_coffee. Over the time, this
script will be the place where novel methods will be integrated

 will be used to run the script *tc\_generic\_method.pl*. The file
tc\_generic\_method.pl is a perl file, automatically generated by
t\_coffee. Over the time this file will make it possible to run all
available methods. You can dump the script using the following command:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee –other\_pg=unpack\_tc\_generic\_method.pl

.. raw:: html

   </div>

Note: If there is a copy of that script in your local directory, that
copy will be used in place of the internal copy of T-Coffee.

.. rubric:: The Mother of All method files…

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6">

\*TC\_METHOD\_FORMAT\_01

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*generic\_method.tc\_method\*\*\*\*\*\*\*\*\*\*\*\*\*

\*

\*       Incorporating new methods in T-Coffee

\*       Cedric Notredame 17/04/05

\*

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

\*This file is a method file

\*Copy it and adapt it to your need so that the method

\*you want to use can be incorporated within T-Coffee

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

\*                  USAGE                              \*

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

\*This file is passed to t\_coffee via –in:

\*

\*      t\_coffee –in Mgeneric\_method.method

\*

\*      The method is passed to the shell using the following

\*call:

\*<EXECUTABLE><IN\_FLAG><seq\_file><OUT\_FLAG><outname><PARAM>

\*

\*Conventions:

\*<FLAG\_NAME> <TYPE>        <VALUE>

\*<VALUE>:     no\_name       <=> Replaced with a space

\*<VALUE>:     &nbsp  <=> Replaced with a space

\*

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

\*                  EXECUTABLE                         \*

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

\*name of the executable

\*passed to the shell: executable

\*     

EXECUTABLE    tc\_generic\_method.pl

\*

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

\*                  ALN\_MODE                           \*

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

\*pairwise   ->all Vs all (no self )[(n2-n)/2aln]

\*m\_pairwise ->all Vs all (no self)[n^2-n]^2

\*s\_pairwise ->all Vs all (self): [n^2-n]/2 + n

\*multiple   ->All the sequences in one go

\*

ALN\_MODE             pairwise

\*

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

\*                  OUT\_MODE                           \*

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

\* mode for the output:

\*External methods:

\* aln -> alignmnent File (Fasta or ClustalW Format)

\* lib-> Library file (TC\_LIB\_FORMAT\_01)

\*Internal Methods:

\* fL -> Internal Function returning a Lib (Librairie)

\* fA -> Internal Function returning an Alignmnent

\*

OUT\_MODE             aln

\*

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

\*                  IN\_FLAG                             \*

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

\*IN\_FLAG

\*flag indicating the name of the in coming sequences

\*IN\_FLAG S no\_name ->no flag

\*IN\_FLAG S &nbsp–in&nbsp -> “ –in “

\*

IN\_FLAG              -infile=

\*

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

\*                  OUT\_FLAG                           \*

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

\*OUT\_FLAG

\*flag indicating the name of the out-coming data

\*same conventions as IN\_FLAG

\*OUT\_FLAG     S no\_name ->no flag

\*

OUT\_FLAG             -outfile=

\*

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

\*                  SEQ\_TYPE                           \*

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

\*G: Genomic, S: Sequence, P: PDB, R: Profile

\*Examples:

\*SEQTYPE      S      sequences against sequences (default)

\*SEQTYPE      S\_P    sequence against structure

\*SEQTYPE      P\_P    structure against structure

\*SEQTYPE      PS     mix of sequences and structure    

\*

SEQ\_TYPE      S

\*

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

\*                  PARAM                              \*

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

\*Parameters sent to the EXECUTABLE

\*If there is more than 1 PARAM line, the lines are

\*concatenated

\*

PARAM  -method clustalw

PARAM   -OUTORDER=INPUT -NEWTREE=core -align -gapopen=-15

\*

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

\*                  END                                \*

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

.. raw:: html

   </div>

.. rubric:: Weighting your Method

By default, the alignment produced by your method will be weighted
according to its percent identity. However, this can be customized via
the WEIGHT parameter.

The WEIGHT parameter supports all the values of the –weight flag. The
only difference is that the –weight value thus declared will only be
applied onto your method.

If needed you can also modify on the fly the WEIGHT value of your
method:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sample\_seq1.fasta –method slow\_pair@WEIGHT@OW2

.. raw:: html

   </div>

Will overweight by a factor 2 the weight of slow\_pair.

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sample\_seq1.fasta –method slow\_pair@WEIGHT@250

.. raw:: html

   </div>

Will cause every pair of slow\_pair to have a weight equal to 250

 

.. rubric:: Plug-Out: Using T-Coffee as a Plug-In

Just because it enjoys enslaving other methods as plug-ins, does not
mean that T-Coffee does not enjoy being incorporated within other
packages. We try to give as much support as possible to anyone who
wishes to incorportae T-Coffee in an alignment pipeline.

If you want to do so, please work out some way to incorporate T-Coffee
in your script . If you need some help along the ways, do not hesitate
to ask, as we will always be happy to either give assistance, or even
modify the package so that it accomodates as many needs as possible.

Once that procedure is over, set aside a couple of input files with the
correct parameterisation and send them to us. These will be included as
a distribution test, to insure that any further distribution remains
compliant with your application.

We currently support:

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6">

Package              Where From

==========================================================

Marna                www.bio.inf.unijena.de/Software/MARNA/download

----------------------------------------------------------

.. raw:: html

   </div>

.. rubric:: Creating Your Own T-Coffee Libraries

If the method you want to use is not integrated, or impossible to
integrate, you can generate your own libraries, either directly or by
turning existing alignments into libraries. You may also want to
precompute your libraries, in order to combine them at your convenience.

.. rubric:: Using Pre-Computed Alignments

If the method you wish to use is not supported, or if you simply have
the alignments, the simplest thing to do is to generate yourself the
pairwise/multiple alignments, in FASTA, ClustalW, msf or Pir format and
feed them into t\_coffee using the *-in* flag:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee –aln=sample\_aln1\_1.aln,sample\_aln1\_2.aln
–outfile=combined\_aln.aln

.. raw:: html

   </div>

.. rubric:: Customizing the Weighting Scheme

The previous integration method forces you to use the same weighting
scheme for each alignment and the rest of the libraries generated on the
fly. This weighting scheme is based on global pairwise sequence
identity. If you want to use a more specific weighting scheme with a
given method, you should either:

generate your own library (cf next section)

convert your aln into a lib, using the –weight flag:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee –aln sample\_aln1.aln –out\_lib=test\_lib.tc\_lib
–lib\_only –weight=sim\_pam250mt

PROMPT: t\_coffee –aln sample\_aln1.aln -lib test\_lib.tc\_lib
–outfile=outaln

PROMPT: t\_coffee –aln=sample\_aln1\_1.aln,sample\_aln1\_2.aln -method=
fast\_pair,lalign\_id\_pair –outfile=out\_aln

.. raw:: html

   </div>

.. rubric:: Generating Your Own Libraries

This is suitable if you have local alignments, or very detailed
information about your potential residue pairs, or if you want to use a
very specific weighting scheme. You will need to generate your own
libraries, using the format described in the last section.

You may also want to pre-compute your libraries in order to save them
for further use. For instance, in the following example, we generate the
local and the global libraries and later re-use them for combination
into a multiple alignment.

 

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sample\_seq1.fasta –method slow\_pair –out\_lib
slow\_pair\_seq1.tc\_lib –lib\_only

PROMPT: t\_coffee sample\_seq1.fasta –method lalign\_id\_pair –out\_lib
lalign\_id\_pair\_seq1.tc\_lib –lib\_only

.. raw:: html

   </div>

 

Once these libraries have been computed, you can then combine tem at
your convenience in a single MSA. Of course you can decide to only use
the local or the global library

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sample\_seq1.fasta –lib
lalign\_id\_pair\_seq1.tc\_lib, slow\_pair\_seq1.tc\_lib

.. raw:: html

   </div>

 

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 3.0pt;
   padding:1.0pt 4.0pt 1.0pt 4.0pt;background:#CCCCCC;margin-left:0cm;margin-right:
   4.25pt">

Frequently Asked Questions

.. raw:: html

   </div>

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.5pt;
   padding:0cm 0cm 0cm 0cm;background:#F7CBB7">

IMPORTANT: All the files mentionned here (sample\_seq...) can be found
in the example directory of the distribution.

.. raw:: html

   </div>

.. rubric:: Abnormal Terminations and Wrong Results

Q: The program keeps crashing when I give my sequences

A: This may be a format problem. Try to reformat your sequences using
any utility (readseq...). We recommend the Fasta format. If the problem
persists, contact us.

A: Your sequences may not be recognized for what they really are.
Normally T-Coffee recognizes the type of your sequences automatically,
but if it fails, use:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sample\_seq1.fasta -type=PROTEIN

.. raw:: html

   </div>

A: Costly computation or data gathered over the net is stored by
T-Coffee in a cache directory. Sometimes, some of these files can be
corrupted and cause an abnormal termination. You can either empty the
cache ( ~/.t\_coffee/cache/) or request T-Coffee to run without using
it:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee –pdb=struc1.pdb,struc2.pdb,struc3.pdb -method
sap\_pair –cache=no

.. raw:: html

   </div>

If you do not want to empty your cache, you may also use –cache=update
that will only update the files corresponding to your data

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee –pdb=struc1.pdb,struc2.pdb,struc3.pdb -method
sap\_pair –cache=update

 

.. raw:: html

   </div>

Q: The default alignment is not good enough

A: see next question

Q: The alignment contains obvious mistakes

A: This happens with most multiple alignment procedures. However, wrong
alignments are sometimes caused by bugs or an implementation mistake.
Please report the most unexpected results to the authors.

Q: The program is crashing

A: If you get the message:

FAILED TO ALLOCATE REQUIRED MEMORY

See the next question.

If the program crashes for some other reason, please check whether you
are using the right syntax and if the problem persists get in touch with
the authors.

Q: I am running out of memory

A: You can use a more accurate, slower and less memory hungry dynamic
programming mode called myers\_miller\_pair\_wise. Simply indicate the
flag:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sample\_seq1.fasta –special\_mode low\_memory

.. raw:: html

   </div>

Note that this mode will be much less time efficient than the default,
although it may be slightly more accurate. In practice the
parameterization associate with special mode turns off every memory
expensive heuristic within T-Coffee. For version 2.11 this amounts to

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee  sample\_seq1.fasta
-method=slow\_pair,lalign\_id\_pair –distance\_matrix\_mode=idscore
-dp\_mode=myers\_miller\_pair\_wise

.. raw:: html

   </div>

If you keep running out of memory, you may also want to lower –maxnseq,
to ensure that t\_coffee\_dpa will be used.

.. rubric:: Input/Output Control

Q: How many Sequences can t\_coffee handle

A: T-Coffee is limited to a maximum of 50 sequences. Above this number,
the program automatically switches to a heuristic mode, named DPA, where
DPA stands for Double Progressive Alignment.

DPA is still in development and the version currently shipped with
T-Coffee is only a beta version.

Q: Can I prevent the Output of all the warnings?

A: Yes, by setting  \ ***–no\_warning***

Q: How many ways to pass parameters to t\_coffee?

A: See the section well behaved parameters

Q: How can I change the default output format?

A: See the -output option, common output formats are:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sample\_seq1.fasta -output=msf,fasta\_aln

.. raw:: html

   </div>

Q: My sequences are slightly different between all the alignments.

A: It does not matter. T-Coffee will reconstruct a set of sequences that
incorporates all the residues potentially missing in some of the
sequences ( see flag -in).

Q: Is it possible to pipe stuff OUT of t\_coffee?

A: Specify stderr or stdout as output filename, the output will be
redirected accordingly. For instance

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sample\_seq1.fasta -outfile=stdout -out\_lib=stdout

.. raw:: html

   </div>

This instruction will output the tree (in new hampshire format) and the
alignment to stdout.

Q: Is it possible to pipe stuff INTO t\_coffee?

A: If as a file name, you specify stdin, the content of this file will
be expected throught pipe:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: cat sample\_seq1.fasta \| t\_coffee -infile=stdin

.. raw:: html

   </div>

will be equivalent to

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sample\_seq1.fasta

.. raw:: html

   </div>

If you do not give any argument to t\_coffee, they will be expected to
come from pipe:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: cat sample\_param\_file.param  \| t\_coffee -parameters=stdin

.. raw:: html

   </div>

For instance:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: echo –seq=sample\_seq1.fasta -method=clustalw\_pair \| t\_coffee
–parameters=stdin

.. raw:: html

   </div>

Q: Can I read my parameters from a file?

A: See the well behaved parameters section.

Q: I want to  decide myself on the name of the output files!!!

A: Use the *-run\_name* flag.

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sample\_seq1.fasta –run\_name=guacamole

.. raw:: html

   </div>

Q: I want to use the sequences in an alignment file

A: Simply fed your alignment, any way you like, but do not forget to
append the prefix S for sequence:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee Ssample\_aln1.aln

PROMPT: t\_coffee -infile=Ssample\_aln1.aln

PROMPT: t\_coffee –seq=sample\_aln1.aln
-method=slow\_pair,lalign\_id\_pair –outfile=outaln

.. raw:: html

   </div>

This means that the gaps will be reset and that the alignment you
provide will not be considered as an alignment, but as a set of
sequences.

Q: I only want to produce a library

A: use the –lib\_only flag

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sample\_seq1.fasta -out\_lib=sample\_lib1.tc\_lib
-lib\_only

.. raw:: html

   </div>

Please, note that the previous usage supersedes the use of the –convert
flag. Its main advantage is to restrict computation time to the actual
library computation.

Q: I want to turn an alignment into a library

A: use the –lib\_only flag

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee –in=Asample\_aln1.aln -out\_lib=sample\_lib1.tc\_lib
-lib\_only

.. raw:: html

   </div>

It is also possible to control the weight associated with this alignment
(see the –weight section).

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee –aln=sample\_aln1.aln -out\_lib=sample\_lib1.tc\_lib
-lib\_only –weight=1000

.. raw:: html

   </div>

Q: I want to concatenate two libraries

A: You cannot concatenate these files on their own. You will have to use
t\_coffee. Assume you want to combine tc\_lib1.tc\_lib and
tc\_lib2.tc\_lib.

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -lib=sample\_lib1.tc\_lib,sample\_lib2.tc\_lib
–lib\_only -out\_lib=sample\_lib3.tc\_lib

.. raw:: html

   </div>

Q: What happens to the gaps when an alignment is fed to T-Coffee

A: An alignment is ALWAYS considered as a library AND a set of
sequences. If you want your alignment to be considered as a library
only, use the S identifier.

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee Ssample\_aln1.aln –outfile=outaln

.. raw:: html

   </div>

It will be seen as a sequence file, even if it has an alignment format
(gaps will be removed).

Q: I cannot print the html graphic display!!!

A: This is a problem that has to do with your browser. Instead of
requesting the score\_html output, request the score\_ps output that can
be read using ghostview:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sample\_seq1.fasta -output=score\_ps

.. raw:: html

   </div>

or     

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sample\_seq2.fasta -output=score\_pdf

.. raw:: html

   </div>

Q: I want to output an html file and a regular file

A: see the next question

Q: I would like to output more than one alignment format at the same
time

A: The flag -output accepts more than one parameter. For instance,

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sample\_seq1.fasta -output=clustalw,html,score\_ps,msf

.. raw:: html

   </div>

This will output founr alignment files in the corresponding formats.
Alignments' names will have the format name as an extension.

Note: you need to have the converter ps2pdf installed on your system
(standard under Linux and cygwin). The latest versions of Internet
Explorer and Netscape now allow the user to print the HTML display *Do
not forget to request Background printing.*

.. rubric:: Alignment Computation 

Q: Is T-Coffee the best? Why Not Using Muscle, or Mafft, or ProbCons???

A: All these packages are good packages and they sometimes outperform
T-Coffee. They also claim to outperform one another... If you have them
installed locally, you can have T-Coffee to generate a consensus
alignment:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sample\_seq1.fasta –method muscle\_msa,probcons\_msa,
mafft\_msa, lalign\_id\_pair,slow\_pair

.. raw:: html

   </div>

Q: Can t\_coffee align Nucleic Acids ???

A: Normally it can, but check in the log that the program recognises the
right type ( In the INPUT SEQ section, Type: xxxx). If this fails, you
will need to manually set the type:

 

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sample\_dnaseq1.fasta –type dna

.. raw:: html

   </div>

Q: I do not want to compute the alignment.

A: use the -convert flag

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sample\_aln1.aln -convert -output=gcg

.. raw:: html

   </div>

This command will read the .aln file and turn it into an .msf alignment.

Q: I would like to force some residues to be aligned.

If you want to brutally force some residues to be aligned, you may use
as a post processing, the force\_aln function of seq\_reformat:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee –other\_pg seq\_reformat –in sample\_aln4.aln –action
+force\_aln seq1 10 seq2 15

PROMPT: t\_coffee –other\_pg seq\_reformat –in sample\_aln4.aln –action
+force\_aln sample\_lib4.tc\_lib02

.. raw:: html

   </div>

sample\_lib4.tc\_lib02 is a T-Coffee library using the tc\_lib02 format:

\*TC\_LIB\_FORMAT\_02

SeqX resY ResY\_index      SeqZ ResZ ResZ\_index

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.5pt;
   padding:0cm 0cm 0cm 0cm;background:#F7CBB7">

The TC\_LIB\_FORMAT\_02 is still experimental and unsupported. It can
only be used in the context of the force\_aln function described here.

.. raw:: html

   </div>

Given more than one constraint, these will be applied one after the
other, in the order they are provided. This greedy procedure means that
the Nth constraint may disrupt the (N-1)th previously imposed
constraint, hence the importance of forcing the constraints in the right
order, with the most important coming last.

We do not recommend imposing hard constraints on an alignment, and it is
much more advisable to use the soft constraints provided by standard
t\_coffee libraries (cf. building your own libraries section)

Q: I would like to use structural alignments.

See the section Using structures in Multiple Sequence Alignments, or see
the question *I want to build my own libraries.*

Q: I want to build my own libraries.

A: Turn your alignment into a library, forcing the residues to have a
very good weight, using structure:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee –aln=sample\_seq1.aln -weight=1000
-out\_lib=sample\_seq1.tc\_lib –lib\_only

.. raw:: html

   </div>

The value 1000 is simply a high value that should make it more likely
for the substitution found in your alignment to reoccur in the final
alignment. This will produce the library *sample\_aln1.tc\_lib* that you
can later use when aligning all the sequences:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee –seq=sample\_seq1.fasta -lib=sample\_seq1.tc\_lib
–outfile sample\_seq1.aln

.. raw:: html

   </div>

If you only want some of these residues to be aligned, or want to give
them individual weights, you will have to edit the library file yourself
or use the –force\_aln option (cf FAQ: I would like to force some
residues to be aligned). A value of N\*N \* 1000 (N being the number of
sequences) usually ensure the respect of a constraint.

Q: I want to use my own tree

A: Use the -usetree=<your own tree> flag.

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sample\_seq1.fasta –usetree=sample\_tree.dnd

.. raw:: html

   </div>

Q: I want to align coding DNA

A: use the fasta\_cdna\_pair method that compares two cDNA using the
best reading frame and taking frameshifts into account.

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee three\_cdna.fasta –method=cdna\_fast\_pair

.. raw:: html

   </div>

Notice that in the resulting alignments, all the gaps are of modulo3,
except one small gap in the first line of sequence hmgl\_trybr. This is
a framshift, made on purpose. You can realign the same sequences while
ignoring their coding potential and treating them like standard DNA:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee three\_cdna.fasta

.. raw:: html

   </div>

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.5pt;
   padding:0cm 0cm 0cm 0cm;background:#F7CBB7">

Note: This method has not yet been fully tested and is only provided
“as-is” with no warranty. Any feedback will be much appreciated.

.. raw:: html

   </div>

Q: I do not want to use all the possible pairs when computing the
library

Q: I only want to use specific pairs to compute the library

A: Simply write in a file the list of sequence groups you want to use:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sample\_seq1.fasta
–method=clustalw\_pair,clustalw\_msa –lib\_list=sample\_list1.lib\_list
–outfile=test

.. raw:: html

   </div>

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6">

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*sample\_list1.lib\_list\*\*\*\*

2 hmgl\_trybr hmgt\_mouse

2 hmgl\_trybr hmgb\_chite

2 hmgl\_trybr hmgl\_wheat

3 hmgl\_trybr hmgl\_wheat hmgl\_mouse

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*sample\_list1.lib\_list\*\*\*\*

 

.. raw:: html

   </div>

Note: Pairwise methods (slow\_pair…) will only be applied to list of
pairs of sequences, while multiple methods (clustalw\_aln) will be
applied to any dataset having more than two sequences.

Q: There are duplicates or quasi-duplicates in my set

A: If you can remove them, this will make the program run faster,
otherwise, the t\_coffee scoring scheme should be able to avoid
over-weighting of over-represented sequences.            

.. rubric:: Using Structures and Profiles

Q: Can I align sequences to a profile with T-Coffee?

A: Yes, you simply need to indicate that your alignment is a profile
with the R tag..

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sample\_seq1.fasta -profile=sample\_aln2.aln –outfile
tacos

.. raw:: html

   </div>

Q: Can I align sequences Two or More Profiles?

A: Yes, you, simply tag your profiles with the letter R and the program
will treat them like standard sequences.

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee -profile=sample\_aln1.fasta,sample\_aln2.aln –outfile
tacos

.. raw:: html

   </div>

Q: Can I align two profiles according to the structures they contain?

A: Yes. As long as the structure sequences are named according to their
PDB identifier

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee  -profile=sample\_profile1.aln,sample\_profile2.aln
–special\_mode=3dcoffee –outfile=aligne\_prf.aln

.. raw:: html

   </div>

Q: T-Coffee becomes very slow when combining sequences and structures

A: This is true. By default the structures are feteched on the net,
using RCSB. The problem arises when T-Coffee looks for the structure of
sequences WITHOUT structures. One solution is to install PDB locally. In
that case you will need to set two environment variables:

          

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

setenv (or export)  PDB\_DIR=”directory containing the pdb structures”
 setenv (or export)  NO\_REMOTE\_PDB\_DIR=1

.. raw:: html

   </div>

Interestingly, the observation that sequences without structures are
those that take the most time to be checked is a reminder of the
strongest rational argument that I know of against torture: any innocent
would require the maximum amount of torture to establish his/her
innocence, which sounds...ahem...strange., and at least inneficient.
Then again I was never struck by the efficiency of the Bush
administration.

Q: Can I use a local installation of PDB?

A: Yes, T-Coffe supports three types of installations:

         -an add-hoc installation where all your structures are in a
directory, under the form pdbid.pdb or pdbid.id.Z or pdbid.pdb.gz. In
that case, all you need to do is set the environement variables
correctly:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

setenv (or export)  PDB\_DIR=”directory containing the pdb structures”
 setenv (or export)  NO\_REMOTE\_PDB\_DIR=1

.. raw:: html

   </div>

-A standard pdb installation using the all section of pdb. In that case,
you must set the variables to:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

setenv (or export)  PDB\_DIR=”<some absolute
path>/data/structures/all/pdb/”
 setenv (or export)  NO\_REMOTE\_PDB\_DIR=1

.. raw:: html

   </div>

-A standard pdb installation using the divided section of pdb:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

setenv (or export)  PDB\_DIR=”<some absolute
path>/data/structures/divided/pdb/”
 setenv (or export)  NO\_REMOTE\_PDB\_DIR=1

.. raw:: html

   </div>

If you need to do more clever things, you should know that all the PDB
manipulation is made in T-Coffee by a perl script named
extract\_from\_pdb. You can extract this script from T-Coffee:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

t\_coffee –other\_pg unpack\_extract\_from\_pdb
 chmod u+x extract\_from\_pdb

.. raw:: html

   </div>

You can then edit the script to suit your needs. T-Coffee will use your
edited version if it is in the current directory. It will issue a
warning that it used a local version.

If you make extensive modifications, I would appreciate you send me the
corrected file so that I can incorporate it in the next distribution.

 

.. rubric:: Alignment Evaluation

Q: How good is my alignment?

A: see what is the color index?

Q: What is that color index?

A: T-Coffee can provide you with a measure of consistency among all the
methods used. You can produce such an output using:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sample\_seq1.fasta -output=html

.. raw:: html

   </div>

This will compute your\_seq.score\_html that you can view using
netscape. An alternative is to use score\_ps or score\_pdf that can be
viewed using ghostview or acroread, score\_ascii will give you an
alignment that can be parsed as a text file.

A book chapter describing the CORE index is available on:

http://www.tcoffee.org/Publications/Pdf/core.pp.pdf

 

Q: Can I evaluate alignments NOT produced with T-Coffee?

A: Yes. You may have an alignment produced from any source you like. To
evaluate it do:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee –infile=sample\_aln1.aln -lib=sample\_aln1.tc\_lib
–special\_mode=evaluate

.. raw:: html

   </div>

If you have no library available, the library will be computed on the
fly using the following command. This can take some time, depending on
your sample size. To monitor the progress in a situation where the
default library is being built, use:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee –infile=sample\_aln1.aln –special\_mode evaluate

.. raw:: html

   </div>

Q: Can I Compare Two Alignments?

A: Yes. You can treat one of your alignments as a library and compare it
with the second alignment:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee –infile=sample\_aln1\_1.aln -aln=sample\_aln1\_2.aln
–special\_mode=evaluate

.. raw:: html

   </div>

If you have no library available, the library will be computed on the
fly using the following command. This can take some time, depending on
your sample size. To monitor the progress in a situation where the
default library is being built, use:

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee –infile=sample\_aln1.aln –special\_mode evaluate

.. raw:: html

   </div>

Q: I am aligning sequences with long regions of very good overlap

A: Increase the ktuple size ( up to 4 or 5 for DNA) and up to 3 for
proteins.

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee sample\_seq1.fasta -ktuple=3

.. raw:: html

   </div>

This will speed up the program. It can be very useful, especially when
aligning ESTs.

Q: Why is T-Coffee changing the names of my sequences!!!!

A: If there is no duplicated name in your sequence set, T-Coffee's
handling of names is consistent with Clustalw, (Cf Sequence Name
Handling in the Format section). If your dataset contains sequences with
identical names, these will automatically be renamed to:

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.0pt;
   mso-border-alt:solid windowtext .5pt;padding:1.0pt 4.0pt 1.0pt 4.0pt;
   background:#E6E6E6">

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

>seq1

>seq1

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

>seq1

>seq1\_1

\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

.. raw:: html

   </div>

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 1.5pt;
   padding:0cm 0cm 0cm 0cm;background:#F7CBB7">

Warning: The behaviour is undefined when this creates two sequence with
a similar names.

.. raw:: html

   </div>

.. rubric:: Improving Your Alignment

Q: How Can I Edit my Alignment Manually?

A: Use jalview, a Java online MSA editor: www.jalview.org

Q: Have I Improved or Not my Alignment?

A: Using structural information is the only way to establish whether you
have improved or not your alignment. The CORE index can also give you
some information.

 

 

 

 

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 3.0pt;
   padding:1.0pt 4.0pt 1.0pt 4.0pt;background:#CCCCCC;margin-left:0cm;margin-right:
   4.25pt">

Addresses and Contacts

.. raw:: html

   </div>

.. rubric:: Contributors

T-coffee is developed, maintained, monitored, used and debugged by a
dedicated team that include:

            Cédric Notredame

Fabrice Armougom

Des Higgins

Sebastien Moretti

Orla O’Sullivan

Eamon O’Toole

Olivier Poirot

Karsten Suhre

Vladimir Keduas

Iain Wallace

.. rubric:: Addresses

We are always very eager to get some user feedback. Please do not
hesitate to drop us a line  at:
`cedric.notredame@europe.com <mailto:cedric.notredame@europe.com>`__ The
latest updates of T-Coffee are always available  on: www.tcoffee.org .
On this address you will also find a link to some of the online T-Coffee
servers, including Tcoffee@igs

 

T-Coffee can be used to automatically check if an updated version is
available, however the program will not update automatically, as this
can cause endless reproducibility problems.

.. raw:: html

   <div style="mso-element:para-border-div;border-top:solid black 1.0pt;
   border-left:none;border-bottom:solid black 1.0pt;border-right:none;mso-border-top-alt:
   solid black .25pt;mso-border-bottom-alt:solid black .25pt;padding:7.0pt 0cm 10.0pt 0cm;
   background:#FFFFCC;margin-left:2.0cm;margin-right:0cm">

PROMPT: t\_coffee –update

.. raw:: html

   </div>

 

.. raw:: html

   <div
   style="mso-element:para-border-div;border:solid windowtext 3.0pt;
   padding:1.0pt 4.0pt 1.0pt 4.0pt;background:#CCCCCC;margin-left:0cm;margin-right:
   4.25pt">

References

.. raw:: html

   </div>

It is important that you cite T-Coffee when you use it. Citing us is
(almost) like giving us money: it helps us convincing our institutions
that what we do is useful and that they should keep paying our salaries
and deliver Donuts to our offices from time to time (Not that they ever
did it, but it would be nice anyway).

 

Cite the server if you used it, otherwise, cite the original paper from
2000 (No, it was never named "T-Coffee 2000").

`Notredame C, Higgins DG, Heringa
J. <http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=pubmed&dopt=Abstract&list_uids=10964570>`__\ 

`Related
Articles, <http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=pubmed&cmd=Display&dopt=pubmed_pubmed&from_uid=10964570>`__
Links

| T-Coffee: A novel method for fast and accurate multiple sequence
  alignment.
|  J Mol Biol. 2000 Sep 8;302(1):205-17.
|  PMID: 10964570 [PubMed - indexed for MEDLINE]

Other useful publications include:

.. rubric:: T-Coffee

`Claude JB, Suhre K, Notredame C, Claverie JM, Abergel
C. <http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=pubmed&dopt=Abstract&list_uids=15215460>`__\ 

`Related
Articles, <http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=pubmed&cmd=Display&dopt=pubmed_pubmed&from_uid=15215460>`__
Links

| CaspR: a web server for automated molecular replacement using homology
  modelling.
|  Nucleic Acids Res. 2004 Jul 1;32(Web Server issue):W606-9.
|  PMID: 15215460 [PubMed - indexed for MEDLINE]

 

`Poirot O, Suhre K, Abergel C, O'Toole E, Notredame
C. <http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=pubmed&dopt=Abstract&list_uids=15215345>`__\ 

`Related
Articles, <http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=pubmed&cmd=Display&dopt=pubmed_pubmed&from_uid=15215345>`__
Links

| 3DCoffee@igs: a web server for combining sequences and structures into
  a multiple sequence alignment.
|  Nucleic Acids Res. 2004 Jul 1;32(Web Server issue):W37-40.
|  PMID: 15215345 [PubMed - indexed for MEDLINE]

 

`O'Sullivan O, Suhre K, Abergel C, Higgins DG, Notredame
C. <http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=pubmed&dopt=Abstract&list_uids=15201059>`__\ 

`Related
Articles, <http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=pubmed&cmd=Display&dopt=pubmed_pubmed&from_uid=15201059>`__
Links

| 3DCoffee: combining protein sequences and structures within multiple
  sequence alignments.
|  J Mol Biol. 2004 Jul 2;340(2):385-95.
|  PMID: 15201059 [PubMed - indexed for MEDLINE]

 

`Poirot O, O'Toole E, Notredame
C. <http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=pubmed&dopt=Abstract&list_uids=12824354>`__\ 

`Related
Articles, <http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=pubmed&cmd=Display&dopt=pubmed_pubmed&from_uid=12824354>`__
Links

| Tcoffee@igs: A web server for computing, evaluating and combining
  multiple sequence alignments.
|  Nucleic Acids Res. 2003 Jul 1;31(13):3503-6.
|  PMID: 12824354 [PubMed - indexed for MEDLINE]

 

`Notredame
C. <http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=pubmed&dopt=Abstract&list_uids=11301309>`__\ 

`Related
Articles, <http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=pubmed&cmd=Display&dopt=pubmed_pubmed&from_uid=11301309>`__
Links

| Mocca: semi-automatic method for domain hunting.
|  Bioinformatics. 2001 Apr;17(4):373-4.
|  PMID: 11301309 [PubMed - indexed for MEDLINE]

 

`Notredame C, Higgins DG, Heringa
J. <http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=pubmed&dopt=Abstract&list_uids=10964570>`__\ 

`Related
Articles, <http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=pubmed&cmd=Display&dopt=pubmed_pubmed&from_uid=10964570>`__
Links

| T-Coffee: A novel method for fast and accurate multiple sequence
  alignment.
|  J Mol Biol. 2000 Sep 8;302(1):205-17.
|  PMID: 10964570 [PubMed - indexed for MEDLINE]

 

`Notredame C, Holm L, Higgins
DG. <http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=pubmed&dopt=Abstract&list_uids=9682054>`__\ 

`Related
Articles, <http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=pubmed&cmd=Display&dopt=pubmed_pubmed&from_uid=9682054>`__
Links

| COFFEE: an objective function for multiple sequence alignments.
|  Bioinformatics. 1998 Jun;14(5):407-22.
|  PMID: 9682054 [PubMed - indexed for MEDLINE]

 

.. rubric:: Mocca

`Notredame
C. <http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=pubmed&dopt=Abstract&list_uids=11301309>`__\ 

`Related
Articles, <http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=pubmed&cmd=Display&dopt=pubmed_pubmed&from_uid=11301309&tool=ExternalSearch>`__
Links

| Mocca: semi-automatic method for domain hunting.
|  Bioinformatics. 2001 Apr;17(4):373-4.
|  PMID: 11301309 [PubMed - indexed for MEDLINE]

.. rubric:: CORE

\ `http://www.tcoffee.org/Publications/Pdf/core.pp.pdf <http://www.tcoffee.org/Publications/Pdf/core.pp.pdf>`__\ 

.. rubric:: Other Contributions

We do not mean to steal code, but we will always try to re-use
pre-existing code whenever that code exists, free of copyright, just
like we expect people to do with our code. However, whenever this
happens, we make a point at properly citing the source of the original
contribution. If ever you recognize a piece of your code improperly
cited, please drop us a note and we will be happy to correct that.

In the mean time, here are some important pieces of code from other
packages that have been incorporated within the T-Coffee package. These
include:

-TM-align package from Zhang, Jeffrey and Skolnik (NAR, 2005, 33:2303)

-The Sim algorithm of Huang and Miller that given two sequences computes
the N best scoring local alignments.

         -The tree reading/computing routines are taken from the
ClustalW Package, courtesy of Julie Thompson, Des Higgins and Toby
Gibson (Thompson, Higgins, Gibson, 1994, 4673-4680,vol. 22, Nucleic Acid
Research).

         -The implementation of the algorithm for aligning two sequences
in linear space was adapted from Myers and Miller, in CABIOS, 1988,
11-17, vol. 1)

         -Various techniques and algorithms have been implemented.
Whenever relevant, the source of the code/algorithm/idea is indicated in
the corresponding function.

         -64 Bits compliance was implemented by Benjamin Sohn,
Performance Computing Center Stuttgart (HLRS), Germany.

**An enormous thanks to these people who believe in free open source
code..**

.. rubric:: Bug Reports and Feedback

         -Prof David Jones (UCL) reported and corrected the PDB1K bug
(now t\_coffee/sap can align PDB sequences longer than 1000 AA).

         -Johan Leckner reported several bugs related to the treatment
of PDB structures, insuring a consistent behavior between version 1.37
and current ones.

 

 

 

.. raw:: html

   </div>

.. |image0| image:: t_coffee_tutorial_files/image002.gif
