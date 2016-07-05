#####################
T-Coffee Installation
#####################

******************************
T-Coffee Standard Installation
******************************

This section describes the installation procedures relevant for a standard use of T-Coffee on the most common operative systems. The procedure is quite straightforward, contact us if you have any problem.


Third party packages and on-demand installations
================================================
T-Coffee is a complex package that interacts with many other third part softwares. If you only want a stand-alone version of T-Coffee, you may install that package on its own with all its dependencies. If you want to use a most sophisticated flavor (3D-Coffee, Expresso, R-Cofeee, etc...), the installer will have to install all the third party packages required.

Note that since version 7.56, T-Coffee uses on-demand installation and install the third party packages it needs. This only works for packages not requiring specific licenses and that can be installed by the regular installer. Please let us know if you would like another third party package to be included.

Whenever on-demand installation or automated installation fails because of unforeseen system specificities, users should install the third party package manually. This documentation gives some tips we have found useful, but users are encouraged to send their feedbacks and share their experiences in order to improve this documentation.


Standard installation of T-Coffee
=================================

Automated installation
----------------------
We now recommend that you use the automated installer provided for Unix or Linux platforms. In any case, you first have to download the T-Coffee installer corresponding to your system from the T-Coffee website and then follow the adequate procedure.


Unix/Linux
^^^^^^^^^^
::

   rpm -i <rpm>


You need to have: gcc, g77, CPAN and an internet connection and your root password. You will also need the two following perl modules: LWP and XML::Simple. These are needed if you want to use the web services provided by the EBI via REST (http://www.ebi.ac.uk/Tools/webservices/tutorials/02_rest)


::

  1.      gunzip t_coffee.tar.gz

  2.      tar -xvf t_coffee.tar

  3.      cd t_coffee

  4.      ./install t_coffee



This installation will only install the stand alone T-Coffee. If you want to install a specific mode of T-Coffee, you may try the following commands that will try to gather all the necessary third party packages. Note that a package already found on your system will not be re-installed.


::

   ./install t_coffee

   ./install mcoffee

   ./install 3dcoffee

   ./install rcoffee

   ./install psicoffee


Or even:


::

   ./install all



All the corresponding executables will be downloaded automatically and installed in:


::

   $HOME/.t_coffee/plugins


- If you executables are in a different location, give it to T-Coffee using the -plugins flag.

- If the installation of any of the package fails, you should install it yourself using the provided link (see below) and following the authors instructions.

- If you have not managed to install SOAP::Lite, you can re-install it later (from anywhere) following steps 1-2.

- This procedure attempts 3 things: installing and compiling T-Coffee (C program), installing and compiling TMalign (Fortran), installing and compiling XML::Simple.

- If you have never installed a perl module before, CPAN will ask you many questions: say Yes to all.

- If everything went well, the procedure has created in the bin directory two executables: t_coffee and TMalign (Make sure these executables are on your $PATH)


Mac OS X
^^^^^^^^
Make sure you have the developer's kit installed (compilers and makefile) and follow the Unix procedure.

::

   Click on the .dmg file and follow the installation procedure


Microsoft Windows (Cygwin)
^^^^^^^^^^^^^^^^^^^^^^^^^^

T-Coffee doesn't run on Windows; in order to be able to run T-Coffee you need to install a Linux environment such as Cygwin:

::

  1.      Install Cygwin

  2.      Download the installer (NOT Cygwin/X)

  3.      Click on view to list ALL the packages

  4.      Select: gcc-core, make, wget

  5.      Optional: ssh, xemacs, nano

  6.      Run mkpasswd in Cygwin (as requested when you start Cygwin)

  7.      Install T-Coffee within Cygwin using the Unix procedure


CLUSTER Installation
^^^^^^^^^^^^^^^^^^^^
In order to run, T-Coffee must have a value for the http_proxy and for the e-mail. In order to do so you can either:

export the following values:

export http_proxy_4_TCOFFEE='proxy' or '' if no proxy

export EMAIL_4_TCOFFEE='your email'

OR

modify the file ~/.t_coffee/t_coffee_env

OR

add to your command line: t_coffee .... -proxy=<proxy> -email=<email>
(if you have no proxy: t_coffee ... -proxy -email=<email>)




*************************
Using BLAST With T-Coffee
*************************

BLAST is a program that searches sequence databases for homologues of a query sequence. It works for proteins and nucleic acids alike. In theory BLAST is just a package like any, but in practice things are a bit more complex. To run well, BLAST requires up-to-date databases (that can be fairly large, like n.r. or UniProt) and a powerful computer.

Fortunately, an increasing number of institutes or companies are now providing BLAST clients that run over the net. It means that all you need is a small program that send your query to the big server and gets the results back. This prevents you from the hassle of installing and maintaining BLAST, but of course it is less private and you rely on the network and the current load of these busy servers.

Thanks to its interaction with BLAST, T-Coffee can gather structures and protein profiles and deliver an alignment significantly more accurate than the default you would get with T-Coffee or any similar method. Let us go through the various modes available for T-Coffee


Why do I need BLAST with T-Coffee?
==================================
The most accurate modes of T-Coffee scan the databases for templates that they use to align the sequences. There are currently two types of templates for proteins: 1) structures (PDB) that can be found by a blastp against the PDB database and 2) profiles that can be constructed using either a BLASTP or a PSIBLAST against n.r. or UniProt. These templates are automatically built if you use the following modes:


::

   t_coffee <yourseq> -mode expresso


that fetches and uses structural templates (PDB), or


::

    t_coffee <your seq> -mode psicoffee


that fetches and uses profile templates, or


::

    t_coffee <your seq> -mode accurate


that does everything and tries to use the best template. Now that you see why it is useful, let's see how to get BLAST up and running, from the easy solution to tailor-made ones.


Using the EBI BLAST client
==========================
This is by far the easiest (and the default mode). The perl clients are already incorporated in T-Coffee and all you need are the proper PERL libraries. In theory, T-Coffee should have already installed these libraries during the standard installation. Yet, this requires having root access. It really is worth the effort, since the EBI is providing one of the best webservice available around, and most notably, the only public PSIBLAST via a web service. Note that because PSIBLAST is time consuming, T-Coffee stores the runs in its cache (./tcoffee/cache) so that it does not need to be re-run. It means that if you re-align your sequences (or add a few extra sequences), things will be considerably faster.


Whenever you use a T-Coffee mode requiring BLAST access, it will ask you for an authentification E-mail. Be Careful! If you provide a fake E-mail, the EBI may suspend the service for all machines associated with your IP address (that could mean your entire lab, entire institute, or even the entire country or, but I doubt it, the whole universe). 



Using the NCBI BLAST client
===========================
The NCBI is the next best alternative. In my hand it was always a bit slower and most of all, it does not incorporate PSIBLAST (as a websevice). A big miss! The NCBI web BLAST client is a small executable that you should install on your system following the instructions given on this link:


::

  ftp://ftp.ncbi.nih.gov/blast/executables/LATEST



Simply go for netbl, download the executable that corresponds to your architecture (Cygwin users should go for the win executable). Despite all the files that come along the executable blastcl3 is a stand alone executable that you can safely move to your $BIN.


All you then need to do is to make sure that T-Coffee uses the right client; when you run T-Coffee, specify the client in the command line with:


::

  -blast_server=NCBI


No need for any E-mail here, but you don't get PSIBLAST, and whenever T-Coffee wants to use it, BLASTP will be used instead.


Using another client
====================
You may have your own client (lucky you). If that is so, all you need is to make sure that this client is complient with the BLAST command line. If your client is named foo.pl, all you need to do is run T-Coffee command line with:


::

  -blast_server=CLIENT_foo.pl



Foo will be called as if it were BLASTPGP, and it is your responsability to make sure it can handle the following command line:


::

  foo.pl -p <method> -d <db> -i <infile> -o <outfile> -m 7



- method can either be blastp or psiblast.


- infile is a FASTA file


- -m 7 triggers the XML output. T-Coffee is able to parse both the EBI XML output and the NCBI XML output.


If foo.pl behaves differently, the easiest will probably be to write a wrapper around it so that wrapped_foo.pl behaves like BLASTPGP.


Using a BLAST local version on Unix
===================================
If you have BLASTPGP installed, you can run it instead of the remote clients by using in your command line:


::

  -blast_server=LOCAL



The documentation for BLASTPGP can be found on:


::

  www.ncbi.nlm.nih.gov/staff/tao/URLAPI/blastpgp.html



and the package is part of the standard BLAST distribution:


::

  ftp://ftp.ncbi.nih.gov/blast/executables/LATEST



Depending on your system, your own skills, your requirements and on more parameters than I have fingers to count, installing a BLAST server suited for your needs can range from a 10 minutes job to an achievement spread over several generations. So at this point, you should roam the NCBI website for suitable information.


If you want to have your own BLAST server to run your own databases, you should know that it is possible to control both the database and the program used by BLAST:


::

  -protein_db: will specify the database used by all the PSIBLAST modes of T-Coffee

  -pdb_db: will specify the database used by the structural modesof T-Coffee


.. tip:: T-Coffee is compliant with BLAST+, the latest NCBI Blast.



Using a BLAST local version on Windows/Cygwin
=============================================

BLAST+
------
BLAST+ is the latest NCBI BLAST. It is easier to install. A default installation should be compliant with a default T-Coffee installation.


Original NCBI BLAST
-------------------
For those of you using Cygwin, be careful. While Cygwin behaves like a Unix system, the BLAST executable required for Cygwin (win32) is expecting Windows paths and not Unix paths. This has three important consequences:


1- the NCBI file declaring the sata directory must be:

::

 C:WINDOWS//ncbi.init [at the root of your WINDOWS]



2- the address mentioned with this file must be WINDOWS formated, for instance, on my system:

::

 Data=C:\cygwin\home\notredame\blast\data


3- the database addresses to BLAST must be in Windows format:

::

 -protein_db='c:/somewhere/somewhereelse/database'



(using the slash (/) or the antislash (\) does not matter on new systems but I would recommend against incorporating white spaces.




******************************
T-Coffee Advanced Installation
******************************

These procedures are not needed for default usage of T-Coffee. You will only need to install/configure these packages for specific purposes. T-Coffee is meant to interact with as many packages as possible, either for aligning or using predictions. If you type:


::

   t_coffee



You will receive a list of supported packages that looks like the next table. In theory, most of these packages can be installed by T-Coffee and we welcome any reasonnable request.


::

  ****** Pairwise Sequence Alignment Methods:

  --------------------------------------------

  fast_pair built_in

  exon3_pair built_in

  exon2_pair built_in

  exon_pair built_in

  slow_pair built_in

  proba_pair built_in

  lalign_id_pair built_in

  seq_pair built_in

  externprofile_pair built_in

  hh_pair built_in

  profile_pair built_in

  cdna_fast_pair built_in

  cdna_cfast_pair built_in

  clustalw_pair ftp://www.ebi.ac.uk/pub/clustalw

  mafft_pair http://www.biophys.kyoto-u.ac.jp/~katoh/programs/align/mafft/

  mafftjtt_pair http://www.biophys.kyoto-u.ac.jp/~katoh/programs/align/mafft/

  mafftgins_pair http://www.biophys.kyoto-u.ac.jp/~katoh/programs/align/mafft/

  dialigntx_pair http://dialign-tx.gobics.de/

  dialignt_pair http://dialign-t.gobics.de/

  poa_pair http://www.bioinformatics.ucla.edu/poa/

  probcons_pair http://probcons.stanford.edu/

  muscle_pair http://www.drive5.com/muscle/

  t_coffee_pair http://www.tcoffee.org

  pcma_pair ftp://iole.swmed.edu/pub/PCMA/

  kalign_pair http://msa.cgb.ki.se

  amap_pair http://bio.math.berkeley.edu/amap/

  proda_pair http://bio.math.berkeley.edu/proda/

  prank_pair http://www.ebi.ac.uk/goldman-srv/prank/

  consan_pair http://selab.janelia.org/software/consan/

  ****** Pairwise Structural Alignment Methods:

  --------------------------------------------

  align_pdbpair built_in

  lalign_pdbpair built_in

  extern_pdbpair built_in

  thread_pair built_in

  fugue_pair http://www-cryst.bioc.cam.ac.uk/fugue/download.html

  pdb_pair built_in

  sap_pair http://www-cryst.bioc.cam.ac.uk/fugue/download.html

  mustang_pair http://www.cs.mu.oz.au/~arun/mustang/

  tmalign_pair http://zhang.bioinformatics.ku.edu/TM-align/

  ****** Multiple Sequence Alignment Methods:

  --------------------------------------------

  clustalw_msa ftp://www.ebi.ac.uk/pub/clustalw

  mafft_msa http://www.biophys.kyoto-u.ac.jp/~katoh/programs/align/mafft/

  mafftjtt_msa http://www.biophys.kyoto-u.ac.jp/~katoh/programs/align/mafft/

  mafftgins_msa http://www.biophys.kyoto-u.ac.jp/~katoh/programs/align/mafft/

  dialigntx_msa http://dialign-tx.gobics.de/

  dialignt_msa http://dialign-t.gobics.de/

  poa_msa http://www.bioinformatics.ucla.edu/poa/

  probcons_msa http://probcons.stanford.edu/

  muscle_msa http://www.drive5.com/muscle/

  t_coffee_msa http://www.tcoffee.org

  pcma_msa ftp://iole.swmed.edu/pub/PCMA/

  kalign_msa http://msa.cgb.ki.se

  amap_msa http://bio.math.berkeley.edu/amap/

  proda_msa http://bio.math.berkeley.edu/proda/

  prank_msa http://www.ebi.ac.uk/goldman-srv/prank/

  ####### Prediction Methods available to generate Templates

  -------------------------------------------------------------

  RNAplfold http://www.tbi.univie.ac.at/~ivo/RNA/

  HMMtop www.enzim.hu/hmmtop/

  GOR4 http://mig.jouy.inra.fr/logiciels/gorIV/

  wublast_client http://www.ebi.ac.uk/Tools/webservices/services/wublast

  blastpgp_client http://www.ebi.ac.uk/Tools/webservices/services/blastpgp

  ==========================================================


In our hands all these packages where very straightforward to compile and install on a standard Cygwin or Linux configuration. Just make sure you have gcc, the C compiler, properly installed.

Once the package is compiled and ready to use, make sure that the executable is on your path, so that t_coffee can find it automatically. Our favorite procedure is to create a bin directory in the home. If you do so, make sure this bin is in your path and fill it with all your executables (this is a standard Unix practice).



Installation of M-Coffee
========================
M-Coffee is a special mode of T-Coffee that makes it possible to combine the output of many Multiple Sequence Alignment packages.


Automated installation
----------------------
In the T-Coffee distribution, type:

::

  ./install mcoffee


In theory, this command should download and install every required package. If, however, it fails, you should switch to the manual installation.


Manual installation
-------------------

By default all the packages will be in the following folder:

::

  $HOME/.t_coffee/plugins


If you want to have these packages in a different directory, you can either set the environement variable:

::

  setenv PLUGINS_4_TCOFFEE=<plugins dir>


or use the command line flag -plugin (overrides every other setting):

::

  t_coffee ... -plugins=<plugins dir>


If for some reason, you do not want this directory to be on your path, or you want to specify a precise directory containing the executables, you can use:

::

   export PLUGINS_4_TCOFFEE=<dir>


If you cannot, or do not want to use a single bin directory, you can set the following environment variables to the absolute path values of the executable you want to use. Whenever they are set, these variables will supersede any other declaration. This is a convenient way to experiment with multiple package versions:

::

  POA_4_TCOFFEE CLUSTALW_4_TCOFFEE TCOFFEE_4_TCOFFEE MAFFT_4_TCOFFEE MUSCLE_4_TCOFFEE
  DIALIGNT_4_TCOFFEE PRANK_4_TCOFFEE DIALIGNTX_4_TCOFFEE 


For three of these packages, you will need to copy some of the files in a special T-Coffee directory:

::

   cp POA_DIR/* ~/.t_coffee/mcoffee/

   cp DIALIGN-T/conf/* ~/.t_coffee/mcoffee

   cp DIALIGN-TX/conf/* ~/.t_coffee/mcoffee


If you would rather have the mcoffee directory in some other location, set the MCOFFEE_4_TCOFFEE environement variable to the propoer directory:

::

   setenv MCOFFEE_4_TCOFFEE <directory containing mcoffee files>
   

Note that the following files are enough for default usage:

::

  BLOSUM.diag_prob_t10 BLOSUM75.scr blosum80_trunc.mat

  dna_diag_prob_100_exp_330000 dna_diag_prob_200_exp_110000

  BLOSUM.scr BLOSUM90.scr dna_diag_prob_100_exp_110000

  dna_diag_prob_100_exp_550000 dna_diag_prob_250_exp_110000

  BLOSUM75.diag_prob_t2 blosum80.mat dna_diag_prob_100_exp_220000

  dna_diag_prob_150_exp_110000 dna_matrix.scr


Configuration for PDB (installed locally)
=========================================
For all the structural modes of T-Coffee (Expresso, 3D-Coffee, tRMSD, iRMSD, etc...), access to structural information is mandatory. You can do so either by having a database installed locally on your own system or by accessing the PDB through the webserver.
If you do not have PDB installed, don't worry, T_Coffee will go and fetch any structure it needs directly from the PDB repository. It will simply be a bit slower than if you had PDB locally. 
If you prefer to have access to a local installation of the PDB in your file system, you have to indicate to T-Coffee their location in your system using the following commands:

::

  setenv (or export) PDB_DIR <abs path>/data/structures/all/pdb/

  OR

  setenv (or export) PDB_DIR <abs path>/structures/divided/pdb/



Installation of tRMSD
=====================
tRMSD comes along with t_coffee but it also requires the package phylip in order to be functional. Phylip can be obtained from:


::

  Package Function

  ===================================================

  ---------------------------------------------------

  Phylip phylogenetic tree computation

  evolution.genetics.washington.edu/phylip.html

  ---------------------------------------------------

  t_coffee -other_pg trmsd


Installation of 3D-Coffee/Expresso
==================================
3D-Coffee/Expresso is a special mode of T-Coffee that makes it possible to combine sequences and structures. The main difference between Expresso and 3D-Coffee is that Expresso fetches the structures itself.


Automated Installation
----------------------
In the T-Coffee distribution, type:


::

  ./install expresso

  OR

  ./install 3dcoffee



In theory, this command should download and install every required package (except fugue). If, however, it fails, you should switch to the manual installation (see next).


Manual Installation
-------------------
In order to make the most out of T-Coffee, you will need to install the following packages (make sure the executable is named as indicated below):


::

  Package Function

  =============================================================

  -------------------------------------------------------------

  wget 3DCoffee
  Automatic downloading of structures
   
  -------------------------------------------------------------
  
  sap structure/structure comparisons
  Obtained from W. Taylor, NIMR-MRC
  
  -------------------------------------------------------------
 
  TMalign zhang.bioinformatics.ku.edu/TM-align/
  
  -------------------------------------------------------------
  
  mustang www.cs.mu.oz.au/~arun/mustang/
  
  -------------------------------------------------------------
  
  wublastclient www.ebi.ac.uk/Tools/webservices/clients/wublast
  
  -------------------------------------------------------------
  
  Blast www.ncbi.nih.nlm.gov
  
  -------------------------------------------------------------

  Fugue protein to structure alignment program
  http://www-cryst.bioc.cam.ac.uk/fugue/download.html

   ***NOT COMPULSORY***
   
  -------------------------------------------------------------


Once the package is installed, make sure make sure that the executable is on your path, so that T-Coffee can find it automatically.


The wublast client makes it possible to run BLAST at the EBI without having to install any database locally. It is an ideal solution if you are only using Expresso occasionally.


Installing Fugue for T-Coffee
-----------------------------
Uses a standard Fugue installation. You only need to install the following packages: joy, melody, fugueali, sstruc, hbond. If you have root privileges, you can install the common data in:

::


 cp fugue/classdef.dat /data/fugue/SUBST/classdef.dat


otherwise:

::


 Setenv MELODY_CLASSDEF=<location>

 Setenv MELODY_SUBST=fugue/allmat.dat


All the other configuration files must be in the right location.


Installation of R-Coffee
========================
R-Coffee is a special mode able to align RNA sequences while taking into account their secondary structure.


Automated installation
----------------------
In the T-Coffee distribution, type:


::

  ./install rcoffee


In theory, this command should download and install every required package (except Consan). If, however, it fails, you should switch to the manual installation (see next).


Manual installation
-------------------
R-Coffee only requires the package Vienna to be installed, in order to compute Multiple Sequence Alignments. To make the best out of it, you should also have all the packages required by M-Coffee.


::

  Package Function

  ===================================================

  ---------------------------------------------------
  
  Consan computes highly accurate pairwise alignments
  selab.janelia.org/software/consan/
  
  ***NOT COMPULSORY***
    
  ---------------------------------------------------
  
  RNAplfold computes RNA secondary structures
  www.tbi.univie.ac.at/~ivo/RNA/
  
  ---------------------------------------------------
  
  ProbConsRNA probcons.stanford.edu/
  
  ---------------------------------------------------

  M-Coffee T-Coffee and the most common MSA Packages
  (cf M-Coffee in this installation guide)

  ---------------------------------------------------
  

Installing ProbConsRNA for R-Coffee
-----------------------------------
Follow the installation procedure, but make sure you rename the probcons executable into probconsRNA.


Installing Consan for R-Coffee
------------------------------
In order to insure a proper interface bewteen Consan and R-Coffee, you must make sure that the file mix80.mod is in the directory ~/.t_coffee/mcoffee or in the mcoffee directory otherwise declared.


