#####################
T-Coffee Installation
#####################
.. warning:: This chapter has been extensively updated in 11/2016. Installation of T-Coffee versions lower than version 9.03 may now be deprecated. Contact us if you need to install an older version !

************
Installation
************
This chapter describes the installation procedures relevant for a standard use of T-Coffee on the most common operative systems: Unix/Linux and Mac OS. T-Coffee can only be installed on windows using virtual box simulating a Linux environement. We maintain two versions of T-Coffee:


 1) the ** `stable version <http://tcoffee.org/Packages/tcoffee/Stable/> `** This is an extensively evaluated version, either linked to the latest paper or benefiting from an important fix- ** recommended **.
 2) the ** `beta version <http://www.tcoffee.org/Packages/Beta/Latest/> `** is the latest release. It passes all the main tests and is usualy of higher quality that the stable. It is recommended to start installing this version if you experience any issue with the stable.
 3) the ** `archives  <http://www.tcoffee.org/Packages/Archives/> `** all previous releases of T-Coffee. If the one you need is missing please contact us. 

T-Coffee is a complex package that interacts with many other third party software and/or servers (such as BLAST, see next section). All available versions on the server (starting from version 1.08) will install on your computer all the third party packages and setup the variables required for the different T-Coffee options to run correctly. Whenever the automated installation fails because of unforeseen system specificities, **don't hesitate to contact us**. If the installation was successfull with the exception of some packages, users should install the third party package manually. This documentation gives some useful tips, but users are encouraged to send their feedbacks and share their experiences in order to improve it.

Unix/Linux Binaries
===================

::

  1. Download the installer package corresponding to your system from:
     http://tcoffee.org/Packages/Stable/Latest/

  2. Grant execution permission to the downloaded file with the following command:
     ##: chmod +x T-COFFEE_installer_"version_x".bin

  3. Launch the installation wizard with:
     ##: ./T-COFFEE_installer_"version_x".bin

  4. Follow the wizard instructions and complete the installation
  
  5. Open a new terminal session to be sure that your environment is updated
  
  6. Type the following command to verify the installation was successful:
     $$: t_coffee -version


MacOS Binaries - Deprecated
===========================

These binaries are not any more supported as such. However, the source installer - see below - will automatically retrive pre-compiled MacOS binaries if compilation fails on your system. 


Installation From Source/Binaries downloader (Mac OSX/Linux)
============================================================ 

The following procedure shows howto install T-Coffee from the pre-packages source distribution. This procedure is recommanded because the version of T-Coffee available there is guarranteed to have passed a set of pre-specified tests. By contrast, the github version is untested, and therefore ongoing work... In most cases the binaries will be downloaded from a repository but you may also trigger local compilation. 

::

  1. Download the desired installer:
     http://www.tcoffee.org/Packages/Stable/Latest/T-COFFEE_distribution_Version_<version>.tar.gz
     http://www.tcoffee.org/Packages/Beta/Latest/T-COFFEE_distribution_Version_<version>.tar.gz
     http://www.tcoffee.org/Packages/Archives/T-COFFEE_distribution_Version_<version>.tar.gz
     https://github.com/cbcrg/tcoffee

  2a. untar the file
     ##: tar -xvf T-COFFEE_distribution.tar.gz
     ##: cd  T-COFFEE_distribution
  
  2b. clone or unzip the github repo
     ##: cd t_coffee/src


  3a. Launch the installer
     ##:  ./install tcoffee
     
  
  3b. If you want to use remote services (psiblast, expresso). You must enter a proxy value and a valid e-mail. These values will be stored in <your home>/.t_coffee/.t_coffee_env and can be edited at any time. The proxy is not compulsory. It is usualy needed when working behind a firewall.
     ##: ./install tcoffee -email=<your email> -proxy=<your proxy>
    

  4. as you can also install any of the T-Coffee mode you are interested in. Get the list with:
     ##: ./install 
     Then install one, for instance
     ##: ./install expresso
     Or install them all
     ##: ./install all
    
  
  5. By default, the installer will not re-install any component already available on your system. If you want to update you must specify
     ##: ./install <any component> -update
     OR
     ##: ./install <any component> -update -nobinaries
     For instance the following will cause T-Coffee to be recompiled on your system
     ##: ./install t_coffee -update -nobinaries
  
  6. When you are done you will need to make this installation permanent by adding the following lines *at the bottom* of your configuration file (typically .bashrc)
     ##: export PATH=<your home>/.t_coffee/bin/<macosx|linux>:$PATH



******************
BLAST and T-Coffee
******************
BLAST is a program that searches databases for homologues of a query sequence. It works for protein and nucleic acid sequences alike. In theory BLAST is just a package like any but in practice things are a bit more complex. To run correctly, BLAST requires up-to-date databases (that can be fairly large, like nr or UniProt) and a powerful computer. Fortunately, an increasing number of institutes or companies are now providing BLAST clients that run over the net. It means that all you need is a small program that send your query to the big server and gets the results back. This prevents you from the hassle of installing and maintaining BLAST, but of course it is less private and you rely on the network and the current load of these busy servers.

**Thanks to its interaction with BLAST, T-Coffee can gather more information and deliver alignments significantly more accurate than the default T-Coffee or any similar method. Let us go through the various modes available for T-Coffee...**


Why do I need BLAST with T-Coffee?
==================================
The most accurate modes of T-Coffee scan the databases for templates that they use to align the sequences. Let's see how to get BLAST up and running, from the easy solution to tailored ones. There are currently two types of templates for proteins: 

 1) **structures**, that can be found by a BLASTP against the PDB database.
 2) **profiles**, constructed using BLASTP or PSI-BLAST against nr or UniProt. 
 
These templates are automatically built by T-Coffee when using one of the following modes:

::

   To fetch and use structural templates:
   ##: t_coffee <yourseq> -mode expresso

   To fetch and use profile templates:
   ##: t_coffee <your seq> -mode psicoffee
   
   To fetch everything possible and get the best templates, structure or profile:
   ##: t_coffee <your seq> -mode accurate
   
   
Using the EBI BLAST client
==========================
This is by far the easiest way and conveniently the default mode of T-Coffee. The PERL clients are already incorporated in T-Coffee and all you need are the proper PERL libraries. In principle, T-Coffee should have already installed these libraries during the standard installation, yet, this requires having root access. It really is worth the effort since the EBI is providing one of the best webservice available around and most notably, the only public PSI-BLAST via a webservice. Note that because PSI-BLAST is time consuming, T-Coffee stores the runs in its cache (**./tcoffee/cache**) so that it does not need to be rerun. It means that if you realign your sequences (or add a few extra sequences), things will be considerably faster.

.. danger:: Whenever you use a T-Coffee mode requiring BLAST access, it will ask you for an authentification e-mail. Be extra careful!!! If you provide a fake e-mail, the EBI may suspend the service for all machines associated with your IP address (that could mean your entire lab, entire institute, even the entire country or, but I doubt it, the whole universe). 

.. tip:: Files in the cache are never erased so remember to empty the cache from time to time otherwise it's just getting bigger and bigger...


Using the NCBI BLAST client
===========================
The NCBI is the next best alternative however in my hands it was always a bit slower and, most of all, it does not incorporate PSI-BLAST as a webservice. A big miss! The NCBI web BLAST client is a small executable that you should install on your system. To do so, you just have to follow the instructions given on this `link <ftp://ftp.ncbi.nih.gov/blast/executables/LATEST>`_. Simply go for netbl, download the executable that corresponds to your architecture (Cygwin users should go for the win executable). Despite all the files that come along the executable blastcl3 is a stand alone executable that you can safely move to your $BIN. All you then need to do is to make sure that T-Coffee uses the right client; when you run T-Coffee, specify the client in the command line with the flag **-blast_server=NCBI**.

.. Attention:: No need for any e-mail here, but you don't get PSI-BLAST. Whenever T-Coffee will need to use it, BLASTP will be used instead.


Using another client
====================
You may have your own client (lucky you). If that is so, all you need is to make sure that this client is complient with the BLAST command line. If your client is named foo.pl, all you need to do is run T-Coffee command line with the flag **-blast_server=CLIENT_foo.pl**. Foo will be called as if it were BLASTPGP, and it is your responsability to make sure it can handle the following command line.

::

  ##: foo.pl -p <method> -d <db> -i <infile> -o <outfile> -m 7

  "method"  : BLAST method for the search ("blastp" or "psiblast")
  "db"      : database used for the search
  "infile"  : input sequence(s) in FASTA format
  "outfile" : name the output file 
  "-m 7"    : triggers the XML output (parses both the EBI & NCBI XML output)

.. tip:: If foo.pl behaves differently, the easiest way will probably be to write a wrapper around it so that wrapped_foo.pl behaves like BLASTPGP.


Using a BLAST local version on Unix
===================================
If you have BLASTPGP installed, you can run it instead of the remote clients by using in your command line the flag **-blast_server=LOCAL**. The documentation for BLASTPGP can be found `here <http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/blastpgp.html>`_ and the package is part of the standard BLAST `distribution <ftp://ftp.ncbi.nih.gov/blast/executables/LATEST>`_. Depending on your system, your own skills, your requirements and on more parameters than I have fingers to count, installing a BLAST server suited for your needs can range from a 10 minutes job to an achievement spread over several generations. So at this point, you should roam the NCBI website for suitable information. If you want to have your own BLAST server to run your own databases, you should know that it is possible to control both the database and the program used by BLAST using T-Coffee flags  **-protein_db** (will specify the database used by all the PSI-BLAST modes) and **-pdb_db** (will specify the database used by the structural modes)

.. tip:: T-Coffee is compliant with BLAST+, the latest NCBI BLAST.


Using a BLAST local version on Windows/Cygwin
=============================================
BLAST+ is the latest NCBI BLAST. It is easier to install and a default installation should be compliant with a default T-Coffee installation. For those of you using Cygwin, be careful!! While Cygwin behaves like a Unix system, the BLAST executable required for Cygwin (win32) is expecting Windows paths and not Unix paths. This has three important consequences:

::

  1. The NCBI file declaring the sata directory must be:
     C:WINDOWS//ncbi.init [at the root of your WINDOWS]

  2. The address mentioned with this file must be WINDOWS formated, for example:
     Data=C:\cygwin\home\notredame\blast\data

  3. The database addresses to BLAST must be in Windows format:
     ##: t_coffee ... -protein_db='c:/somewhere/somewhere else/database'

.. attention:: Using the slash (/) or the antislash (\\) does not matter on new systems but I would recommend against incorporating white spaces.


***************
Troubleshooting
***************

Third party packages
====================
These procedures are not needed for default usage of T-Coffee. You will only need to install/configure these packages for specific purposes. T-Coffee is meant to interact with as many packages as possible, especially for aligning or using predictions. You will receive a list of supported packages that looks like the next table if you simply type **t_coffee**:

::

  Command:
  $$: t_coffee

  Display the list of supported packages:
 
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
  fugue_pair http://mizuguchilab.org/fugue/
  pdb_pair built_in
  sap_pair https://mathbio.crick.ac.uk/wiki/Software#SAP
  mustang_pair http://lcb.infotech.monash.edu.au/mustang/
  tmalign_pair https://zhanglab.ccmb.med.umich.edu/TM-align/

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
  HMMtop http://www.enzim.hu/hmmtop/
  GOR4 http://mig.jouy.inra.fr/logiciels/gorIV/
  wublast_client http://www.ebi.ac.uk/Tools/webservices/services/wublast
  blastpgp_client http://www.ebi.ac.uk/Tools/webservices/services/blastpgp

.. tip:: In our hands all these packages where very straightforward to compile and install on a standard Cygwin or Linux configuration. Just make sure you have gcc, the C compiler, properly installed. Once the package is compiled and ready to use, make sure that the executable is on your path, so that T-Coffee can find it automatically. Our favorite procedure is to create a bin directory in the home. If you do so, make sure this bin is in your path and fill it with all your executables (this is a standard Unix practice).


M-Coffee parameters
===================
M-Coffee is a special mode of T-Coffee that makes it possible to combine the output of many Multiple Sequence Alignment packages. By default all the packages will be in the following folder **$HOME/.t_coffee/plugins/linux/**. If you want to have these packages in a different directory, you can either set the environment variable (option 1) or use the flag **-plugin** (to override every other setting). If for some reason, you do not want this directory to be on your path or you want to specify a precise directory containing the executables, you can use option 2. You can also set the following environment variables to the absolute path of the executable you want to use option 3: whenever they are set these variables will supersede any other declaration. This is a convenient way to experiment with multiple package versions. If you would rather have the mcoffee directory in some other location, set the MCOFFEE_4_TCOFFEE environement variable to the proper directory (option 4).

::

  Option 1: set the environment variable
  ##: setenv PLUGINS_4_TCOFFEE=<plugins dir>
  
  Option 2: specify the directory
  ##: export PLUGINS_4_TCOFFEE=<dir>
  
  Option 3:
  ##: POA_4_TCOFFEE CLUSTALW_4_TCOFFEE TCOFFEE_4_TCOFFEE MAFFT_4_TCOFFEE \
  MUSCLE_4_TCOFFEE DIALIGNT_4_TCOFFEE PRANK_4_TCOFFEE DIALIGNTX_4_TCOFFEE
  
  Option 4:
  ##: setenv MCOFFEE_4_TCOFFEE <directory containing mcoffee files>
  
 
To be able to run M-Coffee, these following files are enough for a default usage:

::

  BLOSUM.diag_prob_t10 BLOSUM75.scr blosum80_trunc.mat
  dna_diag_prob_100_exp_330000 dna_diag_prob_200_exp_110000
  BLOSUM.scr BLOSUM90.scr dna_diag_prob_100_exp_110000
  dna_diag_prob_100_exp_550000 dna_diag_prob_250_exp_110000
  BLOSUM75.diag_prob_t2 blosum80.mat dna_diag_prob_100_exp_220000
  dna_diag_prob_150_exp_110000 dna_matrix.scr


Structural modes (using PDB)
============================
Expresso/3D-Coffee are special modes of T-Coffee that allow to combine sequences and structures to reach more accurate alignments. T-Coffee proposes also other tools (iRMSD/APDB, T-RMSD, etc...) requiring access to structural information. You can do so either by having a database installed locally on your own system or by accessing the PDB through the web server. If you do not have PDB installed, don't worry, T-Coffee will go and fetch any structure it needs directly from the PDB repository, it will simply be a bit slower. If you prefer to have access to a local installation of the PDB in your file system, you have to indicate their location in your system using one of the following commands:

::

  Using a local version of the PDB database:
  ##: setenv (or export) PDB_DIR <PATH>/data/structures/all/pdb/
  ##: setenv (or export) PDB_DIR <PATH>/structures/divided/pdb/

The T-RMSD tools comes along with T_Coffee package in order to build clustering based on structure. In addition to structural information it also requires the package Phylip, containing lots of phylogenetic tree reconstruction tools. If you need more information about the different Phylip tools, information can be obtained `here <http://www.evolution.genetics.washington.edu/phylip.html>`_. 

R-Coffee associated packages
============================
R-Coffee is a special mode able to align RNA sequences while taking into account their secondary structure. R-Coffee only requires the package Vienna to be installed, in order to compute Multiple Sequence Alignments. To make the best out of it, you should also have all the packages required by M-Coffee.

 - `Consan <http://eddylab.org/software/consan/>`_ from Eddy/Riva laboratory.    
 - `RNAplfold <http://www.tbi.univie.ac.at/RNA/>`_ from the Vienna package.
 - `ProbConsRNA <http://probcons.stanford.edu/download.html>`_ from Stanford university.
 
 
.. tip:: Regarding ProbConsRNA, make sure you rename the probcons executable into ProbConsRNA.

.. tip:: In order to insure a proper interface bewteen Consan and R-Coffee, make sure that the file mix80.mod is in the directory **~/.t_coffee/mcoffee** or in the mcoffee directory otherwise declared.

