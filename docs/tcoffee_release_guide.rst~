######################################
T-Coffee release guide
######################################

************
Introduction
************
Before the realease of a new T-Coffee version, many command lines in the documentation are programmatically tested. The purpose of this document is to explain how to make a new release, which command lines are tested, how new command lines can be added to the documentation and how to perform a documentation check before a new release. 


********************
Making a new release
********************

T-Coffee is automatically build by using the [Circle CI](https://circleci.com/gh/cbcrg/tcoffee) service. The procedure starts with an update of the version number followed by a git push that triggers the validation of the new release on circleCI and a posting via an amzon instance hosting the T-Coffee web site. All these operations are detailed below. Depending on the release type the correspinding counter will be incremented by 1 <nmajor>.<stable>.<beta>
	
	
::

  $#: cd t_coffee/src
  $#: make release type=beta(default)|stable|major m='short description of the main feature' test=core(default)|all|full|remote|docs|web|none
  
core  : will test a throrough set of commands
all   : will run all the test commands packaged as dumps in tests
full  : will test all the commands packaged as dumps in the git repo - they are normally in docs and tests
remote: will run all tze tests command dealing with remote access (BLAST and PDB)
docs  : will run all the dumps compiled from the rst files 

Beta Release
============

For a beta, the following needs to go through:

::
  $#: make release type=beta test=core|none

Stable Release
==============

For a Stable, the following needs to go through:

::
  $#: make release type=beta test=core
  $#: make release type=stable teste=none

Major Release
==============

For a Major, the following needs to go through. Note that a major implies some important milestone such as a publication. For a start, you should make sure the documentation compiles well. If needed you should recompile all the available tests


::
  $#: doc2test.pl -play docs/doc2test.rst -data docs/.data/ -dump docs/.dumps/ -update
  $#: doc2test.pl -play tests/doc2test.rst -data tests/.data/ -dump tests/.dumps/ -update


Then you should make sure ALL the tests pass well on your local machine. This is a rsather lengthy procedure that you want to bulletpruf localy before running it on CircleCI

::
  $#: doc2test.pl -replay  . 

If they do, you should the run the following

::
  $#: make release type=beta test=full


If this beta is succesful, then you can pack the major

  $#: make release type=major teste=none
 	

**************************
Release Building Procedure
**************************

This section explains how the release procedure works. Roughly speaking, the release is initiated on a local machine using the makefile (t_coffee/src/makefile). It is this procedure that increments the release version number (lib/version/version_number.version). This procedure makes a dummy distribution on the local machine and eventually pushes all the files in github. Thanks to a CircleCI hook, the push then initiates the building of a virtual machine in CircleCI where the final building procedure is carried out amd followed by a test (as specifiled in makefile test=<test type>. If the test is succesful, the procedure finishes with a publishing of the distribution. 




Triggering a release
====================

- go into t_coffee/src
- type make release=beta|stable|major m='commit comments -- including [ci skip] if needed

This will trigger 
	-an  update of the version number in lib/version/version_number.version: Version_<major>.<stable>.<build>.<github Branch Number>
	-the file .circleci/config.yml will be programmatically edidted by read_program_version.pl to state the kind of release to be done by CircleCI
 	-make distribution will run. Among othe things it will
		-check that the distribution procedure is functional - it is a similar procedure that will be used by CircleCI
		-compile the documentation with sphinx and update docs/.html -- this requires sphinx to be installed on the local machine


CircleCI building
=================

The tcoffee git repoitory is registered to CircleCI via the following hook:
	https://github.com/cbcrg/tcoffee/settings/hooks		

The Circle build is controlled by the [.circleci/circle.yml](circle.yml) file and processes as follows:

1. The T-Coffee [build/build.sh](build/build.sh) script is executed in the [cbcrg/tcoffee-build-box:1.2](docker/Dockerfile.buildbox) container

2. The T-Coffee binaries produced by the build process are created in the folder 
  `~/publish/sandbox/build`. This folder is moved under path `build/tcoffee`
  
3. A new Docker image named `xcoffee` is created copying the content of the folder `build/tcoffee`. 
  The image is built by using this [docker/Dockerfile](docker/Dockerfile). 
  
4. The new `xcoffee` build is tested by running the [docker/run-tests.sh](docker/run-tests.sh) 
  tests suite. 

5. The summary of the tests is available on https://circleci.com/gh/cbcrg/tcoffee/tree/master
  
5. If all tests are passed the `xcoffee` image is pushed to the [Docker Hub](https://hub.docker.com/r/cbcrg/tcoffee/tags/) 
  with the names `cbcrg/tcoffee:latest` and `cbcrg/tcoffee:<version.commit-id>`

7. The environment variable `RELEASE=0|1` is used to mark the build as beta or stable 
  (use the [build/make_release.sh] to trigger a new release build).

Publication
===========
Once the build is complete and all tests are passed the distribution is pushed onto the web along with the associated documentation.



****************************
Compiling and running tests
***************************

The ititility doc2test.pl makes it possible to run a command line on data located in a folder and to produce a dump file that will then be a stand alone tests, that can be transfered and replayed. The dump file is a container in which the input, the stdoin, stdout, environement and output of the run are recorded. 


Tests are systematically carried out on all the command lines that start with the symbol $$. For instance, the following CL will has been tested.

::

  $$: t_coffee

 
Other command lines starting with different symbols are not checked. Two types of command lines identifiers are used:

::

  $#: t_coffee

For a command that could be tested but will not be, either for the sake of time or because it is currently unstable. When a new release needs to be urgently made available because of a critical fix, it is advisable to comment out this way non critical command lines failing the test.

::

  ##: t_coffee

These commands are never tested, either because they contain system dependant information or non programmatic information.

Whenvever adding a new command, input files must be added to the repository directory ./examples/. New commands can be also be built using the existing files, or they can depend on files newly added to the repository.


Using doc2test.pl
=================

The utility doc2test is meant to do the following actions:

The play action will compile dump files. 

 ::

  #$: ./lib/perl/lib/perl4makefile/doc2test.pl -play <*.rst file or *.tests file or list of cmd in txt> -data <folder containing the data> -dumps <target directory for dumps>


The Check action will recursively go through any directory and report the exit status of each dump
::

  #$: ./lib/perl/lib/perl4makefile/doc2test.pl -check <dump directory - recursive>

Note that while compiling cls, you can use check to clean your directory
 
::

  #$: ./lib/perl/lib/perl4makefile/doc2test.pl -check <dump directory - recursive> -clean FAILED or -clean2 <string in dump to be deleted>

You should also use check to make a fresh build and remove all previous dumps:
 
::
  #$: ./lib/perl/lib/perl4makefile/doc2test.pl -check <dump directory - recursive> -clean ALL

And you can the re-run the play that will only regenerate missing dumps:
 

::
#$: ./lib/perl/lib/perl4makefile/doc2test.pl -play <*.rst file or *.tests file or list of cmd in txt> -data <folder containing the data> -dumps <target directory for dumps>

The replay action will re-run each dump and report the status. It will report new warnings.


::
  #$: ./lib/perl/lib/perl4makefile/doc2test.pl -replay <dump directory - recursive>


The unplay option maked it possible to output the inoput files of every dump

::
  #$: ./lib/perl/lib/perl4makefile/doc2test.pl -replay <dump directory - recursive> -outdir <file in which data should be dumped>



Compiling tests from the rst documentation
==========================================

The rst documentation contain lines tagged with two $$ signs. These lines will be selected for validation and used to build dump files. In order to run a fresh collection of dumps for the full documentation, run the foillowing commands:

::

The simplest way to compile all possible tests from the the documentation is to run the following command from the top of the git repo.

 ::
  #$: ./lib/perl/lib/perl4makefile/doc2test.pl -check docs/.dumps -clean ALL
  #$: ./lib/perl/lib/perl4makefile/doc2test.pl -play docs/doc2test.rst -data docs/.data -dump docs/.dumps

The output can then be checked

::

  #$: ./lib/perl/lib/perl4makefile/doc2test.pl -check docs/.dumps

The failing jobs can be removed and re-run once the issue has been identified (usualy a missing file):

::

  #$: ./lib/perl/lib/perl4makefile/doc2test.pl -check docs/.dumps -clean FAILED

Note that in order to figure out an issue you can use the replay debug option

  #$: ./lib/perl/lib/perl4makefile/doc2test.pl -replay <dump file> -debug


You can the re-run all the dumps. This is how the tests will be carried out by CircleCI

::
  #$: ./lib/perl/lib/perl4makefile/doc2test.pl -replay docs/.dumps

Compiling tests commands
========================
tests are special types of pre-compiled commands used to validate specific components of T-Coffee. You can compile and use thsis test in the same way as rst files. 


Creating a new command to be tested  
===================================

In order to generate a new test command, all you need to do is to run your command while setting an environement variable. This will generate a dump file that can be used torerun the same call. Dump files are very easy to produce yout simply need to run

::

  ##:rm -f <your dump file>;export DUMP_4_TCOFFEE=<your dump file>;t_coffee -in seq 
  OR
  ##:rm -f <your dump file>;export DUMP_4_TCOFFEE=<your dump file>;t_coffee -other_pg seq_reformat -in xxxx/xx/s.pep -output fasta_seq > yyy

This command will generate a self contained dumpfile. This dumpfile will contain all the information needed ton reprodduce the runn (i.e file name, path, content and T-Coffee parameters). In order to reproduce the call you need to have T-Coffee installed and run:

::

  ##: ./lib/perl/lib/perl4makefile/doc2test.pl -replay yourdumpfile

Will allow you to check if the command runs and produces similar files. By default the test is "PASSED" if no error is thrown by the call. This can be refined as follows:

  ##: ./lib/perl/lib/perl4makefile/doc2test.pl -replay yourdumpfile -strict

Will report failure whenever an output file is missing or whenever an error is reported. 
  
  ##: ./lib/perl/lib/perl4makefile/doc2test.pl -replay yourdumpfile -very_strict

Will report failure whenever there is a warning OR whenever an output file differs (Note that VERSION and CPU are excluded from the comparison as a consequence different versions will not result in different output files). It is possible to visualize the replayed dump for debugging purposes.
  
  ##: ./lib/perl/lib/perl4makefile/doc2test.pl -replay yourdumpfile -keepreplayed

