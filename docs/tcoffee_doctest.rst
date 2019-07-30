######################################
T-Coffee doc test
######################################

************
Introduction
************
Before the realease of a new T-Coffee version, many command lines in the documentation are programmatically tested. The purpose of this document is to explain which command lines are tested, how new command lines can be added to the documentation and how to perform a documentation check before a new release. 

*********************************
Introducing a new command to test
*********************************

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

**************************
Checking the documentation
**************************

The simplest way to check the documentation is to run the following command from the repository root (other locations are possible, but -ref, -docs must be modified accordingly):

 ::

  ##: ./lib/perl/lib/perl4makefile/doc2test.pl -mode update

This script will extract all tcommand lines from the .rst in docs and will run these command lines. When a reference output is available in /testsuite/docs/ref/ it will produce an other output and compare it with the reference. The new output will be put in /testsuite/docs/latest/. The script will eventually produce a global report that is to be found in /testsuite/docs/log/validation.log. The output of the failed scripts will be put in /testsuite/docs/failed

::

  ##: ./lib/perl/lib/perl4makefile/doc2test.pl -mode new

This mode will only run the new command lines and report (defined as those not having a reference output in /testsuite/docs/ref/). This mode is the default mode

::
  ##: ./lib/perl/lib/perl4makefile/doc2test.pl -mode failed

This mode will only run the command lines that have previously failed, as indicated by the presence of an output in the failed directory.


***************************
Compiling the documentation
***************************

The documentation can be compiled using the following command

::

  ##: ./lib/perl/lib/perl4makefile/doc2test.pl -mode update -stop

It will cause the test to stop whenever a failed is encountered

::

  ##: ./lib/perl/lib/perl4makefile/doc2test.pl -mode reset

Will cause all the reference files to be erased.

*************************
Validating a distribution
*************************

Dump files are safe contained. It is possible to check T-Coffee capacity to reproduce a collection of reference dump files.

::

  ##: ./lib/perl/lib/perl4makefile/doc2test.pl -mode check -ref <dir containing dump files>

It will cause the test to stop whenever a failed is encountered

::

  ##: ./lib/perl/lib/perl4makefile/doc2test.pl -mode reset

Will cause all the reference files to be erased.

*************************
Creating a new command   
*************************

In order to generate a new test command, all you need to do is to run your command while setting an environement variable. This will generate a dump file that can be used to rerun the same call:
Dump files are very easy to produce yout simply need to run

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


