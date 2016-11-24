######################################
T-Coffee doc test
######################################

************
Introduction
************
Before the realease of a new T-Coffee version, all command lines in the documentation are programmatically tested. The purpose of this document is to explain how new command lines can be added to the documentation and which procedure should be carried out in order to perform a test before a new release 

*********************************
Introducing a new command to test
*********************************

Tests are systematically carried out on all the command lines that start with the symbol $$

::

  $$: t_coffee

 
Other command lines starting with different symbols are not systematically checked. Two types of command lines identifiers are used:

::

  $#: t_coffee

For a command that could be tested but will not be, either for the sake of time or because it is currently unstable. When a new release needs to be urgently made available because of a critical fix, it is advisable to comment out this way any command line failing the test.

::

  ##: t_coffee

These commands are never meant to be tested, either because they contain system dependant information or non programmatic information.

Whenvever adding a new command, it is usually advisable to provide the corresponding files. These files can be found in the example directory in the top of the github repository. New commands can be built using the existing files, or they can depend on files newly added to the repository.

**************************
Checking the documentation
**************************

The simplest way to check the documentation is to run the following command:

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

  ##: ./lib/perl/lib/perl4makefile/doc2test.pl -mode update -stop_on_failed

It will cause the test to stop whenever a failed is encountered

::

  ##: ./lib/perl/lib/perl4makefile/doc2test.pl -mode reset

Will cause all the reference files to be erased.

