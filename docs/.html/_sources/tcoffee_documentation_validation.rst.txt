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

The simplest way to check the documentation is to run the following command:

 ::

  ##: ./lib/perl/lib/perl4makefile/doc2test.pl -mode update

This script will extract all tcommand lines from the .rst in docs and will run these command lines. When a reference output is available in /testsuite/docs/ref/ it will produce an other output and compare it with the reference. The new output will be put in /testsuite/docs/latest/. The script will eventually produce a global report that is to be found in /testsuite/docs/log/validation.log. The output of the failed scripts will be put in /testsuite/docs/failed. All the reference output will be added to git repository.

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

::

  ##: ./lib/perl/lib/perl4makefile/doc2test.pl -clean

Will cause all the files in ./examples/ that are not used for any /docs/*.rst commands to be deleted. They will be deleted from both their current locations and the git repository. The remaining files will be added to the git repository.