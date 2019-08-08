T-Coffee
=========

T-Coffee is a collection of tools for Computing, Evaluating and Manipulating 
Multiple Alignments of DNA, RNA, Protein Sequences and Structures.


Prerequisites
--------------
T-Coffee compilation requires the following tools installed on your system ``make``, ``gcc-c++``, ``g77``, ``Perl``, ``sphinx``, ``CPAN``. 


Compile 
--------

Note that the git repository is meant to support the packaging of new releases. 
If you simply want to install T-Coffee on your system please download the latest distribution from
	http://www.tcoffee.org/Projects/tcoffee/index.html#DOWNLOAD

All required installation documentation is available from
	http://www.tcoffee.org/Projects/tcoffee/documentation/index.html#t-coffee-installation

If you nontheless want to compile from git: 

	Clone the git repository on your computer with the following command: 
		git clone git@github.com:cbcrg/tcoffee.git tcoffee
	make sure you have installed all dependencies
	cd tcoffee/t_coffee/src
	make t_coffee
	mv t_coffee <bin location> OR add current dir to your path - not recommended

