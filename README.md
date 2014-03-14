T-Coffee
=========

T-Coffee is a collection of tools for Computing, Evaluating and Manipulating 
Multiple Alignments of DNA, RNA, Protein Sequences and Structures.


Prerequisites
--------------
T-Coffee compilation requires the following tools installed on your system ``make``, ``gcc-c++``, ``g77``, ``Perl`` and ``CPAN``. 


Compile 
--------

Clone the git repository on your computer with the following command: 

    git clone git@github.com:cbcrg/tcoffee.git tcoffee
    
    
Make sure you have installed the required dependencies listed above. 
When done, move in the project root folder named ``tcoffee`` and enter the 
following commands:     
    
    $ cd compile
    $ make t_coffee
    

The binary will be automatically copied to the path specified by the environment 
variable ``$USER_BIN`` (check that it exists before run the make command).

