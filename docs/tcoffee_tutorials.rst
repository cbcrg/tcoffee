*******************
T-Coffee Tutorials
*******************

Using the T-Coffee Package to Build Multiple Sequence Alignments of Proteins, RNAs, DNA 
sequences and 3D-Structures.

Introduction
---------------

T-Coffee is a versatile Multiple Sequence Alignment method suitable for aligning most types 
of biological sequences. The series of protocols presented here show how the package can be 
used to multiply align proteins, DNA and RNA sequences. 

The protein section presents controlled cases for PSI-Coffee the homology extended mode suitable 
for remote homologues, Expresso the structure based multiple aligner and M-Coffee, a meta-version 
able to combine several third party aligners into one. 

We then show how the T-RMSD option can be used to produce a functionally informative structure 
based clustering. RNA alignment procedures are shown for R-Coffee a mode that produces secondary 
structure based MSAs. DNA alignments are illustrated with Pro-Coffee, a multiple aligner specific 
of promoter regions. The last sections presents the many reformatting utilities bundled with T-Coffee. 

The package is an open-source freeware available from `www.tcoffee.org <http://www.tcoffee.org>`_.


Materials
-----------

The list of files (input and output) required by this protocol is available 
from `here <http://www.tcoffee.org/Packages/NatureProtocols/NatureProtocolDataset.tar.gz>`_. 

They can be automatically retrieved using the following command::

    t_coffee -other_pg nature_protocol.pl    

This will create 4 repertories containing the input sequences necessary for the protocols 
we report in this section. For each part, all command lines have been collected into the file README.sh.

Procedures
------------

- Installation
- Protein Multiple Sequence Alignments
- RNA Multiple Sequence Alignments
- Promoter alignments
- Reformat alignments