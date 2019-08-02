##################################
Contacts, Adresses & Contributions
##################################

**********
Contact us
**********
If you have any problem, question, suggestion or else, contact Cedric Notredame by mail at cedric.notredame@crg.eu or at the T-Coffee group tcoffee@googlegroups.com.


********
Adresses
********
The T-Coffee group is currently at the Centre for Genomic Regulation (CRG), programme of Bioinformatics and Genomics (BiG), in Barcelona (PRBB, 88 calle del Doctor Aiguader, 08003, Barcelona). We are always very eager to get some user feedback, so please do not hesitate to drop us a line at cedric.notredame@crg.eu. The latest updates of T-Coffee are always available on our `group webpage <http://www.tcoffee.org>`_ where you will also find links to all the online T-Coffee servers.


*************
Contributions
*************
Main Contributions
==================
T-coffee is developed, maintained, monitored, used and debugged by a dedicated team that include or have included:


Cedric Notredame, Paolo di Tommaso, Jean-François Taly, Cedrik Magis, Fabrice Armougom, Des Higgins, Sebastien Moretti, 
Orla O'Sullivan, Eamon O'Toole, Olivier Poirot, Karsten Suhre, Iain Wallace, Andreas Wilm, Vladimir Keduas.

Other Contributions
===================
We do not mean to steal code, but we will always try to re-use pre-existing code whenever that code exists, free of copyright, 
just like we expect people to do with our code. However, whenever this happens, we make a point at properly citing the source 
of the original contribution. If ever you recognize a piece of your code improperly cited, please drop us a note and we will be 
happy to correct that.

Softwares
---------
In the mean time, here are some important pieces of code from other packages that have been incorporated within the T-Coffee 
package. These include:

 - TMalign package from Zhang, Jeffrey and Skolnik (NAR, 2005, 33:2303).
 - The Sim algorithm of Huang and Miller that given two sequences computes the N best scoring local alignments.
 - The tree reading/computing routines are taken from the ClustalW Package, courtesy of Julie Thompson, Des Higgins and Toby Gibson (Thompson, Higgins, Gibson, 1994, 4673-4680,vol. 22, Nucleic Acid Research).
 - The implementation of the algorithm for aligning two sequences in linear space was adapted from Myers and Miller, in CABIOS, 1988, 11-17, vol. 1).
 - Various techniques and algorithms have been implemented. Whenever relevant, the source of the code/algorithm/idea is indicated in the corresponding function.
 - 64 Bits compliance was implemented by Benjamin Sohn, Performance Computing Center Stuttgart (HLRS), Germany
 - David Mathog (Caltech) provided many fixes and useful feedback for improving the code and making the whole soft behaving more rationally.

An enormous thanks to these people who believe in free open source code.

Bug reports and feedback
------------------------
 - Prof David Jones (UCL) reported and corrected the PDB1K bug (now t_coffee/sap can align PDB sequences longer than 1000 AA).

 - Johan Leckner reported several bugs related to the treatment of PDB structures, insuring a consistent behavior between version 1.37 and current ones.
 
 - Enrico Bonatesta reported several mistakes in the T-Coffee tutorial and provided feedbacks about its readibility.
 
 - ...and to everyone else who send us their feedbacks and suggestions helping us to improve T-Coffee.
 
