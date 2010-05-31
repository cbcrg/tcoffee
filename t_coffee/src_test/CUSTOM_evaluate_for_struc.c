#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>


#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"

#include "dp_lib_header.h"


/*
  23/06/00, Cedric Notredame
  
  
  1-Content of the data structures.
  2-Implementing your own function in pdb_align.
  3-Using that function with T-Coffee (multiple Sequence Alignment).
  4-Syntax rules as defined by Philipp Bucher (19/06/00).
  5-Current Shortcomings
  6-Enquiries.

  1-Content of the data structures

   This file only contains a dummy function to help you create 
   your own matching potential function (Step 2 in the Notations RULES
 
   int evaluate_match_score ( Constraint_list *CL, int A, int i, int B, int j)
 
   returns a score, expected to be between -100 and 100, that corresponds to the matching of 
   A_i with B_j.

   Most needed parameters are included in the data structure CL,
   This Data Structure is declared in util_constraint_list.h
   The following, non exhaustive list explains the most common parameters

   The neighborhood is computed using:
   ((CL->T[A])->pdb_param)->maximum_distance as a radius for the Bubble
   ((CL->T[A])->pdb_param)->n_excluded_nb are excluded around the central residue
                           i.e i-1 and i+1 for n_excluded_nb=1.
			   

   ((CL->T[A])->Bubble)->nb[i][0]     --> Number of residues in the bubble around A_i
   ((CL->T[A])->Bubble)->nb[i][k]=j   --> Index of the kth residue in the bubble
                                         Residues are sorted according to the Ca chain
   ((CL->T[A])->Bubble)->d_nb[i][k]=d --> Distance between A_i and A_j equals d;
				     
   ((CL->T[A])->ca[i]->x -----------> Coordinates of the Ca A_i
   ((CL->T[A])->ca[i]->y
   ((CL->T[A])->ca[i]->z



   ((CL->T[A])->len      -----------> Length of Chain A.
   ((CL->T[A])->n_atom   -----------> n atoms in A.
  
  
   Unspecified parameters can be passed from the command line:

           align_pdb -extra_parameters=10, 10.3, 11, 12.4, my_file
	   
   The values of these parameters can be accessed in:
  
    ((CL->T[A])->pdb_param)->n_extra_param=4
    ((CL->T[A])->pdb_param)->extra_param[0]="10"
    ((CL->T[A])->pdb_param)->extra_param[1]="10.3"
    ((CL->T[A])->pdb_param)->extra_param[2]="11.6"
    ((CL->T[A])->pdb_param)->extra_param[3]="my_file"

   These parameters contain strings! To get the real values, in C, use atoi and atof:
    atoi ( ((CL->T[A])->pdb_param)->extra_param[0])=10;
    atof ( ((CL->T[A])->pdb_param)->extra_param[1])=10.3;
    
    The maximum number of parameters is currently 1000...

   

  2-Implementing your own function

    all you need to do is to edit this file and recompile align_pdb.
    There is no need to prototype any function.

    10 functions holders exist, that correspond to the 10 dummy functions 
    declared in this file:
                  custom_pair_score_function1
		  custom_pair_score_function2
		  custom_pair_score_function3
		  custom_pair_score_function4
		  .....
		  custom_pair_score_function10
    
    Let us imagine, you want to use custom_pair_function1.
    
             1-In CUSTOM_evaluate_for_struc.c, modify custom_pair_function1,
	       so that it computes the score you need.
	
             2-If you need extra parameters, get them from ((CL->T[A])->pdb_param)->extra_param.
	     3-Recompile pdb_align:
	               -put it in your bin
		       -rehash or whatever
	       
	     4-run the program as follows:

	     align_pdb -in <struc1> <struc2> -hasch_mode=custom_pair_score_function1
	               -extra_param=10, 12, 0.4, matrix...

	     5-My advice for a first time: make a very simple dummy function that spits
	     out the content of extra_param.

	     6-Remember it is your responsability to control the number of extra parameters
	     and their proper order, and type. Do not forget to apply atoi and atof to the parameters 
	    
	     7-Remember that the modifications you made to CUSTOM_evaluate_for_sytructure
	     must be preserved by you!!! They may disappear if you update align_pdb, save them 
	     appart as your own patch.

	    
 
    
3-Using that function with T-Coffee (multiple Sequence Alignment).
  
  1- setenv ALIGN_PDB_4_TCOFFEE  <your version of align_pdb>
  
  2- run t_coffee
  To do so, you will NOT NEED to recompile T-Coffee, simply type:
            t_coffee -in <struc1> <struc2> ... custom1_align_pdb_pair
  
	    
  
4-Syntax rules as defined by Philipp Bucher (19/06/00).  
 
  Proposed ascii text notation for align_pdb 

     First, let us summarize the align pdb algorithm in plain 
     english: 

     Given are two protein structures A and B.

        Step 1: For each residue in each structure extract 
           the local structural neighbourhood. A neighbourhood
           is simply a subset of (usually non-consecutive)
           residues from one of the structures. 

        Step 2: For all possible pairs of residues between structures 
           A and B, compute the optimal neighbourhood alignment 
           score. This score, which is also referred to as 
           local neighbourhood similarity (LNS) score indicates
           whether two residues have similar local stuctural
           environemnts.

        Step 3: Generate one (or multiple) optimal structural alignment(s)
           for A and B based on LNS scores plus some gap penalty
           function. 

     Now, some rules for ascii/email notation:

      - Whenever possible use a style which fits on one line (because it 
        is painful to modify formulas that span over several lines). Example: 

        Use: ( a**2 + b**2 )**0.5  
                      ________
                     |  2    2
        instead of: \| a  + b

        Introduce local variables/functions to split long expressions over 
        several lines, e.g. 

        Score = Sum(0<i<I+1) Match(A_i,B_i) where 

                Match(A,B) = ..

      - Pseudosubscript notation for (multiply) indexed variables: 

        A_i A_j_i

        As a general rule, I propose that we always use lower case 
        letters for indices, and that the corresponding upper case letters
        denote the number of indexed objects. 

     Index usage conventions: I propose that we use different indices 
     for different objects: 

        i  a residue of structure A 
        j  a residue of structure B
        k  a residue of a neighbourhood of structure A 
        l  a residue of a neighbourhood of structure B
        m  a residue pair of a neignourhood alignment 
        n  a residue pair of a structure alignment 
     
     Some examples, extensions, and additional conventions: 
     
        I      # of residues in structure A 
        K_i    # of residues of the neighbourhood of residue i of  structure A
        A_k(i) # the kth residue of neighbourhood i. 
        M_i,j  # of residue pairs of an alignment of the neighbourhoods of
               residues i and j. 

     The pseudosubscript notation may not always be optimal in terms clarity. 
     We may occasionally use parenthesis, comma-separated susbscripts, etc.
     instead, e.g. M(i,j) or M_ij. 

     The residues of the structures will be denoted: 

        A_1, A_2 ... A_I
        B_1, B_2 ... B_J
   
     This is for expressing general concepts only. It is of little practical
     importance for the moment since we do not use all residue-related 
     structural information from pdb. Instead we use the C-alpha coordinates 

        C_1, C_2 ... C_I  (for protein C)
        D_1, D_2 ... D_J  (for protein B)

     for all compututations. The D_1 ...  is admittedly not very intuitive
     and I'm open for other suggestions. For the scalar components of the
     C-alpha coordinates I propose that we use 

        C_1 = CX_1, CY_1, CZ_1 = (for example) 7.51, 1.24, 3.01

     For the distance between two C-alpha atoms we write 

        |C_1-C_2| 

     which equals 

        [ (CX_1-CX_2)**2 + (CY_1-CY_2)**2 + (CZ_1-CZ_2)**2 ]**0.5
   
     if I remember correctly from high school.

     Back to the algorithm: 

     > Step 1: For each residue in each structure, extract 
     >         the local structural neighbourhood. A neighbourhood
     >         is simply a subset of (usually non-consecutive)
     >         residues from the same structure. 

     The result is something like: 

       P(i) = P_1(i) .. P_k(i) .. P_K_i(i)  
       Q(i) = Q_1(j) .. Q_l(j) .. Q_L_j(j)  

     These are all ordered integer arrays. The P's and Q's indicate
     residue positions in sequence space. For the C-alpha coordinates,
     we use: 

       C(i) = C_1(i) .. C_k(i) .. C_K_i(i)
       D(i) = D_1(j) .. D_l(j) .. D_L_i(j)

     > Step 2: For all possible pairs of residues between structures 
     >         A and B, compute the optimal neighbourhood alignment 
     >         score. This score, which is also referred to as 
     >         local neighbourhood similarity (NSL) score indicates
     >         whether two residues have similar local stuctural
     >         environemnts.

     We have to define a similarity score: 

       S(i,j) = function[A,B,P(i),Q(j)]

     More specifically, S(i,j) is the score of an opimal alignment between 
     two subsets of C-alpha coordinates from A and B, defined by P(i) and Q(j).
     We use the following notation for an alignment between two neighbourhoods. 

       R = (k_1,l_1) .. (k_m,l_m) .. (k_M, l_M) 

     This is pretty abstract and requires some explanation.

     The alignment consists of M pairs of residues from two neighbourhoods. 
     The paired residues are numbered 1,2...K and 1,2...L, respectively. 
     Obviously M <= K,L. For K=9 and L=7, a possible alignment would
     look as follows: 

       R = (1,2) , (2,3) , (5,4) , (6,5) , (9,7) 
     
     This alignment consists of five paired residues, the first
     residue of neighbourhood P(i) is aligned with with the second residue
     of Q(j), etc.  

     The score of an alignment Z(R) is a function that can be
     defined in many different ways. But independently of its 
     definition: 

        S(i,j) = Z(R*,A,B,P(i),Q(j))
        R* = argmax Z(R,A,B,P(i),Q(j))

     This is just a complicated way of saying that the LNS score
     S(i,j) is an optimal alignment score. A simple alignment
     scoring function would be: 

        Z = Sum(m=1..M) [ 2 - |C_(k_m) - D_l_m)| ] 

     A more complex function could be the sum of the sums of "pair-weights",
     "pair-connection-weights", and unpaired-residue-weights": 

        Z =   Sum(m=1 .. M)  [ PW (i,P_k_m,Q_l_m,C_k_m, D_l_m) ]
            + Sum(m=2 .. M ) [ PCW(j,P_k_m,P_l_m,Q_k_m-1,Q_l_m-1,C_k_m,D_l_m,C_k_m-1,D_j_m-1 ]
            + Sum(over k for all C_i(k) unpaired) UPRW [P_k, C_k ]
            + Sum(over l for all C_i(l) unpaired) UPRW [Q_l, D_l)) ]

     Here, the terms P_k_m ... denote sequence positions, the terms C_k_m ...
     denote coordinates. i and j, the sequence position of the center residues
     of the neighbourhoods under consideration) are included in the argument
     lists of the functions because they are necessary to decide whether
     a residue A_k_m occurs before or after the residue A_i in sequence space.
     We don not want to align a residue A_k_m that occurs before A_i with
     a residue B_j_l that occurs after B_j and vice-versa.
     
     The LNS score could also be defined by a recursive equation system 
     defining a dynamic programming algorithm. However, I find the 
     above formulation more helpful for designing appropriate alignment 
     scoring functions. 

     >       Step 3: Generate one (or multiple) optimal structural alignment(s)r
     >          for A and B based on NLS scores plus some gap penalty
     >          function. 

     This is now pretty simple. We use essentially the same notation as 
     for the neighbourhood alignments. 

       R = (i_1,j_1) .. (i_n,j_n) .. (i_N, j_N) 
     
       X* = X(R*,A,B)
       R* = argmax X(R,A,B)

     The alignment scoring functing X is the sum of the LNS scores 
     of the pairs minus some gap penalties.  
    

5-Current Shortcomings

   At present, it is impossible to use the extra_param flag with T-Coffee. This means that 
   the actual values of your parameters must be HARD-CODED within the custom_pair_score_function
   you are using. 

   On request, I will implement a solution to solve that problem.

6-Contact
  For any enquiry, please send a mail to:
     cedric.notredame@europe.com
 */








int custom_pair_score_function1 (Constraint_list *CL, int s1, int r1, int s2, int r2)
   {
     int score=0;
     int a;
     FILE *fp;

     fp=vfopen ( "test_file", "w");
     for ( a=0; a< ((CL->T[0])->pdb_param)->n_extra_param; a++)
       fprintf (fp, "\n\t%s", ((CL->T[0])->pdb_param)->extra_param[a]);

     fprintf ( fp, "\nTEST OK");
     vfclose ( fp);
     exit (1);

     return score;
     
   }
int custom_pair_score_function2 (Constraint_list *CL, int s1, int r1, int s2, int r2)
   {
     int score=0;
     
     return score;
     
   }
int custom_pair_score_function3 (Constraint_list *CL, int s1, int r1, int s2, int r2)
   {
     int score=0;
     
     return score;
     
   }
int custom_pair_score_function4 (Constraint_list *CL, int s1, int r1, int s2, int r2)
   {
     int score=0;
     
     return score;
     
   }
int custom_pair_score_function5 (Constraint_list *CL, int s1, int r1, int s2, int r2)
   {
     int score=0;
     
     return score;
     
   }
int custom_pair_score_function6 (Constraint_list *CL, int s1, int r1, int s2, int r2)
   {
     int score=0;
     
     return score;
     
   }
int custom_pair_score_function7 (Constraint_list *CL, int s1, int r1, int s2, int r2)
   {
     int score=0;
     
     return score;
     
   }
int custom_pair_score_function8 (Constraint_list *CL, int s1, int r1, int s2, int r2)
   {
     int score=0;
     
     return score;
     
   }

int custom_pair_score_function9 (Constraint_list *CL, int s1, int r1, int s2, int r2)
   {
     int score=0;
     
     return score;
     
   }
int custom_pair_score_function10 (Constraint_list *CL, int s1, int r1, int s2, int r2)
   {
     int score=0;
     
     return score;
     
   }
