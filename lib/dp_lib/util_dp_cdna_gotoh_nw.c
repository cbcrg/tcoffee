#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>


#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"




int gotoh_cdna_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL)
	{
	  /*TREATMENT OF THE TERMINAL GAP PENALTIES*/
	  /*TG_MODE=0---> gop and gep*/
	  /*TG_MODE=1---> ---     gep*/
	  /*TG_MODE=2---> ---     ---*/


	int maximise;
	int l1, l2;	
	/*VARIABLES FOR THE MULTIPLE SEQUENCE ALIGNMENT*/
	int **tot_diag;
	int *diag;
	int ktup;
	static int n_groups;
	static char **group_list;
	int score;
	Alignment *B;
        /********Prepare Penalties******/
	
	
	maximise=CL->maximise;
	ktup=CL->ktup;
	/********************************/
	
	
	
	if ( !group_list)
	   {
	     
	       group_list=make_group_aa (&n_groups, CL->matrix_for_aa_group);
	   }
	l1=strlen (A->seq_al[l_s[0][0]]);
	l2=strlen (A->seq_al[l_s[1][0]]);
	B=dna_aln2_3frame_cdna_aln(A,ns,l_s); 
	tot_diag=evaluate_diagonals_cdna( B, ns, l_s, CL, maximise,n_groups,group_list, ktup);
	diag=extract_N_diag (l1, l2, tot_diag,l1+l2,0);
	score=make_fasta_cdna_pair_wise ( A,B, ns, l_s, CL, diag);
	free_int (tot_diag, -1);
	free_aln (B);
	free (diag);
	return score;
	
	}
     

