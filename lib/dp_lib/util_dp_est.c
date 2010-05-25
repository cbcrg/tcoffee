#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"



int evaluate_est_order (Sequence *S, char *concat, Constraint_list *CL, int ktuple)
    {
	int a;
        static char *alphabet;
	int  *hasched_seq, *hasched_seq1, *hasched_seq2;
	int  *lu_seq, *lu_seq1, *lu_seq2;
	int pos_ktup1, pos_ktup2;
	double score=0;
	int n_ktup;
	int n_dots=0;
	
	if ( !alphabet)alphabet=get_alphabet ( concat, alphabet);
	n_ktup=(int)pow ( (double)alphabet[0]+1, (double)ktuple);
	
	hasch_seq (concat,&hasched_seq, &lu_seq,ktuple, alphabet);
	hasched_seq1=hasched_seq2=hasched_seq;
	lu_seq1=lu_seq2=lu_seq;
	


	for ( a=1; a< n_ktup; a++)
	    {
		pos_ktup1=lu_seq1[a];
		
		while (TRUE)
		    {
		    
		     if (!pos_ktup1)break;
		     pos_ktup2=lu_seq2[a];
		     while (pos_ktup2) 
		            {
			    score+=abs ((int)(pos_ktup1-pos_ktup2));
			    pos_ktup2=hasched_seq2[pos_ktup2];
			    n_dots++;
			    }
		  pos_ktup1=hasched_seq1[pos_ktup1];
		  }
	    }
	
	score=(score/(double)(n_dots*strlen(concat)))*100000;
	vfree ( hasched_seq);
	vfree(lu_seq);


	return score;
    }
	
