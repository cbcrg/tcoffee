#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>


#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"

int ** evaluate_diagonals ( Alignment *A, int *ns, int **l_s, Constraint_list *CL,int maximise,int n_groups, char **group_list, int ktup)
        {
	int **tot_diag;

	if      ( CL->L)
	    tot_diag=evaluate_diagonals_with_clist ( A, ns, l_s, CL, maximise,n_groups,group_list, ktup);
	else if ( CL->use_fragments)
	    tot_diag=evaluate_segments_with_ktup ( A, ns, l_s, CL, maximise,n_groups,group_list, ktup);
	else
	    tot_diag=evaluate_diagonals_with_ktup ( A, ns, l_s, CL, maximise,n_groups,group_list, ktup);
	
	return tot_diag;
	}
int ** evaluate_segments_with_ktup ( Alignment *A, int *ns, int **l_s, Constraint_list *CL,int maximise,int n_groups, char **group_list, int ktup)
    {
   /*
    Reads in an alignmnet A, with two groups of sequences marked.
    1-Turn each group into a conscensus, using the group list identifier.
               -if the group list is left empty original symbols are used
    2-hasSgroupch the two sequences
    3-score each diagonal, sort the list and return it (diag_list)
   */

    char *seq1, *seq2, *alphabet=NULL;
    int a,b,l1, l2, n_ktup,pos_ktup1, pos_ktup2, **pos;
    int *hasched_seq1, *hasched_seq2,*lu_seq1,*lu_seq2;
    int n_diag, **diag, current_diag, **dot_list, n_dots, cost;
    static char *buf;
    
    

    pos=aln2pos_simple ( A,-1, ns, l_s);
    seq1=aln2cons_seq (A, ns[0], l_s[0], n_groups, group_list);
    seq2=aln2cons_seq (A, ns[1], l_s[1], n_groups, group_list);

    

    alphabet=get_alphabet (seq1,alphabet);
    alphabet=get_alphabet (seq2,alphabet);

    

    l1=strlen ( seq1);
    l2=strlen ( seq2);

    n_diag=l1+l2-1;
    diag=declare_int ( n_diag+2, 3);
    n_ktup=(int)pow ( (double)alphabet[0]+1, (double)ktup);
    
    hasch_seq(seq1, &hasched_seq1, &lu_seq1,ktup, alphabet);
    hasch_seq(seq2, &hasched_seq2, &lu_seq2,ktup, alphabet);
    
    
    
    /*EVALUATE THE DIAGONALS*/
    for ( a=0; a<= n_diag; a++)diag[a][0]=a;
    for ( n_dots=0,a=1; a<= n_ktup; a++)
        {
	    pos_ktup1=lu_seq1[a];
	    while (TRUE)
	          {
		  if (!pos_ktup1)break;
		  pos_ktup2=lu_seq2[a];
		  while (pos_ktup2) 
		            {
			    n_dots++;
			    pos_ktup2=hasched_seq2[pos_ktup2];
			    }
		  pos_ktup1=hasched_seq1[pos_ktup1];
		  }
	}

    if ( n_dots==0)
       {
	   if ( !buf)buf=vcalloc (30, sizeof (char));
	   sprintf ( buf, "abcdefghijklmnopqrstuvwxyz");
	   return evaluate_segments_with_ktup ( A,ns,l_s, CL,maximise,1,&buf,1);
       }
	       
    dot_list=declare_int ( n_dots,3);
    
    for ( n_dots=0,a=1; a<= n_ktup; a++)
        {
	    pos_ktup1=lu_seq1[a];
	    while (TRUE)
	          {
		  if (!pos_ktup1)break;
		  pos_ktup2=lu_seq2[a];
		  while (pos_ktup2) 
		            {
			    current_diag=(pos_ktup2-pos_ktup1+l1);
			    dot_list[n_dots][0]=current_diag;
			    dot_list[n_dots][1]=pos_ktup1;
			    dot_list[n_dots][2]=pos_ktup2;			    
			    pos_ktup2=hasched_seq2[pos_ktup2];
			    n_dots++;
			    }
		  pos_ktup1=hasched_seq1[pos_ktup1];
		  }
	}
    
    
    
    hsort_list_array ((void **)dot_list, n_dots, sizeof (int), 3, 0, 3);

    
        
    current_diag= (int)dot_list[0][0];
    for ( b=0; b< ktup; b++)diag[current_diag][2]+=(CL->get_dp_cost) ( A, pos, ns[0], l_s[0], dot_list[0][1]+b-1, pos,ns[1], l_s[1], dot_list[0][2]+b-1, CL);
    
    for ( a=1; a< n_dots; a++)
        {
	    
	    for ( cost=0, b=0; b< ktup; b++)cost+=(CL->get_dp_cost) ( A, pos, ns[0], l_s[0], dot_list[a][1]+b-1, pos,ns[1], l_s[1], dot_list[a][2]+b-1, CL);
	    if ((dot_list[a][0]!=dot_list[a-1][0]) || !((dot_list[a][1]==dot_list[a-1][1]+1)&&(dot_list[a][2]==dot_list[a-1][2]+1)))
	       {
		
		   diag[current_diag][1]=best_of_a_b(diag[current_diag][2], diag[current_diag][1], 1);
		   diag[current_diag][2]=0;
		   current_diag=dot_list[a][0];
	       }
	    diag[current_diag][2]+=cost;
	}
    diag[current_diag][1]=best_of_a_b(diag[current_diag][2], diag[current_diag][1], 1);
    sort_int (diag+1, 3, 1,0, n_diag-1);
    
   
    free (seq1);
    free (seq2);
    free (alphabet);
    free (hasched_seq1);
    free (hasched_seq2);
    free (lu_seq1);
    free (lu_seq2);
    
    free_int (dot_list, -1);
    return diag;
    }
int ** evaluate_diagonals_with_ktup ( Alignment *A, int *ns, int **l_s, Constraint_list *CL,int maximise,int n_groups, char **group_list, int ktup)
    {
   /*
    Reads in an alignmnent A, with two groups of sequences marked.
    1-Turn each group into a conscensus, using the group list identifier.
               -if the group list is left empty original symbols are used
    2-hasch the two sequences
    3-score each diagonal, sort the list and return it (diag_list)
   */

    char *seq1, *seq2, *alphabet=NULL;
    int a,b,l1, l2, n_ktup,pos_ktup1, pos_ktup2, **pos;
    int *hasched_seq1, *hasched_seq2,*lu_seq1,*lu_seq2;
    int n_diag, **diag, current_diag, n_dots;
    static char *buf;
    pos=aln2pos_simple ( A,-1, ns, l_s);
    seq1=aln2cons_seq (A, ns[0], l_s[0], n_groups, group_list);
    seq2=aln2cons_seq (A, ns[1], l_s[1], n_groups, group_list);
    

    alphabet=get_alphabet (seq1,alphabet);
    alphabet=get_alphabet (seq2,alphabet);

    l1=strlen ( seq1);
    l2=strlen ( seq2);

    n_diag=l1+l2-1;
    diag=declare_int ( n_diag+2, 3);
    n_ktup=(int)pow ( (double)alphabet[0]+1, (double)ktup);
    

    hasch_seq(seq1, &hasched_seq1, &lu_seq1,ktup, alphabet);
    hasch_seq(seq2, &hasched_seq2, &lu_seq2,ktup, alphabet);
    
    
    
    /*EVALUATE THE DIAGONALS*/
    for ( a=0; a<= n_diag; a++)diag[a][0]=a;
    for ( n_dots=0,a=1; a<= n_ktup; a++)
        {
	    pos_ktup1=lu_seq1[a];
	    while (TRUE)
	          {
		  if (!pos_ktup1)break;
		  pos_ktup2=lu_seq2[a];
		  while (pos_ktup2) 
		            {
			    current_diag=(pos_ktup2-pos_ktup1+l1);
			    for ( b=0; b< ktup; b++)
			        {
		         	    diag[current_diag][1]+=(CL->get_dp_cost) ( A, pos, ns[0], l_s[0], pos_ktup1+b-1, pos,ns[1], l_s[1], pos_ktup2+b-1, CL);
				    n_dots++;
				}
			    pos_ktup2=hasched_seq2[pos_ktup2];
			    }
		  pos_ktup1=hasched_seq1[pos_ktup1];
		  }
	  
	}
    if ( n_dots==0)
       {
	   if ( !buf)
	       {
	       buf=vcalloc ( 30, sizeof (30));
	       sprintf ( buf, "abcdefghijklmnopqrstuvwxyz");
	       }
	    free ( hasched_seq1);
	    free ( hasched_seq2);
	    free (lu_seq1);
	    free (lu_seq2);
	   return evaluate_diagonals_with_ktup ( A,ns,l_s, CL,maximise,1,&buf,1);
       }
    
   
    /*fprintf ( stderr, "\nDOTS: %.3f", (float)n_dots/(float)(l1*l2));*/
    sort_int (diag+1, 2, 1,0, n_diag-1);
    
    free (seq1);
    free (seq2);
    free (alphabet);
    free ( hasched_seq1);
    free ( hasched_seq2);
    free (lu_seq1);
    free (lu_seq2);
    return diag;
    }


int ** evaluate_diagonals_with_clist ( Alignment *A, int *ns, int **l_s, Constraint_list *CL,int maximise,int n_groups, char **group_list, int ktup)
    {

   /*
    Reads in an alignmnent A, with two groups of sequences marked.
    Weight the diagonals with the values read in the constraint list        
   */

    int l1, l2,n_diag, s1, s2, r1, r2;
    int a, b, c, d;
    int **diag;
    int **code;
    int **pos;
    static int *entry;
    
    if ( !entry)entry=vcalloc ( CL->entry_len, CL->el_size);
    CL=index_constraint_list (CL);
    
    l1=strlen (A->seq_al[l_s[0][0]]);
    l2=strlen (A->seq_al[l_s[1][0]]);

    n_diag=l1+l2-1;
    diag=declare_int ( n_diag+2, 3);
    for ( a=0; a<= n_diag; a++)diag[a][0]=a;

    A->S=CL->S;
    code=seq2aln_pos (A, ns, l_s);
    pos =aln2pos_simple ( A,-1, ns, l_s);

    for (a=0; a<ns[0]; a++)

        {
	s1=A->order[l_s[0][a]][0];
	for (b=0; b<ns[1]; b++)
	    {
	    s2=A->order[l_s[1][b]][0];
	    for (c=CL->start_index[s1][s2], d=0; c<CL->end_index[s1][s2];c++, d++)
	        {
		entry=extract_entry ( entry,c, CL);
		if (s1==entry[SEQ1])
		    {
		    r1=code [s1][entry[R1]];
		    r2=code [s2][entry[R2]];
		    }
		else if ( s2==entry[SEQ1])
		    {
		    r2=code [s2][entry[R1]];
		    r1=code [s1][entry[R2]];
		    }
				

		diag[(r2-r1+l1)][1]+=(CL->get_dp_cost) ( A, pos, ns[0], l_s[0],r1-1, pos,ns[1], l_s[1], r2-1, CL);
		}			
	    }
	}
    

    sort_int (diag+1, 2, 1,0, n_diag-1);
    
    free_int (code,-1);
    
    return diag;
    }

int * flag_diagonals (int l1, int l2, int **sorted_diag, float T)
    {
    int a, b, up, low,current_diag,n_diag;
    int * slopes;
    int *diag_list;
    double mean;
    double sd;
    int window=0;

    

    
    n_diag=l1+l2-1;
    mean=return_mean_int ( sorted_diag, n_diag+1, 1);
    sd  =return_sd_int ( sorted_diag, n_diag+1, 1, (int)mean);
    
    if ( T==0)T=(((double)sorted_diag[n_diag][1]-mean)/sd)/10;
    

    diag_list=vcalloc (l1+l2+1, sizeof (int));
    slopes=vcalloc ( n_diag+1, sizeof (int));
 
    for ( a=n_diag; a>0; a--)
            {
	    current_diag=sorted_diag[a][0];
	    if (((double)sorted_diag[a][1]-mean)/sd>T)
	       {
		   up=MAX(1,current_diag-window);
		   low=MIN(n_diag, current_diag+window);
		   for ( b=up; b<=low; b++)slopes[b]=1;
	       }
	    else break;
		
	    }
    for ( a=1, b=0; a<=n_diag; a++)
        {
	    b+=slopes[a];
	}
/*    fprintf (stderr, "\nN_DIAG=%d [%.3f]\n]", b, (float)b/((float)n_diag));*/
    slopes[1]=1;
    slopes[l1+l2-1]=1;
    slopes[l2]=1;
    for (a=0; a<= (l1+l2-1); a++)
	if ( slopes[a]){diag_list[++diag_list[0]]=a;}

    free (slopes);
    
    return diag_list;
    }
int * extract_N_diag (int l1, int l2, int **sorted_diag, int n_chosen_diag)
    {
    int a, b, up, low,current_diag,n_diag;
    int * slopes;
    int *diag_list;
    int window=0;


    

    
    n_diag=l1+l2-1;
    
    diag_list=vcalloc (l1+l2+1, sizeof (int));
    slopes=vcalloc ( n_diag+1, sizeof (int));
 
   

    for ( a=n_diag; a>0 && a>(n_diag-n_chosen_diag); a--)
            {
	    current_diag=sorted_diag[a][0];
	    up=MAX(1,current_diag-window);
	    low=MIN(n_diag, current_diag+window);
	    for ( b=up; b<=low; b++)slopes[current_diag]=1;
	    }
    
    slopes[1]=1;
    slopes[l1+l2-1]=1;
    slopes[l2]=1;
    for (a=0; a<= (l1+l2-1); a++)
	if ( slopes[a]){diag_list[++diag_list[0]]=a;}

    free (slopes);
    return diag_list;
    }

int cfasta_gotoh_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL)
    {
/*TREATMENT OF THE TERMINAL GAP PENALTIES*/
/*TG_MODE=0---> gop and gep*/
/*TG_MODE=1---> ---     gep*/
/*TG_MODE=2---> ---     ---*/


	int maximise;
		
/*VARIABLES FOR THE MULTIPLE SEQUENCE ALIGNMENT*/
	int **tot_diag;
      
	int *diag;
	int ktup;
	static int n_groups;
	static char **group_list;
	int score, new_score;
        int n_chosen_diag=0;
        int step;
	int max_n_chosen_diag;
	int l1, l2;
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
	
	if ( !CL->fasta_step)
	    {
	    step=MIN(l1,l2);
	    step=(int) log ((double)MAX(step, 1));
	    step=MAX(step, 5);
	    }
	else
	    {
		step=CL->fasta_step;
	    }
	
	

	tot_diag=evaluate_diagonals ( A, ns, l_s, CL, maximise,n_groups,group_list, ktup);

	max_n_chosen_diag=strlen (A->seq_al[l_s[0][0]])+strlen (A->seq_al[l_s[1][0]])-1;
	/*max_n_chosen_diag=(int)log10((double)(l1+l2))*10;*/
	
	n_chosen_diag+=step;
	n_chosen_diag=MIN(n_chosen_diag, max_n_chosen_diag);
	

	diag=extract_N_diag (strlen (A->seq_al[l_s[0][0]]),strlen (A->seq_al[l_s[1][0]]), tot_diag, n_chosen_diag);
	
	score    =make_fasta_gotoh_pair_wise ( A, ns, l_s, CL, diag);
	new_score=0;

	free ( diag);
	
 	while (new_score!=score && n_chosen_diag< max_n_chosen_diag )
	      {

		 score=new_score;
		 ungap_sub_aln ( A, ns[0], l_s[0]);
		 ungap_sub_aln ( A, ns[1], l_s[1]);
		

		 n_chosen_diag+=step;
		 n_chosen_diag=MIN(n_chosen_diag, max_n_chosen_diag);

		 diag     =extract_N_diag (strlen (A->seq_al[l_s[0][0]]),strlen (A->seq_al[l_s[1][0]]), tot_diag, n_chosen_diag); 
		 new_score=make_fasta_gotoh_pair_wise (  A, ns, l_s, CL, diag);
		 free ( diag);
	      }
	
	score=new_score;
	free_int (tot_diag, -1);

	return score;
    }

int fasta_gotoh_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL)
    {
/*TREATMENT OF THE TERMINAL GAP PENALTIES*/
/*TG_MODE=0---> gop and gep*/
/*TG_MODE=1---> ---     gep*/
/*TG_MODE=2---> ---     ---*/


	int maximise;
		
/*VARIABLES FOR THE MULTIPLE SEQUENCE ALIGNMENT*/
	int **tot_diag;
	int *diag;
	int ktup, diagonal_threshold;
	static int n_groups;
	static char **group_list;
	int score;
        /********Prepare Penalties******/
	
	
	maximise=CL->maximise;
	ktup=CL->ktup;
	diagonal_threshold=CL->diagonal_threshold;
	/********************************/
	
	
	
	if ( !group_list)
	   {
	     
	       group_list=make_group_aa (&n_groups, CL->matrix_for_aa_group);
	   }
	
	
	tot_diag=evaluate_diagonals ( A, ns, l_s, CL, maximise,n_groups,group_list, ktup);
	diag=flag_diagonals ( strlen(A->seq_al[l_s[0][0]]), strlen(A->seq_al[l_s[1][0]]), tot_diag,diagonal_threshold);
	score=make_fasta_gotoh_pair_wise ( A, ns, l_s, CL, diag);
	
	free_int (tot_diag, -1);
	free (diag);
	return score;
    }

int make_fasta_gotoh_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL, int *diag)
    {
/*TREATMENT OF THE TERMINAL GAP PENALTIES*/
/*TG_MODE=0---> gop and gep*/
/*TG_MODE=1---> ---     gep*/
/*TG_MODE=2---> ---     ---*/


	int TG_MODE, gop, l_gop, gep,l_gep, maximise;
		
/*VARIABLES FOR THE MULTIPLE SEQUENCE ALIGNMENT*/
	int a, b,c,k, t;
	int l1, l2,eg, ch, sub,score=0, fop,last_i, last_j, i, delta_i, j, pos_j, ala, alb, LEN, n_diag, match1, match2;
	
	int **C, **D, **I, **trace, **pos0, **ddg;
	int lenal[2], len;
	char *buffer, *char_buf;
	char **aln, **al;
			
        /********Prepare Penalties******/
	gop=CL->gop*SCORE_K;
	gep=CL->gep*SCORE_K;
	TG_MODE=CL->TG_MODE;
	maximise=CL->maximise;
	

	/********************************/
	
	
	
		
        n_diag=diag[0];

    
       
       l1=lenal[0]=strlen (A->seq_al[l_s[0][0]]);
       l2=lenal[1]=strlen (A->seq_al[l_s[1][0]]);
       
       

	/*diag:
	  diag[1..n_diag]--> flaged diagonal in order;
	  diag[0]=0--> first diagonal;
	  diag[n_diag+1]=l1+l2-1;
	*/    
	
	/*numeration of the diagonals strats from the bottom right [1...l1+l2-1]*/
	/*sequence s1 is vertical and seq s2 is horizontal*/
	/*D contains the best Deletion  in S2==>comes from diagonal N+1*/
	/*I contains the best insertion in S2=> comes from diagonal N-1*/
	

      
       

       C=declare_int (lenal[0]+lenal[1]+1, n_diag+2);
       D=declare_int (lenal[0]+lenal[1]+1, n_diag+2);
       ddg=declare_int (lenal[0]+lenal[1]+1, n_diag+2);
       I=declare_int (lenal[0]+lenal[1]+1, n_diag+2);
       trace=declare_int (lenal[0]+lenal[1]+1, n_diag+2);
       

       al=declare_char (2,lenal[0]+lenal[1]+lenal[1]+1);  
      
       len= MAX(lenal[0],lenal[1])+1;
       buffer=vcalloc ( 2*len, sizeof (char));	
       char_buf= vcalloc (2*len, sizeof (char));	

       pos0=aln2pos_simple ( A,-1, ns, l_s);
       C[0][0]=0;
     
       t=(TG_MODE==0)?gop:0;	
       for ( j=1; j<= n_diag; j++)
	    {
		l_gop=(TG_MODE==0)?gop:0;
		l_gep=(TG_MODE==2)?0:gep;
		
		

		if ( (diag[j]-lenal[0])<0 )
		    {
		    trace[0][j]=UNDEFINED;  
		    continue;
		    }
		C[0][j]=(diag[j]-lenal[0])*l_gep +l_gop;
		D[0][j]=(diag[j]-lenal[0])*l_gep +l_gop+gop;				
	    }
       D[0][j]=D[0][j-1]+gep;


       t=(TG_MODE==0)?gop:0;	
       for ( i=1; i<=lenal[0]; i++)
           {
	        l_gop=(TG_MODE==0)?gop:0;
		l_gep=(TG_MODE==2)?0:gep;

		C[i][0]=C[i][n_diag+1]=t=t+l_gep;
		I[i][0]=D[i][n_diag+1]=t+    gop;
		
		for ( j=1; j<=n_diag; j++)
		    {
			C[i][j]=C[i][0];
			D[i][j]=I[i][j]=I[i][0];
		    }
			
		for (eg=0, j=1; j<=n_diag; j++)
		    {
			pos_j=diag[j]-lenal[0]+i;
			if (pos_j<=0 || pos_j>l2 )
			    {
			    trace[i][j]=UNDEFINED;
			    continue;
			    }
			sub=(CL->get_dp_cost) ( A, pos0, ns[0], l_s[0], i-1, pos0, ns[1], l_s[1],pos_j-1, CL );
			
		    /*1 identify the best insertion in S2:*/
			l_gop=(i==lenal[0])?((TG_MODE==0)?gop:0):gop;
			l_gep=(i==lenal[0])?((TG_MODE==2)?0:gep):gep;
			len=(j==1)?0:(diag[j]-diag[j-1]);
			if ( a_better_than_b(I[i][j-1], C[i][j-1]+l_gop, maximise))eg++;
			else eg=1;			
			I[i][j]=best_of_a_b (I[i][j-1], C[i][j-1]+l_gop, maximise)+len*l_gep;
			

		    /*2 Identify the best deletion in S2*/
			l_gop=(pos_j==lenal[1])?((TG_MODE==0)?gop:0):gop;
			l_gep=(pos_j==lenal[1])?((TG_MODE==2)?0:gep):gep;

			len=(j==n_diag)?0:(diag[j+1]-diag[j]);			
			delta_i=((i-len)>0)?(i-len):0;

			if ( a_better_than_b(D[delta_i][j+1],C[delta_i][j+1]+l_gop, maximise)){ddg[i][j]=ddg[delta_i][j+1]+1;}
			else {ddg[i][j]=1;}
			D[i][j]=best_of_a_b (D[delta_i][j+1],C[delta_i][j+1]+l_gop, maximise)+len*l_gep;


			/*Identify the best way*/	
			score=C[i][j]=best_int ( 3, maximise, &fop, I[i][j], C[i-1][j]+sub, D[i][j]);
			
			
			fop-=1;
			if ( fop<0)trace[i][j]=fop*eg;
			else if ( fop>0 ) {trace[i][j]=fop*ddg[i][j];}
			else if ( fop==0) trace[i][j]=0;

			last_i=i;
			last_j=j;
		    }
	    }


       /*
	            [0][Positive]
	             ^     ^
	             |    /
                     |   /
                     |  /
                     | /
                     |/
       [Neg]<-------[*]
	*/
        
       
	i=last_i;
	j=last_j;



	ala=alb=0;
	match1=match2=0;
	while (!(match1==l1 && match2==l2))
	      {
		  

		  if ( match1==l1)
		     {
			 len=l2-match2;
			 for ( a=0; a< len; a++)
			     {
			     al[0][ala++]=0;
			     al[1][alb++]=1;
			     match2++;
			     }
			 k=0;
			 break;
			 
			 /*k=-(j-1);*/		
			 
		     }
		  else if ( match2==l2)
		     {
			 len=l1-match1;
			 for ( a=0; a< len; a++)
			     {
			     al[0][ala++]=1;
			     al[1][alb++]=0;
			     match1++;
			     }
			 k=0;
			 break;
			 /*k= n_diag-j;*/
		     }
		  else
		      {
			  k=trace[i][j];
		      }
		  
		
		  if ( k==0)
			     {
				 if ( match2==l2 || match1==l1);
				 else 
				    {
					
				    al[0][ala++]=1;
				    al[1][alb++]=1;
				    i--;
				    match1++;
				    match2++;
				    }
			     }
		  else if ( k>0)
			     {
			     
			     len=diag[j+k]-diag[j];
			     for ( a=0; a<len; a++)
			         {
				     if ( match1==l1)break;
				     al[0][ala++]=1;
				     al[1][alb++]=0;
				     match1++;
				 }
			     i-=len;
			     j+=k;
			     }
		  else if ( k<0)
			     {
			     k*=-1;
			     len=diag[j]-diag[j-k];
			     for ( a=0; a<len; a++)
			         {
				     if ( match2==l2)break;
				     al[0][ala++]=0;
				     al[1][alb++]=1;
				     match2++;
				 }
			    
			     
			     j-=k;			    
			     }
	      }

	LEN=ala;
	c=LEN-1;
	invert_list_char ( al[0], LEN);
	invert_list_char ( al[1], LEN);	
	if ( A->declared_len<=LEN)A=realloc_aln2  ( A,A->max_n_seq, 2*LEN);	
	aln=A->seq_al;

	for ( c=0; c< 2; c++)
	    {
	    for ( a=0; a< ns[c]; a++) 
		{		
		ch=0;
		for ( b=0; b< LEN; b++)
		    {		   
		    if (al[c][b]==1)
			char_buf[b]=aln[l_s[c][a]][ch++];
		    else
			char_buf[b]='-';
		   }
		char_buf[b]='\0';
		sprintf (aln[l_s[c][a]],"%s", char_buf);
	        }
	     }
	
	
	A->len_aln=LEN;
	A->nseq=ns[0]+ns[1];
	
	
	free_int (C, -1);
	free_int (D, -1);
	free_int (I, -1);
	free_int (trace, -1);
	free_int (ddg, -1);
	free_char ( al, -1);
	vfree(buffer);
	vfree(char_buf);
	
	return score;
    }


int hasch_seq(char *seq, int **hs, int **lu,int ktup,char *alp)
    {
	static int a[10];
	
	int i,j,l,limit,code,flag;
	char residue;
	
	int alp_lu[10000];
	int alp_size;
	
	alp_size=alp[0];
	alp++;
	
	

	for ( i=0; i< alp_size; i++)
	    {
		alp_lu[alp[i]]=i;
	    }
	
	
	
	l=strlen (seq);
	limit = (int)   pow((double)(alp_size+1),(double)ktup);
	hs[0]=vcalloc ( l+1,sizeof (int));
	lu[0]=vcalloc ( limit+1, sizeof(int));
	

	if ( l==0)exit(0);
	
	for (i=1;i<=ktup;i++)
           a[i] = (int) pow((double)(alp_size+1),(double)(i-1));
	

	for(i=1;i<=(l-ktup+1);++i) 
	        {
		code=0;
		flag=FALSE;
		for(j=1;j<=ktup;++j) 
		   {
		   if (is_gap(seq[i+j-2])){flag=TRUE;break;}
		   else residue=alp_lu[seq[i+j-2]];
		   code+=residue*a[j];
		   }

		if ( flag)continue;
		++code;
		
		if (lu[0][code])hs[0][i]=lu[0][code];
		lu[0][code]=i;
		}
	return 0;
    }
