#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>


#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"

int fasta_cdna_linear_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL)
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

	tot_diag=evaluate_diagonals_cdna( A, ns, l_s, CL, maximise,n_groups,group_list, ktup);

	
	diag=extract_N_diag (l1, l2, tot_diag,10,3);
	score=make_fasta_cdna_linear_pair_wise ( A, ns, l_s, CL, diag);

	free_int (tot_diag, -1);
	free (diag);
	return score;
    }

int make_fasta_cdna_linear_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL, int *diag)
	{
/*******************************************************************************/
/*                NEEDLEMAN AND WUNSCH (GOTOH)                                 */
/*                                                                             */
/*	makes DP between the the ns[0] sequences and the ns[1]                 */
/*                                                                             */
/*	for MODE, see the function  get_dp_cost                                */
/*******************************************************************************/

	int a, b, l,n;
	

	


/*TREATMENT OF THE TERMINAL GAP PENALTIES*/
/*TG_MODE=0---> gop and gep*/
/*TG_MODE=1---> ---     gep*/
/*TG_MODE=2---> ---     ---*/
	  

	    int TG_MODE;
	    int F_TG_MODE;

	    int gop, gep;
	    int f_gop, f_gep;
	    


/*VARIABLE FOR THE EVALUATION*/
            char *trs1;
            char *trs2;
	    int  **mat;

/*VARIABLES FOR THE MULTIPLE SEQUENCE ALIGNMENT*/
	
	int lenal[2], len, l1, l2;
	int score=0;

/*DP VARIABLES       */
	static int **boundary;
	
	

/*FAST DP VARIABLES       */		
	int pos_j, last_i, last_j, ndiag; 
	int prev_j, j, len_j;
	int prev_i, i, len_i;
/*EVALUATION*/

/*MODEL        VARIABLES       */
	static int nstate;	
	static int **model;
	static int * emission;
	static int * emission_term;
	static int * emission_len;
	static int * emission_len_i;
	static int * emission_len_j;
	
	static int * emission_max_len;
	
	static int * emission_delta_i;
	static int * emission_delta_j;
	
	
	HaschT *Mat;
	HaschT *LMat;
	HaschT *CMat;
	
	static int max_len_diag,len_diag;
	static int max_len,len_al;
	static int mI;
	static int mJ;

	static int M, FM, I, FI, D,FD;
	static int START, END;
	int pc,x,best_pc, t, e, em;

	int Allowed=0;
	int prev_state;
	int prev_score;
	int prev_len;
	int cur_state;

	
	/*DEBUG*/
	int prev[100];
	int state;


	
/********Prepare penalties*******/

	TG_MODE=CL->TG_MODE;
	F_TG_MODE=0;
	gop=CL->gop*SCORE_K;
	gep=CL->gep*SCORE_K;
	
	f_gop=-150;
	f_gep=-100;
	
/********************************/	



	
/*DO MEMORY ALLOCATION FOR DP*/
	ndiag=diag[0];
	


	l1=lenal[0]=strlen (A->seq_al[l_s[0][0]]);
	l2=lenal[1]=strlen (A->seq_al[l_s[1][0]]);
	len= MAX(lenal[0],lenal[1])+1;

	if ( A->cdna_cache) free_int(A->cdna_cache, -1);
	A->cdna_cache=declare_int (1, 2*MAX(l1, l2));

/*PREPARE THE EVALUATION*/      
	if (ns[0]+ns[1]>2)
	  {
	    fprintf ( stderr, "\nERROR: function make_fasta_cdna_pair_wise can only handle two sequences at a time [FATAL:%s]",PROGRAM);
	    crash ("");
	  }
	trs1=translate_dna_seq_on3frame (A->seq_al[l_s[0][0]], 'x', NULL);
	trs2=translate_dna_seq_on3frame (A->seq_al[l_s[1][0]], 'x', NULL);
	mat=CL->M; 
	
	
/*END OF MEMORY ALLOCATION*/
	
	
		/*
		0(s)   +(dd)
		  \      |
		   \     |
		    \    |
		     \   |
		      \  |
		       \ |
		        \|
		-(e)----O
		*/ 
		  
		
	len_al=MAX(lenal[0],lenal[1])+4;
	len_diag=ndiag+4;
	
	max_len=lenal[0]+lenal[1];
	Mat  =hcreate( max_len*2);
	LMat =hcreate( max_len*2);
	CMat =hcreate( max_len*2);

	if ( (len_al>max_len || len_diag>max_len_diag))
	  {
	   
	    
	    free_int (boundary, -1);
	    free_int (model, -1);
	    vfree (emission);
	    vfree (emission_term);
	    vfree (emission_len);
	    vfree (emission_len_i);
	    vfree (emission_len_j);
	    
	    vfree (emission_max_len);
	    
	    vfree (emission_delta_i);
	    vfree (emission_delta_j);
	    nstate=0;
	    max_len_diag=max_len=0;
	 
	  }

	if ( max_len_diag==0 && max_len==0)
	  {	
	    
	    START=nstate++;/*0*/
	    END=nstate++;  /*1*/

	    M=nstate++;    /*2*/
	    I=nstate++;    /*3*/
	    D=nstate++;    /*4*/
	    
	    FM=nstate++;   /*5*/ 
	    FI=nstate++;   /*6*/
	    FD=nstate++;   /*7*/
	    
	    
	    max_len_diag=len_diag;
	    max_len=len_al;
	    
	    mI=max_len*len_diag;
	    mJ=len_diag;
	    
	   
	    
	    
	    model=declare_int (nstate+1, nstate+1);
	    
	    for ( a=0; a<=nstate; a++)
	      for ( b=0; b<= nstate; b++)
		model[a][b]=UNDEFINED;
	    emission=vcalloc ( nstate+1, sizeof (int));
	    emission_term=vcalloc ( nstate+1, sizeof (int));
	    emission_len=vcalloc ( nstate+1, sizeof (int));
	    emission_len_i=vcalloc ( nstate+1, sizeof (int));
	    emission_len_j=vcalloc ( nstate+1, sizeof (int));
	    emission_max_len=vcalloc ( nstate+1, sizeof (int));
	    
	    emission_delta_i=vcalloc ( nstate+1, sizeof (int));
	    emission_delta_j=vcalloc ( nstate+1, sizeof (int));
	    boundary=declare_int (max_len, 2);
	  }

	/*Model Initialization*/
	for ( a=0; a<= nstate; a++)
	  {
	    emission[a]     =UNDEFINED;
	    emission_term[a]=UNDEFINED;
	    emission_len[a]  =1;
	    emission_len_i[a]=1;
	    emission_len_j[a]=1;
	    emission_delta_i[a]=0;
	    emission_delta_j[a]=0;
	  }
	for ( a=0; a<=nstate; a++)
	  for ( b=0; b<= nstate; b++)
	    {
	      model[a][b]=UNDEFINED;
	    }
	

	/*Match*/
	emission_delta_i[M]=-1;
	emission_delta_j[M]= 0;
	emission_len_i[M]=3;
	emission_len_j[M]=3;
	emission[M ]=0;
	emission_max_len[M]=lenal[0]+lenal[1];
	emission_term[M]=0;

	model[M][M]=Allowed;
	model[M][I]=gop;
	model[M][D]=gop;
 	model[M][FM]=f_gop;
        model[M][FI]=f_gop;
        model[M][FD]=f_gop;
	
	model[START][M]=model[M][END]=Allowed;
	/*Insertion*/
	emission_delta_i[I]= 0;
	emission_delta_j[I]=-1;
	emission_len_i[I]=0;
	emission_len_j[I]=3;

	emission[I]=gep;	
	emission_term[I]=(TG_MODE==2)?-emission[I]:0;
	emission_max_len[I]=lenal[0]+lenal[1];


	model[I][I]=Allowed;
	model[I][M]=Allowed;
	model[I][D]=gop;
	model[I][FM]=f_gop;
	model[I][FI]=f_gop;
	model[I][FD]=f_gop;
	
	model[START][I]=model[I][END]=(TG_MODE==0)?0:-model[M][I];
	
	
	
	/*Deletion*/
	emission_delta_i[D]=-1;
	emission_delta_j[D]=+1;
	emission_len_i[D]=3;
	emission_len_j[D]=0;
	
	emission[D ]=gep;
	emission_term[D]=(TG_MODE==2)?-emission[D]:0;
	emission_max_len[D]=lenal[0]+lenal[1];
	
	model[D][D]=Allowed;
	model[D][M]=Allowed;
	model[D][I]=gop;
	model[D][FM]=f_gop;
	model[D][FI]=f_gop;
	model[D][FD]=f_gop;
	
	model[START][D]=model[D][END]=(TG_MODE==0)?0:-model[M][D];

	
	
	/*Frame Match*/
	emission_delta_i[FM]=-1;
	emission_delta_j[FM]= 0;
	emission_len_i[FM]=1;
	emission_len_j[FM]=1;
	emission[FM]=0;
	emission_max_len[FM]=3;
	emission_term[FM]=(F_TG_MODE==2)?-emission[FM]:0;

	model[FM][FM]=Allowed;	 
	model[FM][FI]=Allowed;
 	model[FM][FD]=Allowed; 
	
	model[START][FM]=model[FM][END]=(F_TG_MODE==0)?0:-model[M][FM];

	
	/*Frame Insertion*/
	emission_delta_i[FI]= 0;
	emission_delta_j[FI]=-1;
	emission_len_i[FI]=0;
	emission_len_j[FI]=1;
	emission[FI]=f_gep;	
	emission_max_len[FI]=1;
	emission_term[FI]=(F_TG_MODE==2)?-emission[FI]:0;

	model[FI][I]=gop;
	model[FI][D]=gop;
	model[FI][M]=Allowed;

	model[START][FI]=model[FI][END]=(F_TG_MODE==0)?0:-model[M][FI];

	
	/*Frame Deletion*/
	emission_delta_i[FD]=-1;
	emission_delta_j[FD]=+1;
	emission_len_i[FD]=1;
	emission_len_j[FD]=0;
	emission[FD]=f_gep;
	emission_max_len[FD]=1;
	emission_term[FD]=(F_TG_MODE==2)?-emission[FD]:0;
	
	model[FD][I]=gop;
	model[FD][D]=gop;
	model[FD][M]=Allowed;

	model[START][FD]=model[FD][END]=(F_TG_MODE==0)?0:-model[M][FD];
	


/*Clean The Model*/
	for ( a=0; a< nstate; a++)
	  {
	    emission_len[a]=MAX(emission_len_i[a], emission_len_j[a]);
	    for ( b=0; b< nstate; b++)
	      if( model[a][b]==(UNDEFINED*-1))model[a][b]=UNDEFINED;
	  }

	diag[0]=0;
	for (i=0; i<=lenal[0]; i++)
	  for ( j=0; j<=ndiag; j++)
	    {
	      pos_j=diag[j]-lenal[0]+i;
	      if(pos_j>0 && i>0)continue;
	      else if ( pos_j==0 && i==0)
		{
		  for ( a=0; a< nstate; a++)
		    {
		     hsearch (Mat, 0, a*mI+i*mJ+j, ENTER);
		     hsearch (LMat, 0, a*mI+i*mJ+j, ENTER);
		    }
		}
	      else
		{
		  for ( state=START+1; state<nstate; state++)
		    {
		      pos_j=MAX(pos_j, 0);
		      l=pos_j-i;

		      
		      if (emission_len_i[state] && emission_len_j[state]){continue;}
		      else if( l>0 && emission_len_i[state]){continue;}
		      else if( l<0 && emission_len_j[state]){continue;}
		      if( l==0) continue;
		      else if(FABS(l)>emission_max_len[state])continue;
		      else if (!is_defined_int(4,model[START][state],model[M][state],emission[state],emission_term[state]))continue;
		      else if (l%(emission_len[state]))continue;
		      
		      l=FABS(l);
		      t=model[START][state]+model[M][state];
		      e=emission_term[state]+emission[state];
		      
		      hsearch (Mat,t+e*l/emission_len[state], a*mI+i*mJ+j, ENTER);
		      hsearch (LMat,l, a*mI+i*mJ+j, ENTER);
		    }
		}
	    }

	hstat(Mat);
	fprintf ( stderr, "\n%d %d", lenal[0], j);

/*Forward Pass*/

	

	for (i=1; i<=lenal[0];i++)
	  {
	    hstat(Mat);						
	    for (j=1; j<=ndiag;j++)
	      {
		
		pos_j=diag[j]-lenal[0]+i;
		if (pos_j<=0 || pos_j>l2 )continue;
		last_i=i;
		last_j=j;
		
		
		for (cur_state=START+1; cur_state< nstate; cur_state++)
		  {
		    len_j=len_i=0;
		    prev_i=i;
		    prev_j=j;
		    
		    if ( emission_delta_j[cur_state])
		      {
		      for (; prev_j>0 && prev_j<=ndiag; prev_j+=emission_delta_j[cur_state])
			{
			  prev_i=i+emission_delta_i[cur_state]*FABS((diag[j]-diag[prev_j]));
			  
			  len_i=FABS(i-prev_i);
			  len_j=FABS(((diag[prev_j]-diag[j])+(prev_i-i)));
			  x=MAX(len_i, len_j);
			  

			 
			  if (!x || ((pos_j+emission_delta_j[cur_state]*len_j)<=0))len=UNDEFINED;
			  else if (x%(emission_len[cur_state]))len=UNDEFINED;
			  else
			    {
			      len=x; break;
			    }
			}
		      }
		    else if (!emission_delta_j[cur_state])
		      {
			prev_i=i+emission_delta_i[cur_state]*emission_len_i[cur_state];	
			len_j=len_i=FABS((i-prev_i));
			
			if ((pos_j+emission_delta_j[cur_state]*len_j)<=0)len=UNDEFINED;
			else len=len_j=len_i;
		      }
		    		    
		    if ( cur_state==M)
		      {
			if ( i<3 || pos_j<3)em=UNDEFINED;
			else if ( trs1[i-3]=='x' || trs2[pos_j-3]=='x')em=UNDEFINED;
			else em=(mat[trs1[i-3]-'a'][trs2[pos_j-3]-'a'])*SCORE_K;	
		      }
		    else if ( cur_state==I || cur_state==D)
		      {
			if (i<3 || pos_j<3)em=UNDEFINED;
			
			else if (trs1[i-3]=='x' || trs2[pos_j-3]=='x')em=UNDEFINED;
		      }
		    
		    else 
		      {
			em=emission[cur_state];
		      }

		    for (pc=best_pc=x=UNDEFINED, prev_state=START+1; prev_state< nstate; prev_state++)
		      {

			n=hsearch (CMat,1, prev_state*mI+prev_i*mJ+prev_j, ADD);

			if(prev_i<0 || prev_j<0 ||prev_i>lenal[0] || prev_j>ndiag || len==UNDEFINED)
			  {
			    prev_score=UNDEFINED;
			    prev_len=0;
			  }
			
			else 
			  {
			    prev_score=hsearch ( Mat, 0, prev_state*mI+prev_i*mJ+prev_j, FIND);
			    prev_len=hsearch   ( LMat,0, prev_state*mI+prev_i*mJ+prev_j, FIND);
			  }

			if (n==nstate-2)
			      {
				hsearch (  Mat, 0, prev_state*mI+prev_i*mJ+prev_j,REMOVE);
				hsearch ( LMat, 0, prev_state*mI+prev_i*mJ+prev_j,REMOVE);
				hsearch ( CMat, 0, prev_state*mI+prev_i*mJ+prev_j,REMOVE);
			      }
			/*Set the transition+extension penalties*/

			t=model[prev_state][cur_state];			
			e=em;

			if   (prev_score==UNDEFINED)e=UNDEFINED;			
			else if (len==0 ||len>emission_max_len[cur_state] || e==UNDEFINED)e=UNDEFINED;
			else if (len%(emission_len[cur_state]))e=UNDEFINED;
			else e=e*(len/emission_len[cur_state]);
			
			/*Check that the chain of identical states is not too long*/
			/*Check that ALL the residues have not been matched*/
			if (is_defined_int(3,prev_score,e, t))
			  {
			    if ( prev_state==cur_state && (prev_len+len)>emission_max_len[cur_state])pc=UNDEFINED;
			    else pc=prev_score+t+e;
			  }
			else 
			  pc=UNDEFINED;
			
			/*Identify the best previous score*/
			if (best_pc==UNDEFINED || (pc>best_pc && pc!=UNDEFINED))
			  {
			    prev[cur_state]=prev_state;
			    best_pc=pc;
			  }
		      }
		    hsearch ( Mat,best_pc,cur_state*mI+i*mJ+j,ENTER);

		    if ( best_pc==UNDEFINED)
		      {
			hsearch ( LMat,UNDEFINED,cur_state*mI+i*mJ+j,ENTER);
			continue;
		      }
		    
		    else if ( prev[cur_state]==cur_state)
		      {
			
			n=hsearch ( LMat,0,cur_state*mI+prev_i*mJ+prev_j,FIND);
			hsearch ( LMat,n+len,cur_state*mI+i*mJ+j,ENTER);
		      }
		    else
		      {
			hsearch ( LMat,len,cur_state*mI+i*mJ+j,ENTER);
		      }
		  }
	      }
	  }
	

        i=last_i;
	j=last_j;
	for (pc=best_pc=x=UNDEFINED, state=START+1; state< nstate; state++)
	  {
	    t=model[state][END];
	    e=emission_term[state];
	    l=hsearch ( LMat,0,state*mI+i*mJ+j,FIND);
	    
	    if (!is_defined_int(4,t,e,hsearch(Mat,0,state*mI+i*mJ+j,FIND),l))hsearch(Mat,UNDEFINED,state*mI+i*mJ+j,ENTER);
	    else hsearch(Mat,t+e*(l/emission_len[state]), state*mI+i*mJ+j, ADD);
	    pc=hsearch (Mat, 0,state*mI+i*mJ+j, FIND);
	    if (best_pc==UNDEFINED || (pc>best_pc && pc!=UNDEFINED))
	      {

		best_pc=pc;
	      }
	  }
	 score=best_pc;

	 fprintf ( stderr, "\nSCORE=%d", score);
	 exit (0);

	 return score;
	}

