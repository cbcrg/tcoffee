#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"


int cfasta_cdna_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL)
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
	int l0, l1;
        Alignment *B;
	/********Prepare Penalties******/
		
	maximise=CL->maximise;
	ktup=CL->ktup;

	/********************************/
	if ( !group_list)
	   {
	     
	       group_list=make_group_aa (&n_groups, CL->matrix_for_aa_group);
	   }
	B=dna_aln2_3frame_cdna_aln(A, ns, l_s);
	l0=strlen(B->seq_al[0]);
	l1=strlen(B->seq_al[3]);
	tot_diag=evaluate_diagonals_cdna ( B, ns, l_s, CL, maximise,n_groups,group_list, ktup);
	
	max_n_chosen_diag=100;
	n_chosen_diag=step=10 ;
	n_chosen_diag=MIN(n_chosen_diag, max_n_chosen_diag);
	
	
	diag=extract_N_diag (l0,l1, tot_diag, n_chosen_diag,2);
	score    =make_fasta_cdna_pair_wise ( A,B, ns, l_s, CL, diag);
	
       
	new_score=0;
	vfree ( diag);
	
 	while (new_score!=score && n_chosen_diag< max_n_chosen_diag )
	  {
	    
	    score=new_score;
	    ungap_sub_aln ( A, ns[0], l_s[0]);
	    ungap_sub_aln ( A, ns[1], l_s[1]);
	    
	    
	    n_chosen_diag+=step;
	    n_chosen_diag=MIN(n_chosen_diag, max_n_chosen_diag);
	    
	    diag     =extract_N_diag (l0,l1, tot_diag, n_chosen_diag,3); 
	    new_score=make_fasta_cdna_pair_wise (  A, B,ns, l_s, CL, diag);
	    vfree ( diag);
	  }
	
	score=new_score;
	free_int (tot_diag, -1);
	free_aln(B);
	return score;
    }


int fasta_cdna_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL)
    {
/*TREATMENT OF THE TERMINAL GAP PENALTIES*/
/*TG_MODE=0---> gop and gep*/
/*TG_MODE=1---> ---     gep*/
/*TG_MODE=2---> ---     ---*/


	int maximise;
	int l0, l1;	
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

	B=dna_aln2_3frame_cdna_aln(A, ns, l_s);
	B->nseq=6;
	
	l0=strlen ( B->seq_al[0]);
	l1=strlen ( B->seq_al[3]);
	

	tot_diag=evaluate_diagonals_cdna( B, ns, l_s, CL, maximise,n_groups,group_list, ktup);

	
	diag=extract_N_diag (l0, l1, tot_diag,20,1);
	score=make_fasta_cdna_pair_wise ( A,B, ns, l_s, CL, diag);
	
	free_aln(B);
	free_int (tot_diag, -1);
	vfree (diag);
	return score;
    }

Dp_Model* initialize_dna_dp_model (Constraint_list *CL)
    {
      Dp_Model *M;
      int a, b, c,d;
      int f0, f1;
      int deltaf1, deltaf0,deltatype;
      int type, type1, type0;
            
      M=(Dp_Model*)vcalloc( 1, sizeof (Dp_Model));
     
      for (M->nstate=0,f0=0; f0<3; f0++)
	      for ( f1=0; f1<3; f1++)M->nstate+=3;
      
      M->UM=M->nstate++;
      M->START=M->nstate++;
      M->END  =M->nstate++;

      M->TG_MODE=CL->TG_MODE;
      M->F_TG_MODE=0;
      M->gop=CL->gop*SCORE_K;
      M->gep=CL->gep*SCORE_K;

      

      M->f_gop=CL->f_gop*SCORE_K;
      M->f_gep=CL->f_gep*SCORE_K;
      

      M->bounded_model=declare_int (M->nstate+1, M->nstate+1); 
      M->model=declare_int (M->nstate+1, M->nstate+1); 
      for ( a=0; a<=M->nstate; a++)
	for ( b=0; b<= M->nstate; b++)
	  M->model[a][b]=UNDEFINED;
      M->model_properties=declare_int ( M->nstate, 10);
      
      a=0;     
      M->TYPE=a++;M->F0=a++;M->F1=a++; M->LEN_I=a++; M->LEN_J=a++; M->DELTA_I=a++;M->DELTA_J=a++;M->EMISSION=a++;M->TERM_EMISSION=a++;
      a=M->nstate;
      M->NON_CODING=a++; M->INSERTION=a++; M->DELETION=a++; M->CODING0=a++; M->CODING1=a++;M->CODING2=a++;
      
      
      for ( a=0,f0=0; f0<3; f0++)
	      for ( f1=0; f1<3; f1++, a+=3)
		{
		  M->model_properties[a+0][M->TYPE]=M->CODING0;
		  M->model_properties[a+0][M->F0]=f0;
		  M->model_properties[a+0][M->F1]=f1;
		  M->model_properties[a+0][M->LEN_I]=1;
		  M->model_properties[a+0][M->LEN_J]=1;
		  M->model_properties[a+0][M->DELTA_I]=-1;
		  M->model_properties[a+0][M->DELTA_J]= 0;	
		  M->model_properties[a+0][M->EMISSION]=0;
		  M->model_properties[a+0][M->TERM_EMISSION]=0;
		  
		  M->model_properties[a+1][M->TYPE]=M->DELETION;
		  M->model_properties[a+1][M->F0]=f0;
		  M->model_properties[a+1][M->F1]=f1;
		  M->model_properties[a+1][M->LEN_I]=1;
		  M->model_properties[a+1][M->LEN_J]=0;		 
		  M->model_properties[a+1][M->DELTA_I]=-1;
		  M->model_properties[a+1][M->DELTA_J]=+1;
		  M->model_properties[a+1][M->EMISSION]=M->gep;
		  M->model_properties[a+1][M->TERM_EMISSION]=(M->TG_MODE==2)?0:M->gep;

		  M->model_properties[a+2][M->TYPE]=M->INSERTION;
		  M->model_properties[a+2][M->F0]=f0;
		  M->model_properties[a+2][M->F1]=f1;
		  M->model_properties[a+2][M->LEN_I]=0;
		  M->model_properties[a+2][M->LEN_J]=1;
		  M->model_properties[a+2][M->DELTA_I]= 0;
		  M->model_properties[a+2][M->DELTA_J]=-1;
		  M->model_properties[a+2][M->EMISSION]=M->gep;
		  M->model_properties[a+2][M->TERM_EMISSION]=(M->TG_MODE==2)?0:M->gep;
		}

      /*UM" Unmatched State*/
      M->model_properties[a][M->TYPE]=M->NON_CODING;
      M->model_properties[a][M->F0]=0;
      M->model_properties[a][M->F1]=0;
      M->model_properties[a][M->LEN_I]=1;
      M->model_properties[a][M->LEN_J]=1;
      M->model_properties[a][M->DELTA_I]=-1;
      M->model_properties[a][M->DELTA_J]=0;
      M->model_properties[a][M->EMISSION]=M->f_gep;
      M->model_properties[a][M->TERM_EMISSION]=(M->F_TG_MODE==2)?0:M->f_gep;
      

      M->model_properties[M->START][M->TYPE]=M->NON_CODING;
      M->model_properties[a+1][M->F0]=0;
      M->model_properties[a+1][M->F1]=0;
      M->model_properties[a+1][M->LEN_I]=0;
      M->model_properties[a+1][M->LEN_J]=0;
      M->model_properties[a+1][M->DELTA_I]=0 ;
      M->model_properties[a+1][M->DELTA_J]=0;
      M->model_properties[a+1][M->EMISSION]=0;
      M->model_properties[a+1][M->TERM_EMISSION]=0;
      
      M->model_properties[M->END][M->TYPE]=M->NON_CODING;
      M->model_properties[a+2][M->F0]=0;
      M->model_properties[a+2][M->F1]=0;
      M->model_properties[a+2][M->LEN_I]=0;
      M->model_properties[a+2][M->LEN_J]=0;
      M->model_properties[a+2][M->DELTA_I]=0 ;
      M->model_properties[a+2][M->DELTA_J]=0;
      M->model_properties[a+2][M->EMISSION]=0;
      M->model_properties[a+2][M->TERM_EMISSION]=0;

      /*1: SET THE INDEL PENALTIES*/
     

      for ( a=0; a< M->START; a++)
	      {
		deltaf0=M->model_properties[M->START][M->F0]-M->model_properties[a][M->F0];
		deltaf1=M->model_properties[M->START][M->F1]-M->model_properties[a][M->F1];
		if      ( deltaf0==0 && deltaf1==0)deltatype=0;
		else if ( deltaf0<=0 && deltaf1<=0)deltatype=1;
		else deltatype=-1;
		type=M->model_properties[a][M->TYPE];

		if      ( type==M->NON_CODING)               M->model[a][M->END]=M->model[M->START][a]=(M->F_TG_MODE==0)?M->f_gop:0;
		else if ( type==M->CODING0   && deltatype==0)M->model[a][M->END]=M->model[M->START][a]=ALLOWED;
		else if ( type==M->CODING0   && deltatype==1)M->model[a][M->END]=M->model[M->START][a]=(M->F_TG_MODE==0)?M->f_gop:0;
		else if ( type==M->INSERTION && deltatype==0)M->model[a][M->END]=M->model[M->START][a]=(M->TG_MODE==0)?M->gop:0;
		else if ( type==M->INSERTION && deltatype==1)M->model[a][M->END]=M->model[M->START][a]=(M->TG_MODE==0)?M->gop:0+(M->F_TG_MODE==0)?M->f_gop:0;
		else if ( type==M->DELETION  && deltatype==0)M->model[a][M->END]=M->model[M->START][a]=(M->TG_MODE==0)?M->gop:0;
		else if ( type==M->DELETION  && deltatype==1)M->model[a][M->END]=M->model[M->START][a]=(M->TG_MODE==0)?M->gop:0+(M->F_TG_MODE==0)?M->f_gop:0;
		else  M->model[a][M->END]=M->model[M->START][a]=UNDEFINED;

		/*
		if      (type==M->NON_CODING ||M->model_properties[a][M->F0] ||M->model_properties[a][M->F1]) M->model[M->START][a]=M->model[a][M->END]=(M->F_TG_MODE==0)?M->f_gop:0;
		else if (type==M->INSERTION || type==M->DELETION)M->model[M->START][a]=M->model[a][M->END]=(  M->TG_MODE==0)?M->gop:0;
		else     M->model[M->START][a]=M->model[a][M->END]=ALLOWED;
		*/
		
		for ( b=0; b< M->START; b++)
		  {
		    
		    deltaf0=M->model_properties[a][M->F0]-M->model_properties[b][M->F0];
		    deltaf1=M->model_properties[a][M->F1]-M->model_properties[b][M->F1];
		    type0=M->model_properties[a][M->TYPE];
		    type1=M->model_properties[b][M->TYPE];
		    
		    if      ( deltaf0==0 && deltaf1==0)deltatype=0;
		    else if ( deltaf0<=0 && deltaf1<=0)deltatype=1;
		    else deltatype=-1;
		    
		   
		    
		    if      ( type0==M->NON_CODING && type1==M->NON_CODING                 )M->model[a][b]=UNDEFINED;
		    else if ( type0==M->NON_CODING && type1==M->CODING0                    )M->model[a][b]=ALLOWED ;
		    else if ( type0==M->NON_CODING && type1==M->INSERTION                  )M->model[a][b]=M->gop;
		    else if ( type0==M->NON_CODING && type1==M->DELETION                   )M->model[a][b]=M->gop;
		    
		    else if ( type0==M->CODING0   && type1==M->NON_CODING                 )M->model[a][b]=M->f_gop;
		    else if ( type0==M->CODING0   && type1==M->CODING0   &&  deltatype==0 )M->model[a][b]=ALLOWED;
		    else if ( type0==M->CODING0   && type1==M->CODING0   &&  deltatype==1 )M->model[a][b]=M->f_gop;		    
		    else if ( type0==M->CODING0   && type1==M->INSERTION &&  deltatype==0 )M->model[a][b]=M->gop;
		    else if ( type0==M->CODING0   && type1==M->INSERTION &&  deltatype==1 )M->model[a][b]=M->gop+M->f_gop;		    
		    else if ( type0==M->CODING0   && type1==M->DELETION  &&  deltatype==0 )M->model[a][b]=M->gop;
		    else if ( type0==M->CODING0   && type1==M->DELETION  &&  deltatype==1 )M->model[a][b]=M->gop+M->f_gop;
		    
		    else if ( type0==M->INSERTION && type1==M->NON_CODING                 )M->model[a][b]=M->f_gop;
		    else if ( type0==M->INSERTION && type1==M->CODING0   &&  deltatype==0 )M->model[a][b]=ALLOWED;
		    else if ( type0==M->INSERTION && type1==M->CODING0   &&  deltatype==1 )M->model[a][b]=M->f_gop;		    
		    else if ( type0==M->INSERTION && type1==M->INSERTION &&  deltatype==0 )M->model[a][b]=ALLOWED;
		    else if ( type0==M->INSERTION && type1==M->INSERTION &&  deltatype==1 )M->model[a][b]=M->f_gop;		    
		    else if ( type0==M->INSERTION && type1==M->DELETION  &&  deltatype==0 )M->model[a][b]=M->gop;
		    else if ( type0==M->INSERTION && type1==M->DELETION  &&  deltatype==1 )M->model[a][b]=M->gop+M->f_gop;
		    
		    else if ( type0==M->DELETION  && type1==M->NON_CODING                 )M->model[a][b]=M->f_gop;
		    else if ( type0==M->DELETION  && type1==M->CODING0   &&  deltatype==0 )M->model[a][b]=ALLOWED;
		    else if ( type0==M->DELETION  && type1==M->CODING0   &&  deltatype==1 )M->model[a][b]=M->f_gop;
		    else if ( type0==M->DELETION  && type1==M->INSERTION &&  deltatype==0 )M->model[a][b]=M->gop;
		    else if ( type0==M->DELETION  && type1==M->INSERTION &&  deltatype==1 )M->model[a][b]=M->gop+M->f_gop;
		    else if ( type0==M->DELETION  && type1==M->DELETION  &&  deltatype==0 )M->model[a][b]=ALLOWED;
		    else if ( type0==M->DELETION  && type1==M->DELETION  &&  deltatype==1 )M->model[a][b]=M->f_gop;
		    
		    else {M->model[a][b]=UNDEFINED;}
		    
		  }
	      }
    
	       
      /*2 SET THE FRAMESHIFT PENALTIES
		
      for ( a=0; a< M->START; a++)
	      {
		type=M->model_properties[a][M->TYPE];
		
		for ( b=0; b< M->START; b++)
		  {
		    deltaf0=M->model_properties[a][M->F0]-M->model_properties[b][M->F0];
		    deltaf1=M->model_properties[a][M->F1]-M->model_properties[b][M->F1];
		   
		    
		    
		    
		    if (b==M->UM)                M->model[a][b]+=M->f_gop;
		    else if (a==M->UM)                M->model[a][b]+=ALLOWED;
		    else if (deltaf1==0 && deltaf0==0)M->model[a][b]+=ALLOWED;
		    else if (deltaf1<=0 && deltaf0<=0)M->model[a][b]+=M->f_gop;
		    else M->model[a][b]=UNDEFINED;
		  }
		
	      }         
      M->model[M->UM][M->UM]=UNDEFINED;
      */
     

      for (c=0,a=0, d=0; a< M->START; a++)
	for ( b=0; b<M->START; b++, d++)
	  {
	    if (M->model[a][b]!=UNDEFINED)
	      {
		M->bounded_model[b][1+M->bounded_model[b][0]++]=a;
		c++;
	      }
	  }
      return M;
    }
int make_fasta_cdna_pair_wise (Alignment *B,Alignment *A,int*in_ns, int **l_s,Constraint_list *CL, int *diag)
    {
      int a,c,p,k;
      Dp_Result *DPR;
      static Dp_Model  *M;
      int l0, l1;
      int len_i, len_j;
      int f0=0, f1=0;
      int deltaf0, deltaf1, delta;
      int nr1, nr2;
      int ala, alb, aa0, aa1;
      int type;
      
      char **al;
      int **tl_s;
      int *tns;
      /*DEBUG*/
      int debug_cdna_fasta=0;
      Alignment *DA;
      int score;
      int state,prev_state;
      int t, e;
      int a1, a2;
      
      
      l0=strlen ( B->seq_al[l_s[0][0]]);
      l1=strlen ( B->seq_al[l_s[1][0]]);

      al=declare_char (2, l0+l1+1); 
      B=realloc_aln2 (B,B->nseq,l0+l1+1);


      free_int (B->cdna_cache, -1);
      B->cdna_cache=declare_int(1, l0+l1+1);
      
      if ( !M)M=initialize_dna_dp_model (CL);

     
      M->diag=diag;

      tl_s=(int**)declare_int (2, 2);tns=(int*)vcalloc(2, sizeof(int));tl_s[0][0]=0;tl_s[1][0]=3;tns[0]=tns[1]=1;
      DPR=make_fast_dp_pair_wise (A,tns, tl_s,CL,M);
      vfree(tns);free_int(tl_s, -1);


      
      /*new_trace_back*/
      a=p=0;
      aa0=aa1=ala=alb=0;
      while ( (k=DPR->traceback[a++])!=M->START);
      while ( (k=DPR->traceback[a++])!=M->END)
	{
	  
	  f0=M->model_properties[k][M->F0];
	  f1=M->model_properties[k][M->F1];

	  len_i=M->model_properties[k][M->LEN_I];
	  len_j=M->model_properties[k][M->LEN_J];
	  
	  type=M->model_properties[k][M->TYPE];
	  
	  

	  if (type==M->CODING0)
	    {
	      deltaf0=(aa0*3+f0)-ala;
	      deltaf1=(aa1*3+f1)-alb;

	      delta=MAX(deltaf0, deltaf1);
	      
	      for (nr1=0, nr2=0,c=0; c<delta; c++, nr1++, nr2++,p++)		  
		      {
			if (nr1<deltaf0 && ala<l0)al[0][p]=B->seq_al[l_s[0][0]][ala++];
			else al[0][p]='-';
			
			if (nr2<deltaf1 && alb<l1)al[1][p]=B->seq_al[l_s[1][0]][alb++];
			else al[1][p]='-'; 
			
			B->cdna_cache[0][p]=M->NON_CODING;	
			if ( is_gap(al[1][p]) && is_gap(al[0][p]))p--;
			else if ( debug_cdna_fasta)fprintf (stderr, "\nUM: %c %c",  al[0][p], al[1][p]);
		      } 
	      for ( c=0; c< 3; c++, p++)
		{
		  if ( c==0)B->cdna_cache[0][p]=M->CODING0;
		  else if ( c==1)B->cdna_cache[0][p]=M->CODING1;
		  else if ( c==2)B->cdna_cache[0][p]=M->CODING2;
		  if (ala<l0)al[0][p]=B->seq_al[l_s[0][0]][ala++];
		  else al[0][p]='-';

		  if (alb<l1)al[1][p]=B->seq_al[l_s[1][0]][alb++];
		  else al[1][p]='-';
			
		  if ( is_gap(al[1][p]) && is_gap(al[0][p]))p--;
		  else if ( debug_cdna_fasta)fprintf (stderr, "\n%d: %c %c",k,  al[0][p], al[1][p]);
		}
	    }

	  aa0+=len_i;
	  aa1+=len_j;
	}
      
      deltaf0=(aa0*3+f0)-ala;
      deltaf1=(aa1*3+f1)-alb;
      delta=MAX(deltaf0, deltaf1);
      for (nr1=0, nr2=0,c=0; c<delta; c++, nr1++, nr2++,p++)		  
	{
	  if (nr1<deltaf0 && ala<l0)al[0][p]=B->seq_al[l_s[0][0]][ala++];
	  else al[0][p]='-';
	  
	  if (nr2<deltaf1 && alb<l1)al[1][p]=B->seq_al[l_s[1][0]][alb++];
	  else al[1][p]='-'; 
	  
	  B->cdna_cache[0][p]=M->NON_CODING;	
	  if ( is_gap(al[1][p]) && is_gap(al[0][p]))p--;
	  else if ( debug_cdna_fasta)fprintf (stderr, "\nUM: %c %c",  al[0][p], al[1][p]);
	}
      

      /*End New traceback*/
      



      al[0][p]='\0';
      al[1][p]='\0';


      sprintf( B->seq_al[l_s[0][0]], "%s", al[0]);
      sprintf( B->seq_al[l_s[1][0]], "%s", al[1]);
      B->len_aln=strlen (al[0]);
      B->nseq=2;
     
      
     
      
      if ( debug_cdna_fasta)
	  {
	    fprintf ( stderr, "\nA-A=%d, %d", CL->M['a'-'A']['a'-'A'], CL->M['a'-'A']['a'-'A'] *SCORE_K);
	    for ( a=1; a<diag[0]; a++)
	      {
		fprintf ( stderr, "\nchosen diag: %d", diag[a]);
	      }
	    
	    fprintf ( stderr, "\n  GOP=%d   GEP=%d   TG_MODE=%d", M->gop, M->gep, M->TG_MODE);
	    fprintf ( stderr, "\nF_GOP=%d F_GEP=%d F_TG_MODE=%d", M->gop, M->gep, M->F_TG_MODE);
	    
	    DA=copy_aln (B, NULL);
	    DA=realloc_aln2 (DA,6,(DA->len_aln+1));
	

	    for ( a=0; a<B->len_aln; a++)
	      {

		fprintf ( stderr, "\n%d", DA->cdna_cache[0][a]);
		if (DA->cdna_cache[0][a]>=M->CODING0)DA->seq_al[DA->nseq][a]=DA->cdna_cache[0][a]-M->nstate+'0';
		else DA->seq_al[DA->nseq][a]=DA->cdna_cache[0][a]-M->nstate+'0';

		if (DA->cdna_cache[0][a]==M->CODING0)
		  {
		    DA->seq_al[DA->nseq+1][a]=translate_dna_codon (DA->seq_al[0]+a,'*');
		    DA->seq_al[DA->nseq+2][a]=translate_dna_codon (DA->seq_al[1]+a,'*');
		  }
		else
		  {
		    DA->seq_al[DA->nseq+1][a]='-'; 
		    DA->seq_al[DA->nseq+2][a]='-'; 
		  }
		
	      }
	    DA->nseq+=3;
	    print_aln (DA);
	    
	    free_aln(DA);		      
	    score=0;
	    
	    
	    for (prev_state=M->START,a=0; a< DA->len_aln;)
	      {
		state=DA->cdna_cache[0][a];
		t=M->model[prev_state][state];
		if ( DA->cdna_cache[0][a]==M->CODING0)
		  {
		    a1=translate_dna_codon (A->seq_al[0]+a,'x');
		    a2=translate_dna_codon (A->seq_al[1]+a,'x');
		    
		    if ( a1!='x' && a2!='x')
		      {
			e=CL->M[a1-'A'][a2-'A']*SCORE_K;
		      }
		  }
		else if ( DA->cdna_cache[0][a]>M->CODING0);
		else
		  {
		    e=M->model_properties[B->cdna_cache[0][a]][M->EMISSION];
		  }
		if ( e==UNDEFINED || t==UNDEFINED) fprintf ( stderr, "\nPROBLEM %d\n", a);
		
		fprintf ( stderr, "\n[%c..%c: %d(e)+%d(t)=%d]", A->seq_al[0][a], A->seq_al[1][a], e,t,e+t);
		score+=e+t;
		prev_state=state;
		
		if (B->cdna_cache[0][a]==M->NON_CODING)a++;
		else a+=3;
		
	      }
	    
	  }
      
      for ( a=0; a<B->len_aln; a++)
	{
	  
	  if ( B->cdna_cache[0][a]<M->CODING0)B->cdna_cache[0][a]=0;
	  else B->cdna_cache[0][a]=1;
	}
      
      free_char ( al, -1);
      return DPR->score;
      
    }
		


Dp_Result * make_fast_dp_pair_wise (Alignment *A,int*ns, int **l_s, Constraint_list *CL,Dp_Model *M)
	{
	  
	  /*SIZE VARIABLES*/ 
	  
	  int ndiag;
	  int l0, l1, len_al,len_diag;
	  static int max_len_al, max_len_diag;
	  static int mI, mJ;
	 
	  
	  /*EVALUATION*/
	  int **mat;
	  int a1, a2;
	  
	  /*DP VARIABLES*/
	  static int *Mat, *LMat, *trace;
	  int a, i, j,l;
	  int state, cur_state, prev_state;
	  int pos_i,  pos_j;
	  int last_i=0, last_j=0;
	  int prev_i, prev_j;
	  int len_i, len_j, len;
	  int t, e, em;
	  
	  int prev_score; 
	  int pc, best_pc;
	  
	  int *prev;
	  int model_index;
	  /*TRACEBACK*/
	  Dp_Result *DPR;
	  int k=0, next_k;
	  int new_i, new_j;
	  
	  
	  ndiag=M->diag[0];

	  l0=strlen (A->seq_al[l_s[0][0]]);
	  l1=strlen (A->seq_al[l_s[1][0]]);
	  len_al =l0+l1+1;	
	  len_diag=ndiag+4;
	  
	  if ( (len_al>max_len_al || len_diag>max_len_diag))
	    {
	      
	      vfree (Mat);
	      vfree (LMat);
	      vfree(trace);	    
	      max_len_diag=max_len_al=0;	   
	    }
	  
	  if (max_len_al==0)
	    {
	      max_len_al=len_al;
	      max_len_diag=len_diag;
	      mI=max_len_al*max_len_diag;
	      mJ=max_len_diag;
	      
	      
	      Mat  =(int*)vcalloc ( M->nstate*max_len_al*max_len_diag, sizeof (int));
	      LMat =(int*)vcalloc ( M->nstate*max_len_al*max_len_diag, sizeof (int));
	      trace=(int*)vcalloc ( M->nstate*max_len_al*max_len_diag, sizeof (int));
	      
	    }
	  
	  prev=(int*)vcalloc ( M->nstate, sizeof (int));
	  DPR=( Dp_Result*) vcalloc ( 1, sizeof ( Dp_Result));
	  DPR->traceback=(int*)vcalloc (max_len_al, sizeof (int));
	  
/*PREPARE THE EVALUATION*/      
	  if (ns[0]+ns[1]>2)
	    {
	      fprintf ( stderr, "\nERROR: function make_fasta_cdna_pair_wise can only handle two sequences at a time [FATAL:%s]",PROGRAM);
	      crash ("");
	    }
	  mat=CL->M; 					  		

/*INITIALIZATION OF THE DP MATRICES*/

	for (i=0; i<=l0;i++)
	  {						
	    for (j=0; j<=ndiag;j++)
	      {
		for ( state=0; state<M->nstate; state++)
		  {
		    Mat   [state*mI+i*mJ+j]=UNDEFINED;
		    LMat  [state*mI+i*mJ+j]=UNDEFINED;
		    trace [state*mI+i*mJ+j]=M->START;
		  }
	      }
	  }	

	M->diag[0]=0;

	for (i=0; i<=l0; i++)
	  for ( j=0; j<=ndiag; j++)
	    {
	      pos_j=M->diag[j]-l0+i;
	      pos_i=i;
	      if (!(pos_j==0 || pos_i==0))continue;
	      if ( pos_j<0 || pos_i<0)continue;
	      if ( pos_i==0 && pos_j==0)
		  {
		  for ( a=0; a< M->nstate; a++)
		    {
		     Mat  [a*mI+i*mJ+j]=0;
		     LMat [a*mI+i*mJ+j]=0;
		     trace[a*mI+i*mJ+j]=M->START;
		    }
		}
	      else
		{	
		  l=MAX(pos_i,pos_j);
		  for ( state=0; state<M->START; state++)
		    {		     
		      if (pos_j==0 && M->model_properties[state][M->LEN_J])continue;
		      if (pos_i==0 && M->model_properties[state][M->LEN_I])continue;
		     
		     
		     t=M->model[M->START][state];
		     e=M->model_properties[state][M->TERM_EMISSION];
		     Mat   [state*mI+i*mJ+j]=t+e*l;
		     LMat  [state*mI+i*mJ+j]=l;
		     trace [state*mI+i*mJ+j]=M->START;
		    }
		}
	    }

/*DYNAMIC PROGRAMMING: Forward Pass*/

	

	for (i=1; i<=l0;i++)
	  {						
	    for (j=1; j<=ndiag;j++)
	      {
		pos_j=M->diag[j]-l0+i;
		pos_i=i;
		
		if (pos_j<=0 || pos_j>l1 )continue;
		last_i=i;
		last_j=j;
		
		for (cur_state=0; cur_state<M->START; cur_state++)
		  {
		    if (M->model_properties[cur_state][M->DELTA_J])
		      {
			prev_j=j+M->model_properties[cur_state][M->DELTA_J];
			prev_i=i+M->model_properties[cur_state][M->DELTA_I]*FABS((M->diag[j]-M->diag[prev_j]));			
		      }
		    else
		      {
			prev_j=j;
			prev_i=i+M->model_properties[cur_state][M->DELTA_I];
		      }
		    len_i=FABS((i-prev_i));
		    len_j=FABS((M->diag[prev_j]-M->diag[j]));
		    len=MAX(len_i, len_j);
		    a1=A->seq_al[M->model_properties[cur_state][M->F0]  ][pos_i-1];
		    a2=A->seq_al[M->model_properties[cur_state][M->F1]+3][pos_j-1];
		
		    if (M->model_properties[cur_state][M->TYPE]==M->CODING0)
		      {
			if ( a1=='o' || a2=='o')em=-(mat['w'-'A']['w'-'A'])*SCORE_K;
			else if (a1=='x' || a2=='x')em=UNDEFINED;
			else if ( a1==0 || a2==0)exit (0);
			else 
			  {
			    em=(mat[a1-'A'][a2-'A'])*SCORE_K;
			  }
		      }
		    else
		      {
			em=M->model_properties[cur_state][M->EMISSION];
		      }
		    
		    
		   
		    for (pc=best_pc=UNDEFINED, model_index=1; model_index<=M->bounded_model[cur_state][0]; model_index++)
		      {
			prev_state=M->bounded_model[cur_state][model_index];
			
			if(prev_i<0 || prev_j<0 ||prev_i>l0 || prev_j>ndiag || len==UNDEFINED)prev_score=UNDEFINED;
			else prev_score=Mat[prev_state*mI+prev_i*mJ+prev_j];
			t=M->model[prev_state][cur_state];			
			e=em;
		
			if   (prev_score==UNDEFINED || len==UNDEFINED)e=UNDEFINED;			
			else if (len==0|| e==UNDEFINED)e=UNDEFINED;
			else e=e*len;
			
			if (is_defined_int(3,prev_score,e, t))
			  {
			    pc=prev_score+t+e;
			  }
			else  pc=UNDEFINED;
			
			/*Identify the best previous score*/
			if (best_pc==UNDEFINED || (pc>best_pc && pc!=UNDEFINED))
			  {
			    prev[cur_state]=prev_state;
			    best_pc=pc;
			   
			  }
		      }
		    
		    Mat[cur_state*mI+i*mJ+j]=best_pc;
		   


		    if ( Mat[cur_state*mI+i*mJ+j]==UNDEFINED)
		      {
			LMat[cur_state*mI+i*mJ+j]=UNDEFINED;
			trace[cur_state*mI+i*mJ+j]=UNDEFINED;
			continue;
		      }
		    
		    else if ( prev[cur_state]==cur_state)
		      {
			LMat [cur_state*mI+i*mJ+j]=	LMat [cur_state*mI+prev_i*mJ+prev_j]+len;
			trace[cur_state*mI+i*mJ+j]=     trace[cur_state*mI+prev_i*mJ+prev_j];
		      }
		    else
		      {
			LMat[cur_state*mI+i*mJ+j]=len;
			trace[cur_state*mI+i*mJ+j]=prev[cur_state];
		      }
		  }
	      }
	  }
	
	
        i=last_i;
	j=last_j;
	for (pc=best_pc=UNDEFINED, state=0; state<M->START; state++)
	  {
	    t=M->model[state][M->END];
	    e=M->model_properties[state][M->TERM_EMISSION];
	    l=LMat[state*mI+i*mJ+j];
	    
	   
	    if (!is_defined_int(4,t,e,Mat[state*mI+i*mJ+j],l))Mat[state*mI+i*mJ+j]=UNDEFINED;
	    else Mat[state*mI+i*mJ+j]+=t+e*(l);
	    pc=Mat[state*mI+i*mJ+j];
	    
	   
	    if (best_pc==UNDEFINED || (pc>best_pc && pc!=UNDEFINED))
	      {
		k=state;
		best_pc=pc;
	      }
	  }
	 DPR->score=best_pc;
	
/*TRACEBACK*/ 


	e=0;
	len=0;    
	
	
	while (k!=M->START)
	  {
	    next_k=trace[k*mI+i*mJ+j];
	    new_i=i;
	    new_j=j;
	    l=LMat[k*mI+i*mJ+j];
	   	   
	    for (a=0; a< l; a++)
	      {
		DPR->traceback[len++]=k;
	      }
	   new_i+=M->model_properties[k][M->DELTA_I]*l;
	   
	   
	   if ( M->model_properties[k][M->DELTA_J])
	     {
	       while ( next_k!=M->START && FABS((M->diag[j]-M->diag[new_j]))!=l)new_j+=M->model_properties[k][M->DELTA_J];
	     }

	   i=new_i;
	   j=new_j;
	   k=next_k;
	  }
	DPR->len=len;
	DPR->traceback[DPR->len++]=M->START;
	invert_list_int  (DPR->traceback,DPR->len);
	DPR->traceback[DPR->len]=M->END;
	
	vfree (prev);

	return DPR;
	

	}



int ** evaluate_diagonals_cdna ( Alignment *B, int *ns, int **l_s, Constraint_list *CL,int maximise,int n_groups, char **group_list, int ktup)
        {
	  int f1, f2, c;
	  int **diag;
	  char *s1, *s2;
	  int p1, p2;
	  int **tot_diag;
	  int n_tot_diag;
	  int l0, l1;
	  
	  
	 
	  
	  
	 
	  if ( ns[0]!=1 || ns[1]!=1)
	    {
	      fprintf ( stderr, "\nERROR 2 SEQUENCES ONLY [FATAL:%s", PROGRAM);
	      crash ("");
	    }
	  
	  
	  
	  
	
	l0=strlen ( B->seq_al[0]);
	l1=strlen ( B->seq_al[3]);
	n_tot_diag=(l0+l1-1);

	tot_diag=declare_int ( n_tot_diag+1, 2);
	for ( c=0; c<= n_tot_diag; c++)tot_diag[c][0]=c; 
	  
	for (f1=0; f1< 3; f1++)
	    {
	      for ( f2=0; f2< 3; f2++)
		{
		  s1=B->seq_al[f1];
		  s2=B->seq_al[3+f2];
		  
		  
		  p1=strlen (s1);
		  p2=strlen (s2);


		  diag=evaluate_diagonals_for_two_sequences( s1, s2, maximise,NULL,ktup);
		  for (c=1; c<=(p1+p2-1); c++)
		    { 
		      tot_diag[diag[c][0]][1]+=diag[c][1]*diag[c][1];
		    }
		  free_int (diag, -1);
		  
		}
	    }
	
	

	  sort_int (tot_diag+1, 2, 1,0, n_tot_diag-1);	  
	  
	  return tot_diag;
	  
	}
	  



