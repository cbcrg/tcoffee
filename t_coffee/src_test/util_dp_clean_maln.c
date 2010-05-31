#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"

Alignment *clean_maln ( Alignment *A, Alignment *I, int T, int n_it)
  {
    Alignment *C=NULL;
    int a, b;
    int in_una, in_aln, in_gap, gap, una, aln;
    int Sstart,Rstate, Sstate;
    int n_segment=0;
    int **segment_list;
    
        
    add_warning ( stderr, "\nWARNING: -clean_aln is not supported anymore [PROGRAM:%s]\n", PROGRAM);
    return A;

    
   
    /*Initialization*/
    a=0;
    in_una=a++;in_gap=a++;in_aln=a++; aln=a++;gap=a++;una=a++;       
    segment_list=declare_int ( A->len_aln*A->nseq, 3);
    
    
    /*1: Identify the segments*/
    C=copy_aln(A, C);    
    for ( a=0; a< A->nseq; a++)
      {
	Sstate=in_aln;
	for ( b=0; b<A->len_aln; b++)
	  {
	    if (is_gap(A->seq_al[a][b]))Rstate=gap;
	    else if ( I->seq_al[a][b]<=T){Rstate=una;}
	    else if ( I->seq_al[a][b]==NO_COLOR_RESIDUE)Rstate=una;
	    else Rstate=aln;
	    
	    if (Rstate==una)C->seq_al[a][b]='-';

	    if (Sstate==in_aln)
	      {
		if ( Rstate==gap)
		  {Sstate=in_gap;
		  Sstart=b;
		  }
		else if ( Rstate==una)
		  {
		    Sstate=in_una;
		    Sstart=b;
		  }
		else if ( Rstate==aln)
		  Sstate=in_aln;
	      }
	    else if ( Sstate==in_gap)
	      {
		if ( Rstate==gap);
		else if ( Rstate==una)Sstate=in_una;
		else if ( Rstate==aln)Sstate=in_aln;
	      }
	    else if ( Sstate==in_una)
	      {
		if ( Rstate==gap);
		else if ( Rstate==una);
		else if ( Rstate==aln)
		  {
		    segment_list[n_segment][0]=a;
		    segment_list[n_segment][1]=Sstart;
		    segment_list[n_segment][2]=b-Sstart;
		    Sstate=in_aln;
		    n_segment++;
		  }
	      }
	  }
	if (Sstate==in_una)
	  {
	   segment_list[n_segment][0]=a;
	   segment_list[n_segment][1]=Sstart;
	   segment_list[n_segment][2]=b-Sstart;
	   Sstate=in_aln;
	   n_segment++;
	  }
      }

    /*2 Realign the segments*/
    
    for ( b=0; b< n_it; b++)
      {
	for ( a=0; a< n_segment; a++)
	  {
	    HERE ("1");
	    A=realign_segment ( segment_list[a][0], segment_list[a][1], segment_list[a][2], A, C);
	   
	  }
      }
    free_aln (C);
    free_int ( segment_list, -1);
    make_fast_generic_dp_pair_wise (NULL, NULL, NULL, NULL);
    
    return A;
  }
Alignment *realign_segment (int seq, int start, int len,Alignment *A, Alignment *C)
  {
    Alignment *S1=NULL, *S2=NULL, *S3=NULL;
    int  *ns, **ls;
    int a,b;
    static Constraint_list *CL;
    
    
    /*1 Prepare the Constraint list*/
    if ( !CL)
      {
	CL=vcalloc ( 1, sizeof (Constraint_list));
	CL->extend_jit=0;
	CL->pw_parameters_set=1;
	CL->M=read_matrice ("blosum62mt");
	CL->gop=-20;
	CL->gep=-1;	
	CL->evaluate_residue_pair=evaluate_matrix_score;
	sprintf ( CL->dp_mode, "myers_miller_pair_wise");
      }
   
    S1=copy_aln(A,S1);
    S1=extract_aln (S1,0,start);
    S2=copy_aln(A,S2);
    S2=extract_aln (S2, start, start+len);
    S3=copy_aln(A,S3);
    S3=extract_aln (S3, start+len,A->len_aln);

    
    /*for (a=0; a<S2->nseq; a++){S2->order[a][1]=0;S2->order[a][0]=a;}*/
    
    
    ungap ( S2->seq_al[seq]);   
    CL->S=A->S;/*aln2seq(S2);*/    
    /*3 Prepare Sequence Presentation*/
    ns=vcalloc (2, sizeof (int));
    ls=declare_int (2,S2->nseq);
    
    ns[0]=A->nseq-1;
    for ( a=0,b=0; a< S2->nseq; a++)if (a!=seq)ls[0][b++]=a;
    ns[1]=1;
    ls[1][0]=seq;
    
    pair_wise (S2, ns, ls, CL);
    
    A=realloc_aln (A, strlen (S1->seq_al[0])+ strlen (S2->seq_al[0])+ strlen (S3->seq_al[0])+1);
    for ( a=0; a< A->nseq; a++)
      {
	sprintf ( A->seq_al[a], "%s%s%s", S1->seq_al[a], S2->seq_al[a], S3->seq_al[a]);
      }

    free_aln (S1);
    free_aln (S2);
    free_aln (S3);
    vfree(ns);free_int(ls, -1);
    
    return A;
  }
Alignment *realign_segment_old (int seq, int start, int len,Alignment *A, Alignment *C)
  {
    Alignment *S=NULL;
    int  *ns, **ls;
    char *sub_seq;
    static Dp_Model *M=NULL;
    static Constraint_list *CL=NULL;
    Dp_Result *R=NULL;
    int a,b, c;
 
    
    /*1 Prepare the Constraint list*/
    if ( !CL)
      {
	CL=vcalloc ( 1, sizeof (Constraint_list));
	CL->extend_jit=0;
	CL->pw_parameters_set=1;
	CL->M=read_matrice ("blosum62mt");
	CL->gop=-20;
	CL->gep=-1;	
	CL->evaluate_residue_pair=evaluate_matrix_score;
	
      }
    S=copy_aln(C,S);
    S=extract_aln (S, start, start+len);
    S->len_aln=strlen(S->seq_al[0]);       
    sub_seq=extract_char  (A->seq_al[seq], start, len);   
    
    ungap(sub_seq);
   
    sprintf ( S->seq_al[seq],"%s", sub_seq);
    CL->S=aln2seq(S);    
    
   

    /*2 Prepare the Model*/
    M=initialize_seg2prf_model((start==0)?2:0,(start+len==A->len_aln)?2:0,CL);    
    M->diag=vcalloc ( 2*len+1, sizeof (int));
    M->diag[0]=len+strlen (sub_seq)-1;
    for ( a=1; a<=M->diag[0]; a++)M->diag[a]=a;
    
    /*3 Prepare Sequence Presentation*/
    ns=vcalloc (2, sizeof (int));
    ls=declare_int (2,A->nseq);
    
    ns[0]=A->nseq-1;
    for ( a=0,b=0; a< A->nseq; a++)if (a!=seq)ls[0][b++]=a;
    ns[1]=1;
    ls[1][0]=seq;

    if ( strlen (sub_seq)!=len)
      {

	
	R=make_fast_generic_dp_pair_wise(S, ns, ls, M);

	for (c=0, b=1,a=start; a< start+len; b++,a++)
	  {
	    if (R->traceback[b]==0)
	      {
		A->seq_al[seq][a]=sub_seq[c];
		C->seq_al[seq][a]=sub_seq[c];
		c++;
	      }
	    else 
	      {
		A->seq_al[seq][a]='-';
		C->seq_al[seq][a]='-';
	      }
	  }
      }

    free_dp_model  (M);
    free_aln (S);
    free_dp_result (R);
    vfree(sub_seq);
    vfree(ns);
    free_int (ls, -1);
    free_sequence (CL->S, (CL->S)->nseq);
    
    
    return A;
  }

Dp_Model * initialize_seg2prf_model(int left_tg_mode, int right_tg_mode, Constraint_list *CL)
  {
    
    Dp_Model *M;
    int a, b, c,d;
    
    M=vcalloc ( 1, sizeof (Dp_Model));
    M->nstate=2;
    M->START=M->nstate++;
    M->END  =M->nstate++;
    
    M->TG_MODE=1;
    M->F_TG_MODE=0;
    M->gop=CL->gop*SCORE_K;
    M->gep=CL->gep*SCORE_K;
    
    M->bounded_model=declare_int (M->nstate+1, M->nstate+1); 
    M->model=declare_int (M->nstate+1, M->nstate+1); 
    for ( a=0; a<=M->nstate; a++)
      for ( b=0; b<= M->nstate; b++)
	M->model[a][b]=UNDEFINED;
    
    a=0;     
    M->TYPE=a++;M->LEN_I=a++; M->LEN_J=a++; M->DELTA_I=a++;M->DELTA_J=a++; M->CODING0=a++;M->DELETION=a++;
    M->model_properties=declare_int ( M->nstate, 10);
    
    a=0;
    M->EMISSION=a++;M->TERM_EMISSION=a++;M->START_EMISSION=a++;
    M->model_emission_function=vcalloc(M->nstate, sizeof (int (**)(Alignment*, int **, int, int*, int, int **, int, int*, int, struct Constraint_list *)));
    for ( a=0; a< M->nstate; a++)
       M->model_emission_function[a]=vcalloc(3, sizeof (int (*)(Alignment*, int **, int, int*, int, int **, int, int*, int, struct Constraint_list *)));
    
    
    /*Substitution*/
    M->model_properties[0][M->TYPE]=M->CODING0;
    M->model_properties[0][M->LEN_I]=1;
    M->model_properties[0][M->LEN_J]=1;
    M->model_properties[0][M->DELTA_I]=-1;
    M->model_properties[0][M->DELTA_J]= 0;	

    M->model_emission_function[0][M->EMISSION]      =cw_profile_get_dp_cost;
    M->model_emission_function[0][M->START_EMISSION]=get_start_gep_cost;
    M->model_emission_function[0][M->TERM_EMISSION] =get_start_gep_cost;
   
    /*Deletions*/       
    M->model_properties[1][M->TYPE]=M->DELETION;
    M->model_properties[1][M->LEN_I]=1;
    M->model_properties[1][M->LEN_J]=0;
    M->model_properties[1][M->DELTA_I]=-1;
    M->model_properties[1][M->DELTA_J]=+1;
    M->model_emission_function[1][M->EMISSION]=get_gep_cost;

    if (left_tg_mode ==2)
      M->model_emission_function[1][M->START_EMISSION]=get_start_gep_cost;
    else M->model_emission_function[1][M->START_EMISSION]=get_gep_cost;
    
    if (right_tg_mode ==2)
      M->model_emission_function[1][M->TERM_EMISSION]=get_term_gep_cost;
    else M->model_emission_function[1][M->TERM_EMISSION]=get_gep_cost;
              
    /*Transitions*/
    M->model[0][M->END]=M->model[M->START][0]=ALLOWED;
    M->model[0][1]=M->gop;
    M->model[0][0]=ALLOWED;
    
    M->model[1][M->END]=  (right_tg_mode==0)?0:-M->gop;
    M->model[M->START][1]=( left_tg_mode==0)?M->gop:0;    
    M->model[1][1]=ALLOWED;
    M->model[1][0]=ALLOWED;
    
    
    
    
    /*Prune the model*/

    for (c=0,a=0, d=0; a< M->START; a++)
      for ( b=0; b<M->START; b++, d++)
	{
	  if (M->model[a][b]!=UNDEFINED)
	    {
	      M->bounded_model[b][1+M->bounded_model[b][0]++]=a;
	      c++;
	    }
	}
    M->CL=CL;
   
    return M;
  }

int get_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
{
  return CL->gep*SCORE_K;
}

int get_start_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
{
  return 0;
}
int get_term_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
{
  return CL->gep*SCORE_K*-1;
}




