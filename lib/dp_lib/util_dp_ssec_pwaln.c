#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"

int ssec_pwaln_maln (Alignment *A, int *ns, int **ls, Constraint_list *CL)
  {

    static Dp_Model *M=NULL;
    Dp_Result *R=NULL;
    int a, ndiag;
    int Sa,Sb,St, Da, Db, Dt, Ia, Ib, It;
    int ala, alb, s,b;

    a=0;    
    Sa=a++;Da=a++;Ia=a++;
    Sb=a++;Db=a++;Ib=a++;
    St=a++;Dt=a++;It=a++;
    
    if ( strm (CL->matrices_list[0], "analyse"))
      {
      for ( a=0; a< CL->n_matrices; a++)
	{

	  rescale_two_mat(CL->matrices_list[1],CL->matrices_list[2],1000, 100, AA_ALPHABET);	
	  exit (0);
	}
      }

    

    /*2 Prepare the Model*/
    M=initialize_sseq_model(2,2,CL);    
    ndiag=strlen (A->seq_al[0])+strlen (A->seq_al[1])-1;
    M->diag=(int*)vcalloc (ndiag+1, sizeof (int));
    M->diag[0]=ndiag-1;
    for ( a=1; a<=M->diag[0]; a++)M->diag[a]=a;
    
    /*3 Prepare Sequence Presentation*/
    R=make_fast_generic_dp_pair_wise(A, ns, ls, M);
    
  
    ala=alb=0;

    A=realloc_aln2(A,A->nseq, R->len+1);
    for (b=1; b<R->len;b++)
      {
	if (R->traceback[b]==Sa ||  R->traceback[b]==Sb ||R->traceback[b]==St )
	  {
	    for (s=0; s<ns[0]; s++)
	      A->seq_al[ls[0][s]][b-1]=(CL->S)->seq[A->order[ls[0][s]][0]][ala];
	    ala++;
	    for (s=0; s<ns[1]; s++)
	      A->seq_al[ls[1][s]][b-1]=(CL->S)->seq[A->order[ls[1][s]][0]][alb];
	    alb++;
	  }
	else if ( R->traceback[b]==Da ||  R->traceback[b]==Db ||R->traceback[b]==Dt )
	  {
	    for (s=0; s<ns[0]; s++)
	      A->seq_al[ls[0][s]][b-1]=(CL->S)->seq[A->order[ls[0][s]][0]][ala];
	    ala++;
	    for (s=0; s<ns[1]; s++)
	      A->seq_al[ls[1][s]][b-1]='-';	   
	  }
	else if ( R->traceback[b]==Ia ||  R->traceback[b]==Ib ||R->traceback[b]==It )
	  {
	    for (s=0; s<ns[0]; s++)
	      A->seq_al[ls[0][s]][b-1]='-';
	    
	    for (s=0; s<ns[1]; s++)
	      A->seq_al[ls[1][s]][b-1]=(CL->S)->seq[A->order[ls[1][s]][0]][alb];
	    alb++;
	  }
      }
    for (s=0; s<ns[0]; s++)
      A->seq_al[ls[0][s]][b-1]='\0';
    for (s=0; s<ns[1]; s++)
      A->seq_al[ls[1][s]][b-1]='\0';
    
    A->len_aln=strlen (A->seq_al[ls[0][0]]);
    R->Dp_model=M;
    A->Dp_result=R;
    return A->score;
  }

Dp_Model * initialize_sseq_model(int left_tg_mode, int right_tg_mode, Constraint_list *CL)
  {
    
    Dp_Model *M;
    int a, b, c,d;
    int Sa,Sb,St, Da, Db, Dt, Ia, Ib, It;
    int tgop=CL->gep*3;
    

    
    
    M=(Dp_Model*)vcalloc ( 1, sizeof (Dp_Model));
    
    M->nstate=9;
    M->START=M->nstate++;
    M->END  =M->nstate++;
    
    M->model_comments=declare_char (M->nstate+1, 100);
    M->bounded_model=declare_int (M->nstate+1, M->nstate+1); 
    M->model=declare_int (M->nstate+1, M->nstate+1); 
    for ( a=0; a<=M->nstate; a++)
      for ( b=0; b<= M->nstate; b++)
	M->model[a][b]=UNDEFINED;
    
    
    M->model_properties=declare_int ( M->nstate, 10); 
    
    a=0;     
    M->TYPE=a++;M->LEN_I=a++; M->LEN_J=a++; M->DELTA_I=a++;M->DELTA_J=a++;M->EMISSION=a++;M->TERM_EMISSION=a++;M->START_EMISSION=a++;
    M->CODING0=a++;M->DELETION=a++;
    M->model_properties=declare_int ( M->nstate, 10); 

    a=0;
    M->EMISSION=a++;M->TERM_EMISSION=a++;M->START_EMISSION=a++;
    M->model_emission_function=(int (***)(Alignment*, int **, int, int*, int, int **, int, int*, int, struct Constraint_list *))vcalloc(M->nstate, sizeof (int (**)(Alignment*, int **, int, int*, int, int **, int, int*, int, struct Constraint_list *)));
    for ( a=0; a< M->nstate; a++)
       M->model_emission_function[a]=(int (**)(Alignment*, int **, int, int*, int, int **, int, int*, int, struct Constraint_list *))vcalloc(3, sizeof (int (*)(Alignment*, int **, int, int*, int, int **, int, int*, int, struct Constraint_list *)));
    

        
    a=0;
    
    Sa=a++;Da=a++;Ia=a++;
    Sb=a++;Db=a++;Ib=a++;
    St=a++;Dt=a++;It=a++;
   

    sprintf ( M->model_comments[M->START], "START");
    sprintf ( M->model_comments[M->END], "END");
	      
    /*ALPHA*/
    /*Substitution in Alpha*/
    if (CL->matrices_list[0][0])sprintf ( M->model_comments[Sa], "Substitution %s", CL->matrices_list[0]);
    M->model_properties[Sa][M->TYPE]=Sa;
    M->model_properties[Sa][M->LEN_I]=1;
    M->model_properties[Sa][M->LEN_J]=1;
    M->model_properties[Sa][M->DELTA_I]=-1;
    M->model_properties[Sa][M->DELTA_J]= 0;	

    M->model_emission_function[Sa][M->EMISSION]      =get_alpha_sub_cost;
    M->model_emission_function[Sa][M->START_EMISSION]=get_ssec_no_cost;
    M->model_emission_function[Sa][M->TERM_EMISSION] =get_ssec_no_cost;
   
    /*Deletions*/       
    if (CL->matrices_list[0][0])sprintf ( M->model_comments[Da], "Deletion %s", CL->matrices_list[0]);
    M->model_properties[Da][M->TYPE]=Da;
    M->model_properties[Da][M->LEN_I]=1;
    M->model_properties[Da][M->LEN_J]=0;
    M->model_properties[Da][M->DELTA_I]=-1;
    M->model_properties[Da][M->DELTA_J]=+1;

    
    M->model_emission_function[Da][M->EMISSION]      =get_alpha_gep_cost;
    M->model_emission_function[Da][M->START_EMISSION]=get_alpha_start_gep_cost;
    M->model_emission_function[Da][M->TERM_EMISSION] =get_alpha_term_gep_cost;

        
    /*Insertion*/
    if (CL->matrices_list[0][0])sprintf ( M->model_comments[Ia], "Insertion %s", CL->matrices_list[0]);
    M->model_properties[Ia][M->TYPE]=Ia;
    M->model_properties[Ia][M->LEN_I]=0;
    M->model_properties[Ia][M->LEN_J]=1;
    M->model_properties[Ia][M->DELTA_I]=0;
    M->model_properties[Ia][M->DELTA_J]=-1;
    
    M->model_emission_function[Ia][M->EMISSION]      =get_alpha_gep_cost;
    M->model_emission_function[Ia][M->START_EMISSION]=get_alpha_start_gep_cost;
    M->model_emission_function[Ia][M->TERM_EMISSION] =get_alpha_term_gep_cost;
    
/*BETA*/
    /*Substitution in Beta*/
    if (CL->matrices_list[1][0])sprintf ( M->model_comments[Sb], "Substitution %s", CL->matrices_list[1]);
    M->model_properties[Sb][M->TYPE]=Sb;
    M->model_properties[Sb][M->LEN_I]=1;
    M->model_properties[Sb][M->LEN_J]=1;
    M->model_properties[Sb][M->DELTA_I]=-1;
    M->model_properties[Sb][M->DELTA_J]= 0;	
    
    M->model_emission_function[Sb][M->EMISSION]      =get_beta_sub_cost;
    M->model_emission_function[Sb][M->START_EMISSION]=get_ssec_no_cost;
    M->model_emission_function[Sb][M->TERM_EMISSION] =get_ssec_no_cost;
    
   
    /*Deletions*/       
    if (CL->matrices_list[1][0])sprintf ( M->model_comments[Db], "Deletion %s", CL->matrices_list[1]);
    M->model_properties[Db][M->TYPE]=Db;
    M->model_properties[Db][M->LEN_I]=1;
    M->model_properties[Db][M->LEN_J]=0;
    M->model_properties[Db][M->DELTA_I]=-1;
    M->model_properties[Db][M->DELTA_J]=+1;
    
    M->model_emission_function[Db][M->EMISSION]      =get_beta_gep_cost;
    M->model_emission_function[Db][M->START_EMISSION]=get_beta_start_gep_cost;
    M->model_emission_function[Db][M->TERM_EMISSION] =get_beta_term_gep_cost;
    
    
    /*Insertion*/
    
    if (CL->matrices_list[1][0])sprintf ( M->model_comments[Ib], "Insertion %s", CL->matrices_list[1]);
    M->model_properties[Ib][M->TYPE]=Ib;
    M->model_properties[Ib][M->LEN_I]=0;
    M->model_properties[Ib][M->LEN_J]=1;
    M->model_properties[Ib][M->DELTA_I]=0;
    M->model_properties[Ib][M->DELTA_J]=-1;

    
    
    M->model_emission_function[Ib][M->EMISSION]      =get_beta_gep_cost;
    M->model_emission_function[Ib][M->START_EMISSION]=get_beta_start_gep_cost;
    M->model_emission_function[Ib][M->TERM_EMISSION] =get_beta_term_gep_cost;
    
 /*TURNS*/
    /*Substitution in Turn*/
    if (CL->matrices_list[2][0])sprintf ( M->model_comments[St], "Substitution %s", CL->matrices_list[2]);
    M->model_properties[St][M->TYPE]=St;
    M->model_properties[St][M->LEN_I]=1;
    M->model_properties[St][M->LEN_J]=1;
    M->model_properties[St][M->DELTA_I]=-1;
    M->model_properties[St][M->DELTA_J]= 0;
	
    M->model_emission_function[St][M->EMISSION]      =get_turn_sub_cost;
    M->model_emission_function[St][M->START_EMISSION]=get_ssec_no_cost;
    M->model_emission_function[St][M->TERM_EMISSION] =get_ssec_no_cost;
    
   
    /*Deletions*/       
    if (CL->matrices_list[2][0])sprintf ( M->model_comments[Dt], "Deletion %s", CL->matrices_list[2]);
    M->model_properties[Dt][M->TYPE]=Dt;
    M->model_properties[Dt][M->LEN_I]=1;
    M->model_properties[Dt][M->LEN_J]=0;
    M->model_properties[Dt][M->DELTA_I]=-1;
    M->model_properties[Dt][M->DELTA_J]=+1;
    
    M->model_emission_function[Dt][M->EMISSION]      =get_turn_gep_cost;
    M->model_emission_function[Dt][M->START_EMISSION]=get_turn_start_gep_cost;
    M->model_emission_function[Dt][M->TERM_EMISSION] =get_turn_term_gep_cost;
    /*Insertion*/
    if (CL->matrices_list[2][0])sprintf ( M->model_comments[It], "Insertion %s", CL->matrices_list[2]);
    M->model_properties[It][M->TYPE]=It;
    M->model_properties[It][M->LEN_I]=0;
    M->model_properties[It][M->LEN_J]=1;
    M->model_properties[It][M->DELTA_I]=0;
    M->model_properties[It][M->DELTA_J]=-1;

    M->model_emission_function[It][M->EMISSION]      =get_turn_gep_cost;
    M->model_emission_function[It][M->START_EMISSION]=get_turn_start_gep_cost;
    M->model_emission_function[It][M->TERM_EMISSION] =get_turn_term_gep_cost;


/*Transitions*/

    M->model[M->START][Sa]=ALLOWED;
    M->model[M->START][Sb]=ALLOWED;
    M->model[M->START][St]=ALLOWED;   
    M->model[M->START][Db]=M->model[M->START][Ib]=(CL->TG_MODE==0)?CL->gop*SCORE_K:0;
    M->model[M->START][Da]=M->model[M->START][Ia]=(CL->TG_MODE==0)?CL->gop*SCORE_K:0;
    M->model[M->START][Dt]=M->model[M->START][It]=(CL->TG_MODE==0)?CL->gop*SCORE_K:0;
    
    
    M->model[Sa][M->END]=ALLOWED;
    M->model[Sb][M->END]=ALLOWED;
    M->model[St][M->END]=ALLOWED;
    M->model[Ia][M->END]=M->model[Da][M->END]=(CL->TG_MODE==0)?0:CL->gop*SCORE_K*(-1);
    M->model[Ib][M->END]=M->model[Db][M->END]=(CL->TG_MODE==0)?0:CL->gop*SCORE_K*(-1);
    M->model[It][M->END]=M->model[Dt][M->END]=(CL->TG_MODE==0)?0:CL->gop*SCORE_K*(-1);
    
    for ( a=0; a< M->nstate; a++)M->model[a][a]=ALLOWED;
    
    M->model[Sa][Ia]=M->model[Sa][Da]=CL->gop*SCORE_K;
    M->model[Sa][Ib]=M->model[Sa][Db]=CL->gop*SCORE_K+tgop*SCORE_K;
    M->model[Sa][It]=M->model[Sa][Dt]=CL->gop*SCORE_K+tgop*SCORE_K;
    M->model[Sa][Sb]=M->model[Sa][St]=tgop*SCORE_K;

    M->model[Sb][Ib]=M->model[Sb][Db]=CL->gop*SCORE_K;
    M->model[Sb][Ia]=M->model[Sb][Da]=CL->gop*SCORE_K+tgop*SCORE_K;
    M->model[Sb][It]=M->model[Sb][Dt]=CL->gop*SCORE_K+tgop*SCORE_K;
    M->model[Sb][Sa]=M->model[Sb][St]=tgop*SCORE_K;

    M->model[St][It]=M->model[St][Dt]=CL->gop*SCORE_K;
    M->model[St][Ia]=M->model[St][Da]=CL->gop*SCORE_K+tgop*SCORE_K;
    M->model[St][Ib]=M->model[St][Db]=CL->gop*SCORE_K+tgop*SCORE_K;
    M->model[St][Sa]=M->model[St][Sb]=tgop*SCORE_K;
    
    M->model[Ia][Sa]=M->model[Da][Sa]=ALLOWED;
    M->model[Ia][Sb]=M->model[Da][Sb]=tgop*SCORE_K;
    M->model[Ia][St]=M->model[Da][St]=tgop*SCORE_K;

    M->model[Ib][Sa]=M->model[Db][Sa]=tgop*SCORE_K;
    M->model[Ib][Sb]=M->model[Db][Sb]=ALLOWED;
    M->model[Ib][St]=M->model[Db][St]=tgop*SCORE_K;

    M->model[It][Sa]=M->model[Dt][Sa]=tgop*SCORE_K;
    M->model[It][Sb]=M->model[Dt][Sb]=tgop*SCORE_K;
    M->model[It][St]=M->model[Dt][St]=ALLOWED;
    

        
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





int get_alpha_sub_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
{
  static int **mat;
  int s1, r1, s2, r2;
  float score;
  

  if (!mat && CL->matrices_list[0][0])mat=read_matrice (CL->matrices_list[0]);
  else if ( !CL->matrices_list[0][0])return UNDEFINED;
  

  
  
  s1=A->order[list1[0]][0];
  r1=pos1[list1[0]][col1];
  s2=A->order[list2[0]][0];
  r2=pos1[list2[0]][col2];
  
  if ( r1<0 || r2<0)return 0;

  score=mat[(CL->S)->seq[s1][r1-1]-'A'][(CL->S)->seq[s2][r2-1]-'A']*SCORE_K;
  return (int)score;
  
}
int get_beta_sub_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
{
  static int **mat;
  int s1, r1, s2, r2;
  float score;
  
   if (!mat && CL->matrices_list[1][0])mat=read_matrice (CL->matrices_list[1]);
   else if ( !CL->matrices_list[1][0])return UNDEFINED;
 

  
  s1=A->order[list1[0]][0];
  r1=pos1[list1[0]][col1];
  s2=A->order[list2[0]][0];
  r2=pos1[list2[0]][col2];
  if ( r1<0 || r2<0)return 0;
  
  score=mat[(CL->S)->seq[s1][r1-1]-'A'][(CL->S)->seq[s2][r2-1]-'A']*SCORE_K;
  return (int)score;
  
}
int get_turn_sub_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
{
  static int **mat;
  int s1, r1, s2, r2;
  float score;

  
  
  if (!mat && CL->matrices_list[2][0])mat=read_matrice (CL->matrices_list[2]);
  else if ( !CL->matrices_list[2][0])return UNDEFINED;
  
  
  s1=A->order[list1[0]][0];
  r1=pos1[list1[0]][col1];
  s2=A->order[list2[0]][0];
  r2=pos1[list2[0]][col2];

  
  if ( r1<0 || r2<0)return 0;
  score=mat[(CL->S)->seq[s1][r1-1]-'A'][(CL->S)->seq[s2][r2-1]-'A']*SCORE_K;
  return (int)score;
  
}

int get_turn_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
{
  return (CL->gep) *SCORE_K;
}
int get_turn_start_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
{
  return ((CL->TG_MODE)==2)?0:get_turn_gep_cost(A,pos1, ns1, list1, col1, pos2, ns2, list2, col2, CL);
}
int get_turn_term_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
{
  return ((CL->TG_MODE)==2)?-get_turn_gep_cost(A,pos1, ns1, list1, col1, pos2, ns2, list2, col2, CL):0;
}


int get_alpha_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
{
  return (CL->gep)*SCORE_K;
}
int get_alpha_start_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
{
  return ((CL->TG_MODE)==2)?0:get_alpha_gep_cost(A,pos1, ns1, list1, col1, pos2, ns2, list2, col2, CL);
}
int get_alpha_term_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
{
  return ((CL->TG_MODE)==2)?-get_alpha_gep_cost(A,pos1, ns1, list1, col1, pos2, ns2, list2, col2, CL):0;
}


int get_beta_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
{
  return (CL->gep)*SCORE_K;
}
int get_beta_start_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
{
  return ((CL->TG_MODE)==2)?0:get_beta_gep_cost(A,pos1, ns1, list1, col1, pos2, ns2, list2, col2, CL);
}
int get_beta_term_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
{
  return ((CL->TG_MODE)==2)?-get_beta_gep_cost(A,pos1, ns1, list1, col1, pos2, ns2, list2, col2, CL):0;
}


int get_ssec_no_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
{
  return 0;
}




