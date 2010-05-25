#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "fsa_dp_lib_header.h"



int fsa_nw_pair_wise (Alignment *A, int *ns, int **ls, Constraint_list *CL)
  {
    int a, b, c,d;
    static Fsa_dp_model **M=NULL;
    Fsa_dp_result *R=NULL;
   
   
    /*1 Prepare the Model*/
    if ( !M)
      {
	M=vcalloc (16, sizeof (Fsa_dp_model));
	for ( a=1; a>=0; a--)for(b=1; b>=0; b--)for(c=1; c>=0; c--)for (d=1; d>=0; d--)
	  {M[btoi(4,a,b,c,d)]=initialize_fsa_nw_model(a,b,c,d,CL->TG_MODE,CL);fprintf (stderr, "\n%d: %d %d %d %d",btoi(4,a,b,c,d), a, b, c, d);} 
      }
    


    /*2 Chose The Diagonals*/    
    for ( b=0; b< 16; b++)
      {
	M[b]->diag=vcalloc (strlen (A->seq_al[0])+strlen (A->seq_al[1])+1, sizeof (int));
	M[b]->diag[0]=strlen (A->seq_al[0])+strlen (A->seq_al[1])-1;
	for ( a=1; a<=M[b]->diag[0]; a++)M[b]->diag[a]=a;
      }

    /*3 Get the Trace back*/
    R=fsa_dp_pair_wise(A, ns, ls, M);
    
    /*4 Compute the alignment*/

    A=traceback2aln ( A, ns, ls, M[0], R);
      
    return A->score;
  }

int quadratic_fsa_nw_pair_wise (Alignment *A, int *ns, int **ls, Constraint_list *CL)
  {

    static Fsa_dp_model *M=NULL;
    Fsa_dp_result *R=NULL;
    int a;

    /*1 Prepare the Model*/
    M=initialize_fsa_nw_model(1,1,1,1,CL->TG_MODE,CL);

    /*2 Chose The Diagonals*/    
    M->diag=vcalloc (strlen (A->seq_al[0])+strlen (A->seq_al[1])+1, sizeof (int));
    M->diag[0]=strlen (A->seq_al[0])+strlen (A->seq_al[1])-1;
    for ( a=1; a<=M->diag[0]; a++)M->diag[a]=a;
    

    /*3 Get the Trace back*/
    R=quadratic_fsa_dp_pair_wise(A, ns, ls, M);
    
    /*4 Compute the alignment*/

    A=traceback2aln ( A, ns, ls, M, R);
      
    return A->score;
  }
Fsa_dp_model * initialize_fsa_nw_model(int tl_s1, int tl_s2,int tr_s1, int tr_s2, int tg_mode, Constraint_list *CL)
  {
    
    Fsa_dp_model *M;
    int a, b, c,d;
    int S,D,I;


    
    
    M=vcalloc ( 1, sizeof (Fsa_dp_model));
    
    M->nstate=3;
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
    M->TYPE=a++;M->LEN_I=a++; M->LEN_J=a++; M->DELTA_I=a++;M->DELTA_J=a++;
    a=0;
    M->EMISSION=a++;M->TERM_EMISSION=a++;M->START_EMISSION=a++;
    M->model_properties=declare_int ( M->nstate, 10); 
    
    M->model_emission_function=vcalloc(M->nstate, sizeof (int (**)(Alignment*, int **, int, int*, int, int **, int, int*, int, struct Constraint_list *)));
    for ( a=0; a< M->nstate; a++)
       M->model_emission_function[a]=vcalloc(3, sizeof (int (*)(Alignment*, int **, int, int*, int, int **, int, int*, int, struct Constraint_list *)));
    
    
        
    a=0;    
    S=a++;I=a++;D=a++;

    sprintf ( M->model_comments[M->START], "START");
    sprintf ( M->model_comments[M->END], "END");
	      

    
    sprintf ( M->model_comments[S], "Substitution");
    M->model_properties[S][M->TYPE]=S;
    M->model_properties[S][M->LEN_I]=1;
    M->model_properties[S][M->LEN_J]=1;
    M->model_properties[S][M->DELTA_I]=-1;
    M->model_properties[S][M->DELTA_J]= 0;	

    M->model_emission_function[S][M->EMISSION]      =fsa_nw_get_sub_cost;
    M->model_emission_function[S][M->START_EMISSION]=fsa_nw_get_no_cost;
    M->model_emission_function[S][M->TERM_EMISSION] =fsa_nw_get_no_cost;
   
    /*Deletions*/       
    sprintf ( M->model_comments[D], "Deletion");
    M->model_properties[D][M->TYPE]=D;
    M->model_properties[D][M->LEN_I]=1;
    M->model_properties[D][M->LEN_J]=0;
    M->model_properties[D][M->DELTA_I]=-1;
    M->model_properties[D][M->DELTA_J]=+1;

    
    M->model_emission_function[D][M->EMISSION]      =fsa_nw_get_gep_cost;
    M->model_emission_function[D][M->START_EMISSION]=(tl_s1)?fsa_nw_get_start_gep_cost:fsa_nw_get_gep_cost;
    M->model_emission_function[D][M->TERM_EMISSION] =(tr_s1)?fsa_nw_get_term_gep_cost:fsa_nw_get_gep_cost;
    
    M->model_emission_function[D][M->START_EMISSION]=fsa_nw_get_start_gep_cost;
    M->model_emission_function[D][M->TERM_EMISSION] =fsa_nw_get_term_gep_cost;
    /*Insertion*/

    sprintf ( M->model_comments[I], "Insertion");
    M->model_properties[I][M->TYPE]=I;
    M->model_properties[I][M->LEN_I]=0;
    M->model_properties[I][M->LEN_J]=1;
    M->model_properties[I][M->DELTA_I]=0;
    M->model_properties[I][M->DELTA_J]=-1;
    
    M->model_emission_function[I][M->EMISSION]      =fsa_nw_get_gep_cost;
    M->model_emission_function[I][M->START_EMISSION]=(tl_s2)?fsa_nw_get_start_gep_cost:fsa_nw_get_gep_cost;
    M->model_emission_function[I][M->TERM_EMISSION] =(tr_s2)?fsa_nw_get_term_gep_cost:fsa_nw_get_gep_cost;
    
    M->model_emission_function[I][M->START_EMISSION]=fsa_nw_get_start_gep_cost;
    M->model_emission_function[I][M->TERM_EMISSION] =fsa_nw_get_term_gep_cost;
/*Transitions*/

    M->model[M->START][S]=ALLOWED;
    M->model[M->START][D]=(CL->TG_MODE==0)?CL->gop*SCORE_K:0;
    M->model[M->START][I]=(CL->TG_MODE==0)?CL->gop*SCORE_K:0;
    
    M->model[S][M->END]=ALLOWED;
    
    M->model[D][M->END]=(CL->TG_MODE==0)?0:CL->gop*SCORE_K*(-1);
    M->model[I][M->END]=(CL->TG_MODE==0)?0:CL->gop*SCORE_K*(-1);
    
    for ( a=0; a< M->nstate; a++)M->model[a][a]=ALLOWED;    
    M->model[S][I]=M->model[S][D]=CL->gop*SCORE_K;
    M->model[I][S]=M->model[D][S]=ALLOWED;

            
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





int fsa_nw_get_sub_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
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

  score=mat[(CL->S)->seq[s1][r1-1]-'a'][(CL->S)->seq[s2][r2-1]-'a']*SCORE_K;

  
  return (int)score;
  
}
int fsa_nw_get_no_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
{
  return 0;
}
int fsa_nw_get_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
{
  return (CL->gep) *SCORE_K;
}
int fsa_nw_get_start_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
{
  
  return ((CL->TG_MODE)==2)?0:fsa_nw_get_gep_cost(A,pos1, ns1, list1, col1, pos2, ns2, list2, col2, CL);
}
int fsa_nw_get_term_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
{
  return ((CL->TG_MODE)==2)?-fsa_nw_get_gep_cost(A,pos1, ns1, list1, col1, pos2, ns2, list2, col2, CL):0;
}
