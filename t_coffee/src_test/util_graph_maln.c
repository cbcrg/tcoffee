#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"
int check_link (CL_node ***G, int S1, int r1, int s2, int r2);
void light_nodes (CL_node *A, int va, CL_node*B, int vb, CL_node*C,int vc, char *s);
void print_graph (CL_node *G, Sequence *S);
Sequence *Seq;
CL_node *Start;
static int cycle;



Alignment * add_constraint2aln ( Alignment *A, int s1, int r1, int s2, int r2)
{
  /*Note sCL_node ***G;1 and r1 must be numbered from 0 to n-1*/
  CL_node ***G;
  
  
  G=aln2graph(A);
  G=add_constraint2graph_aln (G,s1,r1, s2, r2);
  A=graph2aln (A,G[0][0],aln2seq(A));
  vfree_graph (G[0][0]);
  return A;
}
Alignment * graph_aln (Alignment *A, Constraint_list *iCL, Sequence *S)
{
  CL_node ***G;
  int a,start;
  CLIST_TYPE *entry=NULL;
  Constraint_list *CL;
  Seq=S;
  
  CL=duplicate_constraint_list (iCL);
  for ( a=0; a<CL->ne; a++)
    {
      CL->L[a*CL->entry_len+WE]=CL->evaluate_residue_pair ( iCL, CL->L[a*CL->entry_len+SEQ1],CL->L[a+CL->entry_len+R1] ,CL->L[a*CL->entry_len+SEQ2],CL->L[a*CL->entry_len+R2]);
     }
   CL=sort_constraint_list_on_n_fields(CL, 0, CL->ne,WE,1);
   G=aln2graph (A);

   
  Start=G[0][0];
  start=0;
  for (a=start; a<CL->ne; a++)
    {
      
      cycle++;
      entry=extract_entry(entry,a, CL);
      fprintf ( stderr, "\r\tCompletion: [%5.2f%%]",(float)(a*100)/CL->ne);
      G=add_constraint2graph_aln (G, entry[SEQ1], entry[R1]-1, entry[SEQ2],entry[R2]-1);
    }
   start=0;
  for (a=start; a<CL->ne; a++)
    {
      
      cycle++;
      entry=extract_entry(entry,a, CL);
      fprintf ( stderr, "\r\tCompletion: [%5.2f%%]",(float)(a*100)/CL->ne);
      G=add_constraint2graph_aln (G, entry[SEQ1], entry[R1]-1, entry[SEQ2],entry[R2]-1);
    }
  
  A=graph2aln (A,G[0][0], S);
  
  return A;

}

void print_graph (CL_node *G, Sequence *S)
{
  Alignment *A=NULL;

  if (S==NULL) S=Seq;
  A=seq2aln (Seq, A, 1);
  A=graph2aln (A,G, Seq);
  print_aln (A);
}

Alignment* graph2aln (Alignment *A, CL_node *G, Sequence *S)
{
  int s, l, a;
  CL_node *Gi;

  /*Rewind G*/
  while ( G->p)G=G->p;
  while ( G->l)G=G->l;
  
  l=s=0;
  Gi=G;
  while (Gi){Gi=Gi->r;l++;}
  Gi=G;
  while (Gi){Gi=Gi->c;s++;}

  A=realloc_alignment (A, l+1);
  
  
  l=0;
  while (G)
    {
      Gi=G;
      s=0;
      while ( G!=NULL)
	{

	  if ( G->res==-1)A->seq_al[s][l]='-';
	  else if ( G->res==-2)A->seq_al[s][l]='*';
	  else if ( G->res==-3)A->seq_al[s][l]='#';
	  else if (G->res>=0)A->seq_al[s][l]=S->seq[G->seq][G->res];
	  
	  G=G->c;
	  s++;
	}
      G=Gi->r;
      l++;
    }
 
  

  for ( a=0;a<s; a++)
    A->seq_al[a][l]='\0';
  A->len_aln=strlen (A->seq_al[0]);
  A->nseq=s;

  
  return A;
}
CL_node ***aln2graph (Alignment *A)
{
  int a=0, b;
  static CL_node ***galn;
  CL_node *N, *iN, *pN;
  int res;
  
  if ( !galn)
    {
      galn=calloc ( A->nseq, sizeof (CL_node**));
      for ( a=0; a< A->nseq; a++)
	{
	  galn[a]=calloc (A->len_aln, sizeof (CL_node*));
	}
    }
  pN=iN=NULL;
  N=declare_cl_nodes(-1, a);
  
  
  for ( a=0; a<A->nseq; a++)
    {
      iN=N;
      for (res=A->order[a][1], b=0; b<A->len_aln; b++)
	{	
	if (b<A->len_aln-1)
	  {
	    N->r=declare_cl_nodes(-1, a);
	    (N->r)->l=N;
	  }
	if ( pN)
	  {
	    N->p=pN;
	    pN->c=N;
	  }
   
	N->seq=A->order[a][0];
	N->res=(is_gap(A->seq_al[a][b]))?-1:res++;
	
	if (N->res!=-1)galn[a][res-1]=N;
	if ( pN)
	  {
	    N->p=pN;
	    pN->c=N;
	    pN=pN->r;
	  }
	N=N->r;
      }
      if ( a<A->nseq-1)iN->c=declare_cl_nodes(-1, a);
      pN=iN;
      N=iN->c;
    }

  return galn;
}
CL_node ***add_constraint2graph_aln (CL_node ***G, int s1, int r1, int s2, int r2)
{
  CL_node *S, *E, *B;
  int d;
  
  
  S=G[s1][r1];
  E=G[s2][r2];
  
    
  d=get_node_distance (S,E);  
  if (d<0){B=S;S=E;E=B;}
  d=(d<0)?-d:d;
 
  
  insert_gap_columns (E,d);
  shift_segment( S,d+1,d);
      
  return G;

}


CL_node *shift_segment ( CL_node *S, int segL, int shiftL)
{
  int a;
  CL_node *E, *G;
  
  
  if ( !shiftL)return S;
  
  /*find segment coordinates*/
  for (E=S, a=1; a< segL; a++)E=E->r;
    
  /*Shift the gaps*/
  
  G=swap_gap_in_graph (S, E);
  for (a=1; a< shiftL; a++)swap_gap_in_graph (S, E);
  
  while (G!=E)G=remove_graph_gap_column (G);
  remove_graph_gap_column (E);
  
  return G;
}
 
int  is_graph_gap_column(CL_node *S)
{
  while (S->p)S=S->p;
  
  while (S)
    {
      if (S->res>=0)return 0;
      S=S->c;
    }
  return 1;
}
CL_node * remove_graph_gap_column (CL_node *S)
{
  CL_node *R,*L, *P, *RV;
  
  RV=S->r;
  while (S->p)
    {
      
      S=S->p;
    }
  
  if ( !is_graph_gap_column (S))return RV;  
  
  
  
  while (S)
    {
      
      R=S->r;
      L=S->l;
      P=S->p;
      
      if (L)L->r=S->r;
      if (R)R->l=S->l;
      
      P=S;
      S=S->c;
      vfree_cl_node (P);
    }
  return RV;
}

CL_node * swap_gap_in_graph ( CL_node*S, CL_node *E)
{
  /*Moves gap AFTER End to BEFORE Start
    SxxxE-
    -xxxxx
   straightens the links in between
  */
  CL_node *G, *N, *iE, *iS, *SP, *SC, *SL;

  
  /*Preserve the E/S values*/
  iE=E;
  iS=S;
 
  
  /*prepare the parent/child links first*/
  
  SP=S->p;
  SC=S->c;
  SL=S->l;

  while ( S!=E->r)
    {
      N=S->r;
      
      S->p=N->p;
      if (N->p)(S->p)->c=S;
      
      S->c=N->c;
      if (N->c)(S->c)->p=S;

      S=S->r;
    }
  
  E=iE;
  S=iS;
  
  /*Remove the gap*/
  G=E->r;
  if ( G->res>=0)fprintf ( stderr, "\nERROR: NOT a GAP");
  
  E->r=G->r;
  if (E->r)(E->r)->l=E;
  
  /*insert the gap*/
  
  G->r=S;
  S->l=G;
  
  G->l=SL;
  if (SL)SL->r=G;
  

  G->p=SP;
  if (SP)SP->c=G;

  G->c=SC;
  if (SC)SC->p=G;

  return G;
 
}

CL_node * declare_cl_nodes ( int len, int seq)
{
  static CL_node **N;
  CL_node *IN;
  static int Nlen;
  int a;

  if (len==-1)
    {
      IN=calloc ( 1, sizeof (CL_node)); 
      IN->res=-1;
      return IN;
    }
      
      

  if ( len>Nlen)
    {
      free (N);
      N=calloc (len, sizeof (CL_node*));
    }

  if ( len==0)return NULL;

  for (a=0; a<len; a++)N[a]=calloc ( 1, sizeof (CL_node)); 
  for (a=0; a<len; a++)
    {
      (N[a])->res=-1;
      (N[a])->seq=seq;
      if (a!=0)(N[a])->l=N[a-1];
      if (a!=len-1)(N[a])->r=N[a+1];
    }
  
  (N[0])->l=N[len-1];
  (N[len-1])->r=N[0];
  
  return N[0];
}
      
CL_node *insert_gap_columns (CL_node *S, int d)
{
  CL_node *Gs,*Ge, *pGs, *Gi, *Si;
  int a;

  if ( d==0)return S;

  pGs=Gi=NULL;
  Si=S;
  while (S->p!=NULL)S=S->p;
  
  while (S!=NULL)
    {
      Gs=declare_cl_nodes(d, S->seq);
      Ge=Gs->l;
      
      Ge->r=S->r;
      if (Ge->r)(Ge->r)->l=Ge;
      
      Gs->l=S;
      S->r=Gs;

      if (pGs)
	{
	  Gi=Gs;
	  for (a=0; a< d; a++)
	    {
	      Gs->p=pGs;
	      pGs->c=Gs;
	      Gs=Gs->r;
	      pGs=pGs->r;
	    }
	  pGs=Gi;
	}
      else
	{
	  pGs=Gs;
	}
      S=S->c;
    }
  return Si;
}

int get_node_distance ( CL_node *S, CL_node *E)
{
  int distance=0;
  CL_node *iS,*B;
  int swap=1;
   
  /*project the two points onto one sequence*/
  if (S->seq>E->seq){B=S;S=E;E=B;swap*=-1;}
  while (S->seq!=E->seq)S=S->c;
    
  /*Walk from E to S */
  iS=S;
  while ( iS->res<0 && iS->r!=NULL){iS=iS->r;}  
  if (iS->res<0 || iS->res>E->res){B=S; S=E; E=B;swap*=-1;}
      
  while ( S!=E)
    {
      S=S->r;
      distance+=swap;
    }
  return distance;
}
  

   
	 
  




int check_graph ( CL_node *S, char *string)
{
  CL_node *iS;
  static int n;
  int lr;
  
  if ( S==NULL)S=Start;
  fprintf ( stderr, "\n\tGRAPH Check %s #%d\n",string, ++n);
  while ( S->p!=NULL)S=S->p;
  while ( S->l!=NULL)S=S->l;
  while ( S)
    {
      iS=S;
      lr=-1;
      while (iS)
	{
	  if (iS->l && (iS->l)->seq!=iS->seq){fprintf ( stderr, "\n\t\tSEq pb");myexit(EXIT_FAILURE);}
	  if (iS->free==1){fprintf ( stderr, "\n\t\tFree Node read");myexit(EXIT_FAILURE);}
	  if (iS->res>0)
	    {
	      if (lr!=-1 && iS->res-lr!=1){fprintf ( stderr, "\n\t\tERROR: lost residues");myexit (EXIT_FAILURE);}
	      lr=iS->res;
	    }
	  if ( iS->r && (iS->r)->l!=iS){fprintf ( stderr, "\n\t\tERROR: left != right: [%d %d][%d %d]", iS->seq, iS->res, (iS->l)->seq, (iS->r)->res);myexit (EXIT_FAILURE);} 
	  if ( iS->p && (iS->p)->c!=iS){fprintf ( stderr, "\n\t\tERROR: parent != child: [%d %d][%d %d]", iS->seq, iS->res, (iS->p)->seq, (iS->p)->res);myexit (EXIT_FAILURE);} 
	  if ( iS->c && (iS->c)->p!=iS){fprintf ( stderr, "\n\t\tERROR: parent != child: [%d %d][%d %d]", iS->seq, iS->res, (iS->c)->seq, (iS->c)->res);myexit (EXIT_FAILURE);} 
	  iS=iS->r;
	}
      S=S->c;      
    }
  return 1;
}

CL_node * vfree_graph ( CL_node *S)
{
  CL_node *Si;
  
  while ( S->p!=NULL)S=S->p;
  while ( S->l!=NULL)S=S->l;

  while ( S)
    {
      Si=S->c;
      while ( S)
	{
	  
	  S=S->r;
	  if (S)vfree_cl_node (S->l);
	}
      S=Si;
    }
  return S;

}
CL_node *vfree_cl_node ( CL_node *N)
{
  if ( N->free==1)crash("freeing free block");
  N->free=1;
  free (N);
  return N;
}


void light_nodes (CL_node *A, int va, CL_node*B, int vb, CL_node*C,int vc, char *string )
{
  int ta=0, tb=0, tc=0;
 
  fprintf ( stderr, "\nCycle %d\n LIGHT NODE: %s", cycle,string);
  if ( A){ta=A->res; A->res=va;fprintf ( stderr, "\nA: seq %d res %d", A->seq, A->res);}
  if ( B){tb=B->res; B->res=vb;fprintf ( stderr, "\nB: seq %d res %d", B->seq, B->res);}
  if ( C){tc=C->res; C->res=vc;fprintf ( stderr, "\nC: seq %d res %d", C->seq, C->res);}
  print_graph (A, 0);
  if ( A){A->res=ta;}
  if ( B){B->res=tb;}
  if ( C){C->res=tc;}
}
int check_link (CL_node ***G, int s1, int r1, int s2, int r2)
{
  CL_node *S;
  CL_node *E;

  S=G[s1][r1];
  E=G[s2][r2];
  while ( S->p)S=S->p;
  while ( S)
    {
      S=S->c;
      if ( S==E)return 1;
    }
  return 0;
}
