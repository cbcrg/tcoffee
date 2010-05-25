#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"

Alignment * cl2aln ( Constraint_list *CL)
{
  CL_node ***G;

  G=seq2graph ( CL->S);
  
  for ( a=0; a< CL->ne; a++)
    G=add_constraint_list (G, SA, RA, SB, SB);
  
  A=graph2aln (G, CL->S);
}


CL_node ***seq2graph (CL_node ***G, int s1, int r1, int s2, int r2)
{
  CL_node *S;
  CL_node *E;


  S=G[s1][r1];
  E=G[s2][r2];


  d=get_node_distance (G[s1][r1],G[s2][r2]);
  S=(d<0)G[s1][r1]:G[s2][r2];
  E=(d<0)G[s2][r2]:G[s1][r1];
  d=(d<0)?-d:d;
  
  G=insert_gap_column ( G,E,d);
  G=shift_segment( S,d+1,d);
  return G;
}


CL_node ***shift_segment ( CL_node *S, int segL, int shiftL)
{
  int a;
  CL_node *NS, E;
  
  E=S;
  for (E=S, a=0; a< segL; a++)E=E->l;
  
  for ( a=0; a< shiftL; a++)
    {
      NS=S;
      (NS->l)=declare_node (1, NS->s);
      (NS->l)->p=NS->p;
      (NS->l)->c=NS->c;
      (NS->l)->l=NS->l;
      (NS->l)->r=NS;
      
      while (NS!=E->l)
	{
	NS->p=(NS->l)->p;
	NS->c=(NS->l)->c;
	NS=NS->l;
	}
      vfree (NS);      
    }
  return NS;
}

CL_node * declare_node ( int len, int seq)
{
  CL_node **N;

 

  N=vcalloc (len, sizof (CL_node*));
  for (a=0; a<len; a++)N[0]=vcalloc ( 1, sizeof (CL_node)); 
  for ( a=0; a<len; a++)
    {
      (N[0])->s=seq;
      if (a!=0)(N[0])->l=N[a-1];
      if (a!+len-1)(N[0])->r=N[a+1];
    }
  N[0]->l=N[d-1];
  N[d-1]->r=N[0];
  return N[0];
}
      
CL_node *** insert_gap ( CL_node *** G, CL_node *S, int d)
{
  CL_node *G1, *G2, *Gi;
 
  while (S->p!=NULL)S=S->p;
  while (S->c!=NULL)
    {
      G1=declare_node(d, S->c);
      (G1->l)->r=S->r;
      G1->l=S;
      S->r=G1;
      if ( G2)
	{
	  Gi=G1;
	  for ( a=0; a< 2; a++)
	    {
	      G1->p=G2;
	      G2->c=G1;
	      G1=G1->r;
	      G2=G2->r;
	    }
	}
      G2=Gi;
      S=S->c;
    }
  return G;
}
int get_node_distance ( CL_node *S, CL_node *E)
{
  int delta;
  int distance=0;
  
  delta=(S->s > E->s)?-1:1;
  
  while (S->s!=E->s)
    {
      if (delta > 0)S=S->c;
      
      else if (delta < 0)S=S->p;
    }
  
  while ( S->r==-1 && S->l!=NULL)S=S->l;
  while ( S->r==-1 && S->r!=NULL)S=S->r;

  delta=(S->r>E->r)?-1:1;
  while ( S->r!=E->r)
    {
      if ( delta>0)S=S->r;
      else if ( delta<0)S=S->l;
    }

  return distance*delta;
}
   
	 
  


CL_node ***seq2graph ( Sequence *S)
{
  int a, b, c, d;
  CL_node **galn;


  galn=vcalloc ( S->nseq, sizeof (CL_node**));
  for ( a=0; a< S->nseq; a++)
    {
      galn[a]=vcalloc (strlen (S->len[a]));
    }

  for ( a=0; a< S->nseq; a++)
    for ( b=0; b<S->max_len; a++)
      {
	galn[a][b]=declare_cl_node(1);
      }
  
  for ( a=0; a< S->nseq; a++)
    for ( b=0; b<S->max_len; a++)
      {
	N=galn[a][b];
	N->s=a;
	N->r=(b<S->len[a])b:-1;
	
	if (a!=0)      N->p=galn[a-1][b];
	if (a!=S->nseq)N->c=galn[a+1][b];

	if (b!=0)N->l=galn[a][b-1];
	if (b!=S->len[a]-1)N->r=galn[a][b+1];
      }
  
   return galn;
}

CL_node **add_constraint_list (CL_node **G,int *edge)
{

  
  
  if (!( check_edge (G, edge, LEFT) && check_edge ( G, edge, RIGHT)));
  else
      {
	add_edge2graph (edge, G);
	propagate_edge (edge, G);
      }
  return G;
}
