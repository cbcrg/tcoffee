#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"


int cl2pair_list_ecl ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in);



/*******************************************************************************/
/*                idscore_pairseq: measure the % id without delivering thze aln*/                                                   
/*                                                                             */
/*	makes DP between the the ns[0] sequences and the ns[1]                 */
/*                                                                             */
/*	for MODE, see the function  get_dp_cost                                */
/*******************************************************************************/
int idscore_pairseq (char *s1, char *s2, int gop, int gep, int **m, char *comp_mode)
{
  int **I, **D, **M, *P;
  int i, j,l1, l2, l,score, id, igop,match;


  l1=strlen (s1); l2=strlen (s2);
  lower_string (s1); lower_string (s2);
  
  I=declare_int (6,l2+1);D=declare_int (6,l2+1);M=declare_int (6,l2+1);
  for (j=0; j<=l2; j++)
    {
      D[0][j]=gep*j;M[0][j]=2*gep*j;D[4][j]=0;
    }
  
  for (i=1; i<=l1; i++)
    {

      I[1][0]=i*gep;
      M[1][0]=2*i*gep;
      
      for (j=1; j<=l2; j++)
	{
	  score=m[s1[i-1]-'a'][s2[j-1]-'a'];	  
	  id=(s1[i-1]==s2[j-1])?1:0;
	  
	  igop=(i==l1 || j==l2)?0:gop;

	  if   ((D[0][j]+gep)>(M[0][j]+igop+gep))   {D[1][j]=D[0][j]+gep;      D[3][j]=D[2][j];        D[5][j]=D[4][j];}
	  else                                      {D[1][j]=M[0][j]+igop+gep; D[3][j]=M[2][j];        D[5][j]=M[4][j];}
	  
	  if ( (I[1][j-1]+gep)>(M[1][j-1]+igop+gep)){I[1][j]=I[1][j-1]+gep;      I[3][j]=I[3][j-1];    I[5][j]=I[5][j-1];}
	  else                                      {I[1][j]=M[1][j-1]+igop+gep; I[3][j]=M[3][j-1];    I[5][j]=M[5][j-1];}
	  
	  match=M[0][j-1]+score;
	  if (I[1][j]>match && I[1][j]>D[1][j])     {M[1][j]=I[1][j]           ; M[3][j]=I[3][j];      M[5][j]=I[5][j];}
	  else if (D[1][j]>match)                   {M[1][j]=D[1][j]           ; M[3][j]=D[3][j];      M[5][j]=D[5][j];}
	  else                                      {M[1][j]=match             ; M[3][j]=M[2][j-1]+id; M[5][j]=M[4][j-1]+1;}
	}
      P=I[0]; I[0]=I[1]; I[1]=P;
      P=I[2]; I[2]=I[3]; I[3]=P;
      P=I[4]; I[4]=I[5]; I[5]=P;
      
      P=D[0]; D[0]=D[1]; D[1]=P;
      P=D[2]; D[2]=D[3]; D[3]=P;
      P=D[4]; D[4]=D[5]; D[5]=P;
      
      P=M[0]; M[0]=M[1]; M[1]=P;
      P=M[2]; M[2]=M[3]; M[3]=P;
      P=M[4]; M[4]=M[5]; M[5]=P;
    }
 

  

  if ( strstr (comp_mode, "sim2"))
    {
       l=MIN(l1,l2);
       score=(l==0)?0:(M[2][l2]*100)/l;
    }
  else if ( strstr (comp_mode, "sim3"))
    {
       l=MAX(l1,l2);
       score=(l==0)?0:(M[2][l2]*100)/l;
    }
  else if ( strstr (comp_mode, "cov"))
    {
      l=MAX(l1,l2);
      score=(l==0)?0:((M[4][l2]*100)/l);
    }
  else
    {
      //default: simple sim
      l=M[4][l2];
      score=(l==0)?0:(M[2][l2]*100)/l;
    }      
  
  free_int (I, -1);
  free_int (D, -1);
  free_int (M, -1);
  
  return score;
}
	      
int test_pair_wise (Alignment *A, int *ns, int **l_s, Constraint_list *CL)
{
  int a,l0, l1, n;
  char buf[VERY_LONG_STRING];
  char *gap, *seq;
  
  l0=strlen (A->seq_al[l_s[0][0]]);
  l1=strlen (A->seq_al[l_s[1][0]]);
  
  n=(l0<5)?l0/2:5;
  gap=generate_null(l1-n);
  for (a=0;a<ns[0]; a++)
    {
      seq=A->seq_al[l_s[0][a]];
      sprintf (buf, "%s%s",seq, gap);
      sprintf (seq, "%s", buf);
    }
  vfree (gap);
  gap=generate_null(l0-n);
  
  for (a=0;a<ns[1]; a++)
    {
      seq=A->seq_al[l_s[1][a]];
      sprintf (buf, "%s%s",seq, gap);
      sprintf (seq, "%s", buf);
    }
  vfree(gap);
  

  A->len_aln=strlen (A->seq_al[l_s[0][0]]);
  A->score=A->score_aln=100;
  return 100;
}  

int idscore_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL)
{
  
  A->score_aln=A->score=idscore_pairseq (A->seq_al[l_s[0][0]], A->seq_al[l_s[1][0]], CL->gop, CL->gep,CL->M, "sim3");
  return A->score_aln;
}
int dp_max (int *trace, int n, ...);
int dp_max (int *trace, int n, ...)
{
  va_list ap;
  int a, v, t, best_v=0;
  
  va_start (ap, n);
  for (a=0; a< n; a++)
    {
      t=va_arg (ap, int);
      v=va_arg (ap, int);

      if (a==0)
	{
	  best_v=v;
	  trace[0]=t;
	}
      else
	{
	  if (v>best_v)
	    {
	      best_v=v;
	      trace[0]=t;
	    }
	}
    }
 
  return best_v;
}
int is_tied (int *trace, int n, ...);
int is_tied(int *trace, int n, ...)
{
  va_list ap;
  int a, v, t, best_v=0;
  int nties=0;
  
  va_start (ap, n);
  for (a=0; a< n; a++)
    {
      t=va_arg (ap, int);
      v=va_arg (ap, int);

      if (a==0)
	{
	  best_v=v;
	  trace[0]=t;
	}
      else
	{
	  if (v>best_v)
	    {
	      best_v=v;
	      trace[0]=t;
	    }
	}
    }
  va_end(ap);
  va_start (ap,n);
  for (a=0; a<n; a++)
    {
      t=va_arg (ap, int);
      v=va_arg (ap, int);
      if (v==best_v && trace[0]!=t)
	nties++;
    }
  va_end (ap);
  return nties;
}

void display_mat (int **M, int l1, int l2, char *title);
void display_mat (int **M, int l1, int l2, char *title)
{
  int a, b;
  
  fprintf ( stdout, "\n\nTitle %s\n", title);
  for ( a=0; a<=l1; a++)
    {
      fprintf ( stdout, "\n");
      for ( b=0; b<=l2; b++)
	fprintf ( stdout, "%3d ", M[a][b]);
    }
}
int glocal_pair_wise ( Alignment *A, int *ns, int **l_s, Constraint_list *CL)
{
  int ***t, ***m;
  int i,j, l1, l2, n, sub, trace,ntrace, a, b, c, score;
  int gop,rgop,tgop, gep, unmatch;
  int M1, M2, I1, D1, LEN;
  char **al, *char_buf, **aln;
  int **pos0;
  

  l1=strlen (A->seq_al[l_s[0][0]]);
  l2=strlen (A->seq_al[l_s[1][0]]);

  n=1;
  M1=n++;D1=n++;I1=n++;M2=n++;
  t=declare_arrayN(3, sizeof (int),n, l1+1, l2+1);
  m=declare_arrayN(3, sizeof (int),n, l1+1, l2+1);
  
  
  gop=CL->gop*SCORE_K;
  gep=CL->gep*SCORE_K;
  tgop=gop;
  unmatch=gep;
  
  pos0=aln2pos_simple ( A,-1, ns, l_s);
 
  
  for (j=1; j<=l2; j++)
    {
      m[D1][0][j]=gep*j;

      m[M1][0][j]=2*gep*j;
      m[M2][0][j]=4*gep*j;
    }
  
  
  for (i=1; i<=l1; i++)
    {
      m[I1][i][0]=i*gep;
      m[M2][i][0]=4*i*gep;
      m[M1][i][0]=2*i*gep;
                 
      for ( j=1; j<=l2; j++)
	{
	  rgop=(i==l1 || j==1)?0:gop;
	  rgop=gop;
	  sub=(CL->get_dp_cost) (A, pos0, ns[0], l_s[0], i-1, pos0, ns[1], l_s[1],j-1,CL);	
	  m[M1][i][j]=dp_max (&trace,4,M1,m[M1][i-1][j-1],I1, m[I1][i-1][j-1],D1,m[D1][i-1][j-1],M2,m[M2][i-1][j-1])+sub;
	  t[M1][i][j]=trace;
	  
	  m[D1][i][j]=dp_max (&trace,3, M1,m[M1][i][j-1]+rgop,D1, m[D1][i][j-1]+gep, M2, m[M2][i][j-1]);
	  t[D1][i][j]=trace;
	  
	  m[I1][i][j]=dp_max (&trace,3, M1,m[M1][i-1][j]+rgop, I1, m[I1][i-1][j]+gep, M2, m[M2][i-1][j]);
	  t[I1][i][j]=trace;

	  m[M2][i][j]=dp_max (&trace,4,M1,m[M1][i-1][j-1]+tgop,I1, m[I1][i-1][j-1]+tgop,D1,m[D1][i-1][j-1]+tgop,M2,m[M2][i-1][j-1])+unmatch;
	  t[M2][i][j]=trace;
	  
	}
	  
    }
  score=dp_max (&trace,4, M1,m[M1][l1][l2],D1,m[D1][l1][l2],I1, m[I1][l1][l2],M2,m[M2][l1][l2]);
  LEN=0;i=l1;j=l2;
  al=declare_char (2, l1+l2+1);
  

  trace=t[trace][i][j];
  while (!(i==0 &&j==0))
    {
  
      ntrace=t[trace][i][j];
      if (i==0)
	{
	  al[0][LEN]=0;
	  al[1][LEN]=1;
	  j--;
	  LEN++;
	}
      else if ( j==0)
	{
	  al[0][LEN]=1;
	  al[1][LEN]=0;
	  i--;
	  LEN++;
	}
      else if ( trace==M1)
	{
	  al[0][LEN]=1;
	  al[1][LEN]=1;
	  i--; j--;
	  LEN++;
	}
      else if ( trace==M2)
	{
	  al[0][LEN]=1;
	  al[1][LEN]=0;
	  LEN++;

	  al[0][LEN]=0;
	  al[1][LEN]=1;
	  LEN++;

	  i--; j--;
	  
	}
      else if ( trace==D1)
	{
	  al[0][LEN]=0;
	  al[1][LEN]=1;
	  j--;
	  LEN++;
	}
      else if ( trace == I1)
	{
	  al[0][LEN]=1;
	  al[1][LEN]=0;
	  i--;
	  LEN++;
	}
      trace=ntrace;     
      
    }
  
  invert_list_char ( al[0], LEN);
  invert_list_char ( al[1], LEN);	
  if ( A->declared_len<=LEN)A=realloc_aln2  ( A,A->max_n_seq, 2*LEN);	
  
  aln=A->seq_al;
  char_buf= vcalloc (LEN+1, sizeof (char));	
  for ( c=0; c< 2; c++)
    {
      for ( a=0; a< ns[c]; a++) 
	{		
	  int ch=0;
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
  free_arrayN((void *)m, 3);
  free_arrayN((void *)t, 3);
  vfree (char_buf);
  free_char (al, -1);
  return score;
}

int ** aln2local_penalties (Alignment *A, int n, int *ls, Constraint_list *CL, int **lg);
int ** aln2local_penalties (Alignment *A, int n, int *ls, Constraint_list *CL, int **lg)
{
  //adapted from gap_count in MAFFT V 5.5
  int p,s,l, c1, c2;
  int gep,gop;
  int open=3, close=4, gap=5;
  
  gop=CL->gop*SCORE_K;
  gep=CL->gep*SCORE_K;
  
  l=strlen (A->seq_al[ls[0]]);
  
  if (!lg)
    {
      lg=declare_int (6, l);
    }
  
  if ( read_array_size_new (lg[0])<l)
    {
      free_int (lg, -1);
      lg=declare_int (6, l);
    }
  
  for( s=0; s<n; s++ ) 
	{
	  c1='x';
	  for (p=0; p<l; p++)
	    {
	      c2=A->seq_al[ls[s]][p];
	      
	      if (c1!='-' && c2=='-')lg[open][p]++;
	      if (c1=='-' && c2!='-')lg[close][p]++;
	      if ( c1=='-')lg[gap][p]++;
	      c1=c2;
	    }
	}
  
  for (p=0; p<l; p++)
    {
      float go, gc, nn;
      nn=n;
      go=lg[open ][p];
      gc=lg[close][p];
     
    
      lg[GOP][p]=0.5*(1-(go/nn))*gop;
      lg[GCP][p]=0.5*(1-(gc/nn))*gop;
      //Checked locacal gep => gives low quality results
      lg[GEP][p]=gep;//(1-((float)lg[gap][p]/(float)n))*gep;
      lg[open][p]=lg[close][p]=lg[gap][p]=0;
      
    }

  return lg;
}
int free_gotoh_pair_wise_lgp()
{
  return gotoh_pair_wise_lgp (NULL, NULL, NULL, NULL);
}
int gotoh_pair_wise_lgp ( Alignment *A, int *ns, int **l_s, Constraint_list *CL)
{
  int i,j, li, lj, n, sub, trace,ntrace, a, b, c, score;
  int I, J;
  int M1, I1, D1, LEN;
  char **al, *char_buf, **aln;
  int **pos0, **pos;
  Alignment *Aln;
  
  int gop[2], gcp[2], gep[2];
  static int ***gpl, ***t, ***m;
  static int max_li, max_lj;
  
  

  //gotoh_pair_wise ( A, ns, l_s,CL);
  //ungap_sub_aln (A, ns[0], l_s[0]);
  //ungap_sub_aln (A, ns[1], l_s[1]);
  
  if (!A)
    {
      free_arrayN((void**)gpl, 3);
      free_arrayN((void**)t, 3);
      free_arrayN((void**)m, 3);
      max_li=max_lj=0;
      return 0;
    }

  I=0;J=1;

  
  li=strlen (A->seq_al[l_s[I][0]]);
  lj=strlen (A->seq_al[l_s[J][0]]);
  
  if ( !gpl)gpl=vcalloc ( 2, sizeof (int**));
  gpl[I]=aln2local_penalties (A,ns[I], l_s[I], CL,gpl[I]);
  gpl[J]=aln2local_penalties (A,ns[J], l_s[J], CL,gpl[J]);
  
  
  n=1;
  M1=n++;D1=n++;I1=n++;
  
  if ( li>max_li ||lj>max_lj )
    {
      free_arrayN((void**)t, 3);
      free_arrayN((void**)m, 3);

     
      max_li=li;
      max_lj=lj;
      t=declare_arrayN(3, sizeof (int),n, max_li+1, max_lj+1);
      m=declare_arrayN(3, sizeof (int),n, max_li+1, max_lj+1);
      
    }
  pos0=aln2pos_simple ( A,-1, ns, l_s);
 
  //Compatibility with Macro
  Aln=A;
  pos=pos0;
  
  for (j=1; j<=lj; j++)
    {
      gep[J]=gpl[J][GEP][j-1];
      m[D1][0][j]=gep[J]*j;
      m[I1][0][j]=m[D1][0][j]-1;
      m[M1][0][j]=m[D1][0][j]-1;
    }
  
  //D1: gap in sequence I
  //I1: gap in sequence J
  
  
  for (i=1; i<=li; i++)
    {
      gep[I]=gpl[I][GEP][i-1];
      gop[I]=gpl[I][GOP][i-1];
      gcp[I]=gpl[I][GCP][i-1];
      
      m[I1][i][0]=i*gep[I];
      m[D1][i][0]= m[I1][i][0]-1;
      m[M1][i][0]= m[I1][i][0]-1;
                 
     
      
      gop[I]=(i==1 || i==li )?0:gop[I];
      gcp[I]=(i==1 || i==li )?0:gcp[I];
      
      
      for ( j=1; j<=lj; j++)
	{
	  
	  gep[J]=gpl[J][GEP][j-1];
	  gop[J]=gpl[J][GOP][j-1];
	  gcp[J]=gpl[J][GCP][j-1];
	  
	  //gep[J]=gep[I]=(gep[J]+gep[I])/2;
	  //gop[J]=gop[I]=(gop[J]+gop[I])/2;
	  //gcp[J]=gcp[I]=(gcp[J]+gcp[I])/2;
	  

	  gop[J]=(j==1 || j==lj )?0:gop[J];
	  gcp[J]=(j==1 || j==lj )?0:gcp[J];
	  
	  
	  //sub=(CL->get_dp_cost) (A, pos0, ns[0], l_s[0], i-1, pos0, ns[1], l_s[1],j-1,CL);	
	  sub=TC_SCORE((i-1), (j-1));
	  
	  m[M1][i][j]=dp_max (&trace,3,M1,m[M1][i-1][j-1],I1, m[I1][i-1][j-1]+gcp[I],D1,m[D1][i-1][j-1]+gcp[J])+sub;
	  t[M1][i][j]=trace;
	  
	  
	  m[D1][i][j]=dp_max (&trace,2, M1,m[M1][i][j-1]+gop[J]+gep[J],D1, m[D1][i][j-1]+gep[J]);
	  t[D1][i][j]=trace;
	  
	  
	  m[I1][i][j]=dp_max (&trace,2, M1,m[M1][i-1][j]+gop[I]+gep[I],I1, m[I1][i-1][j]+gep[I]);
	  t[I1][i][j]=trace;
	  
	}
	  
    }
  score=dp_max (&trace,3, M1,m[M1][li][lj],D1,m[D1][li][lj],I1, m[I1][li][lj]);
  
  LEN=0;i=li;j=lj;
  al=declare_char (2, li+lj);
  

  trace=t[trace][i][j];
  while (!(i==0 &&j==0))
    {
  
      ntrace=t[trace][i][j];
     
      
      if (i==0)
	{
	  al[0][LEN]=0;
	  al[1][LEN]=1;
	  j--;
	  LEN++;
	}
      else if ( j==0)
	{
	  al[0][LEN]=1;
	  al[1][LEN]=0;
	  i--;
	  LEN++;
	}
      else if ( trace==M1)
	{
	  al[0][LEN]=1;
	  al[1][LEN]=1;
	  i--; j--;
	  LEN++;
	}
      else if ( trace==D1)
	{
	  al[0][LEN]=0;
	  al[1][LEN]=1;
	  j--;
	  LEN++;
	}
      else if ( trace == I1)
	{
	  al[0][LEN]=1;
	  al[1][LEN]=0;
	  i--;
	  LEN++;
	}
      trace=ntrace;     
      
    }
  
  invert_list_char ( al[0], LEN);
  invert_list_char ( al[1], LEN);	
  if ( A->declared_len<=LEN)A=realloc_aln  ( A,2*LEN+1);	
  
  aln=A->seq_al;
  char_buf= vcalloc (LEN+1, sizeof (char));	
  for ( c=0; c< 2; c++)
    {
      for ( a=0; a< ns[c]; a++) 
	{		
	  int ch=0;
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
  vfree (char_buf);
  free_char (al, -1);
  free_int (pos0, -1);
  return score;
}
/*******************************************************************************/
/*                GLOCAL 2                                                     */
/*                                                                             */
/*	makes DP between the the ns[0] sequences and the ns[1]                 */
/*                                                                             */
/*	for MODE, see the function  get_dp_cost                                */
/*******************************************************************************/
int glocal2_pair_wise (Alignment *IN,int*ns, int **ls,Constraint_list *CL)
{
  int a, b, s=0;
  Alignment *A, *R,*L;
  char *seq, *buf;
  
  buf=vcalloc (1000, sizeof (char));
  seq=vcalloc (1000, sizeof (char));
  
  A=copy_aln (IN,NULL);
  L=copy_aln (IN,NULL);
  R=copy_aln (IN,NULL);
  
  gotoh_pair_wise_sw (A, ns, ls, CL);
  
  HERE ("1");
  for (a=0; a<2; a++)
    {
      for (b=0; b<ns[a]; b++)
	{
	  s=ls[a][b];
	  sprintf ( seq,"%s", IN->seq_al[s]);
	  
	  seq[A->order[s][2]]='\0';
	  sprintf (L->seq_al[s], "%s", seq);
	  sprintf (R->seq_al[s], "%s", seq+A->order[s][3]+1);
	}
    }
  HERE ("2");
  print_sub_aln (A, ns, ls);
  gotoh_pair_wise(L, ns, ls, CL);
  print_sub_aln (L, ns, ls);
  gotoh_pair_wise(R, ns, ls, CL);
  print_sub_aln (R, ns, ls);
  
  IN=realloc_aln (IN, A->len_aln+L->len_aln+R->len_aln+1);
  for (a=0; a<2; a++)
    {
      for (b=0; b<ns[a]; b++)
	{
	  s=ls[a][b];
	  sprintf (IN->seq_al[s], "%s%s%s",L->seq_al[s], A->seq_al[s], R->seq_al[s]);
	}
    }
  IN->len_aln=strlen (IN->seq_al[s]);
  
  print_sub_aln (IN, ns, ls);
  vfree (seq); vfree (buf);
  free_aln (A); free_aln (L);free_aln (R);
  return IN->score_aln;
}


int gotoh_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL)
	{
/*******************************************************************************/
/*                NEEDLEMAN AND WUNSCH (GOTOH)                                 */
/*                                                                             */
/*	makes DP between the the ns[0] sequences and the ns[1]                 */
/*                                                                             */
/*	for MODE, see the function  get_dp_cost                                */
/*******************************************************************************/
	  

/*TREATMENT OF THE TERMINAL GAP PENALTIES*/
/*TG_MODE=0---> gop and gep*/
/*TG_MODE=1---> ---     gep*/
/*TG_MODE=2---> ---     ---*/


	    int TG_MODE;
	    int l_gop, l_gep;
	    int gop, gep;
	    int maximise;
/*VARIANLES FOR THE MULTIPLE SEQUENCE ALIGNMENT*/
	int a, b, i, j;

	int *cc;	
	int *dd,*ddg;
	int   e, eg;

	int lenal[2], len;
	int t, c=0,s, ch;
	int sub;
	int fop;
	int score=0;
        int **pos0;
	static char **al;
	char **aln;
	int ala, alb,LEN;
	char *buffer;
	char *char_buf;
/*trace back variables       */
	FILE       *long_trace=NULL;
	TRACE_TYPE *buf_trace=NULL;
	static TRACE_TYPE **trace;
	TRACE_TYPE k;
	TRACE_TYPE *tr;
	int long_trace_flag=0;
	int dim;
/********Prepare penalties*******/
	gop=CL->gop*SCORE_K;
	gep=CL->gep*SCORE_K;
	TG_MODE=CL->TG_MODE;
	maximise=CL->maximise;
	
	
/********************************/	
/*CLEAN UP AFTER USE*/
	if ( A==NULL)
	   {
	   free_int (trace,-1);
	   trace=NULL;
	   free_char (al,-1);
	   al=NULL;
	   return 0;
	   }		

/*DO MEMORY ALLOCATION FOR DP*/

	lenal[0]=strlen (A->seq_al[l_s[0][0]]);
	lenal[1]=strlen (A->seq_al[l_s[1][0]]);
	len= MAX(lenal[0],lenal[1])+1;
	buf_trace=vcalloc ( len, sizeof (TRACE_TYPE));	
	buffer=vcalloc ( 2*len, sizeof (char));	
        al=declare_char (2, 2*len);  
	
	char_buf= vcalloc (2*len, sizeof (char));	
	

	dd = vcalloc (len, sizeof (int));
	

	cc = vcalloc (len, sizeof (int));
	ddg=vcalloc (len, sizeof (int));
	

	
	if ( len>=MAX_LEN_FOR_DP)
	    {
	    long_trace_flag=1;
	    long_trace=vtmpfile();
	    }
	else
	    {
	   
	    dim=(trace==NULL)?0:read_size_int ( trace,sizeof (int*));	   
	    trace    =realloc_int ( trace,dim,dim,MAX(0,len-dim), MAX(0,len-dim));
	    }
	
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
		       
	pos0=aln2pos_simple ( A,-1, ns, l_s);


	cc[0]=0;		
	tr=(long_trace_flag)?buf_trace:trace[0];
	tr[0]=(TRACE_TYPE)1;
	for ( j=1; j<=lenal[1]; j++)tr[j]=(TRACE_TYPE)-1;
	if (long_trace_flag)fwrite (buf_trace, sizeof ( TRACE_TYPE),lenal[1]+1, long_trace);
	
	
	t=(TG_MODE==0)?gop:0;
	

	for (cc[0]=0,j=1; j<=lenal[1]; j++)
	    {
	    
	    l_gop=(TG_MODE==0)?gop:0;
	    l_gep=(TG_MODE==2)?0:gep;

	    cc[j]=t=t+l_gep;
	    dd[j]=  t+  gop;
	    }

	t=(TG_MODE==0)?gop:0;	
	
	for (i=1; i<=lenal[0];i++)
			{			
			tr=(long_trace_flag)?buf_trace:trace[i];
			s=cc[0];

			l_gop=(TG_MODE==0)?gop:0;
			l_gep=(TG_MODE==2)?0:gep;
			
			
			
			cc[0]=c=t=t+l_gep;
			e=t+  gop;
			tr[0]=(TRACE_TYPE)1;

			

			for (eg=0,j=1; j<=lenal[1];j++)
				{				   
				 
				  sub=(CL->get_dp_cost) (A, pos0, ns[0], l_s[0], i-1, pos0, ns[1], l_s[1],j-1,CL);	
				      
				/*get the best Insertion*/
				l_gop=(i==lenal[0] || i==1 )?((TG_MODE==0)?gop:0):gop;
				l_gep=(i==lenal[0] || i==1)?((TG_MODE==2)?0:gep):gep;
			

				if ( a_better_than_b ( e,c+l_gop, maximise))eg++;
				else eg=1;	
				e=best_of_a_b (e, c+l_gop, maximise)+l_gep;
				
				/*Get the best deletion*/
				l_gop=(j==lenal[1] || j==1)?((TG_MODE==0)?gop:0):gop;
				l_gep=(j==lenal[1] || j==1)?((TG_MODE==2)?0:gep):gep;
				

				if ( a_better_than_b ( dd[j], cc[j]+l_gop, maximise))ddg[j]++;
				else ddg[j]=1;
				dd[j]=best_of_a_b( dd[j], cc[j]+l_gop,maximise)+l_gep;
				


				c=best_int(3,maximise,&fop, e, s+sub,dd[j]);
				/*Chose Substitution for tie breaking*/
				if ( fop==0 && (s+sub)==e)fop=1;
				else if ( fop==2 && (s+sub)==dd[j])fop=1;
				/*Chose Deletion for tie breaking*/
				else if ( fop==2 && e==dd[j])fop=1;

				fop-=1;
				s=cc[j];
				cc[j]=c;	

	
				if ( fop<0)
					{tr[j]=(TRACE_TYPE)fop*eg;
					}
				else if ( fop>0)
				        {tr[j]=(TRACE_TYPE)fop*ddg[j];
					}
				else if (fop==0)
					{tr[j]=(TRACE_TYPE)0;	
					}					
				fop= -2;
				}
			if (long_trace_flag)
			    {
			    fwrite ( buf_trace, sizeof (TRACE_TYPE), lenal[1]+1, long_trace);
			    }
			}
	
	score=c;
	
        i=lenal[0];
	j=lenal[1];
	ala=alb=0;
	

	while (i>=0 && j>=0 && ((i+j)!=0))
			{
			if ( i==0)
				k=-1;
			else if ( j==0)
				k=1;
			else if ( j==0 && i==0)
				k=1;	
			else
			        {
				if (long_trace_flag)
				   {
				   fseek ( long_trace, sizeof (TRACE_TYPE)*((lenal[1]+1)*(i)+j),SEEK_SET);
				   fread ( &k, sizeof (TRACE_TYPE), 1, long_trace);
				   }
				else
				   {
				   
				   k=trace[i][j];
				   }
				}
				
				
			if (k==0)
				{
				
				al[0][ala++]=1;
				al[1][alb++]=1;
				i--;
				j--;
				}		
			else if (k>0)
				{
				
				for ( a=0; a< k; a++)
					{
					al[0][ala++]=1;
					al[1][alb++]=0;
					i--;
					}
				}
			else if (k<0)
				{
				
				for ( a=0; a>k; a--)
					{
					al[0][ala++]=0;
					al[1][alb++]=1;
					j--;
					}
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
	

	vfree ( cc);
	vfree (dd);		
	vfree (ddg);
	vfree (buffer);
	vfree (char_buf); 
	vfree (buf_trace);
	free_char ( al, -1);
	free_int (pos0, -1);
	if ( long_trace_flag)fclose (long_trace);	



	return score;
	}
     

int get_transition_cost (Alignment *A, int **posi, int ni, int *li, int i, int **posj, int nj, int *lj, int j,Constraint_list *CL);
int gotoh_pair_wise_lgp_sticky ( Alignment *A, int *ns, int **l_s, Constraint_list *CL)
{
  int i,j, li, lj, n, sub, trace,ntrace, a, b, c, score;
  int I, J;
  int M1, I1, D1, LEN;
  char **al, *char_buf, **aln;
  int **pos0;
  
  int gop[2], gcp[2], gep[2];
  static int ***gpl, ***t, ***m;
  static int max_li, max_lj;
  
  

  //gotoh_pair_wise ( A, ns, l_s,CL);
  //ungap_sub_aln (A, ns[0], l_s[0]);
  //ungap_sub_aln (A, ns[1], l_s[1]);
  
  I=0;J=1;

  
  li=strlen (A->seq_al[l_s[I][0]]);
  lj=strlen (A->seq_al[l_s[J][0]]);
  
  if ( !gpl)gpl=vcalloc ( 2, sizeof (int**));
  gpl[I]=aln2local_penalties (A,ns[I], l_s[I], CL,gpl[I]);
  gpl[J]=aln2local_penalties (A,ns[J], l_s[J], CL,gpl[J]);
  
  
  n=1;
  M1=n++;D1=n++;I1=n++;
  
  if ( li>max_li ||lj>max_lj )
    {
      free_arrayN((void**)t, 3);
      free_arrayN((void**)m, 3);

     
      max_li=li;
      max_lj=lj;
      t=declare_arrayN(3, sizeof (int),n, max_li+1, max_lj+1);
      m=declare_arrayN(3, sizeof (int),n, max_li+1, max_lj+1);
      
    }
  pos0=aln2pos_simple ( A,-1, ns, l_s);
 
  
  for (j=1; j<=lj; j++)
    {
      gep[J]=gpl[J][GEP][j-1];
      m[D1][0][j]=gep[J]*j;
      m[I1][0][j]=m[D1][0][j]-1;
      m[M1][0][j]=m[D1][0][j]-1;
    }
  
  //D1: gap in sequence I
  //I1: gap in sequence J
  
  
  for (i=1; i<=li; i++)
    {
      gep[I]=gpl[I][GEP][i-1];
      gop[I]=gpl[I][GOP][i-1];
      gcp[I]=gpl[I][GCP][i-1];
      
      m[I1][i][0]=i*gep[I];
      m[D1][i][0]= m[I1][i][0]-1;
      m[M1][i][0]= m[I1][i][0]-1;
                 
     
      
      gop[I]=(i==1 || i==li )?0:gop[I];
      gcp[I]=(i==1 || i==li )?0:gcp[I];
      
      
      for ( j=1; j<=lj; j++)
	{
	  int transition;
	  
	  gep[J]=gpl[J][GEP][j-1];
	  gop[J]=gpl[J][GOP][j-1];
	  gcp[J]=gpl[J][GCP][j-1];
	  
	  //gep[J]=gep[I]=(gep[J]+gep[I])/2;
	  //gop[J]=gop[I]=(gop[J]+gop[I])/2;
	  //gcp[J]=gcp[I]=(gcp[J]+gcp[I])/2;
	  

	  gop[J]=(j==1 || j==lj )?0:gop[J];
	  gcp[J]=(j==1 || j==lj )?0:gcp[J];
	  
	  
	  sub=(CL->get_dp_cost) (A, pos0, ns[0], l_s[0], i-1, pos0, ns[1], l_s[1],j-1,CL);
	  transition=get_transition_cost (A, pos0, ns[0], l_s[0], i-1, pos0, ns[1], l_s[1],j-1,CL);

	  m[M1][i][j]=dp_max (&trace,3,M1,m[M1][i-1][j-1]+transition,I1, m[I1][i-1][j-1]+gcp[I],D1,m[D1][i-1][j-1]+gcp[J])+sub;
	  t[M1][i][j]=trace;
	  
	  
	  m[D1][i][j]=dp_max (&trace,2, M1,m[M1][i][j-1]+gop[J]+gep[J],D1, m[D1][i][j-1]+gep[J]);
	  t[D1][i][j]=trace;
	  
	  
	  m[I1][i][j]=dp_max (&trace,2, M1,m[M1][i-1][j]+gop[I]+gep[I],I1, m[I1][i-1][j]+gep[I]);
	  t[I1][i][j]=trace;
	  
	}
	  
    }
  score=dp_max (&trace,3, M1,m[M1][li][lj],D1,m[D1][li][lj],I1, m[I1][li][lj]);
  
  LEN=0;i=li;j=lj;
  al=declare_char (2, li+lj);
  

  trace=t[trace][i][j];
  while (!(i==0 &&j==0))
    {
  
      ntrace=t[trace][i][j];
     
      
      if (i==0)
	{
	  al[0][LEN]=0;
	  al[1][LEN]=1;
	  j--;
	  LEN++;
	}
      else if ( j==0)
	{
	  al[0][LEN]=1;
	  al[1][LEN]=0;
	  i--;
	  LEN++;
	}
      else if ( trace==M1)
	{
	  al[0][LEN]=1;
	  al[1][LEN]=1;
	  i--; j--;
	  LEN++;
	}
      else if ( trace==D1)
	{
	  al[0][LEN]=0;
	  al[1][LEN]=1;
	  j--;
	  LEN++;
	}
      else if ( trace == I1)
	{
	  al[0][LEN]=1;
	  al[1][LEN]=0;
	  i--;
	  LEN++;
	}
      trace=ntrace;     
      
    }
  
  invert_list_char ( al[0], LEN);
  invert_list_char ( al[1], LEN);	
  if ( A->declared_len<=LEN)A=realloc_aln  ( A,2*LEN+1);	
  
  aln=A->seq_al;
  char_buf= vcalloc (LEN+1, sizeof (char));	
  for ( c=0; c< 2; c++)
    {
      for ( a=0; a< ns[c]; a++) 
	{		
	  int ch=0;
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
  vfree (char_buf);
  free_char (al, -1);
  free_int (pos0, -1);
  return score;
}
int get_transition_cost (Alignment *A, int **posi, int ni, int *li, int i, int **posj, int nj, int *lj, int j,Constraint_list *CL)
{
  /*counts the number of identical transitions between position i-1, i and j-1..j*/
  float t=0;
  int a,s;
  Sequence *S;
  
  if (i==0 || j==0)return 0;
  
  for (a=0; a<ni; a++)
    {
      s=li[a];
      if (posi[s][i]<0 || posi[s][i-1]<0)continue;
      if (S->seq[li[a]][i-1]==S->seq[li[a]][i-1])t++;
    }
  
  for (a=0; a<nj; a++)
    {
      s=lj[a];
      if (posj[s][j]<0 || posj[s][j-1]<0)continue;
      if (S->seq[li[a]][j-1]==S->seq[li[a]][j-1])t++;
    }

  t=(t*10)/(float)(ni+nj);
  return t;
}
/*******************************************************************************/
/*                idscore_pairseq: measure the % id without delivering thze aln*/                                                   
/*                                                                             */
/*	makes DP between the the ns[0] sequences and the ns[1]                 */
/*                                                                             */
/*	for MODE, see the function  get_dp_cost                                */
/*******************************************************************************/
int cl2pair_list ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in, int mode, int ndiag);
int cl2pair_list_ref ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in);
int cl2pair_list_ecf ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in);
int cl2pair_list_diag ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in, int add);
int cl2list_borders   (Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in);
int cl2diag_cap (Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in); //add one element at the end of each segment so that they can be joined
int** cl2sorted_diagonals   ( Alignment *A, int *ns, int **ls, Constraint_list *CL);
int** cl2sorted_diagonals_mat  ( Alignment *A, int *ns, int **ls, Constraint_list *CL);
int** cl2sorted_diagonals_cs   ( Alignment *A, int *ns, int **ls, Constraint_list *CL);
int list2nodup_list (Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in);
int fill_matrix ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in);
int cl2pair_list ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in, int mode, int ndiag)
{
  int v;
  if (!A)
    {
      free_int (list_in[0], -1); list_in[0]=NULL;
      n_in[0]=0;
    }
  
  cl2list_borders(A, ns, ls, CL, list_in, n_in);
  
  if ( mode==0)
    v=cl2pair_list_ref (A, ns, ls, CL, list_in, n_in);
  else if (mode==1)
    v=cl2pair_list_ecl (A, ns, ls, CL, list_in, n_in);
  else if (mode==2)
    v=cl2pair_list_diag (A, ns, ls, CL, list_in, n_in,ndiag); //add diagonals
  
  cl2diag_cap (A, ns, ls, CL, list_in, n_in);
  //fill_matrix (A, ns, ls, CL, list_in, n_in);//Fill matrix with 0s
  sort_list_int (list_in[0],7, 1, 0, n_in[0]-1);
  list2nodup_list (A, ns, ls, CL, list_in, n_in);
  return v;
  
}
int fill_matrix( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list, int *n)
{
  int a, b, l1, l2, n2=0;
  int score;
  int **pos;
  int max_n;
  if (!A) return 0;
  pos=aln2pos_simple ( A,-1, ns, ls);
  l1=strlen (A->seq_al[ls[0][0]]);
  l2=strlen (A->seq_al[ls[1][0]]);
  
  max_n=read_array_size (list[0], sizeof (int));
  
  for (a=0; a<=l1; a++)
    for (b=0; b<=l2; b++)
      {
	score=0;
	score=(a==0 || b==0)?0:slow_get_dp_cost (A, pos, ns[0], ls[0],a-1, pos, ns[1], ls[1], b-1, CL);
	if ( score>0 && a!=0 && b!=0 && a!=l1 && b!=l2)
	  {
	    if (n[0]==max_n){max_n+=1000;list[0]=vrealloc (list[0], max_n*sizeof (int*));}
	    if (!list[0][n[0]])list[0][n[0]]=vcalloc (7, sizeof (int));
	    list[0][n[0]][0]=a;
	    list[0][n[0]][1]=b;
	    list[0][n[0]][3]=(l1-a)+b;
	    list[0][n[0]][2]=score;
	    if ( a!=0 && b!=0 && a!=l1 && b!=l2)
	      {
		n2++;
	      }
	    n[0]++;
	  }

      }

  return n[0];
  }
int list2nodup_list ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in)
{
  int **list;
  int n, a, b, c;
  
  list=list_in[0];
  n=n_in[0];
  
  if ( !A)return 0;
  
  for (b=a=1; a<n; a++)
    {
      if (list[a][0]==list[b-1][0] && list[a][1]==list[b-1][1])
	{
	  //HERE ("Duplicate");
	   list[b-1][2]=MAX(list[b-1][2],list[a][2]);
	}
      else
	{
	  for (c=0; c<4; c++)list[b][c]=list[a][c];
	  b++;
	}
	
	}
  n_in[0]=b;
  list_in[0]=list;
  return b;
}
int** cl2sorted_diagonals   ( Alignment *A, int *ns, int **ls, Constraint_list *CL)
{
  if ( CL && CL->residue_index)return cl2sorted_diagonals_cs   ( A, ns, ls, CL);
  else return cl2sorted_diagonals_mat   ( A, ns, ls, CL);
}

static int kword;
static char **warray;
int cmp_word ( const int**a, const int**b);
int ** seq2index_list ( Sequence *S, int k);
int** cl2sorted_diagonals_mat  ( Alignment *A, int *ns, int **ls, Constraint_list *CL)
{
  
  int a,b,c,d, comp, k, l1, l2, ndiag;
  int  **diag;
  
  
  static char **alp, alps=5;
  char *buf1, *buf2;
  
  if (!A)return NULL;
  if ( !alp)
    alp=make_group_aa_upgma ("blosum62mt",alps);
  
  
  k=2;
  l1=strlen (A->seq_al[ls[0][0]]);
  buf1=vcalloc ( l1+1, sizeof (char));
  l2=strlen (A->seq_al[ls[1][0]]);
  buf2=vcalloc ( l2+1, sizeof (char));
  
  ndiag=l1+l2;
  diag=declare_int (ndiag+3,2);
  for (a=0; a<=ndiag; a++) diag[a][0]=a;
  vfree (diag[ndiag+1]);
  diag[ndiag+1]=NULL;
 
  for ( a=0; a<ns[0]; a++)
    {
      sprintf (buf1, "%s",A->seq_al[ls[0][a]]);
      lower_string (buf1);
      string_convert (buf1, alps, alp);
      for (b=0; b<ns[1]; b++)
	{
	  sprintf (buf2,"%s", A->seq_al[ls[1][b]]);
	  lower_string (buf2);
	  string_convert (buf2, alps, alp);
	  for (c=0; c<l1-k; c++)
	    {
	      if ( strnrchr (buf1+c, '-', k))continue;
	      for (d=0; d<l2-k; d++)
	      {
		if ( strnrchr (buf2+d, '-', k))continue;
		comp=strncmp (buf1+c,buf2+d, k);
		diag[l1-c+d][1]+=(comp==0)?1:0;
	      }
	    }
	}
    }
  HERE ("dONE");
  /*
   max_len=MAX(l1, l2);
   for (a=1; a<ndiag; a++)
    {
      int l, d, p1_0, p1_1, p2_0, p2_1;
      d=diag[a][0];

      if (d<=l1)l=MIN(d,l2);
      else l=MIN(((l1+l2)-d),l1); 

      diag[a][1]=(diag[a][1]*1000)/max_len;
      diag[a][1]=(float)((float)diag[a][1]*(((float)max_len)/(float)l));
    }
  */
  sort_int_inv (diag, 2, 1,0,ndiag);
  
  return diag;
}







int ** seq2index_list ( Sequence *S, int k)
{
  int **list,**mlist=NULL;
  int a, b,c,ml, n, l, e, s, max=0, nm=0;
  char *cw;
  

  for (ml=0,a=0; a<S->nseq; a++)ml+=strlen (S->seq[a]);
  list=declare_int (ml+1, 2);
  
  for (n=0,a=0; a<S->nseq; a++)
    {
      l=strlen (S->seq[a])-k;
      for ( b=0; b<l; b++, n++)
	{
	  list[n][0]=a;
	  list[n][1]=b;
	}
    }
  list[n][0]=-1;
  
  warray=S->seq;
  kword=k;
  qsort (list, n, sizeof (long**), (int(*)(const void*,const void*))(cmp_word));
    
  cw=NULL;
  e=s=0;
  nm=0;
  for (a=0; a<=n; a++)
    {
      int s1, s2, r1, r2;
      if (!cw ||a==n|| strncmp (warray[list[a][0]]+list[a][1],cw, k)!=0)
	{
	  if (a<n)cw=warray[list[a][0]]+list[a][1];
	  for (b=s; b<a-1; b++)
	    for (c=b+1; c<a; c++)
	      {

		if (list[b][0]<list[c][0])
		  {s1=list[b][0];
		    r1=list[b][1];
		    s2=list[c][0];
		    r2=list[c][1];
		  }
		else
		  {
		    s2=list[b][0];
		    r2=list[b][1];
		    s1=list[c][0];
		    r1=list[c][1];
		  }
		
		if (s1==s2)continue;
		else
		  {
		    if (nm>=max){max+=1000; mlist=vrealloc (mlist, max*sizeof (int*));}
		    mlist[nm]=vcalloc (4, sizeof (int));
		    mlist[nm][0]=s1;
		    mlist[nm][1]=s2;
		    mlist[nm][2]=r1;
		    mlist[nm][3]=r2;
		    nm++;
		  }
		s=a;
	      }
	}
    }
  
  if (nm>=max){max+=1000;mlist=vrealloc (mlist, max*sizeof (int));}
  sort_list_int ( mlist,4,1, 0, nm-1);
  return mlist;
}
int cmp_word ( const int**a, const int**b)
{
  int c;

  c=strncmp (warray[a[0][0]]+a[0][1], warray[b[0][0]]+b[0][1], kword);

  
  if (c) return c;
  else 
    {
      for (c=0; c<2; c++)
	{
	  if ( a[0][c]>b[0][c])return   1;
	  else if (a[0][c]<b[0][c])return   -1;
	}
    }
  return 0;
}
 
int** cl2sorted_diagonals_cs   ( Alignment *A, int *ns, int **ls, Constraint_list *CL)
{
  int p1, p2, s1, s2, r1, r2;
  int a, b, l1, l2;
  int **pos;
  int **diag;
  int ndiag;
  int diag_i;

  
  pos=aln2pos_simple ( A,-1, ns, ls);
  
  l1=strlen (A->seq_al[ls[0][0]]);
  l2=strlen (A->seq_al[ls[1][0]]);
  
  

  ndiag=l1+l2;
  diag=declare_int (ndiag+3, 2);
  
  for (a=1; a<=ndiag; a++)diag[a][0]=a;
    
  for (p1=1; p1<=l1; p1++)
    {
      for (p2=1; p2<=l2; p2++)
	{
	  for (a=0; a<ns[0]; a++)
	    {
	      s1=ls[0][a];
	      r1=pos[s1][p1-1];
	      if (r1<=0)continue;
	      for (b=0; b<ns[1]; b++)
		{
		  s2=ls[1][b];
		  r2=pos[s2][p2-1];
		  if (r2<=0)continue;
		  diag_i=(l1-p1)+p2;
		  
		  diag[diag_i][1]+=residue_pair_extended_list_raw (CL,s1, r1-1,s2, r2-1);
		}
	    }
	}
    }
  
  sort_int_inv (diag, 2, 1,0,ndiag);

  vfree (diag[ndiag+1]);
  diag[ndiag+1]=NULL;
  free_int (pos, -1);
  return diag;
} 
int** cl2sorted_diagonals_cs_old2   ( Alignment *A, int *ns, int **ls, Constraint_list *CL)
{
  int p1, p2, si, s, r, t_s, t_r;
  int a, l1, l2;
  int *sl2, **pos,**inv_pos;
  int **diag;
  int ndiag;
  int diag_i;
  int max_len;
  
  pos=aln2pos_simple ( A,-1, ns, ls);
  inv_pos=vcalloc ((CL->S)->nseq, sizeof (int*));
  for (a=0; a<ns[1]; a++)inv_pos[ls[1][a]] =seq2inv_pos(A->seq_al[ls[1][a]]);
  
  l1=strlen (A->seq_al[ls[0][0]]);
  l2=strlen (A->seq_al[ls[1][0]]);
  sl2=vcalloc ((CL->S)->nseq, sizeof (int));
  
  for (a=0;a<ns[1]; a++)sl2[ls[1][a]]=1;
  
  ndiag=l1+l2;
  diag=declare_int (ndiag+3, 2);
  
  for (a=1; a<=ndiag; a++)diag[a][0]=a;
  for (p1=0; p1<=l1; p1++)
    {
      for (si=0;p1>0 && si<ns[0]; si++)
	{
	  s=ls [0][si];
	  r=pos[s][p1-1];
	  
	  for (a=1; r>0 && a<CL->residue_index[s][r][0];a+=ICHUNK)
	    {
	      
	      t_s=CL->residue_index[s][r][a+SEQ2];
	      t_r=CL->residue_index[s][r][a+R2];
	      
	      if (sl2[t_s])
		{
		  p2=inv_pos[t_s][t_r];
		  diag_i=(l1-p1)+p2;
		  diag[diag_i][1]+=(CL->get_dp_cost) (A, pos, ns[0], ls[0], p1-1, pos, ns[1], ls[1],p2-1,CL);	
		}
	    }
	}
    }
  max_len=MAX(l1, l2);
  for (a=1; a<ndiag; a++)
    {
      int l, d;
      d=diag[a][0];

      if (d<=l1)l=MIN(d,l2);
      else l=MIN(((l1+l2)-d),l1); 

      diag[a][1]/=max_len;
      diag[a][1]=(float)((float)diag[a][1]*(((float)max_len)/(float)l));
    }
  sort_int_inv (diag, 2, 1,0,ndiag);

  vfree (diag[ndiag+1]);
  diag[ndiag+1]=NULL;
  free_int (pos, -1);
  free_int (inv_pos, -1);
  vfree (sl2);
  
  return diag;
}
int cl2list_borders   (Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in)
{
  int maxlen, a,n, p1, p2, l1, l2;
  int **list;
  int **pos;
  if (!A)return 0;
  

  if ( list_in[0] && list_in[0][0]==0)
	  return read_array_size (list_in[0], sizeof (int*));
  

  
  list=list_in[0];
  n=n_in[0];
  if (!list)maxlen=0;
  else maxlen=read_array_size (list, sizeof (int*));
  

  l1=strlen (A->seq_al[ls[0][0]]);
  l2=strlen (A->seq_al[ls[1][0]]);
//    pos=aln2pos_simple ( A,-1, ns, ls);

  for (p1=0; p1<=l1; p1++)
    {
      if (p1==0 || p1==l1)
	{
	  for (p2=0; p2<=l2; p2++)
	    {
	      if (n==maxlen){maxlen+=1000;list=vrealloc (list,maxlen*sizeof (int*));}
	      if (!list[n])list[n]=vcalloc (7, sizeof (int));
	      list[n][0]=p1;
	      list[n][1]=p2;
	      list[n][3]=(l1-(p1))+(p2);
	      //list[n][2]=(p1==0||p2==0)?0:(CL->get_dp_cost) (A, pos, ns[0], ls[0], p1-1, pos, ns[1], ls[1],p2-1,CL);;
	      list[n][2]=(CL->gep)*SCORE_K*p2;
	      n++;
	    }
	}
      else
	{
	  for (a=0; a<2; a++)
	    {
	      p2=(a==0)?0:l2;
	      if (n==maxlen){maxlen+=1000;list=vrealloc (list,maxlen*sizeof (int*));}
	      if (!list[n])list[n]=vcalloc (7, sizeof (int));
	      list[n][0]=p1;
	      list[n][1]=p2;
	      list[n][3]=(l1-(p1))+(p2);
	      //list[n][2]=(p1==0||p2==0)?0:(CL->get_dp_cost) (A, pos, ns[0], ls[0], p1-1, pos, ns[1], ls[1],p2-1,CL);;
	      list[n][2]=(CL->gep)*SCORE_K*p1;
	      n++;
	    }
	}
    }
//   free_int (pos, -1);
  list_in[0]=list;
  n_in[0]=n;
  return read_array_size (list, sizeof (int*));
}

int cl2diag_cap (Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in)
{
  int **list;
  int n, in, a, b, al1, al2;
  int max_n;
  int cap=0;
  
  if (!A) return 0;
  
  al1=strlen (A->seq_al[ls[0][0]]);
  al2=strlen (A->seq_al[ls[1][0]]);
  
  list=list_in[0];
  n=n_in[0];
  max_n=read_array_size (list, sizeof (int*));
  
  
  
  
  for (a=0; a< n; a++)
    {
      b=list[a][3];
      list[a][3]=list[a][0];
      list[a][0]=b;
      
    }
  sort_list_int (list, 4, 1, 0, n-1);
  for (a=0; a< n; a++)
    {
      b=list[a][3];
      list[a][3]=list[a][0];
      list[a][0]=b;
    }
  

  in=n;
  
  for (a=0; a<in; a++)
    {
      int i, j, pi, pj, ni, nj;
      if (list[a][2]==0)continue;
      i=list[a][0];
      j=list[a][1];
      
      if (a==0){pi=-10;pj=-10;}
      else {pi=list[a-1][0];pj=list[a-1][1];}
      
      if (a==in-1){ni=-10; nj=-10;}
      else {ni=list[a+1][0]; nj=list[a+1][1];}
      
      
      if (i==0 || j==0);
      else if ( i==pi || j==pj);
      else if ( i-pi!=1 || j-pj!=1)
	{
	  
	   if (n==max_n){max_n+=1000;list=vrealloc (list, max_n*sizeof (int*));}
	   if (!list[n])list[n]=vcalloc (7, sizeof (int));
	   
	   list[n][0]=i-1;
	   list[n][1]=j-1;
	   list[n][3]=list[a][3];
	   list[n][2]=cap;
	   n++;
	}
    
    
      if (i==al1 || j==al2);
      else if ( i==ni || j==nj);
      else if ( ni-i!=1 || nj-j!=1)
	{
	  
	  if (n==max_n){max_n+=1000;list=vrealloc (list, max_n*sizeof (int*));}
	  if (!list[n])list[n]=vcalloc (7, sizeof (int));
	  
	  list[n][0]=i+1;
	  list[n][1]=j+1;
	  list[n][3]=list[a][3];
	  list[n][2]=cap;
	  n++;
	}
    
    }
  list_in[0]=list;
  n_in[0]=n;
  return n;
}
	  
int cl2pair_list_diag_mat ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in, int add );
int cl2pair_list_diag_cl ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in, int add );	  
int cl2pair_list_diag ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in, int add )
{
  if (CL->residue_index)return cl2pair_list_diag_cl (A, ns, ls, CL, list_in, n_in, add);
  else return cl2pair_list_diag_mat (A, ns, ls, CL, list_in, n_in, add);
}
int cl2pair_list_diag_mat ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in, int add )
{
  int p1, p2, n,d;
  int a, l1, l2;
  int set=0;
  static int **pos;
  int max_n;
  static int **diag;
  int **list;

  if (A==NULL)
    {
      free_int (pos, -1);pos=NULL;
      free_int (diag, -1);diag=NULL;
      //free_int (list_in[0], -1); list_in[0]=NULL;
      return 0;
    }
  
  if ( !pos)
    {
      pos=aln2pos_simple ( A,-1, ns, ls);
      diag=cl2sorted_diagonals (A,ns,ls,CL);
    }
  
  list=list_in[0];
  n=n_in[0];
  max_n=read_array_size (list, sizeof (int**));
  
  
  
  l1=strlen (A->seq_al[ls[0][0]]);
  l2=strlen (A->seq_al[ls[1][0]]);
  
  d=0;
  if ( add)
    {
      while ( diag[d] && diag[d][1]==-1)d++;
      add+=d;
    }
  else
    {
      d=0;
      while (diag[add++]);
    }
  
  HERE ("Add %d diagonals, starts %d N=%d", add, d, n);
  
  for (d=0; d<add && diag[d]; d++)
    {
      int p1_0, p2_0;
      
      set=1;

      HERE ("\t S=%d", diag[d][1]);
      
      p1_0=MAX(0,l1-diag[d][0]);
      p2_0=MAX(0,diag[d][0]-l1);
      diag[d][1]=-1;
      
      for (p1=p1_0, p2=p2_0; p1<=l1 && p2<=l2; p1++,p2++)
	{
	  if (!BORDER(p1,l1, p2,l2) )
	    {

	      if (n==max_n){max_n+=1000;list=vrealloc (list, max_n*sizeof (int*));}
	      if (!list[n])list[n]=vcalloc (7, sizeof (int));
	      
	      list[n][0]=p1;
	      list[n][1]=p2;
	      list[n][3]=(l1-(p1))+(p2);
	      list[n][2]=(CL->get_dp_cost) (A, pos, ns[0], ls[0], p1-1, pos, ns[1], ls[1],p2-1,CL);
	      n++;
	    }
	}
    }
  HERE ("Addition Finished n=%d", n);
  if (!set) return 0;
  sort_list_int (list,4, 1, 0, n-1);
  
  list_in[0]=list;
  n_in[0]=n;
  HERE ("\nN=%d r=%.3f [l1=%d l2=%d]", n, (float)n/(float)(l1*l2), l1, l2);
  return max_n;
}

int cl2pair_list_diag_cl ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in, int add )
{
  int p1, p2,n,d;
  int l1, l2;
  int score, set=0;
  static int **pos;
  int max_n;
  static int **diag;
  int **list;

  

  if (A==NULL)
    {
      free_int (pos, -1);pos=NULL;
      free_int (diag, -1);diag=NULL;
      //free_int (list_in[0], -1); list_in[0]=NULL;
      return 0;
    }
  
  if ( !pos)
    {
      pos=aln2pos_simple ( A,-1, ns, ls);
      diag=cl2sorted_diagonals (A,ns,ls,CL);
    }
  
  list=list_in[0];
  n=n_in[0];
  max_n=read_array_size (list, sizeof (int**));
  
  
  
  l1=strlen (A->seq_al[ls[0][0]]);
  l2=strlen (A->seq_al[ls[1][0]]);
  if ( add==0)add=l1+l2;
  d=0;
  while ( diag[d] && diag[d][1]==-1)d++;
  HERE ("Add %d diagonals, starts %d N=%d", add, d, n);
  add+=d;
  for (d; d<add && diag[d]; d++)
    {
      int p1_0, p2_0;
      if (diag[d][1]==0)continue;
      set=1;

      HERE ("\t S=%d", diag[d][1]);
      
      p1_0=MAX(0,l1-diag[d][0]);
      p2_0=MAX(0,diag[d][0]-l1);
      diag[d][1]=-1;
      
      for (p1=p1_0, p2=p2_0; p1<=l1 && p2<=l2; p1++,p2++)
	{
	  if (!BORDER(p1,l1, p2,l2) && (score=(CL->get_dp_cost) (A, pos, ns[0], ls[0], p1-1, pos, ns[1], ls[1],p2-1,CL))!=0)
	    {

	      if (n==max_n){max_n+=1000;list=vrealloc (list, max_n*sizeof (int*));}
	      if (!list[n])list[n]=vcalloc (7, sizeof (int));
	      
	      list[n][0]=p1;
	      list[n][1]=p2;
	      list[n][3]=(l1-(p1))+(p2);
	      list[n][2]=score;
	      n++;
	    }
	}
    }
  HERE ("Addition Finished n=%d", n);
  if (!set) return 0;
  sort_list_int (list,4, 1, 0, n-1);
  
  list_in[0]=list;
  n_in[0]=n;
  HERE ("\nN=%d r=%.3f [l1=%d l2=%d]", n, (float)n/(float)(l1*l2), l1, l2);
  return max_n;
}
  
int cl2pair_list_ecl_norm ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in);
int cl2pair_list_ecl_raw ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in);
int cl2pair_list_ecl_pc ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in);
int cl2pair_list_ecl_noext_raw ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in);
int cl2pair_list_ecl_rna2 ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in);

int cl2pair_list_ecl_rawquad ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in);
int cl2pair_list_ecl ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in)
{
  int mode=5;
  
  if (atoigetenv ("ECLMODE4TC")){HERE ("****ECLMODE4TC:%d",mode=atoigetenv ("ECLMODE4TC"));}
  
  if ( mode==1)return cl2pair_list_ecl_norm         (A, ns, ls, CL, list_in, n_in);
  else if ( mode==2)return cl2pair_list_ecl_raw     (A, ns, ls, CL, list_in, n_in);
  else if ( mode==3)return cl2pair_list_ecl_rawquad (A, ns, ls, CL, list_in, n_in);
  else if ( mode==4)return cl2pair_list_ecl_noext_raw (A, ns, ls, CL, list_in, n_in);
  else if ( mode==5)return cl2pair_list_ecl_pc     (A, ns, ls, CL, list_in, n_in);
  else if ( mode==6)return cl2pair_list_ecl_ext_pc     (A, ns, ls, CL, list_in, n_in);
  
}
int cl2pair_list_ecl_noext_raw ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in)
{
  int p1, p2, si, s, r, t_s2, t_r2, t_w2, n,n2;
  int a, b, l1, l2;
  int score;
  int **pos;
  int **list;
  int max_n;


  
  int nused;
  int *used_list, **used;
  int *sl2, **inv_pos;

  int filter=10;
  

  
  if ( !A) return 0;
  list=list_in[0];
  n=n_in[0];
  max_n=read_array_size (list, sizeof (int*));
  
 
  n2=0;
  pos=aln2pos_simple ( A,-1, ns, ls);
  inv_pos=vcalloc ((CL->S)->nseq, sizeof (int*));
  for (a=0; a<ns[1]; a++)inv_pos[ls[1][a]] =seq2inv_pos(A->seq_al[ls[1][a]]);
  
  l1=strlen (A->seq_al[ls[0][0]]);
  l2=strlen (A->seq_al[ls[1][0]]);
  sl2=vcalloc ((CL->S)->nseq, sizeof (int));
  
  for (a=0;a<ns[1]; a++)sl2[ls[1][a]]=1;
 

  
  used=declare_int (l2+1,2);
  used_list=vcalloc (l2+1, sizeof (int));
  nused=0;
    
  for (p1=0; p1<=l1; p1++)
    {
      for (nused=0,si=0;p1>0 && si<ns[0]; si++)
	{
	  s=ls [0][si];r=pos[s][p1-1];
	  for (a=1; r>0 && a<CL->residue_index[s][r][0];a+=ICHUNK)
	    {
	      t_s2=CL->residue_index[s][r][a+SEQ2];
	      t_r2=CL->residue_index[s][r][a+R2];
	      t_w2=CL->residue_index[s][r][a+WE];
	      
	      if (sl2[t_s2])
		{
		  p2=inv_pos[t_s2][t_r2];
		  score=t_w2;
		  if (!used[p2][1] && score>0)
		    {
		      used_list[nused++]=p2;
		    }
		  used[p2][0]+=score;
		  used[p2][1]++;
		}
	    }
	}
      for (a=0; a<nused; a++)
	{
	  
	  p2=used_list[a];
	 
	  score=used[p2][0]*SCORE_K;
	  used[p2][0]=used[p2][1]=0;
	  	  
	  if (n==max_n){max_n+=1000;list=vrealloc (list, max_n*sizeof (int*));}
	  if (score>filter && p1!=0 && p2!=0 && p1!=l1 && p2!=l2)
	    {
	      if (!list[n])list[n]=vcalloc (7, sizeof (int));
	      
	      list[n][0]=p1;
	      list[n][1]=p2;
	      list[n][3]=(l1-(p1))+(p2);
	      list[n][2]=score;
	      n++;
	    }
	}
    }
  

  vfree (used);
  vfree (used_list);
  free_int (inv_pos, -1);
  free_int (pos, -1);
  vfree (sl2);

  n_in[0]=n;
  list_in[0]=list;
  return 1;
}

int cl2pair_list_ecl_raw ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in)
{
  int p1, p2, si, s, r, t_s, t_r,t_w, t_s2, t_r2, t_w2, n,tot;
  int a, b, l1, l2;
  int **pos,**list;
  int max_n;
  
  
  int set, raw_max,nscore, score, nused;
  int *used_list, **used;
  int *sl2, **inv_pos;

  int filter1=0, filter2=0, max=0;
  int **nr;
  long tot_score=0, avg;
  int new_n=0;
  
  if ( !A) return 0;
  list=list_in[0];
  n=n_in[0];
  max_n=read_array_size (list, sizeof (int*));
  

  pos=aln2pos_simple ( A,-1, ns, ls);
  inv_pos=vcalloc ((CL->S)->nseq, sizeof (int*));
  for (a=0; a<ns[1]; a++)inv_pos[ls[1][a]] =seq2inv_pos(A->seq_al[ls[1][a]]);
  
  l1=strlen (A->seq_al[ls[0][0]]);
  l2=strlen (A->seq_al[ls[1][0]]);
  sl2=vcalloc ((CL->S)->nseq, sizeof (int));
  nr=declare_int (2, MAX(l1,l2)+1);
  
  for (a=0; a<l1;a++)
    for (b=0; b<ns[0]; b++)
      if (!is_gap(A->seq_al[ls[0][b]][a]))nr[0][a+1]++;
  for (a=0; a<l2;a++)
    for (b=0; b<ns[1]; b++)
      if (!is_gap(A->seq_al[ls[1][b]][a]))nr[1][a+1]++;
  
  for (a=0;a<ns[1]; a++)sl2[ls[1][a]]=1;
  

  
  used=declare_int (l2+1,2);
  used_list=vcalloc (l2+1, sizeof (int));
  nused=0;
    
  for (raw_max=0,p1=0; p1<=l1; p1++)
    {
      for (tot=0,nused=0,si=0;p1>0 && si<ns[0]; si++)
	{
	  s=ls [0][si];r=pos[s][p1-1];
	  for (a=1; r>0 && a<CL->residue_index[s][r][0];a+=ICHUNK)
	    {
	      t_s=CL->residue_index[s][r][a+SEQ2];
	      t_r=CL->residue_index[s][r][a+R2];
	      t_w=CL->residue_index[s][r][a+WE];
	      for (b=0; b<CL->residue_index[t_s][t_r][0];)
		{
		  if (b==0){t_s2=t_s;t_r2=t_r;t_w2=t_w;b++;}
		  else 
		    { 
		      t_s2=CL->residue_index[t_s][t_r][b+SEQ2];
		      t_r2=CL->residue_index[t_s][t_r][b+R2];
		      t_w2=CL->residue_index[t_s][t_r][b+WE];
		      b+=ICHUNK;
		    }
		  
		  if (sl2[t_s2])
		    {
		      p2=inv_pos[t_s2][t_r2];
		      score=MIN(t_w,t_w2);
		      if (score<filter1)score=0;
		      //if ( score<3) continue;
		      if (!used[p2][1] && score>0)
			{
			  used_list[nused++]=p2;
			}
		      tot+=score;
		      used[p2][0]+=score;
		      used[p2][1]++;
		    }
		}
	    }
	}
      //set the threshold to 1/2 of the best normalised score
      
      for (filter2=0,set=0,a=0; a<nused; a++)
	{
	  
	  p2=used_list[a];
	  
	  score=used[p2][0];
	  nscore=(score*100)/tot;
	  if (set==0){filter2=nscore;set=1;}
	  filter2=MAX(nscore,filter2);
	  if ( score<0)HERE ("*********** %d", score);
	}
      filter2/=2;
      max+=nused;
      for (a=0; a<nused; a++)
	{
	  
	  p2=used_list[a];
	 
	  //score=used[p2][0];
	  nscore=(used[p2][0]*100)/tot; //Normalized score used for filtering
	  score =used[p2][0];
	  raw_max=MAX(score, raw_max);
	  
	  used[p2][0]=used[p2][1]=0;
	  
	  if (n==max_n){max_n+=1000;list=vrealloc (list, max_n*sizeof (int*));}
	  if (nscore>=filter2 && p1!=0 && p2!=0 && p1!=l1 && p2!=l2)
	    {
	      if (!list[n])list[n]=vcalloc (7, sizeof (int));
	      
	      list[n][0]=p1;
	      list[n][1]=p2;
	      list[n][3]=(l1-(p1))+(p2);
	      list[n][2]=score;
	      n++;
	      tot_score+=score;
	      new_n++;
	    }
	}
    }
  avg=tot_score/new_n;
  
  //CL->gop=-1*avg*3;CL->gep=0;
  HERE ("FILTER: %d->%d [THR=%d]", max, n-n_in[0], filter2);
 
  vfree (used);
  vfree (used_list);
  free_int (inv_pos, -1);
  free_int (pos, -1);
  vfree (sl2);

  n_in[0]=n;
  list_in[0]=list;
  return 1;
}


/**
 * Calculates scores for diagonal segments.
 * 
 * \param Alignment The sequences.
 * \param ns Number of sequences in each group
 * \param ls sequences in in groups (ls[0][x] sequences in group 1, ls[1][x] squences in group 2).
 * \param CL the constraint list
 * \param list_in the diagonals
 * \param n_in number of sequences?
 */
int cl2pair_list_ecl_pc ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in)
{
  int p1, p2, si, s, r, t_s, t_r,t_w, t_s2, t_r2, t_w2, n;
  int a, b, l1, l2;
  int **pos,**list;
  int max_n;

  int nused;
  int *used_list;
  int *sl2, **inv_pos;

  int **nr;


  float nscore, score, tot, filter, avg=0, new=0;
  float **used;

 
  if ( !A) return 0;
  list=list_in[0];
  n=n_in[0];
  max_n=read_array_size (list, sizeof (int*));

  pos=aln2pos_simple ( A,-1, ns, ls);
  inv_pos=vcalloc ((CL->S)->nseq, sizeof (int*));
  for (a=0; a<ns[1]; a++)inv_pos[ls[1][a]] =seq2inv_pos(A->seq_al[ls[1][a]]);

  l1=strlen (A->seq_al[ls[0][0]]);
  l2=strlen (A->seq_al[ls[1][0]]);
  sl2=vcalloc ((CL->S)->nseq, sizeof (int));
  nr=declare_int (2, MAX(l1,l2)+1);

  for (a=0; a<l1;a++)
    for (b=0; b<ns[0]; b++)
      if (!is_gap(A->seq_al[ls[0][b]][a]))nr[0][a+1]++;
  for (a=0; a<l2;a++)
    for (b=0; b<ns[1]; b++)
      if (!is_gap(A->seq_al[ls[1][b]][a]))nr[1][a+1]++;
  
  for (a=0;a<ns[1]; a++)sl2[ls[1][a]]=1;
  
    
  used=declare_float (l2+1,2);
  used_list=vcalloc (l2+1, sizeof (int));
  nused=0;
  
  for (p1=0; p1<=l1; p1++)
    {
      
      for (tot=0,nused=0,si=0;p1>0 && si<ns[0]; si++)
	{
	  s=ls [0][si];r=pos[s][p1-1];
	  for (a=1; r>0 && a<CL->residue_index[s][r][0];a+=ICHUNK)
	    {
	      t_s=CL->residue_index[s][r][a+SEQ2];
	      t_r=CL->residue_index[s][r][a+R2];
	      t_w=CL->residue_index[s][r][a+WE];
	      for (b=0; b<CL->residue_index[t_s][t_r][0];)
		{
		  if (b==0){t_s2=t_s;t_r2=t_r;t_w2=t_w;b++;}
		  else 
		    { 
		      t_s2=CL->residue_index[t_s][t_r][b+SEQ2];
		      t_r2=CL->residue_index[t_s][t_r][b+R2];
		      t_w2=CL->residue_index[t_s][t_r][b+WE];
		      b+=ICHUNK;
		    }
		  
		  if (sl2[t_s2])
		    {
		      p2=inv_pos[t_s2][t_r2];
		      //score=((float)t_w/(float)NORM_F)*((float)t_w2/(float)NORM_F);
		      score=MIN(((float)t_w/(float)NORM_F),((float)t_w2/(float)NORM_F));
		      
		      if (!used[p2][1] && score>0)
			{
			  used_list[nused++]=p2;
			}
		      
		      tot+=score;
		      used[p2][0]+=score;
		      used[p2][1]++;
		    }
		}
	    }
	}
      //FILTER: Keep in the graph the edges where (p1->p2/(Sum (P1->x))>0.01
      filter=0.01;
      
      for (a=0; a<nused; a++)
	{
	  
	  p2=used_list[a];
	  nscore=used[p2][0]/tot; //Normalized score used for filtering
	  score =used[p2][0];
		  
	  used[p2][0]=used[p2][1]=0;
	  if (n==max_n)
	    {
	      max_n+=10000;list=vrealloc (list, max_n*sizeof (int*));
	    }
	  if (nscore>filter && p1!=0 && p2!=0 && p1!=l1 && p2!=l2)
	    {
	      if (!list[n])
			  list[n]=vcalloc (7, sizeof (int));
	      
	      list[n][0]=p1;
	      list[n][1]=p2;
	      list[n][3]=(l1-(p1))+(p2);
	      score/=(float)((CL->S)->nseq*nr[0][p1]*nr[1][p2]);
	      list[n][2]=(int)((float)score*(float)NORM_F);
	      avg+=(int)((float)score*(float)NORM_F);
	      new++;
	      n++;
	    }
	}
    }
  free_float (used, -1);
  vfree (used_list);
  free_int (inv_pos, -1);
  free_int (pos, -1);
  vfree (sl2);
  free_int (nr, -1);

  n_in[0]=n;
  list_in[0]=list;
  if (new)avg/=new;
  return avg;
}


int cl2pair_list_ecl_rawquad ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in)
{
  int p1, p2, si, s, r, t_s, t_r,t_w, t_s2, t_r2, t_w2,t_s3, t_r3, t_w3, n,n2;
  int a, b, c,l1, l2;
  int score;
  int **pos;
  int **list;
  int max_n;

  int tn;
  
  int nused;
  int *used_list, **used;
  int *sl2, **inv_pos;

  int filter=0;
  int nseq2;
  if ( !A) return 0;
  list=list_in[0];
  n=n_in[0];
  max_n=read_array_size (list, sizeof (int*));
  
 
  n2=0;
  pos=aln2pos_simple ( A,-1, ns, ls);
  inv_pos=vcalloc ((CL->S)->nseq, sizeof (int*));
  for (a=0; a<ns[1]; a++)inv_pos[ls[1][a]] =seq2inv_pos(A->seq_al[ls[1][a]]);
  
  l1=strlen (A->seq_al[ls[0][0]]);
  l2=strlen (A->seq_al[ls[1][0]]);
  sl2=vcalloc ((CL->S)->nseq, sizeof (int));
  
  for (a=0;a<ns[1]; a++)sl2[ls[1][a]]=1;
 
  
  used=declare_int (l2+1,2);
  used_list=vcalloc (l2+1, sizeof (int));
  nused=0;
  nseq2=(CL->S)->nseq*(CL->S)->nseq;
  
  for (p1=0; p1<=l1; p1++)
    {
      for (nused=0,si=0;p1>0 && si<ns[0]; si++)
	{
	  s=ls [0][si];r=pos[s][p1-1];
	  for (a=1; r>0 && a<CL->residue_index[s][r][0];a+=ICHUNK)
	    {
	      t_s=CL->residue_index[s][r][a+SEQ2];
	      t_r=CL->residue_index[s][r][a+R2];
	      t_w=CL->residue_index[s][r][a+WE];
	      for (b=0; b<CL->residue_index[t_s][t_r][0];)
		{
		  if (b==0){t_s2=t_s;t_r2=t_r;t_w2=t_w;b++;}
		  else 
		    { 
		      t_s2=CL->residue_index[t_s][t_r][b+SEQ2];
		      t_r2=CL->residue_index[t_s][t_r][b+R2];
		      t_w2=CL->residue_index[t_s][t_r][b+WE];
		      b+=ICHUNK;
		    }
		  if (sl2[t_s2])
		    {
		      for (c=0; c<CL->residue_index[t_s2][t_r2][0];)
			{
			  if (c==0){t_s3=t_s2;t_r3=t_r2;t_w3=t_w2;c++;}
			  else 
			    { 
			      t_s3=CL->residue_index[t_s2][t_r2][c+SEQ2];
			      t_r3=CL->residue_index[t_s2][t_r2][c+R2];
			      t_w3=CL->residue_index[t_s2][t_r2][c+WE];
			      c+=ICHUNK;
			    }
			  
			  if (sl2[t_s3])
			    {
			      p2=inv_pos[t_s3][t_r3];
			      score=MIN(t_w,t_w2);
			      score=MIN(score,t_w3);
			      if (!used[p2][1] && score>0)
				{
				  used_list[nused++]=p2;
				}
			      used[p2][0]+=score;
			      used[p2][1]++;
			    }
			}
		    }
		}
	    }
	}
      for (a=0; a<nused; a++)
	{
	  
	  p2=used_list[a];
	  score=used[p2][0];
	  
	  used[p2][0]=used[p2][1]=0;
	  	  
	  if (n==max_n){max_n+=1000;list=vrealloc (list, max_n*sizeof (int*));}
	  if (score>filter && p1!=0 && p2!=0 && p1!=l1 && p2!=l2)
	    {
	      if (!list[n])list[n]=vcalloc (7, sizeof (int));
	      
	      list[n][0]=p1;
	      list[n][1]=p2;
	      list[n][3]=(l1-(p1))+(p2);
	      list[n][2]=score;
	      n++;
	    }
	}
    }
  

  vfree (used);
  vfree (used_list);
  free_int (inv_pos, -1);
  free_int (pos, -1);
  vfree (sl2);

  n_in[0]=n;
  list_in[0]=list;
  return 1;
}
int cl2pair_list_ecl_norm ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in)
{
  int p1, p2, si, s, r, t_s, t_r, n,n2;
  int a, b, l1, l2;
  int score;
  int **pos;
  int **list;
  int max_n;


  
  int nused;
  int *used_list, *used;
  int *sl2, **inv_pos;

  int filter=0;
  
  if ( !A) return 0;
  list=list_in[0];
  n=n_in[0];
  max_n=read_array_size (list, sizeof (int*));
  
 
  n2=0;
  pos=aln2pos_simple ( A,-1, ns, ls);
  inv_pos=vcalloc ((CL->S)->nseq, sizeof (int*));
  for (a=0; a<ns[1]; a++)inv_pos[ls[1][a]] =seq2inv_pos(A->seq_al[ls[1][a]]);
  
  l1=strlen (A->seq_al[ls[0][0]]);
  l2=strlen (A->seq_al[ls[1][0]]);
  sl2=vcalloc ((CL->S)->nseq, sizeof (int));
  
  for (a=0;a<ns[1]; a++)sl2[ls[1][a]]=1;
 
  
  used=vcalloc (l2+1, sizeof (int));
  used_list=vcalloc (l2+1, sizeof (int));
  nused=0;
  
  
  
  for (p1=0; p1<=l1; p1++)
    {
      for (si=0;p1>0 && si<ns[0]; si++)
	{
	  s=ls [0][si];
	  r=pos[s][p1-1];
	  
	  for (a=1; r>0 && a<CL->residue_index[s][r][0];a+=ICHUNK)
	    {
	      
	      t_s=CL->residue_index[s][r][a+SEQ2];
	      t_r=CL->residue_index[s][r][a+R2];
	      
	      for (b=0; b<CL->residue_index[t_s][t_r][0];)
		{
		  int t_s2, t_r2;
		  if (b==0){t_s2=t_s;t_r2=t_r;b++;}
		  else { 
		    t_s2=CL->residue_index[t_s][t_r][b+SEQ2];
		    t_r2=CL->residue_index[t_s][t_r][b+R2];
		    b+=ICHUNK;
		  }
		  
		  if (sl2[t_s2])
		    {
		      p2=inv_pos[t_s2][t_r2];
		      if (!used[p2]){used[p2]=1;used_list[nused++]=p2;}
		      else used[p2]++;
		    }
		}
	      
	    }
	}
      if (p1==0 || p1==l1)
	{
	  for (nused=0,p2=0; p2<=l2; p2++)used_list[nused++]=p2;
	}
      else
	{
	  if (!used[0])used_list[nused++]=0;
	  if (!used[l2])used_list[nused++]=l2;
	}
      for (a=0; a<nused; a++)
	{
	  p2=used_list[a];
	  if (p2==0 || p1==0)score=0;
	  
	  else score=(CL->get_dp_cost) (A, pos, ns[0], ls[0], p1-1, pos, ns[1], ls[1],p2-1,CL);	
	  if (score>filter && p1!=0 && p2!=0 && p1!=l1 && p2!=l2)
	    {
	      if (n==max_n){max_n+=1000;list=vrealloc (list, max_n*sizeof (int*));}
	      if (!list[n])list[n]=vcalloc (7, sizeof (int));
	      list[n][0]=p1;
	      list[n][1]=p2;
	      list[n][3]=(l1-(p1))+(p2);
	      list[n][2]=score;
	      n++;
	      if (p1!=0 && p2!=0 && p1!=l1 && p2!=l2)n2++;
	    }
	  used[p2]=0;
	}
      
      nused=0;
    }
  
  vfree (used);
  vfree (used_list);
  free_int (inv_pos, -1);
  free_int (pos, -1);
  vfree (sl2);
  n_in[0]=n;
  list_in[0]=list;
  

   return 1;
}




int cl2pair_list_ref( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list, int *n)
{
  int a, b, l1, l2, n2=0;
  int score;
  int **pos;
  int max_n;
  if (!A) return 0;
  pos=aln2pos_simple ( A,-1, ns, ls);
  l1=strlen (A->seq_al[ls[0][0]]);
  l2=strlen (A->seq_al[ls[1][0]]);
  
  max_n=read_array_size (list[0], sizeof (int));
  
  for (a=0; a<=l1; a++)
    for (b=0; b<=l2; b++)
      {
	score=(a==0 || b==0)?0:slow_get_dp_cost_pc(A, pos, ns[0], ls[0],a-1, pos, ns[1], ls[1], b-1, CL);

	
	if ( score>0 && a!=0 && b!=0 && a!=l1 && b!=l2)
	  {
	    if (n[0]==max_n){max_n+=1000;list[0]=vrealloc (list[0], max_n*sizeof (int*));}
	    if (!list[0][n[0]])list[0][n[0]]=vcalloc (7, sizeof (int));
	    list[0][n[0]][0]=a;
	    list[0][n[0]][1]=b;
	    list[0][n[0]][3]=(l1-a)+b;
	    list[0][n[0]][2]=score;
	    if ( a!=0 && b!=0 && a!=l1 && b!=l2)
	      {
		n2++;
	      }
	    n[0]++;
	  }

      }

  return n[0];
  }

/**
 * Calculates scores for diagonal segments.
 * 
 * \param Alignment The sequences.
 * \param ns Number of sequences in each group
 * \param ls sequences in in groups (ls[0][x] sequences in group 1, ls[1][x] squences in group 2).
 * \param CL the constraint list
 * \param list_in the diagonals
 * \param n_in number of sequences?
 */
int cl2pair_list_ecl_ext_pc ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in)
{
  int **list, **pos, **cache, **reallocT;
  int  *sl2;
  int si,a,b,n, l1, l2, p1, p2,s,r;
  int t_s , t_r , t_w;
  int t_s2, t_r2, t_w2;
  int max_len,nseq, avg;
  int tot=0;
  
  
  if ( !A) return 0;
  list=list_in[0];
  n=n_in[0];
  
  pos=aln2pos_simple ( A,-1, ns, ls);
  l1=strlen (A->seq_al[ls[0][0]]);
  l2=strlen (A->seq_al[ls[1][0]]);
  sl2=vcalloc ((CL->S)->nseq, sizeof (int));
  for (a=0;a<ns[1]; a++)sl2[ls[1][a]]=1;
  
  max_len=(CL->S)->max_len;
  nseq=(CL->S)->nseq;
  
  cache=declare_int ( nseq+1, max_len+1);
  reallocT=declare_int (nseq+1, max_len+1);
  
  for (p1=1; p1<=l1; p1++)
    {
      for (si=0;si<ns[0]; si++)
	{
	  s=ls[0][si];r=pos[s][p1-1];
	  s=name_is_in_list (A->name[s], (CL->S)->name, (CL->S)->nseq, 100);
	  if (r<=0)continue;
	  cache[s][r]=1;
	  for (a=1;a<CL->residue_index[s][r][0];a+=ICHUNK)
	    {
	      
	      
	      t_s=CL->residue_index[s][r][a+SEQ2];
	      t_r=CL->residue_index[s][r][a+R2];
	      t_w=CL->residue_index[s][r][a+WE];
	      
	      cache[t_s][t_r]=1;
	    }
	  for (a=1;a<CL->residue_index[s][r][0];a+=ICHUNK)
	    {
	      
	      
	      t_s=CL->residue_index[s][r][a+SEQ2];
	      t_r=CL->residue_index[s][r][a+R2];
	      t_w=CL->residue_index[s][r][a+WE];
	      
	      for (b=1; b<CL->residue_index[t_s][t_r][0];b+=ICHUNK)
		{
		  t_s2=CL->residue_index[t_s][t_r][b+SEQ2];
		  t_r2=CL->residue_index[t_s][t_r][b+R2];
		  t_w2=CL->residue_index[t_s][t_r][b+WE];
		  
		  if (!cache[t_s2][t_r2])
		    {
		      tot++;
		      int ni;
		      cache[t_s2][t_r2]=1;
		      if (!reallocT[s][r])reallocT[s][r]=CL->residue_index[s][r][0];
		      ni=CL->residue_index[s][r][0];
		      CL->residue_index[s][r][0]+=ICHUNK;
		      
		      CL->residue_index[s][r]=vrealloc ( CL->residue_index[s][r], CL->residue_index[s][r][0]*sizeof (int));
		      CL->residue_index[s][r][ni+SEQ2]=t_s2;
		      CL->residue_index[s][r][ni+R2]=t_r2;
		      CL->residue_index[s][r][ni+WE]=MIN(t_w2,t_w);
		      
		      if (!reallocT[t_s2][t_r2])reallocT[t_s2][t_r2]=CL->residue_index[t_s2][t_r2][0];
		      ni=CL->residue_index[t_s2][t_r2][0];
		      CL->residue_index[t_s2][t_r2][0]+=ICHUNK;
		      CL->residue_index[t_s2][t_r2]=vrealloc ( CL->residue_index[t_s2][t_r2], CL->residue_index[t_s2][t_r2][0]*sizeof (int));
		      CL->residue_index[t_s2][t_r2][ni+SEQ2]=s;
		      CL->residue_index[t_s2][t_r2][ni+R2]=r;
		      CL->residue_index[t_s2][t_r2][ni+WE]=MIN(t_w2,t_w);
		      
		    }
		}
	      
	    }

	  cache[s][r]=0;
	  for (a=1;a<CL->residue_index[s][r][0];a+=ICHUNK)
	    {
	      
	      
	      t_s=CL->residue_index[s][r][a+SEQ2];
	      t_r=CL->residue_index[s][r][a+R2];
	      t_w=CL->residue_index[s][r][a+WE];
	      
	      cache[t_s][t_r]=0;
	    }
	}
    }

  for (p1=1; p1<=l2; p1++)
    {
      for (si=0;si<ns[1]; si++)
	{
	  s=ls[1][si];r=pos[s][p1-1];
	  s=name_is_in_list (A->name[s], (CL->S)->name, (CL->S)->nseq, 100);
	  if (r<=0)continue;
	  cache[s][r]=1;
	  for (a=1;a<CL->residue_index[s][r][0];a+=ICHUNK)
	    {
	      
	      
	      t_s=CL->residue_index[s][r][a+SEQ2];
	      t_r=CL->residue_index[s][r][a+R2];
	      t_w=CL->residue_index[s][r][a+WE];
	      
	      cache[t_s][t_r]=1;
	    }
	  for (a=1;a<CL->residue_index[s][r][0];a+=ICHUNK)
	    {
	      
	      
	      t_s=CL->residue_index[s][r][a+SEQ2];
	      t_r=CL->residue_index[s][r][a+R2];
	      t_w=CL->residue_index[s][r][a+WE];
	      
	      for (b=1; b<CL->residue_index[t_s][t_r][0];b+=ICHUNK)
		{
		  t_s2=CL->residue_index[t_s][t_r][b+SEQ2];
		  t_r2=CL->residue_index[t_s][t_r][b+R2];
		  t_w2=CL->residue_index[t_s][t_r][b+WE];
		  
		  if (!cache[t_s2][t_r2])
		    {
		      tot++;
		      int ni;
		      
		      cache[t_s2][t_r2]=1;
		      if (!reallocT[s][r])reallocT[s][r]=CL->residue_index[s][r][0];
		      ni=CL->residue_index[s][r][0];
		      CL->residue_index[s][r][0]+=ICHUNK;
		      CL->residue_index[s][r]=vrealloc ( CL->residue_index[s][r], CL->residue_index[s][r][0]*sizeof (int));
		      CL->residue_index[s][r][ni+SEQ2]=t_s2;
		      CL->residue_index[s][r][ni +R2]=t_r2;
		      CL->residue_index[s][r][ni+WE]=MIN(t_w2,t_w);
		      
		      if (!reallocT[t_s2][t_r2])reallocT[t_s2][t_r2]=CL->residue_index[t_s2][t_r2][0];
		      ni=CL->residue_index[t_s2][t_r2][0];
		      CL->residue_index[t_s2][t_r2][0]+=ICHUNK;
		      CL->residue_index[t_s2][t_r2]=vrealloc ( CL->residue_index[t_s2][t_r2], CL->residue_index[t_s2][t_r2][0]*sizeof (int));
		      CL->residue_index[t_s2][t_r2][ni+SEQ2]=s;
		      CL->residue_index[t_s2][t_r2][ni+R2]=r;
		      CL->residue_index[t_s2][t_r2][ni+WE]=MIN(t_w2,t_w);
		      
		    }
		}
	      
	    }

	  cache[s][r]=0;
	  for (a=1;a<CL->residue_index[s][r][0];a+=ICHUNK)
	    {
	      
	      
	      t_s=CL->residue_index[s][r][a+SEQ2];
	      t_r=CL->residue_index[s][r][a+R2];
	      t_w=CL->residue_index[s][r][a+WE];
	      
	      cache[t_s][t_r]=0;
	    }
	}
    }

  

  HERE ("TOT=%d", tot);
  avg=cl2pair_list_ecl_pc (A, ns,ls, CL,list_in,n_in);
  
  for (s=0; s<=nseq; s++)
    for (r=0; r<=max_len; r++)
      if ( reallocT[s][r])CL->residue_index[s][r][0]=reallocT[s][r];
  free_int (reallocT, -1);
  free_int (cache, -1);
  return avg;
}

int list2linked_pair_wise ( Alignment *A, int *ns, int **l_s, Constraint_list *CL, int **list, int n);
int two_pass_linked_pair_wise ( Alignment *A, int *ns, int **l_s, Constraint_list *CL)
{
  int n=0, **list=NULL;
  int nscore;
  int mode=2;
  int id;
  
  cl2pair_list (A,ns, l_s, CL, &list, &n, mode, 10);
  nscore=list2linked_pair_wise (A, ns, l_s, CL, list, n);
  id=sub_aln2sim (A, ns, l_s, "idmat_sim");
  
  if (id>50)return nscore;
  ungap_sub_aln ( A, ns[0], l_s[0]);
  ungap_sub_aln ( A, ns[1], l_s[1]);
  cl2pair_list (A,ns, l_s, CL, &list, &n,mode,0);
  nscore=list2linked_pair_wise (A, ns, l_s, CL, list, n);
  cl2pair_list (NULL,ns, l_s, CL, &list, &n, mode, 0);
  return nscore;
}
int clinked_pair_wise ( Alignment *A, int *ns, int **l_s, Constraint_list *CL)
{
  int n=0, **list=NULL;
  int nscore, pscore=0;
  int mode=2;

  int add=0;
  cl2pair_list (A,ns, l_s, CL, &list, &n, mode, 1000);
  nscore=list2linked_pair_wise (A, ns, l_s, CL, list, n);
  HERE ("***********First: %d", nscore);
  if (add)
    {
      while (nscore>pscore)
	{
	  pscore=nscore;
	  ungap_sub_aln ( A, ns[0], l_s[0]);
	  ungap_sub_aln ( A, ns[1], l_s[1]);
	  cl2pair_list (A,ns, l_s, CL, &list, &n, mode, 10);
	  nscore=list2linked_pair_wise (A, ns, l_s, CL, list, n);
	  HERE ("****************New: %d", nscore);
	}
    }
  cl2pair_list (NULL,ns, l_s, CL, &list, &n, mode, 0);
  return nscore;
}
int linked_pair_wise ( Alignment *A, int *nsi, int **lsi, Constraint_list *CL)
{
  int n=0;
  static int **list=NULL;
  int score, a;
  int *ns, **ls;
  int mode=1;//1:ecl, 0:ref

  ns=vcalloc (2, sizeof (int));
  ns[0]=nsi[1]; ns[1]=nsi[0];
  
  ls=declare_int (2, ns[0]+ns[1]);
  for (a=0; a<ns[1]; a++)
    ls[1][a]=lsi[0][a];
  for (a=0; a<ns[0]; a++)
    ls[0][a]=lsi[1][a];
  
  
  cl2pair_list (A,ns, ls, CL, &list, &n, mode, 0);
    
  score=list2linked_pair_wise (A, ns, ls, CL, list, n);
  cl2pair_list (NULL,ns, ls, CL, &list, &n, mode, 0);
 
  return score;
}
 

int list2linked_pair_wise_nolgp ( Alignment *A, int *ns, int **l_s, Constraint_list *CL, int **list, int n);
int list2linked_pair_wise ( Alignment *A, int *ns, int **l_s, Constraint_list *CL, int **list, int n)
 {
   int mode=0;
   if (mode==0)return list2linked_pair_wise_nolgp    (A, ns, l_s, CL, list, n);
   return 0;
 }


#define LIN(a,b,c) a[b*5+c]
int list2linked_pair_wise_nolgp( Alignment *A, int *ns, int **l_s, Constraint_list *CL, int **list, int n)
{
  int a,b,c, i, j, LEN=0, start_trace;
  int pi, pj,ij, delta_i, delta_j, prev_i, prev_j;
  static int **slist;
  static long *MI, *MJ, *MM,*MT2;
  static int *sortseq;
  static int max_size;
  int gop, gep, igop, igep;
  int l1, l2, l, ls;
  char **al;
  char **aln,*char_buf;
  int ni=0, nj=0;
  long score;
  int nomatch;
 
  l1=strlen (A->seq_al[l_s[0][0]]);
  l2=strlen (A->seq_al[l_s[1][0]]);
  al=declare_char (2,l1+l2+1);
  
 
  //Penalties: max score is NORM_F
  //Penalties must be negative
  igop=CL->gop;
  gep=igep=CL->gep;
  
  if (n>max_size)
    {
      max_size=n;
      
      vfree (MI);vfree (MJ); vfree (MM);
      free_int (slist, -1);
     
      slist=declare_int (n,3);
      
      MI=vcalloc (5*n, sizeof (long));
      MJ=vcalloc (5*n, sizeof (long));
      MM=vcalloc (5*n, sizeof (long));
      
    }
  else
    {
      for (a=0; a<n; a++)
	for (b=0; b<5; b++)LIN(MI,a,b)=LIN(MJ,a,b)=LIN(MJ,a,b)=-1000000;
    }
  
  /*New Bit: Start*/
  if (!sortseq) sortseq=vcalloc( 7, sizeof (int));
  sortseq[0]=0; sortseq[1]=1;sortseq[2]=-1;
  sort_list_int2 (list, sortseq,7, 0, n-1);
  for (a=0; a<n; a++)
  {
      slist[a][0]=a;
      list[a][4]=a;
  }
  
  sortseq[0]=1; sortseq[1]=0;sortseq[2]=-1;
  sort_list_int2 (list, sortseq,7, 0, n-1);
  for (a=0; a<n; a++)
  {
	slist[a][1]=list[a][4];
	list[a][5]=a;
  }
 
  sortseq[0]=3; sortseq[1]=0;sortseq[2]=1;sortseq[3]=-1;
  sort_list_int2 (list, sortseq,7, 0, n-1);
  for (a=0; a<n; a++)
    {
      slist[a][2]=list[a][4];
      list[a][6]=a;
    }

  sortseq[0]=0; sortseq[1]=1;sortseq[2]=-1;
  sort_list_int2 (list, sortseq,7, 0, n-1);
 
  /*New Bit: EnD*/
  
  
 
 
  
  
  for (a=0; a<n; a++)
    {

      
      i=list[a][0];
      j=list[a][1];

      
      if (i==l1 || j==l2)gop=0;
      else gop=igop;

      if (i==l1 && j==l2)start_trace=a;
      else if ( i==0 || j==0)
	{
	  LIN(MM,a,0)=-1000000;
	  if (j==0)
	    {
	      
	      LIN(MJ,a,0)=-10000000;
	      LIN(MI,a,0)=gep*i;
	      
	    }
	  else if (i==0)
	    {
	      
	      LIN(MI,a,0)=-10000000;
	      LIN(MJ,a,0)=gep*j;
	      
	    }

	  LIN(MI,a,1)=LIN(MJ,a,1)=-1;
	  LIN(MI,a,2)=LIN(MJ,a,2)=i;
	  LIN(MI,a,3)=LIN(MJ,a,3)=j;
	  continue;
	}
      
      pi=list[a][5];
      pi=slist[pi-1][1];
      
      pj=list[a][4];
      pj=slist[pj-1][0]; 
      
      ij=list[a][6];
      ij=slist[ij-1][2];
      
      
      ij=list[a][6];
      ij=slist[ij-1][2];
      
      
	
     
      
      prev_i=list[pi][0];
      prev_j=list[pj][1];
      
      delta_i=list[a][0]-list[pi][0];
      delta_j=list[a][1]-list[pj][1];
      
      /*Linear Notation*/
      LIN(MI,a,0)=MAX(LIN(MI,pi,0),(LIN(MM,pi,0)+gop))+delta_i*gep;
      LIN(MI,a,1)=pi;
      LIN(MI,a,2)=delta_i;
      LIN(MI,a,3)=0;
      LIN(MI,a,4)=(LIN(MI,pi,0)>=(LIN(MM,pi,0)+gop))?'i':'m';
    
      
      LIN(MJ,a,0)=MAX(LIN(MJ,pj,0),(LIN(MM,pj,0)+gop))+delta_j*gep;
      LIN(MJ,a,1)=pj;
      LIN(MJ,a,2)=0;
      LIN(MJ,a,3)=delta_j;
      
      LIN(MJ,a,4)=(LIN(MJ,pj,0)>=LIN(MM,pj,0)+gop)?'j':'m';
      
    
      if (a>1 && (ls=list[a][0]-list[ij][0])==(list[a][1]-list[ij][1]))
	{
	  LIN(MM,a,0)=MAX3(LIN(MM,ij,0),LIN(MI,ij,0),LIN(MJ,ij,0))+list[a][2]-(ls*CL->nomatch);

	  LIN(MM,a,1)=ij;
	  LIN(MM,a,2)=ls;
	  LIN(MM,a,3)=ls;
	  if ( LIN(MM,ij,0)>=LIN(MI,ij,0) && LIN(MM,ij,0)>=LIN(MJ,ij,0))LIN(MM,a,4)='m';
	  else if ( LIN(MI,ij,0) >= LIN(MJ,ij,0))LIN(MM,a,4)='i';
	  else LIN(MM,a,4)='j';
	  
	}
      else
	{
	  LIN(MM,a,0)=UNDEFINED;
	  LIN(MM,a,1)=-1;
	}  
    }
  
  a=start_trace;
  if (LIN(MM,a,0)>=LIN(MI,a,0) && LIN(MM,a,0) >=LIN(MJ,a,0))MT2=MM;
  else if ( LIN(MI,a,0)>=LIN(MJ,a,0))MT2=MI;
  else MT2=MJ;

  score=MAX3(LIN(MM,a,0), LIN(MI,a,0), LIN(MJ,a,0));
  
  i=l1;
  j=l2;
  
  
  while (!(i==0 &&j==0))
    {
      int next_a;
      l=MAX(LIN(MT2,a,2),LIN(MT2,a,3));
      // HERE ("%c from %c %d %d SCORE=%d [%d %d] [%2d %2d]", T2[a][5],T2[a][4], T2[a][2], T2[a][3], T2[a][0], gop, gep, i, j);
      if (i==0)
	{
	  while ( j>0)
	    {
	      al[0][LEN]=0;
	      al[1][LEN]=1;
	      j--; LEN++;
	    }
	}
      else if (j==0)
	{
	  while ( i>0)
	    {
	      al[0][LEN]=1;
	      al[1][LEN]=0;
	      i--; LEN++;
	    }
	}
      
      else if (l==0) {HERE ("L=0 i=%d j=%d",l, i, j);exit (0);}
      else 
	{
	  for (b=0; b<l; b++, LEN++)
	    {
	      if (LIN(MT2,a,2)){al[0][LEN]=1;i--;ni++;}
	      else al[0][LEN]=0;
	      
	      if (LIN(MT2,a,3)){al[1][LEN]=1;j--;nj++;}
	      else al[1][LEN]=0;
	    }
	  
	  next_a=LIN(MT2,a,1);
	  if (LIN(MT2,a,4)=='m')MT2=MM;
	  else if (LIN(MT2,a,4)=='i')MT2=MI;
	  else if (LIN(MT2,a,4)=='j')MT2=MJ;
	  a=next_a;
	}
    }
  
 
  
  invert_list_char ( al[0], LEN);
  invert_list_char ( al[1], LEN);
  
  	
  if ( A->declared_len<=LEN)A=realloc_aln  ( A,2*LEN+1);
  aln=A->seq_al;
  char_buf= vcalloc (LEN+1, sizeof (char));
	
  for ( c=0; c< 2; c++)
    {
      for ( a=0; a< ns[c]; a++) 
	{		
	  int ch=0;
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
  
  vfree (char_buf);
  free_char (al, -1);
  
  return score;
}
int ** aln2local_penalties4link (Alignment *A, int n, int *ls, Constraint_list *CL, int **lg);
int ** aln2local_penalties4link (Alignment *A, int n, int *ls, Constraint_list *CL, int **lg)
{
  //adapted from gap_count in MAFFT V 5.5
  int p,s,l, c1, c2;
  int gep,gop;
  int open=3, close=4, gap=5;
  
  gop=CL->gop;
  gep=CL->gep;
  
  l=strlen (A->seq_al[ls[0]]);
  
  if (!lg)
    {
      lg=declare_int (6, l);
    }
  
  if ( read_array_size_new (lg[0])<l)
    {
      free_int (lg, -1);
      lg=declare_int (6, l);
    }
  
  for( s=0; s<n; s++ ) 
	{
	  c1='x';
	  for (p=0; p<l; p++)
	    {
	      c2=A->seq_al[ls[s]][p];
	      
	      if (c1!='-' && c2=='-')lg[open][p]++;
	      if (c1=='-' && c2!='-')lg[close][p]++;
	      if ( c1=='-')lg[gap][p]++;
	      c1=c2;
	    }
	}
  
  for (p=0; p<l; p++)
    {
      float go, gc, nn;
      nn=n;
      go=lg[open ][p];
      gc=lg[close][p];
     
    
      lg[GOP][p]=0.5*(1-(go/nn))*gop;
      lg[GCP][p]=0.5*(1-(gc/nn))*gop;
      //Checked locacal gep => gives low quality results
      lg[GEP][p]=gep;//(1-((float)lg[gap][p]/(float)n))*gep;
      lg[open][p]=lg[close][p]=lg[gap][p]=0;
      
    }

  return lg;
}

