#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include "io_lib_header.h"
#include "util_lib_header.h"
#include "dp_lib_header.h"
#include "define_header.h"
/************************************************************************************/
/*                NEW      ANALYZE 2    : SAR                                        */
/************************************************************************************/
float display_prediction_old (int **prediction, int n, Alignment *A, Alignment *S, int field);

float display_prediction (int ***count, Alignment *S, int c, int n);
Alignment * filter_aln4sar0 ( Alignment *A, Alignment *S, int c, int leave, char *mode);
Alignment * filter_aln4sar1 ( Alignment *A, Alignment *S, int c, int leave, char *mode);
Alignment * filter_aln4sar2 ( Alignment *A, Alignment *S, int c, int leave, char *mode);
Alignment * filter_aln4sar3 ( Alignment *A, Alignment *S, int c, int leave, char *mode);
Alignment * filter_aln4sar4 ( Alignment *A, Alignment *S, int c, int leave, char *mode);
Alignment * filter_aln4sar5 ( Alignment *A, Alignment *S, int c, int leave, char *mode);

int **sar2profile ( Alignment *A, Alignment *S, int c, int leave);
int **sar2profile_sim ( Alignment *A, Alignment *S, int **sim, int comp, int leave);
int sar_profile2score ( char *seq, int **profile);
double sar_vs_iseq1( char *sar, int  *seq, float gl, int **sim, char *best_aa);
double sar_vs_seq1 ( char *sar, char *seq, float gl, int **sim, char *best_aa);
double sar_vs_seq2 ( char *sar, char *seq, float ng, int **mat, char *a);
double sar_vs_seq3 ( char *sar, char *seq, float ng, int **mat, char *a);

double sar_vs_iseq4 ( char *sar, int *seq, float ng, int **mat, char *a);//supports an extended alphabet
double sar_vs_seq4 ( char *sar, char *seq, float ng, int **mat, char *a);

double sar_vs_seq5 ( char *sar, char *seq, float ng, int **mat, char *a);
int make_sim_pred ( Alignment *A,Alignment *S, int comp, int seq);

int **sar2profile_sim ( Alignment *A, Alignment *S, int **sim, int comp, int leave)
{

  int a, b, r, c, c1, c2, r1, r2, s, p;
  int ***cache, **profile;


  profile=declare_int (A->len_aln, 26);
  cache=(int***)declare_arrayN (3,sizeof (int),2,A->len_aln, 26);

  for ( a=0; a< A->len_aln; a++)
    for ( b=0; b< A->nseq; b++)
      {
	r=tolower(A->seq_al[b][a]);
	c=( S->seq_al[comp][b]=='I')?1:0;
	if (b==leave || is_gap(r)) continue;
	cache [c][a][r-'a']++;
      }
  for (a=0; a< A->nseq; a++)
    {
      if ( a==leave) continue;
      for ( b=0; b< A->nseq; b++)
	{
	  c1=(S->seq_al[comp][a]=='I')?1:0;
	  c2=(S->seq_al[comp][b]=='I')?1:0;
	  if ( b==leave || b==a || c1!=1 || c1==c2) continue;
	  s=sim[a][b];

	  for (p=0; p<A->len_aln; p++)
	    {
	      r1=tolower(A->seq_al[a][p]);
	      r2=tolower(A->seq_al[b][p]);
	      if ( is_gap(r1) || is_gap(r2) || r1==r2)continue;
	      r1-='a';r2-='a';
	      if (cache[1][p][r2])continue;
	      if ( s<50)continue;
	      profile[p][r2]-=s;
	    }
	}
    }

  free_arrayN((void***)cache,3);
  return profile;

}
int **sar2profile ( Alignment *A, Alignment *S, int comp, int leave)
{

  int a, b,c,r, n, v, npos=0;
  int ***cache, **profile;
  int ncat;
  float n_gap, max_gap;
  profile=declare_int (A->len_aln, 26);
  cache=(int***)declare_arrayN (3,sizeof (int),2,A->len_aln, 26);



  for ( n=0, a=0; a< A->nseq; a++)
    {
      if ( a==leave) continue;
      else n+=(S->seq_al[comp][a]=='I')?1:0;
    }

  for ( a=0; a< A->len_aln; a++)
    for ( b=0; b< A->nseq; b++)
      {
	r=tolower(A->seq_al[b][a]);
	c=( S->seq_al[comp][b]=='I')?1:0;
	if (b==leave) continue;
	else if (is_gap(r))continue;
	r-='a';
	cache [c][a][r]++;
      }

  ncat=15; /*ncat: limit the analysis to columns containing less than ncat categories of aa*/
  max_gap=0.05;
  for (a=0; a< A->len_aln; a++)
    {
      for (n_gap=0,b=0; b< A->nseq; b++)
	n_gap+=(is_gap(A->seq_al[b][a]));
      n_gap/=(float)A->nseq;

      if ( n_gap> max_gap)continue;

      for (v=0,r=0; r< 26; r++)
	{
	  if (cache [0][a][r] || cache[1][a][r])v++;
	}

      for (n=0,r=0; r< 26 && v<ncat; r++)
	{
	  if (cache [0][a][r] && !cache[1][a][r])
	    {
	      n++;
	      profile[a][r]=-cache[0][a][r];
	    }
	}
      if (n) npos++;
    }

  free_arrayN((void***)cache,3);
  return profile;

}
Alignment * filter_aln4sar0 ( Alignment *A, Alignment *S, int comp, int leave, char *mode)
{
  return copy_aln (A,NULL);
}
Alignment * filter_aln4sar1 ( Alignment *inA, Alignment *S, int comp, int leave, char *mode)
{
  Alignment *F, *A;
  int a, b,c, i,r, n0, n1,g, score;
  int ***cache, **list1, **list2;
  int Delta;

  int T1;

  /*Keep only the positions where there are residues ONLY associated with 0 sequences*/

  list1=declare_int ( inA->nseq, 2);
  list2=declare_int ( inA->len_aln, 2);

  cache=(int***)declare_arrayN (3,sizeof (int),inA->len_aln,2, 26);
  F=copy_aln (inA, NULL);

  A=copy_aln (inA, NULL);
  A->nseq=strlen (S->seq_al[comp]);

  strget_param (mode, "_T1_", "5", "%d", &T1);
  for ( a=0; a< A->len_aln; a++)
    {
      n1=n0=g=0;
      for (b=0; b< A->nseq; b++)
	{
	  if ( b==leave) continue;
	  i=(S->seq_al[comp][b]=='I')?1:0;
	  r=tolower(A->seq_al[b][a]);
	  if ( r=='-')continue;
	  cache[a][i][r-'a']++;
	}
    }

  for (a=0; a< A->nseq; a++)
    for ( score=0,b=0; b<A->len_aln; b++)
      {
	r=tolower (A->seq_al[a][b]);
	if ( is_gap(r))continue;
	else if ( cache[b][0][r-'a'] && !cache[b][1][r-'a'])list1[a][0]++;
      }

  for (a=0; a< A->len_aln; a++)
    {
      for ( score=0,b=0; b< A->nseq; b++)
	{
	  r=tolower (A->seq_al[b][a]);
	  if ( r=='-')continue;
	  else r-='a';
	  if ( cache[a][0][r] && !cache[a][1][r])score ++;
	}
      list2[a][0]=a;
      list2[a][1]=score;
    }
  sort_int (list2, 2, 1, 0, F->len_aln-1);

  Delta=A->len_aln/(100/T1);
  for ( a=0; a< F->len_aln-Delta; a++)
    {
      b=list2[a][0];
      for ( c=0; c<F->nseq; c++)
	{
	  F->seq_al[c][b]='-';
	}
    }

  ungap_aln (F);
  free_aln (A);
  free_arrayN ( (void ***)cache, 3);
  free_arrayN ((void**)list1, 2);
  free_arrayN ((void**)list2, 2);

  return F;
}
Alignment * filter_aln4sar2 ( Alignment *inA, Alignment *S, int comp, int leave, char *mode)
{
  Alignment *F, *A;
  int a,b,r,ncat;
  int *cache;
  int max_ncat=10;

  /*Keep Low entropy columns that contain less than ncat categories of different amino acids*/
  /*REmove columns containing 10% or more gaps*/

  cache=(int*)vcalloc ( 500, sizeof (char));
  F=copy_aln (inA, NULL);
  A=copy_aln (inA, NULL);
  A->nseq=strlen (S->seq_al[comp]);
  for ( a=0; a< A->len_aln; a++)
    {
      for (ncat=0,b=0; b< A->nseq; b++)
	{
	  if ( b==leave) continue;

	  r=tolower(A->seq_al[b][a]);
	  if ( !cache[r])ncat++;
	  cache[r]++;
	}

      if ( ncat <max_ncat && ((cache['-']*100)/A->nseq)<10)
	{
	  ;
	}
      else
	{
	  for (b=0; b<F->nseq; b++)
	    {
	      r=tolower(F->seq_al[b][a]);
	      F->seq_al[b][a]='-';
	      cache[r]=0;
	    }
	}
      for (b=0; b<A->nseq; b++)
	  {
	    r=tolower(A->seq_al[b][a]);
	    cache[r]=0;
	  }
    }

  free_aln (A);
  ungap_aln (F);
  vfree (cache);
  return F;
}

Alignment * filter_aln4sar3 ( Alignment *inA, Alignment *S, int comp, int leave, char *mode)
{
  Alignment *F, *rA, *A;
  int a, b,c;
  int **list1;
  char *bufS, *bufA;
  int Delta;
  int T3;

  /*Keep the 10% positions most correlated with the 0/1 pattern*/

  A=copy_aln (inA, NULL);
  A->nseq=strlen (S->seq_al[comp]);
  F=copy_aln (inA, NULL);
  rA=rotate_aln (A, NULL);

  strget_param (mode, "_T3_", "10", "%d", &T3);


  list1=declare_int ( inA->len_aln, 2);
  bufA=(char*)vcalloc ( A->nseq+1, sizeof (char));
  bufS=(char*)vcalloc ( A->nseq+1, sizeof (char));

  sprintf ( bufS, "%s", S->seq_al[comp]);
  splice_out_seg(bufS,leave, 1);


  for (a=0; a< A->len_aln; a++)
    {
      char aa;
      list1[a][0]=a;
      sprintf (bufA, "%s", rA->seq_al[a]);
      splice_out_seg (bufA,leave,1);
      list1[a][1]=(int)sar_vs_seq3 ( bufS, bufA,0,NULL, &aa);
    }

  sort_int (list1, 2, 1, 0, F->len_aln-1);
  Delta=F->len_aln/(100/T3);
  for ( a=0; a< F->len_aln-Delta; a++)
    {
	  b=list1[a][0];

	  for ( c=0; c<F->nseq; c++)
	    {
	      F->seq_al[c][b]='-';
	    }

    }
  F->score_aln=list1[F->len_aln-1][1];
  ungap_aln (F);

  free_aln (rA);
  free_aln(A);
  free_arrayN ((void**)list1, 2);
  vfree (bufS);vfree (bufA);
  return F;
}
Alignment * filter_aln4sar4 ( Alignment *inA, Alignment *S, int comp, int leave, char *mode)
{
  Alignment *F, *A;
  int a, b,c, i,r, n0, n1,g,score;
  int ***cache, **list1, **list2;

  /*Keep only the positions where there are residues ONLY associated with 0 sequences*/

  list1=declare_int ( inA->nseq, 2);
  list2=declare_int ( inA->len_aln, 2);

  cache=(int***)declare_arrayN (3,sizeof (int),inA->len_aln,2, 26);
  F=copy_aln (inA, NULL);
  A=copy_aln (inA, NULL);
  A->nseq=strlen (S->seq_al[comp]);

  for ( a=0; a< A->len_aln; a++)
    {
      n1=n0=g=0;
      for (b=0; b< A->nseq; b++)
	{
	  if ( b==leave) continue;
	  i=(S->seq_al[comp][b]=='I')?1:0;
	  r=tolower(A->seq_al[b][a]);
	  if ( r=='-')continue;
	  cache[a][i][r-'a']++;
	  n1+=i;
	}
    }


  for (a=0; a< A->len_aln; a++)
    {
      for ( score=0,b=0; b< A->nseq; b++)
	{
	  r=tolower (F->seq_al[b][a]);
	  if ( r=='-')continue;
	  else r-='a';
	  if (cache[a][1][r]>=n1/2)score=1;
	}
      list2[a][0]=a;
      list2[a][1]=score;
    }


  for ( a=0; a< F->len_aln; a++)
    {
      if ( list2[a][1]==1);
      else
	{
	  b=list2[a][0];
	  for ( c=0; c<F->nseq; c++)
	    {
	      F->seq_al[c][b]='-';
	    }
	}
    }
  ungap_aln (F);
  free_aln (A);
  free_arrayN ( (void ***)cache, 3);
  free_arrayN ((void**)list1, 2);
  free_arrayN ((void**)list2, 2);

  return F;
}

Alignment * filter_aln4sar5 ( Alignment *inA, Alignment *S, int comp, int leave, char *mode)
{
  Alignment *F, *rA, *A;
  int a, b,c;
  int **list1;
  char *bufS, *bufA;
  int max;
  /*Look for the positions that show the best correlation between the sequence variation and the SAR*/

  A=copy_aln (inA, NULL);
  A->nseq=strlen (S->seq_al[comp]);

  rA=rotate_aln (inA, NULL);
  F=copy_aln (inA, NULL);

  list1=declare_int ( A->len_aln, 2);
  bufA=(char*)vcalloc ( A->nseq+1, sizeof (char));
  bufS=(char*)vcalloc ( A->nseq+1, sizeof (char));



  sprintf ( bufS, "%s", S->seq_al[comp]);
  splice_out_seg(bufS,leave, 1);


  for (a=0; a< A->len_aln; a++)
    {
      char aa;
      list1[a][0]=a;
      sprintf (bufA, "%s", rA->seq_al[a]);
      splice_out_seg (bufA,leave,1);
      list1[a][1]=(int)sar_vs_seq4 ( bufS, bufA,0,NULL, &aa);
    }

  sort_int (list1, 2, 1, 0, F->len_aln-1);
  max=F->score=list1[F->len_aln-1][1];
  max-=(max/10);


  for ( a=0; a< F->len_aln-10; a++)
    {

	  b=list1[a][0];

	  for ( c=0; c<F->nseq; c++)
	    {
	      F->seq_al[c][b]='-';
	    }

    }
  F->score_aln=10;
  ungap_aln (F);
  free_aln (inA);
  free_aln (rA);
  free_arrayN ((void**)list1, 2);
  vfree (bufS);vfree (bufA);
  return F;
}

int sar_profile2score ( char *seq, int **P)
{
  int a,r, l, score;

  l=strlen (seq);
  for ( score=0,a=0; a< l; a++)
    {
      r=seq[a];
      if ( is_gap(r))continue;
      score+=P[a][tolower(r)-'a'];
    }
  return score;
}
int make_sim_pred ( Alignment *A,Alignment *S, int comp, int seq)
{
  int a, b, i, r1, r2;
  static float **cscore;
  static float **tscore;

  if ( !cscore)
    {
      cscore=declare_float (2, 2);
      tscore=declare_float (2, 2);
    }

  for (a=0; a< 2; a++)for (b=0; b<2; b++)cscore[a][b]=tscore[a][b]=0;

  for ( a=0; a<A->len_aln; a++)
    {
      r1=A->seq_al[seq][a];
      if ( r1=='-') continue;
      else
	{
	  for ( b=0; b< A->nseq; b++)
	    {
	      if (b==seq) continue;
	      else
		{
		  r2=A->seq_al[b][a];
		  if (r2=='-')continue;
		  else
		    {

		      i=(S->seq_al[comp][b]=='I')?1:0;
		      cscore[i][0]+=(r1==r2)?1:0;
		      cscore[i][1]++;
		    }
		}
	    }

	  for (i=0; i<2; i++)
	    {
	      cscore[i][0]/=(cscore[i][1]==0)?1:cscore[i][1];
	      tscore[i][0]+=cscore[i][0];tscore[i][1]++;
	      cscore[i][0]=cscore[i][1]=0;
	    }
	}
    }

  fprintf ( stdout, "\nn\t 1: %.2f 0: %.2f", tscore[1][0],tscore[0][0]);
  return ( tscore[1][0]>=tscore[0][0])?1:0;
}


Alignment * sar_analyze (Alignment *inA, Alignment *inS, char *mode)
{
  int ***sim,***glob_results, ***comp_results;
  int *count;
  int a,b,c,m;
  float *tot2;
  Alignment *A=NULL,*S=NULL,*F, *SUBSET;
  char *subset, *target;
  int jack, T, filter;
  filter_func *ff;
  int n_methods=0;
  char *prediction, *reliability;
  int pred_start=0, pred_end, ref_start=0, ref_end;
  int display, CSV=1, NONCSV=0;
  char method[5];

  strget_param (mode, "_METHOD_", "1111", "%s_", method);
  ff=(filter_func*)vcalloc (6,sizeof (filter_func));
  if (method[0]=='1')ff[n_methods++]=filter_aln4sar0;
  if (method[1]=='1')ff[n_methods++]=filter_aln4sar1;
  if (method[2]=='1')ff[n_methods++]=filter_aln4sar2;
  if (method[3]=='1')ff[n_methods++]=filter_aln4sar3;
  /*
    ff[n_methods++]=filter_aln4sar4;
    ff[n_methods++]=filter_aln4sar5;
  */
  sim=(int***)vcalloc (n_methods, sizeof (int**));


  tot2=(float*)vcalloc ( 10, sizeof (float));
  subset=(char*)vcalloc ( 100, sizeof (char));
  target=(char*)vcalloc ( 100, sizeof (char));

  strget_param (mode, "_TARGET_", "no", "%s_", target);
  strget_param (mode, "_SUBSET_", "no", "%s_", subset);
  strget_param (mode, "_JACK_", "0", "%d", &jack);
  strget_param (mode, "_T_", "0", "%d", &T);
  strget_param (mode, "_FILTER_", "11", "%d", &filter);
  strget_param (mode, "_DISPLAY_", "0", "%d", &display);



  if ( !strm (target, "no"))
    {
      Alignment *T;
      T=main_read_aln(target, NULL);
      if ( T->len_aln !=inA->len_aln )
	{
	  printf_exit ( EXIT_FAILURE,stderr, "Error: %s is incompatible with the reference alignment [FATAL:%s]",target,PROGRAM);
	}

      inA=stack_aln (inA, T);

    }

  if ( !strm(subset, "no"))
    {
      SUBSET=main_read_aln (subset, NULL);
      sarset2subsarset ( inA, inS, &A, &S, SUBSET);
    }
  else
    {
      A=inA;
      S=inS;
    }


  prediction=(char*)vcalloc ( n_methods+1, sizeof (char));
  reliability=(char*)vcalloc ( n_methods+1, sizeof (char));

  glob_results=(int***)declare_arrayN(3, sizeof (int), n_methods*2, 2, 2);

  count=(int*)vcalloc (S->nseq, sizeof (int));
  for (a=0; a<S->nseq; a++)
    {
      int l;
      l=strlen (S->seq_al[a]);
      for ( b=0; b<l; b++)
	count[a]+=(S->seq_al[a][b]=='I')?1:0;
    }
  if ( display==CSV)
    {fprintf ( stdout, "\nCompound %s ; Ntargets %d", S->name[a],count[a]);
      pred_start=(strlen (S->seq_al[0])==A->nseq)?0:strlen (S->seq_al[0]);
      pred_end=A->nseq;
      for (a=pred_start; a< pred_end; a++)
	fprintf ( stdout, ";%s", A->name[a]);
      fprintf ( stdout, ";npred;");
    }


  for (a=0; a<S->nseq; a++)
    {
      int n_pred;
      comp_results=(int***)declare_arrayN(3, sizeof (int), n_methods*2, 2, 2);

      pred_start=(strlen (S->seq_al[a])==A->nseq)?0:strlen (S->seq_al[a]);
      pred_end=A->nseq;
      if ( display==CSV)fprintf ( stdout, "\n%s;%d", S->name[a],count[a]);

      for (n_pred=0,b=pred_start; b<pred_end;b++)
	{
	  int t, score=0,pred, real;

	  if ( display==NONCSV)fprintf ( stdout, "\n>%-15s %10s %c ", S->name[a], A->name[b], (pred_start==0)?S->seq_al[a][b]:'?');
	  if (jack || b==pred_start)
	    {
	      for (m=0; m<n_methods; m++)
		{
		  free_int (sim[m], -1);
		  F=(ff[m]) (A,S,a,(jack==0)?-1:b, mode);
		  sim[m]=aln2sim_mat(F, "idmat");
		  free_aln (F);
		}
	    }

	  for (m=0; m<n_methods; m++)
	    {
	      int Nbsim=0,Ybsim=0,bsim=0;
	      ref_start=0;
	      ref_end=strlen (S->seq_al[m]);

	      for (c=ref_start;c<ref_end; c++)
		{
		  if ( b==c) continue;
		  else if ( S->seq_al[a][c]=='O')
		    {
		      Nbsim=MAX(Nbsim,sim[m][b][c]);
		    }
		  else
		    {
		      Ybsim=MAX(Ybsim,sim[m][b][c]);
		    }
		}

	      bsim=(Ybsim>Nbsim)?Ybsim:-Nbsim;
	      pred=(bsim>0)?1:0;
	      real=(S->seq_al[a][b]=='O')?0:1;
	      comp_results[m][pred][real]++;
	      glob_results[m][pred][real]++;
	      score+=pred;
	      prediction[m]=pred+'0';
	      reliability[m]=(FABS((Ybsim-Nbsim))-1)/10+'0';
	    }

	  if ( score>0)n_pred++;
	  prediction[m]=reliability[m]='\0';
	  if (display==NONCSV)fprintf ( stdout, "Compound_Count:%d primary_predictions: %s Total: %d", count[a],prediction, score);
	  else if ( display==CSV)fprintf ( stdout, ";%d", score);
	  for (t=0; t<n_methods; t++)
	    {
	      if (score>t)
		{
		  comp_results[t+n_methods][1][real]++;
		  glob_results[t+n_methods][1][real]++;
		}
	      else
		{
		  comp_results[t+n_methods][0][real]++;
		  glob_results[t+n_methods][0][real]++;
		}
	    }
	}
      if ( display==NONCSV)
	{if ( pred_start==0)display_prediction (comp_results, S,a, n_methods*2);}
      else fprintf (stdout, ";%d;",n_pred);
    }
  if ( display==NONCSV)if (pred_start==0)display_prediction (glob_results, S,-1, n_methods*2);


  myexit (EXIT_SUCCESS);
}
float display_prediction (int ***count, Alignment *S, int c, int n)
{
  float tp,tn,fn,fp,sp,sn,sn2;
  int a, nm;

  nm=n/2;

  for (a=0; a<n; a++)
    {
      tp=count[a][1][1];
      tn=count[a][0][0];
      fp=count[a][1][0];
      fn=count[a][0][1];

      sn2=tp/(tp+fp);
      sn=tp/(tp+fn);
      sp=tn/(tn+fp);
      if ( a<nm)fprintf ( stdout, "\n>#Method %d Compound %15s sp=%.2f sn=%.2f sn2=%.2f",a, (c==-1)?"TOTAL":S->name[c],sp, sn, sn2 );
      else fprintf ( stdout, "\n>#Combined: T=%d Compound %15s sp=%.2f sn=%.2f sn2=%.2f",a-nm, (c==-1)?"TOTAL":S->name[c],sp, sn, sn2 );
    }
  fprintf ( stdout, "\n");
  return 0;
}

float display_prediction_2 (int **prediction, int n,Alignment *A, Alignment *S, int field)
{
  int a, t, T;
  float max_sn, max_sp;

  if ( field==17 || field ==18)
    {
      printf_exit ( EXIT_FAILURE, stderr, "\nERROR: Do not use filed %d in display_prediction", field);
    }

  sort_int_inv ( prediction, 10,field, 0, n-1);
  for (t=0,a=0; a<n; a++)
    {
      t+=prediction[a][3];
      prediction[a][17]=t;
    }

  for (t=0,a=n-1; a>=0; a--)
    {
      prediction[a][18]=t;
      t+=prediction[a][3];
    }

  max_sn=max_sp=T=0;
  for (a=0; a<n; a++)
    {
      float tp, fn, fp, sp, sn;

      tp=prediction[a][17];
      fn=prediction[a][18];
      fp=(a+1)-tp;

      sp=((tp+fp)==0)?0:tp/(tp+fp);
      sn=((tp+fn)==0)?0:tp/(tp+fn);

      if (sp>0.8)
	{
	  if (sn>max_sn)
	    {
	      max_sn=sn;
	      max_sp=sp;

	      T=prediction[a][field];
	    }
	}
    }
  if (max_sn==0)
      fprintf (stdout, "\n T =%d SN=%.2f SP= %.2f",T,max_sn,max_sp);
  else
      fprintf (stdout, "\n T =%d SN=%.2f SP= %.2f",T,max_sn,max_sp);

  return max_sn;
}


/************************************************************************************/
/*                NEW      ANALYZE     : SAR                                        */
/************************************************************************************/
float** cache2pred1 (Alignment *A,int**cache, int *ns, int **ls, Alignment *S, char *compound, char *mode);
float** cache2pred2 (Alignment *A,int**cache, int *ns, int **ls, Alignment *S, char *compound, char *mode);
float** cache2pred3 (Alignment *A,int**cache, int *ns, int **ls, Alignment *S, char *compound, char *mode);
float** cache2pred4 (Alignment *A,int**cache, int *ns, int **ls, Alignment *S, char *compound, char *mode);
float** cache2pred5 (Alignment *A,int**cache, int *ns, int **ls, Alignment *S, char *compound, char *mode);
float** cache2pred_new (Alignment *A,int**cache, int *ns, int **ls, Alignment *S, char *compound, char *mode);

int **sar2cache_adriana ( Alignment *A, int *ns,int **ls, Alignment *S, char *compound, char *mode);
int **sar2cache_proba_old ( Alignment *A, int *ns,int **ls, Alignment *S, char *compound, char *mode);
int **sar2cache_count1 ( Alignment *A, int *ns,int **ls, Alignment *S, char *compound, char *mode);
int **sar2cache_count2 ( Alignment *A, int *ns,int **ls, Alignment *S, char *compound, char *mode);
int **sar2cache_count3 ( Alignment *A, int *ns,int **ls, Alignment *S, char *compound, char *mode);

int **sar2cache_proba_new ( Alignment *A, int *ns,int **ls, Alignment *S, char *compound, char *mode);
int **sar2cache_proba2 ( Alignment *A, int *ns,int **ls, Alignment *S, char *compound, char *mode);
int **analyze_sar_compound1 ( char *name, char *seq, Alignment *A, char *mode);
int **analyze_sar_compound2 ( char *name, char *seq, Alignment *A, char *mode);

int aln2n_comp_col ( Alignment *A, Alignment *S, int ci);




int ***simple_sar_analyze_vot ( Alignment *inA, Alignment *SAR, char *mode);
int ***simple_sar_analyze_col ( Alignment *inA, Alignment *SAR, char *mode);



int sarset2subsarset ( Alignment *A, Alignment *S, Alignment **subA, Alignment **subS, Alignment *SUB);
int benchmark_sar (int v);
int aln2jack_group1 (Alignment *A, int seq, int **l1, int *nl1, int **l2, int *nl2);
int aln2jack_group2 (Alignment *A, int seq, int **l1, int *nl1, int **l2, int *nl2);
int aln2jack_group3 (Alignment *A, char *sar_seq, int **l1, int *nl1, int **l2, int *nl2);
float** jacknife5 (Alignment*A,int **cache, int *ns, int **ls, Alignment *S, char *compound, char *mode);
float** jacknife6 (Alignment*A,int **cache, int *ns, int **ls, Alignment *S, char *compound, char *mode);

int process_cache ( Alignment *A,Alignment *S, int ***Cache, char *mode);
Alignment *analyze_compounds (Alignment *A, Alignment *S, char *mode);

Alignment *analyze_compounds (Alignment *A, Alignment *S, char *mode)
{
  int a, b, c, tot, n;
  int **sim;
  int sar1, sar2;

  sim=aln2sim_mat (A, "idmat");
  for (a=0; a< S->nseq; a++)
    {
      for (n=0, tot=0, b=0; b< A->nseq-1; b++)
	{
	  sar1=(S->seq_al[a][b]=='I')?1:0;
	  for ( c=b+1; c<A->nseq; c++)
	    {
	      sar2=(S->seq_al[a][c]=='I')?1:0;

	      if (sar1 && sar2)
		{
		  tot+=sim[b][c];
		  n++;
		}
	    }
	}
      fprintf ( stdout, ">%-10s   CMPSIM: %.2f\n", S->name[a],(float)tot/(float)n);
    }
  free_int (sim, -1);
  return A;
}

int print_seq_pos ( int pos, Alignment *A, char *seq);
int abl1_evaluation (int p);
int print_seq_pos ( int pos, Alignment *A, char *seq)
{
  int a, b, s;

  s=name_is_in_hlist (seq, A->name, A->nseq);
  fprintf ( stdout, "S=%d", s);

  for (b=0,a=0; a<pos; a++)
    {
      if (!is_gap (A->seq_al[s][a]))b++;
    }
  fprintf ( stdout, "Pos %d SEQ %s: %d ", pos+1, seq, b+246);
  if ( strm ( seq, "ABL1")) fprintf ( stdout , "PT: %d", abl1_evaluation (b+246));
  return 0;
}

int process_cache ( Alignment *A,Alignment *S, int  ***Cache, char *mode)
{
  int a, b;
  int **pos, **pos2;
  int **C;
  int ab1, *ab1_pos;
  int weight_mode;

  strget_param ( mode, "_WEIGHT_", "1", "%d", &weight_mode);
  pos=declare_int(A->len_aln+1,2);
  pos2=declare_int (A->len_aln+1,S->nseq);
  for (a=0; a<S->nseq; a++)
    {
      C=Cache[a];
      for (b=0; b< A->len_aln; b++)
	{
	    pos[b][0]+=C[26][b];
	    if ( C[26][b]>0)
	      {
		pos[b][1]++;
		pos2[b][a]=1;
	      }
	}
    }

  C=Cache[0];
  ab1=name_is_in_hlist ("ABL1", A->name, A->nseq);
  ab1_pos=(int*)vcalloc (A->len_aln+1, sizeof (int));

  for ( b=0,a=0; a< A->len_aln; a++)
    {
      if ( A->seq_al[ab1][a]=='-')ab1_pos[a]=-1;
      else ab1_pos[a]=++b;
    }

  for ( a=0; a< A->len_aln; a++)
    {
      fprintf ( stdout, "\n%4d %5d %5d %5d [%c] [%2d] ALN", a+1, pos[a][0], pos[a][1], ab1_pos[a]+246,A->seq_al[ab1][a],abl1_evaluation (ab1_pos[a]+246));
      for ( b=0; b< S->nseq; b++)fprintf ( stdout, "%d", pos2[a][b]);
    }
  return 1;
}
int abl1_evaluation (int p)
{
  if ( p==248) return 10;
  if ( p==250) return 10;
  if ( p==253) return 10;
  if ( p==254) return 10;
  if ( p==255) return 9;
  if ( p==256) return 10;
  if ( p==257) return 5;
  if ( p==258) return 8;
  if ( p==269) return 8;
  if ( p==291) return 4;
  if ( p==294) return 8;
  if ( p==299) return 10;
  if ( p==306) return 0;
  if ( p==314) return 9;
  if ( p==315) return 10;
  if ( p==318) return 10;

  if ( p==319) return 10;
  if ( p==321) return 10;
  if ( p==323) return 0;
  if ( p==324) return 0;
  if ( p==339) return 0;
  if ( p==340) return 0;
  if ( p==355) return 5;
  if ( p==364) return 10;

  if ( p==366) return 0;
  if ( p==368) return 10;
  if ( p==370) return 10;
  if ( p==372) return 0;
  if ( p==378) return 8;
  if ( p==382) return 10;

  if ( p==384) return 10;
  if ( p==387) return 10;
  if ( p==395) return 8;

  if ( p==398) return 8;
  if ( p==399) return 8;
  if ( p==400) return 8;
  if ( p==403) return 0;
  if ( p==416) return 8;
  if ( p==419) return 5;
  if ( p>400) return 0;
  return -1;
}
float** cache2pred1 (Alignment*A,int **cache, int *ns, int **ls, Alignment *S, char *compound, char *mode)
{
  int s1, s2, seq1, seq2, r1, r2,col, pred, real, ci;
  double score, max, id, m;
  float **R, T;


  int used_col, used_res,is_used_col, n_res=0;
  int weight_mode;
  /*Predict on ns[1] what was trained on ns[0]*/

  strget_param ( mode, "_THR_", "0.09", "%f", &T);
  strget_param ( mode, "_WEIGHT_", "0", "%d", &weight_mode);

  R=declare_float (2, 2);
  ci=name_is_in_hlist ( compound, S->name, S->nseq);



  for (s1=0; s1<ns[1]; s1++)
    {
      int v;
      seq1=ls[1][s1];

      for (max=0,score=0, col=0; col<A->len_aln; col++)
	{
	  int max1;
	  r1=tolower (A->seq_al[seq1][col]);
	  for (max1=0,id=0, m=0,s2=0; s2<ns[0]; s2++)
	    {
	      seq2=ls[0][s2];
	      if ( S->seq_al[ci][seq2]=='O')continue;
	      if ( cache[seq2][col]==0 && !is_gap( A->seq_al[seq2][col]))continue;

	      r2=tolower ( A->seq_al[seq2][col]);
	      if ( is_gap(r2))continue;

	      v=(cache[seq2][col]>0 && weight_mode==1)?cache[seq2][col]:1;

	      max+=v;
	      if ( r2==r1)
		{
		  score+=v;
		}

	    }

	}
      pred=(( score/max) >T)?1:0;
      real=(S->seq_al[ci][seq1]=='I')?1:0;
      R[pred][real]++;

      fprintf ( stdout, "\n>%s %d%d SCORE %.2f C %s [SEQ]\n", A->name[seq1],real, pred, (float)score/(float)max, compound);
    }

  for (used_col=0,used_res=0,col=0; col<A->len_aln; col++)
    {
      for (is_used_col=0,s2=0; s2<ns[0]; s2++)
	{
	  seq2=ls[0][s2];
	  if ( cache[seq2][col]==0 && !is_gap(A->seq_al[seq2][col]))n_res++;
	  else if (is_gap(A->seq_al[seq2][col]));
	  else
	    {
	    is_used_col=1;
	    used_res++;
	    }
	}
      used_col+=is_used_col;
    }
  fprintf ( stdout, "\n>%s USED_POSITIONS: COL: %.2f RES: %.2f COMP\n", S->name[ci],  (float)used_col/(float)A->len_aln, (float)used_res/(float) n_res);

  return R;
}

float** cache2pred2 (Alignment*A,int **cache, int *ns, int **ls, Alignment *S, char *compound, char *mode)
{
  int s1, s2, seq1, seq2, r1, r2,col, pred, real, ci;
  double score, max;
  float **R, T;


  int used_col, used_res,is_used_col, n_res=0;
  /*Predict on ns[1] what was trained on ns[0]*/

  strget_param ( mode, "_THR_", "0.5", "%f", &T);


  R=declare_float (2, 2);
  ci=name_is_in_hlist ( compound, S->name,S->nseq);

  for (s1=0; s1<ns[1]; s1++)
    {
      int v;
      seq1=ls[1][s1];
      fprintf ( stdout, "\n");
      for (max=0,score=0, col=0; col<A->len_aln; col++)
	{
	  int used;

	  r1=tolower (A->seq_al[seq1][col]);
	  for (used=0,s2=0; s2<ns[0]; s2++)
	    {
	      seq2=ls[0][s2];

	      if ( S->seq_al[ci][seq2]=='O')continue;
	      if ( cache[seq2][col]==0 && !is_gap( A->seq_al[seq2][col]))continue;


	      r2=tolower ( A->seq_al[seq2][col]);
	      if ( is_gap(r2))continue;

	      v=cache[seq2][col];
	      if ( r2==r1){score+=v;}
	      used=1;
	      max+=v;
	    }
	  if (used) fprintf ( stdout, "%c", r1);
	}

      pred=(( score/max) >T)?1:0;
      real=(S->seq_al[ci][seq1]=='I')?1:0;
      R[pred][real]++;

      fprintf ( stdout, "PSEQ: %-10s SC: %4d MAX: %4d S: %.2f R: %4d", A->name[seq1],(int)score, (int)max, (float)score/max,real);

    }

  for (used_col=0,used_res=0,col=0; col<A->len_aln; col++)
    {
      for (is_used_col=0,s2=0; s2<ns[0]; s2++)
	{
	  seq2=ls[0][s2];
	  if ( cache[seq2][col]==0 && !is_gap(A->seq_al[seq2][col]))n_res++;
	  else if (is_gap(A->seq_al[seq2][col]));
	  else
	    {
	    is_used_col=1;
	    used_res++;
	    }
	}
      used_col+=is_used_col;
    }
  fprintf ( stdout, "\n>%s USED_POSITIONS: COL: %.2f RES: %.2f COMP\n", S->name[ci],  (float)used_col/(float)A->len_aln, (float)used_res/(float) n_res);

  return R;
}

float** cache2pred3 (Alignment*A,int **cache, int *ns, int **ls, Alignment *S, char *compound, char *mode)
{
  int s1, s2, seq1, seq2, r1, r2,col, pred, real, ci, a, n;
  double score, max;
  float **R, T;



  int tp, tn, fn, fp;
  int best_tp, best_fp;
  int delta, best_delta;
  int **list;

  /*Predict on ns[1] what was trained on ns[0]*/

  strget_param ( mode, "_THR_", "0.5", "%f", &T);


  R=declare_float (2, 2);
  ci=name_is_in_hlist ( compound, S->name, S->nseq);
  list=declare_int ( ns[1],3);

  for (s1=0; s1<ns[1]; s1++)
    {
      int v;
      seq1=ls[1][s1];

      for (max=0,score=0, col=0; col<A->len_aln; col++)
	{
	  int used;

	  r1=tolower (A->seq_al[seq1][col]);
	  for (used=0,s2=0; s2<ns[0]; s2++)
	    {
	      seq2=ls[0][s2];

	      if ( S->seq_al[ci][seq2]=='O')continue;
	      if ( cache[seq2][col]==0 && !is_gap( A->seq_al[seq2][col]))continue;


	      r2=tolower ( A->seq_al[seq2][col]);
	      if ( is_gap(r2))continue;

	      v=cache[seq2][col];
	      if ( r2==r1){score+=v;}
	      used=1;
	      max+=v;
	    }
	}



      pred=(( score/max) >T)?1:0;
      real=(S->seq_al[ci][seq1]=='I')?1:0;

      list[s1][0]=real;
      list[s1][1]=(int)((score/max)*(float)1000);
      list[s1][2]=seq1;



    }
  sort_int_inv (list, 3, 1, 0, ns[1]-1);

  for ( a=0; a<ns[1]; a++)
    {
      seq1=list[a][2];
      fprintf ( stdout, "PSEQ: %-10s SC: %5d R: %4d\n", A->name[seq1],list[a][0], list[a][1]);
    }

  for (n=0, a=0; a<ns[1]; a++)n+=list[a][0];
  for (best_delta=100000,best_tp=0,tp=0,fp=0,best_fp=0,a=0; a<ns[1]; a++)
    {
      tp+=list[a][0];
      fp+=1-list[a][0];
      delta=(n-(tp+fp));
      if (FABS(delta)<best_delta)
	{
	  best_delta=delta;
	  best_tp=tp;
	  best_fp=fp;
	}
    }
  /*R[pred][real]*/
  tp=best_tp;
  fp=best_fp;
  fn=n-tp;
  tn=ns[1]-(tp+fp+fn);
  R[1][1]=tp;
  R[1][0]=fp;
  R[0][1]=fn;
  R[0][0]=tn;
  free_int (list, -1);
  return R;
}
float** cache2pred4 (Alignment*A,int **cache, int *ns, int **ls, Alignment *S, char *compound, char *mode)
{
  int s1, s2, seq1, seq2, ci, a,b, c, n;
  double score;
  float **R;


  int tp, tn, fn, fp;
  int best_tp, best_fp;
  int delta, best_delta;
  int **list;
  int **sim;
  int *ul;
  int nused=0;

  /*Predict on ns[1] what was trained on ns[0]*/
  /*Identify interesting coloumns*/
  ul=(int*)vcalloc ( A->len_aln, sizeof (int));
  for (a=0; a< A->len_aln; a++)
    for ( b=0; b< A->nseq; b++)
      if ( cache[b][a])ul[nused++]=a;

  /*compute the similarity on the used columns*/

  R=declare_float (2, 2);
  sim=declare_int (A->nseq, A->nseq);
  for (a=0; a< A->nseq; a++)
    for ( b=0; b< A->nseq; b++)
      {
	for (c=0; c< nused; c++)
	  {
	    if ( A->seq_al[a][ul[c]]==A->seq_al[b][ul[c]])sim[a][b]++;
	  }
	sim[a][b]=(sim[a][b]*100)/nused;
      }
  vfree (ul);




  ci=name_is_in_hlist ( compound, S->name, S->nseq);
  list=declare_int ( ns[1],2);

  for (s1=0; s1<ns[1]; s1++)
    {

      seq1=ls[1][s1];

      for (score=0,s2=0; s2<ns[0]; s2++)
	{
	  seq2=ls[0][s2];

	  if ( seq1==seq2)continue;
	  if (S->seq_al[ci][seq2]=='I')score=MAX(score, sim[seq1][seq2]);
	}
      list[s1][0]=(S->seq_al[ci][seq1]=='I')?1:0;
      list[s1][1]=(int)score;

    }
  sort_int_inv (list, 2, 1, 0, ns[1]-1);

  for (n=0, a=0; a<ns[1]; a++)n+=list[a][0];
  for (best_delta=100000,best_tp=0,tp=0,fp=0,best_fp=0,a=0; a<ns[1]; a++)
    {
      tp+=list[a][0];
      fp+=1-list[a][0];
      delta=(n-(tp+fp));
      if (FABS(delta)<best_delta)
	{
	  best_delta=delta;
	  best_tp=tp;
	  best_fp=fp;
	}
    }
  /*R[pred][real]*/
  tp=best_tp;
  fp=best_fp;
  fn=n-tp;
  tn=ns[1]-(tp+fp+fn);
  R[1][1]=tp;
  R[1][0]=fp;
  R[0][1]=fn;
  R[0][0]=tn;
  free_int (list, -1);
  free_int (sim, -1);
  return R;
}

float** cache2pred5 (Alignment*A,int **cache, int *ns, int **ls, Alignment *S, char *compound, char *mode)
{
  int s1, s2, seq1, seq2, ci, a, n;
  double score;
  float **R;



  int tp, tn, fn, fp;
  int best_tp, best_fp;
  int delta, best_delta;
  int **list;
  static int **sim;

  /*Predict on ns[1] what was trained on ns[0]*/

  R=declare_float (2, 2);

  if ( sim==NULL)
    sim=aln2sim_mat (A, "idmat");



  ci=name_is_in_hlist ( compound, S->name, S->nseq);
  list=declare_int ( ns[1],2);

  for (s1=0; s1<ns[1]; s1++)
    {

      seq1=ls[1][s1];

      for (score=0,s2=0; s2<ns[0]; s2++)
	{
	  seq2=ls[0][s2];

	  if ( seq1==seq2)continue;
	  if (S->seq_al[ci][seq2]=='I')score=MAX(score, sim[seq1][seq2]);
	}
      list[s1][0]=(S->seq_al[ci][seq1]=='I')?1:0;
      list[s1][1]=(int)score;

    }
  sort_int_inv (list, 2, 1, 0, ns[1]-1);

  for (n=0, a=0; a<ns[1]; a++)n+=list[a][0];
  for (best_delta=100000,best_tp=0,tp=0,fp=0,best_fp=0,a=0; a<ns[1]; a++)
    {
      tp+=list[a][0];
      fp+=1-list[a][0];
      delta=(n-(tp+fp));
      if (FABS(delta)<best_delta)
	{
	  best_delta=delta;
	  best_tp=tp;
	  best_fp=fp;
	}
    }
  /*R[pred][real]*/
  tp=best_tp;
  fp=best_fp;
  fn=n-tp;
  tn=ns[1]-(tp+fp+fn);
  R[1][1]=tp;
  R[1][0]=fp;
  R[0][1]=fn;
  R[0][0]=tn;
  free_int (list, -1);
  return R;
}

float** jacknife5 (Alignment*A,int **cacheIN, int *ns, int **ls, Alignment *S, char *compound, char *mode)
{
  int seq1, ci, a,b, c, n;
  double score, max_score;
  float **R;


  int tp, tn, fn, fp;
  int best_tp, best_fp;
  int delta, best_delta;
  int **list;
  int **cache;

  ci=name_is_in_hlist ( compound, S->name, S->nseq);
  list=declare_int (A->nseq,2);
  R=declare_float (2, 2);


  for ( a=0; a<A->nseq; a++)
    {
      int real, res;

      ns[0]=A->nseq-1;
      ns[1]=1;
      for (c=0,b=0; b<A->nseq; b++)
	if (a!=b)ls[0][c++]=b;
      ls[1][0]=a;


      cache=sar2cache_count1 (A, ns, ls,S, compound, mode);
      for (b=0; b<=26; b++)
	for ( c=0; c< A->len_aln; c++)
	  cacheIN[b][c]+=cache[b][c];

      seq1=a;
      real=(S->seq_al[ci][seq1]=='I')?1:0;
      fprintf ( stdout, ">%-10s %d ", A->name[seq1], real);



      for (max_score=0,b=0; b<A->len_aln; b++)
	max_score+=cache[26][b];

      for (score=0,b=0; b<A->len_aln; b++)
	{
	  res=tolower (A->seq_al[seq1][b]);
	  if ( cache[26][b]==0) continue;
	  if ( !is_gap(res))
	    {
	      score+=cache[res-'a'][b];
	    }
	  /*fprintf ( stdout, "%c[%3d]", res,b);*/
	}
      fprintf ( stdout, " SCORE: %5d SPRED %d RATIO: %.2f \n", (int)score, a, (score*100)/max_score);
      list[a][0]=real;

      if ( strstr (mode, "SIMTEST"))list[a][1]=(score*100)/max_score;
      else list[a][1]=(score*100)/max_score;
      free_int (cache, -1);
    }


  sort_int_inv (list, 2, 1, 0, A->nseq-1);
  for (n=0, a=0; a<A->nseq; a++)
    {
      n+=list[a][0];
    }

  for (best_delta=100000,best_tp=0,tp=0,fp=0,best_fp=0,a=0; a<A->nseq; a++)
    {

      tp+=list[a][0];
      fp+=1-list[a][0];
      delta=(n-(tp+fp));
      if (FABS(delta)<best_delta)
	{
	  best_delta=delta;
	  best_tp=tp;
	  best_fp=fp;
	}
    }

  /*R[pred][real]*/
  tp=best_tp;
  fp=best_fp;
  fn=n-tp;
  tn=A->nseq-(tp+fp+fn);
  R[1][1]=tp;
  R[1][0]=fp;
  R[0][1]=fn;
  R[0][0]=tn;
  free_int (list, -1);

  return R;
}
float** jacknife6 (Alignment*A,int **cache, int *ns, int **ls, Alignment *S, char *compound, char *mode)
{
  int seq1, ci, a,b, c,d,e,f, n;
  double score;
  float **R;


  int tp, tn, fn, fp;
  int best_tp, best_fp;
  int delta, best_delta;
  int **list;

  ci=name_is_in_hlist ( compound, S->name, S->nseq);
  list=declare_int (A->len_aln,2);
  R=declare_float (2, 2);


  for ( a=0; a<A->nseq; a++)
    {
      int sar, res;
      int **new_cache;

      ns[0]=A->nseq-1;
      ns[1]=1;
      for (c=0,b=0; b<A->nseq; b++)
	if (a!=b)ls[0][c++]=b;
      ls[1][0]=a;

      cache=sar2cache_proba_new (A, ns, ls,S, compound, mode);


      new_cache=declare_int (27,A->len_aln);

      for (d=0; d< A->len_aln; d++)
	{
	  int **analyze;
	  if ( cache[26][d]==0)continue;
	  analyze=declare_int (26, 2);

	  for ( e=0; e< ns[0]; e++)
	    {
	      f=ls[0][e];
	      sar=(S->seq_al[ci][f]=='I')?1:0;
	      res=tolower (A->seq_al[f][d]);

	      if ( res=='-') continue;
	      analyze[res-'a'][sar]++;
	    }
	  for (e=0;e<26; e++)
	    {
	      if ( analyze[e][1]){new_cache[26][d]=1;new_cache[e][d]+=cache[e][d];}
	      /*
	      if ( analyze[e][0] && analyze[e][1]){new_cache[26][d]=1;new_cache[e][d]+=analyze[e][1];}
	      else if ( analyze[e][0]){new_cache[26][d]=1;new_cache[e][d]-=analyze[e][0]*10;}
	      else if ( analyze[e][1]){new_cache[26][d]=1;new_cache[e][d]+=analyze[e][1];}
	      else if ( !analyze[e][0] &&!analyze[e][1]);
	      */
	    }
	  free_int (analyze, -1);
	}

      seq1=a;
      sar=(S->seq_al[ci][seq1]=='I')?1:0;
      fprintf ( stdout, ">%-10s %d ", A->name[seq1], sar);

      for (score=0,b=0; b<A->len_aln; b++)
	{
	  res=tolower (A->seq_al[seq1][b]);
	  if ( cache[26][b]==0) continue;
	  if ( !is_gap(res))
	    {
	      score+=new_cache[res-'a'][b];
	    }
	}
      fprintf ( stdout, " SCORE: %5d SPRED\n", (int)score);
      list[seq1][0]=sar;
      list[seq1][1]=(int)score;

      free_int (new_cache, -1);
      free_int (cache, -1);
    }
  sort_int_inv (list, 2, 1, 0, A->nseq-1);
  for (n=0, a=0; a<A->nseq; a++)n+=list[a][0];
  for (best_delta=100000,best_tp=0,tp=0,fp=0,best_fp=0,a=0; a<A->nseq; a++)
    {
      tp+=list[a][0];
      fp+=1-list[a][0];
      delta=(n-(tp+fp));
      if (FABS(delta)<best_delta)
	{
	  best_delta=delta;
	  best_tp=tp;
	  best_fp=fp;
	}
    }

  fprintf ( stderr, "\n%d %d", best_tp, best_fp);
  /*R[pred][real]*/
  tp=best_tp;
  fp=best_fp;
  fn=n-tp;
  tn=A->nseq-(tp+fp+fn);
  R[1][1]=tp;
  R[1][0]=fp;
  R[0][1]=fn;
  R[0][0]=tn;
  free_int (list, -1);


  return R;
}
float** cache2pred_new (Alignment*A,int **cache, int *ns, int **ls, Alignment *S, char *compound, char *mode)
{
  int s1, seq1, ci, a,b, n;
  double score;
  float **R;


  int tp, tn, fn, fp;
  int best_tp, best_fp;
  int delta, best_delta;
  int **list;

  ci=name_is_in_hlist ( compound, S->name, S->nseq);
  list=declare_int ( ns[1],2);
  R=declare_float (2, 2);

  for (s1=0; s1<ns[1]; s1++)
    {
      int res, real;

      seq1=ls[1][s1];
      real=(S->seq_al[ci][seq1]=='I')?1:0;
      fprintf ( stdout, ">%-10s %d ", A->name[seq1], real);
      for (score=0,b=0; b<A->len_aln; b++)
	{
	  res=tolower (A->seq_al[seq1][b]);
	  if ( cache[26][b]==0) continue;
	  if ( !is_gap(res))
	    {
	      score+=cache[res-'a'][b];
	    }
	  fprintf ( stdout, "%c", res);
	}
      fprintf ( stdout, " SCORE: %5d SPRED\n", (int)score);
      list[s1][0]=real;
      list[s1][1]=(int)score;
    }

  sort_int_inv (list, 2, 1, 0, ns[1]-1);

  for (n=0, a=0; a<ns[1]; a++)n+=list[a][0];
  for (best_delta=100000,best_tp=0,tp=0,fp=0,best_fp=0,a=0; a<ns[1]; a++)
    {
      tp+=list[a][0];
      fp+=1-list[a][0];
      delta=(n-(tp+fp));
      if (FABS(delta)<best_delta)
	{
	  best_delta=delta;
	  best_tp=tp;
	  best_fp=fp;
	}
    }



  /*R[pred][real]*/
  tp=best_tp;
  fp=best_fp;
  fn=n-tp;
  tn=ns[1]-(tp+fp+fn);
  R[1][1]=tp;
  R[1][0]=fp;
  R[0][1]=fn;
  R[0][0]=tn;
  free_int (list, -1);


  return R;
}
float** cache2pred_forbiden_res (Alignment*A,int **cache, int *ns, int **ls, Alignment *S, char *compound, char *mode)
{
  int s1,seq1, ci, a,b, c, n;
  double score;
  float **R;


  int tp, tn, fn, fp;
  int best_tp, best_fp;
  int delta, best_delta;
  int **list;
  int **new_cache;
  int **mat;

  mat=read_matrice ( "blosum62mt");

  ci=name_is_in_hlist ( compound, S->name, S->nseq);
  list=declare_int ( ns[1],2);
  R=declare_float (2, 2);

  for (s1=0; s1<ns[1]; s1++)
    {
      int res, real;

      seq1=ls[1][s1];
      real=(S->seq_al[ci][seq1]=='I')?1:0;
      fprintf ( stdout, ">%-10s %d ", A->name[seq1], real);
      for (score=0,b=0; b<A->len_aln; b++)
	{
	  res=tolower (A->seq_al[seq1][b]);
	  if ( cache[26][b]==0) continue;
	  if ( !is_gap(res))
	    {
	      score+=cache[res-'a'][b];
	    }
	  fprintf ( stdout, "%c", res);
	}
      fprintf ( stdout, " SCORE: %5d SPRED\n", (int)score);
      list[s1][0]=real;
      list[s1][1]=(int)score;
    }
  new_cache=declare_int (27,A->len_aln);
  for (a=0; a< A->len_aln; a++)
    {
      int **analyze, real, res, d;
      int *res_type;
      int **sub;
      int *keep;
      keep=(int*)vcalloc ( 26, sizeof (int));
      res_type=(int*)vcalloc ( 26, sizeof (int));
      sub=declare_int (256, 2);

      if ( cache[26][a]==0)continue;
      analyze=declare_int (26, 2);
      for ( b=0; b< ns[0]; b++)
	{
	  seq1=ls[0][b];
	  real=(S->seq_al[ci][seq1]=='I')?1:0;
	  res=tolower (A->seq_al[seq1][a]);

	  if ( res=='-') continue;
	  analyze[res-'a'][real]++;
	}
      fprintf ( stdout, "RSPRED: ");
      for (c=0;c<26; c++)fprintf ( stdout, "%c", c+'a');
      fprintf ( stdout, "\nRSPRED: ");
      for (c=0;c<26; c++)
	{
	  if ( analyze[c][0] && analyze[c][1]){fprintf ( stdout, "1");res_type[c]='1';}
	  else if ( analyze[c][0]){new_cache[26][a]=1;new_cache[c][a]-=analyze[c][0];fprintf ( stdout, "0");res_type[c]='0';}
	  else if ( analyze[c][1]){new_cache[26][a]=1;new_cache[c][a]+=analyze[c][1];fprintf ( stdout, "1");res_type[c]='1';}
	  else if ( !analyze[c][0] &&!analyze[c][1]){fprintf ( stdout, "-");res_type[c]='-';}
	}


      for ( c=0; c<26; c++)
	{
	  for ( d=0; d<26; d++)
	    {

	      if ( res_type[c]==res_type[d])
		{
		  sub[res_type[c]][0]+=mat[c][d];
		  sub[res_type[c]][1]++;
		}
	      if ( res_type[c]!='-' && res_type[d]!='-')
		{
		  sub['m'][0]+=mat[c][d];
		  sub['m'][1]++;
		}
	    }
	}
      for ( c=0; c< 256; c++)
	{
	  if ( sub[c][1])fprintf ( stdout, " %c: %5.2f ", c, (float)sub[c][0]/(float)sub[c][1]);
	}
      fprintf ( stdout, " SC: %d\nRSPRED  ", cache[26][a]);

      for ( c=0; c<26; c++)
	if ( res_type[c]=='1')
	  {
	    for (d=0; d<26; d++)
	      if (mat[c][d]>0)keep[d]++;
	    keep[c]=9;
	  }

      for (c=0; c<26; c++)
	{
	  if ( keep[c]>10)fprintf ( stdout, "9");
	  else fprintf ( stdout, "%d", keep[c]);
	}
      for ( c=0; c<26; c++)
	{
	  if ( keep[c]>8)new_cache[c][a]=10;
	  else new_cache[c][a]=-10;
	}
      fprintf ( stdout, "\n");
      free_int (analyze, -1);
      free_int (sub, -1);
      vfree (res_type);
      vfree (keep);

    }
  for ( a=0; a<25; a++)
    for (b=a+1; b<26; b++)
      {
	int r1, r2;
	r1=a+'a';r2=b+'a';
	if ( strchr("bjoxz", r1))continue;
	if ( strchr("bjoxz",r2))continue;

	if ( mat[a][b]>0 && a!=b)fprintf ( stdout, "\nMATANALYZE %c %c %d", a+'a', b+'a', mat[a][b]);
      }

  for (s1=0; s1<ns[1]; s1++)
    {
      int res, real;

      seq1=ls[1][s1];
      real=(S->seq_al[ci][seq1]=='I')?1:0;
      fprintf ( stdout, ">%-10s %d ", A->name[seq1], real);
      for (score=0,b=0; b<A->len_aln; b++)
	{
	  res=tolower (A->seq_al[seq1][b]);
	  if ( cache[26][b]==0) continue;
	  if ( !is_gap(res))
	    {
	      score+=new_cache[res-'a'][b];
	    }
	  fprintf ( stdout, "%c", res);
	}
      fprintf ( stdout, " SCORE: %5d SPRED\n", (int)score);
      list[s1][0]=real;
      list[s1][1]=(int)score;
    }
  free_int (new_cache, -1);
  sort_int_inv (list, 2, 1, 0, ns[1]-1);


  for (n=0, a=0; a<ns[1]; a++)n+=list[a][0];
  for (best_delta=100000,best_tp=0,tp=0,fp=0,best_fp=0,a=0; a<ns[1]; a++)
    {
      tp+=list[a][0];
      fp+=1-list[a][0];
      delta=(n-(tp+fp));
      if (FABS(delta)<best_delta)
	{
	  best_delta=delta;
	  best_tp=tp;
	  best_fp=fp;
	}
    }



  /*R[pred][real]*/
  tp=best_tp;
  fp=best_fp;
  fn=n-tp;
  tn=ns[1]-(tp+fp+fn);
  R[1][1]=tp;
  R[1][0]=fp;
  R[0][1]=fn;
  R[0][0]=tn;
  free_int (list, -1);


  return R;
}

int **sar2cache_proba_old ( Alignment *A, int *ns,int **ls, Alignment *S, char *compound, char *mode)
{
  int col, s, seq,ms,mseq, res, mres, res1, n,maxn1, maxn2,maxn3, t, ci, a;
  float quant=0;
  int **list;

  int N1msa,N1sar, N, N11, N10, N01,N00, SCORE, COL_INDEX, RES;
  int nfield=0;
  int value;
  float T1, T2, T3, T4;
  int weight_mode;
  int **cache;
  static int **sim;
  int sim_weight, w, sw_thr;
  int train_mode;

  float zscore;

  RES=nfield++;COL_INDEX=nfield++;N1msa=nfield++;N1sar=nfield++;N=nfield++;N11=nfield++;N10=nfield++;N01=nfield++;N00=nfield++;SCORE=nfield++;
  ci=name_is_in_hlist ( compound, S->name, S->nseq);
  cache=declare_int (A->nseq, A->len_aln);

  strget_param ( mode, "_FILTER1_", "0"   , "%f", &T1);
  strget_param ( mode, "_FILTER2_", "1000000", "%f", &T2);
  strget_param ( mode, "_FILTER3_", "0"   , "%f", &T3);
  strget_param ( mode, "_FILTER4_", "1000000", "%f", &T4);
  strget_param ( mode, "_SIMWEIGHT_", "1", "%d", &sim_weight);
  strget_param ( mode, "_SWTHR_", "30", "%d", &sw_thr);
  strget_param (mode, "_TRAIN_","1", "%d", &train_mode);
  strget_param (mode, "_ZSCORE_","0", "%f", &zscore);





  if (sim_weight==1 && !sim) sim=aln2sim_mat(A, "idmat");
  for ( ms=0; ms<ns[0]; ms++)
    {
      mseq=ls[0][ms];
      if ( S->seq_al[ci][mseq]!='I')continue;

      list=declare_int (A->len_aln+1, nfield);
      for (t=0,n=0, col=0; col< A->len_aln; col++)
	{
	  int same_res;

	  mres=tolower(A->seq_al[mseq][col]);
	  list[col][RES]=mres;
	  list[col][COL_INDEX]=col;

	  if ( is_gap(mres))continue;
	  for ( s=0; s<ns[0]; s++)
	    {
	      seq=ls[0][s];
	      res=tolower(A->seq_al[seq][col]);
	      if (is_gap(res))continue;


	      if (sim_weight==1)
		{
		  w=sim[seq][mseq];w=(mres==res)?100-w:w;
		  if (w<sw_thr)w=0;
		}
	      else
		w=1;

	      if ( train_mode==4)
		{
		  if ( S->seq_al[ci][seq]=='I')same_res=1;
		  else same_res=(res==mres)?1:0;
		}
	      else
		same_res=(res==mres)?1:0;

	      list[col][N]+=w;

	      if (S->seq_al[ci][seq]=='I' && same_res)list[col][N11]+=w;
	      else if (S->seq_al[ci][seq]=='I' && same_res)list[col][N10]+=w;
	      else if (S->seq_al[ci][seq]=='O' && same_res)list[col][N01]+=w;
	      else if (S->seq_al[ci][seq]=='O' && same_res)list[col][N00]+=w;

	      if ( S->seq_al[ci][seq]=='I')list[col][N1sar]+=w;
	      if ( same_res)list[col][N1msa]+=w;

	    }

	  list[col][SCORE]=(int)evaluate_sar_score1 (list[col][N], list[col][N11], list[col][N1msa], list[col][N1sar]);

	}

      strget_param ( mode, "_MAXN1_", "5", "%d", &maxn1);
      strget_param ( mode, "_WEIGHT_", "1", "%d", &weight_mode);
      strget_param ( mode, "_QUANT_", "0.0", "%f", &quant);

      sort_int_inv (list,nfield,SCORE,0,A->len_aln-1);
      if ( quant !=0)
	{

	  n=quantile_rank ( list,SCORE, A->len_aln,quant);
	  sort_int (list,nfield,N1msa, 0, n-1);
	  maxn1=MIN(n,maxn1);
	}

      for (a=0; a<maxn1; a++)
	{
	  col=list[a][COL_INDEX];
	  res1=list[a][RES];
	  value=list[a][SCORE];
	  if ( value>T1 && value<T2){cache[mseq][col]= value;}
	}
      free_int (list, -1);
    }

  /*Filter Columns*/
  list=declare_int (A->len_aln+1, nfield);
  for ( col=0; col< A->len_aln; col++)
    {
      list[col][COL_INDEX]=col;
      for ( s=0; s<ns[0]; s++)
	{
	  seq=ls[0][s];
	  list[col][SCORE]+=cache[seq][col];
	}
    }

  /*Filter Columns with a score not between T2 and T3*/

  for (col=0; col< A->len_aln; col++)
    if (list[col][SCORE]<T3 || list[col][SCORE]>T4)
      {
	list[col][SCORE]=0;
	for (s=0; s< A->nseq; s++)
	  if (!is_gap(A->seq_al[s][col]))cache[s][col]=0;
      }

  /*Keep The N Best Columns*/
  if ( zscore!=0)
    {
      double sum=0, sum2=0, z;
      int n=0;
      for (a=0; a< A->len_aln; a++)
	{
	  if ( list[a][SCORE]>0)
	    {
	      sum+=list[a][SCORE];
	      sum2+=list[a][SCORE]*list[a][SCORE];
	      n++;
	    }
	}
      for (a=0; a<A->len_aln; a++)
	{
	  if ( list[a][SCORE]>0)
	    {
	      z=return_z_score (list[a][SCORE], sum, sum2,n);
	      if ((float)z<zscore)
		{
		  col=list[a][COL_INDEX];
		  for (s=0; s<A->nseq; s++)
		    cache [s][col]=0;
		}
	      else
		{
		  fprintf ( stdout, "\nZSCORE: KEEP COL %d SCORE: %f SCORE: %d\n", list[a][COL_INDEX], (float)z, list[a][SCORE]);
		}
	    }
	}
    }
  else
    {
      sort_int_inv (list,nfield,SCORE,0,A->len_aln-1);
      strget_param ( mode, "_MAXN2_", "100000", "%d", &maxn2);

      for (a=maxn2;a<A->len_aln; a++)
	{
	  col=list[a][COL_INDEX];
	  for (s=0; s<A->nseq; s++)
	    cache [s][col]=0;
	}
    }

  /*Get Rid of the N best Columns*/;
  strget_param ( mode, "_MAXN3_", "0", "%d", &maxn3);

  for (a=0; a<maxn3;a++)
    {
      col=list[a][COL_INDEX];
      for (s=0; s<A->nseq; s++)
	cache [s][col]=0;
    }

  return cache;
}
int aln2n_comp_col ( Alignment *A, Alignment *S, int ci)
{
  int  res, seq,sar, col, r;
  int **analyze;

  int tot=0;

  analyze=declare_int (27, 2);
  for ( col=0; col< A->len_aln; col++)
    {
      int n1, n0;


      for ( n1=0, n0=0,seq=0; seq<A->nseq; seq++)
	{
	  res=tolower(A->seq_al[seq][col]);
	  sar=(S->seq_al[ci][seq]=='I')?1:0;
	  n1+=(sar==1)?1:0;
	  n0+=(sar==0)?1:0;
	  if ( res=='-')continue;
	  res-='a';
	  analyze[res][sar]++;
	}

      for (r=0; r<26; r++)
	{
	  int a0,a1;
	  a0=analyze[r][0];
	  a1=analyze[r][1];


	  if ( a1==n1 && a0<n0)
	    {
	      tot++;
	    }
	}
      for ( r=0; r<26; r++)analyze[r][0]=analyze[r][1]=0;
    }

  free_int (analyze, -1);
  return tot;
}
int **sar2cache_count1 ( Alignment *A, int *ns,int **ls, Alignment *S, char *compound, char *mode)
{
  int  maxn2, res, seq,sar, ci, col,s, r;
  int **analyze, **list, **cache;
  static int **mat;

  int a0,a1, w;
  if (!mat) mat=read_matrice ("blosum62mt");


  list=declare_int ( A->len_aln, 2);
  cache=declare_int ( 27, A->len_aln);
  analyze=declare_int (27, 2);

  ci=name_is_in_hlist ( compound, S->name, S->nseq);

  for ( col=0; col< A->len_aln; col++)
    {
      int n1, n0;


      for ( n1=0, n0=0,s=0; s<ns[0]; s++)
	{
	  seq=ls[0][s];
	  res=tolower(A->seq_al[seq][col]);
	  sar=(S->seq_al[ci][seq]=='I')?1:0;
	  n1+=(sar==1)?1:0;
	  n0+=(sar==0)?1:0;
	  if ( res=='-')continue;
	  res-='a';

	  analyze[res][sar]++;
	}

      for (r=0; r<26; r++)
	{

	  a0=analyze[r][0];
	  a1=analyze[r][1];

	  if ( strstr (mode, "SIMTEST"))
	    {
	      w=a1;
	    }
	  else if (a1 )
	    {
	      w=n0-a0;
	    }
	  else w=0;

	  cache[r][col]+=w;
	  cache[26][col]=MAX(w, cache[26][col]);
	}

      for ( r=0; r<26; r++)analyze[r][0]=analyze[r][1]=0;
      list[col][0]=col;
      list[col][1]=cache[26][col];
    }

  free_int (analyze, -1);

  sort_int_inv (list, 2, 1, 0, A->len_aln-1);

  strget_param ( mode, "_MAXN2_", "100000", "%d", &maxn2);

  for ( col=maxn2; col<A->len_aln; col++)
    for ( r=0; r<=26; r++)cache[r][list[col][0]]=0;

  free_int (list, -1);
  return cache;
}


int **sar2cache_count2 ( Alignment *A, int *ns,int **ls, Alignment *S, char *compound, char *mode)
{
  int maxn2, res, seq,sar, ci, col,s, r;
  int **analyze, **list, **cache, **conseq;
  static int **mat;
  int w=0;
  if (!mat) mat=read_matrice ("blosum62mt");


  list=declare_int ( A->len_aln, 2);
  cache=declare_int ( 27, A->len_aln);
  conseq=declare_int ( A->len_aln,3);

  analyze=declare_int (27, 2);

  ci=name_is_in_hlist ( compound, S->name, S->nseq);
  for ( col=0; col< A->len_aln; col++)
    {
      int n1, n0;

      for ( n1=0, n0=0,s=0; s<ns[0]; s++)
	{
	  seq=ls[0][s];
	  res=tolower(A->seq_al[seq][col]);
	  sar=(S->seq_al[ci][seq]=='I')?1:0;
	  n1+=(sar==1)?1:0;
	  n0+=(sar==0)?1:0;
	  if ( res=='-')continue;
	  res-='a';
	  analyze[res][sar]++;
	}
      for (r=0; r<26; r++)
	{
	  int a0,a1;
	  a0=analyze[r][0];
	  a1=analyze[r][1];
	  if ( a1==n1 && a0<n0)
	    {

	      w=n0-a0;
	      conseq[col][0]=r;
	      conseq[col][1]=w;
	    }
	}
      for ( r=0; r<26; r++)analyze[r][0]=analyze[r][1]=0;
    }
  free_int (analyze, -1);

  for (s=0; s<ns[0]; s++)
    {
      int w1, w2;
      seq=ls[0][s];
      for (w1=0,w2=0,col=0; col<A->len_aln; col++)
	{

	  res=tolower(A->seq_al[seq][col]);
	  if ( is_gap(res))continue;
	  else res-='a';

	  if ( conseq[col][1] && res!=conseq[col][0])w1++;
	  if ( conseq[col][1])w2++;
	}
      for (col=0; col<A->len_aln; col++)
	{
	  res=tolower(A->seq_al[seq][col]);
	  if ( is_gap(res))continue;
	  else res-='a';

	  if ( conseq[col][1] && res!=conseq[col][0])conseq[col][2]+=(w2-w1);
	}
    }

  for (col=0; col<A->len_aln; col++)
    {
      r=conseq[col][0];
      w=conseq[col][2];


      cache[r][col]=cache[26][col]=list[col][1]=w;
      list[col][0]=col;
    }
  sort_int_inv (list, 2, 1, 0, A->len_aln-1);
  strget_param ( mode, "_MAXN2_", "100000", "%d", &maxn2);

  for ( col=maxn2; col<A->len_aln; col++)
    for ( r=0; r<=26; r++)cache[r][list[col][0]]=0;


  free_int (list, -1);
  return cache;
}

int **sar2cache_count3 ( Alignment *A, int *ns,int **ls, Alignment *S, char *compound, char *mode)
{
  int  maxn2, res, seq,sar, ci, col,s, r, a1, a0, n1, n0;
  int **analyze, **list, **cache;
  static int **mat;

  if (!mat) mat=read_matrice ("blosum62mt");


  list=declare_int ( A->len_aln, 2);
  cache=declare_int ( 27, A->len_aln);
  analyze=declare_int (27, 2);

  ci=name_is_in_hlist ( compound, S->name, S->nseq);

  for ( col=0; col< A->len_aln; col++)
    {
      double e, g;
      for ( n1=0, n0=0,s=0; s<ns[0]; s++)
	{
	  seq=ls[0][s];
	  res=tolower(A->seq_al[seq][col]);
	  sar=(S->seq_al[ci][seq]=='I')?1:0;
	  n1+=(sar==1)?1:0;
	  n0+=(sar==0)?1:0;
	  if ( res=='-')continue;
	  res-='a';

	  analyze[res][sar]++;
	}

      /*Gap*/
      for (g=0,r=0; r<A->nseq; r++)
	g+=is_gap(A->seq_al[r][col]);
      g=(100*g)/A->nseq;

      /*enthropy
      for (e=0, r=0; r<26; r++)
	{
	  a0=analyze[r][0];
	  a1=analyze[r][1];
	  t=a0+a1;

	  if (t>0)
	    e+= t/(double)A->nseq*log(t/(double)A->nseq);
	}
      e*=-1;
      */
      e=0;
      if (g>10) continue;
      if (e>10) continue;

      if ( strstr ( mode, "SIMTEST"))
	{
	  for (r=0; r<26; r++)
	    {

	      a0=analyze[r][0];
	      a1=analyze[r][1];

	      if (a1)
		{
		  cache[r][col]=a1;
		  cache[26][col]=MAX(cache[26][col],a1);
		}
	    }
	}
      else
	{



	  for (r=0; r<26; r++)
	    {

	      a0=analyze[r][0];
	      a1=analyze[r][1];

	      if (!a1 && a0)
		{
		  cache[r][col]=a0;
		  cache[26][col]=MAX(cache[26][col],a0);
		}
	    }
	}

      for ( r=0; r<26; r++)analyze[r][0]=analyze[r][1]=0;
      list[col][0]=col;
      list[col][1]=cache[26][col];
    }

  free_int (analyze, -1);

  sort_int_inv (list, 2, 1, 0, A->len_aln-1);

  strget_param ( mode, "_MAXN2_", "100000", "%d", &maxn2);

  for ( col=maxn2; col<A->len_aln; col++)
    for ( r=0; r<=26; r++)cache[r][list[col][0]]=0;

  free_int (list, -1);
  return cache;
}


int **sar2cache_proba_new ( Alignment *A, int *ns,int **ls, Alignment *S, char *compound, char *mode)
{
  int col, s, seq,ms,mseq, res, mres, res1, n,maxn1, maxn2,maxn3, t, ci, a,w;

  int **list;

  int N1msa,N1sar, N, N11, N10, N01,N00, SCORE, COL_INDEX, RES;
  int nfield=0;
  int value;


  int **cache;
  static int **sim;
  int sw_thr;
  float zscore;

  RES=nfield++;COL_INDEX=nfield++;N1msa=nfield++;N1sar=nfield++;N=nfield++;N11=nfield++;N10=nfield++;N01=nfield++;N00=nfield++;SCORE=nfield++;
  ci=name_is_in_hlist ( compound, S->name, S->nseq);
  cache=declare_int (27, A->len_aln);

  strget_param ( mode, "_SWTHR_", "30", "%d", &sw_thr);
  strget_param (mode, "_ZSCORE_","0", "%f", &zscore);


  if (!sim)sim=aln2sim_mat(A, "idmat");
  for ( ms=0; ms<ns[0]; ms++)
    {
      mseq=ls[0][ms];
      if ( S->seq_al[ci][mseq]!='I')continue;

      list=declare_int (A->len_aln+1, nfield);
      for (t=0,n=0, col=0; col< A->len_aln; col++)
	{
	  int same_res;

	  mres=tolower(A->seq_al[mseq][col]);
	  if ( is_gap(mres))continue;

	  list[col][RES]=mres;
	  list[col][COL_INDEX]=col;

	  for ( s=0; s<ns[0]; s++)
	    {
	      seq=ls[0][s];
	      res=tolower(A->seq_al[seq][col]);
	      if (is_gap(res))continue;
	      w=sim[seq][mseq];w=(mres==res)?100-w:w;
	      if (w<sw_thr)w=0;
	      same_res=(res==mres)?1:0;

	      list[col][N]+=w;

	      if (S->seq_al[ci][seq]=='I' && same_res)list[col][N11]+=w;
	      else if (S->seq_al[ci][seq]=='I' && same_res)list[col][N10]+=w;
	      else if (S->seq_al[ci][seq]=='O' && same_res)list[col][N01]+=w;
	      else if (S->seq_al[ci][seq]=='O' && same_res)list[col][N00]+=w;

	      if ( S->seq_al[ci][seq]=='I')list[col][N1sar]+=w;
	      if ( same_res)list[col][N1msa]+=w;

	    }

	  list[col][SCORE]=(int)evaluate_sar_score1 (list[col][N], list[col][N11], list[col][N1msa], list[col][N1sar]);

	}
      strget_param ( mode, "_MAXN1_", "5", "%d", &maxn1);
      sort_int_inv (list,nfield,SCORE,0,A->len_aln-1);
      for (a=0; a<maxn1; a++)
	{
	  col=list[a][COL_INDEX];
	  res1=list[a][RES];
	  value=list[a][SCORE];

	  if ( res1!=0)
	    {
	      cache[res1-'a'][col]+= value;
	      cache[26][col]+=value;
	    }
	}
      free_int (list, -1);
    }

  /*Filter Columns*/
  list=declare_int (A->len_aln+1, nfield);
  for ( col=0; col< A->len_aln; col++)
    {
      list[col][COL_INDEX]=col;
      list[col][SCORE]=cache[26][col];
    }
  /*Keep The N Best Columns*/
  if ( zscore!=0)
    {
      double sum=0, sum2=0, z;
      int n=0;
      for (a=0; a< A->len_aln; a++)
	{
	  if ( list[a][SCORE]>0)
	    {
	      sum+=list[a][SCORE];
	      sum2+=list[a][SCORE]*list[a][SCORE];
	      n++;
	    }
	}
      for (a=0; a<A->len_aln; a++)
	{
	  if ( list[a][SCORE]>0)
	    {
	      z=return_z_score (list[a][SCORE], sum, sum2,n);
	      if ((float)z<zscore)
		{
		  col=list[a][COL_INDEX];
		  for (s=0; s<27; s++)
		    cache [s][col]=0;
		}
	      else
		{
		  fprintf ( stdout, "\nZSCORE: KEEP COL %d SCORE: %f SCORE: %d\n", list[a][COL_INDEX], (float)z, list[a][SCORE]);
		}
	    }
	}
    }
  else
    {
      sort_int_inv (list,nfield,SCORE,0,A->len_aln-1);
      strget_param ( mode, "_MAXN2_", "100000", "%d", &maxn2);

      for (a=maxn2;a<A->len_aln; a++)
	{
	  col=list[a][COL_INDEX];
	  for (s=0; s<27; s++)
	    cache [s][col]=0;
	}
    }

  /*Get Rid of the N best Columns*/;
  strget_param ( mode, "_MAXN3_", "0", "%d", &maxn3);

  for (a=0; a<maxn3;a++)
    {
      col=list[a][COL_INDEX];
      for (s=0; s<27; s++)
	cache [s][col]=0;
    }
  return cache;
}
int **sar2cache_adriana ( Alignment *A, int *ns,int **ls, Alignment *S, char *compound, char *mode)
{

  int col,maxn1, s, seq,ms,mseq, res, mres,res1, n, t, ci, a;
  float quant=0;
  int **list;


  int **cache;
  ci=name_is_in_hlist ( compound, S->name, S->nseq);
  cache=declare_int (A->nseq, A->len_aln);


  for ( ms=0; ms<ns[0]; ms++)
    {
      mseq=ls[0][ms];
      if ( S->seq_al[ci][mseq]!='I')continue;

      list=declare_int (A->len_aln+1, 5);
      for (t=0,n=0, col=0; col< A->len_aln; col++)
	{
	  mres=tolower(A->seq_al[mseq][col]);
	  list[col][0]=mres;
	  list[col][1]=col;

	  if ( is_gap(mres))continue;
	  for ( s=0; s<ns[0]; s++)
	    {
	      seq=ls[0][s];
	      res=tolower(A->seq_al[seq][col]);
	      if (is_gap(res))continue;

	      if (S->seq_al[ci][seq]=='I' && res==mres)list[col][3]++;
	      if (res==mres)list[col][2]++;
	    }
	}

      sort_int_inv (list,5,3,0,A->len_aln-1);

      strget_param ( mode, "_MAXN1_", "5", "%d", &maxn1);
      strget_param ( mode, "_QUANT_", "0.95", "%f", &quant);

      n=quantile_rank ( list, 3, A->len_aln,quant);
      sort_int (list, 5, 2, 0, n-1);

      for (a=0; a<maxn1; a++)
	{

	  col=list[a][1];
	  res1=list[a][0];
	  cache[mseq][col]=list[a][3];
	}
      free_int (list, -1);

    }
  return cache;
}
int **sar2cache_proba2 ( Alignment *A, int *ns,int **ls, Alignment *S, char *compound, char *mode)
{
  int col, s, seq,ms,mseq, res, mres,n,maxn1, t, ci, a,b;
  int COL, SCORE;

  float quant=0;
  int **list;

  float T1, T2, T3, T4;

  int **cache;
  cache=declare_int ( A->nseq, A->len_aln);
  ci=name_is_in_hlist ( compound, S->name, S->nseq);

  strget_param ( mode, "_FILTER1_", "0"   , "%f", &T1);
  strget_param ( mode, "_FILTER2_", "1000000", "%f", &T2);
  strget_param ( mode, "_FILTER3_", "0"   , "%f", &T3);
  strget_param ( mode, "_FILTER4_", "1000000", "%f", &T4);

  list=declare_int (A->len_aln+1,A->nseq+2);
  SCORE=A->nseq;
  COL=A->nseq+1;

  for ( ms=0; ms<ns[0]; ms++)
    {
      mseq=ls[0][ms];
      if ( S->seq_al[ci][mseq]!='I')continue;

      for (t=0,n=0, col=0; col< A->len_aln; col++)
	{
	  int N11=0,N10=0,N01=0,N00=0,N1sar=0,N1msa=0,N=0;

	  mres=tolower(A->seq_al[mseq][col]);
	  if ( is_gap(mres))continue;
	  for ( s=0; s<ns[0]; s++)
	    {
	      seq=ls[0][s];
	      res=tolower(A->seq_al[seq][col]);
	      if (is_gap(res))continue;

	      N++;
	      if (S->seq_al[ci][seq]=='I' && res==mres)N11++;
	      else if (S->seq_al[ci][seq]=='I' && res!=mres)N10++;
	      else if (S->seq_al[ci][seq]=='O' && res==mres)N01++;
	      else if (S->seq_al[ci][seq]=='O' && res!=mres)N00++;

	      if ( S->seq_al[ci][seq]=='I')N1sar++;
	      if ( res==mres)N1msa++;
	    }
	  list[col][mseq]=(int)evaluate_sar_score1 (N,N11,N1msa,N1sar);
	  list[col][SCORE]+=list[col][mseq];
	  list[col][COL]=col;
	}
    }

  strget_param ( mode, "_MAXN1_", "5", "%d", &maxn1);
  strget_param ( mode, "_QUANT_", "0.95", "%f", &quant);
  sort_int_inv (list,A->nseq+2,SCORE, 0, A->len_aln-1);
  n=quantile_rank ( list,A->nseq, A->len_aln,quant);
  n=5;


  for (a=0; a<n; a++)
    {
      int value;

      col=list[a][COL];
      for ( b=0; b<A->nseq; b++)
	{
	  value=list[col][b];
	  if ( value>T1 && value<T2){cache[b][col]= value;}
	}
    }

  free_int (list, -1);
  return cache;
}




/************************************************************************************/
/*                ALIGNMENT ANALYZE     : SAR                                            */
/************************************************************************************/
int aln2jack_group3 (Alignment *A,char *comp, int **l1, int *nl1, int **l2, int *nl2)
{
  int **seq_list, **sar_list, nsar=0, nseq=0;
  int a, b, mid;

  vsrand (0);
  sar_list=declare_int (A->nseq, 2);
  seq_list=declare_int (A->nseq, 2);
  for (a=0; a< A->nseq; a++)
    {
      if (comp[a]=='I')
	{
	  sar_list[nsar][0]=a;
	  sar_list[nsar][1]=rand()%100000;
	  nsar++;
	}
      else
	{
	  seq_list[nseq][0]=a;
	  seq_list[nseq][1]=rand()%100000;
	  nseq++;
	}
    }


  l1[0]=(int*)vcalloc (A->nseq, sizeof (int));
  l2[0]=(int*)vcalloc (A->nseq, sizeof (int));
  nl1[0]=nl2[0]=0;

  sort_int (seq_list, 2, 1, 0,nseq-1);
  sort_int (sar_list, 2, 1, 0,nsar-1);
  mid=nsar/2;
  for (a=0; a<mid; a++)
    {
      l1[0][nl1[0]++]=sar_list[a][0];
    }
  for (a=0,b=mid; b<nsar; b++, a++)
    {
      l2[0][nl2[0]++]=sar_list[b][0];
    }

  mid=nseq/2;
  for (a=0; a<mid; a++)
    {
      l1[0][nl1[0]++]=seq_list[a][0];
    }
  for (a=0,b=mid; b<nseq; b++, a++)
    {
      l2[0][nl2[0]++]=seq_list[b][0];
    }


  free_int (seq_list, -1);
  free_int (sar_list, -1);
  return 1;
}

int aln2jack_group2 (Alignment *A, int seq, int **l1, int *nl1, int **l2, int *nl2)
{
  int **list;
  int a, b, mid;


  list=declare_int (A->nseq, 2);
  l1[0]=(int*)vcalloc (A->nseq, sizeof (int));
  l2[0]=(int*)vcalloc (A->nseq, sizeof (int));
  nl1[0]=nl2[0];

  vsrand (0);
  for ( a=0; a< A->nseq; a++)
    {
      list[a][0]=a;
      list[a][1]=rand()%100000;
    }
  sort_int (list, 2, 1, 0,A->nseq-1);
  mid=A->nseq/2;
  for (a=0; a<mid; a++)
    {
      l1[0][nl1[0]++]=list[a][0];
    }
  for (a=0,b=mid; b<A->nseq; b++, a++)
    {
      l2[0][nl2[0]++]=list[b][0];
    }

  free_int (list, -1);
  return 1;
}
int aln2jack_group1 (Alignment *A, int seq, int **l1, int *nl1, int **l2, int *nl2)
{
  int **sim;
  int **list;
  int a, mid;

  list=declare_int ( A->nseq, 3);
  l1[0]=(int*)vcalloc (A->nseq, sizeof (int));
  l2[0]=(int*)vcalloc (A->nseq, sizeof (int));
  nl1[0]=nl2[0];

  sim=aln2sim_mat (A, "idmat");
  for ( a=0; a< A->nseq; a++)
    {
      list[a][0]=seq;
      list[a][1]=a;
      list[a][2]=(a==seq)?100:sim[seq][a];
    }
  sort_int_inv (list, 3, 2, 0, A->nseq-1);
  fprintf ( stderr, "\nJacknife fromsequence %s [%d]\n", A->name[seq], seq);
  mid=A->nseq/2;
  for (a=0; a< mid; a++)
    l1[0][nl1[0]++]=list[a][1];
  for (a=mid; a<A->nseq; a++)
    l2[0][nl2[0]++]=list[a][1];
  return 1;
}


int sarset2subsarset ( Alignment *A, Alignment *S, Alignment **subA, Alignment **subS, Alignment *SUB)
{
  Alignment *rotS, *intS;
  int a,b, *list, nl;

  list=(int*)vcalloc ( SUB->nseq, sizeof (int));
  for (nl=0,a=0; a<SUB->nseq; a++)
    {
      b=name_is_in_hlist(SUB->name[a], A->name, A->nseq);
      if ( b!=-1)list[nl++]=b;
    }

  subA[0]=extract_sub_aln (A, nl, list);
  rotS=rotate_aln (S, NULL);
  intS=extract_sub_aln (rotS, nl, list);

  subS[0]=rotate_aln (intS, NULL);

  for ( a=0; a<S->nseq; a++) sprintf ( (subS[0])->name[a], "%s", S->name[a]);


  return 0;
}

int ***simple_sar_analyze_vot ( Alignment *A, Alignment *SAR, char *mode)
{
  int a, b, c, d;
  int res1, res2, sar1, sar2;
  float s;
  int **sim;
  static float ***result;
  static int ***iresult;
  if (!result)
    {
    result=(float***)declare_arrayN (3,sizeof (float),SAR->nseq, A->len_aln,3);
    iresult=(int***)declare_arrayN (3,sizeof (int),SAR->nseq, A->len_aln,3);
    }

  sim=aln2sim_mat (A, "idmat");


  for (a=0; a<SAR->nseq; a++)
    for (b=0; b<A->len_aln; b++)
      result[a][b][0]=1;

  for ( a=0; a< SAR->nseq; a++)
    for ( b=0; b<A->nseq-1; b++)
      for ( c=b+1; c< A->nseq; c++)
	for ( d=0; d<A->len_aln; d++)
	  {
	    res1=A->seq_al[b][d];
	    res2=A->seq_al[c][d];

	    sar1=(SAR->seq_al[a][b]=='I')?1:0;
	    sar2=(SAR->seq_al[a][c]=='I')?1:0;

	    s=sim[b][c];




	    if ( sar1!=sar2 && res1!=res2)
	      result[a][d][0]*=(1/(100-s));

	    else if ( sar1==sar2 && sar1==1 && res1==res2)
	      result[a][d][0]*=1/s;




	    /*
	    else if ( sar1==sar2 && res1==res2)result[a][d][0]+=(100-s)*(100-s);
	    else if ( sar1==sar2 && res1!=res2)result[a][d][0]-=s*s;
	    else if ( sar1!=sar2 && res1==res2)result[a][d][0]-=(100-s)*(100-s);
	    */

	    result[a][d][1]='a';
	  }
  for ( a=0; a<SAR->nseq; a++)
    for ( b=0; b<A->len_aln; b++)
      {
	fprintf ( stderr, "\n%f", result[a][b][0]);
	iresult[a][b][0]=100*log(1-result[a][b][0]);
      }
  return iresult;
}
int display_simple_sar_analyze_pair_col (Alignment *A, Alignment *SAR, char *mode)
{
  int **r;
  int a, b, n, do_tm;


  Alignment *rA;
  int *nI;


  strget_param (mode, "_TM_", "0", "%d", &do_tm);
  r=simple_sar_analyze_pair_col (A, SAR,mode);
  rA=rotate_aln (A, NULL);
  n=0;

  nI=(int*)vcalloc ( SAR->nseq, sizeof (int));
  for (a=0; a< SAR->nseq; a++)
    for (b=0; b<SAR->len_aln; b++) nI[a]+=(SAR->seq_al[a][b]=='I')?1:0;



  while ( r[n][0]!=-1)
    {
      if (r[n][3]>0)
	{fprintf ( stdout,  "COMP S: %3d %3d %s %20s %2d #\n", r[n][3],0,SAR->seq_al[r[n][0]], SAR->name[r[n][0]], nI[r[n][0]]);
	  fprintf ( stdout, "SEQ1 S: %3d %3d %s %20s %2d #\n", r[n][3],r[n][1],rA->seq_al[r[n][1]], SAR->name[r[n][0]],nI[r[n][0]]);
	  fprintf ( stdout, "SEQ2 S: %3d %3d %s %20s %2d #\n\n", r[n][3],r[n][2],rA->seq_al[r[n][2]], SAR->name[r[n][0]],nI[r[n][0]]);
	}
      n++;
    }
  return 0;
}
int display_simple_sar_analyze_col (Alignment *A, Alignment *SAR, char *mode)
{
  int ***result, **r2, **r3, **r4, **aa;
  int a, b, c, n;
  char *cons;
  int threshold=20;
  int do_tm;
  strget_param (mode, "_TM_", "0", "%d", &do_tm);
  result=simple_sar_analyze_col (A, SAR,mode);
  r2=declare_int (A->len_aln*SAR->nseq, 5);
  r3=declare_int (A->len_aln+1, 5);
  r4=declare_int (A->len_aln+1, SAR->nseq+1);
  aa=declare_int (2, 256);
  cons=(char*)vcalloc (A->len_aln+1, sizeof (char));
  for (a=0; a<A->len_aln; a++){r3[a][0]=a;cons[a]='A';}



  for (n=0,a=0; a< SAR->nseq; a++)
    {
      double sum, sum2;
      for (sum=0, sum2=0,b=0; b<A->len_aln; b++)
	{
	  sum+=result[a][b][0];
	  sum2+=result[a][b][0]*result[a][b][0];
	}

      for (b=0; b<A->len_aln; b++, n++)
	{
	  r2[n][0]=a;//compound
	  r2[n][1]=b;//pos
	  r2[n][2]=result[a][b][1]; //AA
	  r2[n][3]=result[a][b][0]; //Score
	  r2[n][4]=result[a][b][2];  //(int)10*return_z_score ((double)result[a][b][0], sum, sum2, A->len_aln); //ZScore
      }
    }
  sort_int (r2,5, 3, 0, n-1);//sort on Score (3rd field)
  for ( a=0; a< n; a++)
    {
      int comp, pos, bad;

      comp=r2[a][0];
      pos=r2[a][1];
      fprintf ( stdout, "SEQ %5d %5d %5d %s ",r2[a][1]+1,r2[a][3], r2[a][4], (do_tm)?alnpos2hmmtop_pred (A, NULL, r2[a][1], SHORT):"x");
      for (c=0; c<A->nseq; c++)fprintf (stdout, "%c", A->seq_al[c][r2[a][1]]);



      bad=0;
      for (c=0; c< A->nseq; c++)
	{
	  int activity, res;

	  activity=SAR->seq_al[comp][c];
	  res=A->seq_al[c][pos];

	  if (activity=='O')aa[0][res]++;
	  if (activity=='I')aa[1][res]++;
	}

      for (c=0; c< A->nseq; c++)
	{
	  int activity, res;
	  activity=SAR->seq_al[comp][c];
	  res=A->seq_al[c][pos];
	  bad+=(aa[0][res] && aa[1][res])?1:0;
	  aa[0][res]=aa[1][res]=0;
	}
      fprintf ( stdout, " %20s %d |\nCOM %5d %5d %5d %s %s %20s %d |\n\n", SAR->name[r2[a][0]],bad,r2[a][1]+1,r2[a][3],r2[a][4], (do_tm)?alnpos2hmmtop_pred (A, NULL, r2[a][1], SHORT):"x",SAR->seq_al[r2[a][0]], SAR->name[r2[a][0]], bad);


      if (r2[a][4]>threshold)
	{
	  cons[r2[a][1]]++;
	  r3[r2[a][1]][1]++;
	  r3[r2[a][1]][2]+=r2[a][3];
	  r4[r2[a][1]][r2[a][0]]=1;
	}
    }
  sort_int (r3, 3,1,0, A->len_aln-1);

  for (a=0; a<A->len_aln; a++)
    {
      if (r3[a][1]>0)
	{
	  fprintf ( stdout, "\nPOS %4d %4d %4d %c ", r3[a][0]+1, r3[a][1], r3[a][2], cons[r3[a][0]]);
	  for (b=0; b<SAR->nseq; b++)fprintf ( stdout, "%d", r4[r3[a][0]][b]);
	  if (do_tm)fprintf ( stdout, " %s",alnpos2hmmtop_pred (A, NULL, r3[a][0], VERBOSE));
      }
    }
  for (a=0; a< A->nseq; a++)fprintf ( stdout, "\n#MSA >%s\n#MSA %s",A->name[a], A->seq_al[a]);
  fprintf ( stdout, "\n#MSA >cons\n#MSA %s", cons);

  return 0;
}
int ***    simple_sar_predict     (Alignment *A, Alignment *SAR, char *mode)
{
  //This function estimates the z score of every poition with every compound
  //The best Z-score position is then used for the prediction

  int a, b, c, nts, pos, Rscore,Zscore;
  int ***r;
  int ***pred;
  int **aa;


  aa=declare_int (2,256);
  pred=(int***)declare_arrayN (3, sizeof (int),SAR->nseq, A->nseq, 5);


  r=simple_sar_analyze_col (A, SAR, mode);
  nts=SAR->len_aln; //number_trainning_sequences;

  for (a=0; a<SAR->nseq; a++)
    {
      sort_int (r[a],4, 2, 0, A->len_aln-1);

      pos=r[a][A->len_aln-1][3]; //Best Position
      Zscore=r[a][A->len_aln-1][2]; //Best Z-Score
      Rscore=r[a][A->len_aln-1][0]; //Best Z-Score


      for (c=0; c<nts; c++)
	{

	  if (SAR->seq_al[a][c]=='I')aa[1][(int)A->seq_al[c][pos]]++;//Build Positive Alphabet for Compound a
	  if (SAR->seq_al[a][c]=='O')aa[0][(int)A->seq_al[c][pos]]++;//Build Positive Alphabet for Compound a
	}
      for (c=nts; c<A->nseq; c++)
	{
	  pred[a][c][0]=pos;
	  pred[a][c][1]=Zscore;
	  pred[a][c][2]=Rscore;
	  if (aa[1][(int)A->seq_al[c][pos]]>0)
	    {

	      pred[a][c][3]=aa[1][(int)A->seq_al[c][pos]];
	      pred[a][c][4]=aa[0][(int)A->seq_al[c][pos]];
	    }
	}
      for (c=0; c<nts; c++)aa[0][(int)A->seq_al[c][pos]]=aa[1][(int)A->seq_al[c][pos]]=0;
    }

  for ( a=nts; a< A->nseq; a++)
    {
      for ( b=0; b<SAR->nseq; b++)
	{
	  fprintf ( stdout, ">%-25s %-25s Pos %3d ZScore %3d Rscore %3d Activity +: %d -: %d ", A->name [a], SAR->name[b], pred[b][a][0],pred[b][a][1], pred[b][a][2], pred[b][a][3], pred[b][a][4]);
	  if (pred[b][a][4]==0)for (c=0; c<pred[b][a][3]; c++)fprintf ( stdout, "*");
	  fprintf ( stdout, "\n" );
	  for (c=0; c<A->nseq; c++)fprintf ( stdout, "%c", A->seq_al[c][ pred[b][a][0]]);
	  fprintf ( stdout, " %s\n", SAR->name[b]);
	  for (c=0; c<A->nseq-1; c++)fprintf ( stdout, "%c", SAR->seq_al[b][c]);
	  fprintf ( stdout, "  %s\n", SAR->name[b]);	  fprintf ( stdout, "\n");

	}
    }
  return pred;
}
int *pair_seq2seq (int *iseq, char *seq1, char *seq2);
int **simple_sar_analyze_pair_col ( Alignment *inA, Alignment *SAR, char *mode)
{

  int a, b, c, n, n2;
  int *iseq=NULL;
  static int **result, **fresult;
  int sar_mode=1;
  int maxgapratio=0;
  int nresults=10;
  double sum, sum2, score;
  Alignment *A;
  char aa;

  if (!result)
    {
      result=declare_int (inA->len_aln*inA->len_aln,5);

      fresult=declare_int (inA->len_aln*nresults*SAR->nseq, 5);

    }

  A=rotate_aln (inA, NULL);


  for (n2=0,a=0; a<SAR->nseq; a++)
    {

      for (n=0, sum=0, sum2=0,b=0; b<A->nseq-1; b++)
	{
	  for ( c=b+1; c<A->nseq; c++, n++)
	    {

	      iseq=pair_seq2seq (iseq,A->seq_al[b], A->seq_al[c]);
	      if ( sar_mode==1)
		score=sar_vs_iseq1(SAR->seq_al[a],iseq,maxgapratio,NULL,&aa);
	      else if (sar_mode==4)
		score=sar_vs_iseq4(SAR->seq_al[a],iseq,maxgapratio,NULL,&aa);
	      //HERE ("%d", (int)score);
	      result[n][0]=a;//compound;
	      result[n][1]=b; //pos1
	      result[n][2]=c; //pos2
	      result[n][3]=(int)score;

	      sum+=score;
	      sum2+=score*score;
	    }
	}
      for (b=0; b<n; b++)
	result[c][4]=(int)10*return_z_score ((double)result[c][3], sum, sum2,n); //ZScore
      sort_int (result,5,3,0,n-1);
      for ( b=n-nresults; b<n; b++, n2++)
	{

	  for ( c=0; c<5; c++)fresult[n2][c]=result[b][c];
	}
    }
  sort_int (fresult,5,3,0,n2-1);
  fresult[n2][0]=-1; //Terminator;

  return fresult;
}

int *pair_seq2seq (int *iseq, char *seq1, char *seq2)
{
  static int **lu;
  int a;
  if (!lu)
    {
      int a, b, n=0;
      lu=declare_int (256,256);
      for ( n=a=0; a<256; a++)
	for (b=0; b<256; b++)
	  lu[a][b]=n++;
    }

  if (!iseq)iseq=(int*)vcalloc ( strlen (seq1)+1, sizeof (int));
  for ( a=0; a< strlen (seq1); a++)
    {
      if (is_gap(seq1[a]) || is_gap(seq2[a]))iseq[a]=-1;
      else iseq[a]=lu[(int)seq1[a]][(int)seq2[a]];
    }

  return iseq;
}
int ***simple_sar_analyze_col ( Alignment *inA, Alignment *SAR, char *mode)
{
  Alignment *A;
  double score=0, best_score=0;
  int best_pos=0;
  int a, b;

  static int ***result;
  int **sim;
  char aa;
  int sar_mode=3;
  int maxgapratio=0;
  double sum, sum2;


  strget_param (mode, "_METHOD_", "3", "%d", &sar_mode);
  if (!result)
    result=(int***)declare_arrayN (3,sizeof (int),SAR->nseq, inA->len_aln,4);


  sim=aln2sim_mat (inA, "idmat");
  A=rotate_aln (inA, NULL);


  for ( a=0; a<SAR->nseq; a++)
    {
      best_pos=best_score=0;
      for ( sum=0, sum2=0,b=0; b<A->nseq; b++)
	{

	  if ( sar_mode==1)
	    score=sar_vs_seq1(SAR->seq_al[a], A->seq_al[b],maxgapratio, sim, &aa);
	  else if ( sar_mode==2)
	    score=sar_vs_seq2(SAR->seq_al[a], A->seq_al[b],maxgapratio, sim, &aa);
	  else if (sar_mode ==3)
	    score=sar_vs_seq3(SAR->seq_al[a], A->seq_al[b],maxgapratio, sim, &aa);
	  else if (sar_mode ==4)
	    score=sar_vs_seq4(SAR->seq_al[a], A->seq_al[b],maxgapratio, sim, &aa);
	  else if (sar_mode ==5)
	    score=sar_vs_seq5(SAR->seq_al[a], A->seq_al[b],maxgapratio, sim, &aa);


	  result[a][b][0]+=score*10;
	  result[a][b][1]=aa;
	  result[a][b][3]=b;
	  sum+=result[a][b][0];
	  sum2+=result[a][b][0]*result[a][b][0];

	}
      for ( b=0; b< A->nseq; b++)result[a][b][2]=10*return_z_score ((double)result[a][b][0], sum, sum2, A->nseq); //Score
    }

  return result;

 }
int *seq2iseq ( char *seq);
double sar_vs_seq4 ( char *sar, char *seq, float gl, int **sim, char *best_aa)
{

  return sar_vs_iseq4 (sar, seq2iseq(seq), gl, sim, best_aa);
}
double sar_vs_seq1 ( char *sar, char *seq, float gl, int **sim, char *best_aa)
{

  return sar_vs_iseq1 (sar, seq2iseq(seq), gl, sim, best_aa);
}

int *seq2iseq ( char *seq)
{
 static int *iseq, clen;
  int a;

  if (!iseq || clen<strlen (seq))
    {
      clen=strlen (seq);
      vfree (iseq);
      iseq=(int*)vcalloc (clen+1, sizeof (int));
    }
  for (a=0; a<strlen (seq); a++)
    {
      iseq[a]=seq[a];
    }
  return iseq;
}

double sar_vs_iseq1 ( char *sar, int *seq, float gl, int **sim, char *best_aa)
{
  double score=0, return_score=0;
  int RN,N11, Nmsa, Nsar, N, N10, N01, N00;
  int a, b, r, s, res, res1;
  double Ng=0;
  static int *aa, *aal, naa;

  /*measure the E-Value for every amino acid. Returns the best one*/



  N=strlen (sar);
  for (a=0; a<N; a++)
    Ng+=(is_gap(seq[a])|| seq[a]==-1);
  Ng/=N;

  if (Ng>gl) return 0;

  if (!aa)
    {
      aa=(int*)vcalloc (256*256, sizeof(int));
      aal=(int*)vcalloc (N, sizeof (int));
    }
  naa=0;
  for ( a=0; a<N; a++)
    {
      if (!aa[seq[a]])
	{
	  aal[naa++]=seq[a];
	  aa[seq[a]]=1;
	}
    }

  best_aa[0]='-';
  for (a=0; a<naa; a++)
    {

      res=aal[a];

      if (res==-1 || !aa[res]);
      else
	{
	  RN=Nmsa=Nsar=N11=N10=N01=N00=0;

	  for (b=0; b<N; b++)
	    {

	      res1=seq[b];
	      if (res1=='-' || res1==-1)r=0;
	      else r=(res1==res)?1:0;


	      s=(sar[b]=='I')?1:0;

	      Nmsa+=r; Nsar+=s;
	      N11+=(r && s)?1:0;
	      N01+=(!r &&s)?1:0;
	      N10+=(r && !s)?1:0;
	      N00+=(!r && !s)?1:0;
	      RN++;
	    }

	  if (N11)
	    {
	      score=evaluate_sar_score1 ( RN, N11, Nmsa, Nsar);

	    }
	  else
	    {
	      score=0;
	    }

	  if ( score>return_score)
	    {
	      best_aa[0]=res;
	      return_score=score;
	    }
	}
    }
  for ( a=0; a<N; a++)aa[seq[a]]=0;

  return return_score;
}
double sar_vs_seq5 ( char *sar, char *seq, float gl, int **sim, char *best_aa)
{
  double score=0;
  int N11, Nmsa, Nsar, N, N10, N01, N00;
  int a, b, r, s;
  double Ng=0;
  int *aa;

  //measure the E-Value if all the 1AA are considered like alphabet 1*/
  //This function is testing ONE Compound (sar) against one column of MSA (seq) containing N sequences
  //Returns 0 if the active and inactive alphabet overlap
  N=strlen (sar);
  //N is the nmber of sequences

  for (a=0; a<N; a++)
    Ng+=is_gap(seq[a]);
  Ng/=N;

  if (Ng>gl) return 0;

  //Identify all the AA associated with a I (Positive alphabet)
  aa=(int*)vcalloc ( 256, sizeof (int));
  for (b=0; b<N; b++)
    {
      s=(sar[b]=='I')?1:0;
      if (s)aa[(int)seq[b]]=1;
    }
  N11=N10=N01=N00=Nmsa=Nsar=0;


  for (b=0; b<N; b++)
    {

      r=aa[(int)seq[b]];
      s=(sar[b]=='I')?1:0;

      Nmsa+=r; Nsar+=s;
      N11+=(r && s)?1:0;
      N01+=(!r &&s)?1:0;
      N10+=(r && !s)?1:0;
      N00+=(!r && !s)?1:0;
    }
  if (N10>=1 || N01>=1) return 0;
  if (N11)
    {
      score=evaluate_sar_score1 ( N, N11, Nmsa, Nsar);
    }
  else score=0;

  vfree (aa);
  return score;

}



double sar_vs_iseq4 ( char *sar, int *seq, float gl, int **sim, char *best_aa)
{
  int N, Ni, No;
  int a, b,c, r, s;
  double Ng=0;
  static int **aa;

  /*Correlation between AA conservation and Activity*/

  N=strlen (sar);
  for (a=0; a<N; a++)
    Ng+=(is_gap(seq[a]) || seq[a]==-1)?1:0;
  Ng/=N;
  if (gl<1)Ng*=100;

  if (Ng>gl) return 0;


  if (!aa)aa=declare_int(2,257*257);
  for (No=Ni=b=0; b<N; b++)
    {

      s=(sar[b]=='I')?1:0;
      if (s){aa[1][seq[b]]++;Ni++;}
      else  {aa[0][seq[b]]++;No++;}
     }
  for (r=0,b=0; b<N; b++)
    {
      if (aa[1][(int) seq[b]]==Ni)
	{
	  r=(No-aa[0][(int)seq [b]])*100/No;
	  break;
	}

    }
  for (c=0; c<N; c++)aa[0][seq[c]]=aa[1][seq[c]]=0;

  return r;
}

double sar_vs_seq3 ( char *sar, char *seq, float gl, int **sim, char *best_aa)
{
  double score=0;
  int N11, Nmsa, Nsar, N, N10, N01, N00;
  int a, b, r, s;
  double Ng=0;
  int *aa;

  //measure the E-Value if all the 1AA are considered like alphabet 1*/
  //This function is testing ONE Compound (sar) against one column of MSA (seq) containing N sequences

  N=strlen (sar);
  //N is the nmber of sequences

  for (a=0; a<N; a++)
    Ng+=is_gap(seq[a]);
  Ng/=N;

  if (Ng>gl) return 0;

  //Identify all the AA associated with a I (Positive alphabet)
  aa=(int*)vcalloc ( 256, sizeof (int));
  for (b=0; b<N; b++)
    {
      s=(sar[b]=='I')?1:0;
      if (s)aa[(int)seq[b]]=1;
    }
  N11=N10=N01=N00=Nmsa=Nsar=0;


  for (b=0; b<N; b++)
    {

      r=aa[(int)seq[b]];
      s=(sar[b]=='I')?1:0;

      Nmsa+=r; Nsar+=s;
      N11+=(r && s)?1:0;
      N01+=(!r &&s)?1:0;
      N10+=(r && !s)?1:0;
      N00+=(!r && !s)?1:0;
    }

  if (N11)
    {
      score=evaluate_sar_score1 ( N, N11, Nmsa, Nsar);
    }
  else score=0;

  vfree (aa);
  return score;

}

double sar_vs_seq2 ( char *sar, char *seq, float gl, int **sim_mat, char *best_aa)
{
  double score=0, return_score=0;
  int L,N11, Nmsa, Nsar,N10, N01, N;
  int a, b,c,d, r1, s1,r2, s2, res;
  double Ng=0;
  int sim, diff, w;
  char string[5];

  /*Weighted E-Value Similarity*/
  L=strlen (sar);
  for (a=0; a<L; a++)
    Ng+=is_gap(seq[a]);
  Ng/=L;

  if (Ng>gl) return 0;
  for (a=0; a<26; a++)
    {

      N=Nmsa=Nsar=N11=N10=N01=0;
      res='a'+a;
      for (d=0,b=0; b<L; b++)d+=((tolower(seq[b]))==res)?1:0;
      if ( d==0) continue;

      for (b=0; b<L; b++)
	{
	  r1=(tolower(seq[b])==res)?1:0;
	  s1=(sar[b]=='I')?1:0;
	  for ( c=0; c<L; c++)
	    {
	      r2=(tolower(seq[c])==res)?1:0;
	      s2=(sar[c]=='I')?1:0;

	      sprintf ( string, "%d%d%d%d", r1,s1, r2, s2);
	      sim= sim_mat[b][c]/10;
	      diff=10-sim;

	      if (strm (string, "0000"))      {w=diff;N+=2*w;}
	      else if ( strm (string, "0011")){w=sim ;N+=2*w ; N11+=w  ;N10+=0   ;N01+=w   ;Nmsa+=w   ;Nsar+=w;}
	      else if ( strm (string, "1010")){w=diff;N+=2*w ; N11+=0  ;N10+=2*w ;N01+=0   ;Nmsa+=2*w ;Nsar+=0;}
	      else if ( strm (string, "0101")){w=diff;N+=2*w;  N11+=0  ;N10+=0   ;N01+=2*w ;Nmsa+=0   ;Nsar+=2*w;}
	      else if ( strm (string, "1111")){w=diff;N+=2*w;  N11+=2*w;N10+=0   ;N01+=0   ;Nmsa+=2*w ;Nsar+=2*w;}
	      else if ( strm (string, "1001")){w=sim; N+=2*w;  N11+=0  ;N10+=w   ;N01+=w   ;Nmsa+=w;Nsar+=w;}
	      else if ( strm (string, "0110")){w=sim; N+=2*w;  N11+=0  ;N10+=w   ;N01+=w   ;Nmsa+=w;Nsar+=w;}
	    }
	}
      if (N11)
	{

	  score=evaluate_sar_score1 ( N, N11, Nmsa, Nsar);
	}
      return_score=MAX(return_score, score);
    }
  if ( return_score <0)fprintf ( stderr, "\n%.2f", return_score);
  return return_score;
}

float get_sar_sim (char *seq1, char *seq2)
{
  int a, l, s, r;
  int n11=0, n10=0, n01=0, n00=0;


  l=strlen (seq1);
  for ( a=0; a<l; a++)
    {
      s=(seq1[a]=='O')?0:1;
      r=(seq2[a]=='O')?0:1;

      n00+=(!s && !r)?1:0;
      n11+=(s && r)?1:0;
      n01+=(!s && r)?1:0;
      n10+=(s && !r)?1:0;
    }
  if ( n11==0) return 0;
  else return ((float)(n11)*100)/(float)(n11+n10+n01);
}
float get_sar_sim2 (char *seq1, char *seq2)
{
  int a, l, s, r;
  int n11=0, n10=0, n01=0, n00=0, Ns=0, Nr=0;
  float score1;

  l=strlen (seq1);
  for ( a=0; a<l; a++)
    {
      s=(seq1[a]=='O')?0:1;
      r=(seq2[a]=='O')?0:1;

      Ns+=s;
      Nr+=r;

      n00+=(!s && !r)?1:0;
      n11+=(s && r)?1:0;
      n01+=(!s && r)?1:0;
      n10+=(s && !r)?1:0;
    }
  if ( n11==0) return 0;




  score1=evaluate_sar_score1 (l, n11, Ns, Nr);

  return score1;
}
float sar_aln2cor (Alignment *A);
int   sar_aln2ev  (Alignment *A);

int sarseq2ev ( char *s1, char *s2, int MODE);
char*  sarseq2anti_sarseq (char *seq_in, char *seq_out);
Sequence * compare_sar_sequence( Sequence *S1, Sequence *S2, int depth)
{
  int a, b, c, d, e, f;
  char *s1, *s2, *s3,*s4,*s5, *n1, *n2, *n3, *n4, *n5;
  int max=0;;
  int max_depth=5;
  int score;
  char *a1, *a0, *a2, *a3, *a4;
  Alignment *A;

  if ( depth>max_depth)
    {
      printf_exit (EXIT_FAILURE, stderr,"maximum depth: %d", max_depth);
    }
  if ( depth==0) depth=2;
  A=declare_aln2 (strlen (S1->seq[0]),depth);
  a0=A->seq_al[0];
  a1=A->seq_al[1];
  A->len_aln=strlen (S1->seq[0]);

  for (a=0; a< S1->nseq; a++)
    for ( b=0; b<S2->nseq; b++)
      {
	A->nseq=2;
	sprintf (a0, "%s", S1->seq[a]);
	sprintf (a1, "%s", S2->seq[b]);

	if ( strlen (a0)!=strlen (a1))
	  {
	    add_warning (stderr, "WARNING %s (%d) and %s (%d) do not have the same length", S1->name[a], strlen (S1->seq[a]), S2->name[b], strlen (S2->seq[b]));
	    myexit (EXIT_FAILURE);
	  }

	fprintf ( stdout, ">2 %15s %15s CORR: %.3f EVAL: %5d\n",S1->name[a], S2->name[b],  sar_aln2cor (A), sar_aln2ev (A));
	sarseq2anti_sarseq (S1->seq[a],a0);
	fprintf ( stdout, ">2 %15s %15s ANTI: %.3f EVAL: %5d\n", S1->name[a], S2->name[b], sar_aln2cor (A), sar_aln2ev (A));
	if ( depth >=3)
	  {
	    A->nseq=3;
	    a2=A->seq_al[2];
	    for (c=b+1; c<S2->nseq; c++)
	      {
		sprintf (a0, "%s", S1->seq[a]);
		sprintf (a1, "%s", S2->seq[b]);
		sprintf (a2, "%s", S2->seq[c]);
		fprintf ( stdout, ">2 %15s %15s %15s CORR: %.3f EVAL: %5d\n",S1->name[a], S2->name[b],S2->name[c],  sar_aln2cor (A), sar_aln2ev (A));
		sarseq2anti_sarseq (S1->seq[a],a0);
		fprintf ( stdout, ">2 %15s %15s %15s ANTI: %.3f EVAL: %5d\n", S1->name[a], S2->name[b], S2->name[c],sar_aln2cor (A), sar_aln2ev (A));
	      }
	    if ( depth>=4)
	      {
		A->nseq=4;
		a3=A->seq_al[2];
		for (d=c+1; d<S2->nseq; d++)
		  {
		    sprintf (a0, "%s", S1->seq[a]);
		    sprintf (a1, "%s", S2->seq[b]);
		    sprintf (a2, "%s", S2->seq[c]);
		    sprintf (a3, "%s", S2->seq[d]);

		    fprintf ( stdout, ">2 %15s %15s %15s %15s CORR: %.3f EVAL: %5d\n",S1->name[a], S2->name[b],S2->name[c],S2->name[d],  sar_aln2cor (A), sar_aln2ev (A));
		    sarseq2anti_sarseq (S1->seq[a],a0);
		    fprintf ( stdout, ">2 %15s %15s %15s %15s ANTI: %.3f EVAL: %5d\n", S1->name[a], S2->name[b], S2->name[c],S2->name[d],sar_aln2cor (A), sar_aln2ev (A));
		  }
		if (depth>=5)
		  {
		    A->nseq=5;
		    a4=A->seq_al[3];
		    for (e=d+1; e<S2->nseq; e++)
		      {
			sprintf (a0, "%s", S1->seq[a]);
			sprintf (a1, "%s", S2->seq[b]);
			sprintf (a2, "%s", S2->seq[c]);
			sprintf (a3, "%s", S2->seq[d]);
			sprintf (a4, "%s", S2->seq[d]);
		      }
		  }
	      }
	  }
      }

  return S1;
}
char*  sarseq2anti_sarseq (char *seq_in, char *seq_out)
{
  int a;
  if (!seq_out)seq_out=(char*)vcalloc (strlen (seq_in)+1, sizeof (char));
  for (a=0; a<strlen (seq_in); a++)seq_out[a]=(seq_in[a]=='I')?'O':'I';
  return seq_out;
}
float sar_aln2cor (Alignment *A)
{
  float n1, n11, tot_cor;
  int a, b, c,n;

  tot_cor=0;
  for (n=0,a=0; a<A->nseq-1; a++)
    for (b=a+1; b<A->nseq; b++)
      {
	for (n11=n1=0,c=0; c<A->len_aln; c++)
	  {
	    n11+=(A->seq_al[a][c]=='I' && A->seq_al[b][c]=='I');
	    n1+= (A->seq_al[a][c]=='I' || A->seq_al[b][c]=='I');
	  }
	tot_cor+=(n1==0)?0:n11/n1;
	n++;
      }
  tot_cor/=n;
  return tot_cor;
}
int sarseq_pair2ev  ( char *s1, char *s2,int mode);
int sar_aln2ev (Alignment *A)
{
  float n1, n11;
  int a, b, c, tot=0, n=0;

  tot=0;
  for (a=0; a<A->nseq-1; a++)
    for (b=a+1; b<A->nseq; b++)
      {
	tot+=sarseq_pair2ev (A->seq_al[a], A->seq_al[b], 1);
	n++;
      }
  return tot;
}
int sarseq_pair2ev  ( char *s1, char *s2,int mode)
{
  int l, t1, t2, t11,a, n1, n2, s;
  if ( (l=strlen (s1))!=strlen (s2))
    {
      return -1;
    }
  if (mode==2)
    {
      t1=l/2;
      t2=l/2;
      t11=l/2;
    }
  else
    {
      for (t1=t2=t11=0,a=0; a<l; a++)
	{
	  if ( mode==1)//CORRELATED
	    {
	      n1=(s1[a]=='I')?1:0;
	      n2=(s2[a]=='I')?1:0;
	    }
	  else if ( mode==0)//ANTICORRELATED
	    {
	      n1=(s1[a]=='I')?1:0;
	      n2=(s2[a]=='O')?1:0;
	    }
	  t1+=n1;
	  t2+=n2;

	  t11+=(n1&&n2)?1:0;
	}
    }


  s=(int)(100*evaluate_sar_score1 (l, t11, t1, t2));
  return s;
}

double evaluate_sar_score1 ( int N, int n11, int n1msa, int n1sar)
{
  double p, p1, p2;
  int n10, n01;

  n10=n1msa-n11;
  n01=n1sar-n11;

  if ( n11==0)return 0;
  //  if ( (n10)>n11 || n01>n11)return 0;

  p1= M_chooses_Nlog (n1msa, N) + M_chooses_Nlog (n1sar-n11, N-n1msa) + M_chooses_Nlog (n11, n1msa);
  p2=(M_chooses_Nlog (n1msa, N)+  M_chooses_Nlog (n1sar, N));
  p=(p1-p2);

  return -p;

}
double evaluate_sar_score2 ( int N, int n11, int n1msa, int n1sar)
{


  return n11-((n1msa-n11)+(n1sar-n11));

  if ( n11<n1msa) return 0;
  else if ( n11<n1sar) return 0;
  else if ( n11==N)return 0;
  return n11;
}


int benchmark_sar( int value)
{
  static int v[1000];
  static int a;

  if (a==0)
    {
      for (a=0; a< 1000; a++)v[a]=0;
      v[2]=1;
      v[3]=2;
      v[6]=2;
      v[7]=1;
      v[8]=2;
      v[9]=1;
      v[10]=1;
      v[11]=1;
      v[12]=2;
      v[30]=2;
      v[31]=1;
      v[32]=2;
      v[33]=1;
      v[34]=2;
      v[35]=1;
      v[36]=1;
      v[37]=2;
      v[43]=2;
      v[44]=1;
      v[45]=2;
      v[73]=2;
      v[74]=1;
      v[75]=1;
      v[76]=2;
      v[80]=2;
      v[81]=1;
      v[82]=2;
      v[83]=1;
      v[85]=2;
      v[86]=1;
      v[87]=1;
      v[88]=2;
      v[89]=2;
      v[90]=1;
      v[91]=2;
      v[92]=1;
      v[93]=2;
      v[103]=2;
      v[104]=1;
      v[105]=1;
      v[106]=1;
      v[107]=2;
      v[130]=2;
      v[131]=1;
      v[132]=2;
      v[133]=1;
      v[134]=1;
      v[135]=1;
      v[136]=2;
      v[137]=1;
      v[138]=2;
      v[271]=2;
      v[272]=1;
      v[273]=2;
      v[281]=2;
      v[282]=1;
      v[283]=2;
      v[284]=1;
      v[285]=1;
      v[286]=1;
      v[287]=2;
      v[319]=2;
      v[320]=1;
      v[321]=1;
      v[322]=1;
      v[323]=1;
      v[324]=2;
      v[325]=1;
      v[326]=2;
      v[327]=1;
      v[328]=2;
      v[356]=2;
      v[357]=1;
      v[358]=1;
      v[359]=2;
      v[377]=2;
      v[378]=1;
      v[379]=2;
      v[386]=3;
      v[388]=2;
      v[389]=1;
      v[390]=1;
      v[391]=1;
      v[392]=2;
      v[393]=2;
      v[394]=2;
      v[395]=1;
      v[396]=1;
      v[397]=2;
      v[399]=2;
      v[400]=1;
      v[401]=2;
      v[414]=2;
      v[415]=1;
      v[416]=2;
      v[420]=2;
      v[421]=1;
      v[422]=1;
      v[423]=1;
      v[424]=2;
      v[425]=1;
      v[426]=2;
    }
  return v[value];
}

Alignment *weight2sar (Alignment *A, Alignment *SAR, char *weight_file, int limit)
{
  int a, b, c;
  int ***weight;
  char ***list;
  float score;

  weight=(int***)vcalloc (SAR->nseq, sizeof (int**));


  list=file2list (weight_file, " ");

  a=b=0;
  for (a=0; a< SAR->nseq; a++)
    {
      b=c=0;
      while (list[b])
	{
	  if ( strm (list[b][1], SAR->name[a]) && atoi (list[b][3])>0)c++;
	  b++;
	}

      weight[a]=declare_int (c+1, 3);
      fprintf ( stderr, "\n%s %d", SAR->name[a], c);
      b=c=0;
      while (list[b])
	{
	  if ( strm (list[b][1], SAR->name[a]) && atoi (list[b][3])>0)
	    {
	      weight[a][c][0]=atoi(list[b][2])-1;
	      weight[a][c][1]=list[b][5][0];
	      weight[a][c][2]=atoi (list[b][3]);
	      c++;
	    }
	  b++;
	}
      weight[a][c][0]=-1;
    }

  for (a=0; a<A->nseq; a++)
    {
      fprintf ( stdout, ">%s\n", A->name[a]);
      for ( b=0; b< SAR->nseq; b++)
	{
	  score=seq2weighted_sar_score(A->seq_al[a], weight[b]);
	  fprintf ( stdout, "%c", (score>limit)?'I':'O');
	}
      fprintf (stdout, "\n");
    }
  myexit (EXIT_SUCCESS);
  return A;
}

Alignment *display_sar ( Alignment *A, Alignment *SAR, char *compound)
{
  int a,c;
  char name[100];

  c=name_is_in_hlist ( compound, SAR->name, SAR->nseq);
  if ( c==-1)return A;

  for ( a=0; a< A->nseq; a++)
    {
      sprintf (name, "%s", A->name[a]);
      sprintf ( A->name[a], "%c_%s_%s", SAR->seq_al[c][a], name,compound);
    }
  return A;
}
Alignment *aln2weighted_sar_score ( Alignment *A,Alignment *SAR, char *weight_file, char *compound)
{

  int a, b, c=0;
  int **weight;

  int score;
  char reactivity;
  char ***list;


  if ( SAR)
    {
      c=name_is_in_hlist (compound, SAR->name, SAR->nseq);
    }

  list=file2list (weight_file, " ");
  a=b=0;
  while (list[a])
    {
      if (strm (list[a][1], compound))b++;
      a++;
    }
  weight=declare_int ( b+1, 3);


  a=b=0;
  while (list[a])
    {
      if ( !strm (list[a][1], compound) || strm ("TOTPOS", list[a][1]));
      else
	{
	  weight[b][0]=atoi(list[a][2])-1;
	  weight[b][1]=list[a][5][0];
	  weight[b][2]=atoi(list[a][3]);
	  b++;
	}
      a++;
    }
  weight[b][0]=-1;
  for ( a=0; a< A->nseq; a++)
    {
      score=seq2weighted_sar_score (A->seq_al[a], weight);
      reactivity=(!SAR || c==-1)?'U':SAR->seq_al[c][a];

      sprintf (A->seq_comment[a], "Compound %-15s Reactivity %c SAR_SCORE %5d", compound,reactivity, (int) score);

    }
  return A;
}

float seq2weighted_sar_score ( char *seq, int **weight)
{
  int a, p, r, w;
  float score=0;

  a=0;
  while (weight[a][0]!=-1)
    {
      p=weight[a][0];
      r=weight[a][1];
      w=weight[a][2];

      if ( is_gap(seq[p]));
      else if ( tolower(seq[p])==r)score+=w;
      a++;
    }
  return score;
  }

Alignment * sar2simpred (Alignment *A, Alignment *SAR, char *posfile, char *compound, int L1,int L2 )
{
  int a, b, c, c1, c2;
  int **sim, **sim_ref, npred=0;
  float n11, n10, n01, n00;
  float sn, sp;

  int tot_sim=0;
  int N11=1, N01=2, N10=3, NXX=4, SIM=5;
  float ***tot;
  int i1, i2;


  n11=n10=n01=n00=0;
  tot=(float***)declare_arrayN(3,sizeof (float), 10, 6, 2);

  sim_ref=aln2sim_mat (A, "idmat");
  if (strm (posfile, "all"))
    sim=sim_ref;
  else
    {
      Alignment *B;
      B=copy_aln ( A,NULL);
      B=extract_aln3(B,posfile);

      /*if (B->len_aln==0)L1=100;
      else
	L1=((B->len_aln-1)*100)/B->len_aln;

      if (L1<=0)L1=100;
      */
      sim=aln2sim_mat (B, "idmat");
    }

  for (a=0; a< A->nseq-1; a++)
    {
      for ( b=a+1; b< A->nseq; b++)
	{
	  for ( c=0; c<SAR->nseq; c++)
	    {
	      if ( (strm (compound, SAR->name[c]) || strm ( compound, "all")))
		{
		  /*if ( sim_ref[a][b]<30 || sim_ref[a][b]>60)continue;*/
		  i1=0; /*sim_ref[a][b]/10;if (i1==10)i1--;*/

		  i2=sim[a][b];


		  c1=(SAR->seq_al[c][a]=='I')?1:0;
		  c2=(SAR->seq_al[c][b]=='I')?1:0;

		  n11=(c1 && c2)?1:0;
		  n01=(!c1 && c2)?1:0;
		  n10=(c1 && !c2)?1:0;
		  n00=(!c1 && !c2)?1:0;

		  tot[i1][N11][0]+=n11;
		  tot[i1][N01][0]+=n01;
		  tot[i1][N10][0]+=n10;
		  /*tot[i1][N00][0]+=n00;*/
		  tot[i1][NXX][0]++;
		  tot[i1][SIM][0]+=sim_ref[a][b];

		  if ( i2>=L1)
		    {
		      tot[i1][N11][1]+=n11;
		      tot[i1][N01][1]+=n01;
		      tot[i1][N10][1]+=n10;
		      /*tot[i1][N00][1]+=n00;*/
		      tot[i1][NXX][1]++;
		      tot[i1][SIM][1]+=sim_ref[a][b];
		    }
		}
	    }
	}
    }

  for (a=0; a<1; a++)
    {
      sp=(tot[a][N11][0])/(tot[a][N11][0]+tot[a][N10][0]);
      fprintf ( stdout, "\n%15s N11 %5d SP %.2f ",compound, (int)tot[a][N11][0],sp);
      sp=((tot[a][N11][1]+tot[a][N10][1])==0)?1:(tot[a][N11][1])/(tot[a][N11][1]+tot[a][N10][1]);
      sn=(tot[a][N11][0]==0)?1:(tot[a][N11][1]/tot[a][N11][0]);
      fprintf ( stdout, " N11 %5d SP %.2f SN %.2f SIM %.2f", (int)tot[a][N11][1], sp,sn, (tot[a][SIM][1]/tot[a][NXX][1]));
    }

  myexit (EXIT_FAILURE);
  sp=((n11+n01)==0)?1:n11/(n11+n01);
  sn=((n11+n01)==0)?1:n11/(n11+n10);

  fprintf ( stdout, "\nLimit: %d NPRED %d AVGSIM %d SN %.2f   SP %.2f TP %d FP %d FN %d",L1, npred, tot_sim, sn, sp, (int)n11, (int)n01, (int)n10);
  myexit (EXIT_SUCCESS);
  return A;
}

Alignment * sar2simpred2 (Alignment *A, Alignment *SAR, char *seqlist, char *posfile, char *compound, int L )
{
  int a,b, c,c1, c2, p, s;
  float n11, n10, n01, n00, n, sn2, prediction,sp, n1, n0, t, entropy, Delta;
  int *rlist, *tlist, *pred, *npred, tsim, psim;
  int **sim, **sim_ref;
  int nr=0;
  int nrs;
  char *out;
  int delta_max;
  Alignment *B;
  int printall=1;

  out=(char*)vcalloc (A->nseq+1, sizeof (char));
  rlist=(int*)vcalloc ( A->nseq, sizeof (int));
  tlist=(int*)vcalloc ( A->nseq, sizeof (int));
  pred=(int*)vcalloc(2, sizeof (int));
  npred=(int*)vcalloc(2, sizeof (int));

  nrs=0;
  if ( strm (seqlist, "first"))
    {
      for ( a=0; a<SAR->nseq; a++)
	{
	  if ( strm ( compound, SAR->name[a]))
	    {
	      for ( b=0; b<A->nseq; b++)
		{
		  if ( SAR->seq_al[a][b]=='I')
		    {
		      fprintf ( stderr, "COMP: %s REF SEQ: %s\n", A->name[b], compound);
		      rlist[nrs]=b;
		      tlist[rlist[nrs]]=1;
		      nrs++;
		      break;
		    }
		}
	    }
	}
    }
  else if (strm (seqlist, "all"))
    {
      for ( a=0; a< A->nseq; a++)
	{
	  rlist[nrs]=a;
	  tlist[rlist[a]]=1;
	  nrs++;
	}
    }
  else if ((a=name_is_in_hlist ( seqlist, A->name, A->nseq))!=-1)
    {
      rlist[nrs]=a;
      tlist[rlist[nrs]]=1;
      nrs++;
    }
  else
    {
      Alignment *R;
      R=main_read_aln (seqlist, NULL);
      for (a=0; a<R->nseq; a++)
	{
	  rlist[a]=name_is_in_hlist( R->name[a], A->name, A->nseq);
	  tlist[rlist[a]]=1;
	}
      free_aln (R);
    }

  c=name_is_in_hlist ( compound, SAR->name, SAR->nseq);

  sim_ref=aln2sim_mat (A, "idmat");
  if (strm (posfile, "all"))
    {
      sim=sim_ref;
      B=A;
    }
  else
    {
      B=copy_aln ( A,NULL);
      B=extract_aln3(B,posfile);
      sim=aln2sim_mat (B, "idmat");
    }

  n11=n10=n01=n00=n=n1=n0=0;
  delta_max=0;
  for (a=0; a<A->nseq; a++)
    {
      if ( tlist[a] && !strm (seqlist, "all"))
	out[a]=(SAR->seq_al[c][a]=='I')?'Z':'z';/*SAR->seq_al[c][a];*/
      else
	{

	  pred[0]=pred[1]=0;
	  npred[0]=npred[1]=1;
	  c1=(SAR->seq_al[c][a]=='I')?1:0;
	  for (nr=0,tsim=0,psim=0,b=0; b<nrs; b++)
	    {
	      if ( SAR->seq_al[c][rlist[b]]=='o');
	      else
		{
		  c2=(SAR->seq_al[c][rlist[b]]=='I')?1:0;
		  nr+=c2;
		  s=sim[a][rlist[b]];
		  tsim+=sim_ref[a][rlist[b]];
		  psim+=sim[a][rlist[b]];
		  if (s>=L)
		    {
		      pred[c2]+=s;
		      npred[c2]++;
		    }
		}
	    }

	  if (c1==0)n0++;
	  else n1++;
	  t++;


	  Delta=pred[1]-pred[0];

	  if (Delta<-delta_max){p=0;out[a]= (c1==0)?'O':'o';}
	  else if (Delta>delta_max){p=1;out[a]=(c1==1)?'I':'i';}
	  else {p=-1; out[a]=(c1==1)?'U':'u';}

	  if ( p==-1);
	  else if (  p &&  c1)n11++;
	  else if (  p && !c1)n10++;
	  else if ( !p && !c1)n00++;
	  else if ( !p &&  c1)n01++;

	  if (p!=-1)n++;
	  if (printall)fprintf ( stdout, ">%-15s %d %c OVERALL_SIM:%d POSITION_SIM %d\n%s\n", B->name[a], c1, out[a],tsim/nrs,psim/nrs,B->seq_al[a]);
	}
    }
  sp=((n11+n10)==0)?1:n11/(n11+n10);
  sn2=((n1)==0)?1:n11/n1;
  prediction=(n11+n00)/(n1+n0);
  entropy=(float)(M_chooses_Nlog (nr, nrs)/M_chooses_Nlog(nrs/2, nrs));

  fprintf ( stdout, ">%-15s Sp %.2f  Sn %.2f Pred %.2f E %.2f\n", compound,sp, sn2,prediction,entropy );
  fprintf ( stdout, "%s\n", out);

  myexit (EXIT_SUCCESS);
  return A;
}
/************************************************************************************/
/*                ALIGNMENT ANALYZE     : SAR FOR OR                                */
/************************************************************************************/

void display_or_help();
void display_or_help()
{
  fprintf ( stdout, "\nor_sar options:");
  fprintf ( stdout, "\n_ORCL_: Command_line in a file");

  fprintf ( stdout, "\n_ROTATE_  : rotate the sar matrix (if each entry is a compound rather than a sequence)");
  fprintf ( stdout, "\n_JNIT_    : number cycles of Jacknife");
  fprintf ( stdout, "\n_JNSEQ_   : Number of sequences picked up in alignment [0 to keep them all, 10 by default]");
  fprintf ( stdout, "\n_JSAR_    : Number of compounds picked up in the SAR matrix [0 to keep them all => default");
  fprintf ( stdout, "\n_JRSAR_   : Randomization of the SAR file between each JNIT iteration: S, C, R, SC, SR... (S: seq, C, column, R: residue");
  fprintf ( stdout, "\n_JRALN_   : Randomization of the ALN file between each JNIT iteration: S, C, R, SC, SR... (S: seq, C, column, R: residue");


  fprintf ( stdout, "\n_NPOS_    : Number of positions used to make the prediction [4 by default]");
  fprintf ( stdout, "\n_DEPTH_   : Depth of the motif degenerated alphabet   [Default=2]");
  fprintf ( stdout, "\n_POSFILE_ : Predefined list of positions in a file<P: pos score");
  fprintf ( stdout, "\n_RSAR_    : Randomization of the SAR file: S, C, R, SC, SR... (S: seq, C, column, R: residue");
  fprintf ( stdout, "\n_RALN_    : Randomization of the SAR file: S, C, R, SC, SR... (S: seq, C, column, R: residue");

  fprintf ( stdout, "\n_SARTHR_  : E-value Threshold for the filtered displays in comploo mode (def: 3)");



  fprintf ( stdout, "\n_MODE_    : jack, loo (leave one out), pos, sim, predict, self");
  fprintf ( stdout, "\n_COMPLIST_: Compound_list provided in a FASTA file. Each compound MUST correspond to the corresponding column in SAR");
  fprintf ( stdout, "\n_OUTALN_  : print the alignment corresponding to the list of positions");
  fprintf ( stdout, "\n_OUTTREE  : Make the tree (default nj) corrwesponding to _OUTALN_");
  fprintf ( stdout, "\n\n");
  myexit (EXIT_SUCCESS);
}
Alignment * simple_sar_or(Alignment *inA, Alignment *inS, char *param);
void sar2jack (Alignment *A, Alignment *S, int nseq, int sarlen);


void sar2jack (Alignment *A, Alignment *S, int nseq, int sarlen)
{
  A=aln2jacknife (A, nseq, 0);
  S=reorder_aln  (S, A->name, A->nseq);
  S=aln2jacknife (S, 0, sarlen);
}
Alignment *or_scan (Alignment *A,Alignment *S, char *pmode)
{
  int l, a,ax,cx, b;
  char mode[100];
  int start, offset,w;
  int nl, *poslist;


  fprintf ( stdout, "\nPARAMETERS: %s\n", pmode);
  fprintf ( stderr, "\nPARAMETERS: %s\n", pmode);

  strget_param (pmode, "_SW_", "15", "%d",&w);
  strget_param (pmode, "_SMINW_", "5", "%d",&start);
  strget_param (pmode, "_SSTART_", "5", "%d",&offset);
  strget_param (pmode, "_SMODE_", "single", "%s",&mode);


  l=intlen (A->len_aln);
  poslist=(int*)vcalloc ( A->len_aln, sizeof (int));
  nl=0;

  for ( a=0; a< A->len_aln; a++)poslist[nl++]=a+1;

  if ( strm (mode, "single"))
    {
      for ( a=0; a<A->len_aln-2; a++)
	{
	  int c, gap, score;
	  for (gap=0,c=0; c<A->nseq; c++)gap+=(A->seq_al[c][a]=='-');
	  if ( !gap)
	    {
	      Alignment *B;
	      B=extract_aln (A, a+1, a+2);
	      B=or_sar (B, S, pmode, NO_PRINT);
	      score=B->score_aln;
	      free_aln (B);
	    }
	  else
	    {
	      score=0;
	    }
	  fprintf ( stdout, "P: %*d  S: %4d\n",l,a+1, score);
	  //fprintf ( stderr, "P: %*d  S: %4d\n",l,a+1, score);

	}
    }
  else if  strm (mode, "scan")
    {

      for ( ax=0; ax<nl; ax++)
	{
	  int best_score=0;
	  int best_pos=0, best_w=0, best_start=0, best_end=0;
	  a=poslist[ax];
	  for (b=start; b<=w; b++)
	    {
	      Alignment *B;
	      int pstart, pend;
	      pstart=a-b;
	      pend=a+b;

	      if (pstart<1 || pstart>A->len_aln)continue;
	      if (pend<1 || pend>A->len_aln)continue;
	      B=extract_aln (A, a-w, a+w);

	      B=or_sar (B, S,pmode, NO_PRINT);


	      if (B->score_aln>=best_score){best_score=B->score_aln; best_pos=a;best_w=b;best_start=pstart; best_end=pend;}
	      free_aln (B);
	    }
	  fprintf ( stdout, "P: %*d  I: %*d %*d Score: %5d L: %2d\n", l,best_pos, l,best_start+offset, l,best_end, best_score,(best_w*2)+1 );
	}
    }
  else if ( strm ( mode, "scan_comp"))
    {
      Alignment *NS=NULL;
      int *tresults;
      int s, n=0, p1;
      int  nbest;
      int *poscache;
      float sc, eval;
      int **index;
      int **poslist;
      int len=1;
      float best_sc, best_eval;
      int best_pos;
      char *best_word;
      int sth;
      int sev;

      index=aln2pos_simple (A, A->nseq);
      strget_param (pmode, "_NPOS_", "1", "%d", &len);
      strget_param (pmode, "_STH_", "0", "%d", &sth);
      strget_param (pmode, "_SEV_", "0", "%d", &sev);

      if (sev==0)for (a=0, sev=1; a<len; a++)sev*=A->len_aln;

      tresults=(int*)vcalloc ( A->len_aln, sizeof (int));
      poscache=(int*)vcalloc ( A->len_aln, sizeof (int));
      poslist=generate_array_int_list (len, 0, A->len_aln-1,1, &n, NULL);
      for (s=0; s<S->len_aln; s++)
	{

	  char *word;
	  NS=aln2block (S, s+1, s+2, NS);
	  fprintf ( stderr, "\nProcess: %s ...\n", get_compound_name(s, mode));

	  for (best_sc=0,best_pos=0,best_word=NULL,a=0; a<n; a++)
	    {
	      for (b=0; b<len; b++)poscache[poslist[a][b]]=1;
	      word=or_id_evaluate2 (A,NS, pmode,poscache, NO_PRINT, &sc);
	      //if (word)HERE ("%s", word);
	      eval=0;
	      for (b=0; b<len; b++)poscache[poslist[a][b]]=0;
	      if (word && sc>best_sc)
		{
		  best_sc=sc;
		  best_eval=eval;
		  best_word=word;
		  best_pos=a;
		  nbest=1;
		}
	      else if ( word && sc>=best_sc)
		{
		  nbest++;
		}
	      else
		vfree (word);
	    }
	  nbest/=len;
	  for (a=0; a<=A->nseq; a++)
	    {



	      fprintf (stdout, "\n>%s PPPP: S: %s SC: %d EV: %d NBest: %d W: %s",  get_compound_name(s, mode),(a==A->nseq)?"cons":A->name[a],(int)best_sc, (int) best_eval, nbest,best_word);
	      for ( b=0; b<len; b++)
		{
		  if ( a==A->nseq)
		    {
		      p1=poslist[best_pos][b];
		      if (best_sc>sth && nbest<sev)tresults[p1]++;
		      p1++;
		    }
		  else
		    {
		      p1=index[a][poslist[best_pos][b]];
		      if (p1>0)p1++;
		    }
		  fprintf ( stdout, " %d ", p1);

		}
	    }

	}
      for ( p1=0; p1<A->len_aln;p1++)
	fprintf ( stdout, "\n>TOT_P: %d %d", p1+1, tresults[p1]);
    }

  else if ( strm ( mode, "scan_comp_old"))
    {
      Alignment *NS=NULL, *BLOCK=NULL;
      int **results, *tresults;
      int s, n, p1, p2;
      int npos=1;
      int *poscache;

      results=declare_int (A->len_aln*A->len_aln, 3);
      tresults=(int*)vcalloc ( A->len_aln, sizeof (int));
      poscache=(int*)vcalloc ( A->len_aln, sizeof (int));
      for (s=0; s<S->len_aln; s++)
	{
	  int count;
	  NS=aln2block (S, s+1, s+2, NS);
	  fprintf ( stderr, "\nProcess: %s ...", get_compound_name(s, mode));
	  for (n=0,p1=0; p1<A->len_aln-w; p1++, count ++)
	    {
	      if ( count == 50){fprintf ( stderr, "*");count=0;}
	      for ( p2=p1; p2<A->len_aln-w; p2++, n++)
		{
		  poscache[p1]=1;
		  poscache[p2]=1;

		  BLOCK=alnpos2block (A,poscache,BLOCK);

		  if ( aln2ngap(BLOCK)<=0)BLOCK=or_sar (BLOCK,NS, pmode, NO_PRINT);
		  else BLOCK->score_aln=0;

		  //if ( BLOCK->score_aln>0)HERE ("P: %d %d %d", p1, p2,BLOCK->score_aln);
		  results[n][0]=p1;
		  results[n][1]=p2;
		  results[n][2]=BLOCK->score_aln;
		  poscache[p1]=poscache[p2]=0;

		}
	    }
	  sort_int_inv (results, 3, 2, 0, n-1);
	  for (p1=0; p1<npos; p1++)
	    {
	      fprintf ( stdout, "\n>%s PPPP: %d %d SC: %d",  get_compound_name(s, mode), results[p1][0]+1,results[p1][1]+1, results[p1][2]);
	      fprintf ( stderr, "\n>%s PPPP: %d %d SC: %d",  get_compound_name(s, mode), results[p1][0]+1,results[p1][1]+1, results[p1][2]);
	      tresults[results[p1][0]]++;
	      tresults[results[p1][1]]++;
	    }
	}
      for ( p1=0; p1<A->len_aln;p1++)
	if ( tresults[p1])fprintf ( stdout, "\n>TOT_P: %d %d", p1+1, tresults[p1]);
    }

  myexit (EXIT_SUCCESS);
}


ORP * or_sar_compound(Alignment *A, Alignment *S, char *mode, int print);
ORP* set_orp_name ( ORP* P, char *name);
ORP* set_orp_offset ( ORP* P,int offset);
Alignment * or_sar(Alignment *inA, Alignment *inS, char *mode, int print)
{
  char rsar[4],raln[4], sarmode[100];
  int rotate, a, b, start, end;
  ORP **R, *ORS;
  Alignment *A, *S;

  if ( mode && (strstr (mode, "help") || strstr (mode, "HELP")))
   display_or_help();

  strget_param (mode, "_SARMODE_", "single", "%s", &sarmode);//single | all
  strget_param (mode, "_RSAR_", "NO", "%s", rsar);
  strget_param (mode, "_RALN_", "NO", "%s", raln);
  strget_param (mode, "_ROTATE_", "0", "%d", &rotate);

  strget_param (mode, "_START_", "1", "%d", &start);start--;
  strget_param (mode, "_END_", "0", "%d", &end);

  A=copy_aln (inA, NULL);
  S=copy_aln (inS, NULL);

  if ( end==0)end=A->len_aln;
  if ( start!=0 || end!=A->len_aln)
    {
      int c;
      if ( start>=A->len_aln || end<=start || end >A->len_aln)
	printf_exit (EXIT_FAILURE, stderr, "ERROR: _START_%d and _END_%d are incompatible with the aln [len=%d][FATAL:%s]", start, end, A->len_aln, PROGRAM);
      for (b=0; b<A->nseq; b++)
	for (c=0,a=start; a<end; a++, c++)
	  A->seq_al[b][c]=A->seq_al[b][a];
      A->len_aln=c;
      for (b=0; b<A->nseq; b++)A->seq_al[b][c]='\0';
    }


  S=aln2random_aln (S, rsar);
  A=aln2random_aln (A, raln);



  if (rotate)
    {
      Alignment *rS;

      if (S->len_aln!=A->nseq)
	printf_exit ( EXIT_FAILURE,stderr, "ERROR: Alignment and SAR matrix are incompatible [FATAL:%s]", PROGRAM);

      rS=rotate_aln (S, NULL);

      for (a=0; a< A->nseq; a++)sprintf (rS->name[a], "%s", A->name[a]);
      free_aln (S);
      S=rS;
    }

  R=(ORP**)vcalloc ( S->len_aln+2, sizeof (ORP*));
  if (strm (sarmode, "all"))R[0]=or_sar_compound (A, S, mode, print);
  else if ( strm (sarmode, "single"))
    {
      for ( a=0; a<S->len_aln; a++)
	{
	  Alignment *NS=NULL;

	  NS=aln2block (S, a+1, a+2, NS);

	  R[a]=or_sar_compound (A, NS, mode, NO_PRINT);
	  set_orp_name ((R[a]), get_compound_name (a, mode));
	  set_orp_offset(R[a], start);
	  display_or_summary (R[a], mode, stdout, print);

	}
    }



  ORS=combine_n_predictions (R, A, S);

  display_or_summary (ORS, mode, stdout, print);
  inA->score_aln=(int)(ORS->best*(float)1000);

  free_orp_list (R);
  free_orp(ORS);

  return inA;
 }
ORP* set_orp_offset ( ORP* P,int offset)
{
  if (!P) return NULL;
  else
    {
      P->offset=offset;
      return set_orp_offset(P->PR, offset);
    }
}

ORP* set_orp_name ( ORP* P, char *name)
{
  if (!P) return NULL;
  else
    {
      sprintf ( P->name, "%s", name);
      return set_orp_name(P->PR, name);
    }
}

ORP * combine_n_predictions (ORP**R, Alignment *A, Alignment *S)
{
  int a=0;
  ORP*N=NULL;
  while (R[a])
    {
      N=combine_2_predictions (R[a++], N, A, S);
    }
  sprintf ( N->name, "ALL");
  sprintf ( N->mode, "COMBINED");

  return N;
}

ORP *combine_2_predictions ( ORP*IN, ORP *TO,Alignment *A, Alignment *S)
{
  int a;

  if ( !TO)
    {
      TO=declare_or_prediction (IN->ncomp, IN->nseq, IN->len);
      TO->A=A;
      TO->S=S;
      TO->P=copy_aln(S, NULL);
      TO->offset=IN->offset;
      TO->ncomp=0;
    }

  for (a=0; a< IN->len; a++)
    {
      TO->pos[a]+=IN->pos[a];
    }

  TO->fp+=IN->fp;
  TO->fn+=IN->fn;
  TO->tp+=IN->tp;
  TO->tn+=IN->tn;
  rates2sensitivity (TO->tp, TO->tn, TO->fp, TO->fn, &(TO->sp), &(TO->sn), &(TO->sen2), &(TO->best));


  for (a=0; a<(TO->A)->nseq; a++)
    {
      (TO->P)->seq_al[a][TO->ncomp]=(IN->P)->seq_al[a][0];
      //(TO->S)->seq_al[a][TO->ncomp]=(IN->S)->seq_al[a][0];
    }
  TO->ncomp++;
  (TO->P)->len_aln=TO->ncomp;
  (TO->A)->score_aln=TO->best;

  return TO;
}



ORP * display_or_summary (ORP *CP, char *mode, FILE *fp, int print)
{
  int a;
  char *pred;
  char *exp;
  char *motif;
  Alignment *A, *P, *S;



  A=CP->A;
  P=CP->P;
  S=CP->S;



  pred=(char*)vcalloc ( P->nseq*P->len_aln*2, sizeof (char));
  exp=(char*)vcalloc ( P->nseq*P->len_aln*2, sizeof (char));
  motif=(char*)vcalloc (CP->len+1, sizeof (char));



  if (P && S)
    {
      for ( a=0; a<P->nseq; a++)
	{
	  strcat (pred,P->seq_al[a]);
	  strcat (exp, S->seq_al[a]);
	}
      CP->evalue=profile2evalue(pred, exp);
    }
  a=0;
  while ( CP->motif && CP->motif[a] && CP->motif[0][a][0])strcat (motif, CP->motif[0][a++]);

  if ( print==PRINT)
    {
      fprintf (fp, "\n>%-10s Mode: %s Accuracy: %6.2f E-value: %6.2f Motif: ",CP->name,CP->mode, CP->best, CP->evalue);

      if (motif[0])
	{
	  fprintf (fp, " %s",motif);
	}
      for ( a=0; a<CP->len; a++)
	{
	  if ( CP->pos[a]) fprintf ( fp, "\n>%-10s Mode: %s P: cons %3d SC: %4d", CP->name, CP->mode, a+1+CP->offset, CP->pos[a]);
	}
      fprintf ( fp, "\n");
    }
  vfree (pred); vfree(motif); vfree(exp);
  if (CP->PR)display_or_summary (CP->PR, mode, fp, print);
  return CP;
}


ORP * or_sar_compound(Alignment *A, Alignment *S, char *mode, int print)
{
  char rmode[100];
  Alignment *P=NULL;
  ORP *PR=NULL;

  strget_param (mode, "_MODE_", "predict", "%s", rmode);



  if (strm (rmode, "jack"))P=or_jack (A, S, mode);
  else if (strm (rmode, "loo")) PR=or_loo (A, S, mode,NULL, print);
  else if (strm (rmode, "comploo")) P=or_comp_loo (A, S, mode,NULL, print);
  else if ( strm (rmode, "comppos")){or_comp_pos ( A, S, mode,print);myexit (EXIT_FAILURE); return NULL;}
  else if ( strm (rmode, "pos"))P=or_aln2pos_aln ( A, S, mode);
  else if ( strm (rmode, "predict"))P=or_predict ( A, S, mode);
  else if ( strm (rmode, "self"))PR=or_self_predict ( A, S, mode, NULL, print);
  else if ( strm (rmode, "sim"))P=or_sim ( A, S, mode);

  else if ( strm (rmode, "test"))P=or_test ( A, S, mode);
  else
    {
       printf_exit (EXIT_FAILURE, stderr, "ERROR: %s is an unknown mode of or_sar  [FATAL:%s]\n", rmode,PROGRAM);
       return NULL;
    }

  if (!PR)
    {
      PR=(ORP*)vcalloc (1, sizeof (ORP));
      PR->P=P;
    }
  return PR;
}
Alignment * or_test ( Alignment *inA, Alignment *inS, char *mode)
{
  return inA;
}



float or_id_evaluate  ( Alignment *A, Alignment *S, char *mode, int *pos, int print)
{
  char *w;
  float score;

  w=or_id_evaluate2(A,S, mode, pos,print, &score);

  vfree (w);
  return score;
}
char* or_id_evaluate2  ( Alignment *A, Alignment *S, char *mode, int *pos, int print, float *rscore)
{
  static char **words;
  static int *plist;
  char *bword;

  int res, p,nl, w, c, s, exp, pred;
  int tp, tn, fp, fn;
  float sn, sp, sen2, best, best_score;

  if (!A)
    {free_char (words, -1);
      vfree (plist);
      return NULL;
    }
  rscore[0]=0;

  plist=pos2list (pos, A->len_aln, &nl);
  words=declare_char (A->nseq, nl+1);
  bword=(char*)vcalloc (nl+1, sizeof (char));
  for (p=0; p<nl; p++)
    {
      for (s=0; s< A->nseq; s++)
	{
	  res=A->seq_al[s][plist[p]];
	  if (res=='-'){or_id_evaluate2 (NULL, NULL, NULL, NULL, 0, NULL);vfree (bword);return 0;}
	  words[s][p]=res;
	}
    }

  for (best_score=0,w=0; w<A->nseq; w++)
    {
      tp=fp=fn=tn=0;

      for (c=0; c<S->len_aln; c++)
	{
	  for (s=0; s<A->nseq; s++)
	    {
	      exp=(S->seq_al[s][c]=='I')?1:0;
	      pred=strm (words[w], words[s]);
	      if      ( exp &&  pred)tp++;
	      else if ( exp && !pred)fn++;
	      else if (!exp && !pred)tn++;
	      else if (!exp &&  pred)fp++;
	    }
	}
      rates2sensitivity (tp, tn, fp,fn,&sp, &sn, &sen2, &best);
      if ( best>best_score)
	{
	  best_score=best;
	  sprintf (bword, "%s", words[w]);
	}
    }
  rscore[0]=(float)1000*best_score;
  or_id_evaluate2 (NULL, NULL, NULL, NULL, 0, NULL);
  return bword;
}

float or_loo_evaluate2 ( Alignment *A, Alignment *S, char *mode, int *pos, int print)
{
  int c, s, p, res, sar;

  int *plist, nl;
  int tp, tn, fp, fn;
  float sn, sp, sen2, best;
  char **words, **positive, **negative;

  tp=tn=fp=fn=0;
  plist=pos2list (pos, A->len_aln, &nl);

  words=declare_char (A->nseq, nl+1);

  for (p=0; p<nl; p++)
    {
      for (s=0; s< A->nseq; s++)
	{
	  res=A->seq_al[s][plist[p]];
	  if (res=='-'){vfree (plist);free_char (words, -1); return 0;}
	  words[s][p]=res;
	}
    }
  positive=(char**)vcalloc ( A->nseq, sizeof (char*));
  negative=(char**)vcalloc ( A->nseq, sizeof (char*));
  for (c=0; c<S->len_aln; c++)
    {
      //Fill the match matrix
      for (p=0; p<nl; p++)
	{
	  for (s=0; s< A->nseq; s++)
	    {
	      sar=S->seq_al[s][c];
	      if (sar=='I')positive[s]=words[s];
	      else if ( sar=='O')negative[s]=words[s];
	    }
	}

      //Evaluate the scores
      for (s=0; s< A->nseq; s++)
	{
	  int pos=0, neg=0, pred;
	  sar=S->seq_al[s][c];
	  positive[s]=negative[s]=NULL;

	  if ( name_is_in_hlist (words[s], positive, A->nseq)!=-1)
	    pos=1;
	  if ( name_is_in_hlist (words[s], negative, A->nseq)!=-1)
	    neg=1;

	  if (pos & !neg) pred=1;
	  else pred=0;

	  if      ( pred  && sar=='I')tp++;
	  else if (!pred  && sar=='I')fn++;
	  else if (!pred  && sar=='O')tn++;
	  else if ( pred  && sar=='O')fp++;

	  if ( sar=='I')positive[s]=words [s];
	  else negative[s]=words[s];
	}
    }

  vfree (negative); vfree (positive);
  vfree (plist); free_char (words, -1);
  rates2sensitivity (tp, tn, fp,fn,&sp, &sn, &sen2, &best);

  return (float)1000*best;
}
float or_loo_evaluate ( Alignment *A, Alignment *S, char *mode, int *pos, int print)
{
  int c, s, p, res, sar;
  int **matP,**matN;
  int *plist, nl;
  int tp, tn, fp, fn;
  float sn, sp, sen2, best;

  tp=tn=fp=fn=0;
  plist=pos2list (pos, A->len_aln, &nl);
  matP=declare_int (nl, 256);
  matN=declare_int (nl, 256);

  for (c=0; c<S->len_aln; c++)
    {
      //Fill the match matrix
      for (p=0; p<nl; p++)
	{
	  for (s=0; s< A->nseq; s++)
	    {
	      res=A->seq_al[s][plist[p]];
	      sar=S->seq_al[s][c];
	      if (res=='-'){vfree (plist); free_int (matP, -1);free_int (matN, -1); return 0;}
	      if (sar=='I')matP[p][res]++;
	      if (sar=='O')matN[p][res]++;
	    }
	}

      //Evaluate the scores
      for (s=0; s< A->nseq; s++)
	{
	  int scoreP, scoreN;
	  int pred, valP, valN;

	  sar=S->seq_al[s][c];
	  for (scoreN=0,scoreP=0,p=0; p<nl; p++)
	    {
	      res=A->seq_al[s][plist[p]];

	      valP=matP[p][res]-(sar=='I')?1:0;
	      scoreP+=(valP>0)?1:0;

	      valN=matN[p][res]-(sar=='O')?1:0;
	      scoreN+=(valN>0)?1:0;
	    }

	  if ( scoreP==nl && scoreN<nl)pred=1;
	  else pred=0;


	  if      ( pred  && sar=='I')tp++;
	  else if (!pred  && sar=='I')fn++;
	  else if (!pred  && sar=='O')tn++;
	  else if ( pred  && sar=='O')fp++;
	}

      //reset the matrix
      for (p=0; p<nl; p++)
	{
	  for (s=0; s< A->nseq; s++)
	    {
	      res=A->seq_al[s][plist[p]];
	      sar=S->seq_al[s][c];
	      if (sar=='I')matP[p][res]=0;
	      else matN[p][res]=0;
	    }
	}
    }

  vfree (plist); free_int (matP, -1);free_int (matN, -1);
  rates2sensitivity (tp, tn, fp,fn,&sp, &sn, &sen2, &best);

  return (float)1000*best;
}
int* or_comp_pos ( Alignment *inA, Alignment *inS, char *mode,int print)
{
  Alignment *A=NULL, *S=NULL, *inS2=NULL;
  int a, b, c;
  int *main_pos, *pos=NULL;

  set_sar (inA, inS, mode);
  main_pos=(int*)vcalloc ( inA->len_aln, sizeof (int));

  inS2=copy_aln (inS, NULL);
  inS2->len_aln=1;



  //Run every SAR, one at a time
  for ( c=0; c< inS->len_aln; c++)
    {
      int max, p;

      fprintf ( stdout, ">%d\n", c);
      for (a=0; a< inS->nseq; a++)
	{
	  inS2->seq_al[a][0]=inS->seq_al[a][c];
	  inS2->seq_al[a][1]='\0';
	}

      vfree (pos);
      free_aln (S);
      free_aln (A);

      pos=(int*)vcalloc (inA->len_aln, sizeof (int));
      A=copy_aln (inA, NULL);
      S=copy_aln (inS2, NULL);
      set_sar (A,S, mode);
      pos=aln2predictive_positions (A, S, mode,PRINT);

      for (max=0,b=0; b<A->len_aln; b++)
	{
	  main_pos[b]+=pos[b];
	  if (main_pos[b]>max)
	    {
	    max=main_pos[b];
	    p=b+1;
	  }
	}


      for (a=0; a<A->nseq; a++)
	{
	  fprintf ( stdout, "\t");
	  for ( b=0; b<A->len_aln; b++)
	    if ( pos[b]) fprintf ( stdout, "%c", A->seq_al[a][b]);
	  fprintf ( stdout, " %c\n", inS2->seq_al[a][0]);
	}
      fprintf ( stdout, "\n\tBest: %d %d\n", p, max);

    }

  if (print==PRINT)
    {
      for ( a=0; a<inA->len_aln; a++)fprintf ( stdout, "\nP2: cons %4d %4d [FINAL]", a+1, main_pos[a]);
    }
  return main_pos;
}
Alignment * or_comp_loo ( Alignment *inA, Alignment *inS, char *mode, int *pos,int print)
{
  int a, b,c, n;
  char **keep, **remove;
  Alignment *A, *S, *P, *P1, *SEQ, *inS2;
  int **main_pos, *compound_pos, **comp_list;
  int pos_exists=0;
  char *comp_pred, *comp_exp;
  int sar_threshold;

  strget_param (mode, "_SARTHRES_", "3", "%d", &sar_threshold);

  if (pos)pos_exists=1;
  set_sar (inA, inS, mode);
  P=copy_aln (inS, NULL);
  keep=declare_char (inA->nseq, MAXNAMES);
  remove=declare_char (inA->nseq, MAXNAMES);

  main_pos=declare_int ( inA->len_aln,4);
  comp_list=declare_int (inA->len_aln, sizeof (int*));
  inS2=copy_aln (inS, NULL);
  inS2->len_aln=1;


  comp_pred=(char*)vcalloc ( inA->nseq+1, sizeof (int));
  comp_exp=(char*)vcalloc ( inA->nseq+1, sizeof (int));
  compound_pos=NULL;
  //Run every SAR, one at a time
  for ( c=0; c< inS->len_aln; c++)
    {
      for (a=0; a< inS->nseq; a++)
	{
	  inS2->seq_al[a][0]=inS->seq_al[a][c];
	  inS2->seq_al[a][1]='\0';
	}
      vfree (compound_pos);
      compound_pos=(int*)vcalloc (inA->len_aln, sizeof (int));
      for (a=0; a<inA->nseq; a++)
	{
	  char ***motifs;

	  A=copy_aln (inA, NULL);
	  S=copy_aln (inS2, NULL);

	  for (n=0,b=0; b<A->nseq; b++)
	    {
	      if (b!=a)sprintf (keep[n++ ], "%s", A->name [b]);
	    }
	  sprintf ( remove[0], "%s", A->name[a]);
	  reorder_aln (A,keep, A->nseq-1);

	  set_sar (A,S, mode);
	  if (!pos_exists)
	    {
	      pos=aln2predictive_positions (A, S, mode,NO_PRINT);
	    }
	  for (b=0; b<A->len_aln; b++)
	    {
	      compound_pos[b]+=pos[b];
	    }

	  motifs=compounds2motifs (A, S, pos,0, mode, NO_PRINT);

	  SEQ=copy_aln (inA, NULL);
	  SEQ=reorder_aln (SEQ, remove, 1);

	  P1=aln2prediction (SEQ, motifs, pos);
	  comp_pred[a]=P1->seq_al[0][0];
	  comp_exp[a]=inS2->seq_al[a][0];
	  P->seq_al[a][c]=P1->seq_al[0][0];


	  free_aln (SEQ);
	  free_aln (S);
	  free_aln (A);
	  free_aln (P1);
	  free_arrayN( (void *)motifs, 3);
	  if (!pos_exists)vfree (pos);
	}
      if (print==PRINT)
	fprintf ( stdout, ">%-15s SC: %.2f E; %.2f\n%s\n%s\n", get_compound_name(c, mode),profile2sensitivity (comp_pred, comp_exp, NULL, NULL, NULL, NULL),profile2evalue(comp_pred, comp_exp),comp_pred, comp_exp);
      for (b=0; b<A->len_aln; b++)
	    {
	      main_pos[b][2]+=compound_pos[b];
	      if (compound_pos[b])main_pos[b][3]++;
	      if (profile2evalue(comp_pred, comp_exp)>sar_threshold)
		{
		  main_pos[b][0]+=compound_pos[b];
		  if (compound_pos[b])
		    {
		      main_pos[b][1]++;
		      comp_list[b][0]++;
		      comp_list[b]=(int*)vrealloc (comp_list[b], sizeof (int)*(comp_list[b][0]+1));
		      comp_list[b][comp_list[b][0]]=c;
		    }
		}
	    }

    }

  P->score_aln=(int)((float)1000*evaluate_prediction (P, inS, mode,print));

  if (print==PRINT)
    {
      for ( a=0; a<inA->len_aln; a++)fprintf ( stdout, "\nP: cons %4d RS: %4d RC: %5d FC: %4d %4d[FINAL]", a+1, main_pos[a][2], main_pos[a][3], main_pos[a][0], main_pos[a][1]);

      for ( a=0; a<inA->len_aln; a++)
	{
	  fprintf ( stdout, "\nP: cons %4d RS: %4d RC: %5d FC: %4d %4d CLIST: ", a+1, main_pos[a][2], main_pos[a][3], main_pos[a][0], main_pos[a][1]);
	  for ( c=1; c<=comp_list[a][0]; c++)
	    {
	      fprintf ( stdout, "%s ",  get_compound_name(comp_list[a][c], mode));
	    }
	  fprintf ( stdout, " [COMP_LIST]");
	}
    }
  free_int (main_pos, -1);
  free_int (comp_list, -1);
  return P;
}

ORP* or_loo ( Alignment *inA, Alignment *inS, char *mode, int *pos,int print)
{
  int a, b,  n;
  char **keep, **remove;
  Alignment *A, *S, *P, *P1, *SEQ;

  int pos_exists=0;
  ORP *PR;



  if (pos)pos_exists=1;
  set_sar (inA, inS, mode);
  PR=declare_or_prediction (inS->nseq, inA->nseq, inA->len_aln);
  sprintf (PR->mode, "loo ");
  P=copy_aln (inS, NULL);


  keep=declare_char (inA->nseq, MAXNAMES);
  remove=declare_char (inA->nseq, MAXNAMES);


  PR->A=inA;
  PR->P=P;
  PR->S=inS;

  for (a=0; a<inA->nseq; a++)
    {
      char ***motifs;

      A=copy_aln (inA, NULL);
      S=copy_aln (inS, NULL);

      for (n=0,b=0; b<A->nseq; b++)
	{
	  if (b!=a)sprintf (keep[n++ ], "%s", A->name [b]);
	}
      sprintf ( remove[0], "%s", A->name[a]);
      reorder_aln (A,keep, A->nseq-1);

      set_sar (A,S, mode);

      if (!pos_exists)
	{
	  pos=aln2predictive_positions (A, S, mode,print);

	}

      for (b=0; b<A->len_aln; b++)
	{
	  PR->pos[b]+=pos[b];
	}

      motifs=compounds2motifs (A, S, pos,0, mode, print);

      SEQ=copy_aln (inA, NULL);
      SEQ=reorder_aln (SEQ, remove, 1);
      SEQ->nseq=1;

      P1=aln2prediction (SEQ, motifs, pos);


      if (print==PRINT)
	{
	  fprintf ( stdout, "\n%s\nPred: %s\nReal: %s\n", P1->name[0], P1->seq_al[0], inS->seq_al[a]);
	}
      sprintf ( P->seq_al[a], "%s", P1->seq_al[0]);
      free_aln (P1);


      free_aln (SEQ);
      free_aln (S);
      free_aln (A);

      free_arrayN( (void *)motifs, 3);
      if (!pos_exists)vfree (pos);

    }
  free_char (keep, -1);
  free_char (remove, -1);


  PR=new_evaluate_prediction (PR, mode,print);


  PR->PR=or_self_predict(inA, inS, mode,NULL, print);

  if (print==PRINT)for ( a=0; a<inA->len_aln; a++)fprintf ( stdout, "\nP: cons %d %d [FINAL]", a+1, PR->pos[a]);



  return PR;
}



Alignment * or_jack(Alignment *inA, Alignment *inS, char *mode)
{
  int a,b;
  int n_cycles=100;
  int subnseq=10;
  int subsar=0;
  Alignment *A, *S;
  int *main_pos,*pos;
  char jrsar[10], jraln[10];

  strget_param (mode, "_JNIT_", "100", "%d", &n_cycles);
  strget_param (mode, "_JNSEQ_", "10", "%d", &subnseq);
  strget_param (mode, "_JNAR_", "0", "%d", &subsar);

  strget_param (mode, "_JRSAR_", "NO", "%s", jrsar);
  strget_param (mode, "_JRALN_", "NO", "%s", jraln);



  main_pos=(int*)vcalloc ( inA->len_aln, sizeof (int));
  for (a=0; a<n_cycles; a++)
    {

      A=copy_aln (inA, NULL);
      S=copy_aln (inS, NULL);

      S=aln2random_aln (S, jrsar);
      A=aln2random_aln (A, jraln);

      set_sar (A,S, mode);
      sar2jack (A, S,subnseq,subsar);



      pos=aln2predictive_positions (A, S,mode, PRINT);

      for (b=0; b<A->len_aln; b++)main_pos[b]+=pos[b];
      vfree (pos);

    }
  display_pos (A, S, main_pos, mode);


  return inA;
}

Alignment * display_pos (Alignment *A, Alignment *S, int *pos,char *mode)
{
  Alignment *B;
  int a, b;
  int **index;

  int intl;

  intl=intlen (A->len_aln);
  index=aln2pos_simple (A, A->nseq);
  B=copy_aln (A,NULL);
  B->len_aln=0;
  for ( a=0; a<A->len_aln; a++)
    fprintf ( stdout, "\nP: cons %*d %*d S: %4d [DISPLAY_FULL_POS]", intl,a+1,intl, a+2, pos[a]);
  fprintf ( stdout, "\n\n");
  for (a=0; a<A->len_aln; a++)
    {
      if (pos[a])
	{
	  for ( b=0; b<A->nseq; b++)
	    {
	      B->seq_al[b][B->len_aln]=A->seq_al[b][a];
	      if (index[b][a]>0)fprintf ( stdout, "\nP: %s %d %d S: %d [DISPLAY_POS]",A->name[b], index[b][a], index[b][a]+1, pos[a]);
	    }
	  B->len_aln++;
	  fprintf ( stdout, "\nP: cons %d %d S: %d [DISPLAY_POS]", a+1, a+2, pos[a]);
	}
    }
  fprintf ( stdout, "\n");
  for (a=0; a<B->nseq; a++)B->seq_al[a][B->len_aln]='\0';
  return B;
}
Alignment * or_aln2pos_aln (Alignment *A, Alignment *S, char *mode)
{
  Alignment *B;

  int *pos;
  char outaln[100], outtree[100];


  strget_param (mode, "_OUTALN_", "NO", "%s", outaln);
  strget_param (mode, "_OUTTREE_", "NO", "%s", outtree);

  set_sar (A, S, mode);
  pos=aln2predictive_positions (A, S,mode, PRINT);

  B=display_pos (A, S, pos, mode);


  if (!strm(outaln, "NO")) vfclose (output_aln (B, vfopen (outaln, "w")));
  if (!strm(outtree, "NO"))vfclose (print_tree (aln2tree(B), "newick", vfopen (outtree, "w")));

  return B;
}
Alignment * or_sim(Alignment *A, Alignment *S, char *mode)
{
  //Predict all the sequences that are not both in inS and inA
  int *pos;

  set_sar (A, S, mode);
  pos=aln2predictive_positions (A, S,mode, PRINT);
  fprintf ( stdout, "R: %.3f", pos2sim (A,S, pos));

  myexit (EXIT_SUCCESS);
  return A;
}
ORP* or_self_predict(Alignment *A, Alignment *S, char *mode,int *pos, int print)
{
  //Predict all the sequences that are not both in inS and inA
  Alignment *P;
  char ***motifs;


  int a;

  int pre_set_pos=0;
  ORP *PR;


  set_sar (A, S, mode);
  PR=declare_or_prediction (S->nseq, A->nseq, A->len_aln);
  sprintf (PR->mode, "self");
  PR->A=A;
  PR->S=S;

  if (!pos)
    {
      pos=aln2predictive_positions (A, S,mode,print);
      pre_set_pos=0;
    }
  else
    pre_set_pos=1;

  for (a=0; a< A->len_aln; a++)
    PR->pos[a]=pos[a];


  PR->motif=motifs=compounds2motifs (A, S, pos,0, mode, print);
  P=PR->P=aln2prediction (A, motifs, pos);

  if (!pre_set_pos)vfree (pos);

  PR=new_evaluate_prediction (PR, mode,print);
  return PR;
}


Alignment * or_predict(Alignment *inA, Alignment *inS, char *mode)
{
  //Predict all the sequences that are not both in inS and inA
  Alignment *P, *A, *S, *T;
  char ***motifs;
  int *pos;

  int a, b;






  A=copy_aln (inA, NULL);
  S=copy_aln (inS, NULL);
  set_sar (A, S, mode);

  pos=aln2predictive_positions (A, S,mode,PRINT);
  motifs=compounds2motifs (A, S, pos,0, mode, PRINT);
  T=get_prediction_target (inA, inS, mode);


  P=aln2prediction (T, motifs, pos);
  //recall=evaluate_prediction (S, P, mode);
  for ( a=0; a<P->len_aln; a++)
    {
      for (b=0; b<P->nseq; b++)
	{
	  if (tolower(P->seq_al[b][a])=='i')fprintf (stdout, "\n>%20s %20s %c", T->name [0],get_compound_name (a, mode), P->seq_al[b][a]);
	}
    }
  fprintf ( stdout, "\n");
  return P;
}

Alignment *get_prediction_target (Alignment *A, Alignment *S, char *param)
{
  char **name;
  int n, a;
  Alignment *T;

  T=copy_aln (A, NULL);
  name=declare_char (A->nseq, 100);
  for (n=0,a=0; a< A->nseq; a++)
    {
      if ( name_is_in_hlist (A->name[a], S->name, S->nseq)==-1)
	{
	  sprintf (name[n++], "%s", A->name[a]);
	}
    }
  T=reorder_aln (T,name, n);
  return T;
}

Alignment *set_sar (Alignment *A, Alignment *S, char *param)
{
  char **name;
  int n, a;

  name=declare_char (A->nseq, 100);
  for (n=0,a=0; a< A->nseq; a++)
    {
      if ( name_is_in_hlist (A->name[a], S->name, S->nseq)!=-1)
	{
	  sprintf (name[n++], "%s", A->name[a]);
	}
    }
  A=reorder_aln (A,name, n);
  S=reorder_aln (S,name, n);
  free_char (name, -1);
  return S;
}

ORP* new_evaluate_prediction  (ORP *PR, char *mode, int print)
{
  int a,b, i, r, p;
  int tp, tn, fp, fn;
  float sn, sp, sen2, best;
  float tot_best_seq=0;
  float tot_best_comp=0;
  Alignment *P, *R;

  int ns=0;
  float *recall;


  P=PR->P;
  R=PR->S;

  recall=(float*)vcalloc (P->len_aln, sizeof (float));
  if (P->len_aln!=R->len_aln)
    {
      HERE ("Mismatch between number of compounds in prediction and reference");
      myexit (EXIT_FAILURE);
    }
  if (print==PRINT)fprintf ( stdout, "\n");

  for (a=0; a<P->nseq; a++)
    {
      tp=tn=fp=fn=0;
      if ((i=name_is_in_hlist (P->name[a], R->name, R->nseq))!=-1)
	{

	  for (b=0;b<P->len_aln; b++)
	    {
	      r=R->seq_al[i][b];
	      p=P->seq_al[a][b];

	      if ( p=='I' && r=='I')tp++;
	      else if ( p=='I' && r=='O')fp++;
	      else if ( p=='O' && r=='I')fn++;
	      else if ( p=='O' && r=='O')tn++;
	    }
	  rates2sensitivity (tp, tn, fp, fn, &sp, &sn, &sen2, &best);
	  if (print==PRINT)fprintf (stdout, ">%-s sp: %.2f sn: %.2f sn2: %.2f best: %.2f [SEQ]\n",P->name[a], sp, sn, sen2, best);
	  if ( best>0)
	    {
	      ns++;
	      tot_best_seq+=best;
	    }
	}
    }
  if (ns)
    {
      tot_best_seq/=ns;
    }
  if (print==PRINT)fprintf ( stdout, ">TotSeq sp: %.2f N: %d[SEQ]\n",tot_best_seq, ns);

  tot_best_comp=0;
  for (ns=0,b=0; b<P->len_aln; b++)
    {
      tp=tn=fp=fn=0;
      for (a=0; a<P->nseq;a++)
	{
	  if ((i=name_is_in_hlist (P->name[a], R->name, R->nseq))!=-1)
	    {
	      r=R->seq_al[i][b];
	      p=P->seq_al[a][b];

	      if ( p=='I' && r=='I'){PR->tp++;tp++;}
	      else if ( p=='I' && r=='O'){PR->fp++;fp++;}
	      else if ( p=='O' && r=='I'){PR->fn++;fn++;}
	      else if ( p=='O' && r=='O'){PR->tn++;tn++;}
	    }
	}
      rates2sensitivity (tp, tn, fp, fn, &sp, &sn, &sen2, &best);

      if (print==PRINT) fprintf (stdout, ">%-25s sp: %.2f sn: %.2f sen2: %.2f best: %.2f [COMP]\n",get_compound_name (b, mode), PR->sp, PR->sn, PR->sen2,PR->best);
      if ( best>0)
	{
	  ns++;
	  tot_best_comp+=best;
	}
    }

  if (ns)
    {
      tot_best_comp/=ns;
    }
  rates2sensitivity (PR->tp, PR->tn, PR->fp,PR->fn,&(PR->sp), &(PR->sn), &(PR->sen2), &(PR->best));
  if (print==PRINT)fprintf ( stdout, ">FullTot sp: %.2f sn: %.2f sen2: %.2f best: %.2f  N: %d[COMP]\n", PR->sp, PR->sn, PR->sen2,PR->best, ns);
  P->score_aln=(int)((float)1000*(PR->best));
  return PR;
}
float  evaluate_prediction  (Alignment *R, Alignment *P, char *mode, int print)
{
  int a,b, i, r, p;
  int tp, tn, fp, fn;
  int tot_tp, tot_tn, tot_fp, tot_fn;
  float sn, sp, sen2, best;
  float tot_sp=0;
  float tot_sn=0;
  float tot_sen2=0;
  float tot_best_seq=0;
  float tot_best_comp=0;
  float tot_best=0;

  int ns=0;
  float *recall;




  recall=(float*)vcalloc (P->len_aln, sizeof (float));
  if (P->len_aln!=R->len_aln)
    {
      HERE ("Mismatch between number of compounds in prediction and reference");
      myexit (EXIT_FAILURE);
    }
  if (print==PRINT)fprintf ( stdout, "\n");
  for (a=0; a<P->nseq; a++)
    {
      tp=tn=fp=fn=0;
      if ((i=name_is_in_hlist (P->name[a], R->name, R->nseq))!=-1)
	{

	  for (b=0;b<P->len_aln; b++)
	    {
	      r=R->seq_al[i][b];
	      p=P->seq_al[a][b];

	      if ( p=='I' && r=='I')tp++;
	      else if ( p=='I' && r=='O')fp++;
	      else if ( p=='O' && r=='I')fn++;
	      else if ( p=='O' && r=='O')tn++;
	    }
	  rates2sensitivity (tp, tn, fp, fn, &sp, &sn, &sen2, &best);
	  if (print==PRINT)fprintf (stdout, ">%-s sp: %.2f sn: %.2f sn2: %.2f best: %.2f [SEQ]\n",P->name[a], sp, sn, sen2, best);
	  if ( best>0)
	    {
	      ns++;
	      tot_best_seq+=best;
	      tot_sn+=sn;
	      tot_sp+=sp;
	      tot_sen2+=sen2;
	    }
	}
    }
  if (ns)
    {
      tot_best_seq/=ns;
      tot_sn/=ns;
      tot_sp/=ns;
      tot_sen2/=ns;
    }
  if (print==PRINT)fprintf ( stdout, ">Tot sp: %.2f sn: %.2f sen2: %.2f best: %.2f  N: %d[SEQ]\n", tot_sp, tot_sn, tot_sen2,tot_best_seq, ns);

  tot_fp=tot_fn=tot_tp=tot_tn=0;
  tot_sp=tot_sn=tot_sen2=tot_best_comp=0;
  for (ns=0,b=0; b<P->len_aln; b++)
    {
      tp=tn=fp=fn=0;
      for (a=0; a<P->nseq;a++)
	{
	  if ((i=name_is_in_hlist (P->name[a], R->name, R->nseq))!=-1)
	    {
	      r=R->seq_al[i][b];
	      p=P->seq_al[a][b];

	      if ( p=='I' && r=='I'){tot_tp++;tp++;}
	      else if ( p=='I' && r=='O'){tot_fp++;fp++;}
	      else if ( p=='O' && r=='I'){tot_fn++;fn++;}
	      else if ( p=='O' && r=='O'){tot_tn++;tn++;}
	    }
	}
      rates2sensitivity (tp, tn, fp, fn, &sp, &sn, &sen2, &best);

      if (print==PRINT) fprintf (stdout, ">%-25s sp: %.2f sn: %.2f sen2: %.2f best: %.2f [COMP]\n",get_compound_name (b, mode), sp, sn, sen2,best);
      recall[b]=sen2;
      if ( best>0)
	{
	  ns++;
	  tot_best_comp+=best;
	  tot_sn+=sn;
	  tot_sp+=sp;
	  tot_sen2+=sen2;
	}
    }

  if (ns)
    {
      tot_best_comp/=ns;
      tot_sn/=ns;
      tot_sp/=ns;
      tot_sen2/=ns;
    }
  rates2sensitivity (tot_tp, tot_tn, tot_fp,tot_fn,&tot_sp, &tot_sn, &tot_sen2, &tot_best);
  if (print==PRINT)fprintf ( stdout, ">FullTot sp: %.2f sn: %.2f sen2: %.2f best: %.2f  N: %d[COMP]\n", tot_sp, tot_sn, tot_sen2,tot_best, ns);
  return tot_best;
}



Alignment * aln2prediction (Alignment *A,char ***motif, int *pos)
{
  int a, b,nc, nl;
  int *list;
  char **array, **sar;
  Alignment *R;
  Sequence *S;
  nc=read_array_size ((void *)motif, sizeof (char***));


  list=pos2list (pos, A->len_aln, &nl);


  array=declare_char (A->nseq, nl+1);
  sar=declare_char(A->nseq, nc+1);
  for (a=0; a<A->nseq; a++)
    {
      for (b=0; b<nl; b++)
	array[a][b]=A->seq_al[a][list[b]];
    }

  for (a=0; a<nc; a++)
    {
      for (b=0; b<A->nseq; b++)
	{

	  sar[b][a]=(match_motif (array[b], motif[a]))?'I':'O';
	}
    }


  S=fill_sequence_struc (A->nseq,sar,A->name, NULL);
  R=seq2aln (S, NULL, KEEP_GAP);
  free_sequence (S, S->nseq);
  free_char (sar, -1);
  vfree (list);
  free_char (array, -1);
  return R;
}

int *   file2pos_list (Alignment *A, char *posfile)
{
  char ***file;
  int **index;
  int *pos;
  int i, n, p;

  //pos_file:<P: seqname pos score\n>
  //          1  2       3   4


  if ( !check_file_exists (posfile))
    {
      printf_exit ( EXIT_FAILURE, stderr, "ERROR: Could not read posfile %s\n", posfile);
    }

  file=file2list (posfile, " ");

  index=aln2inv_pos (A);
  pos=(int*)vcalloc ( A->len_aln, sizeof (int));

  n=0;
  while (file[n])
    {

      if ( !strm (file[n][1], "P:"));
      else
	{
	  if ( (strm (file[n][2], "cons")))
	    p=atoi(file[n][3])-1;
	  else
	    {
	      i=name_is_in_hlist ( file[n][2], A->name, A->nseq);
	      if (i!=-1)
		p=index[i][atoi(file[n][3])]-1;
	      else p=-1;
	    }
	  if (p!=-1)pos[p]+=atoi(file[n][4]);
	}
      n++;
    }


  free_int (index, -1);
  free_arrayN ( (char **)file, 3);
  return pos;
}
int *   aln2predictive_positions (Alignment *A, Alignment *B, char *mode, int print)
{
  char posmode[100];

  if (!mode) return NULL;
  HERE ("%s", mode);
  strget_param (mode, "_POSMODE_", "scan", "%s", posmode);
  if ( strm (posmode, "mat"))return aln2predictive_positions_mat (A, B, mode, print);
  else if ( strm (posmode, "scan")) return aln2predictive_positions_scan (A, B, mode, print);
  else
    {
      printf_exit (EXIT_FAILURE,stderr, "ERROR: %s is an unknown _POSMODE_ mode",posmode);
      return NULL;
    }
}

int *   aln2predictive_positions_mat (Alignment *A, Alignment *B, char *mode, int print)
  {
    int a, b, c,gap,  res1, res2, sar1, sar2, npos, s, idscore;
    float id1,id2,id3,nid1,nid2,nid3;
    int **pos, *fpos;
    pos=declare_int (A->len_aln,2);
    fpos=(int*)vcalloc ( A->len_aln, sizeof (int));

    strget_param (mode, "_NPOS_", "2", "%d", &npos);
    for ( a=0; a< A->len_aln; a++)
      {
	pos[a][0]=a;
	id1=id2=id3=nid1=nid2=nid3=0;
	for ( gap=0,b=0; b<A->nseq; b++)gap+=(A->seq_al[b][a]=='-');
	if ( gap>0){pos[a][1]=0;continue;}

	for (s=0; s<B->len_aln; s++)
	  {
	    for ( gap=0,b=0; b<A->nseq-1; b++)
	      {
		sar1=B->seq_al[b][s];
		res1=A->seq_al[b][a];

		for ( c=b+1; c<A->nseq; c++)
		  {
		    sar2=B->seq_al[c][s];
		    res2=A->seq_al[c][a];

		    idscore=(res1==res2)?1:0;
		    if ( sar1 == 'I' && sar2=='I'){id1+=idscore;nid1++;}
		    else if ( sar1 =='0' && sar2=='0'){id2+=idscore;nid2++;}
		    else {id3+=idscore; nid3++;}

		  }
	      }
	    id1=(nid1==0)?1:id1/nid1;
	    id2=(nid1==0)?1:id2/nid2;
	    id3=(nid3==0)?1:id3/nid3;
	    pos[a][1]=(int)((float)1000*id1*(1-id3));

	  }
      }

    sort_int (pos, 2,1, 0, A->len_aln-1);
    for ( a=MAX(0,(A->len_aln-npos));a<A->len_aln; a++)
      {
	fpos[pos[a][0]]=1;
      }

    free_int (pos, -1);
    return fpos;
  }
int *   aln2predictive_positions_scan (Alignment *A, Alignment *B, char *mode, int print)
{
  int a, b, c, best_pos,nl, nplist=0, max, posw;
  float best_score, score;
  static int *list, *tpos,**plist,*array;
  int *pos;


  char posfile[100];
  char predmode[100];
  char target_posfile[100];



  if (!A)
    {
      vfree (list);
      vfree (tpos);

      free_int (plist, -1);
      vfree (array);
      return NULL;
    }

  strget_param (mode, "_PREDMODE_", "ID", "%s", predmode);
  strget_param (mode, "_POSW_", "1", "%d", &posw);
  strget_param (mode, "_NPOS_", "2", "%d", &max);
  strget_param (mode, "_POSFILE_", "NO", "%s", posfile);
  strget_param (mode, "_TPOSFILE_", "NO", "%s", target_posfile);

  if ( !strm(posfile, "NO"))return file2pos_list (A,posfile);
  if ( !strm(target_posfile, "NO"))tpos=file2pos_list (A,target_posfile);
  else
    {
      tpos=(int*)vcalloc (A->len_aln, sizeof (int));
      for (a=0; a<A->len_aln; a++)tpos[a]=1;
    }

  //Declare the positions that are going to be scanned


  if (posw==1)
    {
      plist=declare_int (A->len_aln, 2);
      nplist=0;
      for (a=0; a<A->len_aln; a++)
	{
	  if(tpos[a])
	    {
	      plist[nplist][0]=1;
	      plist[nplist][1]=a;
	      nplist++;
	    }
	}
    }
  else if ( posw==2)
    {
      nplist=0;
      plist=declare_int (A->len_aln*A->len_aln, 3);
      for (a=0; a<A->len_aln; a++)
	for (b=0; b<A->len_aln; b++)
	  {
	    plist[nplist][1]=a;
	    plist[nplist][2]=b;
	    plist[nplist][0]=2;
	    nplist++;
	  }
    }
  else if ( posw==3)
    {
      nplist=0;
      plist=declare_int (A->len_aln*A->len_aln*A->len_aln, 3);
      for (a=0; a<A->len_aln; a++)
	for (b=0; b<A->len_aln; b++)
	  {
	    plist[nplist][1]=a;
	    plist[nplist][2]=b;
	    plist[nplist][3]=0;


	    plist[nplist][0]=3;
	    nplist++;
	  }
    }


  pos=(int*)vcalloc ( A->len_aln, sizeof (int));
  if (max==0)max=A->len_aln;
  else if ( max==-1)
    {
      for (a=0; a<A->len_aln; a++)if (tpos[a]){pos[a]=1;}
      aln2predictive_positions_scan (NULL, NULL, NULL, 0);
      return pos;
    }



  pos=(int*)vcalloc ( A->len_aln, sizeof (int));
  list=(int*)vcalloc (A->len_aln, sizeof (int));
  nl=0;



  for (a=0; a< max; a++)
    {
      int previous_best_pos=-1;
      for (best_score=-9999,best_pos=0,b=0; b<nplist; b++)
	{
	  for (c=0; c<nl; c++)pos[list[c]]=1;
	  for (c=1; c<=plist[b][0]; c++)pos[plist[b][c]]=1;

	  if (strm(predmode, "R"))score=sar_aln2r(A,B,pos,0);
	  else if ( strm (predmode, "ID"))
	    {
	      score  =or_id_evaluate (A, B, mode, pos,NO_PRINT);
	    }
	  else if ( strm (predmode, "BP2"))score =or_loo_evaluate2 (A, B, mode, pos,NO_PRINT);
	  else
	    {
	      HERE ("Unknown mode: %s", predmode);
	      myexit (EXIT_FAILURE);
	    }
	  if ( score>best_score)
	    {
	      best_score=score;
	      best_pos=b;
	    }
	  for (c=1; c<=plist[b][0]; c++)pos[plist[b][c]]=0;

	}
      if (best_pos==previous_best_pos)break;
      else previous_best_pos=best_pos;

      //update the best_pos_list
      for (b=1; b<=plist[best_pos][0]; b++)
	list[nl++]=plist[best_pos][b];


      if ( print==PRINT)
	{
	  for (b=0; b<nl; b++) pos[list[b]]=1;
	  fprintf ( stdout, "\nP_current: ");
	  for ( c=1; c<=plist[best_pos][0]; c++)fprintf ( stdout, "%d ",plist[best_pos][c]+1);
	  fprintf ( stdout, " S: %.3f D: %d R:%d", best_score, (int)sar_aln2delta(A,B, pos,0), nl);
	}
      for (b=0; b<nl; b++) pos[list[b]]=0;
    }

  for (a=0; a<nl; a++)
    pos[list[a]]=1;
  if (print==PRINT)fprintf ( stdout, "\nR_best: %.3f with %d pos" ,best_score, nl);

  aln2predictive_positions_scan (NULL, NULL, NULL, 0);

  return pos;
}

char *** compounds2motifs (Alignment *A, Alignment *B, int *pos, int depth, char *mode, int print)
{
  char ***motifs;
  int a;

  motifs=(char***)vcalloc (B->len_aln, sizeof (char**));
  for (a=0; a<B->len_aln; a++)
    {

      motifs[a]=compound2motif (A, B, pos, depth, a, mode, print);
    }

  return motifs;
}
char ** compound2regexp_motif (Alignment *A, Alignment *B, int *pos, int depth, int c, char *mode, int print);
char ** compound2word_motif (Alignment *A, Alignment *B, int *pos, int depth, int c, char *mode, int print);

char ** compound2motif (Alignment *A, Alignment *B, int *pos, int depth, int c, char *mode, int print)
{
  char mmode[100];

  strget_param (mode, "_MOTIFMODE_", "word", "%s", mmode); //words, regexp
  if ( strm (mmode, "regexp"))return compound2regexp_motif (A,B,pos, depth, c, mode, print);
  else if ( strm (mmode, "word"))return compound2word_motif (A,B,pos, depth, c, mode, print);
 else return NULL;}
char ** compound2word_motif (Alignment *A, Alignment *B, int *pos, int depth, int c, char *mode, int print)
{
  int a,l;
  char *word, **motif;
  float score;


  word=or_id_evaluate2 (A, B, mode, pos,print, &score);
  if ( !word) return NULL;
  l=strlen (word);

  motif=declare_char (l+1, 2);
  for (a=0; a<l; a++)motif[a][0]=word[a];

  if (print==PRINT)
    {
      fprintf ( stdout, "\nMotifCompound %25s best: %.2f  motif: ", get_compound_name(c, mode),score);
      for (a=0; a<l; a++)
	{
	  fprintf ( stdout, "[%2s]",motif[a]);
	}
    }
  vfree (word);
  return motif;
}



char ** compound2regexp_motif (Alignment *A, Alignment *B, int *pos, int depth, int c, char *mode, int print)
{
  static Alignment *I;
  static Alignment *O;
  int a, b, o, i;

  float tp,tn,fp,fn,best, sp, sn, sen2, best_sn, best_sp, best_sen2;
  float best_pred=-1;
  int best_motif=0;



  char ***alp=NULL;
  int *alp_size=NULL;

  char *motif_file;
  int n;
  char **m, **m2;

  char *buf=NULL;
  char *best_buf;
  FILE *fpp;

  free_aln (I);
  free_aln (O);

  I=copy_aln(A, NULL);
  O=copy_aln(A, NULL);

  if (depth==0)
    strget_param (mode, "_DEPTH_", "2", "%d", &depth);

  I->nseq=O->nseq=I->len_aln=O->len_aln=0;
  for (a=0; a<A->len_aln; a++)
    {
      if (pos[a])
	{
	  for (i=o=0,b=0; b<A->nseq; b++)
	    {
	      if ( is_gap(A->seq_al[b][a]))return 0;
	      if (B->seq_al[b][c]=='I')I->seq_al[i++][I->len_aln]=A->seq_al[b][a];
	      else O->seq_al[o++][O->len_aln]=A->seq_al[b][a];
	    }
	  I->len_aln++;
	  O->len_aln++;
	}
    }

  if (O->len_aln==0 || I->len_aln==0) return 0;
  O->nseq=o;
  I->nseq=i;
  for (a=0; a<o; a++)O->seq_al[a][O->len_aln]='\0';
  for (a=0; a<i; a++)I->seq_al[a][I->len_aln]='\0';

  if (!I->nseq) return NULL;



  best_pred=best_motif=best_sn=best_sp=best_sen2=0;

  motif_file=vtmpnam (NULL);

  n=0;
  if (depth>0)
    {
       alp=(char***)vcalloc ( sizeof (char**), I->len_aln);
       alp_size= (int*)vcalloc ( I->len_aln, sizeof (int));
       for (a=0; a<I->len_aln; a++)
	 {
	   char *col;
	   alp[a]=string2alphabet ( (col=aln_column2string (I,a)),depth, &alp_size[a]);
	   vfree (col);
	 }
       generate_array_string_list (I->len_aln, alp, alp_size, &n, motif_file, OVERLAP);
    }
  else
    {
      int *used;
      char r;

      used=(int*)vcalloc (256, sizeof (int));
      fpp=vfopen (motif_file,"w");
      for (a=0;a<I->len_aln; a++)
	{
	  for (b=0; b<I->nseq; b++)
	    {
	      r=I->seq_al[b][a];
	      if (!used[(int)r]){fprintf (fpp, "%c", r);used[(int)r]=1;}
	    }
	  for (b=0; b<I->nseq; b++)
	    {
	      r=I->seq_al[b][a];
	      used[(int)r]=0;
	    }
	  fprintf (fpp, " ");
	}
      fprintf (fpp, "\n");
      vfree (used);
      vfclose (fpp);

      n=1;
      depth=I->nseq;
    }

  buf=(char*)vcalloc (2*(I->len_aln*depth)+1, sizeof (char));
  best_buf=(char*)vcalloc (2*(I->len_aln*depth)+1, sizeof (char));
  fpp=vfopen (motif_file, "r");

  for (a=0; a<n; a++)
    {

      buf=vfgets (buf, fpp);
      m2=string2list (buf);
      m2++;


      tp=tn=fp=fn=0;

      for (b=0; b<I->nseq; b++)
	{
	  if (match_motif (I->seq_al[b], m2))tp++;
	  else fn++;
	}
      for (b=0; b<O->nseq; b++)
	{
	  if (match_motif (O->seq_al[b], m2))fp++;
	  else tn++;
	}
      rates2sensitivity (tp, tn, fp, fn, &sp, &sn, &sen2, &best);

      if (best>= best_pred)
	{
	  best_pred=best;
	  best_sp=sp;
	  best_sen2=sen2;
	  best_sn=sn;
	  sprintf (best_buf, "%s", buf);
	}
      m2--;
      free_char (m2, -1);
    }
  vfclose (fpp);
  if (print==PRINT)fprintf ( stdout, "\nMotifCompound %25s sp: %.2f sn: %.2f sen2: %.2f best: %.2f  motif: ", get_compound_name(c, mode), best_sp, best_sn, best_sen2, best_pred);
  m2=string2list (best_buf);
  m=declare_char (I->len_aln+1, depth+1);

  for (a=0; a<I->len_aln; a++)
	{
	  sprintf (m[a], "%s", m2[a+1]);
	  if (print==PRINT) fprintf ( stdout, "[%2s]",m[a]);
	}
  if (print==PRINT)fprintf ( stdout, " N-motifs %d", n);
  free_char (m2, -1);

  if (alp)free_arrayN((void ***) alp, 3);
  if (alp_size)vfree (alp_size);
  vfree (buf); vfree(best_buf);

  return m;
}

double pos2sim (Alignment *A, Alignment *B, int *pos)
{
  return sar_aln2r (A, B,pos, PRINT);
}
double  sar_aln2r (Alignment *A, Alignment *B, int *pos, int print)
{
  int a, b, c, d,r1, r2, n, score, sim;
  double *r, result;
  static double **slist;
  int declare=0;
  static int **M;



  if (!M)M=read_matrice ("blosum62mt");
  if (!slist)
    {
      int maxslist;
      maxslist=A->nseq*A->nseq*10;
      slist=declare_double (maxslist, 2);
    }

  if (pos==NULL)
    {

      declare=1;
      pos=(int*)vcalloc ( A->len_aln+1, sizeof (int));
      for (a=0; a<A->len_aln; a++)pos[a]=1;
      pos[a]=-1;

    }

  for (n=0,a=0; a< A->nseq-1; a++)
    {

      for (b=a+1; b<A->nseq; b++)
	{


	  for (sim=d=0,c=0; c<A->len_aln; c++)
	    {

	      if (pos[c]==0)continue;

	      r1=A->seq_al[a][c];
	      r2=A->seq_al[b][c];
	      if (is_gap(r1) || is_gap(r2))return 0;

	      sim+=M[r1-'A'][r2-'A']*pos[c];
	      d+=MAX((M[r1-'A'][r1-'A']),(M[r2-'A'][r2-'A']));
	    }
	  sim=(d==0)?0:(100*sim)/d;
	  score=(int)get_sar_sim(B->seq_al[a], B->seq_al[b]);
	  slist[n][0]=(double)sim;
	  slist[n][1]=(double)score;
	  if (print==PRINT)fprintf ( stdout, "SIM: %d %d [%s %s]\n", sim, score, A->name[a], A->name[b]);
	  n++;
	}
    }

  r=return_r(slist, n);
  for (a=0; a<n; a++)slist[a][0]=slist[a][1]=0;
  result=r[0];
  vfree (r);
  if (declare) vfree (pos);

  return result;
}


double sar_aln2delta (Alignment *A, Alignment *B, int *pos, int print)
{
  static Alignment *I;
  static Alignment *O;
  int a, b, c, o, i;
  double delta=0;
  if (!I)
    {
      I=copy_aln(A, NULL);
      O=copy_aln(A, NULL);
    }



  for (c=0; c<B->len_aln; c++)
    {

      I->nseq=O->nseq=I->len_aln=O->len_aln=0;
      for (a=0; a<A->len_aln; a++)
	{
	  if (pos[a])
	    {
	      for (i=o=0,b=0; b<B->nseq; b++)
		{
		  if ( is_gap(A->seq_al[b][a]))return 0;
		  if (B->seq_al[b][c]=='I')I->seq_al[i++][I->len_aln]=A->seq_al[b][a];
		  else O->seq_al[o++][O->len_aln]=A->seq_al[b][a];
		}
	      I->len_aln++;
	      O->len_aln++;
	    }
	}
      if (O->len_aln==0 || I->len_aln==0) return 0;
      O->nseq=o;
      I->nseq=i;
      for (a=0; a<o; a++)O->seq_al[a][O->len_aln]='\0';
      for (a=0; a<i; a++)I->seq_al[a][I->len_aln]='\0';

      delta+=aln2sim(I,"blosum62mt")-aln2sim(O, "blosum62mt");

    }

  return delta;
}

char * get_compound_name (int c, char *mode)
{
  static int isset;
  static Alignment *S;
  static char *lname;

  if (!isset)
    {
      char *comp_list;
      isset=1;
      lname=(char*)vcalloc (100, sizeof (char));

      if (!mode);
      else
	{
	  strget_param (mode, "_COMPLIST_", "NO", "%s", comp_list=(char*)vcalloc (100, sizeof (char)));
	  if (strm(comp_list, "NO"));
	  else
	    {
	      S=main_read_aln (comp_list, NULL);
	      vfree (comp_list);
	    }
	}
    }
  if (!S || c>=S->nseq)sprintf (lname, "%d", c);
  else
    {
      sprintf (lname, "%s", S->name [c]);
    }
  return lname;
}
ORP * declare_or_prediction ( int ncomp, int nseq, int len)
{
  ORP *P;
  P=(ORP*)vcalloc ( 1, sizeof (ORP));
  P->ncomp=ncomp;
  P->nseq=nseq;
  P->len=len;
  P->PR=NULL;

  P->pos=(int*)vcalloc (len+1, sizeof (int));

  return P;
}

void free_orp_list ( ORP**P)
{
  int a=0;
  while (P[a])
    {
      free_orp(P[a++]);
    }
}
void free_orp ( ORP*P)
{
  if (!P) return;
  free_aln (P->A);
  free_aln (P->S);
  free_aln (P->P);
  vfree (P->pos);
  free_arrayN((void **)P->motif, 3);
  if (P->PR)free_orp(P->PR);
  vfree (P);
}















