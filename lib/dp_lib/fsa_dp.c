#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"

#define hmm_add(x,y) ((x==UNDEFINED || y==UNDEFINED)?UNDEFINED:(x+y))
#define MAX_EMISSION 256

/*********************************************************************************/
/*                                                                               */
/*                                                                               */
/*                       Procons dp                                              */
/*                                                                               */
/*                                                                               */
/*********************************************************************************/
char alphabetDefault[] = "ARNDCQEGHILKMFPSTWYV";
double emitPairsDefault[20][20] = {
  {0.02373072f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00244502f, 0.01775118f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00210228f, 0.00207782f, 0.01281864f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00223549f, 0.00161657f, 0.00353540f, 0.01911178f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00145515f, 0.00044701f, 0.00042479f, 0.00036798f, 0.01013470f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00219102f, 0.00253532f, 0.00158223f, 0.00176784f, 0.00032102f, 0.00756604f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00332218f, 0.00268865f, 0.00224738f, 0.00496800f, 0.00037956f, 0.00345128f, 0.01676565f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00597898f, 0.00194865f, 0.00288882f, 0.00235249f, 0.00071206f, 0.00142432f, 0.00214860f, 0.04062876f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00114353f, 0.00132105f, 0.00141205f, 0.00097077f, 0.00026421f, 0.00113901f, 0.00131767f, 0.00103704f, 0.00867996f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00318853f, 0.00138145f, 0.00104273f, 0.00105355f, 0.00094040f, 0.00100883f, 0.00124207f, 0.00142520f, 0.00059716f, 0.01778263f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00449576f, 0.00246811f, 0.00160275f, 0.00161966f, 0.00138494f, 0.00180553f, 0.00222063f, 0.00212853f, 0.00111754f, 0.01071834f, 0.03583921f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00331693f, 0.00595650f, 0.00257310f, 0.00252518f, 0.00046951f, 0.00312308f, 0.00428420f, 0.00259311f, 0.00121376f, 0.00157852f, 0.00259626f, 0.01612228f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00148878f, 0.00076734f, 0.00063401f, 0.00047808f, 0.00037421f, 0.00075546f, 0.00076105f, 0.00066504f, 0.00042237f, 0.00224097f, 0.00461939f, 0.00096120f, 0.00409522f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00165004f, 0.00090768f, 0.00084658f, 0.00069041f, 0.00052274f, 0.00059248f, 0.00078814f, 0.00115204f, 0.00072545f, 0.00279948f, 0.00533369f, 0.00087222f, 0.00116111f, 0.01661038f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00230618f, 0.00106268f, 0.00100282f, 0.00125381f, 0.00034766f, 0.00090111f, 0.00151550f, 0.00155601f, 0.00049078f, 0.00103767f, 0.00157310f, 0.00154836f, 0.00046718f, 0.00060701f, 0.01846071f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00631752f, 0.00224540f, 0.00301397f, 0.00285226f, 0.00094867f, 0.00191155f, 0.00293898f, 0.00381962f, 0.00116422f, 0.00173565f, 0.00250962f, 0.00312633f, 0.00087787f, 0.00119036f, 0.00180037f, 0.01346609f, 0.0f, 0.0f, 0.0f, 0.0f},
  {0.00389995f, 0.00186053f, 0.00220144f, 0.00180488f, 0.00073798f, 0.00154526f, 0.00216760f, 0.00214841f, 0.00077747f, 0.00248968f, 0.00302273f, 0.00250862f, 0.00093371f, 0.00107595f, 0.00147982f, 0.00487295f, 0.01299436f, 0.0f, 0.0f, 0.0f},
  {0.00039119f, 0.00029139f, 0.00021006f, 0.00016015f, 0.00010666f, 0.00020592f, 0.00023815f, 0.00038786f, 0.00019097f, 0.00039549f, 0.00076736f, 0.00028448f, 0.00016253f, 0.00085751f, 0.00015674f, 0.00026525f, 0.00024961f, 0.00563625f, 0.0f, 0.0f},
  {0.00131840f, 0.00099430f, 0.00074960f, 0.00066005f, 0.00036626f, 0.00070192f, 0.00092548f, 0.00089301f, 0.00131038f, 0.00127857f, 0.00219713f, 0.00100817f, 0.00054105f, 0.00368739f, 0.00047608f, 0.00102648f, 0.00094759f, 0.00069226f, 0.00999315f, 0.0f},
  {0.00533241f, 0.00169359f, 0.00136609f, 0.00127915f, 0.00119152f, 0.00132844f, 0.00178697f, 0.00194579f, 0.00071553f, 0.01117956f, 0.00914460f, 0.00210897f, 0.00197461f, 0.00256159f, 0.00135781f, 0.00241601f, 0.00343452f, 0.00038538f, 0.00148001f, 0.02075171f}
};

static void DisplayMatState ( MatState *S, char *s);






void check_viterbiL ( Alignment *A,int *ns, int **ls, Constraint_list *CL);
MatState *viterbi2path2 ( double ***Sc, int ***St, Hmm *H, MatState *S, MatState *E);
void testfunc ( MatState *S, char *s);

#ifdef IN_PGROGRESS
/*********************************************************************************/
/*                                                                               */
/*                                                                               */
/*                     MSA Analyzer                                              */
/*                                                                               */
/*                                                                               */
/*********************************************************************************/
Alignment * analyze_alignment ( Alignment *A)
{
  evaluate_alignment (A);
  H=define_msa_model (-100);
  M=seq_viterbi_hmm (A->seq_al[0], H);
  path=seq_viterbi2path ( seq, H, M);
}
    

Hmm* define_msa_model(double penalty)
{
  Hmm *H;
  double freeT=0;
  int n=0;
  HmmState *S;
  
  
  H=declare_hmm(2);
  H->freeT=freeT=0;
 
  H->forbiden=FORBIDEN;
  H->start=START_STATE;
  H->end=END_STATE;
 
  /*define START*/
  S=H->S[n];  
  sprintf (S->name, "START"); S->state=n;
  
  S->DI=0;
  S->DJ=0;
  S->em=freeT;
 
  sprintf ( (S->T[S->nT])->name, "C") ;(S->T[S->nT])->tr=freeT  ;S->nT++;
  sprintf ( (S->T[S->nT])->name, "W");(S->T[S->nT])->tr=freeT  ;S->nT++;
  n++;
  /*define END*/
  S=H->S[n];  
  sprintf (S->name, "END"); S->state=n;
  S->DI=0;
  S->DJ=0;
  S->em=freeT;
  n++;

  /*define Correct*/
  S=H->S[n];  
  sprintf (S->name, "C"); S->state=n;
  S->DI=1;
  S->DJ=0;
  S->em=H->forbiden;
  S->em_func=em_correct_msa;
  
  sprintf ( (S->T[S->nT])->name, "C") ;(S->T[S->nT])->tr=freeT;S->nT++;
  sprintf ( (S->T[S->nT])->name, "W");(S->T[S->nT])->tr=penalty  ;S->nT++;
  sprintf ( (S->T[S->nT])->name, "END");   (S->T[S->nT])->tr=freeT;S->nT++;
  n++;
  
  /*define Wrong*/
  S=H->S[n];  
  sprintf (S->name, "INSERT"); S->state=n;
  S->DI=1;
  S->DJ=0;
  S->em=H->forbiden;
  S->em_func=em_wrong_msa;
  sprintf ( (S->T[S->nT])->name, "C") ; (S->T[S->nT])->tr=penalty;S->nT++;
  sprintf ( (S->T[S->nT])->name, "W"); (S->T[S->nT])->tr=freeT;S->nT++;
  sprintf ( (S->T[S->nT])->name, "END");    (S->T[S->nT])->tr=-gop;S->nT++;
  n++;
  
  /*define LInsert*/
  S=H->S[n];  
  sprintf (S->name, "LINSERT"); S->state=n;
  S->DI=1;
  S->DJ=0;
  S->em=lgep;
  
  sprintf ( (S->T[S->nT])->name, "INSERT")  ;(S->T[S->nT])->tr=freeT;S->nT++;
  sprintf ( (S->T[S->nT])->name, "LINSERT");(S->T[S->nT])->tr=freeT;S->nT++;
  n++;
  
  H=bound_hmm ( H);
  return H;
}
#endif
/*********************************************************************************/
/*                                                                               */
/*                                                                               */
/*                     simple HMM: Viterbi                                       */
/*                                                                               */
/*                                                                               */
/*********************************************************************************/
double pavie_em_func (Hmm*H, HmmState *S, int v);
Hmm* define_full_model(int nstate,char **state_list, char *model_name,Generic_em_func evaluation_func );
char **produce_state_name (int nstate,char **list, char *model_name, Hmm* H);
double** seq_viterbi_hmm (char *seq, Hmm *H);
int * seq_viterbi2path (char *s, Hmm *H, double **M);
double analyze_sequence ( char *seq, Hmm*H);

double pavie_emission (Hmm*H, HmmState *S, int v)
{
  char *n;


  n=S->name;
  
  if ( v==n[0] || ( v=='*' && n[0]=='E')) return H->freeT;
  return H->forbiden;
}				     
Hmm* define_full_model(int nstate, char **list, char *model_name, Generic_em_func emission_function)
{
  /*list: a list of the state names: does not include START or END*/
  /*model_name: a string that will be appended to the names in list*/

  Hmm *H;
  int a,n;
  HmmState *S;
  
  
  H=declare_hmm(nstate+2);
  H->freeT=0;
 
  H->forbiden=FORBIDEN;
  H->start=START_STATE;
  H->end=END_STATE;
 
  list=produce_state_name (nstate,list, model_name, H);
  nstate+=2;
  
  for (n=0; n<nstate; n++)
    {
      
      S=H->S[n];
      S->state=n;
      sprintf ( S->name, "%s", list[n]);
      S->em_func2=emission_function;
      if (n==H->end || n==H->start){S->DI=0;S->DJ=0;}
      else S->DI=1;S->DJ=0;
      
      /*Emmissions*/
      S->em_func2=emission_function;
      for (a=0; a< MAX_EMISSION; a++)S->em2[a]=H->freeT;
      
      for (a=0; a<nstate && n!=H->end; a++)
	{
	  if (a!=H->start && !(n==H->start && a==H->end) )
	  {
	    sprintf ( (S->T[S->nT])->name, "%s", list[a]);
	    (S->T[S->nT])->tr=H->freeT;
	    S->nT++;
	  }
	}
    }
  return H;
}

char **produce_state_name (int nstate,char **list, char *model_name, Hmm* H)
{
  int a,b,c;
  char **new_list;
  nstate+=2;
  
  new_list=declare_char ( nstate, 100);
  for ( a=0, b=0, c=0; a< nstate; a++)
    {
      if ( a==H->start)sprintf ( new_list[a], "START");
      else if ( a==H->end)sprintf ( new_list[a], "END");
      else if ( list==NULL){sprintf ( new_list[a], "%c%s", 'a'+b, (model_name)?model_name:"");b++;}
      else {sprintf ( new_list[a], "%c%s", list[a][c], (model_name)?model_name:"");c++;}
    }
  return new_list;
}

int seq_viterbi_pair_wise (Alignment *A,int*ns, int **ls,Constraint_list *CL)
{
  ungap(A->seq_al[0]);
  analyze_sequence (A->seq_al[0], NULL);
  myexit (EXIT_FAILURE);
  return 1;
}
double analyze_sequence ( char *seq, Hmm *H)
{

  double **M;
  int *path;

  if ( H==NULL)
    {
      H=define_full_model(5, NULL,"_first", pavie_emission);
      H=bound_hmm(H);
      DisplayHmm (H);
    }
  M=seq_viterbi_hmm (seq, H);
  path=seq_viterbi2path (seq, H, M);
  return M[H->end][strlen (seq)];
}
   

double** seq_viterbi_hmm (char *seq, Hmm *H)
{
  /*Given a model H and a sequence seq*/
  double **M;
  double e, v, max;
  int i,pi, bestk, s, k, l1;
  HmmState *S1, *S2;
  
  
  l1=strlen (seq);
  M=declare_double (H->nS*2,l1+2);
  
  /*Handle the start*/
  M[H->start][0]=0;
  for ( i=0; i<=l1; i++)
    {
      for ( s=0; s< H->nS; s++)
	{
	  S1=H->S[s];
	  pi=i-S1->DI;
	  max=H->forbiden;
	  bestk=H->forbiden;
	  if ( pi<0){M[s][i]=H->forbiden;}/*Boundary*/
	  else 
	    {
	      if (pi==0) {max=H->T[(int)H->start][s];bestk=H->start;}/*Start*/
	      else
		{
		  for (k=1; k<=H->fromM[S1->state][0]; k++)
		    {
		      S2=H->S[H->fromM[s][k]];
		      if ( S2->state==H->start || S2->state==H->end)continue;
		      v=hmm_add((M[S2->state][pi]),(H->T[S2->state][S1->state]));
		      if ( v!=H->forbiden && (max==H->forbiden || v>max)){max=v;bestk=S2->state;}
		    }
		} 	
	      if (S1->em2)e=S1->em2[(int)seq[pi]];
	      else e=S1->em_func2(H,S1, (int)seq[pi]);
	      
	      e=hmm_add (e,max);
	      
	      M[s][i]=e;
	      M[s+H->nS][i]=bestk;
	      }
	  }
      }   
  /*Terminate viterbi: connect the path to the END state*/
  max=UNDEFINED;
  bestk=UNDEFINED;
  for (k=0; k< H->nS; k++)
    {
      if (k==H->start || k==H->end);
      else
	{
	  v=(M[k][l1]==H->forbiden || H->T[k][H->end]==H->forbiden)?H->forbiden:M[k][l1]+H->T[k][H->end];
	  if ( max==H->forbiden || v>max){bestk=k;max=v;}
	}
    }
  M[H->end][l1]=max;
  M[H->nS+H->end][l1]=bestk; 
  return M;
}
  
int * seq_viterbi2path (char *s, Hmm *H, double **M)
{
  int i,l,l1;
  int *path;
  HmmState *S1;
  int cs;

  l1=strlen (s);
  path=(int*)vcalloc (l1+1, sizeof (int));
  i=l1;
  l=0;
  cs=M[H->nS+H->end][i];
  
  while (i>0)
    {
      
      S1=H->S[cs];
      path[l++]=cs;

      cs=M[H->nS+cs][i];
      i-=S1->DI;
      /*fprintf ( stderr, "%d", cs);*/
    }
  invert_list_int (path, l);
  path[l++]=H->forbiden;
  
  return path;
}

/*********************************************************************************/
/*                                                                               */
/*                                                                               */
/*                     pairHMM: Viterbi                                              */
/*                                                                               */
/*                                                                               */
/*********************************************************************************/
Hmm* define_mnm_model(Constraint_list *CL);
int viterbi_pair_wise_OLD (Alignment *A,int*ns, int **ls,Constraint_list *CL)
{
  int l1, l2, a;
  double ***M;
  int *path;
  Hmm * H;
  
  A->pos=aln2pos_simple( A, -1, ns, ls);
  
  //  H=define_mnm_model (CL);
  H=define_two_mat_model (CL);
  
  l1=strlen (A->seq_al[ls[0][0]]);
  l2=strlen (A->seq_al[ls[1][0]]);
  M=viterbi_hmm (A, ns, ls, H, CL);
  path=viterbi2path (l1,l2, H,M);
  A=viterbipath2aln (A,ns,ls,path, H); 
  A->score=A->score_aln=M[H->end][l1][l2];
  for ( a=0; a< H->nS*2; a++)free_double (M[a], -1);
  vfree (M);
  free_int (A->pos, -1);
  A->pos=NULL;

  free_Hmm (H);
  vfree (path);
  
  return A->score_aln;
}

Alignment * viterbipath2aln (Alignment *A, int *ns,int **ls,int *tb, Hmm *H)
{
  char **aln;
  char *char_buf;
  int a, b, c, len, ch;
  HmmState *S;
  int l[2];
  
  len=0;while (tb[len]!=H->forbiden)len++;
  
  if ( A->declared_len<=len)A=realloc_aln2  ( A,A->max_n_seq,2*len);
  aln=A->seq_al;
  
  char_buf=(char*)vcalloc (len+1, sizeof (char));
  l[0]=strlen ( A->seq_al[ls[0][0]]);
  l[1]=strlen ( A->seq_al[ls[1][0]]);
  
  for ( c=0; c< 2; c++)
    for ( a=0; a< ns[c]; a++) 
      {
	for (ch=0, b=0; b<len; b++)
	  {
	    S=H->S[tb[b]];
	    if ( (c==0 && S->DI)|| (c==1 && S->DJ) )
	      char_buf[b]=aln[ls[c][a]][ch++];
	    else 
	      char_buf[b]='-';
	  }
	char_buf[b]='\0';
	sprintf (aln[ls[c][a]],"%s", char_buf);
	if ( l[c]!=ch){fprintf (stderr, "\nERROR: Wrong Size Of Alignmnent (Real %d, Observed %d)[FATAL:%s]",l[c], ch, PROGRAM);}
      }
  A->len_aln=len;
  A->nseq=ns[0]+ns[1];
  
  vfree(char_buf);
  return A;
}

double*** viterbi_hmm (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL)
{
  double ***M;
  double e, v, max;
  int a, i,pi, bestk,j,pj, s, k, l1, l2;
  HmmState *S1, *S2;
  
  l1=strlen (A->seq_al[ls[0][0]]);
  l2=strlen (A->seq_al[ls[1][0]]);
  
  M=(double***)vcalloc (H->nS*2, sizeof (double**));
  for ( a=0; a<H->nS*2; a++)M[a]=declare_double (l1+2, l2+2); 
  
  /*Handle the start*/
  
  M[H->start][0][0]=0;
  for ( i=0; i<=l1; i++)
    for ( j=0; j<=l2; j++)
      {
	for ( s=0; s< H->nS; s++)
	  {
	    S1=H->S[s];
	    pi=i-S1->DI;
	    pj=j-S1->DJ;
	    max=H->forbiden;
	    bestk=H->forbiden;
	    if ( pi<0 ||pj<0){M[s][i][j]=H->forbiden;}/*Boundary*/
	    else 
	      {
		if (pi+pj==0) {max=H->T[H->start][s];bestk=H->start;}/*Start*/
		else
		  {
		    for (k=1; k<=H->fromM[S1->state][0]; k++)
		      {
			S2=H->S[H->fromM[s][k]];
			if ( S2->state==H->start || S2->state==H->end)continue;
			v=(M[S2->state][pi][pj]==H->forbiden)?H->forbiden:(M[S2->state][pi][pj]+H->T[S2->state][S1->state]);
			if ( v!=H->forbiden && (max==H->forbiden || v>max)){max=v;bestk=S2->state;}
		      }
		  }
	
		e=(S1->em==H->forbiden)?S1->em_func (A, A->pos, ns[0], ls[0],i-1, A->pos,ns[1], ls[1], j-1, CL):S1->em;
		e=(max==H->forbiden || e==H->forbiden)?H->forbiden:e+max;
		
		M[s][i][j]=e;
		M[s+H->nS][i][j]=bestk;
	      }
	  }
      }
  
  /*Terminate viterbi: connect the path to the END state*/
  max=UNDEFINED;
  bestk=UNDEFINED;
  for (k=0; k< H->nS; k++)
    {
      if (k==H->start || k==H->end);
      else
	{
	  v=(M[k][l1][l2]==H->forbiden || H->T[k][H->end]==H->forbiden)?H->forbiden:M[k][l1][l2]+H->T[k][H->end];
	  if ( max==H->forbiden || v>max){bestk=k;max=v;}
	}
    }
  M[H->end][l1][l2]=max;
  M[H->nS+H->end][l1][l2]=bestk; 
    
  return M;
}

/*********************************************************************************/
/*                                                                               */
/*                                                                               */
/*                      HMM: Decode/Traceback                                       */
/*                                                                               */
/*                                                                               */
/*********************************************************************************/
int * traceback (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL,MatState *S, MatState *E, int **seg_list)

{
  int *path;
  int l=0;
  MatState *N;
  int l1, l2;
  
  l1=strlen (A->seq_al[ls[0][0]]);
  l2=strlen (A->seq_al[ls[1][0]]);
  path=(int*)vcalloc ( l1+l2+1, sizeof (int));
  
  while ( S->st!=H->end)
    {
      DisplayMatState (S, "\n\tTraceback");
      N=S->n;
      if ( N && S && (((N->i-S->i)>1) ||((N->j-S->j)>1)))
	{
	  RviterbiD_hmm (A,ns,ls,H,CL,S,N,seg_list);
	  N=S->n;
	}

      path[l++]=S->st;
      ManageMatState (FREE,S);
      S=N;
    }
  
  path[l]=H->forbiden;
  return path;
}
	
int * viterbi2path (int l1,int l2, Hmm *H, double ***M)
{
  int i, j,l;
  int *path;
  HmmState *S1;
  int cs;

  l=0;
  path=(int*)vcalloc (l1+l2+1, sizeof (int));
  i=l1;j=l2;
  l=0;
  cs=M[H->nS+H->end][i][j];
  
  while (i>0|| j>0)
    {
      
      S1=H->S[cs];
      path[l++]=cs;

      cs=M[H->nS+cs][i][j];
      i-=S1->DI;
      j-=S1->DJ;
      /*fprintf ( stderr, "%d", cs);*/
    }
  invert_list_int (path, l);
  path[l++]=H->forbiden;
  
  return path;
}

/*********************************************************************************/
/*                                                                               */
/*                                                                               */
/*                      HMM Viterbi Linear                                       */
/*                                                                               */
/*                                                                               */
/*********************************************************************************/


int viterbiL_pair_wise (Alignment *A,int*ns, int **ls,Constraint_list *CL)
{
  int l1, l2;
  int *path;
  Hmm * H;
  MatState *Start;
  MatState *End;

  A->pos=aln2pos_simple( A, -1, ns, ls);
  Start=ManageMatState ( DECLARE, NULL);
  End=ManageMatState ( DECLARE, NULL);
  H=define_simple_model (CL);
  l1=strlen (A->seq_al[ls[0][0]]);
  l2=strlen (A->seq_al[ls[1][0]]);

  

  Start->i=0 ;Start->j=0 ; Start->st=H->start;Start->sc=0;
  End->i  =l1;  End->j=l2; End  ->st=H->end;
  Start=RviterbiL_hmm (A, ns, ls, H, CL, Start,End);
  path=traceback (A, ns, ls, H, CL, Start,NULL, NULL);
  
  A=viterbipath2aln (A,ns,ls,path, H); 
  
  free_Hmm (H);
  free_int (A->pos, -1);
  A->pos=NULL;
  
  return A->score_aln;
}


MatState* RviterbiL_hmm (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL,MatState *S, MatState *E)
{
  MatState *MS, *ME;
  MS=S;
  ME=E;
  
  viterbiL_hmm (A,ns, ls,H, CL, S, E);
  
  
  if ( S->n==E)return S;
  if ( E->sc==H->forbiden)
	{
	  DisplayHmm (H);
	  fprintf ( stderr, "\nERROR: The Requested Model (Cf Above) Cannot Produce the Pair-Alignment\nYou must allow extra possible transitions\n[FATAL:%s]", PROGRAM);
	  myexit ( EXIT_FAILURE);
	}
  E=S->n;
  
  while (S!=ME)
    {
      int d1, d2, align;
      d1=MinDeltaMatState(S,E);
      d2=MaxDeltaMatState(S,E);
      align=((d1==1 && d2==1) || ( d1==0))?0:1;
      if (align)RviterbiL_hmm (A,ns, ls,H, CL,S,E);
      S=E;
      E=S->n;
    }
  return MS;
}

#define Dim(i,j) (i*LenJ+j)
MatState* viterbiL_hmm (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL, MatState *S,MatState *E)
{
  int current, previous,row, prow;
  double v;
  int a,i,j,pi,pj, s, k;
  int start_i, start_j, end_i, end_j, l1, l2;
  HmmState *S1, *S2;
  MatState *CC, *PCC,*tS, *tE, *mark=NULL;
  int midpoint;


  static MatState ***M;


  static int LenJ, LenI;
  int MaxDelta=50, DeltaI, DeltaJ;

  

  DisplayMatState (S, "\n\tS");
  DisplayMatState (E, "\n\tE");
  
  
  if ( A==NULL)
    {
      for ( a=0; a<2; a++)memset(M[a],0,LenJ*LenI*sizeof (MatState*));
      free_arrayN((void **)M, 3);M=NULL;
      ManageMatState ( FREE_ALL, NULL);
      return NULL;
    }

  
  if ( MatStateAreIdentical ( S, E))return NULL;
  l1=strlen (A->seq_al[ls[0][0]]);l2=strlen (A->seq_al[ls[1][0]]);
  
  midpoint=S->i+((E->i-S->i)/2);
  DeltaI=E->i-S->i;
  DeltaJ=E->j-S->j;
  
  start_i=S->i;end_i=E->i;start_j=S->j;end_j=E->j;
  current=0;previous=1;
 
  
  if ( !M)
    {
      LenI=l2+1;
      LenJ=H->nS;
      M=(MatState***)declare_arrayN(3, sizeof ( MatState),2,LenI*LenJ,0);
    }
  
  
  /*MAKE THE VITERBI FROM S(tart) to E(nd)*/
  mark=ManageMatState ( MARK, mark);
  for (i=start_i; i<=end_i; i++)
    {
      row=current;
      for ( j=start_j; j<=end_j; j++)
	{
	  DeltaJ=((FABS(j-i))<MaxDelta)?1:MaxDelta+1; /*Keep Pointers on a DealtaMax Band around the main diagonal*/
	  for ( s=H->nS-1;s>=0; s--)
	    {
	      S1=H->S[s];pi=i-S1->DI;prow=S1->DI;pj=j-S1->DJ;
	      	      
	      CC=M[row][Dim(j,s)]=CopyMatState(NULL, M[row][Dim(j,s)]);
	      CC->i=i; CC->j=j; CC->st=s;PCC=NULL;
	      	      
	      if (i==start_i && j==start_j && s==S->st){CC=CopyMatState(S,CC);}
	      else if ( i==end_i && j==end_j && E->st!=H->end && s!=E->st)CC->sc=H->forbiden;
	      else if ( pi<start_i || pj<start_j)        {CC->sc=H->forbiden;}
	      else 
		{
		  for (k=1; k<=H->fromM[S1->state][0]; k++)
		    {
		      S2=H->S[H->fromM[s][k]];
		      PCC=M[prow][Dim((j-S1->DJ),(S2->state))];
		      		    
		      if ( !PCC)PCC=NULL;
		      else if      ( pi+pj!=0 && S2->state==H->start);
		      else if ( !(pi==l1 && pj==l2) && s==H->end);
		      else
			{
			  
			  v=hmm_add(CC->sc,H->T[PCC->st][CC->st]);
			  
			  v=lu_RviterbiD_hmm(A,ns, ls, H, CL,PCC,CC, NULL);
			  if ( v!=H->forbiden && (CC->sc==H->forbiden || v> CC->sc)){CC->sc=v; CC->pst=S2->state;CC->p=PCC;}
			}
		    }
		}
	      if (CC->sc==H->forbiden);
	      else if (i==midpoint || DeltaI<=MaxDelta||DeltaJ<=MaxDelta ||(i==start_i && j==start_j && s==S->st) )
		{
		  CC->m=(CC->p)?(CC->p)->m:NULL;
		  PCC=CopyMatState(CC,NULL);
		  PCC->m=CC->m;CC->m=PCC;
		}
	      else CC->m=(CC->p)?(CC->p)->m:NULL;
	    }
	}
      prow=previous;
      for ( j=start_j; j<=end_j && i!=end_i; j++)
	{
	  for ( s=H->nS-1;s>=0; s--)
	    {
	      
	      CC=(M[prow][Dim(j,s)]);M[prow][Dim(j,s)]=M[row][Dim(j,s)];M[row][Dim(j,s)]=CC;
	      if (M[prow][Dim(j,s)]) M[row ][Dim(j,s)]=CopyMatState ( M[prow][Dim(j,s)], M[row][Dim(j,s)]);
	      
	    }
	} 
      
    }  
  
  mark=ManageMatState ( MARK,mark);
  row=current;
  
  
  if ( E->st==H->end || E->st==H->forbiden){E=CopyMatState ((M[row][Dim(end_j,E->st)]),E);}

  
  

  PCC=CopyMatState (M[row][Dim(end_j,E->st)], NULL);
  
  if ( MatStateAreIdentical(PCC,PCC->m))PCC=PCC->m;
    tS=tE=PCC;
  while (PCC->m)
    {
      tS=CopyMatState (PCC->m,NULL); tS->n=PCC; PCC->p=tS;PCC=tS;
    }
    
  if (tS==tE);
  else
    {
      S->n=tS->n; (S->n)->p=S;
      E->p=tE->p; (E->p)->n=E;
    }
  for ( a=0; a<2; a++)memset(M[a],0,LenJ*LenI*sizeof (MatState*));
  ManageMatState ( FREE_MARK,mark);
  
    
  while (S && S->p!=E){S->m=NULL;S=S->n;}/*Clean the memory of the rturned Cells*/
  return NULL;
}

/*********************************************************************************/
/*                                                                               */
/*                                                                               */
/*                      HMM Viterbi Diagonals                                    */
/*                                                                               */
/*                                                                               */
/*********************************************************************************/
int viterbiD_pair_wise (Alignment *A,int*ns, int **ls,Constraint_list *CL)
{
  int l1, l2;
  int *path;
  Hmm * H;
  MatState *Start;
  MatState *End;
  int **seg_list;
  int a, b, c;
  int main_i;
  int main_j;

  
  A->pos=aln2pos_simple( A, -1, ns, ls);
  
  Start=ManageMatState ( DECLARE, NULL);
  End=ManageMatState ( DECLARE, NULL);
  H=define_simple_model (CL);
  l1=strlen (A->seq_al[ls[0][0]]);
  l2=strlen (A->seq_al[ls[1][0]]);

  main_i=MAX(1,(l2-l1)+1);
  main_j=MAX(1,(l1-l2)+1);
  
  seg_list=(int**)declare_arrayN(2, sizeof (int), l1+l2+3, 3);
  seg_list[0][0]=DIAGONALS;

  
  c=1;
  for ( b=1,a=l1; a>= 1; a--)
    {
      if (a<50 || (b==main_i && a==main_j))
	{
	seg_list[c][0]=a;
	seg_list[c][1]=b;
	seg_list[c][2]=MIN((l1-a), (l2-b));
	c++;
	}
      }
  
  
  for ( b=2,a=1; b<= l2; b++, c++)
    {
      if (b<50 || (b==main_i && a==main_j))
	{
	  seg_list[c][0]=a;
	  seg_list[c][1]=b;
	  seg_list[c][2]=MIN((l1-a), (l2-b));
	}
      }


  seg_list[c][0]=FORBIDEN;
  
  Start->i=0 ;Start->j=0 ; Start->st=H->start;Start->sc=0;
  End->i  =l1;  End->j=l2; End  ->st=H->end;
  Start=RviterbiD_hmm (A, ns, ls, H, CL, Start,End,seg_list);
  
    
  path=traceback (A, ns, ls, H, CL, Start,NULL, NULL);

  

  A=viterbipath2aln (A,ns,ls,path, H); 
  
  viterbiD_hmm (NULL, ns, ls, H, CL, Start,End, seg_list);
  free_Hmm (H);
  free_int (A->pos, -1);
  free_arrayN((void **)seg_list, 2);
    
  A->pos=NULL;
  return A->score_aln;
}


double lu_RviterbiD_hmm (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL,MatState *S, MatState *E, int **seg_list)
{
  HmmState *S1;
  double sc, sc2,e, t;
  static MatState *cS=NULL, *cE=NULL;
  double min, max;
  max=MAX((E->i-S->i), (E->j-S->j));
  min=MIN((E->i-S->i), (E->j-S->j));
  
  
  if ( S->sc==H->forbiden) return H->forbiden;
  else if (min==0)
    {
      e=hmm_add(S->sc,H->T[S->st][E->st]);
      if (  H->T[E->st][E->st]!=H->forbiden)e=hmm_add(e, (max-1)*H->T[E->st][E->st]);
      if ( (H->S[E->st])->em!=H->forbiden)  e=hmm_add(e, max    *(H->S[E->st])->em );
      return e;
    }
  else if ( min>0 && max>1)
    { 
      
      fprintf ( stderr, "\nWarning: Disjoined Diagonals");
      DisplayMatState (S, "\n\tS");
      DisplayMatState (E, "\n\tE");
      
      
      cS=CopyMatState ( S,cS);
      cE=CopyMatState ( E,cE);
      cE->sc=H->forbiden;
      viterbiD_hmm (A,ns,ls, H,CL,cS, cE, NULL);
      sc2=cE->sc;
     
      return sc2;
    }
  else
    {
      S1=H->S[E->st];
      t=H->T[S->st][E->st];
      e=(S1->em==H->forbiden)?S1->em_func (A, A->pos, ns[0], ls[0],E->i-1, A->pos,ns[1], ls[1], E->j-1, CL):S1->em;
      sc=hmm_add(S->sc,t);
      sc=hmm_add(sc,e);
      return sc;
    }
  return H->forbiden;
}


MatState* RviterbiD_hmm (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL,MatState *S, MatState *E, int **seg_list)
{
  MatState *MS, *ME;
  MS=S;
  ME=E;
  
  viterbiD_hmm (A,ns, ls,H, CL, S, E, seg_list);
  
  
  if ( S->n==E)return S;
  if ( E->sc==H->forbiden)
	{
	  DisplayHmm (H);
	  fprintf ( stderr, "\nERROR: The Requested Model (Cf Above) Cannot Produce the Pair-Alignment\nYou must allow extra possible transitions\n[FATAL:%s]", PROGRAM);
	  myexit ( EXIT_FAILURE);
	}
  E=S->n;
  
  while (S!=ME)
    {
      int d1, d2, align;
      d1=MinDeltaMatState(S,E);
      d2=MaxDeltaMatState(S,E);
      align=((d1==1 && d2==1) || ( d1==0))?0:1;
      if (align)RviterbiD_hmm (A,ns, ls,H, CL,S,E, seg_list);
      S=E;
      E=S->n;
    }
  return MS;
}

#define Dim(i,j) (i*LenJ+j)
MatState* viterbiD_hmm (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL, MatState *S,MatState *E, int **seg_list)
{
  int current, previous,row, prow;
  double v;
  int a,b,i,j,pi,pj, s, k;
  int start_i, start_j, end_i, end_j, l1, l2;
  HmmState *S1, *S2;
  MatState *CC, *PCC,*tS, *tE, *mark=NULL;
  int midpoint;
  
  int dj;
  int dc;
  int *jlist=NULL;
  static int **main_jlist;
  static MatState ***M;
  static int *toclean;
  int ntoclean;
  static int LenJ, LenI;
  int MaxDelta=50, DeltaI, DeltaJ;
  int mode;

  DisplayMatState (S, "\n\tS");
  DisplayMatState (E, "\n\tE");
  
  
  if ( A==NULL)
    {
      free_arrayN((void **)main_jlist, 2);main_jlist=NULL;

      for ( a=0; a<2; a++)memset(M[a],0,LenJ*LenI*sizeof (MatState*));
      free_arrayN((void **)M, 3);M=NULL;
      vfree (toclean);
      ManageMatState ( FREE_ALL, NULL);
      return NULL;
    }

  
  if ( MatStateAreIdentical ( S, E))return NULL;
  l1=strlen (A->seq_al[ls[0][0]]);l2=strlen (A->seq_al[ls[1][0]]);
  
  midpoint=S->i+((E->i-S->i)/2);
  DeltaI=E->i-S->i;
  
  
  start_i=S->i;end_i=E->i;start_j=S->j;end_j=E->j;
  current=0;previous=1;
 
  
  if ( !M)
    {
      LenI=l2+1;
      LenJ=H->nS;
      M=(MatState***)declare_arrayN(3, sizeof ( MatState),2,LenI*LenJ,0);
      toclean=(int*)vcalloc ( LenI*LenJ, sizeof (int));
    }
  
  if ( !main_jlist)main_jlist= seglist2table(seg_list, l1, l2);
  

  /*MAKE THE VITERBI FROM S(tart) to E(nd)*/
  mark=ManageMatState ( MARK, mark);
  mode=(!seg_list)?ALL:seg_list[0][0];
  
  for (ntoclean=0,i=start_i; i<=end_i; i++)
    {
      row=current;
      
      if ( mode==ALL)jlist=main_jlist[0];
      else if ( mode==DIAGONALS)jlist=(i==0)?main_jlist[0]:main_jlist[1];
      else if ( mode==SEGMENTS)  jlist=main_jlist[i+2];
      
     
      for ( dj=1; dj<=jlist[0]; dj++)
	{
	  DeltaJ=((FABS(dj-i))<MaxDelta)?1:MaxDelta+1; /*Keep Pointers on a 2*DeltaMax Band around the main diagonal*/

	  dc=(mode==DIAGONALS && dj!=1)?i:0;/*Make sure the diagonal Mode uses the sides of the array*/
	  j=jlist[dj]+dc;
	  
	  
	  if ( j<start_j || j>end_j)continue;
	  for ( s=H->nS-1;s>=0; s--)
	    {
	      S1=H->S[s];pi=i-S1->DI;prow=S1->DI;
	      
	      if ( S1->DI && S1->DJ){pj=j-S1->DJ;}
	      else if ( !S1->DJ)pj=j;
	      else if ( dj>1)pj=jlist[dj-S1->DJ]+dc;
	      else pj=-1;
	      
	      if (!M[row][Dim(j,s)])toclean[ntoclean]=Dim(j,s);
	      
	      CC=M[row][Dim(j,s)]=CopyMatState(NULL, M[row][Dim(j,s)]);
	      CC->i=i; CC->j=j; CC->st=s;PCC=NULL;
	      	      
	      if (i==start_i && j==start_j && s==S->st){CC=CopyMatState(S,CC);}
	      else if ( i==end_i && j==end_j && E->st!=H->end && s!=E->st)CC->sc=H->forbiden;
	      else if ( pi<start_i || pj<start_j)        {CC->sc=H->forbiden;}
	      else 
		{
		  for (k=1; k<=H->fromM[S1->state][0]; k++)
		    {
		      S2=H->S[H->fromM[s][k]];
		      
		      if ( S1->DI && S1->DJ)PCC=M[prow][Dim((j-S1->DJ),(S2->state))];
		      else PCC=M[prow][Dim((jlist[dj-S1->DJ]+dc),(S2->state))];
		    
		      if ( !PCC)PCC=NULL;
		      else if      ( pi+pj!=0 && S2->state==H->start);
		      else if ( !(pi==l1 && pj==l2) && s==H->end);
		      else
			{
			  v=lu_RviterbiD_hmm(A,ns, ls, H, CL,PCC,CC, NULL);
			  if ( v!=H->forbiden && (CC->sc==H->forbiden || v> CC->sc)){CC->sc=v; CC->pst=S2->state;CC->p=PCC;}
			}
		    }
		}
	      if (CC->sc==H->forbiden);
	      else if (i==midpoint || DeltaI<=MaxDelta||DeltaJ<=MaxDelta ||(i==start_i && j==start_j && s==S->st) )
		{
		  CC->m=(CC->p)?(CC->p)->m:NULL;
		  PCC=CopyMatState(CC,NULL);
		  PCC->m=CC->m;CC->m=PCC;
		}
	      else CC->m=(CC->p)?(CC->p)->m:NULL;
	    }
	}
      prow=previous;
      for ( dj=1; dj<=jlist[0] && i!=end_i; dj++)
	{
	  dc=(mode==DIAGONALS && dj!=1)?i:0;
	  j=jlist[dj]+dc;
	  if ( j<start_j || j>end_j)continue;
	  
	  for ( s=H->nS-1;s>=0; s--)
	    {
	      
	      CC=(M[prow][Dim(j,s)]);M[prow][Dim(j,s)]=M[row][Dim(j,s)];M[row][Dim(j,s)]=CC;
	      if (!M[row][Dim(j,s)])toclean[ntoclean++]=Dim(j,s);
	      if (M[prow][Dim(j,s)]) M[row ][Dim(j,s)]=CopyMatState ( M[prow][Dim(j,s)], M[row][Dim(j,s)]);
	      
	    }
	} 
      
    }  
  
  mark=ManageMatState ( MARK,mark);
  row=current;
  
  
  if ( E->st==H->end || E->st==H->forbiden){E=CopyMatState ((M[row][Dim(end_j,E->st)]),E);}

  
  

  PCC=CopyMatState (M[row][Dim(end_j,E->st)], NULL);
  
  if ( MatStateAreIdentical(PCC,PCC->m))PCC=PCC->m;
    tS=tE=PCC;
  while (PCC->m)
    {
      tS=CopyMatState (PCC->m,NULL); tS->n=PCC; PCC->p=tS;PCC=tS;
    }
    
  if (tS==tE);
  else
    {
      S->n=tS->n; (S->n)->p=S;
      E->p=tE->p; (E->p)->n=E;
    }

  ManageMatState ( FREE_MARK,mark);
  
  
  for ( a=0; a<ntoclean; a++)
    {
      for ( b=0; b< 2; b++){M[b][toclean[a]]=NULL;}
      toclean[a]=0;
    }
  
  while (S && S->p!=E){S->m=NULL;S=S->n;}/*Clean the memory of the rturned Cells*/
  return NULL;
}

/*********************************************************************************/
/*                                                                               */
/*                                                                               */
/*                      HMM Viterbi Diagonals  GLOBAL/LOCAL                                  */
/*                                                                               */
/*                                                                               */
/*********************************************************************************/

int viterbiDGL_pair_wise (Alignment *A,int*ns, int **ls,Constraint_list *CL)
{
  int l1, l2;
  int *path;
  Hmm * H;
  MatState *Start;
  MatState *End;
  int **seg_list;
  int a, b, c;
  int main_i;
  int main_j;

  
  A->pos=aln2pos_simple( A, -1, ns, ls);
  
  Start=ManageMatState ( DECLARE, NULL);
  End=ManageMatState ( DECLARE, NULL);
  H=define_simple_model (CL);
  l1=strlen (A->seq_al[ls[0][0]]);
  l2=strlen (A->seq_al[ls[1][0]]);

  main_i=MAX(1,(l2-l1)+1);
  main_j=MAX(1,(l1-l2)+1);
  
  seg_list=(int**)declare_arrayN(2, sizeof (int), l1+l2+3, 3);
  seg_list[0][0]=DIAGONALS;

  
  c=1;
  for ( b=1,a=l1; a>= 1; a--)
    {
      if (a<50 || (b==main_i && a==main_j))
	{
	seg_list[c][0]=a;
	seg_list[c][1]=b;
	seg_list[c][2]=MIN((l1-a), (l2-b));
	c++;
	}
      }
  
  
  for ( b=2,a=1; b<= l2; b++, c++)
    {
      if (b<50 || (b==main_i && a==main_j))
	{
	  seg_list[c][0]=a;
	  seg_list[c][1]=b;
	  seg_list[c][2]=MIN((l1-a), (l2-b));
	}
      }


  seg_list[c][0]=FORBIDEN;
  
  Start->i=0 ;Start->j=0 ; Start->st=H->start;Start->sc=0;
  End->i  =l1;  End->j=l2; End  ->st=H->end;
  Start=RviterbiDGL_hmm (A, ns, ls, H, CL, Start,End,seg_list);
  
    
  path=traceback (A, ns, ls, H, CL, Start,NULL, NULL);

  

  A=viterbipath2aln (A,ns,ls,path, H); 
  
  viterbiD_hmm (NULL, ns, ls, H, CL, Start,End, seg_list);
  free_Hmm (H);
  free_int (A->pos, -1);
  free_arrayN((void **)seg_list, 2);
    
  A->pos=NULL;
  return A->score_aln;
}


double lu_RviterbiDGL_hmm (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL,MatState *S, MatState *E, int **seg_list)
{
  HmmState *S1;
  double sc, sc2,e, t;
  static MatState *cS=NULL, *cE=NULL;
  double min, max;
  max=MAX((E->i-S->i), (E->j-S->j));
  min=MIN((E->i-S->i), (E->j-S->j));
  
  
  
  
  if ( S==NULL || E==NULL || S->sc==H->forbiden) return H->forbiden;
  else if ( S->st==H->start) return 0;
  else if ( E->st==H->end) return S->sc;
  else if (min==0)
    {
      e=hmm_add(S->sc,H->T[S->st][E->st]);
      if (  H->T[E->st][E->st]!=H->forbiden)e=hmm_add(e, (max-1)*H->T[E->st][E->st]);
      if ( (H->S[E->st])->em!=H->forbiden)  e=hmm_add(e, max    *(H->S[E->st])->em );
      return e;
    }
  else if ( min>0 && max>1)
    { 
      
      fprintf ( stderr, "\nWarning: Disjoined Diagonals");
      DisplayMatState (S, "\n\tS");
      DisplayMatState (E, "\n\tE");
      
      
      cS=CopyMatState ( S,cS);
      cE=CopyMatState ( E,cE);
      cE->sc=H->forbiden;
      viterbiD_hmm (A,ns,ls, H,CL,cS, cE, NULL);
      sc2=cE->sc;
     
      return sc2;
    }
  else
    {
      S1=H->S[E->st];
      t=H->T[S->st][E->st];
      e=(S1->em==H->forbiden)?S1->em_func (A, A->pos, ns[0], ls[0],E->i-1, A->pos,ns[1], ls[1], E->j-1, CL):S1->em;
      sc=hmm_add(S->sc,t);
      sc=hmm_add(sc,e);
      return sc;
    }
  return H->forbiden;
}


MatState* RviterbiDGL_hmm (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL,MatState *S, MatState *E, int **seg_list)
{
  MatState *MS, *ME;
  MS=S;
  ME=E;
  
  viterbiDGL_hmm (A,ns, ls,H, CL, S, E, seg_list);
  
  
  if ( S->n==E)return S;
  if ( E->sc==H->forbiden)
	{
	  DisplayHmm (H);
	  fprintf ( stderr, "\nERROR: The Requested Model (Cf Above) Cannot Produce the Pair-Alignment\nYou must allow extra possible transitions\n[FATAL:%s]", PROGRAM);
	  myexit ( EXIT_FAILURE);
	}
  E=S->n;
  
  while (S!=ME)
    {
      int d1, d2, align;
      d1=MinDeltaMatState(S,E);
      d2=MaxDeltaMatState(S,E);
      align=((d1==1 && d2==1) || ( d1==0))?0:1;
      if (align)RviterbiDGL_hmm (A,ns, ls,H, CL,S,E, seg_list);
      S=E;
      E=S->n;
    }
  return MS;
}


#define Dim(i,j) (i*LenJ+j)
MatState* viterbiDGL_hmm (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL, MatState *S,MatState *E, int **seg_list)
{
  int current, previous,row, prow;
  double v;
  int a,i,j,pi,pj, s, k;
  int start_i, start_j, end_i, end_j, l1, l2;
  HmmState *S1, *S2;
  MatState *CC, *PCC,*tS, *tE,*bestE,*bestS, *mark=NULL;
  int midpoint;
  
  int dj;
  int dc;
  int *jlist=NULL;
  static int **main_jlist;
  static MatState ***M;
  static int *toclean;
  int ntoclean;
  static int LenJ, LenI;
  int MaxDelta=50, DeltaI, DeltaJ;
  int mode;

  

  DisplayMatState (S, "\n\tS");
  DisplayMatState (E, "\n\tE");
  
  
  if ( A==NULL)
    {
      free_arrayN((void **)main_jlist, 2);main_jlist=NULL;

      for ( a=0; a<2; a++)memset(M[a],0,LenJ*LenI*sizeof (MatState*));
      free_arrayN((void **)M, 3);M=NULL;
      vfree (toclean);
      ManageMatState ( FREE_ALL, NULL);
      return NULL;
    }

  
  if ( MatStateAreIdentical ( S, E))return NULL;
  l1=strlen (A->seq_al[ls[0][0]]);l2=strlen (A->seq_al[ls[1][0]]);
  
  midpoint=S->i+((E->i-S->i)/2);
  DeltaI=E->i-S->i;
  
  
  start_i=S->i;end_i=E->i;start_j=S->j;end_j=E->j;
  current=0;previous=1;
 
  
  if ( !M)
    {
      LenI=l2+1;
      LenJ=H->nS;
      M=(MatState***)declare_arrayN(3, sizeof ( MatState),2,LenI*LenJ,0);
      toclean=(int*)vcalloc ( LenI*LenJ, sizeof (int));
    }
  
  if ( !main_jlist)main_jlist= seglist2table(seg_list, l1, l2);
  

  /*MAKE THE VITERBI FROM S(tart) to E(nd)*/
  mark=ManageMatState ( MARK, mark);
  mode=(!seg_list)?ALL:seg_list[0][0];
  bestE=CopyMatState (E, NULL);
  bestS=CopyMatState (NULL, NULL);
  for (ntoclean=0,i=start_i; i<=end_i; i++)
    {
      row=current;
      
      if ( mode==ALL)jlist=main_jlist[0];
      else if ( mode==DIAGONALS)jlist=(i==0)?main_jlist[0]:main_jlist[1];
      else if ( mode==SEGMENTS)  jlist=main_jlist[i+2];
      
     
      for ( dj=1; dj<=jlist[0]; dj++)
	{
	  DeltaJ=(FABS(dj-i)<MaxDelta)?1:MaxDelta+1; /*Keep Pointers on a 2*DeltaMax Band around the main diagonal*/

	  dc=(mode==DIAGONALS && dj!=1)?i:0;/*Make sure the diagonal Mode uses the sides of the array*/
	  j=jlist[dj]+dc;
	  
	  
	  if ( j<start_j || j>end_j)continue;
	  for ( s=H->nS-1;s>=0; s--)
	    {
	      if ( s==S->st)continue;
	      S1=H->S[s];pi=i-S1->DI;prow=S1->DI;
	      
	      if ( S1->DI && S1->DJ){pj=j-S1->DJ;}
	      else if ( !S1->DJ)pj=j;
	      else if ( dj>1)pj=jlist[dj-S1->DJ]+dc;
	      else pj=-1;
	      
	      if (!M[row][Dim(j,s)])toclean[ntoclean]=Dim(j,s);
	      
	      CC=M[row][Dim(j,s)]=CopyMatState(NULL, M[row][Dim(j,s)]);
	      CC->i=i; CC->j=j; CC->st=s;PCC=NULL;
	      	      
	      if (i==start_i && j==start_j && s==S->st){CC=CopyMatState(S,CC);}
	      else if ( s==S->st);
	      else if ( i==end_i && j==end_j && E->st!=H->end && s!=E->st)CC->sc=H->forbiden;
	      
	      else if ( pi<start_i || pj<start_j)        {CC->sc=H->forbiden;}
	      else 
		{
		  for (k=1; k<=H->fromM[S1->state][0]; k++)
		    {
		      S2=H->S[H->fromM[s][k]];
		      
		      if ( S1->DI && S1->DJ)PCC=M[prow][Dim((j-S1->DJ),(S2->state))];
		      else PCC=M[prow][Dim((jlist[dj-S1->DJ]+dc),(S2->state))];
		    
		      if ( S2->state==H->start){PCC=bestS;PCC->st=0;PCC->sc=0;PCC->m=PCC->n=PCC->p=NULL;}
		      
		      v=lu_RviterbiDGL_hmm(A,ns, ls, H, CL,PCC,CC, NULL);
		      if ( v!=H->forbiden && (CC->sc==H->forbiden || v> CC->sc)){CC->sc=v; CC->pst=S2->state;CC->p=PCC;}
		    }
		}
	      if ( CC->sc==H->forbiden);
	      else if ( bestE->sc==H->forbiden || bestE->sc>CC->sc) 
		{
		  bestE=CopyMatState(CC, bestE);
		  bestE->m=(CC->p)->m;
		}
	      else if (CC->p && (CC->p)->st==H->start)
		{
		  CC->m=CopyMatState (CC->p, NULL);		  
		}
	      else if (i==midpoint || DeltaI<=MaxDelta||DeltaJ<=MaxDelta ||(i==start_i && j==start_j && s==S->st) )
		{
		  CC->m=(CC->p)?(CC->p)->m:NULL;
		  PCC=CopyMatState(CC,NULL);
		  PCC->m=CC->m;CC->m=PCC;
		}
	      else CC->m=(CC->p)?(CC->p)->m:NULL;
	    }
	}
      prow=previous;
      for ( dj=1; dj<=jlist[0] && i!=end_i; dj++)
	{
	  dc=(mode==DIAGONALS && dj!=1)?i:0;
	  j=jlist[dj]+dc;
	  if ( j<start_j || j>end_j)continue;
	  
	  for ( s=H->nS-1;s>=0; s--)
	    {
	      
	      CC=(M[prow][Dim(j,s)]);M[prow][Dim(j,s)]=M[row][Dim(j,s)];M[row][Dim(j,s)]=CC;
	      /*if (!M[row][Dim(j,s)])toclean[ntoclean++]=Dim(j,s);*/
	      if (M[prow][Dim(j,s)]) M[row ][Dim(j,s)]=CopyMatState ( M[prow][Dim(j,s)], M[row][Dim(j,s)]);
	      
	    }
	} 
      
    }  
  
  mark=ManageMatState ( MARK,mark);
  row=current;
  
  
  if ( E->st==H->end || E->st==H->forbiden){E=CopyMatState ((M[row][Dim(end_j,E->st)]),E);}
  PCC=CopyMatState (bestE, NULL);
  
  if ( MatStateAreIdentical(PCC,PCC->m))PCC=PCC->m;
    tS=tE=PCC;
  while (PCC->m)
    {
      tS=CopyMatState (PCC->m,NULL); tS->n=PCC; PCC->p=tS;PCC=tS;
    }
    
  if (tS==tE);
  else
    {
      CopyMatState ( tS, S);
      CopyMatState ( tE, E);
    }
  ManageMatState ( FREE_MARK,mark);
  
  for ( a=0; a<2; a++)memset(M[a],0,LenJ*LenI*sizeof (MatState*));
  
  while (S && S->p!=E){S->m=NULL;S=S->n;}/*Clean the memory of the rturned Cells*/
  return NULL;
}


/*********************************************************************************/
/*                                                                               */
/*                                                                               */
/*                      HMM Viterbi Diagonals  PROCESSING                        */
/*                                                                               */
/*                                                                               */
/*********************************************************************************/
int **seglist2table ( int **seglist,int l1, int l2)
  {
    int **valuesT;
    int *bvalues;
    int line, a,si, sj, ei, j, c;
    
    /*All: 0*/
    valuesT=(int**)vcalloc ((l1+2)+3, sizeof (int*));
    valuesT[0]=(int*)vcalloc (l2+2, sizeof (int));
    for (a=0; a<=l2; a++)valuesT[0][++valuesT[0][0]]=a;

    if ( !seglist) return valuesT;
    /*Diagonals: 1*/
    valuesT[1]=(int*)vcalloc (l1+l2+2, sizeof (int));
    bvalues=(int*)vcalloc (l1+l2+2, sizeof (int));
    c=1;
    while (seglist[c][0]!=FORBIDEN)
      {
	
	si=seglist[c][0];
	sj=seglist[c][1];
	
	bvalues[(sj-si)+l1]=1;
	c++;
      }
    valuesT[1][++valuesT[1][0]]=0;
    for (a=0; a<=(l1+l2); a++)
      {
	if (bvalues[a])
	  {
	    valuesT[1][++valuesT[1][0]]=a-l1;
	  }
	
      }
    vfree (bvalues);

    /*Segments: 2*/
    valuesT[2]=(int*)vcalloc (l2+2, sizeof (int));
    for (a=0; a<=l2; a++)valuesT[2][++valuesT[2][0]]=a;
    
    bvalues=(int*)vcalloc (l2+2, sizeof (int));
    for ( line=1; line<=l1; line++)
      {
	bvalues[0]=c=0;
	bvalues[++bvalues[0]]=0;
	while (seglist[c][0]!=FORBIDEN)
	  {
		si=seglist[c][0];
		ei=si+seglist[c][2];
		sj=seglist[c][1];
		j=sj+(line-si);
		if ( line<si || line>ei);
		else if (j>=0 && j<=l2 && seglist[c][2])
		  {
		    bvalues[++bvalues[0]]=j;
		  }
		c++;
	  }
	valuesT[line+2]=(int*)vcalloc (bvalues[0]+1, sizeof (int));
	for ( a=0; a<=bvalues[0]; a++) valuesT[line+2][a]=bvalues[a];
      }
    vfree (bvalues);
    return valuesT;
  }


  
/*********************************************************************************/
/*                                                                               */
/*                                                                               */
/*                      HMM modeling                                             */
/*                                                                               */
/*                                                                               */
/*********************************************************************************/

Hmm* declare_hmm(int n)
  {
    Hmm *H;
    int a, b;
    
    H=(Hmm*)vcalloc (1, sizeof (Hmm));
    H->nS=n;
    H->S=(HmmState**)vcalloc (H->nS, sizeof (HmmState*));
    for (a=0; a<H->nS; a++)
      {
	H->S[a]=(HmmState*)vcalloc (1, sizeof (HmmState));
	(H->S[a])->em2=(double*)vcalloc (MAX_EMISSION, sizeof (double));
	
	(H->S[a])->T=(StateTrans**)vcalloc ( H->nS, sizeof (StateTrans*));
	for ( b=0; b< H->nS; b++)
	  (H->S[a])->T[b]=(StateTrans*)vcalloc (1, sizeof (StateTrans));
      }
    return H;
  }

Hmm* free_Hmm(Hmm*H)
  {
    int a, b;
    
    H=(Hmm*)vcalloc (1, sizeof (Hmm));
    free_double (H->T, -1);
    free_int ( H->fromM, -1);
    free_int ( H->toM, -1);
    
    for (a=0; a< H->nS; a++)
      {
	
	for ( b=0; b< H->nS; b++)
	  {
	    vfree ((H->S[a])->em2);
	    vfree((H->S[a])->T[b]);
	  }
	vfree((H->S[a])->T);
	vfree(H->S[a]);
      }
    vfree (H->S);
    vfree (H);
    return NULL;
  }

void DisplayHmm ( Hmm *H)
{
  int a, b;
  HmmState *S1, *S2;
  
  for ( a=0; a< H->nS; a++)
    {
      S1=H->S[a];
      fprintf ( stderr, "\nSTATE %d: %s\n",S1->state,S1->name);
      fprintf ( stderr, "\n\tDI %d", S1->DI);
      fprintf ( stderr, "\n\tDJ %d", S1->DJ);
      fprintf ( stderr, "\n\tE  %f", (float)S1->em);
            
      fprintf ( stderr, "\nReached FROM: ");
      for ( b=1; b<=H->fromM[a][0]; b++)
	{
	  S2=H->S[H->fromM[a][b]];
	  fprintf ( stderr, "[ %s %f ] ", S2->name, H->T[S2->state][S1->state]);
	}
      fprintf ( stderr, "\nGoes TO: ");
      for ( b=1; b<=H->toM[a][0]; b++)
	{
	  S2=H->S[H->toM[a][b]];
	  fprintf ( stderr, "[ %s %f ] ", S2->name, H->T[S1->state][S2->state]);
	}
    }
  return;
}
Hmm * bound_hmm ( Hmm *H)
{
  int a, b, c;
  char **name;
  HmmState *S;
  
  name=declare_char(H->nS, 100);
  H->T=declare_double ( H->nS, H->nS);
  
  for ( a=0; a< H->nS; a++)
    {
      sprintf ( name[a], "%s", (H->S[a])->name);
      H->order=MAX(H->order, (H->S[a])->DI);
      H->order=MAX(H->order, (H->S[a])->DJ);
    }
  
  for ( a=0; a< H->nS; a++)for (b=0; b<H->nS; b++)H->T[a][b]=H->forbiden;
  for (a=0; a< H->nS; a++)
    {
      S=H->S[a];
      for ( b=0; b< S->nT; b++)
	{
	  c=name_is_in_list ((S->T[b])->name, name, H->nS, 100);
	  if ( c!=-1)H->T[a][c]=(S->T[b])->tr;
	}
    }
  
  /*Bound the model:
    bM[state][0]=n_allowed transitions
    bM[state][1]=first allowed transition
  */

  H->toM=declare_int ( H->nS, H->nS);
  H->fromM=declare_int ( H->nS, H->nS);
  
  for ( a=0; a< H->nS; a++)
    for ( b=0; b< H->nS; b++)
      {
	if ( H->T[a][b]!=H->forbiden )
	  {
	    {H->fromM[b][0]++; H->fromM[b][H->fromM[b][0]]=a;}
	    {H->toM[a][0]++; H->toM[a][H->toM[a][0]]=b;}
	  }
      }
  for ( a=0; a< H->nS; a++)
    {
      if (( H->S[a])->em!=H->forbiden)( H->S[a])->em*=SCORE_K;
      for ( b=0; b< H->nS; b++)
	if ( H->T[a][b]!=H->forbiden)H->T[a][b]*=SCORE_K;
    }
  free_arrayN((void**)name, 2);
  return H;
}

/*********************************************************************************/
/*                                                                               */
/*                                                                               */
/*                      Memory Management                                        */
/*                                                                               */
/*                                                                               */
/*********************************************************************************/

MatState * ManageMatState(int Mode, MatState *C)
{
  static MatState *Fheap;
  static MatState *Aheap;
  MatState *Cmark, *Pmark;
  static int alloc, free;
  if (!Fheap || Fheap->Hp==NULL)
    {
      int c=0;
      int extension=1000;
      if (!Fheap){Fheap=(MatState*)vcalloc (1, sizeof (MatState));Fheap->free=1;free++;}
      if (!Aheap)Aheap=(MatState*)vcalloc (1, sizeof (MatState));
      while ( c!=extension)
	{
	  C=(MatState*)vcalloc ( 1, sizeof (MatState));
	  C->free=1;Fheap->Hn=C;C->Hp=Fheap;
	  Fheap=C;
	  c++;
	  free++;
	}
    }
  
  if ( Mode==DECLARE)
    {
     
      C=Fheap;
      Fheap=Fheap->Hp;
      C->Hn=C->Hp=NULL;
      if ( Aheap){Aheap->Hn=C;C->Hp=Aheap;Aheap=C;}
      else Aheap=C;
      alloc++;
      free--;
      C->free=0;
      C=CopyMatState(NULL, C);
      return C;
    }
  else if ( Mode==FREE)
    {
      if ( !C || C->free==1);
      else
	{
	  C=CopyMatState(NULL, C);
	  C->free=1;
	  if (C->Hp==NULL && C==Aheap)crash ("");
	  if (C==Aheap)Aheap=C->Hp;
	  if (C->Hn){(C->Hn)->Hp=C->Hp;}
	  if (C->Hp){(C->Hp)->Hn=C->Hn;}
	  C->Hp=C->Hn=NULL;
	  Fheap->Hn=C;C->Hp=Fheap;
	  Fheap=C;
	  alloc--;
	  free++;
	}
      return NULL;
    }
  else if ( Mode==FREE_ALL)
    {
      while ( Aheap)
	{
	  C=Aheap->Hp;
	  vfree (Aheap);
	  Aheap=C;
	}
      while ( Fheap)
	{
	  C=Fheap->Hp;
	  vfree (Fheap);
	  Fheap=C;
	}
    }
  else if ( Mode==INFO)
    {
      fprintf ( stderr, "\nAllocated: %d Free %d", alloc, free);
    }
  else if ( Mode==MARK)
    {
      
      if (C==NULL);
      else {C->Mn=Aheap;Aheap->Mp=C;}
      
      return Aheap;
    }
  else if ( Mode==UNMARK)
    {
      Pmark=Cmark=NULL;
    }
  else if ( Mode == FREE_MARK)
    {
      Cmark=C;
      Pmark=C->Mp;
      
      if ( Cmark==Pmark)return NULL;
      else if ( Cmark==Aheap)
	{Aheap=Pmark;C=Pmark->Hn;Pmark->Hn=NULL;}
      else
	{
	  (Cmark->Hn)->Hp=Pmark;
	  C=Pmark->Hn;
	  Pmark->Hn=Cmark->Hn;
	}
      
      Fheap->Hn=C;
      C->Hp=Fheap;
      Fheap=Cmark;
      Fheap->Hn=NULL;
  
      C=Fheap;
      while (C && !C->free)
	{
	  free++;alloc--;
	  C->free=1;
	  C=C->Hp;
	}
        
    }
  return NULL;
}


MatState* CopyMatState ( MatState*I, MatState*O)
{
  if (O==NULL || O->free==1) O=ManageMatState(DECLARE, NULL);
  if (I==NULL || I->free==1)I=NULL;
  O->i  =(I)?I->i:0;
  O->j  =(I)?I->j:0;
  O->st =(I)?I->st:FORBIDEN;
  O->pst=(I)?I->pst:FORBIDEN;
  O->sc =(I)?I->sc:FORBIDEN;
  O->n  =(I)?I->n:NULL;
  O->p  =(I)?I->p:NULL;
  O->m  =(I)?I->m:NULL;
  O->s  =(I)?I->m:NULL;
  
  return O;
}

/*********************************************************************************/
/*                                                                               */
/*                                                                               */
/*                      Comparisons                                              */
/*                                                                               */
/*                                                                               */
/*********************************************************************************/
int MaxDeltaMatState (MatState*S, MatState*E)
{
  if ( !S || !E) return -1;
  else return MAX((E->i-S->i),(E->j-S->j));
}
int MinDeltaMatState (MatState*S, MatState*E)
{
  if ( !S || !E) return -1;
  else return MIN((E->i-S->i),(E->j-S->j));
}
int MatStateAreIdentical (MatState*I, MatState*O)
{
  if ( !I || !O)return 0;

  if ( I->i!=O->i)return 0;
  if ( I->j!=O->j)return 0;
  if ( I->st!=O->st)return 0;
  return 1;
}
  


Hmm* define_probcons_model(Constraint_list *CL)
{
  Hmm *H;
  double gop=-10;
  double gep=-1;
  double lgop=-100;
  double lgep=-100;
  double freeT=0;
  int n=0;
  HmmState *S;
  
  
  H=declare_hmm(7);
  H->freeT=freeT=0;
 
  H->forbiden=FORBIDEN;
  H->start=START_STATE;
  H->end=END_STATE;
 
  /*define START*/
  S=H->S[n];  
  sprintf (S->name, "START"); S->state=n;
  
  S->DI=0;
  S->DJ=0;
  S->em=freeT;
 
  sprintf ( (S->T[S->nT])->name, "MATCH") ;(S->T[S->nT])->tr=freeT  ;S->nT++;
  sprintf ( (S->T[S->nT])->name, "INSERT");(S->T[S->nT])->tr=freeT  ;S->nT++;
  sprintf ( (S->T[S->nT])->name, "DELETE");(S->T[S->nT])->tr=freeT  ;S->nT++;
  sprintf ( (S->T[S->nT])->name, "END")   ;(S->T[S->nT])->tr=freeT  ;S->nT++;
  
  n++;
  /*define END*/
  S=H->S[n];  
  sprintf (S->name, "END"); S->state=n;
  S->DI=0;
  S->DJ=0;
  S->em=freeT;
  n++;

  /*define Match*/
  S=H->S[n];  
  sprintf (S->name, "MATCH"); S->state=n;
  S->DI=1;
  S->DJ=1;
  S->em=H->forbiden;
  S->em_func=CL->get_dp_cost;
  
  sprintf ( (S->T[S->nT])->name, "MATCH") ;(S->T[S->nT])->tr=freeT;S->nT++;
  sprintf ( (S->T[S->nT])->name, "INSERT");(S->T[S->nT])->tr=gop  ;S->nT++;
  sprintf ( (S->T[S->nT])->name, "DELETE");(S->T[S->nT])->tr=gop  ;S->nT++;
  sprintf ( (S->T[S->nT])->name, "END");   (S->T[S->nT])->tr=freeT;S->nT++;
  
  n++;
  
  /*define Insert*/
  S=H->S[n];  
  sprintf (S->name, "INSERT"); S->state=n;
  S->DI=1;
  S->DJ=0;
  S->em=gep;
  sprintf ( (S->T[S->nT])->name, "MATCH") ; (S->T[S->nT])->tr=freeT;S->nT++;
  sprintf ( (S->T[S->nT])->name, "INSERT"); (S->T[S->nT])->tr=freeT;S->nT++;
  sprintf ( (S->T[S->nT])->name, "LINSERT");(S->T[S->nT])->tr=lgop ;S->nT++;
  sprintf ( (S->T[S->nT])->name, "END");    (S->T[S->nT])->tr=-gop;S->nT++;
  
  n++;
  
  /*define LInsert*/
  S=H->S[n];  
  sprintf (S->name, "LINSERT"); S->state=n;
  S->DI=1;
  S->DJ=0;
  S->em=lgep;
  
  sprintf ( (S->T[S->nT])->name, "INSERT")  ;(S->T[S->nT])->tr=freeT;S->nT++;
  sprintf ( (S->T[S->nT])->name, "LINSERT");(S->T[S->nT])->tr=freeT;S->nT++;
   
  n++;
  
  
  /*define Delete*/
  S=H->S[n];  
  sprintf (S->name, "DELETE"); S->state=n;
  S->DI=0;
  S->DJ=1;
  S->em=gep;
  
  sprintf ( (S->T[S->nT])->name, "MATCH")   ;(S->T[S->nT])->tr=freeT;S->nT++;
   sprintf ( (S->T[S->nT])->name, "DELETE")  ;(S->T[S->nT])->tr=freeT;S->nT++;
  sprintf ( (S->T[S->nT])->name, "LDELETE") ;(S->T[S->nT])->tr=lgop ;S->nT++;
  sprintf ( (S->T[S->nT])->name, "END")     ;(S->T[S->nT])->tr=-gop;S->nT++;
  
  n++;
  
  /*define LDelete*/
  S=H->S[n];  
  sprintf (S->name, "LDELETE"); S->state=n;
  S->DI=0;
  S->DJ=1;
  S->em=lgep;
  sprintf ( (S->T[S->nT])->name, "DELETE")  ;(S->T[S->nT])->tr=freeT;S->nT++;
  sprintf ( (S->T[S->nT])->name, "LDELETE");(S->T[S->nT])->tr=freeT;S->nT++;
    
  n++;
  
  
  H=bound_hmm ( H);
  return H;
}

Hmm* define_mnm_model(Constraint_list *CL)
{
  Hmm *H;
  double gop=20;



  double freeT=0;
  int n=0;
  HmmState *S;
  
  
  H=declare_hmm(6);
  H->freeT=freeT=0;
 
  H->forbiden=FORBIDEN;
  H->start=START_STATE;
  H->end=END_STATE;
 
  /*define START*/
  S=H->S[n];  
  sprintf (S->name, "START"); S->state=n;
  
  S->DI=0;
  S->DJ=0;
  S->em=freeT;
 
  sprintf ( (S->T[S->nT])->name, "MATCH") ;(S->T[S->nT])->tr=freeT  ;S->nT++;
  sprintf ( (S->T[S->nT])->name, "INSERT");(S->T[S->nT])->tr=freeT  ;S->nT++;
  sprintf ( (S->T[S->nT])->name, "DELETE");(S->T[S->nT])->tr=freeT  ;S->nT++;
  sprintf ( (S->T[S->nT])->name, "NOMATCH");(S->T[S->nT])->tr=freeT  ;S->nT++;
  sprintf ( (S->T[S->nT])->name, "END")   ;(S->T[S->nT])->tr=freeT  ;S->nT++;
  
  n++;
  /*define END*/
  S=H->S[n];  
  sprintf (S->name, "END"); S->state=n;
  S->DI=0;
  S->DJ=0;
  S->em=freeT;
  n++;

  /*define Match*/
  S=H->S[n];  
  sprintf (S->name, "MATCH"); S->state=n;
  S->DI=1;
  S->DJ=1;
  S->em=H->forbiden;
  S->em_func=CL->get_dp_cost;
  
  sprintf ( (S->T[S->nT])->name, "MATCH")  ;(S->T[S->nT])->tr=freeT;S->nT++;
  sprintf ( (S->T[S->nT])->name, "NOMATCH");(S->T[S->nT])->tr=gop  ;S->nT++;
  sprintf ( (S->T[S->nT])->name, "END");   (S->T[S->nT])->tr=freeT;S->nT++;
  
  n++;
  
  /*define NOMatch*/
  S=H->S[n];  
  sprintf (S->name, "NOMATCH"); S->state=n;
  S->DI=1;
  S->DJ=1;
  S->em=freeT;
  S->em_func=NULL;
  
  sprintf ( (S->T[S->nT])->name, "NOMATCH")  ;(S->T[S->nT])->tr=freeT;S->nT++;
  sprintf ( (S->T[S->nT])->name, "MATCH")  ;(S->T[S->nT])->tr=freeT;S->nT++;  
  sprintf ( (S->T[S->nT])->name, "INSERT");(S->T[S->nT])->tr=freeT  ;S->nT++;
  sprintf ( (S->T[S->nT])->name, "DELETE");(S->T[S->nT])->tr=freeT  ;S->nT++;  
  sprintf ( (S->T[S->nT])->name, "END");   (S->T[S->nT])->tr=freeT;S->nT++;
  
  n++;
  /*define Insert*/
  S=H->S[n];  
  sprintf (S->name, "INSERT"); S->state=n;
  S->DI=1;
  S->DJ=0;
  S->em=freeT;
  sprintf ( (S->T[S->nT])->name, "NOMATCH") ; (S->T[S->nT])->tr=freeT;S->nT++;
  sprintf ( (S->T[S->nT])->name, "INSERT")  ; (S->T[S->nT])->tr=freeT;S->nT++;
  sprintf ( (S->T[S->nT])->name, "END");    (S->T[S->nT])->tr=freeT;S->nT++;
  
  n++;
  
  /*define Delete*/
  S=H->S[n];  
  sprintf (S->name, "DELETE"); S->state=n;
  S->DI=0;
  S->DJ=1;
  S->em=freeT;
  
  sprintf ( (S->T[S->nT])->name, "NOMATCH")   ;(S->T[S->nT])->tr=freeT;S->nT++;
  sprintf ( (S->T[S->nT])->name, "DELETE")  ;(S->T[S->nT])->tr=freeT;S->nT++;
  sprintf ( (S->T[S->nT])->name, "END")     ;(S->T[S->nT])->tr=freeT;S->nT++;
  
  n++;
  
  
  H=bound_hmm ( H);
  return H;
}

Hmm* define_simple_model(Constraint_list *CL)
{
  Hmm *H;
  double gop=-10;
  double gep=-1;
  double freeT=0;
  int n=0;
  HmmState *S;
  
  
  H=declare_hmm(5);
  H->freeT=freeT=0;
 
  H->forbiden=FORBIDEN;
  H->start=START_STATE;
  H->end=END_STATE;
 
  /*define START*/
  S=H->S[n];  
  sprintf (S->name, "START"); S->state=n;
  
  S->DI=0;
  S->DJ=0;
  S->em=freeT;
 
  sprintf ( (S->T[S->nT])->name, "MATCH") ;(S->T[S->nT])->tr=freeT  ;S->nT++;
  sprintf ( (S->T[S->nT])->name, "INSERT");(S->T[S->nT])->tr=freeT  ;S->nT++;
  sprintf ( (S->T[S->nT])->name, "DELETE");(S->T[S->nT])->tr=freeT  ;S->nT++;
  sprintf ( (S->T[S->nT])->name, "END")   ;(S->T[S->nT])->tr=freeT  ;S->nT++;
  
  n++;
  /*define END*/
  S=H->S[n];  
  sprintf (S->name, "END"); S->state=n;
  S->DI=0;
  S->DJ=0;
  S->em=freeT;
  n++;

  /*define Match*/
  S=H->S[n];  
  sprintf (S->name, "MATCH"); S->state=n;
  S->DI=1;
  S->DJ=1;
  S->em=H->forbiden;
  S->em_func=CL->get_dp_cost;
  
  sprintf ( (S->T[S->nT])->name, "MATCH") ;(S->T[S->nT])->tr=freeT;S->nT++;
  sprintf ( (S->T[S->nT])->name, "INSERT");(S->T[S->nT])->tr=gop  ;S->nT++;
  sprintf ( (S->T[S->nT])->name, "DELETE");(S->T[S->nT])->tr=gop  ;S->nT++;
  sprintf ( (S->T[S->nT])->name, "END");   (S->T[S->nT])->tr=freeT;S->nT++;
  
  n++;
  
  /*define Insert*/
  S=H->S[n];  
  sprintf (S->name, "INSERT"); S->state=n;
  S->DI=1;
  S->DJ=0;
  S->em=gep;
  sprintf ( (S->T[S->nT])->name, "MATCH") ; (S->T[S->nT])->tr=freeT;S->nT++;
  sprintf ( (S->T[S->nT])->name, "INSERT"); (S->T[S->nT])->tr=freeT;S->nT++;
  sprintf ( (S->T[S->nT])->name, "DELETE"); (S->T[S->nT])->tr=freeT;S->nT++;
  
  sprintf ( (S->T[S->nT])->name, "END");    (S->T[S->nT])->tr=-gop;S->nT++;
  
  n++;
  
    
  /*define Delete*/
  S=H->S[n];  
  sprintf (S->name, "DELETE"); S->state=n;
  S->DI=0;
  S->DJ=1;
  S->em=gep;
  
  sprintf ( (S->T[S->nT])->name, "MATCH")   ;(S->T[S->nT])->tr=freeT;S->nT++;
  sprintf ( (S->T[S->nT])->name, "DELETE")  ;(S->T[S->nT])->tr=freeT;S->nT++;
  sprintf ( (S->T[S->nT])->name, "INSERT"); (S->T[S->nT])->tr=freeT;S->nT++;
  sprintf ( (S->T[S->nT])->name, "END")     ;(S->T[S->nT])->tr=-gop;S->nT++;
  
  n++;
  
  H=bound_hmm ( H);
  return H;
}

Hmm* define_two_mat_model(Constraint_list *CL)
{
  Hmm *H;
  double gop=-15;
  double gep=-2;
  double lgop=-6;
  double lgep=-1;
  double freeT=0;
  int n=0;
  HmmState *S;
  
  
  H=declare_hmm(8);
  H->freeT=freeT=0;
 
  H->forbiden=FORBIDEN;
  H->start=START_STATE;
  H->end=END_STATE;
 
  /*define START*/
  S=H->S[n];  
  sprintf (S->name, "START"); S->state=n;
  
  S->DI=0;
  S->DJ=0;
  S->em=freeT;
 
  sprintf ( (S->T[S->nT])->name, "MATCH1") ;(S->T[S->nT])->tr=freeT  ;S->nT++;
  sprintf ( (S->T[S->nT])->name, "MATCH2") ;(S->T[S->nT])->tr=freeT  ;S->nT++;
  
  sprintf ( (S->T[S->nT])->name, "INSERT");(S->T[S->nT])->tr=freeT  ;S->nT++;
  sprintf ( (S->T[S->nT])->name, "DELETE");(S->T[S->nT])->tr=freeT  ;S->nT++;
  sprintf ( (S->T[S->nT])->name, "END")   ;(S->T[S->nT])->tr=freeT  ;S->nT++;
  
  n++;
  /*define END*/
  S=H->S[n];  
  sprintf (S->name, "END"); S->state=n;
  S->DI=0;
  S->DJ=0;
  S->em=freeT;
  n++;

  /*define Match*/
  S=H->S[n];  
  sprintf (S->name, "MATCH1"); S->state=n;
  S->DI=1;
  S->DJ=1;
  S->em=H->forbiden;
  S->em_func=get_dp_cost_pam_matrix;
  
  sprintf ( (S->T[S->nT])->name, "MATCH1") ;(S->T[S->nT])->tr=freeT;S->nT++;
  sprintf ( (S->T[S->nT])->name, "INSERT");(S->T[S->nT])->tr=gop  ;S->nT++;
  sprintf ( (S->T[S->nT])->name, "DELETE");(S->T[S->nT])->tr=gop  ;S->nT++;
  sprintf ( (S->T[S->nT])->name, "END");   (S->T[S->nT])->tr=freeT;S->nT++;
  
  n++;
  
  /*define Match*/
  S=H->S[n];  
  sprintf (S->name, "MATCH2"); S->state=n;
  S->DI=1;
  S->DJ=1;
  S->em=H->forbiden;
  S->em_func=get_dp_cost_blosum_matrix;
  
  sprintf ( (S->T[S->nT])->name, "MATCH2") ;(S->T[S->nT])->tr=freeT;S->nT++;
  sprintf ( (S->T[S->nT])->name, "INSERT");(S->T[S->nT])->tr=gop  ;S->nT++;
  sprintf ( (S->T[S->nT])->name, "DELETE");(S->T[S->nT])->tr=gop  ;S->nT++;
  sprintf ( (S->T[S->nT])->name, "END");   (S->T[S->nT])->tr=freeT;S->nT++;
  
  n++;

  /*define Insert*/
  S=H->S[n];  
  sprintf (S->name, "INSERT"); S->state=n;
  S->DI=1;
  S->DJ=0;
  S->em=gep;
  sprintf ( (S->T[S->nT])->name, "MATCH2") ; (S->T[S->nT])->tr=freeT;S->nT++;
  sprintf ( (S->T[S->nT])->name, "MATCH1") ; (S->T[S->nT])->tr=freeT;S->nT++;
  
  sprintf ( (S->T[S->nT])->name, "INSERT"); (S->T[S->nT])->tr=freeT;S->nT++;
  sprintf ( (S->T[S->nT])->name, "LINSERT");(S->T[S->nT])->tr=lgop ;S->nT++;
  sprintf ( (S->T[S->nT])->name, "END");    (S->T[S->nT])->tr=-gop;S->nT++;
  
  n++;
  
  /*define LInsert*/
  S=H->S[n];  
  sprintf (S->name, "LINSERT"); S->state=n;
  S->DI=1;
  S->DJ=0;
  S->em=lgep;
  
  sprintf ( (S->T[S->nT])->name, "INSERT")  ;(S->T[S->nT])->tr=freeT;S->nT++;
  sprintf ( (S->T[S->nT])->name, "LINSERT");(S->T[S->nT])->tr=freeT;S->nT++;
    
  n++;
  
  
  /*define Delete*/
  S=H->S[n];  
  sprintf (S->name, "DELETE"); S->state=n;
  S->DI=0;
  S->DJ=1;
  S->em=gep;
  
  sprintf ( (S->T[S->nT])->name, "MATCH2")   ;(S->T[S->nT])->tr=freeT;S->nT++;
  sprintf ( (S->T[S->nT])->name, "MATCH1")   ;(S->T[S->nT])->tr=freeT;S->nT++;
  
  sprintf ( (S->T[S->nT])->name, "DELETE")  ;(S->T[S->nT])->tr=freeT;S->nT++;
  sprintf ( (S->T[S->nT])->name, "LDELETE") ;(S->T[S->nT])->tr=lgop ;S->nT++;
  sprintf ( (S->T[S->nT])->name, "END")     ;(S->T[S->nT])->tr=-gop;S->nT++;
  n++;
  
  /*define LDelete*/
  S=H->S[n];  
  sprintf (S->name, "LDELETE"); S->state=n;
  S->DI=0;
  S->DJ=1;
  S->em=lgep;
  sprintf ( (S->T[S->nT])->name, "DELETE")  ;(S->T[S->nT])->tr=freeT;S->nT++;
  sprintf ( (S->T[S->nT])->name, "LDELETE");(S->T[S->nT])->tr=freeT;S->nT++;
  n++;
  
  if ( n!=H->nS)
    {
      fprintf ( stderr, "\nERROR in HMM definition [FATAL:%s]", PROGRAM);
      myexit (EXIT_FAILURE);
    }
  
  H=bound_hmm ( H);
  return H;
} 
void DisplayMatState ( MatState *S, char *s)
{
  if ( S==NULL)fprintf ( stderr, "%s: Cell is undefined", s);
  else  fprintf ( stderr, "%s: i=%d j=%d st=%d pst=%d sc=%d Free %d", s, S->i, S->j, S->st, S->pst, (int)S->sc, S->free);
}
void testfunc ( MatState *S, char *s)
{
  if ( S==NULL)return;
  fprintf ( stderr, "\n#### %s ", s);
  while ( S){DisplayMatState ( S,"\n\t");S=S->n;}
  fprintf ( stderr, "\n");
}
  
































#ifdef BACKHERE
	  
	  if ( i>0 && j>0)
	  m=emit_pair_default[alphabetDefault[seq1[a]]][alphabetDefault[seq2[a]]];
	  /*Match*/
	  F[M][i][j]=F[M][i-1][j-1];
	  
	  
	
	M[Match][i][j]=m+log_add3( M[Match][i-step_i][j-step_j],M[I][i-step_i][j],M[D][i][j-step_j]);		
	M[D    ][i][j]=log_add3(gep,M[Match][i       ][j-step_j]+gop,M[D][i       ][j-step_j]);
	M[I    ][i][j]=log_add3(gep,M[Match][i-step_i][j       ]+gop,M[I][i-step_i][j       ]);
	
	/*Long gaps
	M[Match][i][j]=log_add3(M[Match][i][j], M[LI][i-step_i][j],M[LD][i][j-step_j]);
	M[LI    ][i][j]=log_add3(lgep, M[I][i-step_i][j       ]+lgop,M[LI][i-step_i][j      ]);	
	M[LD    ][i][j]=log_add3(lgep, M[D][i       ][j-step_j]+lgop,M[LD][i       ][j-step_j]);
	*/
	
	}
    }
  retun M;
MatState* RviterbiL_hmm (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL,MatState *S, MatState *E)
{
  MatState *Mid;
  Mid=viterbiL_hmm (A,ns, ls,H, CL, S, E);
  
  if (!Mid)
    {
      return S;
    }
  else if ( Mid->n)
    {
      return Mid;
    }
      
  else
    {
      Mid->p=S;S->n=Mid;
      Mid->n=E;E->p=Mid;
      RviterbiL_hmm (A,ns, ls,H, CL,S,   Mid);
      RviterbiL_hmm (A,ns, ls,H, CL,Mid, E);
      return S;
    }
}

MatState* viterbiL_hmm (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL, MatState *S,MatState *E)
{
  int current,memory, dim;
  double e, v,t;
  int i,j,pi,pj, s, k;
  int start_i, start_j, end_i, end_j, l1, l2;
  HmmState *S1, *S2;
  static MatState ****M;
  static int maxl;
  MatState *Mid=NULL;

  MatState *CC, *PCC;
  int midpoint;
  int Delta;
  

  if ( MatStateAreIdentical ( S, E))return NULL;
  
  l1=strlen (A->seq_al[ls[0][0]]);
  l2=strlen (A->seq_al[ls[1][0]]);
    
  midpoint=S->i+(E->i-S->i)/2;
  Delta=E->i-S->i;

  start_i=S->i;end_i=E->i;
  start_j=S->j;end_j=E->j;
    
  dim=H->order+2;current=0;memory=H->order+1;
  if (!M || (l2+1)>maxl)
    { free_arrayN((void **)M, 4);
      M=declare_arrayN(4, sizeof ( MatState), dim, maxl=(l2+1), H->nS,1);
    }
      
  /*MAKE THE VITERBI FROM S(tart) to E(nd)*/
  for ( i=start_i; i<=end_i; i++)
    {
      M= (MatState****)recycle  ( (void **)M,H->order+1,1);
      for ( j=start_j; j<=end_j; j++)
	{
	  for ( s=H->nS-1;s>=0; s--)
	    {

	      S1=H->S[s];
	      pi=i-S1->DI;pj=j-S1->DJ;
	      CC=M[current][j][s];
	      CC->i=i; CC->j=j; CC->st=s;CC->sc=H->forbiden;CC->p=CC->n=CC->m=NULL;CC->sc=H->forbiden; 
	      if (i==start_i && j==start_j && s==S->st)  {CopyMatState(S,CC);}
	      else if ( i==end_i && j==end_j && s==E->st && s!=H->end)
		{
		  S2=H->S[E->pst];
		  CopyMatState(E,CC);
		  CC->p=M[S1->DI][j-S1->DJ][S2->state];
		}
	      else if ( pi<start_i || pj<start_j)        {CC->sc=H->forbiden;}
	      else 
		{
		  for (k=1; k<=H->fromM[S1->state][0]; k++)
		    {
		      S2=H->S[H->fromM[s][k]];
		      PCC=M[S1->DI][j-S1->DJ][S2->state];
	
		      if      ( pi+pj!=0 && S2->state==H->start) {t=H->forbiden;}
		      else if ( !(pi==l1 && pj==l2) && s==H->end){t=H->forbiden;}
		      else   t=H->T[S2->state][S1->state];

		      v=hmm_add(t,PCC->sc);
		      if ( v!=H->forbiden && (CC->sc==H->forbiden || v> CC->sc)){CC->sc=v; CC->pst=S2->state;CC->p=PCC;}
		    }
		  
		  e=(S1->em==H->forbiden)?S1->em_func (A, A->pos, ns[0], ls[0],i-1, A->pos,ns[1], ls[1], j-1, CL):S1->em;
		  CC->sc=hmm_add(CC->sc,e);
		}
	      
	      if (i==midpoint)CC->m=CopyMatState(CC, M[memory][j][s]);
	      else if (i>midpoint && CC->sc!=H->forbiden) CC->m=(M[S1->DI][j-S1->DJ][CC->pst])->m;
	    }
	}
    }  
  
  if ( E->st==H->end)CopyMatState ((M[current][end_j][E->st]),E);
  
  if ( Delta>1)
    {
      Mid=CopyMatState ((M[current][end_j][E->st])->m,NULL);
    }
  else if ( Delta==1)
    { 
      CC=M[current][E->j][E->st];
      Mid=E;
      while (!MatStateAreIdentical (CC->p, S) )
	{
	  Mid->p=CopyMatState(CC->p,NULL);
	  (Mid->p)->n=Mid;
	  Mid=Mid->p;CC=CC->p;
	}
      Mid->p=S;
      S->n=Mid;
      Mid=S;
    }

  return Mid;
}
#endif
/*********************************COPYRIGHT NOTICE**********************************/
/*© Centre National de la Recherche Scientifique (CNRS) */
/*and */
/*Please Cite: Notredame*/
/*Mon May 17 20:15:35 MDT 2004. */
/*All rights reserved.*/
/*NOTICE:                                                                                                                                     |*/
/*  This file is an integral part of the */
/*  ALIGN_TWO_SEQ Software. */
/*  Its content is protected and all */
/*  the conditions mentioned in the licensing */
/*  agreement of the software apply to this file.*/
/*...............................................                                                                                      |*/
/*  If you need some more information, or if you */
/*  wish to obtain a full license, please contact: */
/*  cedric.notredame@europe.com*/
/*...............................................                                                                                                                                     |*/
/**/
/**/
/*	*/
/*********************************COPYRIGHT NOTICE**********************************/
