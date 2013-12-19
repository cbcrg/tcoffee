#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <string.h>
#define LOG_ZERO -999999999
#define LOG_UNDERFLOW 0.000000001
#define DELTA 0.000001
#define MINDELTA 0.001
typedef struct
{
  double **M;
  double **NM;
  double **CM;
  int    **PST;
  int      ns;
  int      ne;
  int     *S;
  int      L;
  double **F;
  double **B;
  double **V;
  int    **PTR;
  int     *VT;
  int     *PT;
  double  *Po;
  double   P;
  double  VP;
  double  PP;
  int     nit;
  int     nrounds;
  char    *out;
}Hmm;


#define UNPACK_HMM(name)\
double **M =name->M;\
double **NM =name->NM;\
double **CM =name->CM;\
int    **PST=name->PST;\
int     ns =name->ns;\
int     ne =name->ne;\
int     *S =name->S;\
int      L =name->L;\
double **F =name->F;\
double **B =name->B;\
double **V =name->V;\
int  **PTR =name->PTR;\
int    *VT =name->VT;\
int    *PT =name->PT;\
double *Po =name->Po;\
int     nit=name->nit;\
int nrounds=name->nrounds;\
char   *out=name->out;
double analyze_stability (int **T1,int nr,int L, int n);

//Hmm
double baum_welch_Ntrain (Hmm *H);
double baum_welch_train  (Hmm *H, int n);
double  baum_welch(Hmm *H);
double  viterbi   (Hmm *H);
double  forward   (Hmm *H);
double  backward  (Hmm *H);
double  posterior (Hmm *H);
double  decode    (Hmm *H);
double  evaluationForward (Hmm *H);
//IO

double **undump_model (char *file, Hmm*H);
double **undump_cmodel (char *file, Hmm*H);
int     *undump_seq(char *file, Hmm *L);

int dump_model    (Hmm *H);
int dump_viterbi  (Hmm *H);
int dump_posterior(Hmm *H);
int dump_evalForward (Hmm *H);


//MATH
double  log_multiply4 (double x, double y, double z, double w);
double  log_multiply3 (double x, double y, double z);
double  log_multiply2 (double x, double y);
double  log_divide(double x, double y);
double  log_add (double x, double y);
double  mylog(double x);
double  myexp(double y);

//Conversion
double **model2modelR(Hmm*H);
int model2pruned_model (Hmm*H);
double **model2modelL(Hmm*H);

//Declare
double **declare_double (int x, int y);
int    **declare_int    (int x, int y);
Hmm     *declare_hmm   (Hmm*H);
main (int argc, char *argv[])
{
  
  Hmm * H;

  //argv[0]: mode
  //argv[1]: sequence
  //argv[2]: model
  //argv[3]: nit
  //argv[4]: out
  //

  srand(time(NULL));
  H=calloc (1, sizeof (Hmm));
  
  H->S  =undump_seq (argv[2],H);
  H->M  =undump_model (argv[3],H);
  model2pruned_model(H);
  H->out=argv[4];
  declare_hmm(H);
  
  if (strcmp (argv[1], "decode")==0)
    {
      decode (H);
    }
  else if ( strcmp(argv[1], "bw")==0)
    {
      undump_cmodel (argv[5],H);
      H->nrounds=atoi(argv[6]);
      H->nit=atoi(argv[7]);
      baum_welch_Ntrain (H);
    }
  else if (strcmp(argv[1], "ev")==0)
  	{
	  undump_cmodel (argv[5],H);
	  evaluationForward (H);
  	}
}

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//                                                                         //
//                                Analyze                                  //
//                                                                         //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
double analyze_stability (int **T,int nr,int L, int n)
{
  int a, b,c,d, pa, pb;
  double sc,tot;

  sc=tot=0;
  for (a=0; a<nr-1; a++)
    {
      for (b=a+1; b<nr; b++)
	{
	  for (c=0; c<n; c++)
	    {
	      pa=rand()%L+1;
	      for (d=0; d<n; d++)
		{
		  pb=rand()%L+1;
		  
		  sc+=((T[a][pa]==T[a][pb])==(T[b][pa]==T[b][pb]))?1:0;
		  tot++;
		}
	    }
	}
    }
  return sc/tot;
}

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//                                                                         //
//                                Evaluation                               //
//                                                                         //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
double evaluationForward(Hmm*H)
{
  double P;

  UNPACK_HMM(H);
  model2modelL (H);

  if (isnan(H->P=forward (H)))return LOG_ZERO;
//  fprintf (stderr, "---- SUMMARY: Best: %10.3f STABILITY: %.3f\n", bP,analyze_stability (stab,nrounds, L, 100));
  fprintf (stderr, "******** Result forward: log(P)= %2.3f\n", H->P);
  dump_evalForward (H);
  return P;

}

/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//                                                                         //
//                                Trainning                                //
//                                                                         //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

double baum_welch_Ntrain (Hmm *H)
{
  int **stab;
  int a, b;
  double cP,bP;

  UNPACK_HMM(H);
  stab=declare_int (nrounds,L+1);
  for (a=0; a<nrounds; a++)
    {
      model2modelR(H);
      cP=baum_welch_train(H, a+1);
      viterbi  (H);
      posterior(H);
      for (b=1; b<L; b++)stab[a][b]=H->VT[b];
      if (a==0 || cP>bP)
	{
	  dump_model    (H);
	  dump_viterbi  (H);
	  dump_posterior(H);
	  bP=cP;
	}
    }
  fprintf (stderr, "---- SUMMARY: Best: %10.3f STABILITY: %.3f\n", bP,analyze_stability (stab,nrounds, L, 100));
  H->P=bP;
  return H->P;
}
double baum_welch_train  (Hmm *H, int round)
{
  UNPACK_HMM(H);
  int idle, it, cont,print;
  int maxidle=10;
  float delta;
  double PP, RP=0;
  for (PP=0,it=0,idle=0, cont=1,print=0; it<nit && idle<maxidle && cont==1; it++)
    {
      double sscore;
      
      H->P=baum_welch (H);
      if (H->P==LOG_ZERO)
	{
	  fprintf ( stderr, "\n\n**** OVERFLOW***\n\n");
	  cont=0;
	}
      
      if (PP)
	{
	  delta=H->P-PP;
	  if (delta<MINDELTA)idle++;
	  else idle=0;
	}
      
      if (!RP)RP=H->P;
      fprintf ( stderr, "\r\t--Round: %3d Start %10.3lf IT: %4d Best: %10.3lf ", round, RP, it+1, H->P);
      PP=H->P;
    }
  fprintf ( stderr, "\n");
  return H->P;
}
double baum_welch(Hmm*H)
{
  
  int k,b, l, i;
  int pseudocount=0;
  
  UNPACK_HMM(H);
  
  
  if (isnan(H->P=backward(H)))return LOG_ZERO;
  if (isnan(H->P=forward (H)))return LOG_ZERO;
    
  for(k=0; k<ns; k++)
    for (l=0; l<(ns+ne); l++)
      NM[k][l]=LOG_ZERO;
  

  for (k=0; k<ns; k++)
    {
      int ll=0;
      while ((l=PST[k][ll++])!=-1)
	{
	  for (i=1; i<L; i++)
	    {
	      double fo, ba, tr, em;
	      int symbol=S[i+1];
	      if (M[l][symbol]==LOG_ZERO){continue;}
	      fo=F[i][k];
	      ba=B[i+1][l];
	      tr=M[k][l];
	      em=M[l][symbol];
	      NM[k][l]=log_add(NM[k][l], log_multiply4(fo,tr,em,ba));
	    }
	  NM[k][l]=log_divide (NM[k][l],H->P);
	}
    }


  //update emissions
  for (k=0; k<ns; k++)
    {
      for (b=ns; b<(ns+ne);b++)
	{
	  for (i=1; i<=L; i++)
	    {
	      if (S[i]==b)
		{
		  NM[k][b]=log_add(NM[k][b], log_multiply2(F[i][k],B[i][k]));
		}
	    }
	  NM[k][b]=log_divide(NM[k][b],H->P);
	}
    }

  //modelL2count+psudocounts
  for (k=0; k<ns; k++)
    for (l=0; l<(ns+ne); l++)
      NM[k][l]=myexp(NM[k][l])+pseudocount;
  

  
  
  //update transitions
  for (k=0; k<ns; k++)
    {
      double t=0,f=0;
      for (l=0; l<ns; l++)
	{
	  if (CM[k][l]==LOG_ZERO)t+=NM[k][l];
	  else f+=CM[k][l];
	}
      t/=(f==1)?1:(1-f);
      
      for (l=0; l<ns; l++)
	{
	  if (CM[k][l]==LOG_ZERO)NM[k][l]/=(t<DELTA)?1:t;
	  else NM[k][l]=CM[k][l];
	}
    }

  //update emissions
  for (k=0; k<ns; k++)
    {
      double t=0,f=0;
      for (l=ns; l<(ns+ne); l++)
	{
	  if (CM[k][l]==LOG_ZERO)t+=NM[k][l];
	  else f+=CM[k][l];
	}
      t/=(f==1)?1:(1-f);
      
      for (l=ns; l<(ns+ne); l++)
	{
	  if (CM[k][l]==LOG_ZERO)NM[k][l]/=(t<DELTA)?1:t;
	  else NM[k][l]=CM[k][l];
	}
    }
  
  //moldel2modelL
  for (k=0; k<ns; k++)
    {
      for (l=0; l<(ns+ne); l++)
	{
	  M[k][l]=mylog(NM[k][l]);
	}
    }
  return H->P;
}
  
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//                                                                         //
//                                Decoding                                 //
//                                                                         //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////	  
double decode (Hmm *H)
{
  posterior (H);
  viterbi   (H);
  
  dump_model(H);
  dump_viterbi(H);
  dump_posterior(H);
  return H->PP;
}
double posterior(Hmm *H)
{
  int i, k, bpost_k, symbol;
  double bpost_score=0,p;
  UNPACK_HMM(H);
   
  H->PP=H->P=backward (H);
  H->PP=H->P=forward  (H);
  for (i=1; i<L; i++)
    {
      symbol=S[i];
      bpost_score=0;
      for (k=0; k<ns; k++)
	{
	  p=log_divide (log_multiply2(F[i][k],B[i][k]),H->PP);
	  if (!bpost_score || p>bpost_score)
	    {
	      bpost_score=p;
	      bpost_k=k;
	    }
	}
      PT[i]=bpost_k;
      Po[i]=bpost_score;
    }
  return H->PP;
}
double viterbi  (Hmm *H)
{
  double max_k,v;
  int ptr_k,symbol, i, l, k;
 
  UNPACK_HMM(H);
  for (l=0; l<ns; l++)V[0][l]=0;


for (i=1; i<=L; i++)
    {
      symbol=S[i];
      for (l=0; l<ns; l++)
	{
	  if (M[l][symbol]==LOG_ZERO){V[i][l]=PTR[i][l]=LOG_ZERO;}
	  else
	    {
	      int kk=0;
	      max_k=LOG_ZERO;
	      ptr_k=-1;
	      while ((k=PST[l][kk++])!=-1)
		{
		  v=log_multiply2(V[i-1][k],M[k][l]);
		  if ( v>max_k || max_k==LOG_ZERO)
		    {
		      max_k=v;
		      ptr_k=k;
		    }
		  V[i][l]=log_multiply2(M[l][symbol], max_k);
		  PTR[i][l]=ptr_k;
		}
	    }
	}
    }


  max_k=LOG_ZERO;
  ptr_k=-1;
  for (k=0; k<ns; k++)
    {
      v=V[L][k];
      if (v>max_k || max_k==LOG_ZERO)
	{
	  max_k=v;
	  ptr_k=k;
	}
    }
  for (i=L; i>=1; i--)
    {
      VT[i]=ptr_k;
      ptr_k=PTR[i][ptr_k];
    }
  H->VP=max_k;
  return H->VP;
}
double forward (Hmm *H)
{
  int k, l, i;
  double emit;
  double P;
  
  UNPACK_HMM(H);
  
  for (k=0; k<ns; k++)F[0][k]=0;


 for (i=1; i<=L; i++)
    for (l=0; l<ns; l++)
      {
	F[i][l]=LOG_ZERO;
	emit=M[l][S[i]];
	for (k=0;k<ns; k++)
	  {
	    F[i][l]=log_add(F[i][l],log_multiply2(F[i-1][k],M[k][l]));
	  }
	F[i][l]=log_multiply2(F[i][l],emit);
      } 

  P=LOG_ZERO;
  for (k=0; k<ns; k++)
    {
      P=log_add (P, F[L][k]);
    }
  
  return P;
}

double backward (Hmm *H)
{
  int k, l, i;
  double emission,P;
  
  UNPACK_HMM(H);
  
  for (k=0; k<ns; k++)
    {
      B[L][k]=LOG_ZERO;
      for (l=0; l<ns; l++)
	{
	  B[L][k]=log_add(B[L][k], M[k][l]);
	}
    }
  
  
  for (i=L-1; i>=1; i--)
    {
      for (k=0; k<ns; k++)
	{
	  int ll=0;
	  B[i][k]=LOG_ZERO;
	  while ((l=PST[k][ll++])!=-1)
	    {
	      emission=M[l][S[i+1]];
	      if (emission==LOG_ZERO)continue;
	      B[i][k]=log_add(B[i][k], log_multiply3(M[k][l],emission,B[i+1][l]));
	    }
	}
    }


  P=LOG_ZERO;
  
  for (l=0; l<ns; l++)
    {
      emission=M[l][S[1]];
      P=log_add (P, log_multiply2(B[1][l],emission));
    }
  return P;
}
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//                                                                         //
//                                Conversion                               //
//                                                                         //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

int model2pruned_model (Hmm*H)
{
  int k,l,t=0;
  UNPACK_HMM(H);
  PST=H->PST=declare_int (ns, ns+1);
  for (l=0; l<ns; l++)
    {
      int n=0;
      for (k=0; k<ns; k++)
	{
	  if (M[l][k]!=LOG_ZERO){PST[l][n++]=k;t++;}
	}
      PST[l][n]=-1;
    }
  return t;
}

//model to model random
double **model2modelR(Hmm*H)
{
  
  double ts, te,ft,fe;
  int a, b;
  UNPACK_HMM(H);
   
  for (a=0; a<ns; a++)
    {
      ft=fe=ts=te=0;
      
      for (b=0; b<(ns+ne); b++)
	{
	  M[a][b]=rand()%10000;
	  if (CM[a][b]!=LOG_ZERO)
	    {
	      if (b<ns)ft+=CM[a][b];
	      else fe+=CM[a][b];
	    }
	  else
	    {
	      if (b<ns)ts+=M[a][b];
	      else te+=M[a][b];
	    }
	}
      ts/=(ft==1)?1:(1-ft);
      te/=(fe==1)?1:(1-fe);
      
      for (b=0; b<(ns+ne); b++)
	{
	  if (CM[a][b]!=LOG_ZERO)M[a][b]=CM[a][b]; 
	  else if (b<ns)M[a][b]/=(ts==0)?1:ts;
	  else M[a][b]/=(te==0)?1:te;
	  
	  M[a][b]=mylog (M[a][b]);
	}
    }
  return M;
}

//moldel2modelL
double **model2modelL(Hmm*H)
{
	int k,l;
	UNPACK_HMM(H);

	for (k=0; k<ns; k++)
		{
	      for (l=0; l<(ns+ne); l++)
	      	  {
	    	  	  M[k][l]=mylog(M[k][l]);
	      	  }
	    }

	return M;

}



/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//                                                                         //
//                                IO                                       //
//                                                                         //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

int *undump_seq(char *file, Hmm *H)
{
  FILE *fp;
  int *S;
  int a=1;
  fp=fopen (file, "r");

  
  fscanf (fp, "%d", &(H->L));
  
  H->S=malloc (sizeof (int)*(H->L+1));
  while (fscanf (fp, "%d", &(H->S)[a++])==1);
  
  fclose (fp);
  return H->S;
}

double **undump_cmodel (char *file, Hmm *H)
{
  
  FILE *fp;
  int a, b;
  float p;
  int ns, ne;
  
  H->CM=declare_double (H->ns, H->ns+H->ne);
  if (strcmp(file, "no")==0)
    {
      
      for (a=0; a<H->ns; a++)
	for (b=0; b<(H->ns+H->ne); b++)
	  H->CM[a][b]=LOG_ZERO;
    }
  else
    {
      fp=fopen (file, "r");
      fscanf (fp, "%f %d %d ",&p, &a, &b);
      for (a=0; a<H->ns; a++)
	{
	  for (b=0;b<(H->ns+H->ne); b++)
	    {
	      fscanf (fp, "%lf ", &H->CM[a][b]);
	    }
	}
      fclose (fp);
    }
  return H->CM;
}
double **undump_model (char *file, Hmm *H)
{
  
  FILE *fp;
  int a, b;
  float p;
  int ns, ne;
  
  fp=fopen (file, "r");
  
  fscanf (fp, "%f %d %d ",&p, &H->ns,&H->ne);
  ns=H->ns; ne=H->ne;
  
  H->M=declare_double(ns,(ne+ns));
  for (a=0; a<ns; a++)
    {
      for (b=0;b<(ns+ne); b++)
	{
	  fscanf (fp, "%lf ", &H->M[a][b]);
	}
    }
  fclose (fp);
  return H->M;
}

int display_model (Hmm *H)//for debug only
{
  FILE *fp;
  UNPACK_HMM(H);
  int a, b;
  
  for (a=0; a<ns; a++)
    {
      for (b=0; b<ns; b++)
	{
	  fprintf (stdout, "ST::%d -- ST::%d --> %.5f\n", a, b, myexp(M[a][b]));
	}
      for (b=ns; b<(ns+ne); b++)
	{
	  fprintf (stdout, "ST::%d -- EM::%d --> %.5f\n", a, b-ns, myexp(M[a][b]));
	}
    }
  return 1;
}
int dump_model (Hmm*H)
{
  FILE *fp;
  UNPACK_HMM(H);
  char file[1000];
  
  sprintf (file, "%s.model",H->out);
  fp=fopen (file, "w");
  int a, b;
  
  fprintf (fp, "%lf %d %d ",H->P, ns, ne);
  for (a=0; a<ns; a++)
    for (b=0; b<(ns+ne); b++)
      fprintf ( fp, "%.5f ",M[a][b]);
  return fclose (fp);
}
int dump_posterior (Hmm *H)
{
  FILE *fp;
  UNPACK_HMM(H);
  char file[1000];
  int a;
  
  sprintf (file, "%s.posterior",H->out);
  fp=fopen (file, "w");
  
  fprintf (fp, "%d %.3lf",H->L, H->PP);
  for (a=1; a<=H->L; a++)
    fprintf ( fp, " %d %.5f",H->PT[a],H->Po[a]);
  return fclose (fp);
}
int dump_viterbi    (Hmm *H)
{
  FILE *fp;
  UNPACK_HMM(H);
  char file[1000];
  int a;
  
  sprintf (file, "%s.viterbi",H->out);
  fp=fopen (file, "w");
  
  fprintf (fp, "%d %.3lf ",H->L, H->VP);
  for (a=1; a<=H->L; a++)
    fprintf ( fp, "%d ",H->VT[a]);
  return fclose (fp);
}

int dump_evalForward (Hmm *H)
{
  FILE *fp;
  UNPACK_HMM(H);
  char file[1000];

  sprintf (file, "%s.eval",H->out);
  fp=fopen (file, "w");
//  int a, b;

  fprintf (fp, "%8.5f", H->P);

  return fclose (fp);
}
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//                                                                         //
//                                MATH                                     //
//                                                                         //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
double log_add (double x, double y)
{
  if (x==LOG_ZERO)return y;
  else if (y==LOG_ZERO)return x;
  else if (x>=y)return x+mylog(1+exp(y-x));
  else return y+mylog (1+myexp(x-y));
}
double log_divide(double x, double y)
{
  if (x==LOG_ZERO || y==LOG_ZERO)return LOG_ZERO;
  return x-y;
}
double log_multiply2 (double x, double y)
{
  if (x==LOG_ZERO || y==LOG_ZERO)return LOG_ZERO;
  return x+y;
}

double log_multiply3 (double x, double y, double z)
{
  if (x==LOG_ZERO || y==LOG_ZERO || z==LOG_ZERO)return LOG_ZERO;
  return x+y+z;
}

double log_multiply4 (double x, double y, double z, double w)
{
  if (x==LOG_ZERO || y==LOG_ZERO || z==LOG_ZERO || w==LOG_ZERO)return LOG_ZERO;
  return x+y+z+w;
}
double mylog (double x)
{
  if (x<LOG_UNDERFLOW)return LOG_ZERO;
  else return log (x);
}
double myexp (double x)
{
  double z;
  if ( x==LOG_ZERO)return 0;
  return exp(x);
  }
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//                                                                         //
//                                DECLARE                                  //
//                                                                         //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
double **declare_double (int x, int y)
{
  double **d;
  int a;

  d=malloc (x*sizeof (double*));
  for (a=0; a<x; a++)d[a]=malloc (y*sizeof (double));
  return d;
}    
int **declare_int (int x, int y)
{
  int **d;
  int a;

  d=malloc (x*sizeof (int*));
  for (a=0; a<x; a++)d[a]=malloc (y*sizeof (int));
  return d;
}    
Hmm* declare_hmm (Hmm *H)
{

  H->NM =declare_double (H->ns, H->ns+H->ne);
  H->F  =declare_double (H->L+1, H->ns);
  H->B  =declare_double (H->L+1, H->ns);
  H->V  =declare_double (H->L+1, H->ns);
  
  H->PTR=declare_int    (H->L+1, H->ns);
  
  H->VT =malloc (((H->L)+1)*sizeof (int));
  H->PT =malloc (((H->L)+1)*sizeof (int));
  H->Po =malloc (((H->L)+1)*sizeof (double));
  return H;
}
