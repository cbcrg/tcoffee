#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"


static double mc_delta_matrix ( int ***mat1, int ***mat2, char **alp, int nch);
static double delta_matrix ( int **mat1,int **mat2, char *alp);
static double ***pavie_seq2pavie_fmat (Sequence *S,double *gop, double *gep, char **mat, char *idmat, int id_threshold, int sample_size, int nch, char *param  );
static int **pavie_fmat2pavie_logodd_mat (double **fmat, char *alp);
static double **pavie_aln2fmat(Alignment *A, double **fmat, char *idmat, int id_threshold, int ch, int nch, char *param);
static int pavie_mat2pavie_id_mat ( int **mat,char *in_name, char *alp, char *ignore, char *force,int T, char *out_name);
static double paviemat2gep ( int **mat, char *alp);
static Alignment *align_pavie_sequences (char *seq0,char *seq1,char **mat,double *gop,double *gep,int nch, char *param);
static int pavie_score (char *s0,int p0, char *s1,int p1,char **mat_file, double *gop, double *gep, int nch, float factor, char *param);
static char **seq2pavie_alp (Sequence *S, int nch);
static Sequence * seq2pavie_seq ( Sequence *S, int nch);
static FILE* output_pavie_aln (Alignment *A, int nch, FILE *fp);
static char **output_pavie_mat_list ( int ***current_mat, double *gep, char **alp, int nch,char *prefix,int cycle, char **mat_name);
static float pavie_aln2id ( Alignment *A, int mode);
static int check_pavie_cl ( char *string);
float pavie_aln2delta_age ( Alignment *A,int s0, int s1, int a0, int a1);

static float tgep_factor;
static int id_thres_used_aln;
static int log_odd_mode;
Sequence * pavie_seq2noisy_seq ( Sequence *S, int freq, char *alp)
{
  int a, b, l1, l2;
  
  vsrand(0);
  
  if (alp==NULL)
    {
      char **x;
      x=seq2pavie_alp (S,1);
      alp=x[0];
    }
  
  l2=strlen (alp);
  for (a=0; a< S->nseq; a++)
    {
      l1=strlen (S->seq[a]);
      for ( b=0; b<l1; b++)
	{
	  if ( (rand ()%100+1)<freq)
	    
	    S->seq[a][b]=alp[rand()%l2];
	}
    }
  return S;
}
Sequence * pavie_seq2random_seq ( Sequence *S, char *subst)
{
  int a, b, r, l;

  
  vsrand (0);
  r=subst[0]; subst++;
  l=strlen (subst);
  for ( a=0; a< S->nseq; a++)
    for (b=0; b<S->len[a]; b++)
      if ( S->seq[a][b]==r)S->seq[a][b]=subst[rand()%l];
  return S;
}

double **pavie_seq2pavie_aln(Sequence *S,char *mat, char *mode)
{
  int a, b,c, nch=0;
  char **mat_list;
  char *buf;

  double *gep, *gop;
  Alignment *A;
  char **alp;
  char *pavie_idmat;
  FILE *fp;
  double **dist_mat;
  float score;
  
  check_pavie_cl (mode);
  
  mat_list=declare_char (100, 100);
  
  if ( is_matrix (mat))
    {
      sprintf ( mat_list[nch++], "%s", mat);
    }
  else
    {
      fp=vfopen (mat,"r");
      while ( (c=fgetc(fp))!=EOF)
	{
	  ungetc(c, fp);
	  fscanf (fp, "%s\n",mat_list[nch++]);
	}
      vfclose (fp);
    }
  
  alp=seq2pavie_alp (S, nch);
  S=seq2pavie_seq (S, nch);

  gop=vcalloc (nch, sizeof (double));
  gep=vcalloc (nch, sizeof (double));
 
  for ( a=0; a< nch; a++)
    {
      int **m;
      char *st;
      int v;
      m=read_matrice (mat_list[a]);
      if ((st=vstrstr(mode, "_GEP")))
	{
	  sscanf ( st, "_GEP%d_", &v);
	  gep[a]=v*-1;
	}
      else if ( m[0][GAP_CODE]==0)
      	{
	  gep[a]=paviemat2gep(m,alp[a]);
	}
      else
	{
	  gep[a]=m[0][GAP_CODE];
	}
      free_int (m, -1);
    }
  
  
  if ( (buf=vstrstr (mode, "_TGEPF")))
    {

      sscanf (buf, "_TGEPF%f_", &tgep_factor);
      tgep_factor/=(float)100;
    }
  else
    {
      tgep_factor=0.5;
    }

  pavie_idmat=vtmpnam(NULL);

  pavie_mat2pavie_id_mat (NULL,"idmat", alp[0],"X","",1,pavie_idmat);
  dist_mat=declare_double ( S->nseq, S->nseq);


  
  for ( a=0; a< S->nseq-1; a++)
    {
      for ( b=a+1; b< S->nseq; b++)
        {
	  int a0, a1;
	  float delta_a;

	 
	  if ( ! strstr (mode, "_MSA_"))
	    {
	      A=align_pavie_sequences (S->seq[a],S->seq[b],mat_list,gop,gep,nch, mode);
	      sprintf ( A->name[0], "%s", S->name[a]);
	      sprintf ( A->name[1], "%s", S->name[b]);
	    }
	  else 
	    {

	      A=strings2aln ( 2, S->name[a], S->seq[a], S->name[b], S->seq[b]);
	      sprintf ( A->seq_al[0], "%s", S->seq[a]);
	      sprintf ( A->seq_al[1], "%s", S->seq[b]);
	      A->len_aln=strlen (S->seq[a]);
	      ungap_aln (A);
	    }

	  if (strm (mode, "_ID01_"))A->score=score=pavie_aln2id (A, 1);
	  else if ( vstrstr (mode, "_ID02_"))A->score=score=pavie_aln2id (A, 2);
	  else if ( vstrstr (mode, "_ID04_"))A->score=score=pavie_aln2id (A, 4);
	  else if ( vstrstr (mode, "_ID05_"))A->score=score=pavie_aln2id (A, 5);
	  else if ( vstrstr (mode, "_ID06_"))A->score=score=pavie_aln2id (A, 6);
	  
	  else A->score=score=pavie_aln2id (A, 1);
	  
	  a0=S->seq[a][strlen(S->seq[a])+1];
	  a1=S->seq[b][strlen(S->seq[b])+1];
	  
	  delta_a=pavie_aln2delta_age (A, 0, 1, a0, a1);
	  
	  if ( vstrstr (mode, "_MATDIST_"))    
	    dist_mat[a][b]=dist_mat[b][a]=(double)(vstrstr (mode, "_ID05") || vstrstr (mode, "_ID06"))?-score*100:(100-score);
	  else if ( vstrstr (mode, "_MATSIM_"))
	    dist_mat[a][b]=dist_mat[b][a]=(double)(score);
	  
	  
	  if ( !vstrstr (mode, "_MAT") )
	    {
	      fprintf ( stdout, "#############\nAlignment %s %s: %d %% ID SCORE %d DELTA_AGE %.2f\n", S->name[a], S->name[b], A->score, A->score_aln, delta_a);
	      output_pavie_aln (A,nch, stdout);
	    }
	  free_aln(A);
	}
    }

  if ( vstrstr (mode, "_MAT") && !vstrstr ( mode, "_NOPRINT_"))
    {
      if ( vstrstr (mode, "_MFORMAT2"))
	{
	  int max, n;
	  float *tot,s, bigtot=0;
	  
	  for (max=0, a=0; a< S->nseq; a++)max=MAX(max,(strlen (S->name[a])));
	  tot=vcalloc ( S->nseq, sizeof (float));
	  fprintf (stdout, "# TC_DISTANCE_MATRIX_FORMAT_01\n");
	  for ( a=0; a<S->nseq; a++)
	    fprintf ( stdout, "# SEQ_INDEX %s %d\n",S->name[a],a);
	  fprintf ( stdout, "# PW_SEQ_DISTANCES \n");
	  for (n=0,a=0;a< S->nseq-1; a++)
	    {
	      for ( b=a+1; b<S->nseq; b++, n++)
		{
		  s=dist_mat[a][b];
	   	  
		  fprintf (stdout, "BOT\t %4d %4d\t %5.2f %*s\t %*s\t %5.2f\n", a,b,s,max,S->name[a], max, S->name[b], s);
		  fprintf (stdout, "TOP\t %4d %4d\t %5.2f %*s\t %*s\t %5.2f\n", b,a,s,max,S->name[b], max, S->name[a], s);
		  tot[a]+=s;
		  tot[b]+=s;
		  bigtot+=s;
		}
	    }
	  for ( a=0; a< S->nseq; a++)
	    {
	      fprintf (stdout, "AVG\t %d\t %*s\t %*s\t %5.2f\n", a,max,S->name[a], max, "*", tot[a]/(S->nseq-1));
	    }
	  vfree (tot);
	  fprintf (stdout, "TOT\t %*s\t %*s\t %5.2f\n", max,"TOT", max, "*", bigtot/n);
	  vfclose (stdout);
	}
      else
	{
	  for ( a=0; a<S->nseq; a++)
	    {
	      fprintf ( stdout, "\n%s ", S->name[a]);
	      for ( b=0; b< S->nseq; b++)
		fprintf ( stdout, "%6d ", (int)(dist_mat[a][b]*100));
	    }
	}
    }
  
  return dist_mat;
}

float pavie_aln2delta_age ( Alignment *A,int s0, int s1, int a0, int a1)
{
  int a,r0, r1, g0, g1,  n;
  float delta;
  for (n=0, delta=0, a=0; a< A->len_aln; a++)
    {
      r0=A->seq_al[s0][a];
      r1=A->seq_al[s1][a];
      
      g0=!is_gap(r0);
      g1=!is_gap(r1);
      
      a0+=g0;a1+=g1;
      if ( g0 && g1)
	{
	  delta+=FABS((a0-a1));
	  n++;
	}
    }
  delta/=(float)((n)?n:1);
  return delta;
}

int **pavie_seq2trained_pavie_mat(Sequence *S, char *param)
{
  double ***fmat;
  int ***current_mat;
  int ***previous_mat;
  char **alp;

  char **mat_file;
  double d,delta_min=10;
  double *gep;
  double *gop;
  
  char ignore[100];
  char force [100];
  char pavie_idmat[100];
  int id_threshold;
  int sample_size;
  char *b;
  int a,n=0,nch=1;
  char *buf;
  
  check_pavie_cl (param);
  
  if ( !param)param=vcalloc (1, sizeof (char));

  if ((b=vstrstr(param,"_THR")))sscanf ( b, "_THR%d_", &id_threshold);
  else id_threshold=0;
  
  sample_size=0;
  if ((b=vstrstr(param,"_SAMPLE")))sscanf ( b, "_SAMPLE%d_", &sample_size);
  if ((b=vstrstr(param,"_PARALOGOUS")))
    {
      sscanf ( b, "_PARALOGOUS%d_", &sample_size);
      sample_size*=-1;
    }
  
  if ((b=vstrstr(param,"_CHANNEL")))sscanf ( b, "_CHANNEL%d_", &nch);
  else nch=1;
  
  if ( (buf=vstrstr (param, "_TGEPF")))
    {
      sscanf (buf, "_TGEPF%f_", &tgep_factor);
      tgep_factor/=(float)100;
    }
  else
    {
      tgep_factor=0.5;
    }
  if ( (buf=vstrstr (param, "_PAMLOGODD_")))
    {
      log_odd_mode=1;
    }
  /*Declare Arrays*/
  gep=vcalloc (nch, sizeof (double));
  gop=vcalloc (nch, sizeof (double));
  mat_file=declare_char ( nch, 100);
  current_mat =vcalloc ( nch, sizeof (double**));
  previous_mat=vcalloc ( nch, sizeof (double**));
  

  sprintf (ignore, "X");
  force[0]='\0';
  sprintf ( pavie_idmat, "pavie_idmat");
  
  
  
  alp=seq2pavie_alp (S, nch);
  
  S=seq2pavie_seq (S, nch);
  
  pavie_mat2pavie_id_mat (NULL,"idmat", alp[0],ignore,force,1,pavie_idmat);
  
  for ( a=0; a<nch; a++)sprintf (mat_file[a], "idmat");
  

  fmat=pavie_seq2pavie_fmat ( S,gop,gep,mat_file,pavie_idmat, id_threshold, sample_size, nch, param);
  

  for (a=0; a<nch; a++)
    {
      current_mat[a]=pavie_fmat2pavie_logodd_mat(fmat[a], alp[a]);
      gep[a]=paviemat2gep(current_mat[a], alp[a]);      
    }
  free_arrayN((void*)fmat, 3);
  
  mat_file=output_pavie_mat_list ( current_mat,gep, alp, nch,"", n++, mat_file); 
  
  
  
  fprintf ( stdout, "\n");
  previous_mat=vcalloc ( nch, sizeof (int**));
  while ((d=mc_delta_matrix (previous_mat, current_mat, alp, nch))>delta_min)
    {

      fprintf ( stdout, "\nDelta=%d: ",(int) d);
      for (a=0; a<nch; a++)
	{
	  free_int (previous_mat[a], -1);
	  previous_mat[a]=current_mat[a];
	}
      fprintf ( stdout, "\n");
	
      fmat=pavie_seq2pavie_fmat (S,gop,gep,mat_file, pavie_idmat, id_threshold, sample_size, nch, param);      
      
     
      for (a=0; a< nch; a++)
	{
	  current_mat[a]=pavie_fmat2pavie_logodd_mat(fmat[a], alp[a]);
	  gep[a]=paviemat2gep(current_mat[a], alp[a]);
	}
      
      mat_file=output_pavie_mat_list ( current_mat,gep, alp, nch,"", n, mat_file); 
      free_arrayN((void*)fmat, 3);
      n++;
    }
  
  fprintf ( stdout, "\nDelta=%d Mat: ",(int) d);
  for (a=0; a<nch; a++)
	{
	  fprintf ( stdout, "%s ",mat_file[a]);
	}
  fprintf ( stdout, "\n");
  
  return current_mat[0];
}

double mc_delta_matrix ( int ***mat1, int ***mat2, char **alp, int nch)
{
  int a;
  double delta=0;
  if ( !mat1 || !mat2) return 100000;
  for ( a=0; a< nch; a++)
    delta+=delta_matrix (mat1[a], mat2[a], alp[a]);
  return delta/nch;
}
      
double delta_matrix ( int **mat1,int **mat2, char *alp)
{
  int ns;
  double delta, v;
  int a, b;
  
  if ( mat1==NULL || mat2==NULL) return 100000;
  
  ns=strlen (alp);
  for (delta=0, a=0; a< ns; a++)
    for ( b=0; b< ns; b++)
      {
	v=mat1[alp[a]-'A'][alp[b]-'A']-mat2[alp[a]-'A'][alp[b]-'A'];
	delta+=v*v;
      }
  delta=sqrt(delta);
  
  return delta;
}

double ***pavie_seq2pavie_fmat (Sequence *S,double *gop, double *gep, char **mat, char *idmat, int id_threshold, int sample_size, int nch, char *param  )
{
  int a=0, b, chan;
  double ***fmat=NULL;
  Alignment *A;
  int exclude_id=1;
  static int tot=0;


  id_thres_used_aln=0;
  fmat=vcalloc ( nch, sizeof (double **));
  
  if (sample_size==0)
    {
  
      for (tot=0, a=0; a< S->nseq-exclude_id; a++)
	{

	  output_completion ( stderr,a+1,S->nseq,1, "");
	  
	  for ( b=a+exclude_id; b< S->nseq; b++)
	    {
	      tot++;
	      
	      A=align_pavie_sequences (S->seq[a],S->seq[b],mat,gop,gep,nch,param);
	      
	      for ( chan=0; chan< nch; chan++)
		{
		 
		  fmat[chan]=pavie_aln2fmat (A, fmat[chan], idmat, id_threshold, chan, nch, param);
		}
	      free_aln (A);
	    }
	}
    }
  else 
    {
      int c;
      static int **list;
      
      if ( sample_size>0 && !list)
	{
	  if ( exclude_id==0)sample_size*=3;
	  if (!list)
	    {
	      list=declare_int ((sample_size+1), 2);
	      vsrand(0);
	      tot=0;
	      while (tot<sample_size)
		{
		  a=rand()%(S->nseq);b=rand()%(S->nseq);
		  if ( a!=b)
		    {
		      list[tot][0]=a;list[tot][1]=b;
		      tot++;
		      if ( exclude_id==0)
			{
			  list[tot][0]=a;list[tot][1]=a;
			  tot++;
			  list[tot][0]=b;list[tot][1]=b;
			  tot++;
			}
		    }
		}
	    }
	}
      else if ( sample_size<0 && !list)
	{

	  int **sim;
	  int m;
	  sim=seq2sim_mat (S, "idmat");
	  sample_size*=-1;
	  list=declare_int (S->nseq*S->nseq, 2);
	  
	  m=S->nseq-exclude_id;
	  for (a=0; a<m; a++)
	    for ( b=a+exclude_id; b<S->nseq; b++)
	      {
		if ( sim[a][b]>sample_size)
		  {
		    list[tot][0]=a;
		    list[tot][1]=b;
		    tot++;
		    fprintf ( stderr, "\n%s %s: %d", S->name[a], S->name[b], sim[a][b]);
		    fprintf ( stderr, "\nKeep %s Vs %s : %d%% ID", S->name[a], S->name[b], sim[a][b]);
		  }
	      }
	  free_int(sim, -1);
	}
      
      for (c=0; c<tot; c++)
	{
	  a=list[c][0];b=list[c][1];
	  A=align_pavie_sequences (S->seq[a],S->seq[b],mat,gop,gep,nch, param);
	  for (chan=0; chan< nch; chan++)
	    fmat[chan]=pavie_aln2fmat (A, fmat[chan], idmat, id_threshold,chan, nch, param);
	  
	  free_aln (A);
	  output_completion ( stderr,c,tot,1, "");
	}
    }
  fprintf ( stderr, "\n\tSample_size: %d Used alignments: %d\n", tot, id_thres_used_aln);
  return fmat;
}



int **pavie_fmat2pavie_logodd_mat (double **fmat, char *alp)
{
  int s1, s2,S1, S2;
  double r1, r2;
  int **mat;
  int a, b;
  int ns;
  int logodd=0;
  
  
  ns=strlen (alp);
  mat=declare_int (256, 256);
  
  for ( a=0; a<ns; a++)
    {
      s1=tolower(alp[a]);
      fprintf ( stderr, "\n\tSymbol %c Freq %5.2f %%", s1, ((float)fmat[s1][0]*100)/(float)fmat[0][0]);
    }
  
	
  
  for (a=0; a<ns; a++)
    for (b=0; b<ns; b++)
      {
	s1=tolower(alp[a]);S1=toupper(alp[a]);
	s2=tolower(alp[b]);S2=toupper(alp[b]);

	if ( log_odd_mode==0)
	  {
	    r1=(fmat[s1][s2]+1)/(fmat[s1][s1]+1);
	    r2=(fmat[s2][s1]+1)/(fmat[s2][s2]+1);
	    logodd=(int)(((log(r1)+log(r2))/2)*PAVIE_MAT_FACTOR);
	  }
      	else if ( log_odd_mode==1)
	  {

	    r1=(fmat[s1][s2]+fmat[s2][s1]+1)/(fmat[1][1]+1);
	    r2=((fmat[s1][1]+1)/(fmat[1][1]+1))*((fmat[s2][1]+1)/(fmat[1][1]+1))*2;
	    logodd=(int)(log(r1/r2)*PAVIE_MAT_FACTOR);
	  }
	    
	mat[s1-'A'][s2-'A']=logodd;
	mat[S1-'A'][S2-'A']=logodd;
	mat[S1-'A'][s2-'A']=logodd;
	mat[s1-'A'][S2-'A']=logodd;
	
      }
  return mat;
}
	
double **pavie_aln2fmat(Alignment *A, double **fmat, char *idmat, int id_threshold, int ch, int nch, char *param)
{
  int a;
  int c1, c2;
  int w, id;
  int l,start;
  
  l=(A->len_aln/nch);
  start=l*ch;
  A->len_aln=l;
  A->seq_al[0]+=start;
  A->seq_al[1]+=start;
  


  if ( fmat==NULL)fmat=declare_double(300, 300);
  
  if (  vstrstr (param, "_TWE00_"))w=100;
  else if ( vstrstr (param, "_TWE01_"))w=pavie_aln2id (A, 1);
  else if ( vstrstr (param, "_TWE02_"))w=pavie_aln2id (A, 2);
  else if ( vstrstr (param, "_TWE03_"))w=pavie_aln2id (A, 3);
  else if ( vstrstr (param, "_TWE04_"))w=pavie_aln2id (A, 4);
  else if ( vstrstr (param, "_TWE05_"))w=pavie_aln2id (A, 5);
  else if ( vstrstr (param, "_TWE06_"))w=pavie_aln2id (A, 6);
  
  else w=pavie_aln2id (A, 3);
  
  id=pavie_aln2id(A, 3);

  
  if (id<id_threshold) 
    {
      A->len_aln*=nch;
      A->seq_al[0]-=start;A->seq_al[1]-=start;
      return fmat;
    }
  else
    {
      id_thres_used_aln++;
      for ( a=0; a<A->len_aln; a++)
	{
	  c1=tolower(A->seq_al[0][a]);
	  c2=tolower(A->seq_al[1][a]);
	  fmat[c1][c2]+=w;

	  fmat[c1][0]++;
	  fmat[c1][1]+=w;
	  
	  fmat[c2][0]++;
	  fmat[c1][1]+=w;

	  fmat[0][0]+=2;
	  fmat[1][1]+=2*w;
	}
      A->len_aln*=nch;
      A->seq_al[0]-=start;A->seq_al[1]-=start;
      
      return fmat;
    }
}

int pavie_mat2pavie_id_mat ( int **mat,char *in_name, char *alp, char *ignore, char *force,int T, char *out_name)
{
  int n1, n2, n3;
  int s1, s2, S1, S2;
  int a, b;
  int **idmat;
  
  if      (mat==NULL && in_name==NULL) return 0;
  else if (mat==NULL)
    {
      mat=read_matrice (in_name);
    }
  

  idmat=declare_int ( 256, 256);
  n1=strlen (alp);
  n2=strlen (ignore);
  n3=strlen (force);
  
  for (a=0; a< n1; a++)
    for ( b=0; b<n1; b++)
      {
	s1=tolower(alp[a])-'A';S1=toupper(alp[a])-'A';
	s2=tolower(alp[b])-'A';S2=toupper(alp[b])-'A';
	idmat[s1][s2]=idmat[s1][S2]=idmat[S1][S2]=idmat[S1][s2]=(mat[s1][s2]>=T)?PAVIE_MAT_FACTOR:0;
      }
  for (a=0; a<n3; a++)
    for (b=0; b<n1; b++)
      {
	s1=tolower(force[a])-'A';S1=toupper(force[a])-'A';
	s2=tolower(alp[b])-'A';S2=toupper(alp[b])-'A';
	idmat[s1][s2]=idmat[s1][S2]=idmat[S1][S2]=idmat[S1][s2]=PAVIE_MAT_FACTOR;
      }

  for (a=0; a<n2; a++)
    for (b=0; b<n1; b++)
      {
	s1=tolower(ignore[a])-'A';S1=toupper(ignore[a])-'A';
	s2=tolower(alp[b])-'A';S2=toupper(alp[b])-'A';
	idmat[s1][s2]=idmat[s1][S2]=idmat[S1][S2]=idmat[S1][s2]=0;
      }


  

  output_pavie_mat (idmat, out_name, 0, alp);
  free_int (idmat, -1);
  return 1;
}
double paviemat2gep ( int **mat, char *alp)
{
  int a, b, l;
  int n=0;
  double gep=0;
  l=strlen ( alp);
  
  for (a=0; a<l-1; a++)
    for ( b=a+1; b< l; b++)
      {
	gep+=mat[alp[a]-'A'][alp[b]-'A'];
	n++;
      }
  gep/=n;
 
  return gep;
  
}

Alignment *align_pavie_sequences (char *seq0,char *seq1,char **mat,double *gop,double *gep,int nch, char *param)
{
  double **F; 
  int **T;
  Alignment *A;
  int XL, YL, len;
  int i, j, a, b, c;
  double match, gap_inX, gap_inY, MXY=1, GX=2, GY=3;
  
  char *ax, *ay;
  char *bufx, *bufy, *buf;
  char *x,*y;
  float factor;

  factor=tgep_factor;
  
  
  /*FActor
    terminal gaps are set to gep=gep*factor
  */
  
  if ( strm (seq0, seq1))
    {
      A=declare_aln2 (2,strlen(seq0)+1);
      A->len_aln=strlen (seq0);
      A->nseq=2;
      A->score=A->score_aln=100;
      sprintf ( A->seq_al[0], "%s", seq1);
      sprintf ( A->seq_al[1], "%s", seq0);
      return A;
    }
  

  x=seq0;
  y=seq1;
  
  XL=strlen (x)/nch;
  YL=strlen (y)/nch;
  
  
  ax=vcalloc ( (YL+XL)*nch+1, sizeof (char));
  ay=vcalloc ( (YL+XL)*nch+1, sizeof (char));
  bufx=vcalloc ( (YL+XL)*nch+1, sizeof (char));
  bufy=vcalloc ( (YL+XL)*nch+1, sizeof (char));
  
  F=declare_double (XL+2, YL+2);
  T=declare_int (XL+2, YL+2);
  
  
  /*Fill stage*/
  F[0][0] = 0;
  for(i = 1; i <=XL; i++) 
    {
      
      F[i][0] = F[i-1][0]+pavie_score (x,i-1,NULL,GAP_CODE,mat, gop, gep, nch, factor, param) /*CL->M[x[i-1]-'A'][gap]*/;
      
      T[i][0] = GY;
    }
  
  for(j = 1; j <= YL; j++) 
    {

      F[0][j] = F[0][j-1]+pavie_score (NULL,GAP_CODE,y,j-1,mat, gop, gep, nch, factor, param)/*CL->M[y[j-1]-'A'][gap]*/;
      T[0][j] = GX;
    }

  
  for(i = 1; i <= XL; i++) 
    {
      for(j = 1; j <= YL; j++) 
	{
	  
	  match  = F[i-1][j-1] + /*CL->M[x[i-1]-'A'][y[j-1]-'A']*/pavie_score (x,i-1,y, j-1,mat, gop, gep, nch, 1, param);
	  gap_inY= F[i-1][j] + /*CL->M[x[i-1]-'A'][gap]*/         pavie_score (x,i-1, NULL,GAP_CODE,mat, gop, gep, nch, (j==YL)?factor:1, param); 
	  gap_inX= F[i][j-1] +  /*+ CL->M[y[j-1]-'A'][gap]*/      pavie_score (NULL,GAP_CODE,y, j-1,mat, gop, gep, nch, (i==XL)?factor:1, param);
	  
	  if ( match >= gap_inY && match >=gap_inX){F[i][j]=match; T[i][j]=MXY;}
	  else if ( gap_inX>=gap_inY){F[i][j]=gap_inX; T[i][j]=GX;}
	  else {F[i][j]=gap_inY; T[i][j]=GY;}
	}
    }
  /*Trace back stage*/
  
  
  i = XL; 
  j = YL; 
  len=0;
  while(!(i==0 && j==0)) 
    {
      
      if   (T[i][j]==MXY) 
	{
	  ax[len]=1;i--;
	  ay[len]=1;j--;
	}
      else if ( T[i][j]==GY)
	{
	  ax[len]=1;i--;
	  ay[len]='-';
	}
      else if ( T[i][j]==GX)
	{
	  ax[len]='-';
	  ay[len]=1;j--;
	}
      len++;
    }
  
  for (a=0; a<len; a++)
    for (b=1; b<nch; b++)
      {
	ax[a+len*b]=ax[len];
	ay[a+len*b]=ay[len];
      }
  len=len*nch;
  ax[len]='\0';
  ay[len]='\0';

 
  sprintf ( bufx, "%s", ax);
  sprintf ( bufy, "%s", ay);
  
  for (a=1;a<nch; a++)
    {
      strcat (ax, bufx);
      strcat (ay, bufy);
    }
  
  buf=ax;ax=invert_string (ax);vfree(buf);
  buf=ay;ay=invert_string (ay);vfree(buf);
  
  
  A=declare_aln2 (2,strlen(ax)+1);

  
  A->len_aln=strlen (ax);
  A->nseq=2;
  A->score=A->score_aln=F[XL][YL];
  
  for (a=0, b=0, c=0; a<A->len_aln; a++)
    {
      if (ax[a]==1)ax[a]=seq0[b++];
      if (ay[a]==1)ay[a]=seq1[c++];
    }
  

  
  sprintf ( A->seq_al[0], "%s", ax);
  sprintf ( A->seq_al[1], "%s", ay);
  
  vfree (ax); vfree(ay);vfree (bufx); vfree (bufy);free_double(F, -1); free_int (T, -1);
  return A;
} 


int pavie_score (char *s0,int p0, char *s1,int p1,char **mat_file, double *gop, double *gep, int nch, float factor, char *param)
  
  {
    static char *cmat;
    static int  ***mat;
    static int use_age;
    static int mchscore=-1;
    
    int l0, l1, c0, c1;
    int a, score=0;

    if ( !use_age)
      {
	strget_param ( param, "_AGECHANNEL", "-1", "%d", &use_age);

      }
    if (mchscore==-1)
      {
	strget_param (param, "_MCHSCORE", "0", "%d", &mchscore);
	
      }
    
    if ( !cmat || !mat_file || !strm (cmat, mat_file[0]))
      {
	if ( !cmat)cmat=vcalloc ( 100, sizeof (char));
	sprintf ( cmat, "%s", (mat_file)?mat_file[0]:"idmat");
	if ( !mat)mat=vcalloc ( nch, sizeof (int**));
	for ( a=0; a< nch; a++)
	  {
	    if ( mat[a])free_int (mat[a], -1);
	    mat[a]=read_matrice ((mat_file)?mat_file[a]:"idmat");
	    
	  }
      }
    
    l0=(s0)?strlen (s0)/nch:0;
    l1=(s1)?strlen (s1)/nch:0;

    if (mchscore==0);
    else if (mchscore==1) score=999999;
    else if (mchscore==2)score=-9999999;
    else 
      {
	HERE ("Error: mchscore >2 [FATAL]\n");
	exit (EXIT_FAILURE);
      }
    for ( a=0; a< nch; a++)
      {
	int s;
	c0=(s0)?s0[l0*a+p0]-'A':p0;
	c1=(s1)?s1[l1*a+p1]-'A':p1;
	if ( c0==GAP_CODE)s=(gep[a]!=0)?gep[a]:mat[a][c1][GAP_CODE];
	else if ( c1==GAP_CODE)s=(gep[a]!=0)?gep[a]:mat[a][c0][GAP_CODE];
	else s=mat[a][c0][c1];

	if (mchscore==0)score+=s;
	else if (mchscore==1)score=MIN(s, score);
	else if (mchscore==2)score=MAX(s, score);
	

      }
    
    if ( use_age>0 && s0 && s1)
      {
	
	int a0, a1;
	int s;
	
	a0=s0[strlen(s0)+1];
	a1=s1[strlen(s1)+1];
	
	a0+=p0;
	a1+=p1;
	s=use_age*FABS((a0-a1))*-1;
	
	if (mchscore==0)score+=s;
	else if (mchscore==1)score=MIN(s, score);
	else if (mchscore==2)score=MAX(s, score);
      }
	
	
    score*=factor;
    return score;
  }
Sequence * seq2pavie_seq ( Sequence *S, int nch)
  {
    char *buf, *p;
    int a, b;
    
    S->nseq/=nch;

    for (b=0; b<S->nseq; b++)
      {
	
	buf=vcalloc ((strlen (S->seq[b])*nch)+10, sizeof (char));
	for ( a=0; a< nch; a++)
	  {
	    strcat (buf, S->seq[b+(S->nseq)*a]);
	    vfree ( S->seq[b+(S->nseq)*a]);
	  }
	S->seq[b]=buf;
	/*Code Age on the byte just after the string termination*/
	
	if ((p=strstr (S->seq_comment[b], "FIRSTYEAR")))
	  {
	    sscanf ( p, "FIRSTYEAR%d", (int*)&(S->seq[strlen(buf)+1]));
	  }
	
      }
    return S;
  }
char **seq2pavie_alp (Sequence *S, int nch)
  {
    int a, n;
    char **alp;
    
    n=S->nseq/nch;
    alp=vcalloc (nch, sizeof (char*));
    for ( a=0; a< nch; a++)
      {
	alp[a]=array2alphabet (S->seq+n*a, n, "-.");
      }
    return alp;
  }
FILE *output_pavie_aln (Alignment *A, int nch, FILE *fp)
{
  int a, b, c,d, l, start, end;
  Alignment *B;
  Sequence *S;
  B=declare_aln2(A->nseq*nch+nch, A->len_aln);
  


  l=A->len_aln/nch;

  for ( a=0; a< nch; a++)
    {
      for (b=0; b< A->nseq; b++, B->nseq++)
	{
	  sprintf (B->name[B->nseq], "%s.c%d", A->name[b], a);
	  start=l*a;end=start+l;
	  for (d=0,c=start; c<end; c++, d++)B->seq_al[B->nseq][d]=A->seq_al[b][c];
	  B->seq_al[B->nseq][d]='\0';
	}
      if ( a!=nch-1)
	{
	  B->name[B->nseq][0]='\0';
	  for ( b=0; b<l; b++)B->seq_al[B->nseq][b]='^';
	  B->nseq++;
	}
    }
  
  B->len_aln=l;
  fp=output_Alignment_without_header (B,fp);
  S=free_aln (B);
  free_sequence (S, S->nseq);
  return fp;
  
}
char **output_pavie_mat_list ( int ***current_mat, double *gep, char **alp, int nch,char *prefix,int cycle, char **mat_name)
{
  int a;
  char mat_list_name[100];
  FILE *fp;
  char latest[1000];
  char current[1000];
  char command[1000];
  
  sprintf ( mat_list_name, "pavie_matrix%s.cycle_%d.mat_list", prefix, cycle+1);
  fp=vfopen ( mat_list_name, "w");
  fprintf ( stderr, "\n\tOutput Pavie Matrix: %s", mat_list_name);
  for ( a=0; a< nch; a++)
    {
      sprintf ( mat_name[a], "pavie_matrix%s.ch_%d.cy_%d.pavie_mat", prefix,a+1, cycle+1);
      sprintf  (latest, "pavie_matrix%s.ch_%d.cy_last.pavie_mat",prefix,a+1);
      sprintf  ( current, "matrix.ch%d.pavie_mat", a);
      fprintf ( stderr, "\n\t  Channel %d Matrix: %s",a+1, mat_name[a]);
      output_pavie_mat (current_mat[a],mat_name[a],gep[a], alp[a]);
      sprintf ( command, "cp %s %s", mat_name[a], latest);
      system  (command);
      sprintf ( command, "cp %s %s", latest, current);
      system (command);
      fprintf ( fp, "%s\n", mat_name[a]);
    }
  vfclose (fp);
  return mat_name;
}


int pavie_pair_wise (Alignment *A,int *ns, int **l_s,Constraint_list *CL )
{
  double **F; int **T;
  char *x,*y;
  char *ax, *ay;
  int XL, YL, len;
  int i, j;
  double match, gap_inX, gap_inY, MXY=1, GX=2, GY=3;
  int gap=GAP_CODE;
  char *ix, *iy;
  float factor=0.5;


  /*factor:
    decreases terminal gap penalties with a factor X
    factor=1: terminal gap penalties <=> internal gap penalties
  */



  x=A->seq_al[l_s[0][0]];
  y=A->seq_al[l_s[1][0]];
  XL=strlen (x);
  YL=strlen (y);
  
  ax=vcalloc ( YL+XL+1, sizeof (char));
  ay=vcalloc ( YL+XL+1, sizeof (char));
  
  
  F=declare_double (XL+2, YL+2);
  T=declare_int (XL+2, YL+2);
  
 
  /*Fill stage*/
  F[0][0] = 0;
  for(i = 1; i <=XL; i++) 
    {

      F[i][0] = F[i-1][0]+(CL->M[x[i-1]-'A'][gap]*factor);
      T[i][0] = GY;
    }
  
  for(j = 1; j <= YL; j++) 
    {
      F[0][j] = F[0][j-1]+CL->M[y[j-1]-'A'][gap]*factor;
      T[0][j] = GX;
    }

  
  for(i = 1; i <= XL; i++) 
    {
      for(j = 1; j <= YL; j++) 
	{

	  match  = F[i-1][j-1] + CL->M[x[i-1]-'A'][y[j-1]-'A'];
	  gap_inY= F[i-1][j]   + (CL->M[x[i-1]-'A'][gap]*(j==YL)?factor:1);
	  gap_inX= F[i][j-1]   + (CL->M[y[j-1]-'A'][gap]*(i==XL)?factor:1);
	
	  if ( match >= gap_inY && match >=gap_inX){F[i][j]=match; T[i][j]=MXY;}
	  else if ( gap_inX>=gap_inY){F[i][j]=gap_inX; T[i][j]=GX;}
	  else {F[i][j]=gap_inY; T[i][j]=GY;}
	}
    }
  /*Trace back stage*/
  A->score=A->score_aln=F[XL][YL];
  
  i = XL; 
  j = YL; 
  len=0;
  while(!(i==0 && j==0)) 
    {
      
      if   (T[i][j]==MXY) 
	{
	  ax[len]=x[--i];
	  ay[len]=y[--j];
	}
      else if ( T[i][j]==GY)
	{
	  ax[len]=x[--i];
	  ay[len]='-';
	}
      else if ( T[i][j]==GX)
	{
	  ax[len]='-';
	  ay[len]=y[--j];
	}
      len++;
    }
  ax[len]='\0';
  ay[len]='\0';
  
  ix=invert_string (ax); iy=invert_string(ay);
  A=realloc_aln (A,len+1);
  
  sprintf ( A->seq_al[l_s[0][0]], "%s", ix);
  sprintf ( A->seq_al[l_s[1][0]], "%s", iy);
  A->nseq=2;
  A->len_aln=len;
  
  vfree (ax); vfree(ay);vfree(ix); vfree(iy); free_double(F, -1); free_int (T, -1);
  return A->score;
} 

float pavie_aln2id ( Alignment *A, int mode)
{
  int a, id=0, match=0, l1=0, l2=0, r1, r2, is_res1, is_res2;
  


  for (a=0; a<A->len_aln; a++)
    {
      r1=A->seq_al[0][a];
      r2=A->seq_al[1][a];
      
      
      is_res1=(!is_gap(r1) && r1!='x' && r1!='X')?1:0;
      is_res2=(!is_gap(r2) && r2!='x' && r2!='X')?1:0;

      l1+=is_res1;
      l2+=is_res2;
      
      
      if ( is_res1 && is_res2 )
	{
	  match++;
	  id+=(r1==r2)?1:0;
	}
    }


  if ( mode==1)return (match==0)?0:((id*100)/match);
  else if (mode ==2) return (A->len_aln==0)?0:((id*100)/A->len_aln);
  else if (mode ==3) return (MIN(l1,l2)==0)?0:((id*100)/(MIN(l1,l2)));
  else if (mode ==4) return (MAX(l1,l2)==0)?0:((id*100)/(MAX(l1,l2)));
  else if (mode ==5)return (A->score_aln * -1)/*/PAVIE_MAT_FACTOR*/;
  else if (mode ==6)return ((MAX(l1,l2)==0)?0:((A->score_aln)/(MAX(l1,l2))))*-1;
  else
    {
      fprintf ( stderr, "\nUnknown Mode [pavie_aln2id:FATAL:%s]", PROGRAM);
      myexit (EXIT_FAILURE);
      return EXIT_FAILURE;
    }

}
		   
int check_pavie_cl ( char *string)
{
  if ( !string || string[0]=='\0' ) return 1;
  else if (( string[0]!='_') ||string [strlen (string)-1]!='_')
    {
      fprintf ( stderr, "ERROR: parameters must start and finish with an underscore: _parameters_ [FATAL:%s]\n", PROGRAM);
      myexit (EXIT_FAILURE);
    }
  return 1;
}

Alignment *pavie_seq2pavie_sort ( Sequence *S, char *mat, char *mode)
{
  int a, b, c=0, avg_c;
  int **avg;
  double **dm;
  Alignment *A=NULL;
  char **new_order;
  
  if ( vstrstr (mode, "_IDSORT_") || vstrstr (mode, "_MASTERSORT"))
    {
      char *buf;
      buf=vcat ( mode, "_MATDIST_NOPRINT_");
      dm=pavie_seq2pavie_aln (S, mat,buf);
      avg=declare_int (S->nseq, 2);
      if ( vstrstr (mode,"_IDSORT_"))
	{
	  
	  for ( a=0; a< S->nseq; a++)
	    {
	      avg[a][0]=a;
	      for ( b=0; b<S->nseq; b++)
		if ( b!=a)avg[a][1]+=(int) dm[a][b];
	      avg[a][1]/=S->nseq-1;
	    }
	  sort_int ( avg, 2, 1, 0, S->nseq-1);
	  
	  c=avg[0][0];
	  avg_c=avg[0][1];
	  fprintf ( stderr, "\nAVG\t %s\t %s\t %d",S->name[c],"avg", avg_c);
	}
      else if ( vstrstr (mode, "_MASTERSORT"))
	{
	  char name[100];
	  char *s;
	  s=vstrstr(mode, "_MASTERSORT");
	  mode=substitute ( mode, "_", " ");
	  sscanf (s, " MASTERSORT%s", name);
	  
	  mode=substitute (mode, " ", "_");
	  c=name_is_in_list ( name, S->name, S->nseq, 100);
	  
	  if ( c==-1)
	    {
	      fprintf ( stderr, "\nERROR: Sequence %s is not in the dataset [FATAL:%s]", name, PROGRAM);
	      myexit (EXIT_FAILURE);
	    }
	}
      
      for ( a=0; a<S->nseq; a++)
	{
	  avg[a][0]=a;
	  if ( a!=c)avg[a][1]=dm[c][a];
	  else avg[a][1]=-1;
	}
      
      sort_int ( avg, 2, 1, 0, S->nseq-1);

      new_order=declare_char ( S->nseq, 100);
      sprintf (new_order[0], "%s", S->name[c]);
      for ( a=1; a<S->nseq; a++) 
	{
	  sprintf ( new_order[a], "%s", S->name[avg[a][0]]);
	  
	  fprintf ( stderr, "\nTOP\t %s\t %s\t %d", S->name[c],new_order[a] , avg[a][1]);
	}
      
      fprintf ( stderr, "\n");
      A=seq2aln (S, NULL,RM_GAP);
      A=reorder_aln (A, new_order, A->nseq);
      vfree ( buf);
      free_double(dm, -1);free_int (avg, -1);free_char (new_order, -1);
    }
  else if ( vstrstr ( mode, "_TREESORT_"))
    {
     A=pavie_seq2pavie_msa (S, mat, mode);
    }
  else
    {
      fprintf ( stderr, "\nERROR: pavie_seq2sort <matrix> <_IDSORT_ | _TREESORT_>");
    }
  return A;
  
}
NT_node pavie_seq2pavie_tree (Sequence *S, char *mat, char *mode)
{
  double **dm;
  char *tree_name,*buf;
  
    
  buf=vcat (mode,"_MATDIST_NOPRINT_");  
  dm=pavie_seq2pavie_aln (S, mat,buf);
  dist2nj_tree (dm,S->name, S->nseq,tree_name=vtmpnam (NULL));
  
  free_double(dm, -1);vfree (buf);
  
  return main_read_tree (tree_name);
}

Alignment* pavie_seq2pavie_msa ( Sequence *S, char *mat_in, char *mode)
{
  Constraint_list *CL;
  char **alp, *s;
  Alignment *A;
  NT_node **FT, T;
  int a;
  char mat[100];
  
  
  A=seq2aln (S, NULL, RM_GAP);
  CL=declare_constraint_list (S, NULL, NULL, 0, NULL, NULL);
  sprintf ( CL->dp_mode,   "myers_miller_pair_wise");
  sprintf ( CL->tree_mode, "nj");
  sprintf ( CL->distance_matrix_mode, "idscore");
  CL=choose_extension_mode ("matrix", CL);
  
  if ( !is_matrix (mat_in))
    {
      FILE *fp;
      fp=vfopen ( mat_in, "r");
      fscanf (fp, "%s", mat);
      vfclose (fp);
      add_warning( stderr, "\nWarning: Multiple Channel Not Supported. Used First Channel Only for MSA [Matrix: %s][WARNING:%s]", mat, PROGRAM);
    }
  else
    {
      sprintf ( mat, "%s", mat_in);
    }
  
  CL->M=read_matrice (mat);
  CL->gop=0;
  
  alp=seq2pavie_alp (S, 1);
  CL->gep=paviemat2gep(CL->M, alp[0]);
  
  
  CL->pw_parameters_set=1;
  CL->local_stderr=stderr;

  if ( vstrstr (mode, "_QUICKTREE_"))
    {
      FT=make_tree (A, CL, CL->gop, CL->gep,S, NULL,MAXIMISE);
      T=FT[3][0];
    }
  else if ( (s=vstrstr (mode, "_USETREE")))
    {
      char fname[100];
      mode=substitute ( mode, "_", " ");
      sscanf (s, " USETREE%s", fname);
      mode=substitute (mode, " ", "_");
      T=main_read_tree (fname);
    }
  else
    {
      T=pavie_seq2pavie_tree ( S, mat_in, mode);
    }
  
  for ( a=0; a< A->nseq; a++)ungap (A->seq_al[a]);
  
  tree_aln (T->left,T->right,A,(CL->S)->nseq, CL);
  A=reorder_aln ( A,A->tree_order,A->nseq);
  
  return A;
}
