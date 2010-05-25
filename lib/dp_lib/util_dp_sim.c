#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"
/* extern char name0[], name1[]; */
/* extern int match, mismh; */



static Constraint_list *CL;
static int * ns;
static int **l_s;
static Alignment *Aln;
static int **pos;
static int *seqc0, *seqc1;
static int min0,min1,max0,max1,mins;

static void* sim_vcalloc( size_t nobj, size_t size);
static void sim_free_all (); 
static int sim_reset_static_variable ();
static int big_pass(int *A,int *B,int M,int N,int K, int nseq) ;
static int locate(int *A,int *B,int nseq); 
static int small_pass(int *A,int *B,int count,int nseq);
static int no_cross ();
static int diff_sim( int *A,int *B,int M,int N,int tb,int te); 
int calcons(int *aa0,int n0,int *aa1,int n1,int *res,int *nc,int *nident, Alignment *A, int *ns, int **l_s, Constraint_list *CL);



#define SIM_GAP -1
#define min(x,y) ((x)<=(y) ? (x) : (y))
//#define TC_SCORE_SIM(x,y) TC_SCORE (x,y)

static int q, r;			/* gap penalties */
static int qr;				/* qr = q + r */



typedef struct ONE { int COL ;  struct ONE  *NEXT ;} pair, *pairptr;
pairptr *row, z, z1; 			/* for saving used aligned pairs */


#define PAIRNULL (pairptr)NULL
static int tt;

typedef struct SIM_NODE
	{ int  SIM_SCORE;
	  int  SIM_STARI;
	  int  SIM_STARJ;
	  int  SIM_ENDI;
	  int  SIM_ENDJ;
	  int  SIM_TOP;
	  int  SIM_BOT;
	  int  SIM_LEFT;
	  int  SIM_RIGHT; }  vertex,
#ifdef FAR_PTR
 far *vertexptr;
#else
     *vertexptr;
#endif
		
vertexptr  *LIST;			/* an array for saving k best scores */
vertexptr  low = 0;			/* lowest score node in LIST */
vertexptr  most = 0;			/* latestly accessed node in LIST */
static int numnode;			/* the number of nodes in LIST */

static int *CC, *DD;			/* saving matrix scores */
static int *RR, *SS, *EE, *FF; 		/* saving start-points */
static int *HH, *WW;		 	/* saving matrix scores */
static int *II, *JJ, *XX, *YY; 		/* saving start-points */
static int  m1, mm, n1, nn;		/* boundaries of recomputed area */
static int  rl, cl;			/* left and top boundaries */
static int  lmin;			/* minimum score in LIST */
static int flag;			/* indicate if recomputation necessary*/

/* DIAG() assigns value to x if (ii,jj) is never used before */
#define DIAG(ii, jj, x, value)				\
{ for ( tt = 1, z = row[(ii)]; z != PAIRNULL; z = z->NEXT )	\
    if ( z->COL == (jj) )				\
      { tt = 0; break; }				\
  if ( tt )						\
    x = ( value );					\
}

/* replace (ss1, xx1, yy1) by (ss2, xx2, yy2) if the latter is large */
#define ORDER(ss1, xx1, yy1, ss2, xx2, yy2)		\
{ if ( ss1 < ss2 )					\
    { ss1 = ss2; xx1 = xx2; yy1 = yy2; }		\
  else							\
    if ( ss1 == ss2 )					\
      { if ( xx1 < xx2 )				\
	  { xx1 = xx2; yy1 = yy2; }			\
	else						\
	  if ( xx1 == xx2 && yy1 < yy2 )		\
	    yy1 = yy2;					\
      }							\
}

/* The following definitions are for function diff() */

static int  zero = 0;				/* int type zero        */
#define gap(k)  ((k) <= 0 ? 0 : q+r*(k))	/* k-symbol indel score */

static int *sapp;				/* Current script append ptr */
static int  last;				/* Last script op appended */

static int I, J;				/* current positions of A ,B */
static int no_mat; 				/* number of matches */ 
static int no_mis; 				/* number of mismatches */ 
static int al_len; 				/* length of alignment */
						/* Append "Delete k" op */
#define DEL(k) \
{ I += k;\
  al_len += k;\
  if (last < 0)\
    last = sapp[-1] -= (k);\
  else\
    last = *sapp++ = -(k);\
}
						/* Append "Insert k" op */
#define INS(k) \
{ J += k;\
  al_len += k;\
  if (last < 0)\
    { sapp[-1] = (k); *sapp++ = last; }	\
  else\
    last = *sapp++ = (k);\
}

						/* Append "Replace" op */
#define REP \
{ last = *sapp++ = 0;\
  al_len += 1;\
}


/*
int sim_pair_wise_lalign (Alignment *in_A, int *in_ns, int **in_l_s,Constraint_list *in_CL)
{
  if ( in_ns[0]==1 && in_ns[1]==1)
    return sim_pair_wise_lalign (in_A, in_ns, in_l_s,in_CL);
  else
  */
  
    



int sim_pair_wise_lalign (Alignment *in_A, int *in_ns, int **in_l_s,Constraint_list *in_CL)
/* SIM(A,B,M,N,K,V,Q,R) reports K best non-intersecting alignments of
   the segments of A and B in order of similarity scores, where
   V[a][b] is the score of aligning a and b, and -(Q+R*i) is the score
   of an i-symbol indel. 
*/ 						
{
  int endi, endj, stari, starj;	/* endpoint and startpoint */ 
  int  score;   			/* the max score in LIST */
  int count;				/* maximum size of list */	
  int i;
  int  *S;				/* saving operations for diff */
  int nc, nident;		/* for display */
  vertexptr cur; 			/* temporary pointer */
  vertexptr findmax();	 		/* return the largest score node */
  double percent;
  int t1, t2, g1, g2, r1, r2;
  int a, b, c, d, e;
/*cedric was here 11/2/99*/
  int CEDRIC_MAX_N_ALN=999;
  int CEDRIC_THRESHOLD=50;
  int *A, *B;
  int M, N, K, maxl;
  int nseq;
  int R, Q;
  Alignment *DA;
  

  DA=in_A;

  Aln=copy_aln (in_A, NULL);

  

  l_s=in_l_s;
  ns=in_ns;
  CL=in_CL;
  K=CL->lalign_n_top;
 
  M=strlen (Aln->seq_al[l_s[0][0]]);
  N=strlen (Aln->seq_al[l_s[1][0]]);
  maxl=M+N+1;

  pos=aln2pos_simple (Aln,-1, ns, l_s);
  
  seqc0=(int*)sim_vcalloc (maxl,sizeof (int));
  A=(int*)sim_vcalloc (maxl,sizeof (int));
  for ( a=0; a<maxl; a++){seqc0[a]=A[a]=a;}
  A[M+1]='\0';

  seqc1=(int*)sim_vcalloc (maxl,sizeof (int));
  B=(int*)sim_vcalloc (maxl,sizeof (int));
  for ( a=0; a<maxl; a++){seqc1[a]=B[a]=a;}
  B[N+1]='\0';
  
  nseq=(l_s[0][0]!=l_s[1][0])?2:1;  
  
  
  Q=MAX(CL->gop, -CL->gop)*SCORE_K;
  R=MAX(CL->gep, -CL->gep)*SCORE_K;
  
 
  
  if ( K==CEDRIC_MAX_N_ALN)K--;
  else if ( K<0)
      {
       
       CEDRIC_THRESHOLD=-K; 
       K=CEDRIC_MAX_N_ALN;
      }
  
  /* allocate space for all vectors */
  
  CC = ( int * ) sim_vcalloc(N+1, sizeof(int));
  DD = ( int * ) sim_vcalloc(N+1, sizeof(int));
  RR = ( int * ) sim_vcalloc(N+1, sizeof(int));
  SS = ( int * ) sim_vcalloc(N+1, sizeof(int));
  EE = ( int * ) sim_vcalloc(N+1, sizeof(int));
  FF = ( int * ) sim_vcalloc(N+1, sizeof(int));
  
  HH = ( int * ) sim_vcalloc(M + 1, sizeof(int));
  WW = ( int * ) sim_vcalloc(M + 1, sizeof(int));
  II = ( int * ) sim_vcalloc(M + 1, sizeof(int));
  JJ = ( int * ) sim_vcalloc(M + 1, sizeof(int));
  XX = ( int * ) sim_vcalloc(M + 1, sizeof(int));
  YY = ( int * ) sim_vcalloc(M + 1, sizeof(int));
  S = ( int * )  sim_vcalloc(min(M,N)*5/4+1, sizeof (int));
  row = ( pairptr * ) sim_vcalloc( (M + 1), sizeof(pairptr));


  /* set up list for each row */
  if (nseq == 2) for ( i = 1; i <= M; i++ ) row[i]= PAIRNULL;
  else {
	  z = ( pairptr )sim_vcalloc (M,(int)sizeof(pair));
	  for ( i = 1; i <= M; i++,z++) {
		  row[i] = z;
		  z->COL = i;			
		  z->NEXT = PAIRNULL;
	  }
  }

  
  q = Q;
  r = R;
  qr = q + r;

  LIST = ( vertexptr * ) sim_vcalloc( K, sizeof(vertexptr));
  for ( i = 0; i < K ; i++ )
    LIST[i] = ( vertexptr )sim_vcalloc( 1, sizeof(vertex));
  

  numnode = lmin = 0;
  big_pass(A,B,M,N,K,nseq);
  
 
  
  /* Report the K best alignments one by one. After each alignment is
     output, recompute part of the matrix. First determine the size
     of the area to be recomputed, then do the recomputation         */
  

  for ( count = K - 1; count >= 0; count-- )
    { if ( numnode == 0 )
        {
	  
	  padd_aln (in_A);
	  /*fatal("The number of alignments computed is too large");*/
	  sim_free_all();
	  return 1;
	}
	
      cur = findmax();	/* Return a pointer to a node with max score*/
      score = cur->SIM_SCORE;
      if ( K==CEDRIC_MAX_N_ALN && score<CEDRIC_THRESHOLD)break;
      stari = ++cur->SIM_STARI;
      starj = ++cur->SIM_STARJ;
      endi = cur->SIM_ENDI;
      endj = cur->SIM_ENDJ;
      m1 = cur->SIM_TOP;
      mm = cur->SIM_BOT;
      n1 = cur->SIM_LEFT;
      nn = cur->SIM_RIGHT;
      rl = endi - stari + 1;
      cl = endj - starj + 1;
      I = stari - 1;
      J = starj - 1;
      sapp = S;
      last = 0;
      al_len = 0;
      no_mat = 0;
      no_mis = 0;
      diff_sim(&A[stari]-1, &B[starj]-1,rl,cl,q,q);


      min0 = stari;
      min1 = starj;
      max0 = stari+rl-1;
      max1 = starj+cl-1;
      calcons(A+1,M,B+1,N,S,&nc,&nident, Aln,ns, l_s, CL);
      percent = (double)nident*100.0/(double)nc;
      
      
      
      /*Min0: index of the last residue before the first in a 1..N+1 numerotation*/



      if (!DA->A)DA->A=copy_aln(Aln, DA->A);
      DA->A=realloc_alignment (DA->A,nc+1);
      
 
      DA=DA->A;
      DA->A=NULL;

      for ( c=0; c< 2; c++)
	    {
	    for ( a=0; a< ns[c]; a++) 
		{
		  e=(c==0)?min0:min1;
		  for ( d=0; d<e; d++)
		    {
		      DA->order[l_s[c][a]][1]+=1-is_gap(Aln->seq_al[l_s[c][a]][d]);
		    }
		} 
	    }
      
      
      for ( t1=min0,t2=min1,a=0; a<nc; a++)
	{
	  r1=seqc0[a];
	  r2=seqc1[a];
	  
	  g1=(r1==SIM_GAP || r1>M)?0:1;
	  g2=(r2==SIM_GAP || r2>N)?0:1;
	  t1+=g1;
	  t2+=g2;
	  for (b=0; b<ns[0]; b++)DA->seq_al[l_s[0][b]][a]=(g1)?Aln->seq_al[l_s[0][b]][A[t1-1]]:'-';
	  for (b=0; b<ns[1]; b++)DA->seq_al[l_s[1][b]][a]=(g2)?Aln->seq_al[l_s[1][b]][B[t2-1]]:'-';
	}
      for (b=0; b<ns[0]; b++){DA->seq_al[l_s[0][b]][a]='\0';}
      for (b=0; b<ns[1]; b++){DA->seq_al[l_s[1][b]][a]='\0';}
     
      DA->nseq=ns[0]+ns[1];
      DA->len_aln=nc;
      DA->score=percent;
      DA->score_aln=score;
      fflush(stdout);

      
      if ( count )
	{ flag = 0;
	  locate(A,B,nseq);
	  if ( flag )
	    small_pass(A,B,count,nseq);
	}
    }
  padd_aln (in_A);
  
  sim_free_all();
  free_int (pos, -1);
  free_aln (Aln);

  
 
  return 1;
  
}
int sim_reset_static_variable ()
{
  CC=DD=RR=SS=EE=FF=HH=WW=II=JJ=XX=YY=sapp=NULL;
  min0=min1=max0=max1=mins=q=r=qr=tt=numnode=m1=n1=nn=rl=cl=lmin=flag=zero=last=I=J=no_mat=no_mis=al_len=0;
  most=low=NULL;/*Very important: cause a bug if not reset*/
  LIST=NULL;    /*Very important: cause a bug if not reset*/
  return 0;
}
/* A big pass to compute K best classes */


int big_pass(int *A,int *B,int M,int N,int K, int nseq) 
{ register  int  i, j;			/* row and column indices */
  register  int  c;			/* best score at current point */
  register  int  f;			/* best score ending with insertion */
  register  int  d;			/* best score ending with deletion */
  register  int  p;			/* best score at (i-1, j-1) */
  register  int  ci, cj;		/* end-point associated with c */ 
  register  int  di, dj;		/* end-point associated with d */
  register  int  fi, fj;		/* end-point associated with f */
  register  int  pi, pj;		/* end-point associated with p */
  
  int   addnode();			/* function for inserting a node */

	
	/* Compute the matrix and save the top K best scores in LIST
	   CC : the scores of the current row
	   RR and EE : the starting point that leads to score CC
	   DD : the scores of the current row, ending with deletion
	   SS and FF : the starting point that leads to score DD        */
 	/* Initialize the 0 th row */
	for ( j = 1; j <= N ; j++ )
	  {  CC[j] = 0;
	     RR[j] = 0;
	     EE[j] = j;
	     DD[j] = - (q);
	     SS[j] = 0;
	     FF[j] = j;
	  }
	for ( i = 1; i <= M; i++) 
	  {  c = 0;				/* Initialize column 0 */
	     f = - (q);
	     ci = fi = i;
	     if ( nseq == 2 )
	       { p = 0;
	         pi = i - 1;
	         cj = fj = pj = 0;
	       }
	     else
	       { p = CC[i];
		 pi = RR[i];
		 pj = EE[i];
	         cj = fj = i;
	       }
	     for ( j = (nseq == 2 ? 1 : (i+1)) ; j <= N ; j++ )  
	       {  f = f - r;
		  c = c - qr;
		  ORDER(f, fi, fj, c, ci, cj)
		  c = CC[j] - qr; 
		  ci = RR[j];
		  cj = EE[j];
		  d = DD[j] - r;
		  di = SS[j];
		  dj = FF[j];
		  ORDER(d, di, dj, c, ci, cj)
		  c = 0;
		 
		  DIAG(i, j, c, p+TC_SCORE(A[i-1],B[j-1]))		/* diagonal */
		    
		  if ( c <= 0 )
		    { c = 0; ci = i; cj = j; }
		  else
		    { ci = pi; cj = pj; }
		  ORDER(c, ci, cj, d, di, dj)
		  ORDER(c, ci, cj, f, fi, fj)
		  p = CC[j];
		  CC[j] = c;
		  pi = RR[j];
		  pj = EE[j];
		  RR[j] = ci;
		  EE[j] = cj;
		  DD[j] = d;
		  SS[j] = di;
		  FF[j] = dj;
		  if ( c > lmin )	/* add the score into list */
		    lmin = addnode(c, ci, cj, i, j, K, lmin);
	        }
	  }
return 1;
}

/* Determine the left and top boundaries of the recomputed area */

int locate(int *A,int *B,int nseq) 
{ register  int  i, j;			/* row and column indices */
  register  int  c;			/* best score at current point */
  register  int  f;			/* best score ending with insertion */
  register  int  d;			/* best score ending with deletion */
  register  int  p;			/* best score at (i-1, j-1) */
  register  int  ci, cj;		/* end-point associated with c */ 
  register  int  di=0, dj=0;		/* end-point associated with d */
  register  int  fi, fj;		/* end-point associated with f */
  register  int  pi, pj;		/* end-point associated with p */
  int  cflag, rflag;			/* for recomputation */
  int   addnode();			/* function for inserting a node */
  int  limit;				/* the bound on j */

	/* Reverse pass
	   rows
	   CC : the scores on the current row
	   RR and EE : the endpoints that lead to CC
	   DD : the deletion scores 
	   SS and FF : the endpoints that lead to DD

	   columns
	   HH : the scores on the current columns
	   II and JJ : the endpoints that lead to HH
	   WW : the deletion scores
	   XX and YY : the endpoints that lead to WW
	*/
	for ( j = nn; j >= n1 ; j-- )
          {  CC[j] = 0;
	     EE[j] = j;
	     DD[j] = - (q);
	     FF[j] = j;
	     if ( nseq == 2 || j > mm )
                RR[j] = SS[j] = mm + 1;
	     else
                RR[j] = SS[j] = j;
	  }

        for ( i = mm; i >= m1; i-- )
	  {  c = p = 0;
	     f = - (q);
	     ci = fi = i;
	     pi = i + 1;
	     cj = fj = pj = nn + 1;
	     
	     if ( nseq == 2 || n1 > i )
		limit = n1;
	     else
		limit = i + 1;
	     for ( j = nn; j >= limit ; j-- )  
	       {  f = f - r;
		  c = c - qr;
		  ORDER(f, fi, fj, c, ci, cj)
		  c = CC[j] - qr; 
		  ci = RR[j];
		  cj = EE[j];
		  d = DD[j] - r;
		  di = SS[j];
		  dj = FF[j];
		  ORDER(d, di, dj, c, ci, cj)
		  c = 0;
		  DIAG(i, j, c, p+TC_SCORE(A[i-1],B[j-1]))		/* diagonal */
		  
		  if ( c <= 0 )
		    { c = 0; ci = i; cj = j; }
		  else
		    { ci = pi; cj = pj; }
		  ORDER(c, ci, cj, d, di, dj)
		  ORDER(c, ci, cj, f, fi, fj)
		  p = CC[j];
		  CC[j] = c;
		  pi = RR[j];
		  pj = EE[j];
		  RR[j] = ci;
		  EE[j] = cj;
		  DD[j] = d;
		  SS[j] = di;
		  FF[j] = dj;
		  if ( c > lmin )
		    flag = 1;
	        }
	     if ( nseq == 2 || i < n1 )
	       { HH[i] = CC[n1];
	         II[i] = RR[n1];
	         JJ[i] = EE[n1];
	         WW[i] = DD[n1];
	         XX[i] = SS[n1];
	         YY[i] = FF[n1];
	       }
	  }
      
  for ( rl = m1, cl = n1; ; )
    { for ( rflag = cflag = 1; ( rflag && m1 > 1 ) || ( cflag && n1 > 1 ) ;  )
        { if ( rflag && m1 > 1 )	/* Compute one row */
            { rflag = 0;
	      m1--;
      	      c = p = 0;
	      f = - (q);
	      ci = fi = m1;
	      pi = m1 + 1;
	      cj = fj = pj = nn + 1;

	      for ( j = nn; j >= n1 ; j-- )  
	        { f = f - r;
		  c = c - qr;
		  ORDER(f, fi, fj, c, ci, cj)
		  c = CC[j] - qr; 
		  ci = RR[j];
		  cj = EE[j];
		  d = DD[j] - r;
		  di = SS[j];
		  dj = FF[j];
		  ORDER(d, di, dj, c, ci, cj)
		  c = 0;
		  DIAG(m1, j, c, TC_SCORE(A[m1-1],B[j-1]))		/* diagonal */
		 
		  if ( c <= 0 )
		    { c = 0; ci = m1; cj = j; }
		  else
		    { ci = pi; cj = pj; }
		  ORDER(c, ci, cj, d, di, dj)
		  ORDER(c, ci, cj, f, fi, fj)
		  p = CC[j];
		  CC[j] = c;
		  pi = RR[j];
		  pj = EE[j];
		  RR[j] = ci;
		  EE[j] = cj;
		  DD[j] = d;
		  SS[j] = di;
		  FF[j] = dj;
		  if ( c > lmin )
		     flag = 1;
		  if ( ! rflag && ( (ci > rl && cj > cl) || (di > rl && dj > cl)
				    || (fi > rl && fj > cl) ) )
		      rflag = 1;
	        }
	      HH[m1] = CC[n1];
	      II[m1] = RR[n1];
	      JJ[m1] = EE[n1];
	      WW[m1] = DD[n1];
	      XX[m1] = SS[n1];
	      YY[m1] = FF[n1];
	      if ( ! cflag && ( (ci > rl && cj > cl) || (di > rl && dj > cl)
				|| (fi > rl && fj > cl )) )
	         cflag = 1;
	    }

	  if ( nseq == 1 && n1 == (m1 + 1) && ! rflag )
	     cflag = 0;
	  if ( cflag && n1 > 1 )	/* Compute one column */
	    { cflag = 0;
	      n1--;
	      c = 0;
	      f = - (q);
	      cj = fj = n1;
	      if ( nseq == 2 || mm < n1 )
		{ p = 0;
	          ci = fi = pi = mm + 1;
	          pj = n1 + 1;
		  limit = mm;
		}
	      else
		{ p = HH[n1];
		  pi = II[n1];
		  pj = JJ[n1];
	          ci = fi = n1;
		  limit = n1 - 1;
		}
	      for ( i = limit; i >= m1 ; i-- )  
	        { f = f - r;
		  c = c - qr;
		  ORDER(f, fi, fj, c, ci, cj)
		  c = HH[i] - qr; 
		  ci = II[i];
		  cj = JJ[i];
		  d = WW[i] - r;
		  di = XX[i];
		  dj = YY[i];
		  ORDER(d, di, dj, c, ci, cj)
		  c = 0;
	          DIAG(i, n1, c, p+TC_SCORE(A[i-1], B[n1-1]))
		 
		   
		  
		  if ( c <= 0 )
		    { c = 0; ci = i; cj = n1; }
		  else
		    { ci = pi; cj = pj; }
		  ORDER(c, ci, cj, d, di, dj)
		  ORDER(c, ci, cj, f, fi, fj)
		  p = HH[i];
		  HH[i] = c;
		  pi = II[i];
		  pj = JJ[i];
		  II[i] = ci;
		  JJ[i] = cj;
		  WW[i] = d;
		  XX[i] = di;
		  YY[i] = dj;
		  if ( c > lmin )
		     flag = 1;
	          if ( ! cflag && ( (ci > rl && cj > cl) || (di > rl && dj > cl)
				    || (fi > rl && fj > cl )) )
		     cflag = 1;
	        }
	      CC[n1] = HH[m1];
	      RR[n1] = II[m1];
	      EE[n1] = JJ[m1];
	      DD[n1] = WW[m1];
	      SS[n1] = XX[m1];
	      FF[n1] = YY[m1];
	      if ( ! rflag && ( (ci > rl && cj > cl) || (di > rl && dj > cl)
				|| (fi > rl && fj > cl) ) )
	         rflag = 1;
	    }
	}
      if ( (m1 == 1 && n1 == 1) || no_cross() )
	 break;
   }
  m1--;
  n1--;
return 1;
}

/* recompute the area on forward pass */
int small_pass(int *A,int *B,int count,int nseq)
{ register  int  i, j;			/* row and column indices */
  register  int  c;			/* best score at current point */
  register  int  f;			/* best score ending with insertion */
  register  int  d;			/* best score ending with deletion */
  register  int  p;			/* best score at (i-1, j-1) */
  register  int  ci, cj;		/* end-point associated with c */ 
  register  int  di, dj;		/* end-point associated with d */
  register  int  fi, fj;		/* end-point associated with f */
  register  int  pi, pj;		/* end-point associated with p */
  int   addnode();			/* function for inserting a node */
  int  limit;				/* lower bound on j */

	for ( j = n1 + 1; j <= nn ; j++ )
	  {  CC[j] = 0;
	     RR[j] = m1;
	     EE[j] = j;
	     DD[j] = - (q);
	     SS[j] = m1;
	     FF[j] = j;
	  }
	for ( i = m1 + 1; i <= mm; i++) 
	  {  c = 0;				/* Initialize column 0 */
	     f = - (q);
	     ci = fi = i;
	     
	     if ( nseq == 2 || i <= n1 )
	       { p = 0;
	         pi = i - 1;
	         cj = fj = pj = n1;
		 limit = n1 + 1;
	       }
	     else
	       { p = CC[i];
		 pi = RR[i];
		 pj = EE[i];
	         cj = fj = i;
		 limit = i + 1;
	       }
	     for ( j = limit ; j <= nn ; j++ )  
	       {  f = f - r;
		  c = c - qr;
		  ORDER(f, fi, fj, c, ci, cj)
		  c = CC[j] - qr; 
		  ci = RR[j];
		  cj = EE[j];
		  d = DD[j] - r;
		  di = SS[j];
		  dj = FF[j];
		  ORDER(d, di, dj, c, ci, cj)
		  c = 0;
		  DIAG(i, j, c, p+TC_SCORE(A[i-1], B[j-1]))		/* diagonal */
		    //checked

		  if ( c <= 0 )
		    { c = 0; ci = i; cj = j; }
		  else
		    { ci = pi; cj = pj; }
		  ORDER(c, ci, cj, d, di, dj)
		  ORDER(c, ci, cj, f, fi, fj)
		  p = CC[j];
		  CC[j] = c;
		  pi = RR[j];
		  pj = EE[j];
		  RR[j] = ci;
		  EE[j] = cj;
		  DD[j] = d;
		  SS[j] = di;
		  FF[j] = dj;
		  if ( c > lmin )	/* add the score into list */
		    lmin = addnode(c, ci, cj, i, j, count, lmin);
	        }
	  }
return 1;
}

/* Add a new node into list.  */

int addnode(c, ci, cj, i, j, K, cost)  int c, ci, cj, i, j, K, cost;
{ int found;				/* 1 if the node is in LIST */
  register int d;

  found = 0;
  if ( most != 0 && most->SIM_STARI == ci && most->SIM_STARJ == cj )
    found = 1;
  else
     for ( d = 0; d < numnode ; d++ )
	{ most = LIST[d];
	  if ( most->SIM_STARI == ci && most->SIM_STARJ == cj )
	    { found = 1;
	      break;
	    }
        }
  if ( found )
    { if ( most->SIM_SCORE < c )
        { most->SIM_SCORE = c;
          most->SIM_ENDI = i;
          most->SIM_ENDJ = j;
        }
      if ( most->SIM_TOP > i ) most->SIM_TOP = i;
      if ( most->SIM_BOT < i ) most->SIM_BOT = i;
      if ( most->SIM_LEFT > j ) most->SIM_LEFT = j;
      if ( most->SIM_RIGHT < j ) most->SIM_RIGHT = j;
    }
  else
    { if ( numnode == K )	/* list full */
	 most = low;
      else
         most = LIST[numnode++];
      most->SIM_SCORE = c;
      most->SIM_STARI = ci;
      most->SIM_STARJ = cj;
      most->SIM_ENDI = i;
      most->SIM_ENDJ = j;
      most->SIM_TOP = most->SIM_BOT = i;
      most->SIM_LEFT = most->SIM_RIGHT = j;
    }
  if ( numnode == K )
    { if ( low == most || ! low ) 
        { for ( low = LIST[0], d = 1; d < numnode ; d++ )
            if ( LIST[d]->SIM_SCORE < low->SIM_SCORE )
              low = LIST[d];
	}
      return ( low->SIM_SCORE ) ;
    }
  else
    return cost;
}

/* Find and remove the largest score in list */

vertexptr findmax()
{ vertexptr  cur;
  register int i, j;

  for ( j = 0, i = 1; i < numnode ; i++ )
    if ( LIST[i]->SIM_SCORE > LIST[j]->SIM_SCORE )
       j = i;
  cur = LIST[j];
  if ( j != --numnode )
    { LIST[j] = LIST[numnode];
      LIST[numnode] =  cur;
    }
  most = LIST[0];
  if ( low == cur ) low = LIST[0];
  return ( cur );
}

/* return 1 if no node in LIST share vertices with the area */

int no_cross()
{ vertexptr  cur;
  register int i;

      for ( i = 0; i < numnode; i++ )
	{ cur = LIST[i];
	  if ( cur->SIM_STARI <= mm && cur->SIM_STARJ <= nn && cur->SIM_BOT >= m1-1 && 
	       cur->SIM_RIGHT >= n1-1 && ( cur->SIM_STARI < rl || cur->SIM_STARJ < cl ))
	     { if ( cur->SIM_STARI < rl ) rl = cur->SIM_STARI;
	       if ( cur->SIM_STARJ < cl ) cl = cur->SIM_STARJ;
	       flag = 1;
	       break;
	     }
	}
      if ( i == numnode )
	return 1;
      else
	return 0;
}

/* diff(A,B,M,N,tb,te) returns the score of an optimum conversion between
   A[1..M] and B[1..N] that begins(ends) with a delete if tb(te) is zero
   and appends such a conversion to the current script.                   */

int diff_sim( int *A,int *B,int M,int N,int tb,int te) 

{ int   midi, midj, type;	/* Midpoint, type, and cost */
  int midc;

  {
    register int   i, j;
    register int c, e, d, s;
    int t;
    
    
    /* Boundary cases: M <= 1 or N == 0 */
    
    if (N <= 0)
      { if (M > 0) DEL(M)
	  return - gap(M);
      }
    if (M <= 1)
      { if (M <= 0)
	  { INS(N);
	    return - gap(N);
	  }
	if (tb > te) tb = te;
	midc = - (tb + r + gap(N) );
	midj = 0;
	
	for (j = 1; j <= N; j++)
	  {  for ( tt = 1, z = row[I+1]; z != PAIRNULL; z = z->NEXT )	
              if ( z->COL == j+J )			
		{ tt = 0; break; }		
	    if ( tt )			
	      { c = TC_SCORE (A[0],B[j-1]) - ( gap(j-1) + gap(N-j) );
		//checked
		
		if (c > midc)
		  { midc = c;
		    midj = j;
		  }
	      }
	  }
	if (midj == 0)
	  { INS(N) DEL(1) }
	else
	  { if (midj > 1) INS(midj-1)
	      REP
	      if ( A[1] == B[midj] )
		no_mat += 1;
	      else
		no_mis += 1;
	    /* mark (A[I],B[J]) as used: put J into list row[I] */	
	    I++; J++;
	    
	    
	    z = ( pairptr )sim_vcalloc(1,sizeof(pair));
	    z->COL = J;			
	    z->NEXT = row[I];				
	    row[I] = z;
	    if (midj < N) INS(N-midj)
	      }
	return midc;
      }
    
    /* Divide: Find optimum midpoint (midi,midj) of cost midc */
    
    midi = M/2;			/* Forward phase:                          */
    CC[0] = 0;			/*   Compute C(M/2,k) & D(M/2,k) for all k */
    t = -q;
    for (j = 1; j <= N; j++)
      { CC[j] = t = t-r;
	DD[j] = t-q;
      }
    t = -tb;
    for (i = 1; i <= midi; i++)
      { s = CC[0];
	CC[0] = c = t = t-r;
	e = t-q;
	
	for (j = 1; j <= N; j++)
	  { if ((c = c - qr) > (e = e - r)) e = c;
	    if ((c = CC[j] - qr) > (d = DD[j] - r)) d = c;
	    DIAG(i+I, j+J, c, s+TC_SCORE(A[i-1], B[j-1]))
	      //checked
	      
	      if (c < d) c = d;
	    if (c < e) c = e;
	    s = CC[j];
	    CC[j] = c;
	    DD[j] = d;
	  }
      }
    DD[0] = CC[0];
    
    RR[N] = 0;			/* Reverse phase:                          */
    t = -q;			/*   Compute R(M/2,k) & S(M/2,k) for all k */
    for (j = N-1; j >= 0; j--)
      { RR[j] = t = t-r;
	SS[j] = t-q;
      }
    t = -te;
    for (i = M-1; i >= midi; i--)
      { s = RR[N];
	RR[N] = c = t = t-r;
	e = t-q;
	
	for (j = N-1; j >= 0; j--)
	  { if ((c = c - qr) > (e = e - r)) e = c;
	    if ((c = RR[j] - qr) > (d = SS[j] - r)) d = c;
	    DIAG(i+1+I, j+1+J, c, s+TC_SCORE (A[i],B[j])) /*not -1 on purpose*/
	     
	      if (c < d) c = d;
	    if (c < e) c = e;
	    s = RR[j];
	    RR[j] = c;
	    SS[j] = d;
	  }
      }
    SS[N] = RR[N];
    
    midc = CC[0]+RR[0];		/* Find optimal midpoint */
    midj = 0;
    type = 1;
    for (j = 0; j <= N; j++)
      if ((c = CC[j] + RR[j]) >= midc)
	if (c > midc || (CC[j] != DD[j] && RR[j] == SS[j]))
	  { midc = c;
	    midj = j;
	  }
    for (j = N; j >= 0; j--)
      if ((c = DD[j] + SS[j] + q) > midc)
	{ midc = c;
	  midj = j;
	  type = 2;
	}
  }
  
  /* Conquer: recursively around midpoint */
  
  if (type == 1)
    { diff_sim(A,B,midi,midj,tb,q);
      diff_sim(A+midi,B+midj,M-midi,N-midj,q,te);
    }
  else
    { diff_sim(A,B,midi-1,midj,tb,zero);
      DEL(2);
      diff_sim(A+midi+1,B+midj,M-midi-1,N-midj,zero,te);
    }
  return midc;
}

	


int calcons(int *aa0,int n0,int *aa1,int n1,int *res,int *nc,int *nident, Alignment *A, int *ns, int **l_s, Constraint_list *CL)
{
  int i0, i1;
  int op, nid, lenc, nd;
  int *sp0, *sp1;
  int *rp;
  int a, b, id_col, tot_col, r0, r1;

  min0--; min1--;

  sp0 = seqc0+mins;
  sp1 = seqc1+mins;
  rp = res;
  lenc = nid = op = 0;
  i0 = min0;
  i1 = min1;
  
  while (i0 < max0 || i1 < max1) {
    if (op == 0 && *rp == 0) {
      op = *rp++;
      *sp0 = aa0[i0++];
      *sp1 = aa1[i1++];


      for (id_col=tot_col=0,a=0; a< ns[0]; a++)
	for ( b=0; b< ns[1]; b++)
	  {
	    r0=Aln->seq_al[l_s[0][a]][*sp0-1];
	    r1=Aln->seq_al[l_s[1][a]][*sp1-1];
	    
	    if ( !is_gap(r0) && r1==r0)id_col++;
	    if ( !is_gap(r0) && !is_gap(r1))tot_col++;
	  }
      nid+=(tot_col)?(id_col/tot_col):0;
      lenc++;
      sp0++; sp1++;
    }
    else {
      if (op==0) op = *rp++;
      if (op>0) {
	*sp0++ = SIM_GAP;
	*sp1++ = aa1[i1++];
	op--;
	lenc++;
      }
      else {
	*sp0++ = aa0[i0++];
	*sp1++ = SIM_GAP;
	op++;
	lenc++;
      }
    }
  }

  *nident = nid;
  *nc = lenc;

  nd = 0;
  return mins+lenc+nd;
}

/*Memory management */
struct Mem
    {
    void   *p;
    struct Mem *next;
    };

typedef struct Mem Mem;

Mem *first_mem;
Mem *last_mem;

void *sim_vcalloc ( size_t nobj, size_t size)
{
  void *p;
  Mem *new_mem;

  p=vcalloc (nobj, size);
  
  
  new_mem=vcalloc (1, sizeof (Mem));
  if ( last_mem==NULL)first_mem=last_mem=new_mem;
  else
    {
      last_mem->next=new_mem;
      last_mem=new_mem;
    }
  last_mem->p=p;
  return p;
}

void sim_free_all()
{
  Mem *p1, *p2;
  p1=first_mem;


  while (p1)
    {
      p2=p1->next;
      vfree(p1->p);
      vfree(p1);
      p1=p2;
    }
  first_mem=last_mem=NULL;
  sim_reset_static_variable();
}
  
