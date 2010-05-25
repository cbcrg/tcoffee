#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>


#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"


/*******************************************************************************/
/*                myers and Miller                                             */
/*                                                                             */
/*	makes DP between the the ns[0] sequences and the ns[1]                 */
/*                                                                             */
/*	for MODE, see the function  get_dp_cost                                */
/*******************************************************************************/
	
#define gap(k)  ((k) <= 0 ? 0 : g+h*(k))    /* k-symbol indel cost */
static int *sapp;                /* Current script append ptr */
static int  last;                /* Last script op appended */
                       /* Append "Delete k" op */
#define DEL(k)                \
{ if (last < 0)                \
    last = sapp[-1] -= (k);        \
  else                    \
    last = *sapp++ = -(k);        \
}
                        /* Append "Insert k" op */
#define INS(k)                \
{ if (last < 0)                \
    { sapp[-1] = (k); *sapp++ = last; }    \
  else                    \
    last = *sapp++ = (k);        \
}
 
#define REP { last = *sapp++ = 0; }        /* Append "Replace" op */

int myers_miller_pair_wise (Alignment *A,int *ns, int **l_s,Constraint_list *CL )
	{
	int **pos;
	int a,b, i, j, l,l1, l2, len;
	int *S;
	char ** char_buf;
	int score;

	
	/********Prepare Penalties******/

	/********************************/
	

	pos=aln2pos_simple ( A,-1, ns, l_s);
	

	l1=strlen (A->seq_al[l_s[0][0]]);
	l2=strlen (A->seq_al[l_s[1][0]]);
	S=vcalloc (l1+l2+1, sizeof (int));
	last=0;
	sapp=S;
        score=diff (A,ns, l_s, 0, l1, 0, l2, 0, 0, CL, pos);
	diff (NULL,ns, l_s, 0, l1, 0, l2, 0, 0, CL, pos); 


	i=0;j=0;sapp=S; len=0;
	while (!(i==l1 && j==l2))
	      {
		  if (*sapp==0){i++; j++;len++;}
		  else if ( *sapp<0){i-=*sapp;len-=*sapp;}
		  else if ( *sapp>0){j+=*sapp;len+=*sapp;}		  
		  sapp++;
	      }
	


	A=realloc_aln2  ( A,A->max_n_seq,len+1);
	char_buf=declare_char (A->max_n_seq,len+1);
	
	i=0;j=0;sapp=S; len=0;
	while (!(i==l1 && j==l2))
	      { 

		   if (*sapp==0)
		      {
			  for (b=0; b< ns[0]; b++)
			      char_buf[l_s[0][b]][len]=A->seq_al[l_s[0][b]][i];
			  for (b=0; b< ns[1]; b++)
			      char_buf[l_s[1][b]][len]=A->seq_al[l_s[1][b]][j];
			  i++; j++;len++;
		      }
		  else if ( *sapp>0)
		        {
			    l=*sapp;
			    for ( a=0; a<l; a++, j++, len++)
			        {
				for (b=0; b< ns[0]; b++)
				    char_buf[l_s[0][b]][len]='-';
				for (b=0; b< ns[1]; b++)
				    char_buf[l_s[1][b]][len]=A->seq_al[l_s[1][b]][j];
				}
			}
		  else if ( *sapp<0)
		        {
			    l=-*sapp;
			    for ( a=0; a<l; a++, i++, len++)
			        {
				for (b=0; b< ns[0]; b++)
				    char_buf[l_s[0][b]][len]=A->seq_al[l_s[0][b]][i];;
				for (b=0; b< ns[1]; b++)
				    char_buf[l_s[1][b]][len]='-';
				}		        
			}
		   
		  sapp++;
	      }
	
	
	A->len_aln=len;
	A->nseq=ns[0]+ns[1];
	
	for ( a=0; a< ns[0]; a++){char_buf[l_s[0][a]][len]='\0'; sprintf ( A->seq_al[l_s[0][a]], "%s", char_buf[l_s[0][a]]);}
	for ( a=0; a< ns[1]; a++){char_buf[l_s[1][a]][len]='\0'; sprintf ( A->seq_al[l_s[1][a]], "%s", char_buf[l_s[1][a]]);}

	free (S);
	free_char ( char_buf, -1);
	l1=strlen (A->seq_al[l_s[0][0]]);
	l2=strlen (A->seq_al[l_s[1][0]]);
	if ( l1!=l2) exit(0);
		
	return score;
	}


int diff (Alignment *A, int *ns, int **l_s, int s1, int M,int s2, int N , int tb, int te, Constraint_list *CL, int **pos)
        {
	 static int *CC;
	 static int *DD;
	     /* Forward cost-only vectors */
	 static int *RR;
	 static int *SS;
             /* Reverse cost-only vectors */
         int   midi, midj, type;    /* Midpoint, type, and cost */
         int  midc;
	 
/*TREATMENT OF THE TERMINAL GAP PENALTIES*/
/*TG_MODE=0---> gop and gep*/
/*TG_MODE=1---> ---     gep*/
/*TG_MODE=2---> ---     ---*/







	 if ( !CC)
	    {
		CC=vcalloc ( M+N+1, sizeof (int));
		DD=vcalloc ( M+N+1, sizeof (int));
		RR=vcalloc ( M+N+1, sizeof (int));
		SS=vcalloc ( M+N+1, sizeof (int));
	    }
	 if ( A==NULL)
	    {
		free(CC);
		free(DD);
		free(RR);
		free(SS);
		CC=DD=RR=SS=NULL;
		return 0;
	    }
	 {
	 int   i, j;
	 int   c, e, d, s;
         int t, g,h,m;
	 
	 g=CL->gop*SCORE_K;
	 h=CL->gep*SCORE_K;
	 m=g+h;
	
	 if (N <= 0){if (M > 0) DEL(M);return gap(M);}
	 if (M <= 1)
	       { 

	       if (M <= 0)
	          {INS(N);
		  return gap(N);
		  }

	       
	       if (tb > te) tb = te;
	       midc = (tb+h) + gap(N);
	       midj = 0;
    
	       for (j = 1; j <= N; j++)
	           { c = gap(j-1) +(CL->get_dp_cost) (A, pos, ns[0], l_s[0],s1, pos, ns[1], l_s[1],j-1+s2,CL)+ gap(N-j);
		   if (c > midc)
		      { midc = c;
		        midj = j;
		      }
		   }
	       if (midj == 0)
	          {DEL(1) INS(N)}
	       else
	          {if (midj > 1) INS(midj-1);
		  REP;
		  if (midj < N) INS(N-midj);
		  }
	       
	       return midc;
	       }
/* Divide: Find optimum midpoint (midi,midj) of cost midc */
 

	 midi = M/2;            /* Forward phase:                          */
	 CC[0] = 0;            /*   Compute C(M/2,k) & D(M/2,k) for all k */
	 t = tb;  
	 for (j = 1; j <= N; j++)
	     { CC[j] = t = t+h;
	       DD[j] = t+g;
	     }
	 t = tb;
	 for (i = 1; i <= midi; i++)
	     { 
	     s = CC[0];
	     CC[0] = c = t = t+h;
	     e = t+g;

	     for (j = 1; j <= N; j++)
	         { 
		 if ((c =   c   + m) > (e =   e   + h)) e = c;
		 if ((c = CC[j] + m) > (d = DD[j] + h)) d = c;
		 c = s + (CL->get_dp_cost) (A, pos, ns[0], l_s[0],i-1+s1, pos, ns[1], l_s[1],j-1+s2,CL);
		 if (e > c) c = e;
		 if (d > c) c = d;
		 s = CC[j];
		 CC[j] = c;
		 DD[j] = d;
		 }
	     }
	 DD[0] = CC[0];
	 
	 RR[N] = 0;            /* Reverse phase:                          */
	 t = te;
	 
	 
	 for (j = N-1; j >= 0; j--)
	     { RR[j] = t = t+h;
	     SS[j] = t+g;
	     }
	 t = te;
	 for (i = M-1; i >= midi; i--)
	     { s = RR[N];
	     RR[N] = c = t = t+h;
	     e = t+g;
	     for (j = N-1; j >= 0; j--)
	          { 
		  if ((c =   c   + m) > (e =   e   + h)) e = c;
		  if ((c = RR[j] + m) > (d = SS[j] + h)) d = c;
		  c = s + (CL->get_dp_cost) (A, pos, ns[0], l_s[0],i+s1, pos, ns[1], l_s[1],j+s2,CL);
		  if (e > c) c = e;
		  if (d > c) c = d;
		  s = RR[j];
		  RR[j] = c;
		  SS[j] = d;
		
		  }
	     }
	 SS[N] = RR[N];	 
	 midc = CC[0]+RR[0];        /* Find optimal midpoint */
	 midj = 0;
	 type = 1;
	 for (j = 0; j <= N; j++)
	     if ((c = CC[j] + RR[j]) >= midc)
		 if (c > midc || CC[j] != DD[j] && RR[j] == SS[j])
		    { 
		    midc = c;
		    midj = j;
		    }
	 for (j = N; j >= 0; j--)
	     if ((c = DD[j] + SS[j] - g) > midc)
	         {midc = c;
		 midj = j;
		 type = 2;
		 }
	 }	    
/* Conquer: recursively around midpoint */
 
  if (type == 1)
    { 
	
    diff (A,ns, l_s, s1,midi, s2, midj, tb, CL->gop*SCORE_K, CL, pos); 
    diff (A,ns, l_s, s1+midi,M-midi, s2+midj, N-midj, CL->gop*SCORE_K,te, CL, pos); 
    }
  else
    { 
      diff (A,ns, l_s, s1,midi-1, s2, midj, tb,0, CL, pos); 
      DEL(2);
      diff (A,ns, l_s, s1+midi+1, M-midi-1,s2+midj, N-midj,0,te, CL, pos); 
    }
  return midc;
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
	int *dd, *ddg;
	int lenal[2], len;
	int t, c,s, e,eg, ch;
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
	gop=CL->gop;
	gep=CL->gep;
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
	dd= vcalloc (len, sizeof (int));
	cc= vcalloc (len, sizeof (int));
	ddg=vcalloc (len, sizeof (int));
	

	
	if ( len>=MAX_LEN_FOR_DP)
	    {
	    long_trace_flag=1;
	    long_trace=vtmpfile();
	    }
	else
	    {
	   
	    dim=(trace==NULL)?0:read_size_int ( trace[-1],0);	   
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

	gop=CL->gop*SCORE_K;
	gep=CL->gep*SCORE_K;
	
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
			e=t+gop;
			tr[0]=(TRACE_TYPE)1;

			

			for (eg=0,j=1; j<=lenal[1];j++)
				{				   

				sub=(CL->get_dp_cost) (A, pos0, ns[0], l_s[0], i-1, pos0, ns[1], l_s[1],j-1,CL);	
				
				/*get the best Insertion*/
				l_gop=(i==lenal[0])?((TG_MODE==0)?gop:0):gop;
				l_gep=(i==lenal[0])?((TG_MODE==2)?0:gep):gep;
			

				if ( a_better_than_b ( e,c+l_gop, maximise))eg++;
				else eg=1;	
				e=best_of_a_b (e, c+l_gop, maximise)+l_gep;
				
				/*Get the best deletion*/
				l_gop=(j==lenal[1])?((TG_MODE==0)?gop:0):gop;
				l_gep=(j==lenal[1])?((TG_MODE==2)?0:gep):gep;
				

				if ( a_better_than_b ( dd[j], cc[j]+l_gop, maximise))ddg[j]++;
				else ddg[j]=1;
				dd[j]=best_of_a_b( dd[j], cc[j]+l_gop,maximise)+l_gep;
				
				c=best_int(3,maximise,&fop, e, s+sub,dd[j]);
				/*Chose Substitution for tie breaking*/
				if ( fop==0 && (s+sub)==e)fop=1;
				if ( fop==2 && (s+sub)==dd[j])fop=1;

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
	

	free ( cc);
	free (dd);		
	free (ddg);
	free (buffer);
	free (char_buf); 
	free (buf_trace);
	free_char ( al, -1);
	if ( long_trace_flag)fclose (long_trace);	

	return score;
	}
     

int sw_pair_wise (Alignment *A, int gop,int gep, int scale,int*ns, int **l_s,Constraint_list *CL, int maximise )
	{
/*******************************************************************************/
/*                SMITH AND WATERMAN                                           */
/*                                                                             */
/*	makes DP between the the ns[0] sequences and the ns[1]                 */
/*                                                                             */
/*	for MODE, see the function  get_dp_cost                                */
/*******************************************************************************/
	int a, b, d, i, j;
	int last_i, last_j;

	int *cc;
	int *dd, *ddg;
	int lenal[2], len;
	
	int c,s, e,eg, ch,g,h;
	int sub;
	
	int fop;
	static int **pos0;
	
	char **al;
	char **aln;

	
	int ala, alb,LEN;
	char *buffer;
	char *char_buf;
	
/*trace back variables       */
	int best_i;
	int best_j;
	int best_score;
	
	
	FILE       *long_trace=NULL;
	TRACE_TYPE *buf_trace=NULL;
	static TRACE_TYPE **trace;
	TRACE_TYPE k;
	TRACE_TYPE *tr;
	int long_trace_flag=0;
	int dim;

/********************************/	
/*CLEAN UP AFTER USE*/
	if ( A==NULL)
	   {
	   free_int (trace,-1);
	   
	   return 0;
	   }	   
/*DO MEMORY ALLOCATION FOR SW DP*/
	
	lenal[0]=strlen (A->seq_al[l_s[0][0]]);
	lenal[1]=strlen (A->seq_al[l_s[1][0]]);
	len= (( lenal[0]>lenal[1])?lenal[0]:lenal[1])+1;
	 buf_trace=vcalloc ( len, sizeof (TRACE_TYPE));	
	 buffer=vcalloc ( 2*len, sizeof (char));	
         al=declare_char (2, 2*len);  
	
	 char_buf= vcalloc (2*len, sizeof (char));	
	 dd= vcalloc (len, sizeof (int));
	 cc= vcalloc (len, sizeof (int));
	 ddg=vcalloc (len, sizeof (int));
	
	
	if ( len>=MAX_LEN_FOR_DP)
	    {
	    
	    long_trace_flag=1;
	    long_trace=vtmpfile();	   
	    }
	else
	    {
	    dim=(trace==NULL)?0:read_size_int ( trace[-1],0);
	    trace    =realloc_int ( trace,dim,dim,len-dim, len-dim);
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

	
	g=gop+scale;
	h=gep+scale;

	cc[0]=0;

	best_score=0;
	best_i=0;
	best_j=0;
	
	tr=(long_trace_flag)?buf_trace:trace[0];
	tr[0]=(TRACE_TYPE)UNDEFINED;
	for ( j=1; j<=lenal[1]; j++)
		{
		cc[j]=0;
		dd[j]=0;
		tr[j]=(TRACE_TYPE)UNDEFINED;
		}
	if (long_trace_flag)fwrite (buf_trace, sizeof ( TRACE_TYPE),lenal[1]+1, long_trace);
	


	for (i=1; i<=lenal[0];i++)
			{
			tr=(long_trace_flag)?buf_trace:trace[i];
			s=cc[0];
			cc[0]=c=0;
			e=0;
			tr[0]=(TRACE_TYPE)UNDEFINED;

			for (eg=0,j=1; j<=lenal[1];j++)
				{
				
				sub=get_sw_dp_cost (A, pos0, ns[0], l_s[0], i-1, pos0, ns[1], l_s[1],j-1,CL,gep, scale);
								
				if ( a_better_than_b ( e, c+g, maximise))
					eg++;
				else 
					eg=1;	
					
				e=best_of_a_b (e, c+g, maximise)+h;
				
				if ( a_better_than_b ( dd[j], cc[j]+g, maximise))
					ddg[j]++;
				else 
					ddg[j]=1;
					
				dd[j]=best_of_a_b( dd[j], cc[j]+g, maximise)+h;
				
			
				if ( sub!=UNDEFINED)c=best_int(4,maximise,&fop, e, (s+sub), dd[j],0);
				else 
				   {
				   c=0;
				   fop=3;
				   dd[j]=e=0;
				   eg=ddg[j]=0;
				   }

				if ( c>best_score)
				   {
				   best_i=i;
				   best_j=j;
				   best_score=c;
				   }
				fop-=1;
				s=cc[j];
				cc[j]=c;			  
				
			
				if ( fop==-1)
					{tr[j]=(TRACE_TYPE)fop*eg;
					}
				else if ( fop==1)
				        {tr[j]=(TRACE_TYPE)fop*ddg[j];
					}
				else if (fop==0)
					{tr[j]=(TRACE_TYPE)0;	
					}	
				else if ( fop==2)
				        {
					tr[j]=(TRACE_TYPE)UNDEFINED;
					}
				
				fop= -2;
				}
			if (long_trace_flag)
			    {
			    fwrite ( buf_trace, sizeof (TRACE_TYPE), lenal[1]+1, long_trace);
			    }
			}
	
	

	
	if (best_i==0 ||best_j==0 )
	    {
	    free (buf_trace);
	    free (buffer);
	    free_char ( al,-1);
	    free ( char_buf);
	    free ( dd);
	    free ( cc);		
	    free ( ddg);
	    A->len_aln=0;
	    aln=A->seq_al;
	    for ( c=0; c< 2; c++)
	        {
		for ( a=0; a< ns[c]; a++) 
		    {
		    aln[l_s[c][a]][0]='\0';
		    }
		}
	    if ( long_trace_flag)fclose ( long_trace);
	    return UNDEFINED;
	    }
	else
	    {
	    i=last_i=best_i;
	    j=last_j=best_j;
	    }
	ala=alb=0;
	
	
	while (i>0 && j>0)
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
				last_i=i;
				last_j=j;
				
				}
			else if (k==(TRACE_TYPE)UNDEFINED)
			        {
				i=0;
				j=0;
				
				}
			else if (k>0)
				{
				
				for ( a=0; a< k; a++)
					{
					al[0][ala++]=1;
					al[1][alb++]=0;
					i--;
					}
				last_i=i;
				last_j=j;
				}
			else if (k<0)
				{
				
				for ( a=0; a>k; a--)
					{
					al[0][ala++]=0;
					al[1][alb++]=1;
					j--;
					}
				last_i=i;
				last_j=j;
				}
			}
      
	LEN=ala;	
	c=LEN-1;  
	
	

	invert_list_char ( al[0], LEN);
	invert_list_char ( al[1], LEN);
	
	if ( A->declared_len<=LEN)realloc_alignment  ( A, 2*LEN);	

	
	aln=A->seq_al;
	
	
	for ( c=0; c<2; c++)
	    for ( a=0; a<ns[c]; a++)
	        {
		e=(c==0)?last_i:last_j;
		for ( d=0; d<e; d++)
		    {
		    A->order[l_s[c][a]][1]+=1-is_gap(aln[l_s[c][a]][d]);
		    }
		}
	
		    
	for ( c=0; c< 2; c++)
	    {
	    for ( a=0; a< ns[c]; a++) 
		{
		aln[l_s[c][a]]+=(c==0)?last_i:last_j;
		ch=0;
		for ( b=0; b< LEN; b++)
		    {
		   
		    if (al[c][b]==1)
			char_buf[b]=aln[l_s[c][a]][ch++];
		    else
			char_buf[b]='-';
		   }
		char_buf[b]='\0';
		aln[l_s[c][a]]-=(c==0)?last_i:last_j;
		sprintf (aln[l_s[c][a]],"%s", char_buf);
	        }
	     }
	
	
	A->len_aln=LEN;

	free ( cc);
	free (dd);		
	free (ddg);
	free (buffer);
	free (char_buf); 
	free (buf_trace);
	if ( long_trace_flag)fclose (long_trace);	
	

	return best_score;
	}

/*******************************************************************************/
/*                LOCAL ALN PROCESSING                                         */
/*                                                                             */
/*	                                                                       */
/*                                                                             */
/*	for MODE, see the function  get_dp_cost                                */
/*******************************************************************************/
Alignment * get_best_local_aln ( Alignment *IN,Constraint_list *CL,int gop, int gep, int sw_t, int sw_l, int sw_z, int greedy)
     {
     int a;
     
     int score,best_score;
     int best_a;
     int prev_end;
     static Alignment *A;
   
    
     for (score=0, best_score=0,best_a=0, a=0; a< IN->len_aln;a++)
         {
	 A=copy_aln (IN,A);
	 A=extract_nol_local_aln(A, a, A->len_aln);
	 
	 if ( A->len_aln+a!=prev_end)
	    {
	    
	    score=evaluate_aln(A, 0,A->len_aln,CL,gop,gep,WE);
	 
	    if ( score>=best_score)
	       {
	       best_score=score;
	       best_a=a;
	       }
	    prev_end=a+A->len_aln;
	    }   
	    
	 }
    
    
     return extract_nol_local_aln(IN, best_a, IN->len_aln);  
    
     }     
Alignment * get_best_nol_local_aln ( Alignment *IN, Constraint_list *CL, int gop, int gep,int sw_t,int sw_l, int sw_z, int mode)
     {
     static Alignment *A;
     int a, score, best_score, best_a;
     int end;

       
     return trunkate_local_aln(IN);
     if ( A->len_aln==IN->len_aln)return IN;
     else
         {
	 end =IN->len_aln;
	 for (best_score=0, best_a=0, a=0; a< IN->len_aln;)
             {

	     if ( a==IN->len_aln)break;
	     
	     A=copy_aln( IN, A);
	     A=extract_aln(A, a, A->len_aln);
	     A=trunkate_local_aln(A);
	     if (A->len_aln+a !=end)
	        {
		end=A->len_aln+a;	
		score=evaluate_aln(A, 0, A->len_aln,CL, gop, gep, WE);
		if ( score > best_score)
	           {
		   best_score=score;
		   best_a=a;
		   }
		}
	     a+=A->len_aln;
	     }
	 

	 IN=extract_aln (IN, best_a, IN->len_aln);
	 return trunkate_local_aln(IN);
	 }
     
     }

/*******************************************************************************/
/*                AUTOMATIC GEP+SCALE PENALTY FOR SW                           */
/*                                                                             */
/*	                                                                       */
/*                                                                             */
/*	                                                                       */
/*******************************************************************************/

double compute_penalty   (Constraint_list *CL, char *mode, int len)
     {


     double p;
     int **bin;
     int coor;
     int n;

     bin=bin_list (CL, WE,0);
     n=read_size_int ( bin[-1],0);
     p=return_max_coor_int(bin, n, 3, &coor);     
     
     
     p=bin[coor][3];
     free_int (bin, -1);
     
     return p;
     }



double compute_scale ( Constraint_list *CL,char *mode, int len)
     {
     double p;
     int n;
     int **bin;
     int coor;

    
     bin=bin_list ( CL,WE,0);
     n=read_size_int ( bin[-1],0); 
     p=return_max_coor_int(bin, n, 3, &coor);     
     p=bin[1][3];
  
     free_int (bin, -1);
     return p;
          
     }
	 




int evaluate_penalty (Alignment *A, Constraint_list *CL, int *scale,char *scale_mode, int *penalty, char *penalty_mode, int len_seq)
    {
    Constraint_list *NCL;

    if ( scale[0]!=0 && penalty[0]!=0)return 1;
    else
        {
	if (CL->M)
	    {
	    scale[0]=0;
	    penalty[0]=-0.5*SCORE_K;
	    }
	else
	    {
	    
	    NCL=duplicate_constraint_list ( CL);
	    NCL=mask_list_with_aln_pair ( A, 0, (A==NULL)?0:A->len_aln,NCL, UNDEFINED);
	    scale[0]=compute_scale ( NCL, scale_mode, len_seq)*-SCORE_K;	    
	    penalty[0]=compute_penalty ( NCL, penalty_mode, len_seq)*-SCORE_K;	
	    free_constraint_list (NCL);   
	    }
	return 1;
	}
    }

/*******************************************************************************/
/*                COFFEE LALIGN                                                */
/*                                                                             */
/*	                                                                       */
/*                                                                             */
/*	                                                                       */
/*******************************************************************************/
Alignment ** t_coffee_lalign   (Constraint_list *CL, int scale, int penalty,int maximise,Sequence *S, int sw_t, int sw_l, int sw_z,int *sw_n, int sw_io)
       {


	   Alignment * aln=NULL;
	   Alignment **aln_list=NULL;
	   int n;
	   Constraint_list *NCL;
	   
	   NCL=duplicate_constraint_list(CL);
	   penalty=scale=0;	   
	   evaluate_penalty (aln,NCL, &scale, "default", &penalty, "default",  strlen ( S->seq[0]));			  	      	
	  
	   n=0;
	   while ( (aln==NULL) || ((aln->score_aln>sw_z)|| (sw_n[0]!=0 && sw_n[0]>n)))
	         {
		 n++;
		 if ( aln!=NULL)
		     {
		     aln_list=realloc_aln_array(aln_list,1);
		     aln_list[(aln_list[-1])->nseq-1]=copy_aln(aln,aln_list[(aln_list[-1])->nseq-1] );	 
		     aln=free_aln(aln);		     
		     aln=NULL;
		     }

		 aln=add_seq2aln ( NCL,scale, penalty,aln, maximise,S,sw_t, sw_l, sw_z,sw_n,sw_io);
		 aln=add_seq2aln ( NCL,scale, penalty,aln, maximise,S,sw_t, sw_l, sw_z,sw_n,sw_io);
		 if ( aln==NULL)break; 

		 if ( sw_io)output_aln_with_res_number (aln, stderr);
		 NCL=mask_list_with_aln_pair (aln, 0, aln->len_aln,NCL,UNDEFINED);
		 aln->score_aln=evaluate_seq_z_score (aln, 0, aln->len_aln)/aln->len_aln;
		 fprintf ( stderr, "\nZ_SCORE=%d\n", aln->score_aln);
		 }
       free_constraint_list(NCL);
       return aln_list;
       }

Alignment * add_seq2aln   (Constraint_list *CL,int scale,int penalty, Alignment *IN,int maximise,Sequence  *S,int sw_t, int sw_l, int sw_z,int *sw_n, int sw_io)
           {	
	   int *n_groups;
	   int **group_list;
	   int a;
	   static int series=0;
	   static int ref_score;
	   int score;
	   int pass;
	   int trunkate=1;



	   if ( IN==NULL)
	      {
	      IN=realloc_alignment2(IN, 1, strlen (S->seq[0])+1);
	      IN->S=S;
	      IN->nseq=1;
	      sprintf ( IN->seq_al[0], "%s", S->seq[0]);
	      sprintf (IN->name[0], "%s_%d_1", S->name[0],series);
	      ref_score=0;
	      IN->len_aln=strlen ( IN->seq_al[0]);
	      series++;
	      }	   
	   else
	      {
	      
	      IN=realloc_alignment2 ( IN, IN->nseq+1,MAX(strlen ( S->seq[0])+1, IN->len_aln+1));
	      n_groups=vcalloc ( 2, sizeof (int));
	      group_list=declare_int (2,IN->nseq+1);
	      
	      n_groups[0]=IN->nseq;
	      for ( a=0; a<IN->nseq; a++)group_list[0][a]=a;
	      n_groups[1]=1;
	      group_list[1][0]=IN->nseq;
	      sprintf (IN->name[IN->nseq], "%s_%d_%d",S->name[0],series,IN->nseq+1);
	      sprintf (IN->seq_al[IN->nseq], "%s",S->seq[0]);
	      IN->order[IN->nseq][0]=0;
	      IN->order[IN->nseq][1]=0;
	      IN->nseq++;
	      sw_pair_wise ( IN, penalty,penalty, scale, n_groups, group_list,CL,maximise);		      	     
	      
	      if ( trunkate==1)
	         {
		 IN=trunkate_local_aln(IN);        
		 }
	      else if ( trunkate==0)
	         {
		 IN=get_best_local_aln(IN,CL, scale, penalty, sw_t, sw_l, sw_z,NON_GREEDY);
		 }
	      
	      free (n_groups);
	      free_int ( group_list,-1);	      
	      
	      score=IN->score_aln=evaluate_seq_z_score(IN,0, IN->len_aln);
	      pass =(score>sw_z)&&(score>sw_l)&&((score*100)/MAX(1,ref_score)>sw_t);
	      if ( sw_n[series-1]>0)pass=(sw_n[series-1]>=IN->nseq);

	      
	      if (IN->nseq>2 && !pass)IN->finished=1;
	      else if (IN->nseq<=2 && !pass)IN->finished=2; 
	      else if ( pass)
		  ref_score=score;
	      
	      }	   	   
	   return IN;
	   }

int old_pair_wise (Alignment *A, int*ns, int **l_s,Constraint_list *CL)
	{
/*******************************************************************************/
/*                SMITH AND WATERMAN                                           */
/*                                                                             */
/*	makes DP between the the ns[0] sequences and the ns[1]                 */
/*                                                                             */
/*	for MODE, see the function  get_dp_cost                                */
/*******************************************************************************/
	

	int a, b, i, j;
	int *cc;
	int *dd, *ddg;
	int lenal[2], len;
	int gop, gep, maximise;
	int t, c,s, e,eg, ch,g,h;
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
	gop=CL->gop;
	gep=CL->gep;
	
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
	dd= vcalloc (len, sizeof (int));
	cc= vcalloc (len, sizeof (int));
	ddg=vcalloc (len, sizeof (int));
	

	fprintf ( stderr, "\nUSE THE OLD PAIR WISE");
	
	if ( len>=MAX_LEN_FOR_DP)
	    {
	    long_trace_flag=1;
	    long_trace=vtmpfile();
	    }
	else
	    {
	    
	    dim=(trace==NULL)?0:read_size_int ( trace[-1],0);	   
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
	g=gop;
	h=gep;
	t=h;
	cc[0]=0;
		
	tr=(long_trace_flag)?buf_trace:trace[0];
	tr[0]=(TRACE_TYPE)1;
	for ( j=1; j<=lenal[1]; j++)tr[j]=(TRACE_TYPE)-1;
	if (long_trace_flag)fwrite (buf_trace, sizeof ( TRACE_TYPE),lenal[1]+1, long_trace);

	t=g=gop;
	h=gep;
	for (cc[0]=0,j=1; j<=lenal[1]; j++)
	    {
	    cc[j]=t=t+h;
	    dd[j]=t+g;
	    }
	t=g;
	for (i=1; i<=lenal[0];i++)
			{
			tr=(long_trace_flag)?buf_trace:trace[i];
			s=cc[0];
			cc[0]=c=t=t+h;
			e=t+h;
			tr[0]=(TRACE_TYPE)1;

			for (eg=0,j=1; j<=lenal[1];j++)
				{
				sub=get_dp_cost (A, pos0, ns[0], l_s[0], i-1, pos0, ns[1], l_s[1],j-1,CL);										
				if ( a_better_than_b ( e, c+g, maximise))
					eg++;
				else 
					eg=1;	
					
				e=best_of_a_b (e, c+g, maximise)+h;
				
				if ( a_better_than_b ( dd[j], cc[j]+g, maximise))
					ddg[j]++;
				else 
					ddg[j]=1;
				dd[j]=best_of_a_b( dd[j], cc[j]+g, maximise)+h;
				c=best_int(3,maximise,&fop, e, s+sub,dd[j]);
				fop-=1;
				s=cc[j];
				cc[j]=c;			  
				if ( fop<0)
				        {
				        tr[j]=(TRACE_TYPE)fop*eg;
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

	free ( cc);
	free (dd);		
	free (ddg);
	free (buffer);
	free (char_buf); 
	free (buf_trace);
	if ( long_trace_flag)fclose (long_trace);	

	return score;
	}
