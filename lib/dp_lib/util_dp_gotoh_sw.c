#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"

#include "dp_lib_header.h"

int gotoh_pair_wise_lalign ( Alignment *A, int*ns, int **l_s,Constraint_list *CL)
    {
    Alignment *BUF=NULL;
    Alignment *EA=NULL;
    
    int a;
    BUF=copy_aln (A, BUF);
    
    
    for ( a=0; a<CL->lalign_n_top; a++)
      {
	free_aln (A);

	A=copy_aln (BUF, A);
	
	A->score_aln=gotoh_pair_wise_sw (A, ns, l_s, CL);
	EA=fast_coffee_evaluate_output (A, CL);

	output_format_aln (CL->out_aln_format[0],A,EA,"stdout");
	CL=undefine_sw_aln ( A, CL);
      }
    exit (1);
    return 0;
    }
Constraint_list * undefine_sw_aln ( Alignment *A, Constraint_list *CL)
  {
    int a, b, l;
    int **pos;
    int  r1, rs1;
    int  r2, rs2;
    


    pos=aln2pos_simple ( A,A->nseq);
    
    for ( l=0; l< A->len_aln; l++)
      for ( a=0; a< A->nseq-1; a++)
	{
	  rs1=A->order[a][0];
	  r1 =pos[a][l];
		
	  if ( r1<=0)continue;
	  for ( b=a+1; b< A->nseq;b++)
		  {
		    rs2=A->order[b][0];
		    r2 =pos[b][l];
		    if ( r2<=0)continue;
		    
		    CL=undefine_sw_pair ( CL, rs1, r1, rs2, r2);
		  }
	}
    free_int (pos, -1);
    return CL;
  }
Constraint_list * undefine_sw_pair ( Constraint_list *CL, int s1, int r1, int s2, int r2)
  {
    int a, b;
    
    if ( !CL->forbiden_pair_list)
      {
	
	CL->forbiden_pair_list=vcalloc ( (CL->S)->nseq, sizeof (int ***));
	for ( a=0; a< ((CL->S)->nseq); a++)
	  {
	    CL->forbiden_pair_list[a]=vcalloc ( (CL->S)->nseq, sizeof (int **));
	    for ( b=0; b< ((CL->S)->nseq); b++)
	      CL->forbiden_pair_list[a][b]=vcalloc ( (CL->S)->len[a]+1, sizeof (int *));
	    
	  }
      }
    if ( CL->forbiden_pair_list[s1][s2][r1]==NULL)CL->forbiden_pair_list[s1][s2][r1]=vcalloc ( (CL->S)->len[s2]+1, sizeof (int));
    CL->forbiden_pair_list[s1][s2][r1][r2]=1;
    
    if ( CL->forbiden_pair_list[s2][s1][r2]==NULL)CL->forbiden_pair_list[s2][s1][r2]=vcalloc ( (CL->S)->len[s1]+1, sizeof (int));
    CL->forbiden_pair_list[s2][s1][r2][r1]=1;
    
    return CL;
  }
       
int sw_pair_is_defined ( Constraint_list *CL, int s1, int r1, int s2, int r2)
        {
	  int d;
	  
	  d=(r1-r2);
	  d=(d<0)?-d:d;
	  

	  if ( s1==s2 && d<(CL->sw_min_dist)) return UNDEFINED;
	  else if ( ! CL->forbiden_pair_list) return 1;
	  else if ( CL->forbiden_pair_list[s1][s2][r1]==NULL)return 1;
	  else if ( CL->forbiden_pair_list[s1][s2][r1][r2]==1)return UNDEFINED;
	  else if ( CL->forbiden_pair_list[s1][s2][r1][r2]==0)return 1;
	  
	  else 
	    {
	      crash ("ERROR in function: sw_pair_is_defined\n");
	      return UNDEFINED;
	    }
	    
	}	

  
int gotoh_pair_wise_sw (Alignment *A, int*ns, int **l_s,Constraint_list *CL)
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
	int t;
	int *cc;
	int *dd, *ddg;
	int lenal[2], len;
	
	int c,s, e,eg, ch,g,h, maximise;
	int sub;
	
	int fop;
	static int **pos0;
	
	char **al=NULL;
	char **aln=NULL;

	
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
/********Prepare penalties*******/
	if (CL->moca)
	  {
	    g=(CL->gop+(CL->moca)->moca_scale)*SCORE_K;
	    h=(CL->gep+(CL->moca)->moca_scale)*SCORE_K;
	  }
	else
	  {
	    g=(CL->gop-CL->nomatch)*SCORE_K;
	    h=(CL->gep-CL->nomatch)*SCORE_K;
	  }
	fprintf ( stderr, "\n%d %d", g, h);
	maximise=CL->maximise;
/********************************/	
/*CLEAN UP AFTER USE*/
	if ( A==NULL)
	   {
	   free_int (trace,-1);
	   trace=NULL;
	   if ( al)free_char (al,-1);
	   al=NULL;
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
	    dim=(trace==NULL)?0:read_size_int ( trace,sizeof (int*));
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

	

	cc[0]=0;

	best_score=0;
	best_i=0;
	best_j=0;
	
	tr=(long_trace_flag)?buf_trace:trace[0];
	tr[0]=(TRACE_TYPE)UNDEFINED;

	t=g;
	for ( j=1; j<=lenal[1]; j++)
		{
		cc[j]=t=t+h;
		dd[j]=t+g;
		tr[j]=(TRACE_TYPE)UNDEFINED;
		}
	if (long_trace_flag)fwrite (buf_trace, sizeof ( TRACE_TYPE),lenal[1]+1, long_trace);
	

	t=g;
	for (i=1; i<=lenal[0];i++)
			{		       
			tr=(long_trace_flag)?buf_trace:trace[i];
			s=cc[0];
			cc[0]=c=t=t+h;
			e=t+g;
			tr[0]=(TRACE_TYPE)UNDEFINED;

			for (eg=0,j=1; j<=lenal[1];j++)
				{
				
				sub=(CL->get_dp_cost) (A, pos0, ns[0], l_s[0], i-1, pos0, ns[1], l_s[1],j-1,CL);
				
				/*get the best Insertion*/
				if ( a_better_than_b ( e, c+g, maximise))
					eg++;
				else 
					eg=1;	
				e=best_of_a_b (e, c+g, maximise)+h;
				
				/*Get the best deletion*/
				if ( a_better_than_b ( dd[j], cc[j]+g, maximise))
					ddg[j]++;
				else 
					ddg[j]=1;
				dd[j]=best_of_a_b( dd[j], cc[j]+g, maximise)+h;
				
				/*Chose Substitution for tie breaking*/
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
	    vfree (buf_trace);
	    vfree (buffer);
	    free_char ( al,-1);
	    vfree ( char_buf);
	    vfree ( dd);
	    vfree ( cc);		
	    vfree ( ddg);
	    free_int (pos0, -1);
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
			if ( i==0 || j==0)k=UNDEFINED;
			  /*				k=-1;
		 	else if ( j==0)
			 	k=1;
		 	else if ( j==0 && i==0)
		 	k=1;*/	
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
				
			if ( k==UNDEFINED){i=j=0;}	
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
		  A->order[l_s[c][a]][2]=(c==0)?last_i:last_j;
		  A->order[l_s[c][a]][3]=(c==1)?best_i:best_j;
		  
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
	
	free_int (pos0, -1);
	vfree ( cc);
	vfree (dd);		
	vfree (ddg);
	vfree (buffer);
	vfree (char_buf); 
	vfree (buf_trace);
	if ( long_trace_flag)fclose (long_trace);	
	
		
	return best_score;
	}


/*******************************************************************************/
/*                AUTOMATIC GEP+SCALE PENALTY FOR SW                           */
/*                                                                             */
/*	                                                                       */
/*                                                                             */
/*	                                                                       */
/*******************************************************************************/

Alignment * add_seq2aln   (Constraint_list *CL, Alignment *IN,Sequence  *S)
           {	
	   int *n_groups;
	   int **group_list;
	   int a;
	   static int series=0;
	  
	  
	  
	  
	   int ste; /*sequence to extract, last one if they are packed*/





	   if (CL->packed_seq_lu){ste=S->nseq-1;}
	   else{ste=0;}

	   if ( IN==NULL)
	      {
	      IN=realloc_aln2(IN, 1, strlen (S->seq[ste])+1);
	      IN->S=S;
	      IN->nseq=1;
	      
	      
	      
	      sprintf ( IN->seq_al[0], "%s", S->seq[ste]);
	      sprintf (IN->name[0], "%s_%d_1", S->name[ste],series);
	      IN->order[0][0]=ste;
	      IN->order[0][1]=0;
	      
	      IN->len_aln=strlen ( IN->seq_al[0]);
	      series++;
	    
	      }	   
	   else
	      {
	      
	      IN=realloc_aln2 ( IN, IN->nseq+1,MAX(strlen ( S->seq[ste])+1, IN->len_aln+1));
	      n_groups=vcalloc ( 2, sizeof (int));
	      group_list=declare_int (2,IN->nseq+1);
	      
	      n_groups[0]=IN->nseq;
	      for ( a=0; a<IN->nseq; a++)group_list[0][a]=a;
	      
	      n_groups[1]=1;
	      group_list[1][0]=IN->nseq;
	      sprintf (IN->name[IN->nseq], "%s_%d_%d",S->name[ste],series,IN->nseq+1);
	      sprintf (IN->seq_al[IN->nseq], "%s",S->seq[ste]);
	      IN->order[IN->nseq][0]=ste;
	      IN->order[IN->nseq][1]=0;
	      IN->nseq++;
	      
	     
	      pair_wise ( IN, n_groups, group_list,CL);	
	      
	      }
	  
	   return IN;
	      
	   }
	   

