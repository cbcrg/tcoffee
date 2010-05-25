#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>


#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"

#include "dp_lib_header.h"





int sw_pair_wise (Alignment *A, int gep, int gop,int*ns, int **l_s,int **L, int ne, int maximise )
	{
/*******************************************************************************/
/*	makes DP between the the ns[0] sequences and the ns[1] sequences in A  */
/*                                                                             */
/*	for MODE, see the function  get_dp_cost                                */
/*******************************************************************************/
	
	int a, b, d, f, i, j;
	int *list;
	int *cc;
	int *dd, *ddg;
	int lenal[2], len;
	int nseq;
	int hD, gD, hI, gI,t, c,s, e,eg, k, ch,g,h;
	int sub;
	
	int fop;
	int score=0;
	
	static int **pos1;
	static int **pos0;
	
	
	static int **al;
	static char **aln;
	static int **trace;
      	static int **score_mat;
	static max_len;

	int ala, alb,LEN;
	int *buffer;
	int *char_buf;
	
	nseq=(A->S)->nseq;
	lenal[0]=strlen (A->seq_al[l_s[0][0]]);
	lenal[1]=strlen (A->seq_al[l_s[1][0]]);
	
	len= ( lenal[0]>lenal[1])?lenal[0]:lenal[1];
	
	
	
	
	if ( len>=max_len)
		{
		if ( trace!=NULL)
			{
			free_int(trace,max_len);
			free_int(score_mat, max_len);
			free_char (al, nseq);
			
			}
		trace    =declare_int ( len+1, len+1);
		score_mat=declare_int ( len+1, len+1);
		al=declare_char (nseq, len*2+1);
		max_len=len+1;
		
		}
	buffer=calloc ( 2*max_len+1, sizeof (char));
	char_buf= calloc ( 2*max_len+1, sizeof (char));	
	dd=calloc (max_len+1, sizeof (int));
	cc=calloc (max_len+1, sizeof (int));
	ddg=calloc (max_len+1, sizeof (int));
	
	pos0=aln2pos_simple ( A, nseq);
	
	
	
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
		       
	for (i=0;i<=lenal[0];i++)
		{trace[i][0]=UNDEFINED;
		}
	for (j=0;j<=lenal[1]; j++)
		{trace[0][j]=UNDEFINED;
		}
	
	g=gop;
	h=gep;
	
	trace[0][0]=1;	
	cc[0]=0;
	
	t=0;
	for ( j=1; j<=lenal[1]; j++)
		{
		cc[j]=t=0;
		dd[j]=0;
		}
		
	
	t=0;			
	for (i=1; i<=lenal[0];i++)
			{
			s=cc[0];
			cc[0]=c=t=0;
			e=0;	
			for (eg=0,j=1; j<=lenal[1];j++)
				{
				sub=(CL->get_dp_cost) (A, pos0, ns[0], l_s[0], i-1, pos0, ns[1], l_s[1],j-1,L, ne);
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
				
				c=best_int(4,maximise,&fop, e, s+sub, dd[j],0, maximise);
				fop-=1;
				s=cc[j];
				cc[j]=c;
				
				score_mat[i][j]=c;
				if ( fop==-2)
				   {
				   trace[i][j]==UNDEFINED;
				   }
				else if ( fop<0)
					{trace[i][j]=fop*eg;
					}
				else if ( fop>0)
					{trace[i][j]=fop*ddg[j];
					 }
				else 
					{trace[i][j]=0;	
					 }	
				fop= -2;
				}
			}
	
	
	score=c;	
	
	return_2D_max_coor ( score_mat, lena_al[0], len_al[1], &i, &j);
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
				k=trace[i][j];
				
				
				
			if (k==0)
				{
				al[0][ala++]=i;
				al[1][alb++]=j;
				
				i--;
				j--;
				}
			else if (k>0)
				{
				for ( a=0; a< k; a++)
					{
					al[0][ala++]=i;
					al[1][alb++]=0;
					i--;
					}
				}
			else if (k<0)
				{
				for ( a=0; a>k; a--)
					{
					al[0][ala++]=0;
					al[1][alb++]=j;
					j--;
					}
				}
			else if ( k==UNDEFINED)
			        {
				i=j=0;
				}
			}
	
	
	LEN=ala;	
	c=LEN-1;  
	if ( A->declared_len<=LEN)realloc_alignment ( A, 2*LEN);  
	aln=A->seq_al;
	
	invert_list_int (al[0], LEN);
	invert_list_int (al[1], LEN);
	
	
        for ( c=0; c< 2; c++)
	    {
	    for ( a=0; a< ns[c]; a++) 
		{
		ch=0;
		for ( b=0; b< LEN; b++)
		    {
		   
		    if (al[c][b]>0)
			char_buf[b]=aln[l_s[c][a]][al[c][b]-1];
		    else
			char_buf[b]='-';
		   }
		char_buf[b]='\0';
		sprintf (aln[l_s[c][a]],"%s", char_buf);
	        }
	     }
	
	A->len_aln=LEN;
	
	a=strlen (char_buf);
	
	for ( a=0; a< 2; a++)fprintf ( stderr, "\n%s\n%s",aln[l_s[c][0]], aln[l_s[c][1]]); 

	free ( cc);
	free (dd);		
	free (ddg);
	free (buffer);
	free (char_buf); 
	
	return score;
	}
	
