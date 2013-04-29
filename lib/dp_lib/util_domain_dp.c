#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"

#include "dp_lib_header.h"

int domain_pair_wise (Alignment *A,int*in_ns, int **in_l_s,Constraint_list *CL )
	{
/*******************************************************************************/
/*                SEQ_DOMAIN DP                                                    */
/*                                                                             */
/*	makes DP between the the ns[0] sequences and the ns[1]                 */
/*                                                                             */
/*	for MODE, see the function  get_dp_cost                                */
/*******************************************************************************/
	
	int scale, gop, gep, maximise; 
	int a, b, i, j,l,x;
	int best_j;

	int lenal[2], len;
	int sub;
	int match;
	
	int *ns, **l_s;


	int *f;
	int *pf;

	int *dd;
	int *dd_len;
	
	int * e;
	int *pe;
	int * e_len;
	int *pe_len;
	    
	int fop;
        int **pos0;
	
	int  **al=NULL;
	int  **pos_al=NULL;


	
	int ala,LEN;
	char *buffer;
	char *char_buf;


/*trace back variables       */
	TRACE_TYPE *buf_trace=NULL;
	TRACE_TYPE **trace;
	TRACE_TYPE k;
	TRACE_TYPE *tr;
	int **result_aln;
	int nseq;
/*Test Varaibles*/
	int score;

/*Prepare l_s and ns*/

	
	ns=(int*)vcalloc( 2, sizeof (int));
	l_s=(int**)vcalloc ( 2, sizeof (int*));
	
	ns[0]=in_ns[1];
	ns[1]=in_ns[0];
	  
	l_s[0]=in_l_s[1];
	l_s[1]=in_l_s[0];
	lenal[0]=strlen (A->seq_al[l_s[0][0]]);
	lenal[1]=strlen (A->seq_al[l_s[1][0]]);
	len=lenal[0]+lenal[1]+2;
/********************************/	
gop=(CL->gop)*SCORE_K;
gep=(CL->gep)*SCORE_K;
maximise=CL->maximise;
scale=(lenal[1]*SCORE_K*(CL->moca)->moca_scale);
/*******************************/

	

/*DO MEMORY ALLOCATION FOR DP*/
	

	


	buf_trace=(int*)vcalloc ( len, sizeof (TRACE_TYPE));	
	buffer=(char*)vcalloc ( 2*len, sizeof (char));	
	
        al    =declare_int  (2, 2*len);  
	pos_al=declare_int  (2, 2*len);
	result_aln=declare_int (1,len);
	char_buf=(char*)vcalloc (2*len, sizeof (char));	

	f     =(int*)vcalloc (len, sizeof (int));
       pf     =(int*)vcalloc (len, sizeof (int));
	e     =(int*)vcalloc (len, sizeof (int));
       pe     =(int*)vcalloc (len, sizeof (int));
	e_len =(int*)vcalloc (len, sizeof (int));
       pe_len =(int*)vcalloc (len, sizeof (int));
	dd    =(int*)vcalloc (len, sizeof (int));
	dd_len=(int*)vcalloc (len, sizeof (int));

	
	trace=declare_int (lenal[0]+2, lenal[1]+2);

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

	for ( i=0; i<=lenal[0]+1; i++)
	    {	    
	    tr=trace[i];

	    
	    for ( sub=0,j=0; j<=lenal[1]; j++)
	        {
		if      (i==0 && j==0){tr[j]=1;pe[j]=dd[j]=gop;}
		else if (i==0)        {e[j]=pe[j]=dd[j]=gop;dd_len[j]=e_len[j]=pe_len[j]=f[j]=pf[j]=0;tr[j]=-1;}
		else if (j==0)
		     {
		     for (f[j]=pf[0],best_j=0,a=1; a<=lenal[1]; a++)
		         {
			 if (f[j]!=MAX(pf[a]+scale,f[j]))
			     {
			     f[j]=pf[a]+scale;
			     best_j=a;
			     }			
			 }
		     

		     dd    [j]=e[j]=pe[j]=gop;
		     dd_len[j]=e_len[j]=0;
		     tr    [j]=best_j;
		     
		     }
		else if (i>lenal[0]);
		else
		     {
		     sub=(CL->get_dp_cost) (A, pos0, ns[0], l_s[0], i-1, pos0, ns[1], l_s[1],j-1,CL);
		     
		  
		     
		     match=pf[j-1]+sub;
		     if (a_better_than_b(pf[j]+gop, pe[j]+gep,maximise))
			 {
			 e    [j]=pf[j]+gop;
			 e_len[j]=1;
			 }
		     else 
		         {
			 e    [j]=pe    [j]+gep;
			 e_len[j]=pe_len[j]+1;
			 }
		     
		     
		     if (a_better_than_b(f[j-1]+gop, dd[j-1]+gep,maximise))
		        {
			dd    [j]=f[j-1]+gop;
			dd_len[j]=1;
			}
		     else 
		        {
			dd    [j]=dd    [j-1]+gep;
		        dd_len[j]=dd_len[j-1]+1  ;
			}
		     


		     if ( sub!=UNDEFINED)
			 {			 
			 
			 f[j] =best_int(4,maximise,&fop,e[j],match,dd[j],f[0]);
			 fop-=1;
			 
			 if      (fop==-1)fop= e_len[j]*fop;
			 else if (fop==1 )fop=dd_len[j]*fop;
			 else if (fop==2 )fop=UNDEFINED;
			 
			 
			 
			 }
		     else
		         {
			 dd[j]=e[j]=match=-10000;			 
			 f[j]=f[0];
			 fop=UNDEFINED;
			 }
		     pe    [j]  =e    [j];
		     pf    [j-1]=f    [j-1];
		     pe_len[j]  =e_len[j];
		     tr[j]      =fop;	
		     }
	
		}
	    
	    pf    [j-1]=f    [j-1];
	    pe    [j]  =e    [j];
	    pe_len[j]  =e_len[j];
	    }
	score=f[0];
        i=lenal[0]+1;
	j=0;
	ala=0;
	
	
	while (i!=0)
	      {
	     

	      k=trace[i][j];
	      if (j==0 && i<=lenal[0])
	         {
		 
		 pos_al[0][ala]=i;
		 pos_al[1][ala]=j;
		 al[0][ala]=MATCH;
		 al[1][ala]=UNALIGNED;
		 i--;
		 j=k;
		 ala++;
		 }
	      else if ( j==0 && i>lenal[0])
	         {
		 j=k;
		 i--; 
		 
		 }
	      else if (k==0)
	         {		
		 pos_al[0][ala]=i;
		 pos_al[1][ala]=j;
		 al[0][ala]=MATCH;
		 al[1][ala]=MATCH;
		 i--;
		 j--;
		
		 ala++;
		 }		
	      else if (k!=UNDEFINED && k<0 && i>0)
	         {
		 for (x=0; x< -k && i>0; x++)
		     {
		     pos_al[0][ala]=i;
		     pos_al[1][ala]=j;    
		     al[0][ala]=MATCH;
		     al[1][ala]=GAP;		 
		     i--;
		     ala++;
		     }
		 }
	      else if (k!=UNDEFINED && k>0 && j>0)
	         {
		 for ( x=0; x< k && j>0; x++)
		     {
		     pos_al[0][ala]=i;
		     pos_al[1][ala]=j;    
		     al[0][ala]=GAP;
		     al[1][ala]=MATCH;
		     j--;
		     ala++;
		     }
		 }
	      else if ( k==UNDEFINED){j=0;}
	      
	      }

	
	LEN=ala;	

	invert_list_int ( pos_al[0], LEN);
	invert_list_int ( pos_al[1], LEN);		
	invert_list_int (     al[0], LEN);
	invert_list_int (     al[1], LEN);	

	
	/*O: TARGET SEQUENCE  (long)*/
	/*1: PATTERN SEQUENCE (short)*/

	
	
	for ( b=0; b<lenal[1]; b++)result_aln[0][b]=-1;
	for (l=0, nseq=0, a=0; a< LEN; a++)
	    {
	    i=pos_al[0][a];
	    j=pos_al[1][a];
	    
	    if (al[1][a]==UNALIGNED && l>0)
	             {
		     result_aln=realloc_int ( result_aln, read_size_int ( result_aln,sizeof (int*)),len, 1, 0); 	   
		     nseq++;
		     l=0;
		     for ( b=0; b<lenal[1]; b++)result_aln[nseq][b]=-1;
		     }
	    
	    else if (al[1][a]==MATCH && al[0][a]==MATCH ){l++;result_aln[nseq][j-1]=i-1;}
	    else if (al[1][a]==MATCH && al[0][a]==GAP   ){l++;result_aln[nseq][j-1]=-1;}
	    
	    }
	if ( l>0)nseq++;




	A=domain_match_list2aln ( A,ns,l_s,result_aln,nseq,lenal[1]); 

	
	vfree (f);
	vfree (pf);
	vfree (e);
	vfree (pe);
	vfree (e_len);
	vfree (pe_len);
	vfree (dd_len);
	vfree (dd);
	free_int (pos0, -1);
	vfree (buffer);
	vfree (char_buf); 
	vfree (buf_trace);
	free_int ( pos_al, -1);
        pos_al=NULL;
	  
	free_int (     al, -1);

	free_int (trace,-1);
	free_int ( result_aln, -1);
	return score;
	}


Alignment *domain_match_list2aln ( Alignment *A,int *ns,int **l_s,int **ml, int nseq, int len)
           {

	     /*
	       function documentation: start
	       
	       This function edits the alignment given the results obtained by DP
	       ns: ns[0]->number of sequences serarched (TARGET)
	           ns[1]->number of sequences in the pattern (PATTERN SEQ)
	       l_s:   
		   l_s[0]->list of sequences in the TARGET...
	       
	       nseq: number of occurences of PATTERN in TARGET
	       len:  length of the PATTERN
	       
	       ml: detail of the nseq matches
	           ml[x][y]=k-> residue k of TARGET matches residue y of pattern
	                     -> k=-1 means a gap;
			     
	       NOTE: This implementation can only match ONE target sequence with the PATTERN
	       The Pattern can either be one sequence or a profile.

	       function documentation: end
	     */

	     
	     

	   int a, b, c, d, e;
	   Alignment *B=NULL;
	   int **new_ml;
	   int  *max_ml;
	   int  *start_ml;
	   int tot_nseq;
	   int max_len,seq;
	   char *buf;
	   
	   if ( len==0 || nseq==0)
	      {
	      A->nseq=0;
	      A->len_aln=0;
	      }
	   else
	      {
	      B=copy_aln(A, B);
	      /*1 Extract the sequence used as a pattern, put it on the top*/
	      
	    
	      A=shrink_aln (A, ns[1], l_s[1]);
	      A=realloc_aln2(A, ns[1]+ns[0]*nseq,len+A->len_aln+1);
	      
	      
	      
	      new_ml  =declare_int ( nseq, 3*len);
	      max_ml  =(int*)vcalloc ( nseq, sizeof (int));
	      for ( a=0; a<nseq; a++)max_ml[a]=-1;

	      start_ml=(int*)vcalloc ( nseq, sizeof (int));

	      for ( b=0,a=0; a< len; a++, b+=3)
	          {
		  for ( max_len=0,c=0; c< nseq; c++)
		      {
		      if ( max_ml[c]<0 && ml[c][a]<0)
		         {
			 new_ml[c][b]=ml[c][a];
		         new_ml[c][b+1]=0;
			 }
		      else if ( ml[c][a]>=0)
		         {		       
			 new_ml[c][b]=ml[c][a];
			 if ( max_ml[c]<0) start_ml[c]=max_ml[c]=ml[c][a];
			 for ( d=a+1; d<len; d++){if ( ml[c][d]>=0){ max_ml[c]= ml[c][d];break;}}
			 if (max_ml[c]!=new_ml[c][b])
			    {
			    new_ml[c][b+1]=max_ml[c]-new_ml[c][b]-1;
			    }
			 max_len=MAX( max_len, new_ml[c][b+1]);
			 }
		      else
		         {
			 new_ml[c][b]=ml[c][a];
			 new_ml[c][b+1]=0;
			 }
		      }
	       
		  for ( c=0; c< nseq; c++){new_ml[c][b+2]=max_len;}
		  }
	  
	      tot_nseq=ns[1]+ns[0]*nseq;
	   
	      for ( a=0, b=0; a< len ;a++)
	          {
		 
	          /*1: Place the Match Column*/
		  for ( c=0; c< ns[1]; c++)
		      {
		      A->seq_al[c][b]=B->seq_al[l_s[1][c]][a];
		      A->seq_al[c][b+1]='\0';		   
		      }
		  for ( e=0,c=ns[1]; c<tot_nseq; c+=ns[0],e++)
		      for ( d=0; d< ns[0]; d++)
		          {
			  if ( new_ml[e][3*a]!=-1)
			    A->seq_al[c+d][b]=B->seq_al[l_s[0][d]][new_ml[e][3*a]];
			  else 
			      A->seq_al[c+d][b]='-';
			  A->seq_al[c+d][b+1]='\0';
		       
			  }
		  b++;

	         /*2: Add the Gaps before the next_column*/
		  if ( new_ml[0][3*a+2]>0)
	             {
		     for ( c=0; c< ns[1]; c++)
		         {
			 buf=generate_null(new_ml[0][3*a+2]);
			 strcat ( A->seq_al[c],buf);
			 vfree (buf);
			 }
		     for (e=0,c=ns[1];c< tot_nseq; c+=ns[0], e++)
		         {
			   buf=extract_char (B->seq_al[l_s[0][0]], new_ml[e][3*a]+1, new_ml[e][3*a+1]);
			   strcat ( A->seq_al[c],buf);
			   vfree (buf);
			   buf=generate_null(new_ml[e][3*a+2]-new_ml[e][3*a+1]);
			   strcat ( A->seq_al[c],buf);
			   
			   vfree (buf);
			   
			 }
		     }
		  b+=new_ml[0][3*a+2];
		  }
	      
	      for (e=0,a=ns[1]; a< tot_nseq; a+=ns[0],e++)
	          {
		  for ( b=0; b<ns[0]; b++)
		      {
		      seq=l_s[0][b];
		      A->order[a+b][0]=B->order[seq][0];
		      A->order[a+b][1]=B->order[seq][1];
		      for ( c=0; c<start_ml[e];c++)A->order[a+b][1]+=!is_gap(B->seq_al[seq][c]);
		      sprintf ( A->name[a+b], "Repeat_%d", a+b);
		      }
		  }
	
	      free_aln(B);
	      A->nseq=tot_nseq;
	      A->len_aln=strlen ( A->seq_al[0]);
	      
	      }
	   return A;
	   }
Alignment * domain_seq2domain (Constraint_list *CL,int scale,int gop,int gep,Alignment *SEQ_DOMAIN, Alignment *TARGET) 
	   {	
	   static Alignment *A;
	   int *n_groups;
	   int **group_list;
	   int a,b,c;
	   
	   A=copy_aln (TARGET, A);
	   A=stack_aln( A, SEQ_DOMAIN);

	   
	   n_groups=(int*)vcalloc ( 2, sizeof (int));
	   group_list=declare_int (2, A->nseq);
	   
	   n_groups[0]=TARGET->nseq;
	   n_groups[1]=SEQ_DOMAIN->nseq;
	   for (c=0, a=0; a< 2; a++)
	       {
	       for (b=0; b< n_groups[a]; b++, c++)
	           {
		   group_list[a][b]=c;
		   }
	       }	   
	   A->score_aln=domain_pair_wise (A, n_groups, group_list,CL);
	   
	   SEQ_DOMAIN=copy_aln (A, SEQ_DOMAIN);
	   vfree (n_groups);
	   free_int (group_list,-1);
	   return SEQ_DOMAIN;
	   }
	   

