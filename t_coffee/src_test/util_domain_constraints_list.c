#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
/*********************************************************************/
/*                                                                   */
/*                         MASKING LIST FUNCTIONS                    */
/*                                                                   */
/*                                                                   */
/*********************************************************************/	
Constraint_list * mask_list_with_aln (Alignment *A,int start, int len,Constraint_list *CL, int new_value)
     {
     int a, c,d;
     static int *entry;
     int s1, s2, r1, r2;
     int **pos;
     int **cache;
     int max_nseq;
     int max_len;

     
   
     if ( A==NULL || A->len_aln==0 || A->nseq<=1)return CL;
     if ( entry==NULL) entry=vcalloc (CL->entry_len , CL->el_size);

     cache=declare_int (return_max_int (A->order,read_size_int ( A->order,sizeof (int*)),0)+1,return_max_int (A->order,read_size_int ( A->order,sizeof (int)),1)+A->len_aln+1);

     max_nseq=return_max_int (A->order,read_size_int ( A->order,sizeof (int*)),0)+1;
     max_len =return_max_int (A->order,read_size_int ( A->order,sizeof (int*)),1)+A->len_aln+1;

     pos=aln2pos_simple(A, A->nseq);

     for (a=0; a< A->nseq; a++)
	 for (c=start; c< (start+len); c++)
	     if ( pos[a][c]>0)cache[A->order[a][0]][pos[a][c]]=1;
     
	 
     if (!CL->M)
	  {
	  for ( d=0; d<CL->ne; d++)
		     {
		     s1=vread_clist(CL, d, SEQ1);
		     s2=vread_clist(CL, d, SEQ2);
		     r1=vread_clist(CL, d, R1);
		     r2=vread_clist(CL, d, R2);
		     if ( s1>=max_nseq);
		     else if (r1>=max_len);
		     else if ( cache[s1][r1]==1)
			 {
			 mask_entry(CL,d,new_value);
			 }
		     
		     if ( s2>=max_nseq);		  
		     else if (r2>=max_len);
		     else if ( cache[s2][r2]==1)
		              {
			      mask_entry(CL,d,new_value);
			      }
		     }
    
	  sort_constraint_list_inv (CL,0, CL->ne);    
	  sort_constraint_list (CL,0, CL->ne);         
	  free_int ( cache, -1);
	  }
     else if ( CL->M)
          {
	 
	  for (a=0; a< A->nseq; a++)
	      for (c=start; c< (start+len); c++)
		  if ( pos[a][c]>0)
		      {
		      vwrite_clist(CL,30+A->order[a][0],pos[a][c],UNDEFINED);
		      }
	  }
     free_int (pos, -1);
     return CL;
     }
Constraint_list* mask_list_with_aln_pair (Alignment *A,int start, int len ,Constraint_list *CL,int new_value)
     {
     int a, b, p;
     int *entry;
   

     int l1, l2, r1, r2, s1, s2;
     int x, y;

     if ( A==NULL || A->len_aln==0 || A->nseq<=1)return CL;

     if (CL->M)
        {
	fprintf ( stderr, "\nERROR: AA matrix cannot be masked with  mask_list_with_aln_pair");
	myexit (EXIT_SUCCESS);
	}
     
    
     entry=vcalloc (CL->entry_len, sizeof (int));
     
     for ( a=0; a< A->nseq-1; a++)
         {
	 for (l1=A->order[a][1]+1,p=0    ; p<start        ; p++)l1+=!is_gap(A->seq_al[a][p]);
	 for (r1=l1-1            ,p=start; p<(start+len)  ; p++)r1+=!is_gap(A->seq_al[a][p]);
	 s1=A->order[a][0];

	 for ( b=a+1; b< A->nseq; b++)
	     {
	     for (l2=A->order[b][1]+1,p=0    ; p<start; p++)l2+=!is_gap(A->seq_al[b][p]);
	     for (r2=l2-1            ,p=start; p<(start+len)  ; p++)r2+=!is_gap(A->seq_al[b][p]);
	     s2=A->order[b][0];
	     
	    
	     for ( x=l1; x<=r1; x++)
		 {
		 
		 for ( y=l2; y<=r2; y++)
		     {
		     
		     set_int(entry,4,x,R1,y,R2,s1,SEQ1,s2,SEQ2); 
		     if ( (main_search_in_list_constraint ( entry,&p,4,CL))!=NULL)
			  {			  
			   mask_entry(CL,p,new_value);
			  }
		     set_int(entry,4,y,R1,x,R2,s2,SEQ1,s1,SEQ2); 
		     if ( (main_search_in_list_constraint ( entry,&p,4,CL))!=NULL)
			  {
			  
			  mask_entry(CL,p,new_value); 
			  }
		     }
		 }
	     }
	 }


     vfree(entry);
     
     return CL;
     }

Constraint_list *mask_entry( Constraint_list *CL, int p, int new_value)
     {
     vwrite_clist(CL, p, WE,   new_value);
     vwrite_clist(CL, p, CONS, new_value);
     vwrite_clist(CL, p, MISC, new_value);
     return CL;
     }
/*********************************************************************/
/*                                                                   */
/*                         SEQUENCE CONCATENATION                    */
/*                                                                   */
/*                                                                   */
/*********************************************************************/	
Constraint_list *prepare_list_and_seq4sw(Constraint_list *I, int n_seq, char **seq_name)
        {
	int a, b;
	int len,l;
	char **long_seq=NULL;
	int  **translation;
	char **name;
	int s1, s2;
	Constraint_list *Out=NULL;
	Sequence *S_Out;

	translation=declare_int ( (I->S)->nseq,2);
	name=declare_char (1,STRING);
	long_seq=declare_char(1, STRING);

	for (len=0,a=0; a< (I->S)->nseq; a++)
	    {
	    if((b=name_is_in_list ((I->S)->name[a],seq_name, n_seq, 100))!=-1)
	       {
	       l=strlen((I->S)->seq[a])+1;
	       long_seq[0]=vrealloc(long_seq[0],(len+l+1)*sizeof(char));
	       long_seq[0]=strcat(long_seq[0], (I->S)->seq[a]);
	       long_seq[0]=strcat(long_seq[0], "O");
	       
	       translation[a][0]=b;
	       translation[a][1]=len;
	       len+=l;
	       }
	    else translation[a][0]=-1;
	    }

	long_seq[0][len-1]='\0';
	len--;
	sprintf ( name[0], "concatenat");
	S_Out=fill_sequence_struc(1, long_seq, name);	
	free_char(name, -1);
	free_char(long_seq, -1);

	
	if (!I->M)
	   {
	   if ( I->fp)    Out=declare_constraint_list(S_Out, NULL, NULL, 0, vtmpfile(),NULL);
	   else if ( I->L)Out=declare_constraint_list(S_Out, NULL, NULL, 0, NULL     ,NULL);

	   for (a=0; a<I->ne; a++)
	       {
	       s1=vread_clist(I,a,SEQ1);
	       s2=vread_clist(I,a,SEQ2);

	       if ( translation[s1][0]!=-1 &&  translation[s2][0]!=-1)
		   Out=add_list_entry2list(Out, Out->entry_len, SEQ1, 0, SEQ2, 0, R1,vread_clist(I,a,R1)+translation[s1][1], R2,vread_clist(I,a,R2)+translation[s2][1], WE,vread_clist(I,a,WE),CONS, vread_clist(I,a,CONS), MISC, vread_clist(I,a,MISC));   
	       }

	   for ( a=0; a<(I->S)->nseq; a++)
	       {
	       if (translation[a][0]!=-1 && translation[a][1]!=0)
	          {
		  for (b=1; b<=len; b++)
	              {
		      add_list_entry2list(Out,Out->entry_len, SEQ1, 0, SEQ2, 0, R1, translation[a][1], R2, b, WE, UNDEFINED, CONS, 0,MISC, vread_clist(I,a,MISC));
		      }
		  }
	       }
	   sort_constraint_list (Out, 0, Out->ne);
	   }
	else if (I->M)
	   { 
	   Out=declare_constraint_list(S_Out, NULL, NULL, 0, vtmpfile(),I->M);
	   vfree((Out->M)[30]);
	   (Out->M)[30]=vcalloc ( len+1, sizeof (int));
	   for ( a=0; a<len; a++)
	       {
	       if ( long_seq[0][a]=='O')vwrite_clist(Out, 30, a+1,UNDEFINED);
	       }
	   Out->ne=SIZEOF_AA_MAT;
	   }
	free_int (translation,-1);
	return Out;
	}
/*********************************************************************/
/*                                                                   */
/*                         MISCEANELLOUS                             */
/*                                                                   */
/*                                                                   */
/*********************************************************************/	
int ** get_undefined_list (Constraint_list *CL)
         {
	 int **list;
	 int a;
	 CLIST_TYPE x;
	 
	 list=declare_int ( (CL->S)->nseq+1, (CL->S)->max_len+1);

	 for ( a=0; a< CL->ne; a++)
	     {
	     x=vread_clist(CL, a, WE);
	     list[vread_clist(CL, a, SEQ1)][vread_clist(CL, a, R1)]=(x==UNDEFINED);
	     list[vread_clist(CL, a, SEQ2)][vread_clist(CL, a, R2)]=(x==UNDEFINED);
	     }
	 return list;
	 }
int      is_never_undefined (Constraint_list *CL,int r)
         {
	 int a;
	 for ( a=0; a< CL->ne; a++)
	     {
	     if ( (vread_clist(CL,a,R1)==r || vread_clist(CL,a,R2)==r) && vread_clist(CL,a,WE)==UNDEFINED)return 0;
	     }
	 return 1;
	 }

int* do_analyse_list ( Constraint_list *CL)
        {
	int **seq_score;
	int **seq_score2;
	int  *pos_L;
	int a;
	int n_res;


	int n_it=4;

	double sum, sum2, tot;
	int field=2;
	double z;
	int r1, r2;
	int max_we=0;



	fprintf ( stderr, "\nDO ANALYSE");
	
	n_res=(CL->S)->max_len;
	
	pos_L =vcalloc     (n_res+2, sizeof (int)); 
	pos_L++;
	pos_L[-1]=n_res;
	seq_score=declare_int (n_res+1,n_it+1);
	seq_score2=declare_int(n_res+1,n_it+1);
	


	for ( a=0; a< CL->ne; a++)
	    {
	    r1=vread_clist(CL, a, R1);
	    r2=vread_clist(CL, a, R2);
	   
	    seq_score[r1][0]+=vread_clist(CL, a, WE);
	    seq_score[r2][0]+=vread_clist(CL, a, WE);
	    
	    
	    seq_score[r1][1]+=vread_clist(CL, a, CONS);
	    seq_score[r2][1]+=vread_clist(CL, a, CONS);
	    
	    
	    seq_score[r1][2]+=vread_clist(CL, a, MISC);
	    seq_score[r2][2]+=vread_clist(CL, a, MISC);
	    
	    }
	for ( a=1; a<=n_res; a++)max_we=MAX(seq_score[a][0],max_we);
	for ( a=1; a<=n_res; a++)
	    {
	    if ( a!=n_res && seq_score[a][0]> seq_score[a-1][0]) fprintf ( stderr, "\n%4d %s", a, num2plot(seq_score[a][0],max_we,40));
	    else fprintf ( stderr, "\n");
	    }
	    
	for ( a=0; a< CL->ne; a++)if ( vread_clist(CL, a, MISC)>1000)fprintf ( stderr, "\n%4d %4d %4d", vread_clist(CL, a,R1), vread_clist(CL, a, R2), vread_clist(CL, a, MISC));
	myexit (EXIT_SUCCESS);

	for ( a=0; a<CL->ne; a++)
	    {
	    if ( vread_clist(CL, a, WE)!=UNDEFINED &&vread_clist(CL, a, MISC)>=1 )
		{
		seq_score[vread_clist(CL, a, R1)][0]+=vread_clist(CL, a, CONS);
		seq_score[vread_clist(CL, a, R1)][1]+=vread_clist(CL, a, WE);
		seq_score[vread_clist(CL, a, R1)][2]+=vread_clist(CL, a, MISC);

		seq_score[vread_clist(CL, a, R2)][0]+=vread_clist(CL, a, CONS);
		seq_score[vread_clist(CL, a, R2)][1]+=vread_clist(CL, a, WE);
		seq_score[vread_clist(CL, a, R2)][2]+=vread_clist(CL, a, MISC);		    		
		}
	    }
	
	for (a=1, tot=0,sum=0, sum2=0; a<= n_res  ; a++)
	    {
	    if ( seq_score[a][2]>0)
		{
		sum +=seq_score[a][field];
		sum2+=seq_score[a][field]*seq_score[a][field];
		tot++;
		}
	    }
	
	fprintf ( stderr, "\n");
	


	for (a=1; a<= n_res  ; a++)
	    {
	    z=return_z_score (seq_score[a][field],sum, sum2, tot );	
	    if ( seq_score[a][2]>0)
		{
		
		pos_L[a]=(int)(z*10);
		pos_L[0]++;
		
		}
	    }

	
	free_int (seq_score2,-1);
	free_int (seq_score,-1);
	return pos_L;
	}
