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




      
int evaluate_moca_domain ( Alignment *A, Constraint_list *CL)
        {
	  /*
	    function documentation: start
	    int evaluate_moca_domain ( Alignment *A, Constraint_list *CL)

	    This function evaluates a multiple local alignment
	    If    the alignmnent is to be accepted, return score
	    Else  return UNDEFINED
	    
	    function documentation: end
	  */


	  int score=0;
	  int start, end, a, b;
	  Alignment *B=NULL;
	  char alp[200];
	 

	  score=UNDEFINED;
	  end=0;
	  start=0;
	  
	  sprintf ( alp, "acefghiklmnpqrstuvwy");
	  
	  if ( A->len_aln>0)
	    {
	      score=(int)(output_maln_pval ( "/dev/null", A)*-100);
	      return score;
	    }
	  else
	    return 0;
	  
		
			
	  
	  while ((end+1)!=A->len_aln)
	    {
	      end=get_nol_aln_border (A,start,GO_RIGHT);
	      if ( end==start)break;
	      fprintf ( stderr, "\n**%d %d (%d)",start, end, A->len_aln); 
	      B=copy_aln (A, B);
	      B=extract_aln (B,start,end);
	      for (a=0; a<B->nseq; a++)
		for ( b=0; b<B->len_aln; b++)
		  if ( is_gap (B->seq_al[a][b]))B->seq_al[a][b]=alp[(int)rand()%(strlen (alp))];
	      
	      
	      start=end;
	      fprintf ( stderr, "==>%d",(int)(output_maln_pval ( "/dev/null", B)*-100) );
	      if ( score==UNDEFINED)score=(int)(output_maln_pval ( "/dev/null", B)*-100);
	      else
		score=MAX(score,(int)(output_maln_pval ( "/dev/null", B)*-100));
	      
	      
	    }
	  free_aln (B);
	  return score;
	}


int moca_slow_get_dp_cost ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL)
        {
	  int s;
	  
	  s=slow_get_dp_cost ( A, pos1, ns1, list1, col1, pos2, ns2, list2, col2, CL);
	  
	  
	  if ( s==UNDEFINED)return UNDEFINED;
	  else return s+(CL->moca)->moca_scale;

	}
int moca_evaluate_matrix_score ( Constraint_list *CL, int s1, int r1, int s2, int r2)
{
       /*
	    function documentation: start
	    int moca_residue_pair_extended_list ( Constraint_list *CL, int s1, int r1, int s2, int r2 )

	    THIS FUNCTION RETURNS THE EXTENDED SCORE OF A PAIR OF RESIDUES
	    it is meant to work with local aln pair_wise routines, by using (CL->moca)->forbiden_residues
	    a constant value is substracted from the extended score.
    
	    This function is meant toi be used with omain_dp, therefore, it allows the match of identical residues.

	    function documentation: end
	*/
	
	if (unpack_seq_residues ( &s1, &r1, &s2, &r2, CL->packed_seq_lu)==UNDEFINED)return UNDEFINED;
	else if ( (CL->moca)->forbiden_residues && ((CL->moca)->forbiden_residues[s1][r1]==UNDEFINED ||(CL->moca)->forbiden_residues[s2][r2]==UNDEFINED))return UNDEFINED; 
	else if ( s1==s2 && r1 == r2) return UNDEFINED;
	else return  evaluate_matrix_score(CL, s1, r1, s2, r2);	
	}
  



int moca_residue_pair_extended_list ( Constraint_list *CL, int s1, int r1, int s2, int r2 )
        {
	  /*
	    function documentation: start
	    int moca_residue_pair_extended_list ( Constraint_list *CL, int s1, int r1, int s2, int r2 )

	    THIS FUNCTION RETURNS THE EXTENDED SCORE OF A PAIR OF RESIDUES
	    it is meant to work with local aln pair_wise routines, by using (CL->moca)->forbiden_residues
	    a constant value is substracted from the extended score.
    
	    This function is meant toi be used with omain_dp, therefore, it allows the match of identical residues.

	    function documentation: end
	*/

	if (unpack_seq_residues ( &s1, &r1, &s2, &r2, CL->packed_seq_lu)==UNDEFINED)return UNDEFINED;
	else if ( (CL->moca)->forbiden_residues && ((CL->moca)->forbiden_residues[s1][r1]==UNDEFINED ||(CL->moca)->forbiden_residues[s2][r2]==UNDEFINED))return UNDEFINED; 
	else if ( s1==s2 && r1 == r2) return UNDEFINED;
	else return residue_pair_extended_list (CL, s1, r1, s2, r2);	
	
	
	}

int **cache_cl_with_moca_domain (Alignment *A, Constraint_list *CL)
        {
	  /*
	    function documentation: start
	    int **cache_cl_with_moca_domain (Alignment *A, Constraint_list *CL)
	    
	    Read a multiple alignmnent
	    Given all the residues (CL->S)->seq[x][y] contained in the maln
	    Set (CL->moca)->forbiden_residues[x][y] to UNDEFINED
	    return (CL->moca)->forbiden_residues

	    WARNING
            You must make sure that the evalation strategy uses (CL->moca)->forbiden_residues
	    (CL->moca)->forbiden_residues[0][1]->first residue(1) of First sequence(0) 
	    function documentation: end
	  */

	  int **pos;
	  int a, b;
	  
	  pos=aln2pos_simple(A, A->nseq);
	  
	  if ( !(CL->moca)->forbiden_residues)(CL->moca)->forbiden_residues=declare_int ((CL->S)->nseq, strlen ((CL->S)->seq[0])+1);
	  
	  for ( a=0; a<A->nseq;a++)
	    for ( b=0; b< A->len_aln; b++)
	      (CL->moca)->forbiden_residues[A->order[a][0]][pos[a][b]]=UNDEFINED;

	  free_int (pos, -1);
	  return (CL->moca)->forbiden_residues;
	}
Alignment *make_moca_nol_aln ( Alignment *A, Constraint_list *CL)
{
  
  return A;
}

/*********************************************************************************************/
/*                                                                                           */
/*         DOMAIN Z SCORE EVALUATION                                                         */
/*                                                                                           */
/*********************************************************************************************/

int evaluate_domain_aln_z_score (Alignment *A, int start, int end,Constraint_list *CL, char *alphabet)
    {
    int a;
    static Alignment *B;
    double score, ref_score;
    double N_EVAL=1000;
    double sum=0, sum2=0;


    if ( A==NULL || A->nseq==0 || A->len_aln==0)return 0;
    ref_score=(double)evaluate_domain_aln (A,start, end,CL);
    for (sum=0, sum2=0,a=0;a<N_EVAL; a++)
         {
	 B=make_random_aln ( B, A->nseq, end-start, alphabet);	 	
	 score=(double)evaluate_domain_aln (B,0,B->len_aln,CL);
	 sum+=score;
	 sum2+=score*score;
	 }
     score=(return_z_score(ref_score, sum, sum2, N_EVAL)*100)/A->len_aln;
     
     return(int) score;
     }

int evaluate_domain_aln  ( Alignment *A, int start, int end,Constraint_list *CL)
     {
     int a, b, c;
     int score, c1, c2;
     static int **mat;

     /*
       function documentation: start

       This function uses a pam250 to evaluate the sum of pairs score of A, 
       between position start(included) to position end (exluded), 
       
       the numbering starts 0
       function documentation: end
     */

     if ( !mat)mat=read_matrice ( "pam250mt");
     
     for ( c=start, score=0; c<end; c++)
         {
	 for ( a=0; a< A->nseq-1; a++)
	     for ( b=a+1; b< A->nseq; b++)
	         {
		 c1=tolower(A->seq_al[a][c]);
		 c2=tolower(A->seq_al[b][c]);
		 
		 if ( !is_gap (c1) && !is_gap(c2))score+=mat[c1-'A'][c2-'A'];
		 }
	 }
     return score;
     }

int unpack_seq_residues ( int *s1, int *r1, int *s2, int *r2, int **packed_seq_lu)
        {
	  /* Given a series of sequences concatenated (packed), and the coordinates of two residues
	     This function translates the coordinates into the real ones and allows evaluation
	     Note for this function residues go from [1->N], sequences from [0->N[
	     This is true for in and out comming residues number
	     NOTE: The  sequence cannot be guessed when the residues r1 or r2 are GAPS, therefore UNDEFINED is returned
	     NOTE: Concatenated sequences are separated with X, such residues cause an UNDEFINED to be returned
	  */

	  if ( packed_seq_lu==NULL)return 1;	  
	  else if ( s1[0]!=s2[0])return 1;
	  else if (  r1[0]<=0 || r2[0]<=0)return UNDEFINED;
	  else if (  packed_seq_lu[r1[0]][0]==UNDEFINED || packed_seq_lu[r2[0]][0]==UNDEFINED)return UNDEFINED;
	  else
	    {
	      s1[0]=packed_seq_lu[r1[0]][0];
	      r1[0]=packed_seq_lu[r1[0]][1];
		 	
	      s2[0]=packed_seq_lu[r2[0]][0];
	      r2[0]=packed_seq_lu[r2[0]][1];
	    }
	return 1;
	}

Alignment * unpack_seq_aln ( Alignment *A,Constraint_list *CL)
        {
	  int a, b, r_seq, r_start, r_len;
	  

	  if (!CL->packed_seq_lu) return A;
	  
	  for (a=0; a< A->nseq; a++)
	    {
	      r_seq  =CL->packed_seq_lu[A->order[a][1]+1][0];
	      r_start=CL->packed_seq_lu[A->order[a][1]+1][1];
	      
	      A->order[a][0]=r_seq;
	      A->order[a][1]=r_start-1;

	      for ( r_len=0,b=0; b< A->len_aln; b++)r_len+=!is_gap(A->seq_al[a][b]);
	      sprintf ( A->name[a],"%s_%d_%d", (A->S)->name[r_seq], r_start, r_start+r_len-1);
	    }

	return A;
	}


