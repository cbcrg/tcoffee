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


int commonsextet( int *table, int *pointt );
void makecompositiontable( int *table, int *pointt );
int *code_seq (char *seq, char *type);
int * makepointtable( int *pointt, int *n, int ktup );

static int tsize;

/**
* calculates the number of common tuples
*/
int commonsextet( int *table, int *pointt )
{
	int value = 0;
	int tmp;
	int point;
	static int *memo = NULL;
	static int *ct = NULL;
	static int *cp;

	if( !memo )
	{
		memo =(int*) vcalloc( tsize+1, sizeof( int ) );
		ct =(int*) vcalloc( tsize+1, sizeof( int ) );
	}

	cp = ct;
	while( ( point = *pointt++ ) != END_ARRAY )
	{
	  tmp = memo[point]++;
	  if( tmp < table[point] )
	    value++;
	  if( tmp == 0 )
	    {
	      *cp++ = point;
	    }
	}
	*cp = END_ARRAY;

	cp =  ct;
	while( *cp != END_ARRAY )
		memo[*cp++] = 0;

	return( value );
}

/**
*	calculates how many of each tuple exist
*/
void makecompositiontable( int *table, int *pointt )
{
	int point;

	while( ( point = *pointt++ ) != END_ARRAY )
	  {
	    table[point]++;
	  }
}

int *code_seq (char *seq, char *type)
{
  static int *code;
  static int *aa, ng;
  int a, b, l;


  if (!aa)
    {
      char **gl;
      if ( strm (type, "DNA") || strm (type, "RNA"))
	{
	  gl=declare_char (4,5);
	  sprintf ( gl[ng++], "Aa");
	  sprintf ( gl[ng++], "Gg");
	  sprintf ( gl[ng++], "TtUu");
	  sprintf ( gl[ng++], "Cc");
	}
      else
	{

	  gl=make_group_aa ( &ng, "mafft");
	}
      aa=(int*)vcalloc ( 256, sizeof (int));
      for ( a=0; a<ng; a++)
	{
	  for ( b=0; b< strlen (gl[a]); b++)
	    {
	      aa[(int)gl[a][b]]=a;
	    }
	}
      free_char (gl, -1);
    }


  l=strlen (seq);

  if ( code) code--;

  if ( !code || read_array_size (code, sizeof (int))<(l+2))
    {
      vfree (code);
      code=(int*)vcalloc (l+2, sizeof (int));
    }
  code[0]=ng;
  code++;
  for (a=0; a<l; a++)
    {
      code[a]=aa[(int)seq[a]];
    }

  code[a]=END_ARRAY;
  return code;
}


int * makepointtable( int *pointt, int *n, int ktup )
{
  int point, a, ng;
  register int *p;
  static int *prod;

  ng=n[-1];

  if (!prod)
    {
      prod=(int*)vcalloc ( ktup, sizeof (int));
      for ( a=0; a<ktup; a++)
	{
	  prod[ktup-a-1]=(int)pow(n[-1],a);
	}
    }
  p = n;

  for (point=0,a=0; a<ktup; a++)
    {
      point+= *n++ *prod[a];
    }

  *pointt++ = point;

  while( *n != END_ARRAY )
    {
      point -= *p++ * prod[0];
      point *= ng;
      point += *n++;
      *pointt++ = point;
    }
  *pointt = END_ARRAY;
  return pointt;
}


int ** ktup_dist_mat ( char **seq, int nseq, int ktup, char *type)
{
  //Adapted from MAFFT 5: fast ktup
  int **pointt,*code=NULL, **pscore;
  int i, l, j, minl;
  double **mtx, score0;


  if (!seq || nseq==0)return NULL;
  for (minl=strlen(seq[0]),l=0,i=0;i<nseq; i++)
    {
      int len;
      len=strlen (seq[i]);
      minl=MIN(minl, len);
      l=MAX(l,len);
    }
  ktup=MIN(minl, ktup);
  pointt=declare_int (nseq, l+1);
  mtx=declare_double (nseq, nseq);
  pscore=declare_int ( nseq, nseq);

  for( i=0; i<nseq; i++ )
  {
      makepointtable( pointt[i], code=code_seq (seq[i], type),ktup);
  }
  tsize=(int)pow(code[-1], ktup);

  for ( i=0; i<nseq; i++)
    {
      int *table1;
      table1=(int*)vcalloc ( tsize,sizeof (int));
      makecompositiontable( table1, pointt[i]);
      for (j=i; j<nseq; j++)
	{
	  mtx[i][j] = commonsextet( table1, pointt[j] );
	}
      vfree (table1);
    }
  for( i=0; i<nseq; i++ )
    {
      score0 = mtx[i][i];
      for( j=0; j<nseq; j++ )
	pscore[i][j] = (int)( ( score0 - mtx[MIN(i,j)][MAX(i,j)] ) / score0 * 3 * 10.0 + 0.5 );
    }
  for( i=0; i<nseq-1; i++ )
    for( j=i+1; j<nseq; j++ )
      {
	pscore[i][j] = pscore[j][i]=100-MIN( pscore[i][j], pscore[j][i] );
      }
    return pscore;
}


int ** evaluate_diagonals_with_ktup_1 ( Alignment *A, int *ns, int **l_s, Constraint_list *CL,int maximise,int n_groups, char **group_list, int ktup);
int ** evaluate_diagonals_with_ktup_2 ( Alignment *A, int *ns, int **l_s, Constraint_list *CL,int maximise,int n_groups, char **group_list, int ktup);


int ** evaluate_diagonals_for_two_sequences ( char *seq1, char *seq2,int maximise,Constraint_list *CL,int ktup)
       {

       static int ng;
       static char **gl;
       static int *ns, **l_s;
       Alignment *A;
       int **diag;
       int in_cl;
       char *type;

       if (!CL)
	    {
	      in_cl=0;

	      CL=(Constraint_list*)vcalloc ( 1, sizeof (Constraint_list));
	      CL->maximise=1;
	      sprintf ( CL->matrix_for_aa_group, "vasiliky");
	      CL->M=read_matrice ("blosum62mt");
	      CL->evaluate_residue_pair=evaluate_cdna_matrix_score;
	      CL->get_dp_cost=slow_get_dp_cost;
	      type=get_string_type(seq1);

	      if ( strm (type, "CDNA"))
		   CL->evaluate_residue_pair= evaluate_matrix_score;
	      else if (  strm(type, "PROTEIN"))
		   CL->evaluate_residue_pair=evaluate_matrix_score;
	      else if (  strm (type, "DNA") || strm (type, "RNA"))
		   CL->evaluate_residue_pair= evaluate_matrix_score;
	      vfree(type);
	    }
       else
	    {
	      in_cl=1;
	    }




       if ( !gl)
	 {
	   gl=make_group_aa (&ng, CL->matrix_for_aa_group);
	   ns=(int*)vcalloc (2, sizeof (int));
	   ns[0]=ns[1]=1;
	   l_s=declare_int (2, 2);
	   l_s[0][0]=0;
	   l_s[1][0]=1;
	 }


       A=strings2aln (2, "A",seq1,"B", seq2);
       ungap(A->seq_al[0]);
       ungap(A->seq_al[1]);

       CL->S=A->S;

       diag=evaluate_diagonals ( A,ns, l_s, CL,maximise, ng, gl, ktup);
       free_sequence (A->S, (A->S)->nseq);
       free_aln (A);
       if (!in_cl)
	 {
	  free_int (CL->M, -1);
	  vfree (CL);
	 }


       return diag;
       }


int ** evaluate_diagonals ( Alignment *A, int *ns, int **l_s, Constraint_list *CL,int maximise,int n_groups, char **group_list, int ktup)
        {
	int **tot_diag;



	if      ( CL->residue_index)
	  {
	  tot_diag=evaluate_diagonals_with_clist ( A, ns, l_s, CL, maximise,n_groups,group_list, ktup);
	  }
	else if ( CL->use_fragments)
	    {

	      tot_diag=evaluate_segments_with_ktup ( A, ns, l_s, CL, maximise,n_groups,group_list, ktup);
	    }
	else
	  {

	    tot_diag=evaluate_diagonals_with_ktup ( A, ns, l_s, CL, maximise,n_groups,group_list, ktup);
	  }

	return tot_diag;
	}
int ** evaluate_segments_with_ktup ( Alignment *A, int *ns, int **l_s, Constraint_list *CL,int maximise,int n_groups, char **group_list, int ktup)
    {
   /*
    Reads in an alignmnet A, with two groups of sequences marked.
    1-Turn each group into a conscensus, using the group list identifier.
               -if the group list is left empty original symbols are used
    2-hash groupc the two sequences
    3-score each diagonal, sort the list and return it (diag_list)
   */

    char *seq1, *seq2, *alphabet=NULL;
    int a,b,l1, l2, n_ktup,pos_ktup1, pos_ktup2, **pos;
    int *hasched_seq1, *hasched_seq2,*lu_seq1,*lu_seq2;
    int n_diag, **diag, current_diag, **dot_list, n_dots, cost;
    int l,delta_diag, delta_res;


    pos=aln2pos_simple ( A,-1, ns, l_s);
    seq1=aln2cons_seq (A, ns[0], l_s[0], n_groups, group_list);
    seq2=aln2cons_seq (A, ns[1], l_s[1], n_groups, group_list);



    alphabet=get_alphabet (seq1,alphabet);
    alphabet=get_alphabet (seq2,alphabet);



    l1=strlen ( seq1);
    l2=strlen ( seq2);

    n_diag=l1+l2-1;
    diag=declare_int ( n_diag+2, 3);
    n_ktup=(int)pow ( (double)alphabet[0]+1, (double)ktup);

    hasch_seq(seq1, &hasched_seq1, &lu_seq1,ktup, alphabet);
    hasch_seq(seq2, &hasched_seq2, &lu_seq2,ktup, alphabet);



    /*EVALUATE THE DIAGONALS*/
    for ( a=0; a<= n_diag; a++)diag[a][0]=a;
    for ( n_dots=0,a=1; a<= n_ktup; a++)
        {
	    pos_ktup1=lu_seq1[a];
	    while (TRUE)
	          {
		  if (!pos_ktup1)break;
		  pos_ktup2=lu_seq2[a];
		  while (pos_ktup2)
		            {
			    n_dots++;
			    pos_ktup2=hasched_seq2[pos_ktup2];
			    }
		  pos_ktup1=hasched_seq1[pos_ktup1];
		  }
	}

    if ( n_dots==0)
       {
	    vfree (seq1);
	    vfree (seq2);
	    vfree (alphabet);
	    vfree (hasched_seq1);
	    vfree (hasched_seq2);
	    vfree (lu_seq1);
	    vfree (lu_seq2);
	    free_int (diag, -1);
	   return evaluate_segments_with_ktup (A,ns,l_s,CL,maximise,n_groups, group_list,ktup-1);
       }

    dot_list=declare_int ( n_dots,3);

    for ( n_dots=0,a=1; a<= n_ktup; a++)
        {
	    pos_ktup1=lu_seq1[a];
	    while (TRUE)
	          {
		  if (!pos_ktup1)break;
		  pos_ktup2=lu_seq2[a];
		  while (pos_ktup2)
		            {
			    current_diag=(pos_ktup2-pos_ktup1+l1);
			    dot_list[n_dots][0]=current_diag;
			    dot_list[n_dots][1]=pos_ktup1;
			    dot_list[n_dots][2]=pos_ktup2;
			    pos_ktup2=hasched_seq2[pos_ktup2];
			    n_dots++;
			    }
		  pos_ktup1=hasched_seq1[pos_ktup1];
		  }
	}



    hsort_list_array ((void **)dot_list, n_dots, sizeof (int), 3, 0, 3);
    current_diag= (int)dot_list[0][0];

    for ( b=0; b< ktup; b++)diag[current_diag][2]+=(CL->get_dp_cost) ( A, pos, ns[0], l_s[0], dot_list[0][1]+b-1, pos,ns[1], l_s[1], dot_list[0][2]+b-1, CL);


    for ( l=0,a=1; a< n_dots; a++)
        {

	    delta_diag=dot_list[a][0]-dot_list[a-1][0];
	    delta_res =dot_list[a][1]-dot_list[a-1][1];

	    for ( cost=0, b=0; b< ktup; b++)cost++;

	    /*=(CL->get_dp_cost) ( A, pos, ns[0], l_s[0], dot_list[a][1]+b-1, pos,ns[1], l_s[1], dot_list[a][2]+b-1, CL);*/



	    if (delta_diag!=0 || FABS(delta_res)>5)
	       {

		 l=0;
		 diag[current_diag][1]=best_of_a_b(diag[current_diag][2], diag[current_diag][1], 1);
		 if ( diag[current_diag][2]<0);
		 else diag[current_diag][1]= MAX(diag[current_diag][1],diag[current_diag][2]);
		 diag[current_diag][2]=0;
		 current_diag=dot_list[a][0];
	       }
	    l++;
	    diag[current_diag][2]+=cost;

	}
    diag[current_diag][1]=best_of_a_b(diag[current_diag][2], diag[current_diag][1], 1);
    sort_int (diag+1, 3, 1,0, n_diag-1);


    vfree (seq1);
    vfree (seq2);
    vfree (alphabet);
    vfree (hasched_seq1);
    vfree (hasched_seq2);
    vfree (lu_seq1);
    vfree (lu_seq2);
    free_int (pos, -1);
    free_int (dot_list, -1);
    return diag;
    }





int ** evaluate_diagonals_with_clist ( Alignment *A, int *ns, int **l_s, Constraint_list *CL,int maximise,int n_groups, char **group_list, int ktup)
    {

   /*
    Reads in an alignmnent A, with two groups of sequences marked.
    Weight the diagonals with the values read in the constraint list
   */

    int l1, l2,n_diag, s1, s2, r1=0, r2=0;
    int a, b, c, d;
    int **diag;
    int **code;
    int **pos;
    static int *entry;


    if ( !entry)entry=(int*)vcalloc ( CL->entry_len+1, CL->el_size);
    l1=strlen (A->seq_al[l_s[0][0]]);
    l2=strlen (A->seq_al[l_s[1][0]]);

    n_diag=l1+l2-1;
    diag=declare_int ( n_diag+2, 3);
    for ( a=0; a<= n_diag; a++)diag[a][0]=a;

    A->S=CL->S;
    code=seq2aln_pos (A, ns, l_s);
    pos =aln2pos_simple ( A,-1, ns, l_s);


    for (a=0; a<ns[0]; a++)

        {
	s1=A->order[l_s[0][a]][0];
	for (b=0; b<ns[1]; b++)
	    {
	    s2=A->order[l_s[1][b]][0];
	    for (r1=1; r1<=(A->S)->len[s1]; r1++)
	      {
		int e;
		for (e=1; e<CL->residue_index[s1][r1][0]; e+=ICHUNK)
		  {
		    if (CL->residue_index[s1][r1][e+SEQ2]==s2)
		      {
			r2=CL->residue_index[s1][r1][e+R2];
			diag[(r2-r1+l1)][1]+=(CL->get_dp_cost) ( A, pos, ns[0], l_s[0],r1-1, pos,ns[1], l_s[1], r2-1, CL);
		      }
		  }
	      }
	    }
	}

    sort_int (diag+1, 2, 1,0, n_diag-1);

    free_int (code,-1);
    free_int (pos, -1);
    return diag;
    }

int * flag_diagonals (int l1, int l2, int **sorted_diag, float T, int window)
    {
    int a, b, up, low,current_diag,n_diag;
    int * slopes;
    int *diag_list;
    double mean;
    double sd;
    int use_z_score=1;


    n_diag=l1+l2-1;
    mean=return_mean_int ( sorted_diag, n_diag+1, 1);

    sd  =return_sd_int ( sorted_diag, n_diag+1, 1, (int)mean);

    if ( T==0)
      {
      use_z_score=1;
      T=(((double)sorted_diag[n_diag][1]-mean)/sd)/25;
      }


    diag_list=(int*)vcalloc (l1+l2+1, sizeof (int));
    slopes=(int*)vcalloc ( n_diag+1, sizeof (int));

    for ( a=n_diag; a>0; a--)
            {
	    current_diag=sorted_diag[a][0];


	    if ( !use_z_score && sorted_diag[a][1]>T)
	       {
		   up=MAX(1,current_diag-window);
		   low=MIN(n_diag, current_diag+window);
		   for ( b=up; b<=low; b++)slopes[b]=1;
	       }
	    else if (use_z_score && ((double)sorted_diag[a][1]-mean)/sd>T)
	      {
		up=MAX(1,current_diag-window);
		low=MIN(n_diag, current_diag+window);
		for ( b=up; b<=low; b++)slopes[b]=1;
	      }
	    else break;
	    }

    for ( a=1, b=0; a<=n_diag; a++)
        {
	    b+=slopes[a];
	}

    slopes[1]=1;
    slopes[l1+l2-1]=1;
    slopes[l2]=1;
    for (a=0; a<= (l1+l2-1); a++)
	if ( slopes[a]){diag_list[++diag_list[0]]=a;}

    vfree (slopes);

    return diag_list;
    }
int * extract_N_diag (int l1, int l2, int **sorted_diag, int n_chosen_diag, int window)
    {
    int a, b, up, low,current_diag,n_diag;
    int * slopes;
    int *diag_list;


    n_diag=l1+l2-1;

    diag_list=(int*)vcalloc (l1+l2+1, sizeof (int));
    slopes=(int*)vcalloc ( n_diag+1, sizeof (int));


    for ( a=n_diag; a>0 && a>(n_diag-n_chosen_diag); a--)
            {
	    current_diag=sorted_diag[a][0];
	    up=MAX(1,current_diag-window);
	    low=MIN(n_diag, current_diag+window);

	    for ( b=up; b<=low; b++)slopes[b]=1;
	    }

    /*flag bottom right*/
    up=MAX(1,1-window);low=MIN(n_diag,1+window);
    for ( a=up; a<=low; a++) slopes[a]=1;

    /*flag top left */
    up=MAX(1,(l1+l2-1)-window);low=MIN(n_diag,(l1+l2-1)+window);
    for ( a=up; a<=low; a++) slopes[a]=1;


    /*flag MAIN DIAG SEQ1*/
    up=MAX(1,l1-window);low=MIN(n_diag,l1+window);
    for ( a=up; a<=low; a++) slopes[a]=1;

    /*flag MAIN DIAG SEQ2*/
    up=MAX(1,l2-window);low=MIN(n_diag,l2+window);
    for ( a=up; a<=low; a++) slopes[a]=1;


    for (a=0; a<= (l1+l2-1); a++)
	if ( slopes[a]){diag_list[++diag_list[0]]=a;}

    vfree (slopes);
    return diag_list;
    }




int cfasta_gotoh_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL)
    {
/*TREATMENT OF THE TERMINAL GAP PENALTIES*/
/*TG_MODE=0---> gop and gep*/
/*TG_MODE=1---> ---     gep*/
/*TG_MODE=2---> ---     ---*/


	int maximise;

/*VARIABLES FOR THE MULTIPLE SEQUENCE ALIGNMENT*/
	int **tot_diag;

	int *diag;
	int ktup;
	static int n_groups;
	static char **group_list;
	int score, new_score;
        int n_chosen_diag=20;
        int step;
	int max_n_chosen_diag;
	int l1, l2;
        /********Prepare Penalties******/


	maximise=CL->maximise;
	ktup=CL->ktup;

	/********************************/




	if ( !group_list)
	   {

	       group_list=make_group_aa (&n_groups, CL->matrix_for_aa_group);
	   }

	l1=strlen (A->seq_al[l_s[0][0]]);
	l2=strlen (A->seq_al[l_s[1][0]]);

	if ( !CL->fasta_step)
	    {
	    step=MIN(l1,l2);
	    step=(int) log ((double)MAX(step, 1));
	    step=MAX(step, 20);
	    }
	else
	    {
		step=CL->fasta_step;
	    }


	tot_diag=evaluate_diagonals ( A, ns, l_s, CL, maximise,n_groups,group_list, ktup);


	max_n_chosen_diag=strlen (A->seq_al[l_s[0][0]])+strlen (A->seq_al[l_s[1][0]])-1;
	/*max_n_chosen_diag=(int)log10((double)(l1+l2))*10;*/

	n_chosen_diag+=step;
	n_chosen_diag=MIN(n_chosen_diag, max_n_chosen_diag);


	diag=extract_N_diag (strlen (A->seq_al[l_s[0][0]]),strlen (A->seq_al[l_s[1][0]]), tot_diag, n_chosen_diag, 0);


	score    =make_fasta_gotoh_pair_wise ( A, ns, l_s, CL, diag);

	new_score=0;
	vfree ( diag);


	while (new_score!=score && n_chosen_diag< max_n_chosen_diag    )
	  {


	    score=new_score;

	    ungap_sub_aln ( A, ns[0], l_s[0]);
	    ungap_sub_aln ( A, ns[1], l_s[1]);


	    n_chosen_diag+=step;
	    n_chosen_diag=MIN(n_chosen_diag, max_n_chosen_diag);


	    diag     =extract_N_diag (strlen (A->seq_al[l_s[0][0]]),strlen (A->seq_al[l_s[1][0]]), tot_diag, n_chosen_diag, 0);
	    new_score=make_fasta_gotoh_pair_wise (  A, ns, l_s, CL, diag);

	    vfree ( diag);

	  }

	score=new_score;
	free_int (tot_diag, -1);

	return score;
    }

int fasta_gotoh_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL)
    {
/*TREATMENT OF THE TERMINAL GAP PENALTIES*/
/*TG_MODE=0---> gop and gep*/
/*TG_MODE=1---> ---     gep*/
/*TG_MODE=2---> ---     ---*/


	int maximise;

/*VARIABLES FOR THE MULTIPLE SEQUENCE ALIGNMENT*/
	int **tot_diag;
	int *diag;
	int ktup;
	float diagonal_threshold;
	static int n_groups;
	static char **group_list;
	int score;
        /********Prepare Penalties******/


	maximise=CL->maximise;
	ktup=CL->ktup;
	diagonal_threshold=CL->diagonal_threshold;
	/********************************/



	if ( !group_list)
	   {
	       group_list=make_group_aa (&n_groups, CL->matrix_for_aa_group);
	   }


	tot_diag=evaluate_diagonals ( A, ns, l_s, CL, maximise,n_groups,group_list, ktup);

	if (  !CL->fasta_step)
	  {
	    diag=flag_diagonals (strlen(A->seq_al[l_s[0][0]]),strlen(A->seq_al[l_s[1][0]]), tot_diag,diagonal_threshold,0);
	  }

	else
	  {

	    diag=extract_N_diag (strlen (A->seq_al[l_s[0][0]]),strlen (A->seq_al[l_s[1][0]]), tot_diag,CL->fasta_step,0);

	  }
	score=make_fasta_gotoh_pair_wise ( A, ns, l_s, CL, diag);

	free_int (tot_diag, -1);
	vfree (diag);
	return score;
    }
int very_fast_gotoh_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL)
    {
/*TREATMENT OF THE TERMINAL GAP PENALTIES*/
/*TG_MODE=0---> gop and gep*/
/*TG_MODE=1---> ---     gep*/
/*TG_MODE=2---> ---     ---*/


	int maximise;
/*VARIABLES FOR THE MULTIPLE SEQUENCE ALIGNMENT*/
	int **tot_diag;
	int *diag;
	int ktup;
	static int n_groups;
	static char **group_list;
	int score;
        /********Prepare Penalties******/


	maximise=CL->maximise;
	ktup=CL->ktup;
	/********************************/


	if ( !group_list)
	   {

	       group_list=make_group_aa (&n_groups, CL->matrix_for_aa_group);
	   }

	CL->use_fragments=0;
	tot_diag=evaluate_diagonals ( A, ns, l_s, CL, maximise,n_groups,group_list, ktup);

	/*Note: 20 diagonals. 5 shadows on each side: tunned on Hom39, 2/2/04 */
	diag=extract_N_diag (strlen (A->seq_al[l_s[0][0]]),strlen (A->seq_al[l_s[1][0]]), tot_diag,20,5);
	score=make_fasta_gotoh_pair_wise ( A, ns, l_s, CL, diag);
	free_int (tot_diag, -1);
	vfree (diag);
	return score;
    }
int make_fasta_gotoh_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL, int *diag)
    {
/*TREATMENT OF THE TERMINAL GAP PENALTIES*/
/*TG_MODE=0---> gop and gep*/
/*TG_MODE=1---> ---     gep*/
      /*TG_MODE=2---> ---     ---*/


	int TG_MODE, gop, l_gop, gep,l_gep, maximise;

/*VARIABLES FOR THE MULTIPLE SEQUENCE ALIGNMENT*/
	int a, b,c,k, t;
	int l1, l2,eg, ch, sub,score=0, last_i=0, last_j=0, i, delta_i, j, pos_j, ala, alb, LEN, n_diag, match1, match2;
	int su, in, de, tr;

	int **C, **D, **I, **trace, **pos0, **LD;
	int lenal[2], len;
	char *buffer, *char_buf;
	char **aln, **al;

        /********Prepare Penalties******/
	gop=CL->gop*SCORE_K;
	gep=CL->gep*SCORE_K;
	TG_MODE=CL->TG_MODE;
	maximise=CL->maximise;


	/********************************/


        n_diag=diag[0];



       l1=lenal[0]=strlen (A->seq_al[l_s[0][0]]);
       l2=lenal[1]=strlen (A->seq_al[l_s[1][0]]);

       if ( getenv ("DEBUG_TCOFFEE"))fprintf ( stderr, "\n\tNdiag=%d%%  ", (diag[0]*100)/(l1+l2));

	/*diag:
	  diag[1..n_diag]--> flaged diagonal in order;
	  diag[0]=0--> first diagonal;
	  diag[n_diag+1]=l1+l2-1;
	*/

	/*numeration of the diagonals strats from the bottom right [1...l1+l2-1]*/
	/*sequence s1 is vertical and seq s2 is horizontal*/
	/*D contains the best Deletion  in S2==>comes from diagonal N+1*/
	/*I contains the best insertion in S2=> comes from diagonal N-1*/





       C=declare_int (lenal[0]+lenal[1]+1, n_diag+2);
       D=declare_int (lenal[0]+lenal[1]+1, n_diag+2);
       LD=declare_int (lenal[0]+lenal[1]+1, n_diag+2);
       I=declare_int (lenal[0]+lenal[1]+1, n_diag+2);
       trace=declare_int (lenal[0]+lenal[1]+1, n_diag+2);


       al=declare_char (2,lenal[0]+lenal[1]+lenal[1]+1);

       len= MAX(lenal[0],lenal[1])+1;
       buffer=(char*)vcalloc ( 2*len, sizeof (char));
       char_buf=(char*) vcalloc (2*len, sizeof (char));

       pos0=aln2pos_simple ( A,-1, ns, l_s);
       C[0][0]=0;

       t=(TG_MODE==0)?gop:0;
       for ( j=1; j<= n_diag; j++)
	    {
		l_gop=(TG_MODE==0)?gop:0;
		l_gep=(TG_MODE==2)?0:gep;



		if ( (diag[j]-lenal[0])<0 )
		    {
		    trace[0][j]=UNDEFINED;
		    continue;
		    }
		C[0][j]=(diag[j]-lenal[0])*l_gep +l_gop;
		D[0][j]=(diag[j]-lenal[0])*l_gep +l_gop+gop;
	    }
       D[0][j]=D[0][j-1]+gep;


       t=(TG_MODE==0)?gop:0;
       for ( i=1; i<=lenal[0]; i++)
           {
	        l_gop=(TG_MODE==0)?gop:0;
		l_gep=(TG_MODE==2)?0:gep;

		C[i][0]=C[i][n_diag+1]=t=t+l_gep;
		I[i][0]=D[i][n_diag+1]=t+    gop;

		for ( j=1; j<=n_diag; j++)
		    {
			C[i][j]=C[i][0];
			D[i][j]=I[i][j]=I[i][0];
		    }

		for (eg=0, j=1; j<=n_diag; j++)
		    {

			pos_j=diag[j]-lenal[0]+i;
			if (pos_j<=0 || pos_j>l2 )
			    {
			    trace[i][j]=UNDEFINED;
			    continue;
			    }
			sub=(CL->get_dp_cost) ( A, pos0, ns[0], l_s[0], i-1, pos0, ns[1], l_s[1],pos_j-1, CL );

		    /*1 identify the best insertion in S2:*/
			l_gop=(i==lenal[0])?((TG_MODE==0)?gop:0):gop;
			l_gep=(i==lenal[0])?((TG_MODE==2)?0:gep):gep;
			len=(j==1)?0:(diag[j]-diag[j-1]);
			if ( a_better_than_b(I[i][j-1], C[i][j-1]+l_gop, maximise))eg++;
			else eg=1;
			I[i][j]=best_of_a_b (I[i][j-1], C[i][j-1]+l_gop, maximise)+len*l_gep;

		    /*2 Identify the best deletion in S2*/
			l_gop=(pos_j==lenal[1])?((TG_MODE==0)?gop:0):gop;
			l_gep=(pos_j==lenal[1])?((TG_MODE==2)?0:gep):gep;

			len=(j==n_diag)?0:(diag[j+1]-diag[j]);
			delta_i=((i-len)>0)?(i-len):0;

			if ( a_better_than_b(D[delta_i][j+1],C[delta_i][j+1]+l_gop, maximise)){LD[i][j]=LD[delta_i][j+1]+1;}
			else {LD[i][j]=1;}
			D[i][j]=best_of_a_b (D[delta_i][j+1],C[delta_i][j+1]+l_gop, maximise)+len*l_gep;


			/*Identify the best way*/
			/*
			score=C[i][j]=best_int ( 3, maximise, &fop, I[i][j], C[i-1][j]+sub, D[i][j]);
			fop-=1;
			if ( fop<0)trace[i][j]=fop*eg;
			else if ( fop>0 ) {trace[i][j]=fop*LD[i][j];}
			else if ( fop==0) trace[i][j]=0;
			*/

			su=C[i-1][j]+sub;
			in=I[i][j];
			de=D[i][j];

			/*HERE ("%d %d %d", su, in, de);*/
			if (su>=in && su>=de)
			  {
			    score=su;
			    tr=0;
			  }
			else if (in>=de)
			  {
			    score=in;
			    tr=-eg;
			  }
			else
			  {
			    score=de;
			    tr=LD[i][j];
			  }
			trace[i][j]=tr;
			C[i][j]=score;


			last_i=i;
			last_j=j;
		    }
	    }


       /*
	            [0][Positive]
	             ^     ^
	             |    /
                     |   /
                     |  /
                     | /
                     |/
       [Neg]<-------[*]
	*/


	i=last_i;
	j=last_j;



	ala=alb=0;
	match1=match2=0;
	while (!(match1==l1 && match2==l2))
	      {


		  if ( match1==l1)
		     {
			 len=l2-match2;
			 for ( a=0; a< len; a++)
			     {
			     al[0][ala++]=0;
			     al[1][alb++]=1;
			     match2++;
			     }
			 k=0;
			 break;

			 /*k=-(j-1);*/

		     }
		  else if ( match2==l2)
		     {
			 len=l1-match1;
			 for ( a=0; a< len; a++)
			     {
			     al[0][ala++]=1;
			     al[1][alb++]=0;
			     match1++;
			     }
			 k=0;
			 break;
			 /*k= n_diag-j;*/
		     }
		  else
		      {
			  k=trace[i][j];
		      }


		  if ( k==0)
			     {
				 if ( match2==l2 || match1==l1);
				 else
				    {

				    al[0][ala++]=1;
				    al[1][alb++]=1;
				    i--;
				    match1++;
				    match2++;
				    }
			     }
		  else if ( k>0)
			     {

			     len=diag[j+k]-diag[j];
			     for ( a=0; a<len; a++)
			         {
				     if ( match1==l1)break;
				     al[0][ala++]=1;
				     al[1][alb++]=0;
				     match1++;
				 }
			     i-=len;
			     j+=k;
			     }
		  else if ( k<0)
			     {
			     k*=-1;
			     len=diag[j]-diag[j-k];
			     for ( a=0; a<len; a++)
			         {
				     if ( match2==l2)break;
				     al[0][ala++]=0;
				     al[1][alb++]=1;
				     match2++;
				 }


			     j-=k;
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
		aln[l_s[c][a]]=csprintf (aln[l_s[c][a]],"%s", char_buf);
	        }
	     }


	A->len_aln=LEN;
	A->nseq=ns[0]+ns[1];

	free_int (pos0, -1);
	free_int (C, -1);
	free_int (D, -1);
	free_int (I, -1);
	free_int (trace, -1);
	free_int (LD, -1);
	free_char ( al, -1);
	vfree(buffer);
	vfree(char_buf);


	return score;
    }

int hasch_seq(char *seq, int **hs, int **lu,int ktup,char *alp)
    {
	static int a[10];

	int i,j,l,limit,code,flag;
	char residue;

	int alp_lu[10000];
	int alp_size;

	alp_size=alp[0];
	alp++;



	for ( i=0; i< alp_size; i++)
	    {
	      alp_lu[(int)alp[i]]=i;
	    }



	l=strlen (seq);
	limit = (int)   pow((double)(alp_size+1),(double)ktup);
	hs[0]=(int*)vcalloc ( l+1,sizeof (int));
	lu[0]=(int*)vcalloc ( limit+1, sizeof(int));


	if ( l==0)myexit(EXIT_FAILURE);

	for (i=1;i<=ktup;i++)
           a[i] = (int) pow((double)(alp_size+1),(double)(i-1));


	for(i=1;i<=(l-ktup+1);++i)
	        {
		code=0;
		flag=FALSE;
		for(j=1;j<=ktup;++j)
		   {
		   if (is_gap(seq[i+j-2])){flag=TRUE;break;}
		   else residue=alp_lu[(int)seq[i+j-2]];
		   code+=residue*a[j];
		   }

		if ( flag)continue;
		++code;

		if (lu[0][code])hs[0][i]=lu[0][code];
		lu[0][code]=i;
		}
	return 0;
    }



/*********************************************************************/
/*                                                                   */
/*                         KTUP_DP                                   */
/*                                                                   */
/*                                                                   */
/*********************************************************************/

/**************Hasch DAta Handling*******************************************************/

struct Hasch_data * free_ktup_hasch_data (struct Hasch_data *d);
struct Hasch_data * declare_ktup_hasch_data (struct Hasch_entry *e);
struct Hasch_data * allocate_ktup_hasch_data (struct Hasch_data *e, int action);

struct Hasch_data
{
 int *list;
};
typedef struct Hasch_data Hasch_data;
struct Hasch_data * free_ktup_hasch_data (struct Hasch_data *d)
{
  return allocate_ktup_hasch_data (d, FREE);
}
struct Hasch_data * declare_ktup_hasch_data (struct Hasch_entry *e)
{
  e->data=allocate_ktup_hasch_data (NULL,DECLARE);
  return e->data;
}

struct Hasch_data * allocate_ktup_hasch_data (struct Hasch_data *e, int action)
{
  static struct Hasch_data **heap;
  static int heap_size, free_heap, a;

  if ( action == 100)
    {
      fprintf ( stderr, "\nHeap size: %d, Free Heap: %d", heap_size, free_heap);
      return NULL;
    }
  else if ( action==DECLARE)
    {
      if ( free_heap==0)
	{
	  free_heap=100;
	  heap_size+=free_heap;
	  heap=(Hasch_data**)vrealloc (heap,heap_size*sizeof (struct Hasch_entry *));
	  for ( a=0; a<free_heap; a++)
	    {
	      (heap[a])=(Hasch_data*)vcalloc ( 1, sizeof ( struct Hasch_entry *));
	      (heap[a])->list=(int*)vcalloc ( 10, sizeof (int));
	      (heap[a])->list[0]=10;
	    }
	}
      return heap[--free_heap];
    }
  else if ( action==FREE)
    {
      heap[free_heap++]=e;
      e->list[1]=0;
      return NULL;
    }
  return NULL;
}


/**************Hasch DAta Handling*******************************************************/

int precomputed_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL)
    {
      int l1, l2, a, b, c;
      int nid=0, npos=0, id;
      int r1, r2, s1, s2;

      l1=strlen(A->seq_al[l_s[0][0]]);
      l2=strlen(A->seq_al[l_s[1][0]]);
      if (l1!=l2)
	{
	  fprintf ( stderr, "\nERROR: improper use of the function precomputed pairwise:[FATAL:%s]", PROGRAM);
	  crash ("");
	}
      else if ( l1==0)
	{
	  A->score_aln=A->score=0;
	  return 0;
	}

      for (npos=0, nid=0, a=0; a< ns[0]; a++)
	{
	  s1=l_s[0][a];

	  for (b=0; b< ns[1]; b++)
	    {
	      s2=l_s[1][b];
	      for ( c=0; c<l1; c++)
		{
		r1=A->seq_al[s1][c];
		r2=A->seq_al[s2][c];
		if ( is_gap(r1) || is_gap(r2));
		else
		  {
		    npos++;
		    nid+=(r1==r2);
		  }
		}
	    }
	}
      id=(npos==0)?0:((nid*100)/npos);
      A->score=A->score_aln=id;
      return A->score;
    }
int ktup_comparison_str ( char *seq1, char *seq2, const int ktup);
int ktup_comparison_hasch ( char *i_seq1, char *i_seq2, const int ktup);
int ktup_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL)
    {
      static char **gl;
      static int ng;
      char *seq1;
      char *seq2;

      int min_len=10;



      if ( !gl)
	gl=make_group_aa (&ng, "vasiliky");


      if ( ns[0]>1)seq1=sub_aln2cons_seq_mat (A, ns[0], l_s[0],"blosum62mt");
      else
	{
	  seq1=(char*)vcalloc ( strlen (A->seq_al[l_s[0][0]])+1, sizeof (char));
	  sprintf ( seq1, "%s",A->seq_al[l_s[0][0]]);
	}
      if ( ns[1]>1)seq2=sub_aln2cons_seq_mat (A, ns[1], l_s[1],"blosum62mt");
      else
	{
	  seq2=(char*)vcalloc ( strlen (A->seq_al[l_s[1][0]])+1, sizeof (char));
	  sprintf ( seq2, "%s",A->seq_al[l_s[1][0]]);
	}

      if ( strlen (seq1)<min_len || strlen (seq2)<min_len)
	{
	  Alignment *B;

	  ungap(seq1); ungap(seq2);
	  B=align_two_sequences ( seq1, seq2, "blosum62mt",-10, -1, "myers_miller_pair_wise");
	  A->score=A->score_aln=aln2sim(B, "idmat");
	  free_aln (B);
	  return A->score;
	}
      else
	{

	  string_convert (seq1, ng, gl);
	  string_convert (seq2, ng, gl);
	  A->score=A->score_aln=ktup_comparison (seq1,seq2, CL->ktup);
	}

      vfree (seq1); vfree (seq2);
      return A->score;
    }
int ktup_comparison( char *seq2, char *seq1, const int ktup)
{
  return ktup_comparison_hasch ( seq2, seq1, ktup);
}
int ktup_comparison_str ( char *seq2, char *seq1, const int ktup)
{
  int a,l1, l2,c1, c2, end, start;
  char *s1, *s2;
  double score=0;
  int max_dist=-1;

  if ( max_dist==-1)max_dist=MAX((strlen (seq1)),(strlen (seq2)));
  l1=strlen (seq1)-ktup;
  l2=strlen (seq2);


  for ( a=0; a< l1; a++)
    {
      c1=seq1[a+ktup];seq1[a+ktup]='\0';
      s1=seq1+a;

      start=((a-max_dist)<0)?0:a-max_dist;
      end=((a+max_dist)>=l2)?l2:a+max_dist;

      c2=seq2[end];seq2[end]='\0';
      s2=seq2+start;

      score+=(strstr(s2, s1)!=NULL)?1:0;

      seq1[a+ktup]=c1;
      seq2[end]=c2;
    }
  score/=(l1==0)?1:l1;
  score=((log(0.1+score)-log(0.1))/(log(1.1)-log(0.1)));

  return score*100;

}
int ktup_comparison_hasch ( char *i_seq1, char *i_seq2, const int ktup)
{
  /*Ktup comparison adapted from Rob Edgar, NAR, vol32, No1, 381, 2004*/
  /*1: hasch sequence 1
    2: Count the number of seq2 ktup found in seq1
  */

  char c;
  int key;

  static HaschT*H1;
  static char *pseq;
  Hasch_entry *e;
  char *s;
  int l, ls;
  int p, a, max_dist=-1;
  double score=0;



  if (!strm (i_seq1, pseq))
    {
      if (H1)
	{
	  hdestroy (H1, declare_ktup_hasch_data, free_ktup_hasch_data);
	  string2key (NULL, NULL);
	}
      H1=hasch_sequence ( i_seq1, ktup);
      vfree (pseq);pseq=(char*)vcalloc ( strlen (i_seq1)+1, sizeof (char));
      sprintf ( pseq, "%s", i_seq1);
    }

  ls=l=strlen (i_seq2);
  s=i_seq2;
  p=0;
  while (ls>ktup)
    {
      c=s[ktup];s[ktup]='\0';
      key=string2key (s, NULL);
      e=hsearch (H1,key,FIND, declare_ktup_hasch_data, free_ktup_hasch_data);

      if ( e==NULL);
      else if ( max_dist==-1)score++;
      else
	{
	  for ( a=1; a<=(e->data)->list[1]; a++)
	    if (FABS((p-(e->data)->list[a]))<=max_dist)
	      {score++; break;}
	}
      s[ktup]=c;s++;p++;ls--;
    }
  score/=(l-ktup);
  score=(log(0.1+score)-log(0.1))/(log(1.1)-log(0.1));

  if ( score>100) score=100;
  return (int)(score*100);
}

HaschT* hasch_sequence ( char *seq1, int ktup)
{
  char c;
  int key, offset=0, ls;
  HaschT *H;
  Hasch_entry *e;

  H=hcreate ( strlen (seq1), declare_ktup_hasch_data, free_ktup_hasch_data);
  ls=strlen (seq1);
  while (ls>=(ktup))
    {
      c=seq1[ktup];seq1[ktup]='\0';
      key=string2key (seq1, NULL);
      e=hsearch (H,key,FIND, declare_ktup_hasch_data, free_ktup_hasch_data);

      if (e==NULL)
	{
	 e=hsearch (H,key,ADD,declare_ktup_hasch_data,free_ktup_hasch_data);
	 (e->data)->list[++(e->data)->list[1]+1]=offset;
	}
      else
	{
	  if ((e->data)->list[0]==((e->data)->list[1]+2)){(e->data)->list[0]+=10;(e->data)->list=(int*)vrealloc ((e->data)->list,(e->data)->list[0]*sizeof (int));}
	  (e->data)->list[++(e->data)->list[1]+1]=offset;
	}
       seq1[ktup]=c;seq1++;ls--;
       offset++;
    }
  return H;
}



char *dayhoff_translate (char *seq1)
{
int l, a, c;
l=strlen (seq1);
 for ( a=0; a< l; a++)
  {
    c=tolower(seq1[a]);
    if ( strchr ("agpst", c))seq1[a]='a';
    else if (strchr ("denq", c))seq1[a]='d';
    else if (strchr ("fwy", c))seq1[a]='f';
    else if (strchr ("hkr", c))seq1[a]='h';
    else if (strchr ("ilmv", c))seq1[a]='i';
  }
return seq1;
}

int ** evaluate_diagonals_with_ktup ( Alignment *A, int *ns, int **l_s, Constraint_list *CL,int maximise,int n_groups, char **group_list, int ktup)
{
  /*Ktup comparison as in Rob Edgar, NAR, vol32, No1, 381, 2004*/
  char character;
  int key,ls;
  HaschT*H1, *H2;
  Hasch_entry *e1, *e2;
  char *s, *sb, *seq1, *seq2;
  int l1, l2;
  int score=0;
  int **diag,n_diag, ktup1, ktup2,a,b,c,d, **pos;
  int n_dots=0;

  pos=aln2pos_simple ( A,-1, ns, l_s);

  seq1=aln2cons_maj (A, ns[0], l_s[0], n_groups, group_list);
  seq2=aln2cons_maj (A, ns[1], l_s[1], n_groups, group_list);
  l1=strlen (seq1);
  l2=strlen (seq2);
  n_diag=l1+l2-1;


  diag=declare_int (n_diag+2, 3);
  for ( a=0; a<n_diag+2; a++)diag[a][0]=a;

  H1=hasch_sequence ( seq1, ktup);
  H2=hasch_sequence ( seq2, ktup);
  s=sb=(char*)vcalloc (strlen (seq1)+strlen (seq2)+1, sizeof (char));
  sprintf (s, "%s%s", seq1, seq2);

  ls=strlen(s);
  while (ls>=(ktup))
    {
      character=s[ktup];s[ktup]='\0';
      key=string2key (s, NULL);
      e1=hsearch (H1,key,FIND,declare_ktup_hasch_data, free_ktup_hasch_data);
      e2=hsearch (H2,key,FIND,declare_ktup_hasch_data, free_ktup_hasch_data);
      if ( !e2 || !e1);
      else
	{

	  for (b=2; b<(e1->data)->list[1]+2; b++)
	    for (c=2; c<(e2->data)->list[1]+2; c++)
	      {

		ktup1=(e1->data)->list[b];
		ktup2=(e2->data)->list[c];
		diag[(ktup2-ktup1)+l1][2]++;
		for (score=0, d=0; d<ktup; d++)
		  score+=(CL->get_dp_cost) ( A, pos, ns[0], l_s[0], ktup1+d, pos,ns[1], l_s[1], ktup2+d, CL);
		diag[(ktup2-ktup1)+l1][1]+=score;
		n_dots++;
	      }
	  (e1->data)->list[1]=(e2->data)->list[1]=0;
	}
      s[ktup]=character;s++;ls--;
    }

  sort_int (diag+1, 2, 1,0,n_diag-1);

  hdestroy (H1,declare_ktup_hasch_data, free_ktup_hasch_data); hdestroy (H2,declare_ktup_hasch_data, free_ktup_hasch_data);
  vfree (seq1); vfree (seq2);vfree (sb);free_int (pos, -1);
  return diag;
}
 /*********************************************************************/
/*                                                                   */
/*                         OLD FUNCTIONS                              */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
int ** evaluate_diagonals_with_ktup_1 ( Alignment *A, int *ns, int **l_s, Constraint_list *CL,int maximise,int n_groups, char **group_list, int ktup)
    {
   /*
    Reads in an alignmnent A, with two groups of sequences marked.
    1-Turn each group into a conscensus, using the group list identifier.
               -if the group list is left empty original symbols are used
    2-hasch the two sequences
    3-score each diagonal, sort the list and return it (diag_list)

        diag_list:

   */

    char *seq1, *seq2, *alphabet=NULL;
    int a,b,l1, l2, n_ktup,pos_ktup1, pos_ktup2, **pos;
    int *hasched_seq1, *hasched_seq2,*lu_seq1,*lu_seq2;
    int n_diag, **diag, current_diag, n_dots;
    static char *buf;
    pos=aln2pos_simple ( A,-1, ns, l_s);


    seq1=aln2cons_seq (A, ns[0], l_s[0], n_groups, group_list);
    seq2=aln2cons_seq (A, ns[1], l_s[1], n_groups, group_list);




    alphabet=get_alphabet (seq1,alphabet);
    alphabet=get_alphabet (seq2,alphabet);

    l1=strlen ( seq1);
    l2=strlen ( seq2);

    n_diag=l1+l2-1;
    diag=declare_int ( n_diag+2, 3);
    n_ktup=(int)pow ( (double)alphabet[0]+1, (double)ktup);


    hasch_seq(seq1, &hasched_seq1, &lu_seq1,ktup, alphabet);
    hasch_seq(seq2, &hasched_seq2, &lu_seq2,ktup, alphabet);




    /*EVALUATE THE DIAGONALS*/
    for ( a=0; a<= n_diag; a++)diag[a][0]=a;
    for ( n_dots=0,a=1; a<= n_ktup; a++)
        {
	    pos_ktup1=lu_seq1[a];
	    while (TRUE)
	          {
		  if (!pos_ktup1)break;
		  pos_ktup2=lu_seq2[a];
		  while (pos_ktup2)
		            {
			    current_diag=(pos_ktup2-pos_ktup1+l1);
			    for ( b=0; b< ktup; b++)
			        {
		         	    diag[current_diag][1]+=(CL->get_dp_cost) ( A, pos, ns[0], l_s[0], pos_ktup1+b-1, pos,ns[1], l_s[1], pos_ktup2+b-1, CL);
				    n_dots++;

				}
			    diag[current_diag][2]++;
			    pos_ktup2=hasched_seq2[pos_ktup2];
			    }
		  pos_ktup1=hasched_seq1[pos_ktup1];
		  }

	}
    if ( n_dots==0)
       {
	   if ( !buf)
	       {
	       buf=(char*)vcalloc ( 30, sizeof (30));
	       sprintf ( buf, "abcdefghijklmnopqrstuvwxyz");
	       }
	    vfree ( hasched_seq1);
	    vfree ( hasched_seq2);
	    vfree (lu_seq1);
	    vfree (lu_seq2);
	   return evaluate_diagonals_with_ktup ( A,ns,l_s, CL,maximise,1,&buf,1);
       }


    sort_int (diag+1, 2, 1,0, n_diag-1);
    vfree (seq1);
    vfree (seq2);
    vfree (alphabet);
    vfree ( hasched_seq1);
    vfree ( hasched_seq2);
    vfree (lu_seq1);
    vfree (lu_seq2);
    free_int (pos, -1);
    return diag;
    }
/////////////////////////////////////////////////////////////////

Constraint_list * hasch2constraint_list (Sequence*S, Constraint_list *CL)
{
  int a,b,c, n;
  SeqHasch h,*H=NULL;
  int *entry;
  int ktup=2;


  entry=(int*)vcalloc ( CL->entry_len+1, sizeof (int));

  for (a=0; a<S->nseq; a++)
    {
      H=seq2hasch (a, S->seq[a],ktup,H);
    }

  n=1;
  while (H[n])
    {
      h=H[n];

      for (a=0; a<h->n-2; a+=2)
	{
	  for (b=a+2; b<h->n; b+=2)
	    {

	      if (h->l[a]==h->l[b])continue;
	      else
		{
		  for (c=0; c<ktup; c++)
		    {
		      entry[SEQ1]=h->l[a];
		      entry[SEQ2]=h->l[b];
		      entry[R1]=h->l[a+1]+c;
		      entry[R2]=h->l[b+1]+c;
		      entry[WE]=100;
		      add_entry2list (entry,CL);
		    }
		}
	    }
	}
      n++;
    }

  return CL;
}
SeqHasch *cleanhasch       (SeqHasch *H)
{
  int n=1;
  SeqHasch *N;
  N=(hseq**)vcalloc (2, sizeof (SeqHasch));
  N[0]=H[0];

  while (H[n])
    {
      (H[n])->n=0;
      vfree ((H[n])->l);
      (H[n])->l=NULL;
      n++;
    }
  vfree (H);
  return N;
}
int hasch2sim        (SeqHasch *H, int nseq)
{
  int n=1;

  int a,cs, ps, ns;
  int id=0, tot=0;

  while (H[n])
    {
      for (ps=-1,ns=0,a=0; a<(H[n])->n; a+=2)
	{
	  //HERE ("%d--[%d %d]",n, (H[n])->l[a], (H[n])->l[a+1]);
	  cs=(H[n])->l[a];
	  if (cs!=ps)ns++;
	  ps=cs;
	}
      n++;
      if (ns==nseq)id++;
      tot++;
    }

  return (id*MAXID)/tot;
}
SeqHasch * seq2hasch (int i,char *seq, int ktup, SeqHasch *H)
{
  int a,b,l, n=0;
  SeqHasch h;


  if (!H)
    {
      H=(hseq**)vcalloc (2, sizeof (SeqHasch));
      H[0]=(hseq*)vcalloc (1, sizeof (hseq));
      n=1;
    }
  else
    {
      n=0;
      while (H[++n]);
    }

  l=strlen (seq);
  for (a=0; a<l-ktup; a++)
    {
      h=H[0];
      for (b=a; b<a+ktup; b++)
	{
	  char r;
	  r=seq[b];
	  if (!h->hl[r])  h->hl[r]=(hseq*)vcalloc (1, sizeof (hseq));
	  h=h->hl[r];
	}
      if (!h->l)
	{

	  h->n=2;
	  h->l=(int*)vcalloc (2, sizeof (int));
	  H=(hseq**)vrealloc (H,(n+2)*sizeof (SeqHasch));
	  H[n]=h;
	  n++;
	}
      else
	{
	  h->n+=2;
	  h->l=(int*)vrealloc (h->l, (h->n)*sizeof (int));
	}

      h->l[h->n-2]=i;
      h->l[h->n-1]=a;
    }
  return H;
}

