#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "dp_lib_header.h"
#include "define_header.h"

static int first_seq, last_seq;

static int nseqs;
static char **names;       /* the seq. names */

static double   **tmat;
static double 	*av;
static double 	*left_branch, *right_branch;
static int 	*tkill;


/*!
 *	\file util_make_tree.c
 *	\brief Source code tree algorithms
 */

NT_node ** seq2cw_tree ( Sequence *S, char *tree)
{
  Alignment *A;
  char *aln,command[1000];
  int tot_node;


  
  A=seq2clustalw_aln (S);
  aln=vtmpnam (NULL); if (!tree)tree=vtmpnam (NULL);
  
  output_fasta_aln (aln, A);
  
  sprintf ( command, "clustalw -infile=%s -newtree=%s -tree %s", aln, tree, TO_NULL_DEVICE);
  my_system (command);
  return read_tree (tree, &tot_node, S->nseq, S->name);
}

NT_node ** make_upgma_tree (  Alignment *A,int **distances,int gop, int gep, char **out_seq, char **out_seq_name, int out_nseq, char *tree_file, char *tree_mode)
  {

    if ( distances==NULL && A==NULL)
      {
	fprintf ( stderr, "\nError: make_nj_tree, must provide an alignment or a distance matrix [FATAL:%s]", PROGRAM);
	myexit (EXIT_FAILURE);
      }
    else if ( distances==NULL)
      {

	distances=get_dist_aln_array (A, "idmat");
      }
    return int_dist2upgma_tree (distances,A, A->nseq,tree_file);
  }
NT_node ** make_nj_tree (  Alignment *A,int **distances,int gop, int gep, char **out_seq, char **out_seq_name, int out_nseq, char *tree_file, char *tree_mode)
  {
   
    if ( distances==NULL && A==NULL)
      {
	fprintf ( stderr, "\nError: make_nj_tree, must provide an alignment or a distance matrix [FATAL:%s]", PROGRAM);
	myexit (EXIT_FAILURE);
      }
    else if ( distances==NULL)
      {
	distances=get_dist_aln_array (A, "idmat");
      }
    return int_dist2nj_tree (distances,out_seq_name, out_nseq,tree_file);
  }

NT_node ** int_dist2nj_tree (int **distances, char **out_seq_name, int out_nseq,  char *tree_file) 
        {
	  int a, b;
	  double **d;
	  NT_node **T;

	  d=declare_double( out_nseq, out_nseq);
	  for ( a=0; a< out_nseq; a++)
	    for ( b=0; b< out_nseq; b++)
	      {
		d[a][b]=distances[a][b];
	      }
	  T=dist2nj_tree ( d, out_seq_name, out_nseq, tree_file);
	  
	  free_double (d, -1);
	  return T;
	}
NT_node ** float_dist2nj_tree (float **distances, char **out_seq_name, int out_nseq,  char *tree_file)
        {
	  int a, b;
	  double **d;
	  NT_node **T;
	  
	  d=declare_double( out_nseq, out_nseq);
	  for ( a=0; a< out_nseq; a++)
	    for ( b=0; b< out_nseq; b++)
	      d[a][b]=distances[a][b];
	  T=dist2nj_tree ( d, out_seq_name, out_nseq, tree_file);
	  free_double (d, -1);
	  return T;
	}
NT_node ** dist2nj_tree (double **distances, char **out_seq_name, int out_nseq,  char *tree_file)
	{
	int a, b;
	double **d_dist;
	int tot_node=0;
	NT_node **T;

	if ( !tree_file)tree_file=vtmpnam(NULL);
	d_dist=declare_double( out_nseq+1, out_nseq+1);
	
	for ( a=0; a<out_nseq; a++)
		for ( b=0; b< out_nseq; b++)
			{
			  if ( a!=b)
				d_dist[a+1][b+1]=distances[a][b]/MAXID;
			  else d_dist[a+1][b+1]=0;
			}

	guide_tree ( tree_file, d_dist, out_seq_name, out_nseq);
	free_double (d_dist,-1);
	T=read_tree (tree_file,&tot_node,out_nseq,  out_seq_name);	

	return T;
	}


void guide_tree(char *fname, double **saga_tmat, char **saga_seq_name, int saga_nseq)
/* 
   Routine for producing unrooted NJ trees from seperately aligned
   pairwise distances.  This produces the GUIDE DENDROGRAMS in
   PHYLIP format.
*/
{    

        int i;
        FILE *fp;
        char **standard_tree;


	  
  	nseqs=saga_nseq;
	first_seq=1;
	last_seq=nseqs;

  	names=saga_seq_name;
  	tmat=saga_tmat;
	
        standard_tree   = (char **) vmalloc( (nseqs+1) * sizeof (char *) );
        for(i=0; i<nseqs+1; i++)
                standard_tree[i]  = (char *) vmalloc( (nseqs+1) * sizeof(char));
         	
        	
        nj_tree(standard_tree, nseqs);
        
        fp=fopen ( fname, "w");
        print_phylip_tree(standard_tree,fp,FALSE);

        vfree(left_branch);
        vfree(right_branch);
        vfree(tkill);
        vfree(av);
        
	for (i=0;i<nseqs+1;i++)
                vfree(standard_tree[i]);


        fclose(fp);

        
}


void nj_tree(char **tree_description, int nseq)
{
  if ( nseq<100)
    slow_nj_tree (tree_description);
  else
    fast_nj_tree (tree_description);
  
}



void slow_nj_tree(char **tree_description)
{
	int i;
	int l[4],nude,k;
	int nc,mini,minj,j,ii,jj;
	double fnseqs,fnseqs2=0,sumd;
	double diq,djq,dij,d2r,dr,dio,djo,da;
	double tmin,total,dmin;
	double bi,bj,b1,b2,b3,branch[4];
	int typei,typej;             /* 0 = node; 1 = OTU */
	
	fnseqs = (double)nseqs;

/*********************** First initialisation ***************************/
	
	

	mini = minj = 0;

	left_branch 	= (double *) vcalloc( (nseqs+2),sizeof (double)   );
	right_branch    = (double *) vcalloc( (nseqs+2),sizeof (double)   );
	tkill 		= (int *)    vcalloc( (nseqs+1),sizeof (int) );
	av   		= (double *) vcalloc( (nseqs+1),sizeof (double)   );

	for(i=1;i<=nseqs;++i) 
		{
		  tmat[i][i] = av[i] = 0.0;
		  tkill[i] = 0;
		}

/*********************** Enter The Main Cycle ***************************/

       	/**start main cycle**/
	for(nc=1; nc<=(nseqs-3); ++nc) 
	  {

	    sumd = 0.0;
	    for(j=2; j<=nseqs; ++j)
	      for(i=1; i<j; ++i) 
		{
		tmat[j][i] = tmat[i][j];
		sumd = sumd + tmat[i][j];
		}
	    tmin = 99999.0;

/*.................compute SMATij values and find the smallest one ........*/

	    for(jj=2; jj<=nseqs; ++jj) 
	      if(tkill[jj] != 1) 
		for(ii=1; ii<jj; ++ii)
		  if(tkill[ii] != 1) 
		    {
		      diq = djq = 0.0;
		      
		      for(i=1; i<=nseqs; ++i) 
			{
			  diq = diq + tmat[i][ii];
			  djq = djq + tmat[i][jj];
			}
		      dij = tmat[ii][jj];
		      d2r = diq + djq - (2.0*dij);
		      dr  = sumd - dij -d2r;
		      fnseqs2 = fnseqs - 2.0;
		      total= d2r+ fnseqs2*dij +dr*2.0;
		      total= total / (2.0*fnseqs2);
		      
		      /* commented out to have consistent results with CYGWIN: if(total < tmin)"*/
		      if(total < tmin) 
			{
			  if ( tmin-total>MY_EPS)
			  tmin = total;
			  mini = ii;
			  minj = jj;
			}
		    }


/*.................compute branch lengths and print the results ........*/


	    dio = djo = 0.0;
	    for(i=1; i<=nseqs; ++i) {
	      dio = dio + tmat[i][mini];
	      djo = djo + tmat[i][minj];
	    }
	    
	    dmin = tmat[mini][minj];
	    dio = (dio - dmin) / fnseqs2;
	    djo = (djo - dmin) / fnseqs2;
	    bi = (dmin + dio - djo) * 0.5;
	    bj = dmin - bi;
	    bi = bi - av[mini];
	    bj = bj - av[minj];
	    
	    if( av[mini] > 0.0 )
	      typei = 0;
	    else
	      typei = 1;
	    if( av[minj] > 0.0 )
	      typej = 0;
	    else
	      typej = 1;
	    
	    
	    
/* 
   set negative branch lengths to zero.  Also set any tiny positive
   branch lengths to zero.
*/		
	    if( fabs(bi) < 0.0001) bi = 0.0;
	    if( fabs(bj) < 0.0001) bj = 0.0;

	    


	    left_branch[nc] = bi;
	    right_branch[nc] = bj;
	    
	    for(i=1; i<=nseqs; i++)
	      tree_description[nc][i] = 0;
	    
	    if(typei == 0) { 
	      for(i=nc-1; i>=1; i--)
		if(tree_description[i][mini] == 1) {
		  for(j=1; j<=nseqs; j++)  
		    if(tree_description[i][j] == 1)
		      tree_description[nc][j] = 1;
		  break;
		}
	    }
	    else
	      tree_description[nc][mini] = 1;
	    
	    if(typej == 0) {
	      for(i=nc-1; i>=1; i--) 
		if(tree_description[i][minj] == 1) {
		  for(j=1; j<=nseqs; j++)  
		    if(tree_description[i][j] == 1)
		      tree_description[nc][j] = 1;
		  break;
		}
	    }
	    else
	      tree_description[nc][minj] = 1;
			

/* 
   Here is where the -0.00005 branch lengths come from for 3 or more
   identical seqs.
*/
/*		if(dmin <= 0.0) dmin = 0.0001; */
	    if(dmin <= 0.0) dmin = 0.000001;
	    av[mini] = dmin * 0.5;
	    
 /*........................Re-initialisation................................*/
	    
	    fnseqs = fnseqs - 1.0;
	    tkill[minj] = 1;

	    for(j=1; j<=nseqs; ++j) 
	      if( tkill[j] != 1 ) {
		da = ( tmat[mini][j] + tmat[minj][j] ) * 0.5;
		if( (mini - j) < 0 ) 
		  tmat[mini][j] = da;
		if( (mini - j) > 0)
		  tmat[j][mini] = da;
	      }
	    
	    for(j=1; j<=nseqs; ++j)
	      tmat[minj][j] = tmat[j][minj] = 0.0;
	    
	    

	  }
	/*end main cycle**/
	
/******************************Last Cycle (3 Seqs. left)********************/

	nude = 1;


	for(i=1; i<=nseqs; ++i)
		if( tkill[i] != 1 ) {
			l[nude] = i;
			nude = nude + 1;
		}

	b1 = (tmat[l[1]][l[2]] + tmat[l[1]][l[3]] - tmat[l[2]][l[3]]) * 0.5;
	b2 =  tmat[l[1]][l[2]] - b1;
	b3 =  tmat[l[1]][l[3]] - b1;

	branch[1] = b1 - av[l[1]];
	branch[2] = b2 - av[l[2]];
	branch[3] = b3 - av[l[3]];

/* Reset tiny negative and positive branch lengths to zero */
	if( fabs(branch[1]) < 0.0001) branch[1] = 0.0;
	if( fabs(branch[2]) < 0.0001) branch[2] = 0.0;
	if( fabs(branch[3]) < 0.0001) branch[3] = 0.0;

	left_branch[nseqs-2] = branch[1];
	left_branch[nseqs-1] = branch[2];
	left_branch[nseqs]   = branch[3];

	for(i=1; i<=nseqs; i++)
		tree_description[nseqs-2][i] = 0;

	

	for(i=1; i<=3; ++i) {
	   if( av[l[i]] > 0.0) {
	
		for(k=nseqs-3; k>=1; k--)
			if(tree_description[k][l[i]] == 1) {
				for(j=1; j<=nseqs; j++)
				 	if(tree_description[k][j] == 1)
					    tree_description[nseqs-2][j] = i;
				break;
			}
	   }
	   else  {
	      
		tree_description[nseqs-2][l[i]] = i;
	   }
	   if(i < 3) {
	   }
	}
}

void print_phylip_tree(char **tree_description, FILE *tree, int bootstrap)
{

	fprintf(tree,"(\n");
 
	two_way_split(tree_description, tree, nseqs-2,1,bootstrap);
	fprintf(tree,":%7.5f,\n",left_branch[nseqs-2]);
	two_way_split(tree_description, tree, nseqs-2,2,bootstrap);
	fprintf(tree,":%7.5f,\n",left_branch[nseqs-1]);
	two_way_split(tree_description, tree, nseqs-2,3,bootstrap);

	fprintf(tree,":%7.5f)",left_branch[nseqs]);
	if (bootstrap) fprintf(tree,"TRICHOTOMY");
	fprintf(tree,";\n");
}


void two_way_split
(char **tree_description, FILE *tree, int start_row, int flag, int bootstrap)
{
	int row, new_row, col, test_col=0;
	int single_seq;

	if(start_row != nseqs-2) fprintf(tree,"(\n"); 

	for(col=1; col<=nseqs; col++) {
		if(tree_description[start_row][col] == flag) {
			test_col = col;
			break;
		}
	}

	single_seq = TRUE;
	for(row=start_row-1; row>=1; row--) 
		if(tree_description[row][test_col] == 1) {
			single_seq = FALSE;
			new_row = row;
			break;
		}

	if(single_seq) {
		tree_description[start_row][test_col] = 0;
		fprintf(tree,"%s",names[test_col+0-1]);
	}
	else {
		for(col=1; col<=nseqs; col++) {
		    if((tree_description[start_row][col]==1)&&
		       (tree_description[new_row][col]==1))
				tree_description[start_row][col] = 0;
		}
		two_way_split(tree_description, tree, new_row, (int)1, bootstrap);
	}

	if(start_row == nseqs-2) {
/*		if (bootstrap && (flag==1)) fprintf(tree,"[TRICHOTOMY]");
*/
		return;
	}

	fprintf(tree,":%7.5f,\n",left_branch[start_row]);

	for(col=1; col<=nseqs; col++) 
		if(tree_description[start_row][col] == flag) {
			test_col = col;
			break;
		}
	
	single_seq = TRUE;
	for(row=start_row-1; row>=1; row--) 
		if(tree_description[row][test_col] == 1) {
			single_seq = FALSE;
			new_row = row;
			break;
		}

	if(single_seq) {
		tree_description[start_row][test_col] = 0;
		fprintf(tree,"%s",names[test_col+0-1]);
	}
	else {
		for(col=1; col<=nseqs; col++) {
		    if((tree_description[start_row][col]==1)&&
		       (tree_description[new_row][col]==1))
				tree_description[start_row][col] = 0;
		}
		two_way_split(tree_description, tree, new_row, (int)1, bootstrap);
	}

	fprintf(tree,":%7.5f)\n",right_branch[start_row]);
	
	
}




/****************************************************************************
 * [ Improvement ideas in fast_nj_tree() ] by DDBJ & FUJITSU Limited.
 *						written by Tadashi Koike
 *						(takoike@genes.nig.ac.jp)
 *******************
 * <IMPROVEMENT 1> : Store the value of sum of the score to temporary array,
 *                   and use again and again.
 *
 *	In the main cycle, these are calculated again and again :
 *	    diq = sum of tmat[n][ii]   (n:1 to last_seq-first_seq+1),
 *	    djq = sum of tmat[n][jj]   (n:1 to last_seq-first_seq+1),
 *	    dio = sum of tmat[n][mini] (n:1 to last_seq-first_seq+1),
 *	    djq = sum of tmat[n][minj] (n:1 to last_seq-first_seq+1)
 *		// 'last_seq' and 'first_seq' are both constant values //
 *	and the result of above calculations is always same until 
 *	a best pair of neighbour nodes is joined.
 *
 *	So, we change the logic to calculate the sum[i] (=sum of tmat[n][i]
 *	(n:1 to last_seq-first_seq+1)) and store it to array, before
 *	beginning to find a best pair of neighbour nodes, and after that 
 *	we use them again and again.
 *
 *	    tmat[i][j]
 *	              1   2   3   4   5
 *	            +---+---+---+---+---+
 *	          1 |   |   |   |   |   |
 *	            +---+---+---+---+---+
 *	          2 |   |   |   |   |   |  1) calculate sum of tmat[n][i]
 *	            +---+---+---+---+---+        (n: 1 to last_seq-first_seq+1)
 *	          3 |   |   |   |   |   |  2) store that sum value to sum[i]
 *	            +---+---+---+---+---+
 *	          4 |   |   |   |   |   |  3) use sum[i] during finding a best
 *	            +---+---+---+---+---+     pair of neibour nodes.
 *	          5 |   |   |   |   |   |
 *	            +---+---+---+---+---+
 *	              |   |   |   |   |
 *	              V   V   V   V   V  Calculate sum , and store it to sum[i]
 *	            +---+---+---+---+---+
 *	     sum[i] |   |   |   |   |   |
 *	            +---+---+---+---+---+
 *
 *	At this time, we thought that we use upper triangle of the matrix
 *	because tmat[i][j] is equal to tmat[j][i] and tmat[i][i] is equal 
 *	to zero. Therefore, we prepared sum_rows[i] and sum_cols[i] instead 
 *	of sum[i] for storing the sum value.
 *
 *	    tmat[i][j]
 *	              1   2   3   4   5     sum_cols[i]
 *	            +---+---+---+---+---+     +---+
 *	          1     | # | # | # | # | --> |   | ... sum of tmat[1][2..5]
 *	            + - +---+---+---+---+     +---+
 *	          2         | # | # | # | --> |   | ... sum of tmat[2][3..5]
 *	            + - + - +---+---+---+     +---+
 *	          3             | # | # | --> |   | ... sum of tmat[3][4..5]
 *	            + - + - + - +---+---+     +---+
 *	          4                 | # | --> |   | ... sum of tmat[4][5]
 *	            + - + - + - + - +---+     +---+
 *	          5                     | --> |   | ... zero
 *	            + - + - + - + - + - +     +---+
 *	              |   |   |   |   |
 *	              V   V   V   V   V  Calculate sum , sotre to sum[i]
 *	            +---+---+---+---+---+
 *	sum_rows[i] |   |   |   |   |   |
 *	            +---+---+---+---+---+
 *	              |   |   |   |   |
 *	              |   |   |   |   +----- sum of tmat[1..4][5]
 *	              |   |   |   +--------- sum of tmat[1..3][4]
 *	              |   |   +------------- sum of tmat[1..2][3]
 *	              |   +----------------- sum of tmat[1][2]
 *	              +--------------------- zero
 *
 *	And we use (sum_rows[i] + sum_cols[i]) instead of sum[i].
 *
 *******************
 * <IMPROVEMENT 2> : We manage valid nodes with chain list, instead of
 *                   tkill[i] flag array.
 *
 *	In original logic, invalid(killed?) nodes after nodes-joining
 *	are managed with tkill[i] flag array (set to 1 when killed).
 *	By this method, it is conspicuous to try next node but skip it
 *	at the latter of finding a best pair of neighbor nodes.
 *
 *	So, we thought that we managed valid nodes by using a chain list 
 *	as below:
 *
 *	1) declare the list structure.
 *		struct {
 *		    int n;		// entry number of node.
 *		    void *prev;		// pointer to previous entry.
 *		    void *next;		// pointer to next entry.
 *		}
 *	2) construct a valid node list.
 *
 *       +-----+    +-----+    +-----+    +-----+        +-----+
 * NULL<-|prev |<---|prev |<---|prev |<---|prev |<- - - -|prev |
 *       |  0  |    |  1  |    |  2  |    |  3  |        |  n  |
 *       | next|--->| next|--->| next|--->| next|- - - ->| next|->NULL
 *       +-----+    +-----+    +-----+    +-----+        +-----+
 *
 *	3) when finding a best pair of neighbor nodes, we use
 *	   this chain list as loop counter.
 *
 *	4) If an entry was killed by node-joining, this chain list is
 *	   modified to remove that entry.
 *
 *	   EX) remove the entry No 2.
 *       +-----+    +-----+               +-----+        +-----+
 * NULL<-|prev |<---|prev |<--------------|prev |<- - - -|prev |
 *       |  0  |    |  1  |               |  3  |        |  n  |
 *       | next|--->| next|-------------->| next|- - - ->| next|->NULL
 *       +-----+    +-----+               +-----+        +-----+
 *                             +-----+
 *                       NULL<-|prev |
 *                             |  2  |
 *                             | next|->NULL
 *                             +-----+
 *
 *	By this method, speed is up at the latter of finding a best pair of
 *	neighbor nodes.
 *
 *******************
 * <IMPROVEMENT 3> : Cut the frequency of division.
 *
 * At comparison between 'total' and 'tmin' in the main cycle, total is
 * divided by (2.0*fnseqs2) before comparison.  If N nodes are available, 
 * that division happen (N*(N-1))/2 order.
 *
 * We thought that the comparison relation between tmin and total/(2.0*fnseqs2)
 * is equal to the comparison relation between (tmin*2.0*fnseqs2) and total.
 * Calculation of (tmin*2.0*fnseqs2) is only one time. so we stop dividing
 * a total value and multiply tmin and (tmin*2.0*fnseqs2) instead.
 *
 *******************
 * <IMPROVEMENT 4> : some transformation of the equation (to cut operations).
 *
 * We transform an equation of calculating 'total' in the main cycle.
 *
 */


void fast_nj_tree(char **tree_description)
{
	register int i;
	int l[4],nude,k;
	int nc,mini,minj,j,ii,jj;
	double fnseqs,fnseqs2=0,sumd;
	double diq,djq,dij,dio,djo,da;
	double tmin,total,dmin;
	double bi,bj,b1,b2,b3,branch[4];
	int typei,typej;             /* 0 = node; 1 = OTU */

	/* IMPROVEMENT 1, STEP 0 : declare  variables */
	double *sum_cols, *sum_rows, *join;

	/* IMPROVEMENT 2, STEP 0 : declare  variables */
	int loop_limit;
	typedef struct _ValidNodeID {
	    int n;
	    struct _ValidNodeID *prev;
	    struct _ValidNodeID *next;
	} ValidNodeID;
	ValidNodeID *tvalid, *lpi, *lpj, *lpii, *lpjj, *lp_prev, *lp_next;

	/*
	 * correspondence of the loop counter variables.
	 *   i .. lpi->n,	ii .. lpii->n
	 *   j .. lpj->n,	jj .. lpjj->n
	 */

	fnseqs = (double)last_seq-first_seq+1;

/*********************** First initialisation ***************************/
	

	if (fnseqs == 2) {
		return;
	}

	mini = minj = 0;

	left_branch 	= (double *) ckalloc( (nseqs+2) * sizeof (double)   );
	right_branch    = (double *) ckalloc( (nseqs+2) * sizeof (double)   );
	tkill 		= (int *) ckalloc( (nseqs+1) * sizeof (int) );
	av   		= (double *) ckalloc( (nseqs+1) * sizeof (double)   );

	/* IMPROVEMENT 1, STEP 1 : Allocate memory */
	sum_cols	= (double *) ckalloc( (nseqs+1) * sizeof (double)   );
	sum_rows	= (double *) ckalloc( (nseqs+1) * sizeof (double)   );
	join		= (double *) ckalloc( (nseqs+1) * sizeof (double)   );

	/* IMPROVEMENT 2, STEP 1 : Allocate memory */
	tvalid	= (ValidNodeID *) ckalloc( (nseqs+1) * sizeof (ValidNodeID) );
	/* tvalid[0] is special entry in array. it points a header of valid entry list */
	tvalid[0].n = 0;
	tvalid[0].prev = NULL;
	tvalid[0].next = &tvalid[1];

	/* IMPROVEMENT 2, STEP 2 : Construct and initialize the entry chain list */
	for(i=1, loop_limit = last_seq-first_seq+1,
		lpi=&tvalid[1], lp_prev=&tvalid[0], lp_next=&tvalid[2] ;
		i<=loop_limit ;
		++i, ++lpi, ++lp_prev, ++lp_next)
		{
		tmat[i][i] = av[i] = 0.0;
		tkill[i] = 0;
		lpi->n = i;
		lpi->prev = lp_prev;
		lpi->next = lp_next;

		/* IMPROVEMENT 1, STEP 2 : Initialize arrays */
		sum_cols[i] = sum_rows[i] = join[i] = 0.0;
		}
	tvalid[loop_limit].next = NULL;

	/*
	 * IMPROVEMENT 1, STEP 3 : Calculate the sum of score value that 
	 * is sequence[i] to others.
	 */
	sumd = 0.0;
	for (lpj=tvalid[0].next ; lpj!=NULL ; lpj = lpj->next) {
		double tmp_sum = 0.0;
		j = lpj->n;
		/* calculate sum_rows[j] */
		for (lpi=tvalid[0].next ; lpi->n < j ; lpi = lpi->next) {
			i = lpi->n;
			tmp_sum += tmat[i][j];
			/* tmat[j][i] = tmat[i][j]; */
		}
		sum_rows[j] = tmp_sum;

		tmp_sum = 0.0;
		/* Set lpi to that lpi->n is greater than j */
		if ((lpi != NULL) && (lpi->n == j)) {
			lpi = lpi->next;
		}
		/* calculate sum_cols[j] */
		for( ; lpi!=NULL ; lpi = lpi->next) {
			i = lpi->n;
			tmp_sum += tmat[j][i];
			/* tmat[i][j] = tmat[j][i]; */
		}
		sum_cols[j] = tmp_sum;
	}

/*********************** Enter The Main Cycle ***************************/

	for(nc=1, loop_limit = (last_seq-first_seq+1-3); nc<=loop_limit; ++nc) {

		sumd = 0.0;
		/* IMPROVEMENT 1, STEP 4 : use sum value */
		for(lpj=tvalid[0].next ; lpj!=NULL ; lpj = lpj->next) {
			sumd += sum_cols[lpj->n];
		}

		/* IMPROVEMENT 3, STEP 0 : multiply tmin and 2*fnseqs2 */
		fnseqs2 = fnseqs - 2.0;		/* Set fnseqs2 at this point. */
		tmin = 99999.0 * 2.0 * fnseqs2;


/*.................compute SMATij values and find the smallest one ........*/

		mini = minj = 0;

		/* jj must starts at least 2 */
		if ((tvalid[0].next != NULL) && (tvalid[0].next->n == 1)) {
			lpjj = tvalid[0].next->next;
		} else {
			lpjj = tvalid[0].next;
		}

		for( ; lpjj != NULL; lpjj = lpjj->next) {
			jj = lpjj->n;
			for(lpii=tvalid[0].next ; lpii->n < jj ; lpii = lpii->next) {
				ii = lpii->n;
				diq = djq = 0.0;

				/* IMPROVEMENT 1, STEP 4 : use sum value */
				diq = sum_cols[ii] + sum_rows[ii];
				djq = sum_cols[jj] + sum_rows[jj];
				/*
				 * always ii < jj in this point. Use upper
				 * triangle of score matrix.
				 */
				dij = tmat[ii][jj];

				/*
				 * IMPROVEMENT 3, STEP 1 : fnseqs2 is
				 * already calculated.
				 */
				/* fnseqs2 = fnseqs - 2.0 */

				/* IMPROVEMENT 4 : transform the equation */
  /*-------------------------------------------------------------------*
   * OPTIMIZE of expression 'total = d2r + fnseqs2*dij + dr*2.0'       *
   * total = d2r + fnseq2*dij + 2.0*dr                                 *
   *       = d2r + fnseq2*dij + 2(sumd - dij - d2r)                    *
   *       = d2r + fnseq2*dij + 2*sumd - 2*dij - 2*d2r                 *
   *       =       fnseq2*dij + 2*sumd - 2*dij - 2*d2r + d2r           *
   *       = fnseq2*dij + 2*sumd - 2*dij - d2r                         *
   *       = fnseq2*dij + 2*sumd - 2*dij - (diq + djq - 2*dij)         *
   *       = fnseq2*dij + 2*sumd - 2*dij - diq - djq + 2*dij           *
   *       = fnseq2*dij + 2*sumd - 2*dij + 2*dij - diq - djq           *
   *       = fnseq2*dij + 2*sumd  - diq - djq                          *
   *-------------------------------------------------------------------*/
				total = fnseqs2*dij + 2.0*sumd  - diq - djq;

				/* 
				 * IMPROVEMENT 3, STEP 2 : abbrevlate
				 * the division on comparison between 
				 * total and tmin.
				 */
				/* total = total / (2.0*fnseqs2); */

				if(total < tmin) {
					tmin = total;
					mini = ii;
					minj = jj;
				}
			}
		}

		/* MEMO: always ii < jj in avobe loop, so mini < minj */

/*.................compute branch lengths and print the results ........*/


		dio = djo = 0.0;

		/* IMPROVEMENT 1, STEP 4 : use sum value */
		dio = sum_cols[mini] + sum_rows[mini];
		djo = sum_cols[minj] + sum_rows[minj];

		dmin = tmat[mini][minj];
		dio = (dio - dmin) / fnseqs2;
		djo = (djo - dmin) / fnseqs2;
		bi = (dmin + dio - djo) * 0.5;
		bj = dmin - bi;
		bi = bi - av[mini];
		bj = bj - av[minj];

		if( av[mini] > 0.0 )
			typei = 0;
		else
			typei = 1;
		if( av[minj] > 0.0 )
			typej = 0;
		else
			typej = 1;


/* 
   set negative branch lengths to zero.  Also set any tiny positive
   branch lengths to zero.
*/		if( fabs(bi) < 0.0001) bi = 0.0;
		if( fabs(bj) < 0.0001) bj = 0.0;


	    	left_branch[nc] = bi;
	    	right_branch[nc] = bj;

		for(i=1; i<=last_seq-first_seq+1; i++)
			tree_description[nc][i] = 0;

	     	if(typei == 0) { 
			for(i=nc-1; i>=1; i--)
				if(tree_description[i][mini] == 1) {
					for(j=1; j<=last_seq-first_seq+1; j++)  
					     if(tree_description[i][j] == 1)
						    tree_description[nc][j] = 1;
					break;
				}
		}
		else
			tree_description[nc][mini] = 1;

		if(typej == 0) {
			for(i=nc-1; i>=1; i--) 
				if(tree_description[i][minj] == 1) {
					for(j=1; j<=last_seq-first_seq+1; j++)  
					     if(tree_description[i][j] == 1)
						    tree_description[nc][j] = 1;
					break;
				}
		}
		else
			tree_description[nc][minj] = 1;
			

/* 
   Here is where the -0.00005 branch lengths come from for 3 or more
   identical seqs.
*/
/*		if(dmin <= 0.0) dmin = 0.0001; */
                if(dmin <= 0.0) dmin = 0.000001;
		av[mini] = dmin * 0.5;

/*........................Re-initialisation................................*/

		fnseqs = fnseqs - 1.0;
		tkill[minj] = 1;

		/* IMPROVEMENT 2, STEP 3 : Remove tvalid[minj] from chain list. */
		/* [ Before ]
		 *  +---------+        +---------+        +---------+       
		 *  |prev     |<-------|prev     |<-------|prev     |<---
		 *  |    n    |        | n(=minj)|        |    n    |
		 *  |     next|------->|     next|------->|     next|----
		 *  +---------+        +---------+        +---------+ 
		 *
		 * [ After ]
		 *  +---------+                           +---------+       
		 *  |prev     |<--------------------------|prev     |<---
		 *  |    n    |                           |    n    |
		 *  |     next|-------------------------->|     next|----
		 *  +---------+                           +---------+ 
		 *                     +---------+
		 *              NULL---|prev     |
		 *                     | n(=minj)|
		 *                     |     next|---NULL
		 *                     +---------+ 
		 */
		(tvalid[minj].prev)->next = tvalid[minj].next;
		if (tvalid[minj].next != NULL) {
			(tvalid[minj].next)->prev = tvalid[minj].prev;
		}
		tvalid[minj].prev = tvalid[minj].next = NULL;

		/* IMPROVEMENT 1, STEP 5 : re-calculate sum values. */
		for(lpj=tvalid[0].next ; lpj != NULL ; lpj = lpj->next) {
			double tmp_di = 0.0;
			double tmp_dj = 0.0;
			j = lpj->n;

			/* 
			 * subtrace a score value related with 'minj' from
			 * sum arrays .
			 */
			if (j < minj) {
				tmp_dj = tmat[j][minj];
				sum_cols[j] -= tmp_dj;
			} else if (j > minj) {
				tmp_dj = tmat[minj][j];
				sum_rows[j] -= tmp_dj;
			} /* nothing to do when j is equal to minj. */
			

			/* 
			 * subtrace a score value related with 'mini' from
			 * sum arrays .
			 */
			if (j < mini) {
				tmp_di = tmat[j][mini];
				sum_cols[j] -= tmp_di;
			} else if (j > mini) {
				tmp_di = tmat[mini][j];
				sum_rows[j] -= tmp_di;
			} /* nothing to do when j is equal to mini. */

			/* 
			 * calculate a score value of the new inner node.
			 * then, store it temporary to join[] array.
			 */
			join[j] = (tmp_dj + tmp_di) * 0.5;
		}

		/* 
		 * 1)
		 * Set the score values (stored in join[]) into the matrix,
		 * row/column position is 'mini'.
		 * 2)
		 * Add a score value of the new inner node to sum arrays.
		 */
		for(lpj=tvalid[0].next ; lpj != NULL; lpj = lpj->next) {
			j = lpj->n;
			if (j < mini) {
				tmat[j][mini] = join[j];
				sum_cols[j] += join[j];
			} else if (j > mini) {
				tmat[mini][j] = join[j];
				sum_rows[j] += join[j];
			} /* nothing to do when j is equal to mini. */
		}

		/* Re-calculate sum_rows[mini],sum_cols[mini]. */
		sum_cols[mini] = sum_rows[mini] = 0.0;

		/* calculate sum_rows[mini] */
		da = 0.0;
		for(lpj=tvalid[0].next ; lpj->n < mini ; lpj = lpj->next) {
                      da += join[lpj->n];
		}
		sum_rows[mini] = da;

		/* skip if 'lpj->n' is equal to 'mini' */
		if ((lpj != NULL) && (lpj->n == mini)) {
			lpj = lpj->next;
		}

		/* calculate sum_cols[mini] */
		da = 0.0;
		for( ; lpj != NULL; lpj = lpj->next) {
                      da += join[lpj->n];
		}
		sum_cols[mini] = da;

		/*
		 * Clean up sum_rows[minj], sum_cols[minj] and score matrix
		 * related with 'minj'.
		 */
		sum_cols[minj] = sum_rows[minj] = 0.0;
		for(j=1; j<=last_seq-first_seq+1; ++j)
			tmat[minj][j] = tmat[j][minj] = join[j] = 0.0;


/****/	}						/*end main cycle**/

/******************************Last Cycle (3 Seqs. left)********************/

	nude = 1;

	for(lpi=tvalid[0].next; lpi != NULL; lpi = lpi->next) {
		l[nude] = lpi->n;
		++nude;
	}

	b1 = (tmat[l[1]][l[2]] + tmat[l[1]][l[3]] - tmat[l[2]][l[3]]) * 0.5;
	b2 =  tmat[l[1]][l[2]] - b1;
	b3 =  tmat[l[1]][l[3]] - b1;
 
	branch[1] = b1 - av[l[1]];
	branch[2] = b2 - av[l[2]];
	branch[3] = b3 - av[l[3]];

/* Reset tiny negative and positive branch lengths to zero */
	if( fabs(branch[1]) < 0.0001) branch[1] = 0.0;
	if( fabs(branch[2]) < 0.0001) branch[2] = 0.0;
	if( fabs(branch[3]) < 0.0001) branch[3] = 0.0;

	left_branch[last_seq-first_seq+1-2] = branch[1];
	left_branch[last_seq-first_seq+1-1] = branch[2];
	left_branch[last_seq-first_seq+1]   = branch[3];

	for(i=1; i<=last_seq-first_seq+1; i++)
		tree_description[last_seq-first_seq+1-2][i] = 0;


	for(i=1; i<=3; ++i) {
	   if( av[l[i]] > 0.0) {

		for(k=last_seq-first_seq+1-3; k>=1; k--)
			if(tree_description[k][l[i]] == 1) {
				for(j=1; j<=last_seq-first_seq+1; j++)
				 	if(tree_description[k][j] == 1)
					    tree_description[last_seq-first_seq+1-2][j] = i;
				break;
			}
	   }
	   else  {
		tree_description[last_seq-first_seq+1-2][l[i]] = i;
	   }
	   if(i < 3) {;
	   }
	}
	ckfree(sum_cols);
	ckfree(sum_rows);
	ckfree(join);
	ckfree(tvalid);
}
//////////////////////////////////////////////////////////////////////////////
//
//                              UPGMA_aln
//////////////////////////////////////////////////////////////////////////////

Alignment * upgma_merge_aln_rows (Alignment *A, int *ns, int **ls,int N, int**mat,int *used, int *n,  Constraint_list *CL);
int upgma_pair_wise (Alignment *A, int *ls0, int ns0, int *ls2, int ns2, Constraint_list *CL);


Alignment * upgma_tree_aln  ( Alignment*A, int nseq, Constraint_list *CL)
{
  int a, b,n, *used;
  static int **mat;
  int **ls;
  int *ns;
  nseq=(CL->S)->nseq;
  mat=declare_int (nseq, nseq);
  ls=declare_int  (nseq,nseq);
  ns=vcalloc (nseq,sizeof (int));
  
  for (a=0; a<nseq-1; a++)
    for (b=a+1; b<nseq; b++)
      {
	ns[a]=ns[b]=1;
	ls[a][0]=a;
	ls[b][0]=b;
	mat[a][b]=mat[b][a]=upgma_pair_wise(A,ls[a],ns[a],ls[b],ns[b],CL);
	HERE ("%d %d", a, b, mat[a][b]);
      }
  
  used=vcalloc (nseq, sizeof (int));
  n=nseq;
  while (n>1)
    {
      upgma_merge_aln_rows (A,ns, ls,nseq, mat, used, &n,CL);
    }    
  print_aln (A);
  HERE ("finished");
  free_int ( mat, -1);
  free_int (ls, -1);
  vfree (ns);
  return A;
}

Alignment * upgma_merge_aln_rows (Alignment *A, int *ns, int **ls,int N, int**mat,int *used, int *n,  Constraint_list *CL)
{
  
  int a, b, w, best_a, best_b, set;
  float best_s;
    
  for (set=0,a=0; a<N-1; a++)
    {
      if (used[a])continue;
      for (b=a+1; b<N; b++)
	{
	  if (used[b])continue;
	  w=mat[a][b];
	  if ( !set || w>best_s)
	    {
	      best_s=w;
	      best_a=a;
	      best_b=b;
	      set=1;
	    }
	}
    }
  used[best_b]=1;
  
  //merge a and b
  mat[best_a][best_b]=upgma_pair_wise (A, ls[best_a], ns[best_a], ls[best_b], ns[best_b], CL);
  for (a=0; a<ns[best_b]; a++)ls[best_a][ns[best_a]++]=ls[best_b][a];
  ns[best_b]=0;
  
  //update the a row
  for (a=0; a< A->nseq; a++)
    {
      if (a!=best_a && !used[a])
	mat[best_a][a]=mat[a][best_a]=upgma_pair_wise (A, ls[best_a], ns[best_a], ls[a], ns[a], CL);
    }
  
  HERE ("DONE");
  
  n[0]--;
  return A;
}  
int upgma_pair_wise (Alignment *A, int *ls0, int ns0, int *ls1, int ns1, Constraint_list *CL)
{
  static int **ls;
  static int *ns;
  static int *fl;
  int a, b, n;
  
  if ( !ls )
    {
      ls=vcalloc (2, sizeof (int*));
      ns=vcalloc (2, sizeof (int));
      fl=vcalloc ((CL->S)->nseq, sizeof (int));
    }
  ls[0]=ls0;
  ls[1]=ls1;
  ns[0]=ns0; ns[1]=ns1;

  fprintf ( stderr, "\n");
  for (a=0; a<ns0; a++)
    fprintf ( stderr, "%d", ls0[a]);
  fprintf ( stderr, "\n");
  for (a=0; a<ns1; a++)
    fprintf ( stderr, "%d", ls1[a]);
	      
  ungap_sub_aln (A, ns[0], ls[0]);
  ungap_sub_aln (A, ns[1], ls[1]);
  HERE ("Align");
  pair_wise (A, ns, ls, CL);
  HERE ("Done");
	
  for (n=0,a=0;a<2; a++)
    for (b=0; b<ns[a]; b++)
      fl[n++]=ls[a][b];
  
  return sub_aln2ecl_raw_score (A,CL,n, fl);
}
//////////////////////////////////////////////////////////////////////////////
//
//                              km
///////////////////////////////////////////////////////////////////////////////
static int count;
static Sequence *KS;
static Alignment *KA;
NT_node expand_km_node (NT_node T,int n, char **name, int dim, double **V);


  
    

double **aln2km_vector (Alignment *A, char *mode, int *dim)
{
  double **v;
  int a,b,c;
  int mdim;
  int use_len=0;
  Sequence *S;
  S=aln2seq(A);
  if ( strstr (mode, "diaa"))
    {
      
      mdim=26*26;
      v=declare_double (A->nseq,mdim+4);
      for (a=0; a<A->nseq; a++)v[a]=seq2diaa(S->seq[a], v[a]);
      use_len=1;
    }
  else if ( strstr (mode, "triaa"))
    {
      mdim=26*26*26;
      v=declare_double (A->nseq,mdim+4);
      for (a=0; a<A->nseq; a++) v[a]=seq2triaa(S->seq[a], v[a]);
      use_len=1;
    }
   else if (strstr (mode, "msar"))
    {
      int start=rand()%(A->len_aln-110);
      int end=start+100;
      
      mdim=end-start+1;
      v=declare_double (A->nseq, mdim+4);
      
      for (a=0; a<A->nseq; a++)
	{
	  for (b=0,c=start; c<end;b++,c++)
	    {
	      v[a][b]=(A->seq_al[a][c]=='-')?1:0;
	    }
	}
    }
  else if (strstr (mode, "msa1"))
    {
      int window=10;
      
      mdim=A->len_aln;
      v=declare_double (A->nseq, mdim+4);
      
      for (a=0; a<A->nseq; a++)
	{
	  for (b=0; b<A->len_aln;b++)
	    {
	      v[a][b]=(A->seq_al[a][b]=='-')?1:0;
	    }
	}
    }
   else if (strstr (mode, "msa2"))
    {
      int window=10;
      
      mdim=A->len_aln;
      v=declare_double (A->nseq, mdim+4);
      
      for (a=0; a<A->nseq; a++)
	{
	  for (b=0; b<A->len_aln;b++)
	    {
	      v[a][b]=A->seq_al[a][b];
	    }
	}
    }
   else if (strstr (mode, "msa3"))
    {
      int window=10;
      
      mdim=A->len_aln;
      v=declare_double (A->nseq, mdim+4);
      
      for (a=0; a<A->nseq; a++)
	{
	  
	  for (b=0; b<A->len_aln-window;b++)
	    {
	      for (c=0; c<window; c++)
		v[a][b]+=(A->seq_al[a][b]=='-')?1:0;
	    }
	}
    }
  else if (strstr (mode, "msa4"))
    {
      mdim=A->len_aln;
      v=declare_double (A->nseq, mdim+4);
      
      for (a=0; a<A->nseq; a++)
	{
	  
	  for (b=0; b<A->len_aln;b++)
	    {
	      if (A->seq_al[a][b]=='-')c++;
	      else {v[a][b]=c;c=0;}
	    }
	}
    }
  else
    {
      return aln2km_vector(A,"triaa", dim);
    }
  //Keep only the components with a high variance
  if (!dim[0] || dim[0]>=mdim)
    {
      dim[0]=mdim;
      fprintf ( stderr, "\n\t-- Keep %d components out of %d [mode=%s]\n", dim[0], mdim, mode);
    }
  else
    {
      double **v2, tsd1,tsd2;
      int    **sd;
      sd=declare_int    (mdim,2);
      tsd1=tsd2=0;
      
      for (a=0; a<mdim; a++)
	{
	  double sum, sum2;
	  double avg;
	  sum=sum2=0;
	  
	  for (b=0; b<A->nseq; b++)
	    {
	      sum+=v[b][a];
	      sum2+=v[b][a]*v[b][a];
	    }
	  
	  avg=sum/A->nseq;
	  sd[a][0]=a;
	  sd[a][1]=100000*sqrt ((double)((sum2/A->nseq)-(avg*avg)));
	  tsd1+=(double)sd[a][1];
	}
      
      if (dim[0]<=100)
	{
	  tsd1=(tsd1*(double)dim[0])/(double)100;
	  dim[0]=0;
	  sort_int_inv (sd,2,1, 0,mdim-1);
	 
	  for (a=0; a<mdim; a++)
	    {
	      tsd2+=(double)sd[a][1];
	      if (tsd2<=tsd1){dim[0]++;}
	      else a=mdim;
	    }
	  
	  fprintf ( stderr, "\n\t-- Keep %d components out of %d [mode=%s]\n", dim[0], mdim, mode);
	  v2=declare_double (A->nseq, dim[0]+4);
	  for (a=0; a<A->nseq; a++)
	    for (b=0; b<dim[0]; b++)
	      v2[a][b]=v[a][sd[b][0]];
	  free_double (v, -1);
	  v=v2;
	}
      else
	{
	  fprintf ( stderr, "\n\t-- Keep %d components out of %d [mode=%s]\n", dim[0], mdim, mode);
	}
    }
  
  //Add the len component and scale it so that it is comparable to the other components
  if (use_len)
    {
      double SumL, SumV;
      
      SumL=SumV=0;
      
      for (c=0,a=0; a<A->nseq; a++)
	{
	  for (b=0; b<dim[0]; b++, c++)
	    SumV+=v[a][b];
	  v[a][dim[0]]=strlen(S->seq[a]);
	  SumL+= v[a][dim[0]];
	}
      SumL/=A->nseq;
      SumV/=c;
      
      for (a=0; a<A->nseq; a++)
	{
	  v[a][dim[0]]*=SumV/SumL;
	}
      dim[0]++;
    }
  
  for (a=0; a<A->nseq; a++)
    v[a][dim[0]+2]=a;
  
  return v;
}
static float tid;
static float tpairs;
static int tprint;
static int km_node;
static float km_tbootstrap;
static float km_tnode;
NT_node ** seq2km_tree (Sequence *S, char *file)
{
  int tot_node;
  NT_node T;

  Alignment *A;
  A=seq2aln(S,NULL,RM_GAP);
  T=aln2km_tree (A, "triaa", 1);
  free_aln (A);
  if (!file)file=vtmpnam (NULL);
  vfclose (print_tree (T, "newick", vfopen (file, "w")));
  
  
  return read_tree (file, &tot_node, S->nseq, S->name);
}
  
NT_node    aln2km_tree (Alignment *A, char *mode, int nboot)
{
  NT_node T;
  double **V;
  Sequence *S;
  int dim=60;//Keep all the vector components summing up to 50% of the cumulated sd
  
  KA=A;
  S=KS=aln2seq(A);
  
  if (!nboot)nboot=1;
  fprintf ( stderr, "\n-- Compute vectors\n");
  V=aln2km_vector (A, mode, &dim);
  fprintf ( stderr, "\n-- Estimate Tree (%d boostrap replicates)\n", nboot);
  T=rec_km_tree (A->name,A->nseq,dim,V, nboot);
  
  if (tprint){fprintf ( stderr, "\n---NPAIRS: %d avg id: %.2f %%\n", (int)tpairs, tid/tpairs);}
  
  if (nboot>1)fprintf (stderr, "\n-- %5d tested Nodes -- Average bootstrap: %.2f -- %d Replicates\n", (int)km_tnode, km_tbootstrap/km_tnode, nboot);
  
  return T;
}
NT_node rec_km_tree (char **name,int n,int dim,double **V, int nboot)
{
  NT_node T;
  Alignment *A;
  
  T=new_declare_tree_node ();
  if (n==1)
    {
      
      T->dist=1;
      sprintf (T->name, "%s", name[(int)V[0][dim+2]]);
      T->isseq=1;
    }
  else if (n==2)
    {
      NT_node T0,T1;
     
            
      T0=T->left=new_declare_tree_node ();
      T1=T->right=new_declare_tree_node ();
      
      T0->parent=T1->parent=T;
      T->dist=T0->dist=T1->dist=1;
          
      sprintf (T0->name, "%s", name[(int)V[0][dim+2]]);
      sprintf (T1->name, "%s", name[(int)V[1][dim+2]]);
  
      T0->isseq=1;
      T1->isseq=1;
      if (tprint)
	{
	  Alignment *A;
	  int id;
	  A=align_two_sequences (KS->seq[(int)V[0][dim+2]],KS->seq[(int)V[1][dim+2]],"pam250mt",-10,-2, "myers_miller_pair_wise");
	  id=aln2sim(A, "idmat");
	  tid+=id;
	  tpairs++;
	  fprintf ( stderr, "\nID=%4d L=%5d (%d)", id, A->len_aln,(int)tpairs);
	  free_aln (A);
	}
    }
  else
    {
      double  **V0, **V1;
      int n0,n1,a;
      int d;
      km_node++;
      
      T->bootstrap=(int)km_kmeans_bs (V,n,dim,2,0.001,NULL,nboot);
      
      T->dist=1;
      for (n0=n1=a=0; a<n; a++)(V[a][dim+1]==0)?n0++:n1++;
      fprintf ( stderr, "\t--Resolve Node %5d: (%5d .. %5d) -- %5d Support: %3d\n", km_node, n0, n1,MIN(n1,n0), (int)T->bootstrap);
      km_tbootstrap+=T->bootstrap;
      if (nboot==1)T->bootstrap=0;
      km_tnode++;
      
      V0=vcalloc (n0, sizeof (double **));
      V1=vcalloc (n1, sizeof (double **));
      for (n0=n1=a=0; a<n; a++)
	{
	  if  (V[a][dim+1]==0)V0[n0++]=V[a];
	  else V1[n1++]=V[a];
	}
      if      (n0==0 ){expand_km_node(T,n1,name, dim,V1);}
      else if (n1==0 ){expand_km_node(T,n0,name, dim,V0);}
      else 
	{
	  T->left =rec_km_tree (name, n0,dim,V0,nboot);
	  T->right=rec_km_tree (name, n1,dim,V1,nboot);
	  (T->left)->parent=(T->right)->parent=T;
	}
      vfree(V0);
      vfree(V1);
    }
  return T;
}

NT_node expand_km_node(NT_node T,int n, char ** name, int dim, double **V)
{
  NT_node Root,L, R;
  int a,b;
  
  Root=T;
   
  fprintf ( stderr, "\t**Expand  Node %5d into %5d nodes\n", km_node, n);
  km_node+=n;
  
		
  for (a=0; a<n-1; a++)
    {
      L=T->left=new_declare_tree_node ();
      L->parent=T;
      L->isseq=1;
      L->dist=0;
      sprintf (L->name, "%s", name[(int)V[a][dim+2]]);
      
      R=T->right=new_declare_tree_node ();
      R->parent=T;
      R->dist=0;
      
      T=R;
    }
  T->isseq=1;
  sprintf (T->name, "%s", name[(int)V[n-1][dim+2]]);
  return Root;
}
/*
void km_print_tree(T)
{
  if (!T->isseq)
    {
      fprintf (stdout, "(");
      km_print_tree (T->left);
      km_print_tree (T->right);
      fprintf (stdout, ")");
    }
  else
    {
      fprintf (stdout, "%s ", T->name);
    }
  return T;
}
*/
//////////////////////////////////////////////////////////////////////////////
//
//                              UPGMA
///////////////////////////////////////////////////////////////////////////////


NT_node ** int_dist2upgma_tree (int **mat, Alignment *A, int nseq, char *fname)
{
  NT_node *NL, T;
  int a, n, *used;
  int tot_node;
  if (upgma_node_heap (NULL))
    {
      printf_exit ( EXIT_FAILURE,stderr, "\nERROR: non empty heap in upgma [FATAL]");
    }
  NL=vcalloc (nseq, sizeof (NT_node));
  
  for (a=0; a<A->nseq; a++)
    {
      NL[a]=new_declare_tree_node ();
      upgma_node_heap (NL[a]);
      sprintf (NL[a]->name, "%s", A->name[a]);
      NL[a]->isseq=1;
      NL[a]->leaf=1;
    }
    
  used=vcalloc ( A->nseq, sizeof (int));
  n=A->nseq;
  while (n>1)
    {
      T=upgma_merge (mat, NL,used, &n, A->nseq);
    }    

  vfree (used);
  vfclose (print_tree (T, "newick", vfopen (fname, "w")));
  upgma_node_heap (NULL);
  vfree (NL);
  
  return read_tree (fname,&tot_node,A->nseq,  A->name);
}
NT_node upgma_merge (int **mat, NT_node *NL, int *used, int *n, int N)
{
  NT_node P, LC, RC;
  int a, b, w, best_a, best_b, set;
  float best_s;
  P=new_declare_tree_node();
  upgma_node_heap (P);
  
  for (set=0,a=0; a<N-1; a++)
    {
      if (used[a])continue;
      for (b=a+1; b<N; b++)
	{
	  if (used[b])continue;
	  w=mat[a][b];
	  if ( !set || w<best_s)
	    {
	      best_s=w;
	      best_a=a;
	      best_b=b;
	      set=1;
	    }
	}
    }

  for (a=0; a<N; a++)
    {
      if (!used[a])mat[best_a][a]=mat[a][best_a]=(mat[best_b][a]+mat[best_a][a])/2;
    }
  used[best_b]=1;
  
  LC=NL[best_a];
  RC=NL[best_b];
  best_s/=(float)100;
  LC->dist=RC->dist=best_s;
  LC->parent=RC->parent=P;
  P->left=LC;
  P->right=RC;
  NL[best_a]=P;
  n[0]--;
  return P;

}  


//////////////////////////////////////////////////////////////////////////////
//
//                              SPLIT UPGMA
///////////////////////////////////////////////////////////////////////////////

int upgma_node_heap (NT_node X)
{
  static int n;
  static NT_node *h;
  if ( X==NULL)
    {
      int a,r;
      if (n==0) return 0;
      for (a=0; a<n; a++)
	{
	  free_tree_node (h[a]);
	}
      vfree (h);
      h=NULL;
      r=n;
      n=0;
      return r;
    }
  else
    {
      h=vrealloc (h, sizeof (NT_node)*(n+1));
      h[n++]=X;
    }
  return n;
}
NT_node split2upgma_tree (Split **S, Alignment *A, int nseq, char *fname)
{
  NT_node *NL, T;
  int a, n, *used;
  
  NL=vcalloc (nseq, sizeof (NT_node));
  for (a=0; a<A->nseq; a++)
    {
      
      NL[a]=new_declare_tree_node ();
      NL[a]->lseq2=vcalloc (A->nseq+1, sizeof (int));
      NL[a]->lseq2[a]=1;
      sprintf (NL[a]->name, "%s", A->name[a]);
      NL[a]->isseq=1;
      NL[a]->leaf=1;
      NL[a]->dist=1;
      upgma_node_heap (NL[a]);
    }
  used=vcalloc ( A->nseq, sizeof (int));
  n=A->nseq;
  while (n>1)
    {
      T=split_upgma_merge (A,S, NL,used, &n, A->nseq);
    }    
  vfree (used);
  fprintf ( stdout, "\n");
  vfclose (print_tree (T, "newick", vfopen (fname, "w")));
  upgma_node_heap (NULL);
  return T;
}
NT_node split_upgma_merge (Alignment *A, Split **S, NT_node *NL, int *used, int *n, int N)
{
  NT_node P, LC, RC;
  int a, b, w, best_a, best_b, set;
  float best_s;
  static int **mat;

  if (!mat)
    {
      mat=declare_int (N, N);
      for (a=0; a<N-1; a++)
	for ( b=a+1; b<N; b++)
	  {
	    mat[a][b]=get_split_dist (A, NL[a], NL[b], S);
	  }
    }

  P=new_declare_tree_node();
  upgma_node_heap (P);
  P->lseq2=vcalloc (N, sizeof (int));
  for (set=0,a=0; a<N-1; a++)
    {
      if (used[a])continue;
      for (b=a+1; b<N; b++)
	{
	  if (used[b])continue;
	  w=mat[a][b];
	  if ( !set || w<best_s)
	    {
	      best_s=w;
	      best_a=a;
	      best_b=b;
	      set=1;
	    }
	}
    }

  
  
  
  LC=NL[best_a];
  RC=NL[best_b];
  best_s/=100;

  P->dist=1-best_s;
  LC->parent=RC->parent=P;
  P->left=LC;
  P->right=RC;
  P->bootstrap=best_s*100;
  used[best_b]=1;
  
  n[0]--;
  
  for (a=0; a<A->nseq; a++)
    {
      P->lseq2[a]=(LC->lseq2[a] || RC->lseq2[a])?1:0;
    }
  
  for (a=0; a<N; a++)
    {
      if (!used[a])mat[best_a][a]=mat[a][best_a]=(int)get_split_dist(A,P, NL[a], S);
    }
  NL[best_a]=P;
  
  return P;

}  
float get_split_dist ( Alignment *A, NT_node L, NT_node R, Split **S) 
{
  static char *split;

  
  int n,a;
  
  if (!split)
    {
      split=vcalloc ( A->nseq+1, sizeof (char));
    }
  
  
  
  for ( a=0; a<A->nseq; a++)
    split[a]=((L->lseq2[a] || R->lseq2[a])?1:0)+'0';
  
  n=0;
  while (S[n])
    {
      float score;
      if ( strm (S[n]->split,split))
	{
	  return score=100-S[n]->score;
	}
      n++;
    }
  return 100;
}
