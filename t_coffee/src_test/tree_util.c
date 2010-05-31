#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include "io_lib_header.h"
#include "util_lib_header.h"
#include "dp_lib_header.h"
#include "define_header.h"

#define TOPOLOGY 1
#define WEIGHTED 2
#define LENGTH   3
#define RECODE   4

int distance_tree;
int rooted_tree;
int tot_nseq;
static NT_node compute_fj_tree (NT_node T, Alignment *A, int limit, char *mode);
static NT_node compute_cw_tree (Alignment *A);
static NT_node compute_std_tree (Alignment *A, int n, char **arg);
static NT_node tree2fj_tree (NT_node T);
int tree_contains_duplicates (NT_node T);
int display_tree_duplicates (NT_node T);

static int compare_node1 ( int *b1, int *b2, int n);
static int compare_node2 ( int *b1, int *b2, int n);
static int find_seq_chain (Alignment *A, int **sim,int *used,int seq0,int seq1, int seq2,int chain_length, int limit, int max_chain, int *nseq);
int new_display_tree (NT_node T, int n);
NT_node display_code (NT_node T, int nseq, FILE *fp);
NT_node display_dist (NT_node T, int n, FILE *fp);
/*********************************************************************/
/*                                                                   */
/*                                   dpa_tree_manipulation           */
/*                                                                   */
/*********************************************************************/
static NT_node code_dpa_tree ( NT_node T, int **D);
NT_node collapse_sub_tree ( NT_node T,int nseq, int *list, char *new_name);
NT_node seq2dpa_tree  (Sequence *S, char *mode)
{
  Constraint_list *CL;
  NT_node **T;
  NT_node Tree;
  CL=declare_constraint_list_simple (S);
  CL->local_stderr=NULL;
  
  
  CL->DM=cl2distance_matrix (CL,NOALN,(mode==NULL)?"ktup":mode, NULL, 0);
  
  T=int_dist2nj_tree ( (CL->DM)->similarity_matrix, S->name, S->nseq, vtmpnam (NULL));
  Tree=T[3][0];
  
  Tree=recode_tree (Tree, S);
  Tree=reset_dist_tree (Tree, -1);
  
  Tree=code_dpa_tree (Tree, (CL->DM)->similarity_matrix);
  free_distance_matrix (CL->DM);
  return Tree;
}
  
NT_node tree2dpa_tree (NT_node T, Alignment *A, char *mode)
{
  /*This Function sets the branches with Length values used by DP*/
  /*The tree must be rooted*/
  Sequence *S;
  int **D;

  S=aln2seq (A);
  T=recode_tree     (T, S);
  T=reset_dist_tree (T, -1);
  D=get_sim_aln_array (A,mode);
    
  T=code_dpa_tree (T, D);
  return T;
}

NT_node code_dpa_tree ( NT_node T, int **D)
{
  if ( !T) return T;
  else if ( T->leaf==1)
    {
      T->dist=100;
      return T;
    }
  else
    {
      int nl, *ll;
      int nr, *lr;
      int a, b, min=100;
      float tot, n=0;
      
      nl=(T->left)->nseq;ll=(T->left)->lseq;
      nr=(T->right)->nseq;lr=(T->right)->lseq;

      for (tot=0,n=0, a=0; a< nl; a++)
	for ( b=0; b< nr; b++, n++)
	  {
	    tot+=D[ll[a]][lr[b]];
	    min=MIN(min,(D[ll[a]][lr[b]]));
	  }
      /*      T->dist=(mode==AVERAGE)?(tot/n):min:;*/
      T->dist=(n>0)?tot/n:0;
      T->dist=min;
      code_dpa_tree ( T->right, D);
      code_dpa_tree ( T->left, D);
      return T;
    }
}
static int group_number;
char *tree2Ngroup (Alignment *A, NT_node T, int max_n, char *fname, char *mat)
{
  double top, bot, mid, pmid;
  Sequence *S;
  int n;
  
  

  if (!T)
    {
      char **list;

      list=declare_char ( 2, 100);
      sprintf (list[0], "%s",mat); 
      
      fprintf ( stderr, "\nCompute Phylogenetic tree [Matrix=%s]", mat);
      T=compute_std_tree(A,1, list);
      fprintf ( stderr, "\nCompute dpa tree");
      T=tree2dpa_tree (T,A, mat);
    }
  
  S=tree2seq(T, NULL);
  
  if ( max_n<0)
    {
      max_n*=-1;
      n=tree2group_file (T,S,0, max_n, fname);
      fprintf ( stderr, "\n#TrimTC: Split in %d Groups at a minimum of %d%% ID\n",n, (int)max_n);
      return fname;
      
    }
  else if ( max_n>0)
    {
      if ( max_n>S->nseq)max_n=S->nseq;
            
      top=100; bot=0;
      pmid=0; mid=50;
      n=tree2group_file(T, S,0, (int)mid,fname);
      mid=dichotomy((double)n, (double)max_n,(pmid=mid), &bot, &top);
      while (n!=max_n && (int)pmid!=(int)mid)
	{
	  n=tree2group_file(T, S,0, (int)mid, fname);
	  mid=dichotomy((double)n, (double)max_n,(pmid=mid), &bot, &top);
	}
      fprintf ( stderr, "\nDONE2");
      fprintf ( stderr, "\n#TrimTC: Split in %d Groups at a minimum of %d%% ID\n",n, (int)mid);
      return fname;
    }
  return NULL;
}  
static int group_number;
int tree2group_file ( NT_node T,Sequence *S, int maxnseq, int minsim, char *name)
  {
    FILE *fp;

    
    fp=vfopen (name, "w");
    vfclose (tree2group (T, S,maxnseq,minsim, "tree2ngroup",fp));
   
    return count_n_line_in_file(name);
  }
    
      
FILE * tree2group ( NT_node T,Sequence *S, int maxnseq, int minsim,char *name, FILE *fp)
{
  if ( !T)return fp;
  else
    {
      int m,d;

      m=(maxnseq==0)?S->nseq:maxnseq;
      d=minsim;


      
      if ( T->nseq<=m && T->dist>=d)
	{
	  int a;
	  fprintf ( fp, ">%s_%d ", (name)?name:"", ++group_number);
	  for ( a=0; a< T->nseq; a++)
	    fprintf ( fp, "%s ", S->name[T->lseq[a]]);
	  fprintf (fp, "\n");
	  if (!T->parent)group_number=0;
	  return fp;
	}
      else
	{
	  fp=tree2group (T->right, S, maxnseq, minsim, name,fp);
	  fp=tree2group (T->left, S, maxnseq, minsim, name,fp);
	  if (!T->parent)group_number=0;
	  return fp;
	}
      
    }
}


NT_node  tree2collapsed_tree (NT_node T, int n, char **string)
{
  char ***list;
  Sequence *A;
  int a, *nlist;

  
  A=tree2seq(T, NULL);
  T=recode_tree(T, A);
  list=vcalloc (A->nseq, sizeof (char***));
  nlist=vcalloc (A->nseq, sizeof (int));
  if ( n==0)return T;
  else if (n>1)
    {
      int l;
      char *buf;
      
      for (l=0,a=0; a< n; a++)l+=strlen (string[a]);
      buf=vcalloc ( 2*n+l+1, sizeof (char));
      for (a=0; a< n; a++){buf=strcat (buf,string[a]), buf=strcat ( buf, " ");}     
      list[0]=string2list (buf);
      vfree (buf);
    }
  else if ( file_exists (NULL,string[0]))
    {
      list=read_group (string[0]);
    }
  else
    {
      fprintf (stderr, "\nERROR: file <%s> does not exist [FATAL:%s]\n",string[0], PROGRAM);
      myexit (EXIT_FAILURE);
    }  

  
  a=0;
  while (list[a])
    {
      int i, b;
      n=atoi (list[a][0]);
      for (b=0; b<A->nseq; b++)nlist[b]=0;
      for (b=2; b<n; b++)
	{
	  i=name_is_in_list (list[a][b], A->name, A->nseq, MAXNAMES);
	  nlist[i]=1;
	}
      T=collapse_sub_tree ( T,A->nseq,nlist,list[a][1]);
      free_char (list[a], -1);
      a++;
    }
  vfree (list);
  return T;
}

NT_node collapse_sub_tree ( NT_node T,int nseq, int *list, char *new_name)
{
  if (!T) return T;
  else
    {
      int a=0;
 

      while (a<nseq && list[a]==T->lseq2[a]){a++;} 
      if (a==nseq)
	{
	  sprintf ( T->name, "%s", new_name);
	  T->leaf=T->isseq=1;
	  T->left=T->right=NULL;
	  return T;
	}
      else
	{
	 collapse_sub_tree (T->right, nseq, list, new_name);
	 collapse_sub_tree (T->left, nseq, list, new_name);
	 return T;
	}
    }
}
  
/*********************************************************************/
/*                                                                   */
/*                                   tree pruning                    */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
NT_node remove_leaf ( NT_node T);
NT_node prune_root (NT_node T);
NT_node main_prune_tree ( NT_node T, Sequence *S)
{
  T=prune_tree ( T, S);
  return T;
}

NT_node prune_tree ( NT_node T, Sequence *S)
{
  
  if (!T ) return T;
    
  if (T->leaf && T->isseq && name_is_in_list (T->name,S->name, S->nseq, 100)==-1)
    {
      NT_node C, P, PP;
      
      P=T->parent;
      if ( !P) 
	{
	  int a;
	  for (a=0; a< S->nseq; a++)
	    {
	      HERE ("prune pb ---%s", S->name[a]);
	    }
	  myexit (EXIT_FAILURE);
	}
      C=(P->right==T)?P->left:P->right;
      PP=C->parent=P->parent;
            
      if (PP && PP->right==P)PP->right=C;
      else if (PP)PP->left=C;
      else
	{
	  if (T==P->right)P->right=NULL;
	  else P->left=NULL;
	  T=C;
	  
	}
    }
  else
    {
      prune_tree (T->left, S);
      prune_tree (T->right, S);
    }
  return prune_root(T);
}

NT_node prune_root (NT_node T)
{
  //This function prunes the root if needed (and frees it).
  if (T->parent)return T;
  
  if (!T->right && T->left)
    {
      return prune_root (T->left);
    }
  else if (T->right && !T->left)
     {
      
       return prune_root (T->right);
    }
  else
    {
      return T;
    }
}
/*********************************************************************/
/*                                                                   */
/*                                   tree comparison                 */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
int main_compare_cog_tree (NT_node T1, char *cogfile)
{
  char ***array;
  int a, nbac, n=0, p, c, b;
  Alignment *A;
  
  array=file2list(cogfile, ";\n");
  nbac=atoi(array[0][0])-2;
  
  A=declare_aln2 (nbac+1, 10);
  for (a=0; a<nbac; a++)
    {
      sprintf ( A->name[a], "%s", array[0][a+2]);
      A->seq_al[a][0]='a';
      A->seq_al[a][1]='\0';
      
    }
  sprintf ( A->name[nbac], "cons");
  
  A->nseq=nbac+1;
  A->len_aln=1;
  
  
  n=3;
  while (array[n]!=NULL)
    {
      for (b=0; b<nbac; b++)
	{
	  p=atoi (array[1][b+2]);p=(p==0)?'O':'I';
	  c=atoi (array[n][b+2]);c=(c==0)?'O':'I';
	  A->seq_al[b][0]=p;
	  A->seq_al[b][1]=c;
	  A->seq_al[b][2]='\0';
	 
	}
      sprintf (A->file[0], "%s", array[n][1]);
      A->len_aln=2;
      main_compare_aln_tree (T1, A, stdout);
      n++;
    }
  return n;
}

      
int main_compare_aln_tree (NT_node T1, Alignment *A, FILE *fp)
{
  int n=0;
  
  fprintf ( fp, "\nTOT_CLASH COG %s N %d", A->file[0], compare_aln_tree (T1, A, &n, fp));
  vfclose (fp);
  return n;
}

int compare_aln_tree (NT_node T, Alignment *A, int *n, FILE *fp)
{
  if (T->leaf)
    {
      int i;
      i=name_is_in_list (T->name, A->name, A->nseq, 100);
      T->seqal=A->seq_al[i];
      return 0;
    }
  else
    {
      char *seq1, *seq2;
      if (!(T->left )->seqal)compare_aln_tree (T->left, A,n, fp);
      if (!(T->right)->seqal)compare_aln_tree (T->right, A,n, fp);
      
      seq1=(T->left)->seqal;
      seq2=(T->right)->seqal;
      (T->left)->seqal=(T->right)->seqal=NULL;
      if ( seq1 && seq2)
	{
	  if (strm (seq1, seq2))
	    {
	      T->seqal=seq1;
	      
	    }
	  else
	    {
	   
	      if (seq1[0]!=seq2[0] && seq1[1]!=seq2[1])
		{
		  fprintf ( fp, "\nNODE_CLASH: COG %s (%s,%s):(",A->file[0],seq1,seq2 );
		  display_leaf_below_node (T->left, fp);
		  fprintf ( fp, ");("); 
		  display_leaf_below_node (T->right, fp);
		  fprintf ( fp, ")"); 
		  n[0]++;
		}
	    }
	}
    }
  return n[0];
}
//**********************************************************************

int compare_split (int *s1, int *s2, int l);
int get_split_size (int *s, int l);

int main_compare_splits ( NT_node T1, NT_node T2, char *mode,FILE *fp)
{
  Sequence *S1, *S2, *S;
  int  a, b;

  

  int **sl1, n1;
  int **sl2, n2;
  if ( tree_contains_duplicates (T1))
    {
      display_tree_duplicates (T1);
      printf_exit (EXIT_FAILURE, stderr, "\nFirst Tree Contains Duplicated Sequences [main_compare_trees][FATAL:%s]", PROGRAM);
     
    }
  else if ( tree_contains_duplicates (T2))
    {
      display_tree_duplicates (T2);
      printf_exit (EXIT_FAILURE, stderr, "\nSecond Tree Contains Duplicated Sequences [main_compare_trees]");
      
    }
  
  //Identify the commom Sequence Set  
  S1=tree2seq(T1, NULL);
 
  
  S2=tree2seq(T2, NULL);
 
  
  S=trim_seq ( S1, S2);
  
  //Prune the trees and recode the subtree list
  T1=prune_tree (T1, S);
  T1=recode_tree(T1, S);
  
  T2=prune_tree (T2, S);  
  T2=recode_tree(T2, S);
  HERE ("1");
  sl1=declare_int (10*S->nseq, S->nseq);
  sl2=declare_int (10*S->nseq, S->nseq);
  
  HERE ("2");
  n1=n2=0;
  tree2split_list (T1, S->nseq, sl1, &n1);
  tree2split_list (T2, S->nseq, sl2, &n2);
  
  for (a=0; a<n1; a++)
    {
      int n, best, s;
      
      n=get_split_size (sl1[a], S->nseq);
      for (best=0,b=0; b<n2; b++)
	{
	  s=compare_split (sl1[a], sl2[b], S->nseq);
	  best=MAX(s,best);
	}
      fprintf ( fp, "\n%4d %4d ", MIN(n,(S->nseq)), best);
      for (b=0; b<S->nseq; b++)fprintf ( fp, "%d", sl1[a][b]);
    }
      
  free_sequence (S, -1);
  free_sequence (S1, -1);
  free_sequence (S2, -1);
  myexit (EXIT_SUCCESS);
  return 1;
}  
int compare_split (int *s1, int *s2, int l)
{
  int n1, n2, score1, score2, a;
  n1=get_split_size (s1, l);
  n2=get_split_size (s2, l);

  for (score1=0,a=0; a< l; a++)
    {
      score1+=(s1[a]==1 && s2[a]==1)?1:0;
    }
  score1=(score1*200)/(n1+n2);
  
  for ( score2=0, a=0; a<l; a++)
    {
      score2+=(s1[a]==0 && s2[a]==1)?1:0;
    }
  score2=(score2*200)/((l-n1)+n2);
  return MAX(score1, score2);
}
int get_split_size (int *s, int l)
{
  int a, b;
  for (a=b=0; a<l; a++)b+=s[a];
  return b;
}
//**********************************************************************
  //
  //
  //                       TREEE COMPARISON FUNCTIONS
  //
  //
  //
  //////////////////////////////////////////////////////////////////////

//JM_START
void normalizeScore(float *score, int len)
{
	int i;
	float SCORE_MIN = FLT_MAX;
	float SCORE_MAX = FLT_MIN;
	for(i = 0; i < len; i++)
	{
		if(score[i] < SCORE_MIN)
			SCORE_MIN = score[i];
		if(score[i] > SCORE_MAX)
			SCORE_MAX = score[i];
	}
	for(i = 0; i < len; i++)
		score[i] = (9*(score[i]-SCORE_MIN)/(SCORE_MAX-SCORE_MIN));
}

int new_compare_trees ( NT_node T1, NT_node T2, int nseq, Tree_sim *TS);
NT_node new_search_split (NT_node T, NT_node B, int nseq);
int new_compare_split ( int *b1, int *b2, int n);
Tree_sim*  tree_scan_pos (Alignment *A, int start, int end, char *ptree, NT_node RT);
Tree_sim*  tree_scan_pos_woble (Alignment *A, int center, int max, char *ptree, NT_node RT, int *br, int *bl );
Tree_sim*  tree_scan_pair_pos (Alignment *A, int start, int end, int start2, int end2,char *ptree, NT_node RT);
Tree_sim*  tree_scan_multiple_pos (int *poslist, int *wlist,int nl, Alignment *A, char *ptree, NT_node RT);
NT_node aln2std_tree(Alignment *A, int ipara1, int ipara2, char *mode);

NT_node tree_scan (Alignment *A,NT_node RT, char *pscan, char *ptree)
{
  int l, a,ax, c, cx, b;
  char mode[100];
  int start, w;
  int nl, *poslist;
  char posfile[100];
  
  char *pcFileName = A->file[0];
  char prefix[200] ={0};
  int len = (strrchr(pcFileName,'.')?strrchr(pcFileName,'.')-pcFileName:strlen(pcFileName));
  strncpy(prefix, pcFileName, len);

  float *fascore;
  char out_format[100];

   char *score_csv_file = vcalloc(200, sizeof (char));
   char *score_html_file = vcalloc(200, sizeof (char));
   char *hit_matrix_file = vcalloc(200, sizeof (char));
   char *hit_html_file = vcalloc(200, sizeof (char));
   char *tree_file = vcalloc(200, sizeof (char));

  sprintf(score_csv_file, "%s%s", prefix, ".score_csv");
  sprintf(score_html_file, "%s%s", prefix, ".ts_html");
  sprintf(hit_matrix_file, "%s%s", prefix, ".hit_matrix");
  sprintf(hit_html_file, "%s%s", prefix, ".hit_html");
  sprintf(tree_file, "%s%s", prefix, ".trees_txt");

  if ( pscan && strstr ( pscan, "help"))
    {
      fprintf ( stdout, "\n+tree_scan| _W_     : Window size for the tree computation|STD size in norscan mode");
      fprintf ( stdout, "\n+tree_scan| _MODE_  : Mode for the number of windows (single, double, list, scan, pairscan, norscan, hit, norhit)"); 
      fprintf ( stdout, "\n+tree_scan| _MINW_  : Minimum Window size when using the scan mode (4)");
      fprintf ( stdout, "\n+tree_scan| _OUTTREE_ : specify the format of outputing tree in every position (default: not ouput)");
      myexit (EXIT_SUCCESS);
    }
 
  strget_param (pscan, "_W_", "5", "%d",&w);
  strget_param (pscan, "_MODE_", "single", "%s",mode);
  strget_param (pscan, "_MINW_", "1", "%d",&start);
  strget_param (pscan, "_POSFILE_", "NO", "%s", posfile);
  strget_param (pscan, "_OUTTREE_", "", "%s", &out_format);

	if(strlen(out_format) > 1)
		unlink(tree_file);	

  l=intlen (A->len_aln);
  
  poslist=vcalloc ( A->len_aln, sizeof (int));
  nl=0;
  fascore = vcalloc(A->len_aln, sizeof (float));

  if ( strm (posfile, "NO"))
    {
     
      for ( a=0; a< A->len_aln; a++)poslist[nl++]=a+1;
    }
  else
    {
      int *p; 
      p=file2pos_list (A,posfile);
      poslist=pos2list (p, A->len_aln, &nl);
      for (a=0; a<nl; a++)poslist[a]++;
      vfree (p);
    }
 
//For tree hit
  NT_node *TreeArray = vcalloc(nl, sizeof (NT_node));

  if ( strm (mode, "woble"))
    {
      for  (ax=0; ax<nl; ax++)
	{
	  Tree_sim *TS;
	  int left=0, right=0;
	  a=poslist[ax];
	  
	  TS=tree_scan_pos_woble (A,a,w,ptree, RT, &left, &right);
	  fprintf ( stdout, "P: %*d I: %*d %*d SIM: %6.2f\n", l,a,l,left,l,right,TS->uw);
	  vfree (TS);
	}
    }
  else if ( strm (mode, "single"))
    {
      for (ax=0; ax<nl; ax++)
	{
	  Tree_sim *TS;
	  int pstart, pend;
	  
	  a=poslist[ax];
	  
	  pstart=a-b;
	  pend=a+b;
	  
	  if (pstart<1 || pstart>A->len_aln)continue;
	  if (pend<1 || pend>A->len_aln)continue;
	  TS=tree_scan_pos (A, pstart,pend, ptree, RT);
	  fprintf ( stdout, "P: %*d I: %*d %*d SIM: %6.2f L: %2d\n", l,a,l,pstart,l,pend,TS->uw, (w*2)+1);
	  vfree (TS);
	}
    }
  else if (strm (mode, "scan")||strm (mode, "hit"))
    {      
    	FILE *fp_ts;
    	fp_ts=vfopen (score_csv_file, "w");	
      fprintf ( fp_ts, "Position,Win_Beg,Win_End,Similarity,Win_Len\n");
      for ( ax=0; ax<nl; ax++)
	{
	    float best_score=0;
	    int best_pos=0, best_w=0, best_start, best_end;
	    
	    a=poslist[ax];
	    best_pos = best_start = best_end = a;
	    for (b=start; b<=w; b++)
	      {
		Tree_sim *TS;
		int pstart, pend;
		pstart=a-b;
		pend=a+b;
		
		if (pstart<1 || pstart>A->len_aln)continue;
		if (pend<1 || pend>A->len_aln)continue;
		TS=tree_scan_pos (A, pstart,pend, ptree, RT);
		if (TS->uw>=best_score)
			{best_score=TS->uw;best_w=b;best_start=pstart; best_end=pend;}
		vfree (TS);
	      }
	    fprintf (fp_ts, "%*d,%*d,%*d,%6.2f,%2d\n", l,best_pos, l,best_start, l,best_end, best_score,(best_w*2)+1);
	    fascore[ax]=(float)best_score;
		if(strlen(out_format) > 1)
			vfclose (print_tree (aln2std_tree(A, best_start, best_end, mode), out_format, vfopen (tree_file, "a+")));
		if(strm (mode, "hit"))
			TreeArray[ax] = aln2std_tree(A, best_start, best_end, mode);
	  }
	vfclose(fp_ts);
    }
//tree scan by using normal distribution window
//or
//generate hit matrix
  else if ( strm (mode, "norscan")||strm (mode, "norhit"))
    {
      FILE *fp_ts;
      ptree=vcalloc(100, sizeof (char));
      fp_ts=vfopen (score_csv_file, "w");	
      fprintf ( fp_ts, "Position,Similarity,STD_Len\n");
      for ( ax=0; ax<nl; ax++)
	  {
	    float best_score=DBL_MIN;
	    int best_STD = start;
	    a=poslist[ax];
	    for (b=start; b<=w; b++)
	      {
		Tree_sim *TS;
		sprintf ( ptree, "+aln2tree _COMPARE_nordisaln__STD_%d__CENTER_%d_", b, a); //should be used a or ax
		TS=tree_scan_pos (A, 1, nl, ptree, RT);
		if (TS->uw>=best_score)
			{best_score=TS->uw;best_STD=b;}
		vfree (TS);		
	      }
	      fascore[ax]=best_score;
	      fprintf ( fp_ts, "%*d,%6.2f,%d\n", l,a, fascore[ax], best_STD);
		if(strlen(out_format) > 1)
			vfclose (print_tree (aln2std_tree(A, best_STD, a, mode), out_format, vfopen (tree_file, "a+")));
		if(strm (mode, "norhit"))
			TreeArray[ax] = aln2std_tree(A, best_STD, a, mode);
	  }
	vfclose(fp_ts);
    }
//generate hit matrix
	if (strm (mode, "hit")||strm (mode, "norhit"))
    {
//Compute the pair score of tree scan segqtion
	fprintf (stdout, "[STRAT] Calculate the hit matrix of the tree scan\n");
	float **ffpHitScoreMatrix;

	ffpHitScoreMatrix=vcalloc (nl, sizeof (float*));
	int i, j;
	for(i = 0; i < nl; i++)
      		ffpHitScoreMatrix[i]=vcalloc (nl-i, sizeof (float));

	fprintf (stdout, "Process positions\n", i);
	for(i = 0; i < nl; i++)
	{
		fprintf (stdout, "%d, ", i);
		for(j = i; j < nl; j++)
		{
			Tree_sim *TS;
			TS=tree_cmp (TreeArray[i], TreeArray[j]);
			ffpHitScoreMatrix[i][j-i] = TS->uw;
			vfree (TS);
		}
	}
	vfree(TreeArray);
	fprintf (stdout, "\n");
	output_hit_matrix(hit_matrix_file, ffpHitScoreMatrix, nl);
	fprintf (stdout, "[END]Calculate the hit matrix of the tree scan\n");

//Output Hit Score into color html
	output_hit_color_html  (A, ffpHitScoreMatrix, nl, hit_html_file);
	vfree(ffpHitScoreMatrix);
    }
  else if ( strm (mode, "pairscan"))
    {
      int d, set;
     
      for ( ax=0; ax<nl; ax++)
	  {
	    float best_score=0;
	    int best_pos=0, best_w=0, best_w2=0, best_start, best_end, best_pos2=0, best_start2, best_end2;

	    Tree_sim *TS;
	    int pstart, pend, p2start, p2end;
	    a=poslist[ax];
	    for ( cx=0; cx<nl; cx++)
	      {
		set=0;
		c=poslist[cx];
		for (d=start; d<=w; d++)
		  {
		    for (b=start; b<=w; b++)
		      {
			pstart=a-b;
			pend=a+b;
			p2start=c-d;
			p2end=c+d;
			if (pstart<1 || pstart>A->len_aln)continue;
			if (pend<1 || pend>A->len_aln)continue;
			if (p2start<1 || p2start>A->len_aln)continue;
			if (p2end<1 || p2end>A->len_aln)continue;
			if (pstart<=p2start && pend>=p2start) continue;
			if (pstart<=p2end && pend>=p2end) continue;
			TS=tree_scan_pair_pos (A, pstart,pend,p2start, p2end, ptree, RT);

			if (TS->uw>=best_score){best_score=TS->uw; best_pos=a;best_w=b;best_start=pstart; best_end=pend; best_pos2=c, best_w2=d, best_start2=p2start, best_end2=p2end;set=1;}
			vfree (TS);
		      }
		  }
		if (set)fprintf ( stdout, "P1: %*d  I1: %*d %*d P2: %*d I2: %*d %*d SIM: %6.2f L: %2d\n", l,best_pos, l,best_start, l,best_end, l, best_pos2, l, best_start2, l, best_end2, best_score,(best_w*2)+1 );
		
		set=0;
	      }
	  }
    }
  else if ( strm (mode, "multiplescan"))
    {
      int n, **wlist, best_pos;
      float best_score;
      Tree_sim *TS;
      wlist=generate_array_int_list (nl*2,start, w,1, &n, NULL);
      HERE ("Scan %d Possibilities", n);
      
      for (best_score=best_pos=0,a=0; a<n; a++)
	{
	  TS=tree_scan_multiple_pos (poslist,wlist[a],nl, A, ptree, RT);
	  if (TS && TS->uw>best_score)
	    {
	      best_score=TS->uw;
	      fprintf ( stdout, "\n");
	      for (b=0; b<nl; b++)
		{
		  fprintf ( stdout, "[%3d %3d %3d]", poslist[b], wlist[a][b*2], wlist[a][b*2+1]);
		}
	      fprintf ( stdout, " SCORE: %.2f", best_score);
	    }
	  if (TS)vfree (TS);
	}
    }

//Output Tree Scan core into color html
    normalizeScore(fascore, nl);
    Alignment *ST;
    Sequence *S;
    int i, r1;

    S=A->S;
    ST=copy_aln (A, NULL);
    for (a=0; a<ST->nseq; a++)
    {
	i=name_is_in_list (ST->name[a],S->name, S->nseq, 100);
	if ( i!=-1)
	{
		for (b=0; b<ST->len_aln; b++)
		{
			r1=ST->seq_al[a][b];
			if ( r1!='-')
				r1 = (int)fascore[b] + 48;
			ST->seq_al[a][b]=r1;
		}
	}
    }	    
    output_color_html  ( A, ST, score_html_file);

//free memory
	free_aln(ST);
	vfree(fascore);
  vfree(score_csv_file);
  vfree(score_html_file);
  vfree(hit_matrix_file);
  vfree(hit_html_file);

  myexit(EXIT_SUCCESS);
}

NT_node aln2std_tree(Alignment *A, int ipara1, int ipara2, char *mode)
{
     Alignment *B;
     NT_node T;
     char *cpSet = vcalloc(100, sizeof (char));

	if(strm (mode, "norhit"))    
	{
		B=extract_aln (A, 1, A->len_aln);
		sprintf ( cpSet, "+aln2tree _COMPARE_nordisaln__STD_%d__CENTER_%d_", ipara1, ipara2);
	}
	else
     		B=extract_aln (A, ipara1, ipara2);

	T=compute_std_tree (B, (cpSet)?1:0, (cpSet)?&cpSet:NULL);
	free_aln(B);
	return T;
}

Tree_sim*  tree_scan_multiple_pos (int *poslist, int *wlist,int nl, Alignment *A, char *ptree, NT_node RT)
{
    static Alignment *B;
    static int *pos;
    int a, b, n, s, p, left, right;
    Tree_sim *TS;
    NT_node T=NULL;
    

    //poslist positions come [1..n]
    vfree(pos);
    free_aln (B);

    pos=vcalloc ( A->len_aln+1, sizeof (int));
    B=copy_aln (A, NULL);
    
    for (a=0; a<nl; a++)
      {
	p=poslist[a]-1;
	left =wlist[2*a];
	right=wlist[2*a+1];
	for (b=p-left; b<=p+right; b++)
	  {
	    if (b<1 ||b>A->len_aln) return NULL;
	    else pos[b]++;
	    
	    if (pos[b]>1) return NULL;
	  }
      }
    
    for (s=0; s<A->nseq; s++)
      {
	for (n=0,a=1; a<=A->len_aln; a++)
	  {
	    if (pos[a])B->seq_al[s][n++]=A->seq_al[s][a-1];
	  }
      }
    
    B->len_aln=n;
    for (s=0; s<A->nseq; s++)B->seq_al[s][B->len_aln]='\0';
   
    T=compute_std_tree (B, (ptree)?1:0, (ptree)?&ptree:NULL);
    
    TS=tree_cmp (T, RT);
    
    free_tree(T);
    return TS;
  }
  
Tree_sim*  tree_scan_pair_pos (Alignment *A, int start, int end, int start2, int end2,char *ptree, NT_node RT)
  {
    Tree_sim *TS;
    Alignment *B,*B1, *B2;
    NT_node T=NULL;
    int a;
    
 
    B=copy_aln (A, NULL);
     B1=extract_aln (A,start,end);
     B2=extract_aln (A,start2, end2);

    
     for ( a=0; a< B->nseq;a++)
       sprintf (B->seq_al[a], "%s%s", B1->seq_al[a], B2->seq_al[a]);
     B->len_aln=strlen (B->seq_al[0]);
     
     T=compute_std_tree (B, (ptree)?1:0, (ptree)?&ptree:NULL);
     TS=tree_cmp (T, RT);
     
     free_tree(T);
     free_aln (B);free_aln(B1); free_aln(B2);
     
     return TS;
  }

Tree_sim*  tree_scan_pos (Alignment *A, int start, int end, char *ptree, NT_node RT)
  {
     Tree_sim *TS;
     Alignment *B;
     NT_node T;
     
     if ( start<1 || start>A->len_aln) return NULL;
     if ( end<1 || end>A->len_aln) return NULL;

     B=extract_aln (A,start,end);
     T=compute_std_tree (B, (ptree)?1:0, (ptree)?&ptree:NULL);
     TS=tree_cmp (T, RT);
     free_tree(T);free_aln (B);
     return TS;
  }
Tree_sim*  tree_scan_pos_woble (Alignment *A, int center, int max, char *ptree, NT_node RT, int *br, int *bl )
  {
    Tree_sim *TS,*BTS;


     int left, right;
     float best_score=0;
     int start, end;
     
     br[0]=bl[0]=0;
     BTS=vcalloc (1, sizeof (Tree_sim));
     
     for (left=0; left<max; left++)
       for (right=0; right<max; right++)
	 {
	   start=center-left;
	   end=center+right;
	   TS=tree_scan_pos (A, start, end, ptree, RT);
	   
	   if (TS && TS->uw >best_score)
	     {
	       best_score=TS->uw;
	       BTS[0]=TS[0];
	       br[0]=right; bl[0]=left;
	       vfree(TS);
	     }
	 }
     return BTS;
  }

Tree_sim* tree_cmp( NT_node T1, NT_node T2)
{
  Sequence *S1, *S2, *S;
  int n;
  int a;
  
  Tree_sim *TS1, *TS2;

  if ( tree_contains_duplicates (T1))
    {
      display_tree_duplicates (T1);
      printf_exit (EXIT_FAILURE, stderr, "\nFirst Tree Contains Duplicated Sequences [main_compare_trees][FATAL:%s]", PROGRAM);
     
    }
  else if ( tree_contains_duplicates (T2))
    {
      display_tree_duplicates (T2);
      printf_exit (EXIT_FAILURE, stderr, "\nSecond Tree Contains Duplicated Sequences [main_compare_trees]");
      
    }
  
  //Identify the commom Sequence Set  
  S1=tree2seq(T1, NULL);
  S2=tree2seq(T2, NULL);

  S=trim_seq ( S1, S2);

  if ( S->nseq<=2)
    {
      fprintf ( stderr, "\nERROR: Your two do not have enough common leaf to be compared [FATAL:PROGRAM]");
    }

  //Prune the trees and recode the subtree list
  T1=prune_tree (T1, S);
  T1=recode_tree(T1, S);
  
  T2=prune_tree (T2, S);  
  T2=recode_tree(T2, S);

  TS1=vcalloc (1, sizeof (Tree_sim));
  TS2=vcalloc (1, sizeof (Tree_sim));

  new_compare_trees ( T1, T2, S->nseq, TS1);
  new_compare_trees ( T2, T1, S->nseq, TS2);
  

  TS1->n=tree2nnode (T1);
  TS1->nseq=S->nseq;
  
  TS2->n=tree2nnode (T2);
  /*if (TS1->n !=TS2->n)
    printf_exit (EXIT_FAILURE, stderr,"\nERROR: Different number of Nodes in the two provided trees after prunning [FATAL: %s]", PROGRAM);
  */
  
  free_sequence (S, -1);
  free_sequence (S1, -1);
  free_sequence (S2, -1);

  TS1->uw=(TS1->uw+TS2->uw)*100/(TS1->max_uw+TS2->max_uw);
  TS1->w=(TS1->w+TS2->w)*100/(TS1->max_w+TS2->max_w);
  TS1->d=(TS1->d+TS2->d)*100/(TS1->max_d+TS2->max_d);
  TS1->rf=(TS1->rf+TS2->rf)/2;
  vfree (TS2);
  return TS1;
}

NT_node main_compare_trees ( NT_node T1, NT_node T2, FILE *fp)
{
  Tree_sim *T;
  
  T=tree_cmp (T1, T2);
  fprintf ( fp, "\n#tree_cmp|T: %.f W: %.2f L: %.2f RF: %d N: %d S: %d", T->uw, T->w, T->d, T->rf, T->n, T->nseq);
  fprintf ( fp, "\n#tree_cmp_def|T: ratio of identical nodes");
  fprintf ( fp, "\n#tree_cmp_def|W: ratio of identical nodes weighted with the min Nseq below node");
  fprintf ( fp, "\n#tree_cmp_def|L: average branch length similarity");
  fprintf ( fp, "\n#tree_cmp_def|RF: Robinson and Foulds");
  fprintf ( fp, "\n#tree_cmp_def|N: number of Nodes in T1 [unrooted]");
  fprintf ( fp, "\n#tree_cmp_def|S: number of Sequences in T1\n");
  
  vfree (T);
  return T1;
}  

int new_compare_trees ( NT_node T1, NT_node T2, int nseq, Tree_sim *TS)
{
  int n=0;
  NT_node N;
  float t1, t2;
  
  if (!T1 || !T2) return 0;
  
  n+=new_compare_trees (T1->left,  T2, nseq,TS);
  n+=new_compare_trees (T1->right, T2, nseq,TS);
  
  //Exclude arbitrary splits (dist==0)
  if ((T1->dist==0) && !(T1->parent))return n;
  
  N=new_search_split (T1, T2, nseq);
  t1=FABS(T1->dist);
  t2=(N)?FABS(N->dist):0;
  TS->max_d+=MAX(t1, t2);

  if (!N)TS->rf++;
  if (T1->nseq>1)
    {
      int w;
      w=MIN((nseq-T1->nseq),T1->nseq);
      TS->max_uw++;
      TS->max_w+=w;
      
      if (N)
	{
	  TS->uw++;
	  TS->w+=w;
	  TS->d+=MIN(t1, t2);
	  //T1->dist=T1->nseq;
	}
      else
	{
	  //T1->dist=T1->nseq*-1;
	  ;
	}
    }
  else
    {
      TS->d+=MIN(t1, t2);
      //T1->dist=1;
    }
  return ++n;
}
NT_node new_search_split (NT_node T, NT_node B, int nseq)
{
  NT_node N;
  if (!T || !B) return NULL;
  else if ( new_compare_split (T->lseq2, B->lseq2, nseq)==1)return B;
  else if ( (N=new_search_split (T, B->right, nseq)))return N;
  else return new_search_split (T, B->left, nseq);
}
int new_compare_split ( int *b1, int *b2, int n)
{
  int a, flag;
  
  for (flag=1, a=0; a<n; a++)
    if (b1[a]!=b2[a])flag=0;
  if (flag) return flag;
  
  for (flag=1, a=0; a<n; a++)
    if (b1[a]==b2[a])flag=0;
  return flag;
}
  
float compare_trees ( NT_node T1, NT_node T2, int nseq,int  mode)
{
  /*search each branch of T1 in T2*/
  float n=0;
  
  
  if ( !T1 || !T2)return 0;
  
  if (getenv4debug("DEBUG_TREE_COMPARE"))display_node (T1, "\nNODE ", nseq);
  
  if (T1->parent && T1->nseq>1)n+=search_node ( T1, T2, nseq, mode);
    
  n+=compare_trees ( T1->left, T2, nseq, mode);
  n+=compare_trees ( T1->right, T2, nseq, mode);
  
  return n;
}

float search_node ( NT_node B, NT_node T, int nseq, int mode)
{
  int n=0;
  if ( !B || !T) return -1;
  if (getenv4debug("DEBUG_TREE_COMPARE"))display_node ( T, "\n\t", nseq);
  
  n=compare_node ( B->lseq2, T->lseq2, nseq );
  
  if ( n==1) 
    {
      if (getenv4debug("DEBUG_TREE_COMPARE"))fprintf ( stderr, "[1][%d]", (int)evaluate_node_similarity ( B, T, nseq, mode));
      if (mode==RECODE)B->dist=B->leaf;
      return evaluate_node_similarity ( B, T, nseq, mode);
    }
  else if ( n==-1)
    {
      if (getenv4debug("DEBUG_TREE_COMPARE"))fprintf ( stderr, "[-1]"); 
      if (mode==RECODE)B->dist=-B->leaf;
      return 0;
    }
  else
    {
      if (getenv4debug("DEBUG_TREE_COMPARE"))fprintf ( stderr, "[0]");
      n=search_node ( B, T->left, nseq, mode);   
      if ( n>0) return n;
      n=search_node ( B, T->right, nseq, mode);
      if ( n>0) return n;
      n=search_node ( B, T->bot, nseq, mode);
      if ( n>0) return n;
    }
  return n;
}

float evaluate_node_similarity ( NT_node B, NT_node T, int nseq, int mode)
{
int a, c;
 
 if ( mode==TOPOLOGY || mode ==RECODE)
   {
     for ( a=0; a< nseq; a++)
       if ( B->lseq2[a]!=T->lseq2[a]) return 0;
     return 1;
   }
 else if ( mode == WEIGHTED)
   {
     for (c=0, a=0; a< nseq; a++)
       {
	 if ( B->lseq2[a]!=T->lseq2[a]) return 0; 
	 else c+=B->lseq2[a];
       }      
     return (float)(MIN(c,nseq));
   }
 else if ( mode == LENGTH )
   {
     float d1, d2;
     
     for (c=0, a=0; a< nseq; a++)
       {
	 if ( B->lseq2[a]!=T->lseq2[a]) return 0; 
       }
     d1=FABS((B->dist-T->dist));
     d2=MAX(B->dist, T->dist);
     return (d2>0)?(d1*100)/d2:0;
   }
 else
   {
     return 0;
   }
}
int compare_node ( int *b1, int *b2, int nseq)
{
  int n1, n2;
  
  n1=compare_node1 ( b1, b2, nseq);
  /*fprintf ( stderr, "[%d]", n1);*/
  if ( n1==1) return 1;
  
  n2=compare_node2 ( b1, b2, nseq);
  /* fprintf ( stderr, "[%d]", n2);*/
  if ( n2==1)return 1;
  else if ( n2==-1 && n1==-1) return -1;
  else return 0;
}
int compare_node1 ( int *b1, int *b2, int n)
{
  int a;
  int l1, l2;
  int r=1;
  for ( a=0; a< n; a++)
    {
      l1=b1[a];
      l2=b2[a];
      if ( l1==1 && l2==0) return -1;
      if ( l1!=l2)r=0;
    }
  return r;
}
int compare_node2 ( int *b1, int *b2, int n)
{
  int a;
  int l1, l2;
  int r=1;
  
  for ( a=0; a< n; a++)
    {
      l1=1-b1[a];
      l2=b2[a];
      if ( l1==1 && l2==0) return -1;
      if ( l1!=l2) r=0;
    }
  return r;
}
void display_node (NT_node N, char *string,int nseq)
{
  int a;
  fprintf ( stderr, "%s", string);
  for (a=0; a< nseq; a++)fprintf ( stderr, "%d", N->lseq2[a]);
}


/*********************************************************************/
/*                                                                   */
/*                                   FJ_tree Computation             */
/*                                                                   */
/*                                                                   */
/*********************************************************************/


NT_node tree_compute ( Alignment *A, int n, char ** arg_list)
{
  if (n==0 || strm (arg_list[0], "cw"))
    {
      return compute_cw_tree (A);
    }
  else if ( strm (arg_list[0], "fj"))
    {
      return compute_fj_tree ( NULL, A, (n>=1)?atoi(arg_list[1]):8, (n>=2)?arg_list[2]:NULL);
    }
  
  else if ( ( strm (arg_list[0], "nj")))
    {
      return compute_std_tree (A, n, arg_list);
    }
  else
    return compute_std_tree (A, n, arg_list);
}

NT_node compute_std_tree (Alignment *A, int n, char **arg_list)
{
  return compute_std_tree_2 (A, NULL, list2string (arg_list, n));
}
NT_node compute_std_tree_2 (Alignment *A, int **s, char *cl)
{
  
  NT_node T, **BT=NULL;
  char *tree_name;

  char matrix[100];
  char score [100];
  char compare[100];
  char tmode[100];
  int free_s=0;
  
  tree_name =vtmpnam (NULL);
  
  if (strstr (cl, "help"))
    {
      fprintf ( stdout, "\n+aln2tree| _MATRIX_ : matrix used for the comparison (idmat, sarmat, pam250mt..)\n");
      fprintf ( stdout, "\n+aln2tree| _SCORE_  :  score mode used for the distance (sim, raw)\n");
      fprintf ( stdout, "\n+aln2tree| _COMPARE_:  comparison mode (aln, ktup, align, nordisaln)\n");
      fprintf ( stdout, "\n+aln2tree| _TMODE_  :  tree mode (nj, upgma)\n");
      myexit (EXIT_SUCCESS);
    }
      

  //matrix: idmat, ktup,sarmat, sarmat2 
  strget_param (cl, "_MATRIX_", "idmat", "%s",matrix);
  
  //score: sim, raw
  strget_param (cl, "_SCORE_", "sim", "%s",score);
  
  //compare: aln, ktup, align
  strget_param (cl, "_COMPARE_", "aln", "%s",compare);

  //compare: aln, ktup, align
  strget_param (cl, "_TMODE_", "nj", "%s",tmode);

  int STD, CENTER;
  if ( strm (compare, "nordisaln"))
  {
  	strget_param (cl, "_STD_", "1", "%d", &STD);
  	strget_param (cl, "_CENTER_", "5", "%d", &CENTER);
  }
  //Use external msa2tree methods
  if ( strm (tmode, "cw"))
    {
      free_int (s, -1);
      return compute_cw_tree (A);
    }
  
  //compute distance matrix if needed
  if ( !s)
    {
      free_s=1;
    if ( strm (compare, "ktup"))
	{
	  ungap_array (A->seq_al, A->nseq);
	  s=get_sim_aln_array ( A,cl);
	}
      else if ( strm ( compare, "aln"))
	{
	  if (strm (score, "sim"))
	    s=get_sim_aln_array(A, matrix);
	  else if ( strm (score, "raw"))
	    {
	      s=get_raw_sim_aln_array (A,matrix);
	    }
	}
      else if ( strm ( compare, "nordisaln"))
	{
	  s=get_sim_aln_array_normal_distribution(A, matrix, &STD, &CENTER);
	}
      s=sim_array2dist_array(s, 100);
    }

  //Compute the tree
  if (strm (tmode, "nj"))
    {
     
      BT=int_dist2nj_tree (s, A->name, A->nseq, tree_name);
      T=main_read_tree (tree_name);
      free_read_tree(BT);
    }
  else if (strm (tmode, "upgma"))
    {
      BT=int_dist2upgma_tree (s,A, A->nseq, tree_name);
      T=main_read_tree (tree_name);
      free_read_tree(BT);
    }
  
  if ( strm ( cl, "dpa"))
    {
      s=dist_array2sim_array(s, 100);
      T=code_dpa_tree (T,s);
    }

  if (free_s)free_int (s, -1);
  return T;
}
	  
NT_node similarities_file2tree (char *mat)
{
  int **s;
  Alignment *A;
  char *tree_name;
  NT_node T;

  

  tree_name =vtmpnam (NULL);
  
  s=input_similarities (mat,NULL, NULL);
  
  
  A=similarities_file2aln(mat);
  s=sim_array2dist_array(s, 100);
 
  
  int_dist2nj_tree (s, A->name, A->nseq, tree_name);
  T=main_read_tree(tree_name);
  free_int (s, -1);
  return T;
}
  
NT_node compute_cw_tree (Alignment *A)
{
  char *tmp1, *tmp2, tmp3[1000];
  char command[1000];
  
  tmp1=vtmpnam (NULL);
  tmp2=vtmpnam (NULL);

  sprintf ( tmp3, "%s.ph", tmp1);
  output_clustal_aln (tmp1, A);
  sprintf ( command, "clustalw -infile=%s -tree -newtree=%s %s ", tmp1,tmp3, TO_NULL_DEVICE);
  my_system ( command);
  sprintf ( command, "mv %s %s", tmp3, tmp2);
  my_system ( command);
  return main_read_tree(tmp2);
}

NT_node compute_fj_tree (NT_node T, Alignment *A, int limit, char *mode)
{
  static int in_fj_tree;
  if (!in_fj_tree)fprintf ( stderr, "\nComputation of an NJ tree using conserved positions\n");
  
  in_fj_tree++;
  if (T && T->leaf<=2);
  else
    {
      T=aln2fj_tree(T,A,limit, mode);
      T->right=compute_fj_tree ( T->right, A, limit, mode);
      T->left=compute_fj_tree  ( T->left, A, limit, mode);
    }
  in_fj_tree--;
  return T;
}


NT_node aln2fj_tree(NT_node T, Alignment *A, int limit_in, char *mode)
{
  NT_node NT;
  Sequence *S=NULL;
  Alignment *subA=NULL;
  int fraction_gap;
  int l, limit;
  
  if (T)
    S=tree2seq (T,NULL);
  else
    S=aln2seq (A);

  
  l=0;
  for ( fraction_gap=100; fraction_gap<=100 && l<1; fraction_gap+=10)
    for ( limit=limit_in; limit>0 && l<1; limit--)
      {
	fprintf ( stderr, "\n%d %d", limit, fraction_gap);
	free_aln (subA);
	subA=extract_sub_aln2 (A,S->nseq,S->name);
	subA=filter_aln4tree (subA,limit,fraction_gap, mode);
	l=subA->len_aln;
      }

  /*  while ( subA->len_aln<1)
    {
    subA=extract_sub_aln2 (A,S->nseq,S->name);
    subA=filter_aln4tree (subA,limit,fraction_gap,mode);
    free_aln (subA);
    subA=extract_sub_aln2 (A,S->nseq,S->name);
    subA=filter_aln4tree (subA,--limit,fraction_gap, mode);
    }
  */
  NT=aln2tree (subA);
  NT=tree2fj_tree (NT);

  NT=realloc_tree (NT,A->nseq); 
  fprintf ( stderr, "Limit:%d Gap: %d Columns: %4d Left: %4d Right %4d BL:%4.2f\n",limit,fraction_gap,  subA->len_aln, (NT->right)->leaf,(NT->left)->leaf, (NT->left)->dist+(NT->right)->dist);
  
  if ( T)
    {
      NT->dist=T->dist;
      NT->parent=T->parent;
    }
  free_tree(T);
  free_aln (subA);
  free_sequence (S, -1);
  return NT;
}  
  
Alignment * filter_aln4tree (Alignment *A, int n,int fraction_gap,char *mode)
{
  char *aln_file;
  char *ungaped_aln_file;
  char *scored_aln_file;
  char *filtered_aln_file;
  char command[1000];

  
  aln_file=vtmpnam(NULL);
  ungaped_aln_file=vtmpnam (NULL);
  scored_aln_file=vtmpnam (NULL);
  scored_aln_file=vtmpnam(NULL);
  filtered_aln_file=vtmpnam(NULL);
  
  

  output_clustal_aln (aln_file, A);
  /* 1: remove columns with too many gaps*/
  sprintf ( command, "t_coffee -other_pg seq_reformat -in %s -action +rm_gap %d -output clustalw > %s", aln_file,fraction_gap, ungaped_aln_file);
  my_system ( command);
  /* 2: evaluate the alignment*/
  
  sprintf ( command, "t_coffee -other_pg seq_reformat -in %s -action +evaluate %s -output clustalw > %s", ungaped_aln_file,(mode)?mode:"categories", scored_aln_file);
  my_system ( command);
  
  /*3 extract the high scoring columns*/
  sprintf ( command, "t_coffee -other_pg seq_reformat -in %s -struc_in %s -struc_in_f number_aln -action +use_cons +keep '[%d-8]' +rm_gap -output clustalw > %s", ungaped_aln_file, scored_aln_file,n, filtered_aln_file);
  my_system ( command);
  
  free_aln (A);
  
  A=main_read_aln ( filtered_aln_file, NULL);
  print_aln (A);
  
  return A;
}

NT_node tree2fj_tree (NT_node T)
{
  NT_node L;
  
  return T;

  L=find_longest_branch (T, NULL);
  T=reroot_tree (T, L);
  return T;
}


/*********************************************************************/
/*                                                                   */
/*                                   Tree Filters and MAnipulation   */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
int tree2star_nodes (NT_node R, int n_max)
{
  if ( !R) return 0;
  else if (!R->left && !R->right)
    {
      if (n_max>=1)R->dist=0;
      return 1;
    }
  else 
    {
      int n=0;
      
      n+=tree2star_nodes (R->right, n_max);
      n+=tree2star_nodes (R->left, n_max);
      
      if (n<n_max)R->dist=0;
      return n;
    }
}
  
NT_node aln2tree (Alignment *A)
{
  NT_node **T=NULL;
  

  T=make_nj_tree (A, NULL, 0, 0, A->seq_al, A->name, A->nseq, NULL, NULL);
  tree2nleaf (T[3][0]);
  
  return T[3][0];
}
NT_node realloc_tree ( NT_node R, int n)
{
  if ( !R)return R;
  if ( !R->leaf) 
    {
      R->right=realloc_tree (R->right,n);
      R->left=realloc_tree (R->left,n);
      R->bot=realloc_tree (R->bot,n);
    }
  R->lseq=vrealloc (R->lseq, n*sizeof (int));
  R->lseq2=vrealloc (R->lseq2, n*sizeof (int));
  return R;
}

NT_node reset_boot_tree ( NT_node R, int n)
{
  if ( !R)return R;
  if ( !R->leaf) 
    {
      
      R->right=reset_boot_tree (R->right,n);
      R->left=reset_boot_tree (R->left,n);
      R->bot=reset_boot_tree (R->bot,n);
    }
  R->bootstrap=(float)n;
  	    
  return R;
}
NT_node tree_dist2normalized_tree_dist ( NT_node R, float max)
{
  if (!R)return R;
  else
    {
      tree_dist2normalized_tree_dist ( R->right, max);
      tree_dist2normalized_tree_dist ( R->left, max);
      R->bootstrap=(int)((R->dist*100)/max);
    }
  return R;
}
NT_node reset_dist_tree ( NT_node R, float n)
{
  if ( !R)return R;
  if ( !R->leaf) 
    {
      
      R->right=reset_dist_tree (R->right,n);
      R->left=reset_dist_tree (R->left,n);
      R->bot=reset_dist_tree (R->bot,n);
    }
  if (R->parent && !(R->parent)->parent && !(R->parent)->bot)R->dist=n/2;
  else R->dist=n;

  return R;
}

      
NT_node* free_treelist (NT_node *L)
{
  int n=0;
  while (L[n])free_tree (L[n++]);
  vfree (L);
  return NULL;
}
NT_node free_tree ( NT_node R)
{
  if ( !R)return R;
  
  
    
  if ( R->leaf!=1) 
    {
      R->right=free_tree (R->right);
      R->left=free_tree (R->left);
      R->bot=free_tree (R->bot);
    }
  
  free_tree_node (R);
  return R;
}

NT_node free_tree_node ( NT_node R)
{
  if (!R)return NULL;
  
  vfree (R->seqal);
  vfree (R->idist);
  vfree (R->ldist);
  vfree (R->file);
  vfree ( R->name);
  vfree ( R->lseq); vfree ( R->lseq2);
  vfree (R);
  return NULL;
}

NT_node   rename_seq_in_tree ( NT_node R, char ***list)
{
  if ( !R || !list) return R;

  if ( R->leaf!=1) 
    {
      R->right=rename_seq_in_tree (R->right, list);
      R->left=rename_seq_in_tree (R->left, list);
      R->bot=rename_seq_in_tree (R->bot, list);
    }
  else
    {
      int n=0;
      while ( list[n][0][0])
	{
	  if ( strm (list[n][0], R->name))sprintf (R->name, "%s",list[n][1]);
	  n++;
	}
    }
  return R;
}
Sequence * tree2seq    (NT_node R, Sequence *S)
{

  if ( !R)return S;
  if ( !S)
    {
      S=declare_sequence (10, 10, tree2nseq (R));
      S->nseq=0;
    }
  
  if (R->leaf==1)
    {
      sprintf ( S->name[S->nseq++], "%s", R->name);
    }
  else
    {
      S=tree2seq (R->left, S);
      S=tree2seq (R->right, S);
    }
  return S;
}

NT_node balance_tree (NT_node T)
{
  static int **list;
  NT_node NL[3];
  
  if ( !T) return T;
  else if ( T->leaf<=2)return T;
  else    
    {
      if (!list)list=declare_int (3, 2);
      
      NL[0]=T->left;
      NL[1]=T->right;
      NL[2]=T->bot;
      
      list[0][0]=(T->left)?(T->left)->leaf:0;
      list[0][1]=0;
      list[1][0]=(T->right)?(T->right)->leaf:0;
      list[1][1]=1;
      list[2][0]=(T->bot)?(T->bot)->leaf:0;
      list[2][1]=2;
      
      sort_int (list,2,0,0,2);
      
      T->left=NL[list[2][1]];
      T->right=NL[list[1][1]];
      T->bot=NL[list[0][1]];
      
      T->left=balance_tree (T->left);
      T->right=balance_tree (T->right);
      T->bot=balance_tree (T->bot);
      return T;	
    }
}
FILE * display_tree (NT_node R, int nseq, FILE *fp)
{
  int a;
  
  if ( !R);
  else 
    {
      /*
	if ( R->nseq==1)fprintf (stderr,"\n[%s] ", R->name);
	else fprintf ( stderr, "\n[%d Node] ",R->nseq);
	for ( a=0; a< R->nseq; a++) fprintf ( stderr, "[%d]", R->lseq[a]);
      */
      fprintf (fp, "\n %10s N ", R->name);
      for ( a=0; a< nseq; a++)fprintf (fp, "%d", R->lseq2[a]);
      fprintf (fp, "\n %10s D ", R->name);
      for ( a=0; a< nseq; a++)fprintf (fp, "%d", R->idist[a]);
      
      
      if (R->leaf==1) fprintf (fp, " %s", R->name);
      fprintf (fp, " :%.4f", R->dist);
      HERE ("\nGo Left");fp=display_tree (R->left, nseq, fp);
      HERE ("\nGo Right");fp=display_tree (R->right, nseq, fp);
      HERE ("\nGo Bot");fp=display_tree (R->bot, nseq, fp);
    }
  return fp;
}
int tree2nnode_unresolved (NT_node R, int *l)
{
  if ( !R)return 0;
  else if (R->leaf && R->dist==0){return 1;}
  else
    {
      int n=0;
      n+=tree2nnode_unresolved (R->right, l);
      n+=tree2nnode_unresolved (R->left, l);
      if (R->dist==0)
	{
	  return n;
	}
      else 
	{
	  if (n)l[n]++;
	  return 0;
	}
    }
    
}

int tree2nnode ( NT_node R)
{
  int n;
  if ( !R)n=0;
  else if ( R->leaf==1){R->node=1;n=1;}
  else 
    {
      n=1;
      n+=tree2nnode (R->right);
      n+=tree2nnode (R->left);
      n+=tree2nnode (R->bot);
      R->node=n;
    }
  return n;
}
int tree2nleaf (NT_node R)
{
  if ( !R)return 0;
  else if (R->leaf==1){return 1;}
  else if (R->right==NULL && R->left==NULL && R->bot==NULL){R->leaf=1; return 1;}
  else
    {
      int n=0;
      n+=tree2nleaf (R->right);
      n+=tree2nleaf (R->left);
      n+=tree2nleaf (R->bot);
      
      R->leaf=n;
      return n;
    }
}

int tree2nseq ( NT_node R)
{
 return tree2nleaf(R);
}

int tree_file2nseq (char *fname)
{
  FILE *fp;
  char *string;
  int p, a, b, c, n;
  
  string=vcalloc (count_n_char_in_file(fname)+1, sizeof (char));
  
  fp=vfopen (fname, "r");
  n=0;
  while ( (c=fgetc(fp))!=EOF){if (c=='(' || c==')' || c==',' || c==';') string[n++]=c;}
  vfclose (fp);string[n]='\0';
  
  for (n=0, p=1; p<strlen (string); p++)
    {
      a=string[p-1];
      b=string[p];

      if ( a=='(' && b==',')n++;
      else if ( a==',' && b==')')n++;
      else if ( a==',' && b==',')n++;
    }
  if ( getenv4debug("DEBUG_TREE"))fprintf (stderr, "\n[DEBUG_TREE:tree_file2nseq]%s", string);
  vfree (string);
  return n;
}
  

void clear_tree ( NT_node R)
{
  if (!R) return;
  
  R->visited=0;
  
  if ( R->leaf==1);
  else
    {
      clear_tree ( R->right);
      clear_tree ( R->left);
      clear_tree ( R->bot);
    }
}
int display_leaf_below_node (NT_node T, FILE *fp)
{
  int n=0;
  if ( !T)return 0;
  
  if ( T->leaf==1)
    {
      fprintf (fp, " %s", T->name);
      return 1;
    }
  else
    {
      n+=display_leaf_below_node ( T->right, fp);
      n+=display_leaf_below_node ( T->left, fp);
      return n;
    }
}
int display_leaf ( NT_node T, FILE *fp)
{
  int n=0;
  if ( !T)return 0;
  else if ( T->visited)return 0;
  else T->visited=1;
  
  if ( T->leaf==1)
    {
      fprintf (fp, " %s", T->name);
      return 1;
    }
  else
    {
      n+=display_leaf ( T->right, fp);
      n+=display_leaf ( T->left, fp);
      n+=display_leaf ( T->bot, fp);
      return n;
    }
}

      

      
NT_node find_longest_branch ( NT_node T, NT_node L)
  {
    
    if ( !L || T->dist>L->dist)
      {
	
	L=T;
      }

    if ( T->leaf==1)return L;
    else 
      {
	L=find_longest_branch ( T->right, L);
	L=find_longest_branch ( T->left,  L);
	return L;
      }
  }
int node2side (NT_node N);
int test_print (NT_node T);
NT_node straighten_node (NT_node N);
NT_node EMPTY;
NT_node Previous;
NT_node reroot_tree ( NT_node TREE, NT_node Right)
{
  /*ReRoots the tree between Node R and its parent*/
  NT_node NR;
  int n1, n2;
  
  if (!EMPTY)EMPTY=vcalloc (1, sizeof (NT_node));
  if ( !Right->parent)return Right;
    
  TREE=unroot_tree (TREE);
  if (Right->parent==NULL && Right->bot)
    Right=Right->bot;
  
  n1=tree2nleaf (TREE);
  
  NR=declare_tree_node(TREE->maxnseq);
  
  NR->right=Right;
  NR->left=Right->parent;
  Right->parent=NR;
  
  Right->dist=Right->dist/2;
  
  if ((NR->left)->right==Right)(NR->left)->right=EMPTY;
  else if ( (NR->left)->left==Right) (NR->left)->left=EMPTY;
  
  Previous=NULL;
  
  
  NR->left=straighten_node (NR->left);
  
  
  
  (NR->left)->parent=NR;
  (NR->left)->dist=Right->dist;
  

     
  n2=tree2nleaf(NR);
  
  if ( n1!=n2){fprintf ( stderr, "\n%d %d", n1, n2);myexit (EXIT_FAILURE);}
  return NR;
}

NT_node straighten_node ( NT_node N)
{
  NT_node Child;
  
  
  if ( N->parent)
    {
      if (N->right==EMPTY)N->right=N->parent;
      else if ( N->left==EMPTY) N->left=N->parent;
      
      Child=N->parent;
      if (Child->right==N)
	{
	  Child->right=EMPTY;
	}
      else if (Child->left==N)
	{
	  Child->left=EMPTY;
	}
      
      Previous=N;
      Child=straighten_node (Child);
      Child->parent=N;
      Child->dist=N->dist;
      return N;
    }
  else if ( N->bot && N->bot!=Previous)
    {
      if ( N->right==EMPTY)N->right=N->bot;
      else if ( N->left==EMPTY)N->left=N->bot;
      
      N->bot=NULL;
      return N;
    }
  else
    {
      N->bot=NULL;
      return N;
    }
}
int test_print (NT_node T)
{
  if ( !T) 
    {
      fprintf ( stderr, "\nEMPTY");
    }
  else if ( !T->left && !T->right)
    {
      fprintf ( stderr, "\n%s",T->name);
    }
  else
    {
      fprintf ( stderr, "\nGoing Right");
      test_print (T->right);
      fprintf ( stderr, "\nGoing Left");
      test_print (T->left);
    }
  return 1;
}
int node2side (NT_node C)
{
  if ( !C->parent) return UNKNOWN;
  else if ( (C->parent)->left==C)return LEFT;
  else if ( (C->parent)->right==C)return RIGHT;
  else return UNKNOWN;
}
NT_node straighten_tree ( NT_node P, NT_node C, float new_dist)
{
  float dist;
  
  if ( C==NULL)return NULL;


  dist=C->dist;
  C->dist=new_dist;
  C->bot=NULL;
  
  if (C->left && C->right)
    {
      C->parent=P;
    }
  else if (!C->left)
    {
      C->left=C->parent;
      C->parent=P;
    }

  if ( C->parent==P);
  else if ( C->left==NULL && C->right==NULL)
    {
      C->parent=P;
    }
  else if ( C->right==P)
    {
      C->right=C->parent;
      C->parent=P;

      C=straighten_tree(C, C->right, dist);
    }
  else if ( C->left==P)
    {
      C->left=C->parent;
      C->parent=P;
      C=straighten_tree (C, C->left, dist);
    }
  else if ( C->parent==NULL)
    {
      C->parent=P;
    }
  
  return C;
}


NT_node unroot_tree ( NT_node T)
{
  
  if (!T || T->visited) return T;
  else T->visited=1;
  
  if (T->parent==NULL)
    {      
      
      (T->right)->dist=(T->left)->dist=(T->right)->dist+(T->left)->dist;
      (T->right)->parent=T->left;
      (T->left)->parent=T->right;
      T=T->left;
      T->leaf=0;
      vfree (T->parent);
    }
  else
    {
      T->parent=unroot_tree (T->parent);
      T->right=unroot_tree (T->right);
      T->left=unroot_tree (T->left);
    }
  T->visited=0;
  return T;
}

FILE * print_tree_list ( NT_node *T, char *format,FILE *fp)
{
  int a=0;
  while ( T[a])
    {
      fp=print_tree (T[a], format, fp);
      a++;
    }
  return fp;
}
char * tree2string (NT_node T)
{
  if (!T) return NULL;
  else
    {
      static char *f;
      FILE *fp;
      
      if (!f)f=vtmpnam (NULL);
      fp=vfopen (f, "w");
      print_tree (T, "newick", fp);
      vfclose (fp);
      return file2string (f);
    }
}
char * tree2file (NT_node T, char *name, char *mode)
{
  if (!name)name=vtmpnam (NULL);
  string2file (name, mode, tree2string(T));
  return name;
}
FILE * print_tree ( NT_node T, char *format,FILE *fp)
{
  Sequence *S;
  
  tree2nleaf(T);
  S=tree2seq(T, NULL);
  
  recode_tree (T, S);
  
  free_sequence (S, -1);
  if ( format && strm (format, "binary"))
    fp=display_tree ( T,S->nseq, fp);
  else if ( ! format || strm2 (format, "newick_tree","newick"))
    {
      /*T=balance_tree (T);*/
      fp=rec_print_tree (T, fp);
      fprintf ( fp, ";\n");
    }
  else
    {
      fprintf ( stderr, "\nERROR: %s is an unknown tree format [FATAL:%s]\n", format, PROGRAM);
      myexit (EXIT_FAILURE);
    }
  return fp;
}
int print_newick_tree ( NT_node T, char *name)
{
  FILE *fp;
  fp=vfopen (name, "w");
  fp=rec_print_tree (T,fp);
  fprintf (fp, ";\n");
  vfclose (fp);
  return 1;
}
FILE * rec_print_tree ( NT_node T, FILE *fp)
{
 


  if (!T)return fp;

  if ( T->isseq)
    {
      fprintf ( fp, " %s:%.5f",T->name, T->dist);
    }
  else
    {
      if (T->left && T->right)
	{
	  fprintf ( fp, "(");fp=rec_print_tree ( T->left, fp);
	  fprintf ( fp, ",");fp=rec_print_tree ( T->right, fp);
	  fprintf ( fp, ")");
	  if (T->parent || T->dist)
	    {
	      if ( T->bootstrap!=0)fprintf (fp, " %d", (int)T->bootstrap);
	      fprintf (fp, ":%.5f", T->dist);
	    }
	}
      else if (T->left)fp=rec_print_tree (T->left, fp);
      else if (T->right)fp=rec_print_tree(T->right, fp);
    }

  return fp;
}





/*********************************************************************/
/*                                                                   */
/*                                  Tree Functions                   */
/*                                                                   */
/*                                                                   */
/*********************************************************************/

int ** make_sub_tree_list ( NT_node **T, int nseq, int n_node)
        {


/*This function produces a list of all the sub trees*/
	 

/*  	  /A */
/*  	-* */
/*  	  \      /B */
/*  	   \    / */
/*  	    ---* */
/*                  \ */
/*                   *--C */
/*                    \ */
/*                     \D */
	
/*  	Contains 4 i_nodes */
/*  	       8   nodes (internal nodes +leaves) */
/*  	       8 sub trees: */
/*  	       ABCD */
/*                 1111 */
/*  	       0111 */
/*  	       1000 */
/*  	       0100 */
/*  	       0011 */
/*  	       0001 */
/*  	       0010 */
       
	int **sub_tree_list;
	int a, n=0;
	
	
	if (T)
	     {
		 sub_tree_list=declare_int ( (n_node), nseq);
		 make_all_sub_tree_list (T[3][0],sub_tree_list, &n);

	     }
	else
	     {
		 sub_tree_list=declare_int (nseq, nseq);
		 for ( a=0; a< nseq; a++)sub_tree_list[a][a]=1;
	     }
	
	return sub_tree_list;
	}

void make_all_sub_tree_list ( NT_node N, int **list, int *n)
        {
	make_one_sub_tree_list (N, list[n[0]++]);
	if (N->leaf!=1)
	    {
	    make_all_sub_tree_list (N->left , list, n);
	    make_all_sub_tree_list (N->right, list, n);
	    }
	return;
	}

void make_one_sub_tree_list ( NT_node T,int *list)
        {
	if (T->leaf==1)
	    {
	    
	    list[T->seq]=1;
	    }
	else
	    {
	    make_one_sub_tree_list(T->left , list);
	    make_one_sub_tree_list(T->right, list);
	    }
	return;
	}


NT_node old_main_read_tree(char *treefile)
{
  /*Reads a tree w/o needing the sequence file*/
  NT_node **T;
  T=simple_read_tree (treefile);
  return T[3][0];
}



NT_node** simple_read_tree(char *treefile)
{
  int tot_node=0;
  NT_node **T;
  T=read_tree ( treefile, &tot_node,tree_file2nseq (treefile),NULL);
  return T;
}

void free_read_tree ( NT_node **BT)
{
  int a, s;
  
 
  if (!BT) return;
  
  for (s=0,a=0; a<3; a++)
    {
      vfree (BT[a]);
    }
  free_tree (BT[3][0]);
  vfree (BT);
  return;
}
  
NT_node** read_tree(char *treefile, int *tot_node,int nseq, char **seq_names)
	{

	    /*The Tree Root is in the TREE[3][0]...*/
	    /*TREE[0][ntot]--> pointer to each node and leave*/
  	char ch;
	int  a,b;

  	FILE *fp;
  	int nseq_read = 0;
  	int nnodes = 0;/*Number of Internal Nodes*/
  	int ntotal = 0;/*Number of Internal Nodes + Number of Leaves*/
  	int flag;
  	int c_seq;
	NT_node **lu_ptr;
	NT_node seq_tree, root,p;


	tot_nseq=nseq;
	rooted_tree=distance_tree=TRUE;
	
  	fp = vfopen(treefile, "r");
	fp=skip_space(fp);
  	ch = (char)getc(fp);
  	if (ch != '(')
    		{
    	  	fprintf(stderr, "Error: Wrong format in tree file %s\n", treefile);
      		myexit (EXIT_FAILURE);
    		}
  	rewind(fp);


  	lu_ptr=(NT_node **)vcalloc(4,sizeof(NT_node*));
	lu_ptr[0] = (NT_node *)vcalloc(10*nseq,sizeof(NT_node));
	lu_ptr[1] = (NT_node *)vcalloc(10*nseq,sizeof(NT_node));
	lu_ptr[2] = (NT_node *)vcalloc(10*nseq,sizeof(NT_node));
  	lu_ptr[3] =(NT_node *) vcalloc(1,sizeof(NT_node));
  	
  	seq_tree =(NT_node) declare_tree_node(nseq);
  
  	set_info(seq_tree, NULL, 0, "  ", 0.0, 0);
	

	fp=create_tree(seq_tree,NULL,&nseq_read, &ntotal, &nnodes, lu_ptr, fp);
  	fclose (fp);
	
	
  	if (nseq != tot_nseq)
     		{
        	fprintf(stderr," Error: tree not compatible with alignment (%d sequences in alignment and %d in tree\n", nseq,nseq_read);
         	myexit (EXIT_FAILURE);
     		}
	
  	if (distance_tree == FALSE)
     		{
  		if (rooted_tree == FALSE)
          		{
       	     		fprintf(stderr,"Error: input tree is unrooted and has no distances, cannot align sequences\n");
             		myexit (EXIT_FAILURE);
          		}
          	}

        if (rooted_tree == FALSE)
     		{
		  root = reroot(seq_tree, nseq,ntotal,nnodes, lu_ptr);
		  lu_ptr[1][nnodes++]=lu_ptr[0][ntotal++]=root;
		  
     		}
  	else
     		{
        	root = seq_tree;
     		}    	
	
  	lu_ptr[3][0]=root;
  	tot_node[0]=nnodes;
	
	
	
  	for ( a=0; a< ntotal; a++)
  		{
  		(lu_ptr[0][a])->isseq=(lu_ptr[0][a])->leaf;
  		(lu_ptr[0][a])->dp=(lu_ptr[0][a])->dist;
  		}
  		
  	
  	for ( a=0; a< nseq; a++)
  		{
		if (!seq_names)
		  {
		    flag=1;
		    (lu_ptr[2][a])->order=(lu_ptr[2][a])->seq=a;
		  }
		else
		  {
		    for ( flag=0,b=0; b<nseq; b++)
		      {		       
  			if ( strncmp ( (lu_ptr[2][a])->name, seq_names[b], MAXNAMES)==0)
			  {
			    flag=1;
			    
			    (lu_ptr[2][a])->order=(lu_ptr[2][a])->seq=b;
  				/*vfree ( (lu_ptr[2][a])->name);*/
			    sprintf ((lu_ptr[2][a])->name, "%s", seq_names[b]);
			  }
		      }
		  }
		/*
		if ( flag==0  && (lu_ptr[0][a])->leaf==1)
		      {
  			fprintf ( stderr, "\n%s* not in tree",(lu_ptr[2][a])->name);
  			for ( a=0; a< ntotal; a++)
			  {
			    fprintf ( stderr, "\n%d %s",(lu_ptr[2][a])->leaf, (lu_ptr[2][a])->name);
			  } 
		      }
		*/
  		}
  	
	if (seq_names)
	  {
	    int tnseq;
	    char *s;
	    char **tree_names;
	    int fail_flag=0;
	    tnseq=tree_file2nseq(treefile);
	    tree_names=vcalloc ( tnseq, sizeof (char*));
	    for (a=0; a<tnseq; a++)
	      {
		s=(lu_ptr[2][a])->name;
		tree_names[a]=s;
		if ( name_is_in_list(s, seq_names, nseq, MAXNAMES+1)==-1)
		  {
		    fprintf (stderr, "\nERROR: Sequence %s in the tree [%s] is not in the alignment[FATAL:%s]\n", s, treefile, PROGRAM);
		    fail_flag=1;
		  }
	      }
	    for (a=0; a<nseq; a++)
	      {
		s=seq_names[a];
		if ( name_is_in_list(s, tree_names, nseq, MAXNAMES+1)==-1)
		  {
		    fprintf (stderr, "\nERROR: Sequence %s in the sequences is not in the tree [%s][FATAL:%s]\n", s, treefile, PROGRAM);
		    fail_flag=1;
		  }
	      }
	    vfree (tree_names);
	    if ( fail_flag==1)myexit (EXIT_FAILURE);
	  }
	    
  	for ( a=0; a< nseq; a++)
  		{
  		p=lu_ptr[2][a];
  		c_seq=p->seq;
  		
  		while ( p!=NULL)
  			{
  			p->lseq[p->nseq]=c_seq;
  			p->nseq++;
  			p=p->parent;
  			}
  		}
  	
	
  	return lu_ptr;   	
 	}

FILE * create_linear_tree ( char **name, int n, FILE *fp)
{

  if (!name || n==0 ||!fp) return NULL;
  
  
  if (n==2)
    fprintf ( fp, "(%s,%s);",name[0],name[1]);
  else if ( n==3)
    fprintf ( fp, "((%s,%s),%s);",name[0],name[1], name[2]);
  else
    {
      int a;
      for (a=0; a<n-2; a++)fprintf (fp, "(");
      fprintf (fp, "%s, %s),", name[0], name[1]);
      for ( a=2; a<n-2; a++)fprintf ( fp, "%s),",name[a]);
      fprintf ( fp, "%s,%s);\n", name[n-2], name[n-1]);
    }
  return fp;
}
FILE * create_tree(NT_node ptree, NT_node parent,int *nseq,int  *ntotal,int  *nnodes,NT_node **lu, FILE *fp)
	{

  	int i, type;
  	float dist=0;
	float bootstrap=0;
  	char *name;
  	int ch;
	

	name=vcalloc ( MAXNAMES+1, sizeof (char));
	sprintf ( name, "   ");
  	fp=skip_space(fp);
  	ch = (char)getc(fp);
  	
  	if (ch == '(')
    		{
    		type = NODE;
      		name[0] = '\0';
      		lu[0][ntotal[0]] = lu[1][nnodes[0]] = ptree;
      		ntotal[0]++;
      		nnodes[0]++;
      		create_tree_node(ptree, parent);
      		fp=create_tree(ptree->left, ptree, nseq,ntotal,nnodes,lu,fp);
     		ch = (char)getc(fp);
                if ( ch == ',')
       			{
          		fp=create_tree(ptree->right, ptree,nseq,ntotal,nnodes,lu,fp);
			ch = (char)getc(fp);
          		if ( ch == ',')
            			{
            			
               			ptree = insert_tree_node(ptree);
               			lu[0][ntotal[0]] = lu[1][nnodes[0]] = ptree;
               			ntotal[0]++;
               			nnodes[0]++;
               			fp=create_tree(ptree->right, ptree,nseq,ntotal,nnodes,lu,fp);
               			rooted_tree = FALSE;
				if ( getenv4debug ( "DEBUG_TREE")){fprintf ( stderr, "\n[DEBUG_TREE:create_tree] Unrooted Tree");}
            			}
       			}

      		fp=skip_space(fp);
      		ch = (char)getc(fp);
    		}
 	else
    		{
     	 	type=LEAF;
      		lu[0][ntotal[0]] = lu[2][nseq[0]] = ptree;
      		ntotal[0]++;
      		nseq[0]++;
      		name[0] = ch;
      		i=1;
      		ch = (char)getc(fp);
		if ( name[0]=='\'')
		  {
		    /*This protects names that are between single quotes*/
		    while ( ch!='\'')
		      {
			if (i < MAXNAMES) name[i++] = ch;
          		ch = (char)getc(fp);
		      }
		    if (i < MAXNAMES) name[i++] = ch;
		    while ((ch != ':') && (ch != ',') && (ch != ')'))ch = (char)getc(fp);
		  }
		else
		  {
		    while ((ch != ':') && (ch != ',') && (ch != ')'))
		      {
			if (i < MAXNAMES) name[i++] = ch;
			ch = (char)getc(fp);
		      }
		  }
		
      		name[i] = '\0';
	
      		if ( i>=(MAXNAMES+1)){fprintf (stderr, "\nName is too long");myexit (EXIT_FAILURE);}
      		if (ch != ':' && !isdigit(ch))
         		{
			  /*distance_tree = FALSE*/;
         		}
		}
	if (ch == ':')
     		{
       		fp=skip_space(fp);
       		fscanf(fp,"%f",&dist);
       		fp=skip_space(fp);
		bootstrap=0;
       		}
	/*Tree with Bootstrap information*/
	else if (isdigit (ch))
	  {
	    ungetc(ch,fp);
	    fscanf(fp,"%f",&bootstrap);
	    if ( fscanf(fp,":%f",&dist)==1);
	    else dist=0;
	    fp=skip_space(fp);
	  }
	else
	  {
	    ungetc ( ch, fp);
	    skip_space(fp);
	  }

	set_info(ptree, parent, type, name, dist, bootstrap);
  	
  	
  	vfree (name);
  	return fp;
  	}

NT_node declare_tree_node (int nseq)
	{
	NT_node p;
	
	p= (NT_node)vcalloc (1, sizeof ( Treenode)); 
	p->left = NULL;
   	p->right = NULL;
   	p->parent = NULL;
   	p->dist = 0.0;
   	p->leaf = 0;
   	p->order = 0;
	p->maxnseq=nseq;
   	p->name=(char*)vcalloc (MAXNAMES+1,sizeof (char));
   	p->name[0]='\0';
   	p->lseq=(int*)vcalloc ( nseq, sizeof (int));
   	return p;
	
	}

void set_info(NT_node p, NT_node parent, int pleaf, char *pname, float pdist, float bootstrap)
   	{
 	p->parent = parent;
   	p->leaf = pleaf;
   	p->dist = pdist;
	p->bootstrap=bootstrap;
   	p->order = 0;
   	

   	sprintf (p->name, "%s", pname);
   	
	if (pleaf ==1)
    	 	{
        	p->left = NULL;
        	p->right = NULL;
     		}

   }
NT_node insert_tree_node(NT_node pptr)
	{

   	NT_node newnode;

   	newnode = declare_tree_node( pptr->maxnseq);
   	create_tree_node(newnode, pptr->parent);

   	newnode->left = pptr;
   	pptr->parent = newnode;

   	set_info(newnode, pptr->parent, 0, "", 0.0, 0);

   	return(newnode);
	}

void create_tree_node(NT_node pptr, NT_node parent)
	{
  	pptr->parent = parent;
  	pptr->left =declare_tree_node(pptr->maxnseq) ;
  	(pptr->left)->parent=pptr;
  	
  	pptr->right =declare_tree_node(pptr->maxnseq) ;    
	(pptr->right)->parent=pptr;
	}	

FILE * skip_space(FILE *fp)
	{
  	int   c;
  
  	do
     	c = getc(fp);
  	while(isspace(c));
	if ( c==EOF)
		{
		fprintf ( stderr, "\nEOF");
		myexit (EXIT_FAILURE);
		}
  	ungetc(c, fp);
  	return fp;
	}


NT_node reroot(NT_node ptree, int nseq, int ntotal, int nnodes, NT_node **lu)
	{
	NT_node p, rootnode, rootptr;
 	float   diff, mindiff=0, mindepth = 1.0, maxdist;
	int   i;
	int first = TRUE;
	
	
	
	rootptr = ptree;
	
   	for (i=0; i<ntotal; i++)
     		{
        	p = lu[0][i];
        	if (p->parent == NULL)
           		diff = calc_root_mean(p, &maxdist, nseq, lu);
        	else
           		diff = calc_mean(p, &maxdist, nseq, lu);

        	if ((diff == 0) || ((diff > 0) && (diff < 2 * p->dist)))
         		 {
              		if ((maxdist < mindepth) || (first == TRUE))
                 		{
                    		first = FALSE;
                    		rootptr = p;
                    		mindepth = maxdist;
                    		mindiff = diff;
                 		}
           		}

     		}
	if (rootptr == ptree)
     		{
        	mindiff = rootptr->left->dist + rootptr->right->dist;
        	rootptr = rootptr->right;
     		}
	
   	rootnode = insert_root(rootptr, mindiff);
	diff = calc_root_mean(rootnode, &maxdist, nseq, lu);
	return(rootnode);
	}


float calc_root_mean(NT_node root, float *maxdist, int nseq, NT_node **lu)
	{
   	float dist , lsum = 0.0, rsum = 0.0, lmean,rmean,diff;
   	NT_node p;
   	int i;
   	int nl, nr;
   	int direction;
   	
   
   	dist = (*maxdist) = 0;
   	nl = nr = 0;
   	for (i=0; i< nseq; i++)
   	  {
          p = lu[2][i];
          dist = 0.0;
          while (p->parent != root)
         	{
               	dist += p->dist;
               	p = p->parent;
           	}
         if (p == root->left) direction = LEFT;
         else direction = RIGHT;
         dist += p->dist;

         if (direction == LEFT)
           	{
           	lsum += dist;
             	nl++;
           	}
         else
           	{
             	rsum += dist;
             	nr++;
           	}
        
        if (dist > (*maxdist)) *maxdist = dist;
     	}

      lmean = lsum / nl;
      rmean = rsum / nr;

      diff = lmean - rmean;
      return(diff);
      }

float calc_mean(NT_node nptr, float *maxdist, int nseq,NT_node **lu)
	{
   	float dist , lsum = 0.0, rsum = 0.0, lmean,rmean,diff;
   	NT_node p, *path2root;
   	float *dist2node;
   	int depth = 0, i,j , n;
   	int nl , nr;
   	int direction, found;
	
	
	path2root = (NT_node *)vcalloc(nseq,sizeof(Treenode));
	dist2node = (float *)vcalloc(nseq,sizeof(float));

   	depth = (*maxdist) = dist = 0;
   	nl = nr = 0;
   	p = nptr;
   	while (p != NULL)
     		{
         	path2root[depth] = p;
         	dist += p->dist;
         	dist2node[depth] = dist;
         	p = p->parent;
         	depth++;
     		}
 
/***************************************************************************
   *nl = *nr = 0;
   for each leaf, determine whether the leaf is left or right of the node.
   (RIGHT = descendant, LEFT = not descendant)
****************************************************************************/
   	for (i=0; i< nseq; i++)
     		{
       		p = lu[2][i];
       		if (p == nptr)
         		{
            		direction = RIGHT;
            		dist = 0.0;
         		}
       		else
         		{
         		direction = LEFT;
         		dist = 0.0;
			
        		found = FALSE;
         		n = 0;
         		while ((found == FALSE) && (p->parent != NULL))
           			{
               			for (j=0; j< depth; j++)
                 			if (p->parent == path2root[j])
                    				{ 
                      				found = TRUE;
                      				n = j;
                    				}
               			dist += p->dist;
               			p = p->parent;
           			}
         
         		if (p == nptr) direction = RIGHT;
        
			}
         	if (direction == LEFT)
           		{
             		lsum += dist;
             		lsum += dist2node[n-1];
             		nl++;
           		}
         	else
           		{
             		rsum += dist;
             		nr++;
           		}

        	if (dist > (*maxdist)) *maxdist = dist;
     		}

   vfree(dist2node);
   vfree(path2root);

   

   if ( nl==0 || nr==0)
     {
       myexit (EXIT_FAILURE);
     }
   lmean = lsum / nl;
   rmean = rsum / nr;
   
   diff = lmean - rmean;
   return(diff);
}

NT_node insert_root(NT_node p, float diff)
{
   NT_node newp, prev, q, t;
   float dist, prevdist,td;

   
   newp = declare_tree_node( p->maxnseq);
   t = p->parent;
   
   
   prevdist = t->dist;
   p->parent = newp;

   dist = p->dist;

   p->dist = diff / 2;
   if (p->dist < 0.0) p->dist = 0.0;
   if (p->dist > dist) p->dist = dist;
   
   t->dist = dist - p->dist; 

   newp->left = t;
   newp->right = p;
   newp->parent = NULL;
   newp->dist = 0.0;
   newp->leaf = NODE;
   
   if (t->left == p) t->left = t->parent;
   else t->right = t->parent;

   prev = t;
   q = t->parent;

   t->parent = newp;
   
   while (q != NULL)
     {
        if (q->left == prev)
           {
              q->left = q->parent;
              q->parent = prev;
              td = q->dist;
              q->dist = prevdist;
              prevdist = td;
              prev = q;
              q = q->left;
           }
        else
           {
              q->right = q->parent;
              q->parent = prev;
              td = q->dist;
              q->dist = prevdist;
              prevdist = td;
              prev = q;
              q = q->right;
           }
    }
   
/*
   remove the old root node
*/
   q = prev;
   if (q->left == NULL)
      {
         dist = q->dist;
         q = q->right;
         q->dist += dist;
         q->parent = prev->parent;
         if (prev->parent->left == prev)
            prev->parent->left = q;
         else
            prev->parent->right = q;
         prev->right = NULL;
      }
   else
      {
         dist = q->dist;
         q = q->left;
         q->dist += dist;
         q->parent = prev->parent;
         if (prev->parent->left == prev)
            prev->parent->left = q;
         else
            prev->parent->right = q;
         prev->left = NULL;
      }
   
   return(newp);
}




/*********************************************************************/
/*                                                                   */
/*                                  TrimTC3                          */
/*                                                                   */
/*                                                                   */
/*********************************************************************/

int *aln2seq_chain (Alignment *A, int **sim,int seq1, int seq2, int limit, int max_chain);

Alignment *seq2seq_chain (Alignment *A,Alignment*T, char *arg)
{
  int **sim=NULL;
  int *buf=NULL, *seq2keep, *list, *tname;
  int a, b, c, nl;
  int sim_limit;
  int min_sim=15;
  int max_chain=20;

  /*Estimate Similarity within the incoming sequences*/
  sim=seq2comp_mat (aln2seq(A), "blosum62mt", "sim2");
  
  /*Read and store the list of sequences to keep*/
  seq2keep=vcalloc (A->nseq, sizeof (int));
  tname=vcalloc (T->nseq, sizeof (int));
  for ( a=0; a< T->nseq; a++)
    {
      tname[a]=name_is_in_list ( T->name[a], A->name, A->nseq, 100);
      if (tname[a]>=0)seq2keep[tname[a]]=1;
    }
  
  /*Consider Every Pair of Sequences within the list of sequences to keep*/

  fprintf ( stderr, "\n");
  for ( a=0; a< T->nseq-1; a++)
    {
      if (tname[a]<0) continue;
      for ( b=a+1;b<T->nseq; b++)
	{

	  if (tname[b]<0) continue;

	  buf=NULL;sim_limit=90;
	  while (!buf && sim_limit>min_sim)
	    {
	      buf=aln2seq_chain ( A, sim,tname[a],tname[b],sim_limit, max_chain);
	      sim_limit-=5;
	    }
	  
	  if ( buf)
	    {
	      for (c=0; c< A->nseq; c++)seq2keep[c]+=buf[c];
	      vfree (buf);
	    }
	  else
	    {
	      fprintf ( stderr, "\n#Could Not Find any Intermediate sequence [MAx chain %d MinID %d\n", max_chain, min_sim);
	    }
	}
    }
  
  list=vcalloc (A->nseq, sizeof (int));
  for ( nl=0,a=0; a< A->nseq; a++)
    if ( seq2keep[a])
      list[nl++]=a;
  
  A=extract_sub_aln (A, nl, list);
  
  free_int (sim, -1);
  vfree (list);
  return A;
}
int max_explore=10000000;/*Limits the number of explorations that tends to increase when id is small*/
int n_explore;

int *aln2seq_chain (Alignment *A, int **sim, int seq1, int seq2, int limit, int max_chain)
{
  int *used;
  int **chain;
  char output1[10000];
  char output2[10000];
  int a;
  int *list;
  int n, nseq=0;

  
  output1[0]=output2[0]='\0';
  used=vcalloc (A->nseq, sizeof(int));
  used[seq1]=1;
  
  if (find_seq_chain ( A, sim,used,seq1,seq1, seq2,1,limit, max_chain, &nseq))
    {
      list=vcalloc (A->nseq, sizeof (int));
      chain=declare_int (A->nseq, 2);
      for (n=0, a=0; a< A->nseq; a++)
	{
	  if ( used[a])
	    {
	      chain[n][0]=used[a];
	      chain[n][1]=a;
	      list[used[a]-1]=a;n++;
	    }
	}
      
      sprintf ( output2, "#%s %s N: %d Lower: %d Sim: %d DELTA: %d\n", A->name[list[0]], A->name[list[n-1]],n, limit,sim[list[0]][list[n-1]],limit-sim[list[0]][list[n-1]]);strcat (output1, output2);
      
      sort_int ( chain, 2, 0, 0, n-1);
      sprintf ( output2, "#");strcat(output1, output2);
      
      for ( a=0; a< n-1; a++)
	{
	  sprintf (output2, "%s -->%d -->", A->name[chain[a][1]],sim[chain[a][1]][chain[a+1][1]]);strcat ( output1, output2);
	}
      sprintf ( output2, "%s\n", A->name[chain[n-1][1]]);strcat (output1, output2);
      
      free_int (chain, -1);
      vfree (list);
    }
  else
    {
      vfree (used);
      used=NULL;
    }
  /*  fprintf ( stdout, "%s", output1);*/
  fprintf ( stderr, "%s", output1);
  n_explore=0;
  return used;
}
static int ***pw_sim;
int find_seq_chain (Alignment *A, int **sim,int *used,int seq0,int seq1, int seq2,int chain_length, int limit, int max_chain, int *nseq)
{
  int a,b, seq, seq_sim;
  
  n_explore++;
  if ( n_explore>=max_explore)
    {
      return 0;
    }
  if (!pw_sim)
    {
      pw_sim=declare_arrayN(3, sizeof (int), A->nseq, A->nseq, 3);
      for ( a=0; a< A->nseq; a++)
	{
	  for ( b=0; b<A->nseq; b++)
	    {
	      pw_sim[a][b][0]=b;
	      pw_sim[a][b][1]=sim[a][b];
	      pw_sim[a][b][2]=sim[b][seq2];
	    }
	  sort_int_inv ( pw_sim[a],3, 1, 0, A->nseq-1);
	}
    }
  
  if ( chain_length>max_chain)return 0;
  else if ( sim[seq1][seq2]>=limit)
    {
      used[seq2]=chain_length+1;
      nseq[0]++;
      return 1;
    }
  else
    {
      int delta_seq2;
      for ( a=0; a< A->nseq; a++)
	{
	  seq=pw_sim[seq1][a][0];
	  seq_sim=pw_sim[seq1][a][1];
	  delta_seq2=pw_sim[seq1][a][2]-sim[seq1][seq2];
	  
	  

	  if ( used[seq])continue;
	  else if ( seq_sim<limit)continue;
	  else
	    {
	      used[seq]=chain_length+1;
	      nseq[0]++;
	      if ( find_seq_chain( A, sim,used,seq0,seq, seq2, chain_length+1,limit, max_chain, nseq))
		{return 1;}
	      else
		used[seq]=0;
	    }
	}
      return 0;
    }
  return 0;
}


//
//
//
//
/*********************************************************************/
/*                                                                   */
/*                                   New Tree Parsing                */
/*                                                                   */
/*********************************************************************/
NT_node new_insert_node (NT_node T);
int scan_name_and_dist ( FILE *fp, char *name, float *dist);




NT_node check_tree (NT_node T);
NT_node main_read_tree (char *treefile)
{
  FILE *fp;
  Sequence *S;
  
  NT_node T;



  
  fp=vfopen (remove_charset_from_file (treefile, " \t\n\r"), "r");
  T=new_get_node (NULL,fp);
  vfclose (fp);
    
  S=tree2seq(T, NULL);
  
  T=recode_tree(T, S);
  free_sequence (S,S->nseq);
  vfree (T->file);
  T->file=vcalloc ( strlen (treefile)+1, sizeof (char));
  sprintf ( T->file, "%s", treefile);
  return T;
}

//This function codes the tree into lseq and lseq2
//lseq:  list of the N->nseq child sequences of the node
//lsseq2:Array of size Nseq, with lseq[a]=1 if sequence a is child of node N 
static int node_index;
NT_node index_tree_node    (NT_node T)
{
  if (!T)return T;
  if (!T->parent){node_index=tree2nseq (T)+1;}
  
  index_tree_node(T->left);
  index_tree_node(T->right);

  if (!T->left && !T->right)T->index=T->lseq[0]+1;
  else T->index=node_index++;
  return T;
}
  
 
  
NT_node simple_recode_tree (NT_node T, int nseq)
{

  //recodes atree wher the leafs are already coded
  if (!T) return T;
  
  
 
  
  T->nseq=0;
  
  if ( T->isseq)
    {
      
      ;
      
    }
  else
    {
      NT_node R,L;
      int a;
       vfree (T->lseq); T->lseq=vcalloc (nseq, sizeof (int));
       vfree (T->lseq2); T->lseq2=vcalloc (nseq, sizeof (int));
       vfree (T->idist); T->idist=vcalloc (nseq, sizeof (int));
       vfree (T->ldist); T->ldist=vcalloc (nseq, sizeof (int));
       
       R=simple_recode_tree (T->left,nseq);
      
       L=simple_recode_tree (T->right,nseq);
       
       if (R)for (a=0; a<R->nseq; a++)
	 {
	   T->lseq2[R->lseq[a]]=1;
	 }
       
       if (L)for (a=0; a<L->nseq; a++)
	 {
	   T->lseq2[L->lseq[a]]=1;
	 }
       
       for (a=0; a<nseq; a++)
	 {
	   if (T->lseq2[a])T->lseq[T->nseq++]=a;
	   if (T->lseq2[a])T->idist[a]=(!R)?0:R->idist[a]+((!L)?0:L->idist[a])+1;
	   if (T->lseq2[a])T->ldist[a]=(!R)?0:R->ldist[a]+((!L)?0:L->ldist[a])+(int)(T->dist*10000);
	 }
    }
  return T;
}

NT_node recode_tree (NT_node T, Sequence *S)
{

  
  if (!T) return T;
  
  
  vfree (T->lseq);  T->lseq=vcalloc (S->nseq, sizeof (int));
  vfree (T->lseq2); T->lseq2=vcalloc (S->nseq, sizeof (int));
  vfree (T->idist); T->idist=vcalloc (S->nseq, sizeof (int));
  vfree (T->ldist); T->ldist=vcalloc (S->nseq, sizeof (int));
  T->nseq=0;
  
  if ( T->isseq)
    {
      
      int i;
      i=name_is_in_list (T->name, S->name, S->nseq, -1);
      
      if (i!=-1)
	{
	  T->lseq[T->nseq++]=i;
	  T->lseq2[i]=1;
	  T->idist[i]=1;
	  T->ldist[i]=(int)(T->dist*10000);;
	}
      else
	{
	  printf_exit ( EXIT_FAILURE, stderr, "\nERROR: Sequence %s is in the Tree but Not in the  Sequence dataset [code_lseq][FATAL:%s]", T->name, PROGRAM);
	}
      
    }
  else
    {
      NT_node R,L;
      int a;
      
      R=recode_tree (T->left, S);
      
      L=recode_tree (T->right, S);
      
      if (R)
	for (a=0; a<R->nseq; a++)
	  {
	    T->lseq2[R->lseq[a]]=1;
	  }
      
      if (L)for (a=0; a<L->nseq; a++)
	{
	  T->lseq2[L->lseq[a]]=1;
	}
      
      for (a=0; a<S->nseq; a++)
	{
	  //don't count the root
	  int d;
	  
	  if ( !(T->parent) || !(T->parent)->parent)d=0;
	  else if ( T->dist==0)d=0;
	  else d=1;

	  if (T->lseq2[a])T->lseq[T->nseq++]=a;
	  if (T->lseq2[a])T->idist[a]=(!R)?0:(R->idist[a]+((!L)?0:L->idist[a])+d);
	  if (T->lseq2[a])T->ldist[a]=(!R)?0:R->ldist[a]+((!L)?0:L->ldist[a])+(int)(T->dist*10000);
	  
	}
    }
  return T;
}
int tree2split_list (NT_node T, int ns,int **sl, int* n)
{
  if (!T) return 0;
  if (!sl) return 0;
  
  tree2split_list (T->right, ns, sl, n);
  tree2split_list (T->left , ns, sl, n);
  
  if (!T->right) return 1;
  else if (T->parent && !(T->parent)->parent)return 1;
  else if ( T->dist==0)return 1;
  else
    {
      int t=0,t2=0, c=0, a;

      for (a=0; a< ns; a++)
	{
	  t2+=(a+1)*T->lseq2[a];
	  t+=T->lseq2[a];
	}

      if (t2==0) HERE ("0");
      c=(t>(ns-t))?1:0;
      sl[n[0]][ns]=t2;//Hash value for quick comparison;
     
      for (a=0; a< ns; a++)sl[n[0]][a]=(c==0)?T->lseq2[a]:(1-T->lseq2[a]);
      n[0]++;
      
    }
  return 1;
}

NT_node display_splits (NT_node T,Sequence *S, FILE *fp)
{
  int a;
  if (!T) return T;
  
  if (!S)S=tree2seq (T,NULL);
  
  display_splits (T->right,S, fp);
  display_splits (T->left, S, fp);
  
  
 
  if (!T->right);
  else if (T->parent && !(T->parent)->parent);
  else  
    {
      int t=0;
      for (a=0; a< S->nseq; a++)
	{
	  fprintf (fp, "%d", T->lseq2[a]);
	  t+=T->lseq2[a];
	}
      
      fprintf ( fp, " %5d \n", MIN(t,((S->nseq)-t)));
    }
  return T;
}
NT_node display_leaf_nb (NT_node T, int n, FILE *fp, char * name)
{
  int a;
  if (!T) return T;
  
  
  display_leaf_nb (T->right, n, fp, name);
  display_leaf_nb (T->left, n, fp, name);
  
 
  if (!T->isseq);
  else
    {
      NT_node P;
 
      P=T->parent;
      fprintf (fp, "%s ", T->name);
      for (a=0; a< n; a++)fprintf (fp, "%d", P->lseq2[a]);
      fprintf ( fp," %s\n", name);
    }
  return T;
}
static int root4dc;
NT_node display_code (NT_node T, int n, FILE *fp)
{
  int a, debug=0, t=0;
  if (!T) return T;
  
  if (!T->parent)
    root4dc=0;

  
  
  if (!T->parent && debug) fprintf ( fp, "\nDISPLAY TREE: START");
  display_code (T->right, n, fp);
  display_code (T->left, n, fp);

  fprintf ( fp, "\n");
  if (!T->parent) return T;
  else if ( !(T->parent)->parent && root4dc==1)return T;
  else if ( !(T->parent)->parent && root4dc==0)root4dc=1;
  
  for (a=0; a< n; a++)
    t+=T->lseq2[a];
  if ( t<=n/2)
    for (a=0; a< n; a++)fprintf (fp, "%d", T->lseq2[a]);
  else
    for (a=0; a< n; a++)fprintf (fp, "%d", 1-T->lseq2[a]);
  if (T->isseq && debug)fprintf (fp, "%s", T->name);
  
  if (!T->parent && debug) fprintf (fp, "\nDISPLAY TREE: FINISHED");
  return T;
}
NT_node display_dist (NT_node T, int n, FILE *fp)
{
  int a;
  if (!T) return T;
  
  if (!T->parent)
    root4dc=0;
  
  display_dist (T->right, n, fp);
  display_dist (T->left, n, fp);
  
  fprintf ( stdout, "\n");
  for ( a=0; a< n; a++)
    fprintf ( stdout, " %2d ", T->idist[a]);
  fprintf ( stdout, "\n");
  
  return T;
}
  
NT_node check_tree (NT_node T)
{
  if (T) HERE("CHECK %s", T->name);
  if (!T)
    {
      HERE ("ERROR: Empty Group");
    }
  
  else if (T->isseq)return T;
  else 
    {
      HERE ("R");
      check_tree (T->right);
      HERE ("L");
      check_tree (T->left);
      return NULL;
    }
  return 0;}

NT_node new_reroot_tree( NT_node T)
{
  T=unroot_tree (T);
  return T;
}
NT_node new_get_node (NT_node T, FILE *fp)
{
  NT_node NN;
  int c;
  static int n;

  c=fgetc (fp);
  if (!T)T=declare_tree_node (100);
  
  
  if ( c==';')
    {
     
      if (!T->right)T=T->left;
      else if (!T->left)T=T->right;
      vfree (T->parent);T->parent=NULL;
      return T;
    }
  else if ( c==')')
    {
      --n;
      scan_name_and_dist (fp, T->name, &T->dist);
      return new_get_node (T->parent, fp);
    }
  else if ( c==',')
    {
      return new_get_node (T, fp);
    }
  else
    {
      NN=new_insert_node (T);
      
      if ( c=='(')
	{
	  ++n;
	  return new_get_node (NN, fp);
	}
      else
	{
	  ungetc (c, fp);
	  scan_name_and_dist (fp, NN->name, &NN->dist);
	  
	  NN->leaf=1;
	  NN->isseq=1;
	  return new_get_node (T, fp);
	}
    }
}
int scan_name_and_dist ( FILE *fp, char *name, float *dist)
{
  int a, c;
  char number [1000];
  
  a=0;
  c=fgetc (fp);ungetc (c, fp);
  
  
  if ( c==';')return 0;
    
  while ((c=fgetc(fp))!=':' && c!=EOF && c!=')' && c!=';' && c!=',')
    {
      name[a++]=c;
    }
  name [a]='\0';

  if ( c!=':')
    {
      ungetc (c, fp);
      dist[0]=FLT_MIN;
      return 1;
    }
  a=0;
  while (isdigit((c=fgetc(fp))) || c=='.' || c=='-')
    {
      number[a++]=c;
    }

  ungetc (c, fp);
  number[a]='\0';
  
  dist[0]=atof (number);
  return 2;
}
NT_node new_insert_node (NT_node T)
{
  NT_node NN;


  NN=new_declare_tree_node ();
  NN->parent=T;
  if (!T)
    {
      return NN;
    }
  else if (T->left==NULL)
    {
      T->left=NN;

    }
  else if ( T->right==NULL)
    {
      T->right=NN;
    }
  else
    {
      NT_node NN2;
      NN2=new_declare_tree_node ();
      NN2->left=T->left;
      NN2->right=T->right;
      NN2->parent=T;

      T->left=NN2;
      T->right=NN;
    }

  /*
    
  else
    {
      NN->right=T->right;
      (T->right)->parent=NN;
      
      NN->parent=T;
      T->right=NN;
      NN->left=new_declare_tree_node ();
      (NN->left)->parent=NN;
      return NN->left;
    }
  */
  /*
    This caused a crash when internal undefined nodes, removed 19/02/08
   else
    {
      NT_node P;
      NN->right=T;
      P=NN->parent=T->parent;
      T->parent=NN;
      
      if (P && P->right==T)P->right=NN;
      else if ( P && P->left==T)P->left=NN;
      
      NN->left=new_declare_tree_node ();
      (NN->left)->parent=NN;
      return NN->left;
    }
  */
  return NN;
}

NT_node new_declare_tree_node ()
{
	NT_node p;
	static int node_index;
	p= (NT_node)vcalloc (1, sizeof ( Treenode)); 
	p->left = NULL;
   	p->right = NULL;
   	p->parent = NULL;
   	p->dist = 0.0;
   	p->leaf = 0;
   	p->order = 0;
	p->index=++node_index;
	p->maxnseq=1000;
   	p->name=(char*)vcalloc (MAXNAMES+1,sizeof (char));
   	p->name[0]='\0';
   	return p;
	
	}
int new_display_tree (NT_node T, int n)
{
  int in;

  in=n;
  
  
  if ( T->parent)fprintf (stdout, "\nNode %d: has parents)", in);
  else fprintf (stdout, "\nNode %d: NO parents)", in);

  if ( T->right)
    {
      fprintf (stdout, "\nNode %d has Right Child", in);
      n=new_display_tree (T->right, n+1);
    }
  else fprintf ( stdout, "\nNode %d No Right\n", in);
  
  if ( T->left)
    {
      fprintf (stdout, "\nNode %d has Left Child", in);
      n=new_display_tree (T->left, n+1);
    }
  else fprintf ( stdout, "\nNode %d No Left\n", in);
  
  if ( T->bot)
    {
      fprintf (stdout, "\nNode %d has Bot Child", in);
      n=new_display_tree (T->bot, n+1);
    }
  else fprintf ( stdout, "\nNode %d No Bot\n", in);
	
  
  if (T->isseq)
    {
      fprintf (stdout, "\nNode %d is %s", in, T->name);
      return in;
    }
  else return 0;}
int display_tree_duplicates (NT_node T)
{
  static Sequence *S;
  static int *dup;
  int a, b;
  
  free_sequence (S, -1);
  vfree (dup);
  
  S=tree2seq (T, NULL);
  dup=vcalloc ( S->nseq, sizeof (int));
  
  for (a=0; a< S->nseq-1; a++)
    for ( b=a+1; b<S->nseq; b++)
      {
	if ( strm (S->name[a], S->name[b]))
	  {
	    dup[a]++;
	  }
      }
  for (a=0; a< S->nseq-1; a++)
    for ( b=a+1; b<S->nseq; b++)
      {
	if ( strm (S->name[a], S->name[b]) && dup[a])
	  {
	    fprintf ( stderr, "\nSequence %s is duplicated %d Times in the tree", S->name[a], dup[a]);
	    dup[a]=0;
	  }
      }
  return 0;
}
int tree_contains_duplicates (NT_node T)
{
  static Sequence *S;
  int a, b;
  
  free_sequence (S, -1);

  S=tree2seq (T, NULL);
  for (a=0; a< S->nseq-1; a++)
    for ( b=a+1; b<S->nseq; b++)
      {
	if ( strm (S->name[a], S->name[b]))return 1;
      }
  return 0;
}
 
float display_avg_bootstrap ( NT_node T)
{
  float tot;
  int n;
  
  tot=tree2tot_dist (T, BOOTSTRAP);
  n=tree2n_branches (T, BOOTSTRAP);
  fprintf ( stdout, "\nAVERAGE BOOTSRAP: %.3f on %d Branches\n", (n>0)?tot/n:0, n);
  return (n>0)?tot/n:0;
}
	    

int tree2n_branches(NT_node T, int mode)
{
  int n=0;

  if (!T) return 0;
  if (!T->parent);
  else if  ((T->isseq && mode !=BOOTSTRAP) || !T->isseq)
    {
      n++;
    }
  n+=tree2n_branches(T->right, mode);
  n+=tree2n_branches(T->left, mode);

  return n;
}
  
float tree2tot_dist ( NT_node T, int mode)
{
  float t=0;
  
  
  if ( !T)return 0;
  
  if ( !T->parent);
  else if  ((T->isseq && mode !=BOOTSTRAP) || !T->isseq) 
    {
      t+=T->dist;
    }
  
  t+=tree2tot_dist(T->right, mode);
  t+=tree2tot_dist(T->left, mode);
  return t;
}

//This function displays all the sequences within the tree sorted by node label
int cmp_tree_array ( const void *vp, const void *vq);
int node_sort ( char *name, NT_node T)
{
  NT_node N;
  int nseq;
  int **array, a;
  Sequence *S;
  while (T->parent)T=T->parent;
  
  nseq=tree2nseq (T);
  array=declare_int (nseq, 2);
  N=tree2node (name, T);
  
  if (N==NULL)printf_exit (EXIT_FAILURE, stderr, "ERROR: %s is not in the tree [FATAL:%s]\n", name, PROGRAM);
  array=display_tree_from_node (N,0,0, array);
  qsort ( array, nseq, sizeof (int**), cmp_tree_array);
  S=tree2seq(T, NULL);
  for (a=0; a<nseq; a++)
    fprintf ( stdout, ">%s %d %d\n", S->name[array[a][0]], array[a][1], array[a][2]);
  myexit (EXIT_SUCCESS);
}

NT_node tree2root ( NT_node R)
{
  if (R)while (R->parent)R=R->parent;
  return R;
}

NT_node tree2node (char *name, NT_node T)
{
  NT_node T1, T2;
  if ( !T) return T;
  else if (T->leaf && strm (T->name, name)) return T;
  else 
    {

      T1=tree2node ( name, T->right);
      T2=tree2node ( name, T->left);
      return (T1>T2)?T1:T2;
  }
  
}
NT_node * tree2node_list (NT_node T, NT_node *L)
{
  if (!T) return NULL;
  if (!L) L=vcalloc (T->node+1, sizeof (NT_node));
  
  tree2node_list (T->left, L);
  tree2node_list (T->right, L);
  L[T->index]=T;
  return L;
}
  
  
  
int ** display_tree_from_node (NT_node T, int up, int down, int **array)
{
  
  if (!T || T->visited)return array;
  
  T->visited=1;
  if (T->isseq)
    {
      array[T->lseq[0]][0]=T->lseq[0];
      array[T->lseq[0]][1]=up;
      array[T->lseq[0]][2]=down;
      
    }
  else
    {
      array=display_tree_from_node ( T->left ,up, down+1, array);
      array=display_tree_from_node ( T->right,up, down+1, array);
    }
  array=display_tree_from_node ( T->parent,up+1, 0, array);
  T->visited=0;
  return array;

}

int cmp_tree_array ( const void *p, const void *q)
{
  const int **vp=(const int**)p;
  const int **vq=(const int**)q;
  if (vp[0][1]>vq[0][1])return 1;
  else if ( vp[0][1]<vq[0][1]) return -1;
  else if ( vp[0][2]>vq[0][2]) return 1;
  else if ( vp[0][2]<vq[0][2]) return -1;
  else return 0;
}

NT_node * read_tree_list (Sequence *S)
{
  NT_node *T;
  int a;
  
  T=vcalloc ( S->nseq+1, sizeof (NT_node));
  
  for ( a=0; a<S->nseq; a++)
    {
      char *fname;
      if (S->seq && S->seq[a] && strlen (S->seq[a])<2)
	fname=S->name[a];
      else
	string2file ((fname=vtmpnam(NULL)), "w", S->seq[a]);
      
      T[a]=main_read_tree (fname);
      T[a]->file=vcalloc (strlen (S->name[a])+1, sizeof (char));
      sprintf (T[a]->file, "%s", S->name[a]);
    }
  return T;
}

int treelist2dmat ( Sequence *S)
{
  NT_node *T;
  int n=0, a, b;
  float v;
  Sequence *TS;


  
  n=S->nseq;
  T=read_tree_list (S);
  TS=tree2seq(T[0], NULL);
  fprintf (stdout, "\n%d", S->nseq);
  for (a=0; a<n; a++)
    {
      fprintf ( stdout,"\n%-10s ", S->name[a]); 
      for ( b=0; b<n; b++)
	{
	  v=100-simple_tree_cmp (T[a], T[b], TS, 1);
	  fprintf ( stdout, "%.4f ", v);
	}
      
    }
  
  myexit (EXIT_SUCCESS);
  return 0;
}

int treelist2leafgroup ( Sequence *S, Sequence *TS, char *taxon)
{
  NT_node *T;
  int n=0,nseq, a, c,s;
  
  int *used;
  
  char *split_file, *sorted_split_file;
  char *buf=NULL;
  FILE *fp;
  char *name, *fname, *group, *ref_group, *list;
  

  
  T=read_tree_list (S);
  if (!TS)TS=tree2seq(T[0], NULL);
  
  name=vcalloc (1000, sizeof (char));
  fname=vcalloc (1000, sizeof (char));
  group=vcalloc (TS->nseq*10, sizeof (char));
  ref_group=vcalloc (TS->nseq*10, sizeof (char));
  list=vcalloc (100*S->nseq, sizeof (char));
  split_file=vtmpnam (NULL);
  sorted_split_file =vtmpnam (NULL);
  
  n=S->nseq;
  used=vcalloc (n, sizeof (int));
 
  T=read_tree_list (S);
  if (!TS)TS=tree2seq(T[0], NULL);
  nseq=TS->nseq;
  fp=vfopen (split_file, "w");
  
  for ( a=0; a< S->nseq; a++)
    {
      
      T[a]=prune_tree  (T[a], TS);
      T[a]=recode_tree (T[a], TS);
      display_leaf_nb (T[a], TS->nseq,fp, S->name[a]);
    }
  vfclose (fp);
  
  
  for (s=0; s< TS->nseq; s++)
    {
      int i;
      

      if (taxon && !(strm (taxon, TS->name[s]) ))continue;
      else
	printf_system ( "cat %s | grep %s| sort > %s::IGNORE_FAILURE::", split_file,TS->name[s], sorted_split_file);
     
      vfopen (sorted_split_file, "r");
      ref_group[0]=group[0]='\0';
      
      while ( (c=fgetc (fp))!=EOF)
	{
	  
	  ungetc (c, fp);
	  buf=vfgets (buf, fp);
	  sscanf (buf, "%s %s %s\n", name, group, fname);

	  if ( !ref_group[0]|| !strm (group, ref_group))
	    {
	      if (ref_group[0])
		
		{fprintf (stdout, "%s %6.2f %s",name, (((float)n*100)/(float)S->nseq), ref_group);
		  for (i=0,a=0; a<nseq; a++)
		    if (ref_group[a]=='1' && a!=s)fprintf (stdout, " %s ", TS->name[a]);
		  fprintf ( stdout, "\n");
		  fprintf (stdout, "\nLIST: %s\n", list);
		}
	      list[0]='\0';
	      sprintf ( ref_group, "%s", group);
	      list=strcatf (list, " %s", fname);
	      n=1;
	    }
	  else
	    {

	      list=strcatf (list, " %s", fname);
	      n++;
	    }
	}

      fprintf (stdout, "%s %6.2f %s",name, (((float)n*100)/(float)S->nseq), group);
      for (i=0,a=0; a<nseq; a++)
	if (group[a]=='1' && a!=s)fprintf (stdout, " %s ", TS->name[a]);
      fprintf (stdout, "\nLIST %s\n", list);
      fprintf ( stdout, "\n");
      vfclose (fp);
    }
  
  myexit (0);
}
int count_tree_groups( Sequence *LIST, char *group_file)
{
  NT_node *T;
  Sequence *S;
  int a, b, c,n, w, wo, ng=0;
  int  **list, ***rlist, **blist;
  char ***l;
  int *gs;



  T=read_tree_list (LIST);
  S=tree2seq(T[0], NULL);
  for ( a=0; a< LIST->nseq; a++)
    {
      T[a]=prune_tree  (T[a], S);
      T[a]=recode_tree (T[a], S);
    }


  
  gs=vcalloc (2, sizeof (int));
  list=declare_int (LIST->nseq*S->nseq*2, S->nseq+1);
  
  blist=declare_int (2, S->nseq+1);
  for ( n=0, a=0; a< LIST->nseq; a++)
    {
      int n2=0;
      tree2split_list (T[a], S->nseq, list+n, &n2);
      n+=n2;
      
      
      for (b=0; b<n2; b++)
	{
	  for ( c=0; c<S->nseq; c++)
	    list[n+b][c]=1-list[n-n2+b][c];
	}
      n+=n2;
    }
  
  if ( group_file)
    {
      rlist=declare_arrayN(3, sizeof (int), 2,LIST->nseq*S->nseq, S->nseq+1);
      l=file2list (group_file, " ");
      
      while (l[ng])
	{
	  int i, b, g;
	  if (!strstr (l[ng][1], "group")){ng++;continue;}
	  g=(strm (l[ng][1], "group2"))?0:1;
	  
	  for (b=2; b<atoi (l[ng][0]); b++)
	    {
	      if ((i=name_is_in_list(l[ng][b], S->name, S->nseq, 100))!=-1)rlist[g][gs[g]][i]=1;
	    }
	  gs[g]++;
	  ng++;
	}
    }
  else
    {
      rlist=vcalloc ( 2, sizeof (int**));
      rlist[1]=count_int_strings (list, n, S->nseq);
      gs[1]=read_array_size_new (rlist[1]);
      
      rlist[0]=declare_int (S->nseq, S->nseq);
      gs[0]=S->nseq;
      for ( a=0; a<S->nseq; a++)rlist[0][a][a]=1;
    }
 
  
  for (wo=w=0,a=0; a<gs[0]; a++)
    {
      for (c=0; c< gs[1]; c++)
	{
	  wo=w=0;
	  for (b=0; b<S->nseq; b++)
	    {
	      blist[0][b]=blist[1][b]=rlist[1][c][b];
	      blist[0][b]=(rlist[0][a][b]==1)?1:blist[0][b]; //WITH GROUP 1 
	      blist[1][b]=(rlist[0][a][b]==1)?0:blist[1][b]; //wiTHOUT gROUP 1
	      
	    }
	  for (b=0; b<n; b++)
	    {
	      int x1, x2;
	      x1=(memcmp (blist[0], list[b], sizeof (int)*S->nseq)==0)?1:0;
	      w+=x1;
	      x2=(memcmp (blist[1], list[b], sizeof (int)*S->nseq)==0)?1:0;
	      wo+=x2;
	      
	    }
	  fprintf ( stdout, "\n%d ", MIN(wo, w));
	  fprintf ( stdout, "(");
	  for (b=0; b<S->nseq; b++)if (rlist[1][c][b])fprintf ( stdout, "%s ",S->name[b]);
	  fprintf ( stdout, ") +/- (");
	  for (b=0; b<S->nseq; b++)if (rlist[0][a][b])fprintf ( stdout, "%s ",S->name[b]);
	  
	  fprintf (stdout , ") + %d - %d Delta %d", w, wo, FABS((wo-w)));
	}
    }
  myexit (0);
}
Split * print_split ( int n, int **list, Sequence *LIST, Sequence *S, char *buf, char *file);

NT_node split2tree ( NT_node RT,Sequence *LIST, char *param)
{
  Split **S;
  Alignment *A;
  S=count_splits (RT, LIST, param);
  A=seq2aln ((S[0])->S,NULL, KEEP_GAP);
  
  return split2upgma_tree (S,A, A->nseq, "no");
}

Split** count_splits( NT_node RT,Sequence *LIST, char *param)
{
  NT_node *T, OrderT;
  Sequence *S=NULL;
  int a, b, c, d, n1, n2;
  int  **list1, **list2;
  Split **SL;
  int nb, tlist;
  char *main_buf;
  char *in=NULL,*in2=NULL, *out=NULL, order[100], filter[100];
  FILE *fp, *fp2;
  char *def_param;
  char *cache=NULL;
  //+count_splits _NB_x_FILTER_<file>
  //_<file is a fasta file containing the list of species to keep>
  if (!def_param)def_param=vcalloc ( 10, sizeof (char));
  


  if (!param)param=def_param;
  
 
  strget_param (param, "_NB_", "0", "%d", &nb);
  strget_param (param, "_TLIST_", "0", "%d", &tlist);
  strget_param (param, "_ORDER_", "NO", "%s", order);
  strget_param (param, "_FILTER_", "NO", "%s", filter);
  
  fprintf ( stderr, "\nREAD TREE LIST [%d Trees...", LIST->nseq);
  T=read_tree_list (LIST);
  fprintf ( stderr, "..]");
  
  if ( !(strm (order, "NO")))
    {
      if (is_newick (order))
	{
	  OrderT=main_read_tree (order);
	}
      else 
	{
	  S=main_read_seq (order);
	}
    }
  else
    {
      OrderT=(RT)?RT:T[0];
    }
  fprintf ( stderr, "\nTrees Ordered according to: %s", (strm (order, "NO"))?"First Tree":order);
  
  
  if (!S)S=tree2seq(OrderT, NULL);
  
  for (a=0; a<S->nseq; a++)
    {
      fprintf ( stdout, "\n#ORDER %15s : %3d", S->name[a], a+1);
    }
  if ( !strm (filter, "NO"))
    {
      Sequence *F;
      int i;
      
      F=main_read_seq (filter);
      cache=vcalloc (S->nseq, sizeof (int));
      for ( a=0; a<F->nseq; a++)
	{
	  if ( (i=name_is_in_list (F->name[a], S->name, S->nseq, 100))!=-1)
	    cache[i]=1;
	}
      free_sequence (F, -1);
    }
  
  main_buf=vcalloc ( S->nseq*(STRING+1), sizeof(int));

  list1=declare_int (S->nseq*3, S->nseq+1);
  list2=declare_int (S->nseq*3, S->nseq+1);
  
  for ( a=0; a< LIST->nseq; a++)
    {
      T[a]=prune_tree  (T[a], S);
      T[a]=recode_tree (T[a], S);
    }
  

  
  if (!RT)
    {
      char *buf;
      int i,nl;
      
      in=vtmpnam (NULL);in2=vtmpnam(NULL); out=vtmpnam (NULL);
      
      fp=vfopen (in, "w");
      fp2=vfopen (in2, "w");
      for ( a=0; a< LIST->nseq; a++)
	{
	  n2=0;
	  tree2split_list (T[a], S->nseq, list2, &n2);
	  for ( b=0; b<n2; b++)
	    {
	      for (c=0; c< S->nseq; c++)
		{fprintf (fp, "%d", list2[b][c]);}
	      fprintf (fp, "\n");
	      for (c=0; c< S->nseq; c++)
		{fprintf (fp, "%d", 1-list2[b][c]);}
	      fprintf (fp, "\n");
	      
	      for (c=0; c< S->nseq; c++)
		{fprintf (fp2, "%d", list2[b][c]);}
	      fprintf (fp2, " ");
	      for (c=0; c< S->nseq; c++)
		{fprintf (fp2, "%d", 1-list2[b][c]);}
	      fprintf (fp2, " %s\n",LIST->name[a]);
	    }
	}
      vfclose (fp2);
      vfclose (fp);
     
      count_strings_in_file (in, out);
      nl=count_n_line_in_file(out);
      list1=declare_int (nl+1, S->nseq+2);
      
      fp=vfopen (out, "r");
      n1=0;
      buf=vcalloc (measure_longest_line_in_file (out)+1, sizeof (char));
      while ( fscanf (fp, "%s %d",buf, &i)==2)
	{
	  for (a=0; a<S->nseq; a++)list1[n1][a]=buf[a]-'0';
	  list1[n1++][S->nseq+1]=i;
	}
      vfclose (fp);
      vfree (buf);
    }
  else
    {
     
      RT=prune_tree  (RT, S);
      RT=recode_tree (RT, S);
      n1=0;
      tree2split_list (RT, S->nseq, list1,&n1);
      for ( a=0; a< LIST->nseq; a++)
	{
	  n2=0;
	  tree2split_list (T[a], S->nseq, list2, &n2);
	  for (b=0; b<n1; b++)
	    {
	      for ( c=0; c<n2; c++)
		{
		  int di=0;
		  for (d=0; d<S->nseq; d++)
		    {
		      if (list1[b][d]!=list2[c][d])di++;
		    }
		  list1[b][S->nseq+1]+=(di==0 || di== S->nseq)?1:0;
		}
	      
	    }
	}
    }
  SL=vcalloc ( n1+1, sizeof (Split*));
  
  for (a=0; a<n1; a++)
    {
      int s1, s2;
      int cont=1;
      if (nb)fprintf ( stdout, "\nSPLIT: %d and its Neighborhood +/^- %d\n", a+1, nb);
      
      if (cache)
	for (b=0; b<S->nseq; b++)if (cache[b]!=list1[a][b])cont=0;
      if (!cont) continue;
      
      SL[a]=print_split (a, list1, LIST, S, main_buf, (tlist==1)?in2:NULL);
      for (b=0; b<n1; b++)
	{
	  if ( a==b)continue;
	  else
	    {

	      for (d=s1=s2=0,c=0; c<S->nseq; c++)
		{
		  s1+=list1[b][c];
		  s2+=list1[a][c];
		  d+=(list1[a][c]!=list1[b][c])?1:0;
		}
	      
	    }
	  if (d<=nb &&((s1==s2)|| ((S->nseq-s1)==s2)))print_split (b, list1, LIST, S, main_buf, (tlist==1)?in2:NULL);
	}
    }
  a=0;
  vfree (cache);
  return SL;
}
Split * declare_split (int nseq, int ntrees);
Split* print_split ( int a, int **list1, Sequence *LIST, Sequence *S, char *buf, char *split_file)
  {
    int f1,t,b;
    Split *SP=NULL;
    
    

    SP=declare_split (S->nseq, LIST->nseq);
    
    fprintf ( stdout, "\n>");
    for (t=0,b=0; b<S->nseq; b++){fprintf ( stdout, "%d", list1[a][b]);t+=list1[a][b];SP->split[b]='0'+list1[a][b];}
    fprintf ( stdout, " NumberSplit %5d SplitSize %5d Score %5.2f %s ", list1[a][S->nseq+1],t, (float)(list1[a][S->nseq+1]*100)/LIST->nseq, (buf)?buf:"");
    SP->n= list1[a][S->nseq+1];
    SP->score=(float)(list1[a][S->nseq+1]*100)/LIST->nseq;
    SP->S=S;
    
    for (f1=1,b=0; b< S->nseq; b++)
      {
	
	if (list1[a][b])
	  {
	    if (f1==1)fprintf ( stdout, "(");
	    else fprintf (stdout, ",");
	    f1=0;
	    fprintf ( stdout, "%s", S->name [b]);
	  }
      }
    fprintf ( stdout, ")");
    if (split_file)
      {
	char *buf=NULL;
	FILE *fp;

	char c;
	fp=vfopen (split_file, "r");
	while ( (c=fgetc(fp))!=EOF)
	  {
	    
	    c=ungetc (c, fp);
	    buf=vfgets (buf, fp);
	    if ( strstr (buf, SP->split))
	      {
		char **list;
		list=string2list (buf);
		fprintf ( stdout, "\n\t%s %s", SP->split, list[3]);
		free_char (list, -1);
	      }
	  }
	vfclose (fp);
      }
    
    return SP;
  }
Split * declare_split (int nseq, int ntrees)
{
  Split *S;
  S=vcalloc (1, sizeof (Split));
  S->split=vcalloc ( nseq+1, sizeof (char));
  return S;
}
int treelist2splits( Sequence *S, Sequence *TS)
{
  NT_node *T;
  int n=0,nseq, a, c;

  int *used;

  char *split_file, *sorted_split_file;
  char *buf=NULL, *ref_buf=NULL;
  FILE *fp;
  
  split_file=vtmpnam (NULL);
  sorted_split_file =vtmpnam (NULL);
  
  n=S->nseq;
  used=vcalloc (n, sizeof (int));
 
  T=read_tree_list (S);
  if (!TS)TS=tree2seq(T[0], NULL);
  nseq=TS->nseq;
  fp=vfopen (split_file, "w");
  
  
  for ( a=0; a< S->nseq; a++)
    {
      
      T[a]=prune_tree  (T[a], TS);
      T[a]=recode_tree (T[a], TS);
      display_splits (T[a], TS,fp);
    }
  
  vfclose (fp);
  printf_system ("cp %s split_file::IGORE_FAILURE::", split_file);
  
  printf_system ( "cat %s | grep 1| sort > %s::IGNORE_FAILURE::", split_file, sorted_split_file);
  
  fp=vfopen (sorted_split_file, "r");
  fprintf (stdout, "LEGEND: <#occurences> <coded split> <min group size> <(group1,)> <(group2,>\n");
  
  for ( a=0; a<TS->nseq; a++)fprintf ( stdout, "SEQ_INDEX %d %s\n", a+1, TS->name[a]);
  while ( (c=fgetc (fp))!=EOF)
    {
      
      ungetc (c, fp);
      buf=vfgets (buf, fp);
      buf [strlen(buf)-1]='\0';
      
      if ( ref_buf==NULL)
	{
	  ref_buf=vcalloc (strlen (buf)+1, sizeof (char));
	  sprintf ( ref_buf, "%s", buf);
	  n=1;
	}
      else if ( !strm (buf, ref_buf))
	{
	  int i;
	  fprintf ( stdout, "SPLIT_COUNT %3d %s (", n, ref_buf);
	  for (i=0,a=0; a<nseq; a++)
	    if (ref_buf[a]=='1')
	      {
		if (i==1)fprintf(stdout, ",");
		fprintf (stdout, "%s", TS->name[a]);
		i=1;
	      }
	  fprintf ( stdout, "),(");
	  for (i=0,a=0; a<nseq; a++)
	    if (ref_buf[a]=='0')
	      {
		if (i==1) fprintf ( stdout, ",");
		fprintf (stdout, "%s", TS->name[a]);
		i=1;
	      }
		  
	  fprintf (stdout, ")\n");
	  sprintf ( ref_buf, "%s", buf);
	  n=1;
	}
      else
	{
	  n++;
	}
    }
  vfclose (fp);
      
  
  myexit (0);
}

int treelist2splits_old ( Sequence *S, Sequence *TS)
{
  NT_node *T;
  int n=0,nseq, a,c;
  
  int *used;
  
  char *split_file, *sorted_split_file;
  char *buf=NULL, *ref_buf=NULL;
  FILE *fp;
  
  split_file=vtmpnam (NULL);
  sorted_split_file =vtmpnam (NULL);
  
  n=S->nseq;
  used=vcalloc (n, sizeof (int));
 
  T=read_tree_list (S);
  if (!TS)TS=tree2seq(T[0], NULL);
  nseq=TS->nseq;
  fp=vfopen (split_file, "w");
  
  for ( a=0; a< S->nseq; a++)
    {
      
      T[a]=prune_tree  (T[a], TS);
      T[a]=recode_tree (T[a], TS);
      display_leaf_nb (T[a], TS->nseq,fp, S->name[a]);
    }
  vfclose (fp);
  printf_system ("cp %s split_file::IGNORE_FAILURE::", split_file);myexit (0);

  printf_system ( "cat %s | grep 1| sort > %s::IGNORE_FAILURE::", split_file, sorted_split_file);
  
  vfopen (sorted_split_file, "r");

  while ( (c=fgetc (fp))!=EOF)
    {
      
      ungetc (c, fp);
      buf=vfgets (buf, fp);
      buf [strlen(buf)-1]='\0';
      
      if ( ref_buf==NULL)
	{
	  ref_buf=vcalloc (strlen (buf)+1, sizeof (char));
	  sprintf ( ref_buf, "%s", buf);
	  n=1;
	}
      else if ( !strm (buf, ref_buf))
	{
	  int i;
	  fprintf ( stdout, "%3d %s(", n, ref_buf);
	  for (i=0,a=0; a<nseq; a++)
	    if (ref_buf[a]=='1')
	      {
		if (i==1)fprintf(stdout, ",");
		fprintf (stdout, "%s", TS->name[a]);
		i=1;
	      }
	  fprintf ( stdout, "),(");
	  for (i=0,a=0; a<nseq; a++)
	    if (ref_buf[a]=='0')
	      {
		if (i==1) fprintf ( stdout, ",");
		fprintf (stdout, "%s", TS->name[a]);
		i=1;
	      }
		  
	  fprintf (stdout, ")\n");
	  sprintf ( ref_buf, "%s", buf);
	  n=1;
	}
      else
	{
	  n++;
	}
    }
  vfclose (fp);
      
  
  myexit (0);
}

NT_node *treelist2prune_treelist (Sequence *S, Sequence *TS, FILE *out)
{
  NT_node *T;
  int a, b, c;

  T=read_tree_list (S);
  T=vrealloc (T, (S->nseq+1)*sizeof (NT_node));
  for (b=0,a=0; a<S->nseq; a++)
    {
      T[a]=prune_tree  (T[a], TS);
      if (tree2nleaf(T[a])<TS->nseq)
	{
	  ;
	}
      else 
	{
	  char *s;
	  T[b]=T[a];
	  T[b]=recode_tree (T[b], TS);
	  sprintf ( S->name[b], "%s", S->name[a]);
	  s=tree2string (T[a]);
	  S->seq[b]=vrealloc (S->seq[b], (strlen (s)+1)*sizeof (char));
	  sprintf (S->seq[b], "%s",s);
	  sprintf (S->seq_comment[b], " NSPECIES: %d", TS->nseq); 
	  vfree (s);

	  b++;
	}
      
    }

  S->nseq=b;
  T[S->nseq]=NULL;

  if (out)
    {
      for (a=0; a<S->nseq; a++)
	{
	  print_tree (T[a], "newick", out);
	}
    }
  return T;
}
int** treelist2lti2 ( Sequence *S, Sequence *TS, int ngb, FILE *out);
int treelist2frame (Sequence *S, Sequence *TS)
{
  int n, a, b, c,d, **r, **order;
  Sequence *temp;
  
  temp=duplicate_sequence (S);
  order= treelist2lti (temp, TS,0,stdout);
  
  TS=reorder_seq_2 (TS, order, 0, TS->nseq);
  n=TS->nseq;
  
  for (a=3; a<n; a++)
    {
      NT_node tree;

      TS->nseq=a+1;
      temp=duplicate_sequence (S);
      r=treelist2groups (temp,TS, NULL, NULL);
      fprintf ( stdout, "\n>Tree_%d [%d %%]\n ", a+1,r[0][1]);
      tree=main_read_tree (temp->name[r[0][0]]);
      tree=prune_tree (tree, TS);
      print_tree (tree, "newick",stdout);
      
      free_int (r, -1);
      free_sequence (temp,-1);
    }
  myexit (EXIT_SUCCESS);
}
int** treelist2lti2 ( Sequence *S, Sequence *TS, int ngb, FILE *out)
{
  NT_node *T;
  int a,b, c, d, ****dist, i;
  int **score, **order;

  score=declare_int (TS->nseq, 3);
  order=declare_int (TS->nseq, 2);
  vsrand (0);
  
  for (a=0; a<50; a++)
    {
      Sequence *seq, *trees;
      int **r;
      trees=duplicate_sequence (S);
      seq=duplicate_sequence (TS);
      for (b=0; b<TS->nseq; b++){order[b][0]=b;order[b][1]=rand()%10000;}
      sort_int (order, 2, 1, 0, TS->nseq-1);
      seq=reorder_seq_2(seq, order, 0,5);
      r=treelist2groups (trees,seq, NULL, NULL);
      
      for (b=0; b<5; b++)
	{
	  score[order[b][0]][1]+=r[0][1];
	  score[order[b][0]][2]++;
	}
      HERE ("Score=%d", r[0][1]);
      free_int (r, -1);
      free_sequence (seq, -1);
      free_sequence (trees, -1);
      
    }
  
  for ( a=0; a< TS->nseq; a++)
    {
      score[a][0]=a;
      HERE ("%s => %d [%d]",TS->name[a], score[a][1]/score[a][2], score[a][2]);
      score[a][1]/=(score[a][2])?score[a][2]:1;
    }
  sort_int_inv (score, 3, 1, 0, TS->nseq-1);
  
  return score;
}


int** treelist2lti ( Sequence *S, Sequence *TS, int ngb, FILE *out)
{
  NT_node *T;
  int a,b, c, d, ****dist, i;
  float score0=0, score1=0;
  int **result;
  

  i=S->nseq;
  T=treelist2prune_treelist (S, TS,NULL);
  
  if (!ngb)ngb=TS->nseq*2;
  dist=vcalloc ( S->nseq, sizeof (int****));
  result=declare_int (TS->nseq, 2);
  for (a=0; a<TS->nseq; a++)
    {
      float score_seq=0;
      float n_seq=0;
      for (b=0; b<TS->nseq;b++)
	{
	  float score_pair=0;
	  float n_pair=0;
	  for (c=0; c<S->nseq; c++)
	    {
	      if (!dist[c])dist[c]=tree2dist(T[c], TS, NULL);
	      for (d=0; d<S->nseq; d++)
		{
		  float score, d1, d2;

		  if (!dist[d])dist[d]=tree2dist(T[d], TS, NULL);
		  d1=dist[c][0][a][b];
		  d2=dist[d][0][a][b];
		  score=FABS((d1-d2));
		  if (d1>ngb || d2>ngb);
		  else
		    {
		      score_seq+=score;
		      score_pair+=score;
		      n_seq++;
		      n_pair++;
		    }
		  // if (d1 && d2) HERE ("%d %d", (int)d1, (int)d2);
		}
	    }
	  score_pair=(score_pair*100)/(float)n_pair;
	  if (out)fprintf ( stdout, "\n>%-20s %-20s LTI: %7.3f [Kept %d Trees Out of %d] ", TS->name[a],TS->name[b], score_pair, S->nseq,i);
	}

      score_seq=(score_seq*100)/n_seq;
      result[a][0]=a;
      result[a][1]=(int)(100*score_seq);
      if (out)fprintf ( stdout, "\n>%-20s %-20s LTI: %7.3f [Kept %d Trees Out of %d] ", TS->name[a],"*", score_seq, S->nseq, i);
    }
  sort_int (result,2,1,0, TS->nseq-1);
  return result;
}


int ***tree2dist (NT_node T, Sequence *S, int ***d)
{
  int *l0, *r0,*l1, *r1, a, b;
  
  
  if (!T) return d;
  if (!S)S=tree2seq(T, NULL);
  if (!d)
    {
      d=declare_arrayN (3, sizeof (float),2, S->nseq, S->nseq);
      T=prune_tree(T, S);
      T=recode_tree (T, S);
    }
  
  if (!T->left)return d;
  if (!T->right) return d;
  
  l0=(T->left)->idist;
  r0=(T->right)->idist;

  l1=(T->left)->ldist;
  r1=(T->right)->ldist;
  
  
  
  for (a=0; a< S->nseq; a++)
    for (b=0; b<S->nseq; b++)
      {
	if (l0[a]>0 && r0[b]>0)d[0][a][b]=d[0][b][a]=l0[a]+r0[b];
	if (l0[a]>0 && r0[b]>0)d[1][a][b]=d[1][b][a]=l1[a]+r1[b];
      }
  
  d=tree2dist (T->left, S, d);
  d=tree2dist (T->right, S, d);
  
 
  return d;
}



int **tree2dist_split ( NT_node T, Sequence *S, int **dist)
{
  
  FILE *fp;
  int a, b, c, n=0;
  char *buf=NULL, **list=NULL, *split_file;


  if (!S)S=tree2seq(T, NULL);
  
  T=prune_tree  (T, S);
  T=recode_tree (T, S);
  
  split_file=vtmpnam (NULL);
  fp=vfopen (split_file, "w");
  display_code (T, S->nseq,fp);
  vfclose (fp);
 
  list=declare_char (2*S->nseq, S->nseq+1);
  fp=vfopen (split_file, "r");

  while ((buf=vfgets (buf,fp))!=NULL)
    {
      if (buf[0]=='1' || buf[0]=='0')sprintf (list[n++], "%s", buf);
    }
  vfclose (fp);
  dist=declare_int ( S->nseq, S->nseq);
  for (a=0; a< S->nseq; a++)
    for ( b=0; b<S->nseq; b++)
      for (c=0; c<n; c++)
	if (list[c][a]!=list[c][b])dist[a][b]++;
  
 
  return dist;
}

int** treelist2groups (Sequence *S, Sequence *TS, char *star_node, FILE *out)
   {
   NT_node *T;
   int a, b, tot, n,i;
   int v;
   int *used;
   int ntop;
   int nsn;
   int cov=100;
   int **results;

   
   i=S->nseq;
   T=treelist2prune_treelist (S, TS,NULL);
   nsn=(star_node)?atoi(star_node):0;
   
   results=declare_int (S->nseq+1, 2);
   
   if (nsn)
     {
       for (a=0; a< S->nseq; a++)tree2star_nodes(T[a],nsn);
     }
   
   used=vcalloc (S->nseq, sizeof (int));
   for (ntop=0,a=0; a<S->nseq; a++)
     {
       
       if (used[a]==0)
	 {
	   ntop++;
	   if (out)fprintf ( out, "\nTree %s:",S->name[a]);
	   used[a]=1;
	 }
       else continue;
       tot=1;
       for ( b=0; b<S->nseq; b++)
	 {
	   v=0;
	   
	   v=(int)simple_tree_cmp (T[a], T[b], TS, 1);
	   if ( v==100)
	     {
	       used[b]=1;
	       used[a]++;
	       if (out)fprintf (stdout," %s ", S->name[b]);
	       tot++;
	     }
	 }
       
       if (out)fprintf ( stdout, "__ N=%d\n", tot-1);
     }


   for (n=0,a=0; a<S->nseq; a++)
     {
       if ( used[a]>1)
	 {
	   if (out)fprintf ( out, "\n>%-15s %4d %6.2f TOPOLOGY_LIST\n", S->name[a], used[a]-1, (float)(((float)used[a]-1)*100/(float)S->nseq));
	   if (out)print_tree (T[a], "newick_tree", out);
	   results[n][0]=a;
	   results[n][1]=((used[a]-1)*100)/i;
	   n++;
	 }
     }

   for (a=0; a<S->nseq; a++) free_tree(T[a]);
   vfree (T);

   if (out)fprintf ( stdout, "\nTotal Number of different topologies: %d\n", ntop);
   results[n][0]=-1;
   sort_int_inv (results,2,1,0, n-1);
   for (a=0; a<S->nseq; a++) free_tree(T[a]);
   vfree (T);
   return results;
   }
float simple_tree_cmp (NT_node T1, NT_node T2,Sequence *S, int mode)
{
  Tree_sim *TS1, *TS2;
  float t, w, l, n;
  
  TS1=vcalloc (1, sizeof (Tree_sim));
  TS2=vcalloc (1, sizeof (Tree_sim));
  
  
  T1=recode_tree(T1, S);
  T2=recode_tree(T2, S);
  
  n=new_compare_trees ( T1, T2, S->nseq, TS1);
  new_compare_trees ( T2, T1, S->nseq, TS2);
  


  t=(TS1->uw+TS2->uw)*100/(TS1->max_uw+TS2->max_uw);
  w=(TS1->w+TS2->w)*100/(TS1->max_w+TS2->max_w);
  l=(TS1->d+TS2->d)*100/(TS1->max_d+TS2->max_d);
  
  vfree (TS1); vfree (TS2);
  if ( mode ==1)return t;
  else if (mode ==2) return w;
  else  return l;
}
int treelist2n (NT_node *L)
{
  int n=0;
  while (L[n])n++;
  return n;
}
int **treelist2avg_treecmp (NT_node *L, char *file)
{
  int a, b, n;
  int **score;
  
  if (file) L=read_tree_list (main_read_seq(file));
  n=treelist2n (L);
  
  score=declare_int (n, 2);
  for (a=0; a<n; a++)score[a][0]=a;
  
  for (a=0; a<n-1; a++)
    {
      output_completion (stderr,a,n,1, "Tree Cmp");
      for (b=a+1; b<n; b++)
	{
	  Tree_sim *ts;
	  ts=tree_cmp (L[a],L[b]);
	  score[a][1]+=ts->uw;
	  score[b][1]+=ts->uw;
	  vfree (ts);
	}
    }
  sort_int_inv (score, 2, 1, 0, n-1);
  if (file)free_treelist(L);
  return score;
}

int treelist_file2consense (char *tree_file, char *outtree, char *outfile)
{
  static char *command;
  static char *tmp_outtree;
  static char *tmp_outfile;
  FILE *fp;
  int flag1=0; 
  int flag2=0;
  
  if (!command)
    {
      command=vtmpnam (NULL);
      tmp_outtree=vtmpnam (NULL);
      tmp_outfile=vtmpnam (NULL);
    }
  if (!check_program_is_installed ("consense",NULL,NULL,"www.phylip.com",NO_REPORT))return 0;
  
  fp=vfopen (command, "w");fprintf ( fp, "%s\nY\n", tree_file);fclose (fp);
  if ( check_file_exists ("outtree")){flag1=1;printf_system ("mv outtree %s::IGNORE_FAILURE::", tmp_outtree);}
  if ( check_file_exists ("outfile")){flag2=1;printf_system ("mv outfile %s::IGNORE_FAILURE::", tmp_outfile);}
  printf_system ("consense <%s > /dev/null 2>/dev/null::IGNORE_FAILURE::", command);
  
  if ( outtree)printf_system ("mv outtree %s::IGNORE_FAILURE::", outtree);
  remove ("outtree");
  if ( outfile)printf_system ("mv outfile %s::IGNORE_FAILURE::", outfile);
  remove ("outfile");
  if (flag1)printf_system ("mv %s outtree::IGNORE_FAILURE::", tmp_outtree);
  if (flag2)printf_system ("mv %s outfile::IGNORE_FAILURE::", tmp_outfile);
  return 1;
}
  
  
NT_node treelist2filtered_bootstrap ( NT_node *L,char *file, int **score, float t)
{
  NT_node BT, *L2;
  int n,a;
 
  if (t==1 || t==0 || !score)return treelist2bootstrap (L, file);
 
  if (file)L=read_tree_list (main_read_seq(file));

  n=treelist2n(L)*t;

  if (n==0) return NULL;
  
  L2=vcalloc ( n+1, sizeof (NT_node));
  for (a=0; a<n; a++)
    L2[a]=L[score[a][0]];
    
  BT=treelist2bootstrap (L2, NULL);
  
  vfree (L2);
  if (file)free_treelist(L);
  return BT;
}

NT_node treelist2bootstrap ( NT_node *L, char *file)
{
  char *outfile;
  NT_node T;
  FILE *fp;

  if (!file)
    {
      file=vtmpnam (NULL);
      vfclose (print_tree_list (L,"newick", vfopen (file, "w"))); 
    }
 
  outfile=vtmpnam (NULL);
 
  printf_system( "msa2bootstrap.pl -i %s -o %s -input tree >/dev/null 2>/dev/null", file, outfile);
 
  T=main_read_tree (outfile);
  T=tree_dist2normalized_tree_dist (T,treelist2n(L));
  

  return T;
}


  
Sequence * treelist2seq (Sequence *S)
{
  int a, b, c, n, i;
  char **name;
  NT_node *T;
  Sequence *TS;
  char *fname;
  FILE *fp;
  
  name=vcalloc (1, sizeof (char*));
  fp=vfopen ((fname=vtmpnam (NULL)), "w");
 
  T=read_tree_list (S);
  for (n=0,a=0; a< S->nseq; a++)
    {
      TS=tree2seq(T[a], NULL);
      for (b=0; b<TS->nseq; b++)
	{
	  if ( (i=name_is_in_list (TS->name[b], name, n, 100))==-1)
	    {
	      name[n]=vcalloc (100, sizeof (int));
	      sprintf ( name[n], "%s", TS->name[b]);
	      n++;
	      name=vrealloc (name, (n+1)*sizeof (char*));
	      fprintf ( fp, ">%s\n", TS->name[b]);
	    }
	}
      free_sequence(TS, TS->nseq);
      free_tree (T[a]);
    }
  
  vfclose (fp);
  vfree (T);
  return get_fasta_sequence (fname, NULL);
}


Sequence * treelist2sub_seq ( Sequence *S, int f)
{
  NT_node *T;
  int a,b,c, s, i, n, maxnseq, tot;
  int **count, **grid;
  char *fname;
  Sequence *FS, *TS;
  FILE *fp;
  if (!f)return treelist2seq(S);
  

  //keep as many taxons as possible so that f% of the trees are kept
  //1: count the frequency of each taxon
  
  FS=treelist2seq (S);
  maxnseq=FS->nseq;
  
  count=declare_int (maxnseq, 3);
  grid=declare_int (S->nseq,maxnseq+1);
  T=read_tree_list (S);
  

  
  for (a=0; a<FS->nseq; a++){count[a][0]=a;count[a][2]=1;}
  for (n=0,a=0; a< S->nseq; a++)
    {
      TS=tree2seq(T[a], NULL);
      for (b=0; b<TS->nseq; b++)
	{
	  i=name_is_in_list (TS->name[b], FS->name, FS->nseq, 100);
	  if ( i==-1){myexit (EXIT_FAILURE);}
	  count[i][1]++;
	  grid[a][i]=1;
	}
      free_sequence(TS, TS->nseq);
      free_tree (T[a]);
    }
  vfree (T);
  sort_int ( count,3,1, 0, maxnseq-1);
  
  for (a=0; a<maxnseq; a++)
    {
      count[a][2]=0;
      for ( b=0; b< S->nseq; b++)grid[b][maxnseq]=1;//prepare to keep everything
      for ( tot=S->nseq, b=0; b< S->nseq; b++)
	{
	  for (c=0; c<maxnseq; c++)
	    {
	      s=count[c][0];
	      if (count[c][2] && !grid[b][s])
		{
		  grid[b][maxnseq]=0;
		  tot--;
		  break;
		}
	    }
	}
      tot=(tot*100)/S->nseq;
      if ( tot>=f)break;
    }
  if (tot<f)return NULL;
  
  fname=vtmpnam (NULL);
  fp=vfopen (fname, "w");
  for (a=0; a<maxnseq; a++)
    {
      if (count[a][2])
	{
	  fprintf ( fp, ">%s LIMIT: %d %%\n", FS->name[count[a][0]], f);
	  
	}
    }
  vfclose (fp);
  free_int (grid, -1); free_int (count, -1); 
  free_sequence (FS, FS->nseq);
  
  return get_fasta_sequence (fname, NULL);
}
