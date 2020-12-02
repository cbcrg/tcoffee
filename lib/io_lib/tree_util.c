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


  CL->DM=cl2distance_matrix (CL,NOALN,const_cast<char*>( (mode==NULL)?"ktup":mode), NULL, 0);

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
  
  
  list=(char***)vcalloc (A->nseq, sizeof (char***));
  nlist=(int*)vcalloc (A->nseq, sizeof (int));
  
  if ( n==0)return T;
  else if (n>1)
    {
      int l;
      char *buf;

      for (l=0,a=0; a< n; a++)l+=strlen (string[a]);
      buf=(char*)vcalloc ( 2*n+l+1, sizeof (char));
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

NT_node collapse_tree (NT_node T, Sequence *S, char *string)
{
  char *r, *p;
  int a;
  int collapse;

  if (!T) return NULL;
  if (!S)
    {
      S=tree2seq(T, NULL);
      T=recode_tree (T,S);
    }


  for (a=0; a<T->nseq; a++)
    {
      
      p=strstr((S->name[T->lseq[a]]),string);
	
      if (p && strm (p, string));
      else 
	{
	  collapse_tree (T->left, S, string);
	  collapse_tree (T->right, S, string);
	  return T;
	}
    }
  
  T->isseq=1;
  T->right=T->left=NULL;
  T->nseq=1;
  sprintf ( T->name, "%s",string);
  return T;
}

/*********************************************************************/
/*                                                                   */
/*                                   tree pruning                    */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
NT_node remove_leaf ( NT_node T);
NT_node prune_root (NT_node T);


NT_node prune_tree(NT_node T, Sequence *S)
{
  return (main_read_tree (prune_treeF(T, S, NULL)));
}

int check_treeF_seqF(char *tree, char *seq)
{
  int ret=0;
  
  Sequence*S=get_fasta_sequence (seq,NULL);
  NT_node T =main_read_tree(tree);
  ret=check_tree_seq(T,S);
  free_tree(T);
  free_sequence (S,-1);
  return ret;
}
int check_tree_seq(NT_node T, Sequence *S)
{
  NT_node *L=tree2seqnode_list (T,NULL);
  int a;
  int ntseq=0;
  char **tname=(char**)vcalloc (arrlen(L), sizeof (char*));
  while (L[ntseq])
    {
      if (name_is_in_hlist (L[ntseq]->name,S->name, S->nseq)==-1)
	{
	  HERE ("%s is not in seq", L[ntseq]->name);
	  return 0;
	}
      tname[ntseq]=L[ntseq]->name;
      ntseq++;
    }
  
  
  

  for (a=0; a<S->nseq; a++)
    if (name_is_in_hlist (S->name[a], tname,ntseq)==-1)
	{
	  HERE ("%s is not in seq", S->name[a]);
	  return 0;
	}
  vfree(L); vfree(tname);
  return 1;
}
  
char* prune_treeF (NT_node T, Sequence *S, char *tmp)
{
  NT_node *L;
  NT_node *R;
 
  int a;
  int **r;
  float *d;
  int special=0;
  NT_node B,C,P,GP;
  if (!tmp)tmp=vtmpnam (NULL);

  R=tree2node_list(T,NULL);
  a=0;
  while (R[a])
    {
      R[a]->order=a;
      a++;
    }
  
  r=declare_int (arrlen(R), 3);
  d=(float*)vcalloc (sizeof (float), arrlen(R));
  
  //This stores the original Node connectivity that defines the tree topology
  a=0;
  while (R[a])
    {
      if (R[a]->parent)r[a][0]=(R[a]->parent)->order;
      else r[a][0]=-1;
      
      if (R[a]->left  )r[a][1]=(R[a]->left  )->order;
      else r[a][1]=-1;
      
      if (R[a]->right )r[a][2]=(R[a]->right )->order;
      else r[a][2]=-1;
      
      d[a]=R[a]->dist;
      a++;
    }
  //T will now be pruned
  
  L=tree2seqnode_list (T,NULL);
  
  a=0;
  while (!special && L[a])
    {
      C=L[a];
      C->index=1;
      P=C->parent;
      GP=(P)?P->parent:NULL;
      
      
      if (P->right==C)B=P->left;
      else B=P->right;
      
      if (name_is_in_hlist (C->name,S->name, S->nseq)==-1)
	{
	  C->index=0;
	  if (!GP && B->isseq)
	    {
	      
	      B->parent=NULL;
	      special=1;
	    }
	  else if (!GP)
	    {
	      P->right=B->right;
	      P->left =B->left;
	      (B->right)->parent=P;
	      (B->left )->parent=P;
	      P->dist=0;
	    }
	  else
	    {
	      if (GP->left==P)GP->left=B;
	      else GP->right=B;
	      B->parent=GP;
	      B->dist+=P->dist;
	    }
	}
      a++;
    }
  
  print_newick_tree ((special)?B:T,tmp);
  
  //restore the original Tree
  a=0;
  while (R[a])
    {
      
      R[a]->dist=d[a];
      R[a]->parent=(r[a][0]==-1)?NULL:R[r[a][0]];
      R[a]->left  =(r[a][1]==-1)?NULL:R[r[a][1]];
      R[a]->right =(r[a][2]==-1)?NULL:R[r[a][2]];
      a++;
    }
  vfree(L); vfree(R);vfree(d);
  free_int (r,-1);
  
  return tmp;
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
  
  sl1=declare_int (10*S->nseq, S->nseq+1);
  sl2=declare_int (10*S->nseq, S->nseq+1);

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

   char *score_csv_file =(char*)vcalloc(200, sizeof (char));
   char *score_html_file =(char*) vcalloc(200, sizeof (char));
   char *hit_matrix_file =(char*) vcalloc(200, sizeof (char));
   char *hit_html_file =(char*) vcalloc(200, sizeof (char));
   char *tree_file =(char*) vcalloc(200, sizeof (char));

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

  poslist=(int*)vcalloc ( A->len_aln, sizeof (int));
  nl=0;
  fascore =(float*) vcalloc(A->len_aln, sizeof (float));

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
  NT_node *TreeArray =(NT_node*) vcalloc(nl, sizeof (NT_node));

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
      for (b=0,ax=0; ax<nl; ax++)
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
      ptree=(char*)vcalloc(100, sizeof (char));
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

	ffpHitScoreMatrix=(float**)vcalloc (nl, sizeof (float*));
	int i, j;
	for(i = 0; i < nl; i++)
      		ffpHitScoreMatrix[i]=(float*)vcalloc (nl-i, sizeof (float));

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

  myexit(EXIT_SUCCESS);return NULL;
}

NT_node aln2std_tree(Alignment *A, int ipara1, int ipara2, char *mode)
{
     Alignment *B;
     NT_node T;
     char *cpSet =(char*) vcalloc(100, sizeof (char));

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

    pos=(int*)vcalloc ( A->len_aln+1, sizeof (int));
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
     BTS=(Tree_sim*)vcalloc (1, sizeof (Tree_sim));

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

  TS1=(Tree_sim*)vcalloc (1, sizeof (Tree_sim));
  TS2=(Tree_sim*)vcalloc (1, sizeof (Tree_sim));

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
  TS1->rf=(((TS1->rf+TS2->rf))*10)/2;
  vfree (TS2);
  return TS1;
}
int print_node_list (NT_node T, Sequence *RS)
{
  NT_node *L;
  Sequence *S;
  int *nlseq2;

  S=tree2seq(T, NULL);
  L=tree2node_list (T, NULL);
  if (!RS)RS=S;
  nlseq2=(int*)vcalloc ( RS->nseq, sizeof (int));
  while (L[0])
    {
      int d,b;
      d=MIN(((L[0])->nseq), (S->nseq-(L[0])->nseq));
      fprintf ( stdout, "Bootstrap: %5.2f Depth: %5d Splits: ", (L[0])->bootstrap, d);

      for (b=0; b<RS->nseq; b++)nlseq2[b]='-';
      for (b=0; b<S->nseq; b++)
	{
	  int p;
	  p=name_is_in_list (S->name[b], RS->name, RS->nseq, 100);
	  if (p!=-1)nlseq2[p]=(L[0])->lseq2[b]+'0';
	}
      for (b=0; b<RS->nseq; b++) fprintf ( stdout, "%c", nlseq2[b]);
      fprintf (stdout, "\n");
      L++;
    }
  return 1;
}
Alignment * phylotrim1 (Alignment *A, char *reftree,int narg,char **nargl);


Alignment * phylotrim1 (Alignment *A, char *reftree,int narg,char **narglm, char *name);
Alignment * phylotrim2 (Alignment *A, char *reftree,int *l, int ph, int narglist, char **arglist, char *name);

Alignment * phylotrim (Alignment *A, NT_node RT, char *Ns,  char *treemode, char *tempfile)
{
  int N, a, b;
  char *tmptree;
  char **arglist;
  int narglist=0;
  int free_rt=0;
  char *name;

  name=(char *)vcalloc ( strlen (A->name[0])+strlen(A->name[A->nseq-1])+100, sizeof (char));
  sprintf (name, "%s%s%d", A->name[0], A->name[A->nseq-1], A->nseq);
  
  
  if (!Ns)N=1;
  else if (strstr (Ns, "split"))
    {
      N=-1;
    }
  else if (strstr (Ns, "%"))
    {
      sscanf (Ns, "%d%%", &N);
      N=(A->nseq*100)/N;
    }
  else 
    {
      N=atoi (Ns);
      N=MIN(N, ((A->nseq)-3));
    }
  
  arglist=(char**)vcalloc (2, sizeof (char*));
  arglist[narglist++]=treemode;
  if (tempfile)arglist[narglist++]=tempfile;
  
  tmptree=vtmpnam (NULL);
  if (!RT)
    {
      RT=main_aln2tree(A,narglist, arglist, name);
      free_rt=1;
    }
  
  if (!RT)
    printf_exit (EXIT_FAILURE, stderr, "ERROR: Could Not Compute the reference tree [FATAL:%s]\n", PROGRAM);
  else
    vfclose (print_tree (RT,"newick", vfopen (tmptree, "w")));
  
  
  if (N>0)
    {
      for (a=0; a<N; a++)
	{
	  A=phylotrim1(A, tmptree,narglist, arglist, name);
	}
    }
  else
    {
      Sequence *S;
      int ns=0;
      int **sl;
      FILE *fp;
      char * alnF=vtmpnam (NULL);
      
      S=tree2seq(RT, NULL);
      
      RT=recode_tree (RT, S);
      fp=vfopen (alnF, "w");
      for (a=0; a< S->nseq; a++)
	{
	  int i=name_is_in_list (S->name[a], A->name, A->nseq, 100);
	  if (i>=0)fprintf ( fp, ">%s\n%s\n", A->name[i], A->seq_al[i]);
	  else
	    printf_exit (EXIT_FAILURE, stderr, "ERROR: %s is not in the aln [FATAL:%s]\n",S->name[a], PROGRAM);
	}
      vfclose (fp);
      free_aln (A);
      A=main_read_aln (alnF, NULL);
      sl=declare_int (10*S->nseq, S->nseq+1);
      tree2split_list (RT, S->nseq, sl, &ns);
      
      //add the leafs to the split list
      for (a=0; a<A->nseq; a++)
	{
	  fprintf ( stdout, ">seq %s\n", A->name [a]);
	  for (b=0; b<A->nseq; b++)
	    {
	      sl[ns][b]=(a==b)?1:0;
	    }
	  ns++;
	}

      for (a=0; a< ns; a++)
	{
	  phylotrim2 (A, tmptree, sl[a], 0, narglist, arglist, name);
	  phylotrim2 (A, tmptree, sl[a], 1, narglist, arglist, name);
	}
      
      free_int (sl, -1);
    }
  if (free_rt)free_tree (RT);
  
  myexit (EXIT_SUCCESS);
}


  
Alignment * phylotrim2 (Alignment *A, char *reftree,int *l, int ph, int narg, char **argl, char *name)
{
  Alignment *NA=NULL;
  NT_node RT=NULL,NT=NULL;
  Tree_sim *T;
  int a, b, n, rf;
  FILE *fp;
  
  char *tmpaln=vtmpnam (NULL);
  
  
  for (n=0,a=0; a< A->nseq; a++)
    {
      n+=(l[a]==ph)?1:0;
    }
  
  if (n<3) return A;
  
  
  fp=vfopen (tmpaln, "w");
  for (a=0; a<A->nseq;a++){if (l[a]==ph)fprintf (fp, ">%s\n%s\n", A->name[a], A->seq_al[a]);}
  vfclose (fp);
  
  NA=main_read_aln (tmpaln, NULL);
  if ((NT=main_aln2tree (NA, narg, argl, name)))    
    {
      RT=main_read_tree (reftree);
      T=tree_cmp (RT, NT);
      rf=T->rf;
      vfree(T);
    }
  else
    {
      rf=-1;
    }
  
  fprintf (stdout, ">split ");
  for (a=0; a<A->nseq; a++){fprintf (stdout, "%d",(ph==1)?l[a]:1-l[a]);}
  fprintf ( stdout, " N: %4d R: %4d RF: %4d\n", n,A->nseq-n, rf);
  
  
  
  free_tree (RT);
  free_tree (NT);
  free_aln (NA);
  
  return A;
}
Alignment * phylotrim1 (Alignment *A, char *reftree, int narg, char **argl, char *name)
{
  Alignment *NA;
  NT_node RT,NT;
  Tree_sim *T;
  int a, b, n;
  FILE *fp;
  int seq, brf;
  char *alnF=vtmpnam (NULL);

  seq=brf=-1;
  
  for (a=0; a< A->nseq; a++)
    {
      fp=vfopen (alnF, "w");
      for (b=0; b<A->nseq; b++)
	{
	  if ( a!=b)fprintf (fp, ">%s\n%s\n", A->name[b], A->seq_al[b]);
	}
      vfclose (fp);
      
      NA=main_read_aln (alnF, NULL);
      if (NT=main_aln2tree (NA, narg, argl, name))
	{
	  RT=main_read_tree (reftree);
	  T=tree_cmp (RT, NT);
	  if (T->rf>brf)
	    {
	      brf=T->rf;
	      seq=a;
	    }
	  vfree (T);
	  free_aln (NA);
	  free_tree (RT);
	  free_tree (NT);
	}
      else
	{
	  free_aln (NA);
	}
    }
  fprintf ( stdout, ">%s RF: %d\n", A->name[seq], brf);
  
  fp=vfopen (alnF, "w");
  for (b=0; b<A->nseq; b++)
     {
       if (b!=seq)fprintf (fp, ">%s\n%s\n", A->name[b], A->seq_al[b]);
     }
   vfclose (fp);
   free_aln (A);
   return main_read_aln (alnF, NULL);
}



NT_node main_compare_trees_list ( NT_node RT, Sequence *S, FILE *fp)
{
  Tree_sim *T;
  NT_node *TL;
  Sequence *RS;
  int a;

  T=(Tree_sim*)vcalloc (1, sizeof (Tree_sim));
  RS=tree2seq(RT, NULL);

  TL=read_tree_list (S);

  reset_boot_tree (RT,0.0001);
  for (a=0; a<S->nseq; a++)
    {
      TL[a]=prune_tree(TL[a],RS);
      TL[a]=recode_tree (TL[a],RS);
      new_compare_trees (RT, TL[a], RS->nseq,T);

    }

  vfree (T);
  return RT;
}
NT_node main_compare_trees ( NT_node T1, NT_node T2, FILE *fp)
{
  Tree_sim *T;

  T=tree_cmp (T1, T2);
  fprintf ( fp, "\n#tree_cmp|T: %.2f W: %.2f L: %.2f RF: %.2f N: %d S: %d", T->uw, T->w, T->d, (float)T->rf/10, T->n, T->nseq);
  fprintf ( fp, "\n#Distance obtained by averaging T1 vs T2 and T2 vs T1");
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
	  T1->bootstrap++;
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
NT_node aln2trmsd_tree (Alignment *A, char *temp);
NT_node aln2phyml_tree (Alignment *A);
NT_node main_aln2tree (Alignment *A, int n, char **arglist, char *name)
{
  NT_node T;
  char *CacheAln;
  char *CacheTree;
  
  

  if (!name || strm (name, "nocache"))
    {
      T=tree_compute (A, n, arglist);
      return T;
    }
  else
    {
      char *tmpA=vtmpnam(NULL);
      long int HV;
      
      output_fasta_aln (tmpA,A);
      HV=hash_file (tmpA);
      
      
      CacheAln=(char* )vcalloc (strlen (get_cache_dir ())+100+strlen (arglist[0])+100, sizeof (char));
      CacheTree=(char*)vcalloc (strlen (get_cache_dir ())+100+strlen (arglist[0])+100, sizeof (char));
      sprintf (CacheAln, "%s/%ld.%s.aln", get_cache_dir(),HV, arglist[0]);
      sprintf (CacheTree,"%s/%ld.%s.nwk", get_cache_dir(),HV, arglist[0]);
      
      //HERE ("\nALN: %s (%d)\nTree %s(%d)\nDiff: %d", CacheAln, file_exists (NULL,CacheAln),CacheTree,file_exists (NULL,CacheTree),file2diff (CacheAln, tmpA));
      
      if (!file_exists (NULL,CacheAln) || !file_exists (NULL,CacheTree)|| file2diff (CacheAln, tmpA))
	{
	  //HERE ("Compute");
	  printf_system ("cp %s %s", tmpA, CacheAln);
	  T=tree_compute (A, n, arglist);
	  if (T)vfclose (print_tree (T, "newick", vfopen (CacheTree, "w")));
	  free_tree (T);
	}
      else
	{
	  //HERE ("Cache");
	  ;
	}
      T=main_read_tree (CacheTree);
      vfree (CacheAln);
      vfree (CacheTree);
      return T;
    }
}
NT_node tree_compute ( Alignment *A, int n, char ** arg_list)
{
  if (n==0 || strm (arg_list[0], "cw"))
    {
      return compute_cw_tree (A);
    }
  else if (strm (arg_list[0], "phyml"))
    return aln2phyml_tree (A); 
  else if ( strm (arg_list[0], "fj"))
    {
      return compute_fj_tree ( NULL, A, (n>=1)?atoi(arg_list[1]):8, (n>=2)?arg_list[2]:NULL);
    }
  
  else if (strm (arg_list[0], "kmeans"))return aln2km_tree (A, NULL, 1);
  else if (strm (arg_list[0], "trmsd"))return aln2trmsd_tree (A,arg_list[1]); 
  else if ( ( strm (arg_list[0], "nj")))
    {
      return compute_std_tree (A, n, arg_list);
    }
  else if ( ( strm (arg_list[0], "upgma")))
    {
      return compute_std_tree_2(A,NULL, "_TMODE_upgma");
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

NT_node aln2trmsd_tree (Alignment *A, char *temp)
{
  char *fname=vtmpnam (NULL);
  char *tree=vtmpnam (NULL);
  NT_node T=NULL;


  A->residue_case=KEEP_CASE;
  output_fasta_aln (fname, A);
  
    
  printf_system ("t_coffee -other_pg trmsd -aln %s -template_file %s -outfile %s -quiet %s", fname, temp, tree,  TO_NULL_DEVICE);
  //printf_system ("t_coffee -other_pg trmsd -aln %s -template_file %s -outfile %s ", fname, temp, tree);
  if ( check_file_exists (tree))
    T=main_read_tree (tree);
  if (!T) exit (0);
  
  return T;
}
 
NT_node aln2phyml_tree (Alignment *A)
{
  char *type;
  char *tmpaln=vtmpnam (NULL);
  char *name;
  NT_node RT;

  A=get_aln_type (A);
  type=(A->S)->type;
  output_phylip_aln (tmpaln, A, "w");
  
  /*
  if ( strstr (type, "DNA") || strstr (type, "RNA"))
    {
      printf_system ("phyml -i %s -d nt ", tmpaln);
    }
  else
    printf_system ("phyml -i %s -d aa ", tmpaln);
  */
  


  if ( strstr (type, "DNA") || strstr (type, "RNA"))
    {
      printf_system ("phyml -i %s -d nt --quiet %s ", tmpaln, TO_NULL_DEVICE);
    }
  else
    {
      printf_system ( "phyml -i %s -o tlr --rand_start --n_rand_starts 5 -d aa -b -1 -s BEST --quiet %s", tmpaln, TO_NULL_DEVICE);
      //printf_system ("phyml -i %s -d aa --quiet  %s", tmpaln,TO_NULL_DEVICE);
    }
    
  name=(char *)vcalloc ( strlen (tmpaln)+100, sizeof (char));
  
  sprintf (name,"%s_phyml_stats.txt",tmpaln);
  if ( check_file_exists (name))delete_file (name);
  else printf_exit ( EXIT_FAILURE,stderr, "ERROR: Failed to compute a PhyMl tree [FATAL:%s]\n", PROGRAM);
  
  sprintf (name,"%s_phyml_tree.txt",tmpaln);
  if ( !check_file_exists (name))
    printf_exit ( EXIT_FAILURE,stderr, "ERROR: Failed to compute a PhyMl tree [FATAL:%s]\n", PROGRAM);
  
  RT=main_read_tree (name);
  delete_file (name);
  vfree (name);
  return RT;
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

  aln_file=vtmpnam(NULL);
  ungaped_aln_file=vtmpnam (NULL);
  scored_aln_file=vtmpnam (NULL);
  scored_aln_file=vtmpnam(NULL);
  filtered_aln_file=vtmpnam(NULL);



  output_clustal_aln (aln_file, A);
  /* 1: remove columns with too many gaps*/
  printf_system ("t_coffee -other_pg seq_reformat -in %s -action +rm_gap %d -output clustalw > %s", aln_file,fraction_gap, ungaped_aln_file);

  /* 2: evaluate the alignment*/

  printf_system ("t_coffee -other_pg seq_reformat -in %s -action +evaluate %s -output clustalw > %s", ungaped_aln_file,(mode)?mode:"categories", scored_aln_file);


  /*3 extract the high scoring columns*/
  printf_system("t_coffee -other_pg seq_reformat -in %s -struc_in %s -struc_in_f number_aln -action +use_cons +keep '[%d-8]' +rm_gap -output clustalw > %s", ungaped_aln_file, scored_aln_file,n, filtered_aln_file);


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
void reset_node_count (NT_node root)
{
  //MorrisTraversal: No Recursion, No stack
  NT_node current,pre;
  
  if( root== NULL)
    return; 
 
  current = root;
  root->visited=0;
  
  while(current)
    {                 
      if(!current->left)
	{
	  if (current)current->visited=0;
	  current = current->right;   
	  if ( current) current->visited=0;
	}    
      else
	{
	  pre = current->left;
	  if (pre)pre->visited=0;
	  
	  while(pre->right && pre->right != current)
	    {
	      pre = pre->right;
	      if (pre)pre->visited=0;
	    }
	  
	  
	  if(!pre->right )
	    {
	      pre->right = current;
	      if ( current)current->visited=0;
	      current = current->left;
	      if (current) current->visited=0;
	    }
	  
	  else 
	    {
	      pre->right = NULL;
	      current = current->right;     
	      if (current) current->visited=0;
	    } 
	} 
    } 
}


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
  R->right=realloc_tree (R->right,n);
  R->left=realloc_tree (R->left,n);
  R->bot=realloc_tree (R->bot,n);

  R->lseq=(int*)vrealloc (R->lseq, n*sizeof (int));
  R->lseq2=(int*)vrealloc (R->lseq2, n*sizeof (int));
  return R;
}

NT_node reset_boot_tree ( NT_node R, int n)
{
  if ( !R)return R;
  R->right=reset_boot_tree (R->right,n);
  R->left=reset_boot_tree (R->left,n);
  R->bot=reset_boot_tree (R->bot,n);
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
  R->right=reset_dist_tree (R->right,n);
  R->left=reset_dist_tree (R->left,n);
  R->bot=reset_dist_tree (R->bot,n);

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

NT_node free_tree ( NT_node root)
{
  NT_node stack, current, pre;
  
  if (!root)return NULL;
  reset_node_count (root);
  
  stack=root;
  stack->bot=NULL;
  current = root;
  root->visited=1;
  
  while(current)
    {                 
      if(!current->left)
	{
	  if (current)
	    {
	      if (!current->visited){ current->bot=stack;stack=current;}
	      current->visited=1;
	    }
	  current = current->right;   
	  if ( current) 
	    {
	      if (!current->visited){current->bot=stack;stack=current;}
	      current->visited=1;
	    }
	}    
      else
	{
	  pre = current->left;
	  if (pre)
	    {
	      if (!pre->visited){pre->bot=stack;stack=pre;}
	      pre->visited=1;
	    }
	  
	  while(pre->right && pre->right != current)
	    {
	      pre = pre->right;
	      if (pre)
		 {
		   if (!pre->visited){pre->bot=stack;stack=pre;}
		   pre->visited=1;
		 }
	    }
	  
	  
	  if(!pre->right )
	    {
	      pre->right = current;
	      if ( current)
		{
		  if (!current->visited){current->bot=stack;stack=current;}
		  current->visited=1;
		}
	      current = current->left;
	      if (current)
		{
		  if (!current->visited){current->bot=stack;stack=current;}
		  current->visited=1;
		}
	    }
	  
	  else 
	    {
	      pre->right = NULL;
	      current = current->right;     
	      if (current) 
		{
		  if (!current->visited){current->bot=stack;stack=current;}
		  current->visited=1;
		}
	    } 
	} 
    } 
  while ( stack)
    {
      NT_node tofree=stack;
      stack->visited++;
      stack=stack->bot;
      free_tree_node (tofree);
    }
}
  

NT_node rec_free_tree ( NT_node R)
{
  if ( !R)return R;
  if(R->right)
  	  R->right=rec_free_tree (R->right);
  if(R->left)
  	  R->left=rec_free_tree (R->left);
  if(R->bot)
  	  R->bot=rec_free_tree (R->bot);
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

int decode_seq_in_tree (NT_node R, char **name)
{
  //seq are expected to be named  1--N in the order of **name
  //returns the number of sequences effectively decoded
  int t=0;
  if (!R) return 0;
  if (R->leaf!=1)
    {
      t+=decode_seq_in_tree (R->right,name);
      t+=decode_seq_in_tree (R->left, name);
    }
  else
    {
      int s=atoi (R->name);
      vfree (R->name);R->name=(char*)vcalloc ( strlen (name[s-1])+1, sizeof (char));
      sprintf (R->name, "%s", name[s-1]);
      t=1;

    }
  return t;
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
char * node2seq_file (NT_node N, Sequence *S, int cache, char *file)
{
  int a;
  FILE *fp;
  //cache==0: above, ==1: below
  if (!file) file=vtmpnam(NULL);
  fp=vfopen (file, "w");
  for (a=0; a<S->nseq; a++)
    {
      if (N->lseq2[a]==cache)
	fprintf (fp, ">%s\n%s\n", S->name[a],S->seq[a]);
    }
  vfclose (fp);
  return file;
}
Sequence * tree2seq (NT_node p, Sequence *S)
{
  
  if ( !S)
    {
      S=declare_sequence (10, 10, tree2nseq(p));
      S->nseq=0;
    }
  reset_node_count (p);
  while (p)
    {
      int x=++(p->visited);
      if (!p->isseq)
	{
	  if (x==1)
	    {
	      p=p->right;
	    }
	  else if (x==2)
	    {
	      p=p->left;
	    }
	  else if (x==3 && !p->parent)
	    {
	      p=p->parent;
	    }
	}
      
      if (x>=3 || p->isseq)
	{
	  if (p)
	    {
	      if (p->isseq)sprintf (S->name[S->nseq++], "%s",p->name);
	    }
	  if (x>3)p=NULL;
	  else if (p) p=p->parent;
	}
    }
  return S;
}



int seqindex2seqname4tree (NT_node root, Sequence *S)
{
  //index goes 1 to N
  NT_node current,pre;
  if ( !root)return 0;
  int rv=0;
  reset_node_count (root);
  current = root;
  root->visited=1;
  
  while(current != NULL)
    {                 
      if(!current->left)
	{
	  if (current && current->isseq && !current->visited)
	    {
	      int i=atoi (current->name)-1;
	      if (i<S->nseq)current->name=csprintf (current->name, "%s", S->name[i]);
	      else printf_exit ( EXIT_FAILURE,stderr, "\nERROR: Sequence %s is in tree but not in source sequence file [FATAL]", current->name);
	      rv++;
	    } 
	  current->visited=1;
	  current = current->right;  
	}    
      else
	{
	  pre = current->left;
	  
	  if (!pre->visited)
	    {
	      if (pre && !pre->visited && pre->isseq)
		{
		  int i=atoi (pre->name)-1;
		  if (i<S->nseq)current->name=csprintf (pre->name, "%s", S->name[i]);
		  else printf_exit ( EXIT_FAILURE,stderr, "\nERROR: Sequence %s is in tree but not in source sequence file [FATAL]",pre->name);
		  rv++;
		}
	      pre->visited=1;
	    }
	  while(pre->right && pre->right != current)
	    {
	      pre = pre->right;
	      if (pre)
		{
		  if (pre && !pre->visited && pre->isseq)
		    {
		       int i=atoi (pre->name)-1;
		       if (i<S->nseq)pre->name=csprintf (pre->name, "%s", S->name[i]);
		       else printf_exit ( EXIT_FAILURE,stderr, "\nERROR: Sequence %s is in tree but not in source sequence file [FATAL]",pre->name);
		       rv++;
		    }
		  pre->visited=1;
		}
	    }
	  
	  
	  if(!pre->right )
	    {
	      pre->right = current;
	      current = current->left;
	    }
	  
	  else 
	    {
	      pre->right = NULL;
	      current = current->right;      
	    } 
	} 
    } 
  return rv;
}
int seqname2seqindex4tree (NT_node root, Sequence *S)
{
  //index goes 1 to N
  NT_node current,pre;
  if ( !root)return 0;
  int rv=0;
  reset_node_count (root);
  current = root;
  root->visited=1;
  
  while(current != NULL)
    {                 
      if(!current->left)
	{
	  if (current && current->isseq && !current->visited)
	    {
	      int i=name_is_in_hlist (current->name, S->name, S->nseq);
	      if (i!=-1)current->name=csprintf (current->name, "%d", i+1);
	      else printf_exit ( EXIT_FAILURE,stderr, "\nERROR: Sequence %s is in tree but not in source sequence file [FATAL]", current->name);
	      rv++;
	    } 
	  current->visited=1;
	  current = current->right;  
	}    
      else
	{
	  pre = current->left;
	  
	  if (!pre->visited)
	    {
	      if (pre && !pre->visited && pre->isseq)
		{
		  int i=name_is_in_hlist (pre->name, S->name, S->nseq);
		  if (i!=-1)current->name=csprintf (pre->name, "%d", i+1);
		  else printf_exit ( EXIT_FAILURE,stderr, "\nERROR: Sequence %s is in tree but not in source sequence file [FATAL]", pre->name);
		  rv++;
		}
	      pre->visited=1;
	    }
	  while(pre->right && pre->right != current)
	    {
	      pre = pre->right;
	      if (pre)
		{
		  if (pre && !pre->visited && pre->isseq)
		    {
		       int i=name_is_in_hlist (pre->name, S->name, S->nseq);
		       if (i!=-1)pre->name=csprintf (pre->name, "%d", i+1);
		       else printf_exit ( EXIT_FAILURE,stderr, "\nERROR: Sequence %s is in tree but not in source sequence file [FATAL]", pre->name);
		       
		       rv++;
		    }
		  pre->visited=1;
		}
	    }
	  
	  
	  if(!pre->right )
	    {
	      pre->right = current;
	      current = current->left;
	    }
	  
	  else 
	    {
	      pre->right = NULL;
	      current = current->right;      
	    } 
	} 
    } 
  return rv;
}

int seqname2seqindexname4tree (NT_node root, Sequence *S)
{
  //index goes 1 to N
  NT_node current,pre;
  if ( !root)return 0;
  int rv=0;
  reset_node_count (root);
  current = root;
  root->visited=1;
  
  while(current != NULL)
    {                 
      if(!current->left)
	{
	  if (current && current->isseq && !current->visited)
	    {
	      int i=name_is_in_hlist (current->name, S->name, S->nseq);
	      if (i!=-1)current->name=csprintf (current->name, "%d_%s", i+1, current->name);
	      else printf_exit ( EXIT_FAILURE,stderr, "\nERROR: Sequence %s is in tree but not in source sequence file [FATAL]", current->name);
	      rv++;
	    } 
	  current->visited=1;
	  current = current->right;  
	}    
      else
	{
	  pre = current->left;
	  
	  if (!pre->visited)
	    {
	      if (pre && !pre->visited && pre->isseq)
		{
		  int i=name_is_in_hlist (pre->name, S->name, S->nseq);
		  if (i!=-1)current->name=csprintf (pre->name, "%d_%s", i+1, pre->name);
		  else printf_exit ( EXIT_FAILURE,stderr, "\nERROR: Sequence %s is in tree but not in source sequence file [FATAL]", pre->name);
		  rv++;
		}
	      pre->visited=1;
	    }
	  while(pre->right && pre->right != current)
	    {
	      pre = pre->right;
	      if (pre)
		{
		  if (pre && !pre->visited && pre->isseq)
		    {
		       int i=name_is_in_hlist (pre->name, S->name, S->nseq);
		       if (i!=-1)pre->name=csprintf (pre->name, "%d_%s", i+1, pre->name);
		       else printf_exit ( EXIT_FAILURE,stderr, "\nERROR: Sequence %s is in tree but not in source sequence file [FATAL]", pre->name);
		       
		       rv++;
		    }
		  pre->visited=1;
		}
	    }
	  
	  
	  if(!pre->right )
	    {
	      pre->right = current;
	      current = current->left;
	    }
	  
	  else 
	    {
	      pre->right = NULL;
	      current = current->right;      
	    } 
	} 
    } 
  return rv;
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
int tree2nnode (NT_node p)
{
  int n=0;
  NT_node in=p;
  reset_node_count (p);
  while (p)
    {
      if (!p->visited)n++;
      int x=++(p->visited);

      if (!p->isseq)
	{
	  if (x==1)
	    {
	      p=p->right;
	    }
	  else if (x==2)
	    {
	      p=p->left;
	    }
	  else if (x==3 && !p->parent)
	    {
	      p=p->parent;
	    }
	}
      
      if (x>=3 || p->isseq)
	{
	  if (p&& p->isseq && !p->visited){n++;}
	  if (x>3)p=NULL;
	  else if (p) p=p->parent;
	}
    }
  reset_node_count (in);
  return n;
}

int tree2nleaf (NT_node root)
{
  //MorrisTraversal: No Recursion, No stack
  NT_node current,pre;
  int nleaf=0;
  
  reset_node_count (root);
  
  if( root== NULL)
    return nleaf; 
 
  current = root;
  root->visited=1;
  
  while(current)
    {                 
      if(!current->left)
	{
	  if (!current->visited && current->isseq)nleaf++;
	  current->visited=1;
	  current = current->right;      
	}    
      else
	{
	  pre = current->left;
	  
	  if (!pre->visited && pre->isseq)nleaf++;
	  pre->visited=1;
	  
	  while(pre->right && pre->right != current)
	    {
	      pre = pre->right;
	      if (pre)
		{
		  if (!pre->visited && pre->isseq)nleaf++;
		  pre->visited=1;
		}
	    }
	  
	  
	  if(!pre->right )
	    {
	      pre->right = current;
	      current = current->left;
	    }
	  
	  else 
	    {
	      pre->right = NULL;
	      current = current->right;      
	    } 
	} 
    } 
  return nleaf;
}


int rec_tree2nleaf (NT_node R)
{
  if ( !R)return 0;
  else if (R->isseq==1){R->nseq=1;return 1;}
  else if (R->right==NULL && R->left==NULL && R->bot==NULL){R->leaf=1; return 1;}
  else
    {
      int n=0;
      n+=tree2nleaf (R->right);
      n+=tree2nleaf (R->left);
      n+=tree2nleaf (R->bot);

      R->leaf=n;
      R->nseq=n;
      return n;
    }
}
NT_node tree2deepest_nodeR(NT_node T, NT_node *DN, int *depth);
NT_node tree2deepest_node(NT_node T)
{
  NT_node DN;
  int depth=0;

  return tree2deepest_nodeR(T, &DN, &depth);
}
NT_node tree2deepest_nodeR(NT_node T, NT_node *DN, int *depth)
{
  //tree must be recoded
  if (!T)return NULL;
  int a;
  int n1, n2, n;
  
  n1=n2=n=0;
  for (a=0; a<T->nseq; a++)
    {
      if (T->lseq2[a])n1++;
      else n2++;
    }
  n=MIN(n1, n2);
  if (n>depth[0]){depth[0]=n; DN[0]=T;}
  tree2deepest_nodeR(T->left , DN, depth);
  tree2deepest_nodeR(T->right, DN, depth);
  return DN[0];
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

  string=(char*)vcalloc (count_n_char_in_file(fname)+1, sizeof (char));

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

  if (!EMPTY)EMPTY=(NT_node)vcalloc (1, sizeof (NT_node));
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
FILE * print_tree ( NT_node T, char *format,FILE *fp)
{
  
  tree2file (T,NULL, format, fp);
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
      vfclose (tree2file (T,NULL,"newick", vfopen (f, "w")));
      return file2string (f);
    }
}


FILE * tree2file ( NT_node T, Sequence *S,char *format,FILE *fp)
{
  


  /*
    if (!S)
    {
      S=tree2seq(T, NULL);
      recode_tree (T, S);
      free_sequence (S, -1);
    }
  else
    recode_tree (T, S);
  */
  
  if ( format && strm (format, "binary"))
    fp=display_tree ( T,S->nseq, fp);
  else if ( strm (format, "shuffle_newick"))
    {
      vsrand(0);
      fp=no_rec_print_tree_shuffle (T, fp);
      fprintf ( fp, ";\n");
    }
  else if ( ! format || strm2 (format, "newick_tree","newick"))
    {
      fp=no_rec_print_tree (T, fp);
      fprintf ( fp, ";\n");
    }
  else if ( strm (format, "mafftnewick"))
    {
      if (!S)printf_exit ( EXIT_FAILURE,stderr, "\nERROR: tree2file/mafft requires sequence file [FATAL]\n");
      seqname2seqindex4tree (T, S);
      fp=no_rec_print_tree (T, fp);
      fprintf ( fp, ";\n");
    }
  else
    {
      fprintf ( stderr, "\nERROR: %s is an unknown tree format [FATAL:%s]\n", format, PROGRAM);
      myexit (EXIT_FAILURE);
    }
  return fp;
}


char * print_newick_tree ( NT_node T, char *name)
{
  FILE *fp;
  
  
  fp=vfopen (name, "w");

  
  fp=no_rec_print_tree (T, fp);


  
  fprintf (fp, ";\n");
  vfclose (fp);
  return name;
}

NT_node tree2shuffle (NT_node p)
{
  NT_node *L;
  int a=0;
  L=tree2node_list (p,NULL);
  
  while (L[a])
    {
      int order=rand()%2;
      if (order)
	{
	  NT_node buf=L[a]->right;
	  L[a]->right=L[a]->left;
	  L[a]->left=buf;
	}
      a++;
    }
  vfree (L);
  return p;
}

int newick2random4dpa (Sequence *S, NT_node p,int n, int ntrees)
{
  NT_node *L;
  int a, nseq;
  char **names;
  int **sorted;
  KT_node K;
  float *w;
  int sim;
  Sequence *S2;
  Alignment *A2;
  char * tmp=vtmpnam (NULL);
 

  int nt;
 
  L=tree2seqnode_list (p,NULL);
  nseq=0;
  while (L[nseq])nseq++;
  names=(char**)vcalloc (nseq, sizeof (char*));
  sorted=declare_int (nseq, 2);
  
  w=seq2dpa_weight (S, "longuest");
  for (nt=0; nt<ntrees; nt++)
    {
      for (a=0; a<nseq; a++)
	{
	  names[a]=L[a]->name;
	  sorted[a][0]=a;
	  sorted[a][1]=rand()%nseq+1;
	}
      sort_int (sorted,2,1, 0,nseq-1);
      
      for (a=0; a<nseq; a++)
	{
	  L[a]->name=names[sorted[a][0]];
	}
      p=node2master (p, S, w);
      K=tree2ktree  (p,p, S, n);
      seq2mafft_aln_file (K->seqF, tmp);
      A2=main_read_aln(tmp, NULL);
      sim=aln2sim(A2, "idmat");
      HERE ("Tree: %d -- Sim=%d", nt+1, sim);
      fprintf (stdout, "SIM: %d --", sim);
      no_rec_print_tree (p,stdout);
      fprintf (stdout, ";\n");
      free_ktree(K);
    }
  return 1;
}

FILE * no_rec_print_tree_randomize ( NT_node p, FILE *fp)
{
  NT_node *L;
  int a, nseq;
  char **names;
  int **sorted;
  
  L=tree2seqnode_list (p,NULL);
  nseq=0;
  while (L[nseq])nseq++;
  names=(char**)vcalloc (nseq, sizeof (char*));
  sorted=declare_int (nseq, 2);
  
  for (a=0; a<nseq; a++)
    {
      names[a]=L[a]->name;
      sorted[a][0]=a;
      sorted[a][1]=rand()%nseq+1;
    }
  sort_int (sorted,2,1, 0,nseq-1);
  
  for (a=0; a<nseq; a++)
    {
      L[a]->name=names[sorted[a][0]];
    }
  free_int (sorted, -1);
  vfree (names);
  vfree (L);
  return no_rec_print_tree (p,fp);
}
      
FILE * no_rec_print_tree_shuffle ( NT_node p, FILE *fp)
{
  p=tree2shuffle (p);
  return no_rec_print_tree (p,fp);
}

FILE * no_rec_print_tree ( NT_node p, FILE *fp)
{
  int ns=tree2nleaf(p);
  if ( ns==1)
    {
      
      if (p->name)fprintf (fp, "(%s:1.000)", p->name);
      else if ( p->right)fprintf (fp, "(%s:1.000)", (p->right)->name);
      else if ( p->left)fprintf (fp, "(%s:1.000)", (p->left)->name);
      return fp;
    }

  reset_node_count (p);
  
  while (p)
    {
      int x=++(p->visited);
      
      if (!p->isseq)
	{
	  if (x==1)
	    {
	      fprintf ( fp, "(");
	      p=p->right;
	    }
	  else if (x==2)
	    {
	      fprintf ( fp, ",");
	      p=p->left;
	    }
	  else if (x==3 && !p->parent)
	    {
	      fprintf ( fp, ")");
	      p=p->parent;
	    }
	  else if (x>=3)
	    {
	      fprintf ( fp, ")");
	    }
	}
      
      if (x>=3 || p->isseq)
	{
	  if (p)
	    {
	      if (p->isseq)fprintf ( fp, "%s:%.5f",p->name,p->dist);
	      else 
		{
		  if ( p->bootstrap!=0)fprintf (fp, "%d", (int)p->bootstrap);
		  fprintf ( fp, ":%.5f",p->dist);
		}
	    }
	  if (x>3)p=NULL;
	  else if (p) p=p->parent;
	}
    }
  return fp;
}
void reset_tree_distances( NT_node p, float d)
{
  int ns=tree2nleaf(p);
  if ( ns==1)
    {
      p->dist=d;
      return;
    }

  reset_node_count (p);
  
  while (p)
    {
      int x=++(p->visited);
      
      if (!p->isseq)
	{
	  p->dist=1;
	  if (x==1)
	    {
	      p=p->right;
	      if (p)p->dist=d;
	    }
	  else if (x==2)
	    {
	      p=p->left;
	      if (p)p->dist=d;
	    }
	  else if (x==3 && !p->parent)
	    {
	      p=p->parent;
	      if (p)p->dist=d;
	    }
	  else if (x>=3){;}
	}
      if (x>=3 || p->isseq)
	{
	  if (p)p->dist=d;
	  if (x>3)p=NULL;
	  else if (p) 
	    {
	      p=p->parent;
	      p->dist=1;
	    }
	}
    }
  reset_node_count (p);
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
	    tree_names=(char**)vcalloc ( tnseq, sizeof (char*));
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


	name=(char*)vcalloc ( MAXNAMES+1, sizeof (char));
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
  seq2keep=(int*)vcalloc (A->nseq, sizeof (int));
  tname=(int*)vcalloc (T->nseq, sizeof (int));
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

  list=(int*)vcalloc (A->nseq, sizeof (int));
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
  used=(int*)vcalloc (A->nseq, sizeof(int));
  used[seq1]=1;

  if (find_seq_chain ( A, sim,used,seq1,seq1, seq2,1,limit, max_chain, &nseq))
    {
      list=(int*)vcalloc (A->nseq, sizeof (int));
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
      pw_sim=(int***)declare_arrayN(3, sizeof (int), A->nseq, A->nseq, 3);
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
int scan_name_and_dist ( FILE *fp, char *name, float *dist);
 

Sequence *get_nexus (char*file)
{
  char *fasta=vtmpnam(NULL);
  char *buf=NULL;
  FILE *fpin;
  FILE *fpout;
  int i, max_i=0, nt=0;
  char **lu=NULL;
  
  
  fpin =vfopen (file, "r");
  fpout=vfopen (fasta, "w");
  
  while ((buf=vfgets (buf,fpin)))
    {
      if (strstr (buf, "Translate"))
	{
	  int n=0;
	  while((buf=vfgets (buf, fpin)) && !strstr(buf, ";"))
	    {
	      char *buf2=(char*)vcalloc (strlen(buf)+1, sizeof (char));
	      buf=substitute (buf, ",", "");
	      sscanf (buf, "%d %s", &i, buf2);

	      if (i>=max_i)
		{
		  max_i+=1000;
		  lu=(char**)vrealloc (lu, sizeof (char*)*max_i);
		}
	      lu[i]=buf2;
	      if (i>n)n=i;
	    }
	  
	  while((buf=vfgets (buf, fpin)) && !strstr(buf, "tree"));
	  if (buf)
	    {
	      char *tree=(char*)vcalloc ( strlen (buf)+1, sizeof (char));
	      char index[5];
	      sscanf (buf, "tree PAUP_1 = [&U] %s", tree);
	      
	      for (i=n; i>0; i--)
		{
		  if (lu[i])
		    {
		      sprintf (index, "%d", i);
		      tree=substitute (tree, index, lu[i]);
		      vfree (lu[i]);lu[i]=NULL;
		    }
		}
	      fprintf (fpout, ">tree_%d\n%s", ++nt,tree);
	      vfree (tree);
	    }
	}
    }
  vfclose (fpin);
  vfclose (fpout);
  vfree (lu);
  return get_fasta_tree (fasta, NULL);
}
Sequence * get_treelist (char *fname)
{
  FILE *fpin;
  FILE *fpout;
  int n=0;
  char *buf=NULL;
  char *seq=vtmpnam (NULL);
  
  
  fpin=vfopen (fname, "r");
  fpout=vfopen(seq, "w");
  while ((buf=vfgets(buf, fpin)))
    {
      if (buf[0]!='(' && check_file_exists (buf))
	{
	  char *tree=file2string(buf);
	  fprintf (fpout, ">Tree_%d\n%s", ++n,tree);
	  vfree (tree);
	}
      else
	fprintf (fpout, ">Tree_%d\n%s", ++n,buf);
    }
  vfclose (fpin);
  vfclose (fpout);
  return get_fasta_tree(seq, NULL);
}
Sequence*get_fasta_tree (char *fname, char *comment_out)
{
  Sequence *LS;
    char *buffer;
    FILE *fp;
    int a;

    int   c;
    char *name;
    int clen=0;
    int current=0;
    int p=0;
    int max;
    int max_len_seq=0;
    int min_len_seq=0;
    int nseq=0, l=0;




    int *sub;

    buffer=(char*)vcalloc (1000, sizeof (char));
    name=(char*)vcalloc ( 100, sizeof (char));

    nseq=count_n_char_x_in_file(fname, '>');
    min_len_seq=max=count_n_char_in_file(fname);
    sub=(int*)vcalloc (max+1, sizeof (int));

    fp=vfopen (fname,"r");


    c=fgetc(fp);
    while (c!=EOF)
	 	{
		 if (c=='>')
			{
			fscanf_seq_name (fp,name);
			while ((c=fgetc(fp))!='\n' && c!=EOF);
			while ((c=fgetc(fp))!='>' && c!=EOF)
			  if (isgraph(c))
			    clen++;
			 max_len_seq=(clen> max_len_seq)?clen: max_len_seq;
			 min_len_seq=(clen< min_len_seq)?clen: min_len_seq;
			 clen=0;
			}
		else
		    c=fgetc (fp);

		}

    vfclose (fp);
    LS=declare_sequence (  min_len_seq,  max_len_seq,nseq);

    LS->nseq=nseq;

    fp=vfopen (fname,"r");
    current=0;
    c=fgetc(fp);
    while (c!=EOF)
		{
	 	if (c=='>')
			{

			fscanf_seq_name (fp,LS->name[current]);
			l=strlen ( LS->name[current]);
			if ( LS->name[current][l-1]==','||LS->name[current][l-1]==';')LS->name[current][l-1]='\0';
			LS->name[current]=translate_name ( LS->name[current]);
			a=0;
			while ((c=fgetc(fp))!='\n' && c!=EOF && a<(COMMENT_SIZE-1))LS->seq_comment[current][a++]=c;
			LS->seq_comment[current][a]='\0';


			p=0;
			while ((c=fgetc(fp))!='>' && c!=EOF)
			        {
				  LS->seq[current][p++]=c;
				}

			LS->seq[current][p]='\0';
			LS->len[current]=strlen ( LS->seq[current]);

			current++;

			}

		else
		    c=fgetc ( fp);
		}


    vfclose (fp);


    vfree (sub);
    vfree (name);
    vfree (buffer);

    return LS;
}

void output_fasta_tree (char *fname, Alignment*A)
	{
	int a;
	FILE *fp;
	if ( !A || !A->nseq) return;

	if (A->Tree)return  output_fasta_tree (fname, A->Tree);
	
	fp=vfopen ( fname, "w");
	for ( a=0; a< A->nseq; a++)
	  {
	    A->seq_al[a]=substitute (A->seq_al[a], "\n", "");
	    fprintf ( fp, ">%s %s\n%s\n", A->name[a], A->seq_comment[a], A->seq_al[a]);
	  }
	vfclose (fp);
	}
void output_treelist (char *fname, Alignment*A)
	{
	int a;
	FILE *fp;
	if ( !A || !A->nseq) return;

	if (A->Tree)return  output_treelist (fname, A->Tree);
	
	fp=vfopen ( fname, "w");
	for ( a=0; a< A->nseq; a++)
	  {
	    A->seq_al[a]=substitute (A->seq_al[a], "\n", "");
	    fprintf ( fp, "%s\n",A->seq_al[a]);
	  }
	vfclose (fp);
	}  
  


NT_node check_tree (NT_node T);

NT_node main_read_tree (char *treefile)
{
  return file2tree(treefile);
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


NT_node find_node_in_tree (int *key, int nseq, NT_node T)
{
  if (!T || !key || !nseq)return NULL;
  else
    {
      int yes,a;
      NT_node C;
      int debug=0;
      if (debug)
	{
	  fprintf ( stderr, "\n\n");
	  for (a=0; a<nseq; a++)
	    fprintf (stderr, "%c", (key[a])?key[a]:'*');
	  fprintf ( stderr, "\n");
	  for (a=0; a<nseq; a++)
	    fprintf (stderr, "%d", T->lseq2[a]);
	  fprintf ( stderr, "\n");
	}
      for (yes=1,a=0; a<nseq; a++)
	{
	  if      ( key[a]&& !T->lseq2[a])
	    {
	      return NULL;
	    }
	  else if (!key[a]&&  T->lseq2[a])yes=0;
	}

      if (yes) return T;
      else if ((C=find_node_in_tree(key, nseq, T->right)))return C;
      else if ((C=find_node_in_tree(key, nseq, T->left )))return C;
      else
	{
	  HERE ("NOT FOUND AT ALL");
	  return NULL;
	}
    }
  return NULL;
}


NT_node recode_tree (NT_node T, Sequence *S)
{


  if (!T) return T;


  vfree (T->lseq);  T->lseq=(int*)vcalloc (S->nseq, sizeof (int));
  vfree (T->lseq2); T->lseq2=(int*)vcalloc (S->nseq, sizeof (int));
  vfree (T->idist); T->idist=(int*)vcalloc (S->nseq, sizeof (int));
  vfree (T->ldist); T->ldist=(int*)vcalloc (S->nseq, sizeof (int));
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
float *seq2dpa_weight (Sequence *S, char *mode)
{
  float *w;
  int a, b;
  if (!S) return NULL;
  w=(float*)vcalloc (S->nseq, sizeof (float));
  
  if (!mode || strm (mode, "longuest"))
    {
      for (a=0; a<S->nseq; a++)
	w[a]=strlen (S->seq[a]);
    }
  
  else if (strm (mode, "shortest"))
    {
      for (a=0; a<S->nseq; a++)
	w[a]=-strlen (S->seq[a]);
    }
  else if ( strm (mode, "name"))
    {
      for (a=0; a<S->nseq; a++)
	w[a]=-strlen (S->name[a]);
    }
  else if (mode[0]=='#')
    {
      char *seqf=vtmpnam (NULL);
      char *wf=vtmpnam (NULL);
      output_fasta_seqS (seqf,S);
      vfree (w);
      printf_system ("%s %s > %s", seqf, wf);
      return seq2dpa_weight (S, wf);
    }
  else if ( check_file_exists (mode))
    {
      char ***list=file2list (mode, " ");
      int n=0;
      while (list[n])
	{
	  char *name=list[n][1];
	  if ( name[0]=='>')
	    {
	      int i=name_is_in_list (name+1, S->name, S->nseq, MAXNAMES);
	      if (i!=-1)w[i]=atof (list[n][2]);
	    }
	  n++;
	}
    }
  else if ( strm (mode, "kmeans"))
    {
      int mdim=26*26*26;
      int nk=S->nseq/100;
      int *size;
      
      
      double **v=declare_double (S->nseq+1,mdim+4);
      
      size=(int*)vcalloc (nk+1, sizeof (char));
      for (a=0; a<S->nseq; a++)
	{
	  v[a]=seq2triaa(S->seq[a], v[a]);
	}
      
      v=vector2strip_vector (v, S->nseq, &mdim,0.80);
      
      km_kmeans (v, S->nseq, mdim,100,0.0001, NULL);
      
      for (a=0; a<S->nseq; a++)
	{
	  int c=(int)v[a][mdim+1];
	  size[c]++;
	}
      
      for (a=0; a<S->nseq; a++)
	{
	  int c=(int)v[a][mdim+1];
	  w[a]=(float)size[c];
	}
      free_double (v, -1);
      vfree (size);
    }
  else if (strm (mode, "swl") ||strm (mode, "iswl") || strm (mode, "swa"))
    {
      char **seql;
      int step,a,b, n;
      int **order;
      int invert;
      int mdim;
      double **v;
      float *aw;
      
      invert=(strm (mode, "iswl"))?-1:1;
      mdim=get_int_variable ("swlN");
      if (!mdim)mdim=100;
      if (mdim>S->nseq)mdim=S->nseq;
      aw=(float*)vcalloc (S->nseq, sizeof (float));
      seql=(char**)vcalloc (mdim, sizeof (char *));
      step=(S->nseq/mdim>0)?S->nseq/mdim:1;
      v=(double **)vcalloc (S->nseq+1, sizeof (double**));
      v[S->nseq]=(double*)vcalloc (mdim, sizeof (double));
      
      order=declare_int (S->nseq, 2);
      for (a=0; a<S->nseq; a++){order[a][0]=a; order[a][1]=strlen (S->seq[a]);}
      sort_int_inv (order, 2, 1, 0, S->nseq-1);
      

      for (n=0,a=0; a<S->nseq && n<mdim; a+=step, n++)
	{
	  seql[n]=S->seq[order[a][0]];
	}
      for (a=0; a<S->nseq; a++)
	{
	  output_completion (stderr,a,S->nseq, 100, "swl");
	  v[a]=seq2swr(S->seq[a], seql,mdim);
	  for (b=0; b<mdim; b++)
	    {
	      aw[a]+=(float)v[a][b];
	      v[S->nseq][b]+=v[a][b];
	    }
	}
      for (a=0; a<S->nseq; a++)
	{
	  
	  aw[a]/=(float)mdim;
	  HERE ("%s %.2f", S->name[a], aw[a]);
	}
      if (strm (mode, "swa"))
	{
	  w=aw;
	  vfree (w);
	}
      else
	{
	  vfree (aw);
	  for (b=0; b<mdim; b++)
	    {
	      v[S->nseq][b]/=(double)S->nseq;
	    }
	  
	  for (a=0; a<S->nseq; a++)
	    {
	      
	      double d=0;
	      for (b=0; b<mdim; b++)
		{
		  d+=(v[a][b]-v[S->nseq][b])*(v[a][b]-v[S->nseq][b]);
		}
	      w[a]=(float)invert*(float)sqrt(d);
	      
	    }
	}
      free_double (v, -1);
      free_int (order, -1);
      vfree (seql);
    }
  else if (strm (mode, "diaa") || strm (mode, "triaa")||strm (mode, "idiaa") || strm (mode, "itriaa"))
	{
      int mdim=26*26;
      double **v;
      int invert;
      
      if (mode[0]=='i'){mode++;invert=1;}
      else invert=-1;
      
      if (strm (mode, "triaa"))mdim*=26;
      
      v=declare_double (S->nseq+1,mdim+4);
      for (a=0; a<S->nseq; a++)
	{
	  if ( strm (mode, "diaa"))v[a]=seq2diaa(S->seq[a], v[a]);
	  else v[a]=seq2triaa(S->seq[a], v[a]);
	  
	  for (b=0; b<mdim;b++)v[S->nseq][b]+=v[a][b];
	}
      
      for (b=0; b<mdim; b++)
	{
	  v[S->nseq][b]/=(double)S->nseq;
	}
      
      for (a=0; a<S->nseq; a++)
	{
	  
	  double d=0;
	  for (b=0; b<mdim; b++)
	    {
	      
	      d+=(v[a][b]-v[S->nseq][b])*(v[a][b]-v[S->nseq][b]);
	    }
	  w[a]=(float)invert*(float)sqrt(d);
	  
	}
      free_double (v, -1);
    }
  else 
    {
      myexit (fprintf_error ( stderr, "\nERROR: unknown weight mode : %s [FATAL:%s]", mode, PROGRAM));
    }
  return w;
}

NT_node indextree2nametree (Sequence*S, NT_node root)
{
  NT_node current,pre;
  if ( !root)return root;
  if ( !S)
    {
      int n=tree2nseq (root);
      S=declare_sequence (10, 10, tree2nseq (root));
      S->nseq=0;
    }
  
  reset_node_count (root);
  current = root;
  root->visited=1;
  
  while(current != NULL)
    {                 
      if(!current->left)
	{
	  if (current && current->isseq && !current->visited)
	    {
	      int i=atoi (current->name)-1;
	      current->seq=i;
	      vfree (current->name);
	      current->name=(char*)vcalloc ( strlen (S->name[i])+1, sizeof (char));
	      sprintf (current->name, "%s", S->name[i]);
	      current->isseq=current->nseq=1;
	   } 
	  current->visited=1;
	  current = current->right;  
	  
	}    
      else
	{
	  pre = current->left;
	  
	  if (!pre->visited)
	    {
	      if (pre && !pre->visited && pre->isseq)
		{
		  int i=atoi (pre->name)-1;
		  pre->seq=i;
		  vfree (pre->name);
		  pre->name=(char*)vcalloc ( strlen (S->name[i])+1, sizeof (char));
		  sprintf (pre->name, "%s", S->name[i]);
		  pre->isseq=pre->nseq=1;
		}
	      pre->visited=1;
	    }
	  while(pre->right && pre->right != current)
	    {
	      pre = pre->right;
	      if (pre)
		{
		  if (pre && !pre->visited && pre->isseq)
		    {
		      int i=atoi (pre->name)-1;
		      pre->seq=i;
		      vfree (pre->name);
		      pre->name=(char*)vcalloc ( strlen (S->name[i])+1, sizeof (char));
		      sprintf (pre->name, "%s", S->name[i]);
		      pre->isseq=pre->nseq=1;
		      
		    }
		  pre->visited=1;
		}
	    }
	  
	  
	  if(!pre->right )
	    {
	      pre->right = current;
	      current = current->left;
	    }
	  
	  else 
	    {
	      pre->right = NULL;
	      current = current->right;      
	    } 
	} 
    } 
  return root;
}

NT_node node2master(NT_node T, Sequence *S, float *w)
{
  //select the one with a PDB template OR select the longuest
  if (!T)return T;
  else if (!T->left && !T->right)
    {
      //int i=name_is_in_list (T->name, S->name, S->nseq, MAXNAMES);
      int i=name_is_in_hlist (T->name, S->name, S->nseq);
      
      if (i==-1)HERE ("Could Not Find %s in Tree", T->name);
            
      T->seq=i;
      T->score=w[i];
      T->isseq=1;
      T->nseq=1;
      
    }
  else if (!T->left || !T->right)
    {
      HERE ("incorrectly parsed tree Root: %d", T->parent?1:0);
      exit (0);
    }
  else 
    {
      NT_node L=node2master(T->left,S,w);
      NT_node R=node2master(T->right,S,w);
      
      T->isseq=0;
      T->leaf=T->nseq=L->nseq+R->nseq;
      if      (R->seq==-1 && L->seq==-1)T->seq=-1;
      else if (L->seq==-1)L=R;
      else
	{
	  int pl=(seq2P_template_file(S, L->seq))?1:0;
	  int pr=(seq2P_template_file(S, L->seq))?1:0;
	  
	  int ll=strlen (S->seq[L->seq]);
	  int lr=strlen (S->seq[R->seq]);
	  float wl=w[L->seq];
	  float wr=w[R->seq];
	  
	  if ((!pl &&!pr) || (pl && pr))
	    {
	      if (abs((wr-wl))<0.0001)
		{
		  if (ll<lr)R=L;
		}
	      else if (wl>wr)R=L;
	    }
	  else if (pl) 
	    R=L;
	}
      T->score=R->score;
      T->seq=R->seq;
      
      sprintf (T->name,"%s",R->name);
    }
  
  return T;
}

NT_node * sort_nodelist4dpa (NT_node *L, int n);
NT_node * sort_nodelist4dpa (NT_node *L, int n)
{
  int a;
  float **lu=declare_float ( n, 2);
  NT_node *NL=(NT_node*)vcalloc (n, sizeof (NT_node));
  
  for (a=0; a<n; a++)
    {
      lu[a][0]=a;
      lu[a][1]=(L[a])->score*-1;
    }
  sort_float (lu, 2, 1, 0, n-1);
  for (a=0; a<n; a++)
    {
      NL[a]=L[(int)lu[a][0]];
    }
  vfree (L);
  free_float (lu, -1);
  return NL;
}




  
int tree2clusters   (NT_node T, int *nc,int **cl, double **dist, double Thr, int min)
{
  if (!T) return nc[0];
  int a,b;
  double td=0;
  double n=0;


  
  
  for (a=0; a<T->nseq-1; a++)
    {
      for (b=a+1; b<T->nseq; b++)
	{
	  int s1=T->lseq[a];
	  int s2=T->lseq[b];
	  td+=dist [s1][s2];
	  //fprintf ( stdout, "\n\t %d %d ---> %f", s1, s2,dist[s1][s2]);
	  n++;
	}
    }
  
  if (n>0.00001)td/=n;

  //fprintf (stdout, "\nTest Cluster => %f\n", td);
  //for (a=0; a<T->nseq; a++)
  //  fprintf ( stdout, "%d ", T->lseq[a]);

  
  if (td<Thr && T->nseq>=min)
    {
      for (a=0; a<T->nseq; a++)
	cl[nc[0]][a+1]=T->lseq[a];
      cl[nc[0]][0]=T->nseq;
      nc[0]++;
    }
  else
    {
      tree2clusters (T->left, nc,cl,dist, Thr,min);
      tree2clusters (T->right, nc,cl,dist, Thr, min);
    }
  return nc[0];
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

NT_node display_splits (NT_node T,Sequence *S, FILE *fp, char *name)
{
  int a;
  if (!T) return T;

  if (!S)S=tree2seq (T,NULL);

  display_splits (T->right,S, fp, name);
  display_splits (T->left, S, fp, name);



  if (!T->right);
  else if (!T->left);
  else if (!T->parent );
  
  //else if (T->parent && !(T->parent)->parent);
  //Removed because it prevents some splits to be reported replaced with !T-Parent test.
  else
    {
      int t=0;
      for (a=0; a< S->nseq; a++)
	{
	  fprintf (fp, "%d", T->lseq2[a]);
	  t+=T->lseq2[a];
	}

      fprintf ( fp, " %d", MIN(t,((S->nseq)-t)));
      if (name)fprintf (fp, " %s", name);
      fprintf (fp, "\n");
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

NT_node file2tree (char *fname)
{
  FILE *fp, *fp1, *fp2;
  
  NT_node R,T,N;
  int c, lastc; 
  char *tmp=vtmpnam (NULL);
  char **nl=NULL;
  int  nn=0;
  char *dup=NULL;
  
  if (!fname || !check_file_exists (fname))return NULL;
  fp=vfopen (remove_charset_from_file (fname, " \t\n\r"), "r");
  
  if (!fp)return NULL;
  R=T=new_declare_tree_node ();
  
  while ((c=fgetc(fp))!=';' && c!=EOF)
    {
      if (c=='(')
	{
	  N=new_declare_tree_node ();
	  N->parent=T;
	  if     (!T->right)T->right=N;
	  else if(!T->left)T->left=N;
	  else
	    {
	      N->right=T->right; (T->right)->parent=N;
	      N->left=T->left; (T->left)->parent=N;
	      T->right=N;
	      N=T->left=new_declare_tree_node ();
	      N->parent=T;
	    }
	  
	  T=N;
	  
	  lastc=0;
	}
      else if (c==')')
	{
	  T=T->parent;
	  scan_name_and_dist (fp, T->name, &T->dist);
	  if (T->name && T->name [0])T->bootstrap=atof (T->name);
	  
	  lastc=0;
	}
      else if (c==',')
	{
	  T=T->parent;
	  lastc++;
	}
      else
	{
	  N=new_declare_tree_node ();
	  N->parent=T;
	  if      (T && !T->right)T->right=N;
	  else if (T && !T->left)T->left=N;
	  else 
	    {
	      N->right=T->right; (T->right)->parent=N;
	      N->left=T->left; (T->left)->parent=N;
	      T->right=N;
	      N=T->left=new_declare_tree_node ();
	      N->parent=T;
	    }
	  
	  T=N;
	  	  
	  ungetc (c, fp);
	  scan_name_and_dist (fp, T->name, &T->dist);
	  nl=(char**)vrealloc (nl, (nn+1)*sizeof (char*));
	  
	  nl[nn++]=T->name;
	  T->leaf=1;
	  T->isseq=1;
	  lastc=0;
	}
    }
  T=T->parent;
  vfclose (fp);
  
  if (R!=T)
    {
      HERE ("Parsing Error\n");
    }
  
  if ( (dup=check_hlist_for_dup(nl,nn)))
	{
	  fprintf ( stderr, "ERROR -- Duplicated Sequences %s", dup);
	  myexit(fprintf_error (stderr,"ERROR - file2tree: %s is duplicated in %s  ", fname));
	}
 
  else
    {
      vfree (nl);
    }
  
 if (!T->right && T->left)T=T->left;
 else if (T->right && !T->left)T=T->right;
 T->parent=NULL;
      
  return T;
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
  while (isdigit((c=fgetc(fp))) || c=='.' || c=='-' || c=='e')
    {
      number[a++]=c;
    }

  ungetc (c, fp);
  number[a]='\0';

  dist[0]=atof (number);

  return 2;
}


NT_node new_declare_tree_node ()
{
  //does not declare seqlist --> important when making very large trees
  NT_node p=NULL;
  static int node_index;
	
  p= (NT_node)vcalloc (1, sizeof ( Treenode));
  if (!p)HERE ("Could Not Allocate Node");
	
	
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
void display_tree_lseq2 (NT_node T, int n)
{
  if (!T) return;
  else 
    {
      int a;
      for (a=0; a<n; a++)
	fprintf ( stderr, "%d", T->lseq2[a]);
      fprintf (stderr, "\n");
      display_tree_lseq2(T->left,n);
      display_tree_lseq2(T->right,n);
    }
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
  dup=(int*)vcalloc ( S->nseq, sizeof (int));

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

float newick2avg_bs  (char *string)
{
  return tree2avg_bs (newick_string2tree(string));
}
float tree2avg_bs ( NT_node T)
{
  float tot;
  int n;
  if (!T) return 0;
  tot=tree2tot_dist (T, BOOTSTRAP);
  n=tree2n_branches (T, BOOTSTRAP);
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
      if ( mode == BOOTSTRAP && T->bootstrap!=0)t+=T->bootstrap;
      else t+=T->dist;
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
NT_node * tree2seqnode_list (NT_node p, NT_node *L)
{
  int n=0;
  int nn=tree2nnode(p);
  reset_node_count (p);
  
  if (!L) {L=(NT_node*)vcalloc (nn+1, sizeof (NT_node));}
  while (p)
    {
      if (!p->visited && p->isseq)L[n++]=p;
      
      int x=++(p->visited);
      
      if (!p->isseq)
	{
	  if (x==1)
	    {
	      p=p->right;
	    }
	  else if (x==2)
	    {
	      p=p->left;
	    }
	  else if (x==3 && !p->parent)
	    {
	      p=p->parent;
	    }
	}
      
      if (x>=3 || p->isseq)
	{
	  if (p && p->isseq && !p->visited ){L[n++]=p, p->visited=1;}
	  
	  if (x>3)p=NULL;
	  else if (p) p=p->parent;
	}
    }
  return L;
}

NT_node   duplicate_tree (NT_node T)
{
  NT_node *L1=tree2node_list(T,NULL);
  NT_node *L2=(NT_node*)vcalloc (arrlen(L1)+2, sizeof (NT_node));
  NT_node ROOT;
  int a;
  
  a=0;
  while (L1[a])
    {
      L1[a]->order=a;
      L2[a]=new_declare_tree_node();
      a++;
    }
  a=0;
  while (L1[a])
    {
      NT_node R=L1[a];
      NT_node N=L2[a];
      int p, r, l;
      
      p=(R->parent)?(R->parent)->order:-1;
      l=(R->left )?(R->left )->order:-1;
      r=(R->right)?(R->right)->order:-1;
      
      if (R->parent)N->parent=L2[(R->parent)->order];
      if (R->left  )N->left  =L2[(R->left  )->order];
      if (R->right )N->right =L2[(R->right )->order];
      
      if (!N->parent)ROOT=N;
      if (R->isseq)
	{
	  N->isseq=1;
	  N->name=csprintf (N->name, "%s", R->name);
	}
      a++;
    }
  vfree (L1); vfree (L2);
  return ROOT;
}
      
NT_node * tree2node_list (NT_node p, NT_node *L)
{
  int n=0;
  int nn=tree2nnode(p);
  reset_node_count (p);
  
  if (!L) {L=(NT_node*)vcalloc (nn+2, sizeof (NT_node));}
  while (p)
    {
      if (!p->visited)L[n++]=p;
      
      int x=++(p->visited);
      
      if (!p->isseq)
	{
	  if (x==1)
	    {
	      p=p->right;
	    }
	  else if (x==2)
	    {
	      p=p->left;
	    }
	  else if (x==3 && !p->parent)
	    {
	      p=p->parent;
	    }
	}
      
      if (x>=3 || p->isseq)
	{
	  if (p && p->isseq && !p->visited ){L[n++]=p, p->visited=1;}
	  
	  if (x>3)p=NULL;
	  else if (p) p=p->parent;
	}
    }
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

NT_node newick_string2tree (char *string)
{
  char *tmp=vtmpnam (NULL);
  
  if (!string) return NULL;
  printf_file (tmp,"w", "%s", string);
  
  return main_read_tree (tmp);
}
  
NT_node * read_tree_list (Sequence *S)
{
  NT_node *T;
  int a;

  T=(NT_node*)vcalloc ( S->nseq+1, sizeof (NT_node));

  for ( a=0; a<S->nseq; a++)
    {
      char *fname;
      if (S->seq && S->seq[a] && strlen (S->seq[a])<2)
	fname=S->name[a];
      else
	string2file ((fname=vtmpnam(NULL)), "w", S->seq[a]);

      T[a]=main_read_tree (fname);
      T[a]->file=(char*)vcalloc (strlen (S->name[a])+1, sizeof (char));
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

  name=(char*)vcalloc (1000, sizeof (char));
  fname=(char*)vcalloc (1000, sizeof (char));
  group=(char*)vcalloc (TS->nseq*10, sizeof (char));
  ref_group=(char*)vcalloc (TS->nseq*10, sizeof (char));
  list=(char*)vcalloc (100*S->nseq, sizeof (char));
  split_file=vtmpnam (NULL);
  sorted_split_file =vtmpnam (NULL);

  n=S->nseq;
  used=(int*)vcalloc (n, sizeof (int));

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



  gs=(int*)vcalloc (2, sizeof (int));
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
      rlist=(int***)declare_arrayN(3, sizeof (int), 2,LIST->nseq*S->nseq, S->nseq+1);
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
      rlist=(int***)vcalloc ( 2, sizeof (int**));
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
  static char *def_param;
  char *cache=NULL;
  //+count_splits _NB_x_FILTER_<file>
  //_<file is a fasta file containing the list of species to keep>
  if (!def_param)def_param=(char*)vcalloc ( 10, sizeof (char));



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
      cache=(char*)vcalloc (S->nseq, sizeof (int));
      for ( a=0; a<F->nseq; a++)
	{
	  if ( (i=name_is_in_list (F->name[a], S->name, S->nseq, 100))!=-1)
	    cache[i]=1;
	}
      free_sequence (F, -1);
    }

  main_buf=(char*)vcalloc ( S->nseq*(STRING+1), sizeof(int));

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
      buf=(char*)vcalloc (measure_longest_line_in_file (out)+1, sizeof (char));
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
  SL=(Split**)vcalloc ( n1+1, sizeof (Split*));

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
  S=(Split*)vcalloc (1, sizeof (Split));
  S->split=(char*)vcalloc ( nseq+1, sizeof (char));
  return S;
}
int treelist2splits( Sequence *S, Sequence *TS)
{
  NT_node *T;
  int n=0,nseq, a, c;

  int *used;

  char *split_file, *sorted_split_file;
  char *buf=NULL, *ref_spl=NULL;
  char *spl;
  char *fname;
  char **wl;
  FILE *fp;
  char file_list[100000];
  
  split_file=vtmpnam (NULL);
  sorted_split_file =vtmpnam (NULL);

  n=S->nseq;
  used=(int*)vcalloc (n, sizeof (int));

  T=read_tree_list (S);
  if (!TS)TS=tree2seq(T[0], NULL);
  nseq=TS->nseq;
  fp=vfopen (split_file, "w");


  for ( a=0; a< S->nseq; a++)
    {
      
      T[a]=prune_tree  (T[a], TS);
      T[a]=recode_tree (T[a], TS);
      display_splits (T[a], TS,fp, S->name[a]);
    }
  
  vfclose (fp);
  printf_system ("cp %s split_file::IGNORE_FAILURE::", split_file);
    
  printf_system ( "cat %s | grep 1| sort > %s::IGNORE_FAILURE::", split_file, sorted_split_file);

  fp=vfopen (sorted_split_file, "r");
  
  for ( a=0; a<TS->nseq; a++)fprintf ( stdout, "SEQ_INDEX %d %s\n", a+1, TS->name[a]);
  fprintf(stdout, "#Legend: <1: minimum split size>\t<2: number of occurences>\t<3: coded split>\t<4: File list <file1>####<file2>>\t<5: Group 1 <leaf1>####<leaf2>>\t<6: Group 2 <leaf1>####<leaf2>>\n"); 
  while ( (c=fgetc (fp))!=EOF)
    {
      ungetc (c, fp);
      buf=vfgets (buf, fp);
      buf [strlen(buf)-1]='\0';
      
      wl=string2list (buf);
      spl=wl[1];
      fname=wl[3];
            
      if (!ref_spl)
	{
	  ref_spl=(char*)vcalloc (strlen (spl)+1, sizeof (char));
	  sprintf ( ref_spl, "%s", spl);
	  file_list[0]='\0';
	  n=1;
	}
      else if ( !strm (spl, ref_spl))
	{
	  int i;
	  int n0=0, n1=0;
	  
	  for (i=0; i<nseq; i++)
	    {
	      n0+=(ref_spl[i]=='0')?1:0;
	      n1+=(ref_spl[i]=='1')?1:0;
	    }
	  
	  if (n>1 && strstr (file_list, fname))n--;
	  else
	    {
	      if (n>1)strcat (file_list, "####");
	      strcat (file_list, fname);
	    }
	  
	  fprintf ( stdout, "%d\t%d\t%s\t%s\t",((n0>n1)?n1:n0),n,ref_spl, file_list);
	  
	  for (i=0,a=0; a<nseq; a++)
	    if (ref_spl[a]=='1')
	      {
		if (i==1)fprintf(stdout, "####");
		fprintf (stdout, "%s", TS->name[a]);
		i=1;
	      }
	  fprintf ( stdout, ",GROUP2::,");
	  for (i=0,a=0; a<nseq; a++)
	    if (ref_spl[a]=='0')
	      {
		if (i==1) fprintf ( stdout, "####");
		fprintf (stdout, "%s", TS->name[a]);
		i=1;
	      }
	  fprintf (stdout, "\n");
	  file_list[0]='\0';
	  sprintf ( ref_spl, "%s", spl);
	  n=1;
	}
      else
	{
	  if (n>1)strcat (file_list, "####");
	  strcat (file_list, fname);
	  n++;
	}
      free_char (wl, -1);
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
  used=(int*)vcalloc (n, sizeof (int));

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
	  ref_buf=(char*)vcalloc (strlen (buf)+1, sizeof (char));
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
  T=(NT_node*)vrealloc (T, (S->nseq+1)*sizeof (NT_node));
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
	  S->seq[b]=(char*)vrealloc (S->seq[b], (strlen (s)+1)*sizeof (char));
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
  dist=(int****)vcalloc ( S->nseq, sizeof (int****));
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
      d=(int***)declare_arrayN (3, sizeof (float),2, S->nseq, S->nseq);
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

   used=(int*)vcalloc (S->nseq, sizeof (int));
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

  TS1=(Tree_sim*)vcalloc (1, sizeof (Tree_sim));
  TS2=(Tree_sim*)vcalloc (1, sizeof (Tree_sim));


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

  L2=(NT_node*)vcalloc ( n+1, sizeof (NT_node));
  for (a=0; a<n; a++)
    L2[a]=L[score[a][0]];

  BT=treelist2bootstrap (L2, NULL);

  vfree (L2);
  if (file)free_treelist(L);
  return BT;
}
Alignment *treelist2cons (Alignment *T)
{
  char *infile=vtmpnam (NULL);
  char *outfile=vtmpnam (NULL);
  FILE *fp;
  int a=0;

  if (!T) return T;
  else if (T->Tree)return treelist2cons (T->Tree);
  
  fp=vfopen (infile, "w");
  for (a=0; a<T->nseq; a++)
    {
      fprintf (fp, "%s\n", T->seq_al[a]);
    }
  vfclose (fp);
    
  printf_system( "msa2bootstrap.pl -i %s -o %s -input tree >/dev/null 2>/dev/null", infile, outfile);
  
  vfree (T->seq_al[0]);
  T->seq_al[0]=file2string(outfile);
  
  sprintf (T->name[0], "MajorityConsensusTree");
  sprintf (T->seq_comment[0], "Raplicates: %d Source: consense", T->nseq);
  
  T->nseq=1;
  return T;
}

float group2tree_score (Sequence *G, char *nwtree)
{
  //returns the fraction of G sequences in the node containing at least all the G sequences
  float bs=0;
  NT_node T, *L, *iL;
  Sequence *S;
  int a, b;
  int *l1, *bl2;
  
  if (!nwtree || !nwtree[0]) return -1;
  
  bl2=NULL;
  T=newick_string2tree(nwtree);
  S=tree2seq(T, NULL);
  l1=(int*)vcalloc (S->nseq, sizeof (int));
  
  for (a=0; a<G->nseq; a++)
    {
      if ((b=name_is_in_list (G->name[a], S->name, S->nseq, 100))!=-1)
	l1[b]=1;
    }
  iL=L=tree2node_list (T, NULL);

  while (L && L[0])
    {
      if (L[0]->lseq2)
	{
	  float c, t,cs;
	  int *l2=(L[0])->lseq2;
	  c=t=0;
	  for (a=0; a<S->nseq; a++)
	    {
	      int x=l1[a];
	      int y=l2[a];
	      if (l1[a]+l2[a]==2)c++;
	      else if (l2[a])t++;
	      else if (l1[a]){a=S->nseq; c=0;}
	    }
	  
	  cs=(c>0)?c/(c+t):0;
	  if (cs>bs){bs=cs;bl2=l2;}
	  
	  for (c=0,t=0,a=0; a<S->nseq; a++)l2[a]=1-l2[a];
	  for (c=0,t=0,a=0; a<S->nseq; a++)
	    {
	      if (l1[a]+l2[a]==2)c++;
	      else if (l2[a])t++;
	      else if (l1[a]){a=S->nseq; c=0;}
	    }
	  cs=(c>0)?c/(c+t):0;
	  if (cs>bs){bs=cs;bl2=l2;}
	}
      L++;
    }
  
  
  free_tree (T);
  free_sequence (S, S->nseq);
  vfree (iL);
  vfree (l1);
  return bs;
}

  
Alignment *treelist2node_support_best (Alignment *T)
{
  NT_node *RT, TT;
  Sequence *S=NULL;
  int a, b, best_tree;
  float s=0, best_s=0;
  char *p;
  
  if (T->Tree)return treelist2node_support_best (T->Tree);
  
  RT=(NT_node*)vcalloc (T->nseq, sizeof (NT_node));
  for (a=0; a<T->nseq; a++)
    {
      output_completion (stderr,a,T->nseq,1, "Input Replicate Trees");
      RT[a]=newick_string2tree(T->seq_al[a]);
      if (!S)S=tree2seq(RT[a], NULL);
      RT[a]=prune_tree  (RT[a],S);
      RT[a]=recode_tree (RT[a],S);
    }
  
  for (best_s=best_tree=0,a=0;a<T->nseq; a++)
    {
      output_completion (stderr,a,T->nseq,1, "Scan for the best Supported Tree");
      for (s=0,b=0; b<T->nseq; b++)
	{
	  s+=update_node_support (RT[a],RT[b],S->nseq);
	}
      if (s>best_s){best_s=s; best_tree=a;}
    }
  
  s=(best_s/(float)(T->nseq-1))*100;

  p=tree2string (RT[best_tree]);
  T->seq_al[0]=(char*)vrealloc (T->seq_al[0],(strlen (p)+1)*sizeof (char));
  sprintf (T->seq_al[0], "%s",p);
  
  sprintf (T->name[0], "OriginalDistanceTree");
  sprintf (T->seq_comment[0], "Replicates: %d AverageNodeSupport: %.2f %% Source: seq_reformat", T->nseq, s);
  T->nseq=1;
  T->score=T->score_aln=s*100;

  for (a=0; a<T->nseq; a++)free_tree (RT[a]);
  vfree (RT);
  return T;
}
Alignment *treelist2bs_compare (Alignment *T1, Alignment *T2)
{
  int a;
  
  if (!T1) return T1;
  if (!T2) T2=T1;
  fprintf ( stdout, "#TreeIndex;Node_Count;Node Bootsrap;NormalizedNodeDepth\n");
  for (a=0; a<T1->nseq; a++)
    {
      float bs=0;

      char *t=tree2node_support_simple (T1->seq_al[a], T2, &bs);

      compare_bs (T1->name[a], t, T1->seq_al[a]);
      vfree (t);
    }
  return T1;
}
void compare_bs (char *string, char *T1, char *T2)
{
  NT_node t1, t2;
  Sequence *S;
  t1=newick_string2tree(T1);
  t2=newick_string2tree(T2);
  S=tree2seq(t1, NULL);
  t1=recode_tree (t1,S);
  t2=recode_tree (t2,S);
  
  display_compare_bs (string,t1, t2, S->nseq);
  free_tree (t1);
  free_tree(t2);
  free_sequence (S, S->nseq);
  
}
void display_compare_bs (char *string, NT_node t1, NT_node t2, int nseq)
{
  if (!t1) return;
  
  else
    {
      int a;
      if (t1->left)
	{
	  float depth=0;
	  for (a=0; a<nseq; a++) 
	    {
	      //fprintf ( stdout, "%d", t1->lseq2[a]);
	      depth+=t1->lseq2[a];
	    }
	  
	  depth=MIN(depth, (nseq-depth))/(nseq/2);
	  if (depth<1)fprintf ( stdout, "%s;%.4f;%.4f;%.2f\n", string?string:"",t1->bootstrap, t2->bootstrap, depth);
	}
      display_compare_bs (string, t1->left ,t2->left,nseq );
      display_compare_bs (string, t1->right,t2->right,nseq);
    }
  return;
}
Alignment *treelist2node_support (Alignment *T)
{
  if (!T) return T;
  else if (T->Tree)return treelist2node_support (T->Tree);
  
  return tree2node_support (T->seq_al[0], T);
}
NT_node reset_bs2 (NT_node T);
NT_node reset_bs2 (NT_node T)
{
  if (!T) return T;
  T->bootstrap=0;
  reset_bs2 (T->left);
  reset_bs2 (T->right);
  T->bootstrap=0;
  return T;
}
char*tree2node_support_simple (char *newick_tree, Alignment *T,float *bs)
{
  NT_node RT, TT;
  Sequence *S;
  int a;
  float s=0;
  char *p;
  RT=newick_string2tree(newick_tree);
  
  S=tree2seq(RT, NULL);
  RT=recode_tree (RT,S);
  RT=reset_bs2 (RT);
  
  for (a=1; a<T->nseq; a++)
    {
      TT=newick_string2tree(T->seq_al[a]);
      TT=prune_tree (TT,S);
      TT=recode_tree(TT,S);
      s+=update_node_support (RT, TT,S->nseq);
      free_tree (TT);
    }
  
  
  
  //make sure new tree is not longuer than previous
  if (bs)bs[0]=s;
  p=tree2string (RT);
  free_sequence (S, S->nseq);
  free_tree (RT);
  return p;
}
Alignment *tree2node_support (char *newick_tree, Alignment *T)
{
  char *p;
  float s=0;
  

  p=tree2node_support_simple(newick_tree, T, &s);
  
  T->seq_al[0]=(char*)vrealloc (T->seq_al[0],(strlen (p)+1)*sizeof (char));
  sprintf (T->seq_al[0], "%s",p);
  sprintf (T->name[0], "OriginalDistanceTree");
  sprintf (T->seq_comment[0], "Replicates: %d AverageNodeSupport: %.2f %% Source: seq_reformat", T->nseq, s);
  s=(s/(float)(T->nseq-1))*100;
  T->score=T->score_aln=(float)s*100;
  //T->nseq=1;
  vfree (p);
  return T;
}
float update_node_support (NT_node T1, NT_node T2, int nseq)
{
  //Trees must have been -prune_tree- and -recode_tree- with the same referecne sequence
  //updates node boostrap values in T1
  //Returns the number of nodes in T2 supporting nodes in T1
  
  float sc=0;
  float tot=0;
  NT_node *L1,*iL1,*L2,*iL2;
  int *l1seq, *l2seq, *rl1seq;
  int a;
  if (!T1 || !T2 || !nseq) return 0;
  

 
  L1=iL1=tree2node_list (T1, NULL);
  L2=iL2=tree2node_list (T2, NULL);
  rl1seq=(int*)vcalloc (nseq, sizeof (int));
  
  while (L1 && L1[0])
    {
      NT_node*L2=iL2;
      while (L2 && L2[0])
	{
	  
	  l1seq=(L1[0])->lseq2;
	  l2seq=(L2[0])->lseq2;
	  for (a=0; a<nseq; a++)rl1seq[a]=1-l1seq[a];
	  
	  if (memcmp (rl1seq, l2seq, nseq*sizeof(int)) && memcmp (l1seq, l2seq, sizeof (int)*nseq))L2++;
	  else
	    {
	      sc++;
	      (L1[0])->bootstrap++;
	      L2=NULL;
	    }
	}
      tot++;
      L1++;
    }
  
  vfree (iL1);
  vfree (iL2);
  vfree (rl1seq);
  return (tot==0)?0:(sc/tot);
}
//Quantify the support of every node in a tree against collections of trees
//the collections of trees are provided via B
int **treelist2ns (NT_node T,Sequence *B,char *ref)
  {
    NT_node R, TT, *NL;
    Sequence *S;
    int n, r, a,c,b, t;
    int **support;
    int *max;
    int *tot;
    int *nr;
    int *bsr;
    
    S=tree2seq(T, NULL);
    T=recode_tree(T,S);
    NL=tree2node_list (T, NULL);
    n=tree2nnode(T);
    for (a=0; a<n; a++)(NL[a])->bootstrap=0;
    
    if (ref)
      {
	R=main_read_tree (ref);
	R=recode_tree(R,S);
	update_node_support (T, R, S->nseq);
      }
    
    support=declare_int (B->nseq+1, n);
    max=(int*)vcalloc (B->nseq+1, sizeof (int));
    tot=(int*)vcalloc (B->nseq+1, sizeof (int));
    nr=(int*)vcalloc (B->nseq+1, sizeof (int));
    
    for (a=0; a<n; a++)
      {
	support[0][a]=(int)(NL[a])->bootstrap;
	(NL[a])->bootstrap=0;
      }
    
    bsr=(int*)vcalloc (B->nseq+1, sizeof (int));
    bsr[0]=1;
    for (r=1; r<=B->nseq; r++)
      {
	if (check_file_exists (B->name[r-1]))bsr[r]=1;
	else bsr[r]=0;
      }
    for (r=1; r<=B->nseq; r++)
      {
	int nt;
	char ** list;
	

	HERE ("--- Process %s", B->name[r-1]);
	
	if (!check_file_exists (B->name[r-1]))
	  {
	    printf_exit (EXIT_FAILURE,stderr, "%s is not a valid file. treelist2ns takes as input a fasta-like list of files containing trees", B->name[r-1]);
	  }
	    
	list=file2lines (B->name[r-1]);
	if (list==NULL)myexit (fprintf_error (stderr,"File %s is empty",B->name[r-1]));
	
	nt=atoi(list[0]);
	nr[r]=nt-1;
	for (t=1; t<nt; t++)
	  {
	    float x;
	    TT=newick_string2tree(list[t]);
	    TT=prune_tree (TT,S);
	    TT=recode_tree(TT,S);
	    update_node_support (T, TT,S->nseq);
	    free_tree (TT);
	  } 
	for (a=0; a<n; a++)
	  {
	    support[r][a]=(int)(NL[a])->bootstrap;
	    (NL[a])->bootstrap=0;
	  }
	free_char (list, -1);
      }
    
    fprintf (stdout, "#NodeSupport_01\n#TaxonList:");
    for (a=0; a<S->nseq; a++)
      fprintf (stdout, "%s ", S->name[a]);
    fprintf (stdout, "\n");
    if (ref)fprintf (stdout, "#REF: %s\n", ref);
    else fprintf (stdout, "#REF: none\n");
    
    fprintf (stdout, "#N_REPLICATE: %d\n", B->nseq);
    for (a=0; a<B->nseq; a++)
      fprintf (stdout, "#REPLICATE: %3d %s %4d\n", a+1, B->name[a], nr[a+1]);
    
    fprintf (stdout, "#NODE: <list> <depth> <correct> <support from replicate X..>\n");
	     
    for (a=0; a< n; a++)
      {
	if ((NL[a])->nseq>1 && (NL[a]->nseq<(S->nseq-1)))
	  {
	    int depth=0;
	    for (r=0; r<S->nseq; r++){fprintf ( stdout , "%d", (NL[a]->lseq2[r]));depth+=NL[a]->lseq2[r];}
	    fprintf (stdout, " ");
	    depth=MIN(depth, (S->nseq-depth));
	    fprintf (stdout, "%2d ", depth);
	    for (r=0; r<=B->nseq; r++)
	      {
		if (r>0)
		  {
		    tot[r]+=support[r][a];
		    max[r]++;
		    fprintf (stdout, "%6.4f ", (float)support[r][a]/(float)nr[r]);
		  }
		else
		  fprintf (stdout, "%3d ",support[r][a]);
	      }
	    fprintf (stdout, "\n");
	  }
      }
    fprintf (stdout, "#AVG: ");
    for (r=1; r<=B->nseq; r++)fprintf (stdout, "%6.4f ", ((float)tot[r]/(float)max[r])/(float)nr[r]);
    fprintf (stdout, "\n");
    
    if (ref)
      {
	float yes=0;
	float no=0;
	float totn=0;
	
	for (a=0; a<n; a++)
	  if ((NL[a])->nseq>1 && (NL[a]->nseq<(S->nseq-1)))
	    {
	      if (support[0][a])yes++;
	      else no++;
	      totn++;
	    }
	//B->nseq=1;
	for (c=0; c<=1; c++)
	  {
	    fprintf (stdout, "#AVG_NODE_SUPPORT: Correct: %s %6.4f ", (c==0)?"N":"Y", ((c==1)?yes:no)/totn);
	    for (r=1; r<=B->nseq; r++)
	      {
		float t=0;
		float m=0;
		for (a=0; a<n; a++)
		  if ((NL[a])->nseq>1 && (NL[a]->nseq<(S->nseq-1)) && support[0][a]==c)
		    {
		      //fprintf  (stdout, "\n---- %d %d - %d\n", c, (int)support[r][a], nr[r]);
		      t+=support[r][a];
		      m+=nr[r];
		    }
		if (m)fprintf (stdout, "%6.4f ", t/m);
		else fprintf  (stdout, "------ ");
	      }
	    fprintf (stdout, "\n");
	  }
      }
    fprintf (stdout, "#END\n");
    exit (0);
    return support;
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

  name=(char**)vcalloc (1, sizeof (char*));
  fp=vfopen ((fname=vtmpnam (NULL)), "w");

  T=read_tree_list (S);
  for (n=0,a=0; a< S->nseq; a++)
    {
      TS=tree2seq(T[a], NULL);
      for (b=0; b<TS->nseq; b++)
	{
	  if ( (i=name_is_in_list (TS->name[b], name, n, 100))==-1)
	    {
	      name[n]=(char*)vcalloc (100, sizeof (int));
	      sprintf ( name[n], "%s", TS->name[b]);
	      n++;
	      name=(char**)vrealloc (name, (n+1)*sizeof (char*));
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

NT_node tree2nni (NT_node S, NT_node T)
{

  if (!S)return NULL;
  if (!T)T=S;
  
  
  if (has_nni(S))
    {
      nni (S,0);
      print_newick_tree (T, "stdout");
      nni (S,0);
      
      nni (S,1);
      print_newick_tree (T, "stdout");
      nni (S,1);
    }
  
  tree2nni(S->left,T);
  tree2nni(S->right,T);
  return S;
}
int has_nni (NT_node N)
{
  if (!N)return 0;
  if (!N->left)return 0;
  if (!N->right)return 0;
  if (!(N->left)->left)return 0;
  if (!(N->left)->right)return 0;
  if (!(N->right)->left)return 0;
  if (!(N->right)->right)return 0;
  return 1;
}
NT_node nni (NT_node S, int n)
{
  NT_node I;

  I=(S->right)->right;

  if      (n==0){(S->right)->right=(S->left)->right;(S->left)->right=I;}
  else if (n==1){(S->right)->right=(S->left)->left ;(S->left)->left =I;}

  return S;
}
////////////////////////////////////// Paralel DPA
ALN_node declare_aln_node(int mode);
int ktree2aln_bucketsF(KT_node K,char *fname);


char   *kmsa2msa (KT_node K,Sequence *S, ALNcol***S2,ALNcol*PG);
ALNcol * msa2graph (Alignment *A, Sequence *S, ALNcol***S2, ALNcol*msa,int seq);



NT_node kmsa2dnd(Sequence *S,KT_node *KL, int n);

NT_node tree2dnd4dpa (NT_node T, Sequence *S, int N, char *method)
{
  int n=0;
  char *outname;
  int cn=0;
  KT_node K =tree2ktree  (T,T, S, N);
  KT_node*KL=(KT_node*)vcalloc (K->tot, sizeof (KT_node));
  
  n=ktree2klist(K,KL,&n);
  kseq2kmsa(T,KL,n, method);
  
  return kmsa2dnd  (S,KL,n);
}

//This function pools into a single file all the children sequences until they contain a maximum of N*2-1 sequences
KT_node*pool (KT_node *K1,int n1,int *n2in, int N)
{
  int a,b, cn;
  int pool=0;
  KT_node*K2;
  int n2;
  int **nseq;

  //for ( a=0; a<n1; a++){if (K1[a]->nseq<N)pool=1;}
  //if (!pool)return NULL;
  
 

  nseq=declare_int(n1, 2);
  for ( a=0; a<n1; a++)
    {
      nseq[a][0]=a;
      nseq[a][1]=K1[a]->nseq;
      
    }
  sort_int (nseq, 2, 1, 0,n1-1);
  

  K2=(KT_node*)vcalloc (n1, sizeof (KT_node));
  for (a=0; a<n1; a++)K2[a]=(KT_node)vcalloc (1, sizeof (KTreenode));
  
  for (cn=0,n2=0, b=0; b< n1; b++)
    {
      if ( cn==0)
	{
	  K2[n2]->seqF=vtmpnam (NULL);
	  K2[n2]->msaF=vtmpnam (NULL);
	  K2[n2]->treeF=vtmpnam (NULL);
	  
	  n2++;
	}
      a=nseq[b][0];
      cn+=K1[a]->nseq;
      
      printf_system ("cat %s >> %s", K1[a]->seqF, K2[n2-1]->seqF);
      K2[n2-1]->nseq=cn;
      
      K1[a]->msaF=K2[n2-1]->msaF;
     
      if (cn>=N)cn=0;
      
    }
  
  free_int (nseq, -1);
  n2in[0]=n2;
  return K2;
}


char* tree2msa4dpa (NT_node T, Sequence *S, int N, char *method)
{
  int n=0;
  char *outname;
  int cn=0;
  int dopool=get_int_variable ("reg_pool");
  
  KT_node K =tree2ktree  (T,T, S, N);
  KT_node*KL=(KT_node*)vcalloc (K->tot, sizeof (KT_node));
  KT_node*KL2;
  int n2=0;
  int docons=0;
  int N2;
  
  n=ktree2klist(K,KL,&n);
  ktree2display (K, "1");
  if (getenv ("DUMP_SEQ_BUCKETS") ||getenv ("DUMP_SEQ_BUCKETS_ONLY"))
    {
      ktree2seq_bucketsF(K, "seqdump.1");
      if (getenv ("DUMP_SEQ_BUCKETS_ONLY"))exit (0);
    }  
 
  //if (docons){N2=N*10;pool=1;}
  //else N2=N;
  
  //This is where the slave MSAs are computed, all at once.
  if ( dopool && (KL2=pool(KL, n, &n2,N))!=NULL)
    {
      int a;
         
      kseq2kmsa(T,KL2,n2, method);
      
      for (a=0; a<n; a++)
	{

	  //Be carefull not to overwrite msaF: it is common to all the pooled buckets
	  char *tmp=vtmpnam (NULL);
	  trim_fastaF_big  ( KL[a]->msaF, KL[a]->seqF,tmp, NULL, NULL, NULL);
	  KL[a]->msaF=tmp;
	  ungap_fastaF_big ( KL[a]->msaF, KL[a]->msaF, 100);
	  vfree(KL2[a]);
	}
      vfree (KL2);
    }
  else
    {
       
      kseq2kmsa(T,KL,n, method);
    }
  if (getenv ("DUMP_ALN_BUCKETS") ||getenv ("DUMP_ALN_BUCKETS_ONLY"))
    ktree2aln_bucketsF(K, "alndump.");
  
  outname=kmsa2msa (K,S,NULL,NULL); 

  declare_aln_node (-1);//Free all the nodes declared
  vfree (KL);
  return outname;
}

	
static int nalncol;
ALNcol* declare_alncol ();
ALNcol* declare_alncol ()
{
  static int n;
  ALNcol *p=(ALNcol*)vcalloc (1, sizeof (ALNcol));
  //p->id=++n; //put this back when debugging pointers
  return p;
}

char *kmsa2msa (KT_node K,Sequence *S, ALNcol***S2,ALNcol*start)
{
  int a, s, c;
  Alignment *A=NULL;
  char *out=NULL;
  FILE *fp;
  ALNcol *end;
  char *output;
  ALNcol *msa;
    
  if (!start)
    {
      S2=(ALNcol***)vcalloc (S->nseq, sizeof (ALNcol**));
      for (s=0; s<S->nseq; s++)S2[s]=(ALNcol**)vcalloc (S->len[s], sizeof (ALNcol*));
      A=quick_read_fasta_aln (A,K->msaF);
      
      start=msa2graph(A,S, S2,start, -1);
      out=vtmpnam (NULL);
    }
  
  for (a=0; a<K->nc; a++)
    {
      int i =name_is_in_hlist((K->child[a])->name,S->name, S->nseq);
      A=quick_read_fasta_aln (A,K->child[a]->msaF);
      start=msa2graph (A,S, S2,start,i);
      kmsa2msa (K->child[a], S, S2,start);
    }
  free_aln (A);
  
  //OUT ionly defined in the parent process
  //This is how the recursion stops
  if (!out) return out;
  
  output=get_string_variable ("output");

  msa=start;
  if (!output || strm (output, "fasta_aln"))
    {
      int nn;
      fp=vfopen (out, "w");
      for (nn=0,s=0; s<S->nseq; s++, nn++)
	{
	  int r=0;
	  msa=start;
	  if (nn==1000)
	    {
	      vfclose (fp);
	      fp=vfopen (out, "a");
	      nn=0;
	    }
	  output_completion (stderr,s,S->nseq, 100, "Final MSA");
	  for (c=0; c<S->len[s]; c++)
	    {
	      S2[s][c]->aa=S->seq[s][c];
	    }
	  fprintf (fp, ">%s\n", S->name[s]);
	  
	  r=0;
	  while (msa->next)
	    {
	      if (msa->aa==0)
		{
		  fprintf (fp, "-");
		}
	      else if (msa->aa>0) 
		{
		  fprintf (fp, "%c",msa->aa);
		  msa->aa=0;
		}
	      r++;
	      msa=msa->next;
	    }
	  
	  fprintf (fp, "\n");
	}
      vfclose (fp);
    }
  else if (strm (output, "fastaz_aln"))
    {
      fp=vfopen (out, "w");
      for (s=0; s<S->nseq; s++)
	{
	  int r=0;
	  ALNcol*msa=start;
	  int cg=0;
	  
	  output_completion (stderr,s,S->nseq, 100, "Final MSA");
	  for (c=0; c<S->len[s]; c++)
	    {
	      S2[s][c]->aa=S->seq[s][c];
	    }
	  fprintf (fp, ">%s\n", S->name[s]);
	  
	  while (msa->next)
	    {
	      if (msa->aa==0){cg++;}
	      else if (msa->aa>0) 
		{
		  if (cg){fprintf (fp, "%d",cg);cg=0;}
		  fprintf (fp, "%c",msa->aa);
		  msa->aa=0;
		}
	      msa=msa->next;
	    }
	  if (cg){fprintf (fp, "%d",cg);cg=0;}
	  fprintf (fp, "\n");
	}
      vfclose (fp);
    }
  else
    {
       myexit (fprintf_error (stderr, "-output=%s is not supported when using -reg [FATAL:%s]", output,PROGRAM));
    }

  if (get_string_variable ("homoplasy"))
    {
       ALNcol*msa=start;
       int homoplasy=0;
       int whomoplasy=0;
       int whomoplasy2=0;
       FILE *fp2;
       unsigned long ngap=0;
       unsigned long ngap2=0;
       int len=0;
       
       while (msa->next)
	 {
	   homoplasy+=msa->homoplasy;
	   whomoplasy+=msa->whomoplasy;	   
	   whomoplasy2+=msa->whomoplasy2;
	   
	   ngap+=msa->ngap;
	   ngap2+=msa->ngap*msa->ngap;
	   msa=msa->next;
	   if (msa->aa>=0)len++;
	 }
       
       fp2=vfopen (get_string_variable ("homoplasy"), "w");
       fprintf ( fp2, "HOMOPLASY: %d\n", homoplasy);
       fprintf ( fp2, "WEIGHTED_HOMOPLASY: %d\n", whomoplasy);
       fprintf ( fp2, "WEIGHTED_HOMOPLASY2: %d\n", whomoplasy2);
       
       fprintf ( fp2, "LEN: %d\n", len-2);//substract start and end of msa data structure
       fprintf ( fp2, "NGAP: %lu\n", ngap);
       fprintf ( fp2, "NGAP2: %ul\n",ngap2);
       vfclose (fp2);
    }

  vfree (start);
  for (s=0; s<S->nseq; s++)vfree (S2[s]);
  vfree (S2);
  
  return out;
}
ALNcol * msa2graph (Alignment *A, Sequence *S, ALNcol***S2,ALNcol*msa,int seq)
{
  //Thread msa in the child section
  int s, c,ir, subseq, a;
  int * lu =(int*) vcalloc (A->nseq, sizeof (int));
  int **pos=(int**)declare_int (A->nseq, A->len_aln);
  ALNcol**graph=(ALNcol**)vcalloc ( A->len_aln, sizeof (ALNcol*));
  ALNcol *lp, *p, *test;
  ALNcol *start, *end, *parent, *child, *lchild, *lparent;
  int compact=-1;
  int check_homoplasy=1;
  int *rescount;
  int *gapcount;
  int  nnseq;
  static int tt;
  static int tt2;
  
  if (A->len_aln==0)return msa;
  nnseq=(msa)?(A->nseq+(msa->next)->nseq-1):A->nseq;

  //1-fill up the look up section: pos
  subseq=-1;

  for (s=0; s<A->nseq; s++)
    {
      int r;
      lu[s]=name_is_in_hlist (A->name[s],S->name, S->nseq);
      for (r=0,c=0; c<A->len_aln; c++)
	{
	  if (A->seq_al[s][c]!='-')pos[s][c]=r++;
	  else pos[s][c]=-1;
	}
      //This matches the index of seq in S - the complete sequence dataset and in A, the subMSA
      if (lu[s]==seq)subseq=s;
    }
  
  //The MSA is a linked list of ALNcol having the length of the full MSA
  //S2 is a list of residues with each residue pointing to a position on MSA
  //This way, MSA is only ONE string of Length MSA as opposed to N Strings of length MSA
  //When A is integrated within MSA, for each column, we look for an S2 residue pointing to an ALNcol in msa. 
  //If there is none, then this is an unmatched column and a column of gpas must be inserted in the rest of the sequences
  
  gapcount=(int*)vcalloc ( A->len_aln, sizeof (int));
  rescount=(int*)vcalloc ( A->len_aln, sizeof (int));
  for ( s=0; s<A->nseq; s++)
    if (s!=subseq)
      for (c=0; c<A->len_aln; c++)
	{
	  if (A->seq_al[s][c]=='-')gapcount[c]++;
	  else rescount[c]++;
	}
     
  if (msa && check_homoplasy)
    {
      int len, r;
      int *rpos=(int*)vcalloc (A->len_aln, sizeof (int));
      
      //Homoplasy 1
      for (len=0,c=0; c<A->len_aln; c++)
	{
	  if (A->seq_al[subseq][c]!='-')rpos[len++]=c;
	}
      
      for (r=0; r<len-1; r++)
	{
	  int d1=rpos[r+1]-rpos[r];
	  int d2=(S2[seq][r+1])->index-(S2[seq][r])->index;
	  if (d1>1 && d2>1)
	    {
	      (S2[seq][r])->homoplasy++;
	      (S2[seq][r])->whomoplasy+=MIN(d1,d2);
	    }
	}

      //whomoplasy2
      //Counts the minimum number of indels on each side (parent and subMSA), keep the lowest value
      for (r=0; r<len-1; r++)
	{
	  int c, g,re;
	  int sub, main;
	  //Count gaps and res in child, keep the lowest count -> parsimony
	  g=re=0;
	  for (c=rpos[r]+1;c<rpos[r+1]; c++)
	    {
	      g+=gapcount[c];
	      re+=rescount[c];
	    }
	  sub=MIN(g,re);
	  
	  //Count gaps and res in main, keep lowest count
	  start=(S2[seq][r])->next;
	  end  =(S2[seq][r+1]);
	  g=re=0;
	  while (start!=end)
	    {
	      g+=start->ngap;
	      re+=start->nres;
	      start=start->next;
	    }
	  main=MIN(g,re);
	  (S2[seq][r])->whomoplasy2+=MIN(main,sub);
	}
      vfree (rpos);
    }

  
  
	

  for (c=0; c<A->len_aln; c++)
    {
      p=NULL;
      //scan A->nseq to check if one of the residue is already mapped onto MSA
      //If such a residue is mapped then the whole column is mapped
      for (s=0; s<A->nseq; s++)
	{
	  ir=pos[s][c];

	  if (ir!=-1 && S2[lu[s]][ir] )
	    {
	      p=S2[lu[s]][ir];graph[c]=p;
	      p->nres+=rescount[c];
	      
	      break;
	    }
	}
      
      //If no residue was found to be mapped for msa
      //all residues of this column will then be mapped onto this MSA position
      if (!p)
	{
	  p=graph[c]=(ALNcol*)declare_alncol();
	  p->nres=rescount[c];
	}
      for (s=0; s<A->nseq; s++)
	{	  
	  ir=pos[s][c];
	  if (ir!=-1 && !S2[lu[s]][ir])S2[lu[s]][ir]=p;
	  //p is an empty column it will trigger a insertion in msa in the next step.
	}
    }
 
  if (!msa)//graph gets turned into msa
    {
      msa=start=declare_alncol();
      end=declare_alncol();
      start->aa=-1;
      end->aa=-1;
      start->next=graph[0];

      for (c=0; c<A->len_aln; c++)
	{
	  graph[c]->next=(c<A->len_aln-1)?graph[c+1]:end;
	}
     
    }
  //Insert graph into msa by expanding seq
  //Trick: graph[A->len_aln] contains all the nodes
  //the ones from seq are already in. They can be recognised with aa==1
  //the new ones are the gaps and they need to be connected to the aa==1
  
  else 
    {
      ALNcol*last=msa;
      int len=0;
      
      for (c=0; c<S->len[seq]; c++)(S2[seq][c])->aa=1;
      for (c=0; c<A->len_aln; c++)
	{
	  if (graph[c]->aa)last=graph[c];
	  else 
	    {
	      //This is where new columns of gaps are inserted in the MSA
	      graph[c]->next=last->next;
	      last->next=graph[c];
	      last=last->next;
	    }
	}
      for (c=0; c<S->len[seq]; c++)(S2[seq][c])->aa=0;
    }
  
  //This assigns an index to every column in the linked list MSA
  
  if (check_homoplasy)
    {
      int i=0;
      ALNcol*st=msa;
      while (st)
	{
	  if (st->aa>=0)
	    {
	      st->nseq=nnseq;
	      st->index=i++;
	      st->ngap=st->nseq - st->nres;
	    }
	  st=st->next;
	}
    }
  
  
  vfree (graph);
  vfree (lu);
  free_int (pos,-1);
  vfree (gapcount); vfree(rescount);
  return msa;
}  

int ktree2aln_bucketsF(KT_node K,char *fname)
{

  if (!K)return 0;
  else
    {
      char *nfname=(char*)vcalloc (1000, sizeof (char));
      int a;
      
      for (a=0; a<K->nc; a++)
	{
	  sprintf (nfname, "%s.%d.aln_bucket",fname, a+1);
	  printf_system ("cp %s %s", K->msaF, nfname);
	  sprintf (nfname, "%s.%d",fname, a+1);
	  ktree2seq_bucketsF (K->child[a],nfname); 
	}
      vfree (nfname);
    }
  return 1;
}
int ktree2parent_seq_bucketsF(KT_node K,char *fname)
{
  if (!K)return 0;
  else
    {
      char *nfname=(char*)vcalloc (1000, sizeof (char));
      int a;
      
      printf_system ("cp %s %s", K->seqF, fname);
    }
  return 1;
}
int ktree2display(KT_node K,char *fname)
{

  if (!K)return 0;
  else
    {
      int a;
      fprintf (stderr, "!\t%-10s -- %10d Seq\n", fname, K->nseq);

      for (a=0; a<K->nc; a++)
	{
	  char *nfname=csprintf(NULL, "%s.%d", fname, a+1);
	  ktree2display(K->child[a],nfname); 
	  vfree (nfname);
	}
    }
  return 1;
}
int ktree2seq_bucketsF(KT_node K,char *fname)
{

  if (!K)return 0;
  else
    {
      int a;
      fprintf (stderr, "!DUMP: \t%-10s -- %10d Seq\n", fname, K->nseq);
      printf_system ("cp %s %s", K->seqF, fname);
      for (a=0; a<K->nc; a++)
	{
	  char *nfname=csprintf(NULL, "%s.%d", fname, a+1);
	  ktree2seq_bucketsF (K->child[a],nfname); 
	  vfree (nfname);
	}

    }
  return 1;
}
    
int ktree2klist (KT_node K, KT_node *KL, int *n)
{
  int a;
  if (!K) return n[0];
  KL[n[0]++]=K;
  for (a=0; a<K->nc; a++)
    {
      ktree2klist (K->child[a],KL,n);
    }
  return n[0];
}

Alignment * sorttrim (Alignment *A,int ntrim)
{
  int **lu, *max, *used;
  int s, c;
  FILE *fp;
  char *tmp=vtmpnam (NULL);
  int maxnseq;
  int newn;
  int bin,nbin=100;
  if (!A) return A;
  else if (A->nseq<=ntrim)return A;
  

  

  lu=declare_int     (nbin, A->nseq);
  max =(int*)vcalloc (nbin, sizeof (int));
  used=(int*)vcalloc (nbin, sizeof (int));
    
  for (maxnseq=0,s=1; s<A->nseq; s++)
    {
      int bin, n, t;
      
      for (t=0, n=0,c=0; c<A->len_aln; c++)
	{
	  char c1=A->seq_al[0][c]; 
	  char c2=A->seq_al[s][c]; 
	  

	  if (c1=='-' || c2=='-')continue;
	  if (c1!=c2)n++;
	  t++;
	}
      if (t==0)continue;
      bin=((n*100)/t)%nbin;
      lu[bin][max[bin]++]=s;
      maxnseq++;
    }
  
  fp=vfopen (tmp, "w");
  fprintf (fp, ">%s\n%s\n",A->name[0], A->seq_al[0]);
  
  newn=0;
  while ( newn<maxnseq && newn<ntrim)
    {
      for (bin=0; bin<nbin && newn<ntrim; bin++)
	{
	  if (used[bin]==max[bin])continue;
	  else s=lu[bin][used[bin]++];
	  fprintf (fp, ">%s\n%s\n",A->name[s], A->seq_al[s]);
	  newn++;
	}
    }
  vfclose (fp);
  
  vfree (max);
  vfree (used);
  free_int (lu, -1);
  return quick_read_fasta_aln (A, tmp);
}



Sequence  * regtrim (Sequence *S, NT_node T, int N)
{
  NT_node *CL, *NL;
  int left, right, nc,nn, a,s,terminal;
  FILE *fp;
  int **ordered;
  char *tmp=vtmpnam(NULL);
  float *w;
  int *used;

  //This defines the number of sequences that will be kept
  
  
  if (S->nseq==1)
    {
      fp=vfopen (tmp, "w");
      fprintf ( fp, ">%s\n%s\n", S->name[0], S->seq[0]);
      vfclose (fp);
      return get_fasta_sequence (tmp,NULL);
    }
  //This is typically set to 1 to keep the first sequence
  int keep=get_int_variable ("keep");
  
  if (!T)
    {
      char *tmode=get_string_variable ("treemode");
      if (tmode)T=seq2dnd(S,tmode);
      else T=seq2dnd(S,"codnd");
    }
  
  if (!T)return NULL;
  
  w=seq2dpa_weight (S, "longuest");
  T=node2master (T, S, w);
  T->leaf=0;
    
  CL=(NT_node*)vcalloc (1, sizeof (NT_node));
  nc=0;CL[nc++]=T;
  
  terminal=0;
  while (nc<N && !terminal)
    {
      
      nn=0;
      terminal=1;
      NL=(NT_node*)vcalloc (nc*2, sizeof (NT_node));
      CL=sort_nodelist4dpa (CL, nc);

      for (left=0; left<nc && (nn+(nc-left))<N; left++)
	{
	  
	  NT_node N=CL[left];
	  if (N->isseq){NL[nn++]=N; T->leaf++;}
	  else
	    {
	      NL[nn++]=N->right;
	      NL[nn++]=N->left;
	      terminal=0;
	      T->leaf+=2;
	    }
	}
      for (a=left; a<nc; a++)
	{
	  NL[nn++]=CL[a];
	  if (!CL[a]->isseq)terminal=0;
	  T->leaf++;
	}
      
      vfree (CL);
      CL=NL;
      nc=nn;
    }
  
  fp=vfopen (tmp, "w");
  //used makes sure the keep first sequences are not duplicated
  used=(int*)vcalloc (S->nseq, sizeof (int));
  if (keep)
    {
      for (s=0; s<keep;s++)
	{
	  fprintf (fp, ">%s\n%s\n", S->name[s],S->seq[s]);
	  used[s]=1;
	}
    }
  if (S->nseq>N)
    {
      ordered=declare_int (nc, 1);
      for (a=0; a<nc; a++)ordered[a][0]=CL[a]->seq;
      sort_int (ordered,1,0, 0, nc-1);
      for (a=0; a<nc; a++)
	{
	  int s=ordered[a][0];
	  
	  if (!used[s])fprintf (fp, ">%s\n%s\n", S->name[s],S->seq[s]);
	}
      vfclose (fp);
      free_int (ordered, -1);
      
    }
  else
    {
      for (a=0; a<nc; a++)
	{
	  
	  int s=CL[a]->seq;
	  
	  if (!used[s])fprintf (fp, ">%s\n%s\n", S->name[s],S->seq[s]);
	}
      vfclose (fp);
    }
  vfree (used);
  vfree (w);
  vfree(CL);
  return get_fasta_sequence (tmp,NULL);
  
}

KT_node tree2ktree (NT_node ROOT,NT_node T,Sequence *S, int N)
{
  NT_node *CL, *NL;
  KT_node K;
  int left, right, nc,nn, a,terminal;
  FILE *fp;
  int **ordered;
  float *w;
  if (N<2)N=2;
  //Dynamic is the variable that will be used to increase the bucket size from root to leaf
  int dynamic=get_int_variable ("reg_dynamic");
  
  if (N==0)printf_exit ( EXIT_FAILURE,stderr, "\nERROR: Bucket size set to 0 (N=%d Dynamic=%d) [FATAL:tree2ktree]", N, dynamic);
  if (!T)return NULL;
  
  T->leaf=0;
  //Make sure the node content is declared
  K=(KT_node)vcalloc (1, sizeof (KTreenode));
  
  CL=(NT_node*)vcalloc (1, sizeof (NT_node));
  nc=0;CL[nc++]=T;
  
  terminal=0;
  while (nc<N && !terminal)
    {
      
      nn=0;
      terminal=1;
      NL=(NT_node*)vcalloc (nc*2, sizeof (NT_node));
      CL=sort_nodelist4dpa (CL, nc);

      for (left=0; left<nc && (nn+(nc-left))<N; left++)
	{
	  
	  NT_node N=CL[left];
	  if (N->isseq){NL[nn++]=N; T->leaf++;}
	  else
	    {
	      NL[nn++]=N->right;
	      NL[nn++]=N->left;
	      terminal=0;
	      T->leaf+=2;
	    }
	}
      for (a=left; a<nc; a++)
	{
	  NL[nn++]=CL[a];
	  if (!CL[a]->isseq)terminal=0;
	  T->leaf++;
	}
      
      vfree (CL);
      CL=NL;
      nc=nn;
    }
  K->nseq=nc;
  
  
  K->seqF=vtmpnam(NULL);
  fp=vfopen (K->seqF, "w");

  

  if (S->nseq>N)
    {
      ordered=declare_int (nc, 1);
      for (a=0; a<nc; a++)ordered[a][0]=CL[a]->seq;
      sort_int (ordered,1,0, 0, nc-1);
      for (a=0; a<nc; a++)
	{
	  int s=ordered[a][0];
	  fprintf (fp, ">%s\n%s\n", S->name[s],S->seq[s]);
	}
     
      free_int (ordered, -1);
    }
  else
    {
      for (a=0; a<nc; a++)
	{
	  
	  int s=CL[a]->seq;
	  
	  fprintf (fp, ">%s\n%s\n", S->name[s],S->seq[s]);
	}
    }
  vfclose (fp);
  
  
  
  K->msaF=vtmpnam(NULL);
    
  K->child=(KT_node*)vcalloc (nc, sizeof (KT_node));
  K->tot=1;
    
  for (a=0; a<nc; a++)
    {
      if (!CL[a]->isseq)
	{
	  //The Children will be at max, dynamic times larger than their parent 
	  K->child[K->nc]=tree2ktree (ROOT,CL[a],S, (dynamic>0)?N*dynamic:(-1*N/dynamic));
	  (K->child[K->nc])->name=(CL[a])->name;
	  K->tot+=(K->child[K->nc])->tot;
	  K->nc++;
	}
      else
	{
	  K->tot++;
	}
    }
  vfree (CL);
  return K;
}
KT_node *free_ktree (KT_node K)
{
  int n;
  
  if (!K) return NULL;
  else if (K->nc==0)vfree (K);
  else 
    {
      for (n=0; n<K->nc; n++)
	free_ktree (K->child[n]);
      vfree ((KT_node**)K->child);
      vfree ((KT_node*)K);
      
    }
  return NULL;
}
  
int kseq2kmsa_serial   (NT_node T,KT_node *K, int n, char *method);
int kseq2kmsa_nextflow (NT_node T,KT_node *K, int n, char *method);
int kseq2kmsa_thread   (NT_node T,KT_node *K, int n, char *method);

int kseq2kmsa   (NT_node T,KT_node *K, int n, char *method)
{
  int nproc=get_nproc();
  
    
  if ( strstr (method, "NF_"))
    return kseq2kmsa_nextflow(T,K, n, method);
  
  return kseq2kmsa_thread (T,K, n, method);
}
  

int kseq2kmsa_nextflow   (NT_node T,KT_node *K, int n, char *met)
{
  int a;
  char *in=vtmpnam (NULL);
  char *out=vtmpnam (NULL);
  TC_method *method=method_file2TC_method(method_name2method_file(met+3));
  char *command=make_aln_command (method,"inputfile", "outputfile");
 
 
  FILE *fp1=vfopen (in, "w");
  FILE *fp2=vfopen (out, "w");
  HERE ("%s", command);
  
  for (a=0; a<n; a++)
    {
      fprintf (fp1, "%s %s\n", K[a]->seqF);
      fprintf (fp2, "%s %s\n", K[a]->msaF);
    }
  vfclose (fp1);
  vfclose (fp2);
  
  //Generate a NF pipleline that runs command on each seqF file so as to generate an msaF
  //Deploy the Nextflow command
  vfree (command);
  
  printf_exit (EXIT_FAILURE, stderr,"NF_<method> is currently NOT supported\n");
  return 1;
}

char *tree2child_tree(NT_node T, char *seqF, char *mode)
{
   
  NT_node C;
  Sequence *S;
  char *treeF=NULL;
  static int a;
  if (!mode)treeF=csprintf (treeF, "default");
  else if (strm(mode, "master") || strm(mode, "parent"))
    {
      S=get_fasta_sequence (seqF, NULL);
      treeF=prune_treeF(T,S,NULL);
      free_sequence (S, -1);
    }
  else treeF=csprintf (treeF, "%s",mode);

return treeF;
}

    

int kseq2kmsa_thread   (NT_node T,KT_node *K, int n, char *method)
{
  int nproc=get_nproc();
  
  int * njobs;
  KT_node **KL;
  int npid;
  int failed=0;
  int a, b;
  int nt=0;
  //split the jobs
  
 
  KL=(KT_node**)vcalloc(nproc+1, sizeof (KT_node*));
  for (b=0; b<nproc; b++)
    KL[b]=(KT_node*)vcalloc ((n/nproc)+5,sizeof (KT_node)); 
  
  njobs=(int*)vcalloc (nproc, sizeof (int));
  for (a=0, b=0; a<n; a++, b++)
    {
      if (b==nproc)b=0;
      KL[b][njobs[b]++]=K[a];
    }
  for (a=0; a<nproc; a++)
    {
      if (!njobs[a]){vfree(KL[a]); KL[a]=NULL;}
    }
  //deploy the jobs
  a=0;
  
  while (KL[a]) 
    {
      nt++;
      if (vvfork(NULL)==0)//child process
	{
	  int b=0;
	  while (KL[a][b])
	    {
	      KT_node LK=KL[a][b];
	      initiate_vtmpnam(NULL);//make sure existing tmp are not deleted when exiting.
	      LK->treeF=tree2child_tree(T,LK->seqF,getenv("child_tree_4_TCOFFEE"));
	      reg_seq_file2msa_file (method,LK->nseq,LK->seqF, LK->msaF, LK->treeF);
	      b++;
	      if ( a==0)
		output_completion (stderr, b, njobs[0], 100, method);
	     
	    }
	  myexit (EXIT_SUCCESS);
	}
      else
	{
	  a++;
	}
    }
  
  //collect the jobs
  a=0;
  while (KL[a++])
    {
      vwait(NULL);
    }
  
  for ( a=0; a<n; a++)failed+=(check_file_exists (K[a]->msaF))?0:1;
  if (failed>0)
    printf_exit ( EXIT_FAILURE,stderr, "\nERROR: method %s failed to produce %d out of %d node alignments [FATAL]", method, failed, n);
  else
    fprintf (stderr, "\n!All Jobs collected\n");
  
  vfree (njobs);
  for (a=0; a<nproc; a++)vfree (KL[a]);
  vfree (KL);
  
  return n;
}
char *reg_seq_file2msa_file (char *method,int nseq, char* seqF, char* msaF, char *treeF)
{
  int use_old_constraint_list=0;
  char *com;
 
  
  if (nseq==1)com=csprintf (NULL, "cp %s %s", seqF, msaF);    
  else if (use_old_constraint_list)return seq_file2msa_file (method,seqF, msaF);
  else com=csprintf (NULL, "dynamic.pl -method %s -seq %s -outfile %s -tree %s", method, seqF, msaF, treeF);
 
  printf_system (com);

  if (!isfile(msaF))
    {
      printf_exit (EXIT_FAILURE, stderr, "ERROR: Impossible to run %s [FATAL:%s]\n",com, PROGRAM);
      return NULL;
    }
  else
    {
      vfree(com);
      return msaF;
    }
}

ALN_node* kmsa2graph2 (Sequence *S,KT_node K,Alignment *A0, ALN_node **lu0, ALN_node *list, int *ns, int *done, int max);
int node2gap_len2 (ALN_node n);
ALN_node ** msa2graph2 (Sequence *S, Alignment *A, ALN_node**lu);
ALN_node insert_gap_in_graph2 (ALN_node start);
int link_nodes2 (ALN_node s, ALN_node m);
int align_graph (ALN_node m, int lm, ALN_node s, int ls);
char *graph2cons (ALN_node m, int len);
char **graph2master_seq (ALN_node n, char **seq);
ALN_node n2FirstSeq (ALN_node n);
ALN_node n2FirstAA (ALN_node n);

ALN_node n2top (ALN_node n);
ALN_node n2bot (ALN_node n);
ALN_node n2end (ALN_node n);
ALN_node n2start (ALN_node n);

int node2gap_len (ALN_node n);
int node2len (ALN_node n);
int node2len_left (ALN_node n);
int node2nseq     (ALN_node n);
int check_graph_integrity (ALN_node g);;
void display_graph_top   (ALN_node n, char *text);
ALN_node insert_gap_in_graph (ALN_node start);
int check_node(ALN_node n, char *txt);

void display_right   (ALN_node n, char *text);
void display_parent   (ALN_node n, char *text);
void display_child   (ALN_node n, char *text);
void display_left   (ALN_node n, char *text);


int insert_msa_in_msa (ALN_node s, ALN_node m);
ALN_node insert_node (ALN_node left, ALN_node right, ALN_node parent, ALN_node child, char value);

ALN_node* kmsa2graph_multi (Sequence *S, KT_node K,Alignment *A0, ALN_node **lu0, ALN_node *list, int *ns, int *done, int max);
ALN_node* kmsa2graph (Sequence *S, KT_node K,Alignment *A0, ALN_node **lu0, ALN_node *list, int *ns, int *done, int max);
ALN_node* kmsa2graph_seq (Sequence *S,KT_node K,Alignment *A0, ALN_node **lu0, ALN_node *list, int *ns, int *done, int max);

ALN_node ** msa2graph (Sequence *S, Alignment *A, ALN_node**lu);
void check_seq_graph (Sequence *S, ALN_node iseq);

void check_aln_graph (Sequence *S, ALN_node *aln, int nseq);


char* graph2aln (Sequence *S, ALN_node *aln, int nseq,char *out);

ALN_node* kmsa2graph_multi (Sequence *S, KT_node K,Alignment *A0, ALN_node **lu0, ALN_node *list, int *ns, int *done, int max)
{
  int nproc=get_nproc();
  int a, b, c;
  reset_output_completion ();
  for (a=0; a<K->nc;)
    {
      for (b=0; b<nproc && a<K->nc;)
	{
	  output_completion (stderr,a,K->nc, 100, "Incorporating the MSAs - [Multi Thread]");
	  if (vvfork(NULL)==0)//child process
	    {
	      ALN_node *aln;
	      initiate_vtmpnam(NULL);
	      aln=kmsa2graph (S,K->child[a],A0,lu0,list,ns,done, max);
	      graph2aln (S,aln,ns[0],(K->child[a])->msaF);
	      myexit (EXIT_SUCCESS);
	    }
	  else
	    {
	      a++;
	      b++;
	    }
	  
	}
      for (c=0; c<b; c++)vwait(NULL);
    }
  reset_output_completion ();
  for (a=0; a<K->nc;a++)(K->child[a])->nc=0;
  return kmsa2graph (S,K,A0,lu0,list,ns,done, max);
}
char* graph2aln (Sequence *S, ALN_node *aln, int nseq,char *out)
{
  FILE *fp=vfopen (out, "w");
  short*lu=(short*)vcalloc (S->nseq, sizeof (short));
  int a;
  
  
  fp=vfopen (out, "w");
  for (a=0; a<nseq; a++)
    {
      ALN_node seq=aln[a];
      
      if (lu[seq->seqN]);
      else
	{
	  lu[seq->seqN]=1;
	  fprintf (fp, ">%s\n", S->name[seq->seqN]);
	  while (seq->aa!=')')
	    {
	      ALN_node pseq=seq;
	      
	      if (seq->aa!='(')fprintf (fp, "%c", seq->aa);
	      seq=seq->r;
	    }
	  fprintf (fp, "\n");
	}
    }
  
  vfclose (fp);
  vfree (lu);
  return out;
}
char *getseq4kmsa2dnd (char *seq, char *alseq, ALNcol *start, ALNcol **S2);
NT_node kmsa2dnd (Sequence *S,KT_node *KL, int n)
{
  char *out=vtmpnam (NULL);
  ALNcol ***S2;
  ALNcol *start, *end;
  int a, s, c;
  FILE*fp=vfopen (out,"w");
  int len_aln;
  NT_node tree;
  char *refseq, *newseq;
  int  **ro;
  char **nname;
  int maxl, refseqI;
  NT_node T;
  
  S2=(ALNcol***)vcalloc (S->nseq, sizeof (ALNcol**));
  for (s=0; s<S->nseq; s++)S2[s]=(ALNcol**)vcalloc (S->len[s], sizeof (ALNcol*));
  
  
  for (a=0; a<n; a++)
    {
      
      Alignment *A=quick_read_aln (KL[a]->msaF);
      int * lu =(int*) vcalloc (A->nseq, sizeof (int));
      int **pos=(int**)declare_int (A->nseq, A->len_aln);
      output_completion (stderr,a,n, 100, "Incorporating Children MSA to produce guide tree");
      for (s=0; s<A->nseq; s++)
	{
	  int r;
	  lu[s]=name_is_in_hlist (A->name[s],S->name, S->nseq);
	  for (r=0,c=0; c<A->len_aln; c++)
	    {
	      if (A->seq_al[s][c]!='-')pos[s][c]=r++;
	      else pos[s][c]=-1;
	    }
	}
      
      for (c=0; c<A->len_aln; c++)
	{
	  ALNcol *p=NULL;
	  for (s=0; s<A->nseq; s++)
	    {
	      int ir=pos[s][c];
	      
	      if (ir!=-1)
		{
		  if (S2[lu[s]][ir])
		    {
		      p=S2[lu[s]][ir];
		      break;
		    }
		}
	    }
	  if (!p)p=(ALNcol*)vcalloc (1, sizeof (ALNcol));
	  for (s=0; s<A->nseq; s++)
	    {
	      int ir=pos[s][c];
	      if (ir!=-1)
		{
		  S2[lu[s]][ir]=p;
		  
		}
	    }
	}
      free_aln (A);
      free_int (pos, -1);
      vfree(lu);
    }
  
 
  start=(ALNcol*)vcalloc (1, sizeof (ALNcol));
  end  =(ALNcol*)vcalloc (1, sizeof (ALNcol));
  start->aa=end->aa=-1;
  start->next=end;
  
  for (s=0; s<S->nseq; s++)
    {
      ALNcol *cpos=start;
      output_completion (stderr,s,S->nseq, 100, "Threading Sequences to produce guide tree");
      for (c=0; c<S->len[s]; c++)
	{
	  
	  ALNcol *p=S2[s][c];
	  if (!p->next)
	    {
	      p->next=cpos->next;
	      cpos->next=p;
	      len_aln++;
	    }
	  cpos=p;
	}
    }
  maxl=refseqI=0;
  for (s=0; s<S->nseq; s++)
    if (S->len[a]>maxl)
      {
	maxl=S->len[a];
	refseqI=a;
      }
  refseq=(char*)vcalloc (len_aln+1, sizeof (char));
  newseq=(char*)vcalloc (len_aln+1, sizeof (char));
  ro=declare_int (S->nseq, 2);
  nname=(char**)vcalloc (S->nseq, sizeof (char*));
  
  getseq4kmsa2dnd(S->seq[refseqI],refseq,start,S2[refseqI]);
  for (a=0; a<S->nseq; a++)
    {
      getseq4kmsa2dnd(S->seq[a],newseq,start, S2[a]);
      ro[a][0]=a;
      ro[a][1]=get_seq_sim (refseq, newseq, "-", "sim1");
      
    }
  sort_int (ro, 2, 1, 0, S->nseq-1);
    
  for (a=0;a<S->nseq; a++)
    nname[a]=S->name[ro[a][0]];
  
  T=list2balanced_dnd (nname, S->nseq);
  vfree (refseq);
  vfree (newseq);
  vfree (nname);
  free_int (ro, -1);
  
  return T;
  
}

char *getseq4kmsa2dnd (char *seq, char *alseq, ALNcol *start, ALNcol **S2)
{
  ALNcol *msa=start;
  int l=strlen (seq);
  int c, p, r;
  
  for (c=0; c<l; c++)
    {
      S2[c]->aa=1;
    }
  p=r=0;
  while (msa->next)
    {
      if (!msa->aa){alseq[p++]='-';}
      else if (msa->aa==1) 
	{
	  alseq[p++]=seq[r++];
	  msa->aa=0;
	}
      msa=msa->next;
    }
  return alseq;
}

ALN_node ** msa2graph (Sequence *S, Alignment *A, ALN_node**lu)
{
  int *cc=(int*)vcalloc (A->nseq, sizeof (int));
  ALN_node **aln;
  int s, c;
  
  if (!lu)
    {
      lu=(ALN_node**)vcalloc (A->nseq, sizeof (ALN_node*));
      for (s=0; s<A->nseq; s++)
	{
	  lu[s]=(ALN_node*)vcalloc (A->len_aln+2, sizeof (ALN_node));
	}
    }
  
  aln=(ALN_node**)vcalloc (A->nseq+4,sizeof (ALN_node*));
  for (s=0;s<A->nseq+4; s++)
    {
      aln[s]=(ALN_node*)vcalloc ( A->len_aln+4, sizeof (ALN_node));
      aln[s]+=2;
    }
  aln+=2;
  for (s=-1; s<=A->nseq; s++)
    for (c=-1; c<=A->len_aln; c++)
      aln[s][c]=declare_aln_node(1);

  
  for (s=-1; s<=A->nseq; s++)
    {
      int seqN=(s<0 || s==A->nseq)?-1:name_is_in_hlist (A->name[s], S->name, S->nseq);
      
      for (c=-1; c<=A->len_aln; c++)
	{
	  ALN_node n=aln[s][c];
	  
	  n->seqN =seqN;
	  n->p    =aln[s-1    ][c  ];
	  n->c    =aln[s+1    ][c  ];
	  n->l    =aln[s      ][c-1];
	  n->r    =aln[s      ][c+1];
	  
	  if      (s==-1      ){n->aa='[';}
	  else if (s== A->nseq){n->aa=']';}
	  else if (c==-1      ){n->aa='(';if (n->seqN!=-1)lu[s][0]=n;}
	  else if (c== A->len_aln){n->aa=')';}
	  else
	    {
	      n->aa=A->seq_al[s][c];
	      if (n->aa!='-')
		{
		  cc[s]++;
		  lu[s][cc[s]]=n;
		}
	    }
	}
      
    }
    
  aln-=2;
  for (s=0; s<A->nseq+4; s++)
    {
      vfree (aln[s]-2);
    }
  vfree (aln);
  vfree (cc);
  return lu;
}
  
ALN_node* kmsa2graph (Sequence *S,KT_node K,Alignment *A0, ALN_node **lu0, ALN_node *list, int *ns, int *done, int max)
{
  int a,b,c, n;
  Alignment *A1;
  
  if (!A0)
    {
      A0=quick_read_aln (K->msaF);
      lu0=msa2graph (S,A0, NULL);
      for (a=0; a<A0->nseq; a++)list[ns[0]++]=lu0[a][0];
      reset_output_completion ();
    }  
  
  for (a=0; a<K->nc; a++)
    {
      ALN_node m, s, t;
      A1 =quick_read_aln ((K->child[a])->msaF);
      ALN_node **lu1=msa2graph (S,A1, NULL);
      int i0=name_is_in_list ((K->child[a])->name,A0->name, A0->nseq,MAXNAMES);
      int i1=name_is_in_list ((K->child[a])->name,A1->name, A1->nseq,MAXNAMES);
      
      for (b=0; b<A1->nseq; b++)list[ns[0]++]=lu1[b][0];
            
      n=0;
      while (lu0[i0][n])
	{
	  int g;
	  int doprint=0;
	  m=lu0[i0][n];
	  s=lu1[i1][n];
	  
	  int lm=node2gap_len (m->r);
	  int ls=node2gap_len (s->r);
	  if (lm>0 || ls>0)
	    {
	      int c;
	      int d,e;
	      char *master=graph2cons(m->r, lm);
	      char *slave =graph2cons(s->r, ls);
	      ALN_node tm=n2top(m);
	      ALN_node ts=n2top(s);
	      static Alignment *A;
	      ALN_node buf;
	      
	      
	      A=align_two_streches4dpa (master, slave, "blosum62mt",-4,-1, "myers_miller_pair_wise", A);
	      for (c=0; c<A->len_aln; c++)
		{
		  if (A->seq_al[0][c]=='-')tm=insert_gap_in_graph (tm);
		  else tm=tm->r;
		  if (A->seq_al[1][c]=='-')ts=insert_gap_in_graph (ts);
		  else ts=ts->r;
		}
	      vfree (master); vfree (slave);
	    }
	  n++;
	}
      
      done[0]+=1;
      //output_completion (stderr,done[0],max, 100, "Incorporating MSAs");
      insert_msa_in_msa (lu1[i1][0],lu0[i0][0]);
      list=kmsa2graph(S,K->child[a], A1, lu1, list, ns, done, max);
      
      for (b=0; b<A1->nseq; b++)vfree(lu1[b]);
      vfree (lu1);
      free_aln (A1);
    }
  
  return list;
}


ALN_node* kmsa2graph_seq (Sequence *S,KT_node K,Alignment *A0, ALN_node **lu0, ALN_node *list, int *ns, int *done, int max)
{
  int a,b,c, n;
  Alignment *A1;
  
  if (!A0)
    {
      A0=quick_read_aln (K->msaF);
      lu0=msa2graph (S,A0, NULL);
      for (a=0; a<A0->nseq; a++)list[ns[0]++]=lu0[a][0];
      reset_output_completion ();
    }  
  
  for (a=0; a<K->nc; a++)
    {
      static Alignment *A;
      static char **seq0;
      static char **seq1;
      
      ALN_node m, s, t;
      A1 =quick_read_aln ((K->child[a])->msaF);
      ALN_node **lu1=msa2graph (S,A1, NULL);
      int i0=name_is_in_list ((K->child[a])->name,A0->name, A0->nseq,MAXNAMES);
      int i1=name_is_in_list ((K->child[a])->name,A1->name, A1->nseq,MAXNAMES);
      
      for (b=0; b<A1->nseq; b++)list[ns[0]++]=lu1[b][0];

      seq0=graph2master_seq (lu0[i0][0], seq0);
      seq1=graph2master_seq (lu1[i1][0], seq1);
      
      A=align_two_sequences4dpa (seq0[0],seq0[1],seq1[0], seq1[1],"blosum62mt",-4,-1, "myers_miller_pair_wise", A);
      
      
      m=n2top(lu0[i0][0]);
      s=n2top(lu1[i1][0]);
      
      for (b=0; b<A->len_aln; b++)
	{
	  if (A->seq_al[0][b]=='-')m=insert_gap_in_graph (m);
	  else m=m->r;
	  
	  if (A->seq_al[1][b]=='-')s=insert_gap_in_graph (s);
	  else
	    s=s->r;
	}
                  
      done[0]+=1;
      //output_completion (stderr,done[0],max, 100, "Incorporating MSAs");
      insert_msa_in_msa (lu1[i1][0],lu0[i0][0]);
      list=kmsa2graph(S,K->child[a], A1, lu1, list, ns, done, max);
      
      for (b=0; b<A1->nseq; b++)vfree(lu1[b]);
      vfree (lu1);
      free_aln (A1);
    }
  
  return list;
}

char **graph2master_seq (ALN_node n, char **seq)
  {
    int l=node2len (n);
    static int *aa=(int*)vcalloc(256, sizeof (int));
    ALN_node t;

    if (!seq) seq=(char**)vcalloc (2, sizeof (char*));
    seq[0]=(char*)vrealloc(seq[0], (l+2)*sizeof (char));
    seq[1]=(char*)vrealloc(seq[1], (l+2)*sizeof (char));
    
    t=n2FirstAA(n);
    n=n->r;
    
    l=0;
    while (t->r)
      {
	seq[0][l]=seq[1][l]=n->aa;
	if (n->aa=='-')
	  {
	    int a;
	    ALN_node s=t;
	    int best_naa=0;
	    char best_aa=0;
	    for (a=0; a<256; a++)aa[a]=0;
	    while (s->c)
	      {
		char r=s->aa;
		aa[r]++;
		
		if (aa[r]>best_naa){best_naa=aa[r]; best_aa=r;}
		s=s->c;
	      }
	    seq[0][l]=(best_aa=='-')?'X':best_aa;
	  }
	t=t->r;
	n=n->r;
	l++;
      }
    seq[0][l]=seq[1][l]='\0';
    return seq;
  }
    
char *graph2cons (ALN_node m, int len)
{
  char *master=(char*)vcalloc (len+1, sizeof (char));
  static int *aa=(int*)vcalloc(256, sizeof (int));
  int a;
  int l=0;
  m=n2top(m);
  m=m->c;
  while (l<len)
    {
      ALN_node s=m;
      int best_naa=0;
      char best_aa=0;
      for (a=0; a<256; a++)aa[a]=0;
      
      while (s->c)
	{
	  char r=tolower(s->aa);
	  aa[r]++;
	  
	  if (aa[r]>best_naa){best_naa=aa[r]; best_aa=r;}
	  s=s->c;
	}
      master[l++]=(best_aa=='-')?'x':best_aa;
      //master[l++]='x';
      
      m=m->r;
    }
  master[l]='\0';
  return master;
}
int insert_msa_in_msa (ALN_node s, ALN_node m)
{
  int n=0;
  ALN_node botM, topM, botS, topS,lastM, firstM, lastS, firstS, is, im, aln;
  

  m=n2FirstSeq(m);
  s=n2FirstSeq(s);
  
  lastM=n2bot(m);lastM=lastM->p;
  firstS=n2top(s); firstS=firstS->c;
  lastS=n2bot (s); lastS=lastS->p;
    
  while (lastM)
    {
      ALN_node bot=lastM->c;
      ALN_node newn=firstS;
      
      lastM->c=firstS;
      firstS->p=lastM;
      
      lastS->c=bot;
      bot->p=lastS;
            
      lastM=lastM->r;
      firstS=firstS->r;
      lastS=lastS->r;
    }
  return n;
}

//Node navigation
ALN_node n2FirstAA (ALN_node n)
{
  n=n2top(n);
  n=n2start(n);
  n=n->c;
  n=n->r;
  return n;
}
ALN_node n2FirstSeq (ALN_node n)
{
  n=n2top(n);
  n=n->c;
  return n2start(n);
}
ALN_node n2top (ALN_node n)
{
  while (n->p)n=n->p;
  return n;
}

ALN_node n2bot (ALN_node n)
{
  while (n->c)n=n->c;
  return n;
}
ALN_node n2end (ALN_node n)
{
  while (n->r)n=n->r;
  return n;
}
ALN_node n2start (ALN_node n)
{
  while (n->l)n=n->l;
  return n;
}
int node2nseq     (ALN_node n)
{
  int nseq=0;
  n=n2top(n);
  n=n->c;
  
  while (n->c)
    {
      nseq++;
      n=n->c;
    }
  return nseq;
}
int node2len (ALN_node n)
{
  int len=0;
  
  while (n){len++, n=n->r;}
  return len;
}
int node2len_left (ALN_node n)
{
  int len=0;
  
  while (n){len++, n=n->l;}
  return len;
}
int node2gap_len (ALN_node n)
{
  int len=0;
  
  while (n && n->aa=='-'){len++, n=n->r;}
  return len;
}

ALN_node insert_gap_in_graph (ALN_node start)
{
  ALN_node pn=NULL;
  ALN_node n, top;
  if (!start) return NULL;

  while (start)
    {
      n=declare_aln_node(1);
      
      //parent
      n->p=pn;
      if (pn)pn->c=n;
      else top=n;
      
      //chilren
      n->c=NULL;
      
      //right
      n->r=start->r;
      if (n->r)(n->r)->l=n;
      
      //left
      n->l=start;
      start->r=n;

      if      (!start->p)n->aa='[';
      else if (!start->c)n->aa=']';
      else n->aa='-';
      
      pn=n;
      start=start->c;
    }
  return top;
}

      


  

ALN_node declare_aln_node (int mode)
{

  ALN_node r;
  int chunk=1000;
  static ALN_node *bufbuf;
  static int bb;
  static int available;
  static ALN_node buf;


  
  if (mode==-1)
    {
      int a;
      for (a=0; a<bb; a++)
	{
	  vfree (bufbuf[a]);
	}
      vfree (bufbuf);
      bufbuf=NULL;
      available=0;
      buf=NULL;
      return NULL;
    }
  
  if (!available)
    {
      int s;

      buf=(ALN_node)vmalloc (sizeof (alnnode)*chunk);
      available=chunk;
      
      bufbuf=(ALN_node*)vrealloc (bufbuf,(bb+1)*sizeof (ALN_node));
      bufbuf[bb++]=buf;
    }
  available--;
  r=buf++;
  return r;
}   

    
  
//checking Function
int check_node(ALN_node n, char *txt)
{
  ALN_node b;
  b=n;
  int count;
  
  while (b)
    {
      if (b->p)
	if ((b->p)->c!=b){HERE ("pb1: %s", txt); exit (0);}
      b=b->p;
    }
  b=n;
  count=0;
  while (b)
    {
      count ++;
      if (b->c)
	if ((b->c)->p!=b){HERE ("pb2 (%d): %s [] ", count,txt); exit (0);}
      b=b->c;
    }
  b=n;
  while (b)
    {
      if (b->r)
	if ((b->r)->l!=b){HERE ("pb3: %s", txt); exit (0);}
      b=b->r;
    }
  b=n;
  while (b)
    {
      if (b->l)
	if ((b->l)->r!=b){HERE ("pb4: %s", txt); exit (0);}
      b=b->r;
    }
}
  
int check_node_integrity (ALN_node n, char *txt)
{
  while (n)
    {
      int l1=node2len (n);
      int l2=node2len (n2top(n));
      int l3=node2len (n2bot(n));
      
      if (l1!= l2){HERE ("top %d %d", l1, l2);HERE ("%s", txt);exit (0);}
      if (l1!= l3){HERE ("bot %d %d", l1, l3);HERE ("%s", txt);exit (0);}
      n=n->r;
    }
  return 1;
}


int check_graph_integrity (ALN_node g)
{
  ALN_node t, top, bot, col;
  int l;
  // check all the sequences have the same length
  t=n2top(g);
  l=node2len (g);
  while (t)
    {
      int l2=node2len (g);
      if (l2!=l){HERE ("different length"); exit (0);}
      t=t->c;
    }
  // check the top/bottom
  top=t=n2top(g);
  bot=n2bot(g);
  l=node2len (g);
  while (t)
    {
      col=t;
      bot=n2bot(t);
      top=n2top(t);
      while (col)
	{
	  if      (n2top(col)!=top){HERE ("column pb with top");exit (0);}
	  else if (n2bot(col)!=bot){HERE ("column pb with bot");exit (0);}
	  col=col->c;
	}
      t=t->r;
    }
  return 1;
}


void display_graph_top   (ALN_node n, char *text)
{
  if (text)fprintf (stdout, "\n%s\n", text);
  n=n2start(n);
  n=n2top(n);
  while (n)
    {
      ALN_node s=n;
      while (s)
	{
	  fprintf (stdout, "%c", s->aa);
	  s=s->c;
	}
      fprintf (stdout, "\n");
      n=n->r;
    }
}


void display_right   (ALN_node n, char *text)
{

  while (n)
    {
      fprintf (stdout, "%c", n->aa);
      n=n->r;
    }
  fprintf ( stdout, " ::: %s\n",text);
}


void display_parent   (ALN_node n, char *text)
{
  if (text)fprintf (stdout, "\n%s\n", text);
  while (n)
    {
      fprintf (stdout, "%c", n->aa);
      n=n->p;
    }
  fprintf ( stdout, "\n");
}


void display_child   (ALN_node n, char *text)
{
  if (text)fprintf (stdout, "\n%s\n", text);
  while (n)
    {
      fprintf (stdout, "%c", n->aa);
      n=n->c;
    }
  fprintf ( stdout, "\n");
}

void display_left   (ALN_node n, char *text)
{
  if (text)fprintf (stdout, "\n%s\n", text);
  while (n)
    {
      fprintf (stdout, "%c", n->aa);
      n=n->l;
    }
  fprintf ( stdout, "\n");
}
      
void check_aln_graph (Sequence *S, ALN_node *aln, int nseq)
{
  int a; 
  for (a =0; a<nseq; a++)check_seq_graph (S, aln[a]);
}


void check_seq_graph (Sequence *S, ALN_node iseq)
{
  ALN_node seq=iseq;
  while (seq->aa!=')')
    {
      fprintf ( stdout, "%c", seq->aa);
      seq=seq->r;
    }
  fprintf (stdout, " %s\n",S->name[seq->seqN]);
  
  seq=iseq;
  while (seq->aa!=')')
    {
      if      (!seq->p      )fprintf (stdout, "*");
      else fprintf (stdout, "p");
      seq=seq->r;
    }
  fprintf (stdout, " %s\n",S->name[seq->seqN]);
  
  seq=iseq;
  while (seq->aa!=')')
    {
      if      (!seq->c     )fprintf (stdout, "?");
      else fprintf (stdout, "1");
      seq=seq->r;
    }
  fprintf (stdout, " %s\n\n",S->name[seq->seqN]);
}

////////on disc


unsigned long* kmsa2graph_d (Sequence *S,KT_node K,Alignment *A0, unsigned long **lu0, unsigned long *list, int *ns, int *done, int max);
unsigned long ** msa2graph_d (Sequence *S, Alignment *A, unsigned long**lu);
int node2gap_len_d (unsigned long n);
char *graph2cons_d (unsigned long m, int len);
unsigned long insert_gap_in_graph_d (unsigned long start);
int insert_msa_in_msa_d (unsigned long s, unsigned long m);
unsigned long n2FirstSeq_d (unsigned long n);
unsigned long n2start_d (unsigned long n);

unsigned long n2top_d (unsigned long n);
unsigned long n2bot_d (unsigned long n);

unsigned long rn_c    (unsigned long n);
unsigned long rn_p    (unsigned long n);
unsigned long rn_l    (unsigned long n);
unsigned long rn_r    (unsigned long n);
char       rn_aa   (unsigned long n);
int        rn_seqN (unsigned long n);

void wn_c    (unsigned long n, unsigned long t);
void wn_p    (unsigned long n, unsigned long t);
void wn_l    (unsigned long n, unsigned long t);
void wn_r    (unsigned long n, unsigned long t);
void wn_aa   (unsigned long n, char aa);
void wn_seqN (unsigned long n, int seqN);

ALN_node d2m(unsigned long n);
unsigned long m2d(ALN_node m,unsigned long d );

unsigned long declare_aln_node_d (int i);


char * kmsa2msa_d (Sequence *S,KT_node K, int max, int *cn)
{
  char *out=vtmpnam (NULL);
  int  a;
  int nseq=0;
  FILE*fp;
  unsigned long *aln;
  short *lu;
  int done=1;
  aln=(unsigned long*)vcalloc ((2*S->nseq)+1, sizeof (unsigned long));
  lu=(short*)vcalloc (S->nseq, sizeof (short));
  aln=kmsa2graph_d (S,K, NULL, NULL,aln,&nseq,&done,max);
  
  
  fp=vfopen (out, "w");
  for (a=0; a<nseq; a++)
    {
      unsigned long seq=aln[a];
      int seqN=rn_seqN(seq);
      if (lu[seqN]);
      else
	{
	  char aa;
	  lu[rn_seqN(seq)]=1;
	  fprintf (fp, ">%s\n", S->name[seqN]);
	  while ((aa=rn_aa(seq))!=')')
	    {
	      if (aa!='(')fprintf (fp, "%c",aa);
	      seq=rn_r(seq);
	    }
	  fprintf (fp, "\n");
	}
    }
  
  vfclose (fp);
  vfree (aln);
  vfree (lu);
  return out;
}
unsigned long* kmsa2graph_d (Sequence *S,KT_node K,Alignment *A0, unsigned long **lu0, unsigned long *list, int *ns, int *done, int max)
{
  int a,b,c, n;
  Alignment *A1;
  
  if (!A0)
    {
      A0=quick_read_aln (K->msaF);
      lu0=msa2graph_d (S,A0, NULL);
      for (a=0; a<A0->nseq; a++)list[ns[0]++]=lu0[a][0];
    }  
  
  for (a=0; a<K->nc; a++)
    {
      unsigned long m, s, t;
      A1 =quick_read_aln ((K->child[a])->msaF);
      unsigned long **lu1=msa2graph_d (S,A1, NULL);
      int i0=name_is_in_list ((K->child[a])->name,A0->name, A0->nseq,MAXNAMES);
      int i1=name_is_in_list ((K->child[a])->name,A1->name, A1->nseq,MAXNAMES);
      
      for (b=0; b<A1->nseq; b++)list[ns[0]++]=lu1[b][0];
            
      n=0;
      while (lu0[i0][n])
	{
	  int g;
	  
	  m=lu0[i0][n];
	  s=lu1[i1][n];
	  
	  int lm=node2gap_len_d (rn_r(m));
	  int ls=node2gap_len_d (rn_r(s));
	  if (lm>0 || ls>0)
	    {
	      int c;
	      char *master=graph2cons_d(rn_r(m), lm);
	      char *slave =graph2cons_d(rn_r(s), ls);
	      unsigned long tm=n2top_d(m);
	      unsigned long ts=n2top_d(s);
	      static Alignment *A;
	      ALN_node buf;
	      
	      
	      A=align_two_streches4dpa (master, slave, "blosum62mt",-4,-1, "myers_miller_pair_wise", A);
	      for (c=0; c<A->len_aln; c++)
		{
		  if (A->seq_al[0][c]=='-')tm=insert_gap_in_graph_d (tm);
		  else tm=rn_r(tm);
		  if (A->seq_al[1][c]=='-')ts=insert_gap_in_graph_d (ts);
		  else ts=rn_r(ts);
		}
	      vfree (master); vfree (slave);
	    }
	  n++;
	}
      
      done[0]+=1;
      output_completion (stderr,done[0],max, 100, "Incorporating MSAs");
      insert_msa_in_msa_d (lu1[i1][0],lu0[i0][0]);
      list=kmsa2graph_d(S,K->child[a], A1, lu1, list, ns, done, max);
      
      for (b=0; b<A1->nseq; b++)vfree(lu1[b]);
      vfree (lu1);
      free_aln (A1);
    }
  
  return list;
}
unsigned long ** msa2graph_d (Sequence *S, Alignment *A, unsigned long**lu)
{
  int *cc=(int*)vcalloc (A->nseq, sizeof (int));
  unsigned long **aln;
  int s, c;
  
  if (!lu)
    {
      lu=(unsigned long**)vcalloc (A->nseq, sizeof (unsigned long*));
      for (s=0; s<A->nseq; s++)
	{
	  lu[s]=(unsigned long*)vcalloc (A->len_aln+2, sizeof (unsigned long));
	}
    }
  
  aln=(unsigned long**)vcalloc (A->nseq+4,sizeof (unsigned long*));
  for (s=0;s<A->nseq+4; s++)
    {
      aln[s]=(unsigned long*)vcalloc ( A->len_aln+4, sizeof (unsigned long));
      aln[s]+=2;
    }
  aln+=2;
  for (s=-1; s<=A->nseq; s++)
    for (c=-1; c<=A->len_aln; c++)
      aln[s][c]=declare_aln_node_d(1);

  
  for (s=-1; s<=A->nseq; s++)
    {
      int seqN=(s<0 || s==A->nseq)?-1:name_is_in_hlist (A->name[s], S->name, S->nseq);
      
      for (c=-1; c<=A->len_aln; c++)
	{
	  unsigned long n=aln[s][c];
	  
	  wn_seqN(n,seqN);
	  wn_p(n,aln[s-1    ][c  ]);
	  wn_c(n,aln[s+1    ][c  ]);
	  wn_l(n,aln[s      ][c-1]);
	  wn_r(n,aln[s      ][c+1]);
	  
	  if      (s==-1      ){wn_aa(n,'[');}
	  else if (s== A->nseq){wn_aa(n,']');}
	  else if (c==-1      )
	    {
	      
	      wn_aa(n,'(');
	      if (seqN!=-1)lu[s][0]=n;
	    }
	  else if (c== A->len_aln){wn_aa(n,')');}
	  else
	    {
	      wn_aa(n,A->seq_al[s][c]);
	      if (A->seq_al[s][c]!='-')
		{
		  cc[s]++;
		  lu[s][cc[s]]=n;
		}
	    }
	}
      
    }
  
  aln-=2;
  for (s=0; s<A->nseq+4; s++)
    {
      vfree (aln[s]-2);
    }
  vfree (aln);
  vfree (cc);
  return lu;
}
int node2gap_len_d (unsigned long n)
{
  int len=0;
  
  while (n && rn_aa(n)=='-'){len++, n=rn_r(n);}
  return len;
}

char *graph2cons_d (unsigned long m, int len)
{
  char *master=(char*)vcalloc (len+1, sizeof (char));
  static int *aa=(int*)vcalloc(256, sizeof (int));
  int a;
  int l=0;
  m=n2top_d(m);
  m=rn_c(m);
  while (l<len)
    {
      unsigned long s=m;
      unsigned long c;
      int best_naa=0;
      char best_aa=0;
      for (a=0; a<256; a++)aa[a]=0;
      
      while ((c=rn_c(s))!=NULL)
	{
	  char r=tolower(rn_aa(s));
	  aa[r]++;
	  
	  if (aa[r]>best_naa){best_naa=aa[r]; best_aa=r;}
	  s=c;
	}
      master[l++]=(best_aa=='-')?'x':best_aa;
      m=rn_r(m);
    }
  master[l]='\0';
  return master;
}


unsigned long insert_gap_in_graph_d (unsigned long start)
{
  unsigned long pn=NULL;
  unsigned long n, top;
  if (!start) return NULL;

  while (start)
    {
      unsigned long rr;
      n=declare_aln_node_d(1);
      
      //parent
      wn_p(n,pn);
      if (pn)wn_c(pn,n);
      else top=n;
      
      //chilren
      wn_c(n,NULL);
      
      //right
      wn_r(n, rn_r(start));
      rr=rn_r(n);
      if (rr)wn_l(rr,n);
      
      //left
      wn_l(n,start);
      wn_r(start,n);

      if      (!rn_p(start))wn_aa(n,'[');
      else if (!rn_c(start))wn_aa(n,']');
      else wn_aa(n,'-');
      
      pn=n;
      start=rn_c(start);
    }
  return top;
}

int insert_msa_in_msa_d (unsigned long s, unsigned long m)
{
  int n=0;
  unsigned long botM, topM, botS, topS,lastM, firstM, lastS, firstS, is, im, aln;
  

  m=n2FirstSeq_d(m);
  s=n2FirstSeq_d(s);
  
  lastM =n2bot_d(m); lastM =rn_p(lastM);
  firstS=n2top_d(s); firstS=rn_c(firstS);
  lastS =n2bot_d(s); lastS =rn_p(lastS);
  
  while (lastM)
    {
      unsigned long bot =rn_c(lastM);
            
      wn_c(lastM,firstS);
      wn_p(firstS,lastM);
      
      wn_c(lastS,bot);
      wn_p(bot,lastS);
            
      lastM =rn_r(lastM);
      firstS=rn_r(firstS);
      lastS =rn_r(lastS);
    }
  return n;
}
unsigned long n2FirstSeq_d (unsigned long n)
{
  n=n2top_d(n);
  n=rn_c(n);
  return n2start_d(n);
}
unsigned long n2start_d (unsigned long n)
{
  unsigned long l;
  while ((l=rn_l(n)))n=l;
  return n;
}

unsigned long n2top_d (unsigned long n)
{
  unsigned long p;
  while ((p=rn_p(n)))n=p;
  return n;
}
unsigned long n2bot_d (unsigned long n)
{
  unsigned long c;
  while ((c=rn_c(n)))n=c;
  return n;
}

unsigned long rn_c (unsigned long n)
{
  ALN_node nn=d2m(n);
  return (unsigned long)nn->c;
}
unsigned long rn_p (unsigned long n)
{
  ALN_node nn=d2m(n);
  return (unsigned long)nn->p;
}
unsigned long rn_l (unsigned long n)
{
  ALN_node nn=d2m(n);
  return (unsigned long)nn->l;
}
unsigned long rn_r (unsigned long n)
{
  ALN_node nn=d2m(n);
  return (unsigned long)nn->r;
}
char rn_aa (unsigned long n)
{
  ALN_node nn=d2m(n);
  return nn->aa;
}
int rn_seqN (unsigned long n)
{
  ALN_node nn=d2m(n);
  return nn->seqN;
}

void wn_c (unsigned long n, unsigned long t)
{
  ALN_node nn=d2m(n);
  nn->c=(ALN_node)t;
  m2d(nn,n);
}
void wn_p (unsigned long n, unsigned long t)
{
  ALN_node nn=d2m(n);
  nn->p=(ALN_node)t;
  m2d(nn,n);
}
void wn_l (unsigned long n, unsigned long t)
{
  ALN_node nn=d2m(n);
  nn->l=(ALN_node)t;
  m2d(nn,n);
}
void wn_r (unsigned long n, unsigned long t)
{
  ALN_node nn=d2m(n);
  nn->r=(ALN_node)t;
  m2d(nn,n);
}

void wn_aa (unsigned long n, char t)
{
  ALN_node nn=d2m(n);
  nn->aa=t;
  m2d(nn,n);
}
void wn_seqN (unsigned long n, int t)
{
  ALN_node nn=d2m(n);
  nn->seqN=t;
  m2d(nn,n);
}

ALN_node d2m(unsigned long n)
{
  return (ALN_node)n;
}
unsigned long m2d(ALN_node m,unsigned long d )
{
  return d;
}

unsigned long declare_aln_node_d (int i)
{
  return (unsigned long)declare_aln_node (i);
}
#ifdef ONDISC
FILE *disc;
ALN_node d2m(unsigned long d)
{
  static ALN_node m=declare_aln_node (1);
 
  if (!disc)disc=tmpfile();
  fseek (disc, sizeof(ALNnode)*(d-1), SEEK_SET);
  fread (m   , sizeof(ALNnode), 1, disc);
  
  return m;
}
unsigned long m2d(ALN_node m, unsigned long d)
{
 
  if (!disc)disc=tmpfile();
  fseek  (disc, sizeof(ALNnode)*(d-1), SEEK_SET);
  fwrite (m   , sizeof(ALNnode), 1, disc);
  
  return d;
}
unsigned long declare_aln_node_d_old (int i)
{
  static ALN_node n=declare_aln_node (i);
  static unsigned long d;
  return m2d (n,++d);
}



unsigned long declare_aln_node_d (int mode)
{

  int chunk=100000;
  static int available;
  static int tot;
  static unsigned long allocated;
  
  if (!disc)disc=tmpfile();
  
  
  if (!available)
    {
      static ALN_node *buf=(ALN_node*)vcalloc (chunk, sizeof (alnnode));
      fseek  (disc, sizeof(ALNnode)*(tot-1), SEEK_SET);
      fwrite (buf , sizeof(ALNnode),chunk+1, disc);
      tot+=(chunk+1);
      available=chunk;
      HERE ("DEclared chunk: %d", tot);
    }
  available--;
  return ++allocated;
}   

ALN_node* kmsa2graph2 (Sequence *S,KT_node K,Alignment *A0, ALN_node **lu0, ALN_node *list, int *ns, int *done, int max)
{
  int a,b,c, n;
  Alignment *A1;
  
  if (!A0)
    {
      A0=quick_read_aln (K->msaF);
      lu0=msa2graph (S,A0, NULL);
      for (a=0; a<A0->nseq; a++)list[ns[0]++]=lu0[a][0];
      reset_output_completion ();
    }  
  
  for (a=0; a<K->nc; a++)
    {
      ALN_node m, s, t;
      A1 =quick_read_aln ((K->child[a])->msaF);
      ALN_node **lu1=msa2graph (S,A1, NULL);
      int i0=name_is_in_list ((K->child[a])->name,A0->name, A0->nseq,MAXNAMES);
      int i1=name_is_in_list ((K->child[a])->name,A1->name, A1->nseq,MAXNAMES);
      
      for (b=0; b<A1->nseq; b++)list[ns[0]++]=lu1[b][0];
            
      n=0;
      while (lu0[i0][n])
	{
	  int g;
	  
	  m=lu0[i0][n];
	  s=lu1[i1][n];
	  
	  int lm=node2gap_len (m->r);
	  int ls=node2gap_len (s->r);
	  if (lm>0 || ls>0)
	    {
	      int c;
	      char *master=graph2cons(m->r, lm);
	      char *slave =graph2cons(s->r, ls);
	      ALN_node tm=n2top(m);
	      ALN_node ts=n2top(s);
	      static Alignment *A;
	      ALN_node buf;
	      
	      
	      A=align_two_streches4dpa (master, slave, "blosum62mt",-4,-1, "myers_miller_pair_wise", A);
	      for (c=0; c<A->len_aln; c++)
		{
		  if (A->seq_al[0][c]=='-')tm=insert_gap_in_graph (tm);
		  else tm=tm->r;
		  if (A->seq_al[1][c]=='-')ts=insert_gap_in_graph (ts);
		  else ts=ts->r;
		}
	      vfree (master); vfree (slave);
	    }
	  n++;
	}
      
      done[0]+=1;
      //output_completion (stderr,done[0],max, 100, "Incorporating MSAs");
      insert_msa_in_msa (lu1[i1][0],lu0[i0][0]);
      list=kmsa2graph(S,K->child[a], A1, lu1, list, ns, done, max);
      
      for (b=0; b<A1->nseq; b++)vfree(lu1[b]);
      vfree (lu1);
      free_aln (A1);
    }
  
  return list;
}
#endif
