

typedef struct treesim{
  float w;
  float uw;
  float d;
  
  float max_w;
  float max_uw;
  float max_d;
  
  int rf;
  int n;//n nodes;
  int nseq;// nseq in the common subset
    }Tree_sim;


typedef struct tnode *NT_node;

/**
* Node of a tree
*/
typedef struct tnode{
  int visited;
  char *name;
  char *file;
  
  ///The parent node
  NT_node parent;
  ///Left child node
  NT_node left;
  ///Right child node
  NT_node right;
  NT_node bot;
  /// is leaf?
  int isseq;
  int seq;
  int maxnseq;
  int nseq;
  
  ///contains a list of the sequences
  int *lseq; 
  ///contains a coded version of the node: 10010101
  int *lseq2;
  ///contains distances to the root, in nodes
  int *idist;
  ///contains real distances *1000
  int *ldist;
  float dist;
  float bootstrap;
  float dp;
  int order;
  int aligned;
  ///Number of leave below the considered node
  int leaf;
  ///Number of nodes below the considered node
  int node;
  int group;
  float score;
  int align;
  char *seqal;
  int index;
  int fork;
    }Treenode;

typedef struct split_struc Split;

typedef struct split_struc{
  char *split;
  int n;
  int tot;
  float score;
  char **tlist;//Not used yet
  Sequence *S;
  NT_node *L;
}Split_struc;

NT_node main_prune_tree ( NT_node T, Sequence *S);
NT_node prune_tree ( NT_node T, Sequence *S);
/*********************************************************************/
/*                                                                   */
/*                                   dpa_tree_manipulation           */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
char *tree2Ngroup (Alignment *A, NT_node T, int max_n, char *fname, char *mat4dist);
int tree2group_file ( NT_node T,Sequence *S, int maxnseq, int minsim, char *name);

NT_node seq2dpa_tree  (Sequence *S, char *align_mode);
NT_node tree2dpa_tree (NT_node T, Alignment *A, char *matrix4distance);
FILE * tree2group ( NT_node T,Sequence *S,int maxnseq, int mindist,char *name, FILE *fp);

NT_node  collapse_tree (NT_node T, Sequence *S, char *string);
NT_node  tree2collapsed_tree (NT_node T, int n, char **string);

/*********************************************************************/
/*                                                                   */
/*                                   tree comparison                 */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
int main_compare_cog_tree (NT_node T1, char *cogfile);
int main_compare_aln_tree (NT_node T1, Alignment *A, FILE *fp);
int compare_aln_tree (NT_node T, Alignment *A, int *n, FILE *fp);

int main_compare_splits (NT_node T1, NT_node T2, char *mode, FILE *fp);
Tree_sim * tree_cmp( NT_node T1, NT_node T2);
NT_node tree_scan (Alignment *A,NT_node RT, char *pscan, char *ptree);

int print_node_list (NT_node T, Sequence *S);
NT_node main_compare_trees ( NT_node T1, NT_node T2, FILE *fp);
NT_node main_compare_trees_list ( NT_node T1, Sequence *S, FILE *fp);

float compare_trees ( NT_node T1, NT_node T2, int nseq, int mode);
float search_node ( NT_node B, NT_node T, int nseq, int mode);
float evaluate_node_similarity ( NT_node B, NT_node T, int nseq, int mode);

int compare_node ( int *b1, int *b2, int n);
void display_node (NT_node N, char *string,int nseq);
NT_node index_tree_node    (NT_node T);
NT_node simple_recode_tree (NT_node T, int nseq);
NT_node recode_tree ( NT_node T, Sequence *S);
int compare_branch2 ( int *b1, int *b2, int n);

/*********************************************************************/
/*                                                                   */
/*                                   FJ_tree Computation             */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
NT_node similarities_file2tree (char *mat);
NT_node tree_compute ( Alignment *A, int n, char ** arg_list);
static NT_node compute_std_tree (Alignment *A, int n, char **arg_list);
NT_node compute_std_tree_2 (Alignment *A, int **s, char *arg_list);
NT_node aln2fj_tree(NT_node T, Alignment *A, int limit,char* mode);
Alignment * filter_aln4tree (Alignment *A, int n,int fg,char* mode);

/*********************************************************************/
/*                                                                   */
/*                                   Tree Filters and MAnipulation   */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
int  tree2star_nodes (NT_node R, int n_max);
NT_node aln2tree (Alignment *A);
NT_node reset_boot_tree ( NT_node R, int n);
NT_node tree_dist2normalized_tree_dist ( NT_node R, float max);
NT_node reset_dist_tree ( NT_node R, float n);
NT_node* free_treelist  ( NT_node *R);
NT_node free_tree  ( NT_node R);
NT_node realloc_tree( NT_node R, int n);
NT_node free_tree_node ( NT_node R);

Sequence * tree2seq    (NT_node R, Sequence *S);
NT_node  rename_seq_in_tree ( NT_node R, char ***list);

NT_node balance_tree (NT_node);
int tree2nseq ( NT_node R);
int tree_file2nseq ( char *file);

int tree2nleaf ( NT_node R);
int tree2nnode ( NT_node R);
int tree2_nnode_unresolved (NT_node R, int *l);

FILE* display_tree ( NT_node R, int n, FILE *fp);
void clear_tree (NT_node T);
int display_leaf ( NT_node T, FILE *fp);
int display_leaf_below_node ( NT_node T, FILE *fp);
NT_node display_leaf_nb (NT_node T, int n, FILE *fp, char *name);
NT_node display_splits (NT_node T,Sequence *S, FILE *fp);
int tree2split_list (NT_node T, int nseq, int **split_list, int *n);

NT_node reroot_tree ( NT_node TREE, NT_node T);
NT_node straighten_tree ( NT_node P, NT_node C, float new_dist);
NT_node unroot_tree ( NT_node T);
FILE* print_tree_list ( NT_node *T,char *format, FILE *fp);
FILE* print_tree ( NT_node T,char *format, FILE *fp);
char *tree2string (NT_node T);
char *tree2file   (NT_node T, char *name, char *mode);

int print_newick_tree ( NT_node T, char *name);
FILE * rec_print_tree ( NT_node T, FILE *fp);


NT_node find_longest_branch ( NT_node T, NT_node L);
NT_node shift_root ( NT_node R);

int ** tree2cluster (NT_node T, float thres);
int ** make_sub_tree_list ( NT_node **T, int nseq, int n_node);
void make_all_sub_tree_list ( NT_node N, int **list, int *n);
void make_one_sub_tree_list ( NT_node T, int *list);
NT_node main_read_tree(char *treefile);

NT_node new_read_tree ( char *teefile);
NT_node new_get_node (NT_node T, FILE *fp);


NT_node** simple_read_tree(char *treefile);
void free_read_tree (NT_node **BT);
NT_node** read_tree(char *treefile, int *nnodes,int nseq, char **seq_names);
FILE * create_linear_tree ( char **name, int n, FILE *fp);
FILE * create_tree(NT_node ptree, NT_node parent,int *numseq,int  *ntotal,int  *nnodes,NT_node **lu, FILE *fp);
NT_node declare_tree_node (int nseq);
void set_info(NT_node p, NT_node parent, int pleaf, char *pname, float pdist, float bootstrap);
NT_node insert_tree_node(NT_node pptr);
FILE * skip_space(FILE *fd);
void create_tree_node(NT_node pptr, NT_node parent);
float calc_mean(NT_node nptr, float *maxdist, int nseq,NT_node **lu);
NT_node insert_root(NT_node p, float diff);
float calc_root_mean(NT_node root, float *maxdist, int neq, NT_node **lu);
NT_node reroot(NT_node ptree, int nseq, int ntotal, int nnodes, NT_node **lu);


Alignment *seq2seq_chain (Alignment *A,Alignment *B, char *arg);

float display_avg_bootstrap ( NT_node T);
float tree2tot_dist ( NT_node T, int mode);
int tree2n_branches(NT_node T, int mode);
int **display_tree_from_node (NT_node T, int up, int down, int **array);
NT_node tree2node ( char *name, NT_node T);
NT_node * tree2node_list (NT_node T, NT_node *L);
NT_node tree2root ( NT_node T);
int new_tree_sort ( char *name, NT_node T);


NT_node split2tree ( NT_node RT,Sequence *LIST, char *param);
NT_node * read_tree_list (Sequence *S);

int count_groups( Sequence *S, char *s);

Split ** count_splits( NT_node RT, Sequence *S, char *s);
NT_node *treelist2prune_treelist (Sequence *S, Sequence *TS, FILE *out);
int** treelist2groups (Sequence *S, Sequence *ST, char *depth, FILE *out);
int treelist2splits (Sequence *S, Sequence *ST);
int treelist2leafgroup ( Sequence *S, Sequence *TS, char *taxon);
int ***tree2dist ( NT_node T, Sequence *S, int ***d);
int treelist2frame (Sequence *S, Sequence *TS);
int** treelist2lti ( Sequence *S, Sequence *TS, int nb, FILE *out);

float simple_tree_cmp (NT_node T1, NT_node T2,Sequence *S, int mode);

int treelist2dmat ( Sequence *S);
NT_node new_declare_tree_node ();
int count_tree_groups( Sequence *LIST, char *group_file);
int node_sort ( char *name, NT_node T);
int    treelist2n (NT_node *L);
int ** treelist2avg_treecmp (NT_node *L, char *file);
NT_node treelist2bootstrap ( NT_node *L, char *file);
NT_node treelist2filtered_bootstrap ( NT_node *L, char *file, int **score,float f);

Sequence * treelist2seq ( Sequence *S);
Sequence * treelist2sub_seq ( Sequence *S, int f);

int treelist_file2consense (char *tree_file, char *outtree, char *outfile);
