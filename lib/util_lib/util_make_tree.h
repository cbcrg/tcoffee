char* seq2aln_file (char *pg,char *in, char *out);
char* seq2cw_aln_file (char *in, char *out);
char* seq2co_aln_file (char *in, char *out);
char* seq2mafft_aln_file (char *in, char *out);

char* prf_pair2cw_aln_file (char *prf1,char *prf2, char *out);


NT_node ** make_nj_tree (  Alignment *A,int **distances,int gop, int gep, char **out_seq, char **out_seq_name, int out_nseq, char *tree_file, char *tree_mode);
NT_node ** make_upgma_tree (  Alignment *A,int **distances,int gop, int gep, char **out_seq, char **out_seq_name, int out_nseq, char *tree_file, char *tree_mode);


NT_node ** int_dist2nj_tree (int **distances, char **out_seq_name, int out_nseq,  char *tree_file);
NT_node ** float_dist2nj_tree (float **distances, char **out_seq_name, int out_nseq,  char *tree_file);
NT_node ** dist2nj_tree (double **distances, char **out_seq_name, int out_nseq,  char *tree_file);

//////////////////////////////////////////////////////////////////////////////
//
//                              km
///////////////////////////////////////////////////////////////////////////////
double **vector2strip_vector (double**v, int n, int *dim, float frac);
double **aln2km_vector (Alignment *A, char *mode, int *dim);

NT_node   seq2dnd (Sequence *S, char *mode);
NT_node **seq2km_tree_old (Sequence *S, char *file);
NT_node   seq2cat_dnd (Sequence *S, char *mode);
NT_node   seq2swl_dnd (Sequence *S);
NT_node   seq2km_dnd (Sequence *S);
NT_node   seq2co_dnd (Sequence *S);
NT_node   list2balanced_dnd (char **name, int n);
NT_node   seq2cw_dnd ( Sequence *S);
NT_node   seq2cwquick_dnd ( Sequence *S);

NT_node   seq2parttree_dnd ( Sequence *S);
NT_node   seq2mafft_dnd ( Sequence *S);
NT_node   seq2fftns1_dnd ( Sequence *S);
NT_node   seq2fftns2_dnd ( Sequence *S);


NT_node   seq2dpparttree_dnd ( Sequence *S);
NT_node   seq2fastparttree_dnd ( Sequence *S);
NT_node   seq2cw_tree ( Sequence *S);

//functions for the estimation of a reg_tree
NT_node   seq2reg_tree( Sequence *S);
NT_node   addseq2reg_tree (NT_node T, Sequence *S, int seq, int depth);
float* node2reg_score(NT_node T, Sequence *S, char *s1, float *v, int depth);

NT_node compute_cw_tree (Alignment *A);
NT_node aln2cw_tree     (Alignment *A);
NT_node aln2km_tree (Alignment *A, char *mode, int nboot);
NT_node rec_km_tree (char **name,int n,int dim,double **V, int nboot);
  

NT_node ** int_dist2upgma_tree (int **mat, Alignment *A, int nseq, char *fname);

NT_node int_dist2upgma_tree_new (int **mat, char **name, int nseq);
NT_node **dist2upgma_tree (double **mat, char **name, int nseq, char *fname);
NT_node upgma_merge (int **mat, NT_node *NL, int *used, int *n, int N);

void nj_tree(char **tree_description, int nseq);
void fast_nj_tree(char **tree_description);
void slow_nj_tree(char **tree_description);

void print_phylip_tree(char **tree_description, FILE *tree, int bootstrap);
void two_way_split(char **tree_description, FILE *tree, int start_row, int flag, int bootstrap);
void guide_tree(char *fname, double **saga_tmat, char **sag_seq_name, int saga_nseq);



NT_node split2upgma_tree (Split **S, Alignment *A, int nseq, char *fname);
NT_node split_upgma_merge (Alignment *A, Split **S, NT_node *NL, int *used, int *n, int N);
float get_split_dist ( Alignment *A, NT_node L, NT_node R, Split **S) ;

Alignment * upgma_tree_aln  (Alignment*A, int nseq, Constraint_list *CL);
int ** dist_mat2best_split (int **mat, int nseq);
int upgma_node_heap (NT_node X);
