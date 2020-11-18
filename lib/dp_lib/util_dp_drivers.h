
/****************************************************************************/
/*                                                                          */	
/*                                                                          */	
/*            Alignment Methods                                             */
/*                                                                          */	
/*                                                                          */	
/****************************************************************************/
Alignment * sorted_aln (Alignment *A, Constraint_list *CL);
Alignment * sorted_aln_seq (int seq, Alignment *A, Constraint_list *CL);
Alignment * full_sorted_aln (Alignment *A, Constraint_list *CL);

/******************************************************************/
/*                   MAIN DRIVER                                  */
/*                                                                */
/*                                                                */
/******************************************************************/
int *  set_profile_master (Alignment *A,int *ns, int **ls, Constraint_list *CL);
int *unset_profile_master (Alignment *A,int *ns, int **ls, Constraint_list *CL);

Constraint_list *profile2list     ( Job_TC *job,int nprf);
Constraint_list *seq2list     (Job_TC *Job);
Constraint_list *method2pw_cl (TC_method *M, Constraint_list *CL);
int method_uses_structure(TC_method *M);
int method_uses_profile(TC_method *M);

/******************************************************************/
/*                   MULTIPLE ALIGNMENTS                          */
/*                                                                */
/*                                                                */
/******************************************************************/
Alignment * compute_prrp_aln (Alignment *A, Constraint_list *CL);
Alignment * compute_clustalw_aln (Alignment *A, Constraint_list *CL);
Alignment * compute_tcoffee_aln_quick (Alignment *A, Constraint_list *CL);
Alignment * seq2clustalw_aln (Sequence *S);
Alignment * aln2clustalw_aln (Alignment *A, Constraint_list *CL);
Alignment * realign_block ( Alignment *A, int col1, int col2, char *pg);
char *seq_file2msa_file (char *file, char *seq, char *aln);
int seq_are_duplicated (char *seq);
Alignment *seq2msa (Sequence *S);

/******************************************************************/
/*                  DNA                                           */
/*                                                                */
/*                                                                */
/******************************************************************/
Constraint_list * align_coding_nucleotides (char *seq, char *method, char *weight, char *mem_mode, Constraint_list *CL);
/******************************************************************/
/*                   STRUCTURES                                   */
/*                                                                */
/*                                                                */
/******************************************************************/
Constraint_list * seq_msa (TC_method *M , char *in_seq, Constraint_list *CL);

Constraint_list *align_pdb_pair   (char *seq_in, char *dp_mode,char *evaluate_mode, char *file, Constraint_list *CL, Job_TC *job);
Constraint_list * align_pdb_pair_2 (char *seq, Constraint_list *CL);

Constraint_list * pdb_pair  ( TC_method*M,char *seq, Constraint_list *CL);
Constraint_list * pdbid_pair  ( TC_method*M,char *seq, Constraint_list *CL);
Constraint_list * profile_subset_pair  ( TC_method*M,char *seq, Constraint_list *CL);
Constraint_list * profile_pair             ( TC_method*M,char *seq, Constraint_list *CL);
Constraint_list * profile_pair_decomposed  ( TC_method*M,char *seq, Constraint_list *CL, char *mode);

Constraint_list * thread_pair  ( TC_method*M,char *seq, Constraint_list *CL);
Constraint_list * thread_pair2 ( TC_method *M,int s1, int s2, Constraint_list *CL);
Constraint_list * sap_pair (char *seq, char *weight, Constraint_list *CL);
Constraint_list * lsqman_pair (char *seq, Constraint_list *CL);
Constraint_list * rnapdb_pair (TC_method *M , char *in_seq, Constraint_list *CL);

/******************************************************************/
/*                   GENERIC PAIRWISE METHODS                     */
/*                                                                */
/*                                                                */
/******************************************************************/
Constraint_list *best_pair4prot (Job_TC *job);
Constraint_list *best_pair4rna (Job_TC *job);
Alignment * fast_pair      (Job_TC *job);
Constraint_list *hash_pair (TC_method *M , char *in_seq, Constraint_list *CL);

void toggle_case_in_align_two_sequences(int value);
Alignment * align_two_structures ( Sequence *S, int s1, int s2, char *mode);
Alignment * align_two_sequences ( char *seq1, char *seq2, char *matrix, int gop, int gep, char *align_mode);
Alignment *align_two_streches4dpa ( char *s0, char *s1, char *in_matrix, int gop, int gep, char *in_align_mode, Alignment *A);
Alignment * align_two_sequences4dpa ( char *padded1, char *seq1, char *padded2, char *seq2,char *matrix, int gop, int gep, char *align_mode, Alignment *R);
Alignment * align_two_aln ( Alignment *A1, Alignment  *A2, char *in_matrix, int gop, int gep, char *in_align_mode);
NT_node make_root_tree ( Alignment *A,Constraint_list *CL,int gop, int gep,Sequence *S,  char *tree_file,int maximise);
NT_node ** make_tree ( Alignment *A,Constraint_list *CL,int gop, int gep,Sequence *S,  char *tree_file, int maximise);
int ** get_pw_distances ( Alignment *A,Constraint_list *CL,int gop, int gep, char **out_seq, char **out_seq_name, int out_nseq, char *tree_file, char *tree_mode, int maximise);
Alignment *stack_progressive_nol_aln_with_seq_coor(Constraint_list *CL,int gop, int gep,Sequence *S, int **seq_coor, int nseq);
Alignment *stack_progressive_aln_with_seq_coor (Constraint_list*CL,int gop, int gep, Sequence *S, int **coor, int nseq);
Alignment *stack_progressive_aln(Alignment *A, Constraint_list *CL, int gop, int gep);
Alignment *est_progressive_aln(Alignment *A, Constraint_list *CL, int gop, int gep);
void analyse_seq ( Alignment *A, int s);

char ** list_file2dpa_list_file (char **list_file, int *len,int maxnseq, Sequence *S);
Alignment * seq2aln_group (Alignment *A, int T, Constraint_list *CL);

Alignment *profile_aln (Alignment *A, Constraint_list *CL);
Alignment * iterative_tree_aln (Alignment *A,int n, Constraint_list *CL);
Alignment * iterative_aln ( Alignment*A, int nseq, Constraint_list *CL);
Alignment * seq_aln ( Alignment*A, int nseq, Constraint_list *CL);
Alignment *tsp_aln (Alignment *A, Constraint_list *iCL, Sequence *S);
Alignment *iterate_aln ( Alignment*A, int nit, Constraint_list *CL);
Alignment *realign_aln ( Alignment*A, Constraint_list *CL);
Alignment *very_fast_aln (Alignment*A, int nseq, Constraint_list *CL);
Alignment *simple_progressive_aln (Sequence *S, NT_node **T, Constraint_list *CL, char *mat);
Alignment *frame_aln (Alignment *A, int n,Constraint_list *CL);
Alignment *dpa_aln (Alignment *A, Constraint_list *CL);
Alignment *new_dpa_aln (Alignment *A, Constraint_list *CL);
Alignment * make_delayed_tree_aln (Alignment *A,int n, Constraint_list *CL);

NT_node* delayed_tree_aln ( NT_node LT, NT_node RT, Alignment*A, int nseq, Constraint_list *CL);
int  node2seq_list (NT_node P, int *ns, int *ls);
Alignment* delayed_tree_aln1 ( NT_node P,Alignment*A,Constraint_list *CL, int threshold);
Alignment* delayed_tree_aln2 ( NT_node P,Alignment*A,Constraint_list *CL, int threshold);

NT_node* tree2ao (NT_node LT, NT_node RT,Alignment *A, int nseq,Constraint_list *CL);//tree2align_order
NT_node* tree_aln ( NT_node LT, NT_node RT, Alignment*A, int nseq, Constraint_list *CL);
NT_node* local_tree_aln ( NT_node LT, NT_node RT, Alignment*A, int nseq, Constraint_list *CL);
NT_node* seqan_tree_aln ( NT_node LT, NT_node RT, Alignment*A, int nseq, Constraint_list *CL);


NT_node* tree_realn ( NT_node LT, NT_node RT, Alignment*A, int nseq, Constraint_list *CL);

int split_condition (int nseq, int score, Constraint_list *CL);

int profile_pair_wise (Alignment *A, int n1, int *l1, int n2, int *l2, Constraint_list *CL);
int pair_wise (Alignment *A, int*ns, int **l_s,Constraint_list *CL );
int hh_pair_wise (Alignment *A, int *ns, int **ls, Constraint_list *CL);
int co_pair_wise (Alignment *A, int *ns, int **ls, Constraint_list *CL);
int empty_pair_wise ( Alignment *A, int *ns, int **l_s, Constraint_list *CL, int glocal);


Pwfunc get_pair_wise_function (Pwfunc func, char *dp_mode, int *glocal);


/*********************************************************************/
char *build_consensus ( char *seq1, char *seq2, char *dp_mode);
int fastal (int argv, char **arg);
int tree_aln_N ( NT_node P, Sequence *S, int N, int argc, char **argv);
int updown_tree_aln (NT_node T, Sequence *S, int max, int *n, int argc, char **argv);
