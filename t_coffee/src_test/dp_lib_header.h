struct CL_node
    {
     
      int copy_mode;
      struct CL_node *c;
      struct CL_node *p;
      struct CL_node *l;
      struct CL_node *r;
      int seq;
      int res;
      int free;
};

typedef struct CL_node CL_node;
Alignment * add_constraint2aln ( Alignment *A, int s1, int r1, int s2, int r2);
Alignment * graph_aln (Alignment *A, Constraint_list *CL, Sequence *S);
Alignment* graph2aln (Alignment *A, CL_node *G, Sequence *S);
CL_node ***add_constraint2graph_aln (CL_node ***G, int s1, int r1, int s2, int r2);
CL_node * shift_segment ( CL_node *S, int segL, int shiftL);

int  is_graph_gap_column(CL_node *S);
CL_node * remove_graph_gap_column (CL_node *S);
CL_node * swap_gap_in_graph ( CL_node*S, CL_node *E);

CL_node * declare_cl_nodes ( int len, int seq);
CL_node * insert_gap_columns (CL_node *S, int d);
int get_node_distance ( CL_node *S, CL_node *E);
CL_node ***aln2graph (Alignment *A);
CL_node *vfree_graph (CL_node *S);
CL_node *vfree_cl_node ( CL_node *N);

int gotoh_pair_wise_lalign ( Alignment *A, int*ns, int **l_s,Constraint_list *CL);
Constraint_list * undefine_sw_aln ( Alignment *A, Constraint_list *CL);
Constraint_list * undefine_sw_pair ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int sw_pair_is_defined ( Constraint_list *CL, int s1, int r1, int s2, int r2);


int gotoh_pair_wise_sw (Alignment *A,int*ns, int **l_s,Constraint_list *CL);

Alignment * get_best_local_aln ( Alignment *IN,Constraint_list *CL,int gop, int gep, int sw_t, int sw_l, int sw_z, int greedy);
Alignment * get_best_nol_local_aln ( Alignment *IN, Constraint_list *CL, int gop, int gep,int sw_t,int sw_l, int sw_z, int mode);
double compute_penalty   (Constraint_list *CL, char *mode, int len);
double compute_scale ( Constraint_list *CL,char *mode, int len);
int evaluate_penalty (Alignment *A, Constraint_list *CL, int *scale,char *scale_mode, int *penalty, char *penalty_mode, int len_seq);
Alignment ** t_coffee_lalign   (Constraint_list *CL, int scale, int penalty,int maximise,Sequence *S, int sw_t, int sw_l, int sw_z,int *sw_n, int sw_io);
Alignment * add_seq2aln   (Constraint_list *CL, Alignment *IN,Sequence  *S);







struct Dp_Model
{
  int *diag;

  int TG_MODE;
  int F_TG_MODE;
  int gop;
  int gep;
  int f_gop;
  int f_gep;
  int nstate;
  int START;
  int END;
  
  char**model_comments;
  int **model;
  int **model_properties;
  int **bounded_model;
  int (***model_emission_function)(Alignment*, int **, int, int*, int, int **, int, int*, int, struct Constraint_list *);

  int LEN_I;
  int LEN_J;
  int DELTA_I;
  int DELTA_J;
  int EMISSION;
  int START_EMISSION;
  int TERM_EMISSION;
  
  int ALN_TYPE;
  Constraint_list *CL;
  /*Associated Functions*/
  
  /*To Deprecate*/
  int UM;
  
  int TYPE;
  int F0;
  int F1;
  
  
  int NON_CODING;
  int INSERTION;
  int DELETION;
  int CODING0;
  int CODING1;
  int CODING2;
  

};
typedef struct Dp_Model Dp_Model;

struct Dp_Result
{
  int *traceback;
  int len;
  int score;
  Dp_Model *Dp_model;
};
typedef struct Dp_Result Dp_Result;

Dp_Result * make_fast_generic_dp_pair_wise (Alignment *A, int*ns, int **l_s,Dp_Model *M);

Constraint_list* free_dp_model  (Dp_Model *D);
Dp_Result * free_dp_result (Dp_Result *D );

typedef struct hseq* SeqHasch;

typedef struct hseq
{
  SeqHasch hl[256];
  int  n;
  int  *l;
} hseq;

int ** ktup_dist_mat ( char **seq, int nseq, int ktup, char *type);
int ** evaluate_diagonals_for_two_sequences ( char *seq1, char *seq2,int maximise,Constraint_list *CL,int ktup);
int ** evaluate_diagonals ( Alignment *A, int *ns, int **l_s, Constraint_list *CL,int maximise,int n_groups, char **group_list,int ktup);
int ** evaluate_segments_with_ktup  ( Alignment *A, int *ns, int **l_s, Constraint_list *CL,int maximise,int n_groups, char **group_list,int ktup);
int ** evaluate_diagonals_with_ktup ( Alignment *A, int *ns, int **l_s, Constraint_list *CL,int maximise,int n_groups, char **group_list,int ktup);

int ** evaluate_diagonals_with_clist ( Alignment *A, int *ns, int **l_s, Constraint_list *CL,int maximise,int n_groups, char **group_list,int ktup);

int * flag_diagonals (int l1, int l2, int **sorted_diag,float T, int window);                
int * extract_N_diag (int l1, int l2, int **sorted_diag, int n_chosen_diag, int window);

int hasch_seq(char *seq1, int **hs, int **lu,int ktup, char *alph);
int fasta_gotoh_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL);
int cfasta_gotoh_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL);
int very_fast_gotoh_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL);

int make_fasta_gotoh_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL, int *diag);
/*********************************************************************/
/*                                                                   */
/*                         KTUP_DP                                   */
/*                                                                   */
/*                                                                   */
/*********************************************************************/

int precomputed_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL);
int ktup_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL);
int ktup_comparison ( char *seq1, char *seq2, int ktup);
HaschT* hasch_sequence ( char *seq1, int ktup);

SeqHasch * seq2hasch (int i,char *seq, int ktup, SeqHasch *H);
Constraint_list * hasch2constraint_list (Sequence*S, Constraint_list *CL);
SeqHasch *cleanhasch       (SeqHasch *H);
int hasch2sim        (SeqHasch *H, int nseq);
int cfasta_gotoh_pair_wise_sw (Alignment *A,int*ns, int **l_s,Constraint_list *CL);
int fasta_gotoh_pair_wise_sw (Alignment *A,int*ns, int **l_s,Constraint_list *CL);
int make_fasta_gotoh_pair_wise_sw (Alignment *A,int*ns, int **l_s,Constraint_list *CL, int *diag);

/*pair wise aln implementations*/

int idscore_pairseq (char *s1, char *s2, int gop, int gep, int **m, char *mode);
int idscore_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL);
int gotoh_pair_wise    (Alignment *A,int*ns, int **l_s,Constraint_list *CL);
int glocal_pair_wise ( Alignment *A, int *ns, int **l_s, Constraint_list *CL);
int gotoh_pair_wise_lgp ( Alignment *A, int *ns, int **l_s, Constraint_list *CL);
int test_pair_wise (Alignment *A, int *ns, int **l_s, Constraint_list *CL);

int glocal2_pair_wise (Alignment *A,int*ns, int **ls,Constraint_list *CL);
int linked_pair_wise ( Alignment *A, int *ns, int **l_s, Constraint_list *CL);
int linked_pair_wise_collapse ( Alignment *A, int *ns, int **l_s, Constraint_list *CL);
int cl2pair_list_ecl_ext_pc ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in);
void free_proba_pair_wise();

int subop1_pair_wise ( Alignment *A, int *ns, int **ls, Constraint_list *CL);
int subop2_pair_wise ( Alignment *A, int *ns, int **ls, Constraint_list *CL);
int proba_pair_wise ( Alignment *A, int *ns, int **ls, Constraint_list *CL);
int viterbi_pair_wise ( Alignment *A, int *ns, int **ls, Constraint_list *CL);
int biphasic_pair_wise ( Alignment *A, int *ns, int **l_s, Constraint_list *CL);

  

  
      
  
 

int cfasta_cdna_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL);
int fasta_cdna_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL);
Dp_Model* initialize_dna_dp_model (Constraint_list *CL);
Dp_Result * make_fast_dp_pair_wise (Alignment *A,int*ns, int **l_s, Constraint_list *CL,Dp_Model *M);
int make_fasta_cdna_pair_wise (Alignment *B,Alignment *A,int*ns, int **l_s,Constraint_list *CL, int *diag);



int ** evaluate_diagonals_cdna ( Alignment *A, int *ns, int **l_s, Constraint_list *CL,int maximise,int n_groups, char **group_list, int ktup);
int cfasta_cdna_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL);
Alignment *clean_maln ( Alignment *A, Alignment *I, int T, int n_it);
Alignment *realign_segment (int seq, int start, int len,Alignment *A, Alignment *C);
Dp_Model * initialize_seg2prf_model(int left_tg_mode, int right_tg_mode, Constraint_list *CL);

int get_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_start_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_term_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);

Dp_Model * initialize_sseq_model(int left_tg_mode, int right_tg_mode, Constraint_list *CL);
int ssec_pwaln_maln (Alignment *A, int *ns, int **ls, Constraint_list *CL);


int get_turn_gep_cost       (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_turn_start_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_turn_term_gep_cost  (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);

int get_alpha_gep_cost      (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_alpha_start_gep_cost(Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_alpha_term_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);

int get_beta_gep_cost       (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_beta_start_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_beta_term_gep_cost  (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);

int get_alpha_sub_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_beta_sub_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_turn_sub_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);

int get_ssec_no_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int sim_pair_wise_lalign (Alignment *in_A, int *in_ns, int **in_l_s, Constraint_list *in_CL);

/*pair wise aln implementations*/
int myers_miller_pair_wise (Alignment *A, int *ns, int **l_s,Constraint_list *CL);
int diff (Alignment *A, int *ns, int **ls, int s1, int M,int s2, int N , int tb, int te, Constraint_list *CL, int **pos);
int evaluate_est_order (Sequence *S, char *concat, Constraint_list *CL, int ktuple);
 
Constraint_list *prepare_cl_for_moca ( Constraint_list *CL);
Alignment ** moca_aln ( Constraint_list *CL);
Alignment * extract_domain ( Constraint_list *CL);
Alignment * interactive_domain_extraction ( Constraint_list *CL);
int print_moca_interactive_choices ();

Alignment * approximate_domain ( int min_start, int max_start, int step_start,int min_len, int max_len, int step_len, int *best_start, int *best_len, int *best_score, Constraint_list *CL);
 
int measure_domain_length ( Constraint_list *CL,Alignment *IN, int start, int min_size, int max_size, int step);
Alignment *extract_domain_with_coordinates ( Alignment *RESULT,int start, int len, Constraint_list *CL);
int get_starting_point ( Constraint_list *CL);

Alignment * find_domain_coordinates (Constraint_list *CL, int *start, int *len);
Alignment * extend_domain ( Constraint_list *CL, int *start, int *len, int dstart, int dlen);
Alignment * modify_domain ( Constraint_list *CL, Alignment *IN, int *start, int *len, int dstart, int dlen);

int * analyse_sequence ( Constraint_list *CL);

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
Constraint_list * profile_pair  ( TC_method*M,char *seq, Constraint_list *CL);
Constraint_list * thread_pair  ( TC_method*M,char *seq, Constraint_list *CL);
Constraint_list * thread_pair2 ( TC_method *M,int s1, int s2, Constraint_list *CL);
Constraint_list * sap_pair (char *seq, char *weight, Constraint_list *CL);
Constraint_list * lsqman_pair (char *seq, Constraint_list *CL);
Constraint_list * rna_pair (TC_method *M , char *in_seq, Constraint_list *CL);

/******************************************************************/
/*                   GENERIC PAIRWISE METHODS                     */
/*                                                                */
/*                                                                */
/******************************************************************/
Constraint_list *best_pair4prot (Job_TC *job);
Constraint_list *best_pair4rna (Job_TC *job);
Alignment * fast_pair      (Job_TC *job);
			    
void toggle_case_in_align_two_sequences(int value);
Alignment * align_two_sequences ( char *seq1, char *seq2, char *matrix, int gop, int gep, char *align_mode);
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

int empty_pair_wise ( Alignment *A, int *ns, int **l_s, Constraint_list *CL, int glocal);


Pwfunc get_pair_wise_function (Pwfunc func, char *dp_mode, int *glocal);



char *build_consensus ( char *seq1, char *seq2, char *dp_mode);
int fastal (int argv, char **arg);

int domain_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL );
Alignment *domain_match_list2aln ( Alignment *A,int *ns,int **l_s,int **ml, int nseq, int len);
Alignment * domain_seq2domain (Constraint_list *CL,int scale,int gop,int gep,Alignment *SEQ_DOMAIN, Alignment *TARGET);


int custom_pair_score_function1  (Constraint_list *CL, int s1, int r1, int s2, int r2);
int custom_pair_score_function2  (Constraint_list *CL, int s1, int r1, int s2, int r2);
int custom_pair_score_function3  (Constraint_list *CL, int s1, int r1, int s2, int r2);
int custom_pair_score_function4  (Constraint_list *CL, int s1, int r1, int s2, int r2);
int custom_pair_score_function5  (Constraint_list *CL, int s1, int r1, int s2, int r2);
int custom_pair_score_function6  (Constraint_list *CL, int s1, int r1, int s2, int r2);
int custom_pair_score_function7  (Constraint_list *CL, int s1, int r1, int s2, int r2);
int custom_pair_score_function8  (Constraint_list *CL, int s1, int r1, int s2, int r2);
int custom_pair_score_function9  (Constraint_list *CL, int s1, int r1, int s2, int r2);
int custom_pair_score_function10 (Constraint_list *CL, int s1, int r1, int s2, int r2);
int apdb (int argc, char *argv[]);

Constraint_list * set_constraint_list4align_pdb (Constraint_list *inCL,int seq, char *dp_mode, char *hasch_mode, char *param_file);
int evaluate_ca_trace_sap2_bubble (Constraint_list *CL, int s1, int r1, int s2, int r2);
int evaluate_ca_trace_nb (Constraint_list *CL, int s1, int s2, int r1, int r2);
int evaluate_ca_trace_bubble (Constraint_list *CL, int s1, int s2, int r1, int r2);
int evaluate_ca_trace_sap1_bubble (Constraint_list *CL, int s1, int s2, int r1, int r2);
int evaluate_ca_trace_3D_bubble (Constraint_list *CL, int s1, int s2, int r1, int r2);
int evaluate_ca_trace_transversal (Constraint_list *CL, int s1, int s2, int r1, int r2);
int evaluate_ca_trace_bubble_2 (Constraint_list *CL, int s1, int s2, int r1, int r2);
int evaluate_ca_trace_bubble_3 (Constraint_list *CL, int s1, int s2, int r1, int r2);


/*********************************************************************************************/
/*                                                                                           */
/*         FUNCTIONS FOR COMPARING TWO NEIGHBORHOODS:START                                   */
/*                                                                                           */
/*********************************************************************************************/
float matrix_match (Constraint_list *CL, int s1, int s2, int r1, int r2, Struct_nb *nbs1, Struct_nb *nbs2);
float transversal_match (Constraint_list *CL, int s1, int s2, int r1, int r2, Struct_nb *nbs1, Struct_nb *nbs2);
float neighborhood_match (Constraint_list *CL, int s1, int s2, int r1, int r2,Struct_nb *nbs1, Struct_nb *nbs2);
float sap1_neighborhood_match (Constraint_list *CL, int s1, int s2, int r1, int r2, Struct_nb *nbs1, Struct_nb *nbs2);
float sap2_neighborhood_match (Constraint_list *CL, int s1, int s2, int r1, int r2, Struct_nb *nbs1, Struct_nb *nbs2);


/*********************************************************************************************/
/*                                                                                           */
/*         FUNCTIONS FOR COMPARING TWO NEIGHBORHOODS:END                                     */
/*                                                                                           */
/*********************************************************************************************/
Alignment * analyse_pdb ( Alignment *A, Alignment *ST, char *filename);
Alignment * msa2struc_dist ( Alignment *A, Alignment *ST, char *filename, int gapped, int min_ncol);
float **** analyse_pdb_residues ( Alignment *A, Constraint_list *CL, Pdb_param *pdb_param);

float square_atom ( Atom *X);
Atom* reframe_atom ( Atom *X, Atom*Y, Atom *Z, Atom *IN, Atom *R);
Atom* add_atom ( Atom *A, Atom*B, Atom *R);
Atom* diff_atom ( Atom *A, Atom*B, Atom *R);
Atom * copy_atom ( Atom *A, Atom*R);
/************************************************************************/
/*                                                                      */
/*                            NUSSINOV                                  */
/*                                                                      */
/************************************************************************/
char *nussinov (char *S, int min_dist);
char * rna_struc2rna_lib ( char *seq_name, char *seq, char *name);
int display_rna_ss ( int pos, char *seq, char *struc, FILE *fp);
int evaluate_aln_gibbs   ( Alignment *A, Constraint_list *CL);
int evaluate_moca_domain ( Alignment *A, Constraint_list *CL);
int moca_residue_pair_extended_list ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int moca_evaluate_matrix_score      ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int moca_slow_get_dp_cost ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);

int **cache_cl_with_moca_domain (Alignment *A, Constraint_list *CL);
Alignment *make_moca_nol_aln ( Alignment *A, Constraint_list *CL);
/*********************************************************************************************/
/*                                                                                           */
/*         DOMAIN Z SCORE EVALUATION                                                         */
/*                                                                                           */
/*********************************************************************************************/

int evaluate_domain_aln_z_score (Alignment *A, int start, int end,Constraint_list *CL, char *alphabet);
int evaluate_domain_aln  ( Alignment *A, int start, int end,Constraint_list *CL);



int unpack_seq_residues ( int *s1, int *r1, int *s2, int *r2, int **packed_seq_lu);
Alignment * unpack_seq_aln ( Alignment *A,Constraint_list *C);
typedef struct 
{
  int N_COMPONENT;
  double *double_logB_alpha;
  double *exponant_list;
  double **ALPHA;
  double *DM_Q;
  double *alpha_tot;
  int n_aa;
  int tot_n;
}
Mixture;

double int_logB (int *i, int n);
double float_logB (float *i, int n);
double double_logB (double *x, int n);
double *** make_lup_table ( Mixture *D);
double  double_logB2(int j, double *n,Mixture *D);
double compute_exponant ( double *n, int j, Mixture *D);

double *compute_matrix_p ( double *n,int Nseq);
double* compute_dirichlet_p ( double *n,int Nseq);
void precompute_log_B ( double *table,Mixture *D);
double compute_X (double *n,int i,Mixture *D);
Mixture *read_dirichlet ( char *name);
int dirichlet_code( char aa);

double lgamma2 ( double x);
double lgamma_r(double x, int *signgamp);
/*********************************************************************************************/
/*                                                                                           */
/*         FUNCTIONS FOR EVALUATING THE CONSISTENCY BETWEEN ALN AND CL                       */
/*                                                                                           */
/*********************************************************************************************/
Alignment * overlay_alignment_evaluation     ( Alignment *I, Alignment *O);
Alignment * main_coffee_evaluate_output          ( Alignment *IN,Constraint_list *CL, const char *mode );
int aln2ecl_raw_score (Alignment *A, Constraint_list *C);
int sub_aln2ecl_raw_score (Alignment *A, Constraint_list *CL, int ns, int *ls);
int sub_aln2sub_aln_raw_score ( Alignment *IN,Constraint_list *CL, const char *mode, int *ns, int **ls);
int node2sub_aln_score(Alignment *A,Constraint_list *CL, char *mode, NT_node T);
int sub_aln2sub_aln_score ( Alignment *IN,Constraint_list *CL, const char *mode, int *ns, int **ls);
Alignment *  main_coffee_evaluate_output_sub_aln ( Alignment *IN,Constraint_list *CL, const char *mode, int *ns, int **ls);

Alignment * categories_evaluate_output               ( Alignment *IN,Constraint_list *CL);
Alignment * matrix_evaluate_output               ( Alignment *IN,Constraint_list *CL);
Alignment * sar_evaluate_output ( Alignment *IN,Constraint_list *CL);
Alignment * boxshade_evaluate_output ( Alignment *IN,Constraint_list *CL, int T);

Alignment * fast_coffee_evaluate_output          ( Alignment *IN,Constraint_list *CL);

int slow_coffee_evaluate_sub_aln         ( Alignment *IN,int *ns, int **ls, Constraint_list *CL);
Alignment * slow_coffee_evaluate_output          ( Alignment *IN,Constraint_list *CL);
Alignment * non_extended_t_coffee_evaluate_output( Alignment *IN,Constraint_list *CL);
Alignment * heuristic_coffee_evaluate_output     ( Alignment *IN,Constraint_list *CL);
Alignment *coffee_seq_evaluate_output ( Alignment *IN, Constraint_list *CL);
/*Old Function: To deprecate*/
Alignment * coffee_evaluate_output ( Alignment *IN,Constraint_list *CL);

/*********************************************************************************************/
/*                                                                                           */
/*        PROFILE/PRofile Functions                                                          */
/*                                                                                           */
/*********************************************************************************************/
Profile_cost_func get_profile_mode_function (char *name, Profile_cost_func func);
int generic_evaluate_profile_score     (Constraint_list *CL,Alignment *Prf1, int s1, int r1, Alignment *Prf2, int s2, int r2, Profile_cost_func prf_prf);
int cw_profile_profile     (int *prf1, int *prf2, Constraint_list *CL);
int muscle_profile_profile     (int *prf1, int *prf2, Constraint_list *CL);

/*********************************************************************************************/
/*                                                                                           */
/*         FUNCTIONS FOR GETING THE COST : (Sequences) ->evaluate_residue_pair               */
/*                                                                                           */
/*********************************************************************************************/
int initialize_scoring_scheme (Constraint_list *CL);
int evaluate_blast_profile_score (Constraint_list *CL, int s1, int r1, int s2, int r2);
int evaluate_aln_profile_score   (Constraint_list *CL, int s1, int r1, int s2, int r2);

int evaluate_profile_score                         ( Constraint_list *CL,Alignment *Prf1, int s1, int r1, Alignment *Prf2, int s2, int r2);
int evaluate_cdna_matrix_score                     ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int evaluate_diaa_matrix_score                     ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int evaluate_monoaa_matrix_score                   ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int evaluate_matrix_score                          ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int evaluate_tm_matrix_score                          ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int evaluate_ssp_matrix_score                          ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int evaluate_curvature_score                          ( Constraint_list *CL, int s1, int r1, int s2, int r2);

int evaluate_combined_matrix_score                 ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int evaluate_physico_score                         ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int residue_pair_non_extended_list                 ( Constraint_list *CL, int s1, int r1, int s2, int r2);

int residue_pair_extended_list4rna1                 ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int residue_pair_extended_list4rna2                 ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int residue_pair_extended_list4rna3                 ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int residue_pair_extended_list4rna4                 ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int residue_pair_extended_list4rna                  ( Constraint_list *CL, Constraint_list *R, int s1, int r1, int s2, int r2);

int residue_pair_extended_list_raw                     ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int residue_pair_extended_list_pc                     ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int residue_pair_extended_list_4gp                     ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int residue_pair_extended_list                     ( Constraint_list *CL, int s1, int r1, int s2, int r2);




int residue_pair_extended_list_g_coffee_quadruplet ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int residue_pair_extended_list_g_coffee            ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int residue_pair_extended_list_quadruplet          ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int residue_pair_extended_list_mixt                ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int residue_pair_test_function                     ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int extend_residue_pair                            ( Constraint_list *CL, int s1, int r1, int s2, int r2);

int residue_pair_relative_extended_list ( Constraint_list *CL, int s1, int r1, int s2, int r2 );
/*********************************************************************************************/
/*                                                                                           */
/*         FUNCTIONS FOR GETTING THE PW COST :  CL->get_dp_cost                                     */
/*                                                                                           */
/*********************************************************************************************/
int get_dp_cost_blosum_matrix (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_dp_cost_pam_matrix (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_dp_cost_pw_matrix (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);

int get_cdna_best_frame_dp_cost ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_dp_cost                 ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_dp_cost_quadruplet      ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
      

int very_fast_get_dp_cost       ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);

int cw_profile_get_dp_cost ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int cw_profile_get_dp_cost_window ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);

int consensus_get_dp_cost ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);

int fast_get_dp_cost            ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int fast_get_dp_cost_2          ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int fast_get_dp_cost_3          ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);

int fast_get_dp_cost_quadruplet ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);

int check_fast_profile_mode    (Alignment *A, int ns1,int *list1,int ns2, int *list2, Constraint_list *CL);
int check_fast_mode    (Alignment *A, int ns1,int *list1,int ns2, int *list2, Constraint_list *CL);


int slow_get_dp_cost            ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int slow_get_dp_cost_pc         ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int slow_get_dp_cost_test       ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);

int sw_get_dp_cost              ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_domain_dp_cost          ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2,Constraint_list *CL , int scale , int gop, int gep);

/*********************************************************************************************/
/*                                                                                           */
/*         FUNCTIONS FOR ANALYSING AL OR MATRIX                                              */
/*                                                                                           */
/*********************************************************************************************/
float ** initialise_aa_physico_chemical_property_table (int *n);
int aln2n_res ( Alignment *A, int start, int end);
float get_gop_scaling_factor ( int **matrix,float id, int l1, int l2);
float get_avg_matrix_mm ( int **matrix, char *alphabet);
float get_avg_matrix_match ( int **matrix, char *alphabet);
float get_avg_matrix_diff ( int **matrix1,int **matrix2, char *alphabet);
float measure_matrix_enthropy (int **matrix,char *alphabet);
float measure_matrix_pos_avg (int **matrix,char *alphabet);
float evaluate_random_match (char *matrix, int n, int len,char *alp);
float compare_two_mat (char  *mat1,char*mat2, int n, int len,char *alp);
float compare_two_mat_array (int  **matrix1,int **matrix2, int n, int len,char *alp);
int ** rescale_two_mat (char  *mat1,char*mat2, int n, int len,char *alp);
int ** rescale_matrix ( int **mat, float lambda, char *alp);
int **mat2inverted_mat (int **matrix, char *alp);
void output_matrix_header ( char *name, int **matrix, char *alp);
float evaluate_random_match2 (int **matrix, int n, int len,char *alp);
float measure_lambda2(char  *mat1,char*mat2, int n, int len,char *alp);
float measure_lambda(char  *mat1,char*mat2, int n, int len,char *alp);
Constraint_list * choose_extension_mode ( char *extend_mode, Constraint_list *CL);
int ** combine_two_matrices ( int **mat1, int **mat2);
/*********************************************************************************************/
/*                                                                                           */
/*         OFF THE SHELVES EVALUATION                                              */
/*                                                                                           */
/*********************************************************************************************/
float  sum_pair ( Alignment*A,char *mat_name, int gap_op, int gap_ext);
int lat_sum_pair (Alignment *A, char *mat);
int ** show_pair(int istart, int iend, int jstart, int jend, int *seqlen_array, char **seq_array, int dna_ktup, int dna_window, int dna_wind_gap, int dna_signif,int prot_ktup, int prot_window,int prot_wind_gap,int prot_signif, int nseqs,int dnaflag, int max_aa, int max_aln_length  );
/*By convention, 0: START, 1: END*/

struct Hmm
{
  
  double freeT;           /*Free transition*/
  double forbiden;        /*Forbiden transition*/
  int start;              /*start, by convention: 0*/
  int end;                /*end, by convention: 1*/

  int nS;                 /*Number of states*/
  int order;
  struct HmmState  **S;   /*State List*/

  /*Bounded HMM*/
  double **T;             /*Transition matrix*/
  int **fromM;            /*For any sate s, fromM[0]->number of states leading to s*/
  int **toM;              /*For any sate s, toM[0]->number of s can go to*/
                          /*End and Start are NOT included in toM and FromM*/
  
  
};
typedef struct Hmm Hmm;

struct HmmAln
{
  Hmm *H;
  int *state_list;
};
typedef struct HmmAln HmmAln;

typedef double (*Generic_em_func)(struct Hmm*, struct HmmState*, int);

struct HmmState
{
char name[100];
int state;
int DJ;
int DI;

  /*Pair HMM*/
double em;
Col_cost_func em_func;

  /*Linear HMM*/
double *em2;
Generic_em_func em_func2;
int nT;
struct StateTrans **T;
};
typedef struct HmmState HmmState;

struct StateTrans
{
  char name[101];
  double tr;
};
typedef struct StateTrans StateTrans;

struct MatState
{
  int i;
  int j;
  int st;
  int pst;
  double sc;
  struct MatState *n;
  struct MatState *p;
  struct MatState *m; /*memory*/
  struct MatState *s; /*memory of the start point*/
  /*Heap Mamagement: Never copy*/
  struct MatState *Hn;/*Heap Next*/ 
  struct MatState *Hp;/*Heap Previous*/ 
  
  struct MatState *Mn;/*Marked Heap section*/ 
  struct MatState *Mp;/*Marked Heap Section*/ 
  int free;
};
typedef struct MatState MatState;


/*********************************************************************************/
/*                                                                               */
/*                                                                               */
/*                     simple HMM: Viterbi                                       */
/*                                                                               */
/*                                                                               */
/*********************************************************************************/
int seq_viterbi_pair_wise (Alignment *A,int*ns, int **ls,Constraint_list *CL);

/*********************************************************************************/
/*                                                                               */
/*                                                                               */
/*                     HMM: Viterbi                                              */
/*                                                                               */
/*                                                                               */
/*********************************************************************************/

int viterbi_pair_wise_OLD (Alignment *A,int*ns, int **ls,Constraint_list *CL);
Alignment * viterbipath2aln (Alignment *A, int *ns,int **ls,int *tb, Hmm *H);
double*** viterbi_hmm (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL);
int * viterbi2path (int l1,int l2, Hmm *H, double ***M);
/*********************************************************************************/
/*                                                                               */
/*                                                                               */
/*                      HMM modeling                                             */
/*                                                                               */
/*                                                                               */
/*********************************************************************************/
int viterbiL_pair_wise (Alignment *A,int*ns, int **ls,Constraint_list *CL);
MatState* RviterbiL_hmm (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL,MatState *S, MatState *E);
MatState* viterbiL_hmm (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL, MatState *S, MatState *E);

int viterbiD_pair_wise (Alignment *A,int*ns, int **ls,Constraint_list *CL);
double lu_RviterbiD_hmm (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL,MatState *S, MatState *E, int **seg_list);
MatState* RviterbiD_hmm (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL,MatState *S, MatState *E, int **seg_list);
MatState* viterbiD_hmm (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL, MatState *S, MatState *E, int **seg_list);
int **seglist2table ( int **seglist,int l1, int l2);
int *seglist2line ( int **seglist, int line, int start, int end);
int * traceback (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL,MatState *S, MatState *E, int **seg_list);

int viterbiDGL_pair_wise (Alignment *A,int*ns, int **ls,Constraint_list *CL);
double lu_RviterbiDGL_hmm (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL,MatState *S, MatState *E, int **seg_list);
MatState* RviterbiDGL_hmm (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL,MatState *S, MatState *E, int **seg_list);
MatState* viterbiDGL_hmm (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL,MatState *S, MatState *E, int **seg_list);


/*********************************************************************************/
/*                                                                               */
/*                                                                               */
/*                      HMM modeling                                             */
/*                                                                               */
/*                                                                               */
/*********************************************************************************/
int MatStateAreIdentical (MatState*I, MatState*O);
int MaxDeltaMatState (MatState*I, MatState*O);
int MinDeltaMatState (MatState*I, MatState*O);

MatState * ManageMatState(int Mode, MatState *C);
MatState* CopyMatState ( MatState*I, MatState*O);

Hmm* read_hmm(char *file);
Hmm* declare_hmm(int n);
Hmm* free_Hmm(Hmm*H);
void DisplayHmm ( Hmm *H);

/*********************************************************************************/
/*                                                                               */
/*                                                                               */
/*                      HMM Models                                               */
/*                                                                               */
/*                                                                               */
/*********************************************************************************/
Hmm* define_two_mat_model(Constraint_list *CL);
Hmm* define_probcons_model(Constraint_list *CL);
Hmm* define_simple_model(Constraint_list *CL);

Hmm * bound_hmm ( Hmm *H);
Sequence *pavie_seq2random_seq (Sequence *S, char *subst); 
double **pavie_seq2pavie_aln(Sequence *S,char *mat, char *mode);
Alignment *pavie_seq2pavie_sort ( Sequence *S, char *mat, char *mode);
Alignment* pavie_seq2pavie_msa  ( Sequence *S, char *mat, char *mode);
NT_node pavie_seq2pavie_tree ( Sequence *S, char *mat, char *mode);
int **pavie_seq2trained_pavie_mat(Sequence *S, char *param);
int pavie_pair_wise (Alignment *A,int *ns, int **l_s,Constraint_list *CL );
Sequence *pavie_seq2noisy_seq (Sequence *S, int val,char *subst); 
