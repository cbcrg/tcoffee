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
Alignment * triplet_coffee_evaluate_output          ( Alignment *IN,Constraint_list *CL);
Alignment * fast_coffee_evaluate_output          ( Alignment *IN,Constraint_list *CL);
Alignment * sp3_evaluate4tcoffee (Alignment *RNA, Constraint_list *CLin);
Alignment * distance_evaluate4tcoffee (Alignment *A, Constraint_list *CL, float max, float delta, int enb);
Alignment * strike_evaluate4tcoffee (Alignment *A, Constraint_list *CL,char *matrix);
Alignment * struc_evaluate4tcoffee (Alignment *A, Constraint_list *CL,char *mode,float max, int enb,char *in_matrix_name);
Alignment * struc_evaluate4tcoffee4gt (Alignment *A, Constraint_list *CL,char *mode,float max, int enb,char *in_matrix_name);
Alignment * msa2distances (Alignment *A, Constraint_list *CL, float radius, float threshold, int mon);
Alignment * msa_list2struc_evaluate4tcoffee (Sequence *Species, Sequence *MSAs, Sequence *CL,char *mode,float max, int enb,char *in_matrix_name);


char *shuffle_seq_file(char *file);
int realign_node4hot (Sequence *S, NT_node T,char *pg, int shuffle, char *name, int n);
int hot (Sequence *S,NT_node T, char *pg, int shuffle, char *name, int n);
float shuff (Sequence *S,char *pg,int n);
float hotshot (Sequence *S,NT_node T, char *pg, float *tot, float *n);
float* hotshot2 (Sequence *S,NT_node T, char *pg, float *results);
float realign_node4hotshot (Sequence *S, NT_node T,char *pg, int shuffle);
float* realign_node4hotshot2 (Sequence *S, NT_node T,char *pg,float *results);
Alignment * hotnode2aln  (Sequence *S, NT_node T,char *pg, int shuffle, char *flip);

Alignment * evaluate_tree_group (Alignment *T, Sequence *G);
Alignment *treealn_evaluate4tcoffee (Alignment *A, Sequence *G);

int  sp_triplet_coffee_evaluate_output2  ( Alignment *IN,Constraint_list *CL, char *fname);
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
int get_dp_cost_sankoff_tree  (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_dp_cost_pam_matrix (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_dp_cost_pw_matrix (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);

int get_cdna_best_frame_dp_cost ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int get_dp_cost4dpa             ( Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
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
