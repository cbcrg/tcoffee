typedef struct tnode *NT_node;

/**
* Node of a tree
*/
typedef struct aanode *LAA;
typedef struct aanode{
  int aa;
  LAA bef;
  LAA col;
}Aanode;



typedef struct
    {
      int   p1;
      int   p2;
      int   p3;
      int   p4;
      int   t;
      int   f;
      char  mode[20];//lower, unalign
      char  model[20];//fsa1 fsa2
}
OveralnP;

//RNA
Constraint_list * seq2contacts (Sequence *S, Sequence *T,Constraint_list *CL, char *mode);
Constraint_list * pdb2contacts (Sequence *S, Sequence *T,Constraint_list *CL, char *mode, char *type,float maxD);
int vienna2template_file (char *outfile, Sequence *R, Sequence *ST);
Constraint_list * vienna2tc_lib (char *outfile, Sequence *R, Sequence *ST);


int ** alifold_list2cov_list (Alignment *A, int **list);
int ** update_RNAfold_list (Alignment *A, int **pos, int s, int **l);
int ** vienna2list ( char *seq);
Alignment *compare_RNA_fold ( Alignment *A, Alignment *B);
Alignment *trim_RNA (Alignment *RNA, Sequence *ST, int max);
Alignment *sp3_evaluate (Alignment *A, Sequence *S);
Alignment *alifold2analyze (Alignment *A, Alignment *ST, char *mode);
Alignment *alifold2cov_aln (Alignment *A, int **l, int ug);
Alignment *alifold2cov_stat (Alignment *A, int **l, int ug);
Alignment *alifold2cov_list (Alignment *A, int **l, int ug);
Alignment *alifold2cov_cache (Alignment *inA,  int **l, int ug);


Alignment *add_alifold2aln  (Alignment *A, Alignment *ST);
Alignment *aln2alifold(Alignment *A);

//end
Alignment * aln2bootstrap (Alignment *A, int n);
Alignment * aln2sample    (Alignment *A, int n);
Alignment * aln2random_aln(Alignment *A, char *mode);
Alignment * aln2scale     (Alignment *A, char *offset);
Alignment * aln2case_aln  (Alignment *A, char *upper, char *lower);
Alignment * aln2gap_cache (Alignment *A, int val);
Alignment * score_aln2score_ascii_aln (Alignment *A, Alignment *C);
int **aln2resindex ( Alignment *A, Alignment *B, FILE *fp);
int **index_seq_res      ( Sequence *S1, Sequence *S2, int **name_index);
int **index_seq_name ( Sequence *S1, Sequence *S2);
int *get_name_index (char **l1, int n1, char **l2, int n2);

int* get_res_index (char *seq1, char *seq2);
int * pos2list (int * pos, int len, int *nl);
  int *list2pos (int *list, int nl, int len);


int change_residue_coordinate ( char *in_seq1, char *in_seq2, int v);

int ** minimise_repeat_coor (int **coor, int nseq, Sequence *S);
int ** get_nol_seq( Constraint_list *CL,int **coor, int nseq, Sequence *S);


int compare_pos_column( int **pos1,int p1, int **pos2,int p2, int nseq);
char  * seq2alphabet (Sequence *S);

double *seq2swv   (char *seq, char **seql, int n);
double *seq2swr   (char *seq, char **seql, int n);
double* seq2diaa (char *buf, double *v);
double* seq2triaa(char *buf, double *v);
double* seq2tetraa(char *buf, double *v);

char *aln2alphabet (Alignment *A);
char *array2alphabet (char **array, int n, char *forbiden);






//TM Predictions
char* alnpos2hmmtop_pred (Alignment *A, Alignment *Pred, int pos, int mode);
Alignment * aln2hmmtop_pred (Alignment *A);
char * seq2tmstruc ( char *seq);

void set_blast_default_values();
char      *  seq2pdb   ( Sequence *S);
Sequence *  seq2prf ( Sequence *S);
char* clean_sname(char sname);
Sequence *  seq2blast ( Sequence *S);
Sequence *  seq2blast_thread ( Sequence *S);


Sequence * seq2unique_name_seq ( Sequence *S);
Alignment * aln2unique_name_aln ( Alignment *S);
int name_list2unique_name_list (int n, char **name);
Sequence *seq2clean_seq ( Sequence *S, char *alp);//remove all alp characters from seq
char**gene2exons    (char **seq, int nseq);

int       ** seq2aln_pos      (Alignment *A, int *n, int **ls);


Alignment *padd_aln ( Alignment *A);
char **padd_string ( char **string, int n,char pad);

Alignment *local_maln2global_maln (char *seq, Alignment *A);

Alignment * seq2profile (Sequence *S, int index);

Sequence *remove_empty_sequence (Sequence *S);
Alignment *  aln2profile (Alignment * A);
Alignment * aln2collapsed_aln (Alignment * A, int n, char **string);
Alignment* remove_seq_from_aln (Alignment *A, char *seq);

Alignment* aln2sub_aln_file (Alignment *A, int n, char **string);
Alignment* aln2sub_seq (Alignment *A, int n, char **string);

int       ** aln2inv_pos  (Alignment *A);
int        * seq2inv_pos ( char *seq);
int       ** aln2pos_simple   (Alignment *A, int n_nseq, ...);
int       ** aln2pos_simple_2 (Alignment *A);
Alignment ** split_seq_in_aln_list ( Alignment **aln, Sequence *S, int l_seq, char **seq_list);

Sequence * fill_sequence_struc ( int nseq, char **sequences, char **seq_name, Genomic_info *genome_co);

int  seq_list2in_file ( TC_method *M, Sequence *S, char *list, char *file);
int  seq_list2fasta_file( Sequence *S,  char *list, char *file, char *outmode);
Structure * seq2struc    ( Sequence *S, Structure *ST);
Alignment *strings2aln (int nseq,...);

Alignment * seq2aln      ( Sequence *S, Alignment *A,int rm_gap);
Alignment  *seq_coor2aln ( Sequence *S, Alignment *A, int **coor, int nseq);

Alignment *stack_aln (Alignment *A, Alignment *B);
Alignment *chseqIaln(char *name, int seq_n, int start,int len,Sequence *S, int seqIaln, Alignment *A);


char *dna_aln2cons_seq ( Alignment *A);
char *aln2cons_seq ( Alignment *A, int ns, int *ls, int n_groups, char **group_list);
char *aln2cons_maj ( Alignment *A, int ns, int *ls, int n_groups, char **group_list);
Alignment *aln2conservation ( Alignment *A, int threshold,char *seq);


char *sub_aln2cons_seq_mat ( Alignment *A,int ns, int *ls, char *mat_name);
char *aln2cons_seq_mat ( Alignment*A, char *mat_name);
char *aln2cons_seq_cov ( Alignment*A);

Alignment *aln2short_aln( Alignment *A, char *list, char *nnew, int spacer);
Sequence  *keep_residues_in_seq ( Sequence *S,char *list, char replacement);
Alignment *keep_residues_in_aln ( Alignment *A,char *list, char replacement);
Alignment *filter_keep_residues_in_aln ( Alignment *A,Alignment *ST, int use_cons, int value, char *list, char replacement);

Alignment *aln_convert (Alignment *A, Alignment *ST, int use_cons, int value,int n, ...);
Alignment *aln2number (Alignment *A);
Alignment * filter_aln ( Alignment *A, Alignment *ST, int value);
Alignment * filter_aln_lower_upper ( Alignment *A, Alignment *ST,int use_cons, int value);
Alignment * filter_aln_upper_lower ( Alignment *A, Alignment *ST, int use_cons,int value);
Alignment * filter_aln_switchcase ( Alignment *A, Alignment *ST, int use_cons, int value);

Alignment * STseq2STaln ( Alignment *A, Alignment *ST);
Alignment * merge_annotation   ( Alignment *A, Alignment *ST, char *seq);
Alignment * filter_aln_convert ( Alignment *A, Alignment *ST, int use_cons,int value, int n_symbol,char** symbol_list);
int aln2ngap (Alignment *A);

int  * count_in_aln ( Alignment *A, Alignment *ST, int value, int n_symbol,char **symbol_list, int *table);
void count_misc (Alignment*A, Alignment *B);

char *msaF2fastaF(char *file);
int         trim_fastaF_big (char *in_aln1, char*in_aln2, char *out_aln1, char *out_aln2, long **map1, long **map2);
int         trim_aln_file (char *in_aln1, char*in_aln2, char *out_aln1, char *out_aln2);
Alignment * trim_aln_with_seq ( Alignment *S, Alignment *P);
Alignment * add_align_seq2aln ( Alignment *A, char *seq, char *seq_name);
Alignment * aln2X (Alignment *A, int x);
Sequence  * aln2seq    ( Alignment *A);
Sequence  * aln2seq_main    ( Alignment *A, int mode);
Alignment * thread_profile_files2aln (Alignment *A, char *template_file, Fname *F);
Alignment * expand_aln (Alignment *A);
Alignment * aln2expanded_aln (Alignment *A);
Alignment * expand_number_aln (Alignment *A,Alignment *EA);
Alignment * remove_gap_column ( Alignment *A, char *mode);
Alignment*  ungap_sub_aln        ( Alignment *A, int nseq, int *ls);
Sequence *  ungap_seq       ( Sequence *A);
Alignment * insert_gap_col (Alignment *A, int p, int l);
Alignment * unalign_residues (Alignment *A, int i1, int i2);
Alignment * unalign_aln_2 (Alignment *A, Alignment *C, int t);
Alignment * unalign_aln (Alignment *A, Alignment *C, int t);
Alignment * unalign_aln_pos (Alignment *A, int s, int p, int l);

Alignment *degap_aln (Alignment *A);

Alignment * RmLowerInAln (Alignment *A, char *gap);
char * ungap_fastaF_big        (char *in, char*out, int n);
Alignment* ungap_aln_n        ( Alignment *A, int n);
Alignment * ungap_aln        ( Alignment *A);
void compress_aln     ( Alignment *A);
Alignment* condense_aln (Alignment *A);

Alignment * probabilistic_rm_aa ( Alignment *A, int pos, int len);
Alignment * aln_gap2random_aa(Alignment *A);
Alignment * make_random_aln(Alignment *A,int nseq, int len, char *alphabet);
Alignment * add_random_sequence2aln( Alignment *A, char *alphabet);

int ** trim_aln_borders            ( char **seq1, char **seq2, int nseq);
Sequence * trim_aln_seq      ( Alignment  *A, Alignment *B);
Sequence * trim_aln_seq_name ( Alignment  *A, Alignment *B);
Sequence *get_defined_residues( Alignment *A);


Alignment *thread_defined_residues_on_aln ( Alignment *A, Sequence *S1);
Sequence *seq2number (Sequence *S);
Sequence * merge_seq    ( Sequence *IN, Sequence *OUT);
char * seq_name2coor ( char *s, int *start, int *end, char sep);
Alignment *seq_name2removed_seq_name(Sequence *S, Alignment *NA, float **diff);
int seq_name2index (char *name, Sequence *S);

Sequence *extract_one_seq(char *n,int start, int end, Alignment *S,int keep_name);
Sequence  * extract_sub_seq( Sequence  *COOR, Sequence *S);


Sequence * add_prf2seq  (char *alnfile, Sequence *S);
int prf_in_seq ( Sequence *S);
Sequence * add_sequence ( Sequence *IN, Sequence *OUT, int i);
Sequence  * trim_seq     ( Sequence   *A, Sequence  *B);
Sequence  * reorder_seq  ( Sequence   *A, char **name, int nseq);
Sequence  * reorder_seq_2  ( Sequence   *A, int **name,int field, int nseq);

char * concatenate_seq ( Sequence *S, char *conc, int *order);
Sequence * swap_header ( Sequence *S, Sequence *H);
char *template_file2abs_template_file(char *name);
Alignment *aln2jacknife (Alignment *A, int nseq, int len);
char ** name2random_subset (char **in_name, int n_in, int n_out);
Alignment * aln2random_order   ( Alignment  *A);
Alignment * aln2scramble_seq  ( Alignment  *A);

Alignment * shuffle_aln ( Alignment *A,int N, char *name_i, char *mode);
Alignment * reorder_aln        ( Alignment  *A, char **name, int nseq);
  
char ** rm_name_tag (char **name, int nseq, char *tag);

/******************************************************************************/
/*                          TEMPLATE MANAGEMENENT                             */
/******************************************************************************/
char * string_contains_template_tag (char *string);
Sequence * seq2template_type(Sequence *Seq);

Sequence * vremove_seq_template_files (Sequence *S);
Sequence * display_seq_template_files (Sequence *S);

Sequence * handle_seq_template_file (Sequence *S, char *mode);
int handle_X_template_files ( X_template *T, char *mode);


Sequence * seq2template_seq ( Sequence *S, char *template_file, Fname *F);
char * seq2template_file (Sequence *S, char *file);
int seq2template_file2 (Sequence *S, char *file, char *mode);

Sequence * profile_seq2template_seq ( Sequence *S, char *template_file, Fname *F);
int seq2n_X_template ( Sequence *S, char *type);

struct X_template *fill_X_template (char *name, char *p, char *type);
FILE * display_seq_template (Sequence *S, FILE *io);

char *template_type2type_name (char *type);
char *template_type2short_type_name (char *type);


FILE * display_sequence_templates ( Sequence *S, int i, FILE *io);
FILE * display_profile_templates (Sequence *S,int i, FILE *io);
FILE * display_X_template (struct X_template *X, FILE *io);

struct X_template* free_X_template ( struct X_template *X);

struct X_template *fill_P_template (char *name, char *p, Sequence *S);
int **seq2pdb_index (Sequence *S);
struct X_template *fill_F_template (char *name, char *p, Sequence *S);
struct X_template *fill_S_template ( char *name,char *p, Sequence *S);
struct X_template *fill_R_template (char *name, char *p, Sequence *S);
struct X_template *fill_G_template (char *name, char *p, Sequence *S);
struct X_template *fill_T_template (char *name, char *p, Sequence *S);
struct X_template *fill_E_template (char *name, char *p, Sequence *S);
struct X_template *fill_U_template (char *name, char *p, Sequence *S);

char *seq2T_value ( Sequence *S, int i, char *param_name, char *template_type);
char *profile2P_template_file (Sequence *S, int n);
Alignment * seq2R_template_profile (Sequence *S, int n);
char *seq2P_pdb_id (Sequence *S, int n);
char      * seq2P_template_file (Sequence *S, int n);
char      * seq2T_template_string (Sequence *S, int n);
char      * seq2E_template_string (Sequence *S, int n);
int       * seq2U_template (Sequence *S, int n);


int seq2n_template (Sequence *S, char *type);
struct X_template * seq_has_template ( Sequence *S, int n, char *type);

/******************************************************************************/
/*                          ALIGNMENT MANIPULATION                            */
/******************************************************************************/

char *aln_column2string (Alignment *A, int p);
int **fix_seq_aln (Sequence *S, Alignment*A, int **cache);
int **fix_seq_seq ( Sequence *S1, Sequence *S2);
int **fix_aln_seq_new (Alignment *S1, Sequence *S2);

Alignment * fix_aln_seq  ( Alignment *A, Sequence *S);
Alignment * rotate_aln ( Alignment *A, char *name);
Alignment * invert_aln ( Alignment *A);
Sequence  * invert_seq2 ( Sequence  *A);
int invert_seq_file (char *seq);
int invert_aln_file (char *seq);


char * complement_string (char *s);
Alignment * complement_aln ( Alignment *A);
Alignment * extract_nol_local_aln( Alignment *A, int start, int max_end);
Alignment * aln2block   (Alignment  *A, int start, int end, Alignment *B);
Alignment * alnpos2block   (Alignment  *A, int*pos, Alignment *B);

Alignment * extract_aln          ( Alignment *A, int start, int end);
Alignment * extract_aln2          ( Alignment *A, int start, int end, char *seq_name);
Alignment * extract_aln3          ( Alignment *A, char *filename);
Alignment * alnpos_list2block (Alignment *A, int n, char **in_list);

Alignment * trunkate_local_aln   ( Alignment *A);
int get_nol_aln_border ( Alignment *A, int start, int direction);
Alignment ** trim_local_aln ( Alignment *A, int **List, int ne, int **residue_list, Sequence *S);

Alignment * aln_cat ( Alignment *A, Alignment *B);
Alignment * concatenate_aln ( Alignment *A, Alignment *B, char *sep);
char * extract_defined_seq ( char *in, int in_of, int in_start, int *aa_def, int dir, int *out_start, char *out_seq);
int verify_aln ( Alignment *A, Sequence *S, char * error);
Alignment * remove_end (Alignment *A);
Alignment * orthologous_concatenate_aln (Alignment *A, Sequence *S, char *mode);
Alignment * aln2N_replicate (Alignment *A, char *nn, char *name);
FILE *aln2replicate (Alignment *A, FILE *fp);

Alignment * voronoi_concatenate_aln (Alignment *A, Sequence *S);

Alignment *adjust_est_aln ( Alignment *PW, Alignment *M, int s);
Alignment * rename_seq_in_aln (Alignment *A, char ***list);
Sequence * rename_seq_in_seq (Sequence *A, char ***list);
/********************************************************************/
/*                                                                  */
/*                   FLOAT SIMILARITIES                             */
/*                                                                  */
/*                                                                  */
/*                                                                  */
/********************************************************************/
float get_seq_fsim ( char *string1, char *string2, char *ignore, char *similarity_groups, int **matrix, int mode);
float get_seq_fsim2 ( char *string1, char *string2, char *ignore, char *in_mode);
float ** get_fsim_aln_array ( Alignment *A, char *mode);
/********************************************************************/
/*                                                                  */
/*                   ALIGNMENT ANALYSES                             */
/*                                                                  */
/*                                                                  */
/*                                                                  */
/********************************************************************/
int **sim_array2dist_array ( int **p, int max);
int **dist_array2sim_array ( int **p, int max);
int **normalize_array (int **p, int max, int norm);

int aln2most_similar_sequence ( Alignment *A, char *mode);
int aln2coverage ( Alignment *A, int ref_seq);

double* aln2column_normalized_entropy (Alignment *A);
double* aln2column_entropy (Alignment *A);
double aln2entropy (Alignment *A, int *in_ls, int in_ns, float gap_threshold);
int sub_aln2sim ( Alignment *A, int *ns, int **ls, char *mode);
int sub_aln2max_sim ( Alignment *A, int *ns, int **ls, char *mode);
int aln2sim2    ( Alignment *A);
int aln2sim     ( Alignment *A, char *mode);
int seq2idscore_sim ( char *seq1, char *seq2);

int aln_is_aligned ( Alignment *A);
int* get_cdna_seq_winsim ( int *cache, char *string1, char *string2, char *ignore, char *mode, int *w);
int  get_cdna_seq_sim    ( int *cache, char *string1, char *string2, char *ignore, char *mode);

int seq2aln2sim    (char *seq1, char *seq2, char *mode_aln, char *mode_id);
int* get_seq_winsim( char *string1, char *string2, char *ignore, char *mode, int *w);
int  get_seq_sim ( char *string1, char *string2, char *ignore, char *mode);
int  get_seq_sim_2 ( char *string1, char *string2, char *ignore, char **gr, int ng);
int  get_seq_sim_3 ( char *string1, char *string2, char *ignore, int **mat);


int *** get_winsim_aln_array ( Alignment *A, char *mode, int ***w);
int ** get_sim_master_aln_array ( Alignment *A,int n, char *mode);

int ** seq2sim_mat (Sequence *S, char *mode);
int ** seq2cov_mat (Sequence *S, char *mode);
int ** seq2comp_mat (Sequence *S, char *mode, char *comp_mode);

int logid_score (int sim, int len);
int ** fast_aln2sim_mat (Alignment *A, char *mode);
int ** fast_aln2sim_list (Alignment *A, char *mode, int *ns, int **ls);

int ** aln2dist_mat(Alignment *A);
int ** aln2dist_mat_gap (Alignment *A);
int ** aln2sim_mat_km (Alignment *A, char *mode);
int ** aln2sim_mat    (Alignment *A, char *mode);
int **aln2cov (Alignment *A);
int ** get_dist_aln_array ( Alignment *A, char *mode);
int ** get_raw_sim_aln_array ( Alignment *A, char *mode);
int ** array2sim (char **array, int n, char *mode);
int ** get_sim_aln_array ( Alignment *A, char *mode);
int generic_get_seq_sim  ( char *seq1, char *seq2, int *cache, char *mode);
Alignment * grep_seq (Alignment *S,char *field, char *mode, char *string);
Alignment* modify_seq (Alignment *S,char *field, char *string1, char *string2);

Sequence * seq2filter (Sequence *S_in, int min, int max);
int ** get_cov_aln_array ( Alignment *A, char *mode);
int ** get_cov_master_aln_array ( Alignment *A,int n, char *mode);


int * get_aln_col_weight ( Alignment *A, char *mode);
int analyse_aln_column   ( Alignment *B, int col);

int sub_aln2nseq_prf ( Alignment *A, int ns, int *ls);
int **aln2count_mat   (Alignment *A);
int **sub_aln2count_mat2   (Alignment *A, int ns, int *ls);
int **sub_aln2count_mat3   (char **al, int n);
int **aln2count_mat2   (Alignment *A);
char *aln2random_seq (Alignment *A, int noise1, int noise2, int noise3, int gap_noise);

int * km2centroids (Alignment *A, int k, char *mode,int *keep);
Alignment* km_seq (Alignment *S, int k, char *mode,char*name );
Alignment** seq2kmeans_subset (Alignment *A, int k, int *n, char *mode);
Alignment** seq2id_subset (Alignment *A, int k, int *n, char *mode);

int* seq2kmeans_class  (Alignment *A, int k, char *mode);

int aln2gap_trimmed (Alignment *A, int n, char *alnf, char *seqf);
Alignment *gap_trim (Alignment *A, int f);
Alignment * master_trimseq( Alignment *A, Sequence *S,char *mode);
Alignment * trimseq( Alignment *A, Sequence *S, char *mode);
int *aln2subset (Alignment *A, char *mode, int *n);
Alignment *simple_trimseq (Alignment *A,Alignment*K, char *mode, char *seq, int **sim);
Alignment *sim_filter (Alignment *A, char *in_mode, char *seq_list);

Sequence  * seq_weight2species_weight (Alignment *A, Sequence *S);
Alignment * aln2voronoi_weights (Alignment *A);
float ** get_weight ( Alignment *A, Sequence *S, char *mode);
float **seq2pwsim (	   Alignment *A, Sequence *S, char *mode);
Alignment * trimseq( Alignment *A, Sequence *S,char *mode);
Alignment * tc_trimseq( Alignment *A, Sequence *S,char *mode);
Alignment* seq2subseq3( Alignment *A, Sequence *S,int use_aln, int lower_sim,int upper_sim, int min_nseq, int trim_direction, char *weight_mode, float ***sim_weight, int *seq_list);
Alignment* seq2subseq2( Alignment *A, Sequence *S,int use_aln, int lower_sim,int upper_sim, int max_nseq, int trim_direction, char *weight_mode, float ***weight_table, int *seq_list);
float extreme_seq (int direction, Alignment *A,float **sim_weight,int *seq_list, int *seq_index);


Alignment* seq2subseq1( Alignment *A, Sequence *S,int use_aln, int percent,int max_nseq,int max_diff, char *weight_mode);
/********************************************************************/
/*                                                                  */
/*			AMINO ACID FUNCTIONS                        */
/*                                                                  */
/*                                                                  */
/*                                                                  */
/********************************************************************/
char** string2alphabet (char *string, int depth, int *falp_size);
int is_in_same_group_aa ( char r1, char r2, int n_group, char **gl, char *mode);
int find_group_aa_distribution (char *col, int nseq,int n_group, char **gl,  int *distrib, char *mode );
char** make_group_aa (int *ngroup, char *mode);
char** make_group_aa_upgma (char *mat, int max_size);


char * test_gene2prot (Constraint_list *CL, int s1);
Alignment* gene2prot (Alignment *A);
Alignment * dna_aln2_3frame_cdna_aln(Alignment *A,int *ns,int **l_s);

int ** get_sim_aln_array_normal_distribution ( Alignment *A, char *mode, int *STD, int *CENTER);
double normal(double x, double mean, double std);
int generic_get_seq_sim_normal_distribution ( char *seq1, char *seq2, int*cache, char *mode, int *STD, int *CENTER);
int get_seq_sim_distribution ( char *string1, char *string2, char *ignore, char *in_mode, int *STD, int *CENTER);

Alignment *aln2clean_pw_aln (Alignment *A,OveralnP *F);
char **pw_aln2clean_pw_aln (char ** aln,OveralnP *F);
int  * pw_aln2clean_aln_weight ( char *seq1, char *seq2, int w, OveralnP *F);

float* aln2pred  ( Alignment *A, Alignment*B, char *mode);
float* analyze_overaln ( Alignment *A, Alignment *B, char *mode, int f,int p1,int p2, int p3,int filter);


Alignment * mark_exon_boundaries  (Alignment *A, Alignment *E);
//simple_trimseq2
//Creates Clusters
//In each cluser there is a path between every pair of sequence
//A path is made of edges connecting tow nodes with w>min_sim
int ** simple_trimseq2 (int n, int **sim, int min_sim);
int thread_msa2msa(char *small, char *big, char *seq);
