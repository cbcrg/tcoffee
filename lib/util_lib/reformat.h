#ifndef __REFORMAT_H
#define __REFORMAT_H

typedef struct
    {
      Alignment *A;
      Weights *W;
      Sequence *S;
      int **M;
      Structure *RNA_ST;
      NT_node T;
      Constraint_list *CL;
      char format[100];
      char file[10000];
      int rm_gap;

}Sequence_data_struc;

typedef struct
    {
	char **symbol_list;
        int n_symbol;
        char *coor_file;
        int rm_gap;
        int keep_case;
        int keep_name;
        int use_consensus;
}Action_data_struc;

/*Control of alignment sizes*/
int  set_landscape_msa (int len);
int get_msa_line_length (int line, int aln_len);

int seq_reformat (int argc, char **argv);

Sequence_data_struc *read_data_structure ( char *in_format, char *in_file,Action_data_struc *RAD);
Alignment * read_fasta_aln_noceck ( char *name, Alignment *A);

Alignment * main_read_aln ( char *name, Alignment *A);
int big();
int *list2sample(int nseq, int sample);
char*FastaRecord2name(char*rec);
char*FastaRecord2seq(char*rec);
char *FastaRecord2header (char *record);
char *seq2clean  (char*seq);
char *file2record (char *file, int i, long *map);
char *file2record_it (char *file, int i, long *map);
long *fasta2map(char *file);
int fasta2nseq (char *file);


Sequence  * quick_read_seq ( char *file);
Alignment *reload_aln (Alignment*A);
char *      split_fasta (char *file, int size);
Alignment * quick_read_aln( char *name);//reads onl�y fasta and clustalw
Alignment * quick_read_fasta_aln( Alignment *A,char *name);//reads onl�y fasta and clustalw
Alignment * quick_read_aln_static( char *name);//reads onl�y fasta and clustalw
int read_nameseq (char *file,char ***nam, char ***seq, char ***com);

Sequence  * read_sequences ( char *name);
Sequence  * read_alifold   ( char *name);
Alignment *alifold2aln     ( char *name);
Sequence  * main_read_seq ( char *mname);
int output_format_aln ( char *format, Alignment *A, Alignment *EA,char *name);
int main_output   ( Sequence_data_struc *D1, Sequence_data_struc *D2, Sequence_data_struc *DST, char *out_format, char *out_file);
int is_aligned (int n,char **al);
char * identify_seq_format ( char *file);
char * name2type_name ( char *name);
char identify_format (char **fname);
char **identify_list_format ( char **list, int n);

int type_is_exon_boundaries(char **seq, int n);

int format_is_oligo  ( char *file);
int format_is_msf  ( char *file);
int format_is_fasta( char *file);
int format_is_not_fasta (char *file);
// int format_is_fasta_aln( char *file);
int format_is_fasta_aln ( char *file, int i_know_that_it_not_seq);
int format_is_fasta_seq( char *file);
int is_pir_name (char *name);
int format_is_pir  ( char *file);
int format_is_pir_aln( char *file);
int format_is_pir_seq( char *file);
int pir_name (char *name);
int format_is_conc_aln (char *file);
int format_is_saga  ( char *file);
int format_is_swissprot (char *name);

int is_seq ( char *name);
int is_aln ( char *name);
int has_pdb (char *name);
int is_stockhom_aln ( char *name);
int is_blast_file (char *name);
int is_sap_file (char *name);
int is_pdb_file ( char *name);
int is_simple_pdb_file ( char *name);
char *fix_pdb_file (char *name);
char *atom2pdbF(char*file);
int pdb_has_atom ( char *name);

int seqres_equal_atom (char*fname);
int is_pdb_name ( char *name);
char* get_pdb_id(char *name);
char* get_pdb_struc(char *name, int start, int end);
char*  seq_is_pdb_struc ( Sequence *S, int i);
char* is_pdb_struc ( char *name); /*Returns NULL if not a PDB structure Or a the name of a file containing a PDB structure*/
int is_matrix (char *name);

int is_lib (char *name);
int is_lib_01 (char *name);
int is_lib_02 (char *name);
int is_lib_list ( char *name);
int is_single_seq_weight_file (char *fname);
int is_treelist(char *name);
int is_newick  (char *name);
int is_mafft_newick  (char *name);
int is_mafft_newick  (char *name);
int is_nexus (char *file);
int is_nameseq (char *file);
int is_method ( char *file);

char *format_name2aln_format_name (char *name);
int is_in_format_list ( char *name);
int is_out_format_list ( char *name);
int is_struc_in_format_list ( char *name);
int is_struc_out_format_list ( char *name);
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               INPUT MISC                                               */
/*                                                                                         */
/***************************************************************************************** */

char *** read_rename_file ( char *fname, int mode);
void get_barton_list_tc_seq ( char *in_file);
int process_barton_entry (char *buf, char *name);

Structure *read_rna_struc_number ( Alignment *A, char *fname);
char ** read_lib_list (char *name, int *n);
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               INPUT WEIGHTS                                             */
/*                                                                                         */
/***************************************************************************************** */
Weights* get_amps_sd_scores ( char *fname);
Weights *read_seq_weight (char **name, int nseq, char* seq_weight);
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               INPUT SEQUENCES                                            */
/*                                                                                         */
/***************************************************************************************** */
char ***read_group ( char *file);
int *pdb2atom_pos_list(char *pdb);
Sequence* get_pdb_sequence_from_field   (char *fname, char *field);
Sequence* get_pdb_sequence           ( char *fname);
Sequence* get_struc_gor              ( char *fname);
Sequence* get_dialign_sequence       ( char *fname);
Sequence* get_pima_sequence          ( char *fname);
Sequence* get_sequence_dali          ( char *fname);
Sequence* get_pir_sequence           ( char *fname, char *comment_name);
Sequence* perl_reformat2fasta        ( char *perl_script, char *file);

int get_next_fasta_sequence          (FILE*fp, char **name, char **comment, char **seq);
Sequence *reload_seq(Sequence *A);
Sequence* get_fasta_sequence         ( char *fname, char *comment_name);
Sequence* get_fasta_sequence_num     ( char *fname, char *comment_name);
Sequence* get_fasta_sequence_raw     ( char *fname, char *comment_name);
Sequence *get_file_list ( char *fname);
Sequence *get_tree_file_list ( char *fname);

Sequence* get_gor_sequence           ( char *fname, char *comment_name);
Sequence* get_swissprot_sequence     ( char *fname, char *comment_name);
int  fscanf_seq_name ( FILE *fp, char *sname);

void read_check ( Alignment *A, char *check_file);
void read_stockholm_aln ( char *fname, Alignment *A);
void read_aln ( char *fname, Alignment *A);
void read_number_aln ( char *fname, Alignment *A);
Alignment *read_blast_aln  ( char *fname, Alignment *A);
void read_msf_aln ( char *fname, Alignment *A);
void read_amps_aln ( char *in_file, Alignment *A);
int get_amps_seq_name ( char **name, char* fname);
Alignment *read_gotoh_aln ( char *fname, Alignment *A);

Alignment* undump_msa ( Alignment *A, char *tmp);
Sequence * undump_seq ( Sequence  *A, char *tmp);
char* dump_msa   ( Alignment *A, char *tmp);
char* dump_seq   ( Sequence  *A, char *tmp);

/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               OUTPUT MATRICES                                           */
/*                                                                                         */
/***************************************************************************************** */
int output_freq_mat ( char *outfile, Alignment *A);
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               OUTPUT P-Values                                           */
/*                                                                                         */
/***************************************************************************************** */
float output_maln_pval ( char *outfile, Alignment *A);
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               OUTPUT WEIGHTS                                            */
/*                                                                                         */
/***************************************************************************************** */
void  output_similarities (char *file, Alignment *A, char *mode);
void  output_similarities_pw (char *file, Alignment *A, Alignment *B, char *mode);
Alignment * similarities_file2aln ( char *file);
int** input_similarities (char *file, Alignment *A, char *mode);

void output_statistics (char *file, Alignment *A, char *mode);
void output_pw_weights4saga ( Weights *W, float **w_list, char *wfile);
int  output_seq_weights ( Weights *W, char *wfile);
FILE * display_weights (Weights *W, FILE *fp);
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               OUTPUT SEQ                                                */
/*                                                                                         */
/***************************************************************************************** */
char** clean_seq_names (char **names, int n, int mode);
char *clean_seq_name (char *name, int mode);


void output_pir_seq1 (char *fname, Alignment*A );
void output_pir_seq (char *fname, Alignment*A );
void output_gor_seq (char *fname, Alignment*A );
void output_mult_fasta_seq (char *fname, Alignment*A, int n );

void main_output_fasta_seq ( char *fname, Alignment *A, int header);
void output_fasta_simple   ( char *name, Sequence *S);
void output_fasta_seqS (char *fname, Sequence *S );
void output_fasta_seq1 (char *fname, Alignment*A );
void output_fasta_seq2 (char *fname, Alignment*A );
char *output_fasta_seqX (char *name, char *mode, Sequence *S, Alignment *A, int i);

void output_pir_check (char *fname,int nseq, char **A );
void output_fasta_seq (char *fname, Alignment*A );
void output_gotoh_seq (char *fname, Alignment*A );
void output_est_prf   (char *fname, Alignment *A);
void output_gor_seq (char *fname, Alignment*A );
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               OUTPUT ALN                                                */
/*                                                                                         */
/***************************************************************************************** */
void output_pir_aln    ( char *fname,Alignment*A);
void output_model_aln  ( char *fname,Alignment*A );
char * output_fasta_sub_aln (char *fname, Alignment*A, int ns, int *ls  );
char * output_fasta_sub_aln2 (char *fname, Alignment*A, int *ns, int **ls  );

void ouput_suchard_aln ( char *fname,Alignment*A);
void output_fasta_aln  ( char *fname, Alignment *A);
void output_mfasta_aln  ( char *fname, Alignment *A);

void output_xmfa_aln  ( char *fname, Alignment *A);
void output_msf_aln    ( char *fname,Alignment*B);
FILE * output_generic_interleaved_aln (FILE *fp, Alignment *B, int line, char gap, char *mode);
void output_stockholm_aln (char *file, Alignment *A, Alignment *ST);
void output_clustal_aln( char *name, Alignment*B);
void output_strict_clustal_aln( char *name, Alignment*B);
void output_generic_clustal_aln( char *name, Alignment*B, char *format);
void output_saga_aln   ( char *name, Alignment*B);
void output_rphylip_aln ( char *name, Alignment*B);
void output_phylip_aln ( char *name, Alignment*B, char *mode);
void output_mocca_aln  ( char *name, Alignment*B,Alignment*S);
void output_rnalign    (char *out_file, Alignment*A,Sequence *STRUC);
void output_pw_lib_saga_aln (char *lib_name, Alignment *A );
void output_lib        (char *lib_name, Alignment *A );
void output_compact_aln( char *name, Alignment *B);

void print_sub_aln ( Alignment *B, int *ns, int **ls);
void print_aln ( Alignment *B);
FILE * output_aln( Alignment *B, FILE *fp);


FILE * output_aln_score ( Alignment *B, FILE *fp);
FILE * output_aln_with_res_number ( Alignment *B, FILE *fp);


FILE* output_Alignment ( Alignment *B, FILE *fp);
FILE* output_Alignment_without_header ( Alignment *B, FILE *fp);
FILE * output_Alignment_score ( Alignment *B, FILE *fp);
FILE * output_Alignment_with_res_number ( Alignment *B, FILE *fp);
void output_constraints ( char *fname, char *mode, Alignment *A);

Alignment *input_conc_aln ( char *name, Alignment *A);
void output_conc_aln ( char *name, Alignment *B);
void output_glalign       ( char *name, Alignment *B, Alignment *S);
void output_lalign_header( char *name, Alignment *B);
void output_lalign       ( char *name, Alignment *B);
void output_lalign_aln   ( char *name, Alignment *B);

/**************************************************************************************************/
/*                                                                                                */
/*                                                                                                */
/*                               INPUT/OUTPUT MATRICES                                                  */
/*                                                                                                */
/**************************************************************************************************/
int is_blast_matrix (char *fname);
int is_pavie_matrix (char *fname);
int is_clustalw_matrix (char *fname);

int is_distance_matrix_file (char *name);
int is_similarity_matrix_file (char *name);

void aln2proba_mat (Sequence *S);
void aln2mat (Sequence *S);
void aln2mat_diaa (Sequence *S);
int **seq2latmat ( Sequence *S, char *fname);
int output_mat (int **mat, char *fname, char *alp, int offset);
int ** read_blast_matrix ( char *mat_name);
int output_blast_mat (int **mat, char *fname);
double* mat2cmp (int **mat1, int **mat2);

void output_pavie_mat (int **mat, char *fname, double gep, char *alp);
int ** read_pavie_matrix ( char *mat_name);

/****************************************************************************************************/
/***************************                                    *************************************/
/***************************             PROCESSING 		*************************************/
/***************************                                    *************************************/
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                              THREADING                                                  */
/***************************************************************************************** */




Structure * declare_rna_structure_num (Sequence *SA);

char *thread_aa_seq_on_dna_seq( char *s);
void thread_seq_struc2aln ( Alignment *A, Sequence *ST);
Alignment *thread_dnaseq_on_prot_aln (Sequence *S, Alignment *A);
void cache_id ( Alignment *A);



int process_est_sequence ( Sequence *S, int *cluster_list);
char * invert_seq ( char *seq);
int get_best_match ( char *seq1, char *seq2);
int** extract_m_diag_streches ( int ** m, int l1, int l2,char *seq1, char *seq2, int *n_mdiag);
int is_strech ( char *AA, char *seq1, char *seq2, int len, int x, int y);

int search_for_cluster ( int seq, int cluster_number, int *cluster_list, int T, int nseq, int **S);
int * SHC ( int nseq, int **NST, int **ST);
int mutate_sol (int *sol, int nseq);
int evaluate_sol ( int*sol, int nseq, int **ST, int **NST);



char **make_symbols ( char *name, int *n);
Alignment *code_dna_aln (Alignment *A);
char* back_translate_dna_codon ( char aa, int deterministic);
int translate_dna_codon ( char *seq, char stop);
char* mutate_amino_acid ( char aa, char *mode);
Alignment * mutate_aln ( Alignment *A, char *r);

int extend_seqaln (Sequence *S, Alignment *A);
int unextend_seqaln (Sequence *S, Alignment *A);
char *extend_seq (char *seq);
char *unextend_seq (char *seq);

Sequence * transform_sequence ( Sequence *S, char *mode);
Alignment *translate_splice_dna_aln (Alignment *A,Alignment *ST );
Alignment * mutate_cdna_aln ( Alignment *A);

char *test_dna2gene (char *dna, int *w);
Sequence *dnaseq2geneseq (Sequence *S, int **w);

int ** shift_res_weights ( Sequence *R, int **w, int shift);
int res_weights2min(Sequence *R, int **w);
int res_weights2max(Sequence *R, int **w);
int res_weights2avg(Sequence *R, int **w);
int output_wexons (char *name, Alignment *A);
int scan_res_weights4ac (Sequence *R, int **w, int start, int end, int step);
float *res_weights2accuracy_counts ( Sequence *R, int **w,int T, float *result);
float* genepred_seq2accuracy_counts (Sequence *R, Sequence *T,float *result);
void genepred_seq2accuracy_counts4all (Sequence *R, Sequence *Ts); //JM
float* genepred2accuracy_counts     (char *ref,  char *target , float *result);

char *dna2gene (char *dna, int *w);
char * translate_dna_seq_on3frame (  char *dna_seq, char stop, char *prot);

char * translate_dna_seq ( char *dna_seq, int frame, char stop, char *prot);
int is_stop (char r1, char r2, char r3);
int seq2tblastx_db (char *file,Sequence *S, int strand);
int seq2blastdb (char *out, Sequence *S);
char * back_translate_dna_seq ( char *in_seq,char *out_seq, int mode);
Alignment *back_translate_dna_aln (Alignment *A);
Sequence  *translate_dna_seqS     (Sequence *S, int frame, int stop);
Alignment *translate_dna_aln (Alignment *A, int frame);
char *dna_seq2pep_seq (char *seq, int frame);

Alignment *clean_gdna_aln (Alignment *A);
Alignment *clean_cdna_aln (Alignment *A);
Alignment *clean_est      (Alignment *A);
/**************************************************************************************************/
/********************************                      ********************************************/
/********************************    PROCESSING        ********************************************/
/*************** ****************                      ********************************************/
void modify_data  (Sequence_data_struc *D1, Sequence_data_struc *D2, Sequence_data_struc *DST, char **action_list,int n_actions, Action_data_struc *RAD);

//
// Name MAnipulation
//

Alignment *clean_aln (Alignment *A);
Sequence *clean_sequence ( Sequence *S);
char ** translate_names (int n, char **name);
char * translate_name ( char *name);
char *decode_name (char *name, int mode);
FILE * display_sequences_names (Sequence *S, FILE *fp, int check_pdb_status, int print_templates);
Sequence *add_file2file_list (char *name, Sequence *S);

#endif
