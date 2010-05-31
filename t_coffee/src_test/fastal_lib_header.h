

#define ENLARGEMENT_PER_STEP 50
#define PROFILE_ENLARGEMENT 550

// static char pos2aa[] = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};


/**
 * \brief Struct to save a diagonal
 * 
 */
typedef struct
{
	///Start of diagonal in seq1.
	int x;
	///Start of diagonal in seq2.
	int y;
	///Length of the diagonal.
	int length;
	///expansion at the beginning;
	int front_exp;
	///expansion at the end;
	int end_exp;
}
Diagonal;




/**
 * \brief Struct to save all informations of a profile.
 * 
 */
typedef struct
{
	/// Number of sequences in this profile
	int num_sequences;
	/// number of the profile
	int prf_number;
	///0 = combination of two profiles, 1 = profile of a single sequence -> name1 = seq_name
	int is_leaf;
	///length of the profile
	int length;
	///weight of the sequence
	int weight;
	///saves the amount of allocated memory
	int allocated_memory;	
	///the profile itself [alphabet_size][profile_length]
	int **prf;
	///number_of_sequences
	int number_of_sequences;
}
Fastal_profile;

/**
 * \brief Struct to save all parameters for fastal.
 * 
 */
typedef struct
{
	/// size of alphabet_size
	int alphabet_size;
	/// converting char2position (for profile)
	int char2pos[256];
	/// converting pos2char (for profile)
	char pos2char[256];
	/// gap opening costs
	double gop;
	/// gap extension costs
	double gep;
	/// nomatch???
	int nomatch;
	///method to align profile
	char method[20];
	///scoring Matrix;
	int **M;
	///saves the diag method -> move to sparse!
	char diag_method[30];
}
Fastal_param;

/**
 * \brief Struct to save the parameters and memory for the sparse dynamic programming algorithm.
 */
typedef struct
{
	/// saves the diagonals
	int *diagonals;
	/// saves ....
	int *dig_length;
	/// list of points to be considered during the alignment process
	int **list;
	/// number of points in \a list
	int *list_length;
	char *file_name1;
	char *file_name2;
	// 	static char *file_name1 = vtmpnam(NULL);
// 	static char *file_name2 = vtmpnam(NULL);
}
Sparse_dynamic_param;


typedef struct
{
	int window_length;
	int step_legnth;
	int threshold;
	int oligo_length;
	int distance;
}
Udisc_param;


/**
 * \brief Struct to save the parameters and memory for the needleman-wunsch algorithm.
 */
typedef struct
{
	/// dynamic programming matrix
	double ** dyn_matrix;
	/// length of dimension 1
	int *length1;
	/// length of dimension 2
	int *length2;
	/// summed up version of profile
	int **sumup_prf;
	/// number of entries in \a sumup_prf
	int *sumup_length;
}
Nw_param;



/**
 * \brief Struct to save the parameters and memory for the Gotoh algorithm.
 */
typedef struct
{
	/// dynamic programming matrix
	double ** m_matrix;
	/// dynamic programming matrix
	double ** i_matrix;
	/// dynamic programming matrix
	double ** d_matrix;
	/// length of dimension 1
	int *length1;
	/// length of dimension 2
	int *length2;
	/// summed up version of profile
	int **sumup_prf;
	/// number of entries in \a sumup_prf
	int *sumup_length;
	/// log saver
	double **log_saver;
}
Gotoh_param;




typedef struct
{
	int diagonal;
	int count;
	int start;
	int end;
}
Swift_diagonal;


typedef struct
{
	int diagonal;
	int count;
}
Diagonal_counter;

//tree
void generate_random_tree(int number);


Fastal_profile* make_profile_of_sequence(char *seq_name, char *sequence, int number);



//Definite use

//*********************    input/output    **********************************
void file2profile(FILE* profile_f, Fastal_profile *profile, int prf_number, Fastal_param *param_set);
void file_pos2profile(FILE *seq_file, long off_set, Fastal_profile *profile, int prf_number, Fastal_param *param_set);
void profile2file(Fastal_profile *profile, FILE* prf_f, Fastal_param *param_set);

//index
int make_index_of_file(char *file_name, long **result);


//*********************    pairwise alignment methods     ************************

	//Needleman-Wunsch
	int prf_nw(Fastal_profile *profile1, Fastal_profile *profile2, double **prog_matrix, FILE *edit_file_name, int **sumup_prf, int *sumup_length, Fastal_param *param_set);
	int nw_matrix2edit_file(double **prog_matrix, Fastal_profile *profile1, Fastal_profile *profile2, FILE *edit_f, int **prf_field, int *field_length, Fastal_param *param_set);
	int** sumup_profile(Fastal_profile *profile, int **sumup_prf, Fastal_param *param_set);
	void write2file(int **sumup_prf, int length, FILE *file, int number, int num_sequences, Fastal_param *param_set);

	//Gotoh
	void free_gotoh(Gotoh_param* method_arguments_p, int alphabet_size);
	void fill_arguments_gotoh(Gotoh_param* method_arguments_p, int alphabet_size);
	int prf_gotoh(Fastal_profile *profile1, Fastal_profile *profile2, FILE *edit_file_name, Gotoh_param *arguments, Fastal_param *param_set);
			
	//Sparse dynamic programming
	void free_sparse(Sparse_dynamic_param* method_arguments_p);
	void fill_arguments_sparse(Sparse_dynamic_param* method_arguments_p);
	int **diagonals2int(int *diagonals, int num_diagonals, char *seq1, char *seq2, int *num_points, Fastal_param *param_set);
	int seq_pair2blast_diagonal(char *seq_file_name1, char *seq_file_name2, int **diagonals, int *dig_length, int l1, int l2, int is_dna);
	int sparse_dyn(Fastal_profile **profiles, Fastal_param *param_set, void *method_arguments_p, int is_dna, FILE *edit_file, FILE *prof_file, int number);
	char *profile2consensus(Fastal_profile *profile, Fastal_param *param_set);
	int ** diagonals2int_gap_test(int *diagonals, int num_diagonals, char *seq1, char *seq2, Fastal_profile *profile1, Fastal_profile *profile2, int *num_points, Fastal_param *param_set);
	int ** diagonals2int_euclidf(int *diagonals, int num_diagonals, char *seq1, char *seq2, Fastal_profile *profile1, Fastal_profile *profile2, int *num_points, Fastal_param *param_set);
	int ** diagonals2int_dot(int *diagonals, int num_diagonals, char *seq1, char *seq2, Fastal_profile *profile1, Fastal_profile *profile2, int *num_points, Fastal_param *param_set);
	int seq_pair2blastz_diagonal(char *seq_file_name1, char *seq_file_name2, int **diagonals, int *dig_length, int l1, int l2, int is_dna);
	int list2linked_pair_wise_fastal(Fastal_profile *prf1, Fastal_profile *prf2, Fastal_param *param_set, int **list, int n, FILE *edit_f, FILE *prof_f, int node_number);
	
//edit_files 2 alignment
void edit2alignment(FILE *sequence_file, long *seq_positions, FILE *edit_file, long *edit_positions, int node_number, int number_of_sequences, char *aligned_sequence, int alignment_length, FILE *edit_seq_file, int offset, FILE* alignment_file);
void edit_seq2aligned_seq(char *aligned_sequence, FILE *sequence_file, long sequence_position, FILE *alignment_file);


//main
int fastal(int argc, char **argv);
void alignment2files(Fastal_profile **profiles, Fastal_param *param_set,int **alignment, int alignment_length, FILE *edit_f, FILE *prof_f, int node_number);

//toolbox
double calculate_sum_of_pairs_score_affine(char *alignment_file_name, int **score_matrix, double gop, double gep);
double calculate_sum_of_pairs_score_affine_test(char *alignment_file_name, int **score_matrix, double gop, double gep);
void initiate_profile_files(FILE **profile_files);
void initiate_profiles(Fastal_profile **profiles, Fastal_param *param_set);
void free_fastal_profile(Fastal_profile *profile, int alphabet_size);
double **resize_dyn_matrix(double **dyn_matrix, int old_length1, int old_length2, int length1, int length2);
void free_dyn_matrix(int length1, double **dyn_matrix);
void fill_parameters(int is_dna, Fastal_param *param_set, char *method, char *diag_method);


int seq_pair2diagonal_own(char *seq1, char *seq2, int **diagonals, int *dig_length, int l1, int l2, int is_dna, int word_length);







typedef struct
{
	///field saving the positions [x1,y1,l1,x2,y2,l2,...]
	int *segments;
	/// points to the current used segment
 	int *current_pos;
	/// saves the previous diagonal position.
	int prev_position;
	/// diagonal number
// 	int diagonal_num;
}
Segment;


Segment* extend_diagonals(Diagonal *diagonals, int *num_diagonals, int l1, int l2);
int seq_pair2blast_diagonal2(char *seq_file_name1, char *seq_file_name2, Diagonal **diagonals, int *dig_length, int l1, int l2, int is_dna);

int ** segments2int(Segment *diagonals, int num_diagonals, char *seq1, char *seq2, Fastal_profile *profile1, Fastal_profile *profile2, int *num_points, Fastal_param *param_set);
// #include "string.h"

struct Fastal_arguments
{
// 	char *args[2];
	//IO
	
	
	char *sequence_file;
	char *tree_file;
	char *output_file;
	int tree_out;
	
	char *method;
	int is_dna;
	double gop;
	double gep;
	
	

	char *diag_method;
// 	Tree Computation
// 	int retree;
	int gap_iterate;
	int tree_only;
	char *tree_method;
	int tree_param1;
	int tree_param2;
	
	
	//Scoring
	int agreement_score;
	int evaluate;
	int score;
	int num_ref_aln;
	int num_seq_in_ref;
	char *aln_ref;
	char *aln2test;
};




void tree_parse(struct Fastal_arguments *arguments, char* param);



void arg_parse (int argc, char **argv, struct Fastal_arguments *arguments);



//sum of pairs score
double calculate_sum_of_pairs_score_affine(char *alignment_file_name, int **score_matrix, double gop, double gep);


//compare with reference alignment
void make_ref_alignment(char *seq_file_name, char *tree_file_name, char *ref_aln_name, int num_seq_in_ref);
double agreement_score(char *ref_file_name, char *aln_file_name);


void seq2profile2(char *seq, Fastal_profile *prf, int *char2pos);
void split_set(FILE *aln_file_name, Fastal_profile *gap_prf, Fastal_profile *no_gap_prf, char *seq, int index, int *char2pos, char* split_file_name);
void iterate(Fastal_param *param, void *method_arguments_p, char *aln_file_name, char *out_file_name, int iteration_number);

void edit2seq_pattern(FILE *edit_file, char *seq1, char *seq2);
int *del_gap_from_profile(Fastal_profile *prf, int alphabet_size, int *gap_list, int *gap_list_length, int *num_gaps);

void write_iterated_aln(char* old_aln_file_name, char* new_aln_file_name, char *gap_file_name, char *seq1, int *gap_list1, int num_gap1, char *seq2, int *gap_list2, int num_gap2);
