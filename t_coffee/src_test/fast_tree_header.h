/*!
 * Node of a Tree
 */
typedef struct
{
	///number of this node
	int name;
	///left child
	int left;
	///right child
	int right;
//  ///weight of the node
// 	double weight;
	///position of the leafs during construction
	int *index;
	///position of...
	int *index2;
	///number of leafes below this node
	int num_leafs;
}
Tree_fastal;



/*!
 * Saves all important paramters for the PartTree algorithm
 *
 */
typedef struct
{
	///Size of the ktup
	short ktup;
	///power
	short ng;
	///number of possible ktup
	int tsize;
	///subgroup size
	int subgroup;
	///File with the ktup tables
	FILE *ktup_table_f;
	///Index of the positions in the ktup file
	long *ktup_positions;
	///number of sequences
	int num_sequences;
	///length of the sequences
	int *seq_lengths;
	///threshold for filter
	double threshold;
	/*! \brief An Integer array to save the tree
	*
	* Always three fields in this array belong together. The first two denote the number of the children the third one the number of the current node.
	*/
	Tree_fastal *tree;
	/// Current position in the tree
	int pos_tree;
}
PartTree_param;



//*******************************   METHODS   ***********************************


//*********** ktup-distance ***************
void makepointtable_fast(int *coded_seq, int ktup, int ng, int length);
void makecompositiontable_fastal(FILE* tables_f, int *table, int *pointt, int length);
void get_table(short *table, FILE *tables, long index);



void make_fast_tree(char *file_name, int n, int ktup);
int make_pos_len_index_of_file(char *file_name, char *ktable_f, long **file_positions, int **seq_lengths, int ktup, int is_dna);
int make_pos_len_index_of_file_retree(char *file_name, long **file_positions, int **seq_lengths);
// void ktup_one2all(FILE *seq_f, int* file_positions, int *least_similar);


int **make_distance_matrix(FILE *table_f, long *file_positions, int *sequence_group, int number, int **dist_mat);
int euclidean_dist(FILE* ktup_f, long pos1, long pos2, short *table1, short *table2, int length);
int euclidean_dist_half(FILE* ktup_f, long pos2, short *table1, short *table2, int length);

void partTree_r(PartTree_param *param_set);
void partTree_retree(PartTree_param *param_set);
void make_partTree(char *sequence_f, char *tree_f, int ktup, int subgroup, int is_dna, int retree);

int filter(int *sequence_group, int **dist_mat, int *seed_set_cleaned, PartTree_param *param_set);
int fast_aln2sim_mat2 (char *seq1, char *seq2);







/**
 * \brief Struct to save the features for a set of sequences.
 * 
 */
typedef struct
{
	/// The number of this sequence.
	int seq_number;
	/// The feature to use in qsort.
	int *elem_2_compare;
	/// An element of features.
	double *features;
}
Cluster_info;


int				cluster_recursive(int start, int end, Cluster_info* matrix, int dim, double *mean, double *variance, int *elem_2_compare, FILE* tree_f, int *node);
int*			recode_sequence(char * seq, int seq_length, int *char2value, int alphabet_size, int word_length);
double*			get_features(int *recoded_seq, int seq_length, int k_tup, int max_coding, int max_dist, int num_features);
Cluster_info*	feature_extract(char *seq_file_name, int max_dist, int k_tup, int *char2value, int alphabet_size, int *elem_2_compare, int *num_seq_p, int num_features);
void			compute_oligomer_distance_tree(char *seq_file_name, int* char2value, char *tree_file_name, int max_dist, int word_length, int alphabet_size);

