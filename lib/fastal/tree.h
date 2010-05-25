




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

