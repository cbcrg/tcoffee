

//sum of pairs score
double calculate_sum_of_pairs_score_affine(char *alignment_file_name, int **score_matrix, double gop, double gep);


//compare with reference alignment
void make_ref_alignment(char *seq_file_name, char *tree_file_name, char *ref_aln_name, int num_seq_in_ref);
double agreement_score(char *ref_file_name, char *aln_file_name);
