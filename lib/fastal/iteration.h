

void seq2profile2(char *seq, Fastal_profile *prf, int *char2pos);
void split_set(FILE *aln_file_name, Fastal_profile *gap_prf, Fastal_profile *no_gap_prf, char *seq, int index, int *char2pos, char* split_file_name);
void iterate(Fastal_param *param, void *method_arguments_p, char *aln_file_name, char *out_file_name, int iteration_number);

void edit2seq_pattern(FILE *edit_file, char *seq1, char *seq2);
int *del_gap_from_profile(Fastal_profile *prf, int alphabet_size, int *gap_list, int *gap_list_length, int *num_gaps);

void write_iterated_aln(char* old_aln_file_name, char* new_aln_file_name, char *gap_file_name, char *seq1, int *gap_list1, int num_gap1, char *seq2, int *gap_list2, int num_gap2);
