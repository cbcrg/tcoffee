






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
