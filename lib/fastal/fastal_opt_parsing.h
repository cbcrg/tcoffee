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
	char *mat;
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

