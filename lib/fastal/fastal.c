#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#include "coffee_defines.h"
#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"
#include "fastal_lib_header.h"


// #include <omp.h>
// #define CHUNKSIZE 100
// #define N     1000


//Fastal_param *param_set;


/*!
 *	\file fastal.c
 *	\brief Source code for the fastal algorithm
 *	\author Carsten Kemena
 */


// the following functions are defined in scoring.c
int complete_agreement_score(char *aln_file_name, const char *ref_directory); 
void compute_ref_alignments(char *seq_file_name, char* ref_directory, int num_alignments, int num_seq_in_aln);


//**************************   sparse dynamic aligning **********************************************************


void
fill_arguments_sparse(Sparse_dynamic_param* method_arguments_p)
{
	method_arguments_p->diagonals = (int*)vcalloc(3,sizeof(Diagonal));
	method_arguments_p->dig_length =(int*) vcalloc(1,sizeof(int));
	*method_arguments_p->dig_length = 3;
	method_arguments_p->list = NULL;
	method_arguments_p->list_length = (int*)vcalloc(1,sizeof(int));
	*method_arguments_p->list_length = 0;
 	method_arguments_p->file_name1 = vtmpnam(NULL);
 	method_arguments_p->file_name2 = vtmpnam(NULL);
}

void
free_sparse(Sparse_dynamic_param* method_arguments_p)
{
	vfree(method_arguments_p->diagonals);
	vfree(method_arguments_p->dig_length);
	vfree(method_arguments_p->list_length);
}


/**
 * \brief One run of sparse dynamic programming.
 *
 * \param profiles The profiles.
 * \param param_set The fastal parameters.
 * \param method_arguments_p The method arguments.
 * \param is_dna Sequences are DNA (\a is_dna = 1) or protein.
 * \param edit_file The edit file.
 * \param prof_file the profile file.
 * \param number Number of the parent node.
 * \return The length of the alignment.
 */
int
sparse_dyn(Fastal_profile **profiles,
		   Fastal_param *param_set,
		   void *method_arguments_p,
		   int is_dna,
		   FILE *edit_file,
		   FILE *prof_file,
		   int number)
{
// 	printf("WHAT THE HELL ARE YOU DOING HERE?\n");
	Sparse_dynamic_param *arguments = (Sparse_dynamic_param*)method_arguments_p;
// 	static char *file_name1 = vtmpnam(NULL);
// 	static char *file_name2 = vtmpnam(NULL);
	char *file_name1 = arguments->file_name1;
	char *file_name2 = arguments->file_name2;
	char *seq1, *seq2;
	Fastal_profile *tmp1 = profiles[0];
	Fastal_profile *tmp2 = profiles[1];

	seq1 = profile2consensus(tmp1, param_set);
	seq2 = profile2consensus(tmp2, param_set);


	int **diagonals_p = &(arguments->diagonals);
	int num_diagonals = -1;
	if (!strcmp(param_set->diag_method, "blastz"))
	{
		FILE *cons_f = fopen(file_name1,"w");
		fprintf(cons_f, ">%i\n", tmp1->prf_number);
		fprintf(cons_f, "%s", seq1);
		fprintf( cons_f, "\n");
		fclose(cons_f);
		cons_f = fopen(file_name2,"w");
		fprintf(cons_f, ">%i\n", tmp2->prf_number);
		fprintf(cons_f, "%s", seq2);
		fprintf( cons_f, "\n");
		fclose(cons_f);
		num_diagonals = seq_pair2blastz_diagonal(file_name1, file_name2, diagonals_p, arguments->dig_length, strlen(seq1),strlen(seq2), is_dna);
	}
	else if (!strcmp(param_set->diag_method, "blast"))
	{
		FILE *cons_f = fopen(file_name1,"w");
		fprintf(cons_f, ">%i\n", tmp1->prf_number);
		fprintf(cons_f, "%s", seq1);
		fprintf( cons_f, "\n");
		fclose(cons_f);
		cons_f = fopen(file_name2,"w");
		fprintf(cons_f, ">%i\n", tmp2->prf_number);
		fprintf(cons_f, "%s", seq2);
		fprintf( cons_f, "\n");
		fclose(cons_f);
		int l1 = strlen(seq1);
		int l2 = strlen(seq2);
		num_diagonals = seq_pair2blast_diagonal(file_name1, file_name2, diagonals_p, arguments->dig_length, l1, l2, is_dna);
// 		int *num_p = &num_diagonals;
// 		Segment* seg = extend_diagonals(*diagonals_p, num_p, l1, l2);
// 		printf("A: %i\n", num_diagonals);
	}
	else if (!strcmp(param_set->diag_method, "ktup"))
	{
		num_diagonals = seq_pair2diagonal_own(seq1, seq2, diagonals_p, arguments->dig_length, strlen(seq1),strlen(seq2), is_dna, 3);
//	num_diagonals = seq_pair2diagonal_swift(seq1, seq2, diagonals_p, arguments->dig_length, strlen(seq1),strlen(seq2), is_dna, 3);
	}


// 	arguments->diagonals = diagonals_p[0];
// 	arguments->list = segments2int(seg, *num_p, seq1, seq2, profiles[0], profiles[1], arguments->list_length, param_set);

// t ** segments2int_gap(Segment *diagonals, int num_diagonals, char *seq1, char *seq2, Fastal_profile *profile1, Fastal_profile *profile2, int *num_points, Fastal_param *param_set);

//  arguments->list = diagonals2int_dot(arguments->diagonals, num_diagonals, seq1, seq2, profiles[0], profiles[1], arguments->list_length, param_set);
// 	arguments->list = diagonals2int_euclidf(arguments->diagonals, num_diagonals, seq1, seq2, profiles[0], profiles[1], arguments->list_length, param_set);
	arguments->list = diagonals2int_gap_test(arguments->diagonals, num_diagonals, seq1, seq2, profiles[0], profiles[1], arguments->list_length, param_set);
// 	arguments->list = diagonals2int(arguments->diagonals, num_diagonals, seq1, seq2, arguments->list_length, param_set);
	int alignment_length = list2linked_pair_wise_fastal(profiles[0], profiles[1], param_set, arguments->list, *arguments->list_length, edit_file, prof_file, number);
	int x;

	for (x = 0; x < *arguments->list_length; ++x)
	{
		vfree(arguments->list[x]);
	}
	vfree(arguments->list);
	arguments->list = NULL;
	vfree(seq1);
	vfree(seq2);
	return alignment_length;
}


int
fastal_compare (const void * a, const void * b)
{
	return (*(int*)a - *(int*)b);
}


/**
 * \brief Makes a sorted list out of diagonals.
 *
 * \param diagonals A list of diagonals to use during dynamic programming.
 * \param num_diagonals Number of diagonals.
 * \param seq1 Sequence 1.
 * \param seq2 Sequence 2.
 * \param num_points Number of points in the list
 * \param param_set Fastal parameters.
 * \return A 2-dim array which contains all points needed for the sparse dynamic programming algorithm.
 */
int **
diagonals2int(int *diagonals,
			  int num_diagonals,
			  char *seq1,
			  char *seq2,
			  int *num_points,
			  Fastal_param *param_set)
{

	int l1 = strlen(seq1);
	int l2 = strlen(seq2);
	int gep = param_set->gep;

	int dig_length;
	if (seq1 > seq2)
		dig_length = l1;
	else
		dig_length = l2;

	int current_size = num_diagonals*dig_length + l1 +l2;

	int **list = (int**)vcalloc(current_size, sizeof(int*));
	int *diags = (int*)vcalloc(num_diagonals, sizeof(int));
	int i;
	for (i = 0; i < num_diagonals; ++i)
	{
		diags[i] = l1 - diagonals[i*3] + diagonals[i*3+1];
	}

	qsort (diags, num_diagonals, sizeof(int), fastal_compare);


	int *diagx =(int*) vcalloc(num_diagonals, sizeof(int));
	int *diagy =(int*) vcalloc(num_diagonals, sizeof(int));


	//+1 because diagonals start here at position 1, like in "real" dynamic programming
	int a = -1, b = -1;
	for (i = 0; i < num_diagonals; ++i)
	{
		if (diags[i] < l1)
		{
			diagx[i] = l1 - diags[i];
			diagy[i] = 0;
			a= i;
		}
		else
			break;
	}
	++a;
	b=a-1;
	for (; i < num_diagonals; ++i)
	{
		diagx[i] = 0;
		diagy[i] = diags[i]-l1;
		b = i;
	}

	vfree(diags);
	int tmpy_pos;
	int tmpy_value;
	int **M = param_set->M;
	int *last_y = (int*)vcalloc(l2+1, sizeof(int));
	int *last_x = (int*)vcalloc(l1+1, sizeof(int));
	last_y[0] = 0;

	last_x[0] = 0;
	list[0] =(int*) vcalloc(6, sizeof(int));

	int list_pos = 1;
	int dig_num = l1;
	int tmp_l2 = l2 + 1;

	//left border
	for (; list_pos < tmp_l2; ++list_pos)
	{
		list[list_pos] =(int*) vcalloc(6, sizeof(int));
		list[list_pos][0] = 0;
		list[list_pos][1] = list_pos;
		last_y[list_pos] = list_pos;
		list[list_pos][2] = list_pos*gep;
		list[list_pos][4] = list_pos-1;
	}

	int pos_x = 0;
	int y;
	int tmp_l1 = l1-1;
	while (pos_x < tmp_l1)
	{
		if (list_pos + num_diagonals+2 > current_size)
		{
			current_size += num_diagonals*1000;
			list =(int**) vrealloc(list, current_size * sizeof(int*));
		}
		//upper border
		list[list_pos] = (int*)vcalloc(6, sizeof(int));
		list[list_pos][0] = ++pos_x;
		list[list_pos][1] = 0;
		list[list_pos][2] = pos_x * gep;
		list[list_pos][3] = last_y[0];
		tmpy_value = list_pos;
		tmpy_pos = 0;
		last_x[pos_x] = list_pos;
		++list_pos;

		//diagonals
		for (i = a; i <= b; ++i)
		{
			list[list_pos] =(int*) vcalloc(6, sizeof(int));

			list[list_pos][0] = ++diagx[i];

			list[list_pos][1] = ++diagy[i];
			list[list_pos][3] = last_y[diagy[i]];
			list[list_pos][4] = list_pos-1;
			list[list_pos][5] = last_y[diagy[i]-1];
			list[list_pos][2] = M[toupper(seq1[diagx[i]-1])-'A'][toupper(seq2[diagy[i]-1])-'A'];
			last_y[tmpy_pos] = tmpy_value;
			tmpy_value = list_pos;
			tmpy_pos = diagy[i];

			++list_pos;
		}
		last_y[tmpy_pos] = tmpy_value;


		//lower border
		if (list[list_pos-1][1] != l2)
		{
			list[list_pos] =(int*) vcalloc(6, sizeof(int));
			list[list_pos][0] = pos_x;
			list[list_pos][1] = l2;
			list[list_pos][3] = last_y[l2];

			list[list_pos][2] = -1000;
			list[list_pos][4] = list_pos-1;
			if (pos_x > l2)
				list[list_pos][5] = last_x[pos_x-l2];
			else
				list[list_pos][5] = l2-pos_x;
			last_y[l2] = list_pos;
			++list_pos;

		}


		if ((b >= 0) && (diagy[b] == l2))
			--b;

		if ((a >0) && (diagx[a-1] == pos_x))
			--a;
	}


	dig_num = -1;
	if (list_pos + l2+2 > current_size)
	{
		current_size += list_pos + l2 + 2;
		list =(int**) vrealloc(list, current_size * sizeof(int*));
	}


// 	right border
	list[list_pos] =(int*) vcalloc(6, sizeof(int));
	list[list_pos][0] = l1;
	list[list_pos][1] = 0;
	list[list_pos][3] = last_x[l1-1];
	list[list_pos][2] = -1000;
	++list_pos;



	for (i = 1; i <= l2; ++i)
	{
		list[list_pos] =(int*) vcalloc(6, sizeof(int));
		list[list_pos][0] = l1;
		list[list_pos][1] = i;
		list[list_pos][3] = last_y[i];
		list[list_pos][4] = list_pos-1;
		y = last_y[i-1];
		if ((list[y][0] == l1-1) && (list[y][1] == i-1))
		{
			list[list_pos][5] = y;
			list[list_pos][2] = M[toupper(seq1[l1-1])-'A'][toupper(seq2[i-1])-'A'];
		}
		else
		{
			if (i <= l1)
			{
				list[list_pos][5] = last_x[l1-i];
			}
			else
			{
				list[list_pos][5] = i-l1;
			}
			list[list_pos][2] = -1000;
		}
		++list_pos;
	}

	list[list_pos - l2][2] = -1000;

	*num_points = list_pos;
	vfree(diagx);
	vfree(diagy);


	return list;
}


/**
 * \brief Makes a sorted list out of diagonals.
 *
 * \param diagonals A list of diagonals to use during dynamic programming.
 * \param num_diagonals Number of diagonals.
 * \param seq1 Sequence 1.
 * \param seq2 Sequence 2.
 * \param num_points Number of points in the list
 * \param param_set Fastal parameters.
 * \return A 2-dim array which contains all points needed for the sparse dynamic programming algorithm.
 */
int **
diagonals2int_gap_test(int *diagonals, int num_diagonals, char *seq1, char *seq2, Fastal_profile *profile1, Fastal_profile *profile2, int *num_points, Fastal_param *param_set)
{
	int alphabet_size = param_set->alphabet_size;
	int l1 = strlen(seq1);
	int l2 = strlen(seq2);
	int gep = param_set->gep;

	int current_size = l2+l1;

	int **list =(int**) vcalloc(current_size, sizeof(int*));
	int *diags =(int*) vcalloc(num_diagonals, sizeof(int));
	int i;
	for (i = 0; i < num_diagonals; ++i)
	{
		diags[i] = l1 - diagonals[i*3] + diagonals[i*3+1];
	}
	qsort (diags, num_diagonals, sizeof(int), fastal_compare);


	int *diagx =(int*) vcalloc(num_diagonals, sizeof(int));
	int *diagy =(int*) vcalloc(num_diagonals, sizeof(int));
	int *old_pos =(int*) vcalloc(num_diagonals, sizeof(int));

	//+1 because diagonals start here at position 1, like in "real" dynamic programming
	int a = -1, b = -1;
	for (i = 0; i < num_diagonals; ++i)
	{

		if (diags[i] < l1)
		{
			diagx[i] = l1 - diags[i];
			diagy[i] = 0;
			a= i;
		}
		else
			break;
	}
	++a;
	b=a-1;
	for (; i < num_diagonals; ++i)
	{
		diagx[i] = 0;
		diagy[i] = diags[i]-l1;
		b = i;
	}

	vfree(diags);
	int tmpy_pos;
	int tmpy_value;
	int **M = param_set->M;
	int *last_y =(int*) vcalloc(l2+1, sizeof(int));
	int *last_x =(int*) vcalloc(l1+1, sizeof(int));
	last_y[0] = 0;

	last_x[0] = 0;
	list[0] =(int*) vcalloc(6, sizeof(int));

	int list_pos = 1;
	int dig_num = l1;
	int tmp_l2 = l2 + 1;

	//left border
	for (; list_pos < tmp_l2; ++list_pos)
	{
		list[list_pos] =(int*) vcalloc(6, sizeof(int));
		list[list_pos][0] = 0;
		list[list_pos][1] = list_pos;
		last_y[list_pos] = list_pos;
		list[list_pos][2] = list_pos*gep;
		list[list_pos][4] = list_pos-1;
	}

	int pos_x = 0;
// 	int diags_old = l2;

// 	int tmp = l1;
	int y;
	int tmp_l1 = l1-1;
	while (pos_x < tmp_l1)
	{
		if (list_pos + num_diagonals+2 > current_size)
		{
			current_size += num_diagonals*1000;
			list =(int**) vrealloc(list, current_size * sizeof(int*));
		}
		//upper border
		list[list_pos] =(int*) vcalloc(6, sizeof(int));
		list[list_pos][0] = ++pos_x;
		list[list_pos][1] = 0;
		list[list_pos][2] = pos_x * gep;
		list[list_pos][3] = last_y[0];
		tmpy_value = list_pos;
		tmpy_pos = 0;
		last_x[pos_x] = list_pos;
		++list_pos;

		//diagonals
		for (i = a; i <= b; ++i)
		{
			list[list_pos] =(int*) vcalloc(6, sizeof(int));

			list[list_pos][0] = ++diagx[i];

			list[list_pos][1] = ++diagy[i];
			list[list_pos][3] = last_y[diagy[i]];
			list[list_pos][4] = list_pos-1;
			list[list_pos][5] = last_y[diagy[i]-1];




//SIMPLEGAP
			int num_seq = profile1->number_of_sequences + profile2->number_of_sequences;
			double gap_num = 0;
			int char_c;
			for (char_c = 0; char_c < alphabet_size; ++char_c)
			{

				gap_num += profile1->prf[char_c][diagx[i]-1] + profile2->prf[char_c][diagy[i]-1];
			}

			gap_num /= num_seq;

			list[list_pos][2] = M[toupper(seq1[diagx[i]-1])-'A'][toupper(seq2[diagy[i]-1])-'A'] * gap_num;

// CLUSTAL
// 			int num_seq = profile1->number_of_sequences + profile2->number_of_sequences;
// 			double gap_num = 0;
// 			int char_c, char_c2;
// 			for (char_c = 0; char_c < alphabet_size; ++char_c)
// 				for (char_c2 = 0; char_c2 < alphabet_size; ++char_c2)
// 			{
// 				gap_num += (profile1->prf[char_c][diagx[i]-1]/profile1->number_of_sequences) * (profile2->prf[char_c][diagy[i]-1]/profile2->number_of_sequences) * M[param_set->pos2char[char_c]-'A'][param_set->pos2char[char_c2]-'A'];
// 			}
// 			list[list_pos][2] = gap_num;


			last_y[tmpy_pos] = tmpy_value;
			tmpy_value = list_pos;
			tmpy_pos = diagy[i];

			++list_pos;
		}
		last_y[tmpy_pos] = tmpy_value;


		//lower border
		if (list[list_pos-1][1] != l2)
		{
			list[list_pos] =(int*) vcalloc(6, sizeof(int));
			list[list_pos][0] = pos_x;
			list[list_pos][1] = l2;
			list[list_pos][3] = last_y[l2];

			list[list_pos][2] = -1000;
			list[list_pos][4] = list_pos-1;
			if (pos_x > l2)
				list[list_pos][5] = last_x[pos_x-l2];
			else
				list[list_pos][5] = l2-pos_x;
			last_y[l2] = list_pos;
			++list_pos;

		}


		if ((b >= 0) && (diagy[b] == l2))
			--b;

		if ((a >0) && (diagx[a-1] == pos_x))
			--a;
	}


	dig_num = -1;
	if (list_pos + l2+2 > current_size)
	{
		current_size += list_pos + l2 + 2;
		list =(int**) vrealloc(list, current_size * sizeof(int*));
	}


// 	right border
	list[list_pos] =(int*) vcalloc(6, sizeof(int));
	list[list_pos][0] = l1;
	list[list_pos][1] = 0;
	list[list_pos][3] = last_x[l1-1];
	list[list_pos][2] = -1000;
	++list_pos;



	for (i = 1; i <= l2; ++i)
	{
		list[list_pos] =(int*) vcalloc(6, sizeof(int));
		list[list_pos][0] = l1;
		list[list_pos][1] = i;
		list[list_pos][3] = last_y[i];
		list[list_pos][4] = list_pos-1;
		y = last_y[i-1];
		if ((list[y][0] == l1-1) && (list[y][1] == i-1))
		{
			list[list_pos][5] = y;
			int num_seq = profile1->number_of_sequences + profile2->number_of_sequences;
			double gap_num = 0;
			int char_c;
			for (char_c = 0; char_c < alphabet_size; ++char_c)
			{
				gap_num += profile1->prf[char_c][l1-1] + profile2->prf[char_c][i-1];
			}

			gap_num /= num_seq;

			list[list_pos][2] = M[toupper(seq1[l1-1])-'A'][toupper(seq2[i-1])-'A'] * gap_num;
		}
		else
		{
			if (i <= l1)
			{
				list[list_pos][5] = last_x[l1-i];
			}
			else
			{
				list[list_pos][5] = i-l1;
			}
			list[list_pos][2] = -1000;
		}
		++list_pos;
	}

	list[list_pos - l2][2] = -1000;

	*num_points = list_pos;
	vfree(diagx);
	vfree(diagy);
	vfree(old_pos);

	return list;
}


int **
diagonals2int_euclidf(int *diagonals, int num_diagonals, char *seq1, char *seq2, Fastal_profile *profile1, Fastal_profile *profile2, int *num_points, Fastal_param *param_set)
{
	int alphabet_size = param_set->alphabet_size;
	int l1 = strlen(seq1);
	int l2 = strlen(seq2);
	int gep = param_set->gep;

	int current_size = l2+l1;

	int **list =(int**) vcalloc(current_size, sizeof(int*));
	int *diags =(int*) vcalloc(num_diagonals, sizeof(int));
	int i;
	for (i = 0; i < num_diagonals; ++i)
	{
		diags[i] = l1 - diagonals[i*3] + diagonals[i*3+1];
	}

	qsort (diags, num_diagonals, sizeof(int), fastal_compare);


	int *diagx =(int*) vcalloc(num_diagonals, sizeof(int));
	int *diagy =(int*) vcalloc(num_diagonals, sizeof(int));
	int *old_pos =(int*) vcalloc(num_diagonals, sizeof(int));

	//+1 because diagonals start here at position 1, like in "real" dynamic programming
	int a = -1, b = -1;
	for (i = 0; i < num_diagonals; ++i)
	{

		if (diags[i] < l1)
		{
			diagx[i] = l1 - diags[i];
			diagy[i] = 0;
			a= i;
		}
		else
			break;
	}
	++a;
	b=a-1;
	for (; i < num_diagonals; ++i)
	{
		diagx[i] = 0;
		diagy[i] = diags[i]-l1;
		b = i;
	}

	vfree(diags);
	int tmpy_pos;
	int tmpy_value;
// 	int **M = param_set->M;
	int *last_y =(int*) vcalloc(l2+1, sizeof(int));
	int *last_x =(int*) vcalloc(l1+1, sizeof(int));
	last_y[0] = 0;

	last_x[0] = 0;
	list[0] =(int*) vcalloc(6, sizeof(int));

	int list_pos = 1;
	int dig_num = l1;
	int tmp_l2 = l2 + 1;

	//left border
	for (; list_pos < tmp_l2; ++list_pos)
	{
		list[list_pos] =(int*) vcalloc(6, sizeof(int));
		list[list_pos][0] = 0;
		list[list_pos][1] = list_pos;
		last_y[list_pos] = list_pos;
		list[list_pos][2] = list_pos*gep;
		list[list_pos][4] = list_pos-1;
	}

	int pos_x = 0;
// 	int diags_old = l2;

// 	int tmp = l1;
	int y;
	int tmp_l1 = l1-1;
	while (pos_x < tmp_l1)
	{
		if (list_pos + num_diagonals+2 > current_size)
		{
			current_size += num_diagonals*1000;
			list =(int**) vrealloc(list, current_size * sizeof(int*));
		}
		//upper border
		list[list_pos] =(int*) vcalloc(6, sizeof(int));
		list[list_pos][0] = ++pos_x;
		list[list_pos][1] = 0;
		list[list_pos][2] = pos_x * gep;
		list[list_pos][3] = last_y[0];
		tmpy_value = list_pos;
		tmpy_pos = 0;
		last_x[pos_x] = list_pos;
		++list_pos;

		//diagonals
		for (i = a; i <= b; ++i)
		{
			list[list_pos] =(int*) vcalloc(6, sizeof(int));

			list[list_pos][0] = ++diagx[i];

			list[list_pos][1] = ++diagy[i];
			list[list_pos][3] = last_y[diagy[i]];
			list[list_pos][4] = list_pos-1;
			list[list_pos][5] = last_y[diagy[i]-1];
			int char_c;
			double tmp_score = 0;
			double freq1, freq2;
			for (char_c = 0; char_c < alphabet_size; ++char_c)
			{
				freq1 = (double)profile1->prf[char_c][diagx[i]-1] / profile1->number_of_sequences;

				freq2 = (double)profile2->prf[char_c][diagy[i]-1] / profile2->number_of_sequences;

				tmp_score += ( freq1 - freq2) * (freq1 - freq2);
			}

			list[list_pos][2] = 10 - sqrt(tmp_score);

			last_y[tmpy_pos] = tmpy_value;
			tmpy_value = list_pos;
			tmpy_pos = diagy[i];

			++list_pos;
		}
		last_y[tmpy_pos] = tmpy_value;


		//lower border
		if (list[list_pos-1][1] != l2)
		{
			list[list_pos] =(int*) vcalloc(6, sizeof(int));
			list[list_pos][0] = pos_x;
			list[list_pos][1] = l2;
			list[list_pos][3] = last_y[l2];

			list[list_pos][2] = -1000;
			list[list_pos][4] = list_pos-1;
			if (pos_x > l2)
				list[list_pos][5] = last_x[pos_x-l2];
			else
				list[list_pos][5] = l2-pos_x;
			last_y[l2] = list_pos;
			++list_pos;

		}


		if ((b >= 0) && (diagy[b] == l2))
			--b;

		if ((a >0) && (diagx[a-1] == pos_x))
			--a;
	}


	dig_num = -1;
	if (list_pos + l2+2 > current_size)
	{
		current_size += list_pos + l2 + 2;
		list =(int**) vrealloc(list, current_size * sizeof(int*));
	}


// 	right border
	list[list_pos] =(int*) vcalloc(6, sizeof(int));
	list[list_pos][0] = l1;
	list[list_pos][1] = 0;
	list[list_pos][3] = last_x[l1-1];
	list[list_pos][2] = -1000;
	++list_pos;



	for (i = 1; i <= l2; ++i)
	{
		list[list_pos] =(int*) vcalloc(6, sizeof(int));
		list[list_pos][0] = l1;
		list[list_pos][1] = i;
		list[list_pos][3] = last_y[i];
		list[list_pos][4] = list_pos-1;
		y = last_y[i-1];
		if ((list[y][0] == l1-1) && (list[y][1] == i-1))
		{
			list[list_pos][5] = y;
			int char_c;
			int tmp_score = 0;
			double freq1, freq2;
			for (char_c = 0; char_c < alphabet_size; ++char_c)
			{
				freq1 = profile1->prf[char_c][l1-1] / profile1->number_of_sequences;
				freq2 = profile2->prf[char_c][i-1] / profile2->number_of_sequences;
				tmp_score += ( freq1 - freq2) * (freq2 - freq1);
			}
			list[list_pos][2] = 10 - sqrt(tmp_score);
// 			list[list_pos][2] = M[toupper(seq1[l1-1])-'A'][toupper(seq2[i-1])-'A'];
		}
		else
		{
			if (i <= l1)
			{
				list[list_pos][5] = last_x[l1-i];
			}
			else
			{
				list[list_pos][5] = i-l1;
			}
			list[list_pos][2] = -1000;
		}
		++list_pos;
	}

	list[list_pos - l2][2] = -1000;

	*num_points = list_pos;
	vfree(diagx);
	vfree(diagy);
	vfree(old_pos);

	return list;
}

int **
diagonals2int_dot(int *diagonals, int num_diagonals, char *seq1, char *seq2, Fastal_profile *profile1, Fastal_profile *profile2, int *num_points, Fastal_param *param_set)
{
	int alphabet_size = param_set->alphabet_size;
	int l1 = strlen(seq1);
	int l2 = strlen(seq2);
	int gep = param_set->gep;

	int current_size = l2+l1;

	int **list = (int**)vcalloc(current_size, sizeof(int*));
	int *diags =(int*) vcalloc(num_diagonals, sizeof(int));
	int i;
	for (i = 0; i < num_diagonals; ++i)
	{
		diags[i] = l1 - diagonals[i*3] + diagonals[i*3+1];
	}

	qsort (diags, num_diagonals, sizeof(int), fastal_compare);


	int *diagx =(int*) vcalloc(num_diagonals, sizeof(int));
	int *diagy =(int*) vcalloc(num_diagonals, sizeof(int));
	int *old_pos =(int*) vcalloc(num_diagonals, sizeof(int));

	//+1 because diagonals start here at position 1, like in "real" dynamic programming
	int a = -1, b = -1;
	for (i = 0; i < num_diagonals; ++i)
	{

		if (diags[i] < l1)
		{
			diagx[i] = l1 - diags[i];
			diagy[i] = 0;
			a= i;
		}
		else
			break;
	}
	++a;
	b=a-1;
	for (; i < num_diagonals; ++i)
	{
		diagx[i] = 0;
		diagy[i] = diags[i]-l1;
		b = i;
	}

	vfree(diags);
	int tmpy_pos;
	int tmpy_value;
// 	int **M = param_set->M;
	int *last_y =(int*) vcalloc(l2+1, sizeof(int));
	int *last_x =(int*) vcalloc(l1+1, sizeof(int));
	last_y[0] = 0;

	last_x[0] = 0;
	list[0] =(int*) vcalloc(6, sizeof(int));

	int list_pos = 1;
	int dig_num = l1;
	int tmp_l2 = l2 + 1;

	//left border
	for (; list_pos < tmp_l2; ++list_pos)
	{
		list[list_pos] =(int*) vcalloc(6, sizeof(int));
		list[list_pos][0] = 0;
		list[list_pos][1] = list_pos;
		last_y[list_pos] = list_pos;
		list[list_pos][2] = list_pos*gep;
		list[list_pos][4] = list_pos-1;
	}

	int pos_x = 0;
// 	int diags_old = l2;

// 	int tmp = l1;
	int y;
	int tmp_l1 = l1-1;
	while (pos_x < tmp_l1)
	{
		if (list_pos + num_diagonals+2 > current_size)
		{
			current_size += num_diagonals*1000;
			list =(int**) vrealloc(list, current_size * sizeof(int*));
		}
		//upper border
		list[list_pos] =(int*) vcalloc(6, sizeof(int));
		list[list_pos][0] = ++pos_x;
		list[list_pos][1] = 0;
		list[list_pos][2] = pos_x * gep;
		list[list_pos][3] = last_y[0];
		tmpy_value = list_pos;
		tmpy_pos = 0;
		last_x[pos_x] = list_pos;
		++list_pos;

		//diagonals
		for (i = a; i <= b; ++i)
		{
			list[list_pos] =(int*) vcalloc(6, sizeof(int));

			list[list_pos][0] = ++diagx[i];

			list[list_pos][1] = ++diagy[i];
			list[list_pos][3] = last_y[diagy[i]];
			list[list_pos][4] = list_pos-1;
			list[list_pos][5] = last_y[diagy[i]-1];
			int char_c;
			double tmp_score = 0;
			double freq1, freq2;
			for (char_c = 0; char_c < alphabet_size; ++char_c)
			{
				freq1 = (double)profile1->prf[char_c][diagx[i]-1] / profile1->number_of_sequences;

				freq2 = (double)profile2->prf[char_c][diagy[i]-1] / profile2->number_of_sequences;

				tmp_score += freq1 * freq2;
			}

			list[list_pos][2] = tmp_score * 10;

			last_y[tmpy_pos] = tmpy_value;
			tmpy_value = list_pos;
			tmpy_pos = diagy[i];

			++list_pos;
		}
		last_y[tmpy_pos] = tmpy_value;


		//lower border
		if (list[list_pos-1][1] != l2)
		{
			list[list_pos] =(int*) vcalloc(6, sizeof(int));
			list[list_pos][0] = pos_x;
			list[list_pos][1] = l2;
			list[list_pos][3] = last_y[l2];

			list[list_pos][2] = -1000;
			list[list_pos][4] = list_pos-1;
			if (pos_x > l2)
				list[list_pos][5] = last_x[pos_x-l2];
			else
				list[list_pos][5] = l2-pos_x;
			last_y[l2] = list_pos;
			++list_pos;

		}


		if ((b >= 0) && (diagy[b] == l2))
			--b;

		if ((a >0) && (diagx[a-1] == pos_x))
			--a;
	}


	dig_num = -1;
	if (list_pos + l2+2 > current_size)
	{
		current_size += list_pos + l2 + 2;
		list = (int**)vrealloc(list, current_size * sizeof(int*));
	}


// 	right border
	list[list_pos] =(int*) vcalloc(6, sizeof(int));
	list[list_pos][0] = l1;
	list[list_pos][1] = 0;
	list[list_pos][3] = last_x[l1-1];
	list[list_pos][2] = -1000;
	++list_pos;



	for (i = 1; i <= l2; ++i)
	{
		list[list_pos] =(int*) vcalloc(6, sizeof(int));
		list[list_pos][0] = l1;
		list[list_pos][1] = i;
		list[list_pos][3] = last_y[i];
		list[list_pos][4] = list_pos-1;
		y = last_y[i-1];
		if ((list[y][0] == l1-1) && (list[y][1] == i-1))
		{
			list[list_pos][5] = y;
			int char_c;
			int tmp_score = 0;
			double freq1, freq2;
			for (char_c = 0; char_c < alphabet_size; ++char_c)
			{
				freq1 = profile1->prf[char_c][l1-1] / profile1->number_of_sequences;
				freq2 = profile2->prf[char_c][i-1] / profile2->number_of_sequences;
				tmp_score += freq2 * freq1;
			}
			list[list_pos][2] = tmp_score * 10;
// 			list[list_pos][2] = M[toupper(seq1[l1-1])-'A'][toupper(seq2[i-1])-'A'];
		}
		else
		{
			if (i <= l1)
			{
				list[list_pos][5] = last_x[l1-i];
			}
			else
			{
				list[list_pos][5] = i-l1;
			}
			list[list_pos][2] = -1000;
		}
		++list_pos;
	}

	list[list_pos - l2][2] = -1000;

	*num_points = list_pos;
	vfree(diagx);
	vfree(diagy);
	vfree(old_pos);

	return list;
}


void
combine_profiles2file(int **prf1,
						int **prf2,
						int pos1,
						int pos2,
						Fastal_param *param_set,
						FILE *prof_f,
						char state)
{
	int alphabet_size = param_set->alphabet_size;
	char *pos2aa = &(param_set->pos2char[0]);
	int i;
	int x = 0;
	if (state == 'M')
	{
		for (i = 0; i < alphabet_size; ++i)
			if (prf1[i][pos1] + prf2[i][pos2] > 0)
			{
				if (x)
					fprintf(prof_f," %c%i", pos2aa[i],prf1[i][pos1]+prf2[i][pos2]);
				else
					fprintf(prof_f,"%c%i", pos2aa[i],prf1[i][pos1]+prf2[i][pos2]);
				x = 1;
			}
		fprintf(prof_f,"\n");
	}
	else if (state == 'D')
	{
		for (i = 0; i < alphabet_size; ++i)
			if (prf2[i][pos2] > 0)
		{
			if (x)
				fprintf(prof_f," %c%i", pos2aa[i],prf2[i][pos2]);
			else
				fprintf(prof_f,"%c%i", pos2aa[i],prf2[i][pos2]);
			x = 1;
		}
		fprintf(prof_f,"\n");
	}
	else
	{
		for (i = 0; i < alphabet_size; ++i)
			if (prf1[i][pos1] > 0)
		{
			if (x)
				fprintf(prof_f," %c%i", pos2aa[i],prf1[i][pos1]);
			else
				fprintf(prof_f,"%c%i", pos2aa[i],prf1[i][pos1]);
			x = 1;
		}
		fprintf(prof_f,"\n");
	}
}



#define LIN(a,b,c) a[b*5+c]
/**
 * Calculates a fast and sparse dynamic programming matrix
 *
 * \param prf1 Profile of first sequence.
 * \param prf2 Profile of second sequence.
 * \param param_set The parameter for the alignment.
 * \param list The list of diagonals.
 * \param n number of dots.
 * \param edit_f File to save the edit information.
 * \param prof_f File to save the profile.
 * \param node_number Number of the new profile.
 */
int
list2linked_pair_wise_fastal(Fastal_profile *prf1,
							 Fastal_profile *prf2,
							 Fastal_param *param_set,
							 int **list,
							 int n,
							 FILE *edit_f,
							 FILE *prof_f,
							 int node_number)
{
	int a,b, i, j, LEN=0, start_trace = -1;
	int pi, pj,ij, delta_i, delta_j, prev_i, prev_j;
// 	static int **slist;
	static long *MI, *MJ, *MM,*MT2;
// 	static int *sortseq;
	static int max_size;
	int gop, gep, igop, igep;
	int l1, l2, l, ls;
	char **al;
	int ni=0, nj=0;
	long score;
	int nomatch = param_set->nomatch;

	l1=prf1->length;
	l2=prf2->length;

	al=declare_char (2,l1+l2+1);



	igop=param_set->gop;
	gep=igep=param_set->gep;
	if (n>max_size)
	{
		max_size=n;

		vfree (MI);vfree (MJ); vfree (MM);

		MI=(long int*)vcalloc (5*n, sizeof (long));
		MJ=(long int*)vcalloc (5*n, sizeof (long));
		MM=(long int*)vcalloc (5*n, sizeof (long));

	}
	else
	{
		for (a=0; a<n; a++)
			for (b=0; b<5; b++)
				LIN(MI,a,b)=LIN(MJ,a,b)=LIN(MJ,a,b)=-1000000;
	}

	for (a=0; a<n; a++)
	{
		i=list[a][0];
		j=list[a][1];


		if (i==l1 || j==l2)gop=0;
		else gop=igop;

		if (i==l1 && j==l2)start_trace=a;
		else if ( i==0 || j==0)
		{
			LIN(MM,a,0)=-1000000;
			if (j==0)
			{
				LIN(MJ,a,0)=-10000000;
				LIN(MI,a,0)=gep*i;
			}
			else if (i==0)
			{
				LIN(MI,a,0)=-10000000;
				LIN(MJ,a,0)=gep*j;
			}

			LIN(MI,a,1)=LIN(MJ,a,1)=-1;
			LIN(MI,a,2)=LIN(MJ,a,2)=i;
			LIN(MI,a,3)=LIN(MJ,a,3)=j;
			continue;
		}

		pi = list[a][3];
		pj = list[a][4];
		ij = list[a][5];

		prev_i=list[pi][0];
		prev_j=list[pj][1];

		delta_i=list[a][0]-list[pi][0];
		delta_j=list[a][1]-list[pj][1];

		/*Linear Notation*/
		LIN(MI,a,0)=MAX(LIN(MI,pi,0),(LIN(MM,pi,0)+gop))+delta_i*gep;
		LIN(MI,a,1)=pi;
		LIN(MI,a,2)=delta_i;
		LIN(MI,a,3)=0;
		LIN(MI,a,4)=(LIN(MI,pi,0) >=(LIN(MM,pi,0)+gop))?'i':'m';

		LIN(MJ,a,0)=MAX(LIN(MJ,pj,0),(LIN(MM,pj,0)+gop))+delta_j*gep;
		LIN(MJ,a,1)=pj;
		LIN(MJ,a,2)=0;
		LIN(MJ,a,3)=delta_j;

		LIN(MJ,a,4)=(LIN(MJ,pj,0) >=LIN(MM,pj,0)+gop)?'j':'m';

		if (a>1 && (ls=list[a][0]-list[ij][0])==(list[a][1]-list[ij][1]))
		{
			LIN(MM,a,0)=MAX3(LIN(MM,ij,0),LIN(MI,ij,0),LIN(MJ,ij,0))+list[a][2]-(ls*nomatch);

			LIN(MM,a,1)=ij;
			LIN(MM,a,2)=ls;
			LIN(MM,a,3)=ls;
			if ( LIN(MM,ij,0) >=LIN(MI,ij,0) && LIN(MM,ij,0)>=LIN(MJ,ij,0))LIN(MM,a,4)='m';
			else if ( LIN(MI,ij,0) >= LIN(MJ,ij,0))LIN(MM,a,4)='i';
			else LIN(MM,a,4)='j';

		}
		else
		{
			LIN(MM,a,0)=UNDEFINED;
			LIN(MM,a,1)=-1;
		}
	}

	a=start_trace;
	if (LIN(MM,a,0)>=LIN(MI,a,0) && LIN(MM,a,0) >=LIN(MJ,a,0))MT2=MM;
	else if ( LIN(MI,a,0)>=LIN(MJ,a,0))MT2=MI;
	else MT2=MJ;

	score=MAX3(LIN(MM,a,0), LIN(MI,a,0), LIN(MJ,a,0));

	i=l1;
	j=l2;

	while (!(i==0 &&j==0))
	{
		int next_a;
		l=MAX(LIN(MT2,a,2),LIN(MT2,a,3));
      // HERE ("%c from %c %d %d SCORE=%d [%d %d] [%2d %2d]", T2[a][5],T2[a][4], T2[a][2], T2[a][3], T2[a][0], gop, gep, i, j);
		if (i==0)
		{
			while ( j>0)
			{
				al[0][LEN]=0;
				al[1][LEN]=1;
				j--; LEN++;
			}
		}
		else if (j==0)
		{
			while ( i>0)
			{
				al[0][LEN]=1;
				al[1][LEN]=0;
				i--; LEN++;
			}
		}

// 		else if (l==0) {HERE ("L=0 i=%d j=%d",l, i, j);exit (0);}
		else
		{
			for (b=0; b<l; b++, LEN++)
			{
				if (LIN(MT2,a,2)){al[0][LEN]=1;i--;ni++;}
				else al[0][LEN]=0;

				if (LIN(MT2,a,3)){al[1][LEN]=1;j--;nj++;}
				else al[1][LEN]=0;
			}

			next_a=LIN(MT2,a,1);
			if (LIN(MT2,a,4)=='m')MT2=MM;
			else if (LIN(MT2,a,4)=='i')MT2=MI;
			else if (LIN(MT2,a,4)=='j')MT2=MJ;
			a=next_a;
		}
	}

	invert_list_char ( al[0], LEN);
	invert_list_char ( al[1], LEN);

	fprintf(edit_f, "%i\n%i\n%i\n%i\n",prf1->prf_number, prf2->prf_number, prf1->is_leaf, prf2->is_leaf);
	fprintf(prof_f, "%i\n0\n%i\n1\n%i\n", node_number,LEN, prf1->number_of_sequences+prf2->number_of_sequences);

	char statec[] = {'M','D','I'};
	int num = 0;
	int state = 0;
	i = 0;
	j = 0;

	for ( b=0; b< LEN; b++)
	{
		if ((al[0][b]==1) && (al[1][b]==1))
		{

			combine_profiles2file(prf1->prf, prf2->prf, i, j, param_set, prof_f, 'M');
			++i;
			++j;
			if (state != 0)
			{
				fprintf(edit_f, "%c%i\n",statec[state], num);
				num =1;
				state = 0;
			}
			else
				++num;
		}
		else if (al[0][b]==1)
		{
// 			prf1->prf[param_set->alphabet_size-1] += prf2->num_sequences;
			combine_profiles2file(prf1->prf, prf2->prf, i, j, param_set, prof_f, 'I');
			++i;
			if (state != 2)
			{
				fprintf(edit_f, "%c%i\n",statec[state], num);
				num =1;
				state = 2;
			}
			else
				++num;
		}
		else if (al[1][b]==1)
		{
//			prf2->prf[param_set->alphabet_size-1] += prf1->num_sequences;
			combine_profiles2file(prf1->prf, prf2->prf, i, j, param_set, prof_f, 'D');
			++j;
			if (state != 1)
			{
				fprintf(edit_f, "%c%i\n",statec[state], num);
				num =1;
				state = 1;
			}
			else
				++num;
		}
	}


	fprintf(edit_f, "%c%i\n",statec[state], num);
	num =1;
	state = 1;


	fprintf(edit_f,"*\n");
	fprintf(prof_f,"*\n");
	free_char (al, -1);
//    exit(0);
	return LEN;
}






/**
 * \brief Turns a profile into a consensus sequence.
 *
 * The character with the highest number of occurences is used as consensus. Gaps are not included. For example: 10 '-' and one 'A' would give 'A' as consensus.
 * \param profile The profile.
 * \param file_name Name of the file to save the consensus sequence in.
 * \param param_set The parameter of the fastal algorithm.
 * \return the sequence
 */
char*
profile2consensus(Fastal_profile *profile, Fastal_param *param_set)
{

// 	FILE *cons_f = fopen(file_name,"w");
// 	fprintf(cons_f, ">%i\n", profile->prf_number);
	char* seq =(char*) vcalloc(profile->length+1, sizeof(char));
	int i, j;
	int most_pos = -1, most;
	int alphabet_size = param_set->alphabet_size;
	int **prf = profile->prf;
	char *pos2char = param_set->pos2char;
	for (i = 0; i < profile->length; ++i)
	{
		most = -1;
		for (j = 0; j < alphabet_size; ++j)
		{
			if (prf[j][i] > most)
			{
				most = prf[j][i];
				most_pos = j;
			}
		}
		seq[i] = pos2char[most_pos];
// 		fprintf(cons_f, "%c",pos2char[most_pos]);
	}
	return seq;
}


int
diag_compare (const void * a, const void * b)
{
	return (((Diagonal_counter*)b)->count - ((Diagonal_counter*)a)->count);
}

/**
 * \brief Calculates the diagonals between two sequences.
 *
 * Uses  to calculate the diagonals.
 * \param seq_file1 File with sequence 1.
 * \param seq_file2 File with sequence 2.
 * \param diagonals An array where the diagonal points will be stored.
 * \param dig_length length of \a diagonals .
 * \param num_points Number of points in all diagonals.
 * \return number of diagonals;
 */
int
seq_pair2diagonal_own(char *seq1,
					  char *seq2,
					  int **diagonals,
					  int *dig_length,
					  int l1,
					  int l2,
					  int is_dna,
					  int word_length)
{
// 	word_length = 7;
	int word_number, i;
	int ng;
	if (is_dna)
	{
		word_number = (int)pow(5, word_length);
		ng = 4;
	}
	else
	{
		word_number = (int)pow(24, word_length);
		ng = 24;
	}
	int **word_index =(int**) vcalloc(word_number, sizeof(int*));
	for (i = 0 ; i < word_number; ++i)
	{
		word_index[i] =(int*) vcalloc(20, sizeof(int));
		word_index[i][0] = 2;
		word_index[i][1] = 20;
	}


	//making of k-tup index of seq1

	int *prod=(int*)vcalloc (word_length, sizeof(int));
	for ( i=0; i<word_length; i++)
	{
		prod[word_length-i-1]=(int)pow(ng,i);
	}

	int aa[256];
	if (is_dna)
	{
		aa['A'] = 0;
		aa['C'] = 1;
		aa['G'] = 2;
		aa['T'] = 3;
		aa['U'] = 3;
	}
	else
	{
		aa['A'] = 0;
 		aa['B'] = 20;
		aa['C'] = 1;
		aa['D'] = 2;
		aa['E'] = 3;
		aa['F'] = 4;
		aa['G'] = 5;
		aa['H'] = 6;
		aa['I'] = 7;
		aa['J'] = 20;
		aa['K'] = 8;
		aa['L'] = 9;
		aa['M'] = 10;
		aa['N'] = 11;
		aa['P'] = 12;
		aa['Q'] = 13;
		aa['R'] = 14;
		aa['S'] = 15;
		aa['T'] = 16;
		aa['V'] = 17;
		aa['W'] = 18;
 		aa['X'] = 20;
		aa['Y'] = 19;
		aa['X'] = 20;
	}
	int index = 0;
	for (i = 0; i < word_length; ++i)
	{
		index += aa[(short)seq1[i]] *prod[i];
	}
	word_index[index][2] = 0;
	word_index[index][0] = 3;
	int z = -1;
	int *tmp;
	for (; i < l1; ++i)
	{
		index -= aa[(short)seq1[++z]] * prod[0];
		index *= ng;
		index += aa[(short)seq1[i]];
		tmp = word_index[index];
		if (tmp[0] == tmp[1])
		{
			tmp[1] += 25;
			tmp =(int*) vrealloc(tmp, word_index[index][1] *sizeof(int));
			word_index[index] = tmp;
		}
		tmp[tmp[0]++] = i;
	}



	//counting diagonals
	const int window_length = 14;

	Diagonal_counter *diag_index =(Diagonal_counter*) vcalloc(l1+l2, sizeof(Diagonal_counter));
	int num = l1+l2;
	for (i = 0; i < num; ++i)
	{
		diag_index[i].diagonal = i;
		diag_index[i].count = 0;
	}
	index = 0;

	int j;
	for (i = 0; i < word_length; ++i)
	{
		index += aa[(short)seq2[i]] *prod[i];
		for (j = 2; j < word_index[index][0]; ++j)
		{
			++(diag_index[i - word_index[index][j] + l1].count);
		}
	}

	z = -1;
	int i2 = i-1;
	int second_index = index;
	for (; i < window_length; ++i)
	{
		index -= aa[(short)seq2[++z]] * prod[0];
		index *= ng;
		index += aa[(short)seq2[i]];
		tmp = word_index[index];
		for (j = 2; j < tmp[0]; ++j)
		{
			++(diag_index[i - tmp[j] + l1].count);
		}
	}
	int z2 = -1;
	for (; i < l2; ++i)
	{
		index -= aa[(short)seq2[++z]] * prod[0];
		index *= ng;
		index += aa[(short)seq2[i]];
		second_index -= aa[(short)seq2[++z2]] * prod[0];
		second_index *= ng;
		second_index += aa[(short)seq2[++i2]];

		tmp = word_index[index];
		for (j = 2; j < tmp[0]; ++j)
		{
			++(diag_index[i - tmp[j] + l1].count);
		}


		tmp = word_index[second_index];
		for (j = 2; j < tmp[0]; ++j)
		{
			if (diag_index[i2 - tmp[j] + l1].count > window_length-3)
				diag_index[i2 - tmp[j] + l1].count = window_length+100;
			else
				--diag_index[i2 - tmp[j] + l1].count;
		}
	}


	//choose diagonals
	int *diags = diagonals[0];
	int current_pos = 0;


	qsort (diag_index, num, sizeof(Diagonal_counter*), diag_compare);

	i = 0;
	int y, x;
	while (diag_index[i].count > window_length+10)
	{
		if (current_pos > (*dig_length)-3)
		{
			(*dig_length) += 30;
			diags =(int*) vrealloc(diags, sizeof(int)*(*dig_length));
		}


		y = diag_index[i].diagonal - l1;
		if (y < 0)
		{
			x = y * (-1);
			y = 0;
		}
		else
		{
			x = 0;
		}
		diags[current_pos++] = x;
		diags[current_pos++] = y;
		diags[current_pos++] = 200;
		++i;
	}

	vfree(diag_index);
	for (i = 0; i < word_number; ++i)
		vfree(word_index[i]);
	vfree(word_index);
	diagonals[0] = diags;
	return current_pos/3;
}



int
seq_pair2diagonal_swift(char *seq1,
						char *seq2,
						int **diagonals,
						int *dig_length,
						int l1,
						int l2,
						int is_dna,
						int word_length)
{
	int word_number, i;
	int ng;
	if (is_dna)
	{
		word_number = (int)pow(5, word_length);
		ng = 5;
	}
	else
	{
		word_number = (int)pow(24, word_length);
		ng = 24;
	}
	int **word_index =(int**) vcalloc(word_number, sizeof(int*));
	for (i = 0 ; i < word_number; ++i)
	{
		word_index[i] =(int*) vcalloc(20, sizeof(int));
		word_index[i][0] = 2;
		word_index[i][1] = 20;
	}


	//making of k-tup index of seq1

	int *prod=(int*)vcalloc (word_length, sizeof(int));
	for ( i=0; i<word_length; i++)
	{
		prod[word_length-i-1]=(int)pow(ng,i);
	}

	int aa[256];
	aa['A'] = 0;
	aa['C'] = 1;
	aa['G'] = 2;
	aa['T'] = 3;
	aa['U'] = 4;
	int index = 0;
	for (i = 0; i < word_length; ++i)
	{
		index += aa[(short)seq1[i]] *prod[i];
	}
	word_index[index][2] = 0;
	word_index[index][0] = 3;
	int z = -1;
	int *tmp;
	for (; i < l1; ++i)
	{
		index -= aa[(short)seq1[++z]] * prod[0];
		index *= ng;
		index += aa[(short)seq1[i]];
		tmp = word_index[index];
		if (tmp[0] == tmp[1])
		{
			tmp[1] += 25;
			tmp =(int*) vrealloc(tmp, word_index[index][1] *sizeof(int));
			word_index[index] = tmp;
		}
		tmp[tmp[0]++] = i;
	}


	//counting diagonals
	const int window_length = 14;
	const int threshold = 12;

	Swift_diagonal *diag_index =(Swift_diagonal*) vcalloc(l1+l2, sizeof(Swift_diagonal));
	int num = l1+l2;
	for (i = 0; i < num; ++i)
	{
		diag_index[i].diagonal = i;
		diag_index[i].count = 0;
		diag_index[i].start = -99999;
		diag_index[i].end = -99999;
	}

	index = 0;

	int j;
	for (i = 0; i < word_length; ++i)
	{
		index += aa[(short)seq2[i]] *prod[i];
		for (j = 2; j < word_index[index][0]; ++j)
		{
			++(diag_index[i - word_index[index][j] + l1].count);
		}
	}

	z = -1;
	int tmp_index;
	for (; i < l2; ++i)
	{
		index -= aa[(short)seq2[++z]] * prod[0];
		index *= ng;
		index += aa[(short)seq2[i]];
		tmp = word_index[index];
		for (j = 2; j < tmp[0]; ++j)
		{
			tmp_index = i - tmp[j] + l1;
			if (i - diag_index[tmp_index].end > window_length)
			{
				if (diag_index[tmp_index].count < threshold)
				{
					diag_index[tmp_index].count = 0;
					diag_index[tmp_index].start = i;
					diag_index[tmp_index].end = i + word_length;
				}

			}
			else
			{
				++(diag_index[tmp_index].count);
			}
		}

	}



	// choose diagonals
	int *diags = diagonals[0];
	int current_pos = 0;
	int x, y;
	for (i = 0; i < num; ++i)
	{
		if (diag_index[i].count > threshold)
		{
			if (current_pos > (*dig_length)-3)
			{
				(*dig_length) += 30;
				diags =(int*) vrealloc(diags, sizeof(int)*(*dig_length));
			}
			y = diag_index[i].diagonal - l1;
			if (y < 0)
			{
				x = y * (-1);
				y = 0;
			}
			else
			{
				x = 0;
			}
			diags[current_pos++] = x;
			diags[current_pos++] = y;
			diags[current_pos++] = 200;
		}
	}

	vfree(diag_index);
	for (i = 0; i < word_number; ++i)
		vfree(word_index[i]);
	vfree(word_index);
	diagonals[0] = diags;

	return current_pos/3;
}



/**
 * \brief Calculates the diagonals between two sequences.
 *
 * Uses a k-tup index to choose diagonals.
 * \param seq_file1 File with sequence 1.
 * \param seq_file2 File with sequence 2.
 * \param diagonals An array where the diagonal points will be stored.
 * \param dig_length length of \a diagonals .
 * \param num_points Number of points in all diagonals.
 * \return number of diagonals;
 */
int
seq_pair2blast_diagonal(char *seq_file_name1,
						char *seq_file_name2,
						int **diagonals,
						int *dig_length,
						int l1,
						int l2,
						int is_dna)
{
// 	static int blast_measure[12]={0,0,0,0,0,0,0,0,0,0,0,0};
	int *diag = (int*) vcalloc(l1 + l2, sizeof(int));
	char *out_file = vtmpnam(NULL);
	char blast_command[200];
// 	char blast_command2[200];
// 	if (x)
// 	{
// 		int i;
// 		printf("BLAST-Types:\n");
// 		for (i = 0; i < 11; ++i)
// 		{
// 			printf("Type %i: %i\n", i, blast_measure[i]);
// 		}
// 		return 0;
// 	}
// 	char blast_command2[600];
// 	sprintf(blast_command2, "less %s", out_file);

	if (is_dna)
	{
		sprintf(blast_command, "bl2seq -p blastn -i %s -j %s -D 1 -g F -o %s -S 1 -F F", seq_file_name1, seq_file_name2, out_file);
	}
	else
	{
 		sprintf(blast_command, "bl2seq -p blastp -i %s -j %s -D 1 -g F -o %s -F F -S 1", seq_file_name1, seq_file_name2, out_file);
	}
	system(blast_command);

	int *diags = diagonals[0];
	FILE *diag_f = fopen(out_file,"r");
	char line[300];
	fgets(line, 300, diag_f);
	fgets(line, 300, diag_f);
	fgets(line, 300, diag_f);


	char delims[] = "\t";
	int length, pos_q, pos_d;
	int current_pos = 0;
	while (fgets(line, 300, diag_f) != NULL)
	{
		strtok(line, delims);
		strtok(NULL, delims);
		strtok(NULL, delims);
		length =  atoi(strtok(NULL, delims));
		strtok(NULL, delims);
		strtok(NULL, delims);
		pos_q = atoi(strtok(NULL, delims))-1;
		strtok(NULL, delims);
		pos_d = atoi(strtok(NULL, delims))-1;

		if (current_pos >= *dig_length-20)
		{
			(*dig_length) += 90;
			diags =(int*)  vrealloc(diags, sizeof(int)*(*dig_length));
		}
		if (diag[l1-(pos_q)+pos_d] == 0)
		{
			diag[l1-(pos_q)+pos_d] =1;
			diags[current_pos++] = pos_q;
			diags[current_pos++] = pos_d;
			diags[current_pos++] = length;
		}

	}
	fclose(diag_f);
	int round = 0;
	int e_threshold = 10;
	while ((current_pos == 0) && (round < 10))
	{
		if (is_dna)
		{
			sprintf(blast_command, "bl2seq -p blastn -i %s -j %s -D 1 -g F -o %s -S 1 -F F -W 6 -e %i", seq_file_name1, seq_file_name2, out_file, e_threshold);
		}
		else
		{
		  sprintf(blast_command, "bl2seq -p blastp -i %s -j %s -D 1 -g F -o %s -F F -S 1 -e %i", seq_file_name1, seq_file_name2, out_file, e_threshold);
		}
		system(blast_command);
		e_threshold *= 10;
		FILE *diag_f = fopen(out_file,"r");
		char line[300];
		fgets(line, 300, diag_f);
		fgets(line, 300, diag_f);
		fgets(line, 300, diag_f);


		char delims[] = "\t";
		while (fgets(line, 300, diag_f) != NULL)
		{
			strtok(line, delims);
			strtok(NULL, delims);
			strtok(NULL, delims);
			length =  atoi(strtok(NULL, delims));
			strtok(NULL, delims);
			strtok(NULL, delims);
			pos_q = atoi(strtok(NULL, delims))-1;
			strtok(NULL, delims);
			pos_d = atoi(strtok(NULL, delims))-1;

			if (current_pos >= *dig_length-20)
			{
				(*dig_length) += 90;
				diags =(int*) vrealloc(diags, sizeof(int)*(*dig_length));
			}
			if (diag[l1-(pos_q)+pos_d] == 0)
			{
				diag[l1-(pos_q)+pos_d] =1;
				diags[current_pos++] = pos_q;
				diags[current_pos++] = pos_d;
				diags[current_pos++] = length;
			}
		}
		fclose(diag_f);
		++round;
		if (current_pos < 27)
			current_pos = 0;
	}
// 	++blast_measure[round];

	if (current_pos == 0)
	{
		printf("BLAST NOT SUCCESFULL\n");
		if (l1 < l2)
		{
			int i;
			int diff = l2 - l1 + 10;
			for (i = diff; i > 0; --i)
			{
				if (current_pos >= *dig_length-20)
				{
					(*dig_length) += 90;
					diags = (int*) vrealloc(diags, sizeof(int)*(*dig_length));
				}
				diags[current_pos++] = 0;
				diags[current_pos++] = i;
				diags[current_pos++] = 100;
			}
			diff = 10;
// 			printf("A: %i\n", diff);
			for (i = 0; i < diff; ++i)
			{
				if (current_pos >= *dig_length-20)
				{
					(*dig_length) += 90;
					diags = (int*) vrealloc(diags, sizeof(int)*(*dig_length));
				}
				diags[current_pos++] = i;
				diags[current_pos++] = 0;
				diags[current_pos++] = 100;
			}
		}
		else
		{
			int i;
			int diff = 10;
// 			printf("A: %i\n", diff);
			for (i = diff; i > 0; --i)
			{
				if (current_pos >= *dig_length-20)
				{
					(*dig_length) += 90;
					diags = (int*) vrealloc(diags, sizeof(int)*(*dig_length));
				}
				diags[current_pos++] = 0;
				diags[current_pos++] = i;
				diags[current_pos++] = 100;
			}
			diff = l1 - l2 + 10;
			for (i = 0; i < diff; ++i)
			{
				if (current_pos >= *dig_length-20)
				{
					(*dig_length) += 90;
					diags = (int*) vrealloc(diags, sizeof(int)*(*dig_length));
				}
				diags[current_pos++] = i;
				diags[current_pos++] = 0;
				diags[current_pos++] = 100;
			}

		}
	}


	vfree(diag);

	diagonals[0] = diags;
	return current_pos/3;
}



int
seq_pair2blat_diagonal(char *seq_file_name1,
						char *seq_file_name2,
						int **diagonals,
						int *dig_length,
						int l1,
						int l2,
						int is_dna)
{
	int *diag = (int*) vcalloc(l1 + l2, sizeof(int));
	char *out_file = vtmpnam(NULL);
	char blast_command[200];
// 	char blast_command2[200];
// 	char blast_command2[600];
// 	sprintf(blast_command2, "less %s", out_file);

	if (is_dna)
	{
		sprintf(blast_command, "blat %s %s %s -out=blast8 -q=dna -t=dna -maxGap=0 >/dev/null 2>/dev/null", seq_file_name2, seq_file_name1, out_file);
	}
	else
	{
 		sprintf(blast_command, "blat %s %s %s -out=blast8 -prot -maxGap=0 >/dev/null 2>/dev/null", seq_file_name2, seq_file_name1, out_file);
	}
	system(blast_command);

	int *diags = diagonals[0];
	FILE *diag_f = fopen(out_file,"r");
	char line[300];
// 	fgets(line, 300, diag_f);
// 	fgets(line, 300, diag_f);
// 	fgets(line, 300, diag_f);


	char delims[] = "\t";
	int length, pos_q, pos_d;
	int current_pos = 0;
	while (fgets(line, 300, diag_f) != NULL)
	{
		strtok(line, delims);
		strtok(NULL, delims);
		strtok(NULL, delims);
		length =  atoi(strtok(NULL, delims));
		strtok(NULL, delims);
		strtok(NULL, delims);
		pos_q = atoi(strtok(NULL, delims))-1;
		strtok(NULL, delims);
		pos_d = atoi(strtok(NULL, delims))-1;

		if (current_pos >= *dig_length-20)
		{
			(*dig_length) += 90;
			diags =(int*) vrealloc(diags, sizeof(int)*(*dig_length));
		}
		if (diag[l1-(pos_q)+pos_d] == 0)
		{
			diag[l1-(pos_q)+pos_d] =1;
			diags[current_pos++] = pos_q;
			diags[current_pos++] = pos_d;
			diags[current_pos++] = length;
		}
	}
	if (current_pos == 0)
	{
		printf("BLAT NOT SUCCESFULL\n");
		if (l1 < l2)
		{
			int i;
			int diff = l2 - l1 + 10;
			for (i = diff; i > 0; --i)
			{
				if (current_pos >= *dig_length-20)
				{
					(*dig_length) += 90;
					diags = (int*) vrealloc(diags, sizeof(int)*(*dig_length));
				}
				diags[current_pos++] = 0;
				diags[current_pos++] = i;
				diags[current_pos++] = 100;
			}
			diff = 10;
// 			printf("A: %i\n", diff);
			for (i = 0; i < diff; ++i)
			{
				if (current_pos >= *dig_length-20)
				{
					(*dig_length) += 90;
					diags = (int*) vrealloc(diags, sizeof(int)*(*dig_length));
				}
				diags[current_pos++] = i;
				diags[current_pos++] = 0;
				diags[current_pos++] = 100;
			}
		}
		else
		{
			int i;
			int diff = 10;
// 			printf("A: %i\n", diff);
			for (i = diff; i > 0; --i)
			{
				if (current_pos >= *dig_length-20)
				{
					(*dig_length) += 90;
					diags = (int*) vrealloc(diags, sizeof(int)*(*dig_length));
				}
				diags[current_pos++] = 0;
				diags[current_pos++] = i;
				diags[current_pos++] = 100;
			}
			diff = l1 - l2 + 10;
			for (i = 0; i < diff; ++i)
			{
				if (current_pos >= *dig_length-20)
				{
					(*dig_length) += 90;
					diags = (int*) vrealloc(diags, sizeof(int)*(*dig_length));
				}
				diags[current_pos++] = i;
				diags[current_pos++] = 0;
				diags[current_pos++] = 100;
			}

		}
	}

// 	printf("END\n");
	vfree(diag);
	fclose(diag_f);
	diagonals[0] = diags;
	return current_pos/3;
}



/**
 * \brief Calculates the diagonals between two sequences.
 *
 * Uses blastz to calculate the diagonals.
 * \param seq_file1 File with sequence 1.
 * \param seq_file2 File with sequence 2.
 * \param diagonals An array where the diagonal points will be stored.
 * \param dig_length length of \a diagonals .
 * \param num_points Number of points in all diagonals.
 * \return number of diagonals;
 */
int
seq_pair2blastz_diagonal(char *seq_file_name1,
						 char *seq_file_name2,
						 int **diagonals,
						 int *dig_length,
						 int l1,
						 int l2,
						 int is_dna)
{
	int *diag = (int*) vcalloc(l1 + l2, sizeof(int));
	char *out_file = vtmpnam(NULL);
	char blast_command[200];
// 	char blast_command2[200];
// 	char blast_command2[600];
// 	sprintf(blast_command2, "less %s", out_file);

	if (is_dna)
	{
		sprintf(blast_command, "~/Download/blastz-source/blastz %s %s B=0 K=10000> %s", seq_file_name1, seq_file_name2, out_file);
	}
	else
	{
		printf("SORRY - no BLASTZ with amino accid\n");
		exit(0);
	}
	system(blast_command);

	int *diags = diagonals[0];
	FILE *diag_f = fopen(out_file,"r");
	char line[300];
	char delims[] = " ";
// 	char *result = NULL;
	int length, pos_q, pos_d;
	int current_pos = 0;
	while (fgets(line, 300, diag_f) != NULL)
	{
		if (line[0] == 'a')
		{
			char *line_tmp;
			while (fgets(line, 300, diag_f) != NULL)
			{
				if (line[0] == '}')
					break;

				if (line[2] == 'l')
				{
					line_tmp = &line[4];
					if (current_pos >= *dig_length-20)
					{
						(*dig_length) += 90;
						diags = (int*) vrealloc(diags, sizeof(int)*(*dig_length));
					}
					pos_q = atoi(strtok(line_tmp, delims));
					pos_d = atoi(strtok(NULL, delims));
					length = atoi(strtok(NULL, delims) - pos_q);
					if (diag[l1-(pos_q)+pos_d] == 0)
					{
						diag[l1-(pos_q)+pos_d] =1;
						diags[current_pos++] = pos_q;
						diags[current_pos++] = pos_d;
						diags[current_pos++] = length;
					}
				}
			}
		}
	}


	if (current_pos == 0)
	{
		printf("BLASTZ NOT SUCCESFULL\n");
		if (l1 < l2)
		{
			int i;
			int diff = l2 - l1 + 10;

			for (i = diff; i > 0; --i)
			{
				if (current_pos >= *dig_length-20)
				{
					(*dig_length) += 90;
					diags = (int*) vrealloc(diags, sizeof(int)*(*dig_length));
				}
				diags[current_pos++] = 0;
				diags[current_pos++] = i;
				diags[current_pos++] = 100;
			}
			diff = 10;

			for (i = 0; i < diff; ++i)
			{
				if (current_pos >= *dig_length-20)
				{
					(*dig_length) += 90;
					diags = (int*) vrealloc(diags, sizeof(int)*(*dig_length));
				}
				diags[current_pos++] = i;
				diags[current_pos++] = 0;
				diags[current_pos++] = 100;
			}
		}
		else
		{
			int i;
			int diff = 10;

			for (i = diff; i > 0; --i)
			{
				if (current_pos >= *dig_length-20)
				{
					(*dig_length) += 90;
					diags = (int*) vrealloc(diags, sizeof(int)*(*dig_length));
				}
				diags[current_pos++] = 0;
				diags[current_pos++] = i;
				diags[current_pos++] = 100;
			}
			diff = l1 - l2 + 10;
			for (i = 0; i < diff; ++i)
			{
				if (current_pos >= *dig_length-20)
				{
					(*dig_length) += 90;
					diags = (int*) vrealloc(diags, sizeof(int)*(*dig_length));
				}
				diags[current_pos++] = i;
				diags[current_pos++] = 0;
				diags[current_pos++] = 100;
			}

		}
	}

	vfree(diag);
	fclose(diag_f);
	diagonals[0] = diags;
	return current_pos/3;
}



//**************************   needleman-wunsch aligning **********************************************************


void fill_arguments_nw(Nw_param* method_arguments_p, int alphabet_size)
{
	method_arguments_p-> dyn_matrix =(double**) vcalloc(1,sizeof(double*));
	method_arguments_p->dyn_matrix[0] =(double*) vcalloc(1,sizeof(double));
	method_arguments_p->length1 =(int*) vcalloc(1,sizeof(int));
	method_arguments_p->length2 =(int*) vcalloc(1,sizeof(int));
	*method_arguments_p->length1 = 1;
	*method_arguments_p->length2 = 1;
	method_arguments_p->sumup_prf =(int**) vcalloc(alphabet_size+1,sizeof(int*));
	int i;
	for (i = 0; i < alphabet_size+1; ++i)
		method_arguments_p->sumup_prf[i] = (int*) vcalloc(1,sizeof(int));
	method_arguments_p->sumup_length = (int*) vcalloc(1,sizeof(int));
	*method_arguments_p->sumup_length = 1;
}


void free_nw(Nw_param* method_arguments_p, int alphabet_size)
{
	free_dyn_matrix(*method_arguments_p->length1,method_arguments_p->dyn_matrix);
	int i;
	for (i = 0; i <= alphabet_size; ++i)
	{
		vfree(method_arguments_p->sumup_prf[i]);
	}
	vfree(method_arguments_p->sumup_prf);
	vfree(method_arguments_p->length1);
	vfree(method_arguments_p->length2);
	vfree(method_arguments_p->sumup_length);
}


/**
 * \brief One run of needleman-wunsch dynamic programming.
 *
 * \param profiles The profiles.
 * \param param_set The fastal parameters.
 * \param method_arguments_p The method arguments.
 * \param is_dna Sequences are DNA (\a is_dna = 1) or protein.
 * \param edit_file The edit file.
 * \param prof_file the profile file.
 * \param number Number of the parent node.
 * \return The length of the alignment.
 */
int
nw_dyn(Fastal_profile **profiles, Fastal_param *param_set, void *method_arguments_p, int is_dna, FILE *edit_file, FILE *prof_file, int number)
{
	Nw_param *arguments = (Nw_param*)method_arguments_p;
// 	int old_length1 = *arguments->length1;
// 	int old_length2 = *arguments->length2;
	arguments->dyn_matrix = resize_dyn_matrix(arguments->dyn_matrix, *arguments->length1, *arguments->length2, profiles[0]->length+1, profiles[1]->length+1);
	*arguments->length1 = profiles[0]->length+1;
	*arguments->length2 = profiles[1]->length+1;
	int alignment_length = prf_nw(profiles[0], profiles[1], arguments->dyn_matrix, edit_file, arguments->sumup_prf, arguments->sumup_length, param_set);
	write2file(arguments->sumup_prf, alignment_length, prof_file, number,profiles[0]->number_of_sequences + profiles[1]->number_of_sequences, param_set);
	return alignment_length;
}


/**
 * \brief This method takes a profile and turns it into a sumed up version.
 *
 * Required for NW-algorithm.
 * \param profile The profile to sum up.
 * \param sumup A field where the result will be stored.
 * \param param_set Parameters for the fastal algorithm.
 * \return The new \a sumup.
 */
int**
sumup_profile(Fastal_profile *profile,
			  int **sumup,
			  Fastal_param *param_set)
{

	char *pos2aa = &(param_set->pos2char[0]);
	int alphabet_size = param_set->alphabet_size;
	int **M = param_set->M;
	int prof_length = profile->length;

	int i,j,k;

	for (i = 0; i < prof_length; ++i)
	{
		sumup[alphabet_size][i] = 0;
		for (k = 0; k < alphabet_size; ++k)
		{
			sumup[k][i] = 0;
			sumup[alphabet_size][i] += profile->prf[k][i];
			for (j = 0; j < alphabet_size; ++j)
			{
				sumup[k][i] += profile->weight * profile->prf[j][i] * M[pos2aa[j]-'A'][pos2aa[k]-'A'];
			}
		}
	}

	return sumup;
}


/**
 * \brief Turns the dynamic programming matrix into a editfile and calculates the new profile.
 *
 * Required for NW-algorithm.
 * \param prog_matrix The dynamic programming matrix.
 * \param prf1 The first profile.
 * \param prf2 The second profile.
 * \param edit_f A File object (already opened) to write the edit sequence into.
 * \param prf_field A 2-dim array to save the new profile into.
 * \param field_length Length of the new profile.
 * \param param_set Parameters of the Fastal-Algorithm
 */
int
nw_matrix2edit_file(double **prog_matrix,	//dynamic programming matrix
					Fastal_profile *prf1,	//profile of dim1
					Fastal_profile *prf2,	//profile of dim2
					FILE *edit_f,			//file to safe the edit in
					int **prf_field,		//space to safe the new profile
					int *field_length,
					Fastal_param *param_set)		//length of prf_field
{
// 	int **M = param_set->M;
	int alphabet_size = param_set->alphabet_size;
	double gap_cost = param_set -> gop;
	fprintf(edit_f, "%i\n%i\n%i\n%i\n",prf1->prf_number, prf2->prf_number, prf1->is_leaf, prf2->is_leaf);
	int sum[] = {0,0,0};
	char sumc[] = {'M','I','D'};
	int last = 0;
	int n = 0;
	int m = 0;
	int field_pos = 0;
	int i;
	int prf1_length = prf1->length;
	int prf2_length = prf2->length;
	while ((n < prf1_length) && (m < prf2_length))
	{
		//if necesarry allocate more memory for result
		if ((*field_length)-alphabet_size < field_pos)
		{
			(*field_length) += ENLARGEMENT_PER_STEP;

			for (i = 0; i <alphabet_size+1; ++i)
			{
				prf_field[i] = (int*)vrealloc(prf_field[i], (*field_length)*sizeof(int));
			}
		}

		if (prog_matrix[n][m] == (prog_matrix[n+1][m] +gap_cost))
		{
			for (i = 0; i<alphabet_size; ++i)
			{
				prf_field[i][field_pos] = prf1->prf[i][n];
			}
			++n;
			++ field_pos;

			if (last != 1)
			{
				fprintf(edit_f,"%c%i\n",sumc[last],sum[last]);
				sum[last] = 0;
			}
			last = 1;
			++sum[last];
		}
		else if (prog_matrix[n][m] == (prog_matrix[n][m+1] +gap_cost))
		{

			for (i = 0; i<alphabet_size; ++i)
			{
				prf_field[i][field_pos] = prf2->prf[i][m];
			}
			++m;
			++ field_pos;
			if (last != 2)
			{
				fprintf(edit_f,"%c%i\n",sumc[last],sum[last]);
				sum[last] = 0;
			}
			last = 2;
			++sum[last];
		}
		else
		{
			for (i = 0; i<alphabet_size; ++i)
			{
				prf_field[i][field_pos] = prf1->prf[i][n] + prf2->prf[i][m];
			}
			++n;
			++m;
			++ field_pos;
			if (last != 0)
			{
				fprintf(edit_f,"%c%i\n",sumc[last],sum[last]);
				sum[last] = 0;
			}
			last = 0;
			++sum[last];
		}
	}
	fprintf(edit_f,"%c%i\n",sumc[last],sum[last]);

	//gaps in prf2
	last = 0;
	while (n < prf1_length)
	{
		for (i = 0; i<alphabet_size; ++i)
		{
			prf_field[i][field_pos] = prf1->prf[i][n];
		}
		++n;
		++ field_pos;
		++last;
	}
	if (last > 0)
		fprintf(edit_f,"I%i\n",last);

	//gaps in prf1
	last = 0;
	while (m < prf2_length)
	{
		for (i = 0; i<alphabet_size; ++i)
		{
			prf_field[i][field_pos] = prf2->prf[i][m];
		}
		++m;
		++ field_pos;
		++last;
	}
	if (last > 0)
		fprintf(edit_f,"D%i\n",last);
	fprintf(edit_f,"*\n");
	return field_pos;
}



/**
 * \brief Pairwise alignments of profile is done here.
 *
 * \param profile1 Profile of sequence 1
 * \param profile2 Profile of sequence 2
 * \param prog_matrix Matrix for dynamic programming
 * \param edit_file_name The edit_file_name
 * \param sumup_prf The sumup version of profile 1, which later contains the aligned profile.
 * \param sumup_length Contains the length of the aligned profile.
 * \return length of the aligned profile
 */
int
prf_nw(Fastal_profile *profile1,
	   Fastal_profile *profile2,
	   double **prog_matrix,
	   FILE *edit_file_name,
	   int **sumup_prf,
	   int *sumup_length,
	   Fastal_param *param_set)
{
	int alphabet_size = param_set->alphabet_size;
	double gap_cost = param_set->gop;

	int i;
	if (*sumup_length < profile1->length)
	{
		for (i = 0; i < alphabet_size+1; ++i)
		{
			sumup_prf[i] = (int*)vrealloc(sumup_prf[i], profile1->length*sizeof(int));
		}
		*sumup_length = profile1->length;
	}
	sumup_prf = sumup_profile(profile1, sumup_prf, param_set);



	int j,k;
	int prof1_length = profile1->length;
	int prof2_length = profile2->length;

// 	int** M = param_set->M;
	double match_score;
// 	int amino_counter;
	int residue_pairs = 0;

	for (i = prof2_length; i > 0; --i)
	{
		prog_matrix[prof1_length][i] = gap_cost * (prof2_length-i);
	}

	i = prof1_length-1;
	prog_matrix[prof1_length][prof2_length] = 0.0;
	while (i >=0)
	{
		j = prof2_length-1;

		prog_matrix[i][prof2_length] = gap_cost*(prof1_length-i);
		while (j >=0)
		{
			match_score = 0.0;
			residue_pairs = 0;
			for (k = 0; k < alphabet_size; ++k)
			{
				residue_pairs += profile2->prf[k][j];
				match_score += (profile2->prf[k][j] * sumup_prf[k][i]);
			}
			match_score /= (residue_pairs * sumup_prf[alphabet_size][i]);
			prog_matrix[i][j] = MAX3(prog_matrix[i+1][j+1]+match_score, prog_matrix[i+1][j]+gap_cost, prog_matrix[i][j+1]+gap_cost);

			--j;
		}
		--i;
	}
	return nw_matrix2edit_file(prog_matrix, profile1, profile2, edit_file_name, sumup_prf, sumup_length, param_set);
}


/************** GOTOH ***********************/


void
fill_arguments_gotoh(Gotoh_param* method_arguments_p, int alphabet_size)
{
	method_arguments_p->m_matrix =(double**) vcalloc(1,sizeof(double*));
	method_arguments_p->m_matrix[0] = (double*)vcalloc(1,sizeof(double));
	method_arguments_p->d_matrix = (double**)vcalloc(1,sizeof(double*));
	method_arguments_p->d_matrix[0] = (double*)vcalloc(1,sizeof(double));
	method_arguments_p->i_matrix =(double**) vcalloc(1,sizeof(double*));
	method_arguments_p->i_matrix[0] = (double*)vcalloc(1,sizeof(double));
	method_arguments_p->length1 = (int*)vcalloc(1,sizeof(int));
	method_arguments_p->length2 =(int*) vcalloc(1,sizeof(int));
	method_arguments_p->log_saver =(double**) vcalloc(alphabet_size+1, sizeof(double*));
	*method_arguments_p->length1 = 1;
	*method_arguments_p->length2 = 1;
	method_arguments_p->sumup_prf =(int**) vcalloc(alphabet_size+1,sizeof(int*));
	int i;
	for (i = 0; i < alphabet_size+1; ++i)
	{
		method_arguments_p->sumup_prf[i] =(int*) vcalloc(1,sizeof(int));
		method_arguments_p->log_saver[i] =(double*) vcalloc(1, sizeof(double));
	}
	method_arguments_p->sumup_length =(int*) vcalloc(1,sizeof(int));
	*method_arguments_p->sumup_length = 1;
}


void
free_gotoh(Gotoh_param* method_arguments_p, int alphabet_size)
{
	free_dyn_matrix(*method_arguments_p->length1,method_arguments_p->m_matrix);
	free_dyn_matrix(*method_arguments_p->length1,method_arguments_p->i_matrix);
	free_dyn_matrix(*method_arguments_p->length1,method_arguments_p->d_matrix);

	int i;
	for (i = 0; i <= alphabet_size; ++i)
	{
		vfree(method_arguments_p->sumup_prf[i]);
	}
	vfree(method_arguments_p->sumup_prf);
	vfree(method_arguments_p->length1);
	vfree(method_arguments_p->length2);
	vfree(method_arguments_p->sumup_length);
}


int gotoh_dyn(Fastal_profile **profiles, Fastal_param *param_set, void *method_arguments_p, int is_dna, FILE *edit_file, FILE *prof_file, int number)
{
	Gotoh_param *arguments = (Gotoh_param*)method_arguments_p;
	arguments->m_matrix = resize_dyn_matrix(arguments->m_matrix, *arguments->length1, *arguments->length2, profiles[0]->length+1, profiles[1]->length+1);
	arguments->i_matrix = resize_dyn_matrix(arguments->i_matrix, *arguments->length1, *arguments->length2, profiles[0]->length+1, profiles[1]->length+1);
	arguments->d_matrix = resize_dyn_matrix(arguments->d_matrix, *arguments->length1, *arguments->length2, profiles[0]->length+1, profiles[1]->length+1);
	int i;
	if (profiles[1]->length > *arguments->length2-1)
	{
		for (i = 0; i < param_set->alphabet_size; ++i)
		{
			arguments->log_saver[i] =(double*) vrealloc(arguments->log_saver[i], (profiles[1]->length)*sizeof(double*));
		}
	}
	*arguments->length1 = profiles[0]->length+1;
	*arguments->length2 = profiles[1]->length+1;
	int alignment_length = prf_gotoh(profiles[0], profiles[1], edit_file, arguments, param_set);
	write2file(arguments->sumup_prf, alignment_length, prof_file, number, profiles[0]->number_of_sequences + profiles[1]->number_of_sequences, param_set);
	return alignment_length;
}


int
gotoh_matrix2edit_file(double **m_matrix,		//dynamic programming matrix
					   double **v_matrix,		//dynamic programming matrix
					   double **h_matrix,		//dynamic programming matrix
					   Fastal_profile *prf1,	//profile of dim1
					   Fastal_profile *prf2,	//profile of dim2
					   FILE *edit_f,			//file to safe the edit in
					   int **prf_field,			//space to safe the new profile
					   int *field_length,
					   Fastal_param *param_set)	//length of prf_field
{
	double comp_num = log((double)prf1->number_of_sequences) + log((double)prf2->number_of_sequences);
	int** M = param_set->M;
	int alphabet_size = param_set->alphabet_size;
	double gep = param_set -> gep;
	fprintf(edit_f, "%i\n%i\n%i\n%i\n",prf1->prf_number, prf2->prf_number, prf1->is_leaf, prf2->is_leaf);
	int sum[] = {0,0,0};
	char sumc[] = {'M','I','D'};
	int last = 0;
	int n = 0;
	int m = 0;
	int field_pos = 0;
	int i;
	int prf1_length = prf1->length;
	int prf2_length = prf2->length;
	int current_mode = 0;
	//determine start mode
	char *pos2char = param_set->pos2char;
	if (h_matrix[n][m] == m_matrix[n][m])
	{
		current_mode = 2;
	}
	else
	{

		if (v_matrix[n][m] == m_matrix[n][m])
		{
			current_mode = 1;
		}
		else
		{
			current_mode = 0;
		}
	}
// 	printf("%f %f %f - %i\n",h_matrix[n][m],v_matrix[n][m],m_matrix[n][m], current_mode);
	while ((n < prf1_length) && (m < prf2_length))
	{
		//if necesarry allocate more memory for result
		if ((*field_length)-alphabet_size < field_pos)
		{
			(*field_length) += ENLARGEMENT_PER_STEP;

			for (i = 0; i <alphabet_size+1; ++i)
			{
				prf_field[i] = (int*)vrealloc(prf_field[i], (*field_length)*sizeof(int));
			}
		}


		if (current_mode == 2)
		{
			for (i = 0; i<alphabet_size; ++i)
			{
				prf_field[i][field_pos] = prf2->prf[i][m];
			}
			if (h_matrix[n][m] != (h_matrix[n][m+1]+gep))
			{
				current_mode = 0;
			}
			++m;
			++ field_pos;
			if (last != 2)
			{
				fprintf(edit_f,"%c%i\n",sumc[last],sum[last]);
				sum[last] = 0;
			}
			last = 2;
			++sum[last];
		}
		else
		{
			if (current_mode== 1)
			{
				for (i = 0; i<alphabet_size; ++i)
				{
					prf_field[i][field_pos] = prf1->prf[i][n];
				}
				if (v_matrix[n][m] != (v_matrix[n+1][m]+gep))
				{
					current_mode = 0;
				}
				++n;
				++ field_pos;

				if (last != 1)
				{
					fprintf(edit_f,"%c%i\n",sumc[last],sum[last]);
					sum[last] = 0;
				}
				last = 1;
				++sum[last];
			}
			else
			{
				double match_score = 0.0;
				int char_c, char_c2;
				for (char_c = 0; char_c < alphabet_size; ++char_c)
				{
					for (char_c2 = 0; char_c2 < alphabet_size; ++char_c2)
					{

						if ((log(prf1->prf[char_c][n]) != -1) && ( log(prf2->prf[char_c2][m]) != -1))
						{
							match_score += exp(log((double)prf1->prf[char_c][n]) + log((double)prf2->prf[char_c2][m])-comp_num) * M[pos2char[char_c]-'A'][pos2char[char_c2]-'A'];
						}
					}
				}
				if (m_matrix[n+1][m+1] + match_score != m_matrix[n][m])
				{
					if (m_matrix[n][m] == v_matrix[n][m])
					{
						current_mode = 1;
						continue;
					}
					if (m_matrix[n][m] == h_matrix[n][m])
					{
						current_mode = 2;
						continue;
					}
				}
				for (i = 0; i<alphabet_size; ++i)
				{
					prf_field[i][field_pos] = prf1->prf[i][n] + prf2->prf[i][m];
				}
				++n;
				++m;
				++ field_pos;
				if (last != 0)
				{
					fprintf(edit_f,"%c%i\n",sumc[last],sum[last]);
					sum[last] = 0;
				}
				last = 0;
				++sum[last];
			}
		}

	}
	fprintf(edit_f,"%c%i\n",sumc[last],sum[last]);

	int needed = MAX(prf1_length -n, prf2_length -m);

	if ((*field_length) - needed -10 < field_pos)
	{
		(*field_length) += needed +10;

		for (i = 0; i <alphabet_size+1; ++i)
		{
			prf_field[i] = (int*)vrealloc(prf_field[i], (*field_length)*sizeof(int));
		}
	}
	//gaps in prf2
	last = 0;
	while (n < prf1_length)
	{
		for (i = 0; i<alphabet_size; ++i)
		{
			prf_field[i][field_pos] = prf1->prf[i][n];
		}
		++n;
		++ field_pos;
		++last;
	}
	if (last > 0)
		fprintf(edit_f,"I%i\n",last);

	//gaps in prf1
	last = 0;
	while (m < prf2_length)
	{
		for (i = 0; i<alphabet_size; ++i)
		{
			prf_field[i][field_pos] = prf2->prf[i][m];
		}
		++m;
		++ field_pos;
		++last;
	}
	if (last > 0)
		fprintf(edit_f,"D%i\n",last);

	fprintf(edit_f,"*\n");
	return field_pos;
}




/**
 * \brief The gotoh dynamic programming algorithm.
 *
 * \param profile1 The first profile.
 */
int
prf_gotoh(Fastal_profile *profile1,
		  Fastal_profile *profile2,
		  FILE *edit_file_name,
		  Gotoh_param *arguments,
		  Fastal_param *param_set)
{

// printf("I AM HERE - again\n");
	int **sumup_prf   = arguments->sumup_prf;
	int *sumup_length = arguments->sumup_length;
	int alphabet_size = param_set->alphabet_size;
	double gop = param_set->gop;
	double gep = param_set->gep;

	const int INF = -999999;
	int i;

	double **m_matrix = arguments->m_matrix;
	double **h_matrix = arguments->i_matrix;
	double **v_matrix = arguments->d_matrix;

	int j;
	int prof1_length = profile1->length;
	int prof2_length = profile2->length;

	int** M = param_set->M;
	double match_score;
	for (i = prof2_length; i >= 0; --i)
	{
		m_matrix[prof1_length][i] = gop + gep * (prof2_length-i);
		v_matrix[prof1_length][i] = INF;
		h_matrix[prof1_length][i] = m_matrix[prof1_length][i];
	}

	m_matrix[prof1_length][prof2_length] = 0.0;
	h_matrix[prof1_length][prof2_length] = INF;
	v_matrix[prof1_length][prof2_length] = INF;
	int l;
	double comp_num = log((double)profile1->number_of_sequences) + log((double)profile2->number_of_sequences);
	static double *log_test = NULL;
	if (!log_test)
		log_test =(double*) vcalloc(alphabet_size, sizeof(double));
// 	int k;
	int **prf1 = profile1->prf;
	int **prf2 = profile2->prf;
	double **log_test2 = arguments->log_saver;
	for (l = 0; l < alphabet_size; ++l)
	{
		for (i = 0; i < profile2->length; ++i)
		{
			if (prf2[l][i] > 0)
			{
				log_test2[l][i] = log(prf2[l][i]);
			}
			else
				log_test2[l][i] = -1;
		}
	}

	char *pos2char = param_set->pos2char;
	i = prof1_length-1;
	while (i >=0)
	{
		j = prof2_length-1;

		for (l = 0; l < alphabet_size; ++l)
		{
			if (prf1[l][i] > 0)
				log_test[l] = log((double)prf1[l][i]);
			else
				log_test[l] = -1;
		}
		m_matrix[i][prof2_length] = gop + gep *(prof1_length-i);
		v_matrix[i][prof2_length] = m_matrix[i][prof2_length];
		h_matrix[i][prof2_length] = INF;
		while (j >=0)
		{

			match_score = 0.0;
			v_matrix[i][j] = (MAX(v_matrix[i+1][j], m_matrix[i+1][j] + gop) + gep);
			h_matrix[i][j] = (MAX(h_matrix[i][j+1], m_matrix[i][j+1] + gop) + gep);

			int char_c, char_c2;
			int num = 0;
			for (char_c = 0; char_c < alphabet_size; ++char_c)
			{
				for (char_c2 = 0; char_c2 < alphabet_size; ++char_c2)
				{

					if ((log_test[char_c] != -1) && (log_test2[char_c2][j] != -1))
					{
						match_score += exp(log_test[char_c] + log_test2[char_c2][j]-comp_num) * M[pos2char[char_c]-'A'][pos2char[char_c2]-'A'];
					}
				}
			}

			m_matrix[i][j] = m_matrix[i+1][j+1]+match_score;

			if (m_matrix[i][j] < v_matrix[i][j])
			{
				m_matrix[i][j] = v_matrix[i][j];
			}
			if (m_matrix[i][j] < h_matrix[i][j])
			{
				m_matrix[i][j] = h_matrix[i][j];
			}

			--j;
		}
		--i;
	}
	return gotoh_matrix2edit_file(m_matrix, v_matrix, h_matrix, profile1, profile2, edit_file_name, sumup_prf, sumup_length, param_set);
}


/************* OTHER STUFF ******************/

/**
 * \brief Writes the alignment into the profile file and the edit file.
 *
 * \param profiles The two profiles to combine.
 * \param alignment The alinment information.
 * \param alignment The length of the alignment.
 * \param edit_f The edit file.
 * \param prof_f The profile file.
 * \param node_number the new node number.
 */
void
alignment2files(Fastal_profile **profiles,
				Fastal_param *param_set,
				int **alignment,
				int alignment_length,
				FILE *edit_f,
				FILE *prof_f,
				int node_number)
{
	fprintf(edit_f, "%i\n%i\n%i\n%i\n",profiles[0]->prf_number, profiles[1]->prf_number, profiles[0]->is_leaf, profiles[1]->is_leaf);
	fprintf(prof_f, "%i\n0\n%i\n1\n", node_number, alignment_length);

	int **prf1 = profiles[0]->prf;
	int **prf2 = profiles[1]->prf;
	int i = 0;
	int pos = 0;
	int pos1, pos2;

	char statec[] = {'M','D','I'};
	int num = 0;
	int state = 0;

	while (i < alignment_length)
	{

		pos1 = alignment[0][pos];
		pos2 = alignment[1][pos];
		// match
		if ((pos1 != -1) && (pos2 != -1))
		{

			combine_profiles2file(prf1, prf2, pos1, pos2, param_set, prof_f, 'M');
			if (state != 0)
			{
				fprintf(edit_f, "%c%i\n",statec[state], num);
				num =1;
				state = 0;
			}
			else
				++num;
			++i;
		}
		// insertion in seq 1
		else if (pos1 != -1)
		{
			combine_profiles2file(prf1, prf2, pos1, pos2, param_set, prof_f, 'I');
			if (state != 2)
			{
				fprintf(edit_f, "%c%i\n",statec[state], num);
				num =1;
				state = 2;
			}
			else
				++num;
			++i;
		}
		// deletion in seq 1
		else if (pos2 != -1)
		{
			combine_profiles2file(prf1, prf2, pos1, pos2, param_set, prof_f, 'D');
			if (state != 1)
			{
				fprintf(edit_f, "%c%i\n",statec[state], num);
				num =1;
				state = 1;
			}
			else
				++num;
			++i;
		}
		++pos;
	}
	fprintf(edit_f, "%c%i\n",statec[state], num);

	fprintf(edit_f,"*\n");
	fprintf(prof_f,"*\n");

}


//******************************* OTHER STUFF ***********************

/**
 *	\brief Reads the sequence from a given position in a fasta file and turns it into a profile.
 *
 * \param seq_file The file where the sequence is stored.
 * \param off_set The off_set from the beginning of the file to the position of the sequence name.
 * \param profile The profile where the sequence will be stored into.
 * \param prf_number The number of this profile.
 */
void
file_pos2profile(FILE *seq_file,			//File with sequences
				 long off_set,				//offset of sequence from the beginning of file point to the sequence name, not to the sequence itself
				 Fastal_profile *profile,	//profile to save into
				 int prf_number,			//number of the profile
				 Fastal_param *param_set)
{
	int alphabet_size = param_set->alphabet_size;
	profile->is_leaf = 1;
	profile->number_of_sequences = 1;
	int *aa2pos = &(param_set->char2pos[0]);
	const int LINE_LENGTH = 500;
 	char line[LINE_LENGTH];
	profile->num_sequences = 1;
	profile->prf_number = prf_number;
	fseek (seq_file , off_set , SEEK_SET );

	fgets (line, LINE_LENGTH , seq_file);
	int seq_length = 0;
	int i, j, x;

	while(fgets(line, LINE_LENGTH, seq_file)!=NULL)
	{
		if (line[0] != '>')
		{
			line[LINE_LENGTH-1] = '\n';
			if (seq_length + LINE_LENGTH >= profile->allocated_memory)
			{
				for (i = 0; i < alphabet_size; ++i)
				{
					profile->prf[i] =(int*) vrealloc(profile->prf[i], (profile->allocated_memory+PROFILE_ENLARGEMENT)*sizeof(int));
				}
				profile->allocated_memory += PROFILE_ENLARGEMENT;
			}

			i = 0;
			x = 0;
			while ((line[i] != '\n') && (line[i] != '\0'))
			{
				if (line[i] != '-')
				{
					for(j = 0; j<alphabet_size; ++j )
						profile->prf[j][seq_length+x] = 0;
					profile->prf[aa2pos[toupper((short)line[i])]][seq_length+x] = 1;
					++x;
				}
				++i;
			}
			seq_length += x;

		}
		else
			break;
	}
	profile->length = seq_length;

}



/**
 * \brief Constructs index of fasta_file.
 *
 * The index is of length n (n= number of sequences in the given multi fasta file.). In the order of appearance in the file the position of each sequence in the file is stored.
 * \param file_name The file with the sequences.
 * \param file_positions Array to save the positions in.
 * \return The number of sequences in \a file_name.
 */
int
make_index_of_file(char *file_name, 		//file with sequences
				   long **file_positions)	//array to save the positions
{
	const int LINE_LENGTH = 150;
	(*file_positions) =(long int*) vcalloc(ENLARGEMENT_PER_STEP,  sizeof(long));

	FILE *file = fopen(file_name,"r");

	char line[LINE_LENGTH];

	int num_of_sequences = 0;
	int mem_for_pos = ENLARGEMENT_PER_STEP;


	if (file == NULL)
	{
		printf("FILE NOT FOUND\n");
		exit(1);
	}
	else
	{
		(*file_positions)[num_of_sequences] = ftell(file);
		while(fgets(line, LINE_LENGTH , file)!=NULL)
		{
// 			int length = strlen(line);
			if (line[0] == '>')
			{
				++num_of_sequences;

				if (num_of_sequences == mem_for_pos)
				{
					(*file_positions) = (long int*)vrealloc((*file_positions),(ENLARGEMENT_PER_STEP+mem_for_pos) * sizeof(long));
					mem_for_pos += ENLARGEMENT_PER_STEP;
				}
			}
			(*file_positions)[num_of_sequences] = ftell(file);
		}
	}

	fclose(file);
	return num_of_sequences;
}



/**
 * \brief Reads a profile from a profile file.
 *
 * \param prof A Fastal_profile object to save the profile in.
 * \param profile_f file where the profile is stored.
 * \param position Position of the profile in \a profile_f.
 * \param param_set The parameter set for Fastal
 */

void
profile_file2profile(Fastal_profile *prof,	//structure to save the profile in
					 FILE *profile_f,		//file where the profile is stored
					 long position,			//position in profile_f where the profile is stored
					 Fastal_param *param_set)
{

	int alphabet_size = param_set->alphabet_size;

	int *aa2pos = &(param_set->char2pos[0]);


	fseek(profile_f,position,SEEK_SET);
	const int LINE_LENGTH = 500;
	char line[500];

	fgets(line, LINE_LENGTH, profile_f);

	prof->prf_number = atoi(line);
	fgets(line, LINE_LENGTH, profile_f);
	prof->is_leaf = atoi(line);

	fgets(line, LINE_LENGTH, profile_f);
	prof->length = atoi(line);
	fgets(line, LINE_LENGTH, profile_f);
	prof->weight = atoi(line);
	fgets(line, LINE_LENGTH, profile_f);
	prof->number_of_sequences = atoi(line);
	int i,j;
	if (prof->length > prof->allocated_memory)
		for (i = 0;i < alphabet_size; ++i)
		{
			prof->prf[i] =(int*) vrealloc(prof->prf[i],prof->length*sizeof(int));
		}
	prof->allocated_memory = prof->length;
	char delims[] = " ";
	char *result = NULL;
	char *result_num = NULL;

	int length = prof->length;

	for (i = 0; i < length; ++i)
	{
		for(j = 0; j<alphabet_size; ++j )
			prof->prf[j][i] = 0;
		fgets(line, LINE_LENGTH , profile_f);
		result = strtok( line, delims );

		while( result != NULL)
		{
			result_num = &result[1];
			prof->prf[aa2pos[(short)result[0]]][i] = atoi(result_num);
			result = strtok( NULL, delims );
		}
	}
}



/**
 * \brief Writes a profile into a file
 *
 * \param profile Pointer to the profile which has to be saved.
 * \param file A File object (already opened) to write the profile to.
 * \param param_set The parameters for the fastal algorithm.
 */
void
profile2file(Fastal_profile *profile,	//the profile to save
			 FILE* file,				//file to save in
			 Fastal_param *param_set)
{
	int alphabet_size = param_set->alphabet_size;

	char *pos2aa = &(param_set->pos2char[0]);

	fseek(file,0,SEEK_SET);

	fprintf(file,"%i\n", profile->prf_number);


	fprintf(file,"%i\n", profile->is_leaf);
	fprintf(file,"%i\n", profile->length);
	fprintf(file,"%i\n", profile->weight);
	int i = 0, j = 0;
	int max = profile->length;
	int x= 0;
	--alphabet_size;
	while (i < max)
	{
		for (j = 0; j < alphabet_size; ++j)
			if (profile->prf[j][i] > 0)
			{
				if (x)
					fprintf(file," %c%i", pos2aa[j],profile->prf[j][i]);
				else
					fprintf(file,"%c%i", pos2aa[j],profile->prf[j][i]);
					x = 1;
			}
		if (profile->prf[j][i] > 0)
		{
			if (x)
				fprintf(file," %c%i", pos2aa[j],profile->prf[j][i]);
			else
				fprintf(file,"%c%i", pos2aa[j],profile->prf[j][i]);
			x = 1;
		}
		x = 0;
		fprintf(file,"\n");
		++i;
	}
	fprintf(file,"*\n");
}



/**
*	Reads the profile out of an alignment (NOT IN USE)
*/
// void
// file2profile(FILE* profile_f,			//file to read the profile of
// 			 Fastal_profile *prof,		//profile saved in here
// 			 int prf_number,			//number of the profile
// 			 Fastal_param *param_set)
// {
// 	int alphabet_size = param_set->alphabet_size;
//
// 	int *aa2pos =  &(param_set->char2pos[0]);
//
//
// 	fseek(profile_f,0,SEEK_SET);
// 	const int LINE_LENGTH = 500;
// 	char line[500];
//
// 	fgets(line, LINE_LENGTH, profile_f);
// 	prof->prf_number = atoi(line);
// 	fgets(line, LINE_LENGTH, profile_f);
// 	prof->is_leaf = atoi(line);
//
// 	fgets(line, LINE_LENGTH, profile_f);
// 	prof->length = atoi(line);
//
// 	fgets(line, LINE_LENGTH, profile_f);
// 	prof->weight = atoi(line);
// 	int i,j;
// 	if (prof->length > prof->allocated_memory)
// 		for (i = 0;i < alphabet_size; ++i)
// 		{
// 			prof->prf[i] = vrealloc(prof->prf[i],prof->length*sizeof(int));
// 		}
//
// 	char delims[] = " ";
// 	char *result = NULL;
// 	char *result_num = NULL;
//
// 	int length = prof->length;
//
// 	for (i = 0; i < length; ++i)
// 	{
// 		for(j = 0; j<alphabet_size; ++j )
// 			prof->prf[j][i] = 0;
// 		fgets(line, LINE_LENGTH , profile_f);
// 		result = strtok( line, delims );
//
// 		while( result != NULL)
// 		{
// 			result_num = &result[1];
// 			prof->prf[aa2pos[(short)result[0]]][i] = atoi(result_num);
// 			result = strtok( NULL, delims );
// 		}
// 	}
// }






/**
 * \brief Writes the sequence into the alignment_file.
 *
 * \param aligned_sequence Pattern of aligned sequence.
 * \param sequence_file File with sequences.
 * \param sequence_position Positions of sequences in \a sequence_file.
 * \param alignment_file The file to write the sequence into.
 *
*/
void
edit_seq2aligned_seq(char *aligned_sequence,	//pattern for aligned sequence
					 FILE *sequence_file,		//file with all the sequences
					 long sequence_position,	//position in sequence file with the correct sequence
					 FILE *alignment_file)		//file to write the alignment into
{
	fseek(sequence_file, sequence_position, SEEK_SET);
	const int LINE_LENGTH = 300;
	char line[LINE_LENGTH];
	fgets (line, LINE_LENGTH , sequence_file);
	fprintf(alignment_file,"%s", line);	//writing of sequence name
	int pos = 0;
	int i = 0;
	while(fgets(line, LINE_LENGTH, sequence_file)!=NULL)
	{
		if (line[0] != '>')
		{

			line[LINE_LENGTH-1] = '\n';
			i = 0;
			while ((line[i] != '\n') && (line[i] != '\0'))
			{
				while (aligned_sequence[pos] == '-')
				{

					fprintf(alignment_file,"-");
					++pos;
				}
				if (line[i] != '-')
				{

					fprintf(alignment_file,"%c",line[i]);
					++pos;
				}
				++i;
// 				++pos;
			}
		}
		else
			break;
	}
	while (aligned_sequence[pos] != '\n')
	{
		fprintf(alignment_file,"-");
		++pos;
	}
	fprintf(alignment_file,"\n");
}



/**
 * \brief Recursive function to turn the edit_file into the alignment.
 *
 * \param sequence_file File with all sequences.
 * \param sequence_position The array of sequence positions in \a sequence_file
 * \param edit_file File to safe the edit profiles in.
 * \param edit_positions Array saving the coorespondence between edit profile and position in \a edit_file
 * \param node_number The current node.
 * \param number_of_sequences The number of sequences.
 * \param aligned_sequence The sequence that is edited.
 * \param alignment_length The length of the alignment.
 * \param edit_seq_file File that saves the edited_sequences of the internal nodes.
 * \param offset Saves the size of the edited_sequences.
 * \param alignment_file File where the alignment is saved.
 */
void
edit2alignment(FILE *sequence_file,		//sequence file
			   long *seq_positions,		//sequence positions
			   FILE *edit_file,			//file saving the edit profiles
			   long *edit_positions,	//array saving the correspondence between edit profile and position in edit_file
			   int node_number,			//the current node
			   int number_of_sequences,	//number of sequences
			   char *aligned_sequence,	//the sequence that is edited
			   int alignment_length,	//length of the alignment - and thus of aligned_sequence
			   FILE *edit_seq_file,		//file saving the edited_sequences of the internal nodes
			   int offset,				//saves the size of the edited_sequence
			   FILE* alignment_file)	//file saving the alignments
{
	fseek(edit_file, edit_positions[node_number-number_of_sequences], SEEK_SET);
	const int LINE_LENGTH = 50;
	char line[LINE_LENGTH];
	fgets(line, LINE_LENGTH , edit_file);
	int child1 = atoi(line);
	fgets(line, LINE_LENGTH , edit_file);
	int child2 = atoi(line);
	fgets(line, LINE_LENGTH , edit_file);
	int is_leaf1 = atoi(line);
	fgets(line, LINE_LENGTH , edit_file);
	int is_leaf2 = atoi(line);

// 	static char seq_line[10];
// 	printf("SO EINE VERDAMMTE SCHEISE ABER AUCH\n");
	char x;
	int number;
	int pos = 0;

	//first child
	while(fgets(line, LINE_LENGTH , edit_file)!=NULL)
	{

		x = line[0];
		if (x == '*')
			break;
		number = atoi(&line[1]);
		if (x == 'M')
		{
			while (number > 0)
			{
				if (aligned_sequence[pos] == 'X')
					--number;
				++pos;
			}
		}
		else if (x == 'I')
		{
			while (number > 0)
			{
				if (aligned_sequence[pos] == 'X')
					--number;
				++pos;
			}
		}
		else if (x == 'D')
		{
			while (number > 0)
			{
				if (aligned_sequence[pos] == 'X')
				{
					aligned_sequence[pos] = '-';
					--number;
				}
				++pos;
			}
		}
	}

	if (is_leaf1)
	{
// 		printf("LEAF\n");
		edit_seq2aligned_seq(aligned_sequence, sequence_file, seq_positions[child1], alignment_file);
	}
	else
	{
		fprintf(edit_seq_file, "%s", aligned_sequence);
		edit2alignment(sequence_file, seq_positions, edit_file, edit_positions, child1, number_of_sequences, aligned_sequence, alignment_length, edit_seq_file, offset, alignment_file);
	}

	//second child
	fseek(edit_seq_file, offset, SEEK_CUR);
	fgets(aligned_sequence, alignment_length+3, edit_seq_file);
	fseek(edit_seq_file, offset, SEEK_CUR);

	pos = 0;
	fseek(edit_file, edit_positions[node_number-number_of_sequences], SEEK_SET);
	while(fgets(line, LINE_LENGTH , edit_file)!=NULL)
	{
		x = line[0];
		if (x == '*')
			break;
		number = atoi(&line[1]);
		if (x == 'M')
		{
			while (number > 0)
			{
				if (aligned_sequence[pos] == 'X')
					--number;
				++pos;
			}
		}
		else if (x == 'I')
		{
			while (number > 0)
			{
				if (aligned_sequence[pos] == 'X')
				{
					aligned_sequence[pos] = '-';
					--number;
				}
				++pos;
			}
		}
		else if (x == 'D')
		{
			while (number > 0)
			{
				if (aligned_sequence[pos] == 'X')
					--number;
				++pos;
			}
		}
	}

	if (is_leaf2)
	{
		edit_seq2aligned_seq(aligned_sequence, sequence_file, seq_positions[child2], alignment_file);
	}
	else
	{
		fprintf(edit_seq_file, "%s", aligned_sequence);
		edit2alignment(sequence_file, seq_positions, edit_file, edit_positions, child2, number_of_sequences, aligned_sequence, alignment_length, edit_seq_file, offset, alignment_file);
	}
}




//  * The file has the follwing format (# and text behind are only comments and not included into the file):
// 	* 1		# Number of profile.
// 	* 1		# is leaf.
// 	* 5		# Number of columns in the profile.
// 	* 4A 1C	# In the first column are 4 'A' and 1 'C'
// 	* 3G	# In the second column are 3 'G'
// 	* 5A	# In the third column are 5 'A'
// 	* 2A 3C	# In the fourth column are 2 'A' and 3 'C'
// 	* 5C	# In the fifth column are 5 'C'
// 	* *		# Marks the end of this profile



/**
 * \brief Writes a profile to a file.
 *
 * \param sumup_prf The profile array, not a real profile.
 * \param length The length of the profile. The format can be seen in ./test.txt
 * \param file The FILE object to write the the profile into.
 * \param is_dna The type of sequence.
 * \param number The number of the profile.
 */
void
write2file(int **sumup_prf,
		   int length,
		   FILE *file,
		   int number,
		   int num_sequences,
		   Fastal_param *param_set)
{
	char *pos2aa = &(param_set->pos2char[0]);
	fprintf(file,"%i\n0\n%i\n1\n%i\n",number, length, num_sequences );
	int i, j;
	int alphabet_size = param_set->alphabet_size;

	i = 0;
	int x = 0;
	while (i < length)
	{
		for (j = 0; j < alphabet_size; ++j)
			if (sumup_prf[j][i] > 0)
			{
				if (x)
					fprintf(file," %c%i", pos2aa[j],sumup_prf[j][i]);
				else
					fprintf(file,"%c%i", pos2aa[j],sumup_prf[j][i]);
				x = 1;
			}
		x = 0;
		fprintf(file,"\n");
		++i;
	}
	fprintf(file,"*\n");
}



//*************************************   main function of the fasta algorithm ***********************************************

/**
*	\brief main of the fastal algorithm
*/
int fastal_main(int argc,		//number of arguments
			char **argv)	//arguments first = fastal, second = tree
{

	int i;
	//pointer to arguments
	void * method_arguments_p;
	int (*alignment_function)(Fastal_profile **profiles, Fastal_param *param_set, void *method_arguments_p, int is_dna, FILE *edit_file, FILE *prof_file, int number);

	struct Fastal_arguments arguments;

	arg_parse (argc, argv, &arguments);

	Fastal_param *param_set =(Fastal_param*) vcalloc(1,sizeof(Fastal_param));

	fill_parameters(arguments.is_dna, param_set, arguments.method, arguments.diag_method, arguments.mat);
	param_set->gep = arguments.gep;
	param_set->gop = arguments.gop;

// 	printf("%s",arguments.mat);
	if (arguments.evaluate)
	{
		printf("Calculate Sum of pairs Score\n");
		printf("Score: %f\n", calculate_sum_of_pairs_score_affine(arguments.sequence_file, param_set->M, param_set->gop, param_set->gep));
		vfree(param_set);
		exit(0);
	}


	if (arguments.agreement_score)
	{
		complete_agreement_score(arguments.aln2test, arguments.aln_ref);
		return 0;
	}


	if (arguments.num_ref_aln)
	{
		compute_ref_alignments(arguments.sequence_file, arguments.aln_ref, arguments.num_ref_aln, arguments.num_seq_in_ref);
		return 0;
	}



	int alphabet_size = param_set->alphabet_size;


	//sequence file management
// 	char **seq_name;
	long *file_positions = NULL;
	long **tmp = &file_positions;
	int number_of_sequences = make_index_of_file(arguments.sequence_file, tmp);



	//edit file management

// 	long current_edit_pos;
	long *edit_positions =(long int*) vcalloc(number_of_sequences,sizeof(long));


	//profile management
	Fastal_profile **profiles =(Fastal_profile**) vcalloc(3,sizeof(Fastal_profile*));
	initiate_profiles(profiles, param_set);
	FILE * prof_file = fopen(vtmpnam(NULL),"w+");
	long* profile_positions = (long int*)vcalloc(4,sizeof(long*));
	int max_prof = 4;
	int saved_prof = 0;


	printf("METHOD: %s\n",param_set->method);
	if (strcmp(param_set->method, "fast") == 0)
	{
		method_arguments_p = vcalloc(1,sizeof(Sparse_dynamic_param));
		fill_arguments_sparse((Sparse_dynamic_param*)method_arguments_p);
		alignment_function = sparse_dyn;
	}
	else if (strcmp(param_set->method, "nw") == 0)
	{
		method_arguments_p = vcalloc(1,sizeof(Nw_param));
		fill_arguments_nw((Nw_param*)method_arguments_p, alphabet_size);
		alignment_function = nw_dyn;
	}
	else if (strcmp(param_set->method, "gotoh") == 0)
	{
		method_arguments_p = vcalloc(1,sizeof(Gotoh_param));
		fill_arguments_gotoh((Gotoh_param*)method_arguments_p, alphabet_size);
		alignment_function = gotoh_dyn;
	}
	else if (strcmp(param_set->method, "udisc") == 0)
	{
		method_arguments_p = vcalloc(1,sizeof(Udisc_param));
		fill_arguments_gotoh((Gotoh_param*)method_arguments_p, alphabet_size);
		alignment_function = gotoh_dyn;

	}
	else
	{
		printf("ERROR - METHOD");
		exit(1);
	}


	if (arguments.gap_iterate)
	{
		iterate(param_set, method_arguments_p, arguments.sequence_file, arguments.output_file, arguments.gap_iterate);
		return 0;
	}

	if (arguments.tree_file == NULL)
	{
		arguments.tree_file = vtmpnam(NULL);
		printf("CONSTRUCT TREE\n");
		if (strcmp(arguments.tree_method, "parttree")==0)
		{
			make_partTree(arguments.sequence_file, arguments.tree_file, arguments.tree_param1, arguments.tree_param2, arguments.is_dna, 0);
		}
		else if (strcmp(arguments.tree_method, "oligotree") == 0)
		{
			compute_oligomer_distance_tree(arguments.sequence_file, param_set->char2pos, arguments.tree_file, arguments.tree_param2, arguments.tree_param1, param_set->alphabet_size);
		}

		if (arguments.tree_only == 1)
			return 0;
	}


	if (arguments.tree_out == 1)
	{
		char tree_out_file_name[500];
		sprintf(tree_out_file_name, "%s.tree",arguments.output_file);
		char const LINE_LENGTH = 50;
		char line[LINE_LENGTH];

		FILE* in = fopen(arguments.tree_file, "r");
		FILE* out = fopen(tree_out_file_name, "w");
		while( (fgets(line, LINE_LENGTH, in)) != NULL)
			fprintf(out, "%s", line);
		fclose(in);
		fclose(out);
	}





	FILE *seq_file = fopen(arguments.sequence_file,"r");
// 	FILE *edit_file = fopen(vtmpnam(NULL),"w+");
	FILE *edit_file = fopen("aha","w+");

	printf("CONSTRUCT ALIGNMENT\n");
	FILE *tree_file = fopen(arguments.tree_file,"r");
	const int LINE_LENGTH = 100;
	char line[LINE_LENGTH];
	char delims[] = " ";
	int node[3];
	int alignment_length = -1;
	node[2] = -1;


	//bottom-up traversal
	while(fgets(line, LINE_LENGTH, tree_file)!=NULL)
	{
		//read profiles
		node[0] = atoi(strtok(line,delims));
		node[1] = atoi(strtok(NULL,delims));
		node[2] = atoi(strtok(NULL,delims));

		//getting profile of second child
		if (node[1] < number_of_sequences)
		{
			file_pos2profile(seq_file, file_positions[node[1]], profiles[1], node[1], param_set);	//profile to save into
		}
		else
		{
			profile_file2profile(profiles[1], prof_file, profile_positions[--saved_prof], param_set);
			fseek (prof_file , profile_positions[saved_prof] , SEEK_SET);
		}

		//getting profile of first child
		if (node[0] < number_of_sequences)
		{
			file_pos2profile(seq_file, file_positions[node[0]], profiles[0], node[0], param_set);	//profile to save into
		}
		else
		{
			profile_file2profile(profiles[0], prof_file, profile_positions[--saved_prof], param_set);
			fseek (prof_file , profile_positions[saved_prof] , SEEK_SET);
		}


		if (saved_prof == max_prof)
		{
			max_prof += 5;
			profile_positions = (long int*)vrealloc(profile_positions, max_prof*sizeof(long));
		}

		edit_positions[node[2]-number_of_sequences] = ftell(edit_file);
		profile_positions[saved_prof] = ftell(prof_file);
		++saved_prof;

		//aligning the sequences
		alignment_length = alignment_function(profiles, param_set, method_arguments_p, arguments.is_dna, edit_file, prof_file, node[2]);
	}


//	bottom-down traversal (edit_files --> alignment)
//	tmp_out_file_name = vtmpnam(NULL);

// 	FILE *alignment_file = fopen(tmp_out_file_name, "w");
	FILE *alignment_file = fopen(arguments.output_file, "w");
	FILE *edit_seq_file = fopen(vtmpnam(NULL),"w+");

	char *aligned_sequence =(char*) vcalloc(alignment_length+3, sizeof(char));


	long offset = ftell(edit_seq_file);
	for (i = 0; i < alignment_length; ++i)
	{
		fprintf(edit_seq_file, "X");
		aligned_sequence[i] = 'X';
	}
	aligned_sequence[i]= '\n';
	aligned_sequence[i+1]= '\0';
	fprintf(edit_seq_file, "\n");
	offset = (ftell(edit_seq_file) - offset)*-1;


	edit2alignment(seq_file, file_positions, edit_file, edit_positions, node[2], number_of_sequences, aligned_sequence, alignment_length, edit_seq_file, offset, alignment_file);
	fclose(alignment_file);
	fclose(tree_file);
	fclose(edit_file);
	fclose(seq_file);
	fclose(edit_seq_file);

		//set stuff for the next cycle
// 		arguments.sequence_file	= tmp_out_file_name;


// // 		//DEBUG
// // 		char copy_command[500];
// // 		sprintf(copy_command, "cp %s %s_%i", tmp_out_file_name, arguments.output_file, cycle);
// // 		system(copy_command);
// 		++cycle;
// 	}

// 	printf("HERE_COPY\n");
// 	char copy_command[2000];
// 	sprintf(copy_command, "mv %s %s", tmp_out_file_name, arguments.output_file);
// 	printf("%s\n", copy_command);
// 	int error = system(copy_command);
// 	printf("ERROR %i\n", error);


	//free_memory & close files
	fclose(prof_file);
	free_fastal_profile(profiles[0], alphabet_size);
	free_fastal_profile(profiles[1], alphabet_size);
	vfree(profiles);
	vfree(profile_positions);



// number_of_sequences


	if (arguments.score)
	{
		printf("Calculate Score\n");
		double aln_score = calculate_sum_of_pairs_score_affine(arguments.output_file, param_set->M, param_set->gop, param_set->gep);
		printf("SCORE: %f\n", aln_score);
	}





	if (!strcmp(param_set->method, "fast"))
	{
		free_sparse((Sparse_dynamic_param*)method_arguments_p);
	}
	else if (!strcmp(param_set->method, "nw"))
	{
		free_nw((Nw_param*)method_arguments_p, alphabet_size);
	}
	else if (!strcmp(param_set->method, "gotoh"))
	{
		free_gotoh((Gotoh_param*)method_arguments_p, alphabet_size);
	}

	vfree(param_set);

	//free_memory & close files

	vfree(edit_positions);


	return 0;
}




//******************   toolbox   ***************************


/**
 * \brief Enlargement of the dynamic programming matrix in case it is to small.
 *
 * \param dyn_matrix The dynamic programming matrix.
 * \param old_length1 Current size of dimension 1.
 * \param old_length2 Current size of dimension 2.
 * \param length1 New size of dimension 1.
 * \param length2 New size of dimension 2.
 * \return Pointer to the new array.
 */
double**
resize_dyn_matrix(double **dyn_matrix,	//the dynamic programming matrix
				  int old_length1,		//old length of dimension 1
				  int old_length2,		//old length of dimension 2
				  int length1,			//new minimum length of dimension 1
				  int length2)			//new maximum length of dimension 2
{
	int i;
	if (old_length1 < length1)
	{
		dyn_matrix =(double**) vrealloc(dyn_matrix,length1*sizeof(double*));
		for (i = old_length1; i < length1; ++i)
			dyn_matrix[i] =(double*) vcalloc(old_length2,sizeof(double));
		old_length1 = length1;
	}

	if (old_length2 < length2)
	{
		for (i = 0;i<old_length1; ++i)
			dyn_matrix[i] =(double*) vrealloc(dyn_matrix[i], length2*sizeof(double));
		old_length2 = length2;
	}
	return dyn_matrix;
}



/**
 * \brief Frees the memory of a dynamic programming matrix.
 *
 * \param length1 Size of the first dimension of the matrix.
 * \param dyn_matrix The dynamic programming matrix.
 */
void
free_dyn_matrix(int length1,			//length of first dimension
				double **dyn_matrix)	//dynamic matrix
{
	int i = 0;
	for (; i<length1; ++i)
		vfree(dyn_matrix[i]);
	vfree(dyn_matrix);
}



/**
 * \brief Initialises the profiles with basic values.
 *
 * \param profiles Array of 3 profiles.
 * \param param_set The fastal parameters
 */
void
initiate_profiles(Fastal_profile **profiles,	//profiles pointer
				  Fastal_param *param_set)
{
	int alphabet_size = param_set->alphabet_size;
	int i,j;
	for (i =0; i < 3; ++i)
	{
		profiles[i] =(Fastal_profile*) vcalloc(1,sizeof(Fastal_profile));
		profiles[i]->weight = 1;
		profiles[i]->is_leaf = 1;
		profiles[i]->prf =(int**) vcalloc(alphabet_size, sizeof(int*));
		for (j = 0; j < alphabet_size; ++j)
		{
			profiles[i]->prf[j] =(int*) vcalloc(PROFILE_ENLARGEMENT, sizeof(int));
		}
		profiles[i]->allocated_memory = PROFILE_ENLARGEMENT;
	}
}




/**
 * \brief frees all memory occupied by the profile
 *
 * \param profile The profile to free.
 * \param alphabet_size The alphabet_size.
 */
void
free_fastal_profile(Fastal_profile* profile, int alphabet_size)
{
	--alphabet_size;
	for (;alphabet_size >= 0; --alphabet_size)
		vfree(profile->prf[alphabet_size]);
	vfree(profile->prf);
}


/**
 * \brief Initialize the Fastal parameter set.
 *
 * \param is_dna 1 when sequences are dna, 0 when not
 * \param param_set The fastal parameter set.
 * \param method The method to use in Fastal.
*/
void
fill_parameters(int is_dna, Fastal_param *param_set, char *method, char *diag_method, char *mat)
{
	sprintf(param_set->method,"%s",method);
	sprintf(param_set->diag_method,"%s",diag_method);
	int i;
	printf("%s",mat);
	param_set->M = read_matrice(mat);
	if (is_dna)
	{
		param_set->alphabet_size = 4;
		char tmp1[] = {'A','C','G','T'};

// 		int  tmp2[] = { 0, 1,  1,  0, -1, -1, 2, 0, -1, -1, 3, -1, 1, 0, -1, -1, -1, 0, 1, 3, 4, -1, 3, -1, 1, -1};
		for (i = 0; i<param_set->alphabet_size; ++i)
			param_set->pos2char[i] = tmp1[i];
// 		for (i = 0; i<26; ++i)
// 			param_set->char2pos[i] = tmp2[i];
		param_set->char2pos['A'] = 0;
		param_set->char2pos['B'] = 1;
		param_set->char2pos['C'] = 1;
		param_set->char2pos['D'] = 0;
		param_set->char2pos['G'] = 2;
		param_set->char2pos['H'] = 0;
		param_set->char2pos['K'] = 3;
		param_set->char2pos['M'] = 1;
		param_set->char2pos['N'] = 0;
		param_set->char2pos['R'] = 0;
		param_set->char2pos['S'] = 1;
		param_set->char2pos['T'] = 3;
		param_set->char2pos['U'] = 3;
		param_set->char2pos['W'] = 3;
		param_set->char2pos['Y'] = 1;
// 		param_set->M[0][3] = param_set->M[3][0] = -10;
// 		param_set->M[1][2] = param_set->M[2][1] = -10;
// 		param_set->M[0][1] = param_set->M[0][2] = param_set->M[1][0] = param_set->M[2][0] = -10;
// 		param_set->M[3][1] = param_set->M[3][2] = param_set->M[1][3] = param_set->M[2][3] = -10;
	}
	else
	{
		param_set->alphabet_size = 21;
		char tmp1[] = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','X'};
// 		int tmp2[] = { 0, 20,  1,  5, 16,  4,  2,  6,  7, 21,  8,  9,  10, 11, -1, 12, 13, 14, 15, 3, -1, 17,  18, 22, 19, 23};
		for (i = 0; i<param_set->alphabet_size; ++i)
			param_set->pos2char[i] = tmp1[i];
// 		for (i = 0; i<26; ++i)
// 			param_set->char2pos[i] = tmp2[i];
		param_set->char2pos['A'] = 0;
 		param_set->char2pos['B'] = 20;
		param_set->char2pos['C'] = 1;
		param_set->char2pos['D'] = 2;
		param_set->char2pos['E'] = 3;
		param_set->char2pos['F'] = 4;
		param_set->char2pos['G'] = 5;
		param_set->char2pos['H'] = 6;
		param_set->char2pos['I'] = 7;
		param_set->char2pos['J'] = 20;
		param_set->char2pos['K'] = 8;
		param_set->char2pos['L'] = 9;
		param_set->char2pos['M'] = 10;
		param_set->char2pos['N'] = 11;
		param_set->char2pos['P'] = 12;
		param_set->char2pos['Q'] = 13;
		param_set->char2pos['R'] = 14;
		param_set->char2pos['S'] = 15;
		param_set->char2pos['T'] = 16;
		param_set->char2pos['V'] = 17;
		param_set->char2pos['W'] = 18;
 		param_set->char2pos['X'] = 20;
		param_set->char2pos['Y'] = 19;
		param_set->char2pos['a'] = 0;
		param_set->char2pos['b'] = 20;
		param_set->char2pos['c'] = 1;
		param_set->char2pos['d'] = 2;
		param_set->char2pos['e'] = 3;
		param_set->char2pos['f'] = 4;
		param_set->char2pos['g'] = 5;
		param_set->char2pos['h'] = 6;
		param_set->char2pos['i'] = 7;
		param_set->char2pos['j'] = 20;
		param_set->char2pos['k'] = 8;
		param_set->char2pos['l'] = 9;
		param_set->char2pos['m'] = 10;
		param_set->char2pos['n'] = 11;
		param_set->char2pos['p'] = 12;
		param_set->char2pos['q'] = 13;
		param_set->char2pos['r'] = 14;
		param_set->char2pos['s'] = 15;
		param_set->char2pos['t'] = 16;
		param_set->char2pos['v'] = 17;
		param_set->char2pos['w'] = 18;
		param_set->char2pos['x'] = 20;
		param_set->char2pos['y'] = 19;
	}
}



