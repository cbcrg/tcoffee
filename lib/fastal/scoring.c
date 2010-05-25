
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <sys/stat.h>
#include <dirent.h>


#include "io_lib_header.h"
#include "util_lib_header.h"
#include "fastal_lib_header.h"




/**
 * \brief Function to calculate the sum of pairs score with quasi affine gap costs.
 *
 * \param alignment_file File with the aligned sequences in fasta_format.
 * \param file_index Positions of the sequences in \a alignment_file.
 * \param number_of_sequences The number of sequences in \a alignment_file.
 * \param score_matrix Matrix containing the scores.
 * \param gop The gap opening costs.
 * \param gep The gep extension costs.
 * \return The score of the alignment.
 */
double
calculate_sum_of_pairs_score_affine(char *alignment_file_name,
									int **score_matrix,
									double gop,
									double gep)
{

	const int LINE_LENGTH = 1000;
	int number_of_sequences = 0;
	FILE *alignment_file = fopen(alignment_file_name,"r");
	if (alignment_file == NULL)
	{
		printf("FILE COULD NOT BE OPENED!\n");
		exit(1);
	}
	int alignment_length = 0;
	char line[LINE_LENGTH+1];
	line[0] = '0';

	while (line[0] != '>')
	{
		fgets(line, LINE_LENGTH, alignment_file);
	}
	fgets(line, LINE_LENGTH, alignment_file);
	while (line[0] != '>')
	{
		alignment_length += strlen(line);
		fgets(line, LINE_LENGTH, alignment_file);
	}

	int alphabet_size = 28;
	int **counting = vcalloc(alphabet_size, sizeof(int*));

	int i;
	for (i = 0; i < alphabet_size; ++i)
	{
		counting[i] = vcalloc(alignment_length, sizeof(int));
	}

	alphabet_size -= 2;
	char c;
	int was_gap = 0;
	int pos_line, pos_alignment = -1;
	fseek (alignment_file, 0, SEEK_SET);
	pos_line = 0;

	int gap_line = 0;
	int current_size = 100;
	int gap_size = 0;
	int gap_pos = 0;

	int **gaps = vcalloc(100, sizeof(int*));

	while (fgets(line, LINE_LENGTH, alignment_file) != NULL)
	{
		if (line[0] != '>')
		{
			pos_line = -1;
			while (((c = line[++pos_line]) != '\n') && (c != '\0'))
			{
				if (isalpha(c))
				{
					++counting[toupper(c)-'A'][++pos_alignment];
					if (was_gap)
					{
						gaps[gap_line][gap_pos++] = pos_alignment-1;
					}
					was_gap = 0;
				}
				else
				{
					++counting[alphabet_size][++pos_alignment];
					if (!was_gap)
					{
						++counting[alphabet_size+1][pos_alignment];
						was_gap = 1;
						if (gap_pos >= gap_size-4)
						{
							gap_size += 50;
							gaps[gap_line] = vrealloc(gaps[gap_line], gap_size * sizeof(int));
							gaps[gap_line][gap_pos++] = pos_alignment;
						}
					}
				}
			}
		}
		else
		{
			if (number_of_sequences != 0)
			{
				if ((gap_pos % 2) != 0)
				{
					gaps[gap_line][gap_pos++] = pos_alignment-1;
				}
				gaps[gap_line][gap_pos++] = alignment_length+2;
				gaps[gap_line][gap_pos] = alignment_length+2;
				++gap_line;
			}
			if (current_size == gap_line)
			{
				current_size += 50;
				gaps = vrealloc(gaps, current_size * sizeof(int*));
			}
			gap_size = 2;
			gaps[gap_line] = vcalloc(gap_size, sizeof(int));
			++number_of_sequences;
			pos_alignment = -1;
			was_gap = 0;

			gap_pos = 0;
		}
	}

	gaps[gap_line][gap_pos++] = alignment_length+2;
	gaps[gap_line][gap_pos] = alignment_length+2;


	long double score = 0;
	int j, k, tmp2;
	int non_gap;
	for (i = 0; i < alignment_length; ++i)
	{
		for (j = 0; j < alphabet_size; ++j)
		{
			if ((tmp2 = counting[j][i]) >1)
			{
				while (tmp2 > 1)
				{
					score += score_matrix[j][j] * (--tmp2);
				}
			}
			for (k = j+1; k < alphabet_size; ++k)
			{
				score += score_matrix[j][k] * counting[j][i] * counting[k][i];
			}

		}


		non_gap = number_of_sequences - counting[alphabet_size][i];
		score += counting[alphabet_size][i] * non_gap  * gep;

	}

	++gap_line;
	alignment_length += 2;
	unsigned long gap_open = 0;

	int chunk = 500;

// 	#pragma omp parallel shared(gaps, chunk, alignment_length, gap_open) private(i, j)
	{
		int alignment_length2 = alignment_length;
// 		#pragma omp for schedule(dynamic,chunk) reduction(+:gap_open) nowait
		for (i = 0; i < gap_line; ++i)
		{
			int *gaps1 = gaps[i];
			for (j = i+1; j < gap_line; ++j)
			{
				int *gaps2 = gaps[j];
				int k = 0;
				int l = 0;
				while ((gaps1[k] != alignment_length2) && (gaps2[l] != alignment_length2))
				{
					if (gaps1[k+1] < gaps2[l])
					{
						++gap_open;
						k+=2;
						continue;
					}
					if (gaps2[l+1] < gaps1[k])
					{
						++gap_open;
						l+=2;
						continue;
					}
					if (gaps1[k] < gaps2[l])
					{
						if (gaps1[k+1] < gaps2[l+1])
						{
							++gap_open;
							k+=2;
							continue;
						}
						else
						{
							l+=2;
							continue;
						}
					}
					if (gaps1[k] > gaps2[l])
					{
						if (gaps1[k+1] > gaps2[l+1])
						{
							++gap_open;
							l+=2;
							continue;
						}
						else
						{
							k+=2;
							continue;
						}
					}
					if (gaps1[k] == gaps2[l])
					{
						if (gaps1[k+1] == gaps2[l+1])
						{
							k+=2;
							l+=2;
							continue;
						}
						else
						{
							if (gaps1[k+1] <gaps2[l+1])
							{
								k+=2;
							}
							else
							{
								l+=2;
							}
						}
					}

				}
				while (gaps1[k] != alignment_length2)
				{
					++gap_open;
					k+=2;
				}
				while (gaps2[l] != alignment_length2)
				{
					++gap_open;
					l+=2;
				}
			}
		}
	}

	score += gop * gap_open;


	return score / number_of_sequences;
}


int
compare (const void * a, const void * b)
{
	return (*(int*)a - *(int*)b);
}



int
get_random_value(int a, int b)
{
	return rand() % (b - a + 1) + a;
// 	j = 1 + (int)( 10.0 * rand() / ( RAND_MAX + 1.0 ) );
}





void
compute_ref_alignments(char *seq_file_name, char* ref_directory, int num_alignments, int num_seq_in_aln)
{

	printf("%s %s %i %i\n", seq_file_name, ref_directory, num_alignments, num_seq_in_aln);
	//check if target directory exists, if not make it
	char directory[500];
	if (ref_directory[0]== '~')
	{
		sprintf(directory,"%s%s",getenv("HOME"),ref_directory+1);
	}
	else
	{
		sprintf(directory,"%s",ref_directory);
	}
	struct stat st;
	if(stat(directory, &st) != 0)
	{
		mkdir(directory,0700);
	}

	printf("made\n");

	// get sequences & produce alignments
	//make index
	long *file_positions = NULL;
	long **tmp_pos = &file_positions;
	int number_of_sequences = make_index_of_file(seq_file_name, tmp_pos);
	char *tmp_file_name = vtmpnam(NULL);
	FILE *seq_file = fopen(seq_file_name, "r");
	int *already_taken = vcalloc(num_seq_in_aln, sizeof(int));
	int i, j, k, tmp, already, num_already;
	for (i = 0; i < num_alignments; ++i)
	{
		//choose randomly the sequences

		num_already = 0;
		while (num_already < num_seq_in_aln)
		{
			already = 0;
			tmp = get_random_value(0, number_of_sequences-1);
			for (k = 0; k < num_already; ++k)
			{
				if (already_taken[k] == tmp)
				{
					already = 1;
					break;
				}
			}
			if (!already)
			{
				already_taken[num_already] = tmp;
				++num_already;
			}
		}
		//write sequences into temporary file
		const int LINE_LENGTH = 200;
		char line[LINE_LENGTH];


		FILE* tmp_file = fopen(tmp_file_name, "w");
		for (k = 0; k < num_already; ++k)
		{
			fseek(seq_file, file_positions[already_taken[k]], SEEK_SET);
			fgets(line, LINE_LENGTH, seq_file);
			fprintf(tmp_file, "%s", line);
			while( (fgets(line, LINE_LENGTH, seq_file)) != NULL)
			{
				if (line[0] == '>')
					break;
				fprintf(tmp_file, "%s", line);
			}
		}
		fclose(tmp_file);

		//calculate alignment
		char aln_command[1000];
		sprintf(aln_command,"t_coffee -in %s -outfile %s/%i.fa -output fasta_aln 2>/dev/null >/dev/null", tmp_file_name, directory, i);
		system(aln_command);
	}

	vfree(already_taken);
}



// void
// compute agreement_score(char *alignment_file_name, char *alignment_directory)
// {
//
// }




/**
 * This algorithm choses according to a tree a number of sequenes of which it computes a reference algnment.
 * \param seq_file_name The sequence file.
 * \param tree_file_name The given tree.
 * \param ref_aln_name The file where the reference alignment should be written to.
 * \param num_seq_in_ref The number of sequences to choose.
**/
void
make_ref_alignment(char *seq_file_name, char *tree_file_name, char *ref_aln_name, int num_seq_in_ref)
{
	const int LINE_LENGTH = 200;
	char line[LINE_LENGTH];

	FILE *tree_f = fopen(tree_file_name,"r");


	long *file_positions = NULL;
	long **tmp = &file_positions;
	int num_sequences = make_index_of_file(seq_file_name, tmp);


	int every_x = num_sequences/(num_seq_in_ref);

	int *seq_ids = vcalloc(num_seq_in_ref,sizeof(int));
	char delims[] = " ";

	int pos = -1;
	fseek (tree_f, 0, SEEK_SET);
	int node;
	int dist = every_x;
	while(fgets(line, LINE_LENGTH, tree_f)!=NULL)
	{
		node = atoi(strtok(line,delims));
		if (node < num_sequences)
		{
			if (dist == every_x)
			{

				seq_ids[++pos] = node;
				dist = 0;
			}
			else
			{
				++dist;
			}
		}

		node = atoi(strtok(NULL,delims));
		if (node < num_sequences)
		{
			if (dist == every_x)
			{
				seq_ids[++pos] = node;
				dist = 0;
			}
			else
			{
				++dist;
			}
		}
	}
	fclose(tree_f);


	int i = 0;
	qsort(seq_ids, num_seq_in_ref, sizeof(int), compare);

	char *seq_ref_name = vtmpnam(NULL);
	FILE *seq_ref_f = fopen(seq_ref_name, "w");;
	FILE *seq_f = fopen(seq_file_name, "r");

	for (i = 0; i < num_seq_in_ref; ++i)
	{
		fseek (seq_f, file_positions[seq_ids[i]], SEEK_SET);
		fgets(line, LINE_LENGTH, seq_f);
		fprintf(seq_ref_f, "%s", line);
		while(fgets(line, LINE_LENGTH, seq_f)!=NULL)
		{
			if (line[0] == '>')
				break;
			fprintf(seq_ref_f, "%s", line);
		}
	}
	fclose(seq_ref_f);
	fclose(seq_f);

	char command[500];
	sprintf(command, "mafft --quiet %s > %s", seq_ref_name, ref_aln_name);
	system(command);
}


/**
 * This function calculates the agreement between two alignments in a specified number of sequences.
 * \param ref_file_name The reference alignment.
 * \param aln_file_name The test alignment.
 **/
double
agreement_score(char *ref_file_name, char *aln_file_name)
{
	const int LINE_LENGTH = 200;
	char line[LINE_LENGTH];

	int number_of_sequences = 0;
	FILE *ref_f = fopen(ref_file_name,"r");
	while(fgets(line, LINE_LENGTH, ref_f)!=NULL)
	{
		if (line[0] == '>')
		{
			++number_of_sequences;
		}
	}
	fseek(ref_f, 0, SEEK_SET);
	char **seq_names = vcalloc(number_of_sequences, sizeof(char*));
	int i = 0;
	while(fgets(line, LINE_LENGTH, ref_f)!=NULL)
	{
		if (line[0] == '>')
		{
			seq_names[i] = vcalloc(LINE_LENGTH, sizeof(char));
			sprintf(seq_names[i], "%s",line);
			++i;
		}
	}


	fclose(ref_f);
	FILE *aln_f = fopen(aln_file_name,"r");
	char *tmp_name = vtmpnam(NULL);
	FILE *tmp_f = fopen(tmp_name,"w");
	int x = 0;
	long last_pos = -1;
	while(fgets(line, LINE_LENGTH, aln_f)!=NULL)
	{
		if (line[0] == '>')
		{
			for (i = 0; i < number_of_sequences; ++i)
			{
				if (!strcmp(line, seq_names[i]))
				{
					fprintf(tmp_f, "%s", line);
					last_pos = ftell(aln_f);
					while(fgets(line, LINE_LENGTH, aln_f)!=NULL)
					{
						if (line[0] == '>')
							break;
						fprintf(tmp_f, "%s", line);
// 						printf("ARSCH!\n");
						last_pos = ftell(aln_f);
					}
					++x;
					fseek(aln_f, last_pos, SEEK_SET);
					break;
				}
			}
		}
		if (x == number_of_sequences)
			break;
	}
	fclose(aln_f);
	fclose(tmp_f);
	char command[500];
// 	sprintf(command, "cp %s /users/cn/ckemena/Desktop/1.fa", tmp_name);
// 	system(command);
// 	fp = popen("ls -l", "r");
	sprintf(command, "t_coffee -other_pg aln_compare -al1 %s -al2 %s ", ref_file_name, tmp_name);


	FILE *result;
	result = popen(command, "r");
	fgets(line, LINE_LENGTH, result);
	fgets(line, LINE_LENGTH, result);
	fgets(line, LINE_LENGTH, result);

// 	sprintf(command, "cp %s ~/Destkop/1.fa", tmp_name);
// 	system(command);

	char delims[] = " ";
	char *tmp_str;
	strtok(line, delims);
	for (i = 0; i < 3; ++i)
	{
		while((tmp_str = strtok(NULL, delims))==NULL);
	}
	return(atof(tmp_str));
}


complete_agreement_score(char *aln_file_name, const char *ref_directory)
{
	struct dirent *dp;
	char *name = strrchr(aln_file_name,'/')+1;
// 	printf("%s ",name);
     // enter existing path to directory below
	char directory[500];
	if (ref_directory[0]== '~')
	{
		sprintf(directory,"%s%s",getenv("HOME"),ref_directory+1);
	}
	else
	{
		sprintf(directory,"%s",ref_directory);
	}
	DIR *dir = opendir(directory);

	char ref_file_name[200];
	while ((dp=readdir(dir)) != NULL)
	{
		if ((strcmp(dp->d_name,".")) && (strcmp(dp->d_name,"..")))
		{

			sprintf(ref_file_name, "%s/%s",directory, dp->d_name);
			printf("%s %f\n",name, agreement_score(ref_file_name, aln_file_name));
// 			printf("%f ", agreement_score(ref_file_name, aln_file_name));

		}
	}
	printf("\n");
	closedir(dir);
	return 0;

}


// 	char delims[] = " ";
// 	int node[3];
// 	int alignment_length = -1;
// 	node[2] = -1;


	//bottom-up traversal
// 	while(fgets(line, LINE_LENGTH, tree_file)!=NULL)
// {
// 		//read profiles
// 	node[0] = atoi(strtok(line,delims));
// 	node[1] = atoi(strtok(NULL,delims));
// 	node[2] = atoi(strtok(NULL,delims));
