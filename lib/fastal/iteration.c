
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include <ctype.h>
// #include <string.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "fastal_lib_header.h"


// the following fuctions are defined in fastal.c
int gotoh_dyn(Fastal_profile **profiles, Fastal_param *param_set, void *method_arguments_p, int is_dna, FILE *edit_file, FILE *prof_file, int number);
void free_nw(Nw_param* method_arguments_p, int alphabet_size);

/**
 * Adds a single sequence to a profile.
 * \param seq The sequence.
 * \param prf The profile.
 */
void
seq2profile2(char *seq, Fastal_profile *prf, int *char2pos)
{
// 	printf("SEDQ: %s\n", seq);
	int **prf_in = prf->prf;
	int i = -1;
	--seq;
	while (*(++seq) !='\0')
	{	
		if (*seq != '-')
			++(prf_in[char2pos[*seq]][++i]);
		else
			++i;
// 		++seq;
	}
}




// typedef struct
// {
// 	/// Number of sequences in this profile
// 	int num_sequences;
// 	/// number of the profile
// 	int prf_number;
// 	///0 = combination of two profiles, 1 = profile of a single sequence -> name1 = seq_name
// 	int is_leaf;
// 	///length of the profile
// 	int length;
// 	///weight of the sequence
// 	int weight;
// 	///saves the amount of allocated memory
// 	int allocated_memory;	
// 	///the profile itself [alphabet_size][profile_length]
// 	int **prf;
// 	///number_of_sequences
// 	int number_of_sequences;
// }
// Fastal_profile;





/**
 * Deletes all gap columns from the profile.
 * 
 * \param prf The profile.
 * \param alphabet_size The size of the alphabet.
 * \param gap_list The gap_list.
 * \param gap_list_length The length of the gap list.
 * \param num_gaps The number of gaps.
 * \return The new gap list.
 */
int *
del_gap_from_profile(Fastal_profile *prf, int alphabet_size, int *gap_list, int *gap_list_length, int *num_gaps)
{
	int i;
	int pos = -1;
	int not_gap = 0;
	int gap_pos = 0;
	int write_pos = 0;
// 	int gap_counter = 0;
	int **prf_in = prf->prf;
	int prf_length = prf->length;
	
	while (gap_pos < prf_length)
	{
// 		printf("PRF1: %i %i %i %i\n", prf_in[0][gap_pos], prf_in[1][gap_pos], prf_in[2][gap_pos], prf_in[3][gap_pos]);

		not_gap = 0;
		for (i = 0; i < alphabet_size; ++i)
		{
			if (prf_in[i][gap_pos])
			{
				not_gap = 1;
// 				printf("ARG: %i\n", i);
				break;
			}
		}
		if (not_gap)
		{
			
			if (gap_pos != write_pos)
			{
				for (i = 0; i < alphabet_size; ++i)
				{
					prf_in[i][write_pos] = prf_in[i][gap_pos];
				}
			}
			++write_pos;
		}
		else
		{
			
			if (++pos >= *gap_list_length)
			{
				*gap_list_length += 10;
				gap_list = (int*)vrealloc(gap_list, (*gap_list_length)*sizeof(int));
			}
			gap_list[pos] = gap_pos;
// 			++gap_counter;
		}
		++gap_pos;
	}
// 	printf("DEEEEEEEEEEEEEEEEEL %i %i\n", gap_counter, pos+1);
	
	*num_gaps = pos+1;
	prf->length -= *num_gaps;
	return gap_list;
}


/**
 * Splits a given alignment according to a given position.
 * \param alignment_file The alignment file.
 * \param gap_prf Saves the profile with the gap sequences.
 * \param no_gap_prf Saves the profile of the non gap sequences.
 * \param seq A temporary place to save the sequence.
 * \param index The index where to look for the gaps.
 * \param char2pos The character to position translation array.
 */
void
split_set(FILE *alignment_file, Fastal_profile *gap_prf, Fastal_profile *no_gap_prf, char *seq, int index, int *char2pos, char* split_file_name)
{
	FILE *split_set_file1 = fopen("GAP", "w");
	FILE *split_set_file2 = fopen("NOGAP", "w");
// 	FILE *shorted = fopen("SHORTT.fa", "w");
	printf("INDEX: %i\n",index);
	fseek (alignment_file , 0 , SEEK_SET);
	FILE *split_file = fopen(split_file_name, "w");
	int i,j;
	const int LINE_LENGTH = 200;
	char line[LINE_LENGTH];
	int pos = -1;
	int overall_pos = -1;
	fgets(line, LINE_LENGTH , alignment_file);
// 	fprintf(shorted, "%s", line);
	while(fgets(line, LINE_LENGTH , alignment_file)!=NULL)
	{
		pos = -1;
		if (line[0] == '>')
		{
// 			fprintf(shorted, "%s", line);
			if (seq[index] == '-')
			{
				
				seq[++overall_pos] = '\0';
				fprintf(split_file, "g\n");
				fprintf(split_set_file1,"%s\n",seq);
				seq2profile2(seq, gap_prf, char2pos);
				++(gap_prf->number_of_sequences);
			}
			else
			{
				seq[++overall_pos] = '\0';
				fprintf(split_file, "n\n");
				fprintf(split_set_file2,"%s\n",seq);
				seq2profile2(seq, no_gap_prf, char2pos);
				++(no_gap_prf->number_of_sequences);
			}
			overall_pos = -1;
			
		}
		else
		{
// 			for (j = 35; j < 96; ++j)
// 			{
// 				fprintf(shorted, "%c", line[j]);
// 			}
// 			fprintf(shorted, "\n");
			while ((line[++pos] != '\n') && (line[pos] != '\0'))
			{
				seq[++overall_pos] = line[pos];
			}
		}
	}

	if (seq[index] == '-')
	{
		seq[++overall_pos] = '\0';
		fprintf(split_file, "g\n");
		seq2profile2(seq, gap_prf, char2pos);
		++(gap_prf->number_of_sequences);
	}
	else
	{
		seq[++overall_pos] = '\0';
		fprintf(split_file, "n\n");
		seq2profile2(seq, no_gap_prf, char2pos);
		++(no_gap_prf->number_of_sequences);
	}
	fclose(split_file);
}


Fastal_profile *
enlarge_prof(Fastal_profile *prof, int new_length, int alphabet_size)
{
	if (new_length > prof->allocated_memory)
	{
		int i;
		for (i = 0;i < alphabet_size; ++i)
		{
			prof->prf[i] = (int*)vrealloc(prof->prf[i],new_length*sizeof(int));
		}
		prof->allocated_memory = new_length;
	}
	return prof;
}




void
iterate(Fastal_param *param, void *method_arguments_p, char *aln_file_name, char *out_file_name_end, int iteration_number)
{
// 	calculate alignment length
	
//count gap
	int it_coutner_2 = 0;
	const int LINE_LENGTH = 200;
	char line[LINE_LENGTH];
	char *seq1 = (char*)vcalloc(1,sizeof(char));
	char *seq2 = (char*)vcalloc(1,sizeof(char));
	Fastal_profile **profiles =(Fastal_profile**) vcalloc(3,sizeof(Fastal_profile*));
	initiate_profiles(profiles, param);
	Fastal_profile *gap_prf = profiles[0];
	Fastal_profile *no_gap_prf = profiles[1];
	int alphabet_size = param->alphabet_size;

	int *gap_list_1 = (int*)vcalloc(1, sizeof(int));
	int *gap_list_1_length = (int*)vcalloc(1, sizeof(int));
	*gap_list_1_length = 1;
	int num_gaps_1 = 0;
	int *gap_list_2 = (int*)vcalloc(1, sizeof(int));
	int *gap_list_2_length = (int*)vcalloc(1, sizeof(int));
	*gap_list_2_length = 1;
	int num_gaps_2 = 0;

	int alignment_length = 1;
	//from here repeat!
	int it_counter = 0;
	char *out_file_name = aln_file_name;
	int *gap_profile = (int*)vcalloc(alignment_length, sizeof(int));
	
// 	while (it_counter < alignment_length)
// 	{
		
		
		FILE *alignment_file = fopen(aln_file_name,"r");
		if (alignment_file == NULL)
		{
			printf("Could not open alignment file %s\n", aln_file_name);
			exit(1);
		}
	
	
		
		alignment_length = 0;
		int tmp_len = -1;
		fgets(line, LINE_LENGTH , alignment_file);
		while(fgets(line, LINE_LENGTH , alignment_file)!=NULL)
		{
			tmp_len = -1;
			if (line[0] == '>')
			{
				break;
			}
			while ((line[++tmp_len] != '\n') && (line[tmp_len] != '\0'));
			alignment_length += tmp_len;
		}
	// 	printf("ALN_LENGTH %i\n", alignment_length);
		seq1 =(char*)vrealloc(seq1, (1+alignment_length)*sizeof(char));
		
		gap_profile = (int*)vrealloc(gap_profile, alignment_length * sizeof(int));
		int i;
		for (i = 0; i < alignment_length; ++i)
		{
			gap_profile[i] = 0;
		}
	
		int number_of_sequences = 0;
		int pos = -1;
		fseek (alignment_file , 0 , SEEK_SET);
		while(fgets(line, LINE_LENGTH , alignment_file)!=NULL)
		{
			if (line[0] == '>')
			{
				++number_of_sequences;
				pos = -1;
			}
			else
			{
				i = -1;
				while ((line[++i] != '\n') && (line[i] != '\0'))
				{
					++pos;
					if (line[i] == '-')
					{
						++gap_profile[pos];
					}
				}
			}
		}
	
		double max_entrop = 0;
		int column_index = 0;
		double entrop = 0;
		double last = 0;
		double p;

		
// 		for (i = it_counter; i<=it_counter; ++i)
// 		{
// 			p = gap_profile[i]/(double)number_of_sequences;
// 			if (!p)
// 			{
// 				entrop = 0;
// 			}
// 			else
// 			{
// 				entrop = (-1)*(p*log(p) + (1-p)*log(1-p) ) ;
// 			}
// 			if (entrop > max_entrop)
// 			{
// 				column_index = i;
// 				max_entrop = entrop;
// 			}
// 			last = entrop;
// 		}
// 		++it_counter;
// 		if (max_entrop < 0.6)
// 		{
// 			printf("CONTINUE %f\n",entrop);
// 			continue;
// 		}
		
// 		column_index = 18;//it_counter-1;
// 		if (column_index == 19)
			column_index = 58;
		out_file_name = vtmpnam(NULL);
		
		char *edit_file_name = "EDIT";//vtmpnam(NULL);
		FILE *edit_file = fopen(edit_file_name,"w");
		char *profile_file_name = vtmpnam(NULL);
		FILE *profile_file = fopen(profile_file_name,"w");
		char *split_file_name = "AHA";//vtmpnam(NULL);
		++it_coutner_2;

		gap_prf = enlarge_prof(gap_prf, alignment_length, alphabet_size);
		no_gap_prf = enlarge_prof(no_gap_prf, alignment_length, alphabet_size );
		no_gap_prf->number_of_sequences = 0;
		gap_prf->number_of_sequences = 0;
	
		split_set(alignment_file, gap_prf, no_gap_prf, seq1, column_index, param->char2pos, split_file_name);
		gap_prf -> length = alignment_length;
		no_gap_prf -> length = alignment_length;


		gap_list_1 = del_gap_from_profile(gap_prf, alphabet_size, gap_list_1, gap_list_1_length, &num_gaps_1 );
		gap_list_2 = del_gap_from_profile(no_gap_prf, alphabet_size, gap_list_2, gap_list_2_length, &num_gaps_2 );


		fclose(alignment_file);
		profiles[0] = gap_prf;
		profiles[1] = no_gap_prf;

		alignment_length = gotoh_dyn(profiles, param, method_arguments_p, 0, edit_file, profile_file, 0);
		seq1 =(char*)vrealloc(seq1, (1+alignment_length)*sizeof(char));
		seq2 =(char*)vrealloc(seq2, (1+alignment_length)*sizeof(char));
		
		fclose(edit_file);
		edit_file = fopen(edit_file_name,"r");
		edit2seq_pattern(edit_file, seq1, seq2);
		fclose(edit_file);
		fclose(profile_file);
		
		write_iterated_aln(aln_file_name, out_file_name, split_file_name, seq1, gap_list_1, num_gaps_1, seq2, gap_list_2, num_gaps_2);
		aln_file_name = out_file_name;

// 		if (it_coutner_2 == 2)
// 		break;
// 	}
	char command[1000];
	sprintf(command, "mv %s %s",out_file_name, out_file_name_end);
	system(command);
	
	vfree(seq1);
	vfree(seq2);
	vfree(gap_list_1);
	vfree(gap_list_2);
	if (!strcmp(param->method, "fast"))
	{
		free_sparse((Sparse_dynamic_param*)method_arguments_p);
	}
	else if (!strcmp(param->method, "nw"))
	{
		free_nw((Nw_param*)method_arguments_p, alphabet_size);
	}
	else if (!strcmp(param->method, "gotoh"))
	{
		free_gotoh((Gotoh_param*)method_arguments_p, alphabet_size);
	}
}




void
write_iterated_aln(char* old_aln_file_name, char* new_aln_file_name, char *gap_file_name, char *seq1, int *gap_list1, int num_gap1, char *seq2, int *gap_list2, int num_gap2)
{
	int i;
// 	for (i = 0; i < num_gap1; ++i)
// 	{
// 		printf("%i ", gap_list1[i]);
// 	}
// 	printf("\n-----\n");
// 	for (i = 0; i < num_gap2; ++i)
// 	{
// 		printf("%i ", gap_list2[i]);
// 	}
// 	printf("\n");
// 	
// 	printf("%s\n%s\n",seq1, seq2);
// 	printf("EX: %i\n", gap_list1[0]);
	FILE *gap_file = fopen(gap_file_name, "r");
	FILE *new_aln_file = fopen(new_aln_file_name, "w");
// 	if (new_aln_file == NULL)
// 		printf("AHA\n");
	FILE *old_aln_file = fopen(old_aln_file_name, "r");
	
	const int GAPLINE_LENGTH = 5;
	char gapline[GAPLINE_LENGTH];
	char *seq;
	int *gap_list;
	int num_gap;
	
	const int ALNLINE_LENGTH = 200;
	char alnline[ALNLINE_LENGTH];
	
	int alnline_pos  = -1;
	int overall_aln_pos = -1;
	int pattern_pos = -1;
	int gap_pos = 0;
	
	
// 	fprintf(new_aln_file, "%s", alnline);
	while(fgets(gapline, GAPLINE_LENGTH , gap_file)!=NULL)
	{
		
		//this is the sequence name
		fgets(alnline, ALNLINE_LENGTH , old_aln_file);
		fprintf(new_aln_file,"%s",alnline);
		//getting (part if) the sequence
		fgets(alnline, ALNLINE_LENGTH , old_aln_file);
		if (gapline[0] == 'g')
		{
			seq = seq1;
			gap_list = gap_list1;
			num_gap = num_gap1;
		}
		else
		{
			seq = seq2;
			gap_list = gap_list2;
			num_gap = num_gap2;
		}
		overall_aln_pos = 0;
		pattern_pos = 0;
		gap_pos = 0;
		alnline_pos = -1;
		
		//go along the pattern
		while (seq[pattern_pos] != '\0')
		{
			if (seq[pattern_pos] == '-')
			{
				++pattern_pos;
				fprintf(new_aln_file,"-");
				continue;
			}
			//get next part of the sequence if you are at the end
			if ((alnline[++alnline_pos] == '\n') || (alnline[alnline_pos] == '\0'))
			{
				alnline_pos = 0;
				fgets(alnline, ALNLINE_LENGTH, old_aln_file);
			}
			if ((gap_pos < num_gap) && (overall_aln_pos == gap_list[gap_pos]))
			{
				++gap_pos;
				++overall_aln_pos;
				continue;
			}
			fprintf(new_aln_file,"%c",alnline[alnline_pos]);
			++overall_aln_pos;
			++pattern_pos;
		}
		fprintf(new_aln_file,"\n");
	}
	
	fclose(new_aln_file);
	fclose(old_aln_file);
	fclose(gap_file);
}


void
edit2seq_pattern(FILE *edit_file, char *seq1, char *seq2)
{
	
	fseek(edit_file, 0, SEEK_SET);
	const int LINE_LENGTH = 50;
	char line[LINE_LENGTH];
	fgets(line, LINE_LENGTH , edit_file);
// 	int child1 = atoi(line);
	fgets(line, LINE_LENGTH , edit_file);
// 	int child2 = atoi(line);
	fgets(line, LINE_LENGTH , edit_file);
// 	int is_leaf1 = atoi(line);
	fgets(line, LINE_LENGTH , edit_file);
// 	int is_leaf2 = atoi(line);

	char x;
	int number;
	int pos = -1;
	while(fgets(line, LINE_LENGTH , edit_file)!=NULL)
	{
// 		printf("%s\n",line);
		x = line[0];
		if (x == '*')
			break;
		number = atoi(&line[1]);
		if (x == 'M')
		{
			while (--number >= 0)
			{
				seq1[++pos] = 'X';
				seq2[pos] = 'X';
			}
		}
		else if (x == 'I')
		{
			while (--number >= 0)
			{
				seq1[++pos] = 'X';
				seq2[pos] = '-';
			}
		}
		else if (x == 'D')
		{
			while (--number >= 0)
			{
// 				printf("POS: %i\n", pos);
				seq1[++pos] = '-';
				seq2[pos] = 'X';
			}
		}
	}
	seq1[++pos] = '\0';
	seq2[pos] = '\0';
// 	printf("%s\n%s\n",seq1, seq2);
}
