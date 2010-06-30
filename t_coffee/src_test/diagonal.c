#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "fastal_lib_header.h"



/*!
 *	\file diagonal.c
 *	\brief Source code for the calculation and preparation of diagonals
 */







int 
diagonal_compare (const void * a, const void * b)
{
	
	int tmp = ((((Diagonal*)a)->y - ((Diagonal*)a)->x) - (((Diagonal*)b)->y - ((Diagonal*)b)->x));
	if (tmp != 0)
		return tmp;
	return (((Diagonal*)a)->y - ((Diagonal*)b)->y);
}

int
max(int a, int b)
{
	if (a < b)
		return b;
	else
		return a;
}



Segment*
extend_diagonals(Diagonal *diagonals, int *num_diagonal, int l1, int l2)
{
// 	sort diagonal by diagonal number
	int i;
// 	printf("INPUT_PRE\n");
// 	for (i = 0; i < num_diagonals; ++i)
// 	{
// 		printf("%i %i %i\n",diagonals[i].x, diagonals[i].y, diagonals[i].length);
// 	}
	int num_diagonals = *num_diagonal;
	qsort (diagonals, num_diagonals, sizeof(Diagonal), diagonal_compare);


//	find nearest diagonal and expand
//	make shure that overlapping segments on the same diagonal are merged

	int diff;
// 	printf("INPUT\n");
	for (i = 0; i < num_diagonals; ++i)
	{
		printf("%i %i %i\n",diagonals[i].x, diagonals[i].y, diagonals[i].length);
	}
	for (i = 0; i < num_diagonals; ++i)
	{
		diagonals[i].end_exp = 0;
		diagonals[i].front_exp = 0;
	}
	int first = 0;
	int next = 1;
	while (next < num_diagonals)
	{

		if (!((diagonals[next].y >= diagonals[first].y) && diagonals[next].x <= diagonals[first].x))
		{
			printf("%i %i %i %i %i %i\n",diagonals[first].x, diagonals[first].y, diagonals[next].x, diagonals[next].y, diagonals[first].y - diagonals[first].x, diagonals[next].y - diagonals[next].x);
			diff = diagonals[next].y - (diagonals[first].y + diagonals[first].length);
			if (diff > 0)
			{
				if ((diagonals[first].y - diagonals[first].x) != ((diagonals[next].y - diagonals[next].x)))
				{
					diagonals[first].end_exp = max(diagonals[first].end_exp, diff);
					diagonals[next].front_exp = max(diagonals[next].front_exp, diff);
	// 				if (diagonals[next].x - diagonals[next].front_exp < 0)
	
					while (diagonals[++first].x == -1);
				}
				else
				{
					printf("ICH TUWAS\n");
					diagonals[first].length = diagonals[next].y + diagonals[next].length - diagonals[first].y;
					diagonals[first].end_exp = diagonals[next].end_exp;
					diagonals[next].x = -1;
				}
			}
			diff = diagonals[first].x - (diagonals[next].x + diagonals[next].length);
			if (diff > 0)
			{
	// 			diff = diagonals[next].x + diagonals[next].length - diagonals[first].x;
	
				diagonals[first].front_exp =  max(diagonals[first].front_exp, diff);
				diagonals[next].end_exp = max(diagonals[next].end_exp, diff);
	// 			++num_segments;
	// 			if (diagonals[first].x - diagonals[first].front_exp < 0)
	
				while (diagonals[++first].x == -1);
			}
		} else
			++ first;
		++next;
		
	}

	int num_segments =0;
	Segment *seg = vcalloc(num_diagonals, sizeof(Segment));
	int pos = 0;
// 	int diag_num1, diag_num2;
	int *tmp1;
	Diagonal *tmp2;
	for (i = 0; i < num_diagonals; ++i)
	{
		if (diagonals[i].x != -1)
		{
			++num_segments;
			tmp2 = &diagonals[i];
			seg[pos].segments = vcalloc(6, sizeof(int));
			seg[pos].current_pos = seg[pos].segments;
			tmp1 = seg[pos].segments;
			tmp1[0] = tmp2->x - tmp2->front_exp;
			tmp1[1] = tmp2->y - tmp2->front_exp;
			tmp1[2] = tmp2->front_exp + tmp2->length + tmp2->end_exp;
			tmp1[5] = 1;
			tmp1[4] = l2;
			tmp1[3] = tmp1[0] + (l2 -  tmp1[1]);
			++pos;
		}
	}

	for ( i = 0; i < num_segments; ++i )
	{
		printf("%i %i %i || %i %i %i\n", seg[i].current_pos[0], seg[i].current_pos[1], seg[i].current_pos[2], seg[i].current_pos[3], seg[i].current_pos[4], seg[i].current_pos[5]);
	}
// 	exit(1);
	vrealloc(seg,num_segments*sizeof(Segment));
	
	*num_diagonal = num_segments;
	
// 	vfree(diagonals);
	return seg;
}


int **
segments2int(Segment *diagonals,
			int num_diagonals,
			char *seq1,
			char *seq2,
			Fastal_profile *profile1,
			Fastal_profile *profile2,
			int *num_points,
			Fastal_param *param_set)
{
// 	printf("NUM: %i\n", num_diagonals);
	int l1 = strlen(seq1);
	int l2 = strlen(seq2);
	int gep = param_set->gep;

	int dig_length;
	if (seq1 > seq2)
		dig_length = l1;
	else
		dig_length = l2;

	int current_size = num_diagonals*dig_length + l1 +l2;

	int **list = vcalloc(current_size, sizeof(int*));
// 	int *diags = vcalloc(num_diagonals, sizeof(int));
	int i;
// 	for (i = 0; i < num_diagonals; ++i)
// 	{
// 		diags[i] = l1 - diagonals[i*3] + diagonals[i*3+1];
// 	}

// 	qsort (diags, num_diagonals, sizeof(int), fastal_compare);


// 	int *diagx = vcalloc(num_diagonals, sizeof(int));
// 	int *diagy = vcalloc(num_diagonals, sizeof(int));


	//+1 because diagonals start here at position 1, like in "real" dynamic programming
	int a = -1, b = -1;
	for (i = 0; i < num_diagonals; ++i)
	{
		if (diagonals[i].current_pos[1] - diagonals[i].current_pos[0] < 0)
		{
// 			diagx[i] = l1 - diags[i];
// 			diagy[i] = 0;
			a= i;
		}
		else
			break;
	}
	++a;
	b=a-1;
	for (; i < num_diagonals; ++i)
	{
// 		diagx[i] = 0;
// 		diagy[i] = diags[i]-l1;
		diagonals[i].prev_position = diagonals[i].current_pos[1] - diagonals[i].current_pos[0];
		b = i;
	}
	
// 	vfree(diags);
	int tmpy_pos;
	int tmpy_value;
	int **M = param_set->M;
	int *last_y = vcalloc(l2+1, sizeof(int));
	int *last_x = vcalloc(l1+1, sizeof(int));
	last_y[0] = 0;

	last_x[0] = 0;
	list[0] = vcalloc(6, sizeof(int));

	int list_pos = 1;
	int dig_num = l1;
	int tmp_l2 = l2 + 1;

	//left border
	for (; list_pos < tmp_l2; ++list_pos)
	{
		list[list_pos] = vcalloc(6, sizeof(int));
		list[list_pos][0] = 0;
		list[list_pos][1] = list_pos;
		last_y[list_pos] = list_pos;
		list[list_pos][2] = list_pos*gep;
		list[list_pos][4] = list_pos-1;
	}

	int pos_x = 0;
	int y;
	int tmp_l1 = l1-1;
	int tmpx, tmpy;
	while (pos_x < tmp_l1)
	{
		if (list_pos + num_diagonals+2 > current_size)
		{
			current_size += num_diagonals*1000;
			list = vrealloc(list, current_size * sizeof(int*));
		}
		//upper border
		list[list_pos] = vcalloc(6, sizeof(int));
		list[list_pos][0] = ++pos_x;
		list[list_pos][1] = 0;
		list[list_pos][2] = pos_x * gep;
		list[list_pos][3] = last_y[0];
// 		tmpy_value = list_pos;
// 		tmpy_pos = 0;
		last_x[pos_x] = list_pos;
		++list_pos;

		//diagonals
		for (i = a; i <= b; ++i)
		{
			
			if (pos_x == diagonals[i].current_pos[0])
			{
				if (diagonals[i].current_pos[2] == 0)
				{
					diagonals[i].current_pos+=3;
				}
				list[list_pos] = vcalloc(6, sizeof(int));
				tmpx = diagonals[i].current_pos[0];
				tmpy = diagonals[i].current_pos[1];
				list[list_pos][0] = tmpx;
				list[list_pos][1] = tmpy;
				list[list_pos][3] = last_y[tmpy];
				list[list_pos][4] = list_pos-1;
				list[list_pos][5] = diagonals[i].prev_position;//last_y[tmpy-1];
				diagonals[i].prev_position = list_pos;
// 				printf("A: %i %i\n",tmpx-1, tmpy-1);
				list[list_pos][2] = M[toupper(seq1[tmpx-1])-'A'][toupper(seq2[tmpy-1])-'A'];
				last_y[tmpy] = list_pos;
// 				tmpy_value = list_pos;
// 				tmpy_pos = tmpy;
				--diagonals[i].current_pos[2];
				++diagonals[i].current_pos[0];
				++diagonals[i].current_pos[1];
				++list_pos;
			}
		}
// 		last_y[tmpy_pos] = tmpy_value;

		
		//lower border
		if (list[list_pos-1][1] != l2)
		{
			list[list_pos] = vcalloc(6, sizeof(int));
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


		if ((b >= 0) && (diagonals[b].current_pos[1] > l2))
			--b;
		
		if ((a >0) && (diagonals[a-1].current_pos[0] - diagonals[a-1].current_pos[1] == pos_x))
		{
			--a;
			diagonals[a].prev_position = last_x[pos_x-1];
		}
	}


	dig_num = -1;
	if (list_pos + l2+2 > current_size)
	{
		current_size += list_pos + l2 + 2;
		list = vrealloc(list, current_size * sizeof(int*));
	}


// 	right border	
	list[list_pos] = vcalloc(6, sizeof(int));
	list[list_pos][0] = l1;
	list[list_pos][1] = 0;
	list[list_pos][3] = last_x[l1-1];
	list[list_pos][2] = -1000;
	++list_pos;
	
	

	for (i = 1; i <= l2; ++i)
	{
		list[list_pos] = vcalloc(6, sizeof(int));
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
// 	vfree(diagx);
// 	vfree(diagy);


	return list;
}




int
seq_pair2blast_diagonal2(char *seq_file_name1,
						 char *seq_file_name2,
						 Diagonal **diagonals,
						 int *dig_length,
						 int l1,
						 int l2,
						 int is_dna)
{
	char *out_file = vtmpnam(NULL);
	char blast_command[200];
	char blast_command2[200];

	if (is_dna)
	{
		sprintf(blast_command, "bl2seq -p blastn -i %s -j %s -D 1 -g F -o %s -S 1 -F F", seq_file_name1, seq_file_name2, out_file);
	}
	else
	{
		sprintf(blast_command, "bl2seq -p blastp -i %s -j %s -D 1 -g F -o %s -F F -S 1", seq_file_name1, seq_file_name2, out_file);
	}
	system(blast_command);

	Diagonal *diags = diagonals[0];
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
		pos_q = atoi(strtok(NULL, delims));
		strtok(NULL, delims);
		pos_d = atoi(strtok(NULL, delims));

		if (current_pos >= *dig_length)
		{
			(*dig_length) += 40;
			diags = vrealloc(diags, sizeof(Diagonal)*(*dig_length));
		}
		diags[current_pos].x = pos_q;
		diags[current_pos].y = pos_d;
		diags[current_pos++].length = length;
	}
	fclose(diag_f);
	diagonals[0] = diags;
	return current_pos;
}

