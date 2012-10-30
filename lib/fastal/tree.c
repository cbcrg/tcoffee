#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "ctype.h"


#include "io_lib_header.h"
#include "util_lib_header.h"
// #include "define_header.h"
// #include "dp_lib_header.h"
// #include "fastal_lib_header.h"



#include "fastal_lib_header.h"

/*!
 *	\file tree.c
 *	\brief Source code for the fastal tree algorithm
 */

/*


int
main(int argc, char** argv)
{
	int alphabet_size = 5;
	int aa[256];
	if (alphabet_size == 5)
	{
		aa['A'] = 0;
		aa['B'] = 1;
		aa['C'] = 1;
		aa['D'] = 0;
		aa['G'] = 2;
		aa['H'] = 0;
		aa['K'] = 3;
		aa['M'] = 1;
		aa['N'] = 0;
		aa['R'] = 0;
		aa['S'] = 1;
		aa['T'] = 3;
		aa['U'] = 4;
		aa['W'] = 3;
		aa['Y'] = 1;
	}
	else
	{
		aa['A'] = 0;
		aa['B'] = 1;
		aa['C'] = 2;
		aa['D'] = 3;
		aa['E'] = 4;
		aa['F'] = 5;
		aa['G'] = 6;
		aa['H'] = 7;
		aa['I'] = 8;
		aa['J'] = 9;
		aa['K'] = 10;
		aa['L'] = 11;
		aa['M'] = 12;
		aa['N'] = 13;
		aa['P'] = 14;
		aa['Q'] = 15;
		aa['R'] = 16;
		aa['S'] = 17;
		aa['T'] = 18;
		aa['V'] = 19;
		aa['W'] = 20;
		aa['X'] = 21;
		aa['Y'] = 22;
		aa['Z'] = 23;
	}
	compute_oligomer_distance_tree(argv[1], &aa[0], argv[2], 5, 2, alphabet_size);
	return 0;
}*/


void
compute_oligomer_distance_tree(char *seq_file_name, int* char2value, char *tree_file_name, int max_dist, int word_length, int alphabet_size)
{
	int i;
	int *compare = vcalloc(1,sizeof(int));
	int *num_seq = vcalloc(1,sizeof(int));

	int num_features = (int)pow(alphabet_size,word_length);
	
	num_features *= num_features * (max_dist+1);


	Cluster_info *matrix  = feature_extract(seq_file_name, max_dist, word_length, char2value, alphabet_size, compare, num_seq, num_features);
	
	int *node = vcalloc(1, sizeof(int));
	*node = *num_seq;
	
	FILE *tree_f = fopen(tree_file_name,"w");
	double *mean = vcalloc(num_features, sizeof(double));
	double *variance = vcalloc(num_features, sizeof(double));
	cluster_recursive(0, *num_seq-1, matrix, num_features, mean, variance, compare, tree_f, node);
	
	for (i = 0; i < *num_seq; ++i)
	{
// 		printf("ERG: %i\n", i);
		vfree(matrix[i].features);
	}
	vfree(matrix);
	fclose (tree_f);
	vfree(mean);
	vfree(variance);
	vfree(compare);
	vfree(node);
	vfree(num_seq);
}



/**
 * Recodes a given sequence into a sequence of it's k-mers.
 *
 * \param seq The sequence to encode.
 * \param seq_length The length of \a seq.
 * \param char2value Array giving each character in seq a unique value.
 * \param alphabet_size The size of the alphabet.
 * \param word_length The length of the k-mers to use.
 * 
 * \return The recodes sequence.
 */
int*
recode_sequence(char * seq, int seq_length, int *char2value, int alphabet_size, int word_length)
{

	int *recoded = vcalloc(seq_length - word_length + 1, sizeof(int));
	int i;
	
	if (word_length == 1)
	{
		for (i = 0; i < seq_length; ++i)
		{
			recoded[i] = char2value[(short)seq[i]];
		}
	}
	else
	{
		int *prod=vcalloc (word_length, sizeof(int));
		for ( i=0; i<word_length; ++i)
		{
			prod[word_length-i-1]=(int)pow(alphabet_size,i);
		}

		int index = 0;
		for (i = 0; i < word_length; ++i)
		{
			index += char2value[(short)seq[i]] * prod[i];
		}
		recoded[0] = index;
		int z = 0;
		for (; i < seq_length; ++i)
		{
			index -= char2value[(short)seq[z]] * prod[0];
			index *= alphabet_size;
			index += char2value[(short)seq[i]];
			recoded[++z] = index;
		}
	}

	return recoded;
}


/**
 * Extracts the features from a given sequence.
 *
 * \param recoded_seq The recoded sequence, consisting of k-mer codes.
 * \param seq_length The length of \a recoded_seq.
 * \param k_tup Length of the k-mers.
 * \param max_coding The maximal number used for coding the alphabet.
 * \param max_dist The maximal distance to use between to k-mers.
 * 
 * \return The feature vector of the sequence.
 */
double *
get_features(int *recoded_seq, int seq_length, int k_tup, int max_coding, int max_dist, int num_features)
{

	
	int vector_length = num_features;//max_coding * max_coding*(max_dist+1);
	
	double *feature_vector = vcalloc(vector_length, sizeof(double));
	int i;
	for (i = 0; i < vector_length; ++i)
	{
		feature_vector[i] = 0;
	}
	int j, max_length, indJ, indK, indDist;
	for (i = 0; i < seq_length; ++i)
	{
		if (seq_length <= i + max_dist)
			max_length = seq_length-1;
		else
			max_length = i + max_dist;
		for (j = i; j <= max_length; ++j)
		{
			
			indJ = max_coding * (max_dist+1) * recoded_seq[i];
			indK = (max_dist+1) * recoded_seq[j];
			indDist = j - i;
			++feature_vector[indJ + indK + indDist];
		}
	}

	double normalize = 0;
	for (i = 0; i < vector_length; ++i)
	{
		normalize += (feature_vector[i] * feature_vector[i]);
	}
	normalize = sqrt(normalize);
	for (i = 0; i < vector_length; ++i)
	{
		feature_vector[i] /= normalize; 
	}

	return feature_vector;
}


/**
 * Calculates the feature values for all sequences in a file.
 *
 * The File has to have fasta format.
 * \param seq_file_name The name of the file (in fasta format).
 * \param max_dist The maximal distance to be considered.
 * \param k_tup Size of the k_tup.
 * \param alphabet_size The size of the alphabet.
 * 
 * \return The feature values for every sequence in a Cluster_info object.
 */
Cluster_info *
feature_extract(char *seq_file_name, int max_dist, int k_tup, int *char2value, int alphabet_size, int * elem_2_compare, int *num_seq_p, int num_features)
{
	char line[500];
	FILE *tree_f = fopen(seq_file_name,"r");
// 	fgets(line, 500, tree_f);
	int num_seq = -1;
	int size = 1000;
	int matrix_size = 1000;
	const int STEP = 1000;
	int seq_pos = -1;
	char *seq = vcalloc(size, sizeof(char));
	Cluster_info *matrix = vcalloc(matrix_size, sizeof(Cluster_info));
	char *c_p;
	int max_coding = pow(alphabet_size,k_tup);
	while (fgets(line, 500, tree_f) != NULL)
	{
		if (line[0] == '>')
		{
			if (num_seq >= 0)
			{
				seq[++seq_pos] = '\0';
				if (num_seq > matrix_size -2)
				{
					matrix_size += STEP;
					matrix = vrealloc(matrix, matrix_size * sizeof(Cluster_info));
				}
				int *recoded_seq = recode_sequence(seq, seq_pos, char2value, alphabet_size, k_tup);
				matrix[num_seq].seq_number = num_seq;
				matrix[num_seq].elem_2_compare = elem_2_compare;
				matrix[num_seq].features=get_features(recoded_seq, seq_pos- k_tup +1, k_tup, max_coding, max_dist, num_features);
				vfree(recoded_seq);
				seq_pos = -1;
			}
			++num_seq;
		}
		else
		{
			if (size - seq_pos < 500)
			{
				size += STEP;
				seq = vrealloc(seq, size*sizeof(char));
			}
			c_p = &line[0];
			while ((*c_p != '\n') && (*c_p != '\0'))
			{
				seq[++seq_pos] = toupper(*(c_p));
				++c_p;
			}
		}
	}
	
	fclose(tree_f);
	seq[++seq_pos] = '\0';
	int *recoded_seq = recode_sequence(seq, seq_pos, char2value, alphabet_size, k_tup);
	matrix[num_seq].seq_number = num_seq;
	matrix[num_seq].elem_2_compare = elem_2_compare;
	matrix[num_seq].features=get_features(recoded_seq, seq_pos- k_tup +1, k_tup, max_coding, max_dist, num_features);
	*num_seq_p = ++num_seq;
	vfree(seq);
	vfree(recoded_seq);
	return matrix;
}



/**
 * Compare function for qsort given an array of Cluster_info.
 *
 * The field used for sorting is done according to the element determined in the Cluster_info object.
 * \param a void pointer to an object of Cluster_info.
 * \param b void pointer to an object of Cluster_info.
 *
 * \return Value showing wheather \a a is bigger, smaller or equal to \a b.
 */
int 
cluster_compare (const void * a, const void * b)
{
	Cluster_info *tmp1 = (Cluster_info*)a;
	Cluster_info *tmp2 = (Cluster_info*)b;
	int elem = *(tmp1->elem_2_compare);
	if ((tmp1->features[elem]) > (tmp2->features[elem]))
		return 1;
	else if ((tmp1->features[elem]) < (tmp2->features[elem]))
		return -1;
	else
		return 0;
}



/**
 * Recursive function to build a tree according to a given feature matrix.
 *
 * \param start The begin of this set in \a matrix.
 * \param end The end of this set in \a matrix.
 * \param matrix The feature matrix.
 * \param dim The number of features.
 * \param mean Array of size \a dim.
 * \param variance Array of size \a dim.
 * \param elem_2_compare Pointer to a single in. This will save the feature according to which shall be sorted (biggest variance).
 * \param tree_f The file in which the tree will be written.
 * \param node The current node.
 * 
 * \return The number of the current node.
 */
int
cluster_recursive(int start,
				  int end,
				  Cluster_info* matrix,
				  int dim,
				  double *mean,
				  double *variance,
				  int *elem_2_compare,
				  FILE* tree_f,
				  int *node)
{

	//stop recursion
	if (end-start < 1)
	{
		return matrix[start].seq_number;
	}

	//reset mean/variance
	int i;
	int num_seq = end - start +1;
	for (i = 0; i < dim; ++i)
	{
		mean[i] = 0;
		variance[i] = 0;
	}

	//calculate mean
	Cluster_info *start_p = &matrix[start]-1;
	Cluster_info   *end_p = &matrix[end];
	double *tmp;
	while (start_p < end_p)
	{
		tmp = &(++start_p)->features[0];
		for (i = 0; i < dim; ++i)
		{
			mean[i] += *(tmp++);
		}
	}
	for (i = 0; i < dim; ++i)
			mean[i] /= num_seq;


	//calculate variance
	start_p = &matrix[start] -1;

	while (start_p < end_p)
	{
		tmp = &(++start_p)->features[0];
		for (i = 0; i < dim; ++i)
		{
			variance[i] +=  ((*tmp) - mean[i]) * ((*tmp) - mean[i]);
			++tmp;
		}
	}

	//get maximal variance
	double max  = variance[0];
	int index  = 0;
	for (i = 1; i < dim; ++i)
	{
		if (variance[i] > max)
		{
			max = variance[i];
			index = i;
		}
	}
	*elem_2_compare = index;


	//devide into to different sets and start recursion on these child sets
	qsort (&matrix[start], num_seq, sizeof(Cluster_info), cluster_compare);
	int split = start + (end-start)/2;
	int split_value = matrix[split].features[index];
	while ((matrix[split+1].features[index] == split_value) && (split < end-1))
		++split;

	int left_child = cluster_recursive(start, split, matrix, dim, mean, variance, elem_2_compare, tree_f, node);
	int right_child = cluster_recursive(split+1, end, matrix, dim, mean, variance, elem_2_compare, tree_f, node);
	int node2 = *node;
	++(*node);
	//print node into file
	fprintf(tree_f, "%i %i %i\n", left_child, right_child, node2);

	return node2;
}

