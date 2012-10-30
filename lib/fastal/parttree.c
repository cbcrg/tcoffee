#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
// #include <argp.h>

#include "fastal_lib_header.h"
#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"
#include "dp_lib_header.h"
//TODO: -change kick-out value
//		-pass arrays in partTree_r

/*!
 *	\file parttree.c
 *	\brief Source code for PartTree algorithm.
 */



#define ENLARGEMENT_PER_STEP 50



void
print_fastal_tree(Tree_fastal *tree,
				  int pos,
				  FILE *tree_file,
				  int num_seq)
{
	if (tree[pos].left >= num_seq)
		print_fastal_tree(tree, tree[pos].left-num_seq, tree_file, num_seq);
	if (tree[pos].right >= num_seq)
		print_fastal_tree(tree, tree[pos].right-num_seq, tree_file, num_seq);
	
	fprintf(tree_file, "%i %i %i\n", tree[pos].left, tree[pos].right, tree[pos].name);
}





PartTree_param *param_set;


//**********************   UPGMA *****************************

/**
 * Function to write tree to file in fastal_format. Leafs in \a tree are leafs in the complete tree as well.
 *
 * \param tree The root node of the (sub)tree.
 * \param param_set Parameter of PartTree.
 * \param start_write_here Current node to write into.
 * \return position in tree
 * \see tree_process
 */
int
tree_process_simple(NT_node tree,
				    PartTree_param *param_set,
					int start_write_here)
{
	if (tree->isseq)
	{
// 		printf("T: %s\n", tree->name);
		return atoi(tree->name);
	}
	else
	{
		Tree_fastal *tree_flat = &param_set->tree[start_write_here];
		tree_flat->name = start_write_here +param_set->num_sequences;
		if (start_write_here == param_set->pos_tree)
		{
			++param_set->pos_tree;
		}
		start_write_here = param_set->pos_tree;
		int left = tree_process_simple(tree->left, param_set, start_write_here);
		start_write_here = param_set->pos_tree;
		int right = tree_process_simple(tree->right, param_set, start_write_here);
		tree_flat->index = NULL;
		tree_flat->right = right;
		tree_flat->left = left;
		return tree_flat->name;
	}
}


void
read_sequence_from_position(FILE *file, long position1, char *seq)
{
	fseek(file, position1, SEEK_SET);
	char line[500];
	int pos = -1;
	while ((fgets(line, 500, file) != NULL) && (line[0] != '>'))
	{
		int l = -1;
		while ((line[++l] != '\n') && (line[l] != '\0'))
			seq[++pos] = line[l];
	}
	seq[++pos] = '\0';
}

// char*
// read_sequence_from_position(FILE *file, long position1, long position2)
// {
// 	
// }

/**
* Function to write tree to file in fastal_format. Leafs in \a tree do not need to be leafs in the complete tree.
*
* \param tree The root node of the (sub)tree.
* \param param_set Parameter of PartTree.
* \param clusters Number of sequences in each cluster.
* \param subgroup The sequences for each cluster.
* \param start_write_here Current node to write into.
* \return position in tree
* \see tree_process_simple
*/
int
tree_process(NT_node tree,
			 PartTree_param *param_set,
			 int *clusters,
		     int *subgroup,
			 int start_write_here)
{
	if (tree->isseq)
	{
		int node_num = atoi(tree->name);
		int num_in_sub = clusters[node_num+1] - clusters[node_num];
// 		printf("NUM: %i %i %i %i\n",node_num, num_in_sub, clusters[node_num+1], clusters[node_num]);
		if (num_in_sub > 1)
		{
			Tree_fastal *tree_flat = &param_set->tree[start_write_here];
			tree_flat->name = start_write_here +param_set->num_sequences;
			if (start_write_here == param_set->pos_tree)
			{
				++param_set->pos_tree;
			}
			tree_flat->left  = -1;
			tree_flat->right = -1;
			tree_flat->index = &subgroup[clusters[node_num]];

			tree_flat->num_leafs = num_in_sub;
			return tree_flat->name;
		}
		else
		{
			return(subgroup[clusters[node_num]]);
		}
	}
	else
	{
// 		printf("TREEPOS: %i\n",param_set->pos_tree);
		Tree_fastal *tree_flat = &param_set->tree[start_write_here];
		tree_flat->name = start_write_here +param_set->num_sequences;
		if (start_write_here == param_set->pos_tree)
		{
			++param_set->pos_tree;
		}
		start_write_here = param_set->pos_tree;
		int left = tree_process(tree->left, param_set, clusters, subgroup, start_write_here);
		start_write_here = param_set->pos_tree;
		int right = tree_process(tree->right, param_set, clusters, subgroup, start_write_here);
		tree_flat->index = NULL;
		tree_flat->right = right;
		tree_flat->left = left;
		return tree_flat->name;
	}
}


/**
* \brief Calculates tree out of distance matrix.
* 
* Calculates the upgma tree using a given distance matrix.
* \param mat The distance matrix.
* \param nseq Number of sequences.
* \param fname Filename for temporary storage.
* \param seqnames Names of the sequences.
* \return The calculated UPGMA Tree.
*/
NT_node ** int_dist2upgma_tree_fastal (int **mat, int nseq, char *fname, char **seqnames)
{
	NT_node *NL, T;
	int a, n, *used;
	int tot_node;
	if (upgma_node_heap (NULL))
	{
		printf_exit ( EXIT_FAILURE,stderr, "\nERROR: non empty heap in upgma [FATAL]");
	}
	NL=vcalloc (nseq, sizeof (NT_node));

	for (a=0; a<nseq; a++)
	{
		NL[a]=new_declare_tree_node ();
		upgma_node_heap (NL[a]);
		sprintf (NL[a]->name, "%s", seqnames[a]);
		NL[a]->isseq=1;
		NL[a]->leaf=1;
	}
	used=vcalloc ( nseq, sizeof (int));
	n=nseq;
	while (n>1)
	{
		T=upgma_merge(mat, NL,used, &n, nseq);
	}
	vfree (used);
	vfclose (print_tree (T, "newick", vfopen (fname, "w")));
	upgma_node_heap (NULL);
	vfree (NL);

	return read_tree (fname,&tot_node, nseq, seqnames);
}




//Part_Tree

/*!
* \brief Constructs a guide tree for multiple sequence alignment.
*
* This algorithm is an implementation of the partTree algorithm (PartTree: an algorithm to build an approximate tree from a large number of unaligned
* sequences. Katoh et al. 2007).
* \param sequence_f Filename of file with sequences.
* \param tree_f Filename of file where the tree will be stored.
* \param ktup Size of the ktups.
* \param subgroup Parameter for subgroupsize.
*/
void
make_partTree(char *sequence_f,
			  char *tree_f,
			  int ktup,
			  int subgroup,
			  int is_dna,
			  int retree)
{
	long *file_positions = NULL;
	long **tmp1 = &file_positions;
	int *seq_lengths = NULL;
	int **tmp2 = &seq_lengths;
	int number_of_sequences;
	param_set = vcalloc(1,sizeof(PartTree_param));
	param_set->ktup = ktup;
	param_set->subgroup = subgroup;
	Tree_fastal *tree;

	
	if (!retree)
	{
		//make index
		char *ktup_table = vtmpnam(NULL);
		param_set->num_sequences = number_of_sequences =  make_pos_len_index_of_file(sequence_f, ktup_table, tmp1, tmp2, ktup, is_dna);
		tree = vcalloc(number_of_sequences-1,sizeof(Tree_fastal));
		param_set->tree = tree;
		param_set->ktup_positions = file_positions;
		param_set->seq_lengths = seq_lengths;
		param_set->threshold = 0.01;
		param_set->ktup_table_f = fopen(ktup_table,"r");

	
		partTree_r(param_set);
	}
	else
	{

		//make index
		param_set->num_sequences = number_of_sequences =  make_pos_len_index_of_file_retree(sequence_f, tmp1, tmp2); 
		tree = vcalloc(number_of_sequences-1,sizeof(Tree_fastal));
		param_set->tree = tree;
		param_set->ktup_positions = file_positions;
		param_set->seq_lengths = seq_lengths;
		param_set->threshold = 0.01;
		param_set->ktup_table_f = fopen(sequence_f,"r");

		partTree_retree(param_set);

	}
	
	FILE * tree_file = fopen(tree_f,"w");
	print_fastal_tree(tree, 0, tree_file, number_of_sequences);
	fclose(tree_file);
	vfree(tree);
// 	exit(1);
}


/**
* Filters seed set.
*
* \param sequence_group Sequences to filter.
* \param dist_mat The distance matrix.
* \param seed_set_cleaned ordered_seed_set.
* \param param_set Parameters for PartTree algorithm.
* \return number in the filtered set.
*/
int
filter(int *sequence_group,
	   int **dist_mat,
	   int *seed_set_cleaned,
	   PartTree_param *param_set)
{
	int i, j;
	int num_in_subgroup = param_set->subgroup;
	int *seq_lengths = param_set->seq_lengths;
	int num_in_clean = 0;
	double threshold = param_set->threshold;
// 	printf("threshold: %f\n", threshold);
	double min;
	for (i = 0; i < num_in_subgroup; ++i)
	{
		if (!seed_set_cleaned[i])
			continue;
		for (j = i+1; j < num_in_subgroup; ++j)
		{
			min = MIN(seq_lengths[sequence_group[i]], seq_lengths[sequence_group[j]]);
// 			printf("MINIMUM: %i\n",min);
			min = (threshold * min);
// 			printf("MINIMUM: %f\n",min);
			if (seed_set_cleaned[j] &&(dist_mat[i][j] < min))
			{
				if (seq_lengths[sequence_group[i]] < seq_lengths[sequence_group[j]])
				{
					seed_set_cleaned[i] = 0;
					break;
				}
				else
					seed_set_cleaned[j] = 0;
			}
		}
	}

	for (i = 0; i < num_in_subgroup; ++i)
	{
		num_in_clean += seed_set_cleaned[i];
	}
	int max = num_in_subgroup -1;
	i = 0;
	int tmp;
// 	printf("CLEAN: %i\n", num_in_clean);
	while (i < num_in_clean)
	{
		if (seed_set_cleaned[i])
		{
			++i;
		}
		else
		{
			seed_set_cleaned[i] = seed_set_cleaned[max];
			seed_set_cleaned[max] = 0;
			tmp = sequence_group[i];
			sequence_group[i] = sequence_group[max];
			sequence_group[max] = tmp;
			--max;
		}
	}
	return num_in_clean;
}





/*!
*	\brief Function to create a tree using the PartTree algorithm.
*	
*	\param param_set A \a PartTree_param object containing all necessary parameters and the data.
*	\return The node_number.
*/
void
partTree_r(PartTree_param *param_set)
{
	
	int num_of_tree_nodes = param_set->num_sequences-1;
	int loop_tree_node;

	Tree_fastal *tree = param_set->tree;
// 	int this_node = param_set->pos_tree;
	
	int i;
	int tsize = param_set->tsize;
	
	
	//get some memory
	short *table1 = vcalloc(tsize, sizeof(short));
	short *table2 = vcalloc(tsize, sizeof(short));
	char **names = declare_char(param_set->subgroup, 8);
	int **dist_mat = declare_int(param_set->subgroup, param_set->subgroup);
	int **dist_mat2 = declare_int(param_set->subgroup, param_set->subgroup);
	char * file_name_tmp = vtmpnam(NULL);
	int *seed_set_cleaned = vcalloc(param_set->subgroup, sizeof(int));
	FILE *table_f = param_set->ktup_table_f;
	long *file_positions = param_set->ktup_positions;
	int max_n_group = param_set->subgroup;
	int num_in_subgroup = param_set->subgroup;
	int *seq_lengths = param_set->seq_lengths;
	int *clusters = vcalloc(param_set->subgroup+1, sizeof(int));
	int *min_dist = vcalloc(param_set->num_sequences, sizeof(int));
	int *belongs_to = vcalloc(param_set->num_sequences, sizeof(int));
	
	
	
	
	
	
	//Prepare first node
	
	tree[0].index = vcalloc(param_set->num_sequences,sizeof(int));
	int *index = tree[0].index;
	for (i = 0; i< param_set->num_sequences; ++i)
		index[i] = i;
	tree[0].name = param_set->pos_tree +param_set->num_sequences;
	
	tree[0].num_leafs = param_set->num_sequences;
	int *sequence_group2 = vcalloc(param_set->num_sequences,sizeof(int));
	
	Tree_fastal *current_node;
	for (loop_tree_node = 0; loop_tree_node < num_of_tree_nodes; ++loop_tree_node)
	{
// 		printf("ROUND: %i\n", loop_tree_node);
		current_node = &tree[loop_tree_node];
		index= current_node->index;
		if (current_node->index == NULL)
		{
			continue;
		}
		int num_sequences = current_node->num_leafs;

		//if number of sequences in this group smaller than number subgoup size: make tree, finisch
		if (num_sequences <= max_n_group)
		{
			dist_mat = make_distance_matrix(table_f, file_positions, index, num_sequences, dist_mat);
			for (i = 0; i < num_sequences; ++i)
			{
				sprintf(names[i],"%i", current_node->index[i]);
			}
			NT_node **tree= (int_dist2upgma_tree_fastal (dist_mat, num_sequences, file_name_tmp , names));
			tree_process_simple(tree[0][0], param_set,loop_tree_node);
			continue;
		}

		
		for (i = 0; i < num_in_subgroup; ++i)
		{
			seed_set_cleaned[i] = 1;
		}
		
		//finde longest sequence and put into the first field
		
		int index_longest = 0;
		int length_of_longest = 0;
		
		for(i = 0; i < num_sequences; ++i)
		{
			if (seq_lengths[index[i]] > length_of_longest)
			{
				index_longest = i;
				length_of_longest = seq_lengths[index[i]];
			}
		}
		int tmp = index[index_longest];
		index[index_longest] = index[0];
		index[0] = tmp;

		//distance of longest to rest
		int seq_index = 1;
		int min= euclidean_dist(table_f, file_positions[index[0]], file_positions[index[1]], table1, table2, param_set->tsize);
		for (i = 2; i < num_sequences; ++i)
		{
			tmp = euclidean_dist_half(table_f, file_positions[index[i]], table1, table2, param_set->tsize);
			if (tmp < min)
			{
				min = tmp;
				seq_index = i;
			}
		}

		//get the new seed_set in the first n spaces
		tmp = index[1];
		index[1] = index[seq_index];
		index[seq_index] = tmp;
		int r,j;
		num_in_subgroup = param_set->subgroup;

		
		for (i = 2; i < num_in_subgroup; ++i)
		{
			r = i + rand() / ( RAND_MAX / ( num_sequences-i) + 1 );
// 			printf("RANDOM: %i\n",r);
			tmp = index[r];
			index[r] = index[i];
			index[i] = tmp;
		}

		//Calculate matrix
		dist_mat = make_distance_matrix(table_f, file_positions, index, param_set->subgroup, dist_mat);
	
		//Filter out sequences that are to similar & reorder
		
		NT_node **upgma_tree;
		
		
		int num_in_clean = filter(index, dist_mat, seed_set_cleaned, param_set);

		
		if (num_in_clean ==1)
		{
			num_in_clean = 2;
			seed_set_cleaned[1] = 1;
		}
			//make_tree
		int col = 0;
		int row = 0;
		for (i = 0; i <  num_in_subgroup; ++i)
		{
			if (seed_set_cleaned[i])
			{
				row = col+1;
				for (j = i+1; j <  num_in_subgroup; ++j)
				{
					if (seed_set_cleaned[j])
					{
						dist_mat2[row][col] = dist_mat2[col][row] = dist_mat[i][j];
						++row;
					}
				}
				++col;
			}
		}
		for (i = 0; i < num_in_clean; ++i)
		{
			sprintf(names[i],"%i",i);
		}
		upgma_tree= (int_dist2upgma_tree_fastal (dist_mat2, num_in_clean, file_name_tmp , names));
 

		//cluster
		//calculate distances from n' to N
		get_table(table1, table_f, file_positions[index[0]]);
		for (j = num_in_clean; j < num_sequences; ++j)
		{
			min_dist[j] = euclidean_dist_half(table_f, file_positions[index[j]], table1, table2, param_set->tsize);
			belongs_to[j] = 0;
		}
		for(i = 1; i < num_in_clean; ++i)
		{
			get_table(table1, table_f, file_positions[index[i]]);
			belongs_to[i] = i;
			for (j = num_in_clean; j < num_sequences; ++j)
			{
				tmp = euclidean_dist_half(table_f, file_positions[index[j]], table1, table2, param_set->tsize);
				if (tmp < min_dist[j])
				{
					min_dist[j] = tmp;
					belongs_to[j] = i;
				}
			}
		}

		//how_many sequences has each cluster
		for (j = 0; j <= num_in_subgroup; ++j)
		{
			clusters[j] = 0;
		}
		for (j = 0; j < num_sequences; ++j)
		{
			++clusters[belongs_to[j]];
		}
// 		for (j = 0; j <= num_in_subgroup; ++j)
// 		{
// 			printf("CL: %i ",clusters[j]);
// 		}
// 		printf("\n");
		for(i = 1; i < num_in_clean; ++i)
		{
			clusters[i] += clusters[i-1];
		}
		clusters[num_in_clean] = clusters[num_in_clean-1];

		for (i = 0; i < num_sequences; ++i)
		{
			sequence_group2[--clusters[belongs_to[i]]] = index[i];
		}

		for (i = 0; i < num_sequences; ++i)
		{
			index[i] = sequence_group2[i];
		}

		
		for (i = 0; i < num_in_clean; ++i)
		{
			sprintf(names[i],"%i",i);
		}
		tree_process(upgma_tree[0][0], param_set, clusters, index, loop_tree_node);
		NT_node tmp_tree = upgma_tree[3][0];
		vfree(upgma_tree[0]);
		vfree(upgma_tree[1]);
		vfree(upgma_tree[2]);
		vfree(upgma_tree[3]);
		vfree(upgma_tree);
		free_tree(tmp_tree);
	}
	vfree(min_dist);
	vfree(belongs_to);
	vfree(clusters);
}



/*!
 *	\brief Makes the distance matrix between all sequences.
 *
 *	\param table_file File with the ktup tables
 *	\param file_positions Index of positions where the tabels are stored in \a table_file
 *	\param sequence_group the group of sequences
 *	\param number number of sequences
 *	\param dist_mat distance matrix
 *	\return the distance matrix. (same as \a dist_mat )
*/
int **
make_distance_matrix(FILE *table_f,
					 long *file_positions,
					 int *sequence_group,
					 int number,
					 int **dist_mat)
{
	static short *table1 = NULL;
	static short *table2;
	int tsize = param_set->tsize;
	if (table1 == NULL)
	{
		table1 = vcalloc(tsize, sizeof(short));
		table2 = vcalloc(tsize, sizeof(short));
	}
	int i, j, num = number-1;
	for (i = 0; i < num; ++i)
	{
		j = i+1;
		dist_mat[i][j] = dist_mat[j][i]= euclidean_dist(table_f, file_positions[sequence_group[i]], file_positions[sequence_group[j]], table1, table2, tsize);
		++j;
		for (; j < number; ++j)
		{
			dist_mat[i][j] = dist_mat[j][i] = euclidean_dist_half(table_f, file_positions[sequence_group[j]], table1, table2, tsize);
		}
	}
	return dist_mat;
}


int **
make_distance_matrix_sim(FILE *aln_f,
						 long *file_positions,
						 int *sequence_group,
						 int number,
						 int **dist_mat)
{
	static char *seq1 = NULL;
	static char *seq2;
	char line[500];
	if (seq1 == NULL)
	{
		int length = 0;
		int i;
		fseek(aln_f, file_positions[0], SEEK_SET);
		fgets(line, 500, aln_f);
		while (line[0] != '>')
		{
			i = -1;
			while ((line[++i] != '\n') && (line[i] != '\0'));
			length += i;
			fgets(line, 500, aln_f);
		}
		seq1 = vcalloc(length+1, sizeof(short));
		seq2 = vcalloc(length+1, sizeof(short));
	}


	int i, j, num = number-1, pos;
	for (i = 0; i < num; ++i)
	{

		fseek(aln_f, file_positions[sequence_group[i]], SEEK_SET);
  		
		
		pos = -1;
		while ((fgets(line, 500, aln_f) != NULL) && (line[0] != '>'))
		{
			int l = -1;
			while ((line[++l] != '\n') && (line[l] != '\0'))
				seq1[++pos] = line[l];
		}
		seq1[++pos] = '\0';

		for (j = i+1; j < number; ++j)
		{
			pos = -1;
			fseek(aln_f, file_positions[sequence_group[j]], SEEK_SET);
			while ((fgets(line, 500, aln_f) != NULL) && (line[0] != '>'))
			{
				int l = -1;
				while ((line[++l] != '\n') && (line[l] != '\0'))
					seq2[++pos] = line[l];
			}
			seq2[++pos] = '\0';
			dist_mat[i][j] = dist_mat[j][i] = 100 - fast_aln2sim_mat2(seq1, seq2);
		}
	}
	return dist_mat;
}


/**
* Replaces the coded sequence with coded tuples
*
* \param coded_seq The coded sequence which will be replaced by the tuple number
* \param ktup Size of the ktup
* \param ng Coded alphabet size
* \param length Lengths of coded sequence
*/
void
makepointtable_fast(int *coded_seq,	//sequence
					int ktup,		//ktup size
					int ng,			//hmm...
					int length)		//length of coded_seq
{
	int point, a;
	register int *p;
	static int *prod;

	if (!prod)
	{
		prod=vcalloc ( ktup, sizeof (int));
		for ( a=0; a<ktup; a++)
		{
			prod[ktup-a-1]=(int)pow(ng,a);
		}
	}
	p = coded_seq;

  
	for (point=0,a=0; a<ktup; a++)
	{
		point+= *coded_seq++ *prod[a];
	}

	int i = ktup;
	while (i < length)
	{
		point -= *p * prod[0];
		point *= ng;
		point += *coded_seq;
		*p = point;
		++p;
		++coded_seq;
		++i;
	}
	*p = END_ARRAY;
}


/**
* \brief Calculates the number of occurences for each ktup.
*
* \param tables_f File to save the tables in.
* \param table Table to save the result in.
* \param pointt Array with all ktups listed one after another. Has to end with END_ARRAY.
* \param length length of \a table
*/
void
makecompositiontable_fastal(FILE* tables_f,	//File to save the tables in
							int *table,		//table to calculate in
							int *pointt,	//ktups array
							int length)		//length of the table
{
	int point;
	while(*pointt != END_ARRAY )
	{
		++table[*pointt];
		++pointt;
	}
	for (point = 0; point < length; ++point)
	{
		if (table[point] > 0)
			fprintf(tables_f, "%i %i\n", point, table[point]);
	}
	fprintf(tables_f, "*\n");
}


/** JUST FOR TEST */
void 
make_fast_tree(char *file_name,
			   int n,
			   int ktup)
{
	
	make_partTree(file_name, "TREE_OUT", ktup, n, 1, 0);

}



/**
* \brief Reads ktup_table from file
*
* \param table	Table to save the file content in.
* \param tables_f File in which the tables are stored.
* \param index Position of the table in \a tables_f
*/
void
get_table(short *table,		//Table to save the readings in
		  FILE* tables_f,	//File with tables
		  long index)		//index positin of ktup-tables
{
	fseek(tables_f, index, SEEK_SET);
	const int LINE_LENGTH = 101;
	char line[LINE_LENGTH];
	fgets(line, LINE_LENGTH, tables_f);
	
	char delims[] = " ";
	char *result = NULL;
	int code;

	while (line[0] != '*')
	{
		result = strtok( line, delims );
		code = atoi(result);
		table[code] = atoi(strtok( NULL, delims));
		fgets(line, LINE_LENGTH, tables_f);
	}
}



/**
* \brief calculates the euclidean ktub distance between two sequences
*
* @param ktup_f, ktup_file
* @param pos1 position of sequence 1 in \a ktup_f
* @param pos2 position of sequence 2 in \a ktup_f
* @param table1 Saves the number of occurences for each ktup in sequence 1
* @param table2 Saves the number of occurences for each ktup in sequence 2
*/
int 
euclidean_dist(FILE* ktup_f,	//ktup_file
				 long pos1,		//position of table1
				 long pos2,		//position of table2
				 short *table1,	//table to save ktups in
				 short *table2,	//table to save ktups in
				 int length)	
{
	const int LINE_LENGTH = 101;
	char line[LINE_LENGTH];
	
	
	char delims[] = " ";
	char *result = NULL;
	int code;

	fseek(ktup_f, pos1, SEEK_SET);
	fgets(line, LINE_LENGTH, ktup_f);
	int i;
	for (i = 0; i < length; ++i)
	{
		table1[i] = 0;
		table2[i] = 0;
	}
	while (line[0] != '*')
	{
		result = strtok( line, delims );
		code = atoi(result);
		table1[code] = atoi(strtok( NULL, delims));
		fgets(line, LINE_LENGTH, ktup_f);
	}
	fseek(ktup_f, pos2, SEEK_SET);
	fgets(line, LINE_LENGTH, ktup_f);
	while (line[0] != '*')
	{
		result = strtok( line, delims );
		code = atoi(result);
		table2[code] = atoi(strtok( NULL, delims));
		fgets(line, LINE_LENGTH, ktup_f);
	}

	int dist = 0;
	for (i = 0; i < length; ++i)
	{
		dist += (table1[i]-table2[i])*(table1[i]-table2[i]);
	}
	return dist;
}



/**
 * \brief calculates the euclidean ktub distance between two sequences.
 *
 * The difference to \a euclidean_dist is, that this uses the ktups stored in \a table1
 * @param ktup_f, ktup_file
 * @param pos2 position of sequence 2 in \a ktup_f
 * @param table1 Saves the number of occurences for each ktup in sequence 1
 * @param table2 Saves the number of occurences for each ktup in sequence 2
 * \see euclidean_dist
 */
int 
euclidean_dist_half(FILE* ktup_f,	//ktup_file
					long pos2,		//position of table1
					short *table1,	//table to save ktups in
					short *table2,	//table to save ktups in
					int length)
{
	const int LINE_LENGTH = 101;
	char line[LINE_LENGTH];
	
	
	char delims[] = " ";
	char *result = NULL;
	int code;

	fseek(ktup_f, pos2, SEEK_SET);
	fgets(line, LINE_LENGTH, ktup_f);
	int i;
	for (i = 0; i < length; ++i)
	{
		table2[i] = 0;
	}
	while (line[0] != '*')
	{
		result = strtok( line, delims );
		code = atoi(result);
		table2[code] = atoi(strtok( NULL, delims));
		fgets(line, LINE_LENGTH, ktup_f);
	}

	int dist = 0;
	for (i = 0; i < length; ++i)
	{
		dist += (table1[i]-table2[i])*(table1[i]-table2[i]);
	}
	return dist;
}




/**
*	Makes an index of a file
*/
int
make_pos_len_index_of_file(char *file_name,			//file with sequences
						   char *ktable_f,			//file with the ktup-tables
						   long **file_positions,	//array to save the positions
						   int **seq_lengths,		//array to save the sequence length
						   int ktup,				//length of ktup
						   int is_dna)				//type of the seuqence	
{
	//preparations for recoding sequence
	int *aa;
	int a, b;
	
	int ng = 0;
	char **gl;
	if (is_dna)
	{
		gl=declare_char (5,13);
		sprintf ( gl[ng++], "Aa");
		sprintf ( gl[ng++], "Gg");
		sprintf ( gl[ng++], "TtUu");
		sprintf ( gl[ng++], "Cc");
		sprintf ( gl[ng++], "NnRrYyDdMmWw");
	}
	else
	{
		gl=make_group_aa ( &ng, "mafft");
	}
	aa=vcalloc ( 256, sizeof (int));
	for ( a=0; a<ng; a++)
	{
		for ( b=0; b< strlen (gl[a]); b++) 
		{
			aa[(int)gl[a][b]]=a;
		}
	}
	free_char (gl, -1);

	
	int tsize=(int)pow(ng, ktup);
	param_set->tsize = tsize;
	param_set->ng = ng;
	
	int *table=vcalloc ( tsize,sizeof (int));


	//Reading and recoding squences
	const int LINE_LENGTH = 501;
	int *coded_seq = vcalloc(2*LINE_LENGTH, sizeof(int));
	int allocated_mem = 2*LINE_LENGTH;
			
	(*file_positions) = vcalloc(ENLARGEMENT_PER_STEP,  sizeof(long));
	(*seq_lengths) = vcalloc(ENLARGEMENT_PER_STEP,  sizeof(int));


	FILE *file = fopen(file_name,"r");

	char line[LINE_LENGTH];

	int num_of_sequences = 0;
	int str_len = 0;
	int mem_for_pos = ENLARGEMENT_PER_STEP;

	int real_len;
	int *c_seq;

	FILE *tables_f = fopen(ktable_f, "w");


	if (file == NULL)
	{
		printf("FILE NOT FOUND\n");
		exit(1);
	}
	else
	{

		while(fgets(line, LINE_LENGTH , file)!=NULL)
		{
			if ( str_len >= allocated_mem - LINE_LENGTH)
			{
				allocated_mem += LINE_LENGTH;
				coded_seq = vrealloc(coded_seq, allocated_mem*sizeof(int));
			}
			
			int i;
			
			if (line[0] == '>')
			{
				if (num_of_sequences >0)
				{
					(*seq_lengths)[num_of_sequences-1] = str_len;
// 					printf("len: %i\n", str_len);
					c_seq = coded_seq;
					makepointtable_fast(coded_seq,ktup,ng, str_len);
					
					(*file_positions)[num_of_sequences-1] = ftell(tables_f );
					for (i=0; i < tsize; ++i)
						table[i] = 0;
					makecompositiontable_fastal(tables_f, table, coded_seq,tsize );


				}
				str_len = 0;
				++num_of_sequences;

				if (num_of_sequences == mem_for_pos)
				{
					mem_for_pos += ENLARGEMENT_PER_STEP;
					(*file_positions) = vrealloc((*file_positions), mem_for_pos * sizeof(long));
					(*seq_lengths) = vrealloc((*seq_lengths), mem_for_pos * sizeof(int));
				}
			}
			else
			{
				int i;
				real_len = strlen(line);
				if (line[real_len-1] == '\n')
					--real_len;
				for (i = 0; i < real_len; ++i)
				{
					coded_seq[str_len++] = aa[(short)line[i]];
				}
			}
		}
	}

	(*seq_lengths)[num_of_sequences-1] = str_len;
	c_seq = coded_seq;
	makepointtable_fast(coded_seq,ktup,ng, str_len);
	(*file_positions)[num_of_sequences] = ftell(tables_f );
	makecompositiontable_fastal(tables_f, table, coded_seq,tsize );
	fclose(file);
	fclose(tables_f);
	return num_of_sequences;
}


/**
 *	Makes an index of a file
 */
int
make_pos_len_index_of_file_retree(char *file_name,			//file with sequences
								  long **file_positions,	//array to save the positions
								  int **seq_lengths)		//array to save the sequence length
{

// 	printf("HALLO\n");
	//Reading sequences
	const int LINE_LENGTH = 501;
	(*file_positions) = vcalloc(ENLARGEMENT_PER_STEP,  sizeof(long));
	(*seq_lengths) = vcalloc(ENLARGEMENT_PER_STEP,  sizeof(int));


	FILE *file = fopen(file_name,"r");

	char line[LINE_LENGTH];

	int num_of_sequences = 0;
	int mem_for_pos = ENLARGEMENT_PER_STEP;
// 	fgets(line, LINE_LENGTH , file)
	 int seq_len = 0;
			
			
	if (file == NULL)
	{
		printf("FILE NOT FOUND\n");
		exit(1);
	}
	else
	{
		int i;
		while(fgets(line, LINE_LENGTH , file)!=NULL)
		{
// 			line[LINE_LENGTH -2] = '\n';
			if (line[0] == '>')
			{
				(*file_positions)[num_of_sequences] = ftell(file);
				if (num_of_sequences >0)
				{
					(*seq_lengths)[num_of_sequences-1] = seq_len;
					seq_len = 0;
				}
				++num_of_sequences;

				if (num_of_sequences == mem_for_pos)
				{
					mem_for_pos += ENLARGEMENT_PER_STEP;
					(*file_positions) = vrealloc((*file_positions), mem_for_pos * sizeof(long));
					(*seq_lengths) = vrealloc((*seq_lengths), mem_for_pos * sizeof(int));
				}
			}
			else
			{
				i = -1;
				while ((line[++i] != '\n') && (line[i] != '\0'))
				{
					if (line[i] != '-')
					{
// 						printf("A: %c\n", line[i]);
						++seq_len;
					}
				}
			}
		}
	}
	(*seq_lengths)[num_of_sequences-1] = seq_len;
// 	printf("%i %li\n", (*seq_lengths)[0], (*file_positions)[0]);
// 	printf("%i %li\n", (*seq_lengths)[1], (*file_positions)[1]);
	fclose(file);
	return num_of_sequences;
}



int logid_score2 ( int sim, int len)
{
	float score;
  
	if ( len==0)return (int)(0.33*(float)MAXID);
  
	score=(float)sim/(float)len;
	if (score>0.9) score=1.0;
	else score=-log10 (1.0-score);
  
	score=(score*MAXID);
	return score;
}


int fast_aln2sim_mat2 (char *seq1, char *seq2)
{
	int r1, r2;
	int len = 0;
	int sim = 0;
	int c = -1;
	int simm=100;
// 	printf("SEQ: %s %s\n", seq1, seq2);
	while (seq1[++c] != '\0')
	{
		r1=tolower (seq1[c]);
		r2=tolower (seq2[c]);
		if ((r1 == '-') && (r2 == '-')) continue;
		if (r1==r2) sim++;
		len++;
	}


	simm = logid_score2 ( sim, len);
	return simm;
}



/*!
 *	\brief Function to create a tree using the PartTree algorithm.
 *	
 *	\param param_set A \a PartTree_param object containing all necessary parameters and the data.
 *	\return The node_number.
 */
void
partTree_retree(PartTree_param *param_set)
{
	int num_of_tree_nodes = param_set->num_sequences-1;
	int loop_tree_node;

	Tree_fastal *tree = param_set->tree;
// 	int this_node = param_set->pos_tree;
	
	int i;
// 	int tsize = param_set->tsize;
	
	
	//get some memory
// 	short *table1 = vcalloc(tsize, sizeof(short));
// 	short *table2 = vcalloc(tsize, sizeof(short));
	char **names = declare_char(param_set->subgroup, 8);
	int **dist_mat = declare_int(param_set->subgroup, param_set->subgroup);
	int **dist_mat2 = declare_int(param_set->subgroup, param_set->subgroup);
	char * file_name_tmp = vtmpnam(NULL);
	int *seed_set_cleaned = vcalloc(param_set->subgroup, sizeof(int));
	FILE *aln_f = param_set->ktup_table_f;
	long *file_positions = param_set->ktup_positions;
	int max_n_group = param_set->subgroup;
	int num_in_subgroup = param_set->subgroup;
	int *seq_lengths = param_set->seq_lengths;
	int *clusters = vcalloc(param_set->subgroup+1, sizeof(int));
	int *min_dist = vcalloc(param_set->num_sequences, sizeof(int));
	int *belongs_to = vcalloc(param_set->num_sequences, sizeof(int));

	int aln_length = 1;
	fseek(aln_f, file_positions[0], SEEK_SET);
	char line[500];
	while ((fgets(line, 500, aln_f) != NULL) && (line[0] != '>'))
	{
		i = -1;
		while ((line[++i] != '\n') && (line[i] != '\0'));
		aln_length += i;
	}
	char *seq1 = vcalloc(aln_length, sizeof(char));
	char *seq2 = vcalloc(aln_length, sizeof(char));
	
	//Prepare first node
	
	tree[0].index = vcalloc(param_set->num_sequences,sizeof(int));
	int *index = tree[0].index;
	for (i = 0; i< param_set->num_sequences; ++i)
		index[i] = i;
	tree[0].name = param_set->pos_tree +param_set->num_sequences;
	
	tree[0].num_leafs = param_set->num_sequences;
	int *sequence_group2 = vcalloc(param_set->num_sequences,sizeof(int));
// 	
	Tree_fastal *current_node;
	for (loop_tree_node = 0; loop_tree_node < num_of_tree_nodes; ++loop_tree_node)
	{
// 		printf("ROUND: %i\n", loop_tree_node);
		current_node = &tree[loop_tree_node];
		index= current_node->index;
		if (current_node->index == NULL)
		{
			continue;
		}
		int num_sequences = current_node->num_leafs;

		//if number of sequences in this group smaller than number subgoup size: make tree, finisch
		if (num_sequences <= max_n_group)
		{
			dist_mat = make_distance_matrix_sim(aln_f, file_positions, index, num_sequences, dist_mat);
			for (i = 0; i < num_sequences; ++i)
			{
				sprintf(names[i],"%i", current_node->index[i]);
	 		}
			NT_node **tree= (int_dist2upgma_tree_fastal (dist_mat, num_sequences, file_name_tmp , names));
			tree_process_simple(tree[0][0], param_set,loop_tree_node);
			continue;
		}


		for (i = 0; i < num_in_subgroup; ++i)
		{
			seed_set_cleaned[i] = 1;
		}
		
		//finde longest sequence and put into the first field
		
		int index_longest = 0;
		int length_of_longest = 0;
		
		for(i = 0; i < num_sequences; ++i)
		{
			if (seq_lengths[index[i]] > length_of_longest)
			{
				index_longest = i;
				length_of_longest = seq_lengths[index[i]];
			}
		}
		int tmp = index[index_longest];
		index[index_longest] = index[0];
		index[0] = tmp;

		//distance of longest to rest
		int seq_index = 1;
		read_sequence_from_position(aln_f, file_positions[index[0]], seq1);
		int min = -1;
		
		
		for (i = 1; i < num_sequences; ++i)
		{
			read_sequence_from_position(aln_f, file_positions[index[1]], seq2);
			tmp = 100 - fast_aln2sim_mat2(seq1, seq2);	
			if (tmp < min)
			{
				min = tmp;
				seq_index = i;
			}
		}

		//get the new seed_set in the first n spaces
		tmp = index[1];
		index[1] = index[seq_index];
		index[seq_index] = tmp;
		int r,j;
		num_in_subgroup = param_set->subgroup;

		
		for (i = 2; i < num_in_subgroup; ++i)
		{
			r = i + rand() / ( RAND_MAX / ( num_sequences-i) + 1 );
			tmp = index[r];
			index[r] = index[i];
			index[i] = tmp;
		}

		//Calculate matrix
		dist_mat = make_distance_matrix_sim(aln_f, file_positions, index, param_set->subgroup, dist_mat);
// 	
// 		//Filter out sequences that are to similar & reorder
// 		
		NT_node **upgma_tree;
		
		
		int num_in_clean = filter(index, dist_mat, seed_set_cleaned, param_set);

		
		if (num_in_clean ==1)
		{
			num_in_clean = 2;
			seed_set_cleaned[1] = 1;
		}
			//make_tree
		int col = 0;
		int row = 0;
		for (i = 0; i <  num_in_subgroup; ++i)
		{
			if (seed_set_cleaned[i])
			{
				row = col+1;
				for (j = i+1; j <  num_in_subgroup; ++j)
				{
					if (seed_set_cleaned[j])
					{
						dist_mat2[row][col] = dist_mat2[col][row] = dist_mat[i][j];
						++row;
					}
				}
				++col;
			}
		}
		for (i = 0; i < num_in_clean; ++i)
		{
			sprintf(names[i],"%i",i);
		}
		upgma_tree= (int_dist2upgma_tree_fastal (dist_mat2, num_in_clean, file_name_tmp , names));
 

// 		//cluster
// 		//calculate distances from n' to N
		read_sequence_from_position(aln_f, file_positions[index[0]], seq1);
// 		get_table(table1, table_f, file_positions[index[0]]);
		for (j = num_in_clean; j < num_sequences; ++j)
		{
			read_sequence_from_position(aln_f, file_positions[index[j]], seq2);
			min_dist[j] = 100 - fast_aln2sim_mat2(seq1, seq2);	
// 			min_dist[j] = euclidean_dist_half(table_f, file_positions[index[j]], table1, table2, param_set->tsize);
			belongs_to[j] = 0;
		}
		for(i = 1; i < num_in_clean; ++i)
		{
			read_sequence_from_position(aln_f, file_positions[index[0]], seq1);
// 			get_table(table1, table_f, file_positions[index[i]]);
			belongs_to[i] = i;
			for (j = num_in_clean; j < num_sequences; ++j)
			{
				read_sequence_from_position(aln_f, file_positions[index[j]], seq2);
				tmp = 100 - fast_aln2sim_mat2(seq1, seq2);
// 				tmp = euclidean_dist_half(table_f, file_positions[index[j]], table1, table2, param_set->tsize);
				if (tmp < min_dist[j])
				{
					min_dist[j] = tmp;
					belongs_to[j] = i;
				}
			}
		}
// 
		//how_many sequences has each cluster
		for (j = 0; j <= num_in_subgroup; ++j)
		{
			clusters[j] = 0;
		}
		for (j = 0; j < num_sequences; ++j)
		{
			++clusters[belongs_to[j]];
		}
// 		for (j = 0; j <= num_in_subgroup; ++j)
// 		{
// 			printf("CL: %i ",clusters[j]);
// 		}
// 		printf("\n");
		for(i = 1; i < num_in_clean; ++i)
		{
			clusters[i] += clusters[i-1];
		}
		clusters[num_in_clean] = clusters[num_in_clean-1];

		for (i = 0; i < num_sequences; ++i)
		{
			sequence_group2[--clusters[belongs_to[i]]] = index[i];
		}

		for (i = 0; i < num_sequences; ++i)
		{
			index[i] = sequence_group2[i];
		}

		
		for (i = 0; i < num_in_clean; ++i)
		{
			sprintf(names[i],"%i",i);
		}
		tree_process(upgma_tree[0][0], param_set, clusters, index, loop_tree_node);
		NT_node tmp_tree = upgma_tree[3][0];
		vfree(upgma_tree[0]);
		vfree(upgma_tree[1]);
		vfree(upgma_tree[2]);
		vfree(upgma_tree[3]);
		vfree(upgma_tree);
		free_tree(tmp_tree);
	}
// 	exit(1);
	vfree(min_dist);
	vfree(belongs_to);
	vfree(clusters);
}

