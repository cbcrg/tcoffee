
// #include "km_coffee.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <fcntl.h>
#include <sys/file.h>
#include <sys/types.h>
#include <unistd.h>

#include "io_lib_header.h"
#include "util_lib_header.h"
#include "define_header.h"

#include "km_coffee_header.h"

/**
 * \brief Writes the sequences to the file
 * \param node The node containing the data to write
 * \param vecs The vector set
 * \param seq_set The sequence Set
 * \param name_f The file to write the sequence into.
 */
void write_files(KM_node *node, int *vecs, const SeqSet *seq_set, char *name_f)
{
	size_t end = node->end;
	FILE * tmp_F = fopen(name_f, "w");
	size_t i;
	for (i = node->start; i<end; ++i)
	{
		fprintf(tmp_F, ">%s\n%s\n", seq_set->seqs[vecs[i]]->name, seq_set->seqs[vecs[i]]->seq);
	}
	fclose(tmp_F);
}

/**
 * \brief Frees all memory occupied by this structure
 * \param root The root of the tree to delete.
 */
void
del_tree(KM_node* root)
{
		Stack *to_do = Stack_init();

		Node_pair *tmp = my_malloc(sizeof(Node_pair));
		tmp->node=root;
		tmp->id = 0;
		push(to_do, tmp);

		KM_node* current;
		size_t child;
		while (to_do->size != 0)
		{
			current = ((Node_pair*)to_do->last)->node;
			child = ((Node_pair*)to_do->last)->id;
			if (current->n_children==0)
			{
				free(pop(to_do));
				delKM_node(current);
			}
			else
			{
				if (child == current->n_children)
				{
					delKM_node(current);
					free(pop(to_do));
				}
				else
				{
					++((Node_pair*)to_do->last)->id;
					tmp = my_malloc(sizeof(Node_pair));
					tmp->node=current->children[child];
					tmp->id = 0;
					push(to_do, tmp);
				}
			}
		}
		delStack(to_do);
}


/**
 * \brief Traverses the tree and calls t_coffee to produce the alignments
 * \param root The root of the tree.
 * \param vecs The vector set.
 * \param seq_set The sequence set.
 * \param out_f The output file.
 * \param n_core The number of cores to use (OPEN-MP for kmeans/fork in T-Coffee)
 * \param gapopen The gapopening costs
 * \param gapext The gapextension costs
 * \param method The method to use in the alignment
 */
void
traverse_km_tree(KM_node* root, int *vecs, const SeqSet *seq_set, char *out_f, int n_cores, int gapopen, int gapext, char *method)
{
	Stack *to_do =Stack_init();
	Node_pair *tmp = my_malloc(sizeof(Node_pair));
	tmp->node=root;
	tmp->id = 0;
	push(to_do, tmp);

	KM_node* current;
	size_t child;
	size_t j;
// 	size_t pos;
	char command[1000];
	while (to_do->size != 0)
	{

		current = ((Node_pair*)to_do->last)->node;
		child = ((Node_pair*)to_do->last)->id;
		if (current->n_children==0)
		{
			if(current->end-current->start > 2)
			{
				sprintf(command, "%li",current->id);
				write_files(current, vecs, seq_set, command);
// 				sprintf(command,"clustalo --guidetree-out co_tree.dnd -i %li --force >/dev/null 2>/dev/null", current->id);
// 				if (system(command))
// 				{
// 					printf("%s\n",command);
// 					exit(1);
// 				}
				if (current->id!=0)
					sprintf(command,"t_coffee -in %li -output fasta_aln -outfile %li.fa -n_core %i -gapopen %i -gapext %i -method %s -quiet >/dev/null 2>/dev/null", current->id, current->id, n_cores, gapopen, gapext, method);
				else
					sprintf(command,"t_coffee -in %li -output fasta_aln -outfile %s -n_core %i -gapopen %i -gapext %i -method %s -quiet >/dev/null 2>/dev/null ",  current->id, out_f, n_cores, gapopen, gapext, method);
				if (system(command))
				{
					printf("ERROR when running: %s\n",command);
					exit(1);
				}
			}
			else
				if(current->end-current->start > 1)
			{
				sprintf(command, "%li",current->id);
				write_files(current, vecs, seq_set, command);
				if (current->id!=0)
					sprintf(command,"t_coffee -in %li -output fasta_aln -outfile %li.fa -n_core %i -gapopen %i -gapext %i -method %s -quiet >/dev/null 2>/dev/null", current->id, current->id, n_cores, gapopen, gapext, method);
				else
					sprintf(command,"t_coffee -in %li -output fasta_aln -outfile %s -n_core %i -gapopen %i -gapext %i -method %s -quiet >/dev/null 2>/dev/null ",  current->id, out_f, n_cores, gapopen, gapext, method);

				printf("%s\n", command);
				if (system(command))
				{
					printf("ERROR when running: %s\n",command);
					exit(1);
				}
			}
			else
			{
				sprintf(command, "%li.fa",current->id);
				write_files(current, vecs, seq_set, command);
			}
			tmp = pop(to_do);
		}
		else
		{
			if (child == current->n_children)
			{
				if (current->id!=0)
					sprintf(command, "t_coffee -output fasta_aln -method %s -quiet -outfile %li.fa -n_core %i -gapopen %i -profile FILE::prf.fa ", method, current->id,n_cores, gapopen);
				else
					sprintf(command, "t_coffee -output fasta_aln -method %s -quiet -outfile %s -n_core %i -gapopen %i   -profile FILE::prf.fa ", method, out_f, n_cores, gapopen);
				FILE *prf_F = my_fopen("prf.fa", "w");
				for(j=0; j<current->n_children;++j)
				{
					fprintf(prf_F, "%li.fa\n", current->children[j]->id);
				}
				fclose(prf_F);

				printf("%s\n", command);
				if (system(command))
				{

					printf("ERROR when running: %s\n",command);
					exit(1);
				}
				pop(to_do);
			}
			else
			{
				++((Node_pair*)to_do->last)->id;
				tmp = my_malloc(sizeof(Node_pair));
				tmp->node=current->children[child];
				tmp->id = 0;
				push(to_do, tmp);
			}
		}
	}
	exit(1);
}


int
my_seq_sort (const void *i, const void *j)
{
	return strcmp((*(Seq**)i)->seq,(*(Seq**)j)->seq);
}




void
print_km_tree(KM_node *root, int *vecs, const SeqSet *seq_set, char *out_f)
{
	FILE *out_F = fopen(out_f, "w");
	Stack *to_do =Stack_init();
	Node_pair *tmp = my_malloc(sizeof(Node_pair));
	tmp->node=root;
	tmp->id = 0;
	push(to_do, tmp);

	KM_node* current;
	size_t child;
	size_t k;
	// 	size_t pos;
// 	char command[1000];
	while (to_do->size != 0)
	{

		current = ((Node_pair*)to_do->last)->node;
		child = ((Node_pair*)to_do->last)->id;
		if (current->n_children==0)
		{
			if (current->end-current->start >1)
				fprintf(out_F, "(");
			for (k=current->start; k<current->end; ++k)
			{
				fprintf(out_F, "%s", seq_set->seqs[vecs[k]]->name);
				if (k<current->end-1)
					fprintf(out_F, ",");
			}
			tmp = pop(to_do);
			if (current->end-current->start >1)
				fprintf(out_F, ")");
		}
		else
		{
			if (child == 0)
			{
				fprintf(out_F, "(");
			}
			if (child == current->n_children)
			{
				fprintf(out_F, ")");
				pop(to_do);
			}
			else
			{
				if (child != 0)
					fprintf(out_F, ",");
				++((Node_pair*)to_do->last)->id;
				tmp = my_malloc(sizeof(Node_pair));
				tmp->node=current->children[child];
				tmp->id = 0;
				push(to_do, tmp);
			}
		}
	}
	fprintf(out_F, ";");
}


typedef struct
{
	size_t id;
	float sim;
} sorter;




int
my_val_compare (const void *i, const void *j)
{
	return (*(Vector**)j)->val < (*(Vector**)j)->val;
}

/*
KM_node*
simple_clust(VectorSet *vecSet, unsigned int k)
{
	size_t n_vecs = vecSet->n_vecs;
	size_t i, max_id=0, max_val=0;
	Vector **vecs = vecSet->vecs;
	size_t dim=vecSet->dim;

	// 	determine longest_seq
	for (i=0; i<n_vecs; ++i)
	{
		if (vecs[i]->seq_len > max_val)
		{
			max_val=vecs[i]->seq_len;
			max_id=i;
		}
	}

	size_t num_nodes= ceil(n_vecs*1.0/k);
	KM_node **nodes = (KM_node**)malloc(sizeof(KM_node*)*num_nodes);
	size_t j;

	//determine leaves
	size_t node_id=0;
	for (j=0; j<num_nodes; ++j)
	{
		for (i=0; i<n_vecs; ++i)
			vecs[i]->val=km_sq_dist(vecs[max_id], vecs[i],dim);
		qsort (vecs, n_vecs, sizeof(Vector*), my_val_compare);
		nodes[node_id] = malloc(sizeof(KM_node));
		nodes[node_id]->children = NULL;
		nodes[node_id]->start = node_id*k;
		nodes[node_id]->end = (n_vecs>k )? (node_id+1)*k: (node_id*k)+n_vecs;
		nodes[node_id]->n_children=0;
		nodes[node_id]->id=++node_id;
		n_vecs -= k;
		max_id =0;
		vecs+=k;
	}

	// generate inner nodes
	j=num_nodes;
	KM_node *tmp;
	size_t pos=0;
	size_t overall_pos=0;
// 	printf("1: %li %li\n", num_nodes, node_id);
	while (1)
	{
		for(i=0; i<num_nodes; ++i)
		{
// 			printf("RUN %li\n", num_nodes);
			if ((i%k)==0)
			{
				tmp=malloc(sizeof(KM_node));
				tmp->children=malloc(sizeof(KM_node*)*num_nodes);
				tmp->children[0]=nodes[i];
				tmp->n_children=1;
				tmp->id=++node_id;
				nodes[overall_pos]=tmp;
				pos=0;
				++overall_pos;
			}
			else
			{
				tmp->children[++pos]=nodes[i];
				++tmp->n_children;
			}
		}
		overall_pos=0;
		if (num_nodes==1)
			break;
		num_nodes=ceil(num_nodes*1.0/k);

	}
// 	printf("%li\n", node_id);
	free(nodes);
	tmp->id=0;
	return tmp;
}*/

/*
typedef struct KM_node{
	size_t start;
	size_t end;
	struct KM_node **children;
	size_t id;
	size_t n_children;
	} KM_node;*/



int
km_coffee_align3(char *seq_f, int k, int k_leaf, char *method, char *aln_f, int n_cores, int gapopen, int gapext, char *init)
{
	char *use_as_temp = get_tmp_4_tcoffee();

	#ifdef _OPENMP
		omp_set_num_threads(n_cores);
	#endif

	SeqSet *seq_set = read_fasta(seq_f);
	qsort(seq_set->seqs, seq_set->n_seqs, sizeof(Seq*), my_seq_sort);
	srand(time(0));


	short j = -1;
	short i;
	/****************************************************
	Sequences to vector using k-mers
	*****************************************************/
	short alphabet[256];

	// standard alphabet
	for (i = 65; i < 91; ++i)
		if ((i==66) || (i==74) || (i==79) || (i==88) || (i==90))
			alphabet[i] = 0;
		else
			alphabet[i] = ++j;
	j=-1;
	for (i = 97; i < 123; ++i)
		if ((i==98) || (i==106) || (i==111) || (i==120) || (i==122))
			alphabet[i] = 0;
		else
			alphabet[i] = ++j;

	// shrinked alphabet
//	for (i = 0; i < 256; ++i)
//		alphabet[i] = 0;

// 	char *groups[]={"LlVvIiMmCcAaGgSsTtPpFfYyWw","EeDdNnQqKkRrHh"};
//	char *groups[]={"LlVvIiMmCc","AaGgSsTtPp","FfYyWw","EeDdNnQqKkRrHh"};
//	size_t n_groups = 4;
// 	size_t len;
// 	char *group;
// 	for (i=0; i<n_groups; ++i)
// 	{
// 		group=groups[i];
// 		len=strlen(group);
// 		for (j=0; j<len; ++j)
// 			alphabet[group[j]]=i+1;
// 	}

	VectorSet *vec_set = seqset2vecs_kmer(seq_set, 2, 21, alphabet);


	/****************************************************
		Sequences to vector using distances
	*****************************************************/
	// 	char *groups[]={"LVIMC","AGSTP","FYW","EDNQKRH"};
// 	size_t n_groups = 4;
// 	char *groups[]={"LlVvIiMmCc","AaGgSsTtPp","FfYyWw","EeDdNnQqKkRrHh"};
// 	char *groups[]={"LlVvIiMmCcAaGgSsTtPpFfYyWw","EeDdNnQqKkRrHh"};
// 	size_t n_groups = 2;
// 	VectorSet *vec_set = seqset2vecs_whatever(seq_set, groups, n_groups);






// 	char vec_file[500];
// 	sprintf(vec_file, "%s_2_8_%li_%li.txt", strrchr(seq_f, '/')+1, vec_set->n_vecs, vec_set->dim);

// 	print_vecs(vec_set, &vec_file[0]);
// 	read_vecs(vec_set, "matrix_59");
// 	exit(1);
//	normalize(vec_set);
	KM_node *root = hierarchical_kmeans(vec_set, k, k_leaf, init, 0.001);
// 	KM_node *root = simple_clust(vec_set, k);


	char template[400];
	sprintf(template, "%s/km_coffee_tmp_XXXXXX", use_as_temp);
	char tmp_str[FILENAME_MAX];
	km_cwd = getcwd(tmp_str, FILENAME_MAX);

	km_tmp_dir = my_make_temp_dir(template, "main");
	chdir(km_tmp_dir);
	char out_f[500];
	if (aln_f[0] != '/')
		sprintf(out_f, "%s/%s", km_cwd, aln_f);
	else
		sprintf(out_f, "%s", aln_f);



	size_t n_vecs = seq_set->n_seqs;
	int *assignment = malloc(n_vecs*sizeof(int));
	size_t l;
	for (l = 0; l< n_vecs; ++l)
		assignment[l]=vec_set->vecs[l]->id;

// 	printf("TRAVERSE\n");
	delVecSet(vec_set);
	traverse_km_tree(root, assignment, seq_set, out_f, n_cores, gapopen, gapext, method);
	free( assignment);
	del_tree(root);
	delSeqSet(seq_set);

	free(km_tmp_dir);




	return EXIT_SUCCESS;
}
