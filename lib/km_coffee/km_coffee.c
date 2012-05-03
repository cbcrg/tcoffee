
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

void
traverse_km_tree(KM_node* root, int *vecs, const SeqSet *seq_set, char *out_f, int n_cores, char *method)
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
				sprintf(command,"clustalo --guidetree-out co_tree.dnd -i %li --force >/dev/null 2>/dev/null", current->id);
				if (system(command))
				{
					printf("%s\n",command);
					exit(1);
				}
				if (current->id!=0)
					sprintf(command,"t_coffee -usetree co_tree.dnd -in %li -output fasta_aln -dp_mode gotoh_pair_wise -outfile %li.fa -n_core %i -method %s -quiet >/dev/null 2>/dev/null", current->id, current->id, n_cores, method);
				else
					sprintf(command,"t_coffee -usetree co_tree.dnd -in %li -output fasta_aln  -dp_mode gotoh_pair_wise -outfile %s -n_core %i -method %s -quiet >/dev/null 2>/dev/null ",  current->id, out_f, n_cores, method);
				if (system(command))
				{
					printf("%s\n",command);
					exit(1);
				}
			}
			else if(current->end-current->start > 1)
			{
				sprintf(command, "%li",current->id);
				write_files(current, vecs, seq_set, command);
				if (current->id!=0)
					sprintf(command,"t_coffee -in %li -output fasta_aln -dp_mode gotoh_pair_wise -outfile %li.fa -n_core %i -method %s -quiet >/dev/null 2>/dev/null", current->id, current->id, n_cores, method);
				else
					sprintf(command,"t_coffee -in %li -output fasta_aln  -dp_mode gotoh_pair_wise -outfile %s -n_core %i -method %s -quiet >/dev/null 2>/dev/null ",  current->id, out_f, n_cores, method);
				if (system(command))
				{
					printf("%s\n",command);
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
					sprintf(command, "t_coffee -output fasta_aln  -dp_mode gotoh_pair_wise -quiet -outfile %li.fa -n_core %i -method %s -profile FILE::prf.fa 2>/dev/null >/dev/null", current->id,n_cores, method);
				else
					sprintf(command, "t_coffee -output fasta_aln  -dp_mode gotoh_pair_wise -quiet -outfile %s -n_core %i -method %s -profile FILE::prf.fa 2>/dev/null >/dev/null", out_f, n_cores, method);
				FILE *prf_F = my_fopen("prf.fa", "w");
				for(j=0; j<current->n_children;++j)
				{

					fprintf(prf_F, "%li.fa\n", current->children[j]->id);
				}
				fclose(prf_F);

				if (system(command))
					printf("%s\n",command);
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


int
km_coffee_align3(char *seq_f, int k, char *method, char *aln_f, int n_cores, char *init)
{
	char *use_as_temp = get_tmp_4_tcoffee();

	#ifdef _OPENMP
		omp_set_num_threads(n_cores);
	#endif


	SeqSet *seq_set = read_fasta(seq_f);
	qsort(seq_set->seqs, seq_set->n_seqs, sizeof(Seq*), my_seq_sort);
// 	printf("%s %s %i %i %i\n", seq_f, aln_f, k, n_cores,seq_set->n_seqs);
	srand(time(0));
	short alphabet[127];
	short j = -1;
	short i;
// 	printf("Convert Seqs\n");
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



	VectorSet *vec_set = seqset2vecs_kmer(seq_set, 2, 21, alphabet);
	KM_node *root = hierarchical_kmeans(vec_set, k, init, 0.001);

// 	printf("clustered\n");
	char template[400];
	sprintf(template, "%s/km_coffee_tmp_XXXXXX", use_as_temp);
	char tmp_str[FILENAME_MAX];
	km_cwd = getcwd(tmp_str, FILENAME_MAX);

	km_tmp_dir = my_make_temp_dir(template, "main", "main.c");
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
	traverse_km_tree(root, assignment, seq_set, out_f, n_cores, method);
	free( assignment);
	del_tree(root);
	delSeqSet(seq_set);
	char command[100];

// 	sprintf(command, "rm -r %s", km_tmp_dir);
// 	system(command);
	free(km_tmp_dir);



	return EXIT_SUCCESS;
}
