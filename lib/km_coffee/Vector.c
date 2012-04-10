// #include "Vector.h"
#include "km_coffee_header.h"
Vector *
seq2vec_kmer(const Seq *seq, short k, unsigned int *factor, size_t vec_len, size_t vec_num, short *alphabet, int *used)
{
	Vector *t = my_malloc(sizeof(Vector));
	t->data = my_calloc(vec_len, sizeof(double));
	t->id = vec_num;
	size_t i;
	unsigned int value = 0;
	for (i = 0; i<k; ++i)
		value += factor[i] * alphabet[(int)seq->seq[i]];
	++t->data[used[value]];
	size_t j=0;
	unsigned char c;
	size_t seq_len = seq->size;
	for (i = k; i<seq_len; ++i)
	{
		c = alphabet[(int)seq->seq[j]];
		value -= factor[0]*c;
		value *= factor[k-2];
		value += alphabet[(int)seq->seq[i]];
		++j;
		++t->data[used[value]];
	}

	return t;
}



int*
identify_fields(const SeqSet *seq_set, short k, unsigned int *factor, size_t *vec_len, short *alphabet )
{
	size_t vec_length = *vec_len;
	int *used = my_calloc(vec_length, sizeof(int));
	int *value_arg = my_calloc(vec_length, sizeof(int));
	int *value_test = my_malloc(vec_length * sizeof(int));


	unsigned int value;
	Seq *seq;
	size_t seq_len, j, m;
	size_t vec_num = seq_set->n_seqs;

	short l;

	seq = seq_set->seqs[9];
	value = 0;
	for (l = 0; l<k; ++l)
		value += factor[l] * alphabet[(int)seq->seq[l]];
	++value_arg[value];


	j=0;
	seq_len = seq->size;
	for (m = k; m<seq_len; ++m)
	{
		value -= factor[0]*alphabet[(int)seq->seq[j]];
		value *= factor[k-2];
		value += alphabet[(int)seq->seq[m]];
		++j;
		++value_arg[value];
	}


	size_t i;
	size_t x;
	for (i = 1; i<vec_num; ++i)
	{
		for (x = 0; x < vec_length; ++x)
		{
			value_test[x]=0;
		}
		seq = seq_set->seqs[i];
		value = 0;
		for (l = 0; l<k; ++l)
			value += factor[l] * alphabet[(int)seq->seq[l]];
		value_test[value]=1;


		j=0;
		seq_len = seq->size;
		for (m = k; m<seq_len; ++m)
		{
			value -= factor[0]*alphabet[(int)seq->seq[j]];
			value *= factor[k-2];
			value += alphabet[(int)seq->seq[m]];
			++j;
			++value_test[value];
		}

		for (x = 0; x < vec_length; ++x)
		{
			if (value_test[x] != value_arg[x])
			{
				used[x]=1;
			}
		}

	}
	free(value_test);
	free(value_arg);

	j=0;
	for (i=0; i<vec_length; ++i)
	{
		if (used[i] != 0)
			used[i] = j++;
	}
	*vec_len=j;
	printf("DIM=%li\n", *vec_len);
	return used;
}

int
my_variance_sort (const void *a, const void *b)
{
	double i = *(double*)a;
	double j = *(double*)b;
	if (i<j)
		return 1;
	if (i>j)
		return -1;
	else
		return 0;
}


int*
identify_fields_variance(const SeqSet *seq_set, short k, unsigned int *factor, size_t *vec_len, short *alphabet)
{
	size_t vec_length = *vec_len;
	int *used = my_calloc(vec_length, sizeof(int));
	double *mean = my_calloc(vec_length, sizeof(double));
	double *variance = my_calloc(vec_length, sizeof(double));
	double *value_test = my_malloc(vec_length * sizeof(double));

	size_t i;

	unsigned int value;
	Seq *seq;
	size_t seq_len, j, m;
	size_t n_vecs = seq_set->n_seqs;

// 	calc mean
	short l;
	for (i = 1; i<n_vecs; ++i)
	{
		seq = seq_set->seqs[i];
		value = 0;
		for (l = 0; l<k; ++l)
			value += factor[l] * alphabet[(int)seq->seq[l]];
		++mean[value];

		j=0;
		seq_len = seq->size;
		for (m = k; m<seq_len; ++m)
		{
			value -= factor[0]*alphabet[(int)seq->seq[j]];
			value *= factor[k-2];
			value += alphabet[(int)seq->seq[m]];
			++j;
			++mean[value];
		}
	}

	for (i = 0; i<vec_length; ++i)
		mean[i] /= n_vecs;


	//calc variance
	size_t x;
	for (i = 1; i<n_vecs; ++i)
	{
		for (x = 0; x < vec_length; ++x)
			value_test[x]=0;

		seq = seq_set->seqs[i];
		value = 0;
		for (l = 0; l<k; ++l)
			value += factor[l] * alphabet[(int)seq->seq[l]];


		j=0;
		seq_len = seq->size;
		for (m = k; m<seq_len; ++m)
		{
			value -= factor[0]*alphabet[(int)seq->seq[j]];
			value *= factor[k-2];
			value += alphabet[(int)seq->seq[m]];
			++j;
			++value_test[value];
		}

		for (j = 0; j<vec_length; ++j)
			variance[j] += (value_test[j]-mean[j])*(value_test[j]-mean[j]);
	}
	free(value_test);
	free(mean);

	for (i = 0; i<vec_length; ++i)
		variance[i] /= n_vecs;

	qsort (variance, vec_length, sizeof(double), my_variance_sort );
// 	for (i = 0; i<vec_length; ++i)
// 		printf("%f ", variance[i]);
// 	printf("\n");
	double threshold=variance[399];

	j=0;
	for (i=0; i<vec_length; ++i)
	{
		if (variance[i] >= threshold)
			used[i] = j++;
	}
	free(variance);
	*vec_len=j;
	printf("DIM=%li\n", *vec_len);
	return used;
}

VectorSet*
seqset2vecs_kmer(SeqSet *seq_set, short k, short alphabet_size, short *alphabet)
{
	VectorSet *vec_set = my_malloc(sizeof(VectorSet));
	vec_set->dim = seq_set->seqs[0]->size;
	size_t n_seqs = seq_set->n_seqs;
	unsigned int *factor = my_malloc(k*sizeof(unsigned int));
	factor[k-1] = 1;

	short j;
	for (j=k-2; j>=0; --j)
		factor[j] = factor[j+1] *alphabet_size;
	size_t vec_len = factor[0] *alphabet_size;

	Vector **vecs= malloc(n_seqs*sizeof(Vector*));
	vec_set->n_vecs=n_seqs;
	vec_set->vecs=vecs;
	int *used = identify_fields_variance(seq_set, k, factor, &vec_len, alphabet);
	size_t i;
	for (i = 0; i<n_seqs; ++i)
		vecs[i] = seq2vec_kmer(seq_set->seqs[i], k, factor, vec_len, i, alphabet, used);
	vec_set->dim=vec_len;
	free(used);
	free(factor);
	return vec_set;
}

void
delVecSet(VectorSet* set)
{
	size_t n_vecs = set->n_vecs;
	int i;
	Vector** vecs = set->vecs;
	for (i=0; i < n_vecs; ++i)
	{
		free(vecs[i]->data);
		free(vecs[i]);
	}
	free(vecs);
	free(set);
}



