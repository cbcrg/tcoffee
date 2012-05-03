#ifndef VECTOR_H_
#define VECTOR_H_

#include <stdlib.h>
#include <math.h>
// #include "classes.h"
// #include "km_util.h"

typedef struct
{
	double *data;
	size_t seq_len;
	size_t id;
	size_t assignment;
} Vector;


typedef struct
{
	Vector **vecs;
	size_t n_vecs;
	size_t used;
	size_t dim;
} VectorSet;

Vector *
seq2vec_kmer(const Seq *seq, short k, unsigned int *factor, size_t vec_len, size_t vec_num, short *alphabet, int *used);

VectorSet*
seqset2vecs_kmer(SeqSet *seq_set, short k, short alphabet_size, short *alphabet);

void
delVecSet(VectorSet *set);

#endif
