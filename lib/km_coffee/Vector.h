#ifndef VECTOR_H_
#define VECTOR_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <limits.h>
// #include "classes.h"
// #include "km_util.h"


typedef struct
{
	size_t x;
	double y;
} id_pair;

typedef struct
{
	double *data;
	size_t seq_len;
	size_t id;
	size_t assignment;
	double val;
} Vector;


typedef struct
{
	Vector **vecs;
	size_t n_vecs;
	size_t used;
	size_t dim;
} VectorSet;


/**
 * \brief Turns a sequence into a vector
 * \param seq The sequence to use
 * \param k The kmer size
 * \param factor The multiplication factor for the coding
 * \param vec_len The length of the new vector
 * \param vec_num The vector id to set
 * \param alphabet The alphabet coding to use
 * \param used shifting of the entries
 * \return A pointer to to the vector encoding the given sequence
 */
Vector *
seq2vec_kmer(const Seq *seq, short k, unsigned int *factor, size_t vec_len, size_t vec_num, short *alphabet, int *used);


VectorSet*
seqset2vecs_dist(SeqSet *seq_set, char *groups[], size_t n_groups);

VectorSet*
seqset2vecs_whatever(SeqSet *seq_set, char *groups[], size_t n_groups);

/**
 * \brief Turns as sequence set into a set of vectors
 * \param seq_set The sequence set
 * \param k The size of the kmers
 * \param alphabet_size The size of the alphbet to use
 * \param alphabet The encoding of the alphabet
 * \return A pointer to the created vector set
 */
VectorSet*
seqset2vecs_kmer(SeqSet *seq_set, short k, short alphabet_size, short *alphabet);

/**
 * \brief Calculates the squared distance between two vectors.
 * \param vec1 The first vector.
 * \param vec2 The second vector.
 * \param dim The length of the two vectors.
 */
double
km_sq_dist(const Vector *vec1, const Vector *vec2, size_t dim);

/**
 * \brief Calculated the angle between two vectors
 * \param vec1 The first vector.
 * \param vec2 The second vector.
 * \param dim The length of the two vectors.
 */
double
km_angle_dist(const Vector *vec1, const Vector *vec2, size_t dim);

double
km_muscle_dist(const Vector *vec1, const Vector *vec2, size_t dim, int k);

void
delVecSet(VectorSet *set);

Vector *
new_vec(Vector *vec, int vec_len);

Vector *
new_vec_nodata(Vector *vec, int vec_len);

void
print_vecs(VectorSet *set, char *out_f);

void
read_vecs(VectorSet *set, char *in_f);

SeqSet *
find_distant(SeqSet *seq_set, VectorSet *set);

#endif
