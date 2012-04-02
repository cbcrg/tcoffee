#ifndef CLASSES_H_
#define CLASSES_H_

#include <string.h>
#include <stdio.h>

// #include "km_util.h"

// Single Seq
typedef struct
{
	char *name;
	char *seq;
	size_t size;
	size_t reserved;
}Seq;

Seq*
init_Seq(char *name, int reserve);

void
append(Seq *seq, char *sequence);


// Set of sequences
typedef struct
{
	Seq **seqs;
	size_t  n_seqs;
	size_t reserved;
} SeqSet;


SeqSet*
read_fasta(char *seq_f);

void
delSeqSet(SeqSet *set);

#endif
