// #include "classes.h"
#include "km_coffee_header.h"

// Seq

Seq*
init_Seq(char *name, int reserve)
{
	Seq *tmp_S = my_malloc(sizeof(Seq));
	size_t str_len=strlen(name);
	tmp_S->name=my_malloc((str_len+1)*sizeof(char));
	strncpy(tmp_S->name, name, str_len);
	tmp_S->name[str_len]='\0';
	tmp_S->reserved=reserve;
	tmp_S->seq=my_malloc(reserve*sizeof(char));
	tmp_S->size=0;
	return tmp_S;
}

void
append(Seq *seq, char *sequence)
{
	size_t len = strlen(sequence);
	if (seq->size+len+1>=seq->reserved)
	{
		seq->reserved +=len+1;
		seq->seq=my_realloc(seq->seq, seq->reserved*sizeof(char));
	}
	strncpy(&seq->seq[seq->size], sequence, len);
	seq->size+=len;
	seq->seq[seq->size]='\0';
}


// Seq Set

void
add(Seq *seq, SeqSet *set)
{

	if (set->n_seqs == set->reserved)
	{
		set->reserved += 10;
		set->seqs = my_realloc(set->seqs, set->reserved*sizeof(Seq **));
	}
	set->seqs[set->n_seqs] = seq;
	++set->n_seqs;
}


SeqSet*
read_fasta(char *seq_f)
{

	FILE *seq_F = my_fopen(seq_f, "r");
	SeqSet *set = my_malloc(sizeof(SeqSet));
	set->seqs=NULL;
	set->n_seqs=0;
	set->reserved=0;

	const unsigned int LINE_LENGTH = 500;
	char line[LINE_LENGTH];

	Seq *tmp_S = NULL;
	while (fgets(line, LINE_LENGTH, seq_F) != NULL)
	{
		if (line[0] == '>')
		{
			line[strlen(line)-1] = '\0';
			tmp_S=init_Seq(&line[1], 100);
			add(tmp_S, set);
		}
		else
		{
			if (line[strlen(line)-1] == '\n')
				line[strlen(line)-1] = '\0';
			append(tmp_S, line);
		}
	}
	fclose(seq_F);
	return set;
}


void
delSeqSet(SeqSet *set)
{
	size_t n_seqs = set->n_seqs;
	size_t i;
	for (i=0; i<n_seqs; ++i)
	{
		free(set->seqs[i]->name);
		free(set->seqs[i]->seq);
		free(set->seqs[i]);
	}
	free(set->seqs);
	free(set);
}



