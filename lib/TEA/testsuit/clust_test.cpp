/*
 * clust_test.cpp
 *
 *  Created on: Jan 21, 2012
 *      Author: Carsten Kemena
 *
 * This file is part of BioTools.
 *
 * BioTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BioTools is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BioTools.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <cstdlib>
#include <cstdio>
#include <vector>
#include <ctime>
#include <cstring>


#include "../Clustering/clustering.h"
#include "../Clustering/Vector.h"
#include "../Clustering/kmeans.h"
#include "../Sequence/SequenceSet.h"

using namespace std;
using namespace BioTools::Seq;
using namespace BioTools::Clustering;

typedef boost::shared_ptr<Sequence> Seq_ptr;

void
read_fasta_file(vector<Seq_ptr> &seq_set, const string &seq_f)
{
	const int LINE_LENGTH = 501;
	char line[LINE_LENGTH];
	FILE *seq_F = fopen(seq_f.c_str(), "r");
	if (seq_F == NULL)
		printf("STRANGE %s \n", seq_f.c_str());
	Sequence *tmp = NULL;
	while (fgets (line , LINE_LENGTH , seq_F) != NULL )
	{
		line[strlen(line)-1]='\0';
		if (line[0] == '>')
		{
			tmp = new Sequence(string(&line[1]), "", "");
			seq_set.push_back(Seq_ptr(tmp));
		}
		else
		{
			tmp->append(line);
		}
	}
	fclose(seq_F);



}

void
write_files(KM_node *node, int *vecs, const SequenceSet &seq_set, char *name_f)
{
	size_t end = node->end;
	FILE * tmp_F = fopen(name_f, "w");
	for (size_t i = node->start; i<end; ++i)
	{
		fprintf(tmp_F, ">%s\n%s\n", seq_set[vecs[i]].name().c_str(), seq_set[vecs[i]].sequence().c_str() );
	}
	fclose(tmp_F);

}


void
del_tree(KM_node* root)
{
	stack<pair<KM_node*, size_t> > to_do;
		to_do.push(pair<KM_node*, size_t>(root, 0));
		KM_node* current;
		size_t child;
		while (!to_do.empty())
		{
			current = to_do.top().first;
			child = to_do.top().second;
			if (current->childs.empty())
			{
				to_do.pop();
				delete current;
			}
			else
			{
				if (child == current->childs.size())
				{
					delete current;
					to_do.pop();

				}
				else
				{
					++to_do.top().second;
					to_do.push(pair<KM_node*, size_t>(current->childs[child],0));
				}
			}
		}
}

void
traverse_km_tree(KM_node* root, int *vecs, const SequenceSet &seq_set)
{
	stack<pair<KM_node*, size_t> > to_do;
	to_do.push(pair<KM_node*, size_t>(root, 0));
	KM_node* current;
	size_t child;
	size_t j;
	size_t pos;
	char command[1000];
	while (!to_do.empty())
	{
		current = to_do.top().first;
		child = to_do.top().second;
		if (current->childs.empty())
		{
			if(current->end-current->start > 1)
			{
				sprintf(command, "%li",current->id);
				write_files(current, vecs, seq_set, command);
				sprintf(command,"t_coffee -in %li -output fasta_aln -outfile %li.fa -quiet >/dev/null 2>/dev/null", current->id, current->id);
				system(command);
			}
			else
			{
				sprintf(command, "%li.fa",current->id);
				write_files(current, vecs, seq_set, command);
			}
			to_do.pop();
		}
		else
		{
			if (child == current->childs.size())
			{
				sprintf(command, "t_coffee -output fasta_aln -quiet -outfile %li.fa -profile ", current->id);
				for(j=0; j<current->childs.size();++j)
				{
					pos=strlen(command);
					sprintf(&command[pos], "%li.fa ", current->childs[j]->id);
				}
				pos=strlen(command);
				sprintf(&command[pos], "2>/dev/null >/dev/null");
				system(command);
				to_do.pop();

			}
			else
			{
				++to_do.top().second;
				to_do.push(pair<KM_node*, size_t>(current->childs[child],0));
			}
		}
	}
}


int
main(int argc, char *argv[])
{
	string seq_f = "";
	if (argc == 2)
		seq_f.append(argv[1]);
	srand(time(NULL));
	short alphabet[127];
	short j = -1;
	for (short i = 65; i < 91; ++i)
		if ((i==66) || (i==74) || (i==79) || (i==88) || (i==90))
			alphabet[i] = 0;
		else
			alphabet[i] = ++j;
	j=-1;
	for (short i = 97; i < 123; ++i)
		if ((i==98) || (i==106) || (i==111) || (i==120) || (i==122))
			alphabet[i] = 0;
		else
			alphabet[i] = ++j;

	SequenceSet seq_set;
	seq_set.read_fasta(seq_f);


	vector<Vec_double_ptr>* vecs = seqset2vecs_kmer(seq_set, 3, 21, alphabet);
	printf("TURNED SEQUENCES\n");
	string s = "random";
	KM_node *root = hierarchical_kmeans(*vecs, 100, "first", 0.001);
	printf("Traverse\n");
	exit(1);
	size_t n_vecs = seq_set.num_seqs();
	int *assignment = new int[n_vecs];
	for (size_t i = 0; i< n_vecs; ++i)
		assignment[i]=(*vecs)[i]->id();
	(*vecs).clear();
	delete vecs;
	traverse_km_tree(root, assignment, seq_set);
	delete[] assignment;
	del_tree(root);
//	printf("ERG: %li %li %li %li %li %li", ((*vecs)[0])->assignment(), ((*vecs)[1])->assignment(), ((*vecs)[2])->assignment(), ((*vecs)[3])->assignment(), ((*vecs)[4])->assignment(), ((*vecs)[5])->assignment() );

	return EXIT_SUCCESS;
}


