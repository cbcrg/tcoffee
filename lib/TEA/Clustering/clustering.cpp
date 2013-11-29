/*
 * clustering.cpp
 *
 *  Created on: Dec 30, 2011
 *      Author: Carsten Kemena
 *
 * This file is part of BioTools++.
 *
 * BioTools++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BioTools++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BioTools++.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "clustering.h"



namespace BioTools
{
namespace Clustering
{

typedef boost::shared_ptr<Vector<double> > Vec_double_ptr;
typedef boost::shared_ptr<Vector<unsigned int> > Vec_uint_ptr;
using namespace std;
using namespace BioTools::Seq;
using namespace BioTools::Utils;


Vec_double_ptr
seq2vec_kmer(const Sequence &seq, short k, unsigned int *factor, size_t vec_len, size_t vec_num, short *alphabet, int *used)
{
	Vector<double> *t = new Vector<double>(vec_len, 0, vec_num);
	Vec_double_ptr vec(t);
	unsigned int value = 0;
	short bad=0;
	short bad_val=alphabet[0];
	for (short i = 0; i<k; ++i)
	{
		value += factor[i] * alphabet[static_cast<int>(seq[i])];
		if (alphabet[static_cast<int>(seq[i])] == bad_val)
			++bad;
	}
	if (!bad)
		++(*t)[used[value]];
	size_t j=0;
	unsigned char c;
	size_t seq_len = seq.size();

	for (size_t i = k; i<seq_len; ++i)
	{
		c = alphabet[static_cast<int>(seq[j])];
		if (c == bad_val)
			--bad;
		value -= factor[0]*c;
		value *= factor[k-2];
		value += alphabet[static_cast<int>(seq[i])];
		if (alphabet[static_cast<int>(seq[i])] == bad_val)
			++bad;
		++j;
		if (!bad)
			++(*t)[used[value]];
	}
	return vec;
}



int*
identify_fields(const SequenceSet &seq_set, short k, unsigned int *factor, size_t &vec_len, short *alphabet )
{
	int *used = new int[vec_len];
	size_t i;
	for (i = 0; i<vec_len; used[i++]=0);

	unsigned int value;
	size_t seq_len, j, m;
	size_t vec_num = seq_set.n_seqs();
	short l;
	short bad=0;
	short bad_val=alphabet[0];
	unsigned char c;
	for (i = 0; i<vec_num; ++i)
	{
		bad=0;
		const Sequence &seq = seq_set[i];
		value = 0;
		for (l = 0; l<k; ++l)
		{
			value += factor[l] * alphabet[static_cast<int>(seq[l])];
			if (alphabet[static_cast<int>(seq[l])] == bad_val)
				++bad;
		}
		if (!bad)
			used[value]=true;

		j=0;
		seq_len = seq.size();

		for (m = k; m<seq_len; ++m)
		{
			c = alphabet[static_cast<int>(seq[j])];
			if (c == bad_val)
				--bad;

			value -= factor[0]*c;
			value *= factor[k-2];
			value += alphabet[static_cast<int>(seq[m])];
			if (alphabet[static_cast<int>(seq[m])] == bad_val)
				++bad;
			++j;
			if (!bad)
				used[value]=true;
		}
		/*
		for (m = k; m<seq_len; ++m)
		{
			value -= factor[0]*alphabet[static_cast<int>(seq[j])];
			value *= factor[k-2];
			value += alphabet[static_cast<int>(seq[m])];
			++j;
			++used[value];

		}*/
	}

	j=0;
	for (i=0; i<vec_len; ++i)
	{
	//	printf("%li %i\n", i,used[i]);
		if (used[i])
			used[i] = j++;

	}
	vec_len=j;
	return used;
}



vector<Vec_double_ptr>*
seqset2vecs_kmer(const SequenceSet &seq_set, short k, const string &alphabet)
{
	short *encoded = encode(alphabet);
	short alphabet_size =encoded[0];
	//printf("%i", alphabet_size);
	size_t n_seqs = seq_set.n_seqs();
	unsigned int *factor = new unsigned int [k];
	factor[k-1] = 1;
	for (int i=k-2; i>=0; --i)
		factor[i] = factor[i+1] *alphabet_size;
	size_t vec_len = factor[0] *alphabet_size;
	//printf("%li\n", vec_len);
	vector<Vec_double_ptr> *vec_set= new vector<Vec_double_ptr>;
	int *used = identify_fields(seq_set, k, factor, vec_len, encoded);
	for (size_t i = 0; i<n_seqs; ++i)
		vec_set->push_back(seq2vec_kmer(seq_set[i], k, factor, vec_len, i, encoded, used));
	delete[] used;
	delete[] factor;
	return vec_set;
}





}
}
