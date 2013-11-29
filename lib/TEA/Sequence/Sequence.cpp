/*
 * Sequence.cpp
 *
 *  Created on: Oct 9, 2011
 *      Author: Carsten Kemena
 *
 *   Copyright 2011 Carsten Kemena
 *
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


#include "Sequence.h"


namespace BioTools
{
namespace Seq
{

using namespace std;

Sequence::Sequence(const std::string &seq_name, const std::string &comment_, const std::string &seq, size_t seq_id): _name(seq_name), _comment(comment_), _sequence(seq), _input_id(seq_id), _aln_id(seq_id)
{


}

Sequence::Sequence(const std::string &seq_name, const std::string &comment_, unsigned int seq_length, size_t seq_id):_name(seq_name), _comment(comment_), _input_id(seq_id), _aln_id(seq_id)
{
	_sequence.reserve(seq_length);
}

Sequence::Sequence(const Sequence &seq):_name(seq._name), _comment(seq._comment), _sequence(seq._sequence), _input_id(seq._input_id), _aln_id(seq._aln_id)
{
}


Sequence::~Sequence()
{
	// TODO Auto-generated destructor stub
}


size_t Sequence::ungapped_size() const
{
	size_t length=0, seq_size=this->size();
	for (size_t i=0; i<seq_size; ++i)
		if (_sequence[i] != '-')
			++length;

	return length;
}



void
Sequence::insert_gaps(const vector<pair<unsigned int, unsigned int> > vec)
{
	size_t vec_size=vec.size();
	size_t seq_size = _sequence.size();
	int n_gaps = 0;
	for (size_t i=0; i<vec_size; ++i)
		n_gaps += vec[i].second;
	_sequence.resize(_sequence.size()+n_gaps);

	int seq_pos =_sequence.size();
	int j;
	for (int i = vec_size-1; i>=0; --i)
	{
		while (seq_size > vec[i].first)
			_sequence[--seq_pos]=_sequence[--seq_size];
		n_gaps=vec[i].second;
		for (j=0; j<n_gaps; ++j)
			_sequence[--seq_pos]='-';
	}
}


char
identify_seq_type(const Sequence &seq)
{
	size_t seq_len = seq.size();
	char c;
	for (unsigned int i = 0; i < seq_len; ++i)
	{
		c = tolower(seq[i]);
		if ((c != 'a') && (c != 'c') && (c != 'g') && (c != 't') && (c != 'u'))
			return 'P';
	}
	return 'N';
}


pair<size_t, size_t>
coverage(const Sequence &seq1, const Sequence &seq2)
{
	size_t len=seq1.size();
	if (len != seq2.size())
		return pair<size_t, size_t>(-1,-1);
	size_t pair_len=0, pos=0;
	for (size_t i=0; i<len; ++i)
	{
		if ((seq1[i]!='-') || (seq2[i]!='-'))
			++pair_len;
		if ((seq1[i]!='-') && (seq2[i]!='-'))
			++pos;
	}

	return pair<size_t, size_t>(pos, pair_len);
}

pair<size_t, size_t>
id(const Sequence &seq1, const Sequence &seq2)
{
	size_t len=seq1.size();
	if (len != seq2.size())
		return pair<size_t, size_t>(-1,-1);
	size_t pair_len=0, pos=0;
	for (size_t i=0; i<len; ++i)
	{
		if ((seq1[i]!='-') || (seq2[i]!='-'))
			++pair_len;
		if (seq1[i] == seq2[i])
			++pos;
	}

	return pair<size_t, size_t>(pos, pair_len);
}


bool
seq_check(const Sequence &seq1, const Sequence &seq2)
{
	size_t seq1_l = seq1.size();
	size_t seq2_l = seq2.size();
	size_t i,j=0;
	for (i=0; i < seq1_l; ++i)
	{
		if (seq1[i] == '-')
			continue;
		while ((j<seq2_l) && (seq2[j] == '-'))
			++j;

		if ((j==seq2_l) || (tolower(seq1[i]) != tolower(seq2[j])))
			return 0;
		++j;
	}
	while (j!=seq2_l)
	{
		if (seq2[j] != '-')
			return 0;
		++j;
	}
	return 1;
}



bool
bio_seq(const Sequence &seq)
{
	size_t len = seq.size();
	int c;
	for (size_t i=0; i<len; ++i)
	{
		c = tolower(seq[i]);
		if ((c!=45) && ((c<97) || (c>122)))
			return false;
	}

	return true;
}




} // end namespace Sequence
} // end namespace Biotools



