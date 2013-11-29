/*
 * SequenceSet.cpp
 *
 *  Created on: Jan 28, 2012
 *      Author: Carsten Kemena
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

#include "SequenceSet.h"

namespace BioTools
{
namespace Seq
{

using namespace std;
using namespace Utils;

SequenceSet::SequenceSet() {
	// TODO Auto-generated constructor stub

}

SequenceSet::SequenceSet(const std::string &seq_f)
{
	read(seq_f);
}

SequenceSet::~SequenceSet() {
	// TODO Auto-generated destructor stub
}


double
SequenceSet::avg_size() const
{
	size_t n_seqs = this->n_seqs();
	double length = 0;
	for (size_t i=0; i<n_seqs; ++i)
		length+= this->_seqs[i]->ungapped_size();
	return length/n_seqs;
}


size_t
SequenceSet::max_size() const
{
	size_t max_len = 0;
	size_t n_seqs = this->n_seqs();
	for (size_t i=0; i<n_seqs; ++i)
		if (_seqs[i]->size() > max_len)
			max_len=_seqs[i]->size();
	return max_len;
}

void
SequenceSet::delete_seqs(std::vector<size_t> &indices)
{
	std::sort(indices.begin(), indices.end());
	size_t n_dels = indices.size();
	size_t num_seqs = this->n_seqs();
	size_t index_pos = 0;
	size_t seq_pos = 0;
	for (size_t i = 0; i < num_seqs; ++i)
	{
		if ((index_pos==n_dels) || ((index_pos<n_dels) && (indices[index_pos] != i)))
			_seqs[seq_pos++] = _seqs[i];
		else
			++index_pos;
	}
	_seqs.resize(seq_pos);
}

void
SequenceSet::keep_seqs(std::vector<size_t> &indices)
{
	std::sort(indices.begin(), indices.end());
	size_t n_dels = indices.size();
	size_t num_seqs = this->n_seqs();
	size_t index_pos = 0;
	size_t seq_pos = 0;
	for (size_t i = 0; i < num_seqs; ++i)
	{
		if (indices[index_pos] == i)
		{
			_seqs[seq_pos++] = _seqs[i];
			++index_pos;
		}
		if (index_pos==n_dels)
			break;
	}
	_seqs.resize(seq_pos);
}


void
SequenceSet::delete_seqs(const std::map<std::string,bool> &names)
{
	size_t num_seqs = this->n_seqs();
	size_t seq_pos = 0;
	for (size_t i = 0; i < num_seqs; ++i)
	{
		if (names.count(_seqs[i]->name()))
			_seqs[seq_pos++] = _seqs[i];
	}
	_seqs.resize(seq_pos);
}



bool input_sort(Seq_ptr i, Seq_ptr j)
{
	return (i->id() < j->id());
}

bool seq_sort(Seq_ptr i, Seq_ptr j)
{
	return (*i<*j);
}

bool name_sort(Seq_ptr i, Seq_ptr j)
{
	return (i->name()<j->name());
}



void
SequenceSet::sort(std::string type)
{
	if (type == "input")
	{
		std::sort(_seqs.begin(), _seqs.end(), input_sort);
	}
	else if (type == "seq")
	{
		std::sort(_seqs.begin(), _seqs.end(), seq_sort);
	}
	else if (type == "name")
	{
		std::sort(_seqs.begin(), _seqs.end(), name_sort);
	}
}


bool
check_set(const SequenceSet &set)
{
	size_t n_seqs = set.n_seqs();
	size_t j, len;
	size_t *val_counting = new size_t[256];
	for (j=0; j<256; ++j)
		val_counting[j] = 0;

	// counting occurrences of characters
	for (size_t i=0; i<n_seqs; ++i)
	{
		const Sequence &seq = set[i];
		len=seq.size();
		for (j=0; j<len; ++j)
			++val_counting[static_cast<int>(seq[j])];
	}

	// see if strange character has been found
	for (j=0; j<45; ++j)
	{
		if (val_counting[j] != 0)
			return false;
	}
	for (j=46; j<65; ++j)
	{
		if (val_counting[j] != 0)
			return false;
	}
	for (j=91; j<97; ++j)
	{
		if (val_counting[j] != 0)
			return false;
	}
	for (j=123; j<256; ++j)
	{
		if (val_counting[j] != 0)
			return false;
	}
	return true;
}


}

}
