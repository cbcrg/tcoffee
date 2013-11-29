/*
 * Alignment.cpp
 *
 *  Created on: Oct 9, 2011
 *      Author: Carsten Kemena
 *
 *   Copyright 2011 Carsten Kemena
 *
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


#include "Alignment.h"

using namespace std;

namespace BioTools
{
namespace Seq
{



Alignment::Alignment() : _aln_length(0)
{
}

Alignment::Alignment(Sequence *seq_):_aln_length(seq_->size())
{
	this->addSeq(seq_);
	this->seq_type(identify_seq_type(*seq_));
}

Alignment::Alignment(const SequenceSet &set, size_t id):_aln_length(set[id].size())
{
	this->share(set, id);
	//this->addSeq(new Sequence(set[id]));
	this->seq_type(set.seq_type());
}


Alignment::~Alignment()
{

}


// Manipulation methods

/*
SequenceSet*
Alignment::del_gaps()
{
	SequenceSet *seqSet = new SequenceSet();

	return seqSet;
}*/

void Alignment::_delete_gap_columns()
{
	unsigned int num_seqs = this->n_seqs();
	size_t aln_length = size();
	unsigned int j;
	unsigned int *is_gap_column = new unsigned [aln_length];
	for (j=0; j<aln_length; ++j)
		is_gap_column[j] = 0;
	for (unsigned int i=0; i<num_seqs; ++i)
	{
		const Sequence &tmp_seq = (*this)[i];
		for (j=0; j<aln_length; ++j)
			if (tmp_seq[j] == '-')
				++is_gap_column[j];
	}
	size_t pos = 0;
	for (unsigned int i=0; i<num_seqs; ++i)
	{
		Sequence &tmp_seq = (*this)[i];
		pos=0;
		for (j=0; j<aln_length; ++j)
		{
			if (is_gap_column[j] != num_seqs)
				tmp_seq[pos++] = tmp_seq[j];
		}
		tmp_seq.resize(pos);
	}

	_aln_length = pos;
	delete[] is_gap_column;
}


void
Alignment::delete_seqs(const std::map<std::string,bool> &names)
{
	SequenceSet::delete_seqs(names);
	_delete_gap_columns();
}

void
Alignment::delete_seqs(std::vector<size_t> &indices)
{
	SequenceSet::delete_seqs(indices);
	_delete_gap_columns();
}



void
Alignment::delete_columns(vector<size_t> &col_to_delete)
{
	unsigned int num_seqs = this->n_seqs();
	size_t j;
	if (col_to_delete.size() == _aln_length)
	{
		size_t pos = 0;
		for (unsigned int i=0; i<num_seqs; ++i)
		{
			Sequence &tmp_seq = (*this)[i];
			pos=0;
			for (j=0; j<_aln_length; ++j)
			{
				if (!col_to_delete[j])
					tmp_seq[pos++] = tmp_seq[j];
			}
			tmp_seq.resize(pos);
		}
		_aln_length = pos;
	}
	else
	{
		std::sort(col_to_delete.begin(), col_to_delete.end());
		size_t pos = 0;
		size_t col_pos = 0;
		size_t n_cols_to_del = col_to_delete.size();
		for (unsigned int i=0; i<num_seqs; ++i)
		{
			Sequence &tmp_seq = (*this)[i];
			pos=0;
			col_pos=0;
			for (j=0; j<_aln_length; ++j)
			{
				if ((col_pos >= n_cols_to_del) || (col_to_delete[col_pos] != j))
					tmp_seq[pos++] = tmp_seq[j];
				else
					++col_pos;
			}
			tmp_seq.resize(pos);
		}
		_aln_length = pos;
	}
}

void
Alignment::column_trim(double gap_percentage)
{
	size_t num_seqs = this->n_seqs();
	size_t aln_length = size();

	size_t i,j;
	size_t *n_gaps = new size_t[aln_length];
	for (j=0; j<aln_length; ++j)
		n_gaps[j]=0;
	for (i=0; i<num_seqs; ++i)
	{
		Sequence &tmp_seq = (*this)[i];
		for (j=0; j<aln_length; ++j)
		{
			if (tmp_seq[j] == '-')
				++n_gaps[j];
		}
	}

	vector<size_t> del_columns;
	for (j=0; j<aln_length; ++j)
		if (((double)n_gaps[j]/(double)num_seqs) >= gap_percentage)
			del_columns.push_back(j);
	delete[] n_gaps;
	delete_columns(del_columns);
}



void
Alignment::merge(const vector<pair<unsigned int, unsigned int> > &aln_gaps, Alignment &aln, const vector<pair<unsigned int, unsigned int> > &aln_gaps2)
{
	size_t num_seqs = this->n_seqs();
	for (size_t i=0; i<num_seqs; ++i)
		(*this)[i].insert_gaps(aln_gaps);
	num_seqs=aln.n_seqs();
	for (size_t i=0; i<num_seqs; ++i)
		aln[i].insert_gaps(aln_gaps2);
	_aln_length=(*this)[0].size();

	this->transfer(aln);
}

void
Alignment::add(const vector<pair<unsigned int, unsigned int> > &aln_gaps, const Sequence &new_seq, const vector<pair<unsigned int, unsigned int> > &seq_gaps)
{
	size_t num_seqs = this->n_seqs();
	for (size_t i=0; i<num_seqs; ++i)
		(*this)[i].insert_gaps(aln_gaps);
	_aln_length=(*this)[0].size();
	Sequence *tmp_seq = new Sequence(new_seq);
	tmp_seq->insert_gaps(seq_gaps);
	this->addSeq(tmp_seq);
}


void
Alignment::keep_seqs(std::vector<size_t> &indices)
{
	SequenceSet::keep_seqs(indices);
	_delete_gap_columns();
}


void
replace_char(char c1, char c2, Alignment &aln)
{
	size_t n_seqs = aln.n_seqs();
	size_t j;
	size_t aln_length = aln.size();
	for (size_t i = 0; i< n_seqs; ++i)
	{
		Sequence &seq = aln[i];
		for (j=0; j<aln_length; ++j)
		{
			if (seq[j] == c1)
				seq[j] = c2;
		}
	}
}

} // end namespace Sequence
} // end namespace Biotools


