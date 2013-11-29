/*
 * setAnalyser.h
 *
 *  Created on: Jun 24, 2012
 *      Author: Carsten Kemena
 *
 *   Copyright 2011 Carsten Kemena
 *
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

#ifndef SETANALYSER_H_
#define SETANALYSER_H_


// C header
#include <cstdlib>
#include <cstdio>

// C++ header
#include <string>
#include <vector>

// boost header
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

// other headers
#include <omp.h>

// BioTools headers
#include "../lib/Sequence/Alignment.h"
#include "../lib/Sequence/aln_analysis.h"
#include "../lib/Sequence/SequenceSet.h"
#include "../lib/utils/filesystem.h"
#include "../lib/utils/graph_algorithms.h"
#include "../lib/utils/Matrix.h"
#include "../lib/utils/ScoringMatrix.h"


template <typename DataType>
void
composition(const DataType &aln, size_t *composition)
{
	size_t n_seqs=aln.n_seqs();
	size_t i,j, len;
	for (i=0; i<n_seqs; ++i)
	{
		const BioTools::Seq::Sequence seq = aln[i];
		len=seq.size();
		for (j=0; j<len; ++j)
		{
			if (seq[j] != '-')
				++composition[toupper(seq[j])];
		}
	}
}

template <typename DataType>
void
coverage(const DataType &aln, std::map<std::string, size_t> &genome_len)
{
	size_t n_seqs=aln.n_seqs();
	size_t i,j, len, ungap_len;
	std::map<std::string, size_t>::iterator it_end = genome_len.end();
	for (i=0; i<n_seqs; ++i)
	{
		const BioTools::Seq::Sequence seq = aln[i];
		if (genome_len.find(seq.name())==it_end)
			continue;
		len=seq.size();
		ungap_len=0;
		for (j=0; j<len; ++j)
			if (seq[j] != '-')
				++ungap_len;
		genome_len[seq.name()] += ungap_len;
	}
}


#endif /* SETANALYSER_H_ */
