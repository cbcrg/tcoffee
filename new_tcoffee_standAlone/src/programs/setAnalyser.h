
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
