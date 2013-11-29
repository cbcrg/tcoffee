/*
 * Library.h
 *
 *  Created on: Apr 6, 2012
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
#ifndef LIBRARY_H_
#define LIBRARY_H_


// C++ header
#include <algorithm>
#include <limits>
#include <map>
#include <memory>
#include <utility>
#include <vector>


// BioTools++ header
#include "./Matrix.h"
#include "./fast_math.h"
#include "./filesystem.h"
#include "../Sequence/Sequence.h"
#include "../Sequence/SequenceSet.h"
#include "../Sequence/Alignment.h"

namespace BioTools {
namespace Utils {

typedef std::pair<unsigned int, unsigned int> Match;


struct Match_point
{
	Match_point():x(0),y(0),score(0)
	{}
	Match_point(unsigned int x_, unsigned int y_, float score_):x(x_),y(y_),score(score_)
	{}

	/**
	 * \brief Lesser or equal operator.
	 * \param a The first point.
	 * \param b The second point.
	 */
	friend bool operator <(const Match_point &a, const Match_point &b)
	{
		if (a.x != b.x)
			return (a.x<b.x);
		else
			return (a.y < b.y);
	}

	unsigned int x;
	unsigned int y;
	float score;
};


struct pair_compare
{
	bool
	operator()(const Match &a, const Match &b)
	{
		if (a.first<b.first)
			return true;
		return (a.second<b.second);
	}
};


/**
 * \brief A class to store pairwise matches.
 */
class Library
{

private:
	size_t _n_seqs;
	size_t _max_length;
	std::vector<std::vector<std::vector<Match_point> > > _pairs;
	std::vector<std::vector<Match_point> > _relaxed_pairs;
	omp_lock_t add_lock;
	omp_lock_t relax_lock;

public:

	// Constructors & Destructors

	/**@{*/
	//! \name Constructors & Destructors

	/**
	 * Constructor
	 * \param set The sequence set the library is constructed for.
	 */
	Library(const BioTools::Seq::SequenceSet &set);

	/**
	 * Destructor
	 */
	virtual ~Library();
	/**@}*/


	/**
	 * \brief Adds a matching pare to the library.
	 *
	 * If the pair does not exist yet in the library, the score will be set to 1. If the pair already exists, the score is increased by 1.
	 * \param seq1_id The sequence id of the first sequence.
	 * \param seq2_id The sequence id of the second sequence.
	 * \param pos1 The residue of sequence 1.
	 * \param pos2 The residue of sequence 2.
	 * \param score The score of the pair.
	 */
	void add(unsigned int seq1_id, unsigned int seq2_id, unsigned int pos1, unsigned int pos2, float score=1.0);

	/**
	 * \brief Calculates the scores for each residue pair.
	 * \param[in] seq1_id The id of sequence 1
	 * \param[in] seq2_id The id of sequence 2
	 * \param[out] match_points In this structure, the pairs with the scores will be saved.
	 */
	void get(size_t seq1_id, size_t seq2_id, std::map<Match, int> &match_points) const;

	/**
	 * \brief Calculates the scores for each residue pair.
	 *
	 * This function compares the matching scores for two alignments. The scores between every pair of sequence from the different alignmetns are computed
	 * and converted to fit the aligned sequences.
	 * \param[in] aln1 The first alignment.
	 * \param[in] aln2 The second alignment.
	 * \param[out] match_points In this structure, the pairs with the scores will be saved.
	 *
	 * \sa get
	 */
	void get(const BioTools::Seq::Alignment &aln1, const BioTools::Seq::Alignment &aln2, std::map<Match, int> &match_points) const;

	/**
	 * \brief Calculates the scores for each residue pair.
	 *
	 * The score between each sequence in the alignment and the single sequence is calculated.
	 * \param[in] aln The first alignment.
	 * \param[in] seq The sequence.
	 * \param[out] match_points In this structure, the pairs with the scores will be saved.
	 */
	void get(const BioTools::Seq::Alignment &aln, const BioTools::Seq::Sequence &seq, std::map<Match, int> &match_points) const
	{}

	void print_pairs()
	{

		size_t i,j,k;
		for (i=0; i<_pairs.size(); ++i)
			for (j=0; j<_pairs[i].size(); ++j)
				for (k=0; k<_pairs[i][j].size(); ++k)
					printf("%u %u %f\n",_pairs[i][j][k].x, _pairs[i][j][k].y, _pairs[i][j][k].score);
	}
	/**
	 * \brief Prints the library to a file.
	 *
	 * \param set The sequence set.
	 * \param out_f The file tow rite the library to.
	 */
	void print()//const BioTools::Seq::SequenceSet &set, std::string out_f)
	{
		size_t i,j;
		for (i=0; i<_relaxed_pairs.size(); ++i)
		{
			for (j=0; j<_relaxed_pairs[i].size(); ++j)
				printf("%u %u %f\n",_relaxed_pairs[i][j].x, _relaxed_pairs[i][j].y, _relaxed_pairs[i][j].score);
		}
	}



	void
	relax();


};









} /* namespace Utils */
} /* namespace BioTools */
#endif /* LIBRARY_H_ */
