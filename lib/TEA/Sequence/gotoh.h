/*
 * gotoh.h
 *
 *  Created on: Feb 11, 2012
 *      Author: Carsten Kemena
 *
 *   Copyright 2012 Carsten Kemena
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

/*! \file gotoh.h
 *   \brief Contains several methods related to the gotoh algorithm.
 *
 *   The Gotoh algorithm has been first described in: O. Gotoh: An improved algorithm for matching biological sequences. In: Journal of Molecular Biology. 162, 1982, S. 705-708
 *
 */


//C headers
#include <climits>

// C++ headers
#include <stack>
#include <vector>
#include <map>

// Boost headers
#define BOOST_DISABLE_ASSERTS
#include "boost/multi_array.hpp"

// BioTools++ headers
#include "Sequence.h"
#include "Alignment.h"
#include "../utils/Library.h"
#include "../utils/ScoringMatrix.h"



#ifndef GOTOH_H_
#define GOTOH_H_



namespace BioTools
{


namespace Seq
{


/**@{*/
//! \name Gotoh alignment function

/**
 * \brief Typedef of a boost multi_array
 */
typedef boost::multi_array<int, 3> Gotoh_Int_array;

/**
 * \brief Tyepdef of a boost multi_array
 */
typedef boost::multi_array<short, 3> Gotoh_Trace_array;





/**
 * \brief aligns two sequences using the gotoh algorithm.
 * \param[in] seq1 The first sequence.
 * \param[in] seq2 The second sequence.
 * \param[in] dp_mat The dynamic programming matrix.
 * \param[in] trace_mat The trace matrix.
 * \param[in] score_mat The scoring matrix.
 * \param[in] gapopen The gap opening costs.
 * \param[in] gapext The gap extension cost.
 * \param[out] matches The matches as computed by the algorithm.
 * \pre dp_mat and trace_mat need to have at least the sizes [3][seq1_length+1][seq2_length+2]
 */
void
gotoh(const Sequence &seq1, const Sequence &seq2, Gotoh_Int_array &dp_mat, Gotoh_Trace_array &trace_mat, const BioTools::Utils::Scoring_Matrix &score_mat, int gapopen, int gapext);


/**
 * \brief Aligns two sequences using the gotoh algorithm.
 *
 * The matrices need to be of the correct size. It produces a trace matrix which can be turned into library pairs
 * using the gotoh_trace_all_optimal function.
 * \param[in] seq1 The first sequence.
 * \param[in] seq2 The second sequence.
 * \param[in] dp_mat The dynamic programming matrix.
 * \param[in] trace_mat The trace matrix.
 * \param[in] score_mat The scoring matrix.
 * \param[in] gapopen The gap opening costs.
 * \param[in] gapext The gap extension cost.
 * \param[out] matches The matches as computed by the algorithm.
 *
 */
void
gotoh_all_optimal(const Sequence &seq1, const Sequence &seq2, Gotoh_Int_array &dp_mat, Gotoh_Trace_array &trace_mat, const BioTools::Utils::Scoring_Matrix &score_mat, int gapopen, int gapext);

/**
 * \brief Turns the gotoh trace matrix into values for a library.
 * \param[in] seq_id1 The sequence id of sequence 1.
 * \param[in] seq_id2 The sequence id of sequence 2.
 * \param[in] trace_mat The trace matrix.
 * \param[out] lib The library to which the matches should be added.
 */
void
gotoh_trace(size_t seq_id1, size_t seq_id2, const Gotoh_Trace_array &trace_mat, BioTools::Utils::Library &lib);


/**
 * \brief Turns the gotoh trace matrix as produced by into values for a library.
 *
 * This trace back adds all matches contained in all possible optimal alignments into the library. Each match pair is only added
 * once even if it occurs in several optimal alignments.
 * \param[in] seq1 The sequence id of sequence 1.
 * \param[in] seq2 The sequence id of sequence 2.
 * \param[in] dp_mat The Programming matrix to use.
 * \param[in] trace_mat The trace matrix.
 * \param[out] lib The library to which the matches should be added.
 */
void
gotoh_trace_all_optimal(const Sequence &seq1, const Sequence &seq2, Gotoh_Int_array &dp_mat, Gotoh_Trace_array &trace_mat, BioTools::Utils::Library &lib);

/**
 * \brief Does a gotoh "chaining" to compute an alignment.
 */
void
gotoh_chaining(size_t l_seq1, size_t l_seq2, const std::map<BioTools::Utils::Match, int> &match_points, Gotoh_Int_array &dp_mat, Gotoh_Trace_array &trace_mat, int gapopen, int gapext);

/**
 * \brief Traces back of the chaining algorithm.
 */
void
gotoh_chain_trace(std::vector<std::pair<unsigned int, unsigned int> > &aln_gaps1, std::vector<std::pair<unsigned int, unsigned int> > &aln_gaps2, const Gotoh_Trace_array &trace_mat);



void
all_gotoh_pairs(const SequenceSet &set, BioTools::Utils::Library &lib, const BioTools::Utils::Scoring_Matrix &score_mat, int gapopen, int gapext);


/**@}*/

} // Seq
} // BioTools++


#endif /* GOTOH_H_ */
