/*
 * fw_bw.h
 *
 *  Created on: Apr 12, 2012
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

/*! \file viterbi.h
    \brief The forward and backward algorithms.
*/

#ifndef FW_BW_H_
#define FW_BW_H_



// C header
#include <cstdlib>
#include <cfloat>


// C++ header
#include <stack>
#include <vector>
#include <cmath>

// Boost header
//#define BOOST_DISABLE_ASSERTS
//#include "boost/multi_array.hpp"

// external header
//#include "../../external_src/fmath/fmath.hpp"

// BioTools++ header
#include "HMM.h"
#include "../Sequence/Sequence.h"
#include "../Sequence/SequenceSet.h"
#include "../utils/fast_math.h"
#include "../utils/Library.h"
#include "../utils/Matrix.h"



namespace BioTools {
namespace HMM {



/**
 * Turns the results of the backward and forward algorithms into matches.
 * @param seq1 The first sequence.
 * @param seq2 The second sequence.
 * @param forward_mat The forward matrix.
 * @param backward_mat The backward matrix.
 * @param lib The library to use.
 * @param total_probability The total probability of the alignment.
 */
void
hmm2lib(const Seq::Sequence &seq1, const Seq::Sequence &seq2, const BioTools::Utils::Matrix<float> &forward_mat, const BioTools::Utils::Matrix<float> &backward_mat, BioTools::Utils::Library &lib, float total_probability);
void
hmm2lib(const Seq::Alignment &aln1, const Seq::Alignment &aln2, const Utils::Matrix<float> &forward_mat, const Utils::Matrix<float> &backward_mat, unsigned short n_states, BioTools::Utils::Library &lib, float total_probability);

/**
 * Calculates all pairwise alignments between each pair of sequences.
 * @param set The sequence set.
 * @param lib The library to which the matches should be added.
 */
void
all_hmm_pairs(const BioTools::Seq::SequenceSet &set, BioTools::Utils::Library &lib, BioTools::Utils::Matrix<float> &dist_mat);

void
all_hmm_pairs(const std::vector<BioTools::Seq::Alignment*> &set, BioTools::Utils::Library &lib, BioTools::Utils::Matrix<float> &dist_mat, size_t start, size_t end);

/*
void
all_hmm_pairs(const std::vector<BioTools::Seq::Alignment> &set, BioTools::Utils::Library &lib, BioTools::Utils::Matrix<float> &dist_mat)
{
	all_hmm_pairs(set, lib, dist_mat, 0, set.size());
}
*/

/**
 * \brief Computes the probability using the forward algorithm.
 *
 * The insert matrices are of dimension [2,seq_length]. Memory is reused during computation.
 * @param[in] seq1 The first sequence.
 * @param[in] seq2 The second sequence.
 * @param[in] hmm The hmm to use.
 * @param[out] dp_mat The dynamic programming matrix
 * @param[out] insert_matrices The insert matrices.
 * @return The total probability.
 */
float
hmm_forward(const Seq::Sequence &seq1, const Seq::Sequence &seq2, const HMM &hmm, BioTools::Utils::Matrix<float> &dp_mat, float **insert_matrices);
float
hmm_forward(const HMM &hmm, BioTools::Utils::Matrix<double> &ins_probs1, BioTools::Utils::Matrix<double> &ins_probs2, BioTools::Utils::Matrix<double> &match_probs, BioTools::Utils::Matrix<float> &dp_mat, float **insert_matrices);

/**
 * \brief Computes the probability using the backward algorithm.
 *
 * The insert matrices are of dimension [2,seq_length]. Memory is reused during computation
 * @param[in] seq1 The first sequence.
 * @param[in] seq2 The second sequence.
 * @param[in] hmm The hmm to use.
 * @param[out] dp_mat The dynamic programming matrix
 * @param[out] insert_matrices The insert matrices.
 * @return The total probability.
 */
float
hmm_backward(const Seq::Sequence &seq1, const Seq::Sequence &seq2, const HMM &hmm, BioTools::Utils::Matrix<float> &dp_mat, float **insert_matrices);
float
hmm_backward(const HMM &hmm, BioTools::Utils::Matrix<double> &ins_probs1, BioTools::Utils::Matrix<double> &ins_probs2, BioTools::Utils::Matrix<double> &match_probs, BioTools::Utils::Matrix<float> &dp_mat, float **insert_matrices);

//void
//hmm2lib(const Seq::Sequence &seq1, const Seq::Sequence &seq2, const hmm_dp_array &forward_mat, const hmm_dp_array &backward_mat, unsigned short n_states, BioTools::Utils::Library &lib);


double logsumexp(double *nums, size_t ct);

/*
inline
double logsumexp(double a, double b)
{
 // double max_exp = abs(a)>abs(b)? a : b;
  return 0.0;
 // return (fmath::log(fmath::exp(a - max_exp) + fmath::exp(b-max_exp)) +max_exp);
}*/












} // namespace HMM

} // namespace BioTools




#endif /* FW_BW_H_ */
