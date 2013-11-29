/*
 * vierbi.h
 *
 *  Created on: Dec 15, 2012
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

#ifndef VIERBI_H_
#define VIERBI_H_

#include "HMM.h"
#include "../Sequence/SequenceSet.h"


//typedef boost::multi_array<float, 3> Viterbi_dp_array;
//typedef boost::multi_array<float, 2> hmm_dp_array;
//typedef boost::multi_array<short, 3> Viterbi_trace_array;

/**
 * The viterbi algorithm to calculate the most probable alignment.
 * \param seq1 The first sequence.
 * \param seq2 The second sequence.
 * \param hmm The HMM to use.
 * \param dp_mat The dynamic programming matrix.
 * \param trace_mat The trace matrix.
 * \param n_states The number of states.
 */
//void viterbi_all_optimal(const Seq::Sequence &seq1, const Seq::Sequence &seq2, const HMM &hmm, Viterbi_dp_array &dp_mat, Viterbi_trace_array &trace, unsigned short n_states);

//double
//viterbi_trace_all_optimal(const Seq::Sequence &seq1, const Seq::Sequence &seq2, Viterbi_dp_array &dp_mat, Viterbi_trace_array &trace_mat, BioTools::Utils::Library &lib);

#endif /* VIERBI_H_ */
