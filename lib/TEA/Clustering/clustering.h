/*
 * clustering.h
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

#ifndef CLUSTERING_H_
#define CLUSTERING_H_


/*! \file clustering.h
    \brief Algorithms to turn sequences into a vector.
*/


// C++ header files
#include <vector>

// Boost header
#include <boost/shared_ptr.hpp>

// BioTools header files
#include "Vector.h"
#include "../Sequence/SequenceSet.h"
#include "../utils/utils.h"

namespace BioTools
{
namespace Clustering
{


typedef boost::shared_ptr<Vector<double> > Vec_double_ptr;
typedef boost::shared_ptr<Vector<unsigned int> > Vec_uint_ptr;
typedef boost::shared_ptr<BioTools::Seq::Sequence > Seq_ptr;

/**
 * \brief Turns a sequence into a vector.
 * @param seq The sequence.
 * @param k The size of the k-mers
 * @param factor The factors needed.
 * @param vec_len The length of the vector.
 * @param vec_num The number of vectors.
 * @param alphabet The alphabet to use.
 * @return A vector
 */
Vec_double_ptr
seq2vec_kmer(const BioTools::Seq::Sequence &seq, short k, unsigned int *factor, size_t vec_len, size_t vec_num, short *alphabet );

/**
 * \brief Turns a sequence set into a set of vectors.
 * @param seq_set The sequence set.
 * @param k The size of k-mers.
 * @param alphabet The alphabet to use.
 * @return The set of vectors.
 */
std::vector<Vec_double_ptr>*
seqset2vecs_kmer(const BioTools::Seq::SequenceSet &seq_set, short k, const std::string &alphabet);


void
reduce(std::vector<Vec_double_ptr> &vec_set);


}
}



#endif /* CLUSTERING_H_ */
