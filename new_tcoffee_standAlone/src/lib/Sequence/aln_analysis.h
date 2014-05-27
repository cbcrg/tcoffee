
/*! \file aln_analysis.h
    \brief Contains several methods to analysis an alignment.
*/


#ifndef ALN_ANALYSIS_H_
#define ALN_ANALYSIS_H_

#include <utility>


#include "Alignment.h"
#include "../utils/string_functions.h"
#include "../utils/ScoringMatrix.h"
#include "../utils/Matrix.h"

namespace BioTools
{
namespace Seq
{

//******   Analysis methods   ******
/**@{*/
//! \name Analysis methods

/**
 * \relates Alignment
 * \brief Calculates sequence identity of the alignment.
 *
 * For efficiency purposes the sequences are not compared pairwise. The alignments are turned into a profile and then
 * the correct number of matches and mismatches are calculated from these profiles.
 * @param aln An alignment.
 * @return A pair of values. The first is the number of identical pairs and the second the number of total pairs.
 */
std::pair<size_t, size_t> msaIdentity(const Alignment &aln);


/**
 * \relates Alignment
 * \brief Calculates the Sum-of-Pairs score of an alignment.
 * @param aln An alignment.
 * @param matrix The matrix to be used for the scoring of the alignment.
 * @param gop Gap opening costs.
 * @param gep Gap extension costs.
 * \return The Sum-of-Pairs score of the alignment.
 */
double sum_of_pairs_score(const Alignment &aln, const BioTools::Utils::Scoring_Matrix &matrix, double gop, double gep);
/**@}*/



/**@{*/
//! \name Alignment comparison

/**
 * \relates Alignment
 * \brief Calculates the column or sum-of-pairs similarity of two alignments.
 * @param ref The first alignment
 * @param test The second alignment
 * @param gap_limit The maximum allowed number of gaps in columns to consider.
 * @param mode The mode to do the comparison.
 * @param ignore_missing_seqs Ignores missing sequences in the test alignment.
 * @return The similarity between two alignments.
 */
std::pair<double, double> sim(const Alignment &ref, const Alignment &test, double gap_limit, const std::string &mode, bool ignore_missing_seqs);

/**@}*/



/**@{*/
//! \name Alignment consensus function

/**
 * \relates Alignment
 * \brief Returns a simple consensus derived from the alignment.
 * \param aln The alignment
 * \param threshold The threshold has to be passed to produce a consensus. If the score stays below the threshold the ambiguous character is given.
 * \param ambiguous The characters that should be inserted if the threshold is not passed.
 */
Sequence * simple_consensus(const Alignment &aln, float threshold=0.7, char ambiguous='N');

/**@}*/

}
}

#endif /* ALN_ANALYSIS_H_ */
