/*
 * 5_state_hmm.h
 *
 */

#ifndef HMM_H_
#define HMM_H_


// C header
#include <cctype>
#include <cfloat>

// C++ header
#include <cmath>
#include <unordered_map>

// BioTools header
#include "../utils/Matrix.h"
#include "../Sequence/Alignment.h"

namespace BioTools {
namespace HMM {



/**
 * The HMM class using probabilities from Probalign.
 */
class HMM {

private:
	BioTools::Utils::Matrix<float> _transProb;
	BioTools::Utils::Matrix<float> _insProb;
	BioTools::Utils::Matrix<float> _matchProb;
	float *_initDistr;
	short _num_states;
	short _num_ins_states;


public:
	/**
	 * \brief Constructor
	 * @param type The probabilities to use.
	 */
	HMM(char type);
	virtual ~HMM();

	/**
	 * \brief Returns the number of states.
	 * @return The number of states in the HMM.
	 */
	short num_states() const
	{
		return _num_states;
	}

	/**
	 * \brief The number of insert states.
	 * @return The number of insert states.
	 */
	short num_ins_states() const
	{
		return _num_ins_states;
	}

	/**
	 * Returns the transition probabilities.
	 * @return Transition probabilities.
	 */
	const BioTools::Utils::Matrix<float>&
	trans_probs() const
	{
		return _transProb;
	}

	/**
	 * \brief Returns the insertion probabilities.
	 * @return Insertion probabilities.
	 */
	const BioTools::Utils::Matrix<float>&
	ins_probs() const
	{
		return _insProb;
	}

	/**
	 * \brief Returns the match probabilities.
	 * @return Match probabilities.
	 */
	const BioTools::Utils::Matrix<float>&
	match_probs() const
	{
		return _matchProb;
	}


	/**
	 * \brief Returns initial distribution.
	 * @return Initial distribution
	 */
	const float *
	init_distribution() const
	{
		return _initDistr;
	}


	void
	aln_probs(const BioTools::Seq::Alignment &aln1, const BioTools::Seq::Alignment &aln2, BioTools::Utils::Matrix<double> &ins_probs1, BioTools::Utils::Matrix<double> &ins_probs2, BioTools::Utils::Matrix<double> &match_probs);

	void
	calculate_insertion_probs(const BioTools::Seq::Alignment &aln, BioTools::Utils::Matrix<float> &ins_probs);
	void
	calculate_match_probs(const BioTools::Seq::Alignment &aln1, const BioTools::Seq::Alignment &aln2, BioTools::Utils::Matrix<float> &match_probs);
};




}
}


#endif /* HMM_H_ */
