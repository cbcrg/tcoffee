/*
 * Sequence.h
 *
 *  Created on: Oct 9, 2011
 *      Author: Carsten Kemena
 *
 *   Copyright 2011 Carsten Kemena
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

#ifndef SEQUENCE_H_
#define SEQUENCE_H_

// C++ header
#include <map>
#include <string>
#include <vector>
#include <utility>

// Boost header
#include <boost/shared_ptr.hpp>


// My header
#include "../utils/string_functions.h"

/*! \file Sequence.h
    \brief Contains the Sequence class.
*/


namespace BioTools
{


namespace Seq
{


/**
 * \brief A class to handle biological sequences.
 *
 * It offers basic information storage capabilities.
 */
class Sequence
{
private:
	std::string _name; //Saves the name of the sequence
	std::string _comment; //Saves a comment
	std::string _sequence; //saves the sequence itself
	size_t _input_id;
	size_t _aln_id;


public:
	// Constructors & Destructors
	/**@{*/
	//! \name Constructors & Destructors
	/**
	 * \brief A simple sequence constuctor.
	 * @param seq_name The sequence name.
	 * @param comment_ A sequence comment.
	 * @param seq The sequence.
	 * @param seq_id The sequence id to give to the sequence.
	 */
	Sequence(const std::string &seq_name, const std::string &comment_, const std::string &seq, size_t seq_id=0);
	/**
	 * \brief A sequence constructor which allows the reservation of memory.
	 *
	 * The memory for the actual sequence will be reserved but the actual sequence can be larger without causing problems.
	 * @param seq_name The sequence name.
	 * @param comment_ A sequence comment.
	 * @param seq_length The sequence size for which memory should be reserved.
	 * @param seq_id The sequence id to give to the sequence.
	 */
	Sequence(const std::string &seq_name, const std::string &comment_, unsigned int seq_length, size_t seq_id=0);

	/**
	 * \brief The copy constructor.
	 *
	 * @param seq The sequence of wich a copy should be made.
	 */
	Sequence(const Sequence &seq);

	/**
	 * \brief Destructor
	 */
	virtual ~Sequence();
	/**@}*/



	/**@{*/
	//! \name %Operators

	/**
	 * \brief Access operator.
	 * \param index The index of to access.
	 */
	char &operator[](unsigned int index)
	{
		return _sequence[index];
	}

	/**
	 * \overload
	 */
	const char &operator[](unsigned int index) const
	{
		return _sequence[index];
	}

	/**
	 * \brief Comparison operators.
	 * \param a The first sequence.
	 * \param b The second sequence.
	 */
	friend bool operator ==(const Sequence& a, const Sequence& b)
	{
		return(a.sequence() == b.sequence());
	}

	/**
	 * \brief Comparison operators.
	 * \param a The first sequence.
	 * \param b The second sequence.
	 */
	friend bool operator !=(const Sequence& a, const Sequence& b)
	{
		return(a.sequence() != b.sequence());
	}

	/**
	 * \brief Lesser operator.
	 * \param a The first sequence.
	 * \param b The second sequence.
	 */
	friend bool operator <(const Sequence& a, const Sequence& b)
	{
		return(a.sequence() < b.sequence());
	}

	/**
	 * \brief Greater operator.
	 * \param a The first sequence.
	 * \param b The second sequence.
	 */
	friend bool operator >(const Sequence& a, const Sequence& b)
	{
		return(a.sequence() > b.sequence());
	}

	/**
	 * \brief Lesser or equal operator.
	 * \param a The first sequence.
	 * \param b The second sequence.
	 */
	friend bool operator <=(const Sequence& a, const Sequence& b)
	{
		return(a.sequence() <= b.sequence());
	}

	/**
	 * \brief Greater or equal operator
	 * \param a The first sequence.
	 * \param b The second sequence.
	 */
	friend bool operator >=(const Sequence& a, const Sequence& b)
	{
		return(a.sequence()==b.sequence());
	}
	/**@}*/


	//Functions

	/**@{*/
	//! \name Basic methods

	/**
	 * \brief Returns the name of the sequence
	 * \return The name.
	 */
	const
	std::string	& name() const
	{
		return _name;
	}

	/**
	 * \brief Returns the sequence as a string.
	 * \return The sequence.
	 */
	const
	std::string	& sequence() const
	{
		return _sequence;
	}

	/**
	 * \brief Returns the comment of a string.
	 * \return The comment.
	 */
	const
	std::string	& comment() const
	{
		return _comment;
	}

	/**
	 * \brief The length of the sequnce.
	 * \return Length of the sequence.
	 */
	size_t size() const
	{
		return _sequence.size();
	}

	/**
	 * \brief The size of the sequence without gaps.
	 * \return The size of the ungapped sequence.
	 */
	size_t ungapped_size() const;

	/**
	 * \brief Returns the sequence id of the sequence.
	 *
	 * It represents the order in which the sequences were read.
	 * \return The id of the sequence.
	 */
	size_t id() const
	{
		return _input_id;
	}

	void id(size_t val)
	{
		_input_id=val;
	}

	size_t aln_id() const
	{
		return _aln_id;
	}

	void aln_id(size_t val)
	{
		_aln_id = val;
	}
	/**@}*/



	/**@{*/
	//! \name Manipulation methods


	/**
	 * \brief Appends a string to the sequence.
	 */
	template<class T>
	void append(const T &seq)
	{
		_sequence.append(seq);
	}

	/**
	 * \brief Appends a char to the sequence.
	 */
	void append(char c)
	{
		_sequence.push_back(c);
	}

	/**
	 * \brief Resizes the sequence.
	 * \param new_length The new_length.
	 */
	void
	resize(unsigned int new_length)
	{
		_sequence.resize(new_length);
	}

	/**
	 * \brief Turns all charachters of the sequence to uppercase.
	 */
	void
	to_upper()
	{
		std::string::iterator it, it_end =_sequence.end();
		for (it = _sequence.begin(); it != it_end; ++it)
			*it = toupper(*it);
	}

	/**
	 * \brief Turns all characters of the sequence to lowercase.
	 */
	void
	to_lower()
	{
		std::string::iterator it, it_end =_sequence.end();
		for (it = _sequence.begin(); it != it_end; ++it)
			*it = tolower(*it);
	}


	/**
	 * \brief Inserts gaps into the sequence.
	 * \param vec A vector of pairs. The first value of a pair gives the position of the gap, the second one the length of it.
	 */
	void
	insert_gaps(const std::vector<std::pair<unsigned int, unsigned int> > vec);
	/**@}*/
};




typedef boost::shared_ptr<Sequence> Seq_ptr;



/**@{*/
//! \name Analyze function


/**
 * \relates Sequence
 * \brief Identifies the sequence type.
 * \param seq The sequence.
 * \return N for nucleotide sequence, P for protein sequence.
 */
char identify_seq_type(const Sequence &seq);

/**
 * \relates Sequence
 * \brief Returns the coverage of two sequences.
 * \param seq1 The first sequence.
 * \param seq2 The second sequence.
 * \return (-1, -1) if sequences do not have same length else a pair, first the number of matching/mismatching nucleotide pairs, the second the general number of nucleotide pairs.
 */
std::pair<size_t, size_t> coverage(const Sequence &seq1, const Sequence &seq2);

/**
 * \relates Sequence
 * \brief Returns the identity of two sequence.
 * \param seq1 The first sequence.
 * \param seq2 The second sequence.
 * \return (-1, -1) if sequences do not have same length else a pair, first the number of identical nucleotide pairs, the second the general number of nucleotide pairs.
 */
std::pair<size_t, size_t> id(const Sequence &seq1, const Sequence &seq2);


/**
 * \relates Sequence
 * \brief Compares two sequences, ignoring gaps.
 * @param seq1 The first sequence.
 * @param seq2 The second sequence.
 * @return true, if the gapless sequences are the same.
 */
bool
seq_check(const Sequence &seq1, const Sequence &seq2);


/**
 * \relates Sequence
 * \brief Checks if the sequence is a prober biological sequence
 * \param seq The sequence
 * \return true if it is a proper biological sequence
 */
bool
bio_seq(const Sequence &seq);

/**@}*/
} // end namespace Sequence
} // end namespace Biotools


#endif /* SEQUENCE_H_ */
