/*
 * AlignmentBlock.h
 *
 *  Created on: Feb 12, 2012
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


#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_


// C header
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <cmath>

// C++ header
#include <vector>
#include <map>
#include <exception>


// Boost header
//#include <boost/algorithm/string.hpp>


#include "omp.h"

// BioTools++ header
#include "./Sequence.h"
#include "./SequenceSet.h"
#include "../utils/Matrix.h"
#include "../utils/string_functions.h"
#include "../utils/filesystem.h"


namespace BioTools {
namespace Seq {



/*! \file Alignment.h
    \brief Contains the Alignment class as well as a class for possible exceptions.
*/

/**
 * \brief A class to handle alignments.
 *
 * Alignment offers a method to read and write several different alignment methods as well as allowing some basic analysis and manipulation.
 */
class Alignment : public BioTools::Seq::SequenceSet
{
private:
	// Data
	size_t _aln_length;   		/// The length of the alignment


	// methods
	//deletes columns consisting of gaps only
	void _delete_gap_columns();

	void _write_fasta(FILE *aln_f, unsigned int line_break=60) const;
	void _write_fasta_seq(FILE *aln_F, unsigned int line_break=60) const;
	void _write_clustalw(FILE *aln_f, unsigned int line_break=60) const;
	void _write_msf(FILE *aln_f) const;
	void _write_phylip_sequential(FILE *aln_f) const;
	void _write_phylip_interleaved(FILE *aln_f) const;


	bool _len_check();
public:

	/**@{*/
	//! \name Constructors & Destructors

	/**
	 * \brief Alignment constructor
	 */
	Alignment();

	/**
	 * \brief Alignment constructor
	 */
	Alignment(Sequence *seq);

	/**
	 * \brief Alignment constructor
	 */
	Alignment(const SequenceSet &set, size_t id);

	/**
	 * \brief Alignment destructor
	 */
	virtual ~Alignment();
	/**@}*/


	//******   Methods   ******

	/**@{*/
	//! \name Basic methods


	/**
	 * \brief Returns the length of the alignment.
	 * @return Length of the alignment.
	 */
	size_t size() const throw()
	{
		return _aln_length;
	}


	/**@}*/



	//******   read/write alignment   ******
	/**@{*/
	//! \name Alignment IO

	/**
	 * \brief Reads the alignment
	 *
	 * @param aln_f The alignment file to read.
	 * @param check If true the alignment is checked if it is a proper alignment
	 * @param format The format of the alignment. (-1 enables automatic format detection)
	 */
	void read(const std::string &aln_f, bool check=false, short format =-1)
	{
		return read(aln_f, std::vector<std::string>(), check, format);
	}

	/**
	 * \brief Extracts a subalignment from an alignment.
	 *
	 * Only sequences which are denoted in seq_names are extracted. Names in seq_names not occurring in the alignment are ignored. Columns consisting of gaps only are removed.
	 * @param aln_f The alignment file to read.
	 * @param seq_names The names to read.
	 * @param check If true the alignment is checked if it is a proper alignment
	 * @param format The format of the alignment. (-1 enables automatic format detection)
	 */
	void read(const std::string &aln_f, const std::vector<std::string> &seq_names, bool check=false, short format =-1);

	/**
	 * \brief Writes the alignment into a file.
	 *
	 * This function supports the following formats: FASTA, MSF.
	 * @param aln_f The file to write the alignment to
	 * @param format The format to use (fasta, clustalw, msf, phylip_i, phylip_s)
	 */
	void write(const std::string &aln_f, const std::string format) const;
	/**@}*/



	//******   Manipulation methods   ******
	/**@{*/
	//! \name Manipulation methods

	/**
	 * \brief Deletes columns from the alignment containing more or equal gaps than a given percentage of columns.
	 *
	 * \param gap_percentage Percentage threshold of number of gaps.
	 */
	void column_trim(double gap_percentage);

	/**
	 * \brief Deletes columns from an alignment.
	 *
	 * Two modes exist depending on the length of the vector.
	 * \param col_to_delete A vector of two different types. If the vector has the same length as the alignment it is treated as a bool vector. All columns != 0 are deleted. If the vector is shorter the values denote the columns to delete.
	 *
	 */
	void delete_columns(std::vector<size_t> &col_to_delete);

	/**
	 * \brief Deletes sequences from the alignment
	 * \param names The names of the sequences to delete
	 */
	void delete_seqs(const std::map<std::string,bool> &names);

	/**
	 * \brief Deletes sequences from the alignment
	 * \param indices The indices of the sequences to delete.
	 */
	void delete_seqs(std::vector<size_t> &indices);

	/**
	 * \brief Deletes sequences if they are not in the given list.
	 * \param indices The indices of the sequences to keep.
	 */
	void keep_seqs(std::vector<size_t> &indices);


	/**
	 * \brief Merges two alignments.
	 * \param aln_gaps The gaps to insert into the alignment. The vector needs to be sorted py position.
	 * \param aln2 The alignment no merge into this one.
	 * \param aln_gaps2 The gaps in the second alignment.
	 */
	void
	merge(const std::vector<std::pair<unsigned int, unsigned int> > &aln_gaps, Alignment &aln2, const std::vector<std::pair<unsigned int, unsigned int> > &aln_gaps2);

	/**
	 * \brief Adds a Sequence to the alignment.
	 * \param aln_gaps The gaps to insert into the alignment. The vector needs to be sorted py position.
	 * \param seq The sequence to add to the alignment.
	 * \param seq_gaps The gaps to insert into the sequence.
	 */
	void
	add(const std::vector<std::pair<unsigned int, unsigned int> > &aln_gaps, const Sequence &seq, const std::vector<std::pair<unsigned int, unsigned int> > &seq_gaps);


};

void
replace_char(char c1, char c2, Alignment &aln);


} /* namespace Seq */
} /* namespace BioTools */

#endif /* ALIGNMENTBLOCK_H_ */
