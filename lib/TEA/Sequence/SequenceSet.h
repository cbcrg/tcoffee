/*
 * SequenceSet.h
 *
 *  Created on: Jan 28, 2012
 *      Author: ck
 */

#ifndef SEQUENCESET_H_
#define SEQUENCESET_H_

// C header
#include <cstdio>

// C++ header
#include <map>
#include <string>
#include <vector>

// BioTools header
#include "Sequence.h"
#include "Alignment_excep.h"
#include "../utils/filesystem.h"
#include "../utils/utils.h"
#include "../utils/Uncopyable.h"



namespace BioTools
{
namespace Seq
{


/*! \file SequenceSet.h
    \brief Contains the sequence class as well as a class for possible exceptions.
*/


/**
 * \brief A class to handle sequences.
 *
 * SequenceSet offers basic methods to handle sequences.
 */
class SequenceSet : private BioTools::Utils::Uncopyable
{

private:
	std::vector<Seq_ptr> _seqs; // saves the single sequences
	char _seq_type;				// Sequence type
	std::string _file;			// File the alignment was read from
	int _id;

	// read multiple sequence alignment file
	short _identify_aln_format(FILE *aln_F);
	void _read_fasta_f(FILE *aln_F, const std::map<std::string, short> &seq_names);
	void _read_clustalw_f(FILE *aln_F, const std::map<std::string, short> &seq_names);
	void _read_msf_f(FILE *aln_F, const std::map<std::string, short> &seq_names);
	void _read_stockholm_f(FILE *aln_F, const std::map<std::string, short> &seq_names);
	void _read_codata_f(FILE *aln_F, const std::map<std::string, short> &seq_names);
	void _read_amps_f(FILE *aln_F, const std::map<std::string, short> &seq_names);
	void _read_phylip_f(FILE *aln_F, const std::map<std::string, short> &seq_names);


	// write alignment to file
	void _write_fasta(FILE *aln_f, unsigned int line_break=60) const;



	void _gap_replace();
	bool _seq_check();

public:
	/**@{*/
	//! \name Constructors & Destructors
	/**
	 * \brief Basic constructor for this class.
	 */
	SequenceSet();

	/**
	 * \brief Constructor for this class which reads a file.
	 * \param seq_f The file to read.
	 */
	SequenceSet(const std::string &seq_f);


	/**
	 * \brief Basic destructor.
	 */
	virtual ~SequenceSet();
	/**@}*/

	void addSeq(Sequence *seq)
	{
		_seqs.push_back(Seq_ptr(seq));
	}


	// Operators
	/**@{*/
	//! \name Operators

	/**
	 * \brief Operator to access the sequence
	 * @param index The sequence position to return.
	 * @return Pointer to the sequence.
	 */
	Sequence &operator[](unsigned int index)
	{
		return *_seqs[index];
	}

	/**
	 * \overload
	 */
	const Sequence &operator[](unsigned int index) const
	{
		return *_seqs[index];
	}
	/**@}*/


	/**@{*/
	//! \name Basic methods

	const Sequence* seq(unsigned int index) const
	{
		return &(*(_seqs[index]));
	}


	/**
	 * \brief returns the number of sequences.
	 * \return The number of sequences.
	 */
	size_t n_seqs() const
	{
		return _seqs.size();
	}


	/**
	 * \brief The maximum size of the sequence set.
	 * \return The maximum size.
	 */
	size_t max_size() const;


	/**
	 * \brief The average size of the sequence set.
	 * \return The average size.
	 */
	double avg_size() const;

	/**
	 * \brief Returns true if no sequences are contained in this object.
	 */
	bool empty() const
	{
		return (_seqs.empty());
	}

	/**
	 * \brief Returns the file the sequences were read from.
	 */
	std::string file() const
	{
		return _file;
	}

	/**
	 * \brief Returns type.
	 * @return The sequence type.;
	 */
	char seq_type() const throw()
	{
		return _seq_type;
	}

	void seq_type(char seq_type_) throw()
	{
		_seq_type = seq_type_;
	}


	int id() const throw()
	{
		return _id;
	}

	void id(int val)
	{
		_id=val;
		size_t n_seqs= _seqs.size();
		for (size_t i=0; i<n_seqs; ++i)
			_seqs[i]->aln_id(val);
	}


	/**
	 * \brief Sets everything to 0.
	 */
	void clear()
	{
		_seqs.clear();
		_file.clear();
	}
	/**@}*/

	/**@{*/
	//! \name Read/Write methods


	/**
	 * \brief Reading a fasta file.
	 * @param seq_f The fasta file to read.
	 * @param seq_names Sequence name filter.
	 */
//	void read_fasta(const std::string &seq_f, const std::map<std::string, short> &seq_names);

	/**
	 * \brief Reading a fasta file.
	 * @param seq_f The fasta file to read.
	 */
	//void read_fasta(const std::string &seq_f);
	/**@}*/



	//******   read/write SequenceSets   ******
	/**@{*/
	//! \name SequenceSet IO

	/**
	 * \brief Reads a set of sequences.
	 *
	 * This function can read unaligned sequences in FASTA format as well as aligned sequences in several formats.
	 * @param seq_f The file with the sequences to read.
	 * @param format The format of the alignment. (-1 enables automatic format detection)
	 * @param check Checks if the sequence is a proper biological sequence
	 */
	virtual void read(const std::string &seq_f, bool check=false, short format =-1)
	{
		return read(seq_f, std::vector<std::string>(), check, format);
	}

	/**
	 * \brief Extracts a subalignment from the sequence set.
	 *
	 * Only sequences which are denoted in seq_names are extracted. Names in seq_names not occurring in the alignment are ignored. Columns consisting of gaps only are removed.
	 * @param seq_f The file of the sequences to read.
	 * @param seq_names The names to read.
	 * @param format The format of the alignment. (-1 enables automatic format detection)
	 * @param check Checks if the sequence is a proper biological sequence
	 */
	virtual void read(const std::string &seq_f, const std::vector<std::string> &seq_names, bool check=false, short format =-1);

	/**
	 * \brief Writes the sequences into a file.
	 *
	 * This function supports the following formats: FASTA, MSF.
	 * @param seq_f The file to write the alignment to
	 * @param format The format to use (fasta, clustalw, msf, phylip_i, phylip_s)
	 */
	virtual void write(const std::string &seq_f, const std::string format) const;
	/**@}*/


	//******   Manipulation methods   ******
	/**@{*/
	//! \name Manipulation methods


	/**
	 * \brief Turns all characters to uppercase.
	 */
	void to_upper()
	{
		for (size_t i = 0; i<_seqs.size(); ++i)
			_seqs[i]->to_upper();
	}

	/**
	 * \brief Turns all characters to lowercase.
	 */
	void to_lower()
	{
		for (size_t i = 0; i<_seqs.size(); ++i)
			_seqs[i]->to_lower();
	}



	/**
	 * \brief Deletes sequences from the alignment
	 * \param names The names of the sequences to delete
	 */
	virtual void delete_seqs(const std::map<std::string,bool> &names);

	/**
	 * \brief Deletes sequences from the alignment
	 * \param indices The indices of the sequences to delete.
	 */
	virtual void delete_seqs(std::vector<size_t> &indices);

	/**
	 * \brief Deletes sequences if they are not in the given list.
	 * \param indices The indices of the sequences to keep.
	 */
	virtual void keep_seqs(std::vector<size_t> &indices);


	void
	share(const SequenceSet &set, size_t id)
	{
		_seqs.push_back(Seq_ptr(set._seqs[id]));
	}


	void
	transfer(SequenceSet &set, size_t id)
	{
		_seqs.push_back(Seq_ptr(set._seqs[id]));
		set._seqs.erase(set._seqs.begin()+id);
	}

	void
	transfer(SequenceSet &set)
	{
		size_t n_seqs=set.n_seqs();
		_seqs.reserve(this->n_seqs()+n_seqs);
		for (size_t i=0; i<n_seqs; ++i)
			_seqs.push_back(Seq_ptr(set._seqs[i]));
		set.clear();
	}


	/**@}*/

	/**
	 * \brief Sorts the sequences.
	 * \param type "input" sorts the sequences by order of the input. "name" sorts by sequence name. "seq" sorts the sequences by alphabetical order.
	 */
	void sort(std::string type);
	/**@}*/

};

/**
 * \relates SequenceSet
 * @param set The set to check.
 * @return True if everything is fine.
 */
bool
check_set(const SequenceSet &set);

}
}

#endif /* SEQUENCESET_H_ */
