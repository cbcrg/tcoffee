/*
 * maskHelper.cpp
 *
 *  Created on: May 18, 2012
 *      Author: Carsten Kemena
 *
 * This file is part of BioTools.
 *
 * BioTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BioTools is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BioTools.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


// C header
#include <cstdlib>
#include <cctype>
#include <cstdio>

// C++ header
#include <string>
#include <vector>

// boost header
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

// other headers
#include <omp.h>

// BioTools headers
#include "../lib/Sequence/SequenceSet.h"
#include "../lib/utils/filesystem.h"


using namespace std;
using namespace BioTools::Seq;

namespace po = boost::program_options;
namespace fs = boost::filesystem;


void
merge(SequenceSet &set1, const SequenceSet &set2)//, bool merge_soft, bool merge_hard)
{
	if (set1.n_seqs() != set2.n_seqs())
	{
		fprintf(stderr, "ERROR! %s and %s contain a different number of sequences!\n", set1.file().c_str(), set2.file().c_str());
		exit(EXIT_FAILURE);
	}

	size_t n_seqs=set1.n_seqs();
	size_t seq_length, j;
	size_t i;
	int chunk =1;
	#pragma omp parallel shared(n_seqs) private(i,j, seq_length)
	{
		#pragma omp for schedule(dynamic,chunk) nowait
		for (i=0; i<n_seqs; ++i)
		{
			Sequence &seq1=set1[i];
			const Sequence &seq2=set2[i];
			if ((seq1.name() != seq2.name()) || (seq1.size() != seq2.size()))
			{
				fprintf(stderr, "ERROR! Sequence name/length of %s (%s) and %s (%s) do not agree!\n", seq1.name().c_str(), set1.file().c_str(), seq2.name().c_str(), set2.file().c_str());
				exit(EXIT_FAILURE);
			}
			seq_length = seq1.size();
			for (j=0; j<seq_length; ++j)
			{
				if (seq2[j] == 'N')
					seq1[j] = 'N';
				else if (islower(seq2[j]))
					seq1[j]=seq2[j];
			}
		}
	}

}

void
hard_mask(SequenceSet &set)
{
	size_t n_seqs=set.n_seqs();
	size_t seq_length, j;
	for (size_t i=0; i<n_seqs; ++i)
	{
		Sequence &seq=set[i];
		seq_length = seq.size();
		for (j=0; j<seq_length; ++j)
		{
			if (islower(seq[j]))
				seq[j]='N';
		}
	}
}

void
soft_mask(const SequenceSet &original, SequenceSet &set)
{
	if (original.n_seqs() != set.n_seqs())
	{
		fprintf(stderr, "ERROR! %s and %s contain a different number of sequences!\n", set.file().c_str(), original.file().c_str());
		exit(EXIT_FAILURE);
	}

	size_t n_seqs=set.n_seqs();
	size_t seq_length, j;
	for (size_t i=0; i<n_seqs; ++i)
	{
		Sequence &seq=set[i];
		const Sequence &ori_seq = original[i];
		if ((seq.name() != ori_seq.name()) || (seq.size() != ori_seq.size()))
		{
			fprintf(stderr, "ERROR! Sequence name/length of %s (%s) and %s (%s) do not agree!\n", seq.name().c_str(), set.file().c_str(), ori_seq.name().c_str(), original.file().c_str());
			exit(EXIT_FAILURE);
		}
		seq_length = seq.size();
		for (j=0; j<seq_length; ++j)
		{
			if (seq[j] == 'N')
				seq[j]= tolower(ori_seq[j]);
		}
	}
}


int
main(int argc, char *argv[])
{
	int n_threads = 1;
	vector<string> input_list;
	string out_f ="";
	string original_seqs;
	string soft_mask_f;
	bool do_hard_mask;

	po::options_description general("General options");
	general.add_options()
		("help,h", "Produces this help message")
		("num_threads,n", po::value<int>(&n_threads), "Number of threads to use.")
		("out_file,o", po::value<string>(&out_f), "The file to write the result into.")
		("in,i", po::value<vector<string> >(&input_list)->multitoken(), "A list of input files.")
	;

	po::options_description mask("Mask options");
	mask.add_options()
		("soft_mask,s", po::value<string>(&soft_mask_f), "Turn hardmasked sequences into softmasking using this sequences.")
		("hard_mask,a", po::value<bool>(&do_hard_mask)->default_value(false)->zero_tokens(), "Turns all lower case characters into 'N'")
	;

	po::options_description all("maskHelper v1.0.\n\nAllowed options are displayed below.");
	all.add(general).add(mask);
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(all).run(), vm);
	po::notify(vm);

	if (vm.count("help") || input_list.empty()) {
		cout << all<< "\n";
		return EXIT_SUCCESS;
	}

	omp_set_num_threads(n_threads);
	SequenceSet result, tmp;
	result.read(input_list[0], 1);
	size_t n_sets = input_list.size();
	for (size_t i=1; i<n_sets; ++i)
	{
		tmp.read(input_list[i], 1);
		merge(result, tmp);
	}

	if (!soft_mask_f.empty())
	{
		SequenceSet original;
		original.read(soft_mask_f, 1);
		soft_mask(original, result);
	} else if (do_hard_mask)
		hard_mask(result);

	result.write(out_f,"fasta");


	return EXIT_SUCCESS;
}
