/*
 * hmm_test.h
 *
 *  Created on: 30 Oct 2013
 *      Author: CarstenK
 *		 Email: CarstenK[@]posteo.de
 *	 Copyright: 2013
 *
 *  This file is part of BioTools.
 *
 *  BioTools is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  BioTools is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with BioTools.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef HMM_TEST_H_
#define HMM_TEST_H_



// C header
#include <cstdlib>

// CxxTest header
#include <cxxtest/TestSuite.h>

#include "../src/lib/Sequence/Alignment.h"
#include "../src/lib/HMM/HMM.h"
#include "../src/lib/HMM/fw_bw.h"

using namespace std;
using namespace BioTools;

class SkipList_Test : public CxxTest::TestSuite
{
public:

	void test_calc_match_probs()
	{
		Seq::Alignment aln1(new Seq::Sequence("seq1", "comment", "AAACGTAAA.AAAAT"));
		Seq::Sequence *seq = new Seq::Sequence("seq1", "comment", "AAACGTAAATAAAAA");
		aln1.addSeq(seq);

		Seq::Alignment aln2(new Seq::Sequence("seq1", "comment", "AAACGTAaA.AAAAA"));
		seq = new Seq::Sequence("seq1", "comment", "AAACGTAAATAAAAC");
		aln2.addSeq(seq);
		TS_ASSERT_EQUALS(aln2.n_seqs(), 2);

		HMM::HMM hmm('P');
		const Utils::Matrix<float> &probs= hmm.match_probs();
		BioTools::Utils::Matrix<float> match_probs;
		hmm.calculate_match_probs(aln1, aln2, match_probs);
		TS_ASSERT_EQUALS(match_probs[0][0], probs['A']['A']);
		float field_prob=(probs['A']['A']+probs['T']['A']+probs['T']['C']+probs['A']['C'])/4;
		TS_ASSERT_EQUALS(match_probs[14][14], field_prob);
	}

	void test_calc_ins_probs()
	{
		Seq::Alignment  aln1(new Seq::Sequence("seq1", "comment", "AAAcGTA.A.AAAAT"));
		Seq::Sequence *seq = new Seq::Sequence("seq1", "comment", "AAAcGT.AATAAAAA");
		aln1.addSeq(seq);

		Seq::Alignment aln2(new Seq::Sequence("seq1", "comment", "AAACGTAAA.AAAAA"));
		seq =               new Seq::Sequence("seq1", "comment", "AAACGTAAATAAAAC");
		aln2.addSeq(seq);
		TS_ASSERT_EQUALS(aln2.n_seqs(), 2);

		HMM::HMM hmm('P');
		const Utils::Matrix<float> &probs= hmm.ins_probs();
		BioTools::Utils::Matrix<float> insertion_probs;
		hmm.calculate_insertion_probs(aln1, insertion_probs);
		float field_prob=(probs['A'][0]+probs['A'][0])/2;
		TS_ASSERT_EQUALS(insertion_probs[0][0], field_prob);
		field_prob=(probs['A'][0]+probs['T'][0])/2;
		TS_ASSERT_EQUALS(insertion_probs[14][0], field_prob);
		field_prob=probs['T'][0];
		TS_ASSERT_EQUALS(insertion_probs[9][0], field_prob);


	}


	/*
	Seq::Alignment *
	produce_aln_from_profiles(vector<Seq::Alignment*> &aln_set, const Aln_opts &aln_opts, const Tree_opts tree_opts)
	{
		// pairwise alignments to build library
		if (verbose)
			printf("Calculate pairwise alignments\n");

		size_t n_alns=aln_set.size();

		Seq::Sequence *tmp_seq;
		Seq::SequenceSet consensus_set;
		size_t i;

		for (i=0; i<n_alns; ++i)
			aln_set[i]->write("arg","fasta");

		for (i=0; i<n_alns; ++i)
		{
			tmp_seq = BioTools::Seq::simple_consensus(*(aln_set[i]));
			tmp_seq->id(i);
			consensus_set.addSeq(tmp_seq);
		}

		Utils::Library lib(aln_set);
		Utils::Matrix<float> dist_mat(n_alns, n_alns);

		if (aln_opts.pair_aln_method=="hmm")
			HMM::all_hmm_pairs(aln_set, lib, dist_mat, 0, n_alns-1);
		//else if (aln_opts.pair_aln_method == "gotoh")
		//	all_gotoh_pairs(consensus_set, lib, aln_opts.score_mat, aln_opts.gapopen, aln_opts.gapext);


		// tree construction
		if (verbose)
			printf("Calculate/Read tree\n");

		//Matrix<float> dist_mat(n_seqs, n_seqs);

		Clustering::Tree tree;
		if (!tree_opts.tree_i_f.empty())
		{
			tree.read_newick(tree_opts.tree_i_f);
			tree.match_seq_set(consensus_set);
		}
		else
		{
			// calculate distance matrix if needed
			if (tree_opts.dist_mat_method=="kmer")
				Clustering::make_dist_mat(consensus_set, tree_opts.tree_k, tree_opts.alphabet, dist_mat);

			// calculate tree
			vector<std::string> names;
			for (i=0; i<n_alns; ++i)
				names.push_back(consensus_set[i].name());

			if (tree_opts.tree_method=="nj")
				tree.nj(dist_mat, names);
			else if (tree_opts.tree_method=="upgma")
				tree.upgma(dist_mat, names);
			else
			{
				fprintf(stderr, "Error! Tree method unknown\n");
				exit(EXIT_FAILURE);
			}
		}


		if (!tree_opts.tree_o_f.empty())
			tree.print_newick(tree_opts.tree_o_f);


		if (verbose)
			printf("Relaxing library\n");
		lib.relax();

		// calculate msa using scores from pairwise alignments
		std::map<Utils::Match, int> match_points;
		std::map<Utils::Match, int>::iterator it, it_end;

		Seq::Alignment *result_aln = progressive_aln(tree.root() , lib, aln_set);
		return result_aln;
	}

*/
};


#endif /* HMM_TEST_H_ */
