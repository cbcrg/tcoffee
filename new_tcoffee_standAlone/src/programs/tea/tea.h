

#ifndef TEA_H_
#define TEA_H_

// C headers
#include <cstdlib>
#include <cstdio>

// C++ headers
#include <map>
#include <stack>
#include <string>


// boost header
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

// other headers
#include <omp.h>


// BioTools headers
#include "../../lib/Clustering/Tree.h"
#include "../../lib/Clustering/clustering.h"
#include "../../lib/Clustering/kmeans.h"
#include "../../lib/HMM/fw_bw.h"
#include "../../lib/Sequence/gotoh.h"
#include "../../lib/Sequence/aln_analysis.h"
#include "../../lib/utils/Library.h"
#include "../../lib/utils/Matrix.h"
#include "../../lib/utils/ScoringMatrix.h"



struct Tree_opts
{
	std::string tree_i_f;
	std::string tree_o_f;
	std::string dist_mat_method;
	std::string tree_method;
	std::string alphabet;
	unsigned short tree_k;
	unsigned short cluster_size;

};

struct Aln_opts
{
	std::string pair_aln_method;
	int gapopen;
	int gapext;
	std::string matrix_name;
	char seq_type;
	BioTools::Utils::Scoring_Matrix score_mat;
};

BioTools::Seq::Alignment *
progressive_aln(BioTools::Clustering::TreeNode *root, BioTools::Utils::Library &lib, std::vector<BioTools::Seq::Alignment*> alns);

void
prod_seq_aln(BioTools::Seq::SequenceSet &seq_set, Tree_opts &tree_opts, Aln_opts &aln_opts, BioTools::Seq::Alignment &result_aln);


#endif /* TEA_H_ */
