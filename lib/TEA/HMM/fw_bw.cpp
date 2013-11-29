/*
 * fw_bw.c
 *
 *  Created on: Apr 12, 2012
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

#include "fw_bw.h"


namespace BioTools {
namespace HMM {

using namespace std;
using namespace Seq;
using namespace BioTools::Utils;
using namespace boost;



/*
double logsumexp(double *nums, size_t ct)
{
	//double max_exp = nums[0], sum = 0.0;
	//size_t i;
	//for (i = 1; i < ct ; ++i)
		//if (abs(nums[i]) > abs(max_exp))
			//max_exp = nums[i];

	//for (i = 0; i < ct ; ++i)
//		sum += fmath::exp(nums[i] - max_exp);

//  return fmath::log(sum) + max_exp;
	return 0.0;
}*/



float
hmm_forward(const Seq::Sequence &seq1, const Seq::Sequence &seq2, const HMM &hmm, Matrix<float> &dp_mat, float **insert_matrices)
{
	unsigned short n_states=5;
	size_t l_seq1 = seq1.size()+1;
	size_t l_seq2 = seq2.size()+1;
	size_t i,j;
	float single_c1, single_c2;
	short k;
	const Matrix<float> &matchProbs = hmm.match_probs();
	const Matrix<float> &insProbs = hmm.ins_probs();
	const Matrix<float> &transProbs = hmm.trans_probs();
	const float *initDistr = hmm.init_distribution();

	// 0 matches, 1 insert1, 2 insert2, 3 del1 4 del2
	dp_mat[0][0] =initDistr[0];
	for (k=1; k<n_states; ++k)
		insert_matrices[2*k-2][0]=initDistr[k];

	// initialize y direction
	short c1,c2;
	float tmp1, tmp2;
	for (i=1; i<l_seq2; ++i)
	{
		c2= static_cast<short>(seq2[i-1]);
		single_c2 = insProbs[c2][0];
		tmp1 = insert_matrices[2][i] = LOG_ADD(insert_matrices[2][i-1]+transProbs[2][2], dp_mat[0][i-1]+transProbs[0][2]) + single_c2;
		tmp2 = insert_matrices[6][i] = LOG_ADD(insert_matrices[6][i-1]+transProbs[4][4], dp_mat[0][i-1]+transProbs[0][4]) + single_c2;
		dp_mat[0][i] = LOG_ADD(tmp1, tmp2);
	}

	// transition matrix: 0 match 1 insert1 2 deletion1 3 insert2 4 del2
	short l,m;
	float tmp;

	//State ID: match, insert1, del1, insert2, del2
	for (i=1; i<l_seq1; ++i)
	{
		c1=static_cast<short>(seq1[i-1]);
		single_c1 = insProbs[c1][0];

		// insert matrices
		tmp1 = insert_matrices[1][0] = LOG_ADD(insert_matrices[0][0]+transProbs[1][1], dp_mat[i-1][0]+transProbs[0][1]) + single_c1;
		tmp2 = insert_matrices[5][0] = LOG_ADD(insert_matrices[4][0]+transProbs[3][3], dp_mat[i-1][0]+transProbs[0][3]) + single_c1;
		dp_mat[i][0] = LOG_ADD(tmp1, tmp2);
		for (j=1; j<l_seq2; ++j)
		{
			c2=static_cast<short>(seq2[j-1]);
			single_c2 = insProbs[c2][0];

			// insertion / deletion values
			for (k=1; k<n_states; ++k)
			{
				l = (k&1)?1:0;
				m = (k&1)?0:1;
				tmp = (l) ? single_c1 : single_c2;
				if (((i==1) && (l==1)) || ((j==1) && (m==1)))
					insert_matrices[2*k-1][j] = dp_mat[i-l][j-m] + transProbs[0][k] + tmp;
				else
					insert_matrices[2*k-1][j] = LOG_ADD(dp_mat[i-l][j-m] + transProbs[0][k], insert_matrices[2*k-1-l][j-m] + transProbs[k][k]) + tmp;
			}

			// match value
			tmp = dp_mat[i-1][j-1] + transProbs[0][0];
			if ((i>1) && (j>1))
			{
				for (k=1; k<n_states; ++k)
					LOG_PLUS_EQUALS(&tmp, insert_matrices[2*k-2][j-1] + transProbs[k][0]);
			}
			else if (i>1)
			{
				LOG_PLUS_EQUALS(&tmp, insert_matrices[0][j-1] + transProbs[1][0]);
				LOG_PLUS_EQUALS(&tmp, insert_matrices[4][j-1] + transProbs[3][0]);
			}
			else
			{
				LOG_PLUS_EQUALS(&tmp, insert_matrices[2][j-1] + transProbs[2][0]);
				LOG_PLUS_EQUALS(&tmp, insert_matrices[6][j-1] + transProbs[4][0]);
			}
			dp_mat[i][j] = tmp  + matchProbs[c1][c2];
		}
		swap(insert_matrices[0], insert_matrices[1]);
		swap(insert_matrices[2], insert_matrices[3]);
		swap(insert_matrices[4], insert_matrices[5]);
		swap(insert_matrices[6], insert_matrices[7]);
	}
	float total=dp_mat[l_seq1-1][l_seq2-1];
	LOG_PLUS_EQUALS(&total, insert_matrices[0][l_seq2-1]);
	LOG_PLUS_EQUALS(&total, insert_matrices[2][l_seq2-1]);
	LOG_PLUS_EQUALS(&total, insert_matrices[4][l_seq2-1]);
	LOG_PLUS_EQUALS(&total, insert_matrices[6][l_seq2-1]);
	return total;
}

float
hmm_forward(const HMM &hmm, Matrix<double> &ins_probs1, Matrix<double> &ins_probs2, Matrix<double> &match_probs, Matrix<float> &dp_mat, float **insert_matrices)
{
	unsigned short n_states=5;
	size_t l_seq1 = ins_probs1.dim1()+1;
	size_t l_seq2 = ins_probs2.dim1()+1;
	size_t i,j;
	float single_c1, single_c2;
	short k;
	const Matrix<float> &transProbs = hmm.trans_probs();
	const float *initDistr = hmm.init_distribution();

	// 0 matches, 1 insert1, 2 insert2, 3 del1 4 del2
	dp_mat[0][0] =initDistr[0];
	for (k=1; k<n_states; ++k)
		insert_matrices[2*k-2][0]=initDistr[k];

	// initialize y direction
	float tmp1, tmp2;
	for (i=1; i<l_seq2; ++i)
	{
		//c2= static_cast<short>(seq2[i-1]);
		single_c2 = ins_probs2[i-1][0];
		tmp1 = insert_matrices[2][i] = LOG_ADD(insert_matrices[2][i-1]+transProbs[2][2], dp_mat[0][i-1]+transProbs[0][2]) + single_c2;
		tmp2 = insert_matrices[6][i] = LOG_ADD(insert_matrices[6][i-1]+transProbs[4][4], dp_mat[0][i-1]+transProbs[0][4]) + single_c2;
		dp_mat[0][i] = LOG_ADD(tmp1, tmp2);
	}

	// transition matrix: 0 match 1 insert1 2 deletion1 3 insert2 4 del2
	short l,m;
	float tmp;

	//State ID: match, insert1, del1, insert2, del2
	for (i=1; i<l_seq1; ++i)
	{
		single_c1 = ins_probs1[i-1][0];

		// insert matrices
		tmp1 = insert_matrices[1][0] = LOG_ADD(insert_matrices[0][0]+transProbs[1][1], dp_mat[i-1][0]+transProbs[0][1]) + single_c1;
		tmp2 = insert_matrices[5][0] = LOG_ADD(insert_matrices[4][0]+transProbs[3][3], dp_mat[i-1][0]+transProbs[0][3]) + single_c1;
		dp_mat[i][0] = LOG_ADD(tmp1, tmp2);
		for (j=1; j<l_seq2; ++j)
		{
			//cout << j << endl;
			single_c2 = ins_probs2[j-1][0];

			// insertion / deletion values
			for (k=1; k<n_states; ++k)
			{
				l = (k&1)?1:0;
				m = (k&1)?0:1;
				tmp = (l) ? single_c1 : single_c2;
				if (((i==1) && (l==1)) || ((j==1) && (m==1)))
					insert_matrices[2*k-1][j] = dp_mat[i-l][j-m] + transProbs[0][k] + tmp;
				else
					insert_matrices[2*k-1][j] = LOG_ADD(dp_mat[i-l][j-m] + transProbs[0][k], insert_matrices[2*k-1-l][j-m] + transProbs[k][k]) + tmp;
			}

			// match value
			tmp = dp_mat[i-1][j-1] + transProbs[0][0];
			if ((i>1) && (j>1))
			{
				for (k=1; k<n_states; ++k)
					LOG_PLUS_EQUALS(&tmp, insert_matrices[2*k-2][j-1] + transProbs[k][0]);
			}
			else if (i>1)
			{
				LOG_PLUS_EQUALS(&tmp, insert_matrices[0][j-1] + transProbs[1][0]);
				LOG_PLUS_EQUALS(&tmp, insert_matrices[4][j-1] + transProbs[3][0]);
			}
			else
			{
				LOG_PLUS_EQUALS(&tmp, insert_matrices[2][j-1] + transProbs[2][0]);
				LOG_PLUS_EQUALS(&tmp, insert_matrices[6][j-1] + transProbs[4][0]);
			}
			dp_mat[i][j] = tmp  + match_probs[i-1][j-1];
		}
		swap(insert_matrices[0], insert_matrices[1]);
		swap(insert_matrices[2], insert_matrices[3]);
		swap(insert_matrices[4], insert_matrices[5]);
		swap(insert_matrices[6], insert_matrices[7]);
	}
	float total=dp_mat[l_seq1-1][l_seq2-1];
	LOG_PLUS_EQUALS(&total, insert_matrices[0][l_seq2-1]);
	LOG_PLUS_EQUALS(&total, insert_matrices[2][l_seq2-1]);
	LOG_PLUS_EQUALS(&total, insert_matrices[4][l_seq2-1]);
	LOG_PLUS_EQUALS(&total, insert_matrices[6][l_seq2-1]);
	return total;
}


float
hmm_backward(const Seq::Sequence &seq1, const Seq::Sequence &seq2, const HMM &hmm, Matrix<float> &dp_mat, float **insert_matrices)
{
	unsigned short n_states=5;
	int l_seq1 = seq1.size();
	int l_seq2 = seq2.size();
	int i,j;
	float single_c1, single_c2;
	short k;
	const Matrix<float> &matchProbs = hmm.match_probs();
	const Matrix<float> &insProbs = hmm.ins_probs();
	const Matrix<float> &transProbs = hmm.trans_probs();
	const float *initDistr = hmm.init_distribution();

	// 0 matches, 1 insert1, 2 insert2, 3 del1, 4 del2
	dp_mat[l_seq1][l_seq2] =initDistr[0];
	for (k=1; k<n_states; ++k)
		insert_matrices[2*k-2][l_seq2]=initDistr[k];

	// initialize
	short c1,c2;
	float tmp1, tmp2;
	for (i=l_seq2-1; i>=0; --i)
	{
		c2= static_cast<short>(seq2[i]);
		single_c2 = insProbs[c2][0];
		tmp1 = insert_matrices[2][i] = LOG_ADD(insert_matrices[2][i+1]+transProbs[2][2], dp_mat[l_seq1][i+1]+transProbs[0][2]) + single_c2;
		tmp2 = insert_matrices[6][i] = LOG_ADD(insert_matrices[6][i+1]+transProbs[4][4], dp_mat[l_seq1][i+1]+transProbs[0][4]) + single_c2;
		dp_mat[l_seq1][i] = LOG_ADD(tmp1, tmp2);
	}

	// transition matrix: 0 match 1 insert1 2 deletion1 3 insert2 4 del2
	short l,m;
	float tmp;

	//State ID: match, insert1, del1, insert2, del2
	for (i=l_seq1-1; i>=0; --i)
	{
		c1=static_cast<short>(seq1[i]);
		single_c1 = insProbs[c1][0];

		// insert matrices
		tmp1 = insert_matrices[1][l_seq2] = LOG_ADD(insert_matrices[0][l_seq2]+transProbs[1][1], dp_mat[i+1][l_seq2]+transProbs[0][1]) + single_c1;
		tmp2 = insert_matrices[5][l_seq2] = LOG_ADD(insert_matrices[4][l_seq2]+transProbs[3][3], dp_mat[i+1][l_seq2]+transProbs[0][3]) + single_c1;
		dp_mat[i][l_seq2] = LOG_ADD(tmp1, tmp2);
		for (j=l_seq2-1; j>=0; --j)
		{
			c2=static_cast<short>(seq2[j]);
			single_c2 = insProbs[c2][0];

			// insertion / deletion values
			for (k=1; k<n_states; ++k)
			{
				l = (k&1)?1:0;
				m = (k&1)?0:1;
				tmp = (l) ? single_c1 : single_c2;
				if (((i==(l_seq1-1)) && (l==1)) || ((j==(l_seq2-1)) && (m==1)))
					insert_matrices[2*k-1][j] = dp_mat[i+l][j+m] + transProbs[0][k] + tmp;
				else
					insert_matrices[2*k-1][j] = LOG_ADD(dp_mat[i+l][j+m] + transProbs[0][k], insert_matrices[2*k-1-l][j+m] + transProbs[k][k]) + tmp;
			}

			// match value
			tmp = dp_mat[i+1][j+1] + transProbs[0][0];
			if ((i<l_seq1-1) && (j<l_seq2-1))
			{
				for (k=1; k<n_states; ++k)
					LOG_PLUS_EQUALS(&tmp, insert_matrices[2*k-2][j+1] + transProbs[k][0]);
			}
			else if (i<l_seq1-1)
			{
				LOG_PLUS_EQUALS(&tmp, insert_matrices[0][j+1] + transProbs[1][0]);
				LOG_PLUS_EQUALS(&tmp, insert_matrices[4][j+1] + transProbs[3][0]);
			}
			else
			{
				LOG_PLUS_EQUALS(&tmp, insert_matrices[2][j+1] + transProbs[2][0]);
				LOG_PLUS_EQUALS(&tmp, insert_matrices[6][j+1] + transProbs[4][0]);
			}
			dp_mat[i][j] = tmp  + matchProbs[c1][c2];
		}
		swap(insert_matrices[0], insert_matrices[1]);
		swap(insert_matrices[2], insert_matrices[3]);
		swap(insert_matrices[4], insert_matrices[5]);
		swap(insert_matrices[6], insert_matrices[7]);
	}
	float total=dp_mat[0][0];
	LOG_PLUS_EQUALS(&total, insert_matrices[0][0]);
	LOG_PLUS_EQUALS(&total, insert_matrices[2][0]);
	LOG_PLUS_EQUALS(&total, insert_matrices[4][0]);
	LOG_PLUS_EQUALS(&total, insert_matrices[6][0]);
	return total;
}

float
hmm_backward(const HMM &hmm, Matrix<double> &ins_probs1, Matrix<double> &ins_probs2, Matrix<double> &match_probs, Matrix<float> &dp_mat, float **insert_matrices)
{
	unsigned short n_states=5;
	int l_seq1 = ins_probs1.dim1();
	int l_seq2 = ins_probs2.dim2();
	int i,j;
	float single_c1, single_c2;
	short k;
	const Matrix<float> &transProbs = hmm.trans_probs();
	const float *initDistr = hmm.init_distribution();

	// 0 matches, 1 insert1, 2 insert2, 3 del1, 4 del2
	dp_mat[l_seq1][l_seq2] =initDistr[0];
	for (k=1; k<n_states; ++k)
		insert_matrices[2*k-2][l_seq2]=initDistr[k];

	// initialize
	float tmp1, tmp2;
	for (i=l_seq2-1; i>=0; --i)
	{
		single_c2 = ins_probs2[i][0];
		tmp1 = insert_matrices[2][i] = LOG_ADD(insert_matrices[2][i+1]+transProbs[2][2], dp_mat[l_seq1][i+1]+transProbs[0][2]) + single_c2;
		tmp2 = insert_matrices[6][i] = LOG_ADD(insert_matrices[6][i+1]+transProbs[4][4], dp_mat[l_seq1][i+1]+transProbs[0][4]) + single_c2;
		dp_mat[l_seq1][i] = LOG_ADD(tmp1, tmp2);
	}

	// transition matrix: 0 match 1 insert1 2 deletion1 3 insert2 4 del2
	short l,m;
	float tmp;

	//State ID: match, insert1, del1, insert2, del2
	for (i=l_seq1-1; i>=0; --i)
	{
		single_c1 = ins_probs1[i][0];

		// insert matrices
		tmp1 = insert_matrices[1][l_seq2] = LOG_ADD(insert_matrices[0][l_seq2]+transProbs[1][1], dp_mat[i+1][l_seq2]+transProbs[0][1]) + single_c1;
		tmp2 = insert_matrices[5][l_seq2] = LOG_ADD(insert_matrices[4][l_seq2]+transProbs[3][3], dp_mat[i+1][l_seq2]+transProbs[0][3]) + single_c1;
		dp_mat[i][l_seq2] = LOG_ADD(tmp1, tmp2);
		for (j=l_seq2-1; j>=0; --j)
		{
			single_c2 = ins_probs1[j][0];

			// insertion / deletion values
			for (k=1; k<n_states; ++k)
			{
				l = (k&1)?1:0;
				m = (k&1)?0:1;
				tmp = (l) ? single_c1 : single_c2;
				if (((i==(l_seq1-1)) && (l==1)) || ((j==(l_seq2-1)) && (m==1)))
					insert_matrices[2*k-1][j] = dp_mat[i+l][j+m] + transProbs[0][k] + tmp;
				else
					insert_matrices[2*k-1][j] = LOG_ADD(dp_mat[i+l][j+m] + transProbs[0][k], insert_matrices[2*k-1-l][j+m] + transProbs[k][k]) + tmp;
			}

			// match value
			tmp = dp_mat[i+1][j+1] + transProbs[0][0];
			if ((i<l_seq1-1) && (j<l_seq2-1))
			{
				for (k=1; k<n_states; ++k)
					LOG_PLUS_EQUALS(&tmp, insert_matrices[2*k-2][j+1] + transProbs[k][0]);
			}
			else if (i<l_seq1-1)
			{
				LOG_PLUS_EQUALS(&tmp, insert_matrices[0][j+1] + transProbs[1][0]);
				LOG_PLUS_EQUALS(&tmp, insert_matrices[4][j+1] + transProbs[3][0]);
			}
			else
			{
				LOG_PLUS_EQUALS(&tmp, insert_matrices[2][j+1] + transProbs[2][0]);
				LOG_PLUS_EQUALS(&tmp, insert_matrices[6][j+1] + transProbs[4][0]);
			}
			dp_mat[i][j] = tmp  + match_probs[i][j];
		}
		swap(insert_matrices[0], insert_matrices[1]);
		swap(insert_matrices[2], insert_matrices[3]);
		swap(insert_matrices[4], insert_matrices[5]);
		swap(insert_matrices[6], insert_matrices[7]);
	}
	float total=dp_mat[0][0];
	LOG_PLUS_EQUALS(&total, insert_matrices[0][0]);
	LOG_PLUS_EQUALS(&total, insert_matrices[2][0]);
	LOG_PLUS_EQUALS(&total, insert_matrices[4][0]);
	LOG_PLUS_EQUALS(&total, insert_matrices[6][0]);
	return total;
}



struct hmm_match
{
	float match;
	size_t x;
	size_t y;

	hmm_match(float m, size_t x_, size_t y_):match(m), x(x_), y(y_)
	{}

	friend bool operator <(const hmm_match& a, const hmm_match& b)
	{
		return(a.match > b.match);
	}
};



void
hmm2lib(const Sequence &seq1, const Sequence &seq2, const Matrix<float> &forward_mat, const Matrix<float> &backward_mat, unsigned short n_states, BioTools::Utils::Library &lib, float total_probability)
{
	size_t seq_id1 = seq1.id();
	size_t seq_id2 = seq2.id();
	size_t l_seq1 = seq1.size();
	size_t l_seq2 = seq2.size();
	size_t i,j;
	vector<hmm_match> matches;
	matches.reserve(l_seq1*l_seq2);
	float tmp;

	for (i=0; i<l_seq1; ++i)
	{
		for (j=0; j<l_seq2; ++j)
		{
			if ((tmp=EXP(min(LOG_ONE,(forward_mat[i][j] + backward_mat[i][j] - total_probability)))) >= 0.01)
				matches.push_back(hmm_match(tmp, i,j));
		}
	}

	// sorting

	sort(matches.begin(), matches.end());
	size_t max2 = min(4*min(l_seq1,l_seq2), matches.size());
	float min_score=matches[max2-1].match;
	i=0;
	while ((matches.size()!=i)&&(matches[i].match>=min_score))
	{
		lib.add(seq_id1, seq_id2, matches[i].x, matches[i].y, matches[i].match);
		++i;
	}
}

void
hmm2lib(const Sequence &seq1, int id1, const Sequence &seq2, int id2, const Matrix<float> &forward_mat, const Matrix<float> &backward_mat, unsigned short n_states, BioTools::Utils::Library &lib, float total_probability)
{
	size_t l_seq1 = seq1.size();
	size_t l_seq2 = seq2.size();
	size_t i,j;
	vector<hmm_match> matches;
	matches.reserve(l_seq1*l_seq2);
	float tmp;

	for (i=0; i<l_seq1; ++i)
	{
		for (j=0; j<l_seq2; ++j)
		{
			if ((tmp=EXP(min(LOG_ONE,(forward_mat[i][j] + backward_mat[i][j] - total_probability)))) >= 0.01)
				matches.push_back(hmm_match(tmp, i,j));
		}
	}

	// sorting

	sort(matches.begin(), matches.end());
	size_t max2 = min(4*min(l_seq1,l_seq2), matches.size());
	float min_score=matches[max2-1].match;
	i=0;
	while ((matches.size()!=i)&&(matches[i].match>=min_score))
	{
		lib.add(id1, id2, matches[i].x, matches[i].y, matches[i].match);
		++i;
	}
}


void
hmm2lib(const Alignment &aln1, const Alignment &aln2, const Matrix<float> &forward_mat, const Matrix<float> &backward_mat, unsigned short n_states, BioTools::Utils::Library &lib, float total_probability)
{
	size_t seq_id1 = aln1.id();
	size_t seq_id2 = aln2.id();
	size_t l_seq1 = aln1.size();
	size_t l_seq2 = aln2.size();
	size_t i,j;
	vector<hmm_match> matches;
	matches.reserve(l_seq1*l_seq2);
	float tmp;

	for (i=0; i<l_seq1; ++i)
	{
		for (j=0; j<l_seq2; ++j)
		{
			if ((tmp=EXP(min(LOG_ONE,(forward_mat[i][j] + backward_mat[i][j] - total_probability)))) >= 0.01)
				matches.push_back(hmm_match(tmp, i,j));
		}
	}

	// sorting
	sort(matches.begin(), matches.end());
	size_t max2 = min(4*min(l_seq1,l_seq2), matches.size());
	float min_score=matches[max2-1].match;
	i=0;
	while ((matches.size()!=i)&&(matches[i].match>=min_score))
	{
//		cout << seq_id1 << " " << seq_id2 << " " << matches[i].x << " " << matches[i].y << " " << matches[i].match << endl;
		lib.add(seq_id1, seq_id2, matches[i].x, matches[i].y, matches[i].match);
		++i;
	}
}



void
all_hmm_pairs(const BioTools::Seq::SequenceSet &set, BioTools::Utils::Library &lib, Matrix<float> &dist_mat)
{
	size_t max_len = set.max_size()+1;
	size_t i,j;
	size_t n_seqs = set.n_seqs();
	//HMM hmm(set.seq_type());
	HMM hmm('P');

	for (i=0; i< n_seqs; ++i)
		dist_mat[i][i]=0;

	#pragma omp parallel shared(set, n_seqs, max_len, hmm, dist_mat, lib) private(i, j)
	{
		float bw_p, fw_p=0;
		Matrix<float> forward_mat = Matrix<float>(max_len, max_len);
		Matrix<float> backward_mat = Matrix<float>(max_len, max_len);
		float **insert_matrices = new float*[8];
		for (i=0; i<8; ++i)
			insert_matrices[i] = new float[max_len];

		#pragma omp for nowait
		for (i=0; i< n_seqs; ++i)
		{
			for (j=i+1; j<n_seqs; ++j)
			{
				const Sequence &seq1 = set[i];
				const Sequence &seq2 = set[j];
				fw_p = hmm_forward(seq1, seq2, hmm, forward_mat, insert_matrices);
				bw_p = hmm_backward(seq1, seq2, hmm, backward_mat, insert_matrices);
				hmm2lib(seq1, i, seq2, j, forward_mat, backward_mat, 5, lib, (bw_p+fw_p)/2);
				dist_mat[i][j] = dist_mat[j][i] = 1-EXP(((fw_p+bw_p)/2));
			}
		}
	} //omp
}


void
all_hmm_pairs(const vector<BioTools::Seq::Alignment*> &set, BioTools::Utils::Library &lib, Matrix<float> &dist_mat, size_t start, size_t end)
{
	size_t n_seqs=end-start+1;
	size_t max_len = 0;
	size_t i,j;
	//cout << "N: " << n_seqs << endl;
	for (i=0; i<n_seqs; ++i)
	{
		if (max_len<set[i]->size())
			max_len=set[i]->size();
	}
	++max_len;
	HMM hmm(set[0]->seq_type());


	for (i=0; i< n_seqs; ++i)
		dist_mat[i][i]=0;

//	#pragma omp parallel shared(set, n_seqs, max_len, hmm, dist_mat, lib) private(i, j)
	{
		float bw_p, fw_p=0;
		Matrix<float> forward_mat = Matrix<float>(max_len, max_len);
		Matrix<float> backward_mat = Matrix<float>(max_len, max_len);
		float **insert_matrices = new float*[8];
		for (i=0; i<8; ++i)
			insert_matrices[i] = new float[max_len];
		Matrix<double> ins_probs1(1,1);
		Matrix<double> ins_probs2(1,1);
		Matrix<double> match_probs(1,1);
		#pragma omp for nowait

		for (i=start; i< end; ++i)
		{
			for (j=i+1; j<end; ++j)
			{
				const Alignment *aln1 = set[i];
				const Alignment *aln2 = set[j];
				hmm.aln_probs(*aln1, *aln2, ins_probs1, ins_probs2, match_probs);
				fw_p = hmm_forward(hmm, ins_probs1, ins_probs2, match_probs, forward_mat, insert_matrices);
				bw_p = hmm_backward(hmm, ins_probs1, ins_probs2, match_probs, backward_mat, insert_matrices);
				hmm2lib(*aln1, *aln2, forward_mat, backward_mat, 5, lib, (bw_p+fw_p)/2);

				dist_mat[i-start][j-start] = dist_mat[j-start][i-start] = 1-EXP(((fw_p+bw_p)/2));
			}
		}
	} //omp
}


void
all_hmm_pairs(const SequenceSet &set, BioTools::Utils::Library &lib, Matrix<float> &dist_mat, size_t start, size_t end)
{
	size_t n_seqs=end-start+1;
	size_t max_len = 0;
	size_t i,j;
	for (i=0; i<n_seqs; ++i)
	{
		if (max_len<set[i].size())
			max_len=set[i].size();
	}
	++max_len;
	HMM hmm(set.seq_type());

	for (i=0; i< n_seqs; ++i)
		dist_mat[i][i]=0;

	#pragma omp parallel shared(set, n_seqs, max_len, hmm, dist_mat, lib) private(i, j)
	{
		float bw_p, fw_p=0;
		Matrix<float> forward_mat = Matrix<float>(max_len, max_len);
		Matrix<float> backward_mat = Matrix<float>(max_len, max_len);
		float **insert_matrices = new float*[8];
		for (i=0; i<8; ++i)
			insert_matrices[i] = new float[max_len];

		#pragma omp for nowait
		for (i=start; i< end; ++i)
		{
			for (j=i+1; j<end; ++j)
			{
				const Sequence &seq1 = set[i];
				const Sequence &seq2 = set[j];
				fw_p = hmm_forward(seq1, seq2, hmm, forward_mat, insert_matrices);
				bw_p = hmm_backward(seq1, seq2, hmm, backward_mat, insert_matrices);
				hmm2lib(seq1, seq2, forward_mat, backward_mat, 5, lib, (bw_p+fw_p)/2);
				dist_mat[i-start][j-start] = dist_mat[j-start][i-start] = 1-EXP(((fw_p+bw_p)/2));
			}
		}
	} //omp
}



}
}

