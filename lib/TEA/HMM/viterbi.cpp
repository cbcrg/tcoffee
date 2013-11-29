/*
 * viterbi.cpp
 *
 *  Created on: Dec 15, 2012
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

#include "viterbi.h"

using namespace Seq;
using namespace BioTools::Utils;
using namespace boost;


void
viterbi_all_optimal(const Seq::Sequence &seq1, const Seq::Sequence &seq2, const HMM &hmm, Viterbi_dp_array &dp_mat, Viterbi_trace_array &trace_mat, unsigned short n_states)
{
	size_t l_seq1 = seq1.size()+1;
	size_t l_seq2 = seq2.size()+1;
	size_t i,j;
	short k;
	const Matrix<float> &matchProbs = hmm.match_probs();
	const Matrix<float> &insProbs = hmm.ins_probs();
	const Matrix<float> &transProbs = hmm.trans_probs();
	const float *initDistr = hmm.init_distribution();
	//dp_mat.resize(extents[n_states][l_seq1][l_seq2]);
	// 0 matches, 1 insert1, 2 insert2, 3 del1 4 del2
	for (k=0; k<n_states; ++k)
		dp_mat[k][0][0]=initDistr[k];
	short c1,c2;
	for (i=1; i<l_seq1; ++i)
	{
		c1= static_cast<short>(seq1[i-1]);
		dp_mat[1][i][0] = dp_mat[1][i-1][0]+transProbs[1][1] + insProbs[c1][0];
		dp_mat[2][i][0] = -999999999;
		dp_mat[3][i][0] = dp_mat[3][i-1][0]+transProbs[3][3] + insProbs[c1][0];
		dp_mat[4][i][0] = -999999999;
		dp_mat[0][i][0] = max(dp_mat[1][i][0], dp_mat[3][i][0]);
	}


	for (i=1; i<l_seq2; ++i)
	{
		c2= static_cast<short>(seq2[i-1]);
		dp_mat[1][0][i] = -999999999;
		dp_mat[2][0][i] = dp_mat[2][0][i-1]+transProbs[2][2] + insProbs[c2][0];
		dp_mat[3][0][i] = -999999999;
		dp_mat[4][0][i] = dp_mat[4][0][i-1]+transProbs[4][4] + insProbs[c2][0];
		dp_mat[0][0][i] = max(dp_mat[2][0][i], dp_mat[4][0][i]);
	}

	// transition matrix: 0 match 1 insert1 2 deletion1 3 insert2 4 del2

	float tmp_score1, tmp_score2;
	float values[5];
	short traces[5];
	float tr[5] = {1, 2, 3, 8, 16};
	short l,m;
	float tmp;

	//State ID: match, insert1, del1, insert2, del2
	for (i=1; i<l_seq1; ++i)
	{
		c1=static_cast<short>(seq1[i-1]);
		for (j=1; j<l_seq2; ++j)
		{
			c2=static_cast<short>(seq2[j-1]);

			// insertion / deletion values
			for (k=1; k<n_states; ++k)
			{
				l = (k&1)?1:0;
				m = (k&1)?0:1;
				tmp = ((i==(l_seq1-1)) || (j==(l_seq2-1))) ?initDistr[k] : transProbs[0][k];
				tmp_score1 = dp_mat[0][i-l][j-m] + tmp;
				tmp_score2 = dp_mat[k][i-l][j-m] + transProbs[k][k];
				traces[k] = 0;
				if (tmp_score1 >= tmp_score2)
				{
					++traces[k];
					values[k] = tmp_score1;
				}
				if (tmp_score2 >= tmp_score1)
				{
					traces[k] += tr[k];
					values[k] = tmp_score2;
				}
				tmp = (l) ? insProbs[c1][0] : insProbs[c2][0];
				values[k] += tmp;
				dp_mat[k][i][j] = values[k];
				trace_mat[k][i][j] = traces[k];
			}

			// match value
			tmp_score1 = dp_mat[0][i-1][j-1] + transProbs[0][0] + matchProbs[c1][c2];
			traces[0] = 1;
			for (k=1; k<n_states; ++k)
			{
				if (values[k] > tmp_score1)
				{
					traces[0] = traces[k];
					tmp_score1 = values[k];
				}
				else if (values[k] == tmp_score1)
					traces[0] += traces[k];
			}
			dp_mat[0][i][j] = tmp_score1;
			trace_mat[0][i][j] = traces[0];

		}
	}

}


struct TracePosition
{
	TracePosition(short mat_, short pos1_, short pos2_) : mat(mat_), pos1(pos1_), pos2(pos2_)
	{};
	short mat;
	short pos1;
	short pos2;
};

double
viterbi_trace_all_optimal(const Seq::Sequence &seq1, const Seq::Sequence &seq2, Viterbi_dp_array &dp_mat, Viterbi_trace_array &trace_mat, BioTools::Utils::Library &lib)
{
	size_t seq_id1 = seq1.id();
	size_t seq_id2 = seq2.id();
	stack<TracePosition> to_do;
	size_t i = seq1.size();
	size_t j = seq2.size();
	to_do.push(TracePosition(0, i, j));

	short state=-1;
	int mat = 0;
	double total=0, match=0;
	while (!to_do.empty())
	{
		TracePosition &trace_pos = to_do.top();
		mat=trace_pos.mat;
		i=trace_pos.pos1;
		j=trace_pos.pos2;
		to_do.pop();
		state=1;
		while ((state!=0) && (i!=0) && (j!=0))
		{
			state = trace_mat[mat][i][j];
			if (state&1) // if state is odd then it is a match state
			{
				--trace_mat[mat][i][j];
				if ((state != 1) && (i!=0) && (j!=0))
					to_do.push(TracePosition(mat, i, j));
				if (mat!=0)
				{
					if ((mat==1) || (mat==3))
						--i;
					else
						--j;
					mat = 0;
				}
				else
				{
					--i;
					--j;
					if (dp_mat[0][i][j] != FLT_MAX)
					{
						dp_mat[0][i][j] = FLT_MAX;
						lib.add(seq_id1, seq_id2, i, j);
						if (seq1[i]==seq2[j])
							++match;
						++total;
					}
				}
			}
			else if ((state==2) || (state==6) || (state==10) || (state==14) || (state==18) || (state==22) || (state==30) )
			{
				trace_mat[mat][i][j]-=2;
				if (state !=2)
					to_do.push(TracePosition(mat, i, j));
				if (mat!=1)
					mat=1;
				else
					--i;
			}
			else if ((state==4) || (state==12) || (state==20) || (state==28) )
			{
				trace_mat[mat][i][j]-=4;
				if (state !=4)
					to_do.push(TracePosition(mat, i, j));
				if (mat!=3)
					mat=3;
				else
					--j;
			}
			else if ((state==8) || (state==24) )
			{
				trace_mat[mat][i][j]-=8;
				if (state !=8)
					to_do.push(TracePosition(mat, i, j));
				if (mat!=3)
					mat=3;
				else
					--i;
			}
			else if ((state))
			{
				if (mat!=4)
					mat=4;
				else
					--j;
			}

		}
	}

	if (total >= 0)
		return 100*match/total;
	else
		return 0;
}





