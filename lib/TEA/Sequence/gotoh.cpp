/*
 * pw_align.cpp
 *
 *  Created on: Feb 11, 2012
  *      Author: Carsten Kemena
 *
 *   Copyright 2012 Carsten Kemena
 *
 *
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


#include "gotoh.h"



namespace BioTools
{
namespace Seq
{

using namespace std;
using namespace boost;
using namespace BioTools::Utils;

void
gotoh(const Sequence &seq1, const Sequence &seq2, Gotoh_Int_array &dp_mat, Gotoh_Trace_array &trace_mat, const BioTools::Utils::Scoring_Matrix &score_mat, int gapopen, int gapext)
{
	size_t l_seq1 = seq1.size()+1;
	size_t l_seq2 = seq2.size()+1;
	dp_mat.resize(extents[3][l_seq1][l_seq2]);
	trace_mat.resize(extents[3][l_seq1][l_seq2]);

	// initialization
	// 0=match, 1=insert, 2=deletion
	dp_mat[0][0][0]=0;
	dp_mat[1][0][0]=gapopen;
	dp_mat[2][0][0]=gapopen;
	size_t i,j;
	for (i=1; i<l_seq1; ++i)
	{
		dp_mat[1][i][0]=INT_MIN;
		dp_mat[0][i][0]=dp_mat[2][i][0]=dp_mat[2][i-1][0]+gapext;
	}
	for (j=1; j<l_seq2; ++j)
	{
		dp_mat[0][0][j]=dp_mat[1][0][j]=dp_mat[1][0][j-1]+gapext;
		dp_mat[2][0][j]=INT_MIN;
	}

	// fill matrix
	short c1, c2;
	for (i=1; i<l_seq1; ++i)
	{
		c1 = tolower(seq1[i-1])-97;
		for (j=1; j<l_seq2; ++j)
		{
			//calculate insert value
			if (dp_mat[1][i-1][j] > dp_mat[0][i-1][j] +gapopen)
			{
				trace_mat[1][i][j] = 'i';
				dp_mat[1][i][j] = dp_mat[1][i-1][j];
			}
			else
			{
				trace_mat[1][i][j] = 'm';
				dp_mat[1][i][j] = dp_mat[0][i-1][j]+gapopen;
			}
			dp_mat[1][i][j] += gapext;

			//calculate deletion value
			if (dp_mat[2][i][j-1] > dp_mat[0][i][j-1] +gapopen)
			{
				trace_mat[2][i][j] = 'd';
				dp_mat[2][i][j] = dp_mat[2][i][j-1];
			}
			else
			{
				trace_mat[2][i][j] = 'm';
				dp_mat[2][i][j] = dp_mat[0][i][j-1]+gapopen;
			}
			dp_mat[2][i][j] += gapext;

			//calculate match value
			if (dp_mat[1][i][j] > dp_mat[2][i][j])
			{
				trace_mat[0][i][j] = 'i';
				dp_mat[0][i][j] = dp_mat[1][i][j];
			}
			else
			{
				trace_mat[0][i][j] = 'd';
				dp_mat[0][i][j] = dp_mat[2][i][j];
			}

			c2 = tolower(seq2[j-1])-97;
			if (dp_mat[0][i-1][j-1] + score_mat[c1][c2] >= dp_mat[0][i][j])
			{
				trace_mat[0][i][j] = 'm';
				dp_mat[0][i][j] = dp_mat[0][i-1][j-1] + score_mat[c1][c2];
			}
		}
	}
}



void
gotoh_all_optimal(const Sequence &seq1, const Sequence &seq2, Gotoh_Int_array &dp_mat, Gotoh_Trace_array &trace_mat, const BioTools::Utils::Scoring_Matrix &score_mat, int gapopen, int gapext)
{
	size_t l_seq1 = seq1.size()+1;
	size_t l_s1 = seq1.size();
	size_t l_s2 = seq2.size();
	size_t l_seq2 = seq2.size()+1;


	// initialization
	// 0=match, 1=insert, 2=deletion
	dp_mat[0][0][0]=0;
	dp_mat[1][0][0]=0;
	dp_mat[2][0][0]=0;
	size_t i,j;
	for (i=1; i<l_seq1; ++i)
	{
		dp_mat[1][i][0]=INT_MIN;
		dp_mat[0][i][0]=dp_mat[2][i][0]=dp_mat[2][i-1][0]+gapext;
	}
	for (j=1; j<l_seq2; ++j)
	{
		dp_mat[0][0][j]=dp_mat[1][0][j]=dp_mat[1][i][j-1]+gapext;
		dp_mat[2][0][j]=INT_MIN;
	}

	// fill matrix
	short c1, c2;
	short trace_value;
	int fin_score, ins_score, del_score, tmp_score1, tmp_score2;
	for (i=1; i<l_seq1; ++i)
	{
		c1 = tolower(seq1[i-1])-97;
		for (j=1; j<l_seq2; ++j)
		{
			//calculate insert value
			tmp_score1 = dp_mat[1][i-1][j];
			tmp_score2 = (i == l_s1) ? dp_mat[0][i-1][j] : dp_mat[0][i-1][j] + gapopen;
			trace_value=0;
			if (tmp_score1 >= tmp_score2)
			{
				trace_value = 2;
				ins_score= tmp_score1+gapext;
			}

			if (tmp_score2 >= tmp_score1)
			{
				++trace_value;
				ins_score = tmp_score2+gapext;
			}
			dp_mat[1][i][j] = ins_score;
			trace_mat[1][i][j] = trace_value;

			//calculate deletion value
			tmp_score1 = dp_mat[2][i][j-1];
			tmp_score2 = (j == l_s2) ? dp_mat[0][i][j-1] : dp_mat[0][i][j-1] +gapopen;
		//	if (i==l_seq1)
			//	tmp_score1-=gapext;
			trace_value=0;
			if (tmp_score1 >= tmp_score2)
			{
				trace_value = 4;
				del_score = tmp_score1+gapext;
			}

			if (tmp_score2 >= tmp_score1)
			{
				++trace_value;
				del_score = tmp_score2+gapext;
			}
			dp_mat[2][i][j] = del_score;
			trace_mat[2][i][j] = trace_value;


			//calculate match value
			trace_value =0;
			if (ins_score >= del_score)
			{
				trace_value += 2;
				fin_score = ins_score;
			}

			if (del_score >= ins_score)
			{
				trace_value += 4;
				fin_score = del_score;
			}

			c2 = tolower(seq2[j-1])-97;
			if (dp_mat[0][i-1][j-1] + score_mat[c1][c2] == fin_score)
				++trace_value;
			else if (dp_mat[0][i-1][j-1] + score_mat[c1][c2] > fin_score)
			{
				trace_value=1;
				fin_score = dp_mat[0][i-1][j-1] + score_mat[c1][c2];
			}
			dp_mat[0][i][j] = fin_score;
			trace_mat[0][i][j] = trace_value;

		}


	}


}


void
gotoh_trace(size_t seq_id1, size_t seq_id2, const Gotoh_Trace_array &trace_mat, Library &lib)
{
	size_t i = trace_mat[0].size()-1;
	size_t j = trace_mat[0][0].size()-1;
	char state='c';
	int mat = 0;
	while ((i!=0) && (j!=0))
	{
		state = trace_mat[mat][i][j];
		if (state=='m')
		{
			if (mat!=0)
			{
				if (mat==1)
					--i;
				else
					--j;
				mat = 0;
			}
			else
			{
				lib.add(seq_id1, seq_id2, --i, --j);
				//printf("A %li %li\n", i, j);
			}
		}
		else if (state=='i')
		{
			if (mat==0)
				mat=1;
			else
				--i;
		}
		else if (state=='d')
		{
			if (mat==0)
				mat=2;
			else
				--j;
		}
	}
	//printf("\n");
}


struct TracePosition
{
	TracePosition(short mat_, short pos1_, short pos2_) : mat(mat_), pos1(pos1_), pos2(pos2_)
	{};
	short mat;
	short pos1;
	short pos2;
};



void
gotoh_trace_all_optimal(const Sequence &seq1, const Sequence &seq2, Gotoh_Int_array &dp_mat, Gotoh_Trace_array &trace_mat, Library &lib)
{
	size_t seq_id1 = seq1.id();
	size_t seq_id2 = seq2.id();
	stack<TracePosition> to_do;
	size_t i = seq1.size();
	size_t j = seq2.size();
	to_do.push(TracePosition(0, i, j));

	short state=-1;
	int mat = 0;
//	double total=0, match=0;
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
					if (mat==1)
						--i;
					else
						--j;
					mat = 0;
				}
				else
				{
					--i;
					--j;
					if (dp_mat[0][i][j] != INT_MAX)
					{
						dp_mat[0][i][j] = INT_MAX;
						lib.add(seq_id1, seq_id2, i, j);
				//		if (seq1[i]==seq2[j])
				//			++match;
				//		++total;
					}
				}
			}
			else if ((state==2) || (state==6))
			{
				trace_mat[mat][i][j]-=2;
				if (state !=2)
					to_do.push(TracePosition(mat, i, j));
				if (mat==0)
					mat=1;
				else
					--i;
			}
			else
			{
				if (mat==0)
					mat=2;
				else
					--j;
			}

		}
	}

//	if (total >= 0)
//		return 100*match/total;
//	else
//		return 0;
}




void
gotoh_chaining(size_t l_seq1, size_t l_seq2, const std::map<Match, int> &match_points, Gotoh_Int_array &dp_mat, Gotoh_Trace_array &trace_mat, int gapopen, int gapext)
{
	size_t lseq1=l_seq1+1;
	size_t lseq2=l_seq2+1;
	dp_mat.resize(extents[3][lseq1][lseq2]);
	trace_mat.resize(extents[3][lseq1][lseq2]);

	// initialization
	// 0=match, 1=insert, 2=deletion
	dp_mat[0][0][0]=0;
	dp_mat[1][0][0]=gapopen;
	dp_mat[2][0][0]=gapopen;
	size_t i,j;
	for (i=1; i<lseq1; ++i)
	{
		dp_mat[1][i][0]=INT_MIN;
		dp_mat[0][i][0]=dp_mat[2][i][0]=gapopen;
	}
	for (j=1; j<lseq2; ++j)
	{
		dp_mat[0][0][j]=dp_mat[1][0][j]=gapopen;
		dp_mat[2][0][j]=INT_MIN;
	}

	// fill matrix
	int score;
	std::map<Match, int>::const_iterator it, it_end=match_points.end();

	for (i=1; i<lseq1; ++i)
	{
		for (j=1; j<lseq2; ++j)
		{
			//calculate insert value
			if (dp_mat[1][i-1][j] >= dp_mat[0][i-1][j] +gapopen)
			{
				trace_mat[1][i][j] = 'i';
				dp_mat[1][i][j] = dp_mat[1][i-1][j];
			}
			else
			{
				trace_mat[1][i][j] = 'm';
				dp_mat[1][i][j] = dp_mat[0][i-1][j]+gapopen;
			}
			dp_mat[1][i][j] += gapext;

			//calculate deletion value
			if (dp_mat[2][i][j-1] >= dp_mat[0][i][j-1] +gapopen)
			{
				trace_mat[2][i][j] = 'd';
				dp_mat[2][i][j] = dp_mat[2][i][j-1];
			}
			else
			{
				trace_mat[2][i][j] = 'm';
				dp_mat[2][i][j] = dp_mat[0][i][j-1]+gapopen;
			}
			dp_mat[2][i][j] += gapext;

			//calculate match value
			if (dp_mat[1][i][j] > dp_mat[2][i][j])
			{
				trace_mat[0][i][j] = 'i';
				dp_mat[0][i][j] = dp_mat[1][i][j];
			}
			else
			{
				trace_mat[0][i][j] = 'd';
				dp_mat[0][i][j] = dp_mat[2][i][j];
			}

			// if match has been seen give appropiate score else ignore match
			it = match_points.find(Match(i-1,j-1));
			if (it!=it_end)
			{
//				printf("ARG %li %li %i\n", i-1, j-1, it->second);
				score=it->second;
				if (dp_mat[0][i-1][j-1] + score >= dp_mat[0][i][j])
				{
					trace_mat[0][i][j] = 'm';
					dp_mat[0][i][j] = dp_mat[0][i-1][j-1] + score;
				}
			}
		}
	}
}



void
gotoh_chain_trace(std::vector<std::pair<unsigned int, unsigned int> > &aln_gaps1, std::vector<std::pair<unsigned int, unsigned int> > &aln_gaps2, const Gotoh_Trace_array &trace_mat)
{
	aln_gaps1.clear();
	aln_gaps2.clear();
	size_t i = trace_mat[0].size()-1;
	size_t j = trace_mat[0][0].size()-1;
	char state=trace_mat[0][i][j];
	int mat = 0;
	int gap_counter=0;
	while ((i!=0) && (j!=0))
	{
		if (state != trace_mat[mat][i][j])
		{
			if (state=='i')
				aln_gaps2.push_back(pair<unsigned int, unsigned int>(j,gap_counter));
			if (state=='d')
				aln_gaps1.push_back(pair<unsigned int, unsigned int>(i,gap_counter));
			gap_counter=0;
		}
		state = trace_mat[mat][i][j];
		if (state=='m')
		{
			if (mat!=0)
			{
				if (mat==1)
					--i;
				else
					--j;
				mat = 0;
			}
			else
			{
				--i;
				--j;
			}
		}
		else if (state=='i')
		{
			if (mat==0)
				mat=1;
			else
				--i;
			++gap_counter;
		}
		else if (state=='d')
		{
			if (mat==0)
				mat=2;
			else
				--j;
			++gap_counter;
		}
	}


	if (state != trace_mat[mat][i][j])
	{
		if (state=='i')
			aln_gaps2.push_back(pair<unsigned int, unsigned int>(j,gap_counter-1));
		if (state=='d')
			aln_gaps1.push_back(pair<unsigned int, unsigned int>(i,gap_counter-1));
		gap_counter=0;
	}

	if (i!=0)
		aln_gaps2.push_back(pair<unsigned int, unsigned int>(0, i));
	if (j!=0)
		aln_gaps1.push_back(pair<unsigned int, unsigned int>(0, j));
	sort(aln_gaps1.begin(), aln_gaps1.end());
	sort(aln_gaps2.begin(), aln_gaps2.end());
}


void
all_gotoh_pairs(const SequenceSet &set, BioTools::Utils::Library &lib, const BioTools::Utils::Scoring_Matrix &score_mat, int gapopen, int gapext)
{
	size_t max_len = set.max_size()+1;
	Gotoh_Int_array dp_mat;
	Gotoh_Trace_array trace_mat;
	dp_mat.resize(boost::extents[3][max_len][max_len]);
	trace_mat.resize(boost::extents[3][max_len][max_len]);
	size_t i,j;
	size_t n_seqs = set.n_seqs();
	for (i=0; i< n_seqs; ++i)
	{
		const Sequence &seq1 = set[i];
		for (j=i+1; j<n_seqs; ++j)
		{
			const Sequence &seq2 = set[j];
			vector<pair<unsigned int, unsigned int> > gap1, gap2;
			Library lib2(set);
			std::map<Match, int> match_points;
			gotoh_all_optimal(seq1, seq2, dp_mat, trace_mat, score_mat, gapopen, gapext);
			gotoh_trace_all_optimal(seq1, seq2, dp_mat, trace_mat, lib);
		}
	}
}


} // Seq
} // BioTools


