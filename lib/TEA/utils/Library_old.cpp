/*
 * Library.cpp
 *
 *  Created on: Apr 6, 2012
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


#include "Library.h"

namespace BioTools {
namespace Utils {


using namespace std;
using namespace BioTools::Seq;

Library::Library(const SequenceSet &set):_n_seqs(set.n_seqs()),_max_length(0)
{
	//_pairs.reset(new std::vector<std::vector<std::vector<Match_point> > >);
	_pairs.resize(_n_seqs);
	for (size_t i = 0; i<_n_seqs; ++i)
	{
		_pairs[i].resize(set[i].size());
		if (_max_length<set[i].size())
			_max_length=set[i].size();
	}
	omp_init_lock(&add_lock);
	omp_init_lock(&relax_lock);
}

Library::~Library() {
	 omp_destroy_lock(&add_lock);
	 omp_destroy_lock(&relax_lock);
	// TODO Auto-generated destructor stub
}


void
Library::add(unsigned int seq1_id, unsigned int seq2_id, unsigned int pos1, unsigned int pos2, float score)
{
	omp_set_lock(&add_lock);
	#pragma omp flush (_pairs)
	Match_point match_p = Match_point(seq2_id, pos2, score);
	_pairs[seq1_id][pos1].push_back(match_p);
	match_p.x=seq1_id;
	match_p.y=pos1;
	_pairs[seq2_id][pos2].push_back(match_p);
	//#pragma omp flush (_pairs)
	omp_unset_lock(&add_lock);
}


void
Library::relax()
{

	size_t n_pairs=_n_seqs*(_n_seqs-1)/2;
	_relaxed_pairs.resize(n_pairs);

	size_t i,j;
	for (i=0; i< _n_seqs; ++i)
	{
		for (j=0; j<_pairs[i].size(); ++j)
		{
			sort(_pairs[i][j].begin(), _pairs[i][j].end());
		}

	}

	int chunk=1;
	#pragma omp parallel shared(n_pairs, chunk)
	{
		Matrix<float> hash(_n_seqs, _max_length);
		Match_point match_point;
		size_t id;
		float score;
		unsigned int i,n_residues, r1, s1,s2,r2,x,len1,len2,t_s2,t_r2,a;

		#pragma omp for schedule(dynamic,chunk) nowait
		for (s1=0; s1<_n_seqs; ++s1)
		{
			n_residues=_pairs[s1].size();
			for (r1=0; r1<n_residues; ++r1)
			{
				len1=_pairs[s1][r1].size();
				for (x=0; x< len1; ++x)
					hash[_pairs[s1][r1][x].x][_pairs[s1][r1][x].y] = _pairs[s1][r1][x].score;

				for ( a=0; a<len1; ++a)
				{
					score=0;
					s2 = _pairs[s1][r1][a].x;
					r2 = _pairs[s1][r1][a].y;
					len2=_pairs[s2][r2].size();
					for (x=0; x< len2; ++x)
					{
						t_s2 = _pairs[s2][r2][x].x;
						t_r2 = _pairs[s2][r2][x].y;
						if (t_s2==s1 && t_r2==r1)
							score += (2*_pairs[s2][r2][x].score);
						else if (hash[t_s2][t_r2])
							score += (hash[t_s2][t_r2]*_pairs[s2][r2][x].score);
					}

					if (score !=0)
					{
						score /= _n_seqs;
						score *=100;
						if (s1<s2)
						{
							id = s1*_n_seqs - (s1*(s1+1)/2) +s2-s1-1;
							match_point.x=r1;
							match_point.y=r2;
						}
						else
						{
							id = s2*_n_seqs - (s2*(s2+1)/2) +s1-s2-1;
							match_point.x=r2;
							match_point.y=r1;
						}
						match_point.score=score;
						omp_set_lock(&relax_lock);
							#pragma omp flush (_relaxed_pairs)
							_relaxed_pairs[id].push_back(match_point);
							#pragma omp flush (_relaxed_pairs)
						omp_unset_lock(&relax_lock);
					}
				}
				for (x=1; x< len1; ++x)
					hash[_pairs[s1][r1][x].x][_pairs[s1][r1][x].y]=0;
			}
		}
	}//OPEN_MP
}



void Library::get(const BioTools::Seq::Alignment &aln1, const BioTools::Seq::Alignment &aln2, std::map<Match, int> &match_points) const
{

	match_points.clear();
	size_t n_seq1 = aln1.n_seqs();
	size_t n_seq2 = aln2.n_seqs();

	size_t i,j,k,pos,id, seq1_id, seq2_id;
	map<Match, int> tmp_match_points;
	map<Match, int>::iterator it, it_end;
	vector<Match_point>::const_iterator it_lib, it_lib_end;
	vector<unsigned int> convert1, convert2;
	convert1.resize(aln1.size());
	convert2.resize(aln2.size());
	pair<map<Match, int>::iterator,bool> ret;
	pair<Match, int> tmp_match;

	for (i=0; i<n_seq1; ++i)
	{
		const Sequence &seq1 = aln1[i];
		seq1_id = seq1.id();
		pos=0;

		for (k=0; k<aln1.size(); ++k)
		{
			if (seq1[k] != '-')
				convert1[pos++] = k;
		}

		for (j=0; j<n_seq2; ++j)
		{

			const Sequence &seq2 = aln2[j];
			seq2_id = seq2.id();
			pos=0;

			for (k=0; k<aln2.size(); ++k)
			{
				if (seq2[k] != '-')
					convert2[pos++] = k;
			}

			if (seq1_id < seq2_id)
			{
				id = seq1_id*_n_seqs - (seq1_id*(seq1_id+1)/2) +seq2_id-seq1_id-1;
				it_lib_end=_relaxed_pairs[id].end();
				//it_end=tmp_match_points.end();
				for (it_lib=_relaxed_pairs[id].begin(); it_lib!=it_lib_end; ++it_lib)
				{
					tmp_match.first.first  = convert1[it_lib->x];
					tmp_match.first.second = convert2[it_lib->y];
					tmp_match.second=it_lib->score;
					ret=match_points.insert(tmp_match);
					//printf("1: %u %u %i\n", tmp_match.first.first, tmp_match.first.second, tmp_match.second);
					if (!ret.second)
						ret.first->second+=it_lib->score;
				}
			}
			else
			{
				id = seq2_id*_n_seqs - (seq2_id*(seq2_id+1)/2) +seq1_id-seq2_id-1;
				it_lib_end=_relaxed_pairs[id].end();
				for (it_lib=_relaxed_pairs[id].begin(); it_lib!=it_lib_end; ++it_lib)
				{
					tmp_match.first.first  = convert1[it_lib->y];
					tmp_match.first.second = convert2[it_lib->x];
					tmp_match.second=it_lib->score;
					//printf("1: %u %u %i\n", tmp_match.first.first, tmp_match.first.second, tmp_match.second);
					ret=match_points.insert(tmp_match);
					if (!ret.second)
						ret.first->second+=it_lib->score;
				}
			}
		}
	}
	//it_end = match_points.end();
	//size_t overall=n_seq1*n_seq2; //TODO better *???
	//printf("AHA2\n");
	//for (it = match_points.begin(); it != it_end; ++it )
	//{
//		it->second/=overall;
		//if ((it->first.first < 0) && (it->first.second < 0))
			//printf("ARG\n");
	//	printf("%li %li %i\n", it->first.first, it->first.second, it->second);
	//}
}





} /* namespace Utils */
} /* namespace BioTools */

