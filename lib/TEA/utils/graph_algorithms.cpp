/*
 * max_clique.cpp
 *
 *  Created on: Jun 5, 2012
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

#include "graph_algorithms.h"

using namespace std;


void
greedy_coloring(const vector<size_t> &vertices, const BioTools::Utils::Matrix<bool> &edges, set<pair<size_t, size_t> > &colored_vertices)
{
	colored_vertices.clear();
	colored_vertices.insert(make_pair<size_t, size_t> (1,vertices[0]));
	size_t j, n_vertices=vertices.size();
	vector<size_t> vals(n_vertices,1);
	vector<size_t> used;
	vector<size_t>::iterator it, it_end;
	size_t val;

	for (size_t i=1; i<n_vertices; ++i)
	{
		used.clear();
		for (j=0; j<i; ++j)
			if (edges[vertices[i]][vertices[j]])
				used.push_back(vals[j]);
		it_end=used.end();
		val=0;
		sort(used.begin(), used.end());
		for (it=used.begin(); it!=it_end; ++it)
		{
			if (*it-val>1)
				break;
			else
				val=*it;
		}
		vals[i]=++val;
		colored_vertices.insert(pair<size_t, size_t> (val,vertices[i]));
	}

}

/*
float
calc_score(set<size_t> vertices, const BioTools::Utils::Matrix<float> &edges)
{
	float score =0;
	set<size_t>::iterator it1, it2, it_end=vertices.end();
	for (it1=vertices.begin(); it1!=it_end; ++it1)
	{
		for (it2=it1; it2!=it_end; ++it2)
			score += edges[*it1][*it2];
	}
	return score;
}
*/


void
expand(set<pair<size_t, size_t> > r, set<pair<size_t, size_t>  > &q, set<pair<size_t, size_t>  > &q_max,  const BioTools::Utils::Matrix<bool> &edges)
{
	set<pair<size_t, size_t>  >r_p;
	pair<size_t, size_t> p;
	set<pair<size_t, size_t> >::iterator it_end, it;
	vector<size_t> vertices;
	while (!r.empty())
	{
		p = *r.rbegin();
		if (p.first + q.size() > q_max.size())
		{
			q.insert(p);
			it_end=--r.end();
			vertices.clear();
			for (it=r.begin(); it != it_end; ++it)
				if (edges[p.second][it->second])
					r_p.insert(*it);
			if (!r_p.empty())
				expand(r_p, q, q_max, edges);
			else
				if (q.size()>q_max.size())
					q_max=q;
			q.erase(p);
		}
		else
			return;
		r.erase(p);
	}
}


//#include <cstdio>
void
max_clique(const BioTools::Utils::Matrix<bool> &edges, set<size_t> &result)
{
	size_t n_vertices=edges.dim1();
	set<pair<size_t, size_t>  > q, colors, q_max;
	vector<pair<size_t, size_t> >vertices;
	vertices.reserve(n_vertices);
	size_t i,j;
	for (i=0; i<n_vertices; ++i)
		vertices.push_back(make_pair<size_t,size_t>(0,i));
	for (i=0; i<n_vertices; ++i)
	{
		for (j=0; j<i; ++j)
		{
			if (edges[i][j])
			{
				++vertices[i].first;
				++vertices[j].first;
			}
		}
	}
	std::sort(vertices.begin(), vertices.end());
	vector<size_t> tmp;
	tmp.reserve(n_vertices);
	vector<pair<size_t, size_t> >::reverse_iterator degree_it, degree_it_end=vertices.rend();
	for (degree_it=vertices.rbegin(); degree_it != degree_it_end; ++degree_it)
		tmp.push_back(degree_it->second);
	greedy_coloring(tmp, edges, colors);
	expand(colors, q, q_max, edges);

	set<pair<size_t, size_t>  >::iterator it, it_end =q_max.end();
	for (it = q_max.begin(); it!=it_end; ++it)
		result.insert(it->second);
}



