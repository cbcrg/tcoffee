/*
 * Tree.cpp
 *
 *  Created on: May 14, 2012
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
#include "Tree.h"

namespace BioTools {
namespace Clustering {

using namespace std;
using namespace BioTools::Utils;



Tree::Tree():_root(NULL), _n_species(0)
{}

Tree::~Tree()
{
	clear();
}


void Tree::upgma(const Matrix<float> &dist_mat, const vector<string> &names)
{
	size_t n_taxa = dist_mat.dim1();

	Matrix<float> d_mat(n_taxa, n_taxa);
	float *r_val = new float[n_taxa];
	d_mat=dist_mat;
	size_t i,j, k;

	// Calculate nodes
	float min_val = FLT_MAX;
	TreeNode** nodes = new TreeNode*[n_taxa];
	for (i=0; i<n_taxa; ++i)
		nodes[i]=new TreeNode(i, names[i]);

	if (n_taxa==2)
	{
		_root=new TreeNode();
		_root->children.push_back(nodes[0]);
		_root->children.push_back(nodes[1]);
		_root->edge_length=0;
		_n_species=n_taxa;
		return;
	}


	TreeNode *tmpNode;
	float dist_i_new;
	size_t x_taxa = n_taxa;
	size_t min_i=-1, min_j=-1;
	size_t n_inner_nodes = n_taxa-2;

	for (i=0; i<n_taxa; ++i)
		for (j=i; j<n_taxa; ++j)
			d_mat[i][j] = d_mat[j][i] = dist_mat[i][j];

	float tmp_score;
	for (k=0; k<n_inner_nodes; ++k)
	{
		for (i=0; i<n_taxa; ++i)
			r_val[i]=0;
		for (i=0; i<n_taxa; ++i)
		{
			if (nodes[i]==NULL)
				continue;
			for (j=i+1; j<n_taxa; ++j)
			{
				if (nodes[j]==NULL)
					continue;
				r_val[i]+=d_mat[i][j];
				r_val[j]+=d_mat[i][j];
			}
			r_val[i]*=(1.0/(x_taxa-2));
		}
		min_val = FLT_MAX;

		// Find minimum value
		for (i=0; i<n_taxa; ++i)
		{
			if (nodes[i]==NULL)
				continue;
			for (j=i+1; j<n_taxa; ++j)
			{
				if (nodes[j]==NULL)
					continue;
				tmp_score = d_mat[i][j] - (r_val[i]+r_val[j]);
				if (tmp_score < min_val)
				{
					min_val=tmp_score;
					min_i = i;
					min_j = j;
				}
			}
		}

		// calc unrooted distances for last three nodes
		dist_i_new = (d_mat[min_i][min_j] + r_val[min_i]-r_val[min_j])/2;
		nodes[min_i]->edge_length=dist_i_new;
		nodes[min_j]->edge_length=d_mat[min_i][min_j]-dist_i_new;
		tmpNode = new TreeNode();
		tmpNode->children.push_back(nodes[min_i]);
		tmpNode->children.push_back(nodes[min_j]);
		nodes[min_j]=NULL;
		nodes[min_i]=tmpNode;
		for (i=0; i<n_taxa; ++i)
		{
			if ((i!=min_i) && (i!=min_j) && (nodes[i]!=NULL))
				d_mat[min_i][i]=d_mat[i][min_i]=(d_mat[min_i][i] + d_mat[min_j][i] - d_mat[min_i][min_j])*0.5;
		}
		--x_taxa;
	}

	// Merge the last three taxa
	size_t remaining[3];
	j=0;
	for (i=0; i<n_taxa; ++i)
		if (nodes[i] !=0)
			remaining[j++] = i;

	float e1, e2, e3;
	e1 = nodes[remaining[0]]->edge_length=abs(d_mat[remaining[0]][remaining[1]]-d_mat[remaining[0]][remaining[2]]);
	e2 = nodes[remaining[1]]->edge_length=d_mat[remaining[0]][remaining[1]]-nodes[remaining[0]]->edge_length;
	e3 = nodes[remaining[2]]->edge_length=d_mat[remaining[0]][remaining[2]]-nodes[remaining[0]]->edge_length;

	tmpNode = new TreeNode();
	tmpNode->edge_length=0;
	short last;
	// put together the two first which are the most similar.
	if ((e1>e2) && (e1>e3))
	{
		tmpNode->children.push_back(nodes[remaining[1]]);
		tmpNode->children.push_back(nodes[remaining[2]]);
		last=0;
	}
	else if ((e2>e1) && (e2>e3))
	{
		tmpNode->children.push_back(nodes[remaining[0]]);
		tmpNode->children.push_back(nodes[remaining[2]]);
		last=1;
	}
	else
	{
		tmpNode->children.push_back(nodes[remaining[0]]);
		tmpNode->children.push_back(nodes[remaining[1]]);
		last=2;
	}

	_root=new TreeNode();
	_root->children.push_back(tmpNode);
	_root->children.push_back(nodes[remaining[last]]);
	_root->edge_length=0;
	_n_species=n_taxa;


	delete[] nodes;
	delete[] r_val;


}


void Tree::nj(const Matrix<float> &dist_mat, const vector<string> &names)
{

	size_t n_taxa = dist_mat.dim1();

	Matrix<float> d_mat(n_taxa, n_taxa);
	float *r_val = new float[n_taxa];
	d_mat=dist_mat;
	size_t i,j, k;

	// Calculate nodes
	float min_val = FLT_MAX;
	TreeNode** nodes = new TreeNode*[n_taxa];
	for (i=0; i<n_taxa; ++i)
		nodes[i]=new TreeNode(i, names[i]);

	if (n_taxa==2)
	{
		_root=new TreeNode();
		_root->children.push_back(nodes[0]);
		_root->children.push_back(nodes[1]);
		_root->edge_length=0;
		_n_species=n_taxa;
		return;
	}


	TreeNode *tmpNode;
	float dist_i_new;
	size_t x_taxa = n_taxa;
	size_t min_i=-1, min_j=-1;
	size_t n_inner_nodes = n_taxa-3;


	float tmp_score;
	for (k=0; k<n_inner_nodes; ++k)
	{
		for (i=0; i<n_taxa; ++i)
			r_val[i]=0;
		for (i=0; i<n_taxa; ++i)
		{
			if (nodes[i]==NULL)
				continue;
			for (j=i+1; j<n_taxa; ++j)
			{
				if (nodes[j]==NULL)
					continue;
				r_val[i]+=d_mat[i][j];
				r_val[j]+=d_mat[i][j];
			}
			r_val[i]*=(1.0/(x_taxa-2));
		}
		min_val = FLT_MAX;

		// Find minimum value
		for (i=0; i<n_taxa; ++i)
		{
			if (nodes[i]==NULL)
				continue;
			for (j=i+1; j<n_taxa; ++j)
			{
				if (nodes[j]==NULL)
					continue;
				tmp_score = d_mat[i][j] - (r_val[i]+r_val[j]);
				if (tmp_score < min_val)
				{
					min_val=tmp_score;
					min_i = i;
					min_j = j;
				}
			}
		}

		// calc unrooted distances for last three nodes
		dist_i_new = (d_mat[min_i][min_j] + r_val[min_i]-r_val[min_j])/2;
		nodes[min_i]->edge_length=dist_i_new;
		nodes[min_j]->edge_length=d_mat[min_i][min_j]-dist_i_new;
		tmpNode = new TreeNode();
		tmpNode->children.push_back(nodes[min_i]);
		tmpNode->children.push_back(nodes[min_j]);
		nodes[min_j]=NULL;
		nodes[min_i]=tmpNode;
		for (i=0; i<n_taxa; ++i)
		{
			if ((i!=min_i) && (i!=min_j) && (nodes[i]!=NULL))
				d_mat[min_i][i]=d_mat[i][min_i]=(d_mat[min_i][i] + d_mat[min_j][i] - d_mat[min_i][min_j])*0.5;
		}
		--x_taxa;
	}

	// Merge the last three taxa
	size_t remaining[3];
	j=0;
	for (i=0; i<n_taxa; ++i)
		if (nodes[i] !=0)
			remaining[j++] = i;

	float e1, e2, e3;
	e1 = nodes[remaining[0]]->edge_length=abs(d_mat[remaining[0]][remaining[1]]-d_mat[remaining[0]][remaining[2]]);
	e2 = nodes[remaining[1]]->edge_length=d_mat[remaining[0]][remaining[1]]-nodes[remaining[0]]->edge_length;
	e3 = nodes[remaining[2]]->edge_length=d_mat[remaining[0]][remaining[2]]-nodes[remaining[0]]->edge_length;

	tmpNode = new TreeNode();
	tmpNode->edge_length=0;
	short last;
	// put together the two first which are the most similar.
	if ((e1>e2) && (e1>e3))
	{
		tmpNode->children.push_back(nodes[remaining[1]]);
		tmpNode->children.push_back(nodes[remaining[2]]);
		last=0;
	}
	else if ((e2>e1) && (e2>e3))
	{
		tmpNode->children.push_back(nodes[remaining[0]]);
		tmpNode->children.push_back(nodes[remaining[2]]);
		last=1;
	}
	else
	{
		tmpNode->children.push_back(nodes[remaining[0]]);
		tmpNode->children.push_back(nodes[remaining[1]]);
		last=2;
	}

	_root=new TreeNode();
	_root->children.push_back(tmpNode);
	_root->children.push_back(nodes[remaining[last]]);
	_root->edge_length=0;
	_n_species=n_taxa;


	delete[] nodes;
	delete[] r_val;
}


void
Tree::match_seq_set(const BioTools::Seq::SequenceSet &set)
{
	map<string, size_t> names;
	map<string, size_t>::iterator it, it_end;
	size_t n_seqs = set.n_seqs();
	size_t i;
	for (i=0; i<n_seqs; ++i)
		names.insert(pair<string, size_t>(set[i].name(), set[i].id()));
	it_end=names.end();

	// traverse tree and add id
	stack<pair<TreeNode*,int> > to_do;
	to_do.push(make_pair(_root, 0));
	TreeNode *current;
	size_t id;
	while (!to_do.empty())
	{
		current=to_do.top().first;
		id=to_do.top().second;
		if (!current->children.empty())
		{
			if (id != current->children.size())
			{
				++to_do.top().second;
				to_do.push(make_pair(current->children[id],0));
			}
			else
				to_do.pop();
		}
		else
		{
			it = names.find(current->name);
			if (it!=it_end)
				current->id=it->second;
			//else
			//TODO: throw exception if name not found!

		}
	}
}





void
Tree::read_newick(const string &tree_f)
{
	FILE *tree_F = my_fopen(tree_f, "r");
	string tree_string = "";
	const int LINE_LENGTH = 501;
	char line[LINE_LENGTH];
	while (fgets(line, LINE_LENGTH, tree_F) != NULL)
		tree_string.append(line);
	fclose(tree_F);

	// TODO write rest of this function :)
}

void
Tree::print_newick(const string &tree_f)
{
	typedef pair<TreeNode*, unsigned int> Node;
	stack<Node > to_do;
	TreeNode *current_node;
	to_do.push(Node(_root,0));
	FILE *tree_F = my_fopen(tree_f, "w");
	unsigned int child_id;
	while (!to_do.empty())
	{
		current_node = to_do.top().first;
		child_id = to_do.top().second;
		if (!current_node->children.empty())
		{
			if (child_id == 0)
				fprintf(tree_F, "(");
			else if (current_node->children.size() == child_id)
			{
				fprintf(tree_F, "):%f", current_node->edge_length);
				to_do.pop();
				continue;
			}
			else
			{
				fprintf(tree_F, ",");
			}
			++to_do.top().second;
			to_do.push(Node(current_node->children[child_id],0));
		}
		else
		{
			fprintf(tree_F, "%s:%f", current_node->name.c_str(), current_node->edge_length);
			to_do.pop();
		}
	}
	fprintf(tree_F, ";");
	fclose(tree_F);
}

void
Tree::clear()
{
	unsigned int i;
	stack<TreeNode* > to_do;
	TreeNode *current_node;
	to_do.push(_root);
	while (!to_do.empty())
	{
		current_node = to_do.top();
		to_do.pop();
		for (i=0; i< current_node->children.size(); ++i)
		{
			if (current_node->children[i]->children.empty())
				delete current_node->children[i];
			else
				to_do.push(current_node->children[i]);
		}
		delete current_node;
	}
	_n_species=0;
}


void
make_dist_mat(const BioTools::Seq::SequenceSet &seq_set, short k, const string &alphabet, Matrix<float> &dist_mat)
{
	vector<Vec_double_ptr>* vec_set = seqset2vecs_kmer(seq_set, k, alphabet);
	size_t n_seqs = seq_set.n_seqs();
	size_t i,j;
	for (i=0; i<n_seqs; ++i)
	{
		dist_mat[i][i] = 0;
		for (j=i+1; j<n_seqs; ++j)
		{
			dist_mat[i][j] = dist_mat[j][i] = muscle_dist((*(*vec_set)[i]), (*(*vec_set)[j]));
		}
	}

}

} /* namespace Clustering */
} /* namespace BioTools */
