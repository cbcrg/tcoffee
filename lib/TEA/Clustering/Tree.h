/*
 * Tree.h
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

/*! \file Tree.h
 * \brief Tree management
 */

#ifndef TREE_H_
#define TREE_H_


// C header
#include <cfloat>
#include <cstdio>
#include <cmath>
#include <cstdlib>

// C++ header
#include <map>
#include <stack>
#include <string>
#include <vector>

// BioTools++ header
#include "../Clustering/clustering.h"
#include "../Clustering/Vector.h"
#include "../Sequence/SequenceSet.h"
#include "../utils/filesystem.h"
#include "../utils/Matrix.h"


namespace BioTools {
namespace Clustering {

struct TreeNode
{
	unsigned int id;
	std::string name;
	float edge_length;
	std::vector<TreeNode*> children;

	TreeNode():id(0), name(""), edge_length(-1)
	{}

	TreeNode(unsigned int id_):id(id_), name(""), edge_length(-1)
	{}

	TreeNode(const std::string &name_):id(0), name(name_), edge_length(-1)
	{}

	TreeNode(unsigned int id_, const std::string &name_):id(id_),name(name_), edge_length(-1)
	{}


};


/**
 * \brief This class provides several methods to read, compute and manipulate trees.
 */
class Tree {
private:
	TreeNode *_root;
	size_t _n_species;

public:

	/**@{*/
	//! \name Constructors & Destructors
	/**
	 * Simple Constructor for a tree.
	 */
	Tree();
	virtual ~Tree();
	/**@}*/


	/**@{*/
	//! \name Tree computation methods

	/**
	 * \brief Tree construction using the neighbour joining method.
	 *
	 * This method returns a rooted binary tree. As the method usually returns a unrooted tree
	 * one of the edge lengths of the tree has been set to 0. This method has been first described in:
	 * Saitou N, Nei M. "The neighbor-joining method: a new method for reconstructing phylogenetic trees." Molecular Biology and Evolution, volume 4, issue 4, pp. 406-425, July 1987.
	 * \param dist_mat A matrix giving the distances between the species.
	 * \param names The names of the species.
	 */
	void nj(const BioTools::Utils::Matrix<float> &dist_mat, const std::vector<std::string> &names);

	/**
	 * \brief Tree construction using the UPGMA method.
	 *
	 * The UPGMA (Unweighted Pair Group Method with Arithmetic Mean) method has been first described in:
	 * Sokal R and Michener C (1958). "A statistical method for evaluating systematic relationships". University of Kansas Science Bulletin 38: 1409–1438.
	 * \param dist_mat A matrix giving the distances between the species.
	 * \param names The names of the species.
	 */
	void upgma(const BioTools::Utils::Matrix<float> &dist_mat, const std::vector<std::string> &names);
	/**@}*/


	/**@{*/
	//! \name Input/Output functions for trees


	/**
	 * \brief Reads a tree in newick format from a file.
	 * \param tree_f The file to read a tree from.
	 */
	void read_newick(const std::string &tree_f);

	/**
	 * \brief Writes the tree in newick format into to a file.
	 * \param out_f The file to write the tree to.
	 */
	void print_newick(const std::string &out_f);
	/**@}*/

	/**
	 * \brief Matches the tree to a given sequence set.
	 * @param set The sequence set
	 */
	void match_seq_set(const BioTools::Seq::SequenceSet &set);


	/**
	 * \brief Deletes all nodes from the tree.
	 */
	void clear();

	/**
	 * \brief Returns the number of species in this tree.
	 * \return The number of species
	 */
	size_t n_species()
	{
		return _n_species;
	}


	/**
	 * \brief Returns pointer to the root of the tree.
	 * @return Pointer to the root.
	 */
	TreeNode* root()
	{
		return _root;
	}

	/**
	 * \brief Returns true if no species are in the node.
	 * \return True if empty else false.
	 */
	bool empty()
	{
		return (_n_species==0);
	}
};

void
make_dist_mat(const BioTools::Seq::SequenceSet &seq_set, short k, const std::string &alphabet, BioTools::Utils::Matrix<float> &dist_mat);



} /* namespace Clustering */
} /* namespace BioTools */
#endif /* TREE_H_ */
