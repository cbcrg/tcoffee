/*
 * tree_test.cpp
 *
 *  Created on: May 21, 2012
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

#include <cstdlib>

#include "../src/Tree/Tree.h"
#include "../src/utils/Matrix.h"

using namespace std;
using namespace BioTools::Utils;
using namespace BioTools::Tree;

int
main(int argc, char *argv[])
{
	Matrix<float> mat(4,4);
	vector<string> names;
	names.push_back("A");
	names.push_back("B");
	names.push_back("C");
	names.push_back("D");
	mat[0][0]=0;
	mat[0][1]=7;
	mat[0][2]=11;
	mat[0][3]=14;
	mat[1][0]=7;
	mat[1][1]=0;
	mat[1][2]=6;
	mat[1][3]=9;
	mat[2][0]=11;
	mat[2][1]=6;
	mat[2][2]=0;
	mat[2][3]=7;
	mat[3][0]=14;
	mat[3][1]=9;
	mat[3][2]=7;
	mat[3][3]=0;

	Tree tree;
	tree.nj(mat, names);

	tree.print_newick("test.dnd");
	return EXIT_SUCCESS;
}


