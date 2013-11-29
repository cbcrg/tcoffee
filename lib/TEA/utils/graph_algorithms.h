/*
 * max_clique.h
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

/*! \file graph_algorithms.h
    \brief This file contains some graph algorithms.
*/


#ifndef MAX_CLIQUE_H_
#define MAX_CLIQUE_H_

// C header
#include <cstdlib>

// C++ header
#include <vector>
#include <set>
#include <algorithm>

// BioTools++ header
#include "../utils/Matrix.h"


/**
 * \brief Calculates the maximum clique. Returns the first maximum clique found.
 *
 * This algorithm will return a maximum clique but if several maximum cliques the one returned will depend on the input order of the edges. It is an implementation of the MCQ algorithm described in: Etsuji Tomita, Tatsuya Akutsu, Tsutomu Matsunaga, "Efficient algorithms for finding maximum and maximal cliques: Effective tools for bioinformatics" in "Biomedical Engineering, Trends in Electronics, Communications and Software," Anthony N. Laskovski (Ed.), ISBN: 978-953-307-475-7, InTech, pp.625-640 (2011).
 * \param[in] edges The edges of the graph.
 * \param[out] result The maximum clique found.
 */
void
max_clique(const BioTools::Utils::Matrix<bool> &edges, std::set<size_t> &result);


#endif /* MAX_CLIQUE_H_ */
