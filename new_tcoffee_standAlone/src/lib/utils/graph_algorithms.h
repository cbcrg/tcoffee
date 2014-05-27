
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
#include <utility>

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
