/*
 * kmeans.h
 *
 *  Created on: Jan 17, 2012
 *      Author: Carsten Kemena
 */

#ifndef KMEANS_H_
#define KMEANS_H_

/*! \file kmeans.h
    \brief Contains the algorithms for the kmeans clustering algorithm.
*/


// C headers
#include<cmath>
#include<cfloat>
#include<cstdio>


// C++ headers
#include<vector>
#include<map>
#include<stack>
#include<algorithm>

// Boost header
#include <boost/shared_ptr.hpp>

// BioTools header
#include "Vector.h"

using namespace std;

namespace BioTools {
namespace Clustering {




/**
 * \brief Initializes the centers by using the kmeans++ initialization.
 * \tparam double The type of data stored in the Vector.
 * \param vecs The set of vectors.
 * \param k The number of centers.
 * \param centers The vector where the centers will be saved in.
 * \param start The start of the position to consider.
 * \param end The end position to consider.
 *
 */
void
plusplus_init(const vector<boost::shared_ptr<Vector<double> > > &vecs, unsigned int k, vector<boost::shared_ptr<Vector<double> > > &centers, size_t start, size_t end);


/**
 * \brief Initializes the centers by randomly taking k vectors.
 * \param vecs The set of vectors.
 * \param k The number of centers.
 * \param centers The vector where the centers will be saved in.
 * \param start The start of the position to consider.
 * \param end The end position to consider.
 *
 */
void
random_init(const vector<boost::shared_ptr<Vector<double> > > &vecs, unsigned int k, vector<boost::shared_ptr<Vector<double> > > &centers, size_t start, size_t end);


/**
 * \brief Initializes the centers by taking the first k vectors.
 * \param vecs The set of vectors.
 * \param k The number of centers.
 * \param centers The vector where the centers will be saved in.
 * \param start The start of the position to consider.
 *
 */
void
first_init(const vector<boost::shared_ptr<Vector<double> > > &vecs, unsigned int k, vector<boost::shared_ptr<Vector<double> > > &centers, size_t start);


/**
 * \brief Initializes the centers by maximizing the differences between them.
 *
 * This method has been proposed in:
 * I. Katsavounidis, C.-C. J. Kuo, and Z. Zhang. A new initialization technique for generalized Lloyd iteration. IEEE Signal Processing Letters, 1(10):144–146, 1994.
 *
 * \param vecs The set of vectors.
 * \param k The number of centers.
 * \param centers The vector where the centers will be saved in.
 * \param start The start of the position to consider.
 *
 */
void
kkz_init(const vector<boost::shared_ptr<Vector<double> > > &vecs, unsigned int k, vector<boost::shared_ptr<Vector<double> > > &centers, size_t start);






struct KM_node {
	size_t start;
	size_t end;
	vector<KM_node*> children;
	size_t id;
};

typedef struct KM_node KM_node;


bool
my_assignment_sort (boost::shared_ptr<Vector<double> > i, boost::shared_ptr<Vector<double> > j);



//bool my_assignment_sort (boost::shared_ptr<Vector<unsigned int> > i, boost::shared_ptr<Vector<unsigned int> > j);



/**
 * \brief Repeated clustering until size threshold is met.
 *
 * \param vecs The vectors to cluster.
 * \param k the number of clusters to use.
 * \param init The methods to use for initializing the first centers: "random", "++" or "first"
 * \param error_threshold The change of error to stop kmeans
 */
KM_node*
hierarchical_kmeans(vector<boost::shared_ptr<Vector<double> > > &vecs, unsigned int k, const string &init, double error_threshold);


/**
 * \brief The kmeans clustering algorithm.
 *
 * \param vecs The vectors to cluster.
 * \param k the number of clusters to use.
 * \param init The methods to use for initializing the first centers: "random", "++" or "first"
 * \param error_threshold The change of error to stop kmeans
 */
void
kmeans(const vector<boost::shared_ptr<Vector<double> > > &vecs, unsigned int k, const string &init, double error_threshold);

/**
 * \brief The kmeans clustering algorithm.
 *
 * \param vecs The vectors to cluster.
 * \param k the number of clusters to use.
 * \param init The methods to use for initializing the first centers: "random", "++" or "first".
 * \param error_threshold The change of error to stop kmeans.
 * \param start The startint point for the clustering.
 * \param end The end point for the clustering.
 */
void
kmeans(const vector<boost::shared_ptr<Vector<double> > > &vecs, unsigned int k, const string &init, double error_threshold, size_t start, size_t end);

/**
 * \brief The kmeans clustering algorithm.
 *
 * This kmeans is an implementation of the acceleration as described in:
 * Elkan, C. (2003). "Using the triangle inequality to accelerate k-means". Proceedings of the Twentieth International Conference on Machine Learning (ICML).
 * \param vecs The vectors to cluster.
 * \param k the number of clusters to use.
 * \param init The methods to use for initializing the first centers: "random", "++" or "first"
 * \param error_threshold The change of error to stop kmeans
 */
void
kmeans_elkan(const vector<boost::shared_ptr<Vector<double> > > &vecs, unsigned int k, const string &init, double error_threshold);

/**
 * \brief The kmeans clustering algorithm.
 *
 * This kmeans is an implementation of the acceleration as described in:
 * Elkan, C. (2003). "Using the triangle inequality to accelerate k-means". Proceedings of the Twentieth International Conference on Machine Learning (ICML).
 * \param vecs The vectors to cluster.
 * \param k the number of clusters to use.
 * \param init The methods to use for initializing the first centers: "random", "++" or "first".
 * \param error_threshold The change of error to stop kmeans.
 * \param start The startint point for the clustering.
 * \param end The end point for the clustering.
 */
void
kmeans_elkan(const vector<boost::shared_ptr<Vector<double> > > &vecs, unsigned int k, const string &init, double error_threshold, size_t start, size_t end);



}
} /* namespace BioTools */
#endif /* KMEANS_H_ */
