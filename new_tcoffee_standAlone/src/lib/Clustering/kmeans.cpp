/*
 * kmeans.cpp
 *
 *  Created on: Jan 17, 2012
 *      Author: Carsten Kemena
 */


// C++ Header
#include <vector>
#include <algorithm>
#include <iostream>
#include <cstdlib>

#include "kmeans.h"
#include "../utils/Matrix.h"
#include "../Sequence/SequenceSet.h"
#include "../utils/utils.h"

namespace BioTools {
namespace Clustering {



void
plusplus_init(const vector<boost::shared_ptr<Vector<double> > > &vecs, unsigned int k, vector<boost::shared_ptr<Vector<double> > > &centers, size_t start, size_t end)
{
	centers.clear();
	size_t n_points = vecs.size();
	size_t center = rand() % n_points;
	double *distances = new double[n_points];
	size_t i;
	for (i = 0; i < n_points; ++i)
		distances[i] = sq_dist(*(vecs[center]), *(vecs[i]));

	if (start == end)
		k=4;
	if (k==5)
		k=3;
	//TODO
	delete[] distances;
}



void
random_init(const vector<boost::shared_ptr<Vector<double> > > &vecs, unsigned int k, vector<boost::shared_ptr<Vector<double> > > &centers, size_t start, size_t end)
{
	centers.clear();
	size_t id, j, n_vecs=end-start, found = 0, dimension = vecs[0]->size();
	map<size_t, bool> taken;
	Vector<double> *p, *vec;

	while (found < k)
	{
		id = start + rand() % n_vecs;
		if (taken.count(id))
		{
			taken.insert(pair<size_t, bool>(id,true));
			p = new Vector<double>(dimension);
			centers.push_back(boost::shared_ptr<Vector<double> >(p));
			vec = &(*(vecs[id]));
			for (j=0; j<dimension; ++j)
				(*p)[j] += (*vec)[j];
		}
	}
}


void
first_init(const vector<boost::shared_ptr<Vector<double> > > &vecs, unsigned int k, vector<boost::shared_ptr<Vector<double> > > &centers, size_t start)
{
	centers.clear();
	size_t dimension = vecs[0]->size();
	size_t j,x = start+k;
	Vector<double> *p, *vec;

	for (size_t i = start; i < x; ++i)
	{
		p = new Vector<double>(dimension);
		vec = &(*(vecs[i]));
		centers.push_back(boost::shared_ptr<Vector<double> >(p));
		for (j=0; j<dimension; ++j)
			(*p)[j] += (*vec)[j];
	}
}



void
kkz_init(const vector<boost::shared_ptr<Vector<double> > > &vecs, unsigned int k, vector<boost::shared_ptr<Vector<double> > > &centers, size_t start, size_t end)
{
	centers.clear();
	size_t n_vecs = end-start;
	double tmp_norm;
	double max_norm=euclidean_norm(*vecs[start]);
	size_t i, j;
	size_t id=start;

	// take as a first center the vector with highest norm
	for (i=start+1; i<end; ++i)
	{
		if ((tmp_norm=euclidean_norm(*vecs[i]))>max_norm)
		{
			max_norm=tmp_norm;
			id=i;
		}
	}
	Vector<double> *p = new Vector<double>(*vecs[id]);
	centers.push_back(boost::shared_ptr<Vector<double> >(p));

	// add the most distant vector (from all already existing centers) as a next center
	double *distances = new double[n_vecs];
	for (i=0; i<n_vecs; ++i)
		distances[i]=DBL_MAX;

	double tmp_dist, max_dist;
	for (i=1; i<k; ++i)
	{
		max_dist=-1;
		for(j=0; j<n_vecs; ++j)
		{
			tmp_dist = sq_dist(*centers[i-1], *vecs[start+j]);
			if (tmp_dist<distances[j])
				distances[j]=tmp_dist;
			if (distances[j]>max_dist)
			{
				max_dist=distances[j];
				id=j+start;
			}
		}
		p = new Vector<double>(*vecs[id]);
		centers.push_back(boost::shared_ptr<Vector<double> >(p));
	}

	delete[] distances;
}




void
kmeans_elkan(const vector<boost::shared_ptr<Vector<double> > > &vecs, unsigned int k, const string &init, double error_threshold)
{
	kmeans_elkan(vecs, k, init, error_threshold, 0, vecs.size());
}



void
kmeans(const vector<boost::shared_ptr<Vector<double> > > &vecs, unsigned int k, const string &init, double error_threshold)
{
	kmeans(vecs, k, init, error_threshold, 0, vecs.size());
}

int
my_assignment_sort2 (const void* p1, const void* p2)
{
	return ((*(boost::shared_ptr<Vector<double> > *)p1)->assignment()-(*(boost::shared_ptr<Vector<double> > *)p2)->assignment());
}

bool
my_vecs_sort (boost::shared_ptr<Vector<double> > i, boost::shared_ptr<Vector<double> > j)
{
	if (i->seqsize() != j->seqsize())
		return (i->seqsize() <j->seqsize());
	if (i->seqdiv() != j->seqdiv())
		return (i->seqdiv() > j->seqdiv());
	/*if (i->id() != j->id())
		return (i->id() < j->id());*/
	return false;
		
}

bool
my_vecs_sort2 (boost::shared_ptr<Vector<double> > i, boost::shared_ptr<Vector<double> > j)
{
	if (i->seqdiv() == j->seqdiv())
		return (i->seqsize() > j->seqsize());
	else
		return (i->seqdiv() <j->seqdiv());
}

KM_node*
hierarchical_kmeans(vector<boost::shared_ptr<Vector<double> > > &vecs, const BioTools::Seq::SequenceSet &seq_set, unsigned int k, const string &init, double error_threshold)
{
	stack<KM_node*> todo;
	KM_node *root = new KM_node, *current, *tmp;
	root->start=0;
	root->end=vecs.size();
	root->id=0;
	todo.push(root);
	size_t node_id = 0;
	size_t start, end, i, j, old_index, old_assignment;
	
	for ( i=0; i<vecs.size();++i)
	{
	  vecs[i]->seqsize(vec_sum(*(vecs[i])));
	  vecs[i]->seqdiv(vec_div(*(vecs[i])));
	}
	
	sort(vecs.begin(), vecs.end(), my_vecs_sort);
	/*for ( i=0; i<vecs.size();++i)
	{
		    cout << vecs[i]->id() <<" " << vecs[i]->seqsize() <<" " << vecs[i]->seqdiv() << "\t";
		    const BioTools::Seq::Sequence &seq = seq_set[vecs[i]->id()];
		    size_t len=seq.size();
		    for (j=0; j<len; ++j)
		        cout<< seq[j];
		    cout<< endl;
	}*/
	
	while (!todo.empty())
	{
		current = todo.top();
		todo.pop();
		start = current->start;
		end = current->end;
		unsigned int use_k = 2;// min<unsigned int>( (log10(k)+1), k);
		kmeans(vecs, use_k, init, error_threshold, start, end);
		//qsort(vecs.data()+start, end-start, sizeof(boost::shared_ptr<Vector<double> >), my_assignment_sort2);
	
		sort(vecs.begin()+current->start, vecs.begin()+current->end, my_assignment_sort);
		
		old_index=start;
		old_assignment=vecs[start]->assignment();

		for (i=start; i<end; ++i)
		{
			if (vecs[i]->assignment() != old_assignment)
			{
				tmp=new KM_node;
				current->children.push_back(tmp);
				tmp->start=old_index;
				tmp->end=i;
				tmp->id = ++node_id;
				if (i - old_index > k)
					todo.push(tmp);
				old_index=i;
				old_assignment=vecs[i]->assignment();
			}
		}
		if (vecs[start]->assignment() != vecs[end-1]->assignment())
		{
			tmp=new KM_node;
			current->children.push_back(tmp);
			tmp->start=old_index;
			tmp->end=i;
			tmp->id=++node_id;
			if (i - old_index > k)
				todo.push(tmp);
		}
		/*for ( i=current->start; i<current->end; ++i)
		{
		    cout << vecs[i]->id() <<" " << vecs[i]->assignment() << "\t";
		    const BioTools::Seq::Sequence &seq = seq_set[vecs[i]->id()];
		    size_t len=seq.size();
		    for (j=0; j<len; ++j)
		        cout<< seq[j];
		    cout<< endl;
		}cout<<"\n\n";*/
	}
	return root;
}




void
kmeans(const vector<boost::shared_ptr<Vector<double> > > &vecs, unsigned int k, const string &init, double error_threshold, size_t start, size_t end)
{
	// determine first centers
	vector<boost::shared_ptr<Vector<double> > > centers;
	if (init == "++")
		plusplus_init(vecs, k, centers, start, end);
	else if (init == "random")
		random_init(vecs, k, centers, start, end);
	else if (init == "kkz")
		kkz_init(vecs, k, centers, start, end);
	else //(init == "first")
		first_init(vecs, k, centers, start);

	size_t id, center_id, new_center_id=0;
	double error = DBL_MAX, old_error;
	double tmp_dist, min_dist;
	size_t dimension = vecs[0]->size();
	size_t j;
	boost::shared_ptr<Vector<double> > center_tmp, vec_tmp;
	size_t *nums = new size_t[k];

	int chunk = 30;
	do
	{
		old_error = error;
		error=0;
		for (size_t i = 0; i<k; ++i)
			nums[i] = 0;
		#pragma omp parallel shared(start, end, vecs, chunk, k, centers, nums) private(id, new_center_id, min_dist, center_id, tmp_dist)
		{
			#pragma omp for schedule(dynamic,chunk) nowait
			for (id = start; id < end; ++id)
			{
				min_dist = DBL_MAX;
				new_center_id=0;
				for (center_id = 0; center_id < k; ++center_id)
				{
					tmp_dist = sq_dist(*(vecs[id]), *(centers[center_id]));
					if (tmp_dist < min_dist)
					{
						new_center_id = center_id;
						min_dist = tmp_dist;
					}
				}
				vecs[id]->assignment(new_center_id);
				#pragma omp atomic
				error += min_dist;
				#pragma omp atomic
				++nums[new_center_id];

			}
		}

		chunk=10;
		//boost::shared_ptr<Vector<double> > center_tmp, vec_tmp;
		// set centers to 0
		#pragma omp parallel shared(k, dimension) private(center_id, center_tmp, j)
		{
			#pragma omp for schedule(dynamic, chunk) nowait
			for (center_id = 0; center_id < k; ++center_id)
			{
				center_tmp = centers[center_id];
				for (j=0; j<dimension; (*center_tmp)[j++] =0);
			}
		}

		// calculate new centers
		for (id = start; id < end; ++id)
		{
			vec_tmp = vecs[id];
			center_tmp = centers[vec_tmp->assignment()];
			//++nums[vecs[id]->assignment()];  
			//cout << "|  "<< nums[vecs[id]->assignment()] << "  "<< vecs[id]->assignment() ;
			for (j=0; j<dimension; ++j)
				(*center_tmp)[j] += (*vec_tmp)[j];
		}//cout<<endl<<endl;

		#pragma omp parallel shared(k, dimension, nums) private(center_id, center_tmp, j)
		{
			#pragma omp for schedule(dynamic, chunk) nowait
			for (center_id = 0; center_id < k; ++center_id)
			{
				center_tmp = centers[center_id];
				if (nums[center_id])
				{
					for (j=0; j<dimension; ++j)
						(*center_tmp)[j] /= nums[center_id];
				}
			}
		}
	} while (old_error-error > error_threshold);
	delete[] nums;

}



double
norm_inner_dist(const Vector<double> &vec1, const Vector<double> &vec2)
{
	size_t vec_len = vec1.size();
	double inner_product=0, norm1=0, norm2=0;
	for (size_t i=0; i<vec_len; ++i)
	{
		inner_product += vec1[i]*vec2[i];
		norm1=vec1[i]*vec1[i];
		norm2=vec2[i]*vec2[i];
	}
	return (inner_product/(sqrt(norm1)*sqrt(norm2)));
}


double
norm_inner_criteria()
{
	return 0.0;
}


bool
my_assignment_sort (boost::shared_ptr<Vector<double> > i, boost::shared_ptr<Vector<double> > j)
{
	if (i->assignment() == j->assignment())
		return (i->id() < j->id());
	else
		return (i->assignment() <j->assignment());
}















//TODO not working yet!
void
kmeans_elkan(const vector<boost::shared_ptr<Vector<double> > > &vecs, unsigned int k, const string &init, double error_threshold, size_t start, size_t end)
{
	// determine first center		
}







}
}
