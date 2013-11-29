/*
 * kmeans.cpp
 *
 *  Created on: Jan 17, 2012
 *      Author: Carsten Kemena
 */

#include "kmeans.h"

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
kkz_init(const vector<boost::shared_ptr<Vector<double> > > &vecs, unsigned int k, vector<boost::shared_ptr<Vector<double> > > &centers, size_t start)
{
	centers.clear();
	size_t n_vecs = vecs.size();
	double tmp_norm, max_norm=-1;
	size_t id, i, j;

	// take as a first center the vector with highest norm
	for (i=0; i<n_vecs; ++i)
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
		distances[i]=0;

	double tmp_dist, max_dist;
	for (i=1; i<k; ++i)
	{
		max_dist=0;
		for(j=0; j<n_vecs;++j)
		{
			tmp_dist = sq_dist(*centers[k-1], *vecs[i]);
			if (tmp_dist>distances[i])
				distances[i]=tmp_dist;
			if (distances[i]>max_dist)
			{
				max_dist=0;
				id=i;
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


KM_node*
hierarchical_kmeans(vector<boost::shared_ptr<Vector<double> > > &vecs, unsigned int k, const string &init, double error_threshold)
{
	stack<KM_node*> todo;
	KM_node *root = new KM_node, *current, *tmp;
	root->start=0;
	root->end=vecs.size();
	root->id=0;
	todo.push(root);
	size_t node_id = 0;
	size_t start, end, i, old_index, old_assignment;
	while (!todo.empty())
	{
		current = todo.top();
		todo.pop();
		start = current->start;
		end = current->end;
		//printf("START-STOP %li %li\n", start, end);
		kmeans(vecs, k, init, error_threshold, start, end);
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
	}

	return root;
}



//TODO not working yet!
void
kmeans_elkan(const vector<boost::shared_ptr<Vector<double> > > &vecs, unsigned int k, const string &init, double error_threshold, size_t start, size_t end)
{
	// determine first centers
	vector<boost::shared_ptr<Vector<double> > > centers;
	if (init == "++")
		plusplus_init(vecs, k, centers, start, end);
	else if (init == "random")
		random_init(vecs, k, centers, start, end);
	else //(init == "first")
		first_init(vecs, k, centers, start);

	size_t n_vecs=end-start;
	size_t i,j;


	// Saving from upper and lower bounds for the speadup
	double *upper_bound = new double[n_vecs];
	bool *r = new bool[n_vecs];
	for (j=0; j<n_vecs; ++j)
		r[j]=true;
	float **lower_bound= new float*[k];
	float **center_distances = new float*[k];
	for (i=0; i<k; ++i)
	{
		center_distances[i] = new float[k];
		lower_bound[i] = new float[n_vecs];
		for (j=0; j<n_vecs; ++j)
			lower_bound[i][j]=0;
	}

	// initialize values
	for (i=0; i<k; ++i)
		for (j=i+1; j<k; ++j)
			center_distances[i][j] = center_distances[j][i] = sq_dist(*(centers[i]), *(centers[j]));
	size_t id, center_id, new_center_id;
	size_t bound_id =0;
	double tmp_dist, min_dist;
	for (id = start; id < end; ++id)
	{
		new_center_id = 0;
		min_dist=sq_dist(*(vecs[id]), *(centers[0]));
		lower_bound[0][bound_id] = min_dist;
		for (center_id = 1; center_id < k; ++center_id)
		{
			if (0.5*center_distances[new_center_id][center_id]< min_dist)
			{
				tmp_dist=sq_dist(*(vecs[id]), *(centers[center_id]));
				lower_bound[center_id][bound_id] = min_dist;
				if (tmp_dist < min_dist)
				{
					min_dist=tmp_dist;
					new_center_id = center_id;
				}
			}
		}
		upper_bound[bound_id] =min_dist;
		++bound_id;
	}


	// prepare loop
	double error = DBL_MAX, old_error;
	size_t dimension = vecs[0]->size();
	boost::shared_ptr<Vector<double> > center_tmp, vec_tmp;
	size_t *nums = new size_t[k];

	double *scores = new double[k];
	//for (i=0; i<k; ++i)



	// reassign centers until convergence
	int chunk = 30;
	do
	{
		old_error = error;
		error=0;
		for (i = 0; i<k; ++i)
		{
			nums[i] = 0;
			scores[i] = DBL_MAX;
		}
		for (i=0; i<k; ++i)
		{
			for (j=i+1; j<k; ++j)
			{
				center_distances[i][j] = center_distances[j][i] = sq_dist(*(centers[i]), *(centers[j]));
				if (center_distances[i][j] < scores[i])
					scores[i] = center_distances[i][j];
				if (center_distances[i][j] < scores[j])
					scores[j] = center_distances[i][j];
			}
			scores[i] *= 0.5;
		}
		bound_id=0;
		#pragma omp parallel shared(start, end, vecs, chunk, k, centers, nums) private(id, new_center_id, min_dist, center_id, tmp_dist)
		{
			#pragma omp for schedule(dynamic,chunk) nowait
			for (id = start; id < end; ++id)
			{
				if (upper_bound[id-start] > scores[vecs[id]->assignment()])
				{
					min_dist = DBL_MAX;
					new_center_id=0;
					for (center_id = 0; center_id < k; ++center_id)
					{
						if ((center_id != vecs[id]->assignment()) && (upper_bound[id] > lower_bound[center_id][id]) && (upper_bound[id] > 0.5*center_distances[vecs[id]->assignment()][center_id]))
						{

							if (r[start-id])
							{
								tmp_dist = sq_dist(*(vecs[id]), *(centers[center_id]));
								r[start-id]=false;
							}
							else
								tmp_dist = upper_bound[id-start];


							if (tmp_dist < min_dist)
							{
								new_center_id = center_id;
								min_dist = tmp_dist;
							}
						}
					}
					vecs[id]->assignment(new_center_id);
					#pragma omp atomic
					error += min_dist;
					#pragma omp atomic
					++nums[new_center_id];
				}
			}
		}
		for (j=0; j<n_vecs; ++j)
			r[j]=true;

		chunk=10;
		//set centers to 0
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
			for (j=0; j<dimension; ++j)
				(*center_tmp)[j] += (*vec_tmp)[j];
		}

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


		// calculate lower bounds
		for (id = start; id < end; ++id)
		{
			for (center_id = 0; center_id < k; ++center_id)
				lower_bound[start-id][center_id] = max<double>(lower_bound[start-id][center_id],0.0);
		}

	} while (old_error-error > error_threshold);

	// free memory
	delete[] nums;
	for (i=0; i<k; ++i)
	{
		delete[] lower_bound[i];
		delete[] center_distances[k];
	}
	delete[] center_distances;
	delete[] scores;
	delete[] lower_bound;
	delete[] upper_bound;

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
		// set centers to 0
		//#pragma omp parallel shared(k, dimension) private(center_id, center_tmp, j)
		//{
			//#pragma omp for schedule(dynamic, chunk) nowait
			for (center_id = 0; center_id < k; ++center_id)
			{
				center_tmp = centers[center_id];
				for (j=0; j<dimension; (*center_tmp)[j++] =0);
			}
		//}

		// calculate new centers
		for (id = start; id < end; ++id)
		{
			vec_tmp = vecs[id];
			center_tmp = centers[vec_tmp->assignment()];
			for (j=0; j<dimension; ++j)
				(*center_tmp)[j] += (*vec_tmp)[j];
		}

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
	return (i->assignment() <j->assignment());
}








}
}
