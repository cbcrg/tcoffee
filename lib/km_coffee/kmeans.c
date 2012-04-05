// #include "kmeans.h"

#include "km_coffee_header.h"

/*
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
}*/


/*
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
*/

VectorSet *
distributed_init(const VectorSet *vec_set, unsigned int k, size_t start, size_t end)
{
	VectorSet *centers = my_malloc(sizeof(VectorSet));
	centers->dim=vec_set->dim;
	centers->n_vecs=k;
	size_t dim = vec_set->dim;
	centers->vecs=my_malloc(k*sizeof(Vector*));
	size_t i;
	for (i = 0; i <k; ++i)
	{
		centers->vecs[i]=my_malloc(sizeof(Vector));
		centers->vecs[i]->data=my_calloc(dim,sizeof(double));
	}


	double *p, *vec;
	size_t j,pos=start;
	int step = (end-start)/k;

	for (i = 0; i < k; ++i)
	{
		p = centers->vecs[i]->data;
		vec = vec_set->vecs[pos]->data;
		for (j=0; j<dim; ++j)
			p[j] += vec[j];

		pos+=step;
	}
	return centers;
}



VectorSet *
first_init(const VectorSet *vec_set, unsigned int k, size_t start)
{
	VectorSet *centers = my_malloc(sizeof(VectorSet));
	centers->dim=vec_set->dim;
	centers->n_vecs=k;
	size_t dim = vec_set->dim;
	centers->vecs=my_malloc(k*sizeof(Vector*));
	size_t i;
	for (i = 0; i <k; ++i)
	{
		centers->vecs[i]=my_malloc(sizeof(Vector));
		centers->vecs[i]->data=my_calloc(dim,sizeof(double));
	}

	size_t j,x = start+k;
	double *p, *vec;


	for (i = start; i < x; ++i)
	{
		p = centers->vecs[i-start]->data;
		vec = vec_set->vecs[i]->data;

		for (j=0; j<dim; ++j)
			p[j] += vec[j];
	}
	return centers;
}


void
kmeans(VectorSet *vecs, unsigned int k, const char *init, double error_threshold)
{
	kmeans_sub(vecs, k, init, error_threshold, 0, vecs->n_vecs);
}


int
my_assignment_sort (const void *i, const void *j)
{
	return (*(Vector**)i)->assignment < (*(Vector**)j)->assignment;
}

void
delKM_node(KM_node *node)
{
	free(node->children);
	free(node);
}

KM_node*
hierarchical_kmeans(VectorSet *vecs, unsigned int k, const char *init, double error_threshold)
{
	Stack *todo =Stack_init();
	KM_node *root = my_malloc(sizeof(KM_node));
	root->children=malloc(k*sizeof(KM_node));
	KM_node *current, *tmp;
	root->n_children=0;
	root->start=0;
	root->end=vecs->n_vecs;
	root->id=0;
	push(todo, root);
	size_t node_id = 0;
	size_t start, end, i, old_index, old_assignment;
	while (todo->size != 0)
	{
		current = todo->last;
		pop(todo);
		start = current->start;
		end = current->end;
		kmeans_sub(vecs, k, init, error_threshold, start, end);
		// 		printf("%li %li\n", current->start, current->end-current->start);
		qsort(&vecs->vecs[current->start], current->end-current->start, sizeof(Vector*), my_assignment_sort);
		old_index=start;
		old_assignment=vecs->vecs[start]->assignment;

		for (i=start; i<end; ++i)
		{
			if (vecs->vecs[i]->assignment != old_assignment)
			{
				tmp=my_malloc(sizeof(KM_node));
				tmp->children=malloc(k*sizeof(KM_node*));
				tmp->n_children=0;
				current->children[current->n_children++]=tmp;

				tmp->start=old_index;
				tmp->end=i;
				tmp->id = ++node_id;
				if (i - old_index > k)
					push(todo, tmp);
				old_index=i;
				old_assignment=vecs->vecs[i]->assignment;
			}
		}

		tmp=my_malloc(sizeof(KM_node));
		tmp->children=malloc(k*sizeof(KM_node*));
		tmp->n_children=0;
		current->children[current->n_children++]=tmp;
		tmp->start=old_index;
		tmp->end=i;
		tmp->id=++node_id;
		if ((vecs->vecs[start]->assignment != vecs->vecs[end-1]->assignment) && (i - old_index > k))
			push(todo, tmp);

	}

	delStack(todo);
	return root;
}


void
kmeans_sub(const VectorSet *vecs, unsigned int k, const char *init, double error_threshold, size_t start, size_t end)
{
	// 	printf("%li-%li\n", start, end);
	// determine first centers
	VectorSet *centers;
	// 	if (init == "++")
	// 		plusplus_init(vecs, k, centers, start, end);
	if (!strcmp(init, "distributed"))
		centers = distributed_init(vecs, k, start, end);
	else if (!strcmp(init, "first"))
		centers = first_init(vecs, k, start);
	else
	{
		printf("Unknown initialization value!");
		exit(1);
	}
	size_t id, center_id, new_center_id;
	double error = DBL_MAX, old_error;
	double tmp_dist, min_dist;
	size_t dimension = vecs->dim;
	size_t j;
	Vector *center_tmp;
	size_t *nums = my_malloc(k*sizeof(size_t));

	int chunk = 30;
	size_t i;
	do
	{
		old_error = error;
		error=0;
		for (i = 0; i<k; ++i)
			nums[i] = 0;
		#pragma omp parallel shared(start, end, vecs, chunk, k, centers, nums) private(id, new_center_id, min_dist, center_id, tmp_dist)
		{
			#pragma omp for schedule(dynamic,chunk) reduction(+:error) nowait
			for (id = start; id < end; ++id)
			{
				min_dist = DBL_MAX;
				new_center_id=-1;
				for (center_id = 0; center_id < k; ++center_id)
				{
					tmp_dist = km_sq_dist(vecs->vecs[id]->data, centers->vecs[center_id]->data,dimension);

					if (tmp_dist < min_dist)
					{
						new_center_id = center_id;
						min_dist = tmp_dist;
					}
				}
				vecs->vecs[id]->assignment=new_center_id;
				error += min_dist;
				#pragma omp atomic
				++nums[new_center_id];
			}
		}

		chunk=10;
		// set centers to 0
		double *tmp,*tmp2;
		#pragma omp parallel shared(k, dimension) private(tmp, center_id, center_tmp, j)
		{
			#pragma omp for schedule(dynamic, chunk) nowait
			for (center_id = 0; center_id < k; ++center_id)
			{
				tmp=centers->vecs[center_id]->data;
				for (j=0; j<dimension; tmp[j++] =0);
			}
		}

		// calculate new centers
		for (id = start; id < end; ++id)
		{
			tmp = vecs->vecs[id]->data;
			tmp2 = centers->vecs[vecs->vecs[id]->assignment]->data;
			for (j=0; j<dimension; ++j)
				tmp2[j] += tmp[j];
		}

		#pragma omp parallel shared(k, dimension, nums) private(center_id, center_tmp, j)
		{
			#pragma omp for schedule(dynamic, chunk) nowait
			for (center_id = 0; center_id < k; ++center_id)
			{
				center_tmp = centers->vecs[center_id];
				if (nums[center_id])
				{
					for (j=0; j<dimension; ++j)
						center_tmp->data[j] /= nums[center_id];
				}
			}
		}
	} while (old_error-error > error_threshold);
	free(nums);
	delVecSet(centers);

}

double
km_sq_dist(const double *vec1, const double *vec2, size_t dim)
{
	double dist = 0;
	double tmp;
	size_t i;
	for (i=0; i<dim; ++i)
	{
		tmp = vec1[i]-vec2[i];
		dist += (tmp * tmp);
	}
	return dist;
}





