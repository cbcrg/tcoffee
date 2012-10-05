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

int aha_compare (const void * a, const void * b)
{
	return ( *(int*)a - *(int*)b );
}

/*
VectorSet *
spss_init(const VectorSet *vec_set, unsigned int k, size_t start, size_t end)
{

	printf("init\n");
	VectorSet *centers = my_malloc(sizeof(VectorSet));
	centers->dim=vec_set->dim;
	centers->n_vecs=k;
	size_t dim = vec_set->dim;
	centers->vecs=my_malloc(k*sizeof(Vector*));
	size_t i,j;
	for (i = 0; i <k; ++i)
	{
		centers->vecs[i]=my_malloc(sizeof(Vector));
		centers->vecs[i]->data=my_calloc(dim,sizeof(double));
	}


	// Calculate sum of distances from each point to each other
	double *sum_dists=my_calloc(end-start, sizeof(double));
	double dist;
	for (i=start; i<end; ++i)
	{
		for (j=i+1; j<end; ++j)
		{
			dist = km_sq_dist(vec_set->vecs[i], vec_set->vecs[j], dim);
			sum_dists[i-start] += dist;
			sum_dists[j-start] += dist;
		}
	}

	// Find minimum to every other and set as center
	double min_dist = DBL_MAX;
	size_t n_vecs = end-start, index=0;
	for (i=0; i<n_vecs; ++i)
	{
		if (sum_dists[i]<min_dist)
		{
			index=i;
			min_dist=sum_dists[i];
		}
	}

	index+=start;

	double *p, *vec;
	p = centers->vecs[0]->data;
	vec = vec_set->vecs[index]->data;
	for (j=0; j<dim; ++j)
		p[j] += vec[j];


	//calc y
	for (i=start; i<end; ++i)
		sum_dists[i-start] = km_sq_dist(vec_set->vecs[i], vec_set->vecs[index], dim);
	qsort (sum_dists, end-start, sizeof(double), aha_compare);
	int n_points = ceil(n_vecs/k);
	double y=0;
	for (i=0; i<n_points; ++i)
		y+=sum_dists[i];


	double tmp_dist;
	size_t l, min_index;
	for (i=1; i<k; ++i)
	{
		index = 0;
		dist=0;
		min_index=start;
		while (dist < y)
		{
			min_dist=DBL_MAX;

			for (l=0; l<i; ++l)
			{
				tmp_dist = km_sq_dist(centers->vecs[l], vec_set->vecs[min_index], dim);
				if (tmp_dist < min_dist)
					min_dist = tmp_dist;
			}
			++min_index;
			printf("ARG: %li\n", min_index);
			dist+=min_dist*min_dist;
		}
		printf("%li\n", min_index);
		p = centers->vecs[i]->data;
		vec = vec_set->vecs[min_index]->data;
		for (j=0; j<dim; ++j)
			p[j] += vec[j];
	}*/
//
// 	/*
// 5. For each point xi, set D (xi) to be the distance
// between xi and the nearest point in C
// 6. Find y as the sum of distances of first n/k nearest
// points from the Index
// 7. Find the unique integer i so that
// 8. D(x1)2+D(x2)2+...+D(xi)2> = y>D(x1)2+D(x2)2+...+
// D(x(i-1))2
// 9. Add xi to C
// 10. Repeat steps 5-8 until k centers*/
//
// }

typedef struct
{
	size_t x;
	double y;
} init_pair;

int init_compare (const void * a, const void * b)
{
	return ( ((init_pair*)a)->y - ((init_pair*)b)->y );
}



VectorSet *
spss_like_init(const VectorSet *vec_set, unsigned int k, size_t start, size_t end)
{
	VectorSet *centers = my_malloc(sizeof(VectorSet));
	centers->dim=vec_set->dim;
	centers->n_vecs=k;
	size_t dim = vec_set->dim;
	centers->vecs=my_malloc(k*sizeof(Vector*));
	size_t i,j;
	for (i = 0; i <k; ++i)
	{
		centers->vecs[i]=my_malloc(sizeof(Vector));
		centers->vecs[i]->data=my_calloc(dim,sizeof(double));
	}


	// Calculate sum of distances from each point to each other
	init_pair *dists=my_malloc((end-start)* sizeof(init_pair));
	init_pair *keep = dists;
	int points=end-start;
	for (i=0; i<points; ++i)
	{
		dists[i].x=start+i;
		dists[i].y=0;
	}
	double dist;

	for (i=0; i<points; ++i)
	{
		for (j=i+1; j<points; ++j)
		{
			dist = km_sq_dist(vec_set->vecs[i], vec_set->vecs[j], dim);
			dists[i].y += dist;
			dists[j].y += dist;
		}
	}

	// Find minimum to every other and set as center
	double min_dist = DBL_MAX;
	size_t n_vecs = end-start, index=0;
	for (i=0; i<n_vecs; ++i)
	{
		if (dists[i].y<min_dist)
		{
			index=dists[i].x;
			min_dist=dists[i].y;
		}
	}

// 	index+=start;

	double *p, *vec;
	p = centers->vecs[0]->data;
	vec = vec_set->vecs[index]->data;
	for (j=0; j<dim; ++j)
		p[j] += vec[j];

	for (i=start; i<end; ++i)
		dists[i-start].y = DBL_MAX;
	//calc y
	double tmp_dist;
	size_t l, min_index;
	int n_points = ceil(n_vecs/k);

	for (i=1; i<k; ++i)
	{
		for (j=0; j<points; ++j)
		{
			tmp_dist = km_sq_dist(centers->vecs[i-1], vec_set->vecs[dists[j].x], dim);
			if (tmp_dist < dists[j].y)
				dists[j].y = tmp_dist;
		}
		qsort (dists, points, sizeof(init_pair), init_compare);

		min_index=dists[n_points].x;
// 		printf("T: %li\n", min_index);
		p = centers->vecs[i]->data;
		vec = vec_set->vecs[min_index]->data;
		for (j=0; j<dim; ++j)
			p[j] += vec[j];

		dists=&dists[n_points];
		points-=n_points;
	}

	free(keep);
	return centers;
}


/*
Initialization as proposed int:
	I. Katsavounidis, C.-C. J. Kuo, and Z. Zhang. A new initialization technique for generalized Lloyd iteration. IEEE Signal Processing Letters, 1(10):144?146, 1994*/
VectorSet *
kkz_init(const VectorSet *vec_set, unsigned int k, size_t start, size_t end)
{
	VectorSet *centers = my_malloc(sizeof(VectorSet));
	centers->dim=vec_set->dim;
	centers->n_vecs=k;
	size_t dim = vec_set->dim;
	centers->vecs=my_malloc(k*sizeof(Vector*));
	size_t i,j;
	for (i = 0; i <k; ++i)
	{
		centers->vecs[i]=my_malloc(sizeof(Vector));
		centers->vecs[i]->data=my_calloc(dim,sizeof(double));
		centers->vecs[i]->id=i;
	}

	//get vector with largest l2norm
	double max_norm =-1, tmp;
	size_t index =0;
	for (i=start; i<end; ++i)
	{
		tmp=l2norm(vec_set->vecs[i], dim);
		if (tmp>max_norm)
		{
			max_norm=tmp;
			index=i;
		}
	}

	double *p, *vec;
	p = centers->vecs[0]->data;
	vec = vec_set->vecs[index]->data;
	for (j=0; j<dim; ++j)
		p[j] += vec[j];

	//calculate most distant vectors to existing centers and set as new vector. distance is the distance to the closest centroid.
	double max_dist, tmp_dist,min_dist;
	size_t l, min_index;
	for (i=1; i<k; ++i)
	{
		max_dist = -1;
		index = 0;
		for (j=start; j<end; ++j)
		{
			min_dist=DBL_MAX;
			//TODO Loop can be saved when saving results from previous round!
			for (l=0; l<i; ++l)
			{
				tmp_dist = km_sq_dist(centers->vecs[l], vec_set->vecs[j], dim);
				if (tmp_dist < min_dist)
					min_dist = tmp_dist;
			}
			if (min_dist > max_dist)
			{
				max_dist=min_dist;
				index=j;
			}
		}
// 		printf("INDEX: %i\n", index);


		p = centers->vecs[i]->data;
		vec = vec_set->vecs[index]->data;
		for (j=0; j<dim; ++j)
			p[j] += vec[j];
	}

	return centers;
}


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
	return (*(Vector**)i)->assignment - (*(Vector**)j)->assignment;
}

void
delKM_node(KM_node *node)
{
	free(node->children);
	free(node);
}

KM_node*
hierarchical_kmeans(VectorSet *vecs, unsigned int k, unsigned int k_leaf, const char *init, double error_threshold)
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
	int use_k;
	size_t start, end, i, old_index, old_assignment;
	while (todo->size != 0)
	{
		current = todo->last;
		pop(todo);
		start = current->start;
		end = current->end;
		use_k = k;//(((end-start)/k) >k) ? k : max(2,(end-start)/k);
		kmeans_sub(vecs, use_k, init, error_threshold, start, end);
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
				if (i - old_index > k_leaf)//(k*1.5))
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
		{
			printf("size\n", end-start);
			push(todo, tmp);
		}
	}

	delStack(todo);
	return root;
}


KM_node*
hierarchical_kmeans2(VectorSet *vecs, unsigned int k, unsigned int k_leaf, const char *init, double error_threshold)
{
// 	Stack *todo =Stack_init();
	size_t node_id = 0;
	int use_k;
	size_t n_vecs, i, old_assignment;
	VectorSet *centers;
	KM_node **node_array=malloc(vecs->n_vecs*sizeof(KM_node*));

	// each sequence gets a node
	for (i=0; i<vecs->n_vecs; ++i)
	{
		node_array[i]=my_malloc(sizeof(KM_node));
		node_array[i]->start=i;
		node_array[i]->end=i+1;
		node_array[i]->id = ++node_id;
		node_array[i]->n_children=0;
	}

	// cluster the sequences so that each cluster on average has k sub nodes
	int tmp_node_id;
	KM_node *tmp;
	size_t max_reserve, old_id;
	while (vecs->n_vecs > k)
	{
		n_vecs= vecs->n_vecs;
		use_k = max(2, n_vecs/k);
		centers=kmeans_sub2(vecs, use_k, init, error_threshold, 0, n_vecs);
		qsort(&vecs->vecs[0], n_vecs, sizeof(Vector*), my_assignment_sort);
		old_assignment=-1;
		tmp_node_id=-1;
		old_id=0;
		KM_node **tmp_nodes=malloc(use_k*sizeof(KM_node*));
		VectorSet *tmpSet = my_malloc(sizeof(VectorSet));
		tmpSet->dim=vecs->dim;
		tmpSet->n_vecs=0;
		tmpSet->vecs=my_malloc(use_k*sizeof(Vector*));
		for (i=0; i<vecs->n_vecs; ++i)
		{
			if (vecs->vecs[i]->assignment != old_assignment)
			{
// 				printf("%li\n", i - old_id);
				old_id=i;
				tmp=my_malloc(sizeof(KM_node));
				tmp_nodes[++tmp_node_id]=tmp;
				tmpSet->vecs[tmp_node_id]=centers->vecs[vecs->vecs[i]->assignment];
				tmpSet->vecs[tmp_node_id]->id=tmp_node_id;
				++tmpSet->n_vecs;

				max_reserve=k;
				tmp->children=malloc(max_reserve*sizeof(KM_node*));
				tmp->n_children=0;
				tmp->children[tmp->n_children++]=node_array[vecs->vecs[i]->id];
				tmp->id = ++node_id;
				old_assignment=vecs->vecs[i]->assignment;
			}
			else
			{
				if (tmp->n_children == max_reserve)
				{
					max_reserve += 5;
					tmp->children=realloc(tmp->children, max_reserve*sizeof(KM_node*));
				}
				tmp->children[tmp->n_children++]=node_array[vecs->vecs[i]->id];
			}
		}

		for (i=0; i<=tmp_node_id; ++i)
			node_array[i] = tmp_nodes[i];
		free(tmp_nodes);
		vecs=tmpSet;
		tmpSet->n_vecs=++tmp_node_id;
	}

	// add the root node
	KM_node *root;
	if (vecs->n_vecs==1)
		root=node_array[0];
	else
	{
		root=my_malloc(sizeof(KM_node));
		root->children=malloc(use_k*sizeof(KM_node*));
		root->n_children=0;
		for (i=0; i<vecs->n_vecs; ++i)
			root->children[root->n_children++]=node_array[vecs->vecs[i]->id];
	}
	free(node_array);
	root->id=0;
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
	else if (!strcmp(init, "kkz"))
		centers = kkz_init(vecs, k, start, end);
	else if (!strcmp(init, "spss"))
		centers = spss_like_init(vecs, k, start, end);
	else if (!strcmp(init, "first"))
		centers = first_init(vecs, k, start);
	else
	{
		printf("%s\n", init);
		printf("Unknown initialization value!");
		exit(1);
	}
	size_t *n_inside = malloc(k*sizeof(size_t));
	size_t i;
	for (i=0; i<k; ++i)
		n_inside[i]=1;

	size_t id, center_id, new_center_id;
	double error = DBL_MAX, old_error;
	double tmp_dist, min_dist;
	size_t dimension = vecs->dim;
	size_t j;
	Vector *center_tmp;

	int chunk = 30;
	int round=0;

	do
	{
		++round;
		old_error = error;
		error=0;

		#pragma omp parallel shared(start, end, vecs, chunk, k, centers, n_inside) private(id, new_center_id, min_dist, center_id, tmp_dist)
		{
			#pragma omp for schedule(dynamic,chunk) reduction(+:error) nowait
			for (id = start; id < end; ++id)
			{
				min_dist = DBL_MAX;
				new_center_id=-1;
				for (center_id = 0; center_id < k; ++center_id)
				{
					tmp_dist = sqrt(km_sq_dist(vecs->vecs[id], centers->vecs[center_id],dimension));
// 					tmp_dist = km_angle_dist(vecs->vecs[id]->data, centers->vecs[center_id]->data,dimension);
// 					tmp_dist = km_common(vecs->vecs[id]->data, centers->vecs[center_id]->data,dimension);
					if (tmp_dist < min_dist)
					{
						new_center_id = center_id;
						min_dist = tmp_dist;
					}
				}

				++n_inside[new_center_id];
				if ((round != 1) && (vecs->vecs[id]->assignment != new_center_id))
					--n_inside[vecs->vecs[id]->assignment];
				vecs->vecs[id]->assignment=new_center_id;
				error += min_dist;
			}
		}

		chunk=10;
		// set centers to 0
		double *tmp,*tmp2;
		for (i=0; i<k; ++i)
			n_inside[i]=0;
		#pragma omp parallel shared(k, dimension) private(tmp, center_id, center_tmp, j)
		{
			#pragma omp for schedule(dynamic, chunk) nowait
			for (center_id = 0; center_id < k; ++center_id)
			{
				tmp=centers->vecs[center_id]->data;
				for (j=0; j<dimension; ++j)
					tmp[j]=0;
			}
		}

		// calculate new centers
		for (id = start; id < end; ++id)
		{
			tmp = vecs->vecs[id]->data;
			tmp2 = centers->vecs[vecs->vecs[id]->assignment]->data;
			++n_inside[vecs->vecs[id]->assignment];
			for (j=0; j<dimension; ++j)
				tmp2[j] += tmp[j];
		}

		#pragma omp parallel shared(k, dimension, n_inside) private(center_id, center_tmp, j)
		{
			#pragma omp for schedule(dynamic, chunk) nowait
			for (center_id = 0; center_id < k; ++center_id)
			{
				center_tmp = centers->vecs[center_id];
				if (n_inside[center_id])
				{
					for (j=0; j<dimension; ++j)
						center_tmp->data[j] /= n_inside[center_id];
				}
			}
		}

// 		for (i=0; i<k; ++i)
// 			printf("%li ", n_inside[i]);
// 		printf("\n");
	} while (abs(old_error-error) > error_threshold);
	free(n_inside);
	delVecSet(centers);

}


VectorSet *
kmeans_sub2(const VectorSet *vecs, unsigned int k, const char *init, double error_threshold, size_t start, size_t end)
{
	// 	printf("%li-%li\n", start, end);
	// determine first centers
	VectorSet *centers;
	// 	if (init == "++")
	// 		plusplus_init(vecs, k, centers, start, end);
	if (!strcmp(init, "distributed"))
		centers = distributed_init(vecs, k, start, end);
	else if (!strcmp(init, "kkz"))
		centers = kkz_init(vecs, k, start, end);
		else if (!strcmp(init, "spss"))
			centers = spss_like_init(vecs, k, start, end);
	else if (!strcmp(init, "first"))
		centers = first_init(vecs, k, start);
	else
	{
		printf("%s\n", init);
		printf("Unknown initialization value!");
		exit(1);
	}
	size_t *n_inside = malloc(k*sizeof(double));
	size_t i;
	for (i=0; i<k; ++i)
		n_inside[i]=1;

	size_t id, center_id, new_center_id;
	double error = DBL_MAX, old_error;
	double tmp_dist, min_dist;
	size_t dimension = vecs->dim;
	size_t j;
	Vector *center_tmp;
	size_t *nums = my_malloc(k*sizeof(size_t));

	int chunk = 30;
	int round=0;
	do
	{
		++round;
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
					tmp_dist = sqrt(km_sq_dist(vecs->vecs[id], centers->vecs[center_id],dimension))*(n_inside[center_id]*0.2)-log(n_inside[center_id]);
// 					tmp_dist = km_sq_dist(vecs->vecs[id], centers->vecs[center_id],dimension);
// 					tmp_dist = km_muscle_dist(vecs->vecs[id], centers->vecs[center_id],dimension, 3);
// 					tmp_dist = km_angle_dist(vecs->vecs[id]->data, centers->vecs[center_id]->data,dimension);
// 					tmp_dist = 1-(km_common(vecs->vecs[id]->data, centers->vecs[center_id]->data,dimension)/vecs->vecs[id]->seq_len);
					if (tmp_dist < min_dist)
					{
						new_center_id = center_id;
						min_dist = tmp_dist;
					}
				}
				++n_inside[new_center_id];
				if ((round != 1) && (vecs->vecs[id]->assignment != new_center_id))
					--n_inside[vecs->vecs[id]->assignment];
				vecs->vecs[id]->assignment=new_center_id;
				error += min_dist;
				#pragma omp atomic
				++nums[new_center_id];
			}
		}

		chunk=10;
		// set centers to 0
		double *tmp,*tmp2;
		for (i=0; i<k; ++i)
			n_inside[i]=0;
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
			++n_inside[vecs->vecs[id]->assignment];
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
// 		for (i=0; i<k; ++i)
// 			printf("%li ", n_inside[i]);
// 		printf("\n");
	} while (old_error-error > error_threshold);
	free(nums);
//	delVecSet(centers);
	return centers;
}


double
km_common(const double *vec1, const double *vec2, size_t dim)
{
	double common = 0;
	size_t i;
	for (i=0; i<dim; ++i)
	{
		if ((vec1[i] >0) && (vec2[i] >0))
			++common;
	}
	return common;
}



