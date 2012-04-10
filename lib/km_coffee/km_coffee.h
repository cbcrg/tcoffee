
//#include <stdlib.h>
//#include <time.h>

// #include "classes.h"
// #include "kmeans.h"
// #include "Vector.h"

// #include "Stack.h"

#ifdef _OPENMP
	#include<omp.h>
#endif


static char *km_tmp_dir = NULL;
static char *km_cwd =NULL;

typedef struct
{
	KM_node *node;
	size_t id;
}Node_pair;


int
km_coffee_align3(char *seq_f, int k, char *method, char *aln_f, int n_cores);
