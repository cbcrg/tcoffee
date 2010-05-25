/*insert here structure*/

void * vmalloc ( size_t size);
void *vcalloc ( size_t nobj, size_t size);
void *vrealloc ( void *p, size_t size);
Parameter * declare_parameter ();
void declare_bga (Parameter *PARAM);
void proto_mut ( Chromosome *, Chromosome *, Parameter *);
void declare_dos_1 ( Parameter *L_PARAM);
void declare_dos_2 ( Parameter *L_PARAM);
void reset_dos ( Parameter *PARAM);
Statistic * declare_statistic (Parameter *PARAM);
void reset_statistic (Statistic *STAT);
Population* declare_population (int tot_pop, Parameter *PARAM);
char ** declare_char ( int first, int second);
int ** declare_int ( int first, int second);
short ** declare_short ( int first, int second);
float ** declare_float ( int first,int second);
double ** declare_double ( int first,int second);
void free_int (int **array, int first);
void free_char ( char **array, int first);
void free_double (double **array, int first);
void free_float ( float **array, int first);
void free_short ( short **array, int first);
