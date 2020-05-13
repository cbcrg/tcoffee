#ifndef __UTIL_DECLARE_H
#define __UTIL_DECLARE_H


typedef struct Tmpname Tmpname;
struct Memcontrol
    {
      size_t size;
      size_t size_element;
      char check[3];
      struct Memcontrol *p;
      struct Memcontrol *n;
    };
typedef struct Memcontrol Memcontrol;
void free_pair_wise();//Frees static memory in the pair_wise functions
/************************************************************************/
/*                                                                      */
/*            CONSTRAINT_LIST                                           */
/*                                                                      */
/*                                                                      */
/************************************************************************/
int *** duplicate_residue_index (int ***r);
int *** declare_residue_index (Sequence *S);

Constraint_list *free_constraint_list4lib_computation (Constraint_list *CL);
Constraint_list *duplicate_constraint_list4lib_computation (Constraint_list *CL);
Constraint_list * declare_constraint_list_simple ( Sequence *S);
Constraint_list * declare_constraint_list ( Sequence *S, char *name, int *L, int ne,FILE *fp, int **M);
Constraint_list *cache_dp_value4constraint_list ( char mode[],Constraint_list *CL);
Constraint_list *duplicate_constraint_list_soft (Constraint_list *CL);
Constraint_list *duplicate_constraint_list      (Constraint_list *CL);
Constraint_list *copy_constraint_list      (Constraint_list *CL, int mode);
Sequence        * free_constraint_list (Constraint_list *CL);
Constraint_list * free_constraint_list_full (Constraint_list *CL);
Distance_matrix * free_distance_matrix ( Distance_matrix *DM);
Distance_matrix * duplicate_distance_matrix ( Distance_matrix *DMin);
/************************************************************************/
/*                                                                      */
/*            Blast_param Functions                                            */
/*                                                                      */
/*                                                                      */
/************************************************************************/
Blast_param * duplicate_blast_param ( Blast_param*B);
Blast_param * free_blast_param ( Blast_param*B);
/************************************************************************/
/*                                                                      */
/*            TC_param Functions                                            */
/*                                                                      */
/*                                                                      */
/************************************************************************/
TC_param * duplicate_TC_param ( TC_param*B);
TC_param * free_TC_param ( TC_param*B);
/************************************************************************/
/*                                                                      */
/*            MOCA Functions                                            */
/*                                                                      */
/*                                                                      */
/************************************************************************/
Moca * duplicate_moca ( Moca *m);
Moca * free_moca ( Moca *m);
/************************************************************************/
/*                                                                      */
/*            PDB Functions                                             */
/*                                                                      */
/*                                                                      */
/************************************************************************/
Structure       * declare_structure ( int n, char **array);
Structure       * extend_structure ( Structure *S);
/************************************************************************/
/*                                                                      */
/*            Weights Functions                                         */
/*                                                                      */
/*                                                                      */
/************************************************************************/
Weights* declare_weights ( int nseq);
Weights* duplicate_weights (Weights *W);
Weights* free_weights ( Weights* W);
void mem_profile (char *msg);
FILE* print_mem_usage (FILE *fp, char *comment);
void set_max_mem (int m);
int verify_memory (int s);
int my_assert ( void *p, int index);

void * vmalloc ( size_t size);
void * vcalloc ( size_t nobj, size_t size);
void * vcalloc_nomemset ( size_t nobj, size_t size);


void * sub_vcalloc ( size_t nobj, size_t size, int MODE);

void * vrealloc ( void *p, size_t size);
void * vreallocg ( void *p, size_t size, int set, int resize);
void * vrealloc_nomemset ( void *p, size_t size);
void * vrealloc_nomemset_noresize ( void *p, size_t size);
void * vrealloc_noresize ( void *p, size_t size);

void   vfree2 ( void **p);
void   vfree ( void *p);
void * free_arrayN (void *p, int ndim);
void   vfree_all (void *p);
/*********************************************************************/
/*                                                                   */
/*                          SIZES                                    */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
void write_size_short (int x, short  *array, int offset);
void write_size_char  (int x, char   *array, int offset);
void write_size_int   (int x, int    *array, int offset);
void write_size_float (int x, float  *array, int offset);
void write_size_double(int x, double *array, int offset);

int read_size_short ( void  *array, size_t size  );
int read_size_char  ( void  *array, size_t size );
int read_size_int   ( void  *array, size_t size );
int read_size_float ( void  *array, size_t size );
int read_size_double( void  *array, size_t size );
int arrlen ( void  *array);
int read_array_size_new ( void  *array);
int read_array_size ( void  *array, size_t size );
int read_array_new ( void  *array);
int is_dynamic_memory ( void *array);

/*********************************************************************/
/*                                                                   */
/*                          REALLOCATION                             */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
void **realloc_arrayN(int ndim,void **main_array,size_t size, ...);
void **realloc_arrayN2 ( int ndim, void ** p, int *A, size_t size);


void ** realloc_array (void **array,size_t size, int first, int second, int ext1, int ext2);
short     ** realloc_short     ( short     **array, int first, int second, int ext1, int ext2);
char      ** realloc_char      ( char      **array, int first, int second, int ext1, int ext2);
int       ** realloc_int       ( int       **array, int first, int second, int ext1, int ext2);
float     ** realloc_float     ( float     **array, int first, int second, int ext1, int ext2);
double    ** realloc_double    ( double    **array, int first, int second, int ext1, int ext2);
Alignment ** realloc_aln_array ( Alignment **array, int ext1);
/*The new realloc is recommended*/
short     ** new_realloc_short     ( short     **array, int ext1, int ext2);
char      ** new_realloc_char      ( char      **array, int ext1, int ext2);
int       ** new_realloc_int       ( int       **array, int ext1, int ext2);
float     ** new_realloc_float     ( float     **array, int ext1, int ext2);
double    ** new_realloc_double    ( double    **array, int ext1, int ext2);


void * declare_arrayNnomemset (int ndim, size_t size, ...);
void *declare_arrayN2nomemset ( int ndim, int *A, size_t size);

void * declare_arrayN (int ndim, size_t size, ...);
void *declare_arrayN2 ( int ndim, int *A, size_t size);

char       * resize_string (char *buf);
void      ** declare_array     (int first, int second, size_t size);
short     ** declare_short     ( int first, int second);
char      ** declare_char      ( int first, int second);
int       ** declare_int       ( int first, int second);
int       ** declare_int2       ( int first, int *second, int delta);

float     ** declare_float     ( int first, int second);
double    ** declare_double    ( int first, int second);

void      ** declare_array_nomemset     (int first, int second, size_t size);
short     ** declare_short_nomemset      ( int first, int second);
char      ** declare_char_nomemset       ( int first, int second);
int       ** declare_int_nomemset        ( int first, int second);
float     ** declare_float_nomemset      ( int first, int second);
double    ** declare_double_nomemset     ( int first, int second);


Alignment ** declare_aln_array ( int first);

short     **  free_short    ( short     **array, int first);
int       **  free_int      ( int       **array, int first);
char      **  free_char     ( char      **array, int first);
double    ** free_double    ( double    **array, int first);
float     ** free_float     ( float     **array, int first);
Alignment ** free_aln_array ( Alignment **array);
Alignment *free_data_in_aln (Alignment *A);

long aln_stack (Alignment *A, int mode);

Sequence  *free_Alignment     ( Alignment *A);
Sequence  *free_aln     ( Alignment *A);
Alignment *declare_Alignment  ( Sequence  *S);
Alignment *realloc_Alignment4nseq  ( Alignment *A, int nseq);
Alignment *realloc_alignment  ( Alignment *A, int new_len);
Alignment *realloc_alignment2 ( Alignment *A, int new_nseq, int new_len);

Alignment *declare_aln  ( Sequence  *S);
Alignment *declare_aln2 (int nseq, int len);
Alignment *realloc_aln  ( Alignment *A, int new_len);
Alignment *realloc_aln2 ( Alignment *A, int new_nseq, int new_len);
Alignment *update_aln_random_tag (Alignment *A);

Alignment *copy_aln ( Alignment *A, Alignment *B);
Alignment* extract_sub_aln2 ( Alignment *A, int nseq, char **list);
Alignment* extract_sub_aln ( Alignment *A, int nseq, int *list);
Alignment* shrink_aln      ( Alignment *A, int nseq, int *list);

Profile   *copy_profile   (Profile *P1);
Profile   *declare_profile(char *alphabet, int len);
Profile * free_profile ( Profile *P);

Sequence *  reset_sequence_len (Sequence *S);
Sequence  * declare_sequence ( int min, int max, int nseq);
Sequence * realloc_sequence   (Sequence *OUT, int new_nseq, int max_len);
Sequence * duplicate_sequence (Sequence *S );
Sequence * add_sequence ( Sequence *IN, Sequence *OUT, int i);
void free_sequence ( Sequence *LS, int nseq);



Fname *declare_fname ( int );
Fname *free_fname ( Fname *F);

#endif 
