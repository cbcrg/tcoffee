typedef struct
    {
    char *name;
    char *path;
    char *suffix;
    }
Fname;

/*int    memcmp ( const void *a, const void * b, size_t size);
void * memcpy (       void *a,       void * b, size_t size);
*/
/*********************************************************************/
/*                                                                   */
/*                                   HEAPSORT                        */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
FILE *hsort_file     ( FILE *fp ,int n,int len, size_t size,int first_comp_field, int n_comp_fields,int (*compare)(const void *, const void*,int,int,size_t),void * (*copy)(void *,void*,size_t));
void ** hsort_array ( void **ra,int n,int len, size_t size, int first_comp_field, int n_comp_fields,int (*compare)(const void *, const void*,int,int,size_t),void * (*copy)(void *,void*,size_t));
/**********************************************************************/
/*                                                                    */
/*                         HSORT WRAPPERS                             */
/*                                                                    */
/*                                                                    */
/**********************************************************************/
void **hsort_list_array ( void **L, int len, size_t size, int entry_len,int first_comp_field, int n_comp_fields);  
FILE  *hsort_list_file  ( FILE *fp, int len, size_t size, int entry_len,int first_comp_field, int n_comp_fields);  
int hsort_cmp ( const void *a, const void *b, int first, int clen, size_t size);
void *hsort_cpy(void*to, void *from, size_t size);

void test_hsort_list_array();

/*********************************************************************/
/*                                                                   */
/*                         CEDRIC BSEARCH                            */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
void * bsearch_file ( const void *key,int *p,int comp_first,int comp_len, FILE *fp ,int len, int entry_len,size_t el_size, int (*compare)(const void *, const void*,int, int, size_t));
void * bsearch_array( const void *key,int *p,int comp_first,int comp_len,void**list,int len, int entry_len,size_t el_size, int (*compare)(const void *, const void*,int, int, size_t));

/*********************************************************************/
/*                                                                   */
/*                     MY  B_SEARCH_FILE FUNCTIONS                   */
/*                                                                   */
/*                                                                   */
/*********************************************************************/		 
void **search_in_list_file  ( void *key,int *p, int comp_len,FILE *fp, int len, size_t size, int entry_len);
void **search_in_list_array ( void *key,int *p, int comp_len,void **L , int len, size_t size, int entry_len);

/*********************************************************************/
/*                                                                   */
/*                         SORT/COMPARE/SEARCH FUNCTIONS             */
/*                                                                   */
/*                                                                   */
/*********************************************************************/

int **search_in_list_int ( int *key, int k_len, int **list, int ne);
void sort_float ( float **V,int N_F, int F, int left, int right);
int cmp_float ( const float **a, const float **b);
void sort_int_1D ( int *L, int n);
void sort_int ( int **V,int N_F, int F, int left, int right);
void sort_list_int ( int **V,int N_F, int F, int left, int right);
void sort_int_inv ( int **V,int N_F, int F, int left, int right);
void sort_list_int_inv ( int **V,int N_F, int F, int left, int right);
int cmp_int ( const int**a, const int**b);
int cmp_list_int (const int**a, const int**b);
int name_is_in_list ( char *name, char **name_list, int n_name, int len);
char * check_list_for_dup ( char **list, int ne);
FILE *get_number_list_in_file ( FILE *fp, int *list, int *n, int *max_len);

/*********************************************************************/
/*                                                                   */
/*                         DUPLICATION                               */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
short  *  set_short ( short *, int n,...);
char   *  set_char  ( char  *, int n,...);
int    *  set_int   ( int   *, int n,...);
float  *  set_float ( float *, int n,...);
double *  set_double( double*, int n,...);

short  *  ga_memcpy_short ( short  *array1, short  *array2, int n);
int    *  ga_memcpy_int   ( int    *array1, int    *array2, int n);
float  *  ga_memcpy_float ( float  *array1, float  *array2, int n);
double *  ga_memcpy_double( double *array1, double *array2, int n);

short  ** duplicate_short ( short  **array , int len, int field);
int    ** duplicate_int   ( int    **array , int len, int field);
char   ** duplicate_char  ( char   **array , int len, int field);
char    * duplicate_string ( char *string);
float  ** duplicate_float ( float  **array , int len, int field);
double ** duplicate_double( double **array , int len, int field);

short  ** copy_short ( short  **array1, short  **array2, int len, int number_field);
char   ** copy_char  ( char   **array1, char   **array2, int len, int number_field);
int    ** copy_int   ( int    **array1, int    **array2, int len, int number_field);
float  ** copy_float ( float  **array1, float  **array2, int len, int number_field);
double ** copy_double( double **array1, double **array2, int len, int number_field);

/*********************************************************************/
/*                                                                   */
/*                        CONCATENATION                              */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
Alignment ** cat_aln_list ( Alignment **list_to_cat,int first, int end, Alignment **rec_list);

/*********************************************************************/
/*                                                                   */
/*                         NUMBER ARRAY ANALYSE                      */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
FILE * output_array_int (int  **array, int len, int nf ,FILE *fp);
short  return_max_short (short ** array, int len_array, int field);
char  return_max_char( char  ** array, int len_array, int field);
int  return_max_int( int   ** array, int len_array, int field);
float  return_max_float( float ** array, int len_array, int field);
double return_max_double( double** array, int len_array, int field);

short  return_min_short (short ** array, int len_array, int field);
char  return_min_char( char  ** array, int len_array, int field);
int  return_min_int( int   ** array, int len_array, int field);
float  return_min_float( float ** array, int len_array, int field);
double return_min_double( double** array, int len_array, int field);

short  return_max_coor_short (short ** array, int len_array, int field, int *coor);
char  return_max_coor_char( char  ** array, int len_array, int field, int *coor);
int  return_max_coor_int( int   ** array, int len_array, int field, int *coor);
float  return_max_coor_float( float ** array, int len_array, int field, int *coor);
double return_max_coor_double( double** array, int len_array, int field, int *coor);

short  return_min_coor_short (short ** array, int len_array, int field, int *coor);
char  return_min_coor_char( char  ** array, int len_array, int field, int *coor);
int  return_min_coor_int( int   ** array, int len_array, int field, int *coor);
float  return_min_coor_float( float ** array, int len_array, int field, int *coor);
double return_min_coor_double( double** array, int len_array, int field, int *coor);

short  return_2Dmax_short (short ** array, int start, int len_array, int first_field, int number_field);
char  return_2Dmax_char( char  ** array, int start, int len_array, int first_field, int number_field);
int  return_2Dmax_int( int   ** array, int start, int len_array, int first_field, int number_field);
float  return_2Dmax_float( float ** array, int start, int len_array, int first_field, int number_field);
double return_2Dmax_double( double** array, int start, int len_array, int first_field, int number_field);

short  return_2Dmin_short (short ** array, int start, int len_array, int first_field, int number_field);
char  return_2Dmin_char( char  ** array, int start, int len_array, int first_field, int number_field);
int  return_2Dmin_int( int   ** array, int start, int len_array, int first_field, int number_field);
float  return_2Dmin_float( float ** array, int start, int len_array, int first_field, int number_field);
double return_2Dmin_double( double** array, int start, int len_array, int first_field, int number_field);

short  return_2Dmax_coor_short ( short ** array,int start1, int end1, int start2, int end2, int *i, int *j );
char  return_2Dmax_coor_char( char  ** array, int start1, int end1, int start2, int end2, int *i, int *j);
int  return_2Dmax_coor_int( int   ** array, int start1, int end1, int start2, int end2, int *i, int *j);
float  return_2Dmax_coor_float( float ** array, int start1, int end1, int start2, int end2, int *i, int *j);
double return_2Dmax_coor_double( double** array, int start1, int end1, int start2, int end2, int *i, int *j);

short  return_2Dmin_coor_short ( short ** array, int start1, int end1, int start2, int end2, int *i, int *j);
char  return_2Dmin_coor_char( char  ** array, int start1, int end1, int start2, int end2, int *i, int *j);
int  return_2Dmin_coor_int( int   ** array, int start1, int end1, int start2, int end2, int *i, int *j);
float  return_2Dmin_coor_float( float ** array, int start1, int end1, int start2, int end2, int *i, int *j);
double return_2Dmin_coor_double( double** array, int start1, int end1, int start2, int end2, int *i, int *j);

double return_wmean_short ( short ** array, int len, int wfield, int field);
double return_wmean_char  ( char  ** array, int len, int wfield, int field);
double return_wmean_int   ( int   ** array, int len, int wfield, int field);
double return_wmean_float ( float ** array, int len, int wfield, int field);
double return_wmean_double( double** array, int len, int wfield, int field);

double return_mean_short  ( short ** array, int len, int field);
double return_mean_char   ( char  ** array, int len, int field);
double return_mean_int    ( int   ** array, int len, int field);
double return_mean_float  ( float ** array, int len, int field);
double return_mean_double ( double** array, int len, int field);

short  return_sum_short ( short ** array, int len, int field);
char   return_sum_char  ( char  ** array, int len, int field);
int    return_sum_int   ( int   ** array, int len, int field);
float  return_sum_float ( float ** array, int len, int field);
double return_sum_double( double** array, int len, int field);

short  return_sd_short ( short ** array, int len, int field, short mean);
char   return_sd_char  ( char  ** array, int len, int field, char mean);
int    return_sd_int   ( int   ** array, int len, int field, int mean);
float  return_sd_float ( float ** array, int len, int field, float mean);
double return_sd_double( double** array, int len, int field, double mean);

double return_z_score ( double x, double sum, double sum2, double n);

short*  invert_list_short ( short * array, int len );
char*   invert_list_char  ( char  * array, int len );
int*    invert_list_int   ( int   * array, int len );
float*  invert_list_float ( float * array, int len );
double* invert_list_double( double* array, int len );

void   swap_short ( short * array, short * array2,int len );
void   swap_char  ( char  * array, char  * array2,int len );
void   swap_int   ( int   * array, int   * array2,int len );
void   swap_float ( float * array, float * array2,int len );
void   swap_double( double* array, double* array2,int len );

short  return_max_short_hor (short  ** array, int len_array, int field);
char   return_max_char_hor  (char   ** array, int len_array, int field);
int    return_max_int_hor   (int    ** array, int len_array, int field);
float  return_max_float_hor (float  ** array, int len_array, int field);
double return_max_double_hor(double ** array, int len_array, int field);

short  return_min_short_hor ( short ** array, int len_array, int field);
char   return_min_char_hor  ( char  ** array, int len_array, int field);
int    return_min_int_hor   ( int   ** array, int len_array, int field);
float  return_min_float_hor ( float ** array, int len_array, int field);
double return_min_double_hor( double** array, int len_array, int field);

short  best_short (int n, ...);
int    best_int   (int n, ...);
char   best_char  (int n, ...);
float  best_float (int n, ...);
double best_double(int n, ...);

int  is_defined_short (int n, ...);
int  is_defined_int   (int n, ...);
int  is_defined_char  (int n, ...);
int  is_defined_float (int n, ...);
int  is_defined_double(int n, ...);



int return_maxlen ( char ** array, int number);
int return_minlen ( char ** array, int number);

float return_mean_diff_float ( float **array, int len, int field,float mean);


void inverse_int ( int**array, int len, int field, int max, int min);
void inverse_float ( float**array, int len, int field, int max, int min);
void inverse_2D_float ( float **array, int start, int len, int start_field, int number_field, float max,float min);

/*********************************************************************/
/*                                                                   */
/*                          SIZES                                    */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
void write_size_short (int x, short  *array, int ofset);
void write_size_char  (int x, char   *array, int ofset);
void write_size_int   (int x, int    *array, int ofset);
void write_size_float (int x, float  *array, int ofset);
void write_size_double(int x, double *array, int ofset);

int read_size_short ( short  *array, int ofset);
int read_size_char  ( char   *array, int ofset);
int read_size_int   ( int    *array, int ofset);
int read_size_float ( float  *array, int ofset);
int read_size_double( double *array, int ofset);

/*********************************************************************/
/*                                                                   */
/*                         SHELL INTERFACES                          */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
char* get_env_variable ( const char *var, int mode);
void setenv_func ( char *string_name, char *string_value);
void get_pwd ( char *name);
int pg_is_installed ( char *pg);
/*********************************************************************/
/*                                                                   */
/*                           MISC                                    */  
/*                                                                   */
/*********************************************************************/
char *num2plot (int value, int max, int line_len);
float grep_function ( char *pattern, char *file);
void crash_if ( int val, char *s);
void crash ( char *s);
int ** make_recursive_combination_table ( int tot_n_param, int *n_param, int *nc, int**table, int field);
/*********************************************************************/
/*                                                                   */
/*                         STRING PROCESSING                         */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
FILE *print_array_char (FILE *out, char **array, int n, char *sep);
char *extract_suffixe ( char *array);
Fname* parse_fname ( char *array);

void string_array_convert ( char **array, int n_strings, int ns, char **sl);
void string_convert( char *string, int ns, char **sl);
int convert ( char c, int ns, char **sl);
int convert2 ( char c, char *list);

void string_array_upper ( char **string, int n);
void string_array_lower ( char **string, int n);
char *upper_string ( char *string);
char *lower_string ( char *string);

int str_overlap ( char *string1, char *string2, char x);
int get_string_line ( int start, int n_lines, char *in, char *out);
FILE * output_string_wrap ( int wrap,char *string, FILE *fp);
char * extract_char ( char * array, int first, int len);
char** break_list ( char **argv, int *argc, char *separators);
char ** get_list_of_tokens ( char *string, char *separators, int *n_tokens);
char **ungap_array(char ** array, int n);
void ungap ( char *seq);
char *mark_internal_gaps(char *seq, char symbol);

char **char_array2number ( char ** array, int n);
char *char2number ( char * array);

int isblanc ( char *buf);
void splice_out ( char *seq, char x);
int is_number ( char *buf);
int is_alpha_line ( char *buf);
int is_alnum_line ( char *buf);

int get_string_sim ( char *string1, char *string2, char *ignore);
int is_gap ( char x);
int is_aa  ( char x);
int is_dna ( char x);

char * get_alphabet   ( char *seq, char *alphabet);
int is_in_set ( char r, char *list);
int array_is_in_set (char *array, char *set);
char * generate_void ( int x);
char * generate_null ( int x);
char * generate_string ( int x, char y);
void translate_name ( char *name);
int get_longest_string  (char **array,int n, int *len, int *index);
int get_shortest_string (char **array,int n, int *len, int *index);
/*EDIT STRING*/
char **pad_string_array ( char **array, int n, int len, char pad);
char * crop_string (char *string, int start, int end);
int get_distance2char ( char *x, char *list);

/*********************************************************************/
/*                                                                   */
/*                         TIME        FUNCTIONS                     */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
int  get_time  ();
int get_ctime  ();
int  reset_time();
int increase_ref_time(int increase);
/*********************************************************************/
/*                                                                   */
/*                         SYSTEM CALLS                              */
/*                                                                   */
/*                                                                   */
/*********************************************************************/	
int evaluate_sys_call_io ( char *out_file, char *com, char *fonc);

/*********************************************************************/
/*                                                                   */
/*                         IO FUNCTIONS                              */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
char ** standard_initialisation ( char **in_argv, int *in_argc);
int   count_n_res_in_array  (char *array, int len);
int   count_n_gap_in_array  (char *array, int len);
int count_n_symbol_in_array ( char *array, char *array_list, int len);

int   count_n_char_x_in_file(char *name, char x);
int   count_n_char_in_file(char *name);
int   count_n_line_in_file(char *name);
int measure_longest_line_in_file ( char *name );

char ** get_parameter ( char *para_name, int *np, char *fname);
int get_cl_parameter (int argc, char **argv,FILE **fp, char *para_name, int *set_flag, char *type, int optional, int max_n_val,char *usage, ...); 
char *input_name ();
FILE * vtmpfile();
char * vtmpnam ( char *s);
char *  tmpnam_2 (char *s);
void error_function ();
void   clean_function ( );
void   sig_clean_function ( int x);
FILE * vfopen ( char *name, char *mode);
FILE * vfclose (FILE *fp);
int echo ( char *string, char *fname);

int **get_file_block_pattern (char *fname, int *n_blocks, int max_n_line);

FILE * find_token_in_file_nlines ( char *fname, FILE * fp, char *token, int n_line);
FILE * find_token_in_file ( char *fname, FILE * fp, char *token);

FILE * set_fp_after_char ( FILE *fp, char x);
FILE * set_fp_id ( FILE *fp, char *id);

int check_program_is_installed ( char *program_name, char *current_path, char *path_variable, char *where2getit, int fatal);
int check_file_exists ( char *fname);
void create_file ( char *name);
void delete_file ( char *fname);
int  util_rename ( char* from, char *to);
int  util_copy   ( char* from, char *to);
FILE * output_completion ( FILE *fp,int n, int tot, int n_eports);
void * null_function (int a, ...);
int  btoi ( int nc,...);
#include <stdio.h>

unsigned long linrand(unsigned long r);
unsigned long addrand(unsigned long r);
void addrandinit(unsigned long s);

static unsigned long mult(unsigned long p,unsigned long q);
    

#define SEQ1 0
#define SEQ2 1
#define R1 2
#define R2 3
#define WE 4
#define CONS 5
#define MISC 6
#define LIST_N_FIELDS 7
#define CLIST_TYPE int


/*********************************************************************************************/
/*                                                                                           */
/*         STRUCTURES FOR PDB ANALYSIS                                                       */
/*                                                                                           */
/*********************************************************************************************/

typedef struct
    {
	int   n_excluded_nb;
	
	float similarity_threshold;
	float rmsd_threshold;
	int   distance_on_request;
	char  *comparison_io;
	float maximum_distance;/*Diameter of the bubble used to identify the Calpha Neighborhood*/
	int   N_ca;            /*Number of Calpha to be looked at on both side*/
	float max_delta ;      /*Maximum value for delta to be positive*/ 
	char *hasch_mode;
        int   scale;             /*Value substracted to the pdb score in the bubble mode*/
        int   n_extra_param;
        char **extra_param;
    }
Pdb_param;

typedef struct
    {
	int num;
	int res_num;/*Residue number from 1 to N*/
	char type[4];
	float  x;
	float  y;
	float  z;
    }
Atom;

typedef struct
    {
      
      Atom*CA;
      Atom *C;
      Atom *N;
      Atom *CB;
    }
Amino_acid;


typedef struct
    {
    /*Distances used for the Neighbour mode*/
	int    **nb;       /*Neighbors of each Ca ( sorted by distance) given as atoms*/
	                   /*nb[x][0] contains the number of neighbor atoms*/
	float  **d_nb;     /* contains the distances between atom y=nb[x][5] and Ca x*/
	                   /* !!!d_nb[x][0] is empty, the array starts at +1 to folow nb*/
	int max_nb;        /* Largest neigborhood*/ 
}
Struct_nb;

typedef struct 
    {
        
	int   len;         /*Number of Calpha Carbons*/
	int   n_atom;      /*Number of atoms*/
	char  *name;       /*Name of the sequence*/
	char  *seq;        /*Sequence ( Complete)*/
	Atom  **structure; /*Atoms*/
        Atom  **ca;        /*List of pointers to the Calpha Atoms from 0 to N-1*/
        Amino_acid **peptide_chain;/*List of pointers to the Calpha Atoms from 0 to N-1*/
      
        
        Struct_nb *Chain;
        Struct_nb *Bubble;
        Struct_nb *Transversal;
        
        float ** ca_dist;
	Pdb_param *pdb_param;
}

Ca_trace;
/*********************************************************************************************/
/*                                                                                           */
/*         MOCA: Data structure for domains and alignments                                   */
/*                                                                                           */
/*********************************************************************************************/
struct Moca
{
  /*Normalisation factor: value by which each constraint weight is decreased*/
      int moca_scale;
  /*Functions used for domain extraction:*/
      /*Function for evaluating the score of a domain: returns 0 if not acceptable, value if OK*/
      int (*evaluate_domain)(Alignment*,struct Constraint_list *);
      int moca_threshold;

      /*Function for hiding previously used residues*/
      int  ** (*cache_cl_with_domain)(Alignment*, struct Constraint_list *);
      int  **forbiden_residues; /*List of residues already used for domain construction*/
      
      
      /*Function for trunkating the result into a non-overlapping alignment*/
      Alignment* (*make_nol_aln)(Alignment*, struct Constraint_list *);
     
      /*Parameters Coordinates of the first motif to extract*/
      int moca_start;
      int moca_len;
      int moca_interactive;
      
};
typedef struct Moca Moca;
/*********************************************************************************************/
/*                                                                                           */
/*         CONSTRAINT LISTS                                                                  */
/*                                                                                           */
/*********************************************************************************************/

struct Constraint_list
    {
      /*In Case of Modif, synchronize with:
	util_declare/declare_constraint_list
	util_declare/cache_dp_value4constraint_list
	util_declare/duplicate_constraint_list
	util_declare/free_constraint_list
      */

      Sequence *S;         /*Total sequences*/
      Sequence *STRUC_LIST; /*Name of the sequences with a Structure*/
      Sequence *DO_S;      /*Sequences to align*/
      Weights  *W;         /*Sequence Weights*/
      int *translation;   
      char **  out_aln_format;
      int    n_out_aln_format;

      
      /*Packing Sequence: To use with domain analysis*/
      int **packed_seq_lu;
      
      /*DATA*/
      FILE *fp;           /*File used for i/o if disk being used*/
      int **L;            /*Array used for storing Lib if mem being used*/
      int **M;            /*substitution matrix*/
      
      /*List Information*/      
      int ne;             /*Number of elements in the list*/
      char *list_name;    /*Name of the list*/
      int  entry_len;     /*Size of an entry in el_size*/
      size_t el_size;     /*Size of each elements in an entry in bytes*/
      
      /*Normalisation information*/
      int normalise;
      int max_ext_value;
      int max_value;


      /*Pair wise alignment method*/
      int   pw_parameters_set;
      int   gop;
      int   gep;
      int   f_gop;
      int   f_gep;
      int   nm_gop;
      int   nm_gep;
      
      int   nomatch;
      
      int   TG_MODE;
      int   F_TG_MODE;

      char  dp_mode[100];
      int   maximise;
      char  matrix_for_aa_group[100];
      float diagonal_threshold;
      int ktup;
      int use_fragments;
      int fasta_step;
      int lalign_n_top;
      int sw_min_dist;
      char **matrices_list;
      int n_matrices;
      
      /*Functions used for dynamic programming and Evaluation*/
      /*1 Function for evaluating the cost of a column*/
      int (*get_dp_cost)(Alignment*, int **, int, int*, int, int **, int, int*, int, struct Constraint_list *);
      /*2 Function for evaluating the cost of a pair of residues*/
      int (*evaluate_residue_pair)(struct Constraint_list *, int, int, int, int);
      /*3 Function for making dynamic programming*/
      int (*pair_wise)(Alignment *, int*, int **,struct Constraint_list *);
      int weight_field;


      /*Extend a sequence against itself*/

/*Parameters for domain extraction*/      
      Moca *moca;
      
/*Functions for hiding forbiden pairs of residues*/
      int ****forbiden_pair_list;     /* pair_list[S1][S2][L1][L2]=1 ->forbiden*/
				     /* pair_list[S1][S2][L1][L2]=0 ->allowed*/
				     /* pair_list[S1][S2][L1]=NULL  ->all pairs S1L1, S2 allowed */
                                     /* S-> sequences, 0..N   */
				     /* L-> residues , 1..L-1 */

/*extention properties:  copy*/
      int *seq_for_quadruplet;
      int nseq_for_quadruplet;
      
/*extention properties: Do Not copy*/
      int extend_jit;               /*Extend only on request*/
      int extend_threshold;         /*Do not extend pairs below the Theshold*/
      int do_self;                  /*Extend a sequence against itself*/
      char extend_clean_mode[100];  
      char extend_compact_mode[100];
      


/*Lookup table parameteres*/
/*!!!!!do not copy in duplication*/
      /*Residue Index contains residue_index[nseq][seq_len][0]->number of links*/
      /*[seq][res][x  ]->target seq (0->N-1)*/
      /*[seq][res][x+1]->traget res (1->len*/
      /*[seq][res][x+2]->target weight */
      /*It is automatically recomputed when L residue_indexed is set to 0*/
      int residue_indexed;
      int ***residue_index;
      int residue_field;	

      /*Index of the pairs of sequences within L*/
      int seq_indexed;
      int **start_index;
      int **end_index;
      int max_L_len;
      int chunk;
      
      
            
      /*PDB STRUCTURE ALIGNMENTS*/      
      Ca_trace ** T;	/*This structure contains the PDB trace for sequences with a known Struc T[Nseq]*/

       /*MISC*/
      int cpu;
      FILE *local_stderr;
      
};

typedef struct Constraint_list Constraint_list;
/*********************************************************************/
/*                                                                   */
/*                         PRODUCE IN LIST                             */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
Constraint_list *produce_list ( Constraint_list *CL, Sequence *S, char * method,char *weight,char *mem_mode);	
int parse_method ( char *fname, Sequence *S, char ***aln_c, char ***seq_c, char *seq, char *aln, char *aln_mode, char *out_mode);
int is_in_pre_set_method_list (char *fname);
char *make_aln_command(char *command, char *parafile, char *seq, char *aln);

/*********************************************************************/
/*                                                                   */
/*                         WRITE IN LIST                             */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
int vread_clist ( Constraint_list *CL, int a, int b );
int vwrite_clist ( Constraint_list *CL, int a, int b, CLIST_TYPE x);
Constraint_list *index_constraint_list ( Constraint_list *CL);
Constraint_list *index_res_constraint_list ( Constraint_list *CL, int field);
/*********************************************************************/
/*                                                                   */
/*                         ENTRY MANIPULATION                        */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
Constraint_list * add_list_entry2list (Constraint_list *CL, int n_para, ...);
Constraint_list * evaluate_constraint_list_reference ( Constraint_list *CL);
Constraint_list *add_entry2list ( CLIST_TYPE *entry, Constraint_list *CL);
Constraint_list *insert_entry2list ( CLIST_TYPE *entry, int pos,Constraint_list *CL);
CLIST_TYPE* extract_entry(CLIST_TYPE * entry, int pos, Constraint_list *CL);
/*********************************************************************/
/*                                                                   */
/*                         LIST EXTENTION                            */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
Constraint_list *extend_list_pair (Constraint_list *CLin,char *store_mode, int s1, int s2);
Constraint_list *extend_list (Constraint_list *CLin, char *store_mode,char *clean_mode, char *compact_mode,int do_self, Sequence *SUBSET);
void get_bounds (Constraint_list *CL, int s1, int s2, int *start, int *end);
int ** fill_pos_matrix (Constraint_list *CL, int beg, int end, int slen, int **pos, int *len, int mirrored);

/*********************************************************************/
/*                                                                   */
/*                         SEARCH IN LIST (ARRAY AND FILE)           */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
FILE * compare_list (FILE *OUT, Constraint_list *CL1,Constraint_list *CL2);
CLIST_TYPE **search_in_list_constraint(int *key, int k_len, int **L, int ne, int ***start_index, int ***end_index);
CLIST_TYPE **main_search_in_list_constraint ( int *key,int *p,int k_len,Constraint_list *CL);
Constraint_list *sort_constraint_list_inv (Constraint_list *CL, int start, int len);
Constraint_list *invert_constraint_list (Constraint_list *CL, int start,int len);
Constraint_list * sort_constraint_list (Constraint_list *CL, int start, int len);
Constraint_list * sort_constraint_list_on_n_fields (Constraint_list *CL, int start, int len, int first_field, int n_fields);

/*********************************************************************/
/*                                                                   */
/*                         INPUT/OUTPUT                              */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
Constraint_list* read_n_constraint_list(char **fname,int n_list, char *in_mode,char *mem_mode,char *weight_mode,char *type, FILE *local_stderr);
Constraint_list* read_constraint_list(Constraint_list *CL,char *fname,char *in_mode,char *mem_mode,char *weight_mode);


int        read_cpu_in_n_list(char **fname, int n);

Sequence * read_seq_in_n_list(char **fname, int n, char *type);
int        read_cpu_in_list ( char *fname);
int ** read_list ( char *fname, int **list,int *ne, int *nseq, int *cpu, char ***sequences, char ***seq_name);
Constraint_list * read_constraint_list_file(Constraint_list *CL, char *fname);
/*********************************************************************/
/*                                                                   */
/*                         EXTENDED LIST OUTPUT                      */
/*                                                                   */
/*                                                                   */
/*********************************************************************/	
FILE * save_extended_constraint_list      (  Constraint_list *CL, char *mode, FILE *fp) ;
FILE * save_extended_constraint_list_pair (  Constraint_list *CL, char *mode, char* seq1, char * seq2,FILE *fp);

/*********************************************************************/
/*                                                                   */
/*                         LIST OUTPUT                               */
/*                                                                   */
/*                                                                   */
/*********************************************************************/	
FILE * save_constraint_list ( Constraint_list *CL,int start, int len, char *fname, FILE *fp,char *mode,Sequence *S);
FILE * save_sub_list_header ( FILE *OUT, int n, char **name, Constraint_list *CL);
FILE * save_list_header ( FILE *OUT,Constraint_list *CL);
FILE * save_list_footer (FILE *OUT,Constraint_list *CL);
FILE * save_constraint_list_ascii ( FILE *OUT,Constraint_list *CL, int start,int len, int *translation);
FILE * save_constraint_list_bin   ( FILE *OUT,Constraint_list *CL, int start,int len, int *translation);

/*********************************************************************/
/*                                                                   */
/*                         LIST CONVERTION                           */
/*                                                                   */
/*                                                                   */
/*********************************************************************/		
Constraint_list *aln_file2constraint_list (char *alname, Constraint_list *CL,char *weight_mode);
Constraint_list *aln2constraint_list      (Alignment *A, Constraint_list *CL,char *weight_mode);

double **list2mat (Constraint_list *CL,int s1,int s2, double *min, double *max);
Constraint_list * constraint_list2bin_file(Constraint_list *clist);
FILE * bin_file2constraint_list ( Constraint_list *CL, FILE *fp, char *name);

int **list2residue_total_weight ( Constraint_list *CL);
int **list2residue_total_extended_weight ( Constraint_list *CL);
int **list2residue_partial_extended_weight ( Constraint_list *CL);
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                              clean functions                                            */
/*                                                                                         */
/*                                                                                         */
/*                                                                                         */
/*******************************************************************************************/
Constraint_list *clean ( char *clean_mode,Constraint_list *C,int start, int len);
Constraint_list * clean_shadow ( Constraint_list *CL, int start, int len);

/*********************************************************************/
/*                                                                   */
/*                         LIST FUNCTIONS                            */
/*                                                                   */
/*                                                                   */
/*********************************************************************/	
Constraint_list *modify_weight( Constraint_list *CL,int start, int end,  char *modify_mode);
Constraint_list *compact_list (Constraint_list *CL, int start, int len, char *compact_mode);
Constraint_list *rescale_list_simple (Constraint_list *CL,int start, int len,int new_min, int new_max);
Constraint_list *rescale_list (Constraint_list *CL,int start, int len,int max1, int max2);
Constraint_list* filter_list (Constraint_list *CL, int start, int len,int T);
Constraint_list *undefine_list (Constraint_list *CL);
/*********************************************************************/	
/*          DEBUG                                                    */
/*                                                                   */
/*********************************************************************/	
void print_CL_mem(Constraint_list *CL, char *function);
int constraint_list_is_sorted ( Constraint_list *CL);
void check_seq_pair_in_list(Constraint_list *CL,int seq1, int seq2);
/******************************************************************/
/*                    NEW METHODS                                 */
/*                                                                */
/*                                                                */
/******************************************************************/

Constraint_list * align_coding_nucleotides (char *seq, char *method, char *weight, char *mem_mode, Constraint_list *CL);
/*********************************************************************************************/
/*                                                                                           */
/*         FUNCTION FOR WEIGHTING THE LIST                                                   */
/*                                                                                           */
/*********************************************************************************************/
Constraint_list *weight_constraint_list(Constraint_list * CL, char *seq_weight);
Weights* compute_t_coffee_weight(Constraint_list * CL);
Constraint_list *re_weight_constraint_list(Constraint_list * CL,Weights *W);
int ** minimise_repeat_coor (int **coor, int nseq, Sequence *S);
int ** get_nol_seq( Constraint_list *CL,int **coor, int nseq, Sequence *S);

int compare_pos_column( int **pos1,int p1, int **pos2,int p2, int nseq);
Sequence  *  seq2blast_profile ( Sequence *S, int n); 
int       ** seq2aln_pos      (Alignment *A, int *n, int **ls);
Alignment *local_maln2global_maln (char *seq, Alignment *A);


Alignment *  aln2profile (Alignment * A);
int       ** aln2pos_simple   (Alignment *A, int n_nseq, ...);
int       ** aln2pos_simple_2 (Alignment *A);
Alignment ** split_seq_in_aln_list ( Alignment **aln, Sequence *S, int l_seq, char **seq_list);

Sequence * fill_sequence_struc ( int nseq, char **sequences, char **seq_name);
int  seq_list2fasta_file( Sequence *S,  char *list, char *file);
Structure * seq2struc    ( Sequence *S, Structure *ST);
Alignment *strings2aln (int nseq,...);
Alignment * seq2aln      ( Sequence *S, Alignment *A,int rm_gap);
Alignment  *seq_coor2aln ( Sequence *S, Alignment *A, int **coor, int nseq);

Alignment *stack_aln (Alignment *A, Alignment *B);
Alignment *chseqIaln(char *name, int seq_n, int start,int len,Sequence *S, int seqIaln, Alignment *A);

char *dna_aln2cons_seq ( Alignment *A);
char *aln2cons_seq ( Alignment *A, int ns, int *ls, int n_groups, char **group_list);
char *aln2cons_seq_mat ( Alignment*A, char *mat_name);

Alignment *aln2number (Alignment *A);
Alignment * filter_aln ( Alignment *A, Alignment *ST, int value);
Alignment * filter_aln_lower_upper ( Alignment *A, Alignment *ST, int value);
Alignment * filter_aln_upper_lower ( Alignment *A, Alignment *ST, int value);
Alignment * filter_aln_convert ( Alignment *A, Alignment *ST, int value, int n_symbol,char** symbol_list);
int  * count_in_aln ( Alignment *A, Alignment *ST, int value, int n_symbol,char **symbol_list, int *table);
 
Alignment * add_align_seq2aln ( Alignment *A, char *seq, char *seq_name);
Sequence  * aln2seq    ( Alignment *A);

Alignment * expand_aln (Alignment *A);
Alignment * expand_number_aln (Alignment *A);
Alignment * remove_gap_column ( Alignment *A, char *mode);
Alignment*  ungap_sub_aln        ( Alignment *A, int nseq, int *ls);
Sequence *  ungap_seq       ( Sequence *A);
Alignment * ungap_aln_n        ( Alignment *A, int n);
Alignment * ungap_aln        ( Alignment *A);
void compress_aln     ( Alignment *A);

Alignment * aln_gap2random_aa(Alignment *A);
Alignment * make_random_aln(Alignment *A,int nseq, int len, char *alphabet);
Alignment * add_random_sequence2aln( Alignment *A, char *alphabet);

int ** trim_aln_borders            ( char **seq1, char **seq2, int nseq);
Sequence * trim_aln_seq      ( Alignment  *A, Alignment *B);
Sequence * trim_aln_seq_name ( Alignment  *A, Alignment *B);
Sequence *get_defined_residues( Alignment *A);

Alignment *thread_defined_residues_on_aln ( Alignment *A, Sequence *S1);
Sequence *seq2number (Sequence *S);
Sequence * merge_seq    ( Sequence *IN, Sequence *OUT);
char * seq_name2coor ( char *s, int *start, int *end, char sep);
Sequence *extract_one_seq(char *n,int start, int end, Alignment *S,int keep_name);
Sequence  * extract_sub_seq( Sequence  *COOR, Sequence *S);

Sequence * add_prf2seq  ( Alignment *A, Sequence *S);
Sequence * add_sequence ( Sequence *IN, Sequence *OUT, int i);
Sequence  * trim_seq     ( Sequence   *A, Sequence  *B);
Sequence  * reorder_seq  ( Sequence   *A, char **name, int nseq);
char * concatenate_seq ( Sequence *S, char *conc, int *order);

Alignment * reorder_aln        ( Alignment  *A, char **name, int nseq);

Alignment * fix_aln_seq  ( Alignment *A, Sequence *S);
Alignment * extract_nol_local_aln( Alignment *A, int start, int max_end);
Alignment * extract_aln          ( Alignment *A, int start, int end);
Alignment * trunkate_local_aln   ( Alignment *A);
int get_nol_aln_border ( Alignment *A, int start, int direction);
Alignment ** trim_local_aln ( Alignment *A, int **List, int ne, int **residue_list, Sequence *S);

Alignment * aln_cat ( Alignment *A, Alignment *B);
char * extract_defined_seq ( char *in, int in_of, int in_start, int *aa_def, int dir, int *out_start, char *out_seq);
int verify_aln ( Alignment *A, Sequence *S, char * error);
Alignment * remove_end (Alignment *A);

Alignment *adjust_est_aln ( Alignment *PW, Alignment *M, int s);


/********************************************************************/
/*                                                                  */
/*                   ALIGNMENT ANALYSES                             */
/*                                                                  */
/*                                                                  */
/*                                                                  */
/********************************************************************/
int aln_is_aligned ( Alignment *A);
int* get_cdna_seq_winsim ( int *cache, char *string1, char *string2, char *ignore, char *mode, int *w);
int  get_cdna_seq_sim    ( int *cache, char *string1, char *string2, char *ignore, char *mode);

int* get_seq_winsim ( char *string1, char *string2, char *ignore, char *mode, int *w);
int  get_seq_sim ( char *string1, char *string2, char *ignore, char *mode);
int  get_seq_sim_2 ( char *string1, char *string2, char *ignore, char **gr, int ng);
int  get_seq_sim_3 ( char *string1, char *string2, char *ignore, int **mat);


int *** get_winsim_aln_array ( Alignment *A, char *mode, int ***w);
int ** get_sim_aln_array ( Alignment *A, char *mode);
int * get_aln_col_weight ( Alignment *A, char *mode);
int analyse_aln_column   ( Alignment *B, int col);

int **aln2count_mat   (Alignment *A);
char *aln2random_seq (Alignment *A, int noise1, int noise2, int noise3, int gap_noise);

Alignment * master_trimseq( Alignment *A, Sequence *S,char *mode);
Alignment * trimseq( Alignment *A, Sequence *S, char *mode);

float ** get_weight ( Alignment *A, Sequence *S, char *mode);
float **seq2pwsim (	   Alignment *A, Sequence *S, char *mode);
Alignment * trimseq( Alignment *A, Sequence *S,char *mode);

Alignment* seq2subseq3( Alignment *A, Sequence *S,int use_aln, int lower_sim,int upper_sim, int min_nseq, int trim_direction, char *weight_mode, float ***sim_weight, int *seq_list);
Alignment* seq2subseq2( Alignment *A, Sequence *S,int use_aln, int lower_sim,int upper_sim, int max_nseq, int trim_direction, char *weight_mode, float ***weight_table, int *seq_list);
float extreme_seq (int direction, Alignment *A,float **sim_weight,int *seq_list, int *seq_index);

Alignment* seq2subseq1( Alignment *A, Sequence *S,int use_aln, int percent,int max_nseq,int max_diff, char *weight_mode);
/********************************************************************/
/*                                                                  */
/*			AMINO ACID FUNCTIONS                        */
/*                                                                  */
/*                                                                  */
/*                                                                  */
/********************************************************************/
int is_in_same_group_aa ( char r1, char r2, int n_group, char **gl, char *mode);
int find_group_aa_distribution (char *col, int nseq,int n_group, char **gl,  int *distrib, char *mode );
char** make_group_aa (int *ngroup, char *mode);


char * test_gene2prot (Constraint_list *CL, int s1);
Alignment* gene2prot (Alignment *A);
Alignment * dna_aln2_3frame_cdna_aln(Alignment *A,int *ns,int **l_s);
Constraint_list * mask_list_with_aln (Alignment *A,int start, int len,Constraint_list *CL, int new_value);
Constraint_list* mask_list_with_aln_pair (Alignment *A,int start, int end,Constraint_list *CL,int new_value);
Constraint_list *mask_entry( Constraint_list *CL, int p, int new_value);
Constraint_list *prepare_list_and_seq4sw(Constraint_list *I, int n_seq, char **seq_name);
int ** get_undefined_list (Constraint_list *CL);
int      is_never_undefined (Constraint_list *CL,int r);
int* do_analyse_list ( Constraint_list *CL);
void print_list(Constraint_list *CL);
void print_pair (Constraint_list *CL,int p);
int** bin_list (Constraint_list *CL,int field, int Threshold);
void   save_full_list (Constraint_list *CL, char*fname);
FILE * output_list ( Constraint_list *CL, FILE *fp);
FILE * output_pair (Constraint_list *CL,int p, FILE *fp);
NT_node ** make_nj_tree (  Alignment *A,int **distances,int gop, int gep, char **out_seq, char **out_seq_name, int out_nseq, char *tree_file, char *tree_mode);
NT_node ** int_dist2nj_tree (int **distances, char **out_seq_name, int out_nseq,  char *tree_file);
NT_node ** float_dist2nj_tree (float **distances, char **out_seq_name, int out_nseq,  char *tree_file);
NT_node ** dist2nj_tree (double **distances, char **out_seq_name, int out_nseq,  char *tree_file);

void nj_tree(char **tree_description);
void print_phylip_tree(char **tree_description, FILE *tree, int bootstrap);
void two_way_split(char **tree_description, FILE *tree, int start_row, int flag, int bootstrap);
void guide_tree(char *fname, double **saga_tmat, char **sag_seq_name, int saga_nseq);

/************************************************************************/
/*                                                                      */
/*            CONSTRAINT_LIST                                           */
/*                                                                      */
/*                                                                      */
/************************************************************************/ 
Constraint_list * declare_constraint_list ( Sequence *S, char *name, int **L, int ne,FILE *fp, int **M);
Constraint_list *cache_dp_value4constraint_list ( char mode[],Constraint_list *CL);
Constraint_list *duplicate_constraint_list (Constraint_list *CL);
Sequence        * free_constraint_list (Constraint_list *CL);
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


void * vmalloc ( size_t size);
void * vcalloc ( size_t nobj, size_t size);
void * vrealloc ( void *p, size_t size);
void   vfree ( void *p);

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


short     ** declare_short     ( int first, int second);
char      ** declare_char      ( int first, int second);
int       ** declare_int       ( int first, int second);
float     ** declare_float     ( int first, int second);
double    ** declare_double    ( int first, int second);
Alignment ** declare_aln_array ( int first);

short     **  free_short    ( short     **array, int first);
int       **  free_int      ( int       **array, int first);
char      **  free_char     ( char      **array, int first);
double    ** free_double    ( double    **array, int first);
float     ** free_float     ( float     **array, int first);
Alignment ** free_aln_array ( Alignment **array);

Alignment *free_Alignment     ( Alignment *A);
Alignment *declare_Alignment  ( Sequence  *S);
Alignment *realloc_alignment  ( Alignment *A, int new_len);
Alignment *realloc_alignment2 ( Alignment *A, int new_nseq, int new_len);

Alignment *free_aln     ( Alignment *A);
Alignment *declare_aln  ( Sequence  *S);
Alignment *declare_aln2 (int nseq, int len);
Alignment *realloc_aln  ( Alignment *A, int new_len);
Alignment *realloc_aln2 ( Alignment *A, int new_nseq, int new_len);


Alignment *copy_aln ( Alignment *A, Alignment *B);
Alignment* extract_sub_aln ( Alignment *A, int nseq, int *list);
Alignment* shrink_aln      ( Alignment *A, int nseq, int *list);
Profile   *declare_profile(char *alphabet, int len);
Profile * free_profile ( Profile *P);

Sequence  * declare_sequence ( int min, int max, int nseq);
Sequence * realloc_sequence   (Sequence *OUT, int new_nseq, int max_len);
Sequence * duplicate_sequence (Sequence *S );
Sequence * add_sequence ( Sequence *IN, Sequence *OUT, int i);
void free_sequence ( Sequence *LS, int nseq);



Fname *declare_fname ();
Fname *free_fname ( Fname *F);
/*********************************************************************************************/
/*                                                                                           */
/*         STRUCTURES FOR HSEARCH                                                            */
/*                                                                                           */
/*********************************************************************************************/
#define FIND     0
#define ADD      1
#define DECLARE  1
#define MULTIPLY 2
#define ENTER    3
#define REMOVE   4
#define DELETE   4
#define DESTROY  5
#define DESTROY_STACK 6
struct HaschT
{
  int ne;
  struct Hasch_entry **p;
  int allocated;
  int total;
  int freed;

};
typedef struct HaschT HaschT;

struct Hasch_entry
{
  struct Hasch_entry *next;
  struct Hasch_entry *previous;
  int k;
  double data;
  int tag;
};
typedef struct Hasch_entry Hasch_entry;


HaschT * hcreate ( int n_elements);
HaschT * hdestroy (HaschT *T);
double   hsearch  (HaschT *T,double data, int k, int action);

Hasch_entry * allocate_hasch_entry ( Hasch_entry * p, int action);
char     * normalize_pdb_file  (char *name, char *out_file);
Ca_trace * trim_ca_trace (Ca_trace *st, char *seq );
Ca_trace * read_ca_trace (char *file );
Ca_trace * hasch_ca_trace             ( Ca_trace *T);
Ca_trace * hasch_ca_trace_nb          ( Ca_trace *T);
Ca_trace * hasch_ca_trace_bubble      ( Ca_trace *T);
Ca_trace * hasch_ca_trace_transversal ( Ca_trace *TRACE);

float get_atomic_distance ( Atom *A, Atom*B);
float ** measure_ca_distances(Ca_trace *T);

typedef struct
    {
	Alignment *A;
	Weights *W;
	Sequence *S;
	char format[100];
	char file[100];
	int rm_gap;
}Sequence_data_struc;

typedef struct
    {
	char **symbol_list;
        int n_symbol;
        char *coor_file;
        int rm_gap;
        int keep_case;

}Action_data_struc;

Sequence_data_struc *read_data_structure ( char *in_format, char *in_file,Action_data_struc *RAD); 
Alignment * main_read_aln ( char *name, Alignment *A);
Sequence  * main_read_seq ( char *mname);
int output_format_aln ( char *format, Alignment *A, Alignment *EA,char *name);
void main_output   ( Sequence_data_struc *D1, Sequence_data_struc *D2, Sequence_data_struc *DST, char *out_format, char *out_file);

char * identify_seq_format ( char *file);
char identify_format (char **fname);
char **identify_list_format ( char **list, int n);

int format_is_oligo  ( char *file);
int format_is_msf  ( char *file);
int format_is_fasta( char *file);
int format_is_fasta_aln( char *file);
int format_is_fasta_seq( char *file);
int format_is_pir  ( char *file);
int format_is_pir_aln( char *file);
int format_is_pir_seq( char *file);
int format_is_saga  ( char *file);
int format_is_swissprot (char *name);

int is_seq ( char *name);
int is_aln ( char *name);
int has_pdb (char *name);
int is_blast_file (char *name);
int is_sap_file (char *name);
int is_pdb_file ( char *name);
int is_pdb_name ( char *name, Constraint_list *CL);
char*is_pdb_struc ( char *name, Constraint_list *CL);
int is_matrix (char *name);
int is_lib (char *name);
int is_method ( char *file);

int is_in_format_list ( char *name);
int is_out_format_list ( char *name);
int is_struc_in_format_list ( char *name);
int is_struc_out_format_list ( char *name);


void get_barton_list_tc_seq ( char *in_file);
int process_barton_entry (char *buf, char *name);  

Sequence *read_rna_struc_number ( Alignment *A, char *fname);
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               INPUT WEIGHTS                                             */
/*                                                                                         */
/***************************************************************************************** */
Weights* get_amps_sd_scores ( char *fname);
Weights *read_seq_weight (char **name, int nseq, char* seq_weight);
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               INPUT SEQUENCES                                            */
/*                                                                                         */
/***************************************************************************************** */
Sequence* get_pdb_sequence           ( char *fname);
Sequence* get_struc_gor              ( char *fname);
Sequence* get_dialign_sequence       ( char *fname);
Sequence* get_pima_sequence          ( char *fname);
Sequence* get_sequence_dali          ( char *fname);
Sequence* get_pir_sequence           ( char *fname, char *comment_name);
Sequence* get_fasta_sequence         ( char *fname, char *comment_name);
Sequence* get_fasta_sequence_num     ( char *fname, char *comment_name);
Sequence* get_gor_sequence           ( char *fname, char *comment_name);
Sequence* get_swissprot_sequence     ( char *fname, char *comment_name);

void read_check ( Alignment *A, char *check_file);

void read_aln ( char *fname, Alignment *A);
void read_number_aln ( char *fname, Alignment *A);
void read_blast_aln  ( char *fname, Alignment *A);
void read_msf_aln ( char *fname, Alignment *A);
void read_amps_aln ( char *in_file, Alignment *A);
int get_amps_seq_name ( char **name, char* fname);
Alignment *read_gotoh_aln ( char *fname, Alignment *A);


/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               OUTPUT MATRICES                                           */
/*                                                                                         */
/***************************************************************************************** */
int output_freq_mat ( char *outfile, Alignment *A);
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               OUTPUT P-Values                                           */
/*                                                                                         */
/***************************************************************************************** */	
float output_maln_pval ( char *outfile, Alignment *A);
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               OUTPUT WEIGHTS                                            */
/*                                                                                         */
/***************************************************************************************** */
void output_statistics (char *file, Alignment *A, char *mode);
void output_pw_weights4saga ( Weights *W, float **w_list, char *wfile);
int  output_seq_weights ( Weights *W, char *wfile);
FILE * display_weights (Weights *W, FILE *fp);
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               OUTPUT SEQ                                                */
/*                                                                                         */
/***************************************************************************************** */
FILE * display_sequences_names (Sequence *S, FILE *fp);
void output_pir_seq1 (char *fname, Alignment*A );
void output_pir_seq (char *fname, Alignment*A );
void output_gor_seq (char *fname, Alignment*A );
void output_mult_fasta_seq (char *fname, Alignment*A, int n );
void output_fasta_seq1 (char *fname, Alignment*A );
void output_pir_check (char *fname,int nseq, char **A );
void output_fasta_seq (char *fname, Alignment*A );
void output_gotoh_seq (char *fname, Alignment*A );
void output_est_prf   (char *fname, Alignment *A);
void output_gor_seq (char *fname, Alignment*A );
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               OUTPUT ALN                                                */
/*                                                                                         */
/***************************************************************************************** */
void output_pir_aln    ( char *fname,Alignment*A);
void output_model_aln  ( char *fname,Alignment*A );
void output_fasta_aln  ( char *fname,Alignment*A);
void output_msf_aln    ( char *fname,Alignment*B);
void output_clustal_aln( char *name, Alignment*B);
void output_saga_aln   ( char *name, Alignment*B);
void output_phylip_aln ( char *name, Alignment*B);
void output_mocca_aln  ( char *name, Alignment*B,Alignment*S);
void output_rnalign    (char *out_file, Alignment*A,Sequence *STRUC);
void output_pw_lib_saga_aln (char *lib_name, Alignment *A );
void output_lib        (char *lib_name, Alignment *A );
void output_compact_aln( char *name, Alignment *B);

void print_aln ( Alignment *B);
FILE * output_aln( Alignment *B, FILE *fp);


FILE * output_aln_score ( Alignment *B, FILE *fp);
FILE * output_aln_with_res_number ( Alignment *B, FILE *fp);


FILE* output_Alignment ( Alignment *B, FILE *fp);
FILE* output_Alignment_without_header ( Alignment *B, FILE *fp);
FILE * output_Alignment_score ( Alignment *B, FILE *fp);
FILE * output_Alignment_with_res_number ( Alignment *B, FILE *fp);
void output_constraints ( char *fname, char *mode, Alignment *A);

void output_lalign_header( char *name, Alignment *B);
void output_lalign_aln   ( char *name, Alignment *B);

/****************************************************************************************************/
/***************************                                    *************************************/
/***************************             PROCESSING 		*************************************/
/***************************                                    *************************************/
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                              THREADING                                                  */
/***************************************************************************************** */




Structure * declare_rna_structure_num (Sequence *SA);

char *thread_aa_seq_on_dna_seq( char *s);
void thread_seq_struc2aln ( Alignment *A, Sequence *ST);
Alignment *thread_dnaseq_on_prot_aln (Sequence *S, Alignment *A);
void cache_id ( Alignment *A);



int process_est_sequence ( Sequence *S, int *cluster_list);
char * invert_seq ( char *seq);
int get_best_match ( char *seq1, char *seq2);
int** extract_m_diag_streches ( int ** m, int l1, int l2,char *seq1, char *seq2, int *n_mdiag);
int is_strech ( char *AA, char *seq1, char *seq2, int len, int x, int y);

int search_for_cluster ( int seq, int cluster_number, int *cluster_list, int T, int nseq, int **S);	
int * SHC ( int nseq, int **NST, int **ST);
int mutate_sol (int *sol, int nseq);
int evaluate_sol ( int*sol, int nseq, int **ST, int **NST);	

char **make_symbols ( char *name, int *n);
Alignment *code_dna_aln (Alignment *A);
char* back_translate_dna_codon ( char aa, int deterministic);
int translate_dna_codon ( char *seq, char stop);
char* mutate_amino_acid ( char aa, char *mode);
Alignment * mutate_aln ( Alignment *A, char *r);

Alignment *translate_splice_dna_aln (Alignment *A,Alignment *ST );
Alignment * mutate_cdna_aln ( Alignment *A);

char * translate_dna_seq_on3frame (  char *dna_seq, char stop, char *prot);
char * translate_dna_seq ( char *dna_seq, int frame, char stop, char *prot);

char * back_translate_dna_seq ( char *in_seq,char *out_seq, int mode);     
Alignment *back_translate_dna_aln (Alignment *A);
Alignment *translate_dna_aln (Alignment *A, int frame);
Alignment *clean_gdna_aln (Alignment *A);
Alignment *clean_cdna_aln (Alignment *A);
Alignment *clean_est      (Alignment *A);
/**************************************************************************************************/
/********************************                      ********************************************/
/********************************    PROCESSING        ********************************************/
/*************** ****************                      ********************************************/
/**************************************************************************************************/
/*                                                                                                */
/*                                                                                                */
/*                               OUTPUT MATRICES                                                  */
/*                                                                                                */
/**************************************************************************************************/
void modify_data  (Sequence_data_struc *D1, Sequence_data_struc *D2, Sequence_data_struc *DST, char **action_list,int n_actions, Action_data_struc *RAD);
   
