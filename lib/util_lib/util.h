#include <time.h>
#include <signal.h>
#include <errno.h>
#include <sys/times.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <unistd.h>


typedef struct
    {
    char *name;
    char *path;
    char *suffix;
    char *full;
    }
Fname;

struct Tmpname
    {
    char *name;
    struct Tmpname *next;
    };
/*********************************************************************/
/*                                                                   */
/*                                  SANDBOX
/*                                                                   */
/*                                                                   */
/*********************************************************************/
void sb();
/*********************************************************************/
/*                                                                   */
/*                                  DICHOTOMY                        */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
double dichotomy (double value, double target_value, double middle, double *bottom,double *top);
/*********************************************************************/
/*                                                                   */
/*                                   QSORT                        */
/*                                                                   */
/*                                                                   */
/*********************************************************************/

void qsort(void *base, size_t nmemb, size_t size, int (*compar)(const void *, const void *));
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

void sort_int_1D ( int *L, int n);
char** sort_string_array (char **V, int n);
     
void sort_int            ( int **V,int N_F, int F, int left, int right);
int * flash_sort_int_inv ( int **V,int N_F, int F, int left, int right);
int * flash_sort_int     ( int **V,int N_F, int F, int left, int right);
void sort_list_int ( int **V,int N_F, int F, int left, int right);
void sort_list_int2 ( int **V,int *list,int N_F, int left, int right);

void sort_int_inv    ( int    **V,int N_F, int F, int left, int right);
void sort_double_inv ( double **V,int N_F, int F, int left, int right);
void sort_float_inv  ( float  **V,int N_F, int F, int left, int right);

void sort_list_int_inv ( int **V,int N_F, int F, int left, int right);

int cmp_int    ( const int**a,    const int**b);
int cmp_double ( const double**a, const double**b);
int cmp_float  ( const float **a, const float**b);

int cmp_list_int (const int**a, const int**b);
int cmp_list_int2 (const int**a, const int**b);

int name_is_in_list (const char name[], char **name_list, int n_name, int len);
char * check_list_for_dup ( char **list, int ne);
FILE *get_number_list_in_file ( FILE *fp, int *list, int *n, int *max_len);
/*********************************************************************/
/*                                                                   */
/*                         Non Cherry Pick                           */
/*                                                                   */
/*                                                                   */
/*********************************************************************/

int cherry (int argc, char *argv[]);

/*********************************************************************/
/*                                                                   */
/*                         QUANTILE                                  */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
int quantile ( int argc, char *argv[]);

int quantile_rank (int **list,int field, int n, float p);
/*********************************************************************/
/*                                                                   */
/*                         DUPLICATION                               */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
void *  cmemcpy(void*in, void*out,size_t size );
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

short  ** copy_short ( short  **array1, short  **array2);
char   ** copy_char  ( char   **array1, char   **array2);
int    ** copy_int   ( int    **array1, int    **array2);
float  ** copy_float ( float  **array1, float  **array2);
double ** copy_double( double **array1, double **array2);


short  ** copy_short_new ( short  **array1);
char   ** copy_char_new  ( char   **array1);
int    ** copy_int_new   ( int    **array1);
float  ** copy_float_new ( float  **array1);
double ** copy_double_new( double **array1);
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
int output_array_int (char *name, int **array);
int **input_array_int (char *fname);
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
double* return_r (double **list, int n);
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


int max_int (int*i, ...);

int return_maxlen ( char ** array, int number);
int return_minlen ( char ** array, int number);

float return_mean_diff_float ( float **array, int len, int field,float mean);


void inverse_int ( int**array, int len, int field, int max, int min);
void inverse_float ( float**array, int len, int field, int max, int min);
void inverse_2D_float ( float **array, int start, int len, int start_field, int number_field, float max,float min);




void   **recycle   (void   **A, int l, int cycle);

/*********************************************************************/
/*                                                                   */
/*                         SHELL INTERFACES                          */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
char* getenv4debug ( const char *var);
char* get_env_variable ( const char *var, int mode);
float atofgetenv(const char *x);
int   atoigetenv(const char *x);
void setenv_func ( char *string_name, char *string_value);
char *get_pwd ( char *name);
char *pg2path (char *pg);
int pg_is_installed ( char *pg);
/*********************************************************************/
/*                                                                   */
/*                           MISC                                    */  
/*                                                                   */
/*********************************************************************/
char *num2plot (int value, int max, int line_len);
int   perl_strstr ( char *string, char *pattern);
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
char *strnrchr ( char *s,char x, int n);
int intlen (int n);
char * update_string (char *string1, char *string2);

char* csprintf (char *string1,char *string2, ...);
char* strcatf  (char *string1,char *string2, ...);
char *vcat (char *v1, char *v2);

int strget_param ( char *string, char *param_name, char *param_value, char *format, ...);
int   strim (const char *s1, const char *s2);
char * lstrstr ( char *in, char *token);
char * vstrstr ( char *in, char *token);
char * estrstr ( char *in, char *token,...);
char * festrstr ( char *in, char *token,...);

int strscanf (char *in, char *token, char *format, ...);
int strfcanf (char *file, char *token, char *format, ...);

int match_motif ( char *string, char **motif);

char *after_strstr (char *string, char *token);

char ** push_string (char *val, char **stack, int *nval, int mode);
int vsrand (int val);
int  *randomize_list (int *list, int len, int ncycle);
char **list2random_list (char **name, int n);

int vstrcmp (const char *s1, const char *s2);
int vstrncmp (const char *s1, const char *s2, int n);
FILE *print_array_char (FILE *out, char **array, int n, char *sep);

char *extract_suffixe ( char *array);
char * path2filename ( char *array);
char *filename2path (char *name);
Fname* parse_fname ( char *array);

void string_array_convert ( char **array, int n_strings, int ns, char **sl);
void string_convert( char *string, int ns, char **sl);
int convert ( char c, int ns, char **sl);
int convert2 ( char c, char *list);
char *fname2abs(char*name);
void string_array_upper ( char **string, int n);
void string_array_lower ( char **string, int n);
char *upper_string ( char *string);
char *lower_string ( char *string);
char * substitute_double ( char *string, char *token);
char * substitute ( char *string, char *token, char *replacement);
char * substitute_char ( char *string, char token, char replacement);
char * substituteN ( char *string, char *token, char *replacement, int N);
char * tild_substitute ( char *string, char *token, char *replacement);


char ** clean_string ( int n, char **string);

int str_overlap ( char *string1, char *string2, char x);
int get_string_line ( int start, int n_lines, char *in, char *out);
FILE * output_string_wrap ( int wrap,char *string, FILE *fp);
char * extract_char ( char * array, int first, int len);
int check_cl4t_coffee (int argv, char **argc);

char** break_list ( char **argv, int *argc, char *separators);
char** merge_list ( char **argv, int *argc);
int *name_array2index_array ( char **list1, int n1, char **list2, int n2);
char ** get_list_of_tokens ( char *string, char *separators, int *n_tokens);
char **ungap_array(char ** array, int n);
void ungap ( char *seq);
int seq2len (char *seq, char *pset, char *nset);
int seq2res_len (char *seq);
void remove_charset ( char *seq, char *set);
char *remove_charset_from_file (char *fname, char *set);
char *mark_internal_gaps(char *seq, char symbol);

char *list2string  (char **list, int n);
char *list2string2 (char **list, int n, char* sep);

char ** string2list (char *string);
char ** string2list2(char *string, char *separators);
int  *  string2num_list( char *string);
int  *  string2num_list2( char *string, char *separators);
float*  string2float_list2( char *string, char *separators);
int  *  string2int_list2( char *string, char *separators);
char **char_array2number ( char ** array, int n);
char *char2number ( char * array);
long atop(char *);
char *invert_string (char *string);
char *invert_string2 (char *string);
char *string2inverted_string (char *string);
/* Analyse and Compare Strings*/
int isblanc ( char *buf);
/*int islower (char c);
int isupper (char c);
*/
void splice_out ( char *seq, char x);
char* splice_out_seg ( char *seq,int pos, int len);

int is_number ( char *buf);
int is_alpha_line ( char *buf);
int is_alnum_line ( char *buf);
int case_insensitive_strcmp ( char *string1, char *string2);
int get_string_sim ( char *string1, char *string2, char *ignore);

int is_gap ( char x);
int is_gop (int p, char *s);

int is_aa  ( char x);
int is_dna ( char x);
int is_rna ( char x);
int haslower (char *s);
int hasupper (char *s);


char * get_alphabet   ( char *seq, char *alphabet);
int is_in_set ( char r, char *list);
int array_is_in_set (char *array, char *set);
char * generate_void ( int x);
char * generate_null ( int x);
char * generate_string ( int x, char y);
char *string2null (char *string, int len);

char * translate_string (char *string, char *in, char*out);
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
FILE *print_program_information (FILE *fp, char *comment);
FILE* print_cpu_usage (FILE *fp, char *comment);
void print_exit_success_message ();
void print_exit_failure_message ();

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
int   has_error_lock();
int   is_rootpid();
int is_shellpid(int pid);
int   shift_lock (int from, int to, int from_type,int to_type, int action);
char *lock2name (int pid, int type);
int release_all_locks (int pid);
char *lock (int pid, int type, int action, char *value, ...);
int check_process (const char *com,int pid,int r, int failure_handling);
int max_n_pid();
int assert_pid (pid_t p);
pid_t **declare_pidtable ();
pid_t set_pid (pid_t p);
pid_t vvfork(char *mode);
pid_t vwait (pid_t *p);
int   vwait_npid (int submited, int max, int min);
int kill_child_pid(int pid);

int safe_system (const char * commande);
pid_t  vwaitpid (pid_t p, int *status, int options);

int evaluate_sys_call_io ( char *out_file, char *com, char *fonc);
//char *cvsprintf (char *r,char *format, va_list arg_ptr,... );
void HEREf  (char *file,char *string, ...);
void HERE  (char *string, ...);
void HERE2 (char *string, ...);
void printf_exit  (int exit_code, FILE *fp, char *string, ...);
int printf_file ( char *file, char *mode, char *string, ...);
int printf_fork ( FILE *fp,char *string,...);

int  printf_system (char *string, ...);
char*printf_system2string (char *string, ...);
int  printf_system_direct (char *string, ...);
int  printf_system_direct_check (char *string, ...);

int my_system_cl (int argc, char *argv[]);
int my_system ( char *command);
int unpack_perl_script (char *name, char ***unpacked, int n);
void unpack_all_perl_script (char *script);
/*********************************************************************/
/*                                                                   */
/*                         IO FUNCTIONS                              */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
void set_command_line (char *s);
FILE * print_command_line (FILE *fp );
int getpid_ref();
char ** standard_initialisation       ( char **in_argv, int *in_argc);
char ** standard_initialisation_start ( char **in_argv, int *in_argc);
char ** standard_initialisation_end   ( char **in_argv, int *in_argc);
/*
by default : dir_4_tcoffee: $HOME/.t_coffee
tmp: dir_4_tcoffee/tmp OR TMP_4_TCOFFEE
cache: idem
methods: idem
mcoffee: idem
*/
int set_nproc (int nproc);
int get_nproc ();
char *get_os();
char *get_lockdir_4_tcoffee ();
char *get_plugins_4_tcoffee ();
char *get_home_4_tcoffee();
char *get_dir_4_tcoffee();
char *get_tmp_4_tcoffee();
char *get_cache_4_tcoffee();
char *get_methods_4_tcoffee();
char *get_mcoffee_4_tcoffee();



void myexit (int signal);
FILE *fatal_exit ( FILE *fp, int exit_signal, char *string, ...);
int set_warning_mode ( int mode);
FILE *add_warning (FILE *fp, char *string, ...);
FILE *add_information (FILE *fp, char *string, ...);
int  fprintf_error( FILE *fp, char *string, ...);

void output_warning_list();

int   count_n_res_in_array  (char *array, int len);
int   count_n_gap_in_array  (char *array, int len);
int count_n_symbol_in_array ( char *array, char *array_list, int len);
char* count_strings_in_file ( char *in, char *out);
char** count_strings ( char **array, int len);
int ** count_int_strings ( int **array, int len, int s);

int   get_first_non_white_char (char *name);
int   count_n_char_x_in_file(char *name, char x);
int   count_n_char_in_file(char *name);
int   count_n_line_in_file(char *name);
int measure_longest_line_in_file ( char *name );
int file_cat ( char *fname1, char *fname2);
FILE* display_file_content (FILE *output, char *name);
int   cat_file (char *file1, char *file2);

char **list2expanded_flist (char **list, int *n, char *tag);
char **expand_flist (char *file, char **list,int i,int *n, char *tag);

char ***file2list   (char *name, char *sep);
char ** file2lines  (char *name);
int     file2nlines (char *name);
char *  file2string (char *name);
char    file2lastchar (char *name);
char    file2firstchar (char *name);

int     file2size   (char *name);
float **file2float_mat (char * name, char *sep);
int string2file_direct (char *file, char *mode, char *string,...);
int string2file ( char *file, char *mode, char *string,...);
char *chomp (char *name);
int get_cl_param (int argc, char **argv, FILE **fp, const char para_name[], int *set_flag, const char type[], int optional, int max_n_val,const char usage[], ...);
char ** get_parameter ( char *para_name, int *np, char *fname);

char *get_t_coffee_environement (char *file);
char *set_path_4_plugins (char *);
int add_package2_tcoffee_env (char *package);

char *dump_io_start (char *dump);
void  dump_io (char *dump_file, char *dump_nature);
void  dump_error_file();
void  update_error_dir();

FILE* stack_msg(FILE *fp);
FILE* install_msg(FILE *fp);
FILE* proxy_msg(FILE *fp);
FILE* email_msg(FILE *fp);
FILE* error_msg(FILE *fp);
FILE *proxy_msg(FILE *fp);

char *get_proxy();
char *get_proxy_from_env();
int set_proxy (char *proxy);

char *input_name ();
char *Email4cl(int input_mode, int set_mode);
char *Email(int input_mode, int set_mode);
char *input_email ();
char *get_email_from_env ();
char *get_email ();
int set_email (char *email);
int cputenv4pathFirst (char*);
int cputenv4pathLast (char*);
int cputenv4path (char*);

int   cputenv (char*, ...);
int   fcputenv (char *,char *, char*, ...);
char *file_putenv (char *file);
int check_dir_getenv ( char *string);

char* set_string_variable (char *key, char* v);
char* get_string_variable (char *key);
char* unset_string_variable (char *var);
char* store_string_variable (char *var, char * v, int mode);

int int_variable_isset (const char var[]);
int set_int_variable (char *var, int v);
int get_int_variable (const char var[]);
int unset_int_variable (const char var[]);
int store_int_variable (char *var, int v, int mode);

void check_vtmpnam ();
int flag_file2remove_is_on ();
void set_file2remove_off();
void set_file2remove_on();
char *set_file2remove_extension(char *extension, int mode);
char * add2file2remove_list ( char *name);


FILE *get_stdout1(char *name);
FILE * vtmpfile();
void initiate_vtmpnam (char *s);
int     vtmpnam_size ();
int     isvtmpnam( char *s);
char *  vtmpnam ( char *s);
char *  vtmpnamH( );
char *alp2random_string (char*s);

char *  tmpnam_2 (char *s);
void    safe_remove(char*s);
char *  vremove ( char *s);
char *  vremove2 ( char *s);
char *  vremove3 ( char *s, char *ext);

void error_exit_sigsegv (int x);
void error_exit_sigill (int x);
void error_exit_sigfpe (int x);
void error_exit_sigabrt(int x);

void error_exit (int);
void clean_exit();

void signal_exit_sigterm (int x);
void signal_exit_sigint  (int x);
void signal_exit_sigill (int x);
void signal_exit_sigsegv (int x);

void signal_exit(int);

void main_exit ();
int  log_function (char *fname);
   
void   clean_function ( );
void   sig_clean_function ( int x);
char * prepare_cache ( const char *mode);
char * get_cache_dir();
void update_cache ();
void ignore_cache();


int  register_file4dump (char *name, char *mode);
char *capture_stdin ();
int count_openF();
void valgrind_test();
FILE * vfopen ( char *name, char *mode);
FILE * vfclose (FILE *fp);
int echo ( char *string, char *fname);

int **get_file_block_pattern (char *fname, int *n_blocks, int max_n_line);

int token_is_in_file_n ( char *fname, char *token, int nlines);
int token_is_in_file (char *fname, char *token);

int    check_file_for_token      ( char *file , char *token);
int token_is_in_n_lines ( char *fname, char *token, int n_line);
FILE * find_token_in_file_nlines ( char *fname, FILE * fp, char *token, int n_line);
FILE * find_token_in_file ( char *fname, FILE * fp, char *token);
FILE * quick_find_token_in_file  (FILE *fp, char *token);

char * vfgets (char *buf, FILE *fp);

FILE * set_fp_after_char ( FILE *fp, char x);
FILE * set_fp_id ( FILE *fp, char *id);
FILE * skip_commentary_line_in_file ( char com, FILE *fp);
char * strip_file_from_comments (char *com, char *in_file);

int check_for_update ( char *web_address);
int url2file (char *address, char *out);
int wget (char *address, char *out);
int curl (char *address, char *out);
	

int simple_check_internet_connection  (char *address);
int check_internet_connection  (int mode);
int check_environement_variable_is_set ( char *variable, char *description, int fatal);
int check_program_is_installed ( char *program_name, char *current_path, char *path_variable, char *where2getit, int fatal);
FILE * display_output_filename ( FILE *io, char *type, char *format, char *name, int check_output);
FILE * display_input_filename ( FILE *io, char *type, char *format, char *name, int check_output);
int filename_is_special ( char *fname);
char *check_url_exists ( char *fname);
char *check_file_exists ( char *fname);
int my_mkdir ( char *dir);
int my_rmdir ( char *dir);

int file_is_empty(char *fname);
int file_exists (char *path,char *fname);
int isfile(char *fname);
int isexec (char *fname);
int istmp  (char *name);
int isdir  (char *fname);
int iswdir  (char *fname);

int isdir4path (char *fname);
int rrmdir (char *fname);
char * ls_l(char *path,char *fname);

void create_file ( char *name);
void delete_file ( char *fname);
int  util_rename ( char* from, char *to);
int  util_copy   ( char* from, char *to);
void reset_output_completion ();
FILE * output_completion4halfmat ( FILE *fp,int n, int tot, int n_eports, char *s);
FILE * output_completion ( FILE *fp,int n, int tot, int n_eports, char *s);
void * null_function (int a, ...);
int  btoi ( int nc,...);
/*********************************************************************/
/*                                                                   */
/*                         Geometric FUNCTIONS                    */
/*                                                                   */
/*                                                                   */
/*********************************************************************/

float get_geometric_distance ( float ** matrix, int ncoor, int d1, int d2, char *mode);
/*********************************************************************/
/*                                                                   */
/*                         MATHEMATICAL FUNCTIONS                    */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
double log_addN ( int N, double *L);
double log_add6 (double a, double b, double c, double d, double e, double f );
double log_add5 (double a, double b, double c, double d, double e);
double log_add4 (double a, double b, double c, double d);
double log_add3 (double a, double b, double c);
double log_add2 (double a, double b);

float factorial_log ( int start, int end);
float M_chooses_Nlog ( int m, int N);
double factorial ( int start, int end);
double M_chooses_N ( int m, int N);
float my_int_log(int a);
/*********************************************************************/
/*                                                                   */
/*                         Fast Log Additions (adapted from Probcons)*/
/*                                                                   */
/*                                                                   */
/*********************************************************************/
double EXP(double x);
float LOOKUP (float x);
void LOG_PLUS_EQUALS (float *x, float y);
float LOG_ADD (float x, float y);
float LOG_ADD3 (float x1, float x2, float x3);
float LOG_ADD4 (float x1, float x2, float x3, float x4);
float LOG_ADD5 (float x1, float x2, float x3, float x4, float x5);
float LOG_ADD6 (float x1, float x2, float x3, float x4, float x5, float x6);
float LOG_ADD7 (float x1, float x2, float x3, float x4, float x5, float x6, float x7);
///////////////////////////////////////////////////////////////////////////////////////////
// Hash function
////////////////////////////////////////////////////////////////////////////////////////////
unsigned long hash_file(char* file);  //returns the hash value for key 
int file2diff (char *file1, char *file2);
///////////////////////////////////////////////////////////////////////////////////////////
// Generating lists through recirsive exploration
////////////////////////////////////////////////////////////////////////////////////////////
int **generate_array_int_list (int len, int min, int max, int step, int *n, char *filename);
char ***generate_array_string_list (int len, char ***alp, int *alp_size, int *n, char *file, int mode);
float *display_accuracy (float *count, FILE *fp);
float *counts2accuracy (float *count);

float rates2sensitivity (int tp, int tn, int fp, int fn, float *sp, float *sn, float *sen2, float *best);
float profile2sensitivity (char *pred, char *ref, float *sp, float *sn, float *sen2, float *b);	      
float profile2evalue (char *pred, char *ref);
//isexec lib

///////////////////////////////////////////////////////////////////////////////
//                                                                           //  
//                                                                           //  
//                           Kmeans                                          //  
///////////////////////////////////////////////////////////////////////////////



void      km_output_data ( double **data, int n, int dim, int len,  char *infile, char *outfile);
int       km_kmeans(double **data, int n, int m, int k, double t, double **centroids);
float     km_make_kmeans(double **d, int n, int dim, int k,double t, double **centroids, int nrounds);
double ** km_read_data ( char *file, int n, int dim, int len, char **field_list);
int       km_file2dim     ( char *file, int *n,int *dim, int *len);
double**  km_shuffle_data (double **d, double **sd, int n, int r);
void      km_display_data (double **d, int n, int dim);
double    km_data2evaluate ( double **d, int n, int dim);
double km_kmeans_bs (double **data, int n, int dim, int k,double t, double **centroids, int nrounds);

int mat2process (int n, char *flist[]);
