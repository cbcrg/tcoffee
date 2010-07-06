#include <time.h>
#include <signal.h>
#include <errno.h>
#include <sys/times.h>
#include <sys/types.h>
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
int cmp_float ( const float **a, const float **b);
void sort_int_1D ( int *L, int n);
char** sort_string_array (char **V, int n);
     
void sort_int ( int **V,int N_F, int F, int left, int right);
void sort_list_int ( int **V,int N_F, int F, int left, int right);
void sort_list_int2 ( int **V,int *list,int N_F, int left, int right);
void sort_int_inv ( int **V,int N_F, int F, int left, int right);
void sort_list_int_inv ( int **V,int N_F, int F, int left, int right);
int cmp_int ( const int**a, const int**b);
int cmp_list_int (const int**a, const int**b);
int cmp_list_int2 (const int**a, const int**b);

int name_is_in_list ( char *name, char **name_list, int n_name, int len);
char * check_list_for_dup ( char **list, int ne);
FILE *get_number_list_in_file ( FILE *fp, int *list, int *n, int *max_len);


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
int atoigetenv(const char *x);
void setenv_func ( char *string_name, char *string_value);
void get_pwd ( char *name);
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
char* strcatf  (char *string1,char *string2, ...);
char *vcat (char *v1, char *v2);

int strget_param ( char *string, char *param_name, char *param_value, char *format, ...);
char * lstrstr ( char *in, char *token);
char * vstrstr ( char *in, char *token);
int strscanf (char *in, char *token, char *format, ...);
int match_motif ( char *string, char **motif);

char *after_strstr (char *string, char *token);

char ** push_string (char *val, char **stack, int *nval, int mode);
int vsrand (int val);
int  *randomize_list (int *list, int len, int ncycle);
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
int *  string2num_list( char *string);
int *  string2num_list2( char *string, char *separators);
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
char *lock (int pid, int type, int action, char *value, ...);
int check_process (const char *com,int pid,int r, int failure_handling);
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
void HERE (char *string, ...);
void printf_exit  (int exit_code, FILE *fp, char *string, ...);
int printf_file ( char *file, char *mode, char *string, ...);
int printf_fork ( FILE *fp,char *string,...);
int printf_system (char *string, ...);
int printf_system_direct (char *string, ...);
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
char ***file2list   (char *name, char *sep);
char ** file2lines  (char *name);
char *  file2string (char *name);
int string2file ( char *file, char *mode, char *string,...);
char *chomp (char *name);
int get_cl_param (int argc, char **argv, FILE **fp,char *para_name, int *set_flag, char *type, int optional, int max_n_val,char *usage, ...);
char ** get_parameter ( char *para_name, int *np, char *fname);

char *get_t_coffee_environement (char *file);
char *set_path_4_plugins (char *);
int add_package2_tcoffee_env (char *package);


void dump_error_file();
void update_error_dir();

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
int cputenv4path (char*);
int   cputenv (char*, ...);
int   fcputenv (char *,char *, char*, ...);
char *file_putenv (char *file);
int check_dir_getenv ( char *string);

char* set_string_variable (char *var, char* v);
char* get_string_variable (char *var);
char* unset_string_variable (char *var);
char* store_string_variable (char *var, char * v, int mode);

int int_variable_isset (char *var);
int set_int_variable (char *var, int v);
int get_int_variable (char *var);
int unset_int_variable (char *var);
int store_int_variable (char *var, int v, int mode);

void check_vtmpnam ();
int flag_file2remove_is_on ();
void set_file2remove_off();
void set_file2remove_on();
char *set_file2remove_extension(char *extension, int mode);
char * add2file2remove_list ( char *name);



FILE * vtmpfile();
void initiate_vtmpnam (char *s);
char * vtmpnam ( char *s);
char *  tmpnam_2 (char *s);
void    safe_remove(char*s);
char *  vremove ( char *s);
char *  vremove2 ( char *s);
void error_exit ();
void clean_exit();
void signal_exit();
void main_exit ();
int  log_function (char *fname);
   
void   clean_function ( );
void   sig_clean_function ( int x);
char * prepare_cache ( const char *mode);
char * get_cache_dir();
void update_cache ();
void ignore_cache();

FILE * vfopen ( char *name, char *mode);
FILE * vfclose (FILE *fp);
int echo ( char *string, char *fname);

int **get_file_block_pattern (char *fname, int *n_blocks, int max_n_line);

int token_is_in_file (char *fname, char *token);

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
char *check_file_exists ( char *fname);
int my_mkdir ( char *dir);
int file_is_empty(char *fname);
int file_exists (char *path,char *fname);
int isexec (char *fname);
int isdir  (char *fname);
int isdir4path (char *fname);
int rrmdir (char *fname);
char * ls_l(char *path,char *fname);

void create_file ( char *name);
void delete_file ( char *fname);
int  util_rename ( char* from, char *to);
int  util_copy   ( char* from, char *to);
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
#include <stdio.h>

unsigned long linrand(unsigned long r);
unsigned long addrand(unsigned long r);
void addrandinit(unsigned long s);

unsigned long mult(unsigned long p,unsigned long q);
    

struct Job_TC
    {
      int jobid;
      int status;
      
      struct  Job_TC *c;
      struct  Job_TC *p;
      struct  Job_io_TC *io;
      struct  Job_control_TC *control;
      
      struct  Job_param_TC *param;
        
      /*memory mangement*/
      char **pl;
      int np;
};
typedef struct Job_TC Job_TC;

struct Job_control_TC
    {
      
      struct Job_TC* (*submitF) (struct Job_TC*);
      struct Job_TC* (*retrieveF)(struct Job_TC*);
      char *mode;
};
typedef struct Job_control_TC Job_control_TC;

struct Job_io_TC
    {
      char *in;
      char *out;
      struct Constraint_list *CL;
      struct Alignment *A;
};
typedef struct Job_io_TC Job_io_TC;

struct Job_param_TC
{
  char *method;
  struct TC_method *TCM; 
  char *temp_c;
  char *aln_c;
  char *seq_c;
  char *aln_mode;
};
typedef struct Job_param_TC Job_param_TC;

Job_TC* print_lib_job ( Job_TC *job,char *string, ...);
Job_TC *print_lib_job2 ( Job_TC* job, int n, char **name, char **value);


/*Stack Manipulation*/
Job_TC *free_queue  (Job_TC *job);
Job_TC *free_job  (Job_TC *job);
Job_TC * queue2heap (Job_TC*job);
Job_TC * queue2last (Job_TC*job);
int queue2n (Job_TC*job);
Job_TC * descend_queue (Job_TC*job);
Job_TC *queue_cat  (Job_TC *P, Job_TC *C);
Job_TC *delete_job (Job_TC *job);
/*Job Control*/
struct Job_TC* submit_job ( Job_TC *job);
struct Job_TC* retrieve_job ( Job_TC *job);
Job_TC*** split_job_list (Job_TC *job, int ns);
struct Dps_result
    {
      int njobs;
      struct Dps_job **dps_job;
};
typedef struct Dps_result Dps_result;

struct Dps_job
    {
      int JobId;
      struct Constraint_list *CL;
      char *input_file;
      char *output_file;
};
typedef struct Dps_job Dps_job;

struct Dps_result *seq2list_DPS (struct Constraint_list *CL,char *method, char *aln_command, char *seq_command, char *weight, Dps_result *dps_result);
struct Constraint_list * gather_results_DPS ( Dps_result *DPS, struct Constraint_list *CL);
Dps_result *declare_dps_result ( int naln, Dps_result *dps);

#define SEQ2 0
#define R2   1
#define WE   2
#define CONS 3
#define MISC 4
#define SEQ1 5
#define R1   6
#define INDEX 7

#define ICHUNK 5

#define LIST_N_FIELDS 7
#define CLIST_TYPE int

/*********************************************************************************************/
/*                                                                                           */
/*         FUNCTIONS Typedef                                                                 */
/*                                                                                           */
/*********************************************************************************************/
typedef int (*Profile_cost_func) (int*, int *,struct Constraint_list *);
typedef int (*Col_cost_func)(Alignment*, int **, int, int*, int, int **, int, int*, int, struct Constraint_list *);
typedef int (*Pair_cost_func)(struct Constraint_list *, int, int, int, int);
typedef int (*Pwfunc) (Alignment *, int*, int **,struct Constraint_list *);

/*********************************************************************************************/
/*                                                                                           */
/*         STRUCTURES FOR PDB ANALYSIS                                                       */
/*                                                                                           */
/*********************************************************************************************/
typedef struct 
    {
      char *use_seqan;
}
TC_param;
typedef struct 
    {
      char blast_server[FILENAMELEN+1];
      char db[FILENAMELEN+1];
      int min_cov;
      int min_id;
      int max_id;
}
Blast_param;

typedef struct
    {
	int   n_excluded_nb;
	
	float similarity_threshold;
	float rmsd_threshold;
        float md_threshold;
        int   distance_on_request;
	char  *comparison_io;
        int    print_rapdb;
        float maximum_distance;/*Diameter of the bubble used to identify the Calpha Neighborhood*/
	int   N_ca;            /*Number of Calpha to be looked at on both side*/
	float max_delta ;      /*Maximum value for delta to be positive*/ 
	char *local_mode;
        int   scale;             /*Value substracted to the pdb score in the bubble mode*/
        int   n_extra_param;
        char **extra_param;
      char  *evaluate_mode;
      char  *color_mode;
      float filter;
      int filter_aln;
      int irmsd_graph;
      int nirmsd_graph;
      
      
    }
Pdb_param;

typedef struct
    {
	int num;
	int res_num;/*Residue number from 1 to N*/
        char res[4];
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
struct Distance_matrix
{
  char mode[100];
  char sim_mode[100];
  char nseq;
  int     **similarity_matrix; /*Pairwise ID levels: 1-10000*/ 
  int     **score_similarity_matrix; /*Pairwise ID levels: 1-10000*/ 
  int     **distance_matrix; /*Pairwise ID levels: 1-10000*/ 
};
typedef struct Distance_matrix Distance_matrix;
struct Constraint_list
    {
      /*In Case of Modif, synchronize with:
	util_declare/declare_constraint_list
	util_declare/cache_dp_value4constraint_list
	util_declare/duplicate_constraint_list
	util_declare/free_constraint_list
      */

      //Generic parameters
      TC_param *TC;

      int copy_mode;
      struct Constraint_list *pCL; 
      Sequence *S;         /*Total sequences*/
      Sequence *STRUC_LIST; /*Name of the sequences with a Structure*/
      char align_pdb_param_file[FILENAMELEN+1];
      char align_pdb_hasch_mode[FILENAMELEN+1];
     

      Weights  *W;         /*Sequence Weights*/
      Distance_matrix *DM; /*Accurate Distance Matrix*/
      Distance_matrix *ktupDM; /*Fast Distance Matrix*/
      Fname *RunName;
      
      int *translation;   
      char **  out_aln_format;
      int    n_out_aln_format;

      
      /*Packing Sequence: To use with domain analysis*/
      int **packed_seq_lu;
      
      /*DATA*/
      FILE *fp;           /*File used for i/o if disk being used*/
      //int *L;            /*Array used for storing Lib if mem being used*/
      int **M;            /*substitution matrix*/
      char rna_lib[FILENAMELEN+1];  /*name of a file containing the RNA libraries*/
      
      /*List Information*/      
      int ne;             /*Number of elements in the list*/
      char list_name[1000];    /*Name of the list*/
      int  entry_len;     /*Size of an entry in el_size*/
      size_t el_size;     /*Size of each elements in an entry in bytes*/
      
      /*Normalisation information*/
      int normalise;
      int max_ext_value;
      int max_value;
      int overweight;
      int filter_lib;
      
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

      char  dp_mode[FILENAMELEN+1];
      int   reverse_seq;
      int   maximise;
      char  matrix_for_aa_group[FILENAMELEN+1];
      char  method_matrix[FILENAMELEN+1];
      float diagonal_threshold;
      int ktup;
      int use_fragments;
      int fasta_step;
      int lalign_n_top;
      int sw_min_dist;
      char **matrices_list;
      int n_matrices;
      char tree_mode[FILENAMELEN+1];

      char distance_matrix_mode[FILENAMELEN+1];
      char distance_matrix_sim_mode[FILENAMELEN+1];
      
      Alignment *tree_aln;
      
      /*Functions used for dynamic programming and Evaluation*/
      int no_overaln;
      /*1 Function for evaluating the cost of a column*/
      Col_cost_func get_dp_cost;
      Profile_cost_func profile_mode;
      char profile_comparison [FILENAMELEN+1];
      
      /*2 Function for evaluating the cost of a pair of residues*/
      Pair_cost_func evaluate_residue_pair;
      /*3 Function for making dynamic programming*/
      Pwfunc pair_wise;
      
      /*
      int (*get_dp_cost)(Alignment*, int **, int, int*, int, int **, int, int*, int, struct Constraint_list *);
      int (*evaluate_residue_pair)(struct Constraint_list *, int, int, int, int);
      int (*pair_wise)(Alignment *, int*, int **,struct Constraint_list *);
      */

      int weight_field;
      int max_n_pair; /*maximum number of pairs when aligning two profiles*/

      /*Extend a sequence against itself*/
      
      /*Threading parameters*/
      Blast_param *Prot_Blast;
      Blast_param *Pdb_Blast;
      Blast_param *DNA_Blast;
      /*Split parameters*/
      int split;
      int split_nseq_thres;
      int split_score_thres;
      /*Check Structural Status*/
      int check_pdb_status;
      /*log*/
      char method_log[1000];
      char evaluate_mode[1000];
      char method_evaluate_mode[100];
      /*Gene Prediction Parameter*/      
      char genepred_score[100];
      
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
      int ** freeze;
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
      char  multi_thread[100];
      char  lib_list[FILENAMELEN+1];
};

typedef struct Constraint_list Constraint_list;

struct TC_method 
{
  
  char executable[FILENAMELEN+1];
  char executable2[FILENAMELEN+1];
  char in_flag[FILENAMELEN+1];
  char in_flag2[FILENAMELEN+1];
  char out_flag[FILENAMELEN+1];
  char aln_mode[FILENAMELEN+1];
  char out_mode[FILENAMELEN+1];
  char seq_type[FILENAMELEN+1];
  char weight[FILENAMELEN+1];
  char matrix[FILENAMELEN+1];
  int gop;
  int gep;
  int minid;
  int maxid;
  char param[1000];
  char param1[1000];
  char param2[1000];
  
  Constraint_list *PW_CL;
};
typedef struct TC_method TC_method;

/*********************************************************************/
/*                                                                   */
/*                         PRODUCE IN LIST                             */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
Constraint_list *produce_list ( Constraint_list *CL, Sequence *S, char * method,char *weight,char *mem_mode);	
Job_TC* method2job_list ( char *method, Sequence *S,char *weight, char *lib_list, Distance_matrix *DM, Constraint_list *CL);

Job_TC *job_list2multi_thread_job_list (Job_TC* ojob, char *mt, Constraint_list *CL);
Job_TC *retrieve_lib_job ( Job_TC *job);
Job_TC *submit_lib_job ( Job_TC *job);
int add_method_output2method_log (char *l, char *command,Alignment *A, Constraint_list *CL, char *iofile);

int check_seq_type (TC_method *M, char *slist,Sequence *S);
int check_profile_seq_type (Sequence *S, int i, char t);
char **method_list2method4dna_list ( char **list, int n);
int is_in_pre_set_method_list (char *fname);
char *** display_method_names (char *mode, FILE *fp);

char *method_name2method_file (char *method);
char *make_aln_command(TC_method *m, char *seq, char *aln);
struct TC_method* method_file2TC_method ( char *fname);
char *method_file_tag2value (char *method, char *tag);
int TC_method2method_file( struct TC_method*, char *fname );
/*********************************************************************/
/*                                                                   */
/*                         WRITE IN LIST                             */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
Constraint_list *  unfreeze_constraint_list (Constraint_list *CL);
Constraint_list *  freeze_constraint_list (Constraint_list *CL);
Constraint_list *  undump_constraint_list (Constraint_list *CL, char *file);
int   dump_constraint_list (Constraint_list *CL, char *file,char *mode);
int   safe_dump_constraint_list (Constraint_list *CL, char *file,char *mode, Sequence *RS);
int display_constraint_list (Constraint_list *CL, FILE *fp, char *tag);


Constraint_list *index_constraint_list ( Constraint_list *CL);
Constraint_list *index_res_constraint_list ( Constraint_list *CL, int field);
Constraint_list * progressive_index_res_constraint_list ( Alignment *A, int *ns, int **ls, Constraint_list *CL);
char ** reindex_constraint_list (char **profile, int np,char **list, int *inL, Sequence *S);
/*********************************************************************/
/*                                                                   */
/*                         ENTRY MANIPULATION                        */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
Constraint_list * add_list_entry2list (Constraint_list *CL, int n_para, ...);
Constraint_list * evaluate_constraint_list_reference ( Constraint_list *CL);

int CLisCompacted (Constraint_list *CL, char *t);
int checkCL( Constraint_list *CL, char *t);
Constraint_list *add_entry2list ( CLIST_TYPE *entry, Constraint_list *CL);
Constraint_list *add_entry2list2 ( CLIST_TYPE *entry, Constraint_list *CL);
int *extract_entry (Constraint_list *CL);
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
//CLIST_TYPE *search_in_list_constraint(int *key, int k_len, int **L, int ne, int ***start_index, int ***end_index);
CLIST_TYPE *main_search_in_list_constraint ( int *key,int *p,int k_len,Constraint_list *CL);
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
Constraint_list* read_n_constraint_list(char **fname,int n_list, char *in_mode,char *mem_mode,char *weight_mode,char *type, FILE *local_stderr, Constraint_list *CL, char *seq_source);
Constraint_list* read_constraint_list(Constraint_list *CL,char *fname,char *in_mode,char *mem_mode,char *weight_mode);
Constraint_list * read_constraint_list_raw_file(Constraint_list *CL, char *fname);

int        read_cpu_in_n_list(char **fname, int n);
int read_seq_in_list ( char *fname,  int *nseq, char ***sequences, char ***seq_name);

Sequence * read_seq_in_n_list(char **fname, int n, char *type, char *SeqMode);

int        read_cpu_in_list ( char *fname);
int ** read_list ( char *fname, int **list,int *ne, int *nseq, int *cpu, char ***sequences, char ***seq_name);

char * expand_constraint_list_file ( char *file);
Constraint_list * read_constraint_list_file(Constraint_list *CL, char *fname);
Constraint_list * fast_read_constraint_list_file(Constraint_list *CL, char *fname);

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
int constraint_list2raw_file ( Constraint_list *CL, char *fname, char *fmode);
FILE * save_raw_constraint_list   ( FILE *fp,Constraint_list *CL, int start,int len, int *translation);
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
Constraint_list * shrink_constraint_list_indexed (Constraint_list *CL, int T);
Constraint_list * shrink_constraint_list (Constraint_list *CL);
Constraint_list * extend_constraint_list ( Constraint_list *CL);
Constraint_list * relax_constraint_list (Constraint_list *CL);
Constraint_list * relax_constraint_list_4gp (Constraint_list *CL);

Constraint_list * expand_constraint_list_4gp (Constraint_list *CL, int T);

Constraint_list * filter_constraint_list (Constraint_list *CL, int field, int T);
int constraint_list_is_connected ( Constraint_list *CL);
int constraint_list2avg ( Constraint_list *CL);
float constraint_list2connectivity ( Constraint_list *CL);

int constraint_list2fraction_covered ( Constraint_list *CL);

int *seqpair2weight (int s1, int s2, Alignment *A,Constraint_list *CL, char *weight_mode, int *weight);
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
Constraint_list *merge_constraint_list   ( Constraint_list *SL, Constraint_list *ML, char *mode);
CLIST_TYPE return_max_constraint_list ( Constraint_list *CL, int field);
Constraint_list *modify_weight( Constraint_list *CL,int start, int end,  char *modify_mode);
Constraint_list *compact_list (Constraint_list *CL, char *compact_mode);
Constraint_list *rescale_list_simple (Constraint_list *CL,int start, int len,int new_min, int new_max);
Constraint_list *rescale_list (Constraint_list *CL,int start, int len,int max1, int max2);
Constraint_list* filter_list (Constraint_list *CL, int start, int len,int T);
Constraint_list *undefine_list (Constraint_list *CL);
int ** seq2defined_residues ( Sequence *S, Constraint_list *CL);
int ** aln2defined_residues ( Alignment *A, Constraint_list *CL);
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
/*         FUNCTION FOR PRUNING THE LIST                                                   */
/*                                                                                           */
/*********************************************************************************************/
char * list2prune_list (Sequence *S, int **sm);
/*********************************************************************************************/
/*                                                                                           */
/*         FUNCTION FOR WEIGHTING THE LIST                                                   */
/*                                                                                           */
/*********************************************************************************************/
Constraint_list *weight_constraint_list(Constraint_list * CL, char *seq_weight);
Weights* compute_t_coffee_weight(Constraint_list * CL);
Constraint_list *re_weight_constraint_list(Constraint_list * CL,Weights *W);
Constraint_list *set_weight4constraint_list(Constraint_list * CL,int w);

Distance_matrix *cl2distance_matrix (Constraint_list *CL, Alignment *A,  char *mode, char *sim_mode, int print);
Distance_matrix *seq2distance_matrix (Constraint_list *CL, Alignment *A,  char *mode, char *sim_mode, int print);

/*********************************************************************************************/
/*                                                                                           */
/*         MULTI_THREAD                                                                      */
/*                                                                                           */
/*********************************************************************************************/
int run_multi_thread_file (char *fname, char *config);
/*********************************************************************/
/*                                                                   */
/*                        RNA FUNCTIONS                              */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
char * seq2rna_lib ( Sequence *S, char *name);
Constraint_list *read_rna_lib ( Sequence *S, char *fname);
Constraint_list *rna_lib_extension ( Constraint_list *CL, Constraint_list *R);
char *** produce_method_file ( char *method);
typedef struct
    {
      int   p1;
      int   p2;
      int   p3;
      int   p4;
      int   t;
      int   f;
      char  mode[20];//lower, unalign
      char  model[20];//fsa1 fsa2
}
OveralnP;

//RNA

int **    alifold_list2cov_list (Alignment *A, int **list);
int ** update_RNAfold_list (Alignment *A, int **pos, int s, int **l);
int ** vienna2list ( char *seq);
Alignment *compare_RNA_fold ( Alignment *A, Alignment *B);

Alignment *alifold2analyze (Alignment *A, Alignment *ST, char *mode);
Alignment *alifold2cov_aln (Alignment *A, int **l, int ug);
Alignment *alifold2cov_stat (Alignment *A, int **l, int ug);
Alignment *alifold2cov_list (Alignment *A, int **l, int ug);
Alignment *alifold2cov_cache (Alignment *inA,  int **l, int ug);


Alignment *add_alifold2aln  (Alignment *A, Alignment *ST);
Alignment *aln2alifold(Alignment *A);

//end
Alignment * aln2bootstrap (Alignment *A, int n);
Alignment * aln2sample    (Alignment *A, int n);
Alignment * aln2random_aln (Alignment *A, char *mode);
Alignment *aln2scale (Alignment *A, char *offset);
Alignment* aln2case_aln (Alignment *A, char *upper, char *lower);
Alignment*aln2gap_cache (Alignment *A, int val);
Alignment *score_aln2score_ascii_aln (Alignment *A, Alignment *C);
int **aln2resindex ( Alignment *A, Alignment *B, FILE *fp);
int **index_seq_res      ( Sequence *S1, Sequence *S2, int **name_index);
int **index_seq_name ( Sequence *S1, Sequence *S2);
int *get_name_index (char **l1, int n1, char **l2, int n2);

int* get_res_index (char *seq1, char *seq2);
int * pos2list (int * pos, int len, int *nl);
  int *list2pos (int *list, int nl, int len);


int change_residue_coordinate ( char *in_seq1, char *in_seq2, int v);

int ** minimise_repeat_coor (int **coor, int nseq, Sequence *S);
int ** get_nol_seq( Constraint_list *CL,int **coor, int nseq, Sequence *S);


int compare_pos_column( int **pos1,int p1, int **pos2,int p2, int nseq);



char * seq2alphabet (Sequence *S);
char *aln2alphabet (Alignment *A);
char *array2alphabet (char **array, int n, char *forbiden);

//TM Predictions
char* alnpos2hmmtop_pred (Alignment *A, Alignment *Pred, int pos, int mode);
Alignment * aln2hmmtop_pred (Alignment *A);
char * seq2tmstruc ( char *seq);

char * set_blast_default_values();
char      *  seq2pdb   ( Sequence *S);
Alignment *  seq2blast ( Sequence *S);

Sequence * seq2unique_name_seq ( Sequence *S);
Alignment * aln2unique_name_aln ( Alignment *S);
int name_list2unique_name_list (int n, char **name);
Sequence *seq2clean_seq ( Sequence *S, char *alp);//remove all alp characters from seq
char**gene2exons    (char **seq, int nseq);
 
int       ** seq2aln_pos      (Alignment *A, int *n, int **ls);
Alignment *padd_aln ( Alignment *A);
char **padd_string ( char **string, int n,char pad);

Alignment *local_maln2global_maln (char *seq, Alignment *A);

Alignment * seq2profile (Sequence *S, int index);

Sequence *remove_empty_sequence (Sequence *S);
Alignment *  aln2profile (Alignment * A);
Alignment * aln2collapsed_aln (Alignment * A, int n, char **string);
Alignment* remove_seq_from_aln (Alignment *A, char *seq);

Alignment* aln2sub_aln_file (Alignment *A, int n, char **string);
Alignment* aln2sub_seq (Alignment *A, int n, char **string);

int       ** aln2inv_pos  (Alignment *A);
int        * seq2inv_pos ( char *seq);
int       ** aln2pos_simple   (Alignment *A, int n_nseq, ...);
int       ** aln2pos_simple_2 (Alignment *A);
Alignment ** split_seq_in_aln_list ( Alignment **aln, Sequence *S, int l_seq, char **seq_list);

Sequence * fill_sequence_struc ( int nseq, char **sequences, char **seq_name);

int  seq_list2in_file ( TC_method *M, Sequence *S, char *list, char *file);
int  seq_list2fasta_file( Sequence *S,  char *list, char *file, char *outmode);
Structure * seq2struc    ( Sequence *S, Structure *ST);
Alignment *strings2aln (int nseq,...);

Alignment * seq2aln      ( Sequence *S, Alignment *A,int rm_gap);
Alignment  *seq_coor2aln ( Sequence *S, Alignment *A, int **coor, int nseq);

Alignment *stack_aln (Alignment *A, Alignment *B);
Alignment *chseqIaln(char *name, int seq_n, int start,int len,Sequence *S, int seqIaln, Alignment *A);


char *dna_aln2cons_seq ( Alignment *A);
char *aln2cons_seq ( Alignment *A, int ns, int *ls, int n_groups, char **group_list);
char *aln2cons_maj ( Alignment *A, int ns, int *ls, int n_groups, char **group_list);
Alignment *aln2conservation ( Alignment *A, int threshold,char *seq);

char *sub_aln2cons_seq_mat ( Alignment *A,int ns, int *ls, char *mat_name);
char *aln2cons_seq_mat ( Alignment*A, char *mat_name);
Alignment *aln2short_aln( Alignment *A, char *list, char *new, int spacer);
Sequence  *keep_residues_in_seq ( Sequence *S,char *list, char replacement);
Alignment *keep_residues_in_aln ( Alignment *A,char *list, char replacement);
Alignment *filter_keep_residues_in_aln ( Alignment *A,Alignment *ST, int use_cons, int value, char *list, char replacement);

Alignment *aln_convert (Alignment *A, Alignment *ST, int use_cons, int value,int n, ...);
Alignment *aln2number (Alignment *A);
Alignment * filter_aln ( Alignment *A, Alignment *ST, int value);
Alignment * filter_aln_lower_upper ( Alignment *A, Alignment *ST,int use_cons, int value);
Alignment * filter_aln_upper_lower ( Alignment *A, Alignment *ST, int use_cons,int value);
Alignment * filter_aln_switchcase ( Alignment *A, Alignment *ST, int use_cons, int value);

Alignment * STseq2STaln ( Alignment *A, Alignment *ST);
Alignment * merge_annotation   ( Alignment *A, Alignment *ST, char *seq);
Alignment * filter_aln_convert ( Alignment *A, Alignment *ST, int use_cons,int value, int n_symbol,char** symbol_list);
int aln2ngap (Alignment *A);

int  * count_in_aln ( Alignment *A, Alignment *ST, int value, int n_symbol,char **symbol_list, int *table);
void count_misc (Alignment*A, Alignment *B);

Alignment * trim_aln_with_seq ( Alignment *S, Alignment *P);
Alignment * add_align_seq2aln ( Alignment *A, char *seq, char *seq_name);
Sequence  * aln2seq    ( Alignment *A);
Sequence  * aln2seq_main    ( Alignment *A, int mode);
Alignment * thread_profile_files2aln (Alignment *A, char *template_file, Fname *F);
Alignment * expand_aln (Alignment *A);
Alignment * aln2expanded_aln (Alignment *A);
Alignment * expand_number_aln (Alignment *A,Alignment *EA);
Alignment * remove_gap_column ( Alignment *A, char *mode);
Alignment*  ungap_sub_aln        ( Alignment *A, int nseq, int *ls);
Sequence *  ungap_seq       ( Sequence *A);
Alignment * insert_gap_col (Alignment *A, int p, int l);
Alignment * unalign_residues (Alignment *A, int i1, int i2);
Alignment * unalign_aln_2 (Alignment *A, Alignment *C, int t);
Alignment * unalign_aln (Alignment *A, Alignment *C, int t);
Alignment * unalign_aln_pos (Alignment *A, int s, int p, int l);

Alignment *degap_aln (Alignment *A);

Alignment * ungap_aln_n        ( Alignment *A, int n);
Alignment * ungap_aln        ( Alignment *A);
void compress_aln     ( Alignment *A);
Alignment* condense_aln (Alignment *A);

Alignment * probabilistic_rm_aa ( Alignment *A, int pos, int len);
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
Alignment *seq_name2removed_seq_name(Sequence *S, Alignment *NA, float **diff);
int seq_name2index (char *name, Sequence *S);

Sequence *extract_one_seq(char *n,int start, int end, Alignment *S,int keep_name);
Sequence  * extract_sub_seq( Sequence  *COOR, Sequence *S);


Sequence * add_prf2seq  (char *alnfile, Sequence *S);
int prf_in_seq ( Sequence *S);
Sequence * add_sequence ( Sequence *IN, Sequence *OUT, int i);
Sequence  * trim_seq     ( Sequence   *A, Sequence  *B);
Sequence  * reorder_seq  ( Sequence   *A, char **name, int nseq);
Sequence  * reorder_seq_2  ( Sequence   *A, int **name,int field, int nseq);

char * concatenate_seq ( Sequence *S, char *conc, int *order);
Sequence * swap_header ( Sequence *S, Sequence *H);

Alignment *aln2jacknife (Alignment *A, int nseq, int len);
char ** name2random_subset (char **in_name, int n_in, int n_out);
Alignment * aln2random_order   ( Alignment  *A);
Alignment * aln2scramble_seq  ( Alignment  *A);

Alignment * reorder_aln        ( Alignment  *A, char **name, int nseq);

char ** rm_name_tag (char **name, int nseq, char *tag);

/******************************************************************************/
/*                          TEMPLATE MANAGEMENENT                             */
/******************************************************************************/
char * string_contains_template_tag (char *string);
Sequence * seq2template_type(Sequence *Seq);

Sequence * vremove_seq_template_files (Sequence *S);
Sequence * display_seq_template_files (Sequence *S);
Sequence * handle_seq_template_file (Sequence *S, char *mode);
int handle_X_template_files ( X_template *T, char *mode);


Sequence * seq2template_seq ( Sequence *S, char *template_file, Fname *F);
char * seq2template_file (Sequence *S, char *file);
int seq2template_file2 (Sequence *S, char *file, char *mode);

Sequence * profile_seq2template_seq ( Sequence *S, char *template_file, Fname *F);
int seq2n_X_template ( Sequence *S, char *type);

struct X_template *fill_X_template (char *name, char *p, char *type);
FILE * display_seq_template (Sequence *S, FILE *io);
char *template_type2type_name (char *type);
char *template_type2short_type_name (char *type);


FILE * display_sequence_templates ( Sequence *S, int i, FILE *io);
FILE * display_X_template (struct X_template *X, FILE *io);

struct X_template* free_X_template ( struct X_template *X);

struct X_template *fill_P_template (char *name, char *p, Sequence *S);
struct X_template *fill_F_template (char *name, char *p, Sequence *S);
struct X_template *fill_S_template ( char *name,char *p, Sequence *S);
struct X_template *fill_R_template (char *name, char *p, Sequence *S);
struct X_template *fill_G_template (char *name, char *p, Sequence *S);
struct X_template *fill_T_template (char *name, char *p, Sequence *S);
struct X_template *fill_E_template (char *name, char *p, Sequence *S);
struct X_template *fill_U_template (char *name, char *p, Sequence *S);

char *seq2T_value ( Sequence *S, int i, char *param_name, char *template_type);
char *profile2P_template_file (Sequence *S, int n);
Alignment * seq2R_template_profile (Sequence *S, int n);
char *seq2P_pdb_id (Sequence *S, int n);
char      * seq2P_template_file (Sequence *S, int n);
char      * seq2T_template_string (Sequence *S, int n);
char      * seq2E_template_string (Sequence *S, int n);
int       * seq2U_template (Sequence *S, int n);

struct X_template * seq_has_template ( Sequence *S, int n, char *type);

/******************************************************************************/
/*                          ALIGNMENT MANIPULATION                            */
/******************************************************************************/

char *aln_column2string (Alignment *A, int p);
int **fix_seq_aln (Sequence *S, Alignment*A, int **cache);
int **fix_seq_seq ( Sequence *S1, Sequence *S2);
int **fix_aln_seq_new (Alignment *S1, Sequence *S2);

Alignment * fix_aln_seq  ( Alignment *A, Sequence *S);
Alignment * rotate_aln ( Alignment *A, char *name);
Alignment * invert_aln ( Alignment *A);
char * complement_string (char *s);
Alignment * complement_aln ( Alignment *A);
Alignment * extract_nol_local_aln( Alignment *A, int start, int max_end);
Alignment * aln2block   (Alignment  *A, int start, int end, Alignment *B);
Alignment * alnpos2block   (Alignment  *A, int*pos, Alignment *B);

Alignment * extract_aln          ( Alignment *A, int start, int end);
Alignment * extract_aln2          ( Alignment *A, int start, int end, char *seq_name);
Alignment * extract_aln3          ( Alignment *A, char *filename);
Alignment * alnpos_list2block (Alignment *A, int n, char **in_list);

Alignment * trunkate_local_aln   ( Alignment *A);
int get_nol_aln_border ( Alignment *A, int start, int direction);
Alignment ** trim_local_aln ( Alignment *A, int **List, int ne, int **residue_list, Sequence *S);

Alignment * aln_cat ( Alignment *A, Alignment *B);
Alignment * concatenate_aln ( Alignment *A, Alignment *B, char *sep);
char * extract_defined_seq ( char *in, int in_of, int in_start, int *aa_def, int dir, int *out_start, char *out_seq);
int verify_aln ( Alignment *A, Sequence *S, char * error);
Alignment * remove_end (Alignment *A);

Alignment *adjust_est_aln ( Alignment *PW, Alignment *M, int s);
Alignment * rename_seq_in_aln (Alignment *A, char ***list);
Sequence * rename_seq_in_seq (Sequence *A, char ***list);
/********************************************************************/
/*                                                                  */
/*                   FLOAT SIMILARITIES                             */
/*                                                                  */
/*                                                                  */
/*                                                                  */
/********************************************************************/
float get_seq_fsim ( char *string1, char *string2, char *ignore, char *similarity_groups, int **matrix, int mode);
float get_seq_fsim2 ( char *string1, char *string2, char *ignore, char *in_mode);
float ** get_fsim_aln_array ( Alignment *A, char *mode);
/********************************************************************/
/*                                                                  */
/*                   ALIGNMENT ANALYSES                             */
/*                                                                  */
/*                                                                  */
/*                                                                  */
/********************************************************************/
int **sim_array2dist_array ( int **p, int max);
int **dist_array2sim_array ( int **p, int max);
int **normalize_array (int **p, int max, int norm);

int aln2most_similar_sequence ( Alignment *A, char *mode);
int aln2coverage ( Alignment *A, int ref_seq);

double aln2entropy (Alignment *A, int *in_ls, int in_ns, float gap_threshold);
int sub_aln2sim ( Alignment *A, int *ns, int **ls, char *mode);
int sub_aln2max_sim ( Alignment *A, int *ns, int **ls, char *mode);
int aln2sim     ( Alignment *A, char *mode);
int seq2idscore_sim ( char *seq1, char *seq2);

int aln_is_aligned ( Alignment *A);
int* get_cdna_seq_winsim ( int *cache, char *string1, char *string2, char *ignore, char *mode, int *w);
int  get_cdna_seq_sim    ( int *cache, char *string1, char *string2, char *ignore, char *mode);

int seq2aln2sim    (char *seq1, char *seq2, char *mode_aln, char *mode_id);
int* get_seq_winsim( char *string1, char *string2, char *ignore, char *mode, int *w);
int  get_seq_sim ( char *string1, char *string2, char *ignore, char *mode);
int  get_seq_sim_2 ( char *string1, char *string2, char *ignore, char **gr, int ng);
int  get_seq_sim_3 ( char *string1, char *string2, char *ignore, int **mat);


int *** get_winsim_aln_array ( Alignment *A, char *mode, int ***w);
int ** get_sim_master_aln_array ( Alignment *A,int n, char *mode);

int ** seq2sim_mat (Sequence *S, char *mode);
int ** seq2cov_mat (Sequence *S, char *mode);
int ** seq2comp_mat (Sequence *S, char *mode, char *comp_mode);

int logid_score (int sim, int len);
int ** fast_aln2sim_mat (Alignment *A, char *mode);
int ** fast_aln2sim_list (Alignment *A, char *mode, int *ns, int **ls);

int ** aln2sim_mat (Alignment *A, char *mode); 
int **aln2cov (Alignment *A);
int ** get_dist_aln_array ( Alignment *A, char *mode);
int ** get_raw_sim_aln_array ( Alignment *A, char *mode);
int ** get_sim_aln_array ( Alignment *A, char *mode);
int generic_get_seq_sim  ( char *seq1, char *seq2, int *cache, char *mode);  
Alignment * grep_seq (Alignment *S,char *field, char *mode, char *string);
Alignment* modify_seq (Alignment *S,char *field, char *string1, char *string2);

Sequence * seq2filter (Sequence *S_in, int min, int max);
int ** get_cov_aln_array ( Alignment *A, char *mode);
int ** get_cov_master_aln_array ( Alignment *A,int n, char *mode);


int * get_aln_col_weight ( Alignment *A, char *mode);
int analyse_aln_column   ( Alignment *B, int col);

int sub_aln2nseq_prf ( Alignment *A, int ns, int *ls);
int **aln2count_mat   (Alignment *A);
int **sub_aln2count_mat2   (Alignment *A, int ns, int *ls);
int **sub_aln2count_mat3   (char **al, int n);
int **aln2count_mat2   (Alignment *A);
char *aln2random_seq (Alignment *A, int noise1, int noise2, int noise3, int gap_noise);

Alignment * master_trimseq( Alignment *A, Sequence *S,char *mode);
Alignment * trimseq( Alignment *A, Sequence *S, char *mode);
Alignment *simple_trimseq (Alignment *A,Alignment*K, char *mode, char *seq);
Alignment *sim_filter (Alignment *A, char *in_mode, char *seq_list);

float ** get_weight ( Alignment *A, Sequence *S, char *mode);
float **seq2pwsim (	   Alignment *A, Sequence *S, char *mode);
Alignment * trimseq( Alignment *A, Sequence *S,char *mode);
Alignment * tc_trimseq( Alignment *A, Sequence *S,char *mode);
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
char** string2alphabet (char *string, int depth, int *falp_size);
int is_in_same_group_aa ( char r1, char r2, int n_group, char **gl, char *mode);
int find_group_aa_distribution (char *col, int nseq,int n_group, char **gl,  int *distrib, char *mode );
char** make_group_aa (int *ngroup, char *mode);
char** make_group_aa_upgma (char *mat, int max_size);


char * test_gene2prot (Constraint_list *CL, int s1);
Alignment* gene2prot (Alignment *A);
Alignment * dna_aln2_3frame_cdna_aln(Alignment *A,int *ns,int **l_s);

int ** get_sim_aln_array_normal_distribution ( Alignment *A, char *mode, int *STD, int *CENTER);
double normal(double x, double mean, double std);
int generic_get_seq_sim_normal_distribution ( char *seq1, char *seq2, int*cache, char *mode, int *STD, int *CENTER);
int get_seq_sim_distribution ( char *string1, char *string2, char *ignore, char *in_mode, int *STD, int *CENTER);

Alignment *aln2clean_pw_aln (Alignment *A,OveralnP *F);
char **pw_aln2clean_pw_aln (char ** aln,OveralnP *F);
int  * pw_aln2clean_aln_weight ( char *seq1, char *seq2, int w, OveralnP *F);

float* aln2pred  ( Alignment *A, Alignment*B, char *mode);
float* analyze_overaln ( Alignment *A, Alignment *B, char *mode, int f,int p1,int p2, int p3,int filter);


Alignment * mark_exon_boundaries  (Alignment *A, Alignment *E);

struct orp
{
  char name[100];
  char mode[100];
  int ncomp;
  int nseq;
  int len;
  
  Alignment *A;
  Alignment *P;
  Alignment *S;
  
  int *pos;
  char ***motif;
  float sp;
  float sn;
  float sen2;
  float best;
  int tp;
  int tn;
  int fp;
  int fn;

  int offset;
  float evalue;
  struct orp *PR; 
};

typedef struct orp ORP;

typedef Alignment * (*filter_func) (Alignment *, Alignment*, int,int, char *);
/************************************************************************************/
/*                ALIGNMENT ANALYZE     : SAR                                            */
/************************************************************************************/
int display_simple_sar_analyze_pair_col (Alignment *A, Alignment *SAR, char *mode);
int **simple_sar_analyze_pair_col ( Alignment *inA, Alignment *SAR, char *mode);
int ***simple_sar_predict ( Alignment *inA, Alignment *SAR, char *mode);
int display_simple_sar_analyze_col ( Alignment *inA, Alignment *SAR, char *mode);
Alignment *sar_analyze4  (Alignment *A, Alignment *SAR, char *name);/*28/08/06*/
Alignment *sar_analyze3  (Alignment *A, Alignment *SAR, char *name);
Alignment *sar_analyze2  (Alignment *A, Alignment *SAR, char *name);
Alignment *sar_analyze  (Alignment *A, Alignment *SAR, char *name);
int aln2sar_column_list ( Alignment *A, char *filter);
float get_sar_sim (char *seq1, char *seq2);
float get_sar_sim2 (char *seq1, char *seq2);
Alignment *aln2weighted_sar_score ( Alignment *A,Alignment *B, char *weight_file, char *compound);
float seq2weighted_sar_score ( char *seq, int **weight);

int sarset2subsarset ( Alignment *A, Alignment *S, Alignment **subA, Alignment **subS, Alignment *SUB);
int sar2subsar (Alignment *A, Alignment *S, Alignment **subA, Alignment **subS, char **slist, int nl);
int sar2subsar_file ( Alignment *A, Alignment *S, char *aln, char *sar);

Alignment *weight2sar (Alignment *A, Alignment *SAR, char *weight_file, int limit);
Alignment * sar2simpred (Alignment *A, Alignment *SAR, char *pos, char *compound, int L,int U );
Alignment * sar2simpred2 (Alignment *A, Alignment *SAR, char *seqlist, char *posfile, char *compound, int L1 );

Alignment *display_sar ( Alignment *A, Alignment *SAR, char *compound);
NT_node sar2tree (Alignment *A, char *mode);
/************************************************************************************/
/*                ALIGNMENT ANALYZE     : SAR FOR OR                                */
/************************************************************************************/

Alignment * or_scan (Alignment *A, Alignment *B, char *param);
Alignment * or_sar  (Alignment *A, Alignment *B, char *param, int print);
ORP * or_loo ( Alignment *inA, Alignment *inS, char *mode, int *pos,int print);


ORP * combine_n_predictions (ORP **R,Alignment *A, Alignment *B);
ORP* combine_2_predictions ( ORP *IN, ORP *TO,Alignment *A, Alignment *B);
ORP * display_or_summary (ORP *CP, char *mode, FILE *fp, int print);

Alignment * or_comp_loo ( Alignment *inA, Alignment *inS, char *mode, int *pos,int print);
int * or_comp_pos ( Alignment *inA, Alignment *inS, char *mode,int print);
float or_id_evaluate ( Alignment *A, Alignment *S, char *mode, int *pos, int print);
char* or_id_evaluate2  ( Alignment *A, Alignment *S, char *mode, int *pos, int print, float *score);
float or_loo_evaluate ( Alignment *A, Alignment *S, char *mode, int *pos, int print);
float or_loo_evaluate2 ( Alignment *A, Alignment *S, char *mode, int *pos, int print);

Alignment * or_test ( Alignment *inA, Alignment *inS, char *mode);
Alignment * or_jack(Alignment *A, Alignment *S, char *param);
Alignment * or_predict(Alignment *A, Alignment *S, char *mode);
Alignment * or_aln2pos_aln (Alignment *A, Alignment *S, char *mode);
ORP* or_self_predict(Alignment *inA, Alignment *inS, char *mode, int *pos, int print);
Alignment * or_sim(Alignment *A, Alignment *S, char *mode);

Alignment *display_pos (Alignment *A, Alignment *B, int *pos, char *mode);

float evaluate_prediction  (Alignment *R, Alignment *P, char *mode, int print);
ORP* new_evaluate_prediction  (ORP *P, char *mode, int print);

Alignment * aln2prediction (Alignment *A,char ***motif, int *pos);
int *   aln2predictive_positions (Alignment *A, Alignment *B, char *mode, int print);
int *   aln2predictive_positions_mat  (Alignment *A, Alignment *B, char *mode, int print);
int *   aln2predictive_positions_scan (Alignment *A, Alignment *B, char *mode, int print);
char *** compounds2motifs (Alignment *A, Alignment *B, int *pos, int depth, char *mode, int print);
char ** compound2motif (Alignment *A, Alignment *B, int *pos, int depth, int c, char *mode, int print);
double pos2sim (Alignment *A, Alignment *B, int *pos);
double  sar_aln2r (Alignment *A, Alignment *B, int *pos, int print);
double sar_aln2delta (Alignment *A, Alignment *B, int *pos, int print);
Alignment * jack_sar(Alignment *A, Alignment *S, char *param);
Alignment *set_sar (Alignment *A, Alignment *S, char *param);
char * get_compound_name (int c, char *mode);
Alignment *get_prediction_target (Alignment *A, Alignment *S, char *param);
int *   file2pos_list (Alignment *A, char *posfile);
ORP * declare_or_prediction ( int ncomp, int nseq, int len);
void free_orp_list ( ORP**P);
void free_orp ( ORP*P);
double evaluate_sar_score1 ( int len, int n11, int n1a, int n1b);
double evaluate_sar_score2 ( int len, int n11, int n1a, int n1b);

Sequence * compare_sar_sequence( Sequence *S1, Sequence *S2, int depth);
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
NT_node ** seq2cw_tree ( Sequence *S, char *file);
NT_node ** make_nj_tree (  Alignment *A,int **distances,int gop, int gep, char **out_seq, char **out_seq_name, int out_nseq, char *tree_file, char *tree_mode);
NT_node ** make_upgma_tree (  Alignment *A,int **distances,int gop, int gep, char **out_seq, char **out_seq_name, int out_nseq, char *tree_file, char *tree_mode);

NT_node ** int_dist2nj_tree (int **distances, char **out_seq_name, int out_nseq,  char *tree_file);
NT_node ** float_dist2nj_tree (float **distances, char **out_seq_name, int out_nseq,  char *tree_file);
NT_node ** dist2nj_tree (double **distances, char **out_seq_name, int out_nseq,  char *tree_file);

NT_node ** int_dist2upgma_tree (int **mat, Alignment *A, int nseq, char *fname);
NT_node upgma_merge (int **mat, NT_node *NL, int *used, int *n, int N);

void nj_tree(char **tree_description, int nseq);
void fast_nj_tree(char **tree_description);
void slow_nj_tree(char **tree_description);

void print_phylip_tree(char **tree_description, FILE *tree, int bootstrap);
void two_way_split(char **tree_description, FILE *tree, int start_row, int flag, int bootstrap);
void guide_tree(char *fname, double **saga_tmat, char **sag_seq_name, int saga_nseq);



NT_node split2upgma_tree (Split **S, Alignment *A, int nseq, char *fname);
NT_node split_upgma_merge (Alignment *A, Split **S, NT_node *NL, int *used, int *n, int N);
float get_split_dist ( Alignment *A, NT_node L, NT_node R, Split **S) ;

Alignment * upgma_tree_aln  (Alignment*A, int nseq, Constraint_list *CL);
int ** dist_mat2best_split (int **mat, int nseq);
int upgma_node_heap (NT_node X);typedef struct Tmpname Tmpname;
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

FILE* print_mem_usage (FILE *fp, char *comment);
void set_max_mem (int m);
int verify_memory (int s);
int my_assert ( void *p, int index);

void * vmalloc ( size_t size);
void * vcalloc ( size_t nobj, size_t size);
void * vcalloc_nomemset ( size_t nobj, size_t size);
void * sub_vcalloc ( size_t nobj, size_t size, int MODE);

void * vrealloc ( void *p, size_t size);
void   vfree2 ( void **p);
void   vfree ( void *p);
void * free_arrayN (void *p, int ndim);
void   vfree_all ();
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



Fname *declare_fname ();
Fname *free_fname ( Fname *F);
/*********************************************************************************************/
/*                                                                                           */
/*         STRUCTURES FOR HSEARCH                                                            */
/*                                                                                           */
/*********************************************************************************************/
#define FIND           0
#define ADD            1
#define REMOVE         2
#define DECLARE        3
#define MARK           4
#define UNMARK         5
#define FREE           6
#define FREE_STACK     7
#define FREE_ALL       8
#define FREE_MARK      9
#define INFO           10
     
struct HaschT
{
  int ne;
  struct Hasch_entry **p;
};
typedef struct HaschT HaschT;

struct Hasch_entry
{
  struct Hasch_entry *n;
  struct Hasch_entry *p;
  int k;
  struct Hasch_data  *data;
  struct Hasch_data * (*free_data)(struct Hasch_data *); 
  struct Hasch_data * (*declare_data)(struct Hasch_entry*);
  int tag;
};
typedef struct Hasch_entry Hasch_entry;
struct Char_node
{
 struct Char_node **c;
 int key;
 
};
typedef struct Char_node Char_node;

HaschT * hcreate ( int n_elements,struct Hasch_data * declare_data(struct Hasch_entry *), struct Hasch_data *free_data(struct Hasch_data *) );
HaschT *hdestroy (HaschT *T,struct Hasch_data * declare_data(struct Hasch_entry *), struct Hasch_data *free_data(struct Hasch_data *) );
Hasch_entry* hsearch (HaschT *T, int k, int action, struct Hasch_data * declare_data(struct Hasch_entry *), struct Hasch_data *free_data(struct Hasch_data *) );
Hasch_entry * extract_hasch_entry_from_list (Hasch_entry *e, struct Hasch_data * declare_data(struct Hasch_entry *), struct Hasch_data *free_data(struct Hasch_data *) );
Hasch_entry * insert_hasch_entry_in_list (Hasch_entry *p, Hasch_entry *e, Hasch_entry *n, struct Hasch_data * declare_data(struct Hasch_entry *), struct Hasch_data *free_data(struct Hasch_data *) );
Hasch_entry * allocate_hasch_entry (Hasch_entry *e, int action,struct Hasch_data * declare_data(struct Hasch_entry *), struct Hasch_data *free_data(struct Hasch_data *) );




 
int string2key (char *s, Char_node *n);
Char_node * declare_char_node (int action);
char * process_repeat (char *aln, char *seq, char *pdb);
char     * normalize_pdb_file  (char *name, char *seq,char *out_file);
Ca_trace * trim_ca_trace (Ca_trace *st, char *seq );

Ca_trace * read_ca_trace (char *file, char *seq_field );
Ca_trace * simple_read_ca_trace (char *file );
Ca_trace * hasch_ca_trace             ( Ca_trace *T);
Ca_trace * hasch_ca_trace_nb          ( Ca_trace *T);
Ca_trace * hasch_ca_trace_bubble      ( Ca_trace *T);
Ca_trace * hasch_ca_trace_transversal ( Ca_trace *TRACE);

float get_atomic_distance ( Atom *A, Atom*B);
float ** measure_ca_distances(Ca_trace *T);

float** print_contacts ( char  *file1, char *file2, float T);
char *  map_contacts ( char  *file1, char *file2, float T);
int * identify_contacts (Ca_trace *ST1,Ca_trace *ST2, float T);
Sequence *seq2contacts ( Sequence *S, float T);
char *string2contacts (char *seq,char *name,char *comment, float T);
char **struc2nb (char *name,char *seq, char *comment, float Threshold, char *atom_list, char *output);
char **struclist2nb (char *name,char *seq, char *comment, float Threshold, char *atom_list, char *output);

typedef struct
{   Alignment *A;
    Alignment *B;
    Alignment *sim_A;
    Sequence  *S;
    Structure *ST;
/*PARAMETERS*/
    char ***grep_list;
    int n_greps;

    char *sim_aln;
    char *alignment1_file;
    char *alignment2_file;
    
    char *io_format;

    int n_structure;
    char **struct_file;
    char **struct_format;
    int *n_symbol;
    char ***symbol_list;

/*LIST VARIABLES*/
    int **code_A;
    int **code_B;
    int n_elementsA;
    int n_elementsB;
    
    int **end_index;
    int **start_index;
/*RESULTS_VARIABLES*/
    int **tot_count;
    int **pos_count;
    int ***pw_tot_count;
    int ***pw_pos_count;
    int *glob;
    int **pw_glob;
/*IO VARIABLES*/
    int n_categories;
    char ***category;
    char *category_list;
    int *n_sub_categories;
    char sep_l;
    char sep_r;
/*Sims VARIABLES*/
    float **sim;
    float **sim_param;
    char *sim_matrix;
    
    int sim_n_categories;
    char ***sim_category;
    char *sim_category_list;
    int *sim_n_sub_categories;
}Result;


#define MAX_N_CATEGORIES 100
#define MAX_N_STRUC      100
    

    

int aln_compare (int argc, char *argv[]);
int **analyse_distance ( Alignment *A, int **dis);

Structure * read_structure (char *fname, char *format, Alignment *A,Alignment *B, Structure *ST, int n_symbols, char **symbol_table);


int is_in_struct_category ( int s1, int s2, int r1, int r2, Structure *ST, char **cat, int n_sub_cat);
char * get_structure_residue (int s, int r, Structure *S);
int parse_category_list ( char *category_list, char ***category, int *sub_n_categories);
int struc_matches_pattern ( char *struc, char *pattern);
float **get_aln_compare_sim ( Alignment *A, Structure *S, char **cat, int n_cat, char *matrix);
float **analyse_sim ( Alignment *A, float **dis);

/*Output*/
FILE *output_format (char *iof, FILE *fp, Result *R);
FILE *output_pair_wise_sequence_results (FILE *fp,  Result *R);
FILE *output_sequence_results (FILE *fp,  Result *R);
FILE *output_total_results (FILE *fp,  Result *R);
FILE *output_header (FILE *fp, Result *R);
FILE *output_large_header ( FILE *fp, Result *R);

/*Parameter Checking*/
int is_a_struc_format (char *format);
void get_separating_char ( char s, char *l, char *r);
void output_informations ();

int check_configuration4program();
typedef struct
    {
      Alignment *A;
      Weights *W;
      Sequence *S;
      int **M;
      Structure *RNA_ST;
      NT_node T;
      Constraint_list *CL;
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
        int keep_name;
        int use_consensus;
}Action_data_struc;

/*Control of alignment sizes*/
int  set_landscape_msa (int len);
int get_msa_line_length (int line, int aln_len);

int seq_reformat (int argc, char **argv);

Sequence_data_struc *read_data_structure ( char *in_format, char *in_file,Action_data_struc *RAD); 
Alignment * main_read_aln ( char *name, Alignment *A);
Sequence  * read_sequences ( char *name);
Sequence  * read_alifold   ( char *name);
Alignment *alifold2aln     ( char *name);
Sequence  * main_read_seq ( char *mname);
int output_format_aln ( char *format, Alignment *A, Alignment *EA,char *name);
int main_output   ( Sequence_data_struc *D1, Sequence_data_struc *D2, Sequence_data_struc *DST, char *out_format, char *out_file);

char * identify_seq_format ( char *file);
char * name2type_name ( char *name);
char identify_format (char **fname);
char **identify_list_format ( char **list, int n);

int type_is_exon_boundaries(char **seq, int n);

int format_is_oligo  ( char *file);
int format_is_msf  ( char *file);
int format_is_fasta( char *file);
int format_is_fasta_aln( char *file);
int format_is_fasta_seq( char *file);
int is_pir_name (char *name);
int format_is_pir  ( char *file);
int format_is_pir_aln( char *file);
int format_is_pir_seq( char *file);
int pir_name (char *name);
int format_is_conc_aln (char *file);
int format_is_saga  ( char *file);
int format_is_swissprot (char *name);

int is_seq ( char *name);
int is_aln ( char *name);
int has_pdb (char *name);
int is_stockhom_aln ( char *name);
int is_blast_file (char *name);
int is_sap_file (char *name);
int is_pdb_file ( char *name);
int is_simple_pdb_file ( char *name);
char *fix_pdb_file (char *name);

int is_pdb_name ( char *name);
char* get_pdb_id(char *name);
char* get_pdb_struc(char *name, int start, int end);
char*  seq_is_pdb_struc ( Sequence *S, int i);
char* is_pdb_struc ( char *name); /*Returns NULL if not a PDB structure Or a the name of a file containing a PDB structure*/
int is_matrix (char *name);

int is_lib (char *name);
int is_lib_01 (char *name);
int is_lib_02 (char *name);
int is_lib_list ( char *name);
int is_single_seq_weight_file (char *fname);
int is_newick (char *name);

int is_method ( char *file);

char *format_name2aln_format_name (char *name);
int is_in_format_list ( char *name);
int is_out_format_list ( char *name);
int is_struc_in_format_list ( char *name);
int is_struc_out_format_list ( char *name);
/*******************************************************************************************/
/*                                                                                         */
/*                                                                                         */
/*                               INPUT MISC                                               */
/*                                                                                         */
/***************************************************************************************** */

char *** read_rename_file ( char *fname, int mode);
void get_barton_list_tc_seq ( char *in_file);
int process_barton_entry (char *buf, char *name);  

Structure *read_rna_struc_number ( Alignment *A, char *fname);
char ** read_lib_list (char *name, int *n);
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
char ***read_group ( char *file);
Sequence* get_pdb_sequence           ( char *fname);
Sequence* get_struc_gor              ( char *fname);
Sequence* get_dialign_sequence       ( char *fname);
Sequence* get_pima_sequence          ( char *fname);
Sequence* get_sequence_dali          ( char *fname);
Sequence* get_pir_sequence           ( char *fname, char *comment_name);
Sequence* perl_reformat2fasta        ( char *perl_script, char *file);

Sequence* get_fasta_tree             ( char *fname, char *comment_name);
Sequence* get_fasta_sequence         ( char *fname, char *comment_name);
Sequence* get_fasta_sequence_num     ( char *fname, char *comment_name);
Sequence* get_fasta_sequence_raw     ( char *fname, char *comment_name);
Sequence *get_file_list ( char *fname);
Sequence *get_tree_file_list ( char *fname);

Sequence* get_gor_sequence           ( char *fname, char *comment_name);
Sequence* get_swissprot_sequence     ( char *fname, char *comment_name);
int  fscanf_seq_name ( FILE *fp, char *sname);

void read_check ( Alignment *A, char *check_file);
void read_stockholm_aln ( char *fname, Alignment *A);
void read_aln ( char *fname, Alignment *A);
void read_number_aln ( char *fname, Alignment *A);
Alignment *read_blast_aln  ( char *fname, Alignment *A);
void read_msf_aln ( char *fname, Alignment *A);
void read_amps_aln ( char *in_file, Alignment *A);
int get_amps_seq_name ( char **name, char* fname);
Alignment *read_gotoh_aln ( char *fname, Alignment *A);

void undump_msa ( Alignment *A, char *tmp);
void dump_msa ( char *file,Alignment *A, int nseq, int *lseq);
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
void  output_similarities (char *file, Alignment *A, char *mode);
void  output_similarities_pw (char *file, Alignment *A, Alignment *B, char *mode);
Alignment * similarities_file2aln ( char *file);
int** input_similarities (char *file, Alignment *A, char *mode);

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
char** clean_seq_names (char **names, int n, int mode);
char *clean_seq_name (char *name, int mode);


void output_pir_seq1 (char *fname, Alignment*A );
void output_pir_seq (char *fname, Alignment*A );
void output_gor_seq (char *fname, Alignment*A );
void output_mult_fasta_seq (char *fname, Alignment*A, int n );

void main_output_fasta_seq ( char *fname, Alignment *A, int header);
void output_fasta_tree ( char *fname, Alignment *A);
void output_fasta_seqS (char *fname, Sequence *S );
void output_fasta_seq1 (char *fname, Alignment*A );
char *output_fasta_seqX (char *name, char *mode, Sequence *S, Alignment *A, int i);

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
char * output_fasta_sub_aln (char *fname, Alignment*A, int ns, int *ls  );
char * output_fasta_sub_aln2 (char *fname, Alignment*A, int *ns, int **ls  );

void ouput_suchard_aln ( char *fname,Alignment*A);
void output_fasta_aln  ( char *fname,Alignment*A);
void output_msf_aln    ( char *fname,Alignment*B);
FILE * output_generic_interleaved_aln (FILE *fp, Alignment *B, int line, char gap, char *mode);
void output_stockholm_aln (char *file, Alignment *A, Alignment *ST);
void output_clustal_aln( char *name, Alignment*B);
void output_strict_clustal_aln( char *name, Alignment*B);
void output_generic_clustal_aln( char *name, Alignment*B, char *format);
void output_saga_aln   ( char *name, Alignment*B);
void output_phylip_aln ( char *name, Alignment*B);
void output_mocca_aln  ( char *name, Alignment*B,Alignment*S);
void output_rnalign    (char *out_file, Alignment*A,Sequence *STRUC);
void output_pw_lib_saga_aln (char *lib_name, Alignment *A );
void output_lib        (char *lib_name, Alignment *A );
void output_compact_aln( char *name, Alignment *B);

void print_sub_aln ( Alignment *B, int *ns, int **ls);
void print_aln ( Alignment *B);
FILE * output_aln( Alignment *B, FILE *fp);


FILE * output_aln_score ( Alignment *B, FILE *fp);
FILE * output_aln_with_res_number ( Alignment *B, FILE *fp);


FILE* output_Alignment ( Alignment *B, FILE *fp);
FILE* output_Alignment_without_header ( Alignment *B, FILE *fp);
FILE * output_Alignment_score ( Alignment *B, FILE *fp);
FILE * output_Alignment_with_res_number ( Alignment *B, FILE *fp);
void output_constraints ( char *fname, char *mode, Alignment *A);

Alignment *input_conc_aln ( char *name, Alignment *A);
void output_conc_aln ( char *name, Alignment *B);
void output_glalign       ( char *name, Alignment *B, Alignment *S);
void output_lalign_header( char *name, Alignment *B);
void output_lalign       ( char *name, Alignment *B);
void output_lalign_aln   ( char *name, Alignment *B);

/**************************************************************************************************/
/*                                                                                                */
/*                                                                                                */
/*                               INPUT/OUTPUT MATRICES                                                  */
/*                                                                                                */
/**************************************************************************************************/
int is_blast_matrix (char *fname);
int is_pavie_matrix (char *fname);
int is_clustalw_matrix (char *fname);

int is_distance_matrix_file (char *name);
int is_similarity_matrix_file (char *name);

void aln2mat (Sequence *S);
void aln2mat_diaa (Sequence *S);
int **seq2latmat ( Sequence *S, char *fname);
int output_mat (int **mat, char *fname, char *alp, int offset);
int ** read_blast_matrix ( char *mat_name);
int output_blast_mat (int **mat, char *fname);
double* mat2cmp (int **mat1, int **mat2);

void output_pavie_mat (int **mat, char *fname, double gep, char *alp);
int ** read_pavie_matrix ( char *mat_name);

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


Sequence * transform_sequence ( Sequence *S, char *mode);
Alignment *translate_splice_dna_aln (Alignment *A,Alignment *ST );
Alignment * mutate_cdna_aln ( Alignment *A);

char *test_dna2gene (char *dna, int *w);
Sequence *dnaseq2geneseq (Sequence *S, int **w);

int ** shift_res_weights ( Sequence *R, int **w, int shift);
int res_weights2min(Sequence *R, int **w);
int res_weights2max(Sequence *R, int **w);
int res_weights2avg(Sequence *R, int **w);
int output_wexons (char *name, Alignment *A);
int scan_res_weights4ac (Sequence *R, int **w, int start, int end, int step);
float *res_weights2accuracy_counts ( Sequence *R, int **w,int T, float *result);
float* genepred_seq2accuracy_counts (Sequence *R, Sequence *T,float *result);
void genepred_seq2accuracy_counts4all (Sequence *R, Sequence *Ts); //JM
float* genepred2accuracy_counts     (char *ref,  char *target , float *result);

char *dna2gene (char *dna, int *w);
char * translate_dna_seq_on3frame (  char *dna_seq, char stop, char *prot);

char * translate_dna_seq ( char *dna_seq, int frame, char stop, char *prot);
int is_stop (char r1, char r2, char r3);
int seq2tblastx_db (char *file,Sequence *S, int strand);

char * back_translate_dna_seq ( char *in_seq,char *out_seq, int mode);     
Alignment *back_translate_dna_aln (Alignment *A);
Sequence  *translate_dna_seqS     (Sequence *S, int frame, int stop);
Alignment *translate_dna_aln (Alignment *A, int frame);
char *dna_seq2pep_seq (char *seq, int frame);

Alignment *clean_gdna_aln (Alignment *A);
Alignment *clean_cdna_aln (Alignment *A);
Alignment *clean_est      (Alignment *A);
/**************************************************************************************************/
/********************************                      ********************************************/
/********************************    PROCESSING        ********************************************/
/*************** ****************                      ********************************************/
void modify_data  (Sequence_data_struc *D1, Sequence_data_struc *D2, Sequence_data_struc *DST, char **action_list,int n_actions, Action_data_struc *RAD);
   
//
// Name MAnipulation
//

Alignment *clean_aln (Alignment *A);
Sequence *clean_sequence ( Sequence *S);
char ** translate_names (int n, char **name);
char * translate_name ( char *name);
char *decode_name (char *name, int mode);
FILE * display_sequences_names (Sequence *S, FILE *fp, int check_pdb_status, int print_templates);
Sequence *add_file2file_list (char *name, Sequence *S);
