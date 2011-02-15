/*DEBUGGING*/
/*#include "mshell.h"*/
/*MEMORY MANAGEMENT*/
#include <float.h>
#define MY_EPS 1000*DBL_EPSILON
//Maximum number of tries for interactibve things
#define MAX_N_TRIES 3

//Maximum CACHE and Temporary file size and age (Mb and days, 0: unlimited)
#define TMP_MAX_SIZE  0
#define TMP_MAX_KEEP  10
#define CACHE_MAX_SIZE  2000
#define CACHE_MAX_KEEP  180
#define MAX_N_PID       260000
//Importnat Values Affecting the Program Behavior
#define O2A_BYTE         50
#define SCORE_K          10
#define NORM_F           1000
#define PAVIE_MAT_FACTOR 1000
#define MAXID            100
#define CLEAN_FUNCTION   NULL
#define MINSIM_4_TCOFFEE 25 //The minimum similarity between a sequence and its PDB template
#define MINCOV_4_TCOFFEE 25 //The minimum similarity between a sequence and its PDB template


#define TRACE_TYPE       int
#define MAX_LEN_FOR_DP   600


#define GIVE_MEMORY_BACK 0
#define MEMSET0   1
#define NO_MEMSET0 0
/*OUTPUT DEFINITIONS*/
#define  NO_COLOR_RESIDUE 127
#define  NO_COLOR_GAP 126
#define  CLOSE_HTML_SPAN -1
/*SPECIAL_CODES*/
#define GAP_CODE 60
/*TYPE DEFINITIONS*/

//Formats
#define BLAST_XML 100
#define BLAST_TXT 101

/*SWITCHES*/


#define USED 1
#define UNUSED 2


#define TEMPLATES 1
#define NOTEMPLATES 0

#define EXTEND 1
#define RESIZE 2

#define SEN                0 
#define SPE                1 
#define REC                2 
#define SEN2              2 

#define ALL               1
#define SEGMENTS          2
#define DIAGONALS         3

#define START_STATE       0
#define END_STATE         1

#define KEEP_CASE         2 /*Hard set in several places*/
#define LOWER_CASE        0
#define UPPER_CASE        1
#define CHANGE_CASE       3
#define KEEP_GAP          0
#define RM_GAP            1

#define KEEP_NAME         1

#define CHECK             0
#define NO_CHECK          1
#define FORCE             2
#define STORE             3
#define FLUSH             4
#define DUMP              5


#define ON                8
#define OFF               9
#define LOCKED_ON         10
#define LOCKED_OFF        11

#define YES               12
#define NO                13
#define MAYBE             14

#define NEVER             15
#define ALWAYS            16
#define SOMETIMES         17

#define UPPER             18
#define LOWER             19
#define DELETE            20
#define SWITCHCASE        21 

#define VECTOR            22
#define NON_VECTOR        23
#define NON_PROFILE       24
#define BOOTSTRAP         25

#define HEADER            26
#define NO_HEADER         27

#define VERY_VERBOSE      28
#define VERBOSE           29
#define SHORT             30
#define VERY_SHORT        31

#define OVERLAP           32
#define NO_OVERLAP        33

#define PRINT             34
#define NO_PRINT          35

#define FREE_ALN              36
#define DECLARE_ALN           37
#define EXTRACT_ALN           38
#define CLEAN                 39
#define INTERACTIVE           40
#define NON_INTERACTIVE       41
#define PAD                   42
#define NO_PAD                43

#define SET               44
#define UNSET             45
#define RESET             48
#define ISSET             49
#define GET               50

#define ENV               52
#define LLOCK             53
#define LERROR            54
#define LWARNING          55
#define LSET              56
#define LRESET            57
#define LCHECK            58
#define LREAD             59
#define LRELEASE          60

#define RETURN_ON_FAILURE    61
#define EXIT_ON_FAILURE    62
#define IGNORE_FAILURE     63

#define _START 64
#define _TERM  65

#define GOP               0
#define GCP               1
#define GEP               2

#define BOTTOM             0
#define TOP                1

#define FORWARD            -1
#define BACKWARD            1

#define GO_LEFT            -1
#define GO_RIGHT            1

#define LOCAL            1
#define GLOBAL           2
#define LALIGN           3
#define MOCCA            4

#define TRUE             1
#define FALSE            0

#define NEW              1
#define OLD              0

#define RANDOM           0
#define DETERMINISTIC    1

#define GREEDY           1
#define NON_GREEDY       0

#define IS_FATAL         1
#define IS_NOT_FATAL     0
#define NO_REPORT        2
#define INSTALL          3
#define INSTALL_OR_DIE   4

#define OPTIONAL         1
#define NON_OPTIONAL     0

#define GV_MAXIMISE      1
#define GV_MINIMISE      0

#define MAXIMISE      1
#define MINIMISE      0

#define ALLOWED          0
#define FORBIDEN         -99999999
#define END_ARRAY        -99999990
#define SOFT_COPY 1
#define HARD_COPY 2

#define VERY_SLOW 0
#define SLOW 1
#define FAST 2
#define VERY_FAST 3
#define SUPER_FAST 4
#define ULTRA_FAST 5

#define CODE 1
#define DECODE 2
#define CODELIST 3

/*Identity measure*/
#define UNGAPED_POSITIONS 1
#define ALIGNED_POSITIONS 2
#define AVERAGE_POSITIONS 3
#define NOMATRIX         NULL
#define NOGROUP          NULL
#define NOALN            NULL

/*SIZE DEFINITIONS*/
#define SIZE_OF_INT      10
#define UNDEFINED        FORBIDEN
#define UNDEFINED_INT    UNDEFINED
#define UNDEFINED_FLOAT  UNDEFINED
#define UNDEFINED_DOUBLE UNDEFINED
#define UNDEFINED_CHAR   125
#define UNDEFINED_SHORT  -125
#define UNDEFINED_2      0
#define UNDEFINED_RESIDUE '>'



#define FACTOR           1
#define MAX_N_SEQ        1
#define MAX_N_ALN        1
#define MAX_LEN_ALN      1
#define MAX_N_LIST       100

#define COMMENT_SIZE     1000
#define MAXNAMES         100
#define FILENAMELEN 	 500            /* Max. file name length */
#define MAX_N_PARAM      2000
#define MAX_PARAM_LEN    200
#define MAX_LINE_LENGTH  10000
#define ALN_LINE_LENGTH  50
#define SHORT_STRING     10
#define STRING           300
#define LONG_STRING      1000
#define VERY_LONG_STRING 10000

#define AA_ALPHABET            "acdefghiklmnpqrstvwy-ACDEFGHIKLMNPQRSTVWY"
#define DNA_ALPHABET           "AGCTUNRYMKSWHBVD-agctunrymkswhbvd"
#define RNAONLY_ALPHABET       "Uu"
#define BLAST_AA_ALPHABET      "arndcqeghilkmfpstwyvbzx*"
#define NAMES_ALPHABET         "1234567890ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_|�-!%@&#-+=."

#define SIZEOF_AA_MAT   60
#define GAP_LIST         "-.#*~"
#define SSPACE           "   "

#define MATCH            1
#define UNALIGNED        2
#define GAP              3

#define MNE 3
#define CODE4PROTEINS  10
#define CODE4DNA       20

#define STOCKHOLM_CHAR 'z'
#define STOCKHOLM_STRING "z"


/*CODE SHORT CUTS*/

/*1-COMMAND LINE PROCESSING*/
#define GET_COMMAND_LINE_INFO ((strncmp ( argv[1], "-h",2)==0)||(strncmp ( argv[1], "-man",4)==0)||(strncmp ( argv[1], "-",1)!=0))
#define NEXT_ARG_IS_FLAG ((argc<=(a+1)) ||(( argv[a+1][0]=='-') && !(is_number(argv[a+1]))))


/*UTIL MACROS*/
#define BORDER(p1,l1,p2,l2) ((p1==0 || p2==0 || p1==l1 || p2==l2)?1:0) 
#define GET_CASE(f,c) ((f==UPPER_CASE)?toupper(c):((f==LOWER_CASE)?tolower(c):c))

#define SWAP(x,y) {x=x+y;y=x+y; x=y-x; y=y-2*x;}
#define SWAPP(x,y,tp) {tp=y;y=x;x=tp;}

#define MAX(x, y) (((x) >(y)) ? (x):(y))
#define MAX2(x, y) (((x) >(y)) ? (x):(y))
#define MAX3(x,y,z) (MAX(MAX(x,y),z))
#define MAX4(a,b,c,d) (MAX(MAX(a,b),MAX(c,d)))
#define MAX5(a,b,c,d,e) (MAX2((MAX3(a,b,c)),(MAX2(d,e))))
#define MAX6(a,b,c,d,e,f) (MAX2((MAX3(a,b,c)),(MAX3(c,d,e))))

#define MIN(x, y) (((x) <(y)) ? (x):(y))
#define FABS(x) ((x<0)?(-x):(x))
#define is_defined(x) ((x==UNDEFINED)?0:1)
#define a_better_than_b(x,y,m) ((m==1)?(((x)>(y))?1:0):(((x)<(y))?1:0))
#define is_in_range(x,min,max) ((x>=min && x<=max)?1:0)
/*#define bod_a_b(x,y,m)   ((m==1)?(MAX((x),(y))):(MIN((x),(y))))
#define bo_a_b(x,y,m)    ((x==UNEFINED)?y:((y==UNDEFINED)?x:bod_a_b(y,y,m)))
#define best_of_a_b(x,y,m)   ((x==UNDEFINED && y==UNDEFINED)?(UNDEFINED):(bo_a_b(x,y,m)))
*/


#define DIE(x)  HERE(x);exit(0);
#define best_of_a_b(x,y,m) ((m==1)?(MAX((x),(y))):(MIN((x),(y))))

#define strm(x,y)            ((vstrcmp((x),(y))==0)?1:0)
#define strnm(x,y,n)           ((vstrncmp((x),(y),(n))==0)?1:0)
#define strm2(a,b,c)         (strm(a,b) || strm(a,c))
#define strm3(a,b,c,d)       (strm2(a,b,c) || strm(a,d))
#define strm4(a,b,c,d,e)     (strm2(a,b,c) || strm2(a,d,e))
#define strm5(a,b,c,d,e,f)   (strm2(a,b,c) || strm3(a,d,e,f))
#define strm6(a,b,c,d,e,f,g) (strm3(a,b,c,d) || strm3(a,e,f,g))
#define declare_name(x) (x=vcalloc (MAX(FILENAMELEN,L_tmpnam)+1, sizeof (char))) 
#define is_parameter(x) (x[0]=='-' && !isdigit(x[1])) 

/*Freing functions*/
#define free_2(a, b)            free(a);free(b)
#define free_1(a)               free(a)
#define free_3(a, b, c)         free_2(a,b);free_1(c)
#define free_4(a, b, c,d)       free_2(a,b);free_2(c,d)
#define free_5(a, b, c,d,e)     free_3(a,b,e);free_2(c,d)
#define free_6(a, b, c,d,e,f)   free_3(a,b,e);free_3(c,d,f)
#define free_7(a, b, c,d,e,f,g) free_3(a,b,e);free_4(c,d,f,g)
/*2-FILE PARSING*/
#define SEPARATORS "\n \t,;"
#define LINE_SEPARATOR "\n#TC_LINE_SEPARATOR\n"
#define TC_REC_SEPARATOR "#### TC REC SEPARATOR ###"

/*END 1-*/


/*WIDOWS/UNIX DISTINCTIONS
#if defined(_WIN32) || defined(__WIN32__) ||  defined(__WINDOWS__) || defined(__MSDOS__) || defined(__DOS__) || defined(__NT__) || defined(__WIN32__)
#define WIN32
#define TO_NULL_DEVICE " >nul"
#define    NULL_DEVICE "nul"
#define CWF "/" 
#else
#define TO_NULL_DEVICE " >/dev/null 2>&1"
#define    NULL_DEVICE "/dev/null"
*/

#if defined(_WIN32) || defined(__WIN32__) ||  defined(__WINDOWS__) || defined(__MSDOS__) || defined(__DOS__) || defined(__NT__) || defined(__WIN32__)
#define WIN32
#define TO_NULL_DEVICE " >>t_coffee.log"
#define    NULL_DEVICE "t_coffee.log"
#define CWF "/" /*ClustalW Flag*/
#else
#define TO_NULL_DEVICE " >>/dev/null 2>&1"
#define    NULL_DEVICE "/dev/null"


#define CWF "-" /*ClustaW Flag*/
#endif

/*Generic Data*/
#define EMAIL "cedric.notredame@europe.com"
#define URL "http://www.tcoffee.org"

#define PERL_HEADER "#!/usr/bin/env perl"

//Optimize the Score Computation in DP
#define TC_SCORE_2(x,y) (SCORE_K*CL->M[Aln->seq_al[l_s[0][0]][x]-'A'][Aln->seq_al[l_s[1][0]][y]-'A']-SCORE_K*CL->nomatch) 
#define TC_SCORE_N(x,y) ((CL->get_dp_cost)(Aln, pos, ns[0], l_s[0], x, pos, ns[1], l_s[1], y, CL))
#define TC_SCORE(x,y)  ((CL->get_dp_cost==slow_get_dp_cost && CL->evaluate_residue_pair==evaluate_matrix_score && ns[0]+ns[1]==2 && x>=0 && j>=0)? (TC_SCORE_2(x,y)):(TC_SCORE_N(x,y)))

#define NULL_2 NULL,NULL
#define NULL_3 NULL_2,NULL
#define NULL_4 NULL_2,NULL_2
#define NULL_5 NULL_3,NULL_2
#define NULL_6 NULL_4,NULL_2
#define NULL_7 NULL_5,NULL_2
typedef struct
    {
    char *mode;
    char *comments;
    int nseq;
    char **seq_name;
    float **PW_SD;
    float **PW_ID;
    float *SEQ_W;
    }Weights;

typedef struct
    {
    int **list;
    int tot_list;
    int **stem;
    int tot_stem;
    int n_fields;
    int nseq;
    int *len;
    int ***struc;
    struct Sequence *S;  
    }Structure;

struct Sequence
    {
      char **file;          /* file[Nseq][FILENAMELEN] name of the file that contributed each sequence*/
      char **seq_comment;     /* seq_comment[Nseq][LONG_STRING] comment read in the file */
      char **aln_comment;  /*id*/
      char **seq;          /*seq[Nseq][sequence] sequences*/
      int *len;            /*len[Nseq] length of each sequence*/
      int max_len;         /*Lenght of the longest seq */
      int min_len;         /*Length of the shortest seq*/
      int nseq;            /*nseq*/
      int max_nseq;        /*Maximum number of sequences in the datastruct*/
      char **name;         /*name[Nseq][MAXNAMELEN]*/
      int **dc;         /*coordinates on the disk. Coordinates set if seq[i]==NULL*/
/*Constraint list*/
      struct Constraint_list *CL;
      int contains_gap;   /*set to 1 if gaps are to be kept*/
      char *type;         /*PROTEIN, DNA*/
      Weights *W;         /*Associated weights*/
      char template_file[FILENAMELEN+1];
      struct Template **T;
      
};
typedef struct Sequence Sequence;

//_E_
struct Template
{
  char seq_type[10];
  struct X_template *P;//PDB structure
  struct X_template *F;//RNA secondary structure
  struct X_template *S;//sequence
  struct X_template *R;//Profile
  struct X_template *G;//Genomic structure
  struct X_template *T;//transmembrane
  struct X_template *E;//secondary structure
  struct X_template *U;//Unicode, strings
  
  struct X_template *RB;
};
typedef struct Template Template;
//_E_
struct X_template 
{
  char seq_name[FILENAMELEN+1];
  char template_type[FILENAMELEN+1];
  char template_format[100];
  char template_name[FILENAMELEN+1];
  char template_file[FILENAMELEN+1];
  
  struct P_template *VP; 
  struct F_template *VF;
  struct S_template *VS;
  struct R_template *VR;
  struct G_template *VG;
  struct T_template *VT;
  struct E_template *VE;
  struct U_template *VU;
  
  
};
typedef struct X_template X_template;

//
struct P_template
{
  char pdb_id[100];
};
typedef struct P_template P_template;

//RNA secondary Structure
struct F_template
{
  int l;
};
typedef struct F_template F_template;


struct S_template
{
  Sequence *S;
};
typedef struct S_template S_template;

//Prile associated with a sequence
struct R_template
{
  struct Alignment *A;
};
typedef struct R_template R_template;

//Genomic Information
struct G_template
{
  Sequence *S;
};
typedef struct G_template G_template;


struct T_template
{
  Sequence *S;
};
typedef struct T_template T_template;

//_E_
struct E_template
{
  Sequence *S;
};
typedef struct E_template E_template;

struct U_template
{
  int *list;
};
typedef struct U_template U_template;


typedef struct
    {
    int max_len;
    int alp_size;
    char *alphabet;
    int **count3;  
    int **count;
    int **count2;  
    }Profile;

struct Alignment
    {
/*Size*/
    int max_len;
    int min_len;   
    int *  len;
      //int *weight;  
    int declared_len;
    int max_n_seq;
    int nseq;
    int len_aln;
/*Generic Information*/
      char *generic_comment;
/*Sequence Information*/
    char **file;
    char **seq_comment;
    char **aln_comment;
    char **name;
      
    char **expanded_order;
    char **tree_order;
    char **seq_al;
    
    int  **order;
    Profile *P;
    Sequence *S;
    struct Dp_Result *Dp_result;
    struct Constraint_list *CL;

    int **seq_cache; /*Contains the index of the residues:
		       The sequence Numbering is relative to the sequences, and not to the alignmnent
		       
		       seq_cache[0][1]=3
		       indicates that in the aln residue (0)1 corresponds to [order[0][0]][3]
		       residues: 1...N
		       Sequences 0...M
		     */
    int **cdna_cache; /*Contains the information about wheather a nucleotide is coding or not*/
                     /*Only defined if used */


/*Weight*/   
      float *col_weight;
      float *seq_weight;
      float **res_weight;
/*Score*/
    int *  score_seq;
    int ** score_res;
    int score_aln;
    int score;
    int ibit;  
    int cpu;
    int finished;

/*Input/Output Options*/
    int output_res_num;
    int residue_case; /*1 for lower, 0 for Upper, 2 for keeping unchanged*/
      int expand;
    int output_tm;
/*Must Not be copied*/
     int used;
     int num;     
     int **pos; 
/*For linked lists*/      
    struct Alignment * A;  
      /*Misc*/
    int random_tag;
     
    };
    
typedef struct Alignment Alignment;
typedef struct
    {
    int in_seq;	
    FILE *fp;
    int font;
    int x0;
    int y0;
    int x;
    int y;
    int n_pages;
    int max_line_ppage;
    int n_line;
    int line;
    int eop;
    int in_html_span;
    char previous_html_color[100];
   
    }
FILE_format;

typedef struct
    {
    float r;
    float g;
    float b;
    char html_color[30];
    char html_color_class[30]; 
    int ascii_value;
    }
Color;


Sequence * fill_sequence_struc ( int nseq, char **sequences, char **seq_name);
Sequence * cw_read_sequences ( char *seq_name);
Sequence * get_sequence_type (Sequence *S);
char     * get_array_type (int n, char **s);
Alignment* get_aln_type (Alignment *A);

char     * get_string_type   (char *string);

char *store_mode (char *val);
char *retrieve_mode ();
char *unset_mode ();
char *set_mode (int mode, char *val);

char *store_seq_type (char *val);
char *retrieve_seq_type ();
char *unset_seq_type ();
char *set_seq_type (int mode, char *val);

void get_sequence (char *seq_file,int *NSEQ, char ***SEQ, char ***SN, int **sl, int *min, int *max);

int ** get_matrix   ( char *name, char *format);
int ** read_matrice (char *mat_name);
int **neg_matrix2pos_matrix ( int **matrix);   


void   print_aln ( Alignment *B);

int       output_reliability_ps     ( Alignment *B,Alignment *S, char *name);
int       output_reliability_pdf    ( Alignment *B,Alignment *S, char *name);
int       output_reliability_html   ( Alignment *B,Alignment *S, char *name);
int       output_color_ps     ( Alignment *B,Alignment *S, char *name);
int       output_color_pdf    ( Alignment *B,Alignment *S, char *name);
int       output_color_html   ( Alignment *B,Alignment *S, char *name);
int       output_hit_color_html   (Alignment *B, float **ffPScoreTable, int nl, char *name);	//JM_ADD
void 	  output_hit_matrix(char *fileName, float **ffpHitScoreMatrix, int nl);		//JM_ADD
void      get_rgb_values(int val, Color *C);
int       output_reliability_format ( Alignment *B,Alignment *S, char *name, \
FILE_format *(*vfopen_format)          ( char *),\
FILE_format *(*print_format_string)    ( char * ,Color *, Color *, FILE_format*),\
FILE_format *(*print_format_char)      ( int    ,Color *, Color *, FILE_format*),\
void         (*get_rgb_values_format)  ( int    ,Color *),\
FILE_format* (*vfclose_format)         ( FILE_format *));
int       output_score_format ( Alignment *B,Alignment *S, char *name, \
FILE_format *(*vfopen_format)          ( char *),\
FILE_format *(*print_format_string)    ( char * ,Color *, Color *, FILE_format*),\
FILE_format *(*print_format_char)      ( int    ,Color *, Color *, FILE_format*),\
void         (*get_rgb_values_format)  ( int    ,Color *),\
FILE_format* (*vfclose_format)         ( FILE_format *));


FILE_format * print_ps_string      ( char *s , Color *box, Color *ink, FILE_format *f);
FILE_format * print_ps_char        ( int   c,    Color *box, Color *ink, FILE_format *f);



void get_rgb_values_ps ( int val, Color *C);
FILE_format* vfopen_ps ( char *name);
FILE_format* vfclose_ps ( FILE_format *fps);

FILE_format *print_html_string( char *s, Color *box, Color *ink, FILE_format *fhtml);
FILE_format * print_html_char ( int c, Color *box, Color *ink, FILE_format *f);
void get_rgb_values_html ( int val, Color *C);
FILE_format* vfopen_html ( char *name);
FILE_format* vfclose_html ( FILE_format *fhtml);

int       output_reliability_ascii     ( Alignment *B,Alignment *S, char *name);
int       output_color_ascii           ( Alignment *B,Alignment *S, char *name);

FILE_format *print_ascii_string( char *s, Color *box, Color *ink, FILE_format *fascii);
FILE_format * print_ascii_char ( int c, Color *box, Color *ink, FILE_format *f);
void get_rgb_values_ascii ( int val, Color *C);

FILE_format* vfopen_ascii ( char *name);
FILE_format* vfclose_ascii ( FILE_format *fascii);
int       output_seq_reliability_ascii     ( Alignment *B,Alignment *S, char *name);
/*********************CLUSTALW.H*********************************************/
/****************************************************************************/

   /*
   Main header file for ClustalW.  Uncomment ONE of the following 4 lines
   depending on which compiler you wish to use.
   */

#define VMS 1                 /*VAX or ALPHA VMS */

/*#define MAC 1                 Think_C for MacIntosh */

/*#define MSDOS 1               Turbo C for PC's */

/*#define UNIX 1                Ultrix/Decstation, Gnu C for 
                                Sun, IRIX/SGI, OSF1/ALPHA */

/***************************************************************************/
/***************************************************************************/





#define MAXTITLES		60      /* Title length */

	
#define UNKNOWN   0
#define EMBLSWISS 1
#define PIR 	  2
#define PEARSON   3
#define GDE    	  4
#define CLUSTAL   5	/* DES */
#define MSF       6 /* DES */
#define USER      7	/* DES */

#define PAGE_LEN       22   /* Number of lines of help sent to screen */

#ifdef VMS						/* Defaults for VAX VMS */
#define DIRDELIM ']'		/* Last character before file name in full file 
							   specs */
#define SEQ_MAX_LEN		10000	/* Max Sequence Length */
#define MAXN		500		/* Max Number of Sequences */
#define FSIZE       25000   /* Work space for pairwise alignments */
#define MAXTREE		5000	/* Max Nodes for phylogenetic tree */
#define LINELENGTH    	60  /* Output line length */
#define GCG_LINELENGTH 	50  /* Output line length for GCG output */

#elif MAC
#define DIRDELIM ':'
#define SEQ_MAX_LEN		1000
#define MAXN		30
#define FSIZE       5000
#define MAXTREE		1000
#define LINELENGTH     	50
#define GCG_LINELENGTH 	50


#elif MSDOS
#define DIRDELIM '\\'
#define SEQ_MAX_LEN		1300
#define MAXN		30
#define FSIZE           5000
#define MAXTREE		1000
#define LINELENGTH     	50
#define GCG_LINELENGTH 	50

#elif UNIX
#define DIRDELIM '/'
#define SEQ_MAX_LEN		10000
#define MAXN		500
#define FSIZE           25000
#define MAXTREE		5000
#define LINELENGTH     	60
#define GCG_LINELENGTH 	50
#endif

#define NUMRES 26		/* max size of comparison matrix */

#define INPUT 0
#define ALIGNED 1

#define LEFT 1
#define RIGHT 2

#define NODE 0
#define LEAF 1

#define GAPCOL 32		/* position of gap open penalty in profile */
#define LENCOL 33		/* position of gap extension penalty in profile */

typedef struct node {		/* phylogenetic tree structure */
        struct node *left;
        struct node *right;
        struct node *parent;
        float dist;
        int  leaf;
        int order;
        char name[64];
} stree, *treeptr;

void		*ckalloc(size_t);
void 		* ckvrealloc(void *,size_t);
void 		ckfree(void *);

int readseqs(char *saga_file,char ***SAGA_SEQ, char*** SAGA_NAMES, int ***SAGA_LEN) ;/*first_seq is the #no. of the first seq. to read */


typedef struct treesim{
  float w;
  float uw;
  float d;
  
  float max_w;
  float max_uw;
  float max_d;
  
  int rf;
  int n;//n nodes;
  int nseq;// nseq in the common subset
    }Tree_sim;


typedef struct tnode *NT_node;

/**
* Node of a tree
*/
typedef struct tnode{
  int visited;
  char *name;
  char *file;
  
  ///The parent node
  NT_node parent;
  ///Left child node
  NT_node left;
  ///Right child node
  NT_node right;
  NT_node bot;
  /// is leaf?
  int isseq;
  int seq;
  int maxnseq;
  int nseq;
  
  ///contains a list of the sequences
  int *lseq; 
  ///contains a coded version of the node: 10010101
  int *lseq2;
  ///contains distances to the root, in nodes
  int *idist;
  ///contains real distances *1000
  int *ldist;
  float dist;
  float bootstrap;
  float dp;
  int order;
  int aligned;
  ///Number of leave below the considered node
  int leaf;
  ///Number of nodes below the considered node
  int node;
  int group;
  float score;
  int align;
  char *seqal;
  int index;
  int fork;
    }Treenode;

typedef struct split_struc Split;

typedef struct split_struc{
  char *split;
  int n;
  int tot;
  float score;
  char **tlist;//Not used yet
  Sequence *S;
  NT_node *L;
}Split_struc;

NT_node main_prune_tree ( NT_node T, Sequence *S);
NT_node prune_tree ( NT_node T, Sequence *S);
/*********************************************************************/
/*                                                                   */
/*                                   dpa_tree_manipulation           */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
char *tree2Ngroup (Alignment *A, NT_node T, int max_n, char *fname, char *mat4dist);
int tree2group_file ( NT_node T,Sequence *S, int maxnseq, int minsim, char *name);

NT_node seq2dpa_tree  (Sequence *S, char *align_mode);
NT_node tree2dpa_tree (NT_node T, Alignment *A, char *matrix4distance);
FILE * tree2group ( NT_node T,Sequence *S,int maxnseq, int mindist,char *name, FILE *fp);


NT_node  tree2collapsed_tree (NT_node T, int n, char **string);

/*********************************************************************/
/*                                                                   */
/*                                   tree comparison                 */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
int main_compare_cog_tree (NT_node T1, char *cogfile);
int main_compare_aln_tree (NT_node T1, Alignment *A, FILE *fp);
int compare_aln_tree (NT_node T, Alignment *A, int *n, FILE *fp);

int main_compare_splits (NT_node T1, NT_node T2, char *mode, FILE *fp);
Tree_sim * tree_cmp( NT_node T1, NT_node T2);
NT_node tree_scan (Alignment *A,NT_node RT, char *pscan, char *ptree);

int print_node_list (NT_node T, Sequence *S);
NT_node main_compare_trees ( NT_node T1, NT_node T2, FILE *fp);
NT_node main_compare_trees_list ( NT_node T1, Sequence *S, FILE *fp);

float compare_trees ( NT_node T1, NT_node T2, int nseq, int mode);
float search_node ( NT_node B, NT_node T, int nseq, int mode);
float evaluate_node_similarity ( NT_node B, NT_node T, int nseq, int mode);

int compare_node ( int *b1, int *b2, int n);
void display_node (NT_node N, char *string,int nseq);
NT_node index_tree_node    (NT_node T);
NT_node simple_recode_tree (NT_node T, int nseq);
NT_node recode_tree ( NT_node T, Sequence *S);
int compare_branch2 ( int *b1, int *b2, int n);

/*********************************************************************/
/*                                                                   */
/*                                   FJ_tree Computation             */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
NT_node similarities_file2tree (char *mat);
NT_node tree_compute ( Alignment *A, int n, char ** arg_list);
static NT_node compute_std_tree (Alignment *A, int n, char **arg_list);
NT_node compute_std_tree_2 (Alignment *A, int **s, char *arg_list);
NT_node aln2fj_tree(NT_node T, Alignment *A, int limit,char* mode);
Alignment * filter_aln4tree (Alignment *A, int n,int fg,char* mode);

/*********************************************************************/
/*                                                                   */
/*                                   Tree Filters and MAnipulation   */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
int  tree2star_nodes (NT_node R, int n_max);
NT_node aln2tree (Alignment *A);
NT_node reset_boot_tree ( NT_node R, int n);
NT_node tree_dist2normalized_tree_dist ( NT_node R, float max);
NT_node reset_dist_tree ( NT_node R, float n);
NT_node* free_treelist  ( NT_node *R);
NT_node free_tree  ( NT_node R);
NT_node realloc_tree( NT_node R, int n);
NT_node free_tree_node ( NT_node R);

Sequence * tree2seq    (NT_node R, Sequence *S);
NT_node  rename_seq_in_tree ( NT_node R, char ***list);

NT_node balance_tree (NT_node);
int tree2nseq ( NT_node R);
int tree_file2nseq ( char *file);

int tree2nleaf ( NT_node R);
int tree2nnode ( NT_node R);
int tree2_nnode_unresolved (NT_node R, int *l);

FILE* display_tree ( NT_node R, int n, FILE *fp);
void clear_tree (NT_node T);
int display_leaf ( NT_node T, FILE *fp);
int display_leaf_below_node ( NT_node T, FILE *fp);
NT_node display_leaf_nb (NT_node T, int n, FILE *fp, char *name);
NT_node display_splits (NT_node T,Sequence *S, FILE *fp);
int tree2split_list (NT_node T, int nseq, int **split_list, int *n);

NT_node reroot_tree ( NT_node TREE, NT_node T);
NT_node straighten_tree ( NT_node P, NT_node C, float new_dist);
NT_node unroot_tree ( NT_node T);
FILE* print_tree_list ( NT_node *T,char *format, FILE *fp);
FILE* print_tree ( NT_node T,char *format, FILE *fp);
char *tree2string (NT_node T);
char *tree2file   (NT_node T, char *name, char *mode);

int print_newick_tree ( NT_node T, char *name);
FILE * rec_print_tree ( NT_node T, FILE *fp);


NT_node find_longest_branch ( NT_node T, NT_node L);
NT_node shift_root ( NT_node R);

int ** tree2cluster (NT_node T, float thres);
int ** make_sub_tree_list ( NT_node **T, int nseq, int n_node);
void make_all_sub_tree_list ( NT_node N, int **list, int *n);
void make_one_sub_tree_list ( NT_node T, int *list);
NT_node main_read_tree(char *treefile);

NT_node new_read_tree ( char *teefile);
NT_node new_get_node (NT_node T, FILE *fp);


NT_node** simple_read_tree(char *treefile);
void free_read_tree (NT_node **BT);
NT_node** read_tree(char *treefile, int *nnodes,int nseq, char **seq_names);
FILE * create_linear_tree ( char **name, int n, FILE *fp);
FILE * create_tree(NT_node ptree, NT_node parent,int *numseq,int  *ntotal,int  *nnodes,NT_node **lu, FILE *fp);
NT_node declare_tree_node (int nseq);
void set_info(NT_node p, NT_node parent, int pleaf, char *pname, float pdist, float bootstrap);
NT_node insert_tree_node(NT_node pptr);
FILE * skip_space(FILE *fd);
void create_tree_node(NT_node pptr, NT_node parent);
float calc_mean(NT_node nptr, float *maxdist, int nseq,NT_node **lu);
NT_node insert_root(NT_node p, float diff);
float calc_root_mean(NT_node root, float *maxdist, int neq, NT_node **lu);
NT_node reroot(NT_node ptree, int nseq, int ntotal, int nnodes, NT_node **lu);


Alignment *seq2seq_chain (Alignment *A,Alignment *B, char *arg);

float display_avg_bootstrap ( NT_node T);
float tree2tot_dist ( NT_node T, int mode);
int tree2n_branches(NT_node T, int mode);
int **display_tree_from_node (NT_node T, int up, int down, int **array);
NT_node tree2node ( char *name, NT_node T);
NT_node * tree2node_list (NT_node T, NT_node *L);
NT_node tree2root ( NT_node T);
int new_tree_sort ( char *name, NT_node T);


NT_node split2tree ( NT_node RT,Sequence *LIST, char *param);
NT_node * read_tree_list (Sequence *S);

int count_groups( Sequence *S, char *s);

Split ** count_splits( NT_node RT, Sequence *S, char *s);
NT_node *treelist2prune_treelist (Sequence *S, Sequence *TS, FILE *out);
int** treelist2groups (Sequence *S, Sequence *ST, char *depth, FILE *out);
int treelist2splits (Sequence *S, Sequence *ST);
int treelist2leafgroup ( Sequence *S, Sequence *TS, char *taxon);
int ***tree2dist ( NT_node T, Sequence *S, int ***d);
int treelist2frame (Sequence *S, Sequence *TS);
int** treelist2lti ( Sequence *S, Sequence *TS, int nb, FILE *out);

float simple_tree_cmp (NT_node T1, NT_node T2,Sequence *S, int mode);

int treelist2dmat ( Sequence *S);
NT_node new_declare_tree_node ();
int count_tree_groups( Sequence *LIST, char *group_file);
int node_sort ( char *name, NT_node T);
int    treelist2n (NT_node *L);
int ** treelist2avg_treecmp (NT_node *L, char *file);
NT_node treelist2bootstrap ( NT_node *L, char *file);
NT_node treelist2filtered_bootstrap ( NT_node *L, char *file, int **score,float f);

Sequence * treelist2seq ( Sequence *S);
Sequence * treelist2sub_seq ( Sequence *S, int f);

int treelist_file2consense (char *tree_file, char *outtree, char *outfile);
/* General purpose header file - rf 12/90 */

#ifndef _H_general
#define _H_general



#define pint int			/* cast ints in printf statements as pint */
typedef int Boolean;			/* Is already defined in THINK_C */

#undef TRUE						
#undef FALSE
#define TRUE 1
#define FALSE 0

#define EOS '\0'				/* End-Of-String */
#define MAXLINE 512			/* Max. line length */


#endif /* ifndef _H_general */
