#ifndef COFFEE_DEFINES_H
#define COFFEE_DEFINES_H

/*DEBUGGING*/
/*#include "mshell.h"*/
/*MEMORY MANAGEMENT*/
#include <float.h>
#define MY_EPSILON 1000*DBL_EPSILON
//Maximum number of tries for interactibve things
#define MAX_N_TRIES 3
#define MAX_NSEQ_4_DISPLAY 50
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
#define MEMSET 1
#define NOMEMSET -1
#define RESIZE 1
#define NORESIZE -1

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

#define PROFILE 66
#define SEQUENCE 67

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

#define COMMENT_SIZE      1000
#define MAXNAMES          1000
#define FILENAMELEN 	  500            /* Max. file name length */
#define MAX_N_PARAM       2000
#define MAX_PARAM_LEN     200
#define MAX_LINE_LENGTH   10000
#define ALN_LINE_LENGTH   50
#define SHORT_STRING      10
#define STRING            300
#define LONG_STRING       1000
#define VERY_LONG_STRING  10000
#define SUPER_LONG_STRING 100000
#define MEGA_LONG_STRING  1000000


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
#define declare_name(x) (x=(char*)vcalloc (MAX(FILENAMELEN,L_tmpnam)+1, sizeof (char))) 
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
#define EMAIL "cedric.notredame@gmail.com"
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


#define PATCH_PRF ""
//This variable is set so as to compensate a bug in Clustal-Omega
#endif // -- COFFEE_DEFINES_H
