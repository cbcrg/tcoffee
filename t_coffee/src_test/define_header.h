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
/*               PROGRAM PATH                  */


//ERROR MESSAGES


#define ADDRESS_BUILT_IN "built_in"
#define PROGRAM_BUILT_IN "t_coffee"
#define TEST_WWWSITE_4_TCOFFEE "www.google.com"


//TclinkdbStart

#define TCOFFEE_4_TCOFFEE "t_coffee"
#define TCOFFEE_type "sequence_multiple_aligner"
#define TCOFFEE_ADDRESS "http://www.tcoffee.org"
#define TCOFFEE_language "C"
#define TCOFFEE_language2 "C"
#define TCOFFEE_source "http://www.tcoffee.org/Packages/T-COFFEE_distribution.tar.gz"
#define TCOFFEE_update_action "always"
#define TCOFFEE_mode "tcoffee,mcoffee,rcoffee,expresso,3dcoffee"
#define CLUSTALW2_4_TCOFFEE "clustalw2"
#define CLUSTALW2_type "sequence_multiple_aligner"
#define CLUSTALW2_ADDRESS "http://www.clustal.org"
#define CLUSTALW2_language "C++"
#define CLUSTALW2_language2 "CXX"
#define CLUSTALW2_source "http://www.clustal.org/download/2.0.10/clustalw-2.0.10-src.tar.gz"
#define CLUSTALW2_mode "mcoffee,rcoffee"
#define CLUSTALW2_version "2.0.10"
#define CLUSTALW_4_TCOFFEE "clustalw"
#define CLUSTALW_type "sequence_multiple_aligner"
#define CLUSTALW_ADDRESS "http://www.clustal.org"
#define CLUSTALW_language "C"
#define CLUSTALW_language2 "C"
#define CLUSTALW_source "http://www.clustal.org/download/1.X/ftp-igbmc.u-strasbg.fr/pub/ClustalW/clustalw1.82.UNIX.tar.gz"
#define CLUSTALW_mode "mcoffee,rcoffee"
#define CLUSTALW_version "1.82"
#define DIALIGNT_4_TCOFFEE "dialign-t"
#define DIALIGNT_type "sequence_multiple_aligner"
#define DIALIGNT_ADDRESS "http://dialign-tx.gobics.de/"
#define DIALIGNT_DIR "/usr/share/dialign-tx/"
#define DIALIGNT_language "C"
#define DIALIGNT_language2 "C"
#define DIALIGNT_source "http://dialign-tx.gobics.de/DIALIGN-TX_1.0.2.tar.gz"
#define DIALIGNT_mode "mcoffee"
#define DIALIGNT_binary "dialign-t"
#define DIALIGNT_version "1.0.2"
#define DIALIGNTX_4_TCOFFEE "dialign-tx"
#define DIALIGNTX_type "sequence_multiple_aligner"
#define DIALIGNTX_ADDRESS "http://dialign-tx.gobics.de/"
#define DIALIGNTX_DIR "/usr/share/dialign-tx/"
#define DIALIGNTX_language "C"
#define DIALIGNTX_language2 "C"
#define DIALIGNTX_source "http://dialign-tx.gobics.de/DIALIGN-TX_1.0.2.tar.gz"
#define DIALIGNTX_mode "mcoffee"
#define DIALIGNTX_binary "dialign-tx"
#define DIALIGNTX_version "1.0.2"
#define POA_4_TCOFFEE "poa"
#define POA_type "sequence_multiple_aligner"
#define POA_ADDRESS "http://www.bioinformatics.ucla.edu/poa/"
#define POA_language "C"
#define POA_language2 "C"
#define POA_source "http://downloads.sourceforge.net/poamsa/poaV2.tar.gz"
#define POA_DIR "/usr/share/"
#define POA_FILE1 "blosum80.mat"
#define POA_mode "mcoffee"
#define POA_binary "poa"
#define POA_version "2.0"
#define PROBCONS_4_TCOFFEE "probcons"
#define PROBCONS_type "sequence_multiple_aligner"
#define PROBCONS_ADDRESS "http://probcons.stanford.edu/"
#define PROBCONS_language2 "CXX"
#define PROBCONS_language "C++"
#define PROBCONS_source "http://probcons.stanford.edu/probcons_v1_12.tar.gz"
#define PROBCONS_mode "mcoffee"
#define PROBCONS_binary "probcons"
#define PROBCONS_version "1.12"
#define MAFFT_4_TCOFFEE "mafft"
#define MAFFT_type "sequence_multiple_aligner"
#define MAFFT_ADDRESS "http://align.bmr.kyushu-u.ac.jp/mafft/online/server/"
#define MAFFT_language "C"
#define MAFFT_language "C"
#define MAFFT_source "http://align.bmr.kyushu-u.ac.jp/mafft/software/mafft-6.603-with-extensions-src.tgz"
#define MAFFT_windows "http://align.bmr.kyushu-u.ac.jp/mafft/software/mafft-6.603-mingw.tar"
#define MAFFT_mode "mcoffee,rcoffee"
#define MAFFT_binary "mafft.tar.gz"
#define MAFFT_version "6.603"
#define MUSCLE_4_TCOFFEE "muscle"
#define MUSCLE_type "sequence_multiple_aligner"
#define MUSCLE_ADDRESS "http://www.drive5.com/muscle/"
#define MUSCLE_language "C++"
#define MUSCLE_language2 "GPP"
#define MUSCLE_source "http://www.drive5.com/muscle/downloads3.7/muscle3.7_src.tar.gz"
#define MUSCLE_windows "http://www.drive5.com/muscle/downloads3.7/muscle3.7_win32.zip"
#define MUSCLE_linux "http://www.drive5.com/muscle/downloads3.7/muscle3.7_linux_ia32.tar.gz"
#define MUSCLE_mode "mcoffee,rcoffee"
#define MUSCLE_version "3.7"
#define MUS4_4_TCOFFEE "mus4"
#define MUS4_type "sequence_multiple_aligner"
#define MUS4_ADDRESS "http://www.drive5.com/muscle/"
#define MUS4_language "C++"
#define MUS4_language2 "GPP"
#define MUS4_source "http://www.drive5.com/muscle/muscle4.0_src.tar.gz"
#define MUS4_mode "mcoffee,rcoffee"
#define MUS4_version "4.0"
#define PCMA_4_TCOFFEE "pcma"
#define PCMA_type "sequence_multiple_aligner"
#define PCMA_ADDRESS "ftp://iole.swmed.edu/pub/PCMA/"
#define PCMA_language "C"
#define PCMA_language2 "C"
#define PCMA_source "ftp://iole.swmed.edu/pub/PCMA/pcma.tar.gz"
#define PCMA_mode "mcoffee"
#define PCMA_version "1.0"
#define KALIGN_4_TCOFFEE "kalign"
#define KALIGN_type "sequence_multiple_aligner"
#define KALIGN_ADDRESS "http://msa.cgb.ki.se"
#define KALIGN_language "C"
#define KALIGN_language2 "C"
#define KALIGN_source "http://msa.cgb.ki.se/downloads/kalign/current.tar.gz"
#define KALIGN_mode "mcoffee"
#define KALIGN_version "1.0"
#define AMAP_4_TCOFFEE "amap"
#define AMAP_type "sequence_multiple_aligner"
#define AMAP_ADDRESS "http://bio.math.berkeley.edu/amap/"
#define AMAP_language "C++"
#define AMAP_language2 "CXX"
#define AMAP_source "http://amap-align.googlecode.com/files/amap.2.0.tar.gz"
#define AMAP_mode "mcoffee"
#define AMAP_version "2.0"
#define PRODA_4_TCOFFEE "proda"
#define PRODA_type "sequence_multiple_aligner"
#define PRODA_ADDRESS "http://proda.stanford.edu"
#define PRODA_language "C++"
#define PRODA_language2 "CXX"
#define PRODA_source "http://proda.stanford.edu/proda_1_0.tar.gz"
#define PRODA_mode "mcoffee"
#define PRODA_version "1.0"
#define FSA_4_TCOFFEE "fsa"
#define FSA_type "sequence_multiple_aligner"
#define FSA_ADDRESS "http://fsa.sourceforge.net/"
#define FSA_language "C++"
#define FSA_language2 "CXX"
#define FSA_source "http://sourceforge.net/projects/fsa/files/fsa-1.15.3.tar.gz/download/"
#define FSA_mode "mcoffee"
#define FSA_version "1.15.3"
#define PRANK_4_TCOFFEE "prank"
#define PRANK_type "sequence_multiple_aligner"
#define PRANK_ADDRESS "http://www.ebi.ac.uk/goldman-srv/prank/"
#define PRANK_language "C++"
#define PRANK_language2 "CXX"
#define PRANK_source "http://www.ebi.ac.uk/goldman-srv/prank/src/prank/prank.src.100303.tgz"
#define PRANK_mode "mcoffee"
#define PRANK_version "100303"
#define SAP_4_TCOFFEE "sap"
#define SAP_type "structure_pairwise_aligner"
#define SAP_ADDRESS "http://mathbio.nimr.mrc.ac.uk/wiki/Software"
#define SAP_language "C"
#define SAP_language2 "C"
#define SAP_source "http://mathbio.nimr.mrc.ac.uk/download/sap-1.1.1.tar.gz"
#define SAP_mode "expresso,3dcoffee"
#define SAP_version "1.1.1"
#define TMALIGN_4_TCOFFEE "TMalign"
#define TMALIGN_type "structure_pairwise_aligner"
#define TMALIGN_ADDRESS "http://zhang.bioinformatics.ku.edu/TM-align/TMalign.f"
#define TMALIGN_language "Fortran"
#define TMALIGN_language2 "Fortran"
#define TMALIGN_source "http://zhang.bioinformatics.ku.edu/TM-align/TMalign.f"
#define TMALIGN_linux "http://zhang.bioinformatics.ku.edu/TM-align/TMalign_32.gz"
#define TMALIGN_mode "expresso,3dcoffee"
#define TMALIGN_version "1.0"
#define MUSTANG_4_TCOFFEE "mustang"
#define MUSTANG_type "structure_pairwise_aligner"
#define MUSTANG_ADDRESS "http://www.cs.mu.oz.au/~arun/mustang"
#define MUSTANG_language "C++"
#define MUSTANG_language2 "CXX"
#define MUSTANG_source "http://ww2.cs.mu.oz.au/~arun/mustang/mustang_v3.2.1.tgz"
#define MUSTANG_mode "expresso,3dcoffee"
#define MUSTANG_version "3.2.1"
#define LSQMAN_4_TCOFFEE "lsqman"
#define LSQMAN_type "structure_pairwise_aligner"
#define LSQMAN_ADDRESS "empty"
#define LSQMAN_language "empty"
#define LSQMAN_language2 "empty"
#define LSQMAN_source "empty"
#define LSQMAN_update_action "never"
#define LSQMAN_mode "expresso,3dcoffee"
#define ALIGN_PDB_4_TCOFFEE "align_pdb"
#define ALIGN_PDB_type "structure_pairwise_aligner"
#define ALIGN_PDB_ADDRESS "empty"
#define ALIGN_PDB_language "empty"
#define ALIGN_PDB_language2 "empty"
#define ALIGN_PDB_source "empty"
#define ALIGN_PDB_update_action "never"
#define ALIGN_PDB_mode "expresso,3dcoffee"
#define FUGUE_4_TCOFFEE "fugueali"
#define FUGUE_type "structure_pairwise_aligner"
#define FUGUE_ADDRESS "http://www-cryst.bioc.cam.ac.uk/fugue/download.html"
#define FUGUE_language "empty"
#define FUGUE_language2 "empty"
#define FUGUE_source "empty"
#define FUGUE_update_action "never"
#define FUGUE_mode "expresso,3dcoffee"
#define DALILITEc_4_TCOFFEE "dalilite.pl"
#define DALILITEc_type "structure_pairwise_aligner"
#define DALILITEc_ADDRESS "built_in"
#define DALILITEc_ADDRESS2 "http://www.ebi.ac.uk/Tools/webservices/services/dalilite"
#define DALILITEc_language "Perl"
#define DALILITEc_language2 "Perl"
#define DALILITEc_source "empty"
#define DALILITEc_update_action "never"
#define DALILITEc_mode "expresso,3dcoffee"
#define PROBCONSRNA_4_TCOFFEE "probconsRNA"
#define PROBCONSRNA_type "RNA_multiple_aligner"
#define PROBCONSRNA_ADDRESS "http://probcons.stanford.edu/"
#define PROBCONSRNA_language "C++"
#define PROBCONSRNA_language2 "CXX"
#define PROBCONSRNA_source "http://probcons.stanford.edu/probconsRNA.tar.gz"
#define PROBCONSRNA_mode "mcoffee,rcoffee"
#define PROBCONSRNA_version "1.0"
#define CONSAN_4_TCOFFEE "sfold"
#define CONSAN_type "RNA_pairwise_aligner"
#define CONSAN_ADDRESS "http://selab.janelia.org/software/consan/"
#define CONSAN_language "empty"
#define CONSAN_language2 "empty"
#define CONSAN_source "empty"
#define CONSAN_update_action "never"
#define CONSAN_mode "rcoffee"
#define RNAPLFOLD_4_TCOFFEE "RNAplfold"
#define RNAPLFOLD_type "RNA_secondarystructure_predictor"
#define RNAPLFOLD_ADDRESS "http://www.tbi.univie.ac.at/~ivo/RNA/"
#define RNAPLFOLD_language "C"
#define RNAPLFOLD_language2 "C"
#define RNAPLFOLD_source "http://www.tbi.univie.ac.at/~ivo/RNA/ViennaRNA-1.7.2.tar.gz"
#define RNAPLFOLD_mode "rcoffee,"
#define RNAPLFOLD_version "1.7.2"
#define PHYLIP_4_TCOFFEE "retree"
#define PHYLIP_type "RNA_secondarystructure_predictor"
#define PHYLIP_ADDRESS "http://evolution.gs.washington.edu/phylip/"
#define PHYLIP_language "C"
#define PHYLIP_language2 "C"
#define PHYLIP_source "http://evolution.gs.washington.edu/phylip/download/phylip-3.69.tar.gz"
#define PHYLIP_mode "trmsd,"
#define PHYLIP_version "3.69"
#define HMMTOP_4_TCOFFEE "hmmtop"
#define HMMTOP_type "protein_secondarystructure_predictor"
#define HMMTOP_ADDRESS "www.enzim.hu/hmmtop/"
#define HMMTOP_language "C"
#define HMMTOP_language2 "C"
#define HMMTOP_source "empty"
#define HMMTOP_update_action "never"
#define HMMTOP_mode "tcoffee"
#define GOR4_4_TCOFFEE "gorIV"
#define GOR4_type "protein_secondarystructure_predictor"
#define GOR4_ADDRESS "http://mig.jouy.inra.fr/logiciels/gorIV/"
#define GOR4_language "C"
#define GOR4_language2 "C"
#define GOR4_source "http://mig.jouy.inra.fr/logiciels/gorIV/GOR_IV.tar.gz"
#define GOR4_update_action "never"
#define GOR4_mode "tcoffee"
#define EBIWUBLASTc_4_TCOFFEE "wublast.pl"
#define EBIWUBLASTc_type "protein_homology_predictor"
#define EBIWUBLASTc_ADDRESS "built_in"
#define EBIWUBLASTc_ADDRESS2 "http://www.ebi.ac.uk/Tools/webservices/services/wublast"
#define EBIWUBLASTc_language "Perl"
#define EBIWUBLASTc_language2 "Perl"
#define EBIWUBLASTc_source "empty"
#define EBIWUBLASTc_update_action "never"
#define EBIWUBLASTc_mode "psicoffee,expresso,accurate"
#define EBIBLASTPGPc_4_TCOFFEE "blastpgp.pl"
#define EBIBLASTPGPc_type "protein_homology_predictor"
#define EBIBLASTPGPc_ADDRESS "built_in"
#define EBIBLASTPGPc_ADDRESS2 "http://www.ebi.ac.uk/Tools/webservices/services/blastpgp"
#define EBIBLASTPGPc_language "Perl"
#define EBIBLASTPGPc_language2 "Perl"
#define EBIBLASTPGPc_source "empty"
#define EBIBLASTPGPc_update_action "never"
#define EBIBLASTPGPc_mode "psicoffee,expresso,accurate"
#define NCBIWEBBLAST_4_TCOFFEE "blastcl3"
#define NCBIWEBBLAST_type "protein_homology_predictor"
#define NCBIWEBBLAST_ADDRESS "ftp://ftp.ncbi.nih.gov/blast/executables/LATEST"
#define NCBIWEBBLAST_language "C"
#define NCBIWEBBLAST_language2 "C"
#define NCBIWEBBLAST_source "empty"
#define NCBIWEBBLAST_update_action "never"
#define NCBIWEBBLAST_mode "psicoffee,expresso,3dcoffee"
#define blastall_4_TCOFFEE "blastall"
#define blastall_type "protein_homology_predictor"
#define blastall_ADDRESS "ftp://ftp.ncbi.nih.gov/blast/executables/LATEST"
#define blastall_language "C"
#define blastall_language2 "C"
#define blastall_source "empty"
#define blastall_update_action "never"
#define blastall_mode "psicoffee,expresso,3dcoffee"
#define NCBIBLAST_4_TCOFFEE "legacy_blast.pl"
#define NCBIBLAST_type "protein_homology_predictor"
#define NCBIBLAST_ADDRESS "ftp://ftp.ncbi.nih.gov/blast/executables/LATEST"
#define NCBIBLAST_language "C"
#define NCBIBLAST_language2 "C"
#define NCBIBLAST_source "empty"
#define NCBIBLAST_update_action "never"
#define NCBIBLAST_mode "psicoffee,expresso,3dcoffee"
#define SOAPLITE_4_TCOFFEE "SOAP::Lite"
#define SOAPLITE_type "library"
#define SOAPLITE_ADDRESS "http://cpansearch.perl.org/src/MKUTTER/SOAP-Lite-0.710.08/Makefile.PL"
#define SOAPLITE_language "Perl"
#define SOAPLITE_language2 "Perl"
#define SOAPLITE_source "empty"
#define _update_action "never"
#define SOAPLITE_mode "none"
#define XMLSIMPLE_4_TCOFFEE "XML::Simple"
#define XMLSIMPLE_type "library"
#define XMLSIMPLE_ADDRESS "http://search.cpan.org/~grantm/XML-Simple-2.18/lib/XML/Simple.pm"
#define XMLSIMPLE_language "Perl"
#define XMLSIMPLE_language2 "Perl"
#define XMLSIMPLE_source "empty"
#define XMLSIMPLE_mode "psicoffee,expresso,accurate"
//TclinkdbEnd
/*New Methods*/
/********************************************/
/*            Various Methoids              */
/********************************************/
#define METHODS_4_TCOFFEE  "~/.t_coffee/methods/"
#define METHOD_4_MSA_WEIGHTS "petra_weight"
/********************************************/
/*            SEQAN LIBRARY                  */
/********************************************/
#define SEQAN_TCOFFEE_4_TCOFFEE "seqan_tcoffee"
/********************************************/
/*            REFORMATING  AND UTILITIES                 */
/********************************************/
#define WGET_4_TCOFFEE "wget"
#define WGET_ADDRESS "http://www.gnu.org/software/wget/"

#define CURL_4_TCOFFEE "curl"
#define CURL_ADDRESS "http://curl.haxx.se/"

#define SEQ_REFORMAT_4_TCOFFEE          "seq_reformat"
#define PS2PDF                          "ps2pdf"
#define EXTRACT_FROM_PDB_4_TCOFFEE      "extract_from_pdb"
#define BLAST_ALN2FASTA_ALN             "blast_aln2fasta_aln.pl"
#define FASTA_ALN2FASTA_ALN_UNIQUE_NAME "fasta_aln2fasta_aln_unique_name.pl"
#define MSF_ALN2FASTA_ALN               "msf_aln2fasta_aln.pl"
#define SEQ2MSA_WEIGHT        "seq2msa_weight"
/********************************************/
/*            DEPRECATED DEF                */
/********************************************/
//Deprecated definitions
#define SIB_BLAST_4_TCOFFEE "blastall.remote"
#define LOCAL_BLAST_4_TCOFFEE "blastall"
#define BLAST_DB_4_TCOFFEE "nr"
#define NCBI_BLAST_4_TCOFFEE ""
/********************************************/
/*            PARAMETER_FILE                */
/********************************************/




/*               PARAMETER FILES               */
#define COLOR_FILE         "seq_reformat.color"
/*This file specifies the 10 colors available to seq_reformat.
If the file is not on the system, hard coded defaults will be used.
The format is as follow:

-------------------------------------------------------------------------------------------
<Your comments (as many lines as needed >
*
<number in the range 0-9> <HTML code> <R value (float)> <G value (float)> <B value (float)>
-------------------------------------------------------------------------------------------
the RGB values are used for the post-script generation, the html code is used in html documents.
*/
#define DATE "Tue Jan 11 23:44:37 CET 2011"
#define PROGRAM "T-COFFEE"
#define VERSION "Version_8.98"
#define AUTHOR "Cedric Notredame "
#define DISTRIBUTION_ADDRESS "www.tcoffee.org/Packages/"
