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
