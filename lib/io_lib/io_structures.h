/**
 * \file io_structures.h
 *
 * Defines structs like ::Seqeunce, ::Weights, ::Alignment and ::Template.
 */

#include "define_header.h"

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


typedef struct
{
	char *seg_name; /// Name of chromosome, contig...
	char strand; 	/// The strand of the sequence: +/-
	unsigned int start; /// Start point (always from + point of view) start at 0
	unsigned int end;	///last nucleotide included (always from + point of view)
	unsigned int seg_len; ///The length of the chromosome
}
Genomic_info;



/**
 * Container for several sequences and their attributes.
 *
 * \note Not all members have a documentation
 */
struct Sequence
{
  char **file;          	/**< \c file[\c nseq ][\c filenamelength ] Name of the file where each sequence is taken from */
  char **blastfile;          	/**Contains the name of the blastfiles*/
  
  char **seq_comment;   	/**< \c seq_comment[\c nseq ][\c LONG_STRING ] comment read in the fasta file */
  char **aln_comment;  		/*id*/
  char **seq;          		/**< \c seq[\c nseq ][\c len ] Actual sequences */
  int *len;            		/**< \c len[\c nseq ] length of each sequence */
  int max_len;         		/**< Lenght of the longest sequence */
  int min_len;         		/**< Length of the shortest sequence */
  int nseq;            		/**< Number of sequences */
  int max_nseq;        		/*Maximum number of sequences in the datastruct*/
  char **name;         		/**< \c name[\c nseq ][\c MAXNAMELEN ] Names of the sequences */
  int **dc;         		/*coordinates on the disk. Coordinates set if seq[i]==NULL*/
  
  struct Constraint_list *CL; /**< Points to the ::Constraint_list */
  int contains_gap;   		/**< Set to 1 if gaps should be kept */
  char *type;         		/**< PROTEIN, DNA, RNA */
  Weights *W;         		/**< Associated ::Weight object */
  char template_file[FILENAMELEN+1];
  struct Template **T;		/**< \c T[\c nseq ] Pointer to ::Template for each sequence */
  char *blastdb;
  struct Sequence *blastdbS;
  struct Sequence *MasterS;
  Genomic_info *genome_co; //safes genome_coordinates


};
typedef struct Sequence Sequence;

/**
 * Any sort of Template like PDB structure, Profile or Secondary structure.
 *
 * The Template structure looks a little confusing on first sight. It consists of pointers
 * of type ::X_template to different types of templates, named after their template type (*P, *R, *R etc).
 *
 * The ::X_template again contains pointers to all kinds of templates, but this time with each specified type
 * like ::P_template, ::S_template, ::R_template and so on. These specific template structures contain the
 * actual information, like a PDB identifier, an alignment or another sequence.
 *
 * \attention When alternative templates are given for one sequence, the first one superseeds all the others.
 * \todo Find out whether you can use more than one of these pointers at the same time to specify several different templates.
 * \todo See what happens to the template files you specify in :: ... where?
 */
struct Template
{
  char seq_type[10];	 /**< String containing information on which templates are used.
                              Looks like "P..S.......", for example, where dots are actually white spaces */
  struct X_template *P;  /**< PDB structure */
  struct X_template *F;  /**< RNA secondary structure */
  struct X_template *S;  /**< sequence */
  struct X_template *R;  /**< Profile */
  struct X_template *G;  /**< Genomic structure */
  struct X_template *T;  /**< transmembrane */
  struct X_template *E;  /**< secondary structure */
  struct X_template *U;  /**< Unicode, strings */

  struct X_template *RB;  /**< ? */
};
typedef struct Template Template;


/**
 * See ::Template
 */
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

/**
 * See ::Template
 */
struct P_template
{
  char pdb_id[100]; /**< PDB identifier */
};
typedef struct P_template P_template;

/**
 * RNA Secondary structure, see ::Template
 */
struct F_template
{
  int l;
};
typedef struct F_template F_template;

/**
 * See ::Template
 */
struct S_template
{
  Sequence *S; /**< Sequence object */
};
typedef struct S_template S_template;


/**
 * Profile associated with a sequence, see ::Template
 */
struct R_template
{
  struct Alignment *A; /**< Alignment */
};
typedef struct R_template R_template;



/**
 * Genomic Information, see ::Template
 */
struct G_template
{
  Sequence *S; /**< Sequence object */
};
typedef struct G_template G_template;

/**
 * See ::Template
 */
struct T_template
{
  Sequence *S; /**< Sequence object */
};
typedef struct T_template T_template;

/**
 * See ::Template
 */
struct E_template
{
  Sequence *S; /**< Sequence object */
};
typedef struct E_template E_template;

/**
 * See ::Template
 */
struct U_template
{
  int *list; /**< Int aray */
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

/*Trees*/
      Alignment *Tree;
      char *tname;
      int  **RepColList;//last item set to -1
      char **dmF_list;
/*Score*/
    int       nsites;  
    double ** dm;   
    int   *  score_seq;
    int   ** score_res;
    int      score_aln;
    int      score;
    int      ibit;
    int      cpu;
    int      finished;

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
