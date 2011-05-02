
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
      int *master; //Sequences used as master sequences
      int o2a_byte; // number of one to all provided in one go.
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
      int   reverse_seq;//Used for HoT
      int   extend_seq; //Used for RNA or Promoter Alignments
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
  int extend_seq;
  int reverse_seq;
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
Constraint_list *  empty_constraint_list (Constraint_list *CL);
Constraint_list *  unfreeze_constraint_list (Constraint_list *CL);
Constraint_list *  freeze_constraint_list (Constraint_list *CL);
Constraint_list *  undump_constraint_list (Constraint_list *CL, char *file);
int   dump_constraint_list (Constraint_list *CL, char *file,char *mode);
int   safe_dump_constraint_list (Constraint_list *CL, char *file,char *mode, Sequence *RS);
FILE* display_constraint_list (Constraint_list *CL, FILE *fp, char *tag);


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
Sequence *precompute_blast_db (Sequence *S, char **ml, int n);
int read_seq_in_list ( char *fname,  int *nseq, char ***sequences, char ***seq_name, Genomic_info **genome_co);

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
Constraint_list * expand_constraint_list (Constraint_list *CL, int T);
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
/*********************************************************************/
/*                                                                   */
/*                        SCALED CONSISTENCY                          */
/*                                                                   */
/*                                                                   */
/*********************************************************************/
Constraint_list * plib_msa (Constraint_list *CL);
int cl2worst_seq (Constraint_list *CL, int *list, int n);
Constraint_list *add_seq2cl(int s, Constraint_list *CL);
