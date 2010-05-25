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


   
/*Score*/
    int *  score_seq;
    int ** score_res;
    int score_aln;
    int score;
    	
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
