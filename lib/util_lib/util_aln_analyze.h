
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
