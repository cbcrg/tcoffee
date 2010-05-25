/*By convention, 0: START, 1: END*/

struct Hmm
{
  
  double freeT;           /*Free transition*/
  double forbiden;        /*Forbiden transition*/
  int start;              /*start, by convention: 0*/
  int end;                /*end, by convention: 1*/

  int nS;                 /*Number of states*/
  int order;
  struct HmmState  **S;   /*State List*/

  /*Bounded HMM*/
  double **T;             /*Transition matrix*/
  int **fromM;            /*For any sate s, fromM[0]->number of states leading to s*/
  int **toM;              /*For any sate s, toM[0]->number of s can go to*/
                          /*End and Start are NOT included in toM and FromM*/
  
  
};
typedef struct Hmm Hmm;

struct HmmAln
{
  Hmm *H;
  int *state_list;
};
typedef struct HmmAln HmmAln;

typedef double (*Generic_em_func)(struct Hmm*, struct HmmState*, int);

struct HmmState
{
char name[100];
int state;
int DJ;
int DI;

  /*Pair HMM*/
double em;
Col_cost_func em_func;

  /*Linear HMM*/
double *em2;
Generic_em_func em_func2;
int nT;
struct StateTrans **T;
};
typedef struct HmmState HmmState;

struct StateTrans
{
  char name[101];
  double tr;
};
typedef struct StateTrans StateTrans;

struct MatState
{
  int i;
  int j;
  int st;
  int pst;
  double sc;
  struct MatState *n;
  struct MatState *p;
  struct MatState *m; /*memory*/
  struct MatState *s; /*memory of the start point*/
  /*Heap Mamagement: Never copy*/
  struct MatState *Hn;/*Heap Next*/ 
  struct MatState *Hp;/*Heap Previous*/ 
  
  struct MatState *Mn;/*Marked Heap section*/ 
  struct MatState *Mp;/*Marked Heap Section*/ 
  int free;
};
typedef struct MatState MatState;


/*********************************************************************************/
/*                                                                               */
/*                                                                               */
/*                     simple HMM: Viterbi                                       */
/*                                                                               */
/*                                                                               */
/*********************************************************************************/
int seq_viterbi_pair_wise (Alignment *A,int*ns, int **ls,Constraint_list *CL);

/*********************************************************************************/
/*                                                                               */
/*                                                                               */
/*                     HMM: Viterbi                                              */
/*                                                                               */
/*                                                                               */
/*********************************************************************************/

int viterbi_pair_wise_OLD (Alignment *A,int*ns, int **ls,Constraint_list *CL);
Alignment * viterbipath2aln (Alignment *A, int *ns,int **ls,int *tb, Hmm *H);
double*** viterbi_hmm (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL);
int * viterbi2path (int l1,int l2, Hmm *H, double ***M);
/*********************************************************************************/
/*                                                                               */
/*                                                                               */
/*                      HMM modeling                                             */
/*                                                                               */
/*                                                                               */
/*********************************************************************************/
int viterbiL_pair_wise (Alignment *A,int*ns, int **ls,Constraint_list *CL);
MatState* RviterbiL_hmm (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL,MatState *S, MatState *E);
MatState* viterbiL_hmm (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL, MatState *S, MatState *E);

int viterbiD_pair_wise (Alignment *A,int*ns, int **ls,Constraint_list *CL);
double lu_RviterbiD_hmm (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL,MatState *S, MatState *E, int **seg_list);
MatState* RviterbiD_hmm (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL,MatState *S, MatState *E, int **seg_list);
MatState* viterbiD_hmm (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL, MatState *S, MatState *E, int **seg_list);
int **seglist2table ( int **seglist,int l1, int l2);
int *seglist2line ( int **seglist, int line, int start, int end);
int * traceback (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL,MatState *S, MatState *E, int **seg_list);

int viterbiDGL_pair_wise (Alignment *A,int*ns, int **ls,Constraint_list *CL);
double lu_RviterbiDGL_hmm (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL,MatState *S, MatState *E, int **seg_list);
MatState* RviterbiDGL_hmm (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL,MatState *S, MatState *E, int **seg_list);
MatState* viterbiDGL_hmm (Alignment *A,int *ns, int **ls, Hmm *H, Constraint_list *CL,MatState *S, MatState *E, int **seg_list);


/*********************************************************************************/
/*                                                                               */
/*                                                                               */
/*                      HMM modeling                                             */
/*                                                                               */
/*                                                                               */
/*********************************************************************************/
int MatStateAreIdentical (MatState*I, MatState*O);
int MaxDeltaMatState (MatState*I, MatState*O);
int MinDeltaMatState (MatState*I, MatState*O);

MatState * ManageMatState(int Mode, MatState *C);
MatState* CopyMatState ( MatState*I, MatState*O);

Hmm* read_hmm(char *file);
Hmm* declare_hmm(int n);
Hmm* free_Hmm(Hmm*H);
void DisplayHmm ( Hmm *H);

/*********************************************************************************/
/*                                                                               */
/*                                                                               */
/*                      HMM Models                                               */
/*                                                                               */
/*                                                                               */
/*********************************************************************************/
Hmm* define_two_mat_model(Constraint_list *CL);
Hmm* define_probcons_model(Constraint_list *CL);
Hmm* define_simple_model(Constraint_list *CL);

Hmm * bound_hmm ( Hmm *H);
