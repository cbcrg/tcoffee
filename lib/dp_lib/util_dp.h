/*pair wise aln implementations*/
int myers_miller_pair_wise (Alignment *A, int *ns, int **l_s,Constraint_list *CL);
int diff (Alignment *A, int *ns, int **ls, int s1, int M,int s2, int N , int tb, int te, Constraint_list *CL, int **pos);
int gotoh_pair_wise    (Alignment *A,int*ns, int **l_s,Constraint_list *CL);

int sw_pair_wise (Alignment *A, int gop,int gep, int scale,int*ns, int **l_s,Constraint_list *CL, int maximise );
Alignment ** t_coffee_lalign   (Constraint_list *CL, int scale, int penalty, int maximise,Sequence *S,int sw_t, int sw_l, int sw_z, int *sw_n, int sw_io);


Alignment * get_best_local_aln     ( Alignment *A,Constraint_list *CL,int gop, int gep, int sw_t, int sw_l, int sw_z, int mode);
Alignment *   get_best_ol_local_aln ( Alignment *IN, Constraint_list*CL,int gop, int gep, int sw_t, int sw_l, int sw_z, int mode);
Alignment *  get_best_nol_local_aln ( Alignment *IN, Constraint_list*CL,int gop, int gep, int sw_t, int sw_l, int sw_z, int mode);

int    get_start_point ( Alignment *A,Constraint_list**CL,int gop, int gep, int T, int *first, int *last);

double compute_penalty ( Constraint_list *CL, char *mode, int len);
double compute_scale   ( Constraint_list *CL, char *mode, int len);
int evaluate_penalty (Alignment *A, Constraint_list *CL, int *scale,char *scale_mode, int *penalty, char *penalty_mode, int len_seq); 
Alignment *  add_seq2aln           (Constraint_list *CL,int gop,int gep, Alignment *a_list,int maximise,Sequence  *S, int sw_t, int sw_l, int sw_z, int *sw_n, int sw_io);
