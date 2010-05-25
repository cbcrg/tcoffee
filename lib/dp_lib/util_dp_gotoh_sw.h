
int gotoh_pair_wise_lalign ( Alignment *A, int*ns, int **l_s,Constraint_list *CL);
Constraint_list * undefine_sw_aln ( Alignment *A, Constraint_list *CL);
Constraint_list * undefine_sw_pair ( Constraint_list *CL, int s1, int r1, int s2, int r2);
int sw_pair_is_defined ( Constraint_list *CL, int s1, int r1, int s2, int r2);


int gotoh_pair_wise_sw (Alignment *A,int*ns, int **l_s,Constraint_list *CL);

Alignment * get_best_local_aln ( Alignment *IN,Constraint_list *CL,int gop, int gep, int sw_t, int sw_l, int sw_z, int greedy);
Alignment * get_best_nol_local_aln ( Alignment *IN, Constraint_list *CL, int gop, int gep,int sw_t,int sw_l, int sw_z, int mode);
double compute_penalty   (Constraint_list *CL, char *mode, int len);
double compute_scale ( Constraint_list *CL,char *mode, int len);
int evaluate_penalty (Alignment *A, Constraint_list *CL, int *scale,char *scale_mode, int *penalty, char *penalty_mode, int len_seq);
Alignment ** t_coffee_lalign   (Constraint_list *CL, int scale, int penalty,int maximise,Sequence *S, int sw_t, int sw_l, int sw_z,int *sw_n, int sw_io);
Alignment * add_seq2aln   (Constraint_list *CL, Alignment *IN,Sequence  *S);







