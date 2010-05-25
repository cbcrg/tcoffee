int ** evaluate_diagonals ( Alignment *A, int *ns, int **l_s, Constraint_list *CL,int maximise,int n_groups, char **group_list,int ktup);
int ** evaluate_segments_with_ktup  ( Alignment *A, int *ns, int **l_s, Constraint_list *CL,int maximise,int n_groups, char **group_list,int ktup);
int ** evaluate_diagonals_with_ktup ( Alignment *A, int *ns, int **l_s, Constraint_list *CL,int maximise,int n_groups, char **group_list,int ktup);

int ** evaluate_diagonals_with_clist ( Alignment *A, int *ns, int **l_s, Constraint_list *CL,int maximise,int n_groups, char **group_list,int ktup);

int * flag_diagonals (int l1, int l2, int **sorted_diag,float T);                
int hasch_seq(char *seq1, int **hs, int **lu,int ktup, char *alph);
int fasta_gotoh_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL);
int make_fasta_gotoh_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL, int *diag);
