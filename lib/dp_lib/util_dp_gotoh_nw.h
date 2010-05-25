/*pair wise aln implementations*/

int idscore_pairseq (char *s1, char *s2, int gop, int gep, int **m, char *mode);
int idscore_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL);
int gotoh_pair_wise    (Alignment *A,int*ns, int **l_s,Constraint_list *CL);
int glocal_pair_wise ( Alignment *A, int *ns, int **l_s, Constraint_list *CL);
int gotoh_pair_wise_lgp ( Alignment *A, int *ns, int **l_s, Constraint_list *CL);
int test_pair_wise (Alignment *A, int *ns, int **l_s, Constraint_list *CL);

int glocal2_pair_wise (Alignment *A,int*ns, int **ls,Constraint_list *CL);
int gotoh_pair_wise_lgp_sticky ( Alignment *A, int *ns, int **l_s, Constraint_list *CL);
int linked_pair_wise ( Alignment *A, int *ns, int **l_s, Constraint_list *CL);
int clinked_pair_wise ( Alignment *A, int *ns, int **l_s, Constraint_list *CL);
