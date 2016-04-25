#ifndef __UTIL_DP_GOTOH_NW_H
#define __UTIL_DP_GOTOH_NW_H 
/*pair wise aln implementations*/

int idscore_pairseq (char *s1, char *s2, int gop, int gep, int **m, char *mode);
int idscore_pair_wise (Alignment *A,int*ns, int **l_s,Constraint_list *CL);
int gotoh_pair_wise         (Alignment *A,int*ns, int **l_s,Constraint_list *CL);
int gotoh_pair_wise_test    (Alignment *A,int*ns, int **l_s,Constraint_list *CL);

int glocal_pair_wise ( Alignment *A, int *ns, int **l_s, Constraint_list *CL);
int gotoh_pair_wise_lgp ( Alignment *A, int *ns, int **l_s, Constraint_list *CL);
int test_pair_wise (Alignment *A, int *ns, int **l_s, Constraint_list *CL);

int glocal2_pair_wise (Alignment *A,int*ns, int **ls,Constraint_list *CL);
int linked_pair_wise ( Alignment *A, int *ns, int **l_s, Constraint_list *CL);
int linked_pair_wise_collapse ( Alignment *A, int *ns, int **l_s, Constraint_list *CL);
int cl2pair_list_ecl_ext_pc ( Alignment *A, int *ns, int **ls, Constraint_list *CL, int ***list_in, int *n_in);
int procoffee_pair_wise ( Alignment *A, int *nsi, int **lsi, Constraint_list *CL);


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////                                                                           ///////////////
///////////////                                                                           ///////////////
///////////////                       Profile_PW                                          ///////////////
///////////////                                                                           ///////////////
///////////////                                                                           ///////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

int **aln2prf (Alignment *A, char *matrix);
char *add_sequence2prf (char *seq, char *al,int **prf, int lj, int gop, int gep);

#endif 
