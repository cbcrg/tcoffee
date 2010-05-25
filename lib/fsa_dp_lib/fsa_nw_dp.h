int quadratic_fsa_nw_pair_wise (Alignment *A, int *ns, int **ls, Constraint_list *CL);
int fsa_nw_pair_wise (Alignment *A, int *ns, int **ls, Constraint_list *CL);
Fsa_dp_model * initialize_fsa_nw_model(int tl_s1, int tl_s2,int tr_s1, int tr_s2, int tg_mode, Constraint_list *CL);


int fsa_nw_get_sub_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int fsa_nw_get_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int fsa_nw_get_term_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int fsa_nw_get_start_gep_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
int fsa_nw_get_no_cost (Alignment *A, int**pos1, int ns1, int*list1, int col1, int**pos2, int ns2, int*list2, int col2, Constraint_list *CL);
